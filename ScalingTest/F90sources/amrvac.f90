!================================================================================
!
!    BHAC (The Black Hole Accretion Code) solves the equations of
!    general relativistic magnetohydrodynamics and other hyperbolic systems
!    in curved spacetimes.
!
!    Copyright (C) 2019 Oliver Porth, Hector Olivares, Yosuke Mizuno, Ziri Younsi,
!    Luciano Rezzolla, Elias Most, Bart Ripperda and Fabio Bacchini
!
!    This file is part of BHAC.
!
!    BHAC is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    BHAC is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with BHAC.  If not, see <https://www.gnu.org/licenses/>.
!
!================================================================================

!=============================================================================
program amrvac

  use mod_multigrid_coupling

! AMRVAC solves a set of hyperbolic equations:
!              u  +  Div [f(u)]    =  s
!               t       x    
!  using adaptive mesh refinement.

! following line may avoid mystery problems in the combination
!  where we use ifort compiler and MPT (SGIs MPI implementation)
!DEC$ ATTRIBUTES NOINLINE :: read_snapshot


use mod_amrvacdef

integer :: itin
double precision :: time0, time_in, tin
!-----------------------------------------------------------------------------
call comm_start

time_advance=.false.

time0=MPI_WTIME()
time_bc=zero


! get input parameters
call readcommandline
call readparameters
call initialize_vars
call init_comm_types


! Begin of the code
! -----------------
if (snapshotini/=-1) then
   ! restart from previous file or dat file conversion
   ! get input data from previous VAC/AMRVAC run
   itin=it
   tin=t

   ! order: first physics dependent part then user for defaults.
   call initglobal





   select case (typeparIO)
   case (1,-2)
      call read_snapshot
   case (0,-1)
      call read_snapshot_nopar
   end select



   ! modify globals
   if (changeglobals) call initglobal





   if (convert) then
       
       if (npe/=1.and.(.not.(index(convert_type,'mpi')>=1)) &
          .and. convert_type .ne. 'user')  call mpistop("non-mpi conversion &
          only uses 1 cpu")
       call generate_plotfile

       call comm_finalize
       stop
   end if

   if (itreset) it=itin
   if (treset) t=tin
   ! modify initial condition
   if (firstprocess) call modify_IC
   ! reset AMR grid
   if (resetgrid) call settree
   ! Select active grids
   call selectgrids

   ! set up boundary flux conservation arrays
   if (levmax>levmin) call allocateBflux
 
else
   ! order: first physics dependent part then user for defaults.
   call initglobal



   ! form and initialize all grids at level one
   call initlevelone
   ! set up and initialize finer level grids, if needed
   call settree

   ! set up boundary flux conservation arrays
   if (levmax>levmin) call allocateBflux
   


select case(clean_init_divB)
case ('vecpot')
 !Now that the grid has been set and the flux-conservation arrays are allocated,
   ! re-calculate the magnetic field from the vector potential in a completely
   ! divergence free way.
   if (levmax>levmin .and. ndim .eq. 3) call recalculateB

case ('mg')
   ! Project out the divB using poisson solver.
   ! Due to Jannis Teunissen. Thanks!   
   call mg_setup_multigrid()
   call clean_divb_multigrid()

case ('none')
   ! do nothing
case default
   call mpistop('Unknown method selected for clean_init_divB!')
   ! check what you are doing!
end select




   




end if

if (mype==0) then
   print*,'-----------------------------------------------------------------------------'
   write(*,'(a,f17.3,a)')' Startup phase took : ',MPI_WTIME()-time0,' sec'
   print*,'-----------------------------------------------------------------------------'
end if

time_advance=.true.

! do time integration of all grids on all levels
call timeintegration

if (mype==0) then
   print*,'-----------------------------------------------------------------------------'
   write(*,'(a,f17.3,a)')' Finished BHAC in : ',MPI_WTIME()-time0,' sec'
   print*,'-----------------------------------------------------------------------------'
end if


call comm_finalize

end program amrvac
!=============================================================================
subroutine timeintegration

use mod_timing
use mod_amrvacdef

integer :: level, ifile, fixcount

logical, external :: timetosave, fixgrid

!-----------------------------------------------------------------------------
time_in=MPI_WTIME()
fixcount=1

n_saves(filelog_:fileout_) = snapshotini

do ifile=nfile,1,-1
   if (time_accurate) tsavelast(ifile)=t
   itsavelast(ifile)=it
end do

itmin=it
! the next two are used to keep track of the performance during runtime:
itTimeLast=it
timeLast=MPI_WTIME()

call getbc(t,ps,psCoarse)
!  ------ start of integration loop. ------------------
timeloop0=MPI_WTIME()
timefirstbc=timeloop0-time_in
time_bc=0.d0
if (mype==0) then
   write(*,'(a,f12.3,a)')' BCs before Advance took : ',timefirstbc,' sec'
end if

time_evol : do
   ! Select active grids
   call selectgrids
   call setdt
   if(fixprocess) call process(it,t)

   timeio0=MPI_WTIME()
   do ifile=nfile,1,-1
      if(timetosave(ifile)) call saveamrfile(ifile)
   end do

   timeio_tot=timeio_tot+(MPI_WTIME()-timeio0)

   ! exit time loop criteria
   if (it>=itmax) exit time_evol
   if (time_accurate .and. t>=tmax) exit time_evol

   call advance(it)

   if((.not.time_accurate).or.(residmin>smalldouble)) then
      call getresidual(it)
   endif 
   if (residual<residmin .or. residual>residmax) exit time_evol

   ! resetting of tree BEFORE IO and setdt
   timegr0=MPI_WTIME()
   if(ditregrid>1) then
     if(fixcount<ditregrid) then
       fixcount=fixcount+1
     else
       if (mxnest>1 .and. .not.(fixgrid(0))) call resettree
       fixcount=1
     endif
   else
     if (mxnest>1 .and. .not.(fixgrid(0))) call resettree
   endif
   timegr_tot=timegr_tot+(MPI_WTIME()-timegr0)

   it = it + 1
   if (time_accurate) t = t + dt
   if(addmpibarrier) call MPI_BARRIER(icomm,ierrmpi)
   
   if(it>9000000)then
     it = slowsteps+10
     itmin=0
     itsavelast(:)=0
   end if
   
end do time_evol

timeloop=MPI_WTIME()-timeloop0

if (mype==0) then
   write(*,'(a,f12.3,a)')' Total timeloop took        : ',timeloop,' sec'
   write(*,'(a,f12.3,a)')' Time spent on Regrid+Update: ',timegr_tot,' sec'
   write(*,'(a,f12.2,a)')'                  Percentage: ',100.0*timegr_tot&
      /timeloop,' %'
   write(*,'(a,f12.3,a)')' Time spent on IO in loop   : ',timeio_tot,' sec'
   write(*,'(a,f12.2,a)')'                  Percentage: ',100.0*timeio_tot&
      /timeloop,' %'
   write(*,'(a,f12.3,a)')' Time spent on BC           : ',time_bc,' sec'
   write(*,'(a,f12.2,a)')'                  Percentage: ',100.0*time_bc&
      /timeloop,' %'
end if

!timeio0=MPI_WTIME()
!do ifile=nfile,1,-1
!   if(itsavelast(ifile)<it)call saveamrfile(ifile)
!enddo
!if (mype==0) call MPI_FILE_CLOSE(log_fh,ierrmpi)
!timeio_tot=timeio_tot+(MPI_WTIME()-timeio0)




if (mype==0) then
   write(*,'(a,f12.3,a)')' Total time spent on IO     : ',timeio_tot,' sec'
   write(*,'(a,f12.3,a)')' Total timeintegration took : ',MPI_WTIME()-time_in,&
      ' sec'
end if

end subroutine timeintegration
!=============================================================================
logical function timetosave(ifile)

! Save times are defined by either tsave(isavet(ifile),ifile) or
! itsave(isaveit(ifile),ifile) or dtsave(ifile) or ditsave(ifile)
! Other conditions may be included.

use mod_amrvacdef

integer:: ifile
logical:: oksave
!-----------------------------------------------------------------------------
oksave=.false.
if (it==itsave(isaveit(ifile),ifile)) then
   oksave=.true.
   isaveit(ifile)=isaveit(ifile)+1
end if
if (it==itsavelast(ifile)+ditsave(ifile)) oksave=.true.
if (time_accurate) then
   if (t>=tsave(isavet(ifile),ifile)) then
      oksave=.true.
      isavet(ifile)=isavet(ifile)+1
   end if
   if (t>=tsavelast(ifile)+dtsave(ifile)-smalldouble)then
     oksave=.true.
     n_saves(ifile) = n_saves(ifile) + 1
   endif
end if
if (oksave) then
   if (time_accurate) tsavelast(ifile) =t
   itsavelast(ifile)=it
end if
timetosave=oksave

return
end function timetosave
!=============================================================================
logical function fixgrid(dummy)

! fixing the grid at given times or given iteration is defined by either 
! tfixgrid or itfixgrid
! Other conditions may be included (dummy input integer unused).

use mod_amrvacdef

integer:: dummy
!-----------------------------------------------------------------------------

fixgrid= (t>=tfixgrid .or. it>=itfixgrid)

return
end function fixgrid
!=============================================================================
subroutine initglobal

use mod_amrvacdef
!-----------------------------------------------------------------------------

call initglobaldata
call initglobaldata_usr
call checkglobaldata

end subroutine initglobal
!=============================================================================
