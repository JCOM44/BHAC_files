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
subroutine advance(iit)


use mod_particles, only: handle_particles

use mod_amrvacdef
integer, intent(in) :: iit

integer :: iigrid, igrid, idimsplit

logical :: firstsweep, lastsweep
common/first/firstsweep,lastsweep
!-----------------------------------------------------------------------------

! when computing to steady state, store old solution
if(residmin>smalldouble)then
!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
  do iigrid=1,igridstail; igrid=igrids(iigrid);
     call copy_state(ps(igrid),psres(igrid))
  end do
!$OMP END PARALLEL DO
endif



! split source addition
if(ssplitdivb .or. ssplitresis .or. ssplituser) call addsource_all(.true.)


! old solution values at t_n-1 no longer needed: make copy of w(t_n)
!$OMP PARALLEL DO SCHEDULE(runtime) PRIVATE(igrid)
do iigrid=1,igridstail; igrid=igrids(iigrid);
   call copy_state(ps(igrid),psold(igrid))
end do
!$OMP END PARALLEL DO

! split implicit source addition
if (sourceimpl) then
   if(mype==0.and..false.) print *,'implicit source addition'
   call addsource_impl
endif

firstsweep=.true.
if (dimsplit) then
   if ((iit/2)*2==iit .or. typedimsplit=='xy') then
      ! do the sweeps in order of increasing idim,
      do idimsplit=1,ndim
         lastsweep= idimsplit==ndim
         call advect(idimsplit,idimsplit)
      end do
   else
      ! If the parity of "iit" is odd and typedimsplit=xyyx,
      ! do sweeps backwards
      do idimsplit=ndim,1,-1
         lastsweep= idimsplit==1
         call advect(idimsplit,idimsplit)
      end do
   end if
else
   ! Add fluxes from all directions at once
   lastsweep= .true.
   call advect(1,ndim)
end if


! split source addition
if(ssplitdivb .or. ssplitresis .or. ssplituser) call addsource_all(.false.)



if (use_particles) call handle_particles()


if(residmin>smalldouble)then
!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
  do iigrid=1,igridstail; igrid=igrids(iigrid);
     pwres(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)&
        = pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
        1:nwflux)-pwres(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
        1:nwflux)
     
     pwsres(igrid)%w = pws(igrid)%w-pwsres(igrid)%w
     
  end do
!$OMP END PARALLEL DO
endif

end subroutine advance
!=============================================================================
subroutine advect(idimmin,idimmax)

!  integrate all grids by one step of its delta(t)

! This subroutine is in VAC terminology equivalent to
! `advect' (with the difference that it will `advect' all grids)

use mod_amrvacdef

integer, intent(in) :: idimmin,idimmax

integer :: iigrid, igrid
logical :: firstsweep, lastsweep
common/first/firstsweep,lastsweep
!-----------------------------------------------------------------------------
! copy w instead of wold because of potential use of dimsplit or sourcesplit
!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
do iigrid=1,igridstail; igrid=igrids(iigrid);
   call copy_state(ps(igrid),ps1(igrid))
end do
!$OMP END PARALLEL DO

istep=0

select case (typeadvance)
case ("onestep")
   call advect1(typefull1,one,    idimmin,idimmax,t,          t, ps1,ps,psold)
case ("twostep")
   ! predictor step
   call advect1(typepred1,half,   idimmin,idimmax,t,          t, ps,ps1,psold)
   ! corrector step
   call advect1(typefull1,one,    idimmin,idimmax,t+half*dt,  t, ps1,ps,psold)
case ("ImEx12")
   ! ImEx method, e.g. Palenzuela et al. 2009 for srrmhd and Palenzuela et al,
   ! 2013
   ! and Palenzuela 2013 for grmhd.
   ! This is the first-second pseudo-ImEx scheme from Bucciantini and Del Zanna
   ! 2013

   ! predictor step
   call advect1(typepred1,half,   idimmin,idimmax,t,          t, ps,ps1,psold)
   ! corrector step
   call advect1(typefull1,one,    idimmin,idimmax,t+half*dt,  t, ps1,ps,psold)
case ("ImEx2")
   !This is the second order ImEx scheme IMEX-SSP2 (2,2,2)
   !(Palenzuela et al. 2009)

   call mpistop('This ImEx scheme is not implemented yet')
case("ImEx3")
   !This is the third order ImEx scheme IMEX-SSP3 (4,3,3)
   !(Palenzuela et al. 2009)

   call mpistop('This ImEx scheme is not implemented yet')
case ("threestep")
   ! three step Runge-Kutta in accordance with Gottlieb & Shu 1998
   call advect1(typefull1,one, idimmin,idimmax,t,             t, ps,ps1,psold)

   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call alloc_state(ps2(igrid),ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
         ixGhi3)      
   end do
   
!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw2(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)&
         =0.75d0*pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
         1:nwflux)+0.25d0*pw1(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
         ixGlo3:ixGhi3,1:nwflux)
      if (nw>nwflux) pw2(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
         nwflux+1:nw) = pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
         nwflux+1:nw)
      
      pws2(igrid)%w = 0.75d0*pws(igrid)%w + 0.25d0*pws1(igrid)%w
      
   end do
!$OMP END PARALLEL DO

   call advect1(typefull1,0.25d0, idimmin,idimmax,t+dt, t+dt*0.25d0, ps1,ps2,&
      psold)

!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
   do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)=1.0d0&
         /3.0d0*pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
         1:nwflux)+2.0d0/3.0d0*pw2(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
         ixGlo3:ixGhi3,1:nwflux)
      
      pws(igrid)%w=1.0d0/3.0d0*pws(igrid)%w + 2.0d0/3.0d0*pws2(igrid)%w
      
   end do   
!$OMP END PARALLEL DO
   
   call advect1(typefull1,2.0d0/3.0d0, idimmin,idimmax,t+dt/2.0d0, t+dt&
      /3.0d0,  ps2,ps,psold)

case ("ssprk43")
! Strong stability preserving 4 stage RK 3rd order method by Ruuth and Spiteri
!
! Ruuth & Spiteri
! J. S C, 17 (2002) p. 211 - 220
!
! supposed to be stable up to CFL=2.
! don't use time-dependent sources since i did not bother to set those intermediate times.
! oliver.

! === First step ===
   call advect1(typefull1,0.5d0, idimmin,idimmax,t,           t, ps,ps1,psold)

! === Second step ===
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call alloc_state(ps2(igrid),ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3)
   end do

!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw2(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)&
         =pw1(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)
      if (nw>nwflux) pw2(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
         nwflux+1:nw) = pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
         nwflux+1:nw)
      
      pws2(igrid)%w = pws1(igrid)%w
      
   end do
   !$OMP END PARALLEL DO
   
   call advect1(typefull1,0.5d0, idimmin,idimmax,t,           t+dt, ps1,ps2,&
      psold)

! === Third step ===
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call alloc_state(ps3(igrid),ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3)
   end do

!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw3(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)=2.0d0&
         /3.0d0 * pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
         1:nwflux) + 1.0d0/3.0d0 * pw2(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
         ixGlo3:ixGhi3,1:nwflux)
      if (nw>nwflux) pw3(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
         nwflux+1:nw) = pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
         nwflux+1:nw)
      
      pws3(igrid)%w=2.0d0/3.0d0 * pws(igrid)%w + 1.0d0/3.0d0 * pws2(igrid)%w
      
   end do
   !$OMP END PARALLEL DO
   
   call advect1(typefull1,1.0d0/6.0d0, idimmin,idimmax,t,       t+dt, ps2,ps3,&
      psold)

! === Fourth step ===
!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
   do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)&
         =pw3(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)
      
      pws(igrid)%w=pws3(igrid)%w
      
   end do   
!$OMP END PARALLEL DO
   
   call advect1(typefull1,0.5d0, idimmin,idimmax, t,            t+dt, ps3,ps,&
      psold)



case ("ssprk54")
! Strong stability preserving 5 stage RK 4th order method by Ruuth and Spiteri
!
! SIAM J. NUMER. ANAL.
! Vol. 40, No. 2, pp. 469–491
! c 2002 Society for Industrial and Applied Mathematics
! A NEW CLASS OF OPTIMAL
! HIGH-ORDER STRONG-STABILITY-PRESERVING
! TIME DISCRETIZATION METHODS
!
! E.g. Table A.2
!
! I have the actual coefficients however from the overview article by Gottlieb, JoSC 25 (2005)
! ON HIGH ORDER STRONG STABILITY PRESERVING RUNGE-KUTTA AND MULTI STEP TIME DISCRETIZATIONS
!
! there are slight differences in the coefficients (~8th digit after the .)
! This is SSP till CFL number 1.508 which makes the effective CFL per step ceff=0.377, 
! in contrast to the ceff=0.3333 of the classical RK3.  
!
! coded by oliver on 11/05/2013.  Enjoy!

! === First step ===
   call advect1(typefull1,0.391752226571890d0, idimmin,idimmax, t ,t, ps,ps1,&
      psold)

! === Second step ===
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call alloc_state(ps2(igrid),ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3)
   end do

!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw2(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)&
         =0.444370493651235d0 * pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
         ixGlo3:ixGhi3,1:nwflux) + 0.555629506348765d0 * pw1(igrid)%w&
         (ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)
      if (nw>nwflux) pw2(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
         nwflux+1:nw) = pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
         nwflux+1:nw)
      
      pws2(igrid)%w=0.444370493651235d0 * pws(igrid)%w + 0.555629506348765d0 &
         * pws1(igrid)%w
      
   end do
!$OMP END PARALLEL DO
   
   call advect1(typefull1,0.368410593050371d0, idimmin,idimmax,&
       t+0.391752226571890d0*dt, t+0.2176690962611688d0*dt, ps1,ps2,psold)

! === Third step ===
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call alloc_state(ps3(igrid),ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3)
   end do

!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw3(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)&
         =0.620101851488403d0 * pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
         ixGlo3:ixGhi3,1:nwflux) + 0.379898148511597d0 * pw2(igrid)%w&
         (ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)
      if (nw>nwflux) pw3(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
         nwflux+1:nw) = pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
         nwflux+1:nw)
      
      pws3(igrid)%w=0.620101851488403d0 * pws(igrid)%w + 0.379898148511597d0 &
         * pws2(igrid)%w
      

   end do
!$OMP END PARALLEL DO
   
   call advect1(typefull1,0.251891774271694d0, idimmin,idimmax,&
       t+0.5860796893115398d0*dt, t+0.222650588849706d0*dt, ps2,ps3,psold)

! === Fourth step ===
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call alloc_state(ps4(igrid),ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3)
   end do

!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw4(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)&
         =0.178079954393132d0 * pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
         ixGlo3:ixGhi3,1:nwflux) + 0.821920045606868d0 * pw3(igrid)%w&
         (ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)
      if (nw>nwflux) pw4(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
         nwflux+1:nw) = pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
         nwflux+1:nw)
      
      pws4(igrid)%w=0.178079954393132d0 * pws(igrid)%w + 0.821920045606868d0 &
         * pws3(igrid)%w
      
   end do
!$OMP END PARALLEL DO
   
   call advect1(typefull1,0.544974750228521d0, idimmin,idimmax,&
       t+0.4745423631214d0*dt, t+0.390035880739132d0*dt, ps3,ps4,psold)

! Now recover back the dt*L(u3), store in pw1:
!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw1(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux) &
         = ( pw4(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
         1:nwflux) - (0.178079954393132d0 * pw(igrid)%w(ixGlo1:ixGhi1,&
         ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux) + 0.821920045606868d0 * &
         pw3(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)) ) &
         / 0.544974750228521d0
      
      pws1(igrid)%w = ( pws4(igrid)%w &
           - (0.178079954393132d0 * pws(igrid)%w &
           + 0.821920045606868d0 * pws3(igrid)%w) ) / 0.544974750228521d0
      
   end do
!$OMP END PARALLEL DO

! === Fifth step ===
!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)&
         = 0.517231671970585d0 * pw2(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
         ixGlo3:ixGhi3,1:nwflux) + 0.096059710526147d0 * pw3(igrid)%w&
         (ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux) + &
         0.063692468666290d0 * pw1(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
         ixGlo3:ixGhi3,1:nwflux) + 0.386708617503269d0 * pw4(igrid)%w&
         (ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)
      
      pws(igrid)%w= 0.517231671970585d0 * pws2(igrid)%w &
           + 0.096059710526147d0 * pws3(igrid)%w &
           + 0.063692468666290d0 * pws1(igrid)%w &
           + 0.386708617503269d0 * pws4(igrid)%w
      
   end do
!$OMP END PARALLEL DO
   call advect1(typefull1,0.226007483236906d0, idimmin,idimmax,&
       t+0.935010630967653d0*dt, t+0.710300048096804d0*dt, ps4,ps,psold)



case ("rk4")
   ! classical RK4 four step scheme
   call advect1(typefull1,0.5d0, idimmin,idimmax,t,         t, ps,ps1,psold)
   
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call alloc_state(ps2(igrid),ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3)
   end do

!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw2(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)&
         =pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)
      
      pws2(igrid)%w=pws(igrid)%w
      
   end do
!$OMP END PARALLEL DO
   
   call advect1(typefull1,0.5d0, idimmin,idimmax,t+dt/2d0,   t, ps1,ps2,psold)


   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call alloc_state(ps3(igrid),ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3)
   end do

!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw3(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)&
         =pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)
      
      pws3(igrid)%w=pws(igrid)%w
      
   end do
!$OMP END PARALLEL DO
   
   call advect1(typefull1,one,   idimmin,idimmax,t+dt/2d0,   t, ps2,ps3,psold)


!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw1(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)&
         =(pw1(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
         1:nwflux) +two*pw2(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
         ixGlo3:ixGhi3,1:nwflux) +pw3(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
         ixGlo3:ixGhi3,1:nwflux) -4.0d0*pw(igrid)%w(ixGlo1:ixGhi1,&
         ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux))/3.0d0
      
      pws1(igrid)%w = ( pws1(igrid)%w + two*pws2(igrid)%w &
           + pws3(igrid)%w -4.0d0*pws(igrid)%w )/3.0d0
      
   end do
!$OMP END PARALLEL DO
   
   call advect1(typefull1,1.0d0/6.0d0, idimmin,idimmax,t+dt,  t, ps3,ps,psold)

!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)&
         =pw1(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
         1:nwflux)+pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
         1:nwflux)
      
      pws(igrid)%w = pws1(igrid)%w + pws(igrid)%w
      
   end do
!$OMP END PARALLEL DO


case ("fourstep")
   ! four step scheme, variant Hans De Sterck
   call advect1(typefull1,0.12d0, idimmin,idimmax,t,                t, ps,ps1,&
      psold)
   
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call alloc_state(ps2(igrid),ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3)
   end do

!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw2(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)&
         =pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)
      
      pws2(igrid)%w = pws(igrid)%w
      
   end do
!$OMP END PARALLEL DO
   
   call advect1(typefull1,0.25d0, idimmin,idimmax,t+dt*0.12d0,      t, ps1,&
       ps2,psold)
   
!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw1(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)&
         =pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)
      
      pws1(igrid)%w = pws(igrid)%w
      
   end do
!$OMP END PARALLEL DO
   
   call advect1(typefull1,0.5d0,  idimmin,idimmax,t+dt/4d0   ,t, ps2,ps1,&
      psold)
   call advect1(typefull1,one,    idimmin,idimmax,t+dt/2d0   ,t, ps1,ps,psold)


case ("jameson")
   ! four step scheme, variant jameson
   call advect1(typefull1,0.25d0, idimmin,idimmax,t,                  t, ps,&
      ps1,psold)
   
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call alloc_state(ps2(igrid),ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3)
   end do

!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw2(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)&
         =pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)
      
      pws2(igrid)%w = pws(igrid)%w
      
   end do
!$OMP END PARALLEL DO
   
   call advect1(typefull1,(1.0d0/3.0d0), idimmin,idimmax,t+dt*0.25d0, t, ps1,&
      ps2,psold)

!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw1(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)&
         =pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)
      
      pws1(igrid)%w = pws(igrid)%w
      
   end do
!$OMP END PARALLEL DO
   
   call advect1(typefull1,0.5d0,  idimmin,idimmax,t+dt/3d0,            t, ps2,&
      ps1,psold)
   call advect1(typefull1,one,    idimmin,idimmax,t+dt/2d0,            t, ps1,&
      ps,psold)


case default
   write(unitterm,*) "typeadvance=",typeadvance
   write(unitterm,*) "Error in advect: Unknown time integration method"
   call mpistop("Correct typeadvance")
end select


firstsweep=.false.
end subroutine advect
!=============================================================================
subroutine advect1(method,dtfactor,idimmin,idimmax,qtC,qt,psa,psb,psc)

!  integrate all grids by one partial step

! This subroutine is equivalent to VAC's `advect1', but does
! the advection for all grids
use mod_amrvacdef

integer, intent(in)          :: idimmin,idimmax
double precision, intent(in) :: dtfactor, qtC, qt
character(len=*), intent(in) :: method(nlevelshi)
type(state)                  :: psa(ngridshi), psb(ngridshi), psc(ngridshi)

! .. local ..
double precision             :: qdt
integer                      :: iigrid, igrid, level
logical                      :: do_fluxfix
!-----------------------------------------------------------------------------
istep=istep+1

do_fluxfix = ((time_advance.and.levmax>levmin) .and. (istep==nstep&
   .or.nstep>2))


if (do_fluxfix) call init_comm_fix_conserve(idimmin,idimmax)


!$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(igrid,level,qdt)
do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
   level=node(plevel_,igrid)
   qdt=dtfactor*dt_grid(igrid)

   call process1_grid(method(level),igrid,qdt,ixGlo1,ixGlo2,ixGlo3,ixGhi1,&
      ixGhi2,ixGhi3,idimmin,idimmax,qtC,psa(igrid),qt,psb(igrid),psc(igrid))

end do
!$OMP END PARALLEL DO


if (do_fluxfix) then

   do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      call sendflux(igrid,idimmin,idimmax) ! OMP: not threadsafe!
   end do

   call fix_conserve(psb,idimmin,idimmax)
   
   call fix_edges(psb,idimmin,idimmax)
   
end if



! Now we fill the centers for the staggered variables
!$OMP PARALLEL DO SCHEDULE(runtime) PRIVATE(igrid)
do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
   call faces2centers(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,psb(igrid))
end do
!$OMP END PARALLEL DO


! Fill ghost cells
qdt=dtfactor*dt
call getbc(qt+qdt,psb,psCoarse)

end subroutine advect1
!=============================================================================
subroutine process1_grid(method,igrid,qdt,ixImin1,ixImin2,ixImin3,ixImax1,&
   ixImax2,ixImax3,idimmin,idimmax,qtC,sCT,qt,s,sold)

! This subroutine is equivalent to VAC's `advect1' for one grid

use mod_amrvacdef

character(len=*), intent(in) :: method
integer, intent(in)          :: igrid, ixImin1,ixImin2,ixImin3,ixImax1,&
   ixImax2,ixImax3, idimmin,idimmax
double precision, intent(in) :: qdt, qtC, qt
type(state)                  :: sCT, s, sold
! .. local ..
double precision             :: dx1,dx2,dx3
double precision             :: fC(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nwflux,1:ndim)
logical                      :: do_fluxfix
double precision             :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:3)
!-----------------------------------------------------------------------------
do_fluxfix = ((time_advance.and.levmax>levmin) .and. (istep==nstep&
   .or.nstep>2))


dx1=rnode(rpdx1_,igrid);dx2=rnode(rpdx2_,igrid);dx3=rnode(rpdx3_,igrid);
call set_tmpGlobals(igrid)

call advect1_grid(method,qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   idimmin,idimmax,qtC,sCT,qt,s,sold,fC,fE,dx1,dx2,dx3, px(igrid)%x)


if (do_fluxfix) then
   call storeflux(igrid,fC,idimmin,idimmax)

   call storeedge(igrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,fE,&
      idimmin,idimmax) 
end if

end subroutine process1_grid
!=============================================================================
subroutine advect1_grid(method,qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3,idimmin,idimmax,qtC,sCT,qt,s,sold,fC,fE,dx1,dx2,dx3,x)

!  integrate one grid by one partial step 

use mod_amrvacdef

character(len=*), intent(in) :: method
integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
    idimmin,idimmax
double precision, intent(in) :: qdt, qtC, qt, dx1,dx2,dx3, x(ixImin1:ixImax1,&
   ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
type(state)      :: sCT, s, sold
double precision :: fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nwflux,1:ndim)
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nw) :: wprim
double precision :: fE(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:3)
integer          :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
!integer            :: ixtile^D, ntiles^D ! number of tiles
!integer, parameter :: ntile1=32, ntile2=16, ntile3=4 ! number of cells in tile
!-----------------------------------------------------------------------------

ixOmin1=ixImin1+dixB;ixOmin2=ixImin2+dixB;ixOmin3=ixImin3+dixB
ixOmax1=ixImax1-dixB;ixOmax2=ixImax2-dixB;ixOmax3=ixImax3-dixB;

! Always primitive reconstruction now
if (covariant) myM => mygeo%m
call primitive(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixImin1,&
   ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,sCT%w%w,x) !also applies floors to wCT
wprim = sCT%w%w
call conserve(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixImin1,ixImin2,&
   ixImin3,ixImax1,ixImax2,ixImax3,sCT%w%w,x,patchfalse)

! Now tile it:
!ntiles^D=(ixImax^D-ixImin^D-2*dixB+1)/ntile^D;

!{^D& do ixtile^DB = 1, ntiles^DB\}

!{^D& ixOmin^D = ixImin^D + dixB + (ixtile^D-1)*ntile^D \}
!{^D& ixOmax^D = ixOmin^D + ntile^D - 1 \}

select case (method)
case ('tvdlf','tvdlf1')
   call tvdlf(method,qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idimmin,idimmax,qtC,sCT,&
      qt,s,sold,wprim,fC,fE,dx1,dx2,dx3,x)
case ('tvdlfpos')
   call tvdlfpos(method,qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idimmin,idimmax,qtC,sCT,&
      qt,s,sold,wprim,fC,fE,dx1,dx2,dx3,x)
case ('hll','hll1')
   call hll(method,qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idimmin,idimmax,qtC,sCT,&
      qt,s,sold,wprim,fC,fE,dx1,dx2,dx3,x)
case ('source')
   call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim),ixImin1,ixImin2,&
      ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
      ixOmax3,1,nw,qtC,sCT%w%w,qt,s%w%w,x,.false.)
case ('nul')
   ! There is nothing to do
case default
   write(unitterm,*)'Error in advect1_grid:',method,' is unknown!'
   call mpistop("")
end select

!{end do\}

end subroutine advect1_grid
!=============================================================================
subroutine addsource_all(prior)

use mod_amrvacdef

logical, intent(in) :: prior

integer :: iigrid, igrid
double precision :: qdt, qt
!-----------------------------------------------------------------------------

if ((.not.prior).and.(typesourcesplit=='sf' .or. typesourcesplit&
   =='ssf')) return

if (prior) then
   qt=t
else
   qt=t+dt
end if

!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid,qdt)
      do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
         qdt=dt_grid(igrid)
         if (B0field) then
            myB0_cell => pB0_cell(igrid)
            myB0_face1 => pB0_face1(igrid)
            myB0_face2 => pB0_face2(igrid)
            myB0_face3 => pB0_face3(igrid)
         end if
         call addsource1_grid(igrid,qdt,qt,pw(igrid)%w)
      end do
!$OMP END PARALLEL DO
      
      call getbc(qt,ps,psCoarse)

end subroutine addsource_all
!=============================================================================
subroutine addsource1_grid(igrid,qdt,qt,w)

use mod_amrvacdef

integer, intent(in) :: igrid
double precision, intent(in) :: qdt, qt
double precision, intent(inout) :: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
   ixGlo3:ixGhi3,nw)

double precision :: w1(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw)
!-----------------------------------------------------------------------------

call set_tmpGlobals(igrid)

w1(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nw)=w(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nw)

select case (typesourcesplit)
case ('sf')
   call addsource2(qdt  ,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixMlo1,&
      ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,1,nw,qt,w1,qt,w,px(igrid)%x,.true.)
case ('sfs')
   call addsource2(qdt/2,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixMlo1,&
      ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,1,nw,qt,w1,qt,w,px(igrid)%x,.true.)
case ('ssf')
   call addsource2(qdt/2,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixGlo1,&
      ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,1,nw,qt,w,qt,w1,px(igrid)%x,.true.)
   call addsource2(qdt  ,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixMlo1,&
      ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,1,nw,qt,w1,qt,w,px(igrid)%x,.true.)
case ('ssfss')
   call addsource2(qdt/4,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixGlo1,&
      ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,1,nw,qt,w,qt,w1,px(igrid)%x,.true.)
   call addsource2(qdt/2,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixMlo1,&
      ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,1,nw,qt,w1,qt,w,px(igrid)%x,.true.)
case default
   write(unitterm,*)'No such typesourcesplit=',typesourcesplit
   call mpistop("Error: Unknown typesourcesplit!")
end select

end subroutine addsource1_grid
!=============================================================================
subroutine addsource2(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iwmin,iwmax,qtC,wCT,qt,w,x,&
   qsourcesplit)

! Add source within ixO for iws: w=w+qdt*S[wCT]

use mod_amrvacdef

! differences with VAC is in iw^LIM and in declaration of ranges for wCT,w

integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
    ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, iwmin,iwmax
double precision, intent(in) :: qdt, qtC, qt, x(ixImin1:ixImax1,&
   ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
double precision, intent(inout) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw), w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nw)
logical, intent(in) :: qsourcesplit
!-----------------------------------------------------------------------------

! user defined sources, typically explicitly added
if(qsourcesplit .eqv. ssplituser) call specialsource(qdt,ixImin1,ixImin2,&
   ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
   ixOmax3,iwmin,iwmax,qtC,wCT,qt,w,x)

! physics defined sources, typically explicitly added,
! along with geometrical source additions
call addsource(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
   ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iwmin,iwmax,qtC,wCT,qt,w,x,&
   qsourcesplit)

end subroutine addsource2
!=============================================================================
subroutine process(iit,qt)

use mod_amrvacdef

integer,intent(in)            :: iit
double precision, intent(in)  :: qt
! .. local ..
integer:: iigrid, igrid, level
!-----------------------------------------------------------------------------

!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid,level)
do iigrid=1,igridstail; igrid=igrids(iigrid);

   call set_tmpGlobals(igrid)
   level = node(plevel_,igrid)
   
   call process_grid_usr(igrid,level,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
      ixGhi3,ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,qt,pw(igrid)%w,&
      px(igrid)%x)
   
end do
!$OMP END PARALLEL DO

!opedit: Id like to synchronize in case the user does change the state:
call getbc(qt,ps,psCoarse)

end subroutine process
!=============================================================================
subroutine addsource_impl

use mod_amrvacdef

integer :: iigrid, igrid, icycle, ncycle
double precision :: qdt, qt, sumqdt
!-----------------------------------------------------------------------------

energyonly=.true.
if(sourceimplcycle)then
   ncycle=ceiling(dt/dtimpl)
   if(ncycle<1) then
     ncycle=1
     dtimpl=dt
   endif
else
   ncycle=1
endif

if(mype==0.and..false.) then
  print *,'implicit source addition will subcycle with ',ncycle,' subtimesteps'
  print *,'dt and dtimpl= ',dt,dtimpl,' versus ncycle*dtimpl=',ncycle*dtimpl
endif

qt=t
sumqdt=zero
do icycle=1,ncycle
  if(sourceparasts)then
    qdt=dtimpl/(one+parastsnu+(parastsnu-one)* dcos(half*dpi*(two*dble&
       (icycle)-one)/dble(ncycle)))
  else
    qdt=dt/dble(ncycle)
  endif
!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
  do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
   call addsourceimpl1_grid(igrid,qdt,qt,pw(igrid)%w,px(igrid)%x)
  end do
!$OMP END PARALLEL DO

  sumqdt=sumqdt+qdt
  qt=qt+qdt
  call getbc(qt,ps,psCoarse)
enddo

if(mype==0.and..false.) then
  if(sourceparasts)then
     print *,'qt-t=',qt-t,'versus ncycle2*dtimpl=',ncycle*ncycle*dtimpl,' and sumqdt=',sumqdt
  else
     print *,'qt-t=',qt-t,'versus ncycle*dtimpl=',ncycle*dtimpl,' and sumqdt=',sumqdt
  endif
endif
energyonly=.false.

end subroutine addsource_impl
!=============================================================================
subroutine addsourceimpl1_grid(igrid,qdt,qt,w,x)

use mod_amrvacdef

integer, intent(in) :: igrid
double precision, intent(in) :: qdt, qt, x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
   ixGlo3:ixGhi3,1:ndim)
double precision, intent(inout) :: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
   ixGlo3:ixGhi3,nw)

double precision :: w1(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw)
!-----------------------------------------------------------------------------

call set_tmpGlobals(igrid)

w1(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)=w(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)

call specialsource_impl(qdt,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixMlo1,&
   ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,1,nw,qt,w1,qt,w,x)

end subroutine addsourceimpl1_grid
!=============================================================================
