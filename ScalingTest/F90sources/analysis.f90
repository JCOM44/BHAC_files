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
subroutine write_analysis
! This is an example file how to use the analysis capability.  You can schedule 
! this routine using the slot 5 in itsave, dtsave and ditsave.  
! To use, just copy this file to your working directory and make your modifications.
use constants
use mod_amrvacdef
character(len=20):: userconvert_type
!-----------------------------------------------------------------------------
logical :: fileopen
integer :: iigrid, igrid, level, iw
double precision :: volume(1:nlevelshi),  voltotal
double precision :: re, ke, me, te
double precision :: dvolume(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
    inBubble(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3), lfac_pw&
   (ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3), volumeflat(1:nlevelshi)
integer :: numlevels
integer, dimension(1:nlevelshi) :: isum_send
double precision, dimension(1:4) :: dsum_send, dsum_recv
character(len=80) :: filename
character(len=1024) :: line
logical, save :: opened=.false.
logical, save :: file_exists=.false.
integer :: amode, status(MPI_STATUS_SIZE)
double precision :: trcut
!-----------------------------------------------------------------------------

!!$!!!Selects the bubble volume: !!!!!!!!!!
!!$trcut = 1.0d-3
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$volume(1:mxnest)=zero
!!$volumeflat(1:mxnest)=zero
!!$re = 0.0d0
!!$ke = 0.0d0
!!$me = 0.0d0
!!$te = 0.0d0
!!$
!!$
!!$do iigrid=1,igridstail; igrid=igrids(iigrid);
!!$
!!$call primitive(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x)
!!$! First normalize to cgs units:
!!$   {^D&px(igrid)%x(ixM^T,^D) = normvar(0) * px(igrid)%x(ixM^T,^D)\}
!!$   do iw = 1, nw
!!$      pw(igrid)%w(ixM^T,iw) = normvar(iw) * pw(igrid)%w(ixM^T,iw)
!!$   end do
!!$! Mask out bubble volume:
!!$   call get_lfac(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,lfac_pw)
!!$   where (pw(igrid)%w(ixM^T,lfac_) .le. eqpar(flfac_)*lfac_pw(ixM^T) .and. &
!!$       pw(igrid)%w(ixM^T,tr1_) .ge. trcut )
!!$      inBubble(ixM^T) = one
!!$   elsewhere
!!$      inBubble(ixM^T) = zero
!!$   end where
!!$
!!$   level=node(plevel_,igrid)
!!$   volumeflat(level)=volumeflat(level)+ &
!!$          {(rnode(rpxmax^D_,igrid)-rnode(rpxmin^D_,igrid))|*} &
!!$          * normvar(0)**3.0d0
!!$   if (slab) then
!!$      dvolume(ixM^T)={rnode(rpdx^D_,igrid)|*} &
!!$           * normvar(0)**3.0d0
!!$   else
!!$      dvolume(ixM^T)=pgeo(igrid)%dvolume(ixM^T) &
!!$           * normvar(0)**3.0d0 * 2.0d0* dpi
!!$      volume(level)=volume(level)+sum(dvolume(ixM^T))
!!$   end if
!!$! Now compute the energies:
!!$
!!$   re = re + sum(pw(igrid)%w(ixM^T,lfac_)*pw(igrid)%w(ixM^T,rho_)*CONST_c**2.0d0 & 
!!$        * dvolume(ixM^T) * inBubble(ixM^T))
!!$   ke = ke + sum((pw(igrid)%w(ixM^T,lfac_) - 1.0d0) &
!!$        * pw(igrid)%w(ixM^T,lfac_)*pw(igrid)%w(ixM^T,rho_)*CONST_c**2.0d0 &
!!$        * dvolume(ixM^T) * inBubble(ixM^T))
!!$   me = me + sum((({^C& pw(igrid)%w(ixM^T,b^C_)**2.0d0 |+})/4.0d0/dpi &
!!$        - one/8.0d0/dpi/pw(igrid)%w(ixM^T,lfac_)**2.0d0 * &
!!$        (({^C& pw(igrid)%w(ixM^T,b^C_)**2.0d0 |+}) &
!!$        + ({^C& pw(igrid)%w(ixM^T,u^C_)/CONST_c * pw(igrid)%w(ixM^T,b^C_) |+})**2.0d0 )) &
!!$        * dvolume(ixM^T) * inBubble(ixM^T))
!!$   te = te + sum((4.0d0*pw(igrid)%w(ixM^T,lfac_)**2.0d0*pw(igrid)%w(ixM^T,pp_)  &
!!$        - pw(igrid)%w(ixM^T,pp_) ) &
!!$        * dvolume(ixM^T) * inBubble(ixM^T))
!!$
!!$end do
!!$if (slab) volume(levmin:levmax)=volumeflat(levmin:levmax)
!!$
!!$voltotal=sum(volume(levmin:levmax))
!!$
!!$
!!$dsum_send(1)=re
!!$dsum_send(2)=ke
!!$dsum_send(3)=me
!!$dsum_send(4)=te
!!$
!!$call MPI_REDUCE(dsum_send,dsum_recv,4,MPI_DOUBLE_PRECISION, &
!!$                MPI_SUM,0,icomm,ierrmpi)


!if (mype==0) then

!   re = dsum_recv(1)
!   ke = dsum_recv(2)
!   me = dsum_recv(3)
!   te = dsum_recv(4)

!   if (.not.opened) then
!      ! generate filename
!      write(filename,"(a,a)") TRIM("analysis"),".csv"
!      INQUIRE(FILE=filename, EXIST=file_exists)

!      if (.not. file_exists) then
!         open(unit=unitanalysis,file=filename,status='unknown',access='append')
!         write(unitanalysis,"(a)") trim('# t [years] re [erg] ke [erg] me [erg] te [erg]')
!      else
!         open(unit=unitanalysis,file=filename,status='unknown',access='append')
!      end if
!      opened=.true.
!   end if

!   write(unitanalysis,'(5(es14.6))')t*UNIT_LENGTH/UNIT_VELOCITY/CONST_years,re,ke,me,te

!end if

end subroutine write_analysis
!=============================================================================
