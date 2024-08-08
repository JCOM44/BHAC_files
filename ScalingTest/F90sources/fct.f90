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
subroutine recalculateB

use mod_amrvacdef


! ... local ...
integer :: igrid,iigrid
!-----------------------------------------------------------------------------

call init_comm_fix_conserve(1,3)

!$OMP PARALLEL DO PRIVATE(igrid)
do iigrid=1,igridstail; igrid=igrids(iigrid);
   ! Make zero the magnetic fluxes
   ! Fake advance, storing electric fields at edges
   call fake_advance(igrid,1,3,ps(igrid))

end do
!$OMP END PARALLEL DO

! Do correction

do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
   call sendflux(igrid,1,3) ! OMP: not threadsafe!
end do

call fix_conserve(ps,1,3)

call fix_edges(ps,1,3)

! Now we fill the centers for the staggered variables
!$OMP PARALLEL DO PRIVATE(igrid)
do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
   call faces2centers(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ps(igrid))
end do
!$OMP END PARALLEL DO

end subroutine recalculateB
!=============================================================================
subroutine fake_advance(igrid,idimmin,idimmax,s)

use mod_amrvacdef

integer       :: igrid,idimmin,idimmax
type(state)   :: s
! ... local ...
double precision             :: dx1,dx2,dx3
double precision             :: fC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
   1:nwflux,1:ndim)
double precision             :: fE(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
   1:3)
!-----------------------------------------------------------------------------

dx1=rnode(rpdx1_,igrid);dx2=rnode(rpdx2_,igrid);dx3=rnode(rpdx3_,igrid);
call set_tmpGlobals(igrid)

call fake_update(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,s,fC,fE,dx1,dx2,&
   dx3)

call storeflux(igrid,fC,idimmin,idimmax)
call storeedge(igrid,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,fE,idimmin,&
   idimmax) 

end subroutine fake_advance
!=============================================================================
subroutine fake_update(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,s,fC,&
   fE,dx1,dx2,dx3)
! In reality the same as b_from_vectorpotential for staggered case

use mod_amrvacdef

integer       :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3
type(state)   :: s
double precision             :: fC(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nwflux,1:ndim)
double precision             :: fE(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
   1:3)
double precision             :: dx1,dx2,dx3
! ... local ...
integer                            :: ixIsmin1,ixIsmin2,ixIsmin3,ixIsmax1,&
   ixIsmax2,ixIsmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idir
double precision                   :: xC(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndim), A(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:ndir)
double precision                   :: circ(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndim), dxidir
!-----------------------------------------------------------------------------
associate(ws=>s%ws%w,x=>s%x%x)

A(:,:,:,:)=zero
ws(:,:,:,:)=zero

ixIsmin1=s%ws%ixGmin1;ixIsmin2=s%ws%ixGmin2;ixIsmin3=s%ws%ixGmin3
ixIsmax1=s%ws%ixGmax1;ixIsmax2=s%ws%ixGmax2;ixIsmax3=s%ws%ixGmax3;
ixOmin1=ixImin1+dixB;ixOmin2=ixImin2+dixB;ixOmin3=ixImin3+dixB
ixOmax1=ixImax1-dixB;ixOmax2=ixImax2-dixB;ixOmax3=ixImax3-dixB;

call b_from_vectorpotentialA(ixIsmin1,ixIsmin2,ixIsmin3,ixIsmax1,ixIsmax2,&
   ixIsmax3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,&
   ixOmin3,ixOmax1,ixOmax2,ixOmax3, ws, x, A)

! This is important only in 3D

do idir=1,ndim
   fE(:,:,:,idir) =-A(:,:,:,idir)*dxlevel(idir)
end do

end associate
end subroutine fake_update
!=============================================================================

!=============================================================================

subroutine faces2centers(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,s)

  use mod_amrvacdef

  ! Non-staggered interpolation range
  integer, intent(in)                :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
     ixOmax2,ixOmax3
  type(state)                        :: s

  ! --- local ---
  integer                            :: hxOmin1,hxOmin2,hxOmin3,hxOmax1,&
     hxOmax2,hxOmax3, idim
 !-----------------------------------------------------------------------------
  associate(w=>s%w%w, ws=>s%ws%w, mygeo=>s%geo)

  do idim=1,ndim
     ! Displace index to the left
 !Even if ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3 is the full size of the w arrays, this is ok
     ! because the staggered arrays have an additional place to the left.

     hxOmin1=ixOmin1-kr(idim,1);hxOmin2=ixOmin2-kr(idim,2)
     hxOmin3=ixOmin3-kr(idim,3);hxOmax1=ixOmax1-kr(idim,1)
     hxOmax2=ixOmax2-kr(idim,2);hxOmax3=ixOmax3-kr(idim,3);

     ! Interpolate to cell barycentre using arithmetic average
     ! This might be done better later, to make the method less diffusive.
     select case(idim)
        case(1)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+idim)&
           =(half*mygeo%dx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idim)/mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))*&
             (ws(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                bs0_+idim)*mygeo%surfaceC1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)&
             +ws(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,&
                bs0_+idim)*mygeo%surfaceC1(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                hxOmin3:hxOmax3))
        
case(2)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+idim)&
           =(half*mygeo%dx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idim)/mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))*&
             (ws(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                bs0_+idim)*mygeo%surfaceC2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)&
             +ws(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,&
                bs0_+idim)*mygeo%surfaceC2(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                hxOmin3:hxOmax3))
        
case(3)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+idim)&
           =(half*mygeo%dx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idim)/mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))*&
             (ws(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                bs0_+idim)*mygeo%surfaceC3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)&
             +ws(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,&
                bs0_+idim)*mygeo%surfaceC3(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                hxOmin3:hxOmax3))
        
     end select
  end do

  end associate

end subroutine faces2centers
!=============================================================================
subroutine faces2centers4(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,s)

  use mod_amrvacdef

  ! Non-staggered interpolation range
  ! Fourth order equation.
  integer, intent(in)                :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
     ixOmax2,ixOmax3
  type(state)                        :: s

  ! --- local ---
  integer                            :: gxOmin1,gxOmin2,gxOmin3,gxOmax1,&
     gxOmax2,gxOmax3, hxOmin1,hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3,&
      jxOmin1,jxOmin2,jxOmin3,jxOmax1,jxOmax2,jxOmax3, idim
 !-----------------------------------------------------------------------------
  associate(w=>s%w%w, ws=>s%ws%w, mygeo=>s%geo)

  do idim=1,ndim

     gxOmin1=ixOmin1-2*kr(idim,1);gxOmin2=ixOmin2-2*kr(idim,2)
     gxOmin3=ixOmin3-2*kr(idim,3);gxOmax1=ixOmax1-2*kr(idim,1)
     gxOmax2=ixOmax2-2*kr(idim,2);gxOmax3=ixOmax3-2*kr(idim,3);
     hxOmin1=ixOmin1-kr(idim,1);hxOmin2=ixOmin2-kr(idim,2)
     hxOmin3=ixOmin3-kr(idim,3);hxOmax1=ixOmax1-kr(idim,1)
     hxOmax2=ixOmax2-kr(idim,2);hxOmax3=ixOmax3-kr(idim,3);
     jxOmin1=ixOmin1+kr(idim,1);jxOmin2=ixOmin2+kr(idim,2)
     jxOmin3=ixOmin3+kr(idim,3);jxOmax1=ixOmax1+kr(idim,1)
     jxOmax2=ixOmax2+kr(idim,2);jxOmax3=ixOmax3+kr(idim,3);

     ! Interpolate to cell barycentre using fourth order central formula
     select case(idim)
        case(1)     
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+idim)=(1.0d0&
           /16.0d0/mygeo%surface1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)) &
             *( &
             -ws(gxOmin1:gxOmax1,gxOmin2:gxOmax2,gxOmin3:gxOmax3,&
                bs0_+idim)*mygeo%surfaceC1(gxOmin1:gxOmax1,gxOmin2:gxOmax2,&
                gxOmin3:gxOmax3) &
             +9.0d0*ws(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,&
                bs0_+idim)*mygeo%surfaceC1(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                hxOmin3:hxOmax3) &
             +9.0d0*ws(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                bs0_+idim)*mygeo%surfaceC1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3) &
             -ws(jxOmin1:jxOmax1,jxOmin2:jxOmax2,jxOmin3:jxOmax3,&
                bs0_+idim)*mygeo%surfaceC1(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                jxOmin3:jxOmax3) &
             )
        
case(2)     
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+idim)=(1.0d0&
           /16.0d0/mygeo%surface2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)) &
             *( &
             -ws(gxOmin1:gxOmax1,gxOmin2:gxOmax2,gxOmin3:gxOmax3,&
                bs0_+idim)*mygeo%surfaceC2(gxOmin1:gxOmax1,gxOmin2:gxOmax2,&
                gxOmin3:gxOmax3) &
             +9.0d0*ws(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,&
                bs0_+idim)*mygeo%surfaceC2(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                hxOmin3:hxOmax3) &
             +9.0d0*ws(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                bs0_+idim)*mygeo%surfaceC2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3) &
             -ws(jxOmin1:jxOmax1,jxOmin2:jxOmax2,jxOmin3:jxOmax3,&
                bs0_+idim)*mygeo%surfaceC2(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                jxOmin3:jxOmax3) &
             )
        
case(3)     
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+idim)=(1.0d0&
           /16.0d0/mygeo%surface3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)) &
             *( &
             -ws(gxOmin1:gxOmax1,gxOmin2:gxOmax2,gxOmin3:gxOmax3,&
                bs0_+idim)*mygeo%surfaceC3(gxOmin1:gxOmax1,gxOmin2:gxOmax2,&
                gxOmin3:gxOmax3) &
             +9.0d0*ws(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,&
                bs0_+idim)*mygeo%surfaceC3(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                hxOmin3:hxOmax3) &
             +9.0d0*ws(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                bs0_+idim)*mygeo%surfaceC3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3) &
             -ws(jxOmin1:jxOmax1,jxOmin2:jxOmax2,jxOmin3:jxOmax3,&
                bs0_+idim)*mygeo%surfaceC3(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                jxOmin3:jxOmax3) &
             )
        
     end select
  end do

  end associate

end subroutine faces2centers4
!=============================================================================
subroutine faces2centers6(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,s)

  use mod_amrvacdef

  ! Non-staggered interpolation range
  ! Sixth order equation.
  integer, intent(in)                :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
     ixOmax2,ixOmax3
  type(state)                        :: s

  ! --- local ---
  integer                            :: fxOmin1,fxOmin2,fxOmin3,fxOmax1,&
     fxOmax2,fxOmax3, gxOmin1,gxOmin2,gxOmin3,gxOmax1,gxOmax2,gxOmax3,&
      hxOmin1,hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3, jxOmin1,jxOmin2,&
     jxOmin3,jxOmax1,jxOmax2,jxOmax3, kxOmin1,kxOmin2,kxOmin3,kxOmax1,kxOmax2,&
     kxOmax3, idim
 !-----------------------------------------------------------------------------
  associate(w=>s%w%w, ws=>s%ws%w, mygeo=>s%geo)

  do idim=1,ndim

     fxOmin1=ixOmin1-3*kr(idim,1);fxOmin2=ixOmin2-3*kr(idim,2)
     fxOmin3=ixOmin3-3*kr(idim,3);fxOmax1=ixOmax1-3*kr(idim,1)
     fxOmax2=ixOmax2-3*kr(idim,2);fxOmax3=ixOmax3-3*kr(idim,3);
     gxOmin1=ixOmin1-2*kr(idim,1);gxOmin2=ixOmin2-2*kr(idim,2)
     gxOmin3=ixOmin3-2*kr(idim,3);gxOmax1=ixOmax1-2*kr(idim,1)
     gxOmax2=ixOmax2-2*kr(idim,2);gxOmax3=ixOmax3-2*kr(idim,3);
     hxOmin1=ixOmin1-kr(idim,1);hxOmin2=ixOmin2-kr(idim,2)
     hxOmin3=ixOmin3-kr(idim,3);hxOmax1=ixOmax1-kr(idim,1)
     hxOmax2=ixOmax2-kr(idim,2);hxOmax3=ixOmax3-kr(idim,3);
     jxOmin1=ixOmin1+kr(idim,1);jxOmin2=ixOmin2+kr(idim,2)
     jxOmin3=ixOmin3+kr(idim,3);jxOmax1=ixOmax1+kr(idim,1)
     jxOmax2=ixOmax2+kr(idim,2);jxOmax3=ixOmax3+kr(idim,3);
     kxOmin1=ixOmin1+2*kr(idim,1);kxOmin2=ixOmin2+2*kr(idim,2)
     kxOmin3=ixOmin3+2*kr(idim,3);kxOmax1=ixOmax1+2*kr(idim,1)
     kxOmax2=ixOmax2+2*kr(idim,2);kxOmax3=ixOmax3+2*kr(idim,3);

     ! Interpolate to cell barycentre using sixth order central formula
     select case(idim)
        case(1)     
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+idim)=(1.0d0&
           /256.0d0/mygeo%surface1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)) &
             *( &
             +3.0d0*ws(fxOmin1:fxOmax1,fxOmin2:fxOmax2,fxOmin3:fxOmax3,&
                bs0_+idim)*mygeo%surfaceC1(fxOmin1:fxOmax1,fxOmin2:fxOmax2,&
                fxOmin3:fxOmax3) &
             -25.0d0*ws(gxOmin1:gxOmax1,gxOmin2:gxOmax2,gxOmin3:gxOmax3,&
                bs0_+idim)*mygeo%surfaceC1(gxOmin1:gxOmax1,gxOmin2:gxOmax2,&
                gxOmin3:gxOmax3) &
             +150.0d0*ws(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,&
                bs0_+idim)*mygeo%surfaceC1(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                hxOmin3:hxOmax3) &
             +150.0d0*ws(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                bs0_+idim)*mygeo%surfaceC1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3) &
             -25.0d0*ws(jxOmin1:jxOmax1,jxOmin2:jxOmax2,jxOmin3:jxOmax3,&
                bs0_+idim)*mygeo%surfaceC1(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                jxOmin3:jxOmax3) &
             +3.0d0*ws(kxOmin1:kxOmax1,kxOmin2:kxOmax2,kxOmin3:kxOmax3,&
                bs0_+idim)*mygeo%surfaceC1(kxOmin1:kxOmax1,kxOmin2:kxOmax2,&
                kxOmin3:kxOmax3) &
             )
        
case(2)     
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+idim)=(1.0d0&
           /256.0d0/mygeo%surface2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)) &
             *( &
             +3.0d0*ws(fxOmin1:fxOmax1,fxOmin2:fxOmax2,fxOmin3:fxOmax3,&
                bs0_+idim)*mygeo%surfaceC2(fxOmin1:fxOmax1,fxOmin2:fxOmax2,&
                fxOmin3:fxOmax3) &
             -25.0d0*ws(gxOmin1:gxOmax1,gxOmin2:gxOmax2,gxOmin3:gxOmax3,&
                bs0_+idim)*mygeo%surfaceC2(gxOmin1:gxOmax1,gxOmin2:gxOmax2,&
                gxOmin3:gxOmax3) &
             +150.0d0*ws(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,&
                bs0_+idim)*mygeo%surfaceC2(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                hxOmin3:hxOmax3) &
             +150.0d0*ws(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                bs0_+idim)*mygeo%surfaceC2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3) &
             -25.0d0*ws(jxOmin1:jxOmax1,jxOmin2:jxOmax2,jxOmin3:jxOmax3,&
                bs0_+idim)*mygeo%surfaceC2(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                jxOmin3:jxOmax3) &
             +3.0d0*ws(kxOmin1:kxOmax1,kxOmin2:kxOmax2,kxOmin3:kxOmax3,&
                bs0_+idim)*mygeo%surfaceC2(kxOmin1:kxOmax1,kxOmin2:kxOmax2,&
                kxOmin3:kxOmax3) &
             )
        
case(3)     
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+idim)=(1.0d0&
           /256.0d0/mygeo%surface3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)) &
             *( &
             +3.0d0*ws(fxOmin1:fxOmax1,fxOmin2:fxOmax2,fxOmin3:fxOmax3,&
                bs0_+idim)*mygeo%surfaceC3(fxOmin1:fxOmax1,fxOmin2:fxOmax2,&
                fxOmin3:fxOmax3) &
             -25.0d0*ws(gxOmin1:gxOmax1,gxOmin2:gxOmax2,gxOmin3:gxOmax3,&
                bs0_+idim)*mygeo%surfaceC3(gxOmin1:gxOmax1,gxOmin2:gxOmax2,&
                gxOmin3:gxOmax3) &
             +150.0d0*ws(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,&
                bs0_+idim)*mygeo%surfaceC3(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                hxOmin3:hxOmax3) &
             +150.0d0*ws(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                bs0_+idim)*mygeo%surfaceC3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3) &
             -25.0d0*ws(jxOmin1:jxOmax1,jxOmin2:jxOmax2,jxOmin3:jxOmax3,&
                bs0_+idim)*mygeo%surfaceC3(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                jxOmin3:jxOmax3) &
             +3.0d0*ws(kxOmin1:kxOmax1,kxOmin2:kxOmax2,kxOmin3:kxOmax3,&
                bs0_+idim)*mygeo%surfaceC3(kxOmin1:kxOmax1,kxOmin2:kxOmax2,&
                kxOmin3:kxOmax3) &
             )
        
     end select
  end do

  end associate

end subroutine faces2centers6
!=============================================================================
subroutine show_div_staggered(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
   s)
! For calculating the divergence of a face-allocated vector field.
use mod_amrvacdef

integer, intent(in)           :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
   ixOmax3
type(state)                   :: s
! ... local ...
integer                       :: hxOmin1,hxOmin2,hxOmin3,hxOmax1,hxOmax2,&
   hxOmax3,idim
double precision              :: out(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
   ixGlo3:ixGhi3)
!----------------------------------------------------------------------------
print *, 'In div_stg'

print *, 'x1 min = ', s%x%x(ixMlo1,ixMlo2,ixMlo3,1)
print *, 'x2 min = ', s%x%x(ixMlo1,ixMlo2,ixMlo3,2)

associate(w=>s%w%w, ws=>s%ws%w, mygeo=>s%geo)
out=zero

do idim=1,ndim
   ! Displace index to the left
 !Even if ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3 is the full size of the w arrays, this is ok
   ! because the staggered arrays have an additional place to the left.

   hxOmin1=ixOmin1-kr(idim,1);hxOmin2=ixOmin2-kr(idim,2)
   hxOmin3=ixOmin3-kr(idim,3);hxOmax1=ixOmax1-kr(idim,1)
   hxOmax2=ixOmax2-kr(idim,2);hxOmax3=ixOmax3-kr(idim,3);
 !print *, 'ixOmin1','ixOmin2','ixOmin3','ixOmax1','ixOmax2','ixOmax3', ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
 !print *, 'hxOmin1','hxOmin2','hxOmin3','hxOmax1','hxOmax2','hxOmax3', hxOmin1,hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3
   ! Calculate divergence by taking differences of fluxes
   select case(idim)
   case(1)
   out(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=out(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3)+(ws(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      ixOmin3:ixOmax3,bs0_+idim)*mygeo%surfaceC1(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                        -ws(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,&
                           bs0_+idim)*mygeo%surfaceC1(hxOmin1:hxOmax1,&
                           hxOmin2:hxOmax2,hxOmin3:hxOmax3))&
                        /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                           ixOmin3:ixOmax3)
   
case(2)
   out(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=out(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3)+(ws(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      ixOmin3:ixOmax3,bs0_+idim)*mygeo%surfaceC2(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                        -ws(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,&
                           bs0_+idim)*mygeo%surfaceC2(hxOmin1:hxOmax1,&
                           hxOmin2:hxOmax2,hxOmin3:hxOmax3))&
                        /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                           ixOmin3:ixOmax3)
   
case(3)
   out(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=out(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3)+(ws(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      ixOmin3:ixOmax3,bs0_+idim)*mygeo%surfaceC3(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                        -ws(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,&
                           bs0_+idim)*mygeo%surfaceC3(hxOmin1:hxOmax1,&
                           hxOmin2:hxOmax2,hxOmin3:hxOmax3))&
                        /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                           ixOmin3:ixOmax3)
   
   end select
end do
end associate

call printarray(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,out)

end subroutine show_div_staggered
!---------------------------------------------------------------------------
subroutine div_staggered(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,s,&
   out)
! For calculating the divergence of a face-allocated vector field.
use mod_amrvacdef

integer, intent(in)           :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
   ixOmax3
type(state)                   :: s
double precision              :: out(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   ixOmin3:ixOmax3)
! ... local ...
integer                       :: hxOmin1,hxOmin2,hxOmin3,hxOmax1,hxOmax2,&
   hxOmax3,idim
!----------------------------------------------------------------------------

associate(w=>s%w%w, ws=>s%ws%w, mygeo=>s%geo)
out(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=zero

do idim=1,ndim
   ! Displace index to the left
 !Even if ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3 is the full size of the w arrays, this is ok
   ! because the staggered arrays have an additional place to the left.

   hxOmin1=ixOmin1-kr(idim,1);hxOmin2=ixOmin2-kr(idim,2)
   hxOmin3=ixOmin3-kr(idim,3);hxOmax1=ixOmax1-kr(idim,1)
   hxOmax2=ixOmax2-kr(idim,2);hxOmax3=ixOmax3-kr(idim,3);
 !print *, 'ixOmin1','ixOmin2','ixOmin3','ixOmax1','ixOmax2','ixOmax3', ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
 !print *, 'hxOmin1','hxOmin2','hxOmin3','hxOmax1','hxOmax2','hxOmax3', hxOmin1,hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3
   ! Calculate divergence by taking differences of fluxes
   select case(idim)
   case(1)
   out(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=out(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3)+(ws(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      ixOmin3:ixOmax3,bs0_+idim)*mygeo%surfaceC1(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                        -ws(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,&
                           bs0_+idim)*mygeo%surfaceC1(hxOmin1:hxOmax1,&
                           hxOmin2:hxOmax2,hxOmin3:hxOmax3))&
                        /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                           ixOmin3:ixOmax3)
   
case(2)
   out(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=out(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3)+(ws(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      ixOmin3:ixOmax3,bs0_+idim)*mygeo%surfaceC2(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                        -ws(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,&
                           bs0_+idim)*mygeo%surfaceC2(hxOmin1:hxOmax1,&
                           hxOmin2:hxOmax2,hxOmin3:hxOmax3))&
                        /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                           ixOmin3:ixOmax3)
   
case(3)
   out(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=out(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3)+(ws(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      ixOmin3:ixOmax3,bs0_+idim)*mygeo%surfaceC3(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                        -ws(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,&
                           bs0_+idim)*mygeo%surfaceC3(hxOmin1:hxOmax1,&
                           hxOmin2:hxOmax2,hxOmin3:hxOmax3))&
                        /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                           ixOmin3:ixOmax3)
   
   end select
end do
end associate

 !call printarray(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,out)

end subroutine div_staggered
!=============================================================================
subroutine updatefaces(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qdt,fC,fE,s)

use mod_amrvacdef

integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(in)       :: qdt
type(state)                        :: s
 !double precision, intent(in)       :: w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw)
 !double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
double precision, intent(in)       :: fC(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nwflux,1:ndim)
double precision, intent(inout)    :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:3)
 !double precision, intent(inout)    :: bfaces(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,1:ndir)

! --- local ---
integer                            :: ixmin1(1:ndim),ixmin2(1:ndim),&
   ixmin3(1:ndim),ixmax1(1:ndim),ixmax2(1:ndim),ixmax3(1:ndim)
integer                            :: hxCmin1,hxCmin2,hxCmin3,hxCmax1,hxCmax2,&
   hxCmax3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,ixCpmin1,ixCpmin2,&
   ixCpmin3,ixCpmax1,ixCpmax2,ixCpmax3,jxCmin1,jxCmin2,jxCmin3,jxCmax1,&
   jxCmax2,jxCmax3,ixCmmin1,ixCmmin2,ixCmmin3,ixCmmax1,ixCmmax2,ixCmmax3
integer                            :: idim1,idim2,idir,iwdim1,iwdim2,i,j,k
 !double precision                   :: edge(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndir)
double precision                   :: circ(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndim)
!-----------------------------------------------------------------------------
!associate(mygeo=>s%geo,bfaces=>s%ws%w)
associate(bfaces=>s%ws%w,x=>s%x%x)

! Calculate contribution to FEM of each edge,
! that is, estimate value of line integral of
! electric field in the positive idir direction.

ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
ixCmin1=ixOmin1-1;ixCmin2=ixOmin2-1;ixCmin3=ixOmin3-1;

fE(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndir)=zero

do idim1=1,ndim 
   iwdim1 = b0_+idim1
   do idim2=1,ndim
      iwdim2 = b0_+idim2
      do idir=1,ndir ! Direction of line integral
         ! Allow only even permutations
         if (lvc(idim1,idim2,idir).eq.1) then
            ! Assemble indices
            jxCmin1=ixCmin1+kr(idim1,1);jxCmin2=ixCmin2+kr(idim1,2)
            jxCmin3=ixCmin3+kr(idim1,3);jxCmax1=ixCmax1+kr(idim1,1)
            jxCmax2=ixCmax2+kr(idim1,2);jxCmax3=ixCmax3+kr(idim1,3);
            ixCpmin1=ixCmin1+kr(idim2,1);ixCpmin2=ixCmin2+kr(idim2,2)
            ixCpmin3=ixCmin3+kr(idim2,3);ixCpmax1=ixCmax1+kr(idim2,1)
            ixCpmax2=ixCmax2+kr(idim2,2);ixCpmax3=ixCmax3+kr(idim2,3);
            ! Interpolate to edges
            fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)&
               =qdt*quarter*((fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,iwdim1,idim2)&
            +fC(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,iwdim1,&
               idim2))/dxlevel(idim1)&
            -(fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iwdim2,&
               idim1)+fC(ixCpmin1:ixCpmax1,ixCpmin2:ixCpmax2,&
               ixCpmin3:ixCpmax3,iwdim2,idim1))&
            /dxlevel(idim2))

            if (typeaxial.ne.'slab') then


            ! Catch axis and origin, where the line integral is performed on a
            ! zero lenght element and should be zero.
            ! Remember that staggered quantities are assigned to the index
            ! of the cell at their left.
            where((abs(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               z_)+half*dxlevel(z_)).lt.1.0d-9).or.&
                  (abs(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                     z_)+half*dxlevel(z_)-dpi).lt.1.0d-9))

              fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2)=zero
              fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,r_)=zero
            end where
            where(abs(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               r_)+half*dxlevel(r_)).lt.1.0d-9)
              fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)=zero
            end where
            end if
         end if
      end do
   end do
end do

circ(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)=zero

! Calculate circulation on each face

do idim1=1,ndim ! Coordinate perpendicular to face 
   do idim2=1,ndim
      do idir=1,ndir ! Direction of line integral
        ! Assemble indices
        hxCmin1=ixCmin1-kr(idim2,1);hxCmin2=ixCmin2-kr(idim2,2)
        hxCmin3=ixCmin3-kr(idim2,3);hxCmax1=ixCmax1-kr(idim2,1)
        hxCmax2=ixCmax2-kr(idim2,2);hxCmax3=ixCmax3-kr(idim2,3);
        ! Add line integrals in direction idir
        circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
           =circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
                         +lvc(idim1,idim2,idir)&
                         *(fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                            idir)&
                         -fE(hxCmin1:hxCmax1,hxCmin2:hxCmax2,hxCmin3:hxCmax3,&
                            idir))
      end do
   end do
end do

! Decrease bottom limit by one
do idim1=1, ndim
   ixmax1(idim1)=ixOmax1;ixmax2(idim1)=ixOmax2;ixmax3(idim1)=ixOmax3
   ixmin1=ixOmin1-1;ixmin2=ixOmin2-1;ixmin3=ixOmin3-1; !-kr(1,2,3,idim1);
end do

! Divide by the area of the face to get dB/dt

do idim1=1,ndim
   ixCmin1=ixmin1(idim1);ixCmin2=ixmin2(idim1);ixCmin3=ixmin3(idim1)
   ixCmax1=ixmax1(idim1);ixCmax2=ixmax2(idim1);ixCmax3=ixmax3(idim1);
   select case(idim1)
   case(1)
   where(mygeo%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
      .gt. 1.0d-9*mygeo%dvolume(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      ixCmin3:ixCmax3))
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
         =circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
                 /mygeo%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                    ixCmin3:ixCmax3)
   elsewhere
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)=zero
   end where
   
case(2)
   where(mygeo%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
      .gt. 1.0d-9*mygeo%dvolume(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      ixCmin3:ixCmax3))
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
         =circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
                 /mygeo%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                    ixCmin3:ixCmax3)
   elsewhere
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)=zero
   end where
   
case(3)
   where(mygeo%surfaceC3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
      .gt. 1.0d-9*mygeo%dvolume(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      ixCmin3:ixCmax3))
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
         =circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
                 /mygeo%surfaceC3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                    ixCmin3:ixCmax3)
   elsewhere
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)=zero
   end where
   
    end select

   ! Time update
   ! minus!
   bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
      =bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
      idim1)-circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)
end do

end associate

end subroutine updatefaces
!=============================================================================
subroutine updatefacesuct2(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qdt,vbarC,cbarmin,cbarmax,&
   fE,s)


  
use mod_amrvacdef

integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(in)       :: qdt
type(state)                        :: s
double precision, intent(in)       :: vbarC(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndir,2)
double precision, intent(in)       :: cbarmin(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,ndim)
double precision, intent(in)       :: cbarmax(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,ndim)
double precision, intent(inout)    :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndir)
! --- local ---
double precision                   :: vtilL(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,2)
double precision                   :: vtilR(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,2)
double precision                   :: btilL(s%ws%ixGmin1:s%ws%ixGmax1,&
   s%ws%ixGmin2:s%ws%ixGmax2,s%ws%ixGmin3:s%ws%ixGmax3,ndim)
double precision                   :: btilR(s%ws%ixGmin1:s%ws%ixGmax1,&
   s%ws%ixGmin2:s%ws%ixGmax2,s%ws%ixGmin3:s%ws%ixGmax3,ndim)
double precision                   :: sqrtg(s%ws%ixGmin1:s%ws%ixGmax1,&
   s%ws%ixGmin2:s%ws%ixGmax2,s%ws%ixGmin3:s%ws%ixGmax3)
double precision                   :: cp(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,2)
double precision                   :: cm(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,2)
integer                            :: ixGsmin1,ixGsmin2,ixGsmin3,ixGsmax1,&
   ixGsmax2,ixGsmax3
integer                            :: hxCmin1,hxCmin2,hxCmin3,hxCmax1,hxCmax2,&
   hxCmax3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,ixCpmin1,ixCpmin2,&
   ixCpmin3,ixCpmax1,ixCpmax2,ixCpmax3,jxCmin1,jxCmin2,jxCmin3,jxCmax1,&
   jxCmax2,jxCmax3,ixCmmin1,ixCmmin2,ixCmmin3,ixCmmax1,ixCmmax2,ixCmmax3
integer                            :: idim1,idim2,idir
double precision                   :: circ(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndim), dxidir
integer                            :: i1,i2,i3
!-----------------------------------------------------------------------------

associate(bfaces=>s%ws%w,x=>s%x%x,mygeo=>s%geo)

! Calculate contribution to FEM of each edge,
! that is, estimate value of line integral of
! electric field in the positive idir direction.

ixGsmin1=s%ws%ixGmin1;ixGsmin2=s%ws%ixGmin2;ixGsmin3=s%ws%ixGmin3
ixGsmax1=s%ws%ixGmax1;ixGsmax2=s%ws%ixGmax2;ixGsmax3=s%ws%ixGmax3;

! Loop over components of electric field

! idir: electric field component we need to calculate
! idim1: directions in which we already performed the reconstruction
! idim2: directions in which we perform the reconstruction

fE(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndir)=zero

do idir=7-2*ndim,ndir
  ! Indices
  ! idir: electric field component
  ! idim1: one surface
  ! idim2: the other surface
  ! cyclic permutation: idim1,idim2,idir=1,2,3
  ! Velocity components on the surface
  ! follow cyclic premutations:
  ! Sx(1),Sx(2)=y,z ; Sy(1),Sy(2)=z,x ; Sz(1),Sz(2)=x,y

  ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
  ixCmin1=ixOmin1-1+kr(idir,1);ixCmin2=ixOmin2-1+kr(idir,2)
  ixCmin3=ixOmin3-1+kr(idir,3);

  ! Set indices and directions
  idim1=mod(idir,ndir)+1
  idim2=mod(idir+1,ndir)+1

  jxCmin1=ixCmin1+kr(idim1,1);jxCmin2=ixCmin2+kr(idim1,2)
  jxCmin3=ixCmin3+kr(idim1,3);jxCmax1=ixCmax1+kr(idim1,1)
  jxCmax2=ixCmax2+kr(idim1,2);jxCmax3=ixCmax3+kr(idim1,3);
  ixCpmin1=ixCmin1+kr(idim2,1);ixCpmin2=ixCmin2+kr(idim2,2)
  ixCpmin3=ixCmin3+kr(idim2,3);ixCpmax1=ixCmax1+kr(idim2,1)
  ixCpmax2=ixCmax2+kr(idim2,2);ixCpmax3=ixCmax3+kr(idim2,3);

  ! Interpolate sqrt gamma to the edges

  sqrtg(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=zero

  select case(idim1)
 case(1)
    sqrtg(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=sqrtg&
       (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)+quarter*(&
      mygeo%mSurface1 %sqrtgamma(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)+&
      mygeo%mSurface1 %sqrtgamma(ixCpmin1:ixCpmax1,ixCpmin2:ixCpmax2,&
         ixCpmin3:ixCpmax3))

 case(2)
    sqrtg(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=sqrtg&
       (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)+quarter*(&
      mygeo%mSurface2 %sqrtgamma(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)+&
      mygeo%mSurface2 %sqrtgamma(ixCpmin1:ixCpmax1,ixCpmin2:ixCpmax2,&
         ixCpmin3:ixCpmax3))

 case(3)
    sqrtg(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=sqrtg&
       (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)+quarter*(&
      mygeo%mSurface3 %sqrtgamma(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)+&
      mygeo%mSurface3 %sqrtgamma(ixCpmin1:ixCpmax1,ixCpmin2:ixCpmax2,&
         ixCpmin3:ixCpmax3))

  end select
  
  select case(idim2)
 case(1)
    sqrtg(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=sqrtg&
       (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)+quarter*(&
      mygeo%mSurface1 %sqrtgamma(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)+&
      mygeo%mSurface1 %sqrtgamma(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
         jxCmin3:jxCmax3))

 case(2)
    sqrtg(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=sqrtg&
       (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)+quarter*(&
      mygeo%mSurface2 %sqrtgamma(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)+&
      mygeo%mSurface2 %sqrtgamma(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
         jxCmin3:jxCmax3))

 case(3)
    sqrtg(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=sqrtg&
       (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)+quarter*(&
      mygeo%mSurface3 %sqrtgamma(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)+&
      mygeo%mSurface3 %sqrtgamma(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
         jxCmin3:jxCmax3))

  end select

  ! Reconstruct transverse transport velocities
  call rec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,ixCmin2,&
     ixCmin3,ixCmax1,ixCmax2,ixCmax3,idim2,vbarC(ixImin1:ixImax1,&
     ixImin2:ixImax2,ixImin3:ixImax3,idim1,1),&
           vtilL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,2),&
              vtilR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,2))

  call rec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,ixCmin2,&
     ixCmin3,ixCmax1,ixCmax2,ixCmax3,idim1,vbarC(ixImin1:ixImax1,&
     ixImin2:ixImax2,ixImin3:ixImax3,idim2,2),&
           vtilL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1),&
              vtilR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1))

  ! Reconstruct magnetic fields
  ! Eventhough the arrays are larger, rec works with
  ! the limits ixG.
  call rec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,ixCmin2,&
     ixCmin3,ixCmax1,ixCmax2,ixCmax3,idim2,bfaces(ixImin1:ixImax1,&
     ixImin2:ixImax2,ixImin3:ixImax3,idim1),&
           btilL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim1),&
              btilR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim1))

  call rec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,ixCmin2,&
     ixCmin3,ixCmax1,ixCmax2,ixCmax3,idim1,bfaces(ixImin1:ixImax1,&
     ixImin2:ixImax2,ixImin3:ixImax3,idim2),&
           btilL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim2),&
              btilR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim2))

  ! Take the maximum characteristic

  cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)=max(cbarmin&
     (ixCpmin1:ixCpmax1,ixCpmin2:ixCpmax2,ixCpmin3:ixCpmax3,idim1),&
     cbarmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1))
  cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)=max(cbarmax&
     (ixCpmin1:ixCpmax1,ixCpmin2:ixCpmax2,ixCpmin3:ixCpmax3,idim1),&
     cbarmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1))

  cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2)=max(cbarmin&
     (jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,idim2),&
     cbarmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim2))
  cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2)=max(cbarmax&
     (jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,idim2),&
     cbarmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim2))
 
  ! Calculate elctric field
  if (idir .le. ndim) then
     dxidir = dxlevel(idir)
  else
     dxidir = 1.0d0
  end if

  fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)=qdt * dxidir * &
               sqrtg(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) * (&
               -(cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  1)*vtilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  1)*btilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  idim2) &
               + cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  1)*vtilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  1)*btilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  idim2) &
               -cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  1)*cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  1)*(btilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  idim2)-btilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3,idim2)))&
               /(cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  1) + cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  1)) &
               +(cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  2)*vtilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  2)*btilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  idim1) &
               + cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  2)*vtilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  2)*btilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  idim1) &
               -cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  2)*cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  2)*(btilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  idim1)-btilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3,idim1)))&
               /(cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  2) + cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  2)) &
               )

  if (typeaxial.ne.'slab') then


  ! Catch axis and origin, where the line integral is performed on a
  ! zero lenght element and should be zero.
  ! Remember that staggered quantities are assigned to the index
  ! of the cell at their left.
  where((abs(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
     z_)+half*dxlevel(z_)).lt.1.0d-9).or.&
        (abs(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           z_)+half*dxlevel(z_)-dpi).lt.1.0d-9))

    fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2)=zero
    fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,r_)=zero
  end where
  where(abs(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
     r_)+half*dxlevel(r_)).lt.1.0d-9)
    fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)=zero
  end where
  end if
!  

end do
!--------------------------------------------

circ(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)=zero

! Calculate circulation on each face

do idim1=1,ndim ! Coordinate perpendicular to face 
   ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
   ixCmin1=ixOmin1-kr(idim1,1);ixCmin2=ixOmin2-kr(idim1,2)
   ixCmin3=ixOmin3-kr(idim1,3);
   do idim2=1,ndim
      do idir=1,ndir ! Direction of line integral
         
        ! Assemble indices
        hxCmin1=ixCmin1-kr(idim2,1);hxCmin2=ixCmin2-kr(idim2,2)
        hxCmin3=ixCmin3-kr(idim2,3);hxCmax1=ixCmax1-kr(idim2,1)
        hxCmax2=ixCmax2-kr(idim2,2);hxCmax3=ixCmax3-kr(idim2,3);
        ! Add line integrals in direction idir
        circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
           =circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
                         +lvc(idim1,idim2,idir)&
                         *(fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                            idir)&
                         -fE(hxCmin1:hxCmax1,hxCmin2:hxCmax2,hxCmin3:hxCmax3,&
                            idir))
      end do
   end do
end do

! Divide by the area of the face to get dB/dt
do idim1=1,ndim
   ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
   ixCmin1=ixOmin1-kr(idim1,1);ixCmin2=ixOmin2-kr(idim1,2)
   ixCmin3=ixOmin3-kr(idim1,3);
   select case(idim1)
   case(1)
   where(mygeo%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
      .gt. 1.0d-9*mygeo%dvolume(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      ixCmin3:ixCmax3))
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
         =circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
                 /mygeo%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                    ixCmin3:ixCmax3)
   elsewhere
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)=zero
   end where
   
case(2)
   where(mygeo%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
      .gt. 1.0d-9*mygeo%dvolume(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      ixCmin3:ixCmax3))
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
         =circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
                 /mygeo%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                    ixCmin3:ixCmax3)
   elsewhere
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)=zero
   end where
   
case(3)
   where(mygeo%surfaceC3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
      .gt. 1.0d-9*mygeo%dvolume(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      ixCmin3:ixCmax3))
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
         =circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
                 /mygeo%surfaceC3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                    ixCmin3:ixCmax3)
   elsewhere
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)=zero
   end where
   
    end select

   ! Time update
   ! minus!
   bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
      =bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
      idim1)-circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)
end do

end associate
end subroutine updatefacesuct2
!=============================================================================
subroutine updatefacesuct2av(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qdt,vbarC,cbarmin,cbarmax,&
   fC,fE,s)

use mod_amrvacdef

integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(in)       :: qdt
type(state)                        :: s
double precision, intent(in)       :: vbarC(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndir,2)
double precision, intent(in)       :: cbarmin(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,ndim)
double precision, intent(in)       :: cbarmax(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,ndim)
double precision, intent(in)       :: fC(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nwflux,1:ndim)
double precision, intent(inout)    :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndir)
! --- local ---
double precision                   :: vtilL(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,2)
double precision                   :: vtilR(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,2)
double precision                   :: btilL(s%ws%ixGmin1:s%ws%ixGmax1,&
   s%ws%ixGmin2:s%ws%ixGmax2,s%ws%ixGmin3:s%ws%ixGmax3,ndim)
double precision                   :: btilR(s%ws%ixGmin1:s%ws%ixGmax1,&
   s%ws%ixGmin2:s%ws%ixGmax2,s%ws%ixGmin3:s%ws%ixGmax3,ndim)
double precision                   :: sqrtg(s%ws%ixGmin1:s%ws%ixGmax1,&
   s%ws%ixGmin2:s%ws%ixGmax2,s%ws%ixGmin3:s%ws%ixGmax3)
double precision                   :: cp(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,2)
double precision                   :: cm(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,2)
double precision                   :: fEc(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)
double precision                   :: fEvar(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)
double precision                   :: fEmin(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)
double precision                   :: fEmax(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)
integer                            :: ixGsmin1,ixGsmin2,ixGsmin3,ixGsmax1,&
   ixGsmax2,ixGsmax3
integer                            :: hxCmin1,hxCmin2,hxCmin3,hxCmax1,hxCmax2,&
   hxCmax3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,ixCpmin1,ixCpmin2,&
   ixCpmin3,ixCpmax1,ixCpmax2,ixCpmax3,jxCmin1,jxCmin2,jxCmin3,jxCmax1,&
   jxCmax2,jxCmax3,jxCpmin1,jxCpmin2,jxCpmin3,jxCpmax1,jxCpmax2,jxCpmax3,&
   ixCmmin1,ixCmmin2,ixCmmin3,ixCmmax1,ixCmmax2,ixCmmax3
integer                            :: idim1,idim2,idir,iwdim1,iwdim2
double precision                   :: circ(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndim), dxidir
integer                            :: i1,i2,i3
!-----------------------------------------------------------------------------

associate(bfaces=>s%ws%w,w=>s%w%w,x=>s%x%x,mygeo=>s%geo)
! Calculate contribution to FEM of each edge,
! that is, estimate value of line integral of
! electric field in the positive idir direction.

ixGsmin1=s%ws%ixGmin1;ixGsmin2=s%ws%ixGmin2;ixGsmin3=s%ws%ixGmin3
ixGsmax1=s%ws%ixGmax1;ixGsmax2=s%ws%ixGmax2;ixGsmax3=s%ws%ixGmax3;

! Loop over components of electric field

! idir: electric field component we need to calculate
! idim1: directions in which we already performed the reconstruction
! idim2: directions in which we perform the reconstruction

fE(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndir)=zero

do idir=7-2*ndim,ndir
  ! Indices
  ! idir: electric field component
  ! idim1: one surface
  ! idim2: the other surface
  ! cyclic permutation: idim1,idim2,idir=1,2,3
  ! Velocity components on the surface
  ! follow cyclic premutations:
  ! Sx(1),Sx(2)=y,z ; Sy(1),Sy(2)=z,x ; Sz(1),Sz(2)=x,y

  ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
  ixCmin1=ixOmin1-1+kr(idir,1);ixCmin2=ixOmin2-1+kr(idir,2)
  ixCmin3=ixOmin3-1+kr(idir,3);

  ! Set indices and directions
  idim1=mod(idir,ndir)+1
  idim2=mod(idir+1,ndir)+1

  jxCmin1=ixCmin1+kr(idim1,1);jxCmin2=ixCmin2+kr(idim1,2)
  jxCmin3=ixCmin3+kr(idim1,3);jxCmax1=ixCmax1+kr(idim1,1)
  jxCmax2=ixCmax2+kr(idim1,2);jxCmax3=ixCmax3+kr(idim1,3);
  ixCpmin1=ixCmin1+kr(idim2,1);ixCpmin2=ixCmin2+kr(idim2,2)
  ixCpmin3=ixCmin3+kr(idim2,3);ixCpmax1=ixCmax1+kr(idim2,1)
  ixCpmax2=ixCmax2+kr(idim2,2);ixCpmax3=ixCmax3+kr(idim2,3);
  jxCpmin1=jxCmin1+kr(idim2,1);jxCpmin2=jxCmin2+kr(idim2,2)
  jxCpmin3=jxCmin3+kr(idim2,3);jxCpmax1=jxCmax1+kr(idim2,1)
  jxCpmax2=jxCmax2+kr(idim2,2);jxCpmax3=jxCmax3+kr(idim2,3);

  ! Interpolate sqrt gamma to the edges

  sqrtg(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=zero

  select case(idim1)
 case(1)
    sqrtg(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=sqrtg&
       (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)+quarter*(&
      mygeo%mSurface1 %sqrtgamma(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)+&
      mygeo%mSurface1 %sqrtgamma(ixCpmin1:ixCpmax1,ixCpmin2:ixCpmax2,&
         ixCpmin3:ixCpmax3))

 case(2)
    sqrtg(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=sqrtg&
       (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)+quarter*(&
      mygeo%mSurface2 %sqrtgamma(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)+&
      mygeo%mSurface2 %sqrtgamma(ixCpmin1:ixCpmax1,ixCpmin2:ixCpmax2,&
         ixCpmin3:ixCpmax3))

 case(3)
    sqrtg(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=sqrtg&
       (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)+quarter*(&
      mygeo%mSurface3 %sqrtgamma(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)+&
      mygeo%mSurface3 %sqrtgamma(ixCpmin1:ixCpmax1,ixCpmin2:ixCpmax2,&
         ixCpmin3:ixCpmax3))

  end select
  
  select case(idim2)
 case(1)
    sqrtg(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=sqrtg&
       (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)+quarter*(&
      mygeo%mSurface1 %sqrtgamma(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)+&
      mygeo%mSurface1 %sqrtgamma(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
         jxCmin3:jxCmax3))

 case(2)
    sqrtg(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=sqrtg&
       (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)+quarter*(&
      mygeo%mSurface2 %sqrtgamma(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)+&
      mygeo%mSurface2 %sqrtgamma(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
         jxCmin3:jxCmax3))

 case(3)
    sqrtg(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=sqrtg&
       (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)+quarter*(&
      mygeo%mSurface3 %sqrtgamma(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)+&
      mygeo%mSurface3 %sqrtgamma(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
         jxCmin3:jxCmax3))

  end select

  ! Reconstruct transverse transport velocities
  call rec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,ixCmin2,&
     ixCmin3,ixCmax1,ixCmax2,ixCmax3,idim2,vbarC(ixImin1:ixImax1,&
     ixImin2:ixImax2,ixImin3:ixImax3,idim1,1),&
           vtilL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,2),&
              vtilR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,2))

  call rec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,ixCmin2,&
     ixCmin3,ixCmax1,ixCmax2,ixCmax3,idim1,vbarC(ixImin1:ixImax1,&
     ixImin2:ixImax2,ixImin3:ixImax3,idim2,2),&
           vtilL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1),&
              vtilR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1))

  ! Reconstruct magnetic fields
  ! Eventhough the arrays are larger, rec works with
  ! the limits ixG.
  call rec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,ixCmin2,&
     ixCmin3,ixCmax1,ixCmax2,ixCmax3,idim2,bfaces(ixImin1:ixImax1,&
     ixImin2:ixImax2,ixImin3:ixImax3,idim1),&
           btilL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim1),&
              btilR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim1))

  call rec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,ixCmin2,&
     ixCmin3,ixCmax1,ixCmax2,ixCmax3,idim1,bfaces(ixImin1:ixImax1,&
     ixImin2:ixImax2,ixImin3:ixImax3,idim2),&
           btilL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim2),&
              btilR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim2))

  ! Take the maximum characteristic

  cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)=max(cbarmin&
     (ixCpmin1:ixCpmax1,ixCpmin2:ixCpmax2,ixCpmin3:ixCpmax3,idim1),&
     cbarmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1))
  cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)=max(cbarmax&
     (ixCpmin1:ixCpmax1,ixCpmin2:ixCpmax2,ixCpmin3:ixCpmax3,idim1),&
     cbarmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1))

  cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2)=max(cbarmin&
     (jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,idim2),&
     cbarmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim2))
  cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2)=max(cbarmax&
     (jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,idim2),&
     cbarmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim2))
 
  ! Calculate elctric field
  if (idir .le. ndim) then
     dxidir = dxlevel(idir)
  else
     dxidir = 1.0d0
  end if

  fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)=qdt * dxidir * &
               sqrtg(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) * (&
               -(cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  1)*vtilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  1)*btilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  idim2) &
               + cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  1)*vtilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  1)*btilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  idim2) &
               -cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  1)*cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  1)*(btilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  idim2)-btilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3,idim2)))&
               /(cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  1) + cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  1)) &
               +(cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  2)*vtilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  2)*btilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  idim1) &
               + cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  2)*vtilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  2)*btilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  idim1) &
               -cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  2)*cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  2)*(btilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  idim1)-btilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3,idim1)))&
               /(cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  2) + cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  2)) &
               )

  
  if (typeaxial.ne.'slab') then


  ! Catch axis and origin, where the line integral is performed on a
  ! zero lenght element and should be zero.
  ! Remember that staggered quantities are assigned to the index
  ! of the cell at their left.
  where((abs(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
     z_)+half*dxlevel(z_)).lt.1.0d-9).or.&
        (abs(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           z_)+half*dxlevel(z_)-dpi).lt.1.0d-9))

    fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2)=zero
    fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,r_)=zero
  end where
  where(abs(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
     r_)+half*dxlevel(r_)).lt.1.0d-9)
    fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)=zero
  end where
  end if

  ! Allowed range for the electric field in that direction.
  ! For simplicity, use twice cell-centred values, to allow
  ! some overshooting due to reconstruction
  ! Recall: 1,2,3 --> idir, idim1, idim2

  fEc(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) = myM%beta(idim2)%elem&
     (ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)*w(ixImin1:ixImax1,&
     ixImin2:ixImax2,ixImin3:ixImax3,b0_+idim1) &
             - myM%beta(idim1)%elem(ixImin1:ixImax1,ixImin2:ixImax2,&
                ixImin3:ixImax3)*w(ixImin1:ixImax1,ixImin2:ixImax2,&
                ixImin3:ixImax3,b0_+idim2)

  fEvar(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) = &
     2.0*max(abs(w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
     b0_+idim2) - w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
     b0_+idim1)),&
                     abs(w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
                        b0_+idim1)),&
                     abs(w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
                        b0_+idim2)),1e-6)

  fEmax(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) = &
     fEc(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) + &
     myM%alpha(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) &
             *fEvar(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

  fEmin(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) = &
     fEc(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) - &
     myM%alpha(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) &
             *fEvar(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

 !if (any(fEmax(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3).ne.fEmax(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3))) print *,'fEmax is NaN','it=',it
 !if (any(abs(fEmax(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)).gt.bigdouble)) print *,'fEmax is inf','it=',it
!
 !if (any(fEmin(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3).ne.fEmin(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3))) print *,'fEmin is NaN','it=',it
 !if (any(abs(fEmin(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)).gt.bigdouble)) print *,'fEmin is inf','it=',it



  ! For the edge-allocated electric fields, take the maximum and
  ! minimum of the surrounding cells



  fEmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = qdt * dxidir * &
     sqrtg(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) * &
     max(fEmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3),&
     fEmax(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3),&
     fEmax(ixCpmin1:ixCpmax1,ixCpmin2:ixCpmax2,ixCpmin3:ixCpmax3),&
     fEmax(jxCpmin1:jxCpmax1,jxCpmin2:jxCpmax2,jxCpmin3:jxCpmax3)) 

  fEmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = qdt * dxidir * &
     sqrtg(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) * &
     min(fEmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3),&
     fEmin(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3),&
     fEmin(ixCpmin1:ixCpmax1,ixCpmin2:ixCpmax2,ixCpmin3:ixCpmax3),&
     fEmin(jxCpmin1:jxCpmax1,jxCpmin2:jxCpmax2,jxCpmin3:jxCpmax3)) 


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Fall back to Balsara-Spicer averaging
  ! in case of problems
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (any(.not.(fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
     idir).lt.fEmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
     .and.fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
     idir).gt.fEmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)))) then
!     write(*,*) mype, 'WARNING: Had to fall back to Balsara-Spicer','it=',it
     iwdim1 = b0_+idim1; iwdim2 = b0_+idim2
     where (.not.(fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
        idir).lt.fEmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
        .and.fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
        idir).gt.fEmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)))
        fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)&
           =qdt*quarter*((fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           iwdim1,idim2)&
             +fC(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,iwdim1,&
                idim2))/dxlevel(idim1)&
             -(fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iwdim2,&
                idim1)+fC(ixCpmin1:ixCpmax1,ixCpmin2:ixCpmax2,&
                ixCpmin3:ixCpmax3,iwdim2,idim1))&
             /dxlevel(idim2))
     end where
  end if

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Before clipping the electric field,
  ! forcefully remove NaNs that could be there
  ! Brutalistic (probably never helps
  ! as NaNs will also be in other variables then)
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (any(fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
     idir) .ne. fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
     idir))) then
!     write(*,*) mype, 'WARNING: Had to forecefully set E=0','it=',it
     where (fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
        idir) .ne. fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir))
        fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir) = zero
     end where
  end if
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! If even Balsara-Spicer yields electric
  ! field out of bounds, reset to a fraction
  ! of the allowed field
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (any(.not.(fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
     idir).lt.fEmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
     .and.fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
     idir).gt.fEmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)))) then
!     write(*,*) mype, 'WARNING: Clipping the electric field','it=',it
     where (fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
        idir).gt.fEmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3))
        fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)&
           =fEmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) + 0.8 * &
           (fEmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) - &
           fEmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3))
     end where

     where (fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
        idir).lt.fEmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3))
        fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)&
           =fEmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) + 0.2 * &
           (fEmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) - &
           fEmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3))
     end where

  end if

 !if (any(fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir).ne.fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir))) print *,'fE is NaN','it=',it
 !if (any(abs(fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)).gt.bigdouble)) print *,'fE is inf','it=',it




end do
!--------------------------------------------

circ(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)=zero

! Calculate circulation on each face

do idim1=1,ndim ! Coordinate perpendicular to face 
   ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
   ixCmin1=ixOmin1-kr(idim1,1);ixCmin2=ixOmin2-kr(idim1,2)
   ixCmin3=ixOmin3-kr(idim1,3);
   do idim2=1,ndim
      do idir=1,ndir ! Direction of line integral
         
        ! Assemble indices
        hxCmin1=ixCmin1-kr(idim2,1);hxCmin2=ixCmin2-kr(idim2,2)
        hxCmin3=ixCmin3-kr(idim2,3);hxCmax1=ixCmax1-kr(idim2,1)
        hxCmax2=ixCmax2-kr(idim2,2);hxCmax3=ixCmax3-kr(idim2,3);
        ! Add line integrals in direction idir
        circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
           =circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
                         +lvc(idim1,idim2,idir)&
                         *(fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                            idir)&
                         -fE(hxCmin1:hxCmax1,hxCmin2:hxCmax2,hxCmin3:hxCmax3,&
                            idir))
      end do
   end do
end do

! Divide by the area of the face to get dB/dt
do idim1=1,ndim
   ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
   ixCmin1=ixOmin1-kr(idim1,1);ixCmin2=ixOmin2-kr(idim1,2)
   ixCmin3=ixOmin3-kr(idim1,3);
   select case(idim1)
   case(1)
   where(mygeo%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
      .gt. 1.0d-9*mygeo%dvolume(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      ixCmin3:ixCmax3))
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
         =circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
                 /mygeo%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                    ixCmin3:ixCmax3)
   elsewhere
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)=zero
   end where
   
case(2)
   where(mygeo%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
      .gt. 1.0d-9*mygeo%dvolume(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      ixCmin3:ixCmax3))
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
         =circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
                 /mygeo%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                    ixCmin3:ixCmax3)
   elsewhere
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)=zero
   end where
   
case(3)
   where(mygeo%surfaceC3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
      .gt. 1.0d-9*mygeo%dvolume(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      ixCmin3:ixCmax3))
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
         =circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
                 /mygeo%surfaceC3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                    ixCmin3:ixCmax3)
   elsewhere
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)=zero
   end where
   
    end select

   ! Time update
   ! minus!
    bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
       =bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
       idim1)-circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)

 !if (any(bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1) .ne. bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1))) &
!       write(*,*) mype, 'WARNING: Magnetic field is NaN in uct2av'

    
end do

end associate
end subroutine updatefacesuct2av
!=============================================================================
subroutine updatefacesuct1(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qdt,vbarRC,vbarLC,cbarmin,&
   cbarmax,fE,s)

use mod_amrvacdef

integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(in)       :: qdt
type(state)                        :: s
double precision, intent(in)       :: vbarRC(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndir,2)
double precision, intent(in)       :: vbarLC(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndir,2)
double precision, intent(in)       :: cbarmin(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,ndim)
double precision, intent(in)       :: cbarmax(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,ndim)
double precision, intent(inout)    :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndir)
! --- local ---
double precision                   :: vtilRL(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,2)
double precision                   :: vtilLL(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,2)
double precision                   :: vtilRR(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,2)
double precision                   :: vtilLR(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,2)
double precision                   :: btilL(s%ws%ixGmin1:s%ws%ixGmax1,&
   s%ws%ixGmin2:s%ws%ixGmax2,s%ws%ixGmin3:s%ws%ixGmax3,ndim)
double precision                   :: btilR(s%ws%ixGmin1:s%ws%ixGmax1,&
   s%ws%ixGmin2:s%ws%ixGmax2,s%ws%ixGmin3:s%ws%ixGmax3,ndim)
double precision                   :: sqrtg(s%ws%ixGmin1:s%ws%ixGmax1,&
   s%ws%ixGmin2:s%ws%ixGmax2,s%ws%ixGmin3:s%ws%ixGmax3)
double precision                   :: cp(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,2)
double precision                   :: cm(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,2)
double precision                   :: ELL(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)
double precision                   :: ELR(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)
double precision                   :: ERL(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)
double precision                   :: ERR(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)
integer                            :: ixGsmin1,ixGsmin2,ixGsmin3,ixGsmax1,&
   ixGsmax2,ixGsmax3
integer                            :: hxCmin1,hxCmin2,hxCmin3,hxCmax1,hxCmax2,&
   hxCmax3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,ixCpmin1,ixCpmin2,&
   ixCpmin3,ixCpmax1,ixCpmax2,ixCpmax3,jxCmin1,jxCmin2,jxCmin3,jxCmax1,&
   jxCmax2,jxCmax3,ixCmmin1,ixCmmin2,ixCmmin3,ixCmmax1,ixCmmax2,ixCmmax3
integer                            :: idim1,idim2,idir
double precision                   :: circ(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndim), dxidir

!-----------------------------------------------------------------------------

associate(bfaces=>s%ws%w,x=>s%x%x,mygeo=>s%geo)

! Calculate contribution to FEM of each edge,
! that is, estimate value of line integral of
! electric field in the positive idir direction.

ixGsmin1=s%ws%ixGmin1;ixGsmin2=s%ws%ixGmin2;ixGsmin3=s%ws%ixGmin3
ixGsmax1=s%ws%ixGmax1;ixGsmax2=s%ws%ixGmax2;ixGsmax3=s%ws%ixGmax3;

! Loop over components of electric field

! idir: electric field component we need to calculate
! idim1: directions in which we already performed the reconstruction
! idim2: directions in which we perform the reconstruction

fE(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndir)=zero

do idir=7-2*ndim,ndir
  ! Indices
  ! idir: electric field component
  ! idim1: where the first reconstruction was done
  ! idim2: where we will do the second reconstruction
  ! cyclic permutation: idim1,idim2,idir=1,2,3

  ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
  ixCmin1=ixOmin1-1+kr(idir,1);ixCmin2=ixOmin2-1+kr(idir,2)
  ixCmin3=ixOmin3-1+kr(idir,3);

 ! Set indices and directions
  idim1=mod(idir,ndir)+1
  idim2=mod(idir+1,ndir)+1

  jxCmin1=ixCmin1+kr(idim1,1);jxCmin2=ixCmin2+kr(idim1,2)
  jxCmin3=ixCmin3+kr(idim1,3);jxCmax1=ixCmax1+kr(idim1,1)
  jxCmax2=ixCmax2+kr(idim1,2);jxCmax3=ixCmax3+kr(idim1,3);
  ixCpmin1=ixCmin1+kr(idim2,1);ixCpmin2=ixCmin2+kr(idim2,2)
  ixCpmin3=ixCmin3+kr(idim2,3);ixCpmax1=ixCmax1+kr(idim2,1)
  ixCpmax2=ixCmax2+kr(idim2,2);ixCpmax3=ixCmax3+kr(idim2,3);

  ! Interpolate sqrt gamma to the edges

  sqrtg(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=zero

  select case(idim1)
 case(1)
    sqrtg(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=sqrtg&
       (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)+quarter*(&
      mygeo%mSurface1 %sqrtgamma(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)+&
      mygeo%mSurface1 %sqrtgamma(ixCpmin1:ixCpmax1,ixCpmin2:ixCpmax2,&
         ixCpmin3:ixCpmax3))

 case(2)
    sqrtg(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=sqrtg&
       (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)+quarter*(&
      mygeo%mSurface2 %sqrtgamma(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)+&
      mygeo%mSurface2 %sqrtgamma(ixCpmin1:ixCpmax1,ixCpmin2:ixCpmax2,&
         ixCpmin3:ixCpmax3))

 case(3)
    sqrtg(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=sqrtg&
       (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)+quarter*(&
      mygeo%mSurface3 %sqrtgamma(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)+&
      mygeo%mSurface3 %sqrtgamma(ixCpmin1:ixCpmax1,ixCpmin2:ixCpmax2,&
         ixCpmin3:ixCpmax3))

  end select
  
  select case(idim2)
 case(1)
    sqrtg(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=sqrtg&
       (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)+quarter*(&
      mygeo%mSurface1 %sqrtgamma(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)+&
      mygeo%mSurface1 %sqrtgamma(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
         jxCmin3:jxCmax3))

 case(2)
    sqrtg(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=sqrtg&
       (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)+quarter*(&
      mygeo%mSurface2 %sqrtgamma(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)+&
      mygeo%mSurface2 %sqrtgamma(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
         jxCmin3:jxCmax3))

 case(3)
    sqrtg(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=sqrtg&
       (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)+quarter*(&
      mygeo%mSurface3 %sqrtgamma(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)+&
      mygeo%mSurface3 %sqrtgamma(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
         jxCmin3:jxCmax3))

  end select

  ! Reconstruct velocities
  call rec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,ixCmin2,&
     ixCmin3,ixCmax1,ixCmax2,ixCmax3,idim2,vbarLC(ixImin1:ixImax1,&
     ixImin2:ixImax2,ixImin3:ixImax3,idir,1),&
           vtilLL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1),&
              vtilLR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1)) 
  call rec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,ixCmin2,&
     ixCmin3,ixCmax1,ixCmax2,ixCmax3,idim2,vbarRC(ixImin1:ixImax1,&
     ixImin2:ixImax2,ixImin3:ixImax3,idir,1),&
           vtilRL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1),&
              vtilRR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1))

  call rec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,ixCmin2,&
     ixCmin3,ixCmax1,ixCmax2,ixCmax3,idim2,vbarLC(ixImin1:ixImax1,&
     ixImin2:ixImax2,ixImin3:ixImax3,idir,2),&
           vtilLL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,2),&
              vtilLR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,2)) 
  call rec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,ixCmin2,&
     ixCmin3,ixCmax1,ixCmax2,ixCmax3,idim2,vbarRC(ixImin1:ixImax1,&
     ixImin2:ixImax2,ixImin3:ixImax3,idir,2),&
           vtilRL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,2),&
              vtilRR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,2))

  ! Reconstruct magnetic fields
  ! Eventhough the arrays are larger, rec works with
  ! the limits ixG.
  call rec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,ixCmin2,&
     ixCmin3,ixCmax1,ixCmax2,ixCmax3,idim2,bfaces(ixImin1:ixImax1,&
     ixImin2:ixImax2,ixImin3:ixImax3,idim1),&
           btilL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim1),&
              btilR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim1))
  call rec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,ixCmin2,&
     ixCmin3,ixCmax1,ixCmax2,ixCmax3,idim1,bfaces(ixImin1:ixImax1,&
     ixImin2:ixImax2,ixImin3:ixImax3,idim2),&
           btilL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim2),&
              btilR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim2))

  ! Take the maximum characteristic
  cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)=max(cbarmin&
     (ixCpmin1:ixCpmax1,ixCpmin2:ixCpmax2,ixCpmin3:ixCpmax3,idim1),&
     cbarmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)) 
  cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)=max(cbarmax&
     (ixCpmin1:ixCpmax1,ixCpmin2:ixCpmax2,ixCpmin3:ixCpmax3,idim1),&
     cbarmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1))

  cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2)=max(cbarmin&
     (jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,idim2),&
     cbarmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim2)) 
  cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2)=max(cbarmax&
     (jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,idim2),&
     cbarmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim2))

  ! Calculate component idir of vxB partial
  ELL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)= vtilLL&
     (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2)*btilL&
     (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
       -vtilLL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
          1)*btilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim2)

  ELR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)= vtilLR&
     (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2)*btilR&
     (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
       -vtilLR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
          1)*btilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim2)

  ERL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)= vtilRL&
     (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2)*btilL&
     (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
       -vtilRL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
          1)*btilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim2)

  ERR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)= vtilRR&
     (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2)*btilR&
     (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
       -vtilRR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
          1)*btilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim2)

  ! For 3D, use interval
  if (idir .le. ndim) then
     dxidir = dxlevel(idir)
  else
     dxidir = 1.0d0
  end if
  
  ! Calculate elctric field
  fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)=qdt * dxidir * &
   sqrtg(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) * (&
   (cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)*cp(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3,2)*ELL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      ixCmin3:ixCmax3)&
   +cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)*cm(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3,2)*ELR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      ixCmin3:ixCmax3)&
   +cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)*cp(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3,2)*ERL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      ixCmin3:ixCmax3)&
   +cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)*cm(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3,2)*ERR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      ixCmin3:ixCmax3)&
   )/&
    ((cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)+cm&
       (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1))*&
     (cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2)+cm&
        (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2)))&
   +(cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)*cm(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)*(btilR(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim2)-btilL(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim2)))/&
     (cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1) + &
        cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1))&
   -(cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2)*cm(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3,2)*(btilR(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)-btilL(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)))/&
     (cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2) + &
        cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2)))

  if (typeaxial.ne.'slab') then

  ! Catch axis and origin, where the line integral is performed on a
  ! zero lenght element and should be zero.
  ! Remember that staggered quantities are assigned to the index
  ! of the cell at their left.
  where((abs(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
     z_)+half*dxlevel(z_)).lt.1.0d-9).or.&
        (abs(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           z_)+half*dxlevel(z_)-dpi).lt.1.0d-9))
    fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2)=zero
    fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,r_)=zero
  end where
  where(abs(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
     r_)+half*dxlevel(r_)).lt.1.0d-9)
    fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)=zero
  end where
  end if

end do


!--------------------------------------------

circ(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)=zero

! Calculate circulation on each face

do idim1=1,ndim ! Coordinate perpendicular to face 
   ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
   ixCmin1=ixOmin1-kr(idim1,1);ixCmin2=ixOmin2-kr(idim1,2)
   ixCmin3=ixOmin3-kr(idim1,3);
   do idim2=1,ndim
      do idir=1,ndir ! Direction of line integral
                  
        ! Assemble indices
        hxCmin1=ixCmin1-kr(idim2,1);hxCmin2=ixCmin2-kr(idim2,2)
        hxCmin3=ixCmin3-kr(idim2,3);hxCmax1=ixCmax1-kr(idim2,1)
        hxCmax2=ixCmax2-kr(idim2,2);hxCmax3=ixCmax3-kr(idim2,3);
        ! Add line integrals in direction idir
        circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
           =circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
                         +lvc(idim1,idim2,idir)&
                         *(fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                            idir)&
                         -fE(hxCmin1:hxCmax1,hxCmin2:hxCmax2,hxCmin3:hxCmax3,&
                            idir))
      end do
   end do
end do

! Divide by the area of the face to get dB/dt
do idim1=1,ndim
   ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
   ixCmin1=ixOmin1-kr(idim1,1);ixCmin2=ixOmin2-kr(idim1,2)
   ixCmin3=ixOmin3-kr(idim1,3);
   select case(idim1)
   case(1)
   where(mygeo%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
      .gt. 1.0d-9*mygeo%dvolume(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      ixCmin3:ixCmax3))
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
         =circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
                 /mygeo%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                    ixCmin3:ixCmax3)
   elsewhere
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)=zero
   end where
   
case(2)
   where(mygeo%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
      .gt. 1.0d-9*mygeo%dvolume(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      ixCmin3:ixCmax3))
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
         =circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
                 /mygeo%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                    ixCmin3:ixCmax3)
   elsewhere
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)=zero
   end where
   
case(3)
   where(mygeo%surfaceC3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
      .gt. 1.0d-9*mygeo%dvolume(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      ixCmin3:ixCmax3))
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
         =circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
                 /mygeo%surfaceC3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                    ixCmin3:ixCmax3)
   elsewhere
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)=zero
   end where
   
    end select

   ! Time update
   ! minus!
   bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
      =bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
      idim1)-circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)
end do

end associate
end subroutine updatefacesuct1
!=============================================================================
subroutine storeedge(igrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,fE,&
   idimmin,idimmax)
!use mod_fix_conserve
use mod_amrvacdef

integer, intent(in)          :: igrid, ixImin1,ixImin2,ixImin3,ixImax1,&
   ixImax2,ixImax3, idimmin,idimmax
double precision, intent(in) :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:3)

integer :: idims, idir, iside, i1,i2,i3
integer :: pi1,pi2,pi3, mi1,mi2,mi3, ph1,ph2,ph3, mh1,mh2,mh3 !To detect corners
integer :: ixMcmin1,ixMcmin2,ixMcmin3,ixMcmax1,ixMcmax2,ixMcmax3
!-----------------------------------------------------------------------------

do idims = idimmin,idimmax  !loop over face directions
  !! Loop over block faces
  do iside=1,2 
    i1=kr(1,idims)*(2*iside-3);i2=kr(2,idims)*(2*iside-3)
    i3=kr(3,idims)*(2*iside-3);
     if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle
    select case (neighbor_type(i1,i2,i3,igrid))
    case (4)
      ! The neighbour is finer
      ! Face direction, side (left or right), restrict required?, fE
      call fluxtoedge(igrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         idims,iside,.false.,fE)
    case(2)
      ! The neighbour is coarser
      call fluxtoedge(igrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         idims,iside,.true.,fE)
    case(3)
      ! If the neighbour is at the same level,
      ! check if there are corners
      ! If there is any corner, store the fluxes from that side
      do idir=idims+1,ndim
        pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
        mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
        ph1=pi1-kr(idims,1)*(2*iside-3);ph2=pi2-kr(idims,2)*(2*iside-3)
        ph3=pi3-kr(idims,3)*(2*iside-3);
        mh1=mi1-kr(idims,1)*(2*iside-3);mh2=mi2-kr(idims,2)*(2*iside-3)
        mh3=mi3-kr(idims,3)*(2*iside-3);
        if (neighbor_type(pi1,pi2,pi3,igrid).eq.4&
          .and.neighbor_type(ph1,ph2,ph3,igrid).eq.3) then
          call fluxtoedge(igrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
             ixImax3,idims,iside,.false.,fE)
        end if
        if (neighbor_type(mi1,mi2,mi3,igrid).eq.4&
          .and.neighbor_type(mh1,mh2,mh3,igrid).eq.3) then
          call fluxtoedge(igrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
             ixImax3,idims,iside,.false.,fE)
        end if
      end do
    end select
  end do
end do

end subroutine storeedge
!=============================================================================
subroutine fluxtoedge(igrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   idims,iside,restrict,fE)
use mod_fix_conserve
use mod_amrvacdef

integer                      :: igrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3,idims,iside
logical                      :: restrict
double precision, intent(in) :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:3)
! ... local ...
integer                      :: idir1,idir2
integer                      :: ixEmin1,ixEmin2,ixEmin3,ixEmax1,ixEmax2,&
   ixEmax3,ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3, jxFmin1,jxFmin2,&
   jxFmin3,jxFmax1,jxFmax2,jxFmax3, nx1,nx2,nx3,nxCo1,nxCo2,nxCo3
!-----------------------------------------------------------------------------

nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;nx3=ixMhi3-ixMlo3+1;
nxCo1=nx1/2;nxCo2=nx2/2;nxCo3=nx3/2;

! ixE are the indices on the 'edge' array.
! ixF are the indices on the 'fE' array
! jxF are indices advanced to perform the flux restriction (sum) in 3D
! A line integral of the electric field on the coarse side
! lies over two edges on the fine side. So, in 3D we restrict by summing
! over two cells on the fine side.

do idir1=1,ndim-1
  ! 3D: rotate indices among 1 and 2 to save space 
  idir2=mod(idir1+idims-1,3)+1
 

  if (restrict) then
    ! Set up indices for restriction
    ixFmin1=ixMlo1-1+kr(1,idir2);ixFmin2=ixMlo2-1+kr(2,idir2)
    ixFmin3=ixMlo3-1+kr(3,idir2);
    ixFmax1=ixMhi1-kr(1,idir2);ixFmax2=ixMhi2-kr(2,idir2)
    ixFmax3=ixMhi3-kr(3,idir2);
    
    jxFmin1=ixFmin1+kr(1,idir2);jxFmin2=ixFmin2+kr(2,idir2)
    jxFmin3=ixFmin3+kr(3,idir2);jxFmax1=ixFmax1+kr(1,idir2)
    jxFmax2=ixFmax2+kr(2,idir2);jxFmax3=ixFmax3+kr(3,idir2);

    ixEmin1=0+kr(1,idir2);ixEmin2=0+kr(2,idir2);ixEmin3=0+kr(3,idir2);
    ixEmax1=nxCo1;ixEmax2=nxCo2;ixEmax3=nxCo3;
    select case(idims)
   case(1)
      ixEmin1=1;ixEmax1=1;
      select case(iside)
      case(1)
        ixFmax1=ixFmin1
        jxFmax1=ixFmin1
      case(2)
        ixFmin1=ixFmax1
        jxFmin1=ixFmax1
      end select
   
case(2)
      ixEmin2=1;ixEmax2=1;
      select case(iside)
      case(1)
        ixFmax2=ixFmin2
        jxFmax2=ixFmin2
      case(2)
        ixFmin2=ixFmax2
        jxFmin2=ixFmax2
      end select
   
case(3)
      ixEmin3=1;ixEmax3=1;
      select case(iside)
      case(1)
        ixFmax3=ixFmin3
        jxFmax3=ixFmin3
      case(2)
        ixFmin3=ixFmax3
        jxFmin3=ixFmax3
      end select
   
    end select

  pflux(iside,idims,igrid)%edge(ixEmin1:ixEmax1,ixEmin2:ixEmax2,&
     ixEmin3:ixEmax3,idir1)=&
    fE(ixFmin1:ixFmax1:2,ixFmin2:ixFmax2:2,ixFmin3:ixFmax3:2,idir2) +&
    fE(jxFmin1:jxFmax1:2,jxFmin2:jxFmax2:2,jxFmin3:jxFmax3:2,idir2);

  else
    ! Set up indices for copying 
    ixFmin1=ixMlo1-1+kr(1,idir2);ixFmin2=ixMlo2-1+kr(2,idir2)
    ixFmin3=ixMlo3-1+kr(3,idir2);
    ixFmax1=ixMhi1;ixFmax2=ixMhi2;ixFmax3=ixMhi3;

    ixEmin1=0+kr(1,idir2);ixEmin2=0+kr(2,idir2);ixEmin3=0+kr(3,idir2);
    ixEmax1=nx1;ixEmax2=nx2;ixEmax3=nx3;

    select case(idims)
   case(1)
      ixEmin1=1;ixEmax1=1;
      select case(iside)
      case(1)
        ixFmax1=ixFmin1
      case(2)
        ixFmin1=ixFmax1
      end select
   
case(2)
      ixEmin2=1;ixEmax2=1;
      select case(iside)
      case(1)
        ixFmax2=ixFmin2
      case(2)
        ixFmin2=ixFmax2
      end select
   
case(3)
      ixEmin3=1;ixEmax3=1;
      select case(iside)
      case(1)
        ixFmax3=ixFmin3
      case(2)
        ixFmin3=ixFmax3
      end select
   
    end select

    pflux(iside,idims,igrid)%edge(ixEmin1:ixEmax1,ixEmin2:ixEmax2,&
       ixEmin3:ixEmax3,idir1)=&
    fE(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,idir2)

  end if

end do

end subroutine fluxtoedge
!=============================================================================
subroutine fix_edges(psuse,idimmin,idimmax)
use mod_fix_conserve
use mod_amrvacdef

type(state) :: psuse(ngridshi)
integer, intent(in) :: idimmin,idimmax

integer :: iigrid, igrid, idims, iside, iotherside, i1,i2,i3, ic1,ic2,ic3,&
    inc1,inc2,inc3, ixMcmin1,ixMcmin2,ixMcmin3,ixMcmax1,ixMcmax2,ixMcmax3
integer :: nbuf, ibufnext
integer :: ibufnext_cc
integer :: pi1,pi2,pi3, mi1,mi2,mi3, ph1,ph2,ph3, mh1,mh2,mh3 !To detect corners
integer :: ixEmin1(1:ndir),ixEmin2(1:ndir),ixEmin3(1:ndir),ixEmax1(1:ndir),&
   ixEmax2(1:ndir),ixEmax3(1:ndir), ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,&
   ixtEmax2,ixtEmax3, ixFmin1(1:ndim),ixFmin2(1:ndim),ixFmin3(1:ndim),&
   ixFmax1(1:ndim),ixFmax2(1:ndim),ixFmax3(1:ndim), ixfEmin1(1:ndir),&
   ixfEmin2(1:ndir),ixfEmin3(1:ndir),ixfEmax1(1:ndir),ixfEmax2(1:ndir),&
   ixfEmax3(1:ndir)
integer :: nx1,nx2,nx3, idir, ix, ipe_neighbor, ineighbor
logical :: pcorner(1:ndim),mcorner(1:ndim)
!-----------------------------------------------------------------------------
if (nrecv_ct>0) then
   call MPI_WAITALL(nrecv_ct,recvrequest_stg,recvstatus_stg,ierrmpi)
end if

! Initialise buffer counter again
ibuf=1
ibuf_cc=1

do iigrid=1,igridstail; igrid=igrids(iigrid);
  do idims= idimmin,idimmax
    do iside=1,2
      i1=kr(1,idims)*(2*iside-3);i2=kr(2,idims)*(2*iside-3)
      i3=kr(3,idims)*(2*iside-3);
       if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle
      select case(neighbor_type(i1,i2,i3,igrid))
      case(4)
      ! The first neighbour is finer
      if (.not.neighbor_active(i1,i2,i3,igrid).or.&
          .not.neighbor_active(0,0,0,igrid) ) then
         do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
         inc3=2*i3+ic3
do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
         inc2=2*i2+ic2
do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
         inc1=2*i1+ic1
         ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
         !! When the neighbour is in a different process
         if (ipe_neighbor/=mype) then
            ibufnext=ibuf+isize(idims)
            ibuf=ibufnext
            end if
         end do
end do
end do
         cycle
      end if

      ! Check if there are corners
      pcorner=.false.
      mcorner=.false.
      do idir=1,ndim
        pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
        mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
        ph1=pi1-kr(idims,1)*(2*iside-3);ph2=pi2-kr(idims,2)*(2*iside-3)
        ph3=pi3-kr(idims,3)*(2*iside-3);
        mh1=mi1-kr(idims,1)*(2*iside-3);mh2=mi2-kr(idims,2)*(2*iside-3)
        mh3=mi3-kr(idims,3)*(2*iside-3);
        if (neighbor_type(ph1,ph2,ph3,igrid).eq.4) pcorner(idir)=.true.
        if (neighbor_type(mh1,mh2,mh3,igrid).eq.4) mcorner(idir)=.true.
      end do

      ! Calculate indices range
      call set_ix_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,&
         ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1,&
         ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,ixfEmin1,ixfEmin2,ixfEmin3,&
         ixfEmax1,ixfEmax2,ixfEmax3,igrid,idims,iside,.false.,.false.,0,0,0,&
         pcorner,mcorner)

      ! Remove coarse part of circulation
      call add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,&
         ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1,&
         ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,ixfEmin1,ixfEmin2,ixfEmin3,&
         ixfEmax1,ixfEmax2,ixfEmax3,pflux(iside,idims,igrid)%edge,idims,iside,&
         .false.,psuse(igrid))

      ! Add fine part of the circulation
      do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
         inc3=2*i3+ic3
do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
         inc2=2*i2+ic2
do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
         inc1=2*i1+ic1
         ineighbor=neighbor_child(1,inc1,inc2,inc3,igrid)
         ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
         iotherside=3-iside
         nx1=(ixMhi1-ixMlo1+1)/2;nx2=(ixMhi2-ixMlo2+1)/2
         nx3=(ixMhi3-ixMlo3+1)/2;

         call set_ix_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,&
            ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1,&
            ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,ixfEmin1,ixfEmin2,&
            ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,igrid,idims,iside,.true.,&
            .false.,inc1,inc2,inc3,pcorner,mcorner)
         if (ipe_neighbor.eq.mype) then
         call add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,&
            ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1,&
            ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,ixfEmin1,ixfEmin2,&
            ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,pflux(iotherside,idims,&
            ineighbor)%edge,idims,iside,.true.,psuse(igrid))

         else
         ibufnext=ibuf+isize(idims)
         call add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,&
            ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1,&
            ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,ixfEmin1,ixfEmin2,&
            ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,&
         reshape(source=recvbuffer(ibufnext-isize_stg(idims):ibufnext-1),&
         shape=(/ ixtEmax1-ixtEmin1+1,ixtEmax2-ixtEmin2+1,ixtEmax3-ixtEmin3+1 &
            ,3-1 /)),&
         idims,iside,.true.,psuse(igrid))

         ibuf=ibufnext

         end if
      end do
end do
end do

      case(3)
      ! The first neighbour is at the same level
      ! Check if there are corners
      do idir=idims+1,ndim
        pcorner=.false.
        mcorner=.false.
        pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
        mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
        ph1=pi1-kr(idims,1)*(2*iside-3);ph2=pi2-kr(idims,2)*(2*iside-3)
        ph3=pi3-kr(idims,3)*(2*iside-3);
        mh1=mi1-kr(idims,1)*(2*iside-3);mh2=mi2-kr(idims,2)*(2*iside-3)
        mh3=mi3-kr(idims,3)*(2*iside-3);
        if (neighbor_type(pi1,pi2,pi3,igrid).eq.4&
          .and.neighbor_type(ph1,ph2,ph3,igrid).eq.3&
          .and.neighbor_pole(pi1,pi2,pi3,igrid).eq.0) then
          pcorner(idir)=.true.
        ! Remove coarse part
        ! Set indices
          call set_ix_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,&
             ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1,&
             ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,ixfEmin1,ixfEmin2,&
             ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,igrid,idims,iside,.false.,&
             .true.,0,0,0,pcorner,mcorner)
        ! Remove
          call add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,&
             ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1,&
             ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,ixfEmin1,ixfEmin2,&
             ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,pflux(iside,idims,&
             igrid)%edge,idims,iside,.false.,psuse(igrid))

        ! Add fine part
        ! Find relative position of finer grid
          do ix=1,2

          inc1=kr(idims,1)*3*(iside-1)+3*kr(idir,1)+kr(6-idir-idims,1)*ix
          inc2=kr(idims,2)*3*(iside-1)+3*kr(idir,2)+kr(6-idir-idims,2)*ix
          inc3=kr(idims,3)*3*(iside-1)+3*kr(idir,3)+kr(6-idir-idims,3)*ix;
          ineighbor=neighbor_child(1,inc1,inc2,inc3,igrid)
          ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
          iotherside=3-iside

        ! Set indices
          call set_ix_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,&
             ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1,&
             ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,ixfEmin1,ixfEmin2,&
             ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,igrid,idims,iside,.true.,&
             .true.,inc1,inc2,inc3,pcorner,mcorner)

        ! add

          if (ipe_neighbor.eq.mype) then

          call add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,&
             ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1,&
             ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,ixfEmin1,ixfEmin2,&
             ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,pflux(iotherside,idims,&
             ineighbor)%edge,idims,iside,.true.,psuse(igrid))

          else


          ibufnext_cc=ibuf_cc+isize_stg(idims)
          call add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,&
             ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1,&
             ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,ixfEmin1,ixfEmin2,&
             ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,&
          reshape(source=recvbuffer_cc(ibuf_cc:ibufnext_cc-1),&
          shape=(/ ixtEmax1-ixtEmin1+1,ixtEmax2-ixtEmin2+1,&
             ixtEmax3-ixtEmin3+1 ,3-1 /)),&
          idims,iside,.true.,psuse(igrid))

          ibuf_cc=ibufnext_cc

          end if

          end do
        ! Set CoCorner to false again for next step

          pcorner(idir)=.false.

 !call div_staggered(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,psuse(igrid))
        end if

        if (neighbor_type(mi1,mi2,mi3,igrid).eq.4&
          .and.neighbor_type(mh1,mh2,mh3,igrid).eq.3&
          .and.neighbor_pole(mi1,mi2,mi3,igrid).eq.0) then
          mcorner(idir)=.true.
        ! Remove coarse part
        ! Set indices
          call set_ix_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,&
             ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1,&
             ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,ixfEmin1,ixfEmin2,&
             ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,igrid,idims,iside,.false.,&
             .true.,0,0,0,pcorner,mcorner)
        ! Remove
          call add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,&
             ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1,&
             ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,ixfEmin1,ixfEmin2,&
             ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,pflux(iside,idims,&
             igrid)%edge,idims,iside,.false.,psuse(igrid))

        ! Add fine part
        ! Find relative position of finer grid
          do ix=1,2

          inc1=kr(idims,1)*3*(iside-1)+kr(6-idir-idims,1)*ix
          inc2=kr(idims,2)*3*(iside-1)+kr(6-idir-idims,2)*ix
          inc3=kr(idims,3)*3*(iside-1)+kr(6-idir-idims,3)*ix;
          ineighbor=neighbor_child(1,inc1,inc2,inc3,igrid)
          ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
          iotherside=3-iside

        ! Set indices
          call set_ix_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,&
             ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1,&
             ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,ixfEmin1,ixfEmin2,&
             ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,igrid,idims,iside,.true.,&
             .true.,inc1,inc2,inc3,pcorner,mcorner)

        ! add

          if (ipe_neighbor.eq.mype) then

          call add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,&
             ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1,&
             ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,ixfEmin1,ixfEmin2,&
             ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,pflux(iotherside,idims,&
             ineighbor)%edge,idims,iside,.true.,psuse(igrid))

          else

          ibufnext_cc=ibuf_cc+isize_stg(idims)
          call add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,&
             ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1,&
             ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,ixfEmin1,ixfEmin2,&
             ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,&
          reshape(source=recvbuffer_cc(ibuf_cc:ibufnext_cc-1),&
          shape=(/ ixtEmax1-ixtEmin1+1,ixtEmax2-ixtEmin2+1,&
             ixtEmax3-ixtEmin3+1 ,3-1 /)),&
          idims,iside,.true.,psuse(igrid))

          ibuf_cc=ibufnext_cc

          end if

          end do
        ! Set CoCorner to false again for next step

         mcorner(idir)=.false.
        end if
      end do
      end select
    end do
  end do
 ! Average to centers
 call faces2centers(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,psuse(igrid)) 
end do

if (nrecv>0) deallocate(recvbuffer,recvstatus,recvrequest)
if (nrecv_ct>0) deallocate(recvbuffer_cc,recvstatus_stg,recvrequest_stg)

if (nsend>0) then
   call MPI_WAITALL(nsend,sendrequest,sendstatus,ierrmpi)
   deallocate(sendbuffer,sendstatus,sendrequest)
end if

if (nsend_ct>0) then
   call MPI_WAITALL(nsend_ct,sendrequest_stg,sendstatus_stg,ierrmpi)
   deallocate(sendbuffer_cc,sendstatus_stg,sendrequest_stg)
end if

end subroutine fix_edges
!=============================================================================
subroutine set_ix_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,&
   ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1,ixEmin2,&
   ixEmin3,ixEmax1,ixEmax2,ixEmax3,ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,&
   ixfEmax2,ixfEmax3,igrid,idims,iside,add,CoCorner,inc1,inc2,inc3,pcorner,&
   mcorner)
use mod_amrvacdef
! This routine sets the indexes for the correction
! of the circulation according to several different
! cases, as grids located in different cpus,
! presence of corners, and different relative locations
! of the fine grid respect to the coarse one

integer,intent(in)    :: igrid,idims,iside,inc1,inc2,inc3
logical,intent(in)    :: add,CoCorner
logical,intent(inout) :: pcorner(1:ndim),mcorner(1:ndim)
integer,intent(out)   :: ixFmin1(1:ndim),ixFmin2(1:ndim),ixFmin3(1:ndim),&
   ixFmax1(1:ndim),ixFmax2(1:ndim),ixFmax3(1:ndim),ixtEmin1,ixtEmin2,ixtEmin3,&
   ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1(1:ndir),ixEmin2(1:ndir),ixEmin3(1:ndir),&
   ixEmax1(1:ndir),ixEmax2(1:ndir),ixEmax3(1:ndir),ixfEmin1(1:ndir),&
   ixfEmin2(1:ndir),ixfEmin3(1:ndir),ixfEmax1(1:ndir),ixfEmax2(1:ndir),&
   ixfEmax3(1:ndir) !Indices for faces and edges
! ... local ...
integer               :: icor1,icor2,icor3,idim1,idir,nx1,nx2,nx3,middle1,&
   middle2,middle3
integer               :: ixtfEmin1,ixtfEmin2,ixtfEmin3,ixtfEmax1,ixtfEmax2,&
   ixtfEmax3
!-----------------------------------------------------------------------------
! ixF -> Indices for the _F_aces, and
! depends on the field component
! ixtE -> are the _t_otal range of the 'edge' array
! ixE -> are the ranges of the edge array,
! depending on the component
! ixfE -> are the ranges of the fE array (3D),
! and also depend on the component

! ... General ...
! Assign indices for the size of the E field array

ixtfEmin1=ixMlo1-1;ixtfEmin2=ixMlo2-1;ixtfEmin3=ixMlo3-1;
ixtfEmax1=ixMhi1;ixtfEmax2=ixMhi2;ixtfEmax3=ixMhi3;

if (add) then
nx1=(ixMhi1-ixMlo1+1)/2;nx2=(ixMhi2-ixMlo2+1)/2;nx3=(ixMhi3-ixMlo3+1)/2;
else
nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;nx3=ixMhi3-ixMlo3+1;
end if

do idim1=1,ndim
  ixtEmin1=0;ixtEmin2=0;ixtEmin3=0;
  ixtEmax1=nx1;ixtEmax2=nx2;ixtEmax3=nx3;
  select case(idims)
  case(1)
    ixtEmin1=1;ixtEmax1=1;
    if (iside.eq.1) ixtfEmax1=ixtfEmin1;
    if (iside.eq.2) ixtfEmin1=ixtfEmax1;
  
case(2)
    ixtEmin2=1;ixtEmax2=1;
    if (iside.eq.1) ixtfEmax2=ixtfEmin2;
    if (iside.eq.2) ixtfEmin2=ixtfEmax2;
  
case(3)
    ixtEmin3=1;ixtEmax3=1;
    if (iside.eq.1) ixtfEmax3=ixtfEmin3;
    if (iside.eq.2) ixtfEmin3=ixtfEmax3;
  
  end select
end do

! Assign indices, considering only the face
! (idims and iside)
do idim1=1,ndim
  ixFmin1(idim1)=ixMlo1-kr(idim1,1);ixFmin2(idim1)=ixMlo2-kr(idim1,2)
  ixFmin3(idim1)=ixMlo3-kr(idim1,3);
  ixFmax1(idim1)=ixMhi1;ixFmax2(idim1)=ixMhi2;ixFmax3(idim1)=ixMhi3;
  select case(idims)
  case(1)
     select case(iside)
     case(1)
     ixFmax1(idim1)=ixFmin1(idim1)
     case(2)
     ixFmin1(idim1)=ixFmax1(idim1)
     end select
  
case(2)
     select case(iside)
     case(1)
     ixFmax2(idim1)=ixFmin2(idim1)
     case(2)
     ixFmin2(idim1)=ixFmax2(idim1)
     end select
  
case(3)
     select case(iside)
     case(1)
     ixFmax3(idim1)=ixFmin3(idim1)
     case(2)
     ixFmin3(idim1)=ixFmax3(idim1)
     end select
  
  end select
end do

! ... Relative position ...
! Restrict range using relative position
if (add) then
middle1=(ixMhi1+ixMlo1)/2;middle2=(ixMhi2+ixMlo2)/2;middle3=(ixMhi3+ixMlo3)/2;

if (inc1.eq.1) then
ixFmax1(:)=middle1
ixtfEmax1=middle1
end if
if (inc1.eq.2) then
ixFmin1(:)=middle1+1
ixtfEmin1=middle1
end if


if (inc2.eq.1) then
ixFmax2(:)=middle2
ixtfEmax2=middle2
end if
if (inc2.eq.2) then
ixFmin2(:)=middle2+1
ixtfEmin2=middle2
end if


if (inc3.eq.1) then
ixFmax3(:)=middle3
ixtfEmax3=middle3
end if
if (inc3.eq.2) then
ixFmin3(:)=middle3+1
ixtfEmin3=middle3
end if

end if

! ... Adjust ranges of edges according to direction ...

do idim1=1,ndir
  ixfEmax1(idim1)=ixtfEmax1;ixfEmax2(idim1)=ixtfEmax2
  ixfEmax3(idim1)=ixtfEmax3;
  ixEmax1(idim1)=ixtEmax1;ixEmax2(idim1)=ixtEmax2;ixEmax3(idim1)=ixtEmax3;
  ixfEmin1(idim1)=ixtfEmin1+kr(idim1,1);ixfEmin2(idim1)=ixtfEmin2+kr(idim1,2)
  ixfEmin3(idim1)=ixtfEmin3+kr(idim1,3);
  ixEmin1(idim1)=ixtEmin1+kr(idim1,1);ixEmin2(idim1)=ixtEmin2+kr(idim1,2)
  ixEmin3(idim1)=ixtEmin3+kr(idim1,3);
end do



! ... Corners ...
! 'Coarse' corners
if (CoCorner) then
  do idim1=idims+1,ndim
    if (pcorner(idim1)) then
      do idir=1,ndir !Index arrays have size ndim
        if (idir.eq.6-idim1-idims) then
         !!! Something here has to change
         !!! Array ixfE must have size ndir, while
         !!! ixE must have size ndim
         if (1.eq.idim1) then
            ixfEmin1(idir)=ixfEmax1(idir)
            if (add) then
              ixEmax1(idir) =ixEmin1(idir)
            else
              ixEmin1(idir) =ixEmax1(idir)
            end if
          end if
if (2.eq.idim1) then
            ixfEmin2(idir)=ixfEmax2(idir)
            if (add) then
              ixEmax2(idir) =ixEmin2(idir)
            else
              ixEmin2(idir) =ixEmax2(idir)
            end if
          end if
if (3.eq.idim1) then
            ixfEmin3(idir)=ixfEmax3(idir)
            if (add) then
              ixEmax3(idir) =ixEmin3(idir)
            else
              ixEmin3(idir) =ixEmax3(idir)
            end if
          end if
        else
          ixEmin1(idir)=1;ixEmin2(idir)=1;ixEmin3(idir)=1;
          ixEmax1(idir)=0;ixEmax2(idir)=0;ixEmax3(idir)=0;
          ixfEmin1(idir)=1;ixfEmin2(idir)=1;ixfEmin3(idir)=1;
          ixfEmax1(idir)=0;ixfEmax2(idir)=0;ixfEmax3(idir)=0;
        end if
      end do
    end if
    if (mcorner(idim1)) then
      do idir=1,ndir
        if (idir.eq.6-idim1-idims) then
         if (1.eq.idim1) then
            ixfEmax1(idir)=ixfEmin1(idir)
            if (add) then
              ixEmin1(idir) =ixEmax1(idir)
            else
              ixEmax1(idir) =ixEmin1(idir)
            end if
          end if
if (2.eq.idim1) then
            ixfEmax2(idir)=ixfEmin2(idir)
            if (add) then
              ixEmin2(idir) =ixEmax2(idir)
            else
              ixEmax2(idir) =ixEmin2(idir)
            end if
          end if
if (3.eq.idim1) then
            ixfEmax3(idir)=ixfEmin3(idir)
            if (add) then
              ixEmin3(idir) =ixEmax3(idir)
            else
              ixEmax3(idir) =ixEmin3(idir)
            end if
          end if
        else
          ixEmin1(idir)=1;ixEmin2(idir)=1;ixEmin3(idir)=1;
          ixEmax1(idir)=0;ixEmax2(idir)=0;ixEmax3(idir)=0;
          ixfEmin1(idir)=1;ixfEmin2(idir)=1;ixfEmin3(idir)=1;
          ixfEmax1(idir)=0;ixfEmax2(idir)=0;ixfEmax3(idir)=0;
        end if
      end do
    end if
  end do
else
! Other kinds of corners
!! Crop ranges to account for corners
!! When the fine fluxes are added, we consider 
!! whether they come from the same cpu or from
!! a different one, in order to mimimise the 
!! amount of communication
!
!!!!Case for different processors still not implemented!!!
  
   if((idims.gt.1).and.pcorner(1)) then
      if((.not.add).or.(inc1.eq.2)) then
 !ixFmax1(:)=ixFmax1(:)-kr(1,1);!ixFmax2(:)=ixFmax2(:)-kr(1,2);!ixFmax3(:)=ixFmax3(:)-kr(1,3);
        do idir=1,ndir
          if ((idir.eq.idims).or.(idir.eq.1)) cycle
            ixfEmax1(idir)=ixfEmax1(idir)-1
            ixEmax1(idir)=ixEmax1(idir)-1
        end do
      end if
    end if
if((idims.gt.2).and.pcorner(2)) then
      if((.not.add).or.(inc2.eq.2)) then
 !ixFmax1(:)=ixFmax1(:)-kr(2,1);!ixFmax2(:)=ixFmax2(:)-kr(2,2);!ixFmax3(:)=ixFmax3(:)-kr(2,3);
        do idir=1,ndir
          if ((idir.eq.idims).or.(idir.eq.2)) cycle
            ixfEmax2(idir)=ixfEmax2(idir)-1
            ixEmax2(idir)=ixEmax2(idir)-1
        end do
      end if
    end if
if((idims.gt.3).and.pcorner(3)) then
      if((.not.add).or.(inc3.eq.2)) then
 !ixFmax1(:)=ixFmax1(:)-kr(3,1);!ixFmax2(:)=ixFmax2(:)-kr(3,2);!ixFmax3(:)=ixFmax3(:)-kr(3,3);
        do idir=1,ndir
          if ((idir.eq.idims).or.(idir.eq.3)) cycle
            ixfEmax3(idir)=ixfEmax3(idir)-1
            ixEmax3(idir)=ixEmax3(idir)-1
        end do
      end if
    end if
   if((idims.gt.1).and.mcorner(1)) then
      if((.not.add).or.(inc1.eq.1)) then
 !ixFmin1(:)=ixFmin1(:)+kr(1,1);!ixFmin2(:)=ixFmin2(:)+kr(1,2);!ixFmin3(:)=ixFmin3(:)+kr(1,3);
        do idir=1,ndir
          if ((idir.eq.idims).or.(idir.eq.1)) cycle
            ixfEmin1(idir)=ixfEmin1(idir)+1
            ixEmin1(idir)=ixEmin1(idir)+1
        end do
      end if
    end if
if((idims.gt.2).and.mcorner(2)) then
      if((.not.add).or.(inc2.eq.1)) then
 !ixFmin1(:)=ixFmin1(:)+kr(2,1);!ixFmin2(:)=ixFmin2(:)+kr(2,2);!ixFmin3(:)=ixFmin3(:)+kr(2,3);
        do idir=1,ndir
          if ((idir.eq.idims).or.(idir.eq.2)) cycle
            ixfEmin2(idir)=ixfEmin2(idir)+1
            ixEmin2(idir)=ixEmin2(idir)+1
        end do
      end if
    end if
if((idims.gt.3).and.mcorner(3)) then
      if((.not.add).or.(inc3.eq.1)) then
 !ixFmin1(:)=ixFmin1(:)+kr(3,1);!ixFmin2(:)=ixFmin2(:)+kr(3,2);!ixFmin3(:)=ixFmin3(:)+kr(3,3);
        do idir=1,ndir
          if ((idir.eq.idims).or.(idir.eq.3)) cycle
            ixfEmin3(idir)=ixfEmin3(idir)+1
            ixEmin3(idir)=ixEmin3(idir)+1
        end do
      end if
    end if
end if

end subroutine set_ix_circ
!=============================================================================
subroutine add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,&
   ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1,ixEmin2,&
   ixEmin3,ixEmax1,ixEmax2,ixEmax3,ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,&
   ixfEmax2,ixfEmax3,edge,idims,iside,add,s)
use mod_fix_conserve
use mod_amrvacdef

type(state)        :: s
integer,intent(in) :: idims,iside
integer            :: ixFmin1(1:ndim),ixFmin2(1:ndim),ixFmin3(1:ndim),&
   ixFmax1(1:ndim),ixFmax2(1:ndim),ixFmax3(1:ndim),ixtEmin1,ixtEmin2,ixtEmin3,&
   ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1(1:ndir),ixEmin2(1:ndir),ixEmin3(1:ndir),&
   ixEmax1(1:ndir),ixEmax2(1:ndir),ixEmax3(1:ndir),ixfEmin1(1:ndir),&
   ixfEmin2(1:ndir),ixfEmin3(1:ndir),ixfEmax1(1:ndir),ixfEmax2(1:ndir),&
   ixfEmax3(1:ndir)
double precision   :: edge(ixtEmin1:ixtEmax1,ixtEmin2:ixtEmax2,&
   ixtEmin3:ixtEmax3,1:ndim-1)
logical,intent(in) :: add
! ... local ...
integer            :: idim1,idim2,idir,middle1,middle2,middle3
integer            :: ixfECmin1,ixfECmin2,ixfECmin3,ixfECmax1,ixfECmax2,&
   ixfECmax3,ixECmin1,ixECmin2,ixECmin3,ixECmax1,ixECmax2,ixECmax3
double precision   :: fE(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:3) !!!!!!!!
double precision   :: circ(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:ndim) !!!!!!!!
integer            :: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,hxmin1,hxmin2,&
   hxmin3,hxmax1,hxmax2,hxmax3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,&
   ixCmax3,hxCmin1,hxCmin2,hxCmin3,hxCmax1,hxCmax2,hxCmax3 !Indices for edges
!-----------------------------------------------------------------------------
! ixF -> Indices for the faces, depends on the field component
! ixE -> Total range for the edges
! ixfE -> Edges in fE (3D) array
associate(mygeo=>s%geo,bfaces=>s%ws%w)
! ix,hx,ixC,hxC -> Auxiliary indices

! Assign quantities stored ad edges to make it as similar as 
! possible to the routine updatefaces.

fE(:,:,:,:)=zero

do idim1=1,ndim-1
   ! 3D: rotate indices (see routine fluxtoedge)
  idir=mod(idim1+idims-1,3)+1
  
  ixfECmin1=ixfEmin1(idir);ixfECmin2=ixfEmin2(idir);ixfECmin3=ixfEmin3(idir)
  ixfECmax1=ixfEmax1(idir);ixfECmax2=ixfEmax2(idir);ixfECmax3=ixfEmax3(idir);
  ixECmin1=ixEmin1(idir);ixECmin2=ixEmin2(idir);ixECmin3=ixEmin3(idir)
  ixECmax1=ixEmax1(idir);ixECmax2=ixEmax2(idir);ixECmax3=ixEmax3(idir);
  fE(ixfECmin1:ixfECmax1,ixfECmin2:ixfECmax2,ixfECmin3:ixfECmax3,idir)&
     =edge(ixECmin1:ixECmax1,ixECmin2:ixECmax2,ixECmin3:ixECmax3,idim1)

end do

! Calculate part of circulation needed
circ=zero
do idim1=1,ndim
   do idim2=1,ndim
      do idir=1,ndir
        if (lvc(idim1,idim2,idir).eq.0) cycle
        ! Assemble indices
        ixCmin1=ixFmin1(idim1);ixCmin2=ixFmin2(idim1);ixCmin3=ixFmin3(idim1)
        ixCmax1=ixFmax1(idim1);ixCmax2=ixFmax2(idim1);ixCmax3=ixFmax3(idim1);
        hxCmin1=ixCmin1-kr(idim2,1);hxCmin2=ixCmin2-kr(idim2,2)
        hxCmin3=ixCmin3-kr(idim2,3);hxCmax1=ixCmax1-kr(idim2,1)
        hxCmax2=ixCmax2-kr(idim2,2);hxCmax3=ixCmax3-kr(idim2,3);
        if (idim1.eq.idims) then
        circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
           =circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
                         +lvc(idim1,idim2,idir)&
                         *(fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                            idir)&
                         -fE(hxCmin1:hxCmax1,hxCmin2:hxCmax2,hxCmin3:hxCmax3,&
                            idir))

        else
         
          select case(iside)
          case(2)
          circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
             =circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
                         +lvc(idim1,idim2,idir)&
                         *fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                            idir)

          case(1)
          circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
             =circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
                         -lvc(idim1,idim2,idir)&
                         *fE(hxCmin1:hxCmax1,hxCmin2:hxCmax2,hxCmin3:hxCmax3,&
                            idir)

          end select
        end if
      end do
   end do
end do

! Divide circulation by surface and add
do idim1=1,ndim
   ixCmin1=ixFmin1(idim1);ixCmin2=ixFmin2(idim1);ixCmin3=ixFmin3(idim1)
   ixCmax1=ixFmax1(idim1);ixCmax2=ixFmax2(idim1);ixCmax3=ixFmax3(idim1);
   select case(idim1)
   case(1)
   where(mygeo%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
      .gt. 1.0d-9*mygeo%dvolume(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      ixCmin3:ixCmax3))
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
         =circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
                 /mygeo%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                    ixCmin3:ixCmax3)
   elsewhere
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)=zero
   end where
   
case(2)
   where(mygeo%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
      .gt. 1.0d-9*mygeo%dvolume(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      ixCmin3:ixCmax3))
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
         =circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
                 /mygeo%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                    ixCmin3:ixCmax3)
   elsewhere
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)=zero
   end where
   
case(3)
   where(mygeo%surfaceC3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
      .gt. 1.0d-9*mygeo%dvolume(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      ixCmin3:ixCmax3))
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
         =circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
                 /mygeo%surfaceC3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                    ixCmin3:ixCmax3)
   elsewhere
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)=zero
   end where
   
    end select
   ! Add/subtract to field at face

   if (add) then
      bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
         =bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1)-circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)
   else
      bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
         =bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1)+circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)
   end if

end do

end associate

end subroutine add_sub_circ

!=============================================================================


subroutine b_from_vectorpotential(ixIsmin1,ixIsmin2,ixIsmin3,ixIsmax1,&
   ixIsmax2,ixIsmax3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
    ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, ws, x)

use mod_amrvacdef

integer, intent(in)                :: ixIsmin1,ixIsmin2,ixIsmin3,ixIsmax1,&
   ixIsmax2,ixIsmax3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
    ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(inout)    :: ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,&
   ixIsmin3:ixIsmax3,1:nws)
double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndim)
! ... local ...
double precision                   :: Adummy(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndir)
!----------------------------------------------------------------

call b_from_vectorpotentialA(ixIsmin1,ixIsmin2,ixIsmin3,ixIsmax1,ixIsmax2,&
   ixIsmax3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,&
   ixOmin3,ixOmax1,ixOmax2,ixOmax3, ws, x, Adummy)

end subroutine b_from_vectorpotential
!=============================================================================
subroutine b_from_vectorpotentialA(ixIsmin1,ixIsmin2,ixIsmin3,ixIsmax1,&
   ixIsmax2,ixIsmax3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
    ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, ws, x, A)

use mod_amrvacdef

integer, intent(in)                :: ixIsmin1,ixIsmin2,ixIsmin3,ixIsmax1,&
   ixIsmax2,ixIsmax3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
    ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(inout)    :: ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,&
   ixIsmin3:ixIsmax3,1:nws),A(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:ndir)
double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndim)
! .. local ..
integer                            :: ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,&
   ixCmax3, hxCmin1,hxCmin2,hxCmin3,hxCmax1,hxCmax2,hxCmax3, hxOmin1,hxOmin2,&
   hxOmin3,hxOmax1,hxOmax2,hxOmax3, idim, idim1, idim2, idir
double precision                   :: xC(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndim)
double precision                   :: circ(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndim), dxidir
double precision, parameter        :: smallsurface=1.0d-9
!-----------------------------------------------------------------------------

A(:,:,:,:)=zero
ws(:,:,:,:)=zero

ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
ixCmin1=ixOmin1-1;ixCmin2=ixOmin2-1;ixCmin3=ixOmin3-1; ! Extend range by one
do idir=7-2*ndim,ndir

   do idim=1,ndim
      ! Get edge coordinates
      if (idim.ne.idir) then
         xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim) &
            = x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
            idim) + half * dxlevel(idim)
      else
         xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim) &
            = x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim)
      end if
   end do

   ! Initialise vector potencial at the edge
   call initvecpot_usr(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, xC, A(ixImin1:ixImax1,&
      ixImin2:ixImax2,ixImin3:ixImax3,idir), idir)

end do

! Set NaN to zero (can happen e.g. on axis):
where(A(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1:ndir).ne.A&
   (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1:ndir))
   A(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1:ndir)=zero
end where


! --------------------------------------------------
! Take the curl of the vector potential 
! --------------------------------------------------
circ(:,:,:,:) = zero

! Calculate circulation on each face
do idim1=1,ndim ! Coordinate perpendicular to face 
   ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
   ixCmin1=ixOmin1-kr(idim1,1);ixCmin2=ixOmin2-kr(idim1,2)
   ixCmin3=ixOmin3-kr(idim1,3);
   do idim2=1,ndim
      do idir=1,ndir ! Direction of line integral
        
        ! Assemble indices
        hxCmin1=ixCmin1-kr(idim2,1);hxCmin2=ixCmin2-kr(idim2,2)
        hxCmin3=ixCmin3-kr(idim2,3);hxCmax1=ixCmax1-kr(idim2,1)
        hxCmax2=ixCmax2-kr(idim2,2);hxCmax3=ixCmax3-kr(idim2,3);
        ! Add line integrals in direction idir
        if (idir .le. ndim) then
                dxidir = dxlevel(idir)
        else
                dxidir = 1.0d0
        end if
        circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
           =circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
                         +lvc(idim1,idim2,idir)*dxidir &
                         *(A(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                            idir)&
                         -A(hxCmin1:hxCmax1,hxCmin2:hxCmax2,hxCmin3:hxCmax3,&
                            idir))
      end do
   end do
end do

! Set NaN to zero (should not occur)
where(circ.ne.circ)
   circ=zero
end where


! --------------------------------------------------
! Divide by the area of the face to get B
! --------------------------------------------------
do idim1=1,ndim
   ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
   ixCmin1=ixOmin1-kr(idim1,1);ixCmin2=ixOmin2-kr(idim1,2)
   ixCmin3=ixOmin3-kr(idim1,3);
   select case(idim1)
   case(1)
   where(mygeo%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
      .gt. smallsurface)
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
         =circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
                 /mygeo%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                    ixCmin3:ixCmax3)
   elsewhere
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)=zero
   end where
   
case(2)
   where(mygeo%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
      .gt. smallsurface)
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
         =circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
                 /mygeo%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                    ixCmin3:ixCmax3)
   elsewhere
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)=zero
   end where
   
case(3)
   where(mygeo%surfaceC3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
      .gt. smallsurface)
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
         =circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)&
                 /mygeo%surfaceC3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                    ixCmin3:ixCmax3)
   elsewhere
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)=zero
   end where
   
   end select
   ws(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1) &
      = circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)
end do

end subroutine b_from_vectorpotentialA
!=============================================================================

subroutine printarray(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,array)
! Routine for easily printing an array, just for debugging, erase later

use mod_amrvacdef

integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(in)  :: array(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   ixOmin3:ixOmax3)
!double precision, intent(in)  :: array(:,:)
integer                       :: i,j,k
!----------------------------------------------------------------------------

print *, array


   do i=ixOmin1,ixOmax1




      do k=ixOmin3,ixOmax3
        write (*,'(I2,I2,I2,E17.9,$)') i,j,k,array(i,j,k)
      end do
      print *, ' '
     
   end do

   print *, '---'


write (*,'(A)') '------------------------------------------------------------'
end subroutine printarray
!=============================================================================
