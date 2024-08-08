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
subroutine upwindLR(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixLmin1,&
   ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,&
   ixRmax2,ixRmax3,idims,w,wCT,wLC,wRC,wLp,wRp,x,dxdim)

! Determine the upwinded wLC(ixL) and wRC(ixR) from w. 
! the wCT is only used when PPM is exploited.

use mod_limiter
use mod_amrvacdef

integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
    ixLmin1,ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3, ixRmin1,ixRmin2,ixRmin3,&
   ixRmax1,ixRmax2,ixRmax3, idims
double precision, intent(in) :: dxdim
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nw) :: w, wCT
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nw) :: wLC, wRC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nw) :: wLp, wRp
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:ndim) :: x
! .. local ..
integer :: jxRmin1,jxRmin2,jxRmin3,jxRmax1,jxRmax2,jxRmax3, ixCmin1,ixCmin2,&
   ixCmin3,ixCmax1,ixCmax2,ixCmax3, jxCmin1,jxCmin2,jxCmin3,jxCmax1,jxCmax2,&
   jxCmax3, iw
double precision :: ldw(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
    rdw(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3), dwC(ixImin1:ixImax1,&
   ixImin2:ixImax2,ixImin3:ixImax3)
!-----------------------------------------------------------------------------

select case (typelimiter)

case (limiter_wenoZP)
   call WENOZPlimiter(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixLmin1,&
      ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3,idims,dxdim,w,wLp,wRp)
   
case (limiter_weno5)
   call WENO5limiter(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixLmin1,&
      ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3,idims,w,wLp,wRp)
   
case (limiter_ppm)
   ! Since PPM uses the ordinary grid-index:
   ixCmin1=ixLmin1+kr(1,idims);ixCmin2=ixLmin2+kr(2,idims)
   ixCmin3=ixLmin3+kr(3,idims);
   ixCmax1=ixLmax1;ixCmax2=ixLmax2;ixCmax3=ixLmax3;
   call PPMlimiter(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,&
      ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,w,wCT,wLp,wRp)
   
case (limiter_mp5)
   call MP5limiter(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixLmin1,&
      ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3,idims,w,wLp,wRp)
   
case default
   jxRmin1=ixRmin1+kr(idims,1);jxRmin2=ixRmin2+kr(idims,2)
   jxRmin3=ixRmin3+kr(idims,3);jxRmax1=ixRmax1+kr(idims,1)
   jxRmax2=ixRmax2+kr(idims,2);jxRmax3=ixRmax3+kr(idims,3);
   ixCmax1=jxRmax1;ixCmax2=jxRmax2;ixCmax3=jxRmax3
   ixCmin1=ixLmin1-kr(idims,1);ixCmin2=ixLmin2-kr(idims,2)
   ixCmin3=ixLmin3-kr(idims,3);
   jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
   jxCmin3=ixCmin3+kr(idims,3);jxCmax1=ixCmax1+kr(idims,1)
   jxCmax2=ixCmax2+kr(idims,2);jxCmax3=ixCmax3+kr(idims,3);

   do iw=1,nwflux

      dwC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=w(jxCmin1:jxCmax1,&
         jxCmin2:jxCmax2,jxCmin3:jxCmax3,iw)-w(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)

      call dwlimiter2(dwC,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,typelimiter,&
         ldw,rdw)
      wLp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3,iw)&
         =wLp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3,&
         iw)+half*ldw(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3)
      wRp(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,iw)&
         =wRp(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
         iw)-half*rdw(jxRmin1:jxRmax1,jxRmin2:jxRmax2,jxRmin3:jxRmax3)
      
   end do

end select


! Apply floor model to reconstructed state:
if (tlow>zero) then
   call fixp_usr(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixLmin1,&
      ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3,wLp,x)
   call fixp_usr(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixRmin1,&
      ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,wRp,x)
end if

! Also get the conserved state and auxiliaries
wLC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3,1:nwflux) &
   = wLp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3,1:nwflux)
wRC(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,1:nwflux) &
   = wRp(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,1:nwflux)
call conserve(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixLmin1,ixLmin2,&
   ixLmin3,ixLmax1,ixLmax2,ixLmax3,wLC,x,patchfalse)
call conserve(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixRmin1,ixRmin2,&
   ixRmin3,ixRmax1,ixRmax2,ixRmax3,wRC,x,patchfalse)
if (nwaux .gt. 0) then
   wLp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3,nwflux+1:nwflux+nwaux) &
      = wLC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3,&
      nwflux+1:nwflux+nwaux)
   wRp(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,nwflux+1:nwflux+nwaux) &
      = wRC(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
      nwflux+1:nwflux+nwaux)
end if

end subroutine upwindLR
!============================================================================
subroutine tvdlf(method,qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idimmin,idimmax, qtC,sCT,&
   qt,snew,sold,wprim,fC,fE,dx1,dx2,dx3,x)

! method=='tvdmu'  --> 2nd order (may be 3rd order in 1D) TVD-MUSCL scheme.
! method=='tvdmu1' --> 1st order TVD-MUSCL scheme (upwind per charact. var.)
! method=='tvdlf'  --> 2nd order TVD-Lax-Friedrich scheme.
! method=='tvdlf1' --> 1st order TVD-Lax-Friedrich scheme.

use mod_amrvacdef

character(len=*), intent(in)                           :: method
double precision, intent(in)                           :: qdt, qtC, qt, dx1,&
   dx2,dx3
integer, intent(in)                                    :: ixImin1,ixImin2,&
   ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
   ixOmax3, idimmin,idimmax
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:ndim), intent(in)  :: x
type(state)                                            :: sCT, snew, sold
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nw), intent(in)    :: wprim
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nwflux,1:ndim)     :: fC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:3)               :: fE

double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:ndim)              :: xi
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nw)     :: wLC, wRC, wLp, wRp
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nwflux) :: fLC, fRC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)          :: fadd, fdiss
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)          :: vLC, vRC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)          :: cmaxC, cmaxRC, cmaxLC

! UCT related storage
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:ndir,2) :: vbarC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:ndir,2) :: vbarLC,vbarRC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   ndim)   :: cbarmin,cbarmax
integer                                   :: idimE, idimN

double precision :: dxinv(1:ndim),dxdim(1:ndim)
integer :: idims, iw, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3, hxOmin1,&
   hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3, ixCmin1,ixCmin2,ixCmin3,ixCmax1,&
   ixCmax2,ixCmax3, ixCRmin1,ixCRmin2,ixCRmin3,ixCRmax1,ixCRmax2,ixCRmax3,&
    kxCmin1,kxCmin2,kxCmin3,kxCmax1,kxCmax2,kxCmax3, kxRmin1,kxRmin2,kxRmin3,&
   kxRmax1,kxRmax2,kxRmax3, ixtestmin1,ixtestmin2,ixtestmin3,ixtestmax1,&
   ixtestmax2,ixtestmax3
integer :: ixIsmin1,ixIsmin2,ixIsmin3,ixIsmax1,ixIsmax2,ixIsmax3
logical :: transport(1:nwflux)
!-----------------------------------------------------------------------------
associate(wCT=>sCT%w%w, wnew=>snew%w%w, wold=>sold%w%w)

  

ixIsmin1=sCT%ws%ixGmin1;ixIsmin2=sCT%ws%ixGmin2;ixIsmin3=sCT%ws%ixGmin3;
ixIsmax1=sCT%ws%ixGmax1;ixIsmax2=sCT%ws%ixGmax2;ixIsmax3=sCT%ws%ixGmax3;

staggered : associate(wCTs=>sCT%ws%w, wnews=>snew%ws%w, wolds=>sold%ws%w)



if (covariant .and. typelimited .ne. 'predictor') call mpistop&
('tvdlf: Covariant formulation only implemented with typelimited=predictor')


if (idimmax>idimmin .and. typelimited=='original' .and. method&
   /='tvdlf1')call mpistop("Error in TVDMUSCLF: Unsplit dim. and original is &
   limited")

! The flux calculation contracts by one in the idim direction it is applied.
! The limiter contracts the same directions by one more, so expand ixO by 2.
ixmin1=ixOmin1;ixmin2=ixOmin2;ixmin3=ixOmin3;ixmax1=ixOmax1;ixmax2=ixOmax2
ixmax3=ixOmax3;
do idims= idimmin,idimmax
   ixmin1=ixmin1-2*kr(idims,1);ixmin2=ixmin2-2*kr(idims,2)
   ixmin3=ixmin3-2*kr(idims,3);ixmax1=ixmax1+2*kr(idims,1)
   ixmax2=ixmax2+2*kr(idims,2);ixmax3=ixmax3+2*kr(idims,3);
end do
if (ixImin1>ixmin1.or.ixImin2>ixmin2.or.ixImin3>ixmin3.or.ixImax1<ixmax1&
   .or.ixImax2<ixmax2.or.ixImax3<ixmax3) call mpistop("Error in TVDLF: &
   Nonconforming input limits")


dxinv(1)=-qdt/dx1;dxinv(2)=-qdt/dx2;dxinv(3)=-qdt/dx3;
dxdim(1)=dx1;dxdim(2)=dx2;dxdim(3)=dx3;

do idims= idimmin,idimmax
   
   if (B0field) then
      select case (idims)
      case (1)
         myB0 => myB0_face1
      case (2)
         myB0 => myB0_face2
      case (3)
         myB0 => myB0_face3
      end select
   end if

   if (covariant) then
      select case (idims)
      case (1)
         myM => mygeo%mSurface1
      case (2)
         myM => mygeo%mSurface2
      case (3)
         myM => mygeo%mSurface3
      end select
   end if

   ! ================================
   ! ====== Assemble indices: =======  
   ! ================================
   hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
   hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
   hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);
   kxCmin1=ixImin1;kxCmin2=ixImin2;kxCmin3=ixImin3
   kxCmax1=ixImax1-kr(idims,1);kxCmax2=ixImax2-kr(idims,2)
   kxCmax3=ixImax3-kr(idims,3);
   kxRmin1=kxCmin1+kr(idims,1);kxRmin2=kxCmin2+kr(idims,2)
   kxRmin3=kxCmin3+kr(idims,3);kxRmax1=kxCmax1+kr(idims,1)
   kxRmax2=kxCmax2+kr(idims,2);kxRmax3=kxCmax3+kr(idims,3);

   select case (typeemf)
   case default
      ! ixC is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3; ixCmin1=hxOmin1
      ixCmin2=hxOmin2;ixCmin3=hxOmin3;
   case ('average')
      ! Flux-interpolated constrained transport needs one more layer in the tangential components:
      ixCmax1=ixOmax1+1-kr(idims,1);ixCmax2=ixOmax2+1-kr(idims,2)
      ixCmax3=ixOmax3+1-kr(idims,3); ixCmin1=hxOmin1-1+kr(idims,1)
      ixCmin2=hxOmin2-1+kr(idims,2);ixCmin3=hxOmin3-1+kr(idims,3);
   case ('uct1','uct2','uct2+av')
      ! Upwind constrained transport needs dixB more layers in the tangential components
      ixCmax1=ixOmax1+dixB-dixB*kr(idims,1)
      ixCmax2=ixOmax2+dixB-dixB*kr(idims,2)
      ixCmax3=ixOmax3+dixB-dixB*kr(idims,3)
      ixCmin1=hxOmin1-dixB+dixB*kr(idims,1)
      ixCmin2=hxOmin2-dixB+dixB*kr(idims,2)
      ixCmin3=hxOmin3-dixB+dixB*kr(idims,3);
   end select
   ixCRmin1=ixCmin1;ixCRmin2=ixCmin2;ixCRmin3=ixCmin3;ixCRmax1=ixCmax1
   ixCRmax2=ixCmax2;ixCRmax3=ixCmax3; !indices to reconstruct to
   ! ================================
   ! ====== Done with indices ======
   ! ================================
   
   ! Get interface positions:
   xi(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:ndim) &
      = x(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:ndim)
   xi(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,idims) &
      = half* ( x(kxRmin1:kxRmax1,kxRmin2:kxRmax2,kxRmin3:kxRmax3,&
      idims)+x(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,idims) )
   
   ! ================================
   ! === Start with flat (Godunov): =
   ! ================================
   wRp(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:nwflux+nwaux)&
      =wprim(kxRmin1:kxRmax1,kxRmin2:kxRmax2,kxRmin3:kxRmax3,1:nwflux+nwaux)
   wLp(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:nwflux+nwaux)&
      =wprim(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:nwflux+nwaux)
   ! ================================

   ! ================================
   ! for tvdlf (second order scheme): apply limiting
   ! ================================
   if (method=='tvdlf') then
      select case (typelimited)
      case ('previous')
         call primitive(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,wold,x)
         call upwindLR(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixCRmin1,ixCRmin2,ixCRmin3,ixCRmax1,ixCRmax2,ixCRmax3,ixCRmin1,&
            ixCRmin2,ixCRmin3,ixCRmax1,ixCRmax2,ixCRmax3,idims,wold,wprim,wLC,&
            wRC,wLp,wRp,xi,dxdim(idims))
         call conserve(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,wold,x,patchfalse)
         
         wLC(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,b0_+idims)&
            =wolds(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,&
            idims)
         wRC(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,b0_+idims)&
            =wolds(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,&
            idims)
         wLp(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,b0_+idims)&
            =wolds(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,&
            idims)
         wRp(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,b0_+idims)&
            =wolds(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,&
            idims)
         
      case ('predictor')
         call upwindLR(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixCRmin1,ixCRmin2,ixCRmin3,ixCRmax1,ixCRmax2,ixCRmax3,ixCRmin1,&
            ixCRmin2,ixCRmin3,ixCRmax1,ixCRmax2,ixCRmax3,idims,wprim,wprim,&
            wLC,wRC,wLp,wRp,xi,dxdim(idims))
         
         wLC(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,b0_+idims)&
            =wCTs(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,&
            idims)
         wRC(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,b0_+idims)&
            =wCTs(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,idims)
         wLp(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,b0_+idims)&
            =wCTs(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,&
            idims)
         wRp(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,b0_+idims)&
            =wCTs(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,idims)
         
      case ('original')
         call primitive(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,wnew,x)
         call upwindLR(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixCRmin1,ixCRmin2,ixCRmin3,ixCRmax1,ixCRmax2,ixCRmax3,ixCRmin1,&
            ixCRmin2,ixCRmin3,ixCRmax1,ixCRmax2,ixCRmax3,idims,wnew,wprim,wLC,&
            wRC,wLp,wRp,xi,dxdim(idims))
         call conserve(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,wnew,x,patchfalse)
         
         wLC(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,b0_+idims)&
            =wnews(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,&
            idims)
         wRC(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,b0_+idims)&
            =wnews(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,&
            idims)
         wLp(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,b0_+idims)&
            =wnews(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,&
            idims)
         wRp(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,b0_+idims)&
            =wnews(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,&
            idims)
         
      case default
         call mpistop("Error in TVDMUSCLF: no such base for limiter")
      end select
   else
      wRC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:nwflux+nwaux) &
         = wCT(kxRmin1:kxRmax1,kxRmin2:kxRmax2,kxRmin3:kxRmax3,1:nwflux+nwaux)
      wLC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:nwflux+nwaux) &
         = wCT(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:nwflux+nwaux)
   end if
   ! ==== Done with limiting ========

   
   ! ================================
   ! For the high order Lax-Friedrich TVDLF scheme the limiter is based on
   ! the maximum eigenvalue, it is calculated in advance.
   ! ================================
   call getcmax(wLC,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,cmaxLC,cmaxC,&
      .true.)
   cmaxLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
      = max(cmaxLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3),&
      dabs(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)))
   call getcmax(wRC,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,cmaxRC,cmaxC,&
      .true.)
   cmaxRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
      = max(cmaxRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3),&
      dabs(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)))
   ! now take the maximum of left and right states
   cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=max(cmaxRC&
      (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3),cmaxLC&
      (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3))

   ! ================================

   
   ! ================================
   ! Calculate velocities for transport fluxes
   ! ================================
   call getv(wLp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,&
      ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,vLC)
   call getv(wRp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,&
      ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,vRC)
   ! ================================

   

   ! ================================
   ! Storage for various ways of EMF averaging:
   ! ================================
   select case (typeemf)
   case default
   case ('average')
   case ('uct1')
      ! Store magnitude of characteristics
      cbarmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idims)&
         =max(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3),zero)
      cbarmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idims)&
         =cbarmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idims)
      
      
 !If the current direction will be useful for a component of the electric field
      ! (1 on 2D or 3 in 3D), then store relevant velocities.
      if (idims.le.2*ndim-3) then
         idimN=mod(idims,ndir)+1 ! 'Next' direction
         idimE=mod(idims+1,ndir)+1 ! Electric field direction

         vbarLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idimE,1)&
            =vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
         vbarRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idimE,1)&
            =vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
         ! Calculate additional component of velocity needed
         call getv(wLp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idimN,&
            vbarLC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idimE,2))
         call getv(wRp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idimN,&
            vbarRC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idimE,2))
      end if

   case ('uct2','uct2+av')
      ! Store magnitude of characteristics
      cbarmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idims)&
         =max(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3),zero)
      cbarmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idims)&
         =cbarmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idims)


      idimN=mod(idims,ndir)+1 ! 'Next' direction
      idimE=mod(idims+1,ndir)+1 ! Electric field direction

      call getv(wLp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idimN,&
         vbarLC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idims,1))
      call getv(wRp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idimN,&
         vbarRC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idims,1))

      vbarC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idims,1)&
         =half*(vbarLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idims,&
         1) &
           +vbarRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idims,1))

      call getv(wLp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idimE,&
         vbarLC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idims,2))
      call getv(wRp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idimE,&
         vbarRC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idims,2))

      vbarC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idims,2)&
         =half*(vbarLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idims,&
         2) &
           +vbarRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idims,2))

   end select
   ! ================================

   
   ! ================================
   ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
   ! ================================
      
   call getflux(wLC,wLp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,fLC,transport)
   call getflux(wRC,wRp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,fRC,transport)
   do iw=1,nwflux
      
      if (transport(iw)) then
         fadd(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)  &
            = fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
            iw)+vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)*wLC&
            (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
         fdiss(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
            = fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
            iw)+vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)*wRC&
            (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
      else
         fadd(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)  &
            = fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
         fdiss(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
            = fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
      end if
      
      ! Mean flux:
      fadd(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
         =half*(fadd(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)+fdiss&
         (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3))
      ! Add TVDLF dissipation:
      fdiss(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
         =-tvdlfeps*cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
         *half*(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw))
      ! fLC contains physical+dissipative fluxes
      fadd(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
         =fadd(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)+fdiss&
         (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)

      if (slab) then
         fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,idims)&
            =fadd(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
      else
         select case (idims)
            case (1)
            
 !This where statement catches the axis where inverse metric becomes inf:
            where(mygeo%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3) .gt. 1.0d-9*mygeo%dvolume(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,ixCmin3:ixCmax3))
               fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,1)&
                  =mygeo%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3)*fadd(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3)
            elsewhere
               fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,1) = zero
            end where
            
            
            case (2)
            
 !This where statement catches the axis where inverse metric becomes inf:
            where(mygeo%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3) .gt. 1.0d-9*mygeo%dvolume(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,ixCmin3:ixCmax3))
               fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,2)&
                  =mygeo%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3)*fadd(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3)
            elsewhere
               fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,2) = zero
            end where
            
            
            case (3)
            
 !This where statement catches the axis where inverse metric becomes inf:
            where(mygeo%surfaceC3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3) .gt. 1.0d-9*mygeo%dvolume(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,ixCmin3:ixCmax3))
               fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,3)&
                  =mygeo%surfaceC3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3)*fadd(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3)
            elsewhere
               fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,3) = zero
            end where
            
            
         end select
      end if

   end do ! Next iw
   ! ================================
end do ! Next idims


! Update values at faces
select case (typeemf)
case default
   call mpistop('please specify typeemf for staggered field upwinding')
case ('average')
   call updatefaces(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
      ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qdt,fC,fE,snew)
case ('uct1')
   call updatefacesuct1(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qdt,vbarRC,vbarLC,&
      cbarmin,cbarmax,fE,snew)
case ('uct2')
   call updatefacesuct2(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qdt,vbarC,cbarmin,&
      cbarmax,fE,snew)
case ('uct2+av')
   call updatefacesuct2av(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qdt,vbarC,cbarmin,&
      cbarmax,fC,fE,snew)
end select



! ================================
! Now update the state:
! ================================
do idims= idimmin,idimmax
   
   hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
   hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
   hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);
   
   do iw=1,nwflux

      ! Multiply the fluxes by -dt/dx since Flux fixing expects this
      if (slab) then
         fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw,idims)&
            =dxinv(idims)*fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
            iw,idims)
         wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
            =wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            iw) + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw,&
            idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,iw,&
            idims))
      else
         select case (idims)
         case (1)
            fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw,1)&
               =-qdt*fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw,&
               idims)
           wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
              =wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) &
              + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw,&
                 1)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,iw,1))&
                 /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3)
         case (2)
            fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw,2)&
               =-qdt*fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw,&
               idims)
           wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
              =wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) &
              + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw,&
                 2)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,iw,2))&
                 /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3)
         case (3)
            fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw,3)&
               =-qdt*fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw,&
               idims)
           wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
              =wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) &
              + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw,&
                 3)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,iw,3))&
                 /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3)
         end select
      end if

   end do ! Next iw

end do ! Next idims



call faces2centers(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,snew)

! ================================
! === Done updating state ========
! ================================


if (covariant) myM => mygeo%m
if (.not.slab.and.idimmin==1) then 
   call addgeometry(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,wnew,x)
else
   if(nwaux>0) call getaux(.true.,wnew,x,ixImin1,ixImin2,ixImin3,ixImax1,&
      ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
      'tvdlf_wnew')
end if

call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim),ixImin1,ixImin2,&
   ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
   ixOmax3,1,nw,qtC,wCT,qt,wnew,x,.false.)


end associate staggered


end associate
end subroutine tvdlf
!=============================================================================
subroutine tvdlfpos(method,qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idimmin,idimmax,&
    qtC,sCT,qt,snew,sold,wprim,fC,fE,dx1,dx2,dx3,x)
! method=='tvdlfpos'  --> 2nd order TVD-Lax-Friedrich scheme, positivity preserving on density


use mod_amrvacdef

character(len=*), intent(in)                          :: method
double precision, intent(in)                          :: qdt, qtC, qt, dx1,&
   dx2,dx3
integer, intent(in)                                   :: ixImin1,ixImin2,&
   ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
   ixOmax3, idimmin,idimmax
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:ndim), intent(in) ::  x
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:ndim)             ::  xi
type(state)                                           :: sCT, snew, sold
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nw), intent(in)   :: wprim
! .. local ..
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nwflux,1:ndim)        :: fC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:3)                  :: fE
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nwflux,1:ndim)        :: fClow
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)                        :: theta, theta_tau
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nw)     :: wLC, wRC, wLp, wRp
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nwflux) :: fLC, fRC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)          :: vLC, vRC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)          :: cmaxC, cmaxRC, cmaxLC

double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:ndir,2) :: vbarC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:ndir,2) :: vbarLC,vbarRC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   ndim)   :: cbarmin,cbarmax
integer                                   :: idimE, idimN

double precision :: dxinv(1:ndim),dxdim(1:ndim)
integer :: idims, iw, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3, hxOmin1,&
   hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3, ixCmin1,ixCmin2,ixCmin3,ixCmax1,&
   ixCmax2,ixCmax3, ixCRmin1,ixCRmin2,ixCRmin3,ixCRmax1,ixCRmax2,ixCRmax3,&
    kxCmin1,kxCmin2,kxCmin3,kxCmax1,kxCmax2,kxCmax3, kxRmin1,kxRmin2,kxRmin3,&
   kxRmax1,kxRmax2,kxRmax3
integer :: ixIsmin1,ixIsmin2,ixIsmin3,ixIsmax1,ixIsmax2,ixIsmax3
logical :: transport(1:nwflux)
double precision, parameter             :: epsrho = 1.0d-13, epstau = 1.0d-13
!-----------------------------------------------------------------------------
associate(wCT=>sCT%w%w, wnew=>snew%w%w, wold=>sold%w%w)



ixIsmin1=sCT%ws%ixGmin1;ixIsmin2=sCT%ws%ixGmin2;ixIsmin3=sCT%ws%ixGmin3;
ixIsmax1=sCT%ws%ixGmax1;ixIsmax2=sCT%ws%ixGmax2;ixIsmax3=sCT%ws%ixGmax3;

staggered : associate(wCTs=>sCT%ws%w, wnews=>snew%ws%w, wolds=>sold%ws%w)

  

if (typelimited .ne. 'predictor') call mpistop&
   ('tvdlfpos: Only implemented with typelimited=predictor')

! The flux calculation contracts by one in the idim direction it is applied.
! The limiter contracts the same directions by one more, so expand ixO by 2.
ixmin1=ixOmin1;ixmin2=ixOmin2;ixmin3=ixOmin3;ixmax1=ixOmax1;ixmax2=ixOmax2
ixmax3=ixOmax3;
do idims= idimmin,idimmax
   ixmin1=ixmin1-2*kr(idims,1);ixmin2=ixmin2-2*kr(idims,2)
   ixmin3=ixmin3-2*kr(idims,3);ixmax1=ixmax1+2*kr(idims,1)
   ixmax2=ixmax2+2*kr(idims,2);ixmax3=ixmax3+2*kr(idims,3);
end do


dxinv(1)=-qdt/dx1;dxinv(2)=-qdt/dx2;dxinv(3)=-qdt/dx3;
dxdim(1)=dx1;dxdim(2)=dx2;dxdim(3)=dx3;

!==================================================
! Low order flux for density:
!==================================================
do idims= idimmin,idimmax
   if (B0field) then
      select case (idims)
      case (1)
         myB0 => myB0_face1
      case (2)
         myB0 => myB0_face2
      case (3)
         myB0 => myB0_face3
      end select
   end if

   if (covariant) then
      select case (idims)
      case (1)
         myM => mygeo%mSurface1
      case (2)
         myM => mygeo%mSurface2
      case (3)
         myM => mygeo%mSurface3
      end select
   end if

   call assemble_indices()

   ! Get interface positions:
   xi(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:ndim) &
      = x(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:ndim)
   xi(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,idims) &
      = half* ( x(kxRmin1:kxRmax1,kxRmin2:kxRmax2,kxRmin3:kxRmax3,&
      idims)+x(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,idims) ) 

   
   ! ================================
   ! === Godunov : ==================
   ! ================================
   wRC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:nwflux+nwaux)&
      =wCT(kxRmin1:kxRmax1,kxRmin2:kxRmax2,kxRmin3:kxRmax3,1:nwflux+nwaux)
   wLC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:nwflux+nwaux)&
      =wCT(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:nwflux+nwaux)
   wRp(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:nwflux+nwaux)&
      =wprim(kxRmin1:kxRmax1,kxRmin2:kxRmax2,kxRmin3:kxRmax3,1:nwflux+nwaux)
   wLp(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:nwflux+nwaux)&
      =wprim(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:nwflux+nwaux)
   ! ================================

   !==================================================
   ! Get max Eigenvalue
   !==================================================
   call getcmax(wLC,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,cmaxLC,cmaxC,&
      .true.)
   cmaxLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
      = max(cmaxLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3),&
      dabs(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)))
   call getcmax(wRC,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,cmaxRC,cmaxC,&
      .true.)
   cmaxRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
      = max(cmaxRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3),&
      dabs(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)))
   ! now take the maximum of left and right states
   cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=max(cmaxRC&
      (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3),cmaxLC&
      (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3))
   !==================================================
   ! DONE: Get max Eigenvalue
   !==================================================


   ! Calculate velocities for transport fluxes
   call getv(wLp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,&
      ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,vLC)
   call getv(wRp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,&
      ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,vRC)

   call getflux(wLC,wLp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,fLC,transport)
   call getflux(wRC,wRp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,fRC,transport)
   do iw=1,nwflux
      if (transport(iw)) then
         fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)&
            =fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
            iw)+vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)*wLC&
            (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
         fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)&
            =fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
            iw)+vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)*wRC&
            (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
      end if
      ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
      fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)&
         =half*(fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         iw)+fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw))
      
      ! Add TVDLF dissipation to the flux
      ! To save memory we use fRC to store -cmax*half*(w_R-w_L)
      fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)&
         =-tvdlfeps*cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
         *half*(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw))
      ! fLC contains physical+dissipative fluxes
      fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)&
         =fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         iw)+fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
      
      if (slab) then
         fClow(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,idims)&
            =fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
      else
         select case (idims)
            case (1)
 !This where statement catches the axis where inverse metric becomes inf:
            where(mygeo%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3) .gt. 1.0d-9*mygeo%dvolume(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,ixCmin3:ixCmax3))
               fClow(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,1)&
                  =mygeo%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3)*fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3,iw)
            elsewhere
               fClow(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,1) &
                  = zero
            end where
            
            case (2)
 !This where statement catches the axis where inverse metric becomes inf:
            where(mygeo%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3) .gt. 1.0d-9*mygeo%dvolume(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,ixCmin3:ixCmax3))
               fClow(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,2)&
                  =mygeo%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3)*fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3,iw)
            elsewhere
               fClow(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,2) &
                  = zero
            end where
            
            case (3)
 !This where statement catches the axis where inverse metric becomes inf:
            where(mygeo%surfaceC3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3) .gt. 1.0d-9*mygeo%dvolume(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,ixCmin3:ixCmax3))
               fClow(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,3)&
                  =mygeo%surfaceC3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3)*fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3,iw)
            elsewhere
               fClow(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,3) &
                  = zero
            end where
            
         end select
      end if
   end do ! Next iw

end do ! Next idims
!==================================================
! DONE: low order flux
!==================================================


!==================================================
! High order flux:
!==================================================

do idims= idimmin,idimmax
   if (B0field) then
      select case (idims)
         case (1)
         myB0 => myB0_face1
         case (2)
         myB0 => myB0_face2
         case (3)
         myB0 => myB0_face3
      end select
   end if

   if (covariant) then
      select case (idims)
         case (1)
         myM => mygeo%mSurface1
         case (2)
         myM => mygeo%mSurface2
         case (3)
         myM => mygeo%mSurface3
      end select
   end if

   call assemble_indices()

   ! Get interface positions:
   xi(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:ndim) &
      = x(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:ndim)
   xi(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,idims) &
      = half* ( x(kxRmin1:kxRmax1,kxRmin2:kxRmax2,kxRmin3:kxRmax3,&
      idims)+x(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,idims) )


   ! ================================
   ! === Start with flat (Godunov): =
   ! ================================
   wRp(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:nwflux+nwaux)&
      =wprim(kxRmin1:kxRmax1,kxRmin2:kxRmax2,kxRmin3:kxRmax3,1:nwflux+nwaux)
   wLp(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:nwflux+nwaux)&
      =wprim(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:nwflux+nwaux)
   ! ================================

   call upwindLR(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCRmin1,&
      ixCRmin2,ixCRmin3,ixCRmax1,ixCRmax2,ixCRmax3,ixCRmin1,ixCRmin2,ixCRmin3,&
      ixCRmax1,ixCRmax2,ixCRmax3,idims,wprim,wprim,wLC,wRC,wLp,wRp,x,&
      dxdim(idims))
   
   wLC(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,b0_+idims)&
      =wCTs(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,idims)
   wRC(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,b0_+idims)&
      =wCTs(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,idims)
   wLp(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,b0_+idims)&
      =wCTs(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,idims)
   wRp(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,b0_+idims)&
      =wCTs(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,idims)
   

   !==================================================
   ! Get max Eigenvalue
   !==================================================
   call getcmax(wLC,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,cmaxLC,cmaxC,&
      .true.)
   cmaxLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
      = max(cmaxLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3),&
      dabs(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)))
   call getcmax(wRC,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,cmaxRC,cmaxC,&
      .true.)
   cmaxRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
      = max(cmaxRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3),&
      dabs(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)))
   ! now take the maximum of left and right states
   cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=max(cmaxRC&
      (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3),cmaxLC&
      (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3))
   !==================================================
   ! DONE: Get max Eigenvalue
   !==================================================

   
   ! ================================
   ! Calculate velocities for transport fluxes
   ! ================================
   call getv(wLp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,&
      ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,vLC)
   call getv(wRp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,&
      ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,vRC)
   ! ================================


   ! ================================
   ! Storage for various ways of EMF averaging:
   ! ================================
   select case (typeemf)
   case default
   case ('average')
   case ('uct1')
      ! Store magnitude of characteristics
      cbarmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idims)&
         =max(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3),zero)
      cbarmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idims)&
         =cbarmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idims)
      
      
 !If the current direction will be useful for a component of the electric field
      ! (1 on 2D or 3 in 3D), then store relevant velocities.
      if (idims.le.2*ndim-3) then
         idimN=mod(idims,ndir)+1 ! 'Next' direction
         idimE=mod(idims+1,ndir)+1 ! Electric field direction

         vbarLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idimE,1)&
            =vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
         vbarRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idimE,1)&
            =vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
         ! Calculate additional component of velocity needed
         call getv(wLp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idimN,&
            vbarLC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idimE,2))
         call getv(wRp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idimN,&
            vbarRC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idimE,2))
      end if

   case ('uct2','uct2+av')
      ! Store magnitude of characteristics
      cbarmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idims)&
         =max( cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3),zero)
      cbarmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idims)&
         =cbarmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idims)


      idimN=mod(idims,ndir)+1 ! 'Next' direction
      idimE=mod(idims+1,ndir)+1 ! Electric field direction

      call getv(wLp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idimN,&
         vbarLC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idims,1))
      call getv(wRp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idimN,&
         vbarRC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idims,1))

      vbarC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idims,1)&
         =half*(vbarLC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idims,&
         1) &
           +vbarRC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idims,1))

      call getv(wLp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idimE,&
         vbarLC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idims,2))
      call getv(wRp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idimE,&
         vbarRC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idims,2))

      vbarC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idims,2)&
         =half*(vbarLC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idims,&
         2) &
           +vbarRC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idims,2))

   end select
   ! ================================



   ! ================================
   ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
   ! ================================
   call getflux(wLC,wLp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,fLC,transport)
   call getflux(wRC,wRp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,fRC,transport)
   do iw=1,nwflux
      if (transport(iw)) then
         fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)&
            =fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
            iw)+vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)*wLC&
            (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
         fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)&
            =fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
            iw)+vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)*wRC&
            (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
      end if
      ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
      fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)&
         =half*(fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         iw)+fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw))

      ! Add TVDLF dissipation to the flux
      ! To save memory we use fRC to store -cmax*half*(w_R-w_L)
      fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)&
         =-tvdlfeps*cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
         *half*(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw))
      ! fLC contains physical+dissipative fluxes
      fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)&
         =fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         iw)+fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)

     if (slab) then
         fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,idims)&
            =fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
     else
         select case (idims)
            case (1)
 !This where statement catches the axis where inverse metric becomes inf:
            where(mygeo%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3) .gt. 1.0d-9*mygeo%dvolume(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,ixCmin3:ixCmax3))
               fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,1)&
                  =mygeo%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3)*fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3,iw)
            elsewhere
               fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,1) = zero
            end where
            
            case (2)
 !This where statement catches the axis where inverse metric becomes inf:
            where(mygeo%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3) .gt. 1.0d-9*mygeo%dvolume(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,ixCmin3:ixCmax3))
               fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,2)&
                  =mygeo%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3)*fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3,iw)
            elsewhere
               fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,2) = zero
            end where
            
            case (3)
 !This where statement catches the axis where inverse metric becomes inf:
            where(mygeo%surfaceC3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3) .gt. 1.0d-9*mygeo%dvolume(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,ixCmin3:ixCmax3))
               fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,3)&
                  =mygeo%surfaceC3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3)*fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3,iw)
            elsewhere
               fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,3) = zero
            end where
            
         end select
     end if

   end do ! Next iw
   ! ================================

   
end do ! Next idims

!==================================================
! DONE: High order flux
!==================================================



!==================================================
! Ensure positivity preserving flux in density:
!==================================================
do idims= idimmin,idimmax
   call assemble_indices()

   if (.not.slab) then
      call get_theta(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,&
         ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,epsrho,qdt,&
         mygeo%dvolume(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
         wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,rho_),&
         fClow(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,rho_,idims),&
         fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,rho_,idims),theta)
      call get_theta(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,&
         ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,epstau,qdt,&
         mygeo%dvolume(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
         wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,tau_),&
         fClow(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,tau_,idims),&
         fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,tau_,idims),&
         theta_tau)
   else
      call get_theta_slab(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,epsrho,&
         -dxinv(idims),wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         rho_),fClow(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,rho_,&
         idims),fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,rho_,&
         idims),theta)         
      call get_theta_slab(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,epstau,&
         -dxinv(idims),wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         tau_),fClow(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,tau_,&
         idims),fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,tau_,&
         idims),theta_tau)
   end if

   theta(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = &
      min(theta(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3),&
      theta_tau(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3))
   
   do iw = 1, nwflux
      fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,idims) &
         = theta(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)*fC&
         (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,&
         idims) + (one-theta(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3))*fClow(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3,iw,idims)
   end do

end do ! Next idims
!==================================================
! DONE: Positivity preserving fix
!==================================================



! Update values at faces
select case (typeemf)
case default
   call mpistop('please specify typeemf for staggered field upwinding')
case ('average')
   call updatefaces(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
      ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qdt,fC,fE,snew)
case ('uct1')
   call updatefacesuct1(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qdt,vbarRC,vbarLC,&
      cbarmin,cbarmax,fE,snew)
case ('uct2')
   call updatefacesuct2(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qdt,vbarC,cbarmin,&
      cbarmax,fE,snew)
case ('uct2+av')
   call updatefacesuct2av(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qdt,vbarC,cbarmin,&
      cbarmax,fC,fE,snew)
end select




!==================================================
! Update the state:
!==================================================
do idims= idimmin,idimmax
   hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
   hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
   hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);

   do iw=1,nwflux

      ! Multiply the fluxes by -dt/dx since Flux fixing expects this
      if (slab) then
         fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw,idims)&
            =dxinv(idims)*fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
            iw,idims)
         wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
            =wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            iw) + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw,&
            idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,iw,&
            idims))
      else
         select case (idims)
         case (1)
            fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw,1)&
               =-qdt*fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw,&
               idims)
            wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
               =wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) &
              + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw,&
                 1)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,iw,1))&
                 /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3)
         case (2)
            fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw,2)&
               =-qdt*fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw,&
               idims)
            wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
               =wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) &
              + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw,&
                 2)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,iw,2))&
                 /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3)
         case (3)
            fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw,3)&
               =-qdt*fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw,&
               idims)
            wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
               =wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) &
              + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw,&
                 3)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,iw,3))&
                 /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3)
         end select
      end if

   end do ! Next iw

end do ! Next idims


call faces2centers(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,snew)

!==================================================
! DONE: Update the state
!==================================================


if (covariant) myM => mygeo%m
if (.not.slab.and.idimmin==1) then 
   call addgeometry(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,wnew,x)
else
   if(nwaux>0) call getaux(.true.,wnew,x,ixImin1,ixImin2,ixImin3,ixImax1,&
      ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
      'tvdlfpos_wnew')
end if

call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim),ixImin1,ixImin2,&
   ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
   ixOmax3,1,nw,qtC,wCT,qt,wnew,x,.false.)



end associate staggered


end associate
!=============================================================================
contains 
!=============================================================================
  subroutine assemble_indices()

   hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
   hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
   hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);
   kxCmin1=ixImin1;kxCmin2=ixImin2;kxCmin3=ixImin3
   kxCmax1=ixImax1-kr(idims,1);kxCmax2=ixImax2-kr(idims,2)
   kxCmax3=ixImax3-kr(idims,3);
   kxRmin1=kxCmin1+kr(idims,1);kxRmin2=kxCmin2+kr(idims,2)
   kxRmin3=kxCmin3+kr(idims,3);kxRmax1=kxCmax1+kr(idims,1)
   kxRmax2=kxCmax2+kr(idims,2);kxRmax3=kxCmax3+kr(idims,3);

   select case (typeemf)
   case default
      ! ixC is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3; ixCmin1=hxOmin1
      ixCmin2=hxOmin2;ixCmin3=hxOmin3;
   case ('average')
      ! Flux-interpolated constrained transport needs one more layer in the tangential components:
      ixCmax1=ixOmax1+1-kr(idims,1);ixCmax2=ixOmax2+1-kr(idims,2)
      ixCmax3=ixOmax3+1-kr(idims,3); ixCmin1=hxOmin1-1+kr(idims,1)
      ixCmin2=hxOmin2-1+kr(idims,2);ixCmin3=hxOmin3-1+kr(idims,3);
   case ('uct1','uct2','uct2+av')
      ! Upwind constrained transport needs dixB more layers in the tangential components
      ixCmax1=ixOmax1+dixB-dixB*kr(idims,1)
      ixCmax2=ixOmax2+dixB-dixB*kr(idims,2)
      ixCmax3=ixOmax3+dixB-dixB*kr(idims,3)
      ixCmin1=hxOmin1-dixB+dixB*kr(idims,1)
      ixCmin2=hxOmin2-dixB+dixB*kr(idims,2)
      ixCmin3=hxOmin3-dixB+dixB*kr(idims,3);
   end select
   ixCRmin1=ixCmin1;ixCRmin2=ixCmin2;ixCRmin3=ixCmin3;ixCRmax1=ixCmax1
   ixCRmax2=ixCmax2;ixCRmax3=ixCmax3; !indices to reconstruct to

end subroutine assemble_indices
!=============================================================================
end subroutine tvdlfpos
!=============================================================================
subroutine get_theta_slab(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,eps,lambda,u,flow,&
   fhigh,theta)

use mod_amrvacdef

integer, intent(in)                             :: ixImin1,ixImin2,ixImin3,&
   ixImax1,ixImax2,ixImax3, ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,&
    idims
double precision, intent(in)                    :: eps, lambda
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
    intent(in)  :: u
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
    intent(in)  :: flow, fhigh
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
    intent(out) :: theta
! .. local ..
integer                                         :: ixCpmin1,ixCpmin2,ixCpmin3,&
   ixCpmax1,ixCpmax2,ixCpmax3, ixOpmin1,ixOpmin2,ixOpmin3,ixOpmax1,ixOpmax2,&
   ixOpmax3
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)              :: tmp, thp, thm
!-----------------------------------------------------------------------------
ixCpmin1=ixCmin1+kr(idims,1);ixCpmin2=ixCmin2+kr(idims,2)
ixCpmin3=ixCmin3+kr(idims,3);ixCpmax1=ixCmax1+kr(idims,1)
ixCpmax2=ixCmax2+kr(idims,2);ixCpmax3=ixCmax3+kr(idims,3);
ixOpmin1=ixCmin1;ixOpmin2=ixCmin2;ixOpmin3=ixCmin3; ixOpmax1=ixCpmax1
ixOpmax2=ixCpmax2;ixOpmax3=ixCpmax3;

thm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = one
thp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = one

tmp(ixOpmin1:ixOpmax1,ixOpmin2:ixOpmax2,ixOpmin3:ixOpmax3) = 0.5d0&
   /lambda*(u(ixOpmin1:ixOpmax1,ixOpmin2:ixOpmax2,ixOpmin3:ixOpmax3)-eps)

where (tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) .lt. &
   fhigh(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3))
   thm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = tmp(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3) - flow(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      ixCmin3:ixCmax3)
   thm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = thm(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3) / (-flow(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3) + fhigh(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3))
end where

where (tmp(ixCpmin1:ixCpmax1,ixCpmin2:ixCpmax2,ixCpmin3:ixCpmax3) .lt. &
   -fhigh(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3))
   thp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = tmp&
      (ixCpmin1:ixCpmax1,ixCpmin2:ixCpmax2,ixCpmin3:ixCpmax3) + &
      flow(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
   thp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = thp(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3) / (flow(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3) - fhigh(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3))
end where

theta(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = min(thm&
   (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3),thp(ixCmin1:ixCmax1,&
   ixCmin2:ixCmax2,ixCmin3:ixCmax3))
theta(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = min(max(theta&
   (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3),0.0d0),1.0d0)

end subroutine get_theta_slab
!=============================================================================
subroutine get_theta(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,&
   ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,eps,qdt,dV,u,flow,fhigh,&
   theta)

use mod_amrvacdef

integer, intent(in)                             :: ixImin1,ixImin2,ixImin3,&
   ixImax1,ixImax2,ixImax3, ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,&
    idims
double precision, intent(in)                    :: eps, qdt
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
    intent(in)  :: dV, u
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
    intent(in)  :: flow, fhigh
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
    intent(out) :: theta
! .. local ..
integer                                         :: ixCpmin1,ixCpmin2,ixCpmin3,&
   ixCpmax1,ixCpmax2,ixCpmax3, ixOpmin1,ixOpmin2,ixOpmin3,ixOpmax1,ixOpmax2,&
   ixOpmax3
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)              :: tmp, thp, thm
!-----------------------------------------------------------------------------
ixCpmin1=ixCmin1+kr(idims,1);ixCpmin2=ixCmin2+kr(idims,2)
ixCpmin3=ixCmin3+kr(idims,3);ixCpmax1=ixCmax1+kr(idims,1)
ixCpmax2=ixCmax2+kr(idims,2);ixCpmax3=ixCmax3+kr(idims,3);
ixOpmin1=ixCmin1;ixOpmin2=ixCmin2;ixOpmin3=ixCmin3; ixOpmax1=ixCpmax1
ixOpmax2=ixCpmax2;ixOpmax3=ixCpmax3;

thm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = one
thp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = one

tmp(ixOpmin1:ixOpmax1,ixOpmin2:ixOpmax2,ixOpmin3:ixOpmax3) &
   = 0.5d0*dV(ixOpmin1:ixOpmax1,ixOpmin2:ixOpmax2,ixOpmin3:ixOpmax3)&
   /qdt*(u(ixOpmin1:ixOpmax1,ixOpmin2:ixOpmax2,ixOpmin3:ixOpmax3)-eps)


where (tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) .lt. &
   fhigh(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3))
   thm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = tmp(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3) - flow(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      ixCmin3:ixCmax3)
   thm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = thm(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3) / (-flow(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3) + fhigh(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3))
end where

where (tmp(ixCpmin1:ixCpmax1,ixCpmin2:ixCpmax2,ixCpmin3:ixCpmax3) .lt. &
   -fhigh(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3))
   thp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = tmp&
      (ixCpmin1:ixCpmax1,ixCpmin2:ixCpmax2,ixCpmin3:ixCpmax3) + &
      flow(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
   thp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = thp(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3) / (flow(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3) - fhigh(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3))
end where

theta(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = min(thm&
   (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3),thp(ixCmin1:ixCmax1,&
   ixCmin2:ixCmax2,ixCmin3:ixCmax3))
theta(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = min(max(theta&
   (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3),0.0d0),1.0d0)

end subroutine get_theta
!=============================================================================
subroutine hll(method,qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idimmin,idimmax, qtC,sCT,&
   qt,snew,sold,wprim,fC,fE,dx1,dx2,dx3,x)

! method=='hll'  --> 2nd order HLL scheme.
! method=='hll1' --> 1st order HLL scheme.

use mod_amrvacdef

character(len=*), intent(in)                          :: method
double precision, intent(in)                          :: qdt, qtC, qt, dx1,&
   dx2,dx3
integer, intent(in)                                   :: ixImin1,ixImin2,&
   ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
   ixOmax3, idimmin,idimmax
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:ndim), intent(in) :: x
type(state)                                           :: sCT, snew, sold
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nw), intent(in)   :: wprim
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nwflux,1:ndim)    :: fC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:3)              :: fE

double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:ndim)             :: xi
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nw)     :: wLC, wRC, wLp, wRp
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nwflux) :: fLC, fRC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)          :: vLC, vRC

double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:ndir,2) :: vbarC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:ndir,2) :: vbarLC,vbarRC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   ndim)   :: cbarmin,cbarmax
integer                                 :: idimE,idimN

double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)      :: cmaxC, cmaxRC, cmaxLC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)      :: cminC, cminRC, cminLC
double precision, dimension(1:ndim)     :: dxinv, dxdim
integer, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)           &
       :: patchf
integer :: idims, iw, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3, hxOmin1,&
   hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3, ixCmin1,ixCmin2,ixCmin3,ixCmax1,&
   ixCmax2,ixCmax3, ixCRmin1,ixCRmin2,ixCRmin3,ixCRmax1,ixCRmax2,ixCRmax3,&
    jxCmin1,jxCmin2,jxCmin3,jxCmax1,jxCmax2,jxCmax3, kxCmin1,kxCmin2,kxCmin3,&
   kxCmax1,kxCmax2,kxCmax3, kxRmin1,kxRmin2,kxRmin3,kxRmax1,kxRmax2,kxRmax3
integer :: ixIsmin1,ixIsmin2,ixIsmin3,ixIsmax1,ixIsmax2,ixIsmax3
logical :: transport(1:nwflux), logiB
!-----------------------------------------------------------------------------
associate(wCT=>sCT%w%w, wnew=>snew%w%w, wold=>sold%w%w)


ixIsmin1=sCT%ws%ixGmin1;ixIsmin2=sCT%ws%ixGmin2;ixIsmin3=sCT%ws%ixGmin3;
ixIsmax1=sCT%ws%ixGmax1;ixIsmax2=sCT%ws%ixGmax2;ixIsmax3=sCT%ws%ixGmax3;

staggered : associate(wCTs=>sCT%ws%w, wnews=>snew%ws%w, wolds=>sold%ws%w)


  
logiB=(BnormLF.and.b0_>0)

if (idimmax>idimmin .and. typelimited=='original' .and. method&
   /='hll1')call mpistop("Error in hll: Unsplit dim. and original is limited")

if (covariant .and. typelimited .ne. 'predictor') call mpistop&
   ('hll: Covariant formulation only implemented with typelimited=predictor')

! The flux calculation contracts by one in the idim direction it is applied.
! The limiter contracts the same directions by one more, so expand ixO by 2.
ixmin1=ixOmin1;ixmin2=ixOmin2;ixmin3=ixOmin3;ixmax1=ixOmax1;ixmax2=ixOmax2
ixmax3=ixOmax3;
do idims= idimmin,idimmax
   ixmin1=ixmin1-2*kr(idims,1);ixmin2=ixmin2-2*kr(idims,2)
   ixmin3=ixmin3-2*kr(idims,3);ixmax1=ixmax1+2*kr(idims,1)
   ixmax2=ixmax2+2*kr(idims,2);ixmax3=ixmax3+2*kr(idims,3);
end do
if (ixImin1>ixmin1.or.ixImin2>ixmin2.or.ixImin3>ixmin3.or.ixImax1<ixmax1&
   .or.ixImax2<ixmax2.or.ixImax3<ixmax3) call mpistop("Error in hll : &
   Nonconforming input limits")

dxinv(1)=-qdt/dx1;dxinv(2)=-qdt/dx2;dxinv(3)=-qdt/dx3;
dxdim(1)=dx1;dxdim(2)=dx2;dxdim(3)=dx3;
do idims= idimmin,idimmax
   if (B0field) then
      select case (idims)
      case (1)
         myB0 => myB0_face1
      case (2)
         myB0 => myB0_face2
      case (3)
         myB0 => myB0_face3
      end select
   end if

   if (covariant) then
      select case (idims)
      case (1)
         myM => mygeo%mSurface1
      case (2)
         myM => mygeo%mSurface2
      case (3)
         myM => mygeo%mSurface3
      end select
   end if


   call assemble_indices()

   ! Get interface positions:
   xi(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:ndim) &
      = x(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:ndim)
   xi(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,idims) &
      = half* ( x(kxRmin1:kxRmax1,kxRmin2:kxRmax2,kxRmin3:kxRmax3,&
      idims)+x(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,idims) )


   ! ================================
   ! === Start with flat (Godunov): =
   ! ================================
   wRp(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:nwflux+nwaux)&
      =wprim(kxRmin1:kxRmax1,kxRmin2:kxRmax2,kxRmin3:kxRmax3,1:nwflux+nwaux)
   wLp(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:nwflux+nwaux)&
      =wprim(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:nwflux+nwaux)
   ! ================================

   
   ! ================================
   ! for hll (second order scheme): apply limiting
   ! ================================
   if (method=='hll') then
      select case (typelimited)
      case ('previous')
         call primitive(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,wold,x)
         call upwindLR(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixCRmin1,ixCRmin2,ixCRmin3,ixCRmax1,ixCRmax2,ixCRmax3,ixCRmin1,&
            ixCRmin2,ixCRmin3,ixCRmax1,ixCRmax2,ixCRmax3,idims,wold,wprim,wLC,&
            wRC,wLp,wRp,xi,dxdim(idims))
         call conserve(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,wold,x,patchfalse)
         
         wLC(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,b0_+idims)&
            =wolds(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,&
            idims)
         wRC(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,b0_+idims)&
            =wolds(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,&
            idims)
         wLp(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,b0_+idims)&
            =wolds(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,&
            idims)
         wRp(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,b0_+idims)&
            =wolds(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,&
            idims)
         
      case ('predictor')
         call upwindLR(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixCRmin1,ixCRmin2,ixCRmin3,ixCRmax1,ixCRmax2,ixCRmax3,ixCRmin1,&
            ixCRmin2,ixCRmin3,ixCRmax1,ixCRmax2,ixCRmax3,idims,wprim,wprim,&
            wLC,wRC,wLp,wRp,xi,dxdim(idims))
         
         wLC(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,b0_+idims)&
            =wCTs(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,&
            idims)
         wRC(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,b0_+idims)&
            =wCTs(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,idims)
         wLp(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,b0_+idims)&
            =wCTs(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,&
            idims)
         wRp(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,b0_+idims)&
            =wCTs(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,idims)
         
      case ('original')
         call primitive(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,wnew,x)
         call upwindLR(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixCRmin1,ixCRmin2,ixCRmin3,ixCRmax1,ixCRmax2,ixCRmax3,ixCRmin1,&
            ixCRmin2,ixCRmin3,ixCRmax1,ixCRmax2,ixCRmax3,idims,wnew,wprim,wLC,&
            wRC,wLp,wRp,xi,dxdim(idims))
         call conserve(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,wnew,x,patchfalse)
         
         wLC(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,b0_+idims)&
            =wnews(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,&
            idims)
         wRC(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,b0_+idims)&
            =wnews(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,&
            idims)
         wLp(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,b0_+idims)&
            =wnews(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,&
            idims)
         wRp(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,b0_+idims)&
            =wnews(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,&
            idims)
         
      case default
         call mpistop("Error in HLL: no such base for limiter")
      end select
   else
      wRC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:nwflux+nwaux) &
         = wCT(kxRmin1:kxRmax1,kxRmin2:kxRmax2,kxRmin3:kxRmax3,1:nwflux+nwaux)
      wLC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:nwflux+nwaux) &
         = wCT(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:nwflux+nwaux)
   end if
   ! ==== Done with limiting ========

   
   ! ================================
   ! For the high order hll scheme the limiter is based on
   ! the maximum eigenvalue, it is calculated in advance.
   ! ================================
   call getcmax(wLC,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,cmaxLC,cminLC,&
      .true.)
   call getcmax(wRC,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,cmaxRC,cminRC,&
      .true.)
   ! now take the maximum of left and right states
   cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=max(cmaxRC&
      (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3),cmaxLC&
      (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3))
   cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=min(cminRC&
      (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3),cminLC&
      (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3))
   ! ================================

   patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) =  1
   where(cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) >= zero)
        patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = -2
   elsewhere(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) <= zero)
        patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) =  2
   endwhere

   ! ================================
   ! Calculate velocities for transport fluxes
   ! ================================

   call getv(wLp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,&
      ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,vLC)
   call getv(wRp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,&
      ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,vRC)

   ! ================================
   ! Storage for various ways of EMF averaging:
   ! ================================
   select case (typeemf)
   case default
   case ('average')
   case ('uct1')
      ! Store magnitude of characteristics
      cbarmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idims)&
         =max(-cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3),zero)
      cbarmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idims)&
         =max( cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3),zero)

 !If the current direction will be useful for a component of the electric field
   ! (1 on 2D or 3 in 3D), then store relevant velocities.
      if (idims.le.2*ndim-3) then
         idimN=mod(idims,ndir)+1 ! 'Next' direction
         idimE=mod(idims+1,ndir)+1 ! Electric field direction

         ! Store velocities
         vbarLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idimE,1)&
            =vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
         vbarRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idimE,1)&
            =vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)

         ! and calculate additional components of velocity needed
         call getv(wLp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idimN,&
            vbarLC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idimE,2))
         call getv(wRp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idimN,&
            vbarRC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idimE,2))

      end if
   case ('uct2', 'uct2+av')
      ! Store magnitude of characteristics
      cbarmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idims)&
         =max(-cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3),zero)
      cbarmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idims)&
         =max( cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3),zero)

      idimN=mod(idims,ndir)+1 ! 'Next' direction
      idimE=mod(idims+1,ndir)+1 ! Electric field direction


      call getv(wLp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idimN,&
         vbarLC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idims,1))
      call getv(wRp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idimN,&
         vbarRC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idims,1))

      vbarC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idims,1)&
         =(cbarmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idims)*vbarLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idims,&
         1) &
           +cbarmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
              idims)*vbarRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
              idims,1))&
           /(cbarmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
              idims) + cbarmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
              ixCmin3:ixCmax3,idims))


      call getv(wLp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idimE,&
         vbarLC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idims,2))
      call getv(wRp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idimE,&
         vbarRC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idims,2))

      vbarC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idims,2)&
         =(cbarmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idims)*vbarLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idims,&
         2) &
           +cbarmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
              idims)*vbarRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
              idims,2))&
           /(cbarmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
              idims) + cbarmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
              ixCmin3:ixCmax3,idims))

   end select
   ! ================================



   ! ================================
   ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
   ! ================================
   
     if(any(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
        /= 2).or.(logiB)) then 
        call getflux(wLC,wLp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,fLC,&
           transport)
     end if
     if(any(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
        /=-2).or.(logiB)) then 
        call getflux(wRC,wRp,xi,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,fRC,&
           transport)
     end if

   do iw=1,nwflux
      if (transport(iw)) fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         iw)=fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         iw)+vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)*wLC&
         (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
      if (transport(iw)) fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         iw)=fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         iw)+vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)*wRC&
         (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)

     if (logiB.and.(iw==b0_+idims)) then
         ! flat B norm using tvdlf
         fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)&
            = half*((fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
            iw)+fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
            iw)) -tvdlfeps*max(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
            ixCmin3:ixCmax3),dabs(cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
            ixCmin3:ixCmax3)))*(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
            ixCmin3:ixCmax3,iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
            ixCmin3:ixCmax3,iw)))
     else
       where(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)==1)
         ! Add hll dissipation to the flux
         fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)&
            =(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)*fLC&
            (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
            iw)-cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)*fRC&
            (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
            iw) +tvdlfeps*cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
            ixCmin3:ixCmax3)*cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
            ixCmin3:ixCmax3)*(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
            ixCmin3:ixCmax3,iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
            ixCmin3:ixCmax3,iw)))/(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
            ixCmin3:ixCmax3)-cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
            ixCmin3:ixCmax3))
       elsewhere(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)== 2)
         fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)&
            =fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
       elsewhere(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)==-2)
         fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)&
            =fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
       endwhere
     endif

     if (slab) then
         fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,idims)&
            =fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
     else
         select case (idims)
            case (1)
 !This where statement catches the axis where inverse metric becomes inf:
            where(mygeo%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3) .gt. 1.0d-9*mygeo%dvolume(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,ixCmin3:ixCmax3))
               fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,1)&
                  =mygeo%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3)*fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3,iw)
            elsewhere
               fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,1) = zero
            end where
            
            case (2)
 !This where statement catches the axis where inverse metric becomes inf:
            where(mygeo%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3) .gt. 1.0d-9*mygeo%dvolume(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,ixCmin3:ixCmax3))
               fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,2)&
                  =mygeo%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3)*fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3,iw)
            elsewhere
               fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,2) = zero
            end where
            
            case (3)
 !This where statement catches the axis where inverse metric becomes inf:
            where(mygeo%surfaceC3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3) .gt. 1.0d-9*mygeo%dvolume(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,ixCmin3:ixCmax3))
               fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,3)&
                  =mygeo%surfaceC3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3)*fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3,iw)
            elsewhere
               fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,3) = zero
            end where
            
         end select
     end if

  end do ! Next iw
  ! ================================

  
end do ! Next idims


! Update values at faces
select case (typeemf)
case default
   call mpistop('please specify typeemf for staggered field upwinding')
case ('average')
   call updatefaces(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
      ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qdt,fC,fE,snew)
case ('uct1')
   call updatefacesuct1(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qdt,vbarRC,vbarLC,&
      cbarmin,cbarmax,fE,snew)
case ('uct2')
   call updatefacesuct2(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qdt,vbarC,cbarmin,&
      cbarmax,fE,snew)
case ('uct2+av')
   call updatefacesuct2av(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qdt,vbarC,cbarmin,&
      cbarmax,fC,fE,snew)
end select



! ================================
! Now update the state:
! ================================
do idims= idimmin,idimmax
   
   hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
   hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
   hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);
   
   do iw=1,nwflux

      ! Multiply the fluxes by -dt/dx since Flux fixing expects this
      if (slab) then
         fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw,idims)&
            =dxinv(idims)*fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
            iw,idims)
         wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
            =wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            iw) + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw,&
            idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,iw,&
            idims))
      else
         select case (idims)
         case (1)
            fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw,1)&
               =-qdt*fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw,&
               idims)
            wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
               =wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) &
              + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw,&
                 1)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,iw,1))&
                 /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3)
         case (2)
            fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw,2)&
               =-qdt*fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw,&
               idims)
            wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
               =wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) &
              + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw,&
                 2)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,iw,2))&
                 /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3)
         case (3)
            fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw,3)&
               =-qdt*fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw,&
               idims)
            wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
               =wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) &
              + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw,&
                 3)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,iw,3))&
                 /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3)
         end select
      end if

   end do ! Next iw
   
end do ! Next idims

! Centers must be filled to convert to conserved

call faces2centers(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,snew)

! ================================
! === Done updating state ========
! ================================


if (covariant) myM => mygeo%m
if (.not.slab.and.idimmin==1) then
   call addgeometry(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,wnew,x)
else
   if(nwaux>0) call getaux(.true.,wnew,x,ixImin1,ixImin2,ixImin3,ixImax1,&
      ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
      'hll_wnew')
end if

call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), ixImin1,ixImin2,&
   ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
   ixOmax3,1,nw,qtC,wCT,qt,wnew,x,.false.)


end associate staggered


end associate
!=============================================================================
contains 
!=============================================================================
  subroutine assemble_indices()

   hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
   hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
   hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);
   kxCmin1=ixImin1;kxCmin2=ixImin2;kxCmin3=ixImin3
   kxCmax1=ixImax1-kr(idims,1);kxCmax2=ixImax2-kr(idims,2)
   kxCmax3=ixImax3-kr(idims,3);
   kxRmin1=kxCmin1+kr(idims,1);kxRmin2=kxCmin2+kr(idims,2)
   kxRmin3=kxCmin3+kr(idims,3);kxRmax1=kxCmax1+kr(idims,1)
   kxRmax2=kxCmax2+kr(idims,2);kxRmax3=kxCmax3+kr(idims,3);

   select case (typeemf)
   case default
      ! ixC is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3; ixCmin1=hxOmin1
      ixCmin2=hxOmin2;ixCmin3=hxOmin3;
   case ('average')
      ! Flux-interpolated constrained transport needs one more layer in the tangential components:
      ixCmax1=ixOmax1+1-kr(idims,1);ixCmax2=ixOmax2+1-kr(idims,2)
      ixCmax3=ixOmax3+1-kr(idims,3); ixCmin1=hxOmin1-1+kr(idims,1)
      ixCmin2=hxOmin2-1+kr(idims,2);ixCmin3=hxOmin3-1+kr(idims,3);
   case ('uct1','uct2','uct2+av')
      ! Upwind constrained transport needs dixB more layers in the tangential components
      ixCmax1=ixOmax1+dixB-dixB*kr(idims,1)
      ixCmax2=ixOmax2+dixB-dixB*kr(idims,2)
      ixCmax3=ixOmax3+dixB-dixB*kr(idims,3)
      ixCmin1=hxOmin1-dixB+dixB*kr(idims,1)
      ixCmin2=hxOmin2-dixB+dixB*kr(idims,2)
      ixCmin3=hxOmin3-dixB+dixB*kr(idims,3);
   end select
   ixCRmin1=ixCmin1;ixCRmin2=ixCmin2;ixCRmin3=ixCmin3;ixCRmax1=ixCmax1
   ixCRmax2=ixCmax2;ixCRmax3=ixCmax3; !indices to reconstruct to

 end subroutine assemble_indices
!=============================================================================
end subroutine hll
!=============================================================================
