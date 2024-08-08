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
subroutine refine_grid(child_igrid,child_ipe,igrid,ipe,active)

use mod_amrvacdef

integer, dimension(2,2,2), intent(in) :: child_igrid, child_ipe
integer, intent(in) :: igrid, ipe
logical, intent(in) :: active

integer :: ic1,ic2,ic3
!-----------------------------------------------------------------------------

! allocate solution space for new children
do ic3=1,2
do ic2=1,2
do ic1=1,2
   call alloc_node(child_igrid(ic1,ic2,ic3))
end do
end do
end do

if ((time_advance .and. active).or.convert.or.firstprocess) then
   ! prolong igrid to new children
   call prolong_grid(child_igrid,child_ipe,igrid,ipe)
else
   ! Fill new created children with initial condition
   do ic3=1,2
   do ic2=1,2
   do ic1=1,2
      call set_tmpGlobals(child_igrid(ic1,ic2,ic3))
      call initial_condition(child_igrid(ic1,ic2,ic3))
   end do
   end do
   end do
end if

! remove solution space of igrid
call dealloc_node(igrid)

end subroutine refine_grid
!=============================================================================
subroutine prolong_grid(child_igrid,child_ipe,igrid,ipe)

use mod_amr_fct

use mod_amrvacdef

integer, dimension(2,2,2), intent(in) :: child_igrid, child_ipe
integer, intent(in) :: igrid, ipe

integer :: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3, ichild, ixComin1,&
   ixComin2,ixComin3,ixComax1,ixComax2,ixComax3, ic1,ic2,ic3
double precision :: dxCo1,dxCo2,dxCo3, xComin1,xComin2,xComin3, dxFi1,dxFi2,&
   dxFi3, xFimin1,xFimin2,xFimin3
!-----------------------------------------------------------------------------
if (covariant) myM=>pgeo(igrid)%m

if (typegridfill=="linear") then
   dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
   dxlevel(3)=rnode(rpdx3_,igrid);
   if (amrentropy) then
      ixmin1=ixMlo1-1;ixmin2=ixMlo2-1;ixmin3=ixMlo3-1;ixmax1=ixMhi1+1
      ixmax2=ixMhi2+1;ixmax3=ixMhi3+1;
      call e_to_rhos(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixmin1,ixmin2,&
         ixmin3,ixmax1,ixmax2,ixmax3,pw(igrid)%w,px(igrid)%x)
   else if (prolongprimitive) then
      ixmin1=ixMlo1-1;ixmin2=ixMlo2-1;ixmin3=ixMlo3-1;ixmax1=ixMhi1+1
      ixmax2=ixMhi2+1;ixmax3=ixMhi3+1;
      call primitive(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixmin1,ixmin2,&
         ixmin3,ixmax1,ixmax2,ixmax3,pw(igrid)%w,px(igrid)%x)
   end if

   xComin1=rnode(rpxmin1_,igrid)
   xComin2=rnode(rpxmin2_,igrid)
   xComin3=rnode(rpxmin3_,igrid)
   dxCo1=rnode(rpdx1_,igrid)
   dxCo2=rnode(rpdx2_,igrid)
   dxCo3=rnode(rpdx3_,igrid)
end if

call old_neighbors(child_igrid,child_ipe,igrid,ipe)

do ic3=1,2
do ic2=1,2
do ic1=1,2
   ichild=child_igrid(ic1,ic2,ic3)


   if (covariant) myM=>pgeo(ichild)%m
   
   ixComin1=ixMlo1+(ic1-1)*(ixMhi1-ixMlo1+1)/2
   ixComin2=ixMlo2+(ic2-1)*(ixMhi2-ixMlo2+1)/2
   ixComin3=ixMlo3+(ic3-1)*(ixMhi3-ixMlo3+1)/2
   ixComax1=ixMhi1+(ic1-2)*(ixMhi1-ixMlo1+1)/2
   ixComax2=ixMhi2+(ic2-2)*(ixMhi2-ixMlo2+1)/2
   ixComax3=ixMhi3+(ic3-2)*(ixMhi3-ixMlo3+1)/2

   if (typegridfill=="linear") then
      xFimin1=rnode(rpxmin1_,ichild)
      xFimin2=rnode(rpxmin2_,ichild)
      xFimin3=rnode(rpxmin3_,ichild)
      dxFi1=rnode(rpdx1_,ichild)
      dxFi2=rnode(rpdx2_,ichild)
      dxFi3=rnode(rpdx3_,ichild)

      call prolong_2nd(ps(igrid),px(igrid)%x,ixComin1,ixComin2,ixComin3,&
         ixComax1,ixComax2,ixComax3,ps(ichild),px(ichild)%x, dxCo1,dxCo2,&
         dxCo3,xComin1,xComin2,xComin3,dxFi1,dxFi2,dxFi3,xFimin1,xFimin2,&
         xFimin3,ichild)
   else
      call prolong_1st(ps(igrid),ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,&
         ixComax3,ps(ichild),px(ichild)%x)
   end if
end do
end do
end do
   
if (covariant) myM=>pgeo(igrid)%m
if (typegridfill=="linear") then
   if (amrentropy) then
      call rhos_to_e(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixmin1,ixmin2,&
         ixmin3,ixmax1,ixmax2,ixmax3,pw(igrid)%w,px(igrid)%x)
   else if (prolongprimitive) then
      call conserve(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixmin1,ixmin2,&
         ixmin3,ixmax1,ixmax2,ixmax3,pw(igrid)%w,px(igrid)%x,patchfalse)
   end if
end if

end subroutine prolong_grid
!=============================================================================
subroutine prolong_2nd(sCo,xCo,ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,&
   ixComax3,sFi,xFi,dxCo1,dxCo2,dxCo3,xComin1,xComin2,xComin3,dxFi1,dxFi2,&
   dxFi3,xFimin1,xFimin2,xFimin3,igridFi)
  
  use mod_amr_fct
  
use mod_amrvacdef

integer, intent(in) :: ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,ixComax3,&
    igridFi
double precision, intent(in) :: dxCo1,dxCo2,dxCo3, xComin1,xComin2,xComin3,&
    dxFi1,dxFi2,dxFi3, xFimin1,xFimin2,xFimin3
double precision, intent(in) :: xCo(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
   1:ndim), xFi(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:ndim)
type(state), intent(in)      :: sCo
type(state), intent(inout)   :: sFi

integer :: ixCo1,ixCo2,ixCo3, jxCo1,jxCo2,jxCo3, hxCo1,hxCo2,hxCo3, ixFi1,&
   ixFi2,ixFi3, ix1,ix2,ix3, idim, iw
integer :: ixFimin1,ixFimin2,ixFimin3,ixFimax1,ixFimax2,ixFimax3
double precision :: slopeL, slopeR, slopeC, signC, signR
double precision :: slope(nw,ndim)
double precision :: xCo1,xCo2,xCo3, xFi1,xFi2,xFi3, eta1,eta2,eta3, invdxCo1,&
   invdxCo2,invdxCo3

logical :: fine_min1,fine_min2,fine_min3,fine_max1,fine_max2,fine_max3

!-----------------------------------------------------------------------------
associate(wCo=>sCo%w%w, wFi=>sFi%w%w)
invdxCo1=1.d0/dxCo1;invdxCo2=1.d0/dxCo2;invdxCo3=1.d0/dxCo3;


do ixCo3 = ixComin3,ixComax3
   ! cell-centered coordinates of coarse grid point
   xCo3=xComin3+(dble(ixCo3-dixB)-half)*dxCo3

   ixFi3=2*(ixCo3-ixComin3)+ixMlo3
do ixCo2 = ixComin2,ixComax2
   ! cell-centered coordinates of coarse grid point
   xCo2=xComin2+(dble(ixCo2-dixB)-half)*dxCo2

   ixFi2=2*(ixCo2-ixComin2)+ixMlo2
do ixCo1 = ixComin1,ixComax1
   ! cell-centered coordinates of coarse grid point
   xCo1=xComin1+(dble(ixCo1-dixB)-half)*dxCo1

   ixFi1=2*(ixCo1-ixComin1)+ixMlo1

   do idim=1,ndim
      hxCo1=ixCo1-kr(1,idim)
      hxCo2=ixCo2-kr(2,idim)
      hxCo3=ixCo3-kr(3,idim)
      jxCo1=ixCo1+kr(1,idim)
      jxCo2=ixCo2+kr(2,idim)
      jxCo3=ixCo3+kr(3,idim)

      do iw=1,nw
         slopeL=wCo(ixCo1,ixCo2,ixCo3,iw)-wCo(hxCo1,hxCo2,hxCo3,iw)
         slopeR=wCo(jxCo1,jxCo2,jxCo3,iw)-wCo(ixCo1,ixCo2,ixCo3,iw)
         slopeC=half*(slopeR+slopeL)

         ! get limited slope
         signR=sign(one,slopeR)
         signC=sign(one,slopeC)
         select case(typeprolonglimit)
         case('minmod')
           slope(iw,idim)=signR*max(zero,min(dabs(slopeR), signR*slopeL))
         case('woodward')
           slope(iw,idim)=two*signR*max(zero,min(dabs(slopeR), signR*slopeL,&
              signR*half*slopeC))
         case('mcbeta')
           slope(iw,idim)=signR*max(zero,min(mcbeta*dabs(slopeR),&
               mcbeta*signR*slopeL,signR*slopeC))
         case('koren')
           slope(iw,idim)=signR*max(zero,min(two*signR*slopeL,&
               (dabs(slopeR)+two*slopeL*signR)*third,two*dabs(slopeR)))
         case default
           slope(iw,idim)=signC*max(zero,min(dabs(slopeC), signC*slopeL,&
              signC*slopeR))
         end select
      end do
   end do
   do ix3=ixFi3,ixFi3+1
      ! cell-centered coordinates of fine grid point
      xFi3=xFimin3+(dble(ix3-dixB)-half)*dxFi3
   do ix2=ixFi2,ixFi2+1
      ! cell-centered coordinates of fine grid point
      xFi2=xFimin2+(dble(ix2-dixB)-half)*dxFi2
   do ix1=ixFi1,ixFi1+1
      ! cell-centered coordinates of fine grid point
      xFi1=xFimin1+(dble(ix1-dixB)-half)*dxFi1

      ! normalized distance between fine/coarse cell center
      ! in coarse cell: ranges from -0.5 to 0.5 in each direction
      ! (origin is coarse cell center)
      if (slab) then
         eta1=(xFi1-xCo1)*invdxCo1;eta2=(xFi2-xCo2)*invdxCo2
         eta3=(xFi3-xCo3)*invdxCo3;
      else
         eta1=(xFi1-xCo1)*invdxCo1 *two*(one-pgeo(igridFi)%dvolume(ix1,ix2,&
            ix3) /sum(pgeo(igridFi)%dvolume(ixFi1:ixFi1+1,ix2,ix3))) 
         eta2=(xFi2-xCo2)*invdxCo2 *two*(one-pgeo(igridFi)%dvolume(ix1,ix2,&
            ix3) /sum(pgeo(igridFi)%dvolume(ix1,ixFi2:ixFi2+1,ix3))) 
         eta3=(xFi3-xCo3)*invdxCo3 *two*(one-pgeo(igridFi)%dvolume(ix1,ix2,&
            ix3) /sum(pgeo(igridFi)%dvolume(ix1,ix2,ixFi3:ixFi3+1))) 
      end if

      wFi(ix1,ix2,ix3,1:nw) = wCo(ixCo1,ixCo2,ixCo3,1:nw) + (slope(1:nw,&
         1)*eta1)+(slope(1:nw,2)*eta2)+(slope(1:nw,3)*eta3)
   end do
   end do
   end do
end do
end do
end do

call already_fine(sFi,igridFi,fine_min1,fine_min2,fine_min3,fine_max1,&
   fine_max2,fine_max3)
call prolong_2nd_stg(sCo,sFi,ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,&
   ixComax3,ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,dxCo1,dxCo2,dxCo3,&
   xComin1,xComin2,xComin3,dxFi1,dxFi2,dxFi3,xFimin1,xFimin2,xFimin3,.false.,&
   fine_min1,fine_min2,fine_min3,fine_max1,fine_max2,fine_max3)

if (amrentropy) then
   call rhos_to_e(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixMlo1,ixMlo2,&
      ixMlo3,ixMhi1,ixMhi2,ixMhi3,wFi,xFi)
else if (prolongprimitive) then
   call conserve(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixMlo1,ixMlo2,&
      ixMlo3,ixMhi1,ixMhi2,ixMhi3,wFi,xFi,patchfalse)
end if


end associate
end subroutine prolong_2nd
!=============================================================================
subroutine prolong_1st(sCo,ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,&
   ixComax3,sFi,xFi)

use mod_amrvacdef

integer, intent(in)          :: ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,&
   ixComax3
double precision, intent(in) :: xFi(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
   1:ndim)
type(state), intent(in)      :: sCo
type(state), intent(out)     :: sFi
! .. local ..
integer                      :: ixCo1,ixCo2,ixCo3, ixFi1,ixFi2,ixFi3, iw
integer                      :: ixFimin1,ixFimin2,ixFimin3,ixFimax1,ixFimax2,&
   ixFimax3
!-----------------------------------------------------------------------------
associate(wCo=>sCo%w%w, wFi=>sFi%w%w)
  
do ixCo3 = ixComin3,ixComax3
   ixFi3=2*(ixCo3-ixComin3)+ixMlo3
do ixCo2 = ixComin2,ixComax2
   ixFi2=2*(ixCo2-ixComin2)+ixMlo2
do ixCo1 = ixComin1,ixComax1
   ixFi1=2*(ixCo1-ixComin1)+ixMlo1
   forall(iw=1:nw) wFi(ixFi1:ixFi1+1,ixFi2:ixFi2+1,ixFi3:ixFi3+1,iw)&
      =wCo(ixCo1,ixCo2,ixCo3,iw)
end do
end do
end do

end associate
end subroutine prolong_1st
!=============================================================================
