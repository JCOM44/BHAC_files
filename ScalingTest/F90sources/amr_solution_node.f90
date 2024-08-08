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
integer function getnode(ipe)
use mod_forest, only: igrid_inuse
use mod_amrvacdef

! getnode = get first available igrid on processor ipe

integer, intent(in) :: ipe

integer :: igrid, igrid_available
!----------------------------------------------------------------------------
igrid_available=0

do igrid=1,ngridshi
   if (igrid_inuse(igrid,ipe)) cycle

   igrid_available=igrid
   exit
end do

if (igrid_available==0) then
   write(unitterm,*) " out of nodal space - allowed ",ngridshi," grids"
   call mpistop("")
else
   getnode=igrid_available
   igrid_inuse(igrid,ipe)=.true.
end if

if (ipe==mype) then
   ! initialize nodal block
   node(1:nodehi,getnode) = 0
   rnode(1:rnodehi,getnode) = zero
end if

end function getnode
!=============================================================================
subroutine putnode(igrid,ipe)
use mod_forest
implicit none

! putnode = return igrid node on processor ipe
 
integer, intent(in) :: igrid, ipe
!----------------------------------------------------------------------------
igrid_inuse(igrid,ipe)=.false.

end subroutine putnode
!=============================================================================
subroutine alloc_node(igrid)
use mod_forest
use mod_amrvacdef

integer, intent(in) :: igrid

integer :: level, ig1,ig2,ig3, ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,&
   ixCoGmax2,ixCoGmax3, ixCoCoGmin1,ixCoCoGmin2,ixCoCoGmin3,ixCoCoGmax1,&
   ixCoCoGmax2,ixCoCoGmax3, ix
double precision :: rXmin1,rXmin2,rXmin3, dx1,dx2,dx3
!-----------------------------------------------------------------------------
ixCoGmin1=1;ixCoGmin2=1;ixCoGmin3=1;
ixCoGmax1=ixGhi1/2+dixB;ixCoGmax2=ixGhi2/2+dixB;ixCoGmax3=ixGhi3/2+dixB;

! set level information
level=igrid_to_node(igrid,mype)%node%level
ig1=igrid_to_node(igrid,mype)%node%ig1;ig2=igrid_to_node(igrid,mype)%node%ig2
ig3=igrid_to_node(igrid,mype)%node%ig3;

node(plevel_,igrid)=level
node(pig1_,igrid)=ig1
node(pig2_,igrid)=ig2
node(pig3_,igrid)=ig3

! set dx information
rnode(rpdx1_,igrid)=dx(1,level)
rnode(rpdx2_,igrid)=dx(2,level)
rnode(rpdx3_,igrid)=dx(3,level)

! determine the minimal and maximal corners
rnode(rpxmin1_,igrid)=xprobmin1+dble(ig1-1)*dg1(level)
rnode(rpxmin2_,igrid)=xprobmin2+dble(ig2-1)*dg2(level)
rnode(rpxmin3_,igrid)=xprobmin3+dble(ig3-1)*dg3(level)
rnode(rpxmax1_,igrid)=xprobmax1-dble(ng1(level)-ig1)*dg1(level)
rnode(rpxmax2_,igrid)=xprobmax2-dble(ng2(level)-ig2)*dg2(level)
rnode(rpxmax3_,igrid)=xprobmax3-dble(ng3(level)-ig3)*dg3(level)

! Allocate and fill the position arrays
call alloc_px(px(igrid),ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,1,ndim)
call alloc_px(pxCoarse(igrid),ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,&
   ixCoGmax2,ixCoGmax3,1,ndim)

dx1=rnode(rpdx1_,igrid)
dx2=rnode(rpdx2_,igrid)
dx3=rnode(rpdx3_,igrid)
rXmin1=rnode(rpxmin1_,igrid)-dixB*dx1
rXmin2=rnode(rpxmin2_,igrid)-dixB*dx2
rXmin3=rnode(rpxmin3_,igrid)-dixB*dx3
do ix=ixGlo1,ixGhi1
   px(igrid)%x(ix,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1)=rXmin1+(dble(ix)-half)*dx1
end do
do ix=ixGlo2,ixGhi2
   px(igrid)%x(ixGlo1:ixGhi1,ix,ixGlo3:ixGhi3,2)=rXmin2+(dble(ix)-half)*dx2
end do
do ix=ixGlo3,ixGhi3
   px(igrid)%x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ix,3)=rXmin3+(dble(ix)-half)*dx3
end do
dx1=2.0d0*rnode(rpdx1_,igrid)
dx2=2.0d0*rnode(rpdx2_,igrid)
dx3=2.0d0*rnode(rpdx3_,igrid)
rXmin1=rnode(rpxmin1_,igrid)-dixB*dx1
rXmin2=rnode(rpxmin2_,igrid)-dixB*dx2
rXmin3=rnode(rpxmin3_,igrid)-dixB*dx3
do ix=ixCoGmin1,ixCoGmax1
   pxCoarse(igrid)%x(ix,ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,1)&
      =rXmin1+(dble(ix)-half)*dx1
end do
do ix=ixCoGmin2,ixCoGmax2
   pxCoarse(igrid)%x(ixCoGmin1:ixCoGmax1,ix,ixCoGmin3:ixCoGmax3,2)&
      =rXmin2+(dble(ix)-half)*dx2
end do
do ix=ixCoGmin3,ixCoGmax3
   pxCoarse(igrid)%x(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,ix,3)&
      =rXmin3+(dble(ix)-half)*dx3
end do

! initialize solution space
call alloc_state(psold(igrid),ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3)
call alloc_state(ps(igrid),ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3)
call alloc_state(psCoarse(igrid),ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,&
   ixCoGmax2,ixCoGmax3)

if(residmin>smalldouble) then
   call alloc_state(psres(igrid),ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3)
end if

if (errorestimate==1) then
   ixCoCoGmin1=1;ixCoCoGmin2=1;ixCoCoGmin3=1;
   ixCoCoGmax1=ixCoGmax1/2+dixB;ixCoCoGmax2=ixCoGmax2/2+dixB
   ixCoCoGmax3=ixCoGmax3/2+dixB;
   call alloc_state(psCoCo(igrid),ixCoCoGmin1,ixCoCoGmin2,ixCoCoGmin3,&
      ixCoCoGmax1,ixCoCoGmax2,ixCoCoGmax3)
end if

if (.not.slab) call getgridgeo(igrid)

if (B0field) then
   call alloc_B0_grid(igrid)
   call set_B0_grid(igrid)
end if

end subroutine alloc_node
!=============================================================================
subroutine dealloc_node(igrid)

use mod_amrvacdef

integer, intent(in) :: igrid
!-----------------------------------------------------------------------------
if (igrid==0) then
   call mpistop("trying to delete a non-existing grid in dealloc_node")
end if

! We don't deallocate the solution arrays any more 
! but do some smart checking when trying to re-allocate
! This helps to avoid memory fragmentation.  

! We should extend this strategy to the geometry datastructures.

if (.not.slab) call putgridgeo(igrid)

if (B0field) call dealloc_B0_grid(igrid)

end subroutine dealloc_node
!=============================================================================
