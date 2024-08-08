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
subroutine initialize_vars
  use mod_forest
  use mod_metric, only: init_metric
use mod_amrvacdef

integer :: igrid, level, ipe, ig1,ig2,ig3, i, j, k
logical :: ok
!-----------------------------------------------------------------------------

! Initialize Kronecker delta, and Levi-Civita tensor
do i=0,3
   do j=0,3
      if(i==j)then
         kr(i,j)=1
      else
         kr(i,j)=0
      endif
      if (i .gt. 0 .and. j .gt. 0) then
         do k=1,3
            if(i==j.or.j==k.or.k==i)then
               lvc(i,j,k)=0
            else if(i+1==j.or.i-2==j)then
               lvc(i,j,k)=1
            else
               lvc(i,j,k)=-1
            endif
         enddo
      endif
   enddo
enddo

! set time, time counter
if(.not.treset)t=zero
if(.not.itreset)it=0
dt=zero
dtimpl=zero
itmin=0

if(.not.time_accurate.or.residmin>smalldouble) then
  residual=one
endif 

! set all dt to zero
dt_grid(1:ngridshi)=zero

! check resolution
if (mod(ixGhi1,2)/=0.or.mod(ixGhi2,2)/=0.or.mod(ixGhi3,2)/=0) then
   call mpistop("mesh widths must give even number grid points")
end if
ixMlo1=ixGlo1+dixB;ixMlo2=ixGlo2+dixB;ixMlo3=ixGlo3+dixB;ixMhi1=ixGhi1-dixB
ixMhi2=ixGhi2-dixB;ixMhi3=ixGhi3-dixB;
if (errorestimate==1) then
   if (mod(ixMhi1-ixMlo1+1,4)/=0.or.mod(ixMhi2-ixMlo2+1,4)/=0&
      .or.mod(ixMhi3-ixMlo3+1,4)/=0) then
      call mpistop("mesh widths must be divisible by 4 for Richardson")
   end if
end if

if (nbufferx1>(ixMhi1-ixMlo1+1).or.nbufferx2>(ixMhi2-ixMlo2+1)&
   .or.nbufferx3>(ixMhi3-ixMlo3+1)) then
   write(unitterm,*) 'nbufferx^D bigger than mesh size makes no sense.'
   write(unitterm,*) 'Decrease nbufferx or increase mesh size'
   call mpistop("")
end if

! initialize dx arrays on finer (>1) levels
do level=2,mxnest
   dx(1,level) = dx(1,level-1) * half
   dx(2,level) = dx(2,level-1) * half
   dx(3,level) = dx(3,level-1) * half  ! refine ratio 2
end do

! domain decomposition
! physical extent of a grid block at level 1, per dimension
dg1(1)=dx(1,1)*dble(ixGhi1-2*dixB)
dg2(1)=dx(2,1)*dble(ixGhi2-2*dixB)
dg3(1)=dx(3,1)*dble(ixGhi3-2*dixB)

! number of grid blocks at level 1 in simulation domain, per dimension
ng1(1)=nint((xprobmax1-xprobmin1)/dg1(1))
ng2(1)=nint((xprobmax2-xprobmin2)/dg2(1))
ng3(1)=nint((xprobmax3-xprobmin3)/dg3(1))

! total number of grid blocks at level 1
nglev1=ng1(1)*ng2(1)*ng3(1)

do level=2,mxnest
   dg1(level)=half*dg1(level-1);dg2(level)=half*dg2(level-1)
   dg3(level)=half*dg3(level-1);
   ng1(level)=ng1(level-1)*2;ng2(level)=ng2(level-1)*2
   ng3(level)=ng3(level-1)*2;
end do

! check that specified stepsize correctly divides domain
ok=((abs(dble(ng1(1))*dg1(1)-(xprobmax1-xprobmin1))<=smalldouble)&
   .and.(abs(dble(ng2(1))*dg2(1)-(xprobmax2-xprobmin2))<=smalldouble)&
   .and.(abs(dble(ng3(1))*dg3(1)-(xprobmax3-xprobmin3))<=smalldouble))
if (.not.ok) then
   write(unitterm,*)"domain cannot be divided by meshes of given gridsize"
   call mpistop("domain cannot be divided by meshes of given gridsize")
end if


poleB=.false.
if (.not.slab) call set_pole

do igrid=1,ngridshi
! All nullification already at the declaration stage.   
! Associate the state structure with the pw arrays.  

   psold(igrid) = state(igrid=igrid,w=pwold(igrid),x=px(igrid),geo&
      =pgeo(igrid)  , ws=pwsold(igrid) )

   ps(igrid)        = state(igrid=igrid,w=pw(igrid),x=px(igrid),geo&
      =pgeo(igrid)  , ws=pws(igrid) )

   ps1(igrid)       = state(igrid=igrid,w=pw1(igrid),x=px(igrid),geo&
      =pgeo(igrid)  , ws=pws1(igrid) )

   psCoarse(igrid)  = state(igrid=igrid,w=pwCoarse(igrid),x&
      =pxCoarse(igrid),geo=pgeoCoarse(igrid), is_coarse=.true.  , ws&
      =pwsCoarse(igrid) )
   psCoCo(igrid)    = state(igrid=igrid,w=pwCoCo(igrid),geo&
      =pgeoCoCo(igrid), is_coarse=.true.  , ws=pwsCoCo(igrid) )

   if (nstep>2) then
      ps2(igrid)    = state(igrid=igrid,w=pw2(igrid),x=px(igrid),geo&
         =pgeo(igrid)  , ws=pws2(igrid) )
   end if

   if (nstep>3) then
      ps3(igrid)    = state(igrid=igrid,w=pw3(igrid),x=px(igrid),geo&
         =pgeo(igrid)  , ws=pws3(igrid) )
   end if

   if (nstep>4) then
      ps4(igrid)    = state(igrid=igrid,w=pw4(igrid),x=px(igrid),geo&
         =pgeo(igrid)  , ws=pws4(igrid) )
   end if

   if (residmin>smalldouble) then
      psres(igrid)  = state(igrid=igrid,w=pwres(igrid),x=px(igrid),geo&
         =pgeo(igrid)  , ws=pwsres(igrid) )
   end if

end do

! on each processor, create for later use a default patch array
allocate(patchfalse(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3))
patchfalse(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)=.false.

! initialize connectivity data
igridstail=0

! allocate memory for forest data structures
allocate(level_head(mxnest),level_tail(mxnest))
do level=1,mxnest
   nullify(level_head(level)%node,level_tail(level)%node)
end do

allocate(igrid_to_node(ngridshi,0:npe-1))
do ipe=0,npe-1
   do igrid=1,ngridshi
      nullify(igrid_to_node(igrid,ipe)%node)
   end do
end do

allocate(sfc(1:3,ngridshi*npe))

allocate(igrid_to_sfc(ngridshi))

sfc=0
allocate(Morton_start(0:npe-1),Morton_stop(0:npe-1))
allocate(Morton_sub_start(0:npe-1),Morton_sub_stop(0:npe-1))

allocate(nleafs_level(1:nlevelshi))

allocate(coarsen(ngridshi,0:npe-1),refine(ngridshi,0:npe-1))
coarsen=.false.
refine=.false.
if (nbufferx1/=0.or.nbufferx2/=0.or.nbufferx3/=0) then
   allocate(buffer(ngridshi,0:npe-1))
   buffer=.false.
end if
allocate(igrid_inuse(ngridshi,0:npe-1))
igrid_inuse=.false.

allocate(tree_root(1:ng1(1),1:ng2(1),1:ng3(1)))
do ig3=1,ng3(1)
do ig2=1,ng2(1)
do ig1=1,ng1(1)
   nullify(tree_root(ig1,ig2,ig3)%node)
end do
end do
end do


! default the physical scaling parameters:
UNIT_LENGTH   = ONE
UNIT_DENSITY  = ONE
UNIT_VELOCITY = ONE

! initialize the meta-info for the metric datastructure:
if (covariant) then
   call init_metric
end if


end subroutine initialize_vars
!=============================================================================
subroutine set_tmpGlobals(igrid)

  use mod_amrvacdef

  integer, intent(in)              :: igrid
  !-----------------------------------------------------------------------------

  saveigrid       = igrid
  typelimiter     = type_limiter(node(plevel_,igrid))
  typegradlimiter = type_gradient_limiter(node(plevel_,igrid))

  dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
  dxlevel(3)=rnode(rpdx3_,igrid);
  
  if (.not.slab) mygeo => pgeo(igrid)
  if(covariant)myM => mygeo%m
  if (B0field) then
     myB0_cell => pB0_cell(igrid)
     myB0      => pB0_cell(igrid)
     myB0_face1 => pB0_face1(igrid)
     myB0_face2 => pB0_face2(igrid)
     myB0_face3 => pB0_face3(igrid)
  end if

end subroutine set_tmpGlobals
!=============================================================================
