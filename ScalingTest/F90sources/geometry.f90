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
subroutine set_pole

use mod_amrvacdef
!-----------------------------------------------------------------------------
select case (typeaxial)
case ("spherical") 
  ! For spherical grid, check whether phi-direction is periodic
  if(periodB(ndim)) then
    if(2/=3) call mpistop("set setamrvac -phi=3 and recompile")
    if(mod(ng3(1),2)/=0) &
      call mpistop("Number of meshes in phi-direction should be even!")
    if(abs(xprobmin2)<smalldouble) then
      if(mype==0) write(unitterm,*) &
        "Will apply pi-periodic conditions at northpole!"
      poleB(1,2)=.true.
    else
      if(mype==0) write(unitterm,*) "There is no northpole!"
    end if
    if(abs(xprobmax2-dpi)<smalldouble) then
      if(mype==0) write(unitterm,*) &
        "Will apply pi-periodic conditions at southpole!"
      poleB(2,2)=.true.
    else
      if(mype==0) write(unitterm,*) "There is no southpole!"
    end if
  end if
case ("cylindrical")
  if (1==2.and.periodB(1)) then
    if(mod(ng1(1),2)/=0) &
      call mpistop("Number of meshes in phi-direction should be even!")
    if(abs(xprobmin1)<smalldouble) then
      if(mype==0) write(unitterm,*) "Will apply pi-periodic conditions at r=0"
      poleB(1,1)=.true.
    else
      if(mype==0) write(unitterm,*) "There is no cylindrical axis!"
    end if
  end if
  if (2==2.and.periodB(2)) then
    if(mod(ng2(1),2)/=0) &
      call mpistop("Number of meshes in phi-direction should be even!")
    if(abs(xprobmin1)<smalldouble) then
      if(mype==0) write(unitterm,*) "Will apply pi-periodic conditions at r=0"
      poleB(1,1)=.true.
    else
      if(mype==0) write(unitterm,*) "There is no cylindrical axis!"
    end if
  end if
  if (3==2.and.periodB(3)) then
    if(mod(ng3(1),2)/=0) &
      call mpistop("Number of meshes in phi-direction should be even!")
    if(abs(xprobmin1)<smalldouble) then
      if(mype==0) write(unitterm,*) "Will apply pi-periodic conditions at r=0"
      poleB(1,1)=.true.
    else
      if(mype==0) write(unitterm,*) "There is no cylindrical axis!"
    end if
  end if
end select

end subroutine set_pole
!=============================================================================
subroutine getgridgeo(igrid)

use mod_metric
use mod_amrvacdef

integer, intent(in) :: igrid

integer :: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3, ixCoGmin1,ixCoGmin2,&
   ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3, ixCoMmin1,ixCoMmin2,ixCoMmin3,&
   ixCoMmax1,ixCoMmax2,ixCoMmax3, ixComin1,ixComin2,ixComin3,ixComax1,&
   ixComax2,ixComax3, ixCoCoGmin1,ixCoCoGmin2,ixCoCoGmin3,ixCoCoGmax1,&
   ixCoCoGmax2,ixCoCoGmax3, ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,&
   ixGextmax2,ixGextmax3
double precision :: xmin1,xmin2,xmin3, dx1,dx2,dx3
!-----------------------------------------------------------------------------
!!!!!!!!!!!!!!!! For staggered constrained transport, surfaces must be defined
!!!!!!!!!!!!!!!! also on the faces with index 0. This at the moment is done only
!!!!!!!!!!!!!!!! in fillgeo_covariant
!ix^L=ixM^LL^LADD1;
ixmin1=ixGlo1+1;ixmin2=ixGlo2+1;ixmin3=ixGlo3+1;ixmax1=ixGhi1-1
ixmax2=ixGhi2-1;ixmax3=ixGhi3-1;
if (2*int(dixB/2)==dixB) then
   ixGextmin1=ixGlo1;ixGextmin2=ixGlo2;ixGextmin3=ixGlo3;ixGextmax1=ixGhi1
   ixGextmax2=ixGhi2;ixGextmax3=ixGhi3;
else
   ixGextmin1=ixGlo1-1;ixGextmin2=ixGlo2-1;ixGextmin3=ixGlo3-1
   ixGextmax1=ixGhi1+1;ixGextmax2=ixGhi2+1;ixGextmax3=ixGhi3+1;
end if

allocate(pgeo(igrid)%surfaceC1(ixGlo1-1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
   pgeo(igrid)%surfaceC2(ixGlo1:ixGhi1,ixGlo2-1:ixGhi2,ixGlo3:ixGhi3),&
   pgeo(igrid)%surfaceC3(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3-1:ixGhi3),&
    pgeo(igrid)%surface1(ixmin1-1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3),&
   pgeo(igrid)%surface2(ixmin1:ixmax1,ixmin2-1:ixmax2,ixmin3:ixmax3),&
   pgeo(igrid)%surface3(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3-1:ixmax3),&
    pgeo(igrid)%dvolume(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
   ixGextmin3:ixGextmax3), pgeo(igrid)%dx(ixGextmin1:ixGextmax1,&
   ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,1:ndim), pgeo(igrid)%xbar&
   (ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,1:ndim))
! Covariant coordinates rely on the metric datastructure
if (covariant) then
     call allocate_metric(pgeo(igrid)%m,ixGextmin1,ixGextmin2,ixGextmin3,&
        ixGextmax1,ixGextmax2,ixGextmax3)
      call allocate_metric(pgeo(igrid)%mSurface1,ixGextmin1-1,ixGextmin2-1,&
         ixGextmin3-1,ixGextmax1,ixGextmax2,ixGextmax3,need_derivs=.false.)
       call allocate_metric(pgeo(igrid)%mSurface2,ixGextmin1-1,ixGextmin2-1,&
          ixGextmin3-1,ixGextmax1,ixGextmax2,ixGextmax3,need_derivs=.false.)
       call allocate_metric(pgeo(igrid)%mSurface3,ixGextmin1-1,ixGextmin2-1,&
          ixGextmin3-1,ixGextmax1,ixGextmax2,ixGextmax3,need_derivs=.false.)
end if

dx1=rnode(rpdx1_,igrid);dx2=rnode(rpdx2_,igrid);dx3=rnode(rpdx3_,igrid);
xmin1=rnode(rpxmin1_,igrid);xmin2=rnode(rpxmin2_,igrid)
xmin3=rnode(rpxmin3_,igrid);
call fillgeo(pgeo(igrid),ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixGextmin1,&
   ixGextmin2,ixGextmin3,ixGextmax1,ixGextmax2,ixGextmax3,xmin1,xmin2,xmin3,&
   dx1,dx2,dx3,.false.)

if (errorestimate==1) then
   ixCoGmin1=1;ixCoGmin2=1;ixCoGmin3=1; ixCoGmax1=ixGhi1/2+dixB
   ixCoGmax2=ixGhi2/2+dixB;ixCoGmax3=ixGhi3/2+dixB;
   if (2*int(dixB/2)==dixB) then
      ixGextmin1=ixCoGmin1;ixGextmin2=ixCoGmin2;ixGextmin3=ixCoGmin3
      ixGextmax1=ixCoGmax1;ixGextmax2=ixCoGmax2;ixGextmax3=ixCoGmax3;
   else
      ixGextmin1=ixCoGmin1-1;ixGextmin2=ixCoGmin2-1;ixGextmin3=ixCoGmin3-1
      ixGextmax1=ixCoGmax1+1;ixGextmax2=ixCoGmax2+1;ixGextmax3=ixCoGmax3+1;
   end if
   ixCoMmin1=ixCoGmin1+dixB;ixCoMmin2=ixCoGmin2+dixB;ixCoMmin3=ixCoGmin3+dixB
   ixCoMmax1=ixCoGmax1-dixB;ixCoMmax2=ixCoGmax2-dixB;ixCoMmax3=ixCoGmax3-dixB;
   ixComin1=ixCoMmin1-1;ixComin2=ixCoMmin2-1;ixComin3=ixCoMmin3-1
   ixComax1=ixCoMmax1+1;ixComax2=ixCoMmax2+1;ixComax3=ixCoMmax3+1;
   ixComin1=ixCoGmin1+1;ixComin2=ixCoGmin2+1;ixComin3=ixCoGmin3+1
   ixComax1=ixCoGmax1-1;ixComax2=ixCoGmax2-1;ixComax3=ixCoGmax3-1;

   allocate(pgeoCoarse(igrid)%surfaceC1(ixComin1-1:ixComax1,ixComin2:ixComax2,&
      ixComin3:ixComax3),pgeoCoarse(igrid)%surfaceC2(ixComin1:ixComax1,&
      ixComin2-1:ixComax2,ixComin3:ixComax3),pgeoCoarse(igrid)%surfaceC3&
      (ixComin1:ixComax1,ixComin2:ixComax2,ixComin3-1:ixComax3),&
       pgeoCoarse(igrid)%surface1(ixComin1-1:ixComax1,ixComin2:ixComax2,&
      ixComin3:ixComax3),pgeoCoarse(igrid)%surface2(ixComin1:ixComax1,&
      ixComin2-1:ixComax2,ixComin3:ixComax3),pgeoCoarse(igrid)%surface3&
      (ixComin1:ixComax1,ixComin2:ixComax2,ixComin3-1:ixComax3),&
       pgeoCoarse(igrid)%dvolume(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
      ixGextmin3:ixGextmax3), pgeoCoarse(igrid)%dx(ixGextmin1:ixGextmax1,&
      ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,1:ndim), &
      pgeoCoarse(igrid)%xbar(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
      ixGextmin3:ixGextmax3,1:ndim))

   if (covariant) then
      call allocate_metric(pgeoCoarse(igrid)%m,ixGextmin1,ixGextmin2,&
         ixGextmin3,ixGextmax1,ixGextmax2,ixGextmax3)
       call allocate_metric(pgeoCoarse(igrid)%mSurface1,ixGextmin1-1,&
          ixGextmin2-1,ixGextmin3-1,ixGextmax1,ixGextmax2,ixGextmax3,&
          need_derivs=.false.)
        call allocate_metric(pgeoCoarse(igrid)%mSurface2,ixGextmin1-1,&
           ixGextmin2-1,ixGextmin3-1,ixGextmax1,ixGextmax2,ixGextmax3,&
           need_derivs=.false.)
        call allocate_metric(pgeoCoarse(igrid)%mSurface3,ixGextmin1-1,&
           ixGextmin2-1,ixGextmin3-1,ixGextmax1,ixGextmax2,ixGextmax3,&
           need_derivs=.false.)
   end if

   dx1=two*rnode(rpdx1_,igrid);dx2=two*rnode(rpdx2_,igrid)
   dx3=two*rnode(rpdx3_,igrid);
   call fillgeo(pgeoCoarse(igrid),ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,&
      ixCoGmax2,ixCoGmax3,ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,&
      ixGextmax2,ixGextmax3,xmin1,xmin2,xmin3,dx1,dx2,dx3,.false.)

   ixCoCoGmin1=1;ixCoCoGmin2=1;ixCoCoGmin3=1; ixCoCoGmax1=ixCoGmax1/2+dixB
   ixCoCoGmax2=ixCoGmax2/2+dixB;ixCoCoGmax3=ixCoGmax3/2+dixB;
   if (2*int(dixB/2)==dixB) then
      ixGextmin1=ixCoCoGmin1;ixGextmin2=ixCoCoGmin2;ixGextmin3=ixCoCoGmin3
      ixGextmax1=ixCoCoGmax1;ixGextmax2=ixCoCoGmax2;ixGextmax3=ixCoCoGmax3;
   else
      ixGextmin1=ixCoCoGmin1-1;ixGextmin2=ixCoCoGmin2-1
      ixGextmin3=ixCoCoGmin3-1;ixGextmax1=ixCoCoGmax1+1
      ixGextmax2=ixCoCoGmax2+1;ixGextmax3=ixCoCoGmax3+1;
   end if


   allocate(pgeoCoCo(igrid)%dvolume(ixGextmin1:ixGextmax1,&
      ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3), pgeoCoCo(igrid)%xbar&
      (ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
      1:ndim))
   if (covariant) call allocate_metric(pgeoCoCo(igrid)%m,ixGextmin1,&
      ixGextmin2,ixGextmin3,ixGextmax1,ixGextmax2,ixGextmax3)

   dx1=4.0d0*rnode(rpdx1_,igrid);dx2=4.0d0*rnode(rpdx2_,igrid)
   dx3=4.0d0*rnode(rpdx3_,igrid);
   call fillgeo(pgeoCoCo(igrid),ixCoCoGmin1,ixCoCoGmin2,ixCoCoGmin3,&
      ixCoCoGmax1,ixCoCoGmax2,ixCoCoGmax3,ixGextmin1,ixGextmin2,ixGextmin3,&
      ixGextmax1,ixGextmax2,ixGextmax3,xmin1,xmin2,xmin3,dx1,dx2,dx3,.true.)
else
   ixCoGmin1=1;ixCoGmin2=1;ixCoGmin3=1; ixCoGmax1=ixGhi1/2+dixB
   ixCoGmax2=ixGhi2/2+dixB;ixCoGmax3=ixGhi3/2+dixB;
   if (2*int(dixB/2)==dixB) then
      ixGextmin1=ixCoGmin1;ixGextmin2=ixCoGmin2;ixGextmin3=ixCoGmin3
      ixGextmax1=ixCoGmax1;ixGextmax2=ixCoGmax2;ixGextmax3=ixCoGmax3;
   else
      ixGextmin1=ixCoGmin1-1;ixGextmin2=ixCoGmin2-1;ixGextmin3=ixCoGmin3-1
      ixGextmax1=ixCoGmax1+1;ixGextmax2=ixCoGmax2+1;ixGextmax3=ixCoGmax3+1;
   end if


   ixCoMmin1=ixCoGmin1+dixB;ixCoMmin2=ixCoGmin2+dixB;ixCoMmin3=ixCoGmin3+dixB
   ixCoMmax1=ixCoGmax1-dixB;ixCoMmax2=ixCoGmax2-dixB;ixCoMmax3=ixCoGmax3-dixB;
   ixComin1=ixCoMmin1-1;ixComin2=ixCoMmin2-1;ixComin3=ixCoMmin3-1
   ixComax1=ixCoMmax1+1;ixComax2=ixCoMmax2+1;ixComax3=ixCoMmax3+1;
   ixComin1=ixCoGmin1+1;ixComin2=ixCoGmin2+1;ixComin3=ixCoGmin3+1
   ixComax1=ixCoGmax1-1;ixComax2=ixCoGmax2-1;ixComax3=ixCoGmax3-1;

   allocate(pgeoCoarse(igrid)%surfaceC1(ixCoGmin1-1:ixCoGmax1,&
      ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3),pgeoCoarse(igrid)%surfaceC2&
      (ixCoGmin1:ixCoGmax1,ixCoGmin2-1:ixCoGmax2,ixCoGmin3:ixCoGmax3),&
      pgeoCoarse(igrid)%surfaceC3(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,&
      ixCoGmin3-1:ixCoGmax3),&
        pgeoCoarse(igrid)%surface1(ixComin1-1:ixComax1,ixComin2:ixComax2,&
           ixComin3:ixComax3),pgeoCoarse(igrid)%surface2(ixComin1:ixComax1,&
           ixComin2-1:ixComax2,ixComin3:ixComax3),pgeoCoarse(igrid)%surface3&
           (ixComin1:ixComax1,ixComin2:ixComax2,ixComin3-1:ixComax3),&
        pgeoCoarse(igrid)%dvolume(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
           ixGextmin3:ixGextmax3),&
        pgeoCoarse(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
           ixGextmin3:ixGextmax3,1:ndim),&
        pgeoCoarse(igrid)%xbar(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
           ixGextmin3:ixGextmax3,1:ndim))

   if (covariant) then
      call allocate_metric(pgeoCoarse(igrid)%m,ixGextmin1,ixGextmin2,&
         ixGextmin3,ixGextmax1,ixGextmax2,ixGextmax3)
       call allocate_metric(pgeoCoarse(igrid)%mSurface1,ixGextmin1-1,&
          ixGextmin2-1,ixGextmin3-1,ixGextmax1,ixGextmax2,ixGextmax3,&
          need_derivs=.false.)
 call allocate_metric(pgeoCoarse(igrid)%mSurface2,ixGextmin1-1,ixGextmin2-1,&
    ixGextmin3-1,ixGextmax1,ixGextmax2,ixGextmax3,need_derivs=.false.)
 call allocate_metric(pgeoCoarse(igrid)%mSurface3,ixGextmin1-1,ixGextmin2-1,&
    ixGextmin3-1,ixGextmax1,ixGextmax2,ixGextmax3,need_derivs=.false.)
   end if

   dx1=two*rnode(rpdx1_,igrid);dx2=two*rnode(rpdx2_,igrid)
   dx3=two*rnode(rpdx3_,igrid);
   call fillgeo(pgeoCoarse(igrid),ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,&
      ixCoGmax2,ixCoGmax3,ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,&
      ixGextmax2,ixGextmax3,xmin1,xmin2,xmin3,dx1,dx2,dx3,.false.)


end if

end subroutine getgridgeo
!=============================================================================
subroutine putgridgeo(igrid)

  use mod_metric
  use mod_amrvacdef

integer, intent(in) :: igrid
!-----------------------------------------------------------------------------
if (covariant) then
   call deallocate_metric(pgeo(igrid)%m)
   call deallocate_metric(pgeoCoarse(igrid)%m)
    call deallocate_metric(pgeo(igrid)%mSurface1)
     call deallocate_metric(pgeo(igrid)%mSurface2)
     call deallocate_metric(pgeo(igrid)%mSurface3)
end if

if (errorestimate==1) then
   deallocate(pgeo(igrid)%surfaceC1,pgeo(igrid)%surfaceC2,pgeo&
      (igrid)%surfaceC3,pgeo(igrid)%surface1,pgeo(igrid)%surface2,&
      pgeo(igrid)%surface3,pgeo(igrid)%dvolume,pgeo(igrid)%dx,&
       pgeo(igrid)%xbar, pgeoCoarse(igrid)%surfaceC1,pgeoCoarse&
      (igrid)%surfaceC2,pgeoCoarse(igrid)%surfaceC3,pgeoCoarse&
      (igrid)%surface1,pgeoCoarse(igrid)%surface2,pgeoCoarse(igrid)%surface3,&
      pgeoCoarse(igrid)%dvolume,pgeoCoarse(igrid)%dx,pgeoCoarse(igrid)%xbar,&
      pgeoCoCo(igrid)%dvolume,pgeoCoCo(igrid)%xbar)
   
   if (covariant) then
      call deallocate_metric(pgeoCoCo(igrid)%m)
       call deallocate_metric(pgeoCoarse(igrid)%mSurface1)
        call deallocate_metric(pgeoCoarse(igrid)%mSurface2)
        call deallocate_metric(pgeoCoarse(igrid)%mSurface3)
   end if
   
else
   deallocate(pgeo(igrid)%surfaceC1,pgeo(igrid)%surfaceC2,pgeo&
      (igrid)%surfaceC3,pgeo(igrid)%surface1,pgeo(igrid)%surface2,&
      pgeo(igrid)%surface3,pgeo(igrid)%dvolume,pgeo(igrid)%dx,&
      pgeo(igrid)%xbar,pgeoCoarse(igrid)%dvolume,pgeoCoarse(igrid)%xbar)

   deallocate(pgeoCoarse(igrid)%surfaceC1,pgeoCoarse(igrid)%surfaceC2,&
      pgeoCoarse(igrid)%surfaceC3,pgeoCoarse(igrid)%surface1,&
      pgeoCoarse(igrid)%surface2,pgeoCoarse(igrid)%surface3,&
        pgeoCoarse(igrid)%dx)
   if (covariant) then
       call deallocate_metric(pgeoCoarse(igrid)%mSurface1)
 call deallocate_metric(pgeoCoarse(igrid)%mSurface2)
 call deallocate_metric(pgeoCoarse(igrid)%mSurface3)
   end if

end if

end subroutine putgridgeo
!=============================================================================
subroutine fillgeo(pgeogrid,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
   ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,ixGextmax2,ixGextmax3,xmin1,&
   xmin2,xmin3,dx1,dx2,dx3,need_only_volume)

use mod_metric
use mod_amrvacdef

type(geoalloc) :: pgeogrid
integer, intent(in) :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
    ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,ixGextmax2,ixGextmax3
double precision, intent(in) :: xmin1,xmin2,xmin3, dx1,dx2,dx3
logical, intent(in) :: need_only_volume

integer :: idims, ix, ixMmin1,ixMmin2,ixMmin3,ixMmax1,ixMmax2,ixMmax3, ixmin1,&
   ixmin2,ixmin3,ixmax1,ixmax2,ixmax3, ixCmin1,ixCmin2,ixCmin3,ixCmax1,&
   ixCmax2,ixCmax3
double precision :: x(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
   ixGextmin3:ixGextmax3,ndim)
!-----------------------------------------------------------------------------
ixMmin1=ixGmin1+dixB;ixMmin2=ixGmin2+dixB;ixMmin3=ixGmin3+dixB
ixMmax1=ixGmax1-dixB;ixMmax2=ixGmax2-dixB;ixMmax3=ixGmax3-dixB;
!ix^L=ixM^L^LADD1;
ixmin1=ixGmin1+1;ixmin2=ixGmin2+1;ixmin3=ixGmin3+1;ixmax1=ixGmax1-1
ixmax2=ixGmax2-1;ixmax3=ixGmax3-1;

if (covariant) then
   ! call for the general geometry
   call fillgeo_covariant(pgeogrid,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
      ixGmax3,ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,ixGextmax2,&
      ixGextmax3,xmin1,xmin2,xmin3,dx1,dx2,dx3,need_only_volume)
   return
end if

select case (typecoord)
case ("slabtest")

   pgeogrid%dvolume(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
      ixGextmin3:ixGextmax3) = dx1*dx2*dx3

   if (need_only_volume) return

   ixCmin1=ixmin1-kr(1,1);ixCmin2=ixmin2-kr(2,1);ixCmin3=ixmin3-kr(3,1)
   ixCmax1=ixmax1;ixCmax2=ixmax2;ixCmax3=ixmax3;
   pgeogrid%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
      = dx2*dx3
   pgeogrid%surface1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
      = dx2*dx3
   
   ixCmin1=ixmin1-kr(1,2);ixCmin2=ixmin2-kr(2,2);ixCmin3=ixmin3-kr(3,2)
   ixCmax1=ixmax1;ixCmax2=ixmax2;ixCmax3=ixmax3;
   pgeogrid%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=dx1*dx3
   
   ixCmin1=ixmin1-kr(1,2);ixCmin2=ixmin2-kr(2,2);ixCmin3=ixmin3-kr(3,2)
   ixCmax1=ixmax1;ixCmax2=ixmax2;ixCmax3=ixmax3;
   pgeogrid%surface2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=dx1*dx3
   
   ixCmin1=ixmin1-kr(1,3);ixCmin2=ixmin2-kr(2,3);ixCmin3=ixmin3-kr(3,3)
   ixCmax1=ixmax1;ixCmax2=ixmax2;ixCmax3=ixmax3;
   pgeogrid%surfaceC3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=dx1*dx2
   pgeogrid%surface3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=dx1*dx2

   pgeogrid%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
      ixGextmin3:ixGextmax3,1)=dx1
   pgeogrid%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
      ixGextmin3:ixGextmax3,2)=dx2
   pgeogrid%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
      ixGextmin3:ixGextmax3,3)=dx3;
   ! Not yet correct:
   pgeogrid%xbar(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
      ixGextmin3:ixGextmax3,1)=0.0d0
   pgeogrid%xbar(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
      ixGextmin3:ixGextmax3,2)=0.0d0
   pgeogrid%xbar(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
      ixGextmin3:ixGextmax3,3)=0.0d0;


case ("spherical")

   do idims=1,min(ndim,2)
      select case(idims)
      case(1)
         do ix = ixGextmin1,ixGextmax1
            x(ix,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,1)&
               =xmin1+(dble(ix-dixB)-half)*dx1
         end do
      case(2)
         do ix = ixGextmin2,ixGextmax2
            x(ixGextmin1:ixGextmax1,ix,ixGextmin3:ixGextmax3,2)&
               =xmin2+(dble(ix-dixB)-half)*dx2
         end do
      case(3)
         do ix = ixGextmin3,ixGextmax3
            x(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ix,3)&
               =xmin3+(dble(ix-dixB)-half)*dx3
         end do
      end select
   end do

   if(typespherical==0) then
     pgeogrid%dvolume(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
        ixGextmin3:ixGextmax3)=(x(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
        ixGextmin3:ixGextmax3,1)**2+dx1**2/12.0d0)*dx1  *two*dabs(dsin(x&
        (ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
        2)))*dsin(half*dx2)*dx3
   else
     pgeogrid%dvolume(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
        ixGextmin3:ixGextmax3)=(x(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
        ixGextmin3:ixGextmax3,1)**2)*dx1  *dabs(dsin(x(ixGextmin1:ixGextmax1,&
        ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,2)))*dx2*dx3
   endif

   if (need_only_volume) return

   ! Not yet correct:
   pgeogrid%xbar(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
      ixGextmin3:ixGextmax3,1)=0.0d0
   pgeogrid%xbar(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
      ixGextmin3:ixGextmax3,2)=0.0d0
   pgeogrid%xbar(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
      ixGextmin3:ixGextmax3,3)=0.0d0;

   ixCmin1=ixmin1-kr(1,1);ixCmin2=ixmin2-kr(2,1);ixCmin3=ixmin3-kr(3,1)
   ixCmax1=ixmax1;ixCmax2=ixmax2;ixCmax3=ixmax3;
   if(typespherical==0) then
       pgeogrid%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
          =(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)+half*dx1)&
          **2  *two*dsin(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
          2))*dsin(half*dx2)*dx3
   else
       pgeogrid%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
          =(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)+half*dx1)&
          **2  *dsin(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
          2))*dx2*dx3
   endif

   
   ixCmin1=ixmin1-kr(1,2);ixCmin2=ixmin2-kr(2,2);ixCmin3=ixmin3-kr(3,2)
   ixCmax1=ixmax1;ixCmax2=ixmax2;ixCmax3=ixmax3;
   pgeogrid%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
      =x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)*dx1 &
              *dsin(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                 2)+half*dx2)*dx3

   
   ixCmin1=ixmin1-kr(1,3);ixCmin2=ixmin2-kr(2,3);ixCmin3=ixmin3-kr(3,3)
   ixCmax1=ixmax1;ixCmax2=ixmax2;ixCmax3=ixmax3;
   pgeogrid%surfaceC3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
      =x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)*dx1*dx2

   ixCmin1=ixmin1-kr(1,1);ixCmin2=ixmin2-kr(2,1);ixCmin3=ixmin3-kr(3,1)
   ixCmax1=ixmax1;ixCmax2=ixmax2;ixCmax3=ixmax3;
   if(typespherical==0) then
       pgeogrid%surface1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
          =x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)&
          **2  *two*dsin(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
          2))*dsin(half*dx2)*dx3
   else
      pgeogrid%surface1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
         =x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)&
         **2  *dsin(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         2))*dx2*dx3
   endif

   
   ixCmin1=ixmin1-kr(1,2);ixCmin2=ixmin2-kr(2,2);ixCmin3=ixmin3-kr(3,2)
   ixCmax1=ixmax1;ixCmax2=ixmax2;ixCmax3=ixmax3;
   pgeogrid%surface2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
      =x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)*dx1 &
              *dsin(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2))*dx3

   
   ixCmin1=ixmin1-kr(1,3);ixCmin2=ixmin2-kr(2,3);ixCmin3=ixmin3-kr(3,3)
   ixCmax1=ixmax1;ixCmax2=ixmax2;ixCmax3=ixmax3;
   pgeogrid%surface3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
      =x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)*dx1*dx2

   pgeogrid%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
      ixGextmin3:ixGextmax3,1)=dx1
    pgeogrid%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,2)=x(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,1)*dx2
    pgeogrid%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,3)=x(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,1)*dsin(x(ixGextmin1:ixGextmax1,&
       ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,2))*dx3

case ("cylindrical")

   do ix = ixGextmin1,ixGextmax1
      x(ix,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,1)=xmin1+(dble&
         (ix-dixB)-half)*dx1
   end do

   pgeogrid%dvolume(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
      ixGextmin3:ixGextmax3)=dabs(x(ixGextmin1:ixGextmax1,&
      ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,1))*dx1*dx2*dx3
   pgeogrid%dvolume(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
      ixGextmin3:ixGextmax3)=dabs(half*((x(ixGextmin1:ixGextmax1,&
      ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,1)+half*dx1)&
      **2-(x(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
      ixGextmin3:ixGextmax3,1)-half*dx1)**2))*dx2 *dx3

   if (need_only_volume) return

   ! Not yet correct:
   pgeogrid%xbar(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
      ixGextmin3:ixGextmax3,1)=0.0d0
   pgeogrid%xbar(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
      ixGextmin3:ixGextmax3,2)=0.0d0
   pgeogrid%xbar(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
      ixGextmin3:ixGextmax3,3)=0.0d0;

   ixCmin1=ixmin1-kr(1,1);ixCmin2=ixmin2-kr(2,1);ixCmin3=ixmin3-kr(3,1)
   ixCmax1=ixmax1;ixCmax2=ixmax2;ixCmax3=ixmax3;
   !!pgeogrid%surfaceC1(ixC^S)=(x(ixC^S,1)+half*dx1){^DE&*dx^DE }
   pgeogrid%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
      =dabs((x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
      1)+half*dx1))*dx2 *dx3
   
   ixCmin1=ixmin1-kr(1,2);ixCmin2=ixmin2-kr(2,2);ixCmin3=ixmin3-kr(3,2)
   ixCmax1=ixmax1;ixCmax2=ixmax2;ixCmax3=ixmax3;
   if (3==2) pgeogrid%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      ixCmin3:ixCmax3)=x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
      1)*dx1*dx3
   if (2==2) pgeogrid%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      ixCmin3:ixCmax3)=dx1*dx3
   
   ixCmin1=ixmin1-kr(1,3);ixCmin2=ixmin2-kr(2,3);ixCmin3=ixmin3-kr(3,3)
   ixCmax1=ixmax1;ixCmax2=ixmax2;ixCmax3=ixmax3;
   if (3==3) pgeogrid%surfaceC3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      ixCmin3:ixCmax3)=x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
      1)*dx1*dx2
   if (2==3) pgeogrid%surfaceC3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      ixCmin3:ixCmax3)=dx1*dx2

   ixCmin1=ixmin1-kr(1,1);ixCmin2=ixmin2-kr(2,1);ixCmin3=ixmin3-kr(3,1)
   ixCmax1=ixmax1;ixCmax2=ixmax2;ixCmax3=ixmax3;
   !!pgeogrid%surface1(ixC^S)=x(ixC^S,1){^DE&*dx^DE }
   pgeogrid%surface1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
      =dabs(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1))*dx2 *dx3
   
   ixCmin1=ixmin1-kr(1,2);ixCmin2=ixmin2-kr(2,2);ixCmin3=ixmin3-kr(3,2)
   ixCmax1=ixmax1;ixCmax2=ixmax2;ixCmax3=ixmax3;
   if (3==2) pgeogrid%surface2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      ixCmin3:ixCmax3)=x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
      1)*dx1*dx3
   if (2==2) pgeogrid%surface2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      ixCmin3:ixCmax3)=dx1*dx3
   
   ixCmin1=ixmin1-kr(1,3);ixCmin2=ixmin2-kr(2,3);ixCmin3=ixmin3-kr(3,3)
   ixCmax1=ixmax1;ixCmax2=ixmax2;ixCmax3=ixmax3;
   if (3==3) pgeogrid%surface3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      ixCmin3:ixCmax3)=x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
      1)*dx1*dx2
   if (2==3) pgeogrid%surface3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      ixCmin3:ixCmax3)=dx1*dx2


   pgeogrid%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
      ixGextmin3:ixGextmax3,1)=dx1
    if (2==3) pgeogrid%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,2)=dx2
    if (3==3) pgeogrid%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,3)=dx3
    if (2==2) pgeogrid%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,2)=x(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,1)*dx2
    if (3==2) pgeogrid%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,3)=x(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,1)*dx3

case default

   call mpistop("Sorry, typecoord unknown")
   
end select

end subroutine fillgeo
!=============================================================================
subroutine gradient(q,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idir,gradq)

! Calculate gradient of a scalar q within ixL in direction idir

use mod_amrvacdef

integer :: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3, idir
double precision :: q(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
    gradq(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)

double precision :: qC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),invdx
integer :: jxmin1,jxmin2,jxmin3,jxmax1,jxmax2,jxmax3, hxmin1,hxmin2,hxmin3,&
   hxmax1,hxmax2,hxmax3, ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,&
    jxCmin1,jxCmin2,jxCmin3,jxCmax1,jxCmax2,jxCmax3 

!-----------------------------------------------------------------------------

invdx=1.d0/dxlevel(idir)
if (slab) then

   jxmin1=ixmin1+kr(idir,1);jxmin2=ixmin2+kr(idir,2);jxmin3=ixmin3+kr(idir,3)
   jxmax1=ixmax1+kr(idir,1);jxmax2=ixmax2+kr(idir,2);jxmax3=ixmax3+kr(idir,3);
   hxmin1=ixmin1-kr(idir,1);hxmin2=ixmin2-kr(idir,2);hxmin3=ixmin3-kr(idir,3)
   hxmax1=ixmax1-kr(idir,1);hxmax2=ixmax2-kr(idir,2);hxmax3=ixmax3-kr(idir,3);
   gradq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3) = half*(q(jxmin1:jxmax1,&
      jxmin2:jxmax2,jxmin3:jxmax3)-q(hxmin1:hxmax1,hxmin2:hxmax2,&
      hxmin3:hxmax3))*invdx

else
   hxmin1=ixmin1-kr(idir,1);hxmin2=ixmin2-kr(idir,2);hxmin3=ixmin3-kr(idir,3)
   hxmax1=ixmax1-kr(idir,1);hxmax2=ixmax2-kr(idir,2);hxmax3=ixmax3-kr(idir,3);
   ixCmin1=hxmin1;ixCmin2=hxmin2;ixCmin3=hxmin3;ixCmax1=ixmax1;ixCmax2=ixmax2
   ixCmax3=ixmax3;
   jxCmin1=ixCmin1+kr(idir,1);jxCmin2=ixCmin2+kr(idir,2)
   jxCmin3=ixCmin3+kr(idir,3);jxCmax1=ixCmax1+kr(idir,1)
   jxCmax2=ixCmax2+kr(idir,2);jxCmax3=ixCmax3+kr(idir,3);
   select case(idir)
   case(1)
      qC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=mygeo%surfaceC1&
         (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)*half*(q&
         (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)+q(jxCmin1:jxCmax1,&
         jxCmin2:jxCmax2,jxCmin3:jxCmax3))
      gradq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=(qC(ixmin1:ixmax1,&
         ixmin2:ixmax2,ixmin3:ixmax3)-qC(hxmin1:hxmax1,hxmin2:hxmax2,&
         hxmin3:hxmax3))/mygeo%dvolume(ixmin1:ixmax1,ixmin2:ixmax2,&
         ixmin3:ixmax3)
      ! Substract difference divergence and gradient
      gradq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=gradq(ixmin1:ixmax1,&
         ixmin2:ixmax2,ixmin3:ixmax3)-q(ixmin1:ixmax1,ixmin2:ixmax2,&
         ixmin3:ixmax3) &
                     *(mygeo%surfaceC1(ixmin1:ixmax1,ixmin2:ixmax2,&
                        ixmin3:ixmax3)-mygeo%surfaceC1(hxmin1:hxmax1,&
                        hxmin2:hxmax2,hxmin3:hxmax3)) &
                    /mygeo%dvolume(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3) 
   case(2)
      qC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=mygeo%surfaceC2&
         (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)*half*(q&
         (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)+q(jxCmin1:jxCmax1,&
         jxCmin2:jxCmax2,jxCmin3:jxCmax3))
      gradq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=(qC(ixmin1:ixmax1,&
         ixmin2:ixmax2,ixmin3:ixmax3)-qC(hxmin1:hxmax1,hxmin2:hxmax2,&
         hxmin3:hxmax3))/mygeo%dvolume(ixmin1:ixmax1,ixmin2:ixmax2,&
         ixmin3:ixmax3)
      ! Substract difference divergence and gradient
      gradq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=gradq(ixmin1:ixmax1,&
         ixmin2:ixmax2,ixmin3:ixmax3)-q(ixmin1:ixmax1,ixmin2:ixmax2,&
         ixmin3:ixmax3) &
                     *(mygeo%surfaceC2(ixmin1:ixmax1,ixmin2:ixmax2,&
                        ixmin3:ixmax3)-mygeo%surfaceC2(hxmin1:hxmax1,&
                        hxmin2:hxmax2,hxmin3:hxmax3)) &
                    /mygeo%dvolume(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3) 
   case(3)
      qC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=mygeo%surfaceC3&
         (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)*half*(q&
         (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)+q(jxCmin1:jxCmax1,&
         jxCmin2:jxCmax2,jxCmin3:jxCmax3))
      gradq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=(qC(ixmin1:ixmax1,&
         ixmin2:ixmax2,ixmin3:ixmax3)-qC(hxmin1:hxmax1,hxmin2:hxmax2,&
         hxmin3:hxmax3))/mygeo%dvolume(ixmin1:ixmax1,ixmin2:ixmax2,&
         ixmin3:ixmax3)
      ! Substract difference divergence and gradient
      gradq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=gradq(ixmin1:ixmax1,&
         ixmin2:ixmax2,ixmin3:ixmax3)-q(ixmin1:ixmax1,ixmin2:ixmax2,&
         ixmin3:ixmax3) &
                     *(mygeo%surfaceC3(ixmin1:ixmax1,ixmin2:ixmax2,&
                        ixmin3:ixmax3)-mygeo%surfaceC3(hxmin1:hxmax1,&
                        hxmin2:hxmax2,hxmin3:hxmax3)) &
                    /mygeo%dvolume(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3) 
   end select
end if

end subroutine gradient
!=============================================================================
subroutine upwindGradientS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idir,w,x,var,gradient)

  use mod_limiter
  use mod_amrvacdef
  
  integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idir
  double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:nw)
  double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:ndim)
  double precision, intent(in)    :: var(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3)
  double precision, intent(out)   :: gradient(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3)
  ! .. local ..
  integer                                 :: hxOmin1,hxOmin2,hxOmin3,hxOmax1,&
     hxOmax2,hxOmax3, ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,&
      jxCmin1,jxCmin2,jxCmin3,jxCmax1,jxCmax2,jxCmax3, gxCmin1,gxCmin2,&
     gxCmin3,gxCmax1,gxCmax2,gxCmax3, hxCmin1,hxCmin2,hxCmin3,hxCmax1,hxCmax2,&
     hxCmax3
  double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
     1:ndim)             ::  xi
  double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
     1:nw) :: wLC, wRC, wLp, wRp, wprim
  double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)      &
     :: cmaxC, cmaxRC, cmaxLC
  double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)      &
     :: cminC, cminRC, cminLC
  integer, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)               &
     :: patchf
  double precision,dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3)       :: qL, qR, q, ldq, rdq, dqC
  !-----------------------------------------------------------------------------

  ! Always primitive reconstruction:
  wprim(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw) &
     = w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw)
  call primitive(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixImin1,&
     ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,wprim,x) !also applies floors to wprim


  ! Index ranges:
  hxOmin1=ixOmin1-kr(idir,1);hxOmin2=ixOmin2-kr(idir,2)
  hxOmin3=ixOmin3-kr(idir,3);hxOmax1=ixOmax1-kr(idir,1)
  hxOmax2=ixOmax2-kr(idir,2);hxOmax3=ixOmax3-kr(idir,3);
  ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3; ixCmin1=hxOmin1
  ixCmin2=hxOmin2;ixCmin3=hxOmin3;
  jxCmin1=ixCmin1+kr(idir,1);jxCmin2=ixCmin2+kr(idir,2)
  jxCmin3=ixCmin3+kr(idir,3);jxCmax1=ixCmax1+kr(idir,1)
  jxCmax2=ixCmax2+kr(idir,2);jxCmax3=ixCmax3+kr(idir,3);
  gxCmin1=ixCmin1-kr(idir,1);gxCmin2=ixCmin2-kr(idir,2)
  gxCmin3=ixCmin3-kr(idir,3);gxCmax1=jxCmax1;gxCmax2=jxCmax2;gxCmax3=jxCmax3;
  hxCmin1=gxCmin1+kr(idir,1);hxCmin2=gxCmin2+kr(idir,2)
  hxCmin3=gxCmin3+kr(idir,3);hxCmax1=gxCmax1+kr(idir,1)
  hxCmax2=gxCmax2+kr(idir,2);hxCmax3=gxCmax3+kr(idir,3);
  
  !==================================================
  ! left and right states:
  !==================================================
  wRp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1:nwflux)&
     =wprim(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,1:nwflux)
  wLp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1:nwflux)&
     =wprim(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1:nwflux)
  xi(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir) &
     = half * ( x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
     idir)+x(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,idir) )

  call upwindLR(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,&
     ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,&
     ixCmax2,ixCmax3,idir,wprim,wprim,wLC,wRC,wLp,wRp,x,dxlevel(idir))
  ! get auxilaries for L and R states
  if (nwaux>0) then
     call getaux(.true.,wLC,xi,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
        ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,'upwindGradientS-L')
     call getaux(.true.,wRC,xi,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
        ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,'upwindGradientS-R')
  end if
  call getcmax(wLC,xi,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixCmin1,&
     ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idir,cmaxLC,cminLC,.true.)
  call getcmax(wRC,xi,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixCmin1,&
     ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idir,cmaxRC,cminRC,.true.)
  ! now take the maximum of left and right states
  cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=max(cmaxRC&
     (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3),cmaxLC(ixCmin1:ixCmax1,&
     ixCmin2:ixCmax2,ixCmin3:ixCmax3))
  cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=min(cminRC&
     (ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3),cminLC(ixCmin1:ixCmax1,&
     ixCmin2:ixCmax2,ixCmin3:ixCmax3))

  ! select cases based on fasted waves:
  patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) =  1
  where(cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) >= zero)
     patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = -2
  elsewhere(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) <= zero)
     patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) =  2
  endwhere
  !==================================================
  !==================================================

  
  !==================================================
  ! reconstruct the variable itself, from left and right
  !==================================================
  qR(gxCmin1:gxCmax1,gxCmin2:gxCmax2,gxCmin3:gxCmax3) = var(hxCmin1:hxCmax1,&
     hxCmin2:hxCmax2,hxCmin3:hxCmax3)
  qL(gxCmin1:gxCmax1,gxCmin2:gxCmax2,gxCmin3:gxCmax3) = var(gxCmin1:gxCmax1,&
     gxCmin2:gxCmax2,gxCmin3:gxCmax3)

  select case (typegradlimiter)
     
  case (limiter_ppm)
     
     call PPMlimitervar(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,idir,var,var,qL,qR)

  case default

     dqC(gxCmin1:gxCmax1,gxCmin2:gxCmax2,gxCmin3:gxCmax3)= qR(gxCmin1:gxCmax1,&
        gxCmin2:gxCmax2,gxCmin3:gxCmax3)-qL(gxCmin1:gxCmax1,gxCmin2:gxCmax2,&
        gxCmin3:gxCmax3)

     call dwlimiter2(dqC,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        gxCmin1,gxCmin2,gxCmin3,gxCmax1,gxCmax2,gxCmax3,idir,typegradlimiter,&
        ldq,rdq)

     qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = qL(ixCmin1:ixCmax1,&
        ixCmin2:ixCmax2,ixCmin3:ixCmax3) + half*ldq(ixCmin1:ixCmax1,&
        ixCmin2:ixCmax2,ixCmin3:ixCmax3)
     qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = qR(ixCmin1:ixCmax1,&
        ixCmin2:ixCmax2,ixCmin3:ixCmax3) - half*rdq(jxCmin1:jxCmax1,&
        jxCmin2:jxCmax2,jxCmin3:jxCmax3)

  end select
  !==================================================


  ! Sl>=0, use left interface value:
  where (patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) .eq. -2)
     q(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = qL(ixCmin1:ixCmax1,&
        ixCmin2:ixCmax2,ixCmin3:ixCmax3)     
  ! Sr<=0, use right interface value:
  elsewhere  (patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) .eq. 2)
     q(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = qR(ixCmin1:ixCmax1,&
        ixCmin2:ixCmax2,ixCmin3:ixCmax3)
  ! Outgoing in both directions, use average value:
  elsewhere
     q(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = half * &
        (qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) + &
        qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3))
  end where

  if (slab) then
     gradient(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
        =half*(q(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)-q&
        (hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3))/dxlevel(idir)
  else
     select case(idir)
        case(1)
        gradient(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
           =(q(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)-q&
           (hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3))&
           /mygeo%dx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir) 
        case(2)
        gradient(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
           =(q(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)-q&
           (hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3))&
           /mygeo%dx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir) 
        case(3)
        gradient(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
           =(q(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)-q&
           (hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3))&
           /mygeo%dx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir) 
     end select
  end if

end subroutine upwindGradientS
!=============================================================================
subroutine gradientS(q,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idir,gradq)

! Calculate gradient of a scalar q within ixL in direction idir
! first use limiter to go from cell center to edge

use mod_limiter
use mod_amrvacdef

integer :: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3, idir
double precision :: q(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
    gradq(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
double precision :: dxdim

double precision :: qC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
double precision,dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3):: qL,qR,&
   dqC,ldq,rdq
integer                          :: hxmin1,hxmin2,hxmin3,hxmax1,hxmax2,hxmax3,&
   ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,jxCmin1,jxCmin2,jxCmin3,&
   jxCmax1,jxCmax2,jxCmax3,gxCmin1,gxCmin2,gxCmin3,gxCmax1,gxCmax2,gxCmax3,&
   hxCmin1,hxCmin2,hxCmin3,hxCmax1,hxCmax2,hxCmax3
!-----------------------------------------------------------------------------

hxmin1=ixmin1-kr(idir,1);hxmin2=ixmin2-kr(idir,2);hxmin3=ixmin3-kr(idir,3)
hxmax1=ixmax1-kr(idir,1);hxmax2=ixmax2-kr(idir,2);hxmax3=ixmax3-kr(idir,3);
ixCmin1=hxmin1;ixCmin2=hxmin2;ixCmin3=hxmin3;ixCmax1=ixmax1;ixCmax2=ixmax2
ixCmax3=ixmax3;
jxCmin1=ixCmin1+kr(idir,1);jxCmin2=ixCmin2+kr(idir,2)
jxCmin3=ixCmin3+kr(idir,3);jxCmax1=ixCmax1+kr(idir,1)
jxCmax2=ixCmax2+kr(idir,2);jxCmax3=ixCmax3+kr(idir,3);
gxCmin1=ixCmin1-kr(idir,1);gxCmin2=ixCmin2-kr(idir,2)
gxCmin3=ixCmin3-kr(idir,3);gxCmax1=jxCmax1;gxCmax2=jxCmax2;gxCmax3=jxCmax3;
hxCmin1=gxCmin1+kr(idir,1);hxCmin2=gxCmin2+kr(idir,2)
hxCmin3=gxCmin3+kr(idir,3);hxCmax1=gxCmax1+kr(idir,1)
hxCmax2=gxCmax2+kr(idir,2);hxCmax3=gxCmax3+kr(idir,3);


!==================================================
! reconstruct the variable itself, from left and right
!==================================================
qR(gxCmin1:gxCmax1,gxCmin2:gxCmax2,gxCmin3:gxCmax3) = q(hxCmin1:hxCmax1,&
   hxCmin2:hxCmax2,hxCmin3:hxCmax3)
qL(gxCmin1:gxCmax1,gxCmin2:gxCmax2,gxCmin3:gxCmax3) = q(gxCmin1:gxCmax1,&
   gxCmin2:gxCmax2,gxCmin3:gxCmax3)

select case (typegradlimiter)

case (limiter_ppm)

   call PPMlimitervar(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixMlo1,ixMlo2,&
      ixMlo3,ixMhi1,ixMhi2,ixMhi3,idir,q,q,qL,qR)

case default

   dqC(gxCmin1:gxCmax1,gxCmin2:gxCmax2,gxCmin3:gxCmax3)= qR(gxCmin1:gxCmax1,&
      gxCmin2:gxCmax2,gxCmin3:gxCmax3)-qL(gxCmin1:gxCmax1,gxCmin2:gxCmax2,&
      gxCmin3:gxCmax3)

   call dwlimiter2(dqC,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,gxCmin1,&
      gxCmin2,gxCmin3,gxCmax1,gxCmax2,gxCmax3,idir,typegradlimiter,ldq,rdq)

   qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = qL(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3) + half*ldq(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3)
   qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = qR(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3) - half*rdq(jxCmin1:jxCmax1,&
      jxCmin2:jxCmax2,jxCmin3:jxCmax3)

end select
!==================================================


if (slab) then
   gradq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=half*(qR(ixmin1:ixmax1,&
      ixmin2:ixmax2,ixmin3:ixmax3)-qL(hxmin1:hxmax1,hxmin2:hxmax2,&
      hxmin3:hxmax3))/dxlevel(idir)
else
   select case(idir)
   case(1)
    gradq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=(qR(ixmin1:ixmax1,&
       ixmin2:ixmax2,ixmin3:ixmax3)-qL(hxmin1:hxmax1,hxmin2:hxmax2,&
       hxmin3:hxmax3))/mygeo%dx(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       idir) 
   case(2)
    gradq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=(qR(ixmin1:ixmax1,&
       ixmin2:ixmax2,ixmin3:ixmax3)-qL(hxmin1:hxmax1,hxmin2:hxmax2,&
       hxmin3:hxmax3))/mygeo%dx(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       idir) 
   case(3)
    gradq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=(qR(ixmin1:ixmax1,&
       ixmin2:ixmax2,ixmin3:ixmax3)-qL(hxmin1:hxmax1,hxmin2:hxmax2,&
       hxmin3:hxmax3))/mygeo%dx(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       idir) 
   end select
end if

end subroutine gradientS
!=============================================================================
subroutine divvector(qvec,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divq)

! Calculate divergence of a vector qvec within ixL

use mod_amrvacdef

integer :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,&
   ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision :: qvec(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:ndir),&
    divq(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)

double precision :: qC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
    invdx(1:ndim)
integer :: jxOmin1,jxOmin2,jxOmin3,jxOmax1,jxOmax2,jxOmax3, hxOmin1,hxOmin2,&
   hxOmin3,hxOmax1,hxOmax2,hxOmax3, ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,&
   ixCmax3, jxCmin1,jxCmin2,jxCmin3,jxCmax1,jxCmax2,jxCmax3, idims, ixmin1,&
   ixmin2,ixmin3,ixmax1,ixmax2,ixmax3 
!-----------------------------------------------------------------------------

ixmin1=ixOmin1-1;ixmin2=ixOmin2-1;ixmin3=ixOmin3-1;ixmax1=ixOmax1+1
ixmax2=ixOmax2+1;ixmax3=ixOmax3+1;

if (ixImin1>ixmin1.or.ixImax1<ixmax1.or.ixImin2>ixmin2.or.ixImax2<ixmax2&
   .or.ixImin3>ixmin3.or.ixImax3<ixmax3) call mpistop("Error in divvector: &
   Non-conforming input limits")
invdx=1.d0/dxlevel
divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=zero
do idims=1,ndim
   if (slab) then

     jxOmin1=ixOmin1+kr(idims,1);jxOmin2=ixOmin2+kr(idims,2)
     jxOmin3=ixOmin3+kr(idims,3);jxOmax1=ixOmax1+kr(idims,1)
     jxOmax2=ixOmax2+kr(idims,2);jxOmax3=ixOmax3+kr(idims,3);
     hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
     hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
     hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);
     divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=divq&
        (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)+half*(qvec&
        (jxOmin1:jxOmax1,jxOmin2:jxOmax2,jxOmin3:jxOmax3,idims)-qvec&
        (hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,idims))*invdx(idims)

   else
     hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
     hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
     hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);
     ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmin3=hxOmin3;ixCmax1=ixOmax1
     ixCmax2=ixOmax2;ixCmax3=ixOmax3;
     jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
     jxCmin3=ixCmin3+kr(idims,3);jxCmax1=ixCmax1+kr(idims,1)
     jxCmax2=ixCmax2+kr(idims,2);jxCmax3=ixCmax3+kr(idims,3);
     select case(idims)
     case(1)
        qC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
           =mygeo%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*half*(qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,idims)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
           jxCmin3:jxCmax3,idims))
        divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
           =divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)+(qC&
           (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)-qC&
           (hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3))&
           /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) 
     case(2)
        qC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
           =mygeo%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*half*(qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,idims)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
           jxCmin3:jxCmax3,idims))
        divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
           =divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)+(qC&
           (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)-qC&
           (hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3))&
           /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) 
     case(3)
        qC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
           =mygeo%surfaceC3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*half*(qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,idims)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
           jxCmin3:jxCmax3,idims))
        divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
           =divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)+(qC&
           (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)-qC&
           (hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3))&
           /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) 
      end select
   end if
end do


end subroutine divvector 
!=============================================================================
subroutine curl3(qvec,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
   ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,curlvec)
! To calculate the curl 'curlvec' of a covariant 3-vector 'qvec'
! in the range ixO^L. The total extent ixI^L must be at least 1
! cell bigger than ixO^L in all directions.
! The curl is an array of three components.
! A pointer myM must be associated to the metric.
! Derivatives are calculated with a 2nd order approximation and
! without staggering.
! This routine works only for vectors in coordinate basis.

use mod_amrvacdef

integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(in)  :: qvec(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:3)
double precision, intent(out) :: curlvec(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:3)

! ... local ...

integer :: idir1,idir2,idir3
integer :: jxOmin1,jxOmin2,jxOmin3,jxOmax1,jxOmax2,jxOmax3,hxOmin1,hxOmin2,&
   hxOmin3,hxOmax1,hxOmax2,hxOmax3,ixOpmin1,ixOpmin2,ixOpmin3,ixOpmax1,&
   ixOpmax2,ixOpmax3,ixOmmin1,ixOmmin2,ixOmmin3,ixOmmax1,ixOmmax2,ixOmmax3
double precision :: invdx(1:3)
!-----------------------------------------------------------------------------

invdx=1.d0/dxlevel

curlvec=zero

do idir1=1,3
  ! Directions
  idir2=modulo(idir1,3)+1
  idir3=modulo(idir2,3)+1

  ! Indices
  hxOmin1=ixOmin1-kr(idir2,1);hxOmin2=ixOmin2-kr(idir2,2)
  hxOmin3=ixOmin3-kr(idir2,3);hxOmax1=ixOmax1-kr(idir2,1)
  hxOmax2=ixOmax2-kr(idir2,2);hxOmax3=ixOmax3-kr(idir2,3);
  jxOmin1=ixOmin1+kr(idir2,1);jxOmin2=ixOmin2+kr(idir2,2)
  jxOmin3=ixOmin3+kr(idir2,3);jxOmax1=ixOmax1+kr(idir2,1)
  jxOmax2=ixOmax2+kr(idir2,2);jxOmax3=ixOmax3+kr(idir2,3);
  ixOmmin1=ixOmin1-kr(idir3,1);ixOmmin2=ixOmin2-kr(idir3,2)
  ixOmmin3=ixOmin3-kr(idir3,3);ixOmmax1=ixOmax1-kr(idir3,1)
  ixOmmax2=ixOmax2-kr(idir3,2);ixOmmax3=ixOmax3-kr(idir3,3);
  ixOpmin1=ixOmin1+kr(idir3,1);ixOpmin2=ixOmin2+kr(idir3,2)
  ixOpmin3=ixOmin3+kr(idir3,3);ixOpmax1=ixOmax1+kr(idir3,1)
  ixOpmax2=ixOmax2+kr(idir3,2);ixOpmax3=ixOmax3+kr(idir3,3);

  if (idir2.le.3) curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     idir1) = curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     idir1) + invdx(idir2)*half*(qvec(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
     jxOmin3:jxOmax3,idir3)-qvec(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
     hxOmin3:hxOmax3,idir3))

  if (idir3.le.3) curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     idir1) = curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     idir1) - invdx(idir3)*half*(qvec(ixOpmin1:ixOpmax1,ixOpmin2:ixOpmax2,&
     ixOpmin3:ixOpmax3,idir2)-qvec(ixOmmin1:ixOmmax1,ixOmmin2:ixOmmax2,&
     ixOmmin3:ixOmmax3,idir2))

  !! When covariant, divide between metric element
  if (covariant) curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     idir1)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir1)&
     /myM%sqrtgamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

end do

end subroutine curl3
!=============================================================================
subroutine curlvector(qvec,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,curlvec,idirmin,idirmin0,&
   ndir0)

! Calculate curl of a vector qvec within ixL

use mod_amrvacdef

integer :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,&
   ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,&
   ixmax3,idir,jdir,kdir,hxOmin1,hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3,&
   jxOmin1,jxOmin2,jxOmin3,jxOmax1,jxOmax2,jxOmax3,ndir0,idirmin0
double precision :: qvec(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:ndir0),&
   curlvec(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,idirmin0:3),&
    invdx(1:ndim)
double precision :: tmp(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
   tmp2(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),surface(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2,ixGlo3:ixGhi3),mydx(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
   ixGlo3:ixGhi3)
!-----------------------------------------------------------------------------

ixmin1=ixOmin1-1;ixmin2=ixOmin2-1;ixmin3=ixOmin3-1;ixmax1=ixOmax1+1
ixmax2=ixOmax2+1;ixmax3=ixOmax3+1;

if (ixImin1>ixmin1.or.ixImax1<ixmax1.or.ixImin2>ixmin2.or.ixImax2<ixmax2&
   .or.ixImin3>ixmin3.or.ixImax3<ixmax3) call mpistop("Error in curl: &
   Non-conforming input limits")

! Calculate curl within ixL: CurlV_i=eps_ijk*d_j V_k
! Curl can have components (idirmin0:3)
! Determine exact value of idirmin while doing the loop.

invdx=1.d0/dxlevel
idirmin=4
curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idirmin0:3)=zero

do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
   if(lvc(idir,jdir,kdir)/=0)then
      tmp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=qvec(ixmin1:ixmax1,&
         ixmin2:ixmax2,ixmin3:ixmax3,kdir)
      hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
      hxOmin3=ixOmin3-kr(jdir,3);hxOmax1=ixOmax1-kr(jdir,1)
      hxOmax2=ixOmax2-kr(jdir,2);hxOmax3=ixOmax3-kr(jdir,3);
      jxOmin1=ixOmin1+kr(jdir,1);jxOmin2=ixOmin2+kr(jdir,2)
      jxOmin3=ixOmin3+kr(jdir,3);jxOmax1=ixOmax1+kr(jdir,1)
      jxOmax2=ixOmax2+kr(jdir,2);jxOmax3=ixOmax3+kr(jdir,3);
      if(slab)then


         tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
            =half*(tmp(jxOmin1:jxOmax1,jxOmin2:jxOmax2,jxOmin3:jxOmax3)-tmp&
            (hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3))*invdx(jdir)

      else
         ! approximate formula, reduces to slab case
         ! and avoids staggering

         if (kdir .le. ndim) then 
            mydx(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)&
               =mygeo%dx(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,kdir)
         else 
            mydx(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=one
         end if

         select case(idir)
           case(1)
             surface(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                =mygeo%surface1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)
             tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                =half*(mydx(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                jxOmin3:jxOmax3)*tmp(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                jxOmin3:jxOmax3) &
                              -mydx(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                                 hxOmin3:hxOmax3)*tmp(hxOmin1:hxOmax1,&
                                 hxOmin2:hxOmax2,hxOmin3:hxOmax3)) &
                     /surface(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                        ixOmin3:ixOmax3) 
           case(2)
             surface(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                =mygeo%surface2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)
             tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                =half*(mydx(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                jxOmin3:jxOmax3)*tmp(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                jxOmin3:jxOmax3) &
                              -mydx(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                                 hxOmin3:hxOmax3)*tmp(hxOmin1:hxOmax1,&
                                 hxOmin2:hxOmax2,hxOmin3:hxOmax3)) &
                     /surface(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                        ixOmin3:ixOmax3) 
           case(3)
             surface(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                =mygeo%surface3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)
             tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                =half*(mydx(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                jxOmin3:jxOmax3)*tmp(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                jxOmin3:jxOmax3) &
                              -mydx(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                                 hxOmin3:hxOmax3)*tmp(hxOmin1:hxOmax1,&
                                 hxOmin2:hxOmax2,hxOmin3:hxOmax3)) &
                     /surface(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                        ixOmin3:ixOmax3) 
          end select
      endif
      if(lvc(idir,jdir,kdir)==1)then
         curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir)&
            =curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            idir)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      else
         curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir)&
            =curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            idir)-tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      endif
      if(idir<idirmin)idirmin=idir
   endif
enddo; enddo; enddo;

end subroutine curlvector 
!=============================================================================
subroutine divvectorS(qvec,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divq)

! Calculate divergence of a vector qvec within ixL
! using limited extrapolation to cell edges

use mod_limiter
use mod_amrvacdef

integer :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,&
   ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision :: qvec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:ndir), divq(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

double precision,dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3):: qL,qR,dqC,ldq,rdq
double precision :: dxdim, invdx(1:ndim)

integer :: hxOmin1,hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3,ixCmin1,ixCmin2,&
   ixCmin3,ixCmax1,ixCmax2,ixCmax3,jxCmin1,jxCmin2,jxCmin3,jxCmax1,jxCmax2,&
   jxCmax3,idims,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,gxCmin1,gxCmin2,&
   gxCmin3,gxCmax1,gxCmax2,gxCmax3,hxCmin1,hxCmin2,hxCmin3,hxCmax1,hxCmax2,&
   hxCmax3,idummy
character*79, save :: savetypelimiter,savetypegradlimiter,save2typelimiter
!-----------------------------------------------------------------------------
ixmin1=ixOmin1-2;ixmin2=ixOmin2-2;ixmin3=ixOmin3-2;ixmax1=ixOmax1+2
ixmax2=ixOmax2+2;ixmax3=ixOmax3+2;

if (ixImin1>ixmin1.or.ixImax1<ixmax1.or.ixImin2>ixmin2.or.ixImax2<ixmax2&
   .or.ixImin3>ixmin3.or.ixImax3<ixmax3) call mpistop("Error in divvectorS: &
   Non-conforming input limits")

idummy=0
invdx=1.d0/dxlevel
divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=zero
do idims=1,ndim
   hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
   hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
   hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);
   ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmin3=hxOmin3;ixCmax1=ixOmax1
   ixCmax2=ixOmax2;ixCmax3=ixOmax3;
   jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
   jxCmin3=ixCmin3+kr(idims,3);jxCmax1=ixCmax1+kr(idims,1)
   jxCmax2=ixCmax2+kr(idims,2);jxCmax3=ixCmax3+kr(idims,3);
   gxCmin1=ixCmin1-kr(idims,1);gxCmin2=ixCmin2-kr(idims,2)
   gxCmin3=ixCmin3-kr(idims,3);gxCmax1=jxCmax1;gxCmax2=jxCmax2
   gxCmax3=jxCmax3;
   hxCmin1=gxCmin1+kr(idims,1);hxCmin2=gxCmin2+kr(idims,2)
   hxCmin3=gxCmin3+kr(idims,3);hxCmax1=gxCmax1+kr(idims,1)
   hxCmax2=gxCmax2+kr(idims,2);hxCmax3=gxCmax3+kr(idims,3);

   qR(gxCmin1:gxCmax1,gxCmin2:gxCmax2,gxCmin3:gxCmax3) = qvec(hxCmin1:hxCmax1,&
      hxCmin2:hxCmax2,hxCmin3:hxCmax3,idims)
   qL(gxCmin1:gxCmax1,gxCmin2:gxCmax2,gxCmin3:gxCmax3) = qvec(gxCmin1:gxCmax1,&
      gxCmin2:gxCmax2,gxCmin3:gxCmax3,idims)

  !==================================================
  ! reconstruct the variable itself, from left and right
  !==================================================
  select case (typegradlimiter)
     
  case (limiter_ppm)
     
     dqC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=qvec&
        (ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idims)
     call PPMlimitervar(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idims,dqC,dqC,qL,qR)

  case default

     dqC(gxCmin1:gxCmax1,gxCmin2:gxCmax2,gxCmin3:gxCmax3)= qR(gxCmin1:gxCmax1,&
        gxCmin2:gxCmax2,gxCmin3:gxCmax3)-qL(gxCmin1:gxCmax1,gxCmin2:gxCmax2,&
        gxCmin3:gxCmax3)

     call dwlimiter2(dqC,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        gxCmin1,gxCmin2,gxCmin3,gxCmax1,gxCmax2,gxCmax3,idims,typegradlimiter,&
        ldq,rdq)

     qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = qL(ixCmin1:ixCmax1,&
        ixCmin2:ixCmax2,ixCmin3:ixCmax3) + half*ldq(ixCmin1:ixCmax1,&
        ixCmin2:ixCmax2,ixCmin3:ixCmax3)
     qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = qR(ixCmin1:ixCmax1,&
        ixCmin2:ixCmax2,ixCmin3:ixCmax3) - half*rdq(jxCmin1:jxCmax1,&
        jxCmin2:jxCmax2,jxCmin3:jxCmax3)

  end select
  !==================================================

   if (slab) then
     divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=divq&
        (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)+half*(qR&
        (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)-qL(hxOmin1:hxOmax1,&
        hxOmin2:hxOmax2,hxOmin3:hxOmax3))*invdx(idims)
   else
     select case(idims)
     case(1)
        qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
           =mygeo%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)
        qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
           =mygeo%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)
        divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
           =divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)+(qR&
           (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)-qL&
           (hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3))&
           /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) 
     case(2)
        qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
           =mygeo%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)
        qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
           =mygeo%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)
        divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
           =divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)+(qR&
           (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)-qL&
           (hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3))&
           /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) 
     case(3)
        qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
           =mygeo%surfaceC3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)
        qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)&
           =mygeo%surfaceC3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)
        divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
           =divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)+(qR&
           (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)-qL&
           (hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3))&
           /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) 
      end select
   end if
end do

end subroutine divvectorS
!=============================================================================
subroutine rec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,&
   ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idir,q,qL,qR)

  ! Reconstruct scalar q within ixO^L to 1/2 dx in direction idir
  ! Return both left and right reconstructed values 

use mod_limiter
use mod_amrvacdef

integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3, ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, idir
double precision, intent(in)       :: q(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)
double precision, intent(out)      :: qL(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3), qR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

double precision                   :: qC(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)
double precision,dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)  :: dqC,ldq,rdq
integer                            :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
   ixOmax3,jxCmin1,jxCmin2,jxCmin3,jxCmax1,jxCmax2,jxCmax3,gxCmin1,gxCmin2,&
   gxCmin3,gxCmax1,gxCmax2,gxCmax3,hxCmin1,hxCmin2,hxCmin3,hxCmax1,hxCmax2,&
   hxCmax3
!-----------------------------------------------------------------------------

jxCmin1=ixCmin1+kr(idir,1);jxCmin2=ixCmin2+kr(idir,2)
jxCmin3=ixCmin3+kr(idir,3);jxCmax1=ixCmax1+kr(idir,1)
jxCmax2=ixCmax2+kr(idir,2);jxCmax3=ixCmax3+kr(idir,3);
gxCmin1=ixCmin1-kr(idir,1);gxCmin2=ixCmin2-kr(idir,2)
gxCmin3=ixCmin3-kr(idir,3);gxCmax1=jxCmax1;gxCmax2=jxCmax2;gxCmax3=jxCmax3;
hxCmin1=gxCmin1+kr(idir,1);hxCmin2=gxCmin2+kr(idir,2)
hxCmin3=gxCmin3+kr(idir,3);hxCmax1=gxCmax1+kr(idir,1)
hxCmax2=gxCmax2+kr(idir,2);hxCmax3=gxCmax3+kr(idir,3);

qR(gxCmin1:gxCmax1,gxCmin2:gxCmax2,gxCmin3:gxCmax3) = q(hxCmin1:hxCmax1,&
   hxCmin2:hxCmax2,hxCmin3:hxCmax3)
qL(gxCmin1:gxCmax1,gxCmin2:gxCmax2,gxCmin3:gxCmax3) = q(gxCmin1:gxCmax1,&
   gxCmin2:gxCmax2,gxCmin3:gxCmax3)

select case (typelimiter)
   
case (limiter_ppm)
   ! the ordinary grid-index:
   ixOmin1=ixCmin1+kr(idir,1);ixOmin2=ixCmin2+kr(idir,2)
   ixOmin3=ixCmin3+kr(idir,3);
   ixOmax1=ixCmax1;ixOmax2=ixCmax2;ixOmax3=ixCmax3;
   call PPMlimitervar(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
      ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idir,q,q,qL,qR)
   
case (limiter_mp5)
   call MP5limitervar(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,&
      ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idir,q,qL,qR)
   
case (limiter_weno5)
   call WENO5limitervar(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idir,q,qL,qR)
   
case (limiter_wenoZP)
   call WENOZPlimitervar(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idir,dxlevel(idir),q,qL,&
      qR)
   
case default

   dqC(gxCmin1:gxCmax1,gxCmin2:gxCmax2,gxCmin3:gxCmax3)= qR(gxCmin1:gxCmax1,&
      gxCmin2:gxCmax2,gxCmin3:gxCmax3)-qL(gxCmin1:gxCmax1,gxCmin2:gxCmax2,&
      gxCmin3:gxCmax3)
   
   
   call dwlimiter2(dqC,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      gxCmin1,gxCmin2,gxCmin3,gxCmax1,gxCmax2,gxCmax3,idir,typelimiter,ldq,&
      rdq)
   
   qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = qL(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3) + half*ldq(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3)
   qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = qR(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,ixCmin3:ixCmax3) - half*rdq(jxCmin1:jxCmax1,&
      jxCmin2:jxCmax2,jxCmin3:jxCmax3)
   
end select


end subroutine rec
!=============================================================================
