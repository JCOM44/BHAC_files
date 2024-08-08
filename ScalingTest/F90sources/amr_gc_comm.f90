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

subroutine init_comm_gc(timein,psuse,psuseCo)
!! This routine initializes the ghost cell communications:
!! calculates buffer sizes, allocates buffers
!! and calls the routines 'make_task_' that place
!! all the receive requests.

use mod_comm_gc
use mod_amrvacdef

double precision :: timein
type(state), dimension(ngridshi), target    :: psuse, psuseCo

! ... local ...
! Counters and auxiliaries
integer :: i1,i2,i3, ic1,ic2,ic3, inc1,inc2,inc3, idir, igrid, iigrid,&
    ineighbor, ipe_neighbor, my_neighbor_type
!-------------------------------------------------------------------------------

!! Set pointers to the state structures
time=timein
pstate => psuse
pstateCo => psuseCo

!! Define block extents. ixM and ixG from amrvacdef.f ----
ixCoGmin1=1;ixCoGmin2=1;ixCoGmin3=1;
ixCoGmax1=ixGhi1/2+dixB;ixCoGmax2=ixGhi2/2+dixB;ixCoGmax3=ixGhi3/2+dixB;
ixCoMmin1=ixCoGmin1+dixB;ixCoMmin2=ixCoGmin2+dixB;ixCoMmin3=ixCoGmin3+dixB
ixCoMmax1=ixCoGmax1-dixB;ixCoMmax2=ixCoGmax2-dixB;ixCoMmax3=ixCoGmax3-dixB;

nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;nx3=ixMhi3-ixMlo3+1;
nxCo1=nx1/2;nxCo2=nx2/2;nxCo3=nx3/2;

select case (typeghostfill)
case ("copy")
   interpolation_order=1
case ("linear","unlimit")
   interpolation_order=2
case default
   write (unitterm,*) "Undefined typeghostfill ",typeghostfill
   call mpistop("")
end select
dixBCo=int((dixB+1)/2)

if (dixBCo+interpolation_order-1>dixB) then
   call mpistop("interpolation order for prolongation in getbc too high")
end if

!! -------------- Limits for communications -------------------
!! ... same resolution level ...
!! Cell-centred variables


ixS_srl_min1(-1)=ixMlo1
ixS_srl_max1(-1)=ixMlo1-1+dixB
ixS_srl_min1(0) =ixMlo1
ixS_srl_max1(0) =ixMhi1
ixS_srl_min1(1) =ixMhi1+1-dixB
ixS_srl_max1(1) =ixMhi1

ixR_srl_min1(-1)=1
ixR_srl_max1(-1)=dixB
ixR_srl_min1(0) =ixMlo1
ixR_srl_max1(0) =ixMhi1
ixR_srl_min1(1) =ixMhi1+1
ixR_srl_max1(1) =ixGhi1




ixS_srl_min2(-1)=ixMlo2
ixS_srl_max2(-1)=ixMlo2-1+dixB
ixS_srl_min2(0) =ixMlo2
ixS_srl_max2(0) =ixMhi2
ixS_srl_min2(1) =ixMhi2+1-dixB
ixS_srl_max2(1) =ixMhi2

ixR_srl_min2(-1)=1
ixR_srl_max2(-1)=dixB
ixR_srl_min2(0) =ixMlo2
ixR_srl_max2(0) =ixMhi2
ixR_srl_min2(1) =ixMhi2+1
ixR_srl_max2(1) =ixGhi2




ixS_srl_min3(-1)=ixMlo3
ixS_srl_max3(-1)=ixMlo3-1+dixB
ixS_srl_min3(0) =ixMlo3
ixS_srl_max3(0) =ixMhi3
ixS_srl_min3(1) =ixMhi3+1-dixB
ixS_srl_max3(1) =ixMhi3

ixR_srl_min3(-1)=1
ixR_srl_max3(-1)=dixB
ixR_srl_min3(0) =ixMlo3
ixR_srl_max3(0) =ixMhi3
ixR_srl_min3(1) =ixMhi3+1
ixR_srl_max3(1) =ixGhi3



!! Staggered (face-allocated) variables
do idir=1,3


   ixS_srl_stg_min1(idir,-1)=ixMlo1-kr(idir,1)
   ixS_srl_stg_max1(idir,-1)=ixMlo1-1+dixB
   ixS_srl_stg_min1(idir,0) =ixMlo1-kr(idir,1)
   ixS_srl_stg_max1(idir,0) =ixMhi1
   ixS_srl_stg_min1(idir,1) =ixMhi1-dixB+1-kr(idir,1)
   ixS_srl_stg_max1(idir,1) =ixMhi1
   
   ixR_srl_stg_min1(idir,-1)=1-kr(idir,1)
   ixR_srl_stg_max1(idir,-1)=dixB
   ixR_srl_stg_min1(idir,0) =ixMlo1-kr(idir,1)
   ixR_srl_stg_max1(idir,0) =ixMhi1
   ixR_srl_stg_min1(idir,1) =ixMhi1+1-kr(idir,1)
   ixR_srl_stg_max1(idir,1) =ixGhi1




   ixS_srl_stg_min2(idir,-1)=ixMlo2-kr(idir,2)
   ixS_srl_stg_max2(idir,-1)=ixMlo2-1+dixB
   ixS_srl_stg_min2(idir,0) =ixMlo2-kr(idir,2)
   ixS_srl_stg_max2(idir,0) =ixMhi2
   ixS_srl_stg_min2(idir,1) =ixMhi2-dixB+1-kr(idir,2)
   ixS_srl_stg_max2(idir,1) =ixMhi2
   
   ixR_srl_stg_min2(idir,-1)=1-kr(idir,2)
   ixR_srl_stg_max2(idir,-1)=dixB
   ixR_srl_stg_min2(idir,0) =ixMlo2-kr(idir,2)
   ixR_srl_stg_max2(idir,0) =ixMhi2
   ixR_srl_stg_min2(idir,1) =ixMhi2+1-kr(idir,2)
   ixR_srl_stg_max2(idir,1) =ixGhi2




   ixS_srl_stg_min3(idir,-1)=ixMlo3-kr(idir,3)
   ixS_srl_stg_max3(idir,-1)=ixMlo3-1+dixB
   ixS_srl_stg_min3(idir,0) =ixMlo3-kr(idir,3)
   ixS_srl_stg_max3(idir,0) =ixMhi3
   ixS_srl_stg_min3(idir,1) =ixMhi3-dixB+1-kr(idir,3)
   ixS_srl_stg_max3(idir,1) =ixMhi3
   
   ixR_srl_stg_min3(idir,-1)=1-kr(idir,3)
   ixR_srl_stg_max3(idir,-1)=dixB
   ixR_srl_stg_min3(idir,0) =ixMlo3-kr(idir,3)
   ixR_srl_stg_max3(idir,0) =ixMhi3
   ixR_srl_stg_min3(idir,1) =ixMhi3+1-kr(idir,3)
   ixR_srl_stg_max3(idir,1) =ixGhi3


end do


!! Sizes for srl communications
do i3=-1,1
do i2=-1,1
do i1=-1,1

   ! Cell-centred variables
   sizes_srl_send(i1,i2,i3)=(nwflux+nwaux)*(ixS_srl_max1(i1)-ixS_srl_min1&
      (i1)+1)*(ixS_srl_max2(i2)-ixS_srl_min2(i2)+1)*(ixS_srl_max3&
      (i3)-ixS_srl_min3(i3)+1)
   sizes_srl_recv(i1,i2,i3)=(nwflux+nwaux)*(ixR_srl_max1(i1)-ixR_srl_min1&
      (i1)+1)*(ixR_srl_max2(i2)-ixR_srl_min2(i2)+1)*(ixR_srl_max3&
      (i3)-ixR_srl_min3(i3)+1)
   sizes_srl_send_total(i1,i2,i3)=sizes_srl_send(i1,i2,i3)
   sizes_srl_recv_total(i1,i2,i3)=sizes_srl_recv(i1,i2,i3)



   ! Staggered (face-allocated) variables   
   do idir=1,3
     sizes_srl_send_stg(idir,i1,i2,i3)=(ixS_srl_stg_max1(idir,&
        i1)-ixS_srl_stg_min1(idir,i1)+1)*(ixS_srl_stg_max2(idir,&
        i2)-ixS_srl_stg_min2(idir,i2)+1)*(ixS_srl_stg_max3(idir,&
        i3)-ixS_srl_stg_min3(idir,i3)+1)
     sizes_srl_recv_stg(idir,i1,i2,i3)=(ixR_srl_stg_max1(idir,&
        i1)-ixR_srl_stg_min1(idir,i1)+1)*(ixR_srl_stg_max2(idir,&
        i2)-ixR_srl_stg_min2(idir,i2)+1)*(ixR_srl_stg_max3(idir,&
        i3)-ixR_srl_stg_min3(idir,i3)+1)

     sizes_srl_send_total(i1,i2,i3)=sizes_srl_send_total(i1,i2,&
        i3)+sizes_srl_send_stg(idir,i1,i2,i3)
     sizes_srl_recv_total(i1,i2,i3)=sizes_srl_recv_total(i1,i2,&
        i3)+sizes_srl_recv_stg(idir,i1,i2,i3)

   end do


end do
end do
end do

if (levmin/=levmax) then
!! ... restriction ...

   ixS_r_min1(-1)=ixCoMmin1
   ixS_r_min1(0) =ixCoMmin1
   ixS_r_min1(1) =ixCoMmax1+1-dixB
   ixS_r_max1(-1)=ixCoMmin1-1+dixB
   ixS_r_max1(0) =ixCoMmax1
   ixS_r_max1(1) =ixCoMmax1

   ixR_r_min1(0)=1
   ixR_r_min1(1)=ixMlo1
   ixR_r_min1(2)=ixMlo1+nxCo1
   ixR_r_min1(3)=ixMhi1+1
   ixR_r_max1(0)=dixB
   ixR_r_max1(1)=ixMlo1-1+nxCo1
   ixR_r_max1(2)=ixMhi1
   ixR_r_max1(3)=ixGhi1


   ixS_r_min2(-1)=ixCoMmin2
   ixS_r_min2(0) =ixCoMmin2
   ixS_r_min2(1) =ixCoMmax2+1-dixB
   ixS_r_max2(-1)=ixCoMmin2-1+dixB
   ixS_r_max2(0) =ixCoMmax2
   ixS_r_max2(1) =ixCoMmax2

   ixR_r_min2(0)=1
   ixR_r_min2(1)=ixMlo2
   ixR_r_min2(2)=ixMlo2+nxCo2
   ixR_r_min2(3)=ixMhi2+1
   ixR_r_max2(0)=dixB
   ixR_r_max2(1)=ixMlo2-1+nxCo2
   ixR_r_max2(2)=ixMhi2
   ixR_r_max2(3)=ixGhi2


   ixS_r_min3(-1)=ixCoMmin3
   ixS_r_min3(0) =ixCoMmin3
   ixS_r_min3(1) =ixCoMmax3+1-dixB
   ixS_r_max3(-1)=ixCoMmin3-1+dixB
   ixS_r_max3(0) =ixCoMmax3
   ixS_r_max3(1) =ixCoMmax3

   ixR_r_min3(0)=1
   ixR_r_min3(1)=ixMlo3
   ixR_r_min3(2)=ixMlo3+nxCo3
   ixR_r_min3(3)=ixMhi3+1
   ixR_r_max3(0)=dixB
   ixR_r_max3(1)=ixMlo3-1+nxCo3
   ixR_r_max3(2)=ixMhi3
   ixR_r_max3(3)=ixGhi3



!! Staggered (face-allocated) variables
do idir=1,3

   ixS_r_stg_min1(idir,-1)=ixCoMmin1-kr(idir,1)
   ixS_r_stg_max1(idir,-1)=ixCoMmin1-1+dixB
   ixS_r_stg_min1(idir,0) =ixCoMmin1-kr(idir,1)
   ixS_r_stg_max1(idir,0) =ixCoMmax1
   ixS_r_stg_min1(idir,1) =ixCoMmax1+1-dixB-kr(idir,1)
   ixS_r_stg_max1(idir,1) =ixCoMmax1
 
   ixR_r_stg_min1(idir,0)=1-kr(idir,1)
   ixR_r_stg_max1(idir,0)=dixB
   ixR_r_stg_min1(idir,1)=ixMlo1-kr(idir,1)
   ixR_r_stg_max1(idir,1)=ixMlo1-1+nxCo1
   ixR_r_stg_min1(idir,2)=ixMlo1+nxCo1-kr(idir,1)
   ixR_r_stg_max1(idir,2)=ixMhi1
   ixR_r_stg_min1(idir,3)=ixMhi1+1-kr(idir,1)
   ixR_r_stg_max1(idir,3)=ixGhi1



   ixS_r_stg_min2(idir,-1)=ixCoMmin2-kr(idir,2)
   ixS_r_stg_max2(idir,-1)=ixCoMmin2-1+dixB
   ixS_r_stg_min2(idir,0) =ixCoMmin2-kr(idir,2)
   ixS_r_stg_max2(idir,0) =ixCoMmax2
   ixS_r_stg_min2(idir,1) =ixCoMmax2+1-dixB-kr(idir,2)
   ixS_r_stg_max2(idir,1) =ixCoMmax2
 
   ixR_r_stg_min2(idir,0)=1-kr(idir,2)
   ixR_r_stg_max2(idir,0)=dixB
   ixR_r_stg_min2(idir,1)=ixMlo2-kr(idir,2)
   ixR_r_stg_max2(idir,1)=ixMlo2-1+nxCo2
   ixR_r_stg_min2(idir,2)=ixMlo2+nxCo2-kr(idir,2)
   ixR_r_stg_max2(idir,2)=ixMhi2
   ixR_r_stg_min2(idir,3)=ixMhi2+1-kr(idir,2)
   ixR_r_stg_max2(idir,3)=ixGhi2



   ixS_r_stg_min3(idir,-1)=ixCoMmin3-kr(idir,3)
   ixS_r_stg_max3(idir,-1)=ixCoMmin3-1+dixB
   ixS_r_stg_min3(idir,0) =ixCoMmin3-kr(idir,3)
   ixS_r_stg_max3(idir,0) =ixCoMmax3
   ixS_r_stg_min3(idir,1) =ixCoMmax3+1-dixB-kr(idir,3)
   ixS_r_stg_max3(idir,1) =ixCoMmax3
 
   ixR_r_stg_min3(idir,0)=1-kr(idir,3)
   ixR_r_stg_max3(idir,0)=dixB
   ixR_r_stg_min3(idir,1)=ixMlo3-kr(idir,3)
   ixR_r_stg_max3(idir,1)=ixMlo3-1+nxCo3
   ixR_r_stg_min3(idir,2)=ixMlo3+nxCo3-kr(idir,3)
   ixR_r_stg_max3(idir,2)=ixMhi3
   ixR_r_stg_min3(idir,3)=ixMhi3+1-kr(idir,3)
   ixR_r_stg_max3(idir,3)=ixGhi3


end do


!! ... prolongation

   ixS_p_min1(0)=ixMlo1
   ixS_p_max1(0)=ixMlo1-1+dixBCo+(interpolation_order-1)
   ixS_p_min1(1)=ixMlo1
   ixS_p_max1(1)=ixMlo1-1+nxCo1+dixBCo+(interpolation_order-1)
   ixS_p_min1(2)=ixMlo1+nxCo1-dixBCo-(interpolation_order-1)
   ixS_p_max1(2)=ixMhi1
   ixS_p_min1(3)=ixMhi1+1-dixBCo-(interpolation_order-1)
   ixS_p_max1(3)=ixMhi1

   ixR_p_min1(0)=ixCoMmin1-dixBCo-(interpolation_order-1)
   ixR_p_max1(0)=dixB
   ixR_p_min1(1)=ixCoMmin1
   ixR_p_max1(1)=ixCoMmax1+dixBCo+(interpolation_order-1)
   ixR_p_min1(2)=ixCoMmin1-dixBCo-(interpolation_order-1)
   ixR_p_max1(2)=ixCoMmax1
   ixR_p_min1(3)=ixCoMmax1+1
   ixR_p_max1(3)=ixCoMmax1+dixBCo+(interpolation_order-1)


   ixS_p_min2(0)=ixMlo2
   ixS_p_max2(0)=ixMlo2-1+dixBCo+(interpolation_order-1)
   ixS_p_min2(1)=ixMlo2
   ixS_p_max2(1)=ixMlo2-1+nxCo2+dixBCo+(interpolation_order-1)
   ixS_p_min2(2)=ixMlo2+nxCo2-dixBCo-(interpolation_order-1)
   ixS_p_max2(2)=ixMhi2
   ixS_p_min2(3)=ixMhi2+1-dixBCo-(interpolation_order-1)
   ixS_p_max2(3)=ixMhi2

   ixR_p_min2(0)=ixCoMmin2-dixBCo-(interpolation_order-1)
   ixR_p_max2(0)=dixB
   ixR_p_min2(1)=ixCoMmin2
   ixR_p_max2(1)=ixCoMmax2+dixBCo+(interpolation_order-1)
   ixR_p_min2(2)=ixCoMmin2-dixBCo-(interpolation_order-1)
   ixR_p_max2(2)=ixCoMmax2
   ixR_p_min2(3)=ixCoMmax2+1
   ixR_p_max2(3)=ixCoMmax2+dixBCo+(interpolation_order-1)


   ixS_p_min3(0)=ixMlo3
   ixS_p_max3(0)=ixMlo3-1+dixBCo+(interpolation_order-1)
   ixS_p_min3(1)=ixMlo3
   ixS_p_max3(1)=ixMlo3-1+nxCo3+dixBCo+(interpolation_order-1)
   ixS_p_min3(2)=ixMlo3+nxCo3-dixBCo-(interpolation_order-1)
   ixS_p_max3(2)=ixMhi3
   ixS_p_min3(3)=ixMhi3+1-dixBCo-(interpolation_order-1)
   ixS_p_max3(3)=ixMhi3

   ixR_p_min3(0)=ixCoMmin3-dixBCo-(interpolation_order-1)
   ixR_p_max3(0)=dixB
   ixR_p_min3(1)=ixCoMmin3
   ixR_p_max3(1)=ixCoMmax3+dixBCo+(interpolation_order-1)
   ixR_p_min3(2)=ixCoMmin3-dixBCo-(interpolation_order-1)
   ixR_p_max3(2)=ixCoMmax3
   ixR_p_min3(3)=ixCoMmax3+1
   ixR_p_max3(3)=ixCoMmax3+dixBCo+(interpolation_order-1)



   do idir = 1, 3
   if (idir.eq.1) then
     ! Parallel components

     ixS_p_stg_min1(idir,0)=ixMlo1-1 ! -1 to make redundant 
     ixS_p_stg_max1(idir,0)=ixMlo1-1+dixBCo
     ixS_p_stg_min1(idir,1)=ixMlo1-1 ! -1 to make redundant 
     ixS_p_stg_max1(idir,1)=ixMlo1-1+nxCo1+dixBCo
     ixS_p_stg_min1(idir,2)=ixMhi1-nxCo1-dixBCo
     ixS_p_stg_max1(idir,2)=ixMhi1
     ixS_p_stg_min1(idir,3)=ixMhi1-dixBCo
     ixS_p_stg_max1(idir,3)=ixMhi1

     ixR_p_stg_min1(idir,0)=ixCoMmin1-1-dixBCo
     ixR_p_stg_max1(idir,0)=ixCoMmin1-1
     ixR_p_stg_min1(idir,1)=ixCoMmin1-1 ! -1 to make redundant 
     ixR_p_stg_max1(idir,1)=ixCoMmax1+dixBCo
     ixR_p_stg_min1(idir,2)=ixCoMmin1-1-dixBCo
     ixR_p_stg_max1(idir,2)=ixCoMmax1
     ixR_p_stg_min1(idir,3)=ixCoMmax1+1-1 ! -1 to make redundant 
     ixR_p_stg_max1(idir,3)=ixCoMmax1+dixBCo


     ixS_p_stg_min1(idir,0)=ixMlo1-1 ! -1 to make redundant 
     ixS_p_stg_max1(idir,0)=ixMlo1-1+dixBCo
     ixS_p_stg_min1(idir,1)=ixMlo1-1 ! -1 to make redundant 
     ixS_p_stg_max1(idir,1)=ixMlo1-1+nxCo1+dixBCo
     ixS_p_stg_min1(idir,2)=ixMhi1-nxCo1-dixBCo
     ixS_p_stg_max1(idir,2)=ixMhi1
     ixS_p_stg_min1(idir,3)=ixMhi1-dixBCo
     ixS_p_stg_max1(idir,3)=ixMhi1

     ixR_p_stg_min1(idir,0)=ixCoMmin1-1-dixBCo
     ixR_p_stg_max1(idir,0)=ixCoMmin1-1
     ixR_p_stg_min1(idir,1)=ixCoMmin1-1 ! -1 to make redundant 
     ixR_p_stg_max1(idir,1)=ixCoMmax1+dixBCo
     ixR_p_stg_min1(idir,2)=ixCoMmin1-1-dixBCo
     ixR_p_stg_max1(idir,2)=ixCoMmax1
     ixR_p_stg_min1(idir,3)=ixCoMmax1+1-1 ! -1 to make redundant 
     ixR_p_stg_max1(idir,3)=ixCoMmax1+dixBCo


     ixS_p_stg_min1(idir,0)=ixMlo1-1 ! -1 to make redundant 
     ixS_p_stg_max1(idir,0)=ixMlo1-1+dixBCo
     ixS_p_stg_min1(idir,1)=ixMlo1-1 ! -1 to make redundant 
     ixS_p_stg_max1(idir,1)=ixMlo1-1+nxCo1+dixBCo
     ixS_p_stg_min1(idir,2)=ixMhi1-nxCo1-dixBCo
     ixS_p_stg_max1(idir,2)=ixMhi1
     ixS_p_stg_min1(idir,3)=ixMhi1-dixBCo
     ixS_p_stg_max1(idir,3)=ixMhi1

     ixR_p_stg_min1(idir,0)=ixCoMmin1-1-dixBCo
     ixR_p_stg_max1(idir,0)=ixCoMmin1-1
     ixR_p_stg_min1(idir,1)=ixCoMmin1-1 ! -1 to make redundant 
     ixR_p_stg_max1(idir,1)=ixCoMmax1+dixBCo
     ixR_p_stg_min1(idir,2)=ixCoMmin1-1-dixBCo
     ixR_p_stg_max1(idir,2)=ixCoMmax1
     ixR_p_stg_min1(idir,3)=ixCoMmax1+1-1 ! -1 to make redundant 
     ixR_p_stg_max1(idir,3)=ixCoMmax1+dixBCo

    else

     ! Perpendicular component
     ixS_p_stg_min1(idir,0)=ixMlo1
     ixS_p_stg_max1(idir,0)=ixMlo1-1+dixBCo+(interpolation_order-1)
     ixS_p_stg_min1(idir,1)=ixMlo1
     ixS_p_stg_max1(idir,1)=ixMlo1-1+nxCo1+dixBCo+(interpolation_order-1)
     ixS_p_stg_min1(idir,2)=ixMhi1+1-nxCo1-dixBCo-(interpolation_order-1)
     ixS_p_stg_max1(idir,2)=ixMhi1
     ixS_p_stg_min1(idir,3)=ixMhi1+1-dixBCo-(interpolation_order-1)
     ixS_p_stg_max1(idir,3)=ixMhi1
 
     ixR_p_stg_min1(idir,0)=ixCoMmin1-dixBCo-(interpolation_order-1)
     ixR_p_stg_max1(idir,0)=ixCoMmin1-1
     ixR_p_stg_min1(idir,1)=ixCoMmin1
     ixR_p_stg_max1(idir,1)=ixCoMmax1+dixBCo+(interpolation_order-1)
     ixR_p_stg_min1(idir,2)=ixCoMmin1-dixBCo-(interpolation_order-1)
     ixR_p_stg_max1(idir,2)=ixCoMmax1
     ixR_p_stg_min1(idir,3)=ixCoMmax1+1
     ixR_p_stg_max1(idir,3)=ixCoMmax1+dixBCo+(interpolation_order-1)


     ! Perpendicular component
     ixS_p_stg_min1(idir,0)=ixMlo1
     ixS_p_stg_max1(idir,0)=ixMlo1-1+dixBCo+(interpolation_order-1)
     ixS_p_stg_min1(idir,1)=ixMlo1
     ixS_p_stg_max1(idir,1)=ixMlo1-1+nxCo1+dixBCo+(interpolation_order-1)
     ixS_p_stg_min1(idir,2)=ixMhi1+1-nxCo1-dixBCo-(interpolation_order-1)
     ixS_p_stg_max1(idir,2)=ixMhi1
     ixS_p_stg_min1(idir,3)=ixMhi1+1-dixBCo-(interpolation_order-1)
     ixS_p_stg_max1(idir,3)=ixMhi1
 
     ixR_p_stg_min1(idir,0)=ixCoMmin1-dixBCo-(interpolation_order-1)
     ixR_p_stg_max1(idir,0)=ixCoMmin1-1
     ixR_p_stg_min1(idir,1)=ixCoMmin1
     ixR_p_stg_max1(idir,1)=ixCoMmax1+dixBCo+(interpolation_order-1)
     ixR_p_stg_min1(idir,2)=ixCoMmin1-dixBCo-(interpolation_order-1)
     ixR_p_stg_max1(idir,2)=ixCoMmax1
     ixR_p_stg_min1(idir,3)=ixCoMmax1+1
     ixR_p_stg_max1(idir,3)=ixCoMmax1+dixBCo+(interpolation_order-1)


     ! Perpendicular component
     ixS_p_stg_min1(idir,0)=ixMlo1
     ixS_p_stg_max1(idir,0)=ixMlo1-1+dixBCo+(interpolation_order-1)
     ixS_p_stg_min1(idir,1)=ixMlo1
     ixS_p_stg_max1(idir,1)=ixMlo1-1+nxCo1+dixBCo+(interpolation_order-1)
     ixS_p_stg_min1(idir,2)=ixMhi1+1-nxCo1-dixBCo-(interpolation_order-1)
     ixS_p_stg_max1(idir,2)=ixMhi1
     ixS_p_stg_min1(idir,3)=ixMhi1+1-dixBCo-(interpolation_order-1)
     ixS_p_stg_max1(idir,3)=ixMhi1
 
     ixR_p_stg_min1(idir,0)=ixCoMmin1-dixBCo-(interpolation_order-1)
     ixR_p_stg_max1(idir,0)=ixCoMmin1-1
     ixR_p_stg_min1(idir,1)=ixCoMmin1
     ixR_p_stg_max1(idir,1)=ixCoMmax1+dixBCo+(interpolation_order-1)
     ixR_p_stg_min1(idir,2)=ixCoMmin1-dixBCo-(interpolation_order-1)
     ixR_p_stg_max1(idir,2)=ixCoMmax1
     ixR_p_stg_min1(idir,3)=ixCoMmax1+1
     ixR_p_stg_max1(idir,3)=ixCoMmax1+dixBCo+(interpolation_order-1)

    end if
   if (idir.eq.2) then
     ! Parallel components

     ixS_p_stg_min2(idir,0)=ixMlo2-1 ! -1 to make redundant 
     ixS_p_stg_max2(idir,0)=ixMlo2-1+dixBCo
     ixS_p_stg_min2(idir,1)=ixMlo2-1 ! -1 to make redundant 
     ixS_p_stg_max2(idir,1)=ixMlo2-1+nxCo2+dixBCo
     ixS_p_stg_min2(idir,2)=ixMhi2-nxCo2-dixBCo
     ixS_p_stg_max2(idir,2)=ixMhi2
     ixS_p_stg_min2(idir,3)=ixMhi2-dixBCo
     ixS_p_stg_max2(idir,3)=ixMhi2

     ixR_p_stg_min2(idir,0)=ixCoMmin2-1-dixBCo
     ixR_p_stg_max2(idir,0)=ixCoMmin2-1
     ixR_p_stg_min2(idir,1)=ixCoMmin2-1 ! -1 to make redundant 
     ixR_p_stg_max2(idir,1)=ixCoMmax2+dixBCo
     ixR_p_stg_min2(idir,2)=ixCoMmin2-1-dixBCo
     ixR_p_stg_max2(idir,2)=ixCoMmax2
     ixR_p_stg_min2(idir,3)=ixCoMmax2+1-1 ! -1 to make redundant 
     ixR_p_stg_max2(idir,3)=ixCoMmax2+dixBCo


     ixS_p_stg_min2(idir,0)=ixMlo2-1 ! -1 to make redundant 
     ixS_p_stg_max2(idir,0)=ixMlo2-1+dixBCo
     ixS_p_stg_min2(idir,1)=ixMlo2-1 ! -1 to make redundant 
     ixS_p_stg_max2(idir,1)=ixMlo2-1+nxCo2+dixBCo
     ixS_p_stg_min2(idir,2)=ixMhi2-nxCo2-dixBCo
     ixS_p_stg_max2(idir,2)=ixMhi2
     ixS_p_stg_min2(idir,3)=ixMhi2-dixBCo
     ixS_p_stg_max2(idir,3)=ixMhi2

     ixR_p_stg_min2(idir,0)=ixCoMmin2-1-dixBCo
     ixR_p_stg_max2(idir,0)=ixCoMmin2-1
     ixR_p_stg_min2(idir,1)=ixCoMmin2-1 ! -1 to make redundant 
     ixR_p_stg_max2(idir,1)=ixCoMmax2+dixBCo
     ixR_p_stg_min2(idir,2)=ixCoMmin2-1-dixBCo
     ixR_p_stg_max2(idir,2)=ixCoMmax2
     ixR_p_stg_min2(idir,3)=ixCoMmax2+1-1 ! -1 to make redundant 
     ixR_p_stg_max2(idir,3)=ixCoMmax2+dixBCo


     ixS_p_stg_min2(idir,0)=ixMlo2-1 ! -1 to make redundant 
     ixS_p_stg_max2(idir,0)=ixMlo2-1+dixBCo
     ixS_p_stg_min2(idir,1)=ixMlo2-1 ! -1 to make redundant 
     ixS_p_stg_max2(idir,1)=ixMlo2-1+nxCo2+dixBCo
     ixS_p_stg_min2(idir,2)=ixMhi2-nxCo2-dixBCo
     ixS_p_stg_max2(idir,2)=ixMhi2
     ixS_p_stg_min2(idir,3)=ixMhi2-dixBCo
     ixS_p_stg_max2(idir,3)=ixMhi2

     ixR_p_stg_min2(idir,0)=ixCoMmin2-1-dixBCo
     ixR_p_stg_max2(idir,0)=ixCoMmin2-1
     ixR_p_stg_min2(idir,1)=ixCoMmin2-1 ! -1 to make redundant 
     ixR_p_stg_max2(idir,1)=ixCoMmax2+dixBCo
     ixR_p_stg_min2(idir,2)=ixCoMmin2-1-dixBCo
     ixR_p_stg_max2(idir,2)=ixCoMmax2
     ixR_p_stg_min2(idir,3)=ixCoMmax2+1-1 ! -1 to make redundant 
     ixR_p_stg_max2(idir,3)=ixCoMmax2+dixBCo

    else

     ! Perpendicular component
     ixS_p_stg_min2(idir,0)=ixMlo2
     ixS_p_stg_max2(idir,0)=ixMlo2-1+dixBCo+(interpolation_order-1)
     ixS_p_stg_min2(idir,1)=ixMlo2
     ixS_p_stg_max2(idir,1)=ixMlo2-1+nxCo2+dixBCo+(interpolation_order-1)
     ixS_p_stg_min2(idir,2)=ixMhi2+1-nxCo2-dixBCo-(interpolation_order-1)
     ixS_p_stg_max2(idir,2)=ixMhi2
     ixS_p_stg_min2(idir,3)=ixMhi2+1-dixBCo-(interpolation_order-1)
     ixS_p_stg_max2(idir,3)=ixMhi2
 
     ixR_p_stg_min2(idir,0)=ixCoMmin2-dixBCo-(interpolation_order-1)
     ixR_p_stg_max2(idir,0)=ixCoMmin2-1
     ixR_p_stg_min2(idir,1)=ixCoMmin2
     ixR_p_stg_max2(idir,1)=ixCoMmax2+dixBCo+(interpolation_order-1)
     ixR_p_stg_min2(idir,2)=ixCoMmin2-dixBCo-(interpolation_order-1)
     ixR_p_stg_max2(idir,2)=ixCoMmax2
     ixR_p_stg_min2(idir,3)=ixCoMmax2+1
     ixR_p_stg_max2(idir,3)=ixCoMmax2+dixBCo+(interpolation_order-1)


     ! Perpendicular component
     ixS_p_stg_min2(idir,0)=ixMlo2
     ixS_p_stg_max2(idir,0)=ixMlo2-1+dixBCo+(interpolation_order-1)
     ixS_p_stg_min2(idir,1)=ixMlo2
     ixS_p_stg_max2(idir,1)=ixMlo2-1+nxCo2+dixBCo+(interpolation_order-1)
     ixS_p_stg_min2(idir,2)=ixMhi2+1-nxCo2-dixBCo-(interpolation_order-1)
     ixS_p_stg_max2(idir,2)=ixMhi2
     ixS_p_stg_min2(idir,3)=ixMhi2+1-dixBCo-(interpolation_order-1)
     ixS_p_stg_max2(idir,3)=ixMhi2
 
     ixR_p_stg_min2(idir,0)=ixCoMmin2-dixBCo-(interpolation_order-1)
     ixR_p_stg_max2(idir,0)=ixCoMmin2-1
     ixR_p_stg_min2(idir,1)=ixCoMmin2
     ixR_p_stg_max2(idir,1)=ixCoMmax2+dixBCo+(interpolation_order-1)
     ixR_p_stg_min2(idir,2)=ixCoMmin2-dixBCo-(interpolation_order-1)
     ixR_p_stg_max2(idir,2)=ixCoMmax2
     ixR_p_stg_min2(idir,3)=ixCoMmax2+1
     ixR_p_stg_max2(idir,3)=ixCoMmax2+dixBCo+(interpolation_order-1)


     ! Perpendicular component
     ixS_p_stg_min2(idir,0)=ixMlo2
     ixS_p_stg_max2(idir,0)=ixMlo2-1+dixBCo+(interpolation_order-1)
     ixS_p_stg_min2(idir,1)=ixMlo2
     ixS_p_stg_max2(idir,1)=ixMlo2-1+nxCo2+dixBCo+(interpolation_order-1)
     ixS_p_stg_min2(idir,2)=ixMhi2+1-nxCo2-dixBCo-(interpolation_order-1)
     ixS_p_stg_max2(idir,2)=ixMhi2
     ixS_p_stg_min2(idir,3)=ixMhi2+1-dixBCo-(interpolation_order-1)
     ixS_p_stg_max2(idir,3)=ixMhi2
 
     ixR_p_stg_min2(idir,0)=ixCoMmin2-dixBCo-(interpolation_order-1)
     ixR_p_stg_max2(idir,0)=ixCoMmin2-1
     ixR_p_stg_min2(idir,1)=ixCoMmin2
     ixR_p_stg_max2(idir,1)=ixCoMmax2+dixBCo+(interpolation_order-1)
     ixR_p_stg_min2(idir,2)=ixCoMmin2-dixBCo-(interpolation_order-1)
     ixR_p_stg_max2(idir,2)=ixCoMmax2
     ixR_p_stg_min2(idir,3)=ixCoMmax2+1
     ixR_p_stg_max2(idir,3)=ixCoMmax2+dixBCo+(interpolation_order-1)

    end if
   if (idir.eq.3) then
     ! Parallel components

     ixS_p_stg_min3(idir,0)=ixMlo3-1 ! -1 to make redundant 
     ixS_p_stg_max3(idir,0)=ixMlo3-1+dixBCo
     ixS_p_stg_min3(idir,1)=ixMlo3-1 ! -1 to make redundant 
     ixS_p_stg_max3(idir,1)=ixMlo3-1+nxCo3+dixBCo
     ixS_p_stg_min3(idir,2)=ixMhi3-nxCo3-dixBCo
     ixS_p_stg_max3(idir,2)=ixMhi3
     ixS_p_stg_min3(idir,3)=ixMhi3-dixBCo
     ixS_p_stg_max3(idir,3)=ixMhi3

     ixR_p_stg_min3(idir,0)=ixCoMmin3-1-dixBCo
     ixR_p_stg_max3(idir,0)=ixCoMmin3-1
     ixR_p_stg_min3(idir,1)=ixCoMmin3-1 ! -1 to make redundant 
     ixR_p_stg_max3(idir,1)=ixCoMmax3+dixBCo
     ixR_p_stg_min3(idir,2)=ixCoMmin3-1-dixBCo
     ixR_p_stg_max3(idir,2)=ixCoMmax3
     ixR_p_stg_min3(idir,3)=ixCoMmax3+1-1 ! -1 to make redundant 
     ixR_p_stg_max3(idir,3)=ixCoMmax3+dixBCo


     ixS_p_stg_min3(idir,0)=ixMlo3-1 ! -1 to make redundant 
     ixS_p_stg_max3(idir,0)=ixMlo3-1+dixBCo
     ixS_p_stg_min3(idir,1)=ixMlo3-1 ! -1 to make redundant 
     ixS_p_stg_max3(idir,1)=ixMlo3-1+nxCo3+dixBCo
     ixS_p_stg_min3(idir,2)=ixMhi3-nxCo3-dixBCo
     ixS_p_stg_max3(idir,2)=ixMhi3
     ixS_p_stg_min3(idir,3)=ixMhi3-dixBCo
     ixS_p_stg_max3(idir,3)=ixMhi3

     ixR_p_stg_min3(idir,0)=ixCoMmin3-1-dixBCo
     ixR_p_stg_max3(idir,0)=ixCoMmin3-1
     ixR_p_stg_min3(idir,1)=ixCoMmin3-1 ! -1 to make redundant 
     ixR_p_stg_max3(idir,1)=ixCoMmax3+dixBCo
     ixR_p_stg_min3(idir,2)=ixCoMmin3-1-dixBCo
     ixR_p_stg_max3(idir,2)=ixCoMmax3
     ixR_p_stg_min3(idir,3)=ixCoMmax3+1-1 ! -1 to make redundant 
     ixR_p_stg_max3(idir,3)=ixCoMmax3+dixBCo


     ixS_p_stg_min3(idir,0)=ixMlo3-1 ! -1 to make redundant 
     ixS_p_stg_max3(idir,0)=ixMlo3-1+dixBCo
     ixS_p_stg_min3(idir,1)=ixMlo3-1 ! -1 to make redundant 
     ixS_p_stg_max3(idir,1)=ixMlo3-1+nxCo3+dixBCo
     ixS_p_stg_min3(idir,2)=ixMhi3-nxCo3-dixBCo
     ixS_p_stg_max3(idir,2)=ixMhi3
     ixS_p_stg_min3(idir,3)=ixMhi3-dixBCo
     ixS_p_stg_max3(idir,3)=ixMhi3

     ixR_p_stg_min3(idir,0)=ixCoMmin3-1-dixBCo
     ixR_p_stg_max3(idir,0)=ixCoMmin3-1
     ixR_p_stg_min3(idir,1)=ixCoMmin3-1 ! -1 to make redundant 
     ixR_p_stg_max3(idir,1)=ixCoMmax3+dixBCo
     ixR_p_stg_min3(idir,2)=ixCoMmin3-1-dixBCo
     ixR_p_stg_max3(idir,2)=ixCoMmax3
     ixR_p_stg_min3(idir,3)=ixCoMmax3+1-1 ! -1 to make redundant 
     ixR_p_stg_max3(idir,3)=ixCoMmax3+dixBCo

    else

     ! Perpendicular component
     ixS_p_stg_min3(idir,0)=ixMlo3
     ixS_p_stg_max3(idir,0)=ixMlo3-1+dixBCo+(interpolation_order-1)
     ixS_p_stg_min3(idir,1)=ixMlo3
     ixS_p_stg_max3(idir,1)=ixMlo3-1+nxCo3+dixBCo+(interpolation_order-1)
     ixS_p_stg_min3(idir,2)=ixMhi3+1-nxCo3-dixBCo-(interpolation_order-1)
     ixS_p_stg_max3(idir,2)=ixMhi3
     ixS_p_stg_min3(idir,3)=ixMhi3+1-dixBCo-(interpolation_order-1)
     ixS_p_stg_max3(idir,3)=ixMhi3
 
     ixR_p_stg_min3(idir,0)=ixCoMmin3-dixBCo-(interpolation_order-1)
     ixR_p_stg_max3(idir,0)=ixCoMmin3-1
     ixR_p_stg_min3(idir,1)=ixCoMmin3
     ixR_p_stg_max3(idir,1)=ixCoMmax3+dixBCo+(interpolation_order-1)
     ixR_p_stg_min3(idir,2)=ixCoMmin3-dixBCo-(interpolation_order-1)
     ixR_p_stg_max3(idir,2)=ixCoMmax3
     ixR_p_stg_min3(idir,3)=ixCoMmax3+1
     ixR_p_stg_max3(idir,3)=ixCoMmax3+dixBCo+(interpolation_order-1)


     ! Perpendicular component
     ixS_p_stg_min3(idir,0)=ixMlo3
     ixS_p_stg_max3(idir,0)=ixMlo3-1+dixBCo+(interpolation_order-1)
     ixS_p_stg_min3(idir,1)=ixMlo3
     ixS_p_stg_max3(idir,1)=ixMlo3-1+nxCo3+dixBCo+(interpolation_order-1)
     ixS_p_stg_min3(idir,2)=ixMhi3+1-nxCo3-dixBCo-(interpolation_order-1)
     ixS_p_stg_max3(idir,2)=ixMhi3
     ixS_p_stg_min3(idir,3)=ixMhi3+1-dixBCo-(interpolation_order-1)
     ixS_p_stg_max3(idir,3)=ixMhi3
 
     ixR_p_stg_min3(idir,0)=ixCoMmin3-dixBCo-(interpolation_order-1)
     ixR_p_stg_max3(idir,0)=ixCoMmin3-1
     ixR_p_stg_min3(idir,1)=ixCoMmin3
     ixR_p_stg_max3(idir,1)=ixCoMmax3+dixBCo+(interpolation_order-1)
     ixR_p_stg_min3(idir,2)=ixCoMmin3-dixBCo-(interpolation_order-1)
     ixR_p_stg_max3(idir,2)=ixCoMmax3
     ixR_p_stg_min3(idir,3)=ixCoMmax3+1
     ixR_p_stg_max3(idir,3)=ixCoMmax3+dixBCo+(interpolation_order-1)


     ! Perpendicular component
     ixS_p_stg_min3(idir,0)=ixMlo3
     ixS_p_stg_max3(idir,0)=ixMlo3-1+dixBCo+(interpolation_order-1)
     ixS_p_stg_min3(idir,1)=ixMlo3
     ixS_p_stg_max3(idir,1)=ixMlo3-1+nxCo3+dixBCo+(interpolation_order-1)
     ixS_p_stg_min3(idir,2)=ixMhi3+1-nxCo3-dixBCo-(interpolation_order-1)
     ixS_p_stg_max3(idir,2)=ixMhi3
     ixS_p_stg_min3(idir,3)=ixMhi3+1-dixBCo-(interpolation_order-1)
     ixS_p_stg_max3(idir,3)=ixMhi3
 
     ixR_p_stg_min3(idir,0)=ixCoMmin3-dixBCo-(interpolation_order-1)
     ixR_p_stg_max3(idir,0)=ixCoMmin3-1
     ixR_p_stg_min3(idir,1)=ixCoMmin3
     ixR_p_stg_max3(idir,1)=ixCoMmax3+dixBCo+(interpolation_order-1)
     ixR_p_stg_min3(idir,2)=ixCoMmin3-dixBCo-(interpolation_order-1)
     ixR_p_stg_max3(idir,2)=ixCoMmax3
     ixR_p_stg_min3(idir,3)=ixCoMmax3+1
     ixR_p_stg_max3(idir,3)=ixCoMmax3+dixBCo+(interpolation_order-1)

    end if
  
   end do


!! Sizes for multi-resolution communications
do i3=-1,1
do i2=-1,1
do i1=-1,1

   ! Cell-centred variables
   sizes_r_send(i1,i2,i3)=(nwflux+nwaux)*(ixS_r_max1(i1)-ixS_r_min1(i1)+1)*&
      (ixS_r_max2(i2)-ixS_r_min2(i2)+1)*(ixS_r_max3(i3)-ixS_r_min3(i3)+1)
   sizes_r_send_total(i1,i2,i3)=sizes_r_send(i1,i2,i3)


   ! Staggered (face-allocated) variables
   do idir=1,3

     sizes_r_send_stg(idir,i1,i2,i3)=(ixS_r_stg_max1(idir,i1)-ixS_r_stg_min1&
        (idir,i1)+1)*(ixS_r_stg_max2(idir,i2)-ixS_r_stg_min2(idir,&
        i2)+1)*(ixS_r_stg_max3(idir,i3)-ixS_r_stg_min3(idir,i3)+1)

     sizes_r_send_total(i1,i2,i3)=sizes_r_send_total(i1,i2,&
        i3)+sizes_r_send_stg(idir,i1,i2,i3)

   end do


end do
end do
end do

do i3=0,3
do i2=0,3
do i1=0,3

   ! Cell-centred variables
   sizes_r_recv(i1,i2,i3)=(nwflux+nwaux)*(ixR_r_max1(i1)-ixR_r_min1(i1)+1)*&
      (ixR_r_max2(i2)-ixR_r_min2(i2)+1)*(ixR_r_max3(i3)-ixR_r_min3(i3)+1)
   sizes_p_send(i1,i2,i3)=(nwflux+nwaux)*(ixS_p_max1(i1)-ixS_p_min1(i1)+1)*&
      (ixS_p_max2(i2)-ixS_p_min2(i2)+1)*(ixS_p_max3(i3)-ixS_p_min3(i3)+1)
   sizes_p_recv(i1,i2,i3)=(nwflux+nwaux)*(ixR_p_max1(i1)-ixR_p_min1(i1)+1)*&
      (ixR_p_max2(i2)-ixR_p_min2(i2)+1)*(ixR_p_max3(i3)-ixR_p_min3(i3)+1)

   sizes_r_recv_total(i1,i2,i3)=sizes_r_recv(i1,i2,i3)
   sizes_p_send_total(i1,i2,i3)=sizes_p_send(i1,i2,i3)
   sizes_p_recv_total(i1,i2,i3)=sizes_p_recv(i1,i2,i3)


   ! Staggered (face-allocated) variables
   do idir=1,3

     sizes_r_recv_stg(idir,i1,i2,i3)=(ixR_r_stg_max1(idir,i1)-ixR_r_stg_min1&
        (idir,i1)+1)*(ixR_r_stg_max2(idir,i2)-ixR_r_stg_min2(idir,&
        i2)+1)*(ixR_r_stg_max3(idir,i3)-ixR_r_stg_min3(idir,i3)+1)
     sizes_p_send_stg(idir,i1,i2,i3)=(ixS_p_stg_max1(idir,i1)-ixS_p_stg_min1&
        (idir,i1)+1)*(ixS_p_stg_max2(idir,i2)-ixS_p_stg_min2(idir,&
        i2)+1)*(ixS_p_stg_max3(idir,i3)-ixS_p_stg_min3(idir,i3)+1)
     sizes_p_recv_stg(idir,i1,i2,i3)=(ixR_p_stg_max1(idir,i1)-ixR_p_stg_min1&
        (idir,i1)+1)*(ixR_p_stg_max2(idir,i2)-ixR_p_stg_min2(idir,&
        i2)+1)*(ixR_p_stg_max3(idir,i3)-ixR_p_stg_min3(idir,i3)+1)

     sizes_r_recv_total(i1,i2,i3)=sizes_r_recv_total(i1,i2,&
        i3)+sizes_r_recv_stg(idir,i1,i2,i3)
     sizes_p_send_total(i1,i2,i3)=sizes_p_send_total(i1,i2,&
        i3)+sizes_p_send_stg(idir,i1,i2,i3)
     sizes_p_recv_total(i1,i2,i3)=sizes_p_recv_total(i1,i2,&
        i3)+sizes_p_recv_stg(idir,i1,i2,i3)

   end do


end do
end do
end do

end if !! levmin/=levmax

!! Calculate size of communication buffers
!! Loop over local grids and see what is needed
!! This could be done in build_connectivity, in connectivity.t

nrecv_gc_srl=0
nsend_gc_srl=0
nbuff_gc_send_srl=0
nbuff_gc_recv_srl=0

nrecv_gc_r=0
nsend_gc_r=0
nbuff_gc_send_r=0
nbuff_gc_recv_r=0

nrecv_gc_p=0
nsend_gc_p=0
nbuff_gc_send_p=0
nbuff_gc_recv_p=0

do iigrid=1,igridstail; igrid=igrids(iigrid);
   do i3=-1,1
   do i2=-1,1
   do i1=-1,1
      if (i1==0.and.i2==0.and.i3==0) cycle
      my_neighbor_type=neighbor_type(i1,i2,i3,igrid)
      ineighbor=neighbor(1,i1,i2,i3,igrid)
      ipe_neighbor=neighbor(2,i1,i2,i3,igrid)

      select case(my_neighbor_type)
      case(3) !! same resolution level
        if (ipe_neighbor.ne.mype) then
        nsend_gc_srl=nsend_gc_srl+1
        nrecv_gc_srl=nrecv_gc_srl+1
        nbuff_gc_send_srl=nbuff_gc_send_srl+sizes_srl_send_total(i1,i2,i3)
        nbuff_gc_recv_srl=nbuff_gc_recv_srl+sizes_srl_recv_total(i1,i2,i3)
        end if
      case(2) !! Coarser neighbor
        if (ipe_neighbor.ne.mype) then
          ic1=1+modulo(node(pig1_,igrid)-1,2)
          ic2=1+modulo(node(pig2_,igrid)-1,2)
          ic3=1+modulo(node(pig3_,igrid)-1,2);
          ! This is the local index of the grid.
          ! Depending on it, sometimes communication does not happen.
          ! Consider for example the configuration
          !
          !       F2|
          !       --| C
          !       F1|
          !
          ! The upper corner of F1 does not need to be communicated
          ! to C because there is already F2.

          if ((i1==0.or.i1==2*ic1-3).and.(i2==0.or.i2==2*ic2-3).and.(i3==0&
             .or.i3==2*ic3-3)) then
            nsend_gc_r=nsend_gc_r+1
            nrecv_gc_p=nrecv_gc_p+1
            nbuff_gc_send_r=nbuff_gc_send_r+sizes_r_send_total(i1,i2,i3)
          ! This is the local index of the prolonged ghost zone
            inc1=ic1+i1;inc2=ic2+i2;inc3=ic3+i3;
            nbuff_gc_recv_p=nbuff_gc_recv_p+sizes_p_recv_total(inc1,inc2,inc3)
          end if
        end if
      case(4) !! Finer neighbor
        ! Loop over the local indices of children ic^D
        ! and calculate local indices of ghost zone inc^D.
        do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
           inc3=2*i3+ic3
        do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
           inc2=2*i2+ic2
        do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
           inc1=2*i1+ic1
           ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
           if (ipe_neighbor.ne.mype) then
             nsend_gc_p=nsend_gc_p+1
             nrecv_gc_r=nrecv_gc_r+1
             ! Although the indices change, the size of the send buffer
             ! does not depend on whether there is a pole.
             nbuff_gc_send_p=nbuff_gc_send_p+sizes_p_send_total(inc1,inc2,&
                inc3)
             nbuff_gc_recv_r=nbuff_gc_recv_r+sizes_r_recv_total(inc1,inc2,&
                inc3)
           end if
        end do
        end do
        end do
      end select
   end do
   end do
   end do
end do

!! Allocate communication buffers, initialize MPI requests
!! ... same resolution level ...
if (nrecv_gc_srl>0) then
   allocate(recvbuffer_srl(nbuff_gc_recv_srl))
   allocate(recvstatus_srl(MPI_STATUS_SIZE,nrecv_gc_srl),recvrequest_srl&
      (nrecv_gc_srl))
   recvrequest_srl=MPI_REQUEST_NULL

   !! Make 'task' fill ghostcells srl

   ibuf_recv_srl=1
   irecv_srl=0

   do iigrid=1,igridstail; igrid=igrids(iigrid);
     call make_task_fill_gc_srl(igrid)
   end do

end if

if (nsend_gc_srl>0) then
   allocate(sendbuffer_srl(nbuff_gc_send_srl))
   allocate(sendstatus_srl(MPI_STATUS_SIZE,nsend_gc_srl),sendrequest_srl&
      (nsend_gc_srl))
   sendrequest_srl=MPI_REQUEST_NULL

   ibuf_send_srl=1
   isend_srl=0

end if

! ... restrict ...
if (nrecv_gc_r>0) then
   allocate(recvbuffer_r(nbuff_gc_recv_r))
   allocate(recvstatus_r(MPI_STATUS_SIZE,nrecv_gc_r),recvrequest_r&
      (nrecv_gc_r))
   recvrequest_r=MPI_REQUEST_NULL

   !! Make 'task' fill ghostcells srl

   ibuf_recv_r=1
   irecv_r=0

   do iigrid=1,igridstail; igrid=igrids(iigrid);
     call make_task_fill_gc_r(igrid)
   end do

end if

if (nsend_gc_r>0) then
   allocate(sendbuffer_r(nbuff_gc_send_r))
   allocate(sendstatus_r(MPI_STATUS_SIZE,nsend_gc_r),sendrequest_r&
      (nsend_gc_r))
   sendrequest_r=MPI_REQUEST_NULL

   ibuf_send_r=1
   isend_r=0

end if

! ... prolong ...
if (nrecv_gc_p>0) then
   allocate(recvbuffer_p(nbuff_gc_recv_p))
   allocate(recvstatus_p(MPI_STATUS_SIZE,nrecv_gc_p),recvrequest_p&
      (nrecv_gc_p))
   recvrequest_p=MPI_REQUEST_NULL

   !! Make 'task' fill ghostcells srl

   ibuf_recv_p=1
   irecv_p=0

   do iigrid=1,igridstail; igrid=igrids(iigrid);
     call make_task_fill_gc_p(igrid)
   end do

end if

if (nsend_gc_p>0) then
   allocate(sendbuffer_p(nbuff_gc_send_p))
   allocate(sendstatus_p(MPI_STATUS_SIZE,nsend_gc_p),sendrequest_p&
      (nsend_gc_p))
   sendrequest_p=MPI_REQUEST_NULL

   ibuf_send_p=1
   isend_p=0

end if



! Allocate buffers for reversing phi at poles
allocate(pole_buf%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux+nwaux))

allocate(pole_buf_stg%w(ixGlo1-1:ixGhi1,ixGlo2-1:ixGhi2,ixGlo3-1:ixGhi3,1:3))




end subroutine init_comm_gc
!===============================================================================
subroutine make_task_fill_gc_srl(igrid)
use mod_comm_gc
use mod_amrvacdef

! In a loop, place MPI receive requests, updating the part of the receive
! buffer where ghost cells will be written
! Note: ibuf_recv and irecv_srl were already initialized in init_comm_gc_srl,
! here they are just updated 
integer, intent(in) :: igrid
integer :: i1,i2,i3, ineighbor, ipe_neighbor, my_neighbor_type
!-------------------------------------------------------------------------------
call set_tmpGlobals(igrid)

do i3=-1,1
do i2=-1,1
do i1=-1,1
   if (i1==0.and.i2==0.and.i3==0) cycle
   my_neighbor_type=neighbor_type(i1,i2,i3,igrid)
   select case (my_neighbor_type) !! Could the restrict receive be also here?
   case (3) !!! Case 3 = srl
     ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
       if (ipe_neighbor/=mype) then
         irecv_srl=irecv_srl+1
         itag=(3**3+4**3)*(igrid-1)+(i1+1)*3**(1-1)+(i2+1)*3**(2-1)+(i3+1)*3&
            **(3-1)
         call MPI_IRECV(recvbuffer_srl(ibuf_recv_srl),sizes_srl_recv_total(i1,&
            i2,i3), MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
            recvrequest_srl(irecv_srl),ierrmpi)
         !call add_task_to_list() 
         ibuf_recv_srl=ibuf_recv_srl+sizes_srl_recv_total(i1,i2,i3)
       end if
   end select
end do
end do
end do

end subroutine make_task_fill_gc_srl
!===============================================================================
subroutine fill_gc_srl
use mod_comm_gc
use mod_amrvacdef

integer, save :: i1,i2,i3, igrid, iigrid
integer, save :: ipole, ineighbor, ipe_neighbor, my_neighbor_type
!-------------------------------------------------------------------------------


! In a loop with the same order as that in make_task_fill_gc_srl,
! directly fill the ghost cells that are in the same cpu and 
! unpack the receive buffer to fill those that are not.

!$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(igrid)
do iigrid=1,igridstail; igrid=igrids(iigrid);
   call set_tmpGlobals(igrid)
         
  do i3=-1,1
  do i2=-1,1
  do i1=-1,1
     if (i1==0.and.i2==0.and.i3==0) cycle
     ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
     if (ipe_neighbor.ne.mype) cycle
     
     my_neighbor_type=neighbor_type(i1,i2,i3,igrid)
     select case (my_neighbor_type) !! Could the restrict receive be also here?
     case (3) !!! Case 3 = srl
        call bc_fill_srl_omp(igrid,i1,i2,i3)
     end select
  end do
  end do
  end do
  
end do
!$OMP END PARALLEL DO

! Wait for the receive buffer to be complete
if (nrecv_gc_srl>0) then
   call MPI_WAITALL(nrecv_gc_srl,recvrequest_srl,recvstatus_srl,ierrmpi)
   ibuf_recv_srl=1
end if
if (nsend_gc_srl>0) then
   call MPI_WAITALL(nsend_gc_srl,sendrequest_srl,sendstatus_srl,ierrmpi)
end if

do iigrid=1,igridstail; igrid=igrids(iigrid);
   call set_tmpGlobals(igrid)
   
  do i3=-1,1
  do i2=-1,1
  do i1=-1,1
     if (i1==0.and.i2==0.and.i3==0) cycle
     ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
     if (ipe_neighbor.eq.mype) cycle
     
     my_neighbor_type=neighbor_type(i1,i2,i3,igrid)
     select case (my_neighbor_type) !! Could the restrict receive be also here?
     case (3) !!! Case 3 = srl
        call bc_fill_srl_mpi(igrid,i1,i2,i3)
     end select
  end do
  end do
  end do
  
end do

! Deallocate communication buffers
if (nrecv_gc_srl>0) deallocate(recvbuffer_srl,recvstatus_srl,recvrequest_srl)
if (nsend_gc_srl>0) deallocate(sendbuffer_srl,sendstatus_srl,sendrequest_srl)

contains
!=============================================================================
! internal procedures
!=============================================================================
subroutine bc_fill_srl
integer :: ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3,ixRmin1,ixRmin2,&
   ixRmin3,ixRmax1,ixRmax2,ixRmax3,n_i1,n_i2,n_i3,ixSsyncmin1,ixSsyncmin2,&
   ixSsyncmin3,ixSsyncmax1,ixSsyncmax2,ixSsyncmax3,ixRsyncmin1,ixRsyncmin2,&
   ixRsyncmin3,ixRsyncmax1,ixRsyncmax2,ixRsyncmax3
integer :: ibufnext,idir

double precision, dimension(pstate(igrid)%ws%ixGmin1:pstate(igrid)%ws%ixGmax1,&
   pstate(igrid)%ws%ixGmin2:pstate(igrid)%ws%ixGmax2,pstate&
   (igrid)%ws%ixGmin3:pstate(igrid)%ws%ixGmax3) :: tmp

integer :: idirect, iw
!-----------------------------------------------------------------------------
ineighbor=neighbor(1,i1,i2,i3,igrid)
ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
ipole=neighbor_pole(i1,i2,i3,igrid)
idirect=abs(i1)+abs(i2)+abs(i3)

!! Now the special treatment of the pole is done here, at the receive step
if (ipole.eq.0) then    
  n_i1=-i1;n_i2=-i2;n_i3=-i3;
  ixRmin1=ixR_srl_min1(i1);ixRmin2=ixR_srl_min2(i2);ixRmin3=ixR_srl_min3(i3)
  ixRmax1=ixR_srl_max1(i1);ixRmax2=ixR_srl_max2(i2);ixRmax3=ixR_srl_max3(i3);
  ixSmin1=ixS_srl_min1(n_i1);ixSmin2=ixS_srl_min2(n_i2)
  ixSmin3=ixS_srl_min3(n_i3);ixSmax1=ixS_srl_max1(n_i1)
  ixSmax2=ixS_srl_max2(n_i2);ixSmax3=ixS_srl_max3(n_i3);

  if (ipe_neighbor.eq.mype) then
    !! Just copy from the other block
!$OMP PARALLEL DO SCHEDULE(static,1)
do iw=1,nwflux+nwaux
      pstate(igrid)%w%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,iw)&
         =pstate(ineighbor)%w%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
         ixSmin3:ixSmax3,iw)
end do
!$OMP END PARALLEL DO

   
!$OMP PARALLEL DO SCHEDULE(static,1)
   do idir=1,3
      ixSmin1=ixS_srl_stg_min1(idir,n_i1);ixSmin2=ixS_srl_stg_min2(idir,n_i2)
      ixSmin3=ixS_srl_stg_min3(idir,n_i3);ixSmax1=ixS_srl_stg_max1(idir,n_i1)
      ixSmax2=ixS_srl_stg_max2(idir,n_i2);ixSmax3=ixS_srl_stg_max3(idir,n_i3);
      ixRmin1=ixR_srl_stg_min1(idir,i1);ixRmin2=ixR_srl_stg_min2(idir,i2)
      ixRmin3=ixR_srl_stg_min3(idir,i3);ixRmax1=ixR_srl_stg_max1(idir,i1)
      ixRmax2=ixR_srl_stg_max2(idir,i2);ixRmax3=ixR_srl_stg_max3(idir,i3);
      
      ! opedit: in principle, averaging together would be desirable.
 !when setting pstate(ineighbor) = pstate(igrid) here this leads to a bug though, which is
      ! not fully understood.  
 !So for now we just exchange twice, meaning we keep either the left or the right representation
 !of the staggerd field between blocks. In any case, the states are sychronized this way.
!      if (idirect .eq. 1) then
 !call indices_for_syncing(idir,i1,i2,i3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3,ixRsyncmin1,ixRsyncmin2,ixRsyncmin3,ixRsyncmax1,ixRsyncmax2,ixRsyncmax3,ixSsyncmin1,ixSsyncmin2,ixSsyncmin3,ixSsyncmax1,ixSsyncmax2,ixSsyncmax3) ! Overwrites ixR, ixS
 !pstate(igrid)%ws%w(ixRsyncmin1:ixRsyncmax1,ixRsyncmin2:ixRsyncmax2,ixRsyncmin3:ixRsyncmax3,idir) = &
 !half*(pstate(igrid)%ws%w(ixRsyncmin1:ixRsyncmax1,ixRsyncmin2:ixRsyncmax2,ixRsyncmin3:ixRsyncmax3,idir) + pstate(ineighbor)%ws%w(ixSsyncmin1:ixSsyncmax1,ixSsyncmin2:ixSsyncmax2,ixSsyncmin3:ixSsyncmax3,idir))
!      end if
      
      pstate(igrid)%ws%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
         idir) = pstate(ineighbor)%ws%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
         ixSmin3:ixSmax3,idir)
   end do
!$OMP END PARALLEL DO
   
  else
    !! Unpack the buffer and fill the ghost cells
    ibufnext=ibuf_recv_srl+sizes_srl_recv(i1,i2,i3)
    pstate(igrid)%w%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
       1:nwflux+nwaux)=reshape(source=recvbuffer_srl(ibuf_recv_srl:ibufnext-1)&
       ,shape=shape(pstate(igrid)%w%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
       ixRmin3:ixRmax3,1:nwflux+nwaux)))
  
    ibuf_recv_srl=ibufnext
  
    
    do idir=1,3
       ixSmin1=ixS_srl_stg_min1(idir,n_i1)
       ixSmin2=ixS_srl_stg_min2(idir,n_i2)
       ixSmin3=ixS_srl_stg_min3(idir,n_i3)
       ixSmax1=ixS_srl_stg_max1(idir,n_i1)
       ixSmax2=ixS_srl_stg_max2(idir,n_i2)
       ixSmax3=ixS_srl_stg_max3(idir,n_i3);
       ixRmin1=ixR_srl_stg_min1(idir,i1);ixRmin2=ixR_srl_stg_min2(idir,i2)
       ixRmin3=ixR_srl_stg_min3(idir,i3);ixRmax1=ixR_srl_stg_max1(idir,i1)
       ixRmax2=ixR_srl_stg_max2(idir,i2);ixRmax3=ixR_srl_stg_max3(idir,i3);

       ibufnext=ibuf_recv_srl+sizes_srl_recv_stg(idir,i1,i2,i3)
       
       tmp(ixSmin1:ixSmax1,ixSmin2:ixSmax2,ixSmin3:ixSmax3) = reshape(source&
          =recvbuffer_srl(ibuf_recv_srl:ibufnext-1),shape=shape(pstate&
          (igrid)%ws%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,ixSmin3:ixSmax3,&
          idir)))       

       if (idirect .eq. 1) then
          call indices_for_syncing(idir,i1,i2,i3,ixRmin1,ixRmin2,ixRmin3,&
             ixRmax1,ixRmax2,ixRmax3,ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,&
             ixSmax3,ixRsyncmin1,ixRsyncmin2,ixRsyncmin3,ixRsyncmax1,&
             ixRsyncmax2,ixRsyncmax3,ixSsyncmin1,ixSsyncmin2,ixSsyncmin3,&
             ixSsyncmax1,ixSsyncmax2,ixSsyncmax3) !Overwrites ixR, ixS
          pstate(igrid)%ws%w(ixRsyncmin1:ixRsyncmax1,ixRsyncmin2:ixRsyncmax2,&
             ixRsyncmin3:ixRsyncmax3,idir) = &
               half*(tmp(ixSsyncmin1:ixSsyncmax1,ixSsyncmin2:ixSsyncmax2,&
                  ixSsyncmin3:ixSsyncmax3) + pstate(igrid)%ws%w&
                  (ixRsyncmin1:ixRsyncmax1,ixRsyncmin2:ixRsyncmax2,&
                  ixRsyncmin3:ixRsyncmax3,idir))
       end if

       pstate(igrid)%ws%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
          idir) = tmp(ixSmin1:ixSmax1,ixSmin2:ixSmax2,ixSmin3:ixSmax3)
    
      ibuf_recv_srl=ibufnext
   end do
    
  end if

else ! There is a pole
  ixRmin1=ixR_srl_min1(i1);ixRmin2=ixR_srl_min2(i2);ixRmin3=ixR_srl_min3(i3)
  ixRmax1=ixR_srl_max1(i1);ixRmax2=ixR_srl_max2(i2);ixRmax3=ixR_srl_max3(i3);
  select case (ipole)
  case (1)
     n_i1=i1;n_i2=-i2;n_i3=-i3;
  case (2)
     n_i1=-i1;n_i2=i2;n_i3=-i3;
  case (3)
     n_i1=-i1;n_i2=-i2;n_i3=i3;
  end select
  ixSmin1=ixS_srl_min1(n_i1);ixSmin2=ixS_srl_min2(n_i2)
  ixSmin3=ixS_srl_min3(n_i3);ixSmax1=ixS_srl_max1(n_i1)
  ixSmax2=ixS_srl_max2(n_i2);ixSmax3=ixS_srl_max3(n_i3);

  if (ipe_neighbor==mype) then
    !! Fill ghost cells
    call pole_copy(pstate(igrid)%w,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,&
       ixRmax3,pstate(ineighbor)%w,ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,&
       ixSmax3,ipole,i1,i2,i3)

    
    do idir=1,3
       ixRmin1=ixR_srl_stg_min1(idir,i1);ixRmin2=ixR_srl_stg_min2(idir,i2)
       ixRmin3=ixR_srl_stg_min3(idir,i3);ixRmax1=ixR_srl_stg_max1(idir,i1)
       ixRmax2=ixR_srl_stg_max2(idir,i2);ixRmax3=ixR_srl_stg_max3(idir,i3);
       ixSmin1=ixS_srl_stg_min1(idir,n_i1)
       ixSmin2=ixS_srl_stg_min2(idir,n_i2)
       ixSmin3=ixS_srl_stg_min3(idir,n_i3)
       ixSmax1=ixS_srl_stg_max1(idir,n_i1)
       ixSmax2=ixS_srl_stg_max2(idir,n_i2)
       ixSmax3=ixS_srl_stg_max3(idir,n_i3);
       !! Fill ghost cells
       call pole_copy_stg(pstate(igrid)%ws,ixRmin1,ixRmin2,ixRmin3,ixRmax1,&
          ixRmax2,ixRmax3,pstate(ineighbor)%ws,ixSmin1,ixSmin2,ixSmin3,&
          ixSmax1,ixSmax2,ixSmax3,ipole,i1,i2,i3,idir)
    end do
    

  else
    !! Unpack the buffer and fill an auxiliary array
    ibufnext=ibuf_recv_srl+sizes_srl_recv(i1,i2,i3)
    pole_buf%w=zero
    pole_buf%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,ixSmin3:ixSmax3,&
       1:nwflux+nwaux)=reshape(source=recvbuffer_srl(ibuf_recv_srl:ibufnext-1)&
       ,shape=shape(pstate(igrid)%w%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
       ixSmin3:ixSmax3,1:nwflux+nwaux)))
    ibuf_recv_srl=ibufnext
 
    !! Fill ghost cells
    call pole_copy(pstate(igrid)%w,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,&
       ixRmax3,pole_buf,ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3,ipole,&
       i1,i2,i3)

    
    pole_buf_stg%w=zero
    do idir=1,3
       ixRmin1=ixR_srl_stg_min1(idir,i1);ixRmin2=ixR_srl_stg_min2(idir,i2)
       ixRmin3=ixR_srl_stg_min3(idir,i3);ixRmax1=ixR_srl_stg_max1(idir,i1)
       ixRmax2=ixR_srl_stg_max2(idir,i2);ixRmax3=ixR_srl_stg_max3(idir,i3);
       ixSmin1=ixS_srl_stg_min1(idir,n_i1)
       ixSmin2=ixS_srl_stg_min2(idir,n_i2)
       ixSmin3=ixS_srl_stg_min3(idir,n_i3)
       ixSmax1=ixS_srl_stg_max1(idir,n_i1)
       ixSmax2=ixS_srl_stg_max2(idir,n_i2)
       ixSmax3=ixS_srl_stg_max3(idir,n_i3);

       ibufnext=ibuf_recv_srl+sizes_srl_recv_stg(idir,i1,i2,i3)
       pole_buf_stg%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,ixSmin3:ixSmax3,idir)&
          =reshape(source=recvbuffer_srl(ibuf_recv_srl:ibufnext-1),&
         shape=shape(pstate(igrid)%ws%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
            ixSmin3:ixSmax3,idir)))
       ibuf_recv_srl=ibufnext

       call pole_copy_stg(pstate(igrid)%ws,ixRmin1,ixRmin2,ixRmin3,ixRmax1,&
          ixRmax2,ixRmax3,pole_buf_stg,ixSmin1,ixSmin2,ixSmin3,ixSmax1,&
          ixSmax2,ixSmax3,ipole,i1,i2,i3,idir)

    end do
    


  end if
end if

end subroutine bc_fill_srl
!===============================================================================
subroutine bc_fill_srl_omp(igrid,i1,i2,i3)
  integer :: igrid, i1,i2,i3
  integer :: ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3,ixRmin1,ixRmin2,&
     ixRmin3,ixRmax1,ixRmax2,ixRmax3,n_i1,n_i2,n_i3,ixSsyncmin1,ixSsyncmin2,&
     ixSsyncmin3,ixSsyncmax1,ixSsyncmax2,ixSsyncmax3,ixRsyncmin1,ixRsyncmin2,&
     ixRsyncmin3,ixRsyncmax1,ixRsyncmax2,ixRsyncmax3
  integer :: ibufnext,idir
  
  double precision, dimension(pstate(igrid)%ws%ixGmin1:pstate&
     (igrid)%ws%ixGmax1,pstate(igrid)%ws%ixGmin2:pstate(igrid)%ws%ixGmax2,&
     pstate(igrid)%ws%ixGmin3:pstate(igrid)%ws%ixGmax3) :: tmp
  
  !-----------------------------------------------------------------------------
  ineighbor=neighbor(1,i1,i2,i3,igrid)
  ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
  ipole=neighbor_pole(i1,i2,i3,igrid)

  !! Now the special treatment of the pole is done here, at the receive step
  if (ipole.eq.0) then
     
     n_i1=-i1;n_i2=-i2;n_i3=-i3;
     ixRmin1=ixR_srl_min1(i1);ixRmin2=ixR_srl_min2(i2)
     ixRmin3=ixR_srl_min3(i3);ixRmax1=ixR_srl_max1(i1)
     ixRmax2=ixR_srl_max2(i2);ixRmax3=ixR_srl_max3(i3);
     ixSmin1=ixS_srl_min1(n_i1);ixSmin2=ixS_srl_min2(n_i2)
     ixSmin3=ixS_srl_min3(n_i3);ixSmax1=ixS_srl_max1(n_i1)
     ixSmax2=ixS_srl_max2(n_i2);ixSmax3=ixS_srl_max3(n_i3);

     !! Just copy from the other block
     pstate(igrid)%w%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
        1:nwflux+nwaux)=pstate(ineighbor)%w%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
        ixSmin3:ixSmax3,1:nwflux+nwaux)

     
     do idir=1,3
        ixSmin1=ixS_srl_stg_min1(idir,n_i1)
        ixSmin2=ixS_srl_stg_min2(idir,n_i2)
        ixSmin3=ixS_srl_stg_min3(idir,n_i3)
        ixSmax1=ixS_srl_stg_max1(idir,n_i1)
        ixSmax2=ixS_srl_stg_max2(idir,n_i2)
        ixSmax3=ixS_srl_stg_max3(idir,n_i3);
        ixRmin1=ixR_srl_stg_min1(idir,i1);ixRmin2=ixR_srl_stg_min2(idir,i2)
        ixRmin3=ixR_srl_stg_min3(idir,i3);ixRmax1=ixR_srl_stg_max1(idir,i1)
        ixRmax2=ixR_srl_stg_max2(idir,i2);ixRmax3=ixR_srl_stg_max3(idir,i3);
        
        pstate(igrid)%ws%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
           idir) = pstate(ineighbor)%ws%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
           ixSmin3:ixSmax3,idir)
     end do
     
     
  else ! There is a pole
     
     ixRmin1=ixR_srl_min1(i1);ixRmin2=ixR_srl_min2(i2)
     ixRmin3=ixR_srl_min3(i3);ixRmax1=ixR_srl_max1(i1)
     ixRmax2=ixR_srl_max2(i2);ixRmax3=ixR_srl_max3(i3);
     select case (ipole)
        case (1)
        n_i1=i1;n_i2=-i2;n_i3=-i3;
        case (2)
        n_i1=-i1;n_i2=i2;n_i3=-i3;
        case (3)
        n_i1=-i1;n_i2=-i2;n_i3=i3;
     end select
     ixSmin1=ixS_srl_min1(n_i1);ixSmin2=ixS_srl_min2(n_i2)
     ixSmin3=ixS_srl_min3(n_i3);ixSmax1=ixS_srl_max1(n_i1)
     ixSmax2=ixS_srl_max2(n_i2);ixSmax3=ixS_srl_max3(n_i3);

     !! Fill ghost cells
     call pole_copy(pstate(igrid)%w,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,&
        ixRmax3,pstate(ineighbor)%w,ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,&
        ixSmax3,ipole,i1,i2,i3)

     
     do idir=1,3
        ixRmin1=ixR_srl_stg_min1(idir,i1);ixRmin2=ixR_srl_stg_min2(idir,i2)
        ixRmin3=ixR_srl_stg_min3(idir,i3);ixRmax1=ixR_srl_stg_max1(idir,i1)
        ixRmax2=ixR_srl_stg_max2(idir,i2);ixRmax3=ixR_srl_stg_max3(idir,i3);
        ixSmin1=ixS_srl_stg_min1(idir,n_i1)
        ixSmin2=ixS_srl_stg_min2(idir,n_i2)
        ixSmin3=ixS_srl_stg_min3(idir,n_i3)
        ixSmax1=ixS_srl_stg_max1(idir,n_i1)
        ixSmax2=ixS_srl_stg_max2(idir,n_i2)
        ixSmax3=ixS_srl_stg_max3(idir,n_i3);
        !! Fill ghost cells
        call pole_copy_stg(pstate(igrid)%ws,ixRmin1,ixRmin2,ixRmin3,ixRmax1,&
           ixRmax2,ixRmax3,pstate(ineighbor)%ws,ixSmin1,ixSmin2,ixSmin3,&
           ixSmax1,ixSmax2,ixSmax3,ipole,i1,i2,i3,idir)
     end do
     
  end if
  
end subroutine bc_fill_srl_omp
!=============================================================================
subroutine bc_fill_srl_mpi(igrid,i1,i2,i3)
  integer :: igrid, i1,i2,i3
  integer :: ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3,ixRmin1,ixRmin2,&
     ixRmin3,ixRmax1,ixRmax2,ixRmax3,n_i1,n_i2,n_i3,ixSsyncmin1,ixSsyncmin2,&
     ixSsyncmin3,ixSsyncmax1,ixSsyncmax2,ixSsyncmax3,ixRsyncmin1,ixRsyncmin2,&
     ixRsyncmin3,ixRsyncmax1,ixRsyncmax2,ixRsyncmax3
  integer :: ibufnext,idir
  
  double precision, dimension(pstate(igrid)%ws%ixGmin1:pstate&
     (igrid)%ws%ixGmax1,pstate(igrid)%ws%ixGmin2:pstate(igrid)%ws%ixGmax2,&
     pstate(igrid)%ws%ixGmin3:pstate(igrid)%ws%ixGmax3) :: tmp
  
  integer :: idirect
  !-----------------------------------------------------------------------------
  ineighbor=neighbor(1,i1,i2,i3,igrid)
  ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
  ipole=neighbor_pole(i1,i2,i3,igrid)
  idirect=abs(i1)+abs(i2)+abs(i3)

  !! Now the special treatment of the pole is done here, at the receive step
  if (ipole.eq.0) then

     n_i1=-i1;n_i2=-i2;n_i3=-i3;
     ixRmin1=ixR_srl_min1(i1);ixRmin2=ixR_srl_min2(i2)
     ixRmin3=ixR_srl_min3(i3);ixRmax1=ixR_srl_max1(i1)
     ixRmax2=ixR_srl_max2(i2);ixRmax3=ixR_srl_max3(i3);
     ixSmin1=ixS_srl_min1(n_i1);ixSmin2=ixS_srl_min2(n_i2)
     ixSmin3=ixS_srl_min3(n_i3);ixSmax1=ixS_srl_max1(n_i1)
     ixSmax2=ixS_srl_max2(n_i2);ixSmax3=ixS_srl_max3(n_i3);

     !! Unpack the buffer and fill the ghost cells
     ibufnext=ibuf_recv_srl+sizes_srl_recv(i1,i2,i3)
     pstate(igrid)%w%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
        1:nwflux+nwaux)=reshape(source=recvbuffer_srl&
        (ibuf_recv_srl:ibufnext-1),shape=shape(pstate(igrid)%w%w&
        (ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,1:nwflux+nwaux)))

     ibuf_recv_srl=ibufnext

     
     do idir=1,3
        ixSmin1=ixS_srl_stg_min1(idir,n_i1)
        ixSmin2=ixS_srl_stg_min2(idir,n_i2)
        ixSmin3=ixS_srl_stg_min3(idir,n_i3)
        ixSmax1=ixS_srl_stg_max1(idir,n_i1)
        ixSmax2=ixS_srl_stg_max2(idir,n_i2)
        ixSmax3=ixS_srl_stg_max3(idir,n_i3);
        ixRmin1=ixR_srl_stg_min1(idir,i1);ixRmin2=ixR_srl_stg_min2(idir,i2)
        ixRmin3=ixR_srl_stg_min3(idir,i3);ixRmax1=ixR_srl_stg_max1(idir,i1)
        ixRmax2=ixR_srl_stg_max2(idir,i2);ixRmax3=ixR_srl_stg_max3(idir,i3);

        ibufnext=ibuf_recv_srl+sizes_srl_recv_stg(idir,i1,i2,i3)

        tmp(ixSmin1:ixSmax1,ixSmin2:ixSmax2,ixSmin3:ixSmax3) = reshape(source&
           =recvbuffer_srl(ibuf_recv_srl:ibufnext-1),shape=shape(pstate&
           (igrid)%ws%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,ixSmin3:ixSmax3,&
           idir)))       

        if (idirect .eq. 1) then
           call indices_for_syncing(idir,i1,i2,i3,ixRmin1,ixRmin2,ixRmin3,&
              ixRmax1,ixRmax2,ixRmax3,ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,&
              ixSmax3,ixRsyncmin1,ixRsyncmin2,ixRsyncmin3,ixRsyncmax1,&
              ixRsyncmax2,ixRsyncmax3,ixSsyncmin1,ixSsyncmin2,ixSsyncmin3,&
              ixSsyncmax1,ixSsyncmax2,ixSsyncmax3) !Overwrites ixR, ixS
           pstate(igrid)%ws%w(ixRsyncmin1:ixRsyncmax1,ixRsyncmin2:ixRsyncmax2,&
              ixRsyncmin3:ixRsyncmax3,idir) = &
                half*(tmp(ixSsyncmin1:ixSsyncmax1,ixSsyncmin2:ixSsyncmax2,&
                   ixSsyncmin3:ixSsyncmax3) + pstate(igrid)%ws%w&
                   (ixRsyncmin1:ixRsyncmax1,ixRsyncmin2:ixRsyncmax2,&
                   ixRsyncmin3:ixRsyncmax3,idir))
        end if

        pstate(igrid)%ws%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
           idir) = tmp(ixSmin1:ixSmax1,ixSmin2:ixSmax2,ixSmin3:ixSmax3)

        ibuf_recv_srl=ibufnext
     end do
     

  else ! There is a pole

     ixRmin1=ixR_srl_min1(i1);ixRmin2=ixR_srl_min2(i2)
     ixRmin3=ixR_srl_min3(i3);ixRmax1=ixR_srl_max1(i1)
     ixRmax2=ixR_srl_max2(i2);ixRmax3=ixR_srl_max3(i3);
     select case (ipole)
        case (1)
        n_i1=i1;n_i2=-i2;n_i3=-i3;
        case (2)
        n_i1=-i1;n_i2=i2;n_i3=-i3;
        case (3)
        n_i1=-i1;n_i2=-i2;n_i3=i3;
     end select
     ixSmin1=ixS_srl_min1(n_i1);ixSmin2=ixS_srl_min2(n_i2)
     ixSmin3=ixS_srl_min3(n_i3);ixSmax1=ixS_srl_max1(n_i1)
     ixSmax2=ixS_srl_max2(n_i2);ixSmax3=ixS_srl_max3(n_i3);


     !! Unpack the buffer and fill an auxiliary array
     ibufnext=ibuf_recv_srl+sizes_srl_recv(i1,i2,i3)
     pole_buf%w=zero
     pole_buf%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,ixSmin3:ixSmax3,&
        1:nwflux+nwaux)=reshape(source=recvbuffer_srl&
        (ibuf_recv_srl:ibufnext-1),shape=shape(pstate(igrid)%w%w&
        (ixSmin1:ixSmax1,ixSmin2:ixSmax2,ixSmin3:ixSmax3,1:nwflux+nwaux)))
     ibuf_recv_srl=ibufnext

     !! Fill ghost cells
     call pole_copy(pstate(igrid)%w,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,&
        ixRmax3,pole_buf,ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3,&
        ipole,i1,i2,i3)

     
     pole_buf_stg%w=zero
     do idir=1,3
        ixRmin1=ixR_srl_stg_min1(idir,i1);ixRmin2=ixR_srl_stg_min2(idir,i2)
        ixRmin3=ixR_srl_stg_min3(idir,i3);ixRmax1=ixR_srl_stg_max1(idir,i1)
        ixRmax2=ixR_srl_stg_max2(idir,i2);ixRmax3=ixR_srl_stg_max3(idir,i3);
        ixSmin1=ixS_srl_stg_min1(idir,n_i1)
        ixSmin2=ixS_srl_stg_min2(idir,n_i2)
        ixSmin3=ixS_srl_stg_min3(idir,n_i3)
        ixSmax1=ixS_srl_stg_max1(idir,n_i1)
        ixSmax2=ixS_srl_stg_max2(idir,n_i2)
        ixSmax3=ixS_srl_stg_max3(idir,n_i3);

        ibufnext=ibuf_recv_srl+sizes_srl_recv_stg(idir,i1,i2,i3)
        pole_buf_stg%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,ixSmin3:ixSmax3,idir)&
           =reshape(source=recvbuffer_srl(ibuf_recv_srl:ibufnext-1),&
             shape=shape(pstate(igrid)%ws%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
                ixSmin3:ixSmax3,idir)))
        ibuf_recv_srl=ibufnext

        call pole_copy_stg(pstate(igrid)%ws,ixRmin1,ixRmin2,ixRmin3,ixRmax1,&
           ixRmax2,ixRmax3,pole_buf_stg,ixSmin1,ixSmin2,ixSmin3,ixSmax1,&
           ixSmax2,ixSmax3,ipole,i1,i2,i3,idir)

     end do
     

  end if
  
end subroutine bc_fill_srl_mpi
!=============================================================================
subroutine indices_for_syncing(idir,i1,i2,i3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,&
   ixRmax2,ixRmax3,ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3,&
   ixRsyncmin1,ixRsyncmin2,ixRsyncmin3,ixRsyncmax1,ixRsyncmax2,ixRsyncmax3,&
   ixSsyncmin1,ixSsyncmin2,ixSsyncmin3,ixSsyncmax1,ixSsyncmax2,ixSsyncmax3)
  integer, intent(in)       :: i1,i2,i3,idir
  integer, intent(inout)    :: ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,&
     ixRmax3,ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3
  integer, intent(out)      :: ixRsyncmin1,ixRsyncmin2,ixRsyncmin3,&
     ixRsyncmax1,ixRsyncmax2,ixRsyncmax3,ixSsyncmin1,ixSsyncmin2,ixSsyncmin3,&
     ixSsyncmax1,ixSsyncmax2,ixSsyncmax3
  ! .. local ..
  !-----------------------------------------------------------------------------

  ixRsyncmin1=ixRmin1;ixRsyncmin2=ixRmin2;ixRsyncmin3=ixRmin3
  ixRsyncmax1=ixRmax1;ixRsyncmax2=ixRmax2;ixRsyncmax3=ixRmax3;
  ixSsyncmin1=ixSmin1;ixSsyncmin2=ixSmin2;ixSsyncmin3=ixSmin3
  ixSsyncmax1=ixSmax1;ixSsyncmax2=ixSmax2;ixSsyncmax3=ixSmax3;
  
  
  if (i1 .eq. -1 .and. idir .eq. 1) then
     ixRsyncmin1 = ixRmax1
     ixRsyncmax1 = ixRmax1
     ixSsyncmin1 = ixSmax1
     ixSsyncmax1 = ixSmax1
     ixRmax1 = ixRmax1 - 1
     ixSmax1 = ixSmax1 - 1
  else if (i1 .eq. 1 .and. idir .eq. 1) then
     ixRsyncmin1 = ixRmin1
     ixRsyncmax1 = ixRmin1
     ixSsyncmin1 = ixSmin1
     ixSsyncmax1 = ixSmin1
     ixRmin1 = ixRmin1 + 1
     ixSmin1 = ixSmin1 + 1
  end if
  
  
  if (i2 .eq. -1 .and. idir .eq. 2) then
     ixRsyncmin2 = ixRmax2
     ixRsyncmax2 = ixRmax2
     ixSsyncmin2 = ixSmax2
     ixSsyncmax2 = ixSmax2
     ixRmax2 = ixRmax2 - 1
     ixSmax2 = ixSmax2 - 1
  else if (i2 .eq. 1 .and. idir .eq. 2) then
     ixRsyncmin2 = ixRmin2
     ixRsyncmax2 = ixRmin2
     ixSsyncmin2 = ixSmin2
     ixSsyncmax2 = ixSmin2
     ixRmin2 = ixRmin2 + 1
     ixSmin2 = ixSmin2 + 1
  end if
  
  
  if (i3 .eq. -1 .and. idir .eq. 3) then
     ixRsyncmin3 = ixRmax3
     ixRsyncmax3 = ixRmax3
     ixSsyncmin3 = ixSmax3
     ixSsyncmax3 = ixSmax3
     ixRmax3 = ixRmax3 - 1
     ixSmax3 = ixSmax3 - 1
  else if (i3 .eq. 1 .and. idir .eq. 3) then
     ixRsyncmin3 = ixRmin3
     ixRsyncmax3 = ixRmin3
     ixSsyncmin3 = ixSmin3
     ixSsyncmax3 = ixSmin3
     ixRmin3 = ixRmin3 + 1
     ixSmin3 = ixSmin3 + 1
  end if
  

end subroutine indices_for_syncing
!=============================================================================
end subroutine fill_gc_srl
!=============================================================================
subroutine send_gc_srl(igrid)
 use mod_comm_gc
 use mod_amrvacdef
 ! ibuf_send and isend_srl were already initialized in init_comm_gc_srl,
 ! here they are just updated 
 integer :: igrid, i1,i2,i3, ineighbor, ipe_neighbor, my_neighbor_type
 !-------------------------------------------------------------------------------

call set_tmpGlobals(igrid)

do i3=-1,1
do i2=-1,1
do i1=-1,1
   if (i1==0.and.i2==0.and.i3==0) cycle
   my_neighbor_type=neighbor_type(i1,i2,i3,igrid)
   if (my_neighbor_type==3) call bc_send_srl
end do
end do
end do



contains
!=============================================================================
! internal procedures
!=============================================================================
subroutine bc_send_srl
integer :: n_i1,n_i2,n_i3,idir,ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3
integer :: ibufaux,ibufnext,ipole
 ! Auxialiary array to avoid problems with the preprocessor...
 ! (It keeps splitting the line in the wrong place)
integer, dimension(1) :: sizes
!-----------------------------------------------------------------------------
ineighbor=neighbor(1,i1,i2,i3,igrid)
ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
ipole=neighbor_pole(i1,i2,i3,igrid)

!If the neighbor is in another cpu, ...
if (ipe_neighbor.ne.mype) then

  ! The ghost region of the neighbor changes
  ! if there is a pole
  select case (ipole)
  case(0) ! No pole
     n_i1=-i1;n_i2=-i2;n_i3=-i3;
 case (1) ! Pole in this direction
     n_i1=i1;n_i2=-i2;n_i3=-i3;
 case (2) ! Pole in this direction
     n_i1=-i1;n_i2=i2;n_i3=-i3;
 case (3) ! Pole in this direction
     n_i1=-i1;n_i2=-i2;n_i3=i3;
  end select

  ! fill corresponding part of the send buffer...

  ixSmin1=ixS_srl_min1(i1);ixSmin2=ixS_srl_min2(i2);ixSmin3=ixS_srl_min3(i3)
  ixSmax1=ixS_srl_max1(i1);ixSmax2=ixS_srl_max2(i2);ixSmax3=ixS_srl_max3(i3);
  ibufaux=ibuf_send_srl
  ibufnext=ibufaux+sizes_srl_send(i1,i2,i3)

  sizes=(/sizes_srl_send(i1,i2,i3)/)

  sendbuffer_srl(ibufaux:ibufnext-1)=reshape(pstate(igrid)%w%w&
     (ixSmin1:ixSmax1,ixSmin2:ixSmax2,ixSmin3:ixSmax3,1:nwflux+nwaux),sizes)

  ibufaux=ibufnext

  do idir=1,3
    ixSmin1=ixS_srl_stg_min1(idir,i1);ixSmin2=ixS_srl_stg_min2(idir,i2)
    ixSmin3=ixS_srl_stg_min3(idir,i3);ixSmax1=ixS_srl_stg_max1(idir,i1)
    ixSmax2=ixS_srl_stg_max2(idir,i2);ixSmax3=ixS_srl_stg_max3(idir,i3);
    ibufnext=ibufaux+sizes_srl_send_stg(idir,i1,i2,i3)
    sizes=(/sizes_srl_send_stg(idir,i1,i2,i3)/)
    sendbuffer_srl(ibufaux:ibufnext-1)=&
      reshape(pstate(igrid)%ws%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
         ixSmin3:ixSmax3,idir),sizes)   
    ibufaux=ibufnext
  end do


  ! ...and send
  itag=(3**3+4**3)*(ineighbor-1)+(n_i1+1)*3**(1-1)+(n_i2+1)*3&
     **(2-1)+(n_i3+1)*3**(3-1)
  isend_srl=isend_srl+1
  call MPI_ISEND(sendbuffer_srl(ibuf_send_srl),sizes_srl_send_total(i1,i2,i3),&
      MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,sendrequest_srl&
     (isend_srl),ierrmpi)

  ibuf_send_srl=ibufnext

end if
end subroutine bc_send_srl
end subroutine send_gc_srl
!===============================================================================
subroutine make_task_fill_gc_r(igrid)
use mod_comm_gc
use mod_amrvacdef

! In a loop, place MPI receive requests, updating the part of the receive
! buffer where ghost cells will be written
! Note: ibuf_recv_r and irecv_r were already initialized in init_comm_gc,
! here they are just updated 
integer, intent(in) :: igrid
integer :: i1,i2,i3, ic1,ic2,ic3, inc1,inc2,inc3, ineighbor, ipe_neighbor,&
    my_neighbor_type
!-------------------------------------------------------------------------------
call set_tmpGlobals(igrid)

do i3=-1,1
do i2=-1,1
do i1=-1,1
   if (i1==0.and.i2==0.and.i3==0) cycle
   my_neighbor_type=neighbor_type(i1,i2,i3,igrid)
   if (my_neighbor_type.ne.4) cycle
   ! Loop over the local indices of children ic^D
   ! and calculate local indices of ghost zone inc^D.
   do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
      inc3=2*i3+ic3
   do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
      inc2=2*i2+ic2
   do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
      inc1=2*i1+ic1
      ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
      if (ipe_neighbor.ne.mype) then
        irecv_r=irecv_r+1
        itag=(3**3+4**3)*(igrid-1)+3**3+inc1*4**(1-1)+inc2*4**(2-1)+inc3*4&
           **(3-1)
        call MPI_IRECV(recvbuffer_r(ibuf_recv_r),sizes_r_recv_total(inc1,inc2,&
           inc3), MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
           recvrequest_r(irecv_r),ierrmpi)
        ibuf_recv_r=ibuf_recv_r+sizes_r_recv_total(inc1,inc2,inc3)
      end if
   end do
   end do
   end do
end do
end do
end do

end subroutine make_task_fill_gc_r
!===============================================================================
subroutine fill_gc_r
use mod_comm_gc
use mod_amrvacdef

integer :: i1,i2,i3, igrid, iigrid
integer :: ipole, ineighbor, ipe_neighbor, my_neighbor_type
!-------------------------------------------------------------------------------
! Wait for the receive buffer to be complete

if (nrecv_gc_r>0) then
   call MPI_WAITALL(nrecv_gc_r,recvrequest_r,recvstatus_r,ierrmpi)
   ibuf_recv_r=1
end if

! In a loop with the same order as that in make_task_fill_gc_r,
! directly fill the ghost cells that are in the same cpu and 
! unpack the receive buffer to fill those that are not.

do iigrid=1,igridstail; igrid=igrids(iigrid);
  call set_tmpGlobals(igrid)
     
 do i3=-1,1
 do i2=-1,1
 do i1=-1,1
    if (i1==0.and.i2==0.and.i3==0) cycle
    my_neighbor_type=neighbor_type(i1,i2,i3,igrid)
    if (my_neighbor_type.eq.4) then
      call bc_fill_r
    end if
 end do
 end do
 end do
end do


! Wait for the sends to complete and deallocate communication buffers
if (nrecv_gc_r>0) deallocate(recvbuffer_r,recvstatus_r,recvrequest_r)

if (nsend_gc_r>0) then
   call MPI_WAITALL(nsend_gc_r,sendrequest_r,sendstatus_r,ierrmpi)
   deallocate(sendbuffer_r,sendstatus_r,sendrequest_r)
end if

contains
!=============================================================================
! internal procedures
!=============================================================================
subroutine bc_fill_r
integer :: ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3,ixRmin1,ixRmin2,&
   ixRmin3,ixRmax1,ixRmax2,ixRmax3,ic1,ic2,ic3,inc1,inc2,inc3,n_i1,n_i2,n_i3
integer :: ibufnext,idir
!-----------------------------------------------------------------------------
ipole=neighbor_pole(i1,i2,i3,igrid)

! Check first if there is pole
if (ipole.eq.0) then
! Loop over the children ic^D to and their neighbors inc^D
do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
   inc3=2*i3+ic3
do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
   inc2=2*i2+ic2
do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
   inc1=2*i1+ic1
   n_i1=-i1;n_i2=-i2;n_i3=-i3;
   ineighbor=neighbor_child(1,inc1,inc2,inc3,igrid)
   ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
   if (ipe_neighbor.eq.mype) then ! Same processor
     ixSmin1=ixS_r_min1(n_i1);ixSmin2=ixS_r_min2(n_i2)
     ixSmin3=ixS_r_min3(n_i3);ixSmax1=ixS_r_max1(n_i1)
     ixSmax2=ixS_r_max2(n_i2);ixSmax3=ixS_r_max3(n_i3);
     ixRmin1=ixR_r_min1(inc1);ixRmin2=ixR_r_min2(inc2)
     ixRmin3=ixR_r_min3(inc3);ixRmax1=ixR_r_max1(inc1)
     ixRmax2=ixR_r_max2(inc2);ixRmax3=ixR_r_max3(inc3);
     pstate(igrid)%w%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
        1:nwflux+nwaux)=pstateCo(ineighbor)%w%w(ixSmin1:ixSmax1,&
        ixSmin2:ixSmax2,ixSmin3:ixSmax3,1:nwflux+nwaux)
   
     do idir=1,3
        ixSmin1=ixS_r_stg_min1(idir,n_i1);ixSmin2=ixS_r_stg_min2(idir,n_i2)
        ixSmin3=ixS_r_stg_min3(idir,n_i3);ixSmax1=ixS_r_stg_max1(idir,n_i1)
        ixSmax2=ixS_r_stg_max2(idir,n_i2);ixSmax3=ixS_r_stg_max3(idir,n_i3);
        ixRmin1=ixR_r_stg_min1(idir,inc1);ixRmin2=ixR_r_stg_min2(idir,inc2)
        ixRmin3=ixR_r_stg_min3(idir,inc3);ixRmax1=ixR_r_stg_max1(idir,inc1)
        ixRmax2=ixR_r_stg_max2(idir,inc2);ixRmax3=ixR_r_stg_max3(idir,inc3);
        pstate(igrid)%ws%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
           idir)=pstateCo(ineighbor)%ws%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
           ixSmin3:ixSmax3,idir)
    !    print *, 'idir',idir,'n_i',n_i1,n_i2,n_i3
    !    print *, 'ixS',ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3
    !    print *, 'idir',idir,'inc',inc1,inc2,inc3
    !    print *, 'ixR',ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3
 !call printarray(ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,pstate(igrid)%ws%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,idir))
     end do
   

   else ! Different processor
     ixRmin1=ixR_r_min1(inc1);ixRmin2=ixR_r_min2(inc2)
     ixRmin3=ixR_r_min3(inc3);ixRmax1=ixR_r_max1(inc1)
     ixRmax2=ixR_r_max2(inc2);ixRmax3=ixR_r_max3(inc3);
     !! Unpack the buffer and fill the ghost cells
     ibufnext=ibuf_recv_r+sizes_r_recv(inc1,inc2,inc3)
     pstate(igrid)%w%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
        1:nwflux+nwaux)=reshape(source=recvbuffer_r(ibuf_recv_r:ibufnext-1),&
        shape=shape(pstate(igrid)%w%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
        ixRmin3:ixRmax3,1:nwflux+nwaux)))
   
     ibuf_recv_r=ibufnext

    
     do idir=1,3
       ixRmin1=ixR_r_stg_min1(idir,inc1);ixRmin2=ixR_r_stg_min2(idir,inc2)
       ixRmin3=ixR_r_stg_min3(idir,inc3);ixRmax1=ixR_r_stg_max1(idir,inc1)
       ixRmax2=ixR_r_stg_max2(idir,inc2);ixRmax3=ixR_r_stg_max3(idir,inc3);
       ibufnext=ibuf_recv_r+sizes_r_recv_stg(idir,inc1,inc2,inc3)
       pstate(igrid)%ws%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
          idir)=reshape(source=recvbuffer_r(ibuf_recv_r:ibufnext-1),shape&
          =shape(pstate(igrid)%ws%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
          ixRmin3:ixRmax3,idir)))
       ibuf_recv_r=ibufnext
     end do
    
   end if
end do
end do
end do

else !! There is a pole

do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
   inc3=2*i3+ic3
do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
   inc2=2*i2+ic2
do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
   inc1=2*i1+ic1
   select case(ipole)
  case (1)
     n_i1=i1;n_i2=-i2;n_i3=-i3;
  case (2)
     n_i1=-i1;n_i2=i2;n_i3=-i3;
  case (3)
     n_i1=-i1;n_i2=-i2;n_i3=i3;
   end select
   ineighbor=neighbor_child(1,inc1,inc2,inc3,igrid)
   ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
   if (ipe_neighbor.eq.mype) then ! Same processor
     ixSmin1=ixS_r_min1(n_i1);ixSmin2=ixS_r_min2(n_i2)
     ixSmin3=ixS_r_min3(n_i3);ixSmax1=ixS_r_max1(n_i1)
     ixSmax2=ixS_r_max2(n_i2);ixSmax3=ixS_r_max3(n_i3);
     ixRmin1=ixR_r_min1(inc1);ixRmin2=ixR_r_min2(inc2)
     ixRmin3=ixR_r_min3(inc3);ixRmax1=ixR_r_max1(inc1)
     ixRmax2=ixR_r_max2(inc2);ixRmax3=ixR_r_max3(inc3);
     !! Fill ghost cells
     call pole_copy(pstate(igrid)%w,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,&
        ixRmax3,pstateCo(ineighbor)%w,ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,&
        ixSmax3,ipole,i1,i2,i3)

   
     do idir=1,3
        ixSmin1=ixS_r_stg_min1(idir,n_i1);ixSmin2=ixS_r_stg_min2(idir,n_i2)
        ixSmin3=ixS_r_stg_min3(idir,n_i3);ixSmax1=ixS_r_stg_max1(idir,n_i1)
        ixSmax2=ixS_r_stg_max2(idir,n_i2);ixSmax3=ixS_r_stg_max3(idir,n_i3);
        ixRmin1=ixR_r_stg_min1(idir,inc1);ixRmin2=ixR_r_stg_min2(idir,inc2)
        ixRmin3=ixR_r_stg_min3(idir,inc3);ixRmax1=ixR_r_stg_max1(idir,inc1)
        ixRmax2=ixR_r_stg_max2(idir,inc2);ixRmax3=ixR_r_stg_max3(idir,inc3);
        !! Fill ghost cells
        call pole_copy_stg(pstate(igrid)%ws,ixRmin1,ixRmin2,ixRmin3,ixRmax1,&
           ixRmax2,ixRmax3,pstateCo(ineighbor)%ws,ixSmin1,ixSmin2,ixSmin3,&
           ixSmax1,ixSmax2,ixSmax3,ipole,i1,i2,i3,idir)

     end do
   
   else ! Different processor
     !! Unpack the buffer and fill an auxiliary array
     ixSmin1=ixS_r_min1(n_i1);ixSmin2=ixS_r_min2(n_i2)
     ixSmin3=ixS_r_min3(n_i3);ixSmax1=ixS_r_max1(n_i1)
     ixSmax2=ixS_r_max2(n_i2);ixSmax3=ixS_r_max3(n_i3);
     ixRmin1=ixR_r_min1(inc1);ixRmin2=ixR_r_min2(inc2)
     ixRmin3=ixR_r_min3(inc3);ixRmax1=ixR_r_max1(inc1)
     ixRmax2=ixR_r_max2(inc2);ixRmax3=ixR_r_max3(inc3);

     ibufnext=ibuf_recv_r+sizes_r_recv(inc1,inc2,inc3)
     pole_buf%w=zero
     pole_buf%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
        1:nwflux+nwaux)=reshape(source=recvbuffer_r(ibuf_recv_r:ibufnext-1),&
        shape=shape(pstate(igrid)%w%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
        ixRmin3:ixRmax3,1:nwflux+nwaux)))
     ibuf_recv_r=ibufnext
 
     !! Fill ghost cells
     call pole_copy(pstate(igrid)%w,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,&
        ixRmax3,pole_buf,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,&
        ipole,i1,i2,i3)

     
     pole_buf_stg%w=zero
     do idir=1,3
        ixSmin1=ixS_r_stg_min1(idir,n_i1);ixSmin2=ixS_r_stg_min2(idir,n_i2)
        ixSmin3=ixS_r_stg_min3(idir,n_i3);ixSmax1=ixS_r_stg_max1(idir,n_i1)
        ixSmax2=ixS_r_stg_max2(idir,n_i2);ixSmax3=ixS_r_stg_max3(idir,n_i3);
        ixRmin1=ixR_r_stg_min1(idir,inc1);ixRmin2=ixR_r_stg_min2(idir,inc2)
        ixRmin3=ixR_r_stg_min3(idir,inc3);ixRmax1=ixR_r_stg_max1(idir,inc1)
        ixRmax2=ixR_r_stg_max2(idir,inc2);ixRmax3=ixR_r_stg_max3(idir,inc3);

        ibufnext=ibuf_recv_r+sizes_r_recv_stg(idir,inc1,inc2,inc3)
        pole_buf_stg%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,idir)&
           =reshape(source=recvbuffer_r(ibuf_recv_r:ibufnext-1),&
          shape=shape(pstate(igrid)%ws%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
             ixRmin3:ixRmax3,idir)))
        ibuf_recv_r=ibufnext

        call pole_copy_stg(pstate(igrid)%ws,ixRmin1,ixRmin2,ixRmin3,ixRmax1,&
           ixRmax2,ixRmax3,pole_buf_stg,ixRmin1,ixRmin2,ixRmin3,ixRmax1,&
           ixRmax2,ixRmax3,ipole,i1,i2,i3,idir)

     end do
     


   end if
end do
end do
end do

end if !! ipole == 0


end subroutine bc_fill_r
end subroutine fill_gc_r
!===============================================================================
subroutine gc_restrict(igrid)
! This subroutine fills the coarse representation of the given block.
! 
use mod_comm_gc
use mod_amrvacdef
! ibuf_send and isend_srl were already initialized in init_comm_gc_srl,
! here they are just updated 
integer, intent(in) :: igrid
integer :: i1,i2,i3, ineighbor, ipe_neighbor, my_neighbor_type, ixCoRmin1,&
   ixCoRmin2,ixCoRmin3,ixCoRmax1,ixCoRmax2,ixCoRmax3, ixRmin1,ixRmin2,ixRmin3,&
   ixRmax1,ixRmax2,ixRmax3
!-------------------------------------------------------------------------------

associate(x=>pstate(igrid)%x%x,xCo=>pstateCo(igrid)%x%x,pgeoFi&
   =>pstate(igrid)%geo,pgeoCo=>pstateCo(igrid)%geo)

call set_tmpGlobals(igrid)

if (any(neighbor_type(:,:,:,igrid).eq. 2)) then

   do i3=-1,1
   do i2=-1,1
   do i1=-1,1
      my_neighbor_type=neighbor_type(i1,i2,i3,igrid)      
      ! Restriction is necessary in the physical extent
      ! for sending and in srl ghost regions for later
      ! prolongation
      if (my_neighbor_type .eq. 0 .or. my_neighbor_type .eq. 3) then
      
         ixRmin1=ixR_srl_min1(i1);ixRmin2=ixR_srl_min2(i2)
         ixRmin3=ixR_srl_min3(i3);ixRmax1=ixR_srl_max1(i1)
         ixRmax2=ixR_srl_max2(i2);ixRmax3=ixR_srl_max3(i3);
         ixCoRmin1=int((ixRmin1+dixB+1)/2);ixCoRmin2=int((ixRmin2+dixB+1)/2)
         ixCoRmin3=int((ixRmin3+dixB+1)/2);
         ixCoRmax1=int((ixRmax1+dixB+1)/2);ixCoRmax2=int((ixRmax2+dixB+1)/2)
         ixCoRmax3=int((ixRmax3+dixB+1)/2);
         call coarsen_grid(pstate(igrid),x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
            ixGhi3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,&
            pstateCo(igrid),xCo, ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,&
            ixCoGmax2,ixCoGmax3,ixCoRmin1,ixCoRmin2,ixCoRmin3,ixCoRmax1,&
            ixCoRmax2,ixCoRmax3,pgeoFi,pgeoCo, coarsenprimitive,.true.)

      end if
   end do
   end do
   end do
end if

end associate

end subroutine gc_restrict
!===============================================================================
subroutine gc_prolong(igrid)
use mod_comm_gc
use mod_amrvacdef

integer, intent(in) :: igrid
! ... local ...
integer :: i1,i2,i3,iside,idims,my_neighbor_type

logical,dimension(-1:1,-1:1,-1:1)   :: NeedProlong
integer                      :: idim

!-------------------------------------------------------------------------------
call set_tmpGlobals(igrid)



if (any(neighbor_type(:,:,:,igrid)==2)) then
   ! Fill the array NeedProlong
   ! and prolong cell-centered variables
   NeedProlong=.false.

   do i3=-1,1
do i2=-1,1
do i1=-1,1
      if (i1==0.and.i2==0.and.i3==0) cycle
      my_neighbor_type=neighbor_type(i1,i2,i3,igrid)
      if (my_neighbor_type==2) then
         call bc_prolong
         NeedProlong(i1,i2,i3)=.true.
      end if
   end do
end do
end do


   ! Ghost cell prolongation for staggered variables
   ! must be done in a specific order.
   ! First the first neighbours, which have 2 indices=0 in 3D
   ! or one index=0 in 2D
   do idim=1,3
     i1=0;i2=0;i3=0;
     select case(idim)
    case(1)
       do i1=-1,1,2
         if (NeedProlong(i1,i2,i3)) call bc_prolong_stg(NeedProlong)
       end do
     
case(2)
       do i2=-1,1,2
         if (NeedProlong(i1,i2,i3)) call bc_prolong_stg(NeedProlong)
       end do
     
case(3)
       do i3=-1,1,2
         if (NeedProlong(i1,i2,i3)) call bc_prolong_stg(NeedProlong)
       end do
     
     end select
   end do
   ! Then the second neighbours which have 1 index=0 in 3D
   ! (Only in 3D)
   
   i1=0;
   do i2=-1,1,2
   do i3=-1,1,2
   if (NeedProlong(i1,i2,i3)) call bc_prolong_stg(NeedProlong)
   end do
   end do
   i2=0;
   do i3=-1,1,2
   do i1=-1,1,2
   if (NeedProlong(i1,i2,i3)) call bc_prolong_stg(NeedProlong)
   end do
   end do
   i3=0;
   do i1=-1,1,2
   do i2=-1,1,2
   if (NeedProlong(i1,i2,i3)) call bc_prolong_stg(NeedProlong)
   end do
   end do
  
   ! Finally, the corners, that have no index=0
  do i1=-1,1,2
do i2=-1,1,2
do i3=-1,1,2
   if (NeedProlong(i1,i2,i3)) call bc_prolong_stg(NeedProlong)
 end do
end do
end do

end if



contains
!=============================================================================
! internal procedures
!=============================================================================
subroutine bc_prolong
integer :: ixFimin1,ixFimin2,ixFimin3,ixFimax1,ixFimax2,ixFimax3,ixComin1,&
   ixComin2,ixComin3,ixComax1,ixComax2,ixComax3,ii1,ii2,ii3
double precision :: dxFi1,dxFi2,dxFi3, dxCo1,dxCo2,dxCo3, xFimin1,xFimin2,&
   xFimin3, xComin1,xComin2,xComin3, invdxCo1,invdxCo2,invdxCo3
integer :: ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3
logical :: skip_primitive(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
!--- for debugging ----
integer :: idx1,idx2,idx3
!-----------------------------------------------------------------------------
associate(pgeoCo=>pstateCo(igrid)%geo)


ixFimin1=ixR_srl_min1(i1);ixFimin2=ixR_srl_min2(i2);ixFimin3=ixR_srl_min3(i3)
ixFimax1=ixR_srl_max1(i1);ixFimax2=ixR_srl_max2(i2);ixFimax3=ixR_srl_max3(i3);

dxFi1=rnode(rpdx1_,igrid);dxFi2=rnode(rpdx2_,igrid);dxFi3=rnode(rpdx3_,igrid);
dxCo1=two*dxFi1;dxCo2=two*dxFi2;dxCo3=two*dxFi3;
invdxCo1=1.d0/dxCo1;invdxCo2=1.d0/dxCo2;invdxCo3=1.d0/dxCo3;

xFimin1=rnode(rpxmin1_,igrid)-dble(dixB)*dxFi1
xFimin2=rnode(rpxmin2_,igrid)-dble(dixB)*dxFi2
xFimin3=rnode(rpxmin3_,igrid)-dble(dixB)*dxFi3;
xComin1=rnode(rpxmin1_,igrid)-dble(dixB)*dxCo1
xComin2=rnode(rpxmin2_,igrid)-dble(dixB)*dxCo2
xComin3=rnode(rpxmin3_,igrid)-dble(dixB)*dxCo3;

! moved the physical boundary filling here, to only fill the
! part needed

ixComin1=int((xFimin1+(dble(ixFimin1)-half)*dxFi1-xComin1)*invdxCo1)+1-1
ixComin2=int((xFimin2+(dble(ixFimin2)-half)*dxFi2-xComin2)*invdxCo2)+1-1
ixComin3=int((xFimin3+(dble(ixFimin3)-half)*dxFi3-xComin3)*invdxCo3)+1-1;
ixComax1=int((xFimin1+(dble(ixFimax1)-half)*dxFi1-xComin1)*invdxCo1)+1+1
ixComax2=int((xFimin2+(dble(ixFimax2)-half)*dxFi2-xComin2)*invdxCo2)+1+1
ixComax3=int((xFimin3+(dble(ixFimax3)-half)*dxFi3-xComin3)*invdxCo3)+1+1;


do idims=1,ndim
   do iside=1,2
      ii1=kr(1,idims)*(2*iside-3);ii2=kr(2,idims)*(2*iside-3)
      ii3=kr(3,idims)*(2*iside-3);

      if (neighbor_type(ii1,ii2,ii3,igrid)/=1) cycle

      if  (( (iside==1.and.idims==1.and.ixComin1<ixCoGmin1+dixB)&
         .or.(iside==1.and.idims==2.and.ixComin2<ixCoGmin2+dixB)&
         .or.(iside==1.and.idims==3.and.ixComin3<ixCoGmin3+dixB) ) &
         .or.( (iside==2.and.idims==1.and.ixComax1>ixCoGmax1-dixB)&
         .or. (iside==2.and.idims==2.and.ixComax2>ixCoGmax2-dixB)&
         .or. (iside==2.and.idims==3.and.ixComax3>ixCoGmax3-dixB)))then
        ixBmin1=merge(ixCoGmin1,ixComin1,idims==1)
        ixBmin2=merge(ixCoGmin2,ixComin2,idims==2)
        ixBmin3=merge(ixCoGmin3,ixComin3,idims==3);
        ixBmax1=merge(ixCoGmax1,ixComax1,idims==1)
        ixBmax2=merge(ixCoGmax2,ixComax2,idims==2)
        ixBmax3=merge(ixCoGmax3,ixComax3,idims==3);
        if(.not.slab)mygeo=>pgeoCo
        if(covariant)myM => mygeo%m

        call bc_phys(iside,idims,time,pstateCo(igrid),ixBmin1,ixBmin2,ixBmin3,&
           ixBmax1,ixBmax2,ixBmax3)
      end if
   end do
end do

if(.not.slab)mygeo=>pgeoCo
if(covariant)myM => mygeo%m
if (amrentropy) then
   call e_to_rhos(ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3,&
      ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,ixComax3,&
      pstateCo(igrid)%w%w,pxCoarse(igrid)%x)
else if (prolongprimitive) then
   ! Convert to primitives only in the region with meaningful values.
   skip_primitive=.true.
  do idx3=ixComin3,ixComax3
    do idx2=ixComin2,ixComax2
    do idx1=ixComin1,ixComax1


      if ((((ixComin1+2-interpolation_order.lt.idx1).and.&
         (idx1.lt.ixComax1-2+interpolation_order)).and.&
           ((ixComin2+2-interpolation_order.lt.idx2).and.&
              (idx2.lt.ixComax2-2+interpolation_order))).or.&
          (((ixComin2+2-interpolation_order.lt.idx2).and.&
             (idx2.lt.ixComax2-2+interpolation_order)).and.&
           ((ixComin3+2-interpolation_order.lt.idx3).and.&
              (idx3.lt.ixComax3-2+interpolation_order))).or.&
          (((ixComin3+2-interpolation_order.lt.idx3).and.&
             (idx3.lt.ixComax3-2+interpolation_order)).and.&
           ((ixComin1+2-interpolation_order.lt.idx1).and.&
              (idx1.lt.ixComax1-2+interpolation_order)))) then

      skip_primitive(idx1,idx2,idx3)=.false.
      call primitive(ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,&
         ixCoGmax3,idx1,idx2,idx3,idx1,idx2,idx3,pstateCo(igrid)%w%w,&
         pxCoarse(igrid)%x)

end if

  end do
    end do
    end do
end if

select case (typeghostfill)
case ("linear")
   call interpolation_linear(pstate(igrid),ixFimin1,ixFimin2,ixFimin3,&
      ixFimax1,ixFimax2,ixFimax3,dxFi1,dxFi2,dxFi3,xFimin1,xFimin2,xFimin3,&
       pstateCo(igrid),ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,ixComax3,&
      dxCo1,dxCo2,dxCo3,invdxCo1,invdxCo2,invdxCo3,xComin1,xComin2,xComin3)
case ("copy")
   call interpolation_copy(pstate(igrid)%w,ixFimin1,ixFimin2,ixFimin3,&
      ixFimax1,ixFimax2,ixFimax3,dxFi1,dxFi2,dxFi3,xFimin1,xFimin2,xFimin3,&
       pstateCo(igrid)%w,dxCo1,dxCo2,dxCo3,invdxCo1,invdxCo2,invdxCo3,xComin1,&
      xComin2,xComin3)
case ("unlimit")
   call interpolation_unlimit(pstate(igrid)%w,ixFimin1,ixFimin2,ixFimin3,&
      ixFimax1,ixFimax2,ixFimax3,dxFi1,dxFi2,dxFi3,xFimin1,xFimin2,xFimin3,&
       pstateCo(igrid)%w,dxCo1,dxCo2,dxCo3,invdxCo1,invdxCo2,invdxCo3,xComin1,&
      xComin2,xComin3)
case default
   write (unitterm,*) "Undefined typeghostfill ",typeghostfill
   call mpistop("")
end select

if(.not.slab)mygeo=>pgeoCo
if(covariant)myM => mygeo%m
if (amrentropy) then
    call rhos_to_e(ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,&
       ixCoGmax3,ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,ixComax3,&
       pstateCo(igrid)%w%w,pxCoarse(igrid)%x)
else if (prolongprimitive) then
    call conserve(ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3,&
       ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,ixComax3,&
       pstateCo(igrid)%w%w,pxCoarse(igrid)%x,skip_primitive)
end if

end associate
end subroutine bc_prolong

!=============================================================================
subroutine bc_prolong_stg(NeedProlong)
use mod_amr_fct
logical,dimension(-1:1,-1:1,-1:1) :: NeedProlong
!... local ...
logical                    :: fine_min1in,fine_min2in,fine_min3in,fine_max1in,&
   fine_max2in,fine_max3in
integer                    :: ixFimin1,ixFimin2,ixFimin3,ixFimax1,ixFimax2,&
   ixFimax3,ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,ixComax3
double precision           :: dxFi1,dxFi2,dxFi3,dxCo1,dxCo2,dxCo3,xFimin1,&
   xFimin2,xFimin3,xComin1,xComin2,xComin3,invdxCo1,invdxCo2,invdxCo3
!-----------------------------------------------------------------------------
! Check what is already at the desired level
fine_min1in=.false.;fine_min2in=.false.;fine_min3in=.false.
fine_max1in=.false.;fine_max2in=.false.;fine_max3in=.false.;

if(i1.gt.-1) fine_min1in=(.not.NeedProlong(i1-kr(1,1),i2-kr(1,2),i3-kr(1,&
   3)).and.neighbor_type(i1-kr(1,1),i2-kr(1,2),i3-kr(1,3),igrid)/=1)
if(i1.lt.1)  fine_max1in=(.not.NeedProlong(i1+kr(1,1),i2+kr(1,2),i3+kr(1,&
   3)).and.neighbor_type(i1+kr(1,1),i2+kr(1,2),i3+kr(1,3),igrid)/=1)


if(i2.gt.-1) fine_min2in=(.not.NeedProlong(i1-kr(2,1),i2-kr(2,2),i3-kr(2,&
   3)).and.neighbor_type(i1-kr(2,1),i2-kr(2,2),i3-kr(2,3),igrid)/=1)
if(i2.lt.1)  fine_max2in=(.not.NeedProlong(i1+kr(2,1),i2+kr(2,2),i3+kr(2,&
   3)).and.neighbor_type(i1+kr(2,1),i2+kr(2,2),i3+kr(2,3),igrid)/=1)


if(i3.gt.-1) fine_min3in=(.not.NeedProlong(i1-kr(3,1),i2-kr(3,2),i3-kr(3,&
   3)).and.neighbor_type(i1-kr(3,1),i2-kr(3,2),i3-kr(3,3),igrid)/=1)
if(i3.lt.1)  fine_max3in=(.not.NeedProlong(i1+kr(3,1),i2+kr(3,2),i3+kr(3,&
   3)).and.neighbor_type(i1+kr(3,1),i2+kr(3,2),i3+kr(3,3),igrid)/=1)


ixFimin1=ixR_srl_min1(i1);ixFimin2=ixR_srl_min2(i2);ixFimin3=ixR_srl_min3(i3)
ixFimax1=ixR_srl_max1(i1);ixFimax2=ixR_srl_max2(i2);ixFimax3=ixR_srl_max3(i3);

dxFi1=rnode(rpdx1_,igrid);dxFi2=rnode(rpdx2_,igrid);dxFi3=rnode(rpdx3_,igrid);
dxCo1=two*dxFi1;dxCo2=two*dxFi2;dxCo3=two*dxFi3;
invdxCo1=1.d0/dxCo1;invdxCo2=1.d0/dxCo2;invdxCo3=1.d0/dxCo3;

xFimin1=rnode(rpxmin1_,igrid)-dble(dixB)*dxFi1
xFimin2=rnode(rpxmin2_,igrid)-dble(dixB)*dxFi2
xFimin3=rnode(rpxmin3_,igrid)-dble(dixB)*dxFi3;
xComin1=rnode(rpxmin1_,igrid)-dble(dixB)*dxCo1
xComin2=rnode(rpxmin2_,igrid)-dble(dixB)*dxCo2
xComin3=rnode(rpxmin3_,igrid)-dble(dixB)*dxCo3;

! moved the physical boundary filling here, to only fill the
! part needed

ixComin1=int((xFimin1+(dble(ixFimin1)-half)*dxFi1-xComin1)*invdxCo1)+1-1
ixComin2=int((xFimin2+(dble(ixFimin2)-half)*dxFi2-xComin2)*invdxCo2)+1-1
ixComin3=int((xFimin3+(dble(ixFimin3)-half)*dxFi3-xComin3)*invdxCo3)+1-1;
ixComax1=int((xFimin1+(dble(ixFimax1)-half)*dxFi1-xComin1)*invdxCo1)+1+1
ixComax2=int((xFimin2+(dble(ixFimax2)-half)*dxFi2-xComin2)*invdxCo2)+1+1
ixComax3=int((xFimin3+(dble(ixFimax3)-half)*dxFi3-xComin3)*invdxCo3)+1+1;

call prolong_2nd_stg(pstateCo(igrid),pstate(igrid),ixComin1,ixComin2,ixComin3,&
   ixComax1,ixComax2,ixComax3,ixFimin1,ixFimin2,ixFimin3,ixFimax1,ixFimax2,&
   ixFimax3,dxCo1,dxCo2,dxCo3,xComin1,xComin2,xComin3,dxFi1,dxFi2,dxFi3,&
   xFimin1,xFimin2,xFimin3,.true.,fine_min1in,fine_min2in,fine_min3in,&
   fine_max1in,fine_max2in,fine_max3in)

 !The current region has already been refined, so it doesn t need to be prolonged again

 NeedProlong(i1,i2,i3)=.false. 

end subroutine bc_prolong_stg

!=============================================================================
!=== Prolongation routines (Should be made independent and moved elsewhere) ==
!=============================================================================
subroutine interpolation_linear(psFi,ixFimin1,ixFimin2,ixFimin3,ixFimax1,&
   ixFimax2,ixFimax3,dxFi1,dxFi2,dxFi3,xFimin1,xFimin2,xFimin3, psCo,ixComin1,&
   ixComin2,ixComin3,ixComax1,ixComax2,ixComax3,dxCo1,dxCo2,dxCo3,invdxCo1,&
   invdxCo2,invdxCo3,xComin1,xComin2,xComin3)

integer, intent(in) :: ixFimin1,ixFimin2,ixFimin3,ixFimax1,ixFimax2,ixFimax3
double precision, intent(in) :: dxFi1,dxFi2,dxFi3, xFimin1,xFimin2,xFimin3,&
   dxCo1,dxCo2,dxCo3, invdxCo1,invdxCo2,invdxCo3, xComin1,xComin2,xComin3
type(state) :: psCo, psFi
type(walloc) :: pwCo, pwFi
integer :: ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,ixComax3,ixCo1,ixCo2,&
   ixCo3, jxCo1,jxCo2,jxCo3, hxCo1,hxCo2,hxCo3, ixFi1,ixFi2,ixFi3, ix1,ix2,&
   ix3, iw, idims
double precision :: xCo1,xCo2,xCo3, xFi1,xFi2,xFi3, eta1,eta2,eta3
double precision :: slopeL, slopeR, slopeC, signC, signR, slope(nwflux+nwaux,&
   ndim)
!-----------------------------------------------------------------------------
associate(pgeoFi=>pstate(igrid)%geo)

pwCo=psCo%w
pwFi=psFi%w
do ixFi3 = ixFimin3,ixFimax3
   ! cell-centered coordinates of fine grid point
   xFi3=xFimin3+(dble(ixFi3)-half)*dxFi3

   ! indices of coarse cell which contains the fine cell
   ixCo3=int((xFi3-xComin3)*invdxCo3)+1

   ! cell-centered coordinate for coarse cell
   xCo3=xComin3+(dble(ixCo3)-half)*dxCo3
do ixFi2 = ixFimin2,ixFimax2
   ! cell-centered coordinates of fine grid point
   xFi2=xFimin2+(dble(ixFi2)-half)*dxFi2

   ! indices of coarse cell which contains the fine cell
   ixCo2=int((xFi2-xComin2)*invdxCo2)+1

   ! cell-centered coordinate for coarse cell
   xCo2=xComin2+(dble(ixCo2)-half)*dxCo2
do ixFi1 = ixFimin1,ixFimax1
   ! cell-centered coordinates of fine grid point
   xFi1=xFimin1+(dble(ixFi1)-half)*dxFi1

   ! indices of coarse cell which contains the fine cell
   ixCo1=int((xFi1-xComin1)*invdxCo1)+1

   ! cell-centered coordinate for coarse cell
   xCo1=xComin1+(dble(ixCo1)-half)*dxCo1

   ! normalized distance between fine/coarse cell center
   ! in coarse cell: ranges from -0.5 to 0.5 in each direction
   ! (origin is coarse cell center)
   if (slab) then
      eta1=(xFi1-xCo1)*invdxCo1;eta2=(xFi2-xCo2)*invdxCo2
      eta3=(xFi3-xCo3)*invdxCo3;
   else
      ix1=2*int((ixFi1+ixMlo1)/2)-ixMlo1;ix2=2*int((ixFi2+ixMlo2)/2)-ixMlo2
      ix3=2*int((ixFi3+ixMlo3)/2)-ixMlo3;
      eta1=(xFi1-xCo1)*invdxCo1 *two*(one-pgeoFi%dvolume(ixFi1,ixFi2,ixFi3) &
         /sum(pgeoFi%dvolume(ix1:ix1+1,ixFi2,ixFi3))) 
      eta2=(xFi2-xCo2)*invdxCo2 *two*(one-pgeoFi%dvolume(ixFi1,ixFi2,ixFi3) &
         /sum(pgeoFi%dvolume(ixFi1,ix2:ix2+1,ixFi3))) 
      eta3=(xFi3-xCo3)*invdxCo3 *two*(one-pgeoFi%dvolume(ixFi1,ixFi2,ixFi3) &
         /sum(pgeoFi%dvolume(ixFi1,ixFi2,ix3:ix3+1))) 
   end if

   
   do idims=1,ndim
      hxCo1=ixCo1-kr(1,idims)
      hxCo2=ixCo2-kr(2,idims)
      hxCo3=ixCo3-kr(3,idims)
      jxCo1=ixCo1+kr(1,idims)
      jxCo2=ixCo2+kr(2,idims)
      jxCo3=ixCo3+kr(3,idims)

      do iw=1,nwflux+nwaux
         slopeL=pwCo%w(ixCo1,ixCo2,ixCo3,iw)-pwCo%w(hxCo1,hxCo2,hxCo3,iw)
         slopeR=pwCo%w(jxCo1,jxCo2,jxCo3,iw)-pwCo%w(ixCo1,ixCo2,ixCo3,iw)
         slopeC=half*(slopeR+slopeL)

         ! get limited slope
         signR=sign(one,slopeR)
         signC=sign(one,slopeC)
         select case(typeprolonglimit)
         case('minmod')
           slope(iw,idims)=signR*max(zero,min(dabs(slopeR), signR*slopeL))
         case('woodward')
           slope(iw,idims)=two*signR*max(zero,min(dabs(slopeR), signR*slopeL,&
              signR*half*slopeC))
         case('mcbeta')
           slope(iw,idims)=signR*max(zero,min(mcbeta*dabs(slopeR),&
               mcbeta*signR*slopeL,signR*slopeC))
         case('koren')
           slope(iw,idims)=signR*max(zero,min(two*signR*slopeL,&
               (dabs(slopeR)+two*slopeL*signR)*third,two*dabs(slopeR)))
         case default
           slope(iw,idims)=signC*max(zero,min(dabs(slopeC), signC*slopeL,&
              signC*slopeR))
         end select
      end do
   end do

   ! Interpolate from coarse cell using limited slopes
   pwFi%w(ixFi1,ixFi2,ixFi3,1:nwflux+nwaux)=pwCo%w(ixCo1,ixCo2,ixCo3,&
      1:nwflux+nwaux)+(slope(1:nwflux+nwaux,1)*eta1)+(slope(1:nwflux+nwaux,&
      2)*eta2)+(slope(1:nwflux+nwaux,3)*eta3)

end do
end do
end do

if(.not.slab)mygeo=>pgeoFi
if(covariant)myM => mygeo%m
if (amrentropy) then
   call rhos_to_e(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixFimin1,ixFimin2,&
      ixFimin3,ixFimax1,ixFimax2,ixFimax3,pwFi%w,px(igrid)%x)
else if (prolongprimitive) then
   call conserve(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixFimin1,ixFimin2,&
      ixFimin3,ixFimax1,ixFimax2,ixFimax3,pwFi%w,px(igrid)%x,patchfalse)
end if

end associate
end subroutine interpolation_linear
!=============================================================================
subroutine interpolation_copy(pwFi,ixFimin1,ixFimin2,ixFimin3,ixFimax1,&
   ixFimax2,ixFimax3,dxFi1,dxFi2,dxFi3,xFimin1,xFimin2,xFimin3, pwCo,dxCo1,&
   dxCo2,dxCo3,invdxCo1,invdxCo2,invdxCo3,xComin1,xComin2,xComin3)

integer, intent(in) :: ixFimin1,ixFimin2,ixFimin3,ixFimax1,ixFimax2,ixFimax3
double precision, intent(in) :: dxFi1,dxFi2,dxFi3, xFimin1,xFimin2,xFimin3,&
   dxCo1,dxCo2,dxCo3, invdxCo1,invdxCo2,invdxCo3, xComin1,xComin2,xComin3
type(walloc) :: pwCo, pwFi

integer :: ixCo1,ixCo2,ixCo3, ixFi1,ixFi2,ixFi3
double precision :: xFi1,xFi2,xFi3
!-----------------------------------------------------------------------------
associate(pgeoFi=>pstate(igrid)%geo)

do ixFi3 = ixFimin3,ixFimax3
   ! cell-centered coordinates of fine grid point
   xFi3=xFimin3+(dble(ixFi3)-half)*dxFi3

   ! indices of coarse cell which contains the fine cell
   ixCo3=int((xFi3-xComin3)*invdxCo3)+1
do ixFi2 = ixFimin2,ixFimax2
   ! cell-centered coordinates of fine grid point
   xFi2=xFimin2+(dble(ixFi2)-half)*dxFi2

   ! indices of coarse cell which contains the fine cell
   ixCo2=int((xFi2-xComin2)*invdxCo2)+1
do ixFi1 = ixFimin1,ixFimax1
   ! cell-centered coordinates of fine grid point
   xFi1=xFimin1+(dble(ixFi1)-half)*dxFi1

   ! indices of coarse cell which contains the fine cell
   ixCo1=int((xFi1-xComin1)*invdxCo1)+1

   ! Copy from coarse cell
   pwFi%w(ixFi1,ixFi2,ixFi3,1:nwflux+nwaux)=pwCo%w(ixCo1,ixCo2,ixCo3,&
      1:nwflux+nwaux)

end do
end do
end do



if(.not.slab)mygeo=>pgeoFi
if(covariant)myM => mygeo%m
if (amrentropy) then
   call rhos_to_e(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixFimin1,ixFimin2,&
      ixFimin3,ixFimax1,ixFimax2,ixFimax3,pwFi%w,px(igrid)%x)
else if (prolongprimitive) then
   call conserve(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixFimin1,ixFimin2,&
      ixFimin3,ixFimax1,ixFimax2,ixFimax3,pwFi%w,px(igrid)%x,patchfalse)
end if

end associate
end subroutine interpolation_copy
!=============================================================================
subroutine interpolation_unlimit(pwFi,ixFimin1,ixFimin2,ixFimin3,ixFimax1,&
   ixFimax2,ixFimax3,dxFi1,dxFi2,dxFi3,xFimin1,xFimin2,xFimin3, pwCo,dxCo1,&
   dxCo2,dxCo3,invdxCo1,invdxCo2,invdxCo3,xComin1,xComin2,xComin3)

integer, intent(in) :: ixFimin1,ixFimin2,ixFimin3,ixFimax1,ixFimax2,ixFimax3
double precision, intent(in) :: dxFi1,dxFi2,dxFi3, xFimin1,xFimin2,xFimin3,&
    dxCo1,dxCo2,dxCo3,invdxCo1,invdxCo2,invdxCo3, xComin1,xComin2,xComin3
type(walloc) :: pwCo, pwFi

integer :: ixCo1,ixCo2,ixCo3, jxCo1,jxCo2,jxCo3, hxCo1,hxCo2,hxCo3, ixFi1,&
   ixFi2,ixFi3, ix1,ix2,ix3, idims
double precision :: xCo1,xCo2,xCo3, xFi1,xFi2,xFi3, eta1,eta2,eta3
double precision :: slope(nwflux+nwaux,ndim)
!-----------------------------------------------------------------------------
associate(pgeoFi=>pstate(igrid)%geo)

do ixFi3 = ixFimin3,ixFimax3
   ! cell-centered coordinates of fine grid point
   xFi3=xFimin3+(dble(ixFi3)-half)*dxFi3

   ! indices of coarse cell which contains the fine cell
   ixCo3=int((xFi3-xComin3)*invdxCo3)+1

   ! cell-centered coordinate for coarse cell
   xCo3=xComin3+(dble(ixCo3)-half)*dxCo3
do ixFi2 = ixFimin2,ixFimax2
   ! cell-centered coordinates of fine grid point
   xFi2=xFimin2+(dble(ixFi2)-half)*dxFi2

   ! indices of coarse cell which contains the fine cell
   ixCo2=int((xFi2-xComin2)*invdxCo2)+1

   ! cell-centered coordinate for coarse cell
   xCo2=xComin2+(dble(ixCo2)-half)*dxCo2
do ixFi1 = ixFimin1,ixFimax1
   ! cell-centered coordinates of fine grid point
   xFi1=xFimin1+(dble(ixFi1)-half)*dxFi1

   ! indices of coarse cell which contains the fine cell
   ixCo1=int((xFi1-xComin1)*invdxCo1)+1

   ! cell-centered coordinate for coarse cell
   xCo1=xComin1+(dble(ixCo1)-half)*dxCo1

   ! normalized distance between fine/coarse cell center
   ! in coarse cell: ranges from -0.5 to 0.5 in each direction
   ! (origin is coarse cell center)
   if (slab) then
      eta1=(xFi1-xCo1)*invdxCo1;eta2=(xFi2-xCo2)*invdxCo2
      eta3=(xFi3-xCo3)*invdxCo3;
   else
      ix1=2*int((ixFi1+ixMlo1)/2)-ixMlo1;ix2=2*int((ixFi2+ixMlo2)/2)-ixMlo2
      ix3=2*int((ixFi3+ixMlo3)/2)-ixMlo3;
      eta1=(xFi1-xCo1)*invdxCo1 *two*(one-pgeoFi%dvolume(ixFi1,ixFi2,ixFi3) &
         /sum(pgeoFi%dvolume(ix1:ix1+1,ixFi2,ixFi3))) 
      eta2=(xFi2-xCo2)*invdxCo2 *two*(one-pgeoFi%dvolume(ixFi1,ixFi2,ixFi3) &
         /sum(pgeoFi%dvolume(ixFi1,ix2:ix2+1,ixFi3))) 
      eta3=(xFi3-xCo3)*invdxCo3 *two*(one-pgeoFi%dvolume(ixFi1,ixFi2,ixFi3) &
         /sum(pgeoFi%dvolume(ixFi1,ixFi2,ix3:ix3+1))) 
   end if

   do idims=1,ndim
      hxCo1=ixCo1-kr(1,idims)
      hxCo2=ixCo2-kr(2,idims)
      hxCo3=ixCo3-kr(3,idims)
      jxCo1=ixCo1+kr(1,idims)
      jxCo2=ixCo2+kr(2,idims)
      jxCo3=ixCo3+kr(3,idims)

      ! get centered slope
      slope(1:nwflux+nwaux,idims)=half*(pwCo%w(jxCo1,jxCo2,jxCo3,&
         1:nwflux+nwaux)-pwCo%w(hxCo1,hxCo2,hxCo3,1:nwflux+nwaux))
   end do

   ! Interpolate from coarse cell using centered slopes
   pwFi%w(ixFi1,ixFi2,ixFi3,1:nwflux+nwaux)=pwCo%w(ixCo1,ixCo2,ixCo3,&
      1:nwflux+nwaux)+(slope(1:nwflux+nwaux,1)*eta1)+(slope(1:nwflux+nwaux,&
      2)*eta2)+(slope(1:nwflux+nwaux,3)*eta3)
end do
end do
end do



if(.not.slab)mygeo=>pgeoFi
if(covariant)myM => mygeo%m
if (amrentropy) then
   call rhos_to_e(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixFimin1,ixFimin2,&
      ixFimin3,ixFimax1,ixFimax2,ixFimax3,pwFi%w,px(igrid)%x)
else if (prolongprimitive) then
   call conserve(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixFimin1,ixFimin2,&
      ixFimin3,ixFimax1,ixFimax2,ixFimax3,pwFi%w,px(igrid)%x,patchfalse)
end if

end associate
end subroutine interpolation_unlimit
!=============================================================================


end subroutine gc_prolong
!===============================================================================
subroutine send_gc_r(igrid)
! Before this routine is called,
! the coarse representation of the block must be first filled
use mod_comm_gc
use mod_amrvacdef
! ibuf_send and isend_srl were already initialized in init_comm_gc_srl,
! here they are just updated 
integer, intent(in) :: igrid
integer :: i1,i2,i3, ineighbor, ipe_neighbor, my_neighbor_type
integer :: ipole
!-------------------------------------------------------------------------------
call set_tmpGlobals(igrid)

do i3=-1,1
do i2=-1,1
do i1=-1,1
   if (i1==0.and.i2==0.and.i3==0) cycle
   my_neighbor_type=neighbor_type(i1,i2,i3,igrid)
   if (my_neighbor_type==2) call bc_send_restrict
end do
end do
end do

contains
!=============================================================================
! internal procedures
!=============================================================================
subroutine bc_send_restrict
integer :: ic1,ic2,ic3,inc1,inc2,inc3,n_i1,n_i2,n_i3,n_inc1,n_inc2,n_inc3,&
   ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3,idir,ibufaux,ibufnext
 ! Auxialiary array to avoid problems with the preprocessor...
 ! (It keeps splitting the line in the wrong place)
integer, dimension(1) :: sizes
!-----------------------------------------------------------------------------
ic1=1+modulo(node(pig1_,igrid)-1,2);ic2=1+modulo(node(pig2_,igrid)-1,2)
ic3=1+modulo(node(pig3_,igrid)-1,2);
if (.not.(i1==0.or.i1==2*ic1-3).or..not.(i2==0.or.i2==2*ic2-3)&
   .or..not.(i3==0.or.i3==2*ic3-3)) return

ineighbor=neighbor(1,i1,i2,i3,igrid)
ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
ipole=neighbor_pole(i1,i2,i3,igrid)

if (ipe_neighbor.ne.mype) then

  inc1=i1+ic1;inc2=i2+ic2;inc3=i3+ic3;
  ! The ghost region of the neighbor changes
  ! if there is a pole
  select case (ipole)
  case(0) ! No pole
    n_inc1=-2*i1+ic1;n_inc2=-2*i2+ic2;n_inc3=-2*i3+ic3;
 case (1) ! Pole in this direction
    n_inc1=2*i1+(3-ic1);n_inc2=-2*i2+ic2;n_inc3=-2*i3+ic3;
 case (2) ! Pole in this direction
    n_inc1=-2*i1+ic1;n_inc2=2*i2+(3-ic2);n_inc3=-2*i3+ic3;
 case (3) ! Pole in this direction
    n_inc1=-2*i1+ic1;n_inc2=-2*i2+ic2;n_inc3=2*i3+(3-ic3);
  end select

  ! fill corresponding part of the send buffer...
  ixSmin1=ixS_r_min1(i1);ixSmin2=ixS_r_min2(i2);ixSmin3=ixS_r_min3(i3)
  ixSmax1=ixS_r_max1(i1);ixSmax2=ixS_r_max2(i2);ixSmax3=ixS_r_max3(i3);
  ibufaux=ibuf_send_r
  ibufnext=ibufaux+sizes_r_send(i1,i2,i3)

  sizes=(/sizes_r_send(i1,i2,i3)/)

  sendbuffer_r(ibufaux:ibufnext-1)=reshape(pstateCo(igrid)%w%w&
     (ixSmin1:ixSmax1,ixSmin2:ixSmax2,ixSmin3:ixSmax3,1:nwflux+nwaux),sizes)

  ibufaux=ibufnext

  do idir=1,3
    ixSmin1=ixS_r_stg_min1(idir,i1);ixSmin2=ixS_r_stg_min2(idir,i2)
    ixSmin3=ixS_r_stg_min3(idir,i3);ixSmax1=ixS_r_stg_max1(idir,i1)
    ixSmax2=ixS_r_stg_max2(idir,i2);ixSmax3=ixS_r_stg_max3(idir,i3);
    ibufnext=ibufaux+sizes_r_send_stg(idir,i1,i2,i3)
    sizes=(/sizes_r_send_stg(idir,i1,i2,i3)/)
    sendbuffer_r(ibufaux:ibufnext-1)=&
      reshape(pstateCo(igrid)%ws%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
         ixSmin3:ixSmax3,idir),sizes)   
    ibufaux=ibufnext
  end do


  ! ...and send
  itag=(3**3+4**3)*(ineighbor-1)+3**3+n_inc1*4**(1-1)+n_inc2*4&
     **(2-1)+n_inc3*4**(3-1)
  isend_r=isend_r+1
  call MPI_ISEND(sendbuffer_r(ibuf_send_r),sizes_r_send_total(i1,i2,i3),&
      MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,sendrequest_r(isend_r),&
     ierrmpi)

  ibuf_send_r=ibufnext

end if

end subroutine bc_send_restrict
end subroutine send_gc_r
!===============================================================================
subroutine make_task_fill_gc_p(igrid)
use mod_comm_gc
use mod_amrvacdef

! In a loop, place MPI receive requests, updating the part of the receive
! buffer where ghost cells will be written
! Note: ibuf_recv_p and irecv_p were already initialized in init_comm_gc,
! here they are just updated 
integer, intent(in) :: igrid
integer :: i1,i2,i3, ic1,ic2,ic3, inc1,inc2,inc3, ineighbor, ipe_neighbor,&
    my_neighbor_type
!-------------------------------------------------------------------------------

call set_tmpGlobals(igrid)

do i3=-1,1
do i2=-1,1
do i1=-1,1
   if (i1==0.and.i2==0.and.i3==0) cycle
   my_neighbor_type=neighbor_type(i1,i2,i3,igrid)
   if (my_neighbor_type.ne.2) cycle
   ic1=1+modulo(node(pig1_,igrid)-1,2);ic2=1+modulo(node(pig2_,igrid)-1,2)
   ic3=1+modulo(node(pig3_,igrid)-1,2);
   if ((i1==0.or.i1==2*ic1-3).and.(i2==0.or.i2==2*ic2-3).and.(i3==0&
      .or.i3==2*ic3-3)) then
     ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
     if (ipe_neighbor.ne.mype) then
     inc1=ic1+i1;inc2=ic2+i2;inc3=ic3+i3;
     irecv_p=irecv_p+1
     itag=(3**3+4**3)*(igrid-1)+3**3+inc1*4**(1-1)+inc2*4**(2-1)+inc3*4**(3-1)
!     print *, 'recv tag',itag,'mype',mype,'igrid',igrid
       call MPI_IRECV(recvbuffer_p(ibuf_recv_p),sizes_p_recv_total(inc1,inc2,&
          inc3), MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
          recvrequest_p(irecv_p),ierrmpi)
       !call add_task_to_list() 
       ibuf_recv_p=ibuf_recv_p+sizes_p_recv_total(inc1,inc2,inc3)
     end if
   end if
end do
end do
end do

end subroutine make_task_fill_gc_p
!===============================================================================
subroutine fill_gc_p
use mod_comm_gc
use mod_amrvacdef

integer :: i1,i2,i3, igrid, iigrid
integer :: ipole, ineighbor, ipe_neighbor, my_neighbor_type
!-------------------------------------------------------------------------------
! Wait for the receive buffer to be complete

if (nrecv_gc_p>0) then
   call MPI_WAITALL(nrecv_gc_p,recvrequest_p,recvstatus_p,ierrmpi)
   ibuf_recv_p=1
end if

! In a loop with the same order as that in make_task_fill_gc_p,
! directly fill the ghost cells that are in the same cpu and 
! unpack the receive buffer to fill those that are not.

do iigrid=1,igridstail; igrid=igrids(iigrid);
  call set_tmpGlobals(igrid)
     
 do i3=-1,1
 do i2=-1,1
 do i1=-1,1
    if (i1==0.and.i2==0.and.i3==0) cycle
    my_neighbor_type=neighbor_type(i1,i2,i3,igrid)
    if (my_neighbor_type.eq.2) then
      call bc_fill_p
    end if
 end do
 end do
 end do
end do


! Wait for the sends to complete and deallocate communication buffers
if (nrecv_gc_p>0) deallocate(recvbuffer_p,recvstatus_p,recvrequest_p)

if (nsend_gc_p>0) then
   call MPI_WAITALL(nsend_gc_p,sendrequest_p,sendstatus_p,ierrmpi)
   deallocate(sendbuffer_p,sendstatus_p,sendrequest_p)
end if



!! Deallocate pole buffers, now that all communications have ended
deallocate(pole_buf%w)

deallocate(pole_buf_stg%w)



contains
!=============================================================================
! internal procedures
!=============================================================================
subroutine bc_fill_p
integer :: ic1,ic2,ic3,inc1,inc2,inc3,n_inc1,n_inc2,n_inc3,ixSmin1,ixSmin2,&
   ixSmin3,ixSmax1,ixSmax2,ixSmax3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,&
   ixRmax3,idir,ibufnext
!--- for debugging ----
integer :: idx1,idx2,idx3
!-----------------------------------------------------------------------------
ic1=1+modulo(node(pig1_,igrid)-1,2);ic2=1+modulo(node(pig2_,igrid)-1,2)
ic3=1+modulo(node(pig3_,igrid)-1,2);
if (.not.(i1==0.or.i1==2*ic1-3).or..not.(i2==0.or.i2==2*ic2-3)&
   .or..not.(i3==0.or.i3==2*ic3-3)) return

ineighbor=neighbor(1,i1,i2,i3,igrid)
ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
ipole=neighbor_pole(i1,i2,i3,igrid)

if (ipole.eq.0) then   !! There is no pole 

inc1=ic1+i1;inc2=ic2+i2;inc3=ic3+i3;
ixRmin1=ixR_p_min1(inc1);ixRmin2=ixR_p_min2(inc2);ixRmin3=ixR_p_min3(inc3)
ixRmax1=ixR_p_max1(inc1);ixRmax2=ixR_p_max2(inc2);ixRmax3=ixR_p_max3(inc3);
if (ipe_neighbor.eq.mype) then !! Same processor
  n_inc1=-2*i1+ic1;n_inc2=-2*i2+ic2;n_inc3=-2*i3+ic3;
  ixSmin1=ixS_p_min1(n_inc1);ixSmin2=ixS_p_min2(n_inc2)
  ixSmin3=ixS_p_min3(n_inc3);ixSmax1=ixS_p_max1(n_inc1)
  ixSmax2=ixS_p_max2(n_inc2);ixSmax3=ixS_p_max3(n_inc3);
  pstateCo(igrid)%w%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
     1:nwflux+nwaux) =pstate(ineighbor)%w%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
     ixSmin3:ixSmax3,1:nwflux+nwaux)

  
  do idir=1,3
     ixSmin1=ixS_p_stg_min1(idir,n_inc1);ixSmin2=ixS_p_stg_min2(idir,n_inc2)
     ixSmin3=ixS_p_stg_min3(idir,n_inc3);ixSmax1=ixS_p_stg_max1(idir,n_inc1)
     ixSmax2=ixS_p_stg_max2(idir,n_inc2);ixSmax3=ixS_p_stg_max3(idir,n_inc3);
     ixRmin1=ixR_p_stg_min1(idir,inc1);ixRmin2=ixR_p_stg_min2(idir,inc2)
     ixRmin3=ixR_p_stg_min3(idir,inc3);ixRmax1=ixR_p_stg_max1(idir,inc1)
     ixRmax2=ixR_p_stg_max2(idir,inc2);ixRmax3=ixR_p_stg_max3(idir,inc3);
     pstateCo(igrid)%ws%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
        idir)=pstate(ineighbor)%ws%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
        ixSmin3:ixSmax3,idir)
  end do
  
else !! Different processor
  !! Unpack the buffer and fill the ghost cells
  ibufnext=ibuf_recv_p+sizes_p_recv(inc1,inc2,inc3)
  pstateCo(igrid)%w%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
     1:nwflux+nwaux)=reshape(source=recvbuffer_p(ibuf_recv_p:ibufnext-1),&
     shape=shape(pstateCo(igrid)%w%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
     ixRmin3:ixRmax3,1:nwflux+nwaux)))

  ibuf_recv_p=ibufnext

  
  do idir=1,3
    ixRmin1=ixR_p_stg_min1(idir,inc1);ixRmin2=ixR_p_stg_min2(idir,inc2)
    ixRmin3=ixR_p_stg_min3(idir,inc3);ixRmax1=ixR_p_stg_max1(idir,inc1)
    ixRmax2=ixR_p_stg_max2(idir,inc2);ixRmax3=ixR_p_stg_max3(idir,inc3);
 !ixSmin1=ixS_srl_stg_min1(idir,n_i1);ixSmin2=ixS_srl_stg_min2(idir,n_i2);ixSmin3=ixS_srl_stg_min3(idir,n_i3);ixSmax1=ixS_srl_stg_max1(idir,n_i1);ixSmax2=ixS_srl_stg_max2(idir,n_i2);ixSmax3=ixS_srl_stg_max3(idir,n_i3);

    ibufnext=ibuf_recv_p+sizes_p_recv_stg(idir,inc1,inc2,inc3)

    pstateCo(igrid)%ws%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
       idir)=reshape(source=recvbuffer_p(ibuf_recv_p:ibufnext-1),shape&
       =shape(pstateCo(igrid)%ws%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
       ixRmin3:ixRmax3,idir)))
  
    ibuf_recv_p=ibufnext
  end do
  
end if

else !! There is a pole
inc1=ic1+i1;inc2=ic2+i2;inc3=ic3+i3;
!n_inc^D=inc^D; !! (Hope this is correct...)
select case (ipole)
case (1)
   n_inc1=2*i1+(3-ic1);n_inc2=-2*i2+ic2;n_inc3=-2*i3+ic3;
case (2)
   n_inc1=-2*i1+ic1;n_inc2=2*i2+(3-ic2);n_inc3=-2*i3+ic3;
case (3)
   n_inc1=-2*i1+ic1;n_inc2=-2*i2+ic2;n_inc3=2*i3+(3-ic3);
end select

!print *, 'n_inc',n_inc^D,'<--inc',inc^D


if (ipe_neighbor.eq.mype) then
  
  ixSmin1=ixS_p_min1(n_inc1);ixSmin2=ixS_p_min2(n_inc2)
  ixSmin3=ixS_p_min3(n_inc3);ixSmax1=ixS_p_max1(n_inc1)
  ixSmax2=ixS_p_max2(n_inc2);ixSmax3=ixS_p_max3(n_inc3);
  ixRmin1=ixR_p_min1(inc1);ixRmin2=ixR_p_min2(inc2);ixRmin3=ixR_p_min3(inc3)
  ixRmax1=ixR_p_max1(inc1);ixRmax2=ixR_p_max2(inc2);ixRmax3=ixR_p_max3(inc3);
  
  call pole_copy(pstateCo(igrid)%w,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,&
     ixRmax3,pstate(ineighbor)%w,ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,&
     ixSmax3,ipole,i1,i2,i3)
  
  
  do idir=1,3
     ixSmin1=ixS_p_stg_min1(idir,n_inc1);ixSmin2=ixS_p_stg_min2(idir,n_inc2)
     ixSmin3=ixS_p_stg_min3(idir,n_inc3);ixSmax1=ixS_p_stg_max1(idir,n_inc1)
     ixSmax2=ixS_p_stg_max2(idir,n_inc2);ixSmax3=ixS_p_stg_max3(idir,n_inc3);
     ixRmin1=ixR_p_stg_min1(idir,inc1);ixRmin2=ixR_p_stg_min2(idir,inc2)
     ixRmin3=ixR_p_stg_min3(idir,inc3);ixRmax1=ixR_p_stg_max1(idir,inc1)
     ixRmax2=ixR_p_stg_max2(idir,inc2);ixRmax3=ixR_p_stg_max3(idir,inc3);
  
     call pole_copy_stg(pstateCo(igrid)%ws,ixRmin1,ixRmin2,ixRmin3,ixRmax1,&
        ixRmax2,ixRmax3,pstate(ineighbor)%ws,ixSmin1,ixSmin2,ixSmin3,ixSmax1,&
        ixSmax2,ixSmax3,ipole,i1,i2,i3,idir)
  
  end do
  

else
  ixRmin1=ixR_p_min1(inc1);ixRmin2=ixR_p_min2(inc2);ixRmin3=ixR_p_min3(inc3)
  ixRmax1=ixR_p_max1(inc1);ixRmax2=ixR_p_max2(inc2);ixRmax3=ixR_p_max3(inc3);
  !! Unpack the buffer and fill an auxiliary array
  ibufnext=ibuf_recv_p+sizes_p_recv(inc1,inc2,inc3)
  pole_buf%w=zero
  pole_buf%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,1:nwflux+nwaux)&
     =reshape(source=recvbuffer_p(ibuf_recv_p:ibufnext-1),shape&
     =shape(pstateCo(igrid)%w%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
     ixRmin3:ixRmax3,1:nwflux+nwaux)))
  ibuf_recv_p=ibufnext
 
  !! Fill ghost cells
  call pole_copy(pstateCo(igrid)%w,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,&
     ixRmax3,pole_buf,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,ipole,&
     i1,i2,i3)

  
  pole_buf_stg%w=zero
  do idir=1,3
     ixRmin1=ixR_p_stg_min1(idir,inc1);ixRmin2=ixR_p_stg_min2(idir,inc2)
     ixRmin3=ixR_p_stg_min3(idir,inc3);ixRmax1=ixR_p_stg_max1(idir,inc1)
     ixRmax2=ixR_p_stg_max2(idir,inc2);ixRmax3=ixR_p_stg_max3(idir,inc3);

     ibufnext=ibuf_recv_p+sizes_p_recv_stg(idir,inc1,inc2,inc3)
     pole_buf_stg%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,idir)&
        =reshape(source=recvbuffer_p(ibuf_recv_p:ibufnext-1),&
       shape=shape(pstateCo(igrid)%ws%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
          ixRmin3:ixRmax3,idir)))
     ibuf_recv_p=ibufnext

     call pole_copy_stg(pstateCo(igrid)%ws,ixRmin1,ixRmin2,ixRmin3,ixRmax1,&
        ixRmax2,ixRmax3,pole_buf_stg,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,&
        ixRmax3,ipole,i1,i2,i3,idir)

  end do
  

end if

end if !! ipole == 0

end subroutine bc_fill_p
end subroutine fill_gc_p
!===============================================================================
subroutine send_gc_p(igrid)
use mod_comm_gc
use mod_amrvacdef
! ibuf_send and isend_srl were already initialized in init_comm_gc_srl,
! here they are just updated 
integer :: igrid, i1,i2,i3, ineighbor, ipe_neighbor, my_neighbor_type
!-------------------------------------------------------------------------------
call set_tmpGlobals(igrid)

do i3=-1,1
do i2=-1,1
do i1=-1,1
   if (i1==0.and.i2==0.and.i3==0) cycle
   my_neighbor_type=neighbor_type(i1,i2,i3,igrid)
   if (my_neighbor_type==4) call bc_send_prolong
end do
end do
end do



contains
!=============================================================================
subroutine bc_send_prolong
integer :: ic1,ic2,ic3,inc1,inc2,inc3,n_i1,n_i2,n_i3,n_inc1,n_inc2,n_inc3,ii1,&
   ii2,ii3,idir,ibufaux,ibufnext,ipole
integer :: ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3
 ! Auxialiary array to avoid problems with the preprocessor...
 ! (It keeps splitting the line in the wrong place)
integer, dimension(1) :: sizes
!-----------------------------------------------------------------------------

ipole=neighbor_pole(i1,i2,i3,igrid)


do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
  inc3=2*i3+ic3
do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
  inc2=2*i2+ic2
do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
  inc1=2*i1+ic1

  ineighbor=neighbor_child(1,inc1,inc2,inc3,igrid)
  ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)

  if (ipe_neighbor.ne.mype) then
    select case (ipole)
    case(0) ! No pole
       n_i1=-i1;n_i2=-i2;n_i3=-i3;
       n_inc1=ic1+n_i1;n_inc2=ic2+n_i2;n_inc3=ic3+n_i3;
   case (1) ! Pole in this direction
       n_inc1=inc1;n_inc2=ic2-i2;n_inc3=ic3-i3;
!       print *, 'inc',inc1,inc2,inc3,'-->n_inc',n_inc1,n_inc2,n_inc3
       
   case (2) ! Pole in this direction
       n_inc1=ic1-i1;n_inc2=inc2;n_inc3=ic3-i3;
!       print *, 'inc',inc1,inc2,inc3,'-->n_inc',n_inc1,n_inc2,n_inc3
       
   case (3) ! Pole in this direction
       n_inc1=ic1-i1;n_inc2=ic2-i2;n_inc3=inc3;
!       print *, 'inc',inc1,inc2,inc3,'-->n_inc',n_inc1,n_inc2,n_inc3
       
    end select
  
   ! fill corresponding part of the send buffer...
  
    ixSmin1=ixS_p_min1(inc1);ixSmin2=ixS_p_min2(inc2)
    ixSmin3=ixS_p_min3(inc3);ixSmax1=ixS_p_max1(inc1)
    ixSmax2=ixS_p_max2(inc2);ixSmax3=ixS_p_max3(inc3);
  
    ibufaux=ibuf_send_p
    ibufnext=ibufaux+sizes_p_send(inc1,inc2,inc3)
  
    sizes=(/sizes_p_send(inc1,inc2,inc3)/)
  
    sendbuffer_p(ibufaux:ibufnext-1)=reshape(pstate(igrid)%w%w&
       (ixSmin1:ixSmax1,ixSmin2:ixSmax2,ixSmin3:ixSmax3,1:nwflux+nwaux),sizes)
  
    ibufaux=ibufnext
  
    do idir=1,3
      ixSmin1=ixS_p_stg_min1(idir,inc1);ixSmin2=ixS_p_stg_min2(idir,inc2)
      ixSmin3=ixS_p_stg_min3(idir,inc3);ixSmax1=ixS_p_stg_max1(idir,inc1)
      ixSmax2=ixS_p_stg_max2(idir,inc2);ixSmax3=ixS_p_stg_max3(idir,inc3);
      ibufnext=ibufaux+sizes_p_send_stg(idir,inc1,inc2,inc3)
      sizes=(/sizes_p_send_stg(idir,inc1,inc2,inc3)/)
      sendbuffer_p(ibufaux:ibufnext-1)=&
        reshape(pstate(igrid)%ws%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
           ixSmin3:ixSmax3,idir),sizes)   
      ibufaux=ibufnext
    end do
  

    ! ...and send
    itag=(3**3+4**3)*(ineighbor-1)+3**3+n_inc1*4**(1-1)+n_inc2*4&
       **(2-1)+n_inc3*4**(3-1)

    isend_p=isend_p+1
    call MPI_ISEND(sendbuffer_p(ibuf_send_p),sizes_p_send_total(inc1,inc2,&
       inc3), MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
       sendrequest_p(isend_p),ierrmpi)
  
    ibuf_send_p=ibufnext

  end if

end do
end do
end do

end subroutine bc_send_prolong
end subroutine send_gc_p
!===============================================================================
subroutine fill_boundary(igrid)
! Physical boundary conditions
use mod_comm_gc
use mod_amrvacdef
integer, intent(in) :: igrid
! ... local ...
integer :: idims,iside,i1,i2,i3,kmin1,kmin2,kmin3,kmax1,kmax2,kmax3,ixBmin1,&
   ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3
logical :: isphysbound
!-----------------------------------------------------------------------------

   call set_tmpGlobals(igrid)

   do idims=1,ndim
      ! to avoid using as yet unknown corner info in more than 1D, we
      ! fill only interior mesh ranges of the ghost cell ranges at first,
      ! and progressively enlarge the ranges to include corners later
      kmin1=0; kmax1=0;
      
      
       kmin2=merge(1, 0, idims .lt. 2 .and. neighbor_type(0,-1,0,igrid)==1)
       kmax2=merge(1, 0, idims .lt. 2 .and. neighbor_type(0, 1,0,igrid)==1)
       kmin3=merge(1, 0, idims .lt. 3 .and. neighbor_type(0,0,-1,igrid)==1)
       kmax3=merge(1, 0, idims .lt. 3 .and. neighbor_type(0,0, 1,igrid)==1)
      ixBmin1=ixGlo1+kmin1*dixB;ixBmin2=ixGlo2+kmin2*dixB
      ixBmin3=ixGlo3+kmin3*dixB;
      ixBmax1=ixGhi1-kmax1*dixB;ixBmax2=ixGhi2-kmax2*dixB
      ixBmax3=ixGhi3-kmax3*dixB;
      
      do iside=1,2
         i1=kr(1,idims)*(2*iside-3);i2=kr(2,idims)*(2*iside-3)
         i3=kr(3,idims)*(2*iside-3);
         if (neighbor_type(i1,i2,i3,igrid)/=1) cycle
         
         call bc_phys(iside,idims,time,pstate(igrid),ixBmin1,ixBmin2,ixBmin3,&
            ixBmax1,ixBmax2,ixBmax3)
         
      end do
         
   end do

end subroutine fill_boundary
!=============================================================================
subroutine physbound(i1,i2,i3,igrid,isphysbound)
use mod_forest
use mod_amrvacdef

integer, intent(in)  :: i1,i2,i3, igrid
logical, intent(out) :: isphysbound
type(tree_node_ptr)  :: tree
integer              :: level, ig1,ig2,ig3, ign1,ign2,ign3
!-----------------------------------------------------------------------------
isphysbound = .false.

tree%node => igrid_to_node(igrid,mype)%node
level = tree%node%level
ig1 = tree%node%ig1; ig2 = tree%node%ig2; ig3 = tree%node%ig3;

ign1 = ig1 + i1; ign2 = ig2 + i2; ign3 = ig3 + i3;
if (ign1 .gt. ng1(level) .or. ign1 .lt. 1.or.ign2 .gt. ng2(level) &
   .or. ign2 .lt. 1.or.ign3 .gt. ng3(level) .or. ign3 .lt. 1) isphysbound &
   = .true.

end subroutine physbound
!=============================================================================
subroutine pole_copy(pwrecv,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,&
   pwsend,ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3,ipole,i1,i2,i3)
use mod_amrvacdef

integer, intent(in) :: ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,&
    ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3, ipole, i1,i2,i3
type(walloc) :: pwrecv, pwsend

integer :: iw, iB, iside
!-----------------------------------------------------------------------------
select case (ipole)
case (1)
   iside=int((i1+3)/2)
   iB=2*(1-1)+iside
   do iw=1,nwflux+nwaux
      select case (typeB(iw,iB))
      case ("symm","polefix")
         pwrecv%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,iw) &
            = pwsend%w(ixSmax1:ixSmin1:-1,ixSmin2:ixSmax2,ixSmin3:ixSmax3,iw)
      case ("asymm")
         pwrecv%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,iw) &
            =-pwsend%w(ixSmax1:ixSmin1:-1,ixSmin2:ixSmax2,ixSmin3:ixSmax3,iw)
      case default
         call mpistop("Boundary condition at pole should be symm,&
            asymm or polefix")
      end select
   end do 
case (2)
   iside=int((i2+3)/2)
   iB=2*(2-1)+iside
   do iw=1,nwflux+nwaux
      select case (typeB(iw,iB))
      case ("symm","polefix")
         pwrecv%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,iw) &
            = pwsend%w(ixSmin1:ixSmax1,ixSmax2:ixSmin2:-1,ixSmin3:ixSmax3,iw)
      case ("asymm")
         pwrecv%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,iw) &
            =-pwsend%w(ixSmin1:ixSmax1,ixSmax2:ixSmin2:-1,ixSmin3:ixSmax3,iw)
      case default
         call mpistop("Boundary condition at pole should be symm,&
            asymm or polefix")
      end select
   end do 
case (3)
   iside=int((i3+3)/2)
   iB=2*(3-1)+iside
   do iw=1,nwflux+nwaux
      select case (typeB(iw,iB))
      case ("symm","polefix")
         pwrecv%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,iw) &
            = pwsend%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,ixSmax3:ixSmin3:-1,iw)
      case ("asymm")
         pwrecv%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,iw) &
            =-pwsend%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,ixSmax3:ixSmin3:-1,iw)
      case default
         call mpistop("Boundary condition at pole should be symm,&
            asymm or polefix")
      end select
   end do 
end select

end subroutine pole_copy
!=============================================================================
subroutine pole_copy_stg(pwrecv,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,&
   ixRmax3,pwsend,ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3,ipole,i1,i2,&
   i3,idir)
use mod_amrvacdef

integer, intent(in) :: ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,&
    ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3, ipole, i1,i2,i3, idir
type(walloc) :: pwrecv, pwsend

integer :: iw, iB, iside
!-----------------------------------------------------------------------------

iw=idir+b0_

select case (ipole)
case (1)
   iside=int((i1+3)/2)
   iB=2*(1-1)+iside
   select case (typeB(iw,iB))
      case ("symm","polefix")
         pwrecv%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,idir) &
            = pwsend%w(ixSmax1:ixSmin1:-1,ixSmin2:ixSmax2,ixSmin3:ixSmax3,&
            idir)
      case ("asymm")
         pwrecv%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,idir) &
            =-pwsend%w(ixSmax1:ixSmin1:-1,ixSmin2:ixSmax2,ixSmin3:ixSmax3,&
            idir)
      case default
         call mpistop("Boundary condition at pole should be symm,&
            asymm or polefix")
   end select

case (2)
   iside=int((i2+3)/2)
   iB=2*(2-1)+iside
   select case (typeB(iw,iB))
      case ("symm","polefix")
         pwrecv%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,idir) &
            = pwsend%w(ixSmin1:ixSmax1,ixSmax2:ixSmin2:-1,ixSmin3:ixSmax3,&
            idir)
      case ("asymm")
         pwrecv%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,idir) &
            =-pwsend%w(ixSmin1:ixSmax1,ixSmax2:ixSmin2:-1,ixSmin3:ixSmax3,&
            idir)
      case default
         call mpistop("Boundary condition at pole should be symm,&
            asymm or polefix")
   end select

case (3)
   iside=int((i3+3)/2)
   iB=2*(3-1)+iside
   select case (typeB(iw,iB))
      case ("symm","polefix")
         pwrecv%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,idir) &
            = pwsend%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,ixSmax3:ixSmin3:-1,&
            idir)
      case ("asymm")
         pwrecv%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,idir) &
            =-pwsend%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,ixSmax3:ixSmin3:-1,&
            idir)
      case default
         call mpistop("Boundary condition at pole should be symm,&
            asymm or polefix")
   end select

end select
end subroutine pole_copy_stg
!=============================================================================
subroutine fix_auxiliary(igrid)
use mod_comm_gc
use mod_amrvacdef

integer, intent(in) :: igrid
! ... local ...
! Counters and auxiliaries
integer :: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,i1,i2,i3
!-----------------------------------------------------------------------------

associate(x=>pstate(igrid)%x%x,pgeoFi=>pstate(igrid)%geo)

call set_tmpGlobals(igrid)

saveigrid=igrid
      
do i3=-1,1
do i2=-1,1
do i1=-1,1
! same-level auxiliaries have been sent and are not altered:
   if ((i1==0.and.i2==0.and.i3==0).or.neighbor_type(i1,i2,i3,&
      igrid).eq.3) cycle

   ixmin1=ixR_srl_min1(i1);ixmin2=ixR_srl_min2(i2);ixmin3=ixR_srl_min3(i3)
   ixmax1=ixR_srl_max1(i1);ixmax2=ixR_srl_max2(i2);ixmax3=ixR_srl_max3(i3);
   if(.not.slab)mygeo=>pgeoFi
   if(covariant)myM => mygeo%m
   call getaux(.true.,pstate(igrid)%w%w,x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
      ixGhi3,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,"bc")
end do
end do
end do

end associate

end subroutine fix_auxiliary
!=============================================================================


