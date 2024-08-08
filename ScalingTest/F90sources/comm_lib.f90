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
subroutine comm_start

use mod_amrvacdef
!-----------------------------------------------------------------------------
call MPI_INIT(ierrmpi)
call MPI_COMM_RANK(MPI_COMM_WORLD,mype,ierrmpi)
call MPI_COMM_SIZE(MPI_COMM_WORLD,npe,ierrmpi)

icomm=MPI_COMM_WORLD

end subroutine comm_start
!=============================================================================
subroutine comm_finalize

use mod_amrvacdef
!-----------------------------------------------------------------------------
call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
call MPI_FINALIZE(ierrmpi)

end subroutine comm_finalize
!=============================================================================
subroutine init_comm_types

use mod_amrvacdef

integer, dimension(ndim+1) :: sizes, subsizes, start
!integer :: i^D, ic^D, nx^D, nxCo^D, size_double
integer :: i1,i2,i3, ic1,ic2,ic3, nx1,nx2,nx3, nxCo1,nxCo2,nxCo3
integer(kind=MPI_ADDRESS_KIND):: size_double, lb
 integer :: idir
!-----------------------------------------------------------------------------
nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;nx3=ixMhi3-ixMlo3+1;
nxCo1=nx1/2;nxCo2=nx2/2;nxCo3=nx3/2;

sizes(1)=ixGhi1;sizes(2)=ixGhi2;sizes(3)=ixGhi3;
sizes(ndim+1)=nw
subsizes(1)=nx1;subsizes(2)=nx2;subsizes(3)=nx3;
subsizes(ndim+1)=nwflux
start(1)=ixMlo1-1;start(2)=ixMlo2-1;start(3)=ixMlo3-1;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
   MPI_DOUBLE_PRECISION, type_block,ierrmpi)
call MPI_TYPE_COMMIT(type_block,ierrmpi)

sizes(1)=ixGhi1/2+dixB;sizes(2)=ixGhi2/2+dixB;sizes(3)=ixGhi3/2+dixB;
sizes(ndim+1)=nw
subsizes(1)=nxCo1;subsizes(2)=nxCo2;subsizes(3)=nxCo3;
subsizes(ndim+1)=nw
start(1)=ixMlo1-1;start(2)=ixMlo2-1;start(3)=ixMlo3-1;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
   MPI_DOUBLE_PRECISION, type_coarse_block,ierrmpi)
call MPI_TYPE_COMMIT(type_coarse_block,ierrmpi)

sizes(1)=ixGhi1;sizes(2)=ixGhi2;sizes(3)=ixGhi3;
sizes(ndim+1)=nw
do ic3=1,2
do ic2=1,2
do ic1=1,2
   subsizes(1)=nxCo1;subsizes(2)=nxCo2;subsizes(3)=nxCo3;
   subsizes(ndim+1)=nw
   start(1)=ixMlo1-1+(ic1-1)*nxCo1;start(2)=ixMlo2-1+(ic2-1)*nxCo2
   start(3)=ixMlo3-1+(ic3-1)*nxCo3;
   start(ndim+1)=0
   call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start,&
       MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, type_sub_block(ic1,ic2,ic3),&
      ierrmpi)
   call MPI_TYPE_COMMIT(type_sub_block(ic1,ic2,ic3),ierrmpi)
end do
end do
end do

sizes(1)=ixGhi1;sizes(2)=ixGhi2;sizes(3)=ixGhi3;
sizes(ndim+1)=nw
subsizes(1)=nx1;subsizes(2)=nx2;subsizes(3)=nx3;
subsizes(ndim+1)=nw
start(1)=ixMlo1-1;start(2)=ixMlo2-1;start(3)=ixMlo3-1;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
   MPI_DOUBLE_PRECISION, type_block_io,ierrmpi)
call MPI_TYPE_COMMIT(type_block_io,ierrmpi)


!----------------------------------------
! Types for staggered grid:
!----------------------------------------

sizes(1)=ixGhi1+1;sizes(2)=ixGhi2+1;sizes(3)=ixGhi3+1;
sizes(ndim+1)=nws
subsizes(1)=nx1+1;subsizes(2)=nx2+1;subsizes(3)=nx3+1;
subsizes(ndim+1)=nws
start(1)=ixMlo1-1;start(2)=ixMlo2-1;start(3)=ixMlo3-1;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, &
                              MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                              type_block_staggered_io,ierrmpi)
call MPI_TYPE_COMMIT(type_block_staggered_io,ierrmpi)

!--------------------------------------------------------
! Types for sending coarse blocks for coarsening
sizes(1)=ixGhi1/2+dixB+1;sizes(2)=ixGhi2/2+dixB+1;sizes(3)=ixGhi3/2+dixB+1;
sizes(ndim+1)=nws

do ic3=1,2
do ic2=1,2
do ic1=1,2
do idir=1,3
  subsizes(1)=nxCo1+kr(ic1,1)*kr(idir,1)
  subsizes(2)=nxCo2+kr(ic2,1)*kr(idir,2)
  subsizes(3)=nxCo3+kr(ic3,1)*kr(idir,3);
  subsizes(ndim+1)=1
  start(1)=ixMlo1-kr(ic1,1)*kr(idir,1);start(2)=ixMlo2-kr(ic2,1)*kr(idir,2)
  start(3)=ixMlo3-kr(ic3,1)*kr(idir,3);
  start(ndim+1)=idir-1

  call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, &
                                MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                                type_coarse_block_stg(idir,ic1,ic2,ic3),&
                                   ierrmpi)
  call MPI_TYPE_COMMIT(type_coarse_block_stg(idir,ic1,ic2,ic3),ierrmpi)

end do
end do
end do
end do

! Types for receiving for coarsening
sizes(1)=ixGhi1+1;sizes(2)=ixGhi2+1;sizes(3)=ixGhi3+1;
sizes(ndim+1)=nws
do ic3=1,2
do ic2=1,2
do ic1=1,2
do idir=1,3
   subsizes(1)=nxCo1+kr(ic1,1)*kr(idir,1)
   subsizes(2)=nxCo2+kr(ic2,1)*kr(idir,2)
   subsizes(3)=nxCo3+kr(ic3,1)*kr(idir,3);
   subsizes(ndim+1)=1
   start(1)=ixMlo1-kr(ic1,1)*kr(idir,1)+(ic1-1)*nxCo1
   start(2)=ixMlo2-kr(ic2,1)*kr(idir,2)+(ic2-1)*nxCo2
   start(3)=ixMlo3-kr(ic3,1)*kr(idir,3)+(ic3-1)*nxCo3;
   start(ndim+1)=idir-1
   call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, &
                                 MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                                 type_sub_block_stg(idir,ic1,ic2,ic3),ierrmpi)
   call MPI_TYPE_COMMIT(type_sub_block_stg(idir,ic1,ic2,ic3),ierrmpi)
end do
end do
end do
end do

call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,size_double,ierrmpi)
size_block_stg_io=(nx1+1)*(nx2+1)*(nx3+1)*nws*size_double

!----------------------------------------


!call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,size_double,ierrmpi)
call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,size_double,ierrmpi)
size_block_io=nx1*nx2*nx3*nw*size_double



sizes(1)=ixMhi1-ixMlo1+1;sizes(2)=ixMhi2-ixMlo2+1;sizes(3)=ixMhi3-ixMlo3+1;
sizes(ndim+1)=3
subsizes(1)=sizes(1);subsizes(2)=sizes(2);subsizes(3)=sizes(3);
subsizes(ndim+1)=3
start(1)=0;start(2)=0;start(3)=0;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
   MPI_DOUBLE_PRECISION, type_block_xcc_io,ierrmpi)
call MPI_TYPE_COMMIT(type_block_xcc_io,ierrmpi)

sizes(1)=ixMhi1-ixMlo1+2;sizes(2)=ixMhi2-ixMlo2+2;sizes(3)=ixMhi3-ixMlo3+2;
sizes(ndim+1)=3
subsizes(1)=sizes(1);subsizes(2)=sizes(2);subsizes(3)=sizes(3);
subsizes(ndim+1)=3
start(1)=0;start(2)=0;start(3)=0;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
   MPI_DOUBLE_PRECISION, type_block_xc_io,ierrmpi)
call MPI_TYPE_COMMIT(type_block_xc_io,ierrmpi)

sizes(1)=ixMhi1-ixMlo1+1;sizes(2)=ixMhi2-ixMlo2+1;sizes(3)=ixMhi3-ixMlo3+1;
sizes(ndim+1)=nw+nwauxio
subsizes(1)=sizes(1);subsizes(2)=sizes(2);subsizes(3)=sizes(3);
subsizes(ndim+1)=nw+nwauxio
start(1)=0;start(2)=0;start(3)=0;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
   MPI_DOUBLE_PRECISION, type_block_wcc_io,ierrmpi)
call MPI_TYPE_COMMIT(type_block_wcc_io,ierrmpi)

sizes(1)=ixMhi1-ixMlo1+2;sizes(2)=ixMhi2-ixMlo2+2;sizes(3)=ixMhi3-ixMlo3+2;
sizes(ndim+1)=nw+nwauxio
subsizes(1)=sizes(1);subsizes(2)=sizes(2);subsizes(3)=sizes(3);
subsizes(ndim+1)=nw+nwauxio
start(1)=0;start(2)=0;start(3)=0;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
   MPI_DOUBLE_PRECISION, type_block_wc_io,ierrmpi)
call MPI_TYPE_COMMIT(type_block_wc_io,ierrmpi)

end subroutine init_comm_types
!=============================================================================
subroutine mpistop(message)

use mod_amrvacdef

character(len=*), intent(in) :: message

integer :: ierrcode
!------------------------------------------------------------------------------
write(*,*) "ERROR for processor",mype,":"
write(*,*) message
call MPI_ABORT(icomm,ierrcode,ierrmpi)

end subroutine mpistop
!==============================================================================
