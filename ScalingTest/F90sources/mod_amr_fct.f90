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

module mod_amr_fct
  use mod_indices, only: ngridshi
  implicit none
  
  type facealloc
    double precision, dimension(:,:,:), pointer :: face
  end type facealloc

  type fake_neighbors
    integer :: igrid
    integer :: ipe
  end type fake_neighbors

  type(facealloc), dimension(2,3,ngridshi), save :: pface

  type(fake_neighbors), dimension(2,2,2,3,ngridshi), save :: fine_neighbors

  integer, dimension(2,-1:11,-1:12,-1:13,ngridshi), save :: old_neighbor

  integer, save        :: ibuf_recv, ibuf_send, ibuf_send_next

contains


! -------------------------- subroutines ------------------------
!=============================================================================

subroutine prolong_2nd_stg(sCo,sFi,ixComin1in,ixComin2in,ixComin3in,&
   ixComax1in,ixComax2in,ixComax3in,ixFimin1in,ixFimin2in,ixFimin3in,&
   ixFimax1in,ixFimax2in,ixFimax3in,dxCo1,dxCo2,dxCo3,xComin1,xComin2,xComin3,&
   dxFi1,dxFi2,dxFi3,xFimin1,xFimin2,xFimin3,ghost,fine_min1in,fine_min2in,&
   fine_min3in,fine_max1in,fine_max2in,fine_max3in)
! This subroutine performs a 2nd order prolongation for a staggered field F,
! preserving the divergence of the coarse cell.
! This is useful for preserving DivF=0.
  ! If DivF=f(x), a different algorithm must be used.
use mod_amrvacdef

logical, intent(in)          :: ghost
integer, intent(in)          :: ixComin1in,ixComin2in,ixComin3in,ixComax1in,&
   ixComax2in,ixComax3in, ixFimin1in,ixFimin2in,ixFimin3in,ixFimax1in,&
   ixFimax2in,ixFimax3in
double precision, intent(in) :: dxCo1,dxCo2,dxCo3, xComin1,xComin2,xComin3,&
    dxFi1,dxFi2,dxFi3, xFimin1,xFimin2,xFimin3
type(state), intent(in)      :: sCo
type(state), intent(inout)   :: sFi

! Optional
logical, optional :: fine_min1in,fine_min2in,fine_min3in,fine_max1in,&
   fine_max2in,fine_max3in
logical           :: fine_min1,fine_min2,fine_min3,fine_max1,fine_max2,&
   fine_max3

! ... local ...

double precision :: eta1,eta2,eta3, invdxCo1,invdxCo2,invdxCo3
integer :: ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,ixComax3,ixFimin1,&
   ixFimin2,ixFimin3,ixFimax1,ixFimax2,ixFimax3
integer :: idim1,idim2,ix2,ix3,idim3,ixFismin1,ixFismin2,ixFismin3,ixFismax1,&
   ixFismax2,ixFismax3,ixGsmin1,ixGsmin2,ixGsmin3,ixGsmax1,ixGsmax2,ixGsmax3,&
   ixCosmin1,ixCosmin2,ixCosmin3,ixCosmax1,ixCosmax2,ixCosmax3,ixFisCmin1,&
   ixFisCmin2,ixFisCmin3,ixFisCmax1,ixFisCmax2,ixFisCmax3
integer :: ixCosVmin1(1:3),ixCosVmin2(1:3),ixCosVmin3(1:3),ixCosVmax1(1:3),&
   ixCosVmax2(1:3),ixCosVmax3(1:3),ixFisVmin1(1:3),ixFisVmin2(1:3),&
   ixFisVmin3(1:3),ixFisVmax1(1:3),ixFisVmax2(1:3),ixFisVmax3(1:3)
integer :: hxCosmin1,hxCosmin2,hxCosmin3,hxCosmax1,hxCosmax2,hxCosmax3,&
   jxCosmin1,jxCosmin2,jxCosmin3,jxCosmax1,jxCosmax2,jxCosmax3,ixCosEmin1,&
   ixCosEmin2,ixCosEmin3,ixCosEmax1,ixCosEmax2,ixCosEmax3,ixFisEmin1,&
   ixFisEmin2,ixFisEmin3,ixFisEmax1,ixFisEmax2,ixFisEmax3,hxFisCmin1,&
   hxFisCmin2,hxFisCmin3,hxFisCmax1,hxFisCmax2,hxFisCmax3,jxFisCmin1,&
   jxFisCmin2,jxFisCmin3,jxFisCmax1,jxFisCmax2,jxFisCmax3,ipxFisCmin1,&
   ipxFisCmin2,ipxFisCmin3,ipxFisCmax1,ipxFisCmax2,ipxFisCmax3,ixCosCmin1,&
   ixCosCmin2,ixCosCmin3,ixCosCmax1,ixCosCmax2,ixCosCmax3,imxFisCmin1,&
   imxFisCmin2,imxFisCmin3,imxFisCmax1,imxFisCmax2,imxFisCmax3,jpxFisCmin1,&
   jpxFisCmin2,jpxFisCmin3,jpxFisCmax1,jpxFisCmax2,jpxFisCmax3,jmxFisCmin1,&
   jmxFisCmin2,jmxFisCmin3,jmxFisCmax1,jmxFisCmax2,jmxFisCmax3,hpxFisCmin1,&
   hpxFisCmin2,hpxFisCmin3,hpxFisCmax1,hpxFisCmax2,hpxFisCmax3
integer :: hxFimin1,hxFimin2,hxFimin3,hxFimax1,hxFimax2,hxFimax3,jxFimin1,&
   jxFimin2,jxFimin3,jxFimax1,jxFimax2,jxFimax3,hijxFimin1,hijxFimin2,&
   hijxFimin3,hijxFimax1,hijxFimax2,hijxFimax3,hjixFimin1,hjixFimin2,&
   hjixFimin3,hjixFimax1,hjixFimax2,hjixFimax3,hjjxFimin1,hjjxFimin2,&
   hjjxFimin3,hjjxFimax1,hjjxFimax2,hjjxFimax3
integer :: iihxFimin1,iihxFimin2,iihxFimin3,iihxFimax1,iihxFimax2,iihxFimax3,&
   iijxFimin1,iijxFimin2,iijxFimin3,iijxFimax1,iijxFimax2,iijxFimax3,&
   ijhxFimin1,ijhxFimin2,ijhxFimin3,ijhxFimax1,ijhxFimax2,ijhxFimax3,&
   ijjxFimin1,ijjxFimin2,ijjxFimin3,ijjxFimax1,ijjxFimax2,ijjxFimax3,&
   ihixFimin1,ihixFimin2,ihixFimin3,ihixFimax1,ihixFimax2,ihixFimax3,&
   ijixFimin1,ijixFimin2,ijixFimin3,ijixFimax1,ijixFimax2,ijixFimax3,&
   ihjxFimin1,ihjxFimin2,ihjxFimin3,ihjxFimax1,ihjxFimax2,ihjxFimax3
integer :: jihxFimin1,jihxFimin2,jihxFimin3,jihxFimax1,jihxFimax2,jihxFimax3,&
   jijxFimin1,jijxFimin2,jijxFimin3,jijxFimax1,jijxFimax2,jijxFimax3,&
   jjhxFimin1,jjhxFimin2,jjhxFimin3,jjhxFimax1,jjhxFimax2,jjhxFimax3,&
   jjjxFimin1,jjjxFimin2,jjjxFimin3,jjjxFimax1,jjjxFimax2,jjjxFimax3,&
   jhixFimin1,jhixFimin2,jhixFimin3,jhixFimax1,jhixFimax2,jhixFimax3,&
   jjixFimin1,jjixFimin2,jjixFimin3,jjixFimax1,jjixFimax2,jjixFimax3,&
   jhjxFimin1,jhjxFimin2,jhjxFimin3,jhjxFimax1,jhjxFimax2,jhjxFimax3
double precision :: bfluxCo(sCo%ws%ixGmin1:sCo%ws%ixGmax1,&
   sCo%ws%ixGmin2:sCo%ws%ixGmax2,sCo%ws%ixGmin3:sCo%ws%ixGmax3,nws),&
   bfluxFi(sFi%ws%ixGmin1:sFi%ws%ixGmax1,sFi%ws%ixGmin2:sFi%ws%ixGmax2,&
   sFi%ws%ixGmin3:sFi%ws%ixGmax3,nws)
double precision :: slopes(sCo%ws%ixGmin1:sCo%ws%ixGmax1,&
   sCo%ws%ixGmin2:sCo%ws%ixGmax2,sCo%ws%ixGmin3:sCo%ws%ixGmax3,ndim)

! Directional bias, see pdf
 !These arrays should have ranges ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3
 !Since they are still not known, we give ranges ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,
! even though it is wasteful
double precision :: sigmau(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
   sigmad(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),sigma(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:ndim), alpha(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
   ixGlo3:ixGhi3,1:ndim)
! Auxiliary arrays for magnetic fluxes
double precision :: F1(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
   F2(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),F3(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2,ixGlo3:ixGhi3),F4(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)


!-----------------------------------------------------------------------------


! Note on the indices:
! ixCo  Cells where 
!       divergence-preserving prolongation will be applied.
! ixFi  Fine cells which correspond to that extent.
! ixCoE For 'expanded', faces that need to be used to calculate
!       the slopes for the prolongation
! ixCos
! ixFis
! ixCosV And
! ixFisV For 'variable', since the ranges depend on the direction of the faces
! ixCosC For 'copy',
 !ixFisC For 'copy', to fill the information from the fine grid needed in the internal faces prolongation step.
! At the end, ixFisC is used to copy the assign the magnetic
! fields components to the state structure.

fine_min1=.false.;fine_min2=.false.;fine_min3=.false.;fine_max1=.false.
fine_max2=.false.;fine_max3=.false.;

if(present(fine_min1in))fine_min1=fine_min1in
if(present(fine_min2in))fine_min2=fine_min2in
if(present(fine_min3in))fine_min3=fine_min3in;
if(present(fine_max1in))fine_max1=fine_max1in
if(present(fine_max2in))fine_max2=fine_max2in
if(present(fine_max3in))fine_max3=fine_max3in;

 !When filling ghost cells, ixFimin1,ixFimin2,ixFimin3,ixFimax1,ixFimax2,ixFimax3 are given.
 !When refining the block, ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,ixComax3 are given.
if (ghost) then
ixFimin1=ixFimin1in;ixFimin2=ixFimin2in;ixFimin3=ixFimin3in
ixFimax1=ixFimax1in;ixFimax2=ixFimax2in;ixFimax3=ixFimax3in;
ixComin1=int((ixFimin1+dixB+1)/2);ixComin2=int((ixFimin2+dixB+1)/2)
ixComin3=int((ixFimin3+dixB+1)/2);
ixComax1=int((ixFimax1+dixB+1)/2);ixComax2=int((ixFimax2+dixB+1)/2)
ixComax3=int((ixFimax3+dixB+1)/2);
else
ixComin1=ixComin1in;ixComin2=ixComin2in;ixComin3=ixComin3in
ixComax1=ixComax1in;ixComax2=ixComax2in;ixComax3=ixComax3in;
ixFimin1=ixMlo1;ixFimin2=ixMlo2;ixFimin3=ixMlo3;ixFimax1=ixMhi1
ixFimax2=ixMhi2;ixFimax3=ixMhi3;
end if

! Expanded range for staggered variables
ixCosmin1=ixComin1-1;ixCosmin2=ixComin2-1;ixCosmin3=ixComin3-1;
ixCosmax1=ixComax1;ixCosmax2=ixComax2;ixCosmax3=ixComax3;

ixFismin1=ixFimin1-1;ixFismin2=ixFimin2-1;ixFismin3=ixFimin3-1;
ixFismax1=ixFimax1;ixFismax2=ixFimax2;ixFismax3=ixFimax3;


associate(wCos=>sCo%ws%w, wFis=>sFi%ws%w,wCo=>sCo%w%w, wFi=>sFi%w%w)
! Assemble general indices
ixGsmin1=sFi%ws%ixGmin1;ixGsmin2=sFi%ws%ixGmin2;ixGsmin3=sFi%ws%ixGmin3;
ixGsmax1=sFi%ws%ixGmax1;ixGsmax2=sFi%ws%ixGmax2;ixGsmax3=sFi%ws%ixGmax3;
! -------------------------

do idim1=1,3
ixCosVmin1(idim1)=ixComin1-kr(1,idim1);ixCosVmin2(idim1)=ixComin2-kr(2,idim1)
ixCosVmin3(idim1)=ixComin3-kr(3,idim1);
ixCosVmax1(idim1)=ixComax1;ixCosVmax2(idim1)=ixComax2
ixCosVmax3(idim1)=ixComax3;
ixFisVmin1(idim1)=ixFimin1-kr(1,idim1);ixFisVmin2(idim1)=ixFimin2-kr(2,idim1)
ixFisVmin3(idim1)=ixFimin3-kr(3,idim1);
ixFisVmax1(idim1)=ixFimax1;ixFisVmax2(idim1)=ixFimax2
ixFisVmax3(idim1)=ixFimax3;
end do

! Initialize auxiliary arrays at zero
bfluxCo = zero
bfluxFi = zero 
slopes  = zero 


invdxCo1=1.d0/dxCo1;invdxCo2=1.d0/dxCo2;invdxCo3=1.d0/dxCo3;

! Fill coarse magnetic flux array
do idim1=1,ndim
  ! Fill information from parts already at the fine level
  ! First set up indices
  if (ghost) then
    ixFisEmin1=max(1-kr(1,idim1),ixFisVmin1(idim1)-2*(1-kr(1,idim1)))
    ixFisEmin2=max(1-kr(2,idim1),ixFisVmin2(idim1)-2*(1-kr(2,idim1)))
    ixFisEmin3=max(1-kr(3,idim1),ixFisVmin3(idim1)-2*(1-kr(3,idim1)));
    ixFisEmax1=min(ixGsmax1,ixFisVmax1(idim1)+2*(1-kr(1,idim1)))
    ixFisEmax2=min(ixGsmax2,ixFisVmax2(idim1)+2*(1-kr(2,idim1)))
    ixFisEmax3=min(ixGsmax3,ixFisVmax3(idim1)+2*(1-kr(3,idim1)));
    ixCosEmin1=int((ixFisEmin1+dixB+1)/2)
    ixCosEmin2=int((ixFisEmin2+dixB+1)/2)
    ixCosEmin3=int((ixFisEmin3+dixB+1)/2)
    ixCosEmax1=int((ixFisEmax1+dixB+1)/2)
    ixCosEmax2=int((ixFisEmax2+dixB+1)/2)
    ixCosEmax3=int((ixFisEmax3+dixB+1)/2);
  else
    ixCosEmin1=ixCosVmin1(idim1)-(1-kr(idim1,1))
    ixCosEmin2=ixCosVmin2(idim1)-(1-kr(idim1,2))
    ixCosEmin3=ixCosVmin3(idim1)-(1-kr(idim1,3))
    ixCosEmax1=ixCosVmax1(idim1)+(1-kr(idim1,1))
    ixCosEmax2=ixCosVmax2(idim1)+(1-kr(idim1,2))
    ixCosEmax3=ixCosVmax3(idim1)+(1-kr(idim1,3));
    ixFisEmin1=ixFisVmin1(idim1)-2*(1-kr(idim1,1))
    ixFisEmin2=ixFisVmin2(idim1)-2*(1-kr(idim1,2))
    ixFisEmin3=ixFisVmin3(idim1)-2*(1-kr(idim1,3))
    ixFisEmax1=ixFisVmax1(idim1)+2*(1-kr(idim1,1))
    ixFisEmax2=ixFisVmax2(idim1)+2*(1-kr(idim1,2))
    ixFisEmax3=ixFisVmax3(idim1)+2*(1-kr(idim1,3));
  end if

  ! Convert fine fields to fluxes
  select case(idim1)
 case(1)
  bfluxFi(ixFisEmin1:ixFisEmax1,ixFisEmin2:ixFisEmax2,ixFisEmin3:ixFisEmax3,&
     idim1)=wFis(ixFisEmin1:ixFisEmax1,ixFisEmin2:ixFisEmax2,&
     ixFisEmin3:ixFisEmax3,1)*sFi%geo%surfaceC1(ixFisEmin1:ixFisEmax1,&
     ixFisEmin2:ixFisEmax2,ixFisEmin3:ixFisEmax3)
 
case(2)
  bfluxFi(ixFisEmin1:ixFisEmax1,ixFisEmin2:ixFisEmax2,ixFisEmin3:ixFisEmax3,&
     idim1)=wFis(ixFisEmin1:ixFisEmax1,ixFisEmin2:ixFisEmax2,&
     ixFisEmin3:ixFisEmax3,2)*sFi%geo%surfaceC2(ixFisEmin1:ixFisEmax1,&
     ixFisEmin2:ixFisEmax2,ixFisEmin3:ixFisEmax3)
 
case(3)
  bfluxFi(ixFisEmin1:ixFisEmax1,ixFisEmin2:ixFisEmax2,ixFisEmin3:ixFisEmax3,&
     idim1)=wFis(ixFisEmin1:ixFisEmax1,ixFisEmin2:ixFisEmax2,&
     ixFisEmin3:ixFisEmax3,3)*sFi%geo%surfaceC3(ixFisEmin1:ixFisEmax1,&
     ixFisEmin2:ixFisEmax2,ixFisEmin3:ixFisEmax3)
 
  end select

  
  
  idim2=1+mod(idim1,3)
  idim3=1+mod(idim1+1,3)
 
  bfluxCo(ixCosEmin1:ixCosEmax1,ixCosEmin2:ixCosEmax2,ixCosEmin3:ixCosEmax3,&
     idim1) = zero
  ! Add fine fluxes sharing the same fine face
 do ix2=0,1
do ix3=0,1
  
  ixFisCmin1=ixFisEmin1+ix2*kr(idim2,1);ixFisCmin2=ixFisEmin2+ix2*kr(idim2,2)
  ixFisCmin3=ixFisEmin3+ix2*kr(idim2,3);ixFisCmax1=ixFisEmax1+ix2*kr(idim2,1)
  ixFisCmax2=ixFisEmax2+ix2*kr(idim2,2);ixFisCmax3=ixFisEmax3+ix2*kr(idim2,3);
  ixFisCmin1=ixFisCmin1+ix3*kr(idim3,1);ixFisCmin2=ixFisCmin2+ix3*kr(idim3,2)
  ixFisCmin3=ixFisCmin3+ix3*kr(idim3,3);ixFisCmax1=ixFisCmax1+ix3*kr(idim3,1)
  ixFisCmax2=ixFisCmax2+ix3*kr(idim3,2);ixFisCmax3=ixFisCmax3+ix3*kr(idim3,3);

  bfluxCo(ixCosEmin1:ixCosEmax1,ixCosEmin2:ixCosEmax2,ixCosEmin3:ixCosEmax3,&
     idim1)=&
  bfluxCo(ixCosEmin1:ixCosEmax1,ixCosEmin2:ixCosEmax2,ixCosEmin3:ixCosEmax3,&
     idim1)&
 +bfluxFi(ixFisCmin1:ixFisCmax1:2,ixFisCmin2:ixFisCmax2:2,&
    ixFisCmin3:ixFisCmax3:2,idim1)

 end do
end do

!print *, 'Before'
!print *, 'Fine, idim1: ',idim1
!call printarray(&
!0,!0,!0,12,12,12,&
!0,!0,!0,12,12,12,&
!bfluxFi(:,:,:,idim1))
!print *, 'Coarse, idim1: ',idim1
!call printarray(&
!0,!0,!0,8,8,8,&
!0,!0,!0,8,8,8,&
!bfluxCo(:,:,:,idim1))

end do

!---------------

! Omit indices for already refined edge, if any

if (fine_min1) then
ixCosVmin1(1)=ixCosVmin1(1)+1
ixFisVmin1(1)=ixFisVmin1(1)+2
end if
if (fine_max1) then
ixCosVmax1(1)=ixCosVmax1(1)-1
ixFisVmax1(1)=ixFisVmax1(1)-2
end if


if (fine_min2) then
ixCosVmin2(2)=ixCosVmin2(2)+1
ixFisVmin2(2)=ixFisVmin2(2)+2
end if
if (fine_max2) then
ixCosVmax2(2)=ixCosVmax2(2)-1
ixFisVmax2(2)=ixFisVmax2(2)-2
end if


if (fine_min3) then
ixCosVmin3(3)=ixCosVmin3(3)+1
ixFisVmin3(3)=ixFisVmin3(3)+2
end if
if (fine_max3) then
ixCosVmax3(3)=ixCosVmax3(3)-1
ixFisVmax3(3)=ixFisVmax3(3)-2
end if


!---------------

do idim1=1,ndim
  ixCosEmin1=ixCosVmin1(idim1);ixCosEmin2=ixCosVmin2(idim1)
  ixCosEmin3=ixCosVmin3(idim1);ixCosEmax1=ixCosVmax1(idim1)
  ixCosEmax2=ixCosVmax2(idim1);ixCosEmax3=ixCosVmax3(idim1);
  ! Omit part already refined
  
  if (1.ne.idim1) then
    if ((.not.fine_min1).or.(.not.ghost)) then
     ixCosEmin1=ixCosVmin1(idim1)-1
    end if
    if ((.not.fine_max1).or.(.not.ghost)) then
     ixCosEmax1=ixCosVmax1(idim1)+1
    end if
  end if
  

  if (2.ne.idim1) then
    if ((.not.fine_min2).or.(.not.ghost)) then
     ixCosEmin2=ixCosVmin2(idim1)-1
    end if
    if ((.not.fine_max2).or.(.not.ghost)) then
     ixCosEmax2=ixCosVmax2(idim1)+1
    end if
  end if
  

  if (3.ne.idim1) then
    if ((.not.fine_min3).or.(.not.ghost)) then
     ixCosEmin3=ixCosVmin3(idim1)-1
    end if
    if ((.not.fine_max3).or.(.not.ghost)) then
     ixCosEmax3=ixCosVmax3(idim1)+1
    end if
  end if
  
  ! Fill coarse flux array from coarse field
  select case(idim1)
 case(1)
      bfluxCo(ixCosEmin1:ixCosEmax1,ixCosEmin2:ixCosEmax2,&
         ixCosEmin3:ixCosEmax3,1)=wCos(ixCosEmin1:ixCosEmax1,&
         ixCosEmin2:ixCosEmax2,ixCosEmin3:ixCosEmax3,1)&
      *sCo%geo%surfaceC1(ixCosEmin1:ixCosEmax1,ixCosEmin2:ixCosEmax2,&
         ixCosEmin3:ixCosEmax3)

case(2)
      bfluxCo(ixCosEmin1:ixCosEmax1,ixCosEmin2:ixCosEmax2,&
         ixCosEmin3:ixCosEmax3,2)=wCos(ixCosEmin1:ixCosEmax1,&
         ixCosEmin2:ixCosEmax2,ixCosEmin3:ixCosEmax3,2)&
      *sCo%geo%surfaceC2(ixCosEmin1:ixCosEmax1,ixCosEmin2:ixCosEmax2,&
         ixCosEmin3:ixCosEmax3)

case(3)
      bfluxCo(ixCosEmin1:ixCosEmax1,ixCosEmin2:ixCosEmax2,&
         ixCosEmin3:ixCosEmax3,3)=wCos(ixCosEmin1:ixCosEmax1,&
         ixCosEmin2:ixCosEmax2,ixCosEmin3:ixCosEmax3,3)&
      *sCo%geo%surfaceC3(ixCosEmin1:ixCosEmax1,ixCosEmin2:ixCosEmax2,&
         ixCosEmin3:ixCosEmax3)

  end select

!----------------
!print *, 'After'
!print *, 'Fine, idim1: ',idim1
!call printarray(&
!0,!0,!0,12,12,12,&
!0,!0,!0,12,12,12,&
!bfluxFi(:,:,:,idim1))
!print *, 'Coarse, idim1: ',idim1
!call printarray(&
!0,!0,!0,8,8,8,&
!0,!0,!0,8,8,8,&
!bfluxCo(:,:,:,idim1))
!----------------

end do
! Finished filling coarse flux array


!---------------
! Distribute coarse fluxes among fine fluxes
! There are too many loops here, perhaps optimize later

do idim1=1,ndim
   ixCosmin1=ixCosVmin1(idim1);ixCosmin2=ixCosVmin2(idim1)
   ixCosmin3=ixCosVmin3(idim1);ixCosmax1=ixCosVmax1(idim1)
   ixCosmax2=ixCosVmax2(idim1);ixCosmax3=ixCosVmax3(idim1);
   do idim2=1,ndim
     if(idim1.eq.idim2) cycle
     ! Calculate slope in direction idim2
     ! Set up indices
     jxCosmin1=ixCosmin1+kr(idim2,1);jxCosmin2=ixCosmin2+kr(idim2,2)
     jxCosmin3=ixCosmin3+kr(idim2,3);jxCosmax1=ixCosmax1+kr(idim2,1)
     jxCosmax2=ixCosmax2+kr(idim2,2);jxCosmax3=ixCosmax3+kr(idim2,3);
     hxCosmin1=ixCosmin1-kr(idim2,1);hxCosmin2=ixCosmin2-kr(idim2,2)
     hxCosmin3=ixCosmin3-kr(idim2,3);hxCosmax1=ixCosmax1-kr(idim2,1)
     hxCosmax2=ixCosmax2-kr(idim2,2);hxCosmax3=ixCosmax3-kr(idim2,3);
     slopes(ixCosmin1:ixCosmax1,ixCosmin2:ixCosmax2,ixCosmin3:ixCosmax3,&
        idim2)=&
       0.125d0*(bfluxCo(jxCosmin1:jxCosmax1,jxCosmin2:jxCosmax2,&
          jxCosmin3:jxCosmax3,idim1)-bfluxCo(hxCosmin1:hxCosmax1,&
          hxCosmin2:hxCosmax2,hxCosmin3:hxCosmax3,idim1))


     do idim3=1,ndim
      ! Calculate slope in direction idim3
      ! Set up indices
      jxCosmin1=ixCosmin1+kr(idim3,1);jxCosmin2=ixCosmin2+kr(idim3,2)
      jxCosmin3=ixCosmin3+kr(idim3,3);jxCosmax1=ixCosmax1+kr(idim3,1)
      jxCosmax2=ixCosmax2+kr(idim3,2);jxCosmax3=ixCosmax3+kr(idim3,3);
      hxCosmin1=ixCosmin1-kr(idim3,1);hxCosmin2=ixCosmin2-kr(idim3,2)
      hxCosmin3=ixCosmin3-kr(idim3,3);hxCosmax1=ixCosmax1-kr(idim3,1)
      hxCosmax2=ixCosmax2-kr(idim3,2);hxCosmax3=ixCosmax3-kr(idim3,3);
      slopes(ixCosmin1:ixCosmax1,ixCosmin2:ixCosmax2,ixCosmin3:ixCosmax3,&
         idim3)=&
        0.125d0*(bfluxCo(jxCosmin1:jxCosmax1,jxCosmin2:jxCosmax2,&
           jxCosmin3:jxCosmax3,idim1)-bfluxCo(hxCosmin1:hxCosmax1,&
           hxCosmin2:hxCosmax2,hxCosmin3:hxCosmax3,idim1))
      if(lvc(idim1,idim2,idim3).lt.1) cycle
      do ix2=0,1
       do ix3=0,1
       ixFisCmin1=ixFisVmin1(idim1)&
        +ix2*kr(1,idim2)+ix3*kr(1,idim3);ixFisCmin2=ixFisVmin2(idim1)&
        +ix2*kr(2,idim2)+ix3*kr(2,idim3);ixFisCmin3=ixFisVmin3(idim1)&
        +ix2*kr(3,idim2)+ix3*kr(3,idim3);
       ixFisCmax1=ixFisVmax1(idim1)&
        +ix2*kr(1,idim2)+ix3*kr(1,idim3);ixFisCmax2=ixFisVmax2(idim1)&
        +ix2*kr(2,idim2)+ix3*kr(2,idim3);ixFisCmax3=ixFisVmax3(idim1)&
        +ix2*kr(3,idim2)+ix3*kr(3,idim3);
        bfluxFi(ixFisCmin1:ixFisCmax1:2,ixFisCmin2:ixFisCmax2:2,&
           ixFisCmin3:ixFisCmax3:2,idim1)=&
         quarter*bfluxCo(ixCosmin1:ixCosmax1,ixCosmin2:ixCosmax2,&
            ixCosmin3:ixCosmax3,idim1)
        bfluxFi(ixFisCmin1:ixFisCmax1:2,ixFisCmin2:ixFisCmax2:2,&
           ixFisCmin3:ixFisCmax3:2,idim1)=&
         bfluxFi(ixFisCmin1:ixFisCmax1:2,ixFisCmin2:ixFisCmax2:2,&
            ixFisCmin3:ixFisCmax3:2,idim1)&
        +quarter*(2*ix2-1)*slopes(ixCosmin1:ixCosmax1,ixCosmin2:ixCosmax2,&
           ixCosmin3:ixCosmax3,idim2)
        bfluxFi(ixFisCmin1:ixFisCmax1:2,ixFisCmin2:ixFisCmax2:2,&
           ixFisCmin3:ixFisCmax3:2,idim1)=&
         bfluxFi(ixFisCmin1:ixFisCmax1:2,ixFisCmin2:ixFisCmax2:2,&
            ixFisCmin3:ixFisCmax3:2,idim1)&
        +quarter*(2*ix3-1)*slopes(ixCosmin1:ixCosmax1,ixCosmin2:ixCosmax2,&
           ixCosmin3:ixCosmax3,idim3)

       end do
      end do
     end do

   end do
end do

! Calculate interior fine fluxes

! Directional bias for nonlinear prolongation
do idim1=1,ndim
  do idim2=1,ndim
    do idim3=1,ndim
      if (lvc(idim1,idim2,idim3).lt.1) cycle
      ! Set up indices
      hxFimin1=ixFimin1-kr(idim1,1);hxFimin2=ixFimin2-kr(idim1,2)
      hxFimin3=ixFimin3-kr(idim1,3);hxFimax1=ixFimax1-kr(idim1,1)
      hxFimax2=ixFimax2-kr(idim1,2);hxFimax3=ixFimax3-kr(idim1,3);
      jxFimin1=ixFimin1+kr(idim1,1);jxFimin2=ixFimin2+kr(idim1,2)
      jxFimin3=ixFimin3+kr(idim1,3);jxFimax1=ixFimax1+kr(idim1,1)
      jxFimax2=ixFimax2+kr(idim1,2);jxFimax3=ixFimax3+kr(idim1,3);
      ! ---
      hijxFimin1=hxFimin1+kr(idim3,1);hijxFimin2=hxFimin2+kr(idim3,2)
      hijxFimin3=hxFimin3+kr(idim3,3);hijxFimax1=hxFimax1+kr(idim3,1)
      hijxFimax2=hxFimax2+kr(idim3,2);hijxFimax3=hxFimax3+kr(idim3,3);
      hjixFimin1=hxFimin1+kr(idim2,1);hjixFimin2=hxFimin2+kr(idim2,2)
      hjixFimin3=hxFimin3+kr(idim2,3);hjixFimax1=hxFimax1+kr(idim2,1)
      hjixFimax2=hxFimax2+kr(idim2,2);hjixFimax3=hxFimax3+kr(idim2,3);
      hjjxFimin1=hijxFimin1+kr(idim2,1);hjjxFimin2=hijxFimin2+kr(idim2,2)
      hjjxFimin3=hijxFimin3+kr(idim2,3);hjjxFimax1=hijxFimax1+kr(idim2,1)
      hjjxFimax2=hijxFimax2+kr(idim2,2);hjjxFimax3=hijxFimax3+kr(idim2,3);
      ! ---
      iihxFimin1=ixFimin1-kr(idim3,1);iihxFimin2=ixFimin2-kr(idim3,2)
      iihxFimin3=ixFimin3-kr(idim3,3);iihxFimax1=ixFimax1-kr(idim3,1)
      iihxFimax2=ixFimax2-kr(idim3,2);iihxFimax3=ixFimax3-kr(idim3,3);
      iijxFimin1=ixFimin1+kr(idim3,1);iijxFimin2=ixFimin2+kr(idim3,2)
      iijxFimin3=ixFimin3+kr(idim3,3);iijxFimax1=ixFimax1+kr(idim3,1)
      iijxFimax2=ixFimax2+kr(idim3,2);iijxFimax3=ixFimax3+kr(idim3,3);
      ihixFimin1=ixFimin1-kr(idim2,1);ihixFimin2=ixFimin2-kr(idim2,2)
      ihixFimin3=ixFimin3-kr(idim2,3);ihixFimax1=ixFimax1-kr(idim2,1)
      ihixFimax2=ixFimax2-kr(idim2,2);ihixFimax3=ixFimax3-kr(idim2,3);
      ihjxFimin1=ihixFimin1+kr(idim3,1);ihjxFimin2=ihixFimin2+kr(idim3,2)
      ihjxFimin3=ihixFimin3+kr(idim3,3);ihjxFimax1=ihixFimax1+kr(idim3,1)
      ihjxFimax2=ihixFimax2+kr(idim3,2);ihjxFimax3=ihixFimax3+kr(idim3,3);
      ijixFimin1=ixFimin1+kr(idim2,1);ijixFimin2=ixFimin2+kr(idim2,2)
      ijixFimin3=ixFimin3+kr(idim2,3);ijixFimax1=ixFimax1+kr(idim2,1)
      ijixFimax2=ixFimax2+kr(idim2,2);ijixFimax3=ixFimax3+kr(idim2,3);
      ijhxFimin1=ijixFimin1-kr(idim3,1);ijhxFimin2=ijixFimin2-kr(idim3,2)
      ijhxFimin3=ijixFimin3-kr(idim3,3);ijhxFimax1=ijixFimax1-kr(idim3,1)
      ijhxFimax2=ijixFimax2-kr(idim3,2);ijhxFimax3=ijixFimax3-kr(idim3,3);
      ijjxFimin1=ijixFimin1+kr(idim3,1);ijjxFimin2=ijixFimin2+kr(idim3,2)
      ijjxFimin3=ijixFimin3+kr(idim3,3);ijjxFimax1=ijixFimax1+kr(idim3,1)
      ijjxFimax2=ijixFimax2+kr(idim3,2);ijjxFimax3=ijixFimax3+kr(idim3,3);
      ! ---
      jihxFimin1=jxFimin1-kr(idim3,1);jihxFimin2=jxFimin2-kr(idim3,2)
      jihxFimin3=jxFimin3-kr(idim3,3);jihxFimax1=jxFimax1-kr(idim3,1)
      jihxFimax2=jxFimax2-kr(idim3,2);jihxFimax3=jxFimax3-kr(idim3,3);
      jijxFimin1=jxFimin1+kr(idim3,1);jijxFimin2=jxFimin2+kr(idim3,2)
      jijxFimin3=jxFimin3+kr(idim3,3);jijxFimax1=jxFimax1+kr(idim3,1)
      jijxFimax2=jxFimax2+kr(idim3,2);jijxFimax3=jxFimax3+kr(idim3,3);
      jhixFimin1=jxFimin1-kr(idim2,1);jhixFimin2=jxFimin2-kr(idim2,2)
      jhixFimin3=jxFimin3-kr(idim2,3);jhixFimax1=jxFimax1-kr(idim2,1)
      jhixFimax2=jxFimax2-kr(idim2,2);jhixFimax3=jxFimax3-kr(idim2,3);
      jhjxFimin1=jhixFimin1+kr(idim3,1);jhjxFimin2=jhixFimin2+kr(idim3,2)
      jhjxFimin3=jhixFimin3+kr(idim3,3);jhjxFimax1=jhixFimax1+kr(idim3,1)
      jhjxFimax2=jhixFimax2+kr(idim3,2);jhjxFimax3=jhixFimax3+kr(idim3,3);
      jjixFimin1=jxFimin1+kr(idim2,1);jjixFimin2=jxFimin2+kr(idim2,2)
      jjixFimin3=jxFimin3+kr(idim2,3);jjixFimax1=jxFimax1+kr(idim2,1)
      jjixFimax2=jxFimax2+kr(idim2,2);jjixFimax3=jxFimax3+kr(idim2,3);
      jjhxFimin1=jjixFimin1-kr(idim3,1);jjhxFimin2=jjixFimin2-kr(idim3,2)
      jjhxFimin3=jjixFimin3-kr(idim3,3);jjhxFimax1=jjixFimax1-kr(idim3,1)
      jjhxFimax2=jjixFimax2-kr(idim3,2);jjhxFimax3=jjixFimax3-kr(idim3,3);
      jjjxFimin1=jjixFimin1+kr(idim3,1);jjjxFimin2=jjixFimin2+kr(idim3,2)
      jjjxFimin3=jjixFimin3+kr(idim3,3);jjjxFimax1=jjixFimax1+kr(idim3,1)
      jjjxFimax2=jjixFimax2+kr(idim3,2);jjjxFimax3=jjixFimax3+kr(idim3,3);
      ! ---
      sigmau(ixComin1:ixComax1,ixComin2:ixComax2,ixComin3:ixComax3)=&
        abs(bfluxFi(jhixFimin1:jhixFimax1:2,jhixFimin2:jhixFimax2:2,&
           jhixFimin3:jhixFimax3:2,idim2))&
        + abs(bfluxFi(jhjxFimin1:jhjxFimax1:2,jhjxFimin2:jhjxFimax2:2,&
           jhjxFimin3:jhjxFimax3:2,idim2))&
        + abs(bfluxFi(jjixFimin1:jjixFimax1:2,jjixFimin2:jjixFimax2:2,&
           jjixFimin3:jjixFimax3:2,idim2))&
        + abs(bfluxFi(jjjxFimin1:jjjxFimax1:2,jjjxFimin2:jjjxFimax2:2,&
           jjjxFimin3:jjjxFimax3:2,idim2))&
        + abs(bfluxFi(jihxFimin1:jihxFimax1:2,jihxFimin2:jihxFimax2:2,&
           jihxFimin3:jihxFimax3:2,idim3))&
        + abs(bfluxFi(jijxFimin1:jijxFimax1:2,jijxFimin2:jijxFimax2:2,&
           jijxFimin3:jijxFimax3:2,idim3))&
        + abs(bfluxFi(jjhxFimin1:jjhxFimax1:2,jjhxFimin2:jjhxFimax2:2,&
           jjhxFimin3:jjhxFimax3:2,idim3))&
        + abs(bfluxFi(jjjxFimin1:jjjxFimax1:2,jjjxFimin2:jjjxFimax2:2,&
           jjjxFimin3:jjjxFimax3:2,idim3))

      sigmad(ixComin1:ixComax1,ixComin2:ixComax2,ixComin3:ixComax3)=&
        abs(bfluxFi(ihixFimin1:ihixFimax1:2,ihixFimin2:ihixFimax2:2,&
           ihixFimin3:ihixFimax3:2,idim2))&
        + abs(bfluxFi(ihjxFimin1:ihjxFimax1:2,ihjxFimin2:ihjxFimax2:2,&
           ihjxFimin3:ihjxFimax3:2,idim2))&
        + abs(bfluxFi(ijixFimin1:ijixFimax1:2,ijixFimin2:ijixFimax2:2,&
           ijixFimin3:ijixFimax3:2,idim2))&
        + abs(bfluxFi(ijjxFimin1:ijjxFimax1:2,ijjxFimin2:ijjxFimax2:2,&
           ijjxFimin3:ijjxFimax3:2,idim2))&
        + abs(bfluxFi(iihxFimin1:iihxFimax1:2,iihxFimin2:iihxFimax2:2,&
           iihxFimin3:iihxFimax3:2,idim3))&
        + abs(bfluxFi(iijxFimin1:iijxFimax1:2,iijxFimin2:iijxFimax2:2,&
           iijxFimin3:iijxFimax3:2,idim3))&
        + abs(bfluxFi(ijhxFimin1:ijhxFimax1:2,ijhxFimin2:ijhxFimax2:2,&
           ijhxFimin3:ijhxFimax3:2,idim3))&
        + abs(bfluxFi(ijjxFimin1:ijjxFimax1:2,ijjxFimin2:ijjxFimax2:2,&
           ijjxFimin3:ijjxFimax3:2,idim3))

      sigma(ixComin1:ixComax1,ixComin2:ixComax2,ixComin3:ixComax3,idim1)&
         =sigmau(ixComin1:ixComax1,ixComin2:ixComax2,ixComin3:ixComax3)&
         +sigmad(ixComin1:ixComax1,ixComin2:ixComax2,ixComin3:ixComax3)
      where(sigma(ixComin1:ixComax1,ixComin2:ixComax2,ixComin3:ixComax3,&
         idim1).ne.zero)
        sigma(ixComin1:ixComax1,ixComin2:ixComax2,ixComin3:ixComax3,idim1)&
           =abs(sigmau(ixComin1:ixComax1,ixComin2:ixComax2,&
           ixComin3:ixComax3)-sigmad(ixComin1:ixComax1,ixComin2:ixComax2,&
           ixComin3:ixComax3))/sigma(ixComin1:ixComax1,ixComin2:ixComax2,&
           ixComin3:ixComax3,idim1)
      elsewhere
        sigma(ixComin1:ixComax1,ixComin2:ixComax2,ixComin3:ixComax3,idim1)&
           =zero
      end where

    end do
  end do
end do

! Directional bias for Toth-Roe prolongation

do idim1=1, ndim
  do idim2=1, ndim
    do idim3=1, ndim
      if (lvc(idim1,idim2,idim3).lt.1) cycle
!       Nonlinear
 !alpha(ixComin1:ixComax1,ixComin2:ixComax2,ixComin3:ixComax3,idim1)=sigma(ixComin1:ixComax1,ixComin2:ixComax2,ixComin3:ixComax3,idim2)-sigma(ixComin1:ixComax1,ixComin2:ixComax2,ixComin3:ixComax3,idim3)
!       Homogeneous
 !alpha(ixComin1:ixComax1,ixComin2:ixComax2,ixComin3:ixComax3,idim1)=0.d0
!       Toth-Roe
        alpha(ixComin1:ixComax1,ixComin2:ixComax2,ixComin3:ixComax3,idim1)&
           =(sCo%geo%dx(ixComin1:ixComax1,ixComin2:ixComax2,ixComin3:ixComax3,&
           idim2)-sCo%geo%dx(ixComin1:ixComax1,ixComin2:ixComax2,&
           ixComin3:ixComax3,idim3))/(sCo%geo%dx(ixComin1:ixComax1,&
           ixComin2:ixComax2,ixComin3:ixComax3,idim2)+sCo%geo%dx&
           (ixComin1:ixComax1,ixComin2:ixComax2,ixComin3:ixComax3,idim3))
    end do
  end do
end do

 !call printarray(ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,ixComax3,ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,ixComax3,alpha(ixComin1:ixComax1,ixComin2:ixComax2,ixComin3:ixComax3,1))
 !call printarray(ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,ixComax3,ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,ixComax3,alpha(ixComin1:ixComax1,ixComin2:ixComax2,ixComin3:ixComax3,2))
 !call printarray(ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,ixComax3,ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,ixComax3,alpha(ixComin1:ixComax1,ixComin2:ixComax2,ixComin3:ixComax3,3))


do idim1=1,ndim
  do idim2=1,ndim


    do idim3=1,ndim
      if (lvc(idim1,idim2,idim3).lt.1) cycle
      ! Set up indices
      hxFimin1=ixFimin1-kr(idim1,1);hxFimin2=ixFimin2-kr(idim1,2)
      hxFimin3=ixFimin3-kr(idim1,3);hxFimax1=ixFimax1-kr(idim1,1)
      hxFimax2=ixFimax2-kr(idim1,2);hxFimax3=ixFimax3-kr(idim1,3);
      jxFimin1=ixFimin1+kr(idim1,1);jxFimin2=ixFimin2+kr(idim1,2)
      jxFimin3=ixFimin3+kr(idim1,3);jxFimax1=ixFimax1+kr(idim1,1)
      jxFimax2=ixFimax2+kr(idim1,2);jxFimax3=ixFimax3+kr(idim1,3);
      ! ---
      hijxFimin1=hxFimin1+kr(idim3,1);hijxFimin2=hxFimin2+kr(idim3,2)
      hijxFimin3=hxFimin3+kr(idim3,3);hijxFimax1=hxFimax1+kr(idim3,1)
      hijxFimax2=hxFimax2+kr(idim3,2);hijxFimax3=hxFimax3+kr(idim3,3);
      hjixFimin1=hxFimin1+kr(idim2,1);hjixFimin2=hxFimin2+kr(idim2,2)
      hjixFimin3=hxFimin3+kr(idim2,3);hjixFimax1=hxFimax1+kr(idim2,1)
      hjixFimax2=hxFimax2+kr(idim2,2);hjixFimax3=hxFimax3+kr(idim2,3);
      hjjxFimin1=hijxFimin1+kr(idim2,1);hjjxFimin2=hijxFimin2+kr(idim2,2)
      hjjxFimin3=hijxFimin3+kr(idim2,3);hjjxFimax1=hijxFimax1+kr(idim2,1)
      hjjxFimax2=hijxFimax2+kr(idim2,2);hjjxFimax3=hijxFimax3+kr(idim2,3);
      ! ---
      iihxFimin1=ixFimin1-kr(idim3,1);iihxFimin2=ixFimin2-kr(idim3,2)
      iihxFimin3=ixFimin3-kr(idim3,3);iihxFimax1=ixFimax1-kr(idim3,1)
      iihxFimax2=ixFimax2-kr(idim3,2);iihxFimax3=ixFimax3-kr(idim3,3);
      iijxFimin1=ixFimin1+kr(idim3,1);iijxFimin2=ixFimin2+kr(idim3,2)
      iijxFimin3=ixFimin3+kr(idim3,3);iijxFimax1=ixFimax1+kr(idim3,1)
      iijxFimax2=ixFimax2+kr(idim3,2);iijxFimax3=ixFimax3+kr(idim3,3);
      ihixFimin1=ixFimin1-kr(idim2,1);ihixFimin2=ixFimin2-kr(idim2,2)
      ihixFimin3=ixFimin3-kr(idim2,3);ihixFimax1=ixFimax1-kr(idim2,1)
      ihixFimax2=ixFimax2-kr(idim2,2);ihixFimax3=ixFimax3-kr(idim2,3);
     ihjxFimin1=ihixFimin1+kr(idim3,1);ihjxFimin2=ihixFimin2+kr(idim3,2)
     ihjxFimin3=ihixFimin3+kr(idim3,3);ihjxFimax1=ihixFimax1+kr(idim3,1)
     ihjxFimax2=ihixFimax2+kr(idim3,2);ihjxFimax3=ihixFimax3+kr(idim3,3);
      ijixFimin1=ixFimin1+kr(idim2,1);ijixFimin2=ixFimin2+kr(idim2,2)
      ijixFimin3=ixFimin3+kr(idim2,3);ijixFimax1=ixFimax1+kr(idim2,1)
      ijixFimax2=ixFimax2+kr(idim2,2);ijixFimax3=ixFimax3+kr(idim2,3);
      ijhxFimin1=ijixFimin1-kr(idim3,1);ijhxFimin2=ijixFimin2-kr(idim3,2)
      ijhxFimin3=ijixFimin3-kr(idim3,3);ijhxFimax1=ijixFimax1-kr(idim3,1)
      ijhxFimax2=ijixFimax2-kr(idim3,2);ijhxFimax3=ijixFimax3-kr(idim3,3);
      ijjxFimin1=ijixFimin1+kr(idim3,1);ijjxFimin2=ijixFimin2+kr(idim3,2)
      ijjxFimin3=ijixFimin3+kr(idim3,3);ijjxFimax1=ijixFimax1+kr(idim3,1)
      ijjxFimax2=ijixFimax2+kr(idim3,2);ijjxFimax3=ijixFimax3+kr(idim3,3);
      ! ---
      jihxFimin1=jxFimin1-kr(idim3,1);jihxFimin2=jxFimin2-kr(idim3,2)
      jihxFimin3=jxFimin3-kr(idim3,3);jihxFimax1=jxFimax1-kr(idim3,1)
      jihxFimax2=jxFimax2-kr(idim3,2);jihxFimax3=jxFimax3-kr(idim3,3);
      jijxFimin1=jxFimin1+kr(idim3,1);jijxFimin2=jxFimin2+kr(idim3,2)
      jijxFimin3=jxFimin3+kr(idim3,3);jijxFimax1=jxFimax1+kr(idim3,1)
      jijxFimax2=jxFimax2+kr(idim3,2);jijxFimax3=jxFimax3+kr(idim3,3);
      jhixFimin1=jxFimin1-kr(idim2,1);jhixFimin2=jxFimin2-kr(idim2,2)
      jhixFimin3=jxFimin3-kr(idim2,3);jhixFimax1=jxFimax1-kr(idim2,1)
      jhixFimax2=jxFimax2-kr(idim2,2);jhixFimax3=jxFimax3-kr(idim2,3);
      jhjxFimin1=jhixFimin1+kr(idim3,1);jhjxFimin2=jhixFimin2+kr(idim3,2)
      jhjxFimin3=jhixFimin3+kr(idim3,3);jhjxFimax1=jhixFimax1+kr(idim3,1)
      jhjxFimax2=jhixFimax2+kr(idim3,2);jhjxFimax3=jhixFimax3+kr(idim3,3);
      jjixFimin1=jxFimin1+kr(idim2,1);jjixFimin2=jxFimin2+kr(idim2,2)
      jjixFimin3=jxFimin3+kr(idim2,3);jjixFimax1=jxFimax1+kr(idim2,1)
      jjixFimax2=jxFimax2+kr(idim2,2);jjixFimax3=jxFimax3+kr(idim2,3);
      jjhxFimin1=jjixFimin1-kr(idim3,1);jjhxFimin2=jjixFimin2-kr(idim3,2)
      jjhxFimin3=jjixFimin3-kr(idim3,3);jjhxFimax1=jjixFimax1-kr(idim3,1)
      jjhxFimax2=jjixFimax2-kr(idim3,2);jjhxFimax3=jjixFimax3-kr(idim3,3);
      jjjxFimin1=jjixFimin1+kr(idim3,1);jjjxFimin2=jjixFimin2+kr(idim3,2)
      jjjxFimin3=jjixFimin3+kr(idim3,3);jjjxFimax1=jjixFimax1+kr(idim3,1)
      jjjxFimax2=jjixFimax2+kr(idim3,2);jjjxFimax3=jjixFimax3+kr(idim3,3);
      ! ---
      ! Prolongation formulas
      F1(ixComin1:ixComax1,ixComin2:ixComax2,ixComin3:ixComax3)&
         =bfluxFi(ihixFimin1:ihixFimax1:2,ihixFimin2:ihixFimax2:2,&
         ihixFimin3:ihixFimax3:2,idim2)&
                -bfluxFi(jhixFimin1:jhixFimax1:2,jhixFimin2:jhixFimax2:2,&
                   jhixFimin3:jhixFimax3:2,idim2)&
                -bfluxFi(ijixFimin1:ijixFimax1:2,ijixFimin2:ijixFimax2:2,&
                   ijixFimin3:ijixFimax3:2,idim2)&
                +bfluxFi(jjixFimin1:jjixFimax1:2,jjixFimin2:jjixFimax2:2,&
                   jjixFimin3:jjixFimax3:2,idim2)

      F2(ixComin1:ixComax1,ixComin2:ixComax2,ixComin3:ixComax3)&
         =bfluxFi(ihjxFimin1:ihjxFimax1:2,ihjxFimin2:ihjxFimax2:2,&
         ihjxFimin3:ihjxFimax3:2,idim2)&
                -bfluxFi(jhjxFimin1:jhjxFimax1:2,jhjxFimin2:jhjxFimax2:2,&
                   jhjxFimin3:jhjxFimax3:2,idim2)&
                -bfluxFi(ijjxFimin1:ijjxFimax1:2,ijjxFimin2:ijjxFimax2:2,&
                   ijjxFimin3:ijjxFimax3:2,idim2)&
                +bfluxFi(jjjxFimin1:jjjxFimax1:2,jjjxFimin2:jjjxFimax2:2,&
                   jjjxFimin3:jjjxFimax3:2,idim2)

      F3(ixComin1:ixComax1,ixComin2:ixComax2,ixComin3:ixComax3)&
         =bfluxFi(iihxFimin1:iihxFimax1:2,iihxFimin2:iihxFimax2:2,&
         iihxFimin3:iihxFimax3:2,idim3)&
                -bfluxFi(jihxFimin1:jihxFimax1:2,jihxFimin2:jihxFimax2:2,&
                   jihxFimin3:jihxFimax3:2,idim3)&
                -bfluxFi(iijxFimin1:iijxFimax1:2,iijxFimin2:iijxFimax2:2,&
                   iijxFimin3:iijxFimax3:2,idim3)&
                +bfluxFi(jijxFimin1:jijxFimax1:2,jijxFimin2:jijxFimax2:2,&
                   jijxFimin3:jijxFimax3:2,idim3)

      F4(ixComin1:ixComax1,ixComin2:ixComax2,ixComin3:ixComax3)&
         =bfluxFi(ijhxFimin1:ijhxFimax1:2,ijhxFimin2:ijhxFimax2:2,&
         ijhxFimin3:ijhxFimax3:2,idim3)&
                -bfluxFi(jjhxFimin1:jjhxFimax1:2,jjhxFimin2:jjhxFimax2:2,&
                   jjhxFimin3:jjhxFimax3:2,idim3)&
                -bfluxFi(ijjxFimin1:ijjxFimax1:2,ijjxFimin2:ijjxFimax2:2,&
                   ijjxFimin3:ijjxFimax3:2,idim3)&
                +bfluxFi(jjjxFimin1:jjjxFimax1:2,jjjxFimin2:jjjxFimax2:2,&
                   jjjxFimin3:jjjxFimax3:2,idim3)

      bfluxFi(ixFimin1:ixFimax1:2,ixFimin2:ixFimax2:2,ixFimin3:ixFimax3:2,&
         idim1)=&
         half*(bfluxFi(hxFimin1:ixFimax1:2,hxFimin2:ixFimax2:2,&
            hxFimin3:ixFimax3:2,idim1)+bfluxFi(jxFimin1:ixFimax1:2,&
            jxFimin2:ixFimax2:2,jxFimin3:ixFimax3:2,idim1))&
         +6.25d-2*((3.d0+alpha(ixComin1:ixComax1,ixComin2:ixComax2,&
            ixComin3:ixComax3,idim2))*F1(ixComin1:ixComax1,ixComin2:ixComax2,&
            ixComin3:ixComax3)&
                  +(1.d0-alpha(ixComin1:ixComax1,ixComin2:ixComax2,&
                     ixComin3:ixComax3,idim2))*F2(ixComin1:ixComax1,&
                     ixComin2:ixComax2,ixComin3:ixComax3)&
                  +(3.d0-alpha(ixComin1:ixComax1,ixComin2:ixComax2,&
                     ixComin3:ixComax3,idim3))*F3(ixComin1:ixComax1,&
                     ixComin2:ixComax2,ixComin3:ixComax3)&
                  +(1.d0+alpha(ixComin1:ixComax1,ixComin2:ixComax2,&
                     ixComin3:ixComax3,idim3))*F4(ixComin1:ixComax1,&
                     ixComin2:ixComax2,ixComin3:ixComax3))

      bfluxFi(ijixFimin1:ijixFimax1:2,ijixFimin2:ijixFimax2:2,&
         ijixFimin3:ijixFimax3:2,idim1)=&
         half*(bfluxFi(hjixFimin1:hjixFimax1:2,hjixFimin2:hjixFimax2:2,&
            hjixFimin3:hjixFimax3:2,idim1)+bfluxFi(jjixFimin1:jjixFimax1:2,&
            jjixFimin2:jjixFimax2:2,jjixFimin3:jjixFimax3:2,idim1))&
         +6.25d-2*((3.d0+alpha(ixComin1:ixComax1,ixComin2:ixComax2,&
            ixComin3:ixComax3,idim2))*F1(ixComin1:ixComax1,ixComin2:ixComax2,&
            ixComin3:ixComax3)&
                  +(1.d0-alpha(ixComin1:ixComax1,ixComin2:ixComax2,&
                     ixComin3:ixComax3,idim2))*F2(ixComin1:ixComax1,&
                     ixComin2:ixComax2,ixComin3:ixComax3)&
                  +(1.d0+alpha(ixComin1:ixComax1,ixComin2:ixComax2,&
                     ixComin3:ixComax3,idim3))*F3(ixComin1:ixComax1,&
                     ixComin2:ixComax2,ixComin3:ixComax3)&
                  +(3.d0-alpha(ixComin1:ixComax1,ixComin2:ixComax2,&
                     ixComin3:ixComax3,idim3))*F4(ixComin1:ixComax1,&
                     ixComin2:ixComax2,ixComin3:ixComax3))

      bfluxFi(iijxFimin1:iijxFimax1:2,iijxFimin2:iijxFimax2:2,&
         iijxFimin3:iijxFimax3:2,idim1)=&
         half*(bfluxFi(hijxFimin1:hijxFimax1:2,hijxFimin2:hijxFimax2:2,&
            hijxFimin3:hijxFimax3:2,idim1)+bfluxFi(jijxFimin1:jijxFimax1:2,&
            jijxFimin2:jijxFimax2:2,jijxFimin3:jijxFimax3:2,idim1))&
         +6.25d-2*((1.d0-alpha(ixComin1:ixComax1,ixComin2:ixComax2,&
            ixComin3:ixComax3,idim2))*F1(ixComin1:ixComax1,ixComin2:ixComax2,&
            ixComin3:ixComax3)&
                  +(3.d0+alpha(ixComin1:ixComax1,ixComin2:ixComax2,&
                     ixComin3:ixComax3,idim2))*F2(ixComin1:ixComax1,&
                     ixComin2:ixComax2,ixComin3:ixComax3)&
                  +(3.d0-alpha(ixComin1:ixComax1,ixComin2:ixComax2,&
                     ixComin3:ixComax3,idim3))*F3(ixComin1:ixComax1,&
                     ixComin2:ixComax2,ixComin3:ixComax3)&
                  +(1.d0+alpha(ixComin1:ixComax1,ixComin2:ixComax2,&
                     ixComin3:ixComax3,idim3))*F4(ixComin1:ixComax1,&
                     ixComin2:ixComax2,ixComin3:ixComax3))

      bfluxFi(ijjxFimin1:ijjxFimax1:2,ijjxFimin2:ijjxFimax2:2,&
         ijjxFimin3:ijjxFimax3:2,idim1)=&
         half*(bfluxFi(hjjxFimin1:hjjxFimax1:2,hjjxFimin2:hjjxFimax2:2,&
            hjjxFimin3:hjjxFimax3:2,idim1)+bfluxFi(jjjxFimin1:jjjxFimax1:2,&
            jjjxFimin2:jjjxFimax2:2,jjjxFimin3:jjjxFimax3:2,idim1))&
         +6.25d-2*((1.d0-alpha(ixComin1:ixComax1,ixComin2:ixComax2,&
            ixComin3:ixComax3,idim2))*F1(ixComin1:ixComax1,ixComin2:ixComax2,&
            ixComin3:ixComax3)&
                  +(3.d0+alpha(ixComin1:ixComax1,ixComin2:ixComax2,&
                     ixComin3:ixComax3,idim2))*F2(ixComin1:ixComax1,&
                     ixComin2:ixComax2,ixComin3:ixComax3)&
                  +(1.d0+alpha(ixComin1:ixComax1,ixComin2:ixComax2,&
                     ixComin3:ixComax3,idim3))*F3(ixComin1:ixComax1,&
                     ixComin2:ixComax2,ixComin3:ixComax3)&
                  +(3.d0-alpha(ixComin1:ixComax1,ixComin2:ixComax2,&
                     ixComin3:ixComax3,idim3))*F4(ixComin1:ixComax1,&
                     ixComin2:ixComax2,ixComin3:ixComax3))
    end do

  end do
end do



! Go back to magnetic fields

do idim1=1,ndim
   ixFisCmax1=ixFimax1;ixFisCmax2=ixFimax2;ixFisCmax3=ixFimax3;
   ixFisCmin1=ixFimin1-kr(1,idim1);ixFisCmin2=ixFimin2-kr(2,idim1)
   ixFisCmin3=ixFimin3-kr(3,idim1);
   select case(idim1)
   case(1)
    where(sFi%geo%surfaceC1(ixFisCmin1:ixFisCmax1,ixFisCmin2:ixFisCmax2,&
       ixFisCmin3:ixFisCmax3).ne.zero)
      wFis(ixFisCmin1:ixFisCmax1,ixFisCmin2:ixFisCmax2,ixFisCmin3:ixFisCmax3,&
         1)=bfluxFi(ixFisCmin1:ixFisCmax1,ixFisCmin2:ixFisCmax2,&
         ixFisCmin3:ixFisCmax3,1)/sFi%geo%surfaceC1(ixFisCmin1:ixFisCmax1,&
         ixFisCmin2:ixFisCmax2,ixFisCmin3:ixFisCmax3)
    elsewhere
      wFis(ixFisCmin1:ixFisCmax1,ixFisCmin2:ixFisCmax2,ixFisCmin3:ixFisCmax3,&
         1)=zero
    end where
   case(2)
    where(sFi%geo%surfaceC2(ixFisCmin1:ixFisCmax1,ixFisCmin2:ixFisCmax2,&
       ixFisCmin3:ixFisCmax3).ne.zero)
      wFis(ixFisCmin1:ixFisCmax1,ixFisCmin2:ixFisCmax2,ixFisCmin3:ixFisCmax3,&
         2)=bfluxFi(ixFisCmin1:ixFisCmax1,ixFisCmin2:ixFisCmax2,&
         ixFisCmin3:ixFisCmax3,2)/sFi%geo%surfaceC2(ixFisCmin1:ixFisCmax1,&
         ixFisCmin2:ixFisCmax2,ixFisCmin3:ixFisCmax3)
    elsewhere
      wFis(ixFisCmin1:ixFisCmax1,ixFisCmin2:ixFisCmax2,ixFisCmin3:ixFisCmax3,&
         2)=zero
    end where
   case(3)
    where(sFi%geo%surfaceC3(ixFisCmin1:ixFisCmax1,ixFisCmin2:ixFisCmax2,&
       ixFisCmin3:ixFisCmax3).ne.zero)
      wFis(ixFisCmin1:ixFisCmax1,ixFisCmin2:ixFisCmax2,ixFisCmin3:ixFisCmax3,&
         3)=bfluxFi(ixFisCmin1:ixFisCmax1,ixFisCmin2:ixFisCmax2,&
         ixFisCmin3:ixFisCmax3,3)/sFi%geo%surfaceC3(ixFisCmin1:ixFisCmax1,&
         ixFisCmin2:ixFisCmax2,ixFisCmin3:ixFisCmax3)
    elsewhere
      wFis(ixFisCmin1:ixFisCmax1,ixFisCmin2:ixFisCmax2,ixFisCmin3:ixFisCmax3,&
         3)=zero
    end where
  
   end select
end do

call faces2centers(ixFimin1,ixFimin2,ixFimin3,ixFimax1,ixFimax2,ixFimax3,sFi)

end associate


end subroutine prolong_2nd_stg
!==================================================================
subroutine store_faces

use mod_forest, only: refine 

use mod_amrvacdef

! To achive consistency and thus conservation of divergence,
! when refining a block we take into account the faces of the
! already fine neighbours, if any.
! This routine stores them.

integer :: igrid, iigrid, idims, iside, ineighbor, ipe_neighbor
integer :: nx1,nx2,nx3, i1,i2,i3, ic1,ic2,ic3, inc1,inc2,inc3
!------------------------------------------------------------------
if (npe.gt.1) then
  nsend_fc=0
  nrecv_fc=0
end if

! Size of the block face
nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;nx3=ixMhi3-ixMlo3+1;

do iigrid=1,igridstail; igrid=igrids(iigrid);
  ! Check whether it is necessary to store any block face, i.e.
  ! if any coarser neighbour is going to be refined 
 do iside=1,2
    i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3);i3=kr(3,1)*(2*iside-3);

    if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle
    if (neighbor_type(i1,i2,i3,igrid).eq.2) then
      ineighbor   =neighbor(1,i1,i2,i3,igrid)
      ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
      if (refine(ineighbor,ipe_neighbor)) then
 !print *, 'My neighbor will be refined! igrig, ipe',ineighbor,ipe_neighbor
        !! Allocate
        allocate(pface(iside,1,igrid)%face(1,1:nx2,1:nx3))

        !! Store the faces
        if (iside.eq.1) then !! left side
        pface(iside,1,igrid)%face(1,1:nx2,1:nx3)=&
          ps(igrid)%ws%w(ixMlo1-1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,1)
        else !! right side
        pface(iside,1,igrid)%face(1,1:nx2,1:nx3)=&
          ps(igrid)%ws%w(ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,1)
        end if
       ! print *, 'allocated iside, idims, igrid',iside,1,igrid

        if (ipe_neighbor.ne.mype) nsend_fc(1)=nsend_fc(1)+1

      end if
    end if
  end do
do iside=1,2
    i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3);i3=kr(3,2)*(2*iside-3);

    if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle
    if (neighbor_type(i1,i2,i3,igrid).eq.2) then
      ineighbor   =neighbor(1,i1,i2,i3,igrid)
      ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
      if (refine(ineighbor,ipe_neighbor)) then
 !print *, 'My neighbor will be refined! igrig, ipe',ineighbor,ipe_neighbor
        !! Allocate
        allocate(pface(iside,2,igrid)%face(1:nx1,1,1:nx3))

        !! Store the faces
        if (iside.eq.1) then !! left side
        pface(iside,2,igrid)%face(1:nx1,1,1:nx3)=&
          ps(igrid)%ws%w(ixMlo1:ixMhi1,ixMlo2-1,ixMlo3:ixMhi3,2)
        else !! right side
        pface(iside,2,igrid)%face(1:nx1,1,1:nx3)=&
          ps(igrid)%ws%w(ixMlo1:ixMhi1,ixMhi2,ixMlo3:ixMhi3,2)
        end if
       ! print *, 'allocated iside, idims, igrid',iside,2,igrid

        if (ipe_neighbor.ne.mype) nsend_fc(2)=nsend_fc(2)+1

      end if
    end if
  end do
do iside=1,2
    i1=kr(1,3)*(2*iside-3);i2=kr(2,3)*(2*iside-3);i3=kr(3,3)*(2*iside-3);

    if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle
    if (neighbor_type(i1,i2,i3,igrid).eq.2) then
      ineighbor   =neighbor(1,i1,i2,i3,igrid)
      ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
      if (refine(ineighbor,ipe_neighbor)) then
 !print *, 'My neighbor will be refined! igrig, ipe',ineighbor,ipe_neighbor
        !! Allocate
        allocate(pface(iside,3,igrid)%face(1:nx1,1:nx2,1))

        !! Store the faces
        if (iside.eq.1) then !! left side
        pface(iside,3,igrid)%face(1:nx1,1:nx2,1)=&
          ps(igrid)%ws%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3-1,3)
        else !! right side
        pface(iside,3,igrid)%face(1:nx1,1:nx2,1)=&
          ps(igrid)%ws%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMhi3,3)
        end if
       ! print *, 'allocated iside, idims, igrid',iside,3,igrid

        if (ipe_neighbor.ne.mype) nsend_fc(3)=nsend_fc(3)+1

      end if
    end if
  end do

  ! If a grid is going to be refined,
  ! remember what are its neighbours.
  if (refine(igrid,mype)) then
    ! Initialize neighbour array
    fine_neighbors(:,:,:,:,igrid)%igrid=-1
    fine_neighbors(:,:,:,:,igrid)%ipe=-1
    do idims=1,3
      do iside=1,2
        i1=kr(1,idims)*(2*iside-3);i2=kr(2,idims)*(2*iside-3)
        i3=kr(3,idims)*(2*iside-3);
 if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle
        if (neighbor_type(i1,i2,i3,igrid).eq.4) then
         do ic3=1+int((1+i3)/2),2-int((1-i3)/2)
          inc3=ic3+i3
do ic2=1+int((1+i2)/2),2-int((1-i2)/2)
          inc2=ic2+i2
do ic1=1+int((1+i1)/2),2-int((1-i1)/2)
          inc1=ic1+i1
            ineighbor=neighbor_child(1,inc1,inc2,inc3,igrid)
            ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)

            fine_neighbors(ic1,ic2,ic3,idims,igrid)%igrid= ineighbor
            fine_neighbors(ic1,ic2,ic3,idims,igrid)%ipe=ipe_neighbor

            if (ipe_neighbor.ne.mype) nrecv_fc(idims)=nrecv_fc(idims)+1

           ! print *, 'I will be refined, igrid',igrid
           ! print *, 'icD,ineighbor',ic1,ic2,ic3
           ! print *, 'incD,ineighbor',inc1,inc2,inc3

 !print *, 'ineighbor',fine_neighbors(ic1,ic2,ic3,idims,igrid)%igrid
!            print *, 'ipe_neigh',fine_neighbors(ic1,ic2,ic3,idims,igrid)%ipe
         end do
end do
end do
        end if
      end do
    end do
  end if

end do

end subroutine store_faces
!==================================================================
subroutine comm_faces
!! When refining a coarse block with fine neighbours, it is necessary
!! prolong consistently with the already fine faces.
!! This routine takes care of the communication of such faces.
!! THIS SHOULD BE UPGRAFED TO 'MAKE_TASK' KINS OF ROUTINES,
!! IN THE SAME WAY AS THE FLUX FIXING
use mod_forest, only: refine
use mod_fix_conserve, only: nrecv, nsend, isize,&
                            recvbuffer, sendbuffer
use mod_amrvacdef

integer                   :: iigrid,igrid,ineighbor,ipe_neighbor
integer                   :: idims,iside,i1,i2,i3,ic1,ic2,ic3,inc1,inc2,inc3,&
   nx1,nx2,nx3
integer                   :: recvsize, sendsize

! Communicate the block faces to achieve consistency when refining
!-----------------------------------------------------------------
! Initialize communication structures
!print *,'mype, nrecv_fc, nsend_fc 1',mype,nrecv_fc(1),nsend_fc(1)
!print *,'mype, nrecv_fc, nsend_fc 2',mype,nrecv_fc(2),nsend_fc(2)

nrecv=0
nsend=0
recvsize=0
sendsize=0

do idims=1,3
   select case (idims)
   case (1)
      nrecv=nrecv+nrecv_fc(1)
      nsend=nsend+nsend_fc(1)
      nx1=1;nx2=ixMhi2-ixMlo2+1;nx3=ixMhi3-ixMlo3+1;
      isize(1)=nx1*nx2*nx3
      recvsize=recvsize+nrecv_fc(1)*isize(1)
      sendsize=sendsize+nsend_fc(1)*isize(1)
   
case (2)
      nrecv=nrecv+nrecv_fc(2)
      nsend=nsend+nsend_fc(2)
      nx1=ixMhi1-ixMlo1+1;nx2=1;nx3=ixMhi3-ixMlo3+1;
      isize(2)=nx1*nx2*nx3
      recvsize=recvsize+nrecv_fc(2)*isize(2)
      sendsize=sendsize+nsend_fc(2)*isize(2)
   
case (3)
      nrecv=nrecv+nrecv_fc(3)
      nsend=nsend+nsend_fc(3)
      nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;nx3=1;
      isize(3)=nx1*nx2*nx3
      recvsize=recvsize+nrecv_fc(3)*isize(3)
      sendsize=sendsize+nsend_fc(3)*isize(3)
   
   end select
end do

!print *, 'mype, recvsize, sendsize',mype, recvsize, sendsize

if (nrecv>0) then
! Allocate receive buffer
  allocate(recvbuffer(recvsize),recvstatus(MPI_STATUS_SIZE,nrecv), &
           recvrequest(nrecv))
  recvrequest=MPI_REQUEST_NULL
  ibuf_recv=1
  irecv=0

! Receive
  do iigrid=1,igridstail; igrid=igrids(iigrid);
    if (refine(igrid,mype)) then
   do ic3=1,2
do ic2=1,2
do ic1=1,2
      ! Only one of the sides will be necessary,
      ! so we do the loop only over dimensions, instead of
      ! over dimensions and sizes as in the routines
      ! old_neighbors and already_fine.
      do idims=1,3
        ipe_neighbor=fine_neighbors(ic1,ic2,ic3,idims,igrid)%ipe
        ineighbor   =fine_neighbors(ic1,ic2,ic3,idims,igrid)%igrid
        if (ineighbor.gt.0.and.ipe_neighbor.ne.mype) then

         if (idims.eq.1) iside=ic1
if (idims.eq.2) iside=ic2
if (idims.eq.3) iside=ic3
          !!! Check indices
          i1=kr(1,idims)*(2*iside-3);i2=kr(2,idims)*(2*iside-3)
          i3=kr(3,idims)*(2*iside-3);
   if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle

          inc1=ic1+i1;inc2=ic2+i2;inc3=ic3+i3;

          irecv=irecv+1
          itag=4**3*(igrid-1)+inc1*4**(1-1)+inc2*4**(2-1)+inc3*4**(3-1)

        !  print *, 'recv tag',itag

          call MPI_IRECV(recvbuffer(ibuf_recv),isize(idims), &
                         MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                         icomm,recvrequest(irecv),ierrmpi)
          ibuf_recv=ibuf_recv+isize(idims)

        end if
      end do
   end do
end do
end do
    end if
  end do

end if

if (nsend>0) then
! Allocate send buffer
  allocate(sendbuffer(sendsize),sendstatus(MPI_STATUS_SIZE,nsend),&
     sendrequest(nsend))
  sendrequest=MPI_REQUEST_NULL
  isend=0
  ibuf_send=1
! Send
do iigrid=1,igridstail; igrid=igrids(iigrid);
  ! Check whether it is necessary to store any block face, i.e.
  ! if any coarser neighbour is going to be refined 
 do iside=1,2
    i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3);i3=kr(3,1)*(2*iside-3);

    ! When there is a pole, faces are always zero and this is not necessary
    if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle

    if (neighbor_type(i1,i2,i3,igrid).eq.2) then
      ineighbor   =neighbor(1,i1,i2,i3,igrid)
      ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
      if (refine(ineighbor,ipe_neighbor)) then
        if (ipe_neighbor.ne.mype) then

          ic1=1+modulo(node(pig1_,igrid)-1,2)
          ic2=1+modulo(node(pig2_,igrid)-1,2)
          ic3=1+modulo(node(pig3_,igrid)-1,2);
          inc1=-2*i1+ic1;inc2=ic2;inc3=ic3;
          itag=4**3*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)+inc3*4**(3-1)

         ! print *, 'send tad', itag

          isend=isend+1

          ibuf_send_next=ibuf_send+isize(1)

          sendbuffer(ibuf_send:ibuf_send_next-1)=&
          reshape(pface(iside,1,igrid)%face,(/isize(1)/))

          call MPI_ISEND(sendbuffer(ibuf_send),isize(1), &
                         MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                         icomm,sendrequest(isend),ierrmpi)

          ibuf_send=ibuf_send_next

        end if
      end if
    end if
  end do
do iside=1,2
    i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3);i3=kr(3,2)*(2*iside-3);

    ! When there is a pole, faces are always zero and this is not necessary
    if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle

    if (neighbor_type(i1,i2,i3,igrid).eq.2) then
      ineighbor   =neighbor(1,i1,i2,i3,igrid)
      ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
      if (refine(ineighbor,ipe_neighbor)) then
        if (ipe_neighbor.ne.mype) then

          ic1=1+modulo(node(pig1_,igrid)-1,2)
          ic2=1+modulo(node(pig2_,igrid)-1,2)
          ic3=1+modulo(node(pig3_,igrid)-1,2);
          inc1=ic1;inc2=-2*i2+ic2;inc3=ic3;
          itag=4**3*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)+inc3*4**(3-1)

         ! print *, 'send tad', itag

          isend=isend+1

          ibuf_send_next=ibuf_send+isize(2)

          sendbuffer(ibuf_send:ibuf_send_next-1)=&
          reshape(pface(iside,2,igrid)%face,(/isize(2)/))

          call MPI_ISEND(sendbuffer(ibuf_send),isize(2), &
                         MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                         icomm,sendrequest(isend),ierrmpi)

          ibuf_send=ibuf_send_next

        end if
      end if
    end if
  end do
do iside=1,2
    i1=kr(1,3)*(2*iside-3);i2=kr(2,3)*(2*iside-3);i3=kr(3,3)*(2*iside-3);

    ! When there is a pole, faces are always zero and this is not necessary
    if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle

    if (neighbor_type(i1,i2,i3,igrid).eq.2) then
      ineighbor   =neighbor(1,i1,i2,i3,igrid)
      ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
      if (refine(ineighbor,ipe_neighbor)) then
        if (ipe_neighbor.ne.mype) then

          ic1=1+modulo(node(pig1_,igrid)-1,2)
          ic2=1+modulo(node(pig2_,igrid)-1,2)
          ic3=1+modulo(node(pig3_,igrid)-1,2);
          inc1=ic1;inc2=ic2;inc3=-2*i3+ic3;
          itag=4**3*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)+inc3*4**(3-1)

         ! print *, 'send tad', itag

          isend=isend+1

          ibuf_send_next=ibuf_send+isize(3)

          sendbuffer(ibuf_send:ibuf_send_next-1)=&
          reshape(pface(iside,3,igrid)%face,(/isize(3)/))

          call MPI_ISEND(sendbuffer(ibuf_send),isize(3), &
                         MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                         icomm,sendrequest(isend),ierrmpi)

          ibuf_send=ibuf_send_next

        end if
      end if
    end if
  end do
end do

end if

! Waitalls
if (nrecv>0) then
   call MPI_WAITALL(nrecv,recvrequest,recvstatus,ierrmpi)
   deallocate(recvstatus,recvrequest)
   ibuf_recv=1
end if

if (nsend>0) then
   call MPI_WAITALL(nsend,sendrequest,sendstatus,ierrmpi)
   deallocate(sendbuffer,sendstatus,sendrequest)
end if

end subroutine comm_faces
!==================================================================
subroutine end_comm_faces
use mod_fix_conserve, only: nrecv, recvbuffer
use mod_amrvacdef

!-----------------------------------------------------------------

! Deallocate receive buffer
if (nrecv>0) deallocate(recvbuffer)


end subroutine end_comm_faces
!==================================================================
subroutine deallocateBfaces

use mod_amrvacdef

integer :: igrid, iigrid, iside
!-----------------------------------------------------------------

do iigrid=1,igridstail; igrid=igrids(iigrid);
 do iside=1,2
    if (associated(pface(iside,1,igrid)%face)) then
      deallocate(pface(iside,1,igrid)%face)
     ! print *, 'deallocated iside, idims, igrid',iside,1,igrid
    end if
  end do
do iside=1,2
    if (associated(pface(iside,2,igrid)%face)) then
      deallocate(pface(iside,2,igrid)%face)
     ! print *, 'deallocated iside, idims, igrid',iside,2,igrid
    end if
  end do
do iside=1,2
    if (associated(pface(iside,3,igrid)%face)) then
      deallocate(pface(iside,3,igrid)%face)
     ! print *, 'deallocated iside, idims, igrid',iside,3,igrid
    end if
  end do
end do

end subroutine deallocateBfaces
!==================================================================
subroutine old_neighbors(child_igrid,child_ipe,igrid,ipe)

use mod_amrvacdef

integer, dimension(2,2,2), intent(in) :: child_igrid, child_ipe
integer, intent(in) :: igrid, ipe
! ... local ...
integer :: iside, i1,i2,i3, ic1,ic2,ic3
!------------------------------------------------------------------

do ic3=1,2
do ic2=1,2
do ic1=1,2
  old_neighbor(:,:,:,:,child_igrid(ic1,ic2,ic3))=-1
 do iside=1,2
    if (ic1.eq.iside) then
      i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3);i3=kr(3,1)*(2*iside-3);
      old_neighbor(1,i1,i2,i3,child_igrid(ic1,ic2,ic3))=&
         fine_neighbors(ic1,ic2,ic3,1,igrid)%igrid
      old_neighbor(2,i1,i2,i3,child_igrid(ic1,ic2,ic3))=&
         fine_neighbors(ic1,ic2,ic3,1,igrid)%ipe
    end if
  end do
do iside=1,2
    if (ic2.eq.iside) then
      i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3);i3=kr(3,2)*(2*iside-3);
      old_neighbor(1,i1,i2,i3,child_igrid(ic1,ic2,ic3))=&
         fine_neighbors(ic1,ic2,ic3,2,igrid)%igrid
      old_neighbor(2,i1,i2,i3,child_igrid(ic1,ic2,ic3))=&
         fine_neighbors(ic1,ic2,ic3,2,igrid)%ipe
    end if
  end do
do iside=1,2
    if (ic3.eq.iside) then
      i1=kr(1,3)*(2*iside-3);i2=kr(2,3)*(2*iside-3);i3=kr(3,3)*(2*iside-3);
      old_neighbor(1,i1,i2,i3,child_igrid(ic1,ic2,ic3))=&
         fine_neighbors(ic1,ic2,ic3,3,igrid)%igrid
      old_neighbor(2,i1,i2,i3,child_igrid(ic1,ic2,ic3))=&
         fine_neighbors(ic1,ic2,ic3,3,igrid)%ipe
    end if
  end do
end do
end do
end do

end subroutine old_neighbors
!==================================================================
subroutine already_fine(sFi,ichild,fine_min1,fine_min2,fine_min3,fine_max1,&
   fine_max2,fine_max3)
! This routine fills the fine faces before prolonging.
! It is the face equivalent of fix_conserve 
use mod_forest !!!!!! 
use mod_fix_conserve, only: nrecv,isize,recvbuffer
use mod_amrvacdef

type(tree_node_ptr) :: tree !!!!!!!
type(state) :: sFi
integer, intent(in) :: ichild
logical :: fine_min1,fine_min2,fine_min3,fine_max1,fine_max2,fine_max3

integer :: ineighbor,ipe_neighbor,ibufnext
integer :: iside,iotherside,i1,i2,i3,nx1,nx2,nx3
!-----------------------------------------------------------------
! Size of the block face
nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;nx3=ixMhi3-ixMlo3+1;

! Initialise everything to zero and false
fine_min1=.false.;fine_min2=.false.;fine_min3=.false.;
fine_max1=.false.;fine_max2=.false.;fine_max3=.false.;
sFi%ws%w=zero

!if (nrecv>0) ibuf=1


! This face communication is not needed in 1D
do iside=1,2
  i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3);i3=kr(3,1)*(2*iside-3);


  ! This is not necessary at the pole.
  ! We are inside a loop over the children, so the grid index is ichild
  if (neighbor_pole(i1,i2,i3,ichild)/=0) cycle


! Get old ipe and igrid of neighbour from array fake_neighbor
! Then plug it into the structure pfaces and get the faces 
  ineighbor   =old_neighbor(1,i1,i2,i3,ichild)
  ipe_neighbor=old_neighbor(2,i1,i2,i3,ichild)

  iotherside=3-iside
  if (ineighbor.gt.0) then
    if (iside.eq.1) then ! left side
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! if (mype.eq.0) then 
     !   tree%node => igrid_to_node(ichild,0)%node
 !print *, 'iside:',iside,'recv. from coarse at: ',tree%node%ig1,tree%node%ig2,tree%node%ig3
     ! end if
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      if (ipe_neighbor.eq.mype) then
        sFi%ws%w(ixMlo1-1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,1)=pface(iotherside,1,&
           ineighbor)%face(1,1:nx2,1:nx3)
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!1111
 !if (mype.eq.0) print *, pface(iotherside,1,ineighbor)%face(1,1:nx2,1:nx3)
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      else
     ! print *, 'ibuf',ibuf_recv
        ibufnext=ibuf_recv+isize(1)
        sFi%ws%w(ixMlo1-1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,1)=reshape(&
             source=recvbuffer(ibuf_recv:ibufnext-1),&
             shape=shape(sFi%ws%w(ixMlo1-1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,1)))
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!1111
     ! if (mype.eq.0) print *, recvbuffer(ibuf_recv:ibufnext-1)
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ibuf_recv=ibufnext
      end if
      fine_min1=.true.
    else ! right side
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! if (mype.eq.0) then 
     !   tree%node => igrid_to_node(ichild,0)%node
 !print *, 'iside:',iside,'recv. from coarse at: ',tree%node%ig1,tree%node%ig2,tree%node%ig3
     ! end if
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      if (ipe_neighbor.eq.mype) then
        sFi%ws%w(ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,1)=pface(iotherside,1,&
           ineighbor)%face(1,1:nx2,1:nx3)
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!1111
 !if (mype.eq.0) print *, pface(iotherside,1,ineighbor)%face(1,1:nx2,1:nx3)
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      else
     ! print *, 'ibuf',ibuf_recv
        ibufnext=ibuf_recv+isize(1)
        sFi%ws%w(ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,1)=reshape(&
             source=recvbuffer(ibuf_recv:ibufnext-1),&
             shape=shape(sFi%ws%w(ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,1)))
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!1111
     ! if (mype.eq.0) print *, recvbuffer(ibuf_recv:ibufnext-1)
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ibuf_recv=ibufnext
      end if
      fine_max1=.true.
    end if
  end if

end do
do iside=1,2
  i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3);i3=kr(3,2)*(2*iside-3);


  ! This is not necessary at the pole.
  ! We are inside a loop over the children, so the grid index is ichild
  if (neighbor_pole(i1,i2,i3,ichild)/=0) cycle


! Get old ipe and igrid of neighbour from array fake_neighbor
! Then plug it into the structure pfaces and get the faces 
  ineighbor   =old_neighbor(1,i1,i2,i3,ichild)
  ipe_neighbor=old_neighbor(2,i1,i2,i3,ichild)

  iotherside=3-iside
  if (ineighbor.gt.0) then
    if (iside.eq.1) then ! left side
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! if (mype.eq.0) then 
     !   tree%node => igrid_to_node(ichild,0)%node
 !print *, 'iside:',iside,'recv. from coarse at: ',tree%node%ig1,tree%node%ig2,tree%node%ig3
     ! end if
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      if (ipe_neighbor.eq.mype) then
        sFi%ws%w(ixMlo1:ixMhi1,ixMlo2-1,ixMlo3:ixMhi3,2)=pface(iotherside,2,&
           ineighbor)%face(1:nx1,1,1:nx3)
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!1111
 !if (mype.eq.0) print *, pface(iotherside,2,ineighbor)%face(1:nx1,1,1:nx3)
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      else
     ! print *, 'ibuf',ibuf_recv
        ibufnext=ibuf_recv+isize(2)
        sFi%ws%w(ixMlo1:ixMhi1,ixMlo2-1,ixMlo3:ixMhi3,2)=reshape(&
             source=recvbuffer(ibuf_recv:ibufnext-1),&
             shape=shape(sFi%ws%w(ixMlo1:ixMhi1,ixMlo2-1,ixMlo3:ixMhi3,2)))
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!1111
     ! if (mype.eq.0) print *, recvbuffer(ibuf_recv:ibufnext-1)
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ibuf_recv=ibufnext
      end if
      fine_min2=.true.
    else ! right side
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! if (mype.eq.0) then 
     !   tree%node => igrid_to_node(ichild,0)%node
 !print *, 'iside:',iside,'recv. from coarse at: ',tree%node%ig1,tree%node%ig2,tree%node%ig3
     ! end if
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      if (ipe_neighbor.eq.mype) then
        sFi%ws%w(ixMlo1:ixMhi1,ixMhi2,ixMlo3:ixMhi3,2)=pface(iotherside,2,&
           ineighbor)%face(1:nx1,1,1:nx3)
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!1111
 !if (mype.eq.0) print *, pface(iotherside,2,ineighbor)%face(1:nx1,1,1:nx3)
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      else
     ! print *, 'ibuf',ibuf_recv
        ibufnext=ibuf_recv+isize(2)
        sFi%ws%w(ixMlo1:ixMhi1,ixMhi2,ixMlo3:ixMhi3,2)=reshape(&
             source=recvbuffer(ibuf_recv:ibufnext-1),&
             shape=shape(sFi%ws%w(ixMlo1:ixMhi1,ixMhi2,ixMlo3:ixMhi3,2)))
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!1111
     ! if (mype.eq.0) print *, recvbuffer(ibuf_recv:ibufnext-1)
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ibuf_recv=ibufnext
      end if
      fine_max2=.true.
    end if
  end if

end do
do iside=1,2
  i1=kr(1,3)*(2*iside-3);i2=kr(2,3)*(2*iside-3);i3=kr(3,3)*(2*iside-3);


  ! This is not necessary at the pole.
  ! We are inside a loop over the children, so the grid index is ichild
  if (neighbor_pole(i1,i2,i3,ichild)/=0) cycle


! Get old ipe and igrid of neighbour from array fake_neighbor
! Then plug it into the structure pfaces and get the faces 
  ineighbor   =old_neighbor(1,i1,i2,i3,ichild)
  ipe_neighbor=old_neighbor(2,i1,i2,i3,ichild)

  iotherside=3-iside
  if (ineighbor.gt.0) then
    if (iside.eq.1) then ! left side
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! if (mype.eq.0) then 
     !   tree%node => igrid_to_node(ichild,0)%node
 !print *, 'iside:',iside,'recv. from coarse at: ',tree%node%ig1,tree%node%ig2,tree%node%ig3
     ! end if
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      if (ipe_neighbor.eq.mype) then
        sFi%ws%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3-1,3)=pface(iotherside,3,&
           ineighbor)%face(1:nx1,1:nx2,1)
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!1111
 !if (mype.eq.0) print *, pface(iotherside,3,ineighbor)%face(1:nx1,1:nx2,1)
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      else
     ! print *, 'ibuf',ibuf_recv
        ibufnext=ibuf_recv+isize(3)
        sFi%ws%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3-1,3)=reshape(&
             source=recvbuffer(ibuf_recv:ibufnext-1),&
             shape=shape(sFi%ws%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3-1,3)))
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!1111
     ! if (mype.eq.0) print *, recvbuffer(ibuf_recv:ibufnext-1)
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ibuf_recv=ibufnext
      end if
      fine_min3=.true.
    else ! right side
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! if (mype.eq.0) then 
     !   tree%node => igrid_to_node(ichild,0)%node
 !print *, 'iside:',iside,'recv. from coarse at: ',tree%node%ig1,tree%node%ig2,tree%node%ig3
     ! end if
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      if (ipe_neighbor.eq.mype) then
        sFi%ws%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMhi3,3)=pface(iotherside,3,&
           ineighbor)%face(1:nx1,1:nx2,1)
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!1111
 !if (mype.eq.0) print *, pface(iotherside,3,ineighbor)%face(1:nx1,1:nx2,1)
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      else
     ! print *, 'ibuf',ibuf_recv
        ibufnext=ibuf_recv+isize(3)
        sFi%ws%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMhi3,3)=reshape(&
             source=recvbuffer(ibuf_recv:ibufnext-1),&
             shape=shape(sFi%ws%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMhi3,3)))
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!1111
     ! if (mype.eq.0) print *, recvbuffer(ibuf_recv:ibufnext-1)
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ibuf_recv=ibufnext
      end if
      fine_max3=.true.
    end if
  end if

end do

end subroutine already_fine

end module mod_amr_fct

