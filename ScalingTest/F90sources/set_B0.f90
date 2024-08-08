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
subroutine set_B0_grid(igrid)

use mod_amrvacdef

integer, intent(in) :: igrid
!-----------------------------------------------------------------------------

call set_B0_cell(pB0_cell(igrid)%w,px(igrid)%x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,&
   ixGhi2,ixGhi3)
call set_B0_face(igrid,px(igrid)%x)

end subroutine set_B0_grid
!=============================================================================
subroutine set_B0_cell(wB0,x,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3)

use mod_amrvacdef

integer, intent(in):: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
double precision, intent(inout) :: wB0(ixmin1:ixmax1,ixmin2:ixmax2,&
   ixmin3:ixmax3,1:ndir)
double precision, intent(in) :: x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
   1:ndim)
!-----------------------------------------------------------------------------
wB0(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1:ndir)=zero

! approximate cell-averaged B0 as cell-centered B0
select case (typeaxial)
case ("spherical")
   
   if (abs(Bdip)>smalldouble) then
      wB0(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)=2.0d0*Bdip*cos(x&
         (ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2))/x(ixmin1:ixmax1,&
         ixmin2:ixmax2,ixmin3:ixmax3,1)**3.0d0
      wB0(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2)=Bdip*sin(x&
         (ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2))/x(ixmin1:ixmax1,&
         ixmin2:ixmax2,ixmin3:ixmax3,1)**3.0d0
   end if

   if (abs(Bquad)>smalldouble) then
      wB0(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)=wB0(ixmin1:ixmax1,&
         ixmin2:ixmax2,ixmin3:ixmax3,1) &
           +Bquad*0.5d0*(1.0d0+3.0d0*cos(2.0d0*x(ixmin1:ixmax1,ixmin2:ixmax2,&
              ixmin3:ixmax3,2)))/x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
              1)**4
      wB0(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2)=wB0(ixmin1:ixmax1,&
         ixmin2:ixmax2,ixmin3:ixmax3,2)+Bquad*sin(2.0d0*x(ixmin1:ixmax1,&
         ixmin2:ixmax2,ixmin3:ixmax3,2))/x(ixmin1:ixmax1,ixmin2:ixmax2,&
         ixmin3:ixmax3,1)**4
   end if
   if (abs(Boct)>smalldouble) then
      wB0(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)=wB0(ixmin1:ixmax1,&
         ixmin2:ixmax2,ixmin3:ixmax3,1) &
                   +Boct*(10.0d0*cos(2.0d0*x(ixmin1:ixmax1,ixmin2:ixmax2,&
                      ixmin3:ixmax3,2))-2.0d0) &
                        *cos(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2))&
                           /x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)**5
      wB0(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2)=wB0(ixmin1:ixmax1,&
         ixmin2:ixmax2,ixmin3:ixmax3,2) &
                   +Boct*1.5d0*(3.0d0+5.0d0*cos(2.0d0*x(ixmin1:ixmax1,&
                      ixmin2:ixmax2,ixmin3:ixmax3,2))) &
                        *sin(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2))&
                           /x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)**5
   end if
  
end select

if(dabs(Busr)/=zero) then
   call specialset_B0(ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,ixmin1,ixmin2,&
      ixmin3,ixmax1,ixmax2,ixmax3,x,wB0)
end if

end subroutine set_B0_cell
!=============================================================================
subroutine set_B0_face(igrid,x)

use mod_amrvacdef

integer, intent(in) :: igrid
double precision, intent(in) :: x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
   1:ndim)

double precision :: xC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:ndim),dx1,&
   dx2,dx3,xmin1,xmin2,xmin3,xshift1,xshift2,xshift3
integer :: idims, ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, ix, idims2
!-----------------------------------------------------------------------------
dx1=rnode(rpdx1_,igrid);dx2=rnode(rpdx2_,igrid);dx3=rnode(rpdx3_,igrid);
xmin1=rnode(rpxmin1_,igrid);xmin2=rnode(rpxmin2_,igrid)
xmin3=rnode(rpxmin3_,igrid);

do idims=1,ndim
   ixCmin1=ixMlo1-kr(1,idims);ixCmin2=ixMlo2-kr(2,idims)
   ixCmin3=ixMlo3-kr(3,idims); ixCmax1=ixMhi1;ixCmax2=ixMhi2;ixCmax3=ixMhi3;
   xshift1=half*(one-kr(1,idims));xshift2=half*(one-kr(2,idims))
   xshift3=half*(one-kr(3,idims));
   do idims2=1,ndim
      select case(idims2)
      case(1)
        do ix = ixCmin1,ixCmax1
          xC(ix,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)=xmin1+(dble&
             (ix-dixB)-xshift1)*dx1
        end do
      case(2)
        do ix = ixCmin2,ixCmax2
          xC(ixCmin1:ixCmax1,ix,ixCmin3:ixCmax3,2)=xmin2+(dble&
             (ix-dixB)-xshift2)*dx2
        end do
      case(3)
        do ix = ixCmin3,ixCmax3
          xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ix,3)=xmin3+(dble&
             (ix-dixB)-xshift3)*dx3
        end do
      end select
   end do

   select case(idims)
   case (1)
      call set_B0_cell(pB0_face1(igrid)%w,xC,ixCmin1,ixCmin2,ixCmin3,ixCmax1,&
         ixCmax2,ixCmax3) 
   case (2)
      call set_B0_cell(pB0_face2(igrid)%w,xC,ixCmin1,ixCmin2,ixCmin3,ixCmax1,&
         ixCmax2,ixCmax3) 
   case (3)
      call set_B0_cell(pB0_face3(igrid)%w,xC,ixCmin1,ixCmin2,ixCmin3,ixCmax1,&
         ixCmax2,ixCmax3) 
   end select
end do

end subroutine set_B0_face
!=============================================================================
subroutine alloc_B0_grid(igrid)

use mod_amrvacdef

integer, intent(in) :: igrid

integer :: ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3
!-----------------------------------------------------------------------------

allocate(pB0_cell(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:ndir))
ixCmin1=ixMlo1-kr(1,1);ixCmin2=ixMlo2-kr(2,1);ixCmin3=ixMlo3-kr(3,1)
ixCmax1=ixMhi1;ixCmax2=ixMhi2;ixCmax3=ixMhi3;
allocate(pB0_face1(igrid)%w(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
   1:ndir))
ixCmin1=ixMlo1-kr(1,2);ixCmin2=ixMlo2-kr(2,2);ixCmin3=ixMlo3-kr(3,2)
ixCmax1=ixMhi1;ixCmax2=ixMhi2;ixCmax3=ixMhi3;
allocate(pB0_face2(igrid)%w(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
   1:ndir))
ixCmin1=ixMlo1-kr(1,3);ixCmin2=ixMlo2-kr(2,3);ixCmin3=ixMlo3-kr(3,3)
ixCmax1=ixMhi1;ixCmax2=ixMhi2;ixCmax3=ixMhi3;
allocate(pB0_face3(igrid)%w(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
   1:ndir))

end subroutine alloc_B0_grid
!=============================================================================
subroutine dealloc_B0_grid(igrid)

use mod_amrvacdef

integer, intent(in) :: igrid
!-----------------------------------------------------------------------------

deallocate(pB0_cell(igrid)%w)
deallocate(pB0_face1(igrid)%w)
deallocate(pB0_face2(igrid)%w)
deallocate(pB0_face3(igrid)%w)

end subroutine dealloc_B0_grid
!=============================================================================
