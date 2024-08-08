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
subroutine get_divb(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
   ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,divb)

  ! Calculate div B within ixO

  use mod_amrvacdef

  integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
  double precision, intent(in)       :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:nw)
  double precision, intent(out)      :: divb(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3)
  ! .. local ..
  double precision                   :: bvec(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:ndir)

  integer                            :: ixCmin1,ixCmin2,ixCmin3,ixCmax1,&
     ixCmax2,ixCmax3, idir, idim, ixJpmin1,ixJpmin2,ixJpmin3,ixJpmax1,&
     ixJpmax2,ixJpmax3, ic1,ic2,ic3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
  integer                            :: ixKpmin1,ixKpmin2,ixKpmin3,ixKpmax1,&
     ixKpmax2,ixKpmax3, ixJpKpmin1,ixJpKpmin2,ixJpKpmin3,ixJpKpmax1,&
     ixJpKpmax2,ixJpKpmax3, ixJmmin1,ixJmmin2,ixJmmin3,ixJmmax1,ixJmmax2,&
     ixJmmax3, ixJmKmmin1,ixJmKmmin2,ixJmKmmin3,ixJmKmmax1,ixJmKmmax2,&
     ixJmKmmax3
  double precision                   :: divb_corner(ixImin1:ixImax1,&
     ixImin2:ixImax2,ixImin3:ixImax3), sign
  double precision                   :: aux_vol(ixImin1:ixImax1,&
     ixImin2:ixImax2,ixImin3:ixImax3)
  !-----------------------------------------------------------------------------

  bvec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndir)&
     =w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,b0_+1:b0_+ndir)

  if (typeemf .eq. 'none') then
     ! We don't enforce divB to machine precision:
     select case(typediv)
     case("central")
        call divvector(bvec,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb)
     case("limited")
        call divvectorS(bvec,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb)
     end select
     
  else
     ! Use the FCT (larger stencil) formula for divB:

     ! For fct, we calculate the divB on the corners according to Toth (2000), 
     ! eq. (27) and average to the cell centers for output.

     ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
     ixCmin1=ixOmin1-1;ixCmin2=ixOmin2-1;ixCmin3=ixOmin3-1; !Extend range by one


     ! Get the corner-centered divb:

     divb_corner(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = zero
     aux_vol(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = zero
     do idir = 1, ndim ! idir is the component of the field to consider (j)
        do ic3=0,1
        do ic2=0,1
        do ic1=0,1
        ixmin1=ixCmin1+ic1;ixmin2=ixCmin2+ic2;ixmin3=ixCmin3+ic3
        ixmax1=ixCmax1+ic1;ixmax2=ixCmax2+ic2;ixmax3=ixCmax3+ic3;

        select case(idir)
              
                case(1)
           sign = dble(ic1*2 - 1)
           
                 
                case(2)
           sign = dble(ic2*2 - 1)
           
                 
                case(3)
           sign = dble(ic3*2 - 1)
           
        end select

        if (slab) then
           divb_corner(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
              = divb_corner(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
              ixCmin3:ixCmax3) + sign * bvec(ixmin1:ixmax1,ixmin2:ixmax2,&
              ixmin3:ixmax3,idir)/dxlevel(idir)
        else
           divb_corner(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
              = divb_corner(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
              ixCmin3:ixCmax3) + sign * mygeo%dvolume(ixmin1:ixmax1,&
              ixmin2:ixmax2,ixmin3:ixmax3) * bvec(ixmin1:ixmax1,ixmin2:ixmax2,&
              ixmin3:ixmax3,idir)/dxlevel(idir)
           aux_vol(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
              = aux_vol(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) + &
              mygeo%dvolume(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)
        end if
        end do
        end do
        end do
     end do

     divb_corner(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
        = 2.0d0 * ndim * divb_corner(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
        ixCmin3:ixCmax3) / aux_vol(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
        ixCmin3:ixCmax3)

     if (slab) divb_corner(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
        = divb_corner(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) &
        / 2.0d0**(ndim-1)



     ! Now average back to the cell centers:

     divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = zero
     do ic3=-1,0
     do ic2=-1,0
     do ic1=-1,0
     ixCmin1=ixOmin1+ic1;ixCmin2=ixOmin2+ic2;ixCmin3=ixOmin3+ic3
     ixCmax1=ixOmax1+ic1;ixCmax2=ixOmax2+ic2;ixCmax3=ixOmax3+ic3;

     divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
        = divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + &
        divb_corner(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)

     end do
     end do
     end do
     divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
        = divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
        / 2.0d0**(ndim)

  end if

end subroutine get_divb
!=============================================================================
