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
! mod_physaux to perform some useful computations, e.g. obtain b2, the square 
! of the co-moving magnetic field.
! You could also make a subroutine to calculate the current and other things...
! Oliver Porth, 2016-01-27
!=============================================================================
module mod_physaux
  implicit none

contains

  !=============================================================================
  subroutine get_b2(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,b2,w_is_primitive)
    ! Obtains b**2, the invariant co-moving magnetic field.

    use mod_metric, only: lower3
    use mod_amrvacdef 

    integer, intent(in)                                    :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw), intent(in)      :: w
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,ndim), intent(in)    :: x
    logical, optional, intent(in)                          :: w_is_primitive
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), intent(out)        :: b2
    ! .. local ..
    logical                                                :: is_primitive
    double precision                                       :: &
       bD(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:3)
    double precision                                       :: &
       tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        BdotV(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    !-----------------------------------------------------------------------------

    if(.not.present(w_is_primitive)) then
       is_primitive = .false. !By default assuming conserved variables coming in
    else
       is_primitive = w_is_primitive
    end if

    call lower3(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,myM,w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,b1_:b3_),bD)
    if (is_primitive) then
       BdotV(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
          = bD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          1)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          u1_)+bD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          2)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          u2_)+bD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          3)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,u3_)
       if (useprimitiveRel) BdotV(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = BdotV(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          lfac_)
    else
       BdotV(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
          = (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          s1_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          b1_)+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          s2_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          b2_)+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          s3_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b3_))&
          /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_)
    end if

    tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,b1_)*bD(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,b2_)*bD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,2)+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       b3_)*bD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)

    b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,lfac_)**2 + BdotV(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)**2

  end subroutine get_b2
  !=============================================================================
  subroutine get_u4(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,u4,w_is_primitive)
    ! Obtains the contravariant four-velocity u4
    ! from conservative or primitive variables

    use mod_metric, only: raise3, square3u
    use mod_amrvacdef 

    integer, intent(in)                                    :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw), intent(in)      :: w
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,ndim), intent(in)    :: x
    logical, optional, intent(in)                          :: w_is_primitive
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:ndir), intent(out) :: u4
    ! .. local ..
    logical                                                :: is_primitive
    double precision                                       :: &
       v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)                     :: a, b, c, d
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)                     :: VdotB, sqrB
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)              :: bD, sU
    !-----------------------------------------------------------------------------

    if(.not.present(w_is_primitive)) then
       is_primitive = .false. !By default assuming conserved variables coming in
    else
       is_primitive = w_is_primitive
    end if

    if (is_primitive) then
        u4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) &
           = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           u1_) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)&
           /myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
           *myM%beta(1)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) 
         u4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) &
            = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            u2_) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)&
            /myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
            *myM%beta(2)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3) 
         u4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3) &
            = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            u3_) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)&
            /myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
            *myM%beta(3)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)        
    else
       VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
          = (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          s1_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          b1_)+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          s2_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          b2_)+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          s3_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b3_))&
          /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_)
       call square3u(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
          ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,myM,w(ixImin1:ixImax1,&
          ixImin2:ixImax2,ixImin3:ixImax3,b1_:b3_),sqrB(ixImin1:ixImax1,&
          ixImin2:ixImax2,ixImin3:ixImax3))
       call raise3(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
          ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,myM,w(ixImin1:ixImax1,&
          ixImin2:ixImax2,ixImin3:ixImax3,s1_:s3_),sU)

       
          v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
             =(sU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             1)+VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*w&
             (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+1))/ &
            (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               xi_)+sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
       u4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) &
          = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_) * &
            ( v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) - &
               myM%beta(1)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)/myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) )
       
       
          v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
             =(sU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             2)+VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*w&
             (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+2))/ &
            (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               xi_)+sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
       u4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) &
          = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_) * &
            ( v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) - &
               myM%beta(2)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)/myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) )
       
       
          v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
             =(sU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             3)+VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*w&
             (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+3))/ &
            (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               xi_)+sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
       u4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3) &
          = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_) * &
            ( v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) - &
               myM%beta(3)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)/myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) )
       
    end if

    ! Obtain the time-component from u_alpha * u**alpha = -1

    a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = myM%g(0,&
       0)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 2.0d0*( &
       u4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*myM%g(1,&
       0)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)+ &
       u4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)*myM%g(2,&
       0)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)+ &
       u4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)*myM%g(3,&
       0)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
    c(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 1.0d0 + ( &
       u4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)&
       **2*myM%g(1,1)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)+ u4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2)**2*myM%g(2,2)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)+ u4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       3)**2*myM%g(3,3)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)) + 2.0d0 * ( myM%g(1,2)%elem(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) * u4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1) * u4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,2)  + myM%g(1,3)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*u4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1) * u4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3) + myM%g(2,&
       3)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*u4&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) * &
       u4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))
    d(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = b(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)**2 - 4.0d0*a(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)*c(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    u4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,0) = -half&
       /a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*(b(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)+sqrt(d(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)))
    
  end subroutine get_u4
  !=============================================================================
  subroutine get_b4(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,b4,w_is_primitive)
    ! Obtains the contravariant four-magnetic field in the fluid frame (little b) 
    ! from conservative or primitive variables
    !=============================================================================
    use mod_metric, only: lower3
    use mod_amrvacdef 

    integer, intent(in)                                    :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw), intent(in)      :: w
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,ndim), intent(in)    :: x
    logical, optional, intent(in)                          :: w_is_primitive
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:ndir), intent(out) :: b4
    ! .. local ..
    logical                                                :: is_primitive
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)              :: vd, bD
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:ndir)              :: u4
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)                     :: VdotB, sqrB
    !-----------------------------------------------------------------------------

    if(.not.present(w_is_primitive)) then
       is_primitive = .false. !By default assuming conserved variables coming in
    else
       is_primitive = w_is_primitive
    end if

    if (is_primitive) then
       
       call lower3(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
          ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,myM,w(ixImin1:ixImax1,&
          ixImin2:ixImax2,ixImin3:ixImax3,u1_:u3_),vd)
        vd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) &
           = vd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)&
           /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_) 
         vd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) &
            = vd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)&
            /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_) 
         vd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3) &
            = vd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)&
            /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_) 
       
    else
       
       VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
          = (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          s1_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          b1_)+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          s2_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          b2_)+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          s3_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b3_))&
          /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_)
       call lower3(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
          ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,myM,w(ixImin1:ixImax1,&
          ixImin2:ixImax2,ixImin3:ixImax3,b1_:b3_),bD)
       sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)  &
          = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          b1_)*bD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          1)+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          b2_)*bD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          2)+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          b3_)*bD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)

        
       vd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) &
          = (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          s1_)+VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*bd&
          (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1))/ &
       (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_)+sqrB&
          (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
       
         
       vd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) &
          = (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          s2_)+VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*bd&
          (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))/ &
       (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_)+sqrB&
          (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
       
         
       vd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3) &
          = (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          s3_)+VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*bd&
          (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))/ &
       (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_)+sqrB&
          (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
       
       
    end if ! is_primitive

    b4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,0) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_) * ( w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,b1_)*vd(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)+ w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,b2_)*vd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,2)+ w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       b3_)*vd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))&
       /myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    call get_u4(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,u4,w_is_primitive&
       =is_primitive)
     b4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) &
        = (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        b1_)+myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*b4&
        (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,0)*u4&
        (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1))&
        /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)
      b4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) &
         = (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         b2_)+myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*b4&
         (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,0)*u4&
         (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))&
         /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)
      b4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3) &
         = (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         b3_)+myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*b4&
         (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,0)*u4&
         (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))&
         /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)

    
  end subroutine get_b4
  !=============================================================================
  subroutine get_TEM4_ud(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,tem4,w_is_primitive)
    ! Obtains the electromagnetic part of the four-stress-energy tensor 
    ! One index up, one index down.  
    !=============================================================================
    use mod_metric, only: lower4
    use mod_amrvacdef 

    integer, intent(in)                                    :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw), intent(in)      :: w
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,ndim), intent(in)    :: x
    logical, optional, intent(in)                          :: w_is_primitive
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:ndir,0:ndir), intent(out) :: tem4
    ! .. local ..
    logical                                                :: is_primitive
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:ndir)              :: u4, b4, u4d, b4d
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)                     :: b2
    integer                                                :: imu, inu
    !=============================================================================

    if(.not.present(w_is_primitive)) then
       is_primitive = .false. !By default assuming conserved variables coming in
    else
       is_primitive = w_is_primitive
    end if

    call get_u4(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,u4,w_is_primitive&
       =is_primitive)
    call get_b4(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,b4,w_is_primitive&
       =is_primitive)
    call get_b2(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,b2,w_is_primitive&
       =is_primitive)

    call lower4(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,myM,u4,u4D)
    call lower4(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,myM,b4,b4D)

    do imu = 0, ndir
       do inu = 0, ndir
          tem4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,imu,inu) &
             = b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*u4&
             (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             imu)*u4d(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             inu) + half*b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)*kr(imu,inu) - b4(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,imu)*b4d(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,inu)
       end do
    end do

  end subroutine get_TEM4_ud
  !=============================================================================
  subroutine get_TPAKE4_ud(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,tpake4,&
     w_is_primitive)
    ! Obtains the free particle energy part of the four-stress-energy tensor 
    ! as defined by Mc Kinney et al. (2012), MNRAS 423, p. 3083, Eq. (6)
    ! One index up, one index down.  
    !=============================================================================
    use mod_metric, only: lower4
    use mod_amrvacdef 

    integer, intent(in)                                    :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw), intent(in)      :: w
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,ndim), intent(in)    :: x
    logical, optional, intent(in)                          :: w_is_primitive
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:ndir,0:ndir), intent(out) :: tpake4
    ! .. local ..
    logical                                                :: is_primitive
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:ndir)              :: u4, u4d
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)                     :: rho
    integer                                                :: imu, inu
    !=============================================================================

    if(.not.present(w_is_primitive)) then
       is_primitive = .false. !By default assuming conserved variables coming in
    else
       is_primitive = w_is_primitive
    end if

    if (is_primitive) then
       rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
          = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
    else
       rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
          = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,d_)&
          /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)
    end if

    call get_u4(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,u4,w_is_primitive&
       =is_primitive)
    call lower4(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,myM,u4,u4D)

    do imu = 0, ndir
       do inu = 0, ndir
          tpake4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,imu,inu) &
             = (u4d(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             inu)+one*kr(0,inu)) * rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)*u4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,imu)
       end do
    end do

  end subroutine get_TPAKE4_ud
  !=============================================================================
  subroutine get_TEN4_ud(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,ten4,w_is_primitive)
    ! Obtains the enthalpy part of the four-stress-energy tensor 
    ! as defined by Mc Kinney et al. (2012), MNRAS 423, p. 3083, Eq. (6)
    ! One index up, one index down.  
    !=============================================================================
    use mod_metric, only: lower4
    use mod_amrvacdef 

    integer, intent(in)                                    :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw), intent(in)      :: w
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,ndim), intent(in)    :: x
    logical, optional, intent(in)                          :: w_is_primitive
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:ndir,0:ndir), intent(out) :: ten4
    ! .. local ..
    logical                                                :: is_primitive
    double precision                                       :: &
       wprim(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw)
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:ndir)              :: u4, u4d
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)                     :: p, rhoh, e
    integer                                                :: imu, inu
    !=============================================================================

    if(.not.present(w_is_primitive)) then
       is_primitive = .false. !By default assuming conserved variables coming in
    else
       is_primitive = w_is_primitive
    end if

    call Pressure(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,.not.is_primitive,p)

    ! Most general case to use the enthalpy.  For this we need variables in primitive state:
    wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:nw) &
       = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:nw)
    if (.not.is_primitive) call primitiven(ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wprim,&
       patchfalse)
    call Enthalpy(wprim,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,patchfalse,rhoh)

    ! The thermal energy contribution is rhoh-rho, 
    ! e.g. in ideal-gas: e = rhoh-rho = u+p = govergminone*p
    if (is_primitive) then
       e(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
          rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) - &
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
    else
       e(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
          rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) - &
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,d_)&
          /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)
    end if

    call get_u4(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,u4,w_is_primitive&
       =is_primitive)
    call lower4(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,myM,u4,u4D)

    do imu = 0, ndir
       do inu = 0, ndir
          ten4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,imu,inu) &
             = e(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*u4&
             (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             imu)*u4d(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             inu) + kr(imu,inu) * p(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)
       end do
    end do

  end subroutine get_TEN4_ud
  !=============================================================================
end module mod_physaux
!=============================================================================
