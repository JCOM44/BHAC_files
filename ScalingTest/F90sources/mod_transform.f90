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
! Transformations between coordinate systems
! e.g. anything -> BL coordinates
! KS coordinates <-> "Cartesian" KS coordinates
!
!
! Oliver Porth
! 25.01.2015
!=============================================================================
module mod_transform
  implicit none

contains

  !=============================================================================
  subroutine KSToBL(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xKS,xBL)

    ! Transforms the KS coordinate to the BL coordinate
    ! Only transform space, time would be a separate routine.
    !
    ! rBL     = rKS
    ! thetaBL = thetaKS
    ! phiBL   = phiKS - 1/2 a 1/s ln((x-(1+s))/(x-(1-s)))
    ! s = sqrt(1-a**2)
    !
    ! The transformation yields infinite phi-angle at the outer horizon, so take care!  
    ! Adopted from a subroutine from Ziri Younsi, BHOSS code.
    ! 
    use mod_amrvacdef

    integer,intent(in)                                     :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: xKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(out)  :: xBL
    ! .. local ..
    
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)                     :: s, rmrp, rmrm
   
    !-----------------------------------------------------------------------------


    xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) &
       = xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)
    
    xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3) &
       = xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)
   

    
    s(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = sqrt(1.0d0 - &
       eqpar(a_)**2)
    rmrp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1) - (1.0d0 + s(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
    rmrm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1) - (1.0d0 - s(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
    
    xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) &
       = xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2) - half*eqpar(a_)/s(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) &
         * log(rmrp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
            /rmrm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
    
     where (xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        2).ne.xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))
         xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) = zero
     end where

    

  end subroutine KSToBL
  !=============================================================================
  subroutine BLToKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xBL,xKS)

    ! Transforms the BL coordinate to the KS coordinate
    ! Only transform space, time would be a separate routine.
    !
    ! rKS     = rBL
    ! thetaKS = thetaBL
    ! phiKS   = phiBL + 1/2 a 1/s ln((x-(1+s))/(x-(1-s)))
    ! s = sqrt(1-a**2)
    !
    ! The transformation yields infinite phi-angle at the outer horizon, so take care!  
    ! Adopted from a subroutine from Ziri Younsi, BHOSS code.
    !
    use mod_amrvacdef

    integer,intent(in)                                     :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: xBL
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(out)  :: xKS
    ! .. local ..
    
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)                     :: s, rmrp, rmrm
   
    !-----------------------------------------------------------------------------

    xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) &
       = xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)
    
    xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3) &
       = xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)
   
    
    s(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = sqrt(1.0d0 - &
       eqpar(a_)**2)
    rmrp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1) - (1.0d0 + s(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
    rmrm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1) - (1.0d0 - s(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
    
    xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) &
       = xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2) + half*eqpar(a_)/s(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) &
         * log(rmrp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
            /rmrm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))

    where (xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2).ne.xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))
        xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) = zero
    end where
   

  end subroutine BLToKS
  !=============================================================================
  subroutine u4KStoBL(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xKS,u4KS,u4BL,J)

    ! Transforms the (contravariant) four-velocity u4KS from Kerr-Schild coordinates
    ! to BL coordinates u4BL.
    !
    ! See also Font et al., Mon. Not. R. Astron. Soc. 305, 920±936 (1999) for the
    ! transformations.
    !
    ! utBL       = utBL - 2 M rKS / Delta * urKS
    ! uphiBL     = uphiBL - a/Delta * urKS
    ! urBL       = urKS
    ! uthetaBL   = uthetaKS
    !
    ! Jacobian:
    !
    !( U_BL^0     )      | 1       - 2 M rKS / Delta           0                  0            | ( U_KS^0 )
    !( U_BL^r     ) ____ | 0              1                    0                  0            | ( U_KS^r )
    !( U_BL^theta ) ____ | 0              0                    1                  0            | ( U_KS^theta )
    !( U_BL^phi   )      | 0          - a/Delta                0                  1            | ( U_KS^phi )
    !
    use mod_amrvacdef

    integer, intent(in)                                    :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: xKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(in)   :: u4KS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(out)  :: u4BL
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3), optional, intent(out)  :: J
    ! .. local ..
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)                     :: Delta, twoMr
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3)         :: Jac
    integer                                                :: ix, jx
    !-----------------------------------------------------------------------------

    twoMr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
       = 2.0d0*eqpar(m_)*abs(xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1))
    Delta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
       = xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)&
       **2 - twoMr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + &
       eqpar(a_)**2

    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:,:) = 0.0d0

    ! t column:
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,0,0)       = 1.0d0
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,0,1)       &
       = - twoMr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
       / Delta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    ! r column:
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1,1)       = 1.0d0
    ! theta column:
    
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3,3)     = 1.0d0
   
    
    ! phi column:
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2,1)    = - eqpar(a_)&
       /Delta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2,2) = 1.0d0
   

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Jacobian fully assembled
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! transform contravariant four-vector
    call matvec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,u4KS,Jac,u4BL)
    
    if (present(J)) J = Jac

  end subroutine u4KStoBL
  !=============================================================================
  subroutine u4BLtoKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xBL,u4BL,u4KS,J)

    ! Transforms the (contravariant) four-velocity u4BL from Boyer-Lindquist coordinates
    ! to KS coordinates u4KS.  Often initial conditions are
    ! given in terms of BL coordinates and this routine comes in handy.
    !
    ! See also Font et al., Mon. Not. R. Astron. Soc. 305, 920Â±936 (1999) for the
    ! transformations.
    !
    ! utKS       = utBL + 2 M rBL / Delta * urBL
    ! uphiKS     = uphiBL + a/Delta * urBL
    ! urKS       = urBL
    ! uthetaKS   = uthetaBL
    !
    use mod_amrvacdef

    integer, intent(in)                                    :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: xBL
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(in)   :: u4BL
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(out)  :: u4KS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3), optional, intent(out)  :: J
    ! .. local ..
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)                     :: Delta, twoMr
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3)         :: Jac
    !-----------------------------------------------------------------------------

    twoMr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
       = 2.0d0*eqpar(m_)*abs(xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1))
    Delta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
       = xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)&
       **2 - twoMr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + &
       eqpar(a_)**2


    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:,:) = 0.0d0
    
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,0,0) = 1.0d0
    
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,0,1) &
       = twoMr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
       /Delta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1,1) = 1.0d0

    
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2,2) = 1.0d0
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2,1) = eqpar(a_)&
       /Delta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
   
    
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3,3) = 1.0D0
   
!    u4KS(ixO^S,0) = u4BL(ixO^S,0) + twoMr(ixO^S)/Delta(ixO^S) * u4BL(ixO^S,1)
!    u4KS(ixO^S,1) = u4BL(ixO^S,1)
!    {^IFPHI
!    u4KS(ixO^S,^PHI) = u4BL(ixO^S,^PHI) + eqpar(a_)/Delta(ixO^S) * u4BL(ixO^S,1)
!    }
!    {^IFZ
!    u4KS(ixO^S,^Z) = u4BL(ixO^S,^Z)
!    }

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Jacobian fully assembled
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! transform contravariant four-vector
    call matvec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,u4BL,Jac,u4KS)
    
    if (present(J)) J = Jac
    
  end subroutine u4BLtoKS
  !=============================================================================
  subroutine u4KStoCKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xKS,u4KS,u4CKS,J)

    ! Alejandro Cruz-osorio 12.07.2017
    ! Transforms the (contravariant) four-velocity u4KS from Kerr-Schild spherical coordinates
    ! to the cartesian Kerr-Schild coordinates u4CKS.  
    !
    !( U^0 )      | 1          0                   0                 0  | ( U^0     )
    !( U^x ) ____ | 0  cos(phi)sin(theta)     x cot(theta)          -y  | ( U^r     )
    !( U^y ) ____ | 0  sin(phi)sin(theta)     y cot(theta)           x  | ( U^theta )
    !( U^z )      | 0     cos(theta)       -r_ks sint(theta)         0  | ( U^phi   )
    !
    use mod_amrvacdef

    integer, intent(in)                                    :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: xKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(in)   :: u4KS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(out)  :: u4CKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3), optional, intent(out)  :: J
    ! .. local ..
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3)               :: xCKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)                     :: sth, cth, sph, cph
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3)         :: Jac
    ! ... hectors test ...
    integer :: k,l
    !-----------------------------------------------------------------------------
    call KSToCKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xKS,xCKS)

    
    sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       sin(xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))
    cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       cos(xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))
   
    
    sph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       sin(xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))
    cph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       cos(xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))
   

    ! t:
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:,:) = 0.0d0
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,0,0) = 1.0d0

    ! x:
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1,1) &
       = cph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) * &
       sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    
 !Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1,3) = xcKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) * cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)/sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
!   Using expression below instead to avoid singularity at the pole
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1,3) &
       = (xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)*cph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) - &
       eqpar(a_)*sph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))*cth&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
   
    
    
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1,2) &
       = - xCKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) 
   
   
    
    ! y:
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2,1) &
       = sph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) * &
       sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    
    
 !Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2,3) = xCKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) * cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)/sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
!   Using expression below instead to avoid singularity at the pole
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2,3) &
       = (xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)*sph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + &
       eqpar(a_)*cph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))*cth&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
   
   
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2,2) &
       = xCKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)
   

    ! z:
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3,1) &
       = cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3,3) &
       = - xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1) * sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Jacobian fully assembled
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! transform contravariant four-vector
    call matvec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,u4KS,Jac,u4CKS)
    
    if (present(J)) J = Jac

  end subroutine u4KStoCKS
  !=============================================================================
  subroutine u4CKStoBL(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xCKS,u4CKS,u4BL,J)

    ! 20.07.2017 Oliver Porth
    ! Transforms the (contravariant) four-vector u4CKS from Cartesian Kerr-Schild
    ! to spherical Boyer-Lindquist

    use mod_amrvacdef

    integer, intent(in)                                    :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: xCKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(in)   :: u4CKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(out)  :: u4BL
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3), optional, intent(out)  :: J
    ! .. local ..
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3)               :: xKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3)               :: u4dummy1,u4dummy2
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3)         :: Jac1
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3)         :: Jac2
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3)         :: Jac
    !-----------------------------------------------------------------------------

    ! Get Jacobian KS->BL
    call CKSToKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xCKS,xKS)
    call u4KSToBL(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xKS,u4dummy1,u4dummy2,J=Jac1)
    
    ! Get Jacobian CKS->KS
    call u4CKSToKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xCKS,u4dummy1,u4dummy2,J=Jac2)

    ! Compose Jacobian CKS->BL:
    call compose(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,Jac1,Jac2,Jac)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Jacobian fully assembled
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! transform contravariant four-vector
    call matvec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,u4CKS,Jac,u4BL)

    if (present(J)) J = Jac

  end subroutine u4CKStoBL
  !=============================================================================
  subroutine u4CKStoKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xCKS,u4CKS,u4KS,J)

    ! 20.07.2017 Oliver Porth
    ! Transforms the (contravariant) four-vector u4CKS from Cartesian Kerr-Schild
    ! to spherical Kerr-Schild
    !
    ! Jacobian:
    !
    !( U^0     )      | 1         0                                   0                                   0                     | ( U^0 )
    !( U^r     ) ____ | 0  rks sth^2 x / chi                   rks sth^2 y / chi                (cth (x^2 + y^2))/chi           | ( U^x )
    !( U^theta ) ____ | 0  (cth sth x)/chi                     (cth sth y)/chi                  -((sth^2 (cph x + sph y))/chi)  | ( U^y )
    !( U^phi   )      | 0  -((rks sph sth^3 + cth^2 y)/chi)    (cph rks sth^3 + cth^2 x)/chi    (cth sth (-sph x + cph y))/chi  | ( U^z )
    ! 
    ! chi = cph rks sth^3 x + rks sph sth^3 y + cth^2 (x^2 + y^2)
    !
    use mod_amrvacdef

    integer, intent(in)                                    :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: xCKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(in)   :: u4CKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(out)  :: u4KS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3), optional, intent(out)  :: J
    ! .. local ..
    integer                                                :: ix, jx
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3)               :: xKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)                     :: sth, cth, sph, cph, chi
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3)         :: Jac
    !-----------------------------------------------------------------------------

    call CKSToKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xCKS,xKS)

    
    sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       sin(xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))
    cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       cos(xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))
   
    
    sph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       sin(xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))
    cph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       cos(xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))
   
    
    chi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
       **2*(xCKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)&
       **2 +xCKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)&
       **2)+  cph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*sth&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)**3*xCKS&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*xKS&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)  + &
       sph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*sth&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)**3*xCKS&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)*xKS&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)

    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:,:) = 0.0d0

    ! t:
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,0,0) = 1.0d0
    ! r:
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1,1) &
       = (sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
       **2*xCKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)*xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1))&
       /chi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1,2) &
       = (sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
       **2*xCKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2)*xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1))&
       /chi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
   
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1,3) &
       = (cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*(xCKS&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2  + &
       xCKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)**2))&
       /chi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    ! theta:
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3,1) &
       = (cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*sth&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*xCKS(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,1))/chi(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3,2) &
       = (cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*sth&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*xCKS(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))/chi(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)
   
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3,3) &
       = (-1.0d0*sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
       **2*(cph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*xCKS&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)   + &
       sph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*xCKS&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)))&
       /chi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    ! phi:
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2,1) &
       = (-1.0d0*( cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
       **2*xCKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2)+  sph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*sth&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)**3*xKS&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)))&
       /chi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2,2) &
       = (cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
       **2*xCKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)+  cph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*sth&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)**3*xKS&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1))&
       /chi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2,3) &
       = (cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*sth&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*(-1.0d0*sph&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*xCKS(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)   +cph(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)*xCKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,2)))/chi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Jacobian fully assembled
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! transform contravariant four-vector
    call matvec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,u4CKS,Jac,u4KS)

    if (present(J)) J = Jac

  end subroutine u4CKStoKS
  !=============================================================================
  subroutine d4CKStoKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xCKS,d4CKS,d4KS)

    ! 04.06.2018 Hector Olivares
    ! Transforms the (covariant) four-vector d4CKS from Cartesian Kerr-Schild
    ! to spherical Kerr-Schild
    !
    ! The Jacobian comes from the transformation u4KStoCKS
    !
    use mod_amrvacdef

    integer, intent(in)                                    :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: xCKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(in)   :: d4CKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(out)  :: d4KS
    ! .. local ..
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3)               :: xKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3)         :: Jac
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3)               :: dummy4CKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3)               :: dummy4KS
    !-----------------------------------------------------------------------------

    ! Dummy vectors to call u4KStoCKS and get the Jacobian

    dummy4CKS = zero
    dummy4KS = zero

    call CKSToKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xCKS,xKS)

    ! Get Jacobian

    call u4KStoCKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xKS,dummy4KS,dummy4CKS,Jac)

    ! transform covariant four-vector
    call matrow(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,d4CKS,Jac,d4KS)

  end subroutine d4CKStoKS
  !=============================================================================
  subroutine CKSToKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xCKS,xKS)

    ! Oliver Porth 18.07.2017
    ! Transforms the coordinates from Cartesian "Kerr-Schild" to spheroidal KS coordinates
    ! Based on Convert_KS_CartKS from Ziri Younsi, BHOSS-code
    use mod_amrvacdef

    integer,intent(in)                                     :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: xCKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(out)  :: xKS
    ! .. local ..
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)                     :: rho
    !-----------------------------------------------------------------------------

    rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)   &
       = xCKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)&
       **2  + xCKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)&
       **2  + xCKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)&
       **2 - eqpar(a_)**2

    ! r
    xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) &
       = SQRT(0.5D0*(rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + &
       SQRT(rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
       **2  + 4.D0*eqpar(a_)**2 * xCKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,3)**2)))   

     ! theta
 !This is only valid for 3D, in 2D, we must simulate CKS in the r-phi or x-y plane:
    where (abs(xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)).lt.smalldouble)
      ! If xKS1 == 0, then we are on the equatorial plane. 
      xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3) = dpi*half
    elsewhere
      xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3) &
         = ACOS( xCKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)&
         /xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1))                                
    end where
   

    
     ! phi
    xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) &
       = ATAN2(xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1) * xCKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) &
         - eqpar(a_)*xCKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1),&
         xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            1) * xCKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            1) + eqpar(a_)*xCKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,2))
   

  end subroutine CKSToKS
  !=============================================================================
  subroutine d4KStoCKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xKS,d4KS,d4CKS)

    ! 04.06.2018 Hector Olivares
    ! Transforms the (covariant) four-vector d4KS from spherical Kerr-Schild
    ! to Cartesian Kerr-Schild
    !
    ! The Jacobian comes from the transformation u4CKStoKS
    !
    use mod_amrvacdef

    integer, intent(in)                                    :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: xKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(in)   :: d4KS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(out)  :: d4CKS
    ! .. local ..
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3)               :: xCKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3)         :: Jac
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3)               :: dummy4KS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3)               :: dummy4CKS
    !-----------------------------------------------------------------------------

    ! Dummy vectors to call u4KStoCKS and get the Jacobian

    dummy4KS = zero
    dummy4CKS = zero

    call KSToCKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xKS,xCKS)

    ! Get Jacobian

    call u4CKStoKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xCKS,dummy4CKS,dummy4KS,Jac)

    ! transform covariant four-vector
    call matrow(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,d4KS,Jac,d4CKS)

  end subroutine d4KStoCKS
  !=============================================================================
  subroutine KSToCKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xKS,xCKS)

    ! Alejandro Cruz 12.07.2017
    ! Transforms the coordinates from spheroidal "Kerr-Schild" to Cartesian coordinates

    ! x = (rKS Cos(phi)  - a*Sin(phi) )Sin(theta)
    ! y = (rKS Sin(phi)  + a*Cos(phi) )Sin(theta)
    ! z = rKS Cos(theta)

    ! See pp 56. Introduction to 3+1 Numerical Relativity,
    ! Alcubierre 2008, Oxford University Press Inc., New York

    ! xKS(ixO^S,1)=r
    ! xKS(ixO^S,2)=theta
    ! xKS(ixO^S,3)=phi
    use mod_amrvacdef

    integer,intent(in)                                     :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: xKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(out)  :: xCKS
    ! .. local ..
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)                     :: sth, cth, sph, cph
    !-----------------------------------------------------------------------------

    
    sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       sin(xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))
    cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       cos(xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))
   
    
    sph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       sin(xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))
    cph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       cos(xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))
   

    xCKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) &
       = (xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)*cph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) - &
       eqpar(a_)*sph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)) * &
       sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) !x
    
    
    xCKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) &
       = (xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)*sph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
         + eqpar(a_)*cph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)) * &
            sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) !y
   
    
    xCKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3) &
       = xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)*cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) !z
   

  end subroutine KSToCKS
  !=============================================================================
  subroutine CylKSToKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xCylKS,xKS)

    ! Hector Olivares 27.05.2019
    ! Transforms the coordinates from Cylindrical "Kerr-Schild" to spheroidal KS coordinates

    use mod_amrvacdef

    integer,intent(in)                                     :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: xCylKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(out)  :: xKS
    ! .. local ..
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)                     :: rho
    !-----------------------------------------------------------------------------

    rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)   &
       = xCylKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)&
       **2  + xCylKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)&
       **2 - eqpar(a_)**2

    ! r
    xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) &
       = SQRT(0.5D0*(rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + &
       SQRT(rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
       **2  + 4.D0*eqpar(a_)**2 * xCylKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,3)**2)))   


     ! theta
    where (abs(xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)).lt.smalldouble)
      ! If xKS1 == 0, then we are on the equatorial plane. 
      xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3) = dpi*half
    elsewhere
      xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3) &
         = ACOS(xCylKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)&
         /xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1))                                
    end where
   

    
     ! phi
    xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) &
       = xCylKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)
   

  end subroutine CylKSToKS
  !=============================================================================
  subroutine KSToCylKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xKS,xCylKS)

    ! Hector Olivares 27.05.2019
    ! Transforms the coordinates from spheroidal "Kerr-Schild" to Cylindrical
    ! Kerr-Schild

    ! R      = sqrt(rKS^2 + a^2) sin(th)
    ! z      = rKS*cos(th)
    ! phiCyl = phiKS

    use mod_amrvacdef

    integer,intent(in)                                     :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: xKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(out)  :: xCylKS
    ! .. local ..
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)                     :: sth, cth
    !-----------------------------------------------------------------------------

    
    sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       sin(xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))
    cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       cos(xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))
   

    ! R
    xCylKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)  &
       = sqrt(xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)&
       **2 - eqpar(a_)**2) * sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    
    
    ! z
    xCylKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3) &
       = xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)*cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
   

    
    ! phi
    xCylKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) &
       = xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)
   

  end subroutine KSToCylKS
  !=============================================================================
  subroutine u4CylKStoBL(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xCylKS,u4CylKS,u4BL,J)

    ! 28.05.2019 Hector Olivares
    ! Transforms the (contravariant) four-vector u4CylKS from Cylindrical Kerr-Schild
    ! to spherical Boyer-Lindquist

    use mod_amrvacdef

    integer, intent(in)                                    :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: xCylKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(in)   :: u4CylKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(out)  :: u4BL
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3), optional, intent(out)  :: J
    ! .. local ..
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3)               :: xKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3)               :: u4dummy1,u4dummy2
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3)         :: Jac1
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3)         :: Jac2
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3)         :: Jac
    !-----------------------------------------------------------------------------

    ! Get Jacobian KS->BL
    call CylKSToKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xCylKS,xKS)
    call u4KSToBL(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xKS,u4dummy1,u4dummy2,J=Jac1)
    
    ! Get Jacobian CylKS->KS
    call u4CylKSToKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xCylKS,u4dummy1,u4dummy2,J&
       =Jac2)

    ! Compose Jacobian CylKS->BL:
    call compose(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,Jac1,Jac2,Jac)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Jacobian fully assembled
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! transform contravariant four-vector
    call matvec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,u4CylKS,Jac,u4BL)

    if (present(J)) J = Jac

  end subroutine u4CylKStoBL
  !=============================================================================
  subroutine u4CylKStoKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xCylKS,u4CylKS,u4KS,J)

    ! 28.05.2019 Hector Olivares
    ! Transforms the (contravariant) four-vector u4CylKS from Cylindrical Kerr-Schild
    ! to spherical Kerr-Schild
    !
    ! Jacobian:
    !
    !( U^0    )      | 1  0                                         0                                      0  | ( U^0     )
    !( U^rKS  ) ____ | 0  r R/(r^2+a^2 cos^2(th))                   (a^2+r^2)cos(th)/(r^2+a^2 cos^2(th))   0  | ( U^R     )
    !( U^th   ) ____ | 0  R cos(th)/(r^2+a^2 cos^2(th)) sin(th)     rsin(th)/(r^2+a^2 cos^2(th))           0  | ( U^z     )
    !( U^phKS )      | 0  0                                         0                                      1  | ( U^phCyl )
    ! 
    !
    use mod_amrvacdef

    integer, intent(in)                                    :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: xCylKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(in)   :: u4CylKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(out)  :: u4KS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3), optional, intent(out)  :: J
    ! .. local ..
    integer                                                :: ix, jx
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3)               :: xKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)                     :: sth, cth, det
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3)         :: Jac
    !-----------------------------------------------------------------------------

    call CylKSToKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xCylKS,xKS)

    
    sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       sin(xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))
    cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       cos(xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))
   
    
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:,:) = 0.0d0

    ! The following expression divides all of the nontrivial components. I add a small number to
    ! avoid divisions by zero. This should happen only close to the singularity.
    det(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)*xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1) + eqpar(a_)*eqpar(a_)*cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) + 1.0e-6

    ! t:
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,0,0) = 1.0d0
    ! r:
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1,1)  &
       = xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)*xCylKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)&
       /det(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1,3) &
       = (eqpar(a_)*eqpar(a_) + xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1)*xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1))*cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
       /det(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    ! theta:
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3,1)  &
       = xCylKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)*cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
       /(det(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*sth&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))

    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3,3) &
       = -xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)*sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
       /det(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    ! phi:
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2,2) = 1.0d0

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Jacobian fully assembled
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! transform contravariant four-vector
    call matvec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,u4CylKS,Jac,u4KS)

    if (present(J)) J = Jac

  end subroutine u4CylKStoKS
  !=============================================================================
  subroutine u4KStoCylKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xKS,u4KS,u4CylKS,J)

    ! Hector Olivares 12.07.2017
    ! Transforms the (contravariant) four-velocity u4KS from Kerr-Schild spherical coordinates
    ! to the cylindrical Kerr-Schild coordinates u4CylKS.  
    !
    !( U^0 )      | 1          0                   0                      0  | ( U^0    )
    !( U^R ) ____ | 0  (rKS/R)sin^2(th)   ((rKS^2+a^2)/R)cos(th)sen(th)   0  | ( U^rKS  )
    !( U^z ) ____ | 0       cos(th)         -rKS sin(theta)               0  | ( U^th   )
    !( U^phCyl )  | 0          0                   0                      1  | ( U^phKS )
    !
    use mod_amrvacdef

    integer, intent(in)                                    :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: xKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(in)   :: u4KS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(out)  :: u4CylKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3), optional, intent(out)  :: J
    ! .. local ..
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3)               :: xCylKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)                     :: sth, cth
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3)         :: Jac
    ! ... hectors test ...
    integer :: k,l
    !-----------------------------------------------------------------------------
    call KSToCylKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xKS,xCylKS)


    
    sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       sin(xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))
    cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       cos(xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))
   

    ! t:
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:,:) = 0.0d0
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,0,0) = 1.0d0

    ! R:
    ! Avoid NaNs and keep zeros where R_cyl is zero, since sen(th) makes everything zero
    ! at the axis.
    where (abs(xCylKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)).gt.0.0d0)
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1,1) &
       = (xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)&
       /xCylKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1))*sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*sth&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1,3) &
       = ((xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)&
       **2+ eqpar(a_)**2)/xCylKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1))*cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
   
    end where
    
    
    ! z:
 !where (cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3).eq.cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3,1)  &
       = cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3,3) &
       = -xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)*sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
!    end where
   
   

    
    ! phi:
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2,2) = 1.0d0
   

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Jacobian fully assembled
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    if (any(Jac(ixO^S,1,1).ne.Jac(ixO^S,1,1))) print *, 'NaN at Jac 11'
!    if (any(Jac(ixO^S,2,1).ne.Jac(ixO^S,2,1))) print *, 'NaN at Jac 21'
!    if (any(Jac(ixO^S,3,1).ne.Jac(ixO^S,3,1))) print *, 'NaN at Jac 31'
!    if (any(Jac(ixO^S,1,2).ne.Jac(ixO^S,1,2))) print *, 'NaN at Jac 12'
!    if (any(Jac(ixO^S,2,2).ne.Jac(ixO^S,2,2))) print *, 'NaN at Jac 22'
!    if (any(Jac(ixO^S,3,2).ne.Jac(ixO^S,3,2))) print *, 'NaN at Jac 32'
!    if (any(Jac(ixO^S,1,3).ne.Jac(ixO^S,1,3))) print *, 'NaN at Jac 13'
!    if (any(Jac(ixO^S,2,3).ne.Jac(ixO^S,2,3))) print *, 'NaN at Jac 23'
!    if (any(Jac(ixO^S,3,3).ne.Jac(ixO^S,3,3))) print *, 'NaN at Jac 33'

    ! transform contravariant four-vector
    call matvec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,u4KS,Jac,u4CylKS)
    
    if (present(J)) J = Jac

  end subroutine u4KStoCylKS
  !=============================================================================
  subroutine d4KStoCylKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xKS,d4KS,d4CylKS)

    ! 28.05.2019 Hector Olivares
    ! Transforms the (covariant) four-vector d4KS from spherical Kerr-Schild
    ! to Cylindrical Kerr-Schild
    !
    ! The Jacobian comes from the transformation u4CylKStoKS
    !
    use mod_amrvacdef

    integer, intent(in)                                    :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: xKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(in)   :: d4KS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(out)  :: d4CylKS
    ! .. local ..
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3)               :: xCylKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3)         :: Jac
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3)               :: dummy4KS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3)               :: dummy4CylKS
    !-----------------------------------------------------------------------------

    ! Dummy vectors to call u4KStoCylKS and get the Jacobian

    dummy4KS = zero
    dummy4CylKS = zero

    call KSToCylKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xKS,xCylKS)

    ! Get Jacobian

    call u4CylKStoKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xCylKS,dummy4CylKS,dummy4KS,&
       Jac)

    ! transform covariant four-vector
    call matrow(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,d4KS,Jac,d4CylKS)

  end subroutine d4KStoCylKS
  !=============================================================================
  subroutine d4CylKStoKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xCylKS,d4CylKS,d4KS)

    ! 28.05.2019 Hector Olivares
    ! Transforms the (covariant) four-vector d4CylKS from Cylindrical Kerr-Schild
    ! to spherical Kerr-Schild
    !
    ! The Jacobian comes from the transformation u4KStoCylKS
    !
    use mod_amrvacdef

    integer, intent(in)                                    :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: xCylKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(in)   :: d4CylKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(out)  :: d4KS
    ! .. local ..
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3)               :: xKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3)         :: Jac
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3)               :: dummy4CylKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3)               :: dummy4KS
    !-----------------------------------------------------------------------------

    ! Dummy vectors to call u4KStoCylKS and get the Jacobian

    dummy4CylKS = zero
    dummy4KS = zero

    call CylKSToKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xCylKS,xKS)

    ! Get Jacobian

    call u4KStoCylKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xKS,dummy4KS,dummy4CylKS,Jac)

    ! transform covariant four-vector
    call matrow(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,d4CylKS,Jac,d4KS)

  end subroutine d4CylKStoKS
  !=============================================================================
  subroutine KSToLRNF(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xKS,ehatu)

    ! Oliver Porth, 14.01.2014
    ! Obtains the tetrad transform to the locally non-rotating frame (LRNF)
    ! Needs KS coordinates on input.
    ! Needs to be applied to tensor in (spherical) KS coordinates.
    ! Returns the transformation matrix e^\hat{\mu}_\nu
    ! 
    ! Follows Rohta Takahashi, MNRAS, 383, 1155-1165, 2008
    ! His equations (11) - (14)
    ! Checked with Mathematica.
    
    use mod_amrvacdef

    integer,intent(in)                                     :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: xKS !KS-coordinate
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3), intent(out) :: ehatu
    ! .. local ..
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)        :: sth, cth, rho2
    !-----------------------------------------------------------------------------

    
    sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       sin(xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))
    cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       cos(xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))
   
    rho2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,r_)&
       **2 + eqpar(a_)**2 * cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)**2

    ehatu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:,:) = zero

    ehatu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,0,0) = 1&
       /Sqrt(1+(2.0d0*eqpar(m_)*xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,r_))/rho2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))

    ehatu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,r_,0) &
       = (2.0d0*eqpar(m_)*xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       r_))/((rho2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)  + &
       2*eqpar(m_)*xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       r_))*Sqrt((-2*eqpar(m_)*xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,r_))/(rho2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)  + 2*eqpar(m_)*xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,r_))+1/(1  - (eqpar(a_)**2*sth(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)**2)/(eqpar(a_)**2  + &
       xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,r_)**2))))

    ehatu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,r_,r_) = 1&
       /Sqrt((-2.0d0*eqpar(m_)*xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,r_))/(rho2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)  + 2*eqpar(m_)*xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,r_))+1/(1  - (1.0d0*eqpar(a_)**2*sth(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)**2)/(eqpar(a_)**2  + &
       xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,r_)**2)))

    
    ehatu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,z_,z_) &
       = sqrt(rho2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
   

    
    ehatu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,phi_,0) &
       = (-2.0d0*eqpar(a_)*eqpar(m_)*sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)**2*xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,r_)*(1  &
         + (2*eqpar(m_)*xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            r_))/rho2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)))&
            /((rho2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)  &
         + 2*eqpar(m_)*xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            r_))*Sqrt(sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
            **2*(eqpar(a_)**2  &
         + (2*eqpar(a_)**2*eqpar(m_)*sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)**2*xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,r_))/rho2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)  &
         + xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,r_)**2)))

    ehatu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,phi_,r_) &
       = (-1.0d0*eqpar(a_)*sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)**2*(1+  &
         (2*eqpar(m_)*xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            r_))/rho2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)))&
            /Sqrt(sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
            **2*(eqpar(a_)**2  &
         + (2*eqpar(a_)**2*eqpar(m_)*sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)**2*xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,r_))/rho2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)  &
         + xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,r_)**2))

    ehatu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,phi_,phi_) &
       = sqrt(sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
       **2*(eqpar(a_)**2  &
         + (2.0d0*eqpar(a_)**2*eqpar(m_)*sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)**2*xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,r_))/rho2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3) &
         + xKS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,r_)**2))
   
    
  end subroutine KSToLRNF
  !=============================================================================
  subroutine BLToComoving(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xBL,u4BL,ehatd)

    ! Oliver Porth, 06.05.2019
    ! Obtains the tetrad transform to a local flat fluid frame
    ! Needs BL coordinates on input.
    ! Needs the (contravariant) fluid four-velocity in BL coordinates (u4BL) as input
    ! Needs to be applied to tensor in (covariant) BL coordinates d4BL
    ! Returns the transformation matrix e^\mu_(\nu)
    !
    ! Follows Beckwith et al., ApJ 678, 1180-1199, 2008
    ! See their Appendix A.

    use mod_bl
    use mod_amrvacdef

    integer,intent(in)                                     :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim), intent(in)   :: xBL !BL-coordinate
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:ndir), intent(in)   :: u4BL !Fluid-vel. in BL
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:ndir,0:ndir), intent(out) :: ehatd
    ! .. local ..
    integer                  :: ix1,ix2,ix3
    double precision         :: g(0:ndir,0:ndir), k_r, k_phi, k_theta,&
        u(0:ndir), tmp
    !-----------------------------------------------------------------------------

    ehatd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:,:) = zero

    do ix3=ixOmin3,ixOmax3
    do ix2=ixOmin2,ixOmax2
    do ix1=ixOmin1,ixOmax1

        call get_g4_BL(xBL(ix1,ix2,ix3,1),xBL(ix1,ix2,ix3,2),xBL(ix1,ix2,ix3,&
           3),g)
        u(0:ndir)  = u4BL(ix1,ix2,ix3,0:ndir)

        tmp = u(0)*(2.0d0*g(0,phi_)*u(phi_) + g(0,0)*u(0))
        k_phi   = sqrt( abs( g(phi_,phi_)*u(phi_)**2 + tmp ) )
        k_r     = sqrt( abs( g(phi_,phi_)*u(phi_)**2 + g(r_,r_)*u(r_)&
           **2 + tmp ) )
        k_theta = sqrt( abs( g(phi_,phi_)*u(phi_)**2 + g(r_,r_)*u(r_)&
           **2 + g(z_,z_)*u(z_)**2 + tmp ) )
        
        ehatd(ix1,ix2,ix3,0,0:ndir) = u(0:ndir)

        ehatd(ix1,ix2,ix3,r_,0)     = -1.0d0/(k_theta*k_phi) * sqrt(g(r_,&
           r_))*u(r_)*u(0)
        ehatd(ix1,ix2,ix3,r_,r_)    = -1.0d0/(k_theta*k_phi) * k_phi**2&
           /sqrt(g(r_,r_))
        ehatd(ix1,ix2,ix3,r_,phi_)  = -1.0d0/(k_theta*k_phi) * sqrt(g(r_,&
           r_))*u(r_)*u(phi_)
        
        ehatd(ix1,ix2,ix3,z_,0)    = -1.0d0/(k_theta*k_r) * sqrt(g(z_,&
           z_))*u(z_)*u(0) 
        ehatd(ix1,ix2,ix3,z_,r_)   = -1.0d0/(k_theta*k_r) * sqrt(g(z_,&
           z_))*u(z_)*u(r_) 
        ehatd(ix1,ix2,ix3,z_,z_)   = -1.0d0/(k_theta*k_r) * k_r**2&
           /sqrt(g(z_,z_))
        ehatd(ix1,ix2,ix3,z_,phi_) = -1.0d0/(k_theta*k_r) * sqrt(g(z_,&
           z_))*u(z_)*u(phi_)

        tmp = -1.0d0/( k_phi*sqrt( abs( -g(0,phi_)**2 + g(phi_,phi_)*g(0,&
           0) ) ) )
        ehatd(ix1,ix2,ix3,phi_,0)    = tmp * ( g(phi_,phi_)*u(phi_) + g(0,&
           phi_)*u(0) )
        ehatd(ix1,ix2,ix3,phi_,phi_) = - tmp * ( g(0,phi_)*u(phi_) + g(0,&
           0)*u(0) )
        
    end do
    end do
    end do
    
  end subroutine BLToComoving
  !=============================================================================
  subroutine boost(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,u4,uhatu,ubaru)

    ! Lorentz-boosts a four-vector in an orthogonal (tetrad-) basis.
    ! Written for contravariant vectors u^{\hat{\beta}}
    ! The transformation for covariant vectors follows simply by setting u4->-u4
    ! u4 is the four-velocity of the frame to boost into, also in the orthogonal basis.

    use mod_amrvacdef
    
    integer,intent(in)                                           :: ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,&
       ixOmax1,ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(in)         :: u4
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(in)         :: uhatu
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(out)        :: ubaru
    ! .. local ..
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)                           :: lfac
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3)               :: Lambda
    integer                                                      :: ix, jx
    !-----------------------------------------------------------------------------

    lfac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = sqrt(1.0d0 +  &
       u4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)&
       **2 + u4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)&
       **2 + u4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)**2 )

    Lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,0,0) &
       = lfac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
     Lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,0,1) &
        = -u4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)
      Lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,0,2) &
         = -u4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)
      Lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,0,3) &
         = -u4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)
     Lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1,0) &
        = -u4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)
      Lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2,0) &
         = -u4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)
      Lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3,0) &
         = -u4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)

    do ix = 1,ndir
       do jx = 1,ndir
          Lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ix,jx) &
             = dble(kr(ix,jx)) + one/(one+lfac(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3)) * u4(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,ix) * u4(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,jx)
       end do
    end do

    call matvec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,uhatu,Lambda,ubaru)
    
  end subroutine boost
  !=============================================================================
  subroutine compose(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,Jac1,Jac2,Jac3)

    ! Does a four-matrix multiplication like so:
    ! Jac3(i,j) = Jac1(i,k)* Jac2(k,j)
    use mod_amrvacdef

    integer,intent(in)                                           :: ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,&
       ixOmax1,ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3), intent(in)   :: Jac1, Jac2
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3), intent(out)  :: Jac3
    ! .. local ..
    integer                                                      :: ix, jx, kx
    !-----------------------------------------------------------------------------

    Jac3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:,:) = 0.0d0
    
    do ix=0,ndir
       do jx=0,ndir
          do kx=0,ndir
             Jac3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ix,jx) &
                = Jac3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ix,&
                jx) + Jac1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ix,&
                kx) * Jac2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,kx,&
                jx)
          end do
       end do
    end do

  end subroutine compose
  !=============================================================================
  subroutine transpose(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,Jac1,Jac2)

    ! Transposes a four-matrix like so:
    ! Jac2(i,j) = Jac1(j,i)
    use mod_amrvacdef

    integer,intent(in)                                           :: ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,&
       ixOmax1,ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3), intent(in)   :: Jac1
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3), intent(out)  :: Jac2
    ! .. local ..
    integer                                                      :: ix, jx
    !-----------------------------------------------------------------------------
    
    do ix=0,ndir
       do jx=0,ndir
             Jac2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ix,jx) &
                = Jac1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,jx,ix)
       end do
    end do

  end subroutine transpose
  !=============================================================================
  subroutine matvec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,u4in,Jac,u4out)

    ! Multiplies a four-vector and a matrix like so:
    ! u4out(i) = Jac(i,j) * u4in(j)
    use mod_amrvacdef

    integer,intent(in)                                           :: ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,&
       ixOmax1,ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3), intent(in)   :: Jac
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(in)         :: u4in
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(out)        :: u4out
    ! .. local ..
    integer                                                      :: ix, jx
    !-----------------------------------------------------------------------------

    u4out(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:) = 0.0d0
    do ix=0,ndir
       do jx=0,ndir
             u4out(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ix) &
                = u4out(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                ix) + Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ix,&
                jx) * u4in(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,jx)
       end do
    end do

  end subroutine matvec
  !=============================================================================
  subroutine matrow(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,d4in,Jac,d4out)

    ! Multiplies a row four-vector (1-form) and a matrix like so:
    ! d4out(i) = Jac(j,i) * d4in(j)
    use mod_amrvacdef

    integer,intent(in)                                           :: ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,&
       ixOmax1,ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3), intent(in)   :: Jac
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(in)         :: d4in
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(out)        :: d4out
    ! .. local ..
    integer                                                      :: ix, jx
    !-----------------------------------------------------------------------------

    d4out(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:) = 0.0d0
    do ix=0,ndir
       do jx=0,ndir
             d4out(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ix) &
                = d4out(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                ix) + Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,jx,&
                ix) * d4in(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,jx)
       end do
    end do

  end subroutine matrow
  !=============================================================================
end module mod_transform
!=============================================================================
