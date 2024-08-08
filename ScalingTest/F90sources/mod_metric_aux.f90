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
! Sets up more flexible structures
! Inspired by mpiamrvac.org setup procedures.
! 18/08/2017 Oliver Porth
!=============================================================================
module mod_metric_aux
  implicit none
  public

  procedure(sub_get_gammainv_component_analytic), pointer  :: &
     get_gammainv_component_analytic => null()
  procedure(sub_get_gammainv_component_analytic), pointer  :: &
     get_gammainv_component => null()
  procedure(sub_SPToCoord), pointer                        :: SPToCoord &
     => null()
  procedure(dummy_LRNFu), pointer                          :: LRNFu => null()
  procedure(dummy_u4CoordToKS), pointer                    :: u4CoordToKS &
     => null()
  procedure(dummy_u4CoordToBL), pointer                    :: u4CoordToBL &
     => null()
  procedure(dummy_d4CoordToKS), pointer                    :: d4CoordToKS &
     => null()

  abstract interface
     subroutine sub_get_gammainv_component_analytic(iin,jin,x1,x2,x3,ginv,&
        iszero,dginvdk_iszero,dginvdk,kdir)
       integer, intent(in)                      :: iin,jin
       integer, optional, intent(in)            :: kdir
       double precision, intent(in)             :: x1,x2,x3
       double precision, intent(out)            :: ginv
       logical, optional, intent(out)           :: iszero, dginvdk_iszero
       double precision, optional, intent(out)  :: dginvdk
     end subroutine sub_get_gammainv_component_analytic

     subroutine sub_SPToCoord(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xKS,xCKS)
       integer,intent(in)                                     :: ixImin1,&
          ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,&
          ixOmax1,ixOmax2,ixOmax3
       double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
          ixImin3:ixImax3,1:3), intent(in)   :: xKS
       double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
          ixImin3:ixImax3,1:3), intent(out)  :: xCKS
     end subroutine sub_SPToCoord
     
  end interface
contains

  !=============================================================================
  subroutine dummy_get_gammainv_component_analytic(iin,jin,x1,x2,x3,ginv,&
     iszero,dginvdk_iszero,dginvdk,kdir)
    use mod_amrvacdef
    integer, intent(in)                      :: iin,jin
    integer, optional, intent(in)            :: kdir
    double precision, intent(in)             :: x1,x2,x3
    double precision, intent(out)            :: ginv
    logical, optional, intent(out)           :: iszero, dginvdk_iszero
    double precision, optional, intent(out)  :: dginvdk
    !-----------------------------------------------------------------------------

    ginv = zero ! avoiding compiler warning.
    
    call mpistop&
('get_gammainv_component_analytic: Not implemented for these coordinates!')
    
  end subroutine dummy_get_gammainv_component_analytic
  !=============================================================================
  subroutine dummy_SPToCoord(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xKS,xCKS)
    use mod_amrvacdef
    integer,intent(in)                                     :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: xKS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(out)  :: xCKS
    !-----------------------------------------------------------------------------
    
    xCKS = zero ! avoiding compiler warning.

    call mpistop('SPToCoord: Not implemented for these coordinates!')

  end subroutine dummy_SPToCoord
  !=============================================================================
  subroutine dummy_LRNFu(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,m,ehatu)

    ! Obtains orthonormal tetrad vectors in matrix form.
    ! This takes in already the metric datastructure 'm', I might regret
    ! this choice later on. 
    
    use mod_amrvacdef
    integer, intent(in)                                         :: ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,&
       ixOmax1,ixOmax2,ixOmax3
    type(metric)                                                :: m
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3), intent(out) :: ehatu
    ! .. local ..
    !-----------------------------------------------------------------------------

    ehatu = zero ! avoiding compiler warning

    call mpistop('LRNFu: Not implemented for these coordinates!')
    
  end subroutine dummy_LRNFu
  !=============================================================================
  subroutine dummy_u4CoordToKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xCoord,u4Coord,&
     u4KS,J)

    ! Transforms the current coordinates four vector u4Coord to
    ! the ordinary KS four vector u4KS.
    !
    use mod_amrvacdef

    integer,intent(in)                                     :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: xCoord
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(in)   :: u4Coord
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(out)  :: u4KS
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3), optional, intent(out)  :: J
    !-----------------------------------------------------------------------------
    
    u4KS = zero ! avoiding compiler warning

    call mpistop('u4CoordToKS: Not implemented for these coordinates!')
    
  end subroutine dummy_u4CoordToKS
  !=============================================================================
  subroutine dummy_u4CoordToBL(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xCoord,u4Coord,&
     u4BL,J)

    ! Transforms the current coordinates four vector u4Coord to
    ! the ordinary BL four vector u4BL.
    !
    use mod_amrvacdef

    integer,intent(in)                                     :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: xCoord
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(in)   :: u4Coord
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(out)  :: u4BL
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3), optional, intent(out)  :: J
    !-----------------------------------------------------------------------------
    
    u4BL = zero ! avoiding compiler warning

    call mpistop('u4CoordToBL: Not implemented for these coordinates!')
    
  end subroutine dummy_u4CoordToBL
  !=============================================================================
  subroutine dummy_d4CoordToKS(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xCoord,d4Coord,&
     d4KS)

    ! Transforms the current coordinates four vector u4Coord to
    ! the ordinary KS four vector u4KS.
    !
    use mod_amrvacdef

    integer,intent(in)                                     :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: xCoord
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(in)   :: d4Coord
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(out)  :: d4KS
    !-----------------------------------------------------------------------------
    
    d4KS = zero ! avoiding compiler warning

    call mpistop('d4CoordToKS: Not implemented for these coordinates!')
    
  end subroutine dummy_d4CoordToKS
  !=============================================================================
end module mod_metric_aux
!=============================================================================
