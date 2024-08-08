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
module mod_limiter

  implicit none
  public

  !> radius of the asymptotic region [0.001, 10], larger means more accurate in smooth
  !> region but more overshooting at discontinuities
  double precision :: cada3_radius
  integer, parameter :: limiter_minmod = 1
  integer, parameter :: limiter_woodward = 2
  integer, parameter :: limiter_mcbeta = 3
  integer, parameter :: limiter_superbee = 4
  integer, parameter :: limiter_vanleer = 5
  integer, parameter :: limiter_albada = 6
  integer, parameter :: limiter_koren = 7
  integer, parameter :: limiter_cada = 8
  integer, parameter :: limiter_cada3 = 9
  ! Special cases
  integer, parameter :: limiter_ppm = 10
  integer, parameter :: limiter_mp5 = 11
  integer, parameter :: limiter_weno5 = 12
  integer, parameter :: limiter_wenoZP = 13


contains

  !=============================================================================
  integer function limiter_type(namelim)
    character(len=*), intent(in) :: namelim

    select case (namelim)
    case ('minmod')
       limiter_type = limiter_minmod
    case ('woodward')
       limiter_type = limiter_woodward
    case ('mcbeta')
       limiter_type = limiter_mcbeta
    case ('superbee')
       limiter_type = limiter_superbee
    case ('vanleer')
       limiter_type = limiter_vanleer
    case ('albada')
       limiter_type = limiter_albada
    case ('koren')
       limiter_type = limiter_koren
    case ('cada')
       limiter_type = limiter_cada
    case ('cada3')
       limiter_type = limiter_cada3
    case ('ppm')
       limiter_type = limiter_ppm
    case ('mp5')
       limiter_type = limiter_mp5
    case ('weno5')
       limiter_type = limiter_weno5
    case ('wenoZP')
       limiter_type = limiter_wenoZP
    case default
       limiter_type = -1
       write(*,*) 'Unknown limiter: ', namelim
       call mpistop("No such limiter")
    end select
  end function limiter_type
  !=============================================================================
  !> Limit the centered dwC differences within ixC for iw in direction idim.
  !> The limiter is chosen according to typelimiter.
  !>
  !> Note that this subroutine is called from upwindLR (hence from methods
  !> like tvdlf, hancock, hll(c) etc) or directly from tvd.t,
  !> but also from the gradientS and divvectorS subroutines in geometry.t
  !> Accordingly, the typelimiter here corresponds to one of limiter
  !> or one of gradient_limiter.
  subroutine dwlimiter2(dwC,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,typelim,ldw,rdw)

    use mod_amrvacdef

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, idims
    double precision, intent(in) :: dwC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    integer, intent(in) :: typelim
    !> Result using left-limiter (same as right for symmetric)
    double precision, intent(out), optional :: ldw(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)
    !> Result using right-limiter (same as left for symmetric)
    double precision, intent(out), optional :: rdw(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)

    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    integer :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, hxOmin1,&
       hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3
    double precision, parameter :: qsmall=1.d-12, qsmall2=2.d-12
    double precision, parameter :: eps = sqrt(epsilon(1.0d0))

    ! mcbeta limiter parameter value
    double precision, parameter :: c_mcbeta=1.4d0
    ! cada limiter parameter values
    double precision, parameter :: cadalfa=0.5d0, cadbeta=2.0d0, cadgamma&
       =1.6d0
    ! full third order cada limiter
    double precision :: rdelinv
    double precision :: ldwA(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       ldwB(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmpeta(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, parameter :: cadepsilon=1.d-14, invcadepsilon&
       =1.d14,cada3_radius=0.1d0
    !-----------------------------------------------------------------------------

    ! Contract indices in idim for output.
    ixOmin1=ixCmin1+kr(idims,1);ixOmin2=ixCmin2+kr(idims,2)
    ixOmin3=ixCmin3+kr(idims,3); ixOmax1=ixCmax1;ixOmax2=ixCmax2
    ixOmax3=ixCmax3;
    hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
    hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
    hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);

    ! About the notation: the conventional argument theta (the ratio of slopes)
    ! would be given by dwC(ixO^S)/dwC(hxO^S). However, in the end one
    ! multiplies phi(theta) by dwC(hxO^S), which is incorporated in the
    ! equations below. The minmod limiter can for example be written as:
    ! A:
    ! max(0.0d0, min(1.0d0, dwC(ixO^S)/dwC(hxO^S))) * dwC(hxO^S)
    ! B:
    ! tmp(ixO^S)*max(0.0d0,min(abs(dwC(ixO^S)),tmp(ixO^S)*dwC(hxO^S)))
    ! where tmp(ixO^S)=sign(1.0d0,dwC(ixO^S))

    select case (typelim)
    case (limiter_minmod)
       ! Minmod limiter eq(3.51e) and (eq.3.38e) with omega=1
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
          =sign(one,dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
          =tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)* max(zero,&
          min(abs(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)),&
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*dwC&
          (hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3)))
       if (present(ldw)) ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       if (present(rdw)) rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
    case (limiter_woodward)
       ! Woodward and Collela limiter (eq.3.51h), a factor of 2 is pulled out
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
          =sign(one,dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
          =2*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)* max(zero,&
          min(abs(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)),&
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*dwC&
          (hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3),&
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*quarter*(dwC&
          (hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3)+dwC&
          (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))))
       if (present(ldw)) ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       if (present(rdw)) rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
    case (limiter_mcbeta)
       ! Woodward and Collela limiter, with factor beta
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
          =sign(one,dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
          =tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)* max(zero,&
          min(c_mcbeta*abs(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)),c_mcbeta*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
          hxOmin3:hxOmax3),tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*half*(dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
          hxOmin3:hxOmax3)+dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3))))
       if (present(ldw)) ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       if (present(rdw)) rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
    case (limiter_superbee)
       ! Roes superbee limiter (eq.3.51i)
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
          =sign(one,dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
          =tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)* max(zero,&
          min(2*abs(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)),&
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*dwC&
          (hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3)),&
          min(abs(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)),&
          2*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*dwC&
          (hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3)))
       if (present(ldw)) ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       if (present(rdw)) rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
    case (limiter_vanleer)
       ! van Leer limiter (eq 3.51f), but a missing delta2=1.D-12 is added
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
          =2*max(dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3)*dwC&
          (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),zero) &
          /(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)+dwC&
          (hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3)+qsmall)
       if (present(ldw)) ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       if (present(rdw)) rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
    case (limiter_albada)
       ! Albada limiter (eq.3.51g) with delta2=1D.-12
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
          =(dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3)*(dwC&
          (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
          **2+qsmall)+dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*&
          (dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3)**2+qsmall))&
          /(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
          **2+dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3)**2+qsmall2)
       if (present(ldw)) ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       if (present(rdw)) rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
    case (limiter_koren)
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
          =sign(one,dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
       tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
          =min(2*abs(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)),&
          2*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*dwC&
          (hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3))
       if (present(ldw)) then
          ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
             =tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)* max(zero,&
             min(tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
             (dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3)*tmp&
             (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)+2*abs(dwC&
             (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)))*third))
       end if
       if (present(rdw)) then
          rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
             =tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)* max(zero,&
             min(tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
             (2*dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3)*tmp&
             (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)+abs(dwC&
             (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)))*third))
       end if
    case (limiter_cada)
       ! This limiter has been rewritten in the usual form, and uses a division
       ! of the gradients.
       if (present(ldw)) then
          ! Cada Left variant
          ! Compute theta, but avoid division by zero
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
             =dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3)&
             /(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + &
             sign(eps, dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)))
          tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
             =(2+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))*third
          ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
             = max(zero,min(tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3), max(-cadalfa*tmp(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3), min(cadbeta*tmp&
             (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
              tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
              cadgamma)))) * dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)
       end if

       if (present(rdw)) then
          ! Cada Right variant
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
             =dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
             /(dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3) + &
             sign(eps, dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3)))
          tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
             =(2+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))*third
          rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
             = max(zero,min(tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3), max(-cadalfa*tmp(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3), min(cadbeta*tmp&
             (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
              tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
              cadgamma)))) * dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
             hxOmin3:hxOmax3)
       end if
    case (limiter_cada3)
       rdelinv=one/(cada3_radius*dxlevel(idims))**2
       tmpeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
          =(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
          **2+dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3)**2)*rdelinv

       if (present(ldw)) then
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
             =dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3)&
             /(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + &
             sign(eps, dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)))
          ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
             =(two+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))*third
          ldwB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
             = max(zero,min(ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3), max(-cadalfa*tmp(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3), min(cadbeta*tmp&
             (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
              ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
              cadgamma))))
          where(tmpeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
             <=one-cadepsilon)
             ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                =ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
          elsewhere(tmpeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
             >=one+cadepsilon)
             ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                =ldwB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
          elsewhere
             tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                =(tmpeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)-one)&
                *invcadepsilon
             ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                =half*( (one-tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3))*ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3) +(one+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3))*ldwB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3))
          endwhere
          ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
             =ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) * &
             dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
       end if

       if (present(rdw)) then
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
             =dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
             /(dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3) + &
             sign(eps, dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3)))
          ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
             =(two+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))*third
          ldwB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
             = max(zero,min(ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3), max(-cadalfa*tmp(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3), min(cadbeta*tmp&
             (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
              ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
              cadgamma))))
          where(tmpeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
             <=one-cadepsilon)
             rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                =ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
          elsewhere(tmpeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
             >=one+cadepsilon)
             rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                =ldwB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
          elsewhere
             tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                =(tmpeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)-one)&
                *invcadepsilon
             rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                =half*( (one-tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3))*ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3) +(one+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3))*ldwB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3))
          endwhere
          rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
             =rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) * &
             dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3)
       end if

    case default
       call mpistop("Error in dwLimiter: unknown limiter")
    end select

  end subroutine dwlimiter2
  !=============================================================================
  
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

!============================================================================
subroutine PPMlimitervar(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idims,q,qCT,qLC,qRC)

! references:
! Mignone et al 2005, ApJS 160, 199,
! Miller and Colella 2002, JCP 183, 26
! Fryxell et al. 2000 ApJ, 131, 273 (Flash)
! baciotti Phd (http://www.aei.mpg.de/~baiotti/Baiotti_PhD.pdf)
! version : april 2009
! author: zakaria.meliani@wis.kuleuven.be

use mod_amrvacdef

integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3, idims
double precision, intent(in)    :: q(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3),qCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

double precision, intent(inout) :: qRC(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3),qLC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

double precision,dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)  :: dqC,d2qC,ldq
double precision,dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)  :: qMin,qMax,tmp

integer   :: lxCmin1,lxCmin2,lxCmin3,lxCmax1,lxCmax2,lxCmax3,lxRmin1,lxRmin2,&
   lxRmin3,lxRmax1,lxRmax2,lxRmax3
integer   :: ixLLmin1,ixLLmin2,ixLLmin3,ixLLmax1,ixLLmax2,ixLLmax3,ixLmin1,&
   ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
   ixOmax2,ixOmax3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,ixRRmin1,&
   ixRRmin2,ixRRmin3,ixRRmax1,ixRRmax2,ixRRmax3
integer   :: hxLmin1,hxLmin2,hxLmin3,hxLmax1,hxLmax2,hxLmax3,hxCmin1,hxCmin2,&
   hxCmin3,hxCmax1,hxCmax2,hxCmax3,hxRmin1,hxRmin2,hxRmin3,hxRmax1,hxRmax2,&
   hxRmax3
integer   :: kxLLmin1,kxLLmin2,kxLLmin3,kxLLmax1,kxLLmax2,kxLLmax3,kxLmin1,&
   kxLmin2,kxLmin3,kxLmax1,kxLmax2,kxLmax3,kxCmin1,kxCmin2,kxCmin3,kxCmax1,&
   kxCmax2,kxCmax3,kxRmin1,kxRmin2,kxRmin3,kxRmax1,kxRmax2,kxRmax3,kxRRmin1,&
   kxRRmin2,kxRRmin3,kxRRmax1,kxRRmax2,kxRRmax3
!--------------------------------------------------------------------------
ixOmin1=ixmin1-kr(idims,1);ixOmin2=ixmin2-kr(idims,2)
ixOmin3=ixmin3-kr(idims,3);ixOmax1=ixmax1+kr(idims,1)
ixOmax2=ixmax2+kr(idims,2);ixOmax3=ixmax3+kr(idims,3); !ixO[ixMmin1-1,ixMmax1+1]
ixLmin1=ixOmin1-kr(idims,1);ixLmin2=ixOmin2-kr(idims,2)
ixLmin3=ixOmin3-kr(idims,3);ixLmax1=ixOmax1-kr(idims,1)
ixLmax2=ixOmax2-kr(idims,2);ixLmax3=ixOmax3-kr(idims,3); !ixL[ixMmin1-2,ixMmax1]
ixLLmin1=ixLmin1-kr(idims,1);ixLLmin2=ixLmin2-kr(idims,2)
ixLLmin3=ixLmin3-kr(idims,3);ixLLmax1=ixLmax1-kr(idims,1)
ixLLmax2=ixLmax2-kr(idims,2);ixLLmax3=ixLmax3-kr(idims,3); !ixLL[ixMmin1-3,ixMmax1-1]
ixRmin1=ixOmin1+kr(idims,1);ixRmin2=ixOmin2+kr(idims,2)
ixRmin3=ixOmin3+kr(idims,3);ixRmax1=ixOmax1+kr(idims,1)
ixRmax2=ixOmax2+kr(idims,2);ixRmax3=ixOmax3+kr(idims,3); !ixR=[iMmin1,ixMmax+2]
ixRRmin1=ixRmin1+kr(idims,1);ixRRmin2=ixRmin2+kr(idims,2)
ixRRmin3=ixRmin3+kr(idims,3);ixRRmax1=ixRmax1+kr(idims,1)
ixRRmax2=ixRmax2+kr(idims,2);ixRRmax3=ixRmax3+kr(idims,3); !ixRR=[iMmin1+1,ixMmax+3]

hxCmin1=ixOmin1;hxCmin2=ixOmin2;hxCmin3=ixOmin3;hxCmax1=ixmax1;hxCmax2=ixmax2
hxCmax3=ixmax3; !hxC = [ixMmin-1,ixMmax]
hxLmin1=hxCmin1-kr(idims,1);hxLmin2=hxCmin2-kr(idims,2)
hxLmin3=hxCmin3-kr(idims,3);hxLmax1=hxCmax1-kr(idims,1)
hxLmax2=hxCmax2-kr(idims,2);hxLmax3=hxCmax3-kr(idims,3); !hxL = [ixMmin-2,ixMmax-1]
hxRmin1=hxCmin1+kr(idims,1);hxRmin2=hxCmin2+kr(idims,2)
hxRmin3=hxCmin3+kr(idims,3);hxRmax1=hxCmax1+kr(idims,1)
hxRmax2=hxCmax2+kr(idims,2);hxRmax3=hxCmax3+kr(idims,3); !hxR = [ixMmin,ixMmax+1]

kxCmin1=ixLLmin1;kxCmin2=ixLLmin2;kxCmin3=ixLLmin3; kxCmax1=ixRmax1
kxCmax2=ixRmax2;kxCmax3=ixRmax3; !kxC=[iMmin1-3,ixMmax1+2]
kxLmin1=kxCmin1-kr(idims,1);kxLmin2=kxCmin2-kr(idims,2)
kxLmin3=kxCmin3-kr(idims,3);kxLmax1=kxCmax1-kr(idims,1)
kxLmax2=kxCmax2-kr(idims,2);kxLmax3=kxCmax3-kr(idims,3); !kxL=[iMmin1-4,ixMmax1+1]
kxRmin1=kxCmin1+kr(idims,1);kxRmin2=kxCmin2+kr(idims,2)
kxRmin3=kxCmin3+kr(idims,3);kxRmax1=kxCmax1+kr(idims,1)
kxRmax2=kxCmax2+kr(idims,2);kxRmax3=kxCmax3+kr(idims,3); !kxR=[iMmin1-2,ixMmax1+3]

lxCmin1=ixLLmin1-kr(idims,1);lxCmin2=ixLLmin2-kr(idims,2)
lxCmin3=ixLLmin3-kr(idims,3);lxCmax1=ixRRmax1;lxCmax2=ixRRmax2
lxCmax3=ixRRmax3; !ixC=[iMmin1-4,ixMmax1+3]
lxRmin1=lxCmin1+kr(idims,1);lxRmin2=lxCmin2+kr(idims,2)
lxRmin3=lxCmin3+kr(idims,3);lxRmax1=lxCmax1+kr(idims,1)
lxRmax2=lxCmax2+kr(idims,2);lxRmax3=lxCmax3+kr(idims,3); !lxR=[iMmin1-3,ixMmax1+4]


dqC(lxCmin1:lxCmax1,lxCmin2:lxCmax2,lxCmin3:lxCmax3)=q(lxRmin1:lxRmax1,&
   lxRmin2:lxRmax2,lxRmin3:lxRmax3)-q(lxCmin1:lxCmax1,lxCmin2:lxCmax2,&
   lxCmin3:lxCmax3)
! Eq. 64,  Miller and Colella 2002, JCP 183, 26
d2qC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3)=half*(q(kxRmin1:kxRmax1,&
   kxRmin2:kxRmax2,kxRmin3:kxRmax3)-q(kxLmin1:kxLmax1,kxLmin2:kxLmax2,&
   kxLmin3:kxLmax3))
where(dqC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3)*dqC&
   (kxLmin1:kxLmax1,kxLmin2:kxLmax2,kxLmin3:kxLmax3)>zero)
   ! Store the sign of d2qC in qMin
   qMin(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3)= sign(one,&
      d2qC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3))
   ! Eq. 65,  Miller and Colella 2002, JCP 183, 26
   ldq(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3)= qMin(kxCmin1:kxCmax1,&
      kxCmin2:kxCmax2,kxCmin3:kxCmax3)*min(dabs(d2qC(kxCmin1:kxCmax1,&
      kxCmin2:kxCmax2,kxCmin3:kxCmax3)),2.0d0*dabs(dqC(kxLmin1:kxLmax1,&
      kxLmin2:kxLmax2,kxLmin3:kxLmax3)),2.0d0*dabs(dqC(kxCmin1:kxCmax1,&
      kxCmin2:kxCmax2,kxCmin3:kxCmax3)))
elsewhere
   ldq(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3)=zero
end where

! Eq. 66, Miller and Colella 2002, JCP 183, 26
qLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=qLC(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,ixOmin3:ixOmax3)+half*dqC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   ixOmin3:ixOmax3)+(ldq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)-ldq&
   (ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3))/6.0d0

qRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3)=qRC(ixLmin1:ixLmax1,&
   ixLmin2:ixLmax2,ixLmin3:ixLmax3) -(half*dqC(ixLmin1:ixLmax1,&
   ixLmin2:ixLmax2,ixLmin3:ixLmax3)-(ldq(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
   ixLmin3:ixLmax3)-ldq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))&
   /6.0d0)

! Eq. B9, page 217, Mignone et al 2005, ApJS
where((qRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3)-qCT&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))*(qCT(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,ixOmin3:ixOmax3)-qLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   ixOmin3:ixOmax3))<=zero)
   qRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3)=qCT(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3)
   qLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=qCT(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3)
end where

qMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=(qLC(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,ixOmin3:ixOmax3)-qRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
   ixLmin3:ixLmax3))*(qCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)-&
   (qLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)+qRC(ixLmin1:ixLmax1,&
   ixLmin2:ixLmax2,ixLmin3:ixLmax3))/2.0d0)
qMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=(qLC(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,ixOmin3:ixOmax3)-qRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
   ixLmin3:ixLmax3))**2.0d0/6.0d0
tmp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3)=qRC(ixLmin1:ixLmax1,&
   ixLmin2:ixLmax2,ixLmin3:ixLmax3)

! Eq. B10, page 218, Mignone et al 2005, ApJS
where(qMax(hxRmin1:hxRmax1,hxRmin2:hxRmax2,hxRmin3:hxRmax3)&
   >qMin(hxRmin1:hxRmax1,hxRmin2:hxRmax2,hxRmin3:hxRmax3))
   qRC(hxCmin1:hxCmax1,hxCmin2:hxCmax2,hxCmin3:hxCmax3)= 3.0d0*qCT&
      (hxRmin1:hxRmax1,hxRmin2:hxRmax2,hxRmin3:hxRmax3)-2.0d0*qLC&
      (hxRmin1:hxRmax1,hxRmin2:hxRmax2,hxRmin3:hxRmax3)
end where
! Eq. B11, page 218, Mignone et al 2005, ApJS
where(qMax(hxCmin1:hxCmax1,hxCmin2:hxCmax2,hxCmin3:hxCmax3)&
   <-qMin(hxCmin1:hxCmax1,hxCmin2:hxCmax2,hxCmin3:hxCmax3))
   qLC(hxCmin1:hxCmax1,hxCmin2:hxCmax2,hxCmin3:hxCmax3)= 3.0d0*qCT&
      (hxCmin1:hxCmax1,hxCmin2:hxCmax2,hxCmin3:hxCmax3)-2.0d0*tmp&
      (hxLmin1:hxLmax1,hxLmin2:hxLmax2,hxLmin3:hxLmax3)
end where

end subroutine PPMlimitervar
!============================================================================
subroutine PPMlimiter(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixmin1,&
   ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idims,w,wCT,wLC,wRC)

! references:
! Mignone et al 2005, ApJS 160, 199, 
! Miller and Colella 2002, JCP 183, 26 
! Fryxell et al. 2000 ApJ, 131, 273 (Flash)
! baciotti Phd (http://www.aei.mpg.de/~baiotti/Baiotti_PhD.pdf)
! version : april 2009
! author: zakaria.meliani@wis.kuleuven.be

use mod_amrvacdef

integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3, idims
double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw),wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nw)
double precision, intent(inout) :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw),wLC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nw) 
! .. local ..
double precision,dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)  :: dwC,d2wC,ldw
double precision,dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)  :: wMin,wMax,tmp
integer   :: lxCmin1,lxCmin2,lxCmin3,lxCmax1,lxCmax2,lxCmax3,lxRmin1,lxRmin2,&
   lxRmin3,lxRmax1,lxRmax2,lxRmax3
integer   :: ixLLmin1,ixLLmin2,ixLLmin3,ixLLmax1,ixLLmax2,ixLLmax3,ixLmin1,&
   ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
   ixOmax2,ixOmax3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,ixRRmin1,&
   ixRRmin2,ixRRmin3,ixRRmax1,ixRRmax2,ixRRmax3
integer   :: hxLmin1,hxLmin2,hxLmin3,hxLmax1,hxLmax2,hxLmax3,hxCmin1,hxCmin2,&
   hxCmin3,hxCmax1,hxCmax2,hxCmax3,hxRmin1,hxRmin2,hxRmin3,hxRmax1,hxRmax2,&
   hxRmax3
integer   :: kxLLmin1,kxLLmin2,kxLLmin3,kxLLmax1,kxLLmax2,kxLLmax3,kxLmin1,&
   kxLmin2,kxLmin3,kxLmax1,kxLmax2,kxLmax3,kxCmin1,kxCmin2,kxCmin3,kxCmax1,&
   kxCmax2,kxCmax3,kxRmin1,kxRmin2,kxRmin3,kxRmax1,kxRmax2,kxRmax3,kxRRmin1,&
   kxRRmin2,kxRRmin3,kxRRmax1,kxRRmax2,kxRRmax3
integer   :: iw
!--------------------------------------------------------------------------

ixOmin1=ixmin1-kr(idims,1);ixOmin2=ixmin2-kr(idims,2)
ixOmin3=ixmin3-kr(idims,3);ixOmax1=ixmax1+kr(idims,1)
ixOmax2=ixmax2+kr(idims,2);ixOmax3=ixmax3+kr(idims,3); !ixO[ixMmin1-1,ixMmax1+1]
ixLmin1=ixOmin1-kr(idims,1);ixLmin2=ixOmin2-kr(idims,2)
ixLmin3=ixOmin3-kr(idims,3);ixLmax1=ixOmax1-kr(idims,1)
ixLmax2=ixOmax2-kr(idims,2);ixLmax3=ixOmax3-kr(idims,3); !ixL[ixMmin1-2,ixMmax1]
ixLLmin1=ixLmin1-kr(idims,1);ixLLmin2=ixLmin2-kr(idims,2)
ixLLmin3=ixLmin3-kr(idims,3);ixLLmax1=ixLmax1-kr(idims,1)
ixLLmax2=ixLmax2-kr(idims,2);ixLLmax3=ixLmax3-kr(idims,3); !ixLL[ixMmin1-3,ixMmax1-1]
ixRmin1=ixOmin1+kr(idims,1);ixRmin2=ixOmin2+kr(idims,2)
ixRmin3=ixOmin3+kr(idims,3);ixRmax1=ixOmax1+kr(idims,1)
ixRmax2=ixOmax2+kr(idims,2);ixRmax3=ixOmax3+kr(idims,3); !ixR=[iMmin1,ixMmax+2]
ixRRmin1=ixRmin1+kr(idims,1);ixRRmin2=ixRmin2+kr(idims,2)
ixRRmin3=ixRmin3+kr(idims,3);ixRRmax1=ixRmax1+kr(idims,1)
ixRRmax2=ixRmax2+kr(idims,2);ixRRmax3=ixRmax3+kr(idims,3); !ixRR=[iMmin1+1,ixMmax+3]

hxCmin1=ixOmin1;hxCmin2=ixOmin2;hxCmin3=ixOmin3;hxCmax1=ixmax1;hxCmax2=ixmax2
hxCmax3=ixmax3; !hxC = [ixMmin-1,ixMmax]
hxLmin1=hxCmin1-kr(idims,1);hxLmin2=hxCmin2-kr(idims,2)
hxLmin3=hxCmin3-kr(idims,3);hxLmax1=hxCmax1-kr(idims,1)
hxLmax2=hxCmax2-kr(idims,2);hxLmax3=hxCmax3-kr(idims,3); !hxL = [ixMmin-2,ixMmax-1]
hxRmin1=hxCmin1+kr(idims,1);hxRmin2=hxCmin2+kr(idims,2)
hxRmin3=hxCmin3+kr(idims,3);hxRmax1=hxCmax1+kr(idims,1)
hxRmax2=hxCmax2+kr(idims,2);hxRmax3=hxCmax3+kr(idims,3); !hxR = [ixMmin,ixMmax+1]

kxCmin1=ixLLmin1;kxCmin2=ixLLmin2;kxCmin3=ixLLmin3; kxCmax1=ixRmax1
kxCmax2=ixRmax2;kxCmax3=ixRmax3; !kxC=[iMmin1-3,ixMmax1+2]
kxLmin1=kxCmin1-kr(idims,1);kxLmin2=kxCmin2-kr(idims,2)
kxLmin3=kxCmin3-kr(idims,3);kxLmax1=kxCmax1-kr(idims,1)
kxLmax2=kxCmax2-kr(idims,2);kxLmax3=kxCmax3-kr(idims,3); !kxL=[iMmin1-4,ixMmax1+1]
kxRmin1=kxCmin1+kr(idims,1);kxRmin2=kxCmin2+kr(idims,2)
kxRmin3=kxCmin3+kr(idims,3);kxRmax1=kxCmax1+kr(idims,1)
kxRmax2=kxCmax2+kr(idims,2);kxRmax3=kxCmax3+kr(idims,3); !kxR=[iMmin1-2,ixMmax1+3]

lxCmin1=ixLLmin1-kr(idims,1);lxCmin2=ixLLmin2-kr(idims,2)
lxCmin3=ixLLmin3-kr(idims,3);lxCmax1=ixRRmax1;lxCmax2=ixRRmax2
lxCmax3=ixRRmax3; !ixC=[iMmin1-4,ixMmax1+3]
lxRmin1=lxCmin1+kr(idims,1);lxRmin2=lxCmin2+kr(idims,2)
lxRmin3=lxCmin3+kr(idims,3);lxRmax1=lxCmax1+kr(idims,1)
lxRmax2=lxCmax2+kr(idims,2);lxRmax3=lxCmax3+kr(idims,3); !lxR=[iMmin1-3,ixMmax1+4]

do iw = 1, nwflux

   dwC(lxCmin1:lxCmax1,lxCmin2:lxCmax2,lxCmin3:lxCmax3)=w(lxRmin1:lxRmax1,&
      lxRmin2:lxRmax2,lxRmin3:lxRmax3,iw)-w(lxCmin1:lxCmax1,lxCmin2:lxCmax2,&
      lxCmin3:lxCmax3,iw)
   ! Eq. 64,  Miller and Colella 2002, JCP 183, 26 
   d2wC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3)=half*(w&
      (kxRmin1:kxRmax1,kxRmin2:kxRmax2,kxRmin3:kxRmax3,iw)-w(kxLmin1:kxLmax1,&
      kxLmin2:kxLmax2,kxLmin3:kxLmax3,iw))
   where(dwC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3)*dwC&
      (kxLmin1:kxLmax1,kxLmin2:kxLmax2,kxLmin3:kxLmax3)>zero)
      ! Store the sign of dwC in wMin
      wMin(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3)&
         = sign(one,d2wC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3))
      ! Eq. 65,  Miller and Colella 2002, JCP 183, 26 
      ldw(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3)= &
         wMin(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3)*min(dabs(d2wC&
         (kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3)),&
         2.0d0*dabs(dwC(kxLmin1:kxLmax1,kxLmin2:kxLmax2,kxLmin3:kxLmax3)),&
         2.0d0*dabs(dwC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3)))
   elsewhere
      ldw(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3)=zero
   endwhere

   ! Eq. 66,  Miller and Colella 2002, JCP 183, 26 
   wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)=wLC&
      (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)+half*dwC&
      (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)+(ldw(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3)-ldw(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
      ixRmin3:ixRmax3))/6.0d0

   wRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3,iw)=wRC&
      (ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3,iw)-(half*dwC&
      (ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3)-(ldw(ixLmin1:ixLmax1,&
      ixLmin2:ixLmax2,ixLmin3:ixLmax3)-ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      ixOmin3:ixOmax3))/6.0d0)


   ! Eq. B9, page 217, Mignone et al 2005, ApJS
   where((wRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3,&
      iw)-wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
      iw))*(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
      iw)-wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw))<=zero)
      wRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3,iw)&
         =wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)
      wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
         =wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)
   end where

   wMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=(wLC(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)-wRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
      ixLmin3:ixLmax3,iw))*(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      ixOmin3:ixOmax3,iw)-(wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      ixOmin3:ixOmax3,iw)+wRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3,&
      iw))/2.0d0)
   wMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=(wLC(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)-wRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
      ixLmin3:ixLmax3,iw))**2.0d0/6.0d0
   tmp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3)=wRC(ixLmin1:ixLmax1,&
      ixLmin2:ixLmax2,ixLmin3:ixLmax3,iw)
   ! Eq. B10, page 218, Mignone et al 2005, ApJS
   where(wMax(hxRmin1:hxRmax1,hxRmin2:hxRmax2,hxRmin3:hxRmax3)&
      >wMin(hxRmin1:hxRmax1,hxRmin2:hxRmax2,hxRmin3:hxRmax3))
      wRC(hxCmin1:hxCmax1,hxCmin2:hxCmax2,hxCmin3:hxCmax3,iw)&
         = 3.0d0*wCT(hxRmin1:hxRmax1,hxRmin2:hxRmax2,hxRmin3:hxRmax3,&
         iw)-2.0d0*wLC(hxRmin1:hxRmax1,hxRmin2:hxRmax2,hxRmin3:hxRmax3,iw)
   endwhere
   ! Eq. B11, page 218, Mignone et al 2005, ApJS
   where(wMax(hxCmin1:hxCmax1,hxCmin2:hxCmax2,hxCmin3:hxCmax3)&
      <-wMin(hxCmin1:hxCmax1,hxCmin2:hxCmax2,hxCmin3:hxCmax3))
      wLC(hxCmin1:hxCmax1,hxCmin2:hxCmax2,hxCmin3:hxCmax3,iw)&
         = 3.0d0*wCT(hxCmin1:hxCmax1,hxCmin2:hxCmax2,hxCmin3:hxCmax3,&
         iw)-2.0d0*tmp(hxLmin1:hxLmax1,hxLmin2:hxLmax2,hxLmin3:hxLmax3)
   endwhere

end do

if (flatsh) call mpistop("PPMlimiter: flatsh no longer supported")
if (flatcd) call mpistop("PPMlimiter: flatcd no longer supported")

end subroutine PPMlimiter
!============================================================================
subroutine extremaq(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
   ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,q,nshift,qMax,qMin)

use mod_amrvacdef

integer,intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(in) :: q(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)
integer,intent(in)           :: nshift

double precision, intent(out) :: qMax(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3),qMin(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

integer           :: ixsmin1,ixsmin2,ixsmin3,ixsmax1,ixsmax2,ixsmax3,ixsRmin1,&
   ixsRmin2,ixsRmin3,ixsRmax1,ixsRmax2,ixsRmax3,ixsLmin1,ixsLmin2,ixsLmin3,&
   ixsLmax1,ixsLmax2,ixsLmax3,idims,jdims,kdims,ishift,i,j 
!-------------------------------------------------------------------------
do ishift=1,nshift
 idims=1
 ixsRmin1=ixOmin1+ishift*kr(idims,1);ixsRmin2=ixOmin2+ishift*kr(idims,2)
 ixsRmin3=ixOmin3+ishift*kr(idims,3);ixsRmax1=ixOmax1+ishift*kr(idims,1)
 ixsRmax2=ixOmax2+ishift*kr(idims,2);ixsRmax3=ixOmax3+ishift*kr(idims,3);
 ixsLmin1=ixOmin1-ishift*kr(idims,1);ixsLmin2=ixOmin2-ishift*kr(idims,2)
 ixsLmin3=ixOmin3-ishift*kr(idims,3);ixsLmax1=ixOmax1-ishift*kr(idims,1)
 ixsLmax2=ixOmax2-ishift*kr(idims,2);ixsLmax3=ixOmax3-ishift*kr(idims,3);
 if (ishift==1) then
   qMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=max(q&
      (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),q(ixsRmin1:ixsRmax1,&
      ixsRmin2:ixsRmax2,ixsRmin3:ixsRmax3),q(ixsLmin1:ixsLmax1,&
      ixsLmin2:ixsLmax2,ixsLmin3:ixsLmax3))
   qMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=min(q&
      (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),q(ixsRmin1:ixsRmax1,&
      ixsRmin2:ixsRmax2,ixsRmin3:ixsRmax3),q(ixsLmin1:ixsLmax1,&
      ixsLmin2:ixsLmax2,ixsLmin3:ixsLmax3))
 else
   qMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=max(qMax&
      (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),q(ixsRmin1:ixsRmax1,&
      ixsRmin2:ixsRmax2,ixsRmin3:ixsRmax3),q(ixsLmin1:ixsLmax1,&
      ixsLmin2:ixsLmax2,ixsLmin3:ixsLmax3))
   qMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=min(qMin&
      (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),q(ixsRmin1:ixsRmax1,&
      ixsRmin2:ixsRmax2,ixsRmin3:ixsRmax3),q(ixsLmin1:ixsLmax1,&
      ixsLmin2:ixsLmax2,ixsLmin3:ixsLmax3))
 end if
 
 idims=1
 jdims=idims+1
 do i=-1,1
   ixsmin1=ixOmin1+i*ishift*kr(idims,1);ixsmin2=ixOmin2+i*ishift*kr(idims,2)
   ixsmin3=ixOmin3+i*ishift*kr(idims,3);ixsmax1=ixOmax1+i*ishift*kr(idims,1)
   ixsmax2=ixOmax2+i*ishift*kr(idims,2);ixsmax3=ixOmax3+i*ishift*kr(idims,3);
   ixsRmin1=ixsmin1+ishift*kr(jdims,1);ixsRmin2=ixsmin2+ishift*kr(jdims,2)
   ixsRmin3=ixsmin3+ishift*kr(jdims,3);ixsRmax1=ixsmax1+ishift*kr(jdims,1)
   ixsRmax2=ixsmax2+ishift*kr(jdims,2);ixsRmax3=ixsmax3+ishift*kr(jdims,3);
   ixsLmin1=ixsmin1-ishift*kr(jdims,1);ixsLmin2=ixsmin2-ishift*kr(jdims,2)
   ixsLmin3=ixsmin3-ishift*kr(jdims,3);ixsLmax1=ixsmax1-ishift*kr(jdims,1)
   ixsLmax2=ixsmax2-ishift*kr(jdims,2);ixsLmax3=ixsmax3-ishift*kr(jdims,3);
   qMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=max(qMax&
      (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),q(ixsRmin1:ixsRmax1,&
      ixsRmin2:ixsRmax2,ixsRmin3:ixsRmax3),q(ixsLmin1:ixsLmax1,&
      ixsLmin2:ixsLmax2,ixsLmin3:ixsLmax3))
   qMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=min(qMin&
      (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),q(ixsRmin1:ixsRmax1,&
      ixsRmin2:ixsRmax2,ixsRmin3:ixsRmax3),q(ixsLmin1:ixsLmax1,&
      ixsLmin2:ixsLmax2,ixsLmin3:ixsLmax3))
 end do

 
 idims=1
 jdims=idims+1
 kdims=jdims+1
 do i=-1,1
   ixsmin1=ixOmin1+i*ishift*kr(idims,1);ixsmin2=ixOmin2+i*ishift*kr(idims,2)
   ixsmin3=ixOmin3+i*ishift*kr(idims,3);ixsmax1=ixOmax1+i*ishift*kr(idims,1)
   ixsmax2=ixOmax2+i*ishift*kr(idims,2);ixsmax3=ixOmax3+i*ishift*kr(idims,3);
   do j=-1,1
      ixsmin1=ixOmin1+j*ishift*kr(jdims,1)
      ixsmin2=ixOmin2+j*ishift*kr(jdims,2)
      ixsmin3=ixOmin3+j*ishift*kr(jdims,3)
      ixsmax1=ixOmax1+j*ishift*kr(jdims,1)
      ixsmax2=ixOmax2+j*ishift*kr(jdims,2)
      ixsmax3=ixOmax3+j*ishift*kr(jdims,3);
      ixsRmin1=ixsmin1+ishift*kr(kdims,1);ixsRmin2=ixsmin2+ishift*kr(kdims,2)
      ixsRmin3=ixsmin3+ishift*kr(kdims,3);ixsRmax1=ixsmax1+ishift*kr(kdims,1)
      ixsRmax2=ixsmax2+ishift*kr(kdims,2);ixsRmax3=ixsmax3+ishift*kr(kdims,3);
      ixsLmin1=ixsmin1-ishift*kr(kdims,1);ixsLmin2=ixsmin2-ishift*kr(kdims,2)
      ixsLmin3=ixsmin3-ishift*kr(kdims,3);ixsLmax1=ixsmax1-ishift*kr(kdims,1)
      ixsLmax2=ixsmax2-ishift*kr(kdims,2);ixsLmax3=ixsmax3-ishift*kr(kdims,3);
      qMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
         =max(qMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
         q(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2,ixsRmin3:ixsRmax3),&
         q(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2,ixsLmin3:ixsLmax3))
      qMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
         =min(qMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
         q(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2,ixsRmin3:ixsRmax3),&
         q(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2,ixsLmin3:ixsLmax3))
   end do
 end do

enddo

end subroutine  extremaq
!=============================================================================
subroutine extremaw(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
   ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,nshift,wMax,wMin)

use mod_amrvacdef

integer,intent(in)            :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw)
integer,intent(in)            :: nshift

double precision, intent(out) :: wMax(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nwflux),wMin(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nwflux)

integer          :: ixsmin1,ixsmin2,ixsmin3,ixsmax1,ixsmax2,ixsmax3,ixsRmin1,&
   ixsRmin2,ixsRmin3,ixsRmax1,ixsRmax2,ixsRmax3,ixsLmin1,ixsLmin2,ixsLmin3,&
   ixsLmax1,ixsLmax2,ixsLmax3,idims,jdims,kdims,ishift,i,j
!-------------------------------------------------------------------------
do ishift=1,nshift
 idims=1
 ixsRmin1=ixOmin1+ishift*kr(idims,1);ixsRmin2=ixOmin2+ishift*kr(idims,2)
 ixsRmin3=ixOmin3+ishift*kr(idims,3);ixsRmax1=ixOmax1+ishift*kr(idims,1)
 ixsRmax2=ixOmax2+ishift*kr(idims,2);ixsRmax3=ixOmax3+ishift*kr(idims,3);
 ixsLmin1=ixOmin1-ishift*kr(idims,1);ixsLmin2=ixOmin2-ishift*kr(idims,2)
 ixsLmin3=ixOmin3-ishift*kr(idims,3);ixsLmax1=ixOmax1-ishift*kr(idims,1)
 ixsLmax2=ixOmax2-ishift*kr(idims,2);ixsLmax3=ixOmax3-ishift*kr(idims,3);
 if (ishift==1) then
    wMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:nwflux)&
       = max(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:nwflux),&
       w(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2,ixsRmin3:ixsRmax3,1:nwflux),&
       w(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2,ixsLmin3:ixsLmax3,1:nwflux))
    wMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:nwflux)&
       = min(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:nwflux),&
       w(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2,ixsRmin3:ixsRmax3,1:nwflux),&
       w(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2,ixsLmin3:ixsLmax3,1:nwflux))
 else
    wMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:nwflux)&
       = max(wMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:nwflux),&
       w(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2,ixsRmin3:ixsRmax3,1:nwflux),&
       w(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2,ixsLmin3:ixsLmax3,1:nwflux))
    wMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:nwflux)&
       = min(wMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:nwflux),&
       w(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2,ixsRmin3:ixsRmax3,1:nwflux),&
       w(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2,ixsLmin3:ixsLmax3,1:nwflux))
 end if
 
 idims=1
 jdims=idims+1
 do i=-1,1
   ixsmin1=ixOmin1+i*ishift*kr(idims,1);ixsmin2=ixOmin2+i*ishift*kr(idims,2)
   ixsmin3=ixOmin3+i*ishift*kr(idims,3);ixsmax1=ixOmax1+i*ishift*kr(idims,1)
   ixsmax2=ixOmax2+i*ishift*kr(idims,2);ixsmax3=ixOmax3+i*ishift*kr(idims,3);
   ixsRmin1=ixsmin1+ishift*kr(jdims,1);ixsRmin2=ixsmin2+ishift*kr(jdims,2)
   ixsRmin3=ixsmin3+ishift*kr(jdims,3);ixsRmax1=ixsmax1+ishift*kr(jdims,1)
   ixsRmax2=ixsmax2+ishift*kr(jdims,2);ixsRmax3=ixsmax3+ishift*kr(jdims,3);
   ixsLmin1=ixsmin1-ishift*kr(jdims,1);ixsLmin2=ixsmin2-ishift*kr(jdims,2)
   ixsLmin3=ixsmin3-ishift*kr(jdims,3);ixsLmax1=ixsmax1-ishift*kr(jdims,1)
   ixsLmax2=ixsmax2-ishift*kr(jdims,2);ixsLmax3=ixsmax3-ishift*kr(jdims,3);
   wMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:nwflux)= &
     max(wMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:nwflux),&
        w(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2,ixsRmin3:ixsRmax3,1:nwflux),&
        w(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2,ixsLmin3:ixsLmax3,1:nwflux))
   wMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:nwflux)= &
     min(wMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:nwflux),&
        w(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2,ixsRmin3:ixsRmax3,1:nwflux),&
        w(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2,ixsLmin3:ixsLmax3,1:nwflux))
 end do

 
 idims=1
 jdims=idims+1
 kdims=jdims+1
 do i=-1,1
   ixsmin1=ixOmin1+i*ishift*kr(idims,1);ixsmin2=ixOmin2+i*ishift*kr(idims,2)
   ixsmin3=ixOmin3+i*ishift*kr(idims,3);ixsmax1=ixOmax1+i*ishift*kr(idims,1)
   ixsmax2=ixOmax2+i*ishift*kr(idims,2);ixsmax3=ixOmax3+i*ishift*kr(idims,3);
   do j=-1,1
      ixsmin1=ixOmin1+j*ishift*kr(jdims,1)
      ixsmin2=ixOmin2+j*ishift*kr(jdims,2)
      ixsmin3=ixOmin3+j*ishift*kr(jdims,3)
      ixsmax1=ixOmax1+j*ishift*kr(jdims,1)
      ixsmax2=ixOmax2+j*ishift*kr(jdims,2)
      ixsmax3=ixOmax3+j*ishift*kr(jdims,3);
      ixsRmin1=ixsmin1+ishift*kr(kdims,1);ixsRmin2=ixsmin2+ishift*kr(kdims,2)
      ixsRmin3=ixsmin3+ishift*kr(kdims,3);ixsRmax1=ixsmax1+ishift*kr(kdims,1)
      ixsRmax2=ixsmax2+ishift*kr(kdims,2);ixsRmax3=ixsmax3+ishift*kr(kdims,3);
      ixsLmin1=ixsmin1-ishift*kr(kdims,1);ixsLmin2=ixsmin2-ishift*kr(kdims,2)
      ixsLmin3=ixsmin3-ishift*kr(kdims,3);ixsLmax1=ixsmax1-ishift*kr(kdims,1)
      ixsLmax2=ixsmax2-ishift*kr(kdims,2);ixsLmax3=ixsmax3-ishift*kr(kdims,3);
      wMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:nwflux)= &
       max(wMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:nwflux),&
          w(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2,ixsRmin3:ixsRmax3,1:nwflux),&
          w(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2,ixsLmin3:ixsLmax3,1:nwflux))
      wMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:nwflux)= &
       min(wMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:nwflux),&
          w(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2,ixsRmin3:ixsRmax3,1:nwflux),&
          w(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2,ixsLmin3:ixsLmax3,1:nwflux))
   end do
 end do

enddo

end subroutine  extremaw
!=============================================================================
subroutine extremaa(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
   ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,a,nshift,aMin)

use mod_amrvacdef

integer,intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(in) :: a(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
integer,intent(in)           :: nshift

double precision, intent(out) :: aMin(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
   ixGlo3:ixGhi3)

integer          :: ixsmin1,ixsmin2,ixsmin3,ixsmax1,ixsmax2,ixsmax3,ixsRmin1,&
   ixsRmin2,ixsRmin3,ixsRmax1,ixsRmax2,ixsRmax3,ixsLmin1,ixsLmin2,ixsLmin3,&
   ixsLmax1,ixsLmax2,ixsLmax3,idims,jdims,kdims,ishift,i,j
!-------------------------------------------------------------------------
do ishift=1,nshift
  idims=1
  ixsRmin1=ixOmin1+ishift*kr(idims,1);ixsRmin2=ixOmin2+ishift*kr(idims,2)
  ixsRmin3=ixOmin3+ishift*kr(idims,3);ixsRmax1=ixOmax1+ishift*kr(idims,1)
  ixsRmax2=ixOmax2+ishift*kr(idims,2);ixsRmax3=ixOmax3+ishift*kr(idims,3);
  ixsLmin1=ixOmin1-ishift*kr(idims,1);ixsLmin2=ixOmin2-ishift*kr(idims,2)
  ixsLmin3=ixOmin3-ishift*kr(idims,3);ixsLmax1=ixOmax1-ishift*kr(idims,1)
  ixsLmax2=ixOmax2-ishift*kr(idims,2);ixsLmax3=ixOmax3-ishift*kr(idims,3);
  aMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=min(a&
     (ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2,ixsRmin3:ixsRmax3),&
     a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),a(ixsLmin1:ixsLmax1,&
     ixsLmin2:ixsLmax2,ixsLmin3:ixsLmax3))
  
  idims=1
  jdims=idims+1
  do i=-1,1
    ixsmin1=ixOmin1+i*ishift*kr(idims,1);ixsmin2=ixOmin2+i*ishift*kr(idims,2)
    ixsmin3=ixOmin3+i*ishift*kr(idims,3);ixsmax1=ixOmax1+i*ishift*kr(idims,1)
    ixsmax2=ixOmax2+i*ishift*kr(idims,2);ixsmax3=ixOmax3+i*ishift*kr(idims,3);
    ixsRmin1=ixsmin1+ishift*kr(jdims,1);ixsRmin2=ixsmin2+ishift*kr(jdims,2)
    ixsRmin3=ixsmin3+ishift*kr(jdims,3);ixsRmax1=ixsmax1+ishift*kr(jdims,1)
    ixsRmax2=ixsmax2+ishift*kr(jdims,2);ixsRmax3=ixsmax3+ishift*kr(jdims,3);
    ixsLmin1=ixsmin1-ishift*kr(jdims,1);ixsLmin2=ixsmin2-ishift*kr(jdims,2)
    ixsLmin3=ixsmin3-ishift*kr(jdims,3);ixsLmax1=ixsmax1-ishift*kr(jdims,1)
    ixsLmax2=ixsmax2-ishift*kr(jdims,2);ixsLmax3=ixsmax3-ishift*kr(jdims,3);
    aMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=min(aMin&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),a(ixsRmin1:ixsRmax1,&
       ixsRmin2:ixsRmax2,ixsRmin3:ixsRmax3),a(ixsLmin1:ixsLmax1,&
       ixsLmin2:ixsLmax2,ixsLmin3:ixsLmax3))
  end do
 
  
  idims=1
  jdims=idims+1
  kdims=jdims+1
  do i=-1,1
   ixsmin1=ixOmin1+i*ishift*kr(idims,1);ixsmin2=ixOmin2+i*ishift*kr(idims,2)
   ixsmin3=ixOmin3+i*ishift*kr(idims,3);ixsmax1=ixOmax1+i*ishift*kr(idims,1)
   ixsmax2=ixOmax2+i*ishift*kr(idims,2);ixsmax3=ixOmax3+i*ishift*kr(idims,3);
   do j=-1,1
      ixsmin1=ixOmin1+j*ishift*kr(jdims,1)
      ixsmin2=ixOmin2+j*ishift*kr(jdims,2)
      ixsmin3=ixOmin3+j*ishift*kr(jdims,3)
      ixsmax1=ixOmax1+j*ishift*kr(jdims,1)
      ixsmax2=ixOmax2+j*ishift*kr(jdims,2)
      ixsmax3=ixOmax3+j*ishift*kr(jdims,3);
      ixsRmin1=ixsmin1+ishift*kr(kdims,1);ixsRmin2=ixsmin2+ishift*kr(kdims,2)
      ixsRmin3=ixsmin3+ishift*kr(kdims,3);ixsRmax1=ixsmax1+ishift*kr(kdims,1)
      ixsRmax2=ixsmax2+ishift*kr(kdims,2);ixsRmax3=ixsmax3+ishift*kr(kdims,3);
      ixsLmin1=ixsmin1-ishift*kr(kdims,1);ixsLmin2=ixsmin2-ishift*kr(kdims,2)
      ixsLmin3=ixsmin3-ishift*kr(kdims,3);ixsLmax1=ixsmax1-ishift*kr(kdims,1)
      ixsLmax2=ixsmax2-ishift*kr(kdims,2);ixsLmax3=ixsmax3-ishift*kr(kdims,3);
      aMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
         =min(aMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
         a(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2,ixsRmin3:ixsRmax3),&
         a(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2,ixsLmin3:ixsLmax3))
   end do
  end do
 
end do

end subroutine extremaa
!=============================================================================
  
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

subroutine MP5limiter(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,&
   iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,idims,w,wLC,wRC)
! MP5 limiter from Suresh & Huynh 1997
! Following the convention of Mignone et al. 2010.
! Needs at least three ghost cells.  Set dixB=3.

use mod_amrvacdef

integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3, iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3, idims
double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw)
double precision, intent(inout) :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw),wLC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nw) 

! .. local ..
integer                         :: iLmmin1,iLmmin2,iLmmin3,iLmmax1,iLmmax2,&
   iLmmax3, iLmmmin1,iLmmmin2,iLmmmin3,iLmmmax1,iLmmmax2,iLmmmax3, iLpmin1,&
   iLpmin2,iLpmin3,iLpmax1,iLpmax2,iLpmax3, iLppmin1,iLppmin2,iLppmin3,&
   iLppmax1,iLppmax2,iLppmax3, iLpppmin1,iLpppmin2,iLpppmin3,iLpppmax1,&
   iLpppmax2,iLpppmax3
integer                         :: idmin1,idmin2,idmin3,idmax1,idmax2,idmax3,&
    idpmin1,idpmin2,idpmin3,idpmax1,idpmax2,idpmax3, idppmin1,idppmin2,&
   idppmin3,idppmax1,idppmax2,idppmax3, idmmin1,idmmin2,idmmin3,idmmax1,&
   idmmax2,idmmax3, iemin1,iemin2,iemin3,iemax1,iemax2,iemax3, iemmin1,&
   iemmin2,iemmin3,iemmax1,iemmax2,iemmax3, iepmin1,iepmin2,iepmin3,iepmax1,&
   iepmax2,iepmax3, ieppmin1,ieppmin2,ieppmin3,ieppmax1,ieppmax2,ieppmax3
integer                         :: iw
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nw)  :: f, fmp, fmin, fmax, ful, dm4, d, fmd, flc, flim
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3) :: tmp, tmp2, tmp3, a, b, c
double precision, parameter     :: eps=0.0d0, alpha=4.0d0
!double precision                :: alpha
!----------------------------------------------------------------------------

! Variable alpha:
!alpha = float(nstep)/courantpar - one

! Left side:
! range to process:
!iLmin^D=ixmin^D-kr(idims,^D);iLmax^D=ixmax^D;

!{#IFDEF HALL
   ! For Hall, we need one more reconstructed layer since currents are computed in getflux:
   ! also add one ghost zone!
!   {iL^L=iL^L^LADD1;}
!}

! iL^L holds the indices of interfaces to reconstruct to.  Convention is that a center index holds the _right-side_ interface.  

iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
iLmmin3=iLmin3-kr(idims,3);iLmmax1=iLmax1-kr(idims,1)
iLmmax2=iLmax2-kr(idims,2);iLmmax3=iLmax3-kr(idims,3);
iLmmmin1=iLmmin1-kr(idims,1);iLmmmin2=iLmmin2-kr(idims,2)
iLmmmin3=iLmmin3-kr(idims,3);iLmmmax1=iLmmax1-kr(idims,1)
iLmmmax2=iLmmax2-kr(idims,2);iLmmmax3=iLmmax3-kr(idims,3);
iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
iLpmin3=iLmin3+kr(idims,3);iLpmax1=iLmax1+kr(idims,1)
iLpmax2=iLmax2+kr(idims,2);iLpmax3=iLmax3+kr(idims,3);
iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
iLppmin3=iLpmin3+kr(idims,3);iLppmax1=iLpmax1+kr(idims,1)
iLppmax2=iLpmax2+kr(idims,2);iLppmax3=iLpmax3+kr(idims,3);

f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = 1.0d0&
   /60.0d0 * (2.0d0* w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,iLmmmin3:iLmmmax3,&
   1:nwflux) - 13.0d0* w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iLmmin3:iLmmax3,&
   1:nwflux) + 47.0d0* w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
   1:nwflux) + 27.0d0* w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,&
   1:nwflux) - 3.0d0*  w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
   iLppmin3:iLppmax3,1:nwflux))

! get fmp and ful:
do iw=1,nwflux
   a(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = w(iLpmin1:iLpmax1,&
      iLpmin2:iLpmax2,iLpmin3:iLpmax3,iw)-w(iLmin1:iLmax1,iLmin2:iLmax2,&
      iLmin3:iLmax3,iw)
   b(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = alpha*(w(iLmin1:iLmax1,&
      iLmin2:iLmax2,iLmin3:iLmax3,iw)-w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
      iLmmin3:iLmmax3,iw))
   call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,iLmin2,&
      iLmin3,iLmax1,iLmax2,iLmax3,a,b,tmp)
   fmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,iw) = w(iLmin1:iLmax1,&
      iLmin2:iLmax2,iLmin3:iLmax3,iw) + tmp(iLmin1:iLmax1,iLmin2:iLmax2,&
      iLmin3:iLmax3)
   ful(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,iw) = w(iLmin1:iLmax1,&
      iLmin2:iLmax2,iLmin3:iLmax3,iw) + b(iLmin1:iLmax1,iLmin2:iLmax2,&
      iLmin3:iLmax3)
end do ! iw loop

! get dm4:
idmax1=iLmax1;idmax2=iLmax2;idmax3=iLmax3; idmin1=iLmin1-kr(idims,1)
idmin2=iLmin2-kr(idims,2);idmin3=iLmin3-kr(idims,3);
idmmin1=idmin1-kr(idims,1);idmmin2=idmin2-kr(idims,2)
idmmin3=idmin3-kr(idims,3);idmmax1=idmax1-kr(idims,1)
idmmax2=idmax2-kr(idims,2);idmmax3=idmax3-kr(idims,3);
idpmin1=idmin1+kr(idims,1);idpmin2=idmin2+kr(idims,2)
idpmin3=idmin3+kr(idims,3);idpmax1=idmax1+kr(idims,1)
idpmax2=idmax2+kr(idims,2);idpmax3=idmax3+kr(idims,3);

iemax1=idmax1+kr(idims,1);iemax2=idmax2+kr(idims,2);iemax3=idmax3+kr(idims,3)
iemin1=idmin1;iemin2=idmin2;iemin3=idmin3;
iemmin1=iemin1-kr(idims,1);iemmin2=iemin2-kr(idims,2)
iemmin3=iemin3-kr(idims,3);iemmax1=iemax1-kr(idims,1)
iemmax2=iemax2-kr(idims,2);iemmax3=iemax3-kr(idims,3);
iepmin1=iemin1+kr(idims,1);iepmin2=iemin2+kr(idims,2)
iepmin3=iemin3+kr(idims,3);iepmax1=iemax1+kr(idims,1)
iepmax2=iemax2+kr(idims,2);iepmax3=iemax3+kr(idims,3);

d(iemin1:iemax1,iemin2:iemax2,iemin3:iemax3,1:nwflux) = w(iepmin1:iepmax1,&
   iepmin2:iepmax2,iepmin3:iepmax3,1:nwflux)-2.0d0*w(iemin1:iemax1,&
   iemin2:iemax2,iemin3:iemax3,1:nwflux)+w(iemmin1:iemmax1,iemmin2:iemmax2,&
   iemmin3:iemmax3,1:nwflux)

do iw=1,nwflux
   a(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = 4.0d0*d(idmin1:idmax1,&
      idmin2:idmax2,idmin3:idmax3,iw)-d(idpmin1:idpmax1,idpmin2:idpmax2,&
      idpmin3:idpmax3,iw)
   b(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = 4.0d0*d(idpmin1:idpmax1,&
      idpmin2:idpmax2,idpmin3:idpmax3,iw)-d(idmin1:idmax1,idmin2:idmax2,&
      idmin3:idmax3,iw)
   call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,idmin2,&
      idmin3,idmax1,idmax2,idmax3,a,b,tmp)
   a(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = d(idmin1:idmax1,&
      idmin2:idmax2,idmin3:idmax3,iw)
   b(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = d(idpmin1:idpmax1,&
      idpmin2:idpmax2,idpmin3:idpmax3,iw)
   call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,idmin2,&
      idmin3,idmax1,idmax2,idmax3,a,b,tmp2)
   call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,idmin2,&
      idmin3,idmax1,idmax2,idmax3,tmp,tmp2,tmp3)
   dm4(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3,iw) = tmp3(idmin1:idmax1,&
      idmin2:idmax2,idmin3:idmax3)
end do

! get fmd:
fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = (w(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)+w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
   iLpmin3:iLpmax3,1:nwflux))/2.0d0-dm4(iLmin1:iLmax1,iLmin2:iLmax2,&
   iLmin3:iLmax3,1:nwflux)/2.0d0

!get flc: 
flc(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = half*(3.0d0*w&
   (iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) - w(iLmmin1:iLmmax1,&
   iLmmin2:iLmmax2,iLmmin3:iLmmax3,1:nwflux)) + 4.0d0/3.0d0*dm4&
   (iLmmin1:iLmmax1,iLmmin2:iLmmax2,iLmmin3:iLmmax3,1:nwflux)

fmin(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = &
   max(min(w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),&
   w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux),&
   fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)),&
   min(w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),&
   ful(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),flc(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)))

fmax(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = &
   min(max(w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),&
   w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux),&
   fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)),&
   max(w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),&
   ful(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),flc(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)))

do iw=1,nwflux
   a(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = fmin(iLmin1:iLmax1,&
      iLmin2:iLmax2,iLmin3:iLmax3,iw)
   b(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = f(iLmin1:iLmax1,&
      iLmin2:iLmax2,iLmin3:iLmax3,iw)
   c(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = fmax(iLmin1:iLmax1,&
      iLmin2:iLmax2,iLmin3:iLmax3,iw)
   call median(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,iLmin2,&
      iLmin3,iLmax1,iLmax2,iLmax3,a,b,c,tmp)
   flim(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,iw) = tmp(iLmin1:iLmax1,&
      iLmin2:iLmax2,iLmin3:iLmax3)
end do

! check case
where ((f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)-w(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux))*(f(iLmin1:iLmax1,iLmin2:iLmax2,&
   iLmin3:iLmax3,1:nwflux)-fmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
   1:nwflux)) .le. eps)
   wLC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) &
      = f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)
elsewhere
   wLC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) &
      = flim(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)
end where

! Right side:
! the interpolation from the right is obtained when the left-hand formula is applied to
! data mirrored about the interface.  
! thus substitute: 
! i-2 -> i+3
! i-1 -> i+2
! i   -> i+1
! i+1 -> i
! i+2 -> i-1

iLpppmin1=iLppmin1+kr(idims,1);iLpppmin2=iLppmin2+kr(idims,2)
iLpppmin3=iLppmin3+kr(idims,3);iLpppmax1=iLppmax1+kr(idims,1)
iLpppmax2=iLppmax2+kr(idims,2);iLpppmax3=iLppmax3+kr(idims,3);

f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = 1.0d0&
   /60.0d0 * (2.0d0* w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
   iLpppmin3:iLpppmax3,1:nwflux) - 13.0d0* w(iLppmin1:iLppmax1,&
   iLppmin2:iLppmax2,iLppmin3:iLppmax3,1:nwflux) + 47.0d0* w(iLpmin1:iLpmax1,&
   iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux) + 27.0d0* w(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) - 3.0d0*  w(iLmmin1:iLmmax1,&
   iLmmin2:iLmmax2,iLmmin3:iLmmax3,1:nwflux))

! get fmp and ful:
do iw=1,nwflux
   a(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = w(iLmin1:iLmax1,&
      iLmin2:iLmax2,iLmin3:iLmax3,iw)-w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
      iLpmin3:iLpmax3,iw)
   b(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = alpha*(w(iLpmin1:iLpmax1,&
      iLpmin2:iLpmax2,iLpmin3:iLpmax3,iw)-w(iLppmin1:iLppmax1,&
      iLppmin2:iLppmax2,iLppmin3:iLppmax3,iw))
   call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,iLmin2,&
      iLmin3,iLmax1,iLmax2,iLmax3,a,b,tmp)
   fmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,iw) = w(iLpmin1:iLpmax1,&
      iLpmin2:iLpmax2,iLpmin3:iLpmax3,iw) + tmp(iLmin1:iLmax1,iLmin2:iLmax2,&
      iLmin3:iLmax3)
   ful(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,iw) = w(iLpmin1:iLpmax1,&
      iLpmin2:iLpmax2,iLpmin3:iLpmax3,iw) + b(iLmin1:iLmax1,iLmin2:iLmax2,&
      iLmin3:iLmax3)
end do ! iw loop

! get dm4:
idmax1=iLmax1+kr(idims,1);idmax2=iLmax2+kr(idims,2);idmax3=iLmax3+kr(idims,3)
idmin1=iLmin1;idmin2=iLmin2;idmin3=iLmin3;
idmmin1=idmin1-kr(idims,1);idmmin2=idmin2-kr(idims,2)
idmmin3=idmin3-kr(idims,3);idmmax1=idmax1-kr(idims,1)
idmmax2=idmax2-kr(idims,2);idmmax3=idmax3-kr(idims,3);
idpmin1=idmin1+kr(idims,1);idpmin2=idmin2+kr(idims,2)
idpmin3=idmin3+kr(idims,3);idpmax1=idmax1+kr(idims,1)
idpmax2=idmax2+kr(idims,2);idpmax3=idmax3+kr(idims,3);

iemax1=idmax1;iemax2=idmax2;iemax3=idmax3; iemin1=idmin1-kr(idims,1)
iemin2=idmin2-kr(idims,2);iemin3=idmin3-kr(idims,3);
iemmin1=iemin1-kr(idims,1);iemmin2=iemin2-kr(idims,2)
iemmin3=iemin3-kr(idims,3);iemmax1=iemax1-kr(idims,1)
iemmax2=iemax2-kr(idims,2);iemmax3=iemax3-kr(idims,3);
iepmin1=iemin1+kr(idims,1);iepmin2=iemin2+kr(idims,2)
iepmin3=iemin3+kr(idims,3);iepmax1=iemax1+kr(idims,1)
iepmax2=iemax2+kr(idims,2);iepmax3=iemax3+kr(idims,3);
ieppmin1=iepmin1+kr(idims,1);ieppmin2=iepmin2+kr(idims,2)
ieppmin3=iepmin3+kr(idims,3);ieppmax1=iepmax1+kr(idims,1)
ieppmax2=iepmax2+kr(idims,2);ieppmax3=iepmax3+kr(idims,3);

d(iemin1:iemax1,iemin2:iemax2,iemin3:iemax3,1:nwflux) = w(iemin1:iemax1,&
   iemin2:iemax2,iemin3:iemax3,1:nwflux)-2.0d0*w(iepmin1:iepmax1,&
   iepmin2:iepmax2,iepmin3:iepmax3,1:nwflux)+w(ieppmin1:ieppmax1,&
   ieppmin2:ieppmax2,ieppmin3:ieppmax3,1:nwflux)

do iw=1,nwflux
   a(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = 4.0d0*d(idmin1:idmax1,&
      idmin2:idmax2,idmin3:idmax3,iw)-d(idmmin1:idmmax1,idmmin2:idmmax2,&
      idmmin3:idmmax3,iw)
   b(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = 4.0d0*d(idmmin1:idmmax1,&
      idmmin2:idmmax2,idmmin3:idmmax3,iw)-d(idmin1:idmax1,idmin2:idmax2,&
      idmin3:idmax3,iw)
   call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,idmin2,&
      idmin3,idmax1,idmax2,idmax3,a,b,tmp)
   a(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = d(idmin1:idmax1,&
      idmin2:idmax2,idmin3:idmax3,iw)
   b(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = d(idmmin1:idmmax1,&
      idmmin2:idmmax2,idmmin3:idmmax3,iw)
   call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,idmin2,&
      idmin3,idmax1,idmax2,idmax3,a,b,tmp2)
   call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,idmin2,&
      idmin3,idmax1,idmax2,idmax3,tmp,tmp2,tmp3)
   dm4(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3,iw) = tmp3(idmin1:idmax1,&
      idmin2:idmax2,idmin3:idmax3)
end do

! get fmd:
fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = (w(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)+w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
   iLpmin3:iLpmax3,1:nwflux))/2.0d0-dm4(iLmin1:iLmax1,iLmin2:iLmax2,&
   iLmin3:iLmax3,1:nwflux)/2.0d0

!get flc: 
flc(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = half*(3.0d0*w&
   (iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux) - &
   w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,iLppmin3:iLppmax3,1:nwflux)) + 4.0d0&
   /3.0d0*dm4(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux)

fmin(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = &
   max(min(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux),&
   w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),fmd(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)),min(w(iLpmin1:iLpmax1,&
   iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux),ful(iLmin1:iLmax1,iLmin2:iLmax2,&
   iLmin3:iLmax3,1:nwflux),flc(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
   1:nwflux)))

fmax(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = &
   min(max(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux),&
   w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),fmd(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)),max(w(iLpmin1:iLpmax1,&
   iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux),ful(iLmin1:iLmax1,iLmin2:iLmax2,&
   iLmin3:iLmax3,1:nwflux),flc(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
   1:nwflux)))

do iw=1,nwflux
   a(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = fmin(iLmin1:iLmax1,&
      iLmin2:iLmax2,iLmin3:iLmax3,iw)
   b(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = f(iLmin1:iLmax1,&
      iLmin2:iLmax2,iLmin3:iLmax3,iw)
   c(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = fmax(iLmin1:iLmax1,&
      iLmin2:iLmax2,iLmin3:iLmax3,iw)
   call median(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,iLmin2,&
      iLmin3,iLmax1,iLmax2,iLmax3,a,b,c,tmp)
   flim(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,iw) = tmp(iLmin1:iLmax1,&
      iLmin2:iLmax2,iLmin3:iLmax3)
end do

! check case
where ((f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)-w&
   (iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux))*(f&
   (iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)-fmp(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux))  .le. eps)
   wRC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) &
      = f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)
elsewhere
   wRC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) &
      = flim(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)
end where

end subroutine MP5limiter
!=============================================================================
subroutine MP5limitervar(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,idims,w,wLC,wRC)
! MP5 limiter from Suresh & Huynh 1997
! Following the convention of Mignone et al. 2010.
! Needs at least three ghost cells.  Set dixB=3.

use mod_amrvacdef

integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3, iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3, idims
double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)
double precision, intent(inout) :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3),wLC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) 

! .. local ..
integer                         :: iLmmin1,iLmmin2,iLmmin3,iLmmax1,iLmmax2,&
   iLmmax3, iLmmmin1,iLmmmin2,iLmmmin3,iLmmmax1,iLmmmax2,iLmmmax3, iLpmin1,&
   iLpmin2,iLpmin3,iLpmax1,iLpmax2,iLpmax3, iLppmin1,iLppmin2,iLppmin3,&
   iLppmax1,iLppmax2,iLppmax3, iLpppmin1,iLpppmin2,iLpppmin3,iLpppmax1,&
   iLpppmax2,iLpppmax3
integer                         :: idmin1,idmin2,idmin3,idmax1,idmax2,idmax3,&
    idpmin1,idpmin2,idpmin3,idpmax1,idpmax2,idpmax3, idppmin1,idppmin2,&
   idppmin3,idppmax1,idppmax2,idppmax3, idmmin1,idmmin2,idmmin3,idmmax1,&
   idmmax2,idmmax3, iemin1,iemin2,iemin3,iemax1,iemax2,iemax3, iemmin1,&
   iemmin2,iemmin3,iemmax1,iemmax2,iemmax3, iepmin1,iepmin2,iepmin3,iepmax1,&
   iepmax2,iepmax3, ieppmin1,ieppmin2,ieppmin3,ieppmax1,ieppmax2,ieppmax3
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)  :: f, fmp, fmin, fmax, ful, dm4, d, fmd, flc, flim
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)  :: wRCtmp, wLCtmp
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3) :: tmp, tmp2, tmp3, a, b, c
logical, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)       :: &
   flagL, flagR
double precision, parameter     :: eps=0.0d0, alpha=4.0d0
!double precision                :: alpha
!----------------------------------------------------------------------------

! Variable alpha:
!alpha = float(nstep)/courantpar - one

! Left side:
! range to process:
!iLmin^D=ixmin^D-kr(idims,^D);iLmax^D=ixmax^D;

!{#IFDEF HALL
   ! For Hall, we need one more reconstructed layer since currents are computed in getflux:
   ! also add one ghost zone!
!   {iL^L=iL^L^LADD1;}
!}

! iL^L holds the indices of interfaces to reconstruct to.  Convention is that a center index holds the _right-side_ interface.  

iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
iLmmin3=iLmin3-kr(idims,3);iLmmax1=iLmax1-kr(idims,1)
iLmmax2=iLmax2-kr(idims,2);iLmmax3=iLmax3-kr(idims,3);
iLmmmin1=iLmmin1-kr(idims,1);iLmmmin2=iLmmin2-kr(idims,2)
iLmmmin3=iLmmin3-kr(idims,3);iLmmmax1=iLmmax1-kr(idims,1)
iLmmmax2=iLmmax2-kr(idims,2);iLmmmax3=iLmmax3-kr(idims,3);
iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
iLpmin3=iLmin3+kr(idims,3);iLpmax1=iLmax1+kr(idims,1)
iLpmax2=iLmax2+kr(idims,2);iLpmax3=iLmax3+kr(idims,3);
iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
iLppmin3=iLpmin3+kr(idims,3);iLppmax1=iLpmax1+kr(idims,1)
iLppmax2=iLpmax2+kr(idims,2);iLppmax3=iLpmax3+kr(idims,3);

f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = 1.0d0/60.0d0 * (2.0d0* &
   w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,iLmmmin3:iLmmmax3) - 13.0d0* &
   w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iLmmin3:iLmmax3) + 47.0d0* &
   w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) + 27.0d0* w(iLpmin1:iLpmax1,&
   iLpmin2:iLpmax2,iLpmin3:iLpmax3) - 3.0d0*  w(iLppmin1:iLppmax1,&
   iLppmin2:iLppmax2,iLppmin3:iLppmax3))

! get fmp and ful:
a(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = w(iLpmin1:iLpmax1,&
   iLpmin2:iLpmax2,iLpmin3:iLpmax3)-w(iLmin1:iLmax1,iLmin2:iLmax2,&
   iLmin3:iLmax3)
b(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = alpha*(w(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3)-w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
   iLmmin3:iLmmax3))
call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,iLmin2,&
   iLmin3,iLmax1,iLmax2,iLmax3,a,b,tmp)
fmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = w(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3) + tmp(iLmin1:iLmax1,iLmin2:iLmax2,&
   iLmin3:iLmax3)
ful(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = w(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3) + b(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3)

! get dm4:
idmax1=iLmax1;idmax2=iLmax2;idmax3=iLmax3; idmin1=iLmin1-kr(idims,1)
idmin2=iLmin2-kr(idims,2);idmin3=iLmin3-kr(idims,3);
idmmin1=idmin1-kr(idims,1);idmmin2=idmin2-kr(idims,2)
idmmin3=idmin3-kr(idims,3);idmmax1=idmax1-kr(idims,1)
idmmax2=idmax2-kr(idims,2);idmmax3=idmax3-kr(idims,3);
idpmin1=idmin1+kr(idims,1);idpmin2=idmin2+kr(idims,2)
idpmin3=idmin3+kr(idims,3);idpmax1=idmax1+kr(idims,1)
idpmax2=idmax2+kr(idims,2);idpmax3=idmax3+kr(idims,3);

iemax1=idmax1+kr(idims,1);iemax2=idmax2+kr(idims,2);iemax3=idmax3+kr(idims,3)
iemin1=idmin1;iemin2=idmin2;iemin3=idmin3;
iemmin1=iemin1-kr(idims,1);iemmin2=iemin2-kr(idims,2)
iemmin3=iemin3-kr(idims,3);iemmax1=iemax1-kr(idims,1)
iemmax2=iemax2-kr(idims,2);iemmax3=iemax3-kr(idims,3);
iepmin1=iemin1+kr(idims,1);iepmin2=iemin2+kr(idims,2)
iepmin3=iemin3+kr(idims,3);iepmax1=iemax1+kr(idims,1)
iepmax2=iemax2+kr(idims,2);iepmax3=iemax3+kr(idims,3);

d(iemin1:iemax1,iemin2:iemax2,iemin3:iemax3) = w(iepmin1:iepmax1,&
   iepmin2:iepmax2,iepmin3:iepmax3)-2.0d0*w(iemin1:iemax1,iemin2:iemax2,&
   iemin3:iemax3)+w(iemmin1:iemmax1,iemmin2:iemmax2,iemmin3:iemmax3)

a(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = 4.0d0*d(idmin1:idmax1,&
   idmin2:idmax2,idmin3:idmax3)-d(idpmin1:idpmax1,idpmin2:idpmax2,&
   idpmin3:idpmax3)
b(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = 4.0d0*d(idpmin1:idpmax1,&
   idpmin2:idpmax2,idpmin3:idpmax3)-d(idmin1:idmax1,idmin2:idmax2,&
   idmin3:idmax3)
call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,idmin2,&
   idmin3,idmax1,idmax2,idmax3,a,b,tmp)
a(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = d(idmin1:idmax1,idmin2:idmax2,&
   idmin3:idmax3)
b(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = d(idpmin1:idpmax1,&
   idpmin2:idpmax2,idpmin3:idpmax3)
call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,idmin2,&
   idmin3,idmax1,idmax2,idmax3,a,b,tmp2)
call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,idmin2,&
   idmin3,idmax1,idmax2,idmax3,tmp,tmp2,tmp3)
dm4(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = tmp3(idmin1:idmax1,&
   idmin2:idmax2,idmin3:idmax3)

! get fmd:
fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = (w(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3)+w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
   iLpmin3:iLpmax3))/2.0d0-dm4(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3)&
   /2.0d0

!get flc: 
flc(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = half*(3.0d0*w(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3) - w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
   iLmmin3:iLmmax3)) + 4.0d0/3.0d0*dm4(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
   iLmmin3:iLmmax3)

fmin(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = max(min(w(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3),w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
   iLpmin3:iLpmax3),fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3)),&
   min(w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3),ful(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3),flc(iLmin1:iLmax1,iLmin2:iLmax2,&
   iLmin3:iLmax3)))

fmax(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = min(max(w(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3),w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
   iLpmin3:iLpmax3),fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3)),&
   max(w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3),ful(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3),flc(iLmin1:iLmax1,iLmin2:iLmax2,&
   iLmin3:iLmax3)))

call median(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,iLmin2,&
   iLmin3,iLmax1,iLmax2,iLmax3,fmin,f,fmax,tmp)
flim(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = tmp(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3)

! check case
where ((f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3)-w(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3))*(f(iLmin1:iLmax1,iLmin2:iLmax2,&
   iLmin3:iLmax3)-fmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3)) .le. eps)
   wLC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = f(iLmin1:iLmax1,&
      iLmin2:iLmax2,iLmin3:iLmax3)
elsewhere
   wLC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = flim(iLmin1:iLmax1,&
      iLmin2:iLmax2,iLmin3:iLmax3)
end where

! Right side:
! the interpolation from the right is obtained when the left-hand formula is applied to
! data mirrored about the interface.  
! thus substitute: 
! i-2 -> i+3
! i-1 -> i+2
! i   -> i+1
! i+1 -> i
! i+2 -> i-1

iLpppmin1=iLppmin1+kr(idims,1);iLpppmin2=iLppmin2+kr(idims,2)
iLpppmin3=iLppmin3+kr(idims,3);iLpppmax1=iLppmax1+kr(idims,1)
iLpppmax2=iLppmax2+kr(idims,2);iLpppmax3=iLppmax3+kr(idims,3);

f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = 1.0d0/60.0d0 * (2.0d0* &
   w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,iLpppmin3:iLpppmax3) - 13.0d0* &
   w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,iLppmin3:iLppmax3) + 47.0d0* &
   w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3) + 27.0d0* &
   w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) - 3.0d0*  w(iLmmin1:iLmmax1,&
   iLmmin2:iLmmax2,iLmmin3:iLmmax3))

! get fmp and ful:
a(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = w(iLmin1:iLmax1,iLmin2:iLmax2,&
   iLmin3:iLmax3)-w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3)
b(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = alpha*(w(iLpmin1:iLpmax1,&
   iLpmin2:iLpmax2,iLpmin3:iLpmax3)-w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
   iLppmin3:iLppmax3))
call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,iLmin2,&
   iLmin3,iLmax1,iLmax2,iLmax3,a,b,tmp)
fmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = w(iLpmin1:iLpmax1,&
   iLpmin2:iLpmax2,iLpmin3:iLpmax3) + tmp(iLmin1:iLmax1,iLmin2:iLmax2,&
   iLmin3:iLmax3)
ful(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = w(iLpmin1:iLpmax1,&
   iLpmin2:iLpmax2,iLpmin3:iLpmax3) + b(iLmin1:iLmax1,iLmin2:iLmax2,&
   iLmin3:iLmax3)

! get dm4:
idmax1=iLmax1+kr(idims,1);idmax2=iLmax2+kr(idims,2);idmax3=iLmax3+kr(idims,3)
idmin1=iLmin1;idmin2=iLmin2;idmin3=iLmin3;
idmmin1=idmin1-kr(idims,1);idmmin2=idmin2-kr(idims,2)
idmmin3=idmin3-kr(idims,3);idmmax1=idmax1-kr(idims,1)
idmmax2=idmax2-kr(idims,2);idmmax3=idmax3-kr(idims,3);
idpmin1=idmin1+kr(idims,1);idpmin2=idmin2+kr(idims,2)
idpmin3=idmin3+kr(idims,3);idpmax1=idmax1+kr(idims,1)
idpmax2=idmax2+kr(idims,2);idpmax3=idmax3+kr(idims,3);

iemax1=idmax1;iemax2=idmax2;iemax3=idmax3; iemin1=idmin1-kr(idims,1)
iemin2=idmin2-kr(idims,2);iemin3=idmin3-kr(idims,3);
iemmin1=iemin1-kr(idims,1);iemmin2=iemin2-kr(idims,2)
iemmin3=iemin3-kr(idims,3);iemmax1=iemax1-kr(idims,1)
iemmax2=iemax2-kr(idims,2);iemmax3=iemax3-kr(idims,3);
iepmin1=iemin1+kr(idims,1);iepmin2=iemin2+kr(idims,2)
iepmin3=iemin3+kr(idims,3);iepmax1=iemax1+kr(idims,1)
iepmax2=iemax2+kr(idims,2);iepmax3=iemax3+kr(idims,3);
ieppmin1=iepmin1+kr(idims,1);ieppmin2=iepmin2+kr(idims,2)
ieppmin3=iepmin3+kr(idims,3);ieppmax1=iepmax1+kr(idims,1)
ieppmax2=iepmax2+kr(idims,2);ieppmax3=iepmax3+kr(idims,3);

d(iemin1:iemax1,iemin2:iemax2,iemin3:iemax3) = w(iemin1:iemax1,iemin2:iemax2,&
   iemin3:iemax3)-2.0d0*w(iepmin1:iepmax1,iepmin2:iepmax2,iepmin3:iepmax3)+w&
   (ieppmin1:ieppmax1,ieppmin2:ieppmax2,ieppmin3:ieppmax3)

a(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = 4.0d0*d(idmin1:idmax1,&
   idmin2:idmax2,idmin3:idmax3)-d(idmmin1:idmmax1,idmmin2:idmmax2,&
   idmmin3:idmmax3)
b(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = 4.0d0*d(idmmin1:idmmax1,&
   idmmin2:idmmax2,idmmin3:idmmax3)-d(idmin1:idmax1,idmin2:idmax2,&
   idmin3:idmax3)
call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,idmin2,&
   idmin3,idmax1,idmax2,idmax3,a,b,tmp)
a(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = d(idmin1:idmax1,idmin2:idmax2,&
   idmin3:idmax3)
b(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = d(idmmin1:idmmax1,&
   idmmin2:idmmax2,idmmin3:idmmax3)
call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,idmin2,&
   idmin3,idmax1,idmax2,idmax3,a,b,tmp2)
call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,idmin2,&
   idmin3,idmax1,idmax2,idmax3,tmp,tmp2,tmp3)
dm4(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = tmp3(idmin1:idmax1,&
   idmin2:idmax2,idmin3:idmax3)

! get fmd:
fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = (w(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3)+w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
   iLpmin3:iLpmax3))/2.0d0-dm4(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3)&
   /2.0d0

!get flc: 
flc(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = half*(3.0d0*w&
   (iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3) - w(iLppmin1:iLppmax1,&
   iLppmin2:iLppmax2,iLppmin3:iLppmax3)) + 4.0d0/3.0d0*dm4(iLpmin1:iLpmax1,&
   iLpmin2:iLpmax2,iLpmin3:iLpmax3)

fmin(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = max(min(w(iLpmin1:iLpmax1,&
   iLpmin2:iLpmax2,iLpmin3:iLpmax3),w(iLmin1:iLmax1,iLmin2:iLmax2,&
   iLmin3:iLmax3),fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3)),&
   min(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3),ful(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3),flc(iLmin1:iLmax1,iLmin2:iLmax2,&
   iLmin3:iLmax3)))

fmax(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = min(max(w(iLpmin1:iLpmax1,&
   iLpmin2:iLpmax2,iLpmin3:iLpmax3),w(iLmin1:iLmax1,iLmin2:iLmax2,&
   iLmin3:iLmax3),fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3)),&
   max(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3),ful(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3),flc(iLmin1:iLmax1,iLmin2:iLmax2,&
   iLmin3:iLmax3)))

call median(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,iLmin2,&
   iLmin3,iLmax1,iLmax2,iLmax3,fmin,f,fmax,flim)

! check case
where ((f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3)-w(iLpmin1:iLpmax1,&
   iLpmin2:iLpmax2,iLpmin3:iLpmax3))*(f(iLmin1:iLmax1,iLmin2:iLmax2,&
   iLmin3:iLmax3)-fmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3))  .le. eps)
   wRC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = f(iLmin1:iLmax1,&
      iLmin2:iLmax2,iLmin3:iLmax3)
elsewhere
   wRC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = flim(iLmin1:iLmax1,&
      iLmin2:iLmax2,iLmin3:iLmax3)
end where

end subroutine MP5limitervar
!============================================================================
subroutine MP5limiterL(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,&
   iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,idims,w,wLC)
! MP5 limiter from Suresh & Huynh 1997
! Following the convention of Mignone et al. 2010.
! Needs at least three ghost cells.  Set dixB=3.

use mod_amrvacdef

integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3, iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3, idims
double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw)
double precision, intent(inout) :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw) 

! .. local ..
integer                         :: iLmmin1,iLmmin2,iLmmin3,iLmmax1,iLmmax2,&
   iLmmax3, iLmmmin1,iLmmmin2,iLmmmin3,iLmmmax1,iLmmmax2,iLmmmax3, iLpmin1,&
   iLpmin2,iLpmin3,iLpmax1,iLpmax2,iLpmax3, iLppmin1,iLppmin2,iLppmin3,&
   iLppmax1,iLppmax2,iLppmax3
integer                         :: idmin1,idmin2,idmin3,idmax1,idmax2,idmax3,&
    idpmin1,idpmin2,idpmin3,idpmax1,idpmax2,idpmax3, idppmin1,idppmin2,&
   idppmin3,idppmax1,idppmax2,idppmax3, idmmin1,idmmin2,idmmin3,idmmax1,&
   idmmax2,idmmax3, iemin1,iemin2,iemin3,iemax1,iemax2,iemax3, iemmin1,&
   iemmin2,iemmin3,iemmax1,iemmax2,iemmax3, iepmin1,iepmin2,iepmin3,iepmax1,&
   iepmax2,iepmax3, ieppmin1,ieppmin2,ieppmin3,ieppmax1,ieppmax2,ieppmax3
integer                         :: iw
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nw)  :: f, fmp, fmin, fmax, ful, dm4, d, fmd, flc, flim
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3) :: tmp, tmp2, tmp3, a, b, c
double precision, parameter     :: eps=0.0d0, alpha=4.0d0
!double precision                :: alpha
!----------------------------------------------------------------------------

! Variable alpha:
!alpha = float(nstep)/courantpar - one

! Left side:


iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
iLmmin3=iLmin3-kr(idims,3);iLmmax1=iLmax1-kr(idims,1)
iLmmax2=iLmax2-kr(idims,2);iLmmax3=iLmax3-kr(idims,3);
iLmmmin1=iLmmin1-kr(idims,1);iLmmmin2=iLmmin2-kr(idims,2)
iLmmmin3=iLmmin3-kr(idims,3);iLmmmax1=iLmmax1-kr(idims,1)
iLmmmax2=iLmmax2-kr(idims,2);iLmmmax3=iLmmax3-kr(idims,3);
iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
iLpmin3=iLmin3+kr(idims,3);iLpmax1=iLmax1+kr(idims,1)
iLpmax2=iLmax2+kr(idims,2);iLpmax3=iLmax3+kr(idims,3);
iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
iLppmin3=iLpmin3+kr(idims,3);iLppmax1=iLpmax1+kr(idims,1)
iLppmax2=iLpmax2+kr(idims,2);iLppmax3=iLpmax3+kr(idims,3);

f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = 1.0d0&
   /60.0d0 * (2.0d0* w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,iLmmmin3:iLmmmax3,&
   1:nwflux) - 13.0d0* w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iLmmin3:iLmmax3,&
   1:nwflux) + 47.0d0* w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
   1:nwflux) + 27.0d0* w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,&
   1:nwflux) - 3.0d0*  w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
   iLppmin3:iLppmax3,1:nwflux))

! get fmp and ful:
do iw=1,nwflux
   a(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = w(iLpmin1:iLpmax1,&
      iLpmin2:iLpmax2,iLpmin3:iLpmax3,iw)-w(iLmin1:iLmax1,iLmin2:iLmax2,&
      iLmin3:iLmax3,iw)
   b(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = alpha*(w(iLmin1:iLmax1,&
      iLmin2:iLmax2,iLmin3:iLmax3,iw)-w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
      iLmmin3:iLmmax3,iw))
   call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,iLmin2,&
      iLmin3,iLmax1,iLmax2,iLmax3,a,b,tmp)
   fmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,iw) = w(iLmin1:iLmax1,&
      iLmin2:iLmax2,iLmin3:iLmax3,iw) + tmp(iLmin1:iLmax1,iLmin2:iLmax2,&
      iLmin3:iLmax3)
   ful(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,iw) = w(iLmin1:iLmax1,&
      iLmin2:iLmax2,iLmin3:iLmax3,iw) + b(iLmin1:iLmax1,iLmin2:iLmax2,&
      iLmin3:iLmax3)
end do ! iw loop

! get dm4:
idmax1=iLmax1;idmax2=iLmax2;idmax3=iLmax3; idmin1=iLmin1-kr(idims,1)
idmin2=iLmin2-kr(idims,2);idmin3=iLmin3-kr(idims,3);
idmmin1=idmin1-kr(idims,1);idmmin2=idmin2-kr(idims,2)
idmmin3=idmin3-kr(idims,3);idmmax1=idmax1-kr(idims,1)
idmmax2=idmax2-kr(idims,2);idmmax3=idmax3-kr(idims,3);
idpmin1=idmin1+kr(idims,1);idpmin2=idmin2+kr(idims,2)
idpmin3=idmin3+kr(idims,3);idpmax1=idmax1+kr(idims,1)
idpmax2=idmax2+kr(idims,2);idpmax3=idmax3+kr(idims,3);

iemax1=idmax1+kr(idims,1);iemax2=idmax2+kr(idims,2);iemax3=idmax3+kr(idims,3)
iemin1=idmin1;iemin2=idmin2;iemin3=idmin3;
iemmin1=iemin1-kr(idims,1);iemmin2=iemin2-kr(idims,2)
iemmin3=iemin3-kr(idims,3);iemmax1=iemax1-kr(idims,1)
iemmax2=iemax2-kr(idims,2);iemmax3=iemax3-kr(idims,3);
iepmin1=iemin1+kr(idims,1);iepmin2=iemin2+kr(idims,2)
iepmin3=iemin3+kr(idims,3);iepmax1=iemax1+kr(idims,1)
iepmax2=iemax2+kr(idims,2);iepmax3=iemax3+kr(idims,3);

d(iemin1:iemax1,iemin2:iemax2,iemin3:iemax3,1:nwflux) = w(iepmin1:iepmax1,&
   iepmin2:iepmax2,iepmin3:iepmax3,1:nwflux)-2.0d0*w(iemin1:iemax1,&
   iemin2:iemax2,iemin3:iemax3,1:nwflux)+w(iemmin1:iemmax1,iemmin2:iemmax2,&
   iemmin3:iemmax3,1:nwflux)

do iw=1,nwflux
   a(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = 4.0d0*d(idmin1:idmax1,&
      idmin2:idmax2,idmin3:idmax3,iw)-d(idpmin1:idpmax1,idpmin2:idpmax2,&
      idpmin3:idpmax3,iw)
   b(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = 4.0d0*d(idpmin1:idpmax1,&
      idpmin2:idpmax2,idpmin3:idpmax3,iw)-d(idmin1:idmax1,idmin2:idmax2,&
      idmin3:idmax3,iw)
   call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,idmin2,&
      idmin3,idmax1,idmax2,idmax3,a,b,tmp)
   a(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = d(idmin1:idmax1,&
      idmin2:idmax2,idmin3:idmax3,iw)
   b(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = d(idpmin1:idpmax1,&
      idpmin2:idpmax2,idpmin3:idpmax3,iw)
   call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,idmin2,&
      idmin3,idmax1,idmax2,idmax3,a,b,tmp2)
   call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,idmin2,&
      idmin3,idmax1,idmax2,idmax3,tmp,tmp2,tmp3)
   dm4(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3,iw) = tmp3(idmin1:idmax1,&
      idmin2:idmax2,idmin3:idmax3)
end do

! get fmd:
fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = (w(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)+w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
   iLpmin3:iLpmax3,1:nwflux))/2.0d0-dm4(iLmin1:iLmax1,iLmin2:iLmax2,&
   iLmin3:iLmax3,1:nwflux)/2.0d0

!get flc: 
flc(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = half*(3.0d0*w&
   (iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) - w(iLmmin1:iLmmax1,&
   iLmmin2:iLmmax2,iLmmin3:iLmmax3,1:nwflux)) + 4.0d0/3.0d0*dm4&
   (iLmmin1:iLmmax1,iLmmin2:iLmmax2,iLmmin3:iLmmax3,1:nwflux)

fmin(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = &
   max(min(w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),&
   w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux),&
   fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)),&
   min(w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),&
   ful(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),flc(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)))

fmax(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = &
   min(max(w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),&
   w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux),&
   fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)),&
   max(w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),&
   ful(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),flc(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)))

do iw=1,nwflux
   a(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = fmin(iLmin1:iLmax1,&
      iLmin2:iLmax2,iLmin3:iLmax3,iw)
   b(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = f(iLmin1:iLmax1,&
      iLmin2:iLmax2,iLmin3:iLmax3,iw)
   c(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = fmax(iLmin1:iLmax1,&
      iLmin2:iLmax2,iLmin3:iLmax3,iw)
   call median(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,iLmin2,&
      iLmin3,iLmax1,iLmax2,iLmax3,a,b,c,tmp)
   flim(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,iw) = tmp(iLmin1:iLmax1,&
      iLmin2:iLmax2,iLmin3:iLmax3)
end do


! check case
where ((f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)-w(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux))*(f(iLmin1:iLmax1,iLmin2:iLmax2,&
   iLmin3:iLmax3,1:nwflux)-fmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
   1:nwflux)) .le. eps)
   wLC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) &
      = f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)
elsewhere
   wLC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) &
      = flim(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)
end where


end subroutine MP5limiterL
!============================================================================
subroutine MP5limiterR(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,&
   iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,idims,w,wRC)
! MP5 limiter from Suresh & Huynh 1997
! Following the convention of Mignone et al. 2010.
! Needs at least three ghost cells.  Set dixB=3.

use mod_amrvacdef

integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3, iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3, idims
double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw)
double precision, intent(inout) :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw)

! .. local ..
integer                         :: iLmmin1,iLmmin2,iLmmin3,iLmmax1,iLmmax2,&
   iLmmax3, iLpmin1,iLpmin2,iLpmin3,iLpmax1,iLpmax2,iLpmax3, iLppmin1,&
   iLppmin2,iLppmin3,iLppmax1,iLppmax2,iLppmax3, iLpppmin1,iLpppmin2,&
   iLpppmin3,iLpppmax1,iLpppmax2,iLpppmax3
integer                         :: idmin1,idmin2,idmin3,idmax1,idmax2,idmax3,&
    idpmin1,idpmin2,idpmin3,idpmax1,idpmax2,idpmax3, idppmin1,idppmin2,&
   idppmin3,idppmax1,idppmax2,idppmax3, idmmin1,idmmin2,idmmin3,idmmax1,&
   idmmax2,idmmax3, iemin1,iemin2,iemin3,iemax1,iemax2,iemax3, iemmin1,&
   iemmin2,iemmin3,iemmax1,iemmax2,iemmax3, iepmin1,iepmin2,iepmin3,iepmax1,&
   iepmax2,iepmax3, ieppmin1,ieppmin2,ieppmin3,ieppmax1,ieppmax2,ieppmax3
integer                         :: iw
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nw)  :: f, fmp, fmin, fmax, ful, dm4, d, fmd, flc, flim
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3) :: tmp, tmp2, tmp3, a, b, c
double precision, parameter     :: eps=0.0d0, alpha=4.0d0
!double precision                :: alpha
!----------------------------------------------------------------------------
! Right side:
! the interpolation from the right is obtained when the left-hand formula is applied to
! data mirrored about the interface.  
! thus substitute: 
! i-2 -> i+3
! i-1 -> i+2
! i   -> i+1
! i+1 -> i
! i+2 -> i-1

iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
iLmmin3=iLmin3-kr(idims,3);iLmmax1=iLmax1-kr(idims,1)
iLmmax2=iLmax2-kr(idims,2);iLmmax3=iLmax3-kr(idims,3);
iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
iLpmin3=iLmin3+kr(idims,3);iLpmax1=iLmax1+kr(idims,1)
iLpmax2=iLmax2+kr(idims,2);iLpmax3=iLmax3+kr(idims,3);
iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
iLppmin3=iLpmin3+kr(idims,3);iLppmax1=iLpmax1+kr(idims,1)
iLppmax2=iLpmax2+kr(idims,2);iLppmax3=iLpmax3+kr(idims,3);
iLpppmin1=iLppmin1+kr(idims,1);iLpppmin2=iLppmin2+kr(idims,2)
iLpppmin3=iLppmin3+kr(idims,3);iLpppmax1=iLppmax1+kr(idims,1)
iLpppmax2=iLppmax2+kr(idims,2);iLpppmax3=iLppmax3+kr(idims,3);

f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = 1.0d0&
   /60.0d0 * (2.0d0* w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
   iLpppmin3:iLpppmax3,1:nwflux) - 13.0d0* w(iLppmin1:iLppmax1,&
   iLppmin2:iLppmax2,iLppmin3:iLppmax3,1:nwflux) + 47.0d0* w(iLpmin1:iLpmax1,&
   iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux) + 27.0d0* w(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) - 3.0d0*  w(iLmmin1:iLmmax1,&
   iLmmin2:iLmmax2,iLmmin3:iLmmax3,1:nwflux))

! get fmp and ful:
do iw=1,nwflux
   a(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = w(iLmin1:iLmax1,&
      iLmin2:iLmax2,iLmin3:iLmax3,iw)-w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
      iLpmin3:iLpmax3,iw)
   b(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = alpha*(w(iLpmin1:iLpmax1,&
      iLpmin2:iLpmax2,iLpmin3:iLpmax3,iw)-w(iLppmin1:iLppmax1,&
      iLppmin2:iLppmax2,iLppmin3:iLppmax3,iw))
   call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,iLmin2,&
      iLmin3,iLmax1,iLmax2,iLmax3,a,b,tmp)
   fmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,iw) = w(iLpmin1:iLpmax1,&
      iLpmin2:iLpmax2,iLpmin3:iLpmax3,iw) + tmp(iLmin1:iLmax1,iLmin2:iLmax2,&
      iLmin3:iLmax3)
   ful(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,iw) = w(iLpmin1:iLpmax1,&
      iLpmin2:iLpmax2,iLpmin3:iLpmax3,iw) + b(iLmin1:iLmax1,iLmin2:iLmax2,&
      iLmin3:iLmax3)
end do ! iw loop

! get dm4:
idmax1=iLmax1+kr(idims,1);idmax2=iLmax2+kr(idims,2);idmax3=iLmax3+kr(idims,3)
idmin1=iLmin1;idmin2=iLmin2;idmin3=iLmin3;
idmmin1=idmin1-kr(idims,1);idmmin2=idmin2-kr(idims,2)
idmmin3=idmin3-kr(idims,3);idmmax1=idmax1-kr(idims,1)
idmmax2=idmax2-kr(idims,2);idmmax3=idmax3-kr(idims,3);
idpmin1=idmin1+kr(idims,1);idpmin2=idmin2+kr(idims,2)
idpmin3=idmin3+kr(idims,3);idpmax1=idmax1+kr(idims,1)
idpmax2=idmax2+kr(idims,2);idpmax3=idmax3+kr(idims,3);

iemax1=idmax1;iemax2=idmax2;iemax3=idmax3; iemin1=idmin1-kr(idims,1)
iemin2=idmin2-kr(idims,2);iemin3=idmin3-kr(idims,3);
iemmin1=iemin1-kr(idims,1);iemmin2=iemin2-kr(idims,2)
iemmin3=iemin3-kr(idims,3);iemmax1=iemax1-kr(idims,1)
iemmax2=iemax2-kr(idims,2);iemmax3=iemax3-kr(idims,3);
iepmin1=iemin1+kr(idims,1);iepmin2=iemin2+kr(idims,2)
iepmin3=iemin3+kr(idims,3);iepmax1=iemax1+kr(idims,1)
iepmax2=iemax2+kr(idims,2);iepmax3=iemax3+kr(idims,3);
ieppmin1=iepmin1+kr(idims,1);ieppmin2=iepmin2+kr(idims,2)
ieppmin3=iepmin3+kr(idims,3);ieppmax1=iepmax1+kr(idims,1)
ieppmax2=iepmax2+kr(idims,2);ieppmax3=iepmax3+kr(idims,3);

d(iemin1:iemax1,iemin2:iemax2,iemin3:iemax3,1:nwflux) = w(iemin1:iemax1,&
   iemin2:iemax2,iemin3:iemax3,1:nwflux)-2.0d0*w(iepmin1:iepmax1,&
   iepmin2:iepmax2,iepmin3:iepmax3,1:nwflux)+w(ieppmin1:ieppmax1,&
   ieppmin2:ieppmax2,ieppmin3:ieppmax3,1:nwflux)

do iw=1,nwflux
   a(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = 4.0d0*d(idmin1:idmax1,&
      idmin2:idmax2,idmin3:idmax3,iw)-d(idmmin1:idmmax1,idmmin2:idmmax2,&
      idmmin3:idmmax3,iw)
   b(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = 4.0d0*d(idmmin1:idmmax1,&
      idmmin2:idmmax2,idmmin3:idmmax3,iw)-d(idmin1:idmax1,idmin2:idmax2,&
      idmin3:idmax3,iw)
   call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,idmin2,&
      idmin3,idmax1,idmax2,idmax3,a,b,tmp)
   a(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = d(idmin1:idmax1,&
      idmin2:idmax2,idmin3:idmax3,iw)
   b(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = d(idmmin1:idmmax1,&
      idmmin2:idmmax2,idmmin3:idmmax3,iw)
   call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,idmin2,&
      idmin3,idmax1,idmax2,idmax3,a,b,tmp2)
   call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,idmin2,&
      idmin3,idmax1,idmax2,idmax3,tmp,tmp2,tmp3)
   dm4(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3,iw) = tmp3(idmin1:idmax1,&
      idmin2:idmax2,idmin3:idmax3)
end do

! get fmd:
fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = (w(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)+w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
   iLpmin3:iLpmax3,1:nwflux))/2.0d0-dm4(iLmin1:iLmax1,iLmin2:iLmax2,&
   iLmin3:iLmax3,1:nwflux)/2.0d0

!get flc: 
flc(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = half*(3.0d0*w&
   (iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux) - &
   w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,iLppmin3:iLppmax3,1:nwflux)) + 4.0d0&
   /3.0d0*dm4(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux)

fmin(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = &
   max(min(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux),&
   w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),fmd(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)),min(w(iLpmin1:iLpmax1,&
   iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux),ful(iLmin1:iLmax1,iLmin2:iLmax2,&
   iLmin3:iLmax3,1:nwflux),flc(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
   1:nwflux)))

fmax(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = &
   min(max(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux),&
   w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),fmd(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)),max(w(iLpmin1:iLpmax1,&
   iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux),ful(iLmin1:iLmax1,iLmin2:iLmax2,&
   iLmin3:iLmax3,1:nwflux),flc(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
   1:nwflux)))

do iw=1,nwflux
   a(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = fmin(iLmin1:iLmax1,&
      iLmin2:iLmax2,iLmin3:iLmax3,iw)
   b(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = f(iLmin1:iLmax1,&
      iLmin2:iLmax2,iLmin3:iLmax3,iw)
   c(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = fmax(iLmin1:iLmax1,&
      iLmin2:iLmax2,iLmin3:iLmax3,iw)
   call median(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,iLmin2,&
      iLmin3,iLmax1,iLmax2,iLmax3,a,b,c,tmp)
   flim(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,iw) = tmp(iLmin1:iLmax1,&
      iLmin2:iLmax2,iLmin3:iLmax3)
end do

! check case
where ((f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)-w&
   (iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux))*(f&
   (iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)-fmp(iLmin1:iLmax1,&
   iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux))  .le. eps)
   wRC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) &
      = f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)
elsewhere
   wRC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) &
      = flim(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)
end where

end subroutine Mp5limiterR
!============================================================================
subroutine minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
   ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,a,b,minm)

use mod_amrvacdef

integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(in) :: a(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3), b(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
double precision, intent(out):: minm(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)
!-----------------------------------------------------------------------------

minm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = (sign(one,&
   a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))+sign(one,&
   b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)))/2.0d0 * &
   min(abs(a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)),&
   abs(b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)))

end subroutine minmod
!============================================================================
subroutine median(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
   ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,a,b,c,med)

use mod_amrvacdef

integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(in) :: a(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3), b(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
    c(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
double precision, intent(out):: med(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)
! .. local ..
double precision             :: tmp1(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3),tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
!-----------------------------------------------------------------------------

tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = b(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,ixOmin3:ixOmax3) - a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   ixOmin3:ixOmax3)
tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = c(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,ixOmin3:ixOmax3) - a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   ixOmin3:ixOmax3)

med(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = a(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,ixOmin3:ixOmax3) + (sign(one,tmp1(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,ixOmin3:ixOmax3))+sign(one,tmp2(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,ixOmin3:ixOmax3)))/2.0d0 * min(abs(tmp1(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,ixOmin3:ixOmax3)),abs(tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   ixOmin3:ixOmax3)))

end subroutine median
!============================================================================
  
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
subroutine WENO5limiter(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,idims,w,wLC,wRC)
  ! WENO5 Limiter, See Del Zanna et al. 2007 for a small overview
  ! Needs at least three ghost cells.  Set dixB=3.

  use mod_amrvacdef

  integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3, idims
  double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:nw)
  double precision, intent(inout) :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:nw),wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:nw) 
  ! .. local ..
  integer                         :: iLmmin1,iLmmin2,iLmmin3,iLmmax1,iLmmax2,&
     iLmmax3, iLmmmin1,iLmmmin2,iLmmmin3,iLmmmax1,iLmmmax2,iLmmmax3, iLpmin1,&
     iLpmin2,iLpmin3,iLpmax1,iLpmax2,iLpmax3, iLppmin1,iLppmin2,iLppmin3,&
     iLppmax1,iLppmax2,iLppmax3, iLpppmin1,iLpppmin2,iLpppmin3,iLpppmax1,&
     iLpppmax2,iLpppmax3
  double precision                :: f_array(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:nw,3), d_array(3), beta(ixImin1:ixImax1,&
     ixImin2:ixImax2,ixImin3:ixImax3,1:nw,3), beta_coeff(2),&
      alpha_array(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw,3),&
      tau_5(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw),&
      alpha_sum(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw),&
      usum(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw),&
      flux(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw)
  integer                         :: i, iw
  double precision, parameter     :: weno_eps_machine = 1.0d-42
  !----------------------------------------------------------------------------



  ! iL^L holds the indices of interfaces to reconstruct to.  Convention is that a center index holds the _right-side_ interface.  

  iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
  iLmmin3=iLmin3-kr(idims,3);iLmmax1=iLmax1-kr(idims,1)
  iLmmax2=iLmax2-kr(idims,2);iLmmax3=iLmax3-kr(idims,3);
  iLmmmin1=iLmmin1-kr(idims,1);iLmmmin2=iLmmin2-kr(idims,2)
  iLmmmin3=iLmmin3-kr(idims,3);iLmmmax1=iLmmax1-kr(idims,1)
  iLmmmax2=iLmmax2-kr(idims,2);iLmmmax3=iLmmax3-kr(idims,3);
  iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
  iLpmin3=iLmin3+kr(idims,3);iLpmax1=iLmax1+kr(idims,1)
  iLpmax2=iLmax2+kr(idims,2);iLpmax3=iLmax3+kr(idims,3);
  iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
  iLppmin3=iLpmin3+kr(idims,3);iLppmax1=iLpmax1+kr(idims,1)
  iLppmax2=iLpmax2+kr(idims,2);iLppmax3=iLpmax3+kr(idims,3);



  f_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,1) = 3.0d0&
     /8.0d0 * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,iLmmmin3:iLmmmax3,&
     1:nwflux) - 10.0d0/8.0d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
     iLmmin3:iLmmax3,1:nwflux) + 15.0d0/8.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
     iLmin3:iLmax3,1:nwflux) !first stencil
  f_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,2) = -1.0d0&
     /8.0d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iLmmin3:iLmmax3,&
     1:nwflux) + 6.0d0/8.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
     1:nwflux) + 3.0d0/8.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
     iLpmin3:iLpmax3,1:nwflux)
  f_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,3) = 3.0d0&
     /8.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) + 6.0d0&
     /8.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,&
     1:nwflux) - 1.0d0/8.0d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
     iLppmin3:iLppmax3,1:nwflux)  

  d_array(1:3) = (/ 1.0d0/16.0d0, 10.0d0/16.0d0, 5.0d0/16.0d0 /)


  beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)

  beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,1) &
     = beta_coeff(1) * (w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
     iLmmmin3:iLmmmax3,1:nwflux) + w(iLmin1:iLmax1,iLmin2:iLmax2,&
     iLmin3:iLmax3,1:nwflux) - 2.0d0*w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
     iLmmin3:iLmmax3,1:nwflux))**2 + beta_coeff(2) * (w(iLmmmin1:iLmmmax1,&
     iLmmmin2:iLmmmax2,iLmmmin3:iLmmmax3,1:nwflux) - 4.0d0 * &
     w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iLmmin3:iLmmax3,1:nwflux) + &
     3.0d0*w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux))**2

  beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,2) &
     = beta_coeff(1) * (w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iLmmin3:iLmmax3,&
     1:nwflux) + w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,&
     1:nwflux) - 2.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
     1:nwflux))**2 + beta_coeff(2) * (w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
     iLmmin3:iLmmax3,1:nwflux) - w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
     iLpmin3:iLpmax3,1:nwflux))**2

  beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,3) &
     = beta_coeff(1) * (w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
     1:nwflux) + w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,iLppmin3:iLppmax3,&
     1:nwflux) - 2.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,&
     1:nwflux))**2 + beta_coeff(2) * (3.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
     iLmin3:iLmax3, 1:nwflux) - 4.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
     iLpmin3:iLpmax3,1:nwflux) + w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
     iLppmin3:iLppmax3,1:nwflux))**2


  alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = 0.0d0 
  do i = 1,3
     alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,i) &
        = d_array(i)/(beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,&
        i) + weno_eps_machine)**2
     alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) &
        = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
        1:nwflux) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
        1:nwflux,i)
  end do


  flux(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = 0.0d0
  do i = 1,3
     flux(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) &
        = flux(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) + &
        f_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,&
        i) * alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,&
        i)/(alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux))
  end do

  !left value at right interface
  wLC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = &
     flux(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)

  !now other side

  iLpppmin1=iLppmin1+kr(idims,1);iLpppmin2=iLppmin2+kr(idims,2)
  iLpppmin3=iLppmin3+kr(idims,3);iLpppmax1=iLppmax1+kr(idims,1)
  iLpppmax2=iLppmax2+kr(idims,2);iLpppmax3=iLppmax3+kr(idims,3);

  f_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,1) = 3.0d0&
     /8.0d0 * w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,iLpppmin3:iLpppmax3,&
     1:nwflux) - 10.0d0/8.0d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
     iLppmin3:iLppmax3,1:nwflux) + 15.0d0/8.0d0 * w(iLpmin1:iLpmax1,&
     iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux) !first stencil
  f_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,2) = -1.0d0&
     /8.0d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,iLppmin3:iLppmax3,&
     1:nwflux) + 6.0d0/8.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
     iLpmin3:iLpmax3,1:nwflux) + 3.0d0/8.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
     iLmin3:iLmax3,1:nwflux)
  f_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,3) = 3.0d0&
     /8.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,&
     1:nwflux) + 6.0d0/8.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
     1:nwflux) - 1.0d0/8.0d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
     iLmmin3:iLmmax3,1:nwflux)  

  beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,1) &
     = beta_coeff(1) * (w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
     iLpppmin3:iLpppmax3,1:nwflux) + w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
     iLpmin3:iLpmax3,1:nwflux) - 2.0d0*w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
     iLppmin3:iLppmax3,1:nwflux))**2 + beta_coeff(2) * (w(iLpppmin1:iLpppmax1,&
     iLpppmin2:iLpppmax2,iLpppmin3:iLpppmax3,1:nwflux) - 4.0d0 * &
     w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,iLppmin3:iLppmax3,&
     1:nwflux) + 3.0d0*w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,&
     1:nwflux))**2

  beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,2) &
     = beta_coeff(1) * (w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
     iLppmin3:iLppmax3,1:nwflux) + w(iLmin1:iLmax1,iLmin2:iLmax2,&
     iLmin3:iLmax3,1:nwflux) - 2.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
     iLpmin3:iLpmax3,1:nwflux))**2 + beta_coeff(2) * (w(iLppmin1:iLppmax1,&
     iLppmin2:iLppmax2,iLppmin3:iLppmax3,1:nwflux) - w(iLmin1:iLmax1,&
     iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux))**2

  beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,3) &
     = beta_coeff(1) * (w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,&
     1:nwflux) + w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iLmmin3:iLmmax3,&
     1:nwflux) - 2.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
     1:nwflux))**2 + beta_coeff(2) * (3.0d0 * w(iLpmin1:iLpmax1,&
     iLpmin2:iLpmax2,iLpmin3:iLpmax3, 1:nwflux) - 4.0d0 * w(iLmin1:iLmax1,&
     iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) + w(iLmmin1:iLmmax1,&
     iLmmin2:iLmmax2,iLmmin3:iLmmax3,1:nwflux))**2

  alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = 0.0d0 
  do i = 1,3
     alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,i) &
        = d_array(i)/(beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,&
        i) + weno_eps_machine)**2
     alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) &
        = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
        1:nwflux) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
        1:nwflux,i)
  end do

  flux(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = 0.0d0
  do i = 1,3
     flux(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) &
        = flux(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) + &
        f_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,&
        i) * alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,&
        i)/(alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux))
  end do

  !right value at right interface
  wRC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = &
     flux(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)

end subroutine WENO5limiter
!=============================================================================
subroutine WENO5limitervar(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,idims,w,wLC,wRC)

  use mod_amrvacdef

  integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3, idims
  double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3)

  double precision, intent(inout) :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3),wLC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) 
  ! .. local ..
  integer                         :: iLmmin1,iLmmin2,iLmmin3,iLmmax1,iLmmax2,&
     iLmmax3, iLmmmin1,iLmmmin2,iLmmmin3,iLmmmax1,iLmmmax2,iLmmmax3, iLpmin1,&
     iLpmin2,iLpmin3,iLpmax1,iLpmax2,iLpmax3, iLppmin1,iLppmin2,iLppmin3,&
     iLppmax1,iLppmax2,iLppmax3, iLpppmin1,iLpppmin2,iLpppmin3,iLpppmax1,&
     iLpppmax2,iLpppmax3
  double precision                :: f_array(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,3), d_array(3), beta(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,3), beta_coeff(2), alpha_array(ixImin1:ixImax1,&
     ixImin2:ixImax2,ixImin3:ixImax3,3), tau_5(ixImin1:ixImax1,&
     ixImin2:ixImax2,ixImin3:ixImax3), alpha_sum(ixImin1:ixImax1,&
     ixImin2:ixImax2,ixImin3:ixImax3), usum(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3), flux(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
  integer                         :: i
  double precision, parameter     :: weno_eps_machine = 1.0d-42
  !----------------------------------------------------------------------------



  ! iL^L holds the indices of interfaces to reconstruct to.  Convention is that a center index holds the _right-side_ interface.  

  iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
  iLmmin3=iLmin3-kr(idims,3);iLmmax1=iLmax1-kr(idims,1)
  iLmmax2=iLmax2-kr(idims,2);iLmmax3=iLmax3-kr(idims,3);
  iLmmmin1=iLmmin1-kr(idims,1);iLmmmin2=iLmmin2-kr(idims,2)
  iLmmmin3=iLmmin3-kr(idims,3);iLmmmax1=iLmmax1-kr(idims,1)
  iLmmmax2=iLmmax2-kr(idims,2);iLmmmax3=iLmmax3-kr(idims,3);
  iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
  iLpmin3=iLmin3+kr(idims,3);iLpmax1=iLmax1+kr(idims,1)
  iLpmax2=iLmax2+kr(idims,2);iLpmax3=iLmax3+kr(idims,3);
  iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
  iLppmin3=iLpmin3+kr(idims,3);iLppmax1=iLpmax1+kr(idims,1)
  iLppmax2=iLpmax2+kr(idims,2);iLppmax3=iLpmax3+kr(idims,3);



  f_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1) = 3.0d0&
     /8.0d0 * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,iLmmmin3:iLmmmax3) - &
     10.0d0/8.0d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iLmmin3:iLmmax3) + &
     15.0d0/8.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) !first stencil
  f_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,2) = -1.0d0&
     /8.0d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iLmmin3:iLmmax3) + 6.0d0&
     /8.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) + 3.0d0&
     /8.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3)
  f_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,3) = 3.0d0&
     /8.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) + 6.0d0&
     /8.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3) - 1.0d0&
     /8.0d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,iLppmin3:iLppmax3)  

  d_array(1:3) = (/ 1.0d0/16.0d0, 10.0d0/16.0d0, 5.0d0/16.0d0 /)


  beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)

  beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1) = beta_coeff(1) * &
     (w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,iLmmmin3:iLmmmax3) + &
     w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) - 2.0d0*w(iLmmin1:iLmmax1,&
     iLmmin2:iLmmax2,iLmmin3:iLmmax3))**2 + beta_coeff(2) * &
     (w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,iLmmmin3:iLmmmax3) - 4.0d0 * &
     w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iLmmin3:iLmmax3) + 3.0d0*w&
     (iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3))**2

  beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,2) = beta_coeff(1) * &
     (w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iLmmin3:iLmmax3) + w(iLpmin1:iLpmax1,&
     iLpmin2:iLpmax2,iLpmin3:iLpmax3) - 2.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
     iLmin3:iLmax3))**2 + beta_coeff(2) * (w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
     iLmmin3:iLmmax3) - w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3))**2

  beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,3) = beta_coeff(1) * &
     (w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) + w(iLppmin1:iLppmax1,&
     iLppmin2:iLppmax2,iLppmin3:iLppmax3) - 2.0d0 * w(iLpmin1:iLpmax1,&
     iLpmin2:iLpmax2,iLpmin3:iLpmax3))**2 + beta_coeff(2) * (3.0d0 * &
     w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) - 4.0d0 * w(iLpmin1:iLpmax1,&
     iLpmin2:iLpmax2,iLpmin3:iLpmax3) + w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
     iLppmin3:iLppmax3))**2


  alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = 0.0d0 
  do i = 1,3
     alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,i) = d_array(i)&
        /(beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,i) + &
        weno_eps_machine)**2
     alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = &
        alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) + &
        alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,i)
  end do


  flux(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = 0.0d0
  do i = 1,3
     flux(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = flux(iLmin1:iLmax1,&
        iLmin2:iLmax2,iLmin3:iLmax3) + f_array(iLmin1:iLmax1,iLmin2:iLmax2,&
        iLmin3:iLmax3,i) * alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,&
        iLmin3:iLmax3,i)/(alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
        iLmin3:iLmax3))
  end do

  !left value at right interface
  wLC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = flux(iLmin1:iLmax1,&
     iLmin2:iLmax2,iLmin3:iLmax3)

  !now other side

  iLpppmin1=iLppmin1+kr(idims,1);iLpppmin2=iLppmin2+kr(idims,2)
  iLpppmin3=iLppmin3+kr(idims,3);iLpppmax1=iLppmax1+kr(idims,1)
  iLpppmax2=iLppmax2+kr(idims,2);iLpppmax3=iLppmax3+kr(idims,3);

  f_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1) = 3.0d0&
     /8.0d0 * w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,iLpppmin3:iLpppmax3) &
     - 10.0d0/8.0d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
     iLppmin3:iLppmax3) + 15.0d0/8.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
     iLpmin3:iLpmax3) !first stencil
  f_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,2) = -1.0d0&
     /8.0d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,iLppmin3:iLppmax3) + &
     6.0d0/8.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3) + 3.0d0&
     /8.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3)
  f_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,3) = 3.0d0&
     /8.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3) + 6.0d0&
     /8.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) - 1.0d0&
     /8.0d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iLmmin3:iLmmax3)  

  beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1) = beta_coeff(1) * &
     (w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,iLpppmin3:iLpppmax3) + &
     w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3) - 2.0d0*w&
     (iLppmin1:iLppmax1,iLppmin2:iLppmax2,iLppmin3:iLppmax3))&
     **2 + beta_coeff(2) * (w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
     iLpppmin3:iLpppmax3) - 4.0d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
     iLppmin3:iLppmax3) + 3.0d0*w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
     iLpmin3:iLpmax3))**2

  beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,2) = beta_coeff(1) * &
     (w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,iLppmin3:iLppmax3) + &
     w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) - 2.0d0 * w(iLpmin1:iLpmax1,&
     iLpmin2:iLpmax2,iLpmin3:iLpmax3))**2 + beta_coeff(2) * &
     (w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,iLppmin3:iLppmax3) - &
     w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3))**2

  beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,3) = beta_coeff(1) * &
     (w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3) + w(iLmmin1:iLmmax1,&
     iLmmin2:iLmmax2,iLmmin3:iLmmax3) - 2.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
     iLmin3:iLmax3))**2 + beta_coeff(2) * (3.0d0 * w(iLpmin1:iLpmax1,&
     iLpmin2:iLpmax2,iLpmin3:iLpmax3) - 4.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
     iLmin3:iLmax3) + w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iLmmin3:iLmmax3))**2

  alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = 0.0d0 
  do i = 1,3
     alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,i) = d_array(i)&
        /(beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,i) + &
        weno_eps_machine)**2
     alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = &
        alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) + &
        alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,i)
  end do

  flux(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = 0.0d0
  do i = 1,3
     flux(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = flux(iLmin1:iLmax1,&
        iLmin2:iLmax2,iLmin3:iLmax3) + f_array(iLmin1:iLmax1,iLmin2:iLmax2,&
        iLmin3:iLmax3,i) * alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,&
        iLmin3:iLmax3,i)/(alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
        iLmin3:iLmax3))
  end do

  !right value at right interface
  wRC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = flux(iLmin1:iLmax1,&
     iLmin2:iLmax2,iLmin3:iLmax3)

end subroutine WENO5limitervar
!=============================================================================
  
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
subroutine WENOZPlimiter(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,idims,dxdim,w,wLC,wRC)
  ! WENOZ and WENOZ+ Limiter, see Acker, Borges and Costa (2016)
  ! Needs at least three ghost cells.  Set dixB=3.

  use mod_amrvacdef

  integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3, idims
  double precision, intent(in)    :: dxdim
  double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:nw)
  double precision, intent(inout) :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:nw),wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:nw) 
  ! .. local ..
  integer                         :: iLmmin1,iLmmin2,iLmmin3,iLmmax1,iLmmax2,&
     iLmmax3, iLmmmin1,iLmmmin2,iLmmmin3,iLmmmax1,iLmmmax2,iLmmmax3, iLpmin1,&
     iLpmin2,iLpmin3,iLpmax1,iLpmax2,iLpmax3, iLppmin1,iLppmin2,iLppmin3,&
     iLppmax1,iLppmax2,iLppmax3, iLpppmin1,iLpppmin2,iLpppmin3,iLpppmax1,&
     iLpppmax2,iLpppmax3
  double precision                :: f_array(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:nw,3), d_array(3), beta(ixImin1:ixImax1,&
     ixImin2:ixImax2,ixImin3:ixImax3,1:nw,3)
  double precision                :: beta_coeff(2), alpha_array&
     (ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw,3),&
      tau(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw)
  double precision                :: alpha_sum(ixImin1:ixImax1,&
     ixImin2:ixImax2,ixImin3:ixImax3,1:nw), tmp(ixImin1:ixImax1,&
     ixImin2:ixImax2,ixImin3:ixImax3,1:nw), flux(ixImin1:ixImax1,&
     ixImin2:ixImax2,ixImin3:ixImax3,1:nw)
  double precision                :: lambda
  integer                         :: i, iw
  double precision, parameter     :: weno_eps_machine = 1.0d-42
  double precision, parameter     :: weno_dx_exp = 2.0d0/3.0d0
  !----------------------------------------------------------------------------

  ! iL^L holds the indices of interfaces to reconstruct to.  Convention is that a center index holds the _right-side_ interface.  

  iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
  iLmmin3=iLmin3-kr(idims,3);iLmmax1=iLmax1-kr(idims,1)
  iLmmax2=iLmax2-kr(idims,2);iLmmax3=iLmax3-kr(idims,3);
  iLmmmin1=iLmmin1-kr(idims,1);iLmmmin2=iLmmin2-kr(idims,2)
  iLmmmin3=iLmmin3-kr(idims,3);iLmmmax1=iLmmax1-kr(idims,1)
  iLmmmax2=iLmmax2-kr(idims,2);iLmmmax3=iLmmax3-kr(idims,3);
  iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
  iLpmin3=iLmin3+kr(idims,3);iLpmax1=iLmax1+kr(idims,1)
  iLpmax2=iLmax2+kr(idims,2);iLpmax3=iLmax3+kr(idims,3);
  iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
  iLppmin3=iLpmin3+kr(idims,3);iLppmax1=iLpmax1+kr(idims,1)
  iLppmax2=iLpmax2+kr(idims,2);iLppmax3=iLpmax3+kr(idims,3);
  iLpppmin1=iLppmin1+kr(idims,1);iLpppmin2=iLppmin2+kr(idims,2)
  iLpppmin3=iLppmin3+kr(idims,3);iLpppmax1=iLppmax1+kr(idims,1)
  iLpppmax2=iLppmax2+kr(idims,2);iLpppmax3=iLppmax3+kr(idims,3);

  d_array(1:3) = (/ 1.0d0/16.0d0, 10.0d0/16.0d0, 5.0d0/16.0d0 /)
  beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)

  lambda = dxdim**weno_dx_exp

  ! ==================================================
  ! Left side
  ! ==================================================

  f_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,1) = 3.0d0&
     /8.0d0 * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,iLmmmin3:iLmmmax3,&
     1:nwflux) - 10.0d0/8.0d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
     iLmmin3:iLmmax3,1:nwflux) + 15.0d0/8.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
     iLmin3:iLmax3,1:nwflux) !first stencil
  f_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,2) = -1.0d0&
     /8.0d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iLmmin3:iLmmax3,&
     1:nwflux) + 6.0d0/8.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
     1:nwflux) + 3.0d0/8.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
     iLpmin3:iLpmax3,1:nwflux)
  f_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,3) = 3.0d0&
     /8.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) + 6.0d0&
     /8.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,&
     1:nwflux) - 1.0d0/8.0d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
     iLppmin3:iLppmax3,1:nwflux)  

  
  beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,1) &
     = beta_coeff(1) * (w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
     iLmmmin3:iLmmmax3,1:nwflux) + w(iLmin1:iLmax1,iLmin2:iLmax2,&
     iLmin3:iLmax3,1:nwflux) - 2.0d0*w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
     iLmmin3:iLmmax3,1:nwflux))**2 + beta_coeff(2) * (w(iLmmmin1:iLmmmax1,&
     iLmmmin2:iLmmmax2,iLmmmin3:iLmmmax3,1:nwflux) - 4.0d0 * &
     w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iLmmin3:iLmmax3,1:nwflux) + &
     3.0d0*w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux))**2
  beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,2) &
     = beta_coeff(1) * (w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iLmmin3:iLmmax3,&
     1:nwflux) + w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,&
     1:nwflux) - 2.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
     1:nwflux))**2 + beta_coeff(2) * (w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
     iLmmin3:iLmmax3,1:nwflux) - w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
     iLpmin3:iLpmax3,1:nwflux))**2
  beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,3) &
     = beta_coeff(1) * (w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
     1:nwflux) + w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,iLppmin3:iLppmax3,&
     1:nwflux) - 2.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,&
     1:nwflux))**2 + beta_coeff(2) * (3.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
     iLmin3:iLmax3, 1:nwflux) - 4.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
     iLpmin3:iLpmax3,1:nwflux) + w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
     iLppmin3:iLppmax3,1:nwflux))**2

  tau(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = abs( &
     beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,&
     1) - beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,3) )
  
  alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = 0.0d0 
  do i = 1,3
     tmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) &
        = (tau(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) + &
        weno_eps_machine ) / ( beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
        1:nwflux,i) + weno_eps_machine )
     alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,i) &
        = d_array(i) * (1.0d0 + tmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
        1:nwflux)**2 + lambda/tmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
        1:nwflux) )
     alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) &
        = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
        1:nwflux) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
        1:nwflux,i)
  end do


  flux(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = 0.0d0
  do i = 1,3
     flux(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) &
        = flux(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) + &
        f_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,&
        i) * alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,&
        i)/(alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux))
  end do

  ! left value at right interface
  wLC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = &
     flux(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)

  ! ==================================================
  ! Right side
  ! ==================================================

  f_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,1) = 3.0d0&
     /8.0d0 * w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,iLpppmin3:iLpppmax3,&
     1:nwflux) - 10.0d0/8.0d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
     iLppmin3:iLppmax3,1:nwflux) + 15.0d0/8.0d0 * w(iLpmin1:iLpmax1,&
     iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux) !first stencil
  f_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,2) = -1.0d0&
     /8.0d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,iLppmin3:iLppmax3,&
     1:nwflux) + 6.0d0/8.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
     iLpmin3:iLpmax3,1:nwflux) + 3.0d0/8.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
     iLmin3:iLmax3,1:nwflux)
  f_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,3) = 3.0d0&
     /8.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,&
     1:nwflux) + 6.0d0/8.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
     1:nwflux) - 1.0d0/8.0d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
     iLmmin3:iLmmax3,1:nwflux)  

  beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,1) &
     = beta_coeff(1) * (w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
     iLpppmin3:iLpppmax3,1:nwflux) + w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
     iLpmin3:iLpmax3,1:nwflux) - 2.0d0*w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
     iLppmin3:iLppmax3,1:nwflux))**2 + beta_coeff(2) * (w(iLpppmin1:iLpppmax1,&
     iLpppmin2:iLpppmax2,iLpppmin3:iLpppmax3,1:nwflux) - 4.0d0 * &
     w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,iLppmin3:iLppmax3,&
     1:nwflux) + 3.0d0*w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,&
     1:nwflux))**2
  beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,2) &
     = beta_coeff(1) * (w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
     iLppmin3:iLppmax3,1:nwflux) + w(iLmin1:iLmax1,iLmin2:iLmax2,&
     iLmin3:iLmax3,1:nwflux) - 2.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
     iLpmin3:iLpmax3,1:nwflux))**2 + beta_coeff(2) * (w(iLppmin1:iLppmax1,&
     iLppmin2:iLppmax2,iLppmin3:iLppmax3,1:nwflux) - w(iLmin1:iLmax1,&
     iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux))**2
  beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,3) &
     = beta_coeff(1) * (w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,&
     1:nwflux) + w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iLmmin3:iLmmax3,&
     1:nwflux) - 2.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
     1:nwflux))**2 + beta_coeff(2) * (3.0d0 * w(iLpmin1:iLpmax1,&
     iLpmin2:iLpmax2,iLpmin3:iLpmax3, 1:nwflux) - 4.0d0 * w(iLmin1:iLmax1,&
     iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) + w(iLmmin1:iLmmax1,&
     iLmmin2:iLmmax2,iLmmin3:iLmmax3,1:nwflux))**2


  tau(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = abs( &
     beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,&
     1) - beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,3) )
  
  alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = 0.0d0 
  do i = 1,3
     tmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) &
        = (tau(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) + &
        weno_eps_machine ) / ( beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
        1:nwflux,i) + weno_eps_machine )
     alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,i) &
        = d_array(i) * (1.0d0 + tmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
        1:nwflux)**2 + lambda/tmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
        1:nwflux) )
     alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) &
        = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
        1:nwflux) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
        1:nwflux,i)
  end do

  flux(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = 0.0d0
  do i = 1,3
     flux(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) &
        = flux(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) + &
        f_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,&
        i) * alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux,&
        i)/(alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux))
  end do

  ! right value at right interface
  wRC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = &
     flux(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)

end subroutine WENOZPlimiter
!=============================================================================
subroutine WENOZPlimitervar(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,idims,dxdim,w,wLC,wRC)
  ! WENOZ and WENOZ+ Limiter, see Acker, Borges and Costa (2016)
  ! Needs at least three ghost cells.  Set dixB=3.

  use mod_amrvacdef

  integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3, idims
  double precision, intent(in)    :: dxdim
  double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3)
  double precision, intent(inout) :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3),wLC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) 
  ! .. local ..
  integer                         :: iLmmin1,iLmmin2,iLmmin3,iLmmax1,iLmmax2,&
     iLmmax3, iLmmmin1,iLmmmin2,iLmmmin3,iLmmmax1,iLmmmax2,iLmmmax3, iLpmin1,&
     iLpmin2,iLpmin3,iLpmax1,iLpmax2,iLpmax3, iLppmin1,iLppmin2,iLppmin3,&
     iLppmax1,iLppmax2,iLppmax3, iLpppmin1,iLpppmin2,iLpppmin3,iLpppmax1,&
     iLpppmax2,iLpppmax3
  double precision                :: f_array(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,3), d_array(3), beta(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,3)
  double precision                :: beta_coeff(2), alpha_array&
     (ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,3), tau(ixImin1:ixImax1,&
     ixImin2:ixImax2,ixImin3:ixImax3)
  double precision                :: alpha_sum(ixImin1:ixImax1,&
     ixImin2:ixImax2,ixImin3:ixImax3), tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3), flux(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
  double precision                :: lambda
  integer                         :: i, iw
  double precision, parameter     :: weno_eps_machine = 1.0d-42
  double precision, parameter     :: weno_dx_exp = 2.0d0/3.0d0
  !----------------------------------------------------------------------------

  ! iL^L holds the indices of interfaces to reconstruct to.  Convention is that a center index holds the _right-side_ interface.  

  iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
  iLmmin3=iLmin3-kr(idims,3);iLmmax1=iLmax1-kr(idims,1)
  iLmmax2=iLmax2-kr(idims,2);iLmmax3=iLmax3-kr(idims,3);
  iLmmmin1=iLmmin1-kr(idims,1);iLmmmin2=iLmmin2-kr(idims,2)
  iLmmmin3=iLmmin3-kr(idims,3);iLmmmax1=iLmmax1-kr(idims,1)
  iLmmmax2=iLmmax2-kr(idims,2);iLmmmax3=iLmmax3-kr(idims,3);
  iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
  iLpmin3=iLmin3+kr(idims,3);iLpmax1=iLmax1+kr(idims,1)
  iLpmax2=iLmax2+kr(idims,2);iLpmax3=iLmax3+kr(idims,3);
  iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
  iLppmin3=iLpmin3+kr(idims,3);iLppmax1=iLpmax1+kr(idims,1)
  iLppmax2=iLpmax2+kr(idims,2);iLppmax3=iLpmax3+kr(idims,3);
  iLpppmin1=iLppmin1+kr(idims,1);iLpppmin2=iLppmin2+kr(idims,2)
  iLpppmin3=iLppmin3+kr(idims,3);iLpppmax1=iLppmax1+kr(idims,1)
  iLpppmax2=iLppmax2+kr(idims,2);iLpppmax3=iLppmax3+kr(idims,3);

  d_array(1:3) = (/ 1.0d0/16.0d0, 10.0d0/16.0d0, 5.0d0/16.0d0 /)
  beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)

  lambda = dxdim**weno_dx_exp

  ! ==================================================
  ! Left side
  ! ==================================================

  f_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1) = 3.0d0&
     /8.0d0 * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,iLmmmin3:iLmmmax3) - &
     10.0d0/8.0d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iLmmin3:iLmmax3) + &
     15.0d0/8.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) !first stencil
  f_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,2) = -1.0d0&
     /8.0d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iLmmin3:iLmmax3) + 6.0d0&
     /8.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) + 3.0d0&
     /8.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3)
  f_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,3) = 3.0d0&
     /8.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) + 6.0d0&
     /8.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3) - 1.0d0&
     /8.0d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,iLppmin3:iLppmax3)  

  
  beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1) = beta_coeff(1) * &
     (w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,iLmmmin3:iLmmmax3) + &
     w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) - 2.0d0*w(iLmmin1:iLmmax1,&
     iLmmin2:iLmmax2,iLmmin3:iLmmax3))**2 + beta_coeff(2) * &
     (w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,iLmmmin3:iLmmmax3) - 4.0d0 * &
     w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iLmmin3:iLmmax3) + 3.0d0*w&
     (iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3))**2
  beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,2) = beta_coeff(1) * &
     (w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iLmmin3:iLmmax3) + w(iLpmin1:iLpmax1,&
     iLpmin2:iLpmax2,iLpmin3:iLpmax3) - 2.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
     iLmin3:iLmax3))**2 + beta_coeff(2) * (w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
     iLmmin3:iLmmax3) - w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3))**2
  beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,3) = beta_coeff(1) * &
     (w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) + w(iLppmin1:iLppmax1,&
     iLppmin2:iLppmax2,iLppmin3:iLppmax3) - 2.0d0 * w(iLpmin1:iLpmax1,&
     iLpmin2:iLpmax2,iLpmin3:iLpmax3))**2 + beta_coeff(2) * (3.0d0 * &
     w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) - 4.0d0 * w(iLpmin1:iLpmax1,&
     iLpmin2:iLpmax2,iLpmin3:iLpmax3) + w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
     iLppmin3:iLppmax3))**2

  tau(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = abs( beta(iLmin1:iLmax1,&
     iLmin2:iLmax2,iLmin3:iLmax3,1) - beta(iLmin1:iLmax1,iLmin2:iLmax2,&
     iLmin3:iLmax3,3) )
  
  alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = 0.0d0 
  do i = 1,3
     tmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = (tau(iLmin1:iLmax1,&
        iLmin2:iLmax2,iLmin3:iLmax3) + weno_eps_machine ) / ( &
        beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,i) + weno_eps_machine )
     alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,i) &
        = d_array(i) * (1.0d0 + tmp(iLmin1:iLmax1,iLmin2:iLmax2,&
        iLmin3:iLmax3)**2 + lambda/tmp(iLmin1:iLmax1,iLmin2:iLmax2,&
        iLmin3:iLmax3) )
     alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = &
        alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) + &
        alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,i)
  end do


  flux(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = 0.0d0
  do i = 1,3
     flux(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = flux(iLmin1:iLmax1,&
        iLmin2:iLmax2,iLmin3:iLmax3) + f_array(iLmin1:iLmax1,iLmin2:iLmax2,&
        iLmin3:iLmax3,i) * alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,&
        iLmin3:iLmax3,i)/(alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
        iLmin3:iLmax3))
  end do

  ! left value at right interface
  wLC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = flux(iLmin1:iLmax1,&
     iLmin2:iLmax2,iLmin3:iLmax3)

  ! ==================================================
  ! Right side
  ! ==================================================

  f_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1) = 3.0d0&
     /8.0d0 * w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,iLpppmin3:iLpppmax3) &
     - 10.0d0/8.0d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
     iLppmin3:iLppmax3) + 15.0d0/8.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
     iLpmin3:iLpmax3) !first stencil
  f_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,2) = -1.0d0&
     /8.0d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,iLppmin3:iLppmax3) + &
     6.0d0/8.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3) + 3.0d0&
     /8.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3)
  f_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,3) = 3.0d0&
     /8.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3) + 6.0d0&
     /8.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) - 1.0d0&
     /8.0d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iLmmin3:iLmmax3)  

  beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1) = beta_coeff(1) * &
     (w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,iLpppmin3:iLpppmax3) + &
     w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3) - 2.0d0*w&
     (iLppmin1:iLppmax1,iLppmin2:iLppmax2,iLppmin3:iLppmax3))&
     **2 + beta_coeff(2) * (w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
     iLpppmin3:iLpppmax3) - 4.0d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
     iLppmin3:iLppmax3) + 3.0d0*w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
     iLpmin3:iLpmax3))**2
  beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,2) = beta_coeff(1) * &
     (w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,iLppmin3:iLppmax3) + &
     w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) - 2.0d0 * w(iLpmin1:iLpmax1,&
     iLpmin2:iLpmax2,iLpmin3:iLpmax3))**2 + beta_coeff(2) * &
     (w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,iLppmin3:iLppmax3) - &
     w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3))**2
  beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,3) = beta_coeff(1) * &
     (w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3) + w(iLmmin1:iLmmax1,&
     iLmmin2:iLmmax2,iLmmin3:iLmmax3) - 2.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
     iLmin3:iLmax3))**2 + beta_coeff(2) * (3.0d0 * w(iLpmin1:iLpmax1,&
     iLpmin2:iLpmax2,iLpmin3:iLpmax3) - 4.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
     iLmin3:iLmax3) + w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iLmmin3:iLmmax3))**2


  tau(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = abs( beta(iLmin1:iLmax1,&
     iLmin2:iLmax2,iLmin3:iLmax3,1) - beta(iLmin1:iLmax1,iLmin2:iLmax2,&
     iLmin3:iLmax3,3) )
  
  alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = 0.0d0 
  do i = 1,3
     tmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = (tau(iLmin1:iLmax1,&
        iLmin2:iLmax2,iLmin3:iLmax3) + weno_eps_machine ) / ( &
        beta(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,i) + weno_eps_machine )
     alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,i) &
        = d_array(i) * (1.0d0 + tmp(iLmin1:iLmax1,iLmin2:iLmax2,&
        iLmin3:iLmax3)**2 + lambda/tmp(iLmin1:iLmax1,iLmin2:iLmax2,&
        iLmin3:iLmax3) )
     alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = &
        alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) + &
        alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,i)
  end do

  flux(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = 0.0d0
  do i = 1,3
     flux(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = flux(iLmin1:iLmax1,&
        iLmin2:iLmax2,iLmin3:iLmax3) + f_array(iLmin1:iLmax1,iLmin2:iLmax2,&
        iLmin3:iLmax3,i) * alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,&
        iLmin3:iLmax3,i)/(alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
        iLmin3:iLmax3))
  end do

  ! right value at right interface
  wRC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = flux(iLmin1:iLmax1,&
     iLmin2:iLmax2,iLmin3:iLmax3)

end subroutine WENOZPlimitervar
!=============================================================================

  !=============================================================================

  

end module mod_limiter
!=============================================================================
