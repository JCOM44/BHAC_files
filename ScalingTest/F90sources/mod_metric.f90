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
!
! Metric and geometry for general covariant coordinates.
! To use, set the 'typecoord' keyword to a string containting 'covariant'.  
! The metric is initialised by providing values for alpha, shift-beta and 
! the spatial part of the metric.
! 
! For the dynamics, the metric is assumed to be of the form
! ds^2 = (-\alpha^2+\beta_i \beta^i)dt^2 + 2 \beta_i dt dx^i + \gamma_{ij} dx^i dx^j
! where \beta_i = \gamma_{ij}beta^j.  
! See e.g. the book of Alcubierre, Introduction to 3+1 Numerical Relativity (2008)
! Chapter 2.2
!
! The Valencia forumation is used for the evolution with first class citizens:
! Spatial metric \gamma_{ij}
! Lapse function \alpha
! Contravariant shift vector \beta^i
! and their spatial partial derivatives.  
! 
! It is also possible to initialise the four-metric first and the code will
! automatically obtain lapse alpha and shift beta which are needed for the
! dynamical evolution.
! To do so, set init_from_g4 = .true. in your corresponding mod_coord_ file.
! 
! In the rest of the code, the metric is available via the pointer myM
! or via the geometry datastructure.
! 
! For the definition of the datastructure see 'mod_physicaldata.t'.
! 2014-03-23 by Oliver Porth
! 
!=============================================================================
module mod_metric
  use mod_metric_aux
  implicit none

  integer, save                        :: nnonzero_metric !number of non-zero elements in four-metric
  integer, save                        :: nnonzero_beta !number of non-zero elements in contravariant shift
  integer, save                        :: nnonzero_dalphadj !number of non-zero elements in lapse-derivatives  
  integer, save                        :: nnonzero_dbetaidj !number of non-zero elements in shift-derivatives  
  integer, save                        :: nnonzero_dgdk !number of non-zero elements in metric-derivatives
  logical, save,dimension(0:3,0:3) :: g_is_zero       ! four-metric
  logical, save,dimension(1:3)       :: beta_is_zero !contravariant part, not in metric
  logical, save,dimension(1:3,1:3,1:3) :: dgdk_is_zero !partial derivatives of three-metric
  logical, save,dimension(1:3)       :: dalphadj_is_zero !partial derivatives of lapse
  logical, save,dimension(1:3,1:3) :: dbetaidj_is_zero !partial derivatives of shift 
  logical, save                        :: space_is_diagonal,&
      sqrtgamma_is_analytic

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
  ! Metric components for Cartesian coordinates -coord=cart
  ! This file is included by the preprocessor in mod_metric.t
  ! It will be part of the module mod_metric.
  !  
  ! Oliver Porth
  ! 28.06.2016
  !=============================================================================
  
  
  !-----------------------------------------
  ! Define constants specific to coordinates
  !-----------------------------------------
  
  character*20, parameter              :: coord="cart"
  logical, save                        :: init_from_g4 = .false. !We don't provide the four-metric

  ! Migrating to have coordinate-specific constants here.
  ! a_ and m_ however are in the physics part
  ! so put the MKS parameters, e.g. R0_ and h_ here:
  integer, parameter                   :: ncoordpar = 0
  double precision, save               :: coordpar(ncoordpar) 

contains
!=============================================================================
  subroutine init_coord

    SPToCoord=>BLToCoord
    
  end subroutine init_coord
  !=============================================================================
  subroutine get_sqrtgamma_analytic(x1,x2,x3,sqrtgamma,is_analytic)

    use mod_amrvacdef

    ! You can specify sqrtgamma here.
    ! If you don't want to, set is_analytic = .false.
    ! In that case, get_sqrtgamma() will calculate from the metric,
    ! which can be a major slowdown.
    ! init_metric() will set sqrtgamma_is_analytic from the return-value
    ! of "is_analytic".  

    double precision, intent(in)                     :: x1,x2,x3
    double precision, intent(out)                    :: sqrtgamma
    logical, optional                                :: is_analytic
    !-----------------------------------------------------------------------------


    sqrtgamma = one

    if(present(is_analytic)) is_analytic = .true.

  end subroutine get_sqrtgamma_analytic
  !=============================================================================
  subroutine get_alpha(x1,x2,x3,alpha,iszero,dalphadj_iszero,dalphadj,jdir)

    use mod_amrvacdef

    ! get the lapse.  Optional parameter is true if lapse is
    ! identically zero (does not really make sense)
    ! Optional parameters jdir and dalphadj request derivatives
    ! \partial_j \alpha ; j=jdir
    double precision, intent(in)                     :: x1,x2,x3
    integer, optional, intent(in)                    :: jdir
    double precision, intent(out)                    :: alpha
    logical, optional, intent(out)                   :: iszero,&
        dalphadj_iszero
    double precision, optional, intent(out)          :: dalphadj
    ! .. local ..
    !-----------------------------------------------------------------------------
    if(present(dalphadj) .and. .not. present(jdir) .or. present&
       (dalphadj_iszero) .and. .not. present(jdir)) call mpistop("get_alpha: &
       derivatives requested without direction or output-slot given.")


    alpha = one

    if(present(iszero)) iszero = .false.
    if (present(dalphadj_iszero)) dalphadj_iszero = .true.
    if (present(dalphadj)) dalphadj = zero

  end subroutine get_alpha
  !=============================================================================
  subroutine get_beta(idir,x1,x2,x3,beta,iszero,dbetaidj_iszero,dbetaidj,jdir)

    use mod_amrvacdef

    ! get the (contravariant!!) shift vector.
    ! The optional argument iszero is true if shift-component is 
    ! identically zero.
    ! if requested, dbetaidj is the derivative of the contravariant shift.
    ! \partial_j \beta^i ; i=idir, j=jdir
    integer, intent(in)                      :: idir
    double precision, intent(in)             :: x1,x2,x3
    integer, optional, intent(in)            :: jdir
    double precision, intent(out)            :: beta
    logical, optional, intent(out)           :: iszero, dbetaidj_iszero
    double precision, optional, intent(out)  :: dbetaidj
    ! .. local ..
    !-----------------------------------------------------------------------------
    if(present(dbetaidj) .and. .not. present(jdir) .or. present&
       (dbetaidj_iszero) .and. .not. present(jdir)) call mpistop("get_beta: &
       derivatives requested &without direction or output-slot given.")

    beta = zero

    if(present(iszero)) iszero = .true.
    ! \partial_j \beta^i
    if (present(dbetaidj)) dbetaidj = 0.0d0
    if (present(dbetaidj_iszero)) dbetaidj_iszero = .true.

  end subroutine get_beta
  !=============================================================================
  subroutine get_g_component(iin,jin,x1,x2,x3,g,iszero,dgdk_iszero,dgdk,kdir)

    use mod_amrvacdef

    ! This is at the heart of the scheme: Set the (spatial) metric components here
    ! and only here...
    ! Indices of the metric are down (covariant) g_{ij}
    ! The optional argument iszero is true if the element is identically zero
    ! The optional arguments dgdk and kdir request derivatives of the metric
    ! \partial_k g_{ij} ; i=iin, j=jin, k=kdir
    ! The optional argument dgdk_iszero
    integer, intent(in)                      :: iin,jin
    integer, optional, intent(in)            :: kdir
    double precision, intent(in)             :: x1,x2,x3
    double precision, intent(out)            :: g
    logical, optional, intent(out)           :: iszero, dgdk_iszero
    double precision, optional, intent(out)  :: dgdk
    ! .. local ..
    integer                                  :: i,j
    !-----------------------------------------------------------------------------
    if(present(dgdk) .and. .not. present(kdir) .or. present(dgdk_iszero) &
       .and. .not. present(kdir)) call mpistop("get_g_component: derivatives &
       requested without &direction or output-slot given.")

    ! metric is symmetric: swap indices if needed:
    ! User needs only to provide values for i<=j (upper triangle).  
    if (iin>jin) then
       i=jin; j=iin
    else
       i=iin; j=jin
    end if

    ! Diagonal unity
    if (i==j) then
       if(present(iszero)) iszero = .false.

       g = one 
 
    else
       if(present(iszero)) iszero = .true.

       g = zero

    end if


    if (present(dgdk)) dgdk = zero
    if(present(dgdk_iszero)) dgdk_iszero = .true.

  end subroutine get_g_component
  !=============================================================================
    double precision function outerhorizon()
    
    !-----------------------------------------------------------------------------

    outerhorizon = 0.0d0
    
  end function outerhorizon
  !=============================================================================
  subroutine BLToCoord(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xBL,xCoord)

    ! Assuming zero mass and do the transformation for spherical coordinates
    !
    use mod_amrvacdef

    integer,intent(in)                                     :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: xBL
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(out)  :: xCoord
    ! .. local ..
    !-----------------------------------------------------------------------------

    xCoord(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) &
       = xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)  * cos(xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2))  * sin(xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))

    if (3 .lt. 3) then
       ! Simulation in a plane
       if (2 .le. 3 .and. 2 .gt. 0) then
          
          !Simulation in the r-phi plane:
          xCoord(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) &
             = xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             1)  * sin(xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             2))  * sin(xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             3))
         
       else if (3 .le. 3 .and. 3 .gt. 0) then
          !Simulation in the r-z plane:
          
          xCoord(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) &
             = xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             1)  * cos(xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))
         
       end if
    else
       ! Simulation in 3D
       
       xCoord(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) &
          = xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          1)  * sin(xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          2))  * sin(xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))
       xCoord(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3) &
          = xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          1)  * cos(xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))
       
    end if
    
  end subroutine BLToCoord
  !=============================================================================
  subroutine CoordToBL(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xCoord,xBL)

    !
    use mod_amrvacdef

    integer,intent(in)                                     :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: xCoord
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(out)  :: xBL
    ! .. local ..
    !-----------------------------------------------------------------------------

    xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) &
       = sqrt( xCoord(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)&
       **2 + xCoord(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)&
       **2 + xCoord(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)**2 )

    if (3 .lt. 3) then
       ! Simulation in a plane
       if (2 .le. 3 .and. 2 .gt. 0) then 
          !Simulation in the r-phi plane:
           xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) &
              = atan2(xCoord(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              2),xCoord(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1))
       else if (3 .le. 3 .and. 3 .gt. 0) then
          !Simulation in the r-z plane:
           xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3) &
              = acos(xCoord(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              3)) / xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)
       end if
    else
       ! Simulation in 3D
        xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) &
           = atan2(xCoord(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2),&
           xCoord(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1))
        xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)   &
           = acos(xCoord(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)) &
           / xBL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)
    end if
    
  end subroutine CoordToBL
  !=============================================================================
  subroutine u4BLtoCoord(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,u4BL,u4Coord)

    ! Assuming zero mass and do the transformation for spherical coordinates
    ! Most likely we will not use this ever.  
    !
    use mod_amrvacdef

    integer                                                :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: x
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(in)   :: u4BL
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(out)  :: u4Coord
    ! .. local ..
    !-----------------------------------------------------------------------------

    if (eqpar(m_) .gt. zero) call mpistop("BLToCoord: Cannot transform to &
       Cartesian for non-zero mass")

    call mpistop("u4BLtoCoord: Transformation not implemented")

    u4Coord = u4BL
    
  end subroutine u4BLtoCoord
  !============================================================================= 
  subroutine u3CoordToCart(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,u3Coord,u3Cart,J)

    ! u3Cart = u3Coord
    !-----------------------------------------------------------------------------

    use mod_amrvacdef

    integer, intent(in)                                    :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: x
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: u3Coord
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(out)  :: u3Cart
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3,1:3), optional, intent(out)  :: J
    ! .. local ..
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3,1:3)         :: Jac
    integer                                                :: ix, jx
    !-----------------------------------------------------------------------------

    u3Cart = u3Coord
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Assemble (diagonal) Jacobian:
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (present(J)) then
       Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:,:) = zero
       do ix=1,ndir
          do jx=1,ndir
             Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ix,jx) = one
          end do
       end do
       J = Jac
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  end subroutine u3CoordToCart
  !=============================================================================
  subroutine CoordToCart(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,xCart)

    ! identity
    !-----------------------------------------------------------------------------

    use mod_amrvacdef
    integer, intent(in)                                    :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: x
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(out)  :: xCart
    ! .. local ..
    !-----------------------------------------------------------------------------

    xCart = x

  end subroutine CoordToCart
  !=============================================================================
  subroutine get_dgammainvdk(x1,x2,x3,dgammainvdk,k)
    ! Obtains the derivatives of the inverse _spatial_ metric
    ! \partial_k gamma^{ij} for a given k. 
    ! This is not required when initialising gamma, lapse and shift
    ! e.g. when set init_from_g4 = .false.

    double precision, intent(in)           :: x1,x2,x3
    double precision, intent(out)          :: dgammainvdk(1:3,1:3)
    integer, intent(in)                    :: k
    ! .. local ..
    !-----------------------------------------------------------------------------

    dgammainvdk = 0.0d0
    
  end subroutine get_dgammainvdk
  !=============================================================================
  ! End of coordinate-specific definitions.
  !=============================================================================


!=============================================================================
subroutine init_coord_check

  if (.not. associated(get_gammainv_component_analytic)) &
     get_gammainv_component_analytic => dummy_get_gammainv_component_analytic
  if (.not. associated(get_gammainv_component)) get_gammainv_component &
     => get_gammainv_component_numerical
  if (.not. associated(SPToCoord)) SPToCoord => dummy_SPToCoord
  if (.not. associated(LRNFu)) LRNFu => dummy_LRNFu !locally non-rotation orthonormal tetrad transform
  if (.not. associated(u4CoordToKS)) u4CoordToKS => dummy_u4CoordToKS !Transform four-vector to KS coordinates
  if (.not. associated(u4CoordToBL)) u4CoordToBL => dummy_u4CoordToBL !Transform four-vector to BL coordinates
  if (.not. associated(d4CoordToKS)) d4CoordToKS => dummy_d4CoordToKS !Transform four-covariant vector to KS coordinates

end subroutine init_coord_check
!=============================================================================
subroutine get_g4(x1,x2,x3,g,dgdk)
  ! Return the four-metric at point x^D
  ! This is just a convenience routine and fairly slow.
  ! Can also return metric derivatives when dgdk given.  
  ! Avoid extensive use, expecially when init_from_g4=.false.
  double precision, intent(in)                           :: x1,x2,x3
  double precision, dimension(0:3,0:3),intent(out)   :: g
  double precision, dimension(0:3,0:3,1:3), optional, intent(out) :: dgdk
  ! .. local ..
  integer                                                :: i, j, k
  double precision, dimension(1:3)                     :: betaU, betaD,&
      dalphadj
  double precision, dimension(1:3,1:3)               :: dbetaUdj, dbetaDdj
  double precision                                       :: beta2, alpha,&
      dummy
  !-----------------------------------------------------------------------------

  if (.not. init_from_g4) then

     do j=1,3
        do i=1,3
           if (present(dgdk)) then
              do k = 1, 3
                 call get_g_component(i,j,x1,x2,x3,g(i,j),dgdk&
                    =dgdk(i,j,k),kdir=k)
              end do
           else
              call get_g_component(i,j,x1,x2,x3,g(i,j))
           end if
        end do
     end do

     do i=1,3
        call get_alpha(x1,x2,x3,alpha,dalphadj=dalphadj(i),jdir=i)
        if (present(dgdk)) then
           do j=1,3
              call get_beta(i,x1,x2,x3,betaU(i),dbetaidj=dbetaUdj(i,j),jdir=j)
           end do
        else
           call get_beta(i,x1,x2,x3,betaU(i))   
        end if
     end do
     ! lower the beta:
     betaD(:) = 0.0d0
     do j=1,3
        do i=1,3
           if (g_is_zero(i,j) .or. beta_is_zero(j)) cycle
           betaD(i) = betaD(i) + g(i,j)*betaU(j)
        end do
     end do

     if (present(dgdk)) then
        do j=1,3
           do i=1,3
              dbetaDdj(i,j) =  betaU(1) * dgdk(i,1,j) + betaU(2) * dgdk(i,2,&
                 j) + betaU(3) * dgdk(i,3,j)  +  g(i,1)*dbetaUdj(1,j) + g(i,&
                 2)*dbetaUdj(2,j) + g(i,3)*dbetaUdj(3,j) 
           end do
        end do
     end if

     beta2 = betaU(1)*betaD(1)+betaU(2)*betaD(2)+betaU(3)*betaD(3)

     g(0,0) = -alpha**2 + beta2

     do i=1,3
        g(i,0) = betaD(i)
        g(0,i) = betaD(i)
     end do

     if (present(dgdk)) then
        do i=1,3
           dgdk(0,0,i) = -2.0d0*alpha*dalphadj(i) +  betaD(1)*dbetaUdj(1,&
              i) + betaD(2)*dbetaUdj(2,i) + betaD(3)*dbetaUdj(3,&
              i)  +  betaU(1)*dbetaDdj(1,i) + betaU(2)*dbetaDdj(2,&
              i) + betaU(3)*dbetaDdj(3,i) 
        end do
        do i=1,3
           do j=1,3
              dgdk(0,i,j) = dbetaDdj(i,j)
              dgdk(i,0,j) = dbetaDdj(i,j)
           end do
        end do
     end if

  else ! init_from_g4

     do j=0,3
        do i=0,3
           if (present(dgdk)) then
              do k = 1, 3
                 call get_g_component(i,j,x1,x2,x3,g(i,j),dgdk&
                    =dgdk(i,j,k),kdir=k)
              end do
           else
              call get_g_component(i,j,x1,x2,x3,g(i,j))
           end if
        end do
     end do

  end if ! init_from_g4

end subroutine get_g4
!=============================================================================
subroutine g4inv(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
   ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,m,ginv)
  use mod_amrvacdef
  ! Calculates the full inverse 4-metric using 3+1 quantities in m
  integer, intent(in)                      :: ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
  type(metric), pointer                    :: m
  double precision, intent(out), dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,0:ndir,0:ndir)  :: ginv
  ! .. local ..
  integer                                  :: i,j
  !-----------------------------------------------------------------------------

  ginv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,0,0) = -1.0d0&
     /m%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)**2
  do j=1,3
     do i=0,j
        if (i .eq. 0) then
           ginv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,0,j) &
              = m%beta(j)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)/m%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)**2
        else
           ginv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,i,j) &
              = m%gammainv(i,j)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3) - m%beta(i)%elem(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3)*m%beta(j)%elem(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3) / m%alpha(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3)**2
        end if
        if (i .ne. j) ginv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,j,&
           i) = ginv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,i,j)
     end do
  end do
  
end subroutine g4inv
!=============================================================================
subroutine get_sqrtgamma(x1,x2,x3,sqrtgamma)
    ! Calculate the determinant of the spatial metric at point x
    ! Takes analytic values if provided, otherwise bruteforcing...
    double precision, intent(in)                            :: x1,x2,x3
    double precision, intent(out)                           :: sqrtgamma
    double precision, dimension(1:3,1:3)                :: g
    ! .. local ..
    integer                                                 :: i, j
    !-----------------------------------------------------------------------------
    
    if (sqrtgamma_is_analytic) then

       call get_sqrtgamma_analytic(x1,x2,x3,sqrtgamma)
       
    else

       do j=1,3
          do i=1,j
             call get_g_component(i,j,x1,x2,x3,g(i,j))
             if (i .ne. j) g(j,i) = g(i,j)
          end do
       end do

       if (space_is_diagonal) then 
          ! If space is diagonal, just multiply diagonal entries
          sqrtgamma = 1.0d0
          do j=1,3
              sqrtgamma = sqrtgamma*g(j,j)
          end do
       else
          ! Else, compute the determinant by LU decomposition
          sqrtgamma = determinant(g,3)
       end if

       ! Check that det(gamma) is positive and
       ! send warning message if not.
       ! HO: I commented this since it spams the standard output
!       if (sqrtgamma.lt.0.0d0) &
!          print *, 'WARNING: metric determinant is negative at x^D =',&
!                    x^D,'taking absolute value'

       sqrtgamma = sqrt(abs(sqrtgamma))
    end if

  end subroutine get_sqrtgamma
  !=============================================================================
  subroutine get_gammainv_component_numerical(iin,jin,x1,x2,x3,ginv,iszero,&
     dginvdk_iszero,dginvdk,kdir)

    ! Obtains a component of the inverse spatial metric gamma
    ! Not particularly efficient.
    
    use mod_lu
    use mod_amrvacdef

    integer, intent(in)                      :: iin,jin
    integer, optional, intent(in)            :: kdir
    double precision, intent(in)             :: x1,x2,x3
    double precision, intent(out)            :: ginv
    logical, optional, intent(out)           :: iszero, dginvdk_iszero
    double precision, optional, intent(out)  :: dginvdk
    ! .. local ..
    integer, dimension(1:3)                :: indx
    integer                                  :: i, j, code, d
    double precision, dimension(1:3,1:3) :: g
    double precision, dimension(1:3)       :: b
    double precision                         :: dgammainvdk(1:3,1:3)
    !-----------------------------------------------------------------------------


    if (space_is_diagonal) then 

       if (iin .eq. jin) then
          call get_g_component(iin,iin,x1,x2,x3,g(iin,iin))
          ginv = 1.0d0/g(iin,iin)
       else
          ginv = 0.0d0
       end if

    else ! space is diagonal

       do j=1,3
          do i=1,j
             call get_g_component(i,j,x1,x2,x3,g(i,j))
             if ( i .ne. j) g(j,i) = g(i,j)
          end do
       end do

       
       call LUDCMP(g,3,indx,d,code)

       b(1)=kr(jin,1);b(2)=kr(jin,2);b(3)=kr(jin,3)

       call LUBKSB(g,3,indx,b)

       ginv = b(iin)

    end if ! space is diagonal

    if (present(iszero)) iszero = .false. !At least we cant tell with certainty
    if (present(dginvdk_iszero)) dginvdk_iszero = .false. !At least we cant tell with certainty
    if (present(dginvdk)) then
       call get_dgammainvdk(x1,x2,x3,dgammainvdk,kdir) !this might actually be analytic or numerical
       dginvdk = dgammainvdk(iin,jin)
    end if
    
  end subroutine get_gammainv_component_numerical
  !=============================================================================
  subroutine get_alpha_from_g4(x1,x2,x3,myalpha,dalphadj,jdir)

    use mod_amrvacdef

    ! get the lapse from g4.  
    ! Optional parameters jdir and dalphadj request derivatives
    ! \partial_j \alpha ; j=jdir
    double precision, intent(in)                     :: x1,x2,x3
    integer, optional, intent(in)                    :: jdir
    double precision, intent(out)                    :: myalpha
    double precision, optional, intent(out)          :: dalphadj
    ! .. local ..
    integer                                          :: icomp
    double precision                                 :: betad(1:ndir),&
        betau(1:ndir)
    double precision                                 :: b2, gtt, dgttdjdir
    double precision                                 :: dbetaudjdir(1:ndir),&
        dbetaddjdir(1:ndir)
    double precision                                 :: dummy
    !-----------------------------------------------------------------------------
    if(present(dalphadj) .and. .not. present(jdir)) call mpistop&
       ("get_alpha_from_g4: derivatives requested without direction or &
       output-slot given.")

    do icomp = 1, ndir
       call get_g_component(0,icomp,x1,x2,x3,betad(icomp))
    end do
    call raise3_point(x1,x2,x3,betad,betau)
    b2 = sum(betau*betad)
    call get_g_component(0,0,x1,x2,x3,gtt)
    myalpha = sqrt(max(b2-gtt,0.0d0))

    if (present(jdir) .and. present(dalphadj)) then
       call get_g_component(0,0,x1,x2,x3,dummy,dgdk=dgttdjdir,kdir=jdir)
       do icomp = 1, ndir
          call get_g_component(0,icomp,x1,x2,x3,dummy,dgdk=dbetaddjdir(icomp),&
             kdir=jdir)
          call get_beta_from_g4(icomp,x1,x2,x3,dummy,dbetaidj&
             =dbetaudjdir(icomp),jdir=jdir)
       end do
       dalphadj = 0.5d0 * ( sum(betau*dbetaddjdir + betad*dbetaudjdir) - &
          dgttdjdir ) / myalpha
    end if

  end subroutine get_alpha_from_g4
  !=============================================================================
  subroutine get_beta_from_g4(idir,x1,x2,x3,mybeta,dbetaidj,jdir)

    use mod_amrvacdef

    ! get the (contravariant!!) shift vector from g4.
    ! if requested, dbetaidj is the derivative of the contravariant shift.
    ! \partial_j \beta^i ; i=idir, j=jdir
    integer, intent(in)                      :: idir
    double precision, intent(in)             :: x1,x2,x3
    integer, optional, intent(in)            :: jdir
    double precision, intent(out)            :: mybeta
    double precision, optional, intent(out)  :: dbetaidj
    ! .. local ..
    double precision                         :: betad(ndir), betau(ndir)
    integer                                  :: icomp, inonzero, i,j,k
    double precision                         :: dgammainvdjdir(1:ndir,1:ndir)
    double precision                         :: dbetaDownkDjdir
    double precision                         :: gammainvik
    double precision                         :: dummy
    !-----------------------------------------------------------------------------
    if(present(dbetaidj) .and. .not. present(jdir)) call mpistop&
       ("get_beta_from_g4: derivatives requested &without direction or &
       output-slot given.")

    do icomp = 1, ndir
       call get_g_component(0,icomp,x1,x2,x3,betad(icomp))
    end do
    call raise3_point(x1,x2,x3,betad,betau)
    mybeta = betau(idir)

    if (present(jdir) .and. present(dbetaidj)) then

       call get_dgammainvdk(x1,x2,x3,dgammainvdjdir,jdir)

       dbetaidj = 0.0d0
       do k=1,ndir
          call get_g_component(0,k,x1,x2,x3,dummy,dgdk=dbetaDownkDjdir,kdir&
             =jdir)
          call get_gammainv_component(idir,k,x1,x2,x3,gammainvik)
          dbetaidj = dbetaidj + betad(k)*dgammainvdjdir(idir,&
             k) + dbetaDownkDjdir*gammainvik
       end do

    end if

  end subroutine get_beta_from_g4
  !=============================================================================
  double precision function determinant(g,n)
    ! determinant for up to nxn matrices
    use mod_lu
    integer, intent(in) :: n
    double precision, dimension(1:n,1:n), intent(in)   :: g
    ! .. local ..
    double precision :: det
    integer          :: idx
    integer          :: code,indx(3),d
    !-----------------------------------------------------------------------------

    ! Compute the determinant as the product of the diagonal
    ! entries of U in the LU decompostion
    ! (the diagonal of L just contains ones).
    ! The integer d=+-1 keeps track of the permutations and
    ! changes the sign accortingly.

    call ludcmp_d(g,n,indx,d,code)

    determinant = d
    do idx=1,n
       determinant = determinant*g(idx,idx)
    end do

  end function determinant
  !=============================================================================
  subroutine get_sqrtgammai(idims,x1,x2,x3,s)
    ! Calculate the determinant of the induced metric in direction idims
    
    integer, intent(in)                                 :: idims
    double precision,intent(in)                         :: x1,x2,x3
    double precision, intent(out)                       :: s
    ! .. local ..
    integer                                             :: ii,jj,i,j

    double precision, dimension(1:3-1,1:3-1)        :: gi

    !-----------------------------------------------------------------------------




    i=0
    do ii=1,3
       if (ii.eq.idims) cycle
       i = i+1
       j = 0
       do jj=1,3
          if (jj.eq.idims) cycle
          j = j+1
          call get_g_component(ii,jj,x1,x2,x3,gi(i,j))
       end do
    end do



    s = sqrt(gi(1,1)*gi(2,2) - gi(1,2)*gi(2,1))

  end subroutine get_sqrtgammai
  !=============================================================================
  subroutine init_metric()
    use mod_amrvacdef
    
    ! .. local ..
    integer                            :: i,j,k,inonzero
    double precision                   :: dummy
    !-----------------------------------------------------------------------------
    

    ! Initialize the auxiliary subroutines (better strategy):
    call init_coord
    call init_coord_check
    
    ! Check if we want to use analytic expression for determinant:
    call get_sqrtgamma_analytic(xprobmin1,xprobmin2,xprobmin3,dummy,&
       is_analytic=sqrtgamma_is_analytic)
    
    ! Check which components from the spatial part are identical zero:
    do j=1,3
       do i=1,3
          call get_g_component(i,j,xprobmin1,xprobmin2,xprobmin3,dummy,iszero&
             =g_is_zero(i,j))
       end do
    end do

    ! Time component 00 cannot be zero:
    g_is_zero(0,0) = .false.

    ! Contravariant shifts:
    inonzero = 0
    do i=1,3
       call get_beta(i,xprobmin1,xprobmin2,xprobmin3,dummy,beta_is_zero(i))
       if (.not.beta_is_zero(i)) inonzero = inonzero+1
    end do
    nnonzero_beta = inonzero

    ! Time-space components elements of metric
    ! Check if lowering beta would still result in zero:
    do i=1,3
       g_is_zero(0,i) = .true.
       do j=1,3
          if ( .not.(g_is_zero(i,j) .or. beta_is_zero(j)) ) then
             g_is_zero(0,i) = .false.
          end if
       end do
       g_is_zero(i,0) = g_is_zero(0,i)
    end do
       

    ! Check if space is diagonal and set the flag
    ! also count the number of non-zero elements
    space_is_diagonal = .true.
    inonzero = 0
    do j=0,3
       do i=0,3
          if (.not.g_is_zero(i,j)) inonzero = inonzero + 1
          if (i.eq.j .or. i.eq.0 .or. j.eq.0) cycle
          space_is_diagonal = space_is_diagonal .and. g_is_zero(i,j)
       end do
    end do
    nnonzero_metric = inonzero
    
    ! Derivatives of lapse:
    inonzero = 0
    do j=1,3
       call get_alpha(xprobmin1,xprobmin2,xprobmin3,dummy,dalphadj_iszero&
          =dalphadj_is_zero(j),jdir=j)
       if (.not. dalphadj_is_zero(j)) inonzero = inonzero + 1
    end do
    nnonzero_dalphadj = inonzero
    
    ! Derivatives of shift:
    inonzero = 0
    do j=1,3
       do i=1,3
          call get_beta(i,xprobmin1,xprobmin2,xprobmin3,dummy,dbetaidj_iszero&
             =dbetaidj_is_zero(i,j),jdir=j)
          if (.not. dbetaidj_is_zero(i,j)) inonzero = inonzero + 1
       end do
    end do
    nnonzero_dbetaidj = inonzero
    
    ! Derivatives of three-metric:
    inonzero = 0
    do k=1,3
       do j=1,3
          do i=1,3
             call get_g_component(i,j,xprobmin1,xprobmin2,xprobmin3,dummy,&
                dgdk_iszero=dgdk_is_zero(i,j,k),kdir=k)
             if (.not. dgdk_is_zero(i,j,k)) inonzero = inonzero + 1
          end do
       end do
    end do
    nnonzero_dgdk = inonzero
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Done collecting information on spacetime !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (mype.eq.0) print*,'-----------------------------------------------------------------------------'
    if (mype.eq.0) write(*,*) 'Initialised Meta-info on metric:'
    if (mype.eq.0) print*,'-----------------------------------------'
    if (mype.eq.0) write(*,*) 'Using ', trim(coord), ' coordinates'
    if (mype.eq.0 .and. space_is_diagonal) write(*,*) 'Space is diagonal'
    if (mype.eq.0 .and. .not.space_is_diagonal) write(*,*) &
       'Space is non-diagonal'
    if (mype.eq.0) write(*,*) 'Sqrt(gamma) is analytic:', &
       sqrtgamma_is_analytic
    ! The metric will only be filled as a ^NC x ^NC matrix
    ! For the induced metrics (surfaces), it is assumed that the missing entries
    ! are diagonal 1.  If this is not what you have in mind, you should run with
    ! 3 components.  
    if (mype.eq.0.and.3.ne.3) write(*,*) &
       'Assuming left-out components are identity'

    ! print the shape of the metric to screen:
    if (mype.eq.0) then
       print*,'-----------------------------------------'
       write(*,*) 'These metric elements are zero:'
       do i=0,3
          write(*,*) g_is_zero(i,:)
       end do
       write(*,*) 'Number of non-zero elements:', nnonzero_metric
    end if
    
    ! print the shape of the contravariant shift to screen:
    if (mype.eq.0) then
       print*,'-----------------------------------------'
       write(*,*) 'These shift-components are zero:'
       write(*,*) beta_is_zero(:)
       write(*,*) 'Number of non-zero elements:', nnonzero_beta
    end if

    ! print the shape of the metric derivatives to screen:
    if (mype.eq.0) then
       print*,'-----------------------------------------'
       write(*,*) 'These components of the metric derivative DgijDk are zero:'
       do k=1,3
       print*,'-----------'
       write(*,*) 'k=',k
          do i=1,3
             write(*,*) dgdk_is_zero(i,:,k)
          end do
       end do
       print*,'-----------'
       write(*,*) 'Number of non-zero elements:', nnonzero_dgdk
    end if

    ! print the shape of the shift derivatives:
    if (mype.eq.0) then
       print*,'-----------------------------------------'
       write(*,*) 'These elements of shift derivatives dbetaidj are zero:'
       do i=1,3
          write(*,*) dbetaidj_is_zero(i,:)
       end do
       write(*,*) 'Number of non-zero elements:', nnonzero_dbetaidj
    end if

    if (mype.eq.0) then
       print*,'-----------------------------------------'
       write(*,*) 'These elements of the lapse derivative dalphadj are zero:'
       write(*,*) dalphadj_is_zero(:)
       write(*,*) 'Number of non-zero elements:', nnonzero_dalphadj
    end if

    if (mype.eq.0) print*,'-----------------------------------------------------------------------------'
    
  end subroutine init_metric
  !=============================================================================
  subroutine LU(m)
    ! LU-decomposes the spatial metric and calculates the inverse
    ! much room for optimization. 
    use mod_lu
    ! if this breaks for you, try commenting out and use the old gfortran hack below.
    use, intrinsic :: IEEE_ARITHMETIC, only: ieee_value, ieee_quiet_nan
    !
    use mod_amrvacdef

    type(metric), pointer                    :: m
    ! .. local ..
    integer, dimension(3)                  :: indx
    integer                                  :: inonzero, i, j, code, d, ix1,&
       ix2,ix3
    double precision, dimension(1:3,1:3) :: a
    double precision, dimension(1:3)       :: b
    ! old gfortran (< v5) needed this hack but it breaks on gfortran 10:    
    !    double precision, PARAMETER :: D_QNAN = TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_8)
    double precision                         :: D_QNAN
    D_QNAN = ieee_value(1.,ieee_quiet_nan)
    !-----------------------------------------------------------------------------

    if (space_is_diagonal) then 

       do i=1,3
          m%gammainv(i,i)%elem(m%ixGmin1:m%ixGmax1,m%ixGmin2:m%ixGmax2,&
             m%ixGmin3:m%ixGmax3) = 1.0d0/m%g(i,i)%elem(m%ixGmin1:m%ixGmax1,&
             m%ixGmin2:m%ixGmax2,m%ixGmin3:m%ixGmax3)
       end do

    else ! space is diagonal

       
       ! --------------------------------------------------
       ! Initialize as quiet NaN
       ! --------------------------------------------------
       do j=1,3
          do i=1,j
             m%gammainv(i,j)%elem(m%ixGmin1:m%ixGmax1,m%ixGmin2:m%ixGmax2,&
                m%ixGmin3:m%ixGmax3) = D_QNAN
          end do
       end do
       

       ! --------------------------------------------------
       ! re-arrange elements and do the LU by using LUDCMP
       ! from Numerical recipes
       ! --------------------------------------------------
        do ix3=m%ixGmin3, m%ixGmax3
         do ix2=m%ixGmin2, m%ixGmax2
         do ix1=m%ixGmin1, m%ixGmax1
       a(:,:) = 0.0d0
       do inonzero=1,m%nnonzero
          i = m%nonzero(inonzero)%i
          j = m%nonzero(inonzero)%j
          if (i==0 .or. j==0) cycle

          a(i,j) = m%nonzero(inonzero)%elem(ix1,ix2,ix3)

       end do

       call LUDCMP(a,3,indx,d,code)

       ! --------------------------------------------------
       ! The metric becomes singular on axis (return code=1).
       ! Don't stop in these cases, but keep the NaN values.
       ! Fluxes through the axis are identically set to zero to catch this.
       ! (e.g. in subroutine tvdlf)
       ! These Nan's should thus not filter through to the solution.
       ! --------------------------------------------------
       if (code /= 1) then

          do j=1,3
             ! Directly calculate inverse metric:
             b(1)=kr(j,1);b(2)=kr(j,2);b(3)=kr(j,3)
             call LUBKSB(a,3,indx,b)
             do i=1,j
                m%gammainv(i,j)%elem(ix1,ix2,ix3) = b(i)
             end do
          end do

       end if ! metric non-singular
       
        end do
         end do
         end do

    end if ! space is diagonal

  end subroutine LU
  !=============================================================================
  subroutine raise3(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,m,d,u)
    ! Takes a covariant three-vector with ndir components d (stands for down)
    ! and returns the contravariant one u with ndir components (u stands for up)
    use mod_amrvacdef

    integer, intent(in)                      :: ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3
    double precision, intent(in)             :: d(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndir)
    double precision, intent(out)            :: u(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndir)
    type(metric), pointer                    :: m
    ! .. local ..
    integer                                  :: i, j
    !-----------------------------------------------------------------------------
  
    if (space_is_diagonal) then
       
       do i=1,3
          u(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,i) &
             = m%gammainv(i,i)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * d(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,i)
       end do
       
    else

       u(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:) = 0.0d0
       do j=1,3
          do i=1,3
             u(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,i) &
                = u(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                i) + m%gammainv(i,j)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)*d(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3,j)
          end do
       end do

    end if ! space is diagonal

  end subroutine raise3
  !=============================================================================
  subroutine raise3_point(x1,x2,x3,d,u)
    ! Takes a covariant three-vector at coordinate x^D with ndir components d (stands for down)
    ! and returns the contravariant one u with ndir components (u stands for up)
    ! This is for pointwise-values!

    double precision, dimension(1:3), intent(in)    :: d
    double precision, dimension(1:3), intent(out)   :: u
    double precision, intent(in)                      :: x1,x2,x3
    ! .. local ..
    double precision, dimension(1:3,1:3) :: gammainv
    integer                                  :: i, j
    !-----------------------------------------------------------------------------

    ! Component-wise. Could be more efficient.
    do j=1,3
       do i=1,j
          call get_gammainv_component(i,j,x1,x2,x3,gammainv(i,j))
          if (i .ne. j) gammainv(j,i) = gammainv(i,j)
       end do
    end do

    ! raise:

    if (space_is_diagonal) then
       
       do i=1,3
          u(i) = gammainv(i,i) * d(i)
       end do

    else

       u = 0.0d0
       do i=1,3
          do j=1,3
             u(i) = u(i) + gammainv(i,j) * d(j)
          end do
       end do

    end if

  end subroutine raise3_point
  !=============================================================================
  subroutine lower3(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,m,u,d)
    ! Takes a contravariant three-vector with ndir components u (stands for up)
    ! and returns the covariant one d with ndir components (d stands for down)
    use mod_amrvacdef

    integer, intent(in)                      :: ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3
    double precision, intent(in)             :: u(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndir)
    double precision, intent(out)            :: d(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndir)
    type(metric), pointer                    :: m
    ! .. local ..
    integer                                  :: inonzero, i, j
    !-----------------------------------------------------------------------------

    d(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:) = 0.0d0
    
    do inonzero=1,m%nnonzero
       i = m%nonzero(inonzero)%i
       j = m%nonzero(inonzero)%j
       if (i==0 .or. j==0) cycle

       d(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,i) &
          = d(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          i) + m%nonzero(inonzero)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) * u(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,j)
       
    end do
    
  end subroutine lower3
    !=============================================================================
  subroutine lower3_point(x1,x2,x3,u,d)
    ! Takes a contravariant three-vector at coordinate x^D with ndir components u (stands for up)
    ! and returns the covariant one d with ndir components (d stands for down)
    ! This is for pointwise-values!

    double precision, dimension(1:3), intent(in)    :: u
    double precision, dimension(1:3), intent(out)   :: d
    double precision, intent(in)                      :: x1,x2,x3
    ! .. local ..
    double precision, dimension(1:3,1:3) :: g
    integer                                  :: i, j
    !-----------------------------------------------------------------------------

    ! Component-wise. Could be more efficient.
    do j=1,3
       do i=1,j
          call get_g_component(i,j,x1,x2,x3,g(i,j))
          if (i .ne. j) g(j,i) = g(i,j)
       end do
    end do

    ! lower:

    if (space_is_diagonal) then
       
       do i=1,3
          d(i) = g(i,i) * u(i)
       end do

    else

       d = 0.0d0
       do i=1,3
          do j=1,3
             d(i) = d(i) + g(i,j) * u(j)
          end do
       end do

    end if

  end subroutine lower3_point
  !=============================================================================
  subroutine square3u(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,m,u,u2)
    ! Optimised square b^2=b^i b_i for contravariant input
    ! three-vectors u

    use mod_amrvacdef

    integer, intent(in)                      :: ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3
    double precision, intent(in)             :: u(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndir)
    double precision, intent(out)            :: u2(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)
    type(metric), pointer                    :: m
    ! .. local ..
    integer                                  :: i
    !-----------------------------------------------------------------------------

    u2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = zero
    
    ! Diagonal terms:
    do i = 1, 3
       u2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
          = u2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + &
          u(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,i)&
          **2 * m%g(i,i)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end do
    
    ! Offdiagonal entries:
    
    if (.not. g_is_zero(1,2)) & 
         u2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
            = u2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + 2.0d0 * &
            m%g(1,2)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) * &
            u(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            1)*u(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)
   
    if (.not. g_is_zero(1,3)) & 
         u2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
            = u2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + 2.0d0 * &
            m%g(1,3)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) * &
            u(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            1)*u(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)
    if (.not. g_is_zero(2,3)) & 
         u2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
            = u2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + 2.0d0 * &
            m%g(2,3)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) * &
            u(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            2)*u(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)
   
    
  end subroutine square3u
  !=============================================================================
    subroutine acrossbU(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,m,idir,a,b,res)
      ! Do the crossproduct of two covariant vectors
      ! Result is contravariant
      ! res^idir = eta^{idir,j,k} 1/sqrtgamma a_j b_k
    use mod_amrvacdef

    integer, intent(in)                                        :: ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,&
       ixOmax1,ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)       :: a,b
    integer         , intent(in)                               :: idir
    type(metric), pointer                                      :: m
    double precision, intent(inout)                            :: &
       res(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    ! .. local ..
    integer                                                    :: j,k
    !-----------------------------------------------------------------------------

    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = zero
    do j=1,3
       do k=1,3
          if (j .eq. k .or. j .eq. idir .or. k .eq. idir) cycle
             res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                = res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + &
                1.0d0/m%sqrtgamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)*lvc(idir,j,k)*a(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3,j)*b(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3,k)
       end do
    end do

  end subroutine acrossbU
  !=============================================================================
  subroutine lower4(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,m,u,d)
    ! Takes a contravariant four-vector u (stands for up) 
    ! and returns the covariant one d (stands for down)
    ! Note: index of vector starts at zero!
    use mod_amrvacdef

    integer, intent(in)                      :: ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3
    double precision, intent(in)             :: u(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,0:ndir)
    double precision, intent(out)            :: d(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,0:ndir)
    type(metric), pointer                    :: m
    ! .. local ..
    integer                                  :: inonzero, i, j
    !-----------------------------------------------------------------------------

    d(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:) = 0.0d0

    do inonzero=1,nnonzero_metric
       i = m%nonzero(inonzero)%i
       j = m%nonzero(inonzero)%j

       d(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,i) &
          = d(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          i) + m%nonzero(inonzero)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) * u(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,j)
       
    end do
    
  end subroutine lower4
  !=============================================================================
  subroutine allocate_metric(m,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
     ixGmax3,need_derivs)
    use mod_amrvacdef

    type(metric), pointer, intent(inout)     :: m
    integer, intent(in)                      :: ixGmin1,ixGmin2,ixGmin3,&
       ixGmax1,ixGmax2,ixGmax3
    logical, optional                        :: need_derivs
    ! .. local ..
    integer                                  :: inonzero, i,j,k
    logical                                  :: alloc_derivs
    !-----------------------------------------------------------------------------
    if (.not.present(need_derivs)) then
       alloc_derivs = .true.
    else
       alloc_derivs = need_derivs
    end if

    
    allocate(m)
    allocate(m%nonzero(1:nnonzero_metric))
    allocate(m%nonzeroBeta(1:nnonzero_beta))
    allocate(m%g(0:3,0:3))
    allocate(m%beta(1:3))
    allocate(m%betaD(1:3))
    allocate(m%sqrtgamma(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3))
    allocate(m%alpha(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3))
    ! Derivatives:
    allocate(m%nonzeroDgDk(1:nnonzero_dgdk),m%DgDk(1:3,1:3,1:3))
    allocate(m%nonzeroDalphaDj(1:nnonzero_dalphadj),m%DalphaDj(1:3))
    allocate(m%nonzeroDbetaiDj(1:nnonzero_dbetaidj),m%DbetaiDj(1:3,1:3))
    ! Inverse metric:
    allocate(m%gammainv(1:3,1:3))
    
    ! Store shape of arrays:
    m%ixGmin1=ixGmin1;m%ixGmin2=ixGmin2;m%ixGmin3=ixGmin3;m%ixGmax1=ixGmax1
    m%ixGmax2=ixGmax2;m%ixGmax3=ixGmax3;

    ! we also always allocate one element of zeros
    allocate(m%zero(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3))
    m%zero = 0.0d0

    ! Inverse metric
       do i=1,3
          do j=1,3
             m%gammainv(i,j)%i = i
             m%gammainv(i,j)%j = j
             m%gammainv(i,j)%k = -1
             if (space_is_diagonal) then
                if (i .eq. j) then
                   allocate(m%gammainv(i,i)%elem(ixGmin1:ixGmax1,&
                      ixGmin2:ixGmax2,ixGmin3:ixGmax3))
                else
                   m%gammainv(i,j)%elem => m%zero
                end if
             else
                if (j .ge. i) then 
                   allocate(m%gammainv(i,j)%elem(ixGmin1:ixGmax1,&
                      ixGmin2:ixGmax2,ixGmin3:ixGmax3))
                else
                   m%gammainv(i,j)%elem => m%gammainv(j,i)%elem
                end if
             end if
          end do
       end do

       
    ! Allocate contravariant shifts m%beta
    inonzero = 0
    do j=1,3
       m%beta(j)%i = j
       m%beta(j)%j = j
       m%beta(j)%k = j
       if (beta_is_zero(j)) then 
          m%beta(j)%elem => m%zero
       else
          inonzero = inonzero+1
          allocate(m%beta(j)%elem(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
             ixGmin3:ixGmax3))
          call associateIndexlist(m%nonzeroBeta(inonzero),m%beta(j))
       end if
    end do
    m%nnonzeroBeta = inonzero

    ! Covariant shifts point to metric:
    do j = 1,3
          nullify(m%betaD(j)%elem)
          m%betaD(j)%i = j
          m%betaD(j)%j = j
          m%betaD(j)%k = j
    end do

    ! Point them all to zero if there is at least one
    ! zero-element:
    if (any(g_is_zero(:,:))) then
       do i = 0,3
          do j = 0,3
             m%g(i,j)%i = i
             m%g(i,j)%j = j
             m%g(i,j)%k = -1
             m%g(i,j)%elem => m%zero
          end do
       end do
    else
       do i = 0,3
          do j = 0,3
             m%g(i,j)%i = i
             m%g(i,j)%j = j
             m%g(i,j)%k = -1
             nullify(m%g(i,j)%elem)
          end do
       end do
    end if

    !------------------------------------------------
    ! Initialize indices and allocate the metric derivatives:
    !------------------------------------------------
    inonzero = 0
    do k=1,3
       do i=1,3
          do j=1,3
             m%DgDk(i,j,k)%i = i
             m%DgDk(i,j,k)%j = j
             m%DgDk(i,j,k)%k = k
             nullify(m%DgDk(i,j,k)%elem)
             if (.not.dgdk_is_zero(i,j,k)) then
                inonzero = inonzero + 1
                m%nonzeroDgDk(inonzero)%i = i
                m%nonzeroDgDk(inonzero)%j = j
                m%nonzeroDgDk(inonzero)%k = k
                ! We need only to allocate the upper triangle with j>=i
                ! Link the other components to the transposed DgDk(j,i,k) structure
                if (j .ge. i) then
                   if (alloc_derivs) allocate(m%nonzeroDgDk(inonzero)%elem&
                      (ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3))
                   m%DgDk(i,j,k)%elem => m%nonzeroDgDk(inonzero)%elem
                else
                   m%nonzeroDgDk(inonzero)%elem => m%DgDk(j,i,k)%elem 
                   m%DgDk(i,j,k)%elem => m%DgDk(j,i,k)%elem
                end if
             else
                m%DgDk(i,j,k)%elem => m%zero
             end if
          end do
       end do
    end do
    m%nnonzeroDgDk = inonzero

    !------------------------------------------------
    ! Initialize indices and allocate the shift derivatives:
    !------------------------------------------------
    inonzero = 0
    do j=1,3
       do i=1,3
          m%DbetaiDj(i,j)%i = i
          m%DbetaiDj(i,j)%j = j
          m%DbetaiDj(i,j)%k = -1
          nullify(m%DbetaiDj(i,j)%elem)
          if (.not. dbetaidj_is_zero(i,j)) then
             inonzero = inonzero + 1
             m%nonzeroDbetaiDj(inonzero)%i = i
             m%nonzeroDbetaiDj(inonzero)%j = j
             m%nonzeroDbetaiDj(inonzero)%k = -1
             if (alloc_derivs) allocate(m%nonzeroDbetaiDj(inonzero)%elem&
                (ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3))
             m%DbetaiDj(i,j)%elem => m%nonzeroDbetaiDj(inonzero)%elem
          else
             m%DbetaiDj(i,j)%elem => m%zero
          end if
       end do
    end do
    m%nnonzeroDbetaiDj = inonzero

    !------------------------------------------------
    ! Initialize indices and allocate the lapse derivatives:
    !------------------------------------------------
    inonzero = 0
    do j=1,3
       m%DalphaDj(j)%i = j
       m%DalphaDj(j)%j = j
       m%DalphaDj(j)%k = j
       nullify(m%DalphaDj(j)%elem)
       if (.not. dalphadj_is_zero(j)) then
             inonzero = inonzero + 1
             m%nonzeroDalphaDj(inonzero)%i = j
             m%nonzeroDalphaDj(inonzero)%j = j
             m%nonzeroDalphaDj(inonzero)%k = j
             if (alloc_derivs) allocate(m%nonzeroDalphaDj(inonzero)%elem&
                (ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3))
             m%DalphaDj(j)%elem => m%nonzeroDalphaDj(inonzero)%elem
          else
             m%DalphaDj(j)%elem => m%zero
          end if
    end do
    m%nnonzeroDalphaDj = inonzero

    
    ! ==============================
    ! Allocate and link the non-zero metric components
    ! ==============================
    
    inonzero=0
    if (.not.g_is_zero(0,0)) then
       inonzero = inonzero + 1
       allocate(m%nonzero(inonzero)%elem(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
          ixGmin3:ixGmax3))
       m%nonzero(inonzero)%i = 0
       m%nonzero(inonzero)%j = 0
       m%nonzero(inonzero)%k = -1
       m%g(0,0)%elem => m%nonzero(inonzero)%elem
    end if
    
    
    if (.not.g_is_zero(0,1)) then
       inonzero = inonzero + 1
       allocate(m%nonzero(inonzero)%elem(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
          ixGmin3:ixGmax3))
       m%nonzero(inonzero)%i = 0
       m%nonzero(inonzero)%j = 1
       m%nonzero(inonzero)%k = -1
       m%g(0,1)%elem => m%nonzero(inonzero)%elem
       ! Lower triangle:
       inonzero = inonzero + 1       
       m%nonzero(inonzero)%i = 1
       m%nonzero(inonzero)%j = 0
       m%nonzero(inonzero)%k = -1
       m%nonzero(inonzero)%elem => m%g(0,1)%elem
       m%g(1,0)%elem => m%g(0,1)%elem
    end if
    m%betaD(1)%elem => m%g(0,1)%elem
    
    
    if (.not.g_is_zero(0,2)) then
       inonzero = inonzero + 1
       allocate(m%nonzero(inonzero)%elem(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
          ixGmin3:ixGmax3))
       m%nonzero(inonzero)%i = 0
       m%nonzero(inonzero)%j = 2
       m%nonzero(inonzero)%k = -1
       m%g(0,2)%elem => m%nonzero(inonzero)%elem
       ! Lower triangle:
       inonzero = inonzero + 1       
       m%nonzero(inonzero)%i = 2
       m%nonzero(inonzero)%j = 0
       m%nonzero(inonzero)%k = -1
       m%nonzero(inonzero)%elem => m%g(0,2)%elem
       m%g(2,0)%elem => m%g(0,2)%elem
    end if
    m%betaD(2)%elem => m%g(0,2)%elem
    
    
    if (.not.g_is_zero(0,3)) then
       inonzero = inonzero + 1
       allocate(m%nonzero(inonzero)%elem(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
          ixGmin3:ixGmax3))
       m%nonzero(inonzero)%i = 0
       m%nonzero(inonzero)%j = 3
       m%nonzero(inonzero)%k = -1
       m%g(0,3)%elem => m%nonzero(inonzero)%elem
       ! Lower triangle:
       inonzero = inonzero + 1       
       m%nonzero(inonzero)%i = 3
       m%nonzero(inonzero)%j = 0
       m%nonzero(inonzero)%k = -1
       m%nonzero(inonzero)%elem => m%g(0,3)%elem
       m%g(3,0)%elem => m%g(0,3)%elem
    end if
    m%betaD(3)%elem => m%g(0,3)%elem
    

    
    
    if (.not.g_is_zero(1,1)) then
       inonzero = inonzero + 1
       allocate(m%nonzero(inonzero)%elem(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
          ixGmin3:ixGmax3))
       m%nonzero(inonzero)%i = 1
       m%nonzero(inonzero)%j = 1
       m%nonzero(inonzero)%k = -1
       m%g(1,1)%elem => m%nonzero(inonzero)%elem
    end if
    
    
    if (.not.g_is_zero(1,2)) then
       inonzero = inonzero + 1
       allocate(m%nonzero(inonzero)%elem(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
          ixGmin3:ixGmax3))
       m%nonzero(inonzero)%i = 1
       m%nonzero(inonzero)%j = 2
       m%nonzero(inonzero)%k = -1
       m%g(1,2)%elem => m%nonzero(inonzero)%elem
    end if
    
    
    if (.not.g_is_zero(1,3)) then
       inonzero = inonzero + 1
       allocate(m%nonzero(inonzero)%elem(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
          ixGmin3:ixGmax3))
       m%nonzero(inonzero)%i = 1
       m%nonzero(inonzero)%j = 3
       m%nonzero(inonzero)%k = -1
       m%g(1,3)%elem => m%nonzero(inonzero)%elem
    end if
    
       
    
    
    
    if (.not.g_is_zero(2,2)) then
       inonzero = inonzero + 1
       allocate(m%nonzero(inonzero)%elem(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
          ixGmin3:ixGmax3))
       m%nonzero(inonzero)%i = 2
       m%nonzero(inonzero)%j = 2
       m%nonzero(inonzero)%k = -1
       m%g(2,2)%elem => m%nonzero(inonzero)%elem
    end if
    if (.not.g_is_zero(2,1)) then
       inonzero = inonzero + 1
       m%nonzero(inonzero)%i = 2
       m%nonzero(inonzero)%j = 1
       m%nonzero(inonzero)%k = -1
       m%nonzero(inonzero)%elem => m%g(1,2)%elem
       m%g(2,1)%elem => m%g(1,2)%elem
    end if
    if (.not.g_is_zero(2,3)) then
       inonzero = inonzero + 1
       allocate(m%nonzero(inonzero)%elem(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
          ixGmin3:ixGmax3))
       m%nonzero(inonzero)%i = 2
       m%nonzero(inonzero)%j = 3
       m%nonzero(inonzero)%k = -1
       m%g(2,3)%elem => m%nonzero(inonzero)%elem
       ! Lower triangle
       inonzero = inonzero + 1
       m%nonzero(inonzero)%i = 3
       m%nonzero(inonzero)%j = 2
       m%nonzero(inonzero)%k = -1
       m%nonzero(inonzero)%elem => m%g(2,3)%elem
       m%g(3,2)%elem => m%g(2,3)%elem
    end if
    
    if (.not.g_is_zero(3,3)) then
       inonzero = inonzero + 1
       allocate(m%nonzero(inonzero)%elem(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
          ixGmin3:ixGmax3))
       m%nonzero(inonzero)%i = 3
       m%nonzero(inonzero)%j = 3
       m%nonzero(inonzero)%k = -1
       m%g(3,3)%elem => m%nonzero(inonzero)%elem
    end if

    if (.not.g_is_zero(3,1)) then
       inonzero = inonzero + 1
       m%nonzero(inonzero)%i = 3
       m%nonzero(inonzero)%j = 1
       m%nonzero(inonzero)%k = -1
       m%nonzero(inonzero)%elem => m%g(1,3)%elem
       m%g(3,1)%elem => m%g(1,3)%elem
    end if
    

    m%nnonzero  = inonzero

    if (inonzero .ne. nnonzero_metric) call mpistop&
       ('allocate_metric: wrong number of non-zero elements')
    
  end subroutine allocate_metric
  !=============================================================================
  subroutine associateIndexlist(src,tgt)
    use mod_amrvacdef
    
    type(indexlist), intent(inout)                    :: src,tgt
    !-----------------------------------------------------------------------------

    src%elem => tgt%elem
    src%i = tgt%i
    src%j = tgt%j
    src%k = tgt%k

  end subroutine associateIndexlist
  !=============================================================================
  subroutine deallocate_metric(m)
    use mod_amrvacdef

    type(metric), pointer, intent(inout)   :: m
    integer                                :: i,j, inonzero
    !-----------------------------------------------------------------------------
    
    !----------------------------------------
    ! Free metric (and covariant shift) memory:
    !----------------------------------------
    do inonzero=1,m%nnonzero
       ! Only deallocate upper triangle, the others are just links
       if (m%nonzero(inonzero)%j .ge. m%nonzero(inonzero)%i) then
          deallocate(m%nonzero(inonzero)%elem)
       end if
       nullify(m%nonzero(inonzero)%elem)
    end do
    deallocate(m%nonzero); nullify(m%nonzero)

    ! Nullify covariant shifts (elements of four-metric), lived in m%nonzero
    deallocate(m%betaD); nullify(m%betaD)
    
    ! Nullify four-metric, lived in m%nonzero
    deallocate(m%g); nullify(m%g)
    !----------------------------------------

    
    !----------------------------------------
    ! Metric derivatives:
    !----------------------------------------
    do inonzero=1,m%nnonzeroDgDk
       ! Deallocate only upper triangle in i,j, the others are just links
       if (m%nonzeroDgDk(inonzero)%j .ge. m%nonzeroDgDk(inonzero)%i &
          .and. associated(m%nonzeroDgDk(inonzero)%elem)) then
          deallocate(m%nonzeroDgDk(inonzero)%elem)
       end if
       nullify(m%nonzeroDgDk(inonzero)%elem)
    end do
    deallocate(m%nonzeroDgDk); nullify(m%nonzeroDgDk)
    deallocate(m%DgDk); nullify(m%DgDk)
    !----------------------------------------

    
    !----------------------------------------
    ! Lapse derivatives:
    !----------------------------------------
    do inonzero=1,m%nnonzeroDalphaDj
       if (associated(m%nonzeroDalphaDj(inonzero)%elem)) deallocate&
          (m%nonzeroDalphaDj(inonzero)%elem)
       nullify(m%nonzeroDalphaDj(inonzero)%elem)
    end do
    deallocate(m%nonzeroDalphaDj); nullify(m%nonzeroDalphaDj)
    deallocate(m%DalphaDj); nullify(m%DalphaDj)
    !----------------------------------------

    
    !----------------------------------------
    ! Shift derivatives:
    !----------------------------------------
    do inonzero=1,m%nnonzeroDbetaiDj
       if (associated(m%nonzeroDbetaiDj(inonzero)%elem)) deallocate&
          (m%nonzeroDbetaiDj(inonzero)%elem)
       nullify(m%nonzeroDbetaiDj(inonzero)%elem)
    end do
    deallocate(m%nonzeroDbetaiDj); nullify(m%nonzeroDbetaiDj)
    deallocate(m%DbetaiDj); nullify(m%DbetaiDj)
    !----------------------------------------
    
    
    !----------------------------------------
    ! deallocate or just nullify contravariant shift:
    !----------------------------------------
    do i=1,3
       if (.not.beta_is_zero(i)) then 
          deallocate(m%beta(i)%elem); nullify(m%beta(i)%elem)
       else
          nullify(m%beta(i)%elem)
       end if
    end do
    deallocate(m%nonzeroBeta); nullify(m%nonzeroBeta)
    deallocate(m%beta); nullify(m%beta)
    !----------------------------------------

    
    !----------------------------------------
    ! Lapse and shift:
    !----------------------------------------
    deallocate(m%alpha,m%sqrtgamma); nullify(m%alpha,m%sqrtgamma)
    !----------------------------------------

    
    !----------------------------------------
    ! inverse Metric:
    !----------------------------------------
    if (space_is_diagonal) then
       do i=1,3
          deallocate(m%gammainv(i,i)%elem); nullify(m%gammainv(i,i)%elem)
       end do
    else
       do i=1,3
          do j=1,3
             ! Deallocate only upper triangle in i,j, the others are just links
             if (j .ge. i) then
                deallocate(m%gammainv(i,j)%elem)
             end if
             nullify(m%gammainv(i,j)%elem)
          end do
       end do
    end if
    deallocate(m%gammainv); nullify(m%gammainv)
    !----------------------------------------


    deallocate(m%zero); nullify(m%zero)

    deallocate(m); nullify(m)

  end subroutine deallocate_metric
  !=============================================================================
  subroutine fill_metric(m,ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,&
     ixGextmax2,ixGextmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,&
     need_derivs)
    use mod_amrvacdef

    type(metric), pointer                                   :: m
    integer, intent(in)                                     :: ixGextmin1,&
       ixGextmin2,ixGextmin3,ixGextmax1,ixGextmax2,ixGextmax3, ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, dimension(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,1:ndim), intent(in):: x
    logical, optional                                       :: need_derivs
    ! .. local ..
    integer                                                 :: ix1,ix2,ix3, i,&
        j, k, inonzero
    double precision, dimension(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,1:ndir)            :: betaD !covariant shift
    double precision, dimension(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,1:ndir)            :: betaU !contravariant sift
    double precision, dimension(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3)                   :: beta2
    double precision                                        :: dummy
    double precision, dimension(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,1:3,1:3,1:3) :: dgammainvdk
    double precision, dimension(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,0:3,1:3)       :: dbetaDownjDk
    logical                                                 :: fill_derivs
    !-----------------------------------------------------------------------------
    if (.not.present(need_derivs)) then
       fill_derivs = .true.
    else
       fill_derivs = need_derivs
    end if
    
    !--------------------------------------------------
    if (.not. init_from_g4) then
    !--------------------------------------------------


       ! set the lapse:
        do ix3=ixOmin3, ixOmax3
         do ix2=ixOmin2, ixOmax2
         do ix1=ixOmin1, ixOmax1
       call get_alpha(x(ix1,ix2,ix3,1),x(ix1,ix2,ix3,2),x(ix1,ix2,ix3,3),&
          m%alpha(ix1,ix2,ix3))
        end do 
         end do 
         end do 

       ! set the shift (contravariant):
       
            if (.not. beta_is_zero(1)) then
        do ix3=ixOmin3, ixOmax3
        do ix2=ixOmin2, ixOmax2
        do ix1=ixOmin1, ixOmax1
       call get_beta(1,x(ix1,ix2,ix3,1),x(ix1,ix2,ix3,2),x(ix1,ix2,ix3,3),&
          m%beta(1)%elem(ix1,ix2,ix3))
        end do 
        end do 
        end do 
       end if
       
       
            if (.not. beta_is_zero(2)) then
        do ix3=ixOmin3, ixOmax3
        do ix2=ixOmin2, ixOmax2
        do ix1=ixOmin1, ixOmax1
       call get_beta(2,x(ix1,ix2,ix3,1),x(ix1,ix2,ix3,2),x(ix1,ix2,ix3,3),&
          m%beta(2)%elem(ix1,ix2,ix3))
        end do 
        end do 
        end do 
       end if
       
       
            if (.not. beta_is_zero(3)) then
        do ix3=ixOmin3, ixOmax3
        do ix2=ixOmin2, ixOmax2
        do ix1=ixOmin1, ixOmax1
       call get_beta(3,x(ix1,ix2,ix3,1),x(ix1,ix2,ix3,2),x(ix1,ix2,ix3,3),&
          m%beta(3)%elem(ix1,ix2,ix3))
        end do 
        end do 
        end do 
       end if
       
       
       ! set the spatial part:
       do inonzero=1,nnonzero_metric
          i = m%nonzero(inonzero)%i
          j = m%nonzero(inonzero)%j

          ! Just spatial and only upper triangle (the others are already linking here)
          if (i == 0 .or. j == 0 .or. i .gt. j) cycle

           do ix3=ixOmin3, ixOmax3
            do ix2=ixOmin2, ixOmax2
            do ix1=ixOmin1, ixOmax1
          call get_g_component(i,j,x(ix1,ix2,ix3,1),x(ix1,ix2,ix3,2),x(ix1,&
             ix2,ix3,3),m%nonzero(inonzero)%elem(ix1,ix2,ix3))
           end do 
            end do 
            end do 
       end do

       ! set the space-time part (covariant shifts):
        betaU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) &
           = m%beta(1)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
         betaU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) &
            = m%beta(2)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
         betaU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3) &
            = m%beta(3)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
       call lower3(ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,ixGextmax2,&
          ixGextmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,m,betaU,&
          betaD)
       do i=1,3
          if (.not. g_is_zero(0,i)) m%g(0,i)%elem(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3) = betaD(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,i)
       end do
       
       ! set the 00-component from the shift and lapse
       beta2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
          =  betaD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          1)*betaU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          1) + betaD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          2)*betaU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          2) + betaD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          3)*betaU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3) 
       m%g(0,0)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
          = - m%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
          **2 + beta2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

       ! set the square root of the spatial determinant
        do ix3=ixOmin3, ixOmax3
         do ix2=ixOmin2, ixOmax2
         do ix1=ixOmin1, ixOmax1
       call get_sqrtgamma(x(ix1,ix2,ix3,1),x(ix1,ix2,ix3,2),x(ix1,ix2,ix3,3),&
          m%sqrtgamma(ix1,ix2,ix3))
        end do 
         end do 
         end do 

       ! set the inverse of the three-metric
       call LU(m)

       !------------------------------------------------------------
       ! Fill the derivative data:
       !------------------------------------------------------------
       if (fill_derivs) then

          ! set the lapse:
          do inonzero=1, m%nnonzeroDalphaDj
             j = m%nonzeroDalphaDj(inonzero)%j
              do ix3=ixOmin3, ixOmax3
               do ix2=ixOmin2, ixOmax2
               do ix1=ixOmin1, ixOmax1
             call get_alpha(x(ix1,ix2,ix3,1),x(ix1,ix2,ix3,2),x(ix1,ix2,ix3,&
                3),dummy,dalphadj=m%nonzeroDalphaDj(inonzero)%elem(ix1,ix2,&
                ix3),jdir=j)
              end do 
               end do 
               end do 
          end do

          ! set the shift (contravariant):
          do inonzero=1,m%nnonzeroDbetaiDj
             i = m%nonzeroDbetaiDj(inonzero)%i
             j = m%nonzeroDbetaiDj(inonzero)%j
              do ix3=ixOmin3, ixOmax3
               do ix2=ixOmin2, ixOmax2
               do ix1=ixOmin1, ixOmax1
             call get_beta(i,x(ix1,ix2,ix3,1),x(ix1,ix2,ix3,2),x(ix1,ix2,ix3,&
                3),dummy,dbetaidj=m%nonzeroDbetaiDj(inonzero)%elem(ix1,ix2,&
                ix3),jdir=j)
              end do 
               end do 
               end do 
          end do

          ! set the spatial part:
          do inonzero=1,m%nnonzeroDgDk
             i = m%nonzeroDgDk(inonzero)%i
             j = m%nonzeroDgDk(inonzero)%j
             k = m%nonzeroDgDk(inonzero)%k

             ! Just spatial and only upper triangle (the others are already linking here)
             if (i == 0 .or. j == 0 .or. i .gt. j) cycle

              do ix3=ixOmin3, ixOmax3
               do ix2=ixOmin2, ixOmax2
               do ix1=ixOmin1, ixOmax1
             call get_g_component(i,j,x(ix1,ix2,ix3,1),x(ix1,ix2,ix3,2),x(ix1,&
                ix2,ix3,3),dummy,dgdk=m%nonzeroDgDk(inonzero)%elem(ix1,ix2,&
                ix3),kdir=k)
              end do 
               end do 
               end do 
          end do

       end if

    !--------------------------------------------------
    else ! init_from_g4
    !--------------------------------------------------

       ! set the full four metric:
       do inonzero=1,nnonzero_metric
          i = m%nonzero(inonzero)%i
          j = m%nonzero(inonzero)%j

          ! Only upper triangle (the others are already linking here)
          if (i .gt. j) cycle

           do ix3=ixOmin3, ixOmax3
            do ix2=ixOmin2, ixOmax2
            do ix1=ixOmin1, ixOmax1
          call get_g_component(i,j,x(ix1,ix2,ix3,1),x(ix1,ix2,ix3,2),x(ix1,&
             ix2,ix3,3),m%nonzero(inonzero)%elem(ix1,ix2,ix3))

           end do 
            end do 
            end do 
       end do
       
       ! set the inverse of the three-metric
       call LU(m)

       ! raise the shift vector and fill beta:
       do i=1,3
          betaD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,i) &
             = m%g(0,i)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
       end do
       call raise3(ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,ixGextmax2,&
          ixGextmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,m,betaD,&
          betaU)

        m%beta(1)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
           = betaU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) 
         m%beta(2)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
            = betaU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) 
         m%beta(3)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
            = betaU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3) 

       ! Obtain the lapse from shift and g00:
       beta2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
          =  betaD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          1)*betaU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          1) + betaD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          2)*betaU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          2) + betaD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          3)*betaU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3) 

       ! Freeze regions where the lapse becomes negative
       ! (such region should be protected by a horizon)
       m%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
          = sqrt(max(beta2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)-m%g(0,0)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3),0.0d0))

       ! set the square root of the spatial determinant
        do ix3=ixOmin3, ixOmax3
         do ix2=ixOmin2, ixOmax2
         do ix1=ixOmin1, ixOmax1
       call get_sqrtgamma(x(ix1,ix2,ix3,1),x(ix1,ix2,ix3,2),x(ix1,ix2,ix3,3),&
          m%sqrtgamma(ix1,ix2,ix3))
        end do 
         end do 
         end do 
       
       !------------------------------------------------------------
       ! Fill the derivative data:
       !------------------------------------------------------------

       if (fill_derivs) then

          ! set the spatial part:
          do inonzero=1,m%nnonzeroDgDk
             i = m%nonzeroDgDk(inonzero)%i
             j = m%nonzeroDgDk(inonzero)%j
             k = m%nonzeroDgDk(inonzero)%k

             ! Just spatial and only upper triangle (the others are already linking here)
             if (i == 0 .or. j == 0 .or. i .gt. j) cycle

              do ix3=ixOmin3, ixOmax3
               do ix2=ixOmin2, ixOmax2
               do ix1=ixOmin1, ixOmax1
             call get_g_component(i,j,x(ix1,ix2,ix3,1),x(ix1,ix2,ix3,2),x(ix1,&
                ix2,ix3,3),dummy,dgdk=m%nonzeroDgDk(inonzero)%elem(ix1,ix2,&
                ix3),kdir=k)
              end do 
               end do 
               end do 
          end do


          !------------------------------
          ! Calculate derivatives of lapse and contravariant shift:
          !------------------------------

          ! Numerically compute the derivatives of the inverse three-metric (later not needed by scheme):
          do k=1,3
              do ix3=ixOmin3, ixOmax3
               do ix2=ixOmin2, ixOmax2
               do ix1=ixOmin1, ixOmax1
             call get_dgammainvdk(x(ix1,ix2,ix3,1),x(ix1,ix2,ix3,2),x(ix1,ix2,&
                ix3,3),dgammainvdk(ix1,ix2,ix3,1:3,1:3,k),k)
              end do 
               end do 
               end do 
          end do

          ! Store the derivatives of the covariant shift and gtt (later not needed by scheme):
          do j=0,3
             do k=1,3
                 do ix3=ixOmin3, ixOmax3
                  do ix2=ixOmin2, ixOmax2
                  do ix1=ixOmin1, ixOmax1
                call get_g_component(0,j,x(ix1,ix2,ix3,1),x(ix1,ix2,ix3,2),&
                   x(ix1,ix2,ix3,3),dummy,dgdk=dbetaDownjDk(ix1,ix2,ix3,j,k),&
                   kdir=k)
                 end do 
                  end do 
                  end do 
             end do
          end do

          ! Get derivative of contravariant shift:
          ! d_j beta^i = d_j gamma^{ik}*beta_k
          do inonzero=1,m%nnonzeroDbetaiDj
             i = m%nonzeroDbetaiDj(inonzero)%i
             j = m%nonzeroDbetaiDj(inonzero)%j
             m%nonzeroDbetaiDj(inonzero)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3) = +  betaD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3,1)*dgammainvdk(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3,i,1,j)+ betaD(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)*dgammainvdk&
                (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,i,2,&
                j)+ betaD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                3)*dgammainvdk(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3,i,3,j) +  m%gammainv(i,1)%elem&
                (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                *dbetaDownjDk(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                1,j)+ m%gammainv(i,2)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)*dbetaDownjDk(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3,2,j)+ m%gammainv(i,3)%elem(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3)*dbetaDownjDk(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3,3,j)
          end do

          ! Derivative of lapse:
          ! d_j alpha=1/2*(beta^k*d_j beta_k + beta_k*d_j beta^k - d_j g_{00}) / alpha
          do inonzero=1, m%nnonzeroDalphaDj
             j = m%nonzeroDalphaDj(inonzero)%j

             m%nonzeroDalphaDj(inonzero)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3) = 0.5d0*( betaU(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*dbetaDownjDk&
                (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1,&
                j) + betaU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                2)*dbetaDownjDk(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3,2,j) + betaU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3,3)*dbetaDownjDk(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3,3,j)  +  betaD&
                (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                1)*m%DbetaiDj(1,j)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)+ betaD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3,2)*m%DbetaiDj(2,j)%elem(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3)+ betaD(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)*m%DbetaiDj(3,&
                j)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) - &
                dbetaDownjDk(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                0,j) )/m%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)
          end do

       end if

    !--------------------------------------------------
    end if
    !--------------------------------------------------

  end subroutine fill_metric
  !=============================================================================
  subroutine int_surface(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idims,x,dx1,dx2,dx3,func,&
     integral)
    ! integrates the function func over the cell-surface using 
    ! Simpsons rule. The cell extends over (x-dx/2, x+dx/2)
    use mod_amrvacdef
    integer, intent(in)                                               :: &
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,&
       ixOmin3,ixOmax1,ixOmax2,ixOmax3, idims
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)              :: x
    double precision, intent(in)                                      :: dx1,&
       dx2,dx3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), intent(out)                   :: integral
    double precision                                                  :: func
    ! .. local ..
    integer                :: is1,is2,is3, ix1,ix2,ix3
    double precision       :: tmp, xloc1,xloc2,xloc3
 double precision :: stmp(-1:1)
    logical, save          :: initialize=.true.
     double precision, save :: s(-1:1,-1:1)
    !$OMP threadprivate(s,initialize)
    !-----------------------------------------------------------------------------

    if (initialize) then 
       ! fill matrices for Simpsons rule:
       
        stmp = (/1.0d0,4.0d0,1.0d0/)
       do is2=-1,1
          do is1=-1,1
             s(is1,is2)=stmp(is1)*stmp(is2)
          end do
       end do
       
       initialize = .false.
    end if





    

    select case(idims)
       case(1)
        do ix3 = ixOmin3,ixOmax3
 do ix2 = ixOmin2,ixOmax2
 do ix1 = ixOmin1,ixOmax1
       integral(ix1,ix2,ix3) = 0.0d0
       xloc1 = x(ix1,ix2,ix3,1)
       do is2 = -1,1
          xloc2 = x(ix1,ix2,ix3,2)+dx2/2.0d0*dble(is2)
          do is3 = -1,1
             xloc3 = x(ix1,ix2,ix3,3)+dx3/2.0d0*dble(is3)
             call get_sqrtgamma(xloc1,xloc2,xloc3,tmp)
             integral(ix1,ix2,ix3) = integral(ix1,ix2,ix3) + tmp*func(xloc1,&
                xloc2,xloc3)*s(is2,is3)
          end do
       end do
       integral(ix1,ix2,ix3) = integral(ix1,ix2,ix3) * dx2*dx3/36.0d0
        end do
 end do
 end do

       case(2)
        do ix3 = ixOmin3,ixOmax3
 do ix2 = ixOmin2,ixOmax2
 do ix1 = ixOmin1,ixOmax1
       integral(ix1,ix2,ix3) = 0.0d0
       xloc2 = x(ix1,ix2,ix3,2)
       do is1 = -1,1
          xloc1 = x(ix1,ix2,ix3,1)+dx1/2.0d0*dble(is1)
          do is3 = -1,1
             xloc3 = x(ix1,ix2,ix3,3)+dx3/2.0d0*dble(is3)
             call get_sqrtgamma(xloc1,xloc2,xloc3,tmp)
             integral(ix1,ix2,ix3) = integral(ix1,ix2,ix3) + tmp*func(xloc1,&
                xloc2,xloc3)*s(is1,is3)
          end do
       end do
       integral(ix1,ix2,ix3) = integral(ix1,ix2,ix3) * dx1*dx3/36.0d0
        end do
 end do
 end do

       case(3)
        do ix3 = ixOmin3,ixOmax3
 do ix2 = ixOmin2,ixOmax2
 do ix1 = ixOmin1,ixOmax1
       integral(ix1,ix2,ix3) = 0.0d0
       xloc3 = x(ix1,ix2,ix3,3)
       do is1 = -1,1
          xloc1 = x(ix1,ix2,ix3,1)+dx1/2.0d0*dble(is1)
          do is2 = -1,1
             xloc2 = x(ix1,ix2,ix3,2)+dx2/2.0d0*dble(is2)
             call get_sqrtgamma(xloc1,xloc2,xloc3,tmp)
             integral(ix1,ix2,ix3) = integral(ix1,ix2,ix3) + tmp*func(xloc1,&
                xloc2,xloc3)*s(is1,is2)
          end do
       end do
       integral(ix1,ix2,ix3) = integral(ix1,ix2,ix3) * dx1*dx2/36.0d0
        end do
 end do
 end do
    end select


  end subroutine int_surface
  !=============================================================================
  subroutine int_volume(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,dx1,dx2,dx3,func,&
     integral)
    ! integrates the function func over the cell-volume using
    ! Simpsons rule.  The cell extends over (x-dx/2, x+dx/2)
    integer, intent(in)                                               :: &
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,&
       ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)              :: x
    double precision, intent(in)                                      :: dx1,&
       dx2,dx3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), intent(out)                   :: integral
    double precision                                                  :: func
    ! .. local ..
    integer                    :: is1,is2,is3, ix1,ix2,ix3
    double precision           :: tmp
    double precision           :: s1D(-1:1)
     double precision   ::  s2D(-1:1,-1:1)
    double precision,save      :: s(-1:1,-1:1,-1:1)
    logical, save                         :: initialize=.true.
    !$OMP threadprivate(s,initialize)
    !-----------------------------------------------------------------------------
    if (initialize) then
       ! fill matrices for Simpsons rule:
       s1D = (/1.0d0,4.0d0,1.0d0/)
       
       
       do is2=-1,1
          do is1=-1,1
             
             s2D(is1,is2)=s1D(is1)*s1D(is2)
             
          end do
       end do
      
       do is3=-1,1
       do is2=-1,1
       do is1=-1,1
                s(is1,is2,is3)=s2D(is1,is2)*s1D(is3)
       end do
       end do
       end do
      
       initialize = .false.
    end if

    !--------------------------------
    ! Now integrate the function over the volume
    ! Thus obtain an approximation for
    ! I=\int\int\int sqrt(gamma)*func dx dy dz
    !--------------------------------
    
     do ix3 = ixOmin3,ixOmax3
      do ix2 = ixOmin2,ixOmax2
      do ix1 = ixOmin1,ixOmax1
    integral(ix1,ix2,ix3) = 0.0d0

     do is3=-1,1
      do is2=-1,1
      do is1=-1,1
    call get_sqrtgamma(x(ix1,ix2,ix3,1)+dx1/2.0d0*dble(is1),x(ix1,ix2,ix3,&
       2)+dx2/2.0d0*dble(is2),x(ix1,ix2,ix3,3)+dx3/2.0d0*dble(is3),tmp)
    integral(ix1,ix2,ix3) = integral(ix1,ix2,ix3) + tmp*func(x(ix1,ix2,ix3,&
       1)+dx1/2.0d0*dble(is1),x(ix1,ix2,ix3,2)+dx2/2.0d0*dble(is2),x(ix1,ix2,&
       ix3,3)+dx3/2.0d0*dble(is3))*s(is1,is2,is3)
     end do 
      end do 
      end do 

    integral(ix1,ix2,ix3) = integral(ix1,ix2,ix3) *  dx1/6.0d0* dx2&
       /6.0d0* dx3/6.0d0

     end do
      end do
      end do
    
  end subroutine int_volume
  !=============================================================================
  subroutine complex_derivative(x1,x2,x3,func,jdir,derivative)
    ! Performs complex-step derivative for accurate numerical derivatives 
    ! of an analytic real-valued function func.  
    ! This simply takes advantage of the Cauchy-Riemann relations.  
    ! 
    ! Interface for func must look like this:
    !
    ! double complex function func(x^D)
    ! double complex, intent(in)        :: x^D

    use mod_amrvacdef
    
    double precision, intent(in)        :: x1,x2,x3
    integer, intent(in)                 :: jdir
    double precision, intent(out)       :: derivative

    interface
       double complex function func(x1,x2,x3)
         double complex, intent(in)     :: x1,x2,x3
       end function func
    end interface
    
    !  ..local..
    double complex                      :: xloc1,xloc2,xloc3
    !-----------------------------------------------------------------------------

    xloc1=cmplx(x1 , kr(jdir,1)*smalldouble , kind(1.d0))
    xloc2=cmplx(x2 , kr(jdir,2)*smalldouble , kind(1.d0))
    xloc3=cmplx(x3 , kr(jdir,3)*smalldouble , kind(1.d0));
    
    derivative = aimag(func(xloc1,xloc2,xloc3))/smalldouble
    
  end subroutine complex_derivative
  !=============================================================================
  double precision function ones(x1,x2,x3)
    double precision        :: x1,x2,x3
    !-----------------------------------------------------------------------------
    
    ones = 1.0d0

  end function ones
  !=============================================================================
  
  double precision function coordinate_x1(x1,x2,x3)
    double precision        :: x1,x2,x3
 !-----------------------------------------------------------------------------
    
    coordinate_x1 = x1

  end function coordinate_x1
 !=============================================================================,
  double precision function coordinate_x2(x1,x2,x3)
    double precision        :: x1,x2,x3
 !-----------------------------------------------------------------------------
    
    coordinate_x2 = x2

  end function coordinate_x2
 !=============================================================================,
  double precision function coordinate_x3(x1,x2,x3)
    double precision        :: x1,x2,x3
 !-----------------------------------------------------------------------------
    
    coordinate_x3 = x3

  end function coordinate_x3
 !=============================================================================
  subroutine fillgeo_covariant(pgeogrid,ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
     ixGmax2,ixGmax3,ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,ixGextmax2,&
     ixGextmax3,xmin1,xmin2,xmin3,dx1,dx2,dx3,need_only_volume)
    use mod_amrvacdef

    type(geoalloc) :: pgeogrid
    integer, intent(in) :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
        ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,ixGextmax2,ixGextmax3
    double precision, intent(in) :: xmin1,xmin2,xmin3, dx1,dx2,dx3
    logical, intent(in) :: need_only_volume
    ! .. local ..
    integer :: idims, idir, ixMmin1,ixMmin2,ixMmin3,ixMmax1,ixMmax2,ixMmax3,&
        ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3, ixCmin1,ixCmin2,ixCmin3,&
       ixCmax1,ixCmax2,ixCmax3, ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,&
       ixFmax3, ix1,ix2,ix3, ix, ixGsurfmin1,ixGsurfmin2,ixGsurfmin3,&
       ixGsurfmax1,ixGsurfmax2,ixGsurfmax3, imin1,imin2,imin3
    double precision :: x(ixGextmin1-1:ixGextmax1,ixGextmin2-1:ixGextmax2,&
       ixGextmin3-1:ixGextmax3,1:ndim),xi(ixGextmin1-1:ixGextmax1,&
       ixGextmin2-1:ixGextmax2,ixGextmin3-1:ixGextmax3,1:ndim)
    !-----------------------------------------------------------------------------
    ixMmin1=ixGmin1+dixB;ixMmin2=ixGmin2+dixB;ixMmin3=ixGmin3+dixB
    ixMmax1=ixGmax1-dixB;ixMmax2=ixGmax2-dixB;ixMmax3=ixGmax3-dixB;
    ixmin1=ixGmin1+1;ixmin2=ixGmin2+1;ixmin3=ixGmin3+1;ixmax1=ixGmax1-1
    ixmax2=ixGmax2-1;ixmax3=ixGmax3-1;
    ixGsurfmin1=ixGextmin1-1;ixGsurfmin2=ixGextmin2-1
    ixGsurfmin3=ixGextmin3-1;
    ixGsurfmax1=ixGextmax1;ixGsurfmax2=ixGextmax2;ixGsurfmax3=ixGextmax3;
     imin1 = nint((xmin1-xprobmin1)/dx1) 
      imin2 = nint((xmin2-xprobmin2)/dx2) 
      imin3 = nint((xmin3-xprobmin3)/dx3) 

    !--------------------------------------------------
    ! Cell center positions:
    !--------------------------------------------------
    do idims=1,ndim
       select case(idims)
          case(1)
          do ix = ixGsurfmin1,ixGsurfmax1
             x(ix,ixGsurfmin2:ixGsurfmax2,ixGsurfmin3:ixGsurfmax3,1)&
                =xprobmin1+(dble(imin1+ix-dixB)-half)*dx1
          end do
          case(2)
          do ix = ixGsurfmin2,ixGsurfmax2
             x(ixGsurfmin1:ixGsurfmax1,ix,ixGsurfmin3:ixGsurfmax3,2)&
                =xprobmin2+(dble(imin2+ix-dixB)-half)*dx2
          end do
          case(3)
          do ix = ixGsurfmin3,ixGsurfmax3
             x(ixGsurfmin1:ixGsurfmax1,ixGsurfmin2:ixGsurfmax2,ix,3)&
                =xprobmin3+(dble(imin3+ix-dixB)-half)*dx3
          end do
       end select
    end do

    !--------------------------------------------------
    ! Calculate the cell-volume with Simpsons rule:
    !--------------------------------------------------
    call int_volume(ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,ixGextmax2,&
       ixGextmax3,ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,ixGextmax2,&
       ixGextmax3,x(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,1:ndim),dx1,dx2,dx3,ones,pgeogrid%dvolume)     


    !--------------------------------------------------
    ! Calculate the barycenter, also with Simpsons-rule:
    !--------------------------------------------------
    
    call int_volume(ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,ixGextmax2,&
       ixGextmax3,ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,ixGextmax2,&
       ixGextmax3,x(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,1:ndim),dx1,dx2,dx3,coordinate_x1,&
       pgeogrid%xbar(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,1))
    pgeogrid%xbar(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,1) = pgeogrid%xbar(ixGextmin1:ixGextmax1,&
       ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,1) / pgeogrid%dvolume&
       (ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3)
    
    
    call int_volume(ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,ixGextmax2,&
       ixGextmax3,ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,ixGextmax2,&
       ixGextmax3,x(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,1:ndim),dx1,dx2,dx3,coordinate_x2,&
       pgeogrid%xbar(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,2))
    pgeogrid%xbar(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,2) = pgeogrid%xbar(ixGextmin1:ixGextmax1,&
       ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,2) / pgeogrid%dvolume&
       (ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3)
    
    
    call int_volume(ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,ixGextmax2,&
       ixGextmax3,ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,ixGextmax2,&
       ixGextmax3,x(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,1:ndim),dx1,dx2,dx3,coordinate_x3,&
       pgeogrid%xbar(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,3))
    pgeogrid%xbar(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,3) = pgeogrid%xbar(ixGextmin1:ixGextmax1,&
       ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,3) / pgeogrid%dvolume&
       (ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3)
    


    !--------------------------------------------------
    ! Fill metric at the barycenter position:
    !--------------------------------------------------
    call fill_metric(pgeogrid%m,ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,&
       ixGextmax2,ixGextmax3,ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,&
       ixGextmax2,ixGextmax3,pgeogrid%xbar)


    !--------------------------------------------------
    ! Get out when only volume is required
    !--------------------------------------------------
    if (need_only_volume) return


    !--------------------------------------------------
    ! Fill metric and surfaces at the interface/barycenter position:
    !--------------------------------------------------
    
    do idims=1,ndim
       select case(idims)
          case(1)

          do idir=1,ndim
             if (idir .eq. 1) cycle
             xi(ixGsurfmin1:ixGsurfmax1,ixGsurfmin2:ixGsurfmax2,&
                ixGsurfmin3:ixGsurfmax3,idir) = x(ixGsurfmin1:ixGsurfmax1,&
                ixGsurfmin2:ixGsurfmax2,ixGsurfmin3:ixGsurfmax3,idir)
          end do
          do ix = ixGsurfmin1,ixGsurfmax1
             xi(ix,ixGsurfmin2:ixGsurfmax2,ixGsurfmin3:ixGsurfmax3,1)&
                =xprobmin1+dble(imin1+ix-dixB)*dx1
          end do
          
          
          ixCmin1=ixmin1-kr(1,1);ixCmin2=ixmin2-kr(1,2)
          ixCmin3=ixmin3-kr(1,3); ixCmax1=ixmax1;ixCmax2=ixmax2
          ixCmax3=ixmax3;
          ixFmin1=ixCmin1-1;ixFmin2=ixCmin2-1;ixFmin3=ixCmin3-1
          ixFmax1=ixCmax1+1;ixFmax2=ixCmax2+1;ixFmax3=ixCmax3+1;
          ! Metric at interface xi:
          call fill_metric(pgeogrid%mSurface1,ixGsurfmin1,ixGsurfmin2,&
             ixGsurfmin3,ixGsurfmax1,ixGsurfmax2,ixGsurfmax3,ixGsurfmin1,&
             ixGsurfmin2,ixGsurfmin3,ixGsurfmax1,ixGsurfmax2,ixGsurfmax3,&
             xi(ixGsurfmin1:ixGsurfmax1,ixGsurfmin2:ixGsurfmax2,&
             ixGsurfmin3:ixGsurfmax3,1:ndim),need_derivs=.false.)
          ! Surface at interface xi:
          call int_surface(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,&
             ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,idims,&
             xi(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,1:ndim),dx1,&
             dx2,dx3,ones,pgeogrid%surfaceC1)
          ! Surface at barycenter:
          call int_surface(ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,&
             ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,&
             pgeogrid%xbar(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             1:ndim),dx1,dx2,dx3,ones,pgeogrid%surface1)
          
          case(2)

          do idir=1,ndim
             if (idir .eq. 2) cycle
             xi(ixGsurfmin1:ixGsurfmax1,ixGsurfmin2:ixGsurfmax2,&
                ixGsurfmin3:ixGsurfmax3,idir) = x(ixGsurfmin1:ixGsurfmax1,&
                ixGsurfmin2:ixGsurfmax2,ixGsurfmin3:ixGsurfmax3,idir)
          end do
          do ix = ixGsurfmin2,ixGsurfmax2
             xi(ixGsurfmin1:ixGsurfmax1,ix,ixGsurfmin3:ixGsurfmax3,2)&
                =xprobmin2+dble(imin2+ix-dixB)*dx2
          end do
          
          
          ixCmin1=ixmin1-kr(2,1);ixCmin2=ixmin2-kr(2,2)
          ixCmin3=ixmin3-kr(2,3); ixCmax1=ixmax1;ixCmax2=ixmax2
          ixCmax3=ixmax3;
          ixFmin1=ixCmin1-1;ixFmin2=ixCmin2-1;ixFmin3=ixCmin3-1
          ixFmax1=ixCmax1+1;ixFmax2=ixCmax2+1;ixFmax3=ixCmax3+1;
          ! Metric at interface xi:
          call fill_metric(pgeogrid%mSurface2,ixGsurfmin1,ixGsurfmin2,&
             ixGsurfmin3,ixGsurfmax1,ixGsurfmax2,ixGsurfmax3,ixGsurfmin1,&
             ixGsurfmin2,ixGsurfmin3,ixGsurfmax1,ixGsurfmax2,ixGsurfmax3,&
             xi(ixGsurfmin1:ixGsurfmax1,ixGsurfmin2:ixGsurfmax2,&
             ixGsurfmin3:ixGsurfmax3,1:ndim),need_derivs=.false.)
          ! Surface at interface xi:
          call int_surface(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,&
             ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,idims,&
             xi(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,1:ndim),dx1,&
             dx2,dx3,ones,pgeogrid%surfaceC2)
          ! Surface at barycenter:
          call int_surface(ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,&
             ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,&
             pgeogrid%xbar(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             1:ndim),dx1,dx2,dx3,ones,pgeogrid%surface2)
          
          case(3)

          do idir=1,ndim
             if (idir .eq. 3) cycle
             xi(ixGsurfmin1:ixGsurfmax1,ixGsurfmin2:ixGsurfmax2,&
                ixGsurfmin3:ixGsurfmax3,idir) = x(ixGsurfmin1:ixGsurfmax1,&
                ixGsurfmin2:ixGsurfmax2,ixGsurfmin3:ixGsurfmax3,idir)
          end do
          do ix = ixGsurfmin3,ixGsurfmax3
             xi(ixGsurfmin1:ixGsurfmax1,ixGsurfmin2:ixGsurfmax2,ix,3)&
                =xprobmin3+dble(imin3+ix-dixB)*dx3
          end do
          
          
          ixCmin1=ixmin1-kr(3,1);ixCmin2=ixmin2-kr(3,2)
          ixCmin3=ixmin3-kr(3,3); ixCmax1=ixmax1;ixCmax2=ixmax2
          ixCmax3=ixmax3;
          ixFmin1=ixCmin1-1;ixFmin2=ixCmin2-1;ixFmin3=ixCmin3-1
          ixFmax1=ixCmax1+1;ixFmax2=ixCmax2+1;ixFmax3=ixCmax3+1;
          ! Metric at interface xi:
          call fill_metric(pgeogrid%mSurface3,ixGsurfmin1,ixGsurfmin2,&
             ixGsurfmin3,ixGsurfmax1,ixGsurfmax2,ixGsurfmax3,ixGsurfmin1,&
             ixGsurfmin2,ixGsurfmin3,ixGsurfmax1,ixGsurfmax2,ixGsurfmax3,&
             xi(ixGsurfmin1:ixGsurfmax1,ixGsurfmin2:ixGsurfmax2,&
             ixGsurfmin3:ixGsurfmax3,1:ndim),need_derivs=.false.)
          ! Surface at interface xi:
          call int_surface(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,&
             ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,idims,&
             xi(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,1:ndim),dx1,&
             dx2,dx3,ones,pgeogrid%surfaceC3)
          ! Surface at barycenter:
          call int_surface(ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,&
             ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,&
             pgeogrid%xbar(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             1:ndim),dx1,dx2,dx3,ones,pgeogrid%surface3)
          
       end select
    end do

    !--------------------------------------------------
    ! Fill the dx:
    !--------------------------------------------------
    pgeogrid%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,1)=dx1
    pgeogrid%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,2)=dx2
    pgeogrid%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,3)=dx3;

  end subroutine fillgeo_covariant
  !=============================================================================
end module mod_metric
!=============================================================================
