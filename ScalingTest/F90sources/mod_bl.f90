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
  ! This wraps a module around mod_coord_bl.t to be able to use
  ! the Boyer-Lindquist metric elements for initialization of a problem
  ! to be run in different coordinates.
  !
  ! Usage:
  ! use mod_bl, get_alpha_BL => get_alpha, get_beta_BL => get_beta, &
  ! get_g_component_BL => get_g_component, u4BLtoCoord_BL => u4BLtoCoord
  ! 
  ! Oliver Porth
  ! 2016-01-19
  !=============================================================================

!=============================================================================
module mod_bl

  use mod_metric_aux
  implicit none
  
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
  ! Metric components for Boyer-Lindquist coordinates -coord=bl
  !
  ! Check BL-coordinates.nb for the definitions.
  ! Oliver Porth
  ! 04.12.2015
  !=============================================================================

  !-----------------------------------------
  ! Define constants specific to coordinates
  !-----------------------------------------

  character*20, parameter              :: coord="bl"
  logical, save                        :: init_from_g4 = .false. !We don't provide the four-metric


  ! Migrating to have coordinate-specific constants here.
  ! a_ and m_ however are in the physics part
  ! so put the MKS parameters, e.g. R0_ and h_ here:
  integer, parameter                   :: ncoordpar = 0
  double precision, save               :: coordpar(ncoordpar) 

contains
!=============================================================================
  subroutine init_coord

   get_gammainv_component_analytic => bl_get_gammainv_component_analytic
   get_gammainv_component => get_gammainv_component_analytic

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

    if(present(is_analytic)) is_analytic = .true.

    
    sqrtgamma = sqrt(((x1**2  &
         + eqpar(a_)**2*Cos(x3)**2)*Sin(x3)**2*((eqpar(a_)**2  &
         + x1**2)**2-1.0d0*eqpar(a_)**2*(eqpar(a_)**2  &
         + x1*(-2.0d0*eqpar(m_)+x1))*Sin(x3)**2))/(eqpar(a_)**2  &
         + x1*(-2.0d0*eqpar(m_)+x1)))
   


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

    if(present(iszero)) iszero = .false.

    
    alpha = sqrt((eqpar(a_)**2-2.0d0*eqpar(m_)*x1+x1**2)*(x1**2+  &
         eqpar(a_)**2*Cos(x3)**2))/Sqrt((eqpar(a_)**2+  &
         x1**2)**2-1.0d0*eqpar(a_)**2*(eqpar(a_)**2-  &
         2*eqpar(m_)*x1+x1**2)*Sin(x3)**2)
   

    if (present(jdir)) then

       select case(jdir)

       case(1)
          ! Radial derivative:
          if (present(dalphadj)) then
             
             dalphadj = (sqrt(2.0d0)*eqpar(m_)*(-1.0d0*(eqpar(a_)**2+  &
                  x1**2)**2*(eqpar(a_)**2-2*x1**2+  &
                  eqpar(a_)**2*Cos(2*x3))-  &
                  8.0d0*eqpar(a_)**2*eqpar(m_)*x1**3*Sin(x3)**2))&
                     /(Sqrt((eqpar(a_)**2  &
                  +x1*(-2.0d0*eqpar(m_)+x1))*(x1**2+  &
                  eqpar(a_)**2*Cos(x3)**2))*(eqpar(a_)**4+2.0d0*x1**4+  &
                  eqpar(a_)**2*x1*(2.0d0*eqpar(m_)+3.0d0*x1)+  &
                  eqpar(a_)**2*(eqpar(a_)**2+x1*(-2.0d0*eqpar(m_)+  &
                  x1))*Cos(2.0d0*x3))**1.50d0)
            
          end if
          if (present(dalphadj_iszero)) then
             dalphadj_iszero = .false.
          end if

       case(3)
          ! Theta-derivative:
          
          if (present(dalphadj)) then
             dalphadj = (-2.0d0*sqrt(2.0d0)*eqpar(a_)**2*eqpar(m_)*x1*(eqpar&
                (a_)**2  &
                  +x1**2)*(eqpar(a_)**2+x1*(-2*eqpar(m_)+  &
                  x1))*Sin(2*x3))/(Sqrt((eqpar(a_)**2+x1*(-2*eqpar(m_)+  &
                  x1))*(x1**2+eqpar(a_)**2*Cos(x3)**2))*(eqpar(a_)**4+  &
                  2*x1**4+eqpar(a_)**2*x1*(2*eqpar(m_)+3*x1)+  &
                  eqpar(a_)**2*(eqpar(a_)**2+x1*(-2*eqpar(m_)+  &
                  x1))*Cos(2*x3))**1.50d0)
          end if
          if (present(dalphadj_iszero)) then
             dalphadj_iszero = .false.
          end if
         

       case default
          ! Other cases (phi):
          if (present(dalphadj)) then
             dalphadj = 0.0d0
          end if
          if (present(dalphadj_iszero)) then
             dalphadj_iszero = .true.
          end if
       end select

    end if ! present(jdir)

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

    select case(idir)

    case(2)
       ! Betaphi
       
       beta = (-2.0d0*eqpar(a_)*eqpar(m_)*x1)/((eqpar(a_)**2+  &
            x1**2)**2-eqpar(a_)**2*(eqpar(a_)**2-2*eqpar(m_)*x1+  &
            x1**2)*Sin(x3)**2)
      

       if(present(iszero)) then
          iszero = .false.
       end if

       if (present(jdir)) then
          select case(jdir)

          case(1)
             ! dbetaphidr:
             if (present(dbetaidj)) then
                
                dbetaidj = (eqpar(a_)*eqpar(m_)*(-1.0d0*eqpar(a_)**4+  &
                     3.0d0*eqpar(a_)**2*x1**2+6.0d0*x1**4+  &
                     eqpar(a_)**2*(-1.0d0*eqpar(a_)**2+  &
                     x1**2)*Cos(2.0d0*x3)))/((eqpar(a_)**2+x1**2)**2-  &
                     1.0d0*eqpar(a_)**2*(eqpar(a_)**2+x1*(-2*eqpar(m_)+  &
                     x1))*Sin(x3)**2)**2
               
             end if ! present(dbetaidj)

             if (present(dbetaidj_iszero)) then
                dbetaidj_iszero = .false.
             end if
          case(3)
             ! dbetaphidtheta:
             
             if (present(dbetaidj)) then
                dbetaidj = (-4.0d0*eqpar(a_)**3*eqpar(m_)*x1*(eqpar(a_)**2+  &
                     x1*(-2*eqpar(m_)+  &
                     x1))*Cos(x3)*Sin(x3))/((eqpar(a_)**2+x1**2)**2-  &
                     eqpar(a_)**2*(eqpar(a_)**2+x1*(-2*eqpar(m_)+  &
                     x1))*Sin(x3)**2)**2
             end if
             if (present(dbetaidj_iszero)) then
                dbetaidj_iszero = .false.
             end if
            


          case default
             ! dbbetaphidphi is zero
             if (present(dbetaidj)) then
                dbetaidj = 0.0d0
             end if ! present(dbetaidj)           

             if (present(dbetaidj_iszero)) then
                dbetaidj_iszero = .true.
             end if

          end select
       end if ! present(jdir)


    case default
       ! All zeroes except phi:
       beta = 0.0d0
       if(present(iszero)) iszero = .true.

       ! \partial_j \beta^r
       if (present(dbetaidj)) dbetaidj = 0.0d0
       if (present(dbetaidj_iszero)) dbetaidj_iszero = .true.

    end select

  end subroutine get_beta
  !=============================================================================
  subroutine bl_get_gammainv_component_analytic(iin,jin,x1,x2,x3,ginv,iszero,&
     dginvdk_iszero,dginvdk,kdir)

    use mod_amrvacdef

    ! Indices of the 3-metric are up (contravariant) gammainv^{ij}
    ! The optional argument iszero is true if the element is identically zero
    ! The optional arguments dginvdk and kdir request derivatives of the metric
    ! \partial_k gammainv^{ij} ; i=iin, j=jin, k=kdir
    ! The optional argument dginvdk_iszero

    integer, intent(in)                      :: iin,jin
    integer, optional, intent(in)            :: kdir
    double precision, intent(in)             :: x1,x2,x3
    double precision, intent(out)            :: ginv
    logical, optional, intent(out)           :: iszero, dginvdk_iszero
    double precision, optional, intent(out)  :: dginvdk
    integer                                  :: i,j
    double precision                         :: Delta, Sigma, A
    double precision                         :: d1Delta, d1Sigma, d2Sigma,&
        d1A, d2A
    !-----------------------------------------------------------------------------

    if(present(dginvdk) .and. .not. present(kdir) .or. present&
       (dginvdk_iszero) .and. .not. present(kdir)) call mpistop&
       ("get_gammainv_component: derivatives requested without &direction or &
       output-slot given.")

    ! metric is symmetric: swap indices if needed:
    ! User needs only to provide values for i<=j (upper triangle).  
    if (iin>jin) then
       i=jin; j=iin
    else
       i=iin; j=jin
    end if

    ! Auxiliary quantities
    Delta = x1**2 -2.d0*eqpar(m_)*x1 +eqpar(a_)**2
    
    Sigma = x1**2 + eqpar(a_)**2*cos(x3)**2
    A = (x1**2 + eqpar(a_)**2)**2 - Delta*eqpar(a_)**2*sin(x3)**2
   

    ! Depending on the direction of derivation,
    ! compute derivatives of the auxiliary quantities
    if (present(kdir)) then
       if (kdir .eq. 1) then
          d1Delta = 2.d0*x1 - 2.d0*eqpar(m_)
          d1Sigma = 2.d0*x1
          
          d1A = 4.d0*x1*(eqpar(a_)**2 + x1**2) - eqpar(a_)**2*sin(x3)&
             **2*d1Delta
         
       else if (kdir .eq. 2) then
          
          d2Sigma = -2.d0 * eqpar(a_)**2 * sin(x3) * cos(x3)
          d2A = -2.d0 * Delta * eqpar(a_)**2 * sin(x3) * cos(x3)
         
       end if
    end if
    

    ! Diagonal for now
    if (i==j) then
       if(present(iszero)) iszero = .false.

       if (i .eq. 1) then ! r direction
          ginv = Delta / Sigma
          if (present(kdir)) then
             if (kdir .eq. 1) then
                if (present(dginvdk)) dginvdk = - Delta / Sigma&
                   **2 * d1Sigma + d1Delta / Sigma
                if (present(dginvdk_iszero)) dginvdk_iszero = .false.
             else if (kdir .eq. 3) then
                if (present(dginvdk)) dginvdk = - Delta / Sigma**2 * d2Sigma
                if (present(dginvdk_iszero)) dginvdk_iszero = .false.
             else
                if (present(dginvdk)) dginvdk = 0.0d0
                if (present(dginvdk_iszero)) dginvdk_iszero = .true.
             end if
          end if

       else if (i .eq. 2) then ! phi direction
          
          ginv = Sigma / sin(x3)**2 / A
         

          if (present(kdir)) then
             if (kdir .eq. 1) then
                
                if (present(dginvdk)) dginvdk = (-Sigma/A**2*d1A + d1Sigma&
                   /A) / sin(x3)**2
                if (present(dginvdk_iszero)) dginvdk_iszero = .false.
               
             else if (kdir .eq. 3) then
                
                if (present(dginvdk)) dginvdk = (-Sigma/A**2*d2A + d2Sigma&
                   /A) / sin(x3)**2 &
                                              - 2.d0*cos(x3)/sin(x3)**3*Sigma&
                                                 /A
                if (present(dginvdk_iszero)) dginvdk_iszero = .false.
               
             else
                if (present(dginvdk)) dginvdk = 0.0d0
                if (present(dginvdk_iszero)) dginvdk_iszero = .true.
             end if
          end if

       else if (i .eq. 3) then ! theta direction
          ginv = 1.d0 / Sigma

          if (present(kdir)) then
             if (kdir .eq. 1) then
                if (present(dginvdk)) dginvdk = - d1Sigma / Sigma**2
                if (present(dginvdk_iszero)) dginvdk_iszero = .false.
             else if (kdir .eq. 3) then
                
                if (present(dginvdk)) dginvdk = - d2Sigma / Sigma**2
                if (present(dginvdk_iszero)) dginvdk_iszero = .false.
               
             else
                if (present(dginvdk)) dginvdk = 0.d0
                if (present(dginvdk_iszero)) dginvdk_iszero = .true.
             end if
          end if

       end if
    else
       ginv = 0.0d0

       if(present(iszero)) iszero = .true.

       if (present(dginvdk)) dginvdk = 0.0d0

       if(present(dginvdk_iszero)) dginvdk_iszero = .true.

    end if
    
  end subroutine bl_get_gammainv_component_analytic

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

    ! Diagonal for now
    if (i==j) then
       if(present(iszero)) iszero = .false.

       if (i .eq. 1) then
          !grr:
          
          g = (x1**2+eqpar(a_)**2*Cos(x3)**2)/(eqpar(a_)**2+  &
               x1*(-2.0d0*eqpar(m_)+x1))
         
          if (present(kdir)) then
             ! Derivative info requested

             select case(kdir)

             case (1)
                ! dgrrdr:
                if (present(dgdk)) then
                   
                   dgdk = (2.0d0*(x1*(eqpar(a_)**2-eqpar(m_)*x1)+  &
                        eqpar(a_)**2*(eqpar(m_)-  &
                        x1)*Cos(x3)**2))/(eqpar(a_)**2+x1*(-2*eqpar(m_)+  &
                        x1))**2
                  
                end if
                if (present(dgdk_iszero)) dgdk_iszero = .false.

             case(3)
                !dgrrdtheta:
                
                if (present(dgdk)) then
                   dgdk = (-2.0d0*eqpar(a_)**2*Cos(x3)*Sin(x3))&
                      /(eqpar(a_)**2  &
                        -2*eqpar(m_)*x1+x1**2)
                end if
                if (present(dgdk_iszero)) dgdk_iszero = .false.
               

             case default
                ! phi derivatives zero:
                if (present(dgdk)) dgdk = 0.0d0
                if (present(dgdk_iszero)) dgdk_iszero = .true.

             end select
          end if

       else if (i .eq. 2) then
          ! gphiphi
          
          g = (Sin(x3)**2*((eqpar(a_)**2+x1**2)**2-  &
               1.0d0*eqpar(a_)**2*(eqpar(a_)**2+x1*(-2*eqpar(m_)+  &
               x1))*Sin(x3)**2))/(x1**2+eqpar(a_)**2*Cos(x3)**2)
         

          if (present(kdir)) then
             ! Derivative info requested
             select case(kdir)
             case(1)
                !DgphiphiDr
                
                if (present(dgdk)) dgdk = (2.0d0*x1*(eqpar(a_)**2+x1**2)*(x1&
                   **2+  &
                     eqpar(a_)**2*Cos(2*x3))*Sin(x3)**2+  &
                     2.0d0*eqpar(a_)**2*(x1*(eqpar(a_)**2-eqpar(m_)*x1)+  &
                     eqpar(a_)**2*(eqpar(m_)-  &
                     x1)*Cos(x3)**2)*Sin(x3)**4)/(x1**2+  &
                     eqpar(a_)**2*Cos(x3)**2)**2
                if (present(dgdk_iszero)) dgdk_iszero = .false.
               
             case(3)
                !DgphiphiDtheta
                
                if (present(dgdk)) dgdk = (eqpar(a_)**2+x1*(-2.0d0*eqpar&
                   (m_)+x1)+  &
                     (8.0d0*eqpar(m_)*x1*(eqpar(a_)**2+  &
                     x1**2)**2)/(eqpar(a_)**2+2*x1**2+  &
                     eqpar(a_)**2*Cos(2*x3))**2)*Sin(2.0d0*x3)
                if (present(dgdk_iszero)) dgdk_iszero = .false.
               
             case default
                !DphiphiDphi
                if (present(dgdk)) dgdk = 0.0d0
                if (present(dgdk_iszero)) dgdk_iszero = .true.
             end select
          end if

       else if (i .eq. 3) then
          ! gthetatheta
          
          g = x1**2+eqpar(a_)**2*Cos(x3)**2
         

          if (present(kdir)) then
             ! Derivative info requested
             select case(kdir)

             case(1)
                ! DgthetathetaDr:
                if (present(dgdk)) dgdk = 2.0d0*x1
                if (present(dgdk_iszero)) dgdk_iszero = .false.

             case(3)
                ! DthetathetaDtheta:
                
                if (present(dgdk)) &
                     dgdk = -2.0d0*eqpar(a_)**2*Cos(x3)*Sin(x3)
                if (present(dgdk_iszero)) dgdk_iszero = .false.
               

             case default
                ! DthetathetaDphi
                if (present(dgdk)) dgdk = 0.0d0
                if (present(dgdk_iszero)) dgdk_iszero = .true.

             end select
          end if

       end if
    else
       g = 0.0d0

       if(present(iszero)) iszero = .true.

       if (present(dgdk)) dgdk = 0.0d0

       if(present(dgdk_iszero)) dgdk_iszero = .true.

    end if

  end subroutine get_g_component
  !=============================================================================
  double precision function outerhorizon()

    use mod_amrvacdef
    !-----------------------------------------------------------------------------

    outerhorizon = eqpar(m_) + sqrt(eqpar(m_)**2 - eqpar(a_)**2)

  end function outerhorizon
  !=============================================================================
  subroutine BLToCoord(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xBL,xCoord)

    ! identity
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

    xCoord = xBL

  end subroutine BLToCoord
  !=============================================================================
  subroutine CoordToBL(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xCoord,xBL)

    ! Trivial case
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

    xBL = xCoord

  end subroutine CoordToBL
  !=============================================================================
  subroutine u4BLtoCoord(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,u4BL,u4Coord)

    ! Transforms the (contravariant) four-velocity u4BL from Boyer-Lindquist coordinates
    ! to the current (BL) coordinates u4Coord.  Often initial conditions are
    ! given in terms of BL coordinates and this routine comes in handy.
    !
    ! This should just return the same vector since we are already in BL-coord.
    !
    use mod_amrvacdef

    integer, intent(in)                                    :: ixImin1,ixImin2,&
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

    u4Coord(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,0:3) &
       = u4BL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,0:3)

  end subroutine u4BLtoCoord
  !=============================================================================
  subroutine CoordToCart(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,xCart)

    ! Transforms the coordinates to "Boyer-Lindquist Cartesian coordinates"
    ! x = ra Sin(thetaBL) Cos(phiBL)
    ! y = ra Sin(thetaBL) Sin(phiBL)
    ! z = rBL Cos(thetaBL)
    ! where ra=sqrt(rBL**2+a**2) ! WARNING: I decided to put ra=rBL for now...
    ! We are already in BL-cordinates, so we only have to do one transformation.
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
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)                     :: sth, cth, sph, cph, ra
    integer                                                :: zcart_, ycart_
    !-----------------------------------------------------------------------------

    if (ndim.eq.3) then
       ycart_=2
       zcart_=3
    else if (ndim.eq.2 .and. 3.eq.2) then
       ycart_=1 ! must catch ycart_=1 case and don't use
       zcart_=2
    else if (ndim.eq.2 .and. 2.eq.2) then
       ycart_=2
       zcart_=1 ! must catch zcart_=1 case and don't use
    else if (ndim.eq.1) then
       ycart_=1
       zcart_=1
    else
       call mpistop("CoordToCart: unknown parameter combination of -d,&
           -phi and -z !")
    end if
    !-----------------------------------------------------------------------------

    
    sph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       sin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,phi_))
    cph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       cos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,phi_))
   

    
    sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       sin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,z_))
    cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       cos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,z_))
   

    ra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) !sqrt(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2+eqpar(a_)**2)


    xCart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) &
       = ra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) * &
       sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*cph&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    if (ycart_.ne.1) xCart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       ycart_) = ra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) * &
       sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*sph&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    if (zcart_.ne.1) xCart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       zcart_) = x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)*cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

  end subroutine CoordToCart
  !=============================================================================
  subroutine u4CoordToCart(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,u4Coord,u4Cart,J)

    ! Transforms any contravariant four vector u4Coord in the current coordinates
    ! to a vector in the basis of "Boyer-Lindquist Cartesian coordinates".
    ! t = tBL
    ! x = sqrt(rBL**2+a**2) Sin(thetaBL) Cos(phiBL)
    ! y = sqrt(rBL**2+a**2) Sin(thetaBL) Sin(phiBL)
    ! z = rBL Cos(thetaBL)
    ! We are already in BL-cordinates, so we only have to do one transformation.  
    ! The transformation matrix J=del(t,x,y,z)/del(tBL,rBL,thetaBL,phiBL)
    !
    !     | 1               0                             0                        0                          |
    !     | 0   rBL/ra Sin(thetaBL) Cos(phiBl)    ra Cos(thetaBL) Cos(phiBL)    - ra Sin(thetaBL) Sin(phiBL)  |
    ! J = | 0   rBL/ra Sin(thetaBL) Sin(phiBL)    ra Cos(thetaBL) Sin(phiBL)      ra Sin(thetaBL) Cos(phiBL)  |
    !     | 0   Cos(thetaBL)                    - rBL Sin(thetaBL)                 0                          |
    !
    ! where ra=sqrt(rBL**2+a**2) ! WARNING: I decided to put ra=rBL for now...
    ! 
    ! u4Cart = J u4Coord
    !
    ! 3D simulation, three vector components, (-d=33 -z=2, phi=3 , the standard):
    ! u4Coord = (utBL,urBL,uthetaBL,uphiBL)
    ! u4Cart  = (ut,ux,uy,uz)
    !
    ! r-theta plane (phi=0), three vector components (-d=23 -z=2 -phi=3):
    ! u4Coord = (utBL,urBL,uthetaBL,uphiBL)
    ! u4Cart = (ut,ux,uy,uz)
    !
    ! r-theta plane (phi = 0), two vector components (-d=22 -z=2 -phi=0):
    ! u4Coord = (utBL,urBL,uthetaBL)
    ! u4Cart = (ut,ux,uz)
    !
    ! r-phi plane (theta = pi/2), three vector components (-d=23 -phi=2 -z=3)
    ! u4Coord = (utBL,urBL,uphiBL,uthetaBL)
    ! u4Cart = (ut,ux,uy,uz)
    !
    ! r-phi plane (theta = pi/2), two vector components (-d=22 -phi=2 -z=0)
    ! u4Coord = (utBL,urBL,uphiBL)
    ! uCart = (ut,ux,uy)
    !-----------------------------------------------------------------------------

    use mod_amrvacdef

    integer, intent(in)                                    :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3), intent(in)   :: x
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(in)   :: u4Coord
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3), intent(out)  :: u4Cart
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3), optional, intent(out)  :: J
    ! .. local ..
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,0:3,0:3)         :: Jac
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)                     :: sth, cth, sph, cph, ra
    integer                                                :: zcart_, ycart_,&
        ix,jx
    !-----------------------------------------------------------------------------

    if (ndim.eq.3.and.ndir.eq.3) then
       ycart_=2
       zcart_=3
    else if ((ndim.eq.2.or.ndim.eq.1).and.ndir.eq.3 .and. 2 .eq. 2) then
       ycart_=2
       zcart_=3
    else if ((ndim.eq.2.or.ndim.eq.1).and.ndir.eq.3 .and. 3 .eq. 2) then
       ycart_=3
       zcart_=2
    else if (ndir.eq.2 .and. 3.eq.2) then
       ycart_=0 ! must catch ycart_=0 case and don't use
       zcart_=2
    else if (ndir.eq.2 .and. 2.eq.2) then
       ycart_=2
       zcart_=0 ! must catch zcart_=0 case and don't use
    else if (ndir.eq.1) then
       ycart_=0
       zcart_=0
    else
       call mpistop("u4CoordToCart: unknown parameter combination of -d,&
           -phi and -z !")
    end if


    
    sph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       sin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,phi_))
    cph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       cos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,phi_))
   

    
    sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       sin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,z_))
    cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       cos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,z_))
   

    ra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) !sqrt(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2+eqpar(a_)**2)

    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:,:) = zero
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,0,0) = one

    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1,1) &
       = x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)&
       /ra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) * &
       sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*cph&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1,3) &
       = ra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) * &
       cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*cph&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
   
    
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1,2) &
       = - ra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) * &
       sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*sph&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
   

    if (ycart_.ne.0) then
       Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ycart_,1) &
          = x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)&
          /ra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) * &
          sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*sph&
          (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
       
       Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ycart_,3) &
          = ra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) * &
          cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*sph&
          (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      
       
       Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ycart_,2) &
          = ra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) * &
          sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*cph&
          (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      
    end if

    if(zcart_.ne.0) then
       Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,zcart_,1) &
          = cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
       
       Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,zcart_,3) &
          = - x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          1)*sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      
       
       Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,zcart_,2) = 0.0d0
      
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Jacobian fully assembled, now
    ! transform contravariant four-vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    u4Cart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:) = zero
    do ix=0,ndir
       do jx=0,ndir
          u4Cart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ix) &
             = u4Cart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             ix) + Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ix,&
             jx) * u4Coord(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,jx)
       end do
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if (present(J)) J = Jac

  end subroutine u4CoordToCart
  !=============================================================================
  subroutine u3CoordToCart(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,u3Coord,u3Cart,J)

    ! Transforms any contravariant three vector u3Coord in the current coordinates
    ! to a vector in the basis of "Boyer-Lindquist Cartesian coordinates".
    ! x = sqrt(rBL**2+a**2) Sin(thetaBL) Cos(phiBL)
    ! y = sqrt(rBL**2+a**2) Sin(thetaBL) Sin(phiBL)
    ! z = rBL Cos(thetaBL)
    ! We are already in BL-cordinates, so we only have to do one transformation.  
    ! The transformation matrix J=del(x,y,z)/del(rBL,thetaBL,phiBL)
    ! Since the coordinates don't depend on time, this is the same transformation as for the
    ! four-vector u4CoordToCart()
    !
    !     | rBL/ra Sin(thetaBL) Cos(phiBl)    ra Cos(thetaBL) Cos(phiBL)    - ra Sin(thetaBL) Sin(phiBL)  |
    ! J = | rBL/ra Sin(thetaBL) Sin(phiBL)    ra Cos(thetaBL) Sin(phiBL)      ra Sin(thetaBL) Cos(phiBL)  |
    !     | Cos(thetaBL)                    - rBL Sin(thetaBL)                 0                          |
    !
    ! where ra=sqrt(rBL**2+a**2) ! WARNING: I decided to put ra=rBL for now...
    ! 
    ! u3Cart = J u3Coord
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
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)                     :: sth, cth, sph, cph, ra
    integer                                                :: zcart_, ycart_,&
        ix,jx
    !-----------------------------------------------------------------------------

    if (ndim.eq.3.and.ndir.eq.3) then
       ycart_=2
       zcart_=3
    else if ((ndim.eq.2.or.ndim.eq.1).and.ndir.eq.3 .and. 2 .eq. 2) then
       ycart_=2
       zcart_=3
    else if ((ndim.eq.2.or.ndim.eq.1).and.ndir.eq.3 .and. 3 .eq. 2) then
       ycart_=3
       zcart_=2
    else if (ndir.eq.2 .and. 3.eq.2) then
       ycart_=1 ! must catch ycart_=1 case and don't use
       zcart_=2
    else if (ndir.eq.2 .and. 2.eq.2) then
       ycart_=2
       zcart_=1 ! must catch zcart_=1 case and don't use
    else if (ndir.eq.1) then
       ycart_=1
       zcart_=1
    else
       call mpistop("u3CoordToCart: unknown parameter combination of -d,&
           -phi and -z !")
    end if


    
    sph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       sin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,phi_))
    cph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       cos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,phi_))
   

    
    sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       sin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,z_))
    cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       cos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,z_))
   

    ra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) !sqrt(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2+eqpar(a_)**2)

    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1,1) &
       = x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)&
       /ra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) * &
       sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*cph&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1,3) &
       = ra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) * &
       cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*cph&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
   
    
    Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1,2) &
       = - ra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*sth&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*sph(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)
   

    if (ycart_.ne.1) then
       Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ycart_,1) &
          = x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)&
          /ra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) * &
          sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*sph&
          (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
       
       Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ycart_,3) &
          = ra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) * &
          cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*sph&
          (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      
       
       Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ycart_,2) &
          = ra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*sth&
          (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*cph&
          (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      
    end if

    if(zcart_.ne.1) then
       Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,zcart_,1) &
          = cth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
       
       Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,zcart_,3) &
          = - x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          1)*sth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      
       
       Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,zcart_,2) = 0.0d0
      
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Jacobian fully assembled, now
    ! transform contravariant three-vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    u3Cart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:) = zero
    do ix=1,ndir
       do jx=1,ndir
          u3Cart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ix) &
             = u3Cart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             ix) + Jac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ix,&
             jx) * u3Coord(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,jx)
       end do
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if (present(J)) J = Jac

  end subroutine u3CoordToCart
  !=============================================================================
  ! Dummies
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

    call mpistop("get_dgammainvdk: Not required and not implemented.")

    dgammainvdk = 0.0d0
    
  end subroutine get_dgammainvdk
  !=============================================================================
  ! End of coordinate-specific definitions.
  !=============================================================================


  !=============================================================================
  subroutine get_g4_BL(x1,x2,x3,g)
    ! Return the four-metric at point x^D
    double precision, intent(in)                           :: x1,x2,x3
    double precision, dimension(0:3,0:3),intent(out)   :: g
    ! .. local ..
    integer                                                :: i, j
    double precision, dimension(1:3)                     :: betaU, betaD
    double precision                                       :: beta2
    !-----------------------------------------------------------------------------

    do i=1,3
       do j=1,3
          call get_g_component(i,j,x1,x2,x3,g(i,j))
       end do
    end do
    call get_alpha(x1,x2,x3,g(0,0))
    do i=1,3
       call get_beta(i,x1,x2,x3,betaU(i))
    end do
    ! lower the beta:
    betaD(:) = 0.0d0
    do i=1,3
       do j=1,3
          betaD(i) = betaD(i) + g(i,j)*betaU(j)
       end do
    end do

    
    beta2 = betaU(1)*betaD(1)+betaU(2)*betaD(2)+betaU(3)*betaD(3)
    
    g(0,0) = -g(0,0)**2 + beta2

    do i=1,3
       g(i,0) = betaD(i)
       g(0,i) = betaD(i)
    end do


  end subroutine get_g4_BL
  !=============================================================================

  
end module mod_bl
