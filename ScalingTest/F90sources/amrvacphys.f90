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

  !#############################################################################
  ! Covariant rmhd module, based on the public srmhd version of MPI-AMRVAC
  ! 20.10.2015, Oliver Porth
  
  !=============================================================================
  subroutine getdt(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dtnew,dx1,dx2,dx3,x)

    use mod_amrvacdef

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: dx1,dx2,dx3, x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), dtnew
    !-----------------------------------------------------------------------------

    dtnew=bigdouble

  end subroutine getdt
  !=============================================================================
  subroutine checkglobaldata
    use mod_amrvacdef
    !-----------------------------------------------------------------------------
    minrho = max(zero,smallrho)
    
    
    govergminone =eqpar(gamma_)/(eqpar(gamma_)-one)
    minp  = max(zero,smallp)
    smallxi=minrho+minp*govergminone
    smalltau = minp/(eqpar(gamma_)-one)
    
    

    if (useprimitiveRel .eqv. .false.) call mpistop&
       ('useprimitiveRel = .false. abandonned')

  end subroutine checkglobaldata
  !=============================================================================
  subroutine initglobaldata

    ! place to initialize some globals

    use mod_con2prim, only: setup_con2prim
    use mod_amrvacdef
    !-----------------------------------------------------------------------------
    eqpar(adiab_)=1.0d0

    
    eqpar(gamma_)=4.0d0/3.0d0
    
    

    

    eqpar(a_) = 0.0d0
    eqpar(m_) = 1.0d0

    
    
    limitvalue=smalldouble**2

    call setup_con2prim

  end subroutine initglobaldata
  !=============================================================================
  subroutine conserve(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,patchw)

    ! Transform primitive variables into conservative ones
    ! (rho,v,p,B) ---> (D,S,tau,B,lfac,xi)
    ! v is contravariant, S is covarariant  
    ! call to smallvalues
    ! --> latter only used for correcting procedure in correctaux
    ! --> input array patchw for spatially selective transformation

    use mod_amrvacdef

    integer, intent(in)               :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout)   :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    logical, intent(in)               :: patchw(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3)
    double precision, intent(in)      :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    !-----------------------------------------------------------------------------

    call conserven(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,patchw)

  end subroutine conserve
  !=============================================================================
  subroutine conserven(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,patchw)

    ! Transform primitive variables into conservative ones
    ! (rho,v,p,B) ---> (D,S,tau,B,lfac,xi)
    ! no call to smallvalues
    ! --> latter only used for correcting procedure in correctaux
    ! --> input array patchw for spatially selective transformation

    use mod_metric, only: lower3
    use mod_amrvacdef

    integer, intent(in)               :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout)   :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    logical, intent(in)               :: patchw(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3)

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3):: sqrV,sqrU,sqrB,VdotB,rhoh

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3)  :: vD, bD
    !-----------------------------------------------------------------------------

    call lower3(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,myM,w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,u1_:u3_),vD)
    call lower3(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,myM,w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,b1_:b3_),bD)

    if(useprimitiveRel)then
       ! assumes four velocity computed in primitive (rho u p B) with u=lfac*v
       where(.not.patchw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
          sqrU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
             =w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             u1_)*vD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             1)+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             u2_)*vD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             2)+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             u3_)*vD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3) 
          sqrV(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
             =sqrU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
             /(one+sqrU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
          sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
             =w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             b1_)*bD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             1)+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             b2_)*bD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             2)+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             b3_)*bD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)
       endwhere
    else
       ! assumes velocity in primitive (rho v p B) 
       where(.not.patchw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
          sqrV(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
             =w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             v1_)*vD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             1)+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             v2_)*vD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             2)+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             v3_)*vD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)
          sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
             =w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             b1_)*bD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             1)+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             b2_)*bD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             2)+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             b3_)*bD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3) 
       endwhere
    endif

    ! fill the auxilary variable lfac (lorentz factor)
    if(useprimitiveRel)then
       where(.not.patchw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)&
             =dsqrt(one+sqrU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)) 
          VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)  &
             =( bD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             1)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             u1_)+ bD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             2)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             u2_)+ bD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             3)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,u3_))&
             /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)
       endwhere
    else
       where(.not.patchw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)=one&
             /dsqrt(one-sqrV(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
          VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)  &
             = bD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             1)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             v1_)+ bD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             2)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             v2_)+ bD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             3)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,v3_)
       endwhere
    endif

    ! fill the auxilary variable xi and density D
    ! with enthalpy w: xi= lfac^2 rhoh
    ! density: d = lfac * rho
    call Enthalpy(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,patchw,rhoh)
    where(.not.patchw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_)&
          =w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          lfac_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          lfac_)* rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,d_) &
          =w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          rho_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)
       ! recycle sqrU array for storing temporary positive 
       ! array for use in energy variable
       sqrU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
          =sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*sqrV&
          (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)-VdotB&
          (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)**2
    endwhere
    where(.not.patchw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
       .and.sqrU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)<zero)
       sqrU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=zero
    endwhere

    

    

    

    
    ! fill the vector S
    ! s= (xi + B^2) * v - (v.B) * B
    if(useprimitiveRel)then
       where(.not.patchw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,s1_)&
             =(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             xi_)+sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))*vD&
             (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)&
             /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             lfac_) - VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)*bD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1)
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,s2_)&
             =(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             xi_)+sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))*vD&
             (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)&
             /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             lfac_) - VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)*bD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2)
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,s3_)&
             =(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             xi_)+sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))*vD&
             (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)&
             /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             lfac_) - VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)*bD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3);
       endwhere
    else
       where(.not.patchw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,s1_)&
             =(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             xi_)+sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))*vD&
             (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             1) - VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*bD&
             (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,s2_)&
             =(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             xi_)+sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))*vD&
             (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             2) - VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*bD&
             (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,s3_)&
             =(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             xi_)+sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))*vD&
             (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             3) - VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*bD&
             (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3);
       endwhere
    endif

    !{#IFDEF ENERGY
    ! E = xi - p +B^2/2 + (v^2 B^2 - (v.B)^2)/2 
    ! instead of E use tau= E - D
    where(.not.patchw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,tau_)&
          =w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          xi_) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          pp_) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          d_) +half*(sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + &
          sqrU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
    endwhere
    !}

  end subroutine conserven
  !=============================================================================
  subroutine primitive(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)

    ! Transform conservative variables into primitive ones
    ! (D,S,tau,B)-->(rho,v,p,B,lfac,xi)

    use mod_amrvacdef

    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)

    logical, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)          :: &
       patchw

    !-----------------------------------------------------------------------------
    patchw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = .false.

    ! calculate lorentz factor and xi from conservatives only
    call getaux(.true.,w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,'primitive')

    call primitiven(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,patchw)

    if (tlow>zero) call fixp_usr(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)


    
  end subroutine primitive
  !==============================================================================
  subroutine primitiven(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,patchw)

    !  --> needed in correctaux, avoiding recursive calling by smallvalues
    ! assumes updated xi and lfac
    ! Transform conservative variables into primitive ones
    ! (D,S,tau,B)-->(rho,v,p,B,lfac,xi)

    use mod_metric, only: square3u, raise3
    use mod_amrvacdef

    integer, intent(in)                    :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout)        :: w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:nw)
    logical, intent(in),dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3)   :: patchw
    ! .. local ..
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)     :: SdotB, sqrB,tmpP
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3) :: sU
    !-----------------------------------------------------------------------------

    !{#IFDEF ENERGY
    call Pressuren(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,.true.,tmpP,patchw)
    where(.not.patchw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,pp_)&
          =tmpP(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end where
    !}

    call square3u(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,myM,w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,b1_:b3_),sqrB(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3))
    where(.not.patchw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)&
          =w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,d_)&
          /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)
       SdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
          = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          s1_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          b1_)+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          s2_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          b2_)+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          s3_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b3_)
    end where

    call raise3(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,myM,w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,s1_:s3_),sU)
    if(useprimitiveRel)then   
       where (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_)<smallxi)
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,u1_)=zero
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,u2_)=zero
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,u3_)=zero
       elsewhere
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,u1_)&
             =w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             lfac_)*(sU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             1)+SdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*w&
             (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b1_)&
             /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_))&
             / (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             xi_)+sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,u2_)&
             =w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             lfac_)*(sU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             2)+SdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*w&
             (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b2_)&
             /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_))&
             / (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             xi_)+sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,u3_)&
             =w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             lfac_)*(sU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             3)+SdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*w&
             (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b3_)&
             /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_))&
             / (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             xi_)+sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
       end where
    else
       where (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_)<smallxi)
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,v1_)=zero
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,v2_)=zero
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,v3_)=zero
       elsewhere
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,v1_)&
             =(sU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             1)+SdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*w&
             (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b1_)&
             /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_))&
             /(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             xi_)+sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,v2_)&
             =(sU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             2)+SdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*w&
             (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b2_)&
             /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_))&
             /(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             xi_)+sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,v3_)&
             =(sU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             3)+SdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*w&
             (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b3_)&
             /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_))&
             /(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             xi_)+sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
       endwhere
    endif

    

    
    
    
    
 
  end subroutine primitiven
  !=============================================================================
  subroutine getv(wprim,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idims,v)

    ! Returns the contravariant (velocity*lapse - shift)
    ! Now takes input in primitive form for performance!
    use mod_amrvacdef

    integer, intent(in)                              :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,idims
    double precision, intent(in)                     :: wprim(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:nw)
    double precision, intent(in)                     :: x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), intent(out)  :: v
    !-----------------------------------------------------------------------------

    v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = myM%alpha&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*wprim&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,u0_+idims)&
       /wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       lfac_) - myM%beta(idims)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

  end subroutine getv
  !=============================================================================
  subroutine getcmax(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idims,cmax,cmin,needcmin)

    ! Calculate cmax_idim within ixO^L

    use mod_metric, only: square3u
    use mod_amrvacdef

    logical, intent(in)               :: needcmin
    integer, intent(in)               :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idims
    double precision, intent(inout)   :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)      :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)     :: cmax(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),cmin(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    ! .. local ..
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)  :: sqrB,csound2
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)  :: A,B,C
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)  :: B2
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)  :: sUidims
    integer                             :: j
    !-----------------------------------------------------------------------------
 
    ! C = VdotB
    C(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)= (w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,s1_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,b1_)+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       s2_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       b2_)+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       s3_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b3_))&
       /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_)

    ! squared Eulerian field strength (big B^2)
    call square3u(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,myM,w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,b1_:b3_),B2)

    ! squared comoving field strength (little b^2)
    sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)/w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)**2 + C(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)**2

    ! the square of sound speed using EOS in eos.t
    call getcsound2(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,.true.,csound2)

    ! contravariant momentum used for velocity
    ! get only the required component:
    sUidims(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = zero
    do j=1,3
       sUidims(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
          = sUidims(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + &
          myM%gammainv(idims,j)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          s0_+j)
    end do

    select case(typepoly)

    case('gammie')

       A(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)   &
          = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_)&
          /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)**2 !rhoh
       ! squared Alfven speed
       B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)  &
          = sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
          /(sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)+A&
          (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
       ! equation 72 Del zanna et al. 2007A&A...473...11D
       A(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)  &
          = csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + &
          B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) - &
          csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*B&
          (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

       sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)     = one - one&
          /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)**2 !v2
       csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)  &
          = (sUidims(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + &
          C(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) * &
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+idims)) &
          / ( w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          xi_)+B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) ) !vidim

       ! reuse B for square root:
       B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) =  &
          A(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) * &
          (one-sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)) * ( &
          (one-sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*A&
          (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)) * &
          myM%gammainv(idims,idims)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) - csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)**2 * (one-A(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)) )

       where(B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) .gt. zero)
         B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
            = sqrt(B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)) !the standard case
       elsewhere
         B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = zero !the highly relativistic limit
       endwhere

       ! reuse C for first term:
       C(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
          csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*(one-A&
          (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))

       ! characteristic speed Eq. 76 Del zanna et al. 2007A&A...473...11D
       cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
          = ( C(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + &
          B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) ) &
          / ( one - sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*A&
          (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) )
       cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
          = ( C(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) - &
          B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) ) &
          / ( one - sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*A&
          (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) )

    case default

       call mpistop("getcmax: Unknown polynomial type")

    end select
    
    ! Limit by speed of light:
    cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       max(cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3), - 1.0d0&
       /sqrt(myM%g(idims,idims)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)))
    cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       min(cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),   1.0d0&
       /sqrt(myM%g(idims,idims)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)))
    cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       max(cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3), - 1.0d0&
       /sqrt(myM%g(idims,idims)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)))
    cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       min(cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),   1.0d0&
       /sqrt(myM%g(idims,idims)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)))
    
    ! add shift-part and mutliply by lapse
    if (.not. needcmin) then
       cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
          = max( abs(myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) - myM%beta(idims)%elem(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3)), abs(myM%alpha(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3)*cmin(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3) - myM%beta(idims)%elem&
          (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)) )
    else
       cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
          = myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*cmax&
          (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) - &
          myM%beta(idims)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
          = myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*cmin&
          (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) - &
          myM%beta(idims)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
    end if

  end subroutine getcmax
  !=============================================================================
  subroutine getflux(w,wprim,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idims,f,&
     transport)

    ! Calculate non-transport flux f_idim[iw] within ixO^L.
    ! 
    use mod_metric
    use mod_amrvacdef
    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idims
    double precision, intent(in)       :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw) !conserved variables
    double precision, intent(in)       :: wprim(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,nw) !primitive variables
    double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)      :: f(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nwflux)
    logical, intent(out)               :: transport(1:nwflux)
    ! .. local ..
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) :: VdotB, sqrB, ptot, v
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)  :: alphabjU, alphabjD, vU
    integer                                    :: i, iw
    !-----------------------------------------------------------------------------
    transport=.true.

    ! S is covariant, b is contravariant:
    VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
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
    
    do i=1, 3
       vU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,i) &
          = wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,u0_+i)&
          /wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)
       alphabjU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,i) &
          = myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) * &
          (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+i)&
          /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)&
          **2 + VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*vU&
          (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,i) )
    end do
    call lower3(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,myM,alphabjU,alphabjD)
    ptot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       half*(VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
       **2  + sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
       /(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)&
       **2)) + wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,pp_)
    

    do iw = 1, nwflux
       
       select case(iw)

       case(d_)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)=zero

          
          

          

          
          
          



          case(b1_)    !!! Fi[B_c],Fi[B_c],Fi[B_c] = B_c*v_i - v_c*B_i 
          if (idims==1) then
 !f_i[b_i] should be exactly 0, so we do not use the transport flux
             f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)=zero
             transport(iw)=.false.            
          else
             ! transport velocity in direction 1
             call getv(wprim,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,1,v)
             ! multiply by B in direction idims
             f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) &
                = - v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) * &
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+idims)
          end if
          

          
          case(b2_)    !!! Fi[B_c],Fi[B_c],Fi[B_c] = B_c*v_i - v_c*B_i 
          if (idims==2) then
 !f_i[b_i] should be exactly 0, so we do not use the transport flux
             f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)=zero
             transport(iw)=.false.            
          else
             ! transport velocity in direction 2
             call getv(wprim,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,2,v)
             ! multiply by B in direction idims
             f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) &
                = - v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) * &
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+idims)
          end if
          

          
          case(b3_)    !!! Fi[B_c],Fi[B_c],Fi[B_c] = B_c*v_i - v_c*B_i 
          if (idims==3) then
 !f_i[b_i] should be exactly 0, so we do not use the transport flux
             f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)=zero
             transport(iw)=.false.            
          else
             ! transport velocity in direction 3
             call getv(wprim,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,3,v)
             ! multiply by B in direction idims
             f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) &
                = - v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) * &
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+idims)
          end if
          

          

       case(e_)
          where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_)&
             <smallxi)
             f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)=zero
          elsewhere
             ! f=ptot*v^i
             f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) &
                = ptot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*vU&
                (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idims)
             ! f=ptot*v^i - v.B*B^i
             f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) &
                = f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                iw) - VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3,b0_+idims)
             ! f=(ptot*v^i - v.B*B^i)*alpha
             f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) &
                = f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                iw) * myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)
          end where

       case(s1_) ! momentum is covariant, transport flux done.
          ! f=-bj Bi /lfac

          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) &
             = - alphabjD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             1) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+idims)

          if (idims==1) then !!! add ptot
             f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
                =f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                iw) + ptot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3) * myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)
          end if
          
       case(s2_) ! momentum is covariant, transport flux done.
          ! f=-bj Bi /lfac

          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) &
             = - alphabjD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             2) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+idims)

          if (idims==2) then !!! add ptot
             f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
                =f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                iw) + ptot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3) * myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)
          end if
          
       case(s3_) ! momentum is covariant, transport flux done.
          ! f=-bj Bi /lfac

          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) &
             = - alphabjD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             3) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+idims)

          if (idims==3) then !!! add ptot
             f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
                =f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                iw) + ptot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3) * myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)
          end if
          

       end select

    end do

  end subroutine getflux
  !=============================================================================
  subroutine addgeometry(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)

    ! Add geometrical source terms to w

    use mod_metric, only: lower3, raise3, dalphadj_is_zero, beta_is_zero,&
        dbetaidj_is_zero, dgdk_is_zero, square3u
    use mod_amrvacdef

    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)       :: qdt
    double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    ! .. local ..
    integer                            :: iw, i,j,k, inonzero, jnonzero,&
        hxOmin1,hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3, idims
    double precision,dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)  :: vU, bsU, bsD, sU, wikdjgammaik
    double precision,dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)         :: VdotB, sqrB, tmp, tmp2, wikbetajdjgammaik
    double precision,dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)         :: ptot
    !-----------------------------------------------------------------------------

    VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)= &
       (wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       s1_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       b1_)+wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       s2_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       b2_)+wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       s3_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b3_))&
       /wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_)
    call square3u(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,myM,wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,b1_:b3_),sqrB)
    call raise3(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,myM,wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,s1_:s3_),sU)

    
         vU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)&
            =(sU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            1)+VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*wCT&
            (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+1))/ &
         (wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            xi_)+sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
    bsU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) &
       = wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+1)&
       /wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       lfac_) + wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       lfac_)*VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*vU&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)
    
    
         vU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)&
            =(sU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            2)+VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*wCT&
            (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+2))/ &
         (wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            xi_)+sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
    bsU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) &
       = wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+2)&
       /wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       lfac_) + wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       lfac_)*VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*vU&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)
    
    
         vU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)&
            =(sU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            3)+VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*wCT&
            (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+3))/ &
         (wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            xi_)+sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
    bsU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3) &
       = wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+3)&
       /wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       lfac_) + wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       lfac_)*VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*vU&
       (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)
    
    call lower3(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,myM,bsU,bsD)

    call Pressure(wCT,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,.true.,ptot)
    ptot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       half*(VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
       **2  + sqrB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
       /wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)&
       **2)  + ptot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    ! wikdjgammaik is without the total pressure term
    wikdjgammaik(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:) = zero
    do inonzero = 1, myM%nnonzeroDgDk
       i = myM%nonzeroDgDk(inonzero)%i
       k = myM%nonzeroDgDk(inonzero)%j
       j = myM%nonzeroDgDk(inonzero)%k

            !        + ptot(ixO^S)*myM%gammainv(i,k)%elem(ixO^S) &
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
          = sU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          i)*vU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          k) - bsU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          i)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+k)&
          /wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)
       wikdjgammaik(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,j) &
          = wikdjgammaik(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          j) + tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) * &
          myM%DgDk(i,k,j)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
    end do


    wikbetajdjgammaik(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = zero
    do inonzero = 1, myM%nnonzeroDgDk
       i = myM%nonzeroDgDk(inonzero)%i
       k = myM%nonzeroDgDk(inonzero)%j
       j = myM%nonzeroDgDk(inonzero)%k

       if  (beta_is_zero(j)) cycle

       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
          = sU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          i)*vU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          k) + ptot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
          *myM%gammainv(i,k)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) - bsU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,i)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,b0_+k)/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,lfac_)

       wikbetajdjgammaik(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
          = wikbetajdjgammaik(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) + tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) * myM%beta(j)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) * myM%DgDk(i,k,j)%elem(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end do


    do iw = 1, nwflux
       select case(iw)
          
 !s[s1_] = 1/2 alpha W**ik d1dgammaik + S_i d1beta**i -U d1alpha
               case(s1_)
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
             = half * myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * wikdjgammaik(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1)
          ! Treat the total pressure separately to make discretized gradient:
          hxOmin1=ixOmin1-kr(1,1);hxOmin2=ixOmin2-kr(2,1)
          hxOmin3=ixOmin3-kr(3,1);hxOmax1=ixOmax1-kr(1,1)
          hxOmax2=ixOmax2-kr(2,1);hxOmax3=ixOmax3-kr(3,1);
          idims = 1
          select case(idims)
             case(1)
             tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + &
                myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                *ptot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                  *(mygeo%surfaceC1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)-mygeo%surfaceC1(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2,hxOmin3:hxOmax3)) &
                  /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)
          case(2)
             tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + &
                myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                *ptot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                  *(mygeo%surfaceC2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)-mygeo%surfaceC2(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2,hxOmin3:hxOmax3)) &
                  /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)
          case(3)
             tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + &
                myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                *ptot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                  *(mygeo%surfaceC3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)-mygeo%surfaceC3(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2,hxOmin3:hxOmax3)) &
                  /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)
          end select

          do i = 1, 3
             if (dbetaidj_is_zero(i,1)) cycle
             tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + &
                wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                s0_+i) * myM%dbetaiDj(i,1)%elem(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3)
          end do

          if (.not. dalphadj_is_zero(1)) then
             tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) - &
                (wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                tau_)+wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                d_))*myM%dalphaDj(1)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)
          end if


          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) &
             = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             iw) + qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
          
          
 !s[s2_] = 1/2 alpha W**ik d2dgammaik + S_i d2beta**i -U d2alpha
               case(s2_)
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
             = half * myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * wikdjgammaik(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2)
          ! Treat the total pressure separately to make discretized gradient:
          hxOmin1=ixOmin1-kr(1,2);hxOmin2=ixOmin2-kr(2,2)
          hxOmin3=ixOmin3-kr(3,2);hxOmax1=ixOmax1-kr(1,2)
          hxOmax2=ixOmax2-kr(2,2);hxOmax3=ixOmax3-kr(3,2);
          idims = 2
          select case(idims)
             case(1)
             tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + &
                myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                *ptot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                  *(mygeo%surfaceC1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)-mygeo%surfaceC1(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2,hxOmin3:hxOmax3)) &
                  /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)
          case(2)
             tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + &
                myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                *ptot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                  *(mygeo%surfaceC2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)-mygeo%surfaceC2(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2,hxOmin3:hxOmax3)) &
                  /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)
          case(3)
             tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + &
                myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                *ptot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                  *(mygeo%surfaceC3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)-mygeo%surfaceC3(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2,hxOmin3:hxOmax3)) &
                  /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)
          end select

          do i = 1, 3
             if (dbetaidj_is_zero(i,2)) cycle
             tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + &
                wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                s0_+i) * myM%dbetaiDj(i,2)%elem(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3)
          end do

          if (.not. dalphadj_is_zero(2)) then
             tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) - &
                (wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                tau_)+wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                d_))*myM%dalphaDj(2)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)
          end if


          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) &
             = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             iw) + qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
          
          
 !s[s3_] = 1/2 alpha W**ik d3dgammaik + S_i d3beta**i -U d3alpha
               case(s3_)
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
             = half * myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * wikdjgammaik(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3)
          ! Treat the total pressure separately to make discretized gradient:
          hxOmin1=ixOmin1-kr(1,3);hxOmin2=ixOmin2-kr(2,3)
          hxOmin3=ixOmin3-kr(3,3);hxOmax1=ixOmax1-kr(1,3)
          hxOmax2=ixOmax2-kr(2,3);hxOmax3=ixOmax3-kr(3,3);
          idims = 3
          select case(idims)
             case(1)
             tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + &
                myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                *ptot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                  *(mygeo%surfaceC1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)-mygeo%surfaceC1(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2,hxOmin3:hxOmax3)) &
                  /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)
          case(2)
             tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + &
                myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                *ptot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                  *(mygeo%surfaceC2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)-mygeo%surfaceC2(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2,hxOmin3:hxOmax3)) &
                  /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)
          case(3)
             tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + &
                myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
                *ptot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                  *(mygeo%surfaceC3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)-mygeo%surfaceC3(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2,hxOmin3:hxOmax3)) &
                  /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)
          end select

          do i = 1, 3
             if (dbetaidj_is_zero(i,3)) cycle
             tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + &
                wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                s0_+i) * myM%dbetaiDj(i,3)%elem(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3)
          end do

          if (.not. dalphadj_is_zero(3)) then
             tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) - &
                (wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                tau_)+wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                d_))*myM%dalphaDj(3)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)
          end if


          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) &
             = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             iw) + qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
          

          ! s[tau_] = 1/2 W**ik beta**j dgammaikdj + W**j_i*dbetaidj - S**j dalphadj
       case(tau_)
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
             = half*wikbetajdjgammaik(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)
          do inonzero = 1, myM%nnonzeroDbetaiDj
             i = myM%nonzeroDbetaiDj(inonzero)%i
             j = myM%nonzeroDbetaiDj(inonzero)%j
             tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                = vU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                j)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,s0_+i)
             if (i .eq. j) tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3) = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3) + ptot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)
             tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) - &
                bsD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                i)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+j)&
                /wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)
             !            print*, 'Source:',i,j,maxval(abs(tmp2(ixO^S)))
             tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) + &
                tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)  * &
                myM%nonzeroDbetaiDj(inonzero)%elem(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3)
          end do

          do inonzero = 1, myM%nnonzeroDalphaDj
             j = myM%nonzeroDalphaDj(inonzero)%j
             tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) - &
                sU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                j)*myM%nonzeroDalphaDj(inonzero)%elem(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3)
          end do

          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) &
             = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             iw) + qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

          
       case default

       end select
    end do

    call getaux(.true.,w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,'addgeometry')

  end subroutine addgeometry
  !=============================================================================
  subroutine addsource(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iwmin,iwmax,qtC,wCT,qt,w,&
     x,qsourcesplit)

    ! w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO

    use mod_amrvacdef
    use mod_electronphys

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
        iwmin,iwmax
    double precision, intent(in)    :: qdt, qtC, qt
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    logical, intent(in)             :: qsourcesplit

    double precision :: dx1,dx2,dx3
    !-----------------------------------------------------------------------------

    dx1=dxlevel(1);dx2=dxlevel(2);dx3=dxlevel(3);

    ! Sources related to div B

    

    
    
  end subroutine addsource
  !=============================================================================
  subroutine getcurrent(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,current)

    use mod_metric, only: lower3
    use mod_amrvacdef

    integer, intent(in)           :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(out) :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3)

    ! .. local ..
    double precision              :: bvec(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)
    !-----------------------------------------------------------------------------

    call lower3(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1-1,&
       ixOmin2-1,ixOmin3-1,ixOmax1+1,ixOmax2+1,ixOmax3+1,myM,&
       w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,b1_:b3_),bvec)

    call curl3(bvec,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,current)

  end subroutine getcurrent
  !=============================================================================

  !=============================================================================
  ! just dummies:
  !=============================================================================
  subroutine e_to_rhos(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)

    use mod_amrvacdef
    integer:: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,&
       ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision:: w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,nw)
    double precision, intent(in)      :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    !-------------------------------------------------------------------------

    call mpistop("e to rhos unavailable")

  end subroutine e_to_rhos
  !=============================================================================
  subroutine rhos_to_e(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)

    use mod_amrvacdef

    integer:: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,&
       ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision:: w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,nw)
    double precision, intent(in)      :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    !-----------------------------------------------------------------------------

    call mpistop("e to rhos unavailable")

  end subroutine rhos_to_e
  !=============================================================================
  subroutine ppmflatcd(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ixLmin1,ixLmin2,ixLmin3,&
     ixLmax1,ixLmax2,ixLmax3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,&
     w,d2w,drho,dp)

    use mod_amrvacdef

    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ixLmin1,&
       ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3,ixRmin1,ixRmin2,ixRmin3,&
       ixRmax1,ixRmax2,ixRmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw),d2w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
       1:nwflux)

    double precision, intent(inout) :: drho(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3),dp(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)

    !-----------------------------------------------------------------------------

    call mpistop("PPM with flatsh=.true. not implemented for rmhd")

  end subroutine ppmflatcd
  !=============================================================================
  subroutine ppmflatsh(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ixLLmin1,ixLLmin2,&
     ixLLmin3,ixLLmax1,ixLLmax2,ixLLmax3,ixLmin1,ixLmin2,ixLmin3,ixLmax1,&
     ixLmax2,ixLmax3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,ixRRmin1,&
     ixRRmin2,ixRRmin3,ixRRmax1,ixRRmax2,ixRRmax3,idims,w,drho,dp,dv)

    use mod_amrvacdef

    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ixLLmin1,&
       ixLLmin2,ixLLmin3,ixLLmax1,ixLLmax2,ixLLmax3,ixLmin1,ixLmin2,ixLmin3,&
       ixLmax1,ixLmax2,ixLmax3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,&
       ixRmax3,ixRRmin1,ixRRmin2,ixRRmin3,ixRRmax1,ixRRmax2,ixRRmax3
    integer, intent(in)           :: idims
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)

    double precision, intent(inout) :: drho(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3),dp(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
       dv(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)

    double precision :: v(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
       ptot(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
    !-----------------------------------------------------------------------------

    call mpistop("PPM with flatsh=.true. not implemented for rmhd")

  end subroutine ppmflatsh
  !=============================================================================
  ! end module amrvacphys
  !=============================================================================
