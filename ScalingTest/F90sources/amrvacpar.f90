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
module mod_amrvacpar
  
  implicit none

  character*5,parameter:: typephys='rmhd'            ! VACPHYS module name

  
  
  
  CHARACTER*9,PARAMETER:: eqparname='gamma adiab m a' !Equation parameter names
  
  

  


  ! flow variables
  !=====Conserve variables=====!
  integer,parameter:: d_=1
  integer,parameter:: s0_=d_
  integer,parameter:: s1_=s0_+1,s2_=s0_+2,s3_=s0_+3
  integer,parameter:: e_=s3_+1
  integer,parameter:: tau_=e_
  integer,parameter:: b0_=e_
  integer,parameter:: b1_=b0_+1,b2_=b0_+2,b3_=b0_+3
  integer,parameter:: rhos_=e_

  

  
  
  
  integer,parameter:: Ds_=e_
  integer,parameter:: nwmhd=b3_
  
  
  
  ! Number of variables with tracer
  integer,parameter:: nwfluxtr=nwmhd
  
  
  
  ! Number of variables
  integer,parameter:: nwflux=nwfluxtr
  
  integer,parameter:: lfac_=nwflux+1     ! Lorentz factor
  integer,parameter:: xi_=lfac_+1       ! lfac2,lfac2,lfac2 Enthalpy
  !=============================!

  
  integer,parameter:: nwextra=0
  

  integer,parameter:: nwaux=2
  integer,parameter:: nw=nwflux+nwaux+nwextra

  !=====Primitive variables=====!
  integer,parameter:: rho_=d_
  !---- 3-velocities ----!
  integer,parameter:: v0_=s0_
  integer,parameter:: v1_=v0_+1,v2_=v0_+2,v3_=v0_+3
  !---- 4-velocities ----!
  integer,parameter:: u0_=s0_
  integer,parameter:: u1_=u0_+1,u2_=u0_+2,u3_=u0_+3
  
  
  
  integer,parameter:: pp_=e_
  integer,parameter:: s_=Ds_
  !=============================!
  ! polar variable names
  integer,parameter:: sr_=s0_+1
  integer,parameter:: sphi_=s0_+2
  integer,parameter:: sz_=s0_+3
  integer,parameter:: vr_=v0_+1
  integer,parameter:: vphi_=v0_+2
  integer,parameter:: vz_=v0_+3
  integer,parameter:: uz_=v0_+3
  integer,parameter:: ur_=v0_+1
  integer,parameter:: uphi_=v0_+2
  integer,parameter:: br_=b0_+1
  integer,parameter:: bphi_=b0_+2
  integer,parameter:: bz_=b0_+3
  integer,parameter:: ee_=e_

  integer, parameter :: nvector=2 !No. vector vars
  integer, dimension(nvector), parameter :: iw_vector=(/ s0_, b0_ /)

  integer :: itmp
  
  ! Staggered magnetic field
  integer, parameter :: bs0_ = 0
  integer, parameter :: bs1_=bs0_+1,bs2_=bs0_+2,bs3_=bs0_+3
  integer, parameter :: nws=3
  ! Staggered counterparts for centered variables:
  integer, parameter  :: iws(nwflux+nwaux) = [(itmp-b0_, itmp &
     = 1, nwflux+nwaux)]
  ! The number of staggered components of the vector fields is
  ! the number of dimensions. (In 2D, the 3rd component is centred)
  !integer, parameter :: nws=3
  

  integer,parameter:: fastRW_=1,fastLW_=2,slowRW_=3,slowLW_=4 ! Characteristic
  integer,parameter:: entroW_=5,diverW_=6,alfvRW_=7,alfvLW_=8 ! waves
  integer,parameter:: nworkroe=15


  
  
  
  INTEGER,PARAMETER:: gamma_=1, adiab_=2, m_=3, a_=4, neqpar=4 !equation params
  
  

  

  INTEGER,PARAMETER:: nflag_=nw+1
  INTEGER:: flags(nflag_)
  DOUBLE PRECISION:: wflags(nflag_)

  DOUBLE PRECISION::minp,minrho,smallxi,smalltau , govergminone
  DOUBLE PRECISION::limitvalue

end module mod_amrvacpar
!=============================================================================
