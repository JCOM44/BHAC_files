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
! Assemble quantities used for radiation postprocessing
!
!=============================================================================
subroutine convert_postrad(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,igrid,nwpost,w,x,wpost)

  use mod_physaux
  use mod_amrvacdef
  
  integer, intent(in)              :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, igrid, nwpost
  double precision, intent(in)     :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
     1:ndim)
  double precision, intent(out)    :: wpost(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:nwpost)
  ! .. local ..
  double precision                 :: wprim(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:nw)
  integer, parameter               :: ir_=1, itheta_=2, iphi_=3, irho_&
     =iphi_+1, ivr_=irho_+1, ivtheta_=ivr_+1, ivphi_ = ivtheta_+1, ip_&
     =ivphi_+1, iBr_=ip_+1, iBtheta_=iBr_+1, iBphi_ = iBtheta_+1
  !-----------------------------------------------------------------------------

  if (ibphi_ .gt. nwpost) call mpistop&
     ('convert_postrad: nwpost array too small?!')
  
  wpost = zero
  
  ! First convert to primitive variables:
  wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:nw) &
     = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:nw)
  call primitive(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wprim,x)


  ! Coordinates:
  wpost(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ir_)     &
     = x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,r_)
  
  wpost(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,itheta_) &
     = x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,z_)
 
  wpost(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iphi_)   &
     = x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,phi_)
 

  ! density and pressure:
  wpost(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,irho_) &
     = wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
  wpost(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ip_) &
     = wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,pp_)
  
  ! Fluid three-velocity:
  wpost(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ivr_)     &
     = wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,u1_)&
     /wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)
  
  wpost(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ivtheta_) &
     = wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,u3_)&
     /wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)
 
  wpost(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ivphi_)   &
     = wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,u2_)&
     /wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)
 

  ! Magnetic three-field in the Eulerian frame:
  wpost(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iBr_)     &
     = wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b1_)
  
  wpost(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iBtheta_) &
     = wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b3_)
 
  wpost(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iBphi_)   &
     = wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b2_)
 

end subroutine convert_postrad
!=============================================================================
subroutine convert_postrad_el(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,igrid,nwpost,w,x,wpost)
! prepration of quantities for postrad output in electron thermodynamics case

  use mod_physaux
  use mod_amrvacdef
  
  integer, intent(in)              :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, igrid, nwpost
  double precision, intent(in)     :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
     1:ndim)
  double precision, intent(out)    :: wpost(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:nwpost)
  ! .. local ..
  double precision                 :: wprim(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:nw)
  integer, parameter               :: ir_=1, itheta_=2, iphi_=3, irho_&
     =iphi_+1, ivr_=irho_+1, ivtheta_=ivr_+1, ivphi_ = ivtheta_+1, ip_&
     =ivphi_+1, iBr_=ip_+1, iBtheta_=iBr_+1, iBphi_ = iBtheta_+1
  
  !-----------------------------------------------------------------------------

  if (ibphi_ .gt. nwpost) call mpistop&
     ('convert_postrad: nwpost array too small?!')
  
  wpost = zero
  
  ! First convert to primitive variables:
  wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:nw) &
     = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:nw)
  call primitive(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wprim,x)


  ! Coordinates:
  wpost(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ir_)     &
     = x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,r_)
  
  wpost(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,itheta_) &
     = x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,z_)
 
  wpost(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iphi_)   &
     = x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,phi_)
 

  ! density and pressure:
  wpost(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,irho_) &
     = wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
  wpost(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ip_) &
     = wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,pp_)
  
  ! Fluid three-velocity:
  wpost(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ivr_)     &
     = wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,u1_)&
     /wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)
  
  wpost(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ivtheta_) &
     = wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,u3_)&
     /wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)
 
  wpost(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ivphi_)   &
     = wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,u2_)&
     /wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)
 

  ! Magnetic three-field in the Eulerian frame:
  wpost(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iBr_)     &
     = wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b1_)
  
  wpost(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iBtheta_) &
     = wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b3_)
 
  wpost(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iBphi_)   &
     = wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b2_)
 

  
  
end subroutine convert_postrad_el
!=============================================================================
