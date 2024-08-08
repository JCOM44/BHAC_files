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

!> Module for electroncooling, adopted from Mizuno et al. 2022
module mod_electronphys

  use constants
  implicit none
  
  private

  
  
  ! Set parameter for astrophysical target  
  double precision                   :: mbh_Msun = 4.5d6 !BH mass in solar masses
  double precision                   :: mbh !BH mass in g
  
  double precision                   :: mdot_Edd !Eddington limited accretion rate (mbh)
  double precision                   :: rg !gravitational radius & scale unit
  double precision                   :: tg         ! r_g/c

  ! Set scaling parameter to cgs unit
  double precision                   :: lunit      ! length  unit
  double precision                   :: tunit      ! time unit

  double precision                   :: mdot_Edd_frac = 1.0d-5
  double precision                   :: munit !Define M_unit to scale mass accretion rate target flux (M_sun/year) or Eddington fraction
  double precision                   :: Mdot_Msun_year = 1.d-9 !Obtain M_unit from acc. rate

  double precision                   :: rhounit !Scale cgs density
  double precision                   :: bunit               ! Scale B-field
  double precision                   :: nunit !Electron number density scaling
  double precision                   :: uunit !pressure (internal energy) scaling
  double precision                   :: coolunit

  character*31                       :: typeelheat = 'turbulent'
  character*31                       :: typeelcool = 'allcooling'
  character*31                       :: typemscale = 'Medd'

  double precision, parameter        :: teunit = CONST_me*CONST_c*CONST_c&
     /CONST_kB !m_e c2,c2,c2 / k_B
  double precision, parameter        :: tpunit = CONST_mp*CONST_c*CONST_c&
     /CONST_kB !m_p c2,c2,c2 / k_B 
  double precision, parameter        :: mpme = 1836.15267389
  
contains
  

end module mod_electronphys
