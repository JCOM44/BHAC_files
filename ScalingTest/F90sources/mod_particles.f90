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
!
! Based on cleaned up version of amrvac.org by:
!> By Jannis Teunissen, Bart Ripperda, Oliver Porth, and Fabio Bacchini (2017-2020)
!=============================================================================
module mod_particles
  use mod_particles_base
  use mod_particles_advect
  use mod_particles_user, only: init_particles_user

  implicit none

  contains
    
    !=============================================================================
    subroutine init_tracerparticles()
      !> Initialize particle data and parameters

      use mod_amrvacdef
      !-----------------------------------------------------------------------------
      call init_particles_vars
      call init_particles_com

      select case(physics_type_particles)
      case('advect')
         call advect_init
      case default
         if (mype == 0) then
            print *, "Unknown physics_type_particles", trim(physics_type_particles)
            print *, "Options are: advect"
            call mpistop("Unknown physics_type_particles")
         end if
      end select


      call init_particles_user
      
    end subroutine init_tracerparticles
    !=============================================================================
  end module mod_particles
  !=============================================================================
  !=============================================================================
