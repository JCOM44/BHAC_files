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
module mod_particles_user
  use mod_particles_base
  implicit none



contains
  !=============================================================================
  subroutine init_particles_user()

    initial_particles => initial_particles_user
    add_particles     => add_particles_user
    query_destroy     => destroy_user

  end subroutine init_particles_user
  !=============================================================================
  subroutine initial_particles_user()
    ! add the initial particles

    use constants
    use mod_random
    use mod_amrvacdef

    double precision, dimension(ndir)       :: x
    integer                                 :: igrid_particle, ipe_particle
    integer(i8)                             :: seed(2)
    type(rng_t)                             :: myrand
    double precision, dimension(:), allocatable    :: r1,r2,r3
    double precision                        :: v(1:ndir)
    logical, dimension(:),allocatable       :: follow
    type(particle_ptr), dimension(:), allocatable    :: particles_to_inject
    integer                                 :: ninject, nparticles_added,&
        ipart
    !-----------------------------------------------------------------------------

    allocate(follow(nparticles_init), particles_to_inject(nparticles_init))
    allocate(r1(nparticles_init),r2(nparticles_init),r3(nparticles_init))
    do ipart=1,nparticles_init
       allocate(particles_to_inject(ipart)%self)
       allocate(particles_to_inject(ipart)%payload(npayload))
    end do

    ! initialise the random number generator
    seed = [310952_i8,24378_i8]
    call myrand%set_seed(seed)
     call myrand%unif_01_vec(r1)
      call myrand%unif_01_vec(r2)
      call myrand%unif_01_vec(r3)

    ! flags to follow particles
    follow      = .false.
    follow(1:3)  = .true.


    ninject = 0
    nparticles_added = 0
    x=0
    do while (nparticles_added .lt. nparticles_init)
       nparticles_added = nparticles_added + 1

       x(1) = xprobmin1 + r1(nparticles_added) * (xprobmax1 - xprobmin1)
       x(2) = xprobmin2 + r2(nparticles_added) * (xprobmax2 - xprobmin2)
       x(3) = xprobmin3 + r3(nparticles_added) * (xprobmax3 - xprobmin3)

       call find_particle_ipe(x,igrid_particle,ipe_particle)

       if (ipe_particle == mype) then 
          ninject = ninject + 1

          particles_to_inject(ninject)%self%follow = follow(nparticles_added)
          particles_to_inject(ninject)%self%q      = - CONST_e
          particles_to_inject(ninject)%self%m      =   CONST_me

          particles_to_inject(ninject)%self%t      = t
          particles_to_inject(ninject)%self%dt     = 0.0d0

          particles_to_inject(ninject)%self%x(:)   = x

          call get_vec(igrid_particle,x,0.0d0,v,vp1_,vp3_)
          particles_to_inject(ninject)%self%u(:) = v(:)
          particles_to_inject(ninject)%payload(:) = 666

       end if

    end do

    call inject_particles(particles_to_inject(1:ninject))

  end subroutine initial_particles_user
  !=============================================================================
  subroutine add_particles_user()
    ! add some more particles every dt_inject

    use constants
    use mod_random
    use mod_amrvacdef

    double precision, dimension(ndir)       :: x
    integer                                 :: igrid_particle, ipe_particle
    type(rng_t)                             :: myrand
    integer(i8)                             :: seed(2)
    double precision, dimension(:), allocatable    :: r1,r2,r3
    double precision                        :: v(1:ndir)
    type(particle_ptr), dimension(:), allocatable :: particles_to_inject
    integer                                 :: ninject, nparticles_added,&
        ipart
    !-----------------------------------------------------------------------------

    return

    allocate(particles_to_inject(nparticles_inject))
    allocate(r1(nparticles_inject),r2(nparticles_inject),r3&
       (nparticles_inject))
    do ipart=1,nparticles_init
       allocate(particles_to_inject(ipart)%self)
       allocate(particles_to_inject(ipart)%payload(npayload))
    end do

    seed = [310952_i8,24378_i8] + it_particles
    call myrand%set_seed(seed)
     call myrand%unif_01_vec(r1)
      call myrand%unif_01_vec(r2)
      call myrand%unif_01_vec(r3)

    ninject = 0
    nparticles_added = 0
    x=0
    do while (nparticles_added .lt. nparticles_inject)
       nparticles_added = nparticles_added + 1

       x(1) = xprobmin1 + r1(nparticles_added) * (xprobmax1 - xprobmin1)
       x(2) = xprobmin2 + r2(nparticles_added) * (xprobmax2 - xprobmin2)
       x(3) = xprobmin3 + r3(nparticles_added) * (xprobmax3 - xprobmin3)

       call find_particle_ipe(x,igrid_particle,ipe_particle)

       if (ipe_particle == mype) then 
          ninject = ninject + 1

          particles_to_inject(ninject)%self%follow = .false.
          particles_to_inject(ninject)%self%q      = - CONST_e
          particles_to_inject(ninject)%self%m      =   CONST_me

          particles_to_inject(ninject)%self%t      = t_particles
          particles_to_inject(ninject)%self%dt     = 0.0d0

          particles_to_inject(ninject)%self%x(:)   = x

          call get_vec(igrid_particle,x,0.0d0,v,vp1_,vp3_)
          particles_to_inject(ninject)%self%u(:) = v(:)
          particles_to_inject(ninject)%payload(:) = 777

       end if

    end do

    call inject_particles(particles_to_inject(1:ninject))


  end subroutine add_particles_user
  !=============================================================================
  logical function destroy_user(myparticle)

    type(particle_node), intent(in)          :: myparticle
    !-----------------------------------------------------------------------------
    destroy_user = .false.

  end function destroy_user
  !=============================================================================
end module mod_particles_user
!=============================================================================
