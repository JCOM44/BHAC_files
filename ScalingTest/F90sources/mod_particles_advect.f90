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
!> Tracer for advected particles moving with fluid flows based on a version by
!> By Jannis Teunissen, Bart Ripperda, Oliver Porth, and Fabio Bacchini (2017-2020)
module mod_particles_advect
  use mod_particles_base

  private

  public :: advect_init


contains

  !=============================================================================
  subroutine advect_init()
    !-----------------------------------------------------------------------------

    ngridvars=3
    vp1_=1;vp2_=2;vp3_=3

    fill_gridvars           => fill_gridvars_advect
    integrate_particles     => integrate_particles_advect

  end subroutine advect_init
  !=============================================================================
  subroutine integrate_particles_advect(end_time)
    ! this solves dx/dt=v for particles

    use mod_odeint
    use mod_amrvacdef
    !-----------------------------------------------------------------------------
    double precision, intent(in)        :: end_time
    integer                             :: ipart, iipart
    double precision                    :: dt_p
    double precision, dimension(1:ndir) :: v, x
    integer                             :: igrid
    double precision                    :: rho, rho1, rho2, td, tloc
    !    double precision, dimension(ixG^T,1:nw)   :: w
    ! for odeint:
    !    integer                          :: nok, nbad
    !    double precision                 :: h1
    !    double precision,parameter       :: eps=1.0d-6, hmin=1.0d-8
    !    integer, parameter               :: nvar = ^NC
    external derivs_advect
    !-----------------------------------------------------------------------------

    call set_particles_dt_advect(end_time) ! also outputs individual particles

    do iipart=1,nparticles_active_on_mype
    ipart=particles_active_on_mype(iipart);

       dt_p = particle(ipart)%self%dt
       igrid = particle(ipart)%igrid
       igrid_working = igrid
       tloc = particle(ipart)%self%t
       x(1:ndir) = particle(ipart)%self%x(1:ndir)
       call set_tmpGlobals(igrid)

       ! **************************************************
       ! Position update
       ! **************************************************

       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Simple forward Euler:
       call get_vec(igrid,x,tloc,v,vp1_,vp3_)
       particle(ipart)%self%u(1:ndir) = v(1:ndir)

       particle(ipart)%self%x(1:ndir) = particle(ipart)%self%x(1:ndir) + dt_p &
          * v(1:ndir)
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Adaptive stepwidth RK4:
       !   h1 = dt_p/2.0d0
       !   call odeint(x,nvar,tloc,tloc+dt_p,eps,h1,hmin,nok,nbad,derivs_advect,rkqs)
       !   particle(ipart)%self%x(1:ndir) = x(1:ndir)
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



       ! **************************************************
       ! Velocity update
       ! **************************************************
       call get_vec(igrid,x,tloc+dt_p,v,vp1_,vp3_)
       particle(ipart)%self%u(1:ndir) = v(1:ndir)


       ! **************************************************
       ! Time update
       ! **************************************************
       particle(ipart)%self%t = particle(ipart)%self%t + dt_p


    end do

  end subroutine integrate_particles_advect
  !=============================================================================
  subroutine derivs_advect(t_s,x,dxdt)

    use mod_particles_base, only: get_vec, igrid_working
    use mod_amrvacdef

    double precision                :: t_s, x(*)
    double precision                :: dxdt(*)
    ! .. local ..
    double precision                :: v(1:3)
    !-----------------------------------------------------------------------------

    call get_vec(igrid_working,x(1:3),t_s,v,vp1_,vp3_)
     dxdt(1) = v(1); dxdt(2) = v(2); dxdt(3) = v(3);

  end subroutine derivs_advect
  !=============================================================================
  subroutine set_particles_dt_advect(end_time)

    use mod_particles_base
    use mod_amrvacdef

    double precision, intent(in)    :: end_time
    integer                         :: ipart, iipart
    double precision                :: dt_cfl
    double precision                :: v(1:ndir)
    !-----------------------------------------------------------------------------

    do iipart=1,nparticles_active_on_mype
    ipart=particles_active_on_mype(iipart);

       ! make sure we step only one cell at a time:
       v(1)   = abs(particle(ipart)%self%u(1))
       v(2)   = abs(particle(ipart)%self%u(2))
       v(3)   = abs(particle(ipart)%self%u(3))

       
       dt_cfl = min(rnode(rpdx1_,particle(ipart)%igrid)/v(1),rnode(rpdx2_,&
          particle(ipart)%igrid)/v(2),rnode(rpdx3_,particle(ipart)%igrid)&
          /v(3))
       

       particle(ipart)%self%dt = limit_dt(end_time - particle(ipart)%self%t,&
          dt_cfl*cfl)

       ! **************************************************
       ! Check output for followed particles:
       ! **************************************************
       if (particle(ipart)%self%follow) call check_and_output_particle&
          (particle(ipart)%self,particle(ipart)%payload)

    end do !ipart

  end subroutine set_particles_dt_advect
  !=============================================================================
  subroutine fill_gridvars_advect

    use mod_amrvacdef

    integer                                   :: igrid, iigrid, idir
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
       1:nw)   :: w  
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
       1:ndim) :: x
    !-----------------------------------------------------------------------------

    do iigrid=1,igridstail; igrid=igrids(iigrid);
       call set_tmpGlobals(igrid)

       x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:ndim) &
          = px(igrid)%x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:ndim)

       w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nw)   &
          = pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nw)
       call primitive(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixGlo1,ixGlo2,&
          ixGlo3,ixGhi1,ixGhi2,ixGhi3,w,x)

       gridvars(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
          1:ngridvars) = 0.0d0

       ! Code units fields
       ! fill with Transport-velocity:
       do idir = 1, ndir
          call getv(w,x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixGlo1,&
             ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,idir,gridvars(igrid)%w&
             (ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,vp1_-1+idir))
       end do

       if (time_advance) then

          w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nw)   &
             = pwold(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nw)
          call primitive(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixGlo1,&
             ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,w,x)

          gridvars_old(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
             1:ngridvars) = 0.0d0

          ! Code units fields
          ! fill with Transport-velocity:
          do idir = 1, ndir
             call getv(w,x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixGlo1,&
                ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,idir,gridvars_old(igrid)%w&
                (ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,vp1_-1+idir))
          end do

       end if

    end do

  end subroutine fill_gridvars_advect
  !=============================================================================
end module mod_particles_advect
!=============================================================================
