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
module mod_particles_base
  !
  ! Contains global variables and routines for particles
  !
  use mod_physicaldata
  implicit none

  !> String describing the particle physics type
  character(len=128)               :: physics_type_particles = ""
  !> String describing the type of ensemble (and destroy) IO
  character(len=40)                :: ensemble_type = ""
  !> String describing the particle integrator type
  character(len=128)               :: integrator_type_particles = ""
  integer                          :: nparticleshi, nparticles_per_cpu_hi
  !> Number of total payload variables for a particle
  integer                          :: npayload
  integer                          :: nparticles_init, nparticles_inject
  double precision                 :: t_particles !time of the slowest particle
  double precision                 :: dtinject_particles

  double precision                 :: tmax_particles, tinject_particles
  double precision                 :: t_next_injection, t_next_output
  integer                          :: itmax_particles
  double precision                 :: dtsave_ensemble, dtsave_individual
  !> Output count for ensembles
  integer                          :: n_output_ensemble
  double precision                 :: dtheta
  double precision                 :: cfl
  logical                          :: losses
  logical                          :: write_ensemble, write_individual
  logical                          :: inject
  logical                          :: nocartesian !Don't convert particle position to cartesian

  !> Header string used in CSV files
  character(len=400)               :: csv_header

  integer                          :: nparticles !How many particles we have globally atm
  integer                          :: n_active !How many active particles globally atm
  integer                          :: it_particles !How many iterations were don in the main loop
  integer                          :: index_free !next free index for a particle

  ! these two save the list of neighboring cpus:
  integer, dimension(:), allocatable     :: ipe_neighbor
  integer                                :: npe_neighbors

  integer                                :: type_particle

  integer, parameter                     :: unitparticles=15 

  ! list of particles on the processor:
  integer, dimension(:), allocatable     :: particles_on_mype,&
      particles_active_on_mype
  integer                                :: nparticles_on_mype,&
      nparticles_active_on_mype

  type(walloc),dimension(:), allocatable :: gridvars
  type(walloc),dimension(:), allocatable :: gridvars_old
  integer                                :: ngridvars


  ! I use this variables to set the current igrid and ipart for the particle integrator:
  integer                          :: igrid_working, ipart_working
  !$OMP THREADPRIVATE(igrid_working,ipart_working)

  ! Variable indices for the gridvars array
  integer                                :: vp1_,vp2_,vp3_
  integer                                :: bp1_,bp2_,bp3_
  integer                                :: ep1_,ep2_,ep3_

  ! TODO: this just serves as a reminder to make the corresponding mod_particles_xxx:
  
  

  type particle_node
     logical                               :: follow
     integer                               :: index
     double precision                      :: q, m
     double precision                      :: t, dt
     double precision, dimension(3)      :: x
     double precision, dimension(3)      :: u
  end type particle_node

  type particle_ptr
     type(particle_node), allocatable     :: self
     !> extra information carried by the particle
     double precision, allocatable        :: payload(:)
     integer                              :: igrid, ipe
  end type particle_ptr

  ! Array containing all particles
  type(particle_ptr), dimension(:), allocatable   :: particle

  ! Procedure pointers for various implementations
  procedure(output_ensemble_csv), pointer          :: output_ensemble &
     => null()

  ! Pointers for user-defined stuff
  procedure(sub_fill_gridvars), pointer            :: fill_gridvars => null()
  procedure(sub_integrate), pointer                :: integrate_particles &
     => null()
  procedure(sub_dummy), pointer                    :: initial_particles &
     => null()
  procedure(sub_dummy), pointer                    :: add_particles => null()
  procedure(fun_user_destroy), pointer             :: query_destroy => null()



  !INTERFACES
  abstract interface 

     subroutine sub_fill_gridvars
     end subroutine sub_fill_gridvars

     subroutine sub_integrate(end_time)
       double precision, intent(in) :: end_time
     end subroutine sub_integrate

     subroutine sub_dummy
     end subroutine sub_dummy

     function fun_user_destroy(myparticle)
       import particle_node
       logical                         :: fun_user_destroy
       type(particle_node), intent(in) :: myparticle
     end function fun_user_destroy

  end interface

  !=============================================================================
contains
  !=============================================================================
  subroutine particles_params_read

    use mod_amrvacdef
    !-----------------------------------------------------------------------------

    namelist /particleslist/ physics_type_particles,nparticleshi,&
        nparticles_per_cpu_hi, write_individual, write_ensemble, dtheta,&
        losses, integrator_type_particles, inject, dtsave_ensemble,&
        dtsave_individual, nparticles_init, nparticles_inject, npayload,&
        nocartesian, cfl, ensemble_type

    open(unitpar, file=inifile, status="old")
    read(unitpar, particleslist, end=111)
111 close(unitpar)

  end subroutine particles_params_read
  !=============================================================================
  subroutine init_particles_vars()

    use mod_amrvacdef
    character(len=20)  :: strdata
    integer            :: icomp
    !-----------------------------------------------------------------------------

    ! initialise 
    physics_type_particles    = 'advect'
    integrator_type_particles = 'Boris'
    ensemble_type             = 'csv'
    nparticleshi              = 10000000
    nparticles_per_cpu_hi     = 50000
    nparticles_init           = 100
    nparticles_inject         = 10
    npayload                  = 1
    dtinject_particles        = bigdouble
    t_particles               = 0.0d0
    tmax_particles            = bigdouble
    tinject_particles         = bigdouble
    t_next_injection          = 0.0d0
    t_next_output             = 0.0d0
    itmax_particles           = biginteger
    dtheta                    = 2.0d0*dpi / 60.0d0
    cfl                       = 0.5d0
    write_ensemble            = .false.
    write_individual          = .false.
    dtsave_ensemble           = bigdouble
    dtsave_individual         = bigdouble
    n_output_ensemble         = 0

    inject                    = .false.
    losses                    = .false.
    nocartesian               = .false.

    call particles_params_read

    nparticles                = 0
    n_active                  = 0
    it_particles              = 0
    index_free                = 1


    allocate(particle(1:nparticleshi), particles_on_mype&
       (nparticles_per_cpu_hi), particles_active_on_mype&
       (nparticles_per_cpu_hi))
    allocate(ipe_neighbor(0:npe-1))

    allocate(gridvars(ngridshi))
    if (.not. convert) allocate(gridvars_old(ngridshi))


    nparticles_on_mype   = 0
    particles_on_mype(:) = 0

    nparticles_active_on_mype   = 0
    particles_active_on_mype(:) = 0

    csv_header = "# t, dt, x1, x2, x3, u1, u2, u3,"
    do icomp=1,npayload
       write(strdata,"(a,i0.2,a)") 'payload',icomp,','
       csv_header = trim(csv_header) // ' ' // trim(strdata)
    end do
    csv_header = trim(csv_header) // ' ipe, iteration, index'

    select case(ensemble_type)
    case('csv')
       output_ensemble => output_ensemble_csv
    case('dat')
       output_ensemble => output_ensemble_dat
    case default
       call mpistop('init_particles_vars: unknown ensemble_type')
    end select

  end subroutine init_particles_vars
  !=============================================================================
  logical function exit_condition()
    ! check if we should go out of the integration loop

    !-----------------------------------------------------------------------------

    exit_condition = (n_active == 0 .or. nparticles == 0 .or. it_particles &
       >= itmax_particles )

  end function exit_condition
  !=============================================================================
  double precision function limit_dt(t_left,dt_p)
    ! Make sure we don't miss an output or tmax_particles

    double precision, intent(in) :: t_left,dt_p
    ! .. local ..
    double precision                :: n_steps
    double precision, parameter     :: eps = 1d-10
    !-----------------------------------------------------------------------------

    n_steps = t_left / dt_p

    if (n_steps < 1.0d0+eps) then
       ! Last step
       limit_dt = t_left
    else if (n_steps < 2.0d0-eps) then
       ! Divide time left over two steps, to prevent one tiny time step
       limit_dt = 0.5d0 * t_left
    else
       limit_dt = dt_p
    end if

    !    if (limit_dt < 1.0d-4) print*, 'very small step used!', limit_dt, t_left, dt_p

  end function limit_dt
  !=============================================================================
  subroutine get_new_indices(ninject,new_indices,ninject_total)
    ! get new unique indices for ninject particles (synchronized across MPI tasks)

    use mod_amrvacdef
    integer, intent(in)                      :: ninject !how many particles
    ! this task wants to inject
    integer, intent(out), dimension(ninject) :: new_indices !the new uniqe indices
    integer, intent(out)                     :: ninject_total !add over all processors
    ! .. local ..
    integer                                  :: n_particles_to_inject(0:npe-1)
    integer                                  :: ipe, inew, ninjected
    !-----------------------------------------------------------------------------

    if (npe>0) then
       call MPI_ALLGATHER(ninject,1,MPI_INTEGER,n_particles_to_inject,1,&
          MPI_INTEGER,icomm,ierrmpi)
    else
       n_particles_to_inject(0) = ninject
    end if

    ! go through the processors and fill up the index array:
    ninjected = 0
    ninject_total = 0
    do ipe=0,npe-1
       do inew=1,n_particles_to_inject(ipe)
          if (ipe .eq. mype) then
             ninjected = ninjected + 1
             new_indices(ninjected) = index_free
          end if
          index_free = index_free + 1
          ninject_total = ninject_total + 1
       end do
    end do

  end subroutine get_new_indices
  !=============================================================================
  subroutine inject_particles(particles_inject)
    ! inject ninject new particles on each process

    use mod_timing
    use mod_amrvacdef
    !    include 'mpif.h'

    type(particle_ptr), dimension(:), intent(in)         :: particles_inject
    ! .. local ..
    integer                                              :: ninject
    integer                                              :: ipart, index,&
        ipe_particle, igrid_particle
    integer                                              :: ninject_total !How many particles across all processses
    integer, allocatable, dimension(:)                   :: indices
    !-----------------------------------------------------------------------------

    ninject = size(particles_inject)

    ! first get new indices:
    allocate(indices(1:ninject))
    call get_new_indices(ninject,indices,ninject_total)

    ! update total number of particles:
    nparticles = nparticles + ninject_total

    ! allocate and copy data:
    do ipart=1,ninject
       index = indices(ipart)

       call push_particle_into_particles_on_mype(index)

       if (.not. allocated(particle(index)%self)) allocate(particle&
          (index)%self)
       if (.not. allocated(particle(index)%payload)) allocate(particle&
          (index)%payload(npayload))
       particle(index)%self = particles_inject(ipart)%self
       particle(index)%payload = particles_inject(ipart)%payload
       particle(index)%self%index = index

       ! This is probably overkill:
       call find_particle_ipe(particle(index)%self%x,igrid_particle,&
          ipe_particle)
       particle(index)%ipe = ipe_particle
       particle(index)%igrid = igrid_particle

    end do

    ! to be safe, communicate the particles to where they actually belong:
    tpartc_com0=MPI_WTIME()
    call comm_particles_global()
    tpartc_com=tpartc_com + (MPI_WTIME()-tpartc_com0)

  end subroutine inject_particles
  !=============================================================================
  subroutine finish_particles_vars()
    ! de-allocate some dynamic arrays

    use mod_amrvacdef
    !-----------------------------------------------------------------------------

    deallocate(particle)
    deallocate(particles_on_mype)
    deallocate(particles_active_on_mype)
    deallocate(ipe_neighbor)
    deallocate(gridvars)
    if (time_advance) deallocate(gridvars_old)

  end subroutine finish_particles_vars
  !=============================================================================
  subroutine select_active_particles(end_time)

    use mod_amrvacdef
    double precision, intent(in)                    :: end_time
    integer                                         :: ipart, iipart
    logical                                         :: activate
    double precision                                :: t_min_mype
    !-----------------------------------------------------------------------------

    t_min_mype = bigdouble
    nparticles_active_on_mype = 0
    !$OMP PARALLEL DO PRIVATE(ipart)
    do iipart=1,nparticles_on_mype;ipart=particles_on_mype(iipart);

       activate = particle(ipart)%self%t .lt. end_time
       t_min_mype = min(t_min_mype, particle(ipart)%self%t)

       if (activate) then
          !$OMP CRITICAL(enroll)
          nparticles_active_on_mype = nparticles_active_on_mype + 1
          particles_active_on_mype(nparticles_active_on_mype) = ipart
          !$OMP END CRITICAL(enroll)
       end if

    end do
    !$OMP END PARALLEL DO

    ! Determine total number of active particles
    call MPI_ALLREDUCE(nparticles_active_on_mype, n_active, 1, MPI_INTEGER,&
        MPI_SUM, icomm, ierrmpi)

    ! Determine time of last particle
    call MPI_ALLREDUCE(t_min_mype, t_particles, 1, MPI_DOUBLE_PRECISION,&
        MPI_MIN, icomm, ierrmpi)

  end subroutine select_active_particles
  !=============================================================================
  subroutine locate_particle(index,igrid_particle,ipe_particle)
    ! given the particles unique index, tell me on which cpu and igrid it is
    ! returns -1,-1 if particle was not found
    use mod_amrvacdef
    integer, intent(in)                            :: index
    integer, intent(out)                           :: igrid_particle,&
        ipe_particle
    integer                                        :: iipart,ipart,&
       ipe_has_particle,ipe
    logical                                        :: has_particle(0:npe-1)
    integer,dimension(0:1)                         :: buff
    !-----------------------------------------------------------------------------

    has_particle(:) = .false.
    do iipart=1,nparticles_on_mype;ipart=particles_on_mype(iipart);
       if (particle(ipart)%self%index == index) then
          has_particle(mype) = .true.
          exit
       end if
    end do

    if (has_particle(mype)) then
       buff(0) = particle(ipart)%igrid
       buff(1) = mype
    end if

    if (npe>0) call MPI_ALLGATHER(MPI_IN_PLACE,1,MPI_LOGICAL,has_particle,1,&
       MPI_LOGICAL,icomm,ierrmpi)

    ipe_has_particle = -1
    do ipe=0,npe-1
       if (has_particle(ipe) .eqv. .true.) then 
          ipe_has_particle = ipe
          exit
       end if
    end do

    if (ipe_has_particle .ne. -1) then
       if (npe>0) call MPI_BCAST(buff,2,MPI_INTEGER,ipe_has_particle,icomm,&
          ierrmpi)
       igrid_particle=buff(0)
       ipe_particle = buff(1)
    else
       igrid_particle=-1
       ipe_particle = -1
    end if

  end subroutine locate_particle
  !=============================================================================
  subroutine find_particle_ipe(x,igrid_particle,ipe_particle)

    use mod_forest, only: tree_node_ptr, tree_root
    use mod_slice, only: get_igslice
    use mod_amrvacdef

    double precision, dimension(ndir), intent(inout)   :: x
    integer, intent(out)                            :: igrid_particle,&
        ipe_particle

    integer, dimension(ndir,nlevelshi)              :: ig
    integer                                         :: idim, ic(ndim)
    type(tree_node_ptr)                             :: branch
    !-----------------------------------------------------------------------------

    ! first check if the particle is in the domain
    if (.not. particle_in_domain(x)) then
       igrid_particle = -1
       ipe_particle   = -1
       return
    end if

    ! get the index on each level
    do idim = 1, ndim
       call get_igslice(idim,x(idim),ig(idim,:))
    end do

    ! traverse the tree until leaf is found
    branch=tree_root(ig(1,1),ig(2,1),ig(3,1))
    do while (.not.branch%node%leaf)
       ic(1)=ig(1,branch%node%level+1) - 2 * branch%node%ig1 +2
       ic(2)=ig(2,branch%node%level+1) - 2 * branch%node%ig2 +2
       ic(3)=ig(3,branch%node%level+1) - 2 * branch%node%ig3 +2
       branch%node => branch%node%child(ic(1),ic(2),ic(3))%node
    end do

    igrid_particle = branch%node%igrid
    ipe_particle   = branch%node%ipe

  end subroutine find_particle_ipe
  !=============================================================================
  logical function particle_in_domain(x)

    use mod_amrvacdef
    double precision, dimension(ndim), intent(in)  :: x
    integer                                        :: idim
    !-----------------------------------------------------------------------------

    particle_in_domain = .true.

    do idim=1,ndim
       select case(idim)
          case (1)
          if (x(1) .lt. xprobmin1) then
             particle_in_domain = .false.
             exit
          end if
          if (x(1) .ge. xprobmax1) then
             particle_in_domain = .false.
             exit
          end if
          
          case (2)
          if (x(2) .lt. xprobmin2) then
             particle_in_domain = .false.
             exit
          end if
          if (x(2) .ge. xprobmax2) then
             particle_in_domain = .false.
             exit
          end if
          
          case (3)
          if (x(3) .lt. xprobmin3) then
             particle_in_domain = .false.
             exit
          end if
          if (x(3) .ge. xprobmax3) then
             particle_in_domain = .false.
             exit
          end if
          
       end select
    end do

  end function particle_in_domain
  !=============================================================================
  logical function particle_in_igrid(ipart,igrid)
    ! quick check if particle is still in igrid

    use mod_amrvacdef
    integer, intent(in)                            :: igrid,ipart
    integer                                        :: idim
    !-----------------------------------------------------------------------------

    ! first check if the igrid is still there:
    if (.not. associated(pw(igrid)%w)) then
       particle_in_igrid = .false.
       return
    end if

    particle_in_igrid = .true.

    do idim=1,ndim
       select case(idim)
          case (1)
          if (particle(ipart)%self%x(1) .lt. rnode(rpxmin1_,igrid)) then
             particle_in_igrid = .false.
             exit
          end if
          if (particle(ipart)%self%x(1) .ge. rnode(rpxmax1_,igrid)) then
             particle_in_igrid = .false.
             exit
          end if
          
          case (2)
          if (particle(ipart)%self%x(2) .lt. rnode(rpxmin2_,igrid)) then
             particle_in_igrid = .false.
             exit
          end if
          if (particle(ipart)%self%x(2) .ge. rnode(rpxmax2_,igrid)) then
             particle_in_igrid = .false.
             exit
          end if
          
          case (3)
          if (particle(ipart)%self%x(3) .lt. rnode(rpxmin3_,igrid)) then
             particle_in_igrid = .false.
             exit
          end if
          if (particle(ipart)%self%x(3) .ge. rnode(rpxmax3_,igrid)) then
             particle_in_igrid = .false.
             exit
          end if
          
       end select
    end do

    return
  end function particle_in_igrid
  !=============================================================================
  subroutine init_particles_com()
    ! initialise communicators for particles
    use mod_amrvacdef

    integer(kind=MPI_ADDRESS_KIND)         :: size_int, size_double,&
        size_logical, lb
    integer                                :: oldtypes(0:7), offsets(0:7),&
        blocklengths(0:7)
    !-----------------------------------------------------------------------------

    ! create the MPI-datatype for particle

    call MPI_TYPE_GET_EXTENT(MPI_INTEGER, lb, size_int, ierrmpi)
    call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION, lb, size_double, ierrmpi)
    call MPI_TYPE_GET_EXTENT(MPI_LOGICAL, lb, size_logical, ierrmpi)

    oldtypes(0) = MPI_LOGICAL
    oldtypes(1) = MPI_INTEGER
    oldtypes(2:7) = MPI_DOUBLE_PRECISION

    blocklengths(0:5)=1
    blocklengths(6)=3
    blocklengths(7)=3

    offsets(0) = 0
    offsets(1) = size_logical * blocklengths(0)
    offsets(2) = offsets(1) + size_int * blocklengths(1)
    offsets(3) = offsets(2) + size_double * blocklengths(2)
    offsets(4) = offsets(3) + size_double * blocklengths(3)
    offsets(5) = offsets(4) + size_double * blocklengths(4)
    offsets(6) = offsets(5) + size_double * blocklengths(5)
    offsets(7) = offsets(6) + size_double * blocklengths(6)

    call MPI_TYPE_STRUCT(8,blocklengths,offsets,oldtypes,type_particle,&
       ierrmpi)
    call MPI_TYPE_COMMIT(type_particle,ierrmpi)

  end subroutine init_particles_com
  !=============================================================================
  subroutine finish_particles_com()

    use mod_amrvacdef

    call MPI_TYPE_FREE(type_particle,ierrmpi)

  end subroutine finish_particles_com
  !=============================================================================
  subroutine finish_particles()

    integer             :: ipart, iipart
    !-----------------------------------------------------------------------------

    call destroy_particles(nparticles_on_mype,particles_on_mype&
       (1:nparticles_on_mype))

  end subroutine finish_particles
  !=============================================================================
  subroutine set_neighbor_ipe()

    use mod_amrvacdef
    integer              :: igrid, iigrid,ipe
    logical              :: ipe_is_neighbor(0:npe-1)
    integer              :: my_neighbor_type, i1,i2,i3
    !-----------------------------------------------------------------------------

    ipe_is_neighbor(:) = .false.
    do iigrid=1,igridstail; igrid=igrids(iigrid);

       do i3=-1,1
       do i2=-1,1
       do i1=-1,1
       if (i1==0.and.i2==0.and.i3==0) then
          cycle
       end if
       my_neighbor_type=neighbor_type(i1,i2,i3,igrid)

       select case (my_neighbor_type)
       case (1) ! boundary
          ! do nothing
       case (2) ! fine-coarse
          call ipe_fc(i1,i2,i3,igrid,ipe_is_neighbor)
       case (3) ! same level
          call ipe_srl(i1,i2,i3,igrid,ipe_is_neighbor)
       case (4) ! coarse-fine
          call ipe_cf(i1,i2,i3,igrid,ipe_is_neighbor)
       end select

       end do
       end do
       end do

    end do

    ! remove self as neighbor
    ipe_is_neighbor(mype) = .false.

    ! Now make the list of neighbors
    npe_neighbors = 0
    do ipe=0,npe-1
       if (ipe_is_neighbor(ipe)) then
          npe_neighbors = npe_neighbors + 1
          ipe_neighbor(npe_neighbors) = ipe
       end if
    end do ! ipe
    if (npe_neighbors .gt. 32) call mpistop("set_neighbor_ipe: more than 32 &
       neighbors, how is that possible?")

  end subroutine set_neighbor_ipe
  !=============================================================================
  subroutine ipe_fc(i1,i2,i3,igrid,ipe_is_neighbor)

    use mod_amrvacdef
    integer, intent(in) :: i1,i2,i3, igrid
    logical, intent(inout) :: ipe_is_neighbor(0:npe-1)
    !-----------------------------------------------------------------------------

    ipe_is_neighbor(neighbor(2,i1,i2,i3,igrid)) = .true.

  end subroutine ipe_fc
  !=============================================================================
  subroutine ipe_srl(i1,i2,i3,igrid,ipe_is_neighbor)

    use mod_amrvacdef
    integer, intent(in) :: i1,i2,i3, igrid
    logical, intent(inout) :: ipe_is_neighbor(0:npe-1)
    !-----------------------------------------------------------------------------

    ipe_is_neighbor(neighbor(2,i1,i2,i3,igrid)) = .true.

  end subroutine ipe_srl
  !=============================================================================
  subroutine ipe_cf(i1,i2,i3,igrid,ipe_is_neighbor)

    use mod_amrvacdef
    integer, intent(in)    :: i1,i2,i3, igrid
    logical, intent(inout) :: ipe_is_neighbor(0:npe-1)
    integer                :: ic1,ic2,ic3, inc1,inc2,inc3
    !-----------------------------------------------------------------------------

    do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
    inc3=2*i3+ic3
    
    do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
    inc2=2*i2+ic2
    
    do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
    inc1=2*i1+ic1
    
    ipe_is_neighbor( neighbor_child(2,inc1,inc2,inc3,igrid) ) = .true.

    end do
    end do
    end do

  end subroutine ipe_cf
  !=============================================================================
  subroutine check_and_output_particle(myparticle,payload)

    use mod_amrvacdef

    type(particle_node), intent(in) :: myparticle
    double precision, intent(in)    :: payload(npayload)
    ! .. local ..
    integer                         :: nout, nout_next
    logical                         :: file_exists
    character(len=128)              :: filename
    !-----------------------------------------------------------------------------

    nout      = int(myparticle%t/dtsave_individual)
    nout_next = int((myparticle%t + myparticle%dt)/dtsave_individual)

    if (nout_next - nout .le. 0) return  ! get out if nothing to do

    !$OMP CRITICAL(IO)
    ! Create file and write header
    filename = make_particle_filename(0,myparticle%index,'individual')
    inquire(file=filename, exist=file_exists)
    if (.not. file_exists) then
       open(unit=unitparticles, file=filename)
       write(unitparticles,"(a)") trim(csv_header)
    else
       open(unit=unitparticles, file=filename, access='append')
    end if

    call output_particle(myparticle,payload,mype,unitparticles)

    close(unitparticles)
    !$OMP END CRITICAL(IO)

  end subroutine check_and_output_particle
  !=============================================================================
  subroutine particles_output()

    use mod_amrvacdef

    integer                         :: ipart,iipart
    type(particle_node), dimension(nparticles_per_cpu_hi)  :: send_particles
    double precision, dimension(npayload,nparticles_per_cpu_hi)  :: &
       send_payload
    integer                         :: send_n_particles_for_output
    !-----------------------------------------------------------------------------

    send_n_particles_for_output = 0
    do iipart=1,nparticles_on_mype;ipart=particles_on_mype(iipart);

       ! have to send particle to rank zero for output
       send_n_particles_for_output = send_n_particles_for_output + 1
       send_particles(send_n_particles_for_output) = particle(ipart)%self
       send_payload(1:npayload,send_n_particles_for_output) &
          = particle(ipart)%payload(1:npayload)

    end do

    call output_ensemble(send_n_particles_for_output,send_particles,&
       send_payload,'ensemble')

  end subroutine particles_output
  !=============================================================================
  character(len=128) function make_particle_filename(nout,index,type)

    use mod_amrvacdef

    character(len=*), intent(in)    :: type
    integer, intent(in)             :: nout
    integer, intent(in)             :: index
    integer                         :: mysnapshot
    !-----------------------------------------------------------------------------

    if (snapshotini .ne. -1) then
       mysnapshot = snapshotini
    else
       mysnapshot = 0
    end if

    select case(type)
    case ('ensemble')
       write(make_particle_filename,"(a,i4.4,a,i6.6)") trim(filenameout),&
          mysnapshot,'_ensemble',nout
    case ('destroy')
       write(make_particle_filename,"(a,i4.4,a)") trim(filenameout),&
          mysnapshot,'_destroyed'
    case ('individual')
       write(make_particle_filename,"(a,i4.4,a,i6.6,a)") trim(filenameout),&
          mysnapshot,'_particle',index,'.csv'
    end select

  end function make_particle_filename
  !=============================================================================
  subroutine output_ensemble_csv(send_n_particles_for_output,send_particles,&
     send_payload,typefile)

    use mod_amrvacdef

    integer, intent(in)             :: send_n_particles_for_output
    type(particle_node), dimension(send_n_particles_for_output),&
        intent(in)  :: send_particles
    double precision, dimension(npayload,send_n_particles_for_output),&
        intent(in)  :: send_payload
    character(len=*), intent(in)    :: typefile
    ! .. local ..
    character(len=128)              :: filename
    type(particle_node), dimension(nparticles_per_cpu_hi)  :: &
       receive_particles
    double precision, dimension(npayload,nparticles_per_cpu_hi) :: &
       receive_payload
    integer                         :: status(MPI_STATUS_SIZE)
    integer,dimension(0:npe-1)      :: receive_n_particles_for_output_from_ipe
    integer                         :: ipe, ipart
    logical                         :: file_exists
    !-----------------------------------------------------------------------------
    receive_n_particles_for_output_from_ipe(:) = 0

    call MPI_ALLGATHER(send_n_particles_for_output, 1, MPI_INTEGER,&
        receive_n_particles_for_output_from_ipe, 1, MPI_INTEGER, icomm,&
        ierrmpi)

    ! If there are no particles to be written, skip the output
    if (sum(receive_n_particles_for_output_from_ipe(:)) == 0) return

    if (mype > 0) then
       call MPI_SEND(send_particles,send_n_particles_for_output,type_particle,&
          0,mype,icomm,ierrmpi)
       call MPI_SEND(send_payload,npayload*send_n_particles_for_output,&
          MPI_DOUBLE_PRECISION,0,mype,icomm,ierrmpi)
    else
       ! Create file and write header
       filename = make_particle_filename(n_output_ensemble,0,typefile)
       filename = trim(filename)//'.csv'
       if (typefile == 'ensemble') n_output_ensemble = n_output_ensemble + 1
       inquire(file=filename, exist=file_exists)
       if (.not. file_exists) then
          open(unit=unitparticles, file=filename)
          write(unitparticles,"(a)") trim(csv_header)
       else
          open(unit=unitparticles, file=filename, access='append')
       end if

       ! Write own particles
       do ipart=1,send_n_particles_for_output
          call output_particle(send_particles(ipart),send_payload(1:npayload,&
             ipart),0,unitparticles)
       end do

       ! Write particles from other tasks
       do ipe=1,npe-1
          call MPI_RECV(receive_particles,&
             receive_n_particles_for_output_from_ipe(ipe),type_particle,ipe,&
             ipe,icomm,status,ierrmpi)
          call MPI_RECV(receive_payload,&
             npayload*receive_n_particles_for_output_from_ipe(ipe),&
             MPI_DOUBLE_PRECISION,ipe,ipe,icomm,status,ierrmpi)
          do ipart=1,receive_n_particles_for_output_from_ipe(ipe)
             call output_particle(receive_particles(ipart),&
                receive_payload(1:npayload,ipart),ipe,unitparticles)
          end do ! ipart
       end do ! ipe

       close(unitparticles)
    end if

  end subroutine output_ensemble_csv
  !=============================================================================
  subroutine output_ensemble_dat(send_n_particles_for_output,send_particles,&
     send_payload,typefile)

    use mod_amrvacdef

    integer, intent(in)             :: send_n_particles_for_output
    type(particle_node), dimension(send_n_particles_for_output),&
        intent(in)  :: send_particles
    double precision, dimension(npayload,send_n_particles_for_output),&
        intent(in)  :: send_payload
    character(len=*), intent(in)    :: typefile
    ! .. local ..
    character(len=128)              :: filename
    type(particle_node), dimension(nparticles_per_cpu_hi)  :: &
       receive_particles
    double precision, dimension(npayload,nparticles_per_cpu_hi) :: &
       receive_payload
    integer                         :: status(MPI_STATUS_SIZE)
    integer,dimension(0:npe-1)      :: receive_n_particles_for_output_from_ipe
    integer                         :: ipe, ipart
    logical                         :: file_exists
    !-----------------------------------------------------------------------------
    receive_n_particles_for_output_from_ipe(:) = 0

    call MPI_ALLGATHER(send_n_particles_for_output, 1, MPI_INTEGER,&
        receive_n_particles_for_output_from_ipe, 1, MPI_INTEGER, icomm,&
        ierrmpi)

    ! If there are no particles to be written, skip the output
    if (sum(receive_n_particles_for_output_from_ipe(:)) == 0) return

    if (mype > 0) then
       call MPI_SEND(send_particles,send_n_particles_for_output,type_particle,&
          0,mype,icomm,ierrmpi)
       call MPI_SEND(send_payload,npayload*send_n_particles_for_output,&
          MPI_DOUBLE_PRECISION,0,mype,icomm,ierrmpi)
    else
       ! Create file and write header
       if (mype .eq. 0) then
          filename = make_particle_filename(n_output_ensemble,0,typefile)
          filename = trim(filename)//'.dat'
          if (typefile == 'ensemble') n_output_ensemble = n_output_ensemble + &
             1
          INQUIRE(FILE=filename, EXIST=file_exists)
          if (.not. file_exists) then
             open(unit=unitparticles,file=filename,form='unformatted',status&
                ='new',access='stream')
             write(unitparticles) nparticles,it_particles,npayload
             ! FIXME: ^ this will lead to a bug if destroy and injection is used (as nparticles when the file is written need not correspond to the actual number in the file)
          else
             open(unit=unitparticles,file=filename,form='unformatted',status&
                ='old',access='stream',position='append')
          end if
       end if

       ! Write own particles
       do ipart=1,send_n_particles_for_output
          call append_to_snapshot(send_particles(ipart),send_payload&
             (1:npayload,ipart))
       end do

       ! Write particles from other tasks
       do ipe=1,npe-1
          call MPI_RECV(receive_particles,&
             receive_n_particles_for_output_from_ipe(ipe),type_particle,ipe,&
             ipe,icomm,status,ierrmpi)
          call MPI_RECV(receive_payload,&
             npayload*receive_n_particles_for_output_from_ipe(ipe),&
             MPI_DOUBLE_PRECISION,ipe,ipe,icomm,status,ierrmpi)
          do ipart=1,receive_n_particles_for_output_from_ipe(ipe)
             call append_to_snapshot(receive_particles(ipart),&
                receive_payload(1:npayload,ipart))
          end do ! ipart
       end do ! ipe

       close(unitparticles)
    end if

  end subroutine output_ensemble_dat
  !=============================================================================
  subroutine output_particle(myparticle,payload,ipe,file_unit)

    use mod_metric, only: CoordToCart
    use mod_amrvacdef
    type(particle_node),intent(in)                 :: myparticle
    double precision, intent(in)                   :: payload(npayload)
    integer, intent(in)                            :: ipe
    integer, intent(in)                            :: file_unit
    double precision                               :: x(3)
    double precision                               :: xCart(3)
    double precision                               :: u(3)
    character(len=1024)                            :: line, data
    integer                                        :: icomp
    double precision, parameter                    :: minvalue = 1.0d-99
    double precision                               :: roundoff
    !-----------------------------------------------------------------------------

    if (nocartesian) then
       x = myparticle%x
    else
       ! Compute Cartesian position:
       xCart = 0.0d0 !makes sure the coordinates larger than ndim are initialized (and zero)
       call CoordToCart(0,0,0,0,0,0,0,0,0,0,0,0,myparticle%x(1:ndim),&
          xCart(1:ndim))

       ! normalize the coordinate
       x = xCart*UNIT_LENGTH
    end if

    ! Compute Cartesian velocity components:
    !call u3CoordToCart(0^D&,0^D&,0^D&,0^D&,myparticle%x,myparticle%u,u)
    ! Keep the transport velocity in its coordinate form in u

    line = ''
    write(data,"(es14.6, 3a)")roundoff(myparticle%t,minvalue), ',  '
    line = trim(line)//trim(data)
    write(data,"(es14.6, 3a)")roundoff(myparticle%dt,minvalue), ',  '
    line = trim(line)//trim(data)
    do icomp = 1, 3
       write(data,"(es14.6, 3a)")roundoff(x(icomp),minvalue), ',  '
       line = trim(line)//trim(data)
    end do
    do icomp = 1, 3
       write(data,"(es14.6, 3a)")roundoff(myparticle%u(icomp),minvalue), ',  '
       line = trim(line)//trim(data)
    end do
    do icomp = 1, npayload
       write(data,"(es14.6, 3a)")roundoff(payload(icomp),minvalue), ',  '
       line = trim(line)//trim(data)
    end do
    write(data,"(i9.7, 3a)")ipe, ',  '
    line = trim(line)//trim(data)
    write(data,"(i12.10, 3a)")it_particles, ',  '
    line = trim(line)//trim(data)
    write(data,"(i9.7)")myparticle%index
    line = trim(line)//trim(data)

    write(file_unit,"(a)") trim(line)

  end subroutine output_particle
  !=============================================================================
  subroutine read_from_snapshot()

    integer                         :: index,igrid_particle,ipe_particle
    !-----------------------------------------------------------------------------

    read(unitparticles) index
    if (.not. allocated(particle(index)%self)) allocate(particle(index)%self)
    if (.not. allocated(particle(index)%payload)) allocate(particle&
       (index)%payload(1:npayload))
    particle(index)%self%index = index

    read(unitparticles) particle(index)%self%follow
    read(unitparticles) particle(index)%self%q
    read(unitparticles) particle(index)%self%m
    read(unitparticles) particle(index)%self%t
    read(unitparticles) particle(index)%self%dt
    read(unitparticles) particle(index)%self%x
    read(unitparticles) particle(index)%self%u
    read(unitparticles) particle(index)%payload(1:npayload)

    call find_particle_ipe(particle(index)%self%x,igrid_particle,ipe_particle)
    particle(index)%igrid = igrid_particle
    particle(index)%ipe   = ipe_particle

    call push_particle_into_particles_on_mype(index)

  end subroutine read_from_snapshot
  !=============================================================================
  subroutine comm_particles()

    use mod_amrvacdef

    integer                         :: ipart, iipart, igrid_particle,&
        ipe_particle, ipe, iipe
    integer                         :: index
    integer                         :: tag_send, tag_receive, send_buff,&
        rcv_buff
    integer                         :: status(MPI_STATUS_SIZE)
    integer, dimension(0:npe-1)     :: send_n_particles_to_ipe
    integer, dimension(0:npe-1)    :: receive_n_particles_from_ipe
    type(particle_node), allocatable, dimension(:,:)  :: send_particles
    type(particle_node), allocatable, dimension(:,:)  :: receive_particles
    double precision, allocatable, dimension(:,:,:)  :: send_payload
    double precision, allocatable, dimension(:,:,:)  :: receive_payload
    integer, allocatable, dimension(:,:)      :: &
       particle_index_to_be_sent_to_ipe
    integer, dimension(nparticles_per_cpu_hi) :: &
       particle_index_to_be_destroyed
    integer                                   :: destroy_n_particles_mype
    logical                                   :: BC_applied
    integer, allocatable, dimension(:)        :: sndrqst, rcvrqst
    integer, allocatable, dimension(:)        :: sndrqst_payload,&
        rcvrqst_payload
    integer                                   :: isnd, ircv
    integer, parameter                        :: maxneighbors=32 !maximum case: coarse-fine in 3D
    !-----------------------------------------------------------------------------

    send_n_particles_to_ipe(:)      = 0
    receive_n_particles_from_ipe(:) = 0
    destroy_n_particles_mype        = 0

    allocate( particle_index_to_be_sent_to_ipe(nparticles_per_cpu_hi,&
       0:npe-1) )
    allocate( sndrqst(1:npe-1), rcvrqst(1:npe-1) )
    allocate( sndrqst_payload(1:npe-1), rcvrqst_payload(1:npe-1) )
    sndrqst = MPI_REQUEST_NULL; rcvrqst = MPI_REQUEST_NULL;
    sndrqst_payload = MPI_REQUEST_NULL; rcvrqst_payload = MPI_REQUEST_NULL;

    ! check if and where to send each particle, destroy if necessary
    !    !$OMP PARALLEL DO PRIVATE(ipart)
    do iipart=1,nparticles_on_mype;ipart=particles_on_mype(iipart);

       ! first check if the particle should be destroyed (user defined criterion)
       if ( query_destroy(particle(ipart)%self) ) then
          !          !$OMP CRITICAL(destroy)
          destroy_n_particles_mype  = destroy_n_particles_mype + 1
          particle_index_to_be_destroyed(destroy_n_particles_mype) = ipart
          !          !$OMP END CRITICAL(destroy)
          cycle   
       end if

       ! is my particle still in the same igrid?
       if (.not.particle_in_igrid(ipart,particle(ipart)%igrid)) then
          call find_particle_ipe(particle(ipart)%self%x,igrid_particle,&
             ipe_particle)

          ! destroy particle if out of domain (signalled by return value of -1)
          if (igrid_particle == -1 )then 
             call apply_periodB(particle(ipart)%self,igrid_particle,&
                ipe_particle,BC_applied)
             if (.not. BC_applied .or. igrid_particle == -1) then
                !                !$OMP CRITICAL(destroy2)
                destroy_n_particles_mype  = destroy_n_particles_mype + 1
                particle_index_to_be_destroyed(destroy_n_particles_mype) &
                   = ipart
                !                !$OMP END CRITICAL(destroy2)
                cycle
             end if
          end if

          ! particle still present
          particle(ipart)%igrid = igrid_particle
          particle(ipart)%ipe = ipe_particle

          ! if we have more than one core, is it on another cpu?
          if (npe .gt. 1 .and. particle(ipart)%ipe .ne. mype) then
             !             !$OMP CRITICAL(send)
             send_n_particles_to_ipe(ipe_particle) = send_n_particles_to_ipe&
                (ipe_particle) + 1
             if (send_n_particles_to_ipe(ipe_particle) .gt. &
                nparticles_per_cpu_hi) call mpistop('comm_particles: wanting& &
                to send too many particles, increase nparticles_per_cpu_hi')
             particle_index_to_be_sent_to_ipe(send_n_particles_to_ipe&
                (ipe_particle),ipe_particle) = ipart
             !             !$OMP END CRITICAL(send)
          end if ! ipe_particle

       end if ! particle_in_grid

    end do ! ipart
    !    !$OMP END PARALLEL DO

    call destroy_particles(destroy_n_particles_mype,&
       particle_index_to_be_destroyed(1:destroy_n_particles_mype))

    ! get out when only one core:
    if (npe == 1) return

    ! communicate amount of particles to be sent/received
    isnd = 0; ircv = 0;
    do iipe=1,npe_neighbors;ipe=ipe_neighbor(iipe);
       tag_send    = mype * npe + ipe
       tag_receive = ipe * npe + mype
       isnd = isnd + 1; ircv = ircv + 1
       call MPI_ISEND(send_n_particles_to_ipe(ipe),1,MPI_INTEGER,ipe,tag_send,&
          icomm,sndrqst(isnd),ierrmpi)
       call MPI_IRECV(receive_n_particles_from_ipe(ipe),1,MPI_INTEGER,ipe,&
          tag_receive,icomm,rcvrqst(ircv),ierrmpi)
    end do

    call MPI_WAITALL(isnd,sndrqst,MPI_STATUSES_IGNORE,ierrmpi)
    call MPI_WAITALL(ircv,rcvrqst,MPI_STATUSES_IGNORE,ierrmpi)
    sndrqst = MPI_REQUEST_NULL; rcvrqst = MPI_REQUEST_NULL;

    ! allocate the send and receive buffers,
    ! If this is still not good enough, consider optimizing using lists of lists:
    allocate(send_particles(maxval(send_n_particles_to_ipe),&
       1:min(maxneighbors,npe-1)))
    allocate(receive_particles(maxval(receive_n_particles_from_ipe),&
       1:min(maxneighbors,npe-1)))
    allocate(send_payload(npayload,maxval(send_n_particles_to_ipe),&
       1:min(maxneighbors,npe-1)))
    allocate(receive_payload(npayload,maxval(receive_n_particles_from_ipe),&
       1:min(maxneighbors,npe-1)))

    isnd = 0; ircv = 0;
    ! send and receive the data of the particles
    do iipe=1,npe_neighbors;ipe=ipe_neighbor(iipe);
       tag_send    = mype * npe + ipe
       tag_receive = ipe * npe + mype

       ! should i send some particles to ipe?
       if (send_n_particles_to_ipe(ipe) .gt. 0) then

          ! create the send buffer
          do ipart = 1, send_n_particles_to_ipe(ipe)
             send_particles(ipart,iipe) = particle&
                (particle_index_to_be_sent_to_ipe(ipart,ipe))%self
             send_payload(1:npayload,ipart,iipe) = particle&
                (particle_index_to_be_sent_to_ipe(ipart,ipe))%payload&
                (1:npayload)
          end do ! ipart
          isnd = isnd + 1
          call MPI_ISEND(send_particles(:,iipe),send_n_particles_to_ipe(ipe),&
             type_particle,ipe,tag_send,icomm,sndrqst(isnd),ierrmpi)
          call MPI_ISEND(send_payload(:,:,iipe),&
             npayload*send_n_particles_to_ipe(ipe),MPI_DOUBLE_PRECISION,ipe,&
             tag_send,icomm,sndrqst_payload(isnd),ierrmpi)

       end if ! send .gt. 0

       ! should i receive some particles from ipe?
       if (receive_n_particles_from_ipe(ipe) .gt. 0) then

          ircv = ircv + 1
          call MPI_IRECV(receive_particles(:,iipe),&
             receive_n_particles_from_ipe(ipe),type_particle,ipe,tag_receive,&
             icomm,rcvrqst(ircv),ierrmpi)
          call MPI_IRECV(receive_payload(:,:,iipe),&
             npayload*receive_n_particles_from_ipe(ipe),MPI_DOUBLE_PRECISION,&
             ipe,tag_receive,icomm,rcvrqst_payload(ircv),ierrmpi)

       end if ! receive .gt. 0
    end do ! ipe

    call MPI_WAITALL(isnd,sndrqst,MPI_STATUSES_IGNORE,ierrmpi)
    call MPI_WAITALL(ircv,rcvrqst,MPI_STATUSES_IGNORE,ierrmpi)
    call MPI_WAITALL(isnd,sndrqst_payload,MPI_STATUSES_IGNORE,ierrmpi)
    call MPI_WAITALL(ircv,rcvrqst_payload,MPI_STATUSES_IGNORE,ierrmpi)

    do iipe=1,npe_neighbors;ipe=ipe_neighbor(iipe);
       do ipart = 1, send_n_particles_to_ipe(ipe)
          deallocate(particle(particle_index_to_be_sent_to_ipe(ipart,&
             ipe))%self)
          deallocate(particle(particle_index_to_be_sent_to_ipe(ipart,&
             ipe))%payload)
          call pull_particle_from_particles_on_mype&
             (particle_index_to_be_sent_to_ipe(ipart,ipe))
       end do ! ipart
    end do

    do iipe=1,npe_neighbors;ipe=ipe_neighbor(iipe);
       if (receive_n_particles_from_ipe(ipe) .gt. 0) then
          do ipart = 1, receive_n_particles_from_ipe(ipe)

             index = receive_particles(ipart,iipe)%index
             call push_particle_into_particles_on_mype(index)
             if (.not. allocated(particle(index)%self)) allocate(particle&
                (index)%self)
             particle(index)%self = receive_particles(ipart,iipe)
             if (.not. allocated(particle(index)%payload)) &
                allocate(particle(index)%payload(npayload))
             particle(index)%payload(1:npayload) = receive_payload(1:npayload,&
                ipart,iipe)

             ! since we don't send the igrid, need to re-locate it
             call find_particle_ipe(particle(index)%self%x,igrid_particle,&
                ipe_particle)
             particle(index)%igrid = igrid_particle
             particle(index)%ipe = ipe_particle

          end do ! ipart

       end if ! receive .gt. 0
    end do ! ipe

  end subroutine comm_particles
  !=============================================================================
  subroutine comm_particles_global()

    ! does not destroy particles like comm_particles

    use mod_amrvacdef

    integer                         :: ipart, iipart, igrid_particle,&
        ipe_particle, ipe, iipe
    integer                         :: index
    integer                         :: tag_send, tag_receive, send_buff,&
        rcv_buff
    integer                         :: status(MPI_STATUS_SIZE)
    integer, dimension(0:npe-1)     :: send_n_particles_to_ipe
    integer, dimension(0:npe-1)     :: receive_n_particles_from_ipe
    type(particle_node), dimension(nparticles_per_cpu_hi)  :: send_particles
    type(particle_node), dimension(nparticles_per_cpu_hi)  :: &
       receive_particles
    double precision, dimension(npayload,nparticles_per_cpu_hi)  :: &
       send_payload
    double precision, dimension(npayload,nparticles_per_cpu_hi)  :: &
       receive_payload
    integer, allocatable, dimension(:,:)      :: &
       particle_index_to_be_sent_to_ipe
    logical                                   :: BC_applied
    integer, allocatable, dimension(:)        :: sndrqst, rcvrqst
    integer                                   :: isnd, ircv
    !-----------------------------------------------------------------------------
    send_n_particles_to_ipe(:)      = 0
    receive_n_particles_from_ipe(:) = 0

    allocate( particle_index_to_be_sent_to_ipe(nparticles_per_cpu_hi,&
       0:npe-1) )
    allocate( sndrqst(1:npe-1), rcvrqst(1:npe-1) )
    sndrqst = MPI_REQUEST_NULL; rcvrqst = MPI_REQUEST_NULL;

    ! check if and where to send each particle, relocate all of them (in case grid has changed)
    !    !$OMP PARALLEL DO PRIVATE(ipart)
    do iipart=1,nparticles_on_mype;ipart=particles_on_mype(iipart);

       call find_particle_ipe(particle(ipart)%self%x,igrid_particle,&
          ipe_particle)

       particle(ipart)%igrid = igrid_particle
       particle(ipart)%ipe = ipe_particle

       ! if we have more than one core, is it on another cpu?
       if (particle(ipart)%ipe .ne. mype) then
          !          !$OMP CRITICAL(send)
          send_n_particles_to_ipe(ipe_particle) = send_n_particles_to_ipe&
             (ipe_particle) + 1
          if (send_n_particles_to_ipe(ipe_particle) .gt. &
             nparticles_per_cpu_hi) call mpistop('comm_particles_global: &
             wanting& to send too many particles, increase &
             nparticles_per_cpu_hi')
          particle_index_to_be_sent_to_ipe(send_n_particles_to_ipe&
             (ipe_particle),ipe_particle) = ipart
          !          !$OMP END CRITICAL(send)
       end if ! ipe_particle

    end do ! ipart
    !    !$OMP END PARALLEL DO

    ! get out when only one core:
    if (npe == 1) return


    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! communicate amount of particles to be sent/received
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! opedit: I think this was the main bottleneck.
    ! 31.12.2017: made it nonblocking (easy, least invasive)
    ! If this continues to be problematic, better to 
    ! send particles with block in coarsen/refine/loadbalance. (best). 

    isnd = 0; ircv = 0;
    do ipe=0,npe-1; if (ipe .eq. mype) cycle;

       tag_send = mype * npe + ipe;  tag_receive = ipe * npe + mype;
       isnd = isnd + 1;  ircv = ircv + 1;
       call MPI_ISEND(send_n_particles_to_ipe(ipe),1,MPI_INTEGER, ipe,&
          tag_send,icomm,sndrqst(isnd),ierrmpi)
       call MPI_IRECV(receive_n_particles_from_ipe(ipe),1,MPI_INTEGER, ipe,&
          tag_receive,icomm,rcvrqst(ircv),ierrmpi)
    end do

    call MPI_WAITALL(isnd,sndrqst,MPI_STATUSES_IGNORE,ierrmpi)
    call MPI_WAITALL(ircv,rcvrqst,MPI_STATUSES_IGNORE,ierrmpi)
    deallocate( sndrqst, rcvrqst )


    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! send and receive the data of the particles
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do ipe=0,npe-1; if (ipe .eq. mype) cycle;
       tag_send    = mype * npe + ipe
       tag_receive = ipe * npe + mype

       ! should i send some particles to ipe?
       if (send_n_particles_to_ipe(ipe) .gt. 0) then

          ! create the send buffer
          do ipart = 1, send_n_particles_to_ipe(ipe)
             send_particles(ipart) = particle&
                (particle_index_to_be_sent_to_ipe(ipart,ipe))%self
             send_payload(1:npayload,ipart) = particle&
                (particle_index_to_be_sent_to_ipe(ipart,ipe))%payload&
                (1:npayload)
          end do ! ipart

          call MPI_SEND(send_particles,send_n_particles_to_ipe(ipe),&
             type_particle,ipe,tag_send,icomm,ierrmpi)
          call MPI_SEND(send_payload,npayload*send_n_particles_to_ipe(ipe),&
             MPI_DOUBLE_PRECISION,ipe,tag_send,icomm,ierrmpi)
          do ipart = 1, send_n_particles_to_ipe(ipe)
             deallocate(particle(particle_index_to_be_sent_to_ipe(ipart,&
                ipe))%self)
             deallocate(particle(particle_index_to_be_sent_to_ipe(ipart,&
                ipe))%payload)
             call pull_particle_from_particles_on_mype&
                (particle_index_to_be_sent_to_ipe(ipart,ipe))
          end do ! ipart

       end if ! send .gt. 0

       ! should i receive some particles from ipe?
       if (receive_n_particles_from_ipe(ipe) .gt. 0) then

          call MPI_RECV(receive_particles,receive_n_particles_from_ipe(ipe),&
             type_particle,ipe,tag_receive,icomm,status,ierrmpi)
          call MPI_RECV(receive_payload,npayload*receive_n_particles_from_ipe&
             (ipe),MPI_DOUBLE_PRECISION,ipe,tag_receive,icomm,status,ierrmpi)
          do ipart = 1, receive_n_particles_from_ipe(ipe)

             index = receive_particles(ipart)%index
             if (.not. allocated(particle(index)%self)) allocate(particle&
                (index)%self)
             particle(index)%self = receive_particles(ipart)
             if (.not. allocated(particle(index)%payload)) &
                allocate(particle(index)%payload(npayload))
             particle(index)%payload(1:npayload) = receive_payload(1:npayload,&
                ipart)
             call push_particle_into_particles_on_mype(index)

             ! since we don't send the igrid, need to re-locate it
             call find_particle_ipe(particle(index)%self%x,igrid_particle,&
                ipe_particle)
             particle(index)%igrid = igrid_particle
             particle(index)%ipe = ipe_particle

          end do ! ipart

       end if ! receive .gt. 0
    end do ! ipe

  end subroutine comm_particles_global
  !=============================================================================
  subroutine apply_periodB(particle,igrid_particle,ipe_particle,BC_applied)

    use mod_amrvacdef
    type(particle_node), intent(inout)        :: particle
    integer, intent(inout)                    :: igrid_particle, ipe_particle
    logical,intent(out)                       :: BC_applied
    integer                                   :: idim, iside
    !-----------------------------------------------------------------------------
    BC_applied = .false.

    ! first check out poles:
    if (any(poleB)) then
       do idim = 1,ndim
          do iside = 1,2
             if (.not. poleB(iside,idim)) cycle
             select case(idim)
                case(1)
                if (iside == 1) then
                   if (particle%x(1) < xprobmin1) then
                      particle%x(1)   = 2.0d0*xprobmin1 - particle%x(1)
                      particle%x(phi_) = modulo(particle%x(phi_) + dpi,&
                          2.0d0*dpi)
                      particle%u(1)  = - particle%u(1)
                      BC_applied = .true.
                   end if
                else if (iside == 2) then
                   if (particle%x(1) > xprobmax1) then
                      particle%x(1)   = 2.0d0*xprobmax1 - particle%x(1)
                      particle%x(phi_) = modulo(particle%x(phi_) + dpi,&
                          2.0d0*dpi)
                      particle%u(1)  = - particle%u(1)
                      BC_applied = .true.
                   end if
                end if
                
                case(2)
                if (iside == 1) then
                   if (particle%x(2) < xprobmin2) then
                      particle%x(2)   = 2.0d0*xprobmin2 - particle%x(2)
                      particle%x(phi_) = modulo(particle%x(phi_) + dpi,&
                          2.0d0*dpi)
                      particle%u(2)  = - particle%u(2)
                      BC_applied = .true.
                   end if
                else if (iside == 2) then
                   if (particle%x(2) > xprobmax2) then
                      particle%x(2)   = 2.0d0*xprobmax2 - particle%x(2)
                      particle%x(phi_) = modulo(particle%x(phi_) + dpi,&
                          2.0d0*dpi)
                      particle%u(2)  = - particle%u(2)
                      BC_applied = .true.
                   end if
                end if
                
                case(3)
                if (iside == 1) then
                   if (particle%x(3) < xprobmin3) then
                      particle%x(3)   = 2.0d0*xprobmin3 - particle%x(3)
                      particle%x(phi_) = modulo(particle%x(phi_) + dpi,&
                          2.0d0*dpi)
                      particle%u(3)  = - particle%u(3)
                      BC_applied = .true.
                   end if
                else if (iside == 2) then
                   if (particle%x(3) > xprobmax3) then
                      particle%x(3)   = 2.0d0*xprobmax3 - particle%x(3)
                      particle%x(phi_) = modulo(particle%x(phi_) + dpi,&
                          2.0d0*dpi)
                      particle%u(3)  = - particle%u(3)
                      BC_applied = .true.
                   end if
                end if
                
             end select
          end do
       end do
    end if

    ! get out if we don't have any periodic BC:
    if (.not. any(periodB(1:ndim))) return

    ! go through dimensions and try re-inject the particle at the other side
    do idim=1,ndim

       if (.not. periodB(idim)) cycle

       select case(idim)
          case (1)
          if (particle%x(1) .lt. xprobmin1) then
             particle%x(1) = particle%x(1) + (xprobmax1 - xprobmin1)
             BC_applied = .true.
          end if
          if (particle%x(1) .ge. xprobmax1) then
             particle%x(1) = particle%x(1) - (xprobmax1 - xprobmin1)
             BC_applied = .true.
          end if
          
          case (2)
          if (particle%x(2) .lt. xprobmin2) then
             particle%x(2) = particle%x(2) + (xprobmax2 - xprobmin2)
             BC_applied = .true.
          end if
          if (particle%x(2) .ge. xprobmax2) then
             particle%x(2) = particle%x(2) - (xprobmax2 - xprobmin2)
             BC_applied = .true.
          end if
          
          case (3)
          if (particle%x(3) .lt. xprobmin3) then
             particle%x(3) = particle%x(3) + (xprobmax3 - xprobmin3)
             BC_applied = .true.
          end if
          if (particle%x(3) .ge. xprobmax3) then
             particle%x(3) = particle%x(3) - (xprobmax3 - xprobmin3)
             BC_applied = .true.
          end if
          
       end select

    end do

    call find_particle_ipe(particle%x,igrid_particle,ipe_particle)

  end subroutine apply_periodB
  !=============================================================================
  subroutine destroy_particles(destroy_n_particles_mype,&
     particle_index_to_be_destroyed)
    ! clean up destroyed particles on all cores

    use mod_amrvacdef

    integer, intent(in)                                   :: &
       destroy_n_particles_mype
    integer, dimension(destroy_n_particles_mype), intent(in) :: &
       particle_index_to_be_destroyed
    type(particle_node), dimension(destroy_n_particles_mype):: &
       destroy_particles_mype
    double precision, dimension(npayload,destroy_n_particles_mype):: &
       destroy_payload_mype
    integer                                               :: iipart,ipart,&
       destroy_n_particles
    !-----------------------------------------------------------------------------
    destroy_n_particles             = 0

    ! append the particle to list of destroyed particles
    do iipart=1,destroy_n_particles_mype
    ipart=particle_index_to_be_destroyed(iipart);
       destroy_particles_mype(iipart) = particle(ipart)%self
       destroy_payload_mype(1:npayload,iipart) = particle(ipart)%payload&
          (1:npayload)
    end do

    call output_ensemble(destroy_n_particles_mype,destroy_particles_mype,&
       destroy_payload_mype,'destroy')

    if (npe > 1) then
       call MPI_ALLREDUCE(destroy_n_particles_mype,destroy_n_particles,1,&
          MPI_INTEGER,MPI_SUM,icomm,ierrmpi)
    else
       destroy_n_particles = destroy_n_particles_mype
    end if

    nparticles = nparticles - destroy_n_particles

    do iipart=1,destroy_n_particles_mype
    ipart=particle_index_to_be_destroyed(iipart);

       particle(ipart)%igrid = -1
       particle(ipart)%ipe   = -1

       write(*,*) particle(ipart)%self%t, ': particle',ipart,&
          ' has left at it ',it_particles,' on pe', mype
       write(*,*) particle(ipart)%self%t, ': particles remaining:',nparticles,&
           '(total)', nparticles_on_mype-1, 'on pe', mype

       deallocate(particle(ipart)%self)
       deallocate(particle(ipart)%payload)
       call pull_particle_from_particles_on_mype(ipart)
    end do

  end subroutine destroy_particles
  !=============================================================================
  subroutine push_particle_into_particles_on_mype(ipart)

    use mod_amrvacdef

    integer, intent(in)            :: ipart
    !-----------------------------------------------------------------------------

    nparticles_on_mype = nparticles_on_mype + 1
    if (nparticles_on_mype .gt. nparticles_per_cpu_hi) call &
       mpistop("push_particle_into_particles_on_mype: trying to handle more &
       than nparticles_per_cpu_hi particles !")
    particles_on_mype(nparticles_on_mype) = ipart

  end subroutine push_particle_into_particles_on_mype
  !=============================================================================
  subroutine pull_particle_from_particles_on_mype(ipart)

    implicit none

    integer, intent(in)            :: ipart
    integer                        :: i
    !-----------------------------------------------------------------------------


    do i=1,nparticles_on_mype
       if (particles_on_mype(i) == ipart) then
          particles_on_mype(i) = particles_on_mype(nparticles_on_mype)
          exit
       end if
    end do

    nparticles_on_mype = nparticles_on_mype - 1

  end subroutine pull_particle_from_particles_on_mype
  !=============================================================================
  subroutine read_particles_snapshot()

    use mod_amrvacdef
    logical,save                    :: file_exists=.false.
    character(len=128)              :: filename
    integer                         :: mynpayload, mynparticles, ipart, iipart
    integer, dimension(0:1)         :: buff
    !-----------------------------------------------------------------------------
    ! some initialisations:
    nparticles_on_mype = 0
    mynparticles       = 0
    nparticles         = 0
    it_particles       = 0

    ! open the snapshot file on the headnode
    if (mype .eq. 0) then
       write(filename,"(a,a,i4.4,a)") trim(filenameini),'_particles',&
          snapshotini,'.dat'
       INQUIRE(FILE=filename, EXIST=file_exists)

       if (.not. file_exists) then
          write(*,*) 'File '//trim(filename)//' with particle data does not exist, calling init_particles instead'
          buff(0) = -1
          buff(1) = -1

       else
          open(unit=unitparticles,file=filename,form='unformatted',action&
             ='read',status='unknown',access='stream')
          read(unitparticles) nparticles,it_particles,mynpayload
          if (mynpayload .ne. npayload) call mpistop&
('npayload in restart file does not match npayload in mod_particles')
          buff(0) = nparticles
          buff(1) = it_particles
       end if

    end if

    if (npe>0) call MPI_BCAST(buff,2,MPI_INTEGER,0,icomm,ierrmpi)

    ! check if the particle data was found:
    if (buff(1) .eq. -1) then 

       call initial_particles
       return

    end if

    ! particle data is there, fill variables:
    nparticles   = buff(0)
    it_particles = buff(1)
    index_free   = nparticles + 1

    do while (mynparticles .lt. nparticles)

       if (nparticles_on_mype .ge. nparticles_per_cpu_hi) then
          call mpistop("read_particles_snapshot: max number of particles per &
             cpu reached, increase nparticles_per_cpu_hi !")
       end if

       if (mype .eq. 0) then

          do while (nparticles_on_mype .lt. nparticles_per_cpu_hi &
             .and. mynparticles .lt. nparticles)
             call read_from_snapshot
             mynparticles = mynparticles + 1
          end do

       end if ! mype==0

       if (npe>0) call MPI_BCAST(mynparticles,1,MPI_INTEGER,0,icomm,ierrmpi)

       call comm_particles_global

    end do

    if (mype .eq. 0) close(unit=unitparticles)

  end subroutine read_particles_snapshot
  !=============================================================================
  subroutine write_particles_snapshot()

    use mod_amrvacdef

    character(len=128)                                  :: filename
    type(particle_node), dimension(nparticles_per_cpu_hi)  :: send_particles
    type(particle_node), dimension(nparticles_per_cpu_hi)  :: &
       receive_particles
    double precision, dimension(npayload,nparticles_per_cpu_hi)  :: &
       send_payload
    double precision, dimension(npayload,nparticles_per_cpu_hi)  :: &
       receive_payload
    integer                         :: status(MPI_STATUS_SIZE)
    integer,dimension(0:npe-1)      :: receive_n_particles_for_output_from_ipe
    integer                         :: ipe, ipart, iipart, &
       send_n_particles_for_output
    logical,save                    :: file_exists=.false.
    !-----------------------------------------------------------------------------

    receive_n_particles_for_output_from_ipe(:) = 0

    ! open the snapshot file on the headnode
    if (mype .eq. 0) then
       write(filename,"(a,a,i4.4,a)") trim(filenameout),'_particles',&
          snapshot-1,'.dat'
       INQUIRE(FILE=filename, EXIST=file_exists)
       if (.not. file_exists) then
          open(unit=unitparticles,file=filename,form='unformatted',status&
             ='new',access='stream')
       else
          open(unit=unitparticles,file=filename,form='unformatted',status&
             ='replace',access='stream')
       end if
       write(unitparticles) nparticles,it_particles,npayload
    end if

    if (npe==1) then
       do iipart=1,nparticles_on_mype;ipart=particles_on_mype(iipart);
          call append_to_snapshot(particle(ipart)%self,particle&
             (ipart)%payload)
       end do
       return
    end if

    if (mype .ne. 0) then
       call MPI_SEND(nparticles_on_mype,1,MPI_INTEGER,0,mype,icomm,ierrmpi)
       ! fill the send_buffer
       send_n_particles_for_output = nparticles_on_mype
       do iipart=1,nparticles_on_mype;ipart=particles_on_mype(iipart);
          send_particles(iipart) = particle(ipart)%self
          send_payload(1:npayload,iipart) = particle(ipart)%payload&
             (1:npayload)
       end do
    else
       ! get number of particles on other ipes
       do ipe=1,npe-1
          call MPI_RECV(receive_n_particles_for_output_from_ipe(ipe),1,&
             MPI_INTEGER,ipe,ipe,icomm,status,ierrmpi)
       end do
    end if

    if (mype .ne. 0) then
       call MPI_SEND(send_particles,send_n_particles_for_output,type_particle,&
          0,mype,icomm,ierrmpi)
       call MPI_SEND(send_payload,npayload*send_n_particles_for_output,&
          MPI_DOUBLE_PRECISION,0,mype,icomm,ierrmpi)
    end if

    if (mype==0) then
       ! first output own particles (already in receive buffer)
       do iipart=1,nparticles_on_mype;ipart=particles_on_mype(iipart);
          call append_to_snapshot(particle(ipart)%self,particle&
             (ipart)%payload)
       end do
       ! now output the particles sent from the other ipes
       do ipe=1,npe-1
          call MPI_RECV(receive_particles,&
             receive_n_particles_for_output_from_ipe(ipe),type_particle,ipe,&
             ipe,icomm,status,ierrmpi)
          call MPI_RECV(receive_payload,&
             npayload*receive_n_particles_for_output_from_ipe(ipe),&
             MPI_DOUBLE_PRECISION,ipe,ipe,icomm,status,ierrmpi)
          do ipart=1,receive_n_particles_for_output_from_ipe(ipe)
             call append_to_snapshot(receive_particles(ipart),&
                receive_payload(1:npayload,ipart))
          end do ! ipart
       end do ! ipe
       close(unit=unitparticles)
    end if ! mype == 0

  end subroutine write_particles_snapshot
  !=============================================================================
  subroutine append_to_snapshot(myparticle,mypayload)

    type(particle_node),intent(in) :: myparticle
    double precision, intent(in)   :: mypayload(1:npayload)

    write(unitparticles) myparticle%index
    write(unitparticles) myparticle%follow
    write(unitparticles) myparticle%q
    write(unitparticles) myparticle%m
    write(unitparticles) myparticle%t
    write(unitparticles) myparticle%dt
    write(unitparticles) myparticle%x
    write(unitparticles) myparticle%u
    write(unitparticles) mypayload(1:npayload)

  end subroutine append_to_snapshot
  !=============================================================================
  subroutine time_spent_on_particles()

    use mod_timing
    use mod_amrvacdef

    double precision         :: tpartc_avg, tpartc_int_avg, tpartc_io_avg,&
        tpartc_com_avg, tpartc_grid_avg

    !-----------------------------------------------------------------------------
    call MPI_REDUCE(tpartc,tpartc_avg,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,&
       ierrmpi)
    call MPI_REDUCE(tpartc_int,tpartc_int_avg,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
       0,icomm,ierrmpi)
    call MPI_REDUCE(tpartc_io,tpartc_io_avg,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
       icomm,ierrmpi)
    call MPI_REDUCE(tpartc_com,tpartc_com_avg,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
       0,icomm,ierrmpi)
    call MPI_REDUCE(tpartc_grid,tpartc_grid_avg,1,MPI_DOUBLE_PRECISION,&
       MPI_SUM,0,icomm,ierrmpi)

    if( mype ==0 ) then
       tpartc_avg     = tpartc_avg/dble(npe)
       tpartc_int_avg = tpartc_int_avg/dble(npe)
       tpartc_io_avg  = tpartc_io_avg/dble(npe)
       tpartc_com_avg  = tpartc_com_avg/dble(npe)
       tpartc_grid_avg  = tpartc_grid_avg/dble(npe)
       if (convert) timeloop = tpartc
       write(*,'(a,f12.3,a)')' Particle handling took     : ',tpartc,' sec'
       write(*,'(a,f12.2,a)')'                  Percentage: ',100.0d0*tpartc&
          /timeloop,' %'
       write(*,'(a,f12.3,a)')' Particle IO took           : ',tpartc_io_avg,&
          ' sec'
       write(*,'(a,f12.2,a)')'                  Percentage: ',&
          100.0d0*tpartc_io_avg/timeloop,' %'
       write(*,'(a,f12.3,a)')' Particle COM took          : ',tpartc_com_avg,&
          ' sec'
       write(*,'(a,f12.2,a)')'                  Percentage: ',&
          100.0d0*tpartc_com_avg/timeloop,' %'
       write(*,'(a,f12.3,a)')' Particle integration took  : ',tpartc_int_avg,&
          ' sec'
       write(*,'(a,f12.2,a)')'                  Percentage: ',&
          100.0d0*tpartc_int_avg/timeloop,' %'
       write(*,'(a,f12.3,a)')' Particle init grid took    : ',tpartc_grid_avg,&
          ' sec'
       write(*,'(a,f12.2,a)')'                  Percentage: ',&
          100.0d0*tpartc_grid_avg/timeloop,' %'
    end if


  end subroutine time_spent_on_particles
  !=============================================================================
  subroutine init_gridvars()

    use mod_amrvacdef
    integer                                   :: igrid, iigrid
    !-----------------------------------------------------------------------------

    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       if (.not. associated(gridvars(igrid)%w)) allocate(gridvars(igrid)%w&
          (ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:ngridvars))
       if (time_advance) then
          if (.not. associated(gridvars_old(igrid)%w)) allocate(gridvars_old&
             (igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:ngridvars))
       end if
    end do
    !$OMP END PARALLEL DO

    call fill_gridvars

  end subroutine init_gridvars
  !=============================================================================
  subroutine finish_gridvars()

    use mod_amrvacdef
    integer                                   :: igrid, iigrid
    !-----------------------------------------------------------------------------

    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       if (associated(gridvars(igrid)%w)) deallocate(gridvars(igrid)%w)
       if (time_advance) then
          if (associated(gridvars_old(igrid)%w)) deallocate(gridvars_old&
             (igrid)%w)
       end if
    end do
    !$OMP END PARALLEL DO

  end subroutine finish_gridvars
  !=============================================================================
  subroutine interpolate_var(igrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,gf,x,xloc,gfloc)

    use mod_amrvacdef
    integer, intent(in)                   :: igrid,ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3
    double precision, intent(in)          :: gf(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, intent(in)          :: x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision, intent(in)          :: xloc(1:ndir)
    double precision, intent(out)         :: gfloc
    integer                               :: ic1,ic2,ic3, ic11,ic12,ic13,&
        ic21,ic22,ic23, idir
    double precision                      :: xd1,xd2,xd3
    
    
    double precision                      :: c0, c1, c00, c10, c01, c11
    
    character(len=1024)                   :: line
    !-----------------------------------------------------------------------------

    ! flat interpolation:
    ic1 = int((xloc(1)-rnode(rpxmin1_,igrid))/rnode(rpdx1_,igrid)) + 1 + dixB 
    ic2 = int((xloc(2)-rnode(rpxmin2_,igrid))/rnode(rpdx2_,igrid)) + 1 + dixB 
    ic3 = int((xloc(3)-rnode(rpxmin3_,igrid))/rnode(rpdx3_,igrid)) + 1 + dixB 
    !gfloc = gf(ic^D)

     
         if (ic1.lt.ixImin1 .or. ic1.gt.ixImax1) then
    line = ''
    write(line,"(a)")'Trying to flat-interpolate from out of grid!'
    write(line,"(a,a,i3.2)")trim(line),' direction: ',1
    write(line,"(a,a,3es14.6)")trim(line),' position: ',xloc(1:ndim)
    write(line,"(a,a,i4.3)")trim(line),' index: ', ic1
    write(line,"(a,a,2es14.6)")trim(line),' grid range: ', x(ixImin1,ixImin2,&
       ixImin3,1),x(ixImax1,ixImax2,ixImax3,1)
    call mpistop(line)
 end if
 
      
         if (ic2.lt.ixImin2 .or. ic2.gt.ixImax2) then
    line = ''
    write(line,"(a)")'Trying to flat-interpolate from out of grid!'
    write(line,"(a,a,i3.2)")trim(line),' direction: ',2
    write(line,"(a,a,3es14.6)")trim(line),' position: ',xloc(1:ndim)
    write(line,"(a,a,i4.3)")trim(line),' index: ', ic2
    write(line,"(a,a,2es14.6)")trim(line),' grid range: ', x(ixImin1,ixImin2,&
       ixImin3,2),x(ixImax1,ixImax2,ixImax3,2)
    call mpistop(line)
 end if
 
      
         if (ic3.lt.ixImin3 .or. ic3.gt.ixImax3) then
    line = ''
    write(line,"(a)")'Trying to flat-interpolate from out of grid!'
    write(line,"(a,a,i3.2)")trim(line),' direction: ',3
    write(line,"(a,a,3es14.6)")trim(line),' position: ',xloc(1:ndim)
    write(line,"(a,a,i4.3)")trim(line),' index: ', ic3
    write(line,"(a,a,2es14.6)")trim(line),' grid range: ', x(ixImin1,ixImin2,&
       ixImin3,3),x(ixImax1,ixImax2,ixImax3,3)
    call mpistop(line)
 end if
 

 ! linear interpolation:
 
 if (x(ic1,ic2,ic3,1) .lt. xloc(1)) then
    ic11 = ic1
 else
    ic11 = ic1 -1
 end if
 ic21 = ic11 + 1
 
 
 if (x(ic1,ic2,ic3,2) .lt. xloc(2)) then
    ic12 = ic2
 else
    ic12 = ic2 -1
 end if
 ic22 = ic12 + 1
 
 
 if (x(ic1,ic2,ic3,3) .lt. xloc(3)) then
    ic13 = ic3
 else
    ic13 = ic3 -1
 end if
 ic23 = ic13 + 1
 

  
      if (ic11.lt.ixGlo1 .or. ic21.gt.ixGhi1) then
 line = ''
 write(line,"(a)")'Trying to interpolate from out of grid!'
 write(line,"(a,a,i3.2)")trim(line),' direction: ',1
 write(line,"(a,a,3es14.6)")trim(line),' position: ',xloc(1:ndim)
 write(line,"(a,a,2i4.3)")trim(line),' indices: ', ic11,ic21
 write(line,"(a,a,2es14.6)")trim(line),' grid range: ', x(ixImin1,ixImin2,&
    ixImin3,1),x(ixImax1,ixImax2,ixImax3,1)

 call mpistop(line)
end if

   
      if (ic12.lt.ixGlo2 .or. ic22.gt.ixGhi2) then
 line = ''
 write(line,"(a)")'Trying to interpolate from out of grid!'
 write(line,"(a,a,i3.2)")trim(line),' direction: ',2
 write(line,"(a,a,3es14.6)")trim(line),' position: ',xloc(1:ndim)
 write(line,"(a,a,2i4.3)")trim(line),' indices: ', ic12,ic22
 write(line,"(a,a,2es14.6)")trim(line),' grid range: ', x(ixImin1,ixImin2,&
    ixImin3,2),x(ixImax1,ixImax2,ixImax3,2)

 call mpistop(line)
end if

   
      if (ic13.lt.ixGlo3 .or. ic23.gt.ixGhi3) then
 line = ''
 write(line,"(a)")'Trying to interpolate from out of grid!'
 write(line,"(a,a,i3.2)")trim(line),' direction: ',3
 write(line,"(a,a,3es14.6)")trim(line),' position: ',xloc(1:ndim)
 write(line,"(a,a,2i4.3)")trim(line),' indices: ', ic13,ic23
 write(line,"(a,a,2es14.6)")trim(line),' grid range: ', x(ixImin1,ixImin2,&
    ixImin3,3),x(ixImax1,ixImax2,ixImax3,3)

 call mpistop(line)
end if






xd1 = (xloc(1)-x(ic11,ic12,ic13,1)) / (x(ic21,ic12,ic13,1) - x(ic11,ic12,ic13,&
   1))      
xd2 = (xloc(2)-x(ic11,ic12,ic13,2)) / (x(ic11,ic22,ic13,2) - x(ic11,ic12,ic13,&
   2))      
xd3 = (xloc(3)-x(ic11,ic12,ic13,3)) / (x(ic11,ic12,ic23,3) - x(ic11,ic12,ic13,&
   3))    

c00 = gf(ic11,ic12,ic13) * (1.0d0 - xd1) + gf(ic21,ic12,ic13) * xd1
c10 = gf(ic11,ic22,ic13) * (1.0d0 - xd1) + gf(ic21,ic22,ic13) * xd1
c01 = gf(ic11,ic12,ic23) * (1.0d0 - xd1) + gf(ic21,ic12,ic23) * xd1
c11 = gf(ic11,ic22,ic23) * (1.0d0 - xd1) + gf(ic21,ic22,ic23) * xd1

c0  = c00 * (1.0d0 - xd2) + c10 * xd2
c1  = c01 * (1.0d0 - xd2) + c11 * xd2

gfloc = c0 * (1.0d0 - xd3) + c1 * xd3


end subroutine interpolate_var
!=============================================================================
subroutine get_vec(igrid,x,tloc,var,ibeg,iend)

use mod_amrvacdef

integer,intent(in)                                   :: igrid, ibeg, iend
double precision,dimension(3), intent(in)          :: x
double precision, intent(in)                         :: tloc
double precision,dimension(iend-ibeg+1), intent(out) :: var
double precision,dimension(iend-ibeg+1)              :: e1, e2
integer                                              :: ivar, iloc
double precision                                     :: td
!-----------------------------------------------------------------------------

if (.not.time_advance) then

 do ivar=ibeg,iend
    iloc = ivar-ibeg+1

    call interpolate_var(igrid,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
       ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,gridvars(igrid)%w&
       (ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,ivar),px(igrid)%x&
       (ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:ndim),x,var(iloc))

 end do

else

 td = (tloc/(UNIT_LENGTH/UNIT_VELOCITY) - t) / dt

 do ivar=ibeg,iend
    iloc = ivar-ibeg+1

    call interpolate_var(igrid,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
       ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,gridvars_old(igrid)%w&
       (ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,ivar),px(igrid)%x&
       (ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:ndim),x,e1(iloc))
    call interpolate_var(igrid,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
       ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,gridvars(igrid)%w&
       (ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,ivar),px(igrid)%x&
       (ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:ndim),x,e2(iloc))

    var(iloc) = e1(iloc) * (1.0d0 - td) + e2(iloc) * td

 end do

end if !.not.time_advance

end subroutine get_vec
!=============================================================================
subroutine handle_particles()

use mod_timing

use mod_amrvacdef
double precision       :: t_next_event
!-----------------------------------------------------------------------------

tpartc0 = MPI_WTIME()

call set_neighbor_ipe

tpartc_grid_0=MPI_WTIME()
call init_gridvars
tpartc_grid = tpartc_grid + (MPI_WTIME()-tpartc_grid_0)

tpartc_com0=MPI_WTIME()
call comm_particles_global
tpartc_com=tpartc_com + (MPI_WTIME()-tpartc_com0)


! set the maximum integration time   
if (time_advance) then
 tmax_particles = (t + dt)
else
 tmax_particles = t + (tmax-t)
end if

! main integration loop
do

 t_next_event = min(t_next_output,t_next_injection)

 if (tmax_particles >= t_next_event) then

    if (t_next_output <= t_next_injection) then

       ! integrate to t_next_output
       call advance_particles(t_next_output)

       tpartc_io_0=MPI_WTIME()
       if (write_ensemble) call particles_output
       timeio_tot=timeio_tot+(MPI_WTIME()-tpartc_io_0)
       tpartc_io=tpartc_io+(MPI_WTIME()-tpartc_io_0)

       t_next_output = max(t_next_output, t_particles) + dtsave_ensemble
    else
       ! integrate to t_next_injection
       call advance_particles(t_next_injection)

       ! add new particles (user defined)
       if (inject) call add_particles
       t_next_injection = max(t_next_injection, t_particles) + &
          dtinject_particles
    end if

 else
    ! integrate to tmax_particles
    call advance_particles(tmax_particles)
    exit
 end if

 if (nparticles == 0 .and. convert) then
    tpartc = tpartc + (MPI_WTIME() - tpartc0)
    call time_spent_on_particles
    exit
 end if

end do

! deallocate the interpolation arrays
!    call finish_gridvars

tpartc = tpartc + (MPI_WTIME() - tpartc0)

end subroutine handle_particles
!=============================================================================
subroutine advance_particles(end_time)

use mod_timing

include 'mpif.h'
double precision, intent(in) :: end_time !< Advance at most up to this time
! .. local ..
!-----------------------------------------------------------------------------

particle_evol: do

 call select_active_particles(end_time)

 if (exit_condition() .eqv. .true.) then
    exit particle_evol 
 end if

 tpartc_int_0=MPI_WTIME()
 call integrate_particles(end_time)
 tpartc_int=tpartc_int+(MPI_WTIME()-tpartc_int_0)

 tpartc_com0=MPI_WTIME()
 call comm_particles
 tpartc_com=tpartc_com + (MPI_WTIME()-tpartc_com0)

 it_particles = it_particles + 1

end do particle_evol

end subroutine advance_particles
!=============================================================================
subroutine finish_tracerparticles

call finish_gridvars
call finish_particles
call finish_particles_com
call finish_particles_vars

end subroutine finish_tracerparticles
!=============================================================================
end module mod_particles_base
!=============================================================================
!=============================================================================
