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

module mod_physicaldata
  use mod_indices, only: ngridshi
  implicit none
  save

  !=============================================================================

  type walloc
     double precision, dimension(:,:,:,:), pointer :: w=>Null()
     integer                                      ::  ixGmin1=-1,ixGmax1&
        =-1, ixGmin2=-1,ixGmax2=-1, ixGmin3=-1,ixGmax3=-1,nwmin=-1,nwmax=-1 !allocated index-range
     logical                                      :: allocated=.false.
  end type walloc

  ! General state datastructure:
  type state
     integer                     :: igrid=-1 !comes in handy, with igrid, we know everything about the AMR
     integer                     :: iwpos=0 !location of the w-array, in case its reconstructed, etc.
                                                   ! -1: unknown,
                                                   !  0: w is center (default),
                                                   !  1: w is interface 1
                                                   !  2: w is interface 2
                                                   !  3: w is interface 3
     logical                     :: w_is_primitive=.false. !are we in primitive state?
     logical                     :: is_coarse=.false. !are we a coarse buffer (psCoarse)?
     type(xalloc), pointer       :: x=>Null()      ! Center positions
     type(walloc), pointer       :: w=>Null() !Variables, normally center
     type(walloc), pointer       :: ws=>Null() !Staggered variables, always on interface
     type(walloc), pointer       :: we=>Null() !Edge variables, always on edge
     type(walloc), pointer       :: wc=>Null() !Corner variables, always on corner
     type(geoalloc), pointer     :: geo=>Null() !Points to the geometry for the block, geo%m is center-metric
  end type state

  
  type walloc_sub
     double precision, dimension(:,:,:), pointer :: w=>Null()
  end type walloc_sub
 
  

  ! The centered fluid variables live here:
  type(walloc), target, dimension(ngridshi) :: pw, pwold, pw1, pw2, pw3, pw4,&
      pwres
  type(walloc), target, dimension(ngridshi) :: pwCoarse, pwCoCo

  
  ! The staggered fluid variables live here:
  type(walloc), target, dimension(ngridshi) :: pws, pwsold, pws1, pws2, pws3,&
      pws4, pwsres
  type(walloc), target, dimension(ngridshi) :: pwsCoarse, pwsCoCo
  

  ! The meta-structure state:
  type(state), dimension(ngridshi) :: ps, psold, ps1, ps2, ps3, ps4, psres
  type(state), dimension(ngridshi) :: psCoarse, psCoCo

  type(walloc), dimension(ngridshi) :: pwio
  
  type(walloc), dimension(ngridshi) :: pwsio
  
  type(walloc), dimension(ngridshi), target :: pB0_cell,  pB0_face1,pB0_face2,&
     pB0_face3
  type(walloc), pointer :: myB0_cell=>Null(), myB0_face1=>Null(),myB0_face2&
     =>Null(),myB0_face3=>Null(), myB0=>Null()
  type(walloc_sub), dimension(ngridshi) :: pw_sub !For the center values on the slice
  type(walloc_sub), dimension(ngridshi) :: pwC_sub !For the corner values on the slice
  
  
  
  double precision, dimension(:,:,:), allocatable :: collapsedData
 

  type xalloc
     double precision, dimension(:,:,:,:), pointer :: x=>Null()
     integer                                      ::  ixGmin1=-1,ixGmax1&
        =-1, ixGmin2=-1,ixGmax2=-1, ixGmin3=-1,ixGmax3=-1,ndimmin=-1,ndimmax&
        =-1 !allocated index-range
     logical                                      :: allocated=.false.
  end type xalloc

  
  type xalloc_sub
     double precision, dimension(:,:,:), pointer :: x=>Null()
  end type xalloc_sub
 
  
  type(xalloc), dimension(ngridshi), target :: px, pxCoarse
  type(xalloc_sub), dimension(ngridshi) :: px_sub !For the centerpositions on the slice
  type(xalloc_sub), dimension(ngridshi) :: pxC_sub !For the cornerpositions on the slice

  type indexlist
     integer                                       :: i=-1
     integer                                       :: j=-1
     integer                                       :: k=-1
     double precision, dimension(:,:,:), pointer    :: elem => Null() !Data is stored here, everything else are pointers.
  end type indexlist

  type metric
     integer                                       :: nnonzero&
        =-1, nnonzeroBeta=-1, nnonzeroDgDk=-1, nnonzeroDalphaDj=-1
     integer                                       :: nnonzeroDbetaiDj=-1
     integer                                       ::  ixGmin1=-1,ixGmax1&
        =-1, ixGmin2=-1,ixGmax2=-1, ixGmin3=-1,ixGmax3=-1
     !Derivatives:
     type(indexlist), dimension(:), pointer        :: nonzeroDgDk => Null() !Non-zero elements: Metric derivatives
     type(indexlist), dimension(:), pointer        :: nonzeroDalphaDj &
        => Null() !Non-zero elements: Lapse derivatives
     type(indexlist), dimension(:), pointer        :: nonzeroDbetaiDj &
        => Null() !Non-zero elements: Shift derivatives
     type(indexlist), dimension(:,:,:), pointer    :: DgDk => Null() !Metric derivative
     type(indexlist), dimension(:,:), pointer      :: DbetaiDj => Null() !Shift derivative
     type(indexlist), dimension(:), pointer        :: DalphaDj => Null() !Lapse derivative
     !Primary elements:
     type(indexlist), dimension(:,:), pointer      :: g => Null() !Array of 4-metric elements
     type(indexlist), dimension(:), pointer        :: nonzero => Null() !Coordinatelist of non-zero elements
     type(indexlist), dimension(:), pointer        :: nonzeroBeta => Null() !Non-zero elements: Shift vector
     !Inverse of the three-metric:
     type(indexlist), dimension(:,:), pointer      :: gammainv => Null() !Array of inverse 3-metric elements
     double precision, dimension(:,:,:), pointer    :: alpha => Null() !lapse
     double precision, dimension(:,:,:), pointer    :: zero => Null() !auxillary zeroes to point to
     double precision, dimension(:,:,:), pointer    :: sqrtgamma => Null() !square root of the spatial determinant
     type(indexlist), dimension(:), pointer        :: beta => Null() !shift vector (contravariant)
     type(indexlist), dimension(:), pointer        :: betaD => Null() !shift vector (covariant, as in metric)
  end type metric

  type geoalloc
     double precision, dimension(:,:,:), pointer :: dvolume=>Null()
     double precision, dimension(:,:,:), pointer :: surfaceC1&
        =>Null(),surfaceC2=>Null(),surfaceC3=>Null(),surface1&
        =>Null(),surface2=>Null(),surface3=>Null()
     double precision, dimension(:,:,:,:), pointer :: dx=>Null(), xbar=>Null() !xbar is the barycenter-position
     type(metric), pointer                        :: m=>Null(), mSurface1&
        =>Null(),mSurface2=>Null(),mSurface3=>Null()
  end type geoalloc

  type(geoalloc), dimension(ngridshi), target :: pgeo, pgeoCoarse, pgeoCoCo
  type(geoalloc), pointer                     :: mygeo => Null()
  type(metric), pointer                       :: myM => Null() !points to metric on current face or center.

  ! ------------------------------ OMP directives ------------------------------
 !$OMP THREADPRIVATE(myB0_cell,myB0_face1,myB0_face2,myB0_face3,myB0,mygeo,myM)
  ! ----------------------------------------------------------------------------
  
  !=============================================================================
contains
  !=============================================================================
  subroutine alloc_pw(pw,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     nwmin,nwmax)

    type(walloc)              :: pw
    integer, intent(in)       :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, nwmin, nwmax
    ! .. local ..
    logical                   :: match
    !-----------------------------------------------------------------------------

    ! Did we previously use this structure?
    if (pw%allocated .eqv. .true.) then
       ! Check if the index-ranges match
       match = .true.
       if ( pw%ixGmin1 .ne. ixGmin1 .or. pw%ixGmin2 .ne. ixGmin2 &
          .or. pw%ixGmin3 .ne. ixGmin3 ) match = .false.
       if ( pw%ixGmax1 .ne. ixGmax1 .or. pw%ixGmax2 .ne. ixGmax2 &
          .or. pw%ixGmax3 .ne. ixGmax3 ) match = .false.
       if (pw%nwmin .ne. nwmin .or. pw%nwmax .ne. nwmax) match = .false.

       if (match) then
          ! previously allocated and indices match, return
          return
       else
          ! Indices changed, deallocate and newly allocate below
          call dealloc_pw(pw)
       end if
    end if

    pw%ixGmin1=ixGmin1;pw%ixGmin2=ixGmin2;pw%ixGmin3=ixGmin3;
    pw%ixGmax1=ixGmax1;pw%ixGmax2=ixGmax2;pw%ixGmax3=ixGmax3;
    pw%nwmin=nwmin; pw%nwmax=nwmax
    allocate(pw%w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       nwmin:nwmax))
    pw%allocated = .true.

  end subroutine alloc_pw
  !=============================================================================
  subroutine dealloc_pw(pw)

    type(walloc)              :: pw
    !-----------------------------------------------------------------------------

    deallocate(pw%w)
    pw%allocated = .false.
     pw%ixGmin1=-1;pw%ixGmax1=-1; pw%ixGmin2=-1;pw%ixGmax2=-1; pw%ixGmin3=-1
     pw%ixGmax3=-1
    pw%nwmin=-1
    pw%nwmax=-1 

  end subroutine dealloc_pw
  !=============================================================================
  subroutine alloc_px(px,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ndimmin,ndimmax)

    type(xalloc)              :: px
    integer, intent(in)       :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ndimmin, ndimmax
    ! .. local ..
    logical                   :: match
    !-----------------------------------------------------------------------------

    ! Did we previously use this structure?
    if (px%allocated .eqv. .true.) then
       ! Check if the index-ranges match
       match = .true.
       if ( px%ixGmin1 .ne. ixGmin1 .or. px%ixGmin2 .ne. ixGmin2 &
          .or. px%ixGmin3 .ne. ixGmin3 ) match = .false.
       if ( px%ixGmax1 .ne. ixGmax1 .or. px%ixGmax2 .ne. ixGmax2 &
          .or. px%ixGmax3 .ne. ixGmax3 ) match = .false.
       if (px%ndimmin .ne. ndimmin .or. px%ndimmax .ne. ndimmax) match &
          = .false.

       if (match) then
          ! previously allocated and indices match, return
          return
       else
          ! Indices changed, deallocate and newly allocate below
          call dealloc_px(px)
       end if
    end if

    px%ixGmin1=ixGmin1;px%ixGmin2=ixGmin2;px%ixGmin3=ixGmin3;
    px%ixGmax1=ixGmax1;px%ixGmax2=ixGmax2;px%ixGmax3=ixGmax3;
    px%ndimmin=ndimmin; px%ndimmax=ndimmax
    allocate(px%x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       ndimmin:ndimmax))
    px%allocated = .true.

  end subroutine alloc_px
  !=============================================================================
  subroutine dealloc_px(px)

    type(xalloc)              :: px
    !-----------------------------------------------------------------------------

    deallocate(px%x)
    px%allocated = .false.
    px%ixGmin1=-1;px%ixGmin2=-1;px%ixGmin3=-1;
    px%ixGmax1=-1;px%ixGmax2=-1;px%ixGmax3=-1;
    px%ndimmin=-1; px%ndimmax=-1

  end subroutine dealloc_px
  !=============================================================================
  subroutine alloc_state(ps,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3)
    ! Allocate the solution arrays within the state structure
    ! Using info from the physics module

    use mod_amrvacpar, only: nw, nws
    
    !-------------------------------------------------
    ! DEFINITIONS OF GLOBAL PARAMETERS AND VARIABLES
    !-------------------------------------------------
    INTEGER,PARAMETER:: r_=1, phi_=2, z_=3
    INTEGER,PARAMETER:: pphi_=2, zz_=3
!    include 'amrvacpar.f90'
    !-------------------------------------------------

    type(state)                                   :: ps
    integer, intent(in)                           :: ixGmin1,ixGmin2,ixGmin3,&
       ixGmax1,ixGmax2,ixGmax3 
    ! .. local ..
    
    integer                                       :: ixGsmin1,ixGsmin2,&
       ixGsmin3,ixGsmax1,ixGsmax2,ixGsmax3 
    
    !-----------------------------------------------------------------------------

    call alloc_pw(ps%w,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,1,nw)

    
     ixGsmin1 = ixGmin1-1; ixGsmax1 = ixGmax1; ixGsmin2 = ixGmin2-1
     ixGsmax2 = ixGmax2; ixGsmin3 = ixGmin3-1; ixGsmax3 = ixGmax3
    call alloc_pw(ps%ws,ixGsmin1,ixGsmin2,ixGsmin3,ixGsmax1,ixGsmax2,ixGsmax3,&
       1,nws)
    

  end subroutine alloc_state
  !=============================================================================
  subroutine dealloc_state(ps)
    ! Dellocate the solution arrays within the state structure
    type(state)                                   :: ps
    !-----------------------------------------------------------------------------

    call dealloc_pw(ps%w)
    
    call dealloc_pw(ps%ws)
    

  end subroutine dealloc_state
  !=============================================================================
  subroutine copy_state(a,b)
    ! Copys state a to state b as in b=a

    type(state), intent(in)                        :: a
    type(state), intent(inout)                     :: b
    !-----------------------------------------------------------------------------

    b%igrid = a%igrid
    b%iwpos = a%iwpos
    b%w_is_primitive = a%w_is_primitive

    call copy_pw(a%w,b%w)
    
    call copy_pw(a%ws,b%ws)
    

  end subroutine copy_state
  !=============================================================================
  subroutine copy_pw(a,b)
    ! Copys walloc structure a to new walloc structure b as in b=a
    type(walloc), intent(in)                       :: a
    type(walloc), intent(inout)                    :: b
    ! .. local ..
    !-----------------------------------------------------------------------------

    call alloc_pw(b,a%ixGmin1,a%ixGmin2,a%ixGmin3,a%ixGmax1,a%ixGmax2,&
       a%ixGmax3,a%nwmin,a%nwmax)

    b%w = a%w

  end subroutine copy_pw
  !=============================================================================  
end module mod_physicaldata
!=============================================================================
