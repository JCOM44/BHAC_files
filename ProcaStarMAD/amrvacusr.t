!=============================================================================
! amrvacusr.t.nul
!=============================================================================
!INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacnul/specialimpl.t
!INCLUDE:amrvacnul/usrflags.t
!INCLUDE:amrvacnul/correctaux_usr.t
!=============================================================================
subroutine initglobaldata_usr

  use mod_metric
  use mod_amrvacdef

  double precision, external :: lofr
  
  eqpar(gamma_)    = 4.0d0/3.0d0         ! Adiabatic index

  eqpar(a_)        = 0.0d0            ! Spin parameter
  eqpar(kappa_)    = 1.0d-3              ! entropy (will be re-normalized)
  
  eqpar(rin_)      = 20.0d0              ! inner radius
  eqpar(l_)        = lofr(41.0d0)        ! Angular momentum(rmax)

  eqpar(m_)        = 1.0d0               ! Black hole mass

  eqpar(beta_)     = 1.0d+1              ! Measure for Plasma-beta
  eqpar(rhomin_)   = 1.0d-5              ! Minimum density scale
  eqpar(pmin_)     = 1.0d-7              ! Minimum pressure scale

  {^IFCOORDMKS
  coordpar(h_)        = 0.0d0
  coordpar(R0_)       = 0.0d0
  }

  {^IFCOORDNUM
  coordpar(R0_)       = -1.0d0
  }


  {^IFCOORDRZ
  call dilaton3(0.1d0)
  if (mype .eq. 0) &
       print*, 'The horizon is at r=', coord_r0
  }

  {#IFDEF ISO
  eqpar(adiab_)    = eqpar(kappa_)  ! entropy
  }

end subroutine initglobaldata_usr
!=============================================================================
double precision function lofr(rmax)
  use mod_metric, only: get_g4, BLToCoord
  use mod_amrvacdef
  
  double precision, intent(in)                   :: rmax
  ! .. local ..
  double precision                               :: x(1:ndim), xBL(1:ndim)
  double precision, dimension(0:^NC,0:^NC)       :: g
  double precision, dimension(0:^NC,0:^NC,1:^NC) :: dgdk
  !-----------------------------------------------------------------------------

! Get the metric and its derivative at the desired density maximum: 
  g    = zero
  dgdk = zero  
  xBL(1) = rmax
  {^IFZIN
  xBL(z_) = dpi/2.0d0
  }
  {^IFPHIIN
  xBL(phi_) = 0.0d0
  }
  call BLToCoord(1^D&,1^D&,1^D&,1^D&,xBL,x)
  call get_g4({x(^D)},g,dgdk=dgdk) 
  
  lofr = (-1.0d0*dgdk(phi_,phi_,r_)*g(0,0)*g(0,phi_)+  &
       dgdk(0,phi_,r_)*g(0,phi_)**2+  &
       dgdk(0,phi_,r_)*g(0,0)*g(phi_,phi_)-  &
       1.0d0*dgdk(0,0,r_)*g(0,phi_)*g(phi_,phi_)+  &
       sqrt(dgdk(0,phi_,r_)**2*g(0,phi_)**4-  &
       1.0d0*dgdk(0,0,r_)*dgdk(phi_,phi_,r_)*g(0,phi_)**4-  &
       2.0d0*dgdk(0,phi_,r_)**2*g(0,0)*g(0,phi_)**2*g(phi_,phi_)  &
       +  &
       2.0d0*dgdk(0,0,r_)*dgdk(phi_,phi_,r_)*g(0,0)*g(0,phi_)**2*g(phi_,phi_) &
       +dgdk(0,phi_,r_)**2*g(0,0)**2*g(phi_,phi_)**2-  &
       1.0d0*dgdk(0,0,r_)*dgdk(phi_,phi_,r_)*g(0,0)**2*g(phi_,phi_)**2))/( &
       dgdk(phi_,phi_,r_)*g(0,0)**2  &
       -2.0d0*dgdk(0,phi_,r_)*g(0,0)*g(0,phi_)+  &
       dgdk(0,0,r_)*g(0,phi_)**2)

end function lofr
!=============================================================================
subroutine initonegrid_usr(ixI^L,ixO^L,s)

  use mod_metric,only: get_g_component, get_alpha, get_beta, lower3, get_g4, BLToCoord
  ! initialize one grid within ixO^L
  use mod_amrvacdef

  integer, intent(in)              :: ixI^L, ixO^L
  type(state), intent(inout)       :: s
  ! .. local variables ..
  logical                          :: patchw(ixG^T)
  double precision, dimension(ixG^T)  :: rcyl2, utd, omega, utu, h, d2
  double precision                    :: utd_in, rcyl2_in, gin(0:ndir,0:ndir)
  double precision                    :: xin(1:ndim), gmm1, betad_in(1:ndir),betau_in(1:ndir),betain2,d2in
  double precision                    :: xinBL(1:ndim)

  !-----------------------------------------------------------------------------
  associate(x=>s%x%x,w=>s%w%w{#IFDEF STAGGERED ,ws=>s%ws%w})

    gmm1 = (eqpar(gamma_) - 1.0d0) / eqpar(gamma_)

    w(ixI^S,:) = zero
    w(ixO^S,b3_) = zero

    patchw(ixO^S) = .false.

    rcyl2(ixO^S) = myM%g(0,phi_)%elem(ixO^S)**2-myM%g(0,0)%elem(ixO^S)*myM%g(phi_,phi_)%elem(ixO^S)

    d2(ixO^S) = (myM%g(0,0)%elem(ixO^S)*eqpar(l_)**2 &
         +2.0d0*myM%g(0,phi_)%elem(ixO^S)*eqpar(l_)+myM%g(phi_,phi_)%elem(ixO^S))

    where (d2(ixO^S) .gt. zero .and. rcyl2(ixO^S) .gt. zero)
       utd(ixO^S) = - sqrt(rcyl2(ixO^S)/d2(ixO^S))
    elsewhere
       utd(ixO^S) = -1.0d0
    end where


    omega(ixO^S) = - (myM%g(0,phi_)%elem(ixO^S) + myM%g(0,0)%elem(ixO^S) * eqpar(l_))
    omega(ixO^S) = omega(ixO^S) / (myM%g(phi_,phi_)%elem(ixO^S) +myM%g(0,phi_)%elem(ixO^S) * eqpar(l_))

    utu(ixO^S) = sqrt(-1.0d0/(myM%g(0,0)%elem(ixO^S) &
         + 2.0d0*omega(ixO^S)*myM%g(0,phi_)%elem(ixO^S) &
         + myM%g(phi_,phi_)%elem(ixO^S)*omega(ixO^S)**2 ) )


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Get metric at inner disk radius:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    gin(:,:) = zero
    xinBL(1) = eqpar(rin_)
    {^IFZIN
    xinBL(z_) = dpi/2.0d0
    }
    {^IFPHIIN
    xinBL(phi_) = 0.0d0
    }
    call BLToCoord(1^D&,1^D&,1^D&,1^D&,xinBL,xin)
    call get_g4({xin(^D)},gin) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    rcyl2_in = gin(0,phi_)**2-gin(0,0)*gin(phi_,phi_)
    d2in     = (gin(0,0)*eqpar(l_)**2+2.0d0*gin(0,phi_)*eqpar(l_)+gin(phi_,phi_))
    if (d2in >zero .and. rcyl2_in>zero) then 
       utd_in = - sqrt(rcyl2_in/d2in)
    else
       utd_in = -1.0d0
    end if

    h(ixO^S) = utd_in / utd(ixO^S)

    ! Now fill the primitive variables:

    where (utd_in/utd(ixO^S) .gt. 1.0d0 .and. d2(ixO^S) .gt. zero .and. x(ixO^S,1) .gt. xin(1))

       w(ixO^S,rho_) = (gmm1 * (h(ixO^S)-1.0d0)/eqpar(kappa_))**(1.0d0/(eqpar(gamma_)-1.0d0))
       w(ixO^S,pp_) = eqpar(kappa_) * w(ixO^S,rho_)**eqpar(gamma_)
       w(ixO^S,u0_+phi_) =  utu(ixO^S) * ( omega(ixO^S) + myM%beta(phi_)%elem(ixO^S) )
       w(ixO^S,u1_) =  utu(ixO^S) * myM%beta(1)%elem(ixO^S) 
       {^FL& w(ixO^S,tr^FL_) = one\}

    elsewhere

       w(ixO^S,rho_) = eqpar(rhomin_)*1.d-20
       w(ixO^S,pp_) = eqpar(pmin_)*1.d-20
       {^C&w(ixO^S,v^C_)   = myM%beta(^C)%elem(ixO^S)/myM%alpha(ixO^S)*zero \}
       {^FL& w(ixO^S,tr^FL_) = zero\}

    end where

    ! set atmosphere to floor:
    call set_atmo(ixI^L,ixO^L,w,x)

    call perturb(ixI^L,ixO^L, w, x)

    ! Magnetic field is set up in normalize initial (since it needs rhomax=1)

    {#IFDEF ENTROPY
    w(ixO^S,s_) = w(ixO^S,pp_) * w(ixO^S,rho_)**(-eqpar(gamma_))
    }

    call conserve(ixI^L,ixO^L,w,x,patchw)


  end associate
end subroutine initonegrid_usr
!=============================================================================
subroutine set_atmo(ixI^L,ixO^L,w,x)
  ! The atmosphere treatment

  use mod_metric, only: CoordToBL, lower3
  use mod_amrvacdef

  integer, intent(in)              :: ixI^L, ixO^L
  double precision, intent(in)     :: x(ixI^S,1:ndim)
  double precision, intent(inout)  :: w(ixI^S,1:nw)
  ! .. local ..
  double precision                 :: lfac2(ixI^S), tmp(ixI^S), vD(ixI^S,1:^NC)
  double precision                 :: xBL(ixI^S,1:ndim)
  logical, dimension(ixI^S)        :: patchw
  double precision, parameter      :: lfacmax=20.0d0
  double precision                 :: r0
  !-----------------------------------------------------------------------------

  r0=1.0d0 ! Small number to make the atmosphere bounded at the origin

  patchw(ixO^S) = .false.
  call CoordToBL(ixI^L,ixO^L,x,xBL)

  !! If the simulation includes the origin, it is important to remove the
  !! singularity in the atmosphere profle, with r0.

  !! Set rho to atmosphere where it is too small
  where(w(ixO^S,rho_) .lt. eqpar(rhomin_) * (r0+ xBL(ixO^S,1))**(-1.5d0))
     w(ixO^S,rho_) = eqpar(rhomin_) * (r0+ xBL(ixO^S,1))**(-1.5d0)
     patchw(ixO^S) = .true.
  end where

  !! Set p to atmosphere where it is too small
  where(w(ixO^S,pp_) .lt. eqpar(pmin_) * (r0+ xBL(ixO^S,1))**(-2.5d0))
     {#IFDEF ISO
     w(ixO^S,pp_)   = eqpar(adiab_) * w(ixO^S,rho_)**eqpar(gamma_)
     }{#IFNDEF ISO
     w(ixO^S,pp_)   = eqpar(pmin_) * (r0+ xBL(ixO^S,1))**(-2.5d0)
     }
     patchw(ixO^S) = .true.
  end where

  !! Set to atmosphere and set velocities to zero values where the tracer is small
  where(w(ixO^S,tr1_).lt.0.01d0)
     w(ixO^S,rho_) = eqpar(rhomin_) * (r0+ xBL(ixO^S,1))**(-1.5d0)
     {#IFDEF ISO
     w(ixO^S,pp_)   = eqpar(adiab_) * w(ixO^S,rho_)**eqpar(gamma_)
     }{#IFNDEF ISO
     w(ixO^S,pp_)   = eqpar(pmin_) * (r0+ xBL(ixO^S,1))**(-2.5d0)
     }
    {^C&w(ixO^S,v^C_) = zero\}
     patchw(ixO^S) = .true.
  end where

  call lower3(ixI^L,ixO^L,myM,w(ixI^S,u1_:u^NC_),vD)
  ! Limit the four-velocity
  lfac2(ixO^S) = {^C& w(ixO^S,v^C_)*vD(ixO^S,^C) |+}
  tmp(ixO^S)   = sqrt((lfacmax**2)/(lfac2(ixO^S)))
  where (lfac2(ixO^S) .gt. lfacmax**2)
     {^C&w(ixO^S,v^C_) = w(ixO^S,v^C_) * tmp(ixO^S)\}
  end where

  ! re-calculate lfac:
  call lower3(ixI^L,ixO^L,myM,w(ixI^S,u1_:u^NC_),vD)
  w(ixO^S,lfac_) = sqrt({^C& (w(ixO^S,u^C_)*vD(ixO^S,^C))|+} + 1.0d0)

  ! re-calculate xi:
  w(ixO^S,xi_) = w(ixO^S,lfac_)**2 * (w(ixO^S,rho_) + govergminone*w(ixO^S,pp_))

  {#IFDEF ENTROPY
  ! re-calculate entropy where density or pressure were floored:
  where ( patchw(ixO^S) .eqv. .true. )
     w(ixO^S,s_) = w(ixO^S,pp_) * w(ixO^S,rho_)**(-eqpar(gamma_))
  end where
  }

end subroutine set_atmo
!=============================================================================
subroutine perturb(ixI^L,ixO^L, w, x)
  use mod_amrvacdef

  integer, intent(in)              :: ixI^L, ixO^L
  double precision, intent(in)     :: x(ixI^S,1:ndim)
  double precision, intent(inout)  :: w(ixI^S,1:nw)
  ! .. local ..
  integer                          :: ix^D
  double precision                 :: r(ixO^S)
  double precision                 :: Xr
  double precision, parameter      :: amp=0.04d0
  !-----------------------------------------------------------------------------

  CALL init_random_seed()         ! see example of RANDOM_SEED
  CALL RANDOM_NUMBER(r)

  {do ix^D=ixOmin^D,ixOmax^D\}
  Xr = (r(ix^D)*2.0d0-one)
  w(ix^D,pp_) = w(ix^D,pp_) * (1.0d0 + amp*Xr)
  {enddo\}


end subroutine perturb
!=============================================================================
  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,s)

    ! special boundary types, user defined
    ! user must assign conservative variables in bounderies

    use mod_amrvacdef

    integer, intent(in) :: ixI^L, ixO^L, iB
    double precision, intent(in) :: qt
    type(state), intent(inout)   :: s
    !----------------------------------------------------------------------------
    associate(x=>s%x%x,w=>s%w%w{#IFDEF STAGGERED ,ws=>s%ws%w})
      call mpistop("specialbound not defined")


    end associate
  end subroutine specialbound_usr
  !=============================================================================
subroutine correctaux_usr(ixI^L,ixO^L,w,x,patchierror,subname)

use mod_metric, only: CoordToBL
use mod_amrvacdef

integer, intent(in)            :: ixI^L, ixO^L
integer, intent(inout)         :: patchierror(ixG^T)
character(len=*), intent(in)   :: subname
double precision, intent(inout):: w(ixI^S,1:nw)
double precision, intent(in)   :: x(ixI^S,1:ndim)
! .. local ..
double precision               :: wtmp(ixI^S,1:nw)
double precision               :: xBL(ixI^S,1:ndim)
logical                        :: patchw(ixG^T)
!-----------------------------------------------------------------------------

patchw(ixG^T) = .true.

call CoordToBL(ixI^L,ixO^L,x,xBL)

where (patchierror(ixO^S)/=0)
   w(ixO^S,rho_) = eqpar(rhomin_) * xBL(ixO^S,1)**(-1.5d0)
   {#IFDEF ISO
   w(ixO^S,pp_)   = eqpar(adiab_) * w(ixO^S,rho_)**eqpar(gamma_)
   }{#IFNDEF ISO
   w(ixO^S,pp_)   = eqpar(pmin_) * xBL(ixO^S,1)**(-2.5d0)
   }
   {^C&w(ixO^S,v^C_)  = zero \}
   {#IFDEF TRACER
   w(ixO^S,tr1_)      = 0.0d0
   }
   w(ixO^S,lfac_)     = 1.0d0
   w(ixO^S,xi_)       = w(ixO^S,rho_) + govergminone *  w(ixO^S,pp_)
   patchierror(ixO^S) = 0
   {#IFDEF ENTROPY
   ! re-calculate entropy:
   w(ixO^S,s_) = w(ixO^S,pp_) * w(ixO^S,rho_)**(-eqpar(gamma_))
   }
   patchw(ixO^S)      = .false.
end where

call conserven(ixI^L,ixO^L,w,patchw)

end subroutine correctaux_usr
  !=============================================================================
  subroutine fixp_usr(ixI^L,ixO^L,w,x)
    use mod_physaux, only:get_b2
    use mod_metric, only:lower3
!    use FM, only: lfacmax, set_atmo
    use mod_amrvacdef

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(inout)    :: w(ixI^S,1:nw)
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    ! .. local ..
    double precision                   :: b2(ixI^S), vD(ixI^S,1:^NC)
   !----------------------------------------------------------------------------
    
    call set_atmo(ixI^L,ixO^L,w,x)
    
  end subroutine fixp_usr
  !=============================================================================
  module normalize

    ! rhomax needs to be store for when we set-up the magnetic field in normalize_initial
    ! We initialize rhomax=-1 to prevent the magnetic field from being filled-in before
    ! normalize_initial.

    double precision, save  :: rhomax =-1.0d0
    
  end module normalize
  !=============================================================================
  subroutine normalize_initial()

    use mod_physaux, only:get_b2
    use mod_metric, only: CoordToBL
    use mod_multigrid_coupling
    use normalize, only: rhomax

    use mod_amrvacdef

    ! .. local ..
    integer                   :: iigrid, igrid
    double precision          :: b2max, pmax, norm
    double precision          :: b2(ixG^T)
    double precision          :: xBL(ixG^T,1:ndim)
    !-----------------------------------------------------------------------------

    pmax   = -bigdouble
    b2max  = -bigdouble
    rhomax = -bigdouble

    do iigrid=1,igridstail; igrid=igrids(iigrid);
       associate(x=>ps(igrid)%x%x,w=>ps(igrid)%w%w,s=>ps(igrid){#IFDEF STAGGERED ,ws=>ps(igrid)%ws%w})
         call set_tmpGlobals(igrid)

         call primitive(ixG^LL,ixM^LL,w,x)
         ! get the maximum of rho:      
         rhomax = max(maxval(pw(igrid)%w(ixM^T,rho_)),rhomax)
         ! get the maximum of p:      
         pmax = max(maxval(pw(igrid)%w(ixM^T,pp_)),pmax)
       end associate
    end do

    ! Communicate rhomax:
    call MPI_ALLREDUCE(MPI_IN_PLACE,rhomax,1,MPI_DOUBLE_PRECISION, &
         MPI_MAX,icomm,ierrmpi)
    ! Communicate pmax:
    call MPI_ALLREDUCE(MPI_IN_PLACE,pmax,1,MPI_DOUBLE_PRECISION, &
         MPI_MAX,icomm,ierrmpi)


    norm = 1.0d0/rhomax

    do iigrid=1,igridstail; igrid=igrids(iigrid);
       associate(x=>ps(igrid)%x%x,w=>ps(igrid)%w%w,s=>ps(igrid){#IFDEF STAGGERED ,ws=>ps(igrid)%ws%w})
         call set_tmpGlobals(igrid)

         call CoordToBL(ixG^LL,ixM^LL,x,xBL)
         ! re-normalize outside of floors
          where(w(ixM^T,rho_) .ge. eqpar(rhomin_) * xBL(ixM^T,1)**(-1.5d0))
            w(ixM^T,rho_) = w(ixM^T,rho_)*norm
            w(ixM^T,pp_)  = w(ixM^T,pp_)*norm
            {#IFDEF ENTROPY
            w(ixM^T,s_)   = w(ixM^T,pp_) * w(ixM^T,rho_)**(-eqpar(gamma_))
            }
          end where

         ! Now density is normalized to one.  Set up the magnetic field:
         {#IFNDEF STAGGERED
         call b_from_vectorpotential(ixG^LL,ixM^LL,w,x)
         }{#IFDEF STAGGERED
         call b_from_vectorpotential(s%ws%ixG^L,ixG^LL,ixM^LL,ws,x)
         call faces2centers(ixM^LL,s)
         }

         call conserve(ixG^LL,ixM^LL,w,x,patchfalse)

       end associate
    end do


  {#IFDEF STAGGERED
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! clean the magnetic field:
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  select case(clean_init_divB)
  case ('vecpot')
     ! Now that the grid has been set and the flux-conservation arrays are allocated,
     ! re-calculate the magnetic field from the vector potential in a completely
     ! divergence free way.
     if (levmax>levmin .and. ndim .eq. 3) call recalculateB
     {^NOONED
  case ('mg')
     ! Project out the divB using poisson solver.
     ! Due to Jannis Teunissen. Thanks!   
     call mg_setup_multigrid()
     call clean_divb_multigrid()
     }
  case ('none')
     ! do nothing
  case default
     call mpistop('Unknown method selected for clean_init_divB!')
     ! check what you are doing!
  end select
  }


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Normalize magnetic field:
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! since we changed pressure and density, max pressure is now

  ! print *, 'pmax before norm',pmax

  pmax = pmax * norm
  
  do iigrid=1,igridstail; igrid=igrids(iigrid);
     call set_tmpGlobals(igrid)
     ! get the maximum of b2:
     call get_b2(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,b2,w_is_primitive=.false.)
     b2max = max(maxval(b2(ixM^T)),b2max)
  end do

  ! Communicate b2max:
  call MPI_ALLREDUCE(MPI_IN_PLACE,b2max,1,MPI_DOUBLE_PRECISION, &
       MPI_MAX,icomm,ierrmpi)


  ! Apparently no magnetic field:
  if (b2max .le. 0.0d0) return

  ! Otherwise normalize field:
  norm = sqrt(2.0d0*pmax/b2max/eqpar(beta_))

  ! Now apply the normalization factor:
  do iigrid=1,igridstail; igrid=igrids(iigrid);
     associate(x=>ps(igrid)%x%x,w=>ps(igrid)%w%w{#IFDEF STAGGERED ,ws=>ps(igrid)%ws%w})
       call set_tmpGlobals(igrid)

       call primitive(ixG^LL,ixM^LL,w,x)

{#IFNDEF STAGGERED
       {^D& w(ixM^T,b^D_) = norm*w(ixM^T,b^D_) \}
}{#IFDEF STAGGERED
       {^D& ws(ps(igrid)%ws%ixG^S,bs^D_) = norm*ws(ps(igrid)%ws%ixG^S,bs^D_) \}
       call faces2centers(ixM^LL,ps(igrid))
}
       call conserve(ixG^LL,ixM^LL,w,x,patchfalse)

     end associate
  end do

  end subroutine normalize_initial  
!=============================================================================
subroutine flag_grid_usr(qt,ixG^L,ixO^L,w,x,flag)

use mod_amrvacdef

integer, intent(in)             :: ixG^L, ixO^L
integer, intent(inout)          :: flag
double precision, intent(in)    :: qt
double precision, intent(inout) :: w(ixG^S,1:nw)
double precision, intent(in)    :: x(ixG^S,1:ndim)

! flag=-1 : Treat all cells active, omit deactivation (onentry, default)
! flag=0  : Treat as normal domain
! flag=1  : Treat as passive, but reduce by safety belt
! flag=2  : Always treat as passive

!-----------------------------------------------------------------------------
      
end subroutine flag_grid_usr
!=============================================================================
subroutine initvecpot_usr(ixI^L, ixC^L, xC, A, idir)

! initialize the vectorpotential on the corners
! used by b_from_vectorpotential()

  use normalize, only: rhomax
  use mod_metric, only: get_g4, CoordToBL, BLToCoord
  use mod_amrvacdef

integer, intent(in)                :: ixI^L, ixC^L, idir
double precision, intent(in)       :: xC(ixI^S,1:ndim)
double precision, intent(out)      :: A(ixI^S)
! .. local ..
 double precision, dimension(ixG^T)  :: rcyl2, utd, omega, utu, h, d2, rho
 double precision                    :: utd_in, rcyl2_in, gin(0:ndir,0:ndir)
 double precision                    :: xin(1:ndim), gmm1, betad_in(1:ndir),betau_in(1:ndir),betain2,d2in
 double precision                    :: xinBL(1:ndim)
 double precision                    :: xBL(ixI^S,1:ndim)
!-----------------------------------------------------------------------------

 gmm1 = (eqpar(gamma_) - 1.0d0) / eqpar(gamma_)
 
 rcyl2(ixC^S) = myM%g(0,phi_)%elem(ixC^S)**2-myM%g(0,0)%elem(ixC^S)*myM%g(phi_,phi_)%elem(ixC^S)

 d2(ixC^S) = (myM%g(0,0)%elem(ixC^S)*eqpar(l_)**2 &
      +2.0d0*myM%g(0,phi_)%elem(ixC^S)*eqpar(l_)+myM%g(phi_,phi_)%elem(ixC^S))

 where (d2(ixC^S) .gt. zero .and. rcyl2(ixC^S) .gt. zero)
    utd(ixC^S) = - sqrt(rcyl2(ixC^S)/d2(ixC^S))
 elsewhere
    utd(ixC^S) = -1.0d0
 end where
 

 omega(ixC^S) = - (myM%g(0,phi_)%elem(ixC^S) + myM%g(0,0)%elem(ixC^S) * eqpar(l_))
 omega(ixC^S) = omega(ixC^S) / (myM%g(phi_,phi_)%elem(ixC^S) +myM%g(0,phi_)%elem(ixC^S) * eqpar(l_))

 utu(ixC^S) = sqrt(-1.0d0/(myM%g(0,0)%elem(ixC^S) &
      + 2.0d0*omega(ixC^S)*myM%g(0,phi_)%elem(ixC^S) &
      + myM%g(phi_,phi_)%elem(ixC^S)*omega(ixC^S)**2 ) )

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Get metric at inner disk radius:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  gin(:,:) = zero
  xinBL(1) = eqpar(rin_)
  {^IFZIN
  xinBL(z_) = dpi/2.0d0
  }
  {^IFPHIIN
  xinBL(phi_) = 0.0d0
  }
  call BLToCoord(1^D&,1^D&,1^D&,1^D&,xinBL,xin)
  call get_g4({xin(^D)},gin) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 rcyl2_in = gin(0,phi_)**2-gin(0,0)*gin(phi_,phi_)
 d2in     = (gin(0,0)*eqpar(l_)**2+2.0d0*gin(0,phi_)*eqpar(l_)+gin(phi_,phi_))
 if (d2in >zero .and. rcyl2_in>zero) then 
    utd_in = - sqrt(rcyl2_in/d2in)
 else
    utd_in = -1.0d0
 end if

 h(ixC^S) = utd_in / utd(ixC^S)

 ! Now fill the primitive variables:

 where (utd_in/utd(ixC^S) .gt. 1.0d0 .and. d2(ixC^S) .gt. zero .and. xC(ixC^S,1) .gt. xin(1))
    rho(ixC^S) = (gmm1 * (h(ixC^S)-1.0d0)/eqpar(kappa_))**(1.0d0/(eqpar(gamma_)-1.0d0))
 elsewhere
    rho(ixC^S) = zero
 end where
 
 select case(idir)
 case(1)
    A(ixC^S) = zero
 case(^Z)
    A(ixC^S) = zero
 case(^PHI)
    call CoordToBL(ixI^L,ixC^L,xC,xBL)

    A(ixC^S) = ((xBL(ixC^S,1)*sin(xBL(ixC^S,^Z))/eqpar(rin_))**3)*exp(-xBL(ixC^S,1)/400.d0)
    A(ixC^S) = A(ixC^S)*rho(ixC^S)/rhomax - 0.3d0 !0.25d0 !0.2d0

    where (A(ixC^S) .le. 0.0d0)
       A(ixC^S) = 0.0d0
    end where
 end select

end subroutine initvecpot_usr
  !=============================================================================
  subroutine specialvar_output(ixI^L,ixO^L,nwmax,w,s,normconv)

    ! this subroutine can be used in convert, to add auxiliary variables to the
    ! converted output file, for further analysis using tecplot, paraview, ....
    ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
    !
    ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
    ! corresponding normalization values (default value 1)

    use mod_physaux, only: get_b2, get_u4, get_b4, get_TEM4_ud, get_TPAKE4_ud, get_TEN4_ud
    use mod_metric, only: lower3, lower4, LRNFu
    use mod_transform, only: matvec, boost
    use mod_amrvacdef

    integer, intent(in)                :: ixI^L,ixO^L,nwmax
    double precision, intent(inout)    :: w(ixI^S,nwmax)
    type(state)                        :: s
    double precision                   :: normconv(0:nwmax)
    ! .. local ..
    integer                            :: ix^D
    double precision, parameter        :: thmin=dpi/3.0d0, thmax=2.0d0*dpi/3.0d0
    double precision, dimension(ixI^S) :: inDisk, rhoav, gav, tav, bav, pav, pmagav, betaav
    double precision, dimension(ixI^S) :: rho2av, t2av, b2av, p2av, pmag2av, beta2av, omega, b2
    double precision                   :: bD(ixI^S,1:^NC)
    double precision, dimension(ixI^S,0:ndir) :: u4, u4d, u4hat, u4bar, b4, b4hat, b4bar
    double precision, dimension(ixI^S,0:ndir) :: dx4, dx4hat, dx4bar
    double precision, dimension(ixI^S)        :: rhoh, t_rphi
    double precision, dimension(ixI^S,0:ndir,0:ndir) :: tem, tpake, ten, ehatu, t_tmp
    integer :: flag
    !-----------------------------------------------------------------------------
    associate(x=>s%x%x{#IFDEF STAGGERED ,ws=>s%ws%w})

    if (nwmax-nw .gt. 0) then
       call get_b2(ixI^L,ixO^L,w(ixI^S,1:nw),x,b2,w_is_primitive=.false.)
       w(ixO^S,nw+1) = b2(ixO^S)
    end if
    
    call get_u4(ixI^L,ixO^L,w(ixI^S,1:nw),x,u4,w_is_primitive=.false.)

    if (nwmax-nw .gt. 1) then
       {#IFDEF STAGGERED
       call div_staggered(ixO^L,s,w(ixO^S,nw+2))
       }{#IFNDEF STAGGERED
       ! Reduce output array size, +1 was added for eventual pointdata output
       call get_divb(ixI^L,ixO^L^LSUB1,w(ixI^S,1:nw),w(ixI^S,nw+2))
       }
    end if

    ! Qttys for theta-averaging:
    where(x(ixO^S,z_).ge.thmin .and. x(ixO^S,z_).le.thmax)
       inDisk(ixO^S) = 1.0d0
    elsewhere
       inDisk(ixO^S) = 0.0d0
    end where

    call primitive(ixI^L,ixO^L,w(ixI^S,1:nw),x)

    gav(ixO^S)   = inDisk(ixO^S) * myM%sqrtgamma(ixO^S) * myM%alpha(ixO^S)
    rhoav(ixO^S) = inDisk(ixO^S) * myM%sqrtgamma(ixO^S) * myM%alpha(ixO^S) * w(ixO^S,rho_)

    call lower3(ixI^L,ixO^L,myM,w(ixI^S,b1_:b^NC_),bD)
    bav(ixO^S)   = sqrt({^C&(w(ixO^S,b^C_)*bD(ixO^S,^C))|+})
    bav(ixO^S)   = bav(ixO^S) * inDisk(ixO^S) * myM%sqrtgamma(ixO^S) * myM%alpha(ixO^S)

    tav(ixO^S)   = w(ixO^S,pp_)/w(ixO^S,rho_)
    tav(ixO^S)   = tav(ixO^S) * inDisk(ixO^S) * myM%sqrtgamma(ixO^S) * myM%alpha(ixO^S)

    pav(ixO^S)   = w(ixO^S,pp_)
    pav(ixO^S)   = pav(ixO^S) * inDisk(ixO^S) * myM%sqrtgamma(ixO^S) * myM%alpha(ixO^S)

    pmagav(ixO^S)  = w(ixO^S,nw+1)/2.0d0
    pmagav(ixO^S)  = pmagav(ixO^S) * inDisk(ixO^S) * myM%sqrtgamma(ixO^S) * myM%alpha(ixO^S)

    betaav(ixO^S)  = 2.0d0*w(ixO^S,pp_)/w(ixO^S,nw+1)
    betaav(ixO^S)  = betaav(ixO^S) * inDisk(ixO^S) * myM%sqrtgamma(ixO^S) * myM%alpha(ixO^S)

    if (nwmax-nw .gt. 2) &
         w(ixO^S,nw+3) = gav(ixO^S)
    if (nwmax-nw .gt. 3) &
         w(ixO^S,nw+4) = rhoav(ixO^S)
    if (nwmax-nw .gt. 4) &
         w(ixO^S,nw+5) = tav(ixO^S)
    if (nwmax-nw .gt. 5) &
         w(ixO^S,nw+6) = bav(ixO^S)
    if (nwmax-nw .gt. 6) &
         w(ixO^S,nw+7) = pav(ixO^S)
    if (nwmax-nw .gt. 7) &
         w(ixO^S,nw+8) = pmagav(ixO^S)
    if (nwmax-nw .gt. 8) &
         w(ixO^S,nw+9) = betaav(ixO^S)

    if (nwmax-nw .gt. 9) &
         w(ixO^S,nw+10) = - myM%sqrtgamma(ixO^S) * myM%alpha(ixO^S) * w(ixO^S,rho_) * u4(ixO^S,1)
    if (nwmax-nw .gt. 10) &
         w(ixO^S,nw+11) = half * myM%sqrtgamma(ixO^S) * abs(w(ixO^S,b1_))
    
    if (nwmax-nw .gt. 11) then
       call get_TEM4_ud(ixI^L,ixO^L,w(ixI^S,1:nw),x,tem,w_is_primitive=.true.)
       w(ixO^S,nw+12) = - tem(ixO^S,1,0) * myM%alpha(ixO^S) * myM%sqrtgamma(ixO^S)
    end if

    if (nwmax-nw .gt. 12) then
       call get_TPAKE4_ud(ixI^L,ixO^L,w(ixI^S,1:nw),x,tpake,w_is_primitive=.true.)
       w(ixO^S,nw+13) = - tpake(ixO^S,1,0) * myM%alpha(ixO^S) * myM%sqrtgamma(ixO^S)
    end if

    if (nwmax-nw .gt. 13) then
       call get_TEN4_ud(ixI^L,ixO^L,w(ixI^S,1:nw),x,ten,w_is_primitive=.true.)
       w(ixO^S,nw+14) = - ten(ixO^S,1,0) * myM%alpha(ixO^S) * myM%sqrtgamma(ixO^S)
    end if

    if (nwmax-nw .gt. 14) then
       call Enthalpy(w(ixI^S,1:nw),ixI^L,ixO^L,patchfalse,rhoh)
       call lower4(ixI^L,ixO^L,myM,u4,u4d)
       w(ixO^S,nw+15) = - rhoh(ixO^S)/w(ixO^S,rho_)*u4d(ixO^S,0)
    end if

    if (nwmax-nw .gt. 15) then
       ! Total angular momentum flux T^r_\phi
       call get_TEM4_ud(ixI^L,ixO^L,w(ixI^S,1:nw),x,t_tmp,w_is_primitive=.true.)
       t_rphi(ixO^S) = t_tmp(ixO^S,1,phi_)
       call get_TPAKE4_ud(ixI^L,ixO^L,w(ixI^S,1:nw),x,t_tmp,w_is_primitive=.true.)
       t_rphi(ixO^S) = t_rphi(ixO^S) + t_tmp(ixO^S,1,phi_)
       call get_TEN4_ud(ixI^L,ixO^L,w(ixI^S,1:nw),x,t_tmp,w_is_primitive=.true.)
       t_rphi(ixO^S) = t_rphi(ixO^S) + t_tmp(ixO^S,1,phi_)
       
       w(ixO^S,nw+16) = t_rphi(ixO^S) * myM%alpha(ixO^S) * myM%sqrtgamma(ixO^S)
    end if

    if (nwmax-nw .gt. 16) then
       w(ixO^S,nw+17) = dble(mype)
    end if

    call conserve(ixI^L,ixO^L,w,x,patchfalse)

    if (nwmax-nw .gt. 17) then
      call get_Q_diagonal(ixI^L,ixO^L,w(ixI^S,1:nw),mygeo%xbar(ixI^S,1:ndim),w(ixI^S,nw+18:nw+20))
    end if
    
   if (nwmax-nw .gt. 20) then
     call flag_grid_usr(t,ixI^L,ixO^L,w,x,flag)
     w(ixO^S,nw+21)=flag*1.0d0
   end if

  end associate
  end subroutine specialvar_output
  !=============================================================================
  subroutine get_Qtheta(ixI^L,ixO^L,w,x,m,Q)

    use mod_physaux, only: get_b2, get_u4, get_b4
    use mod_metric, only: LRNFu
    use mod_transform, only: matvec, boost
    use mod_amrvacdef
    
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: w(ixI^S,1:nw)
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    type(metric)                       :: m
    double precision, intent(out)      :: Q(ixI^S)
    ! .. local ..
    double precision, dimension(ixI^S)       :: b2, rhoh, omega
    double precision, dimension(ixI^S,0:^NC) :: u4, b4, u4hat, b4hat, u4bar, b4bar
    double precision, dimension(ixI^S,0:^NC) :: dx4, dx4hat, dx4bar
    double precision, dimension(ixI^S,0:ndir,0:ndir) :: ehatu
    !-----------------------------------------------------------------------------
    
       call get_b2(ixI^L,ixO^L,w(ixI^S,1:nw),x,b2,w_is_primitive=.false.)
       call get_b4(ixI^L,ixO^L,w(ixI^S,1:nw),x,b4,w_is_primitive=.false.)
       call get_u4(ixI^L,ixO^L,w(ixI^S,1:nw),x,u4,w_is_primitive=.false.)
       rhoh(ixO^S) = w(ixO^S,xi_)/w(ixO^S,lfac_)**2
       dx4 = zero
       {^D& dx4(ixO^S,^D) = dxlevel(^D)\}
       ! I am not sure what to do with the phi-grid spacing in 2D.
       ! In principle it should be set to infinity...
       call LRNFu(ixI^L,ixO^L,m,ehatu)
       call matvec(ixI^L,ixO^L,u4,ehatu,u4hat)
       call matvec(ixI^L,ixO^L,b4,ehatu,b4hat)
       call matvec(ixI^L,ixO^L,dx4,ehatu,dx4hat)
       call boost(ixI^L,ixO^L,u4hat,b4hat,b4bar)
       call boost(ixI^L,ixO^L,u4hat,dx4hat,dx4bar)
       omega(ixO^S)   = u4(ixO^S,phi_)/u4(ixO^S,0)
       Q(ixO^S) = two*dpi*abs(b4bar(ixO^S,z_)/(sqrt(rhoh(ixO^S)+b2(ixO^S))*omega(ixO^S)) / dx4bar(ixO^S,z_))
       
  end subroutine get_Qtheta
  !=============================================================================
  subroutine get_Qtheta_general(ixI^L,ixO^L,w,x,Q)

    use mod_physaux, only: get_b2, get_u4, get_b4
    use mod_metric, only: CoordToKS, u4CoordToKS
    use mod_transform, only: matvec, boost, KSToLRNF
    use mod_amrvacdef
    
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: w(ixI^S,1:nw)
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision, intent(out)      :: Q(ixI^S)
    ! .. local ..
    double precision                         :: xKS(ixI^S,1:ndim)
    double precision, dimension(ixI^S)       :: b2, rhoh, omega
    double precision, dimension(ixI^S,0:^NC) :: u4, u4MKS, b4, b4MKS, u4hat, b4hat, b4bar
    double precision, dimension(ixI^S,0:^NC) :: dx4, dx4MKS, dx4hat, dx4bar
    double precision, dimension(ixI^S,0:ndir,0:ndir) :: ehatu
    !-----------------------------------------------------------------------------
    
       call get_b2(ixI^L,ixO^L,w(ixI^S,1:nw),x,b2,w_is_primitive=.false.)
       rhoh(ixO^S) = w(ixO^S,xi_)/w(ixO^S,lfac_)**2

       call get_b4(ixI^L,ixO^L,w(ixI^S,1:nw),x,b4MKS,w_is_primitive=.false.)
       call get_u4(ixI^L,ixO^L,w(ixI^S,1:nw),x,u4MKS,w_is_primitive=.false.)       
       dx4MKS = zero
       {^D& dx4MKS(ixO^S,^D) = dxlevel(^D)\}
       ! I am not sure what to do with the phi-grid spacing in 2D.
       ! In principle it should be set to infinity...

       call u4CoordToKS(ixI^L,ixO^L,x,b4MKS,b4)
       call u4CoordToKS(ixI^L,ixO^L,x,u4MKS,u4)
       call u4CoordToKS(ixI^L,ixO^L,x,dx4MKS,dx4)

       call CoordToKS(ixI^L,ixO^L,x,xKS)
       call KSToLRNF(ixI^L,ixO^L,xKS,ehatu)
       
       call matvec(ixI^L,ixO^L,u4,ehatu,u4hat)
       call matvec(ixI^L,ixO^L,b4,ehatu,b4hat)
       call matvec(ixI^L,ixO^L,dx4,ehatu,dx4hat)
       call boost(ixI^L,ixO^L,u4hat,b4hat,b4bar)
       call boost(ixI^L,ixO^L,u4hat,dx4hat,dx4bar)
       omega(ixO^S)   = u4(ixO^S,phi_)/u4(ixO^S,0)

       Q(ixO^S) = two*dpi*abs(b4bar(ixO^S,z_)/(sqrt(rhoh(ixO^S)+b2(ixO^S))*omega(ixO^S)) / dx4bar(ixO^S,z_))
       
  end subroutine get_Qtheta_general
  
!----------------------------------------------------------------------------------------------------------------------  
   subroutine get_Q_diagonal(ixI^L,ixO^L,w,x,Q)

    use mod_physaux, only: get_b2, get_u4, get_b4
    use mod_metric
    use mod_transform, only: matvec, boost, KSToLRNF
    use mod_amrvacdef
    
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: w(ixI^S,1:nw)
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision, intent(out)      :: Q(ixI^S,1:ndir)
    ! .. local ..
    double precision                         :: xKS(ixI^S,1:ndim)
    double precision, dimension(ixI^S)       :: b2, rhoh, rho, omega
    double precision, dimension(ixI^S,0:^NC) :: u4, u4MKS, b4, b4MKS, u4hat, b4hat, b4bar
    double precision, dimension(ixI^S,0:^NC) :: dx4, dx4MKS, dx4hat, dx4bar, dx4tmp
    double precision, dimension(ixI^S,0:ndir,0:ndir) :: ehatu
    !-----------------------------------------------------------------------------
    
       call get_b2(ixI^L,ixO^L,w(ixI^S,1:nw),x,b2,w_is_primitive=.false.)
!       rhoh(ixO^S) = w(ixO^S,xi_)/w(ixO^S,lfac_)**2
       rho(ixO^S) = w(ixO^S,d_)/w(ixO^S,lfac_)**2

       call get_b4(ixI^L,ixO^L,w(ixI^S,1:nw),x,b4,w_is_primitive=.false.)
       call get_u4(ixI^L,ixO^L,w(ixI^S,1:nw),x,u4,w_is_primitive=.false.)

       {^D&
       dx4(ixO^S,^D) = dxlevel(^D)
       \}

    !   {^D&
    !      dx4MKS = zero
    !      dx4MKS(ixO^S,^D) = dxlevel(^D)
    !      call u4CoordToKS(ixI^L,ixO^L,x,dx4MKS,dx4tmp)
    !      dx4(ixO^S,0:ndir) = dx4(ixO^S,0:ndir) + abs(dx4tmp(ixO^S,0:ndir))
    !   \}
       
    !   call u4CoordToKS(ixI^L,ixO^L,x,b4MKS,b4)
    !   call u4CoordToKS(ixI^L,ixO^L,x,u4MKS,u4)

    !   call CoordToKS(ixI^L,ixO^L,x,xKS)
    !   call KSToLRNF(ixI^L,ixO^L,xKS,ehatu)
    !   
    !   call matvec(ixI^L,ixO^L,u4,ehatu,u4hat)
    !   call matvec(ixI^L,ixO^L,b4,ehatu,b4hat)
    !   call matvec(ixI^L,ixO^L,dx4,ehatu,dx4hat)

       u4hat(ixO^S,0)=u4(ixO^S,0)*myM%alpha(ixO^S)
       u4hat(ixO^S,1)=u4(ixO^S,1)*sqrt(myM%g(1,1)%elem(ixO^S))
       u4hat(ixO^S,2)=u4(ixO^S,2)*sqrt(myM%g(2,2)%elem(ixO^S))
       u4hat(ixO^S,3)=u4(ixO^S,3)*sqrt(myM%g(3,3)%elem(ixO^S))

       b4hat(ixO^S,0)=b4(ixO^S,0)*myM%alpha(ixO^S)
       b4hat(ixO^S,1)=b4(ixO^S,1)*sqrt(myM%g(1,1)%elem(ixO^S))
       b4hat(ixO^S,2)=b4(ixO^S,2)*sqrt(myM%g(2,2)%elem(ixO^S))
       b4hat(ixO^S,3)=b4(ixO^S,3)*sqrt(myM%g(3,3)%elem(ixO^S))

       dx4hat(ixO^S,0)=dx4(ixO^S,0)*myM%alpha(ixO^S)
       dx4hat(ixO^S,1)=dx4(ixO^S,1)*sqrt(myM%g(1,1)%elem(ixO^S))
       dx4hat(ixO^S,2)=dx4(ixO^S,2)*sqrt(myM%g(2,2)%elem(ixO^S))
       dx4hat(ixO^S,3)=dx4(ixO^S,3)*sqrt(myM%g(3,3)%elem(ixO^S))


       call boost(ixI^L,ixO^L,u4hat,b4hat,b4bar)
       call boost(ixI^L,ixO^L,u4hat,dx4hat,dx4bar)
       omega(ixO^S)   = u4(ixO^S,phi_)/u4(ixO^S,0)

       {^C& 
!       Q(ixO^S,^C) = two*dpi*abs(b4bar(ixO^S,^C)/(sqrt(rhoh(ixO^S)+b2(ixO^S))*omega(ixO^S)) / dx4bar(ixO^S,^C)
       Q(ixO^S,^C) = two*dpi*abs(b4bar(ixO^S,^C)/(sqrt(rho(ixO^S))*omega(ixO^S))) / dx4bar(ixO^S,^C)
       \}
       
  end subroutine get_Q_diagonal

  
  !=============================================================================
  subroutine specialvarnames_output

    ! newly added variables need to be concatenated with the varnames/primnames string

    use mod_amrvacdef
    integer                            :: iw
    !-----------------------------------------------------------------------------

    write(primnames,"(a,a)") trim(primnames),' B2'
    write(wnames,"(a,a)") trim(wnames),' B2'

    write(primnames,"(a,a)") trim(primnames),' divB'
    write(wnames,"(a,a)") trim(wnames),' divB'

    write(primnames,"(a,a)") trim(primnames),' gav'
    write(wnames,"(a,a)") trim(wnames),' gav'

    write(primnames,"(a,a)") trim(primnames),' rhoav'
    write(wnames,"(a,a)") trim(wnames),' rhoav'

    write(primnames,"(a,a)") trim(primnames),' tav'
    write(wnames,"(a,a)") trim(wnames),' tav'

    write(primnames,"(a,a)") trim(primnames),' bav'
    write(wnames,"(a,a)") trim(wnames),' bav'

    write(primnames,"(a,a)") trim(primnames),' pav'
    write(wnames,"(a,a)") trim(wnames),' pav'

    write(primnames,"(a,a)") trim(primnames),' pmagav'
    write(wnames,"(a,a)") trim(wnames),' pmagav'

    write(primnames,"(a,a)") trim(primnames),' betaav'
    write(wnames,"(a,a)") trim(wnames),' betaav'

    write(primnames,"(a,a)") trim(primnames),' mdot_int'
    write(wnames,"(a,a)") trim(wnames),' mdot_int'

    write(primnames,"(a,a)") trim(primnames),' phi_int'
    write(wnames,"(a,a)") trim(wnames),' phi_int'

    write(primnames,"(a,a)") trim(primnames),' TEMrt'
    write(wnames,"(a,a)") trim(wnames),' TEMrt'

    write(primnames,"(a,a)") trim(primnames),' TPAKErt'
    write(wnames,"(a,a)") trim(wnames),' TPAKErt'

    write(primnames,"(a,a)") trim(primnames),' TENrt'
    write(wnames,"(a,a)") trim(wnames),' TENrt'

    write(primnames,"(a,a)") trim(primnames),' hut'
    write(wnames,"(a,a)") trim(wnames),' hut'

    write(primnames,"(a,a)") trim(primnames),' ldot_int'
    write(wnames,"(a,a)") trim(wnames),' ldot_int'

    write(primnames,"(a,a)") trim(primnames),' mype'
    write(wnames,"(a,a)") trim(wnames),' mype'
    
    write(primnames,"(a,a)") trim(primnames),' Q_r'
    write(wnames,"(a,a)") trim(wnames),' Q_r'
    
    write(primnames,"(a,a)") trim(primnames),' Q_theta'
    write(wnames,"(a,a)") trim(wnames),' Q_theta'
    
    write(primnames,"(a,a)") trim(primnames),' Q_phi'
    write(wnames,"(a,a)") trim(wnames),' Q_phi'
    
    write(primnames,"(a,a)") trim(primnames),' flag'
    write(wnames,"(a,a)") trim(wnames),' flag'

  end subroutine specialvarnames_output
!=============================================================================
subroutine printlog_special

use mod_amrvacdef
!-----------------------------------------------------------------------------

call mpistop("special log file undefined")

end subroutine printlog_special
!=============================================================================
subroutine process_grid_usr(igrid,level,ixI^L,ixO^L,qt,w,x)

! this subroutine is ONLY to be used for computing auxiliary variables
! which happen to be non-local (like div v), and are in no way used for
! flux computations. As auxiliaries, they are also not advanced

use mod_amrvacdef

integer, intent(in):: igrid,level,ixI^L,ixO^L
double precision, intent(in):: qt,x(ixI^S,1:ndim)
double precision, intent(inout):: w(ixI^S,1:nw)
!-----------------------------------------------------------------------------

end subroutine process_grid_usr
!=============================================================================
subroutine bc_int(level,qt,ixI^L,ixO^L,w,x)

! internal boundary, user defined
!
! This subroutine can be used to artificially overwrite ALL conservative
! variables in a user-selected region of the mesh, and thereby act as
! an internal boundary region. It is called just before external (ghost cell)
! boundary regions will be set by the BC selection. Here, you could e.g.
! want to introduce an extra variable (nwextra, to be distinguished from nwaux)
! which can be used to identify the internal boundary region location.
! Its effect should always be local as it acts on the mesh.
!

use mod_amrvacdef

integer, intent(in) :: ixI^L,ixO^L,level
double precision, intent(in) :: qt
double precision, intent(inout) :: w(ixI^S,1:nw)
double precision, intent(in) :: x(ixI^S,1:ndim)

! .. local ..
logical :: patchw(ixG^T)
!----------------------------------------------------------------------------

{#IFDEF TRACER
    ! Dominate tracer:
    where (w(ixO^S,Dtr1_) .lt. zero)
       w(ixO^S,Dtr1_) = zero
    end where
    where (w(ixO^S,Dtr1_) .gt. w(ixO^S,d_))
       w(ixO^S,Dtr1_) = w(ixO^S,d_)
    end where
}
end subroutine bc_int
!=============================================================================
subroutine init_random_seed()
use iso_fortran_env, only: int64
implicit none
integer, allocatable :: seed(:)
  integer :: i, n, un, istat, dt(8), pid
  integer(int64) :: t

  call random_seed(size = n)
  allocate(seed(n))
  ! First try if the OS provides a random number generator
  open(newunit=un, file="/dev/urandom", access="stream", &
       form="unformatted", action="read", status="old", iostat=istat)
  if (istat == 0) then
     read(un) seed
     close(un)
  else
     ! Fallback to XOR:ing the current time and pid. The PID is
     ! useful in case one launches multiple instances of the same
     ! program in parallel.
     call system_clock(t)
     if (t == 0) then
        call date_and_time(values=dt)
        t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
             + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
             + dt(3) * 24_int64 * 60 * 60 * 1000 &
             + dt(5) * 60 * 60 * 1000 &
             + dt(6) * 60 * 1000 + dt(7) * 1000 &
             + dt(8)
     end if
     pid = 1001 !getpid()
     t = ieor(t, int(pid, kind(t)))
     do i = 1, n
        seed(i) = lcg(t)
     end do
  end if
  call random_seed(put=seed)
contains
  ! This simple PRNG might not be good enough for real work, but is
  ! sufficient for seeding a better PRNG.
  function lcg(s)
    integer :: lcg
    integer(int64) :: s
    if (s == 0) then
       s = 104729
    else
       s = mod(s, 4294967296_int64)
    end if
    s = mod(s * 279470273_int64, 4294967291_int64)
    lcg = int(mod(s, int(huge(0), int64)), kind(0))
  end function lcg
end subroutine init_random_seed
!=============================================================================
subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

! Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
! iw=iwmin...iwmax.  wCT is at time qCT

use mod_amrvacdef

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

! integer :: iw
! double precision :: s(ixG^T)
!-----------------------------------------------------------------------------

! do iw= iw^LIM
!    select case(iw)
!    case(m1_)
!       ! The source is based on the time centered wCT
!       call getmyforce(wCT,ixO^L,s)
!       w(ixO^S,m1_)=w(ixO^S,m1_) + qdt*s(ixO^S)
!    case(e_)
!       call getmyheating(wCT,ixO^L,s)
!       w(ixO^S,e_) =w(ixO^S,e_)  + qdt*s(ixO^S)
!    end select
! end do

end subroutine specialsource
!=============================================================================
subroutine getdt_special(w,ixI^L,ixO^L,dtnew,dx^D,x)

! Limit "dt" further if necessary, e.g. due to the special source terms.
! The getdt_courant (CFL condition) and the getdt subroutine in the AMRVACPHYS
! module have already been called.

use mod_amrvacdef

integer, intent(in) :: ixI^L, ixO^L
double precision, intent(in) :: dx^D, x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw), dtnew
!-----------------------------------------------------------------------------

dtnew=bigdouble

end subroutine getdt_special
!=============================================================================
subroutine specialeta(w,ixI^L,ixO^L,idirmin,x,current,eta)

! Set the "eta" array for resistive MHD based on w or the
! "current" variable which has components between idirmin and 3.

use mod_amrvacdef

integer, intent(in) :: ixI^L, ixO^L, idirmin
double precision, intent(in) :: w(ixI^S,nw), x(ixI^S,1:ndim)

double precision :: current(ixG^T,7-2*ndir:3), eta(ixG^T)
!-----------------------------------------------------------------------------

!  eta(ix^S)=...

call mpistop("con2prim can only handle constant and uniform resistivity at the moment")

end subroutine specialeta
!=============================================================================
subroutine specialrefine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)

! Enforce additional refinement or coarsening
! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

! you must set consistent values for integers refine/coarsen:

! refine = -1 enforce to not refine
! refine =  0 doesn't enforce anything
! refine =  1 enforce refinement

! coarsen = -1 enforce to not coarsen
! coarsen =  0 doesn't enforce anything
! coarsen =  1 enforce coarsen

use mod_metric, only: CoordToBL
use mod_amrvacdef

integer, intent(in) :: igrid, level, ixI^L, ixO^L
double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
integer, intent(inout) :: refine, coarsen
! .. local ..
double precision       :: xBL(ixI^S,1:ndim)
!-----------------------------------------------------------------------------

call CoordToBL(ixI^L,ixO^L,x,xBL)
!if ((any(xBL(ixO^S,1).lt.8.0d0) &
!     .and. (any(xBL(ixO^S,2).lt.dpi/10.0d0) &
!     .or. any(xBL(ixO^S,2).gt.dpi-dpi/10.0d0))) &
!     .or. any(xBL(ixO^S,1).gt.100.0d0)) then
!   coarsen=1
!   refine=-1
if (any(xBL(ixO^S,1).lt.8.0d0)&
    .or. any(xBL(ixO^S,1).gt.1000.0d0))then
   coarsen=1
   refine=-1
else if (any(xBL(ixO^S,2).lt.dpi/10.0d0) &
     .or. any(xBL(ixO^S,2).gt.dpi-dpi/10.0d0)) &
      then
   select case (level)
   case(1)
       coarsen=-1
       refine=1
   case(2)
       refine=-1
       coarsen=-1
   case default
       refine=-1
       coarsen=1
   end select
else if ((any(xBL(ixO^S,1).lt.20.0d0) &
     .and. (any(xBL(ixO^S,2).lt.dpi/4.0d0) &
     .or. any(xBL(ixO^S,2).gt.dpi-dpi/4.0d0)))) then
   refine=0
   coarsen=0
else
   coarsen=0
   refine=1
end if

end subroutine specialrefine_grid
!=============================================================================
subroutine specialvarforerrest(ixI^L,ixO^L,iflag,w,var)

! this is the place to compute a local auxiliary variable to be used
! as refinement criterion for the Lohner error estimator only
!  -->it is then requiring and iflag>nw
! note that ixO=ixI=ixG, hence the term local (gradients need special attention!)

use mod_amrvacdef

integer, intent(in)          :: ixI^L,ixO^L,iflag
double precision, intent(in) :: w(ixI^S,1:nw)
double precision, intent(out):: var(ixG^T)
!-----------------------------------------------------------------------------

if (iflag >nw)call mpistop(' iflag> nw, make change in parfile or in user file')

var(ixI^S) = zero 

end subroutine specialvarforerrest
!=============================================================================
subroutine specialset_B0(ixI^L,ixO^L,x,wB0)

! Here one can add a steady (time-independent) potential background field

use mod_amrvacdef

integer, intent(in)           :: ixI^L,ixO^L
double precision, intent(in)  :: x(ixG^T,1:ndim)
double precision, intent(inout) :: wB0(ixI^S,1:ndir)
!-----------------------------------------------------------------------------
call mpistop(' abs(Busr)> 0, make change in parfile or in user file')

wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)
!!wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)+user defined steady potential field

end subroutine specialset_B0
!=============================================================================
