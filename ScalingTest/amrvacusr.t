  !=============================================================================
  ! amrvacusr.t
  !=============================================================================
!  INCLUDE:amrvacnul/speciallog.t
  INCLUDE:amrvacnul/specialbound.t
  INCLUDE:amrvacnul/specialsource.t
  INCLUDE:amrvacnul/specialimpl.t
  INCLUDE:amrvacnul/usrflags.t
  INCLUDE:amrvacnul/correctaux_usr.t
  !=============================================================================
  subroutine initglobaldata_usr

    use mod_amrvacdef
    !-----------------------------------------------------------------------------

    eqpar(gamma_) = 2.0d0 ! Adiabatic index
    eqpar(swap_)  =+1.0d0 ! muliply x with this to swap Riemann problem.

  end subroutine initglobaldata_usr
  !=============================================================================
  subroutine initonegrid_usr(ixG^L,ix^L,s)

    ! initialize one grid within ix^L

    use mod_amrvacdef

    integer, intent(in) :: ixG^L, ix^L
    type(state)         :: s
    !-----------------------------------------------------------------------------
    associate(x=>s%x%x,w=>s%w%w{#IFDEF STAGGERED ,ws=>s%ws%w})

    {#IFDEF PSI
    w(ix^S,psi_) = 0.0d0
    }

       w(ixG^S,rho_) = 1.0d0
       w(ixG^S,pp_)  = 1.0d0

       w(ixG^S,u1_)  = 0.0d0
       w(ixG^S,u2_)  = 0.0d0       
       w(ixG^S,u3_)  = 0.0d0 

       w(ixG^S,b1_)  = 0.5d0
       w(ixG^S,b2_)  = 1.0d0
       w(ixG^S,b3_)  = 0.0d0

       {#IFDEF STAGGERED
       ws(ixG^S,bs1_)  = 0.5d0
       {^NOONED
       ws(ixG^S,bs2_)  = 1.0d0
       }
       {^IFTHREED
       ws(ixG^S,bs3_)  = 0.0d0       
       }
       }

    call conserve(ixG^L,ix^L,w,x,patchfalse)
 
  end associate
end subroutine initonegrid_usr
!=============================================================================
subroutine initvecpot_usr(ixI^L, ixC^L, xC, A, idir)

  ! initialize the vectorpotential on the corners
  ! used by b_from_vectorpotential()


  use mod_amrvacdef

  integer, intent(in)                :: ixI^L, ixC^L, idir
  double precision, intent(in)       :: xC(ixI^S,1:ndim)
  double precision, intent(out)      :: A(ixI^S)
  !-----------------------------------------------------------------------------

  ! Not needed for this setup.

  A(ixC^S) = zero


end subroutine initvecpot_usr
!=============================================================================
subroutine specialvar_output(ixI^L,ixO^L,nwmax,w,s,normconv)

! this subroutine can be used in convert, to add auxiliary variables to the
! converted output file, for further analysis using tecplot, paraview, ....
! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
!
! the array normconv can be filled in the (nw+1:nw+nwauxio) range with 
! corresponding normalization values (default value 1)

use mod_amrvacdef

integer, intent(in)                :: ixI^L,ixO^L,nwmax
double precision                   :: w(ixI^S,1:nwmax)
type(state)                        :: s
double precision                   :: normconv(0:nwmax)
!-----------------------------------------------------------------------------
associate(x=>s%x%x{#IFDEF STAGGERED ,ws=>s%ws%w})

{#IFDEF STAGGERED
call div_staggered(ixO^L,s,w(ixO^S,nw+1))
}{#IFNDEF STAGGERED
! Reduce output array size, +1 was added for eventual pointdata output
call get_divb(ixI^L,ixO^L^LSUB1,w(ixI^S,1:nw),w(ixI^S,nw+1))
}

end associate
end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables to be concatenated with the primnames/wnames string

use mod_amrvacdef
!-----------------------------------------------------------------------------

primnames= TRIM(primnames)//' '//'divB'
wnames=TRIM(wnames)//' '//'divB'

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
! amrvacusr.t
!=============================================================================
