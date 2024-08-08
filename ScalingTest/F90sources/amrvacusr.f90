  !=============================================================================
  ! amrvacusr.t
  !=============================================================================
!  INCLUDE:amrvacnul/speciallog.t
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
subroutine specialbound_usr(qt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iB,s)

! special boundary types, user defined
! user must assign conservative variables in bounderies

use mod_amrvacdef

integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
    ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, iB
double precision, intent(in) :: qt
type(state), intent(inout)   :: s
!----------------------------------------------------------------------------
associate(x=>s%x%x,w=>s%w%w ,ws=>s%ws%w)
call mpistop("specialbound not defined")

! implement boundary conditions for a side, for example
! set the radial momentum at the inner radial boundary (assuming spherical-
! like coordinates) to zero:

!  select case(iB)
!  case(1)
!  w(ixO^S,s1_) = 0.0d0
!  end select

end associate
end subroutine specialbound_usr
!=============================================================================
subroutine bc_int(level,qt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)

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

integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,level
double precision, intent(in) :: qt
double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw)
double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndim)

! .. local ..
!logical :: patchw(ixG^T)
!----------------------------------------------------------------------------

call mpistop("bc_int not defined")

! just to give an example for relativistic MHD
!  -----------------------------------------
!patchw(ixO^S)=.true.
!where (({^D&x(ixO^S,^D)**2+})<half**2.0d0) 
!    patchw(ixO^S) = .false.
!  ^C&w(ixO^S,v^C_)=zero;
!  ^C&w(ixO^S,b^C_)=zero;
!    w(ixO^S,b3_) = one
!    w(ixO^S,v1_) = 0.99
!    w(ixO^S,rho_) = 1.d0
!    w(ixO^S,pp_)  = 2.0d0
!    w(ixO^S,lfac_)=one/dsqrt(one-({^C&w(ixO^S,v^C_)**2.0d0+}))
!end where
!!if (useprimitiveRel) then
!!  where (({^D&x(ixO^S,^D)**2+})<half**2.0d0) 
!!  {^C&w(ixO^S,u^C_)=w(ixO^S,lfac_)*w(ixO^S,v^C_);\}
!!  end where
!!endif
!call conserve(ixI^L,ixO^L,w,x,patchw)

end subroutine bc_int
!=============================================================================
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
subroutine specialsource(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iwmin,iwmax,qtC,wCT,qt,w,x)

! Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
! iw=iwmin...iwmax.  wCT is at time qCT

use mod_amrvacdef

integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
    ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, iwmin,iwmax
double precision, intent(in) :: qdt, qtC, qt, x(ixImin1:ixImax1,&
   ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
double precision, intent(inout) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw), w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nw)

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
subroutine getdt_special(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dtnew,dx1,dx2,dx3,x)

! Limit "dt" further if necessary, e.g. due to the special source terms.
! The getdt_courant (CFL condition) and the getdt subroutine in the AMRVACPHYS
! module have already been called.

use mod_amrvacdef

integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
    ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(in) :: dx1,dx2,dx3, x(ixImin1:ixImax1,&
   ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw), dtnew
!-----------------------------------------------------------------------------

dtnew=bigdouble

end subroutine getdt_special
!=============================================================================
subroutine specialeta(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,x,current,eta)

! Set the "eta" array for resistive MHD based on w or the
! "current" variable which has components between idirmin and 3.

use mod_amrvacdef

integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
    ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idirmin
double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:ndim)

double precision :: current(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
   7-2*ndir:3), eta(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
!-----------------------------------------------------------------------------

!  eta(ix^S)=...

call mpistop("con2prim can only handle constant and uniform resistivity at &
   the moment")

end subroutine specialeta
!=============================================================================
subroutine specialrefine_grid(igrid,level,ixImin1,ixImin2,ixImin3,ixImax1,&
   ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,w,x,&
   refine,coarsen)

! Enforce additional refinement or coarsening
! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

! you must set consistent values for integers refine/coarsen:

! refine = -1 enforce to not refine
! refine =  0 doesn't enforce anything
! refine =  1 enforce refinement

! coarsen = -1 enforce to not coarsen
! coarsen =  0 doesn't enforce anything
! coarsen =  1 enforce coarsen

use mod_amrvacdef

integer, intent(in) :: igrid, level, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(in) :: qt, w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:ndim)
integer, intent(inout) :: refine, coarsen
!-----------------------------------------------------------------------------

! e.g. refine for negative first coordinate x < 0 as
!
! if (any(x(ix^S,1) < zero)) refine=1

end subroutine specialrefine_grid
!=============================================================================
subroutine specialvarforerrest(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iflag,w,var)

! this is the place to compute a local auxiliary variable to be used
! as refinement criterion for the Lohner error estimator only
!  -->it is then requiring and iflag>nw
! note that ixO=ixI=ixG, hence the term local (gradients need special attention!)

use mod_amrvacdef

integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iflag
double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw)
double precision, intent(out):: var(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
!-----------------------------------------------------------------------------

if (iflag >nw)call mpistop&
   (' iflag> nw, make change in parfile or in user file')

var(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) = zero 

end subroutine specialvarforerrest
!=============================================================================
subroutine specialset_B0(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,wB0)

! Here one can add a steady (time-independent) potential background field

use mod_amrvacdef

integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(in)  :: x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
   1:ndim)
double precision, intent(inout) :: wB0(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndir)
!-----------------------------------------------------------------------------
call mpistop(' abs(Busr)> 0, make change in parfile or in user file')

wB0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:ndir)&
   =wB0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:ndir)
!!wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)+user defined steady potential field

end subroutine specialset_B0
!=============================================================================
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
subroutine specialsource_impl(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iwmin,iwmax,qtC,&
   wCT,qt,w,x)

use mod_amrvacdef

integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
    ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, iwmin,iwmax
double precision, intent(in) :: qdt, qtC, qt, x(ixImin1:ixImax1,&
   ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw), wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nw)
!-----------------------------------------------------------------------------

end subroutine specialsource_impl
!=============================================================================
subroutine getdt_impl(w,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
   ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,dtnew,dx1,dx2,dx3,x)

use mod_amrvacdef

integer, intent(in) :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,ixmin1,&
   ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
double precision, intent(in) :: dx1,dx2,dx3, x(ixGmin1:ixGmax1,&
   ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:ndim)
! note that depending on strictsmall etc, w values may change
double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
   ixGmin3:ixGmax3,1:nw), dtnew
!-----------------------------------------------------------------------------
dtnew=bigdouble

end subroutine getdt_impl
!=============================================================================
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

subroutine fixp_usr(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
   ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
use mod_amrvacdef

integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(inout)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw)
double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndim)
!----------------------------------------------------------------------------


end subroutine fixp_usr
!=============================================================================
subroutine flag_grid_usr(qt,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,flag)

use mod_amrvacdef

integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
   ixGmax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
integer, intent(inout)          :: flag
double precision, intent(in)    :: qt
double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
   ixGmin3:ixGmax3,1:nw)
double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
   ixGmin3:ixGmax3,1:ndim)

! flag=-1 : Treat all cells active, omit deactivation (onentry, default)
! flag=0  : Treat as normal domain
! flag=1  : Treat as passive, but reduce by safety belt
! flag=2  : Always treat as passive

!-----------------------------------------------------------------------------
      
end subroutine flag_grid_usr
!=============================================================================
subroutine correctaux_usr(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,patchierror,subname)

use mod_amrvacdef

integer, intent(in)            :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
integer, intent(inout)         :: patchierror(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
   ixGlo3:ixGhi3)
character(len=*), intent(in)   :: subname
double precision, intent(inout):: w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw)
double precision, intent(in)   :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndim)

! correct solution from analytic case

end subroutine correctaux_usr
!==========================================================================================
  !=============================================================================
  subroutine initglobaldata_usr

    use mod_amrvacdef
    !-----------------------------------------------------------------------------

    eqpar(gamma_) = 2.0d0 ! Adiabatic index
    eqpar(swap_)  =+1.0d0 ! muliply x with this to swap Riemann problem.

  end subroutine initglobaldata_usr
  !=============================================================================
  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,s)

    ! initialize one grid within ix^L

    use mod_amrvacdef

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
        ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    type(state)         :: s
    !-----------------------------------------------------------------------------
    associate(x=>s%x%x,w=>s%w%w ,ws=>s%ws%w)

    

       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,rho_) = 1.0d0
       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,pp_)  = 1.0d0

       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,u1_)  = 0.0d0
       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,u2_)  = 0.0d0       
       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,u3_)  = 0.0d0 

       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,b1_)  = 0.5d0
       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,b2_)  = 1.0d0
       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,b3_)  = 0.0d0

       
       ws(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,bs1_)  = 0.5d0
       
       ws(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,bs2_)  = 1.0d0
      
       
       ws(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,bs3_)  &
          = 0.0d0       
      
       

    call conserve(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,ixmin1,&
       ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x,patchfalse)
 
  end associate
end subroutine initonegrid_usr
!=============================================================================
subroutine initvecpot_usr(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
    ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, xC, A, idir)

  ! initialize the vectorpotential on the corners
  ! used by b_from_vectorpotential()


  use mod_amrvacdef

  integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3, ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, idir
  double precision, intent(in)       :: xC(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:ndim)
  double precision, intent(out)      :: A(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3)
  !-----------------------------------------------------------------------------

  ! Not needed for this setup.

  A(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = zero


end subroutine initvecpot_usr
!=============================================================================
subroutine specialvar_output(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,nwmax,w,s,normconv)

! this subroutine can be used in convert, to add auxiliary variables to the
! converted output file, for further analysis using tecplot, paraview, ....
! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
!
! the array normconv can be filled in the (nw+1:nw+nwauxio) range with 
! corresponding normalization values (default value 1)

use mod_amrvacdef

integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,nwmax
double precision                   :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nwmax)
type(state)                        :: s
double precision                   :: normconv(0:nwmax)
!-----------------------------------------------------------------------------
associate(x=>s%x%x ,ws=>s%ws%w)


call div_staggered(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,s,&
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+1))


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
subroutine process_grid_usr(igrid,level,ixImin1,ixImin2,ixImin3,ixImax1,&
   ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,w,x)

! this subroutine is ONLY to be used for computing auxiliary variables
! which happen to be non-local (like div v), and are in no way used for
! flux computations. As auxiliaries, they are also not advanced

use mod_amrvacdef

integer, intent(in):: igrid,level,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(in):: qt,x(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndim)
double precision, intent(inout):: w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw)
!-----------------------------------------------------------------------------

end subroutine process_grid_usr
!=============================================================================
! amrvacusr.t
!=============================================================================
