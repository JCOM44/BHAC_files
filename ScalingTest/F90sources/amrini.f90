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
subroutine initlevelone

use mod_amrvacdef

integer :: iigrid, igrid
!-----------------------------------------------------------------------------
levmin=1
levmax=1

call init_forest_root

call getigrids
call build_connectivity

! fill solution space of all root grids
!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
do iigrid=1,igridstail; igrid=igrids(iigrid);
   call alloc_node(igrid)
   call initial_condition(igrid)
end do
!$OMP END PARALLEL DO

call selectgrids

end subroutine initlevelone
!=============================================================================
subroutine initial_condition(igrid)

! Need only to set the mesh values (can leave ghost cells untouched)

use mod_amrvacdef

integer, intent(in) :: igrid
! .. local ..

integer                                       :: ixGsmin1,ixGsmin2,ixGsmin3,&
   ixGsmax1,ixGsmax2,ixGsmax3 

external initonegrid_usr
!----------------------------------------------------------------------------
call set_tmpGlobals(igrid)

pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nw)=zero

ixGsmin1=pws(igrid)%ixGmin1;ixGsmax1=pws(igrid)%ixGmax1
ixGsmin2=pws(igrid)%ixGmin2;ixGsmax2=pws(igrid)%ixGmax2
ixGsmin3=pws(igrid)%ixGmin3;ixGsmax3=pws(igrid)%ixGmax3;
pws(igrid)%w(ixGsmin1:ixGsmax1,ixGsmin2:ixGsmax2,ixGsmin3:ixGsmax3,1:nws)=zero


call initonegrid_usr(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixMlo1,ixMlo2,&
   ixMlo3,ixMhi1,ixMhi2,ixMhi3,ps(igrid))

end subroutine initial_condition
!=============================================================================
subroutine modify_IC

use mod_amrvacdef

integer :: iigrid, igrid

external initonegrid_usr
!-----------------------------------------------------------------------------
!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
do iigrid=1,igridstail; igrid=igrids(iigrid);
   call set_tmpGlobals(igrid)
   call initonegrid_usr(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixMlo1,&
      ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,ps(igrid))
end do
!$OMP END PARALLEL DO

end subroutine modify_IC
!=============================================================================
