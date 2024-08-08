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
subroutine getbc(time,psuse,psuseCo)

use mod_amrvacdef

double precision, intent(in)                :: time
type(state), dimension(ngridshi), target    :: psuse, psuseCo

! ... local ...
double precision :: time_bcin
integer :: igrid, iigrid
!-----------------------------------------------------------------------------
time_bcin=MPI_WTIME()

if (internalboundary) then 
   call getintbc(time,psuse)
end if

! srl communications
call init_comm_gc(time,psuse,psuseCo)
do iigrid=1,igridstail; igrid=igrids(iigrid);
  call send_gc_srl(igrid) ! OMP: not threadsafe
end do

call fill_gc_srl ! OMP inside

! Restrict
!$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(igrid)
do iigrid=1,igridstail; igrid=igrids(iigrid);
  call gc_restrict(igrid) ! OMP: imbalance
end do
!$OMP END PARALLEL DO

do iigrid=1,igridstail; igrid=igrids(iigrid);
   call send_gc_r(igrid) ! OMP: not threadsafe
end do

call fill_gc_r

! Prolong
do iigrid=1,igridstail; igrid=igrids(iigrid);
  call send_gc_p(igrid) ! OMP: not threadsafe
end do

call fill_gc_p
!$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(igrid)
do iigrid=1,igridstail; igrid=igrids(iigrid);
  call gc_prolong(igrid) ! OMP: imbalance
end do
!$OMP END PARALLEL DO

! Physical boundary conditions
if(.not.energyonly) then
!$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(igrid)
  do iigrid=1,igridstail; igrid=igrids(iigrid);
    call fill_boundary(igrid) ! OMP: imbalance
  end do
!$OMP END PARALLEL DO
end if

if (nwaux>0) then
!$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call fix_auxiliary(igrid)
   end do
!$OMP END PARALLEL DO
end if
time_bc=time_bc+(MPI_WTIME()-time_bcin)

end subroutine getbc
!=============================================================================

