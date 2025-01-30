!=============================================================================
subroutine settree

use mod_amrvacdef

! create and initialize grids on all levels > 1. On entry, all
! level=1 grids have been formed and initialized. This subroutine
! creates and initializes new level grids

integer :: igrid, iigrid, levnew
type(walloc) :: pwtmp
!----------------------------------------------------------------------------
! when only one level allowed, there is nothing to do anymore
if (mxnest == 1) return

call getbc(t,ps,psCoarse)
do levnew=2,mxnest+1
   if (errorestimate==1.or.errorestimate==2) then
      call setdt
      call advance(0)
   end if

   call errest

   if (errorestimate==1.or.errorestimate==2) then
      do iigrid=1,igridstail; igrid=igrids(iigrid);
         pwtmp%w => pwold(igrid)%w
         pwold(igrid)%w => pw(igrid)%w
         pw(igrid)%w => pwtmp%w
      end do
   end if

   call amr_coarsen_refine
   
   if (.not.resetgrid) then
     ! if no finer level grids created: exit
     if (levmax/=levnew) exit
   end if
end do

end subroutine settree
!=============================================================================
subroutine resettree
{#IFDEF STAGGERED
use mod_amr_fct}

use mod_amrvacdef
!-----------------------------------------------------------------------------
if (levmax>levmin) call deallocateBflux
{#IFDEF STAGGERED
call deallocateBfaces}

      call errest

      call amr_coarsen_refine

! set up boundary flux conservation arrays
if (levmax>levmin) call allocateBflux

end subroutine resettree
!=============================================================================
subroutine resettree_convert

use mod_amrvacdef
integer  :: igrid,iigrid, my_levmin, my_levmax
!-----------------------------------------------------------------------------
if (level_io > 0) then
   my_levmin = level_io
   my_levmax = level_io
else
   my_levmin = max(1,level_io_min)
   my_levmax = min(mxnest,level_io_max)
end if


do while(levmin<my_levmin.or.levmax>my_levmax)
 call getbc(t,ps,psCoarse)
 do iigrid=1,igridstail; igrid=igrids(iigrid);
    call forcedrefine_grid_io(igrid,pw(igrid)%w)
 end do

 call amr_coarsen_refine
end do
! set up boundary flux conservation arrays

end subroutine resettree_convert
!=============================================================================
