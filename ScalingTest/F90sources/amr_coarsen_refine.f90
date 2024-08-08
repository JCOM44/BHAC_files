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
subroutine amr_coarsen_refine
  use mod_forest
  
  use mod_amr_fct
  
use mod_amrvacdef

integer :: iigrid, igrid, ipe, igridCo, ipeCo, level, ic1,ic2,ic3
integer, dimension(2,2,2) :: igridFi, ipeFi
type(tree_node_ptr) :: tree, sibling
logical             :: active
integer, external :: getnode
!-----------------------------------------------------------------------------
if(addmpibarrier) call MPI_BARRIER(icomm,ierrmpi)
call proper_nesting


call store_faces
call comm_faces

! to save memory: first coarsen then refine
irecv=0
isend=0
allocate(recvstatus(MPI_STATUS_SIZE,ngridshi),recvrequest(ngridshi),&
    sendstatus(MPI_STATUS_SIZE,ngridshi),sendrequest(ngridshi))
recvrequest=MPI_REQUEST_NULL
sendrequest=MPI_REQUEST_NULL


allocate(recvstatus_stg(MPI_STATUS_SIZE,ngridshi*3),recvrequest_stg&
   (ngridshi*3), &
         sendstatus_stg(MPI_STATUS_SIZE,ngridshi*3),sendrequest_stg&
            (ngridshi*3))
recvrequest_stg=MPI_REQUEST_NULL
sendrequest_stg=MPI_REQUEST_NULL


do ipe=0,npe-1
   do igrid=1,ngridshi
      if (coarsen(igrid,ipe)) then
         if (.not.associated(igrid_to_node(igrid,ipe)%node)) cycle

         tree%node => igrid_to_node(igrid,ipe)%node%parent%node
         do ic3=1,2
         do ic2=1,2
         do ic1=1,2
            sibling%node => tree%node%child(ic1,ic2,ic3)%node
            ipeFi(ic1,ic2,ic3)=sibling%node%ipe
            igridFi(ic1,ic2,ic3)=sibling%node%igrid
         end do
         end do
         end do

         ipeCo=ipeFi(1,1,1)
         igridCo=getnode(ipeCo)

         call coarsen_tree_leaf(igridCo,ipeCo,igridFi,ipeFi,active)

         call coarsen_grid_siblings(igridCo,ipeCo,igridFi,ipeFi,active)

         ! local coarsening done
         do ic3=1,2
         do ic2=1,2
         do ic1=1,2
            if (ipeFi(ic1,ic2,ic3)==ipeCo) then
               call putnode(igridFi(ic1,ic2,ic3),ipeFi(ic1,ic2,ic3))
               coarsen(igridFi(ic1,ic2,ic3),ipeFi(ic1,ic2,ic3))=.false.
            end if
         end do
         end do
         end do
      end if
   end do
end do

if (irecv>0) call MPI_WAITALL(irecv,recvrequest,recvstatus,ierrmpi)

if (irecv>0) call MPI_WAITALL(irecv,recvrequest_stg,recvstatus_stg,ierrmpi)


if (isend>0) call MPI_WAITALL(isend,sendrequest,sendstatus,ierrmpi)

if (isend>0) call MPI_WAITALL(isend,sendrequest_stg,sendstatus_stg,ierrmpi)


! non-local coarsening done
do ipe=0,npe-1
   do igrid=1,ngridshi
      if (coarsen(igrid,ipe)) then
         if (ipe==mype) call dealloc_node(igrid)
         call putnode(igrid,ipe)
         coarsen(igrid,ipe)=.false.
      end if
   end do
end do

do ipe=0,npe-1
   do igrid=1,ngridshi
      if (refine(igrid,ipe)) then

         do ic3=1,2
         do ic2=1,2
         do ic1=1,2
            igridFi(ic1,ic2,ic3)=getnode(ipe)
            ipeFi(ic1,ic2,ic3)=ipe
         end do
         end do
         end do

         call refine_tree_leaf(igridFi,ipeFi,igrid,ipe,active)

         if (ipe==mype) call refine_grid(igridFi,ipeFi,igrid,ipe,active)

         ! refinement done
         call putnode(igrid,ipe)
         refine(igrid,ipe)=.false.
      end if
   end do
end do


call end_comm_faces


call get_level_range

! Update sfc array: igrid and ipe info in space filling curve
call amr_Morton_order()

call load_balance
deallocate(recvstatus,recvrequest,sendstatus,sendrequest)

deallocate(recvstatus_stg,recvrequest_stg,sendstatus_stg,sendrequest_stg)

! Rebuild tree connectivity
call getigrids
call build_connectivity

! Update the list of active grids
call selectgrids
!  grid structure now complete again.

! since we only filled mesh values, and advance assumes filled
! ghost cells, do boundary filling for the new levels
if (time_advance) then
   call getbc(t+dt,ps,psCoarse)
else
   call getbc(t,ps,psCoarse)
end if



if(addmpibarrier) call MPI_BARRIER(icomm,ierrmpi)
end subroutine amr_coarsen_refine
!=============================================================================
subroutine proper_nesting
use mod_forest
use mod_amrvacdef

!! following alternative
!!logical, dimension(1:ngridshi,0:npe-1):: coarsen2

integer :: iigrid, igrid, level, ic1,ic2,ic3, inp1,inp2,inp3, i1,i2,i3,&
    my_neighbor_type,ipe
logical :: coarsening, pole(ndim)
type(tree_node_ptr) :: tree, p_neighbor, my_parent, sibling, my_neighbor,&
    neighborchild
!-----------------------------------------------------------------------------

! For all grids on all processors, do a check on refinement flags. Make
! sure that neighbors will not differ more than one level of refinement.

if (nbufferx1/=0.or.nbufferx2/=0.or.nbufferx3/=0) then
   call MPI_ALLREDUCE(MPI_IN_PLACE,refine,ngridshi*npe,MPI_LOGICAL,MPI_LOR,&
       icomm,ierrmpi)
else
   call MPI_ALLGATHER(MPI_IN_PLACE,ngridshi,MPI_LOGICAL,refine,ngridshi,&
       MPI_LOGICAL,icomm,ierrmpi)
end if

do level=min(levmax,mxnest-1),levmin+1,-1
   tree%node => level_head(level)%node
   do
      if (.not.associated(tree%node)) exit

      if (refine(tree%node%igrid,tree%node%ipe)) then
         ic1=1+modulo(tree%node%ig1-1,2);ic2=1+modulo(tree%node%ig2-1,2)
         ic3=1+modulo(tree%node%ig3-1,2);
         do inp3=ic3-2,ic3-1
         do inp2=ic2-2,ic2-1
         do inp1=ic1-2,ic1-1
            if (inp1==0.and.inp2==0.and.inp3==0) cycle
            p_neighbor%node => tree%node%parent%node
            if (inp1/=0) then
               p_neighbor%node => p_neighbor%node%neighbor(ic1,1)%node
               if (.not.associated(p_neighbor%node)) cycle
            end if
            if (inp2/=0) then
               p_neighbor%node => p_neighbor%node%neighbor(ic2,2)%node
               if (.not.associated(p_neighbor%node)) cycle
            end if
            if (inp3/=0) then
               p_neighbor%node => p_neighbor%node%neighbor(ic3,3)%node
               if (.not.associated(p_neighbor%node)) cycle
            end if
            if (p_neighbor%node%leaf) then
               refine(p_neighbor%node%igrid,p_neighbor%node%ipe)=.true.
            end if
         end do
         end do
         end do
      end if

      tree%node => tree%node%next%node
   end do
end do

! On each processor locally, check if grids set for coarsening are already
! set for refinement.

do iigrid=1,igridstail; igrid=igrids(iigrid);
   if (refine(igrid,mype).and.coarsen(igrid,mype)) coarsen(igrid,mype)=.false.
end do

! For all grids on all processors, do a check on coarsen flags.

!!coarsen2(1:ngridshi,:)=coarsen(1:ngridshi,:)

!!do ipe=0,npe-1
!! if(ipe/=mype) coarsen2(1:ngridshi,ipe)=.false.
!!enddo

call MPI_ALLGATHER(MPI_IN_PLACE,ngridshi,MPI_LOGICAL,coarsen,ngridshi,&
    MPI_LOGICAL,icomm,ierrmpi)
!!call MPI_ALLREDUCE(MPI_IN_PLACE,coarsen2,ngridshi*npe,MPI_LOGICAL,MPI_LOR, &
!!                      icomm,ierrmpi)

!!coarsen(1:ngridshi,:)=coarsen2(1:ngridshi,:)

do level=levmax,max(2,levmin),-1
   tree%node => level_head(level)%node
   do
      if (.not.associated(tree%node)) exit

      if (coarsen(tree%node%igrid,tree%node%ipe)) then
         coarsening=.true.
         my_parent%node => tree%node%parent%node

         ! are all siblings flagged for coarsen ?
check1:  do ic3=1,2
do ic2=1,2
do ic1=1,2
            sibling%node => my_parent%node%child(ic1,ic2,ic3)%node
            if (sibling%node%leaf) then
               if (coarsen(sibling%node%igrid,sibling%node%ipe)) cycle
            end if
            call unflag_coarsen_siblings
            exit check1
         end do
         end do
         end do check1

         ! Make sure that neighbors will not differ more than one level of
         ! refinement, otherwise unflag all siblings
         if (coarsening) then
check2:     do ic3=1,2
do ic2=1,2
do ic1=1,2
               sibling%node => my_parent%node%child(ic1,ic2,ic3)%node
               do i3=ic3-2,ic3-1
               do i2=ic2-2,ic2-1
               do i1=ic1-2,ic1-1
                  if (i1==0.and.i2==0.and.i3==0) cycle
                  call find_neighbor(my_neighbor,my_neighbor_type, sibling,i1,&
                     i2,i3,pole)
                  select case (my_neighbor_type)
                  case (3)
                     if (refine(my_neighbor%node%igrid, my_neighbor%node%ipe))&
                         then
                        call unflag_coarsen_siblings
                        exit check2
                     else
                        cycle
                     end if
                  case (4)
                     neighborchild%node=>my_neighbor%node%child(1,1,1)%node
                     if (neighborchild%node%leaf) then
                        if (coarsen(neighborchild%node%igrid,&
                            neighborchild%node%ipe)) then
                           cycle
                        end if
                     end if
                     call unflag_coarsen_siblings
                     exit check2
                  end select
               end do
               end do
               end do
            end do
            end do
            end do check2
         end if

      end if

      tree%node => tree%node%next%node
   end do
end do

contains
!=============================================================================
! internal procedures
!=============================================================================
subroutine unflag_coarsen_siblings

integer :: ic1,ic2,ic3
type(tree_node_ptr) :: sibling
!-----------------------------------------------------------------------------

do ic3=1,2
do ic2=1,2
do ic1=1,2
   sibling%node => my_parent%node%child(ic1,ic2,ic3)%node
   if (sibling%node%leaf) then
      coarsen(sibling%node%igrid,sibling%node%ipe)=.false.
   end if
end do
end do
end do
coarsening=.false.

end subroutine unflag_coarsen_siblings
!=============================================================================
! end of internal procedures
!=============================================================================
end subroutine proper_nesting
!=============================================================================
