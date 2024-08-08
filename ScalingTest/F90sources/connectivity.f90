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
subroutine get_level_range
use mod_forest
use mod_amrvacdef

integer :: level
!-----------------------------------------------------------------------------

! determine new finest level
do level=mxnest,1,-1
   if (associated(level_tail(level)%node)) then
      levmax=level
      exit
   end if
end do

! determine coarsest level
do level=1,levmax
   if (associated(level_tail(level)%node)) then
      levmin=level
      exit
   end if
end do

end subroutine get_level_range
!=============================================================================
subroutine getigrids
use mod_indices
use mod_connectivity
use mod_forest
implicit none

integer :: iigrid, igrid
!-----------------------------------------------------------------------------
iigrid=0
do igrid=1,ngridshi
   if (igrid_inuse(igrid,mype)) then
      iigrid=iigrid+1
      igrids(iigrid)=igrid
   end if
end do

igridstail=iigrid

end subroutine getigrids
!=============================================================================
subroutine build_connectivity
use mod_forest
use mod_amrvacdef

integer :: iigrid, igrid, i1,i2,i3, my_neighbor_type
integer :: iside, idim, ic1,ic2,ic3, inc1,inc2,inc3, ih1,ih2,ih3, icdim
type(tree_node_ptr) :: tree, my_neighbor, child
logical, dimension(3) :: pole
logical :: nopole

! Variables to detect special corners
integer :: idir,pi1,pi2,pi3, mi1,mi2,mi3, ph1,ph2,ph3, mh1,mh2,mh3,&
    ipe_neighbor

!-----------------------------------------------------------------------------
nrecv_bc_srl_13=0; nsend_bc_srl_13=0
nrecv_bc_srl_2=0; nsend_bc_srl_2=0
nrecv_bc_r_13=0; nsend_bc_r_13=0
nrecv_bc_r_2=0; nsend_bc_r_2=0
nrecv_bc_r=0; nsend_bc_r=0
nrecv_bc_p=0; nsend_bc_p=0
nrecv_fc=0; nsend_fc=0

nrecv_cc=0; nsend_cc=0

do iigrid=1,igridstail; igrid=igrids(iigrid);
   tree%node => igrid_to_node(igrid,mype)%node

   do i3=-1,1
   do i2=-1,1
   do i1=-1,1
      ! skip the grid itself
      if (i1==0.and.i2==0.and.i3==0) then
         neighbor_type(0,0,0,igrid)=0
         neighbor(1,0,0,0,igrid)=igrid
         neighbor(2,0,0,0,igrid)=mype
      else
         call find_neighbor(my_neighbor,my_neighbor_type,tree,i1,i2,i3,pole)
         nopole=.not.any(pole)

         select case (my_neighbor_type)
         ! adjacent to physical boundary
         case (1)
            neighbor(1,i1,i2,i3,igrid)=0
            neighbor(2,i1,i2,i3,igrid)=-1
         ! fine-coarse transition
         case (2)
            neighbor(1,i1,i2,i3,igrid)=my_neighbor%node%igrid
            neighbor(2,i1,i2,i3,igrid)=my_neighbor%node%ipe
            if (my_neighbor%node%ipe/=mype) then
               ic1=1+modulo(tree%node%ig1-1,2)
               ic2=1+modulo(tree%node%ig2-1,2)
               ic3=1+modulo(tree%node%ig3-1,2);
               if ((i1==0.or.i1==2*ic1-3).and.(i2==0.or.i2==2*ic2-3)&
                  .and.(i3==0.or.i3==2*ic3-3)) then
                  nrecv_bc_p=nrecv_bc_p+1
                  nsend_bc_r=nsend_bc_r+1
                 if (abs(i1)+abs(i2)+abs(i3).eq.2) then
                    nsend_bc_r_2=nsend_bc_r_2+1
                 else
                    nsend_bc_r_13=nsend_bc_r_13+1
                 end if
               end if
            end if
         ! same refinement level
         case (3)
            neighbor(1,i1,i2,i3,igrid)=my_neighbor%node%igrid
            neighbor(2,i1,i2,i3,igrid)=my_neighbor%node%ipe
            if (my_neighbor%node%ipe/=mype) then
               if (abs(i1)+abs(i2)+abs(i3).eq.2) then
                 nrecv_bc_srl_2=nrecv_bc_srl_2+1
                 nsend_bc_srl_2=nsend_bc_srl_2+1
               else
                 nrecv_bc_srl_13=nrecv_bc_srl_13+1
                 nsend_bc_srl_13=nsend_bc_srl_13+1
               end if
            end if
         ! coarse-fine transition
         case (4)
            neighbor(1,i1,i2,i3,igrid)=0
            neighbor(2,i1,i2,i3,igrid)=-1
            do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
               inc3=2*i3+ic3
               if (pole(3)) then
                  ih3=3-ic3
               else
                  ih3=ic3
               end if
            do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
               inc2=2*i2+ic2
               if (pole(2)) then
                  ih2=3-ic2
               else
                  ih2=ic2
               end if
            do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
               inc1=2*i1+ic1
               if (pole(1)) then
                  ih1=3-ic1
               else
                  ih1=ic1
               end if
               child%node => my_neighbor%node%child(ih1,ih2,ih3)%node
               neighbor_child(1,inc1,inc2,inc3,igrid)=child%node%igrid
               neighbor_child(2,inc1,inc2,inc3,igrid)=child%node%ipe
               if (child%node%ipe/=mype) then
                  nrecv_bc_r=nrecv_bc_r+1
                  if (abs(i1)+abs(i2)+abs(i3).eq.2) then
                     nrecv_bc_r_2=nrecv_bc_r_2+1
                  else
                     nrecv_bc_r_13=nrecv_bc_r_13+1
                  end if
                  nsend_bc_p=nsend_bc_p+1
               end if
            end do
            end do
            end do
         end select

         ! flux fix for conservation only for pure directional shifts
         if (abs(i1)+abs(i2)+abs(i3)==1) then
            if (i1/=0) then
               idim=1
               iside=int((i1+3)/2)
            end if
            if (i2/=0) then
               idim=2
               iside=int((i2+3)/2)
            end if
            if (i3/=0) then
               idim=3
               iside=int((i3+3)/2)
            end if
            select case (my_neighbor_type)
            ! only across fine-coarse or coarse-fine boundaries
            case (2)
               if (my_neighbor%node%ipe/=mype) then
                  if (.not.pole(idim)) nsend_fc(idim)=nsend_fc(idim)+1
               end if
            case (4)
               if (pole(idim)) then
                  icdim=iside
               else
                  icdim=3-iside
               end if
               select case (idim)
               case (1)
                  do ic1=icdim,icdim
               do ic2=1,2
               do ic3=1,2
                     child%node => my_neighbor%node%child(ic1,ic2,ic3)%node
                     if (child%node%ipe/=mype) then
                        if (.not.pole(1)) nrecv_fc(1)=nrecv_fc(1)+1
                     end if
                  end do
               end do
               end do 
               case (2)
                  do ic1=1,2
               do ic2=icdim,icdim
               do ic3=1,2
                     child%node => my_neighbor%node%child(ic1,ic2,ic3)%node
                     if (child%node%ipe/=mype) then
                        if (.not.pole(2)) nrecv_fc(2)=nrecv_fc(2)+1
                     end if
                  end do
               end do
               end do 
               case (3)
                  do ic1=1,2
               do ic2=1,2
               do ic3=icdim,icdim
                     child%node => my_neighbor%node%child(ic1,ic2,ic3)%node
                     if (child%node%ipe/=mype) then
                        if (.not.pole(3)) nrecv_fc(3)=nrecv_fc(3)+1
                     end if
                  end do
               end do
               end do 
               end select
            end select
         end if

         neighbor_pole(i1,i2,i3,igrid)=0
         
         if (my_neighbor_type>1) then
            do idim=1,3
               if (pole(idim)) then
                  neighbor_pole(i1,i2,i3,igrid)=idim
                  exit ! there can only be one pole between two meshes
               end if
            end do
         end if
         neighbor_type(i1,i2,i3,igrid)=my_neighbor_type

      end if
   end do
   end do
   end do


!  Now all the neighbour information is known.
!  Check if there are special corners that need to be communicated
!  To determine whether to send/receive, we must check three neighbours
   do i3=-1,1
do i2=-1,1
do i1=-1,1
    if (abs(i1)+abs(i2)+abs(i3)==1) then

      if (neighbor_pole(i1,i2,i3,igrid).ne.0) cycle

     ! Assign value to idim and iside
     if (i1/=0) then
         idim=1
         iside=int((i1+3)/2)
      end if
if (i2/=0) then
         idim=2
         iside=int((i2+3)/2)
      end if
if (i3/=0) then
         idim=3
         iside=int((i3+3)/2)
      end if
      ! Fine block surrounded by coarse blocks
      if (neighbor_type(i1,i2,i3,igrid).eq.2) then
        do idir=idim+1,ndim
          pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
          mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
          ph1=pi1-kr(idim,1)*(2*iside-3);ph2=pi2-kr(idim,2)*(2*iside-3)
          ph3=pi3-kr(idim,3)*(2*iside-3);
          mh1=mi1-kr(idim,1)*(2*iside-3);mh2=mi2-kr(idim,2)*(2*iside-3)
          mh3=mi3-kr(idim,3)*(2*iside-3);

          if (neighbor_type(pi1,pi2,pi3,igrid).eq.2.and.&
              neighbor_type(ph1,ph2,ph3,igrid).eq.2.and.&
              mype.ne.neighbor(2,pi1,pi2,pi3,igrid)) then
       if (neighbor_pole(pi1,pi2,pi3,igrid).eq.0) &
                nsend_cc(idim) = nsend_cc(idim) + 1
          end if

          if (neighbor_type(mi1,mi2,mi3,igrid).eq.2.and.&
              neighbor_type(mh1,mh2,mh3,igrid).eq.2.and.&
              mype.ne.neighbor(2,mi1,mi2,mi3,igrid)) then
       if (neighbor_pole(mi1,mi2,mi3,igrid).eq.0) &
                nsend_cc(idim) = nsend_cc(idim) + 1
          end if
        end do
      end if
      ! Coarse block diagonal to fine block(s)
      if (neighbor_type(i1,i2,i3,igrid).eq.3) then
        do idir=idim+1,ndim
          pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
          mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
          ph1=pi1-kr(idim,1)*(2*iside-3);ph2=pi2-kr(idim,2)*(2*iside-3)
          ph3=pi3-kr(idim,3)*(2*iside-3);
          mh1=mi1-kr(idim,1)*(2*iside-3);mh2=mi2-kr(idim,2)*(2*iside-3)
          mh3=mi3-kr(idim,3)*(2*iside-3);

          if (neighbor_type(pi1,pi2,pi3,igrid).eq.4.and.&
              neighbor_type(ph1,ph2,ph3,igrid).eq.3) then
           ! Loop on children (several in 3D)
    if (neighbor_pole(pi1,pi2,pi3,igrid).eq.0) then
           do ic3=1+int((1-pi3)/2),2-int((1+pi3)/2)
               inc3=2*pi3+ic3
do ic2=1+int((1-pi2)/2),2-int((1+pi2)/2)
               inc2=2*pi2+ic2
do ic1=1+int((1-pi1)/2),2-int((1+pi1)/2)
               inc1=2*pi1+ic1
               if (mype.ne.neighbor_child(2,inc1,inc2,inc3,igrid)) then
                 nrecv_cc(idim) = nrecv_cc(idim) + 1
               end if
           end do
end do
end do
   end if
          end if

          if (neighbor_type(mi1,mi2,mi3,igrid).eq.4.and.&
              neighbor_type(mh1,mh2,mh3,igrid).eq.3) then
           ! Loop on children (several in 3D)
    if (neighbor_pole(mi1,mi2,mi3,igrid).eq.0) then
           do ic3=1+int((1-mi3)/2),2-int((1+mi3)/2)
               inc3=2*mi3+ic3
do ic2=1+int((1-mi2)/2),2-int((1+mi2)/2)
               inc2=2*mi2+ic2
do ic1=1+int((1-mi1)/2),2-int((1+mi1)/2)
               inc1=2*mi1+ic1
               if (mype.ne.neighbor_child(2,inc1,inc2,inc3,igrid)) then
                 nrecv_cc(idim) = nrecv_cc(idim) + 1
               end if
           end do
end do
end do
   end if
          end if

        end do
      end if
    end if
   end do
end do
end do  

end do

end subroutine build_connectivity
!=============================================================================
