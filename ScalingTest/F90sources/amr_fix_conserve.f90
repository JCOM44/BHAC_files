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
subroutine init_comm_fix_conserve(idimmin,idimmax)
use mod_fix_conserve
use mod_amrvacdef

integer, intent(in) :: idimmin,idimmax

integer :: iigrid, igrid, idims, iside, i1,i2,i3, nxCo1,nxCo2,nxCo3
integer :: ic1,ic2,ic3, inc1,inc2,inc3, ipe_neighbor
integer :: recvsize, sendsize

integer :: recvsize_cc, sendsize_cc
integer :: pi1,pi2,pi3,mi1,mi2,mi3,ph1,ph2,ph3,mh1,mh2,mh3,idir
!-----------------------------------------------------------------------------
nsend=0
nrecv=0
recvsize=0
sendsize=0

! Special communication for diagonal 'coarse corners'
! nrecv/send_cc (for 'coarse corners' is a dim=ndim-1 array which
! stores the faces that must be communicated in each direction.
! nrecv/send_ct (for 'corners total' is the total number of
! necessary communications. These special cases have their own
! send and receive buffers (send/recvbuffer_cc), their tags, etc.
nsend_ct=0
nrecv_ct=0
recvsize_cc=0
sendsize_cc=0

do idims= idimmin,idimmax
   select case (idims)
   case (1)
      nrecv=nrecv+nrecv_fc(1)
      nsend=nsend+nsend_fc(1)
      nxCo1=1;nxCo2=ixGhi2/2-dixB;nxCo3=ixGhi3/2-dixB;
      isize(1)=nxCo1*nxCo2*nxCo3*(nwflux)
      recvsize=recvsize+nrecv_fc(1)*isize(1)
      sendsize=sendsize+nsend_fc(1)*isize(1)

      ! This does not consider the 'coarse corner' case
      nxCo1=1;nxCo2=ixGhi2/2-dixB+1;nxCo3=ixGhi3/2-dixB+1;
      isize_stg(1)=nxCo1*nxCo2*nxCo3*(3-1)
      ! To avoid writing IFDEF STAGGERED everywhere,
      ! the whole size is used (cell centred and staggered)
      isize(1)=isize(1)+isize_stg(1)      
      recvsize=recvsize+nrecv_fc(1)*isize_stg(1)
      sendsize=sendsize+nsend_fc(1)*isize_stg(1)
      ! Coarse corner case
      nrecv_ct=nrecv_ct+nrecv_cc(1)
      nsend_ct=nsend_ct+nsend_cc(1)
      recvsize_cc=recvsize_cc+nrecv_cc(1)*isize_stg(1)
      sendsize_cc=sendsize_cc+nsend_cc(1)*isize_stg(1)


   
   case (2)
      nrecv=nrecv+nrecv_fc(2)
      nsend=nsend+nsend_fc(2)
      nxCo1=ixGhi1/2-dixB;nxCo2=1;nxCo3=ixGhi3/2-dixB;
      isize(2)=nxCo1*nxCo2*nxCo3*(nwflux)
      recvsize=recvsize+nrecv_fc(2)*isize(2)
      sendsize=sendsize+nsend_fc(2)*isize(2)

      ! This does not consider the 'coarse corner' case
      nxCo1=ixGhi1/2-dixB+1;nxCo2=1;nxCo3=ixGhi3/2-dixB+1;
      isize_stg(2)=nxCo1*nxCo2*nxCo3*(3-1)
      ! To avoid writing IFDEF STAGGERED everywhere,
      ! the whole size is used (cell centred and staggered)
      isize(2)=isize(2)+isize_stg(2)      
      recvsize=recvsize+nrecv_fc(2)*isize_stg(2)
      sendsize=sendsize+nsend_fc(2)*isize_stg(2)
      ! Coarse corner case
      nrecv_ct=nrecv_ct+nrecv_cc(2)
      nsend_ct=nsend_ct+nsend_cc(2)
      recvsize_cc=recvsize_cc+nrecv_cc(2)*isize_stg(2)
      sendsize_cc=sendsize_cc+nsend_cc(2)*isize_stg(2)


   
   case (3)
      nrecv=nrecv+nrecv_fc(3)
      nsend=nsend+nsend_fc(3)
      nxCo1=ixGhi1/2-dixB;nxCo2=ixGhi2/2-dixB;nxCo3=1;
      isize(3)=nxCo1*nxCo2*nxCo3*(nwflux)
      recvsize=recvsize+nrecv_fc(3)*isize(3)
      sendsize=sendsize+nsend_fc(3)*isize(3)

      ! This does not consider the 'coarse corner' case
      nxCo1=ixGhi1/2-dixB+1;nxCo2=ixGhi2/2-dixB+1;nxCo3=1;
      isize_stg(3)=nxCo1*nxCo2*nxCo3*(3-1)
      ! To avoid writing IFDEF STAGGERED everywhere,
      ! the whole size is used (cell centred and staggered)
      isize(3)=isize(3)+isize_stg(3)      
      recvsize=recvsize+nrecv_fc(3)*isize_stg(3)
      sendsize=sendsize+nsend_fc(3)*isize_stg(3)
      ! Coarse corner case
      nrecv_ct=nrecv_ct+nrecv_cc(3)
      nsend_ct=nsend_ct+nsend_cc(3)
      recvsize_cc=recvsize_cc+nrecv_cc(3)*isize_stg(3)
      sendsize_cc=sendsize_cc+nsend_cc(3)*isize_stg(3)


   
   end select
end do

if (nrecv>0) then
   ! Receive for direct neighbors
   allocate(recvbuffer(recvsize),recvstatus(MPI_STATUS_SIZE,nrecv),&
       recvrequest(nrecv))
   recvrequest=MPI_REQUEST_NULL

   ibuf=1
   irecv=0

   do iigrid=1,igridstail; igrid=igrids(iigrid);
     call make_task_reflux_direct(igrid,idimmin,idimmax)
   end do

end if


if (nrecv_ct.gt.0) then
   ! Receive corners
   allocate(recvbuffer_cc(recvsize_cc),recvstatus_stg(MPI_STATUS_SIZE,&
      nrecv_ct), &
            recvrequest_stg(nrecv_ct))
   recvrequest_stg=MPI_REQUEST_NULL

   ibuf_cc=1
   irecv_cc=0

   do iigrid=1,igridstail; igrid=igrids(iigrid);
     call make_task_reflux_corners(igrid,idimmin,idimmax)
   end do

end if

if (nsend>0) then
   allocate(sendbuffer(sendsize),sendstatus(MPI_STATUS_SIZE,nsend),&
      sendrequest(nsend))
   sendrequest=MPI_REQUEST_NULL
   isend=0
   ibuf_send=1
end if

if (nsend_ct>0) then
   allocate(sendbuffer_cc(sendsize_cc),sendstatus_stg(MPI_STATUS_SIZE,&
      nsend_ct),sendrequest_stg(nsend_ct))
   sendrequest_stg=MPI_REQUEST_NULL
   isend_cc=0
   ibuf_cc_send=1
end if


end subroutine init_comm_fix_conserve
!=============================================================================
subroutine make_task_reflux_direct(igrid,idimmin,idimmax)
use mod_fix_conserve
use mod_amrvacdef

integer, intent(in) :: igrid,idimmin,idimmax
integer             :: iside,idims,i1,i2,i3,ic1,ic2,ic3,inc1,inc2,inc3
integer             :: ipe_neighbor
!-----------------------------------------------------------------------------

do idims= idimmin,idimmax
   do iside=1,2
      i1=kr(1,idims)*(2*iside-3);i2=kr(2,idims)*(2*iside-3)
      i3=kr(3,idims)*(2*iside-3);
       if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle
      if (neighbor_type(i1,i2,i3,igrid)/=4) cycle
      do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
         inc3=2*i3+ic3
      do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
         inc2=2*i2+ic2
      do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
         inc1=2*i1+ic1
         ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
         if (ipe_neighbor/=mype) then
            irecv=irecv+1
            itag=4**3*(igrid-1)+inc1*4**(1-1)+inc2*4**(2-1)+inc3*4**(3-1)
            call MPI_IRECV(recvbuffer(ibuf),isize(idims), &
               MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
               recvrequest(irecv),ierrmpi)
           !call add_task_to_list(itag,recvrequest(irecv),recvstatus(irecv)) 
            ibuf=ibuf+isize(idims)
         end if
      end do
      end do
      end do
   end do
end do

end subroutine make_task_reflux_direct
!=============================================================================

subroutine make_task_reflux_corners(igrid,idimmin,idimmax)
use mod_fix_conserve
use mod_amrvacdef

integer, intent(in) :: igrid,idimmin,idimmax
integer :: iside,idims,i1,i2,i3,ic1,ic2,ic3,inc1,inc2,inc3
integer :: ipe_neighbor
integer :: pi1,pi2,pi3,mi1,mi2,mi3,ph1,ph2,ph3,mh1,mh2,mh3,idir
!-----------------------------------------------------------------------------

do idims= idimmin,idimmax
   do iside=1,2
     i1=kr(1,idims)*(2*iside-3);i2=kr(2,idims)*(2*iside-3)
     i3=kr(3,idims)*(2*iside-3);
     ! Check if there are special corners
     ! (Coarse block diagonal to a fine block)
     ! If there are, receive.
     ! Tags are calculated in the same way as for
     ! normal fluxes, but should not overlap because
     ! inc1,inc2,inc3 are different
     if (neighbor_type(i1,i2,i3,igrid).eq.3) then
       do idir=idims+1,ndim
         pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
         mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
         ph1=pi1-kr(idims,1)*(2*iside-3);ph2=pi2-kr(idims,2)*(2*iside-3)
         ph3=pi3-kr(idims,3)*(2*iside-3);
         mh1=mi1-kr(idims,1)*(2*iside-3);mh2=mi2-kr(idims,2)*(2*iside-3)
         mh3=mi3-kr(idims,3)*(2*iside-3);

  if (neighbor_pole(pi1,pi2,pi3,igrid).eq.0) then
         if (neighbor_type(pi1,pi2,pi3,igrid).eq.4.and.&
             neighbor_type(ph1,ph2,ph3,igrid).eq.3) then
          ! Loop on children (several in 3D)
          do ic3=1+int((1-pi3)/2),2-int((1+pi3)/2)
              inc3=2*pi3+ic3
do ic2=1+int((1-pi2)/2),2-int((1+pi2)/2)
              inc2=2*pi2+ic2
do ic1=1+int((1-pi1)/2),2-int((1+pi1)/2)
              inc1=2*pi1+ic1
              ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
              if (mype.ne.ipe_neighbor) then
                irecv_cc=irecv_cc+1
                itag_cc=4**3*(igrid-1)+inc1*4**(1-1)+inc2*4**(2-1)+inc3*4&
                   **(3-1)
                call MPI_IRECV(recvbuffer_cc(ibuf_cc),isize_stg(idims),&
                               MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,&
                               icomm,recvrequest_stg(irecv_cc),ierrmpi)
 !call add_task_to_list(itag_cc,recvrequest_stg(irecv_cc),recvstatus_stg(irecv_cc))
                ibuf_cc=ibuf_cc+isize_stg(idims)
              end if
          end do
end do
end do
         end if
  end if

  if (neighbor_pole(mi1,mi2,mi3,igrid).eq.0) then
         if (neighbor_type(mi1,mi2,mi3,igrid).eq.4.and.&
             neighbor_type(mh1,mh2,mh3,igrid).eq.3) then
          ! Loop on children (several in 3D)
          do ic3=1+int((1-mi3)/2),2-int((1+mi3)/2)
              inc3=2*mi3+ic3
do ic2=1+int((1-mi2)/2),2-int((1+mi2)/2)
              inc2=2*mi2+ic2
do ic1=1+int((1-mi1)/2),2-int((1+mi1)/2)
              inc1=2*mi1+ic1
              ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
              if (mype.ne.ipe_neighbor) then
                irecv_cc=irecv_cc+1
                itag_cc=4**3*(igrid-1)+inc1*4**(1-1)+inc2*4**(2-1)+inc3*4&
                   **(3-1)
                call MPI_IRECV(recvbuffer_cc(ibuf_cc),isize_stg(idims),&
                               MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,&
                               icomm,recvrequest_stg(irecv_cc),ierrmpi)
 !call add_task_to_list(itag_cc,recvrequest_stg(irecv_cc),recvstatus_stg(irecv_cc))
                ibuf_cc=ibuf_cc+isize_stg(idims)
              end if
          end do
end do
end do
         end if
  end if
       end do
     end if
   end do
end do


end subroutine make_task_reflux_corners
!=============================================================================

subroutine allocateBflux
use mod_fix_conserve
use mod_amrvacdef

integer :: iigrid, igrid, iside, i1,i2,i3, nx1,nx2,nx3, nxCo1,nxCo2,nxCo3

integer :: idir,idim,pi1,pi2,pi3, mi1,mi2,mi3, ph1,ph2,ph3, mh1,mh2,mh3 !To detect corners

!-----------------------------------------------------------------------------
nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;nx3=ixMhi3-ixMlo3+1;
nxCo1=nx1/2;nxCo2=nx2/2;nxCo3=nx3/2;

do iigrid=1,igridstail; igrid=igrids(iigrid);
   ! For every grid,
   ! arrays for the fluxes are allocated for every face direction(^D)
   ! and every side (1=left, 2=right)
   do iside=1,2
      i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3);i3=kr(3,1)*(2*iside-3);
       if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle
      select case (neighbor_type(i1,i2,i3,igrid))
      case (4)
         allocate(pflux(iside,1,igrid)%flux(1,1:nx2,1:nx3,1:nwflux) ,&
             pflux(iside,1,igrid)%edge(1,0:nx2,0:nx3,1:3-1))
      case (2)
         allocate(pflux(iside,1,igrid)%flux(1,1:nxCo2,1:nxCo3,1:nwflux) ,&
             pflux(iside,1,igrid)%edge(1,0:nxCo2,0:nxCo3,1:3-1))
     
      case(3)
      ! If there are staggered variables, it is necessary
      ! to store also some corners.
      ! Detect the corners that must be stored
      ! and store them
      ! This produces a warning when D=ndim, it might be improved.
      idim=1
      do idir=idim+1,ndim
        pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
        mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
        ph1=pi1-kr(1,1)*(2*iside-3);ph2=pi2-kr(1,2)*(2*iside-3)
        ph3=pi3-kr(1,3)*(2*iside-3);
        mh1=mi1-kr(1,1)*(2*iside-3);mh2=mi2-kr(1,2)*(2*iside-3)
        mh3=mi3-kr(1,3)*(2*iside-3);
        if ((neighbor_type(pi1,pi2,pi3,igrid).eq.4&
            .and.neighbor_type(ph1,ph2,ph3,igrid).eq.3)&
            .or.&
            (neighbor_type(mi1,mi2,mi3,igrid).eq.4&
            .and.neighbor_type(mh1,mh2,mh3,igrid).eq.3)) then
          allocate(pflux(iside,1,igrid)%edge(1,0:nx2,0:nx3,1:3-1))
          exit
        end if
      end do
     
      end select
   end do
   do iside=1,2
      i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3);i3=kr(3,2)*(2*iside-3);
       if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle
      select case (neighbor_type(i1,i2,i3,igrid))
      case (4)
         allocate(pflux(iside,2,igrid)%flux(1:nx1,1,1:nx3,1:nwflux) ,&
             pflux(iside,2,igrid)%edge(0:nx1,1,0:nx3,1:3-1))
      case (2)
         allocate(pflux(iside,2,igrid)%flux(1:nxCo1,1,1:nxCo3,1:nwflux) ,&
             pflux(iside,2,igrid)%edge(0:nxCo1,1,0:nxCo3,1:3-1))
     
      case(3)
      ! If there are staggered variables, it is necessary
      ! to store also some corners.
      ! Detect the corners that must be stored
      ! and store them
      ! This produces a warning when D=ndim, it might be improved.
      idim=2
      do idir=idim+1,ndim
        pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
        mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
        ph1=pi1-kr(2,1)*(2*iside-3);ph2=pi2-kr(2,2)*(2*iside-3)
        ph3=pi3-kr(2,3)*(2*iside-3);
        mh1=mi1-kr(2,1)*(2*iside-3);mh2=mi2-kr(2,2)*(2*iside-3)
        mh3=mi3-kr(2,3)*(2*iside-3);
        if ((neighbor_type(pi1,pi2,pi3,igrid).eq.4&
            .and.neighbor_type(ph1,ph2,ph3,igrid).eq.3)&
            .or.&
            (neighbor_type(mi1,mi2,mi3,igrid).eq.4&
            .and.neighbor_type(mh1,mh2,mh3,igrid).eq.3)) then
          allocate(pflux(iside,2,igrid)%edge(0:nx1,1,0:nx3,1:3-1))
          exit
        end if
      end do
     
      end select
   end do
   do iside=1,2
      i1=kr(1,3)*(2*iside-3);i2=kr(2,3)*(2*iside-3);i3=kr(3,3)*(2*iside-3);
       if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle
      select case (neighbor_type(i1,i2,i3,igrid))
      case (4)
         allocate(pflux(iside,3,igrid)%flux(1:nx1,1:nx2,1,1:nwflux) ,&
             pflux(iside,3,igrid)%edge(0:nx1,0:nx2,1,1:3-1))
      case (2)
         allocate(pflux(iside,3,igrid)%flux(1:nxCo1,1:nxCo2,1,1:nwflux) ,&
             pflux(iside,3,igrid)%edge(0:nxCo1,0:nxCo2,1,1:3-1))
     
      case(3)
      ! If there are staggered variables, it is necessary
      ! to store also some corners.
      ! Detect the corners that must be stored
      ! and store them
      ! This produces a warning when D=ndim, it might be improved.
      idim=3
      do idir=idim+1,ndim
        pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
        mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
        ph1=pi1-kr(3,1)*(2*iside-3);ph2=pi2-kr(3,2)*(2*iside-3)
        ph3=pi3-kr(3,3)*(2*iside-3);
        mh1=mi1-kr(3,1)*(2*iside-3);mh2=mi2-kr(3,2)*(2*iside-3)
        mh3=mi3-kr(3,3)*(2*iside-3);
        if ((neighbor_type(pi1,pi2,pi3,igrid).eq.4&
            .and.neighbor_type(ph1,ph2,ph3,igrid).eq.3)&
            .or.&
            (neighbor_type(mi1,mi2,mi3,igrid).eq.4&
            .and.neighbor_type(mh1,mh2,mh3,igrid).eq.3)) then
          allocate(pflux(iside,3,igrid)%edge(0:nx1,0:nx2,1,1:3-1))
          exit
        end if
      end do
     
      end select
   end do
end do

end subroutine allocateBflux
!=============================================================================
subroutine deallocateBflux
use mod_fix_conserve
use mod_amrvacdef

integer :: iigrid, igrid, iside
!-----------------------------------------------------------------------------
do iigrid=1,igridstail; igrid=igrids(iigrid);
   do iside=1,2
      if (associated(pflux(iside,1,igrid)%flux)) &
         deallocate(pflux(iside,1,igrid)%flux)

      if (associated(pflux(iside,1,igrid)%edge)) &
         deallocate(pflux(iside,1,igrid)%edge)

   end do
   do iside=1,2
      if (associated(pflux(iside,2,igrid)%flux)) &
         deallocate(pflux(iside,2,igrid)%flux)

      if (associated(pflux(iside,2,igrid)%edge)) &
         deallocate(pflux(iside,2,igrid)%edge)

   end do
   do iside=1,2
      if (associated(pflux(iside,3,igrid)%flux)) &
         deallocate(pflux(iside,3,igrid)%flux)

      if (associated(pflux(iside,3,igrid)%edge)) &
         deallocate(pflux(iside,3,igrid)%edge)

   end do
end do

end subroutine deallocateBflux
!=============================================================================
subroutine fix_conserve(psuse,idimmin,idimmax)
use mod_fix_conserve
use mod_amrvacdef

type(state) :: psuse(ngridshi)
integer, intent(in) :: idimmin,idimmax

integer :: iigrid, igrid, idims, iside, iotherside, i1,i2,i3, ic1,ic2,ic3,&
    inc1,inc2,inc3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
integer :: nxCo1,nxCo2,nxCo3, iw, ix, ipe_neighbor, ineighbor, ibufnext, nbuf
double precision :: CoFiratio
!-----------------------------------------------------------------------------
if (slab) then
   ! The flux is divided by volume of fine cell. We need, however,
   ! to divide by volume of coarse cell => muliply by volume ratio
   CoFiratio=one/dble(2**ndim)
end if

if (nrecv>0) then
   call MPI_WAITALL(nrecv,recvrequest,recvstatus,ierrmpi)
   ibuf=1
end if
nxCo1=(ixMhi1-ixMlo1+1)/2;nxCo2=(ixMhi2-ixMlo2+1)/2;nxCo3=(ixMhi3-ixMlo3+1)/2;

!if (.false.) then 

! for all grids: perform flux update at Coarse-Fine interfaces
do iigrid=1,igridstail; igrid=igrids(iigrid);
   do idims= idimmin,idimmax
      select case (idims)
      case (1)
         do iside=1,2
            i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3)
            i3=kr(3,1)*(2*iside-3);
             if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle
            if (neighbor_type(i1,i2,i3,igrid)/=4) cycle

! opedit: skip over active/passive interface since flux for passive ones is 
! not computed, keep the buffer counter up to date:
            if (.not.neighbor_active(i1,i2,i3,igrid).or.&
                .not.neighbor_active(0,0,0,igrid) ) then
               do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
               inc3=2*i3+ic3
      do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
               inc2=2*i2+ic2
      do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
               inc1=2*i1+ic1
               ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
               if (ipe_neighbor/=mype) then
                  ibufnext=ibuf+isize(1)
                  ibuf=ibufnext
                  end if
               end do
      end do
      end do
               cycle
            end if
!

            select case (iside)
            case (1)
               ix=ixMlo1
            case (2)
               ix=ixMhi1
            end select

            ! remove coarse flux
            if (slab) then
               psuse(igrid)%w%w(ix,ixMlo2:ixMhi2,ixMlo3:ixMhi3,1:nwflux) &
                  = psuse(igrid)%w%w(ix,ixMlo2:ixMhi2,ixMlo3:ixMhi3,&
                     1:nwflux) &
                   -pflux(iside,1,igrid)%flux(1,:,:,1:nwflux)
            else
               do iw=1,nwflux
                  psuse(igrid)%w%w(ix,ixMlo2:ixMhi2,ixMlo3:ixMhi3,iw)&
                     =psuse(igrid)%w%w(ix,ixMlo2:ixMhi2,ixMlo3:ixMhi3,iw)&
                     -pflux(iside,1,igrid)%flux(1,:,:,iw) &
                     /pgeo(igrid)%dvolume(ix,ixMlo2:ixMhi2,ixMlo3:ixMhi3)
               end do
            end if


            ! add fine flux
            do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
               inc3=2*i3+ic3
      do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
               inc2=2*i2+ic2
      do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
               inc1=2*i1+ic1
               ineighbor=neighbor_child(1,inc1,inc2,inc3,igrid)
               ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
               ixmin1=ix;ixmin2=ixMlo2+(ic2-1)*nxCo2
               ixmin3=ixMlo3+(ic3-1)*nxCo3;
               ixmax1=ix;ixmax2=ixmin2-1+nxCo2;ixmax3=ixmin3-1+nxCo3;
               if (ipe_neighbor==mype) then
                  iotherside=3-iside
                  if (slab) then
                     psuse(igrid)%w%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                        ixmin3:ixmax3,1:nwflux) &
                       = psuse(igrid)%w%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                          ixmin3:ixmax3,1:nwflux) &
                       + pflux(iotherside,1,ineighbor)%flux(:,:,:,1:nwflux)&
                       * CoFiratio
                  else
                     do iw=1,nwflux
                        psuse(igrid)%w%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                           ixmin3:ixmax3,iw)=psuse(igrid)%w%w(ixmin1:ixmax1,&
                           ixmin2:ixmax2,ixmin3:ixmax3,iw) &
                            +pflux(iotherside,1,ineighbor)%flux(:,:,:,iw) &
                            /pgeo(igrid)%dvolume(ixmin1:ixmax1,ixmin2:ixmax2,&
                               ixmin3:ixmax3)
                     end do
                  end if
               else
                  if (slab) then
                     ibufnext=ibuf+isize(1)
                     
                     
                     psuse(igrid)%w%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                        ixmin3:ixmax3,1:nwflux) &
                         = psuse(igrid)%w%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                            ixmin3:ixmax3,1:nwflux)+CoFiratio*reshape(source&
                            =recvbuffer(ibuf:ibufnext-isize_stg(1)-1), &
                                 shape=shape(psuse(igrid)%w%w(ixmin1:ixmax1,&
                                    ixmin2:ixmax2,ixmin3:ixmax3,1:nwflux)))
                     
                     ibuf=ibufnext
                  else
                     ibufnext=ibuf+isize(1)
                     
                     
                     nbuf=(isize(1)-isize_stg(1))/nwflux
                     
                     do iw=1,nwflux
                        psuse(igrid)%w%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                           ixmin3:ixmax3,iw)=psuse(igrid)%w%w(ixmin1:ixmax1,&
                           ixmin2:ixmax2,ixmin3:ixmax3,iw) &
                           +reshape(source=recvbuffer(ibuf:ibuf+nbuf-1), &
                                    shape=shape(psuse(igrid)%w%w&
                                       (ixmin1:ixmax1,ixmin2:ixmax2,&
                                       ixmin3:ixmax3,iw))) &
                           /pgeo(igrid)%dvolume(ixmin1:ixmax1,ixmin2:ixmax2,&
                              ixmin3:ixmax3)
                        ibuf=ibuf+nbuf
                     end do
                     ibuf=ibufnext
                  end if
               end if
            end do
      end do
      end do
         end do
      case (2)
         do iside=1,2
            i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3)
            i3=kr(3,2)*(2*iside-3);
             if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle
            if (neighbor_type(i1,i2,i3,igrid)/=4) cycle

! opedit: skip over active/passive interface since flux for passive ones is 
! not computed, keep the buffer counter up to date:
            if (.not.neighbor_active(i1,i2,i3,igrid).or.&
                .not.neighbor_active(0,0,0,igrid) ) then
               do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
               inc3=2*i3+ic3
      do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
               inc2=2*i2+ic2
      do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
               inc1=2*i1+ic1
               ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
               if (ipe_neighbor/=mype) then
                  ibufnext=ibuf+isize(2)
                  ibuf=ibufnext
                  end if
               end do
      end do
      end do
               cycle
            end if
!

            select case (iside)
            case (1)
               ix=ixMlo2
            case (2)
               ix=ixMhi2
            end select

            ! remove coarse flux
            if (slab) then
               psuse(igrid)%w%w(ixMlo1:ixMhi1,ix,ixMlo3:ixMhi3,1:nwflux) &
                  = psuse(igrid)%w%w(ixMlo1:ixMhi1,ix,ixMlo3:ixMhi3,&
                     1:nwflux) &
                   -pflux(iside,2,igrid)%flux(:,1,:,1:nwflux)
            else
               do iw=1,nwflux
                  psuse(igrid)%w%w(ixMlo1:ixMhi1,ix,ixMlo3:ixMhi3,iw)&
                     =psuse(igrid)%w%w(ixMlo1:ixMhi1,ix,ixMlo3:ixMhi3,iw)&
                     -pflux(iside,2,igrid)%flux(:,1,:,iw) &
                     /pgeo(igrid)%dvolume(ixMlo1:ixMhi1,ix,ixMlo3:ixMhi3)
               end do
            end if


            ! add fine flux
            do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
               inc3=2*i3+ic3
      do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
               inc2=2*i2+ic2
      do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
               inc1=2*i1+ic1
               ineighbor=neighbor_child(1,inc1,inc2,inc3,igrid)
               ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
               ixmin1=ixMlo1+(ic1-1)*nxCo1;ixmin2=ix
               ixmin3=ixMlo3+(ic3-1)*nxCo3;
               ixmax1=ixmin1-1+nxCo1;ixmax2=ix;ixmax3=ixmin3-1+nxCo3;
               if (ipe_neighbor==mype) then
                  iotherside=3-iside
                  if (slab) then
                     psuse(igrid)%w%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                        ixmin3:ixmax3,1:nwflux) &
                       = psuse(igrid)%w%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                          ixmin3:ixmax3,1:nwflux) &
                       + pflux(iotherside,2,ineighbor)%flux(:,:,:,1:nwflux)&
                       * CoFiratio
                  else
                     do iw=1,nwflux
                        psuse(igrid)%w%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                           ixmin3:ixmax3,iw)=psuse(igrid)%w%w(ixmin1:ixmax1,&
                           ixmin2:ixmax2,ixmin3:ixmax3,iw) &
                            +pflux(iotherside,2,ineighbor)%flux(:,:,:,iw) &
                            /pgeo(igrid)%dvolume(ixmin1:ixmax1,ixmin2:ixmax2,&
                               ixmin3:ixmax3)
                     end do
                  end if
               else
                  if (slab) then
                     ibufnext=ibuf+isize(2)
                     
                     
                     psuse(igrid)%w%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                        ixmin3:ixmax3,1:nwflux) &
                         = psuse(igrid)%w%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                            ixmin3:ixmax3,1:nwflux)+CoFiratio*reshape(source&
                            =recvbuffer(ibuf:ibufnext-isize_stg(2)-1), &
                                 shape=shape(psuse(igrid)%w%w(ixmin1:ixmax1,&
                                    ixmin2:ixmax2,ixmin3:ixmax3,1:nwflux)))
                     
                     ibuf=ibufnext
                  else
                     ibufnext=ibuf+isize(2)
                     
                     
                     nbuf=(isize(2)-isize_stg(2))/nwflux
                     
                     do iw=1,nwflux
                        psuse(igrid)%w%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                           ixmin3:ixmax3,iw)=psuse(igrid)%w%w(ixmin1:ixmax1,&
                           ixmin2:ixmax2,ixmin3:ixmax3,iw) &
                           +reshape(source=recvbuffer(ibuf:ibuf+nbuf-1), &
                                    shape=shape(psuse(igrid)%w%w&
                                       (ixmin1:ixmax1,ixmin2:ixmax2,&
                                       ixmin3:ixmax3,iw))) &
                           /pgeo(igrid)%dvolume(ixmin1:ixmax1,ixmin2:ixmax2,&
                              ixmin3:ixmax3)
                        ibuf=ibuf+nbuf
                     end do
                     ibuf=ibufnext
                  end if
               end if
            end do
      end do
      end do
         end do
      case (3)
         do iside=1,2
            i1=kr(1,3)*(2*iside-3);i2=kr(2,3)*(2*iside-3)
            i3=kr(3,3)*(2*iside-3);
             if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle
            if (neighbor_type(i1,i2,i3,igrid)/=4) cycle

! opedit: skip over active/passive interface since flux for passive ones is 
! not computed, keep the buffer counter up to date:
            if (.not.neighbor_active(i1,i2,i3,igrid).or.&
                .not.neighbor_active(0,0,0,igrid) ) then
               do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
               inc3=2*i3+ic3
      do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
               inc2=2*i2+ic2
      do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
               inc1=2*i1+ic1
               ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
               if (ipe_neighbor/=mype) then
                  ibufnext=ibuf+isize(3)
                  ibuf=ibufnext
                  end if
               end do
      end do
      end do
               cycle
            end if
!

            select case (iside)
            case (1)
               ix=ixMlo3
            case (2)
               ix=ixMhi3
            end select

            ! remove coarse flux
            if (slab) then
               psuse(igrid)%w%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ix,1:nwflux) &
                  = psuse(igrid)%w%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ix,&
                     1:nwflux) &
                   -pflux(iside,3,igrid)%flux(:,:,1,1:nwflux)
            else
               do iw=1,nwflux
                  psuse(igrid)%w%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ix,iw)&
                     =psuse(igrid)%w%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ix,iw)&
                     -pflux(iside,3,igrid)%flux(:,:,1,iw) &
                     /pgeo(igrid)%dvolume(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ix)
               end do
            end if


            ! add fine flux
            do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
               inc3=2*i3+ic3
      do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
               inc2=2*i2+ic2
      do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
               inc1=2*i1+ic1
               ineighbor=neighbor_child(1,inc1,inc2,inc3,igrid)
               ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
               ixmin1=ixMlo1+(ic1-1)*nxCo1;ixmin2=ixMlo2+(ic2-1)*nxCo2
               ixmin3=ix;
               ixmax1=ixmin1-1+nxCo1;ixmax2=ixmin2-1+nxCo2;ixmax3=ix;
               if (ipe_neighbor==mype) then
                  iotherside=3-iside
                  if (slab) then
                     psuse(igrid)%w%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                        ixmin3:ixmax3,1:nwflux) &
                       = psuse(igrid)%w%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                          ixmin3:ixmax3,1:nwflux) &
                       + pflux(iotherside,3,ineighbor)%flux(:,:,:,1:nwflux)&
                       * CoFiratio
                  else
                     do iw=1,nwflux
                        psuse(igrid)%w%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                           ixmin3:ixmax3,iw)=psuse(igrid)%w%w(ixmin1:ixmax1,&
                           ixmin2:ixmax2,ixmin3:ixmax3,iw) &
                            +pflux(iotherside,3,ineighbor)%flux(:,:,:,iw) &
                            /pgeo(igrid)%dvolume(ixmin1:ixmax1,ixmin2:ixmax2,&
                               ixmin3:ixmax3)
                     end do
                  end if
               else
                  if (slab) then
                     ibufnext=ibuf+isize(3)
                     
                     
                     psuse(igrid)%w%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                        ixmin3:ixmax3,1:nwflux) &
                         = psuse(igrid)%w%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                            ixmin3:ixmax3,1:nwflux)+CoFiratio*reshape(source&
                            =recvbuffer(ibuf:ibufnext-isize_stg(3)-1), &
                                 shape=shape(psuse(igrid)%w%w(ixmin1:ixmax1,&
                                    ixmin2:ixmax2,ixmin3:ixmax3,1:nwflux)))
                     
                     ibuf=ibufnext
                  else
                     ibufnext=ibuf+isize(3)
                     
                     
                     nbuf=(isize(3)-isize_stg(3))/nwflux
                     
                     do iw=1,nwflux
                        psuse(igrid)%w%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                           ixmin3:ixmax3,iw)=psuse(igrid)%w%w(ixmin1:ixmax1,&
                           ixmin2:ixmax2,ixmin3:ixmax3,iw) &
                           +reshape(source=recvbuffer(ibuf:ibuf+nbuf-1), &
                                    shape=shape(psuse(igrid)%w%w&
                                       (ixmin1:ixmax1,ixmin2:ixmax2,&
                                       ixmin3:ixmax3,iw))) &
                           /pgeo(igrid)%dvolume(ixmin1:ixmax1,ixmin2:ixmax2,&
                              ixmin3:ixmax3)
                        ibuf=ibuf+nbuf
                     end do
                     ibuf=ibufnext
                  end if
               end if
            end do
      end do
      end do
         end do
      end select
   end do
end do

!end if



end subroutine fix_conserve
!=============================================================================
subroutine storeflux(igrid,fC,idimmin,idimmax)
use mod_fix_conserve
use mod_amrvacdef

integer, intent(in)          :: igrid, idimmin,idimmax
double precision, intent(in) :: fC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
   1:nwflux,1:ndim)

integer :: idims, iside, i1,i2,i3, ic1,ic2,ic3, inc1,inc2,inc3, ix1,ix2,ix3,&
    ixCo1,ixCo2,ixCo3, nxCo1,nxCo2,nxCo3, iw
!-----------------------------------------------------------------------------
do idims = idimmin,idimmax
   select case (idims)
   case (1)
      do iside=1,2
         i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3);i3=kr(3,1)*(2*iside-3);
          if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle
         select case (neighbor_type(i1,i2,i3,igrid))
         case (4)
            select case (iside)
            case (1)
               pflux(iside,1,igrid)%flux(1,:,:,1:nwflux) = &
                  -fC(dixB,ixMlo2:ixMhi2,ixMlo3:ixMhi3,1:nwflux,1)
            case (2)
               pflux(iside,1,igrid)%flux(1,:,:,1:nwflux) = &
                  fC(ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,1:nwflux,1)
            end select
         case (2)
            nxCo1=1;nxCo2=ixGhi2/2-dixB;nxCo3=ixGhi3/2-dixB;
            select case (iside)
            case (1)
               do iw=1,nwflux
                  do ixCo3=1,nxCo3
   do ixCo2=1,nxCo2
   do ixCo1=1,nxCo1
                     ix1=dixB;ix2=ixMlo2+2*(ixCo2-1);ix3=ixMlo3+2*(ixCo3-1);
                     pflux(iside,1,igrid)%flux(ixCo1,ixCo2,ixCo3,iw) &
                        = sum(fC(ix1,ix2:ix2+1,ix3:ix3+1,iw,1))
                  end do
   end do
   end do
               end do
            case (2)
               do iw=1,nwflux
                  do ixCo3=1,nxCo3
   do ixCo2=1,nxCo2
   do ixCo1=1,nxCo1
                     ix1=ixMhi1;ix2=ixMlo2+2*(ixCo2-1);ix3=ixMlo3+2*(ixCo3-1);
                     pflux(iside,1,igrid)%flux(ixCo1,ixCo2,ixCo3,iw) &
                        =-sum(fC(ix1,ix2:ix2+1,ix3:ix3+1,iw,1))
                  end do
   end do
   end do
               end do
            end select
         end select
      end do
   case (2)
      do iside=1,2
         i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3);i3=kr(3,2)*(2*iside-3);
          if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle
         select case (neighbor_type(i1,i2,i3,igrid))
         case (4)
            select case (iside)
            case (1)
               pflux(iside,2,igrid)%flux(:,1,:,1:nwflux) = &
                  -fC(ixMlo1:ixMhi1,dixB,ixMlo3:ixMhi3,1:nwflux,2)
            case (2)
               pflux(iside,2,igrid)%flux(:,1,:,1:nwflux) = &
                  fC(ixMlo1:ixMhi1,ixMhi2,ixMlo3:ixMhi3,1:nwflux,2)
            end select
         case (2)
            nxCo1=ixGhi1/2-dixB;nxCo2=1;nxCo3=ixGhi3/2-dixB;
            select case (iside)
            case (1)
               do iw=1,nwflux
                  do ixCo3=1,nxCo3
   do ixCo2=1,nxCo2
   do ixCo1=1,nxCo1
                     ix1=ixMlo1+2*(ixCo1-1);ix2=dixB;ix3=ixMlo3+2*(ixCo3-1);
                     pflux(iside,2,igrid)%flux(ixCo1,ixCo2,ixCo3,iw) &
                        = sum(fC(ix1:ix1+1,ix2,ix3:ix3+1,iw,2))
                  end do
   end do
   end do
               end do
            case (2)
               do iw=1,nwflux
                  do ixCo3=1,nxCo3
   do ixCo2=1,nxCo2
   do ixCo1=1,nxCo1
                     ix1=ixMlo1+2*(ixCo1-1);ix2=ixMhi2;ix3=ixMlo3+2*(ixCo3-1);
                     pflux(iside,2,igrid)%flux(ixCo1,ixCo2,ixCo3,iw) &
                        =-sum(fC(ix1:ix1+1,ix2,ix3:ix3+1,iw,2))
                  end do
   end do
   end do
               end do
            end select
         end select
      end do
   case (3)
      do iside=1,2
         i1=kr(1,3)*(2*iside-3);i2=kr(2,3)*(2*iside-3);i3=kr(3,3)*(2*iside-3);
          if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle
         select case (neighbor_type(i1,i2,i3,igrid))
         case (4)
            select case (iside)
            case (1)
               pflux(iside,3,igrid)%flux(:,:,1,1:nwflux) = &
                  -fC(ixMlo1:ixMhi1,ixMlo2:ixMhi2,dixB,1:nwflux,3)
            case (2)
               pflux(iside,3,igrid)%flux(:,:,1,1:nwflux) = &
                  fC(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMhi3,1:nwflux,3)
            end select
         case (2)
            nxCo1=ixGhi1/2-dixB;nxCo2=ixGhi2/2-dixB;nxCo3=1;
            select case (iside)
            case (1)
               do iw=1,nwflux
                  do ixCo3=1,nxCo3
   do ixCo2=1,nxCo2
   do ixCo1=1,nxCo1
                     ix1=ixMlo1+2*(ixCo1-1);ix2=ixMlo2+2*(ixCo2-1);ix3=dixB;
                     pflux(iside,3,igrid)%flux(ixCo1,ixCo2,ixCo3,iw) &
                        = sum(fC(ix1:ix1+1,ix2:ix2+1,ix3,iw,3))
                  end do
   end do
   end do
               end do
            case (2)
               do iw=1,nwflux
                  do ixCo3=1,nxCo3
   do ixCo2=1,nxCo2
   do ixCo1=1,nxCo1
                     ix1=ixMlo1+2*(ixCo1-1);ix2=ixMlo2+2*(ixCo2-1);ix3=ixMhi3;
                     pflux(iside,3,igrid)%flux(ixCo1,ixCo2,ixCo3,iw) &
                        =-sum(fC(ix1:ix1+1,ix2:ix2+1,ix3,iw,3))
                  end do
   end do
   end do
               end do
            end select
         end select
      end do
   end select
end do

end subroutine storeflux
!=============================================================================
subroutine sendflux(igrid,idimmin,idimmax)
use mod_fix_conserve
use mod_amrvacdef

integer, intent(in) :: idimmin,idimmax
integer             :: igrid, iigrid
integer             :: ibuf_send_next
integer :: idims, iside, i1,i2,i3, ic1,ic2,ic3, inc1,inc2,inc3, ix1,ix2,ix3,&
    ixCo1,ixCo2,ixCo3, nxCo1,nxCo2,nxCo3, iw
integer :: ineighbor, ipe_neighbor

integer             :: idir,ibuf_cc_send_next,pi1,pi2,pi3,ph1,ph2,ph3,mi1,mi2,&
   mi3,mh1,mh2,mh3

!----------------------------------------------------------------------------

do idims = idimmin,idimmax
   select case (idims)
   case (1)
      do iside=1,2
         i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3);i3=kr(3,1)*(2*iside-3);
          if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle
         
         if (neighbor_type(i1,i2,i3,igrid)==2) then

            ineighbor=neighbor(1,i1,i2,i3,igrid)
            ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
            if (ipe_neighbor/=mype) then
               ic1=1+modulo(node(pig1_,igrid)-1,2)
               ic2=1+modulo(node(pig2_,igrid)-1,2)
               ic3=1+modulo(node(pig3_,igrid)-1,2);
               inc1=-2*i1+ic1;inc2=ic2;inc3=ic3;
               itag=4**3*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)+inc3*4&
                  **(3-1)
               isend=isend+1

               ibuf_send_next=ibuf_send+isize(1)
               

               
               sendbuffer(ibuf_send:ibuf_send_next-isize_stg(1)-1)=&
               reshape(pflux(iside,1,igrid)%flux,(/isize(1)-isize_stg(1)/))

               sendbuffer(ibuf_send_next-isize_stg(1):ibuf_send_next-1)=&
               reshape(pflux(iside,1,igrid)%edge,(/isize_stg(1)/))

               call MPI_ISEND(sendbuffer(ibuf_send),isize(1), &
                              MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                              icomm,sendrequest(isend),ierrmpi)

               ibuf_send=ibuf_send_next

            end if
            
            ! If we are in a fine block surrounded by coarse blocks
            do idir=idims+1,ndim
              pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
              mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
              ph1=pi1-kr(idims,1)*(2*iside-3);ph2=pi2-kr(idims,2)*(2*iside-3)
              ph3=pi3-kr(idims,3)*(2*iside-3);
              mh1=mi1-kr(idims,1)*(2*iside-3);mh2=mi2-kr(idims,2)*(2*iside-3)
              mh3=mi3-kr(idims,3)*(2*iside-3);

              if (neighbor_type(pi1,pi2,pi3,igrid).eq.2.and.&
                  neighbor_type(ph1,ph2,ph3,igrid).eq.2.and.&
                  mype.ne.neighbor(2,pi1,pi2,pi3,igrid)) then
         if (neighbor_pole(pi1,pi2,pi3,igrid).eq.0) then
                ! Get relative position in the grid for tags
                ineighbor=neighbor(1,pi1,pi2,pi3,igrid)
                ipe_neighbor=neighbor(2,pi1,pi2,pi3,igrid)
                ic1=1+modulo(node(pig1_,igrid)-1,2)
                ic2=1+modulo(node(pig2_,igrid)-1,2)
                ic3=1+modulo(node(pig3_,igrid)-1,2);
                inc1=-2*pi1+ic1;inc2=-2*pi2+ic2;inc3=-2*pi3+ic3;
                itag_cc=4**3*(ineighbor-1)+inc1*4**(1-1)+inc2*4&
                   **(2-1)+inc3*4**(3-1)
                ! Reshape to buffer and send
                isend_cc=isend_cc+1
                ibuf_cc_send_next=ibuf_cc_send+isize_stg(1)
                sendbuffer_cc(ibuf_cc_send:ibuf_cc_send_next-1)=&
                reshape(pflux(iside,1,igrid)%edge,shape=(/isize_stg(1)/))
                call MPI_ISEND(sendbuffer_cc(ibuf_cc_send),isize_stg(1),&
                               MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,&
                               icomm,sendrequest_stg(isend_cc),ierrmpi)
                ibuf_cc_send=ibuf_cc_send_next
         end if
              end if
 
              if (neighbor_type(mi1,mi2,mi3,igrid).eq.2.and.&
                  neighbor_type(mh1,mh2,mh3,igrid).eq.2.and.&
                  mype.ne.neighbor(2,mi1,mi2,mi3,igrid)) then
         if (neighbor_pole(mi1,mi2,mi3,igrid).eq.0) then
                ! Get relative position in the grid for tags
                ineighbor=neighbor(1,mi1,mi2,mi3,igrid)
                ipe_neighbor=neighbor(2,mi1,mi2,mi3,igrid)
                ic1=1+modulo(node(pig1_,igrid)-1,2)
                ic2=1+modulo(node(pig2_,igrid)-1,2)
                ic3=1+modulo(node(pig3_,igrid)-1,2);
                inc1=-2*pi1+ic1;inc2=-2*pi2+ic2;inc3=-2*pi3+ic3;
                inc1=-2*mi1+ic1;inc2=-2*mi2+ic2;inc3=-2*mi3+ic3;
                itag_cc=4**3*(ineighbor-1)+inc1*4**(1-1)+inc2*4&
                   **(2-1)+inc3*4**(3-1)
                ! Reshape to buffer and send
                isend_cc=isend_cc+1
                ibuf_cc_send_next=ibuf_cc_send+isize_stg(1)
                sendbuffer_cc(ibuf_cc_send:ibuf_cc_send_next-1)=&
                
                reshape(pflux(iside,1,igrid)%edge,shape=(/isize_stg(1)/))
                call MPI_ISEND(sendbuffer_cc(ibuf_cc_send),isize_stg(1),&
                               MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,&
                               icomm,sendrequest_stg(isend_cc),ierrmpi)
                ibuf_cc_send=ibuf_cc_send_next
         end if
              end if
            end do
            
         end if
      end do
   case (2)
      do iside=1,2
         i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3);i3=kr(3,2)*(2*iside-3);
          if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle
         
         if (neighbor_type(i1,i2,i3,igrid)==2) then

            ineighbor=neighbor(1,i1,i2,i3,igrid)
            ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
            if (ipe_neighbor/=mype) then
               ic1=1+modulo(node(pig1_,igrid)-1,2)
               ic2=1+modulo(node(pig2_,igrid)-1,2)
               ic3=1+modulo(node(pig3_,igrid)-1,2);
               inc1=ic1;inc2=-2*i2+ic2;inc3=ic3;
               itag=4**3*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)+inc3*4&
                  **(3-1)
               isend=isend+1

               ibuf_send_next=ibuf_send+isize(2)
               

               
               sendbuffer(ibuf_send:ibuf_send_next-isize_stg(2)-1)=&
               reshape(pflux(iside,2,igrid)%flux,(/isize(2)-isize_stg(2)/))

               sendbuffer(ibuf_send_next-isize_stg(2):ibuf_send_next-1)=&
               reshape(pflux(iside,2,igrid)%edge,(/isize_stg(2)/))

               call MPI_ISEND(sendbuffer(ibuf_send),isize(2), &
                              MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                              icomm,sendrequest(isend),ierrmpi)

               ibuf_send=ibuf_send_next

            end if
            
            ! If we are in a fine block surrounded by coarse blocks
            do idir=idims+1,ndim
              pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
              mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
              ph1=pi1-kr(idims,1)*(2*iside-3);ph2=pi2-kr(idims,2)*(2*iside-3)
              ph3=pi3-kr(idims,3)*(2*iside-3);
              mh1=mi1-kr(idims,1)*(2*iside-3);mh2=mi2-kr(idims,2)*(2*iside-3)
              mh3=mi3-kr(idims,3)*(2*iside-3);

              if (neighbor_type(pi1,pi2,pi3,igrid).eq.2.and.&
                  neighbor_type(ph1,ph2,ph3,igrid).eq.2.and.&
                  mype.ne.neighbor(2,pi1,pi2,pi3,igrid)) then
         if (neighbor_pole(pi1,pi2,pi3,igrid).eq.0) then
                ! Get relative position in the grid for tags
                ineighbor=neighbor(1,pi1,pi2,pi3,igrid)
                ipe_neighbor=neighbor(2,pi1,pi2,pi3,igrid)
                ic1=1+modulo(node(pig1_,igrid)-1,2)
                ic2=1+modulo(node(pig2_,igrid)-1,2)
                ic3=1+modulo(node(pig3_,igrid)-1,2);
                inc1=-2*pi1+ic1;inc2=-2*pi2+ic2;inc3=-2*pi3+ic3;
                itag_cc=4**3*(ineighbor-1)+inc1*4**(1-1)+inc2*4&
                   **(2-1)+inc3*4**(3-1)
                ! Reshape to buffer and send
                isend_cc=isend_cc+1
                ibuf_cc_send_next=ibuf_cc_send+isize_stg(2)
                sendbuffer_cc(ibuf_cc_send:ibuf_cc_send_next-1)=&
                reshape(pflux(iside,2,igrid)%edge,shape=(/isize_stg(2)/))
                call MPI_ISEND(sendbuffer_cc(ibuf_cc_send),isize_stg(2),&
                               MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,&
                               icomm,sendrequest_stg(isend_cc),ierrmpi)
                ibuf_cc_send=ibuf_cc_send_next
         end if
              end if
 
              if (neighbor_type(mi1,mi2,mi3,igrid).eq.2.and.&
                  neighbor_type(mh1,mh2,mh3,igrid).eq.2.and.&
                  mype.ne.neighbor(2,mi1,mi2,mi3,igrid)) then
         if (neighbor_pole(mi1,mi2,mi3,igrid).eq.0) then
                ! Get relative position in the grid for tags
                ineighbor=neighbor(1,mi1,mi2,mi3,igrid)
                ipe_neighbor=neighbor(2,mi1,mi2,mi3,igrid)
                ic1=1+modulo(node(pig1_,igrid)-1,2)
                ic2=1+modulo(node(pig2_,igrid)-1,2)
                ic3=1+modulo(node(pig3_,igrid)-1,2);
                inc1=-2*pi1+ic1;inc2=-2*pi2+ic2;inc3=-2*pi3+ic3;
                inc1=-2*mi1+ic1;inc2=-2*mi2+ic2;inc3=-2*mi3+ic3;
                itag_cc=4**3*(ineighbor-1)+inc1*4**(1-1)+inc2*4&
                   **(2-1)+inc3*4**(3-1)
                ! Reshape to buffer and send
                isend_cc=isend_cc+1
                ibuf_cc_send_next=ibuf_cc_send+isize_stg(2)
                sendbuffer_cc(ibuf_cc_send:ibuf_cc_send_next-1)=&
                
                reshape(pflux(iside,2,igrid)%edge,shape=(/isize_stg(2)/))
                call MPI_ISEND(sendbuffer_cc(ibuf_cc_send),isize_stg(2),&
                               MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,&
                               icomm,sendrequest_stg(isend_cc),ierrmpi)
                ibuf_cc_send=ibuf_cc_send_next
         end if
              end if
            end do
            
         end if
      end do
   case (3)
      do iside=1,2
         i1=kr(1,3)*(2*iside-3);i2=kr(2,3)*(2*iside-3);i3=kr(3,3)*(2*iside-3);
          if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle
         
         if (neighbor_type(i1,i2,i3,igrid)==2) then

            ineighbor=neighbor(1,i1,i2,i3,igrid)
            ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
            if (ipe_neighbor/=mype) then
               ic1=1+modulo(node(pig1_,igrid)-1,2)
               ic2=1+modulo(node(pig2_,igrid)-1,2)
               ic3=1+modulo(node(pig3_,igrid)-1,2);
               inc1=ic1;inc2=ic2;inc3=-2*i3+ic3;
               itag=4**3*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)+inc3*4&
                  **(3-1)
               isend=isend+1

               ibuf_send_next=ibuf_send+isize(3)
               

               
               sendbuffer(ibuf_send:ibuf_send_next-isize_stg(3)-1)=&
               reshape(pflux(iside,3,igrid)%flux,(/isize(3)-isize_stg(3)/))

               sendbuffer(ibuf_send_next-isize_stg(3):ibuf_send_next-1)=&
               reshape(pflux(iside,3,igrid)%edge,(/isize_stg(3)/))

               call MPI_ISEND(sendbuffer(ibuf_send),isize(3), &
                              MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                              icomm,sendrequest(isend),ierrmpi)

               ibuf_send=ibuf_send_next

            end if
            
            ! If we are in a fine block surrounded by coarse blocks
            do idir=idims+1,ndim
              pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
              mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
              ph1=pi1-kr(idims,1)*(2*iside-3);ph2=pi2-kr(idims,2)*(2*iside-3)
              ph3=pi3-kr(idims,3)*(2*iside-3);
              mh1=mi1-kr(idims,1)*(2*iside-3);mh2=mi2-kr(idims,2)*(2*iside-3)
              mh3=mi3-kr(idims,3)*(2*iside-3);

              if (neighbor_type(pi1,pi2,pi3,igrid).eq.2.and.&
                  neighbor_type(ph1,ph2,ph3,igrid).eq.2.and.&
                  mype.ne.neighbor(2,pi1,pi2,pi3,igrid)) then
         if (neighbor_pole(pi1,pi2,pi3,igrid).eq.0) then
                ! Get relative position in the grid for tags
                ineighbor=neighbor(1,pi1,pi2,pi3,igrid)
                ipe_neighbor=neighbor(2,pi1,pi2,pi3,igrid)
                ic1=1+modulo(node(pig1_,igrid)-1,2)
                ic2=1+modulo(node(pig2_,igrid)-1,2)
                ic3=1+modulo(node(pig3_,igrid)-1,2);
                inc1=-2*pi1+ic1;inc2=-2*pi2+ic2;inc3=-2*pi3+ic3;
                itag_cc=4**3*(ineighbor-1)+inc1*4**(1-1)+inc2*4&
                   **(2-1)+inc3*4**(3-1)
                ! Reshape to buffer and send
                isend_cc=isend_cc+1
                ibuf_cc_send_next=ibuf_cc_send+isize_stg(3)
                sendbuffer_cc(ibuf_cc_send:ibuf_cc_send_next-1)=&
                reshape(pflux(iside,3,igrid)%edge,shape=(/isize_stg(3)/))
                call MPI_ISEND(sendbuffer_cc(ibuf_cc_send),isize_stg(3),&
                               MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,&
                               icomm,sendrequest_stg(isend_cc),ierrmpi)
                ibuf_cc_send=ibuf_cc_send_next
         end if
              end if
 
              if (neighbor_type(mi1,mi2,mi3,igrid).eq.2.and.&
                  neighbor_type(mh1,mh2,mh3,igrid).eq.2.and.&
                  mype.ne.neighbor(2,mi1,mi2,mi3,igrid)) then
         if (neighbor_pole(mi1,mi2,mi3,igrid).eq.0) then
                ! Get relative position in the grid for tags
                ineighbor=neighbor(1,mi1,mi2,mi3,igrid)
                ipe_neighbor=neighbor(2,mi1,mi2,mi3,igrid)
                ic1=1+modulo(node(pig1_,igrid)-1,2)
                ic2=1+modulo(node(pig2_,igrid)-1,2)
                ic3=1+modulo(node(pig3_,igrid)-1,2);
                inc1=-2*pi1+ic1;inc2=-2*pi2+ic2;inc3=-2*pi3+ic3;
                inc1=-2*mi1+ic1;inc2=-2*mi2+ic2;inc3=-2*mi3+ic3;
                itag_cc=4**3*(ineighbor-1)+inc1*4**(1-1)+inc2*4&
                   **(2-1)+inc3*4**(3-1)
                ! Reshape to buffer and send
                isend_cc=isend_cc+1
                ibuf_cc_send_next=ibuf_cc_send+isize_stg(3)
                sendbuffer_cc(ibuf_cc_send:ibuf_cc_send_next-1)=&
                
                reshape(pflux(iside,3,igrid)%edge,shape=(/isize_stg(3)/))
                call MPI_ISEND(sendbuffer_cc(ibuf_cc_send),isize_stg(3),&
                               MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,&
                               icomm,sendrequest_stg(isend_cc),ierrmpi)
                ibuf_cc_send=ibuf_cc_send_next
         end if
              end if
            end do
            
         end if
      end do
   end select
end do

end subroutine sendflux
!=============================================================================
