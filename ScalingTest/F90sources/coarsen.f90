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
subroutine coarsen_grid_siblings(igrid,ipe,child_igrid,child_ipe,active)

use mod_amrvacdef

integer, intent(in) :: igrid, ipe
integer, dimension(2,2,2), intent(in) :: child_igrid, child_ipe
logical, intent(in) :: active

integer :: igridFi, ipeFi, ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,&
   ixComax3, ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3,&
    ixCoMmin1,ixCoMmin2,ixCoMmin3,ixCoMmax1,ixCoMmax2,ixCoMmax3, ic1,ic2,ic3

integer :: idir

!-----------------------------------------------------------------------------
if(addmpibarrier) call MPI_BARRIER(icomm,ierrmpi)

if (ipe==mype) call alloc_node(igrid)

! New passive cell, coarsen from initial condition:
if (.not. active) then

   if (ipe == mype) then
      call initial_condition(igrid)
      
      do ic3=1,2
      do ic2=1,2
      do ic1=1,2
      igridFi=child_igrid(ic1,ic2,ic3)
      ipeFi=child_ipe(ic1,ic2,ic3)
      if (ipeFi==mype) then
         ! remove solution space of child      
         call dealloc_node(igridFi)
      end if
      end do
      end do
      end do
      
   end if

   return
end if


do ic3=1,2
do ic2=1,2
do ic1=1,2
   igridFi=child_igrid(ic1,ic2,ic3)
   ipeFi=child_ipe(ic1,ic2,ic3)

   if (ipeFi==mype) then
      dxlevel(1)=rnode(rpdx1_,igridFi);dxlevel(2)=rnode(rpdx2_,igridFi)
      dxlevel(3)=rnode(rpdx3_,igridFi);
      if (ipe==mype) then
         ixComin1=ixMlo1+(ic1-1)*(ixMhi1-ixMlo1+1)/2
         ixComin2=ixMlo2+(ic2-1)*(ixMhi2-ixMlo2+1)/2
         ixComin3=ixMlo3+(ic3-1)*(ixMhi3-ixMlo3+1)/2;
         ixComax1=ixMhi1+(ic1-2)*(ixMhi1-ixMlo1+1)/2
         ixComax2=ixMhi2+(ic2-2)*(ixMhi2-ixMlo2+1)/2
         ixComax3=ixMhi3+(ic3-2)*(ixMhi3-ixMlo3+1)/2;

         call coarsen_grid(ps(igridFi),px(igridFi)%x,ixGlo1,ixGlo2,ixGlo3,&
            ixGhi1,ixGhi2,ixGhi3,ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,&
            ps(igrid),px(igrid)%x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
             ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,ixComax3,&
            pgeo(igridFi),pgeo(igrid),restrictprimitive,.false.)

         ! remove solution space of child
         call dealloc_node(igridFi)
      else
         ixCoGmin1=1;ixCoGmin2=1;ixCoGmin3=1;
        ixCoGmax1=ixGhi1/2+dixB;ixCoGmax2=ixGhi2/2+dixB
        ixCoGmax3=ixGhi3/2+dixB;
         ixCoMmin1=ixCoGmin1+dixB;ixCoMmin2=ixCoGmin2+dixB
         ixCoMmin3=ixCoGmin3+dixB;ixCoMmax1=ixCoGmax1-dixB
         ixCoMmax2=ixCoGmax2-dixB;ixCoMmax3=ixCoGmax3-dixB;
         call coarsen_grid(ps(igridFi),px(igridFi)%x,ixGlo1,ixGlo2,ixGlo3,&
            ixGhi1,ixGhi2,ixGhi3,ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,&
            psCoarse(igridFi),pxCoarse(igridFi)%x, ixCoGmin1,ixCoGmin2,&
            ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3,ixCoMmin1,ixCoMmin2,&
            ixCoMmin3,ixCoMmax1,ixCoMmax2,ixCoMmax3,pgeo(igridFi),&
            pgeoCoarse(igridFi),restrictprimitive,.false.)

         itag=ipeFi*ngridshi+igridFi
         isend=isend+1
         call MPI_ISEND(pwCoarse(igridFi)%w,1,type_coarse_block,ipe,itag,&
             icomm,sendrequest(isend),ierrmpi)

         do idir=1,ndim
         itag_stg=(npe+ipeFi+1)*ngridshi+igridFi*(3-1+idir)
         call MPI_ISEND(psCoarse(igridFi)%ws%w,1,type_coarse_block_stg(idir,&
            ic1,ic2,ic3),ipe,itag_stg, &
                        icomm,sendrequest_stg(isend),ierrmpi)
         end do

         
      end if
   else
      if (ipe==mype) then
         itag=ipeFi*ngridshi+igridFi
         irecv=irecv+1
         call MPI_IRECV(pw(igrid)%w,1,type_sub_block(ic1,ic2,ic3),ipeFi,itag,&
             icomm,recvrequest(irecv),ierrmpi)

         do idir=1,ndim
         itag_stg=(npe+ipeFi+1)*ngridshi+igridFi*(3-1+idir)
         call MPI_IRECV(ps(igrid)%ws%w,1,type_sub_block_stg(idir,ic1,ic2,ic3),&
            ipeFi,itag_stg, &
                        icomm,recvrequest_stg(irecv),ierrmpi)
         end do

         
      end if
   end if
end do
end do
end do

if(addmpibarrier) call MPI_BARRIER(icomm,ierrmpi)
end subroutine coarsen_grid_siblings
!=============================================================================
subroutine coarsen_grid(sFi,xFi,ixFiGmin1,ixFiGmin2,ixFiGmin3,ixFiGmax1,&
   ixFiGmax2,ixFiGmax3,ixFimin1,ixFimin2,ixFimin3,ixFimax1,ixFimax2,ixFimax3,&
   sCo,xCo,ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3,&
   ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,ixComax3,pgeogrid,&
   pgeoCoarsegrid,coarsenprim,keepFi)

! coarsen by 2 in every direction - conservatively
use mod_amrvacdef

integer, intent(in)             :: ixFiGmin1,ixFiGmin2,ixFiGmin3,ixFiGmax1,&
   ixFiGmax2,ixFiGmax3, ixFimin1,ixFimin2,ixFimin3,ixFimax1,ixFimax2,ixFimax3,&
    ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3, ixComin1,&
   ixComin2,ixComin3,ixComax1,ixComax2,ixComax3
double precision, intent(inout) :: xFi(ixFiGmin1:ixFiGmax1,&
   ixFiGmin2:ixFiGmax2,ixFiGmin3:ixFiGmax3,1:ndim)
double precision,intent(inout)  :: xCo(ixCoGmin1:ixCoGmax1,&
   ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,1:ndim)
type(state), intent(inout)      :: sFi, sCo
type(geoalloc)                  :: pgeogrid, pgeoCoarsegrid
logical, intent(in)             :: coarsenprim, keepFi

! .. local ..
integer :: ixCo1,ixCo2,ixCo3, ixFi1,ixFi2,ixFi3, iw
double precision :: CoFiratio
!-----------------------------------------------------------------------------
associate(wFi=>sFi%w%w(ixFiGmin1:ixFiGmax1,ixFiGmin2:ixFiGmax2,&
   ixFiGmin3:ixFiGmax3,1:nw), wCo=>sCo%w%w(ixCoGmin1:ixCoGmax1,&
   ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,1:nw))


staggered : associate(wFis=>sFi%ws%w,wCos=>sCo%ws%w)


if (covariant) myM=>pgeogrid%m
if (amrentropy) then
   call e_to_rhos(ixFiGmin1,ixFiGmin2,ixFiGmin3,ixFiGmax1,ixFiGmax2,ixFiGmax3,&
      ixFimin1,ixFimin2,ixFimin3,ixFimax1,ixFimax2,ixFimax3,wFi,xFi)
else if (coarsenprim) then
   call primitive(ixFiGmin1,ixFiGmin2,ixFiGmin3,ixFiGmax1,ixFiGmax2,ixFiGmax3,&
      ixFimin1,ixFimin2,ixFimin3,ixFimax1,ixFimax2,ixFimax3,wFi,xFi)
end if


if (slab) then
   CoFiratio=one/dble(2**ndim)
   do iw=1,nw
      do ixCo3 = ixComin3,ixComax3
         ixFi3=2*(ixCo3-ixComin3)+ixFimin3
      do ixCo2 = ixComin2,ixComax2
         ixFi2=2*(ixCo2-ixComin2)+ixFimin2
      do ixCo1 = ixComin1,ixComax1
         ixFi1=2*(ixCo1-ixComin1)+ixFimin1
         wCo(ixCo1,ixCo2,ixCo3,iw)=sum(wFi(ixFi1:ixFi1+1,ixFi2:ixFi2+1,&
            ixFi3:ixFi3+1,iw))*CoFiratio
      end do
      end do
      end do
   end do
else
   do iw=1,nw
      do ixCo3 = ixComin3,ixComax3
         ixFi3=2*(ixCo3-ixComin3)+ixFimin3
      do ixCo2 = ixComin2,ixComax2
         ixFi2=2*(ixCo2-ixComin2)+ixFimin2
      do ixCo1 = ixComin1,ixComax1
         ixFi1=2*(ixCo1-ixComin1)+ixFimin1
         wCo(ixCo1,ixCo2,ixCo3,iw)= sum(pgeogrid%dvolume(ixFi1:ixFi1+1,&
            ixFi2:ixFi2+1,ixFi3:ixFi3+1)*wFi(ixFi1:ixFi1+1,ixFi2:ixFi2+1,&
            ixFi3:ixFi3+1,iw)) /pgeoCoarsegrid%dvolume(ixCo1,ixCo2,ixCo3)
      end do
      end do
      end do
   end do

   

   do iw=1,nws

       select case(iw)
       case(1)
       ! Start one layer before
       do ixCo3 = ixComin3-kr(3,1),ixComax3
       ixFi3=2*(ixCo3-ixComin3+kr(3,1))+ixFimin3-kr(3,1)
do ixCo2 = ixComin2-kr(2,1),ixComax2
       ixFi2=2*(ixCo2-ixComin2+kr(2,1))+ixFimin2-kr(2,1)
do ixCo1 = ixComin1-kr(1,1),ixComax1
       ixFi1=2*(ixCo1-ixComin1+kr(1,1))+ixFimin1-kr(1,1)
       
            ! This if statement catches the axis where surface is zero:
            if (pgeoCoarsegrid%surfaceC1(ixCo1,ixCo2,ixCo3) .gt. &
               1.0d-9*pgeoCoarsegrid%dvolume(ixCo1,ixCo2,ixCo3)) then !Normal case
              wCos(ixCo1,ixCo2,ixCo3,iw)= &
                sum(pgeogrid%surfaceC1(ixFi1:ixFi1+1-kr(1,1),&
                   ixFi2:ixFi2+1-kr(1,2),ixFi3:ixFi3+1-kr(1,&
                   3))*wFis(ixFi1:ixFi1+1-kr(1,1),ixFi2:ixFi2+1-kr(1,2),&
                   ixFi3:ixFi3+1-kr(1,3),iw)) &
                 /pgeoCoarsegrid%surfaceC1(ixCo1,ixCo2,ixCo3)
            else ! On axis
               wCos(ixCo1,ixCo2,ixCo3,iw)= zero
            end if

       end do
end do
end do
       case(2)
       ! Start one layer before
       do ixCo3 = ixComin3-kr(3,2),ixComax3
       ixFi3=2*(ixCo3-ixComin3+kr(3,2))+ixFimin3-kr(3,2)
do ixCo2 = ixComin2-kr(2,2),ixComax2
       ixFi2=2*(ixCo2-ixComin2+kr(2,2))+ixFimin2-kr(2,2)
do ixCo1 = ixComin1-kr(1,2),ixComax1
       ixFi1=2*(ixCo1-ixComin1+kr(1,2))+ixFimin1-kr(1,2)
       
            ! This if statement catches the axis where surface is zero:
            if (pgeoCoarsegrid%surfaceC2(ixCo1,ixCo2,ixCo3) .gt. &
               1.0d-9*pgeoCoarsegrid%dvolume(ixCo1,ixCo2,ixCo3)) then !Normal case
              wCos(ixCo1,ixCo2,ixCo3,iw)= &
                sum(pgeogrid%surfaceC2(ixFi1:ixFi1+1-kr(2,1),&
                   ixFi2:ixFi2+1-kr(2,2),ixFi3:ixFi3+1-kr(2,&
                   3))*wFis(ixFi1:ixFi1+1-kr(2,1),ixFi2:ixFi2+1-kr(2,2),&
                   ixFi3:ixFi3+1-kr(2,3),iw)) &
                 /pgeoCoarsegrid%surfaceC2(ixCo1,ixCo2,ixCo3)
            else ! On axis
               wCos(ixCo1,ixCo2,ixCo3,iw)= zero
            end if

       end do
end do
end do
       case(3)
       ! Start one layer before
       do ixCo3 = ixComin3-kr(3,3),ixComax3
       ixFi3=2*(ixCo3-ixComin3+kr(3,3))+ixFimin3-kr(3,3)
do ixCo2 = ixComin2-kr(2,3),ixComax2
       ixFi2=2*(ixCo2-ixComin2+kr(2,3))+ixFimin2-kr(2,3)
do ixCo1 = ixComin1-kr(1,3),ixComax1
       ixFi1=2*(ixCo1-ixComin1+kr(1,3))+ixFimin1-kr(1,3)
       
            ! This if statement catches the axis where surface is zero:
            if (pgeoCoarsegrid%surfaceC3(ixCo1,ixCo2,ixCo3) .gt. &
               1.0d-9*pgeoCoarsegrid%dvolume(ixCo1,ixCo2,ixCo3)) then !Normal case
              wCos(ixCo1,ixCo2,ixCo3,iw)= &
                sum(pgeogrid%surfaceC3(ixFi1:ixFi1+1-kr(3,1),&
                   ixFi2:ixFi2+1-kr(3,2),ixFi3:ixFi3+1-kr(3,&
                   3))*wFis(ixFi1:ixFi1+1-kr(3,1),ixFi2:ixFi2+1-kr(3,2),&
                   ixFi3:ixFi3+1-kr(3,3),iw)) &
                 /pgeoCoarsegrid%surfaceC3(ixCo1,ixCo2,ixCo3)
            else ! On axis
               wCos(ixCo1,ixCo2,ixCo3,iw)= zero
            end if

       end do
end do
end do
      
       end select

   end do
   ! average to fill cell-centred values
   call faces2centers(ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,ixComax3,&
      sCo)
   
   
end if


if (amrentropy) then
   if (keepFi) then
      if (covariant) myM=>pgeogrid%m
      call rhos_to_e(ixFiGmin1,ixFiGmin2,ixFiGmin3,ixFiGmax1,ixFiGmax2,&
         ixFiGmax3,ixFimin1,ixFimin2,ixFimin3,ixFimax1,ixFimax2,ixFimax3,wFi,&
         xFi)
   end if
   if (covariant) myM=>pgeoCoarsegrid%m
   call rhos_to_e(ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3,&
      ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,ixComax3,wCo,xCo)
else if (coarsenprim) then
   if (keepFi) then
      if (covariant) myM=>pgeogrid%m
      call conserve(ixFiGmin1,ixFiGmin2,ixFiGmin3,ixFiGmax1,ixFiGmax2,&
         ixFiGmax3,ixFimin1,ixFimin2,ixFimin3,ixFimax1,ixFimax2,ixFimax3,wFi,&
         xFi,patchfalse)
   end if
   if (covariant) myM=>pgeoCoarsegrid%m
   call conserve(ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3,&
      ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,ixComax3,wCo,xCo,&
      patchfalse)
end if



end associate staggered



end associate
end subroutine coarsen_grid
!=============================================================================
