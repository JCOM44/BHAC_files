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
subroutine setdt

! setdt  - set dt for all levels between levmin and levmax. 
!          dtpar>0  --> use fixed dtpar for all level
!          dtpar<=0 --> determine CFL limited timestep 
!        - set dtimpl

use mod_amrvacdef

integer :: iigrid, igrid, ncycle, ncycle2, ifile
double precision :: dtnew, qdtnew, dtmin_mype, factor, dx1,dx2,dx3, dxmin1,&
   dxmin2,dxmin3
double precision :: qdtimpl,dtimpl_mype

double precision :: dtmax, dxmin
integer,save :: stepflag
!----------------------------------------------------------------------------

if(it==0) stepflag = 0

if(sourceimpl) then
   dtimpl_mype=bigdouble
!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid,qdtnew,qdtimpl,&
!$OMP& dx1,dx2,dx3) REDUCTION(min:dtimpl_mype)
   do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      qdtimpl=bigdouble
      dx1=rnode(rpdx1_,igrid);dx2=rnode(rpdx2_,igrid);dx3=rnode(rpdx3_,igrid);
      saveigrid = igrid
      if (.not.slab) mygeo => pgeo(igrid)
      if (B0field) myB0 => pB0_cell(igrid)
      if (covariant) myM => mygeo%m
      call getdt_impl(pw(igrid)%w,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
         ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,qdtnew,dx1,dx2,dx3,&
         px(igrid)%x)
      qdtimpl=min(qdtimpl,qdtnew)
      dtimpl_mype=min(dtimpl_mype,qdtimpl)
   end do
!$OMP END PARALLEL DO
   call MPI_ALLREDUCE(dtimpl_mype,dtimpl,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
       icomm,ierrmpi)
endif

if (dtpar<=zero) then
   dtmin_mype=bigdouble
   cmax_mype = zero
 !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(igrid,qdtnew,dtnew,dx1,dx2,dx3) REDUCTION(min:dtmin_mype)
   !REDUCTION(max:cmax_mype) ! this does not seem to be doing it
   do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      dtnew=bigdouble
      dx1=rnode(rpdx1_,igrid);dx2=rnode(rpdx2_,igrid);dx3=rnode(rpdx3_,igrid);
      call set_tmpGlobals(igrid)

      call getdt_courant(pw(igrid)%w,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
         ixGhi3,ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,qdtnew,dx1,dx2,dx3,&
         px(igrid)%x)

      dtnew=min(dtnew,qdtnew)
      
      call getdt(pw(igrid)%w,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixMlo1,&
         ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,qdtnew,dx1,dx2,dx3,px(igrid)%x)
      
      dtnew=min(dtnew,qdtnew)
      
      call getdt_special(pw(igrid)%w,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
         ixGhi3,ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,qdtnew,dx1,dx2,dx3,&
         px(igrid)%x)
      
      dtnew=min(dtnew,qdtnew)
      dtmin_mype=min(dtmin_mype,dtnew)
      
      dt_grid(igrid)=dtnew
      
   end do
!$OMP END PARALLEL DO
else
   dtmin_mype=dtpar
end if

if (dtmin_mype<dtmin) then
   write(unitterm,*)"Warning: Time step too small!", dtmin_mype
   write(unitterm,*)"on processor:", mype
   write(unitterm,*)"at time:", t," step:", it
   call mpistop("too small timestep")
end if

if (slowsteps>it-itmin+1) then
   factor=one-(one-dble(it-itmin+1)/dble(slowsteps))**2
   if (time_accurate) then
      dtmin_mype=dtmin_mype*factor
   else
!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
      do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
         dt_grid(igrid)=dt_grid(igrid)*factor
      end do
!$OMP END PARALLEL DO
   end if
end if


if (time_accurate) then

   if( stepflag<1.and.mype==0) then
      if(any(dtsave(1:nfile)<dtmin_mype )) then
         write(unitterm,1001) dtmin_mype, dtsave(1:nfile)
         stepflag = 1     
      endif
   endif   
 
   if (tmaxexact) then
      dtmin_mype=min(dtmin_mype,tmax-t)
   end if

   if(any(dtsave(1:nfile)<bigdouble).or.any(tsave(isavet(1:nfile),1:nfile)&
      <bigdouble))then
      dtmax = minval((int(t/dtsave(1:nfile))+1)*dtsave(1:nfile))-t
      do ifile=1,nfile
         dtmax = min(tsave(isavet(ifile),ifile)-t,dtmax)
      end do
      if(dtmax<dtmin_mype .and. dtmax > smalldouble)then 
        dtmin_mype=min(dtmin_mype,dtmax)
      end if      
   end if

   if (dtpar<=zero) then
      call MPI_ALLREDUCE(dtmin_mype,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN, icomm,&
         ierrmpi)
   else
      dt=dtmin_mype
   end if
   
   if(sourceimpl)then
   if(sourceimplcycle)then
     ncycle2=ceiling(dt/dtimpl)
     if(ncycle2<=1)then
     !  write(unitterm,*)'implicit time step not smaller than CFL limit'
     !  write(unitterm,*)'dt =',dt,' dtimpl=',dtimpl
     !  call mpistop("no need for implicit cycling")
     else
       if(sourceparasts)then
         ncycle=floor(dsqrt(dble(ncycle2)))
         !!if(ncycle*dtimpl>dt)call mpistop("increased timestep!")
         dt=ncycle*dtimpl 
       else
         ncycle=ncycle2
       endif
       if (ncycle>ncyclemax) then
          if(mype==0) print *,'too many subcycles, reducing dt to',ncyclemax,'dt_impl!!'
          dt=ncyclemax*dtimpl
       endif
     endif
   endif
   endif


!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
   do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      dt_grid(igrid)=dt
   end do
!$OMP END PARALLEL DO
     

! global Lax-Friedrich finite difference flux splitting needs fastest wave-speed
! so does GLM: 
call MPI_ALLREDUCE(cmax_mype,cmax_global,1,MPI_DOUBLE_PRECISION,MPI_MAX,icomm,&
   ierrmpi)

end if


1001 format(' Warning: timesteps: ',1x,1pe12.5,' exceeding output intervals ',&
   2(1x,1pe12.5))

end subroutine setdt
!=============================================================================
subroutine getdt_courant(w,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
   ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,dtnew,dx1,dx2,dx3,x)

! compute CFL limited dt (for variable time stepping)

use mod_amrvacdef
 
integer, intent(in) :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
    ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
double precision, intent(in) :: dx1,dx2,dx3, x(ixGmin1:ixGmax1,&
   ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:ndim)
double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
   ixGmin3:ixGmax3,1:nw), dtnew

integer :: idims
double precision :: courantmax, dxinv(1:ndim), courantmaxtot, courantmaxtots
double precision :: cmax(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
    cmaxtot(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3), tmp(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2,ixGlo3:ixGhi3)
!-----------------------------------------------------------------------------
dtnew=bigdouble

courantmax=zero
courantmaxtot=zero
courantmaxtots=zero



dxinv(1)=one/dx1;dxinv(2)=one/dx2;dxinv(3)=one/dx3;

cmaxtot(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=zero

if (nwaux>0) call getaux(.true.,w,x,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
   ixGmax3,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,'getdt_courant')

do idims=1,ndim
   call getcmax(w,x,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,ixmin1,&
      ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idims,cmax,tmp,.false.)
   !$OMP atomic
   cmax_mype = max(cmax_mype,maxval(cmax(ixmin1:ixmax1,ixmin2:ixmax2,&
      ixmin3:ixmax3)))
   if (.not.slab) then
      tmp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=cmax(ixmin1:ixmax1,&
         ixmin2:ixmax2,ixmin3:ixmax3)/mygeo%dx(ixmin1:ixmax1,ixmin2:ixmax2,&
         ixmin3:ixmax3,idims)
      cmaxtot(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=cmaxtot&
         (ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)+tmp(ixmin1:ixmax1,&
         ixmin2:ixmax2,ixmin3:ixmax3)
      courantmax=max(courantmax,maxval(tmp(ixmin1:ixmax1,ixmin2:ixmax2,&
         ixmin3:ixmax3)))
   else
      cmaxtot(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=cmaxtot&
         (ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)+cmax(ixmin1:ixmax1,&
         ixmin2:ixmax2,ixmin3:ixmax3)*dxinv(idims)
      courantmax=max(courantmax,maxval(cmax(ixmin1:ixmax1,ixmin2:ixmax2,&
         ixmin3:ixmax3)*dxinv(idims)))
   end if
   courantmaxtot=courantmaxtot+courantmax
end do


select case (typecourant)
case ('minimum')
   ! courantmax='max(c/dx)'
   if (courantmax>smalldouble)     dtnew=min(dtnew,courantpar/courantmax)
case ('summax')
   ! courantmaxtot='summed max(c/dx)'
   if (courantmaxtot>smalldouble)  dtnew=min(dtnew,courantpar/courantmaxtot)
case ('maxsum')
   ! courantmaxtots='max(summed c/dx)'
   courantmaxtots=max(courantmaxtots,maxval(cmaxtot(ixmin1:ixmax1,&
      ixmin2:ixmax2,ixmin3:ixmax3)))
   if (courantmaxtots>smalldouble) dtnew=min(dtnew,courantpar/courantmaxtots)
case default
   write(unitterm,*)'Unknown typecourant=',typecourant
   call mpistop("Error from getdt_courant: no such typecourant!")
end select

end subroutine getdt_courant
!=============================================================================
subroutine getresidual(iit)

! compute residual for steady state calculations

use mod_amrvacdef
 
integer, intent(in) :: iit

integer :: iigrid,igrid,iw
double precision :: wnrm2_send(1:nwflux),wnrm2_recv(1:nwflux)
double precision :: wnrm2localgrids(1:nwflux),residlocalgrids(1:nwflux)
double precision :: resid_send(1:nwflux),resid_recv(1:nwflux)
!----------------------------------------------------------------------------

select case(typeresid)
  case('relative')
    wnrm2localgrids(1:nwflux)=zero
 !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid,iw) REDUCTION(+:wnrm2localgrids)
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      do iw=1,nwflux
        wnrm2localgrids(iw)=wnrm2localgrids(iw)+sum(pw(igrid)%w(ixMlo1:ixMhi1,&
           ixMlo2:ixMhi2,ixMlo3:ixMhi3,iw)**2)
      enddo
    end do 
!$OMP END PARALLEL DO

    wnrm2_send(1:nwflux)=wnrm2localgrids(1:nwflux)

    call MPI_ALLREDUCE(wnrm2_send,wnrm2_recv,nwflux,MPI_DOUBLE_PRECISION,&
        MPI_SUM,icomm,ierrmpi)

    do iw=1,nwflux
      if(wnrm2_recv(iw)<smalldouble) wnrm2_recv(iw)=one
    enddo
  case('absolute')
    wnrm2_recv(1:nwflux)=one
  case default
    call mpistop('no such typeresid')
end select

residlocalgrids(1:nwflux)=zero
 !$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid,iw) REDUCTION(+:residlocalgrids)
do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
    do iw=1,nwflux
      residlocalgrids(iw)=residlocalgrids(iw) +sum(pwres(igrid)%w&
         (ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,iw)**2)/wnrm2_recv(iw)
    enddo
end do 
!$OMP END PARALLEL DO

resid_send(1:nwflux)=residlocalgrids(1:nwflux)

call MPI_ALLREDUCE(resid_send,resid_recv,nwflux,MPI_DOUBLE_PRECISION, MPI_SUM,&
   icomm,ierrmpi)


residual=sqrt(sum(resid_recv(1:nwflux))/(nwflux))


end subroutine getresidual
!=============================================================================
