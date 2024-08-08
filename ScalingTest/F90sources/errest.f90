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
subroutine errest
 
use mod_forest, only: refine, buffer
use mod_amrvacdef

integer :: igrid, iigrid, ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,&
   ixCoGmax3
double precision :: factor
!-----------------------------------------------------------------------------
if (igridstail==0) return

select case (errorestimate)
case (0) 
   ! all refinement solely based on user routine specialrefine_grid
case (1)
   call mpistop("errest: Richardson no longer supported")
   ! Richardson procedure: compare coarse-integrate versus integrate-coarse
   ! done with low order, dimensionally unsplit scheme typelow1
   ! For error estimate: compare 1st order, coarse 2*dx,2*dt step  with
   ! 1st order dt step, since otherwise wCT can not be filled by interpolation

   ! Note:when there are sources: only unsplit sources are taken into account

   ! Note: the low order solutions are obtained with dimsplit=F.
   !       When overall scheme uses dimsplit=T
   !       and courantpar>0.5, the low order can be unstable. In
   !       that case, simplify this scheme, set low order step of size dt
   !       and compare coarse with available solution at t_n. Enforce this
   !       through `skipfinestep' 

case (2) 
   ! simply compare w_n-1 with w_n and trigger refinement on relative
   ! differences
!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call set_tmpGlobals(igrid)
      call compare1_grid(igrid,pwold(igrid)%w,pw(igrid)%w)
   end do
!$OMP END PARALLEL DO

case (3)
   ! Error estimation is based on Lohner's scheme
!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call set_tmpGlobals(igrid)
      call lohner_grid(igrid)
   end do
!$OMP END PARALLEL DO

case (4)
   ! Error estimation is based on Lohner's original scheme
!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call set_tmpGlobals(igrid)
      call lohner_orig_grid(igrid)
   end do
!$OMP END PARALLEL DO


case default
   call mpistop("Unknown error estimator")
end select

!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
do iigrid=1,igridstail; igrid=igrids(iigrid);
   call forcedrefine_grid(igrid,pw(igrid)%w)
end do
!$OMP END PARALLEL DO

! add the buffer region at the end
if (nbufferx1/=0.or.nbufferx2/=0.or.nbufferx3/=0) then 

   call MPI_ALLGATHER(MPI_IN_PLACE,ngridshi,MPI_LOGICAL,refine,ngridshi,&
       MPI_LOGICAL,icomm,ierrmpi)
   
!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      if (refine(igrid,mype) .and. .not.buffer(igrid,mype)) then
         ! We discard the number of cells for the buffer
         ! (meaning always the entire neighbor block is considered)
         ! so this implementation is a tad silly
         call refinebuffer(igrid,.not. patchfalse)
      end if
   end do
!$OMP END PARALLEL DO

   call MPI_ALLREDUCE(MPI_IN_PLACE,refine,ngridshi*npe,MPI_LOGICAL,MPI_LOR,&
       icomm,ierrmpi)

   buffer=.false.
end if

end subroutine errest
!=============================================================================
subroutine lohner_grid(igrid)
use mod_forest, only: coarsen, refine
use mod_amrvacdef

integer, intent(in) :: igrid

integer :: iiflag, iflag, idims, idims2, level
integer :: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3, hxmin1,hxmin2,hxmin3,&
   hxmax1,hxmax2,hxmax3, jxmin1,jxmin2,jxmin3,jxmax1,jxmax2,jxmax3, h2xmin1,&
   h2xmin2,h2xmin3,h2xmax1,h2xmax2,h2xmax3, j2xmin1,j2xmin2,j2xmin3,j2xmax1,&
   j2xmax2,j2xmax3, ix1,ix2,ix3
double precision :: epsilon, tolerance
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3) :: &
   numerator, denominator, error
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3) :: tmp,&
    tmp1, tmp2
logical, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3) :: refineflag,&
    coarsenflag
!-----------------------------------------------------------------------------
epsilon=1.0d-6
level=node(plevel_,igrid)
ixmin1=ixMlo1-1;ixmin2=ixMlo2-1;ixmin3=ixMlo3-1;ixmax1=ixMhi1+1
ixmax2=ixMhi2+1;ixmax3=ixMhi3+1;

error=zero
do iiflag=1,flags(nflag_); iflag=flags(iiflag);
   numerator=zero
   if (iflag>nw)call specialvarforerrest(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
      ixGhi3,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,iflag,pw(igrid)%w,tmp1)
   do idims=1,ndim
      hxmin1=ixmin1-kr(1,idims);hxmin2=ixmin2-kr(2,idims)
      hxmin3=ixmin3-kr(3,idims);hxmax1=ixmax1-kr(1,idims)
      hxmax2=ixmax2-kr(2,idims);hxmax3=ixmax3-kr(3,idims);
      jxmin1=ixmin1+kr(1,idims);jxmin2=ixmin2+kr(2,idims)
      jxmin3=ixmin3+kr(3,idims);jxmax1=ixmax1+kr(1,idims)
      jxmax2=ixmax2+kr(2,idims);jxmax3=ixmax3+kr(3,idims);
      if (iflag<=nw) then
        if (logflag(iiflag)) then
          tmp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=dlog10(pw(igrid)%w&
             (jxmin1:jxmax1,jxmin2:jxmax2,jxmin3:jxmax3,iflag))-dlog10(pw&
             (igrid)%w(hxmin1:hxmax1,hxmin2:hxmax2,hxmin3:hxmax3,iflag))
        else
          tmp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=pw(igrid)%w&
             (jxmin1:jxmax1,jxmin2:jxmax2,jxmin3:jxmax3,iflag)-pw(igrid)%w&
             (hxmin1:hxmax1,hxmin2:hxmax2,hxmin3:hxmax3,iflag)
        end if
      else
        if (logflag(iiflag)) then
          tmp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=dlog10(tmp1&
             (jxmin1:jxmax1,jxmin2:jxmax2,jxmin3:jxmax3))-dlog10(tmp1&
             (hxmin1:hxmax1,hxmin2:hxmax2,hxmin3:hxmax3))
        else
          tmp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=tmp1(jxmin1:jxmax1,&
             jxmin2:jxmax2,jxmin3:jxmax3)-tmp1(hxmin1:hxmax1,hxmin2:hxmax2,&
             hxmin3:hxmax3)
        end if
      end if
      do idims2=1,ndim
         h2xmin1=ixMlo1-kr(1,idims2);h2xmin2=ixMlo2-kr(2,idims2)
         h2xmin3=ixMlo3-kr(3,idims2);h2xmax1=ixMhi1-kr(1,idims2)
         h2xmax2=ixMhi2-kr(2,idims2);h2xmax3=ixMhi3-kr(3,idims2);
         j2xmin1=ixMlo1+kr(1,idims2);j2xmin2=ixMlo2+kr(2,idims2)
         j2xmin3=ixMlo3+kr(3,idims2);j2xmax1=ixMhi1+kr(1,idims2)
         j2xmax2=ixMhi2+kr(2,idims2);j2xmax3=ixMhi3+kr(3,idims2);
         numerator=numerator+(tmp(j2xmin1:j2xmax1,j2xmin2:j2xmax2,&
            j2xmin3:j2xmax3)-tmp(h2xmin1:h2xmax1,h2xmin2:h2xmax2,&
            h2xmin3:h2xmax3))**2.0d0
      end do
   end do
   denominator=zero
   do idims=1,ndim
      if (iflag<=nw) then
         if (logflag(iiflag)) then
          tmp=dabs(dlog10(pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             ixGlo3:ixGhi3,iflag)))
         else
          tmp=dabs(pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
             iflag))
         end if
      else
         if (logflag(iiflag)) then
          tmp=dabs(dlog10(tmp1(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)))
         else
          tmp=dabs(tmp1(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3))
         end if
      end if
      hxmin1=ixmin1-kr(1,idims);hxmin2=ixmin2-kr(2,idims)
      hxmin3=ixmin3-kr(3,idims);hxmax1=ixmax1-kr(1,idims)
      hxmax2=ixmax2-kr(2,idims);hxmax3=ixmax3-kr(3,idims);
      jxmin1=ixmin1+kr(1,idims);jxmin2=ixmin2+kr(2,idims)
      jxmin3=ixmin3+kr(3,idims);jxmax1=ixmax1+kr(1,idims)
      jxmax2=ixmax2+kr(2,idims);jxmax3=ixmax3+kr(3,idims);
      tmp2(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=tmp(jxmin1:jxmax1,&
         jxmin2:jxmax2,jxmin3:jxmax3)+tmp(hxmin1:hxmax1,hxmin2:hxmax2,&
         hxmin3:hxmax3)
      hxmin1=ixMlo1-2*kr(1,idims);hxmin2=ixMlo2-2*kr(2,idims)
      hxmin3=ixMlo3-2*kr(3,idims);hxmax1=ixMhi1-2*kr(1,idims)
      hxmax2=ixMhi2-2*kr(2,idims);hxmax3=ixMhi3-2*kr(3,idims);
      jxmin1=ixMlo1+2*kr(1,idims);jxmin2=ixMlo2+2*kr(2,idims)
      jxmin3=ixMlo3+2*kr(3,idims);jxmax1=ixMhi1+2*kr(1,idims)
      jxmax2=ixMhi2+2*kr(2,idims);jxmax3=ixMhi3+2*kr(3,idims);
      if (iflag<=nw) then
        if (logflag(iiflag)) then
          tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3)=dabs(dlog10(pw&
             (igrid)%w(jxmin1:jxmax1,jxmin2:jxmax2,jxmin3:jxmax3,&
             iflag))-dlog10(pw(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
             ixMlo3:ixMhi3,iflag))) +dabs(dlog10(pw(igrid)%w(ixMlo1:ixMhi1,&
             ixMlo2:ixMhi2,ixMlo3:ixMhi3,iflag))-dlog10(pw(igrid)%w&
             (hxmin1:hxmax1,hxmin2:hxmax2,hxmin3:hxmax3,iflag)))
        else
           tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3)=dabs(pw(igrid)%w&
              (jxmin1:jxmax1,jxmin2:jxmax2,jxmin3:jxmax3,iflag)-pw(igrid)%w&
              (ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,iflag)) &
              +dabs(pw(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,&
              iflag)-pw(igrid)%w(hxmin1:hxmax1,hxmin2:hxmax2,hxmin3:hxmax3,&
              iflag))
        end if
      else
        if (logflag(iiflag)) then
          tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3)=dabs(dlog10(tmp1&
             (jxmin1:jxmax1,jxmin2:jxmax2,jxmin3:jxmax3))-dlog10(tmp1&
             (ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3))) &
             +dabs(dlog10(tmp1(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
             ixMlo3:ixMhi3))-dlog10(tmp1(hxmin1:hxmax1,hxmin2:hxmax2,&
             hxmin3:hxmax3)))
        else
           tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3)=dabs(tmp1&
              (jxmin1:jxmax1,jxmin2:jxmax2,jxmin3:jxmax3)-tmp1(ixMlo1:ixMhi1,&
              ixMlo2:ixMhi2,ixMlo3:ixMhi3)) +dabs(tmp1(ixMlo1:ixMhi1,&
              ixMlo2:ixMhi2,ixMlo3:ixMhi3)-tmp1(hxmin1:hxmax1,hxmin2:hxmax2,&
              hxmin3:hxmax3))
        end if
      end if
      do idims2=1,ndim
         h2xmin1=ixMlo1-kr(1,idims2);h2xmin2=ixMlo2-kr(2,idims2)
         h2xmin3=ixMlo3-kr(3,idims2);h2xmax1=ixMhi1-kr(1,idims2)
         h2xmax2=ixMhi2-kr(2,idims2);h2xmax3=ixMhi3-kr(3,idims2);
         j2xmin1=ixMlo1+kr(1,idims2);j2xmin2=ixMlo2+kr(2,idims2)
         j2xmin3=ixMlo3+kr(3,idims2);j2xmax1=ixMhi1+kr(1,idims2)
         j2xmax2=ixMhi2+kr(2,idims2);j2xmax3=ixMhi3+kr(3,idims2);
         denominator=denominator +(tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
            ixMlo3:ixMhi3)+amr_wavefilter(level)*(tmp2(j2xmin1:j2xmax1,&
            j2xmin2:j2xmax2,j2xmin3:j2xmax3)+tmp2(h2xmin1:h2xmax1,&
            h2xmin2:h2xmax2,h2xmin3:h2xmax3)))**2
      end do
   end do
   error=error+wflags(iiflag)*dsqrt(numerator/max(denominator,epsilon))
end do

refineflag=.false.
coarsenflag=.false.
tolerance=tol(level)
do ix3=ixMlo3,ixMhi3
do ix2=ixMlo2,ixMhi2
do ix1=ixMlo1,ixMhi1

   if (error(ix1,ix2,ix3) >= tolerance) then
      refineflag(ix1,ix2,ix3) = .true.
   else if (error(ix1,ix2,ix3) <= tolratio(level)*tolerance) then
      coarsenflag(ix1,ix2,ix3) = .true.
   end if
end do
end do
end do


if (any(refineflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3))&
   .and.level<mxnest) refine(igrid,mype)=.true.
if (all(coarsenflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3))&
   .and.level>1) coarsen(igrid,mype)=.true.

end subroutine lohner_grid
!=============================================================================
subroutine lohner_orig_grid(igrid)
use mod_forest, only: coarsen, refine
use mod_amrvacdef

integer, intent(in) :: igrid

integer :: iiflag, iflag, idims, level
integer :: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3, hxmin1,hxmin2,hxmin3,&
   hxmax1,hxmax2,hxmax3, jxmin1,jxmin2,jxmin3,jxmax1,jxmax2,jxmax3, ix1,ix2,&
   ix3
double precision :: epsilon
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3) :: &
   numerator, denominator, error
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3) :: dp,&
    dm, dref, tmp1
logical, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3) :: refineflag,&
    coarsenflag
!-----------------------------------------------------------------------------
epsilon=1.0d-6
level=node(plevel_,igrid)
ixmin1=ixMlo1;ixmin2=ixMlo2;ixmin3=ixMlo3;ixmax1=ixMhi1;ixmax2=ixMhi2
ixmax3=ixMhi3;

error=zero
do iiflag=1,flags(nflag_); iflag=flags(iiflag);
   numerator=zero
   denominator=zero
   if (iflag>nw)call specialvarforerrest(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
      ixGhi3,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,iflag,pw(igrid)%w,tmp1)
   do idims=1,ndim
      hxmin1=ixmin1-kr(1,idims);hxmin2=ixmin2-kr(2,idims)
      hxmin3=ixmin3-kr(3,idims);hxmax1=ixmax1-kr(1,idims)
      hxmax2=ixmax2-kr(2,idims);hxmax3=ixmax3-kr(3,idims);
      jxmin1=ixmin1+kr(1,idims);jxmin2=ixmin2+kr(2,idims)
      jxmin3=ixmin3+kr(3,idims);jxmax1=ixmax1+kr(1,idims)
      jxmax2=ixmax2+kr(2,idims);jxmax3=ixmax3+kr(3,idims);
      if (iflag<=nw) then
        if (logflag(iiflag)) then
          dp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=dlog10(pw(igrid)%w&
             (jxmin1:jxmax1,jxmin2:jxmax2,jxmin3:jxmax3,iflag))-dlog10(pw&
             (igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,iflag))
          dm(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=dlog10(pw(igrid)%w&
             (ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,iflag))-dlog10(pw&
             (igrid)%w(hxmin1:hxmax1,hxmin2:hxmax2,hxmin3:hxmax3,iflag))
          dref(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3)=dabs(dlog10(pw&
             (igrid)%w(jxmin1:jxmax1,jxmin2:jxmax2,jxmin3:jxmax3,&
             iflag)))+ 2.0d0 * dabs(dlog10(pw(igrid)%w(ixMlo1:ixMhi1,&
             ixMlo2:ixMhi2,ixMlo3:ixMhi3,iflag))) + dabs(dlog10(pw(igrid)%w&
             (hxmin1:hxmax1,hxmin2:hxmax2,hxmin3:hxmax3,iflag)))
        else
          dp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=pw(igrid)%w&
             (jxmin1:jxmax1,jxmin2:jxmax2,jxmin3:jxmax3,iflag)-pw(igrid)%w&
             (ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,iflag)
          dp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=pw(igrid)%w&
             (ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,iflag)-pw(igrid)%w&
             (hxmin1:hxmax1,hxmin2:hxmax2,hxmin3:hxmax3,iflag)
          dref(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3)=dabs(pw(igrid)%w&
             (jxmin1:jxmax1,jxmin2:jxmax2,jxmin3:jxmax3,iflag))+2.0d0*dabs(pw&
             (igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,&
             iflag)) +dabs(pw(igrid)%w(hxmin1:hxmax1,hxmin2:hxmax2,&
             hxmin3:hxmax3,iflag))
        end if
      else
        if (logflag(iiflag)) then
          dp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=dlog10(tmp1&
             (jxmin1:jxmax1,jxmin2:jxmax2,jxmin3:jxmax3))-dlog10(tmp1&
             (ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3))
          dm(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=dlog10(tmp1&
             (ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3))-dlog10(tmp1&
             (hxmin1:hxmax1,hxmin2:hxmax2,hxmin3:hxmax3))
          dref(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=dabs(dlog10(tmp1&
             (jxmin1:jxmax1,jxmin2:jxmax2,jxmin3:jxmax3)))+ 2.0d0 * &
             dabs(dlog10(tmp1(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3))) + &
             dabs(dlog10(tmp1(hxmin1:hxmax1,hxmin2:hxmax2,hxmin3:hxmax3)))
        else
          dp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=tmp1(jxmin1:jxmax1,&
             jxmin2:jxmax2,jxmin3:jxmax3)-tmp1(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)
          dm(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=tmp1(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)-tmp1(hxmin1:hxmax1,hxmin2:hxmax2,&
             hxmin3:hxmax3)
          dref(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=dabs(tmp1&
             (jxmin1:jxmax1,jxmin2:jxmax2,jxmin3:jxmax3))+2.0d0*dabs(tmp1&
             (ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)) +dabs(tmp1&
             (hxmin1:hxmax1,hxmin2:hxmax2,hxmin3:hxmax3))
        end if
      end if

      numerator(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3)=numerator+(dp&
         (ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3)-dm(ixMlo1:ixMhi1,&
         ixMlo2:ixMhi2,ixMlo3:ixMhi3))**2.0d0

      denominator(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3)&
         =denominator + (dabs(dp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
         ixMlo3:ixMhi3)) + dabs(dm(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
         ixMlo3:ixMhi3)) + amr_wavefilter(level)*dref(ixMlo1:ixMhi1,&
         ixMlo2:ixMhi2,ixMlo3:ixMhi3))**2.0d0

   end do
   error=error+wflags(iiflag)*dsqrt(numerator/max(denominator,epsilon))
end do

refineflag=.false.
coarsenflag=.false.

do ix3=ixMlo3,ixMhi3
do ix2=ixMlo2,ixMhi2
do ix1=ixMlo1,ixMhi1
   if (error(ix1,ix2,ix3) >= tol(level)) then
      refineflag(ix1,ix2,ix3) = .true.
   else if (error(ix1,ix2,ix3) <= tolratio(level)*tol(level)) then
      coarsenflag(ix1,ix2,ix3) = .true.
   end if
end do
end do
end do

if (any(refineflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3))&
   .and.level<mxnest) refine(igrid,mype)=.true.
if (all(coarsenflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3))&
   .and.level>1) coarsen(igrid,mype)=.true.

end subroutine lohner_orig_grid
!=============================================================================
subroutine compare1_grid(igrid,wold,w)
use mod_forest, only: coarsen, refine
use mod_amrvacdef

integer, intent(in) :: igrid
double precision, intent(in) :: wold(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
   ixGlo3:ixGhi3,1:nw), w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nw)

integer :: ix1,ix2,ix3, iiflag, iflag, level
double precision :: epsilon
double precision :: average, error
double precision :: averages(nflag_)
logical, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3) :: refineflag,&
    coarsenflag
!-----------------------------------------------------------------------------
! identify the points to be flagged in two steps:
!  step I: compare w_n-1 with w_n solution, store flags in auxiliary
!  step II: transfer flags from auxiliary to refine and coarsen

epsilon=1.0d-6

refineflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3) = .false.
coarsenflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3) = .false.
level=node(plevel_,igrid)

do ix3=ixMlo3,ixMhi3 
do ix2=ixMlo2,ixMhi2 
do ix1=ixMlo1,ixMhi1 
   average=zero
   error=zero
   do iiflag=1,flags(nflag_); iflag=flags(iiflag);
      averages(iflag) = w(ix1,ix2,ix3,iflag)-wold(ix1,ix2,ix3,iflag)
      average=average+wflags(iiflag)*abs(averages(iflag))
      if (abs(wold(ix1,ix2,ix3,iflag))<smalldouble)then
         error=error+wflags(iiflag)* abs(averages(iflag))/(abs(wold(ix1,ix2,&
            ix3,iflag))+epsilon)
      else
         error=error+wflags(iiflag)* abs(averages(iflag))/(abs(wold(ix1,ix2,&
            ix3,iflag)))
      end if
   end do
   if (error >= tol(level)) then
      refineflag(ix1,ix2,ix3) = .true.
   else if (error <= tolratio(level)*tol(level)) then
      coarsenflag(ix1,ix2,ix3) = .true.
   end if
end do
end do
end do

if (any(refineflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3))) then
   if (level<mxnest) refine(igrid,mype)=.true.
end if
if (time_advance) then
   if (all(coarsenflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3))&
      .and.level>1) coarsen(igrid,mype)=.true.
end if

end subroutine compare1_grid
!=============================================================================
subroutine createCoarse(ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,&
   ixCoGmax3)
use mod_amrvacdef

integer, intent(in) :: ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,&
   ixCoGmax3

integer :: iigrid, igrid
integer :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3
!-----------------------------------------------------------------------------
ixGmin1=ixCoGmin1;ixGmin2=ixCoGmin2;ixGmin3=ixCoGmin3;ixGmax1=ixCoGmax1
ixGmax2=ixCoGmax2;ixGmax3=ixCoGmax3;
!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
do iigrid=1,igridstail; igrid=igrids(iigrid);
   call createCoarse_grid(igrid,psCoarse(igrid),pxCoarse(igrid),ixGmin1,&
      ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,psold(igrid),px(igrid)%x)
end do
!$OMP END PARALLEL DO

end subroutine createCoarse
!=============================================================================
subroutine createCoarse_grid(igrid,sCo,pxCo,ixCoGmin1,ixCoGmin2,ixCoGmin3,&
   ixCoGmax1,ixCoGmax2,ixCoGmax3,sold,xold)

use mod_amrvacdef

integer, intent(in) :: igrid, ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,&
   ixCoGmax2,ixCoGmax3
type(state)         :: sCo, sold
double precision    :: xold(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:ndim)
type(xalloc)        :: pxCo
! .. local ..
integer             :: ixCoMmin1,ixCoMmin2,ixCoMmin3,ixCoMmax1,ixCoMmax2,&
   ixCoMmax3
!-----------------------------------------------------------------------------
dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
dxlevel(3)=rnode(rpdx3_,igrid);

! now coarsen by 2 in every direction - conservatively
! coarse grid dimension half its size
ixCoMmin1=ixCoGmin1+dixB;ixCoMmin2=ixCoGmin2+dixB;ixCoMmin3=ixCoGmin3+dixB
ixCoMmax1=ixCoGmax1-dixB;ixCoMmax2=ixCoGmax2-dixB;ixCoMmax3=ixCoGmax3-dixB;
call coarsen_grid(sold,xold,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixMlo1,&
   ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,sCo,pxCo%x,ixCoGmin1,ixCoGmin2,&
   ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3,ixCoMmin1,ixCoMmin2,ixCoMmin3,&
   ixCoMmax1,ixCoMmax2,ixCoMmax3, pgeo(igrid),pgeoCoarse(igrid),&
   coarsenprimitive,.true.)

end subroutine createCoarse_grid
!=============================================================================
subroutine advectCoarse(ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,&
   ixCoGmax3,factor)
use mod_amrvacdef

integer :: ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3
double precision, intent(in) :: factor

integer :: iigrid, igrid
!-----------------------------------------------------------------------------
!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid)
do iigrid=1,igridstail; igrid=igrids(iigrid);
   call advectCoarse_grid(igrid,psCoarse(igrid),factor)
end do
!$OMP END PARALLEL DO

end subroutine advectCoarse
!=============================================================================
subroutine advectCoarse_grid(igrid,psCo,factor)

use mod_amrvacdef

integer, intent(in) :: igrid
double precision, intent(in) :: factor
type(state)         :: psCo

double precision :: qdt, dx1,dx2,dx3
!double precision, dimension(:^D&,:), allocatable :: wCo1
type(state)      :: sCo1
double precision :: fC(psCo%w%ixGmin1:psCo%w%ixGmax1,&
   psCo%w%ixGmin2:psCo%w%ixGmax2,psCo%w%ixGmin3:psCo%w%ixGmax3,1:nwflux,&
   1:ndim)
integer :: level
!-----------------------------------------------------------------------------
dxlevel(1)=two*rnode(rpdx1_,igrid);dxlevel(2)=two*rnode(rpdx2_,igrid)
dxlevel(3)=two*rnode(rpdx3_,igrid);
dx1=dxlevel(1);dx2=dxlevel(2);dx3=dxlevel(3);

! here we integrate on the coarse grid
call copy_state(psCo,sCo1)
!allocate(wCo1(ixCoG^S,1:nw))
!wCo1(ixCoG^S,1:nwflux)=pwCo%w(ixCoG^S,1:nwflux)



! 1st order scheme: do coarse time step of
!    size 2*dt starting from t_n-1 solution in pwCoarse
!    to arrive at t_n+1 (n-index from normal uncoarsened grid)
! result in pwCoarse: coarse solution at t_n+1

if (.not.slab) mygeo => pgeoCoarse(igrid)
if(covariant)myM => mygeo%m

qdt=factor*dt_grid(igrid)
level=node(plevel_,igrid)
call advect1_grid(typelow1(level),qdt,psCo%w%ixGmin1,psCo%w%ixGmin2,&
   psCo%w%ixGmin3,psCo%w%ixGmax1,psCo%w%ixGmax2,psCo%w%ixGmax3,1,ndim,t,sCo1,&
   t, psCo,sCo1,fC,dx1,dx2,dx3,pxCoarse(igrid)%x)

!deallocate(wCo1)

call dealloc_state(sCo1)

end subroutine advectCoarse_grid
!=============================================================================
subroutine errest1_grid(igrid,s)

use mod_amrvacdef

integer, intent(in) :: igrid
type(state), intent(in) :: s

integer :: level, ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3
double precision :: dx1,dx2,dx3, qdt
double precision :: fC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux,&
   1:ndim) !, wFi(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nw)
type(state)      :: sFi
!-----------------------------------------------------------------------------
level=node(plevel_,igrid)

dx1=dx(1,level);dx2=dx(2,level);dx3=dx(3,level);
dxlevel(1)=dx1;dxlevel(2)=dx2;dxlevel(3)=dx3;
call copy_state(s,sFi)
!wFi(ixG^T,1:nwflux)=w(ixG^T,1:nwflux)

if (.not.skipfinestep) then
   if (.not.slab) mygeo => pgeo(igrid)
   if(covariant)myM => mygeo%m

   qdt=dt_grid(igrid)
   call advect1_grid(typelow1(level),qdt,s%w%ixGmin1,s%w%ixGmin2,s%w%ixGmin3,&
      s%w%ixGmax1,s%w%ixGmax2,s%w%ixGmax3,1,ndim,t+qdt,s, t+qdt,sFi,s,fC,dx1,&
      dx2,dx3,px(igrid)%x)
end if

ixCoGmin1=1;ixCoGmin2=1;ixCoGmin3=1;
ixCoGmax1=ixGhi1/2+dixB;ixCoGmax2=ixGhi2/2+dixB;ixCoGmax3=ixGhi3/2+dixB;

call flagbadpoints(sFi%w%w,pwCoarse(igrid)%w,ixCoGmin1,ixCoGmin2,ixCoGmin3,&
   ixCoGmax1,ixCoGmax2,ixCoGmax3,igrid,level)

call dealloc_state(sFi)

end subroutine errest1_grid
!=============================================================================
subroutine flagbadpoints(w,wCo,ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,&
   ixCoGmax2,ixCoGmax3,igrid,level)

! compare error between coarse and fine solution in wCo, w 
! We base the comparison on the physical field selected by the index flag_ 
!
! on entry:
!  w:       normal time integration 
!  wCo:     time integration on coarsened grid (2*dx)

use mod_forest, only: coarsen, refine
use mod_amrvacdef

integer, intent(in)         :: igrid, ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,&
   ixCoGmax2,ixCoGmax3, level
double precision,intent(in) :: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
   nw), wCo(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,nw)
double precision            :: specialvar(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
   ixGlo3:ixGhi3), specialvarCo(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,&
   ixCoGmin3:ixCoGmax3)

logical :: needgetaux
integer :: iCo1,iCo2,iCo3, iFi1,iFi2,iFi3, ixCoMmin1,ixCoMmin2,ixCoMmin3,&
   ixCoMmax1,ixCoMmax2,ixCoMmax3, iiflag, iflag
double precision :: average, error
double precision :: averages(nflag_)
logical, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3) :: refineflag,&
    coarsenflag
!-----------------------------------------------------------------------------
ixCoMmin1=ixCoGmin1+dixB;ixCoMmin2=ixCoGmin2+dixB;ixCoMmin3=ixCoGmin3+dixB
ixCoMmax1=ixCoGmax1-dixB;ixCoMmax2=ixCoGmax2-dixB;ixCoMmax3=ixCoGmax3-dixB;

needgetaux=.false.
do iiflag=1,flags(nflag_); iflag=flags(iiflag);
  if (iflag>nwflux) needgetaux=.true.
end do
if (nwaux>0.and.needgetaux) then
   saveigrid=igrid
   if(.not.slab)mygeo=>pgeo(igrid)
      if(covariant)myM => mygeo%m
   call getaux(.true.,w,px(igrid)%x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
      ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,'flagbadpoints')
   call getaux(.true.,wCo,pxCoarse(igrid)%x,ixCoGmin1,ixCoGmin2,ixCoGmin3,&
      ixCoGmax1,ixCoGmax2,ixCoGmax3,ixCoMmin1,ixCoMmin2,ixCoMmin3,ixCoMmax1,&
      ixCoMmax2,ixCoMmax3,'flagbadpointsCo')
end if

! identify the points to be flagged in two steps (needed!):
!  step I: compare coarse with fine solution, store flags in fine auxiliary
!  step II: transfer flags from auxiliary to refine and coarsen

refineflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3) = .false.
coarsenflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3) = .false.


do iiflag=1,flags(nflag_); iflag=flags(iiflag);
   if (iflag>nw) then
      call specialvarforerrest(ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,&
         ixCoGmax2,ixCoGmax3,ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,&
         ixCoGmax2,ixCoGmax3,iflag,wCo,specialvarCo)
      call specialvarforerrest(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
         ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,iflag,w,specialvar)
   end if
end do

iFi3 = ixMlo3
do iCo3 = ixCoMmin3,ixCoMmax3 
iFi2 = ixMlo2
do iCo2 = ixCoMmin2,ixCoMmax2 
iFi1 = ixMlo1
do iCo1 = ixCoMmin1,ixCoMmax1 
   average=zero
   error=zero
   do iiflag=1,flags(nflag_); iflag=flags(iiflag);
      if (slab) then
         if (iflag<=nw) averages(iflag)=sum(w(iFi1:iFi1+1,iFi2:iFi2+1,&
            iFi3:iFi3+1,iflag))/two**ndim
         if (iflag>nw)  averages(iflag)=sum(specialvar(iFi1:iFi1+1,&
            iFi2:iFi2+1,iFi3:iFi3+1))/two**ndim
      else
         if (iflag<=nw) averages(iflag)=sum(pgeo(igrid)%dvolume(iFi1:iFi1+1,&
            iFi2:iFi2+1,iFi3:iFi3+1) *w(iFi1:iFi1+1,iFi2:iFi2+1,iFi3:iFi3+1,&
            iflag))/pgeoCoarse(igrid)%dvolume(iCo1,iCo2,iCo3)
         if (iflag>nw)  averages(iflag)=sum(pgeo(igrid)%dvolume(iFi1:iFi1+1,&
            iFi2:iFi2+1,iFi3:iFi3+1) *specialvar(iFi1:iFi1+1,iFi2:iFi2+1,&
            iFi3:iFi3+1))/pgeoCoarse(igrid)%dvolume(iCo1,iCo2,iCo3)
      end if
      average=average+wflags(iiflag)*abs(averages(iflag))
      if (iflag<=nw) error=error+wflags(iiflag)*abs(averages(iflag)-wCo(iCo1,&
         iCo2,iCo3,iflag))
      if (iflag> nw) error=error+wflags(iiflag)*abs(averages&
         (iflag)-specialvarCo(iCo1,iCo2,iCo3))
   end do
   if (abs(average)>smalldouble) then
      error=error/average
   else
      write(unitterm,*)'Warning from flagbadpoints: zero average:',average
      write(unitterm,*)'   wCo(iCo^D,1:nw):',wCo(iCo1,iCo2,iCo3,1:nw),&
         ' indices:',iCo1,iCo2,iCo3
      write(unitterm,*)'On grid:',igrid,' at level ',level
      write(unitterm,*)'   and grid indices : ',node(pig1_,igrid),node(pig2_,&
         igrid),node(pig3_,igrid)
      write(unitterm,*)'cell indices : ',iCo1,iCo2,iCo3
      call mpistop("")
   end if
   if (error >= tol(level)) then
      refineflag(iFi1:iFi1+1,iFi2:iFi2+1,iFi3:iFi3+1) = .true.
   else if (error <= tolratio(level)*tol(level)) then
      coarsenflag(iFi1:iFi1+1,iFi2:iFi2+1,iFi3:iFi3+1) = .true.
   end if
   iFi1 = iFi1+2
end do
   iFi2 = iFi2+2
end do
   iFi3 = iFi3+2
end do

iFi3 = ixMlo3
do iCo3 = ixCoMmin3,ixCoMmax3 
iFi2 = ixMlo2
do iCo2 = ixCoMmin2,ixCoMmax2 
iFi1 = ixMlo1
do iCo1 = ixCoMmin1,ixCoMmax1 
   if (error >= tol(level)) then
      refineflag(iFi1:iFi1+1,iFi2:iFi2+1,iFi3:iFi3+1) = .true.
   else if (error <= tolratio(level)*tol(level)) then
      coarsenflag(iFi1:iFi1+1,iFi2:iFi2+1,iFi3:iFi3+1) = .true.
   end if
   iFi1 = iFi1+2
end do
   iFi2 = iFi2+2
end do
   iFi3 = iFi3+2
end do

if (any(refineflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3))) then
   if (level<mxnest) refine(igrid,mype)=.true.
end if
if (time_advance) then
   if (all(coarsenflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3))&
      .and.level>1) coarsen(igrid,mype)=.true.
end if

end subroutine flagbadpoints
!=============================================================================
subroutine forcedrefine_grid(igrid,w)
use mod_forest, only: coarsen, refine, buffer
use mod_amrvacdef

integer, intent(in) :: igrid
double precision, intent(in) :: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
   nw)

integer :: level
integer :: my_refine, my_coarsen
double precision :: qt
!-----------------------------------------------------------------------------
level=node(plevel_,igrid)

! initialize to 0
my_refine   = 0
my_coarsen  = 0

if (time_advance) then
   qt=t+dt
else
   if (errorestimate==1.or.errorestimate==2) then
      qt=t+dt
   else
      qt=t
   end if
end if
   
call specialrefine_grid(igrid,level,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
   ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,qt,w,px(igrid)%x, my_refine,&
   my_coarsen)

if (my_coarsen==1) then
   if (level>1) then
      refine(igrid,mype)=.false.
      coarsen(igrid,mype)=.true.
   else
      refine(igrid,mype)=.false.
      coarsen(igrid,mype)=.false.
   end if
endif

if (my_coarsen==-1)then
   coarsen(igrid,mype)=.false.
end if

if (my_refine==1) then
   if (level<mxnest) then
      refine(igrid,mype)=.true.
      coarsen(igrid,mype)=.false.
   else
      refine(igrid,mype)=.false.
      coarsen(igrid,mype)=.false.
   end if
end if

if (my_refine==-1) then
  refine(igrid,mype)=.false.
end if

end subroutine forcedrefine_grid
!=============================================================================
subroutine forcedrefine_grid_io(igrid,w)
use mod_forest, only: coarsen, refine
use mod_amrvacdef

integer, intent(in)          :: igrid
double precision, intent(in) :: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
   nw)

integer                   :: level, my_levmin, my_levmax
logical, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3) :: refineflag
!-----------------------------------------------------------------------------
level=node(plevel_,igrid)

if (level_io > 0) then
   my_levmin = level_io
   my_levmax = level_io
else
   my_levmin = max(1,level_io_min)
   my_levmax = min(mxnest,level_io_max)
end if


if (level>my_levmax) then
      refine(igrid,mype)=.false.
      coarsen(igrid,mype)=.true.
elseif (level<my_levmin) then
      refine(igrid,mype)=.true.
      coarsen(igrid,mype)=.false.
end if

if (level==my_levmin .or. level==my_levmax) then
  refine(igrid,mype)=.false.
  coarsen(igrid,mype)=.false.
end if


if(refine(igrid,mype).and.level>=mxnest)refine(igrid,mype)=.false.
if(coarsen(igrid,mype).and.level<=1)coarsen(igrid,mype)=.false.

end subroutine forcedrefine_grid_io
!=============================================================================
subroutine refinebuffer(igrid,refineflag)
use mod_forest, only: refine, buffer
use mod_amrvacdef

integer, intent(in) :: igrid
logical, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
    intent(in) :: refineflag

integer :: ishiftbuf1,ishiftbuf2,ishiftbuf3, i1,i2,i3, ixmin1,ixmin2,ixmin3,&
   ixmax1,ixmax2,ixmax3, ineighbor, ipe_neighbor, level
!-----------------------------------------------------------------------------
ishiftbuf1=ixMhi1-ixMlo1-nbufferx1+1;ishiftbuf2=ixMhi2-ixMlo2-nbufferx2+1
ishiftbuf3=ixMhi3-ixMlo3-nbufferx3+1;
do i3=-1,1
do i2=-1,1
do i1=-1,1
   ixmin1=max(ixMlo1,ixMlo1+i1*ishiftbuf1)
   ixmin2=max(ixMlo2,ixMlo2+i2*ishiftbuf2)
   ixmin3=max(ixMlo3,ixMlo3+i3*ishiftbuf3);
   ixmax1=min(ixMhi1,ixMhi1+i1*ishiftbuf1)
   ixmax2=min(ixMhi2,ixMhi2+i2*ishiftbuf2)
   ixmax3=min(ixMhi3,ixMhi3+i3*ishiftbuf3);
   if (ixmax1<ixmin1.or.ixmax2<ixmin2.or.ixmax3<ixmin3) cycle
   if (any(refineflag(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3))) then
      select case (neighbor_type(i1,i2,i3,igrid))
      case (2)
         ineighbor=neighbor(1,i1,i2,i3,igrid)
         ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
         if (.not.refine(ineighbor,ipe_neighbor)) then
            buffer(ineighbor,ipe_neighbor)=.true.
            refine(ineighbor,ipe_neighbor)=.true.
         end if
      case (3)
         level=node(plevel_,igrid)
         if (level<mxnest) then
            ineighbor=neighbor(1,i1,i2,i3,igrid)
            ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
            if (.not.refine(ineighbor,ipe_neighbor)) then
               buffer(ineighbor,ipe_neighbor)=.true.
               refine(ineighbor,ipe_neighbor)=.true.
            end if
         end if
      end select
   end if
end do
end do
end do

end subroutine refinebuffer
!=============================================================================
