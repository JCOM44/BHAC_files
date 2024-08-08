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
subroutine bc_phys(iside,idims,time,s,ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,&
   ixBmax3)

use mod_amrvacdef

integer, intent(in)          :: iside, idims, ixBmin1,ixBmin2,ixBmin3,ixBmax1,&
   ixBmax2,ixBmax3
double precision, intent(in) :: time
type(state), intent(inout)   :: s

integer :: iw, iB, ix1,ix2,ix3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3, ixDmin1,ixDmin2,ixDmin3,ixDmax1,ixDmax2,ixDmax3, ixIpmin1,&
   ixIpmin2,ixIpmin3,ixIpmax1,ixIpmax2,ixIpmax3, ixGmin1,ixGmin2,ixGmin3,&
   ixGmax1,ixGmax2,ixGmax3

integer :: idir, is
integer :: ixIsmin1,ixIsmin2,ixIsmin3,ixIsmax1,ixIsmax2,ixIsmax3,hxImin1,&
   hxImin2,hxImin3,hxImax1,hxImax2,hxImax3,jxImin1,jxImin2,jxImin3,jxImax1,&
   jxImax2,jxImax3
double precision :: Q(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
   Qp(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
logical, save    :: zerodiv=.true.
double precision, dimension(:,:,:), pointer :: surf

double precision, dimension(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
   1:nw)    :: wsave
integer, parameter :: dixPolefix=2, dixHotaka=1
logical            :: is_coarse(ndim,2)
!-----------------------------------------------------------------------------
associate(x=>s%x%x,w=>s%w%w ,ws=>s%ws%w)
ixGmin1=s%w%ixGmin1;ixGmin2=s%w%ixGmin2;ixGmin3=s%w%ixGmin3
ixGmax1=s%w%ixGmax1;ixGmax2=s%w%ixGmax2;ixGmax3=s%w%ixGmax3;

! Test if we have a coarse neighbor:
call is_neighbor_coarse(s,is_coarse)

select case (idims)
case (1)
if (iside==2) then

   
   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! maximal boundary
   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      iB=ismax1
      ixImin1=ixBmax1+1-dixB;ixImin2=ixBmin2;ixImin3=ixBmin3;
      ixImax1=ixBmax1;ixImax2=ixBmax2;ixImax3=ixBmax3;

      !!!!!!! Make tangential ranges greater for staggered components
      !!!!!!! Add fill the normal component

      ! Apply on primitive variables?
      if (primitiveB(iside,idims)) then
         ixDmin1=ixBmin1+dixB;ixDmin2=ixBmin2;ixDmin3=ixBmin3;
         ixDmax1=ixBmax1-dixB;ixDmax2=ixBmax2;ixDmax3=ixBmax3;
         wsave(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) &
            = w(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) !make a backup of the state in domain since p2c/c2p might lead to slight changes in w
         call primitive(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
            ixDmin1,ixDmin2,ixDmin3,ixDmax1,ixDmax2,ixDmax3,w,x)
      end if
      
      ! cont/symm/asymm types
      do iw=1,nwflux+nwaux
 ! If the variable has no staggered counterpart
         if ((iws(iw).lt.1).or.(iws(iw).gt.nws)) then

         select case (typeB(iw,iB))
         case ("symm")
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               = w(ixImin1-1:ixImin1-dixB:-1,ixImin2:ixImax2,ixImin3:ixImax3,&
               iw)
         case ("asymm")
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               =-w(ixImin1-1:ixImin1-dixB:-1,ixImin2:ixImax2,ixImin3:ixImax3,&
               iw)
         case ("cont")
            do ix1=ixImin1,ixImax1
               w(ix1,ixImin2:ixImax2,ixImin3:ixImax3,iw) = w(ixImin1-1,&
                  ixImin2:ixImax2,ixImin3:ixImax3,iw)
            end do
         case("noinflow")
            if (iw==1+1)then
              do ix1=ixImin1,ixImax1
                  w(ix1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
                     = max(w(ixImin1-1,ixImin2:ixImax2,ixImin3:ixImax3,iw),&
                     zero)
              end do
            else
              do ix1=ixImin1,ixImax1
                  w(ix1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
                     = w(ixImin1-1,ixImin2:ixImax2,ixImin3:ixImax3,iw)
              end do
            end if
         case("limitinflow")
            if (iw==1+1)then
              do ix1=ixImin1,ixImax1
                  w(ix1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
                     = max(w(ixImin1-1,ixImin2:ixImax2,ixImin3:ixImax3,iw), &
                                           w(ixImin1-1,ixImin2:ixImax2,&
                                              ixImin3:ixImax3,iw)*ratebdflux)
              end do
            else
              do ix1=ixImin1,ixImax1
                  w(ix1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
                     = w(ixImin1-1,ixImin2:ixImax2,ixImin3:ixImax3,iw)
              end do
           end if
        case ("polefix")
           ! First overwrite cells in domain
           ixIpmin1=ixImin1-dixPolefix;ixIpmin2=ixImin2;ixIpmin3=ixImin3;
           ixIpmax1=ixImin1-1;ixIpmax2=ixImax2;ixIpmax3=ixImax3;
           do ix1=ixIpmin1,ixIpmax1
              w(ix1,ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,iw) &
                 = w(ixIpmin1-1,ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,iw)
           end do
           ! Now apply symm boundary
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               = w(ixImin1-1:ixImin1-dixB:-1,ixImin2:ixImax2,ixImin3:ixImax3,&
               iw)
        case ("apolefix")
           ! First interpolate cells in domain down to zero at pole
           ixIpmin1=ixImin1-dixPolefix;ixIpmin2=ixImin2;ixIpmin3=ixImin3;
           ixIpmax1=ixImax1-1;ixIpmax2=ixImax2;ixIpmax3=ixImax3;
           do ix1=ixIpmin1,ixIpmax1
              w(ix1,ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,iw) &
                 = w(ixIpmin1-1,ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,iw) &
                   * (xprobmax1-x(ix1,ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,&
                      1)) &
                   / (xprobmax1-x(ixIpmin1-1,ixIpmin2:ixIpmax2,&
                      ixIpmin3:ixIpmax3,1))
           end do
        case ("hotaka")
           ! First overwrite cells in domain with zero
           ixIpmin1=ixImin1-dixHotaka;ixIpmin2=ixImin2;ixIpmin3=ixImin3;
           ixIpmax1=ixImin1-1;ixIpmax2=ixImax2;ixIpmax3=ixImax3;
           do ix1=ixIpmin1,ixIpmax1
              w(ix1,ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,iw) = zero
           end do
           ! Now apply symm boundary
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               = w(ixImin1-1:ixImin1-dixB:-1,ixImin2:ixImax2,ixImin3:ixImax3,&
               iw)
        case ("ahotaka")
            ! First overwrite cells in domain with zero
           ixIpmin1=ixImin1-dixHotaka;ixIpmin2=ixImin2;ixIpmin3=ixImin3;
           ixIpmax1=ixImax1-1;ixIpmax2=ixImax2;ixIpmax3=ixImax3;
           do ix1=ixIpmin1,ixIpmax1
              w(ix1,ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,iw) = zero
           end do
           ! Now apply asymm boundary
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               =-w(ixImin1-1:ixImin1-dixB:-1,ixImin2:ixImax2,ixImin3:ixImax3,&
               iw)
         case ("special")
            ! skip it here, do AFTER all normal type boundaries are set
         case ("aperiodic")
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               = - w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw)
         case ("periodic")
!            call mpistop("periodic bc info should come from neighbors")
         case default
            write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
               "for variable iw=",iw," and side iB=",iB
         end select


      else !If there is a staggered counterpart
 !At this stage, extrapolation is applied only to the tangential components
        idir=mod(iws(iw)-1,ndir)+1
        if (idir.ne.1) then
        ixIsmax1=ixImax1;ixIsmax2=ixImax2;ixIsmax3=ixImax3;
        ixIsmin1=ixImin1-kr(1,idir);ixIsmin2=ixImin2-kr(2,idir)
        ixIsmin3=ixImin3-kr(3,idir);
        select case(typeB(iw,iB))
        case ("symm")
          ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iws(iw)) &
             = ws(ixIsmin1-1:ixIsmin1-dixB:-1,ixIsmin2:ixIsmax2,&
             ixIsmin3:ixIsmax3,iws(iw))
        case ("asymm")
          ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iws(iw)) &
             =-ws(ixIsmin1-1:ixIsmin1-dixB:-1,ixIsmin2:ixIsmax2,&
             ixIsmin3:ixIsmax3,iws(iw))
        case ("cont")
          do ix1=ixIsmin1,ixIsmax1
             ws(ix1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iws(iw)) &
                = ws(ixIsmin1-1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iws(iw))
          end do
        case ("periodic")
        case ("zero")
           ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iws(iw)) &
              = 0.0d0
        case ("specialT")
           call specialbound_usr(time,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
              ixGmax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,-iB,s)
        case ("special")
           ! skip it here, do AFTER all normal type boundaries are set
        case default
          write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
             "for variable iw=",iw," and side iB=",iB
        end select

                  
      ! Special treatment when we have coarse neighbors in tangential dir:
      do is=1,2
         if (is_coarse(idir,is)) then

            ixIsmin1=ixImin1-kr(1,idir);ixIsmin2=ixImin2-kr(2,idir)
            ixIsmin3=ixImin3-kr(3,idir);
            ixIsmax1=ixImax1;ixIsmax2=ixImax2;ixIsmax3=ixImax3;

            if (is.eq.1) then ! coarse neighbor at beginning of block

               select case(idir)
                  
                  case (1)
                     ixIsmax1=min(ixGlo1+dixB-1,ixIsmax1);
                     surf=>s%geo%surfaceC1
                  

                  case (2)
                     ixIsmax2=min(ixGlo2+dixB-1,ixIsmax2);
                     surf=>s%geo%surfaceC2
                  

                  case (3)
                     ixIsmax3=min(ixGlo3+dixB-1,ixIsmax3);
                     surf=>s%geo%surfaceC3
                  ! I am hating LASY just right now.
               end select

            else ! coarse neighbor at end of block

               select case(idir)
                  
                  case (1)
                     ixIsmin1=max(ixGhi1-dixB,ixIsmin1);
                     surf=>s%geo%surfaceC1
                  

                  case (2)
                     ixIsmin2=max(ixGhi2-dixB,ixIsmin2);
                     surf=>s%geo%surfaceC2
                  

                  case (3)
                     ixIsmin3=max(ixGhi3-dixB,ixIsmin3);
                     surf=>s%geo%surfaceC3
                  ! I am hating LASY just right now.
               end select

            end if


            select case(typeB(iw,iB))
            case ("cont")
               do ix1=ixIsmin1,ixIsmax1
                  where ((surf(ixIsmin1-1,ixIsmin2:ixIsmax2,&
                     ixIsmin3:ixIsmax3) + surf(ixIsmin1-2,ixIsmin2:ixIsmax2,&
                     ixIsmin3:ixIsmax3)) .gt. smalldouble )
                     ws(ix1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iws(iw)) &
                        = (ws(ixIsmin1-1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,&
                        iws(iw))*surf(ixIsmin1-1,ixIsmin2:ixIsmax2,&
                        ixIsmin3:ixIsmax3) &
                          + ws(ixIsmin1-2,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,&
                             iws(iw))*surf(ixIsmin1-2,ixIsmin2:ixIsmax2,&
                             ixIsmin3:ixIsmax3) ) / (surf(ixIsmin1-1,&
                             ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3)+surf&
                             (ixIsmin1-2,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3))
                  elsewhere
                     ws(ix1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iws(iw)) &
                        = zero
                  end where
               end do

            case default
 !Only implemented cont type, symm and asymm will also work if problem respects symmetry
               ! 3D case will need to average four cells together!!!
            end select

         end if
      end do

   end if
end if

      end do

      ! Now that the tangential components are set,
      ! fill the normal components using a prescription for the divergence.
      ! This prescription is given by the typeB for the normal component.
      do iw=1,nws
        ! When using more than one staggered vector field
        idir=mod(iw-1,ndir)+1
        ! Consider only normal direction
        if (idir.eq.1) then
          ixIsmin1=ixImin1;ixIsmin2=ixImin2;ixIsmin3=ixImin3;ixIsmax1=ixImax1
          ixIsmax2=ixImax2;ixIsmax3=ixImax3;
          hxImin1=ixImin1-dixB*kr(1,1);hxImin2=ixImin2-dixB*kr(2,1)
          hxImin3=ixImin3-dixB*kr(3,1);hxImax1=ixImax1-dixB*kr(1,1)
          hxImax2=ixImax2-dixB*kr(2,1);hxImax3=ixImax3-dixB*kr(3,1);

          ! Calculate divergence and partial divergence
          Q=zero
          if (.not.zerodiv) call div_staggered(hxImin1,hxImin2,hxImin3,&
             hxImax1,hxImax2,hxImax3,s,Q(hxImin1:hxImax1,hxImin2:hxImax2,&
             hxImin3:hxImax3))
          select case(typeB(iw+b0_,iB))
          case("symm")
            ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iw)=zero
            do ix1=0,dixB-1
              call div_staggered(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                 ixImax3,s,Qp(ixImin1:ixImax1,ixImin2:ixImax2,&
                 ixImin3:ixImax3))
              ws(ixIsmin1+ix1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iw)=&
              (Q(hxImax1-ix1,hxImin2:hxImax2,hxImin3:hxImax3)*s%geo%dvolume&
                 (hxImax1-ix1,hxImin2:hxImax2,hxImin3:hxImax3)&
             -Qp(ixImin1+ix1,ixImin2:ixImax2,ixImin3:ixImax3)*s%geo%dvolume&
                (ixImin1+ix1,ixImin2:ixImax2,ixImin3:ixImax3))&
              /s%geo%surfaceC1(ixIsmin1+ix1,ixIsmin2:ixIsmax2,&
                 ixIsmin3:ixIsmax3)
            end do
          case("asymm")
            ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iw)=zero
            do ix1=0,dixB-1
              call div_staggered(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                 ixImax3,s,Qp(ixImin1:ixImax1,ixImin2:ixImax2,&
                 ixImin3:ixImax3))
              ws(ixIsmin1+ix1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iw)=&
             (-Q(hxImax1-ix1,hxImin2:hxImax2,hxImin3:hxImax3)*s%geo%dvolume&
                (hxImax1-ix1,hxImin2:hxImax2,hxImin3:hxImax3)&
             -Qp(ixImin1+ix1,ixImin2:ixImax2,ixImin3:ixImax3)*s%geo%dvolume&
                (ixImin1+ix1,ixImin2:ixImax2,ixImin3:ixImax3))&
              /s%geo%surfaceC1(ixIsmin1+ix1,ixIsmin2:ixIsmax2,&
                 ixIsmin3:ixIsmax3)
            end do
          case("cont", "zero")
            ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iw)=zero
            do ix1=0,dixB-1
              call div_staggered(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                 ixImax3,s,Qp(ixImin1:ixImax1,ixImin2:ixImax2,&
                 ixImin3:ixImax3))
              ws(ixIsmin1+ix1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iw)=&
              (Q(hxImax1,hxImin2:hxImax2,hxImin3:hxImax3)*s%geo%dvolume&
                 (hxImax1,hxImin2:hxImax2,hxImin3:hxImax3)&
             -Qp(ixImin1+ix1,ixImin2:ixImax2,ixImin3:ixImax3)*s%geo%dvolume&
                (ixImin1+ix1,ixImin2:ixImax2,ixImin3:ixImax3))&
              /s%geo%surfaceC1(ixIsmin1+ix1,ixIsmin2:ixIsmax2,&
                 ixIsmin3:ixIsmax3)
            end do
          case("periodic")
          end select
        end if
      end do

      ! Fill cell averages
      call faces2centers(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,s)



      ! Apply on primitive variables?
      if (primitiveB(iside,idims)) then
         w(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) &
            = wsave(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) !restore from backup of the state in domain since p2c/c2p might lead to slight changes in w
         call conserve(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
            ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,w,x,patchfalse)
      end if
      
   else
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! minimal boundary
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      iB=ismin1
      ixImin1=ixBmin1;ixImin2=ixBmin2;ixImin3=ixBmin3;
      ixImax1=ixBmin1-1+dixB;ixImax2=ixBmax2;ixImax3=ixBmax3;

      ! Apply on primitive variables?
      if (primitiveB(iside,idims)) then         
         ixDmin1=ixBmin1+dixB;ixDmin2=ixBmin2;ixDmin3=ixBmin3;
         ixDmax1=ixBmax1-dixB;ixDmax2=ixBmax2;ixDmax3=ixBmax3;
         wsave(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) &
            = w(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) !make a backup of the state in domain since p2c/c2p might lead to slight changes in w
         call primitive(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
            ixDmin1,ixDmin2,ixDmin3,ixDmax1,ixDmax2,ixDmax3,w,x)
      end if
 
      ! cont/symm/asymm types
      do iw=1,nwflux+nwaux
 ! If the variable has no staggered counterpart
if ((iws(iw).lt.1).or.(iws(iw).gt.nws)) then
   
   select case (typeB(iw,iB))
   case ("symm")
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         = w(ixImax1+dixB:ixImax1+1:-1,ixImin2:ixImax2,ixImin3:ixImax3,iw)
   case ("asymm")
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         =-w(ixImax1+dixB:ixImax1+1:-1,ixImin2:ixImax2,ixImin3:ixImax3,iw)
   case ("cont")
      do ix1=ixImin1,ixImax1
         w(ix1,ixImin2:ixImax2,ixImin3:ixImax3,iw) = w(ixImax1+1,&
            ixImin2:ixImax2,ixImin3:ixImax3,iw)
      end do
   case("noinflow")
      if (iw==1+1)then
         do ix1=ixImin1,ixImax1
            w(ix1,ixImin2:ixImax2,ixImin3:ixImax3,iw) = min(w(ixImax1+1,&
               ixImin2:ixImax2,ixImin3:ixImax3,iw),zero)
         end do
      else
         do ix1=ixImin1,ixImax1
            w(ix1,ixImin2:ixImax2,ixImin3:ixImax3,iw) = w(ixImax1+1,&
               ixImin2:ixImax2,ixImin3:ixImax3,iw)
         end do
      end if
   case("limitinflow")
      if (iw==1+1)then
         do ix1=ixImin1,ixImax1
            w(ix1,ixImin2:ixImax2,ixImin3:ixImax3,iw) = min(w(ixImax1+1,&
               ixImin2:ixImax2,ixImin3:ixImax3,iw), &
                 w(ixImax1+1,ixImin2:ixImax2,ixImin3:ixImax3,iw)*ratebdflux)
         end do
      else
         do ix1=ixImin1,ixImax1
            w(ix1,ixImin2:ixImax2,ixImin3:ixImax3,iw) = w(ixImax1+1,&
               ixImin2:ixImax2,ixImin3:ixImax3,iw)
         end do
      end if
   case ("polefix")
      ! First overwrite cells in domain
      ixIpmin1=ixImax1+1;ixIpmin2=ixImin2;ixIpmin3=ixImin3;
      ixIpmax1=ixImax1+dixPolefix;ixIpmax2=ixImax2;ixIpmax3=ixImax3;
      do ix1=ixIpmin1,ixIpmax1
         w(ix1,ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,iw) = w(ixIpmax1+1,&
            ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,iw)
      end do
      ! Now apply symm boundary
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         = w(ixImax1+dixB:ixImax1+1:-1,ixImin2:ixImax2,ixImin3:ixImax3,iw)
   case ("apolefix")
      ! First interpolate cells in domain down to zero at pole
      ixIpmin1=ixImax1+1;ixIpmin2=ixImin2;ixIpmin3=ixImin3;
      ixIpmax1=ixImax1+dixPolefix;ixIpmax2=ixImax2;ixIpmax3=ixImax3;
      do ix1=ixIpmin1,ixIpmax1
         w(ix1,ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,iw) = w(ixIpmax1+1,&
            ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,iw) &
              * (x(ix1,ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,1))&
                 /x(ixIpmax1+1,ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,1)
      end do
      ! Now apply asymm boundary
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         =-w(ixImax1+dixB:ixImax1+1:-1,ixImin2:ixImax2,ixImin3:ixImax3,iw)
   case ("hotaka")
      ! First overwrite cells in domain with zero
      ixIpmin1=ixImax1+1;ixIpmin2=ixImin2;ixIpmin3=ixImin3;
      ixIpmax1=ixImax1+dixHotaka;ixIpmax2=ixImax2;ixIpmax3=ixImax3;
      do ix1=ixIpmin1,ixIpmax1
         w(ix1,ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,iw) = zero
      end do
      ! Now apply symm boundary
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         = w(ixImax1+dixB:ixImax1+1:-1,ixImin2:ixImax2,ixImin3:ixImax3,iw)
   case ("ahotaka")
      ! First overwrite cells in domain with zero
      ixIpmin1=ixImax1+1;ixIpmin2=ixImin2;ixIpmin3=ixImin3;
      ixIpmax1=ixImax1+dixHotaka;ixIpmax2=ixImax2;ixIpmax3=ixImax3;
      do ix1=ixIpmin1,ixIpmax1
         w(ix1,ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,iw) = zero
      end do
      ! Now apply asymm boundary
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         =-w(ixImax1+dixB:ixImax1+1:-1,ixImin2:ixImax2,ixImin3:ixImax3,iw)
   case ("special")
      ! skip it here, do AFTER all normal type boundaries are set
   case ("aperiodic")
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         = - w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw)
   case ("periodic")
      !            call mpistop("periodic bc info should come from neighbors")
   case default
      write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
           "for variable iw=",iw," and side iB=",iB
   end select


else !If there is a staggered counterpart
   ! At this stage, extrapolation is applied only to the tangential components
   idir=mod(iws(iw)-1,ndir)+1 ! vector direction
   if (idir.ne.1) then
      ixIsmax1=ixImax1;ixIsmax2=ixImax2;ixIsmax3=ixImax3;
      ixIsmin1=ixImin1-kr(1,idir);ixIsmin2=ixImin2-kr(2,idir)
      ixIsmin3=ixImin3-kr(3,idir);
      select case(typeB(iw,iB))
      case ("symm")
         ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iws(iw)) &
            = ws(ixIsmax1+dixB:ixIsmax1+1:-1,ixIsmin2:ixIsmax2,&
            ixIsmin3:ixIsmax3,iws(iw))
      case ("asymm")
         ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iws(iw)) &
            =-ws(ixIsmax1+dixB:ixIsmax1+1:-1,ixIsmin2:ixIsmax2,&
            ixIsmin3:ixIsmax3,iws(iw))
      case ("cont")
         do ix1=ixIsmin1,ixIsmax1
            ws(ix1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iws(iw)) &
               = ws(ixIsmax1+1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iws(iw))
         end do
      case("periodic")
      case ("zero")
         ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iws(iw)) &
            = 0.0d0
      case ("specialT")
      call specialbound_usr(time,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
         ixGmax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,-iB,s)
      case ("special")
         ! skip it here, do AFTER all normal type boundaries are set
      case default
         write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
              "for variable iw=",iw," and side iB=",iB
      end select


      ! Special treatment when we have coarse neighbors in tangential dir:
      do is=1,2
         if (is_coarse(idir,is)) then

            ixIsmin1=ixImin1-kr(1,idir);ixIsmin2=ixImin2-kr(2,idir)
            ixIsmin3=ixImin3-kr(3,idir);
            ixIsmax1=ixImax1;ixIsmax2=ixImax2;ixIsmax3=ixImax3;

            if (is.eq.1) then ! coarse neighbor at beginning of block

               select case(idir)
                  
                  case (1)
                     ixIsmax1=min(ixGlo1+dixB-1,ixIsmax1);
                     surf=>s%geo%surfaceC1
                  

                  case (2)
                     ixIsmax2=min(ixGlo2+dixB-1,ixIsmax2);
                     surf=>s%geo%surfaceC2
                  

                  case (3)
                     ixIsmax3=min(ixGlo3+dixB-1,ixIsmax3);
                     surf=>s%geo%surfaceC3
                  ! I am hating LASY just right now.
               end select

            else ! coarse neighbor at end of block

               select case(idir)
                  
                  case (1)
                     ixIsmin1=max(ixGhi1-dixB,ixIsmin1);
                     surf=>s%geo%surfaceC1
                  

                  case (2)
                     ixIsmin2=max(ixGhi2-dixB,ixIsmin2);
                     surf=>s%geo%surfaceC2
                  

                  case (3)
                     ixIsmin3=max(ixGhi3-dixB,ixIsmin3);
                     surf=>s%geo%surfaceC3
                  ! I am hating LASY just right now.
               end select

            end if


            select case(typeB(iw,iB))
            case ("cont")
               do ix1=ixIsmin1,ixIsmax1
                  where ( (surf(ixIsmax1+1,ixIsmin2:ixIsmax2,&
                     ixIsmin3:ixIsmax3) + surf(ixIsmax1+2,ixIsmin2:ixIsmax2,&
                     ixIsmin3:ixIsmax3)) .gt. smalldouble )
                     ws(ix1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iws(iw)) &
                        = (ws(ixIsmax1+1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,&
                        iws(iw))*surf(ixIsmax1+1,ixIsmin2:ixIsmax2,&
                        ixIsmin3:ixIsmax3) &
                          + ws(ixIsmax1+2,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,&
                             iws(iw))*surf(ixIsmax1+2,ixIsmin2:ixIsmax2,&
                             ixIsmin3:ixIsmax3) ) / (surf(ixIsmax1+1,&
                             ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3)+surf&
                             (ixIsmax1+2,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3))
                  elsewhere
                     ws(ix1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iws(iw)) &
                        = zero
                  end where
               end do

            case default
 !Only implemented cont type, symm and asymm will also work if problem respects symmetry
               ! 3D case will need to average four cells together!!!
            end select

         end if
      end do

   end if

end if

end do


      ! Now that the tangential components are set,
      ! fill the normal components using a prescription for the divergence.
      ! This prescription is given by the typeB for the normal component.
do iw=1,nws
        ! When using more than one staggered vector field
        idir=mod(iw-1,ndir)+1
        ! Consider only normal direction
        if (idir.eq.1) then
          ixIsmin1=ixImin1-kr(1,1);ixIsmin2=ixImin2-kr(2,1)
          ixIsmin3=ixImin3-kr(3,1);ixIsmax1=ixImax1-kr(1,1)
          ixIsmax2=ixImax2-kr(2,1);ixIsmax3=ixImax3-kr(3,1);
          jxImin1=ixImin1+dixB*kr(1,1);jxImin2=ixImin2+dixB*kr(2,1)
          jxImin3=ixImin3+dixB*kr(3,1);jxImax1=ixImax1+dixB*kr(1,1)
          jxImax2=ixImax2+dixB*kr(2,1);jxImax3=ixImax3+dixB*kr(3,1);

          ! Calculate divergence and partial divergence
          Q=0
          if (.not.zerodiv) call div_staggered(jxImin1,jxImin2,jxImin3,&
             jxImax1,jxImax2,jxImax3,s,Q(jxImin1:jxImax1,jxImin2:jxImax2,&
             jxImin3:jxImax3))
          select case(typeB(iw+b0_,iB))
          case("symm")
            ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iw)=zero
            do ix1=0,dixB-1
              call div_staggered(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                 ixImax3,s,Qp(ixImin1:ixImax1,ixImin2:ixImax2,&
                 ixImin3:ixImax3))
              ws(ixIsmax1-ix1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iw)=&
             -(Q(jxImin1+ix1,jxImin2:jxImax2,jxImin3:jxImax3)*s%geo%dvolume&
                (jxImin1+ix1,jxImin2:jxImax2,jxImin3:jxImax3)&
             -Qp(ixImax1-ix1,ixImin2:ixImax2,ixImin3:ixImax3)*s%geo%dvolume&
                (ixImax1-ix1,ixImin2:ixImax2,ixImin3:ixImax3))&
             /s%geo%surfaceC1(ixIsmax1-ix1,ixIsmin2:ixIsmax2,&
                ixIsmin3:ixIsmax3)
            end do
          case("asymm")
             ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iw)=zero
            do ix1=0,dixB-1
              call div_staggered(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                 ixImax3,s,Qp(ixImin1:ixImax1,ixImin2:ixImax2,&
                 ixImin3:ixImax3))
              ws(ixIsmax1-ix1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iw)=&
            -(-Q(jxImin1+ix1,jxImin2:jxImax2,jxImin3:jxImax3)*s%geo%dvolume&
               (jxImin1+ix1,jxImin2:jxImax2,jxImin3:jxImax3)&
             -Qp(ixImax1-ix1,ixImin2:ixImax2,ixImin3:ixImax3)*s%geo%dvolume&
                (ixImax1-ix1,ixImin2:ixImax2,ixImin3:ixImax3))&
             /s%geo%surfaceC1(ixIsmax1-ix1,ixIsmin2:ixIsmax2,&
                ixIsmin3:ixIsmax3)
            end do
          case("cont", "zero")
            ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iw)=zero
            do ix1=0,dixB-1
              call div_staggered(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                 ixImax3,s,Qp(ixImin1:ixImax1,ixImin2:ixImax2,&
                 ixImin3:ixImax3))
              ws(ixIsmax1-ix1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iw)=&
             -(Q(jxImin1,jxImin2:jxImax2,jxImin3:jxImax3)*s%geo%dvolume&
                (jxImin1,jxImin2:jxImax2,jxImin3:jxImax3)&
             -Qp(ixImax1-ix1,ixImin2:ixImax2,ixImin3:ixImax3)*s%geo%dvolume&
                (ixImax1-ix1,ixImin2:ixImax2,ixImin3:ixImax3))&
              /s%geo%surfaceC1(ixIsmax1-ix1,ixIsmin2:ixIsmax2,&
                 ixIsmin3:ixIsmax3)
            end do
          case("periodic")
          end select
        end if
      end do

      ! Fill cell averages
      call faces2centers(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,s)

      
      ! Apply on primitive variables?
      if (primitiveB(iside,idims)) then
         w(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) &
            = wsave(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) !restore from backup of the state in domain since p2c/c2p might lead to slight changes in w
         call conserve(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
            ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,w,x,patchfalse)
      end if
      
   end if 
case (2)
if (iside==2) then

   
   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! maximal boundary
   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      iB=ismax2
      ixImin1=ixBmin1;ixImin2=ixBmax2+1-dixB;ixImin3=ixBmin3;
      ixImax1=ixBmax1;ixImax2=ixBmax2;ixImax3=ixBmax3;

      !!!!!!! Make tangential ranges greater for staggered components
      !!!!!!! Add fill the normal component

      ! Apply on primitive variables?
      if (primitiveB(iside,idims)) then
         ixDmin1=ixBmin1;ixDmin2=ixBmin2+dixB;ixDmin3=ixBmin3;
         ixDmax1=ixBmax1;ixDmax2=ixBmax2-dixB;ixDmax3=ixBmax3;
         wsave(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) &
            = w(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) !make a backup of the state in domain since p2c/c2p might lead to slight changes in w
         call primitive(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
            ixDmin1,ixDmin2,ixDmin3,ixDmax1,ixDmax2,ixDmax3,w,x)
      end if
      
      ! cont/symm/asymm types
      do iw=1,nwflux+nwaux
 ! If the variable has no staggered counterpart
         if ((iws(iw).lt.1).or.(iws(iw).gt.nws)) then

         select case (typeB(iw,iB))
         case ("symm")
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               = w(ixImin1:ixImax1,ixImin2-1:ixImin2-dixB:-1,ixImin3:ixImax3,&
               iw)
         case ("asymm")
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               =-w(ixImin1:ixImax1,ixImin2-1:ixImin2-dixB:-1,ixImin3:ixImax3,&
               iw)
         case ("cont")
            do ix2=ixImin2,ixImax2
               w(ixImin1:ixImax1,ix2,ixImin3:ixImax3,iw) = w(ixImin1:ixImax1,&
                  ixImin2-1,ixImin3:ixImax3,iw)
            end do
         case("noinflow")
            if (iw==1+2)then
              do ix2=ixImin2,ixImax2
                  w(ixImin1:ixImax1,ix2,ixImin3:ixImax3,iw) &
                     = max(w(ixImin1:ixImax1,ixImin2-1,ixImin3:ixImax3,iw),&
                     zero)
              end do
            else
              do ix2=ixImin2,ixImax2
                  w(ixImin1:ixImax1,ix2,ixImin3:ixImax3,iw) &
                     = w(ixImin1:ixImax1,ixImin2-1,ixImin3:ixImax3,iw)
              end do
            end if
         case("limitinflow")
            if (iw==1+2)then
              do ix2=ixImin2,ixImax2
                  w(ixImin1:ixImax1,ix2,ixImin3:ixImax3,iw) &
                     = max(w(ixImin1:ixImax1,ixImin2-1,ixImin3:ixImax3,iw), &
                                           w(ixImin1:ixImax1,ixImin2-1,&
                                              ixImin3:ixImax3,iw)*ratebdflux)
              end do
            else
              do ix2=ixImin2,ixImax2
                  w(ixImin1:ixImax1,ix2,ixImin3:ixImax3,iw) &
                     = w(ixImin1:ixImax1,ixImin2-1,ixImin3:ixImax3,iw)
              end do
           end if
        case ("polefix")
           ! First overwrite cells in domain
           ixIpmin1=ixImin1;ixIpmin2=ixImin2-dixPolefix;ixIpmin3=ixImin3;
           ixIpmax1=ixImax1;ixIpmax2=ixImin2-1;ixIpmax3=ixImax3;
           do ix2=ixIpmin2,ixIpmax2
              w(ixIpmin1:ixIpmax1,ix2,ixIpmin3:ixIpmax3,iw) &
                 = w(ixIpmin1:ixIpmax1,ixIpmin2-1,ixIpmin3:ixIpmax3,iw)
           end do
           ! Now apply symm boundary
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               = w(ixImin1:ixImax1,ixImin2-1:ixImin2-dixB:-1,ixImin3:ixImax3,&
               iw)
        case ("apolefix")
           ! First interpolate cells in domain down to zero at pole
           ixIpmin1=ixImin1;ixIpmin2=ixImin2-dixPolefix;ixIpmin3=ixImin3;
           ixIpmax1=ixImax1;ixIpmax2=ixImax2-1;ixIpmax3=ixImax3;
           do ix2=ixIpmin2,ixIpmax2
              w(ixIpmin1:ixIpmax1,ix2,ixIpmin3:ixIpmax3,iw) &
                 = w(ixIpmin1:ixIpmax1,ixIpmin2-1,ixIpmin3:ixIpmax3,iw) &
                   * (xprobmax2-x(ixIpmin1:ixIpmax1,ix2,ixIpmin3:ixIpmax3,&
                      2)) &
                   / (xprobmax2-x(ixIpmin1:ixIpmax1,ixIpmin2-1,&
                      ixIpmin3:ixIpmax3,2))
           end do
        case ("hotaka")
           ! First overwrite cells in domain with zero
           ixIpmin1=ixImin1;ixIpmin2=ixImin2-dixHotaka;ixIpmin3=ixImin3;
           ixIpmax1=ixImax1;ixIpmax2=ixImin2-1;ixIpmax3=ixImax3;
           do ix2=ixIpmin2,ixIpmax2
              w(ixIpmin1:ixIpmax1,ix2,ixIpmin3:ixIpmax3,iw) = zero
           end do
           ! Now apply symm boundary
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               = w(ixImin1:ixImax1,ixImin2-1:ixImin2-dixB:-1,ixImin3:ixImax3,&
               iw)
        case ("ahotaka")
            ! First overwrite cells in domain with zero
           ixIpmin1=ixImin1;ixIpmin2=ixImin2-dixHotaka;ixIpmin3=ixImin3;
           ixIpmax1=ixImax1;ixIpmax2=ixImax2-1;ixIpmax3=ixImax3;
           do ix2=ixIpmin2,ixIpmax2
              w(ixIpmin1:ixIpmax1,ix2,ixIpmin3:ixIpmax3,iw) = zero
           end do
           ! Now apply asymm boundary
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               =-w(ixImin1:ixImax1,ixImin2-1:ixImin2-dixB:-1,ixImin3:ixImax3,&
               iw)
         case ("special")
            ! skip it here, do AFTER all normal type boundaries are set
         case ("aperiodic")
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               = - w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw)
         case ("periodic")
!            call mpistop("periodic bc info should come from neighbors")
         case default
            write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
               "for variable iw=",iw," and side iB=",iB
         end select


      else !If there is a staggered counterpart
 !At this stage, extrapolation is applied only to the tangential components
        idir=mod(iws(iw)-1,ndir)+1
        if (idir.ne.2) then
        ixIsmax1=ixImax1;ixIsmax2=ixImax2;ixIsmax3=ixImax3;
        ixIsmin1=ixImin1-kr(1,idir);ixIsmin2=ixImin2-kr(2,idir)
        ixIsmin3=ixImin3-kr(3,idir);
        select case(typeB(iw,iB))
        case ("symm")
          ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iws(iw)) &
             = ws(ixIsmin1:ixIsmax1,ixIsmin2-1:ixIsmin2-dixB:-1,&
             ixIsmin3:ixIsmax3,iws(iw))
        case ("asymm")
          ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iws(iw)) &
             =-ws(ixIsmin1:ixIsmax1,ixIsmin2-1:ixIsmin2-dixB:-1,&
             ixIsmin3:ixIsmax3,iws(iw))
        case ("cont")
          do ix2=ixIsmin2,ixIsmax2
             ws(ixIsmin1:ixIsmax1,ix2,ixIsmin3:ixIsmax3,iws(iw)) &
                = ws(ixIsmin1:ixIsmax1,ixIsmin2-1,ixIsmin3:ixIsmax3,iws(iw))
          end do
        case ("periodic")
        case ("zero")
           ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iws(iw)) &
              = 0.0d0
        case ("specialT")
           call specialbound_usr(time,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
              ixGmax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,-iB,s)
        case ("special")
           ! skip it here, do AFTER all normal type boundaries are set
        case default
          write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
             "for variable iw=",iw," and side iB=",iB
        end select

                  
      ! Special treatment when we have coarse neighbors in tangential dir:
      do is=1,2
         if (is_coarse(idir,is)) then

            ixIsmin1=ixImin1-kr(1,idir);ixIsmin2=ixImin2-kr(2,idir)
            ixIsmin3=ixImin3-kr(3,idir);
            ixIsmax1=ixImax1;ixIsmax2=ixImax2;ixIsmax3=ixImax3;

            if (is.eq.1) then ! coarse neighbor at beginning of block

               select case(idir)
                  
                  case (1)
                     ixIsmax1=min(ixGlo1+dixB-1,ixIsmax1);
                     surf=>s%geo%surfaceC1
                  

                  case (2)
                     ixIsmax2=min(ixGlo2+dixB-1,ixIsmax2);
                     surf=>s%geo%surfaceC2
                  

                  case (3)
                     ixIsmax3=min(ixGlo3+dixB-1,ixIsmax3);
                     surf=>s%geo%surfaceC3
                  ! I am hating LASY just right now.
               end select

            else ! coarse neighbor at end of block

               select case(idir)
                  
                  case (1)
                     ixIsmin1=max(ixGhi1-dixB,ixIsmin1);
                     surf=>s%geo%surfaceC1
                  

                  case (2)
                     ixIsmin2=max(ixGhi2-dixB,ixIsmin2);
                     surf=>s%geo%surfaceC2
                  

                  case (3)
                     ixIsmin3=max(ixGhi3-dixB,ixIsmin3);
                     surf=>s%geo%surfaceC3
                  ! I am hating LASY just right now.
               end select

            end if


            select case(typeB(iw,iB))
            case ("cont")
               do ix2=ixIsmin2,ixIsmax2
                  where ((surf(ixIsmin1:ixIsmax1,ixIsmin2-1,&
                     ixIsmin3:ixIsmax3) + surf(ixIsmin1:ixIsmax1,ixIsmin2-2,&
                     ixIsmin3:ixIsmax3)) .gt. smalldouble )
                     ws(ixIsmin1:ixIsmax1,ix2,ixIsmin3:ixIsmax3,iws(iw)) &
                        = (ws(ixIsmin1:ixIsmax1,ixIsmin2-1,ixIsmin3:ixIsmax3,&
                        iws(iw))*surf(ixIsmin1:ixIsmax1,ixIsmin2-1,&
                        ixIsmin3:ixIsmax3) &
                          + ws(ixIsmin1:ixIsmax1,ixIsmin2-2,ixIsmin3:ixIsmax3,&
                             iws(iw))*surf(ixIsmin1:ixIsmax1,ixIsmin2-2,&
                             ixIsmin3:ixIsmax3) ) / (surf(ixIsmin1:ixIsmax1,&
                             ixIsmin2-1,ixIsmin3:ixIsmax3)+surf&
                             (ixIsmin1:ixIsmax1,ixIsmin2-2,ixIsmin3:ixIsmax3))
                  elsewhere
                     ws(ixIsmin1:ixIsmax1,ix2,ixIsmin3:ixIsmax3,iws(iw)) &
                        = zero
                  end where
               end do

            case default
 !Only implemented cont type, symm and asymm will also work if problem respects symmetry
               ! 3D case will need to average four cells together!!!
            end select

         end if
      end do

   end if
end if

      end do

      ! Now that the tangential components are set,
      ! fill the normal components using a prescription for the divergence.
      ! This prescription is given by the typeB for the normal component.
      do iw=1,nws
        ! When using more than one staggered vector field
        idir=mod(iw-1,ndir)+1
        ! Consider only normal direction
        if (idir.eq.2) then
          ixIsmin1=ixImin1;ixIsmin2=ixImin2;ixIsmin3=ixImin3;ixIsmax1=ixImax1
          ixIsmax2=ixImax2;ixIsmax3=ixImax3;
          hxImin1=ixImin1-dixB*kr(1,2);hxImin2=ixImin2-dixB*kr(2,2)
          hxImin3=ixImin3-dixB*kr(3,2);hxImax1=ixImax1-dixB*kr(1,2)
          hxImax2=ixImax2-dixB*kr(2,2);hxImax3=ixImax3-dixB*kr(3,2);

          ! Calculate divergence and partial divergence
          Q=zero
          if (.not.zerodiv) call div_staggered(hxImin1,hxImin2,hxImin3,&
             hxImax1,hxImax2,hxImax3,s,Q(hxImin1:hxImax1,hxImin2:hxImax2,&
             hxImin3:hxImax3))
          select case(typeB(iw+b0_,iB))
          case("symm")
            ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iw)=zero
            do ix2=0,dixB-1
              call div_staggered(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                 ixImax3,s,Qp(ixImin1:ixImax1,ixImin2:ixImax2,&
                 ixImin3:ixImax3))
              ws(ixIsmin1:ixIsmax1,ixIsmin2+ix2,ixIsmin3:ixIsmax3,iw)=&
              (Q(hxImin1:hxImax1,hxImax2-ix2,hxImin3:hxImax3)*s%geo%dvolume&
                 (hxImin1:hxImax1,hxImax2-ix2,hxImin3:hxImax3)&
             -Qp(ixImin1:ixImax1,ixImin2+ix2,ixImin3:ixImax3)*s%geo%dvolume&
                (ixImin1:ixImax1,ixImin2+ix2,ixImin3:ixImax3))&
              /s%geo%surfaceC2(ixIsmin1:ixIsmax1,ixIsmin2+ix2,&
                 ixIsmin3:ixIsmax3)
            end do
          case("asymm")
            ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iw)=zero
            do ix2=0,dixB-1
              call div_staggered(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                 ixImax3,s,Qp(ixImin1:ixImax1,ixImin2:ixImax2,&
                 ixImin3:ixImax3))
              ws(ixIsmin1:ixIsmax1,ixIsmin2+ix2,ixIsmin3:ixIsmax3,iw)=&
             (-Q(hxImin1:hxImax1,hxImax2-ix2,hxImin3:hxImax3)*s%geo%dvolume&
                (hxImin1:hxImax1,hxImax2-ix2,hxImin3:hxImax3)&
             -Qp(ixImin1:ixImax1,ixImin2+ix2,ixImin3:ixImax3)*s%geo%dvolume&
                (ixImin1:ixImax1,ixImin2+ix2,ixImin3:ixImax3))&
              /s%geo%surfaceC2(ixIsmin1:ixIsmax1,ixIsmin2+ix2,&
                 ixIsmin3:ixIsmax3)
            end do
          case("cont", "zero")
            ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iw)=zero
            do ix2=0,dixB-1
              call div_staggered(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                 ixImax3,s,Qp(ixImin1:ixImax1,ixImin2:ixImax2,&
                 ixImin3:ixImax3))
              ws(ixIsmin1:ixIsmax1,ixIsmin2+ix2,ixIsmin3:ixIsmax3,iw)=&
              (Q(hxImin1:hxImax1,hxImax2,hxImin3:hxImax3)*s%geo%dvolume&
                 (hxImin1:hxImax1,hxImax2,hxImin3:hxImax3)&
             -Qp(ixImin1:ixImax1,ixImin2+ix2,ixImin3:ixImax3)*s%geo%dvolume&
                (ixImin1:ixImax1,ixImin2+ix2,ixImin3:ixImax3))&
              /s%geo%surfaceC2(ixIsmin1:ixIsmax1,ixIsmin2+ix2,&
                 ixIsmin3:ixIsmax3)
            end do
          case("periodic")
          end select
        end if
      end do

      ! Fill cell averages
      call faces2centers(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,s)



      ! Apply on primitive variables?
      if (primitiveB(iside,idims)) then
         w(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) &
            = wsave(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) !restore from backup of the state in domain since p2c/c2p might lead to slight changes in w
         call conserve(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
            ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,w,x,patchfalse)
      end if
      
   else
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! minimal boundary
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      iB=ismin2
      ixImin1=ixBmin1;ixImin2=ixBmin2;ixImin3=ixBmin3;
      ixImax1=ixBmax1;ixImax2=ixBmin2-1+dixB;ixImax3=ixBmax3;

      ! Apply on primitive variables?
      if (primitiveB(iside,idims)) then         
         ixDmin1=ixBmin1;ixDmin2=ixBmin2+dixB;ixDmin3=ixBmin3;
         ixDmax1=ixBmax1;ixDmax2=ixBmax2-dixB;ixDmax3=ixBmax3;
         wsave(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) &
            = w(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) !make a backup of the state in domain since p2c/c2p might lead to slight changes in w
         call primitive(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
            ixDmin1,ixDmin2,ixDmin3,ixDmax1,ixDmax2,ixDmax3,w,x)
      end if
 
      ! cont/symm/asymm types
      do iw=1,nwflux+nwaux
 ! If the variable has no staggered counterpart
if ((iws(iw).lt.1).or.(iws(iw).gt.nws)) then
   
   select case (typeB(iw,iB))
   case ("symm")
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         = w(ixImin1:ixImax1,ixImax2+dixB:ixImax2+1:-1,ixImin3:ixImax3,iw)
   case ("asymm")
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         =-w(ixImin1:ixImax1,ixImax2+dixB:ixImax2+1:-1,ixImin3:ixImax3,iw)
   case ("cont")
      do ix2=ixImin2,ixImax2
         w(ixImin1:ixImax1,ix2,ixImin3:ixImax3,iw) = w(ixImin1:ixImax1,&
            ixImax2+1,ixImin3:ixImax3,iw)
      end do
   case("noinflow")
      if (iw==1+2)then
         do ix2=ixImin2,ixImax2
            w(ixImin1:ixImax1,ix2,ixImin3:ixImax3,iw) = min(w(ixImin1:ixImax1,&
               ixImax2+1,ixImin3:ixImax3,iw),zero)
         end do
      else
         do ix2=ixImin2,ixImax2
            w(ixImin1:ixImax1,ix2,ixImin3:ixImax3,iw) = w(ixImin1:ixImax1,&
               ixImax2+1,ixImin3:ixImax3,iw)
         end do
      end if
   case("limitinflow")
      if (iw==1+2)then
         do ix2=ixImin2,ixImax2
            w(ixImin1:ixImax1,ix2,ixImin3:ixImax3,iw) = min(w(ixImin1:ixImax1,&
               ixImax2+1,ixImin3:ixImax3,iw), &
                 w(ixImin1:ixImax1,ixImax2+1,ixImin3:ixImax3,iw)*ratebdflux)
         end do
      else
         do ix2=ixImin2,ixImax2
            w(ixImin1:ixImax1,ix2,ixImin3:ixImax3,iw) = w(ixImin1:ixImax1,&
               ixImax2+1,ixImin3:ixImax3,iw)
         end do
      end if
   case ("polefix")
      ! First overwrite cells in domain
      ixIpmin1=ixImin1;ixIpmin2=ixImax2+1;ixIpmin3=ixImin3;
      ixIpmax1=ixImax1;ixIpmax2=ixImax2+dixPolefix;ixIpmax3=ixImax3;
      do ix2=ixIpmin2,ixIpmax2
         w(ixIpmin1:ixIpmax1,ix2,ixIpmin3:ixIpmax3,iw) = w(ixIpmin1:ixIpmax1,&
            ixIpmax2+1,ixIpmin3:ixIpmax3,iw)
      end do
      ! Now apply symm boundary
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         = w(ixImin1:ixImax1,ixImax2+dixB:ixImax2+1:-1,ixImin3:ixImax3,iw)
   case ("apolefix")
      ! First interpolate cells in domain down to zero at pole
      ixIpmin1=ixImin1;ixIpmin2=ixImax2+1;ixIpmin3=ixImin3;
      ixIpmax1=ixImax1;ixIpmax2=ixImax2+dixPolefix;ixIpmax3=ixImax3;
      do ix2=ixIpmin2,ixIpmax2
         w(ixIpmin1:ixIpmax1,ix2,ixIpmin3:ixIpmax3,iw) = w(ixIpmin1:ixIpmax1,&
            ixIpmax2+1,ixIpmin3:ixIpmax3,iw) &
              * (x(ixIpmin1:ixIpmax1,ix2,ixIpmin3:ixIpmax3,2))&
                 /x(ixIpmin1:ixIpmax1,ixIpmax2+1,ixIpmin3:ixIpmax3,2)
      end do
      ! Now apply asymm boundary
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         =-w(ixImin1:ixImax1,ixImax2+dixB:ixImax2+1:-1,ixImin3:ixImax3,iw)
   case ("hotaka")
      ! First overwrite cells in domain with zero
      ixIpmin1=ixImin1;ixIpmin2=ixImax2+1;ixIpmin3=ixImin3;
      ixIpmax1=ixImax1;ixIpmax2=ixImax2+dixHotaka;ixIpmax3=ixImax3;
      do ix2=ixIpmin2,ixIpmax2
         w(ixIpmin1:ixIpmax1,ix2,ixIpmin3:ixIpmax3,iw) = zero
      end do
      ! Now apply symm boundary
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         = w(ixImin1:ixImax1,ixImax2+dixB:ixImax2+1:-1,ixImin3:ixImax3,iw)
   case ("ahotaka")
      ! First overwrite cells in domain with zero
      ixIpmin1=ixImin1;ixIpmin2=ixImax2+1;ixIpmin3=ixImin3;
      ixIpmax1=ixImax1;ixIpmax2=ixImax2+dixHotaka;ixIpmax3=ixImax3;
      do ix2=ixIpmin2,ixIpmax2
         w(ixIpmin1:ixIpmax1,ix2,ixIpmin3:ixIpmax3,iw) = zero
      end do
      ! Now apply asymm boundary
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         =-w(ixImin1:ixImax1,ixImax2+dixB:ixImax2+1:-1,ixImin3:ixImax3,iw)
   case ("special")
      ! skip it here, do AFTER all normal type boundaries are set
   case ("aperiodic")
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         = - w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw)
   case ("periodic")
      !            call mpistop("periodic bc info should come from neighbors")
   case default
      write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
           "for variable iw=",iw," and side iB=",iB
   end select


else !If there is a staggered counterpart
   ! At this stage, extrapolation is applied only to the tangential components
   idir=mod(iws(iw)-1,ndir)+1 ! vector direction
   if (idir.ne.2) then
      ixIsmax1=ixImax1;ixIsmax2=ixImax2;ixIsmax3=ixImax3;
      ixIsmin1=ixImin1-kr(1,idir);ixIsmin2=ixImin2-kr(2,idir)
      ixIsmin3=ixImin3-kr(3,idir);
      select case(typeB(iw,iB))
      case ("symm")
         ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iws(iw)) &
            = ws(ixIsmin1:ixIsmax1,ixIsmax2+dixB:ixIsmax2+1:-1,&
            ixIsmin3:ixIsmax3,iws(iw))
      case ("asymm")
         ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iws(iw)) &
            =-ws(ixIsmin1:ixIsmax1,ixIsmax2+dixB:ixIsmax2+1:-1,&
            ixIsmin3:ixIsmax3,iws(iw))
      case ("cont")
         do ix2=ixIsmin2,ixIsmax2
            ws(ixIsmin1:ixIsmax1,ix2,ixIsmin3:ixIsmax3,iws(iw)) &
               = ws(ixIsmin1:ixIsmax1,ixIsmax2+1,ixIsmin3:ixIsmax3,iws(iw))
         end do
      case("periodic")
      case ("zero")
         ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iws(iw)) &
            = 0.0d0
      case ("specialT")
      call specialbound_usr(time,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
         ixGmax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,-iB,s)
      case ("special")
         ! skip it here, do AFTER all normal type boundaries are set
      case default
         write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
              "for variable iw=",iw," and side iB=",iB
      end select


      ! Special treatment when we have coarse neighbors in tangential dir:
      do is=1,2
         if (is_coarse(idir,is)) then

            ixIsmin1=ixImin1-kr(1,idir);ixIsmin2=ixImin2-kr(2,idir)
            ixIsmin3=ixImin3-kr(3,idir);
            ixIsmax1=ixImax1;ixIsmax2=ixImax2;ixIsmax3=ixImax3;

            if (is.eq.1) then ! coarse neighbor at beginning of block

               select case(idir)
                  
                  case (1)
                     ixIsmax1=min(ixGlo1+dixB-1,ixIsmax1);
                     surf=>s%geo%surfaceC1
                  

                  case (2)
                     ixIsmax2=min(ixGlo2+dixB-1,ixIsmax2);
                     surf=>s%geo%surfaceC2
                  

                  case (3)
                     ixIsmax3=min(ixGlo3+dixB-1,ixIsmax3);
                     surf=>s%geo%surfaceC3
                  ! I am hating LASY just right now.
               end select

            else ! coarse neighbor at end of block

               select case(idir)
                  
                  case (1)
                     ixIsmin1=max(ixGhi1-dixB,ixIsmin1);
                     surf=>s%geo%surfaceC1
                  

                  case (2)
                     ixIsmin2=max(ixGhi2-dixB,ixIsmin2);
                     surf=>s%geo%surfaceC2
                  

                  case (3)
                     ixIsmin3=max(ixGhi3-dixB,ixIsmin3);
                     surf=>s%geo%surfaceC3
                  ! I am hating LASY just right now.
               end select

            end if


            select case(typeB(iw,iB))
            case ("cont")
               do ix2=ixIsmin2,ixIsmax2
                  where ( (surf(ixIsmin1:ixIsmax1,ixIsmax2+1,&
                     ixIsmin3:ixIsmax3) + surf(ixIsmin1:ixIsmax1,ixIsmax2+2,&
                     ixIsmin3:ixIsmax3)) .gt. smalldouble )
                     ws(ixIsmin1:ixIsmax1,ix2,ixIsmin3:ixIsmax3,iws(iw)) &
                        = (ws(ixIsmin1:ixIsmax1,ixIsmax2+1,ixIsmin3:ixIsmax3,&
                        iws(iw))*surf(ixIsmin1:ixIsmax1,ixIsmax2+1,&
                        ixIsmin3:ixIsmax3) &
                          + ws(ixIsmin1:ixIsmax1,ixIsmax2+2,ixIsmin3:ixIsmax3,&
                             iws(iw))*surf(ixIsmin1:ixIsmax1,ixIsmax2+2,&
                             ixIsmin3:ixIsmax3) ) / (surf(ixIsmin1:ixIsmax1,&
                             ixIsmax2+1,ixIsmin3:ixIsmax3)+surf&
                             (ixIsmin1:ixIsmax1,ixIsmax2+2,ixIsmin3:ixIsmax3))
                  elsewhere
                     ws(ixIsmin1:ixIsmax1,ix2,ixIsmin3:ixIsmax3,iws(iw)) &
                        = zero
                  end where
               end do

            case default
 !Only implemented cont type, symm and asymm will also work if problem respects symmetry
               ! 3D case will need to average four cells together!!!
            end select

         end if
      end do

   end if

end if

end do


      ! Now that the tangential components are set,
      ! fill the normal components using a prescription for the divergence.
      ! This prescription is given by the typeB for the normal component.
do iw=1,nws
        ! When using more than one staggered vector field
        idir=mod(iw-1,ndir)+1
        ! Consider only normal direction
        if (idir.eq.2) then
          ixIsmin1=ixImin1-kr(1,2);ixIsmin2=ixImin2-kr(2,2)
          ixIsmin3=ixImin3-kr(3,2);ixIsmax1=ixImax1-kr(1,2)
          ixIsmax2=ixImax2-kr(2,2);ixIsmax3=ixImax3-kr(3,2);
          jxImin1=ixImin1+dixB*kr(1,2);jxImin2=ixImin2+dixB*kr(2,2)
          jxImin3=ixImin3+dixB*kr(3,2);jxImax1=ixImax1+dixB*kr(1,2)
          jxImax2=ixImax2+dixB*kr(2,2);jxImax3=ixImax3+dixB*kr(3,2);

          ! Calculate divergence and partial divergence
          Q=0
          if (.not.zerodiv) call div_staggered(jxImin1,jxImin2,jxImin3,&
             jxImax1,jxImax2,jxImax3,s,Q(jxImin1:jxImax1,jxImin2:jxImax2,&
             jxImin3:jxImax3))
          select case(typeB(iw+b0_,iB))
          case("symm")
            ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iw)=zero
            do ix2=0,dixB-1
              call div_staggered(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                 ixImax3,s,Qp(ixImin1:ixImax1,ixImin2:ixImax2,&
                 ixImin3:ixImax3))
              ws(ixIsmin1:ixIsmax1,ixIsmax2-ix2,ixIsmin3:ixIsmax3,iw)=&
             -(Q(jxImin1:jxImax1,jxImin2+ix2,jxImin3:jxImax3)*s%geo%dvolume&
                (jxImin1:jxImax1,jxImin2+ix2,jxImin3:jxImax3)&
             -Qp(ixImin1:ixImax1,ixImax2-ix2,ixImin3:ixImax3)*s%geo%dvolume&
                (ixImin1:ixImax1,ixImax2-ix2,ixImin3:ixImax3))&
             /s%geo%surfaceC2(ixIsmin1:ixIsmax1,ixIsmax2-ix2,&
                ixIsmin3:ixIsmax3)
            end do
          case("asymm")
             ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iw)=zero
            do ix2=0,dixB-1
              call div_staggered(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                 ixImax3,s,Qp(ixImin1:ixImax1,ixImin2:ixImax2,&
                 ixImin3:ixImax3))
              ws(ixIsmin1:ixIsmax1,ixIsmax2-ix2,ixIsmin3:ixIsmax3,iw)=&
            -(-Q(jxImin1:jxImax1,jxImin2+ix2,jxImin3:jxImax3)*s%geo%dvolume&
               (jxImin1:jxImax1,jxImin2+ix2,jxImin3:jxImax3)&
             -Qp(ixImin1:ixImax1,ixImax2-ix2,ixImin3:ixImax3)*s%geo%dvolume&
                (ixImin1:ixImax1,ixImax2-ix2,ixImin3:ixImax3))&
             /s%geo%surfaceC2(ixIsmin1:ixIsmax1,ixIsmax2-ix2,&
                ixIsmin3:ixIsmax3)
            end do
          case("cont", "zero")
            ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iw)=zero
            do ix2=0,dixB-1
              call div_staggered(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                 ixImax3,s,Qp(ixImin1:ixImax1,ixImin2:ixImax2,&
                 ixImin3:ixImax3))
              ws(ixIsmin1:ixIsmax1,ixIsmax2-ix2,ixIsmin3:ixIsmax3,iw)=&
             -(Q(jxImin1:jxImax1,jxImin2,jxImin3:jxImax3)*s%geo%dvolume&
                (jxImin1:jxImax1,jxImin2,jxImin3:jxImax3)&
             -Qp(ixImin1:ixImax1,ixImax2-ix2,ixImin3:ixImax3)*s%geo%dvolume&
                (ixImin1:ixImax1,ixImax2-ix2,ixImin3:ixImax3))&
              /s%geo%surfaceC2(ixIsmin1:ixIsmax1,ixIsmax2-ix2,&
                 ixIsmin3:ixIsmax3)
            end do
          case("periodic")
          end select
        end if
      end do

      ! Fill cell averages
      call faces2centers(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,s)

      
      ! Apply on primitive variables?
      if (primitiveB(iside,idims)) then
         w(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) &
            = wsave(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) !restore from backup of the state in domain since p2c/c2p might lead to slight changes in w
         call conserve(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
            ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,w,x,patchfalse)
      end if
      
   end if 
case (3)
if (iside==2) then

   
   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! maximal boundary
   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      iB=ismax3
      ixImin1=ixBmin1;ixImin2=ixBmin2;ixImin3=ixBmax3+1-dixB;
      ixImax1=ixBmax1;ixImax2=ixBmax2;ixImax3=ixBmax3;

      !!!!!!! Make tangential ranges greater for staggered components
      !!!!!!! Add fill the normal component

      ! Apply on primitive variables?
      if (primitiveB(iside,idims)) then
         ixDmin1=ixBmin1;ixDmin2=ixBmin2;ixDmin3=ixBmin3+dixB;
         ixDmax1=ixBmax1;ixDmax2=ixBmax2;ixDmax3=ixBmax3-dixB;
         wsave(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) &
            = w(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) !make a backup of the state in domain since p2c/c2p might lead to slight changes in w
         call primitive(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
            ixDmin1,ixDmin2,ixDmin3,ixDmax1,ixDmax2,ixDmax3,w,x)
      end if
      
      ! cont/symm/asymm types
      do iw=1,nwflux+nwaux
 ! If the variable has no staggered counterpart
         if ((iws(iw).lt.1).or.(iws(iw).gt.nws)) then

         select case (typeB(iw,iB))
         case ("symm")
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               = w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3-1:ixImin3-dixB:-1,&
               iw)
         case ("asymm")
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               =-w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3-1:ixImin3-dixB:-1,&
               iw)
         case ("cont")
            do ix3=ixImin3,ixImax3
               w(ixImin1:ixImax1,ixImin2:ixImax2,ix3,iw) = w(ixImin1:ixImax1,&
                  ixImin2:ixImax2,ixImin3-1,iw)
            end do
         case("noinflow")
            if (iw==1+3)then
              do ix3=ixImin3,ixImax3
                  w(ixImin1:ixImax1,ixImin2:ixImax2,ix3,iw) &
                     = max(w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3-1,iw),&
                     zero)
              end do
            else
              do ix3=ixImin3,ixImax3
                  w(ixImin1:ixImax1,ixImin2:ixImax2,ix3,iw) &
                     = w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3-1,iw)
              end do
            end if
         case("limitinflow")
            if (iw==1+3)then
              do ix3=ixImin3,ixImax3
                  w(ixImin1:ixImax1,ixImin2:ixImax2,ix3,iw) &
                     = max(w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3-1,iw), &
                                           w(ixImin1:ixImax1,ixImin2:ixImax2,&
                                              ixImin3-1,iw)*ratebdflux)
              end do
            else
              do ix3=ixImin3,ixImax3
                  w(ixImin1:ixImax1,ixImin2:ixImax2,ix3,iw) &
                     = w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3-1,iw)
              end do
           end if
        case ("polefix")
           ! First overwrite cells in domain
           ixIpmin1=ixImin1;ixIpmin2=ixImin2;ixIpmin3=ixImin3-dixPolefix;
           ixIpmax1=ixImax1;ixIpmax2=ixImax2;ixIpmax3=ixImin3-1;
           do ix3=ixIpmin3,ixIpmax3
              w(ixIpmin1:ixIpmax1,ixIpmin2:ixIpmax2,ix3,iw) &
                 = w(ixIpmin1:ixIpmax1,ixIpmin2:ixIpmax2,ixIpmin3-1,iw)
           end do
           ! Now apply symm boundary
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               = w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3-1:ixImin3-dixB:-1,&
               iw)
        case ("apolefix")
           ! First interpolate cells in domain down to zero at pole
           ixIpmin1=ixImin1;ixIpmin2=ixImin2;ixIpmin3=ixImin3-dixPolefix;
           ixIpmax1=ixImax1;ixIpmax2=ixImax2;ixIpmax3=ixImax3-1;
           do ix3=ixIpmin3,ixIpmax3
              w(ixIpmin1:ixIpmax1,ixIpmin2:ixIpmax2,ix3,iw) &
                 = w(ixIpmin1:ixIpmax1,ixIpmin2:ixIpmax2,ixIpmin3-1,iw) &
                   * (xprobmax3-x(ixIpmin1:ixIpmax1,ixIpmin2:ixIpmax2,ix3,&
                      3)) &
                   / (xprobmax3-x(ixIpmin1:ixIpmax1,ixIpmin2:ixIpmax2,&
                      ixIpmin3-1,3))
           end do
        case ("hotaka")
           ! First overwrite cells in domain with zero
           ixIpmin1=ixImin1;ixIpmin2=ixImin2;ixIpmin3=ixImin3-dixHotaka;
           ixIpmax1=ixImax1;ixIpmax2=ixImax2;ixIpmax3=ixImin3-1;
           do ix3=ixIpmin3,ixIpmax3
              w(ixIpmin1:ixIpmax1,ixIpmin2:ixIpmax2,ix3,iw) = zero
           end do
           ! Now apply symm boundary
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               = w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3-1:ixImin3-dixB:-1,&
               iw)
        case ("ahotaka")
            ! First overwrite cells in domain with zero
           ixIpmin1=ixImin1;ixIpmin2=ixImin2;ixIpmin3=ixImin3-dixHotaka;
           ixIpmax1=ixImax1;ixIpmax2=ixImax2;ixIpmax3=ixImax3-1;
           do ix3=ixIpmin3,ixIpmax3
              w(ixIpmin1:ixIpmax1,ixIpmin2:ixIpmax2,ix3,iw) = zero
           end do
           ! Now apply asymm boundary
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               =-w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3-1:ixImin3-dixB:-1,&
               iw)
         case ("special")
            ! skip it here, do AFTER all normal type boundaries are set
         case ("aperiodic")
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               = - w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw)
         case ("periodic")
!            call mpistop("periodic bc info should come from neighbors")
         case default
            write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
               "for variable iw=",iw," and side iB=",iB
         end select


      else !If there is a staggered counterpart
 !At this stage, extrapolation is applied only to the tangential components
        idir=mod(iws(iw)-1,ndir)+1
        if (idir.ne.3) then
        ixIsmax1=ixImax1;ixIsmax2=ixImax2;ixIsmax3=ixImax3;
        ixIsmin1=ixImin1-kr(1,idir);ixIsmin2=ixImin2-kr(2,idir)
        ixIsmin3=ixImin3-kr(3,idir);
        select case(typeB(iw,iB))
        case ("symm")
          ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iws(iw)) &
             = ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,&
             ixIsmin3-1:ixIsmin3-dixB:-1,iws(iw))
        case ("asymm")
          ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iws(iw)) &
             =-ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,&
             ixIsmin3-1:ixIsmin3-dixB:-1,iws(iw))
        case ("cont")
          do ix3=ixIsmin3,ixIsmax3
             ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ix3,iws(iw)) &
                = ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3-1,iws(iw))
          end do
        case ("periodic")
        case ("zero")
           ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iws(iw)) &
              = 0.0d0
        case ("specialT")
           call specialbound_usr(time,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
              ixGmax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,-iB,s)
        case ("special")
           ! skip it here, do AFTER all normal type boundaries are set
        case default
          write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
             "for variable iw=",iw," and side iB=",iB
        end select

                  
      ! Special treatment when we have coarse neighbors in tangential dir:
      do is=1,2
         if (is_coarse(idir,is)) then

            ixIsmin1=ixImin1-kr(1,idir);ixIsmin2=ixImin2-kr(2,idir)
            ixIsmin3=ixImin3-kr(3,idir);
            ixIsmax1=ixImax1;ixIsmax2=ixImax2;ixIsmax3=ixImax3;

            if (is.eq.1) then ! coarse neighbor at beginning of block

               select case(idir)
                  
                  case (1)
                     ixIsmax1=min(ixGlo1+dixB-1,ixIsmax1);
                     surf=>s%geo%surfaceC1
                  

                  case (2)
                     ixIsmax2=min(ixGlo2+dixB-1,ixIsmax2);
                     surf=>s%geo%surfaceC2
                  

                  case (3)
                     ixIsmax3=min(ixGlo3+dixB-1,ixIsmax3);
                     surf=>s%geo%surfaceC3
                  ! I am hating LASY just right now.
               end select

            else ! coarse neighbor at end of block

               select case(idir)
                  
                  case (1)
                     ixIsmin1=max(ixGhi1-dixB,ixIsmin1);
                     surf=>s%geo%surfaceC1
                  

                  case (2)
                     ixIsmin2=max(ixGhi2-dixB,ixIsmin2);
                     surf=>s%geo%surfaceC2
                  

                  case (3)
                     ixIsmin3=max(ixGhi3-dixB,ixIsmin3);
                     surf=>s%geo%surfaceC3
                  ! I am hating LASY just right now.
               end select

            end if


            select case(typeB(iw,iB))
            case ("cont")
               do ix3=ixIsmin3,ixIsmax3
                  where ((surf(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,&
                     ixIsmin3-1) + surf(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,&
                     ixIsmin3-2)) .gt. smalldouble )
                     ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ix3,iws(iw)) &
                        = (ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3-1,&
                        iws(iw))*surf(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,&
                        ixIsmin3-1) &
                          + ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3-2,&
                             iws(iw))*surf(ixIsmin1:ixIsmax1,&
                             ixIsmin2:ixIsmax2,ixIsmin3-2) ) &
                             / (surf(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,&
                             ixIsmin3-1)+surf(ixIsmin1:ixIsmax1,&
                             ixIsmin2:ixIsmax2,ixIsmin3-2))
                  elsewhere
                     ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ix3,iws(iw)) &
                        = zero
                  end where
               end do

            case default
 !Only implemented cont type, symm and asymm will also work if problem respects symmetry
               ! 3D case will need to average four cells together!!!
            end select

         end if
      end do

   end if
end if

      end do

      ! Now that the tangential components are set,
      ! fill the normal components using a prescription for the divergence.
      ! This prescription is given by the typeB for the normal component.
      do iw=1,nws
        ! When using more than one staggered vector field
        idir=mod(iw-1,ndir)+1
        ! Consider only normal direction
        if (idir.eq.3) then
          ixIsmin1=ixImin1;ixIsmin2=ixImin2;ixIsmin3=ixImin3;ixIsmax1=ixImax1
          ixIsmax2=ixImax2;ixIsmax3=ixImax3;
          hxImin1=ixImin1-dixB*kr(1,3);hxImin2=ixImin2-dixB*kr(2,3)
          hxImin3=ixImin3-dixB*kr(3,3);hxImax1=ixImax1-dixB*kr(1,3)
          hxImax2=ixImax2-dixB*kr(2,3);hxImax3=ixImax3-dixB*kr(3,3);

          ! Calculate divergence and partial divergence
          Q=zero
          if (.not.zerodiv) call div_staggered(hxImin1,hxImin2,hxImin3,&
             hxImax1,hxImax2,hxImax3,s,Q(hxImin1:hxImax1,hxImin2:hxImax2,&
             hxImin3:hxImax3))
          select case(typeB(iw+b0_,iB))
          case("symm")
            ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iw)=zero
            do ix3=0,dixB-1
              call div_staggered(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                 ixImax3,s,Qp(ixImin1:ixImax1,ixImin2:ixImax2,&
                 ixImin3:ixImax3))
              ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3+ix3,iw)=&
              (Q(hxImin1:hxImax1,hxImin2:hxImax2,hxImax3-ix3)*s%geo%dvolume&
                 (hxImin1:hxImax1,hxImin2:hxImax2,hxImax3-ix3)&
             -Qp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3+ix3)*s%geo%dvolume&
                (ixImin1:ixImax1,ixImin2:ixImax2,ixImin3+ix3))&
              /s%geo%surfaceC3(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,&
                 ixIsmin3+ix3)
            end do
          case("asymm")
            ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iw)=zero
            do ix3=0,dixB-1
              call div_staggered(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                 ixImax3,s,Qp(ixImin1:ixImax1,ixImin2:ixImax2,&
                 ixImin3:ixImax3))
              ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3+ix3,iw)=&
             (-Q(hxImin1:hxImax1,hxImin2:hxImax2,hxImax3-ix3)*s%geo%dvolume&
                (hxImin1:hxImax1,hxImin2:hxImax2,hxImax3-ix3)&
             -Qp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3+ix3)*s%geo%dvolume&
                (ixImin1:ixImax1,ixImin2:ixImax2,ixImin3+ix3))&
              /s%geo%surfaceC3(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,&
                 ixIsmin3+ix3)
            end do
          case("cont", "zero")
            ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iw)=zero
            do ix3=0,dixB-1
              call div_staggered(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                 ixImax3,s,Qp(ixImin1:ixImax1,ixImin2:ixImax2,&
                 ixImin3:ixImax3))
              ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3+ix3,iw)=&
              (Q(hxImin1:hxImax1,hxImin2:hxImax2,hxImax3)*s%geo%dvolume&
                 (hxImin1:hxImax1,hxImin2:hxImax2,hxImax3)&
             -Qp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3+ix3)*s%geo%dvolume&
                (ixImin1:ixImax1,ixImin2:ixImax2,ixImin3+ix3))&
              /s%geo%surfaceC3(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,&
                 ixIsmin3+ix3)
            end do
          case("periodic")
          end select
        end if
      end do

      ! Fill cell averages
      call faces2centers(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,s)



      ! Apply on primitive variables?
      if (primitiveB(iside,idims)) then
         w(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) &
            = wsave(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) !restore from backup of the state in domain since p2c/c2p might lead to slight changes in w
         call conserve(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
            ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,w,x,patchfalse)
      end if
      
   else
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! minimal boundary
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      iB=ismin3
      ixImin1=ixBmin1;ixImin2=ixBmin2;ixImin3=ixBmin3;
      ixImax1=ixBmax1;ixImax2=ixBmax2;ixImax3=ixBmin3-1+dixB;

      ! Apply on primitive variables?
      if (primitiveB(iside,idims)) then         
         ixDmin1=ixBmin1;ixDmin2=ixBmin2;ixDmin3=ixBmin3+dixB;
         ixDmax1=ixBmax1;ixDmax2=ixBmax2;ixDmax3=ixBmax3-dixB;
         wsave(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) &
            = w(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) !make a backup of the state in domain since p2c/c2p might lead to slight changes in w
         call primitive(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
            ixDmin1,ixDmin2,ixDmin3,ixDmax1,ixDmax2,ixDmax3,w,x)
      end if
 
      ! cont/symm/asymm types
      do iw=1,nwflux+nwaux
 ! If the variable has no staggered counterpart
if ((iws(iw).lt.1).or.(iws(iw).gt.nws)) then
   
   select case (typeB(iw,iB))
   case ("symm")
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         = w(ixImin1:ixImax1,ixImin2:ixImax2,ixImax3+dixB:ixImax3+1:-1,iw)
   case ("asymm")
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         =-w(ixImin1:ixImax1,ixImin2:ixImax2,ixImax3+dixB:ixImax3+1:-1,iw)
   case ("cont")
      do ix3=ixImin3,ixImax3
         w(ixImin1:ixImax1,ixImin2:ixImax2,ix3,iw) = w(ixImin1:ixImax1,&
            ixImin2:ixImax2,ixImax3+1,iw)
      end do
   case("noinflow")
      if (iw==1+3)then
         do ix3=ixImin3,ixImax3
            w(ixImin1:ixImax1,ixImin2:ixImax2,ix3,iw) = min(w(ixImin1:ixImax1,&
               ixImin2:ixImax2,ixImax3+1,iw),zero)
         end do
      else
         do ix3=ixImin3,ixImax3
            w(ixImin1:ixImax1,ixImin2:ixImax2,ix3,iw) = w(ixImin1:ixImax1,&
               ixImin2:ixImax2,ixImax3+1,iw)
         end do
      end if
   case("limitinflow")
      if (iw==1+3)then
         do ix3=ixImin3,ixImax3
            w(ixImin1:ixImax1,ixImin2:ixImax2,ix3,iw) = min(w(ixImin1:ixImax1,&
               ixImin2:ixImax2,ixImax3+1,iw), &
                 w(ixImin1:ixImax1,ixImin2:ixImax2,ixImax3+1,iw)*ratebdflux)
         end do
      else
         do ix3=ixImin3,ixImax3
            w(ixImin1:ixImax1,ixImin2:ixImax2,ix3,iw) = w(ixImin1:ixImax1,&
               ixImin2:ixImax2,ixImax3+1,iw)
         end do
      end if
   case ("polefix")
      ! First overwrite cells in domain
      ixIpmin1=ixImin1;ixIpmin2=ixImin2;ixIpmin3=ixImax3+1;
      ixIpmax1=ixImax1;ixIpmax2=ixImax2;ixIpmax3=ixImax3+dixPolefix;
      do ix3=ixIpmin3,ixIpmax3
         w(ixIpmin1:ixIpmax1,ixIpmin2:ixIpmax2,ix3,iw) = w(ixIpmin1:ixIpmax1,&
            ixIpmin2:ixIpmax2,ixIpmax3+1,iw)
      end do
      ! Now apply symm boundary
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         = w(ixImin1:ixImax1,ixImin2:ixImax2,ixImax3+dixB:ixImax3+1:-1,iw)
   case ("apolefix")
      ! First interpolate cells in domain down to zero at pole
      ixIpmin1=ixImin1;ixIpmin2=ixImin2;ixIpmin3=ixImax3+1;
      ixIpmax1=ixImax1;ixIpmax2=ixImax2;ixIpmax3=ixImax3+dixPolefix;
      do ix3=ixIpmin3,ixIpmax3
         w(ixIpmin1:ixIpmax1,ixIpmin2:ixIpmax2,ix3,iw) = w(ixIpmin1:ixIpmax1,&
            ixIpmin2:ixIpmax2,ixIpmax3+1,iw) &
              * (x(ixIpmin1:ixIpmax1,ixIpmin2:ixIpmax2,ix3,3))&
                 /x(ixIpmin1:ixIpmax1,ixIpmin2:ixIpmax2,ixIpmax3+1,3)
      end do
      ! Now apply asymm boundary
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         =-w(ixImin1:ixImax1,ixImin2:ixImax2,ixImax3+dixB:ixImax3+1:-1,iw)
   case ("hotaka")
      ! First overwrite cells in domain with zero
      ixIpmin1=ixImin1;ixIpmin2=ixImin2;ixIpmin3=ixImax3+1;
      ixIpmax1=ixImax1;ixIpmax2=ixImax2;ixIpmax3=ixImax3+dixHotaka;
      do ix3=ixIpmin3,ixIpmax3
         w(ixIpmin1:ixIpmax1,ixIpmin2:ixIpmax2,ix3,iw) = zero
      end do
      ! Now apply symm boundary
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         = w(ixImin1:ixImax1,ixImin2:ixImax2,ixImax3+dixB:ixImax3+1:-1,iw)
   case ("ahotaka")
      ! First overwrite cells in domain with zero
      ixIpmin1=ixImin1;ixIpmin2=ixImin2;ixIpmin3=ixImax3+1;
      ixIpmax1=ixImax1;ixIpmax2=ixImax2;ixIpmax3=ixImax3+dixHotaka;
      do ix3=ixIpmin3,ixIpmax3
         w(ixIpmin1:ixIpmax1,ixIpmin2:ixIpmax2,ix3,iw) = zero
      end do
      ! Now apply asymm boundary
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         =-w(ixImin1:ixImax1,ixImin2:ixImax2,ixImax3+dixB:ixImax3+1:-1,iw)
   case ("special")
      ! skip it here, do AFTER all normal type boundaries are set
   case ("aperiodic")
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         = - w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw)
   case ("periodic")
      !            call mpistop("periodic bc info should come from neighbors")
   case default
      write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
           "for variable iw=",iw," and side iB=",iB
   end select


else !If there is a staggered counterpart
   ! At this stage, extrapolation is applied only to the tangential components
   idir=mod(iws(iw)-1,ndir)+1 ! vector direction
   if (idir.ne.3) then
      ixIsmax1=ixImax1;ixIsmax2=ixImax2;ixIsmax3=ixImax3;
      ixIsmin1=ixImin1-kr(1,idir);ixIsmin2=ixImin2-kr(2,idir)
      ixIsmin3=ixImin3-kr(3,idir);
      select case(typeB(iw,iB))
      case ("symm")
         ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iws(iw)) &
            = ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,&
            ixIsmax3+dixB:ixIsmax3+1:-1,iws(iw))
      case ("asymm")
         ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iws(iw)) &
            =-ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,&
            ixIsmax3+dixB:ixIsmax3+1:-1,iws(iw))
      case ("cont")
         do ix3=ixIsmin3,ixIsmax3
            ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ix3,iws(iw)) &
               = ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmax3+1,iws(iw))
         end do
      case("periodic")
      case ("zero")
         ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iws(iw)) &
            = 0.0d0
      case ("specialT")
      call specialbound_usr(time,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
         ixGmax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,-iB,s)
      case ("special")
         ! skip it here, do AFTER all normal type boundaries are set
      case default
         write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
              "for variable iw=",iw," and side iB=",iB
      end select


      ! Special treatment when we have coarse neighbors in tangential dir:
      do is=1,2
         if (is_coarse(idir,is)) then

            ixIsmin1=ixImin1-kr(1,idir);ixIsmin2=ixImin2-kr(2,idir)
            ixIsmin3=ixImin3-kr(3,idir);
            ixIsmax1=ixImax1;ixIsmax2=ixImax2;ixIsmax3=ixImax3;

            if (is.eq.1) then ! coarse neighbor at beginning of block

               select case(idir)
                  
                  case (1)
                     ixIsmax1=min(ixGlo1+dixB-1,ixIsmax1);
                     surf=>s%geo%surfaceC1
                  

                  case (2)
                     ixIsmax2=min(ixGlo2+dixB-1,ixIsmax2);
                     surf=>s%geo%surfaceC2
                  

                  case (3)
                     ixIsmax3=min(ixGlo3+dixB-1,ixIsmax3);
                     surf=>s%geo%surfaceC3
                  ! I am hating LASY just right now.
               end select

            else ! coarse neighbor at end of block

               select case(idir)
                  
                  case (1)
                     ixIsmin1=max(ixGhi1-dixB,ixIsmin1);
                     surf=>s%geo%surfaceC1
                  

                  case (2)
                     ixIsmin2=max(ixGhi2-dixB,ixIsmin2);
                     surf=>s%geo%surfaceC2
                  

                  case (3)
                     ixIsmin3=max(ixGhi3-dixB,ixIsmin3);
                     surf=>s%geo%surfaceC3
                  ! I am hating LASY just right now.
               end select

            end if


            select case(typeB(iw,iB))
            case ("cont")
               do ix3=ixIsmin3,ixIsmax3
                  where ( (surf(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,&
                     ixIsmax3+1) + surf(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,&
                     ixIsmax3+2)) .gt. smalldouble )
                     ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ix3,iws(iw)) &
                        = (ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmax3+1,&
                        iws(iw))*surf(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,&
                        ixIsmax3+1) &
                          + ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmax3+2,&
                             iws(iw))*surf(ixIsmin1:ixIsmax1,&
                             ixIsmin2:ixIsmax2,ixIsmax3+2) ) &
                             / (surf(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,&
                             ixIsmax3+1)+surf(ixIsmin1:ixIsmax1,&
                             ixIsmin2:ixIsmax2,ixIsmax3+2))
                  elsewhere
                     ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ix3,iws(iw)) &
                        = zero
                  end where
               end do

            case default
 !Only implemented cont type, symm and asymm will also work if problem respects symmetry
               ! 3D case will need to average four cells together!!!
            end select

         end if
      end do

   end if

end if

end do


      ! Now that the tangential components are set,
      ! fill the normal components using a prescription for the divergence.
      ! This prescription is given by the typeB for the normal component.
do iw=1,nws
        ! When using more than one staggered vector field
        idir=mod(iw-1,ndir)+1
        ! Consider only normal direction
        if (idir.eq.3) then
          ixIsmin1=ixImin1-kr(1,3);ixIsmin2=ixImin2-kr(2,3)
          ixIsmin3=ixImin3-kr(3,3);ixIsmax1=ixImax1-kr(1,3)
          ixIsmax2=ixImax2-kr(2,3);ixIsmax3=ixImax3-kr(3,3);
          jxImin1=ixImin1+dixB*kr(1,3);jxImin2=ixImin2+dixB*kr(2,3)
          jxImin3=ixImin3+dixB*kr(3,3);jxImax1=ixImax1+dixB*kr(1,3)
          jxImax2=ixImax2+dixB*kr(2,3);jxImax3=ixImax3+dixB*kr(3,3);

          ! Calculate divergence and partial divergence
          Q=0
          if (.not.zerodiv) call div_staggered(jxImin1,jxImin2,jxImin3,&
             jxImax1,jxImax2,jxImax3,s,Q(jxImin1:jxImax1,jxImin2:jxImax2,&
             jxImin3:jxImax3))
          select case(typeB(iw+b0_,iB))
          case("symm")
            ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iw)=zero
            do ix3=0,dixB-1
              call div_staggered(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                 ixImax3,s,Qp(ixImin1:ixImax1,ixImin2:ixImax2,&
                 ixImin3:ixImax3))
              ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmax3-ix3,iw)=&
             -(Q(jxImin1:jxImax1,jxImin2:jxImax2,jxImin3+ix3)*s%geo%dvolume&
                (jxImin1:jxImax1,jxImin2:jxImax2,jxImin3+ix3)&
             -Qp(ixImin1:ixImax1,ixImin2:ixImax2,ixImax3-ix3)*s%geo%dvolume&
                (ixImin1:ixImax1,ixImin2:ixImax2,ixImax3-ix3))&
             /s%geo%surfaceC3(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,&
                ixIsmax3-ix3)
            end do
          case("asymm")
             ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iw)=zero
            do ix3=0,dixB-1
              call div_staggered(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                 ixImax3,s,Qp(ixImin1:ixImax1,ixImin2:ixImax2,&
                 ixImin3:ixImax3))
              ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmax3-ix3,iw)=&
            -(-Q(jxImin1:jxImax1,jxImin2:jxImax2,jxImin3+ix3)*s%geo%dvolume&
               (jxImin1:jxImax1,jxImin2:jxImax2,jxImin3+ix3)&
             -Qp(ixImin1:ixImax1,ixImin2:ixImax2,ixImax3-ix3)*s%geo%dvolume&
                (ixImin1:ixImax1,ixImin2:ixImax2,ixImax3-ix3))&
             /s%geo%surfaceC3(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,&
                ixIsmax3-ix3)
            end do
          case("cont", "zero")
            ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,iw)=zero
            do ix3=0,dixB-1
              call div_staggered(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                 ixImax3,s,Qp(ixImin1:ixImax1,ixImin2:ixImax2,&
                 ixImin3:ixImax3))
              ws(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmax3-ix3,iw)=&
             -(Q(jxImin1:jxImax1,jxImin2:jxImax2,jxImin3)*s%geo%dvolume&
                (jxImin1:jxImax1,jxImin2:jxImax2,jxImin3)&
             -Qp(ixImin1:ixImax1,ixImin2:ixImax2,ixImax3-ix3)*s%geo%dvolume&
                (ixImin1:ixImax1,ixImin2:ixImax2,ixImax3-ix3))&
              /s%geo%surfaceC3(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,&
                 ixIsmax3-ix3)
            end do
          case("periodic")
          end select
        end if
      end do

      ! Fill cell averages
      call faces2centers(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,s)

      
      ! Apply on primitive variables?
      if (primitiveB(iside,idims)) then
         w(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) &
            = wsave(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) !restore from backup of the state in domain since p2c/c2p might lead to slight changes in w
         call conserve(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
            ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,w,x,patchfalse)
      end if
      
   end if 
end select


! do special case AFTER all normal cases are set
if (any(typeB(1:nwflux+nwaux,iB)=="special")) then
   call specialbound_usr(time,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
      ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iB,s)
end if

   end associate
end subroutine bc_phys
!=============================================================================
subroutine getintbc(time,psuse)

use mod_amrvacdef

double precision, intent(in)              :: time
type(state), dimension(ngridshi)          :: psuse

! .. local ..
integer :: iigrid, igrid, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
   level
!----------------------------------------------------------------------------
ixOmin1=ixGlo1+dixB;ixOmin2=ixGlo2+dixB;ixOmin3=ixGlo3+dixB
ixOmax1=ixGhi1-dixB;ixOmax2=ixGhi2-dixB;ixOmax3=ixGhi3-dixB;

!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid,level)
do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
   call set_tmpGlobals(igrid)
   level = node(plevel_,igrid)
   call bc_int(level,time,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixOmin1,&
      ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,psuse(igrid)%w%w,px(igrid)%x)
end do
!$OMP END PARALLEL DO

      
end subroutine getintbc
!=============================================================================
subroutine is_neighbor_coarse(s,is_coarse)

use mod_amrvacdef

type(state), intent(in)          :: s
logical, intent(out)             :: is_coarse(ndim,2)
! .. local ..
integer                          :: i1,i2,i3, i
!-----------------------------------------------------------------------------

is_coarse=.false.
! If we are the coarse version of a buffer
! dont consider neighbor as coarser
if (s%is_coarse .eqv. .true.) return 

   do i3=-1,1
   do i2=-1,1
   do i1=-1,1
      if ((i1.eq.0.and.i2.eq.0.and.i3.eq.0).or.(abs(i1)+abs(i2)+abs&
         (i3)).ne.1) cycle
      if (neighbor_type(i1,i2,i3,s%igrid).eq.2) then
         if (i1.eq.-1) is_coarse(1,1) = .true.
         if (i2.eq.-1) is_coarse(2,1) = .true.
         if (i3.eq.-1) is_coarse(3,1) = .true.
         if (i1.eq.+1) is_coarse(1,2) = .true.
         if (i2.eq.+1) is_coarse(2,2) = .true.
         if (i3.eq.+1) is_coarse(3,2) = .true.
      end if
   end do
   end do
   end do

end subroutine is_neighbor_coarse
!=============================================================================
