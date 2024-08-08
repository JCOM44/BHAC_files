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
subroutine getaux(clipping,w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,subname)

! Calculate the auxiliary variables lfac and xi within ixO^L
! clipping is not used (yet) 

use mod_metric, only: raise3, square3u
use mod_con2prim, only: con2prim
use mod_amrvacdef

logical, intent(in)            :: clipping
integer, intent(in)            :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(inout):: w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw)
double precision, intent(in)   :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndim)
character(len=*), intent(in)   :: subname
! .. local ..
integer::          err,ix1,ix2,ix3, i1,i2,i3
double precision :: sold1_,sold2_,sold3_,bold1_,bold2_,bold3_, lfacold, xiold
integer          :: patchierror(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3) :: ssqr, bsqr, sdotb2
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:3) :: sU

!-----------------------------------------------------------------------------


! MM, MB2, BB:
call raise3(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,&
   ixOmin3,ixOmax1,ixOmax2,ixOmax3,myM,w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,s1_:s3_),sU)
call square3u(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,&
   ixOmin3,ixOmax1,ixOmax2,ixOmax3,myM,w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,b1_:b3_),bsqr)

ssqr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)   =  w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,ixOmin3:ixOmax3,s0_+1)*sU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   ixOmin3:ixOmax3,1)+ w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
   s0_+2)*sU(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
   2)+ w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,s0_+3)*sU&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)
sdotb2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = ( w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+1)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   ixOmin3:ixOmax3,s0_+1)+ w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
   b0_+2)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
   s0_+2)+ w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
   b0_+3)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,s0_+3))**2

! we compute auxiliaries lfac,xi from D,S,tau,B
! put the lfac and xi in the auxiliary fields lfac_ and xi_
do ix3= ixOmin3,ixOmax3
do ix2= ixOmin2,ixOmax2
do ix1= ixOmin1,ixOmax1
    lfacold = w(ix1,ix2,ix3,lfac_)
    xiold = w(ix1,ix2,ix3,xi_)
     sold1_=w(ix1,ix2,ix3,s1_); sold2_=w(ix1,ix2,ix3,s2_)
     sold3_=w(ix1,ix2,ix3,s3_);
     bold1_=w(ix1,ix2,ix3,b1_); bold2_=w(ix1,ix2,ix3,b2_)
     bold3_=w(ix1,ix2,ix3,b3_);

    call con2prim(w(ix1,ix2,ix3,lfac_),w(ix1,ix2,ix3,xi_), &
             w(ix1,ix2,ix3,d_),sdotb2(ix1,ix2,ix3),ssqr(ix1,ix2,ix3),bsqr(ix1,&
                ix2,ix3),w(ix1,ix2,ix3,tau_),err)



    
    patchierror(ix1,ix2,ix3)=err
    if (err/=0.and.strictgetaux) then
       write(*,*) 'Getaux error:',err,'ix^D=',ix1,ix2,ix3
       write(*,*) 'input lfac=',lfacold,'s=',sold1_,sold2_,sold3_,'xi=',xiold,&
          'b=',bold1_,bold2_,bold3_
       write(*,*) 'New ','lfac=',w(ix1,ix2,ix3,lfac_),'s=',w(ix1,ix2,ix3,s1_),&
          w(ix1,ix2,ix3,s2_),w(ix1,ix2,ix3,s3_),'xi=',w(ix1,ix2,ix3,xi_)
       write(*,*) 'B=',w(ix1,ix2,ix3,b1_),w(ix1,ix2,ix3,b2_),w(ix1,ix2,ix3,&
          b3_)
       write(*,*) 'D=',w(ix1,ix2,ix3,d_), 'tau=',w(ix1,ix2,ix3,tau_)
       write(*,*) 'igrid =', saveigrid 
       write(*,*) 'level:', node(plevel_,saveigrid)
       write(*,*) 'position  ', x(ix1,ix2,ix3, 1:ndim)
       write(*,*) 'Called from: ',subname
       call mpistop("problem in getaux: retry with strictgetaux=F")
    end if
enddo
enddo
enddo


! first try user-defined method:
if(.not.strictgetaux.and.any(patchierror(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   ixOmin3:ixOmax3)/=0)) call correctaux_usr(ixImin1,ixImin2,ixImin3,ixImax1,&
   ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,&
   patchierror,subname)
! then, if any patchierror was not reset to 0, use default:
if(.not.strictgetaux.and.any(patchierror(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   ixOmin3:ixOmax3)/=0)) call correctaux(ixImin1,ixImin2,ixImin3,ixImax1,&
   ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,&
   patchierror,subname)

end subroutine getaux
!=============================================================================
subroutine correctaux(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
   ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,patchierror,subname)

use mod_amrvacdef

integer, intent(in)            :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
character(len=*), intent(in)   :: subname
double precision, intent(inout):: w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw)
double precision, intent(in)   :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndim)
integer, intent(inout)         :: patchierror(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
   ixGlo3:ixGhi3)

integer        :: iw, kxOmin1,kxOmin2,kxOmin3,kxOmax1,kxOmax2,kxOmax3, ix1,&
   ix2,ix3, i
logical        :: patchw(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
!-----------------------------------------------------------------------------

do ix3= ixOmin3,ixOmax3
do ix2= ixOmin2,ixOmax2
do ix1= ixOmin1,ixOmax1
    if (patchierror(ix1,ix2,ix3)/=0) then
        do i=1,nflatgetaux
           kxOmin1= max(ix1-i,ixOmin1);
            kxOmax1= min(ix1+i,ixOmax1);
           kxOmin2= max(ix2-i,ixOmin2);
            kxOmax2= min(ix2+i,ixOmax2);
           kxOmin3= max(ix3-i,ixOmin3);
            kxOmax3= min(ix3+i,ixOmax3);
           if (any(patchierror(kxOmin1:kxOmax1,kxOmin2:kxOmax2,&
              kxOmin3:kxOmax3)==0)) exit
        end do
        if (any(patchierror(kxOmin1:kxOmax1,kxOmin2:kxOmax2,kxOmin3:kxOmax3)&
           ==0))then
            ! in contrast to srhd case: always switch to primitive
            ! as error can become large when this is not done
            patchw(kxOmin1:kxOmax1,kxOmin2:kxOmax2,kxOmin3:kxOmax3)&
               =(patchierror(kxOmin1:kxOmax1,kxOmin2:kxOmax2,kxOmin3:kxOmax3)&
               /=0)
            call primitiven(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
               kxOmin1,kxOmin2,kxOmin3,kxOmax1,kxOmax2,kxOmax3,w,patchw)
            do iw = 1,nwflux+nwaux
               ! in contrast to srhd: do not alter magnetic field
               if (iw/=b1_.and.iw/=b2_.and.iw/=b3_) then
                   w(ix1,ix2,ix3,iw)=sum(w(kxOmin1:kxOmax1,kxOmin2:kxOmax2,&
                      kxOmin3:kxOmax3,iw),patchierror(kxOmin1:kxOmax1,&
                      kxOmin2:kxOmax2,kxOmin3:kxOmax3)==0)/count(patchierror&
                      (kxOmin1:kxOmax1,kxOmin2:kxOmax2,kxOmin3:kxOmax3)==0)
               end if  
            end do
            ! Check if the interpolated variables obey floors:
            if (tlow>zero) call fixp_usr(ixImin1,ixImin2,ixImin3,ixImax1,&
               ixImax2,ixImax3,ix1,ix2,ix3,ix1,ix2,ix3,w,x)
            patchw(ix1,ix2,ix3)=.false.
            call conserven(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
               kxOmin1,kxOmin2,kxOmin3,kxOmax1,kxOmax2,kxOmax3,w,patchw)
            patchierror(ix1,ix2,ix3) = 0
        else
           ! no cells without error were found in cube of size nflatgetaux
           ! --> point of no recovery
           write(*,*) 'Getaux error:',patchierror(ix1,ix2,ix3),'ix^D=',ix1,&
              ix2,ix3
           write(*,*) 'it=', it, 't=', t
           write(*,*) 'Called from: ',subname
           write(*,*) 'New ','d=',w(ix1,ix2,ix3,d_),'s=',w(ix1,ix2,ix3,s1_),&
              w(ix1,ix2,ix3,s2_),w(ix1,ix2,ix3,s3_),'tau=',w(ix1,ix2,ix3,tau_)
           write(*,*) 'B=',w(ix1,ix2,ix3,b1_),w(ix1,ix2,ix3,b2_),w(ix1,ix2,&
              ix3,b3_)
           write(*,*) 'lfac=',w(ix1,ix2,ix3,lfac_),'xi=',w(ix1,ix2,ix3,xi_)


           write(*,*) 'position  ', x(ix1,ix2,ix3, 1:ndim)
           call mpistop("---------correctaux from getaux----------")
        end if
     end if
enddo
enddo
enddo

end subroutine correctaux
!=============================================================================
! end module correctaux.t
!=============================================================================
