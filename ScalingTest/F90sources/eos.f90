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

!###########################################################################
! module amrvacphys - srmhdeos version april 2009
! This module is developed using the paper Meliani et al 2004
!===========================================================================
subroutine Enthalpy(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
   ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,patchw,rhoh)

!================== IMPORTANT ==================!
!This subroutine is used only with primitive variables on input w
!===============================================!

use mod_amrvacdef

integer, intent(in)                                    :: ixImin1,ixImin2,&
   ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
   ixOmax3
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nw), intent(in)    :: w
logical,          dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
    intent(in)         :: patchw
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
    intent(out)        :: rhoh

!--------------------------------------! 




where(.not.patchw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
 rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
    ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_) + govergminone * w(ixOmin1:ixOmax1,&
    ixOmin2:ixOmax2,ixOmin3:ixOmax3,pp_)
end where




end subroutine Enthalpy
!=============================================================================
subroutine Pressuren(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,varconserve,p,patchw)

!================== IMPORTANT ==================!
!This subroutine is used with conserved variables in w when varconserve=T
!This subroutine is used with primitive variables in w when varconserve=F
!   both cases assume updated auxiliary variables xi_ en lfac_
!===============================================!

use mod_amrvacdef

integer,intent(in)                                :: ixImin1,ixImin2,ixImin3,&
   ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(in),dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw):: w
logical,intent(in)                                :: varconserve
logical,intent(in),dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)       &
           :: patchw
double precision,intent(out), dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)    :: p

!--------------------------------------! 



if (varconserve) then
 where(.not.patchw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
    p(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = ( w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_)/(w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)**2) - w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,d_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,lfac_) ) / govergminone
 end where
else
 where(.not.patchw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
    p(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = ( w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_)/(w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)**2) - w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_) ) / govergminone
 end where




end if

return
end subroutine Pressuren
!============================================================================
subroutine Pressure(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
   ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,varconserve,p)

!================== IMPORTANT ==================!
!This subroutine is used with conserved variables in w when varconserve=T
!This subroutine is used with primitive variables in w when varconserve=F
!   both cases assume updated auxiliary variables xi_ en lfac_
!===============================================!

use mod_amrvacdef

integer,intent(in)                                :: ixImin1,ixImin2,ixImin3,&
   ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(in),dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw):: w
logical,intent(in)                                :: varconserve
double precision,intent(out), dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)    :: p

!--------------------------------------! 



if (varconserve) then
    p(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = (w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_)/(w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)**2)-(w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,d_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,lfac_)))/govergminone
else
    p(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = (w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_)/(w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)**2)-w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_))/govergminone
end if




return
end subroutine Pressure
!============================================================================
subroutine getinvcsound2(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,invcsound2)

!================== IMPORTANT ==================!
!This subroutine is used with conserved variables in w
!===============================================!

use mod_amrvacdef

integer,intent(in)                                    :: ixImin1,ixImin2,&
   ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
   ixOmax3
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nw), intent(in)   :: w
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
    intent(out)       :: invcsound2
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)                    :: h

!--------------------------------------!  

h(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_)/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   ixOmin3:ixOmax3,lfac_)**2)



invcsound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
   =  (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_)&
   /(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_) - &
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)*w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,ixOmin3:ixOmax3,d_))) &
               /(eqpar(gamma_)-one)




end subroutine getinvcsound2
!=============================================================================
subroutine getcsound2(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,varconserve,csound2)

!================== IMPORTANT ==================!
!This subroutine is used with conserved variables in w when varconserve=T
!This subroutine is used with primitive variables in w when varconserve=F
!===============================================!

use mod_amrvacdef

integer,intent(in)                                   :: ixImin1,ixImin2,&
   ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
   ixOmax3
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nw), intent(in)  :: w
logical, intent(in)                                  :: varconserve
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
   intent(out)       :: csound2
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)                   :: h, P

!--------------------------------------!

h(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_)/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   ixOmin3:ixOmax3,lfac_)**2)


call Pressure(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
   ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,varconserve,P)
csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = (eqpar(gamma_) * &
   P(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)) / h(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,ixOmin3:ixOmax3)



end subroutine getcsound2
!=============================================================================
subroutine FuncPressure(xicurrent,lfac,d,ssqr,tau,dlfacdxi,p,dpdxi,dpdx)

! compute pointwise value for pressure p and dpdxi
!
! xi := lfac**2 rho*h
! x  := v**2
! chi := rho*h - rho ??? = rho*(h-1) ???  
!  
use mod_amrvacdef
  
double precision, intent(in)         :: xicurrent,lfac,d,ssqr,tau,dlfacdxi
double precision, intent(out)        :: p,dpdxi,dpdx
! .. local ..
double precision                     :: rho,h

double precision                     :: dpdchi, dpdrho 


!_______________________________________________________________!


h=xicurrent/(lfac**2)  ! h is actually rho*h
rho=d/lfac




p = (h - rho)/govergminone
dpdchi = one/govergminone
dpdrho = zero
dpdx   = (half*lfac*d-xicurrent)/govergminone
dpdxi = dpdchi * one/lfac**2
if (dlfacdxi .ne. 0.0d0) &
     dpdxi = dpdxi  + dpdchi * ((d*lfac-2.0d0*xicurrent)/lfac**3) * dlfacdxi




end subroutine FuncPressure
!=============================================================================!
subroutine smallvaluesEOS
! This is the smallvalues for Synge gas to be precise.  
use mod_amrvacdef

double precision::Lsmallrho,Lsmallp,LsmallE
!--------------------------------------------------
Lsmallrho=(1.0d0 + 10.0d0 * minrho) * minrho
Lsmallp=(1.0d0 + 10.0d0 * minp) * minp
LsmallE=Lsmallp/(eqpar(gamma_)-one)+dsqrt((Lsmallp/(eqpar(gamma_)-one))&
   **2+Lsmallrho**2)

smalltau=half*((eqpar(gamma_)+one)*LsmallE-(eqpar(gamma_)-one)*Lsmallrho**2&
   /LsmallE)-Lsmallp-Lsmallrho
! may need to replace by smallxi above

smallxi=half*((eqpar(gamma_)+one)*LsmallE-(eqpar(gamma_)-one)*Lsmallrho**2&
   /LsmallE)

end Subroutine smallvaluesEOS
!=============================================================================!
subroutine xinoFlow(xi,tau,d,bb)

use mod_amrvacdef

double precision:: xi,tau,d,bb
!_______________________________________________________________!



xi = eqpar(gamma_) *(tau - half * bb)+d


end subroutine xinoFlow

!=============================================================================
subroutine FuncEnthalpy(pcurrent,lfac2inv,d,tau,sqrs,xicurrent,dv2d2p,h,dhdp,&
   ierror)

use mod_amrvacdef
  
double precision, intent(in)             :: pcurrent,lfac2inv,d,tau,sqrs,&
   xicurrent,dv2d2p
double precision, intent(out)            :: h,dhdp
integer, intent(inout)                   :: ierror

double precision                         :: rho

!_______________________________________________________________!
rho=d*dsqrt(lfac2inv)




h = rho + govergminone * pcurrent
dhdp = govergminone + d/sqrt(lfac2inv)*sqrs/xicurrent**3



return
end subroutine FuncEnthalpy
!=========================================================================
subroutine Bisection_Enthalpy(pcurrent,lfac2inv,d,tau,sqrs,xicurrent,h,ierror)

use mod_amrvacdef
  
integer:: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,&
   ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision:: pcurrent,lfac2inv,d,tau,sqrs,xicurrent
double precision:: h
integer::ierror
double precision:: rho

!_______________________________________________________________!
rho=d*dsqrt(lfac2inv)



h = rho + govergminone * pcurrent


return
end subroutine Bisection_Enthalpy

!=============================================================================
subroutine Calcule_Geffect(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,varconserve,Geffect)

!================== IMPORTANT ==================!
!This subroutine is used with conserved variables in w when varconserve=T
!This subroutine is used with primitive variables in w when varconserve=F
!   both cases assume updated auxiliary variables xi_ en lfac_
!===============================================!

use mod_amrvacdef

integer:: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,&
   ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision,intent(in):: w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw)
logical,intent(in)         :: varconserve
double precision,intent(out):: Geffect(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)

!_______________________________________________________________!


Geffect(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = eqpar(gamma_)



end subroutine Calcule_Geffect
!=============================================================================
! end module amrvacphys - srmhdeos
!##############################################################################
