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
module mod_odeint

  implicit none

! This module packages odeint from numerical recipes.  
  integer                  :: MAXSTP, NMAX, KMAXX
  PARAMETER (MAXSTP=10000,NMAX=50,KMAXX=200)
  double precision         :: TINY
  PARAMETER (TINY=1.d-30)

contains
!=============================================================================
      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,rkqs)

      INTEGER nbad,nok,nvar
      double precision eps,h1,hmin,x1,x2,ystart(nvar)
      INTEGER i,nstp
      double precision h,hdid,hnext,x,xsav,dydx(NMAX),y(NMAX),yscal(NMAX)
      integer kmax, kount
      double precision dxsav, yp(NMAX,KMAXX), xp(KMAXX)
      EXTERNAL derivs,rkqs

      COMMON /path/ kmax,kount,dxsav,xp,yp
  
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      do 11 i=1,nvar
         y(i)=ystart(i)
 11   continue
      if (kmax.gt.0) xsav=x-2.*dxsav
      do 16 nstp=1,MAXSTP
         call derivs(x,y,dydx)
         do 12 i=1,nvar
            yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
 12      continue
         if(kmax.gt.0)then
            if(abs(x-xsav).gt.abs(dxsav)) then
               if(kount.lt.kmax-1)then
                  kount=kount+1
                  xp(kount)=x
                  do 13 i=1,nvar
                     yp(i,kount)=y(i)
 13               continue
                  xsav=x
               endif
            endif
         endif
         if((x+h-x2)*(x+h-x1).gt.0.d0) h=x2-x
         call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
         if(hdid.eq.h)then
            nok=nok+1
         else
            nbad=nbad+1
         endif
         if((x-x2)*(x2-x1).ge.0.d0)then
            do 14 i=1,nvar
               ystart(i)=y(i)
 14         continue
            if(kmax.ne.0)then
               kount=kount+1
               xp(kount)=x
               do 15 i=1,nvar
                  yp(i,kount)=y(i)
 15            continue
            endif
            return
         endif
         if(abs(hnext).lt.hmin) call stop("stepsize smaller than minimum in &
            odeint")
         h=hnext
 16   continue
      call stop( 'too many steps in odeint')
      return
   END subroutine odeint
!=============================================================================
      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)

      INTEGER n,NMAX
      double precision eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
!     CU USES derivs,rkck
      INTEGER i
      double precision errmax,h,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW,&
         PSHRNK,ERRCON
      PARAMETER (SAFETY=0.9d0,PGROW=-.2d0,PSHRNK=-.25d0,ERRCON=1.89d-4)

      h=htry
 1    call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
      errmax=0.d0
      do 11 i=1,n
         errmax=max(errmax,abs(yerr(i)/yscal(i)))
 11   continue
      errmax=errmax/eps
      if(errmax.gt.1.d0)then
         h=SAFETY*h*(errmax**PSHRNK)
         if(h.lt.0.1d0*h)then
            h=.1d0*h
         endif
         xnew=x+h
         if(xnew.eq.x)call stop( 'stepsize underflow in rkqs')
         goto 1
      else
         if(errmax.gt.ERRCON)then
            hnext=SAFETY*h*(errmax**PGROW)
         else
            hnext=5.d0*h
         endif
         hdid=h
         x=x+h
         do 12 i=1,n
            y(i)=ytemp(i)
 12      continue
         return
      endif
      END subroutine rkqs
!C (C) Copr. 1986-92 Numerical Recipes Software Vs1&v%1jw#<?4210(9p#.
!=============================================================================
      SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)

      INTEGER n,NMAX
      double precision  h,x,dydx(n),y(n),yerr(n),yout(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
!CU    USES derivs
      INTEGER i
      double precision ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX),&
         ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51, B52,B53,B54,&
         B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
      PARAMETER (A2=.2d0,A3=.3d0,A4=.6d0,A5=1.d0, A6=.875d0,B21=.2d0,B31=3.d0&
         /40.d0,B32=9.d0/40.d0,B41=.3d0,B42=-.9d0,B43=1.2d0,B51=-11.d0&
         /54.d0,B52=2.5d0,B53=-70.d0/27.d0,B54=35.d0/27.d0,B61=1631.d0&
         /55296.d0, B62=175.d0/512.d0, B63=575.d0/13824.d0,B64=44275.d0&
         /110592.d0, B65=253.d0/4096.d0,C1=37.d0/378.d0,C3=250.d0&
         /621.d0,C4=125.d0/594.d0,C6=512.d0/1771.d0,DC1=C1-2825.d0&
         /27648.d0,DC3=C3-18575.d0/48384.d0,DC4=C4-13525.d0&
         /55296.d0,DC5=-277.d0/14336.d0,DC6=C6-.25d0)
      do 11 i=1,n
         ytemp(i)=y(i)+B21*h*dydx(i)
 11   continue
      call derivs(x+A2*h,ytemp,ak2)
      do 12 i=1,n
         ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
 12   continue
      call derivs(x+A3*h,ytemp,ak3)
      do 13 i=1,n
         ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
 13   continue
      call derivs(x+A4*h,ytemp,ak4)
      do 14 i=1,n
         ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
 14   continue
      call derivs(x+A5*h,ytemp,ak5)
      do 15 i=1,n
         ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4&
            (i)+B65*ak5(i))
 15   continue
      call derivs(x+A6*h,ytemp,ak6)
      do 16 i=1,n
         yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
 16   continue
      do 17 i=1,n
         yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*ak6(i))
 17   continue
      return
      END subroutine rkck
!C  (C) Copr. 1986-92 Numerical Recipes Software Vs1&v%1jw#<?4210(9p#.
!=============================================================================
      subroutine stop(text)

        character(len=*), intent(in)   :: text
!-----------------------------------------------------------------------------
        
        print*, text
        STOP
      end subroutine stop
!=============================================================================
end module mod_odeint
!=============================================================================
