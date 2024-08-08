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

module mod_con2prim_vars

  implicit none
  save

  double precision       :: d,ssqr,tau,bsqr,sdotb2
  !$OMP THREADPRIVATE(d,ssqr,tau,bsqr,sdotb2)
  
end module mod_con2prim_vars
!=============================================================================
module mod_con2prim

  implicit none
  save

  double precision       :: lfacmax
  double precision       :: xi_guess, lfac_guess
  double precision       :: xil, xih, dplus
  double precision       :: lfac, xi
  integer                :: ierror, typeinversion_

  logical                :: use_entropy
 !$OMP THREADPRIVATE(lfacmax,xi_guess,lfac_guess,xil,xih,dplus,lfac,xi,ierror,typeinversion_,use_entropy)

  integer, parameter     :: t1d_=1, t2d_=2, t1d2d_=3, t2d1d_=4, t2d1dentropy_&
     =5, t1d2dentropy_=6, tentropy_=7, t1dentropy_=8, t2dentropy_=9

contains


  !=============================================================================
  subroutine setup_con2prim

    use mod_amrvacdef

    !-----------------------------------------------------------------------------
    select case (typeinversion)
    case('1DW', '1dw')
       typeinversion_ = t1d_

    case('2DW', '2dw')
       typeinversion_ = t2d_

    case('1D2D')
       typeinversion_ = t1d2d_

    case('2D1D')
       typeinversion_ = t2d1d_

    case('2D1DEntropy')
        call mpistop("Set #ifdef ENTROPY flag when you want to use Entropy &
           inversion") 
       

    case('1D2DEntropy')
        call mpistop("Set #ifdef ENTROPY flag when you want to use Entropy &
           inversion") 
       

    case('1DEntropy')
        call mpistop("Set #ifdef ENTROPY flag when you want to use Entropy &
           inversion") 
       

    case('2DEntropy')
        call mpistop("Set #ifdef ENTROPY flag when you want to use Entropy &
           inversion") 
       

    case('Entropy')
        call mpistop("Set #ifdef ENTROPY flag when you want to use Entropy &
           inversion") 
       

    case default
       call mpistop('Unknown typeinversion in setup_con2prim')
    end select

    ! Should we use the entropy for backup?
    use_entropy = index(typeinversion,'Entropy')>=1

    ! the maximal allowed Lorentz factor
    lfacmax=one/sqrt(one-(one-dmaxvel)**2)

  end subroutine setup_con2prim
  !=============================================================================
  subroutine con2prim(lfac,xi,myd,mysdotb2,myssqr,mybsqr,mytau,ierror)

    ! (D,S,tau,B) --> compute auxiliaries lfac and xi

    use mod_con2prim_vars
    use mod_amrvacdef

    double precision, intent(in)    :: myd
    double precision, intent(in)    :: mysdotb2, myssqr, mybsqr
    double precision, intent(inout) :: lfac, xi, mytau !taufix changes tau acc. to entropy
    integer, intent(inout)          :: ierror

    ! .. local ..

    double precision:: xi1,xi2,xii,dxi,f,df
    double precision:: temp,vsqr,fl,fh,lfacl,lfach,lfaci
    double precision:: er,er1,xiprev
    double precision:: fv,dfv
    double precision:: lastxiok,lastlfacok,lastf
    double precision:: dlfac,p,dpdxi,dpdx
    integer         :: ni,nit,niiter
    logical         :: finished,tests,testsl,testsh
    !-----------------------------------------------------------------------------

    ! Save the input-state in mod_con2prim_vars :
    d = myd; ssqr = myssqr; tau = mytau; bsqr = mybsqr; sdotb2 = mysdotb2
    

    ierror=0

    ! ierror=1 : error on entry: must have D>=0, tau>=smallp/(gamma-1)
    ! ierror=2 : maxitnr reached without convergence
    ! ierror=3 : no range for solution was find
    ! ierror=4 : final xi <smallxi
    ! ierror=5 : lower b/ound non-negative f(xi) value
    ! ierror=6 : v^2>1
    ! ierror=7 :  solution out from the range?
    ! ierror=8 : nonmonotonic function f?
    ! ierror=13: 2DW: maxitnr reached without convergence
    ! ierror>20: Error in Entropy inversion

    ! Check if guess is close enough:
    call funcd(xi,f,df,lfacl,d,ssqr,tau,bsqr,sdotb2,ierror) 
    if (ierror .eq. 0 .and. abs(f/df)<absaccnr) then
       xii   = xi - f/df
       lfaci = lfacl
       if (xii.gt.smallxi .and. lfaci.lt.lfacmax .and. lfaci.ge.1) then
        xi = xii; lfac = lfaci    
        return
       end if
    else
       ierror = 0
    end if


    if(d<minrho)then
       ierror=1
       return
    endif
    if(tau<smalltau)then
       ierror=1
       
       return
    endif



    ! Hydro case: handle exactly as before (using p to iterate on)
    
    if(bsqr<=limitvalue .and. typeinversion_ .eq. t1d_)then
       call con2primHydro(lfac,xi,d,ssqr,tau,ierror)
       return 
    end if ! Hydro
    

    ! starting point for NR on xi
    ! xi1 is meant to be the lower bound for the bracket to be used in NR
    dplus=d+minp*eqpar(gamma_)/(eqpar(gamma_)-one)
    xi1=dplus

    vsqr = ((ssqr + (xi1*two*sdotb2 + bsqr*sdotb2)/(xi1*xi1))&
       /((xi1+bsqr)*(xi1+bsqr)))


    !=============================================!
    !=== find new xi1, in the case v2(xi1) > maxvel^2 = (1-dmaxvel)^2 
    ! locate xi corresponding to maximal velocity allowed, namely 1-dmaxvel
    niiter=-1
    if(vsqr > (one-dmaxvel)**2 )then

       xiprev=xi1
       LoopVmax:  do ni = 1,maxitnr

          ! v^2(xi1) - maxvel^2
          fv= ((ssqr + (xi1*two*sdotb2 + bsqr*sdotb2)/(xi1*xi1))&
             /((xi1+bsqr)*(xi1+bsqr)))-(one-dmaxvel )**2
          ! d(v^2(xi)-maxvel^2)/dxi
          dfv= -two * (sdotb2*(3.0d0*xi1*(xi1+bsqr)+bsqr*bsqr)+ssqr*xi1**3)&
             / ((xi1*(xi1+bsqr))**3)

          if(fv*dfv==zero) then
             if(fv==zero)then
                exit LoopVmax
             else
                !print *,'stop: dfv becomes zero, non-monotonic function of xi!!'
                ierror=8
                
                return
             endif
          else
             xiprev=xi1
             xi1   =xi1 -fv/dfv
             if(fv*dfv>zero)then
                ! xi-iterate decreased
                ! restrict to left
                xi1=max(xi1,(dplus+xiprev)*half)
             else ! fv*dfv <0 
                ! xi-iterate increased
                ! restrict to right
                xi1=min(xi1,(tau+d+xiprev)*half)
             endif
          endif
          er=dabs(fv/dfv)/xi1
          if((er<tolernr).or.(dabs(fv/dfv)<absaccnr))exit LoopVmax
          niiter=ni
       enddo LoopVmax
    endif
    !=============================================!


    if(niiter==maxitnr)then
       ! could not find consistent value of lower bound for xi compliant with maxvel
       !     print *,' could not find value of lower bound for xi compliant with maxvel'
       !     print *,'xi1=',xi1,'dplus=',dplus,' tau+d=',tau+d
       !     print*, bsqr, sdotb, ssqr, tau, d, lfac, xi
       !     print*, xi,lfac
       ierror=2
       
       return
    endif

    ! we now compute f(xi1) and lfac(xi1)

    call funcd(xi1,fl,df,lfacl,d,ssqr,tau,bsqr,sdotb2,ierror)
    if(ierror /=0) then
       
       return
    end if

    if(fl>zero)then
       ierror=5
       
       return
    end if

    !--------------------------------------------------------!
    ! xi2 is meant to be the maximal bound for the bracket to be used in NR
    ! for this we take the value of E=tau+d, increased by smallp

    xi2= max(tau+d+minp - half *bsqr,10.0d0*xi1)
    niiter=-1

    LoopxiMax : do ni=1,maxitnr
       ! we now compute f(xi2) and lfac(xi2)
       call funcd(xi2,fh,df,lfach,d,ssqr,tau,bsqr,sdotb2,ierror)
       if(ierror /=0) then
          
          return
       end if
       testsh=(xi2>=dplus*lfach.and.lfach<=lfacmax)

       ! maximal bound found when f(xi1) opposite sign of f(xi2)
       !, enforce consistency tests on xi2
       if (testsh.and.fh *fl <=zero) exit LoopxiMax
       !==== Zak 17/05 fast convergence ====!
       xi1=xi2
       fl=fh
       lfacl=lfach
       testsl=testsh
       !====================================!
       xi2=two*xi2
       niiter=ni
    end do   LoopxiMax
    !--------------------------------------------------------!

    if(niiter == maxitnr .or. (fh*fl>zero))then
       ! could not find upper bound for NR on xi
       !     print *,'could not find upper bound for NR on xi'
       !     print *,'niiter=',niiter,' versus maxitnr=',maxitnr
       !     print *,'xi1=',xi1,' fl=',fl,' lfacl=',lfacl,' vs        dplus=',dplus
       !     print *,'xi2=',xi2,' fh=',fh,' lfach=',lfach,' vs tau+d+smallp=',tau+d+smallp
       ierror = 3
       
       return
    end if

    finished=.false.
    if(fl==zero)then
       xii=xi1
       lfaci=lfacl
       finished=(xi1>=dplus*lfacl.and.lfacl<=lfacmax)
    endif

    if(fh==zero)then
       xii=xi2
       lfaci=lfach
       finished=(xi2>=dplus*lfach.and.lfach<=lfacmax)
    endif

    if(finished)then
       xi=xii
       lfac=lfaci
       return
    end if


    xil=xi1
    xih=xi2
    ! opedit: Try to take previous values:
    if (xil .le. xi .and. xi .le. xih) then
       xii = xi
       !     lfaci = lfac
       !   print*, 'taking previous values for initial guess, lfac, xi:', lfaci, xii
    else
       xii=half*(xih+xil)    !Initialize the guess for rootfinder
       !   print*, 'previous values not good for initial guess, lfac, xi:', lfaci, xii
    end if


    ! Consistency check on initial xi value:
    !=============================================!
    !=== Increase xii in case v2>1:
    !=============================================!
    vsqr = ((ssqr + (xii*two*sdotb2 + bsqr*sdotb2)/(xii*xii))&
       /((xii+bsqr)*(xii+bsqr)))
    niiter=-1
    if(vsqr > (one-dmaxvel)**2 )then
       loopxi1: do ni =1, maxitnr

          xii = xii * 10.0d0
          vsqr = ((ssqr + (xii*two*sdotb2 + bsqr*sdotb2)/(xii*xii))&
             /((xii+bsqr)*(xii+bsqr)))
          if (vsqr  .lt. (one-dmaxvel)**2 ) exit

       end do loopxi1
    end if

    call funcd(xii,f,df,lfaci,d,ssqr,tau,bsqr,sdotb2,ierror)

    if(ierror /=0) then
       
       return
    end if

    ! Save the initial guess: 
    xi_guess = xii
    lfac_guess = lfaci

    !----------------------------------------------------
    !--- We got the bracketing and the initial guess ----
    !--- Now do the Newton-Raphson iterations -----------
    !----------------------------------------------------

    select case (typeinversion_)
    case(t1d_)
       call invert_1DW

    case(t2d_)
       call invert_2DW

    case(t1d2d_)
       call invert_1DW
       if (ierror /= 0 ) then
          ierror = 0
          call invert_2DW
       end if

    case(t2d1d_)
       call invert_2DW
       if (ierror /= 0 ) then
          ierror = 0
          call invert_1DW
       end if

    case(t2d1dentropy_)
        call mpistop("Set #ifdef ENTROPY flag when you want to use Entropy &
           inversion") 
       

    case(t1d2dentropy_)
        call mpistop("Set #ifdef ENTROPY flag when you want to use Entropy &
           inversion") 
       

    case(t1dentropy_)
        call mpistop("Set #ifdef ENTROPY flag when you want to use Entropy &
           inversion") 
       

    case(t2dentropy_)
        call mpistop("Set #ifdef ENTROPY flag when you want to use Entropy &
           inversion") 
       

    case(tentropy_)
        call mpistop("Set #ifdef ENTROPY flag when you want to use Entropy &
           inversion") 
       

    case default
       call mpistop('Unknown typeinversion in con2prim')
    end select

    !-----------------------------------------------------------------------------
  contains
    !=============================================================================
    
    !=============================================================================
    subroutine invert_2DW()

      use mod_rtsafe2D
      ! .. local ..
      double precision                  :: x, xl, xh!, xold,yold
      !-----------------------------------------------------------------------------

      xii = xi_guess; lfaci = lfac_guess


      maxit = maxitnr
      x  = 1.0d0 - 1.0d0/lfaci**2
      xl = 0.0d0
      xh = (1.0d0-dmaxvel)**2

      !    xold = x
      !    yold = xii

      call rtsafe2D(f1,f2,xl,xh,xil,xih,absaccnr,x,xii,niiter,ierror,validate&
         =validate_2DW)

      if (ierror .ne. 0) then
         ierror = ierror + 10
         !       print*, 'status:',ierror,niiter,nbisect,xold, xl, xh, yold,xil,xih
         !       print*, xii, x
         !       print*, bsqr, sdotb2, ssqr,tau,d
         !       print*, 'end'
         return
      end if

      !===============================!
      ! final values for auxiliary variables are now passed to w-array
      xi=xii
      lfac=sqrt(1.0d0/(1.0d0-x))
      !===============================!

    end subroutine invert_2DW
    !=============================================================================
    subroutine invert_1DW()

      xii = xi_guess

      er1 = one
      nit = 0
      niiter=-1
      xiprev=xii
      lastxiok=-one
      lastlfacok=-one
      lastf=-one

      !--- Start iteration ---!
      LoopNRRMHD :  do ni=1,maxitnr 
         nit = nit + 1
         if(nit>maxitnr/2)then
            ! mix the last  value for convergence
            xii=half*(xii+xiprev)
            ! relax accuracy requirement
            er1=10.0d0*er1
            ! following avoids decrease of accuracy requirement 
            ! *every* iteration step beyond maxitnr/2
            nit = nit - maxitnr/10
         endif
         call funcd(xii,f,df,lfaci,d,ssqr,tau,bsqr,sdotb2,ierror) 
         if(ierror /=0)return
         tests=(xii>=dplus*lfaci.and.lfaci<=lfacmax)
         if(tests)then
            lastxiok=xii
            lastlfacok=lfaci
            lastf=f
         endif
         if(f*df==zero) then
            if(f==zero) then
               tests=(xii>=dplus*lfaci.and.lfaci<=lfacmax)
               if(tests)then
                  exit ! zero found and consistency checks fullfilled
               else
                  ierror=7
                  return
               endif
            else
               ierror=8
               return
            endif
         else
            xiprev=xii
            xii   =xii -f/df
            if(f*df>zero)then
               ! xi-iterate decreased
               ! restrict to left
               xii=max(xii,xil)
            else
               ! xi-iterate increased
               ! restrict to right
               xii=min(xii,xih)
            endif
         endif
         er=dabs(f/df)/xii
         if((er<tolernr*er1).or.(dabs(f/df)<absaccnr))then
            call funcd(xii,f,df,lfaci,d,ssqr,tau,bsqr,sdotb2,ierror) 
            tests=(xii>=dplus*lfaci.and.lfaci<=lfacmax)
            if(tests)then
               exit LoopNRRMHD ! converged solution with tests ensured
            else
               ierror=7
               return

            endif
         endif
         niiter=ni

      enddo LoopNRRMHD

      if(niiter==maxitnr) then
         ! no convergence of NR for xi, although zero bracketed
         ! opedit: for the time being as prolongprimitive calls con2prim for
         ! uninitialized fields, commenting out the warning.
         ! fix in amr_gc_comm.t
!                print *,'no convergence of NR for xi, although zero bracketed'
!                print *,'er=',er,'tolernr=',tolernr,'er1=',er1,'df=',df,'absaccnr=',absaccnr
!                print *,'xii=',xii,' f=',f
!                print *,'brackets xil=',xil,' and xih=',xih,' with fl fh=',fl,fh
!                print *,'lastxiok=',lastxiok,'lastlfacok=',lastlfacok,'lastf=',lastf
         ierror=2
         return
      endif ! niiter==maxitnr

      !===============================!
      ! final values for auxiliary variables are now passed to w-array
      xi=xii
      lfac=lfaci
      !===============================!

    end subroutine invert_1DW
    !=============================================================================
  end subroutine con2prim
  !=============================================================================
  subroutine validate_2DW(xl,xh,yl,yh,xroot,yroot,xroot_old,yroot_old)
    double precision, intent(in)        :: xl, xh, yl, yh, xroot_old,&
        yroot_old
    double precision, intent(inout)     :: xroot, yroot
    !-----------------------------------------------------------------------------

    !       xroot = abs(xroot)
    !       if (xroot .gt. xh) xroot = xroot_old
    !       if (yroot .lt. yl) yroot = yl
    !       if (yroot .gt. yh) yroot = yh

    if (xroot .gt. xh) xroot = 0.5d0 * (xroot_old + xh)
    if (xroot .lt. xl) xroot = 0.5d0 * (xroot_old + xl)

    if (yroot .gt. yh) yroot = 0.5d0 * (yroot_old + yh)
    if (yroot .lt. yl) yroot = 0.5d0 * (yroot_old + yl)


  end subroutine validate_2DW
  !=============================================================================
  subroutine f1(x,y,f,dfdx,dfdy,return_derivatives)

    use mod_con2prim_vars

    double precision, intent(in)        :: x, y
    double precision, intent(out)       :: f, dfdx, dfdy
    logical, intent(in)                 :: return_derivatives
    !-----------------------------------------------------------------------------

    f    = (y + bsqr)**2 * x - sdotb2 * (2.0d0 * y + bsqr) / y**2 - ssqr

    if (return_derivatives) then
       dfdx = (bsqr + y)**2
       dfdy = - 2.0d0 * sdotb2 / y**2 + 2.0d0 * x * (bsqr+y) + 2.0d0 * sdotb2 &
          * (bsqr + 2.0d0 * y) / y**3
    end if

  end subroutine f1
  !=============================================================================
  subroutine f2(x,y,f,dfdx,dfdy,return_derivatives)

    use mod_con2prim_vars

    double precision, intent(in)        :: x, y
    double precision, intent(out)       :: f, dfdx, dfdy
    logical, intent(in)                 :: return_derivatives
    ! .. local ..
    double precision                    :: mylfac, dlfac
    double precision                    :: p, dpdxi, dpdx
    !-----------------------------------------------------------------------------

    mylfac = sqrt(1.0d0/(1.0d0-x))
    !===== Pressure, calculate using EOS =====!
    call FuncPressure(y,mylfac,d,ssqr,tau,0.0d0,p,dpdxi,dpdx)
    !=========================================!

    f    = y - tau - d - p + 0.5d0 * (1.0d0 + x) * bsqr - 0.5d0 * sdotb2 &
       / y**2

    if (return_derivatives) then
       dfdx = 0.5d0 * bsqr - dpdx
       dfdy = 1.0d0 + sdotb2/y**3  - dpdxi
    end if

  end subroutine f2
  !=============================================================================
  subroutine funcd(xi,F,dF,mylfac,d,ssqr,tau,bsqr,sdotb2,ierror)

    use mod_amrvacdef

    double precision, intent(in)  :: xi,d,ssqr,tau,bsqr,sdotb2
    double precision, intent(out) :: F,dF,mylfac
    integer, intent(inout)        :: ierror

    double precision  :: dlfac
    double precision  :: vsqr,p,dpdxi,dpdx
    !-----------------------------------------------------------------------------

    vsqr = ((ssqr + (xi*two*sdotb2 + bsqr*sdotb2)/(xi*xi))/((xi+bsqr)*&
       (xi+bsqr)))

    if (vsqr<one) then
       mylfac = one/dsqrt(one-vsqr)
       dlfac = -mylfac**3*(sdotb2*(3.0d0*xi*(xi+bsqr)+bsqr*bsqr)+ssqr*xi**3) &
          /((xi*(xi+bsqr))**3)

       !===== Pressure, calculate using EOS =====!
       call FuncPressure(xi,mylfac,d,ssqr,tau,dlfac,p,dpdxi,dpdx)
       !=========================================!
       F  = xi-tau-d+half*bsqr*(one+vsqr)-half*sdotb2/xi**2-p
       dF = one + sdotb2/xi**3  -dpdxi+dlfac*bsqr/mylfac**3

    else 
       ! print *,'Warning: erroneous input to funcd since vsrq=',vsqr,' >=1'
       ! print *,'input values d, ssqr, tau, bsqr, sdotb2:',d,ssqr,tau,bsqr,sdotb2
       !       print*,'ierror ==6 ',it,mype,t,saveigrid 
       ierror =6
       return
    end if


  end subroutine funcd
  !=============================================================================
  
 !=============================================================================
  subroutine con2primHydro(lfac,xi,d,sqrs,tau,ierror)
    !use ieee_arithmetic
    use mod_amrvacdef

 !this is a copy of the HD iteration, where we solve for p via NR, and modified
    ! to give xi on output

    double precision :: xi,lfac
    double precision :: d,sqrs,tau
    integer          :: ierror

    integer          :: ni,niiter
    double precision :: pcurrent,pnew,pL,pR
    double precision :: er,er1,ff,df,dp,v2
    double precision :: pmin,lfac2inv,pLabs,pRabs,pprev
    double precision :: s2overcubeG2rh
    double precision :: xicurrent
    double precision :: oldff1,oldff2
    double precision :: Nff
    double precision :: pleft,pright,pnewi
    integer          ::nit,n2it,ni2,ni3
    double precision :: h,dhdp
 !-----------------------------------------------------------------------------

    ierror=0
    ! ierror=0 : ok
    !
    ! ierror<>0
    !
    ! ierror=1 : error on entry: must have D>=minrho, tau>=smalltau
    ! ierror=2 : maxitnr reached without convergence
    ! ierror=3 : final pressure value < smallp or xi<smallxi during iteration
    ! ierror=4 : final v2=1,v2=1,v2=1 hence problem as lfac=1/0
    ! ierror=5 : nonmonotonic function f?
    ! ierror=7 : stop due to strictnr violation

    if(d<minrho .or. tau<smalltau) then
       ierror=1
       return
    endif

 !incase input pressure is not available or random value: replace by smallp


    ! left and right brackets for p-range
    pmin=dsqrt(sqrs)/(one-dmaxvel)-tau-d
    pLabs=max(minp,pmin)
    pRabs=1.0d99
    ! start value from input
    pcurrent=pLabs

    er1=one
    pprev=pcurrent

    ! Fudge Parameters
    oldff1=1.0d7  ! High number
    oldff2=1.0d9  ! High number bigger then oldff1
    n2it = 0
    nit  = 0



    LoopNR:  do ni=1,maxitnr
       nit = nit + 1
       !============= Controle 1=============!
       if(nit>maxitnr/4)then
          !print *,'ni,er,p',ni,er,pcurrent
          ! mix pressure value for convergence
          pcurrent=half*(pcurrent+pprev)
          ! relax accuracy requirement
          er1=10.*er1
          nit = nit - maxitnr/10
       endif
       !=======================================!

       niiter=ni  
       xicurrent=tau+d+pcurrent

       if(xicurrent<smallxi) then
          !        print*,'!--- amrvacphys/t.srhd-- con2prim ---!'
          !        print *,'stop: too small xi iterate:',xicurrent
          !        print *,'for pressure iterate p',pcurrent
          !        print *,'pressure bracket pLabs pRabs',pLabs,pRabs
          !        print *,'iteration number:',ni
          !        print *,'values for d,s,tau,s2:',d,sqrs,tau,sqrs
          ierror=3 
          return
       endif

       v2=sqrs/xicurrent**2
       lfac2inv=one - v2
       if(lfac2inv>zero) then
          lfac=one/dsqrt(lfac2inv)
       else
          !        print*,'!--- amrvacphys/t.srhd-- con2prim ---!'
          !        print *,'stop: negative or zero factor 1-v2:',lfac2inv
          !        print *,'for pressure iterate p',pcurrent
          !        print *,'pressure bracket pL pR',pL,pR
          !        print *,'absolute pressure bracket pLabs pRabs',pLabs,pRabs
          !        print *,'iteration number:',ni
          !        print *,'values for d,s,tau,s2:',d,sqrs,tau,sqrs
          !        print *,'values for v,xi:',v2,xicurrent
          ierror=4
          return
       endif

       s2overcubeG2rh=sqrs/(xicurrent**3)
       !== ZM calculation done using the EOS ==!
       call FuncEnthalpy(pcurrent,lfac2inv,d,tau,sqrs,xicurrent,&
            s2overcubeG2rh,h,dhdp,ierror)
       !=======================================!   
       ff=-xicurrent*lfac2inv + h 
       df=- two*sqrs/xicurrent**2  + dhdp - lfac2inv

       if (ff*df==zero) then
          if (ff==zero) then
             exit ! zero found
          else
             !           print *,'stop: df becomes zero, non-monotonic f(p)!!'
             ierror=5
             return
          endif
       else 
          pnew=pcurrent-ff/df
          if (ff*df>zero) then
             ! pressure iterate has decreased
             ! restrict to left 
             pnew=max(pnew,pLabs)
          else  ! ff*df<0
             ! pressure iterate has increased
             ! restrict to right 
             pnew=min(pnew,pRabs)
          endif
       endif

       !===============================================!
       dp=pcurrent-pnew
       er=two*dabs(dp)/(pnew+pcurrent)
       if(((er<tolernr*er1).or.(dabs(dp)<absaccnr))) exit LoopNR
       !===============================================!

       ! For very small values of pressure, NR algorithm is not efficient to
       ! find root, use Euler algorithm to find precise value of pressure
       if((dabs(oldff2-ff) < 1.0d-8 .or. niiter >= maxitnr-maxitnr/20).and.&
            ff * oldff1 < zero    .and.  dabs(ff)>absaccnr)then

          n2it=n2it+1
          if(n2it<=3) pcurrent=half*(pnew+pcurrent)
          if(n2it>3)then
             pright =pcurrent
             pleft=pprev
             pcurrent=half*(pleft+pright)
             Dicho:  do ni3=1,maxitnr
                !===================!
                xicurrent=tau+d+pcurrent
                v2=sqrs/xicurrent**2
                lfac2inv=one - v2

                if(lfac2inv>zero)then
                   lfac=one/dsqrt(lfac2inv)
                else
                   ierror=4
                   return
                endif
                !===================!


                !== ZM calculation done using the EOS ==!
                call Bisection_Enthalpy(pnew,lfac2inv,d,&
                     tau,sqrs,xicurrent,h,ierror)
                Nff=-xicurrent*lfac2inv + h 
                !=======================================!
                !==== Iterate ====!
                if(ff * Nff < zero)then
                   pleft=pcurrent
                else
                   pright=pcurrent
                endif

                pcurrent=half*(pleft+pright)
                !==================!

                !=== The iteration converge ===!
                if(2.0d0*dabs(pleft-pright)/(pleft+pright)< absaccnr &
                     .or. dabs(ff)<absaccnr)then
                   pnew=pcurrent
                   exit LoopNR
                endif
                !==============================!

                !=== conserve the last value of Nff ===!
                ff=Nff
                !======================================!
             enddo    Dicho
          endif

       else
          !====== There is no problems, continue the NR iteration ======!
          pprev=pcurrent
          pcurrent=pnew
          !=============================================================!
       endif


       !=== keep the values of the 2 last ff ===!
       oldff2=oldff1
       oldff1=ff
       !========================================!
    enddo LoopNR

    if(niiter==maxitnr)then
       ierror=2
       return
    endif

    if(pcurrent<minp) then
       ierror=3
       return
    endif

    !------------------------------!
    xi=tau+d+pcurrent
    v2=sqrs/xicurrent**2
    lfac2inv=one - v2
    if(lfac2inv>zero) then
       lfac=one/dsqrt(lfac2inv)
    else
       ierror=4
       return
    endif

  end subroutine con2primHydro
  
  !=============================================================================
end module mod_con2prim
!=============================================================================
