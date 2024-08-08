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

!FILE GETPISEN.FOR
!
! 24/11/2013 adopted from Serguei Komisssarov's isen.f
! Inversion for isentropic flow. 
! This is used when -eos=iso in the relativistic mhd module
!
! ----------------------------------------------------------------
        SUBROUTINE GETPISEN(MM,RHO,MB2,BB,KAPPA_I,U0,XI,IFOK)
! ----------------------------------------------------------------
!  Description:
!     Iterative solver which is used to find density pressure
!     and total 4-velocity from energy, rest-mass, total
!     momentum density , and magnetic field B.
!    ( isentropic equation of state ) 

use mod_amrvacdef

! ----------------------------------------------------------------
!  Input:
        double precision, intent(in)             :: MM,RHO,MB2,BB,KAPPA_I
!         MM - momentum density m[i]m[i]; i=1,3
!         RHO - mass density
!         MB2  - scalar product m[i]B[i] m[j]B[j]
!         BB  - scalar product B[i]B[i]
!         KAPPA_I - entropy p/rho**gamma        
! ----------------------------------------------------------------
!   Input/Output:
        double precision, intent(inout)          :: U0,XI
        integer, intent(out)                     :: IFOK
!         U0   - Lorentz factor
!         XI   - Auxilary variable xi = rho h lfac**2
!         It is assumed that when the subroutine is called U0 gives the 
!         initial guess (solution from the previous step). 
! ----------------------------------------------------------------
!   Commons:
!        INCLUDE 'FLOOR.CMN'
!        INCLUDE 'STATE.CMN'
! ----------------------------------------------------------------
!   Constants:
        double precision ALFA
        PARAMETER (ALFA=1.E-10) ! precision parameter
        INTEGER CLIM
        PARAMETER (CLIM=30)  ! maximum number of iterations
       
! ----------------------------------------------------------------
!   Variables:
        INTEGER I
        double precision LOR1,LOR2,F,DFDX
        double precision GAM, D
        double precision U0INIT, XIINIT
! ----------------------------------------------------------------

        U0INIT=U0
        XIINIT=XI
        GAM = eqpar(gamma_)
                 
        IFOK = 0 ! everyting is fine

        LOR1=U0
        DO I=1,CLIM+1    ! Newton iterations
             IF(I.EQ.CLIM)THEN 
                 IFOK= 1
                 RETURN
             ENDIF
             CALL FDFDX(MM,RHO,MB2,BB,LOR1,F,DFDX)
             LOR2=LOR1-F/DFDX
!             print *,i,lor1,F,DFDX
             IF(ABS(LOR2-LOR1).LT.ALFA) GOTO 10
             LOR1=LOR2
         ENDDO

10       U0=LOR2
         IF(U0.LT.1.)IFOK=2
         D=RHO/U0
         XI=RHO*U0 + (GAM/(GAM-1.0d0)) * KAPPA_I*D**GAM * U0*U0

contains
! -------------------------------------------------------------
         SUBROUTINE FDFDX(MM,RHO,MB2,BB,LOR,FX,DFDX)
! -------------------------------------------------------------
implicit none
!   Input:
      double precision MM,RHO,MB2,BB,LOR
!   Output:
      double precision FX,DFDX
!   Commons:
!      INCLUDE 'STATE.CMN'
!   Variables: 
      double precision A,DF,DA,V2,F,F2,F3,F4,P,E,B4

      B4=BB**2
      V2=(LOR**2-1.D0)/LOR**2
      P=KAPPA_I*RHO**GAM
      E=(GAM/(GAM-1.D0))*P
      F=RHO*LOR + E*LOR**(2.D0-GAM)
      F2=F*F
      F3=F2*F
      F4=F2*F2
      A=V2*(F4+B4*F2+2.D0*BB*F3)
      FX=BB*MB2 + MM*F2+2.D0*MB2*F - A
      DF=RHO+E*(2.D0-GAM)*LOR**(1.D0-GAM)
      DA = (F4+B4*F2+2.D0*BB*F3)*2.D0/LOR**3
      DA = DA+(4.D0*F3+2.D0*B4*F+6.D0*BB*F2)*V2*DF 
      DFDX = (2.D0*MM*F+2.D0*MB2)*DF-DA

end subroutine FDFDX

end subroutine GETPISEN

! STATE.CMN
!        INTEGER IEOS
!        COMMON /EOS_TYPE/ IEOS
!        REAL*8 GAM, KAPPA_I
!        COMMON /ISEN_POLY/ GAM,KAPPA_I
