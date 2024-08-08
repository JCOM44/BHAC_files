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
module mod_amrvacdef

  !use mpi ! Well, this seemed like a good idea but does not compile on iboga.
  use mod_indices
  use mod_physicaldata
  use mod_connectivity
  use mod_amrvacpar

  IMPLICIT NONE

  ! DEFINITIONS OF GLOBAL PARAMETERS AND VARIABLES
  ! Parameters:

  ! Indices for cylindrical coordinates FOR TESTS, negative value when not used:
  INTEGER,PARAMETER:: r_=1, phi_=2, z_=3

  ! Indices for cylindrical coordinates FOR INDEXING, always positive
  INTEGER,PARAMETER:: pphi_=2, zz_=3

  !include 'amrvacpar.f90'

  ! xprob: problem box; iprob: problem
  INTEGER:: iprob
  DOUBLE PRECISION:: xprobmin1,xprobmin2,xprobmin3,xprobmax1,xprobmax2,&
     xprobmax3

  INTEGER,PARAMETER:: ndim=3, ndir=3

  include 'amrvacsettings.f90'

  INTEGER,PARAMETER:: filelog_=1,fileout_=2,fileslice_=3,filecollapse_&
     =4,fileanalysis_=5,fileshell_=6,nfile=6 !outputfiles
  INTEGER,PARAMETER:: nslicemax=1000, nshellmax=1000

  INTEGER,PARAMETER:: unitstdin=5,unitterm=6,uniterr=6 ! Unit names.

  ! Units reserved for files:
  INTEGER,PARAMETER:: unitpar=9
  INTEGER,PARAMETER:: unitconvert=10
  INTEGER,PARAMETER:: unitslice=11
  INTEGER,PARAMETER:: unitsnapshot=12
  INTEGER,PARAMETER:: unitcollapse=13
  INTEGER,PARAMETER:: unitanalysis=14
  INTEGER,PARAMETER:: unitshell=15

  INTEGER,PARAMETER:: biginteger=10000000

  ! Note: smalldouble must be above machine precision 
  DOUBLE PRECISION,PARAMETER:: smalldouble=1.D-12, bigdouble=1.D+99
  DOUBLE PRECISION,PARAMETER:: zero=0D0,one=1D0,two=2D0,half=0.5D0,quarter&
     =0.25D0,third=0.33333333333333333333d0
  DOUBLE PRECISION,PARAMETER:: dpi&
     =3.141592653589793238462643383279502884197169399375105d0
  DOUBLE PRECISION,PARAMETER:: de=2.71828182845904523536028747135

  ! Physical scaling parameters:
  DOUBLE PRECISION:: UNIT_LENGTH, UNIT_DENSITY, UNIT_VELOCITY

  include 'amrvacusrpar.f90'

  ! For transform variables and save selected data
  INTEGER :: nwtf
  INTEGER :: neqpartf

  !Kronecker delta and Levi-Civita tensors
  INTEGER:: kr(0:3,0:3),lvc(1:3,1:3,1:3)

  !Equation and method parameters
  DOUBLE PRECISION:: eqpar(neqpar+nspecialpar)

  ! Time step control parameters
  DOUBLE PRECISION :: courantpar, dtpar, dtdiffpar, dtTCpar
  CHARACTER*131 :: typecourant,typeresid
  LOGICAL :: time_accurate, addmpibarrier

  !Time parameters
  INTEGER,PARAMETER:: nsavehi=100       ! maximum No. saves into outputfiles
  ! defined by arrays of tsave or itsave
  DOUBLE PRECISION:: t,tmax,dtmin,residmin,residmax,residual
  DOUBLE PRECISION:: tfixgrid
  DOUBLE PRECISION:: tsave(nsavehi,nfile),tsavelast(nfile),dtsave(nfile),&
     slicecoord(nslicemax),shellcoord(nshellmax)
  LOGICAL:: tmaxexact,treset,itreset,firstprocess,resetgrid,fixprocess,&
     changeglobals,collapse(ndim)
  INTEGER:: it,itmax,itmin,slowsteps
  INTEGER:: itsave(nsavehi,nfile),itsavelast(nfile),ditsave(nfile)
  INTEGER:: isavet(nfile),isaveit(nfile), nslices, nshells,&
      slicedir(nslicemax), collapseLevel, nxShell2,nxShell3
  INTEGER:: n_saves(1:nfile)
  INTEGER:: typeparIO
  INTEGER:: itfixgrid,ditregrid
  INTEGER:: nwauxio
  INTEGER:: istep, nstep

  !Method switches
  CHARACTER*131 :: typeadvance
  CHARACTER*131 :: typelow1(nlevelshi),typelimited,typesourcesplit
  CHARACTER*131 :: typefull1(nlevelshi), typepred1(nlevelshi)
  CHARACTER*131 :: typelimiter1(nlevelshi),typegradlimiter1(nlevelshi)
  INTEGER       :: type_limiter(nlevelshi), type_gradient_limiter(nlevelshi)
  CHARACTER*131 :: typeentropy(nw),typetvd,typeaverage
  CHARACTER*131 :: typedimsplit,typeaxial,typecoord,typepoly,typeinversion,&
     typeinversionresis
  INTEGER:: errorestimate,nxdiffusehllc,typespherical,ncyclemax
  CHARACTER*131 :: typeprolonglimit, clean_init_divB
  INTEGER       :: typelimiter, typegradlimiter
  DOUBLE PRECISION:: entropycoef(nw)
  DOUBLE PRECISION:: tvdlfeps, mcbeta, parastsnu, TCphi, betamin
  LOGICAL:: sourceparasts,sourceimpl,sourceimplcycle,conduction,TCsaturate,&
     energyonly
  LOGICAL:: loglimit(nw),logflag(nw),flathllc,flatcd,flatsh
  LOGICAL:: ssplitdivb,ssplitresis,ssplituser,useprimitive,dimsplit,&
     use_particles
  LOGICAL:: restrictprimitive,prolongprimitive, coarsenprimitive,&
     useprimitiveRel, amrentropy

  LOGICAL:: BnormLF
  DOUBLE PRECISION:: smallT,smallp,smallrho,amr_wavefilter(nlevelshi)
  CHARACTER*131 :: typediv,typegrad,typeemf

  LOGICAL:: strictnr,strictgetaux
  DOUBLE PRECISION::dmaxvel,tolernr,absaccnr,tlow
  INTEGER:: maxitnr,nflatgetaux

  LOGICAL:: convert_nocartesian, slice_nocartesian
  LOGICAL:: writew(nw),writelevel(nlevelshi)
  DOUBLE PRECISION:: writespshift(ndim,2)
  INTEGER:: level_io, level_io_min, level_io_max

  ! local and global fastest wave speed (computed in setdt):
  DOUBLE PRECISION :: cmax_mype, cmax_global

  !Boundary region parameters
  INTEGER,PARAMETER:: nhiB=2*ndim         ! maximum No. boundary sections
  LOGICAL:: periodB(ndim), poleB(2,ndim), aperiodB(ndim), primitiveB(2,ndim)
  CHARACTER*131 :: typeB(nw,nhiB)
  CHARACTER*131 :: typeghostfill,typegridfill
  DOUBLE PRECISION::ratebdflux
  LOGICAL:: internalboundary

  !File parameters
  CHARACTER(len=131) :: inifile,filenameout,filenameini,filenamelog
  CHARACTER*131 :: fileheadout
  CHARACTER*1024 :: wnames,primnames,wnameslog
  CHARACTER*131 :: typefilelog
  INTEGER :: snapshotini
  LOGICAL :: sliceascii


  !Convert parameters
  LOGICAL :: convert,autoconvert,saveprim,uselimiter,endian_swap
  CHARACTER*131 :: convert_type, dxfiletype, collapse_type, slice_type,&
      shell_type
  DOUBLE PRECISION :: normvar(0:nw),normt
  ! --------------------------------------------
  INTEGER:: saveigrid
  ! Stores the memory and load imbalance, to be used in printlog:
  DOUBLE PRECISION :: Xload, Xmemory

  DOUBLE PRECISION:: time_bc

  integer,parameter:: nodehi=3+1
  integer,parameter:: plevel_=1
  integer,parameter:: pig1_=plevel_+1,pig2_=plevel_+2,pig3_=plevel_+3

  integer,parameter:: rnodehi=3*3
  integer,parameter:: rpxmin0_=0
  integer,parameter:: rpxmin1_=rpxmin0_+1,rpxmin2_=rpxmin0_+2,rpxmin3_&
     =rpxmin0_+3 
  integer,parameter:: rpxmax0_=3
  integer,parameter:: rpxmax1_=rpxmax0_+1,rpxmax2_=rpxmax0_+2,rpxmax3_&
     =rpxmax0_+3 
  integer,parameter:: rpdx1_=2*3+1,rpdx2_=2*3+2,rpdx3_=2*3+3

  ! parameters for bc_phys
  integer,parameter:: ismin1=-1+2*1,ismin2=-1+2*2,ismin3=-1+2*3
  integer,parameter:: ismax1=2*1,ismax2=2*2,ismax3=2*3

  include 'mpif.h'

  ! nodal info per grid: 
  ! rnode:            corner coordinates,dx(idim)
  ! node:             dimensions and pointers 

  ! mxnest:           maximal number of levels

  ! tol:              error tolerance used for refinement decision
  ! nbufferx^D:       number of cells as buffer zone

  ! dixB:             number of ghost cells surrounding a grid for boundary cndt.

  double precision :: tol(nlevelshi), tolratio(nlevelshi)
  double precision :: rnode(rnodehi,ngridshi), rnode_sub(rnodehi,ngridshi),&
      dx(ndim,nlevelshi), dt, dtimpl, dt_grid(ngridshi), dxlevel(ndim)
  integer          :: mxnest, dixB, nbufferx1,nbufferx2,nbufferx3
  integer          :: node(nodehi,ngridshi), node_sub(nodehi,ngridshi),&
      levmin, levmax, levmax_sub
  logical          :: skipfinestep, time_advance

  !$OMP THREADPRIVATE(dxlevel, saveigrid, typelimiter, typegradlimiter)
end module mod_amrvacdef
!============================================================================
