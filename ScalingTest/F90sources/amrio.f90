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
subroutine readcommandline

!use M_kracken
use mod_amrvacdef
integer           :: len, ier
logical           :: help

INTEGER :: i,stat
CHARACTER(len=32) :: arg
LOGICAL :: unknown_arg

!----------------------------------------------------------------------------

! =================== Kracken command line reading =================
!
!!defaults and usage:
!call kracken('cmd','-i amrvac.par -if data -restart -1 -slice 0 -collapse 0 '//&
!  '-shell 0 --help .false. -convert .false.')

!! Getting the filename
!      call retrev('cmd_i',inifile,len,ier)
!      call retrev('cmd_if',filenameini,len,ier)
!      snapshotini = iget('cmd_restart')
!      snapshotnext= snapshotini+1
!      slicenext   = iget('cmd_slice')
!      collapseNext   = iget('cmd_collapse')
!      shellNext      = iget('cmd_shell')
!      help = lget('cmd_-help')                    ! --help present?
!      convert = lget('cmd_convert')               ! -convert present?
!
! ==================================================================

! =============== Fortran 2003 command line reading ================

  ! Default command line arguments

  inifile="amrvac.par"

  filenameini="data"

  snapshotini=-1

  snapshotnext=0

  slicenext=0

  collapsenext=0

  shellnext=0

  help=.false.

  convert=.false.

  ! Argument 0 is program name, so we start from one
  i = 1

  unknown_arg=.false.

  DO
    CALL get_command_argument(i, arg)

    IF (LEN_TRIM(arg) == 0) EXIT

    select case(arg)
    case("-i")

        i = i+1
    
        CALL get_command_argument(i, arg)
    
        inifile=TRIM(arg)

    case("-if")

        i = i+1
    
        CALL get_command_argument(i, arg)
    
        filenameini=TRIM(arg)    

    case("-restart")

        i = i+1
    
        CALL get_command_argument(i, arg)
   
        read(arg,*,iostat=stat) snapshotini

        snapshotnext = snapshotini+1
 
    case("-slice")

        i = i+1
    
        CALL get_command_argument(i, arg)
   
        read(arg,*,iostat=stat) slicenext


    case("-collapse")

        i = i+1
    
        CALL get_command_argument(i, arg)
   
        read(arg,*,iostat=stat) collapsenext

    case("-shell")

        i = i+1
    
        CALL get_command_argument(i, arg)
   
        read(arg,*,iostat=stat) shellnext

    case("-convert")

        convert=.true.

    case("--help")

        help=.true.

        EXIT

    case default

        unknown_arg=.true.

        help=.true.

        EXIT

    end select

    i = i+1

  END DO

  if (unknown_arg) then
       print*,"======================================="
       print*,"Error: Command line argument ' ",TRIM(arg)," ' not recognized"
       print*,"======================================="

        help=.true.

  end if

! ==================================================================


if (mype==0) then
   print*,'-----------------------------------------------------------------------------'
   print*,'-----------------------------------------------------------------------------'
   print*,'                         ____  _    _          _____                         '
   print*,'                        |  _ \| |  | |   /\   / ____|                        '
   print*,'                        | |_) | |__| |  /  \ | |                             '
   print*,'                        |  _ <|  __  | / /\ \| |                             '
   print*,'                        | |_) | |  | |/ ____ \ |____                         '
   print*,'                        |____/|_|  |_/_/    \_\_____|                        '
   print*,'                        -----------------------------                        '
   print*,'                        The Black Hole Accretion Code                        '
   print*,'-----------------------------------------------------------------------------'
   print*,'-----------------------------------------------------------------------------'
end if

if (help) then
   if (mype==0) then 
      print*,'calling example:                                         '
      print*,'./bhac -i parameterfile -restart 100 -convert -if datamr/data'
      print*,'default parameterfile is amrvac.par                            '
      print*,'Note that parameterfile parameters overwrite the commandline   '
      print*,'-----------------------------------------------------------------------------'
      print*,'Available options are:'
      print*,'-i parfile'
      print*,'-if filenameout'
      print*,'-restart #snapshot'
      print*,'-convert'
      print*,'-slice #slicenext'
      print*,'-collapse #collapsenext'
      print*,'-shell #shellnext'
      print*,'--help 	Prints the help message'
      print*,'-----------------------------------------------------------------------------' 
      print*,'See COMMAND LINE PARAMETERS entry in the manual for more information'
      print*,'(https://bhac.science/documentation_html/commandline.html).'
      print*,'                                                         '
   endif
   call comm_finalize
   STOP
endif

end subroutine readcommandline
!=============================================================================
subroutine readparameters

use mod_limiter
use mod_amrvacdef


logical :: fileopen
integer :: ifile, iB, isave, iw, level, idim, islice
integer :: nxlone1,nxlone2,nxlone3
double precision :: dxlone1,dxlone2,dxlone3

namelist /filelist/   filenameout,filenameini,filenamelog, snapshotini,&
   typefilelog,firstprocess,resetgrid,changeglobals,snapshotnext, convert,&
   convert_type,slice_type,dxfiletype,saveprim,primnames,typeparIO,uselimiter,&
   nwauxio,convert_nocartesian,slice_nocartesian,addmpibarrier, writew,&
   writelevel,writespshift,endian_swap, normvar,normt,level_io,level_io_min,&
   level_io_max, autoconvert,slicenext,collapseNext,collapse_type, shellNext,&
   shell_type
namelist /savelist/   tsave,itsave,dtsave,ditsave,nslices,slicedir,slicecoord,&
   collapse,collapseLevel,nshells,shellcoord,nxShell2,nxShell3
namelist /stoplist/   itmax,tmax,tmaxexact,dtmin,t,it,treset,itreset,residmin,&
   residmax,typeresid
namelist /methodlist/ wnames,fileheadout,typeadvance, ssplitdivb,ssplitresis,&
   ssplituser,typesourcesplit,sourceimpl,sourceimplcycle,conduction,&
   TCsaturate,TCphi,ncyclemax,sourceparasts,parastsnu,dimsplit,typedimsplit,&
   typeaxial,typecoord,typefull1,typepred1,typelow1,typelimiter1,mcbeta,&
   typegradlimiter1,flatcd,flatsh,typelimited,useprimitive, typetvd,&
   typeentropy,entropycoef,typeaverage, B0field,Bdip,Bquad,Boct,Busr,&
   useprimitiveRel,maxitnr,dmaxvel,typepoly,tvdlfeps,BnormLF,smallT,smallp,&
   smallrho,typegrad,typediv,tolernr,absaccnr,strictnr,strictgetaux,&
   nflatgetaux,nxdiffusehllc,typespherical,fixprocess,flathllc, tlow,nwtf,&
   neqpartf,typeinversion,typeemf,typeinversionresis, clean_init_divB,&
    use_particles, betamin
namelist /boundlist/  dixB,typeB,typeghostfill,typegridfill,ratebdflux,&
   internalboundary,primitiveB
namelist /amrlist/    mxnest,nbufferx1,nbufferx2,nbufferx3,tol,tolratio,&
   errorestimate, amr_wavefilter,nxlone1,nxlone2,nxlone3,dxlone1,dxlone2,&
   dxlone3,iprob,xprobmin1,xprobmin2,xprobmin3,xprobmax1,xprobmax2,xprobmax3,&
    skipfinestep,wflags,flags,restrictprimitive,prolongprimitive,&
   coarsenprimitive, typeprolonglimit, amrentropy,logflag,tfixgrid,itfixgrid,&
   ditregrid
namelist /paramlist/  time_accurate, courantpar, dtpar, dtdiffpar, dtTCpar,&
   typecourant, slowsteps
!----------------------------------------------------------------------------

! defaults for boundary treatments
ratebdflux=one
typeghostfill='linear' 
dixB=2
typeB(1:nw,1:nhiB)='cont'
internalboundary=.false.
primitiveB(1:2,1:ndim) = .false.

! code behavior and testing defaults
addmpibarrier=.false.

! defaults for specific options
fixprocess=.false.
typegrad='central'
typediv='central'
smallT=-one
smallp=-one
smallrho=-one

! relativistic module defaults
useprimitiveRel=.true.
strictnr=.true.
strictgetaux=.false.
nflatgetaux=1
typepoly='gammie'
tolernr=1.0d-13
absaccnr=1.0d-13
maxitnr=100
dmaxvel=1.0d-8
tlow=zero

! defaults for convert behavior
nwauxio=0
convert_nocartesian=.false.
slice_nocartesian=.false.
saveprim=.true.
autoconvert=.false.
endian_swap=.false.
convert_type='vtuBCCmpi'
collapse_type='vti'
dxfiletype='lsb'
writew(1:nw)=.true.
writelevel(1:nlevelshi)=.true.
writespshift(1:ndim,1:2)=zero
level_io=-1
level_io_min=1
level_io_max=nlevelshi

! normalization of primitive variables: only for output
! note that normvar(0) is for length
! this scaling is optional, and must be set consistently if used
normvar(0:nw)=one
normt=one

! residual defaults
residmin=-1.0d0
residmax=bigdouble
typeresid='relative'

! AMR related defaults
mxnest=1
nbufferx1=0;nbufferx2=0;nbufferx3=0;
tol(1:nlevelshi)=0.1d0
tolratio(1:nlevelshi)=1.0d0/8.0d0
typegridfill='linear'
amrentropy=.false.
restrictprimitive=.true.
coarsenprimitive=.true.
prolongprimitive=.true.
typeprolonglimit='default'
errorestimate=3
flags(1:nflag_)=0
wflags(1:nflag_)=zero
flags(nflag_)=1
flags(1)=1
wflags(1)=one
logflag(1:nw)=.false.
amr_wavefilter(1:nlevelshi)=1.0d-2
skipfinestep=.false.
tfixgrid=bigdouble
itfixgrid=biginteger
ditregrid=1

! MHD specific defaults
B0field=.false.
Bdip=zero
Bquad=zero
Boct=zero
Busr=zero

! toggles the particle treatment:
use_particles = .false.

! IO defaults
itmax=biginteger

tmax=bigdouble
tmaxexact=.true.
dtmin=1.0d-10
typeparIO=0
nslices=0
collapse=.false.
collapseLevel=1
sliceascii=.true.
slice_type='vtuCC'


nshells=0
shell_type='vtu'

do ifile=1,nfile
   do isave=1,nsavehi
      tsave(isave,ifile)=bigdouble   ! t  of saves into the output files
      itsave(isave,ifile)=biginteger ! it of saves into the output files
   end do
   dtsave(ifile)=bigdouble           ! time between saves
   ditsave(ifile)=biginteger         ! timesteps between saves
   isavet(ifile)=1                   ! index for saves by t
   isaveit(ifile)=1                  ! index for saves by it
end do
typefilelog='default'
fileheadout = 'AMRVAC'
nwtf=0
neqpartf=0

! defaults for input 
firstprocess=.false.
resetgrid=.false.
changeglobals=.false.
treset=.false.
itreset=.false.
filenameout='data'
filenamelog='amrvac'

! Defaults for discretization methods
typeaverage='default'
tvdlfeps=one
BnormLF=.true.
nxdiffusehllc=0
flathllc=.false.
typeaxial='slab'
typecoord='covariant'
typespherical=1
slowsteps=-1
courantpar=0.8d0
typecourant='minimum'
dimsplit=.false.
typedimsplit='default'
typelimited='predictor'
typeinversion='1DW'
typeinversionresis='4DXIUNR'
mcbeta=1.4d0
betamin=1.d-2
useprimitive=.true.
typetvd='roe'
typeemf='none'
sourceimpl=.false.
sourceimplcycle=.false.
conduction=.false.
TCsaturate=.false.
TCphi=1.d0
energyonly=.false.
ncyclemax=1000
sourceparasts=.false.
parastsnu=0.001d0
ssplitdivb=.false.

ssplitresis=.false.
ssplituser=.false.
clean_init_divB='vecpot'
typeadvance='twostep'
do level=1,nlevelshi
   typefull1(level)='tvdlf'
   typepred1(level)='default'
   typelow1(level)='default'
   typelimiter1(level)='minmod'
   typegradlimiter1(level)='minmod'
enddo
flatcd=.false.
flatsh=.false.
typesourcesplit='sfs'
do iw=1,nw
   typeentropy(iw)='nul'      ! Entropy fix type
end do
dtdiffpar=0.5d0
dtTCpar=0.5d0
dtpar=-1.d0
time_accurate=.true.


! problem setup defaults
dxlone1=zero;dxlone2=zero;dxlone3=zero;
nxlone1=0;nxlone2=0;nxlone3=0;
iprob=1

! end defaults

! MPI reads from a file
open(unitpar,file=inifile,status='old')

! Start reading from standard input
primnames='default'

read(unitpar,filelist)

if(TRIM(primnames)=='default'.and.mype==0) write(uniterr,*) &
   'Warning in ReadParameters: primnames not given!'

if(firstprocess .and. snapshotini<0) call mpistop("Please restart from a &
   snapshot when firstprocess=T")

read(unitpar,savelist)

if (mype.eq.0) then
   do ifile=1,nfile
      if(dtsave(ifile)<bigdouble/2) write(unitterm,&
         '(" DTSAVE  for file",i2," =",g12.5)') ifile,dtsave(ifile)
      if(ditsave(ifile)<biginteger) write(unitterm,&
         '(" DITSAVE for file",i2," =",i10)') ifile,ditsave(ifile)
      if(tsave(1,ifile)==bigdouble.and.itsave(1,ifile)==biginteger&
         .and. dtsave(ifile)==bigdouble.and.ditsave(ifile)==biginteger&
         .and.mype==0) write(uniterr,*)'Warning in ReadParameters: ',&
          'No save condition for file ',ifile
   enddo
end if

! Consistency checks for the slices:
if (index(slice_type,'dat')>=1) sliceascii = .false.

do islice=1,nslices
   if(slicedir(islice) > ndim) write(uniterr,*)'Warning in ReadParameters: ',&
       'Slice ', islice,' direction',slicedir(islice),'larger than ndim=',ndim
   if(slicedir(islice) < 1) write(uniterr,*)'Warning in ReadParameters: ',&
       'Slice ', islice,' direction',slicedir(islice),&
      'too small, should be [',1,ndim,']'
end do



read(unitpar,stoplist)
if (mype.eq.0) then
   if(itmax<biginteger)       write(unitterm,*) 'ITMAX=',itmax
   if(tmax<bigdouble)         write(unitterm,*) 'TMAX=',tmax
   if(dtmin>smalldouble)      write(unitterm,*) 'DTMIN=',dtmin
end if

if(itmax==biginteger .and. tmax==bigdouble.and.mype==0) write(uniterr,&
   *) 'Warning in ReadParameters: itmax or tmax not given!'

if(residmin>=zero) then
   if (mype==0) write(unitterm,*)"SS computation with input value residmin"
   if(residmin<=smalldouble) call mpistop("Provide value for residual above &
      smalldouble")
end if

wnames='default'

read(unitpar,methodlist)

if(TRIM(wnames)=='default') call mpistop("Provide wnames and restart code")
wnameslog=wnames

do level=1,nlevelshi
   !if(typefull1(level)=='tvdlf1'.and.typeadvance=='twostep') &
   !   call mpistop(" tvdlf1 is onestep method, reset typeadvance=onestep!")
   !if(typefull1(level)=='hll1'.and.typeadvance=='twostep') &
   !   call mpistop(" hll1 is onestep method, reset typeadvance=onestep!")
   !if(typefull1(level)=='hllc1'.and.typeadvance=='twostep') &
   !   call mpistop(" hllc1 is onestep method, reset typeadvance=onestep!")
   !if(typefull1(level)=='hllcd1'.and.typeadvance=='twostep') &
   !   call mpistop(" hllcd1 is onestep method, reset typeadvance=onestep!")
   !if(typefull1(level)=='tvdmu1'.and.typeadvance=='twostep') &
   !   call mpistop(" tvdmu1 is onestep method, reset typeadvance=onestep!")
   if(typefull1(level)=='tvd'.and.typeadvance=='twostep') call mpistop(" tvd &
      is onestep method, reset typeadvance=onestep !")
   if(typefull1(level)=='tvd1'.and.typeadvance=='twostep') call mpistop(" &
      tvd1 is onestep method, reset typeadvance=onestep !")
   if(typefull1(level)=='tvd'.or.typefull1(level)=='tvd1')then 
      if(mype==0.and.(.not.dimsplit)) write(unitterm,*) &
         'Warning: setting dimsplit=T for tvd, as used for level=',level
      dimsplit=.true.
   endif

   if (typepred1(level)=='default') then
      select case (typefull1(level))
      case ('fd')
         typepred1(level)='fd'
      case ('tvdlf')
         typepred1(level)='tvdlf'
      case ('tvdlfpos')
         typepred1(level)='tvdlfpos'
      case ('hll')
         typepred1(level)='hll'
      case ('tvdlf1','tvdmu1','tvd1','tvd','hll1','hllc1', 'hlld1','hllcd1',&
         'hlldd1','nul','source')
         typepred1(level)='nul'
      case default
         call mpistop("No default predictor for this full step")
      end select
   end if
end do

if (convert_type == 'particles' .or. convert_type == 'particlesmpi') then
   use_particles = .true.
!   call mpistop('particle convert type requires use_particles=.true. in methodlist')
end if

select case (typeadvance)
case ("onestep")
   nstep=1
case ("twostep")
   nstep=2
case ("ImEx12")
   nstep=2
case ("threestep")
   nstep=3
case ("fourstep","rk4","jameson","ssprk43")
   nstep=4
case ("ssprk54")
   nstep=5
case default
   call mpistop("Unknown typeadvance")
end select


do level=1,nlevelshi
  if (typelow1(level)=='default') then
   select case (typefull1(level))
   case ('fd')
      typelow1(level)='fd'
   case ('tvdlf','tvdlf1','tvdmu','tvdmu1','tvd1','tvd','tvdlfpos')
      typelow1(level)='tvdlf1'
   case ('hll','hll1')
      typelow1(level)='hll1'
   case ('nul')
      typelow1(level)='nul'
   case ('source')
      typelow1(level)='source'
   case default
      call mpistop("No default typelow for this full step")
   end select
  end if
enddo

! Harmonize the parameters for dimensional splitting and source splitting
if(typedimsplit   =='default'.and.     dimsplit)   typedimsplit='xyyx'
if(typedimsplit   =='default'.and..not.dimsplit)   typedimsplit='unsplit'
dimsplit   = typedimsplit   /='unsplit'

if (index(typecoord,'covariant')>=1) then
   covariant=.true.
else
   covariant=.false.
end if

if (typeaxial=="slab".and..not.covariant) then
   slab=.true.
else
   slab=.false.
end if

if (typeaxial=='spherical') then
   if (dimsplit) then
      if(mype==0)print *,'Warning: spherical symmetry needs dimsplit=F, resetting'
      dimsplit=.false.
   end if
end if

if (typecoord=='default') then
   typecoord = typeaxial
end if

if (ndim==1) dimsplit=.false.
if (.not.dimsplit.and.ndim>1) then
   select case (typeadvance)
   case ("ssprk54","ssprk43","fourstep", "rk4", "threestep", "twostep",&
       "ImEx12")
      ! Runge-Kutta needs predictor
      typelimited="predictor"
      if(mype==0)print *,'typelimited to predictor for RK'
   end select
end if



if (typeemf .eq. 'none') typeemf = 'average' !Defaulting to average Balsara Spicer


if (B0field) then
   if(mype==0)print *,'B0+B1 split for MHD'
   if (.not.typephys=='mhd') call mpistop("B0+B1 split for MHD only")
end if

do level=1,nlevelshi
   type_limiter(level) = limiter_type(typelimiter1(level))
   type_gradient_limiter(level) = limiter_type(typegradlimiter1(level))
end do


if (any(typelimiter1(1:nlevelshi)== 'ppm').and.(flatsh.and.typephys&
   =='rho')) then
    call mpistop(" PPM with flatsh=.true. can not be used with typephys='rho' !")
end if
if (any(typelimiter1(1:nlevelshi)== 'ppm').and.(flatsh.and.typephys&
   =='hdadiab')) then
     call mpistop(" PPM with flatsh=.true. can not be used with typephys&
        ='hdadiab' !")
end if
if (any(typelimiter1(1:nlevelshi)== 'ppm').and.(flatcd.and.typephys&
   =='hdadiab')) then
     call mpistop(" PPM with flatcd=.true. can not be used with typephys&
        ='hdadiab' !")
end if
if (any(typelimiter1(1:nlevelshi)== 'ppm').and.(flatsh.and..not.useprimitive))&
    then
     call mpistop(" PPM with flatsh=.true. needs useprimitive=T!")
end if
if (any(typelimiter1(1:nlevelshi)== 'ppm').and.(flatcd.and..not.useprimitive))&
    then
     call mpistop(" PPM with flatcd=.true. needs useprimitive=T!")
end if

read(unitpar,boundlist)
do idim=1,ndim
   periodB(idim)=(any(typeB(:,2*idim-1:2*idim)=='periodic'))
   aperiodB(idim)=(any(typeB(:,2*idim-1:2*idim)=='aperiodic'))
   if (periodB(idim).or.aperiodB(idim)) then
      do iw=1,nw
         if (typeB(iw,2*idim-1) .ne. typeB(iw,2*idim)) call mpistop("Wrong &
            counterpart in periodic boundary")
         if (typeB(iw,2*idim-1) /= 'periodic' .and. typeB(iw,2*idim-1) &
            /= 'aperiodic') call mpistop("Each dimension should either have &
            all or no variables periodic, some can be aperiodic")
      end do
   end if
end do

if (any(typelimiter1(1:nlevelshi)=='ppm').and.(dixB<4)) then
    call mpistop(" PPM works only with dixB>=4 !")
end if

if (any(typelimiter1(1:nlevelshi)=='mp5') .and. (dixB<3)) then
 call mpistop("mp5 needs at at least 3 ghost cells! Set dixB=3 in boundlist.")
end if

if ( ixGhi1-ixGlo1+1-2*dixB .lt. 2*dixB.or. ixGhi2-ixGlo2+1-2*dixB .lt. &
   2*dixB.or. ixGhi3-ixGlo3+1-2*dixB .lt. 2*dixB) then
   ! without AMR you would actually get away with this, by why would you use so small blocks?
   call mpistop("For AMR you will need a physical cell to ghost-cell ratio of &
      >2")
end if


read(unitpar,amrlist)

select case (typeaxial)

case ("spherical")
   xprobmin2=xprobmin2*two*dpi;xprobmin3=xprobmin3*two*dpi
   xprobmax2=xprobmax2*two*dpi;xprobmax3=xprobmax3*two*dpi;

case ("cylindrical")
   if (1==2) then
      xprobmin1=xprobmin1*two*dpi;xprobmax1=xprobmax1*two*dpi;
   end if
   if (2==2) then
      xprobmin2=xprobmin2*two*dpi;xprobmax2=xprobmax2*two*dpi;
   end if
   if (3==2) then
      xprobmin3=xprobmin3*two*dpi;xprobmax3=xprobmax3*two*dpi;
   end if
end select

if(nxlone1>1 .and. mod(nxlone1,2)==0)then
   dxlone1=(xprobmax1-xprobmin1)/dble(nxlone1)
   if (mype==0) then
      write(unitterm,*)'Using ',nxlone1,' cells in dimension ',1
      write(unitterm,*)'level one dx(',1,')=',dxlone1
   end if
end if


if(nxlone2>1 .and. mod(nxlone2,2)==0)then
   dxlone2=(xprobmax2-xprobmin2)/dble(nxlone2)
   if (mype==0) then
      write(unitterm,*)'Using ',nxlone2,' cells in dimension ',2
      write(unitterm,*)'level one dx(',2,')=',dxlone2
   end if
end if


if(nxlone3>1 .and. mod(nxlone3,2)==0)then
   dxlone3=(xprobmax3-xprobmin3)/dble(nxlone3)
   if (mype==0) then
      write(unitterm,*)'Using ',nxlone3,' cells in dimension ',3
      write(unitterm,*)'level one dx(',3,')=',dxlone3
   end if
end if


if(dxlone1*dxlone2*dxlone3<smalldouble)then
   write(unitterm,*)'Wrong value(s) for level one dx:',dxlone1,dxlone2,dxlone3
   call mpistop("Reset nxlone or dxlone!")
endif
dx(1,1)=dxlone1;dx(2,1)=dxlone2;dx(3,1)=dxlone3;
if(mxnest>nlevelshi.or.mxnest<1)then
   write(unitterm,*)'Error: mxnest',mxnest,'>nlevelshi ',nlevelshi
   call mpistop("Reset nlevelshi and recompile!")
endif

if (flags(nflag_)>nw) then
   write(unitterm,*)'Error: flags(nw+1)=',flags(nw+1),'>nw ',nw
   call mpistop("Reset flags(nw+1)!")
end if
if (flags(nflag_)==0) errorestimate=0
if (flags(nflag_)<0) then
   if (mype==0) then
      write(unitterm,*) "flags(",nflag_,") can not be negative"
      call mpistop("")
   end if
end if
select case (errorestimate)
case (0)
   if (mype==0) write(unitterm,*)"Error estimation is user defined"
case (1)
   if (mype==0) write(unitterm,*)"Error estimation is richardson procedure"
case (2)
   if (mype==0) write(unitterm,*)"Error estimation is relative error"
case (3)
   if (mype==0) write(unitterm,*)"Error estimation is Lohner's scheme"
case (4)
   if (mype==0) write(unitterm,*)"Error estimation is Lohner's original &
      scheme"
case default
   call mpistop("Unknown error estimator, change errorestimate")
end select
if (B0field.and.errorestimate==1) then
   call mpistop("No Richardson procedure in combination with B0")
end if

if (tfixgrid<bigdouble/2.0d0) then
   if(mype==0)print*,'Warning, at time=',tfixgrid,'the grid will be fixed'
end if
if (itfixgrid<biginteger/2) then
   if(mype==0)print*,'Warning, at iteration=',itfixgrid,'the grid will be fixed'
end if
if (ditregrid>1) then
   if(mype==0)print*,'Note, Grid is reconstructed once every',ditregrid,'iterations'
end if

do islice=1,nslices
select case(slicedir(islice))
case(1)
   if(slicecoord(islice)<xprobmin1.or.slicecoord(islice)>xprobmax1) &
   write(uniterr,*)'Warning in ReadParameters: ', &
        'Slice ', islice, ' coordinate',slicecoord(islice),&
           'out of bounds for dimension ',slicedir(islice)

case(2)
   if(slicecoord(islice)<xprobmin2.or.slicecoord(islice)>xprobmax2) &
   write(uniterr,*)'Warning in ReadParameters: ', &
        'Slice ', islice, ' coordinate',slicecoord(islice),&
           'out of bounds for dimension ',slicedir(islice)

case(3)
   if(slicecoord(islice)<xprobmin3.or.slicecoord(islice)>xprobmax3) &
   write(uniterr,*)'Warning in ReadParameters: ', &
        'Slice ', islice, ' coordinate',slicecoord(islice),&
           'out of bounds for dimension ',slicedir(islice)

end select
end do

read(unitpar,paramlist)

if (dtpar>zero) time_accurate=.true.

if(.not.time_accurate) then
  if(residmin<=smalldouble .or. residmax==bigdouble) then
   if (mype==0) write(unitterm,*)"Non time_accurate SS computation needs &
      values residmin and residmax"
   call mpistop("Provide values for residual bounds in stoplist")
  end if
end if

close(unitpar)




if (mype==0) then
   print*,'Reading from inifile: ', trim(inifile)
   print*,'snapshotini         : ', snapshotini
   print*,'slicenext           : ', slicenext
   print*,'collapsenext        : ', collapsenext
   print*,'shellnext           : ', shellnext
   print*,'Filenameini         : ', trim(filenameini)
   print*,'Converting?         : ', convert
   print*,'                                                                '
endif

if (mype.eq.0) print*,'-----------------------------------------------------------------------------'

end subroutine readparameters
!=============================================================================
subroutine saveamrfile(ifile)

! following specific for Intel compiler and use on VIC3 with MPT
!DEC$ ATTRIBUTES NOINLINE :: write_snapshot

  use mod_slice, only: write_slice
  use mod_collapse, only: write_collapsed
  use mod_particles, only: write_particles_snapshot

  use mod_shell, only: write_shell 
use mod_amrvacdef
integer:: ifile
!-----------------------------------------------------------------------------
select case (ifile)
case (fileout_)


   if(endian_swap) typeparIO=-1
   select case (typeparIO)
      case (1) ! Parallel read and write
         call write_snapshot
      case (0,-2) ! Parallel read, serial write (-2)
         call write_snapshot_nopar
      case (-1) ! Serial read and write
         call write_snapshot_noparf
   end select



!opedit: now we can also convert directly and will when autoconvert is set in inifile: 
   if (autoconvert) call generate_plotfile
   if (use_particles) call write_particles_snapshot
case (fileslice_)
   call write_slice
case (filecollapse_)
   call write_collapsed
case (fileshell_)

    call write_shell


case (filelog_)
   select case (typefilelog)
   case ('default')
      call printlog_default
   case ('special')
      call printlog_special
   case default
      call mpistop("Error in SaveFile: Unknown typefilelog")
   end select
case (fileanalysis_)
  call write_analysis
case default
   write(*,*) 'No save method is defined for ifile=',ifile
   call mpistop("")
end select

! opedit: Flush stdout and stderr from time to time.
call flush(unitterm)

end subroutine saveamrfile
!=============================================================================
subroutine write_snapshot
use mod_forest
use mod_amrvacdef

integer :: file_handle, amode, igrid, Morton_no, iwrite
integer :: nx1,nx2,nx3
integer(kind=MPI_OFFSET_KIND) :: offset
integer, dimension(ngridshi) :: iorequest  , iorequest_stg
integer, dimension(MPI_STATUS_SIZE,ngridshi) :: iostatus  , iostatus_stg
integer, dimension(MPI_STATUS_SIZE) :: status
character(len=80) :: filename, line
logical, save :: firstsnapshot=.true.
!-----------------------------------------------------------------------------
if (firstsnapshot) then
   snapshot=snapshotnext
   firstsnapshot=.false.
end if

if (snapshot >= 10000) then
   if (mype==0) then
      write(*,*) "WARNING: Number of frames is limited to 10000 (0...9999),"
      write(*,*) "overwriting first frames"
   end if
   snapshot=0
end if

! generate filename
write(filename,"(a,i4.4,a)") TRIM(filenameout),snapshot,".dat"

if(mype==0) then
   open(unit=unitsnapshot,file=filename,status='replace')
   close(unit=unitsnapshot)
end if

amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
call MPI_FILE_OPEN(icomm,filename,amode,MPI_INFO_NULL,file_handle,ierrmpi)
iorequest=MPI_REQUEST_NULL

iorequest_stg=MPI_REQUEST_NULL 

iwrite=0
do Morton_no=Morton_start(mype),Morton_stop(mype)
   igrid=sfc_to_igrid(Morton_no)
   if (nwaux>0) then
      ! extra layer around mesh only for later averaging in convert
      ! set dxlevel value for use in gradient subroutine, 
      ! which might be used in getaux
      call set_tmpGlobals(igrid)
      call getaux(.true.,pw(igrid)%w,px(igrid)%x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,&
         ixGhi2,ixGhi3,ixMlo1-1,ixMlo2-1,ixMlo3-1,ixMhi1+1,ixMhi2+1,ixMhi3+1,&
         "write_snapshot")
   endif
   iwrite=iwrite+1
   offset=int(size_block_io  + size_block_stg_io ,kind=MPI_OFFSET_KIND) &
      *int(Morton_no-1,kind=MPI_OFFSET_KIND)
   call MPI_FILE_IWRITE_AT(file_handle,offset,pw(igrid)%w,1,type_block_io,&
       iorequest(iwrite),ierrmpi)
   
   offset=offset + int(size_block_io,kind=MPI_OFFSET_KIND)
   call MPI_FILE_IWRITE_AT(file_handle,offset,pws(igrid)%w,1,&
      type_block_staggered_io, &
        iorequest_stg(iwrite),ierrmpi)

end do

if (iwrite>0) call MPI_WAITALL(iwrite,iorequest,iostatus,ierrmpi)

if (iwrite>0) call MPI_WAITALL(iwrite,iorequest_stg,iostatus_stg,ierrmpi)


call MPI_FILE_CLOSE(file_handle,ierrmpi)
if (mype==0) then
   amode=ior(MPI_MODE_APPEND,MPI_MODE_WRONLY)
   call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode,MPI_INFO_NULL, file_handle,&
      ierrmpi)

   call write_forest(file_handle)

   nx1=ixMhi1-ixMlo1+1
   call MPI_FILE_WRITE(file_handle,nx1,1,MPI_INTEGER,status,ierrmpi)
   nx2=ixMhi2-ixMlo2+1
   call MPI_FILE_WRITE(file_handle,nx2,1,MPI_INTEGER,status,ierrmpi)
   nx3=ixMhi3-ixMlo3+1
   call MPI_FILE_WRITE(file_handle,nx3,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,eqpar,neqpar+nspecialpar,&
       MPI_DOUBLE_PRECISION,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,nleafs,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,levmax,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,ndim,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,ndir,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,nw,1,MPI_INTEGER,status,ierrmpi)
!{#IFDEF STAGGERED
   call MPI_FILE_WRITE(file_handle,nws,1,MPI_INTEGER,status,ierrmpi)
!}
   call MPI_FILE_WRITE(file_handle,neqpar+nspecialpar,1,MPI_INTEGER,status,&
      ierrmpi)
   call MPI_FILE_WRITE(file_handle,it,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,t,1,MPI_DOUBLE_PRECISION,status,ierrmpi)

   call MPI_FILE_CLOSE(file_handle,ierrmpi)
end if
snapshot=snapshot+1

end subroutine write_snapshot
!=============================================================================
subroutine write_snapshot_nopar
use mod_forest
use mod_amrvacdef


integer :: file_handle, amode, igrid, Morton_no, iwrite
integer :: nx1,nx2,nx3

integer(kind=MPI_OFFSET_KIND) :: offset

integer, allocatable :: iorecvstatus(:,:),ioastatus(:,:)
integer, allocatable :: igrecvstatus(:,:)
integer, allocatable :: igrid_recv(:) 

integer, allocatable :: iorecvstatus_stg(:,:),ioastatus_stg(:,:)

integer, dimension(MPI_STATUS_SIZE) :: status

integer  :: ipe,insend,inrecv,nrecv,nwrite
character(len=80) :: filename, line
logical, save :: firstsnapshot=.true.
!-----------------------------------------------------------------------------

call MPI_BARRIER(icomm,ierrmpi)

if (firstsnapshot) then
   snapshot=snapshotnext
   firstsnapshot=.false.
end if

if (snapshot >= 10000) then
   if (mype==0) then
      write(*,*) "WARNING: Number of frames is limited to 10000 (0...9999),"
      write(*,*) "overwriting first frames"
   end if
   snapshot=0
end if

nrecv=0
inrecv=0
nwrite=0
insend=0
iwrite=0

if (mype /= 0) then
 do Morton_no=Morton_start(mype),Morton_stop(mype)
    igrid=sfc_to_igrid(Morton_no)
    itag=Morton_no

    itag_stg=Morton_no+ngridshi

    insend=insend+1
    if (nwaux>0) then
      ! extra layer around mesh only for later averaging in convert
      ! set dxlevel value for use in gradient subroutine, 
      ! which might be used in getaux
      saveigrid=igrid
      dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
      dxlevel(3)=rnode(rpdx3_,igrid);
      mygeo =>pgeo(igrid)
      if(covariant)myM => mygeo%m
      if (B0field) then
         myB0_cell => pB0_cell(igrid)
         myB0_face1 => pB0_face1(igrid)
         myB0_face2 => pB0_face2(igrid)
         myB0_face3 => pB0_face3(igrid)
      end if
      call getaux(.true.,pw(igrid)%w,px(igrid)%x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,&
         ixGhi2,ixGhi3,ixMlo1-1,ixMlo2-1,ixMlo3-1,ixMhi1+1,ixMhi2+1,ixMhi3+1,&
         "write_snapshot")
    endif
    call MPI_SEND(igrid,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
    call MPI_SEND(pw(igrid)%w,1,type_block_io, 0,itag,icomm,ierrmpi)
    
    call MPI_SEND(pws(igrid)%w,1,type_block_staggered_io, 0,itag_stg,icomm,&
       ierrmpi)
    
 end do
else 
 ! mype==0
 nwrite=(Morton_stop(0)-Morton_start(0)+1)

 ! master processor writes out
 write(filename,"(a,i4.4,a)") TRIM(filenameout),snapshot,".dat"

 open(unit=unitsnapshot,file=filename,status='replace')
 close(unitsnapshot)

 amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
 call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode,MPI_INFO_NULL,file_handle,&
    ierrmpi)


 ! writing his local data first
 do Morton_no=Morton_start(0),Morton_stop(0)
   igrid=sfc_to_igrid(Morton_no)
   iwrite=iwrite+1
   if (nwaux>0) then
      ! extra layer around mesh only for later averaging in convert
      ! set dxlevel value for use in gradient subroutine,
      ! which might be used in getaux
      saveigrid=igrid
      dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
      dxlevel(3)=rnode(rpdx3_,igrid);
      mygeo =>pgeo(igrid)
      if(covariant)myM => mygeo%m
      if (B0field) then
         myB0_cell => pB0_cell(igrid)
         myB0_face1 => pB0_face1(igrid)
         myB0_face2 => pB0_face2(igrid)
         myB0_face3 => pB0_face3(igrid)
      end if
      call getaux(.true.,pw(igrid)%w,px(igrid)%x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,&
         ixGhi2,ixGhi3,ixMlo1-1,ixMlo2-1,ixMlo3-1,ixMhi1+1,ixMhi2+1,ixMhi3+1,&
         "write_snapshot")
   endif

   offset=int(size_block_io  + size_block_stg_io ,kind=MPI_OFFSET_KIND) &
      *int(Morton_no-1,kind=MPI_OFFSET_KIND)
   call MPI_FILE_WRITE_AT(file_handle,offset,pw(igrid)%w,1,type_block_io,&
       MPI_STATUSES_IGNORE,ierrmpi)
   
   offset=offset + int(size_block_io,kind=MPI_OFFSET_KIND)
   call MPI_FILE_WRITE_AT(file_handle,offset,pws(igrid)%w,1,&
      type_block_staggered_io, &
        MPI_STATUSES_IGNORE,ierrmpi)
   
   

 end do
 ! write data communicated from other processors
 ! opedit: take care only to receive if other processors have data
 ! opedit: also allocate nrecv based on actual nleafs instead of Morton_stop(npe-1)
 ! opedit: for the case that last processor has no data
 if(npe>1 .and. nleafs.gt.Morton_stop(0))then
    nrecv = (nleafs-Morton_start(1)+1)
    inrecv=0
    allocate(igrid_recv(nrecv))
    allocate(igrecvstatus(MPI_STATUS_SIZE,nrecv),iorecvstatus(MPI_STATUS_SIZE,&
       nrecv))
    allocate(ioastatus(MPI_STATUS_SIZE,nrecv))
    
    allocate(iorecvstatus_stg(MPI_STATUS_SIZE,nrecv))
    allocate(ioastatus_stg(MPI_STATUS_SIZE,nrecv))
    

    do ipe =1, npe-1
       do Morton_no=Morton_start(ipe),Morton_stop(ipe)
          iwrite=iwrite+1
          itag=Morton_no
          
          itag_stg=Morton_no+ngridshi
          
          inrecv=inrecv+1
          call MPI_RECV(igrid_recv(inrecv),1,MPI_INTEGER, ipe,itag,icomm,&
             igrecvstatus(:,inrecv),ierrmpi)

          allocate(pwio(igrid_recv(inrecv))%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             ixGlo3:ixGhi3,1:nw))
          call MPI_RECV(pwio(igrid_recv(inrecv))%w,1,type_block_io,ipe,itag,&
             icomm,iorecvstatus(:,inrecv),ierrmpi)
          
          allocate(pwsio(igrid_recv(inrecv))%w(ixGlo1-1:ixGhi1,&
             ixGlo2-1:ixGhi2,ixGlo3-1:ixGhi3,1:ndim))
          call MPI_RECV(pwsio(igrid_recv(inrecv))%w,1,type_block_staggered_io,&
             ipe,itag_stg,icomm,&
               iorecvstatus_stg(:,inrecv),ierrmpi)
          
          offset=int(size_block_io  + size_block_stg_io ,kind&
             =MPI_OFFSET_KIND) *int(Morton_no-1,kind=MPI_OFFSET_KIND)
          call MPI_FILE_WRITE_AT(file_handle,offset,pwio(igrid_recv&
             (inrecv))%w,1,type_block_io,ioastatus(:,inrecv),ierrmpi)
          
          offset=offset + int(size_block_io,kind=MPI_OFFSET_KIND)
          call MPI_FILE_WRITE_AT(file_handle,offset,pwsio(igrid_recv&
             (inrecv))%w,1,&
               type_block_staggered_io,ioastatus_stg(:,inrecv),ierrmpi)
          
          deallocate(pwio(igrid_recv(inrecv))%w)
          
          deallocate(pwsio(igrid_recv(inrecv))%w)
          
       end do
    end do
    deallocate(igrecvstatus,iorecvstatus,ioastatus,igrid_recv)
    
    deallocate(iorecvstatus_stg,ioastatus_stg)
    
 end if
end if


call MPI_FILE_CLOSE(file_handle,ierrmpi)


if (mype==0) then
   amode=ior(MPI_MODE_APPEND,MPI_MODE_WRONLY)
   call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode,MPI_INFO_NULL, file_handle,&
      ierrmpi)

   call write_forest(file_handle)

   nx1=ixMhi1-ixMlo1+1
   call MPI_FILE_WRITE(file_handle,nx1,1,MPI_INTEGER,status,ierrmpi)
   nx2=ixMhi2-ixMlo2+1
   call MPI_FILE_WRITE(file_handle,nx2,1,MPI_INTEGER,status,ierrmpi)
   nx3=ixMhi3-ixMlo3+1
   call MPI_FILE_WRITE(file_handle,nx3,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,eqpar,neqpar+nspecialpar,&
       MPI_DOUBLE_PRECISION,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,nleafs,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,levmax,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,ndim,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,ndir,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,nw,1,MPI_INTEGER,status,ierrmpi)
!{#IFDEF STAGGERED
   call MPI_FILE_WRITE(file_handle,nws,1,MPI_INTEGER,status,ierrmpi)
!}
   call MPI_FILE_WRITE(file_handle,neqpar+nspecialpar,1,MPI_INTEGER,status,&
      ierrmpi)
   call MPI_FILE_WRITE(file_handle,it,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,t,1,MPI_DOUBLE_PRECISION,status,ierrmpi)

   call MPI_FILE_CLOSE(file_handle,ierrmpi)
end if

snapshot=snapshot+1

call MPI_BARRIER(icomm,ierrmpi)

end subroutine write_snapshot_nopar
!=============================================================================
subroutine write_snapshot_noparf
use mod_forest
use mod_amrvacdef

integer :: igrid, Morton_no
integer :: nx1,nx2,nx3

integer, allocatable :: iostatus(:,:),iorecvstatus(:,:),ioastatus(:,:)
integer, allocatable :: igrecvstatus(:,:)
integer, allocatable :: iorequest(:),igrid_recv(:) 

integer, allocatable :: iostatus_stg(:,:),iorecvstatus_stg(:,:),&
   ioastatus_stg(:,:)
integer, allocatable :: iorequest_stg(:)

integer  :: ipe,insend,inrecv,nrecv,nwrite
character(len=80) :: filename, line
logical, save :: firstsnapshot=.true.
!-----------------------------------------------------------------------------
call MPI_BARRIER(icomm,ierrmpi)

if (firstsnapshot) then
   snapshot=snapshotnext
   firstsnapshot=.false.
end if

if (snapshot >= 10000) then
   if (mype==0) then
      write(*,*) "WARNING: Number of frames is limited to 10000 (0...9999),"
      write(*,*) "overwriting first frames"
   end if
   snapshot=0
end if

nrecv=0
inrecv=0
nwrite=0
insend=0

if (mype /= 0) then
 do Morton_no=Morton_start(mype),Morton_stop(mype)
    igrid=sfc_to_igrid(Morton_no)
    itag=Morton_no

    itag_stg=Morton_no+ngridshi

    insend=insend+1
    if (nwaux>0) then
      ! extra layer around mesh only for later averaging in convert
      ! set dxlevel value for use in gradient subroutine, 
      ! which might be used in getaux
      saveigrid=igrid
      dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
      dxlevel(3)=rnode(rpdx3_,igrid);
      mygeo =>pgeo(igrid)
      if(covariant)myM => mygeo%m
      if (B0field) then
         myB0_cell => pB0_cell(igrid)
         myB0_face1 => pB0_face1(igrid)
         myB0_face2 => pB0_face2(igrid)
         myB0_face3 => pB0_face3(igrid)
      end if
      call getaux(.true.,pw(igrid)%w,px(igrid)%x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,&
         ixGhi2,ixGhi3,ixMlo1-1,ixMlo2-1,ixMlo3-1,ixMhi1+1,ixMhi2+1,ixMhi3+1,&
         "write_snapshot")
    endif
    call MPI_SEND(igrid,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
    call MPI_SEND(pw(igrid)%w,1,type_block_io, 0,itag,icomm,ierrmpi)
    
    call MPI_SEND(pws(igrid)%w,1,type_block_staggered_io, 0,itag_stg,icomm,&
       ierrmpi)
    
 end do
else 
 ! mype==0
 nwrite=(Morton_stop(0)-Morton_start(0)+1)
 allocate(iorequest(nwrite),iostatus(MPI_STATUS_SIZE,nwrite))

 allocate(iorequest_stg(nwrite),iostatus_stg(MPI_STATUS_SIZE,nwrite))

 iorequest=MPI_REQUEST_NULL

 iorequest_stg=MPI_REQUEST_NULL

 ! master processor writes out
 write(filename,"(a,i4.4,a)") TRIM(filenameout),snapshot,".dat"
 if(endian_swap) then
  
   open(unit=unitsnapshot,file=filename,form='unformatted',access='stream',&
        status='replace',convert='BIG_ENDIAN')
  
  
 else
   open(unit=unitsnapshot,file=filename,form='unformatted',access&
      ='stream',status='replace')
 end if
 ! writing his local data first
 do Morton_no=Morton_start(0),Morton_stop(0)
   igrid=sfc_to_igrid(Morton_no)
   if (nwaux>0) then
      ! extra layer around mesh only for later averaging in convert
      ! set dxlevel value for use in gradient subroutine,
      ! which might be used in getaux
      saveigrid=igrid
      dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
      dxlevel(3)=rnode(rpdx3_,igrid);
      mygeo =>pgeo(igrid)
      if(covariant)myM => mygeo%m
      if (B0field) then
         myB0_cell => pB0_cell(igrid)
         myB0_face1 => pB0_face1(igrid)
         myB0_face2 => pB0_face2(igrid)
         myB0_face3 => pB0_face3(igrid)
      end if
      call getaux(.true.,pw(igrid)%w,px(igrid)%x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,&
         ixGhi2,ixGhi3,ixMlo1-1,ixMlo2-1,ixMlo3-1,ixMhi1+1,ixMhi2+1,ixMhi3+1,&
         "write_snapshot")
   endif
   write(unitsnapshot) pw(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,&
      1:nw)

   write(unitsnapshot) pws(igrid)%w(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,&
      ixMlo3-1:ixMhi3,1:nws)

 end do
 ! write data communicated from other processors
 if(npe>1)then
  nrecv=(Morton_stop(npe-1)-Morton_start(1)+1)
  inrecv=0
  allocate(igrid_recv(nrecv))
  allocate(igrecvstatus(MPI_STATUS_SIZE,nrecv),iorecvstatus(MPI_STATUS_SIZE,&
     nrecv))
  allocate(ioastatus(MPI_STATUS_SIZE,nrecv))
  
  allocate(iorecvstatus_stg(MPI_STATUS_SIZE,nrecv))
  allocate(ioastatus_stg(MPI_STATUS_SIZE,nrecv))
  

  do ipe =1, npe-1
   do Morton_no=Morton_start(ipe),Morton_stop(ipe)
     itag=Morton_no

     itag_stg=Morton_no+ngridshi

     inrecv=inrecv+1
     call MPI_RECV(igrid_recv(inrecv),1,MPI_INTEGER, ipe,itag,icomm,&
        igrecvstatus(:,inrecv),ierrmpi)

     allocate(pwio(igrid_recv(inrecv))%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
        ixGlo3:ixGhi3,1:nw))
     call MPI_RECV(pwio(igrid_recv(inrecv))%w,1,type_block_io,ipe,itag,icomm,&
        iorecvstatus(:,inrecv),ierrmpi)

     allocate(pwsio(igrid_recv(inrecv))%w(ixGlo1-1:ixGhi1,ixGlo2-1:ixGhi2,&
        ixGlo3-1:ixGhi3,1:ndim))
     call MPI_RECV(pwsio(igrid_recv(inrecv))%w,1,type_block_staggered_io,ipe,&
        itag_stg,icomm,&
                   iorecvstatus_stg(:,inrecv),ierrmpi)

     write(unitsnapshot) pwio(igrid_recv(inrecv))%w(ixMlo1:ixMhi1,&
        ixMlo2:ixMhi2,ixMlo3:ixMhi3,1:nw)

     write(unitsnapshot) pwsio(igrid_recv(inrecv))%w(ixMlo1-1:ixMhi1,&
        ixMlo2-1:ixMhi2,ixMlo3-1:ixMhi3,1:nws)

     deallocate(pwio(igrid_recv(inrecv))%w)

     deallocate(pwsio(igrid_recv(inrecv))%w)

   end do
  end do
  deallocate(igrecvstatus,iorecvstatus,ioastatus,igrid_recv)
  
  deallocate(iorecvstatus_stg,ioastatus_stg)
  
 end if
end if

if(nwrite>0)then
  call MPI_WAITALL(nwrite,iorequest,iostatus,ierrmpi) 

  call MPI_WAITALL(nwrite,iorequest_stg,iostatus_stg,ierrmpi)

  if(mype==0) then
    deallocate(iorequest,iostatus)
    
    deallocate(iorequest_stg, iostatus_stg)
  end if
end if

if(mype==0) then
  call write_forest(unitsnapshot)
  nx1=ixMhi1-ixMlo1+1
  write(unitsnapshot) nx1
  nx2=ixMhi2-ixMlo2+1
  write(unitsnapshot) nx2
  nx3=ixMhi3-ixMlo3+1
  write(unitsnapshot) nx3
  write(unitsnapshot) eqpar
  write(unitsnapshot) nleafs
  write(unitsnapshot) levmax 
  write(unitsnapshot) ndim 
  write(unitsnapshot) ndir
  write(unitsnapshot) nw
  write(unitsnapshot) nws !! Staggered variables
  write(unitsnapshot) neqpar+nspecialpar
  write(unitsnapshot) it
  write(unitsnapshot) t
  close(unitsnapshot)
end if

snapshot=snapshot+1

call MPI_BARRIER(icomm,ierrmpi)

end subroutine write_snapshot_noparf
!=============================================================================
subroutine read_snapshot
use mod_forest
use mod_amrvacdef

integer :: file_handle, amode, igrid, Morton_no, iread
integer :: levmaxini, ndimini, ndirini, nwini, nwsini, neqparini, nxini1,&
   nxini2,nxini3

integer(kind=MPI_ADDRESS_KIND) :: size_double, size_int, lb

integer(kind=MPI_OFFSET_KIND) :: offset
integer, dimension(ngridshi) :: iorequest  , iorequest_stg
integer, dimension(MPI_STATUS_SIZE,ngridshi) :: iostatus  , iostatus_stg
integer, dimension(MPI_STATUS_SIZE) :: status
character(len=80) :: filename
logical :: fexist
integer(kind=8) :: file_size
!-----------------------------------------------------------------------------
!!!call MPI_BARRIER(icomm,ierrmpi)
! generate filename
write(filename,"(a,i4.4,a)") TRIM(filenameini),snapshotini,".dat"

if(mype==0) then
  inquire(file=filename,exist=fexist)
  if(.not.fexist) call mpistop(filename//"as an input snapshot file is not found !")
endif

amode=MPI_MODE_RDONLY
call MPI_FILE_OPEN(icomm,filename,amode,MPI_INFO_NULL,file_handle,ierrmpi)


!call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,size_double,ierrmpi)
!call MPI_TYPE_EXTENT(MPI_INTEGER,size_int,ierrmpi)
call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,size_double,ierrmpi)
call MPI_TYPE_GET_EXTENT(MPI_INTEGER,lb,size_int,ierrmpi)

!offset=-int( 8*size_int+size_double,kind=MPI_OFFSET_KIND)

inquire(file=filename,exist=fexist,size=file_size)
offset= file_size-int(8*size_int+size_double)

!call MPI_FILE_SEEK_SHARED(file_handle,offset,MPI_SEEK_END,ierrmpi)
!call MPI_FILE_SEEK(file_handle,offset,MPI_SEEK_END,ierrmpi)

call MPI_FILE_SEEK(file_handle,offset,MPI_SEEK_SET,ierrmpi)

call MPI_FILE_READ_ALL(file_handle,nleafs,1,MPI_INTEGER,status,ierrmpi)
nleafs_active = nleafs
call MPI_FILE_READ_ALL(file_handle,levmaxini,1,MPI_INTEGER,status,ierrmpi)
call MPI_FILE_READ_ALL(file_handle,ndimini,1,MPI_INTEGER,status,ierrmpi)
call MPI_FILE_READ_ALL(file_handle,ndirini,1,MPI_INTEGER,status,ierrmpi)
call MPI_FILE_READ_ALL(file_handle,nwini,1,MPI_INTEGER,status,ierrmpi)
!{#IFDEF STAGGERED
call MPI_FILE_READ_ALL(file_handle,nwsini,1,MPI_INTEGER,status,ierrmpi)
!}
call MPI_FILE_READ_ALL(file_handle,neqparini,1,MPI_INTEGER,status,ierrmpi)
call MPI_FILE_READ_ALL(file_handle,it,1,MPI_INTEGER,status,ierrmpi)
call MPI_FILE_READ_ALL(file_handle,t,1,MPI_DOUBLE_PRECISION,status,ierrmpi)

! check if settings are suitable for restart
if (levmaxini>mxnest) then
   if (mype==0) write(*,*) "number of levels in restart file = ",levmaxini
   if (mype==0) write(*,*) "mxnest = ",mxnest
   call mpistop("mxnest should be at least number of levels in restart file")
end if
if (ndimini/=ndim) then
   if (mype==0) write(*,*) "ndim in restart file = ",ndimini
   if (mype==0) write(*,*) "ndim = ",ndim
   call mpistop("reset ndim to ndim in restart file")
end if
if (ndirini/=ndir) then
   if (mype==0) write(*,*) "ndir in restart file = ",ndirini
   if (mype==0) write(*,*) "ndir = ",ndir
   call mpistop("reset ndir to ndir in restart file")
end if
if (nw/=nwini) then
   if (mype==0) write(*,*) "nw=",nw," and nw in restart file=",nwini
   call mpistop("currently, changing nw at restart is not allowed")
end if

if (nws/=nwsini) then
   if (mype==0) write(*,*) "nws=",nws," and nws in restart file=",nwsini
   call mpistop("currently, changing nws at restart is not allowed")
end if

offset=offset-int(ndimini*size_int+neqparini*size_double,kind=MPI_OFFSET_KIND)
!call MPI_FILE_SEEK_SHARED(file_handle,offset,MPI_SEEK_END,ierrmpi)
!call MPI_FILE_SEEK(file_handle,offset,MPI_SEEK_END,ierrmpi)

call MPI_FILE_SEEK(file_handle,offset,MPI_SEEK_SET,ierrmpi)
call MPI_FILE_READ_ALL(file_handle,nxini1,1,MPI_INTEGER,status,ierrmpi)
call MPI_FILE_READ_ALL(file_handle,nxini2,1,MPI_INTEGER,status,ierrmpi)
call MPI_FILE_READ_ALL(file_handle,nxini3,1,MPI_INTEGER,status,ierrmpi)
if (ixGhi1/=nxini1+2*dixB.or.ixGhi2/=nxini2+2*dixB.or.ixGhi3&
   /=nxini3+2*dixB) then
   if (mype==0) write(*,*) "Error: reset resolution to ",nxini1+2*dixB,&
      nxini2+2*dixB,nxini3+2*dixB
   call mpistop("change with setamrvac")
end if
neqparini=min(neqparini,neqpar+nspecialpar)
call MPI_FILE_READ_ALL(file_handle,eqpar,neqparini, MPI_DOUBLE_PRECISION,&
   status,ierrmpi)


call read_forest(file_handle)

iorequest=MPI_REQUEST_NULL

iorequest_stg=MPI_REQUEST_NULL 

iread=0
do Morton_no=Morton_start(mype),Morton_stop(mype)
   igrid=sfc_to_igrid(Morton_no)
   call alloc_node(igrid)
   iread=iread+1
   offset=int(size_block_io  + size_block_stg_io ,kind=MPI_OFFSET_KIND) &
      *int(Morton_no-1,kind=MPI_OFFSET_KIND)
   call MPI_FILE_IREAD_AT(file_handle,offset,pw(igrid)%w,1,type_block_io,&
       iorequest(iread),ierrmpi)
   
   offset=offset + int(size_block_io,kind=MPI_OFFSET_KIND)
   call MPI_FILE_IREAD_AT(file_handle,offset,pws(igrid)%w,1,&
      type_block_staggered_io, &
        iorequest_stg(iread),ierrmpi)

end do

if (iread>0) call MPI_WAITALL(iread,iorequest,iostatus,ierrmpi)

if (iread>0) call MPI_WAITALL(iread,iorequest_stg,iostatus_stg,ierrmpi)


call MPI_FILE_CLOSE(file_handle,ierrmpi)

!!!call MPI_BARRIER(icomm,ierrmpi)
end subroutine read_snapshot
!=============================================================================
subroutine read_snapshot_nopar
use mod_forest
use mod_amrvacdef

double precision :: wio(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nw)  ,&
   wsio(ixGlo1-1:ixGhi1,ixGlo2-1:ixGhi2,ixGlo3-1:ixGhi3,1:nws)
integer :: file_handle, amode, igrid, Morton_no, iread
integer :: levmaxini, ndimini, ndirini, nwini, nwsini, neqparini, nxini1,&
   nxini2,nxini3

integer(kind=MPI_ADDRESS_KIND) :: size_double, size_int, lb

integer(kind=MPI_OFFSET_KIND) :: offset
integer, dimension(ngridshi) :: iorequest  , iorequest_stg
integer, dimension(MPI_STATUS_SIZE) :: status
integer, dimension(MPI_STATUS_SIZE) :: iostatus  , iostatus_stg

integer, allocatable :: iorecvstatus(:,:)  , iorecvstatus_stg(:,:)
integer :: ipe,inrecv,nrecv
integer ::  sendini(8+3)
character(len=80) :: filename
logical :: fexist
integer(kind=8) :: file_size
!-----------------------------------------------------------------------------
!!!call MPI_BARRIER(icomm,ierrmpi)
! generate filename
write(filename,"(a,i4.4,a)") TRIM(filenameini),snapshotini,".dat"
if (mype==0) then
 inquire(file=filename,exist=fexist,size=file_size)
 if(.not.fexist) call mpistop(filename//"as an input snapshot file is not found !")
 amode=MPI_MODE_RDONLY
 call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode,MPI_INFO_NULL,file_handle,&
    ierrmpi)

 !call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,size_double,ierrmpi)
 !call MPI_TYPE_EXTENT(MPI_INTEGER,size_int,ierrmpi)
 call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,size_double,ierrmpi)
 call MPI_TYPE_GET_EXTENT(MPI_INTEGER,lb,size_int,ierrmpi)

 !offset=-int( 8*size_int+size_double,kind=MPI_OFFSET_KIND)
 offset= file_size-int(8*size_int+size_double)
 !call MPI_FILE_SEEK_SHARED(file_handle,offset,MPI_SEEK_END,ierrmpi)
 !call MPI_FILE_SEEK(file_handle,offset,MPI_SEEK_END,ierrmpi)
 call MPI_FILE_SEEK(file_handle,offset,MPI_SEEK_SET,ierrmpi)

 call MPI_FILE_READ(file_handle,nleafs,1,MPI_INTEGER,status,ierrmpi)
 nleafs_active = nleafs
 call MPI_FILE_READ(file_handle,levmaxini,1,MPI_INTEGER,status,ierrmpi)
 call MPI_FILE_READ(file_handle,ndimini,1,MPI_INTEGER,status,ierrmpi)
 call MPI_FILE_READ(file_handle,ndirini,1,MPI_INTEGER,status,ierrmpi)
 call MPI_FILE_READ(file_handle,nwini,1,MPI_INTEGER,status,ierrmpi)
!{#IFDEF STAGGERED
 call MPI_FILE_READ(file_handle,nwsini,1,MPI_INTEGER,status,ierrmpi)
!}
 call MPI_FILE_READ(file_handle,neqparini,1,MPI_INTEGER,status,ierrmpi)
 call MPI_FILE_READ(file_handle,it,1,MPI_INTEGER,status,ierrmpi)
 call MPI_FILE_READ(file_handle,t,1,MPI_DOUBLE_PRECISION,status,ierrmpi)
 ! check if settings are suitable for restart
 if (levmaxini>mxnest) then
    write(*,*) "number of levels in restart file = ",levmaxini
    write(*,*) "mxnest = ",mxnest
    call mpistop("mxnest should be at least number of levels in restart file")
 end if
 if (ndimini/=ndim) then
    write(*,*) "ndim in restart file = ",ndimini
    write(*,*) "ndim = ",ndim
    call mpistop("reset ndim to ndim in restart file")
 end if
 if (ndirini/=ndir) then
    write(*,*) "ndir in restart file = ",ndirini
    write(*,*) "ndir = ",ndir
    call mpistop("reset ndir to ndir in restart file")
 end if
 if (nw/=nwini) then
    write(*,*) "nw=",nw," and nw in restart file=",nwini
    call mpistop("currently, changing nw at restart is not allowed")
 end if

 if (nws/=nwsini) then
   if (mype==0) write(*,*) "nws=",nws," and nws in restart file=",nwsini
   call mpistop("currently, changing nws at restart is not allowed")
 end if

 offset=offset-int(ndimini*size_int+neqparini*size_double,kind&
    =MPI_OFFSET_KIND)
 !call MPI_FILE_SEEK_SHARED(file_handle,offset,MPI_SEEK_END,ierrmpi)
 !call MPI_FILE_SEEK(file_handle,offset,MPI_SEEK_END,ierrmpi)
 call MPI_FILE_SEEK(file_handle,offset,MPI_SEEK_SET,ierrmpi)

 call MPI_FILE_READ(file_handle,nxini1,1,MPI_INTEGER,status,ierrmpi)
 call MPI_FILE_READ(file_handle,nxini2,1,MPI_INTEGER,status,ierrmpi)
 call MPI_FILE_READ(file_handle,nxini3,1,MPI_INTEGER,status,ierrmpi)
 if (ixGhi1/=nxini1+2*dixB.or.ixGhi2/=nxini2+2*dixB.or.ixGhi3&
    /=nxini3+2*dixB) then
    write(*,*) "Error: reset resolution to ",nxini1+2*dixB,nxini2+2*dixB,&
       nxini3+2*dixB
    call mpistop("change with setamrvac")
 end if
 neqparini=min(neqparini,neqpar+nspecialpar)
 call MPI_FILE_READ(file_handle,eqpar,neqparini, MPI_DOUBLE_PRECISION,status,&
    ierrmpi)
end if

! broadcast the global parameters first
if (npe>1) then
  if (mype==0) then
     sendini=(/nleafs,levmaxini,ndimini,ndirini,nwini,  nwsini, neqparini,it ,&
        nxini1,nxini2,nxini3 /)
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! nws is defined always, not only when staggered is activated. !!!
! This might cause bugs.                                       !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  call MPI_BCAST(sendini,{#IFDEF STAGGERED 1+} 7+^ND,MPI_INTEGER,0,icomm,ierrmpi)
  call MPI_BCAST(sendini,8+3,MPI_INTEGER,0,icomm,ierrmpi)
  nleafs=sendini(1);levmaxini=sendini(2);ndimini=sendini(3);
  ndirini=sendini(4);nwini=sendini(5);

  nwsini=sendini(6);neqparini=sendini(7);
  it=sendini(8);
  nxini1=sendini(8+1);nxini2=sendini(8+2);nxini3=sendini(8+3);

  
  nleafs_active = nleafs
  call MPI_BCAST(t,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
  call MPI_BCAST(eqpar,neqparini,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
end if

call read_forest(file_handle)

if (mype==0)then
   iread=0
!!! HO: What is this for? (iorequest)
   iorequest=MPI_REQUEST_NULL
   
   iorequest_stg=MPI_REQUEST_NULL
   do Morton_no=Morton_start(0),Morton_stop(0)
      igrid=sfc_to_igrid(Morton_no)
      call alloc_node(igrid)
      iread=iread+1
      offset=int(size_block_io  + size_block_stg_io ,kind=MPI_OFFSET_KIND) &
         *int(Morton_no-1,kind=MPI_OFFSET_KIND)
      call MPI_FILE_READ_AT(file_handle,offset,pw(igrid)%w,1,type_block_io,&
          iostatus,ierrmpi)
      
      offset=offset + int(size_block_io,kind=MPI_OFFSET_KIND)
      call MPI_FILE_READ_AT(file_handle,offset,pws(igrid)%w,1,&
         type_block_staggered_io, &
                            iostatus_stg,ierrmpi)
       
   end do
   if (npe>1) then
    do ipe=1,npe-1
     do Morton_no=Morton_start(ipe),Morton_stop(ipe)
       iread=iread+1
       itag=Morton_no
       
       itag_stg=Morton_no+ngridshi
       offset=int(size_block_io  + size_block_stg_io ,kind=MPI_OFFSET_KIND) &
          *int(Morton_no-1,kind=MPI_OFFSET_KIND)
       call MPI_FILE_READ_AT(file_handle,offset,wio,1,type_block_io, iostatus,&
          ierrmpi)
       call MPI_SEND(wio,1,type_block_io, ipe,itag,icomm,ierrmpi)
       
       offset=offset + int(size_block_io,kind=MPI_OFFSET_KIND)
       call MPI_FILE_READ_AT(file_handle,offset,wsio,1,&
          type_block_staggered_io, &
                            iostatus_stg,ierrmpi)
       call MPI_SEND(wsio,1,type_block_staggered_io, ipe,itag_stg,icomm,&
          ierrmpi)
       
     end do
    end do
   end if
   call MPI_FILE_CLOSE(file_handle,ierrmpi)
else
   nrecv=(Morton_stop(mype)-Morton_start(mype)+1)
   allocate(iorecvstatus(MPI_STATUS_SIZE,nrecv))

   allocate(iorecvstatus_stg(MPI_STATUS_SIZE,nrecv))
   inrecv=0
   do Morton_no=Morton_start(mype),Morton_stop(mype)
      igrid=sfc_to_igrid(Morton_no)
      itag=Morton_no
      
      itag_stg=Morton_no+ngridshi
      call alloc_node(igrid)
      inrecv=inrecv+1
      call MPI_RECV(pw(igrid)%w,1,type_block_io,0,itag,icomm,iorecvstatus(:,&
         inrecv),ierrmpi)
      
      call MPI_RECV(pws(igrid)%w,1,type_block_staggered_io,0,itag_stg,icomm,&
                    iorecvstatus_stg(:,inrecv),ierrmpi)
      
   end do
   deallocate(iorecvstatus  , iorecvstatus_stg)
end if

call MPI_BARRIER(icomm,ierrmpi)

end subroutine read_snapshot_nopar
!=============================================================================
subroutine read_snapshot_noparf
use mod_forest
use mod_amrvacdef

double precision :: wio(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nw)  ,&
   wsio(ixGlo1-1:ixGhi1,ixGlo2-1:ixGhi2,ixGlo3-1:ixGhi3,1:nws)
integer :: file_handle, amode, igrid, Morton_no, iread
integer :: levmaxini, ndimini, ndirini, nwini, nwsini, neqparini, nxini1,&
   nxini2,nxini3

integer(kind=MPI_ADDRESS_KIND) :: size_double, size_int, lb

integer(kind=MPI_OFFSET_KIND) :: offset
integer, dimension(ngridshi) :: iorequest  , iorequest_stg
integer, dimension(MPI_STATUS_SIZE) :: status
integer, dimension(MPI_STATUS_SIZE) :: iostatus  , iostatus_stg

integer, allocatable :: iorecvstatus(:,:)  , iorecvstatus_stg(:,:)
integer :: ipe,inrecv,nrecv
integer ::  sendini(8+3)
character(len=80) :: filename
logical :: fexist
integer(kind=8) :: file_size
integer :: idx
!-----------------------------------------------------------------------------
call MPI_BARRIER(icomm,ierrmpi)
! generate filename
write(filename,"(a,i4.4,a)") TRIM(filenameini),snapshotini,".dat"
if (mype==0) then
 inquire(file=filename,exist=fexist,size=file_size)
 if(.not.fexist) call mpistop(filename//"as an input snapshot file is not found !")
 ! Open and read file
 size_int=sizeof(4)
 size_double=sizeof(1.0d0)

 offset=file_size+1-(8*size_int+size_double)
 open(unitsnapshot,file=filename,form='unformatted',access='stream',status&
    ='old',action='read')

 read(unitsnapshot,pos=offset) nleafs
 read(unitsnapshot) levmaxini
 read(unitsnapshot) ndimini
 read(unitsnapshot) ndirini
 read(unitsnapshot) nwini
 read(unitsnapshot) nwsini
 read(unitsnapshot) neqparini

 read(unitsnapshot) it

 read(unitsnapshot) t

 ! check if settings are suitable for restart
 if (levmaxini>mxnest) then
    write(*,*) "number of levels in restart file = ",levmaxini
    write(*,*) "mxnest = ",mxnest
    call mpistop("mxnest should be at least number of levels in restart file")
 end if
 if (ndimini/=ndim) then
    write(*,*) "ndim in restart file = ",ndimini
    write(*,*) "ndim = ",ndim
    call mpistop("reset ndim to ndim in restart file")
 end if
 if (ndirini/=ndir) then
    write(*,*) "ndir in restart file = ",ndirini
    write(*,*) "ndir = ",ndir
    call mpistop("reset ndir to ndir in restart file")
 end if
 if (nw/=nwini) then
    write(*,*) "nw=",nw," and nw in restart file=",nwini
    call mpistop("currently, changing nw at restart is not allowed")
 end if

 if (nws/=nwsini) then
   if (mype==0) write(*,*) "nws=",nws," and nws in restart file=",nwsini
   call mpistop("currently, changing nws at restart is not allowed")
 end if


 ! ---
 offset=offset-(ndimini*size_int+neqparini*size_double)
 read(unitsnapshot,pos=offset+(1-1)*size_int)nxini1
 read(unitsnapshot,pos=offset+(2-1)*size_int)nxini2
 read(unitsnapshot,pos=offset+(3-1)*size_int)nxini3
 if (ixGhi1/=nxini1+2*dixB.or.ixGhi2/=nxini2+2*dixB.or.ixGhi3&
    /=nxini3+2*dixB) then
    write(*,*) "Error: reset resolution to ",nxini1+2*dixB,nxini2+2*dixB,&
       nxini3+2*dixB
    call mpistop("change with setamrvac")
 end if
 neqparini=min(neqparini,neqpar+nspecialpar)
 do idx=1,neqparini
     read(unitsnapshot) eqpar(idx)
 end do

 close(unitsnapshot) 
end if ! mype=0

! Broadcast global parameters
if (npe>1) then
  if (mype==0) then
     sendini=(/nleafs,levmaxini,ndimini,ndirini,nwini,  nwsini, neqparini,it ,&
        nxini1,nxini2,nxini3 /)
  end if

  call MPI_BCAST(sendini,8+3,MPI_INTEGER,0,icomm,ierrmpi)
  nleafs=sendini(1);levmaxini=sendini(2);ndimini=sendini(3);
  ndirini=sendini(4);nwini=sendini(5);

  nwsini=sendini(6);neqparini=sendini(7);
  it=sendini(8);
  nxini1=sendini(8+1);nxini2=sendini(8+2);nxini3=sendini(8+3);

  
  nleafs_active = nleafs
  call MPI_BCAST(t,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
  call MPI_BCAST(eqpar,neqparini,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
end if

open(unitsnapshot,file=filename,form='unformatted',access='stream',status&
   ='old',action='read')

! HO: We pass an MPI file handle to avoid changing the interface,
! even though it is not used
call read_forest(file_handle)

! Read and communicate arrays
if (mype==0)then
   iread=0

   do Morton_no=Morton_start(0),Morton_stop(0)
      igrid=sfc_to_igrid(Morton_no)
      call alloc_node(igrid)
      iread=iread+1
      offset=1+(size_block_io  + size_block_stg_io) *(Morton_no-1)

      read(unitsnapshot,pos=offset) pw(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
         ixMlo3:ixMhi3,:)
      
      offset=offset + size_block_io
      read(unitsnapshot,pos=offset) pws(igrid)%w(ixMlo1-1:ixMhi1,&
         ixMlo2-1:ixMhi2,ixMlo3-1:ixMhi3,:) 
   end do
   if (npe>1) then
    do ipe=1,npe-1
     do Morton_no=Morton_start(ipe),Morton_stop(ipe)
       iread=iread+1
       itag=Morton_no
       
       itag_stg=Morton_no+ngridshi
       offset=1+(size_block_io  + size_block_stg_io) *(Morton_no-1)
       read(unitsnapshot,pos=offset) wio(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
          ixMlo3:ixMhi3,:)
       call MPI_SEND(wio,1,type_block_io, ipe,itag,icomm,ierrmpi)
       
       offset=offset + size_block_io
       read(unitsnapshot,pos=offset) wsio(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,&
          ixMlo3-1:ixMhi3,:) 
       call MPI_SEND(wsio,1,type_block_staggered_io,ipe,itag_stg,icomm,&
          ierrmpi)
       
     end do
    end do
   end if

   close(unitsnapshot)

else
   ! mype /= 0
   nrecv=(Morton_stop(mype)-Morton_start(mype)+1)

   allocate(iorecvstatus(MPI_STATUS_SIZE,nrecv))

   allocate(iorecvstatus_stg(MPI_STATUS_SIZE,nrecv))

   iorequest=MPI_REQUEST_NULL
   
   iorequest_stg=MPI_REQUEST_NULL

   inrecv=0
   do Morton_no=Morton_start(mype),Morton_stop(mype)
      igrid=sfc_to_igrid(Morton_no)
      itag=Morton_no
      
      itag_stg=Morton_no+ngridshi
      call alloc_node(igrid)
      inrecv=inrecv+1
      call MPI_RECV(pw(igrid)%w,1,type_block_io,0,itag,icomm,iorecvstatus(:,&
         inrecv),ierrmpi)
      
      call MPI_RECV(pws(igrid)%w,1,type_block_staggered_io,0,itag_stg,icomm,&
                    iorecvstatus_stg(:,inrecv),ierrmpi)
      
   end do
   deallocate(iorecvstatus  , iorecvstatus_stg)
end if

if(nrecv>0)then
  call MPI_WAITALL(nrecv,iorequest,iostatus,ierrmpi) 

  call MPI_WAITALL(nrecv,iorequest_stg,iostatus_stg,ierrmpi)


end if

end subroutine read_snapshot_noparf
!=============================================================================
subroutine printlog_default

! printlog: calculates volume averaged mean values 
use mod_timing
use mod_forest,only:nleafs,nleafs_active,nleafs_level
use mod_amrvacdef

logical          :: fileopen
integer          :: iigrid, igrid, level, iw, i
double precision :: wmean(1:nw), volume(1:nlevelshi), volprob, voltotal
double precision :: dvolume(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
    volumeflat(1:nlevelshi)
integer          :: numlevels, nx1,nx2,nx3, nc, ncells, dit
double precision :: dtTimeLast, now, cellupdatesPerSecond, &
   activeBlocksPerCore, wctPerCodeTime, timeToFinish
integer, dimension(1:nlevelshi) :: isum_send, isum_recv
double precision, dimension(1:nw+1+nlevelshi) :: dsum_send, dsum_recv
character(len=80) :: filename
character(len=2048) :: line
logical, save :: opened=.false.
integer :: amode, status(MPI_STATUS_SIZE)
!-----------------------------------------------------------------------------

volume(1:mxnest)=zero
volumeflat(1:mxnest)=zero
wmean(1:nw)= zero

do iigrid=1,igridstail; igrid=igrids(iigrid);
   level=node(plevel_,igrid)
   volumeflat(level)=volumeflat(level)+ (rnode(rpxmax1_,igrid)-rnode(rpxmin1_,&
      igrid))*(rnode(rpxmax2_,igrid)-rnode(rpxmin2_,igrid))*(rnode(rpxmax3_,&
      igrid)-rnode(rpxmin3_,igrid))
   if (slab) then
      dvolume(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3)=rnode(rpdx1_,&
         igrid)*rnode(rpdx2_,igrid)*rnode(rpdx3_,igrid)
   else
      dvolume(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3)=pgeo(igrid)%dvolume&
         (ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3)
      volume(level)=volume(level)+sum(dvolume(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
         ixMlo3:ixMhi3))
   end if
   do iw=1,nw
      wmean(iw)=wmean(iw)+sum(dvolume(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
         ixMlo3:ixMhi3)*pw(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,&
         iw))
   end do
end do
if (slab) volume(levmin:levmax)=volumeflat(levmin:levmax)

voltotal=sum(volume(levmin:levmax))

numlevels=levmax-levmin+1
dsum_send(1:nw)=wmean(1:nw)
dsum_send(nw+1)=voltotal
dsum_send(nw+2:nw+1+numlevels)=volumeflat(levmin:levmax)
call MPI_REDUCE(dsum_send,dsum_recv,nw+1+numlevels,MPI_DOUBLE_PRECISION,&
    MPI_SUM,0,icomm,ierrmpi)

if (mype==0) then

! To compute cell updates per second, we do the following:
nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;nx3=ixMhi3-ixMlo3+1;
nc=nx1*nx2*nx3
ncells = nc * nleafs_active
! assumes the number of active leafs haven't changed since last compute.
now        = MPI_WTIME()
dit        = it - itTimeLast
dtTimeLast = now - timeLast
itTimeLast = it
timeLast   = now
cellupdatesPerSecond = dble(ncells) * dble(nstep) * dble(dit) &
   / (dtTimeLast * dble(npe))
! blocks per core:
activeBlocksPerCore = dble(nleafs_active) / dble(npe)
! Wall clock time per code time unit in seconds:
if (dt .gt. zero) then 
   wctPerCodeTime = dtTimeLast / (dble(dit) * dt)
else
   wctPerCodeTime = zero
end if

! Wall clock time to finish in hours:
timeToFinish = (tmax - t) * wctPerCodeTime / 3600.0d0

   wmean(1:nw)=dsum_recv(1:nw)
   voltotal=dsum_recv(nw+1)
   volumeflat(levmin:levmax)=dsum_recv(nw+2:nw+1+numlevels)

   wmean=wmean/voltotal

   ! determine coverage in coordinate space
   volprob=(xprobmax1-xprobmin1)*(xprobmax2-xprobmin2)*(xprobmax3-xprobmin3)
   volumeflat(levmin:levmax)=volumeflat(levmin:levmax)/volprob

   if (.not.opened) then
      ! generate filename
      write(filename,"(a,a)") TRIM(filenamelog),".log"

      amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
      amode=ior(amode,MPI_MODE_APPEND)
      call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode, MPI_INFO_NULL,log_fh,&
         ierrmpi)
      opened=.true.

      call MPI_FILE_WRITE(log_fh,fileheadout,len_trim(fileheadout),&
          MPI_CHARACTER,status,ierrmpi)
      !!call MPI_FILE_WRITE(log_fh,new_line('a'),1,MPI_CHARACTER,status,ierrmpi)
      call MPI_FILE_WRITE(log_fh,achar(10),1,MPI_CHARACTER,status,ierrmpi)

      i=len_trim(wnameslog)-1
      do level=1,mxnest
          i=i+3
          if(level<10) then
            if (i+1<1024) write(wnameslog(i:i+1),"(a,i1)") "c",level
          else
            if (i+2<1024) write(wnameslog(i:i+2),"(a,i2)") "c",level
          endif
      end do

      do level=1,mxnest
          i=i+3
          if(level<10) then
            if (i+1<1024) write(wnameslog(i:i+1),"(a,i1)") "n",level
          else
            if (i+2<1024) write(wnameslog(i:i+2),"(a,i2)") "n",level
          endif
      end do
      if (time_accurate) then
         if(residmin>smalldouble) then
           write(line,'(a15,a1024)')"it   t  dt res ",wnameslog
         else
           write(line,'(a15,a1024)')"it   t   dt    ",wnameslog
         endif
      else
         if(residmin>smalldouble) then
           write(line,'(a7,a1024)')"it res ",wnameslog
         else
           write(line,'(a7,a1024)')"it     ",wnameslog
         endif
      end if

      line=trim(line)//"| Xload Xmemory 'Cell_Updates /second/core'"
      line=trim(line)//" 'Active_Blocks/Core' 'Wct Per Code Time [s]' 'TimeToFinish [hrs]'"

      
      call MPI_FILE_WRITE(log_fh,line,len_trim(line),MPI_CHARACTER, status,&
         ierrmpi)
   end if
   !!call MPI_FILE_WRITE(log_fh,new_line('a'),1,MPI_CHARACTER,status,ierrmpi)
   call MPI_FILE_WRITE(log_fh,achar(10),1,MPI_CHARACTER,status,ierrmpi)

   if (time_accurate) then
      if(residmin>smalldouble) then
         write(line,'(i7,3(es12.4))')it,t,dt,residual
      else
         write(line,'(i7,2(es12.4))')it,t,dt
      endif
   else
      if(residmin>smalldouble) then
         write(line,'(i7,1(es12.4))')it,residual
      else
         write(line,'(i7)')it
      endif
   end if
   call MPI_FILE_WRITE(log_fh,line,len_trim(line), MPI_CHARACTER,status,&
      ierrmpi)
   do iw=1,nw
      write(line,'(es12.4)')wmean(iw)
      call MPI_FILE_WRITE(log_fh,line,len_trim(line), MPI_CHARACTER,status,&
         ierrmpi)
   end do
   do level=1,mxnest
      write(line,'(es12.4)')volumeflat(level)
      call MPI_FILE_WRITE(log_fh,line,len_trim(line), MPI_CHARACTER,status,&
         ierrmpi)
   end do

   do level=1,mxnest
      write(line,'(i8)') nleafs_level(level)
      call MPI_FILE_WRITE(log_fh,line,len_trim(line), MPI_CHARACTER,status,&
         ierrmpi)
   end do

   write(line,'(a3,6(es10.2))') ' | ', Xload, Xmemory, cellupdatesPerSecond,&
       activeBlocksPerCore, wctPerCodeTime, timeToFinish
   call MPI_FILE_WRITE(log_fh,line,len_trim(line), MPI_CHARACTER,status,&
      ierrmpi)

end if

end subroutine printlog_default
!=============================================================================

!=============================================================================
