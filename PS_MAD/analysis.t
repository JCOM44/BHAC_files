
!=============================================================================
subroutine write_analysis


  !-----------------------------------------------------------------------------

  call write_analysis_horizon
  call write_analysis_shells

end subroutine write_analysis
!=============================================================================
subroutine write_analysis_shells
  ! Calculates accretion rate, magnetic flux and energy fluxes.  
  ! You can schedule this routine using the slot 5 in itsave, dtsave and ditsave.  
  ! 2016-08-06 by Oliver Porth
  !-----------------------------------------------------------------------------
  use mod_forest, only: Morton_sub_start, Morton_sub_stop
  use mod_slice, only: select_slice, dealloc_subnode
  use mod_metric, only: CoordToBL
  use mod_amrvacdef
  !-----------------------------------------------------------------------------
  integer, parameter                            :: nqttys = 7
  double precision, allocatable, dimension(:,:) :: dsum_send, dsum_recv

  character(len=80) :: filename
  logical           :: file_exists
  integer           :: amode, status(MPI_STATUS_SIZE), Njgrid, jgrid, ix2, ix3
  double precision                               :: roundoff
  double precision, allocatable, dimension(:)    :: normconv 
  
  double precision                               :: xshell
  integer                                        :: mynshells, is
  double precision, allocatable, dimension(:)    :: mdot, phi, edot_em, edot_pake, edot_en, edot, rshell, ldot
  integer                                        :: mdot_, iphi_, em_, pake_, en_, hut_, Trphi_

  integer, save                                  :: filenr = 0
  double precision, dimension(1^D&:1^D&,1:^ND)   :: xBL, x
  !-----------------------------------------------------------------------------

  mynshells = ng1(1)*(ixMhi1-ixMlo1+1)
  allocate(dsum_send(1:nqttys,1:mynshells),dsum_recv(1:nqttys,1:mynshells),normconv(0:nw+nwauxio))
  allocate(mdot(1:mynshells), phi(1:mynshells), edot_em(1:mynshells), edot_pake(1:mynshells), edot_en(1:mynshells), edot(1:mynshells), rshell(1:mynshells), ldot(1:mynshells))
  
  mdot_     = nw+10   ! Make sure indices match in specialvar_output()
  iphi_     = mdot_+1 ! Make sure indices match in specialvar_output()
  em_       = iphi_+1 ! Make sure indices match in specialvar_output()
  pake_     = em_+1   ! Make sure indices match in specialvar_output()
  en_       = pake_+1 ! Make sure indices match in specialvar_output()
  hut_      = en_+1   ! Make sure indices match in specialvar_output()
  Trphi_    = hut_+1  ! Make sure indices match in specialvar_output()

  
  mdot      = 0.0d0
  phi       = 0.0d0
  edot_em   = 0.0d0
  edot_pake = 0.0d0
  edot_en   = 0.0d0
  edot      = 0.0d0
  ldot      = 0.0d0



  !--------------------------------------
  ! Slice through the different shells
  !--------------------------------------
  do is = 1, mynshells
     xshell = xprobmin1 + (dble(is)-half) * dx(1,1) ! Cell center position on base-level resolution
     rshell(is) = xshell
     call select_slice(1,xshell,.false.,unitslice,normconv)
     ! local number of sub-grids:
     Njgrid = Morton_sub_stop(mype) - Morton_sub_start(mype) + 1
     if (Njgrid>0) then 
        do jgrid=1,Njgrid
           {#IFDEF D2
           do ix2=ixMlo2,ixMhi2
              mdot(is)      = mdot(is) + pw_sub(jgrid)%w(ix2,mdot_) * rnode_sub(rpdx1_,jgrid) * 2.0d0 * dpi
              phi(is)       = phi(is) + pw_sub(jgrid)%w(ix2,iphi_) * rnode_sub(rpdx1_,jgrid) * 2.0d0 * dpi
              edot_em(is)   = edot_em(is) + pw_sub(jgrid)%w(ix2,em_) * rnode_sub(rpdx1_,jgrid) * 2.0d0 * dpi
              edot_pake(is) = edot_pake(is) + pw_sub(jgrid)%w(ix2,pake_) * rnode_sub(rpdx1_,jgrid) * 2.0d0 * dpi
              edot_en(is)   = edot_en(is) + pw_sub(jgrid)%w(ix2,en_) * rnode_sub(rpdx1_,jgrid) * 2.0d0 * dpi
              edot(is)      = edot_em(is) + edot_pake(is) + edot_en(is)
              ldot(is)      = ldot(is) + pw_sub(jgrid)%w(ix2,Trphi_) * rnode_sub(rpdx1_,jgrid) * 2.0d0 * dpi
           end do
           }{#IFDEF D3
           do ix2=ixMlo2,ixMhi2
              do ix3=ixMlo3,ixMhi3
                 mdot(is)      = mdot(is) + pw_sub(jgrid)%w(ix2,ix3,mdot_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid)
                 phi(is)       = phi(is) + pw_sub(jgrid)%w(ix2,ix3,iphi_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid)
                 edot_em(is)   = edot_em(is) + pw_sub(jgrid)%w(ix2,ix3,em_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid)
                 edot_pake(is) = edot_pake(is) + pw_sub(jgrid)%w(ix2,ix3,pake_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid)
                 edot_en(is)   = edot_en(is) + pw_sub(jgrid)%w(ix2,ix3,en_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid)
                 edot(is)      = edot_em(is) + edot_pake(is) + edot_en(is)
                 ldot(is)      = ldot(is) + pw_sub(jgrid)%w(ix2,ix3,Trphi_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid)
              end do
           end do
           }
        end do
        do jgrid=1,Njgrid
           call dealloc_subnode(jgrid)
        end do
     end if ! Njgrid>0
  end do
  !------------------------------
  ! DONE: Slice through the different shells
  !------------------------------


 
  dsum_send(1,:) = mdot(:)
  dsum_send(2,:) = phi(:)
  dsum_send(3,:) = edot_em(:)
  dsum_send(4,:) = edot_pake(:)
  dsum_send(5,:) = edot_en(:)
  dsum_send(6,:) = edot(:)
  dsum_send(7,:) = ldot(:)

  call MPI_REDUCE(dsum_send,dsum_recv,nqttys*mynshells,MPI_DOUBLE_PRECISION, &
       MPI_SUM,0,icomm,ierrmpi)


  if (mype .eq. 0) then

     ! generate filename
     write(filename,"(a,a,i4.4,a)") trim(filenameout),TRIM("_shells"),snapshot,".csv"
     filenr = filenr + 1
     INQUIRE(FILE=filename, EXIST=file_exists)
     if (.not. file_exists) then
        open(unit=unitanalysis,file=filename,status='unknown',position='append')
        write(unitanalysis,"(a,ES14.6,a)",advance='yes') '# t=',t
        write(unitanalysis,"(a)",advance='yes') '# r mdot phi edot_em edot_pake edot_en edot ldot'
     else
        open(unit=unitanalysis,file=filename,status='unknown',position='append')
     end if

  
     ! use roundoff because ascii can give underflows which makes the data unreadable
     do is = 1, mynshells
        x(1^D&,1) = rshell(is)
        {^IFZIN   x(1^D&,^Z)   = dpi/2.0d0 }
        {^IFPHIIN x(1^D&,^PHI) = zero }
        call CoordToBL(1^D&,1^D&,1^D&,1^D&,x,xBL)
        
        rshell(is)    = roundoff(xBL,1.0d-99)
        mdot(is)      = roundoff(dsum_recv(1,is),1.0d-99)
        phi(is)       = roundoff(dsum_recv(2,is),1.0d-99)
        edot_em(is)   = roundoff(dsum_recv(3,is),1.0d-99)
        edot_pake(is) = roundoff(dsum_recv(4,is),1.0d-99)
        edot_en(is)   = roundoff(dsum_recv(5,is),1.0d-99)
        edot(is)      = roundoff(dsum_recv(6,is),1.0d-99)
        ldot(is)      = roundoff(dsum_recv(7,is),1.0d-99)
        write(unitanalysis,'(8(ES14.6))', advance='yes') rshell(is), &
             mdot(is), phi(is), edot_em(is), edot_pake(is), &
             edot_en(is), edot(is), ldot(is)
     end do
     close(unitanalysis)

  end if

  deallocate(dsum_send,dsum_recv,normconv)
  deallocate(mdot, phi, edot_em, edot_pake, edot_en, edot, rshell, ldot)

end subroutine write_analysis_shells
!=============================================================================
subroutine write_analysis_horizon
  ! Calculates accretion rate, magnetic flux and energy fluxes.  
  ! You can schedule this routine using the slot 5 in itsave, dtsave and ditsave.  
  ! 2016-08-06 by Oliver Porth
  !-----------------------------------------------------------------------------
  use mod_forest, only: Morton_sub_start, Morton_sub_stop
!  use mod_metric, only: outerhorizon, CoordToBL
  use mod_metric, only: CoordToBL
  use mod_slice, only: select_slice, dealloc_subnode
  use mod_amrvacdef
  !-----------------------------------------------------------------------------
  integer :: iigrid, igrid, level
  integer, parameter                    :: nqttys = 9, naux=26
  double precision, dimension(1:nqttys) :: dsum_send, dsum_recv

  character(len=80) :: filename
  logical       :: file_exists
  integer       :: amode, status(MPI_STATUS_SIZE), Njgrid, jgrid, ix2, ix3
  logical                                        :: sliceascii_save
  integer                                        :: nwauxio_save
  double precision                               :: roundoff
  
  double precision, dimension(ixG^T,1:nw+naux)   :: w
  double precision,dimension(0:nw+naux)          :: normconv 

  double precision, dimension(ixG^T,1:ndim)      :: xBL
  double precision, dimension(ixG^T)             :: mask
  double precision                               :: xouter
  double precision                               :: mdot, phi, edot_em, edot_pake, edot_en, edot, j, Ldot, rbar1, rbar2
  double precision                               :: tmp, tmp2
  integer                                        :: mdot_, iphi_, em_, pake_, en_, hut_, Trphi_, j_, rbary_, rhosqrtg_
  !-----------------------------------------------------------------------------

  !--------------------------------------
  ! Hack the slicing tool:
  !--------------------------------------
  sliceascii_save = sliceascii
  nwauxio_save    = nwauxio
  sliceascii      = .true.
  nwauxio         = naux   ! Make sure indices match in specialvar_output()
  !--------------------------------------
  ! Make sure to reset globals at the end
  !--------------------------------------

  mdot_     = nw+10   ! Make sure indices match in specialvar_output()
  iphi_     = mdot_+1 ! Make sure indices match in specialvar_output()
  em_       = iphi_+1 ! Make sure indices match in specialvar_output()
  pake_     = em_+1   ! Make sure indices match in specialvar_output()
  en_       = pake_+1 ! Make sure indices match in specialvar_output()
  hut_      = en_+1   ! Make sure indices match in specialvar_output()
  Trphi_    = hut_+1  ! Make sure indices match in specialvar_output()
  j_        = Trphi_+1! Make sure indices match in specialvar_output()
  rbary_    = nw+24   ! Make sure indices match in specialvar_output()
  rhosqrtg_ = nw+25   ! Make sure indices match in specialvar_output()

  
  mdot      = 0.0d0
  phi       = 0.0d0
  edot_em   = 0.0d0
  edot_pake = 0.0d0
  edot_en   = 0.0d0
  edot      = 0.0d0
  j         = 0.0d0
  Ldot      = 0.0d0
  rbar1     = 0.0d0
  rbar2     = 0.0d0
  tmp       = 0.0d0
  tmp2      = 0.0d0



  !--------------------------------------
  ! Slice through the horizon:
  !--------------------------------------
  xouter = log(6.0d0+1.0d0) !outerhorizon()
  call select_slice(1,xouter,.false.,unitslice,normconv)
  ! local number of sub-grids:
  Njgrid = Morton_sub_stop(mype) - Morton_sub_start(mype) + 1
  if (Njgrid>0) then 
     do jgrid=1,Njgrid
        {#IFDEF D2
        do ix2=ixMlo2,ixMhi2
           mdot      = mdot + pw_sub(jgrid)%w(ix2,mdot_) * rnode_sub(rpdx1_,jgrid) * 2.0d0 * dpi
           phi       = phi + pw_sub(jgrid)%w(ix2,iphi_) * rnode_sub(rpdx1_,jgrid) * 2.0d0 * dpi
           edot_em   = edot_em + pw_sub(jgrid)%w(ix2,em_) * rnode_sub(rpdx1_,jgrid) * 2.0d0 * dpi
           edot_pake = edot_pake + pw_sub(jgrid)%w(ix2,pake_) * rnode_sub(rpdx1_,jgrid) * 2.0d0 * dpi
           edot_en   = edot_en + pw_sub(jgrid)%w(ix2,en_) * rnode_sub(rpdx1_,jgrid) * 2.0d0 * dpi
           edot      = edot_em + edot_pake + edot_en
           Ldot      = Ldot + pw_sub(jgrid)%w(ix2,Trphi_) * rnode_sub(rpdx1_,jgrid) * 2.0d0 * dpi
        end do
        }{#IFDEF D3
        do ix2=ixMlo2,ixMhi2
           do ix3=ixMlo3,ixMhi3
              mdot      = mdot + pw_sub(jgrid)%w(ix2,ix3,mdot_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid)
              phi       = phi + pw_sub(jgrid)%w(ix2,ix3,iphi_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid)
              edot_em   = edot_em + pw_sub(jgrid)%w(ix2,ix3,em_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid)
              edot_pake = edot_pake + pw_sub(jgrid)%w(ix2,ix3,pake_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid)
              edot_en   = edot_en + pw_sub(jgrid)%w(ix2,ix3,en_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid)
              edot      = edot_em + edot_pake + edot_en
              Ldot      = Ldot + pw_sub(jgrid)%w(ix2,ix3,Trphi_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid)
           end do
        end do
        }
     end do
     do jgrid=1,Njgrid
        call dealloc_subnode(jgrid)
     end do
  end if ! Njgrid>0
  !------------------------------
  ! DONE: Slice through the horizon
  !------------------------------


  !------------------------------
  ! Volume integrate
  !------------------------------
  do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
     call set_tmpGlobals(igrid)
     call CoordToBL(ixG^LL,ixM^LL,px(igrid)%x,xBL)
     
     w(ixG^T,1:nw)=pw(igrid)%w(ixG^T,1:nw)
     call specialvar_output(ixG^LL,ixM^LL,rhosqrtg_,w,ps(igrid),normconv)

     j = j + sum(w(ixM^T,j_)*mygeo%dVolume(ixM^T))

     ! For the full domain:
     rbar1 = rbar1 + sum(w(ixM^T,rbary_)*mygeo%dVolume(ixM^T))
     
     ! Only inside of domian of interest:
     where (xBL(ixM^T,r_) .le. 50.0d0)
        mask(ixM^T) = 1.0d0
     elsewhere
        mask(ixM^T) = 0.0d0
     end where
     rbar2 = rbar2 + sum(w(ixM^T,rbary_)*mygeo%dVolume(ixM^T)*mask(ixM^T))

     ! Normalization:
     tmp  = tmp + sum(w(ixM^T,rhosqrtg_)*mygeo%dVolume(ixM^T))
     tmp2 = tmp2 + sum(w(ixM^T,rhosqrtg_)*mygeo%dVolume(ixM^T)*mask(ixM^T))
     
  end do
  !------------------------------
  ! Done: Volume integrate
  !------------------------------
  
  
  dsum_send(1) = mdot
  dsum_send(2) = phi
  dsum_send(3) = edot
  dsum_send(4) = Ldot
  dsum_send(5) = j
  dsum_send(6) = rbar1
  dsum_send(7) = rbar2
  dsum_send(8) = tmp
  dsum_send(9) = tmp2

  call MPI_REDUCE(dsum_send,dsum_recv,nqttys,MPI_DOUBLE_PRECISION, &
       MPI_SUM,0,icomm,ierrmpi)


  if (mype .eq. 0) then

     ! use roundoff because ascii can give underflows which makes the data unreadable
     mdot      = roundoff(dsum_recv(1),1.0d-99)
     phi       = roundoff(dsum_recv(2),1.0d-99)
     edot      = roundoff(dsum_recv(3),1.0d-99)
     Ldot      = roundoff(dsum_recv(4),1.0d-99)
     j         = roundoff(dsum_recv(5),1.0d-99)
     rbar1     = roundoff(dsum_recv(6)/dsum_recv(8),1.0d-99)
     rbar2     = roundoff(dsum_recv(7)/dsum_recv(9),1.0d-99)

     ! generate filename
     write(filename,"(a,a,a)") trim(filenameout),TRIM("_codeComparison"),".csv"
     INQUIRE(FILE=filename, EXIST=file_exists)
     if (.not. file_exists) then
        open(unit=unitanalysis,file=filename,status='unknown',position='append')
        write(unitanalysis,"(a,a)",advance='yes') trim('# t mdot phi edot Ldot j rbar rbar_doi')
     else
        open(unit=unitanalysis,file=filename,status='unknown',position='append')
     end if
     write(unitanalysis,'(8(ES14.6))', advance='yes') t, mdot, phi, edot, Ldot, j, rbar1, rbar2
     close(unitanalysis)

  end if


  !--------------------------------------
  ! Hack the slicing tool:
  !--------------------------------------
  sliceascii = sliceascii_save
  nwauxio    = nwauxio_save
  !--------------------------------------
  ! Done resetting globals
  !--------------------------------------

end subroutine write_analysis_horizon
!=============================================================================
