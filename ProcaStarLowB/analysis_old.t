!=============================================================================
subroutine write_analysis
  ! Calculates accretion rate, magnetic flux and energy fluxes.  
  ! You can schedule this routine using the slot 5 in itsave, dtsave and ditsave.  
  ! 2016-08-06 by Oliver Porth
  !-----------------------------------------------------------------------------
  use mod_forest, only: Morton_sub_start, Morton_sub_stop
  use mod_metric, only:  BLToCoord !, outerhorizon
  use mod_slice, only: select_slice, dealloc_subnode
  use mod_amrvacdef
  !-----------------------------------------------------------------------------
  integer :: iigrid, igrid, level
  integer, parameter                    :: nqttys = 35, naux=15
  double precision, dimension(1:nqttys) :: dsum_send, dsum_recv

  character(len=80) :: filename
  logical       :: file_exists
  integer       :: amode, status(MPI_STATUS_SIZE), Njgrid, jgrid, ix2, ix3
  double precision,dimension(0:nw+naux)          :: normconv    ! Make sure indices match in specialvar_output()
  logical                                        :: sliceascii_save, saveprim_save
  integer                                        :: nwauxio_save
  double precision                               :: roundoff

  double precision                               :: xouter
  double precision, dimension(1^D&:1^D&,1:^ND)   :: xBL50, x50
  double precision                               :: mdot, phi, edot_em, edot_pake, edot_en
  double precision                               :: mdot50, phi50, edot_em50, edot_pake50, edot_en50
  double precision                               :: mdot50_j, phi50_j, edot_em50_j, edot_pake50_j, edot_en50_j
  double precision                               :: mdot50_jn, phi50_jn, edot_em50_jn, edot_pake50_jn, edot_en50_jn
  double precision                               :: mdot50_js, phi50_js, edot_em50_js, edot_pake50_js, edot_en50_js
  double precision                               :: mdot50_w, phi50_w, edot_em50_w, edot_pake50_w, edot_en50_w
  double precision                               :: mdot50_disk, phi50_disk, edot_em50_disk, edot_pake50_disk, edot_en50_disk
  {#IFDEF D2
  double precision, dimension(ixGlo2:ixGhi2)     :: inJet, inWind, inDisk, mu
  }{#IFDEF D3
  double precision, dimension(ixGlo2:ixGhi2,ixGlo3:ixGhi3)     :: inJet, inWind, inDisk, mu
  }
  integer                                        :: ib2_, mdot_, iphi_, em_, pake_, en_, hut_
  !--------------------------------------
  ! Thresholds:
  double precision, parameter                    :: muthresh  = 2.0d0 ! Jet material should reach Gamma>muthresh
  double precision, parameter                    :: hutthresh = 1.0d0 ! Unbound wind material has -h u_1 > hutthresh
  !--------------------------------------
  !-----------------------------------------------------------------------------

  !--------------------------------------
  ! Hack the slicing tool:
  !--------------------------------------
  sliceascii_save = sliceascii
  saveprim_save   = saveprim
  nwauxio_save    = nwauxio
  sliceascii      = .true.
  saveprim        = .true.
  nwauxio         = naux   ! Make sure indices match in specialvar_output()
  !--------------------------------------
  ! Make sure to reset globals at the end
  !--------------------------------------

  ib2_      = nw+1    ! Make sure indices match in specialvar_output()
  mdot_     = nw+10   ! Make sure indices match in specialvar_output()
  iphi_     = mdot_+1 ! Make sure indices match in specialvar_output()
  em_       = iphi_+1 ! Make sure indices match in specialvar_output()
  pake_     = em_+1   ! Make sure indices match in specialvar_output()
  en_       = pake_+1 ! Make sure indices match in specialvar_output()
  hut_      = en_+1   ! Make sure indices match in specialvar_output()



  !--------------------------------------
  ! Slice through the horizon:
  !--------------------------------------

  mdot      = 0.0d0
  phi       = 0.0d0
  edot_em   = 0.0d0
  edot_pake = 0.0d0
  edot_en   = 0.0d0

  xouter = 2.0d0 !outerhorizon()

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
        end do
        }{#IFDEF D3
        do ix2=ixMlo2,ixMhi2
           do ix3=ixMlo3,ixMhi3
              mdot      = mdot + pw_sub(jgrid)%w(ix2,ix3,mdot_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid)
              phi       = phi + pw_sub(jgrid)%w(ix2,ix3,iphi_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid)
              edot_em   = edot_em + pw_sub(jgrid)%w(ix2,ix3,em_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid)
              edot_pake = edot_pake + pw_sub(jgrid)%w(ix2,ix3,pake_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid)
              edot_en   = edot_en + pw_sub(jgrid)%w(ix2,ix3,en_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid)
           end do
        end do
        }
     end do
     do jgrid=1,Njgrid
        call dealloc_subnode(jgrid)
     end do
  end if ! Njgrid>0

  !--------------------------------------
  ! Slice through rBL=50:
  !--------------------------------------

  ! total qttys:
  mdot50      = 0.0d0
  phi50       = 0.0d0
  edot_em50   = 0.0d0
  edot_pake50 = 0.0d0
  edot_en50   = 0.0d0

  ! jet, northern and southern qttys:
  mdot50_j      = 0.0d0;   mdot50_jn      = 0.0d0;   mdot50_js      = 0.0d0
  phi50_j       = 0.0d0;   phi50_jn       = 0.0d0;   phi50_js       = 0.0d0
  edot_em50_j   = 0.0d0;   edot_em50_jn   = 0.0d0;   edot_em50_js   = 0.0d0
  edot_pake50_j = 0.0d0;   edot_pake50_jn = 0.0d0;   edot_pake50_js = 0.0d0
  edot_en50_j   = 0.0d0;   edot_en50_jn   = 0.0d0;   edot_en50_js   = 0.0d0

  ! wind qttys:
  mdot50_w      = 0.0d0
  phi50_w       = 0.0d0
  edot_em50_w   = 0.0d0
  edot_pake50_w = 0.0d0
  edot_en50_w   = 0.0d0

  ! disk qttys:
  mdot50_disk      = 0.0d0
  phi50_disk       = 0.0d0
  edot_em50_disk   = 0.0d0
  edot_pake50_disk = 0.0d0
  edot_en50_disk   = 0.0d0

  xBL50(1^D&,1) = 50.0d0
  {^IFZIN   xBL50(1^D&,^Z)   = dpi/2.0d0 }
  {^IFPHIIN xBL50(1^D&,^PHI) = zero }

!  call BLToCoord(1^D&,1^D&,1^D&,1^D&,xBL50,x50)
! Dont have BLToCoord yet for the CMKS metric, approximately it is:
  x50 = xBL50
  x50(1^D&,1) = log(xBL50(1^D&,1))

  call select_slice(1,x50(1^D&,1),.false.,unitslice,normconv)
  ! local number of sub-grids:
  Njgrid = Morton_sub_stop(mype) - Morton_sub_start(mype) + 1
  if (Njgrid>0) then 
     do jgrid=1,Njgrid

        {#IFDEF D2
        !Associate for better readability:
        associate (w=>pw_sub(jgrid)%w)
        mu(ixMlo2:ixMhi2) = - (w(ixMlo2:ixMhi2,em_)+w(ixMlo2:ixMhi2,pake_)+w(ixMlo2:ixMhi2,en_)-w(ixMlo2:ixMhi2,mdot_)) &
             / w(ixMlo2:ixMhi2,mdot_)
        end associate


        ! Set jet mask:
        where (pw_sub(jgrid)%w(ixMlo2:ixMhi2,ib2_)/pw_sub(jgrid)%w(ixMlo2:ixMhi2,rho_) .ge. 1.0d0 .or. &
           mu(ixMlo2:ixMhi2) .ge. muthresh)
           inJet(ixMlo2:ixMhi2) = 1.0d0
        elsewhere
           inJet(ixMlo2:ixMhi2) = 0.0d0
        end where

        ! Set wind mask:
        where (inJet(ixMlo2:ixMhi2) .ne. 1.0d0 .and. pw_sub(jgrid)%w(ixMlo2:ixMhi2,hut_) .gt. hutthresh)
           inWind(ixMlo2:ixMhi2) = 1.0d0
        elsewhere
           inWind(ixMlo2:ixMhi2) = 0.0d0
        end where

        ! Set disk mask:
        where (inJet(ixMlo2:ixMhi2) .ne. 1.0d0 .and. inWind(ixMlo2:ixMhi2) .ne. 1.0d0)
           inDisk(ixMlo2:ixMhi2) = 1.0d0
        elsewhere
           inDisk(ixMlo2:ixMhi2) = 0.0d0
        end where


        do ix2=ixMlo2,ixMhi2
           mdot50      = mdot50 + pw_sub(jgrid)%w(ix2,mdot_) * rnode_sub(rpdx1_,jgrid) * 2.0d0 * dpi
           phi50       = phi50 + pw_sub(jgrid)%w(ix2,iphi_) * rnode_sub(rpdx1_,jgrid) * 2.0d0 * dpi
           edot_em50   = edot_em50 + pw_sub(jgrid)%w(ix2,em_) * rnode_sub(rpdx1_,jgrid) * 2.0d0 * dpi
           edot_pake50 = edot_pake50 + pw_sub(jgrid)%w(ix2,pake_) * rnode_sub(rpdx1_,jgrid) * 2.0d0 * dpi
           edot_en50   = edot_en50 + pw_sub(jgrid)%w(ix2,en_) * rnode_sub(rpdx1_,jgrid) * 2.0d0 * dpi

           mdot50_j      = mdot50_j + pw_sub(jgrid)%w(ix2,mdot_) * rnode_sub(rpdx1_,jgrid) * 2.0d0*dpi * inJet(ix2)
           phi50_j       = phi50_j + pw_sub(jgrid)%w(ix2,iphi_) * rnode_sub(rpdx1_,jgrid) * 2.0d0*dpi * inJet(ix2)
           edot_em50_j   = edot_em50_j + pw_sub(jgrid)%w(ix2,em_) * rnode_sub(rpdx1_,jgrid) * 2.0d0*dpi * inJet(ix2)
           edot_pake50_j = edot_pake50_j + pw_sub(jgrid)%w(ix2,pake_) * rnode_sub(rpdx1_,jgrid) * 2.0d0*dpi * inJet(ix2)
           edot_en50_j   = edot_en50_j + pw_sub(jgrid)%w(ix2,en_) * rnode_sub(rpdx1_,jgrid) * 2.0d0*dpi * inJet(ix2)

           if (px_sub(jgrid)%x(ix2,^Z) .lt. dpi/2.0d0) then 
              ! In the northern hemisphere:
              mdot50_jn      = mdot50_jn + pw_sub(jgrid)%w(ix2,mdot_) * rnode_sub(rpdx1_,jgrid) * 2.0d0*dpi * inJet(ix2)
              phi50_jn       = phi50_jn + pw_sub(jgrid)%w(ix2,iphi_) * rnode_sub(rpdx1_,jgrid) * 2.0d0*dpi * inJet(ix2)
              edot_em50_jn   = edot_em50_jn + pw_sub(jgrid)%w(ix2,em_) * rnode_sub(rpdx1_,jgrid) * 2.0d0*dpi * inJet(ix2)
              edot_pake50_jn = edot_pake50_jn + pw_sub(jgrid)%w(ix2,pake_) * rnode_sub(rpdx1_,jgrid) * 2.0d0*dpi * inJet(ix2)
              edot_en50_jn   = edot_en50_jn + pw_sub(jgrid)%w(ix2,en_) * rnode_sub(rpdx1_,jgrid) * 2.0d0*dpi * inJet(ix2)
           else
              ! In the southern hemisphere:
              mdot50_js      = mdot50_js + pw_sub(jgrid)%w(ix2,mdot_) * rnode_sub(rpdx1_,jgrid) * 2.0d0*dpi * inJet(ix2)
              phi50_js       = phi50_js + pw_sub(jgrid)%w(ix2,iphi_) * rnode_sub(rpdx1_,jgrid) * 2.0d0*dpi * inJet(ix2)
              edot_em50_js   = edot_em50_js + pw_sub(jgrid)%w(ix2,em_) * rnode_sub(rpdx1_,jgrid) * 2.0d0*dpi * inJet(ix2)
              edot_pake50_js = edot_pake50_js + pw_sub(jgrid)%w(ix2,pake_) * rnode_sub(rpdx1_,jgrid) * 2.0d0*dpi * inJet(ix2)
              edot_en50_js   = edot_en50_js + pw_sub(jgrid)%w(ix2,en_) * rnode_sub(rpdx1_,jgrid) * 2.0d0*dpi * inJet(ix2)
           end if
              
           mdot50_w      = mdot50_w + pw_sub(jgrid)%w(ix2,mdot_) * rnode_sub(rpdx1_,jgrid) * 2.0d0*dpi * inWind(ix2)
           phi50_w       = phi50_w + pw_sub(jgrid)%w(ix2,iphi_) * rnode_sub(rpdx1_,jgrid) * 2.0d0*dpi * inWind(ix2)
           edot_em50_w   = edot_em50_w + pw_sub(jgrid)%w(ix2,em_) * rnode_sub(rpdx1_,jgrid) * 2.0d0*dpi * inWind(ix2)
           edot_pake50_w = edot_pake50_w + pw_sub(jgrid)%w(ix2,pake_) * rnode_sub(rpdx1_,jgrid) * 2.0d0*dpi * inWind(ix2)
           edot_en50_w   = edot_en50_w + pw_sub(jgrid)%w(ix2,en_) * rnode_sub(rpdx1_,jgrid) * 2.0d0*dpi * inWind(ix2)

           mdot50_disk      = mdot50_disk + pw_sub(jgrid)%w(ix2,mdot_) * rnode_sub(rpdx1_,jgrid) * 2.0d0*dpi * inDisk(ix2)
           phi50_disk       = phi50_disk + pw_sub(jgrid)%w(ix2,iphi_) * rnode_sub(rpdx1_,jgrid) * 2.0d0*dpi * inDisk(ix2)
           edot_em50_disk   = edot_em50_disk + pw_sub(jgrid)%w(ix2,em_) * rnode_sub(rpdx1_,jgrid) * 2.0d0*dpi * inDisk(ix2)
           edot_pake50_disk = edot_pake50_disk + pw_sub(jgrid)%w(ix2,pake_) * rnode_sub(rpdx1_,jgrid) * 2.0d0*dpi * inDisk(ix2)
           edot_en50_disk   = edot_en50_disk + pw_sub(jgrid)%w(ix2,en_) * rnode_sub(rpdx1_,jgrid) * 2.0d0*dpi * inDisk(ix2)
        end do
        }


        {#IFDEF D3
        !Associate for better readability:
        associate (w=>pw_sub(jgrid)%w)
        mu(ixMlo2:ixMhi2,ixMlo3:ixMhi3) = - (w(ixMlo2:ixMhi2,ixMlo3:ixMhi3,em_)+w(ixMlo2:ixMhi2,ixMlo3:ixMhi3,pake_)+w(ixMlo2:ixMhi2,ixMlo3:ixMhi3,en_)-w(ixMlo2:ixMhi2,ixMlo3:ixMhi3,mdot_)) &
             / w(ixMlo2:ixMhi2,ixMlo3:ixMhi3,mdot_)
        end associate


        ! Set jet mask:
        where (pw_sub(jgrid)%w(ixMlo2:ixMhi2,ixMlo3:ixMhi3,ib2_)/pw_sub(jgrid)%w(ixMlo2:ixMhi2,ixMlo3:ixMhi3,rho_) .ge. 1.0d0 .or. &
           mu(ixMlo2:ixMhi2,ixMlo3:ixMhi3) .ge. muthresh)
           inJet(ixMlo2:ixMhi2,ixMlo3:ixMhi3) = 1.0d0
        elsewhere
           inJet(ixMlo2:ixMhi2,ixMlo3:ixMhi3) = 0.0d0
        end where

        ! Set wind mask:
        where (inJet(ixMlo2:ixMhi2,ixMlo3:ixMhi3) .ne. 1.0d0 .and. pw_sub(jgrid)%w(ixMlo2:ixMhi2,ixMlo3:ixMhi3,hut_) .gt. hutthresh)
           inWind(ixMlo2:ixMhi2,ixMlo3:ixMhi3) = 1.0d0
        elsewhere
           inWind(ixMlo2:ixMhi2,ixMlo3:ixMhi3) = 0.0d0
        end where

        ! Set disk mask:
        where (inJet(ixMlo2:ixMhi2,ixMlo3:ixMhi3) .ne. 1.0d0 .and. inWind(ixMlo2:ixMhi2,ixMlo3:ixMhi3) .ne. 1.0d0)
           inDisk(ixMlo2:ixMhi2,ixMlo3:ixMhi3) = 1.0d0
        elsewhere
           inDisk(ixMlo2:ixMhi2,ixMlo3:ixMhi3) = 0.0d0
        end where


        do ix3=ixMlo3,ixMhi3
           do ix2=ixMlo2,ixMhi2
              mdot50      = mdot50 + pw_sub(jgrid)%w(ix2,ix3,mdot_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid)
              phi50       = phi50 + pw_sub(jgrid)%w(ix2,ix3,iphi_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid)
              edot_em50   = edot_em50 + pw_sub(jgrid)%w(ix2,ix3,em_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid)
              edot_pake50 = edot_pake50 + pw_sub(jgrid)%w(ix2,ix3,pake_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid)
              edot_en50   = edot_en50 + pw_sub(jgrid)%w(ix2,ix3,en_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid)

              mdot50_j      = mdot50_j + pw_sub(jgrid)%w(ix2,ix3,mdot_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid) * inJet(ix2,ix3)
              phi50_j       = phi50_j + pw_sub(jgrid)%w(ix2,ix3,iphi_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid) * inJet(ix2,ix3)
              edot_em50_j   = edot_em50_j + pw_sub(jgrid)%w(ix2,ix3,em_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid) * inJet(ix2,ix3)
              edot_pake50_j = edot_pake50_j + pw_sub(jgrid)%w(ix2,ix3,pake_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid) * inJet(ix2,ix3)
              edot_en50_j   = edot_en50_j + pw_sub(jgrid)%w(ix2,ix3,en_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid) * inJet(ix2,ix3)

              if (px_sub(jgrid)%x(ix2,ix3,^Z) .lt. dpi/2.0d0) then
                 ! In the northern hemisphere:
                 mdot50_jn      = mdot50_jn + pw_sub(jgrid)%w(ix2,ix3,mdot_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid) * inJet(ix2,ix3)
                 phi50_jn       = phi50_jn + pw_sub(jgrid)%w(ix2,ix3,iphi_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid) * inJet(ix2,ix3)
                 edot_em50_jn   = edot_em50_jn + pw_sub(jgrid)%w(ix2,ix3,em_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid) * inJet(ix2,ix3)
                 edot_pake50_jn = edot_pake50_jn + pw_sub(jgrid)%w(ix2,ix3,pake_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid) * inJet(ix2,ix3)
                 edot_en50_jn   = edot_en50_jn + pw_sub(jgrid)%w(ix2,ix3,en_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid) * inJet(ix2,ix3)
              else
                 ! In the southern hemisphere:
                 mdot50_js      = mdot50_js + pw_sub(jgrid)%w(ix2,ix3,mdot_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid) * inJet(ix2,ix3)
                 phi50_js       = phi50_js + pw_sub(jgrid)%w(ix2,ix3,iphi_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid) * inJet(ix2,ix3)
                 edot_em50_js   = edot_em50_js + pw_sub(jgrid)%w(ix2,ix3,em_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid) * inJet(ix2,ix3)
                 edot_pake50_js = edot_pake50_js + pw_sub(jgrid)%w(ix2,ix3,pake_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid) * inJet(ix2,ix3)
                 edot_en50_js   = edot_en50_js + pw_sub(jgrid)%w(ix2,ix3,en_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid) * inJet(ix2,ix3)
              end if


              mdot50_w      = mdot50_w + pw_sub(jgrid)%w(ix2,ix3,mdot_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid) * inWind(ix2,ix3)
              phi50_w       = phi50_w + pw_sub(jgrid)%w(ix2,ix3,iphi_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid) * inWind(ix2,ix3)
              edot_em50_w   = edot_em50_w + pw_sub(jgrid)%w(ix2,ix3,em_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid) * inWind(ix2,ix3)
              edot_pake50_w = edot_pake50_w + pw_sub(jgrid)%w(ix2,ix3,pake_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid) * inWind(ix2,ix3)
              edot_en50_w   = edot_en50_w + pw_sub(jgrid)%w(ix2,ix3,en_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid) * inWind(ix2,ix3)

              mdot50_disk      = mdot50_disk + pw_sub(jgrid)%w(ix2,ix3,mdot_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid) * inDisk(ix2,ix3)
              phi50_disk       = phi50_disk + pw_sub(jgrid)%w(ix2,ix3,iphi_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid) * inDisk(ix2,ix3)
              edot_em50_disk   = edot_em50_disk + pw_sub(jgrid)%w(ix2,ix3,em_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid) * inDisk(ix2,ix3)
              edot_pake50_disk = edot_pake50_disk + pw_sub(jgrid)%w(ix2,ix3,pake_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid) * inDisk(ix2,ix3)
              edot_en50_disk   = edot_en50_disk + pw_sub(jgrid)%w(ix2,ix3,en_) * rnode_sub(rpdx1_,jgrid) * rnode_sub(rpdx2_,jgrid) * inDisk(ix2,ix3)
           end do
        end do
        }
     end do
     do jgrid=1,Njgrid
        call dealloc_subnode(jgrid)
     end do
  end if ! Njgrid>0


  dsum_send(1) = mdot
  dsum_send(2) = phi
  dsum_send(3) = edot_em
  dsum_send(4) = edot_pake
  dsum_send(5) = edot_en

  dsum_send(6) = mdot50
  dsum_send(7) = phi50
  dsum_send(8) = edot_em50
  dsum_send(9) = edot_pake50
  dsum_send(10) = edot_en50

  dsum_send(11) = mdot50_j
  dsum_send(12) = phi50_j
  dsum_send(13) = edot_em50_j
  dsum_send(14) = edot_pake50_j
  dsum_send(15) = edot_en50_j

  dsum_send(16) = mdot50_w
  dsum_send(17) = phi50_w
  dsum_send(18) = edot_em50_w
  dsum_send(19) = edot_pake50_w
  dsum_send(20) = edot_en50_w

  dsum_send(21) = mdot50_disk
  dsum_send(22) = phi50_disk
  dsum_send(23) = edot_em50_disk
  dsum_send(24) = edot_pake50_disk
  dsum_send(25) = edot_en50_disk

  dsum_send(26) = mdot50_jn
  dsum_send(27) = phi50_jn
  dsum_send(28) = edot_em50_jn
  dsum_send(29) = edot_pake50_jn
  dsum_send(30) = edot_en50_jn

  dsum_send(31) = mdot50_js
  dsum_send(32) = phi50_js
  dsum_send(33) = edot_em50_js
  dsum_send(34) = edot_pake50_js
  dsum_send(35) = edot_en50_js

  call MPI_REDUCE(dsum_send,dsum_recv,nqttys,MPI_DOUBLE_PRECISION, &
       MPI_SUM,0,icomm,ierrmpi)


  if (mype .eq. 0) then

     ! use roundoff because ascii can give underflows which makes the data unreadable
     mdot      = roundoff(dsum_recv(1),1.0d-99)
     phi       = roundoff(dsum_recv(2),1.0d-99)
     edot_em   = roundoff(dsum_recv(3),1.0d-99)
     edot_pake = roundoff(dsum_recv(4),1.0d-99)
     edot_en   = roundoff(dsum_recv(5),1.0d-99)

     mdot50      = roundoff(dsum_recv(6),1.0d-99)
     phi50       = roundoff(dsum_recv(7),1.0d-99)
     edot_em50   = roundoff(dsum_recv(8),1.0d-99)
     edot_pake50 = roundoff(dsum_recv(9),1.0d-99)
     edot_en50   = roundoff(dsum_recv(10),1.0d-99) 

     mdot50_j      = roundoff(dsum_recv(11),1.0d-99)
     phi50_j       = roundoff(dsum_recv(12),1.0d-99)
     edot_em50_j   = roundoff(dsum_recv(13),1.0d-99)
     edot_pake50_j = roundoff(dsum_recv(14),1.0d-99)
     edot_en50_j   = roundoff(dsum_recv(15),1.0d-99) 

     mdot50_w      = roundoff(dsum_recv(16),1.0d-99)
     phi50_w       = roundoff(dsum_recv(17),1.0d-99)
     edot_em50_w   = roundoff(dsum_recv(18),1.0d-99)
     edot_pake50_w = roundoff(dsum_recv(19),1.0d-99)
     edot_en50_w   = roundoff(dsum_recv(20),1.0d-99) 

     mdot50_disk      = roundoff(dsum_recv(21),1.0d-99)
     phi50_disk       = roundoff(dsum_recv(22),1.0d-99)
     edot_em50_disk   = roundoff(dsum_recv(23),1.0d-99)
     edot_pake50_disk = roundoff(dsum_recv(24),1.0d-99)
     edot_en50_disk   = roundoff(dsum_recv(25),1.0d-99) 

     mdot50_jn      = roundoff(dsum_recv(26),1.0d-99)
     phi50_jn       = roundoff(dsum_recv(27),1.0d-99)
     edot_em50_jn   = roundoff(dsum_recv(28),1.0d-99)
     edot_pake50_jn = roundoff(dsum_recv(29),1.0d-99)
     edot_en50_jn   = roundoff(dsum_recv(30),1.0d-99) 

     mdot50_js      = roundoff(dsum_recv(31),1.0d-99)
     phi50_js       = roundoff(dsum_recv(32),1.0d-99)
     edot_em50_js   = roundoff(dsum_recv(33),1.0d-99)
     edot_pake50_js = roundoff(dsum_recv(34),1.0d-99)
     edot_en50_js   = roundoff(dsum_recv(35),1.0d-99) 

     ! generate filename
     write(filename,"(a,a,a)") trim(filenameout),TRIM("_analysis"),".csv"
     INQUIRE(FILE=filename, EXIST=file_exists)
     if (.not. file_exists) then
        open(unit=unitanalysis,file=filename,status='unknown',position='append')
        write(unitanalysis,"(a,a)",advance='no') trim('# t mdot phi edot_em edot_pake edot_en'), ' '
        write(unitanalysis,"(a,a)",advance='no') trim('mdot50 phi50 edot_em50 edot_pake50 edot_en50'), ' '
        write(unitanalysis,"(a,a)",advance='no') trim('mdot50_j phi50_j edot_em50_j edot_pake50_j edot_en50_j'), ' '
        write(unitanalysis,"(a,a)",advance='no') trim('mdot50_w phi50_w edot_em50_w edot_pake50_w edot_en50_w'), ' '
        write(unitanalysis,"(a,a)",advance='no') trim('mdot50_disk phi50_disk edot_em50_disk edot_pake50_disk edot_en50_disk'), ' '
        write(unitanalysis,"(a,a)",advance='no') trim('mdot50_jn phi50_jn edot_em50_jn edot_pake50_jn edot_en50_jn'), ' '
        write(unitanalysis,"(a)",advance='yes')  trim('mdot50_js phi50_js edot_em50_js edot_pake50_js edot_en50_js')
     else
        open(unit=unitanalysis,file=filename,status='unknown',position='append')
     end if
     write(unitanalysis,'(11(ES14.6))', advance='no') t, mdot, phi, edot_em, edot_pake, edot_en, mdot50, phi50, edot_em50, edot_pake50, edot_en50
     write(unitanalysis,'(5(ES14.6))', advance='no')  mdot50_j, phi50_j, edot_em50_j, edot_pake50_j, edot_en50_j
     write(unitanalysis,'(5(ES14.6))', advance='no')  mdot50_w, phi50_w, edot_em50_w, edot_pake50_w, edot_en50_w
     write(unitanalysis,'(5(ES14.6))', advance='no')  mdot50_disk, phi50_disk, edot_em50_disk, edot_pake50_disk, edot_en50_disk
     write(unitanalysis,'(5(ES14.6))', advance='no')  mdot50_jn, phi50_jn, edot_em50_jn, edot_pake50_jn, edot_en50_jn
     write(unitanalysis,'(5(ES14.6))', advance='yes')  mdot50_js, phi50_js, edot_em50_js, edot_pake50_js, edot_en50_js
     close(unitanalysis)

  end if


  !--------------------------------------
  ! Hack the slicing tool:
  !--------------------------------------
  sliceascii = sliceascii_save
  saveprim   = saveprim_save
  nwauxio    = nwauxio_save
  !--------------------------------------
  ! Done resetting globals
  !--------------------------------------

end subroutine write_analysis
!=============================================================================
