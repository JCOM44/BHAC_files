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

module mod_shell
  !
  ! Contains routines for building an arbitrary shell and integrate on it
  ! Only for 2D and 3D.
  ! 03.Oct. 2017 Oliver Porth
  !
  use mod_physicaldata
  implicit none


  integer, parameter                               :: nshellshi = 512

  type tshell
     double precision                             :: r=0.0d0, dxShell2&
        =0.0d0,dxShell3=0.0d0
     integer                                      ::  ixGmin2=-1,ixGmax2&
        =-1, ixGmin3=-1,ixGmax3=-1
     integer                                      :: nwmin=-1,nwmax=-1 !allocated index-range
     double precision, dimension(:,:,:,:), allocatable :: w
     double precision, dimension(:,:,:,:), allocatable :: x, xShell
     logical                                       :: allocated=.false.
  end type tshell

  type(walloc), save,dimension(ngridshi)           :: gridvars
  type(tshell), save,dimension(nshellshi)          :: shell
  
 !=============================================================================
contains
 !=============================================================================
  subroutine alloc_shell(is,ixGmax2,ixGmax3,nwmax)

    use mod_amrvacdef

    integer, intent(in)                            :: is,ixGmax2,ixGmax3,nwmax
    ! .. local ..
    integer                                        :: ixGmin2,ixGmin3
 !-----------------------------------------------------------------------------

    ixGmin2=1;ixGmin3=1;
    
    if (.not. shell(is)%allocated) then
       allocate(shell(is)%w(1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:nwmax))
       allocate(shell(is)%x(1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:ndim))
       allocate(shell(is)%xShell(1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:ndim))
       shell(is)%ixGmin2=ixGmin2;shell(is)%ixGmin3=ixGmin3;
       shell(is)%ixGmax2=ixGmax2;shell(is)%ixGmax3=ixGmax3;
       shell(is)%nwmin=1;       shell(is)%nwmax=nwmax
       shell(is)%allocated = .true.
    else
       call mpistop('Trying to allocate already allocated shell')
    end if
    
  end subroutine alloc_shell
 !=============================================================================
  subroutine dealloc_shell(is)

    use mod_amrvacdef
    
    integer, intent(in)                            :: is
 !-----------------------------------------------------------------------------
    
    if (shell(is)%allocated) then
       deallocate(shell(is)%w)
       deallocate(shell(is)%x)
       deallocate(shell(is)%xShell)
       shell(is)%ixGmin2=-1;shell(is)%ixGmin3=-1;
       shell(is)%ixGmax2=-1;shell(is)%ixGmax3=-1;
       shell(is)%nwmin=-1;       shell(is)%nwmax=-1
       shell(is)%allocated = .false.
       shell(is)%r=0.0d0
       shell(is)%dxShell2=0.0d0;shell(is)%dxShell3=0.0d0;
    else
       call mpistop('Trying to deallocate already deallocated shell')
    end if
    
  end subroutine dealloc_shell
 !=============================================================================
  subroutine init_gridvars(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     nwmax)
    
    use mod_amrvacdef

    integer, intent(in)                             :: ixGmin1,ixGmin2,&
       ixGmin3,ixGmax1,ixGmax2,ixGmax3, nwmax
    ! .. local ..
    integer                                         :: igrid, iigrid
    double precision,dimension(0:nwmax)             :: normconv 
 !-----------------------------------------------------------------------------

    if (nwmax .lt. nw) call mpistop('init_gridvars: nwmax needs to be .ge. &
       than nw')
    
    do iigrid=1,igridstail; igrid=igrids(iigrid);         
       call set_tmpGlobals(igrid)
       
       call alloc_pw(gridvars(igrid),ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
          ixGmax3,1,nwmax)
       gridvars(igrid)%w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
          1:nw) = pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nw)

       if(nwmax>nw) call specialvar_output(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
          ixGhi3,ixMlo1-1,ixMlo2-1,ixMlo3-1,ixMhi1+1,ixMhi2+1,ixMhi3+1,nwmax,&
            gridvars(igrid)%w,ps(igrid),normconv)

       if(saveprim)  call primitive(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
          ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
            gridvars(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nw),&
               px(igrid)%x)

    end do

  end subroutine init_gridvars
 !=============================================================================
  subroutine finish_gridvars()

    use mod_amrvacdef

    ! .. local ..
    integer             :: iigrid, igrid
 !-----------------------------------------------------------------------------
    
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       call dealloc_pw(gridvars(igrid))
    end do

  end subroutine finish_gridvars
 !=============================================================================
  subroutine fill_shell(is,rshell,SphToCoord)

    use mod_interpolate, only: find_point_ipe, interpolate_var
    use mod_amrvacdef
    
    integer, intent(in)                                        :: is
    double precision, intent(in)                               :: rshell

    interface
       subroutine SphToCoord(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xKS,xCKS)
         integer,intent(in)                                     :: ixImin1,&
            ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,&
            ixOmax1,ixOmax2,ixOmax3
         double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
            ixImin3:ixImax3,1:3), intent(in)   :: xKS
         double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
            ixImin3:ixImax3,1:3), intent(out)  :: xCKS
       end subroutine SphToCoord
    end interface
    
    ! .. local ..
    integer                                          :: ix2,ix3, N3, N2,&
        igrid, ipe, iw
!    double precision, dimension(:,:,:,:), allocatable :: tmpw
    integer, allocatable, dimension(:)               :: rcv_n_points_from_ipe,&
        ipe_with_data, iipe_of_ipe
    double precision, allocatable, dimension(:,:)    :: sendbuff
    double precision, allocatable, dimension(:,:,:)  :: rcvbuff
    integer                                          :: iipe,&
        n_processes_to_receive_from, max_npoints_to_receive, isnd, ircv
    integer, allocatable, dimension(:)               :: rcvrqst, buffer_offset
    integer                                          :: sndrqst(1),&
        one_if_sending, ipebuffer
    integer, allocatable, dimension(:,:)           :: ipe_of_index
 !-----------------------------------------------------------------------------
    associate (s=>shell(is))

      if (.not. s%allocated) call mpistop('fill_shell: trying to fill &
         unallocated shell')

!      allocate(tmpw,mold=s%w)
      allocate(rcv_n_points_from_ipe(0:npe-1), ipe_with_data(0:npe-1),&
          sendbuff(s%nwmin:s%nwmax,1:(s%ixGmax2-s%ixGmin2+1)*&
         (s%ixGmax3-s%ixGmin3+1)))
      allocate(rcvrqst(1:npe-1),iipe_of_ipe(1:npe-1),ipe_of_index&
         (s%ixGmin2:s%ixGmax2,s%ixGmin3:s%ixGmax3))
      sndrqst = MPI_REQUEST_NULL; rcvrqst = MPI_REQUEST_NULL
      
      s%r = rshell
      s%w = 0.0d0
      
      ! First get the coordinates:
      
      N3 = s%ixGmax3-s%ixGmin3+1
      s%dxshell3   = dpi/dble(N3-1)
     
      N2 = s%ixGmax2-s%ixGmin2+1
      s%dxshell2 = 2.0d0*dpi/dble(N2-1)
           

      do ix2=  s%ixGmin2,s%ixGmax2 
do ix3=  s%ixGmin3,s%ixGmax3 
         
         s%xShell(1,ix2,ix3,3)=(ix3-1)*s%dxShell3
        
         s%xShell(1,ix2,ix3,2)=(ix2-1)*s%dxShell2
        
         s%xShell(1,ix2,ix3,1)=rshell
      end do 
end do 

      call SphToCoord(1,s%ixGmin2,s%ixGmin3,1,s%ixGmax2,s%ixGmax3,1,s%ixGmin2,&
         s%ixGmin3,1,s%ixGmax2,s%ixGmax3,s%xShell,s%x)

      
      isnd = 0; one_if_sending=0
      ircv = 0; rcv_n_points_from_ipe = 0
      do ix2= s%ixGmin2,s%ixGmax2 
do ix3= s%ixGmin3,s%ixGmax3 

         !Finding processor and blocks for each point
         call find_point_ipe(s%x(1,ix2,ix3,1:ndim),igrid,ipe)
         ipe_of_index(ix2,ix3) = ipe
         
         if (ipe.ne.0.and.mype.eq.0) &
              rcv_n_points_from_ipe(ipe) = rcv_n_points_from_ipe(ipe) + 1

         
         if(ipe.eq.mype) then
            
            do iw=s%nwmin,s%nwmax
               call interpolate_var(igrid,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
                  ixGhi3,gridvars(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
                  ixGlo3:ixGhi3,iw),&
                    px(igrid)%x,s%x(1,ix2,ix3,1:ndim),s%w(1,ix2,ix3,&
                       iw))               
            end do
            
            if (ipe.ne.0) then
               isnd      = isnd + 1
               do iw=s%nwmin,s%nwmax
                  sendbuff(iw,isnd) = s%w(1,ix2,ix3,iw)
               end do
            end if

         end if
         
      end do
end do
         

      if (isnd > 0) then
         call MPI_ISEND(sendbuff,(s%nwmax-s%nwmin+1)*isnd,&
            MPI_DOUBLE_PRECISION,0,mype,icomm,sndrqst,ierrmpi)
         one_if_sending = 1
      end if

      if (mype .eq. 0) then
         n_processes_to_receive_from = 0; max_npoints_to_receive = 0
         do ipe=1,npe-1
            if (rcv_n_points_from_ipe(ipe) /= 0) then
               n_processes_to_receive_from = n_processes_to_receive_from+1
               ipe_with_data(n_processes_to_receive_from) = ipe
               if (rcv_n_points_from_ipe(ipe) > max_npoints_to_receive) &
                  max_npoints_to_receive = rcv_n_points_from_ipe(ipe)
            end if
         end do
         allocate(rcvbuff(s%nwmin:s%nwmax,1:max_npoints_to_receive,&
            1:n_processes_to_receive_from))

         do iipe=1,n_processes_to_receive_from; ipe=ipe_with_data(iipe)
            iipe_of_ipe(ipe) = iipe
            ircv = ircv + 1
            call MPI_IRECV(rcvbuff(:,:,iipe),rcv_n_points_from_ipe(ipe)*&
               (s%nwmax-s%nwmin+1),MPI_DOUBLE_PRECISION,ipe,ipe,icomm,&
               rcvrqst(ircv),ierrmpi)
         end do
      end if

      call MPI_WAITALL(one_if_sending,sndrqst,MPI_STATUSES_IGNORE,ierrmpi)
      call MPI_WAITALL(ircv,rcvrqst,MPI_STATUSES_IGNORE,ierrmpi)

      if (mype .eq. 0) then
         allocate(buffer_offset(1:n_processes_to_receive_from))
         buffer_offset = 0
         do ix2= s%ixGmin2,s%ixGmax2 
do ix3= s%ixGmin3,s%ixGmax3 

            ipe = ipe_of_index(ix2,ix3)
            if (ipe.eq.0) cycle
            
            iipe = iipe_of_ipe(ipe)
            buffer_offset(iipe) = buffer_offset(iipe) + 1
            do iw=s%nwmin,s%nwmax
               s%w(1,ix2,ix3,iw) = rcvbuff(iw,buffer_offset(iipe),iipe)
            end do

         end do
end do
      end if
      
      
      ! Reduce to head-node:
!      if (mype==0) then
 !call MPI_REDUCE(MPI_IN_PLACE,s%w,size(s%w),MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)
!      else
 !call MPI_REDUCE(s%w,s%w,size(s%w),MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)
!      end if
!      call MPI_REDUCE(s%w,tmpw,size(tmpw),MPI_DOUBLE_PRECISION,MPI_SUM,0,&
!           icomm,ierrmpi)
!      if (mype==0) s%w=tmpw
!      call MPI_BARRIER(icomm, ierrmpi)
            
    end associate
  end subroutine fill_shell
 !=============================================================================
  subroutine write_shell

    use mod_amrvacdef

    ! Writes a topological sphere
    ! 03.Oct 2017
    integer                                   :: is
    logical, save                             :: firstshell=.true.
 !-----------------------------------------------------------------------------

    if (firstshell) then
       ishell=shellNext
       firstshell=.false.
    end if

    call init_gridvars(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,nw+nwauxio)
    do is=1,nshells
       call put_shell(shellcoord(is))
    end do
    call finish_gridvars

    ishell=ishell+1

  end subroutine write_shell
 !=============================================================================
  subroutine put_shell(rshell)

    use mod_metric, only: SPToCoord
    use mod_amrvacdef

    double precision, intent(in)               :: rshell
 !-----------------------------------------------------------------------------

    call alloc_shell(1, nxShell2,nxShell3, nw+nwauxio)

    call fill_shell(1, rshell, SPToCoord)

    select case(shell_type)
    case('csv')
       call put_shell_csv(shell(1),ishell)
    case('vtu')
       call put_shell_vtu(shell(1),ishell)       
    case default
       call mpistop('put_shell: unknown shell_type')
    end select
    
    call dealloc_shell(1)

  end subroutine put_shell
 !=============================================================================
  subroutine put_shell_csv(shell,isout)

    use mod_tocart
    use mod_amrvacdef

    type(tshell), intent(in)      :: shell  ! the shell structure
    integer, intent(in)           :: isout  ! the output index
    ! .. local ..
    character(len=1024)           :: filename, xlabel
    logical                       :: fileopen
    character(len=10)             :: wnamei(1:nw+nwauxio),xandwnamei&
       (1:ndim+nw+nwauxio)
    character(len=1024)           :: outfilehead
    integer                       :: iw, idir, ix2,ix3
    character(len=1024)           :: line, data
    double precision, dimension(1,shell%ixGmin2:shell%ixGmax2,&
       shell%ixGmin3:shell%ixGmax3,shell%nwmin:shell%nwmax) :: wCart
    double precision, dimension(1,shell%ixGmin2:shell%ixGmax2,&
       shell%ixGmin3:shell%ixGmax3,1:ndim)          :: xCart
    double precision, parameter :: minvalue = 1.0d-99, maxvalue = 1.0d+99
    double precision            :: roundoff_minmax
 !-----------------------------------------------------------------------------

    if (mype==0) then

       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Open the file:
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       inquire(unitshell,opened=fileopen)
       if(.not.fileopen)then
          ! generate filename: 
          write(xlabel,"(D10.3)") shell%r
          if(shell%r>=zero)then
             write(xlabel(1:1),"(a)") "+"
          end if
          write(filename,"(a,i4.4,a)") TRIM(filenameout)//'_r'//trim(xlabel)//'_n',isout,'.csv'
          open(unitshell,file=filename,status='unknown',form='formatted')
       end if

       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Write the header:
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       call getheadernames(wnamei,xandwnamei,outfilehead)
       line=''
       do iw=1,ndim+nw+nwauxio-1
          line = trim(line)//trim(xandwnamei(iw))//', '
       end do
       line = trim(line)//trim(xandwnamei(ndim+nw+nwauxio))

       write(unitshell,'(a)') trim(line)

       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Calculate Cartesian components:
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       call cartesian_covariant(1,shell%ixGmin2,shell%ixGmin3,1,shell%ixGmax2,&
          shell%ixGmax3,1,shell%ixGmin2,shell%ixGmin3,1,shell%ixGmax2,&
          shell%ixGmax3,shell%w,shell%x,&
            wCart,xCart,wCoord_is_primitive=saveprim)

       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Write data to file:
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       do ix2=  shell%ixGmin2,shell%ixGmax2 
do ix3=  shell%ixGmin3,shell%ixGmax3 

          line=''
          do idir=1,ndim
             write(data,"(es14.6,a)") roundoff_minmax(xCart(1,ix2,ix3,idir),&
                minvalue,maxvalue),', '
             line=trim(line)//trim(data)
          end do
          do iw=1,nw+nwauxio-1
             write(data,"(es14.6,a)") roundoff_minmax(wCart(1,ix2,ix3,iw),&
                minvalue,maxvalue),', '
             line=trim(line)//trim(data)
          end do
          write(data,"(es14.6)") roundoff_minmax(wCart(1,ix2,ix3,nw+nwauxio),&
             minvalue,maxvalue)
          line=trim(line)//trim(data)

          write(unitshell,"(a)") trim(line)
       
       end do
end do

       close(unitshell)
       
    end if! mype.eq.0
       
  end subroutine put_shell_csv
 !=============================================================================
  subroutine put_shell_vtu(shell,isout)

    use mod_tocart
    use mod_amrvacdef

    type(tshell), intent(in)      :: shell  ! the shell structure
    integer, intent(in)           :: isout  ! the output index
    ! .. local ..
    character(len=1024)           :: filename, xlabel
    logical                       :: fileopen
    character(len=10)             :: wnamei(1:nw+nwauxio),xandwnamei&
       (1:ndim+nw+nwauxio)
    character(len=1024)           :: outfilehead
    integer                       :: iw, ix2,ix3
    character(len=1024)           :: line, data
    double precision, dimension(1,shell%ixGmin2:shell%ixGmax2,&
       shell%ixGmin3:shell%ixGmax3,shell%nwmin:shell%nwmax) :: wCart
    double precision, dimension(1,shell%ixGmin2:shell%ixGmax2,&
       shell%ixGmin3:shell%ixGmax3,1:ndim)          :: xCart
    double precision, parameter :: minvalue = 1.0d-99, maxvalue = 1.0d+99
    double precision            :: roundoff_minmax
    integer                     :: nx2,nx3, nxC2,nxC3, nc, np, icell
    double precision            :: x_VTK(1:3)
    integer                     :: VTK_type
 !-----------------------------------------------------------------------------

    if (mype==0) then

       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Open the file:
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       inquire(unitshell,opened=fileopen)
       if(.not.fileopen)then
          ! generate filename: 
          write(xlabel,"(D10.3)") shell%r
          if(shell%r>=zero)then
             write(xlabel(1:1),"(a)") "+"
          end if
          write(filename,"(a,i4.4,a)") TRIM(filenameout)//'_r'//trim(xlabel)//'_n',isout,'.vtu'
          open(unitshell,file=filename,status='unknown',form='formatted')
       end if

       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Get the header:
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       call getheadernames(wnamei,xandwnamei,outfilehead)


       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Calculate Cartesian components:
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       call cartesian_covariant(1,shell%ixGmin2,shell%ixGmin3,1,shell%ixGmax2,&
          shell%ixGmax3,1,shell%ixGmin2,shell%ixGmin3,1,shell%ixGmax2,&
          shell%ixGmax3,shell%w,shell%x,&
            wCart,xCart,wCoord_is_primitive=saveprim)


       
       nx2=shell%ixGmax2-shell%ixGmin2;nx3=shell%ixGmax3-shell%ixGmin3; !cells
       nxC2=nx2+1;nxC3=nx3+1;                        ! corners
       nc=nx2*nx3
       np=nxC2*nxC3

       
       ! generate xml header
       write(unitshell,'(a)')'<?xml version="1.0"?>'
       write(unitshell,'(a)',advance='no') '<VTKFile type="UnstructuredGrid"'
       
        write(unitshell,'(a)')' version="0.1" byte_order="LittleEndian">'
       write(unitshell,'(a)')'  <UnstructuredGrid>'
       write(unitshell,'(a)')'<FieldData>'
       write(unitshell,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
            'NumberOfTuples="1" format="ascii">'
       write(unitshell,*) real(t)
       write(unitshell,'(a)')'</DataArray>'
       write(unitshell,'(a)')'</FieldData>'
       
       
       ! write out one VTK PIECE
       write(unitshell,'(a,i7,a,i7,a)') &
            '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'

       
       !==============================
       ! output Pointdata
       !==============================
       write(unitshell,'(a)')'<PointData>'
       do iw=1,nw+nwauxio
          if(iw<=nw) then 
             if(.not.writew(iw)) cycle
          endif
          write(unitshell,'(a,a,a)')&
               '<DataArray type="Float64" Name="',TRIM(wnamei(iw)),'" format&
                  ="ascii">'
          write(unitshell,'(200(1pe14.6))') ((roundoff_minmax(wCart(1,ix2,ix3,&
             iw),minvalue,maxvalue),ix2=shell%ixGmin2,shell%ixGmax2),ix3&
             =shell%ixGmin3,shell%ixGmax3)
          write(unitshell,'(a)')'</DataArray>'
       end do
       write(unitshell,'(a)')'</PointData>'
       !==============================
       ! Done: Output Pointdata
       !==============================

       !==============================
       ! output Cornerpoints
       !==============================
       write(unitshell,'(a)')'<Points>'
       write(unitshell,'(a)')'<DataArray type="Float32" NumberOfComponents&
          ="3" format="ascii">'
 !write cell corner coordinates in a backward dimensional loop, always 3D output
       do ix3=shell%ixGmin3,shell%ixGmax3 
 do ix2=shell%ixGmin2,shell%ixGmax2 
            x_VTK(1:3)=zero;
            x_VTK(1:ndim)=xCart(1,ix2,ix3,1:ndim)
            write(unitshell,'(3(1pe14.6))') x_VTK
      end do 
end do 
      write(unitshell,'(a)')'</DataArray>'
      write(unitshell,'(a)')'</Points>'
      !==============================
      ! Done: output Cornerpoints
      !==============================


      !==============================
      ! cell Metainformation
      !==============================
      write(unitshell,'(a)')'<Cells>'

      ! connectivity part
      write(unitshell,'(a)')'<DataArray type="Int32" Name="connectivity" &
         format="ascii">'

       do ix3=1,nx3
 do ix2=1,nx2
      
      write(unitshell,'(4(i7,1x))')(ix3-1)*nxC2+ix2-1, &
           (ix3-1)*nxC2+ix2,ix3*nxC2+ix2-1,ix3*nxC2+ix2
      end do
 end do
              
      write(unitshell,'(a)')'</DataArray>'

      ! offsets data array
      write(unitshell,'(a)')'<DataArray type="Int32" Name="offsets" format&
         ="ascii">'
      do icell=1,nc
         write(unitshell,'(i7)') icell*(2**(3-1))
      end do
      write(unitshell,'(a)')'</DataArray>'

      ! VTK cell type data array
      write(unitshell,'(a)')'<DataArray type="Int32" Name="types" format&
         ="ascii">'
      ! VTK_LINE=3; VTK_PIXEL=8; VTK_VOXEL=11 -> vtk-syntax
      
       VTK_type=8 
      do icell=1,nc
         write(unitshell,'(i2)') VTK_type
      enddo
      write(unitshell,'(a)')'</DataArray>'
      
      write(unitshell,'(a)')'</Cells>'
      !==============================
      ! Done: cell Metainformation
      !==============================
      write(unitshell,'(a)')'</Piece>'

      
       write(unitshell,'(a)')'</UnstructuredGrid>'
       write(unitshell,'(a)')'</VTKFile>'
       close(unitshell)

    end if

  end subroutine put_shell_vtu
 !=============================================================================
  subroutine shell_integrate(s,iw,gfac,int,imask)

    use mod_amrvacdef
    
    type(tshell), intent(in)      :: s      ! the shell structure
    integer, intent(in)           :: iw     ! the variable number to integrate
    integer, optional, intent(in) :: imask  ! the optional mask
    double precision, intent(out) :: int    ! the integral over the shell

    interface
       subroutine gfac(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xShell,gg) !the geometric factor
         integer,intent(in)                                     :: ixImin1,&
            ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,&
            ixOmax1,ixOmax2,ixOmax3
         double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
            ixImin3:ixImax3,1:3), intent(in)   :: xShell
         double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
            ixImin3:ixImax3), intent(out)        :: gg
       end subroutine gfac
    end interface

    ! .. local ..
    double precision, dimension(1,s%ixGmin2:s%ixGmax2,s%ixGmin3:s%ixGmax3)   &
       :: gg, mask
    integer                                                    :: ix1,ix2,ix3,&
        ip2,ip3, itmp2,itmp3, icnt
    double precision                                           :: tmp
 !-----------------------------------------------------------------------------

    if (mype.ne.0) then

       call mpistop("shell_integrate: called for mype.ne.0,&
           currently only headnode has the full shell !")

    else

       if (present(imask)) then
          mask(:,:,:) = s%w(:,:,:,imask)
       else
          mask(:,:,:) = 1.0d0
       end if
       
       call gfac(1,s%ixGmin2,s%ixGmin3,1,s%ixGmax2,s%ixGmax3,1,s%ixGmin2,&
          s%ixGmin3,1,s%ixGmax2,s%ixGmax3,s%xshell,gg)


       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       int = 0.0d0
       do ix2= s%ixGmin2,s%ixGmax2-1
do ix3= s%ixGmin3,s%ixGmax3-1

       ! Average to center (not very fast, but does it matter?):
       icnt=0
       tmp=0.0d0
       do ip2= 0,1
do ip3= 0,1
       itmp2=ix2+ip2;itmp3=ix3+ip3;
       icnt=icnt+1
       tmp = tmp + s%w(1,itmp2,itmp3,iw) * mask(1,itmp2,itmp3) * gg(1,itmp2,&
          itmp3)
       end do
end do
       tmp = tmp/dble(icnt)

       int = int + tmp * s%dxShell2*s%dxShell3

       end do
end do
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
    end if

  end subroutine shell_integrate
 !=============================================================================
end module mod_shell


!=============================================================================
