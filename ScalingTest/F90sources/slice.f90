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
module mod_slice
  implicit none

contains
  !=============================================================================
  subroutine write_slice
    use mod_amrvacdef
    ! Writes a D-1 slice .dat-file with proper Morton order 
    ! by Oliver Porth
    ! 22.Nov 2011
    integer :: islice
    logical, save :: firstslice=.true.
    !-----------------------------------------------------------------------------
    if (firstslice) then
       slice=slicenext
       firstslice=.false.
    end if

    do islice=1,nslices
       call put_slice(slicedir(islice),slicecoord(islice))
    end do

    slice=slice+1
  end subroutine write_slice
  !=============================================================================
  subroutine put_slice(dir,xslice)
    use mod_forest, only: Morton_sub_start, Morton_sub_stop
    use mod_amrvacdef
    ! Writes a D-1 slice .dat-file with proper Morton order 
    ! Can also write csv files by setting sliceascii=.true.
    ! For ONED simulations, the output will be appended to one csv-file per slice
    ! In the latter two cases, the slices are sensitive to the saveprim switch
    ! Thus csv-files with primitive variables are obtained.  
    ! by Oliver Porth
    ! 22.Nov 2011
    integer, intent(in) :: dir
    double precision, intent(in) :: xslice
    ! .. local ..
    integer :: Njgrid, jgrid
    integer, dimension(ndim-1) :: ixsubGlo, ixsubGhi
    integer, dimension(ndim-1) :: ixsubMlo, ixsubMhi
    integer :: size_subblock_io, nx1,nx2,nx3, slice_fh, nwexpand
    integer, dimension(ndim) :: sizes, subsizes, start
    double precision,dimension(0:nw+nwauxio)          :: normconv 
    integer(kind=MPI_ADDRESS_KIND):: size_double, lb
    !-----------------------------------------------------------------------------
    ! Preamble: 
    call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,size_double,ierrmpi)
    nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;nx3=ixMhi3-ixMlo3+1;
    slice_fh=unitslice

    if (slice_type.ne.'dat') then
       nwexpand = nwauxio
    else
       nwexpand = 0
    end if

    ! Do a last consistency check:
    select case(dir)
       case(1)
       if(xslice<xprobmin1.or.xslice>xprobmax1) &
            call mpistop("slice out of bounds")
       
       case(2)
       if(xslice<xprobmin2.or.xslice>xprobmax2) &
            call mpistop("slice out of bounds")
       
       case(3)
       if(xslice<xprobmin3.or.xslice>xprobmax3) &
            call mpistop("slice out of bounds")
       
    end select

    ! Traverse the forest and fill nodes:
    call select_slice(dir,xslice,.false.,slice_fh,normconv)

    ! Create the MPI-datatype and select indices:
    
    select case(dir)
    case (1)
       sizes(1) = ixGhi2; sizes(2) = ixGhi3;
       ixsubGlo(1) = ixGlo2; ixsubGlo(2) = ixGlo3;
       ixsubGhi(1) = ixGhi2; ixsubGhi(2) = ixGhi3;
       subsizes(1)=nx2;subsizes(2)=nx3;
       start(1)=ixMlo2-1;start(2)=ixMlo3-1;
       size_subblock_io=nx2*nx3*(nw+nwexpand)*size_double
    case (2)
       sizes(1) = ixGhi1; sizes(2) = ixGhi3;
       ixsubGlo(1) = ixGlo1; ixsubGlo(2) = ixGlo3;
       ixsubGhi(1) = ixGhi1; ixsubGhi(2) = ixGhi3;
       subsizes(1)=nx1;subsizes(2)=nx3;
       start(1)=ixMlo1-1;start(2)=ixMlo3-1;
       size_subblock_io=nx1*nx3*(nw+nwexpand)*size_double
    case (3)
       ixsubGlo(1) = ixGlo1; ixsubGlo(2) = ixGlo2;
       ixsubGhi(1) = ixGhi1; ixsubGhi(2) = ixGhi2;
       sizes(1) = ixGhi1; sizes(2) = ixGhi2;
       subsizes(1)=nx1;subsizes(2)=nx2;
       start(1)=ixMlo1-1;start(2)=ixMlo2-1;
       size_subblock_io=nx1*nx2*(nw+nwexpand)*size_double
    case default
       call mpistop("slice direction not clear in put_slice")
    end select
   
    
    

    
    ixsubMlo(2-1) = ixsubGlo(2-1)+dixB;ixsubMlo(3-1) = ixsubGlo(3-1)+dixB;
    ixsubMhi(2-1) = ixsubGhi(2-1)-dixB;ixsubMhi(3-1) = ixsubGhi(3-1)-dixB;
   

    sizes(ndim)=nw+nwexpand
    subsizes(ndim)=nw+nwexpand
    start(ndim)=0

    ! Types for center variables:
    call MPI_TYPE_CREATE_SUBARRAY(ndim,sizes,subsizes,start,&
        MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, type_subblock_io,ierrmpi)
    call MPI_TYPE_COMMIT(type_subblock_io,ierrmpi)

    sizes(ndim)=3
    subsizes(ndim)=3
    start(ndim)=0
    call MPI_TYPE_CREATE_SUBARRAY(ndim,sizes,subsizes,start,&
        MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, type_subblock_x_io,ierrmpi)
    call MPI_TYPE_COMMIT(type_subblock_x_io,ierrmpi)


    ! Types for corner variables:
    subsizes(1:ndim-1) = subsizes(1:ndim-1) + 1
    start(1:ndim-1)    = start(1:ndim-1) - 1 
    sizes(ndim)=nw+nwexpand
    subsizes(ndim)=nw+nwexpand
    start(ndim)=0
    
    call MPI_TYPE_CREATE_SUBARRAY(ndim,sizes,subsizes,start,&
        MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, type_subblockC_io,ierrmpi)
    call MPI_TYPE_COMMIT(type_subblockC_io,ierrmpi)

    sizes(ndim)=3
    subsizes(ndim)=3
    start(ndim)=0
    call MPI_TYPE_CREATE_SUBARRAY(ndim,sizes,subsizes,start,&
        MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, type_subblockC_x_io,ierrmpi)
    call MPI_TYPE_COMMIT(type_subblockC_x_io,ierrmpi)

    
    ! local number of sub-grids:
    Njgrid=Morton_sub_stop(mype)-Morton_sub_start(mype)+1

    ! Now output using various schemes: 
    if (ndim==1) then 
       call put_slice_zerod
    else
       select case(slice_type)
       case ('csv')
          call put_slice_csv
       case ('dat')
          call put_slice_dat
       case ('vtu', 'vtuCC')
          call put_slice_vtu
       end select
    end if

    do jgrid=1,Njgrid
       call dealloc_subnode(jgrid)
    end do

    call MPI_TYPE_FREE(type_subblock_io,ierrmpi)
    call MPI_TYPE_FREE(type_subblock_x_io,ierrmpi)
    call MPI_TYPE_FREE(type_subblockC_io,ierrmpi)
    call MPI_TYPE_FREE(type_subblockC_x_io,ierrmpi)


  contains
    !=============================================================================
    ! internal procedures
    !=============================================================================
    subroutine put_slice_vtu

      character(len=1024) :: filename, xlabel
      character(len=79)   :: xxlabel
      logical             :: fileopen
      character(len=10) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
      character(len=1024) :: outfilehead
      integer :: status(MPI_STATUS_SIZE), ipe
      !-----------------------------------------------------------------------------

      if (mype==0) then

         inquire(slice_fh,opened=fileopen)
         if(.not.fileopen)then
            ! generate filename: 
            write(xlabel,"(D9.2)")xslice
            xxlabel=trim(xlabel)
            if(xslice>=zero)then
               write(xxlabel(1:1),"(a)") "+"
            endif
            write(filename,"(a,i1.1,a,i4.4,a)") TRIM(filenameout)//'_d',dir,'_x'//trim(xxlabel)//'_n',slice,'.vtu'
            open(slice_fh,file=filename,status='unknown',form='formatted')
         end if
         ! get and write the header: 
         call getheadernames(wnamei,xandwnamei,outfilehead)
         ! generate xml header
         write(slice_fh,'(a)')'<?xml version="1.0"?>'
         write(slice_fh,'(a)',advance='no') '<VTKFile type="UnstructuredGrid"'
         
          write(slice_fh,'(a)')' version="0.1" byte_order="LittleEndian">'
         write(slice_fh,'(a)')'  <UnstructuredGrid>'
         write(slice_fh,'(a)')'<FieldData>'
         write(slice_fh,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
            'NumberOfTuples="1" format="ascii">'
         write(slice_fh,*) real(t*normt)
         write(slice_fh,'(a)')'</DataArray>'
         write(slice_fh,'(a)')'</FieldData>'

         ! write to file:
         do jgrid=1, Njgrid
            call write_slice_vtk(jgrid,slice_fh,wnamei)
         end do

         ! create a recv buffer using allocate, will be deallocated at the end of the routine:
         call alloc_subnode(Njgrid+1,dir,nwauxio)

      end if

      
            ! Also communicate the normconv array since processor zero might not have it yet:
      if (npe>1) then
         do ipe=1,npe-1
            do jgrid=1,Morton_sub_stop(ipe)-Morton_sub_start(ipe)+1
               itag=Morton_sub_start(ipe)+jgrid-1
               itag=itag*5
               if (ipe == mype ) then 
                  call MPI_SEND(px_sub(jgrid)%x,1,type_subblock_x_io,0,itag,&
                     icomm,ierrmpi)
                  call MPI_SEND(pw_sub(jgrid)%w,1,type_subblock_io,0,itag+1,&
                     icomm,ierrmpi)
                  call MPI_SEND(pxC_sub(jgrid)%x,1,type_subblockC_x_io,0,&
                     itag+2,icomm,ierrmpi)
                  call MPI_SEND(pwC_sub(jgrid)%w,1,type_subblockC_io,0,itag+3,&
                     icomm,ierrmpi)
                  call MPI_SEND(normconv,nw+nwauxio+1,MPI_DOUBLE_PRECISION,0,&
                     itag+4,icomm,ierrmpi)
               end if
               if (mype == 0) then
                  call MPI_RECV(px_sub(Njgrid+1)%x,1,type_subblock_x_io,ipe,&
                     itag,icomm,status,ierrmpi)
                  call MPI_RECV(pw_sub(Njgrid+1)%w,1,type_subblock_io,ipe,&
                     itag+1,icomm,status,ierrmpi)
                  call MPI_RECV(pxC_sub(Njgrid+1)%x,1,type_subblockC_x_io,ipe,&
                     itag+2,icomm,status,ierrmpi)
                  call MPI_RECV(pwC_sub(Njgrid+1)%w,1,type_subblockC_io,ipe,&
                     itag+3,icomm,status,ierrmpi)
                  call MPI_RECV(normconv,nw+nwauxio+1,MPI_DOUBLE_PRECISION,&
                     ipe,itag+4,icomm,status,ierrmpi)
                  call write_slice_vtk(Njgrid+1,slice_fh,wnamei)
               end if
            end do
         end do
      endif
      
     
      if (mype==0) then
         
         write(slice_fh,'(a)')'</UnstructuredGrid>'
         write(slice_fh,'(a)')'</VTKFile>'
         close(slice_fh)
         call dealloc_subnode(Njgrid+1)
         
      end if
 
    end subroutine put_slice_vtu
    !=============================================================================
    subroutine write_slice_vtk(jgrid,slice_fh,wnamei)

      integer, intent(in)           :: jgrid, slice_fh
      character(len=10), intent(in) :: wnamei(1:nw+nwauxio)
      
      ! .. local ..
      integer                       :: ixCmin1,ixCmin2,ixCmin3,ixCmax1,&
         ixCmax2,ixCmax3, ixCCmin1,ixCCmin2,ixCCmin3,ixCCmax1,ixCCmax2,&
         ixCCmax3, nc, np, iw
      integer                       :: nx1,nx2, nxC1,nxC2, icell, ix1,ix2
      double precision              :: x_VTK(1:3)
      integer                       :: VTK_type
      double precision, parameter   :: minvalue = 1.0d-99, maxvalue = 1.0d+99
      double precision              :: roundoff_minmax
 !-----------------------------------------------------------------------------
      ixCCmin1 = ixsubMlo(1);ixCCmin2 = ixsubMlo(2);
      ixCCmax1 = ixsubMhi(1);ixCCmax2 = ixsubMhi(2);
      ixCmin1  = ixsubMlo(1)-1;ixCmin2  = ixsubMlo(2)-1;
      ixCmax1  = ixsubMhi(1);ixCmax2  = ixsubMhi(2);

      nx1=ixCCmax1-ixCCmin1+1;nx2=ixCCmax2-ixCCmin2+1;
      nxC1=nx1+1;nxC2=nx2+1;
      nc=nx1*nx2      ! Number of cells per subgrid
      np=nxC1*nxC2     ! Number of corner points per subgrid
      
      
      ! we write out every grid as one VTK PIECE
      write(slice_fh,'(a,i7,a,i7,a)') &
           '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'

      !==============================
      ! celldata or pointdata?
      !==============================
      select case(slice_type)


      case('vtuCC') ! celldata
         write(slice_fh,'(a)')'<CellData>'
         do iw=1,nw+nwauxio
            if(iw<=nw) then 
               if(.not.writew(iw)) cycle
            endif
            write(slice_fh,'(a,a,a)')&
                 '<DataArray type="Float64" Name="',TRIM(wnamei(iw)),&
                    '" format="ascii">'
            write(slice_fh,'(200(1pe14.6))') ((roundoff_minmax(pw_sub&
               (jgrid)%w(ix1,ix2,iw)*normconv(iw),minvalue,maxvalue),ix1&
               =ixCCmin1,ixCCmax1),ix2=ixCCmin2,ixCCmax2)
            write(slice_fh,'(a)')'</DataArray>'
         enddo
         write(slice_fh,'(a)')'</CellData>'


      case('vtu') ! pointdata
         write(slice_fh,'(a)')'<PointData>'
         do iw=1,nw+nwauxio
            if(iw<=nw) then 
               if(.not.writew(iw)) cycle
            endif
            write(slice_fh,'(a,a,a)')&
                 '<DataArray type="Float64" Name="',TRIM(wnamei(iw)),&
                    '" format="ascii">'
            write(slice_fh,'(200(1pe14.6))') ((roundoff_minmax(pwC_sub&
               (jgrid)%w(ix1,ix2,iw)*normconv(iw),minvalue,maxvalue),ix1&
               =ixCmin1,ixCmax1),ix2=ixCmin2,ixCmax2)
            write(slice_fh,'(a)')'</DataArray>'
         enddo
         write(slice_fh,'(a)')'</PointData>'
         
      
      end select
      !==============================
      ! Done: celldata or pointdata?
      !==============================


      
      !==============================
      ! output Cornerpoints
      !==============================
      write(slice_fh,'(a)')'<Points>'
      write(slice_fh,'(a)')'<DataArray type="Float32" NumberOfComponents&
         ="3" format="ascii">'
 !write cell corner coordinates in a backward dimensional loop, always 3D output
      do ix2=ixCmin2,ixCmax2 
      do ix1=ixCmin1,ixCmax1 
            x_VTK(1:3)=zero;
            x_VTK(1:ndim)=pxC_sub(jgrid)%x(ix1,ix2,1:ndim)*normconv(0);
            write(slice_fh,'(3(1pe14.6))') x_VTK
      end do 
      end do 
      write(slice_fh,'(a)')'</DataArray>'
      write(slice_fh,'(a)')'</Points>'
      !==============================
      ! Done: output Cornerpoints
      !==============================


      !==============================
      ! cell Metainformation
      !==============================
      write(slice_fh,'(a)')'<Cells>'

      ! connectivity part
      write(slice_fh,'(a)')'<DataArray type="Int32" Name="connectivity" &
         format="ascii">'

       do ix2=1,nx2
       do ix1=1,nx1
      
      write(slice_fh,'(4(i7,1x))')(ix2-1)*nxC1+ix1-1, &
           (ix2-1)*nxC1+ix1,ix2*nxC1+ix1-1,ix2*nxC1+ix1
      end do
       end do
              
      write(slice_fh,'(a)')'</DataArray>'

      ! offsets data array
      write(slice_fh,'(a)')'<DataArray type="Int32" Name="offsets" format&
         ="ascii">'
      do icell=1,nc
         write(slice_fh,'(i7)') icell*(2**(3-1))
      end do
      write(slice_fh,'(a)')'</DataArray>'

      ! VTK cell type data array
      write(slice_fh,'(a)')'<DataArray type="Int32" Name="types" format&
         ="ascii">'
      ! VTK_LINE=3; VTK_PIXEL=8; VTK_VOXEL=11 -> vtk-syntax
      
       VTK_type=8 
      do icell=1,nc
         write(slice_fh,'(i2)') VTK_type
      enddo
      write(slice_fh,'(a)')'</DataArray>'
      
      write(slice_fh,'(a)')'</Cells>'
      !==============================
      ! Done: cell Metainformation
      !==============================
      write(slice_fh,'(a)')'</Piece>'

      
    end subroutine write_slice_vtk
    !=============================================================================
    subroutine put_slice_csv

      character(len=1024)           :: filename, xlabel
      character(len=79)             :: xxlabel
      logical                       :: fileopen
      character(len=10)             :: wnamei(1:nw+nwauxio),&
         xandwnamei(1:ndim+nw+nwauxio)
      character(len=1024)           :: outfilehead
      integer                       :: iw, ipe, itag
      character(len=1024)           :: line
      integer                       :: status(MPI_STATUS_SIZE)
      !-----------------------------------------------------------------------------
      if (mype==0) then
         inquire(slice_fh,opened=fileopen)
         if(.not.fileopen)then
            ! generate filename: 
            write(xlabel,"(D9.2)")xslice
            xxlabel=trim(xlabel)
            if(xslice>=zero)then
               write(xxlabel(1:1),"(a)") "+"
            endif
            write(filename,"(a,i1.1,a,i4.4,a)") TRIM(filenameout)//'_d',dir,'_x'//trim(xxlabel)//'_n',slice,'.csv'
            open(slice_fh,file=filename,status='unknown',form='formatted')
         end if
         ! get and write the header: 
         call getheadernames(wnamei,xandwnamei,outfilehead)
         line=''
         ! distinguish x,y,z and r,theta,phi etc:
         if (slice_nocartesian) then
            do iw = 1, ndim
               line = trim(line)//trim(xandwnamei(iw))//', '
            end do
         else
            line = trim(line)//trim('x')//', '
            if (ndim .gt. 1) line = trim(line)//trim('y')//', '
            if (ndim .gt. 2) line = trim(line)//trim('z')//', '
         end if
         do iw=ndim+1,ndim+nw+nwauxio-1
            line = trim(line)//trim(xandwnamei(iw))//', '
         end do
         line = trim(line)//trim(xandwnamei(ndim+nw+nwauxio))
         write(slice_fh,'(a)')trim(line)
         ! create a recv buffer using allocate, will be deallocated at the end of the routine:
         call alloc_subnode(Njgrid+1,dir,nwauxio)

         ! write to file:
         do jgrid=1, Njgrid
            call put_slice_line(jgrid,slice_fh)
         end do
      end if

      ! Also communicate the normconv array since processor zero might not have it yet:
      if (npe>1) then
         do ipe=1,npe-1
            do jgrid=1,Morton_sub_stop(ipe)-Morton_sub_start(ipe)+1
               itag=Morton_sub_start(ipe)+jgrid-1
               if (ipe == mype ) then 
                  call MPI_SEND(px_sub(jgrid)%x,1,type_subblock_x_io,0,itag,&
                     icomm,ierrmpi)
                  call MPI_SEND(pw_sub(jgrid)%w,1,type_subblock_io,0,itag,&
                     icomm,ierrmpi)
                  call MPI_SEND(normconv,nw+nwauxio+1,MPI_DOUBLE_PRECISION,0,&
                     itag,icomm,ierrmpi)
               end if
               if (mype == 0) then
                  call MPI_RECV(px_sub(Njgrid+1)%x,1,type_subblock_x_io,ipe,&
                     itag,icomm,status,ierrmpi)
                  call MPI_RECV(pw_sub(Njgrid+1)%w,1,type_subblock_io,ipe,&
                     itag,icomm,status,ierrmpi)
                  call MPI_RECV(normconv,nw+nwauxio+1,MPI_DOUBLE_PRECISION,&
                     ipe,itag,icomm,status,ierrmpi)
                  call put_slice_line(Njgrid+1,slice_fh)
               end if
            end do
         end do
      endif

      if (mype==0) then
         close(slice_fh)
         call dealloc_subnode(Njgrid+1)
      end if

    end subroutine put_slice_csv
    !=============================================================================
    subroutine put_slice_line(jout,file_handle)
      integer, intent(in) :: jout, file_handle
      ! .. local ..
      character(len=1024) ::line, data
      integer :: ix1,ix2,ix3,idir,iw
      double precision, parameter :: minvalue = 1.0d-99, maxvalue = 1.0d+99
      double precision            :: roundoff_minmax
      !-----------------------------------------------------------------------------
      
      do ix2=ixsubMlo(2),ixsubMhi(2)
         do ix1=ixsubMlo(1),ixsubMhi(1)
            
            
               ! Format the line:
               line = ''
               do idir=1,ndim
                  
                  write(data,"(es14.6)")roundoff_minmax(px_sub(jout)%x(ix1,&
                     ix2,idir),minvalue,maxvalue)
                 
                  
                  

                  line = trim(line)//trim(data)//', '
               end do
               do iw = 1,nw+nwauxio-1
                  
                  write(data,"(es14.6)")roundoff_minmax(pw_sub(jout)%w(ix1,&
                     ix2,iw)*normconv(iw),minvalue,maxvalue)
                 
                  
                  
                  line = trim(line)//trim(data)//', '
               end do
               
               write(data,"(es14.6)")roundoff_minmax(pw_sub(jout)%w(ix1,ix2,&
                  nw+nwauxio)*normconv(nw+nwauxio),minvalue,maxvalue)
              
               
               line = trim(line)//trim(data)
               write(file_handle,'(a)')trim(line)
               
            
         end do
      end do
      

    end subroutine put_slice_line
    !=============================================================================
    subroutine put_slice_dat

      integer, dimension(ngridshi) :: iorequest
      integer, dimension(MPI_STATUS_SIZE,ngridshi) :: iostatus
      integer(kind=MPI_OFFSET_KIND) :: offset
      integer :: nsubleafs
      character(len=1024) :: filename, xlabel
      character(len=79)   :: xxlabel
      integer :: amode, status(MPI_STATUS_SIZE), iwrite
      !-----------------------------------------------------------------------------

      nsubleafs=Morton_sub_stop(npe-1)
      ! generate filename
      write(xlabel,"(D9.2)")xslice
      xxlabel=trim(xlabel)
      if(xslice>=zero)then
         write(xxlabel(1:1),"(a)") "+"
      endif
      write(filename,"(a,i1.1,a,i4.4,a)") TRIM(filenameout)//'_d',dir,'_x'//trim(xxlabel)//'_n',slice,'.dat'

      if(mype==0) then
         open(unit=slice_fh,file=filename,status='replace')
         close(unit=slice_fh)
      end if

      amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
      call MPI_FILE_OPEN(icomm,filename,amode,MPI_INFO_NULL,slice_fh,ierrmpi)
      iorequest=MPI_REQUEST_NULL
      iwrite=0

      do jgrid=1,Njgrid
         iwrite=iwrite+1
         offset=int(size_subblock_io,kind=MPI_OFFSET_KIND) &
            *int(Morton_sub_start(mype)+jgrid-2,kind=MPI_OFFSET_KIND)
         call MPI_FILE_IWRITE_AT(slice_fh,offset,pw_sub(jgrid)%w,1,&
            type_subblock_io, iorequest(iwrite),ierrmpi)
      end do

      if (iwrite>0) call MPI_WAITALL(iwrite,iorequest,iostatus,ierrmpi)
      call MPI_BARRIER(icomm, ierrmpi)
      call MPI_FILE_CLOSE(slice_fh,ierrmpi)

      if (mype==0) then
         amode=ior(MPI_MODE_APPEND,MPI_MODE_WRONLY)
         call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode,MPI_INFO_NULL,&
             slice_fh,ierrmpi)

         call select_slice(dir,xslice,.true.,slice_fh,normconv)

         call MPI_FILE_WRITE(slice_fh,subsizes(2-1),1,MPI_INTEGER,status,&
            ierrmpi)
         call MPI_FILE_WRITE(slice_fh,subsizes(3-1),1,MPI_INTEGER,status,&
            ierrmpi)
         call MPI_FILE_WRITE(slice_fh,eqpar,neqpar+nspecialpar,&
             MPI_DOUBLE_PRECISION,status,ierrmpi)
         call MPI_FILE_WRITE(slice_fh,nsubleafs,1,MPI_INTEGER,status,ierrmpi)
         call MPI_FILE_WRITE(slice_fh,levmax_sub,1,MPI_INTEGER,status,ierrmpi)
         call MPI_FILE_WRITE(slice_fh,ndim-1,1,MPI_INTEGER,status,ierrmpi)
         call MPI_FILE_WRITE(slice_fh,ndir,1,MPI_INTEGER,status,ierrmpi)
         call MPI_FILE_WRITE(slice_fh,nw,1,MPI_INTEGER,status,ierrmpi)
         call MPI_FILE_WRITE(slice_fh,neqpar+nspecialpar,1,MPI_INTEGER,status,&
            ierrmpi)
         call MPI_FILE_WRITE(slice_fh,it,1,MPI_INTEGER,status,ierrmpi)
         call MPI_FILE_WRITE(slice_fh,t,1,MPI_DOUBLE_PRECISION,status,ierrmpi)

         call MPI_FILE_CLOSE(slice_fh,ierrmpi)
      end if

    end subroutine put_slice_dat
    !=============================================================================
    subroutine put_slice_zerod

      integer::  iw
      character(len=10) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
      character(len=1024) :: outfilehead
      logical, save :: opened=.false.
      character(len=1024) ::line, data
      character(len=1024) :: filename, xlabel
      character(len=79)   :: xxlabel
      integer :: amode, iwrite, status(MPI_STATUS_SIZE)
      !-----------------------------------------------------------------------------
      
    end subroutine put_slice_zerod

  end subroutine put_slice
  !=============================================================================
  subroutine select_slice(dir,xslice,writeonly,file_handle,normconv)
    use mod_forest, only: tree_node_ptr, tree_root, Morton_sub_start,&
        Morton_sub_stop
    use mod_amrvacdef
    integer, intent(in) :: dir
    double precision, intent(in) :: xslice
    integer, intent(in) :: file_handle
    logical, intent(in) :: writeonly
    double precision,dimension(0:nw+nwauxio),intent(out)       :: normconv 
    ! .. local ..
    integer :: ig1,ig2,ig3, jgrid, slice_fh, ipe, mylevmax
    integer, dimension(nlevelshi) :: igslice
    !-----------------------------------------------------------------------------
    jgrid = 0
    mylevmax = 0

    ! Find the global slice index for every level:
    call get_igslice(dir,xslice,igslice)

    ! Traverse forest to find grids indicating slice:
    
    select case(dir)
    case (1)
       ig1 = igslice(1)
       do ig3=1,ng3(1)
          do ig2=1,ng2(1)
             call traverse_slice(tree_root(ig1,ig2,ig3))
          end do
       end do
    case (2)
       ig2 = igslice(1)
       do ig3=1,ng3(1)
          do ig1=1,ng1(1)
             call traverse_slice(tree_root(ig1,ig2,ig3))
          end do
       end do
    case (3)
       ig3 = igslice(1)
       do ig2=1,ng2(1)
          do ig1=1,ng1(1)
             call traverse_slice(tree_root(ig1,ig2,ig3))
          end do
       end do
    case default
       call mpistop("slice direction not clear in select_slice")
    end select
   
    
    

    if (.not.writeonly) then
       ! Synchronize the levmax_sub for output (only rank 0 needs it): 
       levmax_sub = mylevmax
       call MPI_ALLREDUCE(MPI_IN_PLACE,levmax_sub,1,MPI_INTEGER,MPI_MAX,icomm,&
          ierrmpi)

       ! Communicate the subgrid indices according to new Morton sub-sfc:
       Morton_sub_start(:) = 1
       do ipe=0,npe-1
          call MPI_GATHER(jgrid,1,MPI_INTEGER,Morton_sub_stop,1,MPI_INTEGER,&
             ipe,icomm,ierrmpi)
       end do

       do ipe = 0, npe-2
          Morton_sub_start(ipe+1) = Morton_sub_stop(ipe)+Morton_sub_start&
             (ipe+1)
          Morton_sub_stop(ipe+1)  = Morton_sub_stop(ipe)+Morton_sub_stop&
             (ipe+1)
       end do
    end if

  contains
    !=============================================================================
    ! internal procedures
    !=============================================================================
    recursive subroutine traverse_slice(tree)
      implicit none
      type(tree_node_ptr) :: tree
      integer :: ic1,ic2,ic3
      integer, dimension(MPI_STATUS_SIZE) :: status
      !-----------------------------------------------------------------------------
      if (writeonly) then
         call MPI_FILE_WRITE(file_handle,tree%node%leaf,1,MPI_LOGICAL, status,&
            ierrmpi)
      end if

      if (tree%node%leaf) then
         if (tree%node%ipe == mype.and..not.writeonly) then
            mylevmax = max(mylevmax,tree%node%level)
            call fill_subnode(tree%node%igrid,tree%node%active,jgrid,dir,&
               xslice,normconv)
         end if
         return
      end if
      ! We are out for leaves now, continue for branches

      ! Get the correct child:
      select case (dir)
         case (1)
         ic1 = igslice(tree%node%level+1) - 2 * tree%node%ig1 + 2
         
         case (2)
         ic2 = igslice(tree%node%level+1) - 2 * tree%node%ig2 + 2
         
         case (3)
         ic3 = igslice(tree%node%level+1) - 2 * tree%node%ig3 + 2
         
      case default
         call mpistop("slice direction not clear in traverse_slice")
      end select

      ! Recursively descend into the correct child branch:
      
      select case(dir)
      case (1)
         do ic3=1,2
            do ic2=1,2
               call traverse_slice(tree%node%child(ic1,ic2,ic3))
            end do
         end do
      case (2)
         do ic3=1,2
            do ic1=1,2
               call traverse_slice(tree%node%child(ic1,ic2,ic3))
            end do
         end do
      case (3)
         do ic2=1,2
            do ic1=1,2
               call traverse_slice(tree%node%child(ic1,ic2,ic3))
            end do
         end do
      case default
         call mpistop("slice direction not clear in traverse_slice")
      end select
     
      
      

    end subroutine traverse_slice
    !=============================================================================
    ! end of internal procedures
    !=============================================================================
  end subroutine select_slice
  !=============================================================================
  subroutine fill_subnode(igrid,active,jgrid,dir,xslice,normconv)
    use mod_amrvacdef
    integer, intent(in)                                       :: igrid, dir
    integer, intent(inout)                                    :: jgrid
    logical, intent(in)                                       :: active
    double precision, intent(in)                              :: xslice
    double precision,dimension(0:nw+nwauxio),intent(out)      :: normconv 
    ! .. local ..
    double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,&
       ixMlo3-1:ixMhi3,ndim)       :: xC_TMP, xC
    double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,&
       ndim)         :: xCC_TMP, xCC
    double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,&
       ixMlo3-1:ixMhi3,nw+nwauxio) :: wC_TMP
    double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,&
       nw+nwauxio)   :: wCC_TMP
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3)                        :: x_save
    integer                :: ixslice, nwexpand, ixCmin1,ixCmin2,ixCmin3,&
       ixCmax1,ixCmax2,ixCmax3, ixCCmin1,ixCCmin2,ixCCmin3,ixCCmax1,ixCCmax2,&
       ixCCmax3
    logical                :: mask(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)&
       =.false.
    !-----------------------------------------------------------------------------
    ! increase grid-count:
    jgrid=jgrid+1
    ! Allocate subdim solution array:
    if (slice_type.ne.'dat') then
       nwexpand = nwauxio
    else
       nwexpand = 0
    end if
    call alloc_subnode(jgrid,dir,nwexpand)
    call fill_subnode_info(igrid,jgrid,dir)

    mask(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3)=.true.

    ! Now hunt for the index closest to the slice:
    
    
    
    select case (dir)
    case (1)
       ixslice = minloc(dabs(xslice-px(igrid)%x(:,ixMlo2,ixMlo3,dir)),1,&
          mask(:,ixMlo2,ixMlo3))
    case (2)
       ixslice = minloc(dabs(xslice-px(igrid)%x(ixMlo1,:,ixMlo3,dir)),1,&
          mask(ixMlo1,:,ixMlo3))
    case (3)
       ixslice = minloc(dabs(xslice-px(igrid)%x(ixMlo1,ixMlo2,:,dir)),1,&
          mask(ixMlo1,ixMlo2,:))
    case default
       call mpistop("slice direction not clear in fill_subnode")
    end select
   

    call calc_x(igrid,xC,xCC)
    ! Set the coordinate to be exactly on the slice:
    xC(:,:,:,dir)  = xslice
    xCC(:,:,:,dir) = xslice
    call calc_grid(unitslice,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,&
       normconv,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,ixCCmin1,&
       ixCCmin2,ixCCmin3,ixCCmax1,ixCCmax2,ixCCmax3,.true.,slice_nocartesian)
    ! CC stands for CellCenter
    ! C  stands for Corner
    
    ! Fill the subdimensional solution and position array:
    
    
    
    select case (dir)
    case (1)
       pw_sub(jgrid)%w(ixCCmin2:ixCCmax2,ixCCmin3:ixCCmax3,1:nw+nwexpand) = &
            wCC_TMP(ixslice,ixCCmin2:ixCCmax2,ixCCmin3:ixCCmax3,1:nw+nwexpand)
       px_sub(jgrid)%x(ixCCmin2:ixCCmax2,ixCCmin3:ixCCmax3,1:ndim) = &
            xCC_TMP(ixslice,ixCCmin2:ixCCmax2,ixCCmin3:ixCCmax3,1:ndim)
       pwC_sub(jgrid)%w(ixCmin2:ixCmax2,ixCmin3:ixCmax3,1:nw+nwexpand) = &
            wC_TMP(ixslice,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1:nw+nwexpand)
       pxC_sub(jgrid)%x(ixCmin2:ixCmax2,ixCmin3:ixCmax3,1:ndim) = &
            xC_TMP(ixslice,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1:ndim)
    case (2)
       pw_sub(jgrid)%w(ixCCmin1:ixCCmax1,ixCCmin3:ixCCmax3,1:nw+nwexpand) = &
            wCC_TMP(ixCCmin1:ixCCmax1,ixslice,ixCCmin3:ixCCmax3,1:nw+nwexpand)
       px_sub(jgrid)%x(ixCCmin1:ixCCmax1,ixCCmin3:ixCCmax3,1:ndim) = &
            xCC_TMP(ixCCmin1:ixCCmax1,ixslice,ixCCmin3:ixCCmax3,1:ndim)
       pwC_sub(jgrid)%w(ixCmin1:ixCmax1,ixCmin3:ixCmax3,1:nw+nwexpand) = &
            wC_TMP(ixCmin1:ixCmax1,ixslice,ixCmin3:ixCmax3,1:nw+nwexpand)
       pxC_sub(jgrid)%x(ixCmin1:ixCmax1,ixCmin3:ixCmax3,1:ndim) = &
            xC_TMP(ixCmin1:ixCmax1,ixslice,ixCmin3:ixCmax3,1:ndim) 
    case (3)
       pw_sub(jgrid)%w(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,1:nw+nwexpand) = &
            wCC_TMP(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,ixslice,&
               1:nw+nwexpand) 
       px_sub(jgrid)%x(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,1:ndim) = &
            xCC_TMP(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,ixslice,1:ndim)
       pwC_sub(jgrid)%w(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1:nw+nwexpand) = &
            wC_TMP(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixslice,1:nw+nwexpand) 
       pxC_sub(jgrid)%x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1:ndim) = &
            xC_TMP(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixslice,1:ndim) 
    case default
       call mpistop("slice direction not clear in fill_subnode")
    end select
   

  end subroutine fill_subnode
  !=============================================================================
  subroutine alloc_subnode(jgrid,dir,nwexpand)
    use mod_amrvacdef
    integer, intent(in) :: jgrid, dir, nwexpand
    !-----------------------------------------------------------------------------
    ! take care, what comes out is not necessarily a right handed system!
    
    
    
    select case (dir)
    case (1)
       allocate(pw_sub(jgrid)%w(ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nw+nwexpand),&
            px_sub(jgrid)%x(ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:ndim), &
            pwC_sub(jgrid)%w(ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nw+nwexpand),&
            pxC_sub(jgrid)%x(ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:ndim))
    case (2)
       allocate(pw_sub(jgrid)%w(ixGlo1:ixGhi1,ixGlo3:ixGhi3,1:nw+nwexpand),&
            px_sub(jgrid)%x(ixGlo1:ixGhi1,ixGlo3:ixGhi3,1:ndim), &
            pwC_sub(jgrid)%w(ixGlo1:ixGhi1,ixGlo3:ixGhi3,1:nw+nwexpand),&
            pxC_sub(jgrid)%x(ixGlo1:ixGhi1,ixGlo3:ixGhi3,1:ndim))
    case (3)
       allocate(pw_sub(jgrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw+nwexpand),&
            px_sub(jgrid)%x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndim), &
            pwC_sub(jgrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw+nwexpand),&
            pxC_sub(jgrid)%x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndim))
    case default
       call mpistop("slice direction not clear in alloc_subnode")
    end select
   
  end subroutine alloc_subnode
  !=============================================================================
  subroutine dealloc_subnode(jgrid)
    use mod_amrvacdef
    integer, intent(in) :: jgrid
    !-----------------------------------------------------------------------------
    if (jgrid==0) then
       call mpistop("trying to delete a non-existing grid in dealloc_subnode")
    end if

    deallocate(pw_sub(jgrid)%w,px_sub(jgrid)%x,pwC_sub(jgrid)%w,&
       pxC_sub(jgrid)%x)

    ! reset the global node info:
    node_sub(:,jgrid)=0
    rnode_sub(:,jgrid)=zero

  end subroutine dealloc_subnode
  !=============================================================================
  subroutine fill_subnode_info(igrid,jgrid,dir)
    use mod_amrvacdef
    integer, intent(in) :: igrid,jgrid,dir
    !-----------------------------------------------------------------------------

    node_sub(plevel_,jgrid)=node(plevel_,igrid)
    
    select case(dir)
    case (1)
       node_sub(pig1_,jgrid)=node(pig2_,igrid)
       node_sub(pig2_,jgrid)=node(pig3_,igrid)
       rnode_sub(rpdx1_,jgrid)=rnode(rpdx2_,igrid)
       rnode_sub(rpdx2_,jgrid)=rnode(rpdx3_,igrid)
       rnode_sub(rpxmin1_,jgrid)=rnode(rpxmin2_,igrid)
       rnode_sub(rpxmin2_,jgrid)=rnode(rpxmin3_,igrid)
    case (2)
       node_sub(pig1_,jgrid)=node(pig1_,igrid)
       node_sub(pig2_,jgrid)=node(pig3_,igrid)
       rnode_sub(rpdx1_,jgrid)=rnode(rpdx1_,igrid)
       rnode_sub(rpdx2_,jgrid)=rnode(rpdx3_,igrid)
       rnode_sub(rpxmin1_,jgrid)=rnode(rpxmin1_,igrid)
       rnode_sub(rpxmin2_,jgrid)=rnode(rpxmin3_,igrid)
    case (3)
       node_sub(pig1_,jgrid)=node(pig1_,igrid)
       node_sub(pig2_,jgrid)=node(pig2_,igrid)
       rnode_sub(rpdx1_,jgrid)=rnode(rpdx1_,igrid)
       rnode_sub(rpdx2_,jgrid)=rnode(rpdx2_,igrid)
       rnode_sub(rpxmin1_,jgrid)=rnode(rpxmin1_,igrid)
       rnode_sub(rpxmin2_,jgrid)=rnode(rpxmin2_,igrid)
    case default
       call mpistop("slice direction not clear in fill_subnode_info")
    end select
   
    
    

  end subroutine fill_subnode_info
  !=============================================================================
  subroutine get_igslice(dir,x,igslice)
    use mod_amrvacdef
    integer, intent(in) :: dir
    double precision, intent(in) :: x
    integer, dimension(nlevelshi), intent(out) :: igslice
    ! .. local ..
    integer level
    !-----------------------------------------------------------------------------

    if (x.ne.x) call mpistop("get_igslice: your slice position is NaN!")

    select case (dir)
       case (1)
       do level = 1, mxnest
          igslice(level) = int((x-xprobmin1)/dg1(level))+1
 !Gets out of domain when x==xprobmax1, not caught by put_slice, so limit:
          if (x>=xprobmax1) igslice(level) =  int((xprobmax1-xprobmin1)&
             /dg1(level))
          ! This is already caught by control in put_slice, but anyways:
          if (x<=xprobmin1) igslice(level) =  1
       end do
       case (2)
       do level = 1, mxnest
          igslice(level) = int((x-xprobmin2)/dg2(level))+1
 !Gets out of domain when x==xprobmax2, not caught by put_slice, so limit:
          if (x>=xprobmax2) igslice(level) =  int((xprobmax2-xprobmin2)&
             /dg2(level))
          ! This is already caught by control in put_slice, but anyways:
          if (x<=xprobmin2) igslice(level) =  1
       end do
       case (3)
       do level = 1, mxnest
          igslice(level) = int((x-xprobmin3)/dg3(level))+1
 !Gets out of domain when x==xprobmax3, not caught by put_slice, so limit:
          if (x>=xprobmax3) igslice(level) =  int((xprobmax3-xprobmin3)&
             /dg3(level))
          ! This is already caught by control in put_slice, but anyways:
          if (x<=xprobmin3) igslice(level) =  1
       end do
    case default
       call mpistop("slice direction not clear in get_igslice")
    end select
  end subroutine get_igslice
  !=============================================================================
end module mod_slice
!=============================================================================