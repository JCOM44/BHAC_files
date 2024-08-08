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

!=======================================================================
module mod_interpolate
  !
  ! Contains routines to find points and to interpolate to them.
  !
  implicit none

  !=======================================================================
contains
  !=======================================================================
  logical function point_in_domain(x)

    use mod_amrvacdef
    double precision, dimension(ndim), intent(in)  :: x
    integer                                        :: idim
    !----------------------------------------------------------------------

    point_in_domain = .true.

    do idim=1,ndim
       select case(idim)
          case (1)
          if (x(1) .lt. xprobmin1) then
             point_in_domain = .false.
             exit
          end if
          if (x(1) .ge. xprobmax1) then
             point_in_domain = .false.
             exit
          end if
          
          case (2)
          if (x(2) .lt. xprobmin2) then
             point_in_domain = .false.
             exit
          end if
          if (x(2) .ge. xprobmax2) then
             point_in_domain = .false.
             exit
          end if
          
          case (3)
          if (x(3) .lt. xprobmin3) then
             point_in_domain = .false.
             exit
          end if
          if (x(3) .ge. xprobmax3) then
             point_in_domain = .false.
             exit
          end if
          
       end select
    end do

  end function point_in_domain
  !====================================================================
  subroutine find_point_ipe(x,igrid_point,ipe_point)

    use mod_forest, only: tree_node_ptr, tree_root
    use mod_slice, only: get_igslice
    use mod_amrvacdef

    double precision, dimension(ndim), intent(in)   :: x
    integer, intent(out)                            :: igrid_point, ipe_point

    integer, dimension(ndir,nlevelshi)              :: ig
    integer                                         :: idim, ic(ndim)
    type(tree_node_ptr)                             :: branch
    !--------------------------------------------------------------------

    ! first check if the point is in the domain
    if (.not. point_in_domain(x)) then
       igrid_point = -1
       ipe_point   = -1
       return
    end if

    ! get the index on each level
    do idim = 1, ndim
       call get_igslice(idim,x(idim),ig(idim,:))
    end do

    ! traverse the tree until leaf is found
    branch=tree_root(ig(1,1),ig(2,1),ig(3,1))
    do while (.not.branch%node%leaf)
       ic(1)=ig(1,branch%node%level+1) - 2 * branch%node%ig1 +2
       ic(2)=ig(2,branch%node%level+1) - 2 * branch%node%ig2 +2
       ic(3)=ig(3,branch%node%level+1) - 2 * branch%node%ig3 +2
       branch%node => branch%node%child(ic(1),ic(2),ic(3))%node
    end do

    igrid_point = branch%node%igrid
    ipe_point   = branch%node%ipe

  end subroutine find_point_ipe
  !===========================================================================
  subroutine interpolate_var(igrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,gf,x,xloc,gfloc)

    use mod_amrvacdef
    integer, intent(in)                   :: igrid,ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3
    double precision, intent(in)          :: gf(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, intent(in)          :: x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision, intent(in)          :: xloc(1:ndim)
    double precision, intent(out)         :: gfloc
    integer                               :: ic1,ic2,ic3, ic11,ic12,ic13,&
        ic21,ic22,ic23, idir
    double precision                      :: xd1,xd2,xd3
    
    
    double precision                      :: c0, c1, c00, c10, c01, c11
    
    character(len=1024)                   :: line
    !---------------------------------------------------------------------------

    ! flat interpolation:
    ic1 = int((xloc(1)-rnode(rpxmin1_,igrid))/rnode(rpdx1_,igrid)) + 1 + dixB 
    ic2 = int((xloc(2)-rnode(rpxmin2_,igrid))/rnode(rpdx2_,igrid)) + 1 + dixB 
    ic3 = int((xloc(3)-rnode(rpxmin3_,igrid))/rnode(rpdx3_,igrid)) + 1 + dixB 
    !gfloc = gf(ic^D)



    if (ic1.lt.ixImin1 .or. ic1.gt.ixImax1.or.ic2.lt.ixImin2 &
       .or. ic2.gt.ixImax2.or.ic3.lt.ixImin3 .or. ic3.gt.ixImax3) then
       line = ''
       write(line,"(a)") 'Trying to flat-interpolate from out of grid!'
       write(line,"(a,a,3es14.6)") trim(line),' position: ',xloc(1:ndim)
       write(line,"(a,a,i4.3)") trim(line),' index: ', ic1,ic2,ic3
       
       if (ic1.lt.ixImin1 .or. ic1.gt.ixImax1) then
          write(line,"(a,a,i3.2)") trim(line),' Bounds exceeded in direction: &
             ',1
       else
          write(line,"(a,a,i3.2)") trim(line),' Bounds OK in direction: ',1
       end if
       write(line,"(a,a,2es14.6)") trim(line),' grid range: ', x(ixImin1,&
          ixImin2,ixImin3,1), x(ixImax1,ixImax2,ixImax3,1)
       
       
       if (ic2.lt.ixImin2 .or. ic2.gt.ixImax2) then
          write(line,"(a,a,i3.2)") trim(line),' Bounds exceeded in direction: &
             ',2
       else
          write(line,"(a,a,i3.2)") trim(line),' Bounds OK in direction: ',2
       end if
       write(line,"(a,a,2es14.6)") trim(line),' grid range: ', x(ixImin1,&
          ixImin2,ixImin3,2), x(ixImax1,ixImax2,ixImax3,2)
       
       
       if (ic3.lt.ixImin3 .or. ic3.gt.ixImax3) then
          write(line,"(a,a,i3.2)") trim(line),' Bounds exceeded in direction: &
             ',3
       else
          write(line,"(a,a,i3.2)") trim(line),' Bounds OK in direction: ',3
       end if
       write(line,"(a,a,2es14.6)") trim(line),' grid range: ', x(ixImin1,&
          ixImin2,ixImin3,3), x(ixImax1,ixImax2,ixImax3,3)
       
       call mpistop(line)
    end if
    
    ! linear interpolation:
    
    if (x(ic1,ic2,ic3,1) .lt. xloc(1)) then
       ic11 = ic1
    else
       ic11 = ic1 -1
    end if
    ic21 = ic11 + 1
    
    
    if (x(ic1,ic2,ic3,2) .lt. xloc(2)) then
       ic12 = ic2
    else
       ic12 = ic2 -1
    end if
    ic22 = ic12 + 1
    
    
    if (x(ic1,ic2,ic3,3) .lt. xloc(3)) then
       ic13 = ic3
    else
       ic13 = ic3 -1
    end if
    ic23 = ic13 + 1
    
    
    
    if (ic11.lt.ixImin1 .or. ic21.gt.ixImax1) then
       line = ''
       write(line,"(a)") 'Trying to interpolate from out of grid!'
       write(line,"(a,a,i3.2)") trim(line),' direction: ',1
       write(line,"(a,a,3es14.6)") trim(line),' position: ',xloc(1:ndim)
       write(line,"(a,a,2i4.3)") trim(line),' indices: ', ic11,ic21
       write(line,"(a,a,2es14.6)") trim(line),' grid range: ', x(ixImin1,&
          ixImin2,ixImin3,1), x(ixImax1,ixImax2,ixImax3,1)
       call mpistop(line)
    end if
    
    
    if (ic12.lt.ixImin2 .or. ic22.gt.ixImax2) then
       line = ''
       write(line,"(a)") 'Trying to interpolate from out of grid!'
       write(line,"(a,a,i3.2)") trim(line),' direction: ',2
       write(line,"(a,a,3es14.6)") trim(line),' position: ',xloc(1:ndim)
       write(line,"(a,a,2i4.3)") trim(line),' indices: ', ic12,ic22
       write(line,"(a,a,2es14.6)") trim(line),' grid range: ', x(ixImin1,&
          ixImin2,ixImin3,2), x(ixImax1,ixImax2,ixImax3,2)
       call mpistop(line)
    end if
    
    
    if (ic13.lt.ixImin3 .or. ic23.gt.ixImax3) then
       line = ''
       write(line,"(a)") 'Trying to interpolate from out of grid!'
       write(line,"(a,a,i3.2)") trim(line),' direction: ',3
       write(line,"(a,a,3es14.6)") trim(line),' position: ',xloc(1:ndim)
       write(line,"(a,a,2i4.3)") trim(line),' indices: ', ic13,ic23
       write(line,"(a,a,2es14.6)") trim(line),' grid range: ', x(ixImin1,&
          ixImin2,ixImin3,3), x(ixImax1,ixImax2,ixImax3,3)
       call mpistop(line)
    end if
    
    
    
    
    
    
    xd1 = (xloc(1)-x(ic11,ic12,ic13,1)) / (x(ic21,ic12,ic13,1) - x(ic11,ic12,&
       ic13,1))
    xd2 = (xloc(2)-x(ic11,ic12,ic13,2)) / (x(ic11,ic22,ic13,2) - x(ic11,ic12,&
       ic13,2))
    xd3 = (xloc(3)-x(ic11,ic12,ic13,3)) / (x(ic11,ic12,ic23,3) - x(ic11,ic12,&
       ic13,3))

    c00 = gf(ic11,ic12,ic13) * (1.0d0 - xd1) + gf(ic21,ic12,ic13) * xd1
    c10 = gf(ic11,ic22,ic13) * (1.0d0 - xd1) + gf(ic21,ic22,ic13) * xd1
    c01 = gf(ic11,ic12,ic23) * (1.0d0 - xd1) + gf(ic21,ic12,ic23) * xd1
    c11 = gf(ic11,ic22,ic23) * (1.0d0 - xd1) + gf(ic21,ic22,ic23) * xd1

    c0  = c00 * (1.0d0 - xd2) + c10 * xd2
    c1  = c01 * (1.0d0 - xd2) + c11 * xd2

    gfloc = c0 * (1.0d0 - xd3) + c1 * xd3
    

  end subroutine interpolate_var
  !=======================================================================
end module mod_interpolate
!=======================================================================
