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

!> Module to couple the octree-mg library to BHAC. This file uses the VACPP
!> preprocessor, but its use is kept to a minimum.


module mod_multigrid_coupling
  
  
  use m_octree_mg_3d
 

  implicit none
  public

  !> Data structure containing the multigrid tree.
  type(mg_t) :: mg

  !> If defined, this routine is called after a new multigrid tree is
  !> constructed.
  procedure(after_new_tree), pointer :: mg_after_new_tree => null()

  integer :: block_nx1,block_nx2,block_nx3

  interface
     subroutine after_new_tree()
     end subroutine after_new_tree
  end interface

contains

  !> Setup multigrid for usage
  subroutine mg_setup_multigrid()
    use mod_amrvacdef

    block_nx1=ixMhi1-ixMlo1+1;block_nx2=ixMhi2-ixMlo2+1
    block_nx3=ixMhi3-ixMlo3+1;

    if (ndim == 1) &
         error stop "Multigrid not available in 1D"

    if (ndim /= mg_ndim) &
         error stop "Multigrid module was compiled for different ndim"

    select case (typeaxial)
    case ("slab")
       if (ndim == 1) error stop "Multigrid only support 2D, 3D"
    case ("cylindrical")
       if (ndim == 3) error stop "Multigrid does not support cylindrical 3D"
       mg%geometry_type = mg_cylindrical
    case default
       error stop "Multigrid does not support your geometry"
    end select

    if (any([ block_nx1,block_nx2,block_nx3 ] /= block_nx1)) &
         error stop "Multigrid requires all block_nx to be equal"

    ! Use Gauss-Seidel red-black smoother, which should be more robust
    mg%smoother_type = mg_smoother_gsrb

    call mg_comm_init(mg)
    call mg_set_methods(mg)
    call mg_tree_from_bhac(mg)
  end subroutine mg_setup_multigrid

 !> Set multigrid boundary conditions for the solution according to variable iw
  subroutine mg_copy_boundary_conditions(mg, iw)
    use mod_amrvacdef
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: iw
    integer                   :: n

    do n = 1, mg_num_neighbors
       select case (typeB(iw, n))
       case ('symm')
          mg%bc(n, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case ('asymm')
          mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case ('cont')
          mg%bc(n, mg_iphi)%bc_type = mg_bc_continuous
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp ! Not needed
       case ('periodic')
          ! Nothing to do here
       case default
          print *, "Not a standard: ", trim(typeB(iw, n))
          error stop "You have to set a user-defined boundary method"
       end select
    end do
  end subroutine mg_copy_boundary_conditions

  !> If the grid has changed, rebuild the full multigrid tree
  subroutine mg_update_refinement(n_coarsen, n_refine)
    use mod_amrvacdef
    integer, intent(in) :: n_coarsen
    integer, intent(in) :: n_refine

    ! Don't build multigrid tree while doing initial refinement
    if (.not. time_advance) return

    if (.not. mg%is_allocated) then
       call mg_tree_from_bhac(mg)
    else if (n_coarsen + n_refine > 0) then
       call mg_deallocate_storage(mg)
       call mg_tree_from_bhac(mg)
    end if
  end subroutine mg_update_refinement

  !> Copy a variable to the multigrid tree, including a layer of ghost cells
  subroutine mg_copy_to_tree(iw_from, iw_to, restrict, restrict_gc)
    use mod_forest
    use mod_amrvacdef
    integer, intent(in)      :: iw_from !< Variable to use as right-hand side
    integer, intent(in)      :: iw_to      !< Copy to this variable
    logical, intent(in)      :: restrict !< Restrict variable on multigrid tree
    logical, intent(in)      :: restrict_gc !< Fill ghost cells after restrict
    integer                  :: iigrid, igrid, id
    integer                  :: nc, lvl
    type(tree_node), pointer :: pnode

    if (.not. mg%is_allocated) &
         error stop "mg_copy_to_tree: tree not allocated yet"

    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)

       ! Include one layer of ghost cells on grid leaves
       
       
       mg%boxes(id)%cc(0:nc+1, 0:nc+1, 0:nc+1, iw_to) = &
            pw(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, &
            ixMlo3-1:ixMhi3+1, iw_from)
      
    end do

    if (restrict) then
       call mg_restrict(mg, iw_to)
       if (restrict_gc) call mg_fill_ghost_cells(mg, iw_to)
    end if

  end subroutine mg_copy_to_tree

  !> Copy a variable from the multigrid tree
  subroutine mg_copy_from_tree(iw_from, iw_to)
    use mod_forest
    use mod_amrvacdef
    integer, intent(in)      :: iw_from !< Variable to use as right-hand side
    integer, intent(in)      :: iw_to   !< Copy to this variable
    integer                  :: iigrid, igrid, id
    integer                  :: nc, lvl
    type(tree_node), pointer :: pnode

    if (.not. mg%is_allocated) &
         error stop "mg_copy_from_tree: tree not allocated yet"

    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)



       pw(igrid)%w(ixMlo1:ixMhi1, ixMlo2:ixMhi2, ixMlo3:ixMhi3, iw_to) = &
            mg%boxes(id)%cc(1:nc, 1:nc, 1:nc, iw_from)

    end do
  end subroutine mg_copy_from_tree

 !> Copy from multigrid tree with one layer of ghost cells. Corner ghost cells
  !> are not used/set.
  subroutine mg_copy_from_tree_gc(iw_from, iw_to)
    use mod_forest
    use mod_amrvacdef
    integer, intent(in)      :: iw_from !< Variable to use as right-hand side
    integer, intent(in)      :: iw_to   !< Copy to this variable
    integer                  :: iigrid, igrid, id
    integer                  :: nc, lvl
    type(tree_node), pointer :: pnode

    if (.not. mg%is_allocated) &
         error stop "mg_copy_from_tree_gc: tree not allocated yet"

    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)



       pw(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, &
            ixMlo3-1:ixMhi3+1, iw_to) = &
            mg%boxes(id)%cc(0:nc+1, 0:nc+1, 0:nc+1, iw_from)

    end do
  end subroutine mg_copy_from_tree_gc

  !> Generate a multigrid tree that includes the bhac tree, but also contains
  !> coarser grid levels. A number of checks has already been performed in
  !> mg_setup_multigrid, so we don't repeat these checks here.
  subroutine mg_tree_from_bhac(mg)
    use mod_forest
    use mod_amrvacdef
    type(mg_t), intent(inout)        :: mg
    integer                          :: i, n, id, ix(ndim)
    integer                          :: n_boxes_total, i_c, c_id, c_ix(ndim)
    integer                          :: min_lvl, lvl
    integer                          :: nb, nb_ix, nb_dim
    integer                          :: n_finer
    integer                          :: nxlone(ndim)
    type(tree_node), pointer         :: pnode, pnode_ch
    type(tree_node_ptr), allocatable :: id_to_node(:)
    real(dp)                         :: dr_coarse

    ! Estimate number of finer blocks
    n_finer = nparents+nleafs

    if (any([ block_nx1,block_nx2,block_nx3 ] /= block_nx1)) then
       call mpistop("Multigrid required square blocks")
    end if

    nxlone = [ ng1(1),ng2(1),ng3(1) ] * block_nx1

    call mg_build_rectangle(mg, nxlone, block_nx1, &
         dx(:,1), [ xprobmin1,xprobmin2,xprobmin3 ], periodB, n_finer)

    mg%highest_lvl = levmax
    n_boxes_total = mg%n_boxes + n_finer

    ! To link the two trees
    allocate(id_to_node(n_boxes_total))

    ! Link base level
    do i = 1, size(mg%lvls(1)%ids)
       id = mg%lvls(1)%ids(i)
       ix = mg%boxes(id)%ix

       pnode               => tree_root(ix(1),ix(2),ix(3))%node
       pnode%id            =  id
       id_to_node(id)%node => pnode
       mg%boxes(id)%rank   =  pnode%ipe
    end do

    ! Add refinement
    do lvl = 1, mg%highest_lvl
       do i = 1, size(mg%lvls(lvl)%ids)
          id = mg%lvls(lvl)%ids(i)
          pnode => id_to_node(id)%node

          if (.not. pnode%leaf) then
             call mg_add_children(mg, id)

             do i_c = 1, mg_num_children
                c_id = mg%boxes(id)%children(i_c)
                c_ix = mg_child_dix(:, i_c) + 1
                pnode_ch => pnode%child(c_ix(1),c_ix(2),c_ix(3))%node
                id_to_node(c_id)%node => pnode_ch
                pnode_ch%id = c_id
                mg%boxes(c_id)%rank = pnode_ch%ipe
             end do
          end if
       end do

       call mg_set_leaves_parents(mg%boxes, mg%lvls(lvl))

       if (lvl < mg%highest_lvl) then
          call mg_set_next_level_ids(mg, lvl)
          call mg_set_neighbors_lvl(mg, lvl+1)
       end if
    end do

    ! Store boxes with refinement boundaries (from the coarse side)
    do lvl = 1, mg%highest_lvl
       call mg_set_refinement_boundaries(mg%boxes, mg%lvls(lvl))
    end do

    ! Assign boxes to MPI processes
    call mg_load_balance_parents(mg)

    ! Allocate storage for boxes owned by this process
    call mg_allocate_storage(mg)

    if (associated(mg_after_new_tree)) then
       call mg_after_new_tree()
    end if

  end subroutine mg_tree_from_bhac

  subroutine clean_divb_multigrid()
    use mod_forest
    use mod_amrvacdef
    integer                     :: iigrid, igrid, id
    integer                     :: n, nc, lvl, ixmin1,ixmin2,ixmin3,ixmax1,&
       ixmax2,ixmax3, idim
    type(tree_node), pointer    :: pnode
    double precision            :: tmp(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3), grad(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3, ndim)
    double precision, parameter :: residual_reduction = 1d-10
    integer, parameter          :: max_its            = 50
    double precision            :: residual_it(max_its)
    integer                     :: i, j, i0, j0
    real(dp)                    :: phi_grad, max_divb
    
    integer                     :: k, k0
   

    mg%operator_type = mg_laplacian

    ! TODO (?) Set boundary conditions
    mg%bc(:, mg_iphi)%bc_type = mg_bc_dirichlet
    mg%bc(:, mg_iphi)%bc_value = 0.0_dp

    ixmin1=ixMlo1-1;ixmin2=ixMlo2-1;ixmin3=ixMlo3-1;ixmax1=ixMhi1+1
    ixmax2=ixMhi2+1;ixmax3=ixMhi3+1;
    max_divb = 0.0d0

    ! Store divergence of B as right-hand side
    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)

       ! Geometry subroutines expect this to be set
       dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
       dxlevel(3)=rnode(rpdx3_,igrid);

       
       call div_staggered(ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,ps(igrid),&
          tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3))
       
       mg%boxes(id)%cc(1:nc,1:nc,1:nc, mg_irhs) = tmp(ixMlo1:ixMhi1,&
          ixMlo2:ixMhi2,ixMlo3:ixMhi3)
       max_divb = max(max_divb, maxval(abs(tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
          ixMlo3:ixMhi3))))
    end do

    call MPI_ALLREDUCE(MPI_IN_PLACE, max_divb, 1, MPI_DOUBLE_PRECISION, &
         MPI_MAX, icomm, ierrmpi)

    ! Solve laplacian(phi) = divB
    if (mype == 0) print *, "Performing multigrid divB cleaning"
    if (mype == 0) print *, "iteration vs residual"

    do n = 1, max_its
       call mg_fas_fmg(mg, n>1, max_res=residual_it(n))

       ! V-cycles are cheaper but converge less quickly
       ! call mg_fas_vcycle(mg, max_res=residual_it(n))

       if (mype == 0) write(*, "(I4,E11.3)") n, residual_it(n)
       if (residual_it(n) < residual_reduction * max_divb) exit
    end do

    if (mype == 0 .and. n > max_its) then
       print *, "divb_multigrid warning: not fully converged"
       print *, "current amplitude of divb: ", residual_it(max_its)
       print *, "multigrid smallest grid: ", &
            mg%domain_size_lvl(:, mg%lowest_lvl)
       print *, "note: smallest grid ideally has <= 8 cells"
       print *, "multigrid dx/dy/dz ratio: ", mg%dr(:, 1)/mg%dr(1, 1)
       print *, "note: dx/dy/dz should be similar"
    end if

    ! Correct the magnetic field
    do iigrid = 1, igridstail
       igrid = igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       nc    =  mg%box_size

       ! Geometry subroutines expect this to be set
       dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
       dxlevel(3)=rnode(rpdx3_,igrid);

       ! Compute the gradient of phi
       tmp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3) = mg%boxes(id)%cc(:,:,:,&
           mg_iphi)

       
       

       
       ! first dimension
       do k = 1, nc
          do j = 1, nc
             do i = 0, nc
                phi_grad = (mg%boxes(id)%cc(i+1, j, k, mg_iphi) - &
                     mg%boxes(id)%cc(i, j, k, mg_iphi)) / dxlevel(1)
                i0 = i + dixB
                j0 = j + dixB
                k0 = k + dixB
                pws(igrid)%w(i0, j0, k0, bs1_) = &
                     pws(igrid)%w(i0, j0, k0, bs1_) - phi_grad
             end do
          end do
       end do

       ! second dimension
       do k = 1, nc
          do j = 0, nc
             do i = 1, nc
                phi_grad = (mg%boxes(id)%cc(i, j+1, k, mg_iphi) - &
                     mg%boxes(id)%cc(i, j, k, mg_iphi)) / dxlevel(2)
                i0 = i + dixB
                j0 = j + dixB
                k0 = k + dixB
                pws(igrid)%w(i0, j0, k0, bs2_) = &
                     pws(igrid)%w(i0, j0, k0, bs2_) - phi_grad
             end do
          end do
       end do

       ! third dimension
       do k = 0, nc
          do j = 1, nc
             do i = 1, nc
                phi_grad = (mg%boxes(id)%cc(i, j, k+1, mg_iphi) - &
                     mg%boxes(id)%cc(i, j, k, mg_iphi)) / dxlevel(3)
                i0 = i + dixB
                j0 = j + dixB
                k0 = k + dixB
                pws(igrid)%w(i0, j0, k0, bs3_) = &
                     pws(igrid)%w(i0, j0, k0, bs3_) - phi_grad
             end do
          end do
       end do

       call faces2centers(ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,&
           ps(igrid))
   

    end do
    
  end subroutine clean_divb_multigrid

end module mod_multigrid_coupling

