submodule (eu_model) forwardstep
implicit none

contains
! Double Maach Reflection
module subroutine initial_model_eu_fs(funit)
    implicit none
    integer,intent(in) :: funit
  
    type(eu_data_list_t), pointer :: flwd_node
    type(eu_data_t), pointer      :: flwd
    type(eu_info_data_t),pointer  :: info
    type(rvectorfield_t),pointer  :: prim
    type(rvectorfield_t),pointer  :: cons 
  
    class(block_list_t),pointer   :: blk_node 
    class(block_t),pointer        :: blk,blk1,blk2
    class(block_data_t),pointer   :: d_blk,d_blk1,d_blk2

    class(bndy_list_t),pointer    :: bndy_node 
    class(bndy_t),pointer         :: bndy 
    class(bndy_data_t),pointer    :: d_bndy 
  
    integer :: iflwd
    integer :: i,j,iblk,bc_num
    integer :: i0,j0
    integer :: ilo(2),ihi(2),jlo(2),jhi(2)
    double precision :: xc0,yc0,rc0
    double precision :: rho0,u0,v0,p0,ei0
    
    print*,nn,ibf_plane,jbf_plane
	! Create Front step problem blocks
	blk_node => add_blk_to_mesh()
	blk1     => blk_node%blk
	call blk1%create_block(3*nn,5*nn,ibf_plane,jbf_plane,1)
    d_blk1   => blk1%data
    
	blk_node => add_blk_to_mesh()
	blk2     => blk_node%blk
	call blk2%create_block(12*nn,4*nn,ibf_plane,jbf_plane,2)
	d_blk2   => blk2%data
	
	! Step-2: Create and setup boundary conditions for each block
    ! Create boundary condition for block 1
    ! Boundary 1 of block 1
    ilo       = (/1, 1/)
    ihi       = (/d_blk1%ibf, 1/)
    jlo       = (/d_blk1%jbf + 1, 1/)
    jhi       = (/d_blk1%jbf + d_blk1%jsize, 1/)
    bndy_node => blk1%add_bndy(left_bndy_id, ilo, ihi, jlo, jhi)	
    d_bndy    => bndy_node%bndy%data
    call create(d_bndy%bc_infos,2)
    d_bndy%bc_infos%p(1) = inlt11_fs_idx
    d_bndy%bc_infos%p(2) = 1

    ! Boundary 2 of block 1
    ilo       = (/d_blk1%ibf + 1, 1/)
    ihi       = (/d_blk1%ibf + d_blk1%isize, 1/)
    jlo       = (/1, 1/)
    jhi       = (/d_blk1%jbf, 1/)
    bndy_node => blk1%add_bndy(bottom_bndy_id, ilo, ihi, jlo, jhi)	
    d_bndy    => bndy_node%bndy%data
    call create(d_bndy%bc_infos,2)
    d_bndy%bc_infos%p(1) = rflt12_fs_idx
    d_bndy%bc_infos%p(2) = 1
    d_bndy%vector_num    = 1
    allocate(d_bndy%vector(1))
    call create(d_bndy%vector(1), 4)
    d_bndy%vector(1)%p(:) = 1.d0
    d_bndy%vector(1)%p(3) =-1.d0

    ! Boundary 3 of block 1
    ilo       = (/d_blk1%ibf + d_blk1%isize + 1, 1/)
    ihi       = (/d_blk1%ibf + d_blk1%isize + d_blk1%ibf, 1/)
    jlo       = (/d_blk1%jbf + 1, 1/)
    jhi       = (/d_blk1%jbf + nn, 1/)
    bndy_node => blk1%add_bndy(right_bndy_id, ilo, ihi, jlo, jhi)	
    d_bndy    => bndy_node%bndy%data
    call create(d_bndy%bc_infos,2)
    d_bndy%bc_infos%p(1) = rflt13_fs_idx
    d_bndy%bc_infos%p(2) = 1
    d_bndy%vector_num    = 1
    allocate(d_bndy%vector(1))
    call create(d_bndy%vector(1), 4)
    d_bndy%vector(1)%p(:) = 1.d0
    d_bndy%vector(1)%p(2) =-1.d0

    ! Boundary 4 of block 1
    ilo       = (/d_blk1%ibf + d_blk1%isize + 1, d_blk2%ibf + 1/)
    ihi       = (/d_blk1%ibf + d_blk1%isize + d_blk1%ibf, d_blk2%ibf + d_blk2%ibf/)
    jlo       = (/d_blk1%jbf + nn + 1, d_blk2%jbf + 1/)
    jhi       = (/d_blk1%jbf + d_blk1%jsize, d_blk2%jbf + d_blk2%jsize/)
    bndy_node => blk1%add_bndy(right_bndy_id, ilo, ihi, jlo, jhi)	
    d_bndy    => bndy_node%bndy%data
    call create(d_bndy%bc_infos,2)
    d_bndy%bc_infos%p(1) = bc_idx_intf
    d_bndy%bc_infos%p(2) = 2

    ! Boundary 5 of block 1
    ilo       = (/d_blk1%ibf + 1, 1/)
    ihi       = (/d_blk1%ibf + d_blk1%isize, 1/)
    jlo       = (/d_blk1%jbf + d_blk1%jsize + 1, 1/)
    jhi       = (/d_blk1%jbf + d_blk1%jsize + d_blk1%jbf, 1/)
    bndy_node => blk1%add_bndy(top_bndy_id, ilo, ihi, jlo, jhi)	
    d_bndy    => bndy_node%bndy%data
    call create(d_bndy%bc_infos,2)
    d_bndy%bc_infos%p(1) = rflt15_fs_idx
    d_bndy%bc_infos%p(2) = 1
    d_bndy%vector_num    = 1
    allocate(d_bndy%vector(1))
    call create(d_bndy%vector(1), 4)
    d_bndy%vector(1)%p(:) = 1.d0
    d_bndy%vector(1)%p(3) =-1.d0

	
	! Create boundary condition for block 2
    ! Boundary 1 of block 2
    ilo       = (/1, d_blk1%isize + 1/)
    ihi       = (/d_blk2%ibf, d_blk1%isize + d_blk1%ibf/)
    jlo       = (/d_blk2%jbf + 1, d_blk1%jbf + nn + 1/)
    jhi       = (/d_blk2%jbf + d_blk2%jsize, d_blk1%jbf + d_blk1%jsize/)
    bndy_node => blk2%add_bndy(left_bndy_id, ilo, ihi, jlo, jhi)	
    d_bndy    => bndy_node%bndy%data
    call create(d_bndy%bc_infos,2)
    d_bndy%bc_infos%p(1) = bc_idx_intf
    d_bndy%bc_infos%p(2) = 1

    ! Boundary 2 of block 2
    ilo       = (/d_blk2%ibf + 1, 1/)
    ihi       = (/d_blk2%ibf + d_blk2%isize, 1/)
    jlo       = (/1, 1/)
    jhi       = (/d_blk2%jbf, 1/)
    bndy_node => blk2%add_bndy(bottom_bndy_id, ilo, ihi, jlo, jhi)	
    d_bndy    => bndy_node%bndy%data
    call create(d_bndy%bc_infos,2)
    d_bndy%bc_infos%p(1) = rflt22_fs_idx
    d_bndy%bc_infos%p(2) = 1
    d_bndy%vector_num    = 1
    allocate(d_bndy%vector(1))
    call create(d_bndy%vector(1), 4)
    d_bndy%vector(1)%p(:) = 1.d0
    d_bndy%vector(1)%p(3) =-1.d0

    ! Boundary 3 of block 2
    ilo       = (/d_blk2%ibf + d_blk2%isize + 1, 1/)
    ihi       = (/d_blk2%ibf + d_blk2%isize + d_blk2%ibf, 1/)
    jlo       = (/d_blk2%jbf + 1, 1/)
    jhi       = (/d_blk2%jbf + d_blk2%jsize, 1/)
    bndy_node => blk2%add_bndy(right_bndy_id, ilo, ihi, jlo, jhi)	
    d_bndy    => bndy_node%bndy%data
    call create(d_bndy%bc_infos,2)
    d_bndy%bc_infos%p(1) = otlt23_fs_idx
    d_bndy%bc_infos%p(2) = axis_xi
    d_bndy%vector_num    = 1
    allocate(d_bndy%vector(1))
    call create(d_bndy%vector(1), 4)
    d_bndy%vector(1)%p(:) = 0.d0

    ! Boundary 4 of block 2
    ilo       = (/d_blk2%ibf + 1, 1/)
    ihi       = (/d_blk2%ibf + d_blk2%isize, 1/)
    jlo       = (/d_blk2%jbf + d_blk2%jsize + 1, 1/)
    jhi       = (/d_blk2%jbf + d_blk2%jsize + d_blk2%jbf, 1/)
    bndy_node => blk2%add_bndy(top_bndy_id, ilo, ihi, jlo, jhi)	
    d_bndy    => bndy_node%bndy%data
    call create(d_bndy%bc_infos,2)
    d_bndy%bc_infos%p(1) = rflt24_fs_idx
    d_bndy%bc_infos%p(2) = 1
    d_bndy%vector_num    = 1
    allocate(d_bndy%vector(1))
    call create(d_bndy%vector(1), 4)
    d_bndy%vector(1)%p(:) = 1.d0
    d_bndy%vector(1)%p(3) =-1.d0

	! Step-3: Initial block geometric data
	call blk1%initial_block(xa0,ya0,xal,yal)
    call blk2%initial_block(xb0,yb0,xbl,ybl)  
  
    call create_work_que()

    rho0 = gamma
    u0   = 3.d0
    v0   = 0.d0
    p0   = 1.d0
    ei0  = p0/(rho0*(gamma-1.d0))

    do iflwd = 1, flwd_count
      flwd       => flwd_que(iflwd)%flwd
      blk        => flwd%blk
      info       => flwd%info
      cons       => flwd%cons(max_stage)
      info%gamma = gamma
      do j=1,info%jsize
        do i=1,info%isize
          i0 = i + info%ibf
          j0 = j + info%jbf          
  
          cons%p(1,i0,j0) = rho0
          cons%p(2,i0,j0) = rho0*u0
          cons%p(3,i0,j0) = rho0*v0
          cons%p(4,i0,j0) = rho0*(ei0+0.5d0*(u0**2+v0**2))
        enddo
      enddo
      info%eig_xi  = 0.d0
      info%eig_eta = 0.d0
      info%time    = 0.d0
    enddo
  end subroutine initial_model_eu_fs
end submodule forwardstep