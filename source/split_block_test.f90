subroutine split_block_test
	use mesh
	implicit none
	class(block_list_t),pointer :: blk_node 
	class(block_iterator_t),dimension(:),allocatable :: blk_iter
	class(block_t),pointer :: blk, tile
	class(block_data_t),pointer :: d_blk, d_tile
	class(bndy_list_t),pointer :: bndy_node 
	class(bndy_t),pointer :: bndy 
	class(bndy_data_t),pointer :: d_bndy 
	class(block_t),pointer :: blk1, blk2
	class(block_data_t),pointer :: d_blk1, d_blk2
	integer :: nn = 40
	integer :: itile,jtile  
	integer :: tile_num = 0
	integer :: tile_size(2) = (/64,64/)

	write(*,*) "********************************************************************************"
	write(*,*) 
	write(*,*) "						       THE DENSE CODE                                   "
	write(*,*) 
	write(*,*) "                              Two Dimensional                                   "
	write(*,*) "                   Multiple-phase Multi-Material Simulation                     "
	write(*,*) "                                   Program                                      "
	write(*,*) 
	write(*,*) "                                 By Jun Peng                                    "
	write(*,*) "                           Email:pengjun@imech.ac.cn                            "
	write(*,*) 
	write(*,*) "********************************************************************************"
	! Steps to add a block to mesh: 3 steps
	allocate(blk_iter(2))
	! Step-1: add a node to block list and create block data
	! Create Front step problem block 1 
	blk_node        => add_blk_to_mesh()
	blk_iter(1)%blk => blk_node%blk
	call blk_iter(1)%blk%create_block(3*nn,5*nn,3,3,1)

	! Create Front step problem block 2
	blk_node        => add_blk_to_mesh()
	blk_iter(2)%blk => blk_node%blk
	call blk_iter(2)%blk%create_block(12*nn,4*nn,3,3,2)
	
	! Step-2: Create and setup boundary conditions for each block
	blk1           => blk_iter(1)%blk
	blk2           => blk_iter(2)%blk
	d_blk1         => blk1%data
	d_blk2         => blk2%data

	! Create boundary condition for block 1
	bndy_node     => blk1%add_bndy()	
	bndy          => bndy_node%bndy
	bndy%group_id = 1
	d_bndy        => bndy%data
	d_bndy%ilo(1) = 1
	d_bndy%ihi(1) = d_blk1%ibf
	d_bndy%jlo(1) = d_blk1%jbf + 1
	d_bndy%jhi(1) = d_blk1%jbf + d_blk1%jsize
	call create(d_bndy%bc_infos,2)
	d_bndy%bc_infos%p(1) = bc_idx_inlt
	d_bndy%bc_infos%p(2) = 0

	bndy_node     => blk1%add_bndy()	
	bndy          => bndy_node%bndy
	bndy%group_id = 2
	d_bndy        => bndy%data
	d_bndy%ilo(1) = d_blk1%ibf + 1
	d_bndy%ihi(1) = d_blk1%ibf + d_blk1%isize
	d_bndy%jlo(1) = d_blk1%jbf
	d_bndy%jhi(1) = 1
	call create(d_bndy%bc_infos,2)
	d_bndy%bc_infos%p(1) = bc_idx_rflt
	d_bndy%bc_infos%p(2) = 0

	bndy_node     => blk1%add_bndy()	
	bndy          => bndy_node%bndy
	bndy%group_id = 3
	d_bndy        => bndy%data
	d_bndy%ilo(1) = d_blk1%ibf + d_blk1%isize + 1
	d_bndy%ihi(1) = d_blk1%ibf + d_blk1%isize + d_blk1%ibf
	d_bndy%jlo(1) = d_blk1%jbf + 1
	d_bndy%jhi(1) = d_blk1%jbf + nn
	call create(d_bndy%bc_infos,2)
	d_bndy%bc_infos%p(1) = bc_idx_rflt
	d_bndy%bc_infos%p(2) = 0

	bndy_node     => blk1%add_bndy()	
	bndy          => bndy_node%bndy
	bndy%group_id = 3
	d_bndy        => bndy%data
	d_bndy%ilo(1) = d_blk1%ibf + d_blk1%isize + 1
	d_bndy%ihi(1) = d_blk1%ibf + d_blk1%isize + d_blk1%ibf
	d_bndy%jlo(1) = d_blk1%jbf + nn + 1
	d_bndy%jhi(1) = d_blk1%jbf + d_blk1%jsize
	d_bndy%ilo(2) = d_blk2%ibf + 1
	d_bndy%ihi(2) = d_blk2%ibf + d_blk2%ibf
	d_bndy%jlo(2) = d_blk2%jbf + 1 
	d_bndy%jhi(2) = d_blk2%jbf + d_blk2%jsize
	call create(d_bndy%bc_infos,2)
	d_bndy%bc_infos%p(1) = bc_idx_intf
	d_bndy%bc_infos%p(2) = 2

	bndy_node     => blk1%add_bndy()	
	bndy          => bndy_node%bndy
	bndy%group_id = 4
	d_bndy        => bndy%data
	d_bndy%ilo(1) = d_blk1%ibf + 1
	d_bndy%ihi(1) = d_blk1%ibf + d_blk1%isize
	d_bndy%jlo(1) = d_blk1%jbf + d_blk1%jsize + 1
	d_bndy%jhi(1) = d_blk1%jbf + d_blk1%jsize + d_blk1%jbf
	call create(d_bndy%bc_infos,2)
	d_bndy%bc_infos%p(1) = bc_idx_rflt
	d_bndy%bc_infos%p(2) = 0
	
	! Create boundary condition for block 2
	bndy_node     => blk2%add_bndy()	
	bndy          => bndy_node%bndy
	bndy%group_id = 1
	d_bndy        => bndy%data
	d_bndy%ilo(1) = 1
	d_bndy%ihi(1) = d_blk2%ibf
	d_bndy%jlo(1) = d_blk2%jbf + 1
	d_bndy%jhi(1) = d_blk2%jbf + d_blk2%jsize
	d_bndy%ilo(2) = d_blk1%isize + 1
	d_bndy%ihi(2) = d_blk1%isize + d_blk1%ibf
	d_bndy%jlo(2) = d_blk1%jbf + nn + 1
	d_bndy%jhi(2) = d_blk1%jbf + d_blk1%jsize
	call create(d_bndy%bc_infos,2)
	d_bndy%bc_infos%p(1) = bc_idx_intf
	d_bndy%bc_infos%p(2) = 1

	bndy_node     => blk2%add_bndy()	
	bndy          => bndy_node%bndy
	bndy%group_id = 2
	d_bndy        => bndy%data
	d_bndy%ilo(1) = d_blk2%ibf + 1
	d_bndy%ihi(1) = d_blk2%ibf + d_blk2%isize
	d_bndy%jlo(1) = d_blk2%jbf
	d_bndy%jhi(1) = 1
	call create(d_bndy%bc_infos,2)
	d_bndy%bc_infos%p(1) = bc_idx_rflt
	d_bndy%bc_infos%p(2) = 0

	bndy_node     => blk2%add_bndy()	
	bndy          => bndy_node%bndy
	bndy%group_id = 3
	d_bndy        => bndy%data
	d_bndy%ilo(1) = d_blk2%ibf + d_blk2%isize + 1
	d_bndy%ihi(1) = d_blk2%ibf + d_blk2%isize + d_blk2%ibf
	d_bndy%jlo(1) = d_blk2%jbf + 1
	d_bndy%jhi(1) = d_blk2%jbf + d_blk2%jsize
	call create(d_bndy%bc_infos,2)
	d_bndy%bc_infos%p(1) = bc_idx_otlt
	d_bndy%bc_infos%p(2) = axis_xi

	bndy_node     => blk2%add_bndy()	
	bndy          => bndy_node%bndy
	bndy%group_id = 4
	d_bndy        => bndy%data
	d_bndy%ilo(1) = d_blk2%ibf + 1
	d_bndy%ihi(1) = d_blk2%ibf + d_blk2%isize
	d_bndy%jlo(1) = d_blk2%jbf + d_blk2%jsize + 1
	d_bndy%jhi(1) = d_blk2%jbf + d_blk2%jsize + d_blk2%jbf
	call create(d_bndy%bc_infos,2)
	d_bndy%bc_infos%p(1) = bc_idx_rflt
	d_bndy%bc_infos%p(2) = 0

	! Step-3: Initial block geometric data
	call blk1%initial_block(0.d0,0.d0,0.6d0,1.d0)
	call blk2%initial_block(0.6d0,0.2d0,2.4d0,0.8d0)

	print*,"Num of Boundaries of block 1:",blk1%bndy_num
	bndy_node => blk1%bndy_list
	do while(associated(bndy_node))
		bndy   => bndy_node%bndy
		d_bndy => bndy%data
		print*,"  bndy:",bndy%idx, bndy%group_id
		print*,"  --- (ilo,ihi)",d_bndy%ilo(1),d_bndy%ihi(1)
		print*,"  --- (jlo,jhi)",d_bndy%jlo(1),d_bndy%jhi(1)
		print*,"  --- (ilo_d,ihi_d)",d_bndy%ilo(2),d_bndy%ihi(2)
		print*,"  --- (jlo_d,jhi_d)",d_bndy%jlo(2),d_bndy%jhi(2)
		print*,"  --- BC_info",d_bndy%bc_infos%p
		print*,""
		bndy_node => bndy_node%next
	enddo
	print*,""

	print*,"Num of Boundaries of block 2:",blk2%bndy_num
	bndy_node => blk2%bndy_list
	do while(associated(bndy_node))
		bndy   => bndy_node%bndy
		d_bndy => bndy%data
		print*,"  bndy:",bndy%idx, bndy%group_id
		print*,"  --- (ilo,ihi)",d_bndy%ilo(1),d_bndy%ihi(1)
		print*,"  --- (jlo,jhi)",d_bndy%jlo(1),d_bndy%jhi(1)
		print*,"  --- (ilo_d,ihi_d)",d_bndy%ilo(2),d_bndy%ihi(2)
		print*,"  --- (jlo_d,jhi_d)",d_bndy%jlo(2),d_bndy%jhi(2)
		print*,"  --- BC_info",d_bndy%bc_infos%p
		print*,""
		bndy_node => bndy_node%next
	enddo
	print*,""

	call partision_mesh()
	
	open(1111,file='grid-blk.dat') 
	blk_node => base_blk_list
	do while(associated(blk_node))
		blk   => blk_node%blk
		d_blk => blk%data
		write(1111,*) 'ZONE I=',d_blk%isize+1, &
						' J=',d_blk%jsize+1
		call blk%write_block_ascii(1111)
		blk_node => blk_node%next
	enddo
	close(1111)	

	open(2222,file='grid-tiles.dat') 
	do itile = 1, blk_count
		tile   => blk_que(itile)%blk
		d_tile => tile%data
		bndy_node => tile%bndy_list
		print*,"Tile-",tile%idx,itile
		print*,"Range",tile%p_ilo,tile%p_ihi,tile%p_jlo,tile%p_jhi
		print*,"Num of Boundaries of tile:",tile%bndy_num
		do while(associated(bndy_node))
			bndy   => bndy_node%bndy
			d_bndy => bndy%data
			print*,"bndy:",bndy%idx, bndy%group_id
			print*,"--- (ilo,ihi)",d_bndy%ilo(1),d_bndy%ihi(1)
			print*,"--- (jlo,jhi)",d_bndy%jlo(1),d_bndy%jhi(1)
			print*,"--- (ilo_d,ihi_d)",d_bndy%ilo(2),d_bndy%ihi(2)
			print*,"--- (jlo_d,jhi_d)",d_bndy%jlo(2),d_bndy%jhi(2)
			print*,"--- BC_info",d_bndy%bc_infos%p
			bndy_node => bndy_node%next
		enddo
		print*,""
		write(2222,*) 'ZONE I=',d_tile%isize+1, &
						' J=',d_tile%jsize+1
		call tile%write_block_ascii(2222)
	enddo
	close(2222)	
	deallocate(base_blk_list)

end subroutine split_block_test

program dense_multiphase
	use simula
	implicit none

	write(*,*) "********************************************************************************"
	write(*,*) 
	write(*,*) "                              Two dimensional                                   "
	write(*,*) "                   Dense Material Multiple-phase Simulation                     "
	write(*,*) "                                   Program                                      "
	write(*,*) 
	write(*,*) "                                 By Jun Peng                                    "
	write(*,*) "                           Email:pengjun@imech.ac.cn                            "
	write(*,*) 
	write(*,*) "********************************************************************************"

	! call write_simula_inp()
	call create_simula()
	call initial_simula()
	call evolve_simula()
	call destroy_simula()
	
end program dense_multiphase
