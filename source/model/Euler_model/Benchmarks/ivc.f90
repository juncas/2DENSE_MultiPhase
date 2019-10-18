submodule (eu_model) ivc
implicit none

contains
module subroutine initial_model_eu_ivc(funit)
    implicit none
    integer,intent(in) :: funit
  
    type(eu_data_list_t), pointer :: flwd_node
    type(eu_data_t), pointer      :: flwd
    type(eu_info_data_t),pointer  :: info
    type(rvectorfield_t),pointer  :: prim
    type(rvectorfield_t),pointer  :: cons 
  
    class(block_list_t),pointer   :: blk_node 
    class(block_t),pointer        :: blk
    class(block_data_t),pointer   :: d_blk
    class(bndy_list_t),pointer    :: bndy_node 
    class(bndy_t),pointer         :: bndy 
    class(bndy_data_t),pointer    :: d_bndy 
  
    integer :: iflwd
    integer :: i,j,iblk,bc_num
    integer :: i0,j0
    integer :: ilo(2),ihi(2),jlo(2),jhi(2)
    double precision :: xc0,yc0,rc0
    double precision :: rho0,u0,v0,p0,ei0
    double precision,parameter :: pie = 4.d0*atan(1.d0)
  
    ! Create base mesh
      blk_node  => add_blk_to_mesh()
      blk       => blk_node%blk
      call blk%create_block(isize_plane,jsize_plane,ibf_plane,jbf_plane,1)
    d_blk => blk%data
  
      ilo       = (/1, d_blk%isize + 1/)
      ihi       = (/d_blk%ibf, d_blk%isize + d_blk%ibf/)
      jlo       = (/d_blk%jbf + 1, d_blk%jbf + 1/)
      jhi       = (/d_blk%jbf + d_blk%jsize, d_blk%jbf + d_blk%jsize/)
      bndy_node => blk%add_bndy(left_bndy_id, ilo, ihi, jlo, jhi)	
      d_bndy    => bndy_node%bndy%data
      call create(d_bndy%bc_infos,2)
      d_bndy%bc_infos%p(1) = bc_idx_intf
    d_bndy%bc_infos%p(2) = 1  
    
      ilo       = (/d_blk%ibf + 1, d_blk%ibf + 1/)
      ihi       = (/d_blk%ibf + d_blk%isize, d_blk%ibf + d_blk%isize/)
      jlo       = (/1, d_blk%jsize + 1/)
    jhi       = (/d_blk%jbf, d_blk%jsize + d_blk%jbf/)
      bndy_node => blk%add_bndy(bottom_bndy_id, ilo, ihi, jlo, jhi)	
      d_bndy    => bndy_node%bndy%data
      call create(d_bndy%bc_infos,2)
      d_bndy%bc_infos%p(1) = bc_idx_intf
    d_bndy%bc_infos%p(2) = 1
    
      ilo       = (/d_blk%ibf + d_blk%isize + 1, d_blk%ibf + 1/)
      ihi       = (/d_blk%ibf + d_blk%isize + d_blk%ibf, d_blk%ibf + d_blk%ibf/)
      jlo       = (/d_blk%jbf + 1, d_blk%jbf + 1/)
      jhi       = (/d_blk%jbf + d_blk%jsize, d_blk%jbf + d_blk%jsize/)
      bndy_node => blk%add_bndy(right_bndy_id, ilo, ihi, jlo, jhi)	
      d_bndy    => bndy_node%bndy%data
      call create(d_bndy%bc_infos,2)
      d_bndy%bc_infos%p(1) = bc_idx_intf
    d_bndy%bc_infos%p(2) = 1
    
      ilo       = (/d_blk%ibf + 1, d_blk%ibf + 1/)
      ihi       = (/d_blk%ibf + d_blk%isize, d_blk%ibf + d_blk%isize/)
      jlo       = (/d_blk%jbf + d_blk%jsize + 1, d_blk%jbf + 1/)
      jhi       = (/d_blk%jbf + d_blk%jsize + d_blk%jbf, d_blk%jbf + d_blk%jbf/)
      bndy_node => blk%add_bndy(top_bndy_id, ilo, ihi, jlo, jhi)	
      d_bndy    => bndy_node%bndy%data
      call create(d_bndy%bc_infos,2)
      d_bndy%bc_infos%p(1) = bc_idx_intf
    d_bndy%bc_infos%p(2) = 1
    
    call blk%initial_block(x0,y0,xl,yl)
    
    call create_work_que()
  
    do iflwd = 1, flwd_count
      flwd => flwd_que(iflwd)%flwd
      blk  => flwd%blk
      info => flwd%info
      cons => flwd%cons(max_stage)
      info%gamma = gamma
      do j=1,info%jsize
        do i=1,info%isize
          i0 = i + info%ibf
          j0 = j + info%jbf
    
          xc0 = blk%data%xc%p(i0,j0)
          yc0 = blk%data%yc%p(i0,j0)
          rc0 = sqrt((xc0-5.d0)**2+(yc0-5.d0)**2)      
          
          rho0 = (1.d0-(gamma-1.d0)*b0**2*exp(1.d0-rc0**2)/(8.d0*gamma*pie**2))**2.5d0      
          u0   = uinf - b0/(2.d0*pie)*exp((1.d0-rc0**2)/2.d0)*(yc0-5.d0)      
          v0   = vinf + b0/(2.d0*pie)*exp((1.d0-rc0**2)/2.d0)*(xc0-5.d0)
          p0   = rho0**(gamma)
          ei0  = p0/(rho0*(gamma-1.d0))
    
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
  end subroutine initial_model_eu_ivc 
end submodule ivc