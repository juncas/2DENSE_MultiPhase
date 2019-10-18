module eu_model
  use data_type
  use reconstruction
  use gradient
  use time_advance
  use mesh
  use eu_data_type
  use eu_physics
  use eu_solver
  use eu_boundary
  implicit none
  ! Riemann solver indices
  integer,parameter :: benchmark_idx_ivc           = 1
  integer,parameter :: benchmark_idx_shockvortex   = 2
  integer,parameter :: benchmark_idx_doublemach    = 3
  integer,parameter :: benchmark_idx_shockshear    = 4 
  integer,parameter :: benchmark_idx_shockboundary = 5 
  integer,parameter :: benchmark_idx_sedov         = 6 
  integer,parameter :: benchmark_idx_rt            = 7 
  integer,parameter :: benchmark_idx_rm            = 8 
  integer,parameter :: benchmark_idx_forwardstep   = 9 
  integer,parameter :: benchmark_idx_backwardstep  = 10 
  integer,parameter :: benchmark_idx_sod           = 11 
  integer,parameter :: benchmark_idx_lax           = 12 
  integer,parameter :: benchmark_idx_shuosher      = 13 
  integer,parameter :: benchmark_idx_blast         = 14

  character(len=64) :: file_fmt

  integer :: base_flwd_count = 0
  class(eu_data_list_t), pointer :: base_flwd_list, &
                                    base_flwd_list_tail

  integer :: flwd_count
  type(eu_data_iterator_t),dimension(:),allocatable :: flwd_que ! decomposited flwd_que

  ! global model variables
  double precision :: max_eig_xi_g,max_eig_eta_g
  double precision :: dt_g = 1.d0
  
  double precision :: gamma = 1.4d0
  double precision :: Mach_inf = 1.d0
  double precision :: T_inf = 300.d0 ,T_ref = 110.0d0
  integer :: riemann_solver = 0
  integer :: benchmark_idx = 0
  NAMELIST/model_eu_params/ gamma, Mach_inf, T_inf,T_ref, &
                            riemann_solver,               &
                            benchmark_idx,file_fmt
    
  integer :: reconst_xi_idx      = upwind1_id
  integer :: reconst_eta_idx     = upwind1_id
  integer :: grad_xi_idx         = grad_2nd_l1_idx
  integer :: grad_eta_idx        = grad_2nd_l1_idx
  integer :: time_integrator_idx = rk3rdtvd_idx
  integer :: max_stage
  logical :: adaptive_dt = .true.
  double precision :: cfl = 0.5d0
  namelist /numeric_params/ reconst_xi_idx,reconst_eta_idx, &
                            grad_xi_idx,grad_eta_idx, &
                            time_integrator_idx,cfl,dt_g,adaptive_dt

  double precision :: x0,y0,xl,yl
  integer :: isize_plane,jsize_plane
  integer :: ibf_plane,jbf_plane  
  double precision :: uinf = 0.5d0, vinf = 0.d0, b0 = 0.5d0 
  NAMELIST/model_eu_ivc_params/ uinf, vinf, b0,    &
                                x0, y0, xl, yl,           &
                                isize_plane, jsize_plane, &
                                ibf_plane, jbf_plane

  NAMELIST/model_eu_shockvortex_params/ x0, y0, xl, yl,   &
                                        isize_plane,jsize_plane, &
                                        ibf_plane,jbf_plane

  NAMELIST/model_eu_doublemach_params/ x0, y0, xl, yl,   &
                                       isize_plane,jsize_plane, &
                                       ibf_plane,jbf_plane

  NAMELIST/model_eu_shockshear_params/ x0, y0, xl, yl,   &
                                       isize_plane,jsize_plane, &
                                       ibf_plane,jbf_plane
                    
  NAMELIST/model_eu_shockboundary_params/ x0, y0, xl, yl,   &
                                       isize_plane,jsize_plane, &
                                       ibf_plane,jbf_plane
                     
  NAMELIST/model_eu_sedov_params/ x0, y0, xl, yl,   &
                                       isize_plane,jsize_plane, &
                                       ibf_plane,jbf_plane
                   
  NAMELIST/model_eu_rt_params/ x0, y0, xl, yl,   &
                                       isize_plane,jsize_plane, &
                                       ibf_plane,jbf_plane
                        
  NAMELIST/model_eu_rm_params/ x0, y0, xl, yl,   &
                                       isize_plane,jsize_plane, &
                                       ibf_plane,jbf_plane
                               
  double precision :: xa0,ya0,xb0,yb0,xal,yal,xbl,ybl
  integer :: nn
  NAMELIST/model_eu_fs_params/  xa0,ya0,xb0,yb0,    &
                                xal,yal,xbl,ybl,    &
                                nn,                 &
                                ibf_plane,jbf_plane

  NAMELIST/model_eu_bs_params/ x0, y0, xl, yl,   &
                                       isize_plane,jsize_plane, &
                                       ibf_plane,jbf_plane
                    
  NAMELIST/model_eu_sod_params/ x0, y0, xl, yl,   &
                                isize_plane,jsize_plane, &
                                ibf_plane,jbf_plane

  NAMELIST/model_eu_lax_params/ x0, y0, xl, yl,   &
                                isize_plane,jsize_plane, &
                                ibf_plane,jbf_plane

  NAMELIST/model_eu_shuosher_params/ x0, y0, xl, yl,   &
                                isize_plane,jsize_plane, &
                                ibf_plane,jbf_plane

  NAMELIST/model_eu_blast_params/ x0, y0, xl, yl,   &
                                isize_plane,jsize_plane, &
                                ibf_plane,jbf_plane

    
  interface
    module subroutine initial_model_eu_ivc(funit)
      implicit none
      integer,intent(in) :: funit
    end subroutine initial_model_eu_ivc

    module subroutine initial_model_eu_sv(funit)
      implicit none
      integer,intent(in) :: funit
    end subroutine initial_model_eu_sv

    module subroutine initial_model_eu_dm(funit)
      implicit none
      integer,intent(in) :: funit
    end subroutine initial_model_eu_dm

    module subroutine initial_model_eu_ssl(funit)
      implicit none
      integer,intent(in) :: funit
    end subroutine initial_model_eu_ssl

    module subroutine initial_model_eu_sb(funit)
      implicit none
      integer,intent(in) :: funit
    end subroutine initial_model_eu_sb
    
    module subroutine initial_model_eu_sedov(funit)
      implicit none
      integer,intent(in) :: funit
    end subroutine initial_model_eu_sedov

    module subroutine initial_model_eu_rt(funit)
      implicit none
      integer,intent(in) :: funit
    end subroutine initial_model_eu_rt

    module subroutine initial_model_eu_rm(funit)
      implicit none
      integer,intent(in) :: funit
    end subroutine initial_model_eu_rm

    module subroutine initial_model_eu_fs(funit)
      implicit none
      integer,intent(in) :: funit
    end subroutine initial_model_eu_fs

    module subroutine initial_model_eu_bs(funit)
      implicit none
      integer,intent(in) :: funit
    end subroutine initial_model_eu_bs

    module subroutine initial_model_eu_sod(funit)
      implicit none
      integer,intent(in) :: funit
    end subroutine initial_model_eu_sod

    module subroutine initial_model_eu_lax(funit)
      implicit none
      integer,intent(in) :: funit
    end subroutine initial_model_eu_lax

    module subroutine initial_model_eu_shuosher(funit)
      implicit none
      integer,intent(in) :: funit
    end subroutine initial_model_eu_shuosher

    module subroutine initial_model_eu_blast(funit)
      implicit none
      integer,intent(in) :: funit
    end subroutine initial_model_eu_blast
  end interface

contains
subroutine write_model_inp_eu(funit)
  implicit none
  integer,intent(in) :: funit
  write(funit,NML=model_eu_params)
  write(funit,NML=numeric_params)
  write(funit,NML=model_eu_ivc_params)
  write(funit,NML=model_eu_shockvortex_params)
  write(funit,NML=model_eu_doublemach_params)
  write(funit,NML=model_eu_shockshear_params)
  write(funit,NML=model_eu_sod_params)
  write(funit,NML=model_eu_lax_params)
  write(funit,NML=model_eu_shuosher_params)
  write(funit,NML=model_eu_blast_params)
end subroutine write_model_inp_eu

subroutine create_model_eu(funit)
  implicit none
  integer,intent(in) :: funit
  read(funit,NML=model_eu_params)
  read(funit,NML=numeric_params)
  max_stage = set_maxstage(time_integrator_idx)
  select case(benchmark_idx)
    case(benchmark_idx_ivc)
      read(funit,NML=model_eu_ivc_params)
    case(benchmark_idx_shockvortex)
      read(funit,NML=model_eu_shockvortex_params)
    case(benchmark_idx_doublemach)
      read(funit,NML=model_eu_doublemach_params)
    case(benchmark_idx_shockshear)
      read(funit,NML=model_eu_shockshear_params)
    case(benchmark_idx_shockboundary)
    case(benchmark_idx_sedov)
    case(benchmark_idx_rt)
    case(benchmark_idx_rm)
    case(benchmark_idx_forwardstep)
      read(funit,NML=model_eu_fs_params)
    case(benchmark_idx_backwardstep)
    case(benchmark_idx_sod)
      read(funit,NML=model_eu_sod_params)
    case(benchmark_idx_lax)
      read(funit,NML=model_eu_lax_params)
    case(benchmark_idx_shuosher)
      read(funit,NML=model_eu_shuosher_params)
    case(benchmark_idx_blast)
      read(funit,NML=model_eu_blast_params)
    case default
      print*,"Error: create_model_eu()"    
      print*,"Invalid benchmark index: ",benchmark_idx
      stop
  end select
end subroutine create_model_eu

subroutine initial_model_eu(funit)
  implicit none
  integer,intent(in) :: funit
  print*,"Initializing Euler model:"
  select case(benchmark_idx)
  case(benchmark_idx_ivc)
    print*," - Isentropic Vortex Convection."  
    call initial_model_eu_ivc(funit)
  case(benchmark_idx_shockvortex)
    print*," - Shock Vortex Interaction."  
    ! call initial_model_eu_svi(funit)
  case(benchmark_idx_doublemach)
    print*," - Double Mach Reflection."  
    call initial_model_eu_dm(funit)
  ! case(benchmark_idx_shockshear)
  !   print*," - Shock Shear Layer Interaction."  
  !   call initial_model_eu_ss(funit)
  ! case(benchmark_idx_shockboundary)
  !   print*," - Shock Boundary Layer Interaction."  
  !   ! call initial_model_eu_sb(funit)
  ! case(benchmark_idx_sedov)
  !   print*," - Sediv Problem."  
  ! case(benchmark_idx_rt)
  !   print*," - Rayleigh-Taylor Interaction."  
  ! case(benchmark_idx_rm)
  !   print*," - Richtmyer-Meshkov Instability."  
  case(benchmark_idx_forwardstep)
    print*," - Forward Wind Tunnel Step."  
    call initial_model_eu_fs(funit)
  ! case(benchmark_idx_backwardstep)
  !   print*," - Backward Wind Tunnel Step."  
  ! case(benchmark_idx_sod)
  !   print*," - Sod shock tube problem."   
  !   call initial_model_eu_sod(funit)
  ! case(benchmark_idx_lax)
  !   print*," - Lax shock tube problem."   
  !   call initial_model_eu_lax(funit)
  ! case(benchmark_idx_shuosher)
  !   print*," - Shu-Osher problem."   
  !   call initial_model_eu_shuosher(funit)
  ! case(benchmark_idx_blast)
  !   print*," - Blast wave problem."   
  !   call initial_model_eu_blast(funit)
  case default
    print*,"Error: initial_model_eu()"    
    print*,"Invalid benchmark index: ",benchmark_idx
    stop
  end select
end subroutine initial_model_eu

subroutine evolve_model_eu(istep,time,elasp_time)
  implicit none
  integer,intent(inout) :: istep
  double precision,intent(inout) :: time,elasp_time
  type(eu_data_t),pointer :: flwd
  integer :: istage,iflwd
  istage = 1
  do while(istage.le.max_stage) 
    !$OMP PARALLEL DO PRIVATE(iflwd,flwd) &
    !$OMP SHARED(flwd_que)
    do iflwd = 1, flwd_count   
      flwd => flwd_que(iflwd)%flwd
      call flwd%set_stage(istage)
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(iflwd,flwd) &
    !$OMP SHARED(flwd_que)
    do iflwd = 1, flwd_count   
      flwd => flwd_que(iflwd)%flwd
      flwd%info%eig_xi  = 0.d0
      flwd%info%eig_eta = 0.d0
    enddo
    !$OMP END PARALLEL DO
    
    ! Convert conservative variables to primatives for interior grids
    !$OMP PARALLEL DO PRIVATE(iflwd,flwd) &
    !$OMP SHARED(flwd_que)
    do iflwd = 1, flwd_count
      flwd => flwd_que(iflwd)%flwd
      call cmpt_prim_eu(flwd)
    enddo
    !$OMP END PARALLEL DO

    ! Set boundary condition and convert cons to prims or prims to cons
    ! according to bc type
    !$OMP PARALLEL DO PRIVATE(iflwd) &
    !$OMP SHARED(flwd_que)
    do iflwd = 1, flwd_count
      call set_bc_eu(flwd_que,flwd_count,iflwd)
    enddo
    !$OMP END PARALLEL DO

    max_eig_xi_g  = 0.d0
    max_eig_eta_g = 0.d0
    do iflwd = 1, flwd_count     
      flwd => flwd_que(iflwd)%flwd 
      max_eig_xi_g  = max(max_eig_xi_g,flwd%info%eig_xi)
      max_eig_eta_g = max(max_eig_eta_g,flwd%info%eig_eta)
    enddo

    !$OMP PARALLEL DO PRIVATE(iflwd,flwd) &
    !$OMP SHARED(flwd_que)
    do iflwd = 1, flwd_count    
      flwd => flwd_que(iflwd)%flwd   
      flwd%info%eig_xi_g  = max_eig_xi_g
      flwd%info%eig_eta_g = max_eig_eta_g
    enddo
    !$OMP END PARALLEL DO

    ! Convective terms
    !$OMP PARALLEL DO PRIVATE(iflwd,flwd) &
    !$OMP SHARED(flwd_que)
    do iflwd = 1, flwd_count
      flwd => flwd_que(iflwd)%flwd 
      call cmpt_cnvt_eu(flwd,           &
                        riemann_solver, &
                        reconst_xi_idx, &
                        reconst_eta_idx)
    enddo
    !$OMP END PARALLEL DO

    ! Diffusive terms
    ! do iflwd = 1, flwd_count
    !   flwd => flwd_que(iflwd)%flwd 
    !   call cmpt_dfsv_eu(flwd,        &
    !                     grad_xi_idx, &
    !                     grad_eta_idx)
    ! enddo

    ! do iflwd = 1, flwd_count
    !   flwd => flwd_que(iflwd)%flwd 
    !   call cmpt_lhs_eu(flwd)
    ! enddo

    !$OMP PARALLEL DO PRIVATE(iflwd,flwd) &
    !$OMP SHARED(flwd_que)
    do iflwd = 1, flwd_count
      flwd => flwd_que(iflwd)%flwd 
      call cmpt_rhs_eu(flwd)
    enddo
    !$OMP END PARALLEL DO

    if((istage.eq.1) .and. &
       (adaptive_dt)) then
      !$OMP PARALLEL DO PRIVATE(iflwd,flwd) &
      !$OMP SHARED(flwd_que)
      do iflwd = 1, flwd_count
        flwd => flwd_que(iflwd)%flwd 
        call cmpt_dt_eu(flwd)
      enddo
      !$OMP END PARALLEL DO

      dt_g = flwd_que(1)%flwd%info%dt
      do iflwd = 2, flwd_count  
        flwd => flwd_que(iflwd)%flwd  
        dt_g = min(dt_g,flwd%info%dt)
      enddo
      dt_g = cfl*dt_g
      
      !$OMP PARALLEL DO PRIVATE(iflwd,flwd) &
      !$OMP SHARED(flwd_que)
      do iflwd = 1, flwd_count  
        flwd => flwd_que(iflwd)%flwd  
        flwd%info%dt = dt_g
      enddo
      !$OMP END PARALLEL DO
    endif

    !$OMP PARALLEL DO PRIVATE(iflwd,flwd) &
    !$OMP SHARED(flwd_que)
    do iflwd = 1, flwd_count
      flwd => flwd_que(iflwd)%flwd 
      call evlv_mdl_eu(flwd, time_integrator_idx)
    enddo    
    !$OMP END PARALLEL DO

    istage = istage + 1
  enddo
  
  istep = istep + 1
  time  = time + dt_g
  elasp_time  = elasp_time + dt_g

  !$OMP PARALLEL DO PRIVATE(iflwd,flwd) &
  !$OMP SHARED(flwd_que)
  do iflwd = 1, flwd_count
    flwd => flwd_que(iflwd)%flwd 
    flwd%info%time = time
  enddo
  !$OMP END PARALLEL DO
end subroutine evolve_model_eu

subroutine print_status_model_eu(istep,time)
  implicit none
  integer,intent(in) :: istep
  double precision,intent(in) :: time

  write(*,*) "======== istep = ",istep,"time = ",time,"=========="
  write(*,*) "dt = ", dt_g
  write(*,*) 'max eigen value xi', max_eig_xi_g
  write(*,*) 'max eigen value eta', max_eig_eta_g

end subroutine print_status_model_eu

subroutine write_model_eu(istep,time)
  implicit none
  integer,intent(in) :: istep
  double precision,intent(in) :: time  
  type(eu_data_list_t),pointer :: flwd_node

  type(eu_data_t),pointer      :: flwd
  type(block_t),pointer        :: blk
  type(block_data_t),pointer   :: d_blk
  type(eu_info_data_t),pointer :: info
  type(rvectorfield_t),pointer :: prim
  type(rvectorfield_t),pointer :: cons

  type(eu_data_t),pointer      :: flwd_tile
  type(block_t),pointer        :: blk_tile
  type(block_data_t),pointer   :: d_blk_tile
  type(eu_info_data_t),pointer :: info_tile
  type(rvectorfield_t),pointer :: prim_tile
  type(rvectorfield_t),pointer :: cons_tile

  character(len=32) :: sstep,stime 
  character(len=1024) :: varnamelist 
  character(len=1024) :: header 
  integer :: funit = 11
  integer :: iflwd,nvars
  integer :: isize,jsize,ibf,jbf
  integer :: i,j  
  integer :: i0,j0
  integer :: itile,jtile

  type(fscalarshape_t) :: lo_int, hi_int
  
  ! Compute current primative variables for all tiles
  do iflwd = 1, flwd_count
    flwd => flwd_que(iflwd)%flwd
    info => flwd%info
    prim => flwd%prim
    cons => flwd%cons(info%stage_num)

    isize   = info%isize
    jsize   = info%jsize
    ibf     = info%ibf
    jbf     = info%jbf

    lo_int = set_shape(1+ibf,1+jbf)
    hi_int = set_shape(isize+ibf,jsize+jbf)
    call cmpt_cons2prim_range(prim, cons, lo_int, hi_int, info) 
  enddo

  flwd_node => base_flwd_list
  do while(associated(flwd_node))
    flwd => flwd_node%flwd
    blk  => flwd%blk
    prim => flwd%prim
    do jtile = 1, blk%tile_dim(2)
      do itile = 1, blk%tile_dim(1)
        flwd_tile => flwd%tiles(itile,jtile)
        blk_tile  => flwd_tile%blk
        info_tile => flwd_tile%info

        prim_tile => flwd_tile%prim

        isize   = info_tile%isize
        jsize   = info_tile%jsize
        ibf     = info_tile%ibf
        jbf     = info_tile%jbf

        do j = 0, jsize - 1
          do i = 0, isize - 1
            i0 = blk_tile%p_ilo + i
            j0 = blk_tile%p_jlo + j
            prim%p(:,i0,j0) = prim_tile%p(:,ibf + i+1,jbf + j+1)
          enddo
        enddo
      enddo
    enddo    
    flwd_node => flwd_node%next
  enddo

  if(file_fmt.eq.'tec') then
    call write_model_tec_eu(istep,time)
  else if(file_fmt.eq.'hdf') then
    call write_model_hdf5_eu(istep,time)
  else 
    print*,'No such File Format:',file_fmt
  endif  
end subroutine write_model_eu

subroutine write_model_tec_eu(istep,time)
  implicit none
  integer,intent(in) :: istep
  double precision,intent(in) :: time  
  type(eu_data_list_t),pointer :: flwd_node

  type(eu_data_t),pointer      :: flwd
  type(block_t),pointer        :: blk
  type(block_data_t),pointer   :: d_blk
  type(eu_info_data_t),pointer :: info
  type(rvectorfield_t),pointer :: prim
  type(rvectorfield_t),pointer :: cons

  type(eu_data_t),pointer      :: flwd_tile
  type(block_t),pointer        :: blk_tile
  type(block_data_t),pointer   :: d_blk_tile
  type(eu_info_data_t),pointer :: info_tile
  type(rvectorfield_t),pointer :: prim_tile
  type(rvectorfield_t),pointer :: cons_tile

  character(len=32) :: sstep,stime 
  character(len=1024) :: varnamelist 
  character(len=1024) :: header 
  integer :: funit = 11
  integer :: iflwd,nvars
  integer :: isize,jsize,ibf,jbf
  integer :: i,j  
  integer :: i0,j0
  integer :: itile,jtile

  type(fscalarshape_t) :: lo_int, hi_int
  
  write(sstep,'(i16)') istep
  write(stime,'(f12.6)') time
  open(funit,file="DMPeu_"//trim(adjustl(sstep))//"_"//trim(adjustl(stime))//".tec",form='formatted')   

  flwd_node => base_flwd_list
  do while(associated(flwd_node))
    flwd => flwd_node%flwd
    info => flwd%info
    blk  => flwd%blk
    write(funit,*)'title=','"flow time ',time,'"'
    varnamelist = ""
    call get_varnamelist(varnamelist)
    header = 'variables="x","y"'//trim(adjustl(varnamelist))
    write(funit,*) header
    header = "DATAPACKING=BLOCK, VARLOCATION=([3-8]=CELLCENTERED)"
  
    write(funit,*) 'ZONE I=',info%isize+1, &
                    ' J=',info%jsize+1, &
                    trim(adjustl(header))
    call get_numvars(nvars)
    call blk%write_block_ascii(funit)
    call flwd%write_data_fmt(funit)
    flwd_node => flwd_node%next
  enddo
  close(funit)
end subroutine write_model_tec_eu

subroutine write_model_hdf5_eu(istep,time)
  use HDF5
  implicit none
  integer,intent(in) :: istep
  double precision,intent(in) :: time  
  type(eu_data_list_t),pointer :: flwd_node

  type(eu_data_t),pointer      :: flwd
  type(block_t),pointer        :: blk
  type(block_data_t),pointer   :: d_blk
  type(eu_info_data_t),pointer :: info
  type(rvectorfield_t),pointer :: prim
  type(rvectorfield_t),pointer :: cons

  character(len=32) :: sstep,stime 
  character(len=1024) :: varnamelist 
  character(len=1024) :: header 
  integer :: iflwd,nvars
  integer :: isize,jsize,ibf,jbf
  integer :: i,j  
  integer :: i0,j0
  integer :: itile,jtile

  type(fscalarshape_t) :: lo_int, hi_int
  
  flwd_node => base_flwd_list
  do while(associated(flwd_node))
    flwd => flwd_node%flwd
    info => flwd%info
    blk  => flwd%blk
    write(funit,*)'title=','"flow time ',time,'"'
    varnamelist = ""
    call get_varnamelist(varnamelist)
    header = 'variables="x","y"'//trim(adjustl(varnamelist))
    write(funit,*) header
    header = "DATAPACKING=BLOCK, VARLOCATION=([3-8]=CELLCENTERED)"
  
    write(funit,*) 'ZONE I=',info%isize+1, &
                    ' J=',info%jsize+1, &
                    trim(adjustl(header))
    call get_numvars(nvars)
    call blk%write_block_ascii(funit)
    call flwd%write_data_hdf5(funit)
    flwd_node => flwd_node%next
  enddo
end subroutine write_model_hdf5_eu

subroutine save_model_eu(istep,time)
  implicit none
  integer,intent(in) :: istep
  double precision,intent(in) :: time  
  character(len=32) :: sstep,stime 
  integer :: funit = 11
  integer :: iflwd
  type(fscalarshape_t) :: lo_int, hi_int

end subroutine save_model_eu

subroutine load_model_eu(istep,time)
  implicit none
  integer,intent(in) :: istep
  double precision,intent(in) :: time  
  character(len=32) :: sstep,stime 
  character(len=1024) :: varnamelist 
  character(len=1024) :: header 
  integer :: funit = 11
  integer :: iflwd,nvars
  type(fscalarshape_t) :: lo_int, hi_int

end subroutine load_model_eu

subroutine destroy_model_eu()
  implicit none
  integer :: iflwd
  deallocate(base_flwd_list)
  deallocate(flwd_que)
  call destroy_mesh()
end subroutine destroy_model_eu

function add_flwd_to_model() result(flwd_node)
  implicit none
  class(eu_data_list_t),pointer :: flwd_node 
  allocate(flwd_node)
  allocate(flwd_node%flwd)
  flwd_node%prev  => null()
  flwd_node%next => null()
  if(.not.associated(base_flwd_list)) then 
    base_flwd_list      => flwd_node
    base_flwd_list_tail => base_flwd_list
  else
    base_flwd_list_tail%next => flwd_node
    flwd_node%prev           => base_flwd_list_tail
    base_flwd_list_tail      => flwd_node
  endif
end function add_flwd_to_model

subroutine map_flwd_to_base_mesh()
  implicit none
  class(eu_data_list_t),pointer :: flwd_node 
  class(block_list_t),pointer   :: blk_node
  class(block_t),pointer        :: blk 
  class(eu_data_t),pointer      :: flwd 
  integer :: itile,jtile
  integer :: iblk

  blk_node => base_blk_list
  do while(associated(blk_node))
    blk       => blk_node%blk
    flwd_node => add_flwd_to_model()
    flwd      => flwd_node%flwd
    call flwd%config_eu(max_stage,blk)
    call flwd%create_eu()
    
    blk_node => blk_node%next
  enddo
end subroutine map_flwd_to_base_mesh

subroutine create_work_que()
  use omp_lib
  implicit none
  class(eu_data_list_t),pointer               :: flwd_node 
  class(block_t),pointer                      :: blk 
  class(eu_data_t),pointer                    :: flwd 
  class(block_t),pointer                      :: blk_tile
  class(eu_data_t),pointer                    :: flwd_tile
  integer :: itile,jtile
  integer :: max_num_threads

  call partision_mesh()
  call map_flwd_to_base_mesh()
  flwd_count = blk_count
  if(allocated(flwd_que)) deallocate(flwd_que)
  allocate(flwd_que(flwd_count))

  flwd_node => base_flwd_list
  do while(associated(flwd_node))
    flwd      => flwd_node%flwd
    blk       => flwd%blk
    if(blk%partisioned) then
      if(associated(flwd%tiles)) deallocate(flwd%tiles)
      allocate(flwd%tiles(blk%tile_dim(1),blk%tile_dim(2)))
      do jtile = 1, blk%tile_dim(2)
        do itile = 1, blk%tile_dim(1)
          blk_tile  => blk%tiles(itile,jtile)
          flwd_tile => flwd%tiles(itile,jtile)
          call flwd_tile%config_eu(max_stage,blk_tile)
          call flwd_tile%create_eu()
          flwd_tile%parent            => flwd
          flwd_que(blk_tile%idx)%flwd => flwd_tile
        enddo
      enddo
    endif    
    flwd_node => flwd_node%next
  enddo
  print*,"Work flow has been created:",blk_count,"tiles in total."  
  max_num_threads = omp_get_max_threads()
  call omp_set_num_threads(min(max_num_threads,blk_count))
  max_num_threads = omp_get_max_threads()
  print*,"Using ",max_num_threads,"threads."  
end subroutine create_work_que

! ! Shock Shear layer
! subroutine initial_model_eu_ss(funit)
!   implicit none
!   integer,intent(in) :: funit
!   integer,parameter :: blk_num = 1

!   type(eu_data_t), pointer     :: flwd
!   type(block_t), pointer       :: blk
!   type(eu_info_data_t),pointer :: info
!   type(rvectorfield_t),pointer :: prim
!   type(rvectorfield_t),pointer :: cons 
!   type(blk_bndy_t), pointer    :: bndy
!   type(eu_bc_data_t), pointer  :: bc
!   type(ivector_t), pointer     :: bc_info

!   integer :: i,j,iblk,bc_num
!   integer :: i0,j0
!   double precision :: xc0,yc0,rc0
!   double precision :: rho0,u0,v0,p0,ei0
!   double precision,parameter :: pie = 4.d0*atan(1.d0)

!   ! initial mesh data
!   call create_mesh(blk_num,mesh_type_cart)
!   call allocate_block(blks(blk_num),blk_num,   &
!                       isize_plane,jsize_plane, &
!                       ibf_plane,jbf_plane)

!   do iblk = 1, nblk
!     call init_cart_blk(blks(iblk),x0,y0,xl,yl)
!   enddo

!   call cmpt_metrics_cart()

!   blks(1)%nbndy = 4
!   allocate(blks(1)%bndys(blks(1)%nbndy))
!   ! boundary conditions
!   ! left boundary
!   bndy => blks(1)%bndys(1)
!   bndy%idx    = 1

!   bndy%ilo(1)  = 1
!   bndy%ihi(1)  = ibf_plane
!   bndy%jlo(1)  = jbf_plane + 1
!   bndy%jhi(1)  = jbf_plane + jsize_plane
  
!   bndy%bc_num = 1
!   allocate(bndy%bc_infos(bndy%bc_num))
!   bc_info => bndy%bc_infos(1)
!   call create(bc_info,2)
!   bc_info%p(1) = ssl_bc_idx
!   bc_info%p(2) = inlt1_ssl_idx

!   ! bottom boundary
!   bndy => blks(1)%bndys(2)
!   bndy%idx    = 2
  
!   bndy%ilo(1) = ibf_plane + 1
!   bndy%ihi(1) = ibf_plane + isize_plane
!   bndy%jlo(1) = 1
!   bndy%jhi(1) = jbf_plane

!   bndy%ilo(2) = ibf_plane + 1
!   bndy%ihi(2) = ibf_plane + isize_plane
!   bndy%jlo(2) = jbf_plane + jbf_plane
!   bndy%jhi(2) = jbf_plane + 1

!   bndy%bc_num = 1
!   allocate(bndy%bc_infos(bndy%bc_num))
!   bc_info => bndy%bc_infos(1)
!   call create(bc_info,2)
!   bc_info%p(1) = ssl_bc_idx
!   bc_info%p(2) = rflt2_ssl_idx

!   ! right boundary
!   bndy => blks(1)%bndys(3)
!   bndy%idx    = 3

!   bndy%ilo(1)  = isize_plane + ibf_plane + 1
!   bndy%ihi(1)  = isize_plane + ibf_plane + ibf_plane  
!   bndy%jlo(1)  = jbf_plane + 1 
!   bndy%jhi(1)  = jbf_plane + jsize_plane 

!   bndy%bc_num = 1
!   allocate(bndy%bc_infos(bndy%bc_num))
!   bc_info => bndy%bc_infos(1)
!   call create(bc_info,3)
!   bc_info%p(1) = ssl_bc_idx
!   bc_info%p(2) = otlt3_ssl_idx
!   bc_info%p(3) = axis_xi
    
!   ! top boundary
!   bndy => blks(1)%bndys(4)
!   bndy%idx    = 4

!   bndy%ilo(1)  = ibf_plane + 1
!   bndy%ihi(1)  = ibf_plane + isize_plane 
!   bndy%jlo(1)  = jsize_plane + jbf_plane + 1
!   bndy%jhi(1)  = jsize_plane + jbf_plane + jbf_plane

!   bndy%bc_num = 1
!   allocate(bndy%bc_infos(bndy%bc_num))
!   bc_info => bndy%bc_infos(1)
!   call create(bc_info,2)
!   bc_info%p(1) = ssl_bc_idx
!   bc_info%p(2) = inlt4_ssl_idx

!   ! partiiton mesh for parallel computation
!   call mesh_partition()
  
!   ! initial model data
!   allocate(flwd_que(1))
!   call flwd_que(1)%config_eu(max_stage, blks(1))
!   call flwd_que(1)%create_eu
!   flwd => flwd_que(1)
!   blk  => flwd%blk
!   info => flwd%info
!   cons => flwd%cons(max_stage)

!   ! Setup boundary condition 2 data
!   bc => flwd%bcs(2)
!   bc%vs_num = 1
!   allocate(bc%vs(bc%vs_num))
!   call create(bc%vs(1), 4)
!   bc%vs(1)%p(:) = 1.d0
!   bc%vs(1)%p(3) =-1.d0

!   ! Setup boundary condition 3 data
!   bc => flwd%bcs(3)
!   bc%vs_num = 1
!   allocate(bc%vs(bc%vs_num))
!   call create(bc%vs(1), 4)
!   bc%vs(1)%p(:) = 0.d0

!   info%gamma = gamma
!   do j=1,info%jsize
!     do i=1,info%isize
!       i0 = i + info%ibf
!       j0 = j + info%jbf

!       xc0 = blk%xc%p(i0,j0)
!       yc0 = blk%yc%p(i0,j0)

!       if(yc0.ge.0) then
!         rho0 = 1.6374d0
!         u0   = 3.d0
!         v0   = 0.d0
!         p0   = 0.3327d0
!       else
!         rho0 = 0.3626d0
!         u0   = 2.d0
!         v0   = 0.d0
!         p0   = 0.3327d0
!       endif
!       ei0  = p0/(rho0*(gamma-1.d0))

!       cons%p(1,i0,j0) = rho0
!       cons%p(2,i0,j0) = rho0*u0
!       cons%p(3,i0,j0) = rho0*v0
!       cons%p(4,i0,j0) = rho0*(ei0+0.5d0*(u0**2+v0**2))
!     enddo
!   enddo
!   info%eig_xi  = 0.d0
!   info%eig_eta = 0.d0
!   info%time    = 0.d0
! end subroutine initial_model_eu_ss


subroutine get_varnamelist(varnamelist)
  implicit none
  character(len=*),intent(inout) :: varnamelist
  integer :: n,iphase
  character(len=16),dimension(8) :: svarname
  character(len=8)  :: svarindex

  svarname(1) = ',"rho"'
  svarname(2) = ',"u"'
  svarname(3) = ',"v"'
  svarname(4) = ',"p"'
  svarname(5) = ',"c"'
  svarname(6) = ',"t"'

  do n=1, 6
    varnamelist = trim(adjustl(varnamelist))//trim(adjustl(svarname(n)))
  enddo

end subroutine get_varnamelist

subroutine get_numvars(nvars)
  implicit none
  integer,intent(inout) :: nvars  
  nvars = 6
end subroutine get_numvars

end module
