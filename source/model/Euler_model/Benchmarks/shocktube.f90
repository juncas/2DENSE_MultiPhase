submodule (eu_model) shocktube
implicit none

contains

! ! Sod
! module subroutine initial_model_eu_sod(funit)
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

!   bndy%ilo(1)  = ibf_plane
!   bndy%ihi(1)  = 1
!   bndy%jlo(1)  = jbf_plane + 1
!   bndy%jhi(1)  = jbf_plane + jsize_plane
  
!   bndy%bc_num = 1
!   allocate(bndy%bc_infos(bndy%bc_num))
!   bc_info => bndy%bc_infos(1)
!   call create(bc_info,3)
!   bc_info%p(1) = shocktube_bc_idx
!   bc_info%p(2) = otlt1_st_idx
!   bc_info%p(3) = axis_xi

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
!   bc_info%p(1) = shocktube_bc_idx
!   bc_info%p(2) = rflt2_st_idx

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
!   bc_info%p(1) = shocktube_bc_idx
!   bc_info%p(2) = otlt3_st_idx
!   bc_info%p(3) = axis_xi
    
!   ! top boundary
!   bndy => blks(1)%bndys(4)
!   bndy%idx    = 4

!   bndy%ilo(1)  = ibf_plane + 1
!   bndy%ihi(1)  = ibf_plane + isize_plane 
!   bndy%jlo(1)  = jsize_plane + jbf_plane + 1
!   bndy%jhi(1)  = jsize_plane + jbf_plane + jbf_plane

!   bndy%ilo(2)  = ibf_plane + 1
!   bndy%ihi(2)  = ibf_plane + isize_plane 
!   bndy%jlo(2)  = jsize_plane + jbf_plane
!   bndy%jhi(2)  = jsize_plane + 1

!   bndy%bc_num = 1
!   allocate(bndy%bc_infos(bndy%bc_num))
!   bc_info => bndy%bc_infos(1)
!   call create(bc_info,2)
!   bc_info%p(1) = shocktube_bc_idx
!   bc_info%p(2) = rflt4_st_idx

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

!   ! Setup boundary condition 1 data
!   bc => flwd%bcs(1)
!   bc%vs_num = 1
!   allocate(bc%vs(bc%vs_num))
!   call create(bc%vs(1), 4)
!   bc%vs(1)%p(:) = 0.d0

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

!   ! Setup boundary condition 4 data
!   bc => flwd%bcs(4)
!   bc%vs_num = 1
!   allocate(bc%vs(bc%vs_num))
!   call create(bc%vs(1), 4)
!   bc%vs(1)%p(:) = 1.d0
!   bc%vs(1)%p(3) =-1.d0

!   info%gamma = gamma
!   do j=1,info%jsize
!     do i=1,info%isize
!       i0 = i + info%ibf
!       j0 = j + info%jbf

!       xc0 = blk%xc%p(i0,j0)
!       yc0 = blk%yc%p(i0,j0)

!       if(xc0.lt.0.d0) then
!         rho0 = 1.d0
!         u0   = 0.d0
!         v0   = 0.d0
!         p0   = 1.d0
!       else
!         rho0 = 0.125d0
!         u0   = 0.d0
!         v0   = 0.d0
!         p0   = 0.1d0
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
! end subroutine initial_model_eu_sod

! ! Lax
! module subroutine initial_model_eu_lax(funit)
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

!   bndy%ilo(1)  = ibf_plane
!   bndy%ihi(1)  = 1
!   bndy%jlo(1)  = jbf_plane + 1
!   bndy%jhi(1)  = jbf_plane + jsize_plane
  
!   bndy%bc_num = 1
!   allocate(bndy%bc_infos(bndy%bc_num))
!   bc_info => bndy%bc_infos(1)
!   call create(bc_info,3)
!   bc_info%p(1) = shocktube_bc_idx
!   bc_info%p(2) = otlt1_st_idx
!   bc_info%p(3) = axis_xi

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
!   bc_info%p(1) = shocktube_bc_idx
!   bc_info%p(2) = rflt2_st_idx

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
!   bc_info%p(1) = shocktube_bc_idx
!   bc_info%p(2) = otlt3_st_idx
!   bc_info%p(3) = axis_xi
    
!   ! top boundary
!   bndy => blks(1)%bndys(4)
!   bndy%idx    = 4

!   bndy%ilo(1)  = ibf_plane + 1
!   bndy%ihi(1)  = ibf_plane + isize_plane 
!   bndy%jlo(1)  = jsize_plane + jbf_plane + 1
!   bndy%jhi(1)  = jsize_plane + jbf_plane + jbf_plane

!   bndy%ilo(2)  = ibf_plane + 1
!   bndy%ihi(2)  = ibf_plane + isize_plane 
!   bndy%jlo(2)  = jsize_plane + jbf_plane
!   bndy%jhi(2)  = jsize_plane + 1

!   bndy%bc_num = 1
!   allocate(bndy%bc_infos(bndy%bc_num))
!   bc_info => bndy%bc_infos(1)
!   call create(bc_info,2)
!   bc_info%p(1) = shocktube_bc_idx
!   bc_info%p(2) = rflt4_st_idx

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

!   ! Setup boundary condition 1 data
!   bc => flwd%bcs(1)
!   bc%vs_num = 1
!   allocate(bc%vs(bc%vs_num))
!   call create(bc%vs(1), 4)
!   bc%vs(1)%p(:) = 0.d0

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

!   ! Setup boundary condition 4 data
!   bc => flwd%bcs(4)
!   bc%vs_num = 1
!   allocate(bc%vs(bc%vs_num))
!   call create(bc%vs(1), 4)
!   bc%vs(1)%p(:) = 1.d0
!   bc%vs(1)%p(3) =-1.d0

!   info%gamma = gamma
!   do j=1,info%jsize
!     do i=1,info%isize
!       i0 = i + info%ibf
!       j0 = j + info%jbf

!       xc0 = blk%xc%p(i0,j0)
!       yc0 = blk%yc%p(i0,j0)

!       if(xc0.lt.0.d0) then
!         rho0 = 0.445d0
!         u0   = 0.698d0
!         v0   = 0.d0
!         p0   = 3.528d0
!       else
!         rho0 = 0.5d0
!         u0   = 0.d0
!         v0   = 0.d0
!         p0   = 0.571d0
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
! end subroutine initial_model_eu_lax

! ! Shu-Osher
! module subroutine initial_model_eu_shuosher(funit)
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

!   bndy%ilo(1)  = ibf_plane
!   bndy%ihi(1)  = 1
!   bndy%jlo(1)  = jbf_plane + 1
!   bndy%jhi(1)  = jbf_plane + jsize_plane
  
!   bndy%bc_num = 1
!   allocate(bndy%bc_infos(bndy%bc_num))
!   bc_info => bndy%bc_infos(1)
!   call create(bc_info,3)
!   bc_info%p(1) = shocktube_bc_idx
!   bc_info%p(2) = otlt1_st_idx
!   bc_info%p(3) = axis_xi

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
!   bc_info%p(1) = shocktube_bc_idx
!   bc_info%p(2) = rflt2_st_idx

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
!   bc_info%p(1) = shocktube_bc_idx
!   bc_info%p(2) = otlt3_st_idx
!   bc_info%p(3) = axis_xi
    
!   ! top boundary
!   bndy => blks(1)%bndys(4)
!   bndy%idx    = 4

!   bndy%ilo(1)  = ibf_plane + 1
!   bndy%ihi(1)  = ibf_plane + isize_plane 
!   bndy%jlo(1)  = jsize_plane + jbf_plane + 1
!   bndy%jhi(1)  = jsize_plane + jbf_plane + jbf_plane

!   bndy%ilo(2)  = ibf_plane + 1
!   bndy%ihi(2)  = ibf_plane + isize_plane 
!   bndy%jlo(2)  = jsize_plane + jbf_plane
!   bndy%jhi(2)  = jsize_plane + 1

!   bndy%bc_num = 1
!   allocate(bndy%bc_infos(bndy%bc_num))
!   bc_info => bndy%bc_infos(1)
!   call create(bc_info,2)
!   bc_info%p(1) = shocktube_bc_idx
!   bc_info%p(2) = rflt4_st_idx

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

!   ! Setup boundary condition 1 data
!   bc => flwd%bcs(1)
!   bc%vs_num = 1
!   allocate(bc%vs(bc%vs_num))
!   call create(bc%vs(1), 4)
!   bc%vs(1)%p(:) = 0.d0

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

!   ! Setup boundary condition 4 data
!   bc => flwd%bcs(4)
!   bc%vs_num = 1
!   allocate(bc%vs(bc%vs_num))
!   call create(bc%vs(1), 4)
!   bc%vs(1)%p(:) = 1.d0
!   bc%vs(1)%p(3) =-1.d0

!   info%gamma = gamma
!   do j=1,info%jsize
!     do i=1,info%isize
!       i0 = i + info%ibf
!       j0 = j + info%jbf

!       xc0 = blk%xc%p(i0,j0)
!       yc0 = blk%yc%p(i0,j0)

!       if(xc0.lt.(-4.d0)) then
!         rho0 = 27.d0/7.d0
!         u0   = 4.d0*sqrt(35.d0)/9.d0
!         v0   = 0.d0
!         p0   = 31.d0/3.d0
!       else
!         rho0 = 1.d0+0.2*sin(5.d0*xc0)
!         u0   = 0.d0
!         v0   = 0.d0
!         p0   = 1.d0
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
! end subroutine initial_model_eu_shuosher

! ! Blast
! module subroutine initial_model_eu_blast(funit)
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

!   bndy%ilo(2)  = ibf_plane + ibf_plane
!   bndy%ihi(2)  = ibf_plane + 1
!   bndy%jlo(2)  = jbf_plane + 1
!   bndy%jhi(2)  = jbf_plane + jsize_plane
  
!   bndy%bc_num = 1
!   allocate(bndy%bc_infos(bndy%bc_num))
!   bc_info => bndy%bc_infos(1)
!   call create(bc_info,2)
!   bc_info%p(1) = blast_bc_idx
!   bc_info%p(2) = rflt1_blst_idx

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
!   bc_info%p(1) = blast_bc_idx
!   bc_info%p(2) = rflt2_blst_idx

!   ! right boundary
!   bndy => blks(1)%bndys(3)
!   bndy%idx    = 3

!   bndy%ilo(1)  = isize_plane + ibf_plane + 1
!   bndy%ihi(1)  = isize_plane + ibf_plane + ibf_plane  
!   bndy%jlo(1)  = jbf_plane + 1 
!   bndy%jhi(1)  = jbf_plane + jsize_plane 

!   bndy%ilo(2)  = isize_plane + ibf_plane
!   bndy%ihi(2)  = isize_plane + 1
!   bndy%jlo(2)  = jbf_plane + 1 
!   bndy%jhi(2)  = jbf_plane + jsize_plane 

!   bndy%bc_num = 1
!   allocate(bndy%bc_infos(bndy%bc_num))
!   bc_info => bndy%bc_infos(1)
!   call create(bc_info,2)
!   bc_info%p(1) = blast_bc_idx
!   bc_info%p(2) = rflt3_blst_idx
    
!   ! top boundary
!   bndy => blks(1)%bndys(4)
!   bndy%idx    = 4

!   bndy%ilo(1)  = ibf_plane + 1
!   bndy%ihi(1)  = ibf_plane + isize_plane 
!   bndy%jlo(1)  = jsize_plane + jbf_plane + 1
!   bndy%jhi(1)  = jsize_plane + jbf_plane + jbf_plane

!   bndy%ilo(2)  = ibf_plane + 1
!   bndy%ihi(2)  = ibf_plane + isize_plane 
!   bndy%jlo(2)  = jsize_plane + jbf_plane
!   bndy%jhi(2)  = jsize_plane + 1

!   bndy%bc_num = 1
!   allocate(bndy%bc_infos(bndy%bc_num))
!   bc_info => bndy%bc_infos(1)
!   call create(bc_info,2)
!   bc_info%p(1) = blast_bc_idx
!   bc_info%p(2) = rflt4_blst_idx

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

!   ! Setup boundary condition 1 data
!   bc => flwd%bcs(1)
!   bc%vs_num = 1
!   allocate(bc%vs(bc%vs_num))
!   call create(bc%vs(1), 4)
!   bc%vs(1)%p(:) = 1.d0
!   bc%vs(1)%p(2) =-1.d0

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
!   bc%vs(1)%p(:) = 1.d0
!   bc%vs(1)%p(2) =-1.d0

!   ! Setup boundary condition 4 data
!   bc => flwd%bcs(4)
!   bc%vs_num = 1
!   allocate(bc%vs(bc%vs_num))
!   call create(bc%vs(1), 4)
!   bc%vs(1)%p(:) = 1.d0
!   bc%vs(1)%p(3) =-1.d0

!   info%gamma = gamma
!   do j=1,info%jsize
!     do i=1,info%isize
!       i0 = i + info%ibf
!       j0 = j + info%jbf

!       xc0 = blk%xc%p(i0,j0)
!       yc0 = blk%yc%p(i0,j0)

!       if(xc0.lt.(0.1d0)) then
!         rho0 = 1.d0
!         u0   = 0.d0
!         v0   = 0.d0
!         p0   = 1000.d0
!       else if((xc0.ge.0.1d0).and. &
!               (xc0.le.0.9d0)) then
!         rho0 = 1.d0
!         u0   = 0.d0
!         v0   = 0.d0
!         p0   = 0.1
!       else
!         rho0 = 1.d0
!         u0   = 0.d0
!         v0   = 0.d0
!         p0   = 100.d0
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
! end subroutine initial_model_eu_blast
end submodule shocktube