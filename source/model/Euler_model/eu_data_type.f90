module eu_data_type
  use IO
  use data_type
  use boundary
  use mesh
  implicit none

  character(len=64),dimension(6),parameter :: prim_name =(/"rho","u","v","p","c","flag"/)
  type :: eu_info_data_t
    ! General data information
    integer :: idx
    integer :: isize,jsize
    integer :: ibf,jbf
    double precision :: dxi,deta
    double precision :: eig_xi,eig_eta
    double precision :: eig_xi_g,eig_eta_g
    double precision :: dt
    double precision :: time

    double precision :: gamma
    logical :: if_vis
    logical :: if_src

    integer :: stage_cur
    integer :: stage_num

    integer :: cons_num
    integer :: cnvt_num
    integer :: dfsv_num
    integer :: src_num
    integer :: prim_num  
  end type

  type :: eu_convect_data_t
    type(rvectorfield_t) :: cnvt  

    !fields related to convective terms
    type(rvectorfield_t) :: flux_xi
    type(rvectorfield_t) :: flux_xip
    type(rvectorfield_t) :: flux_xin
    type(rvectorfield_t) :: h_xip
    type(rvectorfield_t) :: h_xin
    type(rvectorfield_t) :: diff_xip
    type(rvectorfield_t) :: diff_xin

    type(rvectorfield_t) :: flux_eta
    type(rvectorfield_t) :: flux_etap
    type(rvectorfield_t) :: flux_etan
    type(rvectorfield_t) :: h_etan
    type(rvectorfield_t) :: h_etap 
    type(rvectorfield_t) :: diff_etan
    type(rvectorfield_t) :: diff_etap 
  end type

  type :: eu_diffusive_data_t
    type(rvectorfield_t) :: dfsv
    type(rscalarfield_t),dimension(:),allocatable   :: temp
    type(rscalarfield_t),dimension(:),allocatable   :: temp_grad

    type(rscalarfield_t)                            :: mu
    type(rscalarfield_t),dimension(:,:),allocatable :: uv
    type(rscalarfield_t),dimension(:,:),allocatable :: uv_grad
    
    type(rvectorfield_t),dimension(:),allocatable   :: stress
    type(rvectorfield_t),dimension(:),allocatable   :: stress_grad
  end type

  type :: eu_source_data_t
    type(rvectorfield_t) :: src  
  end type

  ! Assembly data of the Euler model
  TYPE :: eu_data_t
    type(eu_info_data_t) :: info
    type(rvectorfield_t) :: prim
    type(rvectorfield_t),dimension(:),allocatable :: cons 
    type(rvectorfield_t),dimension(:),allocatable :: lhs
    type(rvectorfield_t),dimension(:),allocatable :: rhs

    type(eu_convect_data_t)   :: cnvt
    type(eu_diffusive_data_t) :: dfsv
    type(eu_source_data_t)    :: src  

    type(block_t),pointer     :: blk
    
    logical :: partisioned
    class(eu_data_t),dimension(:,:),pointer :: tiles
    class(eu_data_t), pointer               :: parent
  contains
    procedure, public :: config_eu
    procedure, public :: create_eu
    procedure, public :: set_stage
    procedure, public :: write_data_fmt
    final             :: destroy_eu
  end TYPE  

  type :: eu_data_list_t
    class(eu_data_t), pointer :: flwd
    type(eu_data_list_t), pointer :: prev, next
  contains
    final :: empty_flwd_list
  end type eu_data_list_t

  type :: eu_data_iterator_t
    class(eu_data_t),pointer :: flwd
  end type eu_data_iterator_t

contains
subroutine config_eu(eudata, stage_num, blk)
  implicit none
  integer,intent(in) :: stage_num
  class(block_t),target,intent(in) :: blk
  class(eu_data_t),intent(inout) :: eudata
  
  eudata%info%stage_num = stage_num
  eudata%info%prim_num  = 6
  eudata%info%cons_num  = 4
  eudata%info%cnvt_num  = 4
  eudata%info%dfsv_num  = 4
  eudata%info%src_num   = 4

  eudata%blk         => blk
  eudata%info%idx    = blk%idx
  eudata%info%isize  = blk%data%isize
  eudata%info%jsize  = blk%data%jsize
  eudata%info%ibf    = blk%data%ibf
  eudata%info%jbf    = blk%data%jbf
  eudata%info%dxi    = blk%data%dxi
  eudata%info%deta   = blk%data%deta
  eudata%tiles       => null()
  eudata%parent      => null()

end subroutine config_eu

subroutine create_eu(eudata)
  implicit none  
  class(eu_data_t),intent(inout) :: eudata
  type(fvectorshape_t) :: cons_s
  type(fvectorshape_t) :: prim_s
  integer :: istage
  integer :: stage_num
  integer :: prim_num,cons_num
  integer :: cnvt_num
  integer :: dfsv_num
  integer :: src_num
  integer :: isize, jsize ,ibf, jbf
  integer :: ibc
  
  isize     = eudata%info%isize
  jsize     = eudata%info%jsize
  ibf       = eudata%info%ibf
  jbf       = eudata%info%jbf

  prim_num  = eudata%info%prim_num
  stage_num = eudata%info%stage_num
  cons_num  = eudata%info%cons_num
  cnvt_num  = eudata%info%cnvt_num
  dfsv_num  = eudata%info%dfsv_num
  src_num   = eudata%info%src_num

  ! Allocate primative variables data
  prim_s = set_shape(prim_num,       &
                     isize + 2*ibf,  &
                     jsize + 2*jbf)
  call create(eudata%prim,prim_s)
  eudata%prim%p = 0.d0

  ! Allocate conservative variables data
  allocate(eudata%cons(stage_num))
  cons_s = set_shape(cons_num,       &
                     isize + 2*ibf,  &
                     jsize + 2*jbf)
  do istage = 1, stage_num
    call create(eudata%cons(istage),cons_s)
    eudata%cons(istage)%p = 0.d0
  enddo

  ! Allocate RHS data
  allocate(eudata%rhs(stage_num))
  do istage = 1, stage_num
    call create(eudata%rhs(istage),cons_s)
    eudata%rhs(istage)%p = 0.d0
  enddo

  ! Allocate LHS data
  allocate(eudata%lhs(stage_num))
  do istage = 1, stage_num
    call create(eudata%lhs(istage),cons_s)
    eudata%lhs(istage)%p = 0.d0
  enddo

  call create_cnvt_data(eudata%cnvt,cnvt_num,isize,jsize,ibf,jbf)
  call create_dfsv_data(eudata%dfsv,dfsv_num,isize,jsize,ibf,jbf)
  call create_src_data(eudata%src,src_num,isize,jsize,ibf,jbf)
end subroutine create_eu

subroutine set_stage(eu_data,stage)
  implicit none
  class(eu_data_t),intent(inout) :: eu_data
  integer,intent(in) :: stage
  eu_data%info%stage_cur = stage
end subroutine set_stage

subroutine create_cnvt_data(cnvt, &
                            cnvt_num, &   
                            isize,jsize,ibf,jbf)
  implicit none
  integer,intent(in) :: cnvt_num
  integer,intent(in) :: isize,jsize,ibf,jbf
  type(eu_convect_data_t),intent(inout) :: cnvt
  integer :: icnvt
  type(fvectorshape_t) :: cnvt_s
  type(fvectorshape_t) :: flux_xi_s
  type(fvectorshape_t) :: flux_eta_s
  type(fvectorshape_t) :: h_xi_s
  type(fvectorshape_t) :: h_eta_s
  type(fvectorshape_t) :: diff_xi_s
  type(fvectorshape_t) :: diff_eta_s
  if(cnvt_num.ne.0) then
    ! 1. Allocate convect term data
    if(allocated(cnvt%cnvt%p)) then
      write(*,*) "Error: create_cnvt_data()"
      stop 'pointer::cnvt%cnvt has already been allocated.'
    endif
    cnvt_s = set_shape(cnvt_num,      &
                       isize + 2*ibf, &
                       jsize + 2*jbf)
    call create(cnvt%cnvt,cnvt_s)
    cnvt%cnvt%p = 0.d0

    ! 2.Allocate xi direction flux data of the convect term
    if(allocated(cnvt%flux_xi%p)) then
      write(*,*) "Error: create_cnvt_data()"
      stop 'pointer::cnvt%flux_xi has already been allocated.'
    endif
    flux_xi_s = set_shape(cnvt_num,      &
                          isize + 2*ibf, &
                          jsize)
    call create(cnvt%flux_xi,flux_xi_s)
    cnvt%flux_xi%p = 0.d0

    ! 2.1 Allocate positive/negative splitted xi-direction flux data
    if(allocated(cnvt%flux_xip%p)) then
      write(*,*) "Error: create_cnvt_data()"
      stop 'pointer::cnvt%flux_xip has already been allocated.'
    endif
    call create(cnvt%flux_xip,flux_xi_s)
    cnvt%flux_xip%p = 0.d0

    if(allocated(cnvt%flux_xin%p)) then
      write(*,*) "Error: create_cnvt_data()"
      stop 'pointer::cnvt%flux_xin has already been allocated.'
    endif
    call create(cnvt%flux_xin,flux_xi_s)
    cnvt%flux_xin%p = 0.d0

    ! 2.2 Allocate positive/negative reconstructed xi-direction flux data
    h_xi_s = set_shape(cnvt_num, isize + 1, jsize)
    if(allocated(cnvt%h_xip%p)) then
      write(*,*) "Error: create_cnvt_data()"
      stop 'pointer::cnvt%h_xip has already been allocated.'
    endif
    call create(cnvt%h_xip,h_xi_s)
    cnvt%h_xip%p = 0.d0
    if(allocated(cnvt%h_xin%p)) then
      write(*,*) "Error: create_cnvt_data()"
      stop 'pointer::cnvt%h_xin has already been allocated.'
    endif
    call create(cnvt%h_xin,h_xi_s)
    cnvt%h_xin%p = 0.d0

    ! 2.3 Allocate positive/negative xi-direction derivative data
    diff_xi_s = set_shape(cnvt_num, isize, jsize)
    if(allocated(cnvt%diff_xip%p)) then
      write(*,*) "Error: create_cnvt_data()"
      stop 'pointer::cnvt%diff_xip has already been allocated.'
    endif
    call create(cnvt%diff_xip,diff_xi_s)
    cnvt%diff_xip%p = 0.d0

    if(allocated(cnvt%diff_xin%p)) then
      write(*,*) "Error: create_cnvt_data()"
      stop 'pointer::cnvt%diff_xin has already been allocated.'
    endif
    call create(cnvt%diff_xin,diff_xi_s)
    cnvt%diff_xin%p = 0.d0


    ! 3. Allocate eta direction flux data of the convect term
    if(allocated(cnvt%flux_eta%p)) then
      write(*,*) "Error: create_cnvt_data()"
      stop 'pointer::cnvt%flux_eta has already been allocated.'
    endif
    flux_eta_s = set_shape(cnvt_num,      &
                           jsize + 2*jbf, &
                           isize)
    call create(cnvt%flux_eta,flux_eta_s)
    cnvt%flux_eta%p = 0.d0

    ! 3.1 Allocate positive/negative splitted eta-direction flux data
    flux_eta_s = set_shape(cnvt_num, jsize + 2*jbf, isize)
    if(allocated(cnvt%flux_etap%p)) then
      write(*,*) "Error: create_cnvt_data()"
      stop 'pointer::cnvt%flux_etap has already been allocated.'
    endif
    call create(cnvt%flux_etap,flux_eta_s)
    cnvt%flux_etap%p = 0.d0

    if(allocated(cnvt%flux_etan%p)) then
      write(*,*) "Error: create_cnvt_data()"
      stop 'pointer::cnvt%flux_etan has already been allocated.'
    endif     
    call create(cnvt%flux_etan,flux_eta_s)
    cnvt%flux_etan%p = 0.d0

    ! 3.2 Allocate positive/negative reconstructed eta-direction flux data
    h_eta_s = set_shape(cnvt_num, jsize + 1, isize)
    if(allocated(cnvt%h_etap%p)) then
      write(*,*) "Error: create_cnvt_data()"
      stop 'pointer::cnvt%h_etap has already been allocated.'
    endif
    call create(cnvt%h_etap,h_eta_s)
    cnvt%h_etap%p = 0.d0
    if(allocated(cnvt%h_etan%p)) then
      write(*,*) "Error: create_cnvt_data()"
      stop 'pointer::cnvt%h_etan has already been allocated.'
    endif
    call create(cnvt%h_etan,h_eta_s)
    cnvt%h_etan%p = 0.d0

    ! 3.3 Allocate positive/negative eta-direction derivative data
    diff_eta_s = set_shape(cnvt_num, jsize, isize)
    if(allocated(cnvt%diff_etap%p)) then
      write(*,*) "Error: create_cnvt_data()"
      stop 'pointer::cnvt%diff_etap has already been allocated.'
    endif
    call create(cnvt%diff_etap,diff_eta_s)
    cnvt%diff_etap%p = 0.d0

    if(allocated(cnvt%diff_etan%p)) then
      write(*,*) "Error: create_cnvt_data()"
      stop 'pointer::cnvt%diff_etan has already been allocated.'
    endif      
    call create(cnvt%diff_etan,diff_eta_s)
    cnvt%diff_etan%p = 0.d0
  endif
end subroutine create_cnvt_data

subroutine create_dfsv_data(dfsv,      &
                            dfsv_num,  &    
                            isize,jsize,ibf,jbf)
  implicit none
  integer,intent(in) :: dfsv_num
  integer,intent(in) :: isize,jsize,ibf,jbf
  type(eu_diffusive_data_t),intent(inout) :: dfsv
  type(fvectorshape_t) :: dfsv_s
  type(fscalarshape_t) :: temp_xi_s,temp_eta_s
  type(fscalarshape_t) :: mu_s
  type(fscalarshape_t) :: uv_xi_s,uv_eta_s
  type(fscalarshape_t) :: uv_grad_xi_s,uv_grad_eta_s
  type(fvectorshape_t) :: stress_xi_s,stress_eta_s
  type(fvectorshape_t) :: stress_grad_xi_s, stress_grad_eta_s  

  ! if(dfsv_num.ne.0) then
  !   if(allocated(dfsv%dfsv%p)) then
  !     write(*,*) "Error: create_dfsv_data()"
  !     stop 'pointer::dfsv%dfsv has already been allocated.'
  !   endif      
  !   dfsv_s = set_shape(dfsv_num, isize, jsize)
  !   call create(dfsv%dfsv,dfsv_s)
  !   dfsv%dfsv%p = 0.d0

  !   allocate(dfsv%temp(2))
  !   temp_xi_s = set_shape(isize + 2*ibf, jsize)
  !   call create(dfsv%temp(1),temp_xi_s)
  !   dfsv%temp(1)%p = 0.d0
  !   temp_eta_s = set_shape(jsize + 2*jbf, isize)
  !   call create(dfsv%temp(2),temp_eta_s)
  !   dfsv%temp(2)%p = 0.d0
    
  !   allocate(dfsv%temp_grad(2))
  !   temp_grad_xi_s = set_shape(isize + 2*ibf, jsize)
  !   call create(dfsv%temp_grad(1),temp_grad_xi_s)
  !   dfsv%temp_grad(1)%p = 0.d0
  !   temp_grad_eta_s = set_shape(jsize + 2*jbf, isize)
  !   call create(dfsv%temp_grad(2),temp_grad_eta_s)
  !   dfsv%temp_grad(2)%p = 0.d0
    
  !   mu_s = set_shape(isize + 2*ibf, jsize + 2*jbf)
  !   call create(dfsv%mu,mu_s)
  !   dfsv%mu%p = 0.d0

  !   allocate(dfsv%uv(2,2))
  !   uv_xi_s  = set_shape(isize + 2*ibf, jsize)
  !   uv_eta_s = set_shape(jsize + 2*jbf, isize)
  !   call create(dfsv%uv(1,1),uv_xi_s)
  !   dfsv%uv(1,1)%p = 0.d0
  !   call create(dfsv%uv(2,1),uv_xi_s)
  !   dfsv%uv(2,1)%p = 0.d0
  !   call create(dfsv%uv(1,2),uv_eta_s)
  !   dfsv%uv(1,2)%p = 0.d0
  !   call create(dfsv%uv(2,2),uv_eta_s)
  !   dfsv%uv(2,2)%p = 0.d0

  !   allocate(dfsv%uv_grad(2,2))
  !   uv_grad_xi_s  = set_shape(isize + 2*ibf, jsize)
  !   uv_grad_eta_s = set_shape(jsize + 2*jbf, isize)
  !   call create(dfsv%uv_grad(1,1),uv_grad_xi_s)
  !   dfsv%uv_grad(1,1)%p = 0.d0
  !   call create(dfsv%uv_grad(2,1),uv_grad_xi_s)
  !   dfsv%uv_grad(2,1)%p = 0.d0
  !   call create(dfsv%uv_grad(1,2),uv_grad_eta_s)
  !   dfsv%uv_grad(1,2)%p = 0.d0
  !   call create(dfsv%uv_grad(2,2),uv_grad_eta_s)
  !   dfsv%uv_grad(2,2)%p = 0.d0

  !   allocate(dfsv%stress(2))
  !   stress_xi_s  = set_shape(dfsv_num, isize + 2*ibf, jsize)
  !   stress_eta_s = set_shape(dfsv_num, jsize + 2*jbf, isize)
  !   call create(dfsv%stress(1),stress_xi_s)
  !   dfsv%stress(1)%p = 0.d0
  !   call create(dfsv%stress(2),stress_eta_s)
  !   dfsv%stress(2)%p = 0.d0

  !   allocate(dfsv%stress_grad(2))
  !   stress_grad_xi_s  = set_shape(dfsv_num, isize + 2*ibf, jsize)
  !   stress_grad_eta_s = set_shape(dfsv_num, jsize + 2*jbf, isize)
  !   call create(dfsv%stress_grad(1),stress_grad_xi_s)
  !   dfsv%stress_grad(1)%p = 0.d0
  !   call create(dfsv%stress_grad(2),stress_grad_eta_s)
  !   dfsv%stress_grad(2)%p = 0.d0
  ! endif
end subroutine create_dfsv_data

subroutine create_src_data(src,       &    
                           src_num,   &      
                           isize,jsize,ibf,jbf)
  implicit none
  integer,intent(in) :: src_num
  integer,intent(in) :: isize,jsize,ibf,jbf
  type(eu_source_data_t),intent(inout) :: src
  integer :: isrc
  type(fvectorshape_t) :: src_s
  type(fscalarshape_t) :: z_xi_s,z_eta_s,hz_xi_s,hz_eta_s
  if(src_num.ne.0) then
    if(allocated(src%src%p)) then
      write(*,*) "Error: create_src_data()"
      stop 'pointer::src%src has already been allocated.'
    endif
    src_s = set_shape(src_num,       &
                      isize + 2*ibf, &
                      jsize + 2*jbf)
    call create(src%src,src_s)
    src%src%p = 0.d0
  endif
end subroutine create_src_data

subroutine destroy_eu(eudata)
  implicit none  
  type(eu_data_t),intent(inout) :: eudata
  integer :: ibc, istage
    
  if(associated(eudata%tiles)) then  
    deallocate(eudata%tiles)
  endif

  call destroy(eudata%prim)

  do istage = 1, eudata%info%stage_num    
    call destroy(eudata%cons(istage))
  enddo
  deallocate(eudata%cons)

  do istage = 1, eudata%info%stage_num    
    call destroy(eudata%lhs(istage))
  enddo
  deallocate(eudata%lhs)

  do istage = 1, eudata%info%stage_num    
    call destroy(eudata%rhs(istage))
  enddo
  deallocate(eudata%rhs)

  call destroy_cnvt(eudata%cnvt,eudata%info%cnvt_num)
  call destroy_dfsv(eudata%dfsv,eudata%info%dfsv_num)
  call destroy_src(eudata%src,eudata%info%src_num)

end subroutine destroy_eu

subroutine destroy_cnvt(cnvt, cnvt_num)
  implicit none
  integer,intent(in) :: cnvt_num
  type(eu_convect_data_t),intent(inout) :: cnvt
  integer :: icnvt
  if(cnvt_num.ne.0) then
    call destroy(cnvt%cnvt)

    call destroy(cnvt%flux_xi)
    call destroy(cnvt%flux_xip)
    call destroy(cnvt%flux_xin)
    call destroy(cnvt%h_xip)
    call destroy(cnvt%h_xin)
    call destroy(cnvt%diff_xip)
    call destroy(cnvt%diff_xin)

    call destroy(cnvt%flux_eta)
    call destroy(cnvt%flux_etap)
    call destroy(cnvt%flux_etan)
    call destroy(cnvt%h_etap)
    call destroy(cnvt%h_etan)
    call destroy(cnvt%diff_etap)
    call destroy(cnvt%diff_etan)
  endif
end subroutine destroy_cnvt

subroutine destroy_dfsv(dfsv,dfsv_num)
  implicit none
  integer,intent(in) :: dfsv_num
  type(eu_diffusive_data_t),intent(inout) :: dfsv
end subroutine destroy_dfsv

subroutine destroy_src(src,src_num)
  implicit none
  integer,intent(in) :: src_num
  type(eu_source_data_t),intent(inout) :: src
  if(src_num.ne.0) then
    call destroy(src%src)
  endif
end subroutine destroy_src

recursive subroutine empty_flwd_list(flwd_list)
  implicit none
  type(eu_data_list_t), intent(inout) :: flwd_list
  if(associated(flwd_list%next)) then
    deallocate(flwd_list%next)
  endif
  if(associated(flwd_list%flwd)) then
    deallocate(flwd_list%flwd)
  endif
end subroutine empty_flwd_list

subroutine write_data_fmt(eudata, funit)
  implicit none
  class(eu_data_t),intent(in) :: eudata
  integer,intent(in) :: funit
  
  type(fvectorshape_t) :: lo_prm,hi_prm
  type(fvectorshape_t) :: lo_phs,hi_phs
  type(fvectorshape_t) :: lo_cpnt,hi_cpnt
  integer :: isize,jsize,ibf,jbf
  integer :: prim_num

  isize = eudata%info%isize
  jsize = eudata%info%jsize
  ibf   = eudata%info%ibf
  jbf   = eudata%info%jbf

  prim_num = eudata%info%prim_num
  lo_prm  = set_shape(1,1+ibf,1+jbf)
  hi_prm  = set_shape(prim_num,isize+ibf,jsize+jbf)  
  call write_field_fmt(funit,                 & 
                       eudata%prim,           &
                       lo_prm,                &
                       hi_prm)
  
end subroutine write_data_fmt

subroutine write_data_hdf5(eudata, file_id)
  use HDF5
  implicit none
  class(eu_data_t),intent(in) :: eudata
  integer(HID_T),intent(in)   :: file_id

  type(fvectorshape_t) :: lo_prm,hi_prm
  integer :: isize,jsize,ibf,jbf
  integer :: prim_num

  isize = eudata%info%isize
  jsize = eudata%info%jsize
  ibf   = eudata%info%ibf
  jbf   = eudata%info%jbf

  prim_num = eudata%info%prim_num
  lo_prm  = set_shape(1,1+ibf,1+jbf)
  hi_prm  = set_shape(prim_num,isize+ibf,jsize+jbf)  
  call write_field_hdf5(eudata%prim,   &
                        lo_prm,        &
                        hi_prm,        &
                        prim_name,     &
                        file_id)
  
end subroutine write_data_hdf5
end module eu_data_type