module mesh
  use data_type
  use boundary
  implicit none

  type :: bndy_data_t
    ! index range for each boundary
    ! (1) : local boundary range
    ! (2) : donner data range
    integer,dimension(2) :: ilo, ihi
    integer,dimension(2) :: jlo, jhi

    ! boundary condition type
    type(ivector_t) :: bc_infos

    ! Boundary condition specific values
    integer :: scalar_num = 0
    integer :: vector_num = 0
    integer :: matrix_num = 0
    integer :: tensor_num = 0
    double precision,dimension(:),allocatable :: scalar
    type(rvector_t),dimension(:),allocatable  :: vector
    type(rmatrix_t),dimension(:),allocatable  :: matrix
    type(rtensor_t),dimension(:),allocatable  :: tensor
  end type bndy_data_t

  type :: bndy_t
    integer :: idx
    integer :: group_id ! 1,2,3,4 indicates edge position
    type(bndy_data_t),pointer :: data
  contains
    procedure :: deallocate_d_bndy
  end type bndy_t

  type :: bndy_list_t
    class(bndy_t), pointer :: bndy
    type(bndy_list_t), pointer :: prev, next
  contains
    final :: empty_bndy_list
  end type bndy_list_t

  type :: bndy_iterator_t
    class(bndy_t),pointer :: bndy
  end type bndy_iterator_t

  ! Cartisian mesh block
  TYPE :: block_data_t
    ! Block sizes and buffer sizes
    integer :: ibf, jbf
    integer :: isize, jsize

    ! vortex coord
    type(rscalarfield_t) :: xv, yv

    ! Cell center coord 
    type(rscalarfield_t) :: xc, yc

    ! Face center and face norm of each cell
    ! type(rmatrixfield_t) :: f_cntr
    ! type(rmatrixfield_t) :: f_norm

    double precision :: dxi, deta
  contains
    procedure, private :: allocate_d_blk
    procedure, private :: deallocate_d_blk
    procedure, private :: initial_d_blk
    procedure, public  :: cmpt_cell_geo
    procedure, public  :: cmpt_face_geo
    procedure, public  :: cmpt_metrics
  end TYPE block_data_t

  ! Cartisian mesh block
  TYPE :: block_t    
    integer :: idx         ! global unique index
    class(block_data_t),pointer :: data

    integer :: bndy_num
    integer :: bndy_count 
    type(bndy_list_t),pointer :: bndy_list, &
                                 bndy_list_tail
    type(bndy_iterator_t),dimension(:),allocatable :: bndy_que

    ! Partision parameters
    logical :: partisioned
    ! As root block
    integer :: tile_dim(2)
    integer :: tile_size(2)
    type(block_t),pointer :: tile_head, tile_tail
    type(block_t),dimension(:,:),pointer :: tiles
  
    ! As tile block
    integer :: tile_idx(2) ! local index
    type(block_t),pointer :: parent_tile
    type(block_t),pointer :: prev_tile, next_tile
    type(block_t),pointer :: east, south, west, north, &
                             southeast, southwest, northwest, northeast 
    ! parent block index covered by present tile
    integer :: p_ilo, p_ihi, p_jlo, p_jhi     
    contains
    procedure, public :: create_block
    procedure, public :: initial_block
    procedure, public :: write_block_ascii
    procedure, public :: destroy_block
    procedure, public :: split_block_to_tiles
    procedure, public :: split_bndy_to_tiles
    procedure, public :: reset_block_tiles
    procedure, public :: add_bndy
    procedure, public :: get_bndy
    procedure, public :: rmv_bndy
  end TYPE block_t

  type :: block_list_t
    class(block_t), pointer :: blk
    type(block_list_t), pointer :: prev, next
  contains
    final     :: empty_blk_list
  end type block_list_t

  type :: block_iterator_t
    class(block_t),pointer :: blk
  end type block_iterator_t
  
  integer :: base_blk_num = 0
  ! Store base blocks
  class(block_list_t), pointer :: base_blk_list, &
                                  base_blk_list_tail

  ! Store sorted partisioned blocks
  class(block_list_t), pointer :: sorted_blk
  
  ! Store all blocks needs to be looped
  integer :: blk_count
  class(block_iterator_t),dimension(:),allocatable :: blk_que

  logical :: if_mesh_partisioned = .false.
  integer :: g_tile_size(2) = (/64,64/)
  NAMELIST /mesh_params/ g_tile_size

  integer,parameter :: left_bndy_id   = 1, &
                       bottom_bndy_id = 2, &
                       right_bndy_id  = 3, &
                       top_bndy_id    = 4

contains
subroutine create_mesh(funit)
  implicit none
  integer,intent(in) :: funit
  read(funit,NML=mesh_params)
end subroutine create_mesh

function add_blk_to_mesh() result(blk_node)
  implicit none
  class(block_list_t),pointer :: blk_node  
  class(block_list_t),pointer :: curr_node
  allocate(blk_node)
  allocate(blk_node%blk)
  blk_node%prev  => null()
  blk_node%next => null()
  if(.not.associated(base_blk_list)) then 
    base_blk_list      => blk_node
    base_blk_list_tail => base_blk_list
  else
    base_blk_list_tail%next => blk_node
    blk_node%prev            => base_blk_list_tail
    base_blk_list_tail      => blk_node
  endif
end function add_blk_to_mesh

subroutine partision_mesh()
  implicit none
  print*,"Start partisioning blocks..."
  if(allocated(blk_que)) deallocate(blk_que)
  call reset_blk_in_list(base_blk_list)
  call split_blk_in_list(base_blk_list,blk_count,g_tile_size)
  allocate(blk_que(blk_count))
  call assemble_blk_que_in_list(base_blk_list,blk_que,blk_count)
  print*,"Partision finished!"
end subroutine partision_mesh

subroutine destroy_mesh()
  implicit none
  if(allocated(blk_que)) deallocate(blk_que)
  blk_count          = 0
  deallocate(base_blk_list)
  base_blk_num       = 0
  base_blk_list      => null()
  base_blk_list_tail => null()
end subroutine destroy_mesh

subroutine reset_blk_in_list(blk_list)
  implicit none
  class(block_list_t),pointer :: blk_list
  class(block_list_t),pointer :: curr_node
  class(block_t),pointer :: blk,tile
  integer :: itile,jtile
  curr_node => blk_list
  do while(associated(curr_node))
    blk => curr_node%blk
    call blk%reset_block_tiles()
    curr_node => curr_node%next
  enddo
end subroutine reset_blk_in_list

subroutine split_blk_in_list(blk_list,tile_num,tile_size)
  implicit none
  class(block_list_t),pointer :: blk_list
  integer,intent(inout) :: tile_num
  integer,intent(in) :: tile_size(2)
  class(block_list_t),pointer :: curr
  class(block_t),pointer :: blk
  tile_num = 0  
  !Split block to tiles
  curr => blk_list
  do while(associated(curr))
    blk => curr%blk
    if(blk%partisioned) then
      print*,"Error : split_blk_in_list()"
      print*,"blk",blk%idx,"has been splitted!"
      stop
    else
      call blk%split_block_to_tiles(tile_size,tile_num)
      tile_num = tile_num + blk%tile_dim(1)*blk%tile_dim(2)
    endif
    curr => curr%next
  enddo

  !Set tile boundaries
  curr => blk_list
  do while(associated(curr))
    blk => curr%blk
    if(.not.blk%partisioned) then
      print*,"Error : "
      print*,"blk",blk%idx,"has not been splitted!"
      stop
    else
      call blk%split_bndy_to_tiles()
    endif
    curr => curr%next
  enddo
end subroutine split_blk_in_list

subroutine assemble_blk_que_in_list(blk_list,blk_que,que_size)
  implicit none
  class(block_list_t),pointer :: blk_list
  integer,intent(in) :: que_size
  class(block_iterator_t),dimension(que_size),intent(inout) :: blk_que
  class(block_list_t),pointer :: curr
  class(block_t),pointer :: blk
  class(block_t),pointer :: tile
  integer :: itile,jtile  

  curr => blk_list
  do while(associated(curr))
    blk => curr%blk
    if(blk%partisioned) then
      do jtile = 1, blk%tile_dim(2)
        do itile = 1, blk%tile_dim(1)
          tile => curr%blk%tiles(itile, jtile)
          blk_que(tile%idx)%blk => tile
        enddo
      enddo
    else
      print*,"Error : assemble_blk_que_in_list()"
      print*,"blk",blk%idx,"has not been splitted!"
      stop
    endif
    curr => curr%next
  enddo
end subroutine assemble_blk_que_in_list

function push_back_blk(blk_list) result(blk_node)
  implicit none
  class(block_list_t),pointer,intent(inout) :: blk_list 
  class(block_list_t),pointer :: blk_node  
  class(block_list_t),pointer :: curr_node
  allocate(blk_node)
  allocate(blk_node%blk)
  blk_node%prev => null()
  blk_node%next => null()
  if(.not.associated(blk_list)) then 
    blk_list => blk_node
  else
    curr_node     => blk_list
    do while(associated(curr_node%next))
      curr_node => curr_node%next
    enddo
    curr_node%next => blk_node
    blk_node%prev  => curr_node
  endif
end function push_back_blk

function get_blk(blk_list,blk_idx) result(blk_node)
  implicit none
  class(block_list_t),target,intent(inout) :: blk_list 
  integer,intent(in) :: blk_idx 
  class(block_list_t),pointer :: blk_node
  blk_node => blk_list
  do while(associated(blk_node))
    if(blk_node%blk%idx.eq.blk_idx) return
    blk_node => blk_node%next
  enddo
end function get_blk

subroutine remove_blk(blk)
  implicit none
  class(block_t),target,intent(inout) :: blk 
  class(block_list_t),pointer :: curr,next,prev  

end subroutine remove_blk

recursive subroutine empty_blk_list(blk_list)
  implicit none
  type(block_list_t), intent(inout) :: blk_list
  if(associated(blk_list%next)) then
    deallocate(blk_list%next)
  endif
  call blk_list%blk%destroy_block()
  deallocate(blk_list%blk)
end subroutine empty_blk_list

subroutine create_block(blk,isize,jsize,ibf,jbf,idx)
  implicit none
  integer,intent(in) :: isize,jsize,ibf,jbf
  integer,intent(in),optional :: idx
  class(block_t),target,intent(inout) :: blk

  if(present(idx)) then
    blk%idx = idx
  else
    blk%idx = 0
  endif

  allocate(blk%data)
  call blk%data%allocate_d_blk(isize,jsize, &
                               ibf,jbf)
   
  blk%bndy_num   =  0
  blk%bndy_count =  0
  blk%bndy_list  => null()
  
  blk%partisioned = .false. 
  blk%tile_dim    = 0
  blk%tile_size   = 0
  blk%tile_head   => null()
  blk%tile_tail   => null()
  blk%tiles       => null()

  blk%tile_idx    = 0
  blk%parent_tile => null()
  blk%prev_tile   => null()
  blk%next_tile   => null()
  blk%east        => null()
  blk%south       => null()
  blk%west        => null()
  blk%north       => null()
  blk%southeast   => null()
  blk%southwest   => null()
  blk%northwest   => null()
  blk%northeast   => null() 

  blk%p_ilo = 0
  blk%p_ihi = 0
  blk%p_jlo = 0
  blk%p_jhi = 0     
end subroutine create_block

subroutine initial_block(blk,x0,y0,xl,yl)
  implicit none
  double precision,intent(in) :: x0,y0,xl,yl
  class(block_t),target,intent(inout) :: blk
  class(block_data_t),pointer :: d_blk

  d_blk => blk%data
  call d_blk%initial_d_blk(x0,y0,xl,yl)
  call d_blk%cmpt_cell_geo()
  call d_blk%cmpt_face_geo()
  call d_blk%cmpt_metrics()
end subroutine initial_block

subroutine write_block_ascii(blk,funit)
  implicit none
  class(block_t),target,intent(inout) :: blk
  integer,intent(in) :: funit
  type(block_data_t),pointer :: d_blk
  integer :: i,j

  d_blk => blk%data

  do j= d_blk%jbf + 1, d_blk%jbf + d_blk%jsize + 1
    do i= d_blk%ibf + 1, d_blk%ibf + d_blk%isize + 1
      write(funit,*) d_blk%xv%p(i,j)
    end do
  end do
  
  do j= d_blk%jbf + 1, d_blk%jbf + d_blk%jsize + 1
    do i= d_blk%ibf + 1, d_blk%ibf + d_blk%isize + 1
      write(funit,*) d_blk%yv%p(i,j)
    end do
  end do  
end subroutine write_block_ascii

recursive subroutine destroy_block(blk)
  implicit none
  class(block_t),target,intent(inout) :: blk
  class(block_t),pointer :: tile
  integer :: itile,jtile

  if(blk%partisioned) then
    do jtile = 1, blk%tile_dim(2)
      do itile = 1, blk%tile_dim(1)
        tile => blk%tiles(itile,jtile)
        call tile%destroy_block()        
      enddo
    enddo
    deallocate(blk%tiles)
  endif

  call blk%data%deallocate_d_blk()
  deallocate(blk%data)  
  if(associated(blk%bndy_list)) deallocate(blk%bndy_list)
  blk%bndy_num   =  0
  blk%bndy_count =  0
  blk%bndy_list      => null()
  blk%bndy_list_tail => null()
  if(allocated(blk%bndy_que)) deallocate(blk%bndy_que)
  
  blk%partisioned = .false. 
  blk%tile_dim        = 0
  blk%tile_size       = 0
  blk%tile_head   => null()
  blk%tile_tail   => null()
  blk%tiles       => null()

  blk%tile_idx        = 0
  blk%parent_tile => null()
  blk%prev_tile    => null()
  blk%next_tile   => null()
  blk%east        => null()
  blk%south       => null()
  blk%west        => null()
  blk%north       => null()
  blk%southeast   => null()
  blk%southwest   => null()
  blk%northwest   => null()
  blk%northeast   => null() 

  blk%p_ilo = 0
  blk%p_ihi = 0
  blk%p_jlo = 0
  blk%p_jhi = 0     
end subroutine destroy_block

subroutine split_block_to_tiles(blk, tile_size, idx_base)
  implicit none
  integer,intent(in) :: tile_size(2)
  integer,intent(in) :: idx_base
  class(block_t),target,intent(inout) :: blk
  class(block_data_t),pointer         :: d_blk  
  class(block_t), pointer             :: tile,tile0
  class(block_data_t),pointer         :: d_tile,d_tile0 
  class(bndy_list_t), pointer         :: curr
  class(bndy_t), pointer              :: bndy,bndy0
  class(bndy_data_t), pointer         :: d_bndy,d_bndy0

  integer :: irem, jrem

  integer :: idx
  integer :: itile,jtile
  integer :: isize,jsize
  integer :: ibf,jbf
  integer :: i,j
  integer :: i0,j0

  if(blk%partisioned) then
    print*,"Error: split_block_to_tiles()."
    print*,"- Block",blk%idx,"has been splitted."
    stop
  endif

  d_blk => blk%data
  ibf = d_blk%ibf
  jbf = d_blk%jbf
  
  irem = mod(d_blk%isize,tile_size(1))
  jrem = mod(d_blk%jsize,tile_size(2))
  
  blk%tile_size = tile_size
  if((irem.le.d_blk%ibf)) then  
    blk%tile_dim(1) = max(d_blk%isize/tile_size(1),1)
    irem = irem + tile_size(1)
  else
    blk%tile_dim(1) = d_blk%isize/tile_size(1) + 1
  endif

  if((jrem.le.d_blk%jbf)) then
    blk%tile_dim(2) = max(d_blk%jsize/tile_size(2),1)
    jrem = jrem + tile_size(2)
  else
    blk%tile_dim(2) = d_blk%jsize/tile_size(2) + 1
  endif

  allocate(blk%tiles(blk%tile_dim(1),blk%tile_dim(2)))

  do jtile = 1, blk%tile_dim(2)
    do itile = 1, blk%tile_dim(1)
      tile => blk%tiles(itile,jtile)
      if(itile.eq.blk%tile_dim(1)) then
        isize = irem
      else
        isize = tile_size(1)
      endif

      if(jtile.eq.blk%tile_dim(2)) then
        jsize = jrem
      else
        jsize = tile_size(2)
      endif

      idx = idx_base + (jtile-1)*blk%tile_dim(1) + itile
      call tile%create_block(isize, jsize, &
                             ibf, jbf,     &
                             idx)
            
      tile%tile_idx(1) = itile
      tile%tile_idx(2) = jtile
      ! This is interior points range of a tile in parent block coord
      tile%p_ilo = d_blk%ibf + (itile - 1)*tile_size(1) + 1
      tile%p_ihi = d_blk%ibf + (itile - 1)*tile_size(1) + isize
      tile%p_jlo = d_blk%jbf + (jtile - 1)*tile_size(2) + 1
      tile%p_jhi = d_blk%jbf + (jtile - 1)*tile_size(2) + jsize
      
      tile%partisioned = .false.
      
      tile%parent_tile => blk
      if(.not.associated(blk%tile_head)) then
        blk%tile_head => tile
        blk%tile_tail => tile
        
        tile%prev_tile  => null()
        tile%next_tile  => null()
      else
        blk%tile_tail%next_tile => tile
        tile%prev_tile          => blk%tile_tail
        tile%next_tile          => null()
        blk%tile_tail           => tile
      endif
      
      d_tile => tile%data
      do j = 1, jsize + 2*d_tile%jbf + 1
        do i = 1, isize + 2*d_tile%ibf + 1
          i0 = i - 1 + tile%p_ilo - d_tile%ibf
          j0 = j - 1 + tile%p_jlo - d_tile%jbf
          
          d_tile%xv%p(i,j) = d_blk%xv%p(i0,j0)
          d_tile%yv%p(i,j) = d_blk%yv%p(i0,j0)
        enddo
      enddo
      d_tile%dxi  = d_blk%dxi
      d_tile%deta = d_blk%deta
      call d_tile%cmpt_cell_geo()
      call d_tile%cmpt_face_geo()
      call d_tile%cmpt_metrics()
    enddo
  enddo

  do jtile = 1, blk%tile_dim(2)
    do itile = 1, blk%tile_dim(1)
      tile => blk%tiles(itile,jtile) 
      if((itile - 1).ge.1) then
        tile%west => blk%tiles(itile - 1,jtile)        
      else
        tile%west => null()
      endif

      if(((itile - 1).ge.1).and. &
         ((jtile - 1).ge.1)) then
        tile%southwest => blk%tiles(itile - 1,jtile - 1)        
      else
        tile%southwest => null()
      endif

      if(((itile - 1).ge.1).and. &
         ((jtile + 1).le.blk%tile_dim(2))) then
        tile%northwest => blk%tiles(itile - 1,jtile + 1)        
      else
        tile%northwest => null()
      endif

      if((jtile - 1).ge.1) then 
        tile%south => blk%tiles(itile, jtile - 1)        
      else
        tile%south => null()
      endif
      
      if((jtile + 1).le.blk%tile_dim(2)) then
        tile%north => blk%tiles(itile, jtile + 1)        
      else
        tile%north => null()
      endif
      
      if((itile + 1).le.blk%tile_dim(1)) then
        tile%east => blk%tiles(itile + 1,jtile)        
      else
        tile%east => null()
      endif
      
      if(((itile + 1).le.blk%tile_dim(1)).and. &
         ((jtile - 1).ge.1)) then
        tile%southeast => blk%tiles(itile + 1,jtile - 1)        
      else
        tile%southeast => null()
      endif

      if(((itile + 1).le.blk%tile_dim(1)).and. &
         ((jtile + 1).le.blk%tile_dim(2))) then
        tile%northeast => blk%tiles(itile + 1,jtile + 1)        
      else
        tile%northeast => null()
      endif
    enddo
  enddo

  ! create interior interface condition
  do jtile = 1, blk%tile_dim(2)  
    do itile = 1, blk%tile_dim(1)
      tile   => blk%tiles(itile,jtile) 
      d_tile => tile%data
      if(associated(tile%east)) then
        curr          => tile%add_bndy()
        bndy          => curr%bndy
        d_bndy        => bndy%data
        bndy%group_id = right_bndy_id
        
        tile0   => tile%east
        d_tile0 => tile0%data
        
        d_bndy%ilo(1) = d_tile%ibf + d_tile%isize + 1
        d_bndy%ihi(1) = d_tile%ibf + d_tile%isize + d_tile%ibf
        d_bndy%jlo(1) = d_tile%jbf + 1
        d_bndy%jhi(1) = d_tile%jbf + d_tile%jsize

        d_bndy%ilo(2) = d_tile0%ibf + 1
        d_bndy%ihi(2) = d_tile0%ibf + d_tile0%ibf
        d_bndy%jlo(2) = d_tile0%jbf + 1
        d_bndy%jhi(2) = d_tile0%jbf + d_tile0%jsize

        call create(d_bndy%bc_infos,2)
        d_bndy%bc_infos%p(1) = bc_idx_intf
        d_bndy%bc_infos%p(2) = tile0%idx
      endif

      if(associated(tile%south)) then
        curr          => tile%add_bndy()
        bndy          => curr%bndy
        d_bndy        => bndy%data
        bndy%group_id = bottom_bndy_id

        tile0   => tile%south
        d_tile0 => tile0%data

        d_bndy%ilo(1) = d_tile%ibf + 1
        d_bndy%ihi(1) = d_tile%ibf + d_tile%isize
        d_bndy%jlo(1) = 1
        d_bndy%jhi(1) = d_tile%jbf

        d_bndy%ilo(2) = d_tile0%ibf + 1
        d_bndy%ihi(2) = d_tile0%ibf + d_tile0%isize
        d_bndy%jlo(2) = d_tile0%jsize + 1
        d_bndy%jhi(2) = d_tile0%jsize + d_tile%jbf

        call create(d_bndy%bc_infos,2)
        d_bndy%bc_infos%p(1) = bc_idx_intf
        d_bndy%bc_infos%p(2) = tile0%idx
      endif

      if(associated(tile%west)) then
        curr          => tile%add_bndy()
        bndy          => curr%bndy
        d_bndy        => bndy%data
        bndy%group_id = left_bndy_id

        tile0   => tile%west
        d_tile0 => tile0%data

        d_bndy%ilo(1) = 1
        d_bndy%ihi(1) = d_tile%ibf
        d_bndy%jlo(1) = d_tile%jbf + 1
        d_bndy%jhi(1) = d_tile%jbf + d_tile%jsize

        d_bndy%ilo(2) = d_tile0%isize + 1
        d_bndy%ihi(2) = d_tile0%isize + d_tile0%ibf
        d_bndy%jlo(2) = d_tile0%jbf + 1
        d_bndy%jhi(2) = d_tile0%jbf + d_tile0%jsize

        call create(d_bndy%bc_infos,2)
        d_bndy%bc_infos%p(1) = bc_idx_intf
        d_bndy%bc_infos%p(2) = tile0%idx
      endif

      if(associated(tile%north)) then
        curr          => tile%add_bndy()
        bndy          => curr%bndy
        d_bndy        => bndy%data
        bndy%group_id = top_bndy_id

        tile0   => tile%north
        d_tile0 => tile0%data

        d_bndy%ilo(1) = d_tile%ibf + 1
        d_bndy%ihi(1) = d_tile%ibf + d_tile%isize
        d_bndy%jlo(1) = d_tile%jbf + d_tile%jsize + 1
        d_bndy%jhi(1) = d_tile%jbf + d_tile%jsize + d_tile%jbf

        d_bndy%ilo(2) = d_tile0%ibf + 1
        d_bndy%ihi(2) = d_tile0%ibf + d_tile0%isize
        d_bndy%jlo(2) = d_tile0%jbf + 1
        d_bndy%jhi(2) = d_tile0%jbf + d_tile0%jbf

        call create(d_bndy%bc_infos,2)
        d_bndy%bc_infos%p(1) = bc_idx_intf
        d_bndy%bc_infos%p(2) = tile0%idx
      endif
    enddo
  enddo
  blk%partisioned = .true.
end subroutine split_block_to_tiles

subroutine split_bndy_to_tiles(blk)
  implicit none
  class(block_t),target,intent(inout) :: blk
  class(block_data_t),pointer         :: d_blk 
   
  class(block_list_t),pointer         :: blk_node
  class(block_t),pointer              :: blk0
  class(block_data_t),pointer         :: d_blk0 
  
  class(bndy_list_t),pointer          :: bndy_node  
  class(bndy_t),pointer               :: bndy  
  class(bndy_data_t),pointer          :: d_bndy  

  class(bndy_list_t),pointer          :: bndy_list_split  
  class(bndy_list_t),pointer          :: bndy_node_split 
  class(bndy_t),pointer               :: bndy_split  
  class(bndy_data_t),pointer          :: d_bndy_split  

  class(block_t),pointer              :: tile
  class(block_data_t),pointer         :: d_tile  
  class(bndy_list_t),pointer          :: bndy_node_tile
  class(bndy_t),pointer               :: bndy_tile  
  class(bndy_data_t),pointer          :: d_bndy_tile 

  integer :: donner_idx
  integer :: itile,jtile
  integer :: ilo(2),jlo(2)
  integer :: ihi(2),jhi(2)
  integer :: dlo,dhi
  integer :: ival
  
  d_blk      => blk%data
  bndy_node  => blk%bndy_list
  do while(associated(bndy_node))
    bndy_list_split => null()

    bndy   => bndy_node%bndy
    d_bndy => bndy%data   

    ! Split boundary by adjacent block tiles if any
    if(d_bndy%bc_infos%p(1).eq.bc_idx_intf) then
      blk_node => get_blk(base_blk_list,d_bndy%bc_infos%p(2))
      blk0     => blk_node%blk
      d_blk0   => blk0%data
      select case(bndy%group_id)
        case(left_bndy_id)
          do jtile = 1, blk0%tile_dim(2)  
            tile   => blk0%tiles(blk0%tile_dim(1),jtile)
            d_tile => tile%data
            if((tile%p_jlo.gt.d_bndy%jhi(2)).or. &
               (tile%p_jhi.lt.d_bndy%jlo(2))) then
              cycle
            endif         

            dlo = max(d_bndy%jlo(2) - tile%p_jlo, 0)
            dhi = min(d_bndy%jhi(2) - tile%p_jlo, tile%p_jhi - tile%p_jlo)

            ilo(2) = d_tile%isize + 1
            ihi(2) = d_tile%isize + d_tile%ibf
            jlo(2) = dlo + 1 + d_tile%jbf
            jhi(2) = dhi + 1 + d_tile%jbf

            ilo(1) = d_bndy%ilo(1)
            ihi(1) = d_bndy%ihi(1)
            jlo(1) = d_bndy%jlo(1) - min(d_bndy%jlo(2) - tile%p_jlo, 0)
            jhi(1) = d_bndy%jhi(1) - max(d_bndy%jhi(2) - tile%p_jhi, 0)

            bndy_node_split => push_back_bndy(bndy_list_split)
            bndy_split      => bndy_node_split%bndy
            d_bndy_split    => bndy_split%data

            bndy_split%group_id = left_bndy_id
            d_bndy_split%ilo    = ilo
            d_bndy_split%ihi    = ihi
            d_bndy_split%jlo    = jlo
            d_bndy_split%jhi    = jhi

            call create(d_bndy_split%bc_infos,d_bndy%bc_infos%d)
            d_bndy_split%bc_infos%p(1) = bc_idx_intf
            d_bndy_split%bc_infos%p(2) = tile%idx
          enddo       
        case(bottom_bndy_id)
          do itile = 1, blk0%tile_dim(1)  
            tile   => blk0%tiles(itile,blk0%tile_dim(2))
            d_tile => tile%data
            if((tile%p_ilo.gt.d_bndy%ihi(2)).or. &
               (tile%p_ihi.lt.d_bndy%ilo(2))) then                
              cycle
            endif         
            
            dlo = max(d_bndy%ilo(2) - tile%p_ilo, 0)
            dhi = min(d_bndy%ihi(2) - tile%p_ilo, tile%p_ihi - tile%p_ilo)

            ilo(2) = dlo + 1 + d_tile%ibf
            ihi(2) = dhi + 1 + d_tile%ibf
            jlo(2) = d_tile%jsize + 1
            jhi(2) = d_tile%jsize + d_tile%jbf 

            ilo(1) = d_bndy%ilo(1) - min(d_bndy%ilo(2) - tile%p_ilo, 0)
            ihi(1) = d_bndy%ihi(1) - max(d_bndy%ihi(2) - tile%p_ihi, 0)
            jlo(1) = d_bndy%jlo(1)
            jhi(1) = d_bndy%jhi(1)

            bndy_node_split => push_back_bndy(bndy_list_split)
            bndy_split      => bndy_node_split%bndy
            d_bndy_split    => bndy_split%data

            bndy_split%group_id = bottom_bndy_id
            d_bndy_split%ilo    = ilo
            d_bndy_split%ihi    = ihi
            d_bndy_split%jlo    = jlo
            d_bndy_split%jhi    = jhi

            call create(d_bndy_split%bc_infos,d_bndy%bc_infos%d)
            d_bndy_split%bc_infos%p(1) = bc_idx_intf
            d_bndy_split%bc_infos%p(2) = tile%idx
          enddo 
        case(right_bndy_id)
          do jtile = 1, blk0%tile_dim(2)  
            tile   => blk0%tiles(1,jtile)
            d_tile => tile%data
            if((tile%p_jlo.gt.d_bndy%jhi(2)).or. &
               (tile%p_jhi.lt.d_bndy%jlo(2))) then
              cycle
            endif         
            
            dlo = max(d_bndy%jlo(2) - tile%p_jlo, 0)
            dhi = min(d_bndy%jhi(2) - tile%p_jlo, tile%p_jhi - tile%p_jlo)

            ilo(2) = d_tile%ibf + 1
            ihi(2) = d_tile%ibf + d_tile%ibf
            jlo(2) = d_tile%jbf + 1 + dlo
            jhi(2) = d_tile%jbf + 1 + dhi 

            ilo(1) = d_bndy%ilo(1)
            ihi(1) = d_bndy%ihi(1)
            jlo(1) = d_bndy%jlo(1) - min(d_bndy%jlo(2) - tile%p_jlo, 0)
            jhi(1) = d_bndy%jhi(1) - max(d_bndy%jhi(2) - tile%p_jhi, 0)

            bndy_node_split => push_back_bndy(bndy_list_split)
            bndy_split      => bndy_node_split%bndy
            d_bndy_split    => bndy_split%data

            bndy_split%group_id = right_bndy_id
            d_bndy_split%ilo    = ilo
            d_bndy_split%ihi    = ihi
            d_bndy_split%jlo    = jlo
            d_bndy_split%jhi    = jhi

            call create(d_bndy_split%bc_infos,d_bndy%bc_infos%d)
            d_bndy_split%bc_infos%p(1) = bc_idx_intf
            d_bndy_split%bc_infos%p(2) = tile%idx
          enddo
        case(top_bndy_id)    
          do itile = 1, blk%tile_dim(1)  
            tile   => blk0%tiles(itile,1)
            d_tile => tile%data
            if((tile%p_ilo.gt.d_bndy%ihi(2)).or. &
               (tile%p_ihi.lt.d_bndy%ilo(2))) then
              cycle
            endif         
            
            dlo = max(d_bndy%ilo(2) - tile%p_ilo, 0)
            dhi = min(d_bndy%ihi(2) - tile%p_ilo, tile%p_ihi - tile%p_ilo)

            ilo(2) = dlo + 1 + d_tile%ibf
            ihi(2) = dhi + 1 + d_tile%ibf
            jlo(2) = d_tile%jbf + 1
            jhi(2) = d_tile%jbf + d_tile%jbf 

            ilo(1) = d_bndy%ilo(1) - min(d_bndy%ilo(2) - tile%p_ilo, 0)
            ihi(1) = d_bndy%ihi(1) - max(d_bndy%ihi(2) - tile%p_ihi, 0)
            jlo(1) = d_bndy%jlo(1)
            jhi(1) = d_bndy%jhi(1)

            bndy_node_split => push_back_bndy(bndy_list_split)
            bndy_split      => bndy_node_split%bndy
            d_bndy_split    => bndy_split%data

            bndy_split%group_id = top_bndy_id
            d_bndy_split%ilo    = ilo
            d_bndy_split%ihi    = ihi
            d_bndy_split%jlo    = jlo
            d_bndy_split%jhi    = jhi

            call create(d_bndy_split%bc_infos,d_bndy%bc_infos%d)
            d_bndy_split%bc_infos%p(1) = bc_idx_intf
            d_bndy_split%bc_infos%p(2) = tile%idx
          enddo 
        case default
          print*,"Error: split_bndy_to_tiles()"
          print*,"No such group_id",bndy%group_id, &
                "for bndy",bndy%idx,              &
                "of block",blk%idx
                stop
      end select 
    else
      bndy_node_split => push_back_bndy(bndy_list_split)
      bndy_split      => bndy_node_split%bndy
      d_bndy_split    => bndy_split%data

      bndy_split%group_id = bndy%group_id
      d_bndy_split%ilo    = d_bndy%ilo
      d_bndy_split%ihi    = d_bndy%ihi
      d_bndy_split%jlo    = d_bndy%jlo
      d_bndy_split%jhi    = d_bndy%jhi

      call create(d_bndy_split%bc_infos,d_bndy%bc_infos%d)
      d_bndy_split%bc_infos%p = d_bndy%bc_infos%p

      if(d_bndy%scalar_num.ne.0) then
        d_bndy_split%scalar_num = d_bndy%scalar_num
        allocate(d_bndy_split%scalar(d_bndy_split%scalar_num))
        d_bndy_split%scalar = d_bndy%scalar
      endif

      if(d_bndy%vector_num.ne.0) then
        d_bndy_split%vector_num = d_bndy%vector_num
        allocate(d_bndy_split%vector(d_bndy_split%vector_num))
        do ival = 1, d_bndy_split%vector_num
          call create(d_bndy_split%vector(ival),d_bndy%vector(ival)%d)
          d_bndy_split%vector(ival)%p = d_bndy%vector(ival)%p    
        enddo
      endif

      if(d_bndy%matrix_num.ne.0) then
        d_bndy_split%matrix_num = d_bndy%matrix_num
        allocate(d_bndy_split%matrix(d_bndy_split%matrix_num))
        do ival = 1, d_bndy_split%matrix_num
          call create(d_bndy_split%matrix(ival), &
                      d_bndy%matrix(ival)%d(1), &
                      d_bndy%matrix(ival)%d(2))
          d_bndy_split%matrix(ival)%p = d_bndy%matrix(ival)%p    
        enddo
      endif

      if(d_bndy%tensor_num.ne.0) then
        d_bndy_split%tensor_num = d_bndy%tensor_num
        allocate(d_bndy_split%tensor(d_bndy_split%tensor_num))
        do ival = 1, d_bndy_split%tensor_num
          call create(d_bndy_split%tensor(ival), &
                      d_bndy%tensor(ival)%d(1), &
                      d_bndy%tensor(ival)%d(2), &
                      d_bndy%tensor(ival)%d(3))
          d_bndy_split%tensor(ival)%p = d_bndy%tensor(ival)%p    
        enddo
      endif
    endif
        
    ! Set local tile boundaries
    bndy_node_split => bndy_list_split
    do while(associated(bndy_node_split))
      bndy   => bndy_node_split%bndy
      d_bndy => bndy%data  
      select case(bndy%group_id)
        case(left_bndy_id)
          do jtile = 1, blk%tile_dim(2)  
            tile   => blk%tiles(1,jtile)
            d_tile => tile%data
            if((tile%p_jlo.gt.d_bndy%jhi(1)).or. &
              (tile%p_jhi.lt.d_bndy%jlo(1))) then
              cycle
            endif
            dlo = max(d_bndy%jlo(1) - tile%p_jlo, 0)
            dhi = min(d_bndy%jhi(1) - tile%p_jlo, tile%p_jhi - tile%p_jlo)

            ilo(1) = 1
            ihi(1) = d_tile%ibf
            jlo(1) = d_tile%jbf + 1 + dlo
            jhi(1) = d_tile%jbf + 1 + dhi 

            ilo(2) = d_bndy%ilo(2)
            ihi(2) = d_bndy%ihi(2)
            jlo(2) = d_bndy%jlo(2) - min(d_bndy%jlo(1) - tile%p_jlo, 0)
            jhi(2) = d_bndy%jhi(2) - max(d_bndy%jhi(1) - tile%p_jhi, 0)

            bndy_node_tile     => tile%add_bndy()
            bndy_tile          => bndy_node_tile%bndy
            d_bndy_tile        => bndy_tile%data
            bndy_tile%group_id = left_bndy_id

            d_bndy_tile%ilo    = ilo
            d_bndy_tile%ihi    = ihi
            d_bndy_tile%jlo    = jlo
            d_bndy_tile%jhi    = jhi
            call create(d_bndy_tile%bc_infos,d_bndy%bc_infos%d)
            d_bndy_tile%bc_infos%p = d_bndy%bc_infos%p
            
            if(d_bndy%scalar_num.ne.0) then
              d_bndy_tile%scalar_num = d_bndy%scalar_num
              allocate(d_bndy_tile%scalar(d_bndy_tile%scalar_num))
              d_bndy_tile%scalar = d_bndy%scalar
            endif

            if(d_bndy%vector_num.ne.0) then
              d_bndy_tile%vector_num = d_bndy%vector_num
              allocate(d_bndy_tile%vector(d_bndy_tile%vector_num))
              do ival = 1, d_bndy_tile%vector_num
                call create(d_bndy_tile%vector(ival),d_bndy%vector(ival)%d)
                d_bndy_tile%vector(ival)%p = d_bndy%vector(ival)%p    
              enddo
            endif

            if(d_bndy%matrix_num.ne.0) then
              d_bndy_tile%matrix_num = d_bndy%matrix_num
              allocate(d_bndy_tile%matrix(d_bndy_tile%matrix_num))
              do ival = 1, d_bndy_tile%matrix_num
                call create(d_bndy_tile%matrix(ival), &
                            d_bndy%matrix(ival)%d(1), &
                            d_bndy%matrix(ival)%d(2))
                d_bndy_tile%matrix(ival)%p = d_bndy%matrix(ival)%p    
              enddo
            endif

            if(d_bndy%tensor_num.ne.0) then
              d_bndy_tile%tensor_num = d_bndy%tensor_num
              allocate(d_bndy_tile%tensor(d_bndy_tile%tensor_num))
              do ival = 1, d_bndy_tile%tensor_num
                call create(d_bndy_tile%tensor(ival), &
                            d_bndy%tensor(ival)%d(1), &
                            d_bndy%tensor(ival)%d(2), &
                            d_bndy%tensor(ival)%d(3))
                d_bndy_tile%tensor(ival)%p = d_bndy%tensor(ival)%p    
              enddo
            endif
          enddo       
        case(bottom_bndy_id)
          do itile = 1, blk%tile_dim(1)  
            tile   => blk%tiles(itile,1)
            d_tile => tile%data
            if((tile%p_ilo.gt.d_bndy%ihi(1)).or. &
              (tile%p_ihi.lt.d_bndy%ilo(1))) then
              cycle
            endif
            dlo = max(d_bndy%ilo(1) - tile%p_ilo, 0)
            dhi = min(d_bndy%ihi(1) - tile%p_ilo, tile%p_ihi - tile%p_ilo)


            ilo(1) = d_tile%ibf + 1 + dlo
            ihi(1) = d_tile%ibf + 1 + dhi
            jlo(1) = 1
            jhi(1) = d_tile%jbf 

            ilo(2) = d_bndy%ilo(2) - min(d_bndy%ilo(1) - tile%p_ilo, 0)
            ihi(2) = d_bndy%ihi(2) - max(d_bndy%ihi(1) - tile%p_ihi, 0)
            jlo(2) = d_bndy%jlo(2)
            jhi(2) = d_bndy%jhi(2)

            bndy_node_tile     => tile%add_bndy()
            bndy_tile          => bndy_node_tile%bndy
            d_bndy_tile        => bndy_tile%data
            bndy_tile%group_id = bottom_bndy_id

            d_bndy_tile%ilo    = ilo
            d_bndy_tile%ihi    = ihi
            d_bndy_tile%jlo    = jlo
            d_bndy_tile%jhi    = jhi
            call create(d_bndy_tile%bc_infos,d_bndy%bc_infos%d)
            d_bndy_tile%bc_infos%p = d_bndy%bc_infos%p

            if(d_bndy%scalar_num.ne.0) then
              d_bndy_tile%scalar_num = d_bndy%scalar_num
              allocate(d_bndy_tile%scalar(d_bndy_tile%scalar_num))
              d_bndy_tile%scalar = d_bndy%scalar
            endif

            if(d_bndy%vector_num.ne.0) then
              d_bndy_tile%vector_num = d_bndy%vector_num
              allocate(d_bndy_tile%vector(d_bndy_tile%vector_num))
              do ival = 1, d_bndy_tile%vector_num
                call create(d_bndy_tile%vector(ival),d_bndy%vector(ival)%d)
                d_bndy_tile%vector(ival)%p = d_bndy%vector(ival)%p    
              enddo
            endif

            if(d_bndy%matrix_num.ne.0) then
              d_bndy_tile%matrix_num = d_bndy%matrix_num
              allocate(d_bndy_tile%matrix(d_bndy_tile%matrix_num))
              do ival = 1, d_bndy_tile%matrix_num
                call create(d_bndy_tile%matrix(ival), &
                            d_bndy%matrix(ival)%d(1), &
                            d_bndy%matrix(ival)%d(2))
                d_bndy_tile%matrix(ival)%p = d_bndy%matrix(ival)%p    
              enddo
            endif

            if(d_bndy%tensor_num.ne.0) then
              d_bndy_tile%tensor_num = d_bndy%tensor_num
              allocate(d_bndy_tile%tensor(d_bndy_tile%tensor_num))
              do ival = 1, d_bndy_tile%tensor_num
                call create(d_bndy_tile%tensor(ival), &
                            d_bndy%tensor(ival)%d(1), &
                            d_bndy%tensor(ival)%d(2), &
                            d_bndy%tensor(ival)%d(3))
                d_bndy_tile%tensor(ival)%p = d_bndy%tensor(ival)%p    
              enddo
            endif
          enddo 
        case(right_bndy_id)
          do jtile = 1, blk%tile_dim(2)  
            tile   => blk%tiles(blk%tile_dim(1),jtile)
            d_tile => tile%data
            if((tile%p_jlo.gt.d_bndy%jhi(1)).or. &
              (tile%p_jhi.lt.d_bndy%jlo(1))) then
              cycle
            endif
            dlo = max(d_bndy%jlo(1) - tile%p_jlo, 0)
            dhi = min(d_bndy%jhi(1) - tile%p_jlo, tile%p_jhi - tile%p_jlo)

            ilo(1) = d_tile%ibf + d_tile%isize + 1
            ihi(1) = d_tile%ibf + d_tile%isize + d_tile%ibf
            jlo(1) = d_tile%jbf + 1 + dlo
            jhi(1) = d_tile%jbf + 1 + dhi 

            ilo(2) = d_bndy%ilo(2)
            ihi(2) = d_bndy%ihi(2)
            jlo(2) = d_bndy%jlo(2) - min(d_bndy%jlo(1) - tile%p_jlo, 0)
            jhi(2) = d_bndy%jhi(2) - max(d_bndy%jhi(1) - tile%p_jhi, 0)

            bndy_node_tile     => tile%add_bndy()
            bndy_tile          => bndy_node_tile%bndy
            d_bndy_tile        => bndy_tile%data
            bndy_tile%group_id = right_bndy_id

            d_bndy_tile%ilo    = ilo
            d_bndy_tile%ihi    = ihi
            d_bndy_tile%jlo    = jlo
            d_bndy_tile%jhi    = jhi

            call create(d_bndy_tile%bc_infos,d_bndy%bc_infos%d)
            d_bndy_tile%bc_infos%p = d_bndy%bc_infos%p
            
            if(d_bndy%scalar_num.ne.0) then
              d_bndy_tile%scalar_num = d_bndy%scalar_num
              allocate(d_bndy_tile%scalar(d_bndy_tile%scalar_num))
              d_bndy_tile%scalar = d_bndy%scalar
            endif

            if(d_bndy%vector_num.ne.0) then
              d_bndy_tile%vector_num = d_bndy%vector_num
              allocate(d_bndy_tile%vector(d_bndy_tile%vector_num))
              do ival = 1, d_bndy_tile%vector_num
                call create(d_bndy_tile%vector(ival),d_bndy%vector(ival)%d)
                d_bndy_tile%vector(ival)%p = d_bndy%vector(ival)%p    
              enddo
            endif

            if(d_bndy%matrix_num.ne.0) then
              d_bndy_tile%matrix_num = d_bndy%matrix_num
              allocate(d_bndy_tile%matrix(d_bndy_tile%matrix_num))
              do ival = 1, d_bndy_tile%matrix_num
                call create(d_bndy_tile%matrix(ival), &
                            d_bndy%matrix(ival)%d(1), &
                            d_bndy%matrix(ival)%d(2))
                d_bndy_tile%matrix(ival)%p = d_bndy%matrix(ival)%p    
              enddo
            endif

            if(d_bndy%tensor_num.ne.0) then
              d_bndy_tile%tensor_num = d_bndy%tensor_num
              allocate(d_bndy_tile%tensor(d_bndy_tile%tensor_num))
              do ival = 1, d_bndy_tile%tensor_num
                call create(d_bndy_tile%tensor(ival), &
                            d_bndy%tensor(ival)%d(1), &
                            d_bndy%tensor(ival)%d(2), &
                            d_bndy%tensor(ival)%d(3))
                d_bndy_tile%tensor(ival)%p = d_bndy%tensor(ival)%p    
              enddo
            endif
          enddo
        case(top_bndy_id)    
          do itile = 1, blk%tile_dim(1)  
            tile   => blk%tiles(itile,blk%tile_dim(2))
            d_tile => tile%data
            if((tile%p_ilo.gt.d_bndy%ihi(1)).or. &
              (tile%p_ihi.lt.d_bndy%ilo(1))) then
              cycle
            endif
            dlo = max(d_bndy%ilo(1) - tile%p_ilo, 0)
            dhi = min(d_bndy%ihi(1) - tile%p_ilo, tile%p_ihi - tile%p_ilo)

            ilo(1) = d_tile%ibf + 1 + dlo
            ihi(1) = d_tile%ibf + 1 + dhi
            jlo(1) = d_tile%jbf + d_tile%jsize + 1
            jhi(1) = d_tile%jbf + d_tile%jsize + d_tile%jbf 

            ilo(2) = d_bndy%ilo(2) - min(d_bndy%ilo(1) - tile%p_ilo, 0)
            ihi(2) = d_bndy%ihi(2) - max(d_bndy%ihi(1) - tile%p_ihi, 0)
            jlo(2) = d_bndy%jlo(2)
            jhi(2) = d_bndy%jhi(2)

            bndy_node_tile     => tile%add_bndy()
            bndy_tile          => bndy_node_tile%bndy
            d_bndy_tile        => bndy_tile%data
            bndy_tile%group_id = top_bndy_id

            d_bndy_tile%ilo    = ilo
            d_bndy_tile%ihi    = ihi
            d_bndy_tile%jlo    = jlo
            d_bndy_tile%jhi    = jhi
            call create(d_bndy_tile%bc_infos,d_bndy%bc_infos%d)
            d_bndy_tile%bc_infos%p = d_bndy%bc_infos%p

            if(d_bndy%scalar_num.ne.0) then
              d_bndy_tile%scalar_num = d_bndy%scalar_num
              allocate(d_bndy_tile%scalar(d_bndy_tile%scalar_num))
              d_bndy_tile%scalar = d_bndy%scalar
            endif

            if(d_bndy%vector_num.ne.0) then
              d_bndy_tile%vector_num = d_bndy%vector_num
              allocate(d_bndy_tile%vector(d_bndy_tile%vector_num))
              do ival = 1, d_bndy_tile%vector_num
                call create(d_bndy_tile%vector(ival),d_bndy%vector(ival)%d)
                d_bndy_tile%vector(ival)%p = d_bndy%vector(ival)%p    
              enddo
            endif

            if(d_bndy%matrix_num.ne.0) then
              d_bndy_tile%matrix_num = d_bndy%matrix_num
              allocate(d_bndy_tile%matrix(d_bndy_tile%matrix_num))
              do ival = 1, d_bndy_tile%matrix_num
                call create(d_bndy_tile%matrix(ival), &
                            d_bndy%matrix(ival)%d(1), &
                            d_bndy%matrix(ival)%d(2))
                d_bndy_tile%matrix(ival)%p = d_bndy%matrix(ival)%p    
              enddo
            endif

            if(d_bndy%tensor_num.ne.0) then
              d_bndy_tile%tensor_num = d_bndy%tensor_num
              allocate(d_bndy_tile%tensor(d_bndy_tile%tensor_num))
              do ival = 1, d_bndy_tile%tensor_num
                call create(d_bndy_tile%tensor(ival), &
                            d_bndy%tensor(ival)%d(1), &
                            d_bndy%tensor(ival)%d(2), &
                            d_bndy%tensor(ival)%d(3))
                d_bndy_tile%tensor(ival)%p = d_bndy%tensor(ival)%p    
              enddo
            endif
          enddo 
        case default
          print*,"Error: split_bndy_to_tiles()"
          print*,"No such group_id",bndy%group_id, &
                "for bndy",bndy%idx,              &
                "of block",blk%idx
                stop
      end select 
      bndy_node_split => bndy_node_split%next
    enddo
    deallocate(bndy_list_split)
    bndy_node => bndy_node%next
  enddo
end subroutine split_bndy_to_tiles

subroutine reset_block_tiles(blk)
  implicit none
  class(block_t),target,intent(inout) :: blk
  class(block_t),pointer :: tile
  integer :: itile,jtile
  if(blk%partisioned) then
    do jtile = 1, blk%tile_dim(2)
      do itile = 1, blk%tile_dim(1)
        tile => blk%tiles(itile,jtile)
        call tile%destroy_block()
      enddo
    enddo
    deallocate(blk%tiles)
  endif
end subroutine reset_block_tiles

function add_bndy(blk,group_id,ilo,ihi,jlo,jhi) result(bndy_node)
  implicit none
  class(block_t),target,intent(inout) :: blk
  integer,optional,intent(in) :: group_id
  integer,dimension(2),optional,intent(in) :: ilo,ihi,jlo,jhi
  class(bndy_list_t), pointer :: bndy_node
  class(bndy_list_t), pointer :: tail
  class(bndy_t), pointer :: bndy

  if(associated(blk%bndy_list)) then
    tail => blk%bndy_list_tail
    allocate(bndy_node)
    allocate(bndy_node%bndy)
    allocate(bndy_node%bndy%data)
    blk%bndy_num   = blk%bndy_num + 1
    blk%bndy_count = blk%bndy_count + 1
    bndy => bndy_node%bndy
    bndy%idx       = blk%bndy_count 

    bndy_node%prev  => tail
    tail%next => bndy_node
    bndy_node%next => null()
    blk%bndy_list_tail => bndy_node
  else
    allocate(bndy_node)
    allocate(bndy_node%bndy)
    allocate(bndy_node%bndy%data)
    blk%bndy_num   = 1
    blk%bndy_count = 1
    bndy => bndy_node%bndy
    bndy%idx       = blk%bndy_count
    
    blk%bndy_list      => bndy_node
    blk%bndy_list_tail => bndy_node
    bndy_node%prev  => null()
    bndy_node%next => null()  
  endif

  if(present(group_id)) then
    bndy%group_id  = group_id
  else
    bndy%group_id  = 0
  endif

  if(present(ilo)) then
    bndy%data%ilo  = ilo
  else
    bndy%data%ilo  = 0
  endif

  if(present(ihi)) then
    bndy%data%ihi  = ihi
  else
    bndy%data%ihi  = 0
  endif

  if(present(jlo)) then
    bndy%data%jlo  = jlo
  else
    bndy%data%jlo  = 0
  endif

  if(present(jhi)) then
    bndy%data%jhi  = jhi
  else
    bndy%data%jhi  = 0
  endif
end function add_bndy

function get_bndy(blk,bndy_idx) result(bndy_node)
  implicit none
  class(block_t),target,intent(inout) :: blk
  integer,intent(in) :: bndy_idx
  type(bndy_list_t), pointer :: bndy_node
  
  bndy_node => blk%bndy_list
  do while(associated(bndy_node))
    if(bndy_node%bndy%idx.eq.bndy_idx) return
    bndy_node => bndy_node%next
  enddo
end function get_bndy

subroutine rmv_bndy(blk,bndy_idx)
  implicit none
  class(block_t),target,intent(inout) :: blk
  integer,intent(in) :: bndy_idx
  type(bndy_list_t), pointer :: bndy_node
  type(bndy_list_t), pointer :: prev,next

  bndy_node => blk%bndy_list
  do while(associated(bndy_node))
    if(bndy_node%bndy%idx.eq.bndy_idx) then
      prev => bndy_node%prev
      next => bndy_node%next

      if(associated(prev)) then
        prev%next => next
      endif

      if(associated(next)) then
        next%prev => prev
      endif
      bndy_node%prev => null()
      bndy_node%next => null()
      deallocate(bndy_node)
      return
    endif
    bndy_node => bndy_node%next
  enddo
end subroutine rmv_bndy

subroutine allocate_d_blk(d_blk,isize,jsize,ibf,jbf)
  implicit none
  integer,intent(in) :: isize,jsize,ibf,jbf
  class(block_data_t),intent(inout) :: d_blk
  type(fscalarshape_t) :: vertex_s,cell_s
  type(fmatrixshape_t) :: face_s 
  type(fmatrixshape_t) :: jac_s

  d_blk%ibf   = ibf
  d_blk%jbf   = jbf
  d_blk%isize = isize
  d_blk%jsize = jsize

  vertex_s%d(1) = isize + 2*ibf + 1
  vertex_s%d(2) = jsize + 2*jbf + 1
  call create(d_blk%xv,vertex_s)
  call create(d_blk%yv,vertex_s)

  cell_s%d(1) = isize + 2*ibf
  cell_s%d(2) = jsize + 2*jbf
  call create(d_blk%xc,cell_s)
  call create(d_blk%yc,cell_s)

  ! face_s%d(1) = 2
  ! face_s%d(2) = 4
  ! face_s%d(3) = isize + 2*ibf
  ! face_s%d(4) = jsize + 2*jbf
  ! call create(d_blk%f_cntr, face_s)
  ! call create(d_blk%f_norm, face_s)

  ! if(mesh_type.eq.mesh_type_curv) then
  !   jac_s%d(1) = 2
  !   jac_s%d(2) = 2
  !   jac_s%d(3) = isize + 2*ibf
  !   jac_s%d(4) = jsize + 2*jbf
  !   call create(d_blk%x_xi,cell_s)
  !   call create(d_blk%y_xi,cell_s)
  !   call create(d_blk%x_eta,cell_s)
  !   call create(d_blk%y_eta,cell_s)
  !   call create(d_blk%det_jac,cell_s)
  !   call create(d_blk%jac,jac_s)
  ! endif
end subroutine allocate_d_blk

subroutine deallocate_d_blk(d_blk)
  implicit none
  class(block_data_t),intent(inout) :: d_blk

  d_blk%ibf   = 0
  d_blk%jbf   = 0
  d_blk%isize = 0
  d_blk%jsize = 0

  call destroy(d_blk%xv)
  call destroy(d_blk%yv)
  call destroy(d_blk%xc)
  call destroy(d_blk%yc)
end subroutine deallocate_d_blk

subroutine initial_d_blk(d_blk,x0,y0,xl,yl)
  implicit none
  class(block_data_t),intent(inout) :: d_blk
  double precision,intent(in) :: x0,y0,xl,yl
  integer :: i,j
  d_blk%dxi  = xl/d_blk%isize
  d_blk%deta = yl/d_blk%jsize
  do j=1,d_blk%jsize + 2*d_blk%jbf + 1
    do i=1,d_blk%isize + 2*d_blk%ibf + 1
      d_blk%xv%p(i,j) = x0 + (i-d_blk%ibf-1.d0)*d_blk%dxi
      d_blk%yv%p(i,j) = y0 + (j-d_blk%jbf-1.d0)*d_blk%deta
    enddo
  enddo
end subroutine initial_d_blk

subroutine cmpt_cell_geo(d_blk)
  implicit none
  class(block_data_t),intent(inout) :: d_blk 
  integer :: i,j 

  do j=1,d_blk%jsize+d_blk%jbf*2
    do i=1,d_blk%isize+d_blk%ibf*2
      d_blk%xc%p(i,j) =(d_blk%xv%p(i,j)   + d_blk%xv%p(i+1,j) &
                       +d_blk%xv%p(i,j+1) + d_blk%xv%p(i+1,j+1))*0.25d0
      d_blk%yc%p(i,j) =(d_blk%yv%p(i,j)   + d_blk%yv%p(i+1,j) &                       
                       +d_blk%yv%p(i,j+1) + d_blk%yv%p(i+1,j+1))*0.25d0
    enddo
  enddo
end subroutine cmpt_cell_geo

subroutine cmpt_face_geo(d_blk)
  implicit none
  class(block_data_t),intent(inout) :: d_blk
  integer :: i,j

end subroutine cmpt_face_geo

subroutine cmpt_metrics(d_blk)
  implicit none
  class(block_data_t),intent(inout) :: d_blk
  integer :: i,j
  
end subroutine cmpt_metrics
  
function push_back_bndy(bndy_list) result(bndy_node)
  implicit none
  class(bndy_list_t),pointer,intent(inout) :: bndy_list 
  class(bndy_list_t),pointer :: bndy_node,curr_node 

  allocate(bndy_node)
  allocate(bndy_node%bndy)
  allocate(bndy_node%bndy%data)
  bndy_node%prev  => null()
  bndy_node%next => null()
  if(.not.associated(bndy_list)) then 
    bndy_list => bndy_node
  else
    curr_node => bndy_list
    do while(associated(curr_node%next))
      curr_node => curr_node%next
    enddo
    curr_node%next => bndy_node
    bndy_node%prev  => curr_node
  endif
end function push_back_bndy

recursive subroutine empty_bndy_list(bndy_list)
  implicit none
  type(bndy_list_t), intent(inout) :: bndy_list
  if(associated(bndy_list%next)) then
    deallocate(bndy_list%next)
  endif
  call bndy_list%bndy%deallocate_d_bndy()
  deallocate(bndy_list%bndy)
  bndy_list%prev => null()
  bndy_list%next => null()
end subroutine empty_bndy_list

subroutine deallocate_d_bndy(bndy)
  implicit none
  class(bndy_t),intent(inout) :: bndy
  call destroy(bndy%data%bc_infos)
  deallocate(bndy%data)
end subroutine deallocate_d_bndy

end module mesh
