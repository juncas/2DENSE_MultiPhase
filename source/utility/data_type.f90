module data_type
  implicit none
  !=========================================================
  !               
  !                        Array
  !
  !---------------------------------------------------------
  ! Bool type arrays
  !---------------------------------------------------------
  type :: bvector_t
    logical,dimension(:),allocatable :: p
    integer :: d
  end type bvector_t

  type :: bmatrix_t
    logical,dimension(:,:),allocatable :: p
    integer :: d(2)
  end type bmatrix_t
  
  type :: btensor_t
    logical,dimension(:,:,:),allocatable :: p
    integer :: d(3)
  end type btensor_t
  
  !---------------------------------------------------------
  ! Integer type arrays
  !---------------------------------------------------------
  type :: ivector_t
    integer,dimension(:),allocatable :: p
    integer :: d
  end type ivector_t

  type :: imatrix_t
    integer,dimension(:,:),allocatable :: p
    integer :: d(2)
  end type imatrix_t
  
  type :: itensor_t
    integer,dimension(:,:,:),allocatable :: p
    integer :: d(3)
  end type itensor_t
  
  !---------------------------------------------------------
  ! Real type arrays
  !---------------------------------------------------------
  type :: rvector_t
    double precision,dimension(:),allocatable :: p
    integer :: d
  end type

  type :: rmatrix_t
    double precision,dimension(:,:),allocatable :: p
    integer,dimension(2) :: d
  end type

  type :: rtensor_t
    double precision,dimension(:,:,:),allocatable :: p
    integer,dimension(3) :: d
  end type
  !---------------------------------------------------------
  !               
  !                        End of Array
  !
  !=========================================================

  !=========================================================
  !               
  !                       Field shape
  !
  !---------------------------------------------------------
  type :: fscalarshape_t
    integer,dimension(2) :: d
  end type

  type :: fvectorshape_t
    integer,dimension(3) :: d
  end type

  type :: fmatrixshape_t
    integer,dimension(4) :: d
  end type

  type :: ftensorshape_t
    integer,dimension(5) :: d
  end type
  !---------------------------------------------------------
  !               
  !                   End of Field shape
  !
  !=========================================================

  !=========================================================
  !               
  !                        Field
  !
  !---------------------------------------------------------
  ! Bool type field
  !---------------------------------------------------------
  type :: bscalarfield_t
    logical,dimension(:,:),allocatable :: p
    type(fscalarshape_t) :: s
  end type 

  type :: bvectorfield_t
    logical,dimension(:,:,:),allocatable :: p
    type(fvectorshape_t) :: s
  end type 

  type :: bmatrixfield_t
    logical,dimension(:,:,:,:),allocatable :: p
    type(fmatrixshape_t) :: s
  end type 

  type :: btensorfield_t
    logical,dimension(:,:,:,:,:),allocatable :: p
    type(ftensorshape_t) :: s
  end type 

  !---------------------------------------------------------
  ! Integeer type field
  !---------------------------------------------------------  
  type :: iscalarfield_t
    integer,dimension(:,:),allocatable :: p
    type(fscalarshape_t) :: s
  end type 

  type :: ivectorfield_t
    integer,dimension(:,:,:),allocatable :: p
    type(fvectorshape_t) :: s
  end type 

  type :: imatrixfield_t
    integer,dimension(:,:,:,:),allocatable :: p
    type(fmatrixshape_t) :: s
  end type 

  type :: itensorfield_t
    integer,dimension(:,:,:,:,:),allocatable :: p
    type(ftensorshape_t) :: s
  end type 

  !---------------------------------------------------------
  ! Real type field
  !---------------------------------------------------------
  type :: rscalarfield_t
    double precision,dimension(:,:),allocatable :: p
    type(fscalarshape_t) :: s
  end type 

  type :: rvectorfield_t
    double precision,dimension(:,:,:),allocatable :: p
    type(fvectorshape_t) :: s
  end type 

  type :: rmatrixfield_t
    double precision,dimension(:,:,:,:),allocatable :: p
    type(fmatrixshape_t) :: s
  end type 

  type :: rtensorfield_t
    double precision,dimension(:,:,:,:,:),allocatable :: p
    type(ftensorshape_t) :: s
  end type 
  !---------------------------------------------------------
  !               
  !                       End of Field
  !
  !=========================================================


  !=========================================================
  !               
  !                       Field shape
  !
  !---------------------------------------------------------
  type :: ascalarshape_t
    integer :: d
  end type

  type :: avectorshape_t
    integer,dimension(2) :: d
  end type

  type :: amatrixshape_t
    integer,dimension(3) :: d
  end type

  type :: atensorshape_t
    integer,dimension(4) :: d
  end type
  !---------------------------------------------------------
  !               
  !                   End of Field shape
  !
  !=========================================================

  !=========================================================
  !               
  !                        array
  !
  !---------------------------------------------------------
  ! Bool type array
  !---------------------------------------------------------
  type :: bscalararray_t
    logical,dimension(:),allocatable :: p
    type(ascalarshape_t) :: s
  end type 

  type :: bvectorarray_t
    logical,dimension(:,:),allocatable :: p
    type(avectorshape_t) :: s
  end type 

  type :: bmatrixarray_t
    logical,dimension(:,:,:),allocatable :: p
    type(amatrixshape_t) :: s
  end type 

  type :: btensorarray_t
    logical,dimension(:,:,:,:),allocatable :: p
    type(atensorshape_t) :: s
  end type 

  !---------------------------------------------------------
  ! Integeer type array
  !---------------------------------------------------------  
  type :: iscalararray_t
    integer,dimension(:),allocatable :: p
    type(ascalarshape_t) :: s
  end type 

  type :: ivectorarray_t
    integer,dimension(:,:),allocatable :: p
    type(avectorshape_t) :: s
  end type 

  type :: imatrixarray_t
    integer,dimension(:,:,:),allocatable :: p
    type(amatrixshape_t) :: s
  end type 

  type :: itensorarray_t
    integer,dimension(:,:,:,:),allocatable :: p
    type(atensorshape_t) :: s
  end type 

  !---------------------------------------------------------
  ! Real type array
  !---------------------------------------------------------
  type :: rscalararray_t
    double precision,dimension(:),allocatable :: p
    type(ascalarshape_t) :: s
  end type 

  type :: rvectorarray_t
    double precision,dimension(:,:),allocatable :: p
    type(avectorshape_t) :: s
  end type 

  type :: rmatrixarray_t
    double precision,dimension(:,:,:),allocatable :: p
    type(amatrixshape_t) :: s
  end type 

  type :: rtensorarray_t
    double precision,dimension(:,:,:,:),allocatable :: p
    type(atensorshape_t) :: s
  end type 
  !---------------------------------------------------------
  !               
  !                       End of Field
  !
  !=========================================================

  interface create
    module procedure create_bvector
    module procedure create_bmatrix
    module procedure create_btensor
    module procedure create_ivector
    module procedure create_imatrix
    module procedure create_itensor
    module procedure create_rvector
    module procedure create_rmatrix
    module procedure create_rtensor
    module procedure create_bscalarfield
    module procedure create_bvectorfield
    module procedure create_bmatrixfield
    module procedure create_btensorfield
    module procedure create_iscalarfield
    module procedure create_ivectorfield
    module procedure create_imatrixfield
    module procedure create_itensorfield
    module procedure create_rscalarfield
    module procedure create_rvectorfield
    module procedure create_rmatrixfield
    module procedure create_rtensorfield
  end interface create

  interface destroy
    module procedure destroy_bvector
    module procedure destroy_bmatrix
    module procedure destroy_btensor
    module procedure destroy_ivector
    module procedure destroy_imatrix
    module procedure destroy_itensor
    module procedure destroy_rvector
    module procedure destroy_rmatrix
    module procedure destroy_rtensor
    module procedure destroy_bscalarfield
    module procedure destroy_bvectorfield
    module procedure destroy_bmatrixfield
    module procedure destroy_btensorfield
    module procedure destroy_iscalarfield
    module procedure destroy_ivectorfield
    module procedure destroy_imatrixfield
    module procedure destroy_itensorfield
    module procedure destroy_rscalarfield
    module procedure destroy_rvectorfield
    module procedure destroy_rmatrixfield
    module procedure destroy_rtensorfield
  end interface destroy

  interface assignment (=)
    module procedure fscalarshape_to_fscalarshape
    module procedure fvectorshape_to_fvectorshape
    module procedure fmatrixshape_to_fmatrixshape
    module procedure ftensorshape_to_ftensorshape
    module procedure bvector_to_bvector
    module procedure bmatrix_to_bmatrix
    module procedure btensor_to_btensor
    module procedure ivector_to_ivector
    module procedure imatrix_to_imatrix
    module procedure itensor_to_itensor
    module procedure rvector_to_rvector
    module procedure rmatrix_to_rmatrix
    module procedure rtensor_to_rtensor
    module procedure bscalarfield_to_bscalar
    module procedure bvectorfield_to_bvector
    module procedure bmatrixfield_to_bmatrix
    module procedure btensorfield_to_btensor
    module procedure iscalarfield_to_iscalar
    module procedure ivectorfield_to_ivector
    module procedure imatrixfield_to_imatrix
    module procedure itensorfield_to_itensor
    module procedure rscalarfield_to_rscalar
    module procedure rvectorfield_to_rvector
    module procedure rmatrixfield_to_rmatrix
    module procedure rtensorfield_to_rtensor
    module procedure bscalarfield_to_bscalarfield
    module procedure bvectorfield_to_bvectorfield
    module procedure bmatrixfield_to_bmatrixfield
    module procedure btensorfield_to_btensorfield
    module procedure iscalarfield_to_iscalarfield
    module procedure ivectorfield_to_ivectorfield
    module procedure imatrixfield_to_imatrixfield
    module procedure itensorfield_to_itensorfield
    module procedure rscalarfield_to_rscalarfield
    module procedure rvectorfield_to_rvectorfield
    module procedure rmatrixfield_to_rmatrixfield
    module procedure rtensorfield_to_rtensorfield
  end interface 

  interface operator (.eq.)
    module procedure fscalarshape_eq_fscalarshape
    module procedure fvectorshape_eq_fvectorshape
    module procedure fmatrixshape_eq_fmatrixshape
    module procedure ftensorshape_eq_ftensorshape
  end interface 

  interface set_shape
    module procedure set_fscalarshape
    module procedure set_fvectorshape
    module procedure set_fmatrixshape
    module procedure set_ftensorshape
  end interface 
  
  integer,parameter :: bias_upwind   = 0
  integer,parameter :: bias_downwind = 1
  integer,parameter :: axis_xi  = 0
  integer,parameter :: axis_eta = 1

contains
subroutine create_bvector(bvector,dimn)
  implicit none
  integer,intent(in) :: dimn
  type(bvector_t),intent(inout) :: bvector    
  if(dimn.ne.0) then
    bvector%d = dimn
  else
    write(*,*) "create_bvector() error:"
    write(*,*) "Zero allocation size:",dimn
    stop
  endif

  if(.not.allocated(bvector%p)) then 
    allocate(bvector%p(bvector%d))
  else
    write(*,*) "create_bvector() error:"
    write(*,*) "rbvector_t has been allocated"
    stop
  endif
end subroutine create_bvector

subroutine create_bmatrix(bmatrix,dimm,dimn)
  implicit none
  integer,intent(in) :: dimm,dimn
  type(bmatrix_t),intent(inout) :: bmatrix    
  if((dimm.ne.0).and. &
     (dimn.ne.0)) then
      bmatrix%d(1) = dimm
      bmatrix%d(2) = dimn
  else
    write(*,*) "create_bmatrix() error:"
    write(*,*) "Zero allocation size:",dimm,dimn
    stop
  endif

  if(.not.allocated(bmatrix%p)) then 
    allocate(bmatrix%p(bmatrix%d(1),bmatrix%d(2)))
  else
    write(*,*) "create_bmatrix() error:"
    write(*,*) "bmatrix_t has been allocated"
    stop
  endif
end subroutine create_bmatrix

subroutine create_btensor(btensor,diml,dimm,dimn)
  implicit none
  integer,intent(in) :: diml,dimm,dimn
  type(btensor_t),intent(inout) :: btensor    
  if((diml.ne.0).and. &
     (dimm.ne.0).and. &
     (dimn.ne.0)) then
      btensor%d(1) = diml
      btensor%d(2) = dimm
      btensor%d(3) = dimn
  else
    write(*,*) "create_btensor() error:"
    write(*,*) "Zero allocation size:",diml,dimm,dimn
    stop
  endif

  if(.not.allocated(btensor%p)) then 
    allocate(btensor%p(btensor%d(1),btensor%d(2),btensor%d(3)))
  else
    write(*,*) "create_btensor() error:"
    write(*,*) "rbtensor_t has been allocated"
    stop
  endif
end subroutine create_btensor

subroutine create_ivector(ivector,dimn)
  implicit none
  integer,intent(in) :: dimn
  type(ivector_t),intent(inout) :: ivector    
  if(dimn.ne.0) then
    ivector%d = dimn
  else
    write(*,*) "create_ivector() error:"
    write(*,*) "Zero allocation size:",dimn
    stop
  endif

  if(.not.allocated(ivector%p)) then 
    allocate(ivector%p(ivector%d))
  else
    write(*,*) "create_ivector() error:"
    write(*,*) "ivector_t has been allocated"
    stop
  endif
end subroutine create_ivector

subroutine create_imatrix(imatrix,dimm,dimn)
  implicit none
  integer,intent(in) :: dimm,dimn
  type(imatrix_t),intent(inout) :: imatrix    
  if((dimm.ne.0).and. &
     (dimn.ne.0)) then
      imatrix%d(1) = dimm
      imatrix%d(2) = dimn
  else
    write(*,*) "create_imatrix() error:"
    write(*,*) "Zero allocation size:",dimm,dimn
    stop
  endif

  if(.not.allocated(imatrix%p)) then 
    allocate(imatrix%p(imatrix%d(1),imatrix%d(2)))
  else
    write(*,*) "create_imatrix() error:"
    write(*,*) "imatrix_t has been allocated"
    stop
  endif
end subroutine create_imatrix

subroutine create_itensor(itensor,diml,dimm,dimn)
  implicit none
  integer,intent(in) :: diml,dimm,dimn
  type(itensor_t),intent(inout) :: itensor    
  if((diml.ne.0).and. &
     (dimm.ne.0).and. &
     (dimn.ne.0)) then
      itensor%d(1) = diml
      itensor%d(2) = dimm
      itensor%d(3) = dimn
  else
    write(*,*) "create_itensor() error:"
    write(*,*) "Zero allocation size:",diml,dimm,dimn
    stop
  endif

  if(.not.allocated(itensor%p)) then 
    allocate(itensor%p(itensor%d(1),itensor%d(2),itensor%d(3)))
  else
    write(*,*) "create_itensor() error:"
    write(*,*) "itensor_t has been allocated"
    stop
  endif
end subroutine create_itensor

subroutine create_rvector(rvector,dimn)
  implicit none
  integer,intent(in) :: dimn
  type(rvector_t),intent(inout) :: rvector    
  if(dimn.ne.0) then
    rvector%d = dimn
  else
    write(*,*) "create_rvector() error:"
    write(*,*) "Zero allocation size:",dimn
    stop
  endif

  if(.not.allocated(rvector%p)) then 
    allocate(rvector%p(rvector%d))
  else
    write(*,*) "create_rvector() error:"
    write(*,*) "rvector_t has been allocated"
    stop
  endif
end subroutine create_rvector

subroutine create_rmatrix(rmatrix,dimm,dimn)
  implicit none
  integer,intent(in) :: dimm,dimn
  type(rmatrix_t),intent(inout) :: rmatrix    
  if((dimm.ne.0).and. &
     (dimn.ne.0)) then
      rmatrix%d(1) = dimm
      rmatrix%d(2) = dimn
  else
    write(*,*) "create_rmatrix() error:"
    write(*,*) "Zero allocation size:",dimm,dimn
    stop
  endif

  if(.not.allocated(rmatrix%p)) then 
    allocate(rmatrix%p(rmatrix%d(1),rmatrix%d(2)))
  else
    write(*,*) "create_rmatrix() error:"
    write(*,*) "rmatrix_t has been allocated"
    stop
  endif
end subroutine create_rmatrix

subroutine create_rtensor(rtensor,diml,dimm,dimn)
  implicit none
  integer,intent(in) :: diml,dimm,dimn
  type(rtensor_t),intent(inout) :: rtensor    
  if((diml.ne.0).and. &
     (dimm.ne.0).and. &
     (dimn.ne.0)) then
      rtensor%d(1) = diml
      rtensor%d(2) = dimm
      rtensor%d(3) = dimn
  else
    write(*,*) "create_rtensor() error:"
    write(*,*) "Zero allocation size:",diml,dimm,dimn
    stop
  endif

  if(.not.allocated(rtensor%p)) then 
    allocate(rtensor%p(rtensor%d(1),rtensor%d(2),rtensor%d(3)))
  else
    write(*,*) "create_rtensor() error:"
    write(*,*) "rtensor_t has been allocated"
    stop
  endif
end subroutine create_rtensor

subroutine destroy_bvector(bvector)
  implicit none
  type(bvector_t),intent(inout) :: bvector   
  if(allocated(bvector%p)) then 
    deallocate(bvector%p)
  endif
end subroutine destroy_bvector

subroutine destroy_bmatrix(bmatrix)
  implicit none
  type(bmatrix_t),intent(inout) :: bmatrix  
  if(allocated(bmatrix%p)) then 
    deallocate(bmatrix%p)
  endif
end subroutine destroy_bmatrix

subroutine destroy_btensor(btensor)
  implicit none
  type(btensor_t),intent(inout) :: btensor  
  if(allocated(btensor%p)) then 
    deallocate(btensor%p)
  endif
end subroutine destroy_btensor

subroutine destroy_ivector(ivector)
  implicit none
  type(ivector_t),intent(inout) :: ivector   
  if(allocated(ivector%p)) then 
    deallocate(ivector%p)
  endif
end subroutine destroy_ivector

subroutine destroy_imatrix(imatrix)
  implicit none
  type(imatrix_t),intent(inout) :: imatrix  
  if(allocated(imatrix%p)) then 
    deallocate(imatrix%p)
  endif
end subroutine destroy_imatrix

subroutine destroy_itensor(itensor)
  implicit none
  type(itensor_t),intent(inout) :: itensor  
  if(allocated(itensor%p)) then 
    deallocate(itensor%p)
  endif
end subroutine destroy_itensor

subroutine destroy_rvector(rvector)
  implicit none
  type(rvector_t),intent(inout) :: rvector   
  if(allocated(rvector%p)) then 
    deallocate(rvector%p)
  endif
end subroutine destroy_rvector

subroutine destroy_rmatrix(rmatrix)
  implicit none
  type(rmatrix_t),intent(inout) :: rmatrix  
  if(allocated(rmatrix%p)) then 
    deallocate(rmatrix%p)
  endif
end subroutine destroy_rmatrix

subroutine destroy_rtensor(rtensor)
  implicit none
  type(rtensor_t),intent(inout) :: rtensor  
  if(allocated(rtensor%p)) then 
    deallocate(rtensor%p)
  endif
end subroutine destroy_rtensor

subroutine create_bscalarfield(bfield,s)
  implicit none
  type(fscalarshape_t),intent(in) :: s
  type(bscalarfield_t),intent(inout) :: bfield    
  integer :: err
  if((s%d(1).ne.0).and. &
      (s%d(2).ne.0)) then
    bfield%s = s
  else
    write(*,*) "Allocation error:"
    write(*,*) "Zero allocation size:",s%d(:)
    stop
  endif

  if(.not.allocated(bfield%p)) then 
    allocate(bfield%p(s%d(1),s%d(2)),stat = err)
    if(err.ne.0) write(*,*) "Scalarfield allocation failed."
  else
    write(*,*) "Allocation error:"
    write(*,*) "Field has been allocated"
    stop
  endif
end subroutine create_bscalarfield

subroutine create_bvectorfield(bfield,s)
  implicit none
  type(fvectorshape_t),intent(in) :: s
  type(bvectorfield_t),intent(inout) :: bfield    
  if((s%d(1).ne.0).and. &
      (s%d(2).ne.0).and. &
      (s%d(3).ne.0)) then
    bfield%s = s
  else
    write(*,*) "Allocation error:"
    write(*,*) "Zero allocation size:",s%d(:)
    stop
  endif
  
  if(.not.allocated(bfield%p)) then 
    allocate(bfield%p(s%d(1),s%d(2),s%d(3)))
  else
    write(*,*) "Allocation error:"
    write(*,*) "Field has been allocated"
    stop
  endif
end subroutine create_bvectorfield

subroutine create_bmatrixfield(bfield,s)
  implicit none
  type(fmatrixshape_t),intent(in) :: s
  type(bmatrixfield_t),intent(inout) :: bfield    
  if((s%d(1).ne.0).and. &
      (s%d(2).ne.0).and. &
      (s%d(3).ne.0).and. &
      (s%d(4).ne.0)) then
    bfield%s = s
  else
    write(*,*) "Allocation error:"
    write(*,*) "Zero allocation size:",s%d(:)
    stop
  endif
  
  if(.not.allocated(bfield%p)) then 
    allocate(bfield%p(s%d(1),s%d(2),s%d(3),s%d(4)))
  else
    write(*,*) "Allocation error:"
    write(*,*) "Field has been allocated"
    stop
  endif
end subroutine create_bmatrixfield

subroutine create_btensorfield(bfield,s)
  implicit none
  type(ftensorshape_t),intent(in) :: s
  type(btensorfield_t),intent(inout) :: bfield      
  if((s%d(1).ne.0).and. &
      (s%d(2).ne.0).and. &
      (s%d(3).ne.0).and. &
      (s%d(4).ne.0).and. &
      (s%d(5).ne.0)) then
    bfield%s = s
  else
    write(*,*) "Allocation error:"
    write(*,*) "Zero allocation size:",s%d(:)
    stop
  endif
  
  if(.not.allocated(bfield%p)) then 
    allocate(bfield%p(s%d(1),s%d(2),s%d(3),s%d(4),s%d(5)))
  else
    write(*,*) "Allocation error:"
    write(*,*) "Field has been allocated"
    stop
  endif
end subroutine create_btensorfield

subroutine create_iscalarfield(ifield,s)
  implicit none
  type(fscalarshape_t),intent(in) :: s
  type(iscalarfield_t),intent(inout) :: ifield    
  integer :: err
  if((s%d(1).ne.0).and. &
      (s%d(2).ne.0)) then
    ifield%s = s
  else
    write(*,*) "Allocation error:"
    write(*,*) "Zero allocation size:",s%d(:)
    stop
  endif

  if(.not.allocated(ifield%p)) then 
    allocate(ifield%p(s%d(1),s%d(2)),stat = err)
    if(err.ne.0) write(*,*) "Scalarfield allocation failed."
  else
    write(*,*) "Allocation error:"
    write(*,*) "Field has been allocated"
    stop
  endif
end subroutine create_iscalarfield

subroutine create_ivectorfield(ifield,s)
  implicit none
  type(fvectorshape_t),intent(in) :: s
  type(ivectorfield_t),intent(inout) :: ifield    
  if((s%d(1).ne.0).and. &
      (s%d(2).ne.0).and. &
      (s%d(3).ne.0)) then
    ifield%s = s
  else
    write(*,*) "Allocation error:"
    write(*,*) "Zero allocation size:",s%d(:)
    stop
  endif
  
  if(.not.allocated(ifield%p)) then 
    allocate(ifield%p(s%d(1),s%d(2),s%d(3)))
  else
    write(*,*) "Allocation error:"
    write(*,*) "Field has been allocated"
    stop
  endif
end subroutine create_ivectorfield

subroutine create_imatrixfield(ifield,s)
  implicit none
  type(fmatrixshape_t),intent(in) :: s
  type(imatrixfield_t),intent(inout) :: ifield    
  if((s%d(1).ne.0).and. &
      (s%d(2).ne.0).and. &
      (s%d(3).ne.0).and. &
      (s%d(4).ne.0)) then
    ifield%s = s
  else
    write(*,*) "Allocation error:"
    write(*,*) "Zero allocation size:",s%d(:)
    stop
  endif
  
  if(.not.allocated(ifield%p)) then 
    allocate(ifield%p(s%d(1),s%d(2),s%d(3),s%d(4)))
  else
    write(*,*) "Allocation error:"
    write(*,*) "Field has been allocated"
    stop
  endif
end subroutine create_imatrixfield

subroutine create_itensorfield(ifield,s)
  implicit none
  type(ftensorshape_t),intent(in) :: s
  type(itensorfield_t),intent(inout) :: ifield      
  if((s%d(1).ne.0).and. &
      (s%d(2).ne.0).and. &
      (s%d(3).ne.0).and. &
      (s%d(4).ne.0).and. &
      (s%d(5).ne.0)) then
    ifield%s = s
  else
    write(*,*) "Allocation error:"
    write(*,*) "Zero allocation size:",s%d(:)
    stop
  endif
  
  if(.not.allocated(ifield%p)) then 
    allocate(ifield%p(s%d(1),s%d(2),s%d(3),s%d(4),s%d(5)))
  else
    write(*,*) "Allocation error:"
    write(*,*) "Field has been allocated"
    stop
  endif
end subroutine create_itensorfield

subroutine create_rscalarfield(rfield,s)
  implicit none
  type(fscalarshape_t),intent(in) :: s
  type(rscalarfield_t),intent(inout) :: rfield    
  integer :: err
  if((s%d(1).ne.0).and. &
      (s%d(2).ne.0)) then
    rfield%s = s
  else
    write(*,*) "Allocation error:"
    write(*,*) "Zero allocation size:",s%d(:)
    stop
  endif

  if(.not.allocated(rfield%p)) then 
    allocate(rfield%p(s%d(1),s%d(2)),stat = err)
    if(err.ne.0) write(*,*) "Scalarfield allocation failed."
  else
    write(*,*) "Allocation error:"
    write(*,*) "Field has been allocated"
    stop
  endif
end subroutine create_rscalarfield

subroutine create_rvectorfield(rfield,s)
  implicit none
  type(fvectorshape_t),intent(in) :: s
  type(rvectorfield_t),intent(inout) :: rfield    
  if((s%d(1).ne.0).and. &
      (s%d(2).ne.0).and. &
      (s%d(3).ne.0)) then
    rfield%s = s
  else
    write(*,*) "Allocation error:"
    write(*,*) "Zero allocation size:",s%d(:)
    stop
  endif
  
  if(.not.allocated(rfield%p)) then 
    allocate(rfield%p(s%d(1),s%d(2),s%d(3)))
  else
    write(*,*) "Allocation error:"
    write(*,*) "Field has been allocated"
    stop
  endif
end subroutine create_rvectorfield

subroutine create_rmatrixfield(rfield,s)
  implicit none
  type(fmatrixshape_t),intent(in) :: s
  type(rmatrixfield_t),intent(inout) :: rfield    
  if((s%d(1).ne.0).and. &
      (s%d(2).ne.0).and. &
      (s%d(3).ne.0).and. &
      (s%d(4).ne.0)) then
    rfield%s = s
  else
    write(*,*) "Allocation error:"
    write(*,*) "Zero allocation size:",s%d(:)
    stop
  endif
  
  if(.not.allocated(rfield%p)) then 
    allocate(rfield%p(s%d(1),s%d(2),s%d(3),s%d(4)))
  else
    write(*,*) "Allocation error:"
    write(*,*) "Field has been allocated"
    stop
  endif
end subroutine create_rmatrixfield

subroutine create_rtensorfield(rfield,s)
  implicit none
  type(ftensorshape_t),intent(in) :: s
  type(rtensorfield_t),intent(inout) :: rfield      
  if((s%d(1).ne.0).and. &
      (s%d(2).ne.0).and. &
      (s%d(3).ne.0).and. &
      (s%d(4).ne.0).and. &
      (s%d(5).ne.0)) then
    rfield%s = s
  else
    write(*,*) "Allocation error:"
    write(*,*) "Zero allocation size:",s%d(:)
    stop
  endif
  
  if(.not.allocated(rfield%p)) then 
    allocate(rfield%p(s%d(1),s%d(2),s%d(3),s%d(4),s%d(5)))
  else
    write(*,*) "Allocation error:"
    write(*,*) "Field has been allocated"
    stop
  endif
end subroutine create_rtensorfield

subroutine destroy_bscalarfield(bfield)
  implicit none
  type(bscalarfield_t),intent(inout) :: bfield     
  if(allocated(bfield%p)) then 
    deallocate(bfield%p)
  endif
end subroutine destroy_bscalarfield

subroutine destroy_bvectorfield(bfield)
  implicit none
  type(bvectorfield_t),intent(inout) :: bfield     
  if(allocated(bfield%p)) then 
    deallocate(bfield%p)
  endif
end subroutine destroy_bvectorfield

subroutine destroy_bmatrixfield(bfield)
  implicit none
  type(bmatrixfield_t),intent(inout) :: bfield      
  if(allocated(bfield%p)) then 
    deallocate(bfield%p)
  endif
end subroutine destroy_bmatrixfield

subroutine destroy_btensorfield(bfield)
  implicit none
  type(btensorfield_t),intent(inout) :: bfield      
  if(allocated(bfield%p)) then 
    deallocate(bfield%p)
  endif
end subroutine destroy_btensorfield

subroutine destroy_iscalarfield(ifield)
  implicit none
  type(iscalarfield_t),intent(inout) :: ifield     
  if(allocated(ifield%p)) then 
    deallocate(ifield%p)
  endif
end subroutine destroy_iscalarfield

subroutine destroy_ivectorfield(ifield)
  implicit none
  type(ivectorfield_t),intent(inout) :: ifield     
  if(allocated(ifield%p)) then 
    deallocate(ifield%p)
  endif
end subroutine destroy_ivectorfield

subroutine destroy_imatrixfield(ifield)
  implicit none
  type(imatrixfield_t),intent(inout) :: ifield      
  if(allocated(ifield%p)) then 
    deallocate(ifield%p)
  endif
end subroutine destroy_imatrixfield

subroutine destroy_itensorfield(ifield)
  implicit none
  type(itensorfield_t),intent(inout) :: ifield      
  if(allocated(ifield%p)) then 
    deallocate(ifield%p)
  endif
end subroutine destroy_itensorfield

subroutine destroy_rscalarfield(rfield)
  implicit none
  type(rscalarfield_t),intent(inout) :: rfield     
  if(allocated(rfield%p)) then 
    deallocate(rfield%p)
  endif
end subroutine destroy_rscalarfield

subroutine destroy_rvectorfield(rfield)
  implicit none
  type(rvectorfield_t),intent(inout) :: rfield     
  if(allocated(rfield%p)) then 
    deallocate(rfield%p)
  endif
end subroutine destroy_rvectorfield

subroutine destroy_rmatrixfield(rfield)
  implicit none
  type(rmatrixfield_t),intent(inout) :: rfield      
  if(allocated(rfield%p)) then 
    deallocate(rfield%p)
  endif
end subroutine destroy_rmatrixfield

subroutine destroy_rtensorfield(rfield)
  implicit none
  type(rtensorfield_t),intent(inout) :: rfield      
  if(allocated(rfield%p)) then 
    deallocate(rfield%p)
  endif
end subroutine destroy_rtensorfield

subroutine fscalarshape_to_fscalarshape(left,right)
  implicit none
  type(fscalarshape_t),intent(in) :: right
  type(fscalarshape_t),intent(inout) :: left
  left%d(:) = right%d(:)
end subroutine 

subroutine fvectorshape_to_fvectorshape(left,right)
  implicit none
  type(fvectorshape_t),intent(in) :: right
  type(fvectorshape_t),intent(inout) :: left
  left%d(:) = right%d(:)
end subroutine 

subroutine fmatrixshape_to_fmatrixshape(left,right)
  implicit none
  type(fmatrixshape_t),intent(in) :: right
  type(fmatrixshape_t),intent(inout) :: left
  left%d(:) = right%d(:)
end subroutine 

subroutine ftensorshape_to_ftensorshape(left,right)
  implicit none
  type(ftensorshape_t),intent(in) :: right
  type(ftensorshape_t),intent(inout) :: left
  left%d(:) = right%d(:)
end subroutine 

function fscalarshape_eq_fscalarshape(left,right)
  implicit none
  type(fscalarshape_t),intent(in) :: left,right
  logical :: fscalarshape_eq_fscalarshape
  if((left%d(1).eq.right%d(1)) .and. &
  (left%d(2).eq.right%d(2))) then
    fscalarshape_eq_fscalarshape = .true.
  else
    fscalarshape_eq_fscalarshape = .false.
  endif
end function 

function fvectorshape_eq_fvectorshape(left,right)
  implicit none
  type(fvectorshape_t),intent(in) :: left,right
  logical :: fvectorshape_eq_fvectorshape
  if((left%d(1).eq.right%d(1)) .and. &
  (left%d(2).eq.right%d(2)) .and. &
  (left%d(3).eq.right%d(3))) then
    fvectorshape_eq_fvectorshape = .true.
  else
    fvectorshape_eq_fvectorshape = .false.
  endif
end function 

function fmatrixshape_eq_fmatrixshape(left,right)
  implicit none
  type(fmatrixshape_t),intent(in) :: left,right
  logical :: fmatrixshape_eq_fmatrixshape
  if((left%d(1).eq.right%d(1)) .and. &
  (left%d(2).eq.right%d(2)) .and. &
  (left%d(3).eq.right%d(3)) .and. &
  (left%d(4).eq.right%d(4))) then
    fmatrixshape_eq_fmatrixshape = .true.
  else
    fmatrixshape_eq_fmatrixshape = .false.
  endif
end function 

function ftensorshape_eq_ftensorshape(left,right)
  implicit none
  type(ftensorshape_t),intent(in) :: left,right
  logical :: ftensorshape_eq_ftensorshape
  if((left%d(1).eq.right%d(1)) .and. &
  (left%d(2).eq.right%d(2)) .and. &
  (left%d(3).eq.right%d(3)) .and. &
  (left%d(4).eq.right%d(4)) .and. &
  (left%d(5).eq.right%d(5))) then
    ftensorshape_eq_ftensorshape = .true.
  else
    ftensorshape_eq_ftensorshape = .false.
  endif
end function ftensorshape_eq_ftensorshape

function set_fscalarshape(dim1,dim2)
  implicit none
  integer,intent(in) :: dim1,dim2
  type(fscalarshape_t) :: set_fscalarshape
  set_fscalarshape%d(1) = dim1
  set_fscalarshape%d(2) = dim2
end function set_fscalarshape

function set_fvectorshape(dim1,dim2,dim3)
  implicit none
  integer,intent(in) :: dim1,dim2,dim3
  type(fvectorshape_t) :: set_fvectorshape
  set_fvectorshape%d(1) = dim1
  set_fvectorshape%d(2) = dim2
  set_fvectorshape%d(3) = dim3
end function set_fvectorshape

function set_fmatrixshape(dim1,dim2,dim3,dim4)
  implicit none
  integer,intent(in) :: dim1,dim2,dim3,dim4
  type(fmatrixshape_t) :: set_fmatrixshape
  set_fmatrixshape%d(1) = dim1
  set_fmatrixshape%d(2) = dim2
  set_fmatrixshape%d(3) = dim3
  set_fmatrixshape%d(4) = dim4
end function set_fmatrixshape

function set_ftensorshape(dim1,dim2,dim3,dim4,dim5)
  implicit none
  integer,intent(in) :: dim1,dim2,dim3,dim4,dim5
  type(ftensorshape_t) :: set_ftensorshape
  set_ftensorshape%d(1) = dim1
  set_ftensorshape%d(2) = dim2
  set_ftensorshape%d(3) = dim3
  set_ftensorshape%d(4) = dim4
  set_ftensorshape%d(5) = dim5
end function set_ftensorshape

subroutine bvector_to_bvector(left,right)
  implicit none
  type(bvector_t),intent(in) :: right
  type(bvector_t),intent(inout) :: left
  integer :: mindim
  integer :: i
  mindim = min(left%d,right%d)
  do i = 1, mindim
    left%p(:) = right%p(:)
  enddo
end subroutine 

subroutine bmatrix_to_bmatrix(left,right)
  implicit none
  type(bmatrix_t),intent(in) :: right
  type(bmatrix_t),intent(inout) :: left
  integer :: mindimn,mindimm
  integer :: i,j
  mindimn = min(left%d(1),right%d(1))
  mindimm = min(left%d(2),right%d(2))
  do j = 1, mindimm
    do i = 1, mindimn
      left%p(i,j) = right%p(i,j)
    enddo
  enddo
end subroutine 

subroutine btensor_to_btensor(left,right)
  implicit none
  type(btensor_t),intent(in) :: right
  type(btensor_t),intent(inout) :: left
  integer :: mindimn,mindimm,mindiml
  integer :: i,j,k
  mindimn = min(left%d(1),right%d(1))
  mindimm = min(left%d(2),right%d(2))
  mindiml = min(left%d(3),right%d(3))
  do k = 1, mindiml
    do j = 1, mindimm
      do i = 1, mindimn
        left%p(i,j,k) = right%p(i,j,k)
      enddo
    enddo
  enddo
end subroutine 

subroutine ivector_to_ivector(left,right)
  implicit none
  type(ivector_t),intent(in) :: right
  type(ivector_t),intent(inout) :: left
  integer :: mindim
  integer :: i
  mindim = min(left%d,right%d)
  do i = 1, mindim
    left%p(:) = right%p(:)
  enddo
end subroutine 

subroutine imatrix_to_imatrix(left,right)
  implicit none
  type(imatrix_t),intent(in) :: right
  type(imatrix_t),intent(inout) :: left
  integer :: mindimn,mindimm
  integer :: i,j
  mindimn = min(left%d(1),right%d(1))
  mindimm = min(left%d(2),right%d(2))
  do j = 1, mindimm
    do i = 1, mindimn
      left%p(i,j) = right%p(i,j)
    enddo
  enddo
end subroutine 

subroutine itensor_to_itensor(left,right)
  implicit none
  type(itensor_t),intent(in) :: right
  type(itensor_t),intent(inout) :: left
  integer :: mindimn,mindimm,mindiml
  integer :: i,j,k
  mindimn = min(left%d(1),right%d(1))
  mindimm = min(left%d(2),right%d(2))
  mindiml = min(left%d(3),right%d(3))
  do k = 1, mindiml
    do j = 1, mindimm
      do i = 1, mindimn
        left%p(i,j,k) = right%p(i,j,k)
      enddo
    enddo
  enddo
end subroutine 

subroutine rvector_to_rvector(left,right)
  implicit none
  type(rvector_t),intent(in) :: right
  type(rvector_t),intent(inout) :: left
  integer :: mindim
  integer :: i
  mindim = min(left%d,right%d)
  do i = 1, mindim
    left%p(:) = right%p(:)
  enddo
end subroutine 

subroutine rmatrix_to_rmatrix(left,right)
  implicit none
  type(rmatrix_t),intent(in) :: right
  type(rmatrix_t),intent(inout) :: left
  integer :: mindimn,mindimm
  integer :: i,j
  mindimn = min(left%d(1),right%d(1))
  mindimm = min(left%d(2),right%d(2))
  do j = 1, mindimm
    do i = 1, mindimn
      left%p(i,j) = right%p(i,j)
    enddo
  enddo
end subroutine 

subroutine rtensor_to_rtensor(left,right)
  implicit none
  type(rtensor_t),intent(in) :: right
  type(rtensor_t),intent(inout) :: left
  integer :: mindimn,mindimm,mindiml
  integer :: i,j,k
  mindimn = min(left%d(1),right%d(1))
  mindimm = min(left%d(2),right%d(2))
  mindiml = min(left%d(3),right%d(3))
  do k = 1, mindiml
    do j = 1, mindimm
      do i = 1, mindimn
        left%p(i,j,k) = right%p(i,j,k)
      enddo
    enddo
  enddo
end subroutine 

subroutine bscalarfield_to_bscalar(scalarfield,scalar)
  implicit none
  logical,intent(in) :: scalar
  type(bscalarfield_t),intent(inout) :: scalarfield
  integer :: i,j
  do j = 1, scalarfield%s%d(2)
    do i = 1, scalarfield%s%d(1)
      scalarfield%p(i,j) = scalar
    enddo
  enddo
end subroutine 

subroutine bvectorfield_to_bvector(vectorfield,vector)
  implicit none
  type(bvector_t),intent(in) :: vector
  type(bvectorfield_t),intent(inout) :: vectorfield
  integer :: i,j
  if(vectorfield%s%d(1).eq.vector%d) then
    do j = 1, vectorfield%s%d(3)
      do i = 1, vectorfield%s%d(2)
        vectorfield%p(:,i,j) = vector%p(:)
      enddo
    enddo
  else
    write(*,*) "operator vectorfield_to_vector error:"
    write(*,*) "Unequal rank size:",vectorfield%s%d(1),vector%d
  endif
end subroutine 

subroutine bmatrixfield_to_bmatrix(matrixfield,matrix)
  implicit none
  type(bmatrix_t),intent(in) :: matrix
  type(bmatrixfield_t),intent(inout) :: matrixfield
  integer :: i,j
  if((matrixfield%s%d(1).eq.matrix%d(1)).or. &
      (matrixfield%s%d(2).eq.matrix%d(2))) then
    do j = 1, matrixfield%s%d(4)
      do i = 1, matrixfield%s%d(3)
        matrixfield%p(:,:,i,j) = matrix%p(:,:)
      enddo
    enddo
  else
    write(*,*) "operator matrixfield_to_matrix error:"
    write(*,*) "Unequal rank size:",matrixfield%s%d(1:2),matrix%d
  endif
end subroutine 

subroutine btensorfield_to_btensor(tensorfield,tensor)
  implicit none
  type(btensor_t),intent(in) :: tensor
  type(btensorfield_t),intent(inout) :: tensorfield
  integer :: i,j
  if((tensorfield%s%d(1).eq.tensor%d(1)).or. &
      (tensorfield%s%d(2).eq.tensor%d(2)).or. &
      (tensorfield%s%d(3).eq.tensor%d(3))) then
    do j = 1, tensorfield%s%d(5)
      do i = 1, tensorfield%s%d(4)
        tensorfield%p(:,:,:,i,j) = tensor%p(:,:,:)
      enddo
    enddo
  else
    write(*,*) "operator tensorfield_to_tensor error:"
    write(*,*) "Unequal rank size:",tensorfield%s%d(1:3),tensor%d
  endif
end subroutine 

subroutine iscalarfield_to_iscalar(scalarfield,scalar)
  implicit none
  integer,intent(in) :: scalar
  type(iscalarfield_t),intent(inout) :: scalarfield
  integer :: i,j
  do j = 1, scalarfield%s%d(2)
    do i = 1, scalarfield%s%d(1)
      scalarfield%p(i,j) = scalar
    enddo
  enddo
end subroutine 

subroutine ivectorfield_to_ivector(vectorfield,vector)
  implicit none
  type(ivector_t),intent(in) :: vector
  type(ivectorfield_t),intent(inout) :: vectorfield
  integer :: i,j
  if(vectorfield%s%d(1).eq.vector%d) then
    do j = 1, vectorfield%s%d(3)
      do i = 1, vectorfield%s%d(2)
        vectorfield%p(:,i,j) = vector%p(:)
      enddo
    enddo
  else
    write(*,*) "operator vectorfield_to_vector error:"
    write(*,*) "Unequal rank size:",vectorfield%s%d(1),vector%d
  endif
end subroutine 

subroutine imatrixfield_to_imatrix(matrixfield,matrix)
  implicit none
  type(imatrix_t),intent(in) :: matrix
  type(imatrixfield_t),intent(inout) :: matrixfield
  integer :: i,j
  if((matrixfield%s%d(1).eq.matrix%d(1)).or. &
      (matrixfield%s%d(2).eq.matrix%d(2))) then
    do j = 1, matrixfield%s%d(4)
      do i = 1, matrixfield%s%d(3)
        matrixfield%p(:,:,i,j) = matrix%p(:,:)
      enddo
    enddo
  else
    write(*,*) "operator matrixfield_to_matrix error:"
    write(*,*) "Unequal rank size:",matrixfield%s%d(1:2),matrix%d
  endif
end subroutine 

subroutine itensorfield_to_itensor(tensorfield,tensor)
  implicit none
  type(itensor_t),intent(in) :: tensor
  type(itensorfield_t),intent(inout) :: tensorfield
  integer :: i,j
  if((tensorfield%s%d(1).eq.tensor%d(1)).or. &
      (tensorfield%s%d(2).eq.tensor%d(2)).or. &
      (tensorfield%s%d(3).eq.tensor%d(3))) then
    do j = 1, tensorfield%s%d(5)
      do i = 1, tensorfield%s%d(4)
        tensorfield%p(:,:,:,i,j) = tensor%p(:,:,:)
      enddo
    enddo
  else
    write(*,*) "operator tensorfield_to_tensor error:"
    write(*,*) "Unequal rank size:",tensorfield%s%d(1:3),tensor%d
  endif
end subroutine 

subroutine rscalarfield_to_rscalar(scalarfield,scalar)
  implicit none
  double precision,intent(in) :: scalar
  type(rscalarfield_t),intent(inout) :: scalarfield
  integer :: i,j
  do j = 1, scalarfield%s%d(2)
    do i = 1, scalarfield%s%d(1)
      scalarfield%p(i,j) = scalar
    enddo
  enddo
end subroutine 

subroutine rvectorfield_to_rvector(vectorfield,vector)
  implicit none
  type(rvector_t),intent(in) :: vector
  type(rvectorfield_t),intent(inout) :: vectorfield
  integer :: i,j
  if(vectorfield%s%d(1).eq.vector%d) then
    do j = 1, vectorfield%s%d(3)
      do i = 1, vectorfield%s%d(2)
        vectorfield%p(:,i,j) = vector%p(:)
      enddo
    enddo
  else
    write(*,*) "operator vectorfield_to_vector error:"
    write(*,*) "Unequal rank size:",vectorfield%s%d(1),vector%d
  endif
end subroutine 

subroutine rmatrixfield_to_rmatrix(matrixfield,matrix)
  implicit none
  type(rmatrix_t),intent(in) :: matrix
  type(rmatrixfield_t),intent(inout) :: matrixfield
  integer :: i,j
  if((matrixfield%s%d(1).eq.matrix%d(1)).or. &
      (matrixfield%s%d(2).eq.matrix%d(2))) then
    do j = 1, matrixfield%s%d(4)
      do i = 1, matrixfield%s%d(3)
        matrixfield%p(:,:,i,j) = matrix%p(:,:)
      enddo
    enddo
  else
    write(*,*) "operator matrixfield_to_matrix error:"
    write(*,*) "Unequal rank size:",matrixfield%s%d(1:2),matrix%d
  endif
end subroutine 

subroutine rtensorfield_to_rtensor(tensorfield,tensor)
  implicit none
  type(rtensor_t),intent(in) :: tensor
  type(rtensorfield_t),intent(inout) :: tensorfield
  integer :: i,j
  if((tensorfield%s%d(1).eq.tensor%d(1)).or. &
      (tensorfield%s%d(2).eq.tensor%d(2)).or. &
      (tensorfield%s%d(3).eq.tensor%d(3))) then
    do j = 1, tensorfield%s%d(5)
      do i = 1, tensorfield%s%d(4)
        tensorfield%p(:,:,:,i,j) = tensor%p(:,:,:)
      enddo
    enddo
  else
    write(*,*) "operator tensorfield_to_tensor error:"
    write(*,*) "Unequal rank size:",tensorfield%s%d(1:3),tensor%d
  endif
end subroutine 

subroutine bscalarfield_to_bscalarfield(left,right)
  implicit none
  type(bscalarfield_t),intent(in) :: right
  type(bscalarfield_t),intent(inout) :: left
  integer :: i,j
  if(left%s.eq.right%s) then
    do j = 1, left%s%d(2)
      do i = 1, left%s%d(1)
        left%p(i,j) = right%p(i,j)
      enddo
    enddo
  else
    write(*,*) "operator bscalarfield_to_bscalarfield() error:"
    write(*,*) "Unequal rank :"
    write(*,*) "Left size:",left%s%d
    write(*,*) "Rightsize:",right%s%d
  endif
end subroutine 

subroutine bvectorfield_to_bvectorfield(left,right)
  implicit none
  type(bvectorfield_t),intent(in) :: right
  type(bvectorfield_t),intent(inout) :: left
  integer :: i,j
  if(left%s.eq.right%s) then
    do j = 1, left%s%d(3)
      do i = 1, left%s%d(2)
        left%p(:,i,j) = right%p(:,i,j)
      enddo
    enddo
  else
    write(*,*) "operator bvectorfield_to_bvectorfield() error:"
    write(*,*) "Unequal rank :"
    write(*,*) "Left size:",left%s%d
    write(*,*) "Rightsize:",right%s%d
  endif
end subroutine 

subroutine bmatrixfield_to_bmatrixfield(left,right)
  implicit none
  type(bmatrixfield_t),intent(in) :: right
  type(bmatrixfield_t),intent(inout) :: left
  integer :: i,j
  if(left%s.eq.right%s) then
    do j = 1, left%s%d(4)
      do i = 1, left%s%d(3)
        left%p(:,:,i,j) = right%p(:,:,i,j)
      enddo
    enddo
  else
    write(*,*) "operator bmatrixfield_to_bmatrixfield() error:"
    write(*,*) "Unequal rank :"
    write(*,*) "Left size:",left%s%d
    write(*,*) "Rightsize:",right%s%d
  endif
end subroutine 

subroutine btensorfield_to_btensorfield(left,right)
  implicit none
  type(btensorfield_t),intent(in) :: right
  type(btensorfield_t),intent(inout) :: left
  integer :: i,j
  if(left%s.eq.right%s) then
    do j = 1, left%s%d(5)
      do i = 1, left%s%d(4)
        left%p(:,:,:,i,j) = right%p(:,:,:,i,j)
      enddo
    enddo
  else
    write(*,*) "operator btensorfield_to_btensorfield() error:"
    write(*,*) "Unequal rank :"
    write(*,*) "Left size:",left%s%d
    write(*,*) "Rightsize:",right%s%d
  endif
end subroutine 

subroutine iscalarfield_to_iscalarfield(left,right)
  implicit none
  type(iscalarfield_t),intent(in) :: right
  type(iscalarfield_t),intent(inout) :: left
  integer :: i,j
  if(left%s.eq.right%s) then
    do j = 1, left%s%d(2)
      do i = 1, left%s%d(1)
        left%p(i,j) = right%p(i,j)
      enddo
    enddo
  else
    write(*,*) "operator iscalarfield_to_iscalarfield() error:"
    write(*,*) "Unequal rank :"
    write(*,*) "Left size:",left%s%d
    write(*,*) "Rightsize:",right%s%d
  endif
end subroutine 

subroutine ivectorfield_to_ivectorfield(left,right)
  implicit none
  type(ivectorfield_t),intent(in) :: right
  type(ivectorfield_t),intent(inout) :: left
  integer :: i,j
  if(left%s.eq.right%s) then
    do j = 1, left%s%d(3)
      do i = 1, left%s%d(2)
        left%p(:,i,j) = right%p(:,i,j)
      enddo
    enddo
  else
    write(*,*) "operator ivectorfield_to_ivectorfield() error:"
    write(*,*) "Unequal rank :"
    write(*,*) "Left size:",left%s%d
    write(*,*) "Rightsize:",right%s%d
  endif
end subroutine 

subroutine imatrixfield_to_imatrixfield(left,right)
  implicit none
  type(imatrixfield_t),intent(in) :: right
  type(imatrixfield_t),intent(inout) :: left
  integer :: i,j
  if(left%s.eq.right%s) then
    do j = 1, left%s%d(4)
      do i = 1, left%s%d(3)
        left%p(:,:,i,j) = right%p(:,:,i,j)
      enddo
    enddo
  else
    write(*,*) "operator imatrixfield_to_imatrixfield() error:"
    write(*,*) "Unequal rank :"
    write(*,*) "Left size:",left%s%d
    write(*,*) "Rightsize:",right%s%d
  endif
end subroutine 

subroutine itensorfield_to_itensorfield(left,right)
  implicit none
  type(itensorfield_t),intent(in) :: right
  type(itensorfield_t),intent(inout) :: left
  integer :: i,j
  if(left%s.eq.right%s) then
    do j = 1, left%s%d(5)
      do i = 1, left%s%d(4)
        left%p(:,:,:,i,j) = right%p(:,:,:,i,j)
      enddo
    enddo
  else
    write(*,*) "operator itensorfield_to_itensorfield() error:"
    write(*,*) "Unequal rank :"
    write(*,*) "Left size:",left%s%d
    write(*,*) "Rightsize:",right%s%d
  endif
end subroutine 

subroutine rscalarfield_to_rscalarfield(left,right)
  implicit none
  type(rscalarfield_t),intent(in) :: right
  type(rscalarfield_t),intent(inout) :: left
  integer :: i,j
  if(left%s.eq.right%s) then
    do j = 1, left%s%d(2)
      do i = 1, left%s%d(1)
        left%p(i,j) = right%p(i,j)
      enddo
    enddo
  else
    write(*,*) "operator rscalarfield_to_rscalarfield() error:"
    write(*,*) "Unequal rank :"
    write(*,*) "Left size:",left%s%d
    write(*,*) "Rightsize:",right%s%d
  endif
end subroutine 

subroutine rvectorfield_to_rvectorfield(left,right)
  implicit none
  type(rvectorfield_t),intent(in) :: right
  type(rvectorfield_t),intent(inout) :: left
  integer :: i,j
  if(left%s.eq.right%s) then
    do j = 1, left%s%d(3)
      do i = 1, left%s%d(2)
        left%p(:,i,j) = right%p(:,i,j)
      enddo
    enddo
  else
    write(*,*) "operator rvectorfield_to_rvectorfield() error:"
    write(*,*) "Unequal rank :"
    write(*,*) "Left size:",left%s%d
    write(*,*) "Rightsize:",right%s%d
  endif
end subroutine 

subroutine rmatrixfield_to_rmatrixfield(left,right)
  implicit none
  type(rmatrixfield_t),intent(in) :: right
  type(rmatrixfield_t),intent(inout) :: left
  integer :: i,j
  if(left%s.eq.right%s) then
    do j = 1, left%s%d(4)
      do i = 1, left%s%d(3)
        left%p(:,:,i,j) = right%p(:,:,i,j)
      enddo
    enddo
  else
    write(*,*) "operator rmatrixfield_to_rmatrixfield() error:"
    write(*,*) "Unequal rank :"
    write(*,*) "Left size:",left%s%d
    write(*,*) "Rightsize:",right%s%d
  endif
end subroutine 

subroutine rtensorfield_to_rtensorfield(left,right)
  implicit none
  type(rtensorfield_t),intent(in) :: right
  type(rtensorfield_t),intent(inout) :: left
  integer :: i,j
  if(left%s.eq.right%s) then
    do j = 1, left%s%d(5)
      do i = 1, left%s%d(4)
        left%p(:,:,:,i,j) = right%p(:,:,:,i,j)
      enddo
    enddo
  else
    write(*,*) "operator rtensorfield_to_rtensorfield() error:"
    write(*,*) "Unequal rank :"
    write(*,*) "Left size:",left%s%d
    write(*,*) "Rightsize:",right%s%d
  endif
end subroutine 

end module