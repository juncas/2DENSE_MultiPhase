module boundary
  use data_type
implicit none
interface boundary_interface
  module procedure boundary_interface_scalar
  module procedure boundary_interface_vector
  module procedure boundary_interface_matrix
  module procedure boundary_interface_tensor
end interface boundary_interface

interface boundary_reflect
  module procedure boundary_reflect_scalar
  module procedure boundary_reflect_vector
  module procedure boundary_reflect_matrix
  module procedure boundary_reflect_tensor
end interface boundary_reflect

interface boundary_dirichlet
  module procedure boundary_dirichlet_scalar
  module procedure boundary_dirichlet_vector
  module procedure boundary_dirichlet_matrix
  module procedure boundary_dirichlet_tensor
end interface boundary_dirichlet

interface boundary_neumann
  module procedure boundary_neumann_scalar
  module procedure boundary_neumann_vector
  module procedure boundary_neumann_matrix
  module procedure boundary_neumann_tensor
end interface boundary_neumann

! General boundary condition indice
integer,parameter :: bc_idx_intf = 0
integer,parameter :: bc_idx_inlt = 1
integer,parameter :: bc_idx_otlt = 2
integer,parameter :: bc_idx_rflt = 3
integer,parameter :: bc_idx_spcf = 4

contains
subroutine boundary_interface_scalar(f_r,f_d,lo_r,hi_r,lo_d,hi_d,cnnct)
  implicit none
  type(rscalarfield_t),intent(inout) :: f_r
  type(rscalarfield_t),intent(inout) :: f_d
  type(fscalarshape_t),intent(in) :: lo_r,hi_r
  type(fscalarshape_t),intent(in) :: lo_d,hi_d
  integer,intent(in) :: cnnct
  integer :: irng_r,jrng_r
  integer :: irng_d,jrng_d
  integer :: i,j
  integer :: i0,j0
  integer :: i1,j1
  integer :: incr1,incr2

  irng_r = abs(lo_r%d(1) - hi_r%d(1)) + 1
  jrng_r = abs(lo_r%d(2) - hi_r%d(2)) + 1
  irng_d = abs(hi_d%d(1) - lo_d%d(1)) + 1
  jrng_d = abs(hi_d%d(2) - lo_d%d(2)) + 1

  if(lo_d%d(1).lt.hi_d%d(1)) then
    incr1 = 1
  else
    incr1 =-1
  endif

  if(lo_d%d(2).lt.hi_d%d(2)) then
    incr2 = 1
  else
    incr2 =-1
  endif

  if(cnnct.eq.1) then
    if((irng_r.ne.irng_d).or.(jrng_r.ne.jrng_d)) then
      write(*,*) "Error: boundary_interface_scalar()"
      write(*,*) "Incompatibal interface sizes:", &
                 irng_r,jrng_r,irng_d,jrng_d    
    endif
    do j = 0, jrng_r - 1
      do i = 0, irng_r - 1
        i0 = lo_r%d(1) + i
        j0 = lo_r%d(2) + j
        i1 = lo_d%d(1) + incr1*i
        j1 = lo_d%d(2) + incr2*j
        f_r%p(i0,j0) = f_d%p(i1,j1)
      enddo
    enddo
  elseif(cnnct.eq.2) then
    if((irng_r.ne.jrng_d).or.(jrng_r.ne.irng_d)) then
      write(*,*) "Error: boundary_interface_scalar()"
      write(*,*) "Incompatibal interface sizes:", &
                 irng_r,jrng_r,jrng_d,irng_d    
    endif
    do j = 0, jrng_r - 1
      do i = 0, irng_r - 1
        i0 = lo_r%d(1) + i
        j0 = lo_r%d(2) + j
        i1 = lo_d%d(1) + incr1*j
        j1 = lo_d%d(2) + incr2*i
        f_r%p(i0,j0) = f_d%p(i1,j1)
      enddo
    enddo
  else
    print*,"Error: boundary_interface_scalar()"
    print*,"no such cnnct:",cnnct
  endif
end subroutine boundary_interface_scalar  

subroutine boundary_interface_vector(f_r,f_d,lo_r,hi_r,lo_d,hi_d,cnnct)
  implicit none
  type(rvectorfield_t),intent(inout) :: f_r
  type(rvectorfield_t),intent(inout) :: f_d
  type(fscalarshape_t),intent(in) :: lo_r,hi_r
  type(fscalarshape_t),intent(in) :: lo_d,hi_d
  integer,intent(in) :: cnnct
  integer :: irng_r,jrng_r
  integer :: irng_d,jrng_d
  integer :: i,j
  integer :: i0,j0
  integer :: i1,j1
  integer :: incr1,incr2
  irng_r = abs(lo_r%d(1) - hi_r%d(1)) + 1
  jrng_r = abs(lo_r%d(2) - hi_r%d(2)) + 1
  irng_d = abs(hi_d%d(1) - lo_d%d(1)) + 1
  jrng_d = abs(hi_d%d(2) - lo_d%d(2)) + 1

  if(lo_d%d(1).lt.hi_d%d(1)) then
    incr1 = 1
  else
    incr1 =-1
  endif

  if(lo_d%d(2).lt.hi_d%d(2)) then
    incr2 = 1
  else
    incr2 =-1
  endif

  if(cnnct.eq.1) then
    if((irng_r.ne.irng_d).or.(jrng_r.ne.jrng_d)) then
      write(*,*) "Error: boundary_interface_vector()"
      write(*,*) "Incompatibal interface sizes:", &
                 irng_r,jrng_r,irng_d,jrng_d    
    endif
    do j = 0, jrng_r - 1
      do i = 0, irng_r - 1
        i0 = lo_r%d(1) + i
        j0 = lo_r%d(2) + j
        i1 = lo_d%d(1) + incr1*i
        j1 = lo_d%d(2) + incr2*j
        ! print*,"f_d",f_d%p(:,i1,j1)
        ! print*,"f_r",f_r%p(:,i1,j1)
        f_r%p(:,i0,j0) = f_d%p(:,i1,j1)
      enddo
    enddo
  elseif(cnnct.eq.2) then
    if((irng_r.ne.jrng_d).or.(jrng_r.ne.irng_d)) then
      write(*,*) "Error: boundary_interface_vector()"
      write(*,*) "Incompatibal interface sizes:", &
                 irng_r,jrng_r,jrng_d,irng_d    
    endif
    do j = 0, jrng_r - 1
      do i = 0, irng_r - 1
        i0 = lo_r%d(1) + i
        j0 = lo_r%d(2) + j
        i1 = lo_d%d(1) + incr1*j
        j1 = lo_d%d(2) + incr2*i
        f_r%p(:,i0,j0) = f_d%p(:,i1,j1)
      enddo
    enddo
  else
    print*,"Error: boundary_interface_vector()"
    print*,"no such cnnct:",cnnct
  endif
end subroutine boundary_interface_vector

subroutine boundary_interface_matrix(f_r,f_d,lo_r,hi_r,lo_d,hi_d,cnnct)
  implicit none
  type(rmatrixfield_t),intent(inout) :: f_r
  type(rmatrixfield_t),intent(inout) :: f_d
  type(fscalarshape_t),intent(in) :: lo_r,hi_r
  type(fscalarshape_t),intent(in) :: lo_d,hi_d
  integer,intent(in) :: cnnct
  integer :: irng_r,jrng_r
  integer :: irng_d,jrng_d
  integer :: i,j
  integer :: i0,j0
  integer :: i1,j1
  integer :: incr1,incr2

  irng_r = abs(lo_r%d(1) - hi_r%d(1)) + 1
  jrng_r = abs(lo_r%d(2) - hi_r%d(2)) + 1
  irng_d = abs(hi_d%d(1) - lo_d%d(1)) + 1
  jrng_d = abs(hi_d%d(2) - lo_d%d(2)) + 1

  if(lo_d%d(1).lt.hi_d%d(1)) then
    incr1 = 1
  else
    incr1 =-1
  endif

  if(lo_d%d(2).lt.hi_d%d(2)) then
    incr2 = 1
  else
    incr2 =-1
  endif

  if(cnnct.eq.1) then
    if((irng_r.ne.irng_d).or.(jrng_r.ne.jrng_d)) then
      write(*,*) "Error: boundary_interface_matrix()"
      write(*,*) "Incompatibal interface sizes:", &
                 irng_r,jrng_r,irng_d,jrng_d    
    endif
    do j = 0, jrng_r - 1
      do i = 0, irng_r - 1
        i0 = lo_r%d(1) + i
        j0 = lo_r%d(2) + j
        i1 = lo_d%d(1) + incr1*i
        j1 = lo_d%d(2) + incr2*j
        f_r%p(:,:,i0,j0) = f_d%p(:,:,i1,j1)
      enddo
    enddo
  elseif(cnnct.eq.2) then
    if((irng_r.ne.jrng_d).or.(jrng_r.ne.irng_d)) then
      write(*,*) "Error: boundary_interface_matrix()"
      write(*,*) "Incompatibal interface sizes:", &
                 irng_r,jrng_r,jrng_d,irng_d    
    endif
    do j = 0, jrng_r - 1
      do i = 0, irng_r - 1
        i0 = lo_r%d(1) + i
        j0 = lo_r%d(2) + j
        i1 = lo_d%d(1) + incr1*j
        j1 = lo_d%d(2) + incr2*i
        f_r%p(:,:,i0,j0) = f_d%p(:,:,i1,j1)
      enddo
    enddo
  else
    print*,"Error: boundary_interface_matrix()"
    print*,"no such cnnct:",cnnct
  endif
end subroutine boundary_interface_matrix

subroutine boundary_interface_tensor(f_r,f_d,lo_r,hi_r,lo_d,hi_d,cnnct)
  implicit none
  type(rtensorfield_t),intent(inout) :: f_r
  type(rtensorfield_t),intent(inout) :: f_d
  type(fscalarshape_t),intent(in) :: lo_r,hi_r
  type(fscalarshape_t),intent(in) :: lo_d,hi_d
  integer,intent(in) :: cnnct
  integer :: irng_r,jrng_r
  integer :: irng_d,jrng_d
  integer :: i,j
  integer :: i0,j0
  integer :: i1,j1
  integer :: incr1,incr2

  irng_r = abs(lo_r%d(1) - hi_r%d(1)) + 1
  jrng_r = abs(lo_r%d(2) - hi_r%d(2)) + 1
  irng_d = abs(hi_d%d(1) - lo_d%d(1)) + 1
  jrng_d = abs(hi_d%d(2) - lo_d%d(2)) + 1

  if(lo_d%d(1).lt.hi_d%d(1)) then
    incr1 = 1
  else
    incr1 =-1
  endif

  if(lo_d%d(2).lt.hi_d%d(2)) then
    incr2 = 1
  else
    incr2 =-1
  endif

  if(cnnct.eq.1) then
    if((irng_r.ne.irng_d).or.(jrng_r.ne.jrng_d)) then
      write(*,*) "Error: boundary_interface_tensor()"
      write(*,*) "Incompatibal interface sizes:", &
                 irng_r,jrng_r,irng_d,jrng_d    
    endif
    do j = 0, jrng_r - 1
      do i = 0, irng_r - 1
        i0 = lo_r%d(1) + i
        j0 = lo_r%d(2) + j
        i1 = lo_d%d(1) + incr1*i
        j1 = lo_d%d(2) + incr2*j
        f_r%p(:,:,:,i0,j0) = f_d%p(:,:,:,i1,j1)
      enddo
    enddo
  elseif(cnnct.eq.2) then
    if((irng_r.ne.jrng_d).or.(jrng_r.ne.irng_d)) then
      write(*,*) "Error: boundary_interface_tensor()"
      write(*,*) "Incompatibal interface sizes:", &
                 irng_r,jrng_r,jrng_d,irng_d    
    endif
    do j = 0, jrng_r - 1
      do i = 0, irng_r - 1
        i0 = lo_r%d(1) + i
        j0 = lo_r%d(2) + j
        i1 = lo_d%d(1) + incr1*j
        j1 = lo_d%d(2) + incr2*i
        f_r%p(:,:,:,i0,j0) = f_d%p(:,:,:,i1,j1)
      enddo
    enddo
  else
    print*,"Error: boundary_interface_tensor()"
    print*,"no such cnnct:",cnnct
  endif
end subroutine boundary_interface_tensor  

subroutine boundary_reflect_scalar(f,lo,hi,val,plane_id)
  implicit none
  type(fscalarshape_t),intent(in)    :: lo,hi
  integer,intent(in)                 :: plane_id
  double precision,intent(in)        :: val
  type(rscalarfield_t),intent(inout) :: f

  integer :: irng,jrng
  integer :: irng_d,jrng_d
  integer :: i,j
  integer :: i0,j0
  integer :: i1,j1  
  integer :: incr1,incr2
  type(fscalarshape_t) :: lo_d,hi_d

  lo_d = lo 
  hi_d = hi 
  select case(plane_id)
    case(1)
      lo_d%d(1) = lo%d(1) + abs(lo%d(1) - hi%d(1))*2 + 1
      hi_d%d(1) = hi%d(1) + 1
      incr1 =-1
      incr2 = 1
    case(2)
      lo_d%d(2) = lo%d(2) + abs(lo%d(2) - hi%d(2))*2 + 1
      hi_d%d(2) = hi%d(2) + 1 
      incr1 = 1
      incr2 =-1    
    case(3)
      lo_d%d(1) = lo%d(1) - 1
      hi_d%d(1) = hi%d(1) - abs(lo%d(1) - hi%d(1))*2 - 1 
      incr1 =-1
      incr2 = 1     
    case(4)
      lo_d%d(2) = lo%d(2) - 1 
      hi_d%d(2) = hi%d(2) - abs(lo%d(2) - hi%d(2))*2 - 1  
      incr1 = 1
      incr2 =-1
  end select

  irng   = abs(lo%d(1) - hi%d(1)) + 1
  jrng   = abs(lo%d(2) - hi%d(2)) + 1
  do j = 0, jrng - 1
    do i = 0, irng - 1
      i0 = lo%d(1) + i
      j0 = lo%d(2) + j
      i1 = lo_d%d(1) + incr1*i
      j1 = lo_d%d(2) + incr2*j
      f%p(i0,j0) = val*f%p(i1,j1)
    enddo
  enddo
end subroutine boundary_reflect_scalar

subroutine boundary_reflect_vector(f,lo,hi,val,plane_id)
  implicit none
  type(fscalarshape_t),intent(in)    :: lo,hi
  integer,intent(in)                 :: plane_id
  type(rvector_t),intent(in)        :: val
  type(rvectorfield_t),intent(inout) :: f

  integer :: irng,jrng
  integer :: irng_d,jrng_d
  integer :: i,j
  integer :: i0,j0
  integer :: i1,j1  
  integer :: incr1,incr2

  type(fscalarshape_t) :: lo_d,hi_d

  lo_d = lo 
  hi_d = hi 
  select case(plane_id)
    case(1)
      lo_d%d(1) = lo%d(1) + abs(lo%d(1) - hi%d(1))*2 + 1
      hi_d%d(1) = hi%d(1) + 1
      incr1 =-1
      incr2 = 1
    case(2)
      lo_d%d(2) = lo%d(2) + abs(lo%d(2) - hi%d(2))*2 + 1
      hi_d%d(2) = hi%d(2) + 1 
      incr1 = 1
      incr2 =-1    
    case(3)
      lo_d%d(1) = lo%d(1) - 1
      hi_d%d(1) = hi%d(1) - abs(lo%d(1) - hi%d(1))*2 - 1 
      incr1 =-1
      incr2 = 1     
    case(4)
      lo_d%d(2) = lo%d(2) - 1 
      hi_d%d(2) = hi%d(2) - abs(lo%d(2) - hi%d(2))*2 - 1  
      incr1 = 1
      incr2 =-1
  end select

  irng   = abs(lo%d(1) - hi%d(1)) + 1
  jrng   = abs(lo%d(2) - hi%d(2)) + 1
  do j = 0, jrng - 1
    do i = 0, irng - 1
      i0 = lo%d(1) + i
      j0 = lo%d(2) + j
      i1 = lo_d%d(1) + incr1*i
      j1 = lo_d%d(2) + incr2*j
      f%p(:,i0,j0) = val%p(:)*f%p(:,i1,j1)
    enddo
  enddo
end subroutine boundary_reflect_vector

subroutine boundary_reflect_matrix(f,lo,hi,val,plane_id)
  implicit none
  type(fscalarshape_t),intent(in)    :: lo,hi
  integer,intent(in)                 :: plane_id
  type(rmatrix_t),intent(in)        :: val
  type(rmatrixfield_t),intent(inout) :: f

  integer :: irng,jrng
  integer :: irng_d,jrng_d
  integer :: i,j
  integer :: i0,j0
  integer :: i1,j1  
  integer :: incr1,incr2

  type(fscalarshape_t) :: lo_d,hi_d

  lo_d = lo 
  hi_d = hi 
  select case(plane_id)
    case(1)
      lo_d%d(1) = lo%d(1) + abs(lo%d(1) - hi%d(1))*2 + 1
      hi_d%d(1) = hi%d(1) + 1
      incr1 =-1
      incr2 = 1
    case(2)
      lo_d%d(2) = lo%d(2) + abs(lo%d(2) - hi%d(2))*2 + 1
      hi_d%d(2) = hi%d(2) + 1 
      incr1 = 1
      incr2 =-1    
    case(3)
      lo_d%d(1) = lo%d(1) - 1
      hi_d%d(1) = hi%d(1) - abs(lo%d(1) - hi%d(1))*2 - 1 
      incr1 =-1
      incr2 = 1     
    case(4)
      lo_d%d(2) = lo%d(2) - 1 
      hi_d%d(2) = hi%d(2) - abs(lo%d(2) - hi%d(2))*2 - 1  
      incr1 = 1
      incr2 =-1
  end select

  irng   = abs(lo%d(1) - hi%d(1)) + 1
  jrng   = abs(lo%d(2) - hi%d(2)) + 1
  do j = 0, jrng - 1
    do i = 0, irng - 1
      i0 = lo%d(1) + i
      j0 = lo%d(2) + j
      i1 = lo_d%d(1) + incr1*i
      j1 = lo_d%d(2) + incr2*j
      f%p(:,:,i0,j0) = val%p(:,:)*f%p(:,:,i1,j1)
    enddo
  enddo
end subroutine boundary_reflect_matrix

subroutine boundary_reflect_tensor(f,lo,hi,val,plane_id)
  implicit none
  type(fscalarshape_t),intent(in)    :: lo,hi
  integer,intent(in)                 :: plane_id
  type(rtensor_t),intent(in)        :: val
  type(rtensorfield_t),intent(inout) :: f

  integer :: irng,jrng
  integer :: irng_d,jrng_d
  integer :: i,j
  integer :: i0,j0
  integer :: i1,j1  
  integer :: incr1,incr2
  type(fscalarshape_t) :: lo_d,hi_d

  lo_d = lo 
  hi_d = hi 
  select case(plane_id)
    case(1)
      lo_d%d(1) = lo%d(1) + abs(lo%d(1) - hi%d(1))*2 + 1
      hi_d%d(1) = hi%d(1) + 1
      incr1 =-1
      incr2 = 1
    case(2)
      lo_d%d(2) = lo%d(2) + abs(lo%d(2) - hi%d(2))*2 + 1
      hi_d%d(2) = hi%d(2) + 1 
      incr1 = 1
      incr2 =-1    
    case(3)
      lo_d%d(1) = lo%d(1) - 1
      hi_d%d(1) = hi%d(1) - abs(lo%d(1) - hi%d(1))*2 - 1 
      incr1 =-1
      incr2 = 1     
    case(4)
      lo_d%d(2) = lo%d(2) - 1 
      hi_d%d(2) = hi%d(2) - abs(lo%d(2) - hi%d(2))*2 - 1  
      incr1 = 1
      incr2 =-1
  end select

  irng   = abs(lo%d(1) - hi%d(1)) + 1
  jrng   = abs(lo%d(2) - hi%d(2)) + 1
  do j = 0, jrng - 1
    do i = 0, irng - 1
      i0 = lo%d(1) + i
      j0 = lo%d(2) + j
      i1 = lo_d%d(1) + incr1*i
      j1 = lo_d%d(2) + incr2*j
      f%p(:,:,:,i0,j0) = val%p(:,:,:)*f%p(:,:,:,i1,j1)
    enddo
  enddo
end subroutine boundary_reflect_tensor    

subroutine boundary_dirichlet_scalar(f,lo,hi,scalar)
  implicit none
  double precision,intent(in) :: scalar
  type(fscalarshape_t),intent(in) :: lo,hi
  type(rscalarfield_t),intent(inout) :: f

  integer :: irng,jrng
  integer :: i,j
  integer :: i0,j0
  integer :: incr1,incr2

  irng = abs(hi%d(1) - lo%d(1)) + 1
  jrng = abs(hi%d(2) - lo%d(2)) + 1
  
  if(lo%d(1).lt.hi%d(1)) then
    incr1 = 1
  else
    incr1 =-1
  endif

  if(lo%d(2).lt.hi%d(2)) then
    incr2 = 1
  else
    incr2 =-1
  endif

  do j = 0, jrng - 1
    do i = 0, irng - 1
      i0 = lo%d(1) + incr1*i
      j0 = lo%d(2) + incr2*j
      f%p(i0,j0) = scalar
    enddo
  enddo
end subroutine boundary_dirichlet_scalar  

subroutine boundary_dirichlet_vector(f,lo,hi,vector)
  implicit none
  type(rvector_t),intent(in) :: vector
  type(fscalarshape_t),intent(in) :: lo,hi
  type(rvectorfield_t),intent(inout) :: f

  integer :: irng,jrng
  integer :: i,j
  integer :: i0,j0
  integer :: incr1,incr2

  irng = abs(hi%d(1) - lo%d(1)) + 1
  jrng = abs(hi%d(2) - lo%d(2)) + 1
  
  if(lo%d(1).lt.hi%d(1)) then
    incr1 = 1
  else
    incr1 =-1
  endif

  if(lo%d(2).lt.hi%d(2)) then
    incr2 = 1
  else
    incr2 =-1
  endif

  do j = 0, jrng - 1
    do i = 0, irng - 1
      i0 = lo%d(1) + incr1*i
      j0 = lo%d(2) + incr2*j
      f%p(:,i0,j0) = vector%p(:)
    enddo
  enddo
end subroutine boundary_dirichlet_vector  

subroutine boundary_dirichlet_matrix(f,lo,hi,matrix)
  implicit none
  type(rmatrix_t),intent(in) :: matrix
  type(fscalarshape_t),intent(in) :: lo,hi
  type(rmatrixfield_t),intent(inout) :: f

  integer :: irng,jrng
  integer :: i,j
  integer :: i0,j0
  integer :: incr1,incr2

  irng = abs(hi%d(1) - lo%d(1)) + 1
  jrng = abs(hi%d(2) - lo%d(2)) + 1
  
  if(lo%d(1).lt.hi%d(1)) then
    incr1 = 1
  else
    incr1 =-1
  endif

  if(lo%d(2).lt.hi%d(2)) then
    incr2 = 1
  else
    incr2 =-1
  endif

  do j = 0, jrng - 1
    do i = 0, irng - 1
      i0 = lo%d(1) + incr1*i
      j0 = lo%d(2) + incr2*j
      f%p(:,:,i0,j0) = matrix%p(:,:)
    enddo
  enddo
end subroutine boundary_dirichlet_matrix  

subroutine boundary_dirichlet_tensor(f,lo,hi,tensor)
  implicit none
  type(rtensor_t),intent(in) :: tensor
  type(fscalarshape_t),intent(in) :: lo,hi
  type(rtensorfield_t),intent(inout) :: f

  integer :: irng,jrng
  integer :: i,j
  integer :: i0,j0
  integer :: incr1,incr2

  irng = abs(hi%d(1) - lo%d(1)) + 1
  jrng = abs(hi%d(2) - lo%d(2)) + 1
  
  if(lo%d(1).lt.hi%d(1)) then
    incr1 = 1
  else
    incr1 =-1
  endif

  if(lo%d(2).lt.hi%d(2)) then
    incr2 = 1
  else
    incr2 =-1
  endif

  do j = 0, jrng - 1
    do i = 0, irng - 1
      i0 = lo%d(1) + incr1*i
      j0 = lo%d(2) + incr2*j
      f%p(:,:,:,i0,j0) = tensor%p(:,:,:)
    enddo
  enddo
end subroutine boundary_dirichlet_tensor  

subroutine boundary_neumann_scalar(f,lo,hi,scalar,axis)
  implicit none
  integer,intent(in) :: axis
  double precision,intent(in) :: scalar
  type(fscalarshape_t),intent(in) :: lo,hi
  type(rscalarfield_t),intent(inout) :: f

  integer :: irng,jrng
  integer :: i,j
  integer :: i0,j0
  integer :: i1,j1  
  integer :: i2,j2  
  integer :: incr1,incr2
  integer,dimension(2) :: ishift,jshift,shift

  shift = (/-1,-2/)
  irng = abs(hi%d(1) - lo%d(1)) + 1
  jrng = abs(hi%d(2) - lo%d(2)) + 1
  
  if(lo%d(1).lt.hi%d(1)) then
    incr1 = 1
  else
    incr1 =-1
  endif

  if(lo%d(2).lt.hi%d(2)) then
    incr2 = 1
  else
    incr2 =-1
  endif

  select case(axis)
  case(axis_xi)
    ishift = incr1*shift
    jshift = 0
  case(axis_eta)
    ishift = 0
    jshift = incr2*shift
  end select 

  do j = 0, jrng - 1
    do i = 0, irng - 1
      i0 = lo%d(1) + incr1*i
      j0 = lo%d(2) + incr2*j
      i1 = i0 + ishift(1)
      j1 = j0 + jshift(1)
      i2 = i0 + ishift(2)
      j2 = j0 + jshift(2)
      f%p(i0,j0) = f%p(i1,j1)*4.0/3.0-f%p(i2,j2)/3.0 + scalar*2.d0/3.d0
    enddo
  enddo
end subroutine boundary_neumann_scalar  

subroutine boundary_neumann_vector(f,lo,hi,vector,axis)
  implicit none
  integer,intent(in) :: axis
  type(rvector_t),intent(in) :: vector
  type(fscalarshape_t),intent(in) :: lo,hi
  type(rvectorfield_t),intent(inout) :: f

  integer :: irng,jrng
  integer :: i,j,n
  integer :: i0,j0
  integer :: i1,j1  
  integer :: i2,j2  
  integer :: incr1,incr2
  integer,dimension(2) :: ishift,jshift,shift

  shift = (/-1,-2/)
  irng = abs(hi%d(1) - lo%d(1)) + 1
  jrng = abs(hi%d(2) - lo%d(2)) + 1
  
  if(lo%d(1).lt.hi%d(1)) then
    incr1 = 1
  else
    incr1 =-1
  endif

  if(lo%d(2).lt.hi%d(2)) then
    incr2 = 1
  else
    incr2 =-1
  endif

  select case(axis)
  case(axis_xi)
    ishift = incr1*shift
    jshift = 0
  case(axis_eta)
    ishift = 0
    jshift = incr2*shift
  end select 

  do j = 0, jrng - 1
    do i = 0, irng - 1
      i0 = lo%d(1) + incr1*i
      j0 = lo%d(2) + incr2*j
      i1 = i0 + ishift(1)
      j1 = j0 + jshift(1)
      i2 = i0 + ishift(2)
      j2 = j0 + jshift(2)
      f%p(:,i0,j0) = f%p(:,i1,j1)*4.0/3.0-f%p(:,i2,j2)/3.0 + vector%p(:)*2.d0/3.d0
    enddo
  enddo
end subroutine boundary_neumann_vector  

subroutine boundary_neumann_matrix(f,lo,hi,matrix,axis)
  implicit none
  integer,intent(in) :: axis
  type(rmatrix_t),intent(in) :: matrix
  type(fscalarshape_t),intent(in) :: lo,hi
  type(rmatrixfield_t),intent(inout) :: f

  integer :: irng,jrng
  integer :: i,j,n,m
  integer :: i0,j0
  integer :: i1,j1  
  integer :: i2,j2  
  integer :: incr1,incr2
  integer,dimension(2) :: ishift,jshift,shift

  shift = (/-1,-2/)
  irng = abs(hi%d(1) - lo%d(1)) + 1
  jrng = abs(hi%d(2) - lo%d(2)) + 1
  
  if(lo%d(1).lt.hi%d(1)) then
    incr1 = 1
  else
    incr1 =-1
  endif

  if(lo%d(2).lt.hi%d(2)) then
    incr2 = 1
  else
    incr2 =-1
  endif

  select case(axis)
  case(axis_xi)
    ishift = incr1*shift
    jshift = 0
  case(axis_eta)
    ishift = 0
    jshift = incr2*shift
  end select 

  do j = 0, jrng - 1
    do i = 0, irng - 1
      i0 = lo%d(1) + incr1*i
      j0 = lo%d(2) + incr2*j
      i1 = i0 + ishift(1)
      j1 = j0 + jshift(1)
      i2 = i0 + ishift(2)
      j2 = j0 + jshift(2)
      f%p(:,:,i0,j0) = f%p(:,:,i1,j1)*4.0/3.0-f%p(:,:,i2,j2)/3.0 + matrix%p(:,:)*2.d0/3.d0
    enddo
  enddo
end subroutine boundary_neumann_matrix  

subroutine boundary_neumann_tensor(f,lo,hi,tensor,axis)
  implicit none
  integer,intent(in) :: axis
  type(rtensor_t),intent(in) :: tensor
  type(fscalarshape_t),intent(in) :: lo,hi
  type(rtensorfield_t),intent(inout) :: f

  integer :: irng,jrng
  integer :: i,j,n,m,l
  integer :: i0,j0
  integer :: i1,j1  
  integer :: i2,j2  
  integer :: incr1,incr2
  integer,dimension(2) :: ishift,jshift,shift

  shift = (/-1,-2/)
  irng = abs(hi%d(1) - lo%d(1)) + 1
  jrng = abs(hi%d(2) - lo%d(2)) + 1
  
  if(lo%d(1).lt.hi%d(1)) then
    incr1 = 1
  else
    incr1 =-1
  endif

  if(lo%d(2).lt.hi%d(2)) then
    incr2 = 1
  else
    incr2 =-1
  endif

  select case(axis)
  case(axis_xi)
    ishift = incr1*shift
    jshift = 0
  case(axis_eta)
    ishift = 0
    jshift = incr2*shift
  end select 

  do j = 0, jrng - 1
    do i = 0, irng - 1
      i0 = lo%d(1) + incr1*i
      j0 = lo%d(2) + incr2*j
      i1 = i0 + ishift(1)
      j1 = j0 + jshift(1)
      i2 = i0 + ishift(2)
      j2 = j0 + jshift(2)
      f%p(:,:,:,i0,j0) = f%p(:,:,:,i1,j1)*4.0/3.0 &
                        -f%p(:,:,:,i2,j2)/3.0     &
                        +tensor%p(:,:,:)*2.d0/3.d0 
    enddo
  enddo
end subroutine boundary_neumann_tensor  
end module boundary