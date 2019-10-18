submodule (reconstruction) weno3
implicit none
! third order weno coefficients
double precision,parameter :: a200 =-0.5d0
double precision,parameter :: a201 = 1.5d0
double precision,parameter :: a210 = 0.5d0
double precision,parameter :: a211 = 0.5d0
double precision,parameter :: c20  = 1.d0/3.d0
double precision,parameter :: c21  = 2.d0/3.d0

double precision,parameter :: ss   = 1.0e-20

contains

module subroutine reconst_component_scalar_wenojs3(h,f,lo,hi,padding,bias,axis)
    implicit none
    type(rscalarfield_t),intent(inout) :: h
    type(rscalarfield_t),target,intent(in) :: f
    type(fscalarshape_t),intent(in) :: lo,hi,padding
    integer,intent(in) :: bias,axis
  
    double precision,dimension(:,:),pointer :: f0
    double precision,dimension(:,:),pointer :: f1
    double precision,dimension(:,:),pointer :: f2
    integer :: filo,fihi, &
               fjlo,fjhi
  
    double precision :: is0,is1
    double precision :: aa0,aa1,aa_sum
    double precision :: w0,w1
    double precision :: q20,q21
    integer,dimension(3) :: shift
    integer :: i,j
    if(bias.eq.bias_upwind) then
      shift = (/-1,0,1/)
    else 
      shift = (/2,1,0/)
    endif 
  
    filo = lo%d(1) + padding%d(1)
    fihi = hi%d(1) + padding%d(1)
    fjlo = lo%d(2) + padding%d(2)
    fjhi = hi%d(2) + padding%d(2)
    select case(axis)
    case(axis_xi)
      f0 => f%p(filo+shift(1):fihi+shift(1),fjlo:fjhi)
      f1 => f%p(filo+shift(2):fihi+shift(2),fjlo:fjhi)
      f2 => f%p(filo+shift(3):fihi+shift(3),fjlo:fjhi)
    case(axis_eta)
      f0 => f%p(filo:fihi,fjlo+shift(1):fjhi+shift(1))
      f1 => f%p(filo:fihi,fjlo+shift(2):fjhi+shift(2))
      f2 => f%p(filo:fihi,fjlo+shift(3):fjhi+shift(3))
    case default
      print*,"Error: reconst_component_scalar_wenojs3()"
      print*,"No such ireconst: ",axis
      stop
    end select 
  
    do j=lo%d(2), hi%d(2)
      do i=lo%d(1), hi%d(1)      
        is0 = (f0(i,j)-f1(i,j))**2
        is1 = (f1(i,j)-f2(i,j))**2
  
        aa0    = c20/(ss+is0)**2
        aa1    = c21/(ss+is1)**2
        aa_sum = 1.d0/(aa0+aa1)
      
        w0 = aa0*aa_sum
        w1 = aa1*aa_sum
        
        q20 = a200*f0(i,j) + a201*f1(i,j)
        q21 = a210*f1(i,j) + a211*f2(i,j)
  
        h%p(i,j) = w0*q20 + w1*q21
      enddo
    enddo
  end subroutine reconst_component_scalar_wenojs3
  
  module subroutine reconst_component_vector_wenojs3(h,f,lo,hi,padding,bias,axis)
    implicit none
    type(rvectorfield_t),intent(inout) :: h
    type(rvectorfield_t),target,intent(in) :: f
    type(fvectorshape_t),intent(in) :: lo,hi,padding
    integer,intent(in) :: bias,axis
  
    double precision,dimension(:,:,:),pointer :: f0
    double precision,dimension(:,:,:),pointer :: f1
    double precision,dimension(:,:,:),pointer :: f2
    integer :: fnlo,fnhi, &
               filo,fihi, &
               fjlo,fjhi
  
    double precision,dimension(lo%d(1):hi%d(1)) :: is0,is1
    double precision,dimension(lo%d(1):hi%d(1)) :: aa0,aa1,aa_sum
    double precision,dimension(lo%d(1):hi%d(1)) :: w0,w1
    double precision,dimension(lo%d(1):hi%d(1)) :: q20,q21
    integer :: i,j
    integer :: nlo,nhi
    integer,dimension(3) :: shift
    if(bias.eq.bias_upwind) then
      shift = (/-1,0,1/)
    else 
      shift = (/2,1,0/)
    endif 
    fnlo = lo%d(1) + padding%d(1)
    fnhi = hi%d(1) + padding%d(1)
    filo = lo%d(2) + padding%d(2)
    fihi = hi%d(2) + padding%d(2)
    fjlo = lo%d(3) + padding%d(3)
    fjhi = hi%d(3) + padding%d(3)
    select case(axis)
    case(axis_xi)
      f0 => f%p(fnlo:fnhi,filo+shift(1):fihi+shift(1),fjlo:fjhi)
      f1 => f%p(fnlo:fnhi,filo+shift(2):fihi+shift(2),fjlo:fjhi)
      f2 => f%p(fnlo:fnhi,filo+shift(3):fihi+shift(3),fjlo:fjhi)
    case(axis_eta)
      f0 => f%p(fnlo:fnhi,filo:fihi,fjlo+shift(1):fjhi+shift(1))
      f1 => f%p(fnlo:fnhi,filo:fihi,fjlo+shift(2):fjhi+shift(2))
      f2 => f%p(fnlo:fnhi,filo:fihi,fjlo+shift(3):fjhi+shift(3))
    case default
      print*,"Error: reconst_component_vector_wenojs3()"
      print*,"No such ireconst: ",axis
      stop
    end select 
  
    do j=lo%d(3), hi%d(3)
      do i=lo%d(2), hi%d(2)
        is0(:) = (f0(:,i,j)-f1(:,i,j))**2
        is1(:) = (f1(:,i,j)-f2(:,i,j))**2
  
        aa0(:)    = c20/(ss+is0(:))**2
        aa1(:)    = c21/(ss+is1(:))**2
        aa_sum(:) = 1.d0/(aa0(:)+aa1(:))
      
        w0(:) = aa0(:)*aa_sum(:)
        w1(:) = aa1(:)*aa_sum(:)
        
        q20(:) = a200*f0(:,i,j) + a201*f1(:,i,j)
        q21(:) = a210*f1(:,i,j) + a211*f2(:,i,j)
  
        h%p(lo%d(1):hi%d(1),i,j) = w0*q20 + w1*q21
      enddo
    enddo
  end subroutine reconst_component_vector_wenojs3
  
  module subroutine reconst_component_matrix_wenojs3(h,f,lo,hi,padding,bias,axis)
    implicit none
    type(rmatrixfield_t),intent(inout) :: h
    type(rmatrixfield_t),target,intent(in) :: f
    type(fmatrixshape_t),intent(in) :: lo,hi,padding
    integer,intent(in) :: bias,axis
  
    double precision,dimension(:,:,:,:),pointer :: f0
    double precision,dimension(:,:,:,:),pointer :: f1
    double precision,dimension(:,:,:,:),pointer :: f2
    integer :: fmlo,fmhi, &
               fnlo,fnhi, &
               filo,fihi, &
               fjlo,fjhi
  
    double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2)) :: is0,is1
    double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2)) :: aa0,aa1,aa_sum
    double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2)) :: w0,w1
    double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2)) :: q20,q21
    integer :: i,j
    integer,dimension(3) :: shift
    if(bias.eq.bias_upwind) then
      shift = (/-1,0,1/)
    else 
      shift = (/2,1,0/)
    endif 
  
    fmlo = lo%d(1) + padding%d(1)
    fmhi = hi%d(1) + padding%d(1)
    fnlo = lo%d(2) + padding%d(2)
    fnhi = hi%d(2) + padding%d(2)
    filo = lo%d(3) + padding%d(3)
    fihi = hi%d(3) + padding%d(3)
    fjlo = lo%d(4) + padding%d(4)
    fjhi = hi%d(4) + padding%d(4)
    select case(axis)
    case(axis_xi)
      f0 => f%p(fmlo:fmhi,fnlo:fnhi,filo+shift(1):fihi+shift(1),fjlo:fjhi)
      f1 => f%p(fmlo:fmhi,fnlo:fnhi,filo+shift(2):fihi+shift(2),fjlo:fjhi)
      f2 => f%p(fmlo:fmhi,fnlo:fnhi,filo+shift(3):fihi+shift(3),fjlo:fjhi)
    case(axis_eta)
      f0 => f%p(fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(1):fjhi+shift(1))
      f1 => f%p(fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(2):fjhi+shift(2))
      f2 => f%p(fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(3):fjhi+shift(3))
    case default
      print*,"Error: reconst_component_matrix_wenojs3()"
      print*,"No such ireconst: ",axis
      stop
    end select 
  
    do j = lo%d(4), hi%d(4)
      do i = lo%d(3), hi%d(3)
        is0(:,:) = (f0(:,:,i,j)-f1(:,:,i,j))**2
        is1(:,:) = (f1(:,:,i,j)-f2(:,:,i,j))**2
  
        aa0(:,:)    = c20/(ss+is0(:,:))**2
        aa1(:,:)    = c21/(ss+is1(:,:))**2
        aa_sum(:,:) = 1.d0/(aa0(:,:)+aa1(:,:))
      
        w0(:,:) = aa0(:,:)*aa_sum(:,:)
        w1(:,:) = aa1(:,:)*aa_sum(:,:)
        
        q20(:,:) = a200*f0(:,:,i,j) + a201*f1(:,:,i,j)
        q21(:,:) = a210*f1(:,:,i,j) + a211*f2(:,:,i,j)
  
        h%p(lo%d(1):hi%d(1),lo%d(2):hi%d(2),i,j) = w0*q20 + w1*q21
      enddo
    enddo
  end subroutine reconst_component_matrix_wenojs3
  
  module subroutine reconst_component_tensor_wenojs3(h,f,lo,hi,padding,bias,axis)
    implicit none
    type(rtensorfield_t),intent(inout) :: h
    type(rtensorfield_t),target,intent(in) :: f
    type(ftensorshape_t),intent(in) :: lo,hi,padding
    integer,intent(in) :: bias,axis
  
    double precision,dimension(:,:,:,:,:),pointer :: f0
    double precision,dimension(:,:,:,:,:),pointer :: f1
    double precision,dimension(:,:,:,:,:),pointer :: f2
    integer :: fllo,flhi, &
               fmlo,fmhi, &
               fnlo,fnhi, &
               filo,fihi, &
               fjlo,fjhi
  
    double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2),lo%d(3):hi%d(3)) :: is0,is1
    double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2),lo%d(3):hi%d(3)) :: aa0,aa1,aa_sum
    double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2),lo%d(3):hi%d(3)) :: w0,w1
    double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2),lo%d(3):hi%d(3)) :: q20,q21
    integer :: i,j
    integer,dimension(3) :: shift
    if(bias.eq.bias_upwind) then
      shift = (/-1,0,1/)
    else 
      shift = (/2,1,0/)
    endif 
    fllo = lo%d(1) + padding%d(1)
    flhi = hi%d(1) + padding%d(1)
    fmlo = lo%d(2) + padding%d(2)
    fmhi = hi%d(2) + padding%d(2)
    fnlo = lo%d(3) + padding%d(3)
    fnhi = hi%d(3) + padding%d(3)
    filo = lo%d(4) + padding%d(4)
    fihi = hi%d(4) + padding%d(4)
    fjlo = lo%d(5) + padding%d(5)
    fjhi = hi%d(5) + padding%d(5)
    select case(axis)
    case(axis_xi)
      f0 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo+shift(1):fihi+shift(1),fjlo:fjhi)
      f1 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo+shift(2):fihi+shift(2),fjlo:fjhi)
      f2 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo+shift(3):fihi+shift(3),fjlo:fjhi)
    case(axis_eta)
      f0 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(1):fjhi+shift(1))
      f1 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(2):fjhi+shift(2))
      f2 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(3):fjhi+shift(3))
    case default
      print*,"Error: reconst_component_tensor_wenojs3()"
      print*,"No such ireconst: ",axis
      stop
    end select 
    do j = lo%d(5), hi%d(5)
      do i = lo%d(4), hi%d(4)            
        is0(:,:,:) = (f0(:,:,:,i,j)-f1(:,:,:,i,j))**2
        is1(:,:,:) = (f1(:,:,:,i,j)-f2(:,:,:,i,j))**2
  
        aa0(:,:,:)    = c20/(ss+is0(:,:,:))**2
        aa1(:,:,:)    = c21/(ss+is1(:,:,:))**2
        aa_sum(:,:,:) = 1.d0/(aa0(:,:,:)+aa1(:,:,:))
      
        w0(:,:,:) = aa0(:,:,:)*aa_sum(:,:,:)
        w1(:,:,:) = aa1(:,:,:)*aa_sum(:,:,:)
        
        q20(:,:,:) = a200*f0(:,:,:,i,j) + a201*f1(:,:,:,i,j)
        q21(:,:,:) = a210*f1(:,:,:,i,j) + a211*f2(:,:,:,i,j)
  
        h%p(lo%d(1):hi%d(1),lo%d(2):hi%d(2),lo%d(3):hi%d(3),i,j) = w0*q20 + w1*q21
      enddo
    enddo
  end subroutine reconst_component_tensor_wenojs3
end submodule