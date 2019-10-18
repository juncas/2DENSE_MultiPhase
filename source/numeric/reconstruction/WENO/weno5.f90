submodule (reconstruction) weno5
implicit none
! fifth order weno coefficients
double precision,parameter :: a300 = 1.d0/3.d0
double precision,parameter :: a301 =-7.d0/6.d0
double precision,parameter :: a302 = 11.d0/6.d0
double precision,parameter :: a310 =-1.d0/6.d0
double precision,parameter :: a311 = 5.d0/6.d0
double precision,parameter :: a312 = 1.d0/3.d0
double precision,parameter :: a320 = 1.d0/3.d0
double precision,parameter :: a321 = 5.d0/6.d0
double precision,parameter :: a322 =-1.d0/6.d0
double precision,parameter :: c30  = 0.1d0
double precision,parameter :: c31  = 0.6d0
double precision,parameter :: c32  = 0.3d0

double precision,parameter :: ss   = 1.0e-20

contains
module subroutine reconst_component_scalar_wenojs5(h,f,lo,hi,padding,bias,axis)
  implicit none
  type(rscalarfield_t),intent(inout) :: h
  type(rscalarfield_t),target,intent(in) :: f
  type(fscalarshape_t),intent(in) :: lo,hi,padding
  integer,intent(in) :: bias,axis

  double precision,dimension(:,:),pointer :: f0
  double precision,dimension(:,:),pointer :: f1
  double precision,dimension(:,:),pointer :: f2
  double precision,dimension(:,:),pointer :: f3
  double precision,dimension(:,:),pointer :: f4
  integer :: filo,fihi, &
             fjlo,fjhi

  double precision :: is0,is1,is2
  double precision :: aa0,aa1,aa2,aa_sum
  double precision :: w0,w1,w2
  double precision :: q30,q31,q32
  integer :: i,j
  integer,dimension(5) :: shift

  if(bias.eq.bias_upwind) then
    shift = (/-2,-1,0,1,2/)
  else 
    shift = (/3,2,1,0,-1/)
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
    f3 => f%p(filo+shift(4):fihi+shift(4),fjlo:fjhi)
    f4 => f%p(filo+shift(5):fihi+shift(5),fjlo:fjhi)
  case(axis_eta)
    f0 => f%p(filo:fihi,fjlo+shift(1):fjhi+shift(1))
    f1 => f%p(filo:fihi,fjlo+shift(2):fjhi+shift(2))
    f2 => f%p(filo:fihi,fjlo+shift(3):fjhi+shift(3))
    f3 => f%p(filo:fihi,fjlo+shift(4):fjhi+shift(4))
    f4 => f%p(filo:fihi,fjlo+shift(5):fjhi+shift(5))
  case default
    print*,"Error: reconst_component_scalar_wenojs5()"
    print*,"No such ireconst: ",axis
    stop
  end select 

  do j=lo%d(2), hi%d(2)
    do i=lo%d(1), hi%d(1)
      is0 = (f0(i,j) - 2.d0*f1(i,j) + f2(i,j))**2*13.d0/12.d0 &
           +(f1(i,j) - 4.d0*f2(i,j) + 3.d0*f3(i,j))**2/4.d0

      is1 = (f1(i,j) - 2.d0*f2(i,j)+f3(i,j))**2*13.d0/12.d0 &
           +(f1(i,j) - f3(i,j))**2/4.d0

      is2 = (f2(i,j) - 2.d0*f3(i,j) + f4(i,j))**2*13.d0/12.d0 &
           +(3.d0*f2(i,j) - 4.d0*f3(i,j) + f4(i,j))**2/4.d0

      aa0 = c30/(ss+is0)**2
      aa1 = c31/(ss+is1)**2
      aa2 = c32/(ss+is2)**2
      aa_sum = 1.d0/(aa0+aa1+aa2)
    
      w0 = aa0*aa_sum
      w1 = aa1*aa_sum
      w2 = aa2*aa_sum 
      
      q30 = a300*f0(i,j) + a301*f1(i,j) + a302*f2(i,j)
      q31 = a310*f1(i,j) + a311*f3(i,j) + a312*f3(i,j)
      q32 = a320*f2(i,j) + a321*f3(i,j) + a322*f4(i,j)

      h%p(i,j) = w0*q30 + w1*q31 + w2*q32
    enddo
  enddo
end subroutine reconst_component_scalar_wenojs5

module subroutine reconst_component_vector_wenojs5(h,f,lo,hi,padding,bias,axis)
  implicit none
  type(rvectorfield_t),intent(inout) :: h
  type(rvectorfield_t),target,intent(in) :: f
  type(fvectorshape_t),intent(in) :: lo,hi,padding
  integer,intent(in) :: bias,axis

  double precision,dimension(:,:,:),pointer :: f0
  double precision,dimension(:,:,:),pointer :: f1
  double precision,dimension(:,:,:),pointer :: f2
  double precision,dimension(:,:,:),pointer :: f3
  double precision,dimension(:,:,:),pointer :: f4
  integer :: fnlo,fnhi, &
             filo,fihi, &
             fjlo,fjhi

  double precision,dimension(lo%d(1):hi%d(1)) :: is0,is1,is2
  double precision,dimension(lo%d(1):hi%d(1)):: aa0,aa1,aa2,aa_sum
  double precision,dimension(lo%d(1):hi%d(1)) :: w0,w1,w2
  double precision,dimension(lo%d(1):hi%d(1)) :: q30,q31,q32
  integer :: i,j
  integer,dimension(5) :: shift
  if(bias.eq.bias_upwind) then
    shift = (/-2,-1,0,1,2/)
  else 
    shift = (/3,2,1,0,-1/)
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
    f3 => f%p(fnlo:fnhi,filo+shift(4):fihi+shift(4),fjlo:fjhi)
    f4 => f%p(fnlo:fnhi,filo+shift(5):fihi+shift(5),fjlo:fjhi)
  case(axis_eta)
    f0 => f%p(fnlo:fnhi,filo:fihi,fjlo+shift(1):fjhi+shift(1))
    f1 => f%p(fnlo:fnhi,filo:fihi,fjlo+shift(2):fjhi+shift(2))
    f2 => f%p(fnlo:fnhi,filo:fihi,fjlo+shift(3):fjhi+shift(3))
    f3 => f%p(fnlo:fnhi,filo:fihi,fjlo+shift(4):fjhi+shift(4))
    f4 => f%p(fnlo:fnhi,filo:fihi,fjlo+shift(5):fjhi+shift(5))
  case default
    print*,"Error: reconst_component_vector_wenojs5()"
    print*,"No such ireconst: ",axis
    stop
  end select 

  do j=lo%d(3), hi%d(3)
    do i=lo%d(2), hi%d(2)
      is0(:) = (f0(:,i,j)-2.d0*f1(:,i,j)+f2(:,i,j))**2*13.d0/12.d0 &
              +(f0(:,i,j)-4.d0*f1(:,i,j)+3.d0*f2(:,i,j))**2/4.d0

      is1(:) = (f1(:,i,j)-2.d0*f2(:,i,j)+f3(:,i,j))**2*13.d0/12.d0 &
              +(f1(:,i,j)-f3(:,i,j))**2/4.d0

      is2(:) = (f2(:,i,j)-2.d0*f3(:,i,j)+f4(:,i,j))**2*13.d0/12.d0 &
              +(3.d0*f2(:,i,j)-4.d0*f3(:,i,j)+f4(:,i,j))**2/4.d0

      aa0(:)    = c30/(ss+is0(:))**2
      aa1(:)    = c31/(ss+is1(:))**2
      aa2(:)    = c32/(ss+is2(:))**2
      aa_sum(:) = 1.d0/(aa0(:)+aa1(:)+aa2(:))
    
      w0(:) = aa0(:)*aa_sum(:)
      w1(:) = aa1(:)*aa_sum(:)
      w2(:) = aa2(:)*aa_sum(:)
      
      q30(:) = a300*f0(:,i,j) + a301*f1(:,i,j) + a302*f2(:,i,j)
      q31(:) = a310*f1(:,i,j) + a311*f2(:,i,j) + a312*f3(:,i,j)
      q32(:) = a320*f2(:,i,j) + a321*f3(:,i,j) + a322*f4(:,i,j)

      h%p(lo%d(1):hi%d(1),i,j) = w0(:)*q30(:) + w1(:)*q31(:) + w2(:)*q32(:)
    enddo
  enddo
end subroutine reconst_component_vector_wenojs5

module subroutine reconst_component_matrix_wenojs5(h,f,lo,hi,padding,bias,axis)
  implicit none
  type(rmatrixfield_t),intent(inout) :: h
  type(rmatrixfield_t),target,intent(in) :: f
  type(fmatrixshape_t),intent(in) :: lo,hi,padding
  integer,intent(in) :: bias,axis

  double precision,dimension(:,:,:,:),pointer :: f0
  double precision,dimension(:,:,:,:),pointer :: f1
  double precision,dimension(:,:,:,:),pointer :: f2
  double precision,dimension(:,:,:,:),pointer :: f3
  double precision,dimension(:,:,:,:),pointer :: f4
  integer :: fmlo,fmhi, &
             fnlo,fnhi, &
             filo,fihi, &
             fjlo,fjhi

  double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2)) :: is0,is1,is2
  double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2)) :: aa0,aa1,aa2,aa_sum
  double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2)) :: w0,w1,w2
  double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2)) :: q30,q31,q32
  integer :: i,j
  integer,dimension(5) :: shift
  if(bias.eq.bias_upwind) then
    shift = (/-2,-1,0,1,2/)
  else 
    shift = (/3,2,1,0,-1/)
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
    f3 => f%p(fmlo:fmhi,fnlo:fnhi,filo+shift(4):fihi+shift(4),fjlo:fjhi)
    f4 => f%p(fmlo:fmhi,fnlo:fnhi,filo+shift(5):fihi+shift(5),fjlo:fjhi)
  case(axis_eta)
    f0 => f%p(fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(1):fjhi+shift(1))
    f1 => f%p(fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(2):fjhi+shift(2))
    f2 => f%p(fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(3):fjhi+shift(3))
    f3 => f%p(fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(4):fjhi+shift(4))
    f4 => f%p(fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(5):fjhi+shift(5))
  case default
    print*,"Error: reconst_component_matrix_wenojs5()"
    print*,"No such ireconst: ",axis
    stop
  end select 
  do j = lo%d(4), hi%d(4)
    do i = lo%d(3), hi%d(3)
      is0(:,:) = (f0(:,:,i,j)-2.d0*f2(:,:,i,j)+f3(:,:,i,j))**2*13.d0/12.d0 &
                +(f0(:,:,i,j)-4.d0*f1(:,:,i,j)+3.d0*f2(:,:,i,j))**2/4.d0

      is1(:,:) = (f1(:,:,i,j)-2.d0*f2(:,:,i,j)+f3(:,:,i,j))**2*13.d0/12.d0 &
                +(f1(:,:,i,j)-f3(:,:,i,j))**2/4.d0

      is2(:,:) = (f2(:,:,i,j)-2.d0*f3(:,:,i,j)+f4(:,:,i,j))**2*13.d0/12.d0 &
                +(3.d0*f2(:,:,i,j)-4.d0*f3(:,:,i,j)+f4(:,:,i,j))**2/4.d0

      aa0(:,:)    = c30/(ss+is0(:,:))**2
      aa1(:,:)    = c31/(ss+is1(:,:))**2
      aa2(:,:)    = c32/(ss+is2(:,:))**2
      aa_sum(:,:) = 1.d0/(aa0(:,:)+aa1(:,:)+aa2(:,:))
    
      w0(:,:) = aa0(:,:)*aa_sum(:,:)
      w1(:,:) = aa1(:,:)*aa_sum(:,:)
      w2(:,:) = aa2(:,:)*aa_sum(:,:)
      
      q30(:,:) = a300*f0(:,:,i,j) + a301*f1(:,:,i,j) + a302*f2(:,:,i,j)
      q31(:,:) = a310*f1(:,:,i,j) + a311*f2(:,:,i,j) + a312*f3(:,:,i,j)
      q32(:,:) = a320*f2(:,:,i,j) + a321*f3(:,:,i,j) + a322*f4(:,:,i,j)

      h%p(lo%d(1):hi%d(1),lo%d(2):hi%d(2),i,j) = w0(:,:)*q30(:,:) &
                                                +w1(:,:)*q31(:,:) &
                                                +w2(:,:)*q32(:,:)
    enddo
  enddo
end subroutine reconst_component_matrix_wenojs5

module subroutine reconst_component_tensor_wenojs5(h,f,lo,hi,padding,bias,axis)
  implicit none
  type(rtensorfield_t),intent(inout) :: h
  type(rtensorfield_t),target,intent(in) :: f
  type(ftensorshape_t),intent(in) :: lo,hi,padding
  integer,intent(in) :: bias,axis
  
  double precision,dimension(:,:,:,:,:),pointer :: f0
  double precision,dimension(:,:,:,:,:),pointer :: f1
  double precision,dimension(:,:,:,:,:),pointer :: f2
  double precision,dimension(:,:,:,:,:),pointer :: f3
  double precision,dimension(:,:,:,:,:),pointer :: f4
  integer :: fllo,flhi, &
             fmlo,fmhi, &
             fnlo,fnhi, &
             filo,fihi, &
             fjlo,fjhi

  double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2),lo%d(3):hi%d(3)) :: is0,is1,is2
  double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2),lo%d(3):hi%d(3)) :: aa0,aa1,aa2,aa_sum
  double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2),lo%d(3):hi%d(3)) :: w0,w1,w2
  double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2),lo%d(3):hi%d(3)) :: q30,q31,q32
  integer :: i,j
  integer,dimension(5) :: shift
  if(bias.eq.bias_upwind) then
    shift = (/-2,-1,0,1,2/)
  else 
    shift = (/3,2,1,0,-1/)
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
    f3 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo+shift(4):fihi+shift(4),fjlo:fjhi)
    f4 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo+shift(5):fihi+shift(5),fjlo:fjhi)
  case(axis_eta)
    f0 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(1):fjhi+shift(1))
    f1 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(2):fjhi+shift(2))
    f2 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(3):fjhi+shift(3))
    f3 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(4):fjhi+shift(4))
    f4 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(5):fjhi+shift(5))
  case default
    print*,"Error: reconst_component_tensor_wenojs5()"
    print*,"No such ireconst: ",axis
    stop
  end select 
  do j = lo%d(5), hi%d(5)
    do i = lo%d(4), hi%d(4)
      is0(:,:,:) = (f0(:,:,:,i,j)-2.d0*f1(:,:,:,i,j)+f2(:,:,:,i,j))**2*13.d0/12.d0 &
                  +(f0(:,:,:,i,j)-4.d0*f1(:,:,:,i,j)+3.d0*f2(:,:,:,i,j))**2/4.d0

      is1(:,:,:) = (f1(:,:,:,i,j)-2.d0*f2(:,:,:,i,j)+f3(:,:,:,i,j))**2*13.d0/12.d0 &
                  +(f1(:,:,:,i,j)-f3(:,:,:,i,j))**2/4.d0

      is2(:,:,:) = (f2(:,:,:,i,j)-2.d0*f3(:,:,:,i,j)+f4(:,:,:,i,j))**2*13.d0/12.d0 &
                  +(3.d0*f2(:,:,:,i,j)-4.d0*f3(:,:,:,i,j)+f4(:,:,:,i,j))**2/4.d0

      aa0(:,:,:)    = c30/(ss+is0(:,:,:))**2
      aa1(:,:,:)    = c31/(ss+is1(:,:,:))**2
      aa2(:,:,:)    = c32/(ss+is2(:,:,:))**2
      aa_sum(:,:,:) = 1.d0/(aa0(:,:,:)+aa1(:,:,:)+aa2(:,:,:))
    
      w0(:,:,:) = aa0(:,:,:)*aa_sum(:,:,:)
      w1(:,:,:) = aa1(:,:,:)*aa_sum(:,:,:)
      w2(:,:,:) = aa2(:,:,:)*aa_sum(:,:,:)
      
      q30(:,:,:) = a300*f0(:,:,:,i,j) + a301*f1(:,:,:,i,j) + a302*f2(:,:,:,i,j)
      q31(:,:,:) = a310*f1(:,:,:,i,j) + a311*f2(:,:,:,i,j) + a312*f3(:,:,:,i,j)
      q32(:,:,:) = a320*f2(:,:,:,i,j) + a321*f3(:,:,:,i,j) + a322*f4(:,:,:,i,j)

      h%p(lo%d(1):hi%d(1),lo%d(2):hi%d(2),lo%d(3):hi%d(3),i,j) = w0(:,:,:)*q30(:,:,:) &
                                                                +w1(:,:,:)*q31(:,:,:) &
                                                                +w2(:,:,:)*q32(:,:,:)
    enddo
  enddo
end subroutine reconst_component_tensor_wenojs5

module subroutine reconst_component_scalar_wenoz5(h,f,lo,hi,padding,bias,axis)
  implicit none
  type(rscalarfield_t),intent(inout) :: h
  type(rscalarfield_t),target,intent(in) :: f
  type(fscalarshape_t),intent(in) :: lo,hi,padding
  integer,intent(in) :: bias,axis

  double precision,dimension(:,:),pointer :: f0
  double precision,dimension(:,:),pointer :: f1
  double precision,dimension(:,:),pointer :: f2
  double precision,dimension(:,:),pointer :: f3
  double precision,dimension(:,:),pointer :: f4
  integer :: filo,fihi, &
             fjlo,fjhi

  double precision :: is0,is1,is2
  double precision :: tau5
  double precision :: aa0,aa1,aa2,aa_sum
  double precision :: w0,w1,w2
  double precision :: q30,q31,q32
  integer :: i,j
  integer,dimension(5) :: shift

  if(bias.eq.bias_upwind) then
    shift = (/-2,-1,0,1,2/)
  else 
    shift = (/3,2,1,0,-1/)
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
    f3 => f%p(filo+shift(4):fihi+shift(4),fjlo:fjhi)
    f4 => f%p(filo+shift(5):fihi+shift(5),fjlo:fjhi)
  case(axis_eta)
    f0 => f%p(filo:fihi,fjlo+shift(1):fjhi+shift(1))
    f1 => f%p(filo:fihi,fjlo+shift(2):fjhi+shift(2))
    f2 => f%p(filo:fihi,fjlo+shift(3):fjhi+shift(3))
    f3 => f%p(filo:fihi,fjlo+shift(4):fjhi+shift(4))
    f4 => f%p(filo:fihi,fjlo+shift(5):fjhi+shift(5))
  case default
    print*,"Error: reconst_component_scalar_wenoz5()"
    print*,"No such ireconst: ",axis
    stop
  end select 

  do j=lo%d(2), hi%d(2)
    do i=lo%d(1), hi%d(1)
      is0 = (f0(i,j) - 2.d0*f1(i,j) + f2(i,j))**2*13.d0/12.d0 &
           +(f1(i,j) - 4.d0*f2(i,j) + 3.d0*f3(i,j))**2/4.d0

      is1 = (f1(i,j) - 2.d0*f2(i,j)+f3(i,j))**2*13.d0/12.d0 &
           +(f1(i,j) - f3(i,j))**2/4.d0

      is2 = (f2(i,j) - 2.d0*f3(i,j) + f4(i,j))**2*13.d0/12.d0 &
           +(3.d0*f2(i,j) - 4.d0*f3(i,j) + f4(i,j))**2/4.d0

      tau5 = abs(is2-is0)
      aa0  = c30*(1.d0 + (tau5/(ss+is0))**2)
      aa1  = c31*(1.d0 + (tau5/(ss+is1))**2)
      aa2  = c32*(1.d0 + (tau5/(ss+is2))**2)
      aa_sum = 1.d0/(aa0+aa1+aa2)
    
      w0 = aa0*aa_sum
      w1 = aa1*aa_sum
      w2 = aa2*aa_sum 
      
      q30 = a300*f0(i,j) + a301*f1(i,j) + a302*f2(i,j)
      q31 = a310*f1(i,j) + a311*f3(i,j) + a312*f3(i,j)
      q32 = a320*f2(i,j) + a321*f3(i,j) + a322*f4(i,j)

      h%p(i,j) = w0*q30 + w1*q31 + w2*q32
    enddo
  enddo
end subroutine reconst_component_scalar_wenoz5

module subroutine reconst_component_vector_wenoz5(h,f,lo,hi,padding,bias,axis)
  implicit none
  type(rvectorfield_t),intent(inout) :: h
  type(rvectorfield_t),target,intent(in) :: f
  type(fvectorshape_t),intent(in) :: lo,hi,padding
  integer,intent(in) :: bias,axis

  double precision,dimension(:,:,:),pointer :: f0
  double precision,dimension(:,:,:),pointer :: f1
  double precision,dimension(:,:,:),pointer :: f2
  double precision,dimension(:,:,:),pointer :: f3
  double precision,dimension(:,:,:),pointer :: f4
  integer :: fnlo,fnhi, &
             filo,fihi, &
             fjlo,fjhi

  double precision,dimension(lo%d(1):hi%d(1)) :: is0,is1,is2
  double precision,dimension(lo%d(1):hi%d(1)):: aa0,aa1,aa2,aa_sum
  double precision,dimension(lo%d(1):hi%d(1)):: tau5
  double precision,dimension(lo%d(1):hi%d(1)) :: w0,w1,w2
  double precision,dimension(lo%d(1):hi%d(1)) :: q30,q31,q32
  integer :: i,j
  integer,dimension(5) :: shift
  if(bias.eq.bias_upwind) then
    shift = (/-2,-1,0,1,2/)
  else 
    shift = (/3,2,1,0,-1/)
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
    f3 => f%p(fnlo:fnhi,filo+shift(4):fihi+shift(4),fjlo:fjhi)
    f4 => f%p(fnlo:fnhi,filo+shift(5):fihi+shift(5),fjlo:fjhi)
  case(axis_eta)
    f0 => f%p(fnlo:fnhi,filo:fihi,fjlo+shift(1):fjhi+shift(1))
    f1 => f%p(fnlo:fnhi,filo:fihi,fjlo+shift(2):fjhi+shift(2))
    f2 => f%p(fnlo:fnhi,filo:fihi,fjlo+shift(3):fjhi+shift(3))
    f3 => f%p(fnlo:fnhi,filo:fihi,fjlo+shift(4):fjhi+shift(4))
    f4 => f%p(fnlo:fnhi,filo:fihi,fjlo+shift(5):fjhi+shift(5))
  case default
    print*,"Error: reconst_component_vector_wenoz5()"
    print*,"No such ireconst: ",axis
    stop
  end select 

  do j=lo%d(3), hi%d(3)
    do i=lo%d(2), hi%d(2)
      is0(:) = (f0(:,i,j)-2.d0*f1(:,i,j)+f2(:,i,j))**2*13.d0/12.d0 &
              +(f0(:,i,j)-4.d0*f1(:,i,j)+3.d0*f2(:,i,j))**2/4.d0

      is1(:) = (f1(:,i,j)-2.d0*f2(:,i,j)+f3(:,i,j))**2*13.d0/12.d0 &
              +(f1(:,i,j)-f3(:,i,j))**2/4.d0

      is2(:) = (f2(:,i,j)-2.d0*f3(:,i,j)+f4(:,i,j))**2*13.d0/12.d0 &
              +(3.d0*f2(:,i,j)-4.d0*f3(:,i,j)+f4(:,i,j))**2/4.d0

      tau5(:) = abs(is2(:)-is0(:))
      aa0(:)  = c30*(1.d0 + (tau5(:)/(ss+is0(:)))**2)
      aa1(:)  = c31*(1.d0 + (tau5(:)/(ss+is1(:)))**2)
      aa2(:)  = c32*(1.d0 + (tau5(:)/(ss+is2(:)))**2)
      aa_sum(:) = 1.d0/(aa0(:)+aa1(:)+aa2(:))
    
      w0(:) = aa0(:)*aa_sum(:)
      w1(:) = aa1(:)*aa_sum(:)
      w2(:) = aa2(:)*aa_sum(:)
      
      q30(:) = a300*f0(:,i,j) + a301*f1(:,i,j) + a302*f2(:,i,j)
      q31(:) = a310*f1(:,i,j) + a311*f2(:,i,j) + a312*f3(:,i,j)
      q32(:) = a320*f2(:,i,j) + a321*f3(:,i,j) + a322*f4(:,i,j)

      h%p(lo%d(1):hi%d(1),i,j) = w0(:)*q30(:) + w1(:)*q31(:) + w2(:)*q32(:)
    enddo
  enddo
end subroutine reconst_component_vector_wenoz5

module subroutine reconst_component_matrix_wenoz5(h,f,lo,hi,padding,bias,axis)
  implicit none
  type(rmatrixfield_t),intent(inout) :: h
  type(rmatrixfield_t),target,intent(in) :: f
  type(fmatrixshape_t),intent(in) :: lo,hi,padding
  integer,intent(in) :: bias,axis

  double precision,dimension(:,:,:,:),pointer :: f0
  double precision,dimension(:,:,:,:),pointer :: f1
  double precision,dimension(:,:,:,:),pointer :: f2
  double precision,dimension(:,:,:,:),pointer :: f3
  double precision,dimension(:,:,:,:),pointer :: f4
  integer :: fmlo,fmhi, &
             fnlo,fnhi, &
             filo,fihi, &
             fjlo,fjhi

  double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2)) :: is0,is1,is2
  double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2)) :: aa0,aa1,aa2,aa_sum
  double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2)) :: tau5
  double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2)) :: w0,w1,w2
  double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2)) :: q30,q31,q32
  integer :: i,j
  integer,dimension(5) :: shift
  if(bias.eq.bias_upwind) then
    shift = (/-2,-1,0,1,2/)
  else 
    shift = (/3,2,1,0,-1/)
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
    f3 => f%p(fmlo:fmhi,fnlo:fnhi,filo+shift(4):fihi+shift(4),fjlo:fjhi)
    f4 => f%p(fmlo:fmhi,fnlo:fnhi,filo+shift(5):fihi+shift(5),fjlo:fjhi)
  case(axis_eta)
    f0 => f%p(fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(1):fjhi+shift(1))
    f1 => f%p(fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(2):fjhi+shift(2))
    f2 => f%p(fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(3):fjhi+shift(3))
    f3 => f%p(fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(4):fjhi+shift(4))
    f4 => f%p(fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(5):fjhi+shift(5))
  case default
    print*,"Error: reconst_component_matrix_wenoz5()"
    print*,"No such ireconst: ",axis
    stop
  end select 
  do j = lo%d(4), hi%d(4)
    do i = lo%d(3), hi%d(3)
      is0(:,:) = (f0(:,:,i,j)-2.d0*f2(:,:,i,j)+f3(:,:,i,j))**2*13.d0/12.d0 &
                +(f0(:,:,i,j)-4.d0*f1(:,:,i,j)+3.d0*f2(:,:,i,j))**2/4.d0

      is1(:,:) = (f1(:,:,i,j)-2.d0*f2(:,:,i,j)+f3(:,:,i,j))**2*13.d0/12.d0 &
                +(f1(:,:,i,j)-f3(:,:,i,j))**2/4.d0

      is2(:,:) = (f2(:,:,i,j)-2.d0*f3(:,:,i,j)+f4(:,:,i,j))**2*13.d0/12.d0 &
                +(3.d0*f2(:,:,i,j)-4.d0*f3(:,:,i,j)+f4(:,:,i,j))**2/4.d0

      tau5(:,:) = abs(is2(:,:)-is0(:,:))
      aa0(:,:)  = c30*(1.d0 + (tau5(:,:)/(ss+is0(:,:)))**2)
      aa1(:,:)  = c31*(1.d0 + (tau5(:,:)/(ss+is1(:,:)))**2)
      aa2(:,:)  = c32*(1.d0 + (tau5(:,:)/(ss+is2(:,:)))**2)
      aa_sum(:,:) = 1.d0/(aa0(:,:)+aa1(:,:)+aa2(:,:))
    
      w0(:,:) = aa0(:,:)*aa_sum(:,:)
      w1(:,:) = aa1(:,:)*aa_sum(:,:)
      w2(:,:) = aa2(:,:)*aa_sum(:,:)
      
      q30(:,:) = a300*f0(:,:,i,j) + a301*f1(:,:,i,j) + a302*f2(:,:,i,j)
      q31(:,:) = a310*f1(:,:,i,j) + a311*f2(:,:,i,j) + a312*f3(:,:,i,j)
      q32(:,:) = a320*f2(:,:,i,j) + a321*f3(:,:,i,j) + a322*f4(:,:,i,j)

      h%p(lo%d(1):hi%d(1),lo%d(2):hi%d(2),i,j) = w0(:,:)*q30(:,:) &
                                                +w1(:,:)*q31(:,:) &
                                                +w2(:,:)*q32(:,:)
    enddo
  enddo
end subroutine reconst_component_matrix_wenoz5

module subroutine reconst_component_tensor_wenoz5(h,f,lo,hi,padding,bias,axis)
  implicit none
  type(rtensorfield_t),intent(inout) :: h
  type(rtensorfield_t),target,intent(in) :: f
  type(ftensorshape_t),intent(in) :: lo,hi,padding
  integer,intent(in) :: bias,axis
  
  double precision,dimension(:,:,:,:,:),pointer :: f0
  double precision,dimension(:,:,:,:,:),pointer :: f1
  double precision,dimension(:,:,:,:,:),pointer :: f2
  double precision,dimension(:,:,:,:,:),pointer :: f3
  double precision,dimension(:,:,:,:,:),pointer :: f4
  integer :: fllo,flhi, &
             fmlo,fmhi, &
             fnlo,fnhi, &
             filo,fihi, &
             fjlo,fjhi

  double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2),lo%d(3):hi%d(3)) :: is0,is1,is2
  double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2),lo%d(3):hi%d(3)) :: aa0,aa1,aa2,aa_sum
  double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2),lo%d(3):hi%d(3)) :: tau5
  double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2),lo%d(3):hi%d(3)) :: w0,w1,w2
  double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2),lo%d(3):hi%d(3)) :: q30,q31,q32
  integer :: i,j
  integer,dimension(5) :: shift
  if(bias.eq.bias_upwind) then
    shift = (/-2,-1,0,1,2/)
  else 
    shift = (/3,2,1,0,-1/)
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
    f3 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo+shift(4):fihi+shift(4),fjlo:fjhi)
    f4 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo+shift(5):fihi+shift(5),fjlo:fjhi)
  case(axis_eta)
    f0 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(1):fjhi+shift(1))
    f1 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(2):fjhi+shift(2))
    f2 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(3):fjhi+shift(3))
    f3 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(4):fjhi+shift(4))
    f4 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(5):fjhi+shift(5))
  case default
    print*,"Error: reconst_component_tensor_wenoz5()"
    print*,"No such ireconst: ",axis
    stop
  end select 
  do j = lo%d(5), hi%d(5)
    do i = lo%d(4), hi%d(4)
      is0(:,:,:) = (f0(:,:,:,i,j)-2.d0*f1(:,:,:,i,j)+f2(:,:,:,i,j))**2*13.d0/12.d0 &
                  +(f0(:,:,:,i,j)-4.d0*f1(:,:,:,i,j)+3.d0*f2(:,:,:,i,j))**2/4.d0

      is1(:,:,:) = (f1(:,:,:,i,j)-2.d0*f2(:,:,:,i,j)+f3(:,:,:,i,j))**2*13.d0/12.d0 &
                  +(f1(:,:,:,i,j)-f3(:,:,:,i,j))**2/4.d0

      is2(:,:,:) = (f2(:,:,:,i,j)-2.d0*f3(:,:,:,i,j)+f4(:,:,:,i,j))**2*13.d0/12.d0 &
                  +(3.d0*f2(:,:,:,i,j)-4.d0*f3(:,:,:,i,j)+f4(:,:,:,i,j))**2/4.d0

      tau5(:,:,:) = abs(is2(:,:,:)-is0(:,:,:))
      aa0(:,:,:)  = c30*(1.d0 + (tau5(:,:,:)/(ss+is0(:,:,:)))**2)
      aa1(:,:,:)  = c31*(1.d0 + (tau5(:,:,:)/(ss+is1(:,:,:)))**2)
      aa2(:,:,:)  = c32*(1.d0 + (tau5(:,:,:)/(ss+is2(:,:,:)))**2)
      aa_sum(:,:,:) = 1.d0/(aa0(:,:,:)+aa1(:,:,:)+aa2(:,:,:))
    
      w0(:,:,:) = aa0(:,:,:)*aa_sum(:,:,:)
      w1(:,:,:) = aa1(:,:,:)*aa_sum(:,:,:)
      w2(:,:,:) = aa2(:,:,:)*aa_sum(:,:,:)
      
      q30(:,:,:) = a300*f0(:,:,:,i,j) + a301*f1(:,:,:,i,j) + a302*f2(:,:,:,i,j)
      q31(:,:,:) = a310*f1(:,:,:,i,j) + a311*f2(:,:,:,i,j) + a312*f3(:,:,:,i,j)
      q32(:,:,:) = a320*f2(:,:,:,i,j) + a321*f3(:,:,:,i,j) + a322*f4(:,:,:,i,j)

      h%p(lo%d(1):hi%d(1),lo%d(2):hi%d(2),lo%d(3):hi%d(3),i,j) = w0(:,:,:)*q30(:,:,:) &
                                                                +w1(:,:,:)*q31(:,:,:) &
                                                                +w2(:,:,:)*q32(:,:,:)
    enddo
  enddo
end subroutine reconst_component_tensor_wenoz5

module subroutine reconst_component_scalar_wenoe5(h,f,lo,hi,padding,bias,axis)
    implicit none
    type(rscalarfield_t),intent(inout) :: h
    type(rscalarfield_t),target,intent(in) :: f
    type(fscalarshape_t),intent(in) :: lo,hi,padding
    integer,intent(in) :: bias,axis
  
    double precision,dimension(:,:),pointer :: f0
    double precision,dimension(:,:),pointer :: f1
    double precision,dimension(:,:),pointer :: f2
    double precision,dimension(:,:),pointer :: f3
    double precision,dimension(:,:),pointer :: f4
    integer :: filo,fihi, &
               fjlo,fjhi
  
    double precision :: is0,is1,is2
    double precision :: fpp,fppp
    double precision :: tau5
    double precision :: aa0,aa1,aa2,aa_sum
    double precision :: w0,w1,w2
    double precision :: q30,q31,q32
    integer :: i,j
    integer,dimension(5) :: shift
  
    if(bias.eq.bias_upwind) then
      shift = (/-2,-1,0,1,2/)
    else 
      shift = (/3,2,1,0,-1/)
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
      f3 => f%p(filo+shift(4):fihi+shift(4),fjlo:fjhi)
      f4 => f%p(filo+shift(5):fihi+shift(5),fjlo:fjhi)
    case(axis_eta)
      f0 => f%p(filo:fihi,fjlo+shift(1):fjhi+shift(1))
      f1 => f%p(filo:fihi,fjlo+shift(2):fjhi+shift(2))
      f2 => f%p(filo:fihi,fjlo+shift(3):fjhi+shift(3))
      f3 => f%p(filo:fihi,fjlo+shift(4):fjhi+shift(4))
      f4 => f%p(filo:fihi,fjlo+shift(5):fjhi+shift(5))
    case default
      print*,"Error: reconst_component_scalar_wenoe5()"
      print*,"No such ireconst: ",axis
      stop
    end select 
  
    do j=lo%d(2), hi%d(2)
      do i=lo%d(1), hi%d(1)
        is0 = (f0(i,j) - 2.d0*f1(i,j) + f2(i,j))**2*13.d0/12.d0 &
             +(f1(i,j) - 4.d0*f2(i,j) + 3.d0*f3(i,j))**2/4.d0
  
        is1 = (f1(i,j) - 2.d0*f2(i,j)+f3(i,j))**2*13.d0/12.d0 &
             +(f1(i,j) - f3(i,j))**2/4.d0
  
        is2 = (f2(i,j) - 2.d0*f3(i,j) + f4(i,j))**2*13.d0/12.d0 &
             +(3.d0*f2(i,j) - 4.d0*f3(i,j) + f4(i,j))**2/4.d0
        
        fpp  = (f0(i,j)-4.0d0*f1(i,j)+6.0d0*f2(i,j) -4.0d0*f3(i,j)+f4(i,j))*0.6988d0 &
              +(f1(i,j)-2.0d0*f2(i,j)+f3(i,j))
        
        fppp = 0.5*(-f0(i,j) + 2.0*f1(i,j) - 2.0*f3(i,j) + f4(i,j))

        tau5 = abs(is0 - is2 + 13.0/3.0*fpp*fppp)
        
        aa0  = c30*(1.d0 + (tau5/(ss+is0))**2)
        aa1  = c31*(1.d0 + (tau5/(ss+is1))**2)
        aa2  = c32*(1.d0 + (tau5/(ss+is2))**2)
        aa_sum = 1.d0/(aa0+aa1+aa2)
      
        w0 = aa0*aa_sum
        w1 = aa1*aa_sum
        w2 = aa2*aa_sum 
        
        q30 = a300*f0(i,j) + a301*f1(i,j) + a302*f2(i,j)
        q31 = a310*f1(i,j) + a311*f3(i,j) + a312*f3(i,j)
        q32 = a320*f2(i,j) + a321*f3(i,j) + a322*f4(i,j)
  
        h%p(i,j) = w0*q30 + w1*q31 + w2*q32
      enddo
    enddo
  end subroutine reconst_component_scalar_wenoe5
  
  module subroutine reconst_component_vector_wenoe5(h,f,lo,hi,padding,bias,axis)
    implicit none
    type(rvectorfield_t),intent(inout) :: h
    type(rvectorfield_t),target,intent(in) :: f
    type(fvectorshape_t),intent(in) :: lo,hi,padding
    integer,intent(in) :: bias,axis
  
    double precision,dimension(:,:,:),pointer :: f0
    double precision,dimension(:,:,:),pointer :: f1
    double precision,dimension(:,:,:),pointer :: f2
    double precision,dimension(:,:,:),pointer :: f3
    double precision,dimension(:,:,:),pointer :: f4
    integer :: fnlo,fnhi, &
               filo,fihi, &
               fjlo,fjhi
  
    double precision,dimension(lo%d(1):hi%d(1)) :: is0,is1,is2
    double precision,dimension(lo%d(1):hi%d(1)) :: fpp,fppp
    double precision,dimension(lo%d(1):hi%d(1)) :: aa0,aa1,aa2,aa_sum
    double precision,dimension(lo%d(1):hi%d(1)) :: tau5
    double precision,dimension(lo%d(1):hi%d(1)) :: w0,w1,w2
    double precision,dimension(lo%d(1):hi%d(1)) :: q30,q31,q32
    integer :: i,j
    integer,dimension(5) :: shift
    if(bias.eq.bias_upwind) then
      shift = (/-2,-1,0,1,2/)
    else 
      shift = (/3,2,1,0,-1/)
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
      f3 => f%p(fnlo:fnhi,filo+shift(4):fihi+shift(4),fjlo:fjhi)
      f4 => f%p(fnlo:fnhi,filo+shift(5):fihi+shift(5),fjlo:fjhi)
    case(axis_eta)
      f0 => f%p(fnlo:fnhi,filo:fihi,fjlo+shift(1):fjhi+shift(1))
      f1 => f%p(fnlo:fnhi,filo:fihi,fjlo+shift(2):fjhi+shift(2))
      f2 => f%p(fnlo:fnhi,filo:fihi,fjlo+shift(3):fjhi+shift(3))
      f3 => f%p(fnlo:fnhi,filo:fihi,fjlo+shift(4):fjhi+shift(4))
      f4 => f%p(fnlo:fnhi,filo:fihi,fjlo+shift(5):fjhi+shift(5))
    case default
      print*,"Error: reconst_component_vector_wenoe5()"
      print*,"No such ireconst: ",axis
      stop
    end select 
  
    do j=lo%d(3), hi%d(3)
      do i=lo%d(2), hi%d(2)
        is0(:) = (f0(:,i,j)-2.d0*f1(:,i,j)+f2(:,i,j))**2*13.d0/12.d0 &
                +(f0(:,i,j)-4.d0*f1(:,i,j)+3.d0*f2(:,i,j))**2/4.d0
  
        is1(:) = (f1(:,i,j)-2.d0*f2(:,i,j)+f3(:,i,j))**2*13.d0/12.d0 &
                +(f1(:,i,j)-f3(:,i,j))**2/4.d0
  
        is2(:) = (f2(:,i,j)-2.d0*f3(:,i,j)+f4(:,i,j))**2*13.d0/12.d0 &
                +(3.d0*f2(:,i,j)-4.d0*f3(:,i,j)+f4(:,i,j))**2/4.d0
  
        fpp(:) = (f0(:,i,j) - 4.0d0*f1(:,i,j) + 6.0d0*f2(:,i,j) - 4.0d0*f3(:,i,j)+f4(:,i,j))*0.6988d0 &
                +(f1(:,i,j) - 2.0d0*f2(:,i,j) + f3(:,i,j))

        fppp(:) = 0.5*(-f0(:,i,j) + 2.0*f1(:,i,j) - 2.0*f3(:,i,j) + f4(:,i,j))
  
        tau5(:) = abs(is0(:) - is2(:) + 13.0/3.0*fpp(:)*fppp(:))
        aa0(:)  = c30*(1.d0 + (tau5(:)/(ss+is0(:)))**2)
        aa1(:)  = c31*(1.d0 + (tau5(:)/(ss+is1(:)))**2)
        aa2(:)  = c32*(1.d0 + (tau5(:)/(ss+is2(:)))**2)
        aa_sum(:) = 1.d0/(aa0(:)+aa1(:)+aa2(:))
      
        w0(:) = aa0(:)*aa_sum(:)
        w1(:) = aa1(:)*aa_sum(:)
        w2(:) = aa2(:)*aa_sum(:)
        
        q30(:) = a300*f0(:,i,j) + a301*f1(:,i,j) + a302*f2(:,i,j)
        q31(:) = a310*f1(:,i,j) + a311*f2(:,i,j) + a312*f3(:,i,j)
        q32(:) = a320*f2(:,i,j) + a321*f3(:,i,j) + a322*f4(:,i,j)
  
        h%p(lo%d(1):hi%d(1),i,j) = w0(:)*q30(:) + w1(:)*q31(:) + w2(:)*q32(:)
      enddo
    enddo
  end subroutine reconst_component_vector_wenoe5
  
  module subroutine reconst_component_matrix_wenoe5(h,f,lo,hi,padding,bias,axis)
    implicit none
    type(rmatrixfield_t),intent(inout) :: h
    type(rmatrixfield_t),target,intent(in) :: f
    type(fmatrixshape_t),intent(in) :: lo,hi,padding
    integer,intent(in) :: bias,axis
  
    double precision,dimension(:,:,:,:),pointer :: f0
    double precision,dimension(:,:,:,:),pointer :: f1
    double precision,dimension(:,:,:,:),pointer :: f2
    double precision,dimension(:,:,:,:),pointer :: f3
    double precision,dimension(:,:,:,:),pointer :: f4
    integer :: fmlo,fmhi, &
               fnlo,fnhi, &
               filo,fihi, &
               fjlo,fjhi
  
    double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2)) :: is0,is1,is2
    double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2)) :: fpp,fppp
    double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2)) :: aa0,aa1,aa2,aa_sum
    double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2)) :: tau5
    double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2)) :: w0,w1,w2
    double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2)) :: q30,q31,q32
    integer :: i,j
    integer,dimension(5) :: shift
    if(bias.eq.bias_upwind) then
      shift = (/-2,-1,0,1,2/)
    else 
      shift = (/3,2,1,0,-1/)
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
      f3 => f%p(fmlo:fmhi,fnlo:fnhi,filo+shift(4):fihi+shift(4),fjlo:fjhi)
      f4 => f%p(fmlo:fmhi,fnlo:fnhi,filo+shift(5):fihi+shift(5),fjlo:fjhi)
    case(axis_eta)
      f0 => f%p(fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(1):fjhi+shift(1))
      f1 => f%p(fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(2):fjhi+shift(2))
      f2 => f%p(fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(3):fjhi+shift(3))
      f3 => f%p(fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(4):fjhi+shift(4))
      f4 => f%p(fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(5):fjhi+shift(5))
    case default
      print*,"Error: reconst_component_matrix_wenoe5()"
      print*,"No such ireconst: ",axis
      stop
    end select 
    do j = lo%d(4), hi%d(4)
      do i = lo%d(3), hi%d(3)
        is0(:,:) = (f0(:,:,i,j)-2.d0*f2(:,:,i,j)+f3(:,:,i,j))**2*13.d0/12.d0 &
                  +(f0(:,:,i,j)-4.d0*f1(:,:,i,j)+3.d0*f2(:,:,i,j))**2/4.d0
  
        is1(:,:) = (f1(:,:,i,j)-2.d0*f2(:,:,i,j)+f3(:,:,i,j))**2*13.d0/12.d0 &
                  +(f1(:,:,i,j)-f3(:,:,i,j))**2/4.d0
  
        is2(:,:) = (f2(:,:,i,j)-2.d0*f3(:,:,i,j)+f4(:,:,i,j))**2*13.d0/12.d0 &
                  +(3.d0*f2(:,:,i,j)-4.d0*f3(:,:,i,j)+f4(:,:,i,j))**2/4.d0
  
        fpp(:,:) = (f0(:,:,i,j)-4.0d0*f1(:,:,i,j)+6.0d0*f2(:,:,i,j) -4.0d0*f3(:,:,i,j)+f4(:,:,i,j))*0.6988d0 &
                  +(f1(:,:,i,j)-2.0d0*f2(:,:,i,j)+f3(:,:,i,j))

        fppp(:,:) = 0.5*(-f0(:,:,i,j) + 2.0*f1(:,:,i,j) - 2.0*f3(:,:,i,j) + f4(:,:,i,j))
  
        tau5(:,:) = abs(is0(:,:) - is2(:,:) + 13.0/3.0*fpp(:,:)*fppp(:,:))

        aa0(:,:)  = c30*(1.d0 + (tau5(:,:)/(ss+is0(:,:)))**2)
        aa1(:,:)  = c31*(1.d0 + (tau5(:,:)/(ss+is1(:,:)))**2)
        aa2(:,:)  = c32*(1.d0 + (tau5(:,:)/(ss+is2(:,:)))**2)
        aa_sum(:,:) = 1.d0/(aa0(:,:)+aa1(:,:)+aa2(:,:))
      
        w0(:,:) = aa0(:,:)*aa_sum(:,:)
        w1(:,:) = aa1(:,:)*aa_sum(:,:)
        w2(:,:) = aa2(:,:)*aa_sum(:,:)
        
        q30(:,:) = a300*f0(:,:,i,j) + a301*f1(:,:,i,j) + a302*f2(:,:,i,j)
        q31(:,:) = a310*f1(:,:,i,j) + a311*f2(:,:,i,j) + a312*f3(:,:,i,j)
        q32(:,:) = a320*f2(:,:,i,j) + a321*f3(:,:,i,j) + a322*f4(:,:,i,j)
  
        h%p(lo%d(1):hi%d(1),lo%d(2):hi%d(2),i,j) = w0(:,:)*q30(:,:) &
                                                  +w1(:,:)*q31(:,:) &
                                                  +w2(:,:)*q32(:,:)
      enddo
    enddo
  end subroutine reconst_component_matrix_wenoe5
  
  module subroutine reconst_component_tensor_wenoe5(h,f,lo,hi,padding,bias,axis)
    implicit none
    type(rtensorfield_t),intent(inout) :: h
    type(rtensorfield_t),target,intent(in) :: f
    type(ftensorshape_t),intent(in) :: lo,hi,padding
    integer,intent(in) :: bias,axis
    
    double precision,dimension(:,:,:,:,:),pointer :: f0
    double precision,dimension(:,:,:,:,:),pointer :: f1
    double precision,dimension(:,:,:,:,:),pointer :: f2
    double precision,dimension(:,:,:,:,:),pointer :: f3
    double precision,dimension(:,:,:,:,:),pointer :: f4
    integer :: fllo,flhi, &
               fmlo,fmhi, &
               fnlo,fnhi, &
               filo,fihi, &
               fjlo,fjhi
  
    double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2),lo%d(3):hi%d(3)) :: is0,is1,is2
    double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2),lo%d(3):hi%d(3)) :: fpp,fppp
    double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2),lo%d(3):hi%d(3)) :: aa0,aa1,aa2,aa_sum
    double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2),lo%d(3):hi%d(3)) :: tau5
    double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2),lo%d(3):hi%d(3)) :: w0,w1,w2
    double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2),lo%d(3):hi%d(3)) :: q30,q31,q32
    integer :: i,j
    integer,dimension(5) :: shift
    if(bias.eq.bias_upwind) then
      shift = (/-2,-1,0,1,2/)
    else 
      shift = (/3,2,1,0,-1/)
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
      f3 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo+shift(4):fihi+shift(4),fjlo:fjhi)
      f4 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo+shift(5):fihi+shift(5),fjlo:fjhi)
    case(axis_eta)
      f0 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(1):fjhi+shift(1))
      f1 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(2):fjhi+shift(2))
      f2 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(3):fjhi+shift(3))
      f3 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(4):fjhi+shift(4))
      f4 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(5):fjhi+shift(5))
    case default
      print*,"Error: reconst_component_tensor_wenoe5()"
      print*,"No such ireconst: ",axis
      stop
    end select 
    do j = lo%d(5), hi%d(5)
      do i = lo%d(4), hi%d(4)
        is0(:,:,:) = (f0(:,:,:,i,j)-2.d0*f1(:,:,:,i,j)+f2(:,:,:,i,j))**2*13.d0/12.d0 &
                    +(f0(:,:,:,i,j)-4.d0*f1(:,:,:,i,j)+3.d0*f2(:,:,:,i,j))**2/4.d0
  
        is1(:,:,:) = (f1(:,:,:,i,j)-2.d0*f2(:,:,:,i,j)+f3(:,:,:,i,j))**2*13.d0/12.d0 &
                    +(f1(:,:,:,i,j)-f3(:,:,:,i,j))**2/4.d0
  
        is2(:,:,:) = (f2(:,:,:,i,j)-2.d0*f3(:,:,:,i,j)+f4(:,:,:,i,j))**2*13.d0/12.d0 &
                    +(3.d0*f2(:,:,:,i,j)-4.d0*f3(:,:,:,i,j)+f4(:,:,:,i,j))**2/4.d0
  
        fpp(:,:,:) = (f0(:,:,:,i,j)-4.0d0*f1(:,:,:,i,j)+6.0d0*f2(:,:,:,i,j) -4.0d0*f3(:,:,:,i,j)+f4(:,:,:,i,j))*0.6988d0 &
                    +(f1(:,:,:,i,j)-2.0d0*f2(:,:,:,i,j)+f3(:,:,:,i,j))

        fppp(:,:,:) = 0.5*(-f0(:,:,:,i,j) + 2.0*f1(:,:,:,i,j) - 2.0*f3(:,:,:,i,j) + f4(:,:,:,i,j))
  
        tau5(:,:,:) = abs(is0(:,:,:) - is2(:,:,:) + 13.0/3.0*fpp(:,:,:)*fppp(:,:,:))
        aa0(:,:,:)  = c30*(1.d0 + (tau5(:,:,:)/(ss+is0(:,:,:)))**2)
        aa1(:,:,:)  = c31*(1.d0 + (tau5(:,:,:)/(ss+is1(:,:,:)))**2)
        aa2(:,:,:)  = c32*(1.d0 + (tau5(:,:,:)/(ss+is2(:,:,:)))**2)
        aa_sum(:,:,:) = 1.d0/(aa0(:,:,:)+aa1(:,:,:)+aa2(:,:,:))
      
        w0(:,:,:) = aa0(:,:,:)*aa_sum(:,:,:)
        w1(:,:,:) = aa1(:,:,:)*aa_sum(:,:,:)
        w2(:,:,:) = aa2(:,:,:)*aa_sum(:,:,:)
        
        q30(:,:,:) = a300*f0(:,:,:,i,j) + a301*f1(:,:,:,i,j) + a302*f2(:,:,:,i,j)
        q31(:,:,:) = a310*f1(:,:,:,i,j) + a311*f2(:,:,:,i,j) + a312*f3(:,:,:,i,j)
        q32(:,:,:) = a320*f2(:,:,:,i,j) + a321*f3(:,:,:,i,j) + a322*f4(:,:,:,i,j)
  
        h%p(lo%d(1):hi%d(1),lo%d(2):hi%d(2),lo%d(3):hi%d(3),i,j) = w0(:,:,:)*q30(:,:,:) &
                                                                  +w1(:,:,:)*q31(:,:,:) &
                                                                  +w2(:,:,:)*q32(:,:,:)
      enddo
    enddo
  end subroutine reconst_component_tensor_wenoe5

  module subroutine cmpt_commweight_wenojs5(omega,f,lo,hi,padding,bias,axis)
    implicit none
    type(rvectorfield_t),intent(inout)     :: omega
    type(rscalarfield_t),target,intent(in) :: f
    type(fscalarshape_t),intent(in)        :: lo,hi,padding
    integer,intent(in)                     :: bias,axis
  
    double precision,dimension(:,:),pointer :: f0
    double precision,dimension(:,:),pointer :: f1
    double precision,dimension(:,:),pointer :: f2
    double precision,dimension(:,:),pointer :: f3
    double precision,dimension(:,:),pointer :: f4
    integer :: filo,fihi, &
               fjlo,fjhi
  
    double precision :: is0,is1,is2
    double precision :: aa0,aa1,aa2,aa_sum
    double precision :: w0,w1,w2
    integer :: i,j
    integer,dimension(5) :: shift
  
    if(bias.eq.bias_upwind) then
      shift = (/-2,-1,0,1,2/)
    else 
      shift = (/3,2,1,0,-1/)
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
      f3 => f%p(filo+shift(4):fihi+shift(4),fjlo:fjhi)
      f4 => f%p(filo+shift(5):fihi+shift(5),fjlo:fjhi)
    case(axis_eta)
      f0 => f%p(filo:fihi,fjlo+shift(1):fjhi+shift(1))
      f1 => f%p(filo:fihi,fjlo+shift(2):fjhi+shift(2))
      f2 => f%p(filo:fihi,fjlo+shift(3):fjhi+shift(3))
      f3 => f%p(filo:fihi,fjlo+shift(4):fjhi+shift(4))
      f4 => f%p(filo:fihi,fjlo+shift(5):fjhi+shift(5))
    case default
      print*,"Error: cmpt_commweight_wenojs5()"
      print*,"No such ireconst: ",axis
      stop
    end select 
  
    do j=lo%d(2), hi%d(2)
      do i=lo%d(1), hi%d(1)
        is0 = (f0(i,j) - 2.d0*f1(i,j) + f2(i,j))**2*13.d0/12.d0 &
             +(f1(i,j) - 4.d0*f2(i,j) + 3.d0*f3(i,j))**2/4.d0
  
        is1 = (f1(i,j) - 2.d0*f2(i,j)+f3(i,j))**2*13.d0/12.d0 &
             +(f1(i,j) - f3(i,j))**2/4.d0
  
        is2 = (f2(i,j) - 2.d0*f3(i,j) + f4(i,j))**2*13.d0/12.d0 &
             +(3.d0*f2(i,j) - 4.d0*f3(i,j) + f4(i,j))**2/4.d0
  
        aa0 = c30/(ss+is0)**2
        aa1 = c31/(ss+is1)**2
        aa2 = c32/(ss+is2)**2
        aa_sum = 1.d0/(aa0+aa1+aa2)
      
        w0 = aa0*aa_sum
        w1 = aa1*aa_sum
        w2 = aa2*aa_sum 
        
        omega%p(:,i,j) = (/w0,w1,w2/)
      enddo
    enddo
  end subroutine cmpt_commweight_wenojs5
  
  module subroutine cmpt_commweight_wenoz5(omega,f,lo,hi,padding,bias,axis)
    implicit none
    type(rvectorfield_t),intent(inout)     :: omega
    type(rscalarfield_t),target,intent(in) :: f
    type(fscalarshape_t),intent(in)        :: lo,hi,padding
    integer,intent(in)                     :: bias,axis
  
    double precision,dimension(:,:),pointer :: f0
    double precision,dimension(:,:),pointer :: f1
    double precision,dimension(:,:),pointer :: f2
    double precision,dimension(:,:),pointer :: f3
    double precision,dimension(:,:),pointer :: f4
    integer :: filo,fihi, &
               fjlo,fjhi
  
    double precision :: is0,is1,is2
    double precision :: tau5
    double precision :: aa0,aa1,aa2,aa_sum
    double precision :: w0,w1,w2
    integer :: i,j
    integer,dimension(5) :: shift
  
    if(bias.eq.bias_upwind) then
      shift = (/-2,-1,0,1,2/)
    else 
      shift = (/3,2,1,0,-1/)
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
      f3 => f%p(filo+shift(4):fihi+shift(4),fjlo:fjhi)
      f4 => f%p(filo+shift(5):fihi+shift(5),fjlo:fjhi)
    case(axis_eta)
      f0 => f%p(filo:fihi,fjlo+shift(1):fjhi+shift(1))
      f1 => f%p(filo:fihi,fjlo+shift(2):fjhi+shift(2))
      f2 => f%p(filo:fihi,fjlo+shift(3):fjhi+shift(3))
      f3 => f%p(filo:fihi,fjlo+shift(4):fjhi+shift(4))
      f4 => f%p(filo:fihi,fjlo+shift(5):fjhi+shift(5))
    case default
      print*,"Error: cmpt_commweight_wenoz5()"
      print*,"No such ireconst: ",axis
      stop
    end select 
  
    do j=lo%d(2), hi%d(2)
      do i=lo%d(1), hi%d(1)
        is0 = (f0(i,j) - 2.d0*f1(i,j) + f2(i,j))**2*13.d0/12.d0 &
             +(f1(i,j) - 4.d0*f2(i,j) + 3.d0*f3(i,j))**2/4.d0
  
        is1 = (f1(i,j) - 2.d0*f2(i,j)+f3(i,j))**2*13.d0/12.d0 &
             +(f1(i,j) - f3(i,j))**2/4.d0
  
        is2 = (f2(i,j) - 2.d0*f3(i,j) + f4(i,j))**2*13.d0/12.d0 &
             +(3.d0*f2(i,j) - 4.d0*f3(i,j) + f4(i,j))**2/4.d0
  
        tau5 = abs(is2-is0)
        aa0  = c30*(1.d0 + (tau5/(ss+is0))**2)
        aa1  = c31*(1.d0 + (tau5/(ss+is1))**2)
        aa2  = c32*(1.d0 + (tau5/(ss+is2))**2)
        aa_sum = 1.d0/(aa0+aa1+aa2)
      
        w0 = aa0*aa_sum
        w1 = aa1*aa_sum
        w2 = aa2*aa_sum 
        
        omega%p(:,i,j) = (/w0,w1,w2/)
      enddo
    enddo
  end subroutine cmpt_commweight_wenoz5

  module subroutine cmpt_commweight_wenoe5(omega,f,lo,hi,padding,bias,axis)
    implicit none
    type(rvectorfield_t),intent(inout)     :: omega
    type(rscalarfield_t),target,intent(in) :: f
    type(fscalarshape_t),intent(in)        :: lo,hi,padding
    integer,intent(in)                     :: bias,axis
  
    double precision,dimension(:,:),pointer :: f0
    double precision,dimension(:,:),pointer :: f1
    double precision,dimension(:,:),pointer :: f2
    double precision,dimension(:,:),pointer :: f3
    double precision,dimension(:,:),pointer :: f4
    integer :: filo,fihi, &
               fjlo,fjhi
  
    double precision :: is0,is1,is2
    double precision :: fpp,fppp
    double precision :: tau5
    double precision :: aa0,aa1,aa2,aa_sum
    double precision :: w0,w1,w2
    integer :: i,j
    integer,dimension(5) :: shift
  
    if(bias.eq.bias_upwind) then
      shift = (/-2,-1,0,1,2/)
    else 
      shift = (/3,2,1,0,-1/)
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
      f3 => f%p(filo+shift(4):fihi+shift(4),fjlo:fjhi)
      f4 => f%p(filo+shift(5):fihi+shift(5),fjlo:fjhi)
    case(axis_eta)
      f0 => f%p(filo:fihi,fjlo+shift(1):fjhi+shift(1))
      f1 => f%p(filo:fihi,fjlo+shift(2):fjhi+shift(2))
      f2 => f%p(filo:fihi,fjlo+shift(3):fjhi+shift(3))
      f3 => f%p(filo:fihi,fjlo+shift(4):fjhi+shift(4))
      f4 => f%p(filo:fihi,fjlo+shift(5):fjhi+shift(5))
    case default
      print*,"Error: cmpt_commweight_wenoe5()"
      print*,"No such ireconst: ",axis
      stop
    end select 
  
    do j=lo%d(2), hi%d(2)
      do i=lo%d(1), hi%d(1)
        is0 = (f0(i,j) - 2.d0*f1(i,j) + f2(i,j))**2*13.d0/12.d0 &
             +(f1(i,j) - 4.d0*f2(i,j) + 3.d0*f3(i,j))**2/4.d0
  
        is1 = (f1(i,j) - 2.d0*f2(i,j)+f3(i,j))**2*13.d0/12.d0 &
             +(f1(i,j) - f3(i,j))**2/4.d0
  
        is2 = (f2(i,j) - 2.d0*f3(i,j) + f4(i,j))**2*13.d0/12.d0 &
             +(3.d0*f2(i,j) - 4.d0*f3(i,j) + f4(i,j))**2/4.d0
  
        fpp  = (f0(i,j)-4.0d0*f1(i,j)+6.0d0*f2(i,j) -4.0d0*f3(i,j)+f4(i,j))*0.6988d0 &
              +(f1(i,j)-2.0d0*f2(i,j)+f3(i,j))
        
        fppp = 0.5*(-f0(i,j) + 2.0*f1(i,j) - 2.0*f3(i,j) + f4(i,j))

        tau5 = abs(is0 - is2 + 13.0/3.0*fpp*fppp)
        
        aa0    = c30*(1.d0 + (tau5/(ss+is0))**2)
        aa1    = c31*(1.d0 + (tau5/(ss+is1))**2)
        aa2    = c32*(1.d0 + (tau5/(ss+is2))**2)
        aa_sum = 1.d0/(aa0+aa1+aa2)
      
        w0 = aa0*aa_sum
        w1 = aa1*aa_sum
        w2 = aa2*aa_sum 
        
        omega%p(:,i,j) = (/w0,w1,w2/)
      enddo
    enddo
  end subroutine cmpt_commweight_wenoe5

module subroutine reconst_commweight_scalar_weno5(h,f,omega,lo,hi,padding,bias,axis)
  implicit none
  type(rscalarfield_t),intent(inout)     :: h
  type(rscalarfield_t),target,intent(in) :: f
  type(rvectorfield_t),intent(in)        :: omega
  type(fscalarshape_t),intent(in)        :: lo,hi,padding
  integer,intent(in)                     :: bias,axis

  double precision,dimension(:,:),pointer :: f0
  double precision,dimension(:,:),pointer :: f1
  double precision,dimension(:,:),pointer :: f2
  double precision,dimension(:,:),pointer :: f3
  double precision,dimension(:,:),pointer :: f4
  integer :: filo,fihi, &
             fjlo,fjhi

  double precision :: q30,q31,q32
  integer :: i,j
  integer,dimension(5) :: shift

  if(bias.eq.bias_upwind) then
    shift = (/-2,-1,0,1,2/)
  else 
    shift = (/3,2,1,0,-1/)
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
    f3 => f%p(filo+shift(4):fihi+shift(4),fjlo:fjhi)
    f4 => f%p(filo+shift(5):fihi+shift(5),fjlo:fjhi)
  case(axis_eta)
    f0 => f%p(filo:fihi,fjlo+shift(1):fjhi+shift(1))
    f1 => f%p(filo:fihi,fjlo+shift(2):fjhi+shift(2))
    f2 => f%p(filo:fihi,fjlo+shift(3):fjhi+shift(3))
    f3 => f%p(filo:fihi,fjlo+shift(4):fjhi+shift(4))
    f4 => f%p(filo:fihi,fjlo+shift(5):fjhi+shift(5))
  case default
    print*,"Error: reconst_commweight_scalar_weno5()"
    print*,"No such ireconst: ",axis
    stop
  end select 

  do j=lo%d(2), hi%d(2)
    do i=lo%d(1), hi%d(1)
          
      q30 = a300*f0(i,j) + a301*f1(i,j) + a302*f2(i,j)
      q31 = a310*f1(i,j) + a311*f3(i,j) + a312*f3(i,j)
      q32 = a320*f2(i,j) + a321*f3(i,j) + a322*f4(i,j)

      h%p(i,j) = omega%p(1,i,j)*q30 + omega%p(2,i,j)*q31 + omega%p(3,i,j)*q32
    enddo
  enddo
end subroutine reconst_commweight_scalar_weno5

module subroutine reconst_commweight_vector_weno5(h,f,omega,lo,hi,padding,bias,axis)
  implicit none
  type(rvectorfield_t),intent(inout)     :: h
  type(rvectorfield_t),target,intent(in) :: f
  type(rvectorfield_t),intent(in)        :: omega
  type(fvectorshape_t),intent(in)        :: lo,hi,padding
  integer,intent(in)                     :: bias,axis

  double precision,dimension(:,:,:),pointer :: f0
  double precision,dimension(:,:,:),pointer :: f1
  double precision,dimension(:,:,:),pointer :: f2
  double precision,dimension(:,:,:),pointer :: f3
  double precision,dimension(:,:,:),pointer :: f4
  integer :: fnlo,fnhi, &
             filo,fihi, &
             fjlo,fjhi

  double precision,dimension(lo%d(1):hi%d(1)) :: q30,q31,q32
  integer :: i,j
  integer,dimension(5) :: shift
  if(bias.eq.bias_upwind) then
    shift = (/-2,-1,0,1,2/)
  else 
    shift = (/3,2,1,0,-1/)
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
    f3 => f%p(fnlo:fnhi,filo+shift(4):fihi+shift(4),fjlo:fjhi)
    f4 => f%p(fnlo:fnhi,filo+shift(5):fihi+shift(5),fjlo:fjhi)
  case(axis_eta)
    f0 => f%p(fnlo:fnhi,filo:fihi,fjlo+shift(1):fjhi+shift(1))
    f1 => f%p(fnlo:fnhi,filo:fihi,fjlo+shift(2):fjhi+shift(2))
    f2 => f%p(fnlo:fnhi,filo:fihi,fjlo+shift(3):fjhi+shift(3))
    f3 => f%p(fnlo:fnhi,filo:fihi,fjlo+shift(4):fjhi+shift(4))
    f4 => f%p(fnlo:fnhi,filo:fihi,fjlo+shift(5):fjhi+shift(5))
  case default
    print*,"Error: reconst_commweight_vector_weno5()"
    print*,"No such ireconst: ",axis
    stop
  end select 

  do j=lo%d(3), hi%d(3)
    do i=lo%d(2), hi%d(2)
      
      q30(:) = a300*f0(:,i,j) + a301*f1(:,i,j) + a302*f2(:,i,j)
      q31(:) = a310*f1(:,i,j) + a311*f2(:,i,j) + a312*f3(:,i,j)
      q32(:) = a320*f2(:,i,j) + a321*f3(:,i,j) + a322*f4(:,i,j)

      h%p(lo%d(1):hi%d(1),i,j) =  omega%p(1,i,j)*q30(:) &
                                + omega%p(2,i,j)*q31(:) &
                                + omega%p(3,i,j)*q32(:)
    enddo
  enddo
end subroutine reconst_commweight_vector_weno5

module subroutine reconst_commweight_matrix_weno5(h,f,omega,lo,hi,padding,bias,axis)
  implicit none
  type(rmatrixfield_t),intent(inout)     :: h
  type(rmatrixfield_t),target,intent(in) :: f
  type(rvectorfield_t),intent(in)        :: omega
  type(fmatrixshape_t),intent(in)        :: lo,hi,padding
  integer,intent(in)                     :: bias,axis

  double precision,dimension(:,:,:,:),pointer :: f0
  double precision,dimension(:,:,:,:),pointer :: f1
  double precision,dimension(:,:,:,:),pointer :: f2
  double precision,dimension(:,:,:,:),pointer :: f3
  double precision,dimension(:,:,:,:),pointer :: f4
  integer :: fmlo,fmhi, &
             fnlo,fnhi, &
             filo,fihi, &
             fjlo,fjhi

  double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2)) :: q30,q31,q32
  integer :: i,j
  integer,dimension(5) :: shift
  if(bias.eq.bias_upwind) then
    shift = (/-2,-1,0,1,2/)
  else 
    shift = (/3,2,1,0,-1/)
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
    f3 => f%p(fmlo:fmhi,fnlo:fnhi,filo+shift(4):fihi+shift(4),fjlo:fjhi)
    f4 => f%p(fmlo:fmhi,fnlo:fnhi,filo+shift(5):fihi+shift(5),fjlo:fjhi)
  case(axis_eta)
    f0 => f%p(fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(1):fjhi+shift(1))
    f1 => f%p(fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(2):fjhi+shift(2))
    f2 => f%p(fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(3):fjhi+shift(3))
    f3 => f%p(fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(4):fjhi+shift(4))
    f4 => f%p(fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(5):fjhi+shift(5))
  case default
    print*,"Error: reconst_commweight_matrix_weno5()"
    print*,"No such ireconst: ",axis
    stop
  end select 
  do j = lo%d(4), hi%d(4)
    do i = lo%d(3), hi%d(3)      
      q30(:,:) = a300*f0(:,:,i,j) + a301*f1(:,:,i,j) + a302*f2(:,:,i,j)
      q31(:,:) = a310*f1(:,:,i,j) + a311*f2(:,:,i,j) + a312*f3(:,:,i,j)
      q32(:,:) = a320*f2(:,:,i,j) + a321*f3(:,:,i,j) + a322*f4(:,:,i,j)

      h%p(lo%d(1):hi%d(1),lo%d(2):hi%d(2),i,j) =  omega%p(1,i,j)*q30(:,:) &
                                                + omega%p(2,i,j)*q31(:,:) &
                                                + omega%p(3,i,j)*q32(:,:)
    enddo
  enddo
end subroutine reconst_commweight_matrix_weno5

module subroutine reconst_commweight_tensor_weno5(h,f,omega,lo,hi,padding,bias,axis)
  implicit none
  type(rtensorfield_t),intent(inout)     :: h
  type(rtensorfield_t),target,intent(in) :: f
  type(rvectorfield_t),intent(in)        :: omega
  type(ftensorshape_t),intent(in)        :: lo,hi,padding
  integer,intent(in)                     :: bias,axis
  
  double precision,dimension(:,:,:,:,:),pointer :: f0
  double precision,dimension(:,:,:,:,:),pointer :: f1
  double precision,dimension(:,:,:,:,:),pointer :: f2
  double precision,dimension(:,:,:,:,:),pointer :: f3
  double precision,dimension(:,:,:,:,:),pointer :: f4
  integer :: fllo,flhi, &
             fmlo,fmhi, &
             fnlo,fnhi, &
             filo,fihi, &
             fjlo,fjhi

  double precision,dimension(lo%d(1):hi%d(1),lo%d(2):hi%d(2),lo%d(3):hi%d(3)) :: q30,q31,q32
  integer :: i,j
  integer,dimension(5) :: shift
  if(bias.eq.bias_upwind) then
    shift = (/-2,-1,0,1,2/)
  else 
    shift = (/3,2,1,0,-1/)
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
    f3 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo+shift(4):fihi+shift(4),fjlo:fjhi)
    f4 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo+shift(5):fihi+shift(5),fjlo:fjhi)
  case(axis_eta)
    f0 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(1):fjhi+shift(1))
    f1 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(2):fjhi+shift(2))
    f2 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(3):fjhi+shift(3))
    f3 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(4):fjhi+shift(4))
    f4 => f%p(fllo:flhi,fmlo:fmhi,fnlo:fnhi,filo:fihi,fjlo+shift(5):fjhi+shift(5))
  case default
    print*,"Error: reconst_commweight_tensor_weno5()"
    print*,"No such ireconst: ",axis
    stop
  end select 
  do j = lo%d(5), hi%d(5)
    do i = lo%d(4), hi%d(4)      
      q30(:,:,:) = a300*f0(:,:,:,i,j) + a301*f1(:,:,:,i,j) + a302*f2(:,:,:,i,j)
      q31(:,:,:) = a310*f1(:,:,:,i,j) + a311*f2(:,:,:,i,j) + a312*f3(:,:,:,i,j)
      q32(:,:,:) = a320*f2(:,:,:,i,j) + a321*f3(:,:,:,i,j) + a322*f4(:,:,:,i,j)

      h%p(lo%d(1):hi%d(1),lo%d(2):hi%d(2),lo%d(3):hi%d(3),i,j) =  omega%p(1,i,j)*q30(:,:,:) &
                                                                + omega%p(2,i,j)*q31(:,:,:) &
                                                                + omega%p(3,i,j)*q32(:,:,:)
    enddo
  enddo
end subroutine reconst_commweight_tensor_weno5
end submodule