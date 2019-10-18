module gradient
  use data_type
  implicit none

  integer,parameter :: grad_1st_idx    = 1
  integer,parameter :: grad_1st_l1_idx = 11
  integer,parameter :: grad_1st_l0_idx = 10

  integer,parameter :: grad_2nd_idx    = 2
  integer,parameter :: grad_2nd_l2_idx = 22
  integer,parameter :: grad_2nd_l1_idx = 21
  integer,parameter :: grad_2nd_l0_idx = 20

  integer,parameter :: grad_3rd_idx    = 3
  integer,parameter :: grad_3rd_l3_idx = 33
  integer,parameter :: grad_3rd_l2_idx = 32
  integer,parameter :: grad_3rd_l1_idx = 31
  integer,parameter :: grad_3rd_l0_idx = 30

  integer,parameter :: grad_4th_idx    = 4
  integer,parameter :: grad_4th_l4_idx = 44
  integer,parameter :: grad_4th_l3_idx = 43
  integer,parameter :: grad_4th_l2_idx = 42
  integer,parameter :: grad_4th_l1_idx = 41
  integer,parameter :: grad_4th_l0_idx = 40

  integer,parameter :: grad_5th_idx    = 5
  integer,parameter :: grad_5th_l5_idx = 55
  integer,parameter :: grad_5th_l4_idx = 54
  integer,parameter :: grad_5th_l3_idx = 53
  integer,parameter :: grad_5th_l2_idx = 52
  integer,parameter :: grad_5th_l1_idx = 51
  integer,parameter :: grad_5th_l0_idx = 50

  integer,parameter :: grad_6th_idx    = 6
  integer,parameter :: grad_6th_l6_idx = 66
  integer,parameter :: grad_6th_l5_idx = 65
  integer,parameter :: grad_6th_l4_idx = 64
  integer,parameter :: grad_6th_l3_idx = 63
  integer,parameter :: grad_6th_l2_idx = 62
  integer,parameter :: grad_6th_l1_idx = 61
  integer,parameter :: grad_6th_l0_idx = 60
  
  integer,parameter :: grad_7th_idx    = 7
  integer,parameter :: grad_7th_l7_idx = 77
  integer,parameter :: grad_7th_l6_idx = 76
  integer,parameter :: grad_7th_l5_idx = 75
  integer,parameter :: grad_7th_l4_idx = 74
  integer,parameter :: grad_7th_l3_idx = 73
  integer,parameter :: grad_7th_l2_idx = 72
  integer,parameter :: grad_7th_l1_idx = 71
  integer,parameter :: grad_7th_l0_idx = 70

  integer,parameter :: grad_8th_idx    = 8
  integer,parameter :: grad_8th_l8_idx = 88
  integer,parameter :: grad_8th_l7_idx = 87
  integer,parameter :: grad_8th_l6_idx = 86
  integer,parameter :: grad_8th_l5_idx = 85
  integer,parameter :: grad_8th_l4_idx = 84
  integer,parameter :: grad_8th_l3_idx = 83
  integer,parameter :: grad_8th_l2_idx = 82
  integer,parameter :: grad_8th_l1_idx = 81
  integer,parameter :: grad_8th_l0_idx = 80

  ! Finite difference coefficiences
  ! name (order)_(left_point,right_poitn) 
  ! e.g. a1_10 means first order scheme using f(i-1),...,f(i) 

  ! First order 
  double precision,dimension(2),target :: a1_10 = (/-1.d0,1.d0/)
  double precision,dimension(2),target :: a1_01 = (/-1.d0,1.d0/)

  ! Second order 
  double precision,dimension(3),target :: a2_20 = (/0.5d0,-2.d0,1.5d0/)
  double precision,dimension(3),target :: a2_11 = (/-0.5d0,0.d0,0.5d0/)
  double precision,dimension(3),target :: a2_02 = (/-1.5d0,2.d0,-0.5d0/)
  
  ! Third order 
  double precision,dimension(4),target :: a3_30 = (/1.d0/3.d0,3.d0/2.d0,-3.d0,11.d0/6.d0/)
  double precision,dimension(4),target :: a3_21 = (/1.d0/6.d0,-1.d0,0.5d0,1.d0/3.d0/)
  double precision,dimension(4),target :: a3_12 = (/-1.d0/3.d0,-0.5d0,1.d0,-1.d0/6.d0/)
  double precision,dimension(4),target :: a3_03 = (/-11.d0/6.d0,3.d0,-3.d0/2.d0,1.d0/3.d0/)

  ! Fourth order 
  double precision,dimension(5),target :: a4_40 = (/1.d0/4.d0,-4.d0/3.d0,3.d0,-4.d0,25.d0/12.d0/)
  double precision,dimension(5),target :: a4_31 = (/-1.d0/12.d0,1.d0/2.d0,-3.d0/2.d0,5.d0/6.d0,1.d0/4.d0/)
  double precision,dimension(5),target :: a4_22 = (/1.d0/12.d0,-2.d0/3.d0,0.d0,2.d0/3.d0,-1.d0/12.d0/)
  double precision,dimension(5),target :: a4_13 = (/-1.d0/4.d0,-5.d0/6.d0,3.d0/2.d0,-1.d0/2.d0,1.d0/12.d0/)
  double precision,dimension(5),target :: a4_04 = (/-25.d0/12.d0,4.d0,-3.d0,4.d0/3.d0,-1.d0/4.d0/)
  
  ! Fifth order 
  double precision,dimension(6),target :: a5_50 = (/-1.d0/5.d0,5.d0/4.d0,-10.d0/3.d0,5.d0,-5.d0, &
                                                       137.d0/60.d0/)
  double precision,dimension(6),target :: a5_41 = (/1.d0/20.d0,-1.d0/3.d0,1.d0,-2.d0,13.d0/12.d0, &
                                                       1.d0/5.d0/)
  double precision,dimension(6),target :: a5_32 = (/-1.d0/30.d0,1.d0/4.d0,-1.d0,1.d0/3.d0,1.d0/2.d0, &
                                                       -1.d0/20.d0/)
  double precision,dimension(6),target :: a5_23 = (/1.d0/20.d0,-1.d0/2.d0,-1.d0/3.d0,1.d0,-1.d0/4.d0, &
                                                       1.d0/30.d0/)
  double precision,dimension(6),target :: a5_14 = (/-1.d0/5.d0,-13.d0/12.d0,2.d0,-1.d0,1.d0/3.d0, &
                                                       -1.d0/20.d0/)
  double precision,dimension(6),target :: a5_05 = (/-137.d0/60.d0,5.d0,-5.d0,10.d0/3.d0,-5.d0/4.d0, &
                                                       1.d0/5.d0/)
  
  ! Sixth order 
  double precision,dimension(7),target :: a6_60 = (/1.d0/6.d0,-6.d0/5.d0,15.d0/4.d0,-20.d0/3.d0, &
                                                       15.d0/2.d0,-6.d0,49.d0/20.d0/)
  double precision,dimension(7),target :: a6_51 = (/-1.d0/30.d0,1.d0/4.d0,-5.d0/6.d0,5.d0/3.d0, &
                                                       -5.d0/2.d0,-77.d0/60.d0,1.d0/6.d0/)
  double precision,dimension(7),target :: a6_42 = (/1.d0/60.d0,-2.d0/15.d0,1.d0/2.d0,-4.d0/3.d0, &
                                                       7.d0/12.d0,2.d0/5.d0,-1.d0/30.d0/)
  double precision,dimension(7),target :: a6_33 = (/-1.d0/60.d0,3.d0/20.d0,-3.d0/4.d0,0.d0,  &
                                                       3.d0/4.d0,-3.d0/20.d0,1.d0/60.d0/)
  double precision,dimension(7),target :: a6_24 = (/1.d0/30.d0,-2.d0/5.d0,-7.d0/12.d0,4.d0/3.d0, &
                                                       -1.d0/2.d0,2.d0/15.d0,-1.d0/60.d0/)
  double precision,dimension(7),target :: a6_15 = (/-1.d0/6.d0,-77.d0/60.d0,5.d0/2.d0,-5.d0/3.d0, &
                                                       5.d0/6.d0,-1.d0/4.d0,1.d0/30.d0/)
  double precision,dimension(7),target :: a6_06 = (/-49.d0/20.d0,6.d0,-15.d0/2.d0,20.d0/3.d0, &
                                                       -15.d0/4.d0,6.d0/5.d0,-1.d0/6.d0/)
  
  ! Seventh order 
  double precision,dimension(8),target :: a7_70 = (/-1.d0/7.d0,7.d0/6.d0,-21.d0/5.d0,35.d0/4.d0, &
                                                       -35.d0/3.d0,21.d0/2.d0,-7.d0,363.d0/140.d0/)
  double precision,dimension(8),target :: a7_61 = (/1.d0/42.d0,-1.d0/5.d0,3.d0/4.d0,-5.d0/3.d0, &
                                                       5.d0/2.d0,-3.d0,29.d0/20.d0,1.d0/7.d0/)
  double precision,dimension(8),target :: a7_52 = (/-1.d0/105.d0,1.d0/12.d0,-1.d0/3.d0,5.d0/6.d0, &
                                                       -5.d0/3.d0,47.d0/60.d0,1.d0/3.d0,-1.d0/42.d0/)
  double precision,dimension(8),target :: a7_43 = (/1.d0/140.d0,-1.d0/15.d0,3.d0/10.d0,-1.d0, &
                                                       1.d0/4.d0,3.d0/5.d0,-1.d0/10.d0,1.d0/105.d0/)
  double precision,dimension(8),target :: a7_34 = (/-1.d0/105.d0,1.d0/10.d0,-3.d0/5.d0,-1.d0/4.d0, &
                                                       1.d0,-3.d0/10.d0,1.d0/15.d0,-1.d0/140.d0/)
  double precision,dimension(8),target :: a7_25 = (/1.d0/42.d0,-1.d0/3.d0,-47.d0/60.d0,5.d0/3.d0, &
                                                       -5.d0/6.d0,1.d0/3.d0,-1.d0/12.d0,1.d0/105.d0/)
  double precision,dimension(8),target :: a7_16 = (/-1.d0/7.d0,-29.d0/20.d0,3.d0,-5.d0/2.d0, &
                                                       5.d0/3.d0,-3.d0/4.d0,1.d0/5.d0,-1.d0/42.d0/)
  double precision,dimension(8),target :: a7_07 = (/-363.d0/140.d0,7.d0,-21.d0/2.d0,35.d0/3.d0, &
                                                       -35.d0/4.d0,21.d0/5.d0,-7.d0/6.d0,1.d0/7.d0/)
                                                       
  ! Eighth order 
  double precision,dimension(9),target :: a8_80 = (/1.d0/8.d0,-8.d0/7.d0,14.d0/3.d0,-56.d0/5.d0, &
                                                       35.d0/2.d0,-56.d0/3.d0,14.d0,-8.d0,761.d0/280.d0/)
  double precision,dimension(9),target :: a8_71 = (/-1.d0/56.d0,1.d0/6.d0,-7.d0/10.d0,7.d0/4.d0, &
                                                       -35.d0/12.d0,7.d0/2.d0,-7.d0/2.d0,223.d0/140.d0,1.d0/8.d0/)
  double precision,dimension(9),target :: a8_62 = (/1.d0/168.d0,-2.d0/35.d0,1.d0/4.d0,-2.d0/3.d0, &
                                                       5.d0/4.d0,-2.d0,19.d0/20.d0,2.d0/7.d0,-1.d0/56.d0/)
  double precision,dimension(9),target :: a8_53 = (/-1.d0/280.d0,1.d0/28.d0,-1.d0/6.d0,1.d0/2.d0, &
                                                       -5.d0/4.d0,9.d0/20.d0,1.d0/2.d0,-1.d0/14.d0,1.d0/168.d0/)
  double precision,dimension(9),target :: a8_44 = (/1.d0/280.d0,-4.d0/105.d0,1.d0/5.d0,-4.d0/5.d0, &
                                                       0.d0,4.d0/5.d0,-1.d0/5.d0,4.d0/105.d0,-1.d0/280.d0/)
  double precision,dimension(9),target :: a8_35 = (/-1.d0/168.d0,1.d0/14.d0,-1.d0/2.d0,-9.d0/20.d0, &
                                                       5.d0/4.d0,-1.d0/2.d0,1.d0/6.d0,-1.d0/28.d0,1.d0/280.d0/)
  double precision,dimension(9),target :: a8_26 = (/1.d0/56.d0,-2.d0/7.d0,-19.d0/20.d0,2.d0, &
                                                       -5.d0/4.d0,2.d0/3.d0,-1.d0/4.d0,2.d0/35.d0,-1.d0/168.d0/)
  double precision,dimension(9),target :: a8_17 = (/-1.d0/8.d0,-223.d0/140.d0,7.d0/2.d0,-7.d0/2.d0, &
                                                       35.d0/12.d0,-7.d0/4.d0,7.d0/10.d0,-1.d0/6.d0,1.d0/56.d0/)
  double precision,dimension(9),target :: a8_08 = (/-761.d0/280.d0,8.d0,-14.d0,56.d0/3.d0, &
                                                       -35.d0/2.d0,56.d0/5.d0,-14.d0/3.d0,8.d0/7.d0,-1.d0/8.d0/)

interface compute_grad
  module procedure compute_grad_homo_scalar
  module procedure compute_grad_homo_vector
  module procedure compute_grad_homo_matrix
  module procedure compute_grad_homo_tensor
  module procedure compute_grad_heter_scalar
  module procedure compute_grad_heter_vector
  module procedure compute_grad_heter_matrix
  module procedure compute_grad_heter_tensor
end interface compute_grad

contains
subroutine compute_grad_homo_scalar(df,f,lo,hi,padding,axis,ischeme)
  implicit none
  integer,intent(in) :: axis,ischeme
  type(fscalarshape_t),intent(in) :: lo,hi,padding
  type(rscalarfield_t),intent(in) :: f
  type(rscalarfield_t),intent(out) :: df
  double precision,dimension(:),pointer :: coeff
  integer :: i,j
  integer :: i0,j0
  integer :: order,bias

  order = ischeme/10
  bias  = mod(ischeme,10)
  call select_coeff(coeff,order,bias)
  
  select case(axis)
    case(axis_xi)
      do j = lo%d(2), hi%d(2)
        do i = lo%d(1), hi%d(1)
          i0 = i + padding%d(1)
          j0 = j + padding%d(2)
          df%p(i,j) = dot_product(f%p(i0-bias:i0+order-bias,j0),coeff(:))
        enddo
      enddo
    case(axis_eta)
      do j = lo%d(2), hi%d(2)
        do i = lo%d(1), hi%d(1)
          i0 = i + padding%d(1)
          j0 = j + padding%d(2)
          df%p(i,j) = dot_product(f%p(i0,j0-bias:j0+order-bias),coeff(:))
        enddo
      enddo
    case default
      write(*,*) "Error: compute_grad_homo_scalar():"
      write(*,*) "Illegal axis:",axis
      stop
  end select
end subroutine compute_grad_homo_scalar

subroutine compute_grad_homo_vector(df,f,lo,hi,padding,axis,ischeme)
  implicit none
  integer,intent(in) :: axis,ischeme
  type(fvectorshape_t),intent(in) :: lo,hi,padding
  type(rvectorfield_t),intent(in) :: f
  type(rvectorfield_t),intent(out) :: df
  double precision,dimension(:),pointer :: coeff
  integer :: i,j,n
  integer :: i0,j0,n0
  integer :: order,bias

  order = ischeme/10
  bias  = mod(ischeme,10)
  call select_coeff(coeff,order,bias)
  
  select case(axis)
    case(axis_xi)
      do j = lo%d(3), hi%d(3)
        do i = lo%d(2), hi%d(2)
          do n = lo%d(1), hi%d(1)
            n0 = n + padding%d(1)
            i0 = i + padding%d(2)
            j0 = j + padding%d(3)
            df%p(n,i,j) = dot_product(f%p(n0,i0-bias:i0+order-bias,j0),coeff(:))
          enddo
        enddo
      enddo
    case(axis_eta)
      do j = lo%d(3), hi%d(3)
        do i = lo%d(2), hi%d(2)
          do n = lo%d(1), hi%d(1)
            n0 = n + padding%d(1)
            i0 = i + padding%d(2)
            j0 = j + padding%d(3)
            df%p(n,i,j) = dot_product(f%p(n0,i0,j0-bias:j0+order-bias),coeff(:))        
          enddo
        enddo
      enddo
    case default
      write(*,*) "Error: compute_grad_homo_vector():"
      write(*,*) "Illegal axis:",axis
      stop
  end select
end subroutine compute_grad_homo_vector

subroutine compute_grad_homo_matrix(df,f,lo,hi,padding,axis,ischeme)
  implicit none
  integer,intent(in) :: axis,ischeme
  type(fmatrixshape_t),intent(in) :: lo,hi,padding
  type(rmatrixfield_t),intent(in) :: f
  type(rmatrixfield_t),intent(out) :: df
  double precision,dimension(:),pointer :: coeff
  integer :: i,j,n,m
  integer :: i0,j0,n0,m0
  integer :: order,bias

  order = ischeme/10
  bias  = mod(ischeme,10)
  call select_coeff(coeff,order,bias)
  
  select case(axis)
    case(axis_xi)
      do j = lo%d(4), hi%d(4)
        do i = lo%d(3), hi%d(3)
          do m = lo%d(2), hi%d(2)
            do n = lo%d(1), hi%d(1)
              n0 = n + padding%d(1)
              m0 = m + padding%d(2)
              i0 = i + padding%d(3)
              j0 = j + padding%d(4)
              df%p(n,m,i,j) = dot_product(f%p(n0,m0,i0-bias:i0+order-bias,j0),coeff(:))            
            enddo
          enddo
        enddo
      enddo
    case(axis_eta)
      do j = lo%d(4), hi%d(4)
        do i = lo%d(3), hi%d(3)
          do m = lo%d(2), hi%d(2)
            do n = lo%d(1), hi%d(1)
              n0 = n + padding%d(1)
              m0 = m + padding%d(2)
              i0 = i + padding%d(3)
              j0 = j + padding%d(4)
              df%p(n,m,i,j) = dot_product(f%p(n0,m0,i0,j0-bias:j0+order-bias),coeff(:))            
            enddo
          enddo
        enddo
      enddo
    case default
      write(*,*) "Error: compute_grad_homo_matrix():"
      write(*,*) "Illegal axis:",axis
      stop
  end select
end subroutine compute_grad_homo_matrix

subroutine compute_grad_homo_tensor(df,f,lo,hi,padding,axis,ischeme)
  implicit none
  integer,intent(in) :: axis,ischeme
  type(rtensorfield_t),intent(in) :: f
  type(ftensorshape_t),intent(in) :: lo,hi,padding
  type(rtensorfield_t),intent(out) :: df
  double precision,dimension(:),pointer :: coeff
  integer :: i,j,n,m,l
  integer :: i0,j0,n0,m0,l0
  integer :: order,bias

  order = ischeme/10
  bias  = mod(ischeme,10)
  call select_coeff(coeff,order,bias)
  
  select case(axis)
    case(axis_xi)
      do j = lo%d(5), hi%d(5)
        do i = lo%d(4), hi%d(4)
          do l = lo%d(3), hi%d(3)
            do m = lo%d(2), hi%d(2)
              do n = lo%d(1), hi%d(1)
                n0 = n + padding%d(1)
                m0 = m + padding%d(2)
                l0 = l + padding%d(3)
                i0 = i + padding%d(4)
                j0 = j + padding%d(5)
                df%p(n,m,l,i,j) = dot_product(f%p(n0,m0,l0,i0-bias:i0+order-bias,j0),coeff(:))            
              enddo
            enddo
          enddo
        enddo
      enddo
    case(axis_eta)
      do j = lo%d(5), hi%d(5)
        do i = lo%d(4), hi%d(4)
          do l = lo%d(3), hi%d(3)
            do m = lo%d(2), hi%d(2)
              do n = lo%d(1), hi%d(1)
                n0 = n + padding%d(1)
                m0 = m + padding%d(2)
                l0 = l + padding%d(3)
                i0 = i + padding%d(4)
                j0 = j + padding%d(5)
                df%p(n,m,l,i,j) = dot_product(f%p(n0,m0,l0,i0,j0-bias:j0+order-bias),coeff(:))            
              enddo
            enddo
          enddo
        enddo
      enddo
    case default
      write(*,*) "Error: compute_grad_homo_tensor():"
      write(*,*) "Illegal axis:",axis
      stop
  end select
end subroutine compute_grad_homo_tensor

subroutine compute_grad_heter_scalar(df,f,lo,hi,padding,axis,flap,order)
  implicit none
  integer,intent(in) :: axis,flap,order
  type(rscalarfield_t),intent(in) :: f
  type(fscalarshape_t),intent(in) :: lo,hi,padding
  type(rscalarfield_t),intent(out) :: df

  double precision,dimension(:),pointer :: coeff
  integer :: i,j,pnt,bias
  integer :: i0,j0
  integer :: st,ed,pad,wing
  wing = mod(order,2) + order/2 
  select case(axis)
    case(axis_xi)
      if(flap.eq.0) then
        st = lo%d(1)
        ed = lo%d(1) + wing - 1
      else if(flap.eq.1) then
        st = hi%d(1) - wing + 1
        ed = hi%d(1) 
      else
        write(*,*) "Error: compute_grad_heter_scalar():"
        write(*,*) "Illegal flap:", flap
      endif
      pad = flap*wing - st
      do j = lo%d(2), hi%d(2)
        do i = st, ed
          bias = pad + i 
          i0 = i + padding%d(1) - bias
          j0 = j + padding%d(2)
          call select_coeff(coeff,order,bias)
          df%p(i,j) = dot_product(f%p(i0:i0+order,j0),coeff(:))
        enddo
      enddo
    case(axis_eta)
      if(flap.eq.0) then
        st = lo%d(2)
        ed = lo%d(2) + wing - 1
      else if(flap.eq.1) then
        st = hi%d(2) - wing + 1
        ed = hi%d(2) 
      else
        write(*,*) "Error: compute_grad_heter_scalar():"
        write(*,*) "Illegal flap:", flap
      endif
      pad = flap*wing - st
      do j = st, ed
        do i = lo%d(1), hi%d(1)
          bias = pad + j
          i0 = i + padding%d(1)
          j0 = j + padding%d(2) - bias
          call select_coeff(coeff,order,bias)
          df%p(i,j) = dot_product(f%p(i0,j0:j0+order),coeff(:))
        enddo
      enddo
    case default
      write(*,*) "Error: compute_grad_heter_scalar():"
      write(*,*) "Illegal axis:",axis
      stop
  end select
end subroutine compute_grad_heter_scalar

subroutine compute_grad_heter_vector(df,f,lo,hi,padding,axis,flap,order)
  implicit none
  integer,intent(in) :: axis,flap,order
  type(rvectorfield_t),intent(in) :: f
  type(fvectorshape_t),intent(in) :: lo,hi,padding
  type(rvectorfield_t),intent(out) :: df

  double precision,dimension(:),pointer :: coeff
  integer :: i,j,n
  integer :: bias
  integer :: i0,j0,n0
  integer :: st,ed,pad,wing
  wing = mod(order,2) + order/2 
  select case(axis)
    case(axis_xi)
      if(flap.eq.0) then
        st = lo%d(2)
        ed = lo%d(2) + wing - 1
      else if(flap.eq.1) then
        st = hi%d(2) - wing + 1
        ed = hi%d(2) 
      else
        write(*,*) "Error: compute_grad_heter_vector():"
        write(*,*) "Illegal flap:", flap
      endif
      pad = flap*wing - st
      do j = lo%d(3), hi%d(3)
        do i = st, ed
          bias = pad + i 
          i0 = i + padding%d(2) - bias
          j0 = j + padding%d(3)
          call select_coeff(coeff,order,bias)
          do n = lo%d(1), hi%d(1)
            n0 = n + padding%d(1)
            df%p(n,i,j) = dot_product(f%p(n0,i0:i0+order,j0),coeff(:))
          enddo
        enddo
      enddo
    case(axis_eta)
      if(flap.eq.0) then
        st = lo%d(3)
        ed = lo%d(3) + wing - 1
      else if(flap.eq.1) then
        st = hi%d(3) - wing + 1
        ed = hi%d(3) 
      else
        write(*,*) "Error: compute_grad_heter_vector():"
        write(*,*) "Illegal flap:", flap
      endif
      pad = flap*wing - st
      do j = st, ed
        do i = lo%d(2), hi%d(2)
          bias = pad + j
          i0 = i + padding%d(2)
          j0 = j + padding%d(3) - bias
          call select_coeff(coeff,order,bias)
          do n = lo%d(1), hi%d(1)
            n0 = n + padding%d(1)
            df%p(n,i,j) = dot_product(f%p(n0,i0,j0:j0+order),coeff(:))
          enddo
        enddo
      enddo
    case default
      write(*,*) "Error: compute_grad_heter_vector():"
      write(*,*) "Illegal axis:",axis
      stop
  end select
end subroutine compute_grad_heter_vector

subroutine compute_grad_heter_matrix(df,f,lo,hi,padding,axis,flap,order)
  implicit none
  integer,intent(in) :: axis,flap,order
  type(rmatrixfield_t),intent(in) :: f
  type(fmatrixshape_t),intent(in) :: lo,hi,padding
  type(rmatrixfield_t),intent(out) :: df

  double precision,dimension(:),pointer :: coeff
  integer :: i,j,n,m
  integer :: bias
  integer :: i0,j0,n0,m0
  integer :: st,ed,pad,wing
  wing = mod(order,2) + order/2 
  select case(axis)
    case(axis_xi)
      if(flap.eq.0) then
        st = lo%d(3)
        ed = lo%d(3) + wing - 1
      else if(flap.eq.1) then
        st = hi%d(3) - wing + 1
        ed = hi%d(3) 
      else
        write(*,*) "Error: compute_grad_heter_matrix():"
        write(*,*) "Illegal flap:", flap
      endif
      pad = flap*wing - st
      do j = lo%d(4), hi%d(4)
        do i = st, ed
          bias = pad + i 
          i0 = i + padding%d(3) - bias
          j0 = j + padding%d(4)
          call select_coeff(coeff,order,bias)
          do m = lo%d(2), hi%d(2)
            do n = lo%d(1), hi%d(1)
              n0 = n + padding%d(1)
              m0 = m + padding%d(2)
              df%p(n,m,i,j) = dot_product(f%p(n0,m0,i0:i0+order,j0),coeff(:))
            enddo
          enddo
        enddo
      enddo
    case(axis_eta)
      if(flap.eq.0) then
        st = lo%d(4)
        ed = lo%d(4) + wing - 1
      else if(flap.eq.1) then
        st = hi%d(4) - wing + 1
        ed = hi%d(4) 
      else
        write(*,*) "Error: compute_grad_heter_matrix():"
        write(*,*) "Illegal flap:", flap
      endif
      pad = flap*wing - st
      do j = st, ed
        do i = lo%d(3), hi%d(3)
          bias = pad + j
          i0 = i + padding%d(3)
          j0 = j + padding%d(4) - bias
          call select_coeff(coeff,order,bias)
          do m = lo%d(2), hi%d(2)
            do n = lo%d(1), hi%d(1)
              n0 = n + padding%d(1)
              m0 = m + padding%d(2)
              df%p(n,m,i,j) = dot_product(f%p(n0,m0,i0,j0:j0+order),coeff(:))
            enddo
          enddo
        enddo
      enddo
    case default
      write(*,*) "Error: compute_grad_heter_matrix():"
      write(*,*) "Illegal axis:",axis
      stop
  end select
end subroutine compute_grad_heter_matrix

subroutine compute_grad_heter_tensor(df,f,lo,hi,padding,axis,flap,order)
  implicit none
  integer,intent(in) :: axis,flap,order
  type(rtensorfield_t),intent(in) :: f
  type(ftensorshape_t),intent(in) :: lo,hi,padding
  type(rtensorfield_t),intent(out) :: df

  double precision,dimension(:),pointer :: coeff
  integer :: i,j,n,m,l
  integer :: bias
  integer :: i0,j0,n0,m0,l0
  integer :: st,ed,pad,wing
  wing = mod(order,2) + order/2 
  select case(axis)
    case(axis_xi)
      if(flap.eq.0) then
        st = lo%d(4)
        ed = lo%d(4) + wing - 1
      else if(flap.eq.1) then
        st = hi%d(4) - wing + 1
        ed = hi%d(4) 
      else
        write(*,*) "Error: compute_grad_heter_tensor():"
        write(*,*) "Illegal flap:", flap
      endif
      pad = flap*wing - st
      do j = lo%d(5), hi%d(5)
        do i = st, ed
          bias = pad + i 
          i0 = i + padding%d(4) - bias
          j0 = j + padding%d(5)
          call select_coeff(coeff,order,bias)
          do l = lo%d(3), hi%d(3)
            do m = lo%d(2), hi%d(2)
              do n = lo%d(1), hi%d(1)
                n0 = n + padding%d(1)
                m0 = m + padding%d(2)
                l0 = l + padding%d(3)
                df%p(n,m,l,i,j) = dot_product(f%p(n0,m0,l0,i0:i0+order,j0),coeff(:))
              enddo
            enddo
          enddo
        enddo
      enddo
    case(axis_eta)
      if(flap.eq.0) then
        st = lo%d(5)
        ed = lo%d(5) + wing - 1
      else if(flap.eq.1) then
        st = hi%d(5) - wing + 1
        ed = hi%d(5) 
      else
        write(*,*) "Error: compute_grad_heter_tensor():"
        write(*,*) "Illegal flap:", flap
      endif
      pad = flap*wing - st
      do j = st, ed
        do i = lo%d(4), hi%d(4)
          bias = pad + j
          i0 = i + padding%d(4)
          j0 = j + padding%d(5) - bias
          call select_coeff(coeff,order,bias)
          do l = lo%d(3), hi%d(3)
            do m = lo%d(2), hi%d(2)
              do n = lo%d(1), hi%d(1)
                n0 = n + padding%d(1)
                m0 = m + padding%d(2)
                l0 = l + padding%d(3)
                df%p(n,m,l,i,j) = dot_product(f%p(n0,m0,l0,i0,j0:j0+order),coeff(:))
              enddo
            enddo
          enddo
        enddo
      enddo
    case default
      write(*,*) "Error: compute_grad_heter_tensor():"
      write(*,*) "Illegal axis:",axis
      stop
  end select
end subroutine compute_grad_heter_tensor

subroutine select_coeff(coeff,order,bias)
  implicit none
  integer,intent(in) :: order,bias
  double precision,dimension(:),pointer,intent(inout) :: coeff
  select case(order)
    case (1)
      call select_coeff_1st(coeff,bias)
    case (2)
      call select_coeff_2nd(coeff,bias)
    case (3)
      call select_coeff_3rd(coeff,bias)
    case (4)
      call select_coeff_4th(coeff,bias)
    case (5)
      call select_coeff_5th(coeff,bias)
    case (6)
      call select_coeff_6th(coeff,bias)
    case (7)
      call select_coeff_7th(coeff,bias)
    case (8)
      call select_coeff_8th(coeff,bias)
    case default
      write(*,*) "Error: compute_grad_cd():"
      write(*,*) "No such order",order    
      stop
  end select
end subroutine select_coeff

subroutine select_coeff_1st(coeff,bias)
  implicit none
  integer,intent(in) :: bias
  double precision,dimension(:),pointer,intent(inout) :: coeff
  select case(bias)
  case (1)
    coeff => a1_10
  case (0)
    coeff => a1_01
  case default
    write(*,*) "Error: select_coeff_1st():"
    write(*,*) "No such bias",bias    
  end select
end subroutine select_coeff_1st

subroutine select_coeff_2nd(coeff,bias)
  implicit none
  integer,intent(in) :: bias
  double precision,dimension(:),pointer,intent(inout) :: coeff
  select case(bias)
  case (2)
    coeff => a2_20
  case (1)
    coeff => a2_11
  case (0)
    coeff => a2_02
  case default
    write(*,*) "Error: select_coeff_2nd():"
    write(*,*) "No such bias",bias    
  end select
end subroutine select_coeff_2nd

subroutine select_coeff_3rd(coeff,bias)
  implicit none
  integer,intent(in) :: bias
  double precision,dimension(:),pointer,intent(inout) :: coeff
  select case(bias)
  case (3)
    coeff => a3_30
  case (2)
    coeff => a3_21
  case (1)
    coeff => a3_12
  case (0)
    coeff => a3_03
  case default
    write(*,*) "Error: select_coeff_3rd():"
    write(*,*) "No such bias",bias    
  end select
end subroutine select_coeff_3rd

subroutine select_coeff_4th(coeff,bias)
  implicit none
  integer,intent(in) :: bias
  double precision,dimension(:),pointer,intent(inout) :: coeff
  select case(bias)
  case (4)
    coeff => a4_40
  case (3)
    coeff => a4_31
  case (2)
    coeff => a4_22
  case (1)
    coeff => a4_13
  case (0)
    coeff => a4_04
  case default
    write(*,*) "Error: select_coeff_4th():"
    write(*,*) "No such bias",bias    
  end select
end subroutine select_coeff_4th

subroutine select_coeff_5th(coeff,bias)
  implicit none
  integer,intent(in) :: bias
  double precision,dimension(:),pointer,intent(inout) :: coeff
  select case(bias)
  case (5)
    coeff => a5_50
  case (4)
    coeff => a5_41
  case (3)
    coeff => a5_32
  case (2)
    coeff => a5_23
  case (1)
    coeff => a5_14
  case (0)
    coeff => a5_05
  case default
    write(*,*) "Error: select_coeff_5th():"
    write(*,*) "No such bias",bias    
  end select
end subroutine select_coeff_5th

subroutine select_coeff_6th(coeff,bias)
  implicit none
  integer,intent(in) :: bias
  double precision,dimension(:),pointer,intent(inout) :: coeff
  select case(bias)
  case (6)
    coeff => a6_60
  case (5)
    coeff => a6_51
  case (4)
    coeff => a6_42
  case (3)
    coeff => a6_33
  case (2)
    coeff => a6_24
  case (1)
    coeff => a6_15
  case (0)
    coeff => a6_06
  case default
    write(*,*) "Error: select_coeff_6th():"
    write(*,*) "No such bias",bias    
  end select
end subroutine select_coeff_6th

subroutine select_coeff_7th(coeff,bias)
  implicit none
  integer,intent(in) :: bias
  double precision,dimension(:),pointer,intent(inout) :: coeff
  select case(bias)
  case (7)
    coeff => a7_70
  case (6)
    coeff => a7_61
  case (5)
    coeff => a7_52
  case (4)
    coeff => a7_43
  case (3)
    coeff => a7_34
  case (2)
    coeff => a7_25
  case (1)
    coeff => a7_16
  case (0)
    coeff => a7_07
  case default
    write(*,*) "Error: select_coeff_7th():"
    write(*,*) "No such bias",bias    
  end select
end subroutine select_coeff_7th

subroutine select_coeff_8th(coeff,bias)
  implicit none
  integer,intent(in) :: bias
  double precision,dimension(:),pointer,intent(inout) :: coeff
  select case(bias)
  case (8)
    coeff => a8_80
  case (7)
    coeff => a8_71
  case (6)
    coeff => a8_62
  case (5)
    coeff => a8_53
  case (4)
    coeff => a8_44
  case (3)
    coeff => a8_35
  case (2)
    coeff => a8_26
  case (1)
    coeff => a8_17
  case (0)
    coeff => a8_08
  case default
    write(*,*) "Error: select_coeff_8th():"
    write(*,*) "No such bias",bias    
  end select
end subroutine select_coeff_8th

end module gradient