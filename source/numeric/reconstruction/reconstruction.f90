! Reconstruction is alway performed to the first index
module reconstruction
  use data_type
  implicit none

  ! Linear upwind schemes
  integer,parameter :: upwind1_id  = 101
  integer,parameter :: upwind3_id  = 102
  integer,parameter :: upwind5_id  = 103
  integer,parameter :: upwind7_id  = 104

  ! Linear Central schemes
  integer,parameter :: central2_id = 201
  integer,parameter :: central4_id = 202
  integer,parameter :: central6_id = 203
  integer,parameter :: central8_id = 204

  ! WENO schemes
  integer,parameter :: wenojs3_id    = 301
  integer,parameter :: wenojs5_id    = 302
  integer,parameter :: wenojs7_id    = 303
  integer,parameter :: wenoz3_id     = 304
  integer,parameter :: wenoz5_id     = 305
  integer,parameter :: wenoz7_id     = 306
  integer,parameter :: AdaWENOZ5_id  = 307
  integer,parameter :: wenoe5_id     = 308
  
  integer,parameter :: comm_wenojs5_id    = 401
  integer,parameter :: comm_wenoz5_id     = 402
  integer,parameter :: comm_wenoe5_id     = 403

  interface reconst
    module procedure reconst_component_scalar
    module procedure reconst_component_vector
    module procedure reconst_component_matrix
    module procedure reconst_component_tensor
    module procedure reconst_commweight_scalar
    module procedure reconst_commweight_vector
    module procedure reconst_commweight_matrix
    module procedure reconst_commweight_tensor
    ! module procedure reconst_character
    ! module procedure reconst_adaptive
  end interface reconst

  interface
    module subroutine reconst_component_scalar_uw1(h,f,lo,hi,padding,bias,axis)
      implicit none
      type(rscalarfield_t),intent(inout)     :: h
      type(rscalarfield_t),target,intent(in) :: f
      type(fscalarshape_t),intent(in)        :: lo,hi,padding
      integer,intent(in)                     :: bias,axis
    end subroutine

    module subroutine reconst_component_vector_uw1(h,f,lo,hi,padding,bias,axis)
      implicit none
      type(rvectorfield_t),intent(inout)     :: h
      type(rvectorfield_t),target,intent(in) :: f
      type(fvectorshape_t),intent(in)        :: lo,hi,padding
      integer,intent(in)                     :: bias,axis
    end subroutine


    module subroutine reconst_component_matrix_uw1(h,f,lo,hi,padding,bias,axis)
      implicit none
      type(rmatrixfield_t),intent(inout)     :: h
      type(rmatrixfield_t),target,intent(in) :: f
      type(fmatrixshape_t),intent(in)        :: lo,hi,padding
      integer,intent(in)                     :: bias,axis
    end subroutine

    module subroutine reconst_component_tensor_uw1(h,f,lo,hi,padding,bias,axis)
      implicit none
      type(rtensorfield_t),intent(inout)     :: h
      type(rtensorfield_t),target,intent(in) :: f
      type(ftensorshape_t),intent(in)        :: lo,hi,padding
      integer,intent(in)                     :: bias,axis
    end subroutine

    module subroutine reconst_component_scalar_wenojs3(h,f,lo,hi,padding,bias,axis)
        implicit none
        type(rscalarfield_t),intent(inout)     :: h
        type(rscalarfield_t),target,intent(in) :: f
        type(fscalarshape_t),intent(in)        :: lo,hi,padding
        integer,intent(in)                     :: bias,axis
    end subroutine

    module subroutine reconst_component_vector_wenojs3(h,f,lo,hi,padding,bias,axis)
      implicit none
      type(rvectorfield_t),intent(inout)     :: h
      type(rvectorfield_t),target,intent(in) :: f
      type(fvectorshape_t),intent(in)        :: lo,hi,padding
      integer,intent(in)                     :: bias,axis
    end subroutine

    module subroutine reconst_component_matrix_wenojs3(h,f,lo,hi,padding,bias,axis)
      implicit none
      type(rmatrixfield_t),intent(inout)     :: h
      type(rmatrixfield_t),target,intent(in) :: f
      type(fmatrixshape_t),intent(in)        :: lo,hi,padding
      integer,intent(in)                     :: bias,axis
    end subroutine

    module subroutine reconst_component_tensor_wenojs3(h,f,lo,hi,padding,bias,axis)
      implicit none
      type(rtensorfield_t),intent(inout)     :: h
      type(rtensorfield_t),target,intent(in) :: f
      type(ftensorshape_t),intent(in)        :: lo,hi,padding
      integer,intent(in)                     :: bias,axis
    end subroutine

    module subroutine reconst_component_scalar_wenojs5(h,f,lo,hi,padding,bias,axis)
      implicit none
      type(rscalarfield_t),intent(inout)     :: h
      type(rscalarfield_t),target,intent(in) :: f
      type(fscalarshape_t),intent(in)        :: lo,hi,padding
      integer,intent(in)                     :: bias,axis
    end subroutine

    module subroutine reconst_component_vector_wenojs5(h,f,lo,hi,padding,bias,axis)
      implicit none
      type(rvectorfield_t),intent(inout)     :: h
      type(rvectorfield_t),target,intent(in) :: f
      type(fvectorshape_t),intent(in)        :: lo,hi,padding
      integer,intent(in)                     :: bias,axis
    end subroutine
    
    module subroutine reconst_component_matrix_wenojs5(h,f,lo,hi,padding,bias,axis)
      implicit none
      type(rmatrixfield_t),intent(inout)     :: h
      type(rmatrixfield_t),target,intent(in) :: f
      type(fmatrixshape_t),intent(in)        :: lo,hi,padding
      integer,intent(in)                     :: bias,axis
    end subroutine

    module subroutine reconst_component_tensor_wenojs5(h,f,lo,hi,padding,bias,axis)
      implicit none
      type(rtensorfield_t),intent(inout)     :: h
      type(rtensorfield_t),target,intent(in) :: f
      type(ftensorshape_t),intent(in)        :: lo,hi,padding
      integer,intent(in)                     :: bias,axis
    end subroutine

    module subroutine reconst_component_scalar_wenoz5(h,f,lo,hi,padding,bias,axis)
      implicit none
      type(rscalarfield_t),intent(inout)     :: h
      type(rscalarfield_t),target,intent(in) :: f
      type(fscalarshape_t),intent(in)        :: lo,hi,padding
      integer,intent(in)                     :: bias,axis
    end subroutine

    module subroutine reconst_component_vector_wenoz5(h,f,lo,hi,padding,bias,axis)
      implicit none
      type(rvectorfield_t),intent(inout)     :: h
      type(rvectorfield_t),target,intent(in) :: f
      type(fvectorshape_t),intent(in)        :: lo,hi,padding
      integer,intent(in)                     :: bias,axis
    end subroutine

    module subroutine reconst_component_matrix_wenoz5(h,f,lo,hi,padding,bias,axis)
      implicit none
      type(rmatrixfield_t),intent(inout)     :: h
      type(rmatrixfield_t),target,intent(in) :: f
      type(fmatrixshape_t),intent(in)        :: lo,hi,padding
      integer,intent(in)                     :: bias,axis
    end subroutine

    module subroutine reconst_component_tensor_wenoz5(h,f,lo,hi,padding,bias,axis)
      implicit none
      type(rtensorfield_t),intent(inout)     :: h
      type(rtensorfield_t),target,intent(in) :: f
      type(ftensorshape_t),intent(in)        :: lo,hi,padding
      integer,intent(in)                     :: bias,axis
    end subroutine
   
    module subroutine reconst_component_scalar_wenoe5(h,f,lo,hi,padding,bias,axis)
      implicit none
      type(rscalarfield_t),intent(inout)     :: h
      type(rscalarfield_t),target,intent(in) :: f
      type(fscalarshape_t),intent(in)        :: lo,hi,padding
      integer,intent(in)                     :: bias,axis
    end subroutine

    module subroutine reconst_component_vector_wenoe5(h,f,lo,hi,padding,bias,axis)
      implicit none
      type(rvectorfield_t),intent(inout)     :: h
      type(rvectorfield_t),target,intent(in) :: f
      type(fvectorshape_t),intent(in)        :: lo,hi,padding
      integer,intent(in)                     :: bias,axis
    end subroutine

    module subroutine reconst_component_matrix_wenoe5(h,f,lo,hi,padding,bias,axis)
      implicit none
      type(rmatrixfield_t),intent(inout)     :: h
      type(rmatrixfield_t),target,intent(in) :: f
      type(fmatrixshape_t),intent(in)        :: lo,hi,padding
      integer,intent(in)                     :: bias,axis
    end subroutine

    module subroutine reconst_component_tensor_wenoe5(h,f,lo,hi,padding,bias,axis)
      implicit none
      type(rtensorfield_t),intent(inout)     :: h
      type(rtensorfield_t),target,intent(in) :: f
      type(ftensorshape_t),intent(in)        :: lo,hi,padding
      integer,intent(in)                     :: bias,axis
    end subroutine

    module subroutine cmpt_commweight_wenojs5(omega,f,lo,hi,padding,bias,axis)
      implicit none
      type(rvectorfield_t),intent(inout)     :: omega
      type(rscalarfield_t),target,intent(in) :: f
      type(fscalarshape_t),intent(in)        :: lo,hi,padding
      integer,intent(in)                     :: bias,axis
    end subroutine

    module subroutine cmpt_commweight_wenoz5(omega,f,lo,hi,padding,bias,axis)
      implicit none
      type(rvectorfield_t),intent(inout)     :: omega
      type(rscalarfield_t),target,intent(in) :: f
      type(fscalarshape_t),intent(in)        :: lo,hi,padding
      integer,intent(in)                     :: bias,axis
    end subroutine

    module subroutine cmpt_commweight_wenoe5(omega,f,lo,hi,padding,bias,axis)
      implicit none
      type(rvectorfield_t),intent(inout)     :: omega
      type(rscalarfield_t),target,intent(in) :: f
      type(fscalarshape_t),intent(in)        :: lo,hi,padding
      integer,intent(in)                     :: bias,axis
    end subroutine

    module subroutine reconst_commweight_scalar_weno5(h,f,omega,lo,hi,padding,bias,axis)
      implicit none
      type(rscalarfield_t),intent(inout)     :: h
      type(rscalarfield_t),target,intent(in) :: f
      type(rvectorfield_t),intent(in)        :: omega
      type(fscalarshape_t),intent(in)        :: lo,hi,padding
      integer,intent(in)                     :: bias,axis
    
    end subroutine

    module subroutine reconst_commweight_vector_weno5(h,f,omega,lo,hi,padding,bias,axis)
      implicit none
      type(rvectorfield_t),intent(inout)     :: h
      type(rvectorfield_t),target,intent(in) :: f
      type(rvectorfield_t),intent(in)        :: omega
      type(fvectorshape_t),intent(in)        :: lo,hi,padding
      integer,intent(in)                     :: bias,axis
    
    end subroutine
    
    module subroutine reconst_commweight_matrix_weno5(h,f,omega,lo,hi,padding,bias,axis)
      implicit none
      type(rmatrixfield_t),intent(inout)     :: h
      type(rmatrixfield_t),target,intent(in) :: f
      type(rvectorfield_t),intent(in)        :: omega
      type(fmatrixshape_t),intent(in)        :: lo,hi,padding
      integer,intent(in)                     :: bias,axis
    end subroutine

    module subroutine reconst_commweight_tensor_weno5(h,f,omega,lo,hi,padding,bias,axis)
      implicit none
      type(rtensorfield_t),intent(inout)     :: h
      type(rtensorfield_t),target,intent(in) :: f
      type(rvectorfield_t),intent(in)        :: omega
      type(ftensorshape_t),intent(in)        :: lo,hi,padding
      integer,intent(in)                     :: bias,axis      
    end subroutine

  end interface
contains
  subroutine reconst_component_scalar(h,f,lo,hi,padding,bias,axis,ireconst)
    implicit none
    type(rscalarfield_t),intent(inout) :: h
    type(rscalarfield_t),intent(in) :: f
    type(fscalarshape_t),intent(in) :: lo,hi,padding
    integer,intent(in) :: bias,axis
    integer,intent(in) :: ireconst
    select case (ireconst)
    case(upwind1_id)
      call reconst_component_scalar_uw1(h,f,lo,hi,padding,bias,axis)
    case(upwind3_id)
      ! call reconst_component_scalar_uw3(h,f,lo,hi,padding,bias,axis)
    case(upwind5_id)
      ! call reconst_component_scalar_uw5(h,f,lo,hi,padding,bias,axis)
    case(upwind7_id)
      ! call reconst_component_scalar_uw7(h,f,lo,hi,padding,bias,axis)
    case(central2_id)
      ! call reconst_component_scalar_cn2(h,f,lo,hi,padding,bias,axis)
    case(central4_id)
      ! call reconst_component_scalar_cn4(h,f,lo,hi,padding,bias,axis)
    case(central6_id)
      ! call reconst_component_scalar_cn6(h,f,lo,hi,padding,bias,axis)
    case(central8_id)
      ! call reconst_component_scalar_cn8(h,f,lo,hi,padding,bias,axis)
    case(wenojs3_id)
      call reconst_component_scalar_wenojs3(h,f,lo,hi,padding,bias,axis)
    case(wenojs5_id)
      call reconst_component_scalar_wenojs5(h,f,lo,hi,padding,bias,axis)
    case(wenojs7_id)
      ! call reconst_component_scalar_wenojs7(h,f,lo,hi,padding,bias,axis)
    case(wenoz3_id)
      ! call reconst_component_scalar_wenoz3(h,f,lo,hi,padding,bias,axis)
    case(wenoz5_id)
      call reconst_component_scalar_wenoz5(h,f,lo,hi,padding,bias,axis)
    case(wenoe5_id)
      call reconst_component_scalar_wenoe5(h,f,lo,hi,padding,bias,axis)
    case(wenoz7_id)
      ! call reconst_component_scalar_wenoz7(h,f,lo,hi,padding,bias,axis)
    case default
      print*,"No such ireconst: ",ireconst
    end select
  end subroutine reconst_component_scalar

  subroutine reconst_component_vector(h,f,lo,hi,padding,bias,axis,ireconst)
    implicit none
    type(rvectorfield_t),intent(inout) :: h
    type(rvectorfield_t),intent(in) :: f
    type(fvectorshape_t),intent(in) :: lo,hi,padding
    integer,intent(in) :: bias,axis
    integer,intent(in) :: ireconst
    select case (ireconst)
    case(upwind1_id)
      call reconst_component_vector_uw1(h,f,lo,hi,padding,bias,axis)
    case(wenojs3_id)
      call reconst_component_vector_wenojs3(h,f,lo,hi,padding,bias,axis)
    case(wenojs5_id)
      call reconst_component_vector_wenojs5(h,f,lo,hi,padding,bias,axis)
    case(wenoz5_id)
      call reconst_component_vector_wenoz5(h,f,lo,hi,padding,bias,axis)
    case(wenoe5_id)
      call reconst_component_vector_wenoe5(h,f,lo,hi,padding,bias,axis)
    case(wenoz7_id)
      ! call reconst_component_vector_wenoz7(h,f,lo,hi,padding,bias,axis)
    case(wenojs7_id)
      ! call reconst_component_vector_wenojs7(h,f,lo,hi,padding,bias,axis)
    case default
      print*,"No such ireconst: ",ireconst
    end select
  end subroutine reconst_component_vector

  subroutine reconst_component_matrix(h,f,lo,hi,padding,bias,axis,ireconst)
    implicit none
    type(rmatrixfield_t),intent(inout) :: h
    type(rmatrixfield_t),intent(in) :: f
    type(fmatrixshape_t),intent(in) :: lo,hi,padding
    integer,intent(in) :: bias,axis
    integer,intent(in) :: ireconst
    select case (ireconst)
    case(upwind1_id)
      call reconst_component_matrix_uw1(h,f,lo,hi,padding,bias,axis)
    case(wenojs3_id)
      call reconst_component_matrix_wenojs3(h,f,lo,hi,padding,bias,axis)
    case(wenojs5_id)
      call reconst_component_matrix_wenojs5(h,f,lo,hi,padding,bias,axis)
    case(wenoz5_id)
      call reconst_component_matrix_wenoz5(h,f,lo,hi,padding,bias,axis)
    case(wenoe5_id)
      call reconst_component_matrix_wenoe5(h,f,lo,hi,padding,bias,axis)
    case default
      print*,"No such ireconst: ",ireconst
    end select
  end subroutine reconst_component_matrix

  subroutine reconst_component_tensor(h,f,lo,hi,padding,bias,axis,ireconst)
    implicit none
    type(rtensorfield_t),intent(inout) :: h
    type(rtensorfield_t),intent(in) :: f
    type(ftensorshape_t),intent(in) :: lo,hi,padding
    integer,intent(in) :: bias,axis
    integer,intent(in) :: ireconst
    select case (ireconst)
    case(upwind1_id)
      call reconst_component_tensor_uw1(h,f,lo,hi,padding,bias,axis)
    case(wenojs3_id)
      call reconst_component_tensor_wenojs3(h,f,lo,hi,padding,bias,axis)
    case(wenojs5_id)
      call reconst_component_tensor_wenojs5(h,f,lo,hi,padding,bias,axis)
    case(wenoz5_id)
      call reconst_component_tensor_wenoz5(h,f,lo,hi,padding,bias,axis)
    case(wenoe5_id)
      call reconst_component_tensor_wenoe5(h,f,lo,hi,padding,bias,axis)
    case default
      print*,"No such ireconst: ",ireconst
    end select
  end subroutine reconst_component_tensor

  ! subroutine reconst_character(h,f,l,r,lo,hi,padding,bias,ireconst)
  !   implicit none
  !   type(rvectorfield_t),intent(inout) :: h
  !   type(rvectorfield_t),intent(in) :: f
  !   type(rmatrixfield_t),intent(in) :: l,r
  !   type(fvectorshape_t),intent(in) :: lo,hi,padding
  !   integer,intent(in) :: bias
  !   integer,intent(in) :: ireconst
  ! end subroutine reconst_character

  ! subroutine reconst_adaptive(h,f,l,r,flag,lo,hi,padding,bias,ireconst)
  !   implicit none
  !   type(rvectorfield_t),intent(inout) :: h
  !   type(rvectorfield_t),intent(in) :: f
  !   type(rmatrixfield_t),intent(in) :: l,r
  !   type(rvectorfield_t),intent(in) :: flag
  !   type(fvectorshape_t),intent(in) :: lo,hi,padding
  !   integer,intent(in) :: bias
  !   integer,intent(in) :: ireconst
  ! end subroutine reconst_adaptive

  subroutine reconst_commweight_scalar(h,f,measure,lo,hi,padding,bias,axis,ireconst)
    implicit none
    type(rscalarfield_t),intent(inout) :: h
    type(rscalarfield_t),intent(inout) :: f
    type(rscalarfield_t),intent(in)    :: measure
    type(fscalarshape_t),intent(in)    :: lo,hi,padding
    integer,intent(in)                 :: bias,axis
    integer,intent(in)                 :: ireconst
    
    type(rvectorfield_t) :: omega
    type(fvectorshape_t) :: omega_s
    type(fscalarshape_t) :: omega_lo,omega_hi,omega_padding
    
    omega_lo      = set_shape(lo%d(1),lo%d(2))
    omega_hi      = set_shape(hi%d(1),hi%d(2))
    omega_padding = set_shape(padding%d(1),padding%d(2))

    omega_s       = set_shape(3,h%s%d(1),h%s%d(2))
    call create(omega,omega_s)
    select case (ireconst)
    case(wenojs5_id)
      call cmpt_commweight_wenojs5(omega,measure,omega_lo,omega_hi,omega_padding,bias,axis)
    case(wenoz5_id)
      call cmpt_commweight_wenoz5(omega,measure,omega_lo,omega_hi,omega_padding,bias,axis)
    case(wenoe5_id)
      call cmpt_commweight_wenoe5(omega,measure,omega_lo,omega_hi,omega_padding,bias,axis)
    case default
      print*,"Error: reconst_commweight()"
      print*,"No such ireconst: ",ireconst
      stop
    end select
    call reconst_commweight_scalar_weno5(h,f,omega,lo,hi,padding,bias,axis)
    call destroy(omega)
  end subroutine reconst_commweight_scalar

  subroutine reconst_commweight_vector(h,f,measure,lo,hi,padding,bias,axis,ireconst)
    implicit none
    type(rvectorfield_t),intent(inout) :: h
    type(rvectorfield_t),intent(inout) :: f
    type(rscalarfield_t),intent(in)    :: measure
    type(fvectorshape_t),intent(in)    :: lo,hi,padding
    integer,intent(in)                 :: bias,axis
    integer,intent(in)                 :: ireconst
    
    type(rvectorfield_t) :: omega
    type(fvectorshape_t) :: omega_s
    type(fscalarshape_t) :: omega_lo,omega_hi,omega_padding
    
    omega_lo      = set_shape(lo%d(2),lo%d(3))
    omega_hi      = set_shape(hi%d(2),hi%d(3))
    omega_padding = set_shape(padding%d(2),padding%d(3))
    
    omega_s = set_shape(3,h%s%d(2),h%s%d(3))
    call create(omega,omega_s)
    select case (ireconst)
    case(wenojs5_id)
      call cmpt_commweight_wenojs5(omega,measure,omega_lo,omega_hi,omega_padding,bias,axis)
    case(wenoz5_id)
      call cmpt_commweight_wenoz5(omega,measure,omega_lo,omega_hi,omega_padding,bias,axis)
    case(wenoe5_id)
      call cmpt_commweight_wenoe5(omega,measure,omega_lo,omega_hi,omega_padding,bias,axis)
    case default
      print*,"Error: reconst_commweight()"
      print*,"No such ireconst: ",ireconst
      stop
    end select
    call reconst_commweight_vector_weno5(h,f,omega,lo,hi,padding,bias,axis)
    call destroy(omega)
  end subroutine reconst_commweight_vector

  subroutine reconst_commweight_matrix(h,f,measure,lo,hi,padding,bias,axis,ireconst)
    implicit none
    type(rmatrixfield_t),intent(inout) :: h
    type(rmatrixfield_t),intent(inout) :: f
    type(rscalarfield_t),intent(in)    :: measure
    type(fmatrixshape_t),intent(in)    :: lo,hi,padding
    integer,intent(in)                 :: bias,axis
    integer,intent(in)                 :: ireconst
    
    type(rvectorfield_t) :: omega
    type(fvectorshape_t) :: omega_s
    type(fscalarshape_t) :: omega_lo,omega_hi,omega_padding
    
    omega_lo      = set_shape(lo%d(3),lo%d(4))
    omega_hi      = set_shape(hi%d(3),hi%d(4))
    omega_padding = set_shape(padding%d(3),padding%d(4))
    
    omega_s = set_shape(3,h%s%d(3),h%s%d(4))
    call create(omega,omega_s)
    select case (ireconst)
    case(wenojs5_id)
      call cmpt_commweight_wenojs5(omega,measure,omega_lo,omega_hi,omega_padding,bias,axis)
    case(wenoz5_id)
      call cmpt_commweight_wenoz5(omega,measure,omega_lo,omega_hi,omega_padding,bias,axis)
    case(wenoe5_id)
      call cmpt_commweight_wenoe5(omega,measure,omega_lo,omega_hi,omega_padding,bias,axis)
    case default
      print*,"Error: reconst_commweight()"
      print*,"No such ireconst: ",ireconst
      stop
    end select
    call reconst_commweight_matrix_weno5(h,f,omega,lo,hi,padding,bias,axis)
    call destroy(omega)
  end subroutine reconst_commweight_matrix
  
  subroutine reconst_commweight_tensor(h,f,measure,lo,hi,padding,bias,axis,ireconst)
    implicit none
    type(rtensorfield_t),intent(inout) :: h
    type(rtensorfield_t),intent(inout) :: f
    type(rscalarfield_t),intent(in)    :: measure
    type(ftensorshape_t),intent(in)    :: lo,hi,padding
    integer,intent(in)                 :: bias,axis
    integer,intent(in)                 :: ireconst
    
    type(rvectorfield_t) :: omega
    type(fvectorshape_t) :: omega_s
    type(fscalarshape_t) :: omega_lo,omega_hi,omega_padding
    
    omega_lo      = set_shape(lo%d(4),lo%d(5))
    omega_hi      = set_shape(hi%d(4),hi%d(5))
    omega_padding = set_shape(padding%d(4),padding%d(5))
    
    omega_s = set_shape(3,h%s%d(4),h%s%d(5))
    call create(omega,omega_s)
    select case (ireconst)
    case(wenojs5_id)
      call cmpt_commweight_wenojs5(omega,measure,omega_lo,omega_hi,omega_padding,bias,axis)
    case(wenoz5_id)
      call cmpt_commweight_wenoz5(omega,measure,omega_lo,omega_hi,omega_padding,bias,axis)
    case(wenoe5_id)
      call cmpt_commweight_wenoe5(omega,measure,omega_lo,omega_hi,omega_padding,bias,axis)
    case default
      print*,"Error: reconst_commweight()"
      print*,"No such ireconst: ",ireconst
      stop
    end select
    call reconst_commweight_tensor_weno5(h,f,omega,lo,hi,padding,bias,axis)
    call destroy(omega)
  end subroutine reconst_commweight_tensor
end module reconstruction