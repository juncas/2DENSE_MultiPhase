module IO
  use data_type
  implicit none  
  interface write_field_fmt
    module procedure write_scalarfield_fmt
    module procedure write_vectorfield_fmt
    module procedure write_matrixfield_fmt
    module procedure write_tensorfield_fmt
  end interface write_field_fmt  
  
  interface write_field_unf
    module procedure write_scalarfield_unf
    module procedure write_vectorfield_unf
    module procedure write_matrixfield_unf
    module procedure write_tensorfield_unf
  end interface write_field_unf  
  
  interface write_field_hdf5
    module procedure write_scalarfield_hdf5
    module procedure write_vectorfield_hdf5
    module procedure write_matrixfield_hdf5
    module procedure write_tensorfield_hdf5
  end interface write_field_hdf5  
  
  interface read_field_unf
    module procedure read_scalarfield_unf
    module procedure read_vectorfield_unf
    module procedure read_matrixfield_unf
    module procedure read_tensorfield_unf
  end interface read_field_unf  

  ! interface read_field_hdf5
  !   module procedure read_scalarfield_hdf5
  !   module procedure read_vectorfield_hdf5
  !   module procedure read_matrixfield_hdf5
  !   module procedure read_tensorfield_hdf5
  ! end interface read_field_hdf5  
contains
subroutine write_scalarfield_fmt(funit,f,lo,hi)
  implicit none
  integer,intent(in) :: funit
  type(rscalarfield_t),intent(in) :: f
  type(fscalarshape_t),intent(in) :: lo,hi
  integer :: i,j
  do j = lo%d(2), hi%d(2)
    do i = lo%d(1), hi%d(1)
        write(funit,*) f%p(i,j)
    enddo
  enddo
end subroutine write_scalarfield_fmt

subroutine write_vectorfield_fmt(funit,f,lo,hi)
  implicit none
  integer,intent(in) :: funit
  type(rvectorfield_t),intent(in) :: f
  type(fvectorshape_t),intent(in) :: lo,hi
  integer :: n,i,j
  do n = lo%d(1), hi%d(1)
    do j = lo%d(3), hi%d(3)
      do i = lo%d(2), hi%d(2)
        write(funit,*) f%p(n,i,j)
      enddo
    enddo
  enddo
end subroutine write_vectorfield_fmt

subroutine write_matrixfield_fmt(funit,f,lo,hi)
  implicit none
  integer,intent(in) :: funit
  type(rmatrixfield_t),intent(in) :: f
  type(fmatrixshape_t),intent(in) :: lo,hi
  integer :: m,n,i,j
  do m = lo%d(1), hi%d(1)
    do n = lo%d(2), hi%d(2)
      do j = lo%d(4), hi%d(4)
        do i = lo%d(3), hi%d(3)
          write(funit,*) f%p(m,n,i,j)
        enddo
      enddo
    enddo
  enddo
end subroutine write_matrixfield_fmt

subroutine write_tensorfield_fmt(funit,f,lo,hi)
implicit none
integer,intent(in) :: funit
type(rtensorfield_t),intent(in) :: f
type(ftensorshape_t),intent(in) :: lo,hi
integer :: l,m,n,i,j
do l = lo%d(1), hi%d(1)
  do m = lo%d(2), hi%d(2)
    do n = lo%d(3), hi%d(3)
      do j = lo%d(5), hi%d(5)
        do i = lo%d(4), hi%d(4)
          write(funit,*) f%p(l,m,n,i,j)
        enddo
      enddo
    enddo
  enddo
enddo
end subroutine write_tensorfield_fmt

subroutine write_scalarfield_unf(funit,f,lo,hi)
  implicit none
  integer,intent(in) :: funit
  type(rscalarfield_t),intent(in) :: f
  type(fscalarshape_t),intent(in) :: lo,hi
  integer :: i,j
  do j = lo%d(2), hi%d(2)
    do i = lo%d(1), hi%d(1)
        write(funit) f%p(i,j)
    enddo
  enddo
end subroutine write_scalarfield_unf

subroutine write_vectorfield_unf(funit,f,lo,hi)
  implicit none
  integer,intent(in) :: funit
  type(rvectorfield_t),intent(in) :: f
  type(fvectorshape_t),intent(in) :: lo,hi
  integer :: i,j
  do j = lo%d(3), hi%d(3)
    do i = lo%d(2), hi%d(2)    
      write(funit) f%p(lo%d(1):hi%d(1),i,j)
    enddo
  enddo
end subroutine write_vectorfield_unf

subroutine write_matrixfield_unf(funit,f,lo,hi)
  implicit none
  integer,intent(in) :: funit
  type(rmatrixfield_t),intent(in) :: f
  type(fmatrixshape_t),intent(in) :: lo,hi
  integer :: i,j
  do j = lo%d(4), hi%d(4)
    do i = lo%d(3), hi%d(3)    
      write(funit) f%p(lo%d(1):hi%d(1),lo%d(2):hi%d(2),i,j)
    enddo
  enddo
end subroutine write_matrixfield_unf

subroutine write_tensorfield_unf(funit,f,lo,hi)
implicit none
integer,intent(in) :: funit
type(rtensorfield_t),intent(in) :: f
type(ftensorshape_t),intent(in) :: lo,hi
integer :: i,j
do j = lo%d(5), hi%d(5)
  do i = lo%d(4), hi%d(4)
    write(funit) f%p(lo%d(1):hi%d(1),lo%d(2):hi%d(2),lo%d(3):hi%d(3),i,j)
  enddo
enddo
end subroutine write_tensorfield_unf

subroutine write_scalarfield_hdf5(f,lo,hi,name,file_id)
  use HDF5
  implicit none
  type(rscalarfield_t),target,intent(in) :: f
  type(fscalarshape_t),intent(in) :: lo,hi
  CHARACTER(LEN=*),intent(in)     :: name
  integer(HID_T),intent(in)       :: file_id
 
  integer(HID_T)                  :: dspace_id
  integer(HID_T)                  :: dset_id
  integer(HSIZE_T), DIMENSION(2)  :: dim 
  integer                         :: error
  double precision,pointer,dimension(:,:) :: p

  dim(1) = abs(hi%d(1) - lo%d(1)) + 1
  dim(2) = abs(hi%d(2) - lo%d(2)) + 1  
  
  p => f%p(lo%d(1):hi%d(1),lo%d(2):hi%d(2))
  CALL h5screate_simple_f(2, dim, dspace_id, error)
  CALL h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, dspace_id, &
                   dset_id, error)
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, p, dim, error)
  CALL h5dclose_f(dset_id, error)
  CALL h5sclose_f(dspace_id, error)

end subroutine write_scalarfield_hdf5

subroutine write_vectorfield_hdf5(f,lo,hi,name,file_id)
  use HDF5
  implicit none
  type(rvectorfield_t),target,intent(in)   :: f
  type(fvectorshape_t),intent(in)          :: lo,hi
  CHARACTER(LEN=*),dimension(*),intent(in) :: name
  integer(HID_T),intent(in)                :: file_id
 
  integer(HID_T)                          :: dspace_id
  integer(HID_T)                          :: dset_id
  integer(HSIZE_T), DIMENSION(3)          :: dim 
  integer                                 :: error
  integer                                 :: i
  double precision,pointer,dimension(:,:) :: p

  dim(1) = abs(hi%d(2) - lo%d(2)) + 1 
  dim(2) = abs(hi%d(3) - lo%d(3)) + 1  
  do i = lo%d(1), hi%d(1)  
    p => f%p(i,lo%d(2):hi%d(2),lo%d(3):hi%d(3))
    CALL h5screate_simple_f(2, dim, dspace_id, error)
    CALL h5dcreate_f(file_id, name(i), H5T_NATIVE_DOUBLE, dspace_id, &
                    dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, p, dim, error)
    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
  enddo

end subroutine write_vectorfield_hdf5

subroutine write_matrixfield_hdf5(f,lo,hi,name,file_id)
  use HDF5
  implicit none
  type(rmatrixfield_t),target,intent(in) :: f
  type(fmatrixshape_t),intent(in) :: lo,hi
  CHARACTER(LEN=*),intent(in)     :: name
  integer(HID_T),intent(in)       :: file_id
 
  integer(HID_T)                  :: dspace_id
  integer(HID_T)                  :: dset_id
  integer(HSIZE_T), DIMENSION(2)  :: dim 
  integer                         :: error
  integer                         :: i, j
  double precision,pointer,dimension(:,:) :: p

  dim(1) = abs(hi%d(3) - lo%d(3)) + 1
  dim(2) = abs(hi%d(4) - lo%d(4)) + 1 
  
  do j = lo%d(2), hi%d(2)  
    do i = lo%d(1), hi%d(1)  
      p => f%p(i,j,lo%d(3):hi%d(3),lo%d(4):hi%d(4))
      CALL h5screate_simple_f(2, dim, dspace_id, error)
      CALL h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, p, dim, error)
      CALL h5dclose_f(dset_id, error)
      CALL h5sclose_f(dspace_id, error)
    enddo
  enddo

end subroutine write_matrixfield_hdf5

subroutine write_tensorfield_hdf5(f,lo,hi,name,file_id)
  use HDF5
  implicit none
  type(rtensorfield_t),target,intent(in) :: f
  type(ftensorshape_t),intent(in) :: lo,hi
  CHARACTER(LEN=*),intent(in)     :: name
  integer(HID_T),intent(in)       :: file_id

  integer(HID_T)                  :: dspace_id
  integer(HID_T)                  :: dset_id
  integer(HSIZE_T), DIMENSION(2)  :: dim 
  integer                         :: error
  integer                         :: i, j, k
  double precision,pointer,dimension(:,:) :: p

  dim(1) = abs(hi%d(4) - lo%d(4)) + 1
  dim(2) = abs(hi%d(5) - lo%d(5)) + 1
  do k = lo%d(3), hi%d(3)  
    do j = lo%d(2), hi%d(2)  
      do i = lo%d(1), hi%d(1)  
        p => f%p(i,j,k,lo%d(4):hi%d(4),lo%d(5):hi%d(5))

        CALL h5screate_simple_f(2, dim, dspace_id, error)
        CALL h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, dspace_id, &
                        dset_id, error)
        CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, &
                        p, &
                        dim, error)
        CALL h5dclose_f(dset_id, error)
        CALL h5sclose_f(dspace_id, error)
      enddo
    enddo
  enddo
end subroutine write_tensorfield_hdf5

subroutine read_scalarfield_unf(funit,f,lo,hi)
  implicit none
  integer,intent(in) :: funit
  type(fscalarshape_t),intent(in) :: lo,hi
  type(rscalarfield_t),intent(inout) :: f
  integer :: i,j
  do j = lo%d(2), hi%d(2)
    do i = lo%d(1), hi%d(1)
        read(funit) f%p(i,j)
    enddo
  enddo
end subroutine read_scalarfield_unf

subroutine read_vectorfield_unf(funit,f,lo,hi)
  implicit none
  integer,intent(in) :: funit
  type(fvectorshape_t),intent(in) :: lo,hi
  type(rvectorfield_t),intent(inout) :: f
  integer :: i,j
  do j = lo%d(3), hi%d(3)
    do i = lo%d(2), hi%d(2)
      read(funit) f%p(lo%d(1):hi%d(1),i,j)
    enddo
  enddo
end subroutine read_vectorfield_unf

subroutine read_matrixfield_unf(funit,f,lo,hi)
  implicit none
  integer,intent(in) :: funit
  type(fmatrixshape_t),intent(in) :: lo,hi
  type(rmatrixfield_t),intent(inout) :: f
  integer :: i,j
  do j = lo%d(4), hi%d(4)
    do i = lo%d(3), hi%d(3)
      read(funit) f%p(lo%d(1):hi%d(1),lo%d(2):hi%d(2),i,j)
    enddo
  enddo
end subroutine read_matrixfield_unf

subroutine read_tensorfield_unf(funit,f,lo,hi)
implicit none
integer,intent(in) :: funit
type(ftensorshape_t),intent(in) :: lo,hi
type(rtensorfield_t),intent(inout) :: f
integer :: i,j
do j = lo%d(5), hi%d(5)
  do i = lo%d(4), hi%d(4)
    read(funit) f%p(lo%d(1):hi%d(1),lo%d(2):hi%d(2),lo%d(3):hi%d(3),i,j)
  enddo
enddo
end subroutine read_tensorfield_unf
end module IO
