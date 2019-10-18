submodule (reconstruction) uw1
implicit none

contains
    module subroutine reconst_component_scalar_uw1(h,f,lo,hi,padding,bias,axis)
        implicit none
        type(rscalarfield_t),intent(inout) :: h
        type(rscalarfield_t),target,intent(in) :: f
        type(fscalarshape_t),intent(in) :: lo,hi,padding
        integer,intent(in) :: bias,axis
        integer :: i,j
        integer :: i0,j0
        integer :: ishift,jshift
        select case(axis)
        case(axis_xi)
            ishift = padding%d(1) + bias
            jshift = padding%d(2)
        case(axis_eta)
            ishift = padding%d(1)
            jshift = padding%d(2) + bias
        case default
            print*,"Error: reconst_component_scalar_uw1()"
            print*,"No such ireconst: ",axis
            stop
        end select 
        do j=lo%d(2), hi%d(2)
            do i=lo%d(1), hi%d(1)
            i0 = i + ishift
            j0 = j + jshift
            h%p(i,j) = f%p(i0,j0)
            enddo
        enddo
    end subroutine reconst_component_scalar_uw1

    module subroutine reconst_component_vector_uw1(h,f,lo,hi,padding,bias,axis)
        implicit none
        type(rvectorfield_t),intent(inout) :: h
        type(rvectorfield_t),target,intent(in) :: f
        type(fvectorshape_t),intent(in) :: lo,hi,padding
        integer,intent(in) :: bias,axis
        integer :: n,i,j
        integer :: n0,i0,j0
        integer :: ishift,jshift,nshift
        select case(axis)
        case(axis_xi)
            nshift = padding%d(1)
            ishift = padding%d(2) + bias
            jshift = padding%d(3)
        case(axis_eta)
            nshift = padding%d(1)
            ishift = padding%d(2)
            jshift = padding%d(3) + bias
        case default
            print*,"Error: reconst_component_vector_uw1()"
            print*,"No such ireconst: ",axis
            stop
        end select 
        do j=lo%d(3), hi%d(3)
            do i=lo%d(2), hi%d(2)
            do n=lo%d(1), hi%d(1)
                i0 = i + ishift
                j0 = j + jshift
                n0 = n + nshift
                h%p(n,i,j) = f%p(n0,i0,j0)
            enddo
            enddo
        enddo
    end subroutine reconst_component_vector_uw1

    module subroutine reconst_component_matrix_uw1(h,f,lo,hi,padding,bias,axis)
        implicit none
        type(rmatrixfield_t),intent(inout) :: h
        type(rmatrixfield_t),target,intent(in) :: f
        type(fmatrixshape_t),intent(in) :: lo,hi,padding
        integer,intent(in) :: bias,axis
        integer :: n,m,i,j
        integer :: n0,m0,i0,j0
        integer :: ishift,jshift,nshift,mshift
        select case(axis)
        case(axis_xi)
            nshift = padding%d(1)
            mshift = padding%d(2)
            ishift = padding%d(3) + bias
            jshift = padding%d(4)
        case(axis_eta)
            nshift = padding%d(1)
            mshift = padding%d(2)
            ishift = padding%d(3)
            jshift = padding%d(4) + bias
        case default
            print*,"Error: reconst_component_matrix_uw1()"
            print*,"No such ireconst: ",axis
            stop
        end select 
        do j=lo%d(4), hi%d(4)
            do i=lo%d(3), hi%d(3)
            do m=lo%d(2), hi%d(2)
                do n=lo%d(1), hi%d(1)
                i0 = i + ishift
                j0 = j + jshift
                m0 = m + mshift
                n0 = n + nshift

                h%p(n,m,i,j) = f%p(n0,m0,i0,j0)
                enddo
            enddo
            enddo
        enddo
    end subroutine reconst_component_matrix_uw1

    module subroutine reconst_component_tensor_uw1(h,f,lo,hi,padding,bias,axis)
        implicit none
        type(rtensorfield_t),intent(inout) :: h
        type(rtensorfield_t),target,intent(in) :: f
        type(ftensorshape_t),intent(in) :: lo,hi,padding
        integer,intent(in) :: bias,axis
        integer :: l,m,n,i,j
        integer :: l0,m0,n0,i0,j0
        integer :: ishift,jshift,nshift,mshift,lshift
        select case(axis)
        case(axis_xi)
            nshift = padding%d(1)
            mshift = padding%d(2)
            lshift = padding%d(3)
            ishift = padding%d(4) + bias
            jshift = padding%d(5)
        case(axis_eta)
            nshift = padding%d(1)
            mshift = padding%d(2)
            lshift = padding%d(3)
            ishift = padding%d(4)
            jshift = padding%d(5) + bias
        case default
            print*,"Error: reconst_component_tensor_uw1()"
            print*,"No such ireconst: ",axis
            stop
        end select 
        do j=lo%d(5), hi%d(5)
            do i=lo%d(4), hi%d(4)
            do l=lo%d(3), hi%d(3)
                do m=lo%d(2), hi%d(2)
                do n=lo%d(1), hi%d(1)
                    i0 = i + ishift
                    j0 = j + jshift
                    l0 = l + lshift
                    m0 = m + mshift
                    n0 = n + nshift
                    h%p(n,m,l,i,j) = f%p(n0,m0,l0,i0,j0)
                enddo
                enddo
            enddo
            enddo
        enddo
    end subroutine reconst_component_tensor_uw1
end submodule