module time_advance
	use data_type
implicit none

integer,parameter :: euler1st_idx = 1
integer,parameter :: rk2nd_idx    = 2
integer,parameter :: rk3rdtvd_idx = 3
integer,parameter :: irk3rd_idx   = 4

contains
integer function set_maxstage(time_advance_idx)
	implicit none
	integer,intent(in) :: time_advance_idx
	select case(time_advance_idx)
	case(euler1st_idx)
		set_maxstage = 1
	case(rk2nd_idx)
		set_maxstage = 2
	case(rk3rdtvd_idx)
		set_maxstage = 3
	case default
		print*,"Error: compute_time_advance_stage()"
		print*,"No such time method:", time_advance_idx
		set_maxstage = 1
		stop
	end select 
end function set_maxstage

subroutine cmpt_time_advn_stg(u,                  &
	                            lhs,rhs,            &
															dt,                 &
															lo,hi,              &
															istage,             &
															max_stage,          &
															time_integrator_idx)
implicit none
integer,intent(in) :: time_integrator_idx
integer,intent(in) :: max_stage
type(fvectorshape_t),intent(in) :: lo,hi
double precision,intent(in) :: dt
integer,intent(inout) :: istage
type(rvectorfield_t),dimension(max_stage),intent(inout) :: u,lhs,rhs
select case(time_integrator_idx)
case(euler1st_idx)
	call compute_euler1st_stage(u,rhs,   &
															dt,      &
															lo,hi,   &
															istage,  &
															max_stage)
case(rk2nd_idx)
	call compute_rk2ndssp_stage(u,rhs,   &
															dt,      &
															lo,hi,   &
															istage,  &
															max_stage)
case(rk3rdtvd_idx)
	call compute_rk3rdssp_stage(u,rhs,   &
															dt,      &
															lo,hi,   &
															istage,  &
															max_stage)
case default
	print*,"Error: compute_time_advance_stage()"
	print*,"No such time method:", time_integrator_idx
end select 
end subroutine cmpt_time_advn_stg

subroutine compute_euler1st_stage(u,rhs,   &
																	dt,      &
																	lo,hi,   &
																	istage,  &
																	max_stage)
implicit none
double precision,intent(in) :: dt
type(fvectorshape_t),intent(in) :: lo,hi
integer,intent(in) :: max_stage
integer,intent(inout) :: istage
type(rvectorfield_t),dimension(max_stage),intent(inout) :: u,rhs
integer :: i,j,n
select case(istage)
case(1)
  do j = lo%d(3), hi%d(3)
    do i = lo%d(2), hi%d(2)
			do n = lo%d(1), hi%d(1)
				u(istage)%p(n,i,j) = u(max_stage)%p(n,i,j) + dt*rhs(istage)%p(n,i,j)
			enddo
		enddo
	enddo
case default
	print*,"Error: compute_euler1st_stage()"
	print*,"Out of steges of rk3tve:",istage," of 1"
end select
end subroutine compute_euler1st_stage

subroutine compute_rk2ndssp_stage(u,rhs,   &
																	dt,      &
																	lo,hi,   &
																	istage,  &
																	max_stage)
implicit none
double precision,intent(in) :: dt
type(fvectorshape_t),intent(in) :: lo,hi
integer,intent(in) :: max_stage
integer,intent(inout) :: istage
type(rvectorfield_t),dimension(max_stage),intent(inout) :: u,rhs
integer :: i,j,n
select case(istage)
case(1)
  do j = lo%d(3), hi%d(3)
    do i = lo%d(2), hi%d(2)
			do n = lo%d(1), hi%d(1)
				u(istage)%p(n,i,j) = u(max_stage)%p(n,i,j) + dt*rhs(istage)%p(n,i,j)
			enddo
		enddo
	enddo
case(2)
  do j = lo%d(3), hi%d(3)
    do i = lo%d(2), hi%d(2)
			do n = lo%d(1), hi%d(1)
				u(istage)%p(n,i,j) = 0.5d0*u(max_stage)%p(n,i,j) + 0.5d0*(u(istage-1)%p(n,i,j) + dt*rhs(istage)%p(n,i,j))
			enddo
		enddo
	enddo
case default
	print*,"Error: compute_rk2ndssp_stage()"
	print*,"Out of steges of rk3tve:",istage," of 2"
end select
end subroutine compute_rk2ndssp_stage

subroutine compute_rk3rdssp_stage(u,rhs,   &
																	dt,      &
																	lo,hi,   &
																	istage,  &
																	max_stage)
implicit none
double precision,intent(in) :: dt
type(fvectorshape_t),intent(in) :: lo,hi
integer,intent(in) :: max_stage
integer,intent(inout) :: istage
type(rvectorfield_t),dimension(max_stage),intent(inout) :: u,rhs
integer :: i,j,n
select case(istage)
case(1)
  do j = lo%d(3), hi%d(3)
		do i = lo%d(2), hi%d(2)
			do n = lo%d(1), hi%d(1)
				u(istage)%p(n,i,j) = u(max_stage)%p(n,i,j) + dt*rhs(istage)%p(n,i,j)
			enddo
		enddo
	enddo
case(2)
  do j = lo%d(3), hi%d(3)
    do i = lo%d(2), hi%d(2)
			do n = lo%d(1), hi%d(1)
				u(istage)%p(n,i,j) = 3.d0/4.d0*u(max_stage)%p(n,i,j) + 1.d0/4.d0*(u(istage-1)%p(n,i,j) + dt*rhs(istage)%p(n,i,j))
			enddo
		enddo
	enddo
case(3)
  do j = lo%d(3), hi%d(3)
    do i = lo%d(2), hi%d(2)
			do n = lo%d(1), hi%d(1)
				u(istage)%p(n,i,j) = 1.d0/3.d0*u(max_stage)%p(n,i,j) + 2.d0/3.d0*(u(istage-1)%p(n,i,j) + dt*rhs(istage)%p(n,i,j))
			enddo
		enddo
	enddo
case default
	print*,"Out of steges of rk3rdssp:",istage," of 3"
	stop
end select
end subroutine compute_rk3rdssp_stage
end module
