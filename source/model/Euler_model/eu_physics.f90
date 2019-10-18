module eu_physics
  use data_type
  use eu_data_type
  use mesh
  implicit none
contains
subroutine cmpt_prim_eu(flwd)
  implicit none
  type(eu_data_t),target,intent(inout) :: flwd

  type(eu_info_data_t),pointer :: info
  type(rvectorfield_t),pointer :: prim
  type(rvectorfield_t),pointer :: cons 

  integer :: isize, jsize, ibf, jbf
  type(fscalarshape_t) :: lo_int,hi_int
  integer :: pre_stage

  info => flwd%info

  if(info%stage_cur.eq.1) then
    pre_stage = info%stage_num
  else
    pre_stage = info%stage_cur - 1
  endif

  prim => flwd%prim
  cons => flwd%cons(pre_stage)

  isize = info%isize
  jsize = info%jsize
  ibf   = info%ibf
  jbf   = info%jbf

  ! interior
  lo_int = set_shape(1+ibf,1+jbf)
  hi_int = set_shape(isize+ibf,jsize+jbf)
  call cmpt_cons2prim_range(prim, cons, lo_int, hi_int, info) 

end subroutine cmpt_prim_eu

! convert conservative variables to primative variables
subroutine cmpt_cons2prim_range(prim, cons, lo, hi, info)
  implicit none
  type(fscalarshape_t),intent(in) :: lo,hi
  type(rvectorfield_t),target,intent(in) :: cons
  type(eu_info_data_t),intent(inout) :: info
  type(rvectorfield_t),target,intent(inout) :: prim

  double precision,dimension(:),pointer :: cons_elem 
  double precision,dimension(:),pointer :: prim_elem 

  double precision :: ru,rv,ret
  double precision :: r,u,v,et
  double precision :: t,p
  double precision :: c2
  double precision :: ei
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

  do j0 = 0, jrng - 1
    do i0 = 0, irng - 1
      i = lo%d(1) + incr1*i0
      j = lo%d(2) + incr2*j0
      cons_elem => cons%p(:,i,j)
      prim_elem => prim%p(:,i,j)
      
      r  = cons_elem(1)
      u  = cons_elem(2)/r
      v  = cons_elem(3)/r
      et = cons_elem(4)/r
      ei = et - 0.5d0*(u**2 + v**2)
      p  = r*ei*(info%gamma-1.d0)
      t  = ei
      c2 = info%gamma*p/r
      
      if(c2<0.d0) then
        print*,"negative sound speed",r,u,v,c2,i,j,info%idx
        p   = 1.0e-6
        ei  = p/(r*(info%gamma-1.d0))
        et  = r*(ei+0.5d0*(u**2+v**2))
        cons_elem(4) = r*et
        c2 = info%gamma*p/r
      endif
    
      prim_elem(1) = r
      prim_elem(2) = u
      prim_elem(3) = v
      prim_elem(4) = p
      prim_elem(5) = sqrt(c2)
      prim_elem(6) = info%idx
           
      info%eig_xi  = max(info%eig_xi,abs(u)+sqrt(c2))
      info%eig_eta = max(info%eig_eta,abs(v)+sqrt(c2))
    enddo
  enddo
end subroutine cmpt_cons2prim_range

! convert primative variables to conservative variables
subroutine cmpt_prim2cons_range(cons, prim, lo, hi, info)
  implicit none
  type(fscalarshape_t),intent(in) :: lo,hi
  type(rvectorfield_t),target,intent(in) :: prim
  type(eu_info_data_t),intent(inout) :: info
  type(rvectorfield_t),target,intent(inout) :: cons

  double precision,dimension(:),pointer :: cons_elem 
  double precision,dimension(:),pointer :: prim_elem 

  double precision :: ru,rv,ret
  double precision :: r,u,v,et
  double precision :: t,p
  double precision :: c2
  double precision :: ei
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

  do j0 = 0, jrng - 1
    do i0 = 0, irng - 1
      i = lo%d(1) + incr1*i0
      j = lo%d(2) + incr2*j0
      cons_elem => cons%p(:,i,j)
      prim_elem => prim%p(:,i,j)
      
      r  = prim_elem(1)
      u  = prim_elem(2)
      v  = prim_elem(3)
      p  = prim_elem(4)

      ei = p/(r*(info%gamma-1.d0))
      et = ei + 0.5d0*(u**2 + v**2)
          
      cons_elem(1) = r
      cons_elem(2) = r*u
      cons_elem(3) = r*v
      cons_elem(4) = r*et      
    enddo
  enddo
end subroutine cmpt_prim2cons_range
end module eu_physics