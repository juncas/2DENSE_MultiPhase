module eu_boundary
  use data_type
  use boundary
  use mesh
  use eu_data_type
  use eu_physics
  implicit none
  integer,parameter :: sv_bc_idx         = 2
  integer,parameter :: dm_bc_idx         = 3
  integer,parameter :: ssl_bc_idx        = 4 
  integer,parameter :: sb_bc_idx         = 5 
  integer,parameter :: sdv_bc_idx        = 6 
  integer,parameter :: rt_bc_idx         = 7 
  integer,parameter :: rm_bc_idx         = 8 
  integer,parameter :: fs_bc_idx         = 9 
  integer,parameter :: bs_bc_idx         = 10 
  integer,parameter :: shocktube_bc_idx = 11 
  integer,parameter :: blast_bc_idx      = 14

  integer,parameter :: otlt1_sv_idx = 201
  integer,parameter :: rflt2_sv_idx = 202
  integer,parameter :: otlt3_sv_idx = 203
  integer,parameter :: rflt4_sv_idx = 204

  integer,parameter :: inlt1_dm_idx = 301
  integer,parameter :: rflt2_dm_idx = 302
  integer,parameter :: otlt3_dm_idx = 303
  integer,parameter :: inlt4_dm_idx = 304

  integer,parameter :: inlt1_ssl_idx = 401
  integer,parameter :: rflt2_ssl_idx = 402
  integer,parameter :: otlt3_ssl_idx = 403
  integer,parameter :: inlt4_ssl_idx = 404

  integer,parameter :: inlt1_sb_idx = 501
  integer,parameter :: rflt2_sb_idx = 502
  integer,parameter :: otlt3_sb_idx = 503
  integer,parameter :: inlt4_sb_idx = 504

  integer,parameter :: inlt1_sedov_idx = 601
  integer,parameter :: rflt2_sedov_idx = 602
  integer,parameter :: otlt3_sedov_idx = 603
  integer,parameter :: inlt4_sedov_idx = 604

  integer,parameter :: inlt1_rt_idx = 701
  integer,parameter :: rflt2_rt_idx = 702
  integer,parameter :: otlt3_rt_idx = 703
  integer,parameter :: inlt4_rt_idx = 704

  integer,parameter :: inlt1_rm_idx = 801
  integer,parameter :: rflt2_rm_idx = 802
  integer,parameter :: otlt3_rm_idx = 803
  integer,parameter :: inlt4_rm_idx = 804

  integer,parameter :: inlt11_fs_idx = 901
  integer,parameter :: rflt12_fs_idx = 902
  integer,parameter :: rflt13_fs_idx = 903
  integer,parameter :: intf14_fs_idx = 904
  integer,parameter :: rflt15_fs_idx = 905
  integer,parameter :: intf21_fs_idx = 906
  integer,parameter :: rflt22_fs_idx = 907
  integer,parameter :: otlt23_fs_idx = 908
  integer,parameter :: rflt24_fs_idx = 909

  integer,parameter :: inlt1_bs_idx = 1001
  integer,parameter :: rflt2_bs_idx = 1002
  integer,parameter :: otlt3_bs_idx = 1003
  integer,parameter :: inlt4_bs_idx = 1004

  integer,parameter :: otlt1_st_idx = 1101
  integer,parameter :: rflt2_st_idx = 1102
  integer,parameter :: otlt3_st_idx = 1103
  integer,parameter :: rflt4_st_idx = 1104

  integer,parameter :: rflt1_blst_idx = 1401
  integer,parameter :: rflt2_blst_idx = 1402
  integer,parameter :: rflt3_blst_idx = 1403
  integer,parameter :: rflt4_blst_idx = 1404

contains
subroutine set_bc_eu(flwds, nflwd, iflwd)
  implicit none
  integer,intent(in) :: nflwd, iflwd
  type(eu_data_iterator_t),target,dimension(nflwd),intent(inout) :: flwds

  type(eu_data_t),pointer     :: flwd
  type(block_t),pointer       :: blk
  type(block_data_t),pointer  :: d_blk
  type(bndy_list_t),pointer   :: bndy_node
  type(bndy_t),pointer        :: bndy
  type(bndy_data_t),pointer   :: d_bndy
  type(ivector_t),pointer     :: bndy_bc_info
  integer :: bc_group

  type(eu_info_data_t),pointer :: info

  integer :: pre_stage

  integer :: ibndy
  integer :: bndy_bc_idx

  flwd => flwds(iflwd)%flwd
  blk  => flwd%blk
  info => flwd%info

  bndy_node => blk%bndy_list
  do while(associated(bndy_node))
    bndy         => bndy_node%bndy
    bndy_bc_info => bndy%data%bc_infos
    bc_group     = bndy_bc_info%p(1)

    ! General boundary conditions
    select case(bc_group)
      case(bc_idx_intf)
        call set_intf_bc(flwds, nflwd, iflwd, bndy)
      case(bc_idx_inlt)
        ! call set_inlt_bc(flwds, nflwd, iflwd, bndy)
      case(bc_idx_otlt)
        ! call set_otlt_bc(flwds, nflwd, iflwd, bndy)
      case(bc_idx_rflt)
        ! call set_rflt_bc(flwds, nflwd, iflwd, bndy)
    end select

    ! Specific boundary conditions
    select case(bc_group/100)
      case(sv_bc_idx)
        ! call set_sv_bc(flwds, nflwd, iflwd, bndy)
      case(dm_bc_idx)
        call set_dm_bc(flwds, nflwd, iflwd, bndy)
      case(ssl_bc_idx)
        ! call set_ssl_bc(flwds, nflwd, iflwd, bndy)
      case(sb_bc_idx)
        ! call set_sb_bc(flwds, nflwd, iflwd, bndy)
      case(sdv_bc_idx)
        ! call set_sdv_bc(flwds, nflwd, iflwd, bndy)
      case(rt_bc_idx)
        ! call set_rt_bc(flwds, nflwd, iflwd, bndy)
      case(rm_bc_idx)
        ! call set_rm_bc(flwds, nflwd, iflwd, bndy)
      case(fs_bc_idx)
        call set_fs_bc(flwds, nflwd, iflwd, bndy)
      case(bs_bc_idx)
        ! call set_bs_bc(flwds, nflwd, iflwd, bndy)
      case(shocktube_bc_idx)
        ! call set_st_bc(flwds, nflwd, iflwd, bndy) 
      case(blast_bc_idx)
        ! call set_blst_bc(flwds, nflwd, iflwd, bndy)
      ! case default
      !   print*,"Error set_bc_eu():"
      !   print*,"No such BC group index:",bc_group, &
      !          "boundary index:",bndy%idx,         &
      !          "flwd index:", iflwd
      !   stop
    end select
    bndy_node => bndy_node%next
  enddo
end subroutine set_bc_eu

subroutine set_dm_bc(flwds, nflwd,      &
                     iflwd,             & 
                     bndy)
  implicit none
  integer,intent(in) :: nflwd
  type(eu_data_iterator_t),target,dimension(nflwd),intent(inout) :: flwds
  integer,intent(in) :: iflwd
  type(bndy_t),target,intent(in) :: bndy
  type(bndy_data_t),pointer      :: d_bndy
  type(ivector_t),pointer        :: bndy_bc_info

  type(eu_data_t),pointer      :: flwd
  type(eu_info_data_t),pointer :: info
  type(rvectorfield_t),pointer :: cons
  type(rvectorfield_t),pointer :: prim
  type(block_t),pointer        :: blk
  type(block_data_t),pointer   :: d_blk

  integer :: iflwd_d
  type(eu_data_t),pointer      :: flwd_d
  type(eu_info_data_t),pointer :: info_d
  type(rvectorfield_t),pointer :: cons_d

  type(fscalarshape_t) :: lo,hi
  type(fscalarshape_t) :: lo_d,hi_d

  integer :: istage

  integer :: bc_type
  type(rvector_t),pointer :: val
  integer :: axis
  integer :: cnnct
  integer :: irng,jrng,i0,j0,i,j
  integer :: incr1,incr2

  double precision :: xc0,yc0
  double precision :: rho0,u0,v0,p0,ei0
  
  flwd  => flwds(iflwd)%flwd
  blk   => flwd%blk
  d_blk => blk%data
  info  => flwd%info

  if(info%stage_cur.eq.1) then
    istage = info%stage_num
  else
    istage = info%stage_cur - 1
  endif  

  cons => flwd%cons(istage)

  d_bndy       => bndy%data
  bndy_bc_info => d_bndy%bc_infos  
  bc_type      =  bndy_bc_info%p(1)
  
  lo = set_shape(d_bndy%ilo(1), d_bndy%jlo(1))
  hi = set_shape(d_bndy%ihi(1), d_bndy%jhi(1))
  
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
  
  select case(bc_type)
  case(inlt1_dm_idx)
    do j = 0, jrng - 1
      do i = 0, irng - 1
        i0 = lo%d(1) + incr1*i
        j0 = lo%d(2) + incr2*j

        rho0 = 8.d0
        u0   = 8.25*sqrt(3.d0)*0.5
        v0   =-8.25*0.5d0
        p0   = 116.5d0      
        ei0  = p0/(rho0*(info%gamma-1.d0))  

        cons%p(1,i0,j0) = rho0
        cons%p(2,i0,j0) = rho0*u0
        cons%p(3,i0,j0) = rho0*v0
        cons%p(4,i0,j0) = rho0*(ei0+0.5d0*(u0**2+v0**2))
      enddo
    enddo
  case(rflt2_dm_idx)         
    val => d_bndy%vector(1) 
    call boundary_reflect(cons, lo, hi, val, bndy%group_id)
    do j = 0, jrng - 1
      do i = 0, irng - 1
        i0 = lo%d(1) + incr1*i
        j0 = lo%d(2) + incr2*j
        xc0 = d_blk%xc%p(i0,j0)
        if(xc0.le.1.d0/6.d0) then
          rho0 = 8.d0
          u0   = 8.25*sqrt(3.d0)*0.5
          v0   =-8.25*0.5d0
          p0   = 116.5d0      
          ei0  = p0/(rho0*(info%gamma-1.d0))  
          cons%p(1,i0,j0) = rho0
          cons%p(2,i0,j0) = rho0*u0
          cons%p(3,i0,j0) = rho0*v0
          cons%p(4,i0,j0) = rho0*(ei0+0.5d0*(u0**2+v0**2))
        endif
      enddo
    enddo
  case(otlt3_dm_idx)
    axis = bndy_bc_info%p(2)
    val  => d_bndy%vector(1) 
    call boundary_neumann(cons, lo, hi, val, axis)
  case(inlt4_dm_idx)    
    do j = 0, jrng - 1
      do i = 0, irng - 1
        i0  = lo%d(1) + incr1*i
        j0  = lo%d(2) + incr2*j
        xc0 = d_blk%xc%p(i0,j0)
        yc0 = d_blk%yc%p(i0,j0)
        if(xc0.lt.(1.d0/6.d0+(yc0 + 20.d0*info%time)/sqrt(3.d0))) then
            rho0 = 8.d0
            u0   = 8.25*sqrt(3.d0)*0.5
            v0   =-8.25*0.5d0
            p0   = 116.5d0
        else
            rho0 = 1.4d0
            u0   = 0.d0
            v0   = 0.d0
            p0   = 1.d0
        endif       
        ei0  = p0/(rho0*(info%gamma-1.d0))  
        cons%p(1,i0,j0) = rho0
        cons%p(2,i0,j0) = rho0*u0
        cons%p(3,i0,j0) = rho0*v0
        cons%p(4,i0,j0) = rho0*(ei0+0.5d0*(u0**2+v0**2))
      enddo
    enddo
  case default
    print*,"Error set_dm_bc():"
    print*,"No such BC type index:",bc_type, &
           "boundary index:",bndy%idx,         &
           "flwd index:", iflwd
    stop
  end select

  prim => flwd%prim
  call cmpt_cons2prim_range(prim, cons, lo, hi, info) 
end subroutine set_dm_bc

subroutine set_fs_bc(flwds, nflwd,      &
                     iflwd,             & 
                     bndy)
  implicit none
  integer,intent(in) :: nflwd
  type(eu_data_iterator_t),target,dimension(nflwd),intent(inout) :: flwds
  integer,intent(in) :: iflwd
  type(bndy_t),target,intent(in) :: bndy
  type(bndy_data_t),pointer      :: d_bndy
  type(ivector_t),pointer        :: bndy_bc_info

  type(eu_data_t),pointer      :: flwd
  type(eu_info_data_t),pointer :: info
  type(rvectorfield_t),pointer :: cons
  type(rvectorfield_t),pointer :: prim
  type(block_t),pointer        :: blk
  type(block_data_t),pointer   :: d_blk

  type(fscalarshape_t) :: lo,hi

  integer :: istage

  integer :: bc_type
  type(rvector_t),pointer :: val
  integer :: axis
  integer :: cnnct
  integer :: irng,jrng,i0,j0,i,j
  integer :: incr1,incr2

  double precision :: xc0,yc0
  double precision :: rho0,u0,v0,p0,ei0
  
  flwd  => flwds(iflwd)%flwd
  blk   => flwd%blk
  d_blk => blk%data
  info  => flwd%info

  if(info%stage_cur.eq.1) then
    istage = info%stage_num
  else
    istage = info%stage_cur - 1
  endif  

  cons => flwd%cons(istage)

  d_bndy       => bndy%data
  bndy_bc_info => d_bndy%bc_infos  
  bc_type      =  bndy_bc_info%p(1)
  
  lo = set_shape(d_bndy%ilo(1), d_bndy%jlo(1))
  hi = set_shape(d_bndy%ihi(1), d_bndy%jhi(1))
  
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
  
  select case(bc_type)
  case(inlt11_fs_idx)  
    do j = 0, jrng - 1
      do i = 0, irng - 1
        i0 = lo%d(1) + incr1*i
        j0 = lo%d(2) + incr2*j

        rho0 = info%gamma
        u0   = 3.d0
        v0   = 0.d0
        p0   = 1.d0      
        ei0  = p0/(rho0*(info%gamma-1.d0))  

        cons%p(1,i0,j0) = rho0
        cons%p(2,i0,j0) = rho0*u0
        cons%p(3,i0,j0) = rho0*v0
        cons%p(4,i0,j0) = rho0*(ei0+0.5d0*(u0**2+v0**2))
      enddo
    enddo        
  case(rflt12_fs_idx)   
    val => d_bndy%vector(1) 
    call boundary_reflect(cons, lo, hi, val, bndy%group_id) 
  case(rflt13_fs_idx)   
    val => d_bndy%vector(1) 
    call boundary_reflect(cons, lo, hi, val, bndy%group_id) 
  case(rflt15_fs_idx)    
    val => d_bndy%vector(1) 
    call boundary_reflect(cons, lo, hi, val, bndy%group_id)
  case(rflt22_fs_idx)    
    val => d_bndy%vector(1) 
    call boundary_reflect(cons, lo, hi, val, bndy%group_id)
  case(otlt23_fs_idx)   
    axis = bndy_bc_info%p(2)
    val  => d_bndy%vector(1) 
    call boundary_neumann(cons, lo, hi, val, axis) 
  case(rflt24_fs_idx)    
    val => d_bndy%vector(1) 
    call boundary_reflect(cons, lo, hi, val, bndy%group_id)
  case default
    print*,"Error set_fs_bc():"
    print*,"No such BC type index:",bc_type, &
           "boundary index:",bndy%idx,         &
           "flwd index:", iflwd
    stop
  end select

  prim => flwd%prim
  call cmpt_cons2prim_range(prim, cons, lo, hi, info) 
end subroutine set_fs_bc

! subroutine set_ssl_bc(flwds,             &
!                       nflwd, iflwd,      &
!                       bndy,              &
!                       istage)
!   implicit none
!   integer,intent(in) :: nflwd,iflwd
!   type(blk_bndy_t),target,intent(in) :: bndy
!   integer,intent(in) :: istage
!   type(eu_data_t),target,dimension(nflwd),intent(inout) :: flwds
  
!   type(block_t),pointer        :: blk
!   type(eu_data_t),pointer      :: flwd
!   type(eu_info_data_t),pointer :: info
  
!   type(eu_bc_data_t),pointer   :: bc
!   type(ivector_t),pointer      :: bndy_bc_info
!   type(rvectorfield_t),pointer :: cons
!   type(rvectorfield_t),pointer :: prim

!   integer :: iflwd_d
!   type(block_t),pointer        :: blk_d
!   type(eu_data_t),pointer      :: flwd_d
!   type(eu_info_data_t),pointer :: info_d
!   type(rvectorfield_t),pointer :: cons_d

!   type(fscalarshape_t) :: lo,hi
!   type(fscalarshape_t) :: lo_d,hi_d
!   type(fscalarshape_t) :: lo_updt,hi_updt

!   integer :: bc_type

!   type(rvector_t),pointer :: val
!   integer :: axis
!   integer :: cnnct
!   integer :: irng,jrng,i0,j0,i,j
!   integer :: incr1,incr2

!   double precision :: xc0,yc0
!   double precision :: rho0,u0,v0,p0,ei0
!   double precision,parameter :: pie = 4.d0*atan(1.d0)
  
!   flwd => flwds(iflwd)
!   blk  => flwd%blk

!   info => flwd%info
!   bc   => flwd%bcs(bndy%idx)

!   cons => flwd%cons(istage)

!   lo = set_shape(bndy%ilo(1), bndy%jlo(1))
!   hi = set_shape(bndy%ihi(1), bndy%jhi(1))

!   bndy_bc_info => bndy%bc_infos(1)  
!   bc_type = bndy_bc_info%p(2)
!   select case(bc_type)
!   case(inlt1_ssl_idx)
!     irng = abs(hi%d(1) - lo%d(1)) + 1
!     jrng = abs(hi%d(2) - lo%d(2)) + 1    
!     if(lo%d(1).lt.hi%d(1)) then
!       incr1 = 1
!     else
!       incr1 =-1
!     endif  
!     if(lo%d(2).lt.hi%d(2)) then
!       incr2 = 1
!     else
!       incr2 =-1
!     endif  
!     do j = 0, jrng - 1
!       do i = 0, irng - 1
!         i0 = lo%d(1) + incr1*i
!         j0 = lo%d(2) + incr2*j
!         xc0 = blk%xc%p(i0,j0)
!         yc0 = blk%yc%p(i0,j0)

!         if(yc0.ge.0.d0) then
!           rho0 = 1.6374d0
!           u0   = 2.5d0+0.5d0*tanh(2.d0*yc0)
!           v0   = 0.05d0*cos(2.d0*pie*info%time/(30.d0/2.68d0))*exp(-yc0**2/10.d0) &
!                 +0.05d0*cos(2.d0*pie*2.d0*info%time/(30.d0/2.68d0)+pie/2.d0)*exp(-yc0**2/10.d0)
!           p0   = 0.3327d0
!         else
!           rho0 = 0.3626d0
!           u0   = 2.5d0+0.5d0*tanh(2.d0*yc0)
!           v0   = 0.05d0*cos(2.d0*pie*info%time/(30.d0/2.68d0))*exp(-yc0**2/10.d0) &
!                 +0.05d0*cos(2.d0*pie*2.d0*info%time/(30.d0/2.68d0)+pie/2.d0)*exp(-yc0**2/10.d0)
!           p0   = 0.3327d0
!         endif
    
!         ei0  = p0/(rho0*(info%gamma-1.d0))  
!         cons%p(1,i0,j0) = rho0
!         cons%p(2,i0,j0) = rho0*u0
!         cons%p(3,i0,j0) = rho0*v0
!         cons%p(4,i0,j0) = rho0*(ei0+0.5d0*(u0**2+v0**2))
!       enddo
!     enddo
!   case(rflt2_ssl_idx)         
!     lo_d = set_shape(bndy%ilo(2), bndy%jlo(2))
!     hi_d = set_shape(bndy%ihi(2), bndy%jhi(2))   
!     val  => bc%vs(1) 
!     call boundary_reflect(cons, lo, hi, lo_d, hi_d, val)
!   case(otlt3_ssl_idx)
!     axis = bndy_bc_info%p(3)
!     val => bc%vs(1)
!     call boundary_neumann(cons, lo, hi, val, axis)
!   case(inlt4_ssl_idx)    
!     irng = abs(hi%d(1) - lo%d(1)) + 1
!     jrng = abs(hi%d(2) - lo%d(2)) + 1    
!     if(lo%d(1).lt.hi%d(1)) then
!       incr1 = 1
!     else
!       incr1 =-1
!     endif  
!     if(lo%d(2).lt.hi%d(2)) then
!       incr2 = 1
!     else
!       incr2 =-1
!     endif  
!     do j = 0, jrng - 1
!       do i = 0, irng - 1
!         i0  = lo%d(1) + incr1*i
!         j0  = lo%d(2) + incr2*j
        
!         rho0 = 2.1101d0
!         u0   = 2.9709d0
!         v0   =-0.1367d0
!         p0   = 0.4754
        
!         ei0  = p0/(rho0*(info%gamma-1.d0))  
!         cons%p(1,i0,j0) = rho0
!         cons%p(2,i0,j0) = rho0*u0
!         cons%p(3,i0,j0) = rho0*v0
!         cons%p(4,i0,j0) = rho0*(ei0+0.5d0*(u0**2+v0**2))
!       enddo
!     enddo
!   case(intf5_ssl_idx)
!     iflwd_d = bndy_bc_info%p(3)
!     flwd_d => flwds(iflwd_d)
!     cons_d => flwd_d%cons(istage)
  
!     lo_d = set_shape(bndy%ilo(2), bndy%jlo(2))
!     hi_d = set_shape(bndy%ihi(2), bndy%jhi(2))    
    
!     cnnct = bndy_bc_info%p(4)
!     call boundary_interface(cons, cons_d,     &
!                             lo, hi, lo_d, hi_d, &
!                             cnnct)
!   case default
!     print*,"Error set_ssl_bc():"
!     print*,"No such BC type index:",bc_type, &
!            "boundary index:",bndy%idx,         &
!            "flwd index:", iflwd
!     stop
!   end select

!   prim => flwd%prim
!   call cmpt_cons2prim_range(prim, cons, lo, hi, info) 
! end subroutine set_ssl_bc

! subroutine set_st_bc(flwds,             &
!                      nflwd, iflwd,      &
!                      bndy,              &
!                      istage)
!   implicit none
!   integer,intent(in) :: nflwd,iflwd
!   type(blk_bndy_t),target,intent(in) :: bndy
!   integer,intent(in) :: istage
!   type(eu_data_t),target,dimension(nflwd),intent(inout) :: flwds
  
!   type(block_t),pointer        :: blk
!   type(eu_data_t),pointer      :: flwd
!   type(eu_info_data_t),pointer :: info
  
!   type(eu_bc_data_t),pointer   :: bc
!   type(ivector_t),pointer      :: bndy_bc_info
!   type(rvectorfield_t),pointer :: cons
!   type(rvectorfield_t),pointer :: prim

!   integer :: iflwd_d
!   type(block_t),pointer        :: blk_d
!   type(eu_data_t),pointer      :: flwd_d
!   type(eu_info_data_t),pointer :: info_d
!   type(rvectorfield_t),pointer :: cons_d

!   type(fscalarshape_t) :: lo,hi
!   type(fscalarshape_t) :: lo_d,hi_d
!   type(fscalarshape_t) :: lo_updt,hi_updt

!   integer :: bc_type

!   type(rvector_t),pointer :: val
!   integer :: axis
!   integer :: cnnct
!   integer :: irng,jrng,i0,j0,i,j
!   integer :: incr1,incr2

!   double precision :: xc0,yc0
!   double precision :: rho0,u0,v0,p0,ei0
  
!   flwd => flwds(iflwd)
!   blk  => flwd%blk

!   info => flwd%info
!   bc   => flwd%bcs(bndy%idx)

!   cons => flwd%cons(istage)

!   lo = set_shape(bndy%ilo(1), bndy%jlo(1))
!   hi = set_shape(bndy%ihi(1), bndy%jhi(1))

!   bndy_bc_info => bndy%bc_infos(1)  
!   bc_type = bndy_bc_info%p(2)
!   select case(bc_type)
!   case(otlt1_st_idx)
!     axis = bndy_bc_info%p(3)
!     val => bc%vs(1)
!     call boundary_neumann(cons, lo, hi, val, axis)
!   case(rflt2_st_idx)         
!     lo_d = set_shape(bndy%ilo(2), bndy%jlo(2))
!     hi_d = set_shape(bndy%ihi(2), bndy%jhi(2))   
!     val  => bc%vs(1) 
!     call boundary_reflect(cons, lo, hi, lo_d, hi_d, val)
!   case(otlt3_st_idx)
!     axis = bndy_bc_info%p(3)
!     val => bc%vs(1)
!     call boundary_neumann(cons, lo, hi, val, axis)
!   case(rflt4_st_idx)    
!     lo_d = set_shape(bndy%ilo(2), bndy%jlo(2))
!     hi_d = set_shape(bndy%ihi(2), bndy%jhi(2))   
!     val  => bc%vs(1) 
!     call boundary_reflect(cons, lo, hi, lo_d, hi_d, val)
!   case(intf4_st_idx)
!     iflwd_d = bndy_bc_info%p(3)
!     flwd_d => flwds(iflwd_d)
!     cons_d => flwd_d%cons(istage)
  
!     lo_d = set_shape(bndy%ilo(2), bndy%jlo(2))
!     hi_d = set_shape(bndy%ihi(2), bndy%jhi(2))    
    
!     cnnct = bndy_bc_info%p(4)
!     call boundary_interface(cons, cons_d,     &
!                             lo, hi, lo_d, hi_d, &
!                             cnnct)
!   case default
!     print*,"Error set_st_bc():"
!     print*,"No such BC type index:",bc_type, &
!            "boundary index:",bndy%idx,         &
!            "flwd index:", iflwd
!     stop
!   end select

!   prim => flwd%prim
!   call cmpt_cons2prim_range(prim, cons, lo, hi, info) 
! end subroutine set_st_bc

! subroutine set_blst_bc(flwds,             &
!   nflwd, iflwd,      &
!   bndy,              &
!   istage)
!   implicit none
!   integer,intent(in) :: nflwd,iflwd
!   type(blk_bndy_t),target,intent(in) :: bndy
!   integer,intent(in) :: istage
!   type(eu_data_t),target,dimension(nflwd),intent(inout) :: flwds

!   type(block_t),pointer        :: blk
!   type(eu_data_t),pointer      :: flwd
!   type(eu_info_data_t),pointer :: info

!   type(eu_bc_data_t),pointer   :: bc
!   type(ivector_t),pointer      :: bndy_bc_info
!   type(rvectorfield_t),pointer :: cons
!   type(rvectorfield_t),pointer :: prim

!   integer :: iflwd_d
!   type(block_t),pointer        :: blk_d
!   type(eu_data_t),pointer      :: flwd_d
!   type(eu_info_data_t),pointer :: info_d
!   type(rvectorfield_t),pointer :: cons_d

!   type(fscalarshape_t) :: lo,hi
!   type(fscalarshape_t) :: lo_d,hi_d
!   type(fscalarshape_t) :: lo_updt,hi_updt

!   integer :: bc_type

!   type(rvector_t),pointer :: val
!   integer :: axis
!   integer :: cnnct
!   integer :: irng,jrng,i0,j0,i,j
!   integer :: incr1,incr2

!   double precision :: xc0,yc0
!   double precision :: rho0,u0,v0,p0,ei0

!   flwd => flwds(iflwd)
!   blk  => flwd%blk

!   info => flwd%info
!   bc   => flwd%bcs(bndy%idx)

!   cons => flwd%cons(istage)

!   lo = set_shape(bndy%ilo(1), bndy%jlo(1))
!   hi = set_shape(bndy%ihi(1), bndy%jhi(1))

!   bndy_bc_info => bndy%bc_infos(1)  
!   bc_type = bndy_bc_info%p(2)
!   select case(bc_type)
!     case(rflt1_blst_idx)
!       lo_d = set_shape(bndy%ilo(2), bndy%jlo(2))
!       hi_d = set_shape(bndy%ihi(2), bndy%jhi(2))   
!       val  => bc%vs(1) 
!       call boundary_reflect(cons, lo, hi, lo_d, hi_d, val)
!     case(rflt2_blst_idx)         
!       lo_d = set_shape(bndy%ilo(2), bndy%jlo(2))
!       hi_d = set_shape(bndy%ihi(2), bndy%jhi(2))   
!       val  => bc%vs(1) 
!       call boundary_reflect(cons, lo, hi, lo_d, hi_d, val)
!     case(rflt3_blst_idx)
!       lo_d = set_shape(bndy%ilo(2), bndy%jlo(2))
!       hi_d = set_shape(bndy%ihi(2), bndy%jhi(2))   
!       val  => bc%vs(1) 
!       call boundary_reflect(cons, lo, hi, lo_d, hi_d, val)
!     case(rflt4_blst_idx)    
!       lo_d = set_shape(bndy%ilo(2), bndy%jlo(2))
!       hi_d = set_shape(bndy%ihi(2), bndy%jhi(2))   
!       val  => bc%vs(1) 
!       call boundary_reflect(cons, lo, hi, lo_d, hi_d, val)
!     case(intf5_blst_idx)
!       iflwd_d = bndy_bc_info%p(3)
!       flwd_d => flwds(iflwd_d)
!       cons_d => flwd_d%cons(istage)

!       lo_d = set_shape(bndy%ilo(2), bndy%jlo(2))
!       hi_d = set_shape(bndy%ihi(2), bndy%jhi(2))    

!       cnnct = bndy_bc_info%p(4)
!       call boundary_interface(cons, cons_d,     &
!               lo, hi, lo_d, hi_d, &
!               cnnct)
!     case default
!       print*,"Error set_blst_bc():"
!       print*,"No such BC type index:",bc_type, &
!       "boundary index:",bndy%idx,         &
!       "flwd index:", iflwd
!       stop
!   end select

!   prim => flwd%prim
!   call cmpt_cons2prim_range(prim, cons, lo, hi, info) 
! end subroutine set_blst_bc

subroutine set_intf_bc(flwds, nflwd, iflwd, bndy)
  implicit none
  integer,intent(in) :: nflwd
  type(eu_data_iterator_t),target,dimension(nflwd),intent(inout) :: flwds
  integer,intent(in) :: iflwd
  type(bndy_t),target,intent(in) :: bndy
  type(bndy_data_t),pointer      :: d_bndy
  type(ivector_t),pointer        :: bndy_bc_info

  type(eu_data_t),pointer      :: flwd
  type(eu_info_data_t),pointer :: info
  type(rvectorfield_t),pointer :: cons
  type(rvectorfield_t),pointer :: prim

  integer :: iflwd_d
  type(eu_data_t),pointer      :: flwd_d
  type(eu_info_data_t),pointer :: info_d
  type(rvectorfield_t),pointer :: cons_d

  type(fscalarshape_t) :: lo,hi
  type(fscalarshape_t) :: lo_d,hi_d

  integer :: istage

  flwd => flwds(iflwd)%flwd
  info => flwd%info

  if(info%stage_cur.eq.1) then
  istage = info%stage_num
  else
  istage = info%stage_cur - 1
  endif  

  cons => flwd%cons(istage)

  d_bndy       => bndy%data
  bndy_bc_info => d_bndy%bc_infos  
  lo = set_shape(d_bndy%ilo(1), d_bndy%jlo(1))
  hi = set_shape(d_bndy%ihi(1), d_bndy%jhi(1))

  iflwd_d = bndy_bc_info%p(2)
  flwd_d  => flwds(iflwd_d)%flwd
  cons_d  => flwd_d%cons(istage)

  lo_d = set_shape(d_bndy%ilo(2), d_bndy%jlo(2))
  hi_d = set_shape(d_bndy%ihi(2), d_bndy%jhi(2))  

  call boundary_interface(cons, cons_d,     &
        lo, hi,           &
        lo_d, hi_d,       &
        1)

  prim => flwd%prim
  call cmpt_cons2prim_range(prim, cons, lo, hi, info) 
end subroutine set_intf_bc
end module eu_boundary