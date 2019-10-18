module eu_solver
  use data_type
  use eu_data_type
  use reconstruction  
  use time_advance
  use mesh
  implicit none
  ! Riemann solver indices
  integer,parameter :: flux_llf_idx     = 100
  integer,parameter :: flux_glf_idx     = 200
  integer,parameter :: flux_sw_idx      = 300
  integer,parameter :: flux_vl_idx      = 400
  integer,parameter :: flux_ausmfvs_idx = 500
  integer,parameter :: flux_roe_idx     = 600
  integer,parameter :: flux_hllc_idx    = 700

  integer,parameter :: pres_measure_idx = 1
  integer,parameter :: ent_measure_idx  = 2
  integer,parameter :: g_measure_idx    = 3
  integer,parameter :: sw_measure_idx   = 4  

contains
subroutine cmpt_cnvt_eu(flwd, &
                        riemann_solver, &
                        reconst_xi_idx, &
                        reconst_eta_idx)
  implicit none
  integer,intent(in) :: riemann_solver, &
                        reconst_xi_idx, &
                        reconst_eta_idx
  type(eu_data_t),target,intent(inout) :: flwd
  type(eu_info_data_t),pointer :: info
  type(eu_convect_data_t),pointer :: cnvt
  
  select case(riemann_solver)
    case(flux_llf_idx)
      call cmpt_flux_llf_eu(flwd)
      if(mod(riemann_solver,100).eq.0) then
        call cmpt_rnst_component_eu(flwd,             &
                                    reconst_xi_idx,   &
                                    reconst_eta_idx)
      else 
        call cmpt_rnst_commweno_eu(flwd,                     &
                                   mod(riemann_solver,100),  &
                                   reconst_xi_idx,           &
                                   reconst_eta_idx)
      endif
      call cmpt_cnvt_diff(flwd) 
    case(flux_glf_idx)
      call cmpt_flux_glf_eu(flwd)
      if(mod(riemann_solver,100).eq.0) then
        call cmpt_rnst_component_eu(flwd,             &
                                    reconst_xi_idx,   &
                                    reconst_eta_idx)
      else 
        call cmpt_rnst_commweno_eu(flwd,                     &
                                   mod(riemann_solver,100),  &
                                   reconst_xi_idx,           &
                                   reconst_eta_idx)
      endif
      call cmpt_cnvt_diff(flwd) 
    case(flux_sw_idx)
      call cmpt_flux_sw_eu(flwd)
      if(mod(riemann_solver,100).eq.0) then
        call cmpt_rnst_component_eu(flwd,             &
                                    reconst_xi_idx,   &
                                    reconst_eta_idx)
      else 
        call cmpt_rnst_commweno_eu(flwd,                     &
                                   mod(riemann_solver,100),  &
                                   reconst_xi_idx,           &
                                   reconst_eta_idx)
      endif
      call cmpt_cnvt_diff(flwd) 
    case(flux_vl_idx)
      call cmpt_flux_vl_eu(flwd)
      if(mod(riemann_solver,100).eq.0) then
        call cmpt_rnst_component_eu(flwd,             &
                                    reconst_xi_idx,   &
                                    reconst_eta_idx)
      else 
        call cmpt_rnst_commweno_eu(flwd,                     &
                                   mod(riemann_solver,100),  &
                                   reconst_xi_idx,           &
                                   reconst_eta_idx)
      endif
      call cmpt_cnvt_diff(flwd) 
    case(flux_ausmfvs_idx)
      call cmpt_flux_ausmfvs_eu(flwd)
      if(mod(riemann_solver,100).eq.0) then
        call cmpt_rnst_component_eu(flwd,             &
                                    reconst_xi_idx,   &
                                    reconst_eta_idx)
      else 
        call cmpt_rnst_commweno_eu(flwd,                     &
                                   mod(riemann_solver,100),  &
                                   reconst_xi_idx,           &
                                   reconst_eta_idx)
      endif
      call cmpt_cnvt_diff(flwd) 
    case(flux_hllc_idx)
      call cmpt_rnst_prim_eu(flwd,             &
                             reconst_xi_idx,   &
                             reconst_eta_idx)
      ! call cmpt_flux_hllc_eu(flwd,istage,     &
      !                        reconst_xi_idx,  &
      !                        reconst_eta_idx)
    case(flux_roe_idx)
      ! call cmpt_flux_roe_eu(flwd,istage,     &
      !                        reconst_xi_idx,  &
      !                        reconst_eta_idx)
    case default
      write(*,*) "Error: cmpt_cnvt()"
      write(*,*) "No such Riemann Solver Index:", riemann_solver
      stop
  end select
end subroutine cmpt_cnvt_eu

subroutine cmpt_dfsv_eu(flwd,        &
                        grad_xi_idx, &
                        grad_eta_idx)
  implicit none
  integer,intent(in) :: grad_xi_idx, &
                        grad_eta_idx
  type(eu_data_t),target,intent(inout) :: flwd

end subroutine cmpt_dfsv_eu

subroutine cmpt_src_eu(flwd)
  implicit none
  type(eu_data_t),target,intent(inout) :: flwd
  
end subroutine cmpt_src_eu

subroutine cmpt_rhs_eu(flwd)
  implicit none
  type(eu_data_t),target,intent(inout) :: flwd

  type(eu_info_data_t),pointer :: info
  type(rvectorfield_t),pointer :: cnvt
  type(rvectorfield_t),pointer :: dfsv
  type(rvectorfield_t),pointer :: src
  type(rvectorfield_t),pointer :: rhs

  integer :: isize,jsize,ibf,jbf
  integer :: i,j
  integer :: i0,j0

  info => flwd%info

  isize = info%isize
  jsize = info%jsize
  ibf   = info%ibf
  jbf   = info%jbf

  rhs   => flwd%rhs(info%stage_cur)  
  rhs%p = 0.d0  

  cnvt => flwd%cnvt%cnvt  
  do j = 1, jsize
    do i = 1, isize
      i0 = i + ibf
      j0 = j + jbf
      rhs%p(:,i0,j0) = cnvt%p(:,i0,j0)
    enddo
  enddo
end subroutine cmpt_rhs_eu 

subroutine cmpt_lhs_eu(flwd)
  implicit none
  type(eu_data_t),intent(inout) :: flwd
end subroutine cmpt_lhs_eu 

subroutine cmpt_dt_eu(flwd)
  implicit none
  type(eu_data_t),target,intent(inout) :: flwd

  type(eu_info_data_t),pointer :: info
  double precision :: dxi,deta
  double precision :: dt_xi,dt_eta
  info    => flwd%info
  dxi     = info%dxi
  deta    = info%deta
  dt_xi   = dxi/info%eig_xi_g
  dt_eta  = deta/info%eig_eta_g
  info%dt = min(dt_xi,dt_eta)
end subroutine cmpt_dt_eu  

subroutine evlv_mdl_eu(flwd, time_integrator_idx)
  implicit none 
  integer,intent(in) :: time_integrator_idx
  type(eu_data_t),target,intent(inout) :: flwd
  type(eu_info_data_t),pointer :: info
  type(fvectorshape_t) :: hi,lo
  integer :: isize,jsize,ibf,jbf
  integer :: cons_num
  info     => flwd%info
  cons_num = info%cons_num

  isize = info%isize
  jsize = info%jsize
  ibf   = info%ibf
  jbf   = info%jbf

  lo = set_shape(1, 1 + ibf, 1 + jbf)
  hi = set_shape(cons_num,isize+ibf,jsize+jbf)  
  call cmpt_time_advn_stg(flwd%cons,         &
                          flwd%lhs,          &
                          flwd%rhs,          &
                          info%dt,           &
                          lo,hi,             &
                          info%stage_cur,    &
                          info%stage_num,    &
                          time_integrator_idx)
end subroutine evlv_mdl_eu

subroutine cmpt_rnst_component_eu(flwd,           &
                                  reconst_xi_idx, &
                                  reconst_eta_idx)
  implicit none
  integer,intent(in) :: reconst_xi_idx, &
                        reconst_eta_idx
  type(eu_data_t),target,intent(inout) :: flwd

  type(eu_info_data_t),pointer    :: info
  type(eu_convect_data_t),pointer :: cnvt

  type(rvectorfield_t),pointer :: flux
  type(rvectorfield_t),pointer :: h
  type(rvectorfield_t),pointer :: diff
  type(rvectorfield_t),pointer :: diffp,diffn
  type(rvectorfield_t),pointer :: diff_xi,diff_eta

  double precision :: invdxi,invdeta
  integer :: isize,jsize,ibf,jbf
  integer :: i,j
  integer :: i0,j0

  type(fvectorshape_t) :: lo_xi,hi_xi,padding_xi
  type(fvectorshape_t) :: lo_eta,hi_eta,padding_eta
  integer :: cnvt_num

  type(rvectorfield_t) :: measure
  type(fvectorshape_t) :: measure_s

  info => flwd%info
  cnvt => flwd%cnvt
  
  isize    = info%isize
  jsize    = info%jsize
  ibf      = info%ibf
  jbf      = info%jbf
  cnvt_num = info%cnvt_num
  
  ! Xi direction 
  h     => cnvt%h_xip
  flux  => cnvt%flux_xip
  lo_xi = set_shape(1, 1, 1)
  hi_xi = set_shape(cnvt_num, &
                    h%s%d(2), &
                    h%s%d(3))
  padding_xi = set_shape(0, ibf-1, 0)
  call reconst(h,           &
               flux,        &
               lo_xi,       &
               hi_xi,       &
               padding_xi,  &
               bias_upwind, &
               axis_xi,     &
               reconst_xi_idx)

  h     => cnvt%h_xin
  flux  => cnvt%flux_xin
  lo_xi = set_shape(1, 1, 1)
  hi_xi = set_shape(cnvt_num, &
                    h%s%d(2), &
                    h%s%d(3))
  padding_xi = set_shape(0, ibf-1, 0)
  call reconst(h,             &
               flux,          &
               lo_xi,         &
               hi_xi,         &
               padding_xi,    &
               bias_downwind, &
               axis_xi,       &
               reconst_xi_idx)

  ! Eta direction 
  h      => cnvt%h_etap
  flux   => cnvt%flux_etap
  lo_eta = set_shape(1, 1, 1)
  hi_eta = set_shape(cnvt_num, &
                      h%s%d(2), &
                      h%s%d(3))
  padding_eta = set_shape(0, jbf-1, 0)
  call reconst(h,           &
               flux,        &
               lo_eta,      &
               hi_eta,      &
               padding_eta, &
               bias_upwind, &
               axis_xi,     &
               reconst_eta_idx)

  h      => cnvt%h_etan
  flux   => cnvt%flux_etan
  lo_eta = set_shape(1, 1, 1)
  hi_eta = set_shape(cnvt_num, &
                      h%s%d(2), &
                      h%s%d(3))
  padding_eta = set_shape(0, jbf-1, 0)
  call reconst(h,             &
               flux,          &
               lo_eta,        &
               hi_eta,        &
               padding_eta,   &
               bias_downwind, &
               axis_xi,       &
               reconst_eta_idx)  
end subroutine cmpt_rnst_component_eu

subroutine cmpt_rnst_commweno_eu(flwd,           &
                                 measure_idx,    &
                                 reconst_xi_idx, &
                                 reconst_eta_idx)
  implicit none
  integer,intent(in) :: measure_idx,    &
                        reconst_xi_idx, &
                        reconst_eta_idx
  type(eu_data_t),target,intent(inout) :: flwd

  ! type(eu_info_data_t),pointer    :: info
  ! type(eu_convect_data_t),pointer :: cnvt

  ! type(rvectorfield_t),pointer :: prim
  ! type(rvectorfield_t),pointer :: flux
  ! type(rvectorfield_t),pointer :: h
  ! type(rvectorfield_t),pointer :: diff
  ! type(rvectorfield_t),pointer :: diffp,diffn
  ! type(rvectorfield_t),pointer :: diff_xi,diff_eta
  
  ! double precision :: invdxi,invdeta
  ! integer :: isize,jsize,ibf,jbf
  ! integer :: i,j
  ! integer :: i0,j0

  ! type(fvectorshape_t) :: lo_xi,hi_xi,padding_xi
  ! type(fvectorshape_t) :: lo_eta,hi_eta,padding_eta
  ! integer :: cnvt_num

  ! type(rscalarfield_t) :: measure
  ! type(ff) :: measure_s

  ! info => flwd%info
  ! prim => flwd%prim
  ! cnvt => flwd%cnvt
  
  ! isize    = info%isize
  ! jsize    = info%jsize
  ! ibf      = info%ibf
  ! jbf      = info%jbf
  ! cnvt_num = info%cnvt_num
  
  ! ! Xi direction 
  ! h     => cnvt%h_xip
  ! flux  => cnvt%flux_xip
  ! lo_xi = set_shape(1, 1, 1)
  ! hi_xi = set_shape(cnvt_num, &
  !                   h%s%d(2), &
  !                   h%s%d(3))
  ! padding_xi = set_shape(0, ibf-1, 0)
  ! measure_s  = set_shape(flux%s%d(2),flux%s%d(3))
  ! call create(measure, measure_s)
  ! call cmpt_measure(measure, prim, measure_idx, lo_xi, hi_xi, axis_xi)
  ! call reconst(h,           &
  !              flux,        &
  !              measure,     &
  !              lo_xi,       &
  !              hi_xi,       &
  !              padding_xi,  &
  !              bias_upwind, &
  !              axis_xi,     &
  !              reconst_xi_idx)
  ! call destroy(measure)
  ! h     => cnvt%h_xin
  ! flux  => cnvt%flux_xin
  ! lo_xi = set_shape(1, 1, 1)
  ! hi_xi = set_shape(cnvt_num, &
  !                   h%s%d(2), &
  !                   h%s%d(3))
  ! padding_xi = set_shape(0, ibf-1, 0)
  ! measure_s  = set_shape(flux%s%d(2),flux%s%d(3))
  ! call create(measure, measure_s)
  ! call cmpt_measure(measure, prim, measure_idx, lo_xi, hi_xi, axis_xi)
  ! call reconst(h,             &
  !              flux,          &
  !              measure,       &
  !              lo_xi,         &
  !              hi_xi,         &
  !              padding_xi,    &
  !              bias_downwind, &
  !              axis_xi,       &
  !              reconst_xi_idx)
  ! call destroy(measure)

  ! ! Eta direction 
  ! h      => cnvt%h_etap
  ! flux   => cnvt%flux_etap
  ! lo_eta = set_shape(1, 1, 1)
  ! hi_eta = set_shape(cnvt_num, &
  !                     h%s%d(2), &
  !                     h%s%d(3))
  ! padding_eta = set_shape(0, jbf-1, 0)
  ! measure_s  = set_shape(flux%s%d(2),flux%s%d(3))
  ! call create(measure, measure_s)
  ! call cmpt_measure(measure, prim, measure_idx, lo_xi, hi_xi, axis_xi)
  ! call reconst(h,           &
  !              flux,        &
  !              measure,     &
  !              lo_eta,      &
  !              hi_eta,      &
  !              padding_eta, &
  !              bias_upwind, &
  !              axis_xi,     &
  !              reconst_eta_idx)
  ! call destroy(measure)

  ! h      => cnvt%h_etan
  ! flux   => cnvt%flux_etan
  ! lo_eta = set_shape(1, 1, 1)
  ! hi_eta = set_shape(cnvt_num, &
  !                     h%s%d(2), &
  !                     h%s%d(3))
  ! padding_eta = set_shape(0, jbf-1, 0)
  ! measure_s  = set_shape(flux%s%d(2),flux%s%d(3))
  ! call create(measure, measure_s)
  ! call cmpt_measure(measure, prim, measure_idx, lo_xi, hi_xi, axis_xi)
  ! call reconst(h,             &
  !              flux,          &
  !              measure,       &
  !              lo_eta,        &
  !              hi_eta,        &
  !              padding_eta,   &
  !              bias_downwind, &
  !              axis_xi,       &
  !              reconst_eta_idx)  
  ! call destroy(measure)
end subroutine cmpt_rnst_commweno_eu

subroutine cmpt_rnst_prim_eu(flwd,           &
                             reconst_xi_idx, &
                             reconst_eta_idx)
  implicit none
  integer,intent(in) :: reconst_xi_idx, &
                        reconst_eta_idx
  type(eu_data_t),target,intent(inout) :: flwd

  type(eu_info_data_t),pointer    :: info
  type(eu_convect_data_t),pointer :: cnvt

  type(rvectorfield_t),pointer :: prim
  type(rvectorfield_t),pointer :: prim_l,prim_r

  integer :: ibf,jbf

  type(fvectorshape_t) :: lo_xi,hi_xi,padding_xi
  type(fvectorshape_t) :: lo_eta,hi_eta,padding_eta

  info => flwd%info
  cnvt => flwd%cnvt

  ibf      = info%ibf
  jbf      = info%jbf

  ! Xi direction 
  prim   => cnvt%flux_xi
  prim_l => cnvt%h_xip
  lo_xi = set_shape(1, 1, 1)
  hi_xi = set_shape(4, &
                    prim_l%s%d(2), &
                    prim_l%s%d(3))
  padding_xi = set_shape(0, ibf-1, 0)
  call reconst(prim_l,      &
               prim,        &
               lo_xi,       &
               hi_xi,       &
               padding_xi,  &
               bias_upwind, &
               axis_xi,     &
               reconst_xi_idx)

  prim_r => cnvt%h_xin
  lo_xi = set_shape(1, 1, 1)
  hi_xi = set_shape(4, &
                    prim_r%s%d(2), &
                    prim_r%s%d(3))
  padding_xi = set_shape(0, ibf-1, 0)
  call reconst(prim_r,        &
               prim,          &
               lo_xi,         &
               hi_xi,         &
               padding_xi,    &
               bias_downwind, &
               axis_xi,       &
               reconst_xi_idx)

  ! Eta direction 
  prim   => cnvt%flux_eta
  prim_l => cnvt%h_etap
  lo_eta = set_shape(1, 1, 1)
  hi_eta = set_shape(4, prim_l%s%d(2), prim_l%s%d(3))
  padding_eta = set_shape(0, jbf-1, 0)
  call reconst(prim_l,           &
               prim,        &
               lo_eta,      &
               hi_eta,      &
               padding_eta, &
               bias_upwind, &
               axis_xi,     &
               reconst_eta_idx)

  prim_r => cnvt%h_etap
  lo_eta = set_shape(1, 1, 1)
  hi_eta = set_shape(4, prim_r%s%d(2), prim_r%s%d(3))
  padding_eta = set_shape(0, jbf-1, 0)
  call reconst(prim_r,             &
               prim,          &
               lo_eta,        &
               hi_eta,        &
               padding_eta,   &
               bias_downwind, &
               axis_xi,       &
               reconst_eta_idx)  
end subroutine cmpt_rnst_prim_eu

subroutine cmpt_flux_llf_eu(flwd)
  implicit none
  type(eu_data_t),target,intent(inout) :: flwd
  
  type(eu_info_data_t),pointer :: info

  type(rvectorfield_t),pointer :: cons
  type(rvectorfield_t),pointer :: prim

  type(eu_convect_data_t),pointer :: cnvt 
  type(rvectorfield_t),pointer    :: flux
  type(rvectorfield_t),pointer    :: fluxp,fluxn
  
  double precision,dimension(:),pointer :: cons_elem
  double precision,dimension(:),pointer :: prim_elem
  double precision,dimension(:),pointer :: flux_elem
  double precision,dimension(:),pointer :: fluxp_elem
  double precision,dimension(:),pointer :: fluxn_elem

  integer :: ibf,jbf
  integer :: i,j
  integer :: i0,j0
  double precision :: u,v,p
  double precision :: max_eig
  integer :: pre_stage
  
  info => flwd%info
  if(info%stage_cur.eq.1) then
    pre_stage = info%stage_num
  else
    pre_stage = info%stage_cur - 1
  endif 
  ibf       = info%ibf
  jbf       = info%jbf

  prim => flwd%prim
  cons => flwd%cons(pre_stage)
  cnvt => flwd%cnvt

  ! xi direction
  flux => cnvt%flux_xi
  do j=1,flux%s%d(3)
    do i=1,flux%s%d(2)
      i0 = i
      j0 = j + jbf
      prim_elem=>prim%p(:,i0,j0)
      cons_elem=>cons%p(:,i0,j0)
      flux_elem=>flux%p(:,i,j)
      u = prim_elem(2)
      p = prim_elem(4)
      flux_elem(:) = cons_elem(:)*u
      flux_elem(2) = flux_elem(2) + p
      flux_elem(4) = flux_elem(4) + u*p
    enddo
  enddo  

  fluxp => cnvt%flux_xip
  fluxn => cnvt%flux_xin
  do j=1,flux%s%d(3)
    do i=1,flux%s%d(2)
      i0 = i
      j0 = j + jbf
      prim_elem  => prim%p(:,i0,j0)
      cons_elem  => cons%p(:,i0,j0)
      flux_elem  => flux%p(:,i,j)
      fluxp_elem => fluxp%p(:,i,j)
      fluxn_elem => fluxn%p(:,i,j)

      max_eig = abs(prim_elem(2)) + prim_elem(5)
      fluxp_elem(:) = 0.5d0*(flux_elem(:) + max_eig*cons_elem(:))
      fluxn_elem(:) = 0.5d0*(flux_elem(:) - max_eig*cons_elem(:))
    enddo
  enddo

  ! j direction
  flux => cnvt%flux_eta
  do j=1,flux%s%d(3)
    do i=1,flux%s%d(2)
      i0 = j + ibf
      j0 = i
      prim_elem => prim%p(:,i0,j0)
      cons_elem => cons%p(:,i0,j0)
      flux_elem => flux%p(:,i,j)
      v = prim_elem(3)
      p = prim_elem(4)
      flux_elem(:) = cons_elem(:)*v
      flux_elem(3) = flux_elem(3) + p
      flux_elem(4) = flux_elem(4) + v*p
    enddo
  enddo  

  fluxp => cnvt%flux_etap
  fluxn => cnvt%flux_etan
  do j=1,flux%s%d(3)
    do i=1,flux%s%d(2)
      i0 = j + ibf
      j0 = i
      cons_elem  => cons%p(:,i0,j0)
      prim_elem  => prim%p(:,i0,j0)
      flux_elem  => flux%p(:,i,j)
      fluxp_elem => fluxp%p(:,i,j)
      fluxn_elem => fluxn%p(:,i,j)
      max_eig = abs(prim_elem(3)) + prim_elem(5)
      fluxp_elem(:) = 0.5d0*(flux_elem(:) + max_eig*cons_elem(:))
      fluxn_elem(:) = 0.5d0*(flux_elem(:) - max_eig*cons_elem(:))
    enddo
  enddo
end subroutine cmpt_flux_llf_eu

subroutine cmpt_flux_glf_eu(flwd)
  implicit none
  type(eu_data_t),target,intent(inout) :: flwd

  type(eu_info_data_t),pointer :: info

  type(rvectorfield_t),pointer :: cons
  type(rvectorfield_t),pointer :: prim

  type(eu_convect_data_t),pointer :: cnvt 
  type(rvectorfield_t),pointer    :: flux
  type(rvectorfield_t),pointer    :: fluxp,fluxn
  
  double precision,dimension(:),pointer :: flux_elem
  double precision,dimension(:),pointer :: cons_elem
  double precision,dimension(:),pointer :: prim_elem

  integer :: ibf,jbf
  integer :: i,j
  integer :: i0,j0
  double precision :: u,v,p
  integer :: pre_stage
  
  info=>flwd%info
  if(info%stage_cur.eq.1) then
    pre_stage = info%stage_num
  else
    pre_stage = info%stage_cur - 1
  endif 
  ibf       = info%ibf
  jbf       = info%jbf

  prim => flwd%prim
  cons => flwd%cons(pre_stage)
  cnvt => flwd%cnvt

  ! xi direction
  flux => cnvt%flux_xi
  do j=1,flux%s%d(3)
    do i=1,flux%s%d(2)
      i0 = i
      j0 = j + jbf
      prim_elem=>prim%p(:,i0,j0)
      cons_elem=>cons%p(:,i0,j0)
      flux_elem=>flux%p(:,i,j)
      u = prim_elem(2)
      p = prim_elem(4)
      flux_elem(:) = cons_elem(:)*u
      flux_elem(2) = flux_elem(2) + p
      flux_elem(4) = flux_elem(4) + u*p
    enddo
  enddo  

  fluxp => cnvt%flux_xip
  fluxn => cnvt%flux_xin
  do j=1,flux%s%d(3)
    do i=1,flux%s%d(2)
      i0 = i
      j0 = j + jbf
      cons_elem => cons%p(:,i0,j0)
      flux_elem => flux%p(:,i,j)
      fluxp%p(:,i,j) = 0.5d0*(flux_elem(:) + info%eig_xi_g*cons_elem(:))
      fluxn%p(:,i,j) = 0.5d0*(flux_elem(:) - info%eig_xi_g*cons_elem(:))
    enddo
  enddo

  ! j direction
  flux => cnvt%flux_eta
  do j=1,flux%s%d(3)
    do i=1,flux%s%d(2)
      i0 = j + ibf
      j0 = i
      prim_elem => prim%p(:,i0,j0)
      cons_elem => cons%p(:,i0,j0)
      flux_elem => flux%p(:,i,j)
      v = prim_elem(3)
      p = prim_elem(4)
      flux_elem(:) = cons_elem(:)*v
      flux_elem(3) = flux_elem(3) + p
      flux_elem(4) = flux_elem(4) + v*p
    enddo
  enddo  

  fluxp => cnvt%flux_etap
  fluxn => cnvt%flux_etan
  do j=1,flux%s%d(3)
    do i=1,flux%s%d(2)
      i0 = j + ibf
      j0 = i
      cons_elem => cons%p(:,i0,j0)
      flux_elem => flux%p(:,i,j)
      fluxp%p(:,i,j) = 0.5d0*(flux_elem(:) + info%eig_eta_g*cons_elem(:))
      fluxn%p(:,i,j) = 0.5d0*(flux_elem(:) - info%eig_eta_g*cons_elem(:))
    enddo
  enddo
end subroutine cmpt_flux_glf_eu

subroutine cmpt_flux_sw_eu(flwd)
  implicit none
  type(eu_data_t),target,intent(inout) :: flwd

  type(eu_info_data_t),pointer :: info

  type(rvectorfield_t),pointer :: cons
  type(rvectorfield_t),pointer :: prim

  type(eu_convect_data_t),pointer :: cnvt 
  type(rvectorfield_t),pointer    :: flux
  type(rvectorfield_t),pointer    :: fluxp,fluxn
  
  double precision,dimension(:),pointer :: flux_elem
  double precision,dimension(:),pointer :: cons_elem
  double precision,dimension(:),pointer :: prim_elem
  double precision,dimension(:),pointer :: fluxp_elem,fluxn_elem

  integer :: ibf,jbf
  integer :: i,j
  integer :: i0,j0
  integer :: pre_stage

  double precision :: rho,u,v,e,p,alpha
  double precision :: uu2,wii
  double precision :: lambda(3),lambda_pm(3)
  double precision :: gamma 
  double precision,parameter :: epsl = 0.125d0
  
  info => flwd%info
  if(info%stage_cur.eq.1) then
    pre_stage = info%stage_num
  else
    pre_stage = info%stage_cur - 1
  endif 

  ibf   = info%ibf
  jbf   = info%jbf
  gamma = info%gamma 
  prim  => flwd%prim
  cons  => flwd%cons(pre_stage)
  cnvt  => flwd%cnvt
  
  ! i direction
  flux  => cnvt%flux_xi
  fluxp => cnvt%flux_xip
  fluxn => cnvt%flux_xin
  do j=1,flux%s%d(3)
    do i=1,flux%s%d(2)
      i0 = i
      j0 = j + jbf

      prim_elem  => prim%p(:,i0,j0)
      cons_elem  => cons%p(:,i0,j0)
      fluxp_elem => fluxp%p(:,i,j)
      fluxn_elem => fluxn%p(:,i,j)

      rho   = prim_elem(1)
      u     = prim_elem(2)
      v     = prim_elem(3)
      p     = prim_elem(4)
      alpha = prim_elem(5) 

      uu2       = u**2+v**2
      lambda(1) = u
      lambda(2) = u-alpha
      lambda(3) = u+alpha

      lambda_pm(1) = (lambda(1) + sqrt(lambda(1)**2+epsl**2))/2.0
      lambda_pm(2) = (lambda(2) + sqrt(lambda(2)**2+epsl**2))/2.0
      lambda_pm(3) = (lambda(3) + sqrt(lambda(3)**2+epsl**2))/2.0

      wii=(3.d0-gamma)*(lambda_pm(2)+lambda_pm(3))*alpha**2/(2.d0*(gamma-1.d0))

      fluxp_elem(1) =(lambda_pm(1)*2.d0*(gamma-1.d0) &
                     +lambda_pm(2) &
                     +lambda_pm(3))*rho/(2.d0*gamma)

      fluxp_elem(2) =(lambda_pm(1)*2.d0*(gamma-1.d0)*u &
                     +lambda_pm(2)*(u-alpha) &
                     +lambda_pm(3)*(u+alpha))*rho/(2.d0*gamma)

      fluxp_elem(3) =(lambda_pm(1)*2.d0*(gamma-1.d0)*v &
                     +lambda_pm(2)*v &
                     +lambda_pm(3)*v)*rho/(2.d0*gamma)

      fluxp_elem(4) =(lambda_pm(1)*(gamma-1.d0)*uu2 &
                     +lambda_pm(2)*0.5d0*((u-alpha)**2+v**2) &
                     +lambda_pm(3)*0.5d0*((u+alpha)**2+v**2) &
                     +wii)*rho/(2.d0*gamma)
     
      lambda_pm(1) = (lambda(1)-sqrt(lambda(1)**2+epsl**2))/2.0
      lambda_pm(2) = (lambda(2)-sqrt(lambda(2)**2+epsl**2))/2.0
      lambda_pm(3) = (lambda(3)-sqrt(lambda(3)**2+epsl**2))/2.0

      wii=(3.d0-gamma)*(lambda_pm(2)+lambda_pm(3))*alpha**2/(2.d0*(gamma-1.d0))

      fluxn_elem(1) =(lambda_pm(1)*2.d0*(gamma-1.d0) &
                     +lambda_pm(2) &
                     +lambda_pm(3))*rho/(2.d0*gamma)

      fluxn_elem(2) =(lambda_pm(1)*2.d0*(gamma-1.d0)*u &
                     +lambda_pm(2)*(u-alpha) &
                     +lambda_pm(3)*(u+alpha))*rho/(2.d0*gamma)

      fluxn_elem(3) =(lambda_pm(1)*2.d0*(gamma-1.d0)*v &
                     +lambda_pm(2)*v &
                     +lambda_pm(3)*v)*rho/(2.d0*gamma)

      fluxn_elem(4) =(lambda_pm(1)*(gamma-1.d0)*uu2 &
                     +lambda_pm(2)*0.5d0*((u-alpha)**2+v**2) &
                     +lambda_pm(3)*0.5d0*((u+alpha)**2+v**2) &
                     +wii)*rho/(2.d0*gamma)
     
    enddo
  enddo  

  ! j direction
  flux  => cnvt%flux_eta
  fluxp => cnvt%flux_etap
  fluxn => cnvt%flux_etan
  do j=1,flux%s%d(3)
    do i=1,flux%s%d(2)
      i0 = j + ibf
      j0 = i 

      prim_elem  => prim%p(:,i0,j0)
      cons_elem  => cons%p(:,i0,j0)
      fluxp_elem => fluxp%p(:,i,j)
      fluxn_elem => fluxn%p(:,i,j)

      rho   = prim_elem(1)
      u     = prim_elem(2)
      v     = prim_elem(3)
      p     = prim_elem(4)
      alpha = prim_elem(5)

      uu2       = u**2+v**2
      lambda(1) = v
      lambda(2) = v - alpha
      lambda(3) = v + alpha

      lambda_pm(1) = (lambda(1) + sqrt(lambda(1)**2+epsl**2))/2.0
      lambda_pm(2) = (lambda(2) + sqrt(lambda(2)**2+epsl**2))/2.0
      lambda_pm(3) = (lambda(3) + sqrt(lambda(3)**2+epsl**2))/2.0

      wii=(3.d0-gamma)*(lambda_pm(2)+lambda_pm(3))*alpha**2/(2.d0*(gamma-1.d0))

      fluxp_elem(1) =(lambda_pm(1)*2.d0*(gamma-1.d0) &
                     +lambda_pm(2) &
                     +lambda_pm(3))*rho/(2.d0*gamma)

      fluxp_elem(2) =(lambda_pm(1)*2.d0*(gamma-1.d0)*u &
                     +lambda_pm(2)*u &
                     +lambda_pm(3)*u)*rho/(2.d0*gamma)

      fluxp_elem(3) =(lambda_pm(1)*2.d0*(gamma-1.d0)*v &
                     +lambda_pm(2)*(v-alpha) &
                     +lambda_pm(3)*(v+alpha))*rho/(2.d0*gamma)

      fluxp_elem(4) =(lambda_pm(1)*(gamma-1.d0)*uu2 &
                     +lambda_pm(2)*0.5d0*((v-alpha)**2+u**2) &
                     +lambda_pm(3)*0.5d0*((v+alpha)**2+u**2) &
                     +wii)*rho/(2.d0*gamma)
     

      lambda_pm(1) = (lambda(1)-sqrt(lambda(1)**2+epsl**2))/2.0
      lambda_pm(2) = (lambda(2)-sqrt(lambda(2)**2+epsl**2))/2.0
      lambda_pm(3) = (lambda(3)-sqrt(lambda(3)**2+epsl**2))/2.0

      wii=(3.d0-gamma)*(lambda_pm(2)+lambda_pm(3))*alpha**2/(2.d0*(gamma-1.d0))

      fluxn_elem(1) =(lambda_pm(1)*2.d0*(gamma-1.d0) &
                     +lambda_pm(2)                   &
                     +lambda_pm(3))*rho/(2.d0*gamma)

      fluxn_elem(2) =(lambda_pm(1)*2.d0*(gamma-1.d0)*u &
                     +lambda_pm(2)*u                   &
                     +lambda_pm(3)*u)*rho/(2.d0*gamma)

      fluxn_elem(3) =(lambda_pm(1)*2.d0*(gamma-1.d0)*v &
                     +lambda_pm(2)*(v-alpha)           &
                     +lambda_pm(3)*(v+alpha))*rho/(2.d0*gamma)

      fluxn_elem(4) =(lambda_pm(1)*(gamma-1.d0)*uu2          &
                     +lambda_pm(2)*0.5d0*((v-alpha)**2+u**2) &
                     +lambda_pm(3)*0.5d0*((v+alpha)**2+u**2) &
                     +wii)*rho/(2.d0*gamma)     
    enddo
  enddo  
end subroutine cmpt_flux_sw_eu

subroutine cmpt_flux_vl_eu(flwd)
  implicit none
  type(eu_data_t),target,intent(inout) :: flwd

  type(eu_info_data_t),pointer :: info

  type(rvectorfield_t),pointer :: cons
  type(rvectorfield_t),pointer :: prim

  type(eu_convect_data_t),pointer :: cnvt 
  type(rvectorfield_t),pointer    :: flux
  type(rvectorfield_t),pointer    :: fluxp,fluxn
  
  double precision,dimension(:),pointer :: flux_elem
  double precision,dimension(:),pointer :: cons_elem
  double precision,dimension(:),pointer :: prim_elem
  double precision,dimension(:),pointer :: fluxp_elem,fluxn_elem

  integer :: ibf,jbf
  integer :: i,j
  integer :: i0,j0
  integer :: pre_stage

  double precision :: rho,u,v,e,p,alpha,mach
  double precision :: gamma
  
  info => flwd%info
  if(info%stage_cur.eq.1) then
    pre_stage = info%stage_num
  else
    pre_stage = info%stage_cur - 1
  endif 
  gamma     = info%gamma
  ibf       = info%ibf
  jbf       = info%jbf

  prim => flwd%prim
  cons => flwd%cons(pre_stage)
  cnvt => flwd%cnvt
  
  ! i direction
  flux  => cnvt%flux_xi
  fluxp => cnvt%flux_xip
  fluxn => cnvt%flux_xin
  do j=1,flux%s%d(3)
    do i=1,flux%s%d(2)
      i0 = i
      j0 = j + jbf

      prim_elem  => prim%p(:,i0,j0)
      cons_elem  => cons%p(:,i0,j0)
      fluxp_elem => fluxp%p(:,i,j)
      fluxn_elem => fluxn%p(:,i,j)

      rho   = prim_elem(1)
      u     = prim_elem(2)
      v     = prim_elem(3)
      p     = prim_elem(4)
      alpha = prim_elem(5)
      mach  = u/alpha

      if(abs(mach).le.1.d0) then
        fluxp_elem(1) = rho*alpha*0.25d0*(mach+1)**2
        fluxp_elem(2) = fluxp_elem(1)*((gamma-1.d0)*u+2.d0*alpha)/gamma
        fluxp_elem(3) = fluxp_elem(1)*v
        fluxp_elem(4) = fluxp_elem(1)*(((gamma-1.d0)*u+2.d0*alpha)**2/(2.d0*(gamma**2-1.d0))+0.5d0*v**2)

        fluxn_elem(1) =-rho*alpha*0.25d0*(mach-1)**2
        fluxn_elem(2) = fluxn_elem(1)*((gamma-1.d0)*u-2.d0*alpha)/gamma
        fluxn_elem(3) = fluxn_elem(1)*v
        fluxn_elem(4) = fluxn_elem(1)*(((gamma-1.d0)*u-2.d0*alpha)**2/(2.d0*(gamma**2-1.d0))+0.5d0*v**2)
      else if (mach.gt.1.d0) then
        fluxp_elem(:) = cons_elem(:)*u      
        fluxp_elem(2) = fluxp_elem(2) + p
        fluxp_elem(4) = fluxp_elem(4) + u*p
  
        fluxn_elem(:) = 0.d0  
      else
        fluxp_elem(:) = 0.d0  

        fluxn_elem(:) = cons_elem(:)*u      
        fluxn_elem(2) = fluxn_elem(2) + p
        fluxn_elem(4) = fluxn_elem(4) + u*p
      endif  
    enddo
  enddo  

  ! j direction
  flux  => cnvt%flux_eta
  fluxp => cnvt%flux_etap
  fluxn => cnvt%flux_etan
  do j=1,flux%s%d(3)
    do i=1,flux%s%d(2)
      i0 = j + ibf
      j0 = i 

      prim_elem  => prim%p(:,i0,j0)
      cons_elem  => cons%p(:,i0,j0)
      fluxp_elem => fluxp%p(:,i,j)
      fluxn_elem => fluxn%p(:,i,j)

      rho   = prim_elem(1)
      u     = prim_elem(2)
      v     = prim_elem(3)
      p     = prim_elem(4)
      alpha = prim_elem(5)
      mach  = v/alpha

      if(abs(mach).le.1.d0) then
        fluxp_elem(1) = rho*alpha*0.25d0*(mach+1)**2
        fluxp_elem(2) = fluxp_elem(1)*u
        fluxp_elem(3) = fluxp_elem(1)*((gamma-1.d0)*v+2.d0*alpha)/gamma
        fluxp_elem(4) = fluxp_elem(1)*(((gamma-1.d0)*v+2.d0*alpha)**2/(2.d0*(gamma**2-1.d0))+0.5d0*u**2)

        fluxn_elem(1) =-rho*alpha*0.25d0*(mach-1)**2
        fluxn_elem(2) = fluxn_elem(1)*u
        fluxn_elem(3) = fluxn_elem(1)*((gamma-1.d0)*v-2.d0*alpha)/gamma
        fluxn_elem(4) = fluxn_elem(1)*(((gamma-1.d0)*v-2.d0*alpha)**2/(2.d0*(gamma**2-1.d0))+0.5d0*u**2)
      else if (mach.gt.1.d0) then
        fluxp_elem(:) = cons_elem(:)*v      
        fluxp_elem(3) = fluxp_elem(3) + p
        fluxp_elem(4) = fluxp_elem(4) + v*p
  
        fluxn_elem(:) = 0.d0  
      else
        fluxp_elem(:) = 0.d0  

        fluxn_elem(:) = cons_elem(:)*v      
        fluxn_elem(3) = fluxn_elem(3) + p
        fluxn_elem(4) = fluxn_elem(4) + v*p
      endif
      
    enddo
  enddo  
end subroutine cmpt_flux_vl_eu

subroutine cmpt_flux_ausmfvs_eu(flwd)
  implicit none
  type(eu_data_t),target,intent(inout) :: flwd

  type(eu_info_data_t),pointer :: info

  type(rvectorfield_t),pointer :: cons
  type(rvectorfield_t),pointer :: prim

  type(eu_convect_data_t),pointer :: cnvt 
  type(rvectorfield_t),pointer    :: flux
  type(rvectorfield_t),pointer    :: fluxp,fluxn
  
  double precision,dimension(:),pointer :: flux_elem
  double precision,dimension(:),pointer :: cons_elem
  double precision,dimension(:),pointer :: prim_elem
  double precision,dimension(:),pointer :: fluxp_elem,fluxn_elem

  integer :: ibf,jbf
  integer :: i,j
  integer :: i0,j0
  integer :: pre_stage

  double precision :: rho,u,v,e,p,alpha,mach
  double precision :: machp,machn,p_p,p_n
  
  info => flwd%info
  if(info%stage_cur.eq.1) then
    pre_stage = info%stage_num
  else
    pre_stage = info%stage_cur - 1
  endif 

  ibf       = info%ibf
  jbf       = info%jbf

  prim => flwd%prim
  cons => flwd%cons(pre_stage)
  cnvt => flwd%cnvt
  
  ! i direction
  flux  => cnvt%flux_xi
  fluxp => cnvt%flux_xip
  fluxn => cnvt%flux_xin
  do j=1,flux%s%d(3)
    do i=1,flux%s%d(2)
      i0 = i
      j0 = j + jbf

      prim_elem  => prim%p(:,i0,j0)
      cons_elem  => cons%p(:,i0,j0)
      fluxp_elem => fluxp%p(:,i,j)
      fluxn_elem => fluxn%p(:,i,j)

      rho   = prim_elem(1)
      u     = prim_elem(2)
      p     = prim_elem(4)
      alpha = prim_elem(5)
      mach  = u/alpha
      machp = get_machp(mach)
      machn = get_machn(mach)
      p_p   = get_pp(p,mach)
      p_n   = get_pn(p,mach)
      
      fluxp_elem(:) = cons_elem(:)*alpha*machp      
      fluxp_elem(2) = fluxp_elem(2) + p_p
      fluxp_elem(4) = fluxp_elem(4) + alpha*machp*p

      fluxn_elem(:) = cons_elem(:)*alpha*machn  
      fluxn_elem(2) = fluxn_elem(2) + p_n
      fluxn_elem(4) = fluxn_elem(4) + alpha*machn*p
    enddo
  enddo  

  ! j direction
  flux  => cnvt%flux_eta
  fluxp => cnvt%flux_etap
  fluxn => cnvt%flux_etan
  do j=1,flux%s%d(3)
    do i=1,flux%s%d(2)
      i0 = j + ibf
      j0 = i 

      prim_elem  => prim%p(:,i0,j0)
      cons_elem  => cons%p(:,i0,j0)
      fluxp_elem => fluxp%p(:,i,j)
      fluxn_elem => fluxn%p(:,i,j)

      rho   = prim_elem(1)
      v     = prim_elem(3)
      p     = prim_elem(4)
      alpha = prim_elem(5)
      mach  = v/alpha
      machp = get_machp(mach)
      machn = get_machn(mach)
      p_p   = get_pp(p,mach)
      p_n   = get_pn(p,mach)
      
      fluxp_elem(:) = cons_elem(:)*alpha*machp    
      fluxp_elem(3) = fluxp_elem(3) + p_p
      fluxp_elem(4) = fluxp_elem(4) + alpha*machp*p

      fluxn_elem(:) = cons_elem(:)*alpha*machn  
      fluxn_elem(3) = fluxn_elem(3) + p_n
      fluxn_elem(4) = fluxn_elem(4) + alpha*machn*p
    enddo
  enddo  
end subroutine cmpt_flux_ausmfvs_eu

subroutine cmpt_flux_ausm_eu(flwd)
  type(eu_data_t),intent(inout) :: flwd
end subroutine cmpt_flux_ausm_eu

subroutine cmpt_flux_hllc_eu(flwd, reconst_xi_idx, reconst_eta_idx)
  implicit none
  integer,intent(in) :: reconst_xi_idx, reconst_eta_idx
  type(eu_data_t),intent(inout) :: flwd
  ! type(eu_info_data_t),pointer :: info

  ! type(rvectorfield_t),pointer :: cons
  ! type(rvectorfield_t),pointer :: h
  ! type(rvectorfield_t),pointer :: prim,prim_xi,prim_eta
  ! type(rvectorfield_t),pointer :: priml,primr

  ! type(eu_convect_data_t),pointer :: cnvt 
  ! type(rvectorfield_t),pointer    :: flux
  ! type(rvectorfield_t),pointer    :: fluxp,fluxn
  
  ! double precision,dimension(:),pointer :: flux_elem
  ! double precision,dimension(:),pointer :: cons_elem
  ! double precision,dimension(:),pointer :: prim_elem
  ! double precision,dimension(:),pointer :: prim_elem
  ! double precision,dimension(:),pointer :: priml_elem,primr_elem
  ! double precision,dimension(:),pointer :: fluxp_elem,fluxn_elem

  ! integer :: ibf,jbf
  ! integer :: i,j
  ! integer :: i0,j0
  ! integer :: pre_stage

  ! double precision :: rho,u,v,e,p,alpha,mach
  ! double precision :: rho_l,u_l,v_l,e_l,p_l,alpha_l,mach_l
  ! double precision :: rho_r,u_r,v_r,e_r,p_r,alpha_r,mach_r
  ! double precision :: eta,u_bar,d_bar
  ! double precision :: s_l,s_r,s_star
  ! double precision :: Ustar(4)
  ! double precision :: gamma
  
  ! info => flwd%info
  ! if(istage.eq.1) then
  !   pre_stage = info%stage_num
  ! else
  !   pre_stage = istage - 1
  ! endif 
  ! gamma     = info%gamma
  ! ibf       = info%ibf
  ! jbf       = info%jbf

  ! prim => flwd%prim
  ! cnvt => flwd%cnvt
  
  ! flux  => cnvt%flux_xi
  ! do j=1,flux%s%d(3)
  !   do i=1,flux%s%d(2)
  !     i0 = i
  !     j0 = j + jbf
      
  !     prim_elem  => prim%p(:,i0,j0)
  !     flux_elem  => flux%p(:,i,j)

  !     flux_elem(:) = prim_elem(1:4)
  !   enddo
  ! enddo

  ! flux  => cnvt%flux_eta
  ! do j=1,flux%s%d(3)
  !   do i=1,flux%s%d(2)
  !     i0 = j + ibf
  !     j0 = i 
      
  !     prim_elem  => prim%p(:,i0,j0)
  !     flux_elem  => flux%p(:,i,j)

  !     flux_elem(:) = prim_elem(1:4)
  !   enddo
  ! enddo

  ! call cmpt_rnst_prim_eu(flwd, istage,     &
  !                        reconst_xi_idx,   &
  !                        reconst_eta_idx)
 
  ! h      => cnvt%h_xip
  ! priml  => cnvt%h_xip
  ! primr  => cnvt%h_xin
  ! do j=1,h%s%d(3)
  !   do i=1,h%s%d(2)
  !     flux_elem  => h(:,i,j)
  !     priml_elem => priml(:,i,j)
  !     primr_elem => primr(:,i,j)

  !     rho_l = priml_elem(1)
  !     u_l   = priml_elem(2)
  !     v_l   = priml_elem(3)
  !     p_l   = priml_elem(4)

  !     rho_r = primr_elem(1)
  !     u_r   = primr_elem(2)
  !     v_r   = primr_elem(3)
  !     p_r   = primr_elem(4)
      
  !     alpha_l=sqrt(gamma*p_l/rho_l)
  !     alpha_r=sqrt(gamma*p_r/rho_r)

  !     eta = 0.5*sqrt(rho_l)*sqrt(rho_r)/(sqrt(rho_l)+sqrt(rho_r))**2
  !     u_bar = (sqrt(rho_l)*u_l+sqrt(rho_r)*u_r)/(sqrt(rho_l)+sqrt(rho_r))
  !     d_bar = sqrt((sqrt(rho_l)*c_l**2+sqrt(rho_r)*c_r**2)/(sqrt(rho_l)+sqrt(rho_r))+eta*(u_l-u_r)**2)

  !     S_l = u_bar-d_bar
  !     S_r = u_bar+d_bar
  !     S_star = (p_l-p_r-rho_l*u_l*(S_l-u_l)+rho_r*u_r*(S_r-u_r))/(rho_l*u_l-rho_r*u_r+S_r*rho_r-S_l*rho_l)


  !     if (0.0d0.le.S_l) then
  !       flux_elem(1) = u_l*rho_l
  !       flux_elem(2) = u_l*rho_l*u_l+p_l
  !       flux_elem(2) = u_l*rho_l*v_l
  !       flux_elem(3) = u_l*(0.5d0*rho_l*(u_l**2+v_l**2)+p_l/(gamma-1.d0))+u_l*p_l
  !     elseif((S_l.lt.0.0d0).and.(0.0.le.S_star)) then    
  !       U_star(1) = rho_l*(S_l-u_l)/(S_l-S_star)
  !       U_star(2) = rho_l*(S_l-u_l)/(S_l-S_star)*S_star
  !       U_star(3) = rho_l*(S_l-u_l)/(S_l-S_star)*(0.5*u_l*u_l+p_l/(gamma-1)/rho_l+(S_star-u_l)*(S_star+p_l/rho_l/(S_l-u_l)))         
  !       p_star    = p_l+rho_l*(S_l-u_l)*(S_star-u_l)

  !       flux_elem(1) = S_star*U_star(1)
  !       flux_elem(2) = S_star*U_star(2)+p_star
  !       flux_elem(3) = S_star*U_star(3)+S_star*p_star
  !     elseif((S_star.lt.0.0).and.(0.0.le.S_r)) then  
  !       U_star(1) = rho_r*(S_r-u_r)/(S_r-S_star)
  !       U_star(2) = rho_r*(S_r-u_r)/(S_r-S_star)*S_star
  !       U_star(3) = rho_r*(S_r-u_r)/(S_r-S_star)*(0.5*u_r*u_r+p_r/(gamma-1)/rho_r+(S_star-u_r)*(S_star+p_r/rho_r/(S_r-u_r)))         
  !       p_star    = p_r+rho_r*(S_r-u_r)*(S_star-u_r)

  !       flux_elem(1) = S_star*U_star(1)
  !       flux_elem(2) = S_star*U_star(2)+p_star
  !       flux_elem(3) = S_star*U_star(3)+S_star*p_star
  !     elseif((0.0.gt.S_r)) then                
  !       flux_elem(1)= u_r*rho_r
  !       flux_elem(2)= u_r*rho_r*u_r+p_r
  !       flux_elem(2)= u_r*rho_r*v_r
  !       flux_elem(3)= u_r*(0.5d0*rho_r*(u_r**2+v_r**2)+p_r/(gamma-1.d0))+u_r*p_r
  !     endif    
  !   enddo
  ! enddo

  ! h      => cnvt%h_etap
  ! priml  => cnvt%h_etap
  ! primr  => cnvt%h_etan
  ! do j=1,h%s%d(3)
  !   do i=1,h%s%d(2)
  !     flux_elem  => h(:,i,j)
  !     priml_elem => priml(:,i,j)
  !     primr_elem => primr(:,i,j)
      
  !     rho_l = priml_elem(1)
  !     u_l   = priml_elem(2)
  !     v_l   = priml_elem(3)
  !     p_l   = priml_elem(4)

  !     rho_r = primr_elem(1)
  !     u_r   = primr_elem(2)
  !     v_r   = primr_elem(3)
  !     p_r   = primr_elem(4)

  !     alpha_l=sqrt(gamma*p_l/rho_l)
  !     alpha_r=sqrt(gamma*p_r/rho_r)

  !     eta = 0.5*sqrt(rho_l)*sqrt(rho_r)/(sqrt(rho_l)+sqrt(rho_r))**2
  !     v_bar = (sqrt(rho_l)*v_l+sqrt(rho_r)*v_r)/(sqrt(rho_l)+sqrt(rho_r))
  !     d_bar = sqrt((sqrt(rho_l)*c_l**2+sqrt(rho_r)*c_r**2)/(sqrt(rho_l)+sqrt(rho_r))+eta*(v_l-v_r)**2)

  !     S_l = u_bar-d_bar
  !     S_r = u_bar+d_bar
  !     S_star = (p_l-p_r-rho_l*u_l*(S_l-u_l)+rho_r*u_r*(S_r-u_r))/(rho_l*u_l-rho_r*u_r+S_r*rho_r-S_l*rho_l)

  !     if (0.0d0.le.S_l) then
  !       flux_elem(1)= u_l*rho_l
  !       flux_elem(2)= u_l*rho_l*u_l+p_l
  !       flux_elem(3)= u_l*rho_l*v_l
  !       flux_elem(4)= u_l*(0.5d0*rho_l*(u_l**2+v_l**2)+p_l/(gamma-1.d0))+u_l*p_l
  !     elseif((S_l.lt.0.0d0).and.(0.0.le.S_star)) then    
  !       U_star(1)=rho_l*(S_l-u_l)/(S_l-S_star)
  !       U_star(2)=rho_l*(S_l-u_l)/(S_l-S_star)*S_star
  !       U_star(3)=rho_l*(S_l-u_l)/(S_l-S_star)*(0.5*u_l*u_l+p_l/(gamma-1)/rho_l+(S_star-u_l)*(S_star+p_l/rho_l/(S_l-u_l)))         
  !       p_star = p_l+rho_l*(S_l-u_l)*(S_star-u_l)

  !       h(i,1)= S_star*U_star(1)
  !       h(i,2)= S_star*U_star(2)+p_star
  !       h(i,3)= S_star*U_star(3)+S_star*p_star
  !     elseif((S_star.lt.0.0).and.(0.0.le.S_r)) then  
  !       U_star(1)=rho_r*(S_r-u_r)/(S_r-S_star)
  !       U_star(2)=rho_r*(S_r-u_r)/(S_r-S_star)*S_star
  !       U_star(3)=rho_r*(S_r-u_r)/(S_r-S_star)*(0.5*u_r*u_r+p_r/(gamma-1)/rho_r+(S_star-u_r)*(S_star+p_r/rho_r/(S_r-u_r)))         
  !       p_star = p_r+rho_r*(S_r-u_r)*(S_star-u_r)

  !       h(i,1)= S_star*U_star(1)
  !       h(i,2)= S_star*U_star(2)+p_star
  !       h(i,3)= S_star*U_star(3)+S_star*p_star
  !     elseif((0.0.gt.S_r)) then                
  !       flux_elem(1)= u_r*rho_r
  !       flux_elem(2)= u_r*rho_r*u_r+p_r
  !       flux_elem(2)= u_r*rho_r*v_r
  !       flux_elem(3)= u_r*(0.5d0*rho_r*(u_r**2+v_r**2)+p_r/(gamma-1.d0))+u_r*p_r
  !     endif    
  !   enddo
  ! enddo

end subroutine cmpt_flux_hllc_eu

subroutine cmpt_cnvt_diff(flwd)
  implicit none
  type(eu_data_t),target,intent(inout) :: flwd

  type(eu_info_data_t),pointer    :: info
  type(eu_convect_data_t),pointer :: cnvt

  type(rvectorfield_t),pointer :: h
  type(rvectorfield_t),pointer :: diff
  type(rvectorfield_t),pointer :: diffp,diffn
  type(rvectorfield_t),pointer :: diff_xi,diff_eta
  
  double precision :: invdxi,invdeta
  integer :: isize,jsize,ibf,jbf
  integer :: i,j,i0,j0

  info => flwd%info
  cnvt => flwd%cnvt
  
  isize    = info%isize
  jsize    = info%jsize
  ibf      = info%ibf
  jbf      = info%jbf
  
  invdxi   = 1.d0/info%dxi
  invdeta  = 1.d0/info%deta
  ! Derivatives in xi direction
  h    => cnvt%h_xip
  diff => cnvt%diff_xip
  do j = 1, jsize
    do i = 1, isize
      diff%p(:,i,j) = (h%p(:,i+1,j)-h%p(:,i,j))*invdxi
    enddo
  enddo

  h    => cnvt%h_xin
  diff => cnvt%diff_xin
  do j = 1, jsize
    do i = 1, isize
      diff%p(:,i,j) = (h%p(:,i+1,j)-h%p(:,i,j))*invdxi
    enddo
  enddo

  ! Derivatives in eta direction  
  h    => cnvt%h_etap
  diff => cnvt%diff_etap
  do i = 1, isize
    do j = 1, jsize
      diff%p(:,j,i) = (h%p(:,j+1,i)-h%p(:,j,i))*invdeta 
    enddo
  enddo

  h    => cnvt%h_etan
  diff => cnvt%diff_etan
  do i = 1, isize
    do j = 1, jsize
      diff%p(:,j,i) = (h%p(:,j+1,i)-h%p(:,j,i))*invdeta 
    enddo
  enddo

  ! Assemble derivatives
  diffp => cnvt%diff_xip
  diffn => cnvt%diff_xin
  do j = 1, jsize
    do i = 1, isize
      diffp%p(:,i,j) = diffp%p(:,i,j) + diffn%p(:,i,j) 
    enddo
  enddo
      
  diffp => cnvt%diff_etap
  diffn => cnvt%diff_etan
  do i = 1, isize
    do j = 1, jsize
      diffp%p(:,j,i) = diffp%p(:,j,i) + diffn%p(:,j,i)
    enddo
  enddo

  diff_xi  => cnvt%diff_xip
  diff_eta => cnvt%diff_etap
  do i = 1, isize
    do j = 1, jsize
      i0 = i+ibf
      j0 = j+jbf
      cnvt%cnvt%p(:,i0,j0)=-(diff_xi%p(:,i,j) + diff_eta%p(:,j,i))
    enddo
  enddo
end subroutine 

subroutine cmpt_measure(measure, flux, measure_idx, lo, hi, axis)
  implicit none
  type(rscalarfield_t),target,intent(inout) :: measure
  type(rvectorfield_t),target,intent(inout) :: flux
  integer,intent(in) :: measure_idx
  type(fscalarshape_t),intent(in) :: lo, hi
  integer,intent(in) :: axis


end subroutine

double precision function  get_machp(mach)
implicit none
double precision,intent(in) :: mach
if(abs(mach).le.1.d0) then
  get_machp = 0.25*(mach+1.d0)**2*(1.d0+0.5d0*(mach-1.d0)**2)
else
  get_machp = 0.5d0*(mach+abs(mach))
endif
end function get_machp

double precision function get_machn(mach)
implicit none
double precision,intent(in) :: mach
if(abs(mach).le.1.d0) then
    get_machn =-0.25*(mach-1.d0)**2*(1.d0+0.5d0*(mach+1.d0)**2)
else
    get_machn = 0.5d0*(mach-abs(mach))
endif
end function get_machn

double precision function get_pp(p,mach)
implicit none
double precision,intent(in) :: mach,p
if(abs(mach).le.1.d0) then
    get_pp = 0.25*(mach+1.d0)**2*((2.d0-mach)+3.d0*mach*0.25d0*(mach-1.d0)**2)*p
else
    get_pp = 0.5d0*p*(mach+abs(mach))/mach
endif
end function get_pp

double precision function get_pn(p,mach)
implicit none
double precision,intent(in) :: mach,p
if(abs(mach).le.1.d0) then
    get_pn =-0.25*(mach-1.d0)**2*((-2.d0-mach)+3.d0*mach*0.25d0*(mach+1.d0)**2)*p
else
    get_pn = 0.5d0*p*(mach-abs(mach))/mach
endif
end function

end module eu_solver