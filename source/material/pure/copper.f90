module material_copper
  use EOS_CC
  implicit none
  integer,parameter :: material_type = 0

  double precision :: a     = 1.4567e11
  double precision :: b     = 1.4775e11
  double precision :: eps1  = 2.99
  double precision :: eps2  = 1.99
  double precision :: cv    = 3.498e6
  double precision :: gamma = 2.0d0
  double precision :: rho0  = 8900.d0
  double precision :: eref0 = 0.d0
  NAMELIST/material_copper_params/ a, b, gamma, cv, eps1, &
                                  eps2, rho0, eref0  
  
  logical :: if_initialized = .false.
contains
subroutine write_material_copper_inp(funit)
  implicit none
  integer,intent(in) :: funit
  write(funit,NML=material_copper_params)
end subroutine write_material_copper_inp

subroutine initial_material_copper(funit)
  implicit none
  integer,intent(in) :: funit
  if(.not.if_initialized) then
    read(funit,NML=material_copper_params)
  endif
end subroutine initial_material_copper

subroutine scale_material_copper(rho_rv,u_rv,temp_rv,time_rv, &
                          p_rv,e_rv,l_rv,cv_rv)
  implicit none
  double precision,intent(in) :: rho_rv,u_rv,temp_rv,time_rv, &
                                  p_rv,e_rv,l_rv,cv_rv
  
    if(.not.if_initialized) then
      !-------Reactant
      a     = a/p_rv
      b     = b/p_rv
      cv    = cv/cv_rv
      rho0  = rho0/rho_rv
      
      eref0 = eref0/p_rv
      
      print*,"scaled Copper parameters:"
      print*,"Reactant:",a,b,cv,rho0,eref0
      if_initialized = .true.
    endif
end subroutine scale_material_copper

subroutine compute_eref_copper(ei_ref,rho)
  implicit none
  double precision,intent(in) :: rho 
  double precision,intent(inout) :: ei_ref 
  call compute_eiref_eos_cc(ei_ref,rho,a,b,eps1,eps2,rho0)
  ei_ref = ei_ref - eref0
end subroutine compute_eref_copper

subroutine compute_pref_copper(p_ref,rho)
  implicit none
  double precision,intent(in) :: rho 
  double precision,intent(inout) :: p_ref 
  call compute_pref_eos_cc(p_ref,rho,a,b,eps1,eps2,rho0)
end subroutine compute_pref_copper

! compute the internal energy of reactant
subroutine compute_internal_energy_copper(ei,p,r)
  implicit none
  double precision,intent(in) :: p,r 
  double precision,intent(inout) :: ei 
  double precision :: p_ref, ei_ref
  call compute_pref_eos_cc(p_ref,r,a,b,eps1,eps2,rho0)
  call compute_eiref_eos_cc(ei_ref,r,a,b,eps1,eps2,rho0)
  ei = (p - p_ref)/(r*gamma) + ei_ref - eref0
end subroutine compute_internal_energy_copper
  
! compute the sound speed
subroutine compute_sound_speed2_copper(c2,p,rho)
  implicit none
  double precision,intent(in) :: p,rho
  double precision,intent(inout) :: c2
  double precision :: dedp,dedrho
  call compute_dedp_eos_cc(dedp,rho,gamma)
  call compute_dedr_eos_cc(dedrho,p,rho,a,b,eps1,eps2,rho0,gamma)
  c2 = (p/rho**2-dedrho)/dedp
end subroutine compute_sound_speed2_copper

subroutine compute_xi_copper(xi,r)
  implicit none
  double precision :: r
  double precision,intent(inout) :: xi
  double precision :: gg
  xi = 1.d0/gamma
end subroutine compute_xi_copper

subroutine compute_temperature_copper(t,cv_o,p,rho)
  implicit none
  double precision,intent(in) :: p,rho
  double precision,intent(inout) :: t,cv_o
  double precision :: pref
  call compute_pref_eos_cc(pref,rho,a,b,eps1,eps2,rho0)
  t = (p - pref)/(rho*gamma*cv)
  cv_o = cv
end subroutine compute_temperature_copper


subroutine compute_rhogamma_copper(rhogamma,rho)
  implicit none
  double precision,intent(in) :: rho
  double precision,intent(inout) :: rhogamma 
  rhogamma = rho*gamma
end subroutine compute_rhogamma_copper

subroutine compute_pressure_copper(p,ei,rho)
  implicit none
  double precision,intent(in) :: rho,ei 
  double precision,intent(inout) :: p 
  double precision :: pref,eref
  call compute_pref_eos_cc(pref,rho,a,b,eps1,eps2,rho0)
  call compute_eiref_eos_cc(eref,rho,a,b,eps1,eps2,rho0)
  p = pref + gamma*rho*(ei - eref + eref0)
end subroutine compute_pressure_copper

subroutine compute_source_copper(src0, src_num,                     &
  comm_prim, comm_prim_num,          &
  phs_prim, phs_prim_num,            &
  phs_cpnt_prim, phs_cpnt_prim_num)
  implicit none
  integer,intent(in) :: src_num, comm_prim_num, &
                        phs_prim_num, phs_cpnt_prim_num
  double precision,dimension(comm_prim_num),intent(in) :: comm_prim
  double precision,dimension(phs_prim_num),intent(in) :: phs_prim
  double precision,dimension(phs_cpnt_prim_num),intent(in) :: phs_cpnt_prim
  double precision,dimension(src_num),intent(inout) :: src0
end subroutine compute_source_copper
end module