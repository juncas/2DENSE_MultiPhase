module material_nm
use EOS_CC
implicit none
integer,parameter :: material_type = 0
double precision :: a     = 0.00819e11
double precision :: b     = 0.0151e11
double precision :: eps1  = 4.53
double precision :: eps2  = 1.42
double precision :: cv    = 1714.d0*1134.d0
double precision :: gamma = 1.19d0
double precision :: rho0  = 1134.d0
double precision :: q0    = 4.48e6*1134.d0
double precision :: cc    = 2.6e9
double precision :: ta    = 11350.d0
NAMELIST/material_nm_params/ a, b, gamma, cv, eps1, &
                                eps2, rho0, q0, cc, ta 

logical :: if_initialized = .false.
contains
subroutine write_material_nm_inp(funit)
  implicit none
  integer,intent(in) :: funit
  write(funit,NML=material_nm_params)
end subroutine write_material_nm_inp

subroutine initial_material_nm(funit)
  implicit none
  integer,intent(in) :: funit
  if(.not.if_initialized) then
  read(funit,NML=material_nm_params)
endif
end subroutine initial_material_nm

subroutine scale_material_nm(rho_rv,u_rv,temp_rv,time_rv, &
                            p_rv,e_rv,l_rv,cv_rv)
  implicit none
  double precision,intent(in) :: rho_rv,u_rv,temp_rv,time_rv, &
                                  p_rv,e_rv,l_rv,cv_rv
  
    if(.not.if_initialized) then
      !-------Reactant
  a    = a/p_rv
  b    = b/p_rv
  cv   = cv/cv_rv
  rho0 = rho0/rho_rv
  q0   = q0/p_rv
  cc   = cc*time_rv
  ta   = ta/temp_rv
  if_initialized = .true.
endif
end subroutine scale_material_nm

! compute the internal energy of reactant
subroutine compute_internal_energy_nm(ei,p,r)
  implicit none
  double precision,intent(in) :: p,r 
  double precision,intent(inout) :: ei 
  double precision :: p_ref, ei_ref
  call compute_pref_eos_cc(p_ref,r,a,b,eps1,eps2,rho0)
  call compute_eiref_eos_cc(ei_ref,r,a,b,eps1,eps2,rho0)
  ei = (p - p_ref)/(r*gamma) + ei_ref
end subroutine compute_internal_energy_nm
    
! compute the sound speed
subroutine compute_sound_speed2_nm(c2,p,rho)
  implicit none
  double precision,intent(in) :: p,rho
  double precision,intent(inout) :: c2
  double precision :: dedp,dedrho
  call compute_dedp_eos_cc(dedp,rho,gamma)
  call compute_dedr_eos_cc(dedrho,p,rho,a,b,eps1,eps2,rho0,gamma)
  c2 = (p/rho**2-dedrho)/dedp
end subroutine compute_sound_speed2_nm

subroutine compute_xi_nm(xi,r)
  implicit none
  double precision :: r
  double precision,intent(inout) :: xi
  double precision :: gg
  xi = 1.d0/gamma
end subroutine compute_xi_nm

subroutine compute_temperature_nm(t,cv_o,p,rho)
  implicit none
  double precision,intent(in) :: p,rho
  double precision,intent(inout) :: t,cv_o
  double precision :: pref
  call compute_pref_eos_cc(pref,rho,a,b,eps1,eps2,rho0)
  t = (p - pref)/(rho*gamma*cv)
  cv_o = cv
end subroutine compute_temperature_nm

subroutine compute_eref_nm(ei_ref,rho)
  implicit none
  double precision,intent(in) :: rho 
  double precision,intent(inout) :: ei_ref 
  call compute_eiref_eos_cc(ei_ref,rho,a,b,eps1,eps2,rho0)
end subroutine compute_eref_nm

subroutine compute_pref_nm(p_ref,rho)
  implicit none
  double precision,intent(in) :: rho 
  double precision,intent(inout) :: p_ref 
  call compute_pref_eos_cc(p_ref,rho,a,b,eps1,eps2,rho0)
end subroutine compute_pref_nm

subroutine compute_rhogamma_nm(rhogamma,rho)
  implicit none
  double precision,intent(in) :: rho
  double precision,intent(inout) :: rhogamma 
  rhogamma = rho*gamma
end subroutine compute_rhogamma_nm

subroutine compute_pressure_nm(p,ei,rho)
  implicit none
  double precision,intent(in) :: rho,ei 
  double precision,intent(inout) :: p 
  double precision :: pref,eref
  call compute_pref_eos_cc(pref,rho,a,b,eps1,eps2,rho0)
  p = pref + gamma*rho*(ei-eref)
end subroutine compute_pressure_nm
    
! compute the reaction rate of IG model
subroutine compute_source_nm(src0, src_num,                     &
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
  double precision :: lbd,rho,t,z
  t   = comm_prim(6)
  z   = phs_prim(1)
  rho = phs_prim(2)
  src0(1) =-lbd*cc*exp(-ta/t)*rho*z
end subroutine compute_source_nm
end module