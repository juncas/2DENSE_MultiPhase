module material_water
  use EOS_JWL
  implicit none
  integer,parameter :: material_type = 0
  double precision :: a     = 15.85e11
  double precision :: b     =-0.0467e11
  double precision :: r1    = 8.94
  double precision :: r2    = 1.45
  double precision :: cv    = 3.498e6
  double precision :: gamma = 2.0d0
  double precision :: rho0  = 8900.d0
  NAMELIST/material_water_params/ a, b, gamma, cv, r1, &
                                  r2, rho0  

  
  logical :: if_initialized = .false.
contains
subroutine write_material_water_inp(funit)
  implicit none
  integer,intent(in) :: funit
  write(funit,NML=material_water_params)
end subroutine write_material_water_inp

subroutine initial_material_water(funit)
  implicit none
  integer,intent(in) :: funit
  if(.not.if_initialized) then
  read(funit,NML=material_water_params)
endif
end subroutine initial_material_water

subroutine scale_material_water(rho_rv,u_rv,temp_rv,time_rv, &
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
  if_initialized = .true.
endif
end subroutine scale_material_water

! compute the internal energy of reactant
subroutine compute_internal_energy_water(ei,p,r)
  implicit none
  double precision,intent(in) :: p,r 
  double precision,intent(inout) :: ei 
  double precision :: p_ref, ei_ref
  call compute_pref_eos_jwl(p_ref,r,a,b,r1,r2,rho0)
  call compute_eiref_eos_jwl(ei_ref,r,a,b,r1,r2,rho0)
  ei = (p - p_ref)/(r*gamma) + ei_ref
end subroutine compute_internal_energy_water
  
! compute the sound speed
subroutine compute_sound_speed2_water(c2,p,rho)
  implicit none
  double precision,intent(in) :: p,rho
  double precision,intent(inout) :: c2
  double precision :: dedp,dedrho
  call compute_dedp_eos_jwl(dedp,p,rho,gamma)
  call compute_dedr_eos_jwl(dedrho,p,rho,a,b,r1,r2,rho0,gamma)
  c2 = (p/rho**2-dedrho)/dedp
end subroutine compute_sound_speed2_water

subroutine compute_xi_water(xi,r)
  implicit none
  double precision :: r
  double precision,intent(inout) :: xi
  double precision :: gg
  xi = 1.d0/gamma
end subroutine compute_xi_water

subroutine compute_temperature_water(t,cv_o,p,rho)
  implicit none
  double precision,intent(in) :: p,rho
  double precision,intent(inout) :: t,cv_o
  double precision :: pref
  call compute_pref_eos_jwl(pref,rho,a,b,r1,r2,rho0)
  t = (p - pref)/(rho*gamma*cv)
  cv_o = cv
end subroutine compute_temperature_water

subroutine compute_eref_water(ei_ref,rho)
  implicit none
  double precision,intent(in) :: rho 
  double precision,intent(inout) :: ei_ref 
  call compute_eiref_eos_jwl(ei_ref,rho,a,b,r1,r2,rho0)
end subroutine compute_eref_water

subroutine compute_pref_water(p_ref,rho)
  implicit none
  double precision,intent(in) :: rho 
  double precision,intent(inout) :: p_ref 
  call compute_pref_eos_jwl(p_ref,rho,a,b,r1,r2,rho0)
end subroutine compute_pref_water

subroutine compute_rhogamma_water(rhogamma,rho)
  implicit none
  double precision,intent(in) :: rho
  double precision,intent(inout) :: rhogamma 
  rhogamma = rho*gamma
end subroutine compute_rhogamma_water

subroutine compute_pressure_water(p,ei,rho)
  implicit none
  double precision,intent(in) :: rho,ei 
  double precision,intent(inout) :: p 
  double precision :: pref,eref
  call compute_pref_eos_jwl(pref,rho,a,b,r1,r2,rho0)
  p = pref + gamma*rho*(ei-eref)
end subroutine compute_pressure_water  

! compute the reaction rate of IG model
subroutine compute_source_water(src0, src_num,                     &
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
end subroutine compute_source_water
end module