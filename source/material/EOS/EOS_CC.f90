module EOS_CC
  implicit none
contains
subroutine compute_pref_eos_cc(p_ref,rho,a,b,esp1,esp2,rho0)
  implicit none
  double precision,intent(in) :: rho
  double precision,intent(in) :: a,b,esp1,esp2,rho0
  double precision,intent(inout) :: p_ref
  double precision :: v
  v = rho0/rho
  p_ref = a*v**(-esp1) - b*v**(-esp2)
end subroutine compute_pref_eos_cc

subroutine compute_dpref_eos_cc(dp_ref,rho,a,b,esp1,esp2,rho0)
  implicit none
  double precision,intent(in) :: rho
  double precision,intent(in) :: a,b,esp1,esp2,rho0
  double precision,intent(inout) :: dp_ref
  double precision :: v
  v = rho0/rho
  dp_ref = esp1*a*v**(1.d0-esp1)/rho0 &
          -esp2*b*v**(1.d0-esp2)/rho0
end subroutine compute_dpref_eos_cc
  
subroutine compute_dpref2_eos_cc(dp_ref,p_ref,rho,a,b,esp1,esp2,rho0)
  implicit none
  double precision,intent(in) :: rho
  double precision,intent(in) :: a,b,esp1,esp2,rho0
  double precision,intent(inout) :: p_ref,dp_ref
  double precision :: v,vpow1,vpow2
  v = rho0/rho
  vpow1  = v**(-esp1)
  vpow2  = v**(-esp2)
  p_ref  = a*vpow1 &
          -b*vpow2
  dp_ref = esp1*a*vpow1*v/rho0 &
          -esp2*b*vpow2*v/rho0
end subroutine compute_dpref2_eos_cc
  
subroutine compute_eiref_eos_cc(ei_ref,rho,a,b,esp1,esp2,rho0)
  implicit none
  double precision,intent(in) :: rho
  double precision,intent(in) :: a,b,esp1,esp2,rho0
  double precision,intent(inout) :: ei_ref 
  double precision :: v
  v = rho0/rho
  ei_ref =-a*(v**(1.d0-esp1)-1.d0)/(rho0*(1.d0-esp1)) &
          +b*(v**(1.d0-esp2)-1.d0)/(rho0*(1.d0-esp2))
end subroutine compute_eiref_eos_cc

subroutine compute_deiref_eos_cc(dei_ref,rho,a,b,esp1,esp2,rho0)
  implicit none
  double precision,intent(in) :: rho
  double precision,intent(in) :: a,b,esp1,esp2,rho0
  double precision,intent(inout) :: dei_ref 
  double precision :: v
  v = rho0/rho
  dei_ref = a*v**(2.d0-esp1)/rho0**2 &
           -b*v**(2.d0-esp2)/rho0**2
end subroutine compute_deiref_eos_cc

subroutine compute_deiref2_eos_cc(dei_ref,ei_ref,rho,a,b,esp1,esp2,rho0)
  implicit none
  double precision,intent(in) :: rho
  double precision,intent(in) :: a,b,esp1,esp2,rho0
  double precision,intent(inout) :: ei_ref,dei_ref 
  double precision :: v,vpow1,vpow2
  v = rho0/rho
  vpow1 = v**(1.d0-esp1)
  vpow2 = v**(1.d0-esp2)
  ei_ref =-a*(vpow1-1.d0)/(rho0*(1.d0-esp1)) &
          +b*(vpow2-1.d0)/(rho0*(1.d0-esp2))
  dei_ref = a*vpow1*v/rho0**2 &
           -b*vpow2*v/rho0**2
end subroutine compute_deiref2_eos_cc

subroutine compute_gamma_eos_cc(gamma,rho,gamma0)
  implicit none
  double precision,intent(in) :: rho,gamma0
  double precision,intent(inout) :: gamma 
  gamma = gamma0
end subroutine compute_gamma_eos_cc

subroutine compute_xi_eos_cc(xi,rho,gamma0)
  implicit none
  double precision,intent(in) :: rho,gamma0
  double precision,intent(inout) :: xi
  double precision :: gamma
  call compute_gamma_eos_cc(gamma,rho,gamma0)
  xi = 1.d0/gamma
end subroutine compute_xi_eos_cc

! compute the internal energy of product
subroutine compute_ei_eos_cc(ei,p,rho,a,b,esp1,esp2,rho0,gamma0,q0)
  implicit none
  double precision,intent(in) :: p,rho,a,b,esp1,esp2,rho0,gamma0,q0
  double precision,intent(inout) :: ei 
  double precision :: p_ref, ei_ref,gamma
  call compute_pref_eos_cc(p_ref,rho,a,b,esp1,esp2,rho0)
  call compute_eiref_eos_cc(ei_ref,rho,a,b,esp1,esp2,rho0)
  call compute_gamma_eos_cc(gamma,rho,gamma0)
  ei = (p - p_ref)/(rho*gamma) + ei_ref - q0
end subroutine compute_ei_eos_cc

! compute the dedrho
subroutine compute_dedr_eos_cc(dedr,p,rho,a,b,esp1,esp2,rho0,gamma0)
  implicit none
  double precision,intent(in) :: p,rho
  double precision,intent(in) :: a,b,esp1,esp2,rho0,gamma0
  double precision,intent(inout) :: dedr
  double precision :: p_ref, gamma
  double precision :: dprefdr
  double precision :: derefdr
  double precision :: v,vpow1,vpow2
  v = rho0/rho
  vpow1  = v**(-esp1)
  vpow2  = v**(-esp2)
  p_ref  = a*vpow1 &
          -b*vpow2
  dprefdr = esp1*a*vpow1*v/rho0 &
           -esp2*b*vpow2*v/rho0
  derefdr = a*vpow1*v**2/rho0**2 &
           -b*vpow2*v**2/rho0**2
  call compute_gamma_eos_cc(gamma,rho,gamma0)
  dedr = (-dprefdr*gamma*rho-(p-p_ref)*gamma)/(rho*gamma)**2 + derefdr
end subroutine compute_dedr_eos_cc

! compute dedp
subroutine compute_dedp_eos_cc(dedp,rho,gamma0)
  implicit none
  double precision,intent(in) :: rho,gamma0
  double precision,intent(inout) :: dedp
  double precision :: gamma
  call compute_gamma_eos_cc(gamma,rho,gamma0)
  dedp = 1.0/(gamma*rho)
end subroutine compute_dedp_eos_cc

! compute the sound speed
subroutine compute_sound2_eos_cc(c2,p,rho,a,b,esp1,esp2,rho0,gamma0)
  implicit none
  double precision,intent(in) :: p,rho,a,b,esp1,esp2,rho0,gamma0
  double precision,intent(inout) :: c2
  double precision :: dedp,dedr
  call compute_dedr_eos_cc(dedr,p,rho,a,b,esp1,esp2,rho0,gamma0)
  call compute_dedp_eos_cc(dedp,rho,gamma0)
  c2 = (p/rho**2-dedr)/dedp
end subroutine compute_sound2_eos_cc

end module