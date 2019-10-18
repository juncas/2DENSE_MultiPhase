module EOS_JWL
  implicit none

contains
subroutine compute_pref_eos_jwl(p_ref,rho,a,b,r1,r2,rho0)
  implicit none
  double precision,intent(in) :: rho
  double precision,intent(in) :: a,b,r1,r2,rho0
  double precision,intent(inout) :: p_ref
  double precision :: v
  v = rho0/rho
  p_ref = a*exp(-r1*v) + b*exp(-r2*v)
end subroutine compute_pref_eos_jwl

subroutine compute_dpref_eos_jwl(dpref,rho,a,b,r1,r2,rho0)
  implicit none
  double precision,intent(in) :: rho
  double precision,intent(in) :: a,b,r1,r2,rho0
  double precision,intent(inout) :: dpref
  double precision :: v
  v = rho0/rho
  dpref = a*exp(-r1*v)*r1*rho0/rho**2 &
         +b*exp(-r2*v)*r2*rho0/rho**2
end subroutine compute_dpref_eos_jwl

subroutine compute_dpref2_eos_jwl(dpref,pref,rho,a,b,r1,r2,rho0)
  implicit none
  double precision,intent(in) :: rho
  double precision,intent(in) :: a,b,r1,r2,rho0
  double precision,intent(inout) :: dpref,pref
  double precision :: v,expr1v,expr2v
  v = rho0/rho
  expr1v = exp(-r1*v)
  expr2v = exp(-r2*v)
  pref   = a*expr1v + b*expr2v
  dpref  = a*expr1v*r1*rho0/rho**2 &
          +b*expr2v*r2*rho0/rho**2
end subroutine compute_dpref2_eos_jwl

subroutine compute_eiref_eos_jwl(ei_ref,rho,a,b,r1,r2,rho0)
  implicit none
  double precision,intent(in) :: rho
  double precision,intent(in) :: a,b,r1,r2,rho0
  double precision,intent(inout) :: ei_ref 
  double precision :: v
  v = rho0/rho
  ei_ref = a*exp(-r1*v)/(r1*rho0) &
          +b*exp(-r2*v)/(r2*rho0)
end subroutine compute_eiref_eos_jwl

subroutine compute_deiref_eos_jwl(deref,rho,a,b,r1,r2,rho0)
  implicit none
  double precision,intent(in) :: rho
  double precision,intent(in) :: a,b,r1,r2,rho0
  double precision,intent(inout) :: deref 
  double precision :: v
  v = rho0/rho
  deref = a*exp(-r1*v)/rho**2 &
         +b*exp(-r2*v)/rho**2
end subroutine compute_deiref_eos_jwl

subroutine compute_deiref2_eos_jwl(deref,eref,rho,a,b,r1,r2,rho0)
  implicit none
  double precision,intent(in) :: rho
  double precision,intent(in) :: a,b,r1,r2,rho0
  double precision,intent(inout) :: deref,eref
  double precision :: v,expr1v,expr2v
  v = rho0/rho
  expr1v = exp(-r1*v)
  expr2v = exp(-r2*v)
  eref   = a*expr1v/(r1*rho0) &
          +b*expr2v/(r2*rho0)
  deref  = a*expr1v/rho**2 &
          +b*expr2v/rho**2
end subroutine compute_deiref2_eos_jwl

subroutine compute_gamma_eos_jwl(gamma,rho,gamma0)
  implicit none
  double precision,intent(in) :: rho,gamma0
  double precision,intent(inout) :: gamma 
  gamma = gamma0
end subroutine

subroutine compute_xi_eos_jwl(xi,rho,gamma0)
  implicit none
  double precision,intent(in) :: rho,gamma0
  double precision,intent(inout) :: xi
  double precision :: gamma
  call compute_gamma_eos_jwl(gamma,rho,gamma0)
  xi = 1.d0/gamma
end subroutine

! compute the internal energy of product
subroutine compute_ei_eos_jwl(ei,p,rho,a,b,r1,r2,rho0,gamma0,q0)
  implicit none
  double precision,intent(in) :: p,rho,a,b,r1,r2,rho0,gamma0,q0
  double precision,intent(inout) :: ei 
  double precision :: p_ref, ei_ref,gamma
  call compute_pref_eos_jwl(p_ref,rho,a,b,r1,r2,rho0)
  call compute_eiref_eos_jwl(ei_ref,rho,a,b,r1,r2,rho0)
  call compute_gamma_eos_jwl(gamma,rho,gamma0)
  ei = (p - p_ref)/(rho*gamma) + ei_ref - q0
end subroutine compute_ei_eos_jwl

! compute de/drho
subroutine compute_dedr_eos_jwl(dedr,p,rho,a,b,r1,r2,rho0,gamma0)
  implicit none
  double precision,intent(in) :: p,rho,a,b,r1,r2,rho0,gamma0
  double precision,intent(inout) :: dedr
  double precision :: p_ref, gamma
  double precision :: expr1v,expr2v
  double precision :: dprefdr
  double precision :: derefdr
  double precision :: v
  v = rho0/rho
  call compute_gamma_eos_jwl(gamma,rho,gamma0)
  expr1v = exp(-r1*v)
  expr2v = exp(-r2*v)
  p_ref    = a*expr1v + b*expr2v
  dprefdr  = a*expr1v*r1*rho0/rho**2 &
            +b*expr2v*r2*rho0/rho**2
  derefdr = a*expr1v/rho**2 + b*expr2v/rho**2
  dedr    = (-dprefdr*gamma*rho-(p-p_ref)*gamma)/(rho*gamma)**2 + derefdr
end subroutine compute_dedr_eos_jwl

! compute de/dp
subroutine compute_dedp_eos_jwl(dedp,p,rho,gamma0)
  implicit none
  double precision,intent(in) :: p,rho,gamma0
  double precision,intent(inout) :: dedp
  double precision :: gamma
  call compute_gamma_eos_jwl(gamma,rho,gamma0)
  dedp = 1.0/(gamma*rho)
end subroutine compute_dedp_eos_jwl

! compute the sound speed
subroutine compute_sound_eos_jwl(c2,p,rho,a,b,r1,r2,rho0,gamma0)
  implicit none
  double precision,intent(in) :: p,rho,a,b,r1,r2,rho0,gamma0
  double precision,intent(inout) :: c2
  double precision :: dedp,dedr
  call compute_dedp_eos_jwl(dedp,p,rho,gamma0)
  call compute_dedr_eos_jwl(dedr,p,rho,a,b,r1,r2,rho0,gamma0)
  c2 = (p/rho**2-dedr)/dedp
end subroutine compute_sound_eos_jwl 
end module