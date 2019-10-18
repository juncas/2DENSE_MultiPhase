module material_c4
use EOS_JWL
implicit none
integer,parameter :: material_type = 1
!-------Reactant
double precision :: a_s     = 778.10e11
double precision :: b_s     =-0.05031e11
double precision :: gamma_s = 0.8938d0
double precision :: cv_s    = 2.4870e6
double precision :: r1_s    = 11.3d0
double precision :: r2_s    = 1.13d0
double precision :: rho0_s  = 1601.d0
NAMELIST/material_c4rtt_params/ a_s, b_s, gamma_s, cv_s, r1_s, &
                                r2_s, rho0_s

!-------Products
double precision :: a_g     = 6.0977e11
double precision :: b_g     = 0.1295e11
double precision :: gamma_g = 0.25d0
double precision :: cv_g    = 1.0e6
double precision :: r1_g    = 4.5d0
double precision :: r2_g    = 1.4d0
double precision :: rho0_g  = 1601.d0
double precision :: q_ref   = 0.09e11
NAMELIST/material_c4pdt_params/ a_g, b_g, gamma_g, cv_g, r1_g, &
                                r2_g, q_ref

!-------Reaction
double precision :: i_r         = 4.0e12
double precision :: lambda_i_r  = 0.022d0
double precision :: a_r         = 0.0367d0
double precision :: b_r         = 0.667d0
double precision :: x_r         = 7.0d0
double precision :: g1_r        = 140.0e6
double precision :: lambda_g1_r = 1.0d0
double precision :: c_r         = 0.667d0
double precision :: d_r         = 0.333d0
double precision :: y_r         = 2.d0
double precision :: g2_r        = 0.d0
double precision :: lambda_g2_r = 0.0d0
double precision :: e_r         = 0.667d0
double precision :: g_r         = 0.667d0
double precision :: z_r         = 3.d0
NAMELIST/material_c4rtn_params/ i_r, a_r, b_r, x_r, c_r,  &
                                d_r, y_r, g1_r, e_r, g_r, &
                                z_r, g2_r, lambda_i_r,    &
                                lambda_g1_r, lambda_g2_r

logical :: if_initialized = .false.

contains
subroutine write_material_c4_inp(funit)
    implicit none
    integer,intent(in) :: funit
    write(funit,NML=material_c4rtt_params)
    write(funit,NML=material_c4pdt_params)
    write(funit,NML=material_c4rtn_params)
end subroutine write_material_c4_inp

subroutine initial_material_c4(funit)
    implicit none
    integer,intent(in) :: funit
    if(.not.if_initialized) then
    read(funit,NML=material_c4rtt_params)
    read(funit,NML=material_c4pdt_params)
    read(funit,NML=material_c4rtn_params)
endif
end subroutine initial_material_c4

subroutine scale_material_c4(rho_rv,u_rv,temp_rv,time_rv, &
                          p_rv,e_rv,l_rv,cv_rv)
    implicit none
    double precision,intent(in) :: rho_rv,u_rv,temp_rv,time_rv, &
                                p_rv,e_rv,l_rv,cv_rv
    
    if(.not.if_initialized) then
        !-------Reactant
    a_s     = a_s/p_rv
    b_s     = b_s/p_rv
    cv_s    = cv_s/cv_rv
    rho0_s  = rho0_s/rho_rv
    !-------Products
    a_g    = a_g/p_rv
    b_g    = b_g/p_rv
    cv_g   = cv_g/cv_rv
    rho0_g = rho0_g/rho_rv
    q_ref  = q_ref/p_rv

    !-------Reaction
    i_r  = i_r*time_rv
    g1_r = g1_r*((1.0e11/p_rv)**(-y_r)/(1.0/time_rv))
    g2_r = g2_r*((1.0e11/p_rv)**(-z_r)/(1.0/time_rv))
    if_initialized = .true.
endif
end subroutine scale_material_c4

subroutine compute_pref_mixt1_c4(p_ref,rho_s,rho_g,lbd)
    implicit none
    double precision,intent(in) :: rho_s,rho_g,lbd 
    double precision,intent(inout) :: p_ref 
    double precision :: pref_s,pref_g,pref_sg,rhogamma_sg
    call compute_pref_eos_jwl(pref_s,rho_s,a_s,b_s,r1_s,r2_s,rho0_s)
    call compute_pref_eos_jwl(pref_g,rho_g,a_g,b_g,r1_g,r2_g,rho0_g)
    
    pref_sg     = (1.d0-lbd)*pref_s/(rho_s*gamma_s) &
                    +lbd*pref_g/(rho_g*gamma_g)
    rhogamma_sg = (1.d0-lbd)/(rho_s*gamma_s) &
                    +lbd/(rho_g*gamma_g)     
    
    p_ref = pref_sg/rhogamma_sg
end subroutine compute_pref_mixt1_c4

subroutine compute_pref_mixt2_c4(p_ref,dp_ref,rho_s,rho_g,lbd,drho_s,drho_g)
    implicit none
    double precision,intent(in) :: rho_s,rho_g,lbd 
    double precision,intent(in) :: drho_s,drho_g
    double precision,intent(inout) :: p_ref,dp_ref
    double precision :: pref_s,pref_g
    double precision :: dpref_s,dpref_g
    double precision :: pref_num,rhogamma_dom
    double precision :: dpref_num,drhogamma_dom
    double precision :: lbd_s,lbd_g
    call compute_dpref2_eos_jwl(pref_s,dpref_s,rho_s,a_s,b_s,r1_s,r2_s,rho0_s)
    call compute_dpref2_eos_jwl(pref_g,dpref_g,rho_g,a_g,b_g,r1_g,r2_g,rho0_g)
    lbd_s = (1.d0-lbd)
    lbd_g = lbd
    pref_num     = lbd_s*pref_s/(rho_s*gamma_s) &
                    +lbd_g*pref_g/(rho_g*gamma_g)
    rhogamma_dom = lbd_s/(rho_s*gamma_s) &
                    +lbd_g/(rho_g*gamma_g)        
    p_ref = pref_num/rhogamma_dom

    dpref_num     = (lbd_s*dpref_s*rho_s*gamma_s - lbd_s*pref_s*gamma_s)/(rho_s*gamma_s)**2*drho_s &
                   +(lbd_g*dpref_g*rho_g*gamma_g - lbd_g*pref_g*gamma_g)/(rho_g*gamma_g)**2*drho_g
    drhogamma_dom =-lbd_s*gamma_s/(rho_s*gamma_s)**2*drho_s &
                   -lbd_g*gamma_g/(rho_g*gamma_g)**2*drho_g       
    dp_ref = (dpref_num*rhogamma_dom - pref_num*drhogamma_dom)/rhogamma_dom**2
end subroutine compute_pref_mixt2_c4

subroutine compute_eref_mixt1_c4(ei_ref,rho_s,rho_g,lbd)
    implicit none
    double precision,intent(in) :: rho_s,rho_g,lbd 
    double precision,intent(inout) :: ei_ref 
    double precision :: eref_s,eref_g
    call compute_eiref_eos_jwl(eref_s,rho_s,a_s,b_s,r1_s,r2_s,rho0_s)
    call compute_eiref_eos_jwl(eref_g,rho_g,a_g,b_g,r1_g,r2_g,rho0_g)
    ei_ref = (1.d0-lbd)*eref_s + lbd*(eref_g - q_ref)
end subroutine compute_eref_mixt1_c4

subroutine compute_eref_mixt2_c4(ei_ref,dei_ref,rho_s,rho_g,lbd,drho_s,drho_g)
    implicit none
    double precision,intent(in) :: rho_s,rho_g,lbd 
    double precision,intent(in) :: drho_s,drho_g
    double precision,intent(inout) :: ei_ref,dei_ref
    double precision :: eref_s,eref_g
    double precision :: deref_s,deref_g
    call compute_deiref2_eos_jwl(eref_s,deref_s,rho_s,a_s,b_s,r1_s,r2_s,rho0_s)
    call compute_deiref2_eos_jwl(eref_g,deref_g,rho_g,a_g,b_g,r1_g,r2_g,rho0_g)
    ei_ref  = (1.d0-lbd)*eref_s + lbd*(eref_g - q_ref)
    dei_ref = (1.d0-lbd)*deref_s*drho_s + lbd*deref_g*drho_g
end subroutine compute_eref_mixt2_c4

subroutine compute_rhogamma_mixt1_c4(rhogamma_mixt,rho_s,rho_g,lbd)
    implicit none
    double precision,intent(in) :: rho_s,rho_g,lbd
    double precision,intent(inout) :: rhogamma_mixt 
    double precision :: rhogamma_sg 
    rhogamma_sg   = (1.d0-lbd)/(rho_s*gamma_s) + lbd/(rho_g*gamma_g)
    rhogamma_mixt = 1.d0/rhogamma_sg
end subroutine compute_rhogamma_mixt1_c4

subroutine compute_rhogamma_mixt2_c4(rhogamma_mixt,drhogamma_mixt,rho_s,rho_g,lbd,drho_s,drho_g)
    implicit none
    double precision,intent(in) :: rho_s,rho_g,lbd
    double precision,intent(in) :: drho_s,drho_g
    double precision,intent(inout) :: rhogamma_mixt,drhogamma_mixt
    double precision :: rhogamma_dom,drhogamma_dom
    double precision :: lbd_s,lbd_g
    lbd_s = (1.d0-lbd)
    lbd_g = lbd
    rhogamma_dom   = lbd_s/(rho_s*gamma_s) + lbd_g/(rho_g*gamma_g)
    rhogamma_mixt  = 1.d0/rhogamma_dom
    drhogamma_dom  =-lbd_s*gamma_s/(rho_s*gamma_s)**2*drho_s &
                    -lbd_g*gamma_g/(rho_g*gamma_g)**2*drho_g
    drhogamma_mixt =-drhogamma_dom/rhogamma_dom**2
end subroutine compute_rhogamma_mixt2_c4

subroutine compute_pressure_mixt1_c4(p,e,rho_s,rho_g,lbd)
    implicit none
    double precision,intent(in) :: e,rho_s,rho_g,lbd
    double precision,intent(inout) :: p 
    double precision :: pref_s,pref_g,eref_s,eref_g
    double precision :: numerator,denominator
    call compute_eiref_eos_jwl(eref_s,rho_s,a_s,b_s,r1_s,r2_s,rho0_s)
    call compute_eiref_eos_jwl(eref_g,rho_g,a_g,b_g,r1_g,r2_g,rho0_g)
    call compute_pref_eos_jwl(pref_s,rho_s,a_s,b_s,r1_s,r2_s,rho0_s)
    call compute_pref_eos_jwl(pref_g,rho_g,a_g,b_g,r1_g,r2_g,rho0_g)
    numerator = e - (1.d0-lbd)*eref_s - lbd*(eref_g - q_ref) &
               +(1.d0-lbd)*pref_s/(rho_s*gamma_s)    &
               +lbd*pref_g/(rho_g*gamma_g)
    denominator = (1.d0-lbd)/(rho_s*gamma_s)+lbd/(rho_g*gamma_g)
    p = numerator/denominator
end subroutine compute_pressure_mixt1_c4

subroutine compute_pressure_mixt2_c4(p,rho_s,rho_g,lbd)
    implicit none
    double precision,intent(in) :: rho_s,rho_g,lbd
    double precision,intent(inout) :: p 
    double precision :: pref_s,pref_g,rgcv_s,rgcv_g
    double precision :: numerator,denominator
    call compute_pref_eos_jwl(pref_s,rho_s,a_s,b_s,r1_s,r2_s,rho0_s)
    call compute_pref_eos_jwl(pref_g,rho_g,a_g,b_g,r1_g,r2_g,rho0_g)
    rgcv_s = 1.d0/(rho_s*gamma_s*cv_s)
    rgcv_g = 1.d0/(rho_g*gamma_g*cv_g)
    numerator   = pref_s*rgcv_s - pref_g*rgcv_g
    denominator = rgcv_s - rgcv_g
    p = numerator/denominator
end subroutine compute_pressure_mixt2_c4

subroutine compute_pressure_mixt3_c4(p,numerator,denominator,rho_s,rho_g,lbd)
    implicit none
    double precision,intent(in) :: rho_s,rho_g,lbd
    double precision,intent(inout) :: p,numerator,denominator 
    double precision :: pref_s,pref_g,rgcv_s,rgcv_g
    call compute_pref_eos_jwl(pref_s,rho_s,a_s,b_s,r1_s,r2_s,rho0_s)
    call compute_pref_eos_jwl(pref_g,rho_g,a_g,b_g,r1_g,r2_g,rho0_g)
    rgcv_s = 1.d0/(rho_s*gamma_s*cv_s)
    rgcv_g = 1.d0/(rho_g*gamma_g*cv_g)
    numerator = pref_s*rgcv_s - pref_g*rgcv_g
    denominator = rgcv_s - rgcv_g
    p = numerator/denominator
end subroutine compute_pressure_mixt3_c4

subroutine compute_pressure_mixt4_c4(num,dom,dnum,ddom,rho_s,rho_g,lbd,drho_s,drho_g)
    implicit none
    double precision,intent(in) :: rho_s,rho_g,lbd
    double precision,intent(in) :: drho_s,drho_g
    double precision,intent(inout) :: num,dom,dnum,ddom
    double precision :: pref_s,pref_g
    double precision :: dpref_s,dpref_g
    double precision :: rgcv_s,rgcv_g
    double precision :: drgcv_s,drgcv_g
    call compute_dpref2_eos_jwl(pref_s,dpref_s,rho_s,a_s,b_s,r1_s,r2_s,rho0_s)
    call compute_dpref2_eos_jwl(pref_g,dpref_g,rho_g,a_g,b_g,r1_g,r2_g,rho0_g)
    rgcv_s = 1.d0/(rho_s*gamma_s*cv_s)
    rgcv_g = 1.d0/(rho_g*gamma_g*cv_g)
    num = pref_s*rgcv_s - pref_g*rgcv_g
    dom = rgcv_s - rgcv_g

    drgcv_s =-gamma_s*cv_s/(rho_s*gamma_s*cv_s)**2
    drgcv_g =-gamma_g*cv_g/(rho_g*gamma_g*cv_g)**2
    dnum = (dpref_s*rgcv_s + pref_s*drgcv_s)*drho_s &
          -(dpref_g*rgcv_g + pref_g*drgcv_g)*drho_g
    ddom = drgcv_s*drho_s - drgcv_g*drho_g
end subroutine compute_pressure_mixt4_c4

! compute the internal energy of reactant
subroutine compute_ei_reac_c4(ei,p,r)
    implicit none
    double precision,intent(in) :: p,r 
    double precision,intent(inout) :: ei 
    double precision :: p_ref,pr_ref,p1_ref
    double precision :: ei_ref,eir_ref,ei1_ref
    call compute_pref_eos_jwl(pr_ref,r,a_s,b_s,r1_s,r2_s,rho0_s)
    call compute_pref_eos_jwl(p1_ref,rho0_s,a_s,b_s,r1_s,r2_s,rho0_s)
    call compute_eiref_eos_jwl(eir_ref,r,a_s,b_s,r1_s,r2_s,rho0_s)
    call compute_eiref_eos_jwl(ei1_ref,rho0_s,a_s,b_s,r1_s,r2_s,rho0_s)
    p_ref  = pr_ref - p1_ref
    ei_ref = eir_ref - ei1_ref
    ei = (p - p_ref)/(r*gamma_s) + ei_ref
end subroutine compute_ei_reac_c4

! compute the internal energy of product
subroutine compute_ei_prod_c4(ei,p,r)
    implicit none
    double precision,intent(in) :: p,r 
    double precision,intent(inout) :: ei 
    double precision :: p_ref, ei_ref
    call compute_pref_eos_jwl(p_ref,r,a_g,b_g,r1_g,r2_g,rho0_g)
    call compute_eiref_eos_jwl(ei_ref,r,a_g,b_g,r1_g,r2_g,rho0_g)
    ei = (p - p_ref)/(r*gamma_g) + ei_ref - q_ref
end subroutine compute_ei_prod_c4

! compute the internal energy of the mixture
subroutine compute_internal_energy_mixt_c4(ei,eia,eib,p,ra,rb,lbd)
    implicit none
    double precision,intent(in) :: p,ra,rb,lbd
    double precision,intent(inout) :: ei,eia,eib 
    double precision :: p_ref, ei_ref
    call compute_ei_reac_c4(eia,p,ra)
    call compute_ei_prod_c4(eib,p,rb)
    ei = (1.d0-lbd)*eia+lbd*eib
end subroutine compute_internal_energy_mixt_c4    

! compute the sound speed of mixture
subroutine compute_sound_speed2_c4(c2,p,rho,rho_s,rho_g,lambda)
    implicit none
    double precision,intent(in) :: p,rho,rho_s,rho_g,lambda
    double precision,intent(inout) :: c2
    double precision :: desdp,desdrs,degdp,degdrg
    double precision :: dersdrs,dergdrg
    double precision :: pref_s,pref_g,dprsdrs,dprgdrg
    double precision :: drsdr,drgdr
    double precision :: dedp,dedrho
    call compute_dedp_eos_jwl(desdp,p,rho_s,gamma_s)
    call compute_dedp_eos_jwl(degdp,p,rho_g,gamma_g)
    dedp  = (1.d0-lambda)*desdp + lambda*degdp

    call compute_dedr_eos_jwl(desdrs,p,rho_s,a_s,b_s,r1_s,r2_s,rho0_s,gamma_s)
    call compute_dedr_eos_jwl(degdrg,p,rho_g,a_g,b_g,r1_g,r2_g,rho0_g,gamma_g)

    drsdr  = rho_s**2/(rho**2*(1.d0-lambda))
    drgdr  = rho_g**2/(rho**2*lambda)
    dedrho = (1.d0-lambda)*desdrs*drsdr + lambda*degdrg*drgdr

    c2 = (p/rho**2-dedrho)/dedp
end subroutine

subroutine compute_xi_c4(xi,r,ra,rb,lbd)
    implicit none
    double precision,intent(in) :: r,ra,rb,lbd
    double precision,intent(inout) :: xi
    double precision :: ggs,ggg,gg

    ggs = gamma_s + 1.d0
    ggg = gamma_g + 1.d0

    gg = ((1.d0 - lbd)*cv_s*ggs + lbd*cv_g*ggg)/((1.d0 - lbd)*cv_s + lbd*cv_g)

    xi = 1.d0/(gg - 1.d0)
end subroutine compute_xi_c4

subroutine compute_temperature_mixt_c4(t,cv,p,rho_s,rho_g,lbd)
    implicit none
    double precision,intent(in) :: p,rho_s,rho_g,lbd
    double precision,intent(inout) :: t,cv
    double precision :: t_s,t_g,pref_s,pref_g
    call compute_pref_eos_jwl(pref_s,rho_s,a_s,b_s,r1_s,r2_s,rho0_s)
    call compute_pref_eos_jwl(pref_g,rho_g,a_g,b_g,r1_g,r2_g,rho0_g)
    t_s = (p - pref_s)/(rho_s*gamma_s*cv_s)
    t_g = (p - pref_g)/(rho_g*gamma_g*cv_g)
    t  = (1.d0-lbd)*t_s+lbd*t_g
    cv = (1.d0-lbd)*cv_s+lbd*cv_g
end subroutine compute_temperature_mixt_c4


! compute the reaction rate of IG model
subroutine compute_source_c4(src0, src_num,                     &
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
    double precision :: rho,p,lbd,z
    double precision :: h0,h1,h2,hv,hp
    p   = comm_prim(4)
    z   = phs_prim(1)
    rho = phs_prim(2)
    lbd = phs_cpnt_prim(3)
    call heaviside_c4(hv,rho/rho0_s-1.d0-a_r)
    call heaviside_c4(hp,p)
    call heaviside_c4(h0,lambda_i_r-lbd)
    call heaviside_c4(h1,lambda_g1_r-lbd)
    call heaviside_c4(h2,lbd-lambda_g2_r)
    src0(1) = i_r*(1.d0-lbd)**b_r*(rho/rho0_s-1.d0-a_r)**x_r*h0*hv &
          +g1_r*(1.d0-lbd)**c_r*lbd**d_r*p**y_r*h1*hp &
          +g2_r*(1.d0-lbd)**e_r*lbd**g_r*p**z_r*h2*hp
    src0(1) = z*src0(1)*rho
end subroutine compute_source_c4

subroutine heaviside_c4(h,x)
    implicit none
    double precision,intent(in) :: x
    double precision,intent(inout) :: h
    if(x.lt.0.d0) then
        h=0.d0
    else
        h=1.d0
    endif
end subroutine heaviside_c4
end module

    