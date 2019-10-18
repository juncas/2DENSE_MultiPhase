! This is the data of LX-17 condensed phase explosive 
! taken from:
! Kapila A K , Schwendeman D W , Bdzil J B , et al. 
! A study of detonation diffraction in the ignition-and-growth model[J]. 
! Combustion Theory and Modelling, 2007, 11(5):781-822.
module material_lx17
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
double precision :: rho0_s  = 1905.d0
NAMELIST/material_lx17rtt_params/ a_s, b_s, gamma_s, cv_s, r1_s, &
                                r2_s, rho0_s

!-------Products
double precision :: a_g     = 14.81050e11
double precision :: b_g     = 0.63790e11
double precision :: gamma_g = 0.5d0
double precision :: cv_g    = 1.0e6
double precision :: r1_g    = 6.2d0
double precision :: r2_g    = 2.2d0
double precision :: rho0_g  = 1905.d0
double precision :: q_ref   = 0.069e11  
NAMELIST/material_lx17pdt_params/ a_g, b_g, gamma_g, cv_g, r1_g, &
                                  r2_g, rho0_g, q_ref

double precision :: lbd_tol = 1.0e-4, z_tol = 0.2
NAMELIST/material_lx17mixt_params/ lbd_tol, z_tol

!-------Reaction
double precision :: i_r         = 4.0e12
double precision :: b_r         = 0.667d0
double precision :: a_r         = 0.22d0
double precision :: x_r         = 7.0d0
double precision :: c_r         = 0.667d0
double precision :: d_r         = 1.d0
double precision :: y_r         = 3.d0
double precision :: g1_r        = 4500.0e6
double precision :: e_r         = 0.667d0
double precision :: g_r         = 0.667d0
double precision :: z_r         = 1.d0
double precision :: g2_r        = 30.0e6
double precision :: lambda_i_r  = 0.02d0
double precision :: lambda_g1_r = 0.8d0
double precision :: lambda_g2_r = 0.8d0
NAMELIST/material_lx17rtn_params/ i_r, a_r, b_r, x_r, c_r,  &
                                  d_r, y_r, g1_r, e_r, g_r, &
                                  z_r, g2_r, lambda_i_r,    &
                                  lambda_g1_r, lambda_g2_r

logical,private :: if_initialized = .false.

contains
subroutine write_material_lx17_inp(funit)
    implicit none
    integer,intent(in) :: funit
    write(funit,NML=material_lx17rtt_params)
    write(funit,NML=material_lx17pdt_params)
    write(funit,NML=material_lx17mixt_params)
    write(funit,NML=material_lx17rtn_params)
end subroutine write_material_lx17_inp

subroutine initial_material_lx17(funit)
    implicit none
    integer,intent(in) :: funit
    if(.not.if_initialized) then
        read(funit,NML=material_lx17rtt_params)
        read(funit,NML=material_lx17pdt_params)
        read(funit,NML=material_lx17mixt_params)
        read(funit,NML=material_lx17rtn_params)
    endif
end subroutine initial_material_lx17

subroutine scale_material_lx17(rho_rv,u_rv,temp_rv,time_rv, &
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
        g1_r = g1_r*((1.0e11/p_rv)**(-y_r))*time_rv
        g2_r = g2_r*((1.0e11/p_rv)**(-z_r))*time_rv

        print*,"scaled LX-17 parameters:"
        print*,"Reactant:",a_s,b_s,cv_s,rho0_s
        print*,"Product:",a_g,b_g,cv_g,rho0_g,q_ref
        if_initialized = .true.
    endif
end subroutine scale_material_lx17

subroutine compute_pref_mixt1_lx17(p_ref,rho_s,rho_g,lbd)
    implicit none
    double precision,intent(in) :: rho_s,rho_g,lbd 
    double precision,intent(inout) :: p_ref 
    double precision :: pref_s,pref_g,pref_sg,rhogamma_sg
    if(lbd.lt.lbd_tol) then
        call compute_pref_eos_jwl(p_ref,rho_s,a_s,b_s,r1_s,r2_s,rho0_s)
    else if(lbd.gt.(1.d0-lbd_tol)) then    
        call compute_pref_eos_jwl(p_ref,rho_g,a_g,b_g,r1_g,r2_g,rho0_g)
    else
        call compute_pref_eos_jwl(pref_s,rho_s,a_s,b_s,r1_s,r2_s,rho0_s)
        call compute_pref_eos_jwl(pref_g,rho_g,a_g,b_g,r1_g,r2_g,rho0_g)
        
        pref_sg     = (1.d0-lbd)*pref_s/(rho_s*gamma_s) &
                     +       lbd*pref_g/(rho_g*gamma_g)
        rhogamma_sg = (1.d0-lbd)/(rho_s*gamma_s) &
                     +       lbd/(rho_g*gamma_g)     
        
        p_ref       = pref_sg/rhogamma_sg
    endif    
end subroutine compute_pref_mixt1_lx17

subroutine compute_pref_mixt2_lx17(p_ref,dp_ref,rho_s,rho_g,lbd,drho_s,drho_g)
    implicit none
    double precision,intent(in) :: rho_s,rho_g,lbd 
    double precision,intent(in) :: drho_s,drho_g
    double precision,intent(inout) :: p_ref,dp_ref
    double precision :: pref_s,pref_g
    double precision :: dpref_s,dpref_g
    double precision :: pref_num,rhogamma_dom
    double precision :: dpref_num,drhogamma_dom
    double precision :: lbd_s,lbd_g
    if(lbd.lt.lbd_tol) then
        call compute_dpref2_eos_jwl(dp_ref,p_ref,rho_s,a_s,b_s,r1_s,r2_s,rho0_s)
    else if(lbd.gt.(1.d0-lbd_tol)) then 
        call compute_dpref2_eos_jwl(dp_ref,p_ref,rho_g,a_g,b_g,r1_g,r2_g,rho0_g)
    else
        call compute_dpref2_eos_jwl(dpref_s,pref_s,rho_s,a_s,b_s,r1_s,r2_s,rho0_s)
        call compute_dpref2_eos_jwl(dpref_g,pref_g,rho_g,a_g,b_g,r1_g,r2_g,rho0_g)
        lbd_s = (1.d0-lbd)
        lbd_g = lbd
        pref_num     = lbd_s*pref_s/(rho_s*gamma_s) &
                      +lbd_g*pref_g/(rho_g*gamma_g)
        rhogamma_dom = lbd_s/(rho_s*gamma_s) &
                      +lbd_g/(rho_g*gamma_g)        
        p_ref        = pref_num/rhogamma_dom
    
        dpref_num     = (lbd_s*dpref_s*rho_s*gamma_s - lbd_s*pref_s*gamma_s)/(rho_s*gamma_s)**2*drho_s &
                       +(lbd_g*dpref_g*rho_g*gamma_g - lbd_g*pref_g*gamma_g)/(rho_g*gamma_g)**2*drho_g
        drhogamma_dom =-lbd_s*gamma_s/(rho_s*gamma_s)**2*drho_s &
                       -lbd_g*gamma_g/(rho_g*gamma_g)**2*drho_g       
        dp_ref = (dpref_num*rhogamma_dom - pref_num*drhogamma_dom)/rhogamma_dom**2
    endif

end subroutine compute_pref_mixt2_lx17

subroutine compute_eref_mixt1_lx17(ei_ref,rho_s,rho_g,lbd)
    implicit none
    double precision,intent(in) :: rho_s,rho_g,lbd 
    double precision,intent(inout) :: ei_ref 
    double precision :: eref_s,eref_g
    if(lbd.lt.lbd_tol) then
        call compute_eiref_eos_jwl(ei_ref,rho_s,a_s,b_s,r1_s,r2_s,rho0_s)
    else if(lbd.gt.(1.d0-lbd_tol)) then         
        call compute_eiref_eos_jwl(eref_g,rho_g,a_g,b_g,r1_g,r2_g,rho0_g)
        ei_ref = eref_g - q_ref
    else
        call compute_eiref_eos_jwl(eref_s,rho_s,a_s,b_s,r1_s,r2_s,rho0_s)
        call compute_eiref_eos_jwl(eref_g,rho_g,a_g,b_g,r1_g,r2_g,rho0_g)
        ei_ref = (1.d0-lbd)*eref_s + lbd*(eref_g - q_ref)
    endif
end subroutine compute_eref_mixt1_lx17

subroutine compute_eref_mixt2_lx17(ei_ref,dei_ref,rho_s,rho_g,lbd,drho_s,drho_g)
    implicit none
    double precision,intent(in) :: rho_s,rho_g,lbd 
    double precision,intent(in) :: drho_s,drho_g
    double precision,intent(inout) :: ei_ref,dei_ref
    double precision :: eref_s,eref_g
    double precision :: deref_s,deref_g
    if(lbd.lt.lbd_tol) then
        call compute_deiref2_eos_jwl(deref_s,eref_s,rho_s,a_s,b_s,r1_s,r2_s,rho0_s)
        ei_ref  = eref_s
        dei_ref = deref_s
    else if(lbd.gt.(1.d0-lbd_tol)) then 
        call compute_deiref2_eos_jwl(deref_g,eref_g,rho_g,a_g,b_g,r1_g,r2_g,rho0_g)
        ei_ref  = (eref_g - q_ref)
        dei_ref = deref_g        
    else
        call compute_deiref2_eos_jwl(deref_s,eref_s,rho_s,a_s,b_s,r1_s,r2_s,rho0_s)
        call compute_deiref2_eos_jwl(deref_g,eref_g,rho_g,a_g,b_g,r1_g,r2_g,rho0_g)
        ei_ref  = (1.d0-lbd)*eref_s + lbd*(eref_g - q_ref)
        dei_ref = (1.d0-lbd)*deref_s*drho_s + lbd*deref_g*drho_g
    endif
end subroutine compute_eref_mixt2_lx17

subroutine compute_rhogamma_mixt1_lx17(rhogamma_mixt,rho_s,rho_g,lbd)
    implicit none
    double precision,intent(in) :: rho_s,rho_g,lbd
    double precision,intent(inout) :: rhogamma_mixt 
    double precision :: rhogamma_sg 
    if(lbd.lt.lbd_tol) then
        rhogamma_mixt = rho_s*gamma_s
    else if(lbd.gt.(1.d0-lbd_tol)) then 
        rhogamma_mixt = rho_g*gamma_g
    else
        rhogamma_sg   = (1.d0-lbd)/(rho_s*gamma_s) + lbd/(rho_g*gamma_g)
        rhogamma_mixt = 1.d0/rhogamma_sg
    endif
end subroutine compute_rhogamma_mixt1_lx17

subroutine compute_rhogamma_mixt2_lx17(rhogamma_mixt,drhogamma_mixt,rho_s,rho_g,lbd,drho_s,drho_g)
    implicit none
    double precision,intent(in) :: rho_s,rho_g,lbd
    double precision,intent(in) :: drho_s,drho_g
    double precision,intent(inout) :: rhogamma_mixt,drhogamma_mixt
    double precision :: rhogamma_dom,drhogamma_dom
    double precision :: lbd_s,lbd_g
    if(lbd.lt.lbd_tol) then
        rhogamma_mixt  = rho_s*gamma_s
        drhogamma_mixt = gamma_s
    else if(lbd.gt.(1.d0-lbd_tol)) then 
        rhogamma_mixt  = rho_g*gamma_g
        drhogamma_mixt = gamma_g     
    else
        lbd_s = (1.d0-lbd)
        lbd_g = lbd
        rhogamma_dom   = lbd_s/(rho_s*gamma_s) + lbd_g/(rho_g*gamma_g)
        rhogamma_mixt  = 1.d0/rhogamma_dom
        drhogamma_dom  =-lbd_s*gamma_s/(rho_s*gamma_s)**2*drho_s &
                        -lbd_g*gamma_g/(rho_g*gamma_g)**2*drho_g
        drhogamma_mixt =-drhogamma_dom/rhogamma_dom**2
    endif
end subroutine compute_rhogamma_mixt2_lx17

subroutine compute_pressure_mixt1_lx17(p,e,rho_s,rho_g,lbd)
    implicit none
    double precision,intent(in) :: e,rho_s,rho_g,lbd
    double precision,intent(inout) :: p 
    double precision :: pref_s,pref_g,eref_s,eref_g
    double precision :: eia,eib
    double precision :: numerator,denominator
    if(lbd.lt.lbd_tol) then
        call compute_eiref_eos_jwl(eref_s,rho_s,a_s,b_s,r1_s,r2_s,rho0_s)
        call compute_pref_eos_jwl(pref_s,rho_s,a_s,b_s,r1_s,r2_s,rho0_s)
        p = (rho_s*gamma_s)*(e - eref_s) + pref_s
    else if(lbd.gt.(1.d0-lbd_tol)) then 
        call compute_eiref_eos_jwl(eref_g,rho_g,a_g,b_g,r1_g,r2_g,rho0_g)
        call compute_pref_eos_jwl(pref_g,rho_g,a_g,b_g,r1_g,r2_g,rho0_g)
        p = (rho_g*gamma_g)*(e - eref_g + q_ref) + pref_g        
    else
        call compute_eiref_eos_jwl(eref_s,rho_s,a_s,b_s,r1_s,r2_s,rho0_s)
        call compute_pref_eos_jwl(pref_s,rho_s,a_s,b_s,r1_s,r2_s,rho0_s)
        call compute_eiref_eos_jwl(eref_g,rho_g,a_g,b_g,r1_g,r2_g,rho0_g)
        call compute_pref_eos_jwl(pref_g,rho_g,a_g,b_g,r1_g,r2_g,rho0_g)
        numerator = e                                 &
                   -(1.d0-lbd)*eref_s                 &
                   -       lbd*(eref_g - q_ref)       &
                   +(1.d0-lbd)*pref_s/(rho_s*gamma_s) &
                   +       lbd*pref_g/(rho_g*gamma_g)
        denominator = (1.d0-lbd)/(rho_s*gamma_s)+lbd/(rho_g*gamma_g)
        p = numerator/denominator
    endif
end subroutine compute_pressure_mixt1_lx17

subroutine compute_pressure_mixt2_lx17(p,rho_s,rho_g,lbd)
    implicit none
    double precision,intent(in) :: rho_s,rho_g,lbd
    double precision,intent(inout) :: p 
    double precision :: pref_s,pref_g,rgcv_s,rgcv_g
    double precision :: numerator,denominator
    call compute_pref_eos_jwl(pref_s,rho_s,a_s,b_s,r1_s,r2_s,rho0_s)
    call compute_pref_eos_jwl(pref_g,rho_g,a_g,b_g,r1_g,r2_g,rho0_g)
    rgcv_s      = 1.d0/(rho_s*gamma_s*cv_s)
    rgcv_g      = 1.d0/(rho_g*gamma_g*cv_g)
    numerator   = pref_s*rgcv_s - pref_g*rgcv_g
    denominator = rgcv_s - rgcv_g
    p           = numerator/denominator
end subroutine compute_pressure_mixt2_lx17

subroutine compute_pressure_mixt3_lx17(p,numerator,denominator,rho_s,rho_g,lbd)
    implicit none
    double precision,intent(in) :: rho_s,rho_g,lbd
    double precision,intent(inout) :: p,numerator,denominator 
    double precision :: pref_s,pref_g,rgcv_s,rgcv_g
    call compute_pref_eos_jwl(pref_s,rho_s,a_s,b_s,r1_s,r2_s,rho0_s)
    call compute_pref_eos_jwl(pref_g,rho_g,a_g,b_g,r1_g,r2_g,rho0_g)
    rgcv_s      = 1.d0/(rho_s*gamma_s*cv_s)
    rgcv_g      = 1.d0/(rho_g*gamma_g*cv_g)
    numerator   = pref_s*rgcv_s - pref_g*rgcv_g
    denominator = rgcv_s - rgcv_g
    p           = numerator/denominator
end subroutine compute_pressure_mixt3_lx17

subroutine compute_pressure_mixt4_lx17(num,dom,dnum,ddom,rho_s,rho_g,lbd,drho_s,drho_g)
    implicit none
    double precision,intent(in) :: rho_s,rho_g,lbd
    double precision,intent(in) :: drho_s,drho_g
    double precision,intent(inout) :: num,dom,dnum,ddom
    double precision :: pref_s,pref_g
    double precision :: dpref_s,dpref_g
    double precision :: rgcv_s,rgcv_g
    double precision :: drgcv_s,drgcv_g
    call compute_dpref2_eos_jwl(dpref_s,pref_s,rho_s,a_s,b_s,r1_s,r2_s,rho0_s)
    call compute_dpref2_eos_jwl(dpref_g,pref_g,rho_g,a_g,b_g,r1_g,r2_g,rho0_g)
    rgcv_s = 1.d0/(rho_s*gamma_s*cv_s)
    rgcv_g = 1.d0/(rho_g*gamma_g*cv_g)
    num    = pref_s*rgcv_s - pref_g*rgcv_g
    dom    = rgcv_s - rgcv_g

    drgcv_s =-gamma_s*cv_s/(rho_s*gamma_s*cv_s)**2
    drgcv_g =-gamma_g*cv_g/(rho_g*gamma_g*cv_g)**2
    dnum = (dpref_s*rgcv_s + pref_s*drgcv_s)*drho_s &
          -(dpref_g*rgcv_g + pref_g*drgcv_g)*drho_g
    ddom = drgcv_s*drho_s - drgcv_g*drho_g
end subroutine compute_pressure_mixt4_lx17

! compute the internal energy of reactant
subroutine compute_ei_reac_lx17(ei,p,r)
    implicit none
    double precision,intent(in) :: p,r 
    double precision,intent(inout) :: ei 
    double precision :: p_ref
    double precision :: ei_ref
    call compute_pref_eos_jwl(p_ref,r,a_s,b_s,r1_s,r2_s,rho0_s)
    call compute_eiref_eos_jwl(ei_ref,r,a_s,b_s,r1_s,r2_s,rho0_s)
    ei = (p - p_ref)/(r*gamma_s) + ei_ref
end subroutine compute_ei_reac_lx17

! compute the internal energy of product
subroutine compute_ei_prod_lx17(ei,p,r)
    implicit none
    double precision,intent(in) :: p,r 
    double precision,intent(inout) :: ei 
    double precision :: p_ref, ei_ref
    call compute_pref_eos_jwl(p_ref,r,a_g,b_g,r1_g,r2_g,rho0_g)
    call compute_eiref_eos_jwl(ei_ref,r,a_g,b_g,r1_g,r2_g,rho0_g)
    ei = (p - p_ref)/(r*gamma_g) + ei_ref - q_ref
end subroutine compute_ei_prod_lx17

! compute the internal energy of the mixture
subroutine compute_internal_energy_mixt_lx17(ei,eia,eib,p,ra,rb,lbd)
    implicit none
    double precision,intent(in) :: p,ra,rb,lbd
    double precision,intent(inout) :: ei,eia,eib 
    double precision :: p_ref, ei_ref
    if(lbd.lt.lbd_tol) then
        call compute_ei_reac_lx17(eia,p,ra)
        call compute_ei_prod_lx17(eib,p,rb)
        ei = eia
    else if(lbd.gt.(1.d0-lbd_tol)) then 
        call compute_ei_reac_lx17(eia,p,ra)
        call compute_ei_prod_lx17(eib,p,rb)
        ei = eib
    else
        call compute_ei_reac_lx17(eia,p,ra)
        call compute_ei_prod_lx17(eib,p,rb)
        ei = (1.d0-lbd)*eia+lbd*eib
    endif
end subroutine compute_internal_energy_mixt_lx17    

! compute the sound speed of mixture
subroutine compute_sound_speed2_lx17(c2,p,rho,rho_s,rho_g,lbd)
    implicit none
    double precision,intent(in) :: p,rho,rho_s,rho_g,lbd
    double precision,intent(inout) :: c2
    double precision :: desdp,desdrs,degdp,degdrg
    double precision :: dersdrs,dergdrg
    double precision :: pref_s,pref_g,dprsdrs,dprgdrg
    double precision :: drsdr,drgdr,dlbddr
    double precision :: dedp,dedrho
    double precision :: eia,eib
    double precision :: ss2_s,ss2_g
    if(lbd.lt.lbd_tol) then
        call compute_dedp_eos_jwl(desdp,p,rho_s,gamma_s)   
        dedp  = desdp    
        call compute_dedr_eos_jwl(desdrs,p,rho_s,a_s,b_s,r1_s,r2_s,rho0_s,gamma_s)        
        dedrho = desdrs

        c2 = (p/rho**2-dedrho)/dedp
    else if(lbd.gt.(1.d0-lbd_tol)) then    
        call compute_dedp_eos_jwl(degdp,p,rho_g,gamma_g)
        dedp  = degdp
    
        call compute_dedr_eos_jwl(degdrg,p,rho_g,a_g,b_g,r1_g,r2_g,rho0_g,gamma_g)
        
        dedrho = degdrg
        c2 = (p/rho**2-dedrho)/dedp   
    else
        call compute_dedp_eos_jwl(desdp,p,rho_s,gamma_s)    
        call compute_dedp_eos_jwl(degdp,p,rho_g,gamma_g)
        dedp  = (1.d0-lbd)*desdp + lbd*degdp
    
        call compute_dedr_eos_jwl(desdrs,p,rho_s,a_s,b_s,r1_s,r2_s,rho0_s,gamma_s)
        call compute_dedr_eos_jwl(degdrg,p,rho_g,a_g,b_g,r1_g,r2_g,rho0_g,gamma_g)
                
        drsdr  = rho_s**2/(rho**2*(1.d0 - lbd))
        drgdr  = rho_g**2/(rho**2*lbd)
        dedrho = (1.d0-lbd)*desdrs*drsdr + lbd*degdrg*drgdr
        c2 = (p/rho**2-dedrho)/dedp
    endif
    ! if(c2.lt.0.d0) then
    !     print*,"LX-17, negative sound speed:",c2
    !     print*,drsdr,drgdr,dedrho,dedp
    ! endif
end subroutine

subroutine compute_xi_lx17(xi,r,ra,rb,lbd)
    implicit none
    double precision,intent(in) :: r,ra,rb,lbd
    double precision,intent(inout) :: xi
    double precision :: ggs,ggg,gg
    if(lbd.lt.lbd_tol) then
        xi = 1.d0/gamma_s
    else if(lbd.gt.(1.d0-lbd_tol)) then 
        xi = 1.d0/gamma_g     
    else
        ggs = gamma_s + 1.d0
        ggg = gamma_g + 1.d0
        gg = ((1.d0 - lbd)*cv_s*ggs + lbd*cv_g*ggg)/((1.d0 - lbd)*cv_s + lbd*cv_g)
        xi = 1.d0/(gg - 1.d0)
    endif
    if(ISNAN(xi)) then
        print*,gg,lbd,r,ra,rb
    endif
end subroutine compute_xi_lx17

subroutine compute_temperature_mixt_lx17(t,cv,p,rho_s,rho_g,lbd)
    implicit none
    double precision,intent(in) :: p,rho_s,rho_g,lbd
    double precision,intent(inout) :: t,cv
    double precision :: t_s,t_g,pref_s,pref0_s,pref_g
    ! if(lbd.lt.lbd_tol) then
    !     call compute_pref_eos_jwl(pref_s,rho_s,a_s,b_s,r1_s,r2_s,rho0_s)
    !     t_s = (p - pref_s)/(rho_s*gamma_s*cv_s)
    !     t  = t_s
    !     cv = cv_s
    ! else if(lbd.gt.(1.d0-lbd_tol)) then 
    !     call compute_pref_eos_jwl(pref_g,rho_g,a_g,b_g,r1_g,r2_g,rho0_g)
    !     t_g = (p - pref_g)/(rho_g*gamma_g*cv_g)
    !     t  = t_g
    !     cv = cv_g        
    ! else
        call compute_pref_eos_jwl(pref_s,rho_s,a_s,b_s,r1_s,r2_s,rho0_s)
        call compute_pref_eos_jwl(pref_g,rho_g,a_g,b_g,r1_g,r2_g,rho0_g)
        t_s = (p - pref_s)/(rho_s*gamma_s*cv_s)
        t_g = (p - pref_g)/(rho_g*gamma_g*cv_g)
        t  = (1.d0-lbd)*t_s+lbd*t_g
        cv = (1.d0-lbd)*cv_s+lbd*cv_g
    ! endif
end subroutine compute_temperature_mixt_lx17

! compute the reaction rate of IG model
subroutine compute_source_lx17(src0, src_num,                     &
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
    p   = comm_prim(5)
    z   = phs_prim(1)
    rho = phs_prim(2)
    lbd = phs_cpnt_prim(3)
    if(z.gt.z_tol) then
        call heaviside_lx17(hv,rho/rho0_s-1.d0-a_r)
        call heaviside_lx17(hp,p)
        call heaviside_lx17(h0,lambda_i_r-lbd)
        call heaviside_lx17(h1,lambda_g1_r-lbd)
        call heaviside_lx17(h2,lbd-lambda_g2_r)
        src0(1) =  (i_r*(1.d0-lbd)**b_r*(rho/rho0_s-1.d0-a_r)**x_r*h0*hv &
                  +g1_r*(1.d0-lbd)**c_r*lbd**d_r*p**y_r*h1*hp &
                  +g2_r*(1.d0-lbd)**e_r*lbd**g_r*p**z_r*h2*hp)*rho
    else
        src0(1) = 0.d0
    endif
    if(isnan(src0(1))) then
        print*,rho,z,p,lbd
        stop 'rate is a NAN'
    endif
end subroutine compute_source_lx17

subroutine heaviside_lx17(h,x)
    implicit none
    double precision,intent(in) :: x
    double precision,intent(inout) :: h
    if(x.lt.0.d0) then
        h=0.d0
    else
        h=1.d0
    endif
end subroutine heaviside_lx17

integer function getstatus_lx17(lbd)
    implicit none
    double precision,intent(in) :: lbd
    if((lbd.ge.0.d0).and.(lbd.lt.lbd_tol)) then
        getstatus_lx17 = 0
    else if((lbd.ge.lbd_tol).and.(lbd.lt.(1.d0-lbd_tol))) then 
        getstatus_lx17 = 1
    else if((lbd.ge.(1.d0-lbd_tol)).and.(lbd.le.1.d0)) then 
        getstatus_lx17 = 2
    else 
        getstatus_lx17 = -1
    endif
end function

end module

    