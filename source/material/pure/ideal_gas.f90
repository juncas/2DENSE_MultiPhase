module material_idealgas
  use EOS_JWL
  implicit none
  integer,parameter :: material_type = 0
  double precision :: cv    = 718.d0*1.225
  double precision :: gamma = 0.4d0
  double precision :: rho0  = 1.225d0
  NAMELIST/material_idealgas_params/ gamma, cv, rho0  

  logical :: if_initialized = .false.
  contains
  subroutine write_material_idealgas_inp(funit)
    implicit none
    integer,intent(in) :: funit
    write(funit,NML=material_idealgas_params)
  end subroutine write_material_idealgas_inp
  
  subroutine initial_material_idealgas(funit)
    implicit none
    integer,intent(in) :: funit
    if(.not.if_initialized) then
    read(funit,NML=material_idealgas_params)
  endif
  end subroutine initial_material_idealgas
  
  subroutine scale_material_idealgas(rho_rv,u_rv,temp_rv,time_rv, &
                            p_rv,e_rv,l_rv,cv_rv)
    implicit none
    double precision,intent(in) :: rho_rv,u_rv,temp_rv,time_rv, &
                                  p_rv,e_rv,l_rv,cv_rv
                        
    if(.not.if_initialized) then          
      cv    = cv/cv_rv
      rho0  = rho0/rho_rv
      if_initialized = .true.
    endif
  end subroutine scale_material_idealgas
  
  ! compute the internal energy of reactant
  subroutine compute_internal_energy_idealgas(ei,p,r)
    implicit none
    double precision,intent(in) :: p,r 
    double precision,intent(inout) :: ei 
    ei = p/(r*gamma)
  end subroutine compute_internal_energy_idealgas
    
  ! compute the sound speed
  subroutine compute_sound_speed2_idealgas(c2,p,rho)
    implicit none
    double precision,intent(in) :: p,rho
    double precision,intent(inout) :: c2
    double precision :: dedp,dedrho
    c2 = (gamma+1.d0)*p/rho
  end subroutine compute_sound_speed2_idealgas
  
  subroutine compute_xi_idealgas(xi,r)
    implicit none
    double precision :: r
    double precision,intent(inout) :: xi
    double precision :: gg
    xi = 1.d0/gamma
  end subroutine compute_xi_idealgas
  
  subroutine compute_temperature_idealgas(t,cv_o,p,rho)
    implicit none
    double precision,intent(in) :: p,rho
    double precision,intent(inout) :: t,cv_o
    t    = p/(rho*gamma*cv)
    cv_o = cv
  end subroutine compute_temperature_idealgas
  
  subroutine compute_eref_idealgas(ei_ref,rho)
    implicit none
    double precision,intent(in) :: rho 
    double precision,intent(inout) :: ei_ref 
    ei_ref = 0.d0
  end subroutine compute_eref_idealgas
  
  subroutine compute_pref_idealgas(p_ref,rho)
    implicit none
    double precision,intent(in) :: rho 
    double precision,intent(inout) :: p_ref 
    p_ref = 0.d0
  end subroutine compute_pref_idealgas
  
  subroutine compute_rhogamma_idealgas(rhogamma,rho)
    implicit none
    double precision,intent(in) :: rho
    double precision,intent(inout) :: rhogamma 
    rhogamma = rho*gamma
  end subroutine compute_rhogamma_idealgas
  
  subroutine compute_pressure_idealgas(p,ei,rho)
    implicit none
    double precision,intent(in) :: rho,ei 
    double precision,intent(inout) :: p 
    p =  gamma*rho*ei
  end subroutine compute_pressure_idealgas  
  
subroutine compute_source_idealgas(src0, src_num,                     &
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
end subroutine compute_source_idealgas
end module