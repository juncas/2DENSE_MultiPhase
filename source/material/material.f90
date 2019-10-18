module material
implicit none
!------Reference Values
double precision :: rho_rv  = 1906.d0
double precision :: u_rv    = 7679.9473d0
double precision :: temp_rv = 300.d0
double precision :: time_rv = 1.0e-6
double precision :: p_rv    = 7679.9473d0**2*1906.d0
double precision :: e_rv    = 7679.9473d0**2
double precision :: l_rv    = 7679.9473d0*1.0e-6
double precision :: cv_rv   = 7679.9473d0**2*1906.d0/298.d0
NAMELIST/material_rfrn_params/ rho_rv, u_rv, temp_rv, time_rv

interface cmpt_pref
  module procedure compute_pref_mixt1
  module procedure compute_pref_mixt2
  module procedure compute_pref_pure
end interface cmpt_pref

interface cmpt_eref
  module procedure compute_eref_mixt1
  module procedure compute_eref_mixt2
  module procedure compute_eref_pure
end interface cmpt_eref

interface cmpt_rg
  module procedure compute_rhogamma_mixt1
  module procedure compute_rhogamma_mixt2
  module procedure compute_rhogamma_pure
end interface cmpt_rg

interface cmpt_pres
  module procedure compute_pressure_mixt1
  module procedure compute_pressure_mixt2
  module procedure compute_pressure_mixt3
  module procedure compute_pressure_mixt4
  module procedure compute_pressure_pure
end interface cmpt_pres

interface cmpt_temp
  module procedure compute_temperature_mixt
  module procedure compute_temperature_pure
end interface cmpt_temp

interface cmpt_ei
  module procedure compute_internal_energy_mixt
  module procedure compute_internal_energy_pure
end interface cmpt_ei

interface cmpt_ss2
  module procedure compute_sound_speed2_mixt
  module procedure compute_sound_speed2_pure
end interface cmpt_ss2

interface cmpt_xi
  module procedure compute_xi_mixt
  module procedure compute_xi_pure
end interface cmpt_xi

integer,parameter :: mat_type_id_pure = 0
integer,parameter :: mat_type_id_mixt = 1

integer,parameter :: mat_id_lx17     = 1
integer,parameter :: mat_id_pbx9502  = 2
integer,parameter :: mat_id_c4       = 3
integer,parameter :: mat_id_copper   = 4
integer,parameter :: mat_id_idealgas = 5
integer,parameter :: mat_id_nm       = 6
integer,parameter :: mat_id_water    = 7
integer,parameter :: mat_id_mxnm     = 8

contains

subroutine write_material_inp(funit)
  use material_lx17,only : write_lx17 => write_material_lx17_inp
  use material_pbx9502,only : write_pbx9502 => write_material_pbx9502_inp
  use material_c4,only : write_c4 => write_material_c4_inp
  use material_copper,only : write_copper => write_material_copper_inp
  use material_idealgas,only : write_idealgas => write_material_idealgas_inp
  use material_nm,only : write_nm => write_material_nm_inp
  use material_water,only : write_water => write_material_water_inp
  use material_mxnm,only : write_mxnm => write_material_mxnm_inp
  implicit none
  integer,intent(in) :: funit
  write(funit,NML=material_rfrn_params)
  call write_lx17(funit)
  call write_pbx9502(funit)
  call write_c4(funit)
  call write_copper(funit)
  call write_idealgas(funit)
  call write_nm(funit)
  call write_water(funit)
  call write_mxnm(funit)
end subroutine write_material_inp

subroutine create_material(funit)
  implicit none
  integer,intent(in) :: funit
  read(funit,NML=material_rfrn_params)

  p_rv    = u_rv**2*rho_rv
  e_rv    = u_rv**2
  l_rv    = u_rv*time_rv
  cv_rv   = u_rv**2*rho_rv/temp_rv

  print*,u_rv,rho_rv,p_rv,e_rv,cv_rv
end subroutine create_material

subroutine initial_material(funit,mat_id)
  use material_lx17,only : scale_lx17 => scale_material_lx17, &
                           initial_lx17 => initial_material_lx17
  use material_pbx9502,only : scale_pbx9502 => scale_material_pbx9502, &
                              initial_pbx9502 => initial_material_pbx9502
  use material_c4,only : scale_c4 => scale_material_c4, &
                         initial_c4 => initial_material_c4  
  use material_copper,only : scale_copper => scale_material_copper, &
                            initial_copper => initial_material_copper  
  use material_idealgas,only : scale_idealgas => scale_material_idealgas, &
                               initial_idealgas => initial_material_idealgas  
  use material_nm,only : scale_nm => scale_material_nm, &
                         initial_nm => initial_material_nm  
  use material_water,only : scale_water => scale_material_water, &
                            initial_water => initial_material_water  
  use material_mxnm,only : scale_mxnm => scale_material_mxnm, &
                           initial_mxnm => initial_material_mxnm  
  implicit none
  integer,intent(in) :: funit,mat_id
  select case(mat_id)
  case(mat_id_lx17)
    call initial_lx17(funit)
    call scale_lx17(rho_rv,u_rv,temp_rv,time_rv, &
                    p_rv,e_rv,l_rv,cv_rv)
  case(mat_id_pbx9502)
    call initial_pbx9502(funit)
    call scale_pbx9502(rho_rv,u_rv,temp_rv,time_rv, &
                    p_rv,e_rv,l_rv,cv_rv)
  case(mat_id_c4)
    call initial_c4(funit)
    call scale_c4(rho_rv,u_rv,temp_rv,time_rv, &
                  p_rv,e_rv,l_rv,cv_rv)
  case(mat_id_copper)
    call initial_copper(funit)
    call scale_copper(rho_rv,u_rv,temp_rv,time_rv, &
                      p_rv,e_rv,l_rv,cv_rv)
  case(mat_id_idealgas)
    call initial_idealgas(funit)
    call scale_idealgas(rho_rv,u_rv,temp_rv,time_rv, &
                        p_rv,e_rv,l_rv,cv_rv)
  case(mat_id_nm)
    call initial_nm(funit)
    call scale_nm(rho_rv,u_rv,temp_rv,time_rv, &
                  p_rv,e_rv,l_rv,cv_rv)
  case(mat_id_water)
    call initial_water(funit)
    call scale_water(rho_rv,u_rv,temp_rv,time_rv, &
                    p_rv,e_rv,l_rv,cv_rv)                    
  case(mat_id_mxnm)
    call initial_mxnm(funit)
    call scale_mxnm(rho_rv,u_rv,temp_rv,time_rv, &
                    p_rv,e_rv,l_rv,cv_rv)
  case default
    print*,"Error: initial_material()"
    print*,"No such mat_id:",mat_id
    stop
  end select
end subroutine initial_material

subroutine destroy_material()
  implicit none
  print*,"Material has been destroyed."
end subroutine destroy_material

subroutine get_mat_type(mat_type_id,mat_id)  
  use material_lx17,only : mat_type_lx17 => material_type  
  use material_pbx9502,only : mat_type_pbx9502 => material_type 
  use material_c4,only : mat_type_c4 => material_type 
  use material_copper,only : mat_type_copper => material_type 
  use material_idealgas,only : mat_type_idealgas => material_type 
  use material_nm,only : mat_type_nm => material_type 
  use material_water,only : mat_type_water => material_type 
  use material_mxnm,only : mat_type_mxnm => material_type 
  implicit none
  integer,intent(in) :: mat_id
  integer,intent(out) :: mat_type_id
  select case(mat_id)
  case(mat_id_lx17)
    mat_type_id = mat_type_lx17
  case(mat_id_pbx9502)
    mat_type_id = mat_type_pbx9502
  case(mat_id_c4)
    mat_type_id = mat_type_c4
  case(mat_id_copper)
    mat_type_id = mat_type_copper
  case(mat_id_idealgas)
    mat_type_id = mat_type_idealgas
  case(mat_id_nm)
    mat_type_id = mat_type_nm
  case(mat_id_water)
    mat_type_id = mat_type_water
  case(mat_id_mxnm)
    mat_type_id = mat_type_mxnm
  case default
    print*,"Error: get_mat_type()"
    print*,"No such mat_id:",mat_id
    stop
  end select
end subroutine get_mat_type

subroutine compute_pref_mixt1(pref,ra,rb,lbd,mixt_id)
  use material_lx17,only : pref_lx17 => compute_pref_mixt1_lx17  
  use material_pbx9502,only : pref_pbx9502 => compute_pref_mixt1_pbx9502  
  use material_c4,only : pref_c4 => compute_pref_mixt1_c4  
  use material_mxnm,only : pref_mxnm => compute_pref_mixt1_mxnm  
  implicit none
  integer,intent(in) :: mixt_id
  double precision,intent(in) :: ra,rb,lbd
  double precision,intent(inout) :: pref
  select case(mixt_id)
  case(mat_id_lx17)
    call pref_lx17(pref,ra,rb,lbd)
  case(mat_id_pbx9502)
    call pref_pbx9502(pref,ra,rb,lbd)
  case(mat_id_c4)
    call pref_c4(pref,ra,rb,lbd)
  case(mat_id_mxnm)
    call pref_mxnm(pref,ra,rb,lbd)
  case default
    print*,"Error: compute_pref_mixt1()"
    print*,"No such mat_id:",mixt_id
    stop
  end select
end subroutine compute_pref_mixt1

subroutine compute_pref_mixt2(pref,dpref,ra,rb,lbd,dra,drb,mixt_id)
  use material_lx17,only : pref_lx17 => compute_pref_mixt2_lx17  
  use material_pbx9502,only : pref_pbx9502 => compute_pref_mixt2_pbx9502  
  use material_c4,only : pref_c4 => compute_pref_mixt2_c4  
  use material_mxnm,only : pref_mxnm => compute_pref_mixt2_mxnm  
  implicit none
  integer,intent(in) :: mixt_id
  double precision,intent(in) :: ra,rb,lbd
  double precision,intent(in) :: dra,drb
  double precision,intent(inout) :: pref,dpref
  select case(mixt_id)
  case(mat_id_lx17)
    call pref_lx17(pref,dpref,ra,rb,lbd,dra,drb)
  case(mat_id_pbx9502)
    call pref_pbx9502(pref,dpref,ra,rb,lbd,dra,drb)
  case(mat_id_c4)
    call pref_c4(pref,dpref,ra,rb,lbd,dra,drb)
  case(mat_id_mxnm)
    call pref_mxnm(pref,dpref,ra,rb,lbd,dra,drb)
  case default
    print*,"Error: compute_pref_mixt2()"
    print*,"No such mat_id:",mixt_id
    stop
  end select
end subroutine compute_pref_mixt2

subroutine compute_pref_pure(pref,r,mat_id)
  use material_copper,only : pref_copper => compute_pref_copper
  use material_idealgas,only : pref_idealgas => compute_pref_idealgas
  use material_nm,only : pref_nm => compute_pref_nm
  use material_water,only : pref_water => compute_pref_water
  implicit none
  integer,intent(in) :: mat_id
  double precision,intent(in) :: r
  double precision,intent(inout) :: pref
  select case(mat_id)
  case(mat_id_copper)
    call pref_copper(pref,r)
  case(mat_id_idealgas)
    call pref_idealgas(pref,r)
  case(mat_id_nm)
    call pref_nm(pref,r)
  case(mat_id_water)
    call pref_water(pref,r)
  case default
    print*,"Error: compute_pref_pure()"
    print*,"No such mat_id:",mat_id
    stop
  end select
end subroutine compute_pref_pure

subroutine compute_eref_mixt1(eref,ra,rb,lbd,mixt_id)
  use material_lx17,only : eref_lx17 => compute_eref_mixt1_lx17
  use material_pbx9502,only : eref_pbx9502 => compute_eref_mixt1_pbx9502
  use material_c4,only : eref_c4 => compute_eref_mixt1_c4
  use material_mxnm,only : eref_mxnm => compute_eref_mixt1_mxnm
  implicit none
  integer,intent(in) :: mixt_id
  double precision,intent(in) :: ra,rb,lbd
  double precision,intent(inout) :: eref
  select case(mixt_id)
  case(mat_id_lx17)
    call eref_lx17(eref,ra,rb,lbd)
  case(mat_id_pbx9502)
    call eref_pbx9502(eref,ra,rb,lbd)
  case(mat_id_c4)
    call eref_c4(eref,ra,rb,lbd)
  case(mat_id_mxnm)
    call eref_mxnm(eref,ra,rb,lbd)
  case default
    print*,"Error: compute_eref_mixt1()"
    print*,"No such mat_id:",mixt_id
    stop
  end select
end subroutine compute_eref_mixt1

subroutine compute_eref_mixt2(eref,derer,ra,rb,lbd,dra,drb,mixt_id)
  use material_lx17,only : eref_lx17 => compute_eref_mixt2_lx17
  use material_pbx9502,only : eref_pbx9502 => compute_eref_mixt2_pbx9502
  use material_c4,only : eref_c4 => compute_eref_mixt2_c4
  use material_mxnm,only : eref_mxnm => compute_eref_mixt2_mxnm
  implicit none
  integer,intent(in) :: mixt_id
  double precision,intent(in) :: ra,rb,lbd
  double precision,intent(in) :: dra,drb
  double precision,intent(inout) :: eref,derer
  select case(mixt_id)
  case(mat_id_lx17)
    call eref_lx17(eref,derer,ra,rb,lbd,dra,drb)
  case(mat_id_pbx9502)
    call eref_pbx9502(eref,derer,ra,rb,lbd,dra,drb)
  case(mat_id_c4)
    call eref_c4(eref,derer,ra,rb,lbd,dra,drb)
  case(mat_id_mxnm)
    call eref_mxnm(eref,derer,ra,rb,lbd,dra,drb)
  case default
    print*,"Error: compute_eref_mixt2()"
    print*,"No such mat_id:",mixt_id
    stop
  end select
end subroutine compute_eref_mixt2

subroutine compute_eref_pure(eref,r,mat_id)
  use material_copper,only : eref_copper => compute_eref_copper
  use material_idealgas,only : eref_idealgas => compute_eref_idealgas
  use material_nm,only : eref_nm => compute_eref_nm
  use material_water,only : eref_water => compute_eref_water
  implicit none
  integer,intent(in) :: mat_id
  double precision,intent(in) :: r
  double precision,intent(inout) :: eref
  select case(mat_id)
  case(mat_id_copper)
    call eref_copper(eref,r)
  case(mat_id_idealgas)
    call eref_idealgas(eref,r)
  case(mat_id_nm)
    call eref_nm(eref,r)
  case(mat_id_water)
    call eref_water(eref,r)
  case default
    print*,"Error: compute_eref_pure()"
    print*,"No such mat_id:",mat_id
    stop
  end select
end subroutine compute_eref_pure  

subroutine compute_rhogamma_mixt1(gamma,ra,rb,lbd,mixt_id)
  use material_lx17,only : rhogamma_lx17 => compute_rhogamma_mixt1_lx17
  use material_pbx9502,only : rhogamma_pbx9502 => compute_rhogamma_mixt1_pbx9502
  use material_c4,only : rhogamma_c4 => compute_rhogamma_mixt1_c4
  use material_mxnm,only : rhogamma_mxnm => compute_rhogamma_mixt1_mxnm
  implicit none
  integer,intent(in) :: mixt_id
  double precision,intent(in) :: ra,rb,lbd
  double precision,intent(inout) :: gamma
  select case(mixt_id)
  case(mat_id_lx17)
    call rhogamma_lx17(gamma,ra,rb,lbd)
  case(mat_id_pbx9502)
    call rhogamma_pbx9502(gamma,ra,rb,lbd)
  case(mat_id_c4)
    call rhogamma_c4(gamma,ra,rb,lbd)
  case(mat_id_mxnm)
    call rhogamma_mxnm(gamma,ra,rb,lbd)
  case default
    print*,"Error: compute_gamma_mixt()"
    print*,"No such mat_id:",mixt_id
    stop
  end select
end subroutine compute_rhogamma_mixt1

subroutine compute_rhogamma_mixt2(gamma,dgamma,ra,rb,lbd,dra,drb,mixt_id)
  use material_lx17,only : rhogamma_lx17 => compute_rhogamma_mixt2_lx17
  use material_pbx9502,only : rhogamma_pbx9502 => compute_rhogamma_mixt2_pbx9502
  use material_c4,only : rhogamma_c4 => compute_rhogamma_mixt2_c4
  use material_mxnm,only : rhogamma_mxnm => compute_rhogamma_mixt2_mxnm
  implicit none
  integer,intent(in) :: mixt_id
  double precision,intent(in) :: ra,rb,lbd
  double precision,intent(in) :: dra,drb
  double precision,intent(inout) :: gamma,dgamma
  select case(mixt_id)
  case(mat_id_lx17)
    call rhogamma_lx17(gamma,dgamma,ra,rb,lbd,dra,drb)
  case(mat_id_pbx9502)
    call rhogamma_pbx9502(gamma,dgamma,ra,rb,lbd,dra,drb)
  case(mat_id_c4)
    call rhogamma_c4(gamma,dgamma,ra,rb,lbd,dra,drb)
  case(mat_id_mxnm)
    call rhogamma_mxnm(gamma,dgamma,ra,rb,lbd,dra,drb)
  case default
    print*,"Error: compute_rhogamma_mixt2()"
    print*,"No such mat_id:",mixt_id
    stop
  end select
end subroutine compute_rhogamma_mixt2

subroutine compute_rhogamma_pure(gamma,r,mat_id)
  use material_copper,only : rhogamma_copper => compute_rhogamma_copper
  use material_idealgas,only : rhogamma_idealgas => compute_rhogamma_idealgas
  use material_nm,only : rhogamma_nm => compute_rhogamma_nm
  use material_water,only : rhogamma_water => compute_rhogamma_water
  implicit none
  integer,intent(in) :: mat_id
  double precision,intent(in) :: r
  double precision,intent(inout) :: gamma
  select case(mat_id)
  case(mat_id_copper)
    call rhogamma_copper(gamma,r)
  case(mat_id_idealgas)
    call rhogamma_idealgas(gamma,r)
  case(mat_id_nm)
    call rhogamma_nm(gamma,r)
  case(mat_id_water)
    call rhogamma_water(gamma,r)
  case default
    print*,"Error: compute_rhogamma_pure()"
    print*,"No such mat_id:",mat_id
    stop
  end select
end subroutine compute_rhogamma_pure

subroutine compute_pressure_mixt1(p_mixt,e,ra,rb,lbd,mixt_id)
  use material_lx17,only : pmixt_lx17 => compute_pressure_mixt1_lx17
  use material_pbx9502,only : pmixt_pbx9502 => compute_pressure_mixt1_pbx9502
  use material_c4,only : pmixt_c4 => compute_pressure_mixt1_c4
  use material_mxnm,only : pmixt_mxnm => compute_pressure_mixt1_mxnm
  implicit none
  integer,intent(in) :: mixt_id
  double precision,intent(in) :: e,ra,rb,lbd
  double precision,intent(inout) :: p_mixt
  select case(mixt_id)
  case(mat_id_lx17)
    call pmixt_lx17(p_mixt,e,ra,rb,lbd)
  case(mat_id_pbx9502)
    call pmixt_pbx9502(p_mixt,e,ra,rb,lbd)
  case(mat_id_c4)
    call pmixt_c4(p_mixt,e,ra,rb,lbd)
  case(mat_id_mxnm)
    call pmixt_mxnm(p_mixt,e,ra,rb,lbd)
  case default
    print*,"Error: compute_pressure_mixt()"
    print*,"No such mat_id:",mixt_id
    stop
  end select
end subroutine compute_pressure_mixt1

subroutine compute_pressure_mixt2(p_mixt,ra,rb,lbd,mixt_id)
  use material_lx17,only : pmixt_lx17 => compute_pressure_mixt2_lx17
  use material_pbx9502,only : pmixt_pbx9502 => compute_pressure_mixt2_pbx9502
  use material_c4,only : pmixt_c4 => compute_pressure_mixt2_c4
  use material_mxnm,only : pmixt_mxnm => compute_pressure_mixt2_mxnm
  implicit none
  integer,intent(in) :: mixt_id
  double precision,intent(in) :: ra,rb,lbd
  double precision,intent(inout) :: p_mixt
  select case(mixt_id)
  case(mat_id_lx17)
    call pmixt_lx17(p_mixt,ra,rb,lbd)
  case(mat_id_pbx9502)
    call pmixt_pbx9502(p_mixt,ra,rb,lbd)
  case(mat_id_c4)
    call pmixt_c4(p_mixt,ra,rb,lbd)
  case(mat_id_mxnm)
    call pmixt_mxnm(p_mixt,ra,rb,lbd)
  case default
    print*,"Error: compute_pressure_mixt()"
    print*,"No such mat_id:",mixt_id
    stop
  end select
end subroutine compute_pressure_mixt2

subroutine compute_pressure_mixt3(p_mixt,numerator,denominator,ra,rb,lbd,mixt_id)
  use material_lx17,only : pmixt_lx17 => compute_pressure_mixt3_lx17
  use material_pbx9502,only : pmixt_pbx9502 => compute_pressure_mixt3_pbx9502
  use material_c4,only : pmixt_c4 => compute_pressure_mixt3_c4
  use material_mxnm,only : pmixt_mxnm => compute_pressure_mixt3_mxnm
  implicit none
  integer,intent(in) :: mixt_id
  double precision,intent(in) :: ra,rb,lbd
  double precision,intent(inout) :: p_mixt,numerator,denominator
  select case(mixt_id)
  case(mat_id_lx17)
    call pmixt_lx17(p_mixt,numerator,denominator,ra,rb,lbd)
  case(mat_id_pbx9502)
    call pmixt_pbx9502(p_mixt,numerator,denominator,ra,rb,lbd)
  case(mat_id_c4)
    call pmixt_c4(p_mixt,numerator,denominator,ra,rb,lbd)
  case(mat_id_mxnm)
    call pmixt_mxnm(p_mixt,numerator,denominator,ra,rb,lbd)
  case default
    print*,"Error: compute_pressure_mixt()"
    print*,"No such mat_id:",mixt_id
    stop
  end select
end subroutine compute_pressure_mixt3

subroutine compute_pressure_mixt4(num,dom,dnum,ddom,ra,rb,lbd,dra,drb,mixt_id)
  use material_lx17,only : pmixt_lx17 => compute_pressure_mixt4_lx17
  use material_pbx9502,only : pmixt_pbx9502 => compute_pressure_mixt4_pbx9502
  use material_c4,only : pmixt_c4 => compute_pressure_mixt4_c4
  use material_mxnm,only : pmixt_mxnm => compute_pressure_mixt4_mxnm
  implicit none
  integer,intent(in) :: mixt_id
  double precision,intent(in) :: ra,rb,lbd,dra,drb
  double precision,intent(inout) :: num,dom,dnum,ddom
  select case(mixt_id)
  case(mat_id_lx17)
    call pmixt_lx17(num,dom,dnum,ddom,ra,rb,lbd,dra,drb)
  case(mat_id_pbx9502)
    call pmixt_pbx9502(num,dom,dnum,ddom,ra,rb,lbd,dra,drb)
  case(mat_id_c4)
    call pmixt_c4(num,dom,dnum,ddom,ra,rb,lbd,dra,drb)
  case(mat_id_mxnm)
    call pmixt_mxnm(num,dom,dnum,ddom,ra,rb,lbd,dra,drb)
  case default
    print*,"Error: compute_pressure_mixt()"
    print*,"No such mat_id:",mixt_id
    stop
  end select
end subroutine compute_pressure_mixt4

subroutine compute_pressure_pure(p_pure,e,r,mat_id)
  use material_copper,only : pressure_copper => compute_pressure_copper
  use material_idealgas,only : pressure_idealgas => compute_pressure_idealgas
  use material_nm,only : pressure_nm => compute_pressure_nm
  use material_water,only : pressure_water => compute_pressure_water
  implicit none
  integer,intent(in) :: mat_id
  double precision,intent(in) :: e,r
  double precision,intent(inout) :: p_pure
  select case(mat_id)
  case(mat_id_copper)
    call pressure_copper(p_pure,e,r)
  case(mat_id_idealgas)
    call pressure_idealgas(p_pure,e,r)
  case(mat_id_nm)
    call pressure_nm(p_pure,e,r)
  case(mat_id_water)
    call pressure_water(p_pure,e,r)
  case default
    print*,"Error: compute_pressure_pure()"
    print*,"No such mat_id:",mat_id
    stop
  end select
end subroutine compute_pressure_pure

subroutine compute_temperature_mixt(t,cv,p,ra,rb,lbd,mixt_id)
  use material_lx17,only : temperature_lx17 => compute_temperature_mixt_lx17
  use material_pbx9502,only : temperature_pbx9502 => compute_temperature_mixt_pbx9502
  use material_c4,only : temperature_c4 => compute_temperature_mixt_c4
  use material_mxnm,only : temperature_mxnm => compute_temperature_mixt_mxnm
  implicit none
  integer,intent(in) :: mixt_id
  double precision,intent(in) :: p,ra,rb,lbd
  double precision,intent(inout) :: t,cv
  select case(mixt_id)
  case(mat_id_lx17)
    call temperature_lx17(t,cv,p,ra,rb,lbd)
  case(mat_id_pbx9502)
    call temperature_pbx9502(t,cv,p,ra,rb,lbd)
  case(mat_id_c4)
    call temperature_c4(t,cv,p,ra,rb,lbd)
  case(mat_id_mxnm)
    call temperature_mxnm(t,cv,p,ra,rb,lbd)
  case default
    print*,"Error: compute_temperature_mixt()"
    print*,"No such mat_id:",mixt_id
    stop
  end select
end subroutine compute_temperature_mixt

subroutine compute_temperature_pure(t,cv,p,r,mat_id)
  use material_copper,only : temperature_copper => compute_temperature_copper
  use material_idealgas,only : temperature_idealgas => compute_temperature_idealgas
  use material_nm,only : temperature_nm => compute_temperature_nm
  use material_water,only : temperature_water => compute_temperature_water
  implicit none
  integer,intent(in) :: mat_id
  double precision,intent(in) :: p,r
  double precision,intent(inout) :: t,cv
  select case(mat_id)
  case(mat_id_copper)
    call temperature_copper(t,cv,p,r)
  case(mat_id_idealgas)
    call temperature_idealgas(t,cv,p,r)
  case(mat_id_nm)
    call temperature_nm(t,cv,p,r)
  case(mat_id_water)
    call temperature_water(t,cv,p,r)
  case default
    print*,"Error: compute_temperature_pure()"
    print*,"No such mat_id:",mat_id
    stop
  end select
end subroutine compute_temperature_pure

subroutine compute_internal_energy_mixt(ei,eia,eib,p,ra,rb,lbd,mixt_id)
  use material_lx17,only : ei_lx17 => compute_internal_energy_mixt_lx17
  use material_pbx9502,only : ei_pbx9502 => compute_internal_energy_mixt_pbx9502
  use material_c4,only : ei_c4 => compute_internal_energy_mixt_c4
  use material_mxnm,only : ei_mxnm => compute_internal_energy_mixt_mxnm
  implicit none
  integer,intent(in) :: mixt_id
  double precision,intent(in) :: p,ra,rb,lbd
  double precision,intent(inout) :: ei,eia,eib
  select case(mixt_id)
  case(mat_id_lx17)
    call ei_lx17(ei,eia,eib,p,ra,rb,lbd)
  case(mat_id_pbx9502)
    call ei_pbx9502(ei,eia,eib,p,ra,rb,lbd)
  case(mat_id_c4)
    call ei_c4(ei,eia,eib,p,ra,rb,lbd)
  case(mat_id_mxnm)
    call ei_mxnm(ei,eia,eib,p,ra,rb,lbd)
  case default
    print*,"Error: compute_internal_energy_mixt()"
    print*,"No such mat_id:",mixt_id
    stop
  end select
end subroutine compute_internal_energy_mixt

subroutine compute_internal_energy_pure(ei,p,r,mat_id)
  use material_copper,only : ei_copper => compute_internal_energy_copper
  use material_idealgas,only : ei_idealgas => compute_internal_energy_idealgas
  use material_nm,only : ei_nm => compute_internal_energy_nm
  use material_water,only : ei_water => compute_internal_energy_water
  implicit none
  integer,intent(in) :: mat_id
  double precision,intent(in) :: p,r
  double precision,intent(inout) :: ei
  select case(mat_id)
  case(mat_id_copper)
    call ei_copper(ei,p,r)
  case(mat_id_idealgas)
    call ei_idealgas(ei,p,r)
  case(mat_id_nm)
    call ei_nm(ei,p,r)
  case(mat_id_water)
    call ei_water(ei,p,r)
  case default
    print*,"Error: compute_internal_energy_pure()"
    print*,"No such mat_id:",mat_id
    stop
  end select
end subroutine compute_internal_energy_pure

subroutine compute_sound_speed2_mixt(c2,p,r,ra,rb,lbd,mixt_id)
  use material_lx17,only : sound_speed2_lx17 => compute_sound_speed2_lx17
  use material_pbx9502,only : sound_speed2_pbx9502 => compute_sound_speed2_pbx9502
  use material_c4,only : sound_speed2_c4 => compute_sound_speed2_c4
  use material_mxnm,only : sound_speed2_mxnm => compute_sound_speed2_mxnm
  implicit none
  integer,intent(in) :: mixt_id
  double precision,intent(in) :: p,r,ra,rb,lbd
  double precision,intent(inout) :: c2
  select case(mixt_id)
  case(mat_id_lx17)
    call sound_speed2_lx17(c2,p,r,ra,rb,lbd)
  case(mat_id_pbx9502)
    call sound_speed2_pbx9502(c2,p,r,ra,rb,lbd)
  case(mat_id_c4)
    call sound_speed2_c4(c2,p,r,ra,rb,lbd)
  case(mat_id_mxnm)
    call sound_speed2_mxnm(c2,p,r,ra,rb,lbd)
  case default
    print*,"Error: compute_sound_speed2_mixt()"
    print*,"No such mat_id:",mixt_id
    stop
  end select
end subroutine compute_sound_speed2_mixt

subroutine compute_sound_speed2_pure(c2,p,r,mat_id)
  use material_copper,only : sound_speed2_copper => compute_sound_speed2_copper
  use material_idealgas,only : sound_speed2_idealgas => compute_sound_speed2_idealgas
  use material_nm,only : sound_speed2_nm => compute_sound_speed2_nm
  use material_water,only : sound_speed2_water => compute_sound_speed2_water
  implicit none
  integer,intent(in) :: mat_id
  double precision,intent(in) :: p,r
  double precision,intent(inout) :: c2
  select case(mat_id)
  case(mat_id_copper)
    call sound_speed2_copper(c2,p,r)
  case(mat_id_idealgas)
    call sound_speed2_idealgas(c2,p,r)
  case(mat_id_nm)
    call sound_speed2_nm(c2,p,r)
  case(mat_id_water)
    call sound_speed2_water(c2,p,r)
  case default
    print*,"Error: compute_sound_speed2_pure()"
    print*,"No such mat_id:",mat_id
    stop
  end select
end subroutine compute_sound_speed2_pure

subroutine compute_xi_mixt(xi,r,ra,rb,lbd,mixt_id)
  use material_lx17,only : xi_lx17 => compute_xi_lx17
  use material_pbx9502,only : xi_pbx9502 => compute_xi_pbx9502
  use material_c4,only : xi_c4 => compute_xi_c4
  use material_mxnm,only : xi_mxnm => compute_xi_mxnm
  implicit none
  integer,intent(in) :: mixt_id
  double precision,intent(in) :: r,ra,rb,lbd
  double precision,intent(inout) :: xi
  select case(mixt_id)
  case(mat_id_lx17)
    call xi_lx17(xi,r,ra,rb,lbd)
  case(mat_id_pbx9502)
    call xi_pbx9502(xi,r,ra,rb,lbd)
  case(mat_id_c4)
    call xi_c4(xi,r,ra,rb,lbd)
  case(mat_id_mxnm)
    call xi_mxnm(xi,r,ra,rb,lbd)
  case default
    print*,"Error: compute_xi_mixt()"
    print*,"No such mat_id:",mixt_id
    stop
  end select
end subroutine compute_xi_mixt

subroutine compute_xi_pure(xi,r,mat_id)
  use material_copper,only : xi_copper => compute_xi_copper
  use material_idealgas,only : xi_idealgas => compute_xi_idealgas
  use material_nm,only : xi_nm => compute_xi_nm
  use material_water,only : xi_water => compute_xi_water
  implicit none
  integer,intent(in) :: mat_id
  double precision,intent(in) :: r
  double precision,intent(inout) :: xi
  select case(mat_id)
  case(mat_id_copper)
    call xi_copper(xi,r)
  case(mat_id_idealgas)
    call xi_idealgas(xi,r)
  case(mat_id_nm)
    call xi_nm(xi,r)
  case(mat_id_water)
    call xi_water(xi,r)
  case default
    print*,"Error: compute_xi_pure()"
    print*,"No such mat_id:",mat_id
    stop
  end select
end subroutine compute_xi_pure

subroutine cmpt_src(src0, src_num,                     &
                    comm_prim, comm_prim_num,          &
                    phs_prim, phs_prim_num,            &
                    phs_cpnt_prim, phs_cpnt_prim_num,  &
                    mat_id)
  use material_lx17,only : source_lx17 => compute_source_lx17
  use material_pbx9502,only : source_pbx9502 => compute_source_pbx9502
  use material_c4,only : source_c4 => compute_source_c4
  use material_mxnm,only : source_mxnm => compute_source_mxnm
  use material_copper,only : source_copper => compute_source_copper
  use material_idealgas,only : source_idealgas => compute_source_idealgas
  use material_nm,only : source_nm => compute_source_nm
  use material_water,only : source_water => compute_source_water
  implicit none
  integer,intent(in) :: mat_id
  integer,intent(in) :: src_num, comm_prim_num, &
                        phs_prim_num, phs_cpnt_prim_num
  double precision,dimension(comm_prim_num),intent(in) :: comm_prim
  double precision,dimension(phs_prim_num),intent(in) :: phs_prim
  double precision,dimension(phs_cpnt_prim_num),intent(in) :: phs_cpnt_prim
  double precision,dimension(src_num),intent(inout) :: src0
  select case(mat_id)
  case(mat_id_lx17)
    call source_lx17(src0, src_num,                     &
                     comm_prim, comm_prim_num,          &
                     phs_prim, phs_prim_num,            &
                     phs_cpnt_prim, phs_cpnt_prim_num)
  case(mat_id_pbx9502)
    call source_pbx9502(src0, src_num,                     &
                        comm_prim, comm_prim_num,          &
                        phs_prim, phs_prim_num,            &
                        phs_cpnt_prim, phs_cpnt_prim_num)
  case(mat_id_c4)
    call source_c4(src0, src_num,                     &
                   comm_prim, comm_prim_num,          &
                   phs_prim, phs_prim_num,            &
                   phs_cpnt_prim, phs_cpnt_prim_num)
  case(mat_id_mxnm)
    call source_mxnm(src0, src_num,                     &
                     comm_prim, comm_prim_num,          &
                     phs_prim, phs_prim_num,            &
                     phs_cpnt_prim, phs_cpnt_prim_num)
  case(mat_id_copper)
    call source_copper(src0,  src_num,                    &
                       comm_prim, comm_prim_num,          &
                       phs_prim, phs_prim_num,            &
                       phs_cpnt_prim, phs_cpnt_prim_num)
  case(mat_id_idealgas)
    call source_idealgas(src0, src_num,                     &
                         comm_prim, comm_prim_num,          &
                         phs_prim, phs_prim_num,            &
                         phs_cpnt_prim, phs_cpnt_prim_num)
  case(mat_id_nm)
    call source_nm(src0, src_num,                     &
                   comm_prim, comm_prim_num,          &
                   phs_prim, phs_prim_num,            &
                   phs_cpnt_prim, phs_cpnt_prim_num)
  case(mat_id_water)
    call source_water(src0, src_num,                     &
                      comm_prim, comm_prim_num,          &
                      phs_prim, phs_prim_num,            &
                      phs_cpnt_prim, phs_cpnt_prim_num)
  case default
    print*,"Error: cmpt_src()"
    print*,"No such mat_id:",mat_id
    stop
  end select
end subroutine cmpt_src

integer function getstatus(lbd,mat_id)
  use material_lx17,only : status_lx17 => getstatus_lx17
  implicit none 
  double precision,intent(in) :: lbd
  integer,intent(in) :: mat_id
  select case(mat_id)
  case(mat_id_lx17)
    getstatus = status_lx17(lbd)
  ! case(mat_id_pbx9502)
  ! case(mat_id_c4)
  ! case(mat_id_mxnm)
  case default
    print*,"Error: cmpt_src()"
    print*,"No such mat_id:",mat_id
    stop
  end select
end function 

end module