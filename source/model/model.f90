module model
  implicit none

  integer :: model_type = 1
  ! model indices
  integer,parameter :: model_type_idx_eu = 1
  integer,parameter :: model_type_idx_ns = 2
  integer,parameter :: model_type_idx_ig = 3

  NAMELIST/model_params/ model_type

contains
subroutine write_model_inp(funit)
  use eu_model,only: write_model_inp_eu
  ! use ig_model,only: write_model_inp_ig
  implicit none
  integer,intent(in) :: funit
  write(funit,NML=model_params)
  select case(model_type)
    case(model_type_idx_eu)
      call write_model_inp_eu(funit)
    case(model_type_idx_ig)
      ! call write_model_inp_ig(funit)
    case default
      print*,"Error: write_model_inp()"    
      print*,"Invalid model_type index: ",model_type
      stop
  end select
end subroutine write_model_inp

subroutine create_model(funit)
  use eu_model,only: create_model_eu
  ! use ig_model,only: create_model_ig
  implicit none
  integer,intent(in) :: funit
  read(funit,NML=model_params)
  select case(model_type)
    case(model_type_idx_eu)
      call create_model_eu(funit)
    case(model_type_idx_ig)
      ! call create_model_ig(funit)
    case default
      print*,"Error: create_model()"    
      print*,"Invalid model_type index: ",model_type
      stop
  end select
end subroutine create_model

subroutine initial_model(funit)
  use eu_model,only: initial_model_eu
  ! use ig_model,only: initial_model_ig
  implicit none
  integer,intent(in) :: funit
  select case(model_type)
    case(model_type_idx_eu)
      call initial_model_eu(funit)
    case(model_type_idx_ig)
      ! call initial_model_ig(funit)
    case default
      print*,"Error: initial_model()"    
      print*,"Invalid model_type index: ",model_type
      stop
  end select
end subroutine

subroutine evolve_model(istep,time,elasp_time)
  use eu_model,only: evolve_model_eu
  ! use ig_model,only: evolve_model_ig
  implicit none
  integer,intent(inout) :: istep
  double precision,intent(inout) :: time,elasp_time
  select case(model_type)
    case(model_type_idx_eu)
      call evolve_model_eu(istep,time,elasp_time)
    case(model_type_idx_ig)
      ! call evolve_model_ig(istep,time,elasp_time)
    case default
      print*,"Error: evolve_model()"    
      print*,"Invalid model_type index: ",model_type
      stop
  end select
end subroutine 

subroutine destroy_model()
  use eu_model,only: destroy_model_eu
  ! use ig_model,only: destroy_model_ig
  implicit none
  select case(model_type)
    case(model_type_idx_eu)
      call destroy_model_eu()
    case(model_type_idx_ig)
      ! call destroy_model_ig()
    case default
      print*,"Error: destroy_model()"    
      print*,"Invalid model_type index: ",model_type
      stop
  end select
end subroutine 

subroutine print_status_model(istep,time)
  use eu_model,only: print_status_model_eu
  ! use ig_model,only: print_status_model_ig
  implicit none
  integer,intent(in) :: istep
  double precision,intent(in) :: time
  select case(model_type)
    case(model_type_idx_eu)
      call print_status_model_eu(istep,time)
    case(model_type_idx_ig)
      ! call print_status_model_ig(istep,time)
    case default
      print*,"Error: print_status_model()"    
      print*,"Invalid model_type index: ",model_type
      stop
  end select
end subroutine 

subroutine write_model(istep,time)
  use eu_model,only: write_model_eu
  ! use ig_model,only: write_model_ig
  implicit none
  integer,intent(in) :: istep
  double precision,intent(in) :: time
  select case(model_type)
    case(model_type_idx_eu)
      call write_model_eu(istep,time)
    case(model_type_idx_ig)
      ! call write_model_ig(istep,time)
    case default
      print*,"Error: write_model()"    
      print*,"Invalid model_type index: ",model_type
      stop
  end select
end subroutine 

subroutine load_model(istep,time)
  use eu_model,only: load_model_eu
  ! use ig_model,only: load_model_ig
  implicit none
  integer,intent(in) :: istep
  double precision,intent(in) :: time
  select case(model_type)
    case(model_type_idx_eu)
      call load_model_eu(istep,time)
    case(model_type_idx_ig)
      ! call load_model_ig(istep,time)
    case default
      print*,"Error: load_model()"    
      print*,"Invalid model_type index: ",model_type
      stop
  end select
end subroutine 
end module model