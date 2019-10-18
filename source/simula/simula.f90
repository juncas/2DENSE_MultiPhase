module simula
  use mesh
  use model
  use omp_lib
  implicit none

  integer :: istep,nstep,save_step
  double precision :: time,final_time,save_dtime
  integer :: max_thread_num = 8
  logical :: dynamic_threads = .false.

  namelist/simula_params/ nstep,final_time,save_step,save_dtime, &
                          max_thread_num,dynamic_threads
  integer,parameter :: input_id = 10  
  integer,parameter :: input_example_id = 20  
contains
subroutine write_simula_inp()
  implicit none
  logical :: if_exist
  open(input_example_id,file="input_example.in",form='formatted')
  write(input_example_id,NML = simula_params)
  call write_model_inp(input_example_id)
  close(input_example_id)
  stop
end subroutine write_simula_inp

subroutine create_simula()
  implicit none
  logical :: if_exist
  inquire(FILE ='input.in',EXIST = if_exist )
  if(if_exist) then
    open(input_id,file="input.in",form='formatted')
    read(input_id,NML = simula_params)
    call omp_set_dynamic(dynamic_threads)
    call omp_set_num_threads(max_thread_num) 
    call create_mesh(input_id)
    call create_model(input_id)
    close(input_id)
  else
    print*,"Error: create_simula():"
    print*,"File input.in does not exist."
    stop
  endif 
  istep = 0
  time  = 0.d0
	write(*,*) "----- Simulation has been created!"
end subroutine create_simula

subroutine initial_simula()
  implicit none  
  logical :: if_exist
  inquire(FILE ='input.in',EXIST = if_exist )
  if(if_exist) then
    open(input_id,file="input.in",form='formatted')
    call initial_model(input_id)
    close(input_id)
  else
    print*,"Error: initial_simula():"
    print*,"File input.in does not exist."
    stop
  endif 
	write(*,*) "----- Simulation has been initialized!"
  call write_simula(0,0.d0)
end subroutine initial_simula

subroutine evolve_simula()
  implicit none
  double precision :: t0,t1,cputime,elasp_time
  integer :: max_num_threads
	write(*,*) "----- Start evolving simulation..."
  do while((istep.lt.nstep).and. &
            (time.lt.final_time))
    t0 = omp_get_wtime()
    ! call cpu_time(t0) 
    call evolve_model(istep,time,elasp_time)
    ! call cpu_time(t1) 
    t1 = omp_get_wtime()
    cputime = t1 - t0
    call print_status_model(istep,time)
    max_num_threads = omp_get_max_threads()
    write(*,*) "CPU time = ",cputime,max_num_threads
    if(checkpoint_step(istep)) then
      call write_simula(istep,time)
    endif
    if(checkpoint_time(elasp_time)) then
      call write_simula(0,time)
      elasp_time = 0.d0
    endif
  end do
	write(*,*) "----- Finished ..."
  call write_simula(istep,time)
end subroutine evolve_simula  

subroutine destroy_simula()
  implicit none
  call destroy_model()
	write(*,*) "----- Simulation has been destroyed."
end subroutine destroy_simula  

subroutine write_simula(istep,time)
  implicit none
  integer,intent(in) :: istep
  double precision,intent(in) :: time
  call write_model(istep,time)
	write(*,*) "----- Simulation has been written."
end subroutine write_simula  

logical function checkpoint_step(istep)
  implicit none
  integer,intent(in) :: istep  
  checkpoint_step = .false.

  if(mod(istep,save_step).eq.0) then
    checkpoint_step = .true.
  endif
end function checkpoint_step

logical function checkpoint_time(elasp_time)
  implicit none
  double precision,intent(in) :: elasp_time
  
  checkpoint_time = .false.
  if (elasp_time.ge.save_dtime) then
    checkpoint_time = .true.
  endif
end function checkpoint_time
end module simula
