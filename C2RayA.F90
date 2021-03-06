! Adaptive timestep scheme (Hydrogen+Helium version)
Program C2Ray

  use precision, only: dp
  use clocks, only: setup_clocks, update_clocks, report_clocks
  use astroconstants, only: YEAR
  use my_mpi, only: mpi_setup, mpi_end, rank
  use output_module, only: setup_output, output, close_down
  use grid, only: grid_ini
  use radiation, only: rad_ini, sourcetype
  use cosmology, only: cosmology_init, redshift_evol, cosmo_evol
  use c2ray_parameters, only: cosmological
  use material, only: mat_ini, zred00, isothermal
  use radiative_cooling, only: setup_cool
  use sourceprops, only: source_properties_ini
  use evolve, only: evolve1D
  use file_admin, only: stdinput, logf, flag_for_file_input
  use times, only: time_ini, end_time, output_time
  use timestep, only: time_step
  use sizes, only: mesh
  use LTE, only: LTE_calculation

#ifdef XLF
  ! Modules for the xlf (IBM) compiler
  USE XLFUTILITY, only: iargc, getarg, flush => flush_
#endif

  implicit none

  ! CPU time variables
  real(kind=dp) :: front_sim_time                       !< actual time seen by I-front
  real(kind=dp) :: next_output_time               	!< time of next output
  real(kind=dp) :: front_actual_dt 	                !< timestep taken by I-front
  character(len=512) :: inputfile	                !< Input file
  character(len=3) :: powerlaw_string !< temporary
  real(kind=dp) :: powerlaw !< temporary
  
  ! Adative timestep variables
  integer, parameter :: time_partition = 5              !< I-front takes this number of timestep to proceed
  integer :: i_partition                                !< dummy variable of time partition
  real(kind=dp), dimension(1:time_partition) :: xfinal  !< targeting ionized fractions to be achieved 
  real(kind=dp) :: front_dt                             !< time step taken by the I-front
  integer :: front_pos                                  !< position of I-front
  integer :: front_index                                !< count the number of steps taken by I-front in one cell
  integer :: inside_pos                                 !< position of the verge of HII zone
  integer :: equilibrium_pos                            !< position of the verge of equilibrium zone (until next output)
  logical :: output_time_is_reached                     !< is true when an output is made

  real(kind=dp) :: background_HII

  ! Initialize clocks (cpu and wall)
  call setup_clocks

  ! Set up MPI structure
  call mpi_setup()			

  ! Set up input stream (either standard input or from file given
  ! by first argument)
  if (rank == 0) then
     write(logf,*) "screen input or file input?"
     flush(logf)
     if (COMMAND_ARGUMENT_COUNT () > 0) then
        call GET_COMMAND_ARGUMENT(1,inputfile)
        write(logf,*) "reading input from ",trim(adjustl(inputfile))
        open(unit=stdinput,file=inputfile)
        call flag_for_file_input(.true.)
        call GET_COMMAND_ARGUMENT(2,powerlaw_string) !< temporary
        read(powerlaw_string,*) powerlaw !< temporary
     else
        write(logf,*) "reading input from command line"
     endif
     flush(logf)
  endif

  ! Initialize output
  call setup_output()     	

  ! Initialize grid	 
  call grid_ini()         	

  ! Initialize the material properties
  call mat_ini (background_HII)  	                

  ! Setup cooling
  if (.not.isothermal) call setup_cool()

  ! Initialize photo-ionization calculation 
  call rad_ini(powerlaw)         		 

  ! Initialize source property
  call source_properties_ini (sourcetype)

  ! Initialize time step parameters
  call time_ini()        		

  ! Initialize actual time to zero
  front_sim_time=0.0         

  ! An output will be made by that time       	
  next_output_time=0.0          

  ! Initialize cosmology         
  if (cosmological) then  		    
     call cosmology_init(zred00,front_sim_time)
     call redshift_evol(front_sim_time)
     call cosmo_evol( )
  endif

  ! Update clock counters (cpu + wall, to avoid overflowing the counter)
  call update_clocks ()

  ! Report clocks (cpu and wall)
  call report_clocks ()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                    !!
!!                Loop until end time is reached                      !!
!!                                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  inside_pos=0
  front_pos=1
  front_index=1
  equilibrium_pos=0
  output_time_is_reached=.false.

  ! Assignment of expected ionization fraction at the I-front
  !do i_partition=1,time_partition
  !  xfinal(i_partition) = real(i_partition)/real(time_partition+1)
  !enddo

  ! Assignment of expected ionization fraction at the I-front
  do i_partition=1,time_partition
    xfinal(i_partition) = background_HII + (1-background_HII)*real(i_partition)/real(time_partition+1)
  enddo

  !write(*,*)'xfinal ',xfinal
  do  

     ! This happens once only
     if ( (inside_pos.ge.1) .and. (equilibrium_pos.eq.0) ) equilibrium_pos=1

     ! Reset equilibrium position and output_time_is_reached after each output
     if ( (inside_pos.ge.1) .and. (output_time_is_reached.eqv..true.) ) then
        equilibrium_pos=1
        output_time_is_reached=.false.
     endif

     ! Write output
     if (abs(front_sim_time-next_output_time) .le. 1e-6*front_sim_time) then
        call output(front_sim_time,front_dt,end_time,powerlaw)
        next_output_time=next_output_time+output_time
     endif

     ! Reset front_index to one and front_pos goes further out
     if (front_index.gt.time_partition) then
       front_index=1
       front_pos=front_pos+1
     endif
 
     ! If I-front is still inside the simulation box
     if (front_pos.le.mesh) then

       ! To get the adaptive timestep for the front cell
       !call time_step(front_dt, front_pos, front_index, xfinal)

       ! When the timestep is too large that traverse the output time
       !if ((next_output_time-front_sim_time).le.front_dt) then 

         ! Assign timestep to the front so that it hits the output time 
         !front_actual_dt= next_output_time-front_sim_time 
         !output_time_is_reached=.true. 

       ! When the timestep is smaller so it does not traverse the output time
   	   !else 

         !Assign the adaptive timestep to the front
         !front_actual_dt=front_dt
         !front_index=front_index+1
       !endif
front_actual_dt= next_output_time-front_sim_time
output_time_is_reached=.true.
     ! If I-front is out of the simulation box
     else 

       ! Assign timestep to all the cells so that it hits the output time 
       front_actual_dt=next_output_time-front_sim_time 
       output_time_is_reached=.true.

     endif
!front_actual_dt= 0.5*front_actual_dt     
      ! Report time and time step
      write(logf,'(A,2(1pe10.3,1x),A)') 'Time, dt:', &
        front_sim_time/YEAR,front_actual_dt/YEAR,' (years)'

      ! For cosmological simulations evolve proper quantities
      if (cosmological) then
        call redshift_evol(front_sim_time+0.5*front_actual_dt)
        call cosmo_evol()
      endif

      ! Evolve cells from 1 to inside_pos using inside dt
      if ( (equilibrium_pos.le.inside_pos) .and. (inside_pos.ge.1) .and. (equilibrium_pos.ge.1)) then
        call evolve1D(next_output_time-front_sim_time,min(mesh,equilibrium_pos),min(mesh,inside_pos)) 
        equilibrium_pos=inside_pos+1
      endif
     
      ! Evolve cells from inside_pos+1 to mesh using front dt
       
      call evolve1D(front_actual_dt,inside_pos+1,mesh)

      ! Check if the simulation finishes
      if (abs(front_sim_time-end_time) .lt. 1e-6*end_time) exit  

      ! Update time
      front_sim_time=front_sim_time+front_actual_dt 

      ! Find out the position of local thermal equilirium
      ! If in equilibrium, update 
      !if (front_index.gt.size(xfinal) .and. front_pos.gt.1) then
      !  call LTE_calculation(inside_pos,front_pos,next_output_time-front_sim_time)
      !endif
 
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                    !!
!!                       End time is reached                          !!
!!                                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Scale to the current redshift
  if (cosmological) then
     call redshift_evol(front_sim_time)
     call cosmo_evol()
  endif

  ! Clean up some stuff
  call close_down ()

  ! Update clock counters (cpu + wall, to avoid overflowing the counter)
  call update_clocks ()

  ! Report clocks (cpu and wall)
  call report_clocks ()

  ! End the run
  call mpi_end () 

end Program C2Ray
