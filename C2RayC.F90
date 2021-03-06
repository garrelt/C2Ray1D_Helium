! Constant timestep scheme (Hydrogen+Helium version)
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
  use times, only: time_ini, end_time, max_dt, output_time
  use timestep, only: time_step
  use sizes, only: mesh

#ifdef XLF
  ! Modules for the xlf (IBM) compiler
  USE XLFUTILITY, only: iargc, getarg, flush => flush_
#endif

  implicit none

  ! CPU time variables
  real(kind=dp) :: sim_time                     !< actual time seen by I-front
  real(kind=dp) :: next_output_time             !< time of next output
  real(kind=dp) :: actual_dt 	                !< timestep taken by I-front
  character(len=512) :: inputfile	        !< Input file

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
  call rad_ini() 

  ! Initialize source property        		 
  call source_properties_ini (sourcetype) 

  ! Initialize time step parameters
  call time_ini()        		

  ! Initialize actual time to zero
  sim_time=0.0     

  ! An output will be made by that time           		
  next_output_time=0.0 

  ! Initialize cosmology 
  if (cosmological) then  		    
     call cosmology_init(zred00,sim_time)
     call redshift_evol(sim_time)
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
  
  open(unit=33,file='hihi.dat',form='formatted',status='unknown')

  do  

     ! Write output
     if (abs(sim_time-next_output_time) .le. 1e-6*sim_time) then
       call output(sim_time,max_dt,end_time)
       next_output_time=next_output_time+output_time
     endif

     ! Assignment of timestep
     actual_dt=min(next_output_time-sim_time,max_dt)

     ! Report time and time step
     write(logf,'(A,2(1pe10.3,1x),A)') 'Time, dt:', &
           sim_time/YEAR,actual_dt/YEAR,' (years)'

     ! For cosmological simulations evolve proper quantities
     if (cosmological) then
       call redshift_evol(sim_time+0.5*actual_dt)
       call cosmo_evol()
     endif

     ! Take one time step
     call evolve1D(actual_dt,1,mesh) 

     ! Check if the simulation finishes
     if (abs(sim_time-end_time).lt.1e-6*end_time) exit     

     ! Update time
     sim_time=sim_time+actual_dt 

  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                    !!
!!                       End time is reached                          !!
!!                                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Scale to the current redshift
  if (cosmological) then 
     call redshift_evol(sim_time)
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

  close(33)

end Program C2Ray
