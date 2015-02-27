!>
!! \brief This module contains routines for file output

module output_module
  
  ! This file contains routines having to do with the output
  ! of Ifront programs.
  
  ! setup_out : open files
  ! close_down : close files
  ! output : write output

  use precision, only: dp
  use my_mpi, only: rank
  use file_admin, only: stdinput, results_dir, file_input, logf
  use tped
  use abundances, only: abu_h, abu_he
  use radiation, only: L_0
  use cgsconstants, only: hplanck, temphe
  use cgsphotoconstants, only: ion_freq_HI,ion_freq_HeI,ion_freq_HeII
  use radiation, only: NumFreqBnd,NumBndin1,NumBndin2,NumBndin3,&
                       sigma_HI,sigma_HeI,sigma_HeII,freq_max,freq_min	,pl_input_flux			
use grid, only: dr, x, r_in
		  use times, only: end_time

  implicit none
  
  integer,parameter :: max_input_streams=5 !< maximum number of output streams
  integer,dimension(max_input_streams) :: streams !< flags for output streams
  
contains
  !----------------------------------------------------------------------------

  !> Initializes global output (photonstatistics, analytical solution)
  !!
  !! Analytical: ASCII table
  !! \li time (s)
  !! \li Numerical position (cm)
  !! \li Analytical position (cm)
  !! \li Relative error between analytical and numerical
  !! \li Uncertainty due to front width
  !! \li Error due to finite cell size.
  !!
  !! PhotonCounts: ASCII table
  !! \li time
  !! \li Number of (ionizations + recombinations) / photons during time step
  !! \li Number of ionizations /(ionizations + recombinations)
  !! \li Number of recombinations /(ionizations + recombinations)
  !! \li Number of (ionizations + recombinations) / photons 
  !! \li Number of (ionizations + recombinations) / photons since t=0

  subroutine setup_output ()
    
    ! Sets up output
    
    ! photon statistics
    use photonstatistics, only: do_photonstatistics, &
         initialize_photonstatistics

    if (rank == 0) then
       ! Open files
       if (do_photonstatistics) open(unit=90,file=trim(adjustl(results_dir))//'PhotonCounts.out', &
            form='formatted',status='unknown')
       open(unit=91,file=trim(adjustl(results_dir))//'Analytical.out',form='formatted',status='unknown')
       write(91,*) 'Time - Numerical Solution - Analytical solution -', &
            'Fractional Error'

       open(unit=50,file='Outputtimes.out',status='unknown')
       open(unit=52,file='Inputnumbers.out',status='unknown')

    endif
    if (do_photonstatistics) call initialize_photonstatistics ()
    
  end subroutine setup_output
  
  !-----------------------------------------------------------------------------

  !> Closes global output files which have been open the entire run

  subroutine close_down ()
    
    ! Closes down
    
    if (rank == 0) then
       close(90)
       close(91)
       close(50)
       close(52)
    endif

  end subroutine close_down
  
  !----------------------------------------------------------------------------

  !> produces output for a time frame. See below for format
  !!
  !! Output format: ASCII table
  !! \li r-position
  !! \li neutral hydrogen fraction
  !! \li ionized hydrogen fraction
  !! \li temperature
  !! \li density
  !!
  !! No time information is available in the output, but the
  !! file names are Ifront_xx.xxx.dat, where xx.xxx is the
  !! fraction of the simulation time passed (initial condition 0.000,
  !! last output 1.000).
 
  subroutine output (time,dt,end_time,powerlaw) !* (step time dt end_time)

    ! Simple output routine.

    !      Output format:
    !      The output stream (see setup_output)
    !      stream 1: r-position
    !                neutral hydrogen fraction
    !                ionized hydrogen fraction
    !                temperature
    !                density
    !                neutral helium fraction   
    !                single ionized helium fraction
    !                double ionized helium fraction  
    !      Data from subsequent time steps are stacked without
    !      any separators. No time information is available in
    !      the output.
    
    !      Analytical: time (s)
    !                  Numerical position (cm)
    !                  Analytical position (cm)
    !                  Relative error between analytical and numerical
    !                  Uncertainty due to front width
    !                  Error due to finite cell size.
    
    !      PhotonCounts: time
    !                    Number of (ionizations + recombinations) / photons 
    !                     during time step
    !                    Number of ionizations /(ionizations + recombinations)
    !                    Number of recombinations /(ionizations + recombinations)
    !                    Number of (ionizations + recombinations) / photons 
    !                    Number of (ionizations + recombinations) / photons 
    !                     since t=0

    use sizes, only: mesh
    use grid, only: x,dr
    use material
    use photonstatistics
    use radiation, only: r_star,l_star,s_star,pl_photo_flux_wanted,source_ionzing_photon_rate
	use evolve, only: alpha_HI_B
  use mathconstants, only: pi
  use abundances, only: abu_h
  
  implicit none
  
    real(kind=dp),intent(in) :: time !< current simulation time
    real(kind=dp),intent(in) :: dt !< time step taken
    real(kind=dp),intent(in) :: end_time !< end simulation at this time
		    real(kind=dp), intent(in) :: powerlaw
    !integer, intent(in)      :: step, number_timesteps   
    real(kind=dp),dimension(1:mesh) :: R_HII_array
	real(kind=dp) :: R_HII
	real(kind=dp) ::alphaB
	real(kind=dp) :: meanT
	integer :: number
    integer :: i,j,k,ns,v
    character(len=6) :: zred_str
    character(len=40) :: file0,file1,file2,file3,file4,file5,file6
    real(kind=dp) :: totalsrc,tau
    logical crossing,recording_photonstats
	real(kind=dp) :: max_temperature
    real(kind=dp) :: ana_front,num_front,num_front1,num_front2
	real(kind=dp) :: column_HI_this_0100,column_HeI_this_0100,column_HeII_this_0100
	real(kind=dp) :: column_HI_this_0200,column_HeI_this_0200,column_HeII_this_0200
	real(kind=dp) :: column_HI_this_0300,column_HeI_this_0300,column_HeII_this_0300
	real(kind=dp) :: column_HI_this_0400,column_HeI_this_0400,column_HeII_this_0400	
	real(kind=dp) :: column_HI_this_0500,column_HeI_this_0500,column_HeII_this_0500
	real(kind=dp) :: column_HI_this_0700,column_HeI_this_0700,column_HeII_this_0700
	real(kind=dp) :: column_HI_this_0900,column_HeI_this_0900,column_HeII_this_0900
	real(kind=dp) :: sigma_HI_this,sigma_HeI_this,sigma_HeII_this	
	real(kind=dp) :: tau_HI_this_0100,tau_HeI_this_0100,tau_HeII_this_0100
	real(kind=dp) :: tau_HI_this_0200,tau_HeI_this_0200,tau_HeII_this_0200
	real(kind=dp) :: tau_HI_this_0300,tau_HeI_this_0300,tau_HeII_this_0300
	real(kind=dp) :: tau_HI_this_0400,tau_HeI_this_0400,tau_HeII_this_0400
	real(kind=dp) :: tau_HI_this_0500,tau_HeI_this_0500,tau_HeII_this_0500
	real(kind=dp) :: tau_HI_this_0700,tau_HeI_this_0700,tau_HeII_this_0700
	real(kind=dp) :: tau_HI_this_0900,tau_HeI_this_0900,tau_HeII_this_0900
	real(kind=dp) :: a,b
	logical :: no_match
	integer :: HeIII_zone_pos
	real(kind=dp) :: analytical_HeIII_zone_size
	real(kind=dp) :: improved_analytical_HeIII_zone_size
	integer :: eight_Mpc_pos
	real(kind=dp) :: eight_Mpc_HI_column_density
	real(kind=dp) :: eight_Mpc_HI_tau
	real(kind=dp) :: eight_Mpc_HeI_column_density
	real(kind=dp) :: eight_Mpc_HeI_tau
	real(kind=dp) :: HeIII_size_HeII_column_density
	real(kind=dp) :: HeIII_size_HeII_tau
	real(kind=dp) :: lambda
	real(kind=dp) :: alphaHeIIIA, alphaHeIIIB
	real(kind=dp) :: a_factor, b_factor, ndot
	
	
    if (rank == 0) then
       ! Stream 1
       write(file1,'(f6.3)') real(time)/real(end_time)
       write(file2,'(f6.3)') real(time)/real(end_time)
       !write(*,'(f6.3)') real(time)/real(end_time)

	   !do i=1,mesh
		!   if (HI_photoionization_rate(i).ge.gamma_uvb(1)) then
		!	   HI_photoionization_rate(i)=0.0
		 !  else 
		!	   HI_photoionization_rate(i)=1.0
		 !  endif
		!enddo

		!R_HII = ((3.0*source_ionzing_photon_rate)/(4.0*pi*alpha_HI_B*(ndens(mesh,1,1)*abu_h)**2.0))**(1.0/3.0)
		!write(*,*)'flux is ',pl_photo_flux_wanted
		do i=1,mesh
			if (x(i).ge.R_HII) then
				R_HII_array(i)=1.0
			else
				R_HII_array(i)=0.0
			endif
		enddo
		!write(*,*) R_HII/x(1)
if (abs(real(time)/real(end_time)).le.1.0001 .and. abs(real(time)/real(end_time)).ge.0.9999) then
	
	
	
! write the tau of first cell
!write(*,*) 'col den of DLA is (10^20)',ndens(1,1,1)*xHI(1)*abu_h*dr(1)/1e20_dp
!write(*,*) 'xHI = ', xHI(1)
!write(*,*) 'xHeIII = ', xHeIII(1)
file0=trim(adjustl(results_dir))//'spectrum.dat'
open(unit=50,file=file0,form='formatted',status='unknown')
do v = 1,10000
	
	no_match=.true.
	! to get the sigma of HI, HeI and HeII
	do j = 1,NumFreqBnd-1
		if (ion_freq_HI*v .ge. freq_min(j) .and. ion_freq_HI*v.lt.freq_max(j)) then
  		  a = ion_freq_HI*v - freq_min(j)
		  b = freq_max(j) - freq_min(j)
		  sigma_HI_this = (a*sigma_HI(j+1) + b*sigma_HI(j))/(a+b)
		  sigma_HeI_this = (a*sigma_HeI(j+1) + b*sigma_HeI(j))/(a+b)
		  sigma_HeII_this = (a*sigma_HeII(j+1) + b*sigma_HeII(j))/(a+b)
		  		  
		  no_match = .false.	
	  	endif	
	enddo
	
 	if (ion_freq_HI*v.ge.freq_max(NumFreqBnd-1) .and. no_match.eqv..true.) then
 	   a = ion_freq_HI*v - freq_min(NumFreqBnd-1)
 	   b = ion_freq_HI*v - freq_max(NumFreqBnd-1)
	   sigma_HI_this = max(0.0,((a+b)*sigma_HI(NumFreqBnd)-b*sigma_HI(NumFreqBnd-1))/a)
	   sigma_HeI_this = max(0.0,((a+b)*sigma_HeI(NumFreqBnd)-b*sigma_HeI(NumFreqBnd-1))/a)
	   sigma_HeII_this = max(0.0,((a+b)*sigma_HeII(NumFreqBnd)-b*sigma_HeII(NumFreqBnd-1))/a)	   	     				
	endif	

	column_HI_this_0100=0.0
	column_HI_this_0200=0.0
	column_HI_this_0300=0.0
	column_HI_this_0400=0.0
	column_HI_this_0500=0.0
	column_HI_this_0700=0.0
	column_HI_this_0900=0.0
	
	column_HeI_this_0100=0.0
	column_HeI_this_0200=0.0
	column_HeI_this_0300=0.0
	column_HeI_this_0400=0.0
	column_HeI_this_0500=0.0
	column_HeI_this_0700=0.0
	column_HeI_this_0900=0.0
	
	column_HeII_this_0100=0.0
	column_HeII_this_0200=0.0
	column_HeII_this_0300=0.0
	column_HeII_this_0400=0.0
	column_HeII_this_0500=0.0
	column_HeII_this_0700=0.0
	column_HeII_this_0900=0.0
			
	do i=1,mesh
		
		if (i.le.100) then
			column_HI_this_0100=column_HI_this_0100+ndens(i,1,1)*xHI(i)*abu_h*dr(1)
			column_HeI_this_0100=column_HeI_this_0100+ndens(i,1,1)*xHI(i)*abu_he*dr(1)
			column_HeII_this_0100=column_HeII_this_0100+ndens(i,1,1)*xHI(i)*abu_he*dr(1)
			tau_HI_this_0100=column_HI_this_0100*sigma_HI_this
			tau_HeI_this_0100=column_HeI_this_0100*sigma_HeI_this
			tau_HeII_this_0100=column_HeII_this_0100*sigma_HeII_this				
		endif
		
		if (i.le.200) then
			column_HI_this_0200=column_HI_this_0200+ndens(i,1,1)*xHI(i)*abu_h*dr(1)
			column_HeI_this_0200=column_HeI_this_0200+ndens(i,1,1)*xHI(i)*abu_he*dr(1)
			column_HeII_this_0200=column_HeII_this_0200+ndens(i,1,1)*xHI(i)*abu_he*dr(1)
			tau_HI_this_0200=column_HI_this_0200*sigma_HI_this
			tau_HeI_this_0200=column_HeI_this_0200*sigma_HeI_this
			tau_HeII_this_0200=column_HeII_this_0200*sigma_HeII_this				
		endif
		
		if (i.le.300) then
			column_HI_this_0300=column_HI_this_0300+ndens(i,1,1)*xHI(i)*abu_h*dr(1)
			column_HeI_this_0300=column_HeI_this_0300+ndens(i,1,1)*xHI(i)*abu_he*dr(1)
			column_HeII_this_0300=column_HeII_this_0300+ndens(i,1,1)*xHI(i)*abu_he*dr(1)
			tau_HI_this_0300=column_HI_this_0300*sigma_HI_this
			tau_HeI_this_0300=column_HeI_this_0300*sigma_HeI_this
			tau_HeII_this_0300=column_HeII_this_0300*sigma_HeII_this				
		endif
		
		if (i.le.400) then
			column_HI_this_0400=column_HI_this_0400+ndens(i,1,1)*xHI(i)*abu_h*dr(1)
			column_HeI_this_0400=column_HeI_this_0400+ndens(i,1,1)*xHI(i)*abu_he*dr(1)
			column_HeII_this_0400=column_HeII_this_0400+ndens(i,1,1)*xHI(i)*abu_he*dr(1)
			tau_HI_this_0400=column_HI_this_0400*sigma_HI_this
			tau_HeI_this_0400=column_HeI_this_0400*sigma_HeI_this
			tau_HeII_this_0400=column_HeII_this_0400*sigma_HeII_this				
		endif
		
		if (i.le.500) then
			column_HI_this_0500=column_HI_this_0500+ndens(i,1,1)*xHI(i)*abu_h*dr(1)
			column_HeI_this_0500=column_HeI_this_0500+ndens(i,1,1)*xHI(i)*abu_he*dr(1)
			column_HeII_this_0500=column_HeII_this_0500+ndens(i,1,1)*xHI(i)*abu_he*dr(1)
			tau_HI_this_0500=column_HI_this_0500*sigma_HI_this
			tau_HeI_this_0500=column_HeI_this_0500*sigma_HeI_this
			tau_HeII_this_0500=column_HeII_this_0500*sigma_HeII_this				
		endif
		
		if (i.le.700) then
			column_HI_this_0700=column_HI_this_0700+ndens(i,1,1)*xHI(i)*abu_h*dr(1)
			column_HeI_this_0700=column_HeI_this_0700+ndens(i,1,1)*xHI(i)*abu_he*dr(1)
			column_HeII_this_0700=column_HeII_this_0700+ndens(i,1,1)*xHI(i)*abu_he*dr(1)
			tau_HI_this_0700=column_HI_this_0700*sigma_HI_this
			tau_HeI_this_0700=column_HeI_this_0700*sigma_HeI_this
			tau_HeII_this_0700=column_HeII_this_0700*sigma_HeII_this				
		endif
		
		if (i.le.900) then
			column_HI_this_0900=column_HI_this_0900+ndens(i,1,1)*xHI(i)*abu_h*dr(1)
			column_HeI_this_0900=column_HeI_this_0900+ndens(i,1,1)*xHI(i)*abu_he*dr(1)
			column_HeII_this_0900=column_HeII_this_0900+ndens(i,1,1)*xHI(i)*abu_he*dr(1)
			tau_HI_this_0900=column_HI_this_0900*sigma_HI_this
			tau_HeI_this_0900=column_HeI_this_0900*sigma_HeI_this
			tau_HeII_this_0900=column_HeII_this_0900*sigma_HeII_this				
		endif
		
	enddo
	
write(50,*)v*ion_freq_HI,&
L_0*(real(v)**(-powerlaw)),&
L_0*(real(v)**(-powerlaw))*exp(-tau_HI_this_0100-tau_HeI_this_0100-tau_HeII_this_0100),&
L_0*(real(v)**(-powerlaw))*exp(-tau_HI_this_0200-tau_HeI_this_0200-tau_HeII_this_0200),&
L_0*(real(v)**(-powerlaw))*exp(-tau_HI_this_0300-tau_HeI_this_0300-tau_HeII_this_0300),&
L_0*(real(v)**(-powerlaw))*exp(-tau_HI_this_0400-tau_HeI_this_0400-tau_HeII_this_0400),&
L_0*(real(v)**(-powerlaw))*exp(-tau_HI_this_0500-tau_HeI_this_0500-tau_HeII_this_0500),&
L_0*(real(v)**(-powerlaw))*exp(-tau_HI_this_0700-tau_HeI_this_0700-tau_HeII_this_0700),&
L_0*(real(v)**(-powerlaw))*exp(-tau_HI_this_0900-tau_HeI_this_0900-tau_HeII_this_0900)
enddo
endif

       file1=trim(adjustl(results_dir))//'Ifront1_'//trim(adjustl(file1))//'.dat'
       open(unit=51,file=file1,form='formatted',status='unknown')
       do i=1,mesh
          write(51,'(8(1pe11.4,1x))') x(i),xHI(i),xHII(i), temper(i),ndens(i,1,1), &
          HI_photoionization_rate(i),recombination_rate(i),xHeIII(i)
       enddo

       file2=trim(adjustl(results_dir))//'python_'//trim(adjustl(file2))//'.dat'
       open(unit=53,file=file2,form='formatted',status='unknown')
	   
	   ! for IDL
       !do i=1,mesh
       !   write(53,'(8(1pe11.4,1x))') x(i),xHI(i),xHII(i),temper(i),ndens(i,1,1),&
	!	     HI_photoionization_rate(i),recombination_time(i),xHeIII(i)
    !   enddo
	   
	   ! for me
       do i=1,mesh
          write(53,'(8(1pe11.4,1x))') x(i),xHI(i),temper(i),ndens(i,1,1)*(1-abu_he),&
		     xHeIII(i),ndens(i,1,1)*abu_he
       enddo
	   	   
	   if (abs(time-end_time)/end_time .le. 0.0001) then
		   
		   meanT=0.0
		   number=0
		   do i=1,mesh
		   if (xHeIII(i).ge.0.5) then
			   meanT=meanT+temper(i)
			   number=number+1
			 endif  
		 enddo
		 meanT=meanT/real(number)
	   write(*,*) 'mean T is ', meanT
	   write(*,*) 'max T is ', maxval(temper), ' at position ', number-1
	   write(*,*) 'max T is ', (x(number-1)+x(number))*3.24077929e-25/2.0, 'pMpc'   
	   ! to report HeIII zone size in pMpc
	   do i=1,mesh-1
		   if (xHeIII(i).ge.0.5 .and. xHeIII(i+1).le.0.5) then
			 HeIII_zone_pos = i
			 exit
		   endif
       enddo
	   
write(*,*) 'HeIII zone size is ', (x(HeIII_zone_pos)+x(HeIII_zone_pos+1))*3.24077929e-25/2.0, 'pMpc'
! J1148 - 6.5e57
! J1030 - 8.0e57
! UKIDS4 - 0.4e57
ndot = pl_input_flux



 analytical_HeIII_zone_size = (3.0*ndot*4.0**(-powerlaw+1)*end_time/ &
     (4.0*3.1416*ndens(i,1,1)*abu_he))**(1.0/3.0)
 analytical_HeIII_zone_size = analytical_HeIII_zone_size*3.24077929e-25
 
write(*,*) 'anal. HeIII zone size is ', analytical_HeIII_zone_size, 'pMpc' 

lambda=2.0_dp*(temphe(1)/meanT) 
alphaHeIIIA=2.538e-13*lambda**1.503_dp/(1.0_dp+(lambda/0.522_dp)**0.470_dp)**1.923_dp
alphaHeIIIB=5.5060e-14_dp*lambda**1.5_dp/(1.0_dp+(lambda/2.740_dp)**0.407_dp)**2.242_dp

a_factor = ndot*4.0**(-powerlaw+1)/(4.0*3.1416*ndens(i,1,1)*abu_he)
b_factor = (ndens(i,1,1)*abu_h+2.0*ndens(i,1,1)*abu_he)*alphaHeIIIA/3.0

improved_analytical_HeIII_zone_size = (a_factor/b_factor)**(1.0/3.0)* &
     (1-exp(-3.0*b_factor*end_time))**(1.0/3.0)
improved_analytical_HeIII_zone_size = improved_analytical_HeIII_zone_size* &
              3.24077929e-25


write(*,*) 'im. anal. HeIII zone size is ', improved_analytical_HeIII_zone_size, 'pMpc' 	 
	 
	 
	   ! to report column density of HI up to HeIII zone size
	   eight_Mpc_pos = ((3.0/3.24077929e-25)-r_in)/dr(1)+0.5
	   eight_Mpc_HI_column_density = 0
	   eight_Mpc_HeI_column_density = 0
	   
       do i=1,eight_Mpc_pos
		   eight_Mpc_HI_column_density = eight_Mpc_HI_column_density + ndens(i,1,1)*xHI(i)*abu_h*dr(1)
		   eight_Mpc_HeI_column_density = eight_Mpc_HeI_column_density + ndens(i,1,1)*xHeI(i)*abu_he*dr(1)   
       enddo
	   eight_Mpc_HI_tau = eight_Mpc_HI_column_density*sigma_HI(1)
	   eight_Mpc_HeI_tau = eight_Mpc_HeI_column_density*sigma_HeI(2)
	   
write(*,*) 'eight_Mpc_HI_tau is ',eight_Mpc_HI_tau
write(*,*) 'eight_Mpc_HeI_tau is ',eight_Mpc_HeI_tau
			  
       HeIII_size_HeII_column_density = 0
	   do i=1,(HeIII_zone_pos)
		   HeIII_size_HeII_column_density = HeIII_size_HeII_column_density +  &
		   ndens(i,1,1)*xHeII(i)*abu_he*dr(1)
	   enddo 
	   HeIII_size_HeII_tau = HeIII_size_HeII_column_density*sigma_HeII(28)
write(*,*) 'HeIII_size_HeII_tau is ',HeIII_size_HeII_tau

   endif
       ! Check if we are tracking photon conservation
       !if (do_photonstatistics .and. time > 0.0) then
       !   ! Photon Statistics
       !   total_ion=total_ion!+photon_loss
       !   grtotal_ion=grtotal_ion+total_ion-totcollisions
       !   if (time > 0.0) then
       !      write(90,'(7(1pe13.5))') &
       !            real(time)/end_time, &
       !           (total_ion-totcollisions)/(s_star*dt), &
       !            total_ion, &
       !            dh0/dt, &
       !            dh1/dt, &
       !            dhe0/dt, &
       !            dhe1/dt, &
       !            dhe2/dt, &
       !            photon_loss, &
       !            (2.0*dhe0+dhe1+dh0+totrec)/dt
 !                  (dh0+dh1+dhe0+dhe1+dhe2)/dt
 !                  totcollisions
 !                 grtotal_ion/(s_star*time)
		   
                   

        !  endif
       !endif

       ! Compare to analytical solution (ejr, 04072004)
!       if (time > 0.0) then     
!          call calc_ana_front(ana_front,time)
!          call calc_num_front(num_front1,0.1d0)
!          call calc_num_front(num_front,0.5d0)
!          call calc_num_front(num_front2,0.9d0)
!          if (num_front.eq.0.0) num_front=1.0
!          if (ana_front.eq.0.0) ana_front=1.0
!          write(91,'(6(1pe10.3,1x))') time,num_front, &
!               ana_front, &
!               (num_front-ana_front)/ana_front, &
!               (num_front1-num_front2)/num_front, &
!               dr/num_front
!       endif
      
    endif
    return
  end subroutine output

  !---------------------------------------------------------------------------

  !> Calculate the location of the ionization front analytically\n
  !! \b Authors: Erik-Jan Rijkhorst, Ilian Iliev, Garrelt Mellema\n
  !! \b Date: 20-Oct-2004 (5-Aug-2004)\n
  !! This routine needs several external functions:
  !!  - LambertW(x): the LambertW function.(real argument and ouput)
  !!  - expint(n,x): the exponentional integral of order n

  subroutine calc_ana_front(ana_front,time)
    !     
    !     Calculate the location of the ionization front.
    !     
    !     Authors: Erik-Jan Rijkhorst, Ilian Iliev, Garrelt Mellema
    !     Date: 20-Oct-2004 (5-Aug-2004)
    
    !     This routine needs several external functions:
    !     LambertW(x): the LambertW function.(real argument and ouput)
    !     expint(n,x): the exponentional integral of order n

    use cosmology
    use material
    use radiation
    use grid
    
    real(kind=dp),intent(in) :: time !< current simulation time
    real(kind=dp),intent(out) :: ana_front !< front position
    
    real(kind=dp) :: StromgrenRadius
    real(kind=dp) :: typical_dens
    real(kind=dp) :: L,K,a,S1,t_reccore
    real(kind=dp) :: tratio
    !real(kind=dp) :: LambertW
    !real(kind=dp) :: expint     

    ! The analytical solution depends on the input:
    ! test 1: uniform density
    ! test 2: 1/r density
    ! test 3: constant density for r<r_core, 1/r^2 for r>r_core,
    !          and L close to 0 (see below).
    ! test 4: uniform density in expanding universe


    select case (testnum)
    case(1) ! constant density solution
       typical_dens=ndens(1,1,1)
       StromgrenRadius = (3.0*S_star/ &
            (4.0*pi*typical_dens*typical_dens* &
            clumping*bh00))**(1.0/3.0)
       ana_front = StromgrenRadius*(1.0- &
            exp(-typical_dens*clumping*bh00*time))**(1.0/3.0)
    case(2)
       L=S_star/(4.0*pi*dens_core*r_core)
       K=dens_core*r_core*clumping*bh00
       ana_front=L/K*(1.0+LambertW(-exp(-K*K*time/L-1.0)))
    case(3)
       !write(*,*) S_star/(4.0*pi*dens_core*r_core*r_core), &
       !     4./3.*dens_core*r_core*clumping*bh00
       L=S_star/(4.0*pi*dens_core*r_core*r_core) &
            -4./3.*dens_core*r_core*clumping*bh00
       K=dens_core*r_core*r_core*clumping*bh00
       t_reccore=1./(dens_core*clumping*bh00)
       if(abs(L)/(4./3.*dens_core*r_core*clumping*bh00) < 1e-3) then
          ana_front=r_core*sqrt(1.+2.*time/t_reccore)
       else
          write(*,*) 'No analytical solution implemented for', &
               ' these parameters.'
          write(*,*) abs(L)/(4./3.*dens_core*r_core*clumping*bh00)
          ana_front=0.0
          ! The following analytical solution seems to be wrong??
          ! ana_front=K/L*(1.0+ &
          ! LambertW(-(L*r_core/K+1.)/exp(L*L*time/K+L*r_core/K+1.0d0)))
              
       endif
       typical_dens=ndens(1,1,1)
       StromgrenRadius = (3.0*S_star/ &
            (4.0*pi*typical_dens*typical_dens* &
            clumping*bh00))**(1.0/3.0)
       if (time < -t_reccore*log(1.-(r_core/StromgrenRadius)**3)) then
          ana_front = StromgrenRadius*(1.0- &
               exp(-typical_dens*clumping*bh00*time))**(1.0/3.0)
       endif

    case(4)
       StromgrenRadius = (3.0*S_star/ &
            (4.0*pi*dens_core*dens_core* &
            clumping*bh00))**(1.0/3.0) !comoving r_S 
       ! (Note: not necessarily reached!!)
       tratio=t0/(t0+time)
       ana_front = StromgrenRadius*(eta/(1.+zred_t0)**3* & !exp(eta*tratio)* 
            (expint(2,eta*tratio,tratio*eta)/tratio &
            -expint(2,eta,tratio*eta)))**(1./3.) &
            /(1.+zred)
    case default
       write(*,*) 'Unknown test problem'
       ana_front=0.0
    end select

  end subroutine calc_ana_front
      
  !---------------------------------------------------------------------------
      
  !> finds the the location of the numerical ionization front.\n
  !! \b Author: Erik-Jan Rijkhorst\n
  !! \b Date: 5-Aug-2004

  subroutine calc_num_front(num_front,xlimit)
    !     
    !     Calculate the location of the ionization front.
    !     
    !     Author: Erik-Jan Rijkhorst
    !     Date: 5-Aug-2004

    use grid, only: x,dr
    use material, only: xHII

    !> ionization fraction that defines front location
    real(kind=dp),intent(in) :: xlimit 
    real(kind=dp),intent(out) :: num_front !< front position
      
    integer :: i, i1, i2

    num_front = x(1)-0.5*dr(1)
    i=1
    i1=i
    i2=i
    do while(xHII(i).ge.xlimit)
       i1 = i
       i2 = i+1
       i=i+1
    enddo
       if(xHII(i1).eq.0.0.and.xHII(i2).eq.0.0) then
       num_front = x(1)-0.5*dr(1)
    else
!remember: problems if too less ionizing photons
       num_front = (xlimit-xHII(i1)) * (x(i1)-x(i2)) / &
            (xHII(i1)-xHII(i2)) + x(i1)
    endif
  end subroutine calc_num_front
  
  !---------------------------------------------------------------------------

  !> Calculates LambertW function for z
  !
  !      /* Lambert W function. 
  !      Was ~/C/LambertW.c written K M Briggs Keith dot Briggs at 
  !      bt dot com 97 May 21.  
  !      Revised KMB 97 Nov 20; 98 Feb 11, Nov 24, Dec 28; 99 Jan 13; 00 Feb 23; 
  !      01 Apr 09
  !      Translated to Fortran 77 by Garrelt Mellema, 04 Sep 20.
  
  !      Computes Lambert W function, principal branch.
  !      See LambertW1.c for -1 branch.
  
  !      Returned value W(z) satisfies W(z)*exp(W(z))=z
  !      test data...
  !      W(1)= 0.5671432904097838730
  !      W(2)= 0.8526055020137254914
  !      W(20)=2.2050032780240599705
  !      To solve (a+b*R)*exp(-c*R)-d=0 for R, use
  !      R=-(b*W(-exp(-a*c/b)/b*d*c)+a*c)/b/c
  
  !      Test: 
  !      gcc -DTESTW LambertW.c -o LambertW -lm && LambertW
  !      Library:
  !      gcc -O3 -c LambertW.c */
  
  !      double LambertW(const double z);
  !      const int dbgW=0;
  
  function LambertW(z)
    
    real(kind=dp) :: LambertW
    real(kind=dp),intent(in) :: z !< input value, > -0.367879
    integer i
    real(kind=dp) :: eps,em1
    parameter(eps=4.0e-16)
    parameter(em1=0.3678794411714423215955237701614608)
    real(kind=dp) :: p,e,t,w
    real(kind=dp) :: q,r,q2,q3
    if (z.lt.-em1) then
       write(*,*) "LambertW: bad argument ",z," exiting"
       return
    endif
    if (z.eq.0.0) then
       LambertW=0.0
    elseif (z.lt.-em1+1e-4) then !series near -em1 in sqrt(q)
       q=z+em1
       r=sqrt(q)
       q2=q*q
       q3=q2*q
       LambertW=-1.0 &
            +2.331643981597124203363536062168*r &
            -1.812187885639363490240191647568*q &
            +1.936631114492359755363277457668*r*q &
            -2.353551201881614516821543561516*q2 &
            +3.066858901050631912893148922704*r*q2 &
            -4.175335600258177138854984177460*q3 &
            +5.858023729874774148815053846119*r*q3 &
            -8.401032217523977370984161688514*q3*q !! error approx 1e-16
       
       ! initial approx for iteration... */
    else
       if (z.lt.1.0) then    !! /* series near 0 */
          p=sqrt(2.0* &
               (2.7182818284590452353602874713526625*z+1.0))
          w=-1.0+p*(1.0+p*(-0.333333333333333333333+ &
               p*0.152777777777777777777777))
       else 
          w=log(z)          !! /* asymptotic */
       endif
       if (z.gt.3.0) w=w-log(w) !! /* useful? */
       do i=0,19             !! /* Halley iteration */
          e=exp(w)
          t=w*e-z
          p=w+1.0
          t=t/(e*p-0.5*(p+1.0)*t/p)
          w=w-t
          if (abs(t)<eps*(1.0+abs(w))) then
             LambertW=w    !! /* rel-abs error */
             goto 100
          endif
          !     /* should never get here */
       enddo
       write(*,*) "LambertW: No convergence at z= ",z," exiting."
       write(*,*) abs(t),eps*(1.0+abs(w))
    endif
      
100 continue
  end function LambertW
  

  !> calculates exponential integral of order n for x

  FUNCTION expint(n,x,etatratio)
    

    REAL(kind=dp) :: expint

    INTEGER,intent(in) :: n !< order
    REAL(kind=dp),intent(in) :: x !< argument
    REAL(kind=dp),intent(in) ::etatratio !< asymptotic value

    integer,parameter :: MAXIT=100
    REAL(kind=dp),parameter :: EPS=1.e-7,FPMIN=1.e-30,EULER=.5772156649

    INTEGER :: i,ii,nm1
    REAL(kind=dp) :: a,b,c,d,del,fact,h,psi

    nm1=n-1
    ! print*,'expint called',n,x,etatratio
    ! print*,'args',n,x
    if(n.lt.0.or.x.lt.0..or.(x.eq.0..and.(n.eq.0.or.n.eq.1)))then
       pause 'bad arguments in expint'
    else if(n.eq.0)then
       expint=exp(-x)/x
    else if(x.eq.0.)then
       expint=1./nm1
    else if(x.gt.1.)then
       b=x+n
       c=1./FPMIN
       d=1./b
       h=d
       do i=1,MAXIT
          a=-i*(nm1+i)
          b=b+2.
          d=1./(a*d+b)
          c=b+a/c
          del=c*d
          h=h*del
          if(abs(del-1.).lt.EPS)then
             expint=h*exp(-x+etatratio)
             ! print*,'expint done',expint,x,etatratio
             return
          endif
       enddo
       pause 'continued fraction failed in expint'
    else
       if(nm1.ne.0)then
          expint=1./nm1
       else
          expint=-log(x)-EULER
       endif
       fact=1.
       do i=1,MAXIT
          fact=-fact*x/i
          if(i.ne.nm1)then
             del=-fact/(i-nm1)
          else
             psi=-EULER
             do ii=1,nm1
                psi=psi+1./ii
             enddo
             del=fact*(-log(x)+psi)
          endif
          expint=expint+del
          ! if(abs(del).lt.abs(expint)*EPS) return
          if(abs(del).lt.abs(expint)*EPS)then
             expint=expint*exp(etatratio)
             return
          end if
       enddo
       pause 'series failed in expint'
    endif
    return
  END FUNCTION expint
  !  (C) Copr. 1986-92 Numerical Recipes Software 

end module output_module
