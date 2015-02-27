!     This module contains data and routines which deal with radiative
!     effects. Its main part deal with photo-ionizing radiation, but it
!     also initializes other radiative properties, such as cooling (which
!     are contained in different modules).
!     It can be used in hydrodynamic or stand-alone radiative transfer 
!     calculations.

module radiation
  
  use precision, only: dp
  use my_mpi
  use file_admin, only: logf
  use mathconstants, only: pi
  use cgsconstants, only: sigma_SB, &                    ! Stefan-Boltzmann constant
                          hplanck, &                     ! Planck constant
                          k_B, &                         ! Boltzmann constant
                          two_pi_over_c_square           ! two times pi over c aquare
  use cgsphotoconstants, only: ion_freq_HI,&             ! HI ionization energy in frequency
                               ion_freq_HeI,&            ! HeI ionization energy in frequency
                               ion_freq_HeII,&           ! HeII ionization energy in frequency
                               sigma_HI_at_ion_freq,&    ! HI cross section at its ionzing frequency
                               sigma_HeI_at_ion_freq,&   ! HeI cross section at its ionzing frequency
                               sigma_HeII_at_ion_freq    ! HeII cross section at its ionzing frequency
  use astroconstants, only: R_SOLAR, &                   ! Solar radius
                            L_SOLAR                      ! Solar luminosity
  use romberg, only: scalar_romberg, &                   ! 1D integration function
                     scalar_romberg_for_pl_scaling, &                   ! 1D integration function
                     vector_romberg, &                   ! 1D integration subroutine
                     romberg_initialisation              ! Romberg initialisation procedure
  use c2ray_parameters, only: T_eff_nominal,&            ! Black body  effective temperature for for nominal SED
                              S_star_nominal, &          ! Ionizing photon rate for for nominal SED
                              pl_index_nominal,&         ! Power law index for for nominal SED
                              EddLeff_nominal,&          ! Eddington efficiency for for nominal SED
                              EddLum                     ! Eddington luminosity for for nominal SED
  use material, only: isothermal

  use abundances, only: abu_h,abu_he
  implicit none

  integer,parameter :: NumFreq = 512      ! Number of integration points in each of the frequency bins
  integer,parameter :: NumFreq_for_pl_scaling = 4096
  integer,parameter :: NumTau = 2000      ! Number of table points for the optical depth
  integer,parameter :: NumBndin1 = 1      ! Number of frequency sub-bins in interval 1 
  integer,parameter :: NumBndin2 = 26     ! Number of frequency sub-bins in interval 2
  integer,parameter :: NumBndin3 = 20     ! Number of frequency sub-bins in interval 3
  integer,parameter :: NumFreqBnd=NumBndin1+NumBndin2+NumBndin3       ! Total number of frequency bins
  integer,parameter :: NumheatBin=NumBndin1+NumBndin2*2+NumBndin3*3   ! Total number of heating bins 

  ! Optical depths at the entrance of the grid.
  ! It can be used if radiation enters the simulation volume from the outside.
  real(kind=dp) :: boundary_tauHI = 0.0!4.0
  real(kind=dp) :: boundary_tauHeI = 0.0!4.0*(sigma_HeI_at_ion_freq/sigma_HI_at_ion_freq)*(abu_he/abu_h)
  real(kind=dp) :: boundary_tauHeII = 0.0

  ! Parameters defining the optical depth entries in the table.
  real(kind=dp),parameter :: minlogtau = -20.0                             ! Table position starts at log10(minlogtau) 
  real(kind=dp),parameter :: maxlogtau = 4.0                               ! Table position ends at log10(maxlogtau) 
  real(kind=dp),parameter :: dlogtau = (maxlogtau-minlogtau)/real(NumTau)  ! dlogtau is the step size in log10(tau)

  ! Some boring variables  
  real(kind=dp), dimension(1:3) :: CR1, CR2, bR1, dR1, aR2, bR2, y1R, y2R
  real(kind=dp) :: xeb


! for printing out quasar spectrum
real(kind=dp) :: L_0 



  ! Number of photon of the power law source
  real(kind=dp) :: pl_input_flux = 0.0

  ! Logical that determines the use of grey opacities
  logical,parameter :: grey = .false. 

  ! Type of source, B=black body, P=power law source
  Character :: sourcetype 

  ! Stellar properties
  real(kind=dp) :: T_eff        ! Black body effective temperature
  real(kind=dp) :: R_star       ! Black body radius
  real(kind=dp) :: L_star       ! Black body luminosity
  real(kind=dp) :: L_star_ion   ! Black body ionizing luminosity
  real(kind=dp) :: S_star       ! Black body ionizing photons rate

  ! Power law source properties
  real(kind=dp) :: pl_index             ! Power law index
  real(kind=dp) :: pl_minfreq           ! Minimum frequency for integration of total power
  real(kind=dp) :: pl_maxfreq           ! Maximum frequency for integration of total power
  real(kind=dp) :: pl_scaling           ! The scaling of the flux
  real(kind=dp) :: Edd_Efficiency       ! Eddinton efficieny
  real(kind=dp) :: source_ionzing_photon_rate  ! The rate of ionizing photon generated from the source
  real(kind=dp) :: pl_photo_flux_wanted
  real(kind=dp), dimension(:), allocatable :: delta_freq      ! Frequency width of integration 
  real(kind=dp), dimension(:), allocatable :: freq_max        ! Maximum freqeucny of integration 
  real(kind=dp), dimension(:), allocatable :: freq_min        ! Minimum freqeucny of integration

  ! Power law fit parameter for frequency range 1:3
  real(kind=dp), dimension(:), allocatable :: pl_index_cross_section_HI    ! Power law index of cross section of HI
  real(kind=dp), dimension(:), allocatable :: pl_index_cross_section_HeI   ! Power law index of cross section of HeI
  real(kind=dp), dimension(:), allocatable :: pl_index_cross_section_HeII  ! Power law index of cross section of HeII

  ! Cross section of atoms
  real(kind=dp), dimension(:), allocatable :: sigma_HI       ! Cross section of HI
  real(kind=dp), dimension(:), allocatable :: sigma_HeI      ! Cross section of HeI
  real(kind=dp), dimension(:), allocatable :: sigma_HeII     ! Cross section of HeII

  ! Parameters related to fraction of ionization and heating from different species
  real(kind=dp), dimension(:), allocatable :: f1ion_HI       ! Parameters related to ionization
  real(kind=dp), dimension(:), allocatable :: f1ion_HeI      ! Parameters related to ionization
  real(kind=dp), dimension(:), allocatable :: f1ion_HeII     ! Parameters related to ionization
  real(kind=dp), dimension(:), allocatable :: f2ion_HI       ! Parameters related to ionization
  real(kind=dp), dimension(:), allocatable :: f2ion_HeI      ! Parameters related to ionization
  real(kind=dp), dimension(:), allocatable :: f2ion_HeII     ! Parameters related to ionization
  real(kind=dp), dimension(:), allocatable :: f2heat_HI      ! Parameters related to heating
  real(kind=dp), dimension(:), allocatable :: f2heat_HeI     ! Parameters related to heating
  real(kind=dp), dimension(:), allocatable :: f2heat_HeII    ! Parameters related to heating
  real(kind=dp), dimension(:), allocatable :: f1heat_HI      ! Parameters related to heating
  real(kind=dp), dimension(:), allocatable :: f1heat_HeI     ! Parameters related to heating
  real(kind=dp), dimension(:), allocatable :: f1heat_HeII    ! Parameters related to heating
  
  ! Integrands ( frequency, optical depth )
  real(kind=dp),dimension(:,:), allocatable :: bb_photo_thick_integrand
  real(kind=dp),dimension(:,:), allocatable :: bb_photo_thin_integrand
  real(kind=dp),dimension(:,:), allocatable :: pl_photo_thick_integrand
  real(kind=dp),dimension(:,:), allocatable :: pl_photo_thin_integrand
  real(kind=dp),dimension(:,:), allocatable :: bb_heat_thick_integrand_HI
  real(kind=dp),dimension(:,:), allocatable :: bb_heat_thick_integrand_HeI
  real(kind=dp),dimension(:,:), allocatable :: bb_heat_thick_integrand_HeII
  real(kind=dp),dimension(:,:), allocatable :: bb_heat_thin_integrand_HI
  real(kind=dp),dimension(:,:), allocatable :: bb_heat_thin_integrand_HeI
  real(kind=dp),dimension(:,:), allocatable :: bb_heat_thin_integrand_HeII
  real(kind=dp),dimension(:,:), allocatable :: pl_heat_thick_integrand_HI
  real(kind=dp),dimension(:,:), allocatable :: pl_heat_thick_integrand_HeI
  real(kind=dp),dimension(:,:), allocatable :: pl_heat_thick_integrand_HeII
  real(kind=dp),dimension(:,:), allocatable :: pl_heat_thin_integrand_HI
  real(kind=dp),dimension(:,:), allocatable :: pl_heat_thin_integrand_HeI
  real(kind=dp),dimension(:,:), allocatable :: pl_heat_thin_integrand_HeII

  ! Integration table ( optical depth, sub-bin )
  real(kind=dp),dimension(:,:), target, allocatable :: bb_photo_thick_table  
  real(kind=dp),dimension(:,:), target, allocatable :: bb_photo_thin_table
  real(kind=dp),dimension(:,:), target, allocatable :: pl_photo_thick_table 
  real(kind=dp),dimension(:,:), target, allocatable :: pl_photo_thin_table
  real(kind=dp),dimension(:,:), target, allocatable :: bb_heat_thick_table  
  real(kind=dp),dimension(:,:), target, allocatable :: bb_heat_thin_table
  real(kind=dp),dimension(:,:), target, allocatable :: pl_heat_thick_table  
  real(kind=dp),dimension(:,:), target, allocatable :: pl_heat_thin_table

  ! photrates contains all the photo-ionization rates and heating rates
  type photrates    
     real(kind=dp) :: photo_cell_HI          ! HI photoionization rate of the cell    
     real(kind=dp) :: photo_cell_HeI         ! HeI photoionization rate of the cell    
     real(kind=dp) :: photo_cell_HeII        ! HeII photoionization rate of the cell    
     real(kind=dp) :: heat_cell_HI           ! HI heating rate of the cell       
     real(kind=dp) :: heat_cell_HeI          ! HeI heating rate of the cell    
     real(kind=dp) :: heat_cell_HeII         ! HeII heating rate of the cell          
     real(kind=dp) :: photo_in_HI            ! HI photoionization rate incoming to the cell    
     real(kind=dp) :: photo_in_HeI           ! HeI photoionization rate incoming to the cell
     real(kind=dp) :: photo_in_HeII          ! HeII photoionization rate incoming to the cell
     real(kind=dp) :: heat_in_HI             ! HI heating rate incoming to the cell
     real(kind=dp) :: heat_in_HeI            ! HeI heating rate incoming to the cell
     real(kind=dp) :: heat_in_HeII           ! HeII heating rate incoming to the cell 
     real(kind=dp) :: photo_out_HI           ! HI photoionization rate outgoing from the cell
     real(kind=dp) :: photo_out_HeI          ! HeI photoionization rate outgoing from the cell
     real(kind=dp) :: photo_out_HeII         ! HeII photoionization rate outgoing from the cell 
     real(kind=dp) :: heat_out_HI            ! HI heating rate outgoing from the cell
     real(kind=dp) :: heat_out_HeI           ! HeI heating rate outgoing from the cell
     real(kind=dp) :: heat_out_HeII          ! HeII heating rate outgoing from the cell
     real(kind=dp) :: heat                   ! Total heating rate of the cell
     real(kind=dp) :: photo_in               ! Total photoionization rate incoming to the cell
  end type photrates

  ! tablepos helps to locate correct position of the photoionization and heating tables
  type tablepos
    real(kind=dp), dimension(NumFreqBnd) :: tau            
    real(kind=dp), dimension(NumFreqBnd) :: odpos          
    real(kind=dp), dimension(NumFreqBnd) :: residual       
    integer, dimension(NumFreqBnd)       :: ipos           
    integer, dimension(NumFreqBnd)       :: ipos_p1        
  end type tablepos 

  logical,parameter :: do_secondary_ionization = .true.
  !logical,parameter :: do_secondary_ionization = .false.
#ifdef MPI       
    integer,private :: mympierror
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	  subroutine rad_ini (powerlaw)

	    real(kind=dp), intent(in) :: powerlaw
	                                                                                         !        
    ! Ask for the parameters of the spectrum                                                 !
    call spectrum_parms (powerlaw)                                                                   !      
                                                                                             ! 
    ! Initializes constants and tables for radiation processes                               ! 
    ! (heating, cooling and ionization)                                                      !
    call setup_scalingfactors ()                                                             !
                                                                                             !
    ! Initialize integration routines                                                        !
    call romberg_initialisation (NumFreq,NumFreq_for_pl_scaling)                                                    !
                                                                                             !
    ! Determine spectrum diagnostics                                                         !
    call spec_diag (powerlaw)                                                                        !
                                                                                             !
    ! Generate photoionization tables and heating tables                                     !
    call spec_integration ()                                                                 !                            
                                                                                             !
  end subroutine rad_ini                                                                     !                        
                                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Ask for the parameters of the spectrum
  subroutine spectrum_parms(powerlaw)

    use file_admin, only: stdinput, file_input
    
    integer :: i_choice                ! option number
    real(kind=dp) :: bb_total_flux    ! black body total flux
    real(kind=dp), intent(in) :: powerlaw
    ! Set total ionizing luminosity to zero (only changed by option 3)
    L_star_ion=0.0d0
    
    ! Set Eddington luminosity efficiency to zero
    Edd_Efficiency=0.0_dp

    ! Set power law index to 1.0
    pl_index=1.0_dp

    ! Ask for the input if you are processor 0 and the
    ! spectral parameters are not set in the c2ray_parameters
    ! Note that it is assumed that if teff_nominal is set, 
    ! S_star_nominal is ALSO set.
    if (rank == 0 .and. T_eff_nominal == 0.0) then

    	do while((sourcetype.ne.'B').and.(sourcetype.ne.'P'))

    	  if (.not.file_input) write(*,'(A,$)') 'Blackbody source (B) or power law source (P)'
    	  read(stdinput,*) sourcetype     ! Read source type, either blackbody or power law source

    	enddo

        !In blackbody case, ask for effective temperature of the source. 
        !The temperature should be bounded below by 2000 and above by 1000000.
        !And then ask for some parameters of the blackbody.
    	if (sourcetype.eq.'B') then

           write(logf,*) 'Black body source'
           T_eff=0.0
           do while (T_eff < 2000.0 .or. T_eff > 1000000.) 
              if (.not.file_input) write(*,'(A,$)') 'Give black body effective temperature (2000 <= T <= 1000000): '
              read(stdinput,*) T_eff      ! Read temperature of black body
              write(logf,*) 'Temperature of the black body : ', T_eff
              if (T_eff < 2000.0 .or. T_eff > 1000000.) then
                 write(*,*) 'Error: Effective temperature out of range. Try again'
              endif
           enddo

           ! Find total flux of blackbody(Stefan-Boltzmann law)
           bb_total_flux = sigma_SB*T_eff*T_eff*T_eff*T_eff

           ! Ask for radius, luminosity, ionizing luminosity or ionizing photon rate?
           if (.not.file_input) then
           write(*,'(A)') 'You can specify' 
           write(*,'(A)') ' 1) a stellar radius'
           write(*,'(A)') ' 2) a total luminosity'
           write(*,'(A)') ' 3) Total ionizing luminosity'
           write(*,'(A)') ' 4) Total number of ionizing photons'
           endif

       	   i_choice=0

           ! Report error if options are not 1, 2, 3 and 4
           do while (i_choice <= 0 .or. i_choice > 4)
              if (.not.file_input) write(*,'(A,$)') 'Preferred option (1, 2, 3 or 4): '
              read(stdinput,*) i_choice       ! Read option from the input, 1 to 4
              if (i_choice <= 0 .or. i_choice > 4) then
                 write(*,*) 'Error: Choose between 1 2 3 or 4'
              endif
           enddo

           select case (i_choice)
              
              case (1)
                write(logf,*) 'A stellar radius is specified'
          	if (.not.file_input) write(*,'(A,$)') 'Give radius in solar radius: '
          	read(stdinput,*) R_star      ! Read radius of the black body
                write(logf,*) 'The radius is ', R_star, ' solar radius'
          	R_star=R_star*r_solar
          	L_star=4.0d0*pi*R_star*R_star*bb_total_flux
          	S_star=0.0  ! Number of photo-ionizing photons is set to zero here but will be reset in spec_diag routine

              case (2)
                write(logf,*) 'A total luminosity is specified'
          	if (.not.file_input) write(*,'(A,$)') 'Give total luminosity in solar luminosity: '
          	read(stdinput,*) L_star      ! Read luminosity of the black body
                write(logf,*) 'The luminosity is ', L_star, ' solar luminosity'
          	L_star=L_star*l_solar
          	R_star=dsqrt(L_star/(4.0d0*pi*bb_total_flux))
          	S_star=0.0   ! Number of photo-ionizing photons is set to zero here but will be reset in spec_diag routine

              case (3)
                write(logf,*) 'Total ionizing luminosity is specified'
          	if (.not.file_input) write(*,'(A,$)') 'Give total ionizing luminosity in solar luminosity: '
          	read(stdinput,*) L_star_ion   ! Read ionizing luminosity of the black body
                write(logf,*) 'Total ionizing luminosity is ', L_star_ion, ' solar luminosty'
                ! Assign some fiducial values, these are overwritten in routine spec_diag
          	L_star_ion=L_star_ion*l_solar         	
          	R_star=r_solar                   
          	L_star=4.0d0*pi*R_star*R_star*bb_total_flux
          	S_star=0.0   ! Number of photo-ionizing photons is set to zero here but will be reset in spec_diag routine

              case (4)
                write(logf,*) 'Total number of ionzing photon is specified'
          	if (.not.file_input) write(*,'(A,$)') 'Give the number of ionizing photons per second: '
          	read(stdinput,*) S_star
                write(logf,*) 'The number of photons per second is ', S_star
          	! Assign some fiducial values for R_star and L_star, 
                ! these are scaled to correspond to S_star in routine spec_diag
          	R_star=r_solar
          	L_star=4.0d0*pi*R_star*R_star*bb_total_flux

           end select

        ! In power-law case, we ask for number of ionizing photons per second or Eddinton luminosity efficiency 
    	 elseif(sourcetype.eq.'P') then

           write(logf,*) 'Power law source'
           if (.not.file_input) then
           write(*,'(A)') 'You can specify'
     	   write(*,'(A)') ' 1)  Number of ionizing photons per second '
     	   write(*,'(A)') ' 2)  Efficiency parameter assuming a 1e6 solar mass BH' 
           endif
       	   i_choice=0

           ! Report error if options are not 1 and 2
           do while (i_choice <= 0 .or. i_choice > 2)
             if (.not.file_input) write(*,'(A,$)') 'Preferred option (1 or 2): '
             read(stdinput,*) i_choice
             if (i_choice <= 0 .or. i_choice > 2) then
               write(*,*) 'Error: Choose between 1 or 2'
         	endif
           enddo
 	      
           select case (i_choice)

             case (1)
               write(logf,*) 'Rate of ionizing photons is specified'
               if (.not.file_input) write(*,'(A,$)') 'give number of ionizing photons per second '                        	      
               read(stdinput,*) pl_input_flux          ! Read ionizing photons per second
               write(logf,*) 'The rate is ', pl_input_flux

               !Set the Eddinton luminosity efficiency by the default value 1.0
               Edd_Efficiency=EddLeff_nominal

             case (2)
               write(logf,*) 'Efficiency paramter is specified'
               if (.not.file_input)write(*,'(A,$)') 'give efficiency parameter '
               read(stdinput,*) Edd_Efficiency         ! Read Eddington efficiency
               write(logf,*) 'The efficiency parameter is ', Edd_Efficiency
      
           end select
          	
           if (.not.file_input) write(*,'(A,$)') 'power law index (for number of photons, not energy) '
           read(stdinput,*) pl_index      ! Read power law index, this number equal to one plus that of energy 
		              pl_index = powerlaw
           write(logf,*) 'Power law index is ', pl_index
           if (.not.file_input) write(*,'(A,$)') 'give lower and upper frequency limits in eV '
           read(stdinput,*) pl_minfreq,pl_maxfreq     ! Read lower and upper frequency limits in eV	
           write(logf,*) 'The lower energy limit is ', pl_minfreq, ' eV'
           write(logf,*) 'The upper energy limit is ', pl_maxfreq, ' eV'
           if (.not.file_input) write(*,'(A)') 'However, this is not implemented right now '          	

           ! set some fiducial values for the BB source here, though they are not useful
           R_star=r_solar
           S_star=0.0
           L_star=0.0
           T_eff=1.0e5

         endif

    ! In case neither blackbody nor power-law source
    else
       ! T_eff and S_star are assumed to have been set in the c2ray_parameter module
       T_eff=T_eff_nominal
       S_star=S_star_nominal
       bb_total_flux=sigma_SB*T_eff*T_eff*T_eff*T_eff   
       ! Assign some fiducial values, these are scaled to correspond 
       ! to S_star in routine spec_diag
       R_star=r_solar
       L_star=R_star*R_star*(4.0d0*pi*bb_total_flux)
       Edd_Efficiency=EddLeff_nominal
       pl_index=pl_index_nominal
    endif

#ifdef MPI       
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(T_eff,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         mympierror)
    call MPI_BCAST(R_star,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         mympierror)
    call MPI_BCAST(L_star,1,MPI_DOUBLE_PRECISION,0, & 
         MPI_COMM_NEW,mympierror)
    call MPI_BCAST(L_star_ion,1,MPI_DOUBLE_PRECISION,0, & 
         MPI_COMM_NEW,mympierror)
    call MPI_BCAST(S_star,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         mympierror)
#endif

  end subroutine spectrum_parms

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Determine spectrum diagnostics   
  subroutine spec_diag (powerlaw)

	  real(kind=dp), intent(in) :: powerlaw
    integer :: i_freq
    real(kind=dp) :: h_over_kT, freq_step 
    real(kind=dp) :: bb_photo_flux, pl_photo_flux
    real(kind=dp) :: S_star_unscaled, S_scaling
    real(kind=dp) :: S_star_band1, S_star_band2, S_star_band3
    real(kind=dp), dimension(0:NumFreq) :: frequency, weight
    real(kind=dp), dimension(0:NumFreq) :: bb_photon, bb_energy
    real(kind=dp), dimension(0:NumFreq) :: pl_photon, pl_energy

    real(kind=dp) :: freq_step_for_pl_scaling
    real(kind=dp), dimension(0:NumFreq_for_pl_scaling) :: frequency_for_pl_scaling
    real(kind=dp), dimension(0:NumFreq_for_pl_scaling) :: weight_for_pl_scaling	
    real(kind=dp), dimension(0:NumFreq_for_pl_scaling) :: pl_photon_for_pl_scaling
	real(kind=dp), dimension(0:NumFreq_for_pl_scaling) :: pl_energy_for_pl_scaling
		
    real(kind=dp), dimension(1:NumFreqBnd) :: S_star_subband
   
   
	L_0 = hplanck*ion_freq_HI*(powerlaw-1)*pl_input_flux
	write(*,*)'L_0',L_0
	
	
	
	
	
    ! This is h/kT
    h_over_kT = hplanck/(k_B*T_eff)

    ! Frequency step width
    freq_step = (freq_max(NumFreqBnd)-freq_min(1))/real(NumFreq)

    ! Fill the arrays for integration(frequency, weight, spectrum)
    do i_freq=0,NumFreq
      frequency(i_freq) = freq_min(1)+freq_step*real(i_freq)
      weight(i_freq) = freq_step
    enddo

    ! Frequency step width for pl scaling
    freq_step_for_pl_scaling = (freq_max(NumFreqBnd)-freq_min(1))/real(NumFreq_for_pl_scaling)

    ! Fill the arrays for integration(frequency, weight, spectrum) for pl scaling
    do i_freq=0,NumFreq_for_pl_scaling
      frequency_for_pl_scaling(i_freq) = freq_min(1)+freq_step_for_pl_scaling*real(i_freq)
      weight_for_pl_scaling(i_freq) = freq_step_for_pl_scaling
    enddo
	
    ! Case L_star_ion is provided. Find R_star, L_star and blackblody flux
    if (L_star_ion .ne. 0.0d0) then

      do i_freq=0,NumFreq
        if (frequency(i_freq)*h_over_kT .le. 709.0_dp) then
          ! this blackbody is in energy sense
          bb_energy(i_freq) = two_pi_over_c_square*frequency(i_freq)*frequency(i_freq)*frequency(i_freq)* &
                              hplanck/(exp(frequency(i_freq)*h_over_kT)-1.0)
          ! when the argument of the exponential function gets too high
        else
          bb_energy(i_freq) = two_pi_over_c_square*frequency(i_freq)*frequency(i_freq)*frequency(i_freq)*hplanck/ &
                              (exp((frequency(i_freq)*h_over_kT)/2.0_dp))/(exp((frequency(i_freq)*h_over_kT)/2.0_dp))
        endif
      enddo

        ! Black-body flux (photon sense)
        bb_photo_flux = scalar_romberg(bb_energy,weight,NumFreq,NumFreq,0) 
        ! Find radius from the scaled and specified ionizing luminosities
        R_star = dsqrt(L_star_ion/(4.0_dp*pi*bb_photo_flux))
        ! Find total luminosity from Stefan-Boltzmann law
        L_star = R_star*R_star*4.0_dp*pi*sigma_SB*T_eff*T_eff*T_eff*T_eff

    ! Case L_star_ion is not provided.  
    ! Now we know R_star and L_star (except in the case where S_star is provided).
    ! So we can continue to find blackbody flux 
    else
      do i_freq=0,NumFreq
        if (frequency(i_freq)*h_over_kT .le. 709.0_dp)then
          ! this blackbody is in number of photon sense
          bb_photon(i_freq) = two_pi_over_c_square*frequency(i_freq)*frequency(i_freq)/ &
                              (exp(frequency(i_freq)*h_over_kT)-1.0_dp)  
          ! when the argument of the exponential function gets too high
        else
          bb_photon(i_freq) = two_pi_over_c_square*frequency(i_freq)*frequency(i_freq)/ &
                              (exp((frequency(i_freq)*h_over_kT)/2.0_dp))/(exp((frequency(i_freq)*h_over_kT)/2.0_dp))
        endif     
      enddo
      ! Black-body flux (photon sense)
      bb_photo_flux = scalar_romberg(bb_photon,weight,NumFreq,NumFreq,0) 

    endif

    !Find out the power law scaling factor for the case Eddinton luminosity efficiency is provided.
    do i_freq=0,NumFreq_for_pl_scaling
       ! this power-law is in energy sense
       pl_energy_for_pl_scaling(i_freq) = &
	      frequency_for_pl_scaling(i_freq)**(-pl_index)*hplanck*frequency_for_pl_scaling(i_freq)
       ! this power-law is in number of photon sense
       pl_photon_for_pl_scaling(i_freq) = frequency_for_pl_scaling(i_freq)**(-pl_index)       
    enddo

    !Total power-law flux and power law scaling if Eddinton luminosity efficiency is provided.
    if (rank .eq. 0 .and. pl_input_flux .le. 0.001 ) then    
       pl_photo_flux = scalar_romberg_for_pl_scaling(pl_energy_for_pl_scaling,weight_for_pl_scaling,&
	         NumFreq_for_pl_scaling,NumFreq_for_pl_scaling,0) 
			 ! This power-law flux is not normalized (energy sense).
       pl_photo_flux_wanted = EddLum*Edd_Efficiency
	   !write(*,*) pl_photo_flux_wanted 
       pl_scaling = pl_photo_flux_wanted/pl_photo_flux
    endif 
 
    ! If power-law flux is provided, the power-law scaling is as follows.
    if (rank .eq. 0 .and. pl_input_flux .ge. 0.0 .and. T_eff_nominal .eq. 0.0) then
      pl_photo_flux = scalar_romberg_for_pl_scaling(pl_photon_for_pl_scaling,weight_for_pl_scaling,&
	  NumFreq_for_pl_scaling,NumFreq_for_pl_scaling,0) ! This power-law flux is not normalized (photon sense).
      pl_photo_flux_wanted = pl_input_flux
      pl_scaling = pl_photo_flux_wanted/pl_photo_flux
	  !write(*,*)'pl_scaling ',pl_scaling
    endif

    ! Find out S_star with the given radius.
    S_star_unscaled = 4.0*pi*R_star*R_star*bb_photo_flux 

    ! If S_star is zero, it is set here.
    if (S_star == 0.0) then
      S_star=S_star_unscaled

    ! If S_star is sepcified, then R_star and L_star have to be tuned accordingly.
    else
      S_scaling = S_star/S_star_unscaled
      R_star = sqrt(S_scaling)*R_star
      L_star = S_scaling*L_star
    endif

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Find out the number of photons in band 1
    freq_step=(freq_max(NumBndin1)-freq_min(NumBndin1))/real(NumFreq)

    do i_freq=0,NumFreq
      frequency(i_freq)=ion_freq_HI+freq_step*real(i_freq)
      weight(i_freq)=freq_step
      if (frequency(i_freq)*h_over_kT .le. 709.0_dp) then
        ! this blackbody is in number of photon sense
        bb_photon(i_freq) = two_pi_over_c_square*frequency(i_freq)*frequency(i_freq)/ &
                            (exp(frequency(i_freq)*h_over_kT)-1.0_dp)  
      ! when the argument of the exponential function gets too high
      else
        bb_photon(i_freq) = two_pi_over_c_square*frequency(i_freq)*frequency(i_freq)/ &
                            (exp((frequency(i_freq)*h_over_kT)/2.0_dp))/ &
                            (exp((frequency(i_freq)*h_over_kT)/2.0_dp))
      endif
      ! this power-law is in number of photon sense
      pl_photon(i_freq)=frequency(i_freq)**(-pl_index)
    enddo

    bb_photo_flux = scalar_romberg(bb_photon,weight,NumFreq,NumFreq,0)
    pl_photo_flux = scalar_romberg(pl_photon,weight,NumFreq,NumFreq,0)

    if (sourcetype .eq. 'B') S_star_band1 = 4.0*pi*R_star*R_star*bb_photo_flux
    if (sourcetype .eq. 'P') S_star_band1 = pl_photo_flux*pl_scaling
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Find out the number of photons in band 2
    freq_step = (freq_max(NumBndin1+NumBndin2)-freq_min(NumBndin1+1))/real(NumFreq)

    do i_freq=0,NumFreq
      frequency(i_freq) = ion_freq_HeI+freq_step*real(i_freq)
      weight(i_freq) = freq_step
      if (frequency(i_freq)*h_over_kT .le. 709.0_dp) then
        ! this blackbody is in number of photon sense
        bb_photon(i_freq) = two_pi_over_c_square*frequency(i_freq)*frequency(i_freq)/ &
                            (exp(frequency(i_freq)*h_over_kT)-1.0_dp)  
        ! when the argument of the exponential function gets too high
      else
        bb_photon(i_freq) = two_pi_over_c_square*frequency(i_freq)*frequency(i_freq)/ &
                            (exp((frequency(i_freq)*h_over_kT)/2.0_dp))/ &
                            (exp((frequency(i_freq)*h_over_kT)/2.0_dp))
      endif
      ! this power-law is in number of photon sense
      pl_photon(i_freq) = frequency(i_freq)**(-pl_index)
    enddo

    bb_photo_flux = scalar_romberg(bb_photon,weight,NumFreq,NumFreq,0)
    pl_photo_flux = scalar_romberg(pl_photon,weight,NumFreq,NumFreq,0)
    if (sourcetype .eq. 'B') S_star_band2 = 4.0*pi*R_star*R_star*bb_photo_flux
    if (sourcetype .eq. 'P') S_star_band2 = pl_photo_flux*pl_scaling
 
    ! Find out the number of photons in band 3
    freq_step = (freq_max(NumBndin1+NumBndin2+NumBndin3)-freq_min(NumBndin1+NumBndin2+1))/real(NumFreq)
    do i_freq=0,NumFreq
      frequency(i_freq) = ion_freq_HeII+freq_step*real(i_freq)
      weight(i_freq) = freq_step
      if (frequency(i_freq)*h_over_kT .le. 709.0_dp) then
        ! this blackbody is in number of photon sense
        bb_photon(i_freq) = two_pi_over_c_square*frequency(i_freq)*frequency(i_freq)/ &
                            (exp(frequency(i_freq)*h_over_kT)-1.0_dp)  
        ! when the argument of the exponential function gets too high
      else
        bb_photon(i_freq) = two_pi_over_c_square*frequency(i_freq)*frequency(i_freq)/ &
                            (exp((frequency(i_freq)*h_over_kT)/2.0_dp))/ &
                            (exp((frequency(i_freq)*h_over_kT)/2.0_dp))
      endif
       ! this power-law is in number of photon sense
       pl_photon(i_freq) = frequency(i_freq)**(-pl_index)
    enddo

    bb_photo_flux = scalar_romberg(bb_photon,weight,NumFreq,NumFreq,0)
    pl_photo_flux = scalar_romberg(pl_photon,weight,NumFreq,NumFreq,0)
    if (sourcetype .eq. 'B') S_star_band3 = 4.0*pi*R_star*R_star*bb_photo_flux 
    if (sourcetype .eq. 'P') S_star_band3 = pl_photo_flux*pl_scaling

    ! Report back to the log file
    if (rank == 0) then

      if(sourcetype.eq.'B') then
        write(logf,'(/a)')           'Using a black body with'
        write(logf,'(a,1pe10.3,a)')   ' Teff =       ', T_eff, ' K'
        write(logf,'(a,1pe10.3,a)')   ' Radius =     ', R_star/r_solar, ' R_solar'
        write(logf,'(a,1pe10.3,a)')   ' Luminosity = ', L_star/l_solar, ' L_solar'
        write(logf,*) 'ionzing photon production rate = ', S_star, ' s^-1'
      endif

      if(sourcetype.eq.'P') then
        write(logf,'(/a)')           'Using a power law source with'
        write(logf,'(a,1pe10.3)')   ' power law index = ', pl_index
        write(logf,'(a,1pe10.3)')   ' efficiency parameter = ', Edd_Efficiency
      endif

      write(logf,'(A,(1PE12.5),A//)') ' Number of photons in band 1: ', S_star_band1, ' s^-1'
      write(logf,'(A,(1PE12.5),A//)') ' Number of photons in band 2: ', S_star_band2, ' s^-1'
      write(logf,'(A,(1PE12.5),A//)') ' Number of photons in band 3: ', S_star_band3, ' s^-1'
      write(logf,'(A,(1PE12.5),A//)') ' Total number of ionizing photons: ', &
                                               S_star_band1+S_star_band2+S_star_band3, 's^-1'            
    endif
 
	write(*,*) 'photon rate (band1+band2+band3)',S_star_band1+S_star_band2+S_star_band3
	
  end subroutine spec_diag
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Generate photoionization tables and heating tables   
  subroutine spec_integration ()

    integer :: i_freq, i_tau, i_subband
    real(kind=dp) :: R_star2, h_over_kT
    real(kind=dp), dimension(0:NumTau) :: tau
    real(kind=dp), dimension(0:NumTau) :: answer
    real(kind=dp), dimension(0:NumFreq) :: exponent_HI
    real(kind=dp), dimension(0:NumFreq) :: exponent_HeI
    real(kind=dp), dimension(0:NumFreq) :: exponent_HeII
    real(kind=dp), dimension(0:NumFreq) :: frequency
    real(kind=dp), dimension(0:NumFreq, 0:NumTau) :: weight

    ! Photoionization integrand as a function of frequency and tau
    allocate(bb_photo_thick_integrand(0:NumFreq, 0:NumTau))    
    allocate(bb_photo_thin_integrand(0:NumFreq, 0:NumTau)) 
    allocate(pl_photo_thick_integrand(0:NumFreq, 0:NumTau))    
    allocate(pl_photo_thin_integrand(0:NumFreq, 0:NumTau)) 

    ! Heating integrand as a function of frequency and tau
    if (.not.isothermal) then
    	allocate(bb_heat_thick_integrand_HI(0:NumFreq, 0:NumTau))  
    	allocate(bb_heat_thick_integrand_HeI(0:NumFreq, 0:NumTau))  
    	allocate(bb_heat_thick_integrand_HeII(0:NumFreq, 0:NumTau))      
    	allocate(bb_heat_thin_integrand_HI(0:NumFreq, 0:NumTau))   
    	allocate(bb_heat_thin_integrand_HeI(0:NumFreq, 0:NumTau))   
    	allocate(bb_heat_thin_integrand_HeII(0:NumFreq, 0:NumTau))   
    	allocate(pl_heat_thick_integrand_HI(0:NumFreq, 0:NumTau))   
    	allocate(pl_heat_thick_integrand_HeI(0:NumFreq, 0:NumTau))   
    	allocate(pl_heat_thick_integrand_HeII(0:NumFreq, 0:NumTau))      
    	allocate(pl_heat_thin_integrand_HI(0:NumFreq, 0:NumTau))   
    	allocate(pl_heat_thin_integrand_HeI(0:NumFreq, 0:NumTau))   
    	allocate(pl_heat_thin_integrand_HeII(0:NumFreq, 0:NumTau))   
    endif

    ! Photoionization table as a function of photo sub-bin and tau
    allocate(bb_photo_thick_table(0:NumTau, 1:NumFreqBnd))
    allocate(bb_photo_thin_table(0:NumTau, 1:NumFreqBnd))
    allocate(pl_photo_thick_table(0:NumTau, 1:NumFreqBnd))
    allocate(pl_photo_thin_table(0:NumTau, 1:NumFreqBnd))

    ! Heating table as a function of heating sub-bin and tau
    if (.not.isothermal) then
	allocate(bb_heat_thick_table(0:NumTau, 1:NumheatBin))
    	allocate(bb_heat_thin_table(0:NumTau, 1:NumheatBin))
	allocate(pl_heat_thick_table(0:NumTau, 1:NumheatBin))
    	allocate(pl_heat_thin_table(0:NumTau, 1:NumheatBin))
    endif

    ! This is h/kT
    h_over_kT=hplanck/(k_B*T_eff)

    ! This is R_star^2
    R_star2=R_star*R_star

    ! fill the optical depth array used to fill the tables 
    ! it is filled in NumTau logarithmic steps 
    ! from minlogtau to maxlogtau
    do i_tau = 1,NumTau
       tau(i_tau) = 10.0**(minlogtau+dlogtau*real(i_tau-1))
    enddo

    ! Position zero corresponds to zero optical depth
    tau(0)=0.0

    ! Warn about grey opacities:
    if (grey .and. rank == 0) write(logf,*) 'WARNING: Using grey opacities'

    ! In frequency band 1, fill in integrands and make tables

    ! Go through all the sub-bin in band 1
    do i_subband=1,NumBndin1  
     
      ! Assign values to exponent array
      do i_freq=0,NumFreq

        frequency(i_freq) = freq_min(i_subband)+delta_freq(i_subband)*real(i_freq)
        if (grey) then    
          exponent_HI(i_freq) = 1.0
        else
          exponent_HI(i_freq) = ((frequency(i_freq)/freq_min(i_subband))**(-pl_index_cross_section_HI(i_subband)))        
        endif

      enddo

      ! Loop through the tau partition
      do i_tau=0,NumTau 

        ! Loop through the frequency partition
        do i_freq=0,NumFreq
          weight(i_freq,i_tau) = delta_freq(i_subband)  

          ! Assign values to the photo integrands
          if (tau(i_tau)*exponent_HI(i_freq) < 700.0) then   
            bb_photo_thick_integrand(i_freq,i_tau) = 4.0_dp*pi*R_star2*two_pi_over_c_square*frequency(i_freq)* &
                                                     frequency(i_freq)*exp(-tau(i_tau)*exponent_HI(i_freq))/ &
                                                     (exp(frequency(i_freq)*h_over_kT)-1.0)
            bb_photo_thin_integrand(i_freq,i_tau) = 4.0_dp*pi*R_star2*two_pi_over_c_square*frequency(i_freq)* &
                                                    frequency(i_freq)*exponent_HI(i_freq)*exp(-tau(i_tau)* &
                                                    exponent_HI(i_freq))/(exp(frequency(i_freq)*h_over_kT)-1.0)
            pl_photo_thick_integrand(i_freq,i_tau) = pl_scaling*frequency(i_freq)**(-pl_index)* &
                                                     exp(-tau(i_tau)*exponent_HI(i_freq))
            pl_photo_thin_integrand(i_freq,i_tau) = pl_scaling*frequency(i_freq)**(-pl_index)*exponent_HI(i_freq) &
                                                      *exp(-tau(i_tau)*exponent_HI(i_freq))
          else
            bb_photo_thick_integrand(i_freq,i_tau) = 0.0
            bb_photo_thin_integrand(i_freq,i_tau) = 0.0
            pl_photo_thick_integrand(i_freq,i_tau) = 0.0
            pl_photo_thin_integrand(i_freq,i_tau) = 0.0
          endif

          ! Assign values to the heating integrands
          if (.not.isothermal) then
            bb_heat_thick_integrand_HI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HI)* &
                                                       bb_photo_thick_integrand(i_freq,i_tau)
            bb_heat_thin_integrand_HI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HI)* &
                                                      bb_photo_thin_integrand(i_freq,i_tau)
            pl_heat_thick_integrand_HI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HI)* &
                                                       pl_photo_thick_integrand(i_freq,i_tau)
            pl_heat_thin_integrand_HI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HI)* &
                                                      pl_photo_thin_integrand(i_freq,i_tau)
          endif

        enddo

      enddo

      ! Make photo tables
      call vector_romberg (bb_photo_thick_integrand,weight,NumFreq,NumFreq,NumTau,answer)
      bb_photo_thick_table(:,1) = answer
      call vector_romberg (bb_photo_thin_integrand,weight,NumFreq,NumFreq,NumTau,answer)
      bb_photo_thin_table(:,1) = answer
      call vector_romberg (pl_photo_thick_integrand,weight,NumFreq,NumFreq,NumTau,answer)
      pl_photo_thick_table(:,1) = answer
      call vector_romberg (pl_photo_thin_integrand,weight,NumFreq,NumFreq,NumTau,answer)
      pl_photo_thin_table(:,1) = answer

      ! Make heating tables
      if (.not.isothermal) then
        call vector_romberg (bb_heat_thick_integrand_HI,weight,NumFreq,NumFreq,NumTau,answer)
        bb_heat_thick_table(:,1) = answer
        call vector_romberg (bb_heat_thin_integrand_HI,weight,NumFreq,NumFreq,NumTau,answer)
        bb_heat_thin_table(:,1) = answer
        call vector_romberg (pl_heat_thick_integrand_HI,weight,NumFreq,NumFreq,NumTau,answer)
        pl_heat_thick_table(:,1) = answer
        call vector_romberg (pl_heat_thin_integrand_HI,weight,NumFreq,NumFreq,NumTau,answer)
        pl_heat_thin_table(:,1) = answer
      endif

    enddo

    ! In frequency band 2, fill in integrands and make tables

    ! Go through all the sub-bin in band 2
    do i_subband=NumBndin1+1,NumBndin1+NumBndin2
       
      ! Assign values to exponent array
      do i_freq=0,NumFreq

        frequency(i_freq) = freq_min(i_subband)+delta_freq(i_subband)*real(i_freq)
        if (grey) then
	  exponent_HeI(i_freq) = 1.0
        else
          exponent_HeI(i_freq) = ((frequency(i_freq)/freq_min(i_subband))**(-pl_index_cross_section_HeI(i_subband)))
        endif

      enddo
         
      ! Loop through the tau partition
      do i_tau=0,NumTau 

        ! Loop through the frequency partition
        do i_freq=0,NumFreq
          weight(i_freq,i_tau) = delta_freq(i_subband)  
              
          ! Assign values to the photo integrands
          if (tau(i_tau)*exponent_HeI(i_freq) < 700.0) then 
            bb_photo_thick_integrand(i_freq,i_tau) = 4.0_dp*pi*R_star2*two_pi_over_c_square*frequency(i_freq)* &
                                                     frequency(i_freq)*exp(-(tau(i_tau)*exponent_HeI(i_freq)))/ &
                                                     (exp(frequency(i_freq)*h_over_kT)-1.0)
            bb_photo_thin_integrand(i_freq,i_tau) = 4.0_dp*pi*R_star2*two_pi_over_c_square*frequency(i_freq)* &
                                                    frequency(i_freq)* exponent_HeI(i_freq)*exp(-(tau(i_tau)* &
                                                    exponent_HeI(i_freq)))/(exp(frequency(i_freq)*h_over_kT)-1.0)
            pl_photo_thick_integrand(i_freq,i_tau) = pl_scaling*frequency(i_freq)**(-pl_index)* &
                                                     exp(-tau(i_tau)*exponent_HeI(i_freq))
            pl_photo_thin_integrand(i_freq,i_tau) = pl_scaling*frequency(i_freq)**(-pl_index)*exponent_HeI(i_freq)* &
                                                    exp(-tau(i_tau)*exponent_HeI(i_freq))
          else
            bb_photo_thick_integrand(i_freq,i_tau) = 0.0  
            bb_photo_thin_integrand(i_freq,i_tau) = 0.0   
            pl_photo_thick_integrand(i_freq,i_tau) = 0.0  
            pl_photo_thin_integrand(i_freq,i_tau) = 0.0  
          endif            

          ! Assign values to the heating integrands
          if (.not.isothermal) then
            bb_heat_thick_integrand_HI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HI)* &
                                                       bb_photo_thick_integrand(i_freq,i_tau)
            bb_heat_thick_integrand_HeI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HeI)* &
                                                        bb_photo_thick_integrand(i_freq,i_tau)
            bb_heat_thin_integrand_HI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HI)*  &
                                                      bb_photo_thin_integrand(i_freq,i_tau)
            bb_heat_thin_integrand_HeI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HeI)* &
                                                       bb_photo_thin_integrand(i_freq,i_tau)
            pl_heat_thick_integrand_HI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HI)* &
                                                       pl_photo_thick_integrand(i_freq,i_tau)
            pl_heat_thick_integrand_HeI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HeI)* &
                                                        pl_photo_thick_integrand(i_freq,i_tau)
            pl_heat_thin_integrand_HI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HI)* &
                                                      pl_photo_thin_integrand(i_freq,i_tau)
            pl_heat_thin_integrand_HeI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HeI)* &
                                                       pl_photo_thin_integrand(i_freq,i_tau)
          endif 

        enddo 
        
      enddo   

      ! Make photo tables
      call vector_romberg (bb_photo_thick_integrand,weight,NumFreq,NumFreq,NumTau,answer)
      bb_photo_thick_table(:,i_subband) = answer
      call vector_romberg (bb_photo_thin_integrand,weight,NumFreq,NumFreq,NumTau,answer)
      bb_photo_thin_table(:,i_subband) = answer
      call vector_romberg (pl_photo_thick_integrand,weight,NumFreq,NumFreq,NumTau,answer)
      pl_photo_thick_table(:,i_subband) = answer
      call vector_romberg (pl_photo_thin_integrand,weight,NumFreq,NumFreq,NumTau,answer)
      pl_photo_thin_table(:,i_subband) = answer 
 
      ! Make heating tables
      if (.not.isothermal) then
        call vector_romberg (bb_heat_thick_integrand_HI,weight,NumFreq,NumFreq,NumTau,answer)
        bb_heat_thick_table(:,i_subband*2-2) = answer
        call vector_romberg (bb_heat_thick_integrand_HeI,weight,NumFreq,NumFreq,NumTau,answer)
        bb_heat_thick_table(:,i_subband*2-1) = answer
        call vector_romberg (bb_heat_thin_integrand_HI,weight,NumFreq,NumFreq,NumTau,answer)
        bb_heat_thin_table(:,i_subband*2-2) = answer
        call vector_romberg (bb_heat_thin_integrand_HeI,weight,NumFreq,NumFreq,NumTau,answer)
        bb_heat_thin_table(:,i_subband*2-1) = answer
        call vector_romberg (pl_heat_thick_integrand_HI,weight,NumFreq,NumFreq,NumTau,answer)
        pl_heat_thick_table(:,i_subband*2-2) = answer
        call vector_romberg (pl_heat_thick_integrand_HeI,weight,NumFreq,NumFreq,NumTau,answer)
        pl_heat_thick_table(:,i_subband*2-1) = answer
        call vector_romberg (pl_heat_thin_integrand_HI,weight,NumFreq,NumFreq,NumTau,answer)
        pl_heat_thin_table(:,i_subband*2-2) = answer
        call vector_romberg (pl_heat_thin_integrand_HeI,weight,NumFreq,NumFreq,NumTau,answer)
        pl_heat_thin_table(:,i_subband*2-1) = answer
      endif

    enddo            

    ! In frequency band 3, fill in integrands and make tables

    ! Go through all the sub-bin in band 3
    do i_subband=NumBndin1+NumBndin2+1,NumBndin1+NumBndin2+NumBndin3

      ! Assign values to exponent array
      do i_freq=0,NumFreq

        frequency(i_freq) = freq_min(i_subband)+delta_freq(i_subband)*real(i_freq)
        if (grey) then
	  exponent_HeII(i_freq) = 1.0
        else
          exponent_HeII(i_freq) = ((frequency(i_freq)/freq_min(i_subband))**(-pl_index_cross_section_HeII(i_subband)))
        endif 

      enddo

      ! Loop through the tau partition
      do i_tau=0,NumTau 
      
        ! Loop through the frequency partition 
        do i_freq=0,NumFreq
           weight(i_freq,i_tau) = delta_freq(i_subband)  

          ! Assign values to the photo integrands
          if (tau(i_tau)*exponent_HeII(i_freq) < 700.0) then  
            bb_photo_thick_integrand(i_freq,i_tau) = 4.0_dp*pi*R_star2*two_pi_over_c_square*frequency(i_freq)* &
                                                     frequency(i_freq)*exp(-tau(i_tau)*exponent_HeII(i_freq))/ &
                                                     (exp(frequency(i_freq)*h_over_kT)-1.0)  
            bb_photo_thin_integrand(i_freq,i_tau) = 4.0_dp*pi*R_star2*two_pi_over_c_square*frequency(i_freq)* &
                                                    frequency(i_freq)*exponent_HeII(i_freq)*exp(-tau(i_tau)* &
                                                    exponent_HeII(i_freq))/(exp(frequency(i_freq)*h_over_kT)-1.0)  
            pl_photo_thick_integrand(i_freq,i_tau) = pl_scaling*frequency(i_freq)**(-pl_index)* &
                                                     exp(-tau(i_tau)*exponent_HeII(i_freq))
            pl_photo_thin_integrand(i_freq,i_tau) = pl_scaling*frequency(i_freq)**(-pl_index)*exponent_HeII(i_freq)* &
                                                    exp(-tau(i_tau)*exponent_HeII(i_freq))
	  else
            bb_photo_thick_integrand(i_freq,i_tau) = 0.0
            bb_photo_thin_integrand(i_freq,i_tau) = 0.0
            pl_photo_thick_integrand(i_freq,i_tau) = 0.0
            pl_photo_thin_integrand(i_freq,i_tau) = 0.0
	  endif 

          ! Assign values to the heating integrands
          if (.not.isothermal) then
            bb_heat_thick_integrand_HI(i_freq,i_tau) = 4.0_dp*pi*hplanck*(frequency(i_freq)-ion_freq_HI)* &
                                                       bb_photo_thick_integrand(i_freq,i_tau)
            bb_heat_thick_integrand_HeI(i_freq,i_tau) = 4.0_dp*pi*hplanck*(frequency(i_freq)-ion_freq_HeI)* &
                                                        bb_photo_thick_integrand(i_freq,i_tau)
            bb_heat_thick_integrand_HeII(i_freq,i_tau) = 4.0_dp*pi*hplanck*(frequency(i_freq)-ion_freq_HeII)* &
                                                         bb_photo_thick_integrand(i_freq,i_tau)
            bb_heat_thin_integrand_HI(i_freq,i_tau) = 4.0_dp*pi*hplanck*(frequency(i_freq)-ion_freq_HI)* &
                                                      bb_photo_thin_integrand(i_freq,i_tau)
            bb_heat_thin_integrand_HeI(i_freq,i_tau) = 4.0_dp*pi*hplanck*(frequency(i_freq)-ion_freq_HeI)* &
                                                       bb_photo_thin_integrand(i_freq,i_tau)
            bb_heat_thin_integrand_HeII(i_freq,i_tau) = 4.0_dp*pi*hplanck*(frequency(i_freq)-ion_freq_HeII)* &
                                                        bb_photo_thin_integrand(i_freq,i_tau)
            pl_heat_thick_integrand_HI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HI)*  &
                                                       pl_photo_thick_integrand(i_freq,i_tau)
            pl_heat_thick_integrand_HeI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HeI)* &
                                                        pl_photo_thick_integrand(i_freq,i_tau)
            pl_heat_thick_integrand_HeII(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HeII)* &
                                                         pl_photo_thick_integrand(i_freq,i_tau)
            pl_heat_thin_integrand_HI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HI)* &
                                                      pl_photo_thin_integrand(i_freq,i_tau)
            pl_heat_thin_integrand_HeI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HeI)* &
                                                       pl_photo_thin_integrand(i_freq,i_tau)
            pl_heat_thin_integrand_HeII(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HeII)* &
                                                        pl_photo_thin_integrand(i_freq,i_tau)
          endif 

        enddo

      enddo  

      ! Make photo tables
      call vector_romberg (bb_photo_thick_integrand,weight,NumFreq,NumFreq,NumTau,answer)
      bb_photo_thick_table(:,i_subband) = answer
      call vector_romberg (bb_photo_thin_integrand,weight,NumFreq,NumFreq,NumTau,answer)
      bb_photo_thin_table(:,i_subband) = answer
      call vector_romberg (pl_photo_thick_integrand,weight,NumFreq,NumFreq,NumTau,answer)
      pl_photo_thick_table(:,i_subband) = answer
      call vector_romberg (pl_photo_thin_integrand,weight,NumFreq,NumFreq,NumTau,answer)
      pl_photo_thin_table(:,i_subband) = answer  

      ! Make heating tables
      if (.not.isothermal) then
        call vector_romberg (bb_heat_thick_integrand_HI,weight,NumFreq,NumFreq,NumTau,answer)
        bb_heat_thick_table(:,i_subband*3-NumBndin2-4) = answer
        call vector_romberg (bb_heat_thick_integrand_HeI,weight,NumFreq,NumFreq,NumTau,answer)
        bb_heat_thick_table(:,i_subband*3-NumBndin2-3) = answer
        call vector_romberg (bb_heat_thick_integrand_HeII,weight,NumFreq,NumFreq,NumTau,answer)
        bb_heat_thick_table(:,i_subband*3-NumBndin2-2) = answer
        call vector_romberg (bb_heat_thin_integrand_HI,weight,NumFreq,NumFreq,NumTau,answer)
        bb_heat_thin_table(:,i_subband*3-NumBndin2-4) = answer
        call vector_romberg (bb_heat_thin_integrand_HeI,weight,NumFreq,NumFreq,NumTau,answer)
        bb_heat_thin_table(:,i_subband*3-NumBndin2-3) = answer
        call vector_romberg (bb_heat_thin_integrand_HeII,weight,NumFreq,NumFreq,NumTau,answer)
        bb_heat_thin_table(:,i_subband*3-NumBndin2-2) = answer
        call vector_romberg (pl_heat_thick_integrand_HI,weight,NumFreq,NumFreq,NumTau,answer)
        pl_heat_thick_table(:,i_subband*3-NumBndin2-4) = answer
        call vector_romberg (pl_heat_thick_integrand_HeI,weight,NumFreq,NumFreq,NumTau,answer)
        pl_heat_thick_table(:,i_subband*3-NumBndin2-3) = answer
        call vector_romberg (pl_heat_thick_integrand_HeII,weight,NumFreq,NumFreq,NumTau,answer)
        pl_heat_thick_table(:,i_subband*3-NumBndin2-2) = answer
        call vector_romberg (pl_heat_thin_integrand_HI,weight,NumFreq,NumFreq,NumTau,answer)
        pl_heat_thin_table(:,i_subband*3-NumBndin2-4) = answer
        call vector_romberg (pl_heat_thin_integrand_HeI,weight,NumFreq,NumFreq,NumTau,answer)
        pl_heat_thin_table(:,i_subband*3-NumBndin2-3) = answer
        call vector_romberg (pl_heat_thin_integrand_HeII,weight,NumFreq,NumFreq,NumTau,answer)
        pl_heat_thin_table(:,i_subband*3-NumBndin2-2) = answer
      endif

    enddo           

    ! deallocate the useless photo integrand
    deallocate(bb_photo_thick_integrand)
    deallocate(bb_photo_thin_integrand)
    deallocate(pl_photo_thick_integrand)
    deallocate(pl_photo_thin_integrand)

    ! deallocate the useless heating integrand
    if (.not.isothermal) then
      deallocate(bb_heat_thick_integrand_HI)
      deallocate(bb_heat_thick_integrand_HeI)
      deallocate(bb_heat_thick_integrand_HeII)
      deallocate(bb_heat_thin_integrand_HI)
      deallocate(bb_heat_thin_integrand_HeI)
      deallocate(bb_heat_thin_integrand_HeII)
      deallocate(pl_heat_thick_integrand_HI)
      deallocate(pl_heat_thick_integrand_HeI)
      deallocate(pl_heat_thick_integrand_HeII)
      deallocate(pl_heat_thin_integrand_HI)
      deallocate(pl_heat_thin_integrand_HeI)
      deallocate(pl_heat_thin_integrand_HeII)
    endif

  end subroutine spec_integration

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! this subroutine calculates photo-ionization rates at a particular sets of column density
  subroutine photoion (phi,colum_in_HI,colum_out_HI,colum_in_HeI,colum_out_HeI, &
	       	       colum_in_HeII,colum_out_HeII,vol,nsrc,i_state)

    use sourceprops, only: NormFlux,srcMass,NormFluxPL
    use cgsphotoconstants

    ! result of the routine
    type(photrates), intent(out) :: phi

    ! Incoming and outgoing HI column density
    real(kind=dp), intent(in) :: colum_in_HI, colum_out_HI

    ! Incoming and outgoing HeI column density
    real(kind=dp), intent(in) :: colum_in_HeI, colum_out_HeI

    ! Incoming and outgoing HeII column density
    real(kind=dp), intent(in) :: colum_in_HeII, colum_out_HeII

    ! Volume of shell cell
    real(kind=dp), intent(in) :: vol

    real(kind=dp), intent(in) :: i_state

    ! Number of the source
    integer, intent(in) :: nsrc 

    integer :: i_tau, i_subband
    real(kind=dp) :: colum_cell_HI
    real(kind=dp) :: colum_cell_HeI
    real(kind=dp) :: colum_cell_HeII
    real(kind=dp) :: NFlux
    real(kind=dp), dimension(1:NumFreqBnd) :: tau_in_all
    real(kind=dp), dimension(1:NumFreqBnd) :: tau_out_all
    real(kind=dp), dimension(1:NumFreqBnd) :: tau_cell_HI
    real(kind=dp), dimension(1:NumFreqBnd) :: tau_cell_HeI
    real(kind=dp), dimension(1:NumFreqBnd) :: tau_cell_HeII
    type(tablepos) :: tau_pos_in, tau_pos_out

    ! The column densities (HI, HeI, HeII) at current cell
    colum_cell_HI = colum_out_HI-colum_in_HI
    colum_cell_HeI = colum_out_HeI-colum_in_HeI
    colum_cell_HeII = colum_out_HeII-colum_in_HeII 

    ! The optical depths (HI, HeI, HeII) at current cell
    do i_subband=1,NumFreqBnd
      tau_cell_HI(i_subband) = colum_cell_HI*sigma_HI(i_subband)
      tau_cell_HeI(i_subband) = colum_cell_HeI*sigma_HeI(i_subband)
      tau_cell_HeII(i_subband) = colum_cell_HeII*sigma_HeII(i_subband)
    enddo      
              
    ! total tau_in (including HI, HeI, HeII)
    do i_subband=1,NumFreqBnd
      tau_in_all(i_subband) = colum_in_HI*sigma_HI(i_subband)+ &
                              colum_in_HeI*sigma_HeI(i_subband)+ &
                              colum_in_HeII*sigma_HeII(i_subband)
    enddo  

    ! total tau_out (including HI, HeI, HeII)
    do i_subband=1,NumFreqBnd
      tau_out_all(i_subband) = colum_out_HI*sigma_HI(i_subband)+ &
                               colum_out_HeI*sigma_HeI(i_subband)+ &
                               colum_out_HeII*sigma_HeII(i_subband)
    enddo

    ! find the table positions for the optical depth (ingoing)
    do i_subband=1,NumFreqBnd  
      tau_pos_in%tau(i_subband) = log10(max(1.0e-20_dp,tau_in_all(i_subband)))
      tau_pos_in%odpos(i_subband) = min(real(NumTau,dp),max(0.0_dp,1.0+ &
                                      (tau_pos_in%tau(i_subband)-minlogtau)/dlogtau))
      tau_pos_in%ipos(i_subband) = int(tau_pos_in%odpos(i_subband))
      tau_pos_in%residual(i_subband) = tau_pos_in%odpos(i_subband)-real(tau_pos_in%ipos(i_subband),dp)
      tau_pos_in%ipos_p1(i_subband) = min(NumTau,tau_pos_in%ipos(i_subband)+1)
    enddo

    ! find the table positions for the optical depth (outgoing)
    do i_subband=1,NumFreqBnd  
      tau_pos_out%tau(i_subband) = log10(max(1.0e-20_dp,tau_out_all(i_subband)))
      tau_pos_out%odpos(i_subband) = min(real(NumTau,dp),max(0.0_dp,1.0+ &
                                     (tau_pos_out%tau(i_subband)-minlogtau)/dlogtau))
      tau_pos_out%ipos(i_subband) = int(tau_pos_out%odpos(i_subband))
      tau_pos_out%residual(i_subband) = tau_pos_out%odpos(i_subband)-real(tau_pos_out%ipos(i_subband),dp)
      tau_pos_out%ipos_p1(i_subband) = min(NumTau,tau_pos_out%ipos(i_subband)+1)
    enddo 

    if (NormFluxPL(nsrc) .le. 1.0e-100_dp) then  
      NFlux=NormFlux(nsrc) 
      sourcetype='B'
    else
      NFlux=NormFluxPL(nsrc)
      sourcetype='P'
    endif

    call lookuptable(tau_pos_in,tau_pos_out,phi,tau_in_all,tau_out_all, &
                     tau_cell_HI,tau_cell_HeI,tau_cell_HeII,NFlux,sourcetype, &
                     vol,i_state,colum_cell_HI,colum_cell_HeI, colum_cell_HeII )

  end subroutine photoion
 
  ! find out the correct position in the photo and heating tables
  subroutine lookuptable(tau_pos_in,tau_pos_out,phi,tau_in_all,tau_out_all, &
                         tau_cell_HI,tau_cell_HeI,tau_cell_HeII,NFlux,sourcetype, &
                         vol,i_state,colum_cell_HI,colum_cell_HeI,colum_cell_HeII)

    use cgsphotoconstants
    use material, only: gamma_uvb
	
    type(photrates), intent(out) :: phi
    type(tablepos), intent(in) :: tau_pos_in, tau_pos_out
    real(kind=dp), intent(in) :: NFlux, vol, i_state
    real(kind=dp), intent(in) :: colum_cell_HI, colum_cell_HeI, colum_cell_HeII
    real(kind=dp), dimension(NumFreqBnd), intent(in) :: tau_in_all, tau_out_all
    character,intent(in) :: sourcetype
    real(kind=dp), dimension(NumFreqBnd),intent(in) :: tau_cell_HI, tau_cell_HeI, tau_cell_HeII

    integer ::  n, i_subband, i
    real(kind=dp) :: phi_heat_HI, phi_heat_HeI, phi_heat_HeII
    real(kind=dp) :: f_heat, f_ion_HI, f_ion_HeI
    real(kind=dp) :: phi_photo_in_all, phi_photo_out_all, phi_photo_all
    real(kind=dp) :: fra_sum1, fra_sum2, fra_sum3, fra_sum4
    real(kind=dp), pointer, dimension(:,:) :: photo_thick_table, photo_thin_table
    real(kind=dp), pointer, dimension(:,:) :: heat_thick_table, heat_thin_table
    real(kind=dp), dimension(NumBndin1+1:NumBndin1+NumBndin2+NumBndin3) :: scaling_HI
    real(kind=dp), dimension(NumBndin1+1:NumBndin1+NumBndin2+NumBndin3) :: scaling_HeI
    real(kind=dp), dimension(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+NumBndin3) :: scaling_HeII
    real(kind=dp), dimension(1:NumheatBin) :: scaling
    real(kind=dp) :: test1, test2
    real(kind=dp) :: tau_photo_limit = 1.0e-7 
    real(kind=dp) :: tau_heat_limit = 1.0e-4

    ! pointers point to some variables
    if (sourcetype.eq.'B') then 
      photo_thick_table => bb_photo_thick_table
      photo_thin_table => bb_photo_thin_table
      heat_thick_table => bb_heat_thick_table
      heat_thin_table => bb_heat_thin_table
    elseif (sourcetype.eq.'P') then
      photo_thick_table => pl_photo_thick_table
      photo_thin_table => pl_photo_thin_table
      heat_thick_table => pl_heat_thick_table
      heat_thin_table => pl_heat_thin_table
    endif

    ! initialization
    phi%photo_cell_HI = 0.0_dp
    phi%photo_cell_HeI = 0.0_dp
    phi%photo_cell_HeII = 0.0_dp
    phi_heat_HI = 0.0_dp
    phi_heat_HeI = 0.0_dp
    phi_heat_HeII = 0.0_dp
    phi%heat_cell_HI = 0.0_dp
    phi%heat_cell_HeI = 0.0_dp
    phi%heat_cell_HeII = 0.0_dp
    phi%photo_in = 0.0_dp
    f_heat = 0.0_dp
    f_ion_HI = 0.0_dp
    f_ion_HeI = 0.0_dp

    ! loop through all frequency band
	
	
	
    do i_subband=1,NumFreqBnd
    !do i_subband=NumBndin1+6,NumFreqBnd
		
		
		
      ! Incoming, outcoming, current cell total photoionization rate
      phi_photo_in_all = (photo_thick_table(tau_pos_in%ipos(i_subband),i_subband)+ &
                         (photo_thick_table(tau_pos_in%ipos_p1(i_subband),i_subband)- &
                         photo_thick_table(tau_pos_in%ipos(i_subband),i_subband))* &
                         tau_pos_in%residual(i_subband))*NFlux
      phi%photo_in = phi%photo_in+phi_photo_in_all

      ! When current cell is optically thick
      if (abs(tau_out_all(i_subband)-tau_in_all(i_subband)) .gt. tau_photo_limit) then
        phi_photo_out_all = (photo_thick_table(tau_pos_out%ipos(i_subband),i_subband)+ &
                            (photo_thick_table(tau_pos_out%ipos_p1(i_subband),i_subband)- &
                            photo_thick_table(tau_pos_out%ipos(i_subband),i_subband))* &
                            tau_pos_out%residual(i_subband))*NFlux 
        phi_photo_all = phi_photo_in_all-phi_photo_out_all

      ! When current cell is optically thin
      else
        phi_photo_all = ((photo_thin_table(tau_pos_in%ipos(i_subband),i_subband)+ &
                        (photo_thin_table(tau_pos_in%ipos_p1(i_subband),i_subband)- &
                        photo_thin_table(tau_pos_in%ipos_p1(i_subband),i_subband))* &
                        tau_pos_in%residual(i_subband))* &
                        (tau_out_all(i_subband)-tau_in_all(i_subband)))*NFlux
        phi_photo_out_all = phi_photo_in_all-phi_photo_all
      endif

      ! Current cell individual photoionization rate of HI, HeI, HeII
      select case (i_subband) 

      ! band 1
      case (NumBndin1) 
    
        phi%photo_cell_HI = phi_photo_all/vol

      ! band 2
      case (NumBndin1+1:NumBndin1+NumBndin2)
        
        call scale_int2(scaling_HI(i_subband),scaling_HeI(i_subband),colum_cell_HI,colum_cell_HeI, i_subband)

        phi%photo_cell_HI = phi%photo_cell_HI+scaling_HI(i_subband)*phi_photo_all/vol 
        phi%photo_cell_HeI = phi%photo_cell_HeI+scaling_HeI(i_subband)*phi_photo_all/vol

      ! band 3
      case (NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+NumBndin3)

        call scale_int3(scaling_HI(i_subband),scaling_HeI(i_subband),scaling_HeII(i_subband), &
                        colum_cell_HI,colum_cell_HeI,colum_cell_HeII,i_subband)
        phi%photo_cell_HI = phi%photo_cell_HI+scaling_HI(i_subband)*phi_photo_all/vol
        phi%photo_cell_HeI = phi%photo_cell_HeI+scaling_HeI(i_subband)*phi_photo_all/vol
        phi%photo_cell_HeII = phi%photo_cell_HeII+scaling_HeII(i_subband)*phi_photo_all/vol

      end select
     
    enddo

    if(.not.isothermal) then

      ! in general, I'm following Ricotti et al 2002
      CR1 = (/0.3908_dp, 0.0554_dp, 1.0_dp/)
      bR1 = (/0.4092_dp, 0.4614_dp, 0.2663_dp/)
      dR1 = (/1.7592_dp, 1.6660_dp, 1.3163_dp/)
      CR2 = (/0.6941_dp,0.0984_dp,3.9811_dp/)
      aR2 = (/0.2_dp,0.2_dp,0.4_dp/)
      bR2 = (/0.38_dp,0.38_dp,0.34_dp/)
      test1 = 0.0_dp
      test2 = 0.0_dp

      do i=1,3
        y1R(i) = CR1(i)*(1.0_dp-i_state**bR1(i))**dR1(i)
        xeb = 1.0_dp-i_state**bR2(i) 
        y2R(i) = CR2(i)*i_state**aR2(i)*xeb*xeb
      enddo

      ! Current cell individual heating rates of HI, HeI, HeII
	  
      do i_subband=1,NumFreqBnd 
      !do i_subband=NumBndin1+6,NumFreqBnd
       
	   
	    phi_heat_HI = 0.0_dp
        phi_heat_HeI = 0.0_dp
        phi_heat_HeII = 0.0_dp
      
        select case (i_subband)

        ! Incoming, outcoming, current cell HI heating rate at band 1
        case (NumBndin1)

          phi%heat_in_HI = (heat_thick_table(tau_pos_in%ipos(1),1)+(heat_thick_table(tau_pos_in%ipos_p1(1),1)- &
                           heat_thick_table(tau_pos_in%ipos(1),1))*tau_pos_in%residual(1))*NFlux 

          ! When current cell is HI optically thick     
          if (abs(tau_cell_HI(i_subband)) .gt. tau_heat_limit) then
            phi%heat_out_HI = (heat_thick_table(tau_pos_out%ipos(1),1)+(heat_thick_table(tau_pos_out%ipos_p1(1),1)- &
                              heat_thick_table(tau_pos_out%ipos(1),1))*tau_pos_out%residual(1))*NFlux 
            phi_heat_HI = (phi%heat_in_HI-phi%heat_out_HI)/vol
 
          ! When current cell is HI optically thin
          else
            phi_heat_HI = (heat_thin_table(tau_pos_in%ipos(1),1)+(heat_thin_table(tau_pos_in%ipos_p1(1),1)- &
                          heat_thin_table(tau_pos_in%ipos(1),1))*tau_pos_in%residual(1))* &   
                          (tau_out_all(1)-tau_in_all(1))*NFlux 
            phi%heat_out_HI = phi%heat_in_HI+phi_heat_HI
            phi_heat_HI = phi_heat_HI/vol
          endif

          f_heat = phi_heat_HI

        ! Incoming, outcoming, current cell HI, HeI heating rate at band 2      
        case (NumBndin1+1:NumBndin1+NumBndin2)
       
          phi%heat_in_HI = (heat_thick_table(tau_pos_in%ipos(i_subband),2*i_subband-2)+ &
                           (heat_thick_table(tau_pos_in%ipos_p1(i_subband),2*i_subband-2)- &
                           heat_thick_table(tau_pos_in%ipos(i_subband),2*i_subband-2))* &
                           tau_pos_in%residual(i_subband))*NFlux

          phi%heat_in_HeI = (heat_thick_table(tau_pos_in%ipos(i_subband),2*i_subband-1)+ &
                            (heat_thick_table(tau_pos_in%ipos_p1(i_subband),2*i_subband-1)- &
                            heat_thick_table(tau_pos_in%ipos(i_subband),2*i_subband-1))* &
                            tau_pos_in%residual(i_subband))*NFlux

          ! When current cell is HI optically thick  
          if (abs(tau_cell_HI(i_subband)) .gt. tau_heat_limit) then
            phi%heat_out_HI = (heat_thick_table(tau_pos_out%ipos(i_subband),2*i_subband-2)+ &
                              (heat_thick_table(tau_pos_out%ipos_p1(i_subband),2*i_subband-2)- &
                              heat_thick_table(tau_pos_out%ipos(i_subband),2*i_subband-2))* &
                              tau_pos_out%residual(i_subband))*NFlux
            phi_heat_HI = scaling_HI(i_subband)*(phi%heat_in_HI-phi%heat_out_HI)/vol 

          ! When current cell is HI optically thin
          else
            phi_heat_HI = scaling_HI(i_subband)*(((heat_thin_table(tau_pos_in%ipos(i_subband),2*i_subband-2)+ &
                          (heat_thin_table(tau_pos_in%ipos_p1(i_subband),2*i_subband-2)- &
                          heat_thin_table(tau_pos_in%ipos_p1(i_subband),2*i_subband-2))* &
                          tau_pos_in%residual(i_subband))*(tau_out_all(i_subband)-tau_in_all(i_subband))))*NFlux
            phi%heat_out_HI = phi%heat_in_HI+phi_heat_HI
            phi_heat_HI = phi_heat_HI/vol
          endif

          ! When current cell is HeI optically thick      
          if (abs(tau_cell_HeI(i_subband)) .gt. tau_heat_limit) then
            phi%heat_out_HeI = (heat_thick_table(tau_pos_out%ipos(i_subband),2*i_subband-1)+ &
                               (heat_thick_table(tau_pos_out%ipos_p1(i_subband),2*i_subband-1)- &
                               heat_thick_table(tau_pos_out%ipos(i_subband),2*i_subband-1))* &
                               tau_pos_out%residual(i_subband))*NFlux
            phi_heat_HeI = scaling_HeI(i_subband)*(phi%heat_in_HeI-phi%heat_out_HeI)/vol

          ! When current cell is HeI optically thin
          else 
            phi_heat_HeI = scaling_HeI(i_subband)*((heat_thin_table(tau_pos_in%ipos(i_subband),2*i_subband-1)+&
                           (heat_thin_table(tau_pos_in%ipos_p1(i_subband),2*i_subband-1)- &
                           heat_thin_table(tau_pos_in%ipos_p1(i_subband),2*i_subband-1))* &
                           tau_pos_in%residual(i_subband))*(tau_out_all(i_subband)-tau_in_all(i_subband)))*NFlux
            phi%heat_out_HeI=phi%heat_in_HeI+phi_heat_HeI
            phi_heat_HeI = phi_heat_HeI/vol
          endif

          fra_sum1 = f1ion_HI(i_subband)*phi_heat_HI+f1ion_HeI(i_subband)*phi_heat_HeI
          fra_sum2 = f2ion_HI(i_subband)*phi_heat_HI+f2ion_HeI(i_subband)*phi_heat_HeI
          fra_sum3 = f1heat_HI(i_subband)*phi_heat_HI+f1heat_HeI(i_subband)*phi_heat_HeI
          fra_sum4 = f2heat_HI(i_subband)*phi_heat_HI+f2heat_HeI(i_subband)*phi_heat_HeI

     if (do_secondary_ionization.eqv..true.)then
          ! These are all cumulative
          f_ion_HeI = f_ion_HeI+y1R(2)*fra_sum1-y2R(2)*fra_sum2  
          f_ion_HI = f_ion_HI+y1R(1)*fra_sum1-y2R(1)*fra_sum2
          f_heat = f_heat+phi_heat_HI+phi_heat_HeI-y1R(3)*fra_sum3+y2R(3)*fra_sum4 
	  else
		 f_heat = f_heat+phi_heat_HI+phi_heat_HeI
	 endif 
		  
        ! Incoming, outcoming, current cell HI, HeI, HeII heating rate at band 3   
        case (NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+NumBndin3)
          phi%heat_in_HI = (heat_thick_table(tau_pos_in%ipos(i_subband),3*i_subband-NumBndin2-4)+ &
                           (heat_thick_table(tau_pos_in%ipos_p1(i_subband),3*i_subband-NumBndin2-4)- &
                           heat_thick_table(tau_pos_in%ipos(i_subband),3*i_subband-NumBndin2-4))* &
                           tau_pos_in%residual(i_subband))*NFlux 
          phi%heat_in_HeI = (heat_thick_table(tau_pos_in%ipos(i_subband),3*i_subband-NumBndin2-3)+ &
                            (heat_thick_table(tau_pos_in%ipos_p1(i_subband),3*i_subband-NumBndin2-3)- &
                            heat_thick_table(tau_pos_in%ipos(i_subband),3*i_subband-NumBndin2-3))* &
                            tau_pos_in%residual(i_subband))*NFlux
          phi%heat_in_HeII = (heat_thick_table(tau_pos_in%ipos(i_subband),3*i_subband-NumBndin2-2)+ &
                             (heat_thick_table(tau_pos_in%ipos_p1(i_subband),3*i_subband-NumBndin2-2)- &
                             heat_thick_table(tau_pos_in%ipos(i_subband),3*i_subband-NumBndin2-2))* &
                             tau_pos_in%residual(i_subband))*NFlux

          ! When current cell is HI optically thick 
          if (abs(tau_cell_HI(i_subband)) .gt. tau_heat_limit) then
            phi%heat_out_HI = (heat_thick_table(tau_pos_out%ipos(i_subband),3*i_subband-NumBndin2-4)+ &
                              (heat_thick_table(tau_pos_out%ipos_p1(i_subband),3*i_subband-NumBndin2-4)- &
                              heat_thick_table(tau_pos_out%ipos(i_subband),3*i_subband-NumBndin2-4))* &
                              tau_pos_out%residual(i_subband))*NFlux
            phi_heat_HI = scaling_HI(i_subband)*(phi%heat_in_HI-phi%heat_out_HI)/vol

          ! When current cell is HI optically thin
          else
            phi_heat_HI = scaling_HI(i_subband)*((heat_thin_table(tau_pos_in%ipos(i_subband),3*i_subband-NumBndin2-4)+ &
                          (heat_thin_table(tau_pos_in%ipos_p1(i_subband),3*i_subband-NumBndin2-4)- &
                          heat_thin_table(tau_pos_in%ipos_p1(i_subband),3*i_subband-NumBndin2-4))* &
                          tau_pos_in%residual(i_subband))*(tau_out_all(i_subband)-tau_in_all(i_subband)))*NFlux
            phi%heat_out_HI = phi%heat_in_HI+phi_heat_HI
            phi_heat_HI = phi_heat_HI/vol
          endif

          ! When current cell is HeI optically thick   
          if (abs(tau_cell_HeI(i_subband)) .gt. tau_heat_limit) then 
            phi%heat_out_HeI = (heat_thick_table(tau_pos_out%ipos(i_subband),3*i_subband-NumBndin2-3)+ &
                               (heat_thick_table(tau_pos_out%ipos_p1(i_subband),3*i_subband-NumBndin2-3)- &
                               heat_thick_table(tau_pos_out%ipos(i_subband),3*i_subband-NumBndin2-3))* &
                               tau_pos_out%residual(i_subband))*NFlux
            phi_heat_HeI = scaling_HeI(i_subband)*(phi%heat_in_HeI-phi%heat_out_HeI)/vol

          ! When current cell is HeI optically thin
          else
            phi_heat_HeI = scaling_HeI(i_subband)*((heat_thin_table(tau_pos_in%ipos(i_subband),3*i_subband-NumBndin2-3)+ &
                           (heat_thin_table(tau_pos_in%ipos_p1(i_subband),3*i_subband-NumBndin2-3)- &
                           heat_thin_table(tau_pos_in%ipos_p1(i_subband),3*i_subband-NumBndin2-3))* &
                           tau_pos_in%residual(i_subband))*(tau_out_all(i_subband)-tau_in_all(i_subband)))*NFlux
            phi%heat_out_HeI = phi%heat_in_HeI+phi_heat_HeI
            phi_heat_HeI = phi_heat_HeI/vol
          endif

          ! When current cell is HeII optically thick      
          if (abs(tau_cell_HeII(i_subband)) .gt. tau_heat_limit) then
            phi%heat_out_HeII = (heat_thick_table(tau_pos_out%ipos(i_subband),3*i_subband-NumBndin2-2)+ &
                                (heat_thick_table(tau_pos_out%ipos_p1(i_subband),3*i_subband-NumBndin2-2)- &
                                heat_thick_table(tau_pos_out%ipos(i_subband),3*i_subband-NumBndin2-2))* &
                                tau_pos_out%residual(i_subband))*NFlux
            phi_heat_HeII = scaling_HeII(i_subband)*(phi%heat_in_HeII-phi%heat_out_HeII)/vol

          ! When current cell is HeII optically thin
          else 
            phi_heat_HeII = scaling_HeII(i_subband)*((heat_thin_table(tau_pos_in%ipos(i_subband),3*i_subband-NumBndin2-2)+ &
                            (heat_thin_table(tau_pos_in%ipos_p1(i_subband),3*i_subband-NumBndin2-2)- &
                            heat_thin_table(tau_pos_in%ipos_p1(i_subband),3*i_subband-NumBndin2-2))* &
                            tau_pos_in%residual(i_subband))*(tau_out_all(i_subband)- tau_in_all(i_subband)))*NFlux
            phi%heat_out_HeII = phi%heat_in_HeII+phi_heat_HeII
            phi_heat_HeII = phi_heat_HeII/vol
          endif
          
          fra_sum1 = f1ion_HI(i_subband)*phi_heat_HI+f1ion_HeI(i_subband)*phi_heat_HeI+f1ion_HeII(i_subband)*phi_heat_HeII
          fra_sum2 = f2ion_HI(i_subband)*phi_heat_HI+f2ion_HeI(i_subband)*phi_heat_HeI+f2ion_HeII(i_subband)*phi_heat_HeII
          fra_sum3 = f1heat_HI(i_subband)*phi_heat_HI+f1heat_HeI(i_subband)*phi_heat_HeI+f1heat_HeII(i_subband)*phi_heat_HeII
          fra_sum4 = f2heat_HI(i_subband)*phi_heat_HI+f2heat_HeI(i_subband)*phi_heat_HeI+f2heat_HeII(i_subband)*phi_heat_HeII
          
      if (do_secondary_ionization.eqv..true.) then
          ! These are all cumulative
          f_ion_HeI = f_ion_HeI+y1R(2)*fra_sum1-y2R(2)*fra_sum2
          f_ion_HI = f_ion_HI+y1R(1)*fra_sum1-y2R(1)*fra_sum2
          f_heat = f_heat+phi_heat_HI+phi_heat_HeI+phi_heat_HeII-y1R(3)*fra_sum3+y2R(3)*fra_sum4   
	  else
		  f_heat = f_heat+phi_heat_HI+phi_heat_HeI+phi_heat_HeII
	  endif
        end select
	
      enddo 

      !Total heating rate on current cell
      phi%heat = f_heat 
      !Final HI photoionization rate modified by secondary ionization
      phi%photo_cell_HI = phi%photo_cell_HI+f_ion_HI/(ion_freq_HI*hplanck) 
      !Final HeI photoionization rate modified by secondary ionization
      phi%photo_cell_HeI = phi%photo_cell_HeI+f_ion_HeI/(ion_freq_HeI*hplanck)  

    endif 
	
  end subroutine lookuptable

  ! give scalings of species for division of photoionization and heating to species
  subroutine scale_int2(scaling_HI,scaling_HeI,colum_cell_HI,colum_cell_HeI,i_subband)

    real(kind=dp),intent(in) :: colum_cell_HI, colum_cell_HeI
    integer,intent(in) :: i_subband
    real(kind=dp),intent(out):: scaling_HI, scaling_HeI
    real(kind=dp) :: forscaleing

    forscaleing = 1.0_dp/(sigma_HI(i_subband)*colum_cell_HI+sigma_HeI(i_subband)*colum_cell_HeI)
    scaling_HI = sigma_HI(i_subband)*colum_cell_HI*forscaleing
    scaling_HeI = sigma_HeI(i_subband)*colum_cell_HeI*forscaleing

  end subroutine scale_int2  

  ! give scalings of species for division of photoionization and heating to species
  subroutine scale_int3(scaling_HI, scaling_HeI, scaling_HeII, colum_cell_HI, colum_cell_HeI, colum_cell_HeII, i_subband)

    real(kind=dp),intent(in) :: colum_cell_HI ,colum_cell_HeI, colum_cell_HeII
    integer,intent(in) :: i_subband
    real(kind=dp),intent(out):: scaling_HI, scaling_HeI, scaling_HeII
    real(kind=dp) :: forscaleing

    forscaleing = 1.0_dp/(sigma_HI(i_subband)*colum_cell_HI+sigma_HeI(i_subband)*colum_cell_HeI+ &
                  sigma_HeII(i_subband)*colum_cell_HeII) 
    scaling_HI = colum_cell_HI*sigma_HI(i_subband) *forscaleing
    scaling_HeI = colum_cell_HeI*sigma_HeI(i_subband)*forscaleing
    scaling_HeII = colum_cell_HeII*sigma_HeII(i_subband)*forscaleing 

  end subroutine scale_int3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! set up some scaling factors arrays
  subroutine setup_scalingfactors

    integer :: i_subband

    ! Allocate size of arrays
    allocate(delta_freq(1:NumFreqBnd))  
    allocate(freq_max(1:NumFreqBnd))  
    allocate(freq_min(1:NumFreqBnd))  
    allocate(pl_index_cross_section_HI(1:NumFreqBnd))
    allocate(pl_index_cross_section_HeI(1:NumFreqBnd))
    allocate(pl_index_cross_section_HeII(1:NumFreqBnd))
    allocate(sigma_HI(1:NumFreqBnd))
    allocate(sigma_HeI(1:NumFreqBnd))
    allocate(sigma_HeII(1:NumFreqBnd))

    ! Allocate size of arrays of heating parameters
    if (.not.isothermal) then
      allocate(f1ion_HI(NumBndin1+1:NumFreqBnd))  
      allocate(f1ion_HeI(NumBndin1+1:NumFreqBnd))  
      allocate(f1ion_HeII(NumBndin1+1:NumFreqBnd))  
      allocate(f2ion_HI(NumBndin1+1:NumFreqBnd)) 
      allocate(f2ion_HeI(NumBndin1+1:NumFreqBnd)) 
      allocate(f2ion_HeII(NumBndin1+1:NumFreqBnd)) 
      allocate(f2heat_HI(NumBndin1+1:NumFreqBnd)) 
      allocate(f2heat_HeI(NumBndin1+1:NumFreqBnd)) 
      allocate(f2heat_HeII(NumBndin1+1:NumFreqBnd)) 
      allocate(f1heat_HI(NumBndin1+1:NumFreqBnd)) 
      allocate(f1heat_HeI(NumBndin1+1:NumFreqBnd))
      allocate(f1heat_HeII(NumBndin1+1:NumFreqBnd))
    endif

    ! Assignment of maximum frequency in the sub-bin partition.
    select case (NumBndin1)

      case (1)

        freq_max(NumBndin1) = ion_freq_HeI

    end select

    select case (NumBndin2)

    case (26)
    
      freq_max(NumBndin1+1:26) = ion_freq_HeI* &
                                 (/1.02_dp, 1.05_dp, 1.07_dp, 1.10_dp, 1.15_dp, 1.20_dp, &
		  	           1.25_dp, 1.30_dp, 1.35_dp, 1.40_dp, 1.45_dp, 1.50_dp, &
		                   1.55_dp, 1.60_dp, 1.65_dp, 1.70_dp, 1.75_dp, 1.80_dp, &
		                   1.85_dp, 1.90_dp, 1.95_dp, 2.00_dp, 2.05_dp, 2.10_dp, &
                                   2.15_dp/)
      freq_max(NumBndin1+NumBndin2) = ion_freq_HeII

    case (10)
	
      freq_max(NumBndin1+1:10) = ion_freq_HeI* &
                                 (/1.10_dp, 1.20_dp, 1.30_dp, 1.40_dp, 1.50_dp, &
                                   1.60_dp, 1.70_dp, 1.80_dp, 1.90_dp/)
      freq_max(NumBndin1+NumBndin2) = ion_freq_HeII

    case (6)

      freq_max(NumBndin1+1:6) = ion_freq_HeI* &
                                (/1.15_dp, 1.30_dp, 1.50_dp, 1.70_dp, 1.9557_dp/)
      freq_max(NumBndin1+NumBndin2) = ion_freq_HeII

    case (3)

      freq_max(NumBndin1+1:3) = ion_freq_HeI* &                           
                                (/1.3_dp,1.7_dp/)
      freq_max(NumBndin1+NumBndin2) = ion_freq_HeII

    case (2)   
  
      freq_max(NumBndin1+1) = ion_freq_HeI*1.5_dp
      freq_max(NumBndin1+NumBndin2) = ion_freq_HeII

    case (1)  
      freq_max(NumBndin1+NumBndin2) = ion_freq_HeII

    end select

    select case (NumBndin3)

      case (20)
        freq_max(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = ion_freq_HeII* &
                                (/1.05_dp, 1.10_dp, 1.20_dp, 1.40_dp, 1.70_dp, 2.00_dp, &
			          2.50_dp, 3.00_dp, 4.00_dp, 5.00_dp, 7.00_dp, 10.00_dp, &
			          15.00_dp, 20.00_dp, 30.00_dp, 40.00_dp, 50.00_dp, 70.00_dp, &
                                  90.00_dp, 100.00_dp/)  

      case (16)
        freq_max(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+16) = ion_freq_HeII* &
			        (/1.05_dp, 1.10_dp, 1.20_dp, 1.40_dp, 1.70_dp, 2.00_dp, &
                                  3.00_dp, 5.00_dp, 7.00_dp, 10.00_dp, 15.00_dp, 20.00_dp, &
                                  30.00_dp, 50.00_dp, 70.00_dp, 100.00_dp/)

      case (11)
        freq_max(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+11) = ion_freq_HeII* &
                                (/1.10_dp, 1.20_dp, 1.50_dp, 2.00_dp, 3.00_dp, 4.00_dp, &
                                  7.00_dp, 10.00_dp, 20.00_dp, 50.00_dp, 100.0_dp/)

      case (9)
        freq_max(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+9) = ion_freq_HeII* & 
	                        (/1.50_dp, 2.00_dp, 3.00_dp, 4.00_dp, 7.00_dp, 10.00_dp, &
                                  20.00_dp, 50.00_dp, 100.00_dp/)

      case (4)
        freq_max(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+4) = ion_freq_HeII* &
                                (/2.00_dp, 4.00_dp, 10.00_dp, 100.0_dp/)

      case (1)
        freq_max(NumBndin1+NumBndin2+1) = ion_freq_HeII * 100.00_dp

    end select

    ! Assignment of minimum frequency in the sub-bin partition.
       
    freq_min(NumBndin1) = ion_freq_HI
	
    do i_subband=2,NumFreqBnd
      freq_min(i_subband) = freq_max(i_subband-1)
    enddo

    ! calculate the width of frequency sub-bin
    do i_subband=1,NumFreqBnd 
      delta_freq(i_subband) = (freq_max(i_subband)-freq_min(i_subband))/real(NumFreq)
    enddo

    ! Assign f_ion and f_heat for secondary ionization
    if (.not.isothermal) then

      select case (NumBndin2)

        case (26)

          f1ion_HI(NumBndin1+1:NumBndin1+26) = (/ 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                                                  1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                                                  1.0000_dp, 1.0000_dp/) 

          f1ion_HeI(NumBndin1+1:NumBndin1+26) = (/0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 1.0000_dp/) 

          f1ion_HeII(NumBndin1+1:NumBndin1+26) = (/0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp/) 

          f2ion_HI(NumBndin1+1:NumBndin1+26) = (/0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                 0.9971_dp, 0.9802_dp, 0.9643_dp, 0.9493_dp, &
                                                 0.9350_dp, 0.9215_dp, 0.9086_dp, 0.8964_dp, &
                                                 0.8847_dp, 0.8735_dp/)

          f2ion_HeI(NumBndin1+1:NumBndin1+26) = (/0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 0.9960_dp/) 

          f2ion_HeII(NumBndin1+1:NumBndin1+26) = (/0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp/) 

          f1heat_HI(NumBndin1+1:NumBndin1+26) = (/0.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                                                    1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                                                    1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                                                    1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                                                    1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                                                    1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                                                    1.0000_dp, 1.0000_dp/)

          f1heat_HeI(NumBndin1+1:NumBndin1+26) = (/0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp, 0.0000_dp, 1.0000_dp, &
                                                   1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                                                   1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                                                   1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                                                   1.0000_dp, 1.0000_dp /)

          f1heat_HeII(NumBndin1+1:NumBndin1+26) = (/0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                    0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                    0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                    0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                    0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                    0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                    0.0000_dp, 0.0000_dp /)

          f2heat_HI(NumBndin1+1:NumBndin1+26) = (/0.0000_dp, 0.9704_dp, 0.9290_dp, 0.9037_dp, &
                                                  0.8687_dp, 0.8171_dp, 0.7724_dp, 0.7332_dp, &
                                                  0.6985_dp, 0.6675_dp, 0.6397_dp, 0.6145_dp, &
                                                  0.5916_dp, 0.5707_dp, 0.5514_dp, 0.5337_dp, &
                                                  0.5173_dp, 0.5021_dp, 0.4879_dp, 0.4747_dp, &
                                                  0.4623_dp, 0.4506_dp, 0.4397_dp, 0.4293_dp, &
                                                  0.4196_dp, 0.4103_dp/) 

          f2heat_HeI(NumBndin1+1:NumBndin1+26) = (/0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp, 0.0000_dp, 0.9959_dp, &
                                                   0.9250_dp, 0.8653_dp, 0.8142_dp, 0.7698_dp, &
                                                   0.7309_dp, 0.6965_dp, 0.6657_dp, 0.6380_dp, &
                                                   0.6130_dp, 0.5903_dp, 0.5694_dp, 0.5503_dp, &
                                                   0.5327_dp, 0.5164_dp/)

          f2heat_HeII(NumBndin1+1:NumBndin1+26) = (/0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                    0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                    0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                    0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                    0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                    0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                    0.0000_dp, 0.0000_dp/)  

      end select

      select case (NumBndin3)

        case(20)

          f1ion_HI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                        (/1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp /)

           f1ion_HeI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                        (/1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp /)

           f1ion_HeII(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                        (/0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp /)

           f2ion_HI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) =&
                        (/0.8600_dp, 0.8381_dp, 0.8180_dp, 0.7824_dp, 0.7249_dp, 0.6607_dp, &
                          0.6128_dp, 0.5542_dp, 0.5115_dp, 0.4518_dp, 0.4110_dp, 0.3571_dp, &
                          0.3083_dp, 0.2612_dp, 0.2325_dp, 0.1973_dp, 0.1757_dp, 0.1606_dp, &
                          0.1403_dp, 0.1269_dp /)

           f2ion_HeI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                        (/0.9750_dp, 0.9415_dp, 0.9118_dp, 0.8609_dp, 0.7831_dp, 0.7015_dp, &
                          0.6436_dp, 0.5755_dp, 0.5273_dp, 0.4619_dp, 0.4182_dp, 0.3615_dp, &
                          0.3109_dp, 0.2627_dp, 0.2334_dp, 0.1979_dp, 0.1761_dp, 0.1609_dp, &
                          0.1405_dp, 0.1270_dp /)

           f2ion_HeII(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                        (/0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.8841_dp, &
                          0.7666_dp, 0.6518_dp, 0.5810_dp, 0.4940_dp, 0.4403_dp, 0.3744_dp, &
                          0.3183_dp, 0.2668_dp, 0.2361_dp, 0.1993_dp, 0.1771_dp, 0.1616_dp, &
                          0.1409_dp, 0.1273_dp /)

           f1heat_HI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                        (/1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp /)

           f1heat_HeI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                        (/1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp /)

           f1heat_HeII(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                        (/0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp /)

           f2heat_HI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                        (/0.3994_dp, 0.3817_dp, 0.3659_dp, 0.3385_dp, 0.2961_dp, 0.2517_dp, &
                          0.2207_dp, 0.1851_dp, 0.1608_dp, 0.1295_dp, 0.1097_dp, 0.0858_dp, &
                          0.0663_dp, 0.0496_dp, 0.0405_dp, 0.0304_dp, 0.0248_dp, 0.0212_dp, &
                          0.0167_dp, 0.0140_dp /)

           f2heat_HeI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                        (/0.4974_dp, 0.4679_dp, 0.4424_dp, 0.4001_dp, 0.3389_dp, 0.2796_dp, &
                          0.2405_dp, 0.1977_dp, 0.1697_dp, 0.1346_dp, 0.1131_dp, 0.0876_dp, &
                          0.0673_dp, 0.0501_dp, 0.0408_dp, 0.0305_dp, 0.0249_dp, 0.0213_dp, &
                          0.0168_dp, 0.0140_dp /)

           f2heat_HeII(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                        (/0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.6202_dp, 0.4192_dp, &
                          0.3265_dp, 0.2459_dp, 0.2010_dp, 0.1513_dp, 0.1237_dp, 0.0932_dp, &
                          0.0701_dp, 0.0515_dp, 0.0416_dp, 0.0309_dp, 0.0251_dp, 0.0214_dp, &
                          0.0169_dp, 0.0141_dp /)

      end select 

    endif

    ! Assign value sigma of HI, HeI and HeII at different frequencies 
    select case (NumBndin1)

      case (1)

        sigma_HI(1) = sigma_HI_at_ion_freq
        sigma_HeI(1) = 0.0_dp
        sigma_HeII(1) = 0.0_dp  

      end select

    select case (NumBndin2)

      case (26) 

        sigma_HI(NumBndin1+1:NumBndin1+26) = (/1.239152e-18_dp, 1.171908e-18_dp, &
             1.079235e-18_dp, 1.023159e-18_dp, 9.455687e-19_dp, 8.329840e-19_dp, &
             7.374876e-19_dp, 6.559608e-19_dp, 5.859440e-19_dp, 5.254793e-19_dp, &
             4.729953e-19_dp, 4.272207e-19_dp, 3.874251e-19_dp, 3.521112e-19_dp, &
             3.209244e-19_dp, 2.932810e-19_dp, 2.686933e-19_dp, 2.467523e-19_dp, &
             2.271125e-19_dp, 2.094813e-19_dp, 1.936094e-19_dp, 1.792838e-19_dp, &
             1.663215e-19_dp, 1.545649e-19_dp, 1.438778e-19_dp, 1.341418e-19_dp/)
        sigma_HeI(NumBndin1+1:NumBndin1+26) = (/7.434699e-18_dp, 7.210641e-18_dp, &
              6.887151e-18_dp, 6.682491e-18_dp, 6.387263e-18_dp, 5.931487e-18_dp, &
              5.516179e-18_dp, 5.137743e-18_dp, 4.792724e-18_dp, 4.477877e-18_dp, &
              4.190200e-18_dp, 3.926951e-18_dp, 3.687526e-18_dp, 3.465785e-18_dp, &
              3.261781e-18_dp, 3.073737e-18_dp, 2.900074e-18_dp, 2.739394e-18_dp, &
              2.590455e-18_dp, 2.452158e-18_dp, 2.323526e-18_dp, 2.203694e-18_dp, &
              2.091889e-18_dp, 1.987425e-18_dp, 1.889687e-18_dp, 1.798126e-18_dp/)
        sigma_HeII(NumBndin1+1:NumBndin1+26) = 0.0_dp 

      case (10) 

        sigma_HI(NumBndin1+1:NumBndin1+10) = (/1.239152e-18_dp, 9.455687e-19_dp, &
             7.374876e-19_dp, 5.859440e-19_dp, 4.729953e-19_dp, 3.874251e-19_dp, &
             3.209244e-19_dp, 2.686933e-19_dp, 2.271125e-19_dp, 1.936094e-19_dp/)
        sigma_HeI(NumBndin1+1:NumBndin1+10) = (/7.434699e-18_dp, 6.387263e-18_dp, &
              5.516179e-18_dp, 4.792724e-18_dp, 4.190200e-18_dp, 3.687526e-18_dp, &
              3.261781e-18_dp, 2.900074e-18_dp, 2.590455e-18_dp, 2.323526e-18_dp/)
        sigma_HeII(NumBndin1+1:NumBndin1+10) = 0.0_dp 

      case (6)
        sigma_HI(NumBndin1+1:NumBndin1+6) = (/1.164e-18_dp, 8.33e-19_dp,5.859e-19_dp, &
                                              3.874e-19_dp,2.687e-19_dp,1.777e-19_dp/)
        sigma_HeI(NumBndin1+1:NumBndin1+6) = (/sigma_HeI_at_ion_freq, 5.9315e-18_dp, &
                                                       4.7927e-18_dp, 3.6875e-18_dp, &
                                                       2.9001e-18_dp, 2.1906e-18_dp/)
        sigma_HeII(NumBndin1+1:NumBndin1+6) = 0.0_dp 

      case (3) 
        sigma_HI(NumBndin1+1:NumBndin1+3) = (/1.239e-18_dp, 5.86e-19_dp, 2.69e-19_dp/)
        sigma_HeI(NumBndin1+1:NumBndin1+3) = (/sigma_HeI_at_ion_freq, 4.793e-18_dp,2.90e-18_dp/)  
        sigma_HeII(NumBndin1+1:NumBndin1+3) = 0.0_dp 

      case (2)
        sigma_HI(NumBndin1+1:NumBndin1+2) = (/1.239e-18_dp, 3.87e-19_dp/)
        sigma_HeI(NumBndin1+1:NumBndin1+2) = (/sigma_HeI_at_ion_freq, 3.688e-18_dp/)  
        sigma_HeII(NumBndin1+1:NumBndin1+2) = 0.0_dp 

      case (1) 
        sigma_HI(NumBndin1+1) = 1.239e-18_dp
        sigma_HeI(NumBndin1+1) = sigma_HeI_at_ion_freq
        sigma_HeII(NumBndin1+1) = 0.0_dp 

    end select

    select case (NumBndin3)

      case (20)

        sigma_HI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                 (/1.230696e-19_dp, 1.063780e-19_dp, 9.253883e-20_dp, &
                   7.123014e-20_dp, 4.464019e-20_dp, 2.465533e-20_dp, &
                   1.492667e-20_dp, 7.446712e-21_dp, 4.196728e-21_dp, &
                   1.682670e-21_dp, 8.223247e-22_dp, 2.763830e-22_dp, &
                   8.591126e-23_dp, 2.244684e-23_dp, 8.593853e-24_dp, &
                   2.199718e-24_dp, 8.315674e-25_dp, 3.898672e-25_dp, &
                   1.238718e-25_dp, 5.244957e-26_dp/)
        sigma_HeI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                 (/1.690781e-18_dp, 1.521636e-18_dp, 1.373651e-18_dp, &
                   1.128867e-18_dp, 7.845096e-19_dp, 4.825331e-19_dp, &
                   3.142134e-19_dp, 1.696228e-19_dp, 1.005051e-19_dp, &
                   4.278712e-20_dp, 2.165403e-20_dp, 7.574790e-21_dp, &
                   2.429426e-21_dp, 6.519748e-22_dp, 2.534069e-22_dp, &
                   6.599821e-23_dp, 2.520412e-23_dp, 1.189810e-23_dp, &
                   3.814490e-24_dp, 1.624492e-24_dp/)
        sigma_HeII(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                 (/1.587280e-18_dp, 1.391911e-18_dp, 1.227391e-18_dp, &
                   9.686899e-19_dp, 6.338284e-19_dp, 3.687895e-19_dp, &
                   2.328072e-19_dp, 1.226873e-19_dp, 7.214988e-20_dp, &
                   3.081577e-20_dp, 1.576429e-20_dp, 5.646276e-21_dp, &
                   1.864734e-21_dp, 5.177347e-22_dp, 2.059271e-22_dp, &
                   5.526508e-23_dp, 2.151467e-23_dp, 1.029637e-23_dp, &
                   3.363164e-24_dp, 1.450239e-24_dp/)

      case (16) 
	
        sigma_HI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+16) = &
                (/ 1.230696e-19_dp, 1.063780e-19_dp, 9.253883e-20_dp, &
                   7.123014e-20_dp, 4.464019e-20_dp, 2.465533e-20_dp, &	
                   1.492667e-20_dp, 4.196728e-21_dp, 8.223247e-22_dp, &
                   2.763830e-22_dp, 8.591126e-23_dp, 2.244684e-23_dp, &	
                   8.593853e-24_dp, 2.199718e-24_dp, 3.898672e-25_dp, &
                   1.238718e-25_dp/)

        sigma_HeI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+16) = &
                 (/1.690781e-18_dp, 1.521636e-18_dp, 1.373651e-18_dp, &
                   1.128867e-18_dp, 7.845096e-19_dp, 4.825331e-19_dp, &
                   3.142134e-19_dp, 1.005051e-19_dp, 2.165403e-20_dp, &
                   7.574790e-21_dp, 2.429426e-21_dp, 6.519748e-22_dp, &
                   2.534069e-22_dp, 6.599821e-23_dp, 1.189810e-23_dp, &
                   3.814490e-24_dp/)

        sigma_HeII(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+16) = &
                 (/1.587280e-18_dp, 1.391911e-18_dp, 1.227391e-18_dp, &
                   9.686899e-19_dp, 6.338284e-19_dp, 3.687895e-19_dp, &
                   2.328072e-19_dp, 7.214988e-20_dp, 1.576429e-20_dp, &
                   5.646276e-21_dp, 1.864734e-21_dp, 5.177347e-22_dp, &
                   2.059271e-22_dp, 5.526508e-23_dp, 1.029637e-23_dp, &
                   3.363164e-24_dp/)

      case (11) 

        sigma_HI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+11) = &
                       (/1.2307e-19_dp, 9.2539e-20_dp, 7.1230e-20_dp, &
                         3.6176e-20_dp, 1.4927e-20_dp, 4.1967e-21_dp, &
                         1.6827e-21_dp, 2.7638e-22_dp, 8.5911e-23_dp, &
                         8.5939e-24_dp,3.8987e-25_dp/)
        sigma_HeI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+11) = &
                       (/1.6908e-18_dp, 1.3737e-18_dp, 1.1289e-18_dp, &
                         6.6238e-19_dp, 3.1421e-19_dp, 1.0051e-19_dp, &
                         4.2787e-20_dp, 7.5748e-21_dp, 2.4294e-21_dp, &
                         2.5341e-22_dp, 1.1898e-23_dp/)
        sigma_HeII(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+11) = &
                       (/1.5873e-18_dp, 1.2274e-18_dp, 9.6869e-19_dp, &
                         5.2339e-19_dp, 2.3281e-19_dp, 7.2150e-20_dp, &
                         3.0816e-20_dp, 5.6463e-21_dp, 1.8647e-21_dp, &
                         2.0593e-22_dp, 1.0296e-23_dp/)

      case (9) 

        sigma_HI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+9) = &
                 (/1.230696e-19_dp, 3.617600e-20_dp, 1.492667e-20_dp, &
                   4.196728e-21_dp, 1.682670e-21_dp, 2.763830e-22_dp, &
                   8.591126e-23_dp, 8.593853e-24_dp, 3.898672e-25_dp/) 
        sigma_HeI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+9) = &
                 (/1.690781e-18_dp, 6.623773e-19_dp, 3.142134e-19_dp, &
                   1.005051e-19_dp, 4.278712e-20_dp, 7.574790e-21_dp, &
                   2.429426e-21_dp, 2.534069e-22_dp, 1.189810e-23_dp/)
        sigma_HeII(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+9) = &
           (/sigma_HeII_at_ion_freq,5.233870e-19_dp, 2.328072e-19_dp, &
                   7.214988e-20_dp, 3.081577e-20_dp, 5.646276e-21_dp, &
                   1.864734e-21_dp, 2.059271e-22_dp, 1.029637e-23_dp/)

      case (4)

        sigma_HI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+4) = &
          (/1.2307e-19_dp, 1.4927e-20_dp, 1.6827e-21_dp, 8.5900e-23_dp/)
        sigma_HeI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+4) = &
          (/1.6908e-18_dp, 3.1421e-19_dp, 4.2787e-20_dp, 2.4294e-21_dp/)
        sigma_HeII(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+4) = &
          (/1.5873e-18_dp, 2.3280e-19_dp, 3.0816e-20_dp, 1.1865e-21_dp/)

      case (1)

        sigma_HI(NumBndin1+NumBndin2+1) = 1.2300e-19_dp
        sigma_HeI(NumBndin1+NumBndin2+1) = 1.691e-18_dp
        sigma_HeII(NumBndin1+NumBndin2+1) = sigma_HeII_at_ion_freq

    end select


    ! Assign power-law index of HI, HeI and HeII at different frequencies (about absorption)
    select case (NumBndin1)

      case (1)

        pl_index_cross_section_HI(1) = 2.761_dp

    end select

    select case (NumBndin2)

      case (26) 

        pl_index_cross_section_HI(NumBndin1+1:NumBndin1+26) = (/2.8277_dp, 2.8330_dp, 2.8382_dp, &
                            2.8432_dp, 2.8509_dp, 2.8601_dp, 2.8688_dp, 2.8771_dp, &
                            2.8850_dp, 2.8925_dp, 2.8997_dp, 2.9066_dp, 2.9132_dp, &
                            2.9196_dp, 2.9257_dp, 2.9316_dp, 2.9373_dp, 2.9428_dp, &
                            2.9481_dp, 2.9532_dp, 2.9582_dp, 2.9630_dp, 2.9677_dp, &
                            2.9722_dp, 2.9766_dp, 2.9813_dp/)
        pl_index_cross_section_HeI(NumBndin1+1:NumBndin1+26) = (/1.5509_dp, 1.5785_dp, 1.6047_dp, &
                             1.6290_dp, 1.6649_dp, 1.7051_dp, 1.7405_dp, 1.7719_dp, &
                             1.8000_dp, 1.8253_dp, 1.8486_dp, 1.8701_dp, 1.8904_dp, &
                             1.9098_dp, 1.9287_dp, 1.9472_dp, 1.9654_dp, 1.9835_dp, &
                             2.0016_dp, 2.0196_dp, 2.0376_dp, 2.0557_dp, 2.0738_dp, &
                             2.0919_dp, 2.1099_dp, 2.1302_dp/)

      case (10) 

        pl_index_cross_section_HI(NumBndin1+1:NumBndin1+10) = (/2.8360_dp, 2.8554_dp, 2.8729_dp, 2.8887_dp, &
                                2.9031_dp,2.9164_dp, 2.9287_dp,2.9400_dp,2.9507_dp,2.9701_dp/)
        pl_index_cross_section_HeI(NumBndin1+1:NumBndin1+10) = (/1.5932_dp, 1.6849_dp, 1.7561_dp, 1.8126_dp, &
                             1.8592_dp, 1.9000_dp, 1.9379_dp, 1.9744_dp, 2.0105_dp, 2.0840_dp/)

      case (6)

        pl_index_cross_section_HI(NumBndin1+1:NumBndin1+6) = (/2.8408_dp, 2.8685_dp, 2.8958_dp, &
                                                 2.9224_dp, 2.9481_dp, 2.9727_dp/)
        pl_index_cross_section_HeI(NumBndin1+1:NumBndin1+6) = (/1.6168_dp, 1.7390_dp, 1.8355_dp, &
                                                  1.9186_dp, 2.0018_dp, 2.0945_dp/)

      case (3) 

        pl_index_cross_section_HI(NumBndin1+1:NumBndin1+3) = (/2.8542_dp, 2.9086_dp, 2.9600_dp/)
        pl_index_cross_section_HeI(NumBndin1+1:NumBndin1+3) = (/1.6770_dp, 1.8758_dp, 2.0458_dp/)

      case (2)

        pl_index_cross_section_HI(NumBndin1+1:NumBndin1+2) = (/2.8697_dp, 2.9486_dp/)
        pl_index_cross_section_HeI(NumBndin1+1:NumBndin1+2) = (/1.7385_dp, 2.0061_dp/)

      case (1) 

        pl_index_cross_section_HI(NumBndin1+1) = 2.9118_dp
        pl_index_cross_section_HeI(NumBndin1+1) = 1.8832_dp

    end select

    select case (NumBndin3)

      case (20)

        pl_index_cross_section_HI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                            (/2.9884_dp, 2.9970_dp, 3.0088_dp, 3.0298_dp, 3.0589_dp, &
                              3.0872_dp, 3.1166_dp, 3.1455_dp, 3.1773_dp, 3.2089_dp, &
                              3.2410_dp, 3.2765_dp, 3.3107_dp, 3.3376_dp, 3.3613_dp, &
                              3.3816_dp, 3.3948_dp, 3.4078_dp, 3.4197_dp, 3.4379_dp/)
        pl_index_cross_section_HeI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                            (/2.1612_dp, 2.2001_dp, 2.2564_dp, 2.3601_dp, 2.5054_dp, &
                              2.6397_dp, 2.7642_dp, 2.8714_dp, 2.9700_dp, 3.0528_dp, &
                              3.1229_dp, 3.1892_dp, 3.2451_dp, 3.2853_dp, 3.3187_dp, &
                              3.3464_dp, 3.3640_dp, 3.3811_dp, 3.3967_dp, 3.4203_dp/)
        pl_index_cross_section_HeII(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                            (/2.6930_dp, 2.7049_dp, 2.7213_dp, 2.7503_dp, 2.7906_dp, &
                              2.8300_dp, 2.8711_dp, 2.9121_dp, 2.9577_dp, 3.0041_dp, &
                              3.0522_dp, 3.1069_dp, 3.1612_dp, 3.2051_dp, 3.2448_dp, &
                              3.2796_dp, 3.3027_dp, 3.3258_dp, 3.3472_dp, 3.3805_dp/)

      case (16) 
	
        pl_index_cross_section_HI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+16) = &
                            (/2.9884_dp, 2.9970_dp, 3.0088_dp, 3.0298_dp, 3.0589_dp, &
                              3.0872_dp, 3.1303_dp, 3.1920_dp, 3.2410_dp, 3.2765_dp, &
                              3.3107_dp, 3.3376_dp, 3.3613_dp, 3.3878_dp, 3.4078_dp, &
                              3.4343_dp/)
        pl_index_cross_section_HeI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+16) = &
                            (/2.1612_dp, 2.2001_dp, 2.2564_dp, 2.3601_dp, 2.5054_dp, &
                              2.6397_dp, 2.8157_dp, 3.0093_dp, 3.1229_dp, 3.1892_dp, &
                              3.2451_dp, 3.2853_dp, 3.3187_dp, 3.3546_dp, 3.3811_dp, &
                              3.4157_dp/)
        pl_index_cross_section_HeII(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+16) = &
                            (/2.6930_dp, 2.7049_dp, 2.7213_dp, 2.7503_dp, 2.7906_dp, &
                              2.8300_dp, 2.8904_dp, 2.9793_dp, 3.0522_dp, 3.1069_dp, &
                              3.1612_dp, 3.2051_dp, 3.2448_dp, 3.2904_dp, 3.3258_dp, &
                              3.3740_dp/)

      case (11) 

        pl_index_cross_section_HI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+11) = &
                            (/2.9926_dp, 3.0088_dp, 3.0357_dp, 3.0777_dp, 3.1303_dp, &
                              3.1773_dp, 3.2292_dp, 3.2765_dp, 3.3230_dp, 3.3775_dp, &
                              3.4155_dp/)
        pl_index_cross_section_HeI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+11) = &
                            (/2.1803_dp, 2.2564_dp, 2.3901_dp, 2.5951_dp, 2.8157_dp, &
                              2.9700_dp, 3.0976_dp, 3.1892_dp, 3.2636_dp, 3.3407_dp, &
                              3.3913_dp/)
        pl_index_cross_section_HeII(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+11) = &
                            (/2.6989_dp, 2.7213_dp, 2.7585_dp, 2.8167_dp, 2.8904_dp, &
                              2.9577_dp, 3.0345_dp, 3.1069_dp, 3.1811_dp, 3.2727_dp, &
                              3.3397_dp/)

      case (9) 

        pl_index_cross_section_HI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+9) = &
                            (/3.0207_dp, 3.0777_dp, 3.1303_dp, 3.1773_dp, 3.2292_dp, &
                              3.2765_dp, 3.3230_dp, 3.3775_dp, 3.4155_dp/)
        pl_index_cross_section_HeI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+9) = &
                            (/2.3157_dp, 2.5951_dp, 2.8157_dp, 2.9700_dp, 3.0976_dp,&
                              3.1892_dp, 3.2636_dp, 3.3407_dp, 3.3913_dp/)
        pl_index_cross_section_HeII(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+9) = &
                            (/2.7377_dp, 2.8167_dp, 2.8904_dp, 2.9577_dp, 3.0345_dp,&
                              3.1069_dp, 3.1811_dp, 3.2727_dp, 3.3397_dp/)

      case (4)

        pl_index_cross_section_HI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+4) = &
                                       (/3.0465_dp, 3.1516_dp, 3.2501_dp, 3.3833_dp/)
        pl_index_cross_section_HeI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+4) = &
                                       (/2.4431_dp, 2.8878_dp, 3.1390_dp, 3.3479_dp/)
        pl_index_cross_section_HeII(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+4) = &
                                       (/2.7735_dp, 2.9209_dp, 3.0663_dp, 3.2833_dp/)

      case (1)

        pl_index_cross_section_HI(NumBndin1+NumBndin2+1) = 3.3369_dp
        pl_index_cross_section_HeI(NumBndin1+NumBndin2+1) = 3.2681_dp
        pl_index_cross_section_HeII(NumBndin1+NumBndin2+1) = 3.2082_dp

    end select

  end subroutine setup_scalingfactors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module radiation
