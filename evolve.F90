! This module contains routines having to do with the calculation of
! the ionization evolution of the entire grid (1D).

module evolve

  use precision, only: dp
  use file_admin, only: logf
  use my_mpi 
  use sizes, only: Ndim, mesh
  use grid, only: vol, dr
  use material, only: ndens, xHI, xHII, xHeI, xHeII, xHeIII, temper, HI_photoionization_rate,&
                      recombination_time,recombination_rate
  use photonstatistics, only: state_before, calculate_photon_statistics, photon_loss
  use abundances, only: abu_he,abu_h
  use astroconstants, only: YEAR
  use cgsconstants, only: ev2erg
  
  implicit none

  public :: evolve1D 

  ! HI column density to the end of the cell
  real(kind=dp), dimension(mesh), save :: coldens_out_HI
  ! HeI column density to the end of the cell
  real(kind=dp), dimension(mesh), save :: coldens_out_HeI
  ! HeII column density to the end of the cell
  real(kind=dp), dimension(mesh), save :: coldens_out_HeII
  real(kind=dp) :: alpha_HI_B
  integer, save :: report_counter = 0
  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Calculates the evolution of HI, HeI and HeII ionization state
  ! and temperature for the whole 1D grid
  subroutine evolve1D (dt,starting_pos,end_pos)

    ! The time step 
    real(kind=dp),intent(in) :: dt
    ! Evolve starts from starting_pos to end_pos
    integer,intent(in) :: starting_pos, end_pos
    ! Loop variables
    integer :: i_pos
    ! become true if photons remained are too less, stop further evolution
    logical,save :: too_less_photon
	
    ! Initial state (for photon statistics)
    call state_before ()

    ! reset photon loss counter
    photon_loss = 0.0

    ! Photons at the beginning may not be too less
    too_less_photon = .false.

    ! Loop through grid
    do i_pos = starting_pos,end_pos

      if (too_less_photon.eqv..true.) exit
      call evolve0D(dt,i_pos,too_less_photon)

    end do

    ! Calculate photon statistics
    call calculate_photon_statistics (dt)

    return

  end subroutine evolve1D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Calculates the evolution of the hydrogen ionization state in a single cell
  subroutine evolve0D(dt,pos,too_less_photon)

    use mathconstants, only: pi
    use tped, only: electrondens
    use doric_module, only: doric, coldens, &
                            coldens_bndry_HI, coldens_bndry_HeI,coldens_bndry_HeII
    use radiation, only: photoion, photrates, source_ionzing_photon_rate
    use c2ray_parameters, only: minimum_fractional_change, minimum_fraction_of_atoms, minium_fraction_of_photons
    use material, only: isothermal, ionstates, gamma_uvb, avg_ex_ph_eV
    use thermalevolution, only: thermal
    use cgsconstants,only : ini_rec_colion_factors
    use cgsphotoconstants, only: sigma_H_heth, sigma_H_heLya, sigma_He_heLya, sigma_HeI_at_ion_freq, &
                                 sigma_He_he2, sigma_HeII_at_ion_freq, sigma_H_he2

    implicit none

    ! The time step
    real(kind=dp), intent(in) :: dt
    ! mesh position of cell 
    integer,intent(in) :: pos
    ! Indicator if too less photon reach this cell
    logical,intent(inout) :: too_less_photon
    ! Column density for stopping chemisty
    real(kind=dp), parameter :: max_coldens_HI = 2.0e26 
    ! Counter of number of iterations before convergence
    integer :: i_iteration
    ! Incoming column densities
    real(kind=dp) :: coldens_in_HI
    real(kind=dp) :: coldens_in_HeI
    real(kind=dp) :: coldens_in_HeII
    ! Intermidiate column densities of the cell during iterations
    real(kind=dp) :: intermediate_coldens_cell_HI
    real(kind=dp) :: intermediate_coldens_cell_HeI
    real(kind=dp) :: intermediate_coldens_cell_HeII
    ! Intermidiate outgoing column densities during iterations
    real(kind=dp) :: intermediate_coldens_out_HI
    real(kind=dp) :: intermediate_coldens_out_HeI
    real(kind=dp) :: intermediate_coldens_out_HeII
    ! Intermediate averaged ionized fractions
    real(kind=dp) :: intermediate_ion_avg_HI
    real(kind=dp) :: intermediate_ion_avg_HeI
    real(kind=dp) :: intermediate_ion_avg_HeII
    ! Intermediate end time ionized fractions
    real(kind=dp) :: intermediate_ion_end_HI
    real(kind=dp) :: intermediate_ion_end_HII
    real(kind=dp) :: intermediate_ion_end_HeI
    real(kind=dp) :: intermediate_ion_end_HeII
    real(kind=dp) :: intermediate_ion_end_HeIII
    ! Ionized fractions and temperature of the previous run, be used in convergence test
    real(kind=dp) :: prev_temper
    real(kind=dp) :: prev_avg_HI
    real(kind=dp) :: prev_avg_HII
    real(kind=dp) :: prev_avg_HeI
    real(kind=dp) :: prev_avg_HeII
    real(kind=dp) :: prev_avg_HeIII
    ! Ionized fraction of the cell
    type(ionstates) :: ion
    ! Temperature of the cell
    real(kind=dp) :: begin_temper
    real(kind=dp) :: end_temper
    real(kind=dp) :: avg_temper
    ! Photoionization rate and heating rate
    type(photrates) :: phi
    ! Size of a cell
    real(kind=dp) :: path
    ! Volume of the cell
    real(kind=dp) :: vol_cell
    ! Electron density of the cell
    real(kind=dp) :: ndens_electron
    ! Number density of atoms of the cell
    real(kind=dp) :: ndens_atom
    ! Parameters to solve coupled ODEs
    real(kind=dp) :: tau_H_heth, tau_H_heLya, tau_He_heth, tau_He_heLya, y, z
    real(kind=dp) :: y2a, y2b, tau_He2_he2th, tau_He_he2th, tau_H_he2th
    !logical, parameter :: do_artificial_temperature = .true.
    logical, parameter :: do_artificial_temperature = .false.
	
    ! copy the ionization states from the previous time-step
    ion%begin_HI = xHI(pos)
    ion%avg_HI = xHI(pos)
    ion%end_HI = xHI(pos)
    ion%begin_HII = xHII(pos)
    ion%avg_HII = xHII(pos)
    ion%end_HII = xHII(pos)
    ion%begin_HeI = xHeI(pos)
    ion%avg_HeI = xHeI(pos)
    ion%end_HeI = xHeI(pos)
    ion%begin_HeII = xHeII(pos)
    ion%avg_HeII = xHeII(pos)
    ion%end_HeII = xHeII(pos)
    ion%begin_HeIII = xHeIII(pos)
    ion%avg_HeIII = xHeIII(pos)
    ion%end_HeIII = xHeIII(pos)

    ! copy the temperature of the previous time-step
    begin_temper = temper(pos)
    avg_temper = temper(pos)
    end_temper = temper(pos)

    ! copy number density of the previous time-step
    ndens_atom=ndens(pos,1,1)

    ! copy outgoing column density of previous cell 
    ! as the incoming column density of current cell
    if (pos.eq.1) then
       coldens_in_HI = coldens_bndry_HI()
       coldens_in_HeI = coldens_bndry_HeI()
       coldens_in_HeII = coldens_bndry_HeII()
    else
       coldens_in_HI = coldens_out_HI(pos-1)
       coldens_in_HeI = coldens_out_HeI(pos-1)
       coldens_in_HeII = coldens_out_HeII(pos-1)
    endif

    ! Cell size
    path = dr(1)

    ! Find the volume of the shell this cell is part of (dilution factor)
    vol_cell = vol(pos)

    ! Initialize iteration counter to zero
    i_iteration = 0

    if (coldens_in_HI .le. max_coldens_HI) then

      do

          i_iteration = i_iteration+1

          ! Save the value of the previous iteration
          prev_avg_HI = ion%avg_HI
          prev_avg_HII = ion%avg_HII
          prev_avg_HeI = ion%avg_HeI
          prev_avg_HeII = ion%avg_HeII
          prev_avg_HeIII = ion%avg_HeIII
          prev_temper = end_temper

          ! Update current cell average HI, HeI and HeII column densities 
          intermediate_coldens_cell_HI = coldens(path,ion%avg_HI,ndens_atom,(1.0_dp-abu_he))
          intermediate_coldens_cell_HeI = coldens(path,ion%avg_HeI,ndens_atom,abu_he)
          intermediate_coldens_cell_HeII = coldens(path,ion%avg_HeII,ndens_atom,abu_he)

          ! Update outgoing average column density of HI, HeI and HeII
          intermediate_coldens_out_HI = coldens_in_HI+intermediate_coldens_cell_HI
          intermediate_coldens_out_HeI = coldens_in_HeI+intermediate_coldens_cell_HeI  
          intermediate_coldens_out_HeII = coldens_in_HeII+intermediate_coldens_cell_HeII     

          ! Obtain photoionizing rate and heating rate on the current cell
          call photoion(phi,coldens_in_HI,intermediate_coldens_out_HI,coldens_in_HeI, &
                        intermediate_coldens_out_HeI,coldens_in_HeII,intermediate_coldens_out_HeII, &
                        vol_cell,1,ion%avg_HII)

          ! Check if the photon rate is too low
          if (phi%photo_in/source_ionzing_photon_rate .le. minium_fraction_of_photons) then
            too_less_photon = .true.
          endif

		  
			  avg_ex_ph_eV(pos) = phi%photo_cell_HI+phi%photo_cell_HeI+phi%photo_cell_HeII
			  avg_ex_ph_eV(pos) =  phi%heat/avg_ex_ph_eV(pos)
			  avg_ex_ph_eV(pos) = avg_ex_ph_eV(pos)/ev2erg ! change to eV
			  !write(*,*) 'averaged excess photon enegy in eV is ', avg_ex_ph_eV
			  

          ! True photoionization rate
          phi%photo_cell_HI = phi%photo_cell_HI/(ndens_atom*(1.0_dp-abu_he)*ion%avg_HI)
          phi%photo_cell_HeI = phi%photo_cell_HeI/(ndens_atom*abu_he*ion%avg_HeI)
          phi%photo_cell_HeII = phi%photo_cell_HeII/(ndens_atom*abu_he*ion%avg_HeII)

          ! Add the UV background
          phi%photo_cell_HI = phi%photo_cell_HI + gamma_uvb(1)
          phi%photo_cell_HeI = phi%photo_cell_HeI + gamma_uvb(2)
          phi%photo_cell_HeII = phi%photo_cell_HeII + gamma_uvb(3)

          ! Update current cell average electron density
          ndens_electron = electrondens(ndens_atom,ion%avg_HII,ion%avg_HeII,ion%avg_HeIII)

          ! Update the collisional ion and recombination rates for new temperature
          if (.not.isothermal) call ini_rec_colion_factors(avg_temper) 

          ! Update current cell avergae column density of HI, HeI and HeII
          intermediate_coldens_cell_HI = coldens(path,ion%avg_HI, ndens_atom,(1.0_dp-abu_he)) 
          intermediate_coldens_cell_HeI = coldens(path,ion%avg_HeI,ndens_atom, abu_he)
          intermediate_coldens_cell_HeII = coldens(path,ion%avg_HeII,ndens_atom,abu_he)

          ! These are the specific optical depth in order to find out y,z,y2a,y2b
          tau_H_heth = intermediate_coldens_cell_HI*sigma_H_heth 
          tau_He_heth = intermediate_coldens_cell_HeI*sigma_HeI_at_ion_freq 
          tau_H_heLya = intermediate_coldens_cell_HI*sigma_H_heLya 
          tau_He_heLya = intermediate_coldens_cell_HeI*sigma_He_heLya 
          tau_He2_he2th = intermediate_coldens_cell_HeII*sigma_HeII_at_ion_freq
          tau_He_he2th = intermediate_coldens_cell_HeI* sigma_He_he2
          tau_H_he2th = intermediate_coldens_cell_HI*sigma_H_he2

          ! Parameters required to define the coupled ODEs
          y = tau_H_heth /(tau_H_heth +tau_He_heth)
          z = tau_H_heLya/(tau_H_heLya+tau_He_heLya)
          y2a = tau_He2_he2th /(tau_He2_he2th +tau_He_he2th+tau_H_he2th)
          y2b = tau_He_he2th /(tau_He2_he2th +tau_He_he2th+tau_H_he2th)

          ! Calculate the analytical solution of ionized fractions
          call doric(dt,ndens_electron,ion,phi,y,z,y2a,y2b,alpha_HI_B,pos,avg_temper)! 

          ! Update current cell average electron density
          ndens_electron = electrondens(ndens_atom,ion%avg_HII,ion%avg_HeII,ion%avg_HeIII)

          ! Update current cell avergae column density of HI, HeI and HeII
          intermediate_coldens_cell_HI = coldens(path,ion%avg_HI, ndens_atom,(1.0_dp-abu_he)) 
          intermediate_coldens_cell_HeI = coldens(path,ion%avg_HeI,ndens_atom, abu_he)
          intermediate_coldens_cell_HeII = coldens(path,ion%avg_HeII,ndens_atom,abu_he)
 
          ! These are the specific optical depth in order to find out y,z,y2a,y2b
          tau_H_heth = intermediate_coldens_cell_HI *sigma_H_heth
          tau_He_heth = intermediate_coldens_cell_HeI*sigma_HeI_at_ion_freq  
          tau_H_heLya = intermediate_coldens_cell_HI*sigma_H_heLya 
          tau_He_heLya= intermediate_coldens_cell_HeI*sigma_He_heLya 
          tau_He2_he2th = intermediate_coldens_cell_HeII*sigma_HeII_at_ion_freq
          tau_He_he2th = intermediate_coldens_cell_HeI* sigma_He_he2
          tau_H_he2th = intermediate_coldens_cell_HI*sigma_H_he2

          ! Parameters required to define the coupled ODEs
          y = tau_H_heth /(tau_H_heth +tau_He_heth)
          z = tau_H_heLya/(tau_H_heLya+tau_He_heLya)
          y2a = tau_He2_he2th /(tau_He2_he2th +tau_He_he2th+tau_H_he2th)
          y2b = tau_He_he2th /(tau_He2_he2th +tau_He_he2th+tau_H_he2th)

          ! Save the ionization results of the previous doric    
          intermediate_ion_end_HI = ion%end_HI
          intermediate_ion_end_HII = ion%end_HII
          intermediate_ion_end_HeI = ion%end_HeI
          intermediate_ion_end_HeII = ion%end_HeII
          intermediate_ion_end_HeIII = ion%end_HeIII
          intermediate_ion_avg_HI = ion%avg_HI
          intermediate_ion_avg_HeI = ion%avg_HeI
          intermediate_ion_avg_HeII = ion%avg_HeII

          ! Calculate the analytical solution of ionized fractions
          call doric(dt,ndens_electron,ion,phi, y, z,y2a,y2b,alpha_HI_B,pos,avg_temper)

          ! Average the ionization results of two dorics.
          ion%end_HI = 0.5*(ion%end_HI+intermediate_ion_end_HI)
          ion%end_HII = 0.5*(ion%end_HII+intermediate_ion_end_HII)
          ion%end_HeI = 0.5*(ion%end_HeI+intermediate_ion_end_HeI)
          ion%end_HeII = 0.5*(ion%end_HeII+intermediate_ion_end_HeII)
          ion%end_HeIII = 0.5*(ion%end_HeIII+intermediate_ion_end_HeIII)
          ion%avg_HI = 0.5*(ion%avg_HI+intermediate_ion_avg_HI)
          ion%avg_HeI = 0.5*(ion%avg_HeI+intermediate_ion_avg_HeI)
          ion%avg_HeII = 0.5*(ion%avg_HeII+intermediate_ion_avg_HeII)

          ! Update current cell average electron density
          ndens_electron = electrondens(ndens_atom,ion%avg_HII,ion%avg_HeII,ion%avg_HeIII)

          ! Reset the final temperature to the beginning temperature
          end_temper = begin_temper 

          ! Calculate the temperature of the cell
          if (.not.isothermal) call thermal(dt,end_temper,avg_temper,ndens_electron,ndens_atom,ion,phi,pos)


!!!!!!!!!! for the artificial thermal profile !!!!!!!!!!!!
          if (do_artificial_temperature) then
            if (pos.le.84) then
              end_temper = 40000
              avg_temper = 40000
            else
              	end_temper = 10000
              	avg_temper = 10000
            endif		
          endif
!!!!!!!!!! for the artificial thermal profile !!!!!!!!!!!!


!if (isnan(end_temper) .and. report_counter.eq.-1) then
!    write(*,*)'at position ',pos
!    write(*,*) 'avg xHI: ',ion%avg_HI,prev_avg_HI
!    write(*,*) 'avg xHII: ',ion%avg_HII,prev_avg_HII
!    write(*,*) 'avg xHeI: ',ion%avg_HeI, prev_avg_HeI
!    write(*,*) 'avg xHeII: ',ion%avg_HeII, prev_avg_HeII
!    write(*,*) 'avg xHeIII: ',ion%avg_HeIII, prev_avg_HeIII
!    write(*,*) 'end xHI: ',ion%end_HI
!    write(*,*) 'end xHII: ',ion%end_HII
!    write(*,*) 'end xHeI: ',ion%end_HeI
!    write(*,*) 'end xHeII: ',ion%end_HeII
!    write(*,*) 'end xHeIII: ',ion%end_HeIII
!    write(*,*) 'temper: ',end_temper,prev_temper
!    report_counter = report_counter+1
!endif

            ! Convergence test
            if( ( (abs(ion%avg_HI-prev_avg_HI)/ion%avg_HI .lt. minimum_fractional_change) &
               .or. &
               (ion%avg_HI .lt. minimum_fraction_of_atoms) ) &
               .and. &
               ( (abs(ion%avg_HeI-prev_avg_HeI)/ion%avg_HeI .lt. minimum_fractional_change) &
               .or. &
               (ion%avg_HeI .lt. minimum_fraction_of_atoms) ) &
               .and. &
               ( (abs(ion%avg_HeIII-prev_avg_HeIII)/ion%avg_HeIII .lt. minimum_fractional_change) &
               .or. &
               (ion%avg_HeIII .lt. minimum_fraction_of_atoms) ) &   
               .and. &
               abs(end_temper-prev_temper)/end_temper .lt. minimum_fractional_change ) then
               HI_photoionization_rate(pos) = phi%photo_cell_HI
			   recombination_time(pos) = 1.0/(alpha_HI_B*ion%avg_HII*ndens(pos,1,1)*abu_h*YEAR)
			   recombination_rate(pos) = alpha_HI_B
			   if (isnan(end_temper)) write(*,*)'at position ',pos
               exit

            endif

            ! when after 4000 iterations still not converge..
            if (i_iteration.gt.500) then
               write(logf,*) 'convcrits of HI ',abs(ion%avg_HI-prev_avg_HI)/ion%avg_HI,ion%avg_HI 
               write(logf,*) 'convcrits of HeI ', abs(ion%avg_HeI-prev_avg_HeI)/ion%avg_HeI,ion%avg_HeI 
               write(logf,*) 'convcrits of HeII ',abs(ion%avg_HeII-prev_avg_HeII)/ion%avg_HeII, ion%avg_HeII
               write(logf,*) 'convcrits of HeIII ',abs(ion%avg_HeIII-prev_avg_HeIII)/ion%avg_HeIII, ion%avg_HeIII
               write(logf,*) 'convcrits of T ', abs(end_temper-prev_temper)/end_temper, end_temper
               write(33,*) 'over 4000 iterations...'
               HI_photoionization_rate(pos) = phi%photo_cell_HI
               recombination_time(pos) = 1.0/(alpha_HI_B*ion%avg_HII*ndens(pos,1,1)*abu_h*YEAR)
               exit
            endif

      enddo ! end of iteration

    ! when incoming optical depth is too big that no ionizing photons can reach the cell
    else

       phi%photo_cell_HI = 0.0_dp
       phi%photo_out_HI = 0.0_dp
       phi%photo_cell_HeI = 0.0_dp
       phi%photo_cell_HeII = 0.0_dp
       phi%photo_out_HeI = 0.0_dp
       phi%photo_out_HeII = 0.0_dp

    endif

    !Copy the final ionization results back
    xHI(pos) = ion%end_HI
    xHII(pos) = ion%end_HII
    xHeI(pos) = ion%end_HeI    
    xHeII(pos) = ion%end_HeII  
    xHeIII(pos) = ion%end_HeIII  

    ! Copy the final temperature back
    temper(pos) = end_temper

    ! Copy the averaged outgoing column densities back
    coldens_out_HI(pos)  =coldens_in_HI+ &
                         coldens(path,ion%avg_HI,ndens_atom,(1.0_dp-abu_he))
    coldens_out_HeI(pos) =coldens_in_HeI+ &
                         coldens(path,ion%avg_HeI,ndens_atom,abu_he)
    coldens_out_HeII(pos)=coldens_in_HeII+ &
                         coldens(path,ion%avg_HeII,ndens_atom,abu_he)

    ! Photon statistics: register number of photons leaving the grid
    if (pos.eq.mesh) photon_loss = 0.0_dp

  end subroutine evolve0D

end module evolve
