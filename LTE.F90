! This module checks if a particular cell is in local thermal equilibrium
! until the next output time

module LTE

  use precision, only: dp
  use grid, only: vol, dr
  use material, only: ndens, xHI, xHII, xHeI, xHeII, xHeIII, temper, HI_photoionization_rate,&
			   recombination_time
  use abundances, only: abu_he
  use tped, only: electrondens
  use doric_module, only: doric, coldens, &
                          coldens_bndry_HI, coldens_bndry_HeI, coldens_bndry_HeII
  use radiation, only: photoion, photrates
  use c2ray_parameters, only: minimum_fraction_of_atoms, LTE_fractional_change
  use material, only: isothermal, ionstates, gamma_uvb, avg_ex_ph_eV
  use thermalevolution, only: thermal
  use cgsconstants,only : ini_rec_colion_factors
  use cgsphotoconstants, only: sigma_H_heth, sigma_H_heLya, sigma_He_heLya, sigma_HeI_at_ion_freq, &
                               sigma_He_he2, sigma_HeII_at_ion_freq, sigma_H_he2
  use evolve, only: coldens_out_HI, coldens_out_HeI, coldens_out_HeII
  use abundances, only: abu_he,abu_h
    use astroconstants, only: YEAR
  use cgsconstants, only: ev2erg	
	
  implicit none

contains
   
  subroutine LTE_calculation (inside_pos,front_pos,B_dt)

    ! furthest position of cell in LTE from the source 
    integer, intent(inout) :: inside_pos
    ! I-front position
    integer, intent(in) :: front_pos
    ! time elapsed until the next output time
    real(kind=dp), intent(in) :: B_dt
    ! position of cell in  LTE 
    integer :: maybe_inside_pos
    ! logical which indicates if LTE is achieved
    logical :: B_convergence
    ! position of cell being investigated
    integer :: pos
    ! counter of iteration
    integer :: B_iteration

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
    ! Size of a cell
    real(kind=dp) :: path
    ! Volume of the cell
    real(kind=dp) :: vol_cell
    ! Electron density of the cell
    real(kind=dp) :: ndens_electron
    ! Number density of atoms of the cell
    real(kind=dp) :: ndens_atom
    real(kind=dp) :: intermediate_ion_avg_HI
    real(kind=dp) :: intermediate_ion_avg_HeI
    real(kind=dp) :: intermediate_ion_avg_HeII
    ! Parameters to solve coupled ODEs
    real(kind=dp) :: tau_H_heth, tau_H_heLya, tau_He_heth, tau_He_heLya, y, z
    real(kind=dp) :: y2a, y2b, tau_He2_he2th, tau_He_he2th, tau_H_he2th
    ! Temperature of the cell
    real(kind=dp) :: begin_temper
    real(kind=dp) :: end_temper
    real(kind=dp) :: avg_temper
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
    type(ionstates) :: ion
    type(photrates) :: phi
	real(kind=dp) :: alpha_HI_B
		
    ! Initialize the LTE logical to true
    B_convergence = .true.
    ! This cell may be in LTE, check now..
    maybe_inside_pos = inside_pos+1

    ! The LTE check on a cell is performed if this cell is behind from the I-front
    ! and on the previous cell it is in LTE
    do while (B_convergence.eqv..true. .and. maybe_inside_pos.lt.front_pos)

    ! copy the ionization states of the cell
    ion%begin_HI = xHI(maybe_inside_pos)
    ion%avg_HI = xHI(maybe_inside_pos)
    ion%end_HI = xHI(maybe_inside_pos)
    ion%begin_HII = xHII(maybe_inside_pos)
    ion%avg_HII = xHII(maybe_inside_pos)
    ion%end_HII = xHII(maybe_inside_pos)
    ion%begin_HeI = xHeI(maybe_inside_pos)
    ion%avg_HeI = xHeI(maybe_inside_pos)
    ion%end_HeI = xHeI(maybe_inside_pos)
    ion%begin_HeII = xHeII(maybe_inside_pos)
    ion%avg_HeII = xHeII(maybe_inside_pos)
    ion%end_HeII = xHeII(maybe_inside_pos)
    ion%begin_HeIII = xHeIII(maybe_inside_pos)
    ion%avg_HeIII = xHeIII(maybe_inside_pos)
    ion%end_HeIII = xHeIII(maybe_inside_pos)

    ! copy the temperature of the cell
    begin_temper = temper(maybe_inside_pos)
    avg_temper = temper(maybe_inside_pos)
    end_temper = temper(maybe_inside_pos)

    !copy number density of the cell
    ndens_atom = ndens(maybe_inside_pos,1,1)

    ! copy outgoing column density of previous cell 
    ! as the incoming column density of current cell
    if (pos.eq.1) then
       coldens_in_HI = coldens_bndry_HI()
       coldens_in_HeI = coldens_bndry_HeI()
       coldens_in_HeII = coldens_bndry_HeII()
    else
       coldens_in_HI = coldens_out_HI(maybe_inside_pos-1)
       coldens_in_HeI = coldens_out_HeI(maybe_inside_pos-1)
       coldens_in_HeII = coldens_out_HeII(maybe_inside_pos-1)
    endif

    ! Cell size
    path = dr(1)

    ! Find the volume of the shell this cell is part of (dilution factor)
    vol_cell = vol(maybe_inside_pos)

    ! Initialize iteration counter to zero
    B_iteration = 0

    do

      B_iteration = B_iteration+1

      ! Save the value of the previous iteration
      prev_avg_HI = ion%avg_HI
      prev_avg_HII = ion%avg_HII
      prev_avg_HeI = ion%avg_HeI
      prev_avg_HeII = ion%avg_HeII
      prev_avg_HeIII = ion%avg_HeIII
      prev_temper = end_temper

      ! Update current cell average column density of HI, HeI and HeII
      intermediate_coldens_cell_HI = coldens(path,ion%avg_HI,ndens_atom,(1.0_dp-abu_he))
      intermediate_coldens_cell_HeI = coldens(path,ion%avg_HeI,ndens_atom,abu_he)
      intermediate_coldens_cell_HeII = coldens(path,ion%avg_HeII,ndens_atom,abu_he)

      ! Update outgoing average column density of HI, HeI and HeII
      intermediate_coldens_out_HI = coldens_in_HI+intermediate_coldens_cell_HI
      intermediate_coldens_out_HeI = coldens_in_HeI+intermediate_coldens_cell_HeI  
      intermediate_coldens_out_HeII = coldens_in_HeII+intermediate_coldens_cell_HeII     

      ! Obtain photoionizing rate and heating rate on the current cell
      call photoion(phi,coldens_in_HI,intermediate_coldens_out_HI,coldens_in_HeI, &
                    intermediate_coldens_out_HeI,coldens_in_HeII, intermediate_coldens_out_HeII, &
                    vol_cell,1,ion%avg_HII)

avg_ex_ph_eV(maybe_inside_pos) = phi%photo_cell_HI+phi%photo_cell_HeI+phi%photo_cell_HeII
avg_ex_ph_eV(maybe_inside_pos) =  phi%heat/avg_ex_ph_eV(maybe_inside_pos)
avg_ex_ph_eV(maybe_inside_pos) = avg_ex_ph_eV(maybe_inside_pos)/ev2erg ! change to eV
				  
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

      !Update the collisional ion and recombination rates for new temperature
      if (.not.isothermal) call ini_rec_colion_factors(avg_temper) 

      ! Update current cell avergae column density of HI, HeI and HeII
      intermediate_coldens_cell_HI = coldens(path,ion%avg_HI,ndens_atom,(1.0_dp-abu_he)) 
      intermediate_coldens_cell_HeI = coldens(path,ion%avg_HeI,ndens_atom, abu_he)
      intermediate_coldens_cell_HeII = coldens(path,ion%avg_HeII,ndens_atom,abu_he)

      !These are the specific optical depth in order to find out y,z,y2a,y2b
      tau_H_heth = intermediate_coldens_cell_HI*sigma_H_heth 
      tau_He_heth = intermediate_coldens_cell_HeI*sigma_HeI_at_ion_freq 
      tau_H_heLya = intermediate_coldens_cell_HI*sigma_H_heLya 
      tau_He_heLya= intermediate_coldens_cell_HeI*sigma_He_heLya 
      tau_He2_he2th = intermediate_coldens_cell_HeII*sigma_HeII_at_ion_freq
      tau_He_he2th = intermediate_coldens_cell_HeI*sigma_He_he2
      tau_H_he2th = intermediate_coldens_cell_HI*sigma_H_he2

      ! Parameters required to define the coupled ODEs
      y = tau_H_heth /(tau_H_heth +tau_He_heth)
      z = tau_H_heLya/(tau_H_heLya+tau_He_heLya)
      y2a = tau_He2_he2th /(tau_He2_he2th +tau_He_he2th+tau_H_he2th)
      y2b = tau_He_he2th /(tau_He2_he2th +tau_He_he2th+tau_H_he2th)

      ! Calculate the analytical solution of ionized fractions
      call doric(B_dt,ndens_electron,ion,phi,y,z,y2a,y2b,alpha_HI_B,maybe_inside_pos,end_temper)! 

      ! Update current cell average electron density
      ndens_electron = electrondens(ndens_atom,ion%avg_HII,ion%avg_HeII,ion%avg_HeIII)

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
      tau_He_he2th = intermediate_coldens_cell_HeI*sigma_He_he2
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
      call doric(B_dt,ndens_electron,ion,phi, y, z,y2a,y2b,alpha_HI_B,maybe_inside_pos,end_temper)

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
      if (.not.isothermal) call thermal(B_dt,end_temper,avg_temper,ndens_electron,ndens_atom,ion,phi,pos)

      ! Convergence test
      if ( ( (abs(ion%avg_HI-prev_avg_HI)/ion%avg_HI .lt. LTE_fractional_change) &
             .or. &
             (ion%avg_HI.lt.minimum_fraction_of_atoms) &
             ) &
             .and. &
             ( (abs(ion%avg_HeI-prev_avg_HeI)/ion%avg_HeI .lt. LTE_fractional_change) &
             .or. &
             (ion%avg_HeI.lt.minimum_fraction_of_atoms) &
             ) &
             .and. &
             ( (abs(ion%avg_HeIII-prev_avg_HeIII)/ion%avg_HeIII .lt. LTE_fractional_change) &
             .or. &
             (ion%avg_HeIII.lt.minimum_fraction_of_atoms) &
             ) &   
             .and. &
             abs(end_temper-prev_temper)/end_temper .lt. LTE_fractional_change &
             ) then
             HI_photoionization_rate(maybe_inside_pos) = phi%photo_cell_HI
		     recombination_time(maybe_inside_pos) = &
			 1.0/(alpha_HI_B*ion%avg_HI*ndens(maybe_inside_pos,1,1)*abu_h*YEAR)
            exit

       endif

     enddo

     ! LTE test
     if  ( (ion%avg_HI .lt. LTE_fractional_change) &
           .and. &
           (ion%avg_HeI .lt. minimum_fraction_of_atoms) &
           .and. &
           ((abs(ion%avg_HeIII-xHeIII(maybe_inside_pos))/xHeIII(maybe_inside_pos)) .lt. LTE_fractional_change) &
           .and. &
           ((abs(end_temper-temper(maybe_inside_pos))/temper(maybe_inside_pos)) .lt. LTE_fractional_change) &
         ) then

         maybe_inside_pos = maybe_inside_pos+1

      else

         B_convergence = .false.

      endif

      
    enddo

    ! new assignment of maybe inside position
    inside_pos = maybe_inside_pos-1

  end subroutine LTE_calculation

end module LTE

