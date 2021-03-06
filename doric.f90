  ! This module contains routines having to do with the calculation of
  ! the ionization evolution of a single grid point.

module doric_module

  use precision, only: dp

  implicit none

  integer, save :: report_counter = 0
  
  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  ! calculates time dependent ionization state for hydrogen and helium
  subroutine doric (dt,ndens_electron,ion,phi,y,z,y2a,y2b,alpha_HI_B,pos,T)
    
    use cgsconstants, only: bh00, bhe00, bhe10, albpow, alcpow, colh0, colhe, temphe, temph0, &
                            ev2fr,brech0,areche0,breche0, oreche0,breche1, treche1, n_el_crit, &
                            colli_HI, colli_HeI, colli_HeII, ini_rec_colion_factors, &
                            arech0, areche1, v
    use material, only: clumping, ionstates
    use abundances, only: abu_he
    use radiation, only: photrates
    use c2ray_parameters, only: epsilon

    ! time step
    real(kind=dp), intent(in) :: dt     
    ! electron density
    real(kind=dp), intent(inout) :: ndens_electron 
    ! ionzation states
    type(ionstates) :: ion 
    ! photo-ionization rates
    type(photrates), intent(in)  :: phi    
    real(kind=dp), intent(inout) :: y
    real(kind=dp), intent(inout) :: z   
    real(kind=dp), intent(inout) :: y2a
    real(kind=dp), intent(inout) :: y2b
    real(kind=dp), intent(out) :: alpha_HI_B
	integer, intent(in) :: pos
    real(kind=dp), intent(in) :: T
		
    real(kind=dp) :: alpha_HI_A
    real(kind=dp) :: alpha_HeI_1, alpha_HeI_B, alpha_HeI_A
    real(kind=dp) :: alpha_HeII_B, alpha_HeII_A, alpha_HeII_2, alpha_HeII_1 
    real(kind=dp) :: p, f, w, l, m
    real(kind=dp) :: lambda_1, lambda_2, lambda_3
    real(kind=dp) :: lambda_1_dt, lambda_2_dt, lambda_3_dt
    real(kind=dp) :: exp_lambda_1_dt, exp_lambda_2_dt, exp_lambda_3_dt
    real(kind=dp) :: x2_1, x2_2, x3_1, x3_2
    real(kind=dp) :: c_1, c_2, c_3
    real(kind=dp) :: Scoef, Rcoef, Tcoef, Kcoef 
    real(kind=dp) :: Bcoef, Bminus, Bplus
    real(kind=dp) :: avg_factor_1, avg_factor_2, avg_factor_3
    real(kind=dp) :: A_11, A_12, A_13, A_22, A_23, A_32, A_33, two_A_32
    real(kind=dp) :: xp_1, xp_2, xp_3
    real(kind=dp) :: g_1, g_2
    real(kind=dp) :: heliumfraction, normfac
    !logical, parameter :: caseB = .true.
    logical, parameter :: caseB = .false.
	
    logical, parameter :: do_coupled = .true.
	!logical, parameter :: do_coupled = .false.
    
	!if (pos.le.414) then
	!	clumping = 2.0
	!else
	!	clumping = 2.0
    !endif
		
	!clumping = 3-(T-10000.0)/10000.0
	!clumping = 2
	
    ! see Osterbrock, 1989
    p = 0.96_dp 
    ! fraction of photons from 2-photon decay, energetic enough to ionize hydrogen
    l = 1.425_dp 
    ! fraction of photons from 2-photon decay, energetic enough to ionize neutral helium
    m = 0.737_dp 
    heliumfraction = abu_he/(1.0_dp-abu_he)    
    f = max(min(10.0_dp*ion%end_HI,1.0_dp),0.01_dp) 
    ! These numbers come from Flower & Perinotto (1980)
    w = (l-m)+m*y  

    ! absorption coefficients
    alpha_HI_B = clumping*brech0
    alpha_HI_A = clumping*arech0
    alpha_HeI_1 = clumping*oreche0
    alpha_HeI_B = clumping*breche0    
    alpha_HeI_A = clumping*areche0
    alpha_HeII_B = clumping*breche1
    alpha_HeII_A = clumping*areche1
    alpha_HeII_2 = clumping*treche1
    alpha_HeII_1 = alpha_HeII_A-alpha_HeII_B

    ! elements of matrix A and vector g
    !g_1 = max(phi%photo_cell_HI+ndens_electron*colli_HI,1.0e-200_dp)
    !g_2 = max(phi%photo_cell_HeI+ndens_electron*colli_HeI,1.0e-200_dp)
    !A_32 = max(phi%photo_cell_HeII+ndens_electron*colli_HeII,1.0e-200_dp)
    !A_11 = -(g_1+ndens_electron*alpha_HI_B)
	

    !A_12 = (y*alpha_HeI_1+p*alpha_HeI_B)*ndens_electron*heliumfraction
    !A_13 = ((f*z*(1.0_dp-v) +v*w)*alpha_HeII_B +alpha_HeII_2+ &
    !       (1.0_dp-y2a-y2b)*alpha_HeII_1)*heliumfraction*ndens_electron  
    !A_22 = -g_2-A_32-ndens_electron*(alpha_HeI_A-(1.0_dp-y)*alpha_HeI_1)
    !A_33 = -ndens_electron*(alpha_HeII_A-y2a*alpha_HeII_1)            
    !A_23 = -g_2+ndens_electron*alpha_HeII_B*(f*(1.0_dp-z)*(1.0_dp-v)+ &
    !       v*(1.425_dp-w))-A_33+alpha_HeII_1*y2b*ndens_electron

 if (do_coupled) then

     ! elements of matrix A and vector g
     g_1 = max(phi%photo_cell_HI+ndens_electron*colli_HI,1.0e-200_dp)
     g_2 = max(phi%photo_cell_HeI+ndens_electron*colli_HeI,1.0e-200_dp)
     A_32 = max(phi%photo_cell_HeII+ndens_electron*colli_HeII,1.0e-200_dp)
	 
	 if (caseB) then
	 
            A_11 = -(g_1+ndens_electron*alpha_HI_B)
		    A_12 = (y*alpha_HeI_1+p*alpha_HeI_B)*ndens_electron*heliumfraction
		    A_13 = ((f*z*(1.0_dp-v) +v*w)*alpha_HeII_B +alpha_HeII_2+ &
		           (1.0_dp-y2a-y2b)*alpha_HeII_1)*heliumfraction*ndens_electron
		    A_22 = -g_2-A_32-ndens_electron*(alpha_HeI_A-(1.0_dp-y)*alpha_HeI_1)
		    A_33 = -ndens_electron*(alpha_HeII_A-y2a*alpha_HeII_1)
		    A_23 = -g_2+ndens_electron*alpha_HeII_B*(f*(1.0_dp-z)*(1.0_dp-v)+ &
		           v*(1.425_dp-w))-A_33+alpha_HeII_1*y2b*ndens_electron

	 else
		 
            A_11 = -(g_1+ndens_electron*alpha_HI_A)
	        A_12 = 0.0_dp
	        A_13 = 0.0_dp
	        A_22 = -g_2-A_32-ndens_electron*alpha_HeI_A
	        A_33 = -ndens_electron*alpha_HeII_A
	        A_23 = -g_2-A_33	 
			!if (pos.eq.64) write(*,*) 'g_1', g_1
			!if (pos.eq.64) write(*,*) 'alpha_HI_A', alpha_HI_A
			!if (pos.eq.64) write(*,*) 'A_11', A_11
			!if (pos.eq.64) write(*,*) 'A_13', A_13
			!if (pos.eq.64) write(*,*) 'A_22', A_22
			!if (pos.eq.64) write(*,*) 'A_33', A_33
			!if (pos.eq.64) write(*,*) 'A_23', A_23
	 endif
	 
  else

      ! elements of matrix A and vector g
      g_1 = max(phi%photo_cell_HI,1.0e-200_dp)
      g_2 = max(phi%photo_cell_HeI,1.0e-200_dp)
      A_32 = max(phi%photo_cell_HeII,1.0e-200_dp)
	  
	 if (caseB) then
		 	  
            A_11 = -(g_1+ndens_electron*alpha_HI_B)
		    A_12 = 0.0_dp
		    A_13 = 0.0_dp
		    A_22 = -(g_2+A_32+ndens_electron*alpha_HeI_B)
		    A_33 = -ndens_electron*alpha_HeII_B
		    A_23 = -g_2-A_33

	 else
		 
            A_11 = -(g_1+ndens_electron*alpha_HI_A)
	        A_12 = 0.0_dp
	        A_13 = 0.0_dp
	        A_22 = -(g_2+A_32+ndens_electron*alpha_HeI_A)
	        A_33 = -ndens_electron*alpha_HeII_A
	        A_23 = -g_2-A_33
		
	 endif
	 		 
  endif
			
    ! this is just abbrivations
    two_A_32 = 2.0_dp*A_32

    ! These are just abbrivations 
    Bcoef = A_33-A_22
    Scoef = sqrt(Bcoef*Bcoef+4.0_dp*A_32*A_23)
    Kcoef = 1.0_dp/(A_23*A_32-A_33*A_22)
    Bminus = Bcoef-Scoef
    Bplus = Bcoef+Scoef

    ! eigenvalues of the matrix A
    lambda_1 = A_11
    lambda_2 = 0.5_dp*(A_33+A_22-Scoef)
    lambda_3 = 0.5_dp*(A_33+A_22+Scoef)

    ! elements of particular solution vector
    xp_1 = -1.0_dp/A_11*(g_1+(A_12*A_33-A_13*A_32)*(g_2*Kcoef))   
    xp_2 = g_2*(A_33*Kcoef)                        
    xp_3 = -g_2*(A_32*Kcoef)                      
 
    ! elements of some column vectors
    x2_1 = -A_13/(A_11-lambda_2)+(A_12/two_A_32)*Bplus/(A_11-lambda_2)  
    x2_2 = (-Bplus)/(two_A_32)
    x3_1 = (-two_A_32 *A_13+A_12*(Bminus))/(two_A_32*(A_11-lambda_3))  
    x3_2 = (-Bminus)/(two_A_32)

    ! These are just abbrivations
    Rcoef = two_A_32 *(xp_2-ion%begin_HeII)
    Tcoef = xp_3-ion%begin_HeIII

    ! coefficients of homogeneous solutions
    c_1 = -xp_1+ (x3_1-x2_1)*(Rcoef/(2.0_dp*Scoef))+ &
          Tcoef*((Bplus*x3_1/(2.0_dp*Scoef)-Bminus*x2_1/(2.0_dp*Scoef)))+ &
          ion%begin_HII
    c_2 = (Rcoef+(Bminus)*Tcoef)/(2.0_dp*Scoef)
    c_3 = -(Rcoef+(Bplus)*Tcoef)/(2.0_dp*Scoef)

    ! arguments of exponential functions
    lambda_1_dt = dt*lambda_1
    lambda_2_dt = dt*lambda_2
    lambda_3_dt = dt*lambda_3

    ! exponential functions of homogeneous solutions
    exp_lambda_1_dt = exp(lambda_1_dt)
    exp_lambda_2_dt = exp(lambda_2_dt)
    exp_lambda_3_dt = exp(lambda_3_dt)

    ! ionization fractions at the end of the time-step
    ion%end_HII = c_1*exp_lambda_1_dt+c_2*exp_lambda_2_dt*x2_1+c_3*exp_lambda_3_dt*x3_1+xp_1 
    ion%end_HeII = c_2*exp_lambda_2_dt*x2_2+c_3*exp_lambda_3_dt*x3_2+xp_2
    ion%end_HeIII = c_2*exp_lambda_2_dt+c_3*exp_lambda_3_dt+xp_3
    ion%end_HI = 1.0_dp-ion%end_HII
    ion%end_HeI = 1.0_dp-ion%end_HeII-ion%end_HeIII

    !if (isnan(ion%end_HII)) write(*,*)'at position ',pos,' end HII'
	!if (isnan(ion%end_HeII)) write(*,*)'at position ',pos,' end HeII'
	!if (isnan(ion%end_HeIII)) write(*,*)'at position ',pos,' end HeIII'
	
    ! Make sure none of the ionization fractions less than epsilon
    if (ion%end_HI .lt. epsilon) then
      ion%end_HI = epsilon
      ion%end_HII = 1.0_dp-epsilon
    endif
 
    if (ion%end_HII .lt. epsilon) then
      ion%end_HII = epsilon
       ion%end_HI = 1.0_dp-epsilon
    endif
 
    if ((ion%end_HeI.le.epsilon) .or. (ion%end_HeII.le.epsilon) .or. (ion%end_HeIII.le.epsilon))then

      if (ion%end_HeI < epsilon) then
        ion%end_HeI = epsilon
      endif

      if (ion%end_HeII .lt. epsilon) then
        ion%end_HeII = epsilon      
      endif
	
      if (ion%end_HeIII .lt. epsilon) then
        ion%end_HeIII = epsilon
      endif
   
      normfac = ion%end_HeI+ion%end_HeII+ion%end_HeIII
      ion%end_HeI = ion%end_HeI/normfac
      ion%end_HeII = ion%end_HeII/normfac
      ion%end_HeIII = ion%end_HeIII/normfac

    endif
    
if (pos.eq.1 .or. pos.eq. 100) then
	!write(*,*) 'ne at pos ',pos,' is ', ndens_electron
	!write(*,*) 'f at pos ',pos,' is ', f
	!write(*,*) 'z at pos ',pos,' is ', z
	!write(*,*) 'v at pos ',pos,' is ', v
	!write(*,*) 'w at pos ',pos,' is ', w
	!write(*,*) 'y at pos ',pos,' is ', y
	!write(*,*) 'p at pos ',pos,' is ', p
	!write(*,*) 'y2a at pos ',pos,' is ', y2a
	!write(*,*) 'y2b at pos ',pos,' is ', y2b
	!write(*,*) 'alpha_HI_B at pos ',pos,' is ', alpha_HI_B
	!write(*,*) 'alpha_HeII_B at pos ',pos,' is ', alpha_HeII_B
	!write(*,*) 'alpha_HeII_2 at pos ',pos,' is ', alpha_HeII_2
	!write(*,*) 'alpha_HeI_B at pos',pos,' is ',alpha_HeI_B
	!write(*,*) 'alpha_HeI_1 at pos',pos,' is ',alpha_HeI_1
	!write(*,*) 'alpha_HeII_1 at pos ',pos,' is ', alpha_HeII_1
	!write(*,*) 'UHI at pos ',pos,' is ',g_1			
endif
	
! to check some constants
!if (pos .eq. 200) then
	!write(*,*) 'at pos ',pos
	!write(*,*) 'f = ',f
	!write(*,*) 'z = ',z	
	!write(*,*) 'v = ',v
	!write(*,*) 'w = ',w
	!write(*,*) 'y2a = ', y2a
	!write(*,*) 'y2b = ', y2b
	!write(*,*) 'alpha_HI_B = ',alpha_HI_B 		
	!write(*,*) 'alpha_HeII_B = ',alpha_HeII_B
	!write(*,*) 'alpha_HeII_2 = ',alpha_HeII_2
	!write(*,*) 'alpha_HeII_1 = ',alpha_HeII_1
	!write(*,*) '1st term = ', alpha_HI_B*ndens_electron
	!write(*,*) '2nd term = ', heliumfraction*ndens_electron*((f*z*(1-v)+v*w)*alpha_HeII_B)
	!write(*,*) '3rd term = ', heliumfraction*ndens_electron*(alpha_HeII_2)
	!write(*,*) '4th term = ', heliumfraction*ndens_electron*(1-y2a-y2b)*alpha_HeII_1
!endif
		
    ! Some coefficients required to calculate the average ionization fractions
    if (abs(lambda_1_dt) .lt. 1.0e-8) then
      avg_factor_1 = c_1
    else
      avg_factor_1 = c_1*(exp_lambda_1_dt-1.0_dp)/lambda_1_dt
    endif

    if (abs(lambda_2_dt) .lt. 1.0e-8) then
      avg_factor_2 = c_2
    else
      avg_factor_2 = c_2*(exp_lambda_2_dt-1.0_dp)/lambda_2_dt      
    endif

    if (abs(lambda_3_dt) .lt. 1.0e-8) then
      avg_factor_3=c_3
    else
      avg_factor_3 = c_3*(exp_lambda_3_dt-1.0_dp)/lambda_3_dt
    endif

    ! average ionization fractions of the time-step
    ion%avg_HII = xp_1 + avg_factor_1 + x2_1*avg_factor_2 + x3_1*avg_factor_3
    ion%avg_HeII = xp_2 + x2_2*avg_factor_2 + x3_2*avg_factor_3
    ion%avg_HeIII = xp_3 + avg_factor_2 + avg_factor_3
    ion%avg_HI = 1.0_dp-ion%avg_HII
    ion%avg_HeI = 1.0_dp-ion%avg_HeII-ion%avg_HeIII

    !if (isnan(ion%avg_HII)) write(*,*)'at position ',pos,' avg HII'
	!if (isnan(ion%avg_HeII)) write(*,*)'at position ',pos,' avg HeII'
	!if (isnan(ion%avg_HeIII)) write(*,*)'at position ',pos,' avg HeIII'

   ! Make sure none of the ionization fractions less than epsilon
    if (ion%avg_HII .lt. epsilon) then
       ion%avg_HII = epsilon
       ion%avg_HI = 1.0_dp-epsilon
    endif

    if (ion%avg_HI .lt. epsilon) then
      ion%avg_HI = epsilon
      ion%avg_HII = 1.0_dp-epsilon
    endif
 
    if ((ion%avg_HeI.le.epsilon) .or. (ion%avg_HeII.le.epsilon) .or. (ion%avg_HeIII.le.epsilon)) then
      if (ion%avg_HeII < epsilon) ion%avg_HeII = epsilon   
      if (ion%avg_HeIII < epsilon) ion%avg_HeIII = epsilon
      if (ion%avg_HeI < epsilon) ion%avg_HeI = epsilon
      normfac = ion%avg_HeI+ion%avg_HeII+ion%avg_HeIII
      ion%avg_HeI = ion%avg_HeI/normfac
      ion%avg_HeII = ion%avg_HeII/normfac
      ion%avg_HeIII = ion%avg_HeIII/normfac  
    endif
    
	if (isnan(ion%avg_HI) .and. report_counter.eq.0) then
		report_counter=report_counter+1
		write(*,*) 'alpha_HI_A', alpha_HI_A
		write(*,*) 'A_11', A_11
		write(*,*) 'A_13', A_13
		write(*,*) 'A_22', A_22
		write(*,*) 'A_33', A_33
		write(*,*) 'A_23', A_23
		write(*,*) 'A_32', A_32
		write(*,*) 'A_33', A_33
		write(*,*) 'two_A_32',two_A_32
		write(*,*) 'Bcoef',Bcoef
		write(*,*) 'Scoef',Scoef
		write(*,*) 'Kcoef',Kcoef
		write(*,*) 'Bminus',Bminus
		write(*,*) 'Bplus',Bplus
		write(*,*) 'lambda_1',lambda_1
		write(*,*) 'lambda_2',lambda_2
		write(*,*) 'lambda_3',lambda_3
		write(*,*) 'xp_1',xp_1
		write(*,*) 'xp_2',xp_3
		write(*,*) 'xp_3',xp_3
	endif
	
    return

  end subroutine doric
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Calculates the column density
  function coldens(path,neufrac,ndens,abundance)
    
    real(kind=dp) :: coldens
    real(kind=dp),intent(in) :: path     ! path length through cell
    real(kind=dp),intent(in) :: neufrac  ! neutral H/He fraction
    real(kind=dp),intent(in) :: ndens    ! total number density   
    real(kind=dp),intent(in) :: abundance ! the abundance of the element

    ! Column density over a cell
    coldens = neufrac*ndens*path*abundance
    
  end function coldens
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Sets the boundary condition for the HI column density
  function coldens_bndry_HI()
    
    use cgsphotoconstants, only: sigma_HI_at_ion_freq
    use radiation, only: boundary_tauHI
    
    real(kind=dp):: coldens_bndry_HI
    
    ! Column density at the left of the first cell
    coldens_bndry_HI = boundary_tauHI/sigma_HI_at_ion_freq
    
  end function coldens_bndry_HI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Sets the boundary condition for the HeI column density
  function coldens_bndry_HeI()
    use cgsphotoconstants, only: sigma_HeI_at_ion_freq
    use radiation, only: boundary_tauHeI

    real(kind=dp):: coldens_bndry_HeI  

    ! Column density at the left of the first cell
    coldens_bndry_HeI=boundary_tauHeI/sigma_HeI_at_ion_freq

  end function coldens_bndry_HeI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Sets the boundary condition for the HeII column density
  function coldens_bndry_HeII()
    use cgsphotoconstants, only: sigma_HeII_at_ion_freq
    use radiation, only: boundary_tauHeII

    real(kind=dp):: coldens_bndry_HeII  

    ! Column density at the left of the first cell
    coldens_bndry_HeII = boundary_tauHeII/sigma_HeII_at_ion_freq

  end function coldens_bndry_HeII

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module doric_module
  
