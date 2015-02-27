! This module solves the photoionization equation 
! and calculates the timestep for evolution

module timeequation

  use precision, only: dp

  implicit none

contains


  subroutine time_equation(dt,ndens_electron,yHI,yHII,yHeII,yHeIII,phi,I_internal,xfinal,y,z,y2a,y2b,I_pos)

    use cgsconstants, only: bh00, bhe00, bhe10, albpow, alcpow, colh0, colhe, temphe, temph0, &
                            ev2fr,brech0,areche0,breche0, oreche0,breche1, treche1, n_el_crit, &
                            colli_HI, colli_HeI, colli_HeII, ini_rec_colion_factors, &
                            arech0, areche1, v
    use material, only: clumping, ionstates, ndens
    use abundances, only: abu_he
    use radiation, only: photrates

    ! time step
    real(kind=dp), intent(out) :: dt  
    ! electron density 
    real(kind=dp), intent(in) :: ndens_electron
    ! ionized states     
    real(kind=dp), intent(in) :: yHI, yHII, yHeII, yHeIII
    ! the I-front state
    integer, intent(inout) :: I_internal
    ! The timestep will be such that the cell will be evolved to xfinal
    real(kind=dp), dimension(:), intent(in) :: xfinal
    real(kind=dp), intent(inout) :: y
    real(kind=dp), intent(inout) :: z   
    real(kind=dp), intent(inout) :: y2a
    real(kind=dp), intent(inout) :: y2b
    integer, intent(in) :: I_pos
    ! Photoionization rate and heating rate
    type(photrates), intent(in)  :: phi 

    real(kind=dp) :: alpha_HI_B, alpha_HI_A
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
    real(kind=dp) :: A_11, A_12, A_13, A_22, A_23, A_32, A_33, two_A_32
    real(kind=dp) :: xp_1, xp_2, xp_3
    real(kind=dp) :: g_1, g_2
    real(kind=dp) :: heliumfraction
    real(kind=dp) :: bigdt, smalldt, averagedt
    real(kind=dp) :: resultant_xHII
    real(kind=dp) :: range_xfinal


    smalldt = 0.0
    bigdt = 1.0/(arech0*ndens(I_pos,1,1))
    averagedt = 0.5*(smalldt+bigdt)
    dt=0.0 !averagedt ! some value 

    ! see Osterbrock, 1989
    p = 0.96_dp
    ! fraction of photons from 2-photon decay, energetic enough to ionize hydrogen
    l = 1.425_dp
    ! fraction of photons from 2-photon decay, energetic enough to ionize neutral helium
    m = 0.737_dp 
    heliumfraction = abu_he/(1.0_dp-abu_he) 
    f = max(min(10.0_dp*yHI,1.0_dp),0.01_dp) 
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
    g_1 = max(phi%photo_cell_HI+ndens_electron*colli_HI,1.0e-200_dp)
    g_2 = max(phi%photo_cell_HeI+ndens_electron*colli_HeI,1.0e-200_dp)
    A_32 = max(phi%photo_cell_HeII+ndens_electron*colli_HeII,1.0e-200_dp)     
    A_11 = -(g_1+ndens_electron*alpha_HI_B)
    A_12 = (y*alpha_HeI_1+p*alpha_HeI_B)*ndens_electron*heliumfraction
    A_13 = ((f*z*(1.0_dp-v) +v*w)*alpha_HeII_B +alpha_HeII_2 &
           +(1.0_dp-y2a-y2b)*alpha_HeII_1)*heliumfraction*ndens_electron  
    A_22 = -g_2-A_32-ndens_electron*(alpha_HeI_A-(1.0_dp-y)*alpha_HeI_1)
    A_33 = -ndens_electron*(alpha_HeII_A-y2a*alpha_HeII_1)            
    A_23 = -g_2+ndens_electron*alpha_HeII_B*(f*(1.0_dp-z)*(1.0_dp-v)+ &
           v*(1.425_dp-w))-A_33+alpha_HeII_1*y2b*ndens_electron

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
    Rcoef = two_A_32 *(xp_2-yHeII)
    Tcoef = xp_3-yHeIII

    ! coefficients of homogeneous solutions
    c_1 = -xp_1+ (x3_1-x2_1)*(Rcoef/(2.0_dp*Scoef))+ &
          Tcoef*((Bplus*x3_1/(2.0_dp*Scoef)-Bminus*x2_1/(2.0_dp*Scoef)))+ &
          yHII
    c_2 = (Rcoef+(Bminus)*Tcoef)/(2.0_dp*Scoef)
    c_3 = -(Rcoef+(Bplus)*Tcoef)/(2.0_dp*Scoef)

   ! When the ionized fraction of the I-front laready surpasses the targeted ionized fraction
   ! plus one to the I-front state
   do while (yHII.gt.xfinal(I_internal).and.I_internal.le.size(xfinal))
      I_internal = I_internal+1
   enddo

   ! the target ionized fraction is a range
   ! range = [xfinal, xfinal+range_xfinal] 
   range_xfinal = 0.2/(size(xfinal)+1)

   
   if (I_internal.le.size(xfinal)) then

     ! Iteration to guess the timestep required to evolve a cell to the desired ionized fraction
     do 
	   !write(*,*)'dt is ',dt
       ! arguments of exponential functions
       lambda_1_dt = dt*lambda_1
       lambda_2_dt = dt*lambda_2
       lambda_3_dt = dt*lambda_3

       ! exponential functions of homogeneous solutions
       exp_lambda_1_dt = exp(lambda_1_dt)
       exp_lambda_2_dt = exp(lambda_2_dt)
       exp_lambda_3_dt = exp(lambda_3_dt)

       ! ionization fractions at the end of the time-step
       resultant_xHII = c_1*exp_lambda_1_dt+c_2*exp_lambda_2_dt*x2_1+c_3*exp_lambda_3_dt*x3_1+xp_1 

       ! check if the time-step reaches the requirement
       if (resultant_xHII .ge. xfinal(I_internal) .and. &
           resultant_xHII .le. (xfinal(I_internal)+range_xfinal) ) then
         exit
       endif

       ! if the time-step is underestimated
       if (resultant_xHII .le. xfinal(I_internal)) then
         smalldt = averagedt
         averagedt = 0.5*(smalldt+bigdt)
         dt = averagedt
       endif

       ! if the time-step is overestimated
       if (resultant_xHII .ge. (xfinal(I_internal)+range_xfinal)) then
         bigdt = averagedt
         averagedt = 0.5*(smalldt+bigdt)
         dt = averagedt
       endif    

     enddo

   endif

   ! If the I-front state is larger than the partition size
   ! then zero time-step is assigned so that IGM will evolve nothing
   ! but the I-front proceeds to the next cell
   if (I_internal.gt.size(xfinal)) dt = 0.0

  end subroutine time_equation
  
end module timeequation
  
