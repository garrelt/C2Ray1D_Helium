!>
!! \brief This module contains data and routines for handling
!! the material properties on the grid (1D)
!!
!! These properties are; density, temperature, clumping, ionization fractions
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date:
!!
!! \b Version: 1D test problems\n
!! Problem 1: constant density (Strömgren problem)\n
!! Problem 2: 1/r density\n
!! Problem 3: 1/r^2 density (with core of radius r_core)\n
!! Problem 4: cosmological constant density (Shapiro & Giroux problem)\n


module material

  ! This module contains the grid data and routines for initializing them.
  ! These are
  ! - mat_ini : initializes temperature and ionization fractions at start

  ! Version: 1D test problems

  ! Problem 1: constant density (Strömgren problem)
  ! Problem 2: 1/r density
  ! Problem 3: 1/r^2 density (with core of radius r_core)
  ! Problem 4: cosmological constant density (Shapiro & Giroux problem)

    use mathconstants, only: pi
  use precision, only: dp,si
  use cgsconstants, only: bh00,bhe00,bhe10,ini_rec_colion_factors,m_p
  use cgsconstants, only: brech0,breche0,breche1,arech0,areche0,areche1
  use astroconstants, only: YEAR, mpc
  use sizes, only: mesh
  use file_admin, only: stdinput, file_input
  use my_mpi
  use grid, only: x,vol,dr
  use c2ray_parameters, only: epsilon
use abundances, only: abu_h,abu_he,mu
  use cosmology_parameters, only: Omega0, H0, rho_crit_0,Omega_B
  !use cosmology, only: cosmology_init,H0,t0,zred_t0

  implicit none

  ! ndens - number density (cm^-3) of a cell
  ! temper - temperature (K) of a cell
  ! xh - ionization fractions for one cell
  real(kind=dp) :: ndens(mesh,1,1) !< number density (cm^-3) of a cell
  real(kind=dp) :: temper(mesh) !< temperature (K) of a cell
  real(kind=dp) :: xHI(mesh) 
  real(kind=dp) :: xHII(mesh) 
  real(kind=dp) :: xHeI(mesh) 
  real(kind=dp) :: xHeII(mesh) 
  real(kind=dp) :: xHeIII(mesh) 
  real(kind=dp) :: HI_photoionization_rate(mesh)
  real(kind=dp) :: recombination_time(mesh)
  real(kind=dp) :: recombination_rate(mesh)
  real(kind=dp) :: avg_ex_ph_eV(mesh)
  real(kind=dp) :: clumping !< global clumping factor
  real(kind=dp) :: r_core !< core radius (for problems 2 and 3)
  real(kind=dp) :: dens_core !< core density (for problems 2 and 3)
  integer :: testnum !< number of test problem (1 to 4)
  logical :: isothermal !< is the run isothermal?
  real(kind=dp),dimension(3) :: gamma_uvb !< UV background for HI, HeI, HeII
  ! needed for analytical solution of cosmological Ifront
  real(kind=dp) :: t1 !< parameter for analytical solution of test 4
  real(kind=dp) :: t0_t !< parameter for analytical solution of test 4
  real(kind=dp) :: eta !< parameter for analytical solution of test 4
  real(kind=dp) :: zred00
  real(kind=dp),public :: n_LLS ! just because cosmology needs it in the 3D version
  ! b) Model Songaila & Cowie (2010)
  real(kind=dp) :: y_LLS
real(kind=dp) :: avg_ndens

!*TEST******************************************************
   type ionstates
     real(kind=dp) :: end_HI 
     real(kind=dp) :: end_HII
     real(kind=dp) :: end_HeI
     real(kind=dp) :: end_HeII
     real(kind=dp) :: end_HeIII
     real(kind=dp) :: avg_HI
     real(kind=dp) :: avg_HII
     real(kind=dp) :: avg_HeI
     real(kind=dp) :: avg_HeII
     real(kind=dp) :: avg_HeIII
     real(kind=dp) :: begin_HI
     real(kind=dp) :: begin_HII
     real(kind=dp) :: begin_HeI
     real(kind=dp) :: begin_HeII
     real(kind=dp) :: begin_HeIII

     !real(kind=dp) :: h(0:1) !< H ionization fractions
     !real(kind=dp) :: he(0:2) !< He ionization fractions
     !real(kind=dp) :: h_av(0:1) !< average H ionization fractions
     !real(kind=dp) :: he_av(0:2) !< average He ionization fractions
     !real(kind=dp) :: h_old(0:1) !< H ionization fractions from last time step
     !real(kind=dp) :: he_old(0:2) !< He ionization fractions from last time step
  end type ionstates
 


!***********************************************************

#ifdef MPI
  integer,private :: ierror !< MPI error flag
#endif

contains

  ! ============================================================================

  !> Initializes material properties on grid\n
  !! \b Author: Garrelt Mellema\n
  !! \b Date: 20-Aug-2006 (f77 21-May-2005 (derives from mat_ini_cosmo2.f))\n
  !! \b Version:
  !! - 1D\n
  !! - Four different test problems\n
  !! - Initially completely neutral\n

  subroutine mat_ini (background_HII)

    ! Initializes material properties on grid

    ! Author: Garrelt Mellema

    ! Date: 20-Aug-2006 (f77 21-May-2005 (derives from mat_ini_cosmo2.f))

    ! Version:
    ! - 1D
    ! - Four different test problems
    ! - Initially completely neutral

real(kind=dp), intent(out) :: background_HII
    integer :: i,n ! loop counters
    real(kind=dp) :: dens_val
    real(kind=dp) :: temper_val
    real(kind=dp) :: alpha
    real(kind=dp) :: z_quasar
    real(kind=dp) :: overdensity
    real(kind=dp) :: rb
    real(kind=dp) :: xMpc
    real(kind=dp), parameter :: nh_0 = 1.933e-7
    character(len=1) :: answer
 real(kind=dp) :: dr1,dr2	
    real(kind=dp),dimension(3) :: xions
   character(len=40) :: file3
    type(ionstates) :: ion
    real(kind=dp) :: normfac
    ! Ask for input
    if (rank == 0) then
if (.not.file_input) then
write(*,'(A,$)') 'Which test? (1-5): '
       endif
read(stdinput,*) testnum
write(logf,*) 'Test problem ', testnum
#ifdef MPI
       ! Distribute the input parameters to the other nodes
       call MPI_BCAST(testnum,1,MPI_INTEGER,0,MPI_COMM_NEW,ierror)
#endif
endif

    ! Set alpha according to test problem
    select case (testnum)
    case(1,4)
       alpha=0.0
    case(2)
       alpha=-1.0
    case(3)
       alpha=-2.0
    end select

if (rank == 0) then
  if (testnum == 1) then
    if (.not.file_input) write(*,'(A,$)') 'Enter density (cm^-3): '
    read(stdinput,*) dens_val
    write(logf,*) 'Density = ', dens_val
  elseif (testnum == 2 .or. testnum == 3) then
    if (.not.file_input) write(*,'(A,$)') 'Enter reference (core) radius (cm): '
    read(stdinput,*) r_core
    write(logf,*) 'Core radius = ', r_core 
    if (.not.file_input) write(*,'(A,$)') 'Enter density at reference (core) radius(cm^-3): '
    read(stdinput,*) dens_val
    write(logf,*) 'Density at core = ', dens_val
  elseif (testnum == 4 ) then 
    if (.not.file_input) write(*,'(A,$)') 'Enter redshift of the quasar: '
    read(stdinput,*) z_quasar
    write(logf,*) 'redshift = ', z_quasar
  elseif (testnum == 5 ) then 
    if (.not.file_input) write(*,'(A,$)') 'Enter redshift of the quasar: '
    read(stdinput,*) z_quasar
    write(logf,*) 'redshift = ', z_quasar
    if (.not.file_input) write(*,'(A,$)') 'Enter overdensity of the quasar: '
    read(stdinput,*) overdensity
    write(logf,*) 'overdensity = ', overdensity
  endif

if (.not.file_input) write(*,'(A,$)') 'Enter clumping factor: '
       read(stdinput,*) clumping
       write(logf,*) 'Clumping factor = ', clumping 
       if (.not.file_input) write(*,'(A,$)') 'Enter initial temperature (K): '
       read(stdinput,*) temper_val
       write(logf,*) 'Initial temperature = ', temper_val
       if (.not.file_input) write(*,'(A,$)') 'Isothermal? (y/n): '
       read(stdinput,*) answer
       write(logf,*) 'Isothermal? (y/n): ', answer
       ! Isothermal?
       if (answer == 'y' .or. answer == 'Y') then
          isothermal=.true.
       else
          isothermal=.false.
       endif
if (.not.file_input) write(*,'(A,$)') 'Ionizing background (HI,HeI,HeII) (s^-1): '
       read(stdinput,*) gamma_uvb
       write(logf,*) 'Ionizing background of HI (s^-1):', gamma_uvb(1)
       write(logf,*) 'Ionizing background of HeI (s^-1):', gamma_uvb(2)
       write(logf,*) 'Ionizing background of HeII (s^-1):', gamma_uvb(3)
       call ini_rec_colion_factors(temper_val) !initialize the collisional ion and recomb rates for inital temp

endif
#ifdef MPI
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(dens_val,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
    if (testnum == 2.or.testnum == 3) &
         call MPI_BCAST(r_core,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
    call MPI_BCAST(clumping,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
    call MPI_BCAST(temper_val,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
    call MPI_BCAST(isothermal,1,MPI_LOGICAL,0,MPI_COMM_NEW,ierror)
    call MPI_BCAST(gamma_uvb,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
#endif

    ! For test problem 4: cosmological parameters
    !if (testnum == 4) then
       ! Ask for cosmological parameters
    !   if (rank == 0) then
!write(*,'(A,$)') 'Initial redshift?'
    !      read(stdinput,*) zred00
    !      write(*,*) 'redshift=', zred00
    !   endif
!#ifdef MPI
       ! Distribute the input parameters to the other nodes
 !      call MPI_BCAST(zred00,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
!#endif

! call cosmology_init (.true.)
! else
! call cosmology_init (.false.)
    !endif

    ! Assign density and temperature to grid
    
    select case (testnum)
    case(1)
       do i=1,mesh
          ndens(i,1,1)=dens_val
          temper(i)=temper_val
       enddo
	   
	   !ndens(1,1,1)=ndens(1,1,1)*700.0




	
	 
	   !ndens(460,1,1)=10*ndens(50,1,1)
	   
	   !ndens(100,1,1)=20*ndens(100,1,1)
	   	   
	   !ndens(200,1,1)=10*ndens(200,1,1)
	   
	   !ndens(300,1,1)=10*ndens(300,1,1)	   		   
    case(2,3)
       dens_core=dens_val
       do i=1,mesh
          ! This is an attempt to make the initial conditions more
          ! correct: the density of a cell is the integrated density
          ! (mass) divided by the volume. Seems to give a worse fit
          ! to the analytical solution (slightly)
          ! rl=r(i)-0.5*dr
          ! rr=r(i)+0.5*dr
          ! ndens(i)=4.0*pi*dens_val*r_core**(-alpha)/vol(i)*
          ! $ (rr**(3.0+alpha)-rl**(3.0+alpha))/(3.0+alpha)
          
          ! This is just a straight sampling of the density distribution
          ! using the value at the cell centre.
          if (testnum == 3.and.x(i) <= r_core) then
             ! Flat core for test 3
             ndens(i,1,1)=dens_val
          else
             ndens(i,1,1)=dens_val*(x(i)/r_core)**alpha
          endif
             temper(i)=temper_val
       enddo
       
    case(4)
       ! For cosmological simulations, mean IGM
       !dens_core=dens_val ! save initial density value
       ! Parameters needed for the analytical solution
       ! Recombination time at z0
       !t1 = 1./(bh00*clumping*dens_core)
       !t0_t = 2.0_dp*(1.+zred00)**(-1.5)/(3.*H0*sqrt(Omega0))
       !eta = t0_t/t1*(1.+zred00)**3
       
	   !dens_val = 2.039943884e-7 * (1.0+z_quasar)**3
	   dens_val = Omega_B*(rho_crit_0/m_p)*(1.0/(1.0+3.0*abu_he))*(1.0+z_quasar)**3
	   write(*,*) 'rho_bar',rho_crit_0*Omega_B
	   write(*,*) 'density is ',dens_val
       ! Set the density to the comoving value
       ! (as is the spatial coordinate)
       ! evol_cosmo will set it to proper values.
       do i=1,mesh
          ndens(i,1,1)=dens_val
          temper(i)=temper_val
       enddo
       !scale below would be wrong
    case(5)
       rb = 0.4/(1.0+z_quasar)
       do i=1,mesh
          xMpc = x(i)/mpc
		  ndens(i,1,1)=0.8*nh_0*(1+z_quasar)**3*( 1.0+ overdensity*((1.0**2+rb**2)/(xMpc**2+rb**2))**0.7 )
          !ndens(i,1,1)=nh_0*(1+z_quasar)**3*( 1.0+ overdensity*((1.0**2+rb**2)/(xMpc**2+rb**2))**0.7 )
          temper(i)=temper_val
       enddo
	   
	   !ndens(1,1,1)=5.086e-04*419.013
	   !ndens(2,1,1)=6.901e-05*419.013
	   !ndens(3,1,1)=1.797e-06*419.013
	   !ndens(31,1,1)=1.388e-06*419.013
	   !ndens(32,1,1)=1.384e-06*419.013
	   !ndens(33,1,1)=1.195e-06*419.013
	   write(*,*)'flat density = ',ndens(mesh,1,1)
	   !do i=410,mesh
	!do i=320,mesh	   
     !      ndens(i,1,1)=2.3*ndens(1,1,1)
      ! enddo		  

	   !do i=330,mesh
       !    ndens(i,1,1)=2.35*ndens(1,1,1)
       !    temper(i)=temper_val
       !enddo	

	   !do i=335,mesh
       !    ndens(i,1,1)=2.5*ndens(1,1,1)
       !    temper(i)=temper_val
       !enddo	
	   	      		  
    end select

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! read file !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (1<0)then
    file3='ntot_6.483_xneg_1.txt'
	!file3='ntot_6.483_xpos_1.txt'
	!file3='ntot_6.483_yneg_1.txt'
	!file3='ntot_6.483_ypos_1.txt'
	!file3='ntot_6.483_zneg_1.txt'
	!file3='ntot_6.483_zpos_1.txt'
    !file3='ntot_6.483_xneg_2.txt'
	!file3='ntot_6.483_xpos_2.txt'
	!file3='ntot_6.483_yneg_2.txt'
	!file3='ntot_6.483_ypos_2.txt'
	!file3='ntot_6.483_zneg_2.txt'
	!file3='ntot_6.483_zpos_2.txt'
			
	avg_ndens = 0.0
	    open(unit=71,file=file3,form='formatted',status='unknown')
	   	do i=1,mesh		
	       read(71,*) x(i), ndens(i,1,1)
    		   x(i)=x(i)*Mpc
		   ndens(i,1,1)=ndens(i,1,1)*419.013
		   temper(i)=temper_val
		   avg_ndens = avg_ndens + ndens(i,1,1)
	   	enddo 	
		
		write(*,*)'The avg density of file ',file3,' is ',avg_ndens/real(mesh)
	    dr1= (x(mesh)-x(i))/real(mesh-1)
	    do i=1,mesh
	       x(i)=x(i)+0.5*dr1
	       vol(i)=4.0*pi/3.0*((x(i)+0.5*dr1)**3-(x(i)-0.5*dr1)**3)		   
	    enddo
   endif
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! Assign ionization fractions
    ! Use Gamma_UVB_H for this if it is not zero
    if (gamma_uvb(1) > 0.0) then
      do i=1,mesh
          call find_ionfractions_from_uvb(i,ndens(i,1,1), xions)
          xHI(i)=xions(1)
          xHII(i)=1.0-xHI(i)
          xHeI(i)=xions(2)
          xHeII(i)=xions(3)
          xHeIII(i)=1.0-(xHeI(i)+xHeII(i))

      enddo
    else
      do i=1,mesh
          xHI(i) =1.0_dp!1.0_dp-epsilon!-1.0e-14
          xHII(i) =0.0_dp!epsilon!1.0e-14
          xHeI(i)=1.0_dp-2.0_dp*epsilon!1.0_dp-epsilon!1.0_dp-2.0e-14 !1.0_dp-2.4e-9
          xHeII(i)=epsilon!1.0e-14 !1.2e-9
          xHeIII(i)=epsilon!1.0e-14 !1.2e-9
      enddo
    endif
	
    ! Make sure none of the ionization fractions less than epsilon
do i=1,mesh	
     if (xHI(i) .lt. epsilon) then
        xHI(i) = epsilon
        xHII(i) = 1.0_dp-epsilon
     endif

     if (xHII(i) .lt. epsilon) then
        xHII(i) = epsilon
        xHI(i) = 1.0_dp-epsilon
     endif
 
     if ((xHeI(i).le.epsilon) .or. (xHeII(i).le.epsilon) .or. (xHeIII(i).le.epsilon)) then
       if (xHeI(i) < epsilon) xHeI(i) = epsilon   
       if (xHeII(i) < epsilon) xHeII(i) = epsilon 
       if (xHeIII(i) < epsilon) xHeIII(i) = epsilon 
       normfac = xHeI(i)+xHeII(i)+xHeIII(i)
       xHeI(i) = xHeI(i)/normfac
       xHeII(i) = xHeII(i)/normfac
	   xHeIII(i) = xHeIII(i)/normfac  
     endif
enddo
	 	
    write(*,*) 'xHI is ', xHI(1)
    write(*,*) 'xHII is ', xHII(1)	
    write(*,*) 'xHeI is ',xHeI(1)
    write(*,*) 'xHeII is ', xHeII(1)
    write(*,*) 'xHeIII is ', xHeIII(1)		
	
background_HII = xHII(1)

    ! Report recombination time scale (in case of screen input)
    if (.not.file_input) write(*,'(A,1pe10.3,A)') 'Recombination time scale: ', &
         1.0/(dens_val*clumping*bh00*YEAR),' years'

    
  end subroutine mat_ini

  subroutine find_ionfractions_from_uvb (ii,nnd,xions)

    real(kind=dp),parameter :: convergence=0.01
    integer,intent(in) :: ii
    real(kind=dp),intent(in) :: nnd
    real(kind=dp),dimension(3),intent(out) :: xions

    real(kind=dp) :: rech2
    real(kind=dp) :: reche2
    real(kind=dp) :: reche3
    real(kind=dp) :: fe
    real(kind=dp) :: fe_prev

    rech2 = nnd * clumping * brech0
    reche2 = nnd * clumping* breche0
    reche3 = nnd * clumping* breche1
    fe=1.0
    
    ! Iterate to find the proper electron density (fe)
    do
       xions(1)=fe*rech2/(gamma_uvb(1)+fe*rech2)
       xions(2)=fe*reche2/(gamma_uvb(2)*(1.0+gamma_uvb(3)/(fe*reche3)) + &
            fe*reche2)
       xions(3)=(1.0-xions(2))/(1.0+gamma_uvb(3)/(fe*reche3))
       fe_prev=fe
       fe=abu_h*(1.0-xions(1))+abu_he*(2.0-(2.0*xions(2)+xions(3)))
       if (abs(fe-fe_prev)/fe_prev < convergence) exit
    enddo

  end subroutine find_ionfractions_from_uvb

end module material
