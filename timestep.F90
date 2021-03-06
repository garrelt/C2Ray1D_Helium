! to get the suitable timestep

module timestep

  use precision, only: dp
  use grid, only: vol, dr
  use material, only: ndens, xHI, xHII, xHeI, xHeII, xHeIII, temper 
  use abundances, only: abu_h, abu_he
  use tped, only: electrondens
  use doric_module, only: coldens, coldens_bndry_HI, coldens_bndry_HeI, coldens_bndry_HeII
  use timeequation, only: time_equation
  use radiation, only: photoion, photrates
  use evolve, only: coldens_out_HI, coldens_out_HeI, coldens_out_HeII
  use material, only: gamma_uvb
  use cgsconstants, only : ini_rec_colion_factors
  use cgsphotoconstants, only: sigma_H_heth, sigma_H_heLya,sigma_He_heLya,sigma_HeI_at_ion_freq,&
                               sigma_He_he2,sigma_HeII_at_ion_freq,sigma_He_he2,sigma_H_he2

  implicit none

  public :: time_step 

contains

  subroutine time_step(dt,I_pos,I_internal,xfinal)

    ! Timestep is calculated here for the evolution
    real(kind=dp), intent(out) :: dt 
    ! the position of cell that the I-front situates current
    integer, intent(inout) :: I_pos
    ! the I-front state
    integer, intent(inout) :: I_internal
    ! The timestep will be such that the cell will be evolved to xfinal
    real(kind=dp), dimension(:), intent(in) :: xfinal

    ! Size of a cell
    real(kind=dp) :: path
    ! Parameters to solve coupled ODEs
    real(kind=dp) :: tau_H_heth, tau_H_heLya, tau_He_heth, tau_He_heLya, y, z
    real(kind=dp) :: y2a, y2b, tau_He2_he2th, tau_He_he2th, tau_H_he2th
    ! Volume of the cell
    real(kind=dp) :: vol_cell
    ! Ionized fractions of I-front
    real(kind=dp) :: yHI, yHII, yHeI, yHeII, yHeIII
    ! Temperature of I-front
    real(kind=dp) :: temper_cell
    ! Atom density of I-front
    real(kind=dp) :: ndens_atom_cell
    ! Electron density of I-front
    real(kind=dp) :: ndens_electron_cell
    ! Incoming column densities of I-front
    real(kind=dp) :: coldensHI_in
    real(kind=dp) :: coldensHeI_in
    real(kind=dp) :: coldensHeII_in
    ! Intermidiate column densities of I-front
    real(kind=dp) :: coldensHI_cell
    real(kind=dp) :: coldensHeI_cell
    real(kind=dp) :: coldensHeII_cell
    ! Photoionization rate and heating rate
    type(photrates) :: phi

    ! copy the cell ionization fractions of I-front
    yHI = xHI(I_pos)
    yHII = xHII(I_pos)
    yHeI = xHeI(I_pos)
    yHeII = xHeII(I_pos)
    yHeIII = xHeIII(I_pos)
   
    ! copy the cell temperature of I-front
    temper_cell = temper(I_pos)  

    ! copy the cell atomic number density of I-front    
    ndens_atom_cell = ndens(I_pos,1,1)

    ! copy outgoing column density of cell behind I-front
    ! as the incoming column density of current cell
    if (I_pos.eq.1) then
       coldensHI_in = coldens_bndry_HI()
       coldensHeI_in = coldens_bndry_HeI()
       coldensHeII_in = coldens_bndry_HeII()
    else
       coldensHI_in = coldens_out_HI(I_pos-1)
       coldensHeI_in = coldens_out_HeI(I_pos-1)
       coldensHeII_in = coldens_out_HeII(I_pos-1)
    endif

    ! Cell size
    path = dr(1)

    ! Find the volume of the shell this cell is part of (dilution factor)
    vol_cell = vol(I_pos)

    ! Average HI, HeI and HeII column densities 
    coldensHI_cell = coldens(path,yHI,ndens_atom_cell,abu_h)
    coldensHeI_cell = coldens(path,yHeI,ndens_atom_cell,abu_he)
    coldensHeII_cell = coldens(path,yHeII,ndens_atom_cell,abu_he)

    ! Obtain photoionizing rate and heating rate on the I-front
    call photoion(phi,coldensHI_in,coldensHI_in+coldensHI_cell, &
                  coldensHeI_in,coldensHeI_in+coldensHeI_cell, &
                  coldensHeII_in,coldensHeII_in+coldensHeII_cell, &
                  vol_cell, 1, yHII)

    ! True HI photoionization rate
    phi%photo_cell_HI = phi%photo_cell_HI/(ndens_atom_cell*abu_h*yHI)

    ! Add the HI UV background
    phi%photo_cell_HI = phi%photo_cell_HI + gamma_uvb(1)

    ! Calculate (mean) electron density
    ndens_electron_cell = electrondens(ndens_atom_cell,yHII, yHeII, yHeIII)
   
    ! These are the specific optical depth in order to find out y,z,y2a,y2b
    tau_H_heth = coldensHI_cell*sigma_H_heth 
    tau_He_heth = coldensHeI_cell*sigma_HeI_at_ion_freq 
    tau_H_heLya = coldensHI_cell*sigma_H_heLya
    tau_He_heLya= coldensHeI_cell*sigma_He_heLya
    tau_He2_he2th = coldensHeII_cell*sigma_HeII_at_ion_freq
    tau_He_he2th = coldensHeI_cell*sigma_He_he2
    tau_H_he2th = coldensHI_cell*sigma_H_he2

    !Parameters required to define the coupled ODEs
    y = tau_H_heth/(tau_H_heth +tau_He_heth)
    z = tau_H_heLya/(tau_H_heLya+tau_He_heLya)
    y2a = tau_He2_he2th/(tau_He2_he2th +tau_He_he2th+tau_H_he2th)
    y2b = tau_He_he2th/(tau_He2_he2th +tau_He_he2th+tau_H_he2th)

    ! Update the collisional ion and recombination rates for new temperature
    call ini_rec_colion_factors(temper_cell) 

    ! Obtain the evolution timestep
    call time_equation(dt,ndens_electron_cell,yHI,yHII,yHeII,yHeIII,phi, &
                       I_internal,xfinal,y,z,y2a,y2b,I_pos)
       
  end subroutine time_step

end module timestep
