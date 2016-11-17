module donnees

  implicit none

  integer, parameter :: PR=Selected_real_kind(12,60) !12,60

  ! Maillage
  integer   :: Nx, Ny, Nz
  integer   :: Niter, Ndisplay
  real(PR)  :: Lx, Ly, Lz, dx, dy, dz, x, y, z
  real(PR)  :: time, dt, tmax

  ! Thermique
  real(PR)  :: Text, Tad
  !real(PR)  :: rho, cp, lambda
  real(PR),dimension(:),allocatable   :: rho, rhocp, lambda
  real(PR)  :: h, flux
  real(PR)  :: Ea, R, k0, Q, coef

  ! Modèle
  real(PR), dimension(:), allocatable :: U, U0, rhs, Fm, eta, chi, eta0
  real(PR)                            :: epsilon
  real(PR), dimension(:), allocatable :: Cd
  real(PR), dimension(:), allocatable :: Cx, Cy, Cz

  integer :: nb_fichiers

  ! Materiau
  real(PR),dimension(:,:),allocatable :: fraction_vol ! (/Si,N2,fibre/)
  real(PR)                          :: rho_Si, rho_Si3N4, rho_N2, rho_fibre
  real(PR),dimension(:,:),allocatable :: cp_Si, cp_Si3N4, cp_N2, cp_fibre
  real(PR),dimension(:,:),allocatable :: lambda_Si, lambda_Si3N4, lambda_N2, lambda_fibre

  ! fibre
  integer, dimension(:,:,:),allocatable :: Poro
  integer, dimension(:,:,:,:),allocatable :: Orien
  integer :: porox, poroy, poroz

  integer :: iteration

contains

end module donnees
