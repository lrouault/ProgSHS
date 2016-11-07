module donnees

  implicit none

  integer, parameter :: PR=Selected_real_kind(12,60) !12,60

  ! Maillage
  integer   :: N, Nx, Ny
  integer   :: Niter, Nmax, Ndisplay
  real(PR)  :: Lx, Ly, dx, dy, x, y
  real(PR)  :: time, dt, tmax

  ! Thermique
  real(PR)  :: Text, Tad
  real(PR)  :: rho, cp, lambda
  real(PR)  :: h, flux
  real(PR)  :: Ea, R, k0, Q, coef

  ! Modèle  
  real(PR), dimension(:), allocatable :: U, U0, rhs, Fm, eta, chi, eta0
  real(PR)                            :: D, a, b, c, epsilon
  real(PR), dimension(:), allocatable :: Cd
  real(PR), dimension(:), allocatable :: Cx
  real(PR), dimension(:), allocatable :: Cy

  integer :: nb_fichiers

contains

end module donnees
