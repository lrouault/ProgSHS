module donnees

  implicit none

  integer, parameter :: PR=Selected_real_kind(12,60) !12,60
  
  ! Maillage
  integer :: N, Nx, Ny
  integer :: Niter, Nmax, Ndisplay
  real(PR)  :: Lx, Ly, dx, dy, x, y
  real(PR)  :: t, dt, Tmax

  ! Thermique
  real(PR)  :: Text, Tad
  real(PR)  :: rho, cp, lambda
  real(PR)  :: h
  real(PR)  :: Ea, R, k0, Q, coef

  ! Modele  
  real(PR), dimension(:), allocatable :: U, U0, V, Fm, eta, chi, eta0
  real(PR)                            :: D, a, b, c, epsilon

contains

end module donnees
