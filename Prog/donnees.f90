module donnees

  implicit none

  ! Maillage
  integer :: N, Nx, Ny
  integer :: Niter, Nmax, Ndisplay
  real*8  :: Lx, Ly, dx, dy, x, y
  real*8  :: t, dt, Tmax

  ! Thermique
  real*8  :: Text, Tad
  real*8  :: rho, cp, lambda
  real*8  :: h
  real*8  :: Ea, R, k0, Q, coef, tmp

  ! Modele  
  real*8, dimension(:), allocatable :: U, U0, V, Fm, eta, chi, eta0
  real*8                            :: D, a, b, c, epsilon

contains

end module donnees
