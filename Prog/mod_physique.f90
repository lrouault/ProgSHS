module mod_physique
  use donnees

  implicit none

contains

  subroutine creation_matrice()

    D = lambda / (rho*cp)
    a = 1. + (2.*D*dt)/(dx*dx) + (2.*D*dt)/(dy*dy)
    b = -(D*dt) / (dy*dy)
    c = -(D*dt) / (dx*dx)

  end subroutine creation_matrice

end module mod_physique
