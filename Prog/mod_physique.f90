module mod_physique
  use donnees
  use mod_fonction
  implicit none

contains

  subroutine creation_matrice()

    Cd = rho*cp + (2.*lambda*dt)/(dx**2) + (2.*lambda*dt)/(dy**2)
    Cx = -(lambda*dt)/(dx**2)
    Cy = -(lambda*dt)/(dy**2)

  end subroutine creation_matrice

  subroutine second_membre()
    integer :: i,j,k
    real(PR),dimension(Nx*Ny)::phi

    rhs = 0.
    phi = eq_arrhenius()

    do k=1,Nx*Ny
       rhs(k) = rho*Cp*U0(k)+rho*Q*phi(k)
    end do

    do i=1,Nx
       k=i !bas
       rhs(k) = rhs(k) - Cy(1)*(U(k)-h*dy*(U(k)-Text)/lambda)

       k=(Ny-1)*Nx+i !haut
       rhs(k) = rhs(k) - Cy(1)*(U(k)-h*dy*(U(k)-Text)/lambda)
    end do

    do j=1,Ny
       k=(j-1)*Nx+1 !gauche
       rhs(k) = rhs(k) - Cx(1)*(U(k)-h*dx*(U(k)-Text)/lambda)

       k=(j-1)*Nx+Nx !droite
       rhs(k) = rhs(k) - Cx(1)*(U(k)-h*dx*(U(k)-Text)/lambda)
    end do

    rhs((Ny/2)*Nx+1) = rhs((Ny/2)*Nx+1) - Cx(1)*flux*dx/lambda
    
  end subroutine second_membre

  function eq_arrhenius()
    integer  :: i
    real(PR),dimension(Nx*Ny)::eq_arrhenius

    eta0 = eta

    do i = 1,Nx*Ny
       coef   = -Ea/(R*U0(i))
       !if (coef<-700) coef=-700

       eta(i) = (eta(i)+dt*k0*exp(coef)) / (1+k0*dt*exp(coef))

       eq_arrhenius(i) = eta(i)-eta0(i)
    end do

  end function eq_arrhenius
end module mod_physique
