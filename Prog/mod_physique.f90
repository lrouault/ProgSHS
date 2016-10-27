module mod_physique
  use donnees
  use mod_fonction
  implicit none

contains

  subroutine creation_matrice()

    !a = 1. + (2.*D*dt)/(dx*dx) + (2.*D*dt)/(dy*dy)
    !b = -(D*dt) / (dy*dy)
    !c = -(D*dt) / (dx*dx)
    alpha = rho*cp + (2.*lambda*dt)/(dx**2) + (2.*lambda*dt)/(dy**2)
    beta  = -(lambda*dt)/(dx**2)
    gamma = -(lambda*dt)/(dy**2)

    
  end subroutine creation_matrice

  subroutine cd_limite()

    integer :: i,j
    V = 0
    do j = 1,Ny
       V(bij(1,j,Ny)) =2000*(lambda*dt)/(dx**2)
       !   V(bij(1,j,Ny)) = -c*Tad + V(bij(1,j,Ny))!-c*min(Tad,Tad*(500*t/Tmax))+V(bij(1,j,Ny))
       !   V(bij(Nx,j,Ny)) = 1000  *b*dy !-c*Text + V(bij(Nx,j,Ny))!-h*dt/(rho*cp*Lx)*(U0(bij(Nx,j,Ny))-Text)
    end do

    do i = 1,Nx
       !   V(bij(i,1,Ny))  =- 1000  
       !   V(bij(i,Ny,Ny))  =- 1000  *c*dy
       !V(bij(i,1,Ny)) = +V(bij(i,1,Ny)) - 10*(U(bij(i,1,Ny))-Text)!+V(bij(i,1,Ny))!-b*Text+V(bij(i,1,Ny))
       !V(bij(i,Ny,Ny)) = +V(bij(i,Ny,Ny)) - 10*(U(bij(i,Ny,Ny))-Text)!+V(bij(i,Ny,Ny))-b*Text+V(bij(i,Ny,Ny))
    end do

  end subroutine cd_limite

  subroutine eq_arrhenius()

    real(PR) :: tmp
    integer  :: i

    eta0 = eta
    tmp = 0.
    do i = 1,Nx*Ny
       coef   = -Ea/(R*U0(i))
       if (coef<-700) coef=-700
       !print*,"coeff, k0",coef,k0,dt !!!!!!!!!!!!!!!!!!!!!!
       eta(i) = (eta(i)+dt*k0*exp(coef)) / (1+k0*dt*exp(coef))
       !print*,"eta",eta(i)
       tmp    = max(tmp, eta(i)-eta0(i))
       !..!print*, "tmp", tmp
       chi(i) = Q/Cp * (eta(i)-eta0(i)) !k0*(1-eta(i))*exp(coef)/(1,77e+1)
    end do

  end subroutine eq_arrhenius
end module mod_physique
