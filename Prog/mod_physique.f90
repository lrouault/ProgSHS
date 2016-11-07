module mod_physique
  use donnees
  use mod_fonction
  implicit none

contains

  subroutine creation_matrice()
    
    Cd = rho*cp + (2.*lambda*dt)/(dx**2) + (2.*lambda*dt)/(dy**2)
    Cx  = -(lambda*dt)/(dx**2)
    Cy = -(lambda*dt)/(dy**2)
    
  end subroutine creation_matrice

  subroutine cd_limite()
    integer :: i,j,k
    
    V = 0

    do i=1,Nx
       do j=1,Ny
          k=(j-1)*Nx+i
          
          if(i==1) then
             V(k) = V(k) + min(Tad, 2*Text+Tad*t)*(lambda*dt)/(dx**2)
          else if(i==Nx) then
             V(k) = V(k) + U(k)*(lambda*dt)/(dx**2)
          end if

          if(j==1) then
             v(k) = V(k) + U(k)*(lambda*dt)/(dy**2)
          else if(j==Ny) then
             V(k) = V(k) + U(k)*(lambda*dt)/(dy**2)
          end if
          
       end do
    end do

    do i=1,Nx
       k=(j-1)*Nx+1
       V(k
    end do

    do j=1,Ny
       k=(j-1)*Nx+i
       V(k) = V(k) + min(Tad, 2*Text+Tad*t)*(lambda*dt)/(dx**2)
    end do

  end subroutine cd_limite

  subroutine eq_arrhenius()
    integer  :: i

    eta0 = eta
    
    do i = 1,Nx*Ny
       coef   = -Ea/(R*U0(i))
       !if (coef<-700) coef=-700
       
       eta(i) = (eta(i)+dt*k0*exp(coef)) / (1+k0*dt*exp(coef))
       
       chi(i) = rho * Q * (eta(i)-eta0(i)) 
    end do

  end subroutine eq_arrhenius
end module mod_physique
