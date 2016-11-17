module mod_physique
  use donnees
  use mod_fonction
  implicit none

contains

  subroutine creation_matrice()
    !Cd = rho*cp + (2.*lambda*dt)/(dx**2) + (2.*lambda*dt)/(dy**2)
    !Cx = -(lambda*dt)/(dx**2)
    !Cy = -(lambda*dt)/(dy**2)
    integer :: i,j,k
    integer :: num

    call c_rho()
    call c_rhocp()
    call c_lambda()

    do i=1,Nx
      do j=1,Ny
        do k=1,Nz
          num=(k-1)*Nx*Ny + (j-1)*Nx + i

          Cd(num) = rhocp(num) !! + (2.*lambda*dt)/(dx**2) + (2.*lambda*dt)/(dy**2)
          if(i/=1)  Cd(num) = Cd(num) + dt/(dx**2)*f_lambda(num-1,num)
          if(i/=Nx) Cd(num) = Cd(num) + dt/(dx**2)*f_lambda(num,num+1)
          if(j/=1)  Cd(num) = Cd(num) + dt/(dy**2)*f_lambda(num-Nx,num)
          if(j/=Ny) Cd(num) = Cd(num) + dt/(dy**2)*f_lambda(num,num+Nx)
          if(k/=1)  Cd(num) = Cd(num) + dt/(dz**2)*f_lambda(num-Nx*Ny,num)
          if(k/=Nz) Cd(num) = Cd(num) + dt/(dz**2)*f_lambda(num,num+Nx*Ny)

          if(i/=Nx)then
            Cx(num-j+1)  = -dt/(dx**2)*f_lambda(num,num+1)
          end if

          if(j/=Ny)then
            Cy(num) = -dt/(dy**2)*f_lambda(num,num+Nx)
          end if

          if(k/=Nz)then
            Cz(num) = -dt/(dz**2)*f_lambda(num,num+Nx*Ny)
          end if
        end do
      end do
    enddo
  end subroutine creation_matrice

  function cd_lim()
    integer :: i,j,k, num
    real(PR),dimension(Nx*Ny) :: cd_lim

    cd_lim = 0.

    do z=
    do i=1,Nx
      num=i !bas
      num=(k-1)*Nx*Ny + (j-1)*Nx + i
      cd_lim(num) = cd_lim(num) + Cy(num)*h*dy*(U(num)-Text)/lambda(num)

      num=(Ny-1)*Nx+i !haut
      cd_lim(num) = cd_lim(num) + Cy(num)*h*dy*(U(num)-Text)/lambda(num)
    end do

    do j=1,Ny
      num=(j-1)*Nx+1 !gauche
      cd_lim(num) = cd_lim(num) + Cx(num)*h*dx*(U(num)-Text)/lambda(num)

      num=(j-1)*Nx+Nx !droite
      cd_lim(num) = cd_lim(num) + Cx(num)*h*dx*(U(num)-Text)/lambda(num)
    end do

  end function cd_lim

  function chauffage()
    real(PR),dimension(Nx*Ny) :: chauffage
    integer::j

    if(eta((Ny/2)*Nx+2)==1.)then
      flux = 0.
    end if

    chauffage = 0.
    !chauffage((Ny/2)*Nx+1) = -Cx((Ny/2)*Nx+1)*flux*dx/lambda((Ny/2)*Nx+1)
    do j=1,Ny
      ! Répartition du chauffage à gauche -> 100% au milieu, 0% en y=0 et y=Ly
      ! Fonction quadratique -(y/Ly)**2+(y/Ly)
      !chauffage((j-1)*Nx+1) = -Cx((j-1)*Nx+1)*flux*(-(j*dy/Ly)**2+j*dy/Ly)*dx/lambda((j-1)*Nx+1)

      ! Fonction linéaire -|y-Ly/2|/(Ly/2) + 1 (valeur absolue)
      chauffage((j-1)*Nx+1) = -Cx((j-1)*Nx+1)*flux*(-2*abs(j*dy-Ly/2)/Ly+1)*dx/lambda((j-1)*Nx+1)
    end do

  end function chauffage

  function eq_arrhenius()
    integer :: i
    !integer :: p
    real(PR),dimension(Nx*Ny)::eq_arrhenius

    !p = 1
    eta0 = eta

    do i = 1,Nx*Ny
      coef   = -Ea/(R*U(i))
      !if (coef<-700) coef=-700

      !implicite
      eta(i) = (eta(i)+dt*k0*exp(coef)) / (1+k0*dt*exp(coef))

      !explicite, ordre p
      !eta(i) = eta(i) + dt*k0*exp(-Ea/(R*U(i)))*(1-eta(i))**p

      eq_arrhenius(i) = eta(i)-eta0(i)
    end do

  end function eq_arrhenius


  ! Calcul des differentes proprietes des materiaux
  subroutine c_rho()
    ! Calule rhocp pour chaque maille
    integer  :: i

    do i = 1,Nx*Ny
      rho(i) = fraction_vol(i,1) * ((1.-eta(i))*rho_Si + eta(i)*rho_Si3N4) &
      + fraction_vol(i,2) * rho_N2     &
      + fraction_vol(i,3) * rho_fibre
    end do
  end subroutine c_rho


  subroutine c_rhocp()
    ! Calule rhocp pour chaque maille
    integer  :: i

    do i = 1,Nx*Ny
      rhocp(i) = fraction_vol(i,1)* &
      ((1.-eta(i))*rho_Si*interp(cp_Si,U(i)) + eta(i)*rho_Si3N4*interp(cp_Si3N4,U(i))) &
      + fraction_vol(i,2) * rho_N2    * interp(cp_N2,U(i)) &
      + fraction_vol(i,3) * rho_fibre * interp(cp_fibre,U(i))
    end do
  end subroutine c_rhocp


  subroutine c_lambda()
    ! Calule lambda pour chaque maille (isotrope pour le moment)
    integer  :: i

    do i = 1,Nx*Ny
      lambda(i) = fraction_vol(i,1)/ &
      ((1.-eta(i))*interp(lambda_Si,U(i)) + eta(i)*interp(lambda_Si3N4,U(i))) &
      + fraction_vol(i,2)/interp(lambda_N2,U(i)) &
      + fraction_vol(i,3)/interp(lambda_fibre,U(i))
      lambda(i) = 1./lambda(i)
    end do
  end subroutine c_lambda


  function f_lambda(i,j)
    ! Calcul du flux en lambda entre caseles mailles i et j
    integer  :: i,j
    real(PR) :: f_lambda

    f_lambda = 1. / (1./lambda(i) + 1./lambda(j))
  end function f_lambda

end module mod_physique
