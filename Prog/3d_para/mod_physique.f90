!> Gere la physique du modele

module mod_physique
  use donnees
  use mod_fonction
  implicit none

contains

  !> @brief Creer la matrice de discretisation
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
        do k=k1,kN
           num = (k-1)*Nx*Ny + (j-1)*Nx + i

          !! f_lambda renvoie le moyenne en lambda entre les deux mailles
          Cd(num) = rhocp(num) !! + (2.*lambda*dt)/(dx**2) + (2.*lambda*dt)/(dy**2)
          if(i/=1)  Cd(num) = Cd(num) + dt/(dx**2)*f_lambda(num-1,num)
          if(i/=Nx) Cd(num) = Cd(num) + dt/(dx**2)*f_lambda(num,num+1)
          if(j/=1)  Cd(num) = Cd(num) + dt/(dy**2)*f_lambda(num,num)!-Nx)
          if(j/=Ny) Cd(num) = Cd(num) + dt/(dy**2)*f_lambda(num,num)!+Nx)
          if(k/=1)  Cd(num) = Cd(num) + dt/(dz**2)*f_lambda(num,num)!-Nx*Ny)
          if(k/=Nz) Cd(num) = Cd(num) + dt/(dz**2)*f_lambda(num,num)!+Nx*Ny)

          if(i/=Nx)then
            Cx(num-(j-1)-(k-1)*Ny)  = -dt/(dx**2)*f_lambda(num,num+1)
          end if

          if(j/=Ny)then
            Cy(num-(k-1)*Nx) = -dt/(dy**2)*f_lambda(num,num)!+Nx)
          end if

          if(k/=Nz)then
            Cz(num) = -dt/(dz**2)*f_lambda(num,num)!+Nx*Ny)
          end if
        end do
      end do
    enddo
  end subroutine creation_matrice

  !> @brief Construit les conditions limites
  function cd_lim()
    integer :: i,j,k, num
    real(PR),dimension(num1:numN) :: cd_lim

    cd_lim = 0.

    ! ! Conditon Bas / Haut
    ! do k=1,Ny
    !   do i=1,Nx
    !     num=(k-1)*Nx*Ny + i !bas
    !     ! cd_lim(num) = cd_lim(num) + Cy(num-(k-1)*Nx)*h*dy*(U(num)-Text)/lambda(num)
    !
    !     num=(k-1)*Nx*Ny + (Ny-1)*Nx+i !haut
    !     ! cd_lim(num) = cd_lim(num) + Cy(num-(k-1)*Nx)*h*dy*(U(num)-Text)/lambda(num)
    !   end do
    ! end do
    !
    ! ! Condition Gauche / Droite
    ! do k=1,Ny
    !   do j=1,Ny
    !     num=(k-1)*Nx*Ny + (j-1)*Nx + 1 !gauche
    !     ! cd_lim(num) = cd_lim(num) + Cx(num-(j-1)-(k-1)*Ny)*h*dx*(U(num)-Text)/lambda(num)
    !
    !     num=(k-1)*Nx*Ny + j*Nx !droite
    !     ! cd_lim(num) = cd_lim(num) + Cx(num-(j-1)-(k-1)*Ny)*h*dx*(U(num)-Text)/lambda(num)
    !   end do
    ! end do
    !
    ! ! Condition Devant / Derriere
    ! do i=1,Nx
    !   do j=1,Ny
    !     num=(j-1)*Nx + i ! Devant
    !     ! cd_lim(num) = cd_lim(num) + Cz(num)*h*dx*(U(num)-Text)/lambda(num)
    !
    !     num=Nx*Ny*(Nz-1) + (j-1)*Nx + i ! Derriere
    !     ! cd_lim(num) = cd_lim(num) + Cz(num)*h*dx*(U(num)-Text)/lambda(num)
    !   end do
    ! end do
  end function cd_lim

  !> @brief Implique une condition de chauffe
  function chauffage()
    real(PR),dimension(num1:numN) :: chauffage
    integer::i,j,k

    !if(eta((Ny/2)*Nx+2)==1.)then
      !flux = 0.
    !end if

    chauffage = 0.
    !chauffage((Ny/2)*Nx+1) = -Cx((Ny/2)*Nx+1)*flux*dx/lambda((Ny/2)*Nx+1)
    ! do i=1,Nx
    !   do j=1,Ny
    !     ! Répartition du chauffage à gauche -> 100% au milieu, 0% en y=0 et y=Ly
    !     ! Fonction quadratique -(y/Ly)**2+(y/Ly)
    !     !2D  face xy
    !     chauffage((j-1)*Nx + i) = -Cx((j-1)*Nx+i)*flux*dx/lambda((j-1)*Nx+i) &
    !     *16*(j*dy)*(j*dy-Ly)*(i*dx)*(i*dx-Lx)/(Lx**2*Ly**2)
    !     ! Fonction linéaire -|y-Ly/2|/(Ly/2) + 1 (valeur absolue)
    !     !2D chauffage((j-1)*Nx+1) = -Cx((j-1)*Nx+1)*flux*(-2*abs(j*dy-Ly/2)/Ly+1)*dx/lambda((j-1)*Nx+1)
    !   end do
    ! end do
    do i=1,Nx
      do k=k1,kN !!!!===!!! ATTENTION kN NORMALEMENT
        !2D  face xz
        ! chauffage((k-1)*Nx*Ny + i) = -Cx((k-1)*Nx*Ny + i)*flux*dx/lambda((k-1)*Nx*Ny + i) &
        ! *16*(k*dz)*(k*dz-Lz)*(i*dx)*(i*dx-Lx)/(Lx**2*Lz**2)
        ! chauffe pleine
        !if(Me==Np-1) print*,((k-1)*Nx*Ny + i),num1,numN
        chauffage((k-1)*Nx*Ny + i) = -Cx((k-1)*Nx*Ny + i)*flux*dx/lambda((k-1)*Nx*Ny + i)

      end do
    end do

  end function chauffage

  !> @brief Calcule le terme d'energie de la reaction chimique
  function eq_arrhenius()
    integer :: num
    !integer :: p
    real(PR),dimension(Num1:NumN)::eq_arrhenius

    !p = 1
    eta0 = eta

    do Num = Num1,NumN
      coef   = -Ea/(R*U(Num))
      !if (coef<-700) coef=-700

      !implicite
      eta(Num) = (eta(Num)+dt*k0*exp(coef)) / (1+k0*dt*exp(coef))

      !explicite, ordre p
      !eta(i) = eta(i) + dt*k0*exp(-Ea/(R*U(i)))*(1-eta(i))**p

      eq_arrhenius(Num) = eta(Num)-eta0(Num)
    end do

  end function eq_arrhenius


  !> @brief Calcul des differentes proprietes des materiaux
  subroutine c_rho()
    ! Calule rhocp pur chaque maille
    integer  :: num

    do num = num1,numN
      rho(num) = fraction_vol(num,1) * ((1.-eta(num))*rho_Si + eta(num)*rho_Si3N4) &
      + fraction_vol(num,2) * rho_N2     &
      + fraction_vol(num,3) * rho_fibre
    end do
  end subroutine c_rho

  !> @brief Calule rhocp pour chaque maille
  subroutine c_rhocp()
    integer  :: num

    do num = num1,numN
      rhocp(num) = fraction_vol(num,1)* &
      ((1.-eta(num))*rho_Si*interp(cp_Si,U(num)) + eta(num)*rho_Si3N4*interp(cp_Si3N4,U(num))) &
      + fraction_vol(num,2) * rho_N2    * interp(cp_N2,U(num)) &
      + fraction_vol(num,3) * rho_fibre * interp(cp_fibre,U(num))
    end do
  end subroutine c_rhocp

  !> @brief Calul lambda pour chaque maille (isotrope pour le moment)
  subroutine c_lambda()
    integer  :: num

    do num = num1,numN
      lambda(num) = fraction_vol(num,1)/ &
      ((1.-eta(num))*interp(lambda_Si,U(num)) + eta(num)*interp(lambda_Si3N4,U(num))) &
      + fraction_vol(num,2)/interp(lambda_N2,U(num)) &
      + fraction_vol(num,3)/interp(lambda_fibre,U(num))
      lambda(num) = 1./lambda(num)
    end do
  end subroutine c_lambda

  !> @brief Calcul du flux en lambda entre caseles mailles i et j
  function f_lambda(i,j)
    integer  :: i,j
    real(PR) :: f_lambda

    f_lambda = 2. / (1./lambda(i) + 1./lambda(j))
  end function f_lambda

end module mod_physique
