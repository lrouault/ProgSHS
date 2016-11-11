module mod_physique
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR
  use donnees
  use mod_fonction
  implicit none

contains

  subroutine creation_matrice()
    !Cd = rho*cp + (2.*lambda*dt)/(dx**2) + (2.*lambda*dt)/(dy**2)
    !Cx = -(lambda*dt)/(dx**2)
    !Cy = -(lambda*dt)/(dy**2)
    integer :: i,j,k

    call c_rho()
    call c_rhocp()
    call c_lambda()

    do i=1,Nx
      do j=1,Ny
        k=(j-1)*Nx+i

        Cd(k) = rhocp(k) !!  + (2.*lambda*dt)/(dx**2) + (2.*lambda*dt)/(dy**2)
        if(i/=1)  Cd(k) = Cd(k) + dt/(dx**2)*f_lambda(k-1,k)
        if(i/=Nx) Cd(k) = Cd(k) + dt/(dx**2)*f_lambda(k,k+1)
        if(j/=1)  Cd(k) = Cd(k) + dt/(dy**2)*f_lambda(k-Nx,k)
        if(j/=Ny) Cd(k) = Cd(k) + dt/(dy**2)*f_lambda(k,k+Nx)


        if(i/=Nx)then
          Cx(k-j+1)  = -dt/(dx**2)*f_lambda(k,k+1)
        end if

        if(j/=Ny)then
          Cy(k) = -dt/(dy**2)*f_lambda(k,k+Nx)
        end if
      end do
    end do
  end subroutine creation_matrice

  function cd_lim()
    integer :: i,j,k
    real(PR),dimension(Nx*Ny) :: cd_lim

    cd_lim = 0.

    !do i=1,Nx
    !   k=i !bas
    !   cd_lim(k) = cd_lim(k) - Cy(1)*(U(k)-h*dy*(U(k)-Text)/lambda)

    !   k=(Ny-1)*Nx+i !haut
    !   cd_lim(k) = cd_lim(k) - Cy(1)*(U(k)-h*dy*(U(k)-Text)/lambda)
    !  end do

    !do j=1,Ny
    !   k=(j-1)*Nx+1 !gauche
    !   cd_lim(k) = cd_lim(k) - Cx(1)*(U(k)-h*dx*(U(k)-Text)/lambda)

    !     k=(j-1)*Nx+Nx !droite
    !   cd_lim(k) = cd_lim(k) - Cx(1)*(U(k)-h*dx*(U(k)-Text)/lambda)
    !end do

    cd_lim((Ny/2)*Nx+1) = cd_lim((Ny/2)*Nx+1) &
    - Cx((Ny/2)*Nx+1)*flux*dx/lambda((Ny/2)*Nx+1)

  end function cd_lim

  function eq_arrhenius()
    integer  :: i
    real(PR),dimension(Nx*Ny)::eq_arrhenius

    eta0 = eta

    do i = 1,Nx*Ny
      coef   = -Ea/(R*U(i))
      !if (coef<-700) coef=-700

      eta(i) = (eta(i)+dt*k0*exp(coef)) / (1+k0*dt*exp(coef))

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


  subroutine fillFibre(F_NAME,F_NAME2)
    character(len=*),intent(in)                   :: F_NAME,F_NAME2
    CHARACTER(KIND=C_CHAR) :: data2
    integer :: i,j,k

    ! syntaxe complète
    open( 101 , file=F_NAME, form="formatted", action="read")
    read(101,*) porox,poroy,poroz
    close(101)
    allocate(Poro(porox,poroy,poroz))

    open( 102 , file=F_NAME2, form="unformatted", action="read", &
    access="stream")
    1	FORMAT ( Z4 )
    do k=1,poroz
      do j=1,poroy
        do i=1,porox
          read(102) data2
          Poro(i,j,k) = ichar(data2)
        enddo
      enddo
    enddo
    close(102)
  end subroutine fillFibre

end module mod_physique
