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

    call c_rho()
    call c_rhocp()
    call c_lambda()
    do i=1,Nx
       do j=j1,jN
          k=(j-1)*Nx+i

          Cd(k) = rhocp(k) !!  + (2.*lambda*dt)/(dx**2) + (2.*lambda*dt)/(dy**2)
          if(i/=1)  Cd(k) = Cd(k) + dt/(dx**2)*f_lambda(k,k-1)
          if(i/=Nx) Cd(k) = Cd(k) + dt/(dx**2)*f_lambda(k,k+1)
          if(j/=1)  Cd(k) = Cd(k) + dt/(dy**2)*f_lambda(k,k)!-Nx)
          if(j/=Ny) Cd(k) = Cd(k) + dt/(dy**2)*f_lambda(k,k)!+Nx)


          if(i/=Nx) Cx(k-j+1) = -dt/(dx**2)*f_lambda(k,k+1)
          if(j/=Ny) Cy(k) = -dt/(dy**2)*f_lambda(k,k)!+Nx)

       end do
    end do
    
    if(Me/=Np-1)then
       call MPI_SEND(Cy(numN-Nx+1:numN),Nx,MPI_REAL8,Me+1,101,MPI_COMM_WORLD,statinfo)
    end if
    if(Me/=0)then
       call MPI_RECV(Cy(num1-Nx:num1-1),Nx,MPI_REAL8,Me-1,101,MPI_COMM_WORLD,status,statinfo)
    end if
    
  end subroutine creation_matrice

  function cd_lim()
    integer :: i,j,k
    real(PR),dimension(num1:numN) :: cd_lim

    cd_lim = 0.
    if(j1==1) then 
       do i = 1,Nx
          k = i !bas
          cd_lim(k) = cd_lim(k) + Cy(k)*h*dy*(U(k)-Text)/lambda(k)
       end do
    end if
    if (jN==Ny) then 
       do i = 1,Nx
          k = (Ny-1)*Nx + i !haut
          cd_lim(k) = cd_lim(k) + Cy(k-Nx)*h*dy*(U(k)-Text)/lambda(k)
       end do
    end if
    
    do j = j1,jN
       k = (j-1)*Nx + 1 !gauche
       cd_lim(k) = cd_lim(k) + Cx(k-j+1)*h*dx*(U(k)-Text)/lambda(k)
       
       k = (j-1)*Nx + Nx !droite
       cd_lim(k) = cd_lim(k) + Cx(k-j)*h*dx*(U(k)-Text)/lambda(k)
    end do

  end function cd_lim

  function chauffage()
    real(PR),dimension(num1:numN) :: chauffage
    integer  :: j,k
    real(PR) :: y
    real(PR) :: norm_eta

    !Arrêt du chauffage si la réaction est lancée
    call MPI_ALLREDUCE(dot_product(eta,eta),norm_eta,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,statinfo)
    !if(norm_eta>0.16 .and. flux/=0.) flux = 0.
    
    chauffage = 0.
    !chauffage((Ny/2)*Nx+1) = -Cx((Ny/2)*Nx+1)*flux*dx/lambda((Ny/2)*Nx+1)
    do j=j1,jN
       k = (j-1)*Nx + 1
       y = j*dy

       ! Répartition du chauffage à gauche -> 100% au milieu, 0% en y=0 et y=Ly
       ! Fonction quadratique (y/Ly)*(1-y/Ly)
       chauffage(k) = -Cx(k-j+1)*flux*(y*(1-y/Ly)/Ly)*dx/lambda(k)

       ! Fonction linéaire -|y-Ly/2|/(Ly/2) + 1 (valeur absolue)
       !chauffage(k) = -Cx(k-j+1)*flux*(-2*abs(y-Ly/2)/Ly+1)*dx/lambda(k)

       ! Chauffage uniforme
       !chauffage(k) = -Cx(k-j+1)*flux*dx/lambda(k)
    end do
    
  end function chauffage

  function eq_arrhenius()
    integer :: k
    !integer :: p
    real(PR),dimension(num1:numN)::eq_arrhenius

    !p = 1
    eta0 = eta

    do k = num1,numN
       coef   = -Ea/(R*U(k))
       if (coef<-700) coef=-700

       !implicite
       eta(k) = (eta(k)+dt*k0*exp(coef)) / (1+k0*dt*exp(coef))

       !explicite, ordre p
       !eta(k) = eta(k) + dt*k0*exp(-Ea/(R*U(k)))*(1-eta(k))**p

       eq_arrhenius(k) = eta(k)-eta0(k)
    end do

  end function eq_arrhenius


  ! Calcul des differentes proprietes des materiaux
  subroutine c_rho()
    ! Calule rhocp pour chaque maille
    integer  :: k
    
    do k = num1,numN
       rho(k) = fraction_vol(k,1) * ((1.-eta(k))*rho_Si + eta(k)*rho_Si3N4) &
            + fraction_vol(k,2) * rho_N2     &
            + fraction_vol(k,3) * rho_fibre
    end do
  end subroutine c_rho
  
  
  subroutine c_rhocp()
    ! Calule rhocp pour chaque maille
    integer  :: k
    
    do k = num1,numN
       rhocp(k) = fraction_vol(k,1) * & !fraction_vol(i,k)* &
            ((1.-eta(k))*rho_Si*interp(cp_Si,U(k)) + eta(k)*rho_Si3N4*interp(cp_Si3N4,U(k))) &
            + fraction_vol(k,2) * rho_N2    * interp(cp_N2,U(k)) &
            + fraction_vol(k,3) * rho_fibre * interp(cp_fibre,U(k))
    end do
  end subroutine c_rhocp

  
  subroutine c_lambda()
    ! Calule lambda pour chaque maille (isotrope pour le moment)
    integer  :: k

    do k = num1,numN
       lambda(k) = fraction_vol(k,1)/ &
            ((1.-eta(k))*interp(lambda_Si,U(k)) + eta(k)*interp(lambda_Si3N4,U(k))) &
            + fraction_vol(k,2)/interp(lambda_N2,U(k)) &
            + fraction_vol(k,3)/interp(lambda_fibre,U(k))
       lambda(k) = 1./lambda(k)
    end do
  end subroutine c_lambda
  

  function f_lambda(i,j)
    ! Calcul du flux en lambda entre caseles mailles i et j
    integer  :: i,j,tmp
    real(PR) :: f_lambda

    if(i<num1 .or. i>numN)then
       tmp = i
       i = j
       j = tmp
    end if       

    if(j>=num1 .and. j<=numN)then
       f_lambda = 2. / (1./lambda(i) + 1./lambda(j))
    else if(j<num1)then ! lambda proc-1
       
    else !  j>numN        lambda proc+1

    end if
    
  end function f_lambda
  
end module mod_physique
