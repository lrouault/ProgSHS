!> Gere l'initialisation du probleme.

module mod_preprocess
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR,C_INT32_T,C_FLOAT
  use donnees
  use mod_fonction

  implicit none

contains

  !> @brief Initialise toutes les variables
  subroutine initialisation(filename)
    character(len=*), intent(in) :: filename
    character(len=3)             :: bfr ! Variable poubelle
    integer                      :: i,j,k , num

    open(11, file=filename, action="read", status="old")
    read(11,'(A3,I6)')    bfr, Nx     ! "Nx="
    read(11,'(A3,I6)')    bfr, Ny     ! "Ny="
    read(11,'(A3,I6)')    bfr, Nz     ! "Ny="
    read(11,'(A3,F4.6)')  bfr, Lx     ! "Lx="
    read(11,'(A3,F4.6)')  bfr, Ly     ! "Ly="
    read(11,'(A3,F4.6)')  bfr, Lz     ! "Ly="
    !read(11,'(A4,F10.6)') bfr, rho    ! "rho="
    !read(11,'(A3,F10.6)') bfr, cp     ! "cp="
    !read(11,'(A7,F4.6)')  bfr, lambda ! "lambda"

    tmax    = 500.
    Niter   = 10000
    epsilon = 1.e-7
    epsilon = epsilon**2
    Text    = 298.
    Tad     = 850 !2300.0
    h       = 10.0d+0
    Ea      = 372.e3 !266547.
    R       = 8.3144621d+0
    k0      = 1.227e7 !6.2e17 !2.0e+4
    Q       = 372.e3 !287.e+3 !626e+3

    flux = 11.e+6

    dt = 0.5
    dx = Lx/(Nx+1)
    dy = Ly/(Ny+1)
    dz = Lz/(Nz+1)

    allocate(U(Nx*Ny*Nz), rhs(Nx*Ny*Nz), U0(Nx*Ny*Nz), eta(Nx*Ny*Nz))
    allocate(Cd(Nx*Ny*Nz), Cx((Nx-1)*Ny*Nz), Cy(Nx*(Ny-1)*Nz), Cz(Nx*Ny*(Nz-1)))

    rhs = 0.
    eta = 0.

    ! Conditions initiales

    U0 = 298.
    U=U0

    nb_fichiers = 100

    ! MATERIAU
    allocate(rho(Nx*Ny*Nz),   rhocp(Nx*Ny*Nz),    lambda(Nx*Ny*Nz))
    allocate(cp_Si(1,2)    , cp_Si3N4(1,2)    , cp_N2(1,2)    , cp_fibre(1,2)    )
    allocate(lambda_Si(1,2), lambda_Si3N4(1,2), lambda_N2(1,2), lambda_fibre(1,2))

    rho_Si    = 1600.
    cp_Si    = reshape((/298., 228./),(/1,2/))
    lambda_Si    = reshape((/298., 22./),(/1,2/))

    rho_Si3N4 = 1600.
    cp_Si3N4 = reshape((/298., 228./),(/1,2/))
    lambda_Si3N4 = reshape((/298., 22./),(/1,2/))

    rho_N2    = 1600.
    cp_N2    = reshape((/298., 228./),(/1,2/))
    lambda_N2    = reshape((/298., 22./),(/1,2/))

    ! rho_fibre = 1600.
    ! cp_fibre = reshape((/298., 228./),(/1,2/))
    ! lambda_fibre = reshape((/298., 22./),(/1,2/))
    ! Laine de verre
    rho_fibre = 1000.
    cp_fibre = reshape((/298., 500./),(/1,2/))
    lambda_fibre = reshape((/298., 10./),(/1,2/))
    ! ! Laine de verre
    ! rho_fibre = 35.
    ! cp_fibre = reshape((/298., 1030./),(/1,2/))
    ! lambda_fibre = reshape((/298., 0.039/),(/1,2/))




    ! allocate(cp_Si(10,2)    , cp_Si3N4(1,2)    , cp_N2(1,2)    , cp_fibre(1,2)    )
    ! allocate(lambda_Si(1,2), lambda_Si3N4(1,2), lambda_N2(1,2), lambda_fibre(1,2))
    !
    ! rho_Si    = 1600.
    ! cp_Si    = reshape((/300., 228. , 300., 228. ,500., 228. , 700., 228. , 900., 228. , 1100., 228. , 1300., 228. , 1500., 228. , 1700., 228. , 1900., 228. , /),(/1,2/),(/0./),(/2,1/))
    ! lambda_Si    = reshape((/298., 22./),(/1,2/))
    !
    ! rho_Si3N4 = 1600.
    ! cp_Si3N4 = reshape((/298., 228./),(/1,2/))
    ! lambda_Si3N4 = reshape((/298., 22./),(/1,2/))
    !
    ! rho_N2    = 1600.
    ! cp_N2    = reshape((/298., 228./),(/1,2/))
    ! lambda_N2    = reshape((/298., 22./),(/1,2/))
    !
    ! rho_fibre = 1000.
    ! cp_fibre = reshape((/298., 500./),(/1,2/))
    ! lambda_fibre = reshape((/298., 10./),(/1,2/))





    allocate(fraction_vol(Nx*Ny*Nz,3))! (/Si,N2,fibre/)
    do i = 1,Nx
      do j = 1,Ny
        do k = 1,Nz
          num=(k-1)*Nx*Ny + (j-1)*Nx + i
          ! ! Partie superieur en fibre
          ! if (i<=Nx/2)then
          !   fraction_vol(num,3) = 1.
          ! else
          !   fraction_vol(num,3) = 0.
          ! end if
          !
          ! ! 5 bandes
          ! if (i<=Nx/5)then
          !   fraction_vol(num,3) = 1.
          ! elseif (i<=2*Nx/5)then
          !   fraction_vol(num,3) = 0.
          ! elseif (i<=3*Nx/5)then
          !   fraction_vol(num,3) = 1.
          ! elseif (i<=4*Nx/5)then
          !   fraction_vol(num,3) = 0.
          ! else
          !   fraction_vol(num,3) = 1.
          ! end if
          !
          ! 5 bandes *2
          if (i<=Nx/5 .and. k<=Nz/2)then
            fraction_vol(num,3) = 1.
          elseif (i<=2*Nx/5.and. k<=Nz/2)then
            fraction_vol(num,3) = 0.
          elseif (i<=3*Nx/5.and. k<=Nz/2)then
            fraction_vol(num,3) = 1.
          elseif (i<=4*Nx/5.and. k<=Nz/2)then
            fraction_vol(num,3) = 0.
          elseif (i<=Nx.and. k<=Nz/2)then
            fraction_vol(num,3) = 1.
          elseif (i<=Nx/5)then
            fraction_vol(num,3) = 0.
          elseif (i<=2*Nx/5)then
            fraction_vol(num,3) = 1.
          elseif (i<=3*Nx/5)then
            fraction_vol(num,3) = 0.
          elseif (i<=4*Nx/5)then
            fraction_vol(num,3) = 1.
          else
            fraction_vol(num,3) = 0.
          end if


          ! Fibre
          ! fraction_vol(num,3) = Poro(i,j,k)/135.

          ! 256 Pour image crop et 135 pour tex1_
        end do
      end do
    end do

    fraction_vol(:,1) = 0.76*(1-fraction_vol(:,3)) ! 0.76 CFC
    fraction_vol(:,2) = 1.-fraction_vol(:,1)-fraction_vol(:,3)


  end subroutine initialisation

!*******************************************************!
!*******************************************************!

  !> @brief Désalloue les différents tableaux
  subroutine fin()

    deallocate(U, rhs, U0, eta)
    deallocate(rho,   rhocp,    lambda)
    deallocate(cp_Si,     cp_Si3N4 ,    cp_N2 ,    cp_fibre )
    deallocate(lambda_Si, lambda_Si3N4, lambda_N2, lambda_fibre)
    deallocate(fraction_vol) ! (/Si,N2,fibre/)

  end subroutine fin

  !*******************************************************!
  !*******************************************************!

  !> @brief Remplit la porosite du domaine (fibres)
  subroutine fillPoro(F_NAME,F_NAME2)
    character(len=*),intent(in)                   :: F_NAME,F_NAME2
    !CHARACTER(KIND=C_CHAR)  :: data2
    character(kind=C_CHAR)      :: data2
    integer :: i,j,k

    ! syntaxe complète
    open( 101 , file=F_NAME, form="formatted", action="read")
    read(101,*) porox,poroy,poroz
    close(101)
    allocate(Poro(porox,poroy,poroz))

    open( 102 , file=F_NAME2, form="unformatted", action="read", &
    access="stream")
    ! 1	FORMAT ( Z4 ) ! Format pour ecrire les caractere en hexadecimal
    do k=1,poroz
      do j=1,poroy
        do i=1,porox
          read(102) data2
          Poro(i,j,k) = ichar(data2)
        enddo
      enddo
    enddo
    close(102)
    print*, "Remplissage des porosite OK"
  end subroutine fillPoro

  subroutine fillPoro2(F_NAME,sizex,sizey,sizez)
    character(len=*),intent(in)                   :: F_NAME
    character(kind=C_CHAR)      :: data2
    integer,intent(in)         ::sizex,sizey,sizez
    integer :: i,j,k

    ! syntaxe complète
    porox=sizex
    poroy=sizey
    poroz=sizez
    allocate(Poro(porox,poroy,poroz))

    open( 102 , file=F_NAME, form="unformatted", action="read", &
    access="stream")
    ! 1	FORMAT ( Z4 ) ! Format pour ecrire les caractere en hexadecimal
    do k=1,poroz
      do j=1,poroy
        do i=1,porox
          read(102) data2
          Poro(i,j,k) = ichar(data2)
        enddo
      enddo
    enddo
    close(102)
    print*, "Remplissage des porosite OK"
  end subroutine fillPoro2

  !> @brief Remplit l'orientation des fibres
  subroutine fillOrientation(F_NAME)
    character(len=*),intent(in)                   :: F_NAME
    integer(KIND=C_INT32_T) :: data_int
    real(kind=C_FLOAT)      :: data_float

    integer :: i,j,k

    ! syntaxe complète
    open( 102 , file=F_NAME, form="unformatted", action="read", &
    access="stream")
    1	FORMAT ( Z4 )

    read(102) data_int
    if ( porox/=data_int ) print*, "ERREUR DIM X ORIENTATION"
    read(102) data_int
    if ( poroy/=data_int ) print*, "ERREUR DIM Y ORIENTATION"
    read(102) data_int
    if ( poroz/=data_int ) print*, "ERREUR DIM Z ORIENTATION"

    allocate(Orien(porox,poroy,poroz,3))

    do k=1,poroz
      do j=1,poroy
        do i=1,porox
          read(102) data_float
          Orien(i,j,k,1) = data_float
          read(102) data_float
          Orien(i,j,k,2) = data_float
          read(102) data_float
          Orien(i,j,k,3) = data_float
        enddo
      enddo
    enddo
    close(102)
    print*, "Remplissage des orientations OK"
  end subroutine fillOrientation


end module mod_preprocess
