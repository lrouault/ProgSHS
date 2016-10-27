module mod_preprocess
  use donnees
  use mod_fonction

  implicit none

contains

  subroutine initialisation(filename)
    character(len=*), intent(in) :: filename
    character(len=3)             :: bfr ! Variable poubelle
    integer                      :: i

    open(11, file=filename, action="read", status="old")
    read(11,'(A3,I6)')    bfr, Nx     ! "Nx="
    read(11,'(A3,I6)')    bfr, Ny     ! "Ny="
    read(11,'(A3,F4.6)')  bfr, Lx     ! "Lx="
    read(11,'(A3,F4.6)')  bfr, Ly     ! "Ly="
    read(11,'(A4,F10.6)') bfr, rho    ! "rho="
    read(11,'(A3,F10.6)') bfr, cp     ! "cp="
    read(11,'(A7,F4.6)')  bfr, lambda ! "lambda"

    N       = Nx*Ny
    Tmax    = 1000.
    Niter   = 10000
    Nmax    = n+1
    epsilon = 1.e-7
    Text    = 298.
    Tad     = 2300.0
    h       = 10.0d+0
    Ea      = 266547.
    R       = 8.3144621d+0
    k0      = 2.0e+4
    Q       = 626.e+3 !!626e+3

    dt = real(Tmax)/Niter
    dx = Lx/(Nx+1)
    dy = Ly/(Ny+1)

    allocate(U(N), V(N), U0(N), eta(N), chi(N))
    allocate(alpha(Nx*Ny), beta((Nx-1)*Ny), gamma(Nx*(Ny-1)))

    V   = 0
    eta = 0

    ! Conditions initiales

    U0 = 1000
    !do i = 1,Ny
    !   U0(bij(1,i,Ny)) = (3.0/4) * Tad
    !end do
    U=U0

  end subroutine initialisation

!*******************************************************!
!*******************************************************!

  subroutine fin()

    deallocate(U, V, U0, eta, chi)

  end subroutine fin

end module mod_preprocess
