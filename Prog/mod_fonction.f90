module mod_fonction
  use donnees

  implicit none

contains

  function bij(i, j, Ny)

    integer,intent(in) :: i, j, Ny
    integer            :: bij

    bij=Ny*(i-1)+j

  end function bij

  !*******************************************************!
  !*******************************************************!

  subroutine gradc(a, b, c, epsilon, Niter, X, V)

    real(PR), intent(in)                      :: a, b, c, epsilon
    integer, intent(in)                     :: Niter
    real(PR), dimension(:), intent(in)        :: V
    real(PR), dimension(size(V)), intent(out) :: X
    real(PR), dimension(size(V))              :: R, D, W
    real(PR)                                  :: alpha, beta
    integer                                 :: iter

    X = 10
    call prod_mat_vec(a, b, c, X, R)
    R    = R - V
    D    = R
    iter = 1

    do while ((norme(R)>epsilon).AND.(iter<Niter))

       call prod_mat_vec(a, b, c, D, W)
       alpha = dot_product(D,R) / dot_product(D,W)
       X=X-alpha*D

       beta = 1 / (norme(R)**2)
       R    = R - alpha*W
       beta = beta * norme(R)**2

       D = R + beta*D
       iter = iter + 1
       !print*,'iter',iter
       
    end do

  end subroutine gradc

  !*******************************************************!
  !*******************************************************!

  subroutine prod_mat_vec(a, b, c, x, y)


    real(PR), intent(in)::a,b,c
    real(PR), intent(in), dimension(Nx*Ny)::x
    real(PR), intent(inout), dimension(Nx*Ny)::y
    integer::i,j
    
    ! 1er bloc
    y(bij(1,1,Ny)) = (a-b)*x(bij(1,1,Ny))+b*x(bij(1,2,Ny))+c*x(bij(2,1,Ny))
    do j = 2, Ny-1
       y(bij(1,j,Ny)) = b*x(bij(1,j-1,Ny))+a*x(bij(1,j,Ny))+b*x(bij(1,j+1,Ny))+c*x(bij(2,j,Ny))
    end do
    y(bij(1,Ny,Ny)) = b*x(bij(1,Ny-1,Ny))+(a-b)*x(bij(1,Ny,Ny))+c*x(bij(2,Ny,Ny))
    ! Milieu Matrice
    do i = 2, Nx-1
       y(bij(i,1,Ny)) = c*x(bij(i-1,1,Ny))+(a-b)*x(bij(i,1,Ny))+b*x(bij(i,2,Ny))+c*x(bij(i+1,1,Ny))
       do j = 2, Ny-1
          y(bij(i,j,Ny)) = c*x(bij(i-1,j,Ny))+b*x(bij(i,j-1,Ny))+a*x(bij(i,j,Ny))+b*x(bij(i,j+1,Ny))+c*x(bij(i+1,j,Ny))
       end do
       y(bij(i,Ny,Ny)) = c*x(bij(i-1,Ny,Ny))+b*x(bij(i,Ny-1,Ny))+(a-b)*x(bij(i,Ny,Ny))+c*x(bij(i+1,Ny,Ny))
    end do
    ! Dernier bloc
    y(bij(Nx,1,Ny)) = c*x(bij(Nx-1,1,Ny))+(a-b)*x(bij(Nx,1,Ny))+b*x(bij(Nx,2,Ny))
    do j = 2, Ny-1
       y(bij(Nx,j,Ny)) = c*x(bij(Nx-1,j,Ny))+b*x(bij(Nx,j-1,Ny))+a*x(bij(Nx,j,Ny))+b*x(bij(Nx,j+1,Ny))
    end do
    y(bij(Nx,Ny,Ny)) = c*x(bij(Nx-1,Ny,Ny))+b*x(bij(Nx,Ny-1,Ny))+(a-b)*x(bij(Nx,Ny,Ny))
  end subroutine prod_mat_vec

  !*******************************************************!
  !*******************************************************!

  function norme(R) ! norme euclidienne

    implicit none
    real(PR), intent(in), dimension(:) :: R
    real(PR)                           :: norme, S
    integer                          ::i

    S = 0
    do i = 1, size(R)
       S = S + R(i)*R(i)
    end do
    norme = sqrt(S)

  end function norme

  !*******************************************************!
  !*******************************************************!

  subroutine printmat(A,n,ch)

    integer,intent(in)             :: n
    real,dimension(n,n),intent(in) :: A
    character(len=*),intent(in)    :: ch
    integer                        :: i

    print*, ch
    do i=1,n
       print*, A(i,:)
    end do

  end subroutine printmat

  !*******************************************************!
  !*******************************************************!

  subroutine  writeVtk(U, Nx, Ny, dx, dy, N)

    integer,intent(in)                  :: Nx, Ny, N
    real(PR),dimension(Nx*Ny),intent(in) :: U
    real(PR), intent(in)                 :: dx, dy
    integer                             :: i, j
    character(len=20)                   :: F_NAME
    character(len=40)                   :: s1, s2, s3, s4, s5, s

    if (N<10) then
       F_NAME='fichier/T'
       write(F_NAME (10:10),'(I1)') N
       F_NAME(11:14)= '.vtk'

    elseif ((N>=10).and.(N<100)) then
       F_NAME='fichier/T'
       write(F_NAME (10:11),'(I2)') N
       F_NAME(12:15)= '.vtk'

    elseif ((N>=100).and.(N<1000)) then
       F_NAME='fichier/T'
       write(F_NAME (10:12),'(I3)') N
       F_NAME(13:16)= '.vtk'

    else if ((N>=1000).and.(N<10000)) then
       F_NAME='fichier/T'
       write(F_NAME (10:13),'(I4)') N
       F_NAME(14:17)= '.vtk'

    else if ((N>=10000).and.(N<100000)) then
       F_NAME='fichier/T'
       write(F_NAME (10:14),'(I5)') N
       F_NAME(15:18)= '.vtk'
    end if

    open(unit=2, file=F_NAME, action="write")

    write(2,'(a)')       "# vtk DataFile Version 2.0"
    write(2,'(a)')       "Titre"
    write(2,'(a)')       "ASCII"
    write(2,'(a)')       "DATASET STRUCTURED_POINTS"
    write(s1,'(I4)')     Nx
    write(s2,'(I4)')     Ny
    write(2,'(a)')       "DIMENSIONS " // adjustl(trim(s1)) // adjustl(trim(s2)) // "1"
    write(2,'(a)')       "ORIGIN 0.0 0.0 0.0"
    write(s3,'(F20.10)') dx
    write(s4,'(F20.10)') dy
    write(2,'(a)')       "SPACING " // adjustl(trim(s3)) // adjustl(trim(s4)) // "0.0"
    write(s5,'(I8)')     (Nx*Ny)
    write(2,'(a)')       "POINT_DATA " // adjustl(trim(s5))
    write(2,'(a)')       "SCALARS Temperature double"
    write(2,'(a)')       "LOOKUP_TABLE default"

    do i = 1,Ny
       do j = 1,Nx

          if(modulo(j,10) == 0) then
             write(2,'(/)',advance='no')
          end if
          s = ""
          write(s,'(F20.10)')         U(bij(i,j,Ny))
          write(2,'(a)',advance='no') adjustl(trim(s)) // " "

       end do
    end do

    close(2) 

  end subroutine writeVtk

  !*******************************************************!
  !*******************************************************!

  subroutine printvector(U, Num)

    integer,intent(in)                  ::  Num
    real(PR),dimension(Nx*Ny),intent(in) :: U
    integer                             :: i, j
    character(len=20)                   :: F_NAME

    if (Num<10) then
       F_NAME='fichier/T'
       write(F_NAME (10:10),'(I1)') Num
       F_NAME(11:14)= '.dat'

    elseif ((Num>=10).and.(Num<100)) then
       F_NAME='fichier/T'
       write(F_NAME (10:11),'(I2)') Num
       F_NAME(12:15)= '.dat'

    elseif ((Num>=100).and.(Num<1000)) then
       F_NAME='fichier/T'
       write(F_NAME (10:12),'(I3)') Num
       F_NAME(13:16)= '.dat'

    else if ((N>=1000).and.(N<10000)) then
       F_NAME='fichier/T'
       write(F_NAME (10:13),'(I4)') N
       F_NAME(14:17)= '.dat'

    else if ((Num>=10000).and.(Num<100000)) then
       F_NAME='fichier/T'
       write(F_NAME (10:14),'(I5)') Num
       F_NAME(15:18)= '.dat'
    end if

    open(unit=2, file=F_NAME, action="write")

    do i=1,Nx
       do j=1,Ny
          write(2,*) i*dx, j*dy, U(bij(i,j,Ny))
       end do
       write(2,*)
    end do

    close(2)

    print*, 'SVG numero: ',Num

  end subroutine printvector

end module mod_fonction
