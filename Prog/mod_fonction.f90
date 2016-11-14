module mod_fonction
  use donnees

  implicit none

contains

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
          write(s,'(F20.10)')         U((j-1)*Nx+i)
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

    else if ((Num>=1000).and.(Num<10000)) then
       F_NAME='fichier/T'
       write(F_NAME (10:13),'(I4)') Num
       F_NAME(14:17)= '.dat'

    else if ((Num>=10000).and.(Num<100000)) then
       F_NAME='fichier/T'
       write(F_NAME (10:14),'(I5)') Num
       F_NAME(15:18)= '.dat'
    end if

    open(unit=2, file=F_NAME, action="write")

    do i=1,Nx
       do j=1,Ny
          write(2,*) i*dx, j*dy, U((j-1)*Nx+i), eta((j-1)*Nx+i)
       end do
       write(2,*)
    end do

    close(2)

    print*, 'SVG numero: ',Num

  end subroutine printvector


  !*******************************************************!
  !*******************************************************!
  !*******************************************************!
  !*******************************************************!


  subroutine Gradient_conjugue(X,b,eps)
    implicit none
    real(PR),dimension(Nx*Ny) :: x
    real(PR),dimension(Nx*Ny) :: b
    real(PR)                  :: eps

    real(PR),dimension(Nx*Ny) :: residu,p,Ap
    real(PR)                  :: norm_prec_c,norm_c,norm_rhs,App,err
    integer                   :: n,nb_iter

    ! Initialisation des variables
    n=Nx*Ny

    residu=b-Matmula(x)
    norm_c = dot_product(residu,residu)
    norm_prec_c=norm_c
    norm_rhs = dot_product(rhs,rhs)
    p=residu

    err = eps*norm_rhs !norm_c/norm_rhs >= eps

    ! Boucle
    nb_iter=0
    do while ((norm_c >=  err).and.(nb_iter<10000))

       Ap = Matmula(p)
       App = dot_product(Ap,p)
       x=x+(norm_c/App)*p
       residu=residu-(norm_c/App)*Ap
       norm_c = dot_product(residu,residu)
       p=residu+(norm_c/norm_prec_c)*p

       norm_prec_c=norm_c
       nb_iter=nb_iter+1
       iteration = iteration + 1
    end do

  end subroutine Gradient_conjugue


  function Matmula(U)
    implicit none
    real(PR),dimension(Nx*Ny) :: U,Matmula
    integer:: i,j,k

    Matmula=0
    do j= 1,Ny
       do i=1,Nx
          k=(j-1)*Nx+i
          ! Cy gauche
          if(j/=1)then
             Matmula(k)=Matmula(k)+Cy(k-Nx)*U(k-Nx)  !!
          end if
          ! Cx gauche
          if(i/=1)then
             Matmula(k)=Matmula(k)+Cx(k-j)*U(k-1)    !!
          end if
          ! Cd
          Matmula(k)=Matmula(k)+Cd(k)*U(k)
          ! Cx droite
          if(i/=Nx)then
             Matmula(k)=Matmula(k)+Cx(k-j+1)*U(k+1)    !!
          end if
          ! Cy droite
          if(j/=Ny)then
             Matmula(k)=Matmula(k)+Cy(k)*U(k+Nx)  !!
          end if
       end do
    end do
  end function Matmula

  function interp(Tab,val)
    ! Interpole la valeur (col 2) de tab
    ! en fonction de la la valeur voulue (col 1)
    implicit none
    real(PR), dimension(:,:) :: Tab
    real(PR)                 :: val, interp
    integer                  :: i,taille

    taille=size(Tab,1)

    if    (val < Tab(1,1)     )then
       interp = Tab(1,2)
    elseif(val >= Tab(taille,1))then
       interp = Tab(taille,2)
    else
       do i=1,taille-1
          if(val<Tab(i+1,1))then
             interp = (val - Tab(i,1)) / (Tab(i+1,1)-Tab(i,1))
             interp = interp * (Tab(i+1,2)-Tab(i,2)) + Tab(i,2)
             exit
          end if
       end do
    end if

  end function interp

end module mod_fonction
