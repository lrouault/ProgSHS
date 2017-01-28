!> Déclaration de fonctions de base utilisees par le programme

module mod_fonction
  use donnees

  implicit none

contains

  !*******************************************************!
  !*******************************************************!


  !> Calcul la norme d'un vecteur.
  !! @param[in] R
  !! @return La norme du vecteur
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

  !> @brief Permet d'afficher une matrice carree n*n ligne par ligne.
  !!
  !! L'impression de la matrice est precede par l'affichage du parametre
  !! d'entree ch. Il s'en suit une lecture de la matrice ligne aprés ligne sans
  !! formatage spécial.
  !! @param[in] A matrice de taille n*n
  !! @param[in] n taille de la matrice
  !! @param[in] ch chaine de caractere imprimer avant la matrice
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

  !> @brief Ecrit au format vtk ...
  !!
  !!..........
  !! @param[in] U
  !! @param[in] Nx,Ny
  !! @param[in] dx,dy
  !! @param[in] N
  subroutine  writeVtk(U, Nx, Ny, dx, dy, N)

    integer,intent(in)                   :: Nx, Ny, N
    real(PR),dimension(Nx*Ny),intent(in) :: U
    real(PR), intent(in)                 :: dx, dy
    integer                              :: i, j
    character(len=20)                    :: F_NAME
    character(len=40)                    :: s1, s2, s3, s4, s5, s



    F_NAME='fichier/T'
    write(F_NAME (10:13),'(I4.4)') N
    F_NAME(14:17)= '.vtk'


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

    do j = 1,Ny
      do i = 1,Nx

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

  !> @brief Ecrit au format vtk ...
  !!
  !!..........
  !! @param[in] U
  !! @param[in] Numx,Numy,Numz
  !! @param[in] deltax,deltay,deltaz
  !! @param[in] N
  subroutine  write3dVtk(U, Numx,Numy,Numz, deltax,deltay,deltaz, N)

    integer,intent(in)                            :: Numx,Numy,Numz, N
    real(PR)                                      :: deltax,deltay,deltaz
    real(PR),dimension(num1:numN),intent(in) :: U
    integer                             :: i, j, k, num, ideb, ifin
    character(len=20)                   :: F_NAME
    character(len=40)                   :: s1, s2, s3, s4, s5, s
    real(PR),dimension(Numx*Numy*Numz) :: X,chimie


    if(Me/=0)then
       call MPI_SEND(U(num1:numN),numN-num1+1,MPI_REAL8,0,100,MPI_COMM_WORLD,statinfo)
       !call MPI_SEND(eta(num1:numN),numN-num1+1,MPI_REAL8,0,102,MPI_COMM_WORLD,statinfo)
       chimie(num1:numN) = eta(num1:numN)*fraction_vol(num1:numN,1)/0.76
       call MPI_SEND(chimie(num1:numN),numN-num1+1,MPI_REAL8,0,102,MPI_COMM_WORLD,statinfo)
    end if
    if(Me==0)then
       X(num1:numN) = U
       chimie(num1:numN) = eta(num1:numN)*fraction_vol(num1:numN,1)/0.76
       do i=1,Np-1
          call charge(Nz,Np,i,ideb,ifin)
          call MPI_RECV(X((ideb-1)*Nx*Ny+1:ifin*Nx*Ny),&
               Nx*Ny*(ifin-ideb+1),MPI_REAL8,i,100,MPI_COMM_WORLD,status,statinfo)
          call MPI_RECV(chimie((ideb-1)*Nx*Ny+1:ifin*Nx*Ny),&
               Nx*Ny*(ifin-ideb+1),MPI_REAL8,i,102,MPI_COMM_WORLD,status,statinfo)
       end do
    end if

    if(Me==0)then
       F_NAME='fichier/T'
       write(F_NAME (10:13),'(I4.4)') N
       F_NAME(14:17)= '.vtk'

       open(unit=2, file=F_NAME, action="write")

       write(2,'(a)')       "# vtk DataFile Version 2.0"
       write(2,'(a)')       "Titre"
       write(2,'(a)')       "ASCII"
       write(2,'(a)')       "DATASET STRUCTURED_POINTS"
       write(s1,'(I4)')     Numx
       write(s2,'(I4)')     Numy
       write(s3,'(I4)')     Numz
       write(2,'(a)')       "DIMENSIONS " // adjustl(trim(s1)) // adjustl(trim(s2)) // adjustl(trim(s3))
       write(2,'(a)')       "ORIGIN 0.0 0.0 0.0"
       write(s1,'(F20.10)')  deltax !dx
       write(s2,'(F20.10)')  deltay !dx
       write(s3,'(F20.10)')  deltaz
       write(2,'(a)')       "SPACING " // adjustl(trim(s1)) // adjustl(trim(s2)) // adjustl(trim(s3))
       write(s5,'(I8)')     (Numx*Numy*Numz)
       write(2,'(a)')       "POINT_DATA " // adjustl(trim(s5))
       write(2,'(a)')       "SCALARS Temperature double"
       write(2,'(a)')       "LOOKUP_TABLE default"


       do k = 1,Numz
          do j = 1,Numy
             do i = 1,Numx
                num = (k-1)*Nx*Ny + (j-1)*Nx + i

                if(modulo(j,10) == 0) then
                   write(2,'(/)',advance='no')
                end if
                s = ""
                write(s,'(F20.10)')         X(num)
                write(2,'(a)',advance='no') adjustl(trim(s)) // " "
             enddo
          end do
       end do

       write(2,'(a)')       "SCALARS Eta double"
       write(2,'(a)')       "LOOKUP_TABLE default"

       do k = 1,Numz
          do j = 1,Numy
             do i = 1,Numx
                num = (k-1)*Nx*Ny + (j-1)*Nx + i

                if(modulo(j,10) == 0) then
                   write(2,'(/)',advance='no')
                end if
                s = ""
                write(s,'(F20.10)')         chimie(num)
                write(2,'(a)',advance='no') adjustl(trim(s)) // " "
             enddo
          end do
       end do


       close(2)

       print*, "SVG -> ", F_NAME
    end if
  end subroutine write3dVtk

  !> @brief Ecrit au format vtk ...
  !!
  !!..........
  !! @param[in] U
  !! @param[in] Nx,Ny
  !! @param[in] dx,dy
  !! @param[in] N
  subroutine  writeFibreVtk(U, porox,poroy,poroz)

    integer,intent(in)                  :: porox,poroy,poroz
    integer,dimension(porox,poroy,poroz),intent(in) :: U
    integer                             :: i, j, k
    character(len=20)                   :: F_NAME
    character(len=40)                   :: s1, s2, s3, s4, s5, s


    F_NAME='fichier/fibre.vtk'

    open(unit=2, file=F_NAME, action="write")

    write(2,'(a)')       "# vtk DataFile Version 2.0"
    write(2,'(a)')       "Titre"
    write(2,'(a)')       "ASCII"
    write(2,'(a)')       "DATASET STRUCTURED_POINTS"
    write(s1,'(I4)')     porox
    write(s2,'(I4)')     poroy
    write(s3,'(I4)')     poroz
    write(2,'(a)')       "DIMENSIONS " // adjustl(trim(s1)) // adjustl(trim(s2)) // adjustl(trim(s3))
    write(2,'(a)')       "ORIGIN 0.0 0.0 0.0"
    write(s4,'(F20.10)') 1. !dx
    write(2,'(a)')       "SPACING " // adjustl(trim(s4)) // adjustl(trim(s4)) // adjustl(trim(s4))
    write(s5,'(I8)')     (porox*poroy*poroz)
    write(2,'(a)')       "POINT_DATA " // adjustl(trim(s5))
    write(2,'(a)')       "SCALARS Fibre integer"
    write(2,'(a)')       "LOOKUP_TABLE default"

    do k = 1,poroz
      do j = 1,poroy
        do i = 1,porox

          if(modulo(j,10) == 0) then
            write(2,'(/)',advance='no')
          end if
          s = ""
          write(s,'(I4)')         U(i,j,k)
          write(2,'(a)',advance='no') adjustl(trim(s)) // " "
        enddo
      end do
    end do

    close(2)

  end subroutine writeFibreVtk

  !*******************************************************!
  !*******************************************************!

  !> @brief Ecrit au format dat...
  !!
  !!..........
  !! @param[in] U
  !! @param[in] Nx,Ny
  !! @param[in] dx,dy
  !! @param[in] N
  subroutine printvector(X, Y, N)

    integer,intent(in)                               :: N
    real(PR),dimension((k1-1)*Nx*Ny+1:kN*Nx*Ny),intent(in) :: X,Y
    integer                                          :: i, j
    character(len=20)                                :: F_NAME

    if (N<10) then
       F_NAME='fichier/T'
       write(F_NAME (10:10),'(I1)') N
       F_NAME(11:11)='P'
       write(F_NAME (12:12),'(I1)') Me
       F_NAME(13:16)= '.dat'

    elseif ((N>=10).and.(N<100)) then
       F_NAME='fichier/T'
       write(F_NAME (10:11),'(I2)') N
       F_NAME(12:12)='P'
       write(F_NAME (13:13),'(I1)') Me
       F_NAME(14:17)= '.dat'

    elseif ((N>=100).and.(N<1000)) then
       F_NAME='fichier/T'
       write(F_NAME (10:12),'(I3)') N
       F_NAME(13:13)='P'
       write(F_NAME (14:14),'(I1)') Me
       F_NAME(15:18)= '.dat'

    else if ((N>=1000).and.(N<10000)) then
       F_NAME='fichier/T'
       write(F_NAME (10:13),'(I4)') N
       F_NAME(14:14)='P'
       write(F_NAME (15:15),'(I1)') Me
       F_NAME(16:19)= '.dat'

    else if ((N>=10000).and.(N<100000)) then
       F_NAME='fichier/T'
       write(F_NAME (10:14),'(I5)') N
       F_NAME(15:15)='P'
       write(F_NAME (16:16),'(I1)') Me
       F_NAME(17:20)= '.dat'
    end if

    open(unit=2, file=F_NAME, action="write")

    do i=1,Nx
       do j=1,Ny
          write(2,*) i*dx, j*dy, X((j-1)*Nx+i), Y((j-1)*Nx+i)
       end do
       write(2,*)
    end do

    close(2)

    if(Me==0) then
       print*, 'SVG numero: ',N
    end if

  end subroutine printvector


  !*******************************************************!
  !*******************************************************!
  !*******************************************************!
  !*******************************************************!

  !> @brief Resout le systeme lineaire Ax=b
  !!
  !!..........
  !! @param[in] U
  !! @param[in] Nx,Ny
  !! @param[in] dx,dy
  !! @param[in] N
  subroutine Gradient_conjugue(X,b,eps)
    implicit none
    real(PR),dimension(num1:numN) :: x
    real(PR),dimension(num1:numN) :: b
    real(PR)                     :: eps

    integer,parameter :: NitMax = 100000

    real(PR),dimension(num1:numN) :: residu,Ap
    real(PR),dimension(num1-Nx*Ny:numN+Nx*Ny) :: p,x_r
    real(PR)                     :: norm_r0,norm_r,norm_rtotal,norm_b,norm_btot,App,App_total,err
    integer                      :: n,nb_iter
    integer,dimension(MPI_STATUS_SIZE) :: status


    x_r = 0.

    if(Me<Np-1) then ! pas optimisé
       call MPI_SEND(x(numN-Nx*Ny+1:numN),Nx*Ny,mpi_real8,Me+1,101,MPI_COMM_WORLD,statinfo)
    end if
    if(Me>0) then
       call MPI_RECV(x_r(num1-Nx*Ny:num1-1),Nx*Ny,mpi_real8,Me-1,101,MPI_COMM_WORLD,status,statinfo)
       call MPI_SEND(x(num1:num1+Nx*Ny-1),Nx*Ny,mpi_real8,Me-1,101,MPI_COMM_WORLD,statinfo)
    end if
    if(Me<Np-1) then
       call MPI_RECV(x_r(numN+1:numN+Nx*Ny),Nx*Ny,mpi_real8,me+1,101,MPI_COMM_WORLD,status,statinfo)
    end if

    x_r(num1:numN) = x
    residu  = b - Matmula(x_r)
    p(num1:numN) = residu

    !Calcul de la norme et l envoie a tout le monde
    norm_r = dot_product(residu,residu)
    call MPI_ALLREDUCE(norm_r,norm_rtotal,1,mpi_real8,MPI_SUM,MPI_COMM_WORLD,statinfo)
    norm_r0 = norm_rtotal
    !norm_b  = norm_rtotal
    norm_b = dot_product(b,b)
    call MPI_ALLREDUCE(norm_b,norm_btot,1,mpi_real8,MPI_SUM,MPI_COMM_WORLD,statinfo)
    err = norm_btot*eps**2

    ! Boucle
    nb_iter=1
    do while ((norm_rtotal >=  err).and.(nb_iter<NitMax))

       !Reconstruction de p(i1-Nx:in+Nx)
       if(Me<Np-1) then
          call MPI_SEND(p(numN-Nx*Ny+1:numN),Nx*Ny,mpi_real8,Me+1,101,MPI_COMM_WORLD,statinfo)
       end if
       if(Me>0) then
          call MPI_RECV(p(num1-Nx*Ny:num1-1),Nx*Ny,mpi_real8,Me-1,101,MPI_COMM_WORLD,status,statinfo)
          call MPI_SEND(p(num1:num1+Nx*Ny-1),Nx*Ny,mpi_real8,Me-1,101,MPI_COMM_WORLD,statinfo)
       end if
       if(Me<Np-1) then
          call MPI_RECV(p(numN+1:numN+Nx*Ny),Nx*Ny,mpi_real8,me+1,101,MPI_COMM_WORLD,status,statinfo)
       end if

       Ap = Matmula(p)

       ! Calcul de App_total
       App = dot_product(Ap,p(num1:numN))
       call MPI_ALLREDUCE(App,App_total,1,mpi_real8,MPI_SUM,MPI_COMM_WORLD,statinfo)

       x = x + (norm_rtotal/App_total)*p(num1:numN)

       residu = residu-(norm_rtotal/App_total)*Ap

       ! Calcul de la norme du residu
       norm_r = dot_product(residu,residu)!Norme par partie
       call MPI_ALLREDUCE(norm_r,norm_rtotal,1,mpi_real8,MPI_SUM,MPI_COMM_WORLD,statinfo)

       p(num1:numN) = residu + (norm_rtotal/norm_r0)*p(num1:numN)

       norm_r0 = norm_rtotal
       nb_iter = nb_iter + 1
       iteration = iteration + 1
    end do

  end subroutine Gradient_conjugue

  !> @brief Effectue le produit matrice vecteur Au
  !!
  !!..........
  !! @param[in] U
  function Matmula(U)
    implicit none
    real(PR),dimension(num1-Nx*Ny:numN+Nx*Ny) :: U
    real(PR),dimension(num1:numN)             :: Matmula
    integer:: i,j,k, num


    Matmula=0
    do k= k1,kN
      do j= 1,Ny
        do i=1,Nx
          num=(k-1)*Nx*Ny + (j-1)*Nx + i
          ! print*, "---",Nx,Ny,Nz
          ! print*, i,j,k,num

          ! Cz gauche
          if(k/=1)then
            Matmula(num) = Matmula(num) + Cz(num-Nx*Ny)*U(num-Nx*Ny)  !!
          end if
          ! Cy gauche
          if(j/=1)then
            Matmula(num) = Matmula(num) + Cy(num-k*Nx)*U(num-Nx)  !!
          end if
          ! Cx gauche
          if(i/=1)then
            Matmula(num) = Matmula(num) + Cx(num-j-(k-1)*Ny)*U(num-1)    !!
          end if
          ! Cd
          Matmula(num)   = Matmula(num) + Cd(num)*U(num)
          ! Cx droite
          if(i/=Nx)then
            Matmula(num)=Matmula(num)+Cx(num-(j-1)-(k-1)*Ny)*U(num+1)    !!
          end if
          ! Cy droite
          if(j/=Ny)then
            Matmula(num) = Matmula(num) + Cy(num-(k-1)*Nx)*U(num+Nx)  !!
          end if
          ! Cz droite
          if(k/=Nz)then
            Matmula(num) = Matmula(num) + Cz(num)*U(num+Nx*Ny)  !!
          end if
        end do
      end do
    end do
  end function Matmula

  !> @brief Interpole la valeur (col 2) de tab en fonction de la la valeur voulue (col 1)
  !!
  !!..........
  function interp(Tab,val)
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




  ! **** POUR CALCULER UN TENSEUR LOCAL DE CONDUCTIVITE QUI SUIT LES FIBRES ****
  ! D = P. Delta . inv(P)
  ! où Delta est la matrice qui suit la direction des fibres
  ! P est construite sur la direction des fibres "dirfib" plus deux autres vecteurs qui lui sont orthogonaux.
  ! Delta = diag (Dperp, Dperp, Dparall)
  !

  !> @brief Recupere la direction des fibres pour le point (i,j,k)
  !!
  !!..........
  subroutine getdirfib(pos,dirfib)
    integer,dimension(3),intent(in) :: pos ! (i,j,k)
    integer,dimension(3),intent(out) :: dirfib

    dirfib(1) = Orien(pos(1),pos(2),pos(3),1)
    dirfib(2) = Orien(pos(1),pos(2),pos(3),2)
    dirfib(3) = Orien(pos(1),pos(2),pos(3),3)
  end subroutine getdirfib

  !> @brief Tenseur local de conductivite selon (x,y,z)
  subroutine difglo(difpa, difpe, dirfib, MatD)
    real(PR)                           :: difpa, difpe
    real(PR), dimension(3), intent(inout) :: dirfib
    real(PR), dimension(3,3), intent(out) :: MatD

    real(PR), dimension(3,3) :: MatDelta, MatInv, MatP, MatM
    real(PR), dimension(3)   :: vecA, vecB

    if (difpa == difpe) then ! Si lambda isotrope
      MatD = 0.
      MatD(1,1) = difpe
      MatD(2,2) = difpe
      MatD(3,3) = difpa
    else
      !  /* 1. faire la matrice de passage */
      call vec_norm(dirfib)
      call cree_mat_pass(dirfib,vecA,vecB) ! /* trois vecteurs V1,V2,V3 */
      call coltomat(vecA,vecB,dirfib,MatM)

      ! /* 2. obtenir l'inverse de cette matrice */
      call inv_mat_pass(vecA,vecB,dirfib,MatInv)

      ! /* 3. fabriquer la matrice diagonale */
      MatDelta = 0.
      MatDelta(1,1) = difpe
      MatDelta(2,2) = difpe
      MatDelta(3,3) = difpa

      ! /* 4. calculer P = Delta.INV */
      call prodmat(MatDelta,MatInv,MatP);

      ! /* 5. calculer D = M.Delta.inv(M) */
      call prodmat(MatM,MatP,MatD);

      ! if (D->X.x == 0.) then print*, "Aie ! X, difpa, difpe :",difpa,difpe
    end if
  end subroutine difglo


  !> @brief Creation d'une matrice a partir de trois vecteurs colonnes
  subroutine coltomat(vec1,vec2,vec3,Mat)
    real(PR),dimension(3,3), intent(out) :: Mat
    real(PR),dimension(3), intent(in) :: vec1, vec2, vec3

    Mat(1,1) = vec1(1)
    Mat(2,1) = vec1(2)
    Mat(3,1) = vec1(3)
    Mat(1,2) = vec2(1)
    Mat(2,2) = vec2(2)
    Mat(3,2) = vec2(3)
    Mat(1,3) = vec3(1)
    Mat(2,3) = vec3(2)
    Mat(3,3) = vec3(3)
  end subroutine coltomat


  !> @brief cree une bond (base) A, B, V a partir d'une direction V donnee normee
  subroutine cree_mat_pass(dir,per1,per2)
    real(PR),dimension(3),intent(in)  :: dir
    real(PR),dimension(3),intent(out) :: per1, per2
    real(PR),dimension(3)             :: vec

    if ( dir(3)==1. ) then ! dir=(0,0,1)
      per1(1) = 1.; per1(2) = 0.; per1(3) = 0.;
      per2(1) = 0.; per2(2) = 1.; per2(3) = 0.;
    elseif ( dir(3)==-1. ) then ! dir=(0,0,-1)
      per1(1) = 1.; per1(2) = 0. ; per1(3) = 0.;
      per2(1) = 0.; per2(2) = -1.; per2(3) = 0.;
    else
      vec(1) = 0.; vec(2) = 0.; vec(3) = 1.;
      call vec_crotz(dir,vec,per1);  ! A = V x Z
      call vec_norm(per1);
      call vec_crotz(dir,per1,per2);  ! B = V x A
      call vec_norm(per2);
    end if

  end subroutine cree_mat_pass


  !> @brief inverse d'une matrice de determinant unite formee par 3 vecteurs unites
  subroutine inv_mat_pass(vec1,vec2,vec3,inv)
    real(PR),dimension(3,3), intent(out) :: inv
    real(PR),dimension(3), intent(in) :: vec1, vec2, vec3

    inv(1,1) = vec2(2)*vec3(3) - vec2(3)*vec3(2)
    inv(1,2) = vec2(3)*vec3(1) - vec2(1)*vec3(3)
    inv(1,3) = vec2(1)*vec3(2) - vec2(2)*vec3(1)

    inv(2,1) = vec1(3)*vec3(2) - vec1(2)*vec3(3)
    inv(2,2) = vec1(1)*vec3(3) - vec1(3)*vec3(1)
    inv(2,3) = vec1(2)*vec3(1) - vec1(1)*vec3(2)

    inv(3,1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
    inv(3,2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
    inv(3,3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
  end subroutine inv_mat_pass


  !> @brief Produit vectoriel AxB
  subroutine vec_crotz(vec1,vec2, res)
    real(PR),dimension(3) :: res
    real(PR),dimension(3) :: vec1, vec2

    res(1) =  vec1(2)*vec2(3) - vec2(2)*vec1(3)
    res(2) =  vec1(3)*vec2(1) - vec2(3)*vec1(1)
    res(3) =  vec1(1)*vec2(2) - vec2(1)*vec1(2)
  end subroutine vec_crotz


  !> @brief Retourne le vecteur norme
  subroutine vec_norm(vec)
    real(PR), dimension(3), intent(inout) :: vec
    real(PR)                              :: norme

    norme = vec(1)**2 + vec(2)**2 + vec(3)**2
    norme = sqrt(norme)

    vec(1) = vec(1) / norme
    vec(2) = vec(2) / norme
    vec(3) = vec(3) / norme
  end subroutine vec_norm


  !> @brief Produit de deux matrices 3*3
  subroutine prodmat(A,B,res)
    real(PR), dimension(3,3), intent(in) :: A,B
    real(PR), dimension(3,3), intent(out) :: res

    res=matmul(A,B)
  end subroutine prodmat

  subroutine charge(n,Np,me,k1,kN)
    integer,intent(in)::n,Np,me
    integer,intent(out)::k1,kN
    integer::q
    q=n/Np
    if(me<mod(n,Np))then
       k1=me*(q+1)+1
       kN=(me+1)*(1+q)
    else
       k1=mod(n,Np)+1+me*q
       kN=k1+q-1
    end if
  end subroutine charge



end module mod_fonction
