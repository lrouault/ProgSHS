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

    integer,intent(in)                   :: Nx, Ny, N
    real(PR),dimension(Nx*Ny),intent(in) :: U
    real(PR), intent(in)                 :: dx, dy
    integer                              :: i, j
    character(len=20)                    :: F_NAME
    character(len=40)                    :: s1, s2, s3, s4, s5, s



    F_NAME='fichier/T'
    write(F_NAME (10:14),'(I4.4)') N
    F_NAME(15:18)= '.vtk'


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




  ! **** POUR CALCULER UN TENSEUR LOCAL DE CONDUCTIVITE QUI SUIT LES FIBRES ****
  ! D = P. Delta . inv(P)
  ! où Delta est la matrice qui suit la direction des fibres
  ! P est construite sur la direction des fibres "dirfib" plus deux autres vecteurs qui lui sont orthogonaux.
  ! Delta = diag (Dperp, Dperp, Dparall)
  !

  ! Recupere la direction des fibres pour le point (i,j,k)
  subroutine getdirfib(pos,dirfib)
    integer,dimension(3),intent(in) :: pos ! (i,j,k)
    integer,dimension(3),intent(out) :: dirfib

    dirfib(1) = Orien(pos(1),pos(2),pos(3),1)
    dirfib(2) = Orien(pos(1),pos(2),pos(3),2)
    dirfib(3) = Orien(pos(1),pos(2),pos(3),3)
  end subroutine getdirfib

  ! Tenseur local de conductivite selon (x,y,z)
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


  ! /* Creation d'une matrice a partir de trois vecteurs colonnes */
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


  ! /* cree une bond (base) A, B, V a partir d'une direction V donnee normee */
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


  ! /* inverse d'une matrice de determinant unite formee par 3 vecteurs unites */
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


  ! /* Produit vectoriel AxB */
  subroutine vec_crotz(vec1,vec2, res)
    real(PR),dimension(3) :: res
    real(PR),dimension(3) :: vec1, vec2

    res(1) =  vec1(2)*vec2(3) - vec2(2)*vec1(3)
    res(2) =  vec1(3)*vec2(1) - vec2(3)*vec1(1)
    res(3) =  vec1(1)*vec2(2) - vec2(1)*vec1(2)
  end subroutine vec_crotz


  ! /* Retourne le vecteur norme */
  subroutine vec_norm(vec)
    real(PR), dimension(3), intent(inout) :: vec
    real(PR)                              :: norme

    norme = vec(1)**2 + vec(2)**2 + vec(3)**2
    norme = sqrt(norme)

    vec(1) = vec(1) / norme
    vec(2) = vec(2) / norme
    vec(3) = vec(3) / norme
  end subroutine vec_norm


  ! /* Produit de deux matrices 3*3 */
  subroutine prodmat(A,B,res)
    real(PR), dimension(3,3), intent(in) :: A,B
    real(PR), dimension(3,3), intent(out) :: res

    res=matmul(A,B)
  end subroutine prodmat



end module mod_fonction
