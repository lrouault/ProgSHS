!> @mainpage Accueil Doxygen
!!
!! @section Introduction
!! Ce document décrit les sources du code resolvant le probleme ..
!!
!! @section Information
!! Pour plus d'informations regarder les pages documentations du site.


program main
  use donnees
  use mod_preprocess
  use mod_physique
  use mod_fonction
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR,C_INT32_T,C_FLOAT,C_SIGNED_CHAR

  implicit none
  integer  :: iter



  ! character(kind=C_SIGNED_CHAR)      :: data2
  ! character(kind=C_CHAR)      :: data3
  ! integer :: i,j,k
  !
  ! open( 102 , file="image/tex1_249x249x100_vol_uchar.raw", form="unformatted", action="read", &
  ! access="stream")
  ! 1	FORMAT ( Z4 ) ! Format pour ecrire les caractere en hexadecimal
  ! do k=1,249*249*100
  !   ! do j=1,poroy
  !   !   do i=1,porox
  !       read(102) data3
  !       if (ichar(data3)/=0)  print*,ichar(data3)
  !   !   enddo
  !   ! enddo
  ! enddo
  ! close(102)
  !
  !


  !
  ! call fillPoro("image/IMAGE_crop.mat","image/IMAGE_crop2.mat") ! Remplit fibre(porox/y/z)

  call fillPoro2("image/tex1_249x249x100_vol_uchar.raw",249,249,100) ! Remplit fibre(porox/y/z)
  call writeFibreVtk(Poro,porox,poroy,poroz)
  print*, porox,poroy,poroz

  !
  ! call fillOrientation("IMAGE_crop.or")
  !
  call initialisation("param.dat")
  call write3dVtk(U,Nx,Ny,Nz,  dx,dy,dz, 0)

  call creation_matrice()

  !
  ! print *, "rho cp lambda dt Tad",rho,cp,lambda,dt,Tad
  !
  ! call printvector(U0,0)



  ! print*, interp(cp_Si,550._PR)
  ! print*, interp(cp_N2,550._PR)
  ! print*, interp(cp_Si3N4,550._PR)
  ! print*, interp(cp_fibre,550._PR)
  ! print*, interp(lambda_Si,550._PR)
  ! print*, interp(lambda_N2,550._PR)
  ! print*, interp(lambda_Si3N4,550._PR)
  ! print*, interp(lambda_fibre,550._PR)



  write(*,*) dx*(Nx+1),dx,Nx,dt
  !*************Marche en temps*********************
  iter=0
  Niter = nint(tmax/dt)
  iteration = 0
  do while(time<tmax)
     time = time + dt
     iter=iter+1

     call creation_matrice() ! Construit rhocp, rho, lambda, Cd Cx Cy

     rhs = rhocp*U + rho*Q*eq_arrhenius()*fraction_vol(:,1) + cd_lim() + chauffage()
    ! print*,"iter: ",iter," U1:",U(1)
     call Gradient_conjugue(U,rhs,epsilon) !V+chi

     if(modulo(iter*nb_fichiers,Niter) == 0) then
        ! call printvector(U, iter*nb_fichiers/Niter)
        ! call writeVtk(U, Nx, Ny, dx, dy, iter*nb_fichiers/Niter)
        call write3dVtk(U, Nx,Ny,Nz, dx,dy,dz, iter*nb_fichiers/Niter)
        print*, "temps : ",time," T0.5 : ",U((Nz/2-1)*Nx*Ny + Nx/2)
     end if

  end do

  call fin()

  write(*,*) "iterations : ", iteration

end program main
