program main
  use donnees
  use mod_preprocess
  use mod_physique
  use mod_fonction

  implicit none

  integer  :: iter

  ! call fillPoro("IMAGE_crop.mat","IMAGE_crop2.mat") ! Remplit fibre(porox/y/z)
  ! call writeFibreVtk(Poro,porox,poroy,poroz)
  !
  ! call fillOrientation("IMAGE_crop.or")
  !
  call initialisation("param.dat")
  call write3dVtk(U,Nx,Ny,Nz,  dx,dy,dz, 0)

  call creation_matrice()


  !print *, "rho cp lambda dt Tad",rho,cp,lambda,dt,Tad

  ! call printvector(U0,0)

  write(*,*) dx*(Nx+1),dx,Nx
  !*************Marche en temps*********************
  iter=0
  Niter = nint(tmax/dt)
  iteration = 0
  do while(time<tmax)
     time = time + dt
     iter=iter+1

     call creation_matrice() ! Construit rhocp, rho, lambda, Cd Cx Cy

     rhs = rhocp*U + rho*Q*eq_arrhenius() + cd_lim() + chauffage()

     call Gradient_conjugue(U,rhs,epsilon) !V+chi

     if(modulo(iter*nb_fichiers,Niter) == 0) then
        ! call printvector(U, iter*nb_fichiers/Niter)
        ! call writeVtk(U, Nx, Ny, dx, dy, iter*nb_fichiers/Niter)
        call write3dVtk(U, Nx,Ny,Nz, dx,dy,dz, iter*nb_fichiers/Niter)
     end if

  end do

  call fin()

  write(*,*) "iterations : ", iteration

end program main
