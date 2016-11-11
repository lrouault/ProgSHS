program main
  use donnees
  use mod_preprocess
  use mod_physique
  use mod_fonction
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR

  implicit none
  integer  :: iter


  call fillPoro("IMAGE_crop.mat","IMAGE_crop2.mat") ! Remplit fibre(porox/y/z)
  call writeFibreVtk(Poro,porox,poroy,poroz)


  call initialisation("param.dat")

  call creation_matrice()
  call printvector(U0,0)


  flux=12000000.*2

  write(*,*) dx*(Nx+1),dx,Nx
  !*************Marche en temps*********************
  iter=0
  Niter = nint(tmax/dt)

  do while(time<tmax)
     time = time + dt
     iter=iter+1

     call creation_matrice() ! Construit rhocp, rho, lambda, Cd Cx Cy

     rhs = rhocp*U + rho*Q*eq_arrhenius() + cd_lim()

     call Gradient_conjugue(U,rhs,epsilon) !V+chi

     if(eta((Ny/2)*Nx+2)==1.)then
        flux = 0.
     end if

     if(modulo(iter*nb_fichiers,Niter) == 0) then
        call printvector(U, iter*nb_fichiers/Niter)
        call writeVtk(U, Nx, Ny, dx, dy, iter*nb_fichiers/Niter)
     end if

  end do

  call fin()
!  deallocate(Poro)

end program main
