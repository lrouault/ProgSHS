program main
  use donnees
  use mod_preprocess
  use mod_physique
  use mod_fonction

  implicit none

  integer  :: i,j,k,iter

  call initialisation("param.dat")
  call creation_matrice()

  print *, "rho cp lambda dt Tad",rho,cp,lambda,dt,Tad

  call printvector(U0,0)

  flux=12000000.*2

  write(*,*) dx*(Nx+1),dx,Nx
  !*************Marche en temps*********************
  iter=0
  Niter = nint(tmax/dt)
  !do k = 1,Niter
  do while(time<tmax)  
     time = time + dt
     iter=iter+1

     call second_membre()
     call Gradient_conjugue(U,rhs,epsilon) !V+chi

     if(eta((Ny/2)*Nx+2)==1.)then
        flux = 0.
     end if

     if(modulo(iter*nb_fichiers,Niter) == 0) then
        call printvector(U, iter*nb_fichiers/Niter)
     end if

     U0 = U

  end do
  
  call fin()

end program main
