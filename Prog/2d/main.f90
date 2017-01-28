program main
  use donnees
  use mod_preprocess
  use mod_physique
  use mod_fonction

  implicit none

  integer  :: iter
  real(PR) :: t1,t2,diff
  t1 = MPI_WTIME()
  call initialisation("param.dat")
!!$  if(me==0)then
!!$     call fillPoro("IMAGE_crop.mat","IMAGE_crop2.mat") ! Remplit fibre(porox/y/z)
!!$     call writeFibreVtk(Poro,porox,poroy,poroz)
!!$
!!$     call fillOrientation("IMAGE_crop.or")
!!$  end if

  
  call creation_matrice()

  !print *, "rho cp lambda dt Tad",rho,cp,lambda,dt,Tad

  call printvector(U,eta,0)

  !write(*,*) dx*(Nx+1),dx,Nx
  !*************Marche en temps*********************
  iter=0
  Niter = nint(tmax/dt)
  iteration = 0
  do while(time<tmax)
     time = time + dt
     iter = iter + 1

     call creation_matrice() ! Construit rhocp, rho, lambda, Cd Cx Cy

     rhs = rhocp*U + rho*Q*eq_arrhenius() + cd_lim() + chauffage()

     call Gradient_conjugue(U,rhs,epsilon) !V+chi

     if(modulo(iter*nb_fichiers,Niter) == 0) then
        !call printvector(U, eta, iter*nb_fichiers/Niter)
        !call writeVtk(U, Nx, Ny, dx, dy, iter*nb_fichiers/Niter)
     end if

  end do

!!$  if(me==0)then
!!$     print*, iteration, "iterations"
!!$  end if
  t2 = MPI_WTIME()

  call MPI_ALLREDUCE(t2-t1,diff,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,statinfo)
  if(me==0)write(*,*) diff
  call fin()

end program main
