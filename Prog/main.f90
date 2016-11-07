program main
  use donnees
  use mod_preprocess
  use mod_physique
  use mod_fonction

  implicit none

  integer  :: i,j,k


  call initialisation("param.dat")
  call creation_matrice()

  print *, "rho cp lambda dt Tad",rho,cp,lambda,dt,Tad

  call printvector(U0,0)

  !h=h*dt/(Ly*rho*cp) ! h=1.0e-2
  !print*,'h',h

  write(*,*) dx*(Nx+1),dx,Nx
  !*************Marche en temps*********************
  do k = 1,Niter
     t = k*dt

     !do j=1,Ny
     !   U0(bij(1,j,Ny)) = min(Tad, Text+Tad*t)
     !end do


     call cd_limite()
     call eq_arrhenius()
     call Gradient_conjugue(U, rho*cp*U0+V+chi,epsilon) !V+chi

     !~ do i=1,size(U)
     !~ 	if(U(i)>=Tad)then      
     !~ 		U(i)=U(i)
     !~ 	end if
     !~ end do

     if(modulo(k,100) == 0) then
        call printvector(U, k/100)
     end if

     U0 = U

  end do

  call fin()

end program main
