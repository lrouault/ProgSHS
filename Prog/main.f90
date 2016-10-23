program main
  use donnees
  use mod_preprocess
  use mod_physique
  use mod_fonction

  implicit none

  integer::i,j,k

  call initialisation("param.dat")
  call creation_matrice()

  print *, "rho cp lambda dt Tad",rho,cp,lambda,dt,Tad
  print*,'D',D

  call printvector(U0,0)

  h=h*dt/(Ly*rho*cp)
  ! h=1.0e-2
  print*,'h',h


  !*************Marche en temps*********************
  do k = 1,Niter
    !..! print*, "BOUCLE : ",k !!!!!!!!!!!!!!!
     t = k*dt
     do j=1,Ny
        U0(bij(1,j,Ny)) = min(Tad, Tad*(5*t/Tmax)+Text)
     end do
    !..! print*,"Calc cd U0 : OK" !!!!!!!!!!!!
     ! Conditions limites :
     !******************définition de V****************

     V = 0
     do j = 1,Ny
        V(bij(1,j,Ny)) = -c*Tad + V(bij(1,j,Ny))!-c*min(Tad,Tad*(500*t/Tmax))+V(bij(1,j,Ny))
        V(bij(Nx,j,Ny)) = -c*Text + V(bij(Nx,j,Ny))!-h*dt/(rho*cp*Lx)*(U0(bij(Nx,j,Ny))-Text)
     end do

     do i = 1,Nx
        V(bij(i,1,Ny)) = +V(bij(i,1,Ny)) - h*(U(bij(i,1,Ny))-Text)!+V(bij(i,1,Ny))!-b*Text+V(bij(i,1,Ny))
        V(bij(i,Ny,Ny)) = +V(bij(i,Ny,Ny)) - h*(U(bij(i,Ny,Ny))-Text)!+V(bij(i,Ny,Ny))-b*Text+V(bij(i,Ny,Ny))
     end do
     !..!print*,"Calc cd lim : OK" !!!!!!!!!!!!
     !    Chimie  :
     !*******************définition de eta/chi**************

     eta0 = eta
     tmp = 0.
     do i = 1,size(eta)
        coef   = -Ea/(R*U0(i))
        !..!print*,"coeff, k0",coef,k0,dt !!!!!!!!!!!!!!!!!!!!!!
        eta(i) = (eta(i)+dt*k0*exp(coef)) / (1+k0*dt*exp(coef))
        !..!print*,"eta",eta(i)
        tmp    = max(tmp, eta(i)-eta0(i))
        !..!print*, "tmp", tmp
        chi(i) = Q/Cp * (eta(i)-eta0(i))!k0*(1-eta(i))*exp(coef)/(1,77e+1)
     end do
     !..!print*,"Calc eta et chi : OK" !!!!!!!!!!!!
     !print*,'tmp=',tmp

     !**************************************************

     call gradc(a, b, c, epsilon, Nmax, U, U0+V+chi)
     !..!print*,"Calc gradc : OK"
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

end program main
