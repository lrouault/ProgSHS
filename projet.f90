program projet
  use sequentiel
  implicit none

integer::i,j,k

  call initialisation("param.dat")

  call creation_matrice()
 
 print *, "rho cp lambda dt Tad",rho,cp,lambda,dt,Tad

  print*,'D',D

  call printvector(U0,Nx,Ny,dx,dy,0)

  !Conditions initiales ???

  h=h*dt/(Ly*rho*cp)
  ! h=1.0e-2
  print*,'h',h


  !*************Marche en temps*************************************************************************
  do k = 1,Niter
     t = k*dt


     do j=1,Ny
	U0(bij(1,j,Ny))=min(Tad,Tad*(5*t/Tmax)+Text)
     end do


     ! Conditions limites :
     !******************définition de V****************
     V=0
     do j = 1,Ny		
        V(bij(1,j,Ny))=-c*Tad+V(bij(1,j,Ny))!-c*min(Tad,Tad*(500*t/Tmax))+V(bij(1,j,Ny))
        V(bij(Nx,j,Ny))=-c*Text+V(bij(Nx,j,Ny))!-h*dt/(rho*cp*Lx)*(U0(bij(Nx,j,Ny))-Text)
     end do
     do i = 1,Nx
        V(bij(i,1,Ny))=+V(bij(i,1,Ny))-h*(U(bij(i,1,Ny))-Text)!+V(bij(i,1,Ny))!-b*Text+V(bij(i,1,Ny))
        V(bij(i,Ny,Ny))=+V(bij(i,Ny,Ny))-h*(U(bij(i,Ny,Ny))-Text)!+V(bij(i,Ny,Ny))-b*Text+V(bij(i,Ny,Ny))
     end do
     !**************************************************

     !~ print*,'t=',t

     !~ print*,'V='
     !~ do i=1,size(V)                     
     !~ print*,V(i)
     !~ if(modulo(i,Ny)==0) then
     !~ print*,' '
     !~ end if
     !~ end do
     !~ print*,'-------------------------------------------------------------'



     !    Chimie  :
     !*******************définition de eta/chi**************
     eta0=eta
     tmp=0.
     do i=1,size(eta)
        coef=-Ea/(R*U0(i))

        eta(i)=(eta(i)+dt*k0*exp(coef))/(1+k0*dt*exp(coef))
        tmp=max(tmp,eta(i)-eta0(i))
        chi(i)=Q/Cp*(eta(i)-eta0(i))!k0*(1-eta(i))*exp(coef)/(1,77e+1)
     end do

     !print*,'tmp=',tmp

     !**************************************************

     call gradc(a,b,c,epsilon,Nmax,Nx,Ny,U,U0+V+chi,h)


     !~ do i=1,size(U)
     !~ 	if(U(i)>=Tad)then      
     !~ 		U(i)=U(i)
     !~ 	end if
     !~ end do


     if(modulo(k,100)==0) then
        call printvector(chi,Nx,Ny,dx,dy,k/100)
     end if
     U0 = U
  end do

  !******************************************************************************************************
  
end program projet
