program projet

  use sequentiel
  implicit none
  integer::N,i,j,Nx,Ny,Niter,Nmax,Ndisplay,statinfo,me,np,i1,im,k
  real*8::Lx, Ly, D, dt, Tmax, t, dx, dy, a, b, c, epsilon, x, y,Text,Tad,rho,cp,lambda,h,Ea,R,k0,Q,coef,tmp
  real*8, dimension(:), allocatable::U,U0,V,Fm,eta,chi,eta0


  call read_param("param.dat", Nx, Ny, Lx, Ly,rho,cp,lambda)


  N = Nx*Ny
  Tmax = 1000
  Niter = 10000
  Nmax = n+1
  epsilon = 1e-7
  Text=298
  Tad=2300.0
  h=10.0d+0
  Ea=266547
  R=8.3144621d+0
  k0=2.0e+4
  Q=626e+3

  call compute_vars(Nx, Ny, Lx, Ly,rho,cp,lambda,Tmax,Niter, dt, dx, dy, a, b, c,D)
  print *, rho,cp,lambda,dt,Tad



  print*,'D',D

  allocate(U(Nx*Ny),V(Nx*Ny),U0(Nx*Ny),eta(Nx*Ny),chi(Nx*Ny))

  !Conditions initiales
  !******************définition de U0***************
  U0=Text
  do i=1,Ny
     U0(bij(1,i,Ny))=(3.0/4)*Tad

  end do
  call printvector(U0,Nx,Ny,dx,dy,0)
  !*************************************************	

  V=0
  eta=0
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
  deallocate(U,U0,V)
end program projet
