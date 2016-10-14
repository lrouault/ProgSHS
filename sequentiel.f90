module sequentiel

contains


	subroutine prod_mat_vec(a,b,c,x,y,Nx,Ny,h)
		implicit none
		real*8, intent(in)::a,b,c,h
		real*8, intent(in), dimension(Nx*Ny)::x
		real*8, intent(inout), dimension(Nx*Ny)::y
		integer,intent(in)::Nx,Ny
		integer::i,j,k

		
		
		
		
		
		! 1er bloc
		y(bij(1,1,Ny)) = (a-b)*x(bij(1,1,Ny))+b*x(bij(1,2,Ny))+c*x(bij(2,1,Ny))
		do j = 2, Ny-1
			y(bij(1,j,Ny)) = b*x(bij(1,j-1,Ny))+a*x(bij(1,j,Ny))+b*x(bij(1,j+1,Ny))+c*x(bij(2,j,Ny))
		end do
		y(bij(1,Ny,Ny)) = b*x(bij(1,Ny-1,Ny))+(a-b)*x(bij(1,Ny,Ny))+c*x(bij(2,Ny,Ny))
		! Milieu Matrice
		do i = 2, Nx-1
			y(bij(i,1,Ny)) = c*x(bij(i-1,1,Ny))+(a-b)*x(bij(i,1,Ny))+b*x(bij(i,2,Ny))+c*x(bij(i+1,1,Ny))
			do j = 2, Ny-1
				y(bij(i,j,Ny)) = c*x(bij(i-1,j,Ny))+b*x(bij(i,j-1,Ny))+a*x(bij(i,j,Ny))+b*x(bij(i,j+1,Ny))+c*x(bij(i+1,j,Ny))
			end do
			y(bij(i,Ny,Ny)) = c*x(bij(i-1,Ny,Ny))+b*x(bij(i,Ny-1,Ny))+(a-b)*x(bij(i,Ny,Ny))+c*x(bij(i+1,Ny,Ny))
		end do
		! Dernier bloc
		y(bij(Nx,1,Ny)) = c*x(bij(Nx-1,1,Ny))+(a-b)*x(bij(Nx,1,Ny))+b*x(bij(Nx,2,Ny))
		do j = 2, Ny-1
			y(bij(Nx,j,Ny)) = c*x(bij(Nx-1,j,Ny))+b*x(bij(Nx,j-1,Ny))+a*x(bij(Nx,j,Ny))+b*x(bij(Nx,j+1,Ny))
		end do
		y(bij(Nx,Ny,Ny)) = c*x(bij(Nx-1,Ny,Ny))+b*x(bij(Nx,Ny-1,Ny))+(a-b)*x(bij(Nx,Ny,Ny))
	end subroutine

	function norme(R) ! norme euclidienne
		implicit none
		real*8, intent(in), dimension(:)::R
		real*8 norme
		integer::i
		real*8::S
		S = 0
		do i = 1, size(R)
			S = S + R(i)*R(i)
		end do
		norme = sqrt(S)
	end function

	subroutine printmat(A,n,ch)
		implicit none
		real,dimension(n,n),intent(in)::A
		character(len=*),intent(in)::ch
		integer,intent(in)::n
		integer::i
		print*, ch
		do i=1,n
		   print*, A(i,:)
		end do
	end subroutine printmat
	
	function bij(i,j,Ny)
		implicit none
		integer,intent(in)::i,j,Ny
		integer::bij
		bij=Ny*(i-1)+j
	end function bij
	
	subroutine  writeVtk(U,Nx,Ny,dx,dy,N)
		implicit none
		real*8,dimension(Nx*Ny),intent(in)::U
		integer,intent(in)::Nx,Ny,N
		real*8, intent(in)::dx,dy
		integer::i,j
		character(len=20)::F_NAME
		character(len=40)::s1,s2,s3,s4,s5,s
		
		if (N<10) then
			F_NAME='fichier/T'
			write(F_NAME (10:10),'(I1)') N
			F_NAME(11:14)= '.vtk'

		elseif ((N>=10).and.(N<100)) then
			F_NAME='fichier/T'
			write(F_NAME (10:11),'(I2)') N
			F_NAME(12:15)= '.vtk'

		elseif ((N>=100).and.(N<1000)) then
			F_NAME='fichier/T'
			write(F_NAME (10:12),'(I3)') N
			F_NAME(13:16)= '.vtk'

		else if ((N>=1000).and.(N<10000)) then
			F_NAME='fichier/T'
			write(F_NAME (10:13),'(I4)') N
			F_NAME(14:17)= '.vtk'

		else if ((N>=10000).and.(N<100000)) then
			F_NAME='fichier/T'
			write(F_NAME (10:14),'(I5)') N
			F_NAME(15:18)= '.vtk'
		end if

		open(unit=2, file=F_NAME, action="write")
		
		write(2,'(a)') "# vtk DataFile Version 2.0"
		write(2,'(a)') "Titre"
		write(2,'(a)') "ASCII"
		write(2,'(a)') "DATASET STRUCTURED_POINTS"
		write(s1,'(I4)') Nx
		write(s2,'(I4)') Ny
		write(2,'(a)') "DIMENSIONS " // adjustl(trim(s1)) // adjustl(trim(s2)) // "1"
		write(2,'(a)') "ORIGIN 0.0 0.0 0.0"
		write(s3,'(F20.10)') dx
		write(s4,'(F20.10)') dy
		write(2,'(a)') "SPACING " // adjustl(trim(s3)) // adjustl(trim(s4)) // "0.0"
		write(s5,'(I8)') (Nx*Ny)
		write(2,'(a)') "POINT_DATA " // adjustl(trim(s5))
		write(2,'(a)') "SCALARS Temperature double"
		write(2,'(a)') "LOOKUP_TABLE default"
		
		do i=1,Ny
			do j=1,Nx
				if(modulo(j,10)==0) then
					write(2,'(/)',advance='no')
				end if
				s = ""
				write(s,'(F20.10)') U(bij(i,j,Ny))
				write(2,'(a)',advance='no') adjustl(trim(s)) // " "
			end do
		end do
		close(2)
	end subroutine

	subroutine printvector(U,Nx,Ny,dx,dy,N)
		real*8,dimension(Nx*Ny),intent(in)::U
		integer,intent(in)::Nx,Ny,N
		real*8, intent(in)::dx,dy
		integer::i,j
		character(len=20)::F_NAME

		if (N<10) then
			F_NAME='fichier/T'
			write(F_NAME (10:10),'(I1)') N
			F_NAME(11:14)= '.dat'

		elseif ((N>=10).and.(N<100)) then
			F_NAME='fichier/T'
			write(F_NAME (10:11),'(I2)') N
			F_NAME(12:15)= '.dat'

		elseif ((N>=100).and.(N<1000)) then
			F_NAME='fichier/T'
			write(F_NAME (10:12),'(I3)') N
			F_NAME(13:16)= '.dat'

		else if ((N>=1000).and.(N<10000)) then
			F_NAME='fichier/T'
			write(F_NAME (10:13),'(I4)') N
			F_NAME(14:17)= '.dat'

		else if ((N>=10000).and.(N<100000)) then
			F_NAME='fichier/T'
			write(F_NAME (10:14),'(I5)') N
			F_NAME(15:18)= '.dat'
		end if

		open(unit=2, file=F_NAME, action="write")


		do i=1,Nx
			do j=1,Ny


				write(2,*),i*dx,j*dy,U(bij(i,j,Ny))

			end do
			write(2,*)
		end do
		close(2)
		print*, 'N=',N
	end subroutine

	subroutine gradc(a,b,c,epsilon,Niter,Nx,Ny,X,V,h)
		implicit none
		real*8,intent(in)::a,b,c,epsilon,h
		integer,intent(in)::Niter,Nx,Ny
		real*8,dimension(:),intent(in)::V
		real*8,dimension(size(V)),intent(out)::X
		real*8,dimension(size(V))::R,D,W
		real*8::alpha,beta
		integer::iter

		X = 10
		call prod_mat_vec(a,b,c,X,R,Nx,Ny,h)
		R=R-V
		D=R
		iter = 1
	
		do while ((norme(R)>epsilon).AND.(iter<Niter))
			call prod_mat_vec(a,b,c,D,W,Nx,Ny,h)

			alpha=dot_product(D,R)/dot_product(D,W)

			X=X-alpha*D

			beta=1/(norme(R)**2)
			R=R-alpha*W
			beta=beta*norme(R)**2

			D=R+beta*D
			iter=iter+1
			!print*,'iter',iter
		end do
	end subroutine



	subroutine read_param(filename, Nx, Ny, Lx, Ly, rho,cp,lambda)
		implicit none
		character(len=*), intent(in)::filename
		integer, intent(out)::Nx,Ny
		real*8, intent(out)::Lx, Ly, rho,cp,lambda
		character(len=3)::bfr
		open(11, file=filename, action="read", status="old")
		read(11,'(A3,I6)') bfr, Nx ! "Nx="
		read(11,'(A3,I6)') bfr, Ny ! "Ny="
		read(11,'(A3,F4.6)') bfr, Lx ! "Lx="
		read(11,'(A3,F4.6)') bfr, Ly ! "Ly="
		read(11,'(A4,F10.6)') bfr, rho ! "rho="
		read(11,'(A3,F10.6)') bfr, cp ! "cp="
		read(11,'(A7,F4.6)') bfr, lambda ! "lambda"
	end subroutine

	subroutine compute_vars(Nx, Ny, Lx, Ly,rho,cp,lambda,Tmax,Niter,dt, dx, dy, a, b, c,D)
		implicit none
		integer, intent(in)::Nx, Ny, Niter
		real*8, intent(in)::Lx, Ly,Tmax,rho,cp,lambda
		real*8::dx, dy, a, b, c,dt,D
		dt = real(Tmax)/Niter
		dx = Lx/(Nx+1)
		dy = Ly/(Ny+1)
		D=lambda/(rho*cp)
		a = (1.+(2.*D*dt)/(dx*dx)+(2.*D*dt)/(dy*dy))
		b = -(D*dt)/(dy*dy)
		c = -(D*dt)/(dx*dx)
	end subroutine
end module
