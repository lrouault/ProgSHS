module parallel
	use sequentiel
contains


	subroutine prod_mat_vec_p(a,b,c,x,y,Nx,Ny,i1,im)
		implicit none
		real*8, intent(in)::a,b,c
		real*8, intent(in), dimension(Nx*Ny)::x
		real*8, intent(inout), dimension(Nx*Ny)::y
		integer,intent(in)::Nx,Ny,i1,im
		integer::i,j
		
		do i = i1, im
			if(i <= Ny) then
				if(modulo(i,Ny)==1) then
					y(i) = a*x(i)+b*x(i+1)+c*x(i+Ny)
				else if(modulo(i,Ny)==0) then
					y(i) = b*x(i-1)+a*x(i)+c*x(i+Ny)
				else
					y(i) = b*x(i-1)+a*x(i)+b*x(i+1)+c*x(i+Ny)
				end if
			else if(i >= Ny*(Nx-1)+1) then
				if(modulo(i,Ny)==1) then
					y(i) = c*x(i-Ny)+a*x(i)+b*x(i+1)
				else if(modulo(i,Ny)==0) then
					y(i) = c*x(i-Ny)+b*x(i-1)+a*x(i)
				else
					y(i) = c*x(i-Ny)+b*x(i-1)+a*x(i)+b*x(i+1)
				end if
			else
				if(modulo(i,Ny)==1) then
					y(i) = c*x(i-Ny)+a*x(i)+b*x(i+1)+c*x(i+Ny)
				else if(modulo(i,Ny)==0) then
					y(i) = c*x(i-Ny)+b*x(i-1)+a*x(i)+c*x(i+Ny)
				else
					y(i) = c*x(i-Ny)+b*x(i-1)+a*x(i)+b*x(i+1)+c*x(i+Ny)
				end if
			end if
			print *,i,y(i)
		end do
	end subroutine
	
	function dot_prod_p(A,B,i1,im)
		implicit none
		real*8, dimension(:), intent(in)::A,B
		integer,intent(in)::i1,im
		real*8::dot_prod_p
		do i=i1,im
			dot_prod_p = dot_prod_p + A(i)*B(i)
		end do
	end function

	subroutine gradc_p(a,b,c,epsilon,Niter,Nx,Ny,X,V,i1,im)
		implicit none
		real*8,intent(in)::a,b,c,epsilon
		integer,intent(in)::Niter,Nx,Ny,i1,im
		real*8,dimension(:),intent(in)::V
		real*8,dimension(size(V)),intent(out)::X
		real*8,dimension(size(V))::R,D,W
		real*8::alpha,beta,nr
		integer::iter

		X = 10
		call prod_mat_vec(a,b,c,X,R,Nx,Ny)
		R=R-V
		D=R
		iter = 1
	
		do while ((norme(R)>epsilon).AND.(iter<Niter))
			call prod_mat_vec(a,b,c,D,W,Nx,Ny)

			alpha=dot_product(D,R)/dot_product(D,W)

			X=X-alpha*D

			beta=1/(norme(R)**2)
			R=R-alpha*W
			beta=beta*norme(R)**2

			D=R+beta*D
			iter=iter+1
			!print*,'iter',iter
			
			nr = dot_prod_p(R(i1:im),R(i1:im),i1,im)
			
		end do
	end subroutine
	
	subroutine charge(n,np,me,i1,im)  
		implicit none
		integer,intent(in):: n,np,me
		integer,intent(out)::i1,im
		integer :: rep,q

		q=n/np
		rep=n-q*np

		if (me<rep) then
			i1=(q+1)*me+1
			im=(me+1)*(q+1)
		else
			i1=q*me+1+rep
			im=i1+q-1
		end if
	end subroutine
end module
