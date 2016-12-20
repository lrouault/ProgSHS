!> Déclaration de toutes las variables globales utilisées dans le programme

module donnees

  implicit none
  include "mpif.h"

  integer, parameter :: PR=Selected_real_kind(12,60) !12,60

  ! Maillage
  integer   :: Nx !< Nombre de points dans la direction x
  integer   :: Ny !< Nombre de points dans la direction y
  integer   :: Nz !< Nombre de points dans la direction z
  integer   :: Niter
  integer   :: Ndisplay !< Frequence de sauvegarde des solutions
  real(PR)  :: Lx !< Dimension selon x
  real(PR)  :: Ly !< Dimension selon y
  real(PR)  :: Lz !< Dimension selon z
  real(PR)  :: dx !< Pas d'espace principale en x
  real(PR)  :: dy !< Pas d'espace principale en y
  real(PR)  :: dz !< Pas d'espace principale en z
  real(PR)  :: x
  real(PR)  :: y
  real(PR)  :: z
  real(PR)  :: time !< Temps du calcul
  real(PR)  :: dt !< Pas de temps pour les calculs
  real(PR)  :: tmax !< Condition d'arret du calcul en temps

  ! Thermique
  real(PR)  :: Text !< Temperature initiale
  real(Pr)  :: Tad
  !real(PR)  :: rho, cp, lambda
  real(PR),dimension(:),allocatable   :: rho, rhocp, lambda
  real(PR)  :: h, flux
  real(PR)  :: Ea, R, k0, Q, coef

  ! Modèle
  real(PR), dimension(:), allocatable :: U, U0, rhs, Fm, eta, chi, eta0
  real(PR)                            :: epsilon
  real(PR), dimension(:), allocatable :: Cd
  real(PR), dimension(:), allocatable :: Cx, Cy, Cz

  integer :: nb_fichiers

  ! Materiau
  real(PR),dimension(:,:),allocatable :: fraction_vol ! (/Si,N2,fibre/)
  real(PR)                            :: rho_Si, rho_Si3N4, rho_N2, rho_fibre
  real(PR),dimension(:,:),allocatable :: cp_Si, cp_Si3N4, cp_N2, cp_fibre
  real(PR),dimension(:,:),allocatable :: lambda_Si, lambda_Si3N4, lambda_N2, lambda_fibre

  ! fibre
  integer, dimension(:,:,:),allocatable :: Poro
  integer, dimension(:,:,:,:),allocatable :: Orien
  integer :: porox, poroy, poroz

  integer :: iteration

  !Parallèle
  integer :: statinfo,Np,Me,num1,numN,k1,kN
  integer,dimension(MPI_STATUS_SIZE) :: status

contains

end module donnees
