! TESIS CODE 2.1
! 
! Name: Gomez Yanez Uzmar de Jesus
! Account Number: 308287141
! Email: uzmar.gomez@ciencias.unam.mx
! Tutor: Miguel Alcubierre Moya
! Git repository: http
! 
!******************************************************************************
! 
! The following program seeks the resolution of the spherically 
! symmetric Vlasov equation, in the case where all particles have
! same angular momentum L, for the relativistic case. We will use 
! two kind of coordinates for the Schwarzschild geometry (a spherically
! symmetric Black Hole), the Schwarzschild coordinates and the 
! Kerr-Schild coordinates.
! 
! In order to run this program, gnuplot must be installed.
! Also, file "VariablesTesis.txt" has to be in the same carpet as the
! program itself.
! 
!******************************************************************************

module integrals
contains
!------------------------------------------------------------------------------
!                                 Integrals 
!------------------------------------------------------------------------------
! 
! Integral of a one dimensional function of one variable
! 
subroutine int1d1v(f,Nx,dx,integral)
   real(8),dimension(:),allocatable:: f
   real(8):: integral
   real(8):: dx
   real(8):: auxiliar1,auxiliar2
   integer:: Nx
   integer:: j
   
   integral=0.0d0
   auxiliar1=0.0d0
   auxiliar2=0.0d0   

   do j=1,(Nx/2)-1
      auxiliar1=auxiliar1+f(2*j)
   end do
   do j=1,Nx/2
      auxiliar2=auxiliar2+f(2*j-1)
   end do
   integral = (dx/3.0d0)*(f(0)+2.0d0*auxiliar1+4.0d0*auxiliar2+f(Nx))
   
end subroutine


! Integral of a two dimensional function in one variable

subroutine int1d(f,Nx,Ny,dy,integral)
   real(8),dimension(:,:),allocatable:: f
   real(8),dimension(:),allocatable:: integral
   real(8):: dy
   real(8),allocatable,dimension(:):: auxiliar1,auxiliar2
   integer:: Nx,Ny
   integer:: i,j
   
   allocate(auxiliar1(0:Nx))
   allocate(auxiliar2(0:Nx))

   integral=0.0d0
   auxiliar1=0.0d0
   auxiliar2=0.0d0   

   do i=0,Nx
      do j=1,(Ny/2)-1
         auxiliar1(i)=auxiliar1(i)+f(i,2*j)
      end do
      do j=1,Ny/2
         auxiliar2(i)=auxiliar2(i)+f(i,2*j-1)
      end do
   integral(i) = (dy/3.0d0)*(f(i,0)+2.0d0*auxiliar1(i)+4.0d0*auxiliar2(i)+f(i,Ny))
   end do

   deallocate(auxiliar1)
   deallocate(auxiliar2)
end subroutine


! Integral of a two dimensional function in two variables

subroutine int2d(f,Nx,Ny,dx,dy,integral)

   real(8),dimension(:,:),allocatable:: f
   real(8):: dx
   real(8):: dy
   real(8):: integral,aux1,aux2,aux3
   integer:: Nx
   integer:: Ny
   integer:: i,j

   integral=0.0d0
   aux1=0.0d0
   aux2=0.0d0
   aux3=0.0d0

   do i=1,Nx-1
      aux1=aux1+f(i,0)+f(i,Ny)
   end do
   do j=1,Ny-1
      aux2=aux2+f(0,j)+f(Nx,j)
   end do
   do i=1,Nx-1
      do j=1,Ny-1
         aux3=aux3+f(i,j)
      end do
   end do
   integral=integral+0.25d0*dx*dy*(f(0,0)+f(Nx,0)+f(0,Ny)+f(Nx,Ny)+ 2.0d0*(aux1+aux2)+ 4.0d0*aux3) 
end subroutine


!-------------------------------------------------------------------------------
!                         Find boundary Fluxes
! 
! For the integrals of the fluxes at the boundaries we
! must remember that the fluxes in a given direction are
! evaluated at the half integer posisitions, so that to
! get the correct values at the boundaries we must interpolate.
! 
!-------------------------------------------------------------------------------
subroutine boundary(dNdt,fluxr,fluxp,Nr,Np,dr,dp,rmin)

   implicit none

   real(8),dimension(:,:),allocatable:: fluxr,fluxp
   real(8):: dr,dp
   real(8):: dNdt,intflux_rmin,intflux_rmax,intflux_pmin,intflux_pmax
   real(8):: rmin

   integer:: Np,Nr
   integer:: i,j

! Integrated flux on rmin.  Do notice that if rmin=0 then
! there isn't really a boundary there, it is just the origin.

   if (rmin==0.d0) then
      intflux_rmin = 0.d0
   else
      intflux_rmin = 0.25d0*(fluxr(0,0 ) + fluxr(1,0 ) &
                           + fluxr(0,Np) + fluxr(1,Np))
      do j=1,Np-1
         intflux_rmin = intflux_rmin + 0.5d0*(fluxr(0,j) + fluxr(1,j))
      end do
      intflux_rmin = intflux_rmin*dp
   end if


! Integrated flux on rmax.

   intflux_rmax = 0.25d0*(fluxr(Nr-1,0 ) + fluxr(Nr,0 ) &
                        + fluxr(Nr-1,Np) + fluxr(Nr,Np))

   do j=1,Np-1
   intflux_rmax = intflux_rmax + 0.5d0*(fluxr(Nr-1,j) + fluxr(Nr,j))
   end do

  intflux_rmax = intflux_rmax*dp


! Integrated flux on pmin.

   intflux_pmin = 0.25d0*(fluxp(0 ,0) + fluxp(0 ,1) &
                        + fluxp(Nr,0) + fluxp(Nr,1))

   do i=0,Nr-1
      intflux_pmin = intflux_pmin + 0.5d0*(fluxp(i,0) + fluxp(i,1))
   end do

   intflux_pmin = intflux_pmin*dr


! Integrated flux on pmax.

   intflux_pmax = 0.25d0*(fluxp(0 ,Np-1) + fluxp(0 ,Np) &
                        + fluxp(Nr,Np-1) + fluxp(Nr,Np))

   do i=0,Nr-1
      intflux_pmax = intflux_pmax + 0.5d0*(fluxp(i,Np-1) + fluxp(i,Np))
   end do

   intflux_pmax = intflux_pmax*dr


! Total flux of particles through all boundaries.
! The sign convention here is that dNdt is positive
! if the particles enter the domain, and negative
! if they leave it.

   dNdt = (intflux_rmin - intflux_rmax) &
        + (intflux_pmin - intflux_pmax)
end subroutine
end module 




module flux
contains
!------------------------------------------------------------------------------
!                           Calculate Fluxes 
! 
! Calculation of fluxes for the Vlasov Equation:
!    
! flux_r=f_{*}(dr/dt)
! flux_p=f_{*}(dp/dt).
! 
! Because of the conservation form of Vlasov Equation, and to avoid 
! spurious oscillations, we should use a High Resolution Method. First, we need
! to write the above fluxes (for simplicity, we should express a flux in general
! as F in this little paragraph) as F=F_H-(1-g(s))(F_H-F_L),
! where F_L is a low order flux in smooth regions (upwind first order), and F_H
! is a high order flux from a monotonic method near a discontinuity (Lax Wendroff).
! This kind of methods are known as flux "limiters", and we should use them to 
! calculate flux_r and flux_p
!.
!------------------------------------------------------------------------------
subroutine fluxes(fluxr,fluxp,f,rho_r,rho_p,drdt,dpdt,Np,Nr,lim)
   implicit none
   real(8),dimension(:,:),allocatable:: f             ! Initial density function
   integer:: Nr                                       ! Number of grid points in r
   integer:: Np                                       ! Number of grid points in p_r
   real(8),dimension(:,:),allocatable:: dpdt          ! dp_r/dt
   real(8),dimension(:,:),allocatable:: drdt          ! dr/dt    
   real(8),dimension(:,:),allocatable:: fluxr         ! flux in r direction:  dr/dt*f   
   real(8),dimension(:,:),allocatable:: fluxp         ! flux in p direction:  dp/dt*f
   integer:: j,k                                      ! Counters
   character(50):: lim                                ! Limiter
   real(8):: slope1                                   ! f_m - f_{m-1}
   real(8):: slope2                                   ! f_{m+1} - f_m
   real(8):: slopecorrect
   real(8):: rho_r                                    ! Courant factor for r
   real(8):: rho_p                                    ! Courant factor for p
      fluxr=0.0d0
      fluxp=0.0d0
      do j=0,Nr
         do k=0,Np
            if (drdt(j,k)>=0) then
               if (j==0) then
                  fluxr(j,k)=drdt(j,k)*f(j,k)+ 0.5d0*drdt(j,k)*(1.0d0-drdt(j,k)*rho_r)*(f(j+1,k)-f(j,k))
               else if (j==Nr) then
                  fluxr(j,k)=drdt(j,k)*f(j,k)+ 0.5d0*drdt(j,k)*(1.0d0-drdt(j,k)*rho_r)*(f(j,k)-f(j-1,k))
               else
                  slope1 = (f(j+1,k)-f(j,k))
                  slope2 = (f(j,k)-f(j-1,k))
                  call limiter(slope1,slope2,slopecorrect,lim)
                  fluxr(j,k)=drdt(j,k)*f(j,k)+ 0.5d0*drdt(j,k)*(1.0d0-drdt(j,k)*rho_r)*slopecorrect
               end if
            else
               if (j==Nr) then
                  fluxr(j,k)=drdt(j,k)*f(j,k)+ 0.5d0*drdt(j,k)*(-1.0d0-drdt(j,k)*rho_r)*(f(j,k)-f(j-1,k))
               else if (j==0) then
                  fluxr(j,k)=drdt(j,k)*f(j,k)+ 0.5d0*drdt(j,k)*(-1.0d0-drdt(j,k)*rho_r)*(f(j+1,k)-f(j,k))
               else 
                  slope1 = (f(j+1,k)-f(j,k))
                  slope2 = (f(j,k)-f(j-1,k))
                  call limiter(slope1,slope2,slopecorrect,lim)
                  fluxr(j,k)=drdt(j,k)*f(j,k)+ 0.5d0*drdt(j,k)*(-1.0d0-drdt(j,k)*rho_r)*slopecorrect
               end if
            end if
         end do
      end do
      do j=0,Nr
         do k=0,Np
            if (dpdt(j,k)>=0) then
               if (k==Np) then
                  fluxp(j,k)=dpdt(j,k)*f(j,k)+ 0.5d0*dpdt(j,k)*(1.0d0-dpdt(j,k)*rho_p)*(f(j,k)-f(j,k-1))
               else if (k==0) then
                  fluxp(j,k)=dpdt(j,k)*f(j,k)+ 0.5d0*dpdt(j,k)*(1.0d0-dpdt(j,k)*rho_p)*(f(j,k+1)-f(j,k))
               else
                  slope1 = (f(j,k+1)-f(j,k))
                  slope2 = (f(j,k)-f(j,k-1))
                  call limiter(slope1,slope2,slopecorrect,lim)
                  fluxp(j,k)=dpdt(j,k)*f(j,k)+ 0.5d0*dpdt(j,k)*(1.0d0-dpdt(j,k)*rho_p)*slopecorrect
               end if
            else
               if (k==Np) then
                  fluxp(j,k)=dpdt(j,k)*f(j,k)+ 0.5d0*dpdt(j,k)*(-1.0d0-dpdt(j,k)*rho_p)*(f(j,k)-f(j,k-1))
               else if (k==0) then
                  fluxp(j,k)=dpdt(j,k)*f(j,k)+ 0.5d0*dpdt(j,k)*(-1.0d0-dpdt(j,k)*rho_p)*(f(j,k+1)-f(j,k))
               else
                  slope1 = (f(j,k+1)-f(j,k))
                  slope2 = (f(j,k)-f(j,k-1))
                  call limiter(slope1,slope2,slopecorrect,lim)
                  fluxp(j,k)=dpdt(j,k)*f(j,k)+ 0.5d0*dpdt(j,k)*(-1.0d0-dpdt(j,k)*rho_p)*slopecorrect
               end if
            end if
         end do
      end do
return
end subroutine fluxes
end module flux


!******************************************************************************
!*                                                                            *
!*                              MAIN PROGRAM                                  *
!*                                                                            *
!******************************************************************************

program progtesis

!------------------------------------------------------------------------------
!                         Parameters and constants

   use flux
   use integrals

   implicit none

   real(8),dimension(:),allocatable:: r               ! Position
   real(8),dimension(:),allocatable:: p               ! Momentum p_r
   real(8):: L                                        ! Angular momentum
    
   real(8):: gaussian,compact                         ! Kind of density function
   real(8),dimension(:,:),allocatable:: f             ! Density function
   real(8),dimension(:,:),allocatable:: g             ! Auxiliar
   real(8),dimension(:,:),allocatable:: f_old         ! Old value of f

   real(8):: tmax,rmin,rmax,pmin,pmax 
   real(8):: dt                                       ! Step size in t
   real(8):: dr                                       ! Step size in r
   real(8):: dp                                       ! Step size in p_r
   real(8):: rho_r                                    ! Courant factor for r
   real(8):: rho_p                                    ! Courant factor for p
    
   real(8):: rc,pc                                    ! Initial conditions for f
   real(8):: N0                                       ! Number of particles
   real(8):: sr,sp                                    ! Constants for initial f
    
   real(8):: Msch,plusr,plusp,plusL                   ! Schwarzschild mass,+pc,+rc,+Lc
   real(8):: m                                        ! Mass of the shell
    
   integer:: Nt                                       ! Number of grid points in t
   integer:: Nr                                       ! Number of grid points in r
   integer:: Np                                       ! Number of grid points in p_r
    
   real(8),dimension(:,:),allocatable:: p0            ! Momentum temporal coordinate p^0
   real(8),dimension(:,:),allocatable:: dpdt          ! dp_r/dt
   real(8),dimension(:,:),allocatable:: drdt          ! dr/dt

!*************************** 3+1 Quantities ****************************** 

   real(8),dimension(:),allocatable:: A
   real(8),dimension(:),allocatable:: dAdr
   real(8),dimension(:),allocatable:: alpha
   real(8),dimension(:),allocatable:: dalphadr
   real(8),dimension(:),allocatable:: betasubr
   real(8),dimension(:),allocatable:: betasupr
   real(8),dimension(:),allocatable:: dbetadr
   real(8):: B,psi
 
!*************************************************************************

   real(8),dimension(:,:),allocatable:: fluxr         ! flux in r direction:  dr/dt*f   
   real(8),dimension(:,:),allocatable:: fluxp         ! flux in p direction:  dp/dt*f
    
   integer:: i,j,n                                    ! Counters
   integer:: stepr,stepp,stept                        ! Steps for storing data 

   real(8):: dF,dG,dNdt
   real(8),dimension(:,:),allocatable:: k1,k2,k3,k4               
   real(8),dimension(:,:),allocatable:: aux1,aux2,aux3

   character(50):: lim                                ! Limiter (minmod,mclimiter,superbee)
   character(50):: coord                              ! Coordinates for Schwarzschild
   character(50):: initial                            ! Kind of density function (gaussian, compact)

   real(8),dimension(:),allocatable:: density         ! Density   
   real(8):: Npart,Npart_old                          ! Number of particles

   character(100):: generaldatacarpet                 ! Data     
   character(30)::  filenameevolution				      ! Name of the multiple evolution data files
   character(30)::  filenamedensity
   character(30)::  evolutiondatacarpet				   ! Evolution data carpet
   character(30)::  densitydatacarpet
   character(30)::  convergencecarpet                 ! Convergence carpet

   real(8):: S                                        ! For normalization
   real:: start, finish                               ! For execution time of the program

!******************** Stability testing (Virial Theorem) *************************
   real(8),dimension(:,:),allocatable:: T             ! Kinetic Energy
   real(8),dimension(:,:),allocatable:: V             ! Virial
   real(8):: tenergy,venergy                                 

!*************************** Convergence *****************************************
   integer(8)::percentage                             ! convnumb*dt, % of time elapsed to check convergence

   real(8):: pi
    
   pi=acos(-1.0d0)

   call cpu_time(start)

!------------------------------------------------------------------------------
!                               Reading Data

   open(10,file="Variables.txt")
   read(10,*)
   read(10,*)
   read(10,*)
   read(10,*)
   read(10,*)
   read(10,*)
   read(10,*)
   read(10,*)
   read(10,*)
   read(10,*)
   read(10,*)
   read(10,*)
   read(10,*)
   read(10,*) lim
   read(10,*) coord
   read(10,*) initial
   read(10,*)
   read(10,*) tmax
   read(10,*) rmax
   read(10,*) rmin
   read(10,*) pmax
   read(10,*) pmin
   read(10,*)
   read(10,*) dt
   read(10,*) rho_r   
   read(10,*) rho_p
   read(10,*)
   read(10,*)
   read(10,*) stept
   read(10,*) stepr
   read(10,*) stepp
   read(10,*)
   read(10,*)
   read(10,*) percentage
   read(10,*)
   read(10,*)
   read(10,*) N0
   read(10,*) Msch
   read(10,*) m
   read(10,*)
   read(10,*) plusr
   read(10,*) plusp
   read(10,*) plusL
   read(10,*) sr
   read(10,*) sp
   close(10)   
   
!------------------------------------------------------------------------------   
!                              Create Carpets

   call execute_command_line("mkdir -p Data/Graphs/GraphsEvolution")
   call execute_command_line("mkdir -p Data/Graphs/GraphsDensity")   
   call execute_command_line("mkdir -p Data/EvolutionDataFiles")
   call execute_command_line("mkdir -p Data/ConvergenceFiles")
   call execute_command_line("mkdir -p Data/DensityDataFiles")
   call execute_command_line("cp Variables.txt Data/Variables.txt")
   write(generaldatacarpet,*) "Data/"
   write(evolutiondatacarpet,*)  "Data/EvolutionDataFiles/"
   write(densitydatacarpet,*) "Data/DensityDataFiles/"
   write(convergencecarpet,*)  "Data/ConvergenceFiles/"
   
!------------------------------------------------------------------------------

   Nt=int(tmax/dt)
   dr=dt/rho_r
   dp=dt/rho_p

   Nr=int((rmax-rmin)/dr)
   Np=int((pmax-pmin)/dp)
   
!------------------------------------------------------------------------------
!                   Quantities for stable circular orbits
   if (coord=='sch') then

      L=plusL+sqrt(12.0d0)*Msch
      if (L>sqrt(12.0d0)*Msch) then
         rc=plusr+(L**2 + sqrt(L**4 -12.0d0*(L*Msch)**2))/(2.0d0*Msch)
         pc=plusp+(1.0d0/sqrt(1.0d0-3.0d0*Msch/rc))*(2.0d0*Msch/rc)
      else if (L==sqrt(12.0d0)*Msch) then
         rc=plusr+6.0d0*Msch
         pc=plusp+(1.0d0/sqrt(1.0d0-3.0d0*Msch/rc))*(2.0d0*Msch/rc)
      else if (L<sqrt(12.0d0)*Msch) then
         rc=plusr
         pc=plusp
      end if

      rmin=2.1d0*Msch    

   else if (coord=='ks') then

      L=plusL+sqrt(12.0d0)*Msch
      if (L>=sqrt(12.0d0)*Msch) then
         rc=plusr+(L**2 + sqrt(L**4 -12.0d0*(L*Msch)**2))/(2.0d0*Msch)
         pc=plusp+(1.0d0/sqrt(1.0d0-3.0d0*Msch/rc))*(2.0d0*Msch/rc)
      else if (L==sqrt(12.0d0)*Msch) then
         rc=plusr+6.0d0*Msch
         pc=plusp+(1.0d0/sqrt(1.0d0-3.0d0*Msch/rc))*(2.0d0*Msch/rc)
      else if (L<sqrt(12.0d0)*Msch) then
         rc=plusr
         pc=plusp
      end if
      
      rmin=Msch/4.0d0

   else
      print*
      print*, 'Unknown kind of coordinates'
      print*, 'Aborting ...'
      print*
      stop
   end if

!------------------------------------------------------------------------------
!                             Allocate arrays

   allocate(r(0:Nr))
   allocate(p(0:Np))
   allocate(f(0:Nr,0:Np))
   allocate(g(0:Nr,0:Np))
   allocate(f_old(0:Nr,0:Np))
   allocate(T(0:Nr,0:Np))
   allocate(p0(0:Nr,0:Np))
   allocate(dpdt(0:Nr,0:Np))
   allocate(drdt(0:Nr,0:Np))
   allocate(A(0:Nr))
   allocate(dAdr(0:Nr))
   allocate(betasubr(0:Nr))
   allocate(betasupr(0:Nr))
   allocate(dbetadr(0:Nr))
   allocate(alpha(0:Nr))
   allocate(dalphadr(0:Nr))
   allocate(fluxr(0:Nr,0:Np))
   allocate(fluxp(0:Nr,0:Np))
   allocate(density(0:Nr))
   allocate(V(0:Nr,0:Np))
   allocate(k1(0:Nr,0:Np))
   allocate(k2(0:Nr,0:Np))      
   allocate(k3(0:Nr,0:Np))
   allocate(k4(0:Nr,0:Np))
   allocate(aux1(0:Nr,0:Np))
   allocate(aux2(0:Nr,0:Np))
   allocate(aux3(0:Nr,0:Np))

   r=         0.0d0
   p=         0.0d0
   f=         0.0d0
   g=         0.0d0
   f_old=     0.0d0
   T=         0.0d0
   p0=        0.0d0
   dpdt=      0.0d0
   drdt=      0.0d0
   A=         0.0d0
   dAdr=      0.0d0
   betasubr=  0.0d0
   betasupr=  0.0d0
   dbetadr=   0.0d0
   alpha=     0.0d0
   dalphadr=  0.0d0
   V=         0.0d0   
   density=   0.0d0


!------------------------------------------------------------------------------
!                           Filling arrays r,p_r

   if (rmin==0.0d0) then
      do i=0,Nr
         r(i)=-0.5d0*dr+dble(i)*dr
      end do
   else
      do i=0,Nr
         r(i)=rmin+dble(i)*dr
      end do
   end if

   do i=0,Np
      p(i)=pmin+dble(i)*dp
   end do
    
    
!------------------------------------------------------------------------------
!                         Initial density function f

   if (initial=='gaussian') then
      do i=0,Nr
         do j=0,Np
            g(i,j)=gaussian(r(i),p(j),rc,pc,sr,sp)
         end do
      end do
   else if (initial=='compact') then
      do i=0,Nr
         do j=0,Np
            g(i,j)=compact(r(i),p(j),rc,pc,sr,sp)
         end do
      end do
   end if

   call int2d(g,Nr,Np,dr,dp,S)
   do i=0,Nr
      do j=0,Np
         f_old(i,j)=(N0/S)*g(i,j)
      end do
   end do

!------------------------------------------------------------------------------
!          Spherical coordinate choice (Schwarzschild or Kerr-Schild)

if (coord=='sch') then
   do i=0,Nr
      alpha(i)=sqrt(1.0d0-2.0d0*Msch/r(i))
      dalphadr(i)=Msch/(sqrt(1.0d0-2.0d0*Msch/r(i))*r(i)**2)
      A(i)=(1.0d0-2.0d0*Msch/r(i))**(-1)
      dAdr(i)=-2.0d0*Msch/((r(i)-2.0d0*Msch)**2)
      betasubr(i)=0.0d0
      betasupr(i)=0.0d0
      dbetadr(i)=0.0d0
   end do
   B=1.0d0
   psi=1.0d0
   
else if (coord=='ks') then
   do i=0,Nr
      alpha(i)=(1.0d0+2.0d0*Msch/r(i))**(-0.5)
      dalphadr(i)=(Msch/r(i)**2)*(1.0d0+2.0d0*Msch/r(i))**(-1.5)
      A(i)=1.0d0+2.0d0*Msch/r(i)
      dAdr(i)=-2.0d0*Msch/r(i)**2
      betasubr(i)=2.0d0*Msch/r(i)
      betasupr(i)=(2.0d0*Msch/r(i))/(1.0d0+2.0d0*Msch/r(i))
      dbetadr(i)=-2.0d0*Msch/(2.0d0*Msch+r(i))**2
   end do
   B=1.0d0
   psi=1.0d0
end if      


!------------------------------------------------------------------------------
!                          p^0, dp/dt and dr/dt

do i=0,Nr
   do j=0,Np
      p0(i,j)=(1.0d0/alpha(i))*sqrt(m**2 + (1.0d0/(A(i)*psi**4))*p(j)**2 &
             & + (L**2)/(r(i)**2*B*psi**4))
   end do
end do  
do i=0,Nr
   do j=0,Np
      dpdt(i,j)=-alpha(i)*dalphadr(i)*p0(i,j) + dbetadr(i)*p(j) &
               & + (dAdr(i)*p(j)**2)/(2.0d0*p0(i,j)*(A(i)**2)*(psi**4))&
               & + (L**2/(B*(r(i)**3)*(psi**4)*p0(i,j)))
   end do
end do
do i=0,Nr
   do j=0,Np
      drdt(i,j)=p(j)/(A(i)*p0(i,j)*psi**4)-betasupr(i)
   end do
end do

print*, '   ---------------------------------------------'
print*, '   | SOME USEFUL DATA TO BE TAKEN INTO ACCOUNT |'
print*, '   ---------------------------------------------'
print*
print*, 'Steps in:           r=',Nr 
print*, '                    p=',Np
print*, '                    t=',Nt
print*
print*, 'Max velocity:       max(drdt)=',maxval(drdt)
print*, '                    max(dpdt)=',maxval(dpdt)
print*
print*, 'Courant number:     rho_r=',rho_r
print*, '                    rho_p=',rho_p
print*
print*, 'Grid spacing:       dr=',dr
print*, '                    dp=',dp
print*, '                    dt=',dt
print*
print*, 'CFL:                rho_r*max(drdt)=',rho_r*maxval(drdt)
print*, '                    rho_p*max(dpdt)=',rho_p*maxval(dpdt)
if (max(abs(rho_r*maxval(drdt)),abs(rho_p*maxval(dpdt)))<=1.0d0/sqrt(2.0d0)) then
   print*
   print*, 'CFL condition satisfied!'
else
   print*
   print*, 'CFL condition not satisfied! Try put different values of rho_p and rho_r'
   print*, 'Aborting...'
   stop
end if
print*
call sleep(2)

!-----------------------------------------------------------------------------
!                      Save and graph function f at t=0


   open(16,file=trim(adjustl(generaldatacarpet))//'vlasov.dat',status='unknown')
      do j=0,Np
         do i=0,Nr
            if (mod(i,stepr)==0) then
               write(16,*) r(i),p(j),f_old(i,j)
            end if
         end do
         write(16,*) 
      end do
    
   open(17,file=trim(adjustl(generaldatacarpet))//'fold.dat',status='unknown')
      do j=0,Np
         do i=0,Nr
            write(17,*) r(i),p(j),f_old(i,j)
         end do 
      end do
   close(17)

!-----------------------------------------------------------------------------
!                               Save dp/dt

   open(12,file=trim(adjustl(generaldatacarpet))//'dpdt.dat',status='unknown')
      do i=0,Nr
         do j=0,Np
            write(12,*) r(i),p(j),dpdt(i,j)
         end do
      end do
   close(12)
  
!-----------------------------------------------------------------------------
!                               Save dr/dt

   open(13,file=trim(adjustl(generaldatacarpet))//'drdt.dat',status='unknown')
      do i=0,Nr
         do j=0,Np
            write(13,*) r(i),p(j),drdt(i,j)
         end do
      end do
   close(13)

!-----------------------------------------------------------------------------
!                                Save p0 

   open(14,file=trim(adjustl(generaldatacarpet))//'p0.dat',status='unknown')
      do i=0,Nr
         do j=0,Np
            write(14,*) r(i),p(j),p0(i,j)
         end do
      end do
   close(14)

   call execute_command_line("touch graphs")
   open(22,file='graphs')
   write(22,*) "set terminal png"
   write(22,*) "set nokey" 
   write(22,*) "set output 'Data/Graphs/initialf.png'"
   write(22,*) "set xlabel 'r'"
   write(22,*) "set ylabel 'p'"
   write(22,*) "splot 'Data/fold.dat' with lines lt rgb 'blue'"
   close(22)
   call execute_command_line("gnuplot graphs")
   call execute_command_line("rm graphs")


!------------------------------------------------------------------------------
!                  Initial density and number of particles

   call int1d(f_old,Nr,Np,dp,density)
   call int2d(f_old,Nr,Np,dr,dp,Npart)
   
   open(18,file=trim(adjustl(generaldatacarpet))//'density.dat',status='unknown')
      do j=0,Nr
         if (mod(j,stepr)==0) then
            write(18,*) r(j),density(j)/(r(j)**2)
         end if
      end do
      write(18,*)
      write(18,*)

   open(19,file=trim(adjustl(generaldatacarpet))//'numpart.dat',status='unknown')
      write(19,*) 0,Npart

      Npart_old=Npart

   open(23,file=trim(adjustl(generaldatacarpet))//'tenergy.dat',status='unknown')
   open(24,file=trim(adjustl(generaldatacarpet))//'venergy.dat',status='unknown')


!******************************************************************************
!                         Start Main Evolution Loop
! 
! For time integration I will use 4th Order Runge-Kutta. Method of lines is 
! also used to separate time and space.
! 
!******************************************************************************  
!******************************************************************************
!                                                   EVOLUTION
do n=1,Nt
        k1=0.0d0
        k2=0.0d0
        k3=0.0d0   
        k4=0.0d0      
        aux1=0.0d0
        aux2=0.0d0
        aux3=0.0d0
        fluxr=0.0d0
        fluxp=0.0d0
        density=0.0d0   

        call int1d(f_old,Nr,Np,dp,density)
        call int2d(f_old,Nr,Np,dr,dp,Npart)
        call fluxes(fluxr,fluxp,f_old,rho_r,rho_p,drdt,dpdt,Np,Nr,lim)
        call boundary(dNdt,fluxr,fluxp,Nr,Np,dr,dp,rmin)

        do j = 0,Np
                f_old(1,j) = f_old(2,Np-j)
                f_old(0,j) = f_old(3,Np-j)
        end do
        do i=0,Nr
                do j=0,Np					
                        if (drdt(i,j)>=0) then
                                if (i==0) then
                                        dF=drdt(i+1,j)*f_old(i+1,j)-drdt(i,j)*f_old(i,j)
                                else
                                        dF=fluxr(i,j)-fluxr(i-1,j)
                                end if
                        else
                                if (i==0) then
                                        dF=fluxr(i+1,j)-fluxr(i,j)
						else if (i==Nr) then
							dF=drdt(i,j)*f_old(i,j)-drdt(i-1,j)*f_old(i-1,j)
						else
							dF=fluxr(i+1,j)-fluxr(i,j)
						end if
					end if
					if (dpdt(i,j)>=0) then
						if (j==0) then
							dG=dpdt(i,j+1)*f_old(i,j+1)-dpdt(i,j)*f_old(i,j)
						else
							dG=fluxp(i,j)-fluxp(i,j-1)
						end if
					else
						if (j==0) then
							dG=fluxp(i,j+1)-fluxp(i,j)
						else if (j==Np) then
							dG=dpdt(i,j)*f_old(i,j)-dpdt(i,j-1)*f_old(i,j-1)
						else
							dG=fluxp(i,j+1)-fluxp(i,j)
						end if
					end if
					k1(i,j)=-(1.0d0/dr)*dF-(1.0d0/dp)*dG
					aux1(i,j) = f_old(i,j) + 0.5d0*dt*k1(i,j)
				end do
			end do
			call fluxes(fluxr,fluxp,aux1,rho_r,rho_p,drdt,dpdt,Np,Nr,lim)
			do i=0,Nr
				do j=0,Np					
					if (drdt(i,j)>=0) then
						if (i==0) then
							dF=drdt(i+1,j)*aux1(i+1,j)-drdt(i,j)*aux1(i,j)
						else
							dF=fluxr(i,j)-fluxr(i-1,j)
						end if
					else
						if (i==0) then
							dF=fluxr(i+1,j)-fluxr(i,j)
						else if (i==Nr) then
							dF=drdt(i,j)*aux1(i,j)-drdt(i-1,j)*aux1(i-1,j)
						else
							dF=fluxr(i+1,j)-fluxr(i,j)
						end if
					end if
					if (dpdt(i,j)>=0) then
						if (j==0) then
							dG=dpdt(i,j+1)*aux1(i,j+1)-dpdt(i,j)*aux1(i,j)
						else
							dG=fluxp(i,j)-fluxp(i,j-1)
						end if
					else
						if (j==0) then
							dG=fluxp(i,j+1)-fluxp(i,j)
						else if (j==Np) then
							dG=dpdt(i,j)*aux1(i,j)-dpdt(i,j-1)*aux1(i,j-1)
						else
							dG=fluxp(i,j+1)-fluxp(i,j)
						end if
					end if
					k2(i,j)=-(1.0d0/dr)*dF-(1.0d0/dp)*dG
					aux2(i,j)=f_old(i,j)+dt*(2.0d0*k2(i,j)-k1(i,j))
				end do
			end do
			call fluxes(fluxr,fluxp,aux2,rho_r,rho_p,drdt,dpdt,Np,Nr,lim)
			do i=0,Nr
				do j=0,Np					
					if (drdt(i,j)>=0) then
						if (i==0) then
							dF=drdt(i+1,j)*aux2(i+1,j)-drdt(i,j)*aux2(i,j)
						else
							dF=fluxr(i,j)-fluxr(i-1,j)
						end if
					else
						if (i==0) then
							dF=fluxr(i+1,j)-fluxr(i,j)
						else if (i==Nr) then
							dF=drdt(i,j)*aux2(i,j)-drdt(i-1,j)*aux2(i-1,j)
						else
							dF=fluxr(i+1,j)-fluxr(i,j)
						end if
					end if
					if (dpdt(i,j)>=0) then
						if (j==0) then
							dG=dpdt(i,j+1)*aux2(i,j+1)-dpdt(i,j)*aux2(i,j)
						else
							dG=fluxp(i,j)-fluxp(i,j-1)
						end if
					else
						if (j==0) then
							dG=fluxp(i,j+1)-fluxp(i,j)
						else if (j==Np) then
							dG=dpdt(i,j)*aux2(i,j)-dpdt(i,j-1)*aux2(i,j-1)
						else
							dG=fluxp(i,j+1)-fluxp(i,j)
						end if
					end if
					k3(i,j)=-(1.0d0/dr)*dF-(1.0d0/dp)*dG
				end do
			end do

			do i=0,Nr
				do j=0,Np
					if (i==Nr) then
                                                if (p(j)<=0) then
                                                        f(i,j)=0.0d0
                                                else
                                                f(i,j)=f_old(i,j)+(dt/6.0d0)*(k1(i,j)+4.0D0*k2(i,j)+k3(i,j))
                  e nd if
					else if (i==0) then
						if (p(j)>=0) then
							f(i,j)=0.0d0
						else
							f(i,j)=f_old(i,j)+(dt/6.0d0)*(k1(i,j)+4.0D0*k2(i,j)+k3(i,j))
						end if
					else
						f(i,j)=f_old(i,j)+(dt/6.0d0)*(k1(i,j)+4.0D0*k2(i,j)+k3(i,j))
					end if
					if (j==Np) then
						f(i,j)=0.0d0
					end if
					if (f(i,j)<0.0d0) then
						f(i,j)=0.0d0
					end if
				end do
			end do
        Npart_old=Npart
        f_old=f
   end do

   close(18)
   close(19)
   close(16)   
   close(23)
   close(24)
   call cpu_time(finish)

   print* 
   call elapsedtime(int(finish-start))

!------------------------------------------------------------------------------
!                                Graphs files

   call execute_command_line("touch graphs")
   open(15,file='graphs')
   write(15,*) "set terminal png"
   write(15,*) "set nokey" 
   write(15,*) "set output 'Data/Graphs/Npart.png'"
   write(15,*) "set xlabel 'Time (s)'"
   write(15,*) "set ylabel 'Number of particles'"
   write(15,*) "plot 'Data/numpart.dat' with lines lt rgb 'blue'"
   close(15)
   call execute_command_line("gnuplot graphs")
   call execute_command_line("rm graphs")

   call execute_command_line("touch graphs")
   open(15,file='graphs')
   write(15,*) "set terminal png"
   write(15,*) "set nokey" 
   write(15,*) "set output 'Data/Graphs/tenergy.png'"
   write(15,*) "set xlabel 'Time (s)'"
   write(15,*) "set ylabel 'Kinetic'"
   write(15,*) "plot 'Data/tenergy.dat' with lines lt rgb 'blue'"
   close(15)
   call execute_command_line("gnuplot graphs")
   call execute_command_line("rm graphs")
   
   call execute_command_line("touch graphs")
   open(15,file='graphs')
   write(15,*) "set terminal png"
   write(15,*) "set nokey" 
   write(15,*) "set output 'Data/Graphs/venergy.png'"
   write(15,*) "set xlabel 'Time (s)'"
   write(15,*) "set ylabel 'vential Energy'"
   write(15,*) "plot 'Data/venergy.dat' with lines lt rgb 'blue'"
   close(15)
   call execute_command_line("gnuplot graphs")
   call execute_command_line("rm graphs")

   call execute_command_line("touch graphs")
   open(15,file='graphs')
   write(15,*) "set terminal png"
   write(15,*) "set nokey" 
   write(15,*) "set output 'Data/Graphs/initialf.png'"
   write(15,*) "set xlabel 'r'"
   write(15,*) "set ylabel 'p'"
   write(15,*) "splot 'Data/fold.dat' with lines lt rgb 'blue'"
   close(15)
   call execute_command_line("gnuplot graphs")
   call execute_command_line("rm graphs")

   call execute_command_line("touch graphs")
   open(15,file='graphs')
   write(15,*) "set title 'Contour Evolution'"
   write(15,*) "set terminal png nocrop large size 1024,768"
   write(15,*) "set style fill solid 1.00 noborder"
   write(15,*) "set xrange [",r(0),":",r(nr),"]"
   write(15,*) "set mytics"
   write(15,*) "set yrange [",pmin,":",pmax,"]"
   write(15,*) "set xlabel 'r'"
   write(15,*) "set ylabel 'p'"
   write(15,*) "set isosamples 50"
   write(15,*) "set pm3d"
   write(15,*) "set cbrange [0.0:]"
   write(15,*) "unset surface"
   write(15,*) "set view map"
   write(15,*) "set palette rgbformulae 34,35,36"
   write(15,*) "set key outside"
   write(15,*) "nc","=", Nt/stept
   write(15,*) "do for [i=1:nc]{"
   write(15,*) "outfile = sprintf('Data/Graphs/GraphsEvolution/evolution%05.0f.png',i)"
   write(15,*) "set output outfile"
   write(15,*) "splot 'Data/vlasov.dat' using 1:2:3 index (i-1) title sprintf('t=%i',i*",stept,"*0)"
   write(15,*) "}"
   close(15)
   call execute_command_line("gnuplot graphs")
   call execute_command_line("rm graphs")

   call execute_command_line("touch graphs")
   open(15,file='graphs')
   write(15,*) "set title 'Density'"
   write(15,*) "set style data lines"
   write(15,*) "set terminal png nocrop large size 1024,768"
   write(15,*) "stats 'Data/density.dat' nooutput"
   write(15,*) "fmax=STATS_max_y"
   write(15,*) "set xrange [",r(3),":",r(nr),"]"
   write(15,*) "set yrange [",0,":fmax]"
   write(15,*) "set xlabel 'r'"
   write(15,*) "set ylabel 'f(r)'"
   write(15,*) "nc=", Nt/stept
   write(15,*) "do for [i=1:nc]{"
   write(15,*) "outfile = sprintf('Data/Graphs/GraphsDensity/density%05.0f.png',i)"
   write(15,*) "set output outfile"
   write(15,*) "plot 'Data/density.dat' index (i-1) lt 1 lc 1"
   write(15,*) "}"
   close(15)
   call execute_command_line("gnuplot graphs")
   call execute_command_line("rm graphs")
   
   
!------------------------------------------------------------------------------

   call execute_command_line("rm integrals.mod")
   call execute_command_line("rm flux.mod")
   call execute_command_line("rm Data/drdt.dat Data/dpdt.dat Data/p0.dat")

   print*
   print*
   print*, '**********************************'
   print*, '****** Program has finished ******'
   print*, '**********************************'
   print*
   print*, '   Data obtain is in "Data" carpet'
   print*
   print*, '            !See you!'
   print*
   print*
end program




!------------------------------------------------------------------------------
!                      Initial Data for the Density Function
!------------------------------------------------------------------------------ 
real(8) function gaussian(r,p,rc,pc,sr,sp)
   implicit none
    
   real(8):: r,p
   real(8):: rc,pc
   real(8):: sr,sp

   gaussian=(exp((-(r-rc)**2)/(2.0d0*sr**2)+(-(p-pc)**2)/(2.0d0*sp**2)) &
            +exp((-(r+rc)**2)/(2.0d0*sr**2)+(-(p+pc)**2)/(2.0d0*sp**2)))
end function 


real(8) function compact(r,p,rc,pc,sr,sp)
   implicit none
    
   real(8):: r,p
   real(8):: rc,pc
   real(8):: sr,sp
   real(8):: pi

   pi=acos(-1.0d0)
   
   if (abs(r-rc)<=sr*0.5d0 .and. abs(p-pc)<=sp*0.5d0) then
      compact=(1.0d0+cos((pi/sr)*(r-rc)))*(1+cos((pi/sp)*(p-pc)))
   else
      compact=0.0d0
   end if
end function




!-----------------------------------------------------------------------------
!                               Limiter Functions
!-----------------------------------------------------------------------------

subroutine limiter(slope1,slope2,slopecorrect,lim)
   implicit none

   real(8):: slope1                      ! f_m - f_{m-1}
   real(8):: slope2                      ! f_{m+1} - f_m
   real(8):: ss,slopecorrect

   character(50):: lim

! If the slopes have different signs then we set the slope
! correction to zero and return. Otherwise we remember the
! common sign.

   ss=0.0d0

   if ((slope1*slope2)<=0.0d0) then
      slopecorrect = 0.0d0
   else if (slope1>0.0d0) then
      ss = +1.0d0
   else
      ss = -1.0d0
   end if

   if (lim=='minmod') then
      slopecorrect = ss*min(abs(slope1),abs(slope2))

   else if (lim=='mclimiter') then
      slopecorrect = ss*min(abs(2.0d0*slope1),abs(2.0d0*slope2),abs(0.5d0*(slope1+slope2)))

   else if (lim=='superbee') then
      slopecorrect = max(0.0d0,min(2.0d0*slope1,slope2),min(slope1,2.0d0*slope2))
   
   else
      print*
      print*, 'Unknown slope limiter.'
      print*, 'Aborting ...'
      print*
      stop
      
   end if
return
end subroutine




!-----------------------------------------------------------------------------
!           Remove Empty Spaces (Used for the different data files)
!-----------------------------------------------------------------------------
subroutine removespaces(string)

   character(len=*)::      string
   
   integer::               stringlen 
   integer::               last, actual

   stringlen = len (string)
   last = 1
   actual = 1

   do while (actual < stringlen)
      if (string(last:last) == ' ') then
         actual = actual + 1
         string(last:last) = string(actual:actual)
         string(actual:actual) = ' '
      else
         last = last + 1
         if (actual < last) &
            actual = last
      end if
   end do

end subroutine




!----------------------------------------------------------------------------
!                         Calculate elapsed time
!----------------------------------------------------------------------------
!Got from https://rosettacode.org/wiki/Convert_seconds_to_compound_duration#Fortran
subroutine elapsedtime(t)
   
   implicit none

   integer::         t               !Time in seconds
   integer::         ntypes          !Types of time
   parameter (ntypes=5)
   integer::         usize(ntypes)     !Size of each time unit
   character*3::     uname(ntypes)     !Name of each time unit
   parameter (usize=(/7*24*60*60,24*60*60,60*60,60,1/))
   parameter (uname=(/"w","d","h","m","s"/))
   character(28)::   text
   integer::         i,l,n,s

   s=t
   l=0
   
   do i=1,ntypes
      n=s/usize(i)
      if (n.gt.0) then
         s=s-n*usize(i)
         if (l.gt.0) then
            l=l+2
            text(l-1:l)=","
         end if
         write(text(L+1:),1) n,uname(i)
1        format (i0,1x,a)
         l=len_trim(text)
      end if
   end do
   write(6,*) "Time elapsed: ",text(1:l)

end subroutine
