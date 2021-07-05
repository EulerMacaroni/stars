!
!*********************************************************	   
! Module - one propagation step by the Runge-Kutta method
!
      module RungeKutta
      implicit none	  
      integer, parameter :: Dk = 4  ! kind=Dk
	  real(kind = Dk) ht

	  contains
!-------------------------------
! Solution of Nv 1st order differential equarions dy_i/dt, i=1,...Nv
! Input: 
! Nv - number of first order equations to solve
! ht - r step
! Yi(Nv) - initial vector (solution at ri=tf-ht) 
! Yi1(Nv) - solution 
! FuntY(Nv,t,Yin,FH) - function FH_i=dy_i/dt, i=1,..Nv (defined by user)

	   subroutine RKstep(tf,ht,Yi,Yi1,Nv,FuntY)
	   integer, intent(in):: Nv
	   real(kind = Dk), intent(in):: tf,ht,Yi(Nv)
	   real(kind = Dk), intent(out):: Yi1(Nv)
	   real (kind = Dk) k1(Nv),k2(Nv),k3(Nv),k4(Nv),FH(Nv),ri
	   external FuntY
	   
! k1:	 
       ri=tf-ht
	   call FuntY(Nv,ri,Yi,FH)  
	   k1(:)=FH(:)
! k2:	   
	     Yi1(:)=Yi(:)+ht*k1(:)/2.
	     call FuntY(Nv,ri+ht/2.,Yi1,FH) 
		 k2(:)=FH(:)
! k3:		 
	   	   Yi1(:)=Yi(:)+ht*k2(:)/2.
		   call FuntY(Nv,ri+ht/2.,Yi1,FH) 
	       k3(:)=FH(:)	
! k4:
	   	   	 Yi1(:)=Yi(:)+ht*k3(:)
			 call FuntY(Nv,ri+ht,Yi1,FH) 
	         k4(:)=FH(:)
			 
       Yi1(:)=Yi(:)+ht/6.*(k1(:)+2.*k2(:) &
          +2.*k3(:)+k4(:))

	   return
	   end subroutine
!---------------------------	
       end module RungeKutta
!--------------------------------------------------------		   
!	   
!*********************************************************	
! Main program 
!  for the solution	of ToV equation for Fermi star with EoS of Fermi gas


       program MainRK
	   use RungeKutta	
 	   integer, parameter :: N1=100000
	   integer, parameter :: Nv=2
	   integer it
	   real(kind = Dk) :: ri,Yi(Nv),Yi1(Nv),Yan,P0
        
		real :: Patmo,Rstar,Mass
		real(kind = Dk) rmin,rmax,G,PI,rho0
		real :: dr,dM,hcq,Pf,Pfe
	    real :: Rho(N1),En(N1),P(N1)
		real :: Rho1(N1),En1(N1),P1(N1)
        real :: Ean(N1),Pan(N1)
	    Common/inputDHO/ rmin,rmax,N,G,PI,rho0
		Common/energydens/ ener1
		common/EoS/ En1,P1,Rho1
	    external FunDHO,comp
		
       PI=4.*Atan(1.)		
!--->
! Output file with the r evolution and final accuracy:	   
	   open(unit=10,file='1outP.dat',status='unknown')
	   open(unit=19,file='1outM.dat',status='unknown')	
	   open(unit=11,file='1outP0_R.dat',status='unknown')	
	   open(unit=12,file='1outM_R.dat',status='unknown')		   
!-------------------------
! Initial conditions	  
        dm=0.938
        hqc=0.197	
		G=1.
		const2= dm**4/hqc**3   ! epsilon0
		
	   rmin=0.001 ! minimal r
	   RRmax=20.  ! max R considered for investigation in units r'
	   N=1000  ! number of points in r
	   ht=(RRmax-rmin)/real(N)

!--------------- EoS: 
	   call fermi(N1,Rho1,En1,P1)    ! numetical
!	   call EoSFanal(N1,Pan,Ean) ! analitic S-B - works only for pf>>m_N
	   do i1=1,N1
	   EpsilR=En1(i1)/rho1(i1)
	   rhoPr=rho1(i1)/Const2
	   write(554,*) En1(i1),P1(i1),Rho1(i1) !,EpsilR  ! , e', P',rho
!	   write(555,*) Ean(i1),Pan(i1) !,EpsilR  ! e', P',e/rho
	   end do

!--------------ToV:---------------------------------------------
        P0min=P1(2)   !1E-8 !P1(2)
		P0max=10. !P1(N1)
		
        NPmax=100
        hP0=(P0max-P0min)/Npmax

	    hP0=Alog((P0max/P0min))/real(Npmax)   ! log step in P0'
		write(*,*) 'Pmin,Pmax:',P1(1),P1(N1),hp0
		Rhomin=Rho1(2) !  1E-16
	   
	   	   do ip0=1,NPmax ! --- loop over P0
! 	   P0=P0min+hP0*real(ip0-1)	  ! linear step
	   P0=P0min*exp(hP0*real(ip0-1))	  ! log step

! Set initial conditions for RK for given P0
	   Yi(2)=P0    ! in P' 
	   call P_En1(P0,E0,Rho0)    ! interpolation program to find corresponding rho0 - needed to set the limits to the EoS fermiRL:
	   write(*,*) 'P0,E0,rho0=',ip0,P0,E0,Rho0
		
!	   Vol=4.*Pi*(rmin+ht)**3/3.
	   AMan=1E-8  ! - small value E0*Vol
       Yi(1)=AMan

! r  loop at given P0:	 
 	    do it=1,N   ! ---- loop over r for given P0 - N steps
		ri=rmin+ht*real(it)

!		Write(80,*) '1:',ri,Yi
        call RKstep(ri,ht,Yi,Yi1,Nv,FunDHO) ! solution for mass and pressure
!		Write(80,*) '2:',ri,Yi1

! Write output in file '1outRKTask.dat':	

		write(10,555) ri,Yi1(2)
		write(19,555) ri,Yi1(1) !,Mass(ri)   ! Mass: r, Mnumeric, Manalyric
		
		if(Yi1(2).lt.0.) GO TO 55  ! P(ri)=0 ==> ri=Rmax for given P0 --->>>>>> go out of the loop!

		 Yi(:)=Yi1(:) ! reset Y_i <-- Y_{i+1} and go to the text rime step		

! Compute C:	
    	 Compri=comp(Nv,Yi,ri)
! 		 write(13,*) ri,compri
		 
		  end do  ! loop over ri
		 
55      continue		! <<<<<----  

        rmax=ri  ! maximum radius for given P0

		write(11,555) rmax,P0 !, Pan
		write(12,555) rmax,Yi1(1) ! mass
 		  
		 write(10,777)	
		 write(19,777)
				 
		
		end do   ! loop over P0
		      
 777    format(/)

555     format(10(2XE12.5))
553     format(10(2X,E12.5))
       stop
       end 

   
!------------------------------------------------------------------------
! Input: ri = ri
! Pi(2) - function P'_i(r)=y'(r)
!
       subroutine FunDHO(Nv,ri,Yi,FH)
	   
	   USE RungeKutta
	   integer Nv 
       real(kind = Dk), intent(in):: ri,Yi(Nv)
	   real(kind = Dk), intent(out):: FH(Nv)

		real(kind = Dk) rmin,rmax,G,PI,rho0,AA1	
		integer N
	   Common/inputDHO/ rmin,rmax,N,G,PI,rho0
       Common/energydens/ener1	   
	   
	   FH=0.
	   if(Yi(2).lt.0) return
       AMm=Yi(1)
	   call P_En(Yi(2),Eneri,Rhoi)  ! Obtain Energi=E' for given P', needed for RK

	   AA1=(1.+4.*PI*ri**3*Yi(2)/AMm) 
 
       FH(1)= 4.*PI*ri**2*Eneri ! f1   mass
	   FH(2)= -AMm*Eneri/ri**2 *(1.+Yi(2)/Eneri)*AA1/(1.-2.*AMm/ri) ! f2 - pressure
	   
       return
       end	  

!*************************************************************************


!*************************************************
       real function comp(Nv,Yi,ri)
	   USE RungeKutta
	   real :: ri
	   real(kind = Dk), intent(in):: Yi(Nv) 
	   real(kind = Dk) rmin,rmax,G,PI,rho0	
	   Common/inputDHO/ rmin,rmax,N,G,PI,rho0	   
!	   G=1.
	   comp= G*Yi(1)/ri
	   return
	   end function
	   
!******************************************************
       subroutine integral(Pf,EE,PP)
! Input:
! N, Pf [GeV/c] - upper limit of integration int_0^Pf dp 
! Output: EE-energy density, PP-pressure in units GeV/fm^3
       integer N
       real :: dp,dM,pmom,E1,Pf
       real,intent(out) :: EE,PP
       integer :: i
	   
       dM = 0.938 ! GeV
	   N=10000   ! define grid for integration, used ONLY inside this subroutine!
       dp = (Pf-0.)/real(N)  ! GeV/c

       EE = 0.
       PP = 0.
       do i=1,N
         pmom = real(i)*dp  ! momentum GeV/c
         E1 = sqrt(pmom**2 + dM**2)	 ! energy GeV
         EE = EE + dp*E1*pmom**2   ! GeV^4  => GeV/hc**3=GeV/fm^3
         PP = PP + dp*pmom**4/(3.*E1)
!		 write(66,*) i,pmom,E1,EE,PP
       end do	
!555    format(10(2XE12.5))	  
        return
       end subroutine	


!*******************************************************
       subroutine fermi(N,Rho,En,P)
! calculate  rho, e' vs P' in a grid with N steps (defined in main)
	   real :: dr,dM,hcq,Pf,Pfe
	   real :: Rho(N),En(N),P(N) 
	   real :: PI,g,EE,PP
	   integer :: i

	   
	   PI=4.*Atan(1.)	
	   hqc = 0.197   ! unit conversion h_bar*c  [GeV*fm]
	   dM = 0.938   ! mass of a neutron in 10**3 MeV
	   g=2.     ! fermions	  
	   
	    PFac=g/2./pi**2
	    Const= 1./hqc**3*PFac
		const2= dm**4/hqc**3   ! epsilon0

! limits in density for generation of the grid P'(E')
	    Rhomin= 1.E-14 ! 1E-8 !0. !fm^-3
   	    Rhomax= 1000. ! 10.
		
!	    dr = (Rhomax-Rhomin)/real(Np1)  ! fm^-3  - linear step
	    dr = Alog(Rhomax/Rhomin)/real(N)  ! fm^-3	- for log scale	
		write(*,*) 'P-grid:',N,rhomin,rhomax	
		
	     do i=1,N
!		 Rho(i) = rhomin+real(i-1)*dr  ! fm^-3
 		 Rho(i) = rhomin*exp(real(i-1)*dr)  ! fm^-3 - log step

		 Pf= (6.*PI**2*Rho(i)/g)**(1./3.)    !fermi-momentum fm^-1
		 Pfe= Pf*hqc   ! converting units ! GeV/c
!	     write(71,*) i,rho(i),pfe,pf,dr
		 
		 call integral(Pfe,EE,PP)
		 
!		 Rho(i)=Rho(i)  ! as input in fm^-3
		 En(i) = EE *Const/const2  ! dimentionless e'
		 P(i) = PP *Const/const2   ! dimentionless P'

!		 write(77,*) i,Rho(i),En(i),P(i),Pf
		 
		 end do 

!555    format(10(2XE12.5))
	
       end subroutine	

!------------------------------------------------
       subroutine EoSFanal(N,Pan,Ean)
! analitic solution for P', E' for Fermi gas EoS - cf. S-B. paper	   
	   real :: Ean(N),Pan(N) 
	   
	   PI=4.*Atan(1.)
       dm=0.938
	   
	   pfmin=1E-10 !dm
	   pfmax=20.
	   
	   hp=(pfmax-pfmin)/real(N)
	   do ip=1,N
	   pf=pfmin+hp*real(ip-1)
	   
	   z=pf/dm
	   Z1=sqrt(1.+z**2)
	   SZ=1./sinh(z)
	   Ean(ip)=((2.*z**3+z)*Z1-SZ)/8./Pi**2
	   Pan(ip)=((2.*z**3-3.*z)*Z1+3.*SZ)/24./Pi**2
	   end do
	   return
	   end subroutine	
	   
!------------------------------------------------
        subroutine P_En(Pc,Ec,Rhoc)
! Linear interpolation between grid points		
! find E' for given P'
! Output: Ec for given Pc
 	    integer, parameter :: N=100000
	    real :: Ein(N),Pin(N),Rhoin(N) 
		common/EoS/ Ein,Pin,Rhoin
	   
!	     write(*,*) N,Pc,Ec
         Ec=0.
	   	 do i=1,N-1
!		 write(*,*) i,Ein(i),Pin(i)

        if((Pc.lt.Pin(1)) .or. (Pc.gt.Pin(N))) then
		write(*,*) '2:P grid is too small',Pc,Pin(1),Pin(N)
		stop
		end if
		 if( (Pc.gt.Pin(i)) .and.(Pc.le.Pin(i+1))) then ! find position of P' in the grid
		 Ec = (Ein(i)+Ein(i+1))/2.        ! interpolation E'
		 Rhoc = (Rhoin(i)+Rhoin(i+1))/2.  ! interpolation rho
		 go to 321
		 end if
		end do 
 321     continue		
!!		write(433,*) Pc,Ec,i

	   return
	   end subroutine	
!------------------------------------------------
	   	   