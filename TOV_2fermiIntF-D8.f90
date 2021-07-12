!
!*********************************************************	   
! Module - one propagation step by the Runge-Kutta method
!
      module RungeKutta

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
	   IMPLICIT DOUBLE PRECISION (A-H,O-Z)	
	   real(8):: Yi(Nv),Yi1(Nv)
	   real(8):: k1(Nv),k2(Nv),k3(Nv),k4(Nv),FH(Nv)
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
! for the solution	of ToV equation for 2 steps Fermi star with interaction
! with EoS of Fermi gas + DM


        program MainRK
	    use RungeKutta	
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)		   
 	    integer, parameter :: N1=10000  ! -  number of points for EoS
	    integer, parameter :: Nv=2

 	    real(8) :: Yi(Nv),Yi1(Nv)
		real(8) :: En1(N1),P1(N1),PF1(N1)
		real(8) :: En2(N1),P2(N1),PF2(N1)
		common/EoS_wInt/  En2,P2  ! with int
		common/EoS_woInt/ En1,P1  ! witout int
		
	    Common/inputDHO/ G,PI
		Common/energydens/ ener1

		Common/RadStep2/Rstep2,IEoS
		Common/Param_EOS_woInt/dMwoInt
		Common/Param_EOS_wInt/dMwInt,YY

	    external FunDHO,comp
		
	
!--->
! Output files:
	   open(unit=110,file='1EoSwoInt.dat',status='unknown')
	   open(unit=120,file='1EoSwInt.dat',status='unknown')   
	   
	   open(unit=10,file='1outP.dat',status='unknown')
	   open(unit=19,file='1outM.dat',status='unknown')	
	   open(unit=11,file='1outP0_R.dat',status='unknown')	
	   open(unit=12,file='1outM_R.dat',status='unknown')
	   open(unit=16,file='1outCr.dat',status='unknown')	 
	   open(unit=17,file='1outC_R.dat',status='unknown')	
!-------------------------
! Initial conditions	
        PI=4.*Atan(1.)	 
        dm=0.938
        hqc=0.197	
		G=1.
!
! Input flag to define EoS : 
!       IEoS = 1  ! 1-step EoS wo interaction
!		IEoS = 10 ! 1-step EoS with interaction
!		
!		IEoS = 2  ! 2- step EoS: 1st step - EoS with interaction;  2nd step - EoS wo interaction
!		IEoS = 20 ! 2- step EoS: 1st step - EoS wo interaction;    2nd step - EoS with interaction	
! select IEoS:
        IEoS=2
	 	
		
! 1) Input parameters for EoS wo interaction: ---
! Set the mass of the fermions
        dMwoInt=dm
! 2) Input parameters for EoS with interaction: ---
! Set the mass of the fermions
        dMwInt=dm    ! Ex1
!		dMwInt=dm/2. ! Ex2
! add Interaction term: eq. (36) astro-ph/0605724
! YY - interaction scale - defined in main!		
        YY=10.  ! Ex1 - strenght of interaction		
!		YY=100. ! Ex2 - strenght of interaction	
!		
		Rstep2=2.   ! Ex1 radius for starting "step 2"  - second EoS in r' 
!		Rstep2=10.  ! Ex2 radius for starting "step 2"  - second EoS in r' 		
!-------		
	   rmin=0.001 ! minimal r
! rmax:=50 for IEoS=1;  for IEoS=2: 100 for Y=10.	   
	   rmax=100.   !  Ex1  ! max R considered for investigation in units r'
!	   rmax=1000.  ! Ex2  ! max R considered for investigation in units r'
	   
	   N=10000  ! number of points in r
	   ht=(rmax-rmin)/real(N)

!--------------- EoS: 
   
       Call EoS_Fermi_woInt(N1,PF1,En1,P1)	  !  EoS wo inter
       Call EoS_Fermi_wInt(N1,PF2,En2,P2)     !   EoS w int
!	   
	   do i1=1,N1
	   write(110,*) En1(i1),P1(i1),PF1(i1)  ! Eos1: dimentionless e', P' and Pf'
	   write(120,*) En2(i1),P2(i1),PF2(i1)  ! EoS2
	   end do
	   write(*,*) 'EoS is stored in a grid'
	   
!--------------ToV:---------------------------------------------
        P0min=P1(3)   ! minimal P0 - take from EoS with interaction 
		P0max=100.    ! maximal P0
		
        NPmax=100
        hP0=(P0max-P0min)/Npmax

	    hP0=log((P0max/P0min))/real(Npmax)   ! log step in P0'
		write(*,*) 'Pmin,Pmax:',P1(1),P1(N1),hp0
	   
	   	   do ip0=1,NPmax ! --- loop over P0
! 	   P0=P0min+hP0*real(ip0-1)	  ! linear step
	   P0=P0min*exp(hP0*real(ip0-1))	  ! log step

! Set initial conditions for RK for given P0
	   Yi(2)=P0    ! in P' 
	     AMan=1E-8  ! - small value 
       Yi(1)=AMan

! r  loop at given P0:	 
 	    do it=1,N   ! ---- loop over r for given P0 - N steps
		ri=rmin+ht*real(it)

!		Write(80,*) '1:',ri,Yi
        call RKstep(ri,ht,Yi,Yi1,Nv,FunDHO) ! solution for mass and pressure
!		Write(80,*) '2:',ri,Yi1

! Write output in file '1outRKTask.dat':	

!       write(*,555) ri,Yi1
		write(10,555) ri,Yi1(2)
		write(19,555) ri,Yi1(1) !,Mass(ri)   ! Mass: r, Mnumeric, Manalyric
		
! Compute Compactness:	
    	 Compri=comp(Nv,Yi,ri)
 		 write(16,555) ri,compri	
		 
		if(Yi1(2).lt.0.) GO TO 55  ! P(ri)=0 ==> ri=Rmax for given P0 --->>>>>> go out of the loop!

		 Yi(:)=Yi1(:) ! reset Y_i <-- Y_{i+1} and go to the text rime step		

	 
		  end do  ! loop over ri
		 
55      continue		! <<<<<----  

        rmax=ri  ! maximum radius for given P0
		dMax=Yi1(1)  ! mass at rmax
		
        write(*,*) ip0, ri,Yi1		
		
		if(dMax.gt.Aman) then
		write(11,555) rmax,P0     ! Pressure
		write(12,555) rmax,Yi1(1) ! mass
		write(17,555) rmax,Compri       ! mass
		end if
 		  
		 write(10,777)	
		 write(19,777)
		 write(16,777)					 
		
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
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)		   
	   integer Nv 
       real(8) :: ri,Yi(Nv)
	   real(8) :: FH(Nv)

	   Common/inputDHO/ G,PI
       Common/energydens/ener1	
	   Common/RadStep2/Rstep2,IEoS
	   
	   
	   FH=0.
	   if(Yi(2).lt.0) return
       AMm=Yi(1)
	   
! Obtain Energi=E' for given P', needed for RK	

! Now select scenario for EoS:
       If(IEoS.eq.1)      then         ! 1-step EoS wo interaction
 	     call P_En_woInt(Yi(2),Eneri)    	
		 
       else if(IEoS.eq.10) then        ! 1-step EoS with interaction
 	     call P_En_wInt(Yi(2),Eneri)  
		 
       else if(IEoS.eq.2)  then        ! 2-step EoS: 1st step - EoS with interaction;  2nd step - EoS wo interaction
	     if(ri.le.Rstep2)  call P_En_wInt(Yi(2),Eneri)
	     if(ri.gt.Rstep2)  call P_En_woInt(Yi(2),Eneri)		 
		 
       else if(IEoS.eq.20) then        ! 2-step EoS: 1st step - EoS wo interaction;    2nd step - EoS with interaction	
	   	 if(ri.le.Rstep2)  call P_En_woInt(Yi(2),Eneri)
	     if(ri.gt.Rstep2)  call P_En_wInt(Yi(2),Eneri)	
       end if		 

	   AA1=(1.+4.*PI*ri**3*Yi(2)/AMm) 
 
       FH(1)= 4.*PI*ri**2*Eneri ! f1   mass
	   FH(2)= -AMm*Eneri/ri**2 *(1.+Yi(2)/Eneri)*AA1/(1.-2.*AMm/ri) ! f2 - pressure
	   
       return
       end	  

!*************************************************************************

!*******************************************************
       subroutine EoS_Fermi_woInt(N,Pf,En,P)
	   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
! EoS without interaction	   
! calculate   e' vs P' connected by pf (Fermi momentum in a grid with N steps (defined in main)
	   real(8) :: Rho(N),En(N),P(N),Pf(N) 
	   Common/Param_EOS_woInt/dM   ! mass of ferion - input parameter for EoS1

	   
	   PI=4.*Atan(1.)	
	   hqc = 0.197   ! unit conversion h_bar*c  [GeV*fm]

! limits in Femi momenta for generation of the grid P'(E')
	    Pfmin= 1.E-8 !  wo units
   	    Pfmax= 5000. !
		
!	    dpf = (Pfmax-Pfmin)/real(N)  ! - linear step
	    dpf = Dlog(Pfmax/Pfmin)/real(N)  ! 	- for log scale	
!		write(*,*) 'P-grid:',N,pmin,pmax,dpf
		
	     do i=1,N  ! ---> loop over Pf
 		 Pf(i) = Pfmin*Dexp(real(i-1)*dpf)  ! - Pfermi - log step

! Calculate the Integral for energy density e' and pressure p':
  	      EE = 0.
          PP = 0.
		  Nint=1000
		  dp=Pf(i)/real(Nint)
            do ip=1,Nint
         pmom = real(ip)*dp  ! momentum GeV/c
         E1 = sqrt(pmom**2 + dM**2)	 ! energy GeV
         EE = EE + dp*E1*pmom**2   ! GeV^4  => GeV/hc**3=GeV/fm^3
		 PP = PP + dp*pmom**4/E1 
           end do		 
		 
		 En(i) = EE/PI**2/dM**4         ! dimentionless e'
		 P(i) = PP/(3.*PI**2)/dM**4     ! dimentionless P'

		 end do   ! -- loop ovre Pf

!555    format(10(2XE12.5))
	
       end subroutine
!*******************************************************
       subroutine EoS_Fermi_wInt(N,Pf,En,P)
	   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
! EoS without interaction	   
! calculate   e' vs P' connected by pf (Fermi momentum in a grid with N steps (defined in main)
	   real(8) :: Rho(N),En(N),P(N),Pf(N) 
	   Common/Param_EOS_wInt/dM,YY   ! mass of ferion - input parameter for EoS1

	   
	   PI=4.*Atan(1.)	
	   hqc = 0.197   ! unit conversion h_bar*c  [GeV*fm]

! limits in Femi momenta for generation of the grid P'(E')
!	    Pfmin= 1.E-8 !  Ex1 -wo units
		Pfmin= 1.E-16 ! Ex2 - wo units
   	    Pfmax= 5000. !
		
!	    dpf = (Pfmax-Pfmin)/real(N)  ! - linear step
	    dpf = log(Pfmax/Pfmin)/real(N)  ! 	- for log scale	
!		write(*,*) 'P-grid:',N,pmin,pmax,dpf		

	     do i=1,N  ! ---> loop over Pf
 		 Pf(i) = Pfmin*exp(real(i-1)*dpf)  ! - Pfermi - log step

! Calculate the Integral for energy density e' and pressure p':
  	      EE = 0.
          PP = 0.
		  Nint=1000
		  dp=Pf(i)/real(Nint)
            do ip=1,Nint
         pmom = real(ip)*dp  ! momentum GeV/c
         E1 = sqrt(pmom**2 + dM**2)	 ! energy GeV
         EE = EE + dp*E1*pmom**2   ! GeV^4  => GeV/hc**3=GeV/fm^3
		 PP = PP + dp*pmom**4/E1 
           end do		 
		 
		 EnWOI = EE/PI**2/dM**4         ! dimentionless e'
		 PWOI = PP/(3.*PI**2)/dM**4     ! dimentionless P'
		 
! add Interaction term: eq. (36) astro-ph/0605724
! YY - interaction scale - defined in main!
		 
         ZZ=Pf(i)/dM  ! dimentionless Fermi momentum
		 ENint=(1./3./Pi**2)**2*YY**2*ZZ**6
		 Pint=ENint
		 
         En(i)=EnWOI+ENint ! dimentionless e' with int.
		 P(i)=PWOI+Pint    ! dimentionless P' with int.
		 
!		 write(66,*) i,pf(i),En(i),P(i),ENint
		 end do   ! -- loop ovre Pf

555    format(10(2XE12.5))
	
       end subroutine	   


!------------------------------------------------
        subroutine P_En_woInt(Pc,Ec)   ! for EoS wo int
! Linear interpolation between grid points		
! find E' for given P'
! Output: Ec for given Pc
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 	    integer, parameter :: N=10000
	    real(8) :: Ein(N),Pin(N)
		common/EoS_woInt/ Ein,Pin  ! wo int

         Ec=0.
	   	 do i=1,N-1
        if((Pc.lt.Pin(1)) .or. (Pc.gt.Pin(N))) then
		write(*,*) 'wo:P-E grid is too small',Pc,Pin(1),Pin(N)
		stop
		end if
		 if( (Pc.gt.Pin(i)) .and.(Pc.le.Pin(i+1))) then ! find position of P' in the grid
		 Ec = (Ein(i)+Ein(i+1))/2.        ! interpolation E'
		 go to 321
		 end if
		end do 
 321     continue		
!!		write(433,*) Pc,Ec,i

	   return
	   end subroutine	
!------------------------------------------------
        subroutine P_En_wInt(Pc,Ec)  ! for EoS with int
! Linear interpolation between grid points		
! find E' for given P'
! Output: Ec for given Pc
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 	    integer, parameter :: N=10000
	    real(8) :: Ein(N),Pin(N) 
		common/EoS_wInt/ Ein,Pin ! with int

         Ec=0.
	   	 do i=1,N-1
        if((Pc.lt.Pin(1)) .or. (Pc.gt.Pin(N))) then
		write(*,*) 'w: P-E grid is too small',Pc,Pin(1),Pin(N)
		stop
		end if
		 if( (Pc.gt.Pin(i)) .and.(Pc.le.Pin(i+1))) then ! find position of P' in the grid
		 Ec = (Ein(i)+Ein(i+1))/2.        ! interpolation E'
		 go to 321
		 end if
		end do 
 321     continue		
!!		write(433,*) Pc,Ec,i

	   return
	   end subroutine
	   
!*************************************************
       DOUBLE PRECISION function comp(Nv,Yi,ri)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)	
	   real(8):: Yi(Nv) 
	   Common/inputDHO/ G,PI
!	   G=1.
	   comp= G*Yi(1)/ri
	   return
	   end function
	   