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
			
       Yi1(:)=Yi(:)+ht/6.*(k1(:)+2.*k2(:) +2.*k3(:)+k4(:))

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
 	    integer, parameter :: N1=50000  ! -  number of points for EoS
	    integer, parameter :: Nv=2
		integer, parameter :: Nrm=500
		real(8) :: rmaxV(Nrm),dMaxR(Nrm)

 	    real(8) :: Yi(Nv),Yi1(Nv)
		real(8) :: En1(N1),P1(N1),PF1(N1)
		real(8) :: En2(N1),P2(N1),PF2(N1)
		real(8) :: En3(N1),P3(N1),PF3(N1)
		
		common /EoS_wInt/  En2,P2  ! with int
		common /EoS_woInt/ En1,P1  ! without int
		common /EoS_Bag/   En3,P3, NBag  ! with int - MIT Bag model
		
	    Common/inputDHO/ G,PI
		Common/RadStep2/Rstep2,IEoS
		Common/Param_EOS_woInt/dMwoInt
		Common/Param_EOS_wInt/dMwInt,YY
        Common/Param_EOS_Bag/dMBag
	    external FunDHO,comp,EnergiIncomp
		
	
!--->
! Output files:
	   open(unit=110,file='1EoSwoInt.dat',status='unknown')
	   open(unit=120,file='1EoSwInt.dat',status='unknown')
	   open(unit=121,file='1EoSBag.dat',status='unknown')
	   
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
!		IEoS = 21 ! 2-step EoS: 1st step -  EoS with interaction - Quark core: Bag MIT for quarks;    
!                               2nd step -  EoS with interactionfermions - nucleons
!		IEoS = 22 ! 2-step EoS: 1st step -  EoS with interactionfermions - DM fermions
!                               2nd step -  EoS with interaction -  Bag MIT for quarks                           
!
!       IEoS = 30 ! 2- step EoS: 1st step - EoS incompressible fluid;    2nd step - EoS with interaction
! select IEoS:
        IEoS=21
		
! 1) Input parameters for EoS wo interaction: ---
! Set the mass of the fermions
        dMwoInt=dm
! 2) Input parameters for EoS with interaction: ---
! Set the mass of the fermions
        dMwInt=dm    ! Ex1
!		dMwInt=dm/2. ! Ex2
        dMwInt=1E-6 ! dm/1.E6  ! Ex3 axino star warm dark matter mass
!        dMwInt=dm*100   ! Ex4   neutralino star : cold dark matter
		
! 3) Input parameters for EoS with MIT Bag interaction for quark stars: ---
! Set the mass of the quarks
        dMBag=5.E-3  ! mass of the bare u,d quark =5-10 MeV (all must be in units GeV!)
		
! add Interaction term: eq. (36) astro-ph/0605724
! YY - interaction scale - defined in main!		
        YY=10.  ! Ex1 - strenght of interaction	for fermions	
!		YY=100. ! Ex2 - strenght of interaction	
!       YY=1. ! Ex3 - strenght of interaction
!		
!		Rstep2=2.   ! Ex1 radius for starting "step 2"  - second EoS in r'
!		Rstep2=10.  ! Ex2 radius for starting "step 2"  - second EoS in r'
!        Rstep2=(1.E+13)/2.
        Rstep2= 10.  !7.		! radius for change of EoS
!-------		
	   rmin=0.01 !0.001 ! minimal r
! rmax:=50 for IEoS=1;  for IEoS=2: 100 for Y=10.	
!	   rmax=100.   !  Ex1  ! max R considered for investigation in units r'
	   rmax=50.  ! Ex2  ! max R considered for investigation in units r'
	
	   N=1000 !2000, 500  ! number of points in r   !!! 1000
	   ht=(rmax-rmin)/real(N)

!--------------- EoS:

       Call EoS_Fermi_woInt(N1,PF1,En1,P1)	  !  EoS wo interaction
       Call EoS_Fermi_wInt(N1,PF2,En2,P2)     !   EoS w int
	   Call EoS_Fermi_Bag(N1,PF3,En3,P3,NBag)  
!	
	   do i1=1,N1
	   write(110,*) En1(i1),P1(i1),PF1(i1)  ! Eos1: wo int  dimentionless e', P' and Pf'
	   write(120,*) En2(i1),P2(i1),PF2(i1)  ! EoS2: with int
!	   write(121,*) En3(i1),P3(i1),PF3(i1)  ! EoS2: with int
	   end do
	   	   do i1=1,N1 !NBag
	   write(121,*) En3(i1),P3(i1),PF3(i1)  ! EoS2: with int with Bag MIT quarks - Nbag points!
	   end do
	   write(*,*) 'EoS1,2 is stored on a grid',EN1(1),EN1(N1),P1(1),P1(N1)
	   write(*,*) 'EoS-Bag is stored on a grid',Nbag,EN3(1),EN3(N1)
     &            ,P3(1),P3(N1)

!        stop
!---------------------------------------------------------------	
!--------------ToV:---------------------------------------------
        write(*,*) '--> go to RK with iEoS=',IEoS
!		
        P0min= P1(3)   ! minimal P0 - take from EoS without interaction
		P0min= P3(3)   ! minimal P0 - take from EoS without interaction MIT Bag
!       P0min= P2(2)   ! case EX4 
        P0min=1E-10
		P0max= 50.     !10. !0.01   ! maximal P0
!        P0max=P2(N1)
!		P0max=P3(NBag-1) !*1E-12
        NPmax=1000  !1000      ! number of steps for TOV ! 1000 ! 5000
        hP0=(P0max-P0min)/Npmax    ! grid step in P0

	    hP0=log((P0max/P0min))/real(Npmax)   ! log step in P0'
		write(*,*) 'P0min,P0max,h:',P0min,P0max,hp0,Nbag
		
		irm=0
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
	   	   do ip0=1,NPmax ! --- loop over P0 (initial pressure)
! 	   P0=P0min+hP0*real(ip0-1)	  ! linear step
	   P0=P0min*exp(hP0*real(ip0-1))	  ! log step

! Set initial conditions for RK for given P0
	   Yi(2)=P0    ! in P'
	     AMan=1E-8  ! - small value
!		AMan=4.*Pi*rmin**2*En3(1)  ! - small value
!		write(*,*) 'Mmin=',Aman
       Yi(1)=AMan
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! r  loop at given P0:	

 	    do it=1,N   ! ---- loop over r for given P0 - N steps
		ri=rmin+ht*real(it)

!		Write(80,*) '1:',ri,Yi
        call RKstep(ri,ht,Yi,Yi1,Nv,FunDHO) ! solution for mass and pressure by Runge-Kutta
!		Write(80,*) '2:',ri,Yi1

! Write output in file '1outRKTask.dat':	

!        write(*,555) ri,Yi1       !----------------- pressure screen
		write(10,555) ri,Yi1(2)    !----------------- file "1outM"
		write(19,555) ri,Yi1(1)    !----------------- file "1outP" !,Mass(ri)   ! Mass: r, Mnumeric, Manalyric
	   
! Compute Compactness:	
    	 Compri=comp(Nv,Yi1,ri)  
 		 write(16,555) ri,compri	 !----------------- file "1outC_r"
		
! Special output:	E(r)
 !        go to 981  ! to jump over special output
		 Pc=Yi1(2)
		 Eneric=0.
       if(IEoS.eq.30 .and. Pc.gt.0) then        ! 2-step EoS: 1st step - EoS incompressible fluid;    2nd step - EoS with interaction
      	 if(ri.le.Rstep2)   Eneric=EnergiIncomp(ri) 
	     if(ri.gt.Rstep2)  call P_En_wInt(Pc,Eneric) 
	   end if
       if(IEoS.eq.21 .and. Pc.gt.0) then 
	   	 if(ri.le.Rstep2)  call P_En_Bag(Pc,Eneric)
	     if(ri.gt.Rstep2)  call P_En_wInt(Pc,Eneric) 
		end if 
		 write(48,555) ri, Eneric,Pc	
981     continue		
!		 
!---pressure becomes negatine -->> maximum radius!		 
		if(Yi1(2).lt.0.) GO TO 55  ! P(ri)=0 ==> ri=Rmax for given P0 --->>>>>> go out of the loop!

		 Yi(:)=Yi1(:) ! reset Y_i <-- Y_{i+1} and go to the text rime step	
	
		  end do  ! loop over ri
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx		
55      continue		! <<<<<----
	
        rmax=ri  ! maximum radius for given P0
		dMax=Yi1(1)  ! mass at rmax
		
		irm=irm+1
		rmaxV(irm)=rmax
		dMaxR(irm)=dMax
		
 !       write(*,*) ip0, ri,Yi1		
		
		if(dMax.gt.Aman) then !--->>
		
	     IF(ri.gt.Rstep2) then !--- EOS=30! select ONLY rmax > Rc!!!		
		write(11,555) rmax,P0     ! Pressure
		write(12,555) rmax,Yi1(1) ! mass
		write(17,555) rmax,Compri       ! mass
		 END IF    !--- end if r>Rc		
		 
! Calculate derivatives of mass_max over Rmax  
		 Deltar=(rmaxV(irm)-rmaxV(irm-1))
		 dmrm=(dMaxR(irm)-dMaxR(irm-1))/Deltar
		 rmaxR=rmaxV(irm)
!		 write(18,555) rmaxR,dmrm  !,dMaxR(irm),dMaxR(irm-1),deltar
! Polytrops:
        Gam=4./3.
!		Gam=5./3.
!		Gam=2.
        ConstPoly=dMax**(2.-Gam)*rmax**(3.*Gam-4.)		
!		write(20,555) rmaxR,ConstPoly
		end if ! <<--
		
		 write(10,777)	
		 write(19,777)
		 write(16,777)	
		 write(48,777)	 
		
		end do   ! loop over P0
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
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
	   Common/RadStep2/Rstep2,IEoS
	   External EnergiIncomp,P_En_Bag	
	
	   FH=0.
	   if(Yi(2).lt.0) return             ! P<0
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
		 
       else if(IEoS.eq.21) then        ! 2-step EoS: 1st step -  EoS with interaction - Quark core: Bag MIT for quarks;    2nd step - EoS with interactionfermions - nucleons
	   	 if(ri.le.Rstep2)  call P_En_Bag(Yi(2),Eneri)
	     if(ri.gt.Rstep2)  call P_En_wInt(Yi(2),Eneri)   
		 
       else if(IEoS.eq.22) then        ! 2-step EoS: 1st step -  EoS with interaction  DM fermions; 2nd step - EoS with interaction- Quark core: Bag MIT for quarks
	     if(ri.le.Rstep2)  call P_En_wInt(Yi(2),Eneri)  	
	   	 if(ri.gt.Rstep2)  call P_En_Bag(Yi(2),Eneri)
!		 write(*,*) '2:',ri,Yi(2),Eneri
		 		 
       else if(IEoS.eq.30) then        ! 2-step EoS: 1st step - EoS incompressible fluid;    2nd step - EoS with interaction
       	 if(ri.le.Rstep2)   Eneri=EnergiIncomp(ri) 
	     if(ri.gt.Rstep2)  call P_En_wInt(Yi(2),Eneri) 
       end if
		
	   AA1=(1.+4.*PI*ri**3*Yi(2)/AMm)

       FH(1)= 4.*PI*ri**2*Eneri ! f1   mass
	   FH(2)= -AMm*Eneri/ri**2 *(1.+Yi(2)/Eneri)*AA1/(1.-2.*AMm/ri) ! f2 - pressure

!	   write(80,555) ri,Eneri,Yi
555     format(10(2XE12.5))	   
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
!	    Pfmin= 1.E-16 !  wo units  Ex1,2
		Pfmin= 1.E-10 ! Ex3 - wo units
   	    Pfmax= 1.E+4                      !10. !5000. !
		
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

		 end do   ! -- loop over Pf

!555    format(10(2XE12.5))
	
       end subroutine
!*******************************************************
       subroutine EoS_Fermi_wInt(N,Pf,En,P)
	   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
! EoS without interaction	
! calculate   e' vs P' connected by pf (Fermi momentum in a grid with N steps (defined in main)
	   real(8) :: Rho(N),En(N),P(N),Pf(N)
	   Common/Param_EOS_wInt/dM,YY   ! mass of ferion - input parameter for EoS1

	
	    write(*,*) 'dmf=',dm
	   PI=4.*Atan(1.)	
	   hqc = 0.197   ! unit conversion h_bar*c  [GeV*fm]

! limits in Femi momenta for generation of the grid P'(E')
!	    Pfmin= 1.E-8 !  Ex1 -wo units
!		Pfmin= 1.E-16 ! Ex2 - wo units
!		Pfmin= 1.E-30 ! Ex3 - wo units
        Pfmin= 1.E-10 ! Ex3 - wo units
   	    Pfmax= 100. !1.E+4                          !10. !5000. !
		
!	    dpf = (Pfmax-Pfmin)/real(N)  ! - linear step
	    dpf = log(Pfmax/Pfmin)/real(N)  ! 	- for log scale	
!		write(*,*) 'P-grid:',N,pmin,pmax,dpf		

	     do i=1,N  ! ---> loop over Pf
 		 Pf(i) = Pfmin*exp(real(i-1)*dpf)  ! - Pfermi - log step

! Calculate the Integral for energy density e' and pressure p':
  	      EE = 0.
          PP = 0.
		  Nint=10000
		  dp=Pf(i)/real(Nint)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            do ip=1,Nint
         pmom = real(ip)*dp  ! momentum GeV/c
         E1 = sqrt(pmom**2 + dM**2)	 ! energy GeV
         EE = EE + dp*E1*pmom**2   ! GeV^4  => GeV/hc**3=GeV/fm^3
		 PP = PP + dp*pmom**4/E1
           end do		
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++		
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
!		 write(66,*) i,pf(i),En(i),P(i),ENint
		 end do   ! -- loop over Pf

555    format(10(2XE12.5))
	
       end subroutine	

!------------------------------------------------
        subroutine P_En_woInt(Pc,Ec)   ! for EoS wo int
! Linear interpolation between grid points		
! find E' for given P'
! Output: Ec for given Pc
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 	    integer, parameter :: N=50000
	    real(8) :: Ein(N),Pin(N)
		common/EoS_woInt/ Ein,Pin  ! wo int

         Ec=0.
	   	 do i=1,N-1
        if((Pc.lt.Pin(1)) .or. (Pc.gt.Pin(N))) then
		write(*,*) 'wo:P-E grid is too small, wo int',Pc,Pin(1),Pin(N)
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
 	    integer, parameter :: N=50000
	    real(8) :: Ein(N),Pin(N)
		common/EoS_wInt/ Ein,Pin ! with int

         Ec=0.		 
c         WRITE(*,*)'PC,Pin(1),Pin(N)', PC,Pin(1),Pin(N)
	   	 do i=1,N-1
        if((Pc.lt.Pin(1)) .or. (Pc.gt.Pin(N))) then
		write(*,*) 'w: P-E grid is too small , with int',Pc,Pin(1),Pin(N)
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
        subroutine P_En_Bag(Pc,Ec)  ! for EoS with int MIT Bag model
! Linear interpolation between grid points		
! find E' for given P'
! Output: Ec for given Pc
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 	    integer, parameter :: N=50000
	    real(8) :: Ein(N),Pin(N)
		common/EoS_Bag/ Ein,Pin, Nbag ! with int

         Ec=0.		 
 !        WRITE(80,*)'PC,Pin(1),Pin(N)', PC,Pin(1),Pin(N),Nbag
	   	 do i=1,N-1 !Nbag-1
!        if((Pc.lt.Pin(1)) .or. (Pc.gt.Pin(Nbag))) then
        if((Pc.lt.Pin(1)) .or. (Pc.gt.Pin(N))) then
		write(*,*) 'Bag MIT: P-E grid is too small',Pc,Pin(1),Pin(Nbag)
		  if (Pc.lt.Pin(1)) then
		     write(*,*) 'Bag MIT: too small : P_center < Pin1 '
		  else if (Pc.gt.Pin(N)) then
		     write(*,*) 'Bag MIT: too big : P_center > Pin(N) '
		  end if 
		stop
		end if
		 if( (Pc.ge.Pin(i)) .and.(Pc.lt.Pin(i+1))) then ! find position of P' in the grid
		 Ec = (Ein(i)+Ein(i+1))/2.        ! interpolation E'
		 go to 321
		 end if
		end do
 321     continue		
		write(433,*) Pc,Ec,i

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
	   
!*******************************************************************************
       DOUBLE PRECISION function EnergiIncomp(r)
! E=const - incompressuble fluid	   
	   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	   
!	   EnergiIncomp=0.1
!	   EnergiIncomp=2.6569E-5
!	  	EnergiIncomp=0.0001 ! -case 1
 	  	EnergiIncomp=0.01   ! case 2
       end function	
	   
!*******************************************************
       subroutine EoS_Fermi_Bag(N,Pf,En,P,Nbag)
	   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
! EoS without interaction for Qurk stars - MIT bag model
! calculate   e' vs P' connected by pf (Fermi momentum in a grid with N steps (defined in main)
	   real(8) :: Rho(N),En(N),P(N),Pf(N)
	   Common/Param_EOS_Bag/dM   ! mass of fermion - input parameter for EoS1

	   PI=4.*Atan(1.)	
	   hqc = 0.197   ! unit conversion h_bar*c  [GeV*fm]
	    write(*,*) 'dmq=',dm
		
         En=0.    
		 P=0.   
		 Pf=0.
		 
! limits in Femi momenta for generation of the grid P'(E')
        Pfmin= 1E-10 !0.5 !1.E-1 ! Ex3 - wo units
   	    Pfmax= 100. !100. !0.8                           !10. !5000. !
!   	    Pfmax= 1.E+4 
		
!	    dpf = (Pfmax-Pfmin)/real(N)  ! - linear step
	    dpf = log(Pfmax/Pfmin)/real(N)  ! 	- for log scale	
!		write(*,*) 'P-grid:',N,pmin,pmax,dpf		

         iBag=0
	     do i=1,N  ! ---> loop over Pf
 		 Pfi = Pfmin*exp(real(i-1)*dpf)  ! - Pfermi - log step

! Calculate the Integral for energy density e' and pressure p':
  	      EE = 0.
          PP = 0.
		  Nint=10000
		  dp=Pfi/real(Nint)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            do ip=1,Nint
         pmom = real(ip)*dp  ! momentum GeV/c
         E1 = sqrt(pmom**2 + dM**2)	 ! energy GeV
         EE = EE + dp*E1*pmom**2   ! GeV^4  => GeV/hc**3=GeV/fm^3
		 PP = PP + dp*pmom**4/E1
           end do		
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++		
		 EnWOI = EE/PI**2        
		 PWOI = PP/(3.*PI**2)
		
! add Interaction term: eq. (46) SB lecture - MIT Bag Model
! YY - interaction scale - defined in main!
		 Bag=0.145**4  ! GeV
		 Evac0=4.*Bag
		 DeltaE=EnWOI-Evac0
		 
!         EnP=EnWOI /Evac0              ! dimentionless e'   --???
		 EnP=DeltaE /Evac0              ! dimentionless e' 
		 Pp=DeltaE/3. /Evac0        ! dimentionless P' with int. - Bag MIT		 
		 
	     if(DeltaE .gt.0.) then
		 iBag=iBag+1
		 
!        En(iBag)=EnP             ! dimentionless e' 
!		 P(iBag)=Pp    ! dimentionless P' with int. - Bag MIT		 
!		 Pf(iBag)=Pfi
		 
         En(i)=EnP             ! dimentionless e' 
		 P(i)=Pp    ! dimentionless P' with int. - Bag MIT		 
		 Pf(i)=Pfi		 
		 end if		
		 NBag=iBag

		 end do   ! -- loop over Pf

555    format(10(2XE12.5))
	
       end subroutine	
	   