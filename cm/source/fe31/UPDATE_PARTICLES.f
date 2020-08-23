      SUBROUTINE UPDATE_PARTICLES(NBJ,NEELEM,NORD,NPNE,NVJE,nx,NXI,CE,
     & XP, ERROR, *)       

C#### Subroutine: UPDATE_PARTICLES
C###  Description:
C###    

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M),NORD(5,NE_R_M),
     & NPNE(NNM,NBFM,NEM), NVJE(NNM,NBFM,NJM,NEM),
     & NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 CE(NMM,NEM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,ne,noelem,np1,np2,np3,
     & nv1, nv2,ne0,nx,gen0,gen
      REAL*8 k,T,mu,dp,Cc,l,Patm,u,D,rad,rhop,dpmm,
     & A1,A2,A3,lambda,temp,Db,Da,Dc,Loss,Q,Qe,tempA,
     & Area,Area_i, eta, etaO, etaI,etaS,etaD,Stk,theta,phi_g,
     & Tau, Delta, Gravity, ex,ey,ez,vx,vy,vz, infinity
      LOGICAL SEARCH, 
     !     Functions
      REAL*8 LENGTH_1D
      
      CALL ENTERS('UPDATE_PARTICLES',*9999)

C Total diffusion = Da+Db
C Apparent diffusion coefficient, Da = 0.167*u*l
C Brownian diffusion, Db = Cc*k*T/(3PI*mu*dp)
C   u = velocity [from XP & INLET_FLOW]
C   l = length of element [from XP]
C   Cc == Cunningham correction factor [constant]
C   k == Boltzmann's coefficient [constant]
C   T == absolute temperature [constant]
C   mu == dynamic viscosity of air [constant]
C   dp == particle diameter [from PULMAT[1]]
C   Q == flow rate [mm^3/s]
C   rho == particle density [from PULMAT[2]]


!----------------------------------------------
c       infinity =1.d0/0.0d0
       
       dp=PULMAT(1)
       dpmm=dp*1.0d6 ! particle diameter in micro-m
!	   write(*,*) PULMAT(1), PULMAT(2)
! Choi and Kim 2007 used rho_p = 1gm/cm^3 for all simulations
!  which is = 1000kg/m^3	   
       rhop=1000.0d0  !=PULMAT(2)???????
  
!   Add constant values
       k=1.3806505d-23 !Boltzman's constant (J/K)
       T=310d0 !Temperature (K)
!      mu=0.7978d-3      !Viscosity Pa.s
       mu=1.8324d-5 !Viscosity Pa.s
       Gravity=9.81d0 ! Gravitaion const      
!     Cunningham consts
       A1=1.275d0
       A2=0.4d0
       A3=0.55d0
       Patm=101300d0 ! atmospheric pressure (Pa=N.m^-2)
       
! Get the inlet flow
       ne=NEELEM(1) !element number
       nb=NBJ(1,ne) !basis type 
       np1=NPNE(1,nb,ne) !node 1
       nv1=NVJE(1,nb,nj_flow,ne)
 !INLET_FLOW(nx)   is in mm^3/s
       Q=INLET_FLOW(nx)*XP(1,nv1,nj_flow,np1)/1.0d9 !Inlet Flowrate (m^3/s)    
       Q=SQRT(Q*Q) ! make u always positive    
       rad=XP(1,nv1,nj_radius,np1)/1000.0d0 !Radius (m)
       Area_i=PI*rad*rad

c       write(*,*) 'for each element 1..',NEELEM(0)
       DO noelem=1,NEELEM(0)
         ne=NEELEM(noelem) !element number
         nb=NBJ(1,ne) !basis type 
         np1=NPNE(1,nb,ne) !node 1
         np2=NPNE(2,nb,ne) !node 2
         nv1=NVJE(1,nb,nj_loss,ne) !version number
         nv2=NVJE(2,nb,nj_loss,ne) !version number
         
         IF(ne.eq.1)THEN ! first element = extrathroasic airways
 ! NORD(1,ne)=0
c           XP(1,1,3,np1)=150.0d0 !set increase length (mm)
           l=LENGTH_1D(NBJ,ne,NPNE,NVJE,XP) 
           XP(1,nv1,nj_radius,np1)=SQRT(70.0d3/(PI*l)) !set radius (mm)
           Qe=INLET_FLOW(nx)*XP(1,nv1,nj_flow,np1) !Elem Flowrate (mm^3/s)
           Qe= 60.0d0*Qe*1.0d-6 ! (L/min)  
           Qe=SQRT(Qe*Qe) !Qe always positive 
 ! Calculate Cunningham Const
 ! new version of the cunningham correction factor from Nazridoust.2008
           lambda=0.0673d0; ! mean free path for air (microm)
           Cc=1.0d0+(2.0d0*(lambda/dpmm))*
     &       (1.165d0 +0.483d0*exp(-(0.997d0*dpmm/(2.0d0*lambda))))
           
 ! Brownian diffusion         
           Db=Cc*k*T/(3.0d0*PI*mu*dp) 
           
           
           
 ! oral cavity efficiency
           eta=1.0d0-exp(-20.4d0*((Db*10000.0d0)**0.66d0)*(Qe**(
     &       -0.31d0))-0.000241d0*(dpmm**2.0d0)*Qe );
           
! Advective mixing
           Da=0.167d0*(u*1000.0d0)*(l*1000.0d0) !mm^2/s
! Total diffusion
           D=Da+Db*1.0d6 !Update units for Db from m^2/s => mm^2/s
           
         else ! not in oral cavity
           
           rad=XP(1,nv1,nj_radius,np1)/1000.0d0 !Radius (m)
           Area=PI*rad*rad !Elem Area (m^2)
           l=LENGTH_1D(NBJ,ne,NPNE,NVJE,XP)/1000.0d0 !Elem Length (m)
           Qe=INLET_FLOW(nx)*XP(1,nv1,nj_flow,np1)/1.0d9 !Elem Flowrate (m^3/s)    
           Qe=SQRT(Qe*Qe) !Qe always positive               
           u=Qe/Area !Velocity (m/s)
           
! new version of the cunningham correction factor from Nazridoust.2008
           lambda=0.0673d0; ! mean free path for air (microm)
           Cc=1.0d0+(2.0d0*(lambda/dpmm))*
     &       (1.165d0 +0.483d0*exp(-(0.997d0*dpmm/(2.0d0*lambda))))
           
 ! Brownian diffusion         
           Db=Cc*k*T/(3.0d0*PI*mu*dp) !(m^2/s)
! Advective mixing
           Da=0.167d0*(u*1000.0d0)*(l*1000.0d0) !(mm^2/s)
! Total diffusion
           D=Da+Db*1.0d6 !(mm^2/s)
! Stokes number
           Stk=rhop*dp*dp*Cc*u/(36.0d0*mu*rad) !(dimless)! Get element unit vectors
! unit vector=|{ex,ey,ez}|=1 		
           ex=XP(1,2,1,np2) ! x component of unit vector
           ey=XP(1,2,2,np2) ! y component of unit vector
           ez=XP(1,2,3,np2) ! z component of unit vector
!get the generation of the element
           gen =NORD(1,ne)
           tempA=1.0d0
           
! get unit vector of previous element.
           IF(gen.EQ.1) THEN ! First element
             theta=0.0d0 ! Branching Angle
! gen of elem 1 is 1
           ELSE !previous searh to find previous generation element
             SEARCH=.false. !Init search
             theta=0.0d0
             ne0 =NXI(-1,1,ne) !previous element
             gen0=  NORD(1,ne0) !generation of ne0
             IF(gen0.EQ.gen)THEN !Same generation?
               SEARCH=.TRUE. !Then search
             ENDIF
             DO WHILE (SEARCH) !search backwards to find previous bifurcation
               SEARCH=.FALSE. !Stop searching unless...?
               ne0 =NXI(-1,1,ne0) !previous element
               gen0=NORD(1,ne0) !generation of ne0
               IF(gen0.EQ.gen)THEN !Same generation?
                 SEARCH=.TRUE. !Then search
               ENDIF
             ENDDO
             np3 =NPNE(1,nb,ne0) !node 1 of ne0
 ! unit vector=|{vx,vy,vz}|=1 	of ne0	
             vx=XP(1,2,1,np3) ! x component of unit vector
             vy=XP(1,2,2,np3) ! y component of unit vector
             vz=XP(1,2,3,np3) ! z component of unit vector
 !COS(theta)= (v dot u)
             tempA=0.0d0
             tempA=(vx*ex)+(vy*ey)+(vz*ez)
             if ((abs(tempA)-1.0d0).LE.ZERO_TOL) then
               theta=acos((tempA-ZERO_TOL))
             else
               theta=ACOS(tempA)
             endif
             
           ENDIF !end find branching angle	
! Gravity Angle 
! gravity unit vector ={0,0,1}, therefore
           phi_g=ACOS(ez)
           
           
!-----------------------------------
           
           IF(INLET_FLOW(nx).GE.0) THEN ! Inspiration
             etaI=1.0d0-(1.0d0/(1.0d0+14.01d0*Stk**(1.977d0))) ! impaction efficiency (Insp)
           ELSE ! Exhalation
             etaI=1.0d0-(1.0d0/(1.0d0+12.57d0*((DSQRT(Stk)*
     &         SIN(theta))**2.5d0))) ! impaction efficiency (Expir)
           ENDIF ! END Insp/Expir
           
           
           Tau=rhop*dp*dp/(18d0*mu) ! Relaxation time (s)    
           temp=2.0d0*Gravity*Tau*l*SIN(phi_g)*DSQRT((Area/PI))*Qe
           etaS=1.0d0-EXP(-temp)
           Delta = PI*Db/(4.0d0*Qe)
           
           
           etaD=1.0d0-0.819d0*EXP(-14.63d0*Delta)
     &       -0.0975d0*EXP(-89.225d0*Delta)
     &       -0.0325d0*EXP(-228d0*Delta)
     &       -0.0509d0*EXP(-125.9d0*(Delta**(2.0d0/3.0d0)))         
           
!Update units for eta from 1/m => 1/mm 
           etaI=etaI
           etaS=etaS
 ! etaD=etaD*1.0d-6
           
! Calculate the total deposition efficiency 
           eta=etaI+etaS+etaD+etaI*etaS*etaD-etaI*etaS+etaS*etaD+etaI
     &       *etaD
           
           
         ENDIF
         
         CE(1,ne)=D
         XP(1,nv1,nj_loss,np1)=eta
         XP(1,nv2,nj_loss,np2)=eta

c         write(*,*) 'eta',eta
         
       ENDDO ! End loop elements
       
       CALL EXITS('UPDATE_PARTICLES')
       RETURN
 9999  CALL ERRORS('UPDATE_PARTICLES',ERROR)
       CALL EXITS('UPDATE_PARTICLES')
       RETURN 1
       END
      
