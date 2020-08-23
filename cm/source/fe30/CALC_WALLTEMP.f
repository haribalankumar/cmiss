      SUBROUTINE CALC_WALLTEMP(np,nr,nx,beta,Clumen,dH,flux,Kd,L,Lpvc,
     &  tissue_depth,Tlumen,Tlumen0,Twall,CP,Vevap,XP,TUBE,ERROR,*)

C#### Subroutine: CALC_WALLTEMP
C###  Description:
C###    CALC_WALLTEMP finds the temperature at the Mucus-Air-Interface
C###    (MAI) in an airway, given the temperature and water vapour
C###    concentration at the airway centre, and the wall temperature.
C***  Created by Merryn Howatson Tawhai, May 1998

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER nr,np,nx
      REAL*8 beta,Clumen,CP(NMM),dH,Kd,L,Lpvc,tissue_depth,
     &  Tlumen,Tlumen0,Twall,XP(NKM,NVM,NJM)
      LOGICAL CELL,TUBE
      CHARACTER ERROR*(*)
!    Local Variables
      INTEGER i,N,nl,nl1,nl2
      REAL*8 beta0,Cmai,Cmax,Cp_w,Cp_t,Cp_pvc,depth,depth0,
     '  devap,DIVIDE,dTdr,A(3,3),B(3),GAMMA(3),heat_term,Ke,Kl,Kpvc,
     &  Kt,Kw,K_dtdr,Me1,Me2,Rm,R1,R2,Rw,rho_pvc,rho_t,rho_w,
     '  temp_vol,thet,Tmai,Tmai0,T_1,T_2,Vevap,vol_limit,Vsource,X(3),
     '  CONC_EVAL,gam,flux,MAX_EVAP,MAX_COND,gland_flux

      CALL ENTERS('CALC_WALLTEMP',*9999)

      depth=CP(2)
      depth0=CP(2)
      dely = tissue_depth
      RH_MAI = 1.d0
      gam=MAX(2.d0,beta) ! used for the humidity
      beta0=beta
      IF(INLET_FLOW(nx).GE.0.d0)THEN !inspiration
        gam  =MAX(2.d0,WHT(1)*beta) !used for humidity
        beta0=MAX(2.d0,WHT(2)*beta) !used for temperature
      ELSE
        gam  =MAX(2.d0,WHT(3)*beta) !used for temperature
        beta0=MAX(2.d0,WHT(4)*beta) !used for temperature
      ENDIF

      Tmai0=CP(4)
      T_1=CP(5)
      T_2=CP(6)
      Rm=CP(1)-depth0
      R1=CP(1)
      R2=CP(1)+Lpvc
      Rw=CP(1)+Lpvc+dely

      Kw=0.58d-6    !thermal cond. through ASL=water (kJ/mm/s/K)
      Kt=0.2d-6    !tissue, from Ozen et al. Burns; 34(1):45-49, 2008
      Kpvc=0.92d-7 !thermal conductivity pvc (from F&PH,confirmed)
      Kl=PULMAT(5) !thermal conductivity of air (user specified)
      rho_w=0.1d-2     !density of water (ASL) [g.mm-3]
      rho_t=0.1d-2     !same as water (Perry) [g.mm-3]
      rho_pvc=0.171d-2 !density of pvc (from F&P, confirmed)
      Cp_w=0.4184d-2   !specific heat of ASL==water (Perry) [kJ/g/K]
      Cp_t=0.35d-2     !Valvano(1984)=K_t/(rho_t*alpha)
      Cp_pvc=0.105d-2  !specific heat (from F&PH, confirmed)

      vol_limit=PI*(R1**2.d0-Rm**2.d0)*L
      temp_vol=vol_limit
      Cmax=CONC_EVAL(Tlumen)
      Cmai=RH_MAI*CONC_EVAL(CP(4))
      IF(Clumen.GT.Cmax) Clumen=Cmax

      IF(INLET_FLOW(nx).GE.0.d0)THEN
        ! max conc in the lumen - mean concentration in the lumen
        MAX_EVAP = DABS((Cmai-(Cmai-Clumen))*PI*R1**2*L/1.d9/rho_w) !max mm^3 to evap
        Vevap=Kd*gam*(Cmai-Clumen)*2.d0*PI*L*DT/rho_w/1.d9 !mm^3
        Vevap=MIN(Vevap,MAX_EVAP)
      ELSE
        IF(Cmai.GT.Clumen)THEN
          Vevap=0.d0 !not allowing evaporation during expiration.
                     !volume is minimal but causes numerical errors
        ELSE
          Vevap=Kd*gam*(Cmai-Clumen)*2.d0*PI*L*DT/rho_w/1.d9 !mm^3
        ENDIF
      ENDIF

      flux = Vevap/(2.d0*PI*L*DT)/Rm
c      if(np.eq.1) write(*,*) 'Vevap',Vevap,Kd,gam,Cmai,Clumen

c      IF(TUBE)THEN
c        Vsource=0.d0
c      ELSE
      IF(ITYP12(nr,1).EQ.2)THEN
        ! note that the ODEMODEL code uses micro-metres/sec

c        CALL ODEMODEL(depth,DT,flux,XP,ERROR,*9999)
c        devap=R1-(R1**2.d0-(vol_limit-Vevap)/(PI*L))**0.5d0


c        IF(np.EQ.1)THEN
c          WRITE(OP_STRING,'(4D12.4)') depth0,depth0-Vevap/(2.d0*PI*R1
c     &      *L),flux,depth
c          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c        ENDIF
      ELSE
c      IF(TUBE)THEN
c        Vsource=0.d0
c      ELSE
c        IF(devap.LT.0.006d0)THEN

        Vsource=2.d0*PI*R1*L*DT*(SF_RATE/6000.d0)
        devap=R1-(R1**2.d0-(vol_limit-Vevap+Vsource)/(PI*L))**0.5d0
        IF(devap.GT.0.007d0)THEN
          Vsource=(R1**2.d0-Rm**2.d0)*PI*L+Vevap-vol_limit
        ENDIF
      ENDIF !cell
C ENDIF ! tube     
      
c        ELSE IF(devap.GT.0.006d0)THEN
c          Vsource=-2.d0*PI*R1*L*DT*(SF_RATE/6000.d0)
c          devap=R1-(R1**2.d0-(vol_limit-Vevap+Vsource)/(PI*L))**0.5d0
c          IF(devap.LT.0.006d0)THEN
c            Vsource=(R1**2.d0-Rm**2.d0)*PI*L+Vevap-vol_limit
c          ENDIF
c        ELSE
c          Vsource=0.d0
c        ENDIF
c      ENDIF !TUBE
c      IF((vol_limit-Vevap+Vsource).LE.ZERO_TOL)THEN !all evaporates
c        heat_term=(vol_limit+Vsource)*dH*rho_w/(2.d0*PI*Rm*L*DT)
c        temp_vol=0.d0
c        depth=0.d0
c      ELSE !some evaporation or condensation

      heat_term=Vevap*dH*rho_w/(2.d0*PI*Rm*L*DT)
      
      IF(ITYP12(nr,1).EQ.1)THEN
        temp_vol=temp_vol-Vevap+Vsource
        depth=R1-(R1**2.d0-temp_vol/(PI*L))**0.5d0
      ENDIF

      DO nl1=1,3
        DO nl2=1,3
          A(nl1,nl2)=0.d0
        ENDDO !nl2
        B(nl1)=0.d0
      ENDDO !nl1
      thet=2.d0/3.d0
      nl=0
      IF(depth0.GT.ZERO_TOL)THEN
        nl=nl+1
        ! Ke=2*PI*Kt*(Rm+(Rn-Rm)/2)/(Rn-Rm) *[1 -1;-1 1]
        Ke=2.d0*PI*Kw/(R1-Rm)*(Rm+0.5d0*(R1-Rm))
        ! Me=2*PI*rho*Cp*(Rn-Rm)*etc
        Me1=2.d0*PI*rho_w*Cp_w*(R1-Rm)*Rm/6.d0
        Me2=2.d0*PI*rho_w*Cp_w*(R1-Rm)*(R1-Rm)/12.d0
        A(nl,nl)=Me1*2.d0+Me2+DT*thet*Ke
        A(nl,nl+1)=Me1+Me2-DT*thet*Ke
        A(nl+1,nl)=Me1+Me2-DT*thet*Ke
        A(nl+1,nl+1)=Me1*2.d0+Me2*3.d0+DT*thet*Ke
        B(nl)=-Ke*Tmai0+Ke*T_1
        B(nl+1)=Ke*Tmai0-Ke*T_1
        K_dtdr=Kw
      ENDIF
c      IF(Lpvc.GT.ZERO_TOL)THEN
c        nl=nl+1
c        Ke=2.d0*PI*Kpvc/(R2-R1)*(R1+0.5d0*(R2-R1))
c        Me1=2.d0*PI*rho_pvc*Cp_pvc*(R2-R1)*R1/6.d0
c        Me2=2.d0*PI*rho_pvc*Cp_pvc*(R2-R1)*(R2-R1)/12.d0
c        A(nl,nl)=A(nl,nl)+Me1*2.d0+Me2+DT*thet*Ke
c        A(nl,nl+1)=A(nl,nl+1)+Me1+Me2-DT*thet*Ke
c        A(nl+1,nl)=A(nl+1,nl)+Me1+Me2-DT*thet*Ke
c        A(nl+1,nl+1)=A(nl+1,nl+1)+Me1*2.d0+Me2*3.d0+DT*thet*Ke
c        B(nl)=B(nl)-Ke*T_1+Ke*T_2
c        B(nl+1)=B(nl+1)+Ke*T_1-Ke*T_2
c        IF(nl.EQ.1) K_dtdr=Kpvc
c      ENDIF

c      nl=nl+1
c      Ke=2.d0*PI*Kt/(Rw-R2)*(R2+0.5d0*(Rw-R2))
c      Me1=2.d0*PI*rho_t*Cp_t*(Rw-R2)*R2/6.d0
c      Me2=2.d0*PI*rho_t*Cp_t*(Rw-R2)*(Rw-R2)/12.d0
c      A(nl,nl)=A(nl,nl)+Me1*2.d0+Me2+DT*thet*Ke
c      A(nl,nl+1)=A(nl,nl+1)+Me1+Me2-DT*thet*Ke
c      A(nl+1,nl)=A(nl+1,nl)+Me1+Me2-DT*thet*Ke
c      A(nl+1,nl+1)=A(nl+1,nl+1)+Me1*2.d0+Me2*3.d0+DT*thet*Ke
c      B(nl)=B(nl)-Ke*T_1+Ke*Twall
c      B(nl+1)=B(nl+1)+Ke*T_1-Ke*Twall

      nl=nl+1
      Ke=2.d0*PI*Kt/(Rw-R2)*(R2+0.5d0*(Rw-R2))
      Me1=2.d0*PI*rho_t*Cp_t*(Rw-R2)*R2/6.d0
      Me2=2.d0*PI*rho_t*Cp_t*(Rw-R2)*(Rw-R2)/12.d0
      A(nl,nl)=A(nl,nl)+Me1*2.d0+Me2+DT*thet*Ke

      B(nl)=B(nl)-Ke*T_2+Ke*Twall
      
      IF(nl.EQ.1) K_dtdr=Kt

      dTdr=(-Kl*beta0*(Tmai0-Tlumen0)/Rm-heat_term)/K_dtdr

      B(1)=B(1)+dTdr*2.d0*PI*Rm*K_dtdr !flux b.c. at MAI
      N=nl

      IF(N.EQ.1)THEN
        X(1)=B(1)/A(1,1)
      ELSE IF(N.EQ.2)THEN
        X(2)=(A(1,1)*B(2)-A(2,1)*B(1))/
     '    (A(1,1)*A(2,2)-A(1,2)*A(2,1))
        X(1)=(B(1)-A(1,2)*X(2))/A(1,1)
      ELSE IF(N.GE.3)THEN
        DIVIDE=A(1,2)
        X(1)=B(1)/DIVIDE
        DO nl=2,N
          GAMMA(nl)=A(nl-1,3)/DIVIDE
          DIVIDE=A(nl,2)-A(nl,1)*GAMMA(nl)
          X(nl)=(B(nl)-A(nl,1)*X(nl-1))/DIVIDE
        ENDDO
        DO nl=N-1,1,-1
          X(nl)=X(nl)-GAMMA(nl+1)*X(nl+1)
        ENDDO
      ENDIF !N

      CP(8)=depth

      IF(depth0.GT.ZERO_TOL.AND.Lpvc.GT.ZERO_TOL)THEN
        Tmai=Tmai0+X(1)*DT
        CP(9)=0.d0
        CP(10)=CP(4)+X(1)*DT
        CP(11)=CP(5)+X(2)*DT
        CP(12)=CP(6)+X(3)*DT
      ELSE IF(depth0.GT.ZERO_TOL.AND.Lpvc.LE.ZERO_TOL)THEN
        Tmai=Tmai0+X(1)*DT
        CP(9)=CONC_EVAL(Tmai)*RH_MAI
        CP(10)=CP(4)+X(1)*DT
        CP(11)=CP(5)+X(2)*DT
        CP(12)=CP(11)
      ELSE IF(depth0.LE.ZERO_TOL.AND.Lpvc.GT.ZERO_TOL)THEN
        Tmai=Tmai0+X(1)*DT
        CP(9)=0.d0
        CP(11)=CP(5)+X(1)*DT
        CP(10)=CP(11)
        CP(12)=CP(6)+X(2)*DT
      ELSE IF(depth0.LE.ZERO_TOL.AND.Lpvc.LE.ZERO_TOL)THEN
        Tmai=Tmai0+X(1)*DT
        CP(9)=CONC_EVAL(Tlumen0)*RH_MAI
        CP(12)=CP(6)+X(1)*DT
        CP(10)=CP(12)
        CP(11)=CP(12)
      ENDIF
   
c      DO i=10,12
c !i.e. for et tube with time-varying inspired temperature
c !the wall T must be able to be lower than the inspired T
c        IF(FLOW.GT.0.d0.AND.CP(i).LT.Tlumen0)THEN
c          CP(i)=Tlumen0
c          IF(CP(9).NE.0.d0) CP(9)=CONC_EVAL(Tlumen0)*RH_MAI
c        ENDIF
c        IF(FLOW.LT.0.d0.AND.CP(i).LT.CP(i-6)) CP(i)=CP(i-6)
c        IF(CP(i).GT.TMAX)THEN
c          CP(i)=TMAX
c          IF(CP(9).NE.0.d0) CP(9)=CONC_EVAL(TMAX)*RH_MAI
c        ENDIF
c      ENDDO

      CALL EXITS('CALC_WALLTEMP')
      RETURN
 9999 CALL ERRORS('CALC_WALLTEMP',ERROR)
      CALL EXITS('CALC_WALLTEMP')
      RETURN 1
      END



