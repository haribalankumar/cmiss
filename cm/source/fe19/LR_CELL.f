      SUBROUTINE LR_CELL(TIME,Y,DY,CONTROL,MODEL,SIZES,VARIANT,DERIVED,
     '  PARAM,PROTOCOL,AII,AIO,ARI,ARO,ERR_CODE)

C#### Subroutine: LR_CELL
C###  Description:
C###    Computes RHS of Luo-Rudy equations.
C###    Voltage is in mV & time in ms.

      IMPLICIT NONE

      INCLUDE 'cell_lr.inc'
      INCLUDE 'cell_reserved.inc'

!     Parameter List
      INTEGER SIZES(11)
      INTEGER CONTROL(*),MODEL(*),VARIANT,AII(*),AIO(*),ERR_CODE
      REAL*8 TIME(*),Y(*),PARAM(*),PROTOCOL(*),ARI(*),DY(*),
     '  DERIVED(*),ARO(*)
      !common blocks
      REAL*8        I(19),ICa !DPN 3/10/97 - extending I()
      COMMON /curr/ I,ICa
!     Local Variables
      REAL*8 VRatio_MJ, VRatio_MN, flux
c      REAL*8 Iv_temp
c      INTEGER count,ISWTCH(50)

C *** DPN 22 June 1999
      REAL*8 CURRENTS(27),RATES(21)

c      SINGLE_CELL = .TRUE.

      PARAM(RTONF) = (PARAM(Temp)+273.d0)*8314.472d0/96485.341d0 !mV

      VRatio_MJ = PARAM(VolJSR)/PARAM(VolMyo)
      VRatio_MN = PARAM(VolNSR)/PARAM(VolMyo)

C *** Calculate the stimulus current
CC     fist the stuff for distributed models - should only be greater
CC     than zero for distributed models
C      DERIVED(Stimulus)=-(PROTOCOL(PseudoIs)+PROTOCOL(Is1current))/
C     '  (PARAM(Cm)*PARAM(Am))
CC     Then the rest of it - should have no contribution to distributed
CC     models
C      IF ((T.GE.PARAM(TPS)).AND.
C     '  (T.LE.PARAM(TPS)+PARAM(TP))) THEN
C        DERIVED(Stimulus)=DERIVED(Stimulus)+PARAM(Istim)
C      ELSE
C        IF (T.GT.PARAM(TPS)+PARAM(TP)) THEN
C          IF (CONTROL(NUM_APPLIED_STIMULII).LT.
C     '      INT(PARAM(numStim))) THEN
C            PARAM(TPS) = PARAM(TPS) + PARAM(dtp)
C            CONTROL(NUM_APPLIED_STIMULII)=
C     '        CONTROL(NUM_APPLIED_STIMULII)+1
C          ENDIF
C        ENDIF
CC ***   This can go now ???
Cc        DERIVED(Stimulus)=0.d0
C      ENDIF

      CURRENTS(Stimulus)=0.0d0
      CURRENTS(Stimulus)=-PROTOCOL(PseudoIs)
      IF(TIME(TCell).GE.PROTOCOL(Is1start).AND.
     '  TIME(TCell).LT.PROTOCOL(Is1stop))
     '  CURRENTS(Stimulus)=CURRENTS(Stimulus)+PROTOCOL(Is1current)
      IF(TIME(TCell).GE.PROTOCOL(Is2start).AND.
     '  TIME(TCell).LT.PROTOCOL(Is2stop))
     '  CURRENTS(Stimulus)=CURRENTS(Stimulus)+PROTOCOL(Is2current)
      IF(PROTOCOL(IsFreqPeriod).GT.1.0d-6) THEN
        IF(DMOD(TIME(TCell),PROTOCOL(IsFreqPeriod)).LT.
     '    PROTOCOL(IsFreqDuration))
     '    CURRENTS(Stimulus)=CURRENTS(Stimulus)+PROTOCOL(IsFreqMag)
      ENDIF
      CURRENTS(Stimulus)=CURRENTS(Stimulus)/PARAM(Am)
c      Istim=1000.0d0*Istim/PARAM(N98_Am) !Convert to nA.mm^-2

C *** Compute rate constants
      CALL LR_RATES_CELL(TIME,Y,DY,CONTROL,MODEL,SIZES,VARIANT,CURRENTS,
     '  PARAM,PROTOCOL,AII,AIO,ARI,RATES,ERR_CODE)

C *** Compute ionic currents
      CALL LR_CURRENTS_CELL(TIME,Y,DY,CONTROL,MODEL,SIZES,VARIANT,
     '  CURRENTS,PARAM,PROTOCOL,AII,AIO,ARI,RATES,ERR_CODE)

      flux = CURRENTS(Irel)*VRatio_MJ+
     '  (-CURRENTS(Iup)+CURRENTS(Ileak))*VRatio_MN

! Compute o.d.e. RHS
C *** DPN 15/04/98 - adding switch to each current in the total
c      Iv_temp = ISWTCH(14)*LR_IK1+ISWTCH(21)*LR_IKp+
c     +  ISWTCH(18)*LR_INaK
c     +  +ISWTCH(26)*LR_IpCa+ISWTCH(2)*LR_ICab+
c     +  ISWTCH(4)*LR_INab
c      Iv_temp = ISWTCH(24)*Iv_temp ! turns all time indep. currents off
      DY(Vm) =  -(CURRENTS(Stimulus)+CURRENTS(INa)+CURRENTS(ICaL)+
     '  CURRENTS(IK)+CURRENTS(INaCa)+CURRENTS(InsCa)+CURRENTS(Iv))/
     '  PARAM(Cm)
      DY(m) = RATES(alpha_m)*(1.d0-Y(m)) - RATES(beta_m)*Y(m)
      DY(h) = RATES(alpha_h)*(1.d0-Y(h)) - RATES(beta_h)*Y(h)
      DY(j) = RATES(alpha_j)*(1.d0-Y(j)) - RATES(beta_j)*Y(j)
      DY(d) = RATES(alpha_d)*(1.d0-Y(d)) - RATES(beta_d)*Y(d)
      DY(f) = RATES(alpha_f)*(1.d0-Y(f)) - RATES(beta_f)*Y(f)
      DY(x) = RATES(alpha_x)*(1.d0-Y(x)) - RATES(beta_x)*Y(x)
      DY(Nai) =  -(CURRENTS(I_Na)*PARAM(Acap))/
     '  (PARAM(VolMyo)*Farad)
      DY(Cai) =  -(CURRENTS(I_Ca)*PARAM(Acap))/
     '  (PARAM(VolMyo)*2.d0*Farad)+flux
      DY(Ki) =  -(CURRENTS(I_K)*PARAM(Acap)) /
     '  (PARAM(VolMyo)*Farad)
      DY(CaJSR) =  -CURRENTS(F_JSR)
      DY(CaNSR) =  -CURRENTS(F_NSR)
      DY(Cab) =  -((CURRENTS(ICaLCa)+CURRENTS(IpCa)+CURRENTS(ICab))*
     '  1.d-3)/2.d0
      DY(Cao) = 0.0d0 ![Ca]o constant

      IF(CONTROL(RETURN_CURRENTS).NE.0) THEN
        DERIVED(Stimulus)=CURRENTS(Stimulus)
        DERIVED(INa)=CURRENTS(INa)
        DERIVED(ICaLCa)=CURRENTS(ICaLCa)
        DERIVED(ICaLK)=CURRENTS(ICaLK)
        DERIVED(ICaLNa)=CURRENTS(ICaLNa)
        DERIVED(ICaL)=CURRENTS(ICaL)
        DERIVED(IK)=CURRENTS(IK)
        DERIVED(IK1)=CURRENTS(IK1)
        DERIVED(IKp)=CURRENTS(IKp)
        DERIVED(INaCa)=CURRENTS(INaCa)
        DERIVED(INaK)=CURRENTS(INaK)
        DERIVED(InsK)=CURRENTS(InsK)
        DERIVED(InsNa)=CURRENTS(InsNa)
        DERIVED(InsCa)=CURRENTS(InsCa)
        DERIVED(IpCa)=CURRENTS(IpCa)
        DERIVED(ICab)=CURRENTS(ICab)
        DERIVED(INab)=CURRENTS(INab)
        DERIVED(Iv)=CURRENTS(Iv)
        DERIVED(Irel)=CURRENTS(Irel)
        DERIVED(Ileak)=CURRENTS(Ileak)
        DERIVED(Iup)=CURRENTS(Iup)
        DERIVED(I_Na)=CURRENTS(I_Na)
        DERIVED(I_Ca)=CURRENTS(I_Ca)
        DERIVED(I_K)=CURRENTS(I_K)
        DERIVED(F_NSR)=CURRENTS(F_NSR)
        DERIVED(F_JSR)=CURRENTS(F_JSR)
        DERIVED(Itr)=CURRENTS(Itr)
      ENDIF

      RETURN
      END


