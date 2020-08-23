      SUBROUTINE LR_CURRENTS_CELL(TIME,Y,DY,CONTROL,MODEL,SIZES,VARIANT,
     '  CURRENTS,PARAM,PROTOCOL,AII,AIO,ARI,RATES,ERR_CODE)

C#### Subroutine: LR_CURRENTS_CELL
C###  Description:
C###    Computes ionic currents for Luo-Rudy model

C *** total current for each ions
C *** I_Na  = sum of currents containing Na ions
C *** I_Ca  = sum of currents containing Ca ions
C *** I_K   = sum of currents containing K ions
C *** F_JSR = sum of fluxes moving into JSR
C *** F_NSR = sum of fluxes moving into NSR

      IMPLICIT NONE

      INCLUDE 'cell_lr.inc'
      INCLUDE 'cell_reserved.inc'

!     Parameter List
      INTEGER SIZES(11)
      INTEGER CONTROL(*),MODEL(*),VARIANT,AII(*),AIO(*),ERR_CODE
      REAL*8 TIME(*),Y(*),PARAM(*),PROTOCOL(*),ARI(*),DY(*),
     '  CURRENTS(27),RATES(21)
!     Local Variables
      REAL*8 ENa
      REAL*8 EK,GGK,Xi
      REAL*8 EK1,GGK1,K10
      REAL*8 EKp,Kp
      REAL*8 fNaK,sigma
      REAL*8 EnsCa,Vns,IInsNa,IInsK
      REAL*8 ECaN
      REAL*8 ENaN
      REAL*8 CaJSRbuffer
      REAL*8 Grel,GGrel
      REAL*8 t_CICR
      REAL*8 Kleak
      REAL*8 VRatio_JN
      REAL*8       CSQN
      REAL*8        KmCa
      LOGICAL DEBUG

      DEBUG = .FALSE.

! Initialising the volume parameters
      VRatio_JN= PARAM(VolNSR)/PARAM(VolJSR)

! Ca buffers in the myoplasm
      !??? DPN - move to common block ???
      !TRPNTRPN = 70.d-3               !mmol/L
      !CMDNCMDN = 50.d-3               !mmol/L
      !KmTRPN   = 0.5d-3               !mmol/L
      !KmCMDN   = 2.38d-3              !mmol/L
c      TRPN     = TRPNTRPN*(LR_Cai/(LR_Cai+KmTRPN))
c      CMDN     = CMDNCMDN*(LR_Cai/(LR_Cai+KmCMDN))


! Calcium concentration after taking the buffer TRPN and CMDN into
! account
C dpn 03/06/98 - uncomment
!      Cabuffer = Cai-TRPN-CMDN
!      IF(t.GE.Tstim)THEN
!        IF(Cabuffer.GT.1.d-6)THEN
!          Cai = Cabuffer
!        ELSE
!          Cai = 1.d-6
!        ENDIF
!      ENDIF
c      Cabuffer = Cai-TRPN-CMDN
c      IF(t.GE.TPS)THEN
c        IF(Cabuffer.GT.1.d-6)THEN
c          Cai = Cabuffer
c        ELSE
c          Cai = 1.d-6
c        ENDIF
c      ENDIF
C dpn 03/06/98

! Ca buffer in JSR and CSQN
      CSQN      = PARAM(CSQNCSQN)*(Y(CaJSR)/(Y(CaJSR)+
     '  PARAM(KmCSQN)))


! Ca flux in the SR after taking CSQN buffering into account
C dpn 03/06/98 - uncomment
!      CaJSRbuffer = CaJSR-CSQN
!      IF(t.GE.Tstim)THEN
!        IF(CaJSRbuffer.GT.1.d-3)THEN
!          CaJSR = CaJSRbuffer
!        ELSE
!          CaJSR = 1.d-3
!        ENDIF
!      ENDIF
      CaJSRbuffer = Y(CaJSR)-CSQN
      IF(TIME(TCell).GE.PARAM(TPS))THEN
        IF(CaJSRbuffer.GT.1.d-3)THEN
          Y(CaJSR) = CaJSRbuffer
        ELSE
          Y(CaJSR) = 1.d-3
        ENDIF
      ENDIF
C dpn 03/06/98

C *** Ionic currents in the sarcolemma
      !Fast Na current INa
      !Na reversal potential (mV)
      ENa  = PARAM(RTONF)*dlog(PARAM(Nao)/Y(Nai))
C *** DPN 29 April 1998 - grab max INa
c      I(1) = GGNa*m*m*m*h*j*(V-ENa)
      RATES(IINa) = PARAM(GGNa) * (Y(Vm)-ENa)
      CURRENTS(INa) = Y(m)*Y(m)*Y(m) * Y(h) * Y(j) * RATES(IINa)

! L-type Ca current ICa_t
      KmCa = PARAM(KmCa_L) !0.6d-3    !mmol/L
      RATES(fCa) = 1.d0/(1.d0+(Y(Cai)/KmCa)**2)
      !Ca-dep. inact.n gate of L-type Ca

      RATES(IICa) = PARAM(PCa)*4.d0*(Y(Vm)/PARAM(RTONF))*
     '  Farad*
     '  (PARAM(gCai)*Y(Cai)*dexp(2.d0*Y(Vm)/PARAM(RTONF))
     '  -PARAM(gCao)*Y(Cao))/(dexp(2.d0*Y(Vm)/PARAM(RTONF))
     '  -1.d0)

      RATES(IICaNa) = PARAM(PNa)*(Y(Vm)/PARAM(RTONF))*Farad*
     '  (PARAM(gNai)*Y(Nai)*dexp(Y(Vm)/PARAM(RTONF))-
     '  PARAM(gNao)*PARAM(Nao))/(dexp(Y(Vm)/
     '  PARAM(RTONF))-1.d0)

      RATES(IICaK)  = PARAM(PK)*(Y(Vm)/PARAM(RTONF))*Farad*
     '  (PARAM(gKi)*Y(Ki)*dexp(Y(Vm)/PARAM(RTONF))-
     '  PARAM(gKo)*PARAM(Ko))/(dexp(Y(Vm)/PARAM(RTONF))
     '  -1.d0)

      CURRENTS(ICaLCa) = Y(d)*Y(f)*RATES(fCa)*RATES(IICa)
      CURRENTS(ICaLK)  = Y(d)*Y(f)*RATES(fCa)*RATES(IICaK)
      CURRENTS(ICaLNa) = Y(d)*Y(f)*RATES(fCa)*RATES(IICaNa)

      CURRENTS(ICaL)=CURRENTS(ICaLCa)+CURRENTS(ICaLK)+CURRENTS(ICaLNa)

! Time-dep. K current IK
      Xi   = 1.d0/(1.d0+dexp((Y(Vm)-56.26d0)/32.1d0))
                                      !V-dep. K-inactivation gate
      GGK  = 0.00282d0*dsqrt(PARAM(Ko)/5.4d0)  !millisiemens/mm^2
      EK   = PARAM(RTONF)*dlog((PARAM(Ko)+PARAM(PNaK)*
     '  PARAM(Nao))/(Y(Ki)+PARAM(PNaK)*Y(Nai)))
                                      !K reversal potential
C *** DPN 29 April 1998 - grab the max I(K)
c      I(3) = GGK*Xi*X*X*(V-EK)
      RATES(IIK) = GGK*(Y(Vm)-EK)
      CURRENTS(IK) = RATES(IIK) * Xi * Y(x) * Y(x)

! Time-indep. K current IK1
      GGK1 = 0.0075d0*dsqrt(PARAM(Ko)/5.4d0)   !millisiemens/mm^2
      EK1  = PARAM(RTONF)*dlog(PARAM(Ko)/Y(Ki))
      K10  = RATES(alpha_K1)/(RATES(alpha_K1)+RATES(beta_K1))
C *** DPN 29 April 1998 - grab max I(K1)
c      I(4) = GGK1*K10*(V-EK1)
      RATES(IIK1) = GGK1 * (Y(Vm) - EK1)
      CURRENTS(IK1) = RATES(IIK1) * K10

! Plateau K current IKp
      EKp  = EK1                      !Kp reversal potential
      Kp   = 1.d0/(1.d0+dexp((7.488d0-Y(Vm))/5.98d0))
      CURRENTS(IKp) = PARAM(GGKp)*Kp*(Y(Vm)-EKp)

! Na-Ca exchanger current INaCa
      !DPN 7/10/97 - moved to common block
      !kNaCa = 2000.d0   !scaling factor of INaCa (uA/uF)
      !KmNa  = 87.5d0    !half sat.n conc of Na channel (mmol/L)
      KmCa  = PARAM(KmCa_NaCa) !1.38d0
!half sat.n conc of Ca channel (mmol/L)
      CURRENTS(INaCa)  = PARAM(kNaCa)*(dexp(PARAM(eta)*
     '  Y(Vm)/PARAM(RTONF))*Y(Nai)**3*Y(Cao)-
     '  dexp((PARAM(eta)-1.d0)*
     '  Y(Vm)/PARAM(RTONF))*PARAM(Nao)**3*Y(Cai))
     '  /((PARAM(KmNa)**3+PARAM(Nao)**3)*
     '  (KmCa+Y(Cao))*(1.d0+PARAM(ksat)*
     '  dexp((PARAM(eta)-1.d0)*Y(Vm)/PARAM(RTONF))))

! Na-K pump INaK
      sigma = (1.d0/7.d0)*(dexp(PARAM(Nao)/67.3d0)-1.d0)
                                     ![Na]o-dependence factor of fNak
      fNaK  = 1.d0/(1.d0+0.1245d0*dexp(-0.1d0*Y(Vm)/PARAM(RTONF))
     '        +0.0365d0*sigma*dexp(-Y(Vm)/PARAM(RTONF)))
      CURRENTS(INaK) = PARAM(IINaK)*fNaK*PARAM(Ko)/
     '  ((PARAM(Ko)+PARAM(KmKo))*(1.d0+
     '  dsqrt((PARAM(KmNai)/Y(Nai))**3)))

! Non-specific Ca-activated current InsCa
      EnsCa  = PARAM(RTONF)*dlog((PARAM(Ko)+PARAM(Nao))/
     '  (Y(Ki)+Y(Nai)))
      Vns    = Y(Vm)-EnsCa
      IInsK  = PARAM(PnsCa)*(Vns/PARAM(RTONF))*Farad*
     '  (PARAM(gNai)*Y(Nai)*dexp(Vns/PARAM(RTONF))
     '  -PARAM(gNao)*PARAM(Nao))/(dexp(Vns/
     '  PARAM(RTONF))-1.d0)
      CURRENTS(InsK) = IInsK /(1.d0+(PARAM(KmnsCa)/Y(Cai))**3)
      IInsNa = PARAM(PnsCa)*(Vns/PARAM(RTONF))*Farad*
     '  (PARAM(gKi)*Y(Ki)*dexp(Vns/PARAM(RTONF))
     '  -PARAM(gKo)*PARAM(Ko))/(dexp(Vns/PARAM(RTONF))
     '  -1.d0)
      CURRENTS(InsNa) = IInsNa/(1.d0+(PARAM(KmnsCa)/Y(Cai))**3)
      CURRENTS(InsCa) = CURRENTS(InsK) + CURRENTS(InsNa)

! Sarcolemmal Ca pump IpCa
      CURRENTS(IpCa)   = PARAM(IIpCa)*(Y(Cai)/
     '  (PARAM(KmpCa)+Y(Cai)))

! Ca background current ICab
      ECaN = 0.5d0*PARAM(RTONF)*dlog(Y(Cao)/Y(Cai))!Nernst potential of Ca
      CURRENTS(ICab) = PARAM(GGCab)*(Y(Vm)-ECaN)

! Na background current INab
      ENaN   = ENa                    !Nernst potential of Na
      CURRENTS(INab) = PARAM(GGNab)*(Y(Vm)-ENaN)

! Total time independent current
      CURRENTS(Iv)  = CURRENTS(IK1)+CURRENTS(IKp)+CURRENTS(INaK)+
     '  CURRENTS(IpCa)+CURRENTS(ICab)+CURRENTS(INab)

! Compute Cai change at 2 ms after onset of stimulus
      !IF((t.GE.TPS).AND.(t.LT.0.1d0+TPS))THEN
      IF((TIME(TCell).GE.PARAM(TPS)).AND.
     '  (TIME(TCell).LT.PARAM(TABT)+PARAM(TPS)))THEN
        PARAM(Cai_on) = Y(Cab)
        IF(DEBUG) write (*,*) 'Cai_on = ',PARAM(Cai_on)
      ENDIF

      IF((TIME(TCell).GE.2.d0+PARAM(TPS)).AND.
     '  (TIME(TCell).LT.2.d0+PARAM(TABT)+PARAM(TPS)))THEN
        PARAM(Cai2) = Y(Cab)-PARAM(Cai_on)
        IF(DEBUG) write (*,*) 'Cai2 = ',PARAM(Cai2)
      ENDIF

! Ca induced Ca release of JSR
      IF(PARAM(Cai2).GT.PARAM(Caith)) THEN
        GGrel = PARAM(GGrel_)
      ELSE
        GGrel = 0.d0
      ENDIF
      !time of Ca-induced Ca-release
c dpn 03/06/98      t_CICR  = t-TPS
      t_CICR  = TIME(TCell)-(PARAM(TPS)+2.0d0)
      Grel  = GGrel*((PARAM(Cai2)-PARAM(Caith))/
     '  (PARAM(Kmrel)+PARAM(Cai2)-PARAM(Caith)))
     '  *(1.d0-dexp(-t_CICR/PARAM(Tau_on)))*
     '  dexp(-t_CICR/PARAM(Tau_off))
      CURRENTS(Irel) = Grel*(Y(CaJSR)-Y(Cai))

! Ca release of JSR under Ca-overload conditions
!      CSQNth  = 0.7d0
!      IF(CSQN.GE.CSQNth) THEN
!         GGrel = 4.d0                !ms-1
!      ELSE
!         GGrel = 0.d0                !ms-1
!      ENDIF
!      Grel    = GGrel*(1.d0-dexp(-t/Tau_on))*dexp(-t/Tau_off)
!      I(13)   = Grel*(CaJSR-Cai)


! Ca uptake and leakage of NSR Iup and Ileak
      Kleak      = PARAM(IIup)/PARAM(CaNSRCaNSR)        !ms-1
      CURRENTS(Ileak) = Kleak*Y(CaNSR)            !mmol/L/ms
      CURRENTS(Iup)   = PARAM(IIup)*(Y(Cai)/
     '  (Y(Cai)+PARAM(Kmup))) !mmol/L/ms


! Translocation of Ca ions from NSR to JSR Itr
      CURRENTS(Itr)  = (Y(CaNSR)-Y(CaJSR))/PARAM(Tau_tr) !mmol/L/ms

! Computes current carrying the five ions
      CURRENTS(I_Na) = CURRENTS(INa)+CURRENTS(ICaLNa)+3.0d0*
     '  CURRENTS(INaCa)+
     '  3.0d0*CURRENTS(INaK)+CURRENTS(InsNa)+CURRENTS(INab)
      CURRENTS(I_Ca) = CURRENTS(ICaLCa)-CURRENTS(INaCa)+CURRENTS(IpCa)+
     '  CURRENTS(ICab)
      CURRENTS(I_K)  = CURRENTS(ICaLK)+CURRENTS(IK)+CURRENTS(IK1)+
     '  CURRENTS(IKp)-2.0d0*CURRENTS(INaK)+CURRENTS(InsK)
      CURRENTS(F_NSR) = -CURRENTS(Iup)+CURRENTS(Ileak)+CURRENTS(Itr)
      CURRENTS(F_JSR) = CURRENTS(Irel)-CURRENTS(Itr)*VRatio_JN

C *** No error
      ERR_CODE = 0

      RETURN
      END


