      SUBROUTINE JRW_CELL(TIME,Y,DY,CONTROL,MODEL,SIZES,VARIANT,
     '  DERIVED,PARAM,PROTOCOL,AII,AIO,ARI,ARO,ERR_CODE)

C#### Subroutine: JRW_CELL
C###  Description:
C###    Calculates the RHS of the Jafri-Rice-Winslow ionic
C###    current model.

      IMPLICIT NONE

      INCLUDE 'cell_reserved.inc'
      INCLUDE 'cell_jrw.inc'

      INTEGER SIZES(11)
      INTEGER CONTROL(*),MODEL(*),VARIANT,AII(*),AIO(*),ERR_CODE
      REAL*8 TIME(*),Y(*),DY(*),DERIVED(*),PARAM(*),PROTOCOL(*),
     '  ARI(*),ARO(*)

      !local variables
      REAL*8 Istim, !Stimulus current (uA.mm^-2)
     '  INa, !Na+ current (uA.mm^-2)
     '  ICa, !L-type channel Ca++ current (uA.mm^-2)
     '  IK,
     '  IK1, !time independent K+ current (uA.mm^-2)
     '  IKp, !plateau K+ current (uA.mm^-2)
     '  INaCa, !Na+/Ca++ exchanger current (uA.mm^-2)
     '  INaK, !Na+/K+ pump current (uA.mm^-2)
     '  InsCa, !non-specific Ca-activated total current (uA.mm^-2)
     '  IpCa, !sarcolemmal Ca++ pump current (uA.mm^-2)
     '  ICaK, !L-type channel K+ current (uA.mm^-2)
     '  ICab, !Ca++ background current (uA.mm^-2)
     '  INab, !Na+ background current (uA.mm^-2)
     '  InsNa, !non-specific Ca-activated Na+ current (uA.mm^-2)
     '  InsK, !non-specific Ca-activated K+ current (uA.mm^-2)
     '  Jrel, !RyR channel flux (nmol.mm^-3.ms^-1)=(mM.ms^-1)
     '  Jleak, !leak flux from NSR to myoplasm (nmol.mm^-3.ms^-1)
     '  Jup, !SR Ca-ATPase pump flux (nmol.mm^-3.ms^-1)
     '  Jtr, !transfer flux from NSR to JSR (nmol.mm^-3.ms^-1)
     '  Jxfer, !transfer flux from SS to myoplasm (nmol.mm^-3.ms^-1)
     '  Jtrpn !Ca++ binding flux (nmol.mm^-3.ms^-1)

      REAL*8 Bi,BSS,BJSR,ENa,EK,EK1,EKp,ECaN,ENaN,
     '  alfa_m,beta_m,alfa_h,alfa_j,beta_h,beta_j,
     '  alfa_x,beta_x,alfa_K1,beta_K1,alfa_Ca,beta_Ca,
     '  RTF,VRTF,TVRTF,GK,Xi,GK1,K1inf,Kp,sigma,fNaK,InsNa_max,InsK_max,
     '  ICa_max,PK_dash,ICaK_max,y_inf,tau_y,alfa_Ca_dash,beta_Ca_dash,
     '  gama_Ca

C     Physical constants
      REAL*8 R,FC
      PARAMETER(R=8314.472d0) !Gas constant (pJ.nmol^-1.K^-1)
      PARAMETER(FC=96485.341d0) !Faradays constant (nC.nmol^-1)

C *** Check array sizes
      IF(SIZES(NUM_PARAM).LT.117) THEN
        ERR_CODE=1
      ELSE IF(SIZES(NUM_MODEL).LT.1) THEN
        ERR_CODE=2
      ELSE IF(SIZES(NUM_CONTROL).LT.1) THEN
        ERR_CODE=3
      ELSE IF(SIZES(NUM_EQN).LT.25) THEN
        ERR_CODE=4
      ELSE IF(SIZES(NUM_PROTOCOL).LT.10) THEN
        ERR_CODE=5
      ELSE
        IF(CONTROL(RETURN_CURRENTS).NE.0.AND.
     '    SIZES(NUM_DERIVED).LT.37) THEN
          ERR_CODE=6
        ELSE
          ERR_CODE=0

C *** Set stimulus current
          Istim=PROTOCOL(PseudoIs)
          IF(TIME(TCell).GE.PROTOCOL(Is1start).AND.
     '      TIME(TCell).LT.PROTOCOL(Is1stop))
     '      Istim=Istim+PROTOCOL(Is1current)
          IF(TIME(TCell).GE.PROTOCOL(Is2start).AND.
     '      TIME(TCell).LT.PROTOCOL(Is2stop))
     '      Istim=Istim+PROTOCOL(Is2current)
          IF(PROTOCOL(IsFreqPeriod).GT.1.0d-6) THEN
            IF(DMOD(TIME(TCell),PROTOCOL(IsFreqPeriod)).LT.
     '        PROTOCOL(IsFreqDuration))
     '        Istim=Istim+PROTOCOL(IsFreqMag)
          ENDIF
          !convert from volume to area current (uA.mm^-2)
          Istim=Istim/PARAM(Am)

          RTF = R*PARAM(TEMP)/FC !RT/F (mV)
          VRTF = Y(Vm)/RTF !VF/RT (dimensionless)
          TVRTF = 2.d0*VRTF

C ***     Calculate equilibrium potentials
          ENa  = RTF*DLOG(PARAM(Nao)/Y(Nai)) !mV
          EK   = RTF*DLOG((PARAM(Ko)+PARAM(PNaK)*
     '      PARAM(Nao))/(Y(Ki)+PARAM(PNaK)*Y(Nai))) !mV
          EK1  = RTF*DLOG(PARAM(Ko)/Y(Ki)) !mV
          EKp  = EK1 !mV
          ECaN = 0.5d0*RTF*DLOG(Y(Cao)/Y(Cai)) !Ca background (mV)
          ENaN = ENa !Na background (mV)

C ***     Calculate buffer scale factors
          !Myoplasmic buffer (calmodulin)
          Bi   = 1.d0/(1.d0+PARAM(CMDN_tot)*PARAM(km_CMDN)/
     '      (PARAM(km_CMDN)+Y(Cai))**2) !(dimensionless)
          !Sub-space buffer (calmodulin)
          BSS  = 1.d0/(1.d0+PARAM(CMDN_tot)*PARAM(km_CMDN)/
     '      (PARAM(km_CMDN)+Y(CaSS))**2) !(dimensionless)
          !JSR buffer (calsequestrin)
          BJSR = 1.d0/(1.d0+PARAM(CSQN_tot)*PARAM(km_CSQN)/
     '      (PARAM(km_CSQN)+Y(CaJSR))**2) !(dimensionless)

C ***     Calculate Ca++ fluxes
          !RyR channel flux
          Jrel  = PARAM(v1)*(Y(PO1)+Y(PO2))*(Y(CaJSR)-
     '      Y(CaSS)) !(nmol.L^-1.ms^-1)
          !leak flux from NSR to myoplasm
          Jleak = PARAM(v2)*(Y(CaNSR)-Y(Cai)) !(mM.ms^-1)
          !SR Ca-ATPase pump flux
          Jup   = PARAM(v3)*Y(Cai)**2/(PARAM(km_up)**2+
     '      Y(Cai)**2) !(mM.ms^-1)
          !transfer flux from NSR to JSR
          Jtr   = (Y(CaNSR)-Y(CaJSR))/PARAM(tau_tr) !(mM.ms^-1)
          !transfer flux from SS to myoplasm
          Jxfer = (Y(CaSS)-Y(Cai))/PARAM(tau_xfer) !(mM.ms^-1)

C ***     Ca++ binding to low-affinity tronponin
          DY(LTRPNCa) = PARAM(kp_ltrpn)*Y(Cai)*(
     '      PARAM(LTRPN_tot)-Y(LTRPNCa))-PARAM(km_ltrpn)*
     '      Y(LTRPNCa) !(mM.ms^-1)
          !Ca++ binding flux
          Jtrpn = PARAM(kp_htrpn)*Y(Cai)*(PARAM(HTRPN_tot)-
     '      Y(HTRPNCa)) - PARAM(km_htrpn)*Y(HTRPNCa) +
     '      DY(LTRPNCa) !(mM.ms^-1)

C ***     Calculate the gating rate constants
          !INa channel
          alfa_m = 0.32d0*(Y(Vm)+47.13d0)/(1.d0-DEXP(-0.1d0*(Y(Vm)+
     '      47.13d0)))
          beta_m = 0.08d0*DEXP(-Y(Vm)/11.d0)
          IF(Y(Vm).GE.-40.d0) THEN
            alfa_h = 0.d0
            beta_h = 1.d0/(0.13d0*(1.d0+DEXP(-(Y(Vm)+10.66d0)/11.1d0)))
            alfa_j = 0.d0
            beta_j = 0.3d0*DEXP(-2.535d-7*Y(Vm))
     '        /(1.d0+DEXP(-0.1d0*(Y(Vm)+32.d0)))
          ELSE
            alfa_h = 0.135d0*DEXP(-(80.d0+Y(Vm))/6.8d0)
            beta_h = 3.56d0*DEXP(0.079d0*Y(Vm))+3.1d5*DEXP(0.35d0*Y(Vm))
            alfa_j = (-1.27140d5*DEXP(0.2444d0*Y(Vm))
     '        -3.474d-5*DEXP(-0.04391d0*Y(Vm)))
     '        *(Y(Vm)+37.78d0)/(1.d0+DEXP(0.311d0*(Y(Vm)+79.23d0)))
            beta_j = 0.1212d0*DEXP(-0.01052d0*Y(Vm))
     '        /(1.d0+DEXP(-0.1378d0*(Y(Vm)+40.14d0)))
          ENDIF !Y(Vm)

          !IK channel
          alfa_x = 7.19d-5*(Y(Vm)+30.d0)/( 1.d0-DEXP(-0.148d0*(Y(Vm)+
     '      30.d0)))
          beta_x = 1.31d-4*(Y(Vm)+30.d0)/(-1.d0+DEXP(0.0687d0*(Y(Vm)+
     '      30.d0)))

          !IK1 channel
          alfa_k1 = 1.02d0/(1.d0+DEXP(0.2385d0*(Y(Vm)-EK1-59.215d0)))
          beta_k1 = (0.4912d0*DEXP(0.08032d0*(Y(Vm)-EK1+5.476d0))
     '      +DEXP(0.06175d0*(Y(Vm)-EK1-594.31d0)))
     '      /(1.d0+DEXP(-0.5143d0*(Y(Vm)-EK1+4.753d0)))

          !ICaL channel
          alfa_Ca = 0.4d0*DEXP((Y(Vm)+12.d0)/10.d0)
          beta_Ca = 0.05d0*DEXP(-(Y(Vm)+12.d0)/13.d0)

C ***     Ionic currents
          ! fast Na+ current - INa
          INa = PARAM(GNa_max)*Y(mNa)*Y(mNa)*Y(mNa)*Y(hNa)*
     '      Y(jNa)*(Y(Vm)-ENa)
          ! K+ current - IK
          GK = 0.1128d0*DSQRT(PARAM(Ko)/5.4d0)
          Xi = 1.d0/(1.d0+DEXP((Y(Vm)-56.26d0)/32.1d0))
          IK = GK*Xi*Y(xK)*Y(xK)*(Y(Vm)-EK)
          ! time-independent K+ current - IK1
          GK1 = 0.75d0*DSQRT(PARAM(Ko)/5.4d0)
          K1inf = alfa_K1/(alfa_K1+beta_K1)
          IK1 = GK1*K1inf*(Y(Vm)-EK1)
          ! plateau K+ current - IKp
          Kp = 1.d0/(1.d0+DEXP((7.488d0-Y(Vm))/5.98d0))
          IKp = PARAM(GKp_max)*Kp*(Y(Vm)-EKp)
          ! Na+ Ca++ exchanger current - INaCa
          INaCa = PARAM(kNaCa)*(DEXP(PARAM(eta)*VRTF)*Y(Nai)*
     '      Y(Nai)*Y(Nai)*Y(Cao)-DEXP((PARAM(eta)-1.d0)*VRTF)*
     '      PARAM(Nao)*PARAM(Nao)*PARAM(Nao)*Y(Cai))/
     '      ((PARAM(km_Na)**3+PARAM(Nao)**3)*
     '      (PARAM(km_Ca)+Y(Cao))*(1.d0+PARAM(ksat)*
     '      DEXP((PARAM(eta)-1.d0)*VRTF)))
          ! Na+ K+ pump current - INaK
          sigma = (DEXP(PARAM(Nao)/67.3d0)-1.d0)/7.d0
          fNaK = 1.d0/(1.d0+0.1245d0*DEXP(-0.1d0*VRTF)
     '      +0.0365d0*sigma*DEXP(-VRTF))
          INaK = PARAM(INaK_max)*fNaK*PARAM(Ko)/
     '      (PARAM(Ko)+PARAM(km_Ko))/(1.d0+(
     '      PARAM(km_Nai)/Y(Nai))**1.5d0)
          ! Nonspecific Ca++ activated current - InsCa
          InsNa_max = PARAM(P_ns)*VRTF*FC*0.75d0*(Y(Nai)*
     '      DEXP(VRTF)-PARAM(Nao))/(DEXP(VRTF)-1.d0)
          InsNa = InsNa_max/(1.d0+(PARAM(km_ns)/Y(Cai))**3)
          InsK_max = PARAM(P_ns)*VRTF*FC*0.75d0*(Y(Ki)*DEXP(VRTF)-
     '      PARAM(Ko))/(DEXP(VRTF)-1.d0)
          InsK = InsK_max/(1.d0+(PARAM(km_ns)/Y(Cai))**3)
          InsCa = InsNa+InsK
          ! Sarcolemmal Ca++ pump current - IpCa
          IpCa = PARAM(ICap_max)*Y(Cai)/
     '      (PARAM(km_Cap)+Y(Cai))
          ! Ca++ background current - ICab
          ICab = PARAM(GCab_max)*(Y(Vm)-ECaN)
          ! Na+ background current - INab
          INab = PARAM(GNab_max)*(Y(Vm)-ENaN)
          ! L-type Ca++ currents
          ICa_max = PARAM(PCa)*4.0d0*VRTF*FC*(0.001d0*DEXP(TVRTF)-
     '      0.341d0*Y(Cao))/(DEXP(TVRTF)-1.0d0)
          ICa = ICa_max*Y(yL)*(Y(O)+Y(OCa))
          PK_dash = PARAM(PK)/(1.d0+ICa_max/PARAM(ICa_half))
          ICaK_max = PK_dash*VRTF*FC*(Y(Ki)*DEXP(VRTF)-PARAM(Ko))/
     '      (DEXP(VRTF)-1.d0)
          ICaK = ICaK_max*Y(yL)*(Y(O)+Y(OCa))

C ***     Calculate the RHS of the O.D.E.'s
          ! Membrane potential - d(Vm)/dt
          DY(Vm) = -(Istim+INa+ICa+IK+IK1+IKp+INaCa+INaK+InsCa+IpCa+
     '      ICaK+ICab+INab)/PARAM(Cm)
          ! INa activation gate - d(mNa)/dt
          DY(mNa) = alfa_m*(1.d0-Y(mNa))-beta_m*Y(mNa)
          ! INa inactivation gate - d(hNa)/dt
          DY(hNa) = alfa_h*(1.d0-Y(hNa))-beta_h*Y(hNa)
          ! INa slow inactivation gate - d(jNa)/dt
          DY(jNa) = alfa_j*(1.d0-Y(jNa))-beta_j*Y(jNa)
          ! IK inactivation gate - d(xK)/dt
          DY(xK) = alfa_x*(1.d0-Y(xK))-beta_x*Y(xK)
          ! Intracellular [Na+] - d(Nai)/dt
          DY(Nai) = -(INa+INab+InsNa+3.d0*INaCa+3.d0*INaK)*
     '      PARAM(Acap)/(PARAM(V_myo)*FC)
          ! Intracellular [K+] - d(Ki)/dt
          DY(Ki) = -(IK+IK1+IKp+InsK+ICaK-2.d0*INaK)*PARAM(Acap)/
     '      (PARAM(V_myo)*FC)
          ! Intracellular [Ca++] - d(Cai)/dt
          DY(Cai) = Bi*(Jleak+Jxfer-Jup-Jtrpn-(ICab-2.d0*INaCa+
     '      IpCa)*PARAM(Acap)/(2.d0*PARAM(V_myo)*FC))
          ! Sub-space [Ca++] - d(CaSS)/dt
          DY(CaSS) = BSS*(Jrel*PARAM(V_JSR)/PARAM(V_SS)-
     '      Jxfer*PARAM(V_myo)/PARAM(V_SS)-
     '      ICa*PARAM(Acap)/(2.d0*PARAM(V_SS)*FC))
          ! JSR [Ca++] - d(CaJSR)/dt
          DY(CaJSR) = BJSR*(Jtr-Jrel)
          ! NSR [Ca++] - d(CaNSR)/dt
          DY(CaNSR) = (Jup-Jleak)*PARAM(V_myo)/PARAM(V_NSR)
          ! [Ca++] bound to high-affinity troponin - d(HTRPNCa)/dt
          DY(HTRPNCa) = PARAM(kp_htrpn)*Y(Cai)*
     '      (PARAM(HTRPN_tot)-Y(HTRPNCa))-PARAM(km_htrpn)*
     '      Y(HTRPNCa)
          ! [Ca++] bound to low-affinity troponin - d(LTRPNCa)/dt
c         moved up to where Jtrpn is evaluated
          ! Extracellular [Ca++] - d(Cao)/dt
          DY(Cao) = 0.0d0
          ! Fraction of RyR's in closed state 1 - d(PC1)/dt
          DY(PC1) = -PARAM(kap)*Y(CaSS)**PARAM(n)*Y(PC1)+
     '      PARAM(kam)*Y(PO1)
          ! Fraction of RyR's in open state 1 - d(PO1)/dt
          DY(PO1) = PARAM(kap)*Y(CaSS)**PARAM(n)*Y(PC1)-
     '      PARAM(kam)*Y(PO1)-PARAM(kbp)*
     '      Y(CaSS)**PARAM(m)*Y(PO1)+PARAM(kbm)*Y(PO2)-
     '      PARAM(kcp)*Y(PO1)+PARAM(kcm)*Y(PC2)
          ! Fraction of RyR's in open state 2 - d(PO2)/dt
          DY(PO2) = PARAM(kbp)*DY(CaSS)**PARAM(m)*Y(PO1)-
     '      PARAM(kbm)*Y(PO2)
          ! Fraction of RyR's in closed state 2 - d(PC2)/dt
          DY(PC2) = PARAM(kcp)*Y(PO1)-PARAM(kcm)*Y(PC2)
          ! Voltage-dependent, ICaL inactivation gate - d(yL)/dt
          y_inf = 1.0d0/(1.d0+DEXP((Y(Vm)+55.d0)/7.5d0))+
     '      0.1d0/(1.d0+DEXP((-Y(Vm)+21.d0)/6.d0))
          tau_y = 20.d0+600.d0/(1.d0+DEXP((Y(Vm)+30.d0)/9.5d0))
          DY(yL) = (y_inf-Y(yL))/tau_y
          ! Fraction of ICaL channels in closed state 0, mode normal
          alfa_Ca_dash = alfa_Ca*PARAM(a)
          beta_Ca_dash = beta_Ca/PARAM(b)
          gama_Ca = 0.1875d0*Y(CaSS)
          DY(C0) = beta_Ca*Y(C1)+PARAM(w)*Y(Cca0)-(4.d0*alfa_Ca+
     '      gama_Ca)*Y(C0)
          ! Fraction of ICaL channels in closed state 1, mode normal
          DY(C1) = 4.d0*alfa_Ca*Y(C0)+2.d0*beta_Ca*Y(C2)+PARAM(w)/
     '      PARAM(b)*Y(Cca1)-(beta_Ca+3.d0*alfa_Ca+gama_Ca*
     '      PARAM(a))*Y(C1)
          ! Fraction of ICaL channels in closed state 2, mode normal
          DY(C2) = 3.d0*alfa_Ca*Y(C1)+3.d0*beta_Ca*Y(C3)+PARAM(w)/
     '      PARAM(b)**2*Y(Cca2)-(2.d0*beta_Ca+2.d0*alfa_Ca+
     '      gama_Ca*PARAM(a)**2)*Y(C2)
          ! Fraction of ICaL channels in closed state 3, mode normal
          DY(C3) = 2.d0*alfa_Ca*Y(C2)+4.d0*beta_Ca*Y(C4)+PARAM(w)/
     '      PARAM(b)**3*Y(Cca3)-(3.d0*beta_Ca+alfa_Ca+gama_Ca*
     '      PARAM(a)**3)*Y(C3)
          ! Fraction of ICaL channels in closed state 4, mode normal
          DY(C4) = alfa_Ca*Y(C3)+PARAM(g)*Y(O)+PARAM(w)/
     '      PARAM(b)**4*Y(Cca4)-(4.d0*beta_Ca+PARAM(f)+
     '      gama_Ca*PARAM(a)**4)*Y(C4)
          ! Fraction of ICaL channels in the open state, mode normal
          DY(O) = PARAM(f)*Y(C4)-PARAM(g)*Y(O)
          ! Fraction of ICaL channels in closed state 0, mode Ca
          DY(Cca0) = beta_Ca_dash*Y(Cca1)+gama_Ca*Y(C0)-(4.d0*
     '      alfa_Ca_dash+PARAM(w))*Y(Cca0)
          ! Fraction of ICaL channels in closed state 1, mode Ca
          DY(Cca1) = 4.d0*alfa_Ca_dash*Y(Cca0)+2.d0*beta_Ca_dash*
     '      Y(Cca2)+gama_Ca*PARAM(a)*Y(C1)-(beta_Ca_dash+3.d0*
     '      alfa_Ca_dash+PARAM(w)/PARAM(b))*Y(Cca1)
          ! Fraction of ICaL channels in closed state 2, mode Ca
          DY(Cca2) = 3.d0*alfa_Ca_dash*Y(Cca1)+3.d0*beta_Ca_dash*
     '      Y(Cca3)+gama_Ca*PARAM(a)**2*Y(C2)-(2.d0*beta_Ca_dash+
     '      2.d0*alfa_Ca_dash+PARAM(w)/PARAM(b)**2)*Y(Cca2)
          ! Fraction of ICaL channels in closed state 3, mode Ca
          DY(Cca3) = 2.d0*alfa_Ca_dash*Y(Cca2)+4.d0*beta_Ca_dash*
     '      Y(Cca4)+gama_Ca*PARAM(a)**3*Y(C3)-(3.d0*beta_Ca_dash+
     '      alfa_Ca_dash+PARAM(w)/PARAM(b)**3)*Y(Cca3)
          ! Fraction of ICaL channels in closed state 4, mode Ca
          DY(Cca4) = alfa_Ca_dash*Y(Cca3)+PARAM(g_dash)*Y(Oca)+
     '      gama_Ca*PARAM(a)**4*Y(C4)-(4.d0*beta_Ca_dash+
     '      PARAM(f_dash)+PARAM(w)/PARAM(b)**4)*Y(Cca4)
          ! Fraction of ICaL channels in the open state, mode Ca
          DY(Oca) = PARAM(f_dash)*Y(Cca4)-PARAM(g_dash)*Y(Oca)



          IF(CONTROL(RETURN_CURRENTS).NE.0) THEN
            DERIVED(DIstim)=Istim
            DERIVED(DINa)=INa
            DERIVED(DICa)=ICa
            DERIVED(DIK)=IK
            DERIVED(DIK1)=IK1
            DERIVED(DIKp)=IKp
            DERIVED(DINaCa)=INaCa
            DERIVED(DINaK)=INaK
            DERIVED(DInsCa)=InsCa
            DERIVED(DIpCa)=IpCa
            DERIVED(DICaK)=ICaK
            DERIVED(DICab)=ICab
            DERIVED(DINab)=INab
            DERIVED(DInsNa)=InsNa
            DERIVED(DInsK)=InsK
            DERIVED(DJrel)=Jrel
            DERIVED(DJleak)=Jleak
            DERIVED(DJup)=Jup
            DERIVED(DJtr)=Jtr
            DERIVED(DJxfer)=Jxfer
            DERIVED(DJtrpn)=Jtrpn
          ENDIF
        ENDIF
      ENDIF
      RETURN
      END


