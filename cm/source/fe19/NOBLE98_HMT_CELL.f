      SUBROUTINE NOBLE98_HMT_CELL(TIME,Y,DY,CONTROL,MODEL,SIZES,VARIANT,
     '  DERIVED,PARAM,PROTOCOL,AII,AIO,ARI,ARO,ERR_CODE)

C#### Subroutine: NOBLE98_HMT_CELL
C###  Description:
C###    Calculates the RHS for a coupled Noble 98 and HMT model

      IMPLICIT NONE

      INCLUDE 'cell_reserved.inc'
      INCLUDE 'cell_n98_hmt.inc'
      INCLUDE 'tol00.cmn'

      INTEGER SIZES(11)
      INTEGER CONTROL(*),MODEL(*),VARIANT,AII(*),AIO(*),ERR_CODE
      REAL*8 TIME(*),Y(*),DY(*),DERIVED(*),PARAM(*),PROTOCOL(*),
     '  ARI(*),ARO(*)

      !local variables
C *** Noble '98
      REAL*8 IK1, !Time-independent (background) K+ current (uA.mm^-2)
     '  Ito, !Transient outward K+ current (uA.mm^-2)
     '  IKr1, !Time-dependent (delayed) fast K+ current (uA.mm^-2)
     '  IKr2, !Time-dependent (delayed) fast K+ current (uA.mm^-2)
     '  IKs, !Time-dependent (delayed) slow K+ current (uA.mm^-2)
     '  IKNa, !Na+-dependent K+ current (uA.mm^-2)
     '  IbK, !Background K+ current (uA.mm^-2)
     '  IKATP, !ATP-dependent K+ current (uA.mm^-2)
     '  IKACh,  !ACh-dependent K+ current (uA.mm^-2)
     '  INa, !Fast Na+ current (uA.mm^-2)
     '  IbNa, !Background Na+ current (uA.mm^-2)
     '  IpNa, !Persistent Na+ current (uA.mm^-2)
     '  ICaLK, !L-type K+ current (uA.mm^-2)
     '  ICaLNa, !L-type Na+ current (uA.mm^-2)
     '  ICaLCa, !L-type Ca++ current (uA.mm^-2)
     '  ICaLKDS, !Diadic Space L-type K+ current (uA.mm^-2)
     '  ICaLNaDS, !Diadic Space L-type Na+ current (uA.mm^-2)
     '  ICaLCaDS, !Diadic Space L-type Ca++ current (uA.mm^-2)
     '  IbCa, !Background Ca++ current (uA.mm^-2)
     '  INaK, !Na+/K+ exchange current (uA.mm^-2)
     '  INaCa, !Na+/Ca++ exchange current (uA.mm^-2)
     '  INaCaDS, !Na+/Ca++ Diadic space exchange current (uA.mm^-2)
     '  IK_stretch, !K+ stretch activated current (uA.mm^-2)
     '  INa_stretch, !Na+ stretch activated current (uA.mm^-2)
     '  ICa_stretch, !Ca++ stretch activated current (uA.mm^-2)
     '  INs_stretch, !Ns stretch activated current (uA.mm^-2)
     '  IAn_stretch, !An stretch activated current (uA.mm^-2)
     '  Istim, !Stimulus current (uA.mm^-2)
     '  IK_tot, !Total K+ current (uA.mm^-2)
     '  INa_tot, !Total Na+ current (uA.mm^-2)
     '  ICa_tot, !Total Ca+ current (uA.mm^-2)
     '  ICaDS_tot !Total DS Ca+ current (uA.mm^-2)
      REAL*8 jup, !SR Ca++ uptake flux (mmol.L^-1.ms^-1)
     '  jtr, !SR Ca++ translocation flux (mmol.L^-1.ms^-1)
     '  jrel, !SR Ca++ release flux (mmol.L^-1.ms^-1)
     '  jleak, !SR Ca++ leakage flux (mmol.L^-1.ms^-1)
     '  jdecay !DS Ca++ decay flux (mmol.L^-1.ms^-1)
      REAL*8 RTF, !RT/F (mV)
     '  ENa, !Na reversal potential (mV)
     '  EK, !K reversal potential (mV)
     '  EKs, !K (IKs) reversal potential (mV)
     '  ECa, !Ca reversal potential (mV)
     '  Emh !Na (INa) reversal potential (mV)
      REAL*8 alpha_xr1, !IKr x1 alpha rate (ms^-1)
     '  beta_xr1, !IKr x1 beta rate (ms^-1)
     '  alpha_xr2, !IKr x2 alpha rate (ms^-1)
     '  beta_xr2, !IKr x2 beta rate (ms^-1)
     '  alpha_xs, !IKs x alpha rate (ms^-1)
     '  beta_xs, !IKs x beta rate (ms^-1)
     '  alpha_sto, !Ito s alpha rate (ms^-1)
     '  beta_sto, !Ito s beta rate (ms^-1)
     '  alpha_mNa, !INa m alpha rate (ms^-1)
     '  beta_mNa, !INa m beta rate (ms^-1)
     '  alpha_hNa, !INa h alpha rate (ms^-1)
     '  beta_hNa, !INa h beta rate (ms^-1)
     '  alpha_dCa, !ILCa d alpha rate (ms^-1)
     '  beta_dCa, !ILCa d beta rate (ms^-1)
     '  alpha_fCa, !ILCa f alpha rate (ms^-1)
     '  beta_fCa, !ILCa f beta rate (ms^-1)
     '  alpha_xACh1, !IKACh x1 alpha rate (ms^-1)
     '  beta_xACh1, !IKACh x1 beta rate (ms^-1)
     '  alpha_xACh2, !IKACh x2 alpha rate (ms^-1)
     '  beta_xACh2, !IKACh x2 beta rate (ms^-1)
     '  alpha_act, !ICaL activation rate (ms^-1)
     '  alpha_inact !ICaL inactivation rate (ms^-1)
      REAL*8 c1,c2,FractBackSRSites,FractCaUptakeSites,f_stretch,
     '  OpenReleaseChannelFract,PrecursorFract,RegulatoryBindingSite,
     '  VRTF,z1,z2

C *** HMT
      REAL*8 Cb_norm,C50,dPHI1,dPHI2,dPHI3,LENGTH_DIFF,Q,
     '  TTo,Tm_n,Tm_p50

C***  pH dependence
      REAL*8 pHi_IK1, !multiplier for IK1 conductance variation with pHi ()
     '  pHi_ICaL,  !multiplier for ICaL conductance variation with pHi ()
     '  pHi_INaCa !multiplier for INaCa conductance variation with pHi ()


C     Physical constants
      REAL*8 R,F,nNaK,nNaCa
c      PARAMETER(R=8314.41d0) !Gas constant (mJ.mol^-1.K^-1)
c      PARAMETER(F=96485.0d0) !Faradays constant (C.mol^-1)
      PARAMETER(R=8314.472d0) !Gas constant (pJ.nmol^-1.K^-1)
      PARAMETER(F=96485.341d0) !Faradays constant (nC.nmol^-1)
      PARAMETER(nNaK=1.5d0) !Na+/K+ exchange stoichometry ratio ()
      PARAMETER(nNaCa=3.0d0) !Na+/Ca++ exchange stoichometry ratio ()

      IF(SIZES(NUM_PARAM).LT.141) THEN
        ERR_CODE=1
      ELSE IF(SIZES(NUM_MODEL).LT.1) THEN
        ERR_CODE=2
      ELSE IF(SIZES(NUM_CONTROL).LT.1) THEN
        ERR_CODE=3
      ELSE IF(SIZES(NUM_EQN).LT.29) THEN
        ERR_CODE=4
      ELSE IF(SIZES(NUM_PROTOCOL).LT.10) THEN
        ERR_CODE=5
      ELSE
        IF(CONTROL(RETURN_CURRENT).NE.0.AND.SIZES(NUM_DERIVED).LT.37)
     '    THEN
          ERR_CODE=6
        ELSE
          ERR_CODE=0

C *** DPN 27 September 2000 - need to check that the extension ratio has
C         not exceeded the physiological limit, and that the
C         proportion of actin sites does not fall outside the range 0 <=
C         z < 1. But still solve the model even if the values are out of
C         range, simply set them to the limits - this is to allow
C         distributed models to exceed the limits while trying to
C         converge.
C *** ?? Up to the user to check if the converged solution still exceeds
C         the limits ??
          IF ((Y(ExtensionRatio)-PARAM(MaxExtRatio)).GT.ZERO_TOL) THEN
            WRITE(*,
     '        '(''WARNING: Extension ratio above upper limit '',E12.5)')
     '        Y(ExtensionRatio)
            Y(ExtensionRatio) = PARAM(MaxExtRatio)
          ELSEIF ((PARAM(MinExtRatio)-Y(ExtensionRatio)).GT.ZERO_TOL)
     '        THEN
            WRITE(*,
     '        '(''WARNING: Extension ratio below lower limit '',E12.5)')
     '        Y(ExtensionRatio)
            Y(ExtensionRatio) = PARAM(MinExtRatio)
          ENDIF
          IF (Y(z).LT.0.0d0) THEN
            WRITE(*,'(''WARNING: Z below lower limit '',E12.5)') Y(z)
            Y(z) = 0.0d0
          ELSEIF (Y(z).GE.1.0d0) THEN
            WRITE(*,'(''WARNING: Z above upper limit '',E12.5)') Y(z)
            Y(z) = 1.0d0
          ENDIF
          IF(CONTROL(ODE).EQ.1) THEN !evaluating DY for ODE's
            Istim=PROTOCOL(PseudoIs)
            IF(TIME(TCell).GE.PROTOCOL(Is1start).AND.
     '        TIME(TCell).LT.PROTOCOL(Is1stop))
     '        Istim=Istim+PROTOCOL(Is1current)
            IF(TIME(TCell).GE.PROTOCOL(Is2start).AND.
     '        TIME(TCell).LT.PROTOCOL(Is2stop))
     '        Istim=Istim+PROTOCOL(Is2current)
            IF(PROTOCOL(IsFreqPeriod).GT.1.0d-6) THEN
              IF(DMOD(TIME(TCell),PROTOCOL(IsFreqPeriod)).LT.
     '          PROTOCOL(IsFreqDuration))
     '          Istim=Istim+PROTOCOL(IsFreqMag)
            ENDIF
            !convert from volume to area current (uA.mm^-2)
            Istim=Istim/PARAM(Am)

            RTF = R*PARAM(TEMP)/F

C           Equilibrium potentials
            ENa = RTF*DLOG(PARAM(Nao)/Y(Nai)) !(mV) Equation 1
            EK = RTF*DLOG(PARAM(Ko)/Y(Ki)) !(mV) Equation 2
            EKs = RTF*DLOG((PARAM(Ko)+PARAM(PKNa)*
     '        PARAM(Nao))/(Y(Ki)+PARAM(PKNa)*Y(Nai))) !(mV) Equation 3
            ECa = 0.5d0*RTF*DLOG(PARAM(Cao)/Y(Cai)) !(mV) Equation 4
            Emh = RTF*DLOG((PARAM(Nao)+PARAM(PNaK)*
     '        PARAM(Ko))/(Y(Nai)+PARAM(PNaK)*Y(Ki))) !(mV) Equation 5

C pHi dependant modifiers of channel conductance

          IF(MODEL(USE_pHi).EQ.1) THEN
            pHi_IK1=PARAM(KnpHiK1)/(1+(Y(Hi)/
     '        PARAM(KdpHiK1))**2.52d0)
            pHi_ICaL=PARAM(KnpHiCaL)/(1+(Y(Hi)/
     '        PARAM(KdpHiCaL)))
            pHi_INaCa=PARAM(KnpHiNaCa)/(1+(Y(Hi)/
     '        PARAM(KdpHiNaCa))**0.8d0)
          ELSE
            pHi_IK1=1.0d0
            pHi_ICaL=1.0d0
            pHi_INaCa=1.0d0
          ENDIF


C           Time-independent (background) K+ current, IK1
            IK1 = pHi_IK1*PARAM(gK1)*PARAM(Ko)/(PARAM(Ko)+
     '        PARAM(KmK1))*(Y(Vm)-EK)/(1.0d0+
     '        DEXP(PARAM(STEEP_IK1)*(Y(Vm)-EK+10.0d0-
     '        PARAM(SHIFT_IK1))/RTF)) !(uA.mm^-2) Equation 6

C           Transient outward current, Ito
            alpha_sto = PARAM(epsilon_sto)*0.000033d0*
     '        DEXP(-(Y(Vm)-PARAM(SHIFT_sto))/
     '        (2.125d0*PARAM(STEEP_sto))) !(ms^-1) Equation 10
            beta_sto =  PARAM(epsilon_sto)*0.033d0/
     '        (1.0d0+DEXP(-(Y(Vm)+10.0d0-PARAM(SHIFT_sto))/
     '        PARAM(STEEP_sto))) !(ms^-1) Equation 11
            Ito = PARAM(gto)*(PARAM(gtos)+Y(sto)*
     '        (1.0d0-PARAM(gtos)))*Y(rto)*(Y(Vm)-EK) !(uA.mm^-2) Equation 7
            DY(rto) = PARAM(epsilon_rto)*PARAM(alpha_rto)*
     '        (1.0d0/(1.0d0+DEXP(-(Y(Vm)+4.0d0-PARAM(SHIFT_rto))/
     '        PARAM(STEEP_rto)))-Y(rto)) !(ms^-1) Equation 8
            DY(sto) = alpha_sto*(1.0d0-Y(sto))-beta_sto*Y(sto) !(ms^-1) Equation 9

C           Time-dependent (delayed) rapid K+ current, IKr1
            alpha_xr1 = PARAM(epsilon_alpha_xr1)*0.05d0/
     '        (1.0d0+DEXP(-(Y(Vm)-5.0d0-PARAM(SHIFT_alpha_xr1))/
     '        9.0d0)) !(ms^-1) Equation 14
            beta_xr1 = PARAM(epsilon_beta_xr1)*0.00005d0*
     '        DEXP(-(Y(Vm)-20.0d0-PARAM(SHIFT_beta_xr1))/
     '        15.0d0) !(ms^-1) Equation 15
            IKr1 = PARAM(gKr1)*Y(xr1)*(Y(Vm)-EK)/(1.0d0+
     '        DEXP((Y(Vm)+9.0d0-PARAM(SHIFT_IKr))/
     '        22.4d0)) !(uA.mm^-2) Equation 12
            DY(xr1) = alpha_xr1*(1.0d0-Y(xr1))-beta_xr1*Y(xr1) !(ms^-1) Equation 13

C           Time-dependent (delayed) rapid K+ current, IKr2
            alpha_xr2 = PARAM(epsilon_alpha_xr2)*0.05d0/
     '        (1.0d0+DEXP(-(Y(Vm)-5.0d0-PARAM(SHIFT_alpha_xr2))/
     '        9.0d0)) !(ms^-1) Equation 18
            beta_xr2 = PARAM(epsilon_beta_xr2)*0.0004d0*
     '        DEXP(-((Y(Vm)+30.0d0-PARAM(SHIFT_beta_xr2))/
     '        30.0d0)**3) !(ms^-1) Equation 19
            IKr2 = PARAM(gKr2)*Y(xr2)*(Y(Vm)-EK)/(1.0d0+
     '        DEXP((Y(Vm)+9.0d0+PARAM(SHIFT_IKr))/
     '        22.4d0)) !(uA.mm^-2) Equation 16
            DY(xr2) = alpha_xr2*(1.0d0-Y(xr2))-beta_xr2*Y(xr2) !(ms^-1) Equation 17

C           Time-dependent (delayed) slow K+ current, IKs
            alpha_xs = PARAM(epsilon_alpha_xs)*0.014d0/
     '        (1.0d0+DEXP(-(Y(Vm)-40.0d0-PARAM(SHIFT_alpha_xs))/
     '        9.0d0)) !(ms^-1) Equation 22
            beta_xs = PARAM(epsilon_beta_xs)*0.001d0*
     '        DEXP(-(Y(Vm)-PARAM(SHIFT_beta_xs))/45.0d0) !(ms^-1) Equation 23
            IKs = PARAM(gKs)*Y(xs)**2*(Y(Vm)-EKs) !(uA.mm^-2) Equation 20
            DY(xs) = alpha_xs*(1.0d0-Y(xs))-beta_xs*Y(xs) !(ms^-1) Equation 21

C           Na+ dependent K+ current, IKNa
            IKNa = PARAM(gKNa)*Y(Nai)/(Y(Nai)+PARAM(KmgKNa))*
     '        (Y(Vm)-EK) !(uA.mm^-2) Equation 24

C           Background Potassium Current, IbK
            IbK = PARAM(gbK)*(Y(Vm)-EK) !(uA.mm^-2) Equation 25

C           ATP-dependent K+ current, IKATP
            IKATP = PARAM(gKATP)*(Y(Vm)-PARAM(EKATPrev))/
     '        (1.0d0+(PARAM(ATP)/PARAM(KmATP))**2) !(uA.mm^-2) Equation 26

C           ACh-dependent K+ current (Mark Boyett Formulation), IKACh
            alpha_xACh1 = 0.003684211d0 !(ms^-1) Equation 30
            beta_xACh1 = 0.00582d0/(1.0d0+DEXP(-(Y(Vm)+
     '        PARAM(SHIFT_ACh1))/15.0d0)) !(ms^-1) Equation 31
            alpha_xACh2 = 0.07309924d0 !(ms^-1) Equation 32
            beta_xACh2 = 0.12d0/(1.0d0+DEXP(-(Y(Vm)+
     '        PARAM(SHIFT_ACh2))/15.0d0)) !(ms^-1) Equation 33
            IKACh = PARAM(gKACh)*(PARAM(ACh)**1.4969d0/
     '        (PARAM(ACh)**1.4969d0+PARAM(KmACh)**1.4969d0))*
     '        (PARAM(Ko)/(PARAM(Ko)+PARAM(KmKACh)))*
     '        Y(xACh1)*Y(xACh2)*(Y(Vm)-EK)/(1.0d0+DEXP(0.4d0*(Y(Vm)-
     '        EK-140.0d0)/RTF)) !(uA.mm^-2) Equation 27
            DY(xACh1) = alpha_xACh1*(1.0d0-Y(xACh1))-
     '        beta_xACh1*Y(xACh1) !(ms^-1) Equation 28
            DY(xACh2) = alpha_xACh2*(1.0d0-Y(xACh2))-
     '        beta_xACh2*Y(xACh2) !(ms^-1) Equation 29

C           Fast Na+ current, INa
            IF(DABS(Y(Vm)+41.0d0-PARAM(SHIFT_mNa)).LT.0.00001d0) THEN
              alpha_mNa = 2.0d0 !(ms^-1) Equation 37
            ELSE
              alpha_mNa = 0.2d0*(Y(Vm)+41.d0-PARAM(SHIFT_mNa))/
     '          (1.0d0-DEXP(-(Y(Vm)+41.0d0-PARAM(SHIFT_mNa))/
     '          10.0d0)) !(ms^-1) Equation 37
            ENDIF
            beta_mNa = 8.0d0*DEXP(-(Y(Vm)+66.d0-PARAM(SHIFT_mNa))/
     '        18.0d0) !(ms^-1) Equation 38
            alpha_hNa = 0.02d0*DEXP(-(Y(Vm)+75.d0-
     '        PARAM(SHIFT_hNa))/8.0d0) !(ms^-1) Equation 39
            beta_hNa = 2.0d0/(1.d0+320.0d0*DEXP(-(Y(Vm)+
     '        75.0d0-PARAM(SHIFT_hNa))/10.0d0)) !(ms^-1) Equation 40
            INa = PARAM(gNa)*Y(mNa)**3*Y(hNa)*(Y(Vm)-Emh) !(uA.mm^-2) Equation 34
            DY(mNa) = alpha_mNa*(1.0d0-Y(mNa))-beta_mNa*Y(mNa) !(ms^-1) Equation 35

            DY(hNa) = alpha_hNa*(1.0d0-Y(hNa))-beta_hNa*Y(hNa) !(ms^-1) Equation 36

C           Background Na+ current, IbNa
            IbNa = PARAM(gbNa)*(Y(Vm)-ENa) !(uA.mm^-2) Equation 41

C           Persistent Na+ current, IpNa
            IpNa = PARAM(gpNa)*(Y(Vm)-ENa)/(1.0d0+DEXP(-(Y(Vm)+
     '        52.0d0)/8.0d0)) !(uA.mm^-2) Equation 42

C           Background Ca++ current, IbCa
            IbCa  = PARAM(gbCa)*(Y(Vm)-ECa) !(uA.mm^-2) Equation 43

C           L-type Ca++ channel, ICaL
            IF(DABS(Y(Vm)+24.0d0-PARAM(SHIFT_dCa)).LT.0.001d0) THEN
              alpha_dCa = 0.12d0*PARAM(epsilon_d) !(ms^-1) Equation 54
              beta_dCa = 0.12d0*PARAM(epsilon_d) !(ms^-1) Equation 55
            ELSE
              alpha_dCa = 0.03d0*PARAM(epsilon_d)*(Y(Vm)+24.0d0-
     '          PARAM(SHIFT_dCa))/(1.0d0-DEXP(-(Y(Vm)+24.d0-
     '          PARAM(SHIFT_dCa))/
     '          PARAM(STEEP_dCa))) !(ms^-1) Equation 54
              beta_dCa = -0.012d0*PARAM(epsilon_d)*
     '          (Y(Vm)+24.0d0-PARAM(SHIFT_dCa))/
     '          (1.0d0-DEXP((Y(Vm)+24.0d0-PARAM(SHIFT_dCa))/
     '          (2.5d0*PARAM(STEEP_dCa)))) !(ms^-1) Equation 54
            ENDIF
            IF(DABS(Y(Vm)+34.0d0-PARAM(SHIFT_fCa)).LT.0.001d0) THEN
              alpha_fCa = 0.00625d0*PARAM(epsilon_f)*
     '          PARAM(STEEP_fCa) !(ms^-1) Equation 56
            ELSE
              alpha_fCa = -0.00625d0*PARAM(epsilon_f)*(Y(Vm)+
     '          34.0d0-PARAM(SHIFT_fCa))/(1.0d0-DEXP((Y(Vm)+
     '          34.0d0-PARAM(SHIFT_fCa))/
     '          PARAM(STEEP_fCa))) !(ms^-1) Equation 56
            ENDIF
            beta_fCa = 0.012d0*PARAM(epsilon_f)/(1.0d0+
     '        DEXP(-(Y(Vm)+34.0d0-PARAM(SHIFT_fCa))/
     '        PARAM(STEEP_fCa))) !(ms^-1) Equation 57

            z1 = PARAM(PCa)*(1.0d0-PARAM(ACh)/
     '        (PARAM(ACh)+PARAM(KmAChICaL)))
            z2 = DEXP(PARAM(Esurf)/RTF)*Y(Ki)-DEXP(-(Y(Vm)-
     '        PARAM(Esurf))/RTF)*PARAM(Ko)
            ICaLK = pHi_ICaL*(1.0d0-PARAM(ICaLFract))*PARAM(PCaK)*
     '        Y(dCa)*Y(fCa)*Y(f2Ca)*z1*z2 !(uA.mm^-2) Equation 44
            ICaLKDS = pHi_ICaL*PARAM(ICaLFract)*PARAM(PCaK)*Y(dCa)*
     '        Y(fCa)*Y(f2CaDS)*z1*z2 !(uA.mm^-2) Equation 47
            z2 = DEXP(PARAM(Esurf)/RTF)*Y(Nai)-DEXP(-(Y(Vm)-
     '        PARAM(Esurf))/RTF)*PARAM(Nao)

            ICaLNa = pHi_ICaL*(1.0d0-PARAM(ICaLFract))*PARAM(PCaNa)*
     '        Y(dCa)*Y(fCa)*Y(f2Ca)*z1*z2 !(uA.mm^-2) Equation 45
            ICaLNaDS = pHi_ICaL*PARAM(ICaLFract)*PARAM(PCaNa)*Y(dCa)*
     '        Y(fCa)*Y(f2CaDS)*z1*z2 !(uA.mm^-2) Equation 48
            z1 = PARAM(PCa)*(1.0d0-PARAM(ACh)/
     '        (PARAM(ACh)+PARAM(KmAChICaL)))*
     '        (Y(Vm)-PARAM(Esurf))/(RTF*(1.0d0-
     '        DEXP(-2.0d0*(Y(Vm)-PARAM(ESurf))/RTF)))
            z2 = DEXP(2.0d0*PARAM(Esurf)/RTF)*Y(Cai)-
     '        DEXP(-2.0d0*(Y(Vm)-PARAM(Esurf))/RTF)*PARAM(Cao)
            ICaLCa = pHi_ICaL*(1.0d0-PARAM(ICaLFract))*4.0d0*
     '        Y(dCa)*Y(fCa)*Y(f2Ca)*z1*z2 !(uA.mm^-2) Equation 46
            z2 = DEXP(2.0d0*PARAM(Esurf)/RTF)*Y(CaDS)-
     '        DEXP(-2.0d0*(Y(Vm)-PARAM(Esurf))/RTF)*PARAM(Cao)
            ICaLCaDS = pHi_ICaL*PARAM(ICaLFract)*4.0d0*
     '        Y(dCa)*Y(fCa)*Y(f2CaDS)*z1*z2 !(uA.mm^-2) Equation 49

            DY(dCa) = alpha_dCa*(1.0d0-Y(dCa))-beta_dCa*Y(dCa) !(ms^-1) Equation 50

            DY(fCa) = alpha_fCa*(1.0d0-Y(fCa))-beta_fCa*Y(fCa) !(ms^-1) Equation 51

            DY(f2Ca) = PARAM(alpha_CaInact)*(1.0d0-
     '        Y(Cai)/(Y(Cai)+PARAM(KmCaInact))-Y(f2Ca)) !(ms^-1) Equation 52

            DY(f2CaDS) = PARAM(alpha_CaDSInact)*(1.0d0-
     '        Y(CaDS)/(Y(CaDS)+PARAM(KmCaDSInact))-
     '        Y(f2CaDS)) !(ms^-1) Equation 53

C           Na+/K+ Exchange, INaK (Ip)
            INaK = PARAM(INaKmax)*PARAM(Ko)/
     '        (PARAM(Ko)+PARAM(KmK))*
     '        Y(Nai)/(Y(Nai)+PARAM(KmNa)) !(uA.mm^-2) Equation 58

C           Na+/Ca++ exchange current, INaCa
            VRTF = Y(Vm)/RTF
            c1 = DEXP(PARAM(gamma_NaCa)*(nNaCa-2.0d0)*VRTF)
            c2 = DEXP((PARAM(gamma_NaCa)-1.0d0)*(nNaCa-2.0d0)*VRTF)
            z1 = c1*Y(Nai)**nNaCa*PARAM(Cao)-
     '        c2*PARAM(Nao)**nNaCa*Y(Cai)
            z2 = 1.0d0+PARAM(dNaCa)*(PARAM(Nao)**nNaCa*
     '        Y(Cai)+Y(Nai)**nNaCa*PARAM(Cao))
            INaCa = pHi_INaCa*(1.0d0-PARAM(INaCaFract))*PARAM(INaCamax)*
     '        (z1/z2)*1.0d0/(1.0d0+Y(Cai)/0.0069d0) !(uA.mm^-2) Equation 59

C           Diadic-space Na+/Ca++ exchange current, INaCaDS
            z1 = c1*Y(Nai)**nNaCa*PARAM(Cao)-
     '        c2*PARAM(Nao)**nNaCa*Y(CaDS)
            z2 = 1.0d0+PARAM(dNaCa)*(PARAM(Nao)**nNaCa*
     '        Y(CaDS)+Y(Nai)**nNaCa*PARAM(Cao))
            INaCaDS = pHi_INaCa*PARAM(INaCaFract)*PARAM(INaCamax)*
     '        (z1/z2)*1.0d0/(1.0d0+Y(CaDS)/0.0069d0) !(uA.mm^-2) Equation 60

            IF(MODEL(USE_SAC).EQ.1) THEN
C             Stretch activated channels

              f_stretch = 1.0d0/(1.0d0+DEXP(-2.0d0*
     '          PARAM(gamma_SACSL)*PARAM(SLRef)*
     '          (Y(ExtensionRatio)-PARAM(SLHST)))) !() Equation 61


C             K+ stretch activated channel
              IK_stretch = PARAM(gK_stretch)*f_stretch*
     '          (Y(Vm)-EK) !(uA.mm^-2) Equation 62

C             Na+ stretch activated channel
              INa_stretch = PARAM(gNa_stretch)*f_stretch*
     '          (Y(Vm)-ENa) !(uA.mm^-2) Equation 63

C             Ca++ stret ch activated channel
              ICa_stretch = PARAM(gCa_stretch)*f_stretch*
     '          (Y(Vm)-ECa) !(uA.mm^-2) Equation 64

C             INs stretch activated channel
              INs_stretch = PARAM(gNs_stretch)*f_stretch*
     '          (Y(Vm)-PARAM(ENs_stretch)) !(uA.mm^-2) Equation 65

C             IAn stretch activated channel
              IAn_stretch = PARAM(gAn_stretch)*f_stretch*
     '          (Y(Vm)-PARAM(EAn_stretch)) !(uA.mm^-2) Equation 66
            ELSE
              INa_stretch = 0.0d0
              IK_stretch = 0.0d0
              ICa_stretch = 0.0d0
              INs_stretch = 0.0d0
              IAn_stretch = 0.0d0
            ENDIF

C           SR Ca++ uptake flux, jup
            c1 = PARAM(KcyCa)*PARAM(Kxcs)/PARAM(KSRCa)
            c2 = Y(Cai)+Y(Caup)*c1+PARAM(KcyCa)*PARAM(Kxcs)+
     '        PARAM(KcyCa)
            FractCaUptakeSites = Y(Cai)/c2 !() Equation 68
            FractBackSRSites = Y(Caup)*c1/c2 !() Equation 69
            jup = PARAM(alpha_up)*FractCaUptakeSites-
     '        PARAM(beta_up)*FractBackSRSites !(mmol.L^-1.ms^-1) Equation 67

C           SR Ca++ translocation flux, jtr
            jtr = PARAM(alpha_tr)*(Y(Caup)-Y(Carel)) !(mmol.L^-1.ms^-1) Equation 70

C           SR Ca++ leakage flux, jleak
            jleak = DEXP(PARAM(gamma_SRSL)*PARAM(SL))*
     '        PARAM(alpha_SRLeak)*Y(Carel) !(mmol.L^-1.ms^-1) Equation 71

C           SR Ca++ release flux, jrel
            OpenReleaseChannelFract=(Y(ActivatorFract)/
     '        (Y(ActivatorFract)+0.25d0))**2 !() Equation 72
            jrel = PARAM(alpha_rel)*OpenReleaseChannelFract*
     '        Y(Carel) !(mmol.L^-1.ms^-1) Equation 71
            RegulatoryBindingSite=(Y(Cai)/(Y(Cai)+
     '        PARAM(KmCa))+(1.0d0-Y(Cai)/(Y(Cai)+
     '        PARAM(KmCa)))*Y(CaDS)/(Y(CaDS)+
     '        PARAM(KmCaDS)))**2 !() Equation 78
            PrecursorFract=1.0d0-Y(ActivatorFract)-Y(ProductFract) !() Equation 75
            alpha_act = 0.5d0*RegulatoryBindingSite !(ms^-1) Equation 76
            alpha_inact = 0.06d0+0.5d0*RegulatoryBindingSite !(ms^-1) Equation 77
            IF(Y(Vm).LT.-50.0d0) THEN
              DY(ActivatorFract) = PARAM(epsilon_sr)*
     '          (alpha_act*PrecursorFract-alpha_inact*Y(ActivatorFract)) !(ms^-1) Equation 73
              DY(ProductFract) = PARAM(epsilon_sr)*
     '          (alpha_inact*Y(ActivatorFract)-0.001d0*Y(ProductFract)) !(ms^-1) Equation 74
            ELSE
              DY(ActivatorFract) = alpha_act*PrecursorFract-
     '          alpha_inact*Y(ActivatorFract) !(ms^-1) Equation 73
              DY(ProductFract) = alpha_inact*Y(ActivatorFract)-
     '          0.001d0*Y(ProductFract) !(ms^-1) Equation 74
            ENDIF

C           DS Ca++ decay flux, jdecay
            jdecay = PARAM(alpha_CaDSdecay)*Y(CaDS)

C           Membrane potential change
            DY(Vm) = (IStim-IK1-Ito-IKr1-IKr2-IKs-IKNa-IbK-IKATP-
     '        IKACh-INa-IbNa-IpNa-ICaLK-ICaLNa-ICaLCa-
     '        ICaLKDS-ICaLNaDS-ICaLCaDS-IbCa-INaK-INaCa-INaCaDS-
     '        IK_stretch-INa_stretch-ICa_stretch-INs_stretch-
     '        IAn_stretch)/PARAM(Cm) !(mV.ms^-1) Equation 79

C           Concentration changes
C           Intracellular K+, Ki
            DY(Ki) = -PARAM(Am)/(PARAM(Vi)*F)*
     '        (IK1+Ito+IKr1+IKr2+IKs+IKNa+IbK+IKATP+IKACh+ICaLK+
     '        ICaLKDS-1.0d0/(nNaK-1.0d0)*INaK+IK_stretch) !(mmol.L^-1.ms^-1) Equation 80

C           Intracellular Na+, Nai
            DY(Nai) = -PARAM(Am)/(PARAM(Vi)*F)*
     '        (INa+IbNa*PARAM(Nao)/140.0d0+IpNa+ICaLNa+ICaLNaDS+
     '        nNaK/(nNaK-1.0d0)*INaK+
     '        nNaCa/(nNaCa-2.0d0)*(INaCa+INaCaDS)) !(mmol.L^-1.ms^-1) Equation 81

C           Ca++ bound to calmodulin, CaCALM
            DY(CaCALM) = PARAM(alpha_CALM)*(PARAM(CALM)-
     '        Y(CaCALM))*Y(Cai)-PARAM(beta_CALM)*
     '        Y(CaCALM) !(mmol.L^-1.ms^-1) Equation 86

C *** HMT
C ***       Troponin kinetics - Ca++ binding to troponin, CaTROP
            IF (DABS(Y(IsometricTension)).LT.1.0d-7) THEN
              TTo = 0.0d0
            ELSE
              TTo = Y(Tension)/(PARAM(gamma)*Y(IsometricTension))
            ENDIF
            DY(CaTROP) = PARAM(Rho0)*Y(Cai)*(PARAM(CaTropMax)-
     '        Y(CaTrop))-PARAM(Rho1)*(1.0d0-TTo)*Y(CaTrop)
            !(mmol.L^-1.ms^-1)

C ***       Tropomyosin Kinetics
            !length dependence for n
            Tm_n   = PARAM(Tm_n_0)*(1.d0+PARAM(beta1)*
     '        (Y(ExtensionRatio)-1.d0))
            !length dependence for p50
            Tm_p50 = PARAM(Tm_p50_0)*(1.d0+PARAM(beta2)*
     '        (Y(ExtensionRatio)-1.d0))
            C50 = 10**(3.d0-Tm_p50) !mM

C DPN 05/08/98 - scale dz/dt by normalised [Ca]b
            Cb_norm = Y(CaTROP)/PARAM(CaTropMax)

            DY(z)=PARAM(alfa0)*(((Y(CaTROP)/C50)*Cb_norm)**Tm_n*
     '        (1.d0-Y(z)) - Y(z))

C           Ca++ bound to indicator, CaIND
            DY(CaIND) = PARAM(alpha_IND)*(PARAM(IND)-
     '        Y(CaIND))*Y(Cai)-PARAM(beta_IND)*
     '        Y(CaIND) !(mmol.L^-1.ms^-1) Equation 88

C           Intracellular Ca++, Cai
            DY(Cai) = -PARAM(Am)/(2.0d0*PARAM(Vi)*F)*
     '        (ICaLCa+IbCa-2.0d0/(nNaCa-2.0d0)*INaCa+ICa_stretch)+
     '        jrel*PARAM(Vrel)/PARAM(Vi)+
     '        jleak*PARAM(Vrel)/PARAM(Vi)-jup+
     '        jdecay*PARAM(VDS)/PARAM(Vi)-DY(CaCALM)-
     '        DY(CaTROP)-DY(CaIND) !(mmol.L^-1.ms^-1) Equation 82

C           Diadic space Ca++, CaDS
            DY(CaDS) = -PARAM(Am)/(2.0d0*PARAM(VDS)*F)*
     '        (ICaLCaDS-2.0d0/(nNaCa-2.0d0)*INaCaDS)-jdecay !(mmol.L^-1.ms^-1) Equation 83

C           SR uptake Ca++, Caup
            DY(Caup) = jup*PARAM(Vi)/PARAM(Vup)-
     '        jtr !(mmol.L^-1.ms^-1) Equation 84

C           SR release Ca++, Carel
            DY(Carel) = jtr*PARAM(Vup)/PARAM(Vrel)-
     '        jrel-jleak !(mmol.L^-1.ms^-1) Equation 85

            IF(CONTROL(RETURN_CURRENT).NE.0) THEN
              DERIVED(DIK1) = IK1
              DERIVED(DIto) = Ito
              DERIVED(DIKr1) = IKr1
              DERIVED(DIKr2) = IKr2
              DERIVED(DIKs) = IKs
              DERIVED(DIKNa) = IKNa
              DERIVED(DIbK) = IbK
              DERIVED(DIKATP) = IKATP
              DERIVED(DIKACh) = IKACh
              DERIVED(DINa) = INa
              DERIVED(DIbNa) = IbNa
              DERIVED(DIpNa) = IpNa
              DERIVED(DICaLK) = ICaLK
              DERIVED(DICaLNa) = ICaLNa
              DERIVED(DICaLCa) = ICaLCa
              DERIVED(DICaLKDS) = ICaLKDS
              DERIVED(DICaLNaDS) = ICaLNaDS
              DERIVED(DICaLCaDS) = ICaLCaDS
              DERIVED(DIbCa) = IbCa
              DERIVED(DINaK) = INaK
              DERIVED(DINaCa) = INaCa
              DERIVED(DINaCaDS) = INaCaDS
              DERIVED(DIK_stretch) = IK_stretch
              DERIVED(DINa_stretch) = INa_stretch
              DERIVED(DICa_stretch) = ICa_stretch
              DERIVED(DINs_stretch) = INs_stretch
              DERIVED(DIAn_stretch) = IAn_stretch
              DERIVED(Djup) = jup
              DERIVED(Djtr) = jtr
              DERIVED(Djrel) = jrel
              DERIVED(Djleak) = jleak
              DERIVED(Djdecay) = jdecay
              DERIVED(DIK) = IKr1+IKr2+IKs
              IK_tot = IK1+Ito+IKr1+IKr2+IKs+IKNa+IbK+IKATP+IKACh+ICaLK+
     '          ICaLKDS+INaK+IK_stretch
              DERIVED(DIK_tot) = IK_tot
              INa_tot = INa+IbNa+IpNa+ICaLNa+ICaLNaDS+INaK+INaCa+INaCaDS
              DERIVED(DINa_tot) = INa_tot
              ICa_tot = ICaLCa+IbCa+INaCa+ICa_stretch
              DERIVED(DICa_tot) = ICa_tot
              ICaDS_tot = ICaLCaDS+INaCaDS
              DERIVED(DICaDS_tot) = ICaDS_tot
            ENDIF
          ELSE !evaluating non-ODE state variables

C *** HMT *********

C ***       Tension-length-pCa - Isometric Tension
            Y(IsometricTension) = PARAM(Tref)*(1.d0+PARAM(beta0)*
     '        (Y(ExtensionRatio)-1.d0))*Y(z)

C ***       X-bridge kinetics - Active tension
c           Fading memory model...
            LENGTH_DIFF = Y(ExtensionRatio)-Y(ExtensionRatio_prev)

C           1st Fading Memory term
            dPHI1=0.5d0*LENGTH_DIFF*(1.d0+DEXP(-PARAM(alpha1)*
     '        TIME(DTCell)))
            Y(phi1)=DEXP(-PARAM(alpha1)*TIME(DTCell))*Y(phi1)+dPHI1

C           2nd Fading Memory term
            dPHI2=0.5d0*LENGTH_DIFF*(1.d0+DEXP(-PARAM(alpha2)*
     '        TIME(DTCell)))
            Y(phi2)=DEXP(-PARAM(alpha2)*TIME(DTCell))*Y(phi2)+dPHI2

C           3rd Fading Memory term
            dPHI3=0.5d0*LENGTH_DIFF*(1.d0+DEXP(-PARAM(alpha3)*
     '        TIME(DTCell)))
            Y(phi3)=DEXP(-PARAM(alpha3)*TIME(DTCell))*Y(phi3)+dPHI3

C *** DPN 06 October 2000 - Having trouble exporting the values of phi
C           to an exelem file when the values are really small...
            IF (DABS(Y(phi1)).LT.ZERO_TOL) Y(phi1) = 0.0d0
            IF (DABS(Y(phi2)).LT.ZERO_TOL) Y(phi2) = 0.0d0
            IF (DABS(Y(phi3)).LT.ZERO_TOL) Y(phi3) = 0.0d0

            Q = PARAM(A1)*Y(phi1)+PARAM(A2)*Y(phi2)+PARAM(A3)*Y(phi3)

C *** DPN 20 July 2000 - Need to fix HMT for lengthening experiments
            IF(Q.LT.0.00d0) THEN
              Y(Tension)=Y(IsometricTension)*(1.d0+PARAM(a)*Q)/(1.d0-Q)
            ELSE
              Y(Tension) = Y(IsometricTension)
            ENDIF

C *** DPN 02 October 2000 - Adding some viscous damping
            Y(Tension) = Y(Tension) - (PARAM(DampingCoeff)
     '        *DABS(LENGTH_DIFF)/TIME(DTCell))

C           Update previous extension ratio
            Y(ExtensionRatio_prev) = Y(ExtensionRatio)

          ENDIF !ODE/non-ODE
        ENDIF !Errors II
      ENDIF !Errors I
      RETURN
      END
