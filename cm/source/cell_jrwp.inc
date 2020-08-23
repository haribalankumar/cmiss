C  cell_jrwp.inc
C
C  Array elements for the Jafri-Rice-Winslow-Physiome ionic current model.
C
C  NOTE: The base units to be used in cell models are:
C            Area - mm^2
C   Concentration - mMol
C    Conductivity - mS/mm
C     Capacitance - uF
C         Current - uA
C          Length - mm
C       Potential - mV
C     Temperature - K
C            Time - ms
C          Volume - mm^3 (uL)

C *** CONTROL array
      INTEGER RETURN_CURRENT
      PARAMETER(RETURN_CURRENT = 1)

C *** MODEL array
      INTEGER FORCE_TYPE
      PARAMETER(FORCE_TYPE = 1)

C *** ARO - additional real out array

C *** PARAMETERS array
      INTEGER TEMP,Acap,Vmyo,VJSR,VNSR,VSS,Ko,Nao,Cao,GKr,GKs,Gto1,
     +  GK1,GKp,GNa,kNaCa,KmNa,KmCa,KmK1,Ksat,eta,INaKmax,KmNai,KmKo,
     +  IpCamax,KmpCa,GCab,GNab,K1_inf_const,v1,Kfb,Krb,KSR,Nfb,Nrb,
     +  vmaxf,vmaxr,tautr,tauxfer,kaplus,kaminus,fL,gL,bL,aL,omega,
     +  Pscale,PCa,PK,ICahalf,LTRPNtot,HTRPNtot,
     +  khtrpn_plus,khtrpn_minus,kltrpn_plus,kltrpn_minus,
     +  CMDNtot,CSQNtot,EGTAtot,KmCMDN,KmCSQN,KmEGTA,
     +  Fisotonic,La,Lm_prime,Lb,Lz,Lm,SL0,SLmax,SLmin,SLset,Vmax,
     +  smallk1,smallk_1,KCa,f,g0,g0_off,f01,f12,f23,g0_01,g0_12,
     +  g0_23,g0_off_01,small_k_norm,pathsum,fnormnum,fnormmax,
     +  fnormnum2,fnormmax2,GCl_Ca,Kd_Cl_Ca,N_Cl_Ca,Vrev_Cl_Ca,
     +  xNa_sustain
      REAL*8 Farad

C     PARAMETER(          Cm =  1 ) !membrane cap. per unit area (uF/mm^2) given in cell_reserved.inc
C     PARAMETER(          Am =  2 ) !fibre surface to volume ratio - 1/mm  given in cell_reserved.inc
      PARAMETER(        TEMP =  3 ) !temperature (K)
      PARAMETER(        Acap =  4 ) !capacitive membrane area (mm^2)
      PARAMETER(        Vmyo =  5 ) !myoplasmic volume (uL)
      PARAMETER(        VJSR =  6 ) !junctional SR volume (uL)
      PARAMETER(        VNSR =  7 ) !network SR volume (uL)
      PARAMETER(         VSS =  8 ) !subspace volume (uL)
      PARAMETER(          Ko =  9 ) !extracellular [K+] (mMol)
      PARAMETER(         Nao = 10 ) !extracellular [Na+] (mMol)
      PARAMETER(         Cao = 11 ) !extracellular [Ca++] (mMol)
      PARAMETER(         GKr = 12 ) !peak IKr conductance (mS/mm^2)
      PARAMETER(         GKs = 13 ) !peak IKs conductance (mS)
      PARAMETER(        Gto1 = 14 ) !peak Ito1 conductance (mS)
      PARAMETER(         GK1 = 15 ) !peak IK1 conductance (mS)
      PARAMETER(         GKp = 16 ) !peak IKp conductance (mS)
      PARAMETER(         GNa = 17 ) !peak INa conductance (mS)
      PARAMETER(       kNaCa = 18 ) !scaling factor of Na+/Ca++ exchange (uA)
      PARAMETER(        KmNa = 19 ) !Na+ half saturation constant for Na+/Ca++ exchange (mMol)
      PARAMETER(        KmCa = 20 ) !Ca++ half saturation constant for Na+/Ca++ exchange (mMol)
      PARAMETER(        KmK1 = 21 ) !Ca++ half saturation constant for IK1 (mMol)
      PARAMETER(        Ksat = 22 ) !Na+/Ca++ exchange saturation factor at negative potentials (???)
      PARAMETER(         eta = 23 ) !controls voltage dependence of Na+/Ca++ exchange
      PARAMETER(     INaKmax = 24 ) !maximum Na+/K+ pump current (uA)
      PARAMETER(       KmNai = 25 ) !Na+ half saturation constant for Na+/K+ pump (mMol)
      PARAMETER(        KmKo = 26 ) !K+ half saturation constant for Na+/K+ pump (mMol)
      PARAMETER(     IpCamax = 27 ) !maximum sarcolemmal Ca++ pump current (uA)
      PARAMETER(       KmpCa = 28 ) !half saturation constant for sarcolemmal Ca++ pump (mMol)
      PARAMETER(        GCab = 29 ) !maximum background Ca++ current conductance (mS)
      PARAMETER(        GNab = 30 ) !maximum background Na+ current conductance (mS)
      PARAMETER(K1_inf_const = 31 ) !
      PARAMETER(          v1 = 32 ) !maximum RyR channel Ca++ flux (1/ms ???)
      PARAMETER(         Kfb = 33 ) !forward half saturation constant for Ca++ ATPase (mMol)
      PARAMETER(         Krb = 34 ) !backward half saturation constant for Ca++ ATPase (mMol)
      PARAMETER(         KSR = 35 ) !scaling factor for Ca++ ATPase (???)
      PARAMETER(         Nfb = 36 ) !forward cooperativity constant for Ca++ ATPase (???)
      PARAMETER(         Nrb = 37 ) !backward cooperativity constant for Ca++ ATPase (???)
      PARAMETER(       vmaxf = 38 ) !Ca++ ATPase forward rate parameter (mMol/ms)
      PARAMETER(       vmaxr = 39 ) !Ca++ ATPase backward rate parameter (mMol/ms)
      PARAMETER(       tautr = 40 )   !time constant for transfer from NSR to JSR (ms)
      PARAMETER(     tauxfer = 41 ) !time constant for transfer from SS to myoplasm (ms)
      PARAMETER(      kaplus = 42 ) !
      PARAMETER(     kaminus = 43 ) !
      PARAMETER(          fL = 44 ) !transition rate into open state (1/ms)
      PARAMETER(          gL = 45 ) !transition rate out of open state (1/ms)
      PARAMETER(          bL = 46 ) !mode transition parameter (???)
      PARAMETER(          aL = 47 ) !mode transition parameter (???)
      PARAMETER(       omega = 48 ) !mode transition parameter (1/ms)
      PARAMETER(      Pscale = 49 ) !scaling factor for both PK and PCa
      PARAMETER(         PCa = 50 ) !L-type Ca++ channel permeability to Ca++ (cm/s)
      PARAMETER(          PK = 51 ) !L-type Ca++ channel permeability to K+ (cm/s)
      PARAMETER(     ICahalf = 52 ) !ICa level that reduces PK by half (uA)
      PARAMETER(    LTRPNtot = 53 ) !total troponin low affinity site conc (mMol)
      PARAMETER(    HTRPNtot = 54 ) !total troponin high affinity site conc (mMol)
      PARAMETER( khtrpn_plus = 55 ) !Ca++ on rate for troponin high affinity sites (1/[mMol*ms])
      PARAMETER(khtrpn_minus = 56 ) !Ca++ off rate for troponin high affinity sites (1/ms)
      PARAMETER( kltrpn_plus = 57 ) !Ca++ on rate for troponin low affinity sites (1/[mM*ms])
      PARAMETER(kltrpn_minus = 58 ) !Ca++ off rate for troponin low affinity sites (1/ms)
      PARAMETER(     CMDNtot = 59 ) !total myoplasmic calmodulin concentration (mMol)
      PARAMETER(     CSQNtot = 60 ) !total NSR calsequestrin concentration (mMol)
      PARAMETER(     EGTAtot = 61 ) !total myoplasmic EGTA concentration (mMol)
      PARAMETER(      KmCMDN = 62 ) !Ca++ half sat constant for calmodulin (mMol)
      PARAMETER(      KmCSQN = 63 ) !Ca++ half sat constant for calsequestrin (mMol)
      PARAMETER(      KmEGTA = 64 ) !Ca++ half sat constant for EGTA (mMol)
      PARAMETER(   Fisotonic = 65 ) !
      PARAMETER(          La = 66 ) !
      PARAMETER(    Lm_prime = 67 ) !
      PARAMETER(          Lb = 68 ) !unused ? (PJH)
      PARAMETER(          Lz = 69 ) !
      PARAMETER(          Lm = 70 ) !
      PARAMETER(         SL0 = 71 ) !
      PARAMETER(       SLmax = 72 ) !
      PARAMETER(       SLmin = 73 ) !
      PARAMETER(       SLset = 74 ) !
      PARAMETER(        Vmax = 75 ) !
      PARAMETER(     smallk1 = 76 ) !unused ? (PJH)
      PARAMETER(    smallk_1 = 77 ) !unused ? (PJH)
      PARAMETER(         KCa = 78 ) !
      PARAMETER(           f = 79 ) !unused ? (PJH)
      PARAMETER(          g0 = 80 ) !unused ? (PJH)
      PARAMETER(      g0_off = 81 ) !unused ? (PJH)
      PARAMETER(         f01 = 82 ) !
      PARAMETER(         f12 = 83 ) !
      PARAMETER(         f23 = 84 ) !
      PARAMETER(       g0_01 = 85 ) !
      PARAMETER(       g0_12 = 86 ) !
      PARAMETER(       g0_23 = 87 ) !
      PARAMETER(   g0_off_01 = 88 ) !
      PARAMETER(small_k_norm = 89 ) !
      PARAMETER(     pathsum = 90 ) !unused ? (PJH)
      PARAMETER(    fnormnum = 91 ) !unused ? (PJH)
      PARAMETER(    fnormmax = 92 ) !
      PARAMETER(   fnormnum2 = 93 ) !unused ? (PJH)
      PARAMETER(   fnormmax2 = 94 ) !unused ? (PJH)
      PARAMETER(      GCl_Ca = 95 ) !Ca activited Cl Conductance (mS/mm^2)
      PARAMETER(    Kd_Cl_Ca = 96 ) !Ca act Cl 1/2 activation constant (mMol)
      PARAMETER(     N_Cl_Ca = 97 ) !Ca act Cl Hill coeff
      PARAMETER(  Vrev_Cl_Ca = 98 ) !Ca activated Cl reversal potential (mV)
      PARAMETER(  xNa_sustain= 99 ) !Scaling factor for Na_sus ()

C *** DERIVED array
      INTEGER INa,IKr,IKs,Ito1,IK1,IKp,InsK,InsNa,INaCa,INaCa_SS,INaK,
     +  IpCa,ICab,INab,ICa,ICaK,ICl_Ca,Jup,Jrel,Jtrpn,Jtr,Jxfer,Istim

      PARAMETER(         INa =  1 ) !Na+ current
      PARAMETER(         IKr =  2 ) !rapid activating delayed rectifier K+ current
      PARAMETER(         IKs =  3 ) !slow activating delayed rectifier K+ current
      PARAMETER(        Ito1 =  4 ) !transient outward K+ current
      PARAMETER(         IK1 =  5 ) !time independent K+ current
      PARAMETER(         IKp =  6 ) !plateau K+ current
      PARAMETER(        InsK =  7 ) !non-specific Ca-activated K+ current
      PARAMETER(       InsNa =  8 ) !non-specific Ca-activated Na+ current
      PARAMETER(       INaCa =  9 ) !Na+/Ca++ exchanger current (in cytoplasm)
      PARAMETER(    INaCa_SS = 10 ) !Na+/Ca++ exchanger current (in subspace???)
      PARAMETER(        INaK = 11 ) !Na+/K+ pump current
      PARAMETER(        IpCa = 12 ) !sarcolemmal Ca++ pump current
      PARAMETER(        ICab = 13 ) !Ca++ background current
      PARAMETER(        INab = 14 ) !Na+ background current
      PARAMETER(         ICa = 15 ) !L-type channel Ca++ current
      PARAMETER(        ICaK = 16 ) !L-type channel K+ current
      PARAMETER(      ICl_Ca = 17 ) !Ca activated Cl current
      PARAMETER(         Jup = 18 ) !SERCA pump flux
      PARAMETER(        Jrel = 19 ) !JSR Ca++ release flux
      PARAMETER(       Jtrpn = 20 ) !troponin-Ca++ binding/unbinding flux
      PARAMETER(         Jtr = 21 ) !Ca++ flux from NSR to JST
      PARAMETER(       Jxfer = 22 ) !Ca++ flux from SS to myoplasm
      PARAMETER(       Istim = 23 ) !Stimulus current

C *** Y and DY arrays
      INTEGER mNa,hNa,jNa,mto1,y1to1,y2to1,Nai,Ki,Cai,CaNSR,CaSS,
     +  CaJSR,C1_RyR,O1_RyR,O2_RyR,C2_RyR,xKr,xKs,C0,C1,C2,C3,C4,Open,
     +  CCa0,CCa1,CCa2,CCa3,CCa4,yCa,LTRPNCa,HTRPNCa,xto1,yto1,SL,A,T1,
     +  U,T2,T3,Force,CAMK,
     +  rC1_RyR,rC2_RyR,rO3_RyR,rO4_RyR,rC5_RyR,rC6_RyR,yNa2

C     PARAMETER(          Vm =  1 ) !Membrane potential (mV) given in cell_reserved.inc
      PARAMETER(         mNa =  2 ) !INa activation variable ()
      PARAMETER(         hNa =  3 ) !INa fast inactivation variable ()
      PARAMETER(         jNa =  4 ) !INa slow inactivation variable ()
      PARAMETER(        mto1 =  5 ) !Ito1 (JLG) activation variable ()
      PARAMETER(       y1to1 =  6 ) !Ito1 (JLG) inactivation state y1 variable ()
      PARAMETER(       y2to1 =  7 ) !Ito1 (JLG) inactivation state y2 variable ()
      PARAMETER(         Nai =  8 ) !Intracellular [Na+] (mMol)
      PARAMETER(          Ki =  9 ) !Intraceullar [K+]  (mMol)
      PARAMETER(         Cai = 10 ) !Intraceullar [Ca++] (mMol)
      PARAMETER(       CaNSR = 11 ) !Network SR [Ca++] (mMol)
      PARAMETER(        CaSS = 12 ) !Subspace [Ca++] (mMol) ???
      PARAMETER(       CaJSR = 13 ) !Junctional SR [Ca++] (mMol)
      PARAMETER(      C1_RyR = 14 ) !RyR channel closed state 1 variable ()
      PARAMETER(      O1_RyR = 15 ) !RyR channel open state 1 variable ()
      PARAMETER(      O2_RyR = 16 ) !RyR channel open state 2 variable ()
      PARAMETER(      C2_RyR = 17 ) !RyR channel closed state 2 variable ()
      PARAMETER(         xKr = 18 ) !IKr activation variable ()
      PARAMETER(         xKs = 19 ) !IKs activation variable ()
      PARAMETER(          C0 = 20 ) !L-type Ca channel: normal mode, closed state 0 variable ()
      PARAMETER(          C1 = 21 ) !L-type Ca channel: normal mode, closed state 1 variable ()
      PARAMETER(          C2 = 22 ) !L-type Ca channel: normal mode, closed state 2 variable ()
      PARAMETER(          C3 = 23 ) !L-type Ca channel: normal mode, closed state 3 variable ()
      PARAMETER(          C4 = 24 ) !L-type Ca channel: normal mode, closed state 4 variable ()
      PARAMETER(        Open = 25 ) !L-type Ca channel: normal mode, open state variable ()
      PARAMETER(        CCa0 = 26 ) !L-type Ca channel: calcium mode, closed state 0 variable ()
      PARAMETER(        CCa1 = 27 ) !L-type Ca channel: calcium mode, closed state 1 variable ()
      PARAMETER(        CCa2 = 28 ) !L-type Ca channel: calcium mode, closed state 2 variable ()
      PARAMETER(        CCa3 = 29 ) !L-type Ca channel: calcium mode, closed state 3 variable ()
      PARAMETER(        CCa4 = 30 ) !L-type Ca channel: calcium mode, closed state 4 variable ()
      PARAMETER(         yCa = 31 ) !L-type Ca channel: inactivation variable
      PARAMETER(     LTRPNCa = 32 ) !Tn binding: Ca fraction bound to low affinity site ()
      PARAMETER(     HTRPNCa = 33 ) !Tn binding: Ca fraction bound to high affinity site ()
      PARAMETER(        xto1 = 34 ) !Ito1 (RLW) activation variable ()
      PARAMETER(        yto1 = 35 ) !Ito1 (RLW) inactivation variable ()
      PARAMETER(          SL = 36 ) !
      PARAMETER(           A = 37 ) !
      PARAMETER(          T1 = 38 ) !
      PARAMETER(           U = 39 ) !
      PARAMETER(          T2 = 40 ) !
      PARAMETER(          T3 = 41 ) !
      PARAMETER(       Force = 42 ) !
      PARAMETER(        CAMK = 43 ) !Fraction of activated CAM kinase ()
      PARAMETER(     rC1_RyR = 44 ) !
      PARAMETER(     rC2_RyR = 45 ) !
      PARAMETER(     rO3_RyR = 46 ) !
      PARAMETER(     rO4_RyR = 47 ) !
      PARAMETER(     rC5_RyR = 48 ) !
      PARAMETER(     rC6_RyR = 49 ) !
      PARAMETER(        yNa2 = 50 ) !Sustained Na inactivation gate ()


C *** Need to do this sometime !!!
c     Stim_tv   CONTROL(1)
c     Stim_id   18
c     Stim_tv_current_position AIO(1)
