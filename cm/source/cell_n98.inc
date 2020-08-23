!cell_n98.inc

C     Dependent variables
      INTEGER Ki,Nai,Cai,CaDS,Caup,Carel,CaCALM,CaTROP,CaIND,
     '  rto,sto,xr1,xr2,xs,xACh1,xACh2,mNa,hNa,dCa,fCa,f2Ca,f2CaDS,
     '  ActivatorFract,ProductFract
      PARAMETER(Ki=2) !Intracellular [K+] (nmol.mm^-3)
      PARAMETER(Nai=3) !Intracellular [Na+] (nmol.mm^-3)
      PARAMETER(Cai=4) !Intracellular [Ca++] (nmol.mm^-3)
      PARAMETER(CaDS=5) !Diadic Space [Ca++] (nmol.mm^-3)
      PARAMETER(Caup=6) ![Ca++] uptake (nmol.mm^-3)
      PARAMETER(Carel=7) ![Ca++] released (nmol.mm^-3)
      PARAMETER(CaCALM=8) ![Ca++] bound to calmodulin (nmol.mm^-3)
      PARAMETER(CaTROP=9) ![Ca++] bound to troponin (nmol.mm^-3)
      PARAMETER(CaIND=10) ![Ca++] bound to an indicator (nmol.mm^-3)
      PARAMETER(rto=11) !Ito inactivation variable ()
      PARAMETER(sto=12) !Ito activation variable ()
      PARAMETER(xr1=13) !Fast Kr channel, fast gate variable ()
      PARAMETER(xr2=14) !Fast Kr channel, slow gate variable ()
      PARAMETER(xs=15) !Slow Ks channel gate ()
      PARAMETER(xACh1=16) !K ACh channel gate variable 1 ()
      PARAMETER(xACh2=17) !K ACh channel gate variable 2 ()
      PARAMETER(mNa=18) !INa activation gate ()
      PARAMETER(hNa=19) !INa inactivation gate ()
      PARAMETER(dCa=20) !ICaL activation gate ()
      PARAMETER(fCa=21) !ICaL voltage dependent inactivation gate ()
      PARAMETER(f2Ca=22) !ICaL Cai dependent inactivation gate ()
      PARAMETER(f2CaDS=23) !ICaLDS CaDS dependent inactivation gate ()
      PARAMETER(ActivatorFract=24) !Ca SR release activator fraction ()
      PARAMETER(ProductFract=25) !Ca SR release product fraction ()

C     Model parameters
      INTEGER TEMP,PKNa,PNaK,PCaK,PCaNa,PCaCa,PCa,gK1,gto,gKr1,
     '  gKr2,gKs,gKNa,gbK,gKATP,gKACh,gNa,gbNa,gpNa,gbCa,gK_stretch,
     '  gNa_stretch,gCa_stretch,gNs_stretch,gAn_stretch,INaKmax,
     '  INaCamax,EKATPrev,Esurf,ENs_stretch,EAn_stretch,alpha_rto,
     '  alpha_CaInact,alpha_CaDSInact,alpha_tr,alpha_rel,
     '  alpha_CaDSdecay,alpha_CALM,alpha_TROP,alpha_IND,alpha_SRleak,
     '  beta_CALM,beta_TROP,beta_IND,alpha_up,beta_up,KmK1,KmgKNa,
     '  KmATP,KmACh,KmKACh,KmCaInact,KmCaDSInact,KmK,KmNa,KmCa,KmCaDS,
     '  KmAChICaL,STEEP_IK1,STEEP_rto,STEEP_sto,STEEP_dCa,STEEP_fCa,
     '  SHIFT_IK1,SHIFT_rto,SHIFT_sto,SHIFT_IKr,SHIFT_alpha_xr1,
     '  SHIFT_beta_xr1,SHIFT_alpha_xr2,SHIFT_beta_xr2,SHIFT_alpha_xs,
     '  SHIFT_beta_xs,SHIFT_ACh1,SHIFT_ACh2,SHIFT_mNa,SHIFT_hNa,
     '  SHIFT_dCa,SHIFT_fCa,gtos,ICaLFract,
     '  INaCaFract,Ko,Nao,Cao,ATP,ACh,CALM,TROP,IND,dNaCa,SL,
     '  SLref,KcyCa,KSRCa,Kxcs,gamma_NaCa,gamma_SACSL,gamma_SRSL,
     '  gamma_TROPSL,epsilon_rto,epsilon_sto,epsilon_d,epsilon_f,
     '  epsilon_alpha_xr1,epsilon_beta_xr1,epsilon_alpha_xr2,
     '  epsilon_beta_xr2,epsilon_alpha_xs,epsilon_beta_xs,epsilon_sr,
     '  Vi,VDS,Vup,Vrel,Vo
      !PARAMETER(Cm=1) !Membrane capactiance per unit area (uF.mm^-2)
      !PARAMETER(Am=2) !Fibre surface to volume ratio (mm^-1)
      PARAMETER(TEMP=3) !Temperature (K)
      PARAMETER(PKNa=4) !Na+ permeability of IKs ()
      PARAMETER(PNaK=5) !K+ permeability of INa ()
      PARAMETER(PCaK=6) !ICaL permeability ratio of K+:Ca++ ()
      PARAMETER(PCaNa=7) !ICaL permeability ratio of Na+:Ca++ ()
      PARAMETER(PCa=8) !ICaL permeability to Ca++ (uA.mMol^-1.L)
      PARAMETER(gK1=9) !IK1 channel conductivity per unit area (mS.mm^-2)
      PARAMETER(gto=10) !Ito channel conductivity per unit area (mS.mm^-2)
      PARAMETER(gKr1=11) !IKr1 channel conductivity per unit area (mS.mm^-2)
      PARAMETER(gKr2=12) !IKr2 channel conductivity per unit area (mS.mm^-2)
      PARAMETER(gKs=13) !IKs channel conductivity per unit area (mS.mm^-2)
      PARAMETER(gKNa=14) !IKNa channel conductivity per unit area (mS.mm^-2)
      PARAMETER(gbK=15) !IbK channel conductivity per unit area (mS.mm^-2)
      PARAMETER(gKATP=16) !IKATP channel conductivity per unit area (mS.mm^-2)
      PARAMETER(gKACh=17) !IKACh channel conductivity per unit area (mS.mm^-2)
      PARAMETER(gNa=18) !INa channel conductivity  per unit area (mS.mm^-2)
      PARAMETER(gbNa=19) !IbNa channel conductivity per unit area (mS.mm^-2)
      PARAMETER(gpNa=20) !IpNa channel conductivity per unit area (mS.mm^-2)
      PARAMETER(gbCa=21) !IbCa channel conductivity per unit area (mS.mm^-2)
      PARAMETER(gK_stretch=22) !IK_stretch channel conductvity per unit area (mS.mm^-2)
      PARAMETER(gNa_stretch=23) !INa_stretch channel conductivity per unit area (mS.mm^-2)
      PARAMETER(gCa_stretch=24) !ICa_stretch channel conductivity per unit area (mS.mm^-2)
      PARAMETER(gNs_stretch=25) !INs_stretch channel conductivity per unit area (mS.mm^-2)
      PARAMETER(gAn_stretch=26) !IAn_stretch channel conductivity per unit area (mS.mm^-2)
      PARAMETER(INaKmax=27) !Na+/K+ exchange maximum current per unit area (uA.mm^-2)
      PARAMETER(INaCamax=28) !Na+/Ca++ exchange maximum current per unit area (uA.mm^-2.mmol^-4.L^4)
      PARAMETER(EKATPrev=29) !IKATP reversal potential (mV)
      PARAMETER(Esurf=30) !Surface charge potential on L-type Ca++ channel (mV)
      PARAMETER(ENs_stretch=31) !INs_stretch channel reversal potential (mV)
      PARAMETER(EAn_stretch=32) !IAn_stretch channel reversal potential (mV)
      PARAMETER(alpha_rto=33) !rto activation gate rate (ms^-1)
      PARAMETER(alpha_CaInact=34) !Rate of Cai inactivation of ICaL (ms^-1)
      PARAMETER(alpha_CaDSInact=35) !Rate of CaDS inactivation of ICaL (ms^-1)
      PARAMETER(alpha_tr=36) !SR translocation rate (nmol.mm^-3.ms^-1)
      PARAMETER(alpha_rel=37) !SR release rate (nmol.mm^-3.ms^-1)
      PARAMETER(alpha_CaDSdecay=38) !DS to intracellular Ca++ decay rate (nmol.mm^-3.ms^-1)
      PARAMETER(alpha_CALM=39) !Ca++-Calmodulin binding rate (ms^-1)
      PARAMETER(alpha_TROP=40) !Ca++-Troponin binding rate (ms^-1)
      PARAMETER(alpha_IND=41) !Ca++-Indicator binding rate (ms^-1)
      PARAMETER(alpha_SRleak=42) !SR leakage rate (ms^-1)
      PARAMETER(beta_CALM=43) !Ca++-Calmodulin unbinding rate (ms^-1)
      PARAMETER(beta_TROP=44) !Ca++-Troponin unbinding rate (ms^-1)
      PARAMETER(beta_IND=45) !Ca++-Indicator unbinding rate (ms^-1)
      PARAMETER(alpha_up=46) !SR uptake flux forward rate (nmol.mm^-3.ms^-1)
      PARAMETER(beta_up=47) !SR uptake flux reverse rate (nmol.mm^-3.ms^-1)
      PARAMETER(KmK1=48) !IK1 Ki Km (nmol.mm^-3)
      PARAMETER(KmgKNa=49) !IKNa Nai Km (nmol.mm^-3)
      PARAMETER(KmATP=50) !IKATP ATP Km (nmol.mm^-3)
      PARAMETER(KmACh=51) !IKACh ACh Km (nmol.mm^-3)
      PARAMETER(KmKACh=52) !IKACh K Km (nmol.mm^-3)
      PARAMETER(KmCaInact=53) !Cai Km for the inactivation of ICaL (nmol.mm^-3)
      PARAMETER(KmCaDSInact=54) !CaDS Km for the inactivation of ICaL (nmol.mm^-3)
      PARAMETER(KmK=55) !INaK Ki Km (nmol.mm^-3)
      PARAMETER(KmNa=56) !INaK Nai Km (nmol.mm^-3)
      PARAMETER(KmCa=57) !Cai Km for binding to SR release sites (nmol.mm^-3)
      PARAMETER(KmCaDS=58) !CaDS Km for binding to SR release sites (nmol.mm^-3)
      PARAMETER(KmAChICaL=59) !ICaL Km for ACh binding (nmol.mm^-3)
      PARAMETER(STEEP_IK1=60) !IK1 slope factor ()
      PARAMETER(STEEP_rto=61) !Ito activation gate slope factor (mV)
      PARAMETER(STEEP_sto=62) !Ito inactivation gate slope factor (mV)
      PARAMETER(STEEP_dCa=63) !ICaL activation gate slope factor (???)
      PARAMETER(STEEP_fCa=64) !ICaL voltage dependent inactivation gate slope factor (???)
      PARAMETER(SHIFT_IK1=65) !IK1 shift factor (mV)
      PARAMETER(SHIFT_rto=66) !rto shift factor (mV)
      PARAMETER(SHIFT_sto=67) !sto shift factor (mV)
      PARAMETER(SHIFT_IKr=68) !IKr1 and IKr2 shift factor (mV)
      PARAMETER(SHIFT_alpha_xr1=69) !IKr1 gate forward rate shift factor (mV)
      PARAMETER(SHIFT_beta_xr1=70) !IKr1 gate backward rate shift factor (mV)
      PARAMETER(SHIFT_alpha_xr2=71) !IKr2 gate forward rate shift factor (mV)
      PARAMETER(SHIFT_beta_xr2=72) !IKr2 gate backward rate shift factor (mV)
      PARAMETER(SHIFT_alpha_xs=73) !IKs gate forward rate shift factor (mV)
      PARAMETER(SHIFT_beta_xs=74) !IKs gate backward rate shift factor (mV)
      PARAMETER(SHIFT_ACh1=75) !IKACh first relaxation gate shift factor (mV)
      PARAMETER(SHIFT_ACh2=76) !IKACh second relaxation gate shift factor (mV)
      PARAMETER(SHIFT_mNa=77) !INa activation gate shift factor (mV)
      PARAMETER(SHIFT_hNa=78) !INa inactivation gate shift factor (mV)
      PARAMETER(SHIFT_dCa=79) !ICaL activation gate shift factor (mV)
      PARAMETER(SHIFT_fCa=80) !ICaL voltage dependent inactivation gate shift factor (mV)
      PARAMETER(gtos=81) !Sustained component of Ito as a fraction ()
      PARAMETER(ICaLFract=82) !Fraction of ICaL going into the diadic space ()
      PARAMETER(INaCaFract=83) !Fraction of INaCa going into the diadic space ()
      PARAMETER(Ko=84) !Extracellular K+ concentration (nmol.mm^-3)
      PARAMETER(Nao=85) !Extracellular Na+ concentration (nmol.mm^-3)
      PARAMETER(Cao=86) !Extracellular Ca++ concentration (nmol.mm^-3)
      PARAMETER(ATP=87) !Intracellular ATP concentration (nmol.mm^-3)
      PARAMETER(ACh=88) !Intracellular ACh concentration (nmol.mm^-3)
      PARAMETER(CALM=89) !Intracellular Calmodulin concentration (nmol.mm^-3)
      PARAMETER(TROP=90) !Intracellular Troponin concentration (nmol.mm^-3)
      PARAMETER(IND=91) !Intracellular Indicator concentration (nmol.mm^-3)
      PARAMETER(dNaCa=92) !INaCa denominator fraction (nmol^{-4}.mm^12)
      PARAMETER(SL=93) !Sarcomere length (um)
      PARAMETER(SLref=94) !Reference arcomere length (um)
      PARAMETER(KcyCa=95) !Hilgemann factor (nmol.mm^-3)
      PARAMETER(Kxcs=96) !Hilgemann factor ()
      PARAMETER(KSRCa=97) !Hilgemann (nmol.mm^-3)
      PARAMETER(gamma_NaCa=98) !Energy barrier position for Na+:Ca++ exchanger ()
      PARAMETER(gamma_SACSL=99) !??? (pum)
      PARAMETER(gamma_SRSL=100) !??? (pum)
      PARAMETER(gamma_TROPSL=101) !??? (pum)
      PARAMETER(epsilon_rto=102) !Ito activation gate speed factor ()
      PARAMETER(epsilon_sto=103) !Ito inactivation gate speed factor ()
      PARAMETER(epsilon_alpha_xr1=104) !IKr1 gate forward rate speed factor ()
      PARAMETER(epsilon_beta_xr1=105) !IKr1 gate backward rate speed factor ()
      PARAMETER(epsilon_alpha_xr2=106) !IKr2 gate forward rate speed factor ()
      PARAMETER(epsilon_beta_xr2=107) !IKr2 gate backward rate speed factor ()
      PARAMETER(epsilon_alpha_xs=108) !IKs gate forward rate speed factor ()
      PARAMETER(epsilon_beta_xs=109) !IKs gate backward rate speed factor ()
      PARAMETER(epsilon_d=110) !ICaL activation gate speed factor ()
      PARAMETER(epsilon_f=111) !ICaL inactivation gate speed factor ()
      PARAMETER(epsilon_sr=112) !SR speed factor ()
      PARAMETER(Vi=113) !Fraction of cell volume taken up by intracellular space ()
      PARAMETER(VDS=114) !Fraction of cell volume taken up by diadic space ()
      PARAMETER(Vup=115) !Fraction of cell volume taken up by SR uptake space ()
      PARAMETER(Vrel=116) !Fraction of cell volume taken up by SR release space ()
      PARAMETER(Vo=117) !Fraction of cell volume taken up by extracellular space ()

C     Protocol variables
      INTEGER PSTIMULUS,STIMULUS1_ON,STIMULUS1_OFF,
     '  STIMULUS1_MAG,STIMULUS2_ON,STIMULUS2_OFF,
     '  STIMULUS2_MAG,STIMULUS_FREQ_PERIOD,
     '  STIMULUS_FREQ_DURATION,STIMULUS_FREQ_MAG
      PARAMETER(PSTIMULUS=1) !Pseudo stimulus current (uA/mm^3)
      PARAMETER(STIMULUS1_ON=2) !Stimulus 1 time on (ms)
      PARAMETER(STIMULUS1_OFF=3) !Stimulus 1 time off (ms)
      PARAMETER(STIMULUS1_MAG=4) !Stimulus 1 magnitude (uA/mm^3)
      PARAMETER(STIMULUS2_ON=5) !Stimulus 2 time on (ms)
      PARAMETER(STIMULUS2_OFF=6) !Stimulus 2 time off (ms)
      PARAMETER(STIMULUS2_MAG=7) !Stimulus 2 magnitude (uA/mm^3)
      PARAMETER(STIMULUS_FREQ_PERIOD=8) !Stimulus frequency period (ms)
      PARAMETER(STIMULUS_FREQ_DURATION=9) !Stimulus frequency duration (ms)
      PARAMETER(STIMULUS_FREQ_MAG=10) !Stimulus frequency stimulus magnitude (uA/mm^3)

C     Control variables
      INTEGER RETURN_CURRENT
      PARAMETER(RETURN_CURRENT=1)

C     Model variable
      INTEGER USE_SAC
      PARAMETER(USE_SAC=1)

C     Derived variables
      INTEGER DIK1,DIto,DIKr1,DIKr2,DIKs,DIKNa,DIbK,DIKATP,DIKACh,DINa,
     '  DIbNa,DIpNa,DICaLK,DICaLNa,DICaLCa,DICaLKDS,DICaLNaDS,DICaLCaDS,
     '  DIbCa,DINaK,DINaCa,DINaCaDS,DIK_stretch,DINa_stretch,
     '  DICa_stretch,DINs_stretch,DIAn_stretch,Djup,Djtr,Djrel,Djleak,
     '  Djdecay,DIK,DIK_tot,DINa_tot,DICa_tot,DICaDS_tot,DIStim

      PARAMETER(DIK1=1) !Time-independent (background) K+ current (uA.mm^-2)
      PARAMETER(DIto=2) !Transient outward K+ current (uA.mm^-2)
      PARAMETER(DIKr1=3) !Time-dependent (delayed) fast K+ current (uA.mm^-2)
      PARAMETER(DIKr2=4) !Time-dependent (delayed) fast K+ current (uA.mm^-2)
      PARAMETER(DIKs=5) !Time-dependent (delayed) slow K+ current (uA.mm^-2)
      PARAMETER(DIKNa=6) !Na+-dependent K+ current (uA.mm^-2)
      PARAMETER(DIbK=7) !Background K+ current (uA.mm^-2)
      PARAMETER(DIKATP=8) !ATP-dependent K+ current (uA.mm^-2)
      PARAMETER(DIKACh=9) !ACh-dependent K+ current (uA.mm^-2)
      PARAMETER(DINa=10) !Fast Na+ current (uA.mm^-2)
      PARAMETER(DIbNa=11) !Background Na+ current (uA.mm^-2)
      PARAMETER(DIpNa=12) !Persistent Na+ current (uA.mm^-2)
      PARAMETER(DICaLK=13) !K+-dependent L-type Ca++ current (uA.mm^-2)
      PARAMETER(DICaLNa=14) !Na+-dependet L-type Ca++ current (uA.mm^-2)
      PARAMETER(DICaLCa=15) !Ca++-dependent L-type Ca++ current (uA.mm^-2)
      PARAMETER(DICaLKDS=16) !K+-dependent Diadic Space L-type Ca++ current (uA.mm^-2)
      PARAMETER(DICaLNaDS=17) !Na+-dependent Diadic Space L-type Ca++ current (uA.mm^-2)
      PARAMETER(DICaLCaDS=18) !Ca++-dependent Diadic Space L-type Ca++ current (uA.mm^-2)
      PARAMETER(DIbCa=19) !Background Ca++ current (uA.mm^-2)
      PARAMETER(DINaK=20) !Na+/K+ exchange current (uA.mm^-2)
      PARAMETER(DINaCa=21) !Na+/Ca++ exchange current (uA.mm^-2)
      PARAMETER(DINaCaDS=22) !Na+/Ca++ Diadic space exchange current (uA.mm^-2)
      PARAMETER(DIK_stretch=23) !K+ stretch activated current (uA.mm^-2)
      PARAMETER(DINa_stretch=24) !Na+ stretch activated current (uA.mm^-2)
      PARAMETER(DICa_stretch=25) !Ca++ stretch activated current (uA.mm^-2)
      PARAMETER(DINs_stretch=26) !Ns stretch activated current (uA.mm^-2)
      PARAMETER(DIAn_stretch=27) !An stretch activated current (uA.mm^-2)
      PARAMETER(Djup=28) !SR Ca++ uptake flux (nmol.mm^-3.ms^-1)
      PARAMETER(Djtr=29) !SR Ca++ translocation flux (nmol.mm^-3.ms^-1)
      PARAMETER(Djrel=30) !SR Ca++ release flux (nmol.mm^-3.ms^-1)
      PARAMETER(Djleak=31) !DS Ca++ leakage flux (nmol.mm^-3.ms^-1)
      PARAMETER(Djdecay=32) !DS Ca++ decay flux (nmol.mm^-3.ms^-1)
      PARAMETER(DIK=33) !IKr1+IKr2+IKs (uA.mm^-2)
      PARAMETER(DIK_tot=34) !Total K+ current (uA.mm^-2)
      PARAMETER(DINa_tot=35) !Total K+ current (uA.mm^-2)
      PARAMETER(DICa_tot=36) !Total K+ current (uA.mm^-2)
      PARAMETER(DICaDS_tot=37) !Total K+ current (uA.mm^-2)
      PARAMETER(DIStim=38) !Stimulus current
