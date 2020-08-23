C#### Comment: JRW Parameters
C###  Description:
C###    <html>
C###    <b>cell_jrw.inc<b> - Parameter indices for the
C###    Jafri-Rice-Winslow ionic current model.
C###    </html>

C *** State array
      INTEGER mNa,hNa,jNa,xK,Nai,Ki,Cai,CaSS,CaJSR,CaNSR,HTRPNCa,
     '  LTRPNCa,PO1,PO2,PC1,PC2,yL,C0,C1,C2,C3,C4,O,Cca0,Cca1,Cca2,
     '  Cca3,Cca4,Oca,Cao

C      CELL_RESERVED.INC
c      PARAMETER(Vm=1) Membrane potential (mV)
C      JRW-specific
      PARAMETER(mNa=2) !Na channel activation gate ()
      PARAMETER(hNa=3) !Na channel inactivation gate ()
      PARAMETER(jNa=4) !Na channel slow inactivation gate ()
      PARAMETER(xK=5) !K channel gate ()
      PARAMETER(Nai=6) !Intracellular [Na+] (nmol.mm^-3)=mM
      PARAMETER(Ki=7) !Intracellular [K+] (nmol.mm^-3)
      PARAMETER(Cai=8) !Intracelllular [Ca++] (nmol.mm^-3)
      PARAMETER(CaSS=9) !Sub-space [Ca++] (nmol.mm^-3)
      PARAMETER(CaJSR=10) !JSR [Ca++] (nmol.mm^-3)
      PARAMETER(CaNSR=11) !NSR [Ca++] (nmol.mm^-3)
      PARAMETER(HTRPNCa=12) ![Ca++] bound to high affinity toponin sites (nmol.mm^-3)
      PARAMETER(LTRPNCa=13) ![Ca++] bound to low affinity troponin sites (nmol.mm^-3)
      PARAMETER(Cao=14) !Extracellular [Ca++] (nmol.mm^-3)
      PARAMETER(PO1=15) !Fraction of RyR channels in open state 1 ()
      PARAMETER(PO2=16) !Fraction of RyR channels in open state 2 ()
      PARAMETER(PC1=17) !Fraction of RyR channels in closed state 1 ()
      PARAMETER(PC2=18) !Fraction of RyR channels in closed state 2 ()
      PARAMETER(yL=19) !ICaL voltage-dependent inactivation gate ()
      PARAMETER(C0=20) !Fraction of ICaL channels in closed state 0, mode normal ()
      PARAMETER(C1=21) !Fraction of ICaL channels in closed state 1, mode normal ()
      PARAMETER(C2=22) !Fraction of ICaL channels in closed state 2, mode normal ()
      PARAMETER(C3=23) !Fraction of ICaL channels in closed state 3, mode normal ()
      PARAMETER(C4=24) !Fraction of ICaL channels in closed state 4, mode normal ()
      PARAMETER(O=25) !Fraction of ICaL channels in open state, mode normal ()
      PARAMETER(Cca0=26) !Fraction of ICaL channels in closed state 0, mode Ca ()
      PARAMETER(Cca1=27) !Fraction of ICaL channels in closed state 1, mode Ca ()
      PARAMETER(Cca2=28) !Fraction of ICaL channels in closed state 2, mode Ca ()
      PARAMETER(Cca3=29) !Fraction of ICaL channels in closed state 3, mode Ca ()
      PARAMETER(Cca4=30) !Fraction of ICaL channels in closed state 4, mode Ca ()
      PARAMETER(Oca=31) !Fraction of ICaL channels in open state, mode Ca ()

C *** Derived array
      INTEGER DIstim,DINa,DICa,DIK,DIK1,DIKp,DINaCa,DINaK,DInsCa,DIpCa,
     '  DICaK,DICab,DINab,DInsNa,DInsK,DJrel,DJleak,DJup,DJtr,DJxfer,
     '  DJtrpn

      PARAMETER(DIstim=1) !Stimulus current (uA.mm^-2)
      PARAMETER(DINa=2) !Na+ current (uA.mm^-2)
      PARAMETER(DICa=3) !L-type channel Ca++ current (uA.mm^-2)
      PARAMETER(DIK=4)
      PARAMETER(DIK1=5) !time independent K+ current (uA.mm^-2)
      PARAMETER(DIKp=6) !plateau K+ current (uA.mm^-2)
      PARAMETER(DINaCa=7) !Na+/Ca++ exchanger current (uA.mm^-2)
      PARAMETER(DINaK=8) !Na+/K+ pump current (uA.mm^-2)
      PARAMETER(DInsCa=9) !non-specific Ca-activated total current (uA.mm^-2)
      PARAMETER(DIpCa=10) !sarcolemmal Ca++ pump current (uA.mm^-2)
      PARAMETER(DICaK=11) !L-type channel K+ current (uA.mm^-2)
      PARAMETER(DICab=12) !Ca++ background current (uA.mm^-2)
      PARAMETER(DINab=13) !Na+ background current (uA.mm^-2)
      PARAMETER(DInsNa=14) !non-specific Ca-activated Na+ current (uA.mm^-2)
      PARAMETER(DInsK=15) !non-specific Ca-activated K+ current (uA.mm^-2)
      PARAMETER(DJrel=16) !RyR channel flux (nmol.mm^-3.ms^-1)=(mM.ms^-1)
      PARAMETER(DJleak=17) !leak flux from NSR to myoplasm (nmol.mm^-3.ms^-1)
      PARAMETER(DJup=18) !SR Ca-ATPase pump flux (nmol.mm^-3.ms^-1)
      PARAMETER(DJtr=19) !transfer flux from NSR to JSR (nmol.mm^-3.ms^-1)
      PARAMETER(DJxfer=20) !transfer flux from SS to myoplasm (nmol.mm^-3.ms^-1)
      PARAMETER(DJtrpn=21) !Ca++ binding flux (nmol.mm^-3.ms^-1)


C *** Model array


C *** Control array
      INTEGER RETURN_CURRENTS

C     CELL_RESERVED.INC
c      PARAMETER(ODE=1) !==1 if evaluating ODEs in RHS Routine
C     JRW-specific
      PARAMETER(RETURN_CURRENTS=2)


C *** Additional Integer Input array


C *** Additional Integer Output array


C *** Parameter array
      INTEGER Nao,Ko,PNaK,CMDN_tot,km_CMDN,CSQN_tot,km_CSQN,v1,v2,v3,
     '  km_up,tau_xfer,tau_tr,kp_ltrpn,LTRPN_tot,km_ltrpn,kp_htrpn,
     '  km_htrpn,HTRPN_tot,GNa_max,GKp_max,kNaCa,eta,km_Na,km_Ca,ksat,
     '  INaK_max,km_Nai,km_Ko,P_ns,km_ns,Temp,ICap_max,km_Cap,GCab_max,
     '  GNab_max,PCa,PK,ICa_half,V_myo,Acap,V_SS,V_JSR,V_NSR,kap,kam,n,
     '  kbp,kbm,kcp,kcm,m,f,g,f_dash,g_dash,a,b,w

C     CELL_RESERVED.INC
c      PARAMETER(Cm=1) Membrane capacitance (uF.mm^-2)
c      PARAMETER(Am=1) Surface to volume ratio (mm^-1)
C     JRW-specific
      PARAMETER(PNaK=3) !Na+ permeability of K+ channels ()
      PARAMETER(CMDN_tot=4) !Total myoplasmic [calmodulin] (mM)
      PARAMETER(km_CMDN=5) !Ca++ half saturation constant for CMDN (mM)
      PARAMETER(CSQN_tot=6) !Total NSR [calsequestrin] (mM)
      PARAMETER(km_CSQN=7) !Ca++ half saturation constant for CSQN (mM)
      PARAMETER(v1=8) !Maximum RyR channel Ca++ permeability (ms^-1)
      PARAMETER(v2=9) !Ca++ leak rate constant from the NSR (ms^-1)
      PARAMETER(v3=10) !SR Ca-ATPase pump rate (mM.ms^-1)
      PARAMETER(km_up=11) !Half saturation constant for SR Ca-ATPase pump (mM)
      PARAMETER(tau_xfer=12) !Time constant for Ca++ transfer from subspace to myoplasm (ms)
      PARAMETER(tau_tr=13) !Time constant for Ca++ tansfer from NSR to JSR (ms)
      PARAMETER(kp_ltrpn=14) !Ca++ on rate constant for troponin low-affinity sites (mM^-1.ms^-1)
      PARAMETER(LTRPN_tot=15) !Total myoplasmic troponin low-affinity site concentration (mM)
      PARAMETER(km_ltrpn=16) !Ca++ off rate constant for troponin low-affinity sites (ms^-1)
      PARAMETER(kp_htrpn=17) !Ca++ on rate constant for troponin high-affinity sites (mM^-1.ms^-1)
      PARAMETER(km_htrpn=18) !Ca++ off rate constant for troponin high-affinity sites (ms^-1)
      PARAMETER(HTRPN_tot=19) !Total myoplasmic troponin high-affinity site concentration (mM)
      PARAMETER(GNa_max=20) !Maximum Na+ channel conductance (mS.mm^-2)
      PARAMETER(GKp_max=21) !Maximum plateau K+ channel conductance (mS.mm^-2)
      PARAMETER(kNaCa=22) !Scaling factor for Na+-Ca++ exchange (uA.mm^-2)
      PARAMETER(eta=23) !Controls voltage dependence of Na+-Ca++ excahange ()
      PARAMETER(km_Na=24) !Na+ half saturation constant for Na+-Ca++ exchange (mM)
      PARAMETER(km_Ca=25) !Ca++ half saturation constant for Na+-Ca++ exchange (mM)
      PARAMETER(ksat=26) !Na+-Ca++ exchange saturation at very negative potentials ()
      PARAMETER(INaK_max=27) !Maxium Na+-K+ pump current (uA.mm^-2)
      PARAMETER(km_Nai=28) !Na+ half saturation constant for Na+-K+ pump (mM)
      PARAMETER(km_Ko=29) !K+ half saturation constant for Na+-K+ pump (mM)
      PARAMETER(P_ns=30) !Nonspecific Ca++ activated current channel permeability (mm.ms^-1)
      PARAMETER(km_ns=31) !Ca++ half saturation constant for nonspecific Ca++ activated current (mM)
      PARAMETER(Temp=32) !Absolute temperature (K)
      PARAMETER(ICap_max=33) !Maximum sarcolemmal Ca++ pump current (uA.mm^-2)
      PARAMETER(km_Cap=34) !Ca++ half saturation constant for sarcolemmal Ca++ pump (mM)
      PARAMETER(GCab_max=35) !Maximum background Ca++ current conductance (mS.mm^-2)
      PARAMETER(GNab_max=36) !Maximum background Na+ current conductance (mS.mm^-2)
      PARAMETER(PCa=37) !L-type Ca++ channel permeability to Ca++ (mm.ms^-1)
      PARAMETER(PK=38) !L-type Ca++ channel permeability to K+ (mm.ms^-1)
      PARAMETER(ICa_half=39) !ICa level which reduces PK by half (uA.mm^-2)
      PARAMETER(V_myo=40) !Myoplasmic volume (mm^3)
      PARAMETER(Acap=41) !Capacitive membrane area (mm^2)
      PARAMETER(V_SS=42) !Calcium sub-space volume (mm^3)
      PARAMETER(V_JSR=43) !JSR volume (mm^3)
      PARAMETER(V_NSR=44) !NSR volume (mm^3)
      PARAMETER(kap=45) !RyR PC1-PO1 rate constant (mM^-4.ms^-1)
      PARAMETER(kam=46) !RyR PO1-PC1 rate constnat (ms^-1)
      PARAMETER(n=47) !RyR Ca++ cooperativity parameter PC1-PO2 ()
      PARAMETER(kbp=48) !RyR PO1-PO2 rate constant (mM^-3.ms^-1)
      PARAMETER(kbm=49) !RyR PO2-PO1 rate constant (ms^-1)
      PARAMETER(kcp=50) !RyR PO1-PC2 rate constant (ms^-1)
      PARAMETER(kcm=51) !RyR PC2-PO1 rate constant (ms^-1)
      PARAMETER(m=52) !RyR Ca++ cooperativity parameter PO1-PO2 ()
      PARAMETER(f=53) !L-type Ca++ channel transition rate into open state, mode normal (ms^-1)
      PARAMETER(g=54) !L-tpye Ca++ channel transition rate out of open state, mode normal (ms^-1)
      PARAMETER(f_dash=55) !L-type Ca++ channel transition rate into open state, mode Ca (ms^-1)
      PARAMETER(g_dash=56) !L-tpye Ca++ channel transition rate out of open state, mode Ca (ms^-1)
      PARAMETER(a=57) !L-type Ca++ channel mode transition parameter ()
      PARAMETER(b=58) !L-type Ca++ channel mode transition parameter ()
      PARAMETER(w=59) !L-type Ca++ channel mode transition parameter (ms^-1)
      PARAMETER(Nao=60) !Extracellular [Na+] (mM)
      PARAMETER(Ko=61) !Extracellular [K+] (mM)

C *** Protocol array

C     CELL_RESERVED.INC
c      PARAMETER(      PseudoIs =  1 ) Pseudo stimulus current (uA/mm^3)
c      PARAMETER(      Is1start =  2 ) Stimulus 1 time on (ms)
c      PARAMETER(       Is1stop =  3 ) Stimulus 1 time off (ms)
c      PARAMETER(    Is1current =  4 ) Stimulus 1 magnitude (uA/mm^3)
c      PARAMETER(      Is2start =  5 ) Stimulus 2 time on (ms)
c      PARAMETER(       Is2stop =  6 ) Stimulus 2 time off (ms)
c      PARAMETER(    Is2current =  7 ) Stimulus 2 magnitude (uA/mm^3)
c      PARAMETER(  IsFreqPeriod =  8 ) Stimulus frequency period (ms)
c      PARAMETER(IsFreqDuration =  9 ) Stimulus frequency duration (ms)
c      PARAMETER(     IsFreqMag = 10 ) Stimulus frequency stimulus magnitude (uA/mm^3)
C     JRW-specific
