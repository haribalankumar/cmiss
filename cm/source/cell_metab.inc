! cell_hmt.inc
!  Array elements for a coupled HMT -SPM  model.

C     Dependent variables
C *** SPM

      INTEGER ADP,AMP,ATP,B1,B1h,B2,B2h,Be,Behe,Ca,CaCal,Cae,
     '  Cal,CaTrop,CO2,CO2e,Cr,Glyc,h,h2phos,h2phose,hADP,hAMP,
     '  hATP,hco3,hco3e,
     '  he,hphos,hphose,K,Ke,Lac,Lace,MCT,MCTe,
     '  Mg,MgADP,MgATP,Mghphos,MgPCr,Na,Nae,PCr,R1,R1A,R2,R2A,Trop,
     '  Cl,Cle

C      PARAMETER(Vm=1) !membrain potential (units ??)
      PARAMETER(ADP=2) !ADP concentration (M)
      PARAMETER(AMP=3) !AMP concentration (M)
      PARAMETER(ATP=4) !ATP concentration (M)
      PARAMETER(B1=5)  !Buffer 1 concentration
      PARAMETER(B1h=6)
      PARAMETER(B2=7)
      PARAMETER(B2h=8)
      PARAMETER(Be=9)
      PARAMETER(Behe=10)
      PARAMETER(Ca=11)
      PARAMETER(CaCal=12)
      PARAMETER(Cae=13)
      PARAMETER(Cal=14)
      PARAMETER(CaTrop=15)
      PARAMETER(CO2=16)
      PARAMETER(CO2e=17)
      PARAMETER(Cr=18)
      PARAMETER(Glyc=19)
      PARAMETER(h=20)
      PARAMETER(h2phos=21)
      PARAMETER(h2phose=22)
      PARAMETER(hADP=23)
      PARAMETER(hAMP=24)
      PARAMETER(hATP=25)
      PARAMETER(hco3=26)
      PARAMETER(hco3e=27)
      PARAMETER(he=28)
      PARAMETER(hphos=29)
      PARAMETER(hphose=30)
      PARAMETER(K=31)
      PARAMETER(Ke=32)
      PARAMETER(Lac=33)
      PARAMETER(Lace=34)
      PARAMETER(MCT=35)
      PARAMETER(MCTe=36)
      PARAMETER(Mg=37)
      PARAMETER(MgADP=38)
      PARAMETER(MgATP=39)
      PARAMETER(Mghphos=40)
      PARAMETER(MgPCr=41)
      PARAMETER(Na=42)
      PARAMETER(Nae=43)
      PARAMETER(PCr=44)
      PARAMETER(R1=45)
      PARAMETER(R1A=46)
      PARAMETER(R2=47)
      PARAMETER(R2A=48)
      PARAMETER(Trop=49)
      PARAMETER(Cl=50)
      PARAMETER(Cle=51)

C *** PARAMETERS array
C     Material parameters
      INTEGER Ve,Vi,rho,k1,Voxphos,Koxphos,Vglyco,Kglyco,F,R,Temp,
     '  knapump1f,knapump1b,knapump2f,knapump2b,knapump3f,knapump3b,
     '  knapump4f,knapump4b,conapump,Kmct_h,Kmct_lac,kmct1,kmct1b,kmct2,
     '  kmct2b,V_NaPi_symport,K_NaPi_symport,V_Piefflux,K_Piefflux,
     '  P_co2trans,Vncef,Knce_nae,Knce_cae,Knhe,nnhe,Conhe,knhe1,
     '  knhe1b,knhe2,knhe2b,Knhe_Nae,Knhe_Na,Knhe_he,Knhe_h,Co_nbc,
     '  nnbc,Knbc_hco3,Knbc_Na,knbc_1,knbc_1b,knbc_2,n_ae,Kae,Coae,
     '  kae_1,kae_1b,kae_2,Kae_Cle,Kae_Cl,Kae_hco3e,Kae_hco3,Co_che,
     '  kche_1,
     '  kche_1b,kche_2,kche_2b,Kche_Cle,Kche_Cl,Kche_ohe,Kche_oh,
     '  kfh2co3_f,kfh2co3_r,kfh2co3e_f,kfh2co3e_r,Kkatp,Khatp,
     '  kfmgatp_r,Kkadp,Khadp,kfmgadp_r,Kkhphos,Khphos,
     '  KmgPCr,kfmgPCr_r,kfhatp_r,kfhadp_r,
     '  Kkamp,Khamp,kfhamp_r,kfh2phos_r,kintrinsicbuffer1_f,
     '  kintrinsicbuffer1_r,kintrinsicbuffer2_f,kintrinsicbuffer2_r,
     '  kexternalproteinbuffer_f,kexternalproteinbuffer_r,ktrop_f,
     '  ktrop_r,kcal_f,kcal_r,kfmghphos_r

      PARAMETER(Ve=1) !Extracellular volume
      PARAMETER(Vi=2) !Intracellular volume
      PARAMETER(rho=3) !??
      PARAMETER(k1=4)  !??
      PARAMETER(Voxphos=5) !??
      PARAMETER(Koxphos=6) !??
      PARAMETER(Vglyco=7) !??
      PARAMETER(Kglyco=8) !??
      PARAMETER(F=9) ! Faridays constant
      PARAMETER(R=10) ! Gas constant
      PARAMETER(Temp=11) ! Temparature in Kelvin
      PARAMETER(knapump1f=12) !
      PARAMETER(knapump1b=13) !
      PARAMETER(knapump2f=14) !
      PARAMETER(knapump2b=15) !
      PARAMETER(knapump3f=16) !
      PARAMETER(knapump3b=17) !
      PARAMETER(knapump4f=18) !
      PARAMETER(knapump4b=19) !
      PARAMETER(conapump=20) !
      PARAMETER(Kmct_h=21) !
      PARAMETER(Kmct_lac=22) !
      PARAMETER(kmct1=23) !
      PARAMETER(kmct1b=24) !
      PARAMETER(kmct2=25) !
      PARAMETER(kmct2b=26) !
      PARAMETER(V_NaPi_symport=27) !
      PARAMETER(K_NaPi_symport=28) !
      PARAMETER(V_Piefflux=29) !
      PARAMETER(K_Piefflux=30) !
      PARAMETER(P_co2trans=31) !
      PARAMETER(Vncef=32) !
      PARAMETER(Knce_nae=33) !
      PARAMETER(Knce_cae=34) !
      PARAMETER(Knhe=35) !
      PARAMETER(nnhe=36) !
      PARAMETER(Conhe=37) !
      PARAMETER(knhe1=38) !
      PARAMETER(knhe1b=39) !
      PARAMETER(knhe2=40) !
      PARAMETER(knhe2b=41) !
      PARAMETER(Knhe_Nae=42) !
      PARAMETER(Knhe_Na=43) !
      PARAMETER(Knhe_he=44) !
      PARAMETER(Knhe_h=45) !
      PARAMETER(Co_nbc=46)
      PARAMETER(nnbc=47)
      PARAMETER(Knbc_hco3=48)
      PARAMETER(Knbc_Na=49)
      PARAMETER(knbc_1=50)
      PARAMETER(knbc_1b=51)
      PARAMETER(knbc_2=52)
      PARAMETER(n_ae=53)
      PARAMETER(Kae=54)
      PARAMETER(Coae=55)
      PARAMETER(kae_1=56)
      PARAMETER(kae_1b=57)
      PARAMETER(kae_2=58)
      PARAMETER(Kae_Cle=59)
      PARAMETER(Kae_Cl=60)
      PARAMETER(Kae_hco3e=61)
      PARAMETER(Kae_hco3=62)
      PARAMETER(Co_che=63)
      PARAMETER(kche_1=64)
      PARAMETER(kche_1b=65)
      PARAMETER(kche_2=66)
      PARAMETER(kche_2b=67)
      PARAMETER(Kche_Cle=68)
      PARAMETER(Kche_Cl=69)
      PARAMETER(Kche_ohe=70)
      PARAMETER(Kche_oh=71)
      PARAMETER(kfh2co3_f=72)
      PARAMETER(kfh2co3_r=73)
      PARAMETER(kfh2co3e_f=74)
      PARAMETER(kfh2co3e_r=75)
      PARAMETER(Kkatp=76)
      PARAMETER(Khatp=77)
      PARAMETER(kfmgatp_r=78)
      PARAMETER(Kkadp=79)
      PARAMETER(Khadp=80)
      PARAMETER(kfmgadp_r=81)
      PARAMETER(Kkhphos=82)
      PARAMETER(kfmghphos_r=83)
      PARAMETER(Khphos=84)
      PARAMETER(KmgPCr=85)
      PARAMETER(kfmgPCr_r=86)
      PARAMETER(kfhatp_r=87)
      PARAMETER(kfhadp_r=88)
      PARAMETER(Kkamp=89)
      PARAMETER(Khamp=90)
      PARAMETER(kfhamp_r=91)
      PARAMETER(kfh2phos_r=92)
      PARAMETER(kintrinsicbuffer1_f=93)
      PARAMETER(kintrinsicbuffer1_r=94)
      PARAMETER(kintrinsicbuffer2_f=95)
      PARAMETER(kintrinsicbuffer2_r=96)
      PARAMETER(kexternalproteinbuffer_f=97)
      PARAMETER(kexternalproteinbuffer_r=98)
      PARAMETER(ktrop_f=99)
      PARAMETER(ktrop_r=100)
      PARAMETER(kcal_f=101)
      PARAMETER(kcal_r=102)


C *** PROTOCOL array

C     Model variable
      LOGICAL USE_ClOSED_SYS
      PARAMETER(USE_ClOSED_SYS=.TRUE.)

C     Derived variables
      INTEGER DVOXPHOS1,rv_ak,rvfhadp,rvfmgadp,rvfhamp,rvfhatp,
     '  rvfmgatp

      PARAMETER(DVOXPHOS1=1) !rate of oxidative phosphorylation
      PARAMETER(rv_ak=2) !
      PARAMETER(rvfhadp=3) !
      PARAMETER(rvfmgadp=4) !
      PARAMETER(rvfhamp=5) !
      PARAMETER(rvfhatp=6) !
      PARAMETER(rvfmgatp=7) !

