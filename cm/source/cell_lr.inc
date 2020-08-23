! cell_lr.inc
!  Array elements for the Luo-Rudy ionic current model.

C *** CONTROL array
C     DPN 14 February 2000 - unused yet
c      INTEGER NUM_APPLIED_STIMULII
c      PARAMETER(NUM_APPLIED_STIMULII=1)

      INTEGER RETURN_CURRENTS

      !PARAMETER(           ODE=1) !in cell_reserved.inc
      PARAMETER(RETURN_CURRENTS=1) !==1 if currents should be returned
                                   !in the derived array

C *** ARO - additional real out array
      INTEGER alpha_m,beta_m,alpha_h,beta_h,alpha_j,beta_j,alpha_d,
     '  beta_d,alpha_f,beta_f,alpha_x,beta_x,alpha_K1,beta_K1,
     '  IICa,IICaK,IICaNa,IINa,IIK,IIK1,fca

      PARAMETER(     alpha_m =  1 ) !
      PARAMETER(      beta_m =  2 ) !
      PARAMETER(     alpha_h =  3 ) !
      PARAMETER(      beta_h =  4 ) !
      PARAMETER(     alpha_j =  5 ) !
      PARAMETER(      beta_j =  6 ) !
      PARAMETER(     alpha_d =  7 ) !
      PARAMETER(      beta_d =  8 ) !
      PARAMETER(     alpha_f =  9 ) !
      PARAMETER(      beta_f = 10 ) !
      PARAMETER(     alpha_x = 11 ) !
      PARAMETER(      beta_x = 12 ) !
      PARAMETER(    alpha_K1 = 13 ) !
      PARAMETER(     beta_K1 = 14 ) !
      PARAMETER(        IICa = 15 ) !
      PARAMETER(       IICaK = 16 ) !
      PARAMETER(      IICaNa = 17 ) !
      PARAMETER(        IINa = 18 ) !
      PARAMETER(         IIK = 19 ) !
      PARAMETER(        IIK1 = 20 ) !
      PARAMETER(         fca = 21 ) !

C *** PARAMETERS array
      INTEGER TIMESCALE,TSTART,TEND,DT,TABT,TPS,TP,GGCab,KmCa_L,PCa,
     '  gCao,kNaCa,KmNa,KmCa_NaCa,eta,ksat,KmpCa,IIpCa,KmnsCa,PnsCa,
     '  LR_Cm,Acap,Temp,VolMyo,CSQNCSQN,KmCSQN,Kmrel,Tau_on,Tau_off,
     '  Caith,
     '  GGrel_,VolJSR,Tau_tr,Kmup,CaNSRCaNSR,IIup,VolNSR,gCai,TRPNTRPN,
     '  KmTRPN,CMDNCMDN,KmCMDN,GGKp,PK,gNai,gKi,gNao,gKo,PNaK,IINaK,
     '  KmNai,KmKo,GGNab,GGNa,PNa,Istim,Ko,Nao,tps2,numStim,dtp,RTONF,
     '  Cai_on,Cai2
      REAL*8 Farad

      PARAMETER(   TIMESCALE = 64 ) !
      PARAMETER(      TSTART = 65 ) !
      PARAMETER(        TEND =  3 ) !
      PARAMETER(          DT =  4 ) !
      PARAMETER(        TABT =  5 ) !
      PARAMETER(         TPS =  6 ) !
      PARAMETER(          TP =  7 ) !
      PARAMETER(       GGCab =  8 ) !
      PARAMETER(      KmCa_L =  9 ) !
      PARAMETER(         PCa = 10 ) !
      PARAMETER(        gCao = 11 ) !
      PARAMETER(       kNaCa = 12 ) !
      PARAMETER(        KmNa = 13 ) !
      PARAMETER(   KmCa_NaCa = 14 ) !
      PARAMETER(         eta = 15 ) !
      PARAMETER(        ksat = 16 ) !
      PARAMETER(       KmpCa = 17 ) !
      PARAMETER(       IIpCa = 18 ) !
      PARAMETER(      KmnsCa = 19 ) !
      PARAMETER(       PnsCa = 20 ) !
      PARAMETER(       LR_Cm = 21 ) !
      PARAMETER(        Acap = 22 ) !
      PARAMETER(        Temp = 23 ) !
      PARAMETER(      VolMyo = 24 ) !
      PARAMETER(    CSQNCSQN = 25 ) !
      PARAMETER(      KmCSQN = 26 ) !
      PARAMETER(       Kmrel = 27 ) !
      PARAMETER(      Tau_on = 28 ) !
      PARAMETER(     Tau_off = 29 ) !
      PARAMETER(       Caith = 30 ) !
      PARAMETER(      GGrel_ = 31 ) !
      PARAMETER(      VolJSR = 32 ) !
      PARAMETER(      Tau_tr = 33 ) !
      PARAMETER(        Kmup = 34 ) !
      PARAMETER(  CaNSRCaNSR = 35 ) !
      PARAMETER(        IIup = 36 ) !
      PARAMETER(      VolNSR = 37 ) !
      PARAMETER(        gCai = 38 ) !
      PARAMETER(    TRPNTRPN = 39 ) !
      PARAMETER(      KmTRPN = 40 ) !
      PARAMETER(    CMDNCMDN = 41 ) !
      PARAMETER(      KmCMDN = 42 ) !
      PARAMETER(        GGKp = 43 ) !
      PARAMETER(          PK = 44 ) !
      PARAMETER(        gNai = 45 ) !
      PARAMETER(         gKi = 46 ) !
      PARAMETER(        gNao = 47 ) !
      PARAMETER(         gKo = 48 ) !
      PARAMETER(        PNaK = 49 ) !
      PARAMETER(       IINaK = 50 ) !
      PARAMETER(       KmNai = 51 ) !
      PARAMETER(        KmKo = 52 ) !
      PARAMETER(       GGNab = 53 ) !
      PARAMETER(        GGNa = 54 ) !
      PARAMETER(         PNa = 55 ) !
      PARAMETER(       Istim = 56 ) !
      PARAMETER(          Ko = 57 ) !
      PARAMETER(         Nao = 58 ) !
      PARAMETER(        tps2 = 58 ) !
      PARAMETER(     numStim = 59 ) !
      PARAMETER(         dtp = 60 ) !
      PARAMETER(       RTONF = 61 ) !
      PARAMETER(      Cai_on = 62 ) !
      PARAMETER(        Cai2 = 63 ) !
      PARAMETER(       Farad = 96.5d3 ) !C/mol

C *** DERIVED array
      INTEGER Stimulus,INa,ICaLCa,ICaLK,ICaLNa,ICaL,IK,IK1,IKp,INaCa,
     '  INaK,InsK,InsNa,InsCa,IpCa,ICab,INab,Iv,Irel,Ileak,Iup,I_Na,
     '  I_Ca,I_K,F_NSR,F_JSR,Itr

      PARAMETER(    Stimulus =  1 ) !
      PARAMETER(         INa =  2 ) !
      PARAMETER(      ICaLCa =  3 ) !
      PARAMETER(       ICaLK =  4 ) !
      PARAMETER(      ICaLNa =  5 ) !
      PARAMETER(        ICaL =  6 ) !
      PARAMETER(          IK =  7 ) !
      PARAMETER(         IK1 =  8 ) !
      PARAMETER(         IKp =  9 ) !
      PARAMETER(       INaCa = 10 ) !
      PARAMETER(        INaK = 11 ) !
      PARAMETER(        InsK = 12 ) !
      PARAMETER(       InsNa = 13 ) !
      PARAMETER(       InsCa = 14 ) !
      PARAMETER(        IpCa = 15 ) !
      PARAMETER(        ICab = 16 ) !
      PARAMETER(        INab = 17 ) !
      PARAMETER(          Iv = 18 ) !
      PARAMETER(        Irel = 19 ) !
      PARAMETER(       Ileak = 20 ) !
      PARAMETER(         Iup = 21 ) !
      PARAMETER(        I_Na = 22 ) !
      PARAMETER(        I_Ca = 23 ) !
      PARAMETER(         I_K = 24 ) !
      PARAMETER(       F_NSR = 25 ) !
      PARAMETER(       F_JSR = 26 ) !
      PARAMETER(         Itr = 27 ) !

C *** Y and DY arrays
      INTEGER Nai,Cai,x,d,f,h,m,j,Ki,CaNSR,CaJSR,Cab,Cao

c     PARAMETER(          Vm =  1 ) !given in cell_reserved.inc
      PARAMETER(         Nai =  2 ) !
      PARAMETER(         Cai =  3 ) !
      PARAMETER(           x =  4 ) !
      PARAMETER(           d =  5 ) !
      PARAMETER(           f =  6 ) !
      PARAMETER(           h =  7 ) !
      PARAMETER(           m =  8 ) !
      PARAMETER(           j =  9 ) !
      PARAMETER(          Ki = 10 ) !
      PARAMETER(       CaNSR = 11 ) !
      PARAMETER(       CaJSR = 12 ) !
      PARAMETER(         Cab = 13 ) !
      PARAMETER(         Cao = 14 ) !

C *** Need to do this sometime !!!
c     Stim_tv   CONTROL(1)
c     Stim_id   18
c     Stim_tv_current_position AIO(1)
