! cell_hmt.inc
!  Array elements for a coupled HMT -SPM  model.

C     Dependent variables
C *** SPM

      INTEGER Cai,Caf,Cas,CaCSQ,Cas1,Cas2
      !PARAMETER(Vm=1) [reserved]
      PARAMETER(Cai=2) !Quiesent free [Ca++], cytosol (M)
      PARAMETER(Caf=3) !free Ca++ concentration in the fuzzy space
      PARAMETER(Cas=4) !free Ca++ concentration in the SR
      PARAMETER(CaCSQ=5)!Ca++ concentration bound to CSQ
      PARAMETER(Cas1=6)!fractional concentration of bound activating site
      PARAMETER(Cas2=7)!fractional concentration of bound inactivating site


C *** PARAMETERS array
C **** SPM
      INTEGER Rt,Vf,Vs,Vc,K1,K2,kon1,kon2,n1,n2,ks,kl,period,VmaxNaCaX,
     '  KmNaCaX,nNaCaX,Vmaxs,Kms,ns,kf,Cao,Bmaxc,Kbc,dyei,Kbdye,Bmaxf1,
     '  kbf1,Bmaxf2,Kbf2,Bmaxs,Kbs,kons,koffs,koff1,koff2
      PARAMETER(Rt=1) !SRRC concentration (M)
      PARAMETER(Vf=2) !Volume fraction, fuzzy space
      PARAMETER(Vs=3) !Volume fraction SR
      PARAMETER(Vc=4) !Volume fraction, cytosol
      PARAMETER(K1=5) !koff1/kon1 SRRC, fast site (M)
      PARAMETER(K2=6) !!koff2/kon2 SRRC, slow site (M)
      PARAMETER(kon1=7) !Binding rate constant, SRRC, slow site (M.s)-1
      PARAMETER(kon2=8) !Binding rate constant, SRRC, fast site (M.s)-1
      PARAMETER(n1=9) !Hill coeffiecent, fast site
      PARAMETER(n2=10) !Hill coeffiecent, slow site
      PARAMETER(ks=11) !Release rate constant, SRRC (s-1)
      PARAMETER(kl=12) !Net L-type rate constant (s-1)
      PARAMETER(period=13) !Period of Ca++ influx
      PARAMETER(VmaxNaCaX=14) !Maximum rate of Na+/Ca++ exchange (M.s-1)
      PARAMETER(KmNaCaX=15) !Dissociation constant, Na+/Ca++ exchange (M)
      PARAMETER(nNaCaX=16) !Hill coefficent, Na+/Ca++ exchange
      PARAMETER(Vmaxs=17) !Maximum uptake rate, SR Ca++-ATPase
      PARAMETER(Kms=18) !Dissociation constant, SR Ca++-ATPase
      PARAMETER(ns=19) !Hill coeffient, SR Ca++-ATPase
      PARAMETER(kf=20) !Diffusion rate constant, fuzzy space/cytosol (s-1)
      PARAMETER(Cao=21) !Extracellular free Ca++ concentration (M)
      PARAMETER(Bmaxc=22) !Maximum free fast buffers, cytosol (M)
      PARAMETER(Kbc=23) !Buffer dissociation constant (fast), cytosol (M)
      PARAMETER(dyei=24) !Concentration of dye indicator, cytosol (M)
      PARAMETER(Kbdye=25) !Dye indicator dissociation constant, cytosol (M)
      PARAMETER(Bmaxf1=26) !Maximum free fast buffers, cytosol (M)
      PARAMETER(kbf1=27) ! Buffer dissociation constant, fuzzy space 1 (M)
      PARAMETER(Bmaxf2=28) !Maximum free buffer, fuzzy space site 1 (M)
      PARAMETER(Kbf2=29) !Buffer dissociation constant, space site 2 (M)
      PARAMETER(Bmaxs=30) !Maximum free buffer, fuzzy space site 2 (M)
      PARAMETER(Kbs=31) !Buffer dissociation constant, SR (M)
      PARAMETER(kons=32) !Buffer binding rate constant (M-1.s-1)


C *** PROTOCOL array
      INTEGER IsFreqStart,IsFreqStop
      PARAMETER(IsFreqStart=11) ! Start time for frequency stimuluation (ms)
      PARAMETER(IsFreqStop=12) ! Stop time for frequency stimuluation (ms)

C *** AII array

C *** AIO array

C *** ARI array

C *** ARO array
