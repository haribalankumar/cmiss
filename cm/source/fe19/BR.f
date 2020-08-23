      SUBROUTINE BR(TIME,Y,DY,CONTROL,MODEL,SIZES,VARIANT,DERIVED,
     '  PARAM,PROTOCOL,AII,AIO,ARI,ARO,ERR)

C#### Subroutine: BR
C###  Description:
C###    BR is the Beeler-Reuter model, with the option of sodium
C###    kinetics determined by Ebihara-Johnson or Drouhard-Roberge
C###    according to KTYP33
C**** (re)Written by Martin Buist, 7 May 1999

C***  SIZES(1) is 8 equations (Vm,m,h,j,d,f1,x1,Cai).
C***  SIZES(2) is not used.
C***  SIZES(3) is not used.
C***  SIZES(4) is not used.
C***  SIZES(5) is 7 parameters.
C***  PARAM(1) is the membrane capacitance Cm.
C***  PARAM(2) is the surface to volume ration Am.
C***  PARAM(3) is resting potential.
C***  PARAM(4) is sodium reversal potential.
C***  PARAM(5) is sodium conductance.
C***  PARAM(6) is steady state sodium conductance.
C***  PARAM(7) is slow current conductance.
C***  SIZES(6) is 4 protocols.
C***  PROTOCOL(1) is pseudo stimulus current.
C***  PROTOCOL(4) is applied stimulus current.
C***  SIZES(7) is not used.
C***  SIZES(8) is not used.
C***  SIZES(9) is not used.
C***  SIZES(10) is not used.
C***  SIZES(11) is not used.

      IMPLICIT NONE

      INCLUDE 'cell_br.inc'
      INCLUDE 'cell_reserved.inc'
      INCLUDE 'ktyp30.cmn'

!     Parameter List
      INTEGER CONTROL(*),MODEL(*),VARIANT,AII(*),AIO(*),ERR,SIZES(11)
      REAL*8 TIME(*),Y(*),DY(*),DERIVED(*),PARAM(*),PROTOCOL(*),ARI(*),
     '  ARO(*)
C      INTEGER SIZES(10)
C      INTEGER CONTROL(SIZES(NUM_CONTROL)),MODEL(SIZES(NUM_MODEL)),
C     '  VARIANT,AII(SIZES(NUM_AII)),AIO(SIZES(NUM_AIO)),ERR
C      REAL*8 T,Y(SIZES(NUM_EQN)),DY(SIZES(NUM_EQN)),
C     '  DERIVED(SIZES(NUM_DERIVED)),PARAM(SIZES(NUM_PARAM)),
C     '  PROTOCOL(SIZES(NUM_PROTOCOL)),ARI(SIZES(NUM_ARI)),
C     '  ARO(SIZES(NUM_ARO))
!     Local Variables
      REAL*8 Ad,Af,Ah,Aj,Alm,Ax1,Bd,Bf,Bh,Bj,Bm,Bx1,IK1,INa,Is,Ix1,V,Vs
      REAL*8 alpha,beta,gamma

      INa=0.d0
      IF(KTYP33.EQ.1) THEN
        !Beeler-Reuter sodium model
        INa = PARAM(GNa)*Y(m)*Y(m)*Y(m)*Y(h)*Y(j)+PARAM(GNaC)
      ELSE IF(KTYP33.EQ.2.OR.KTYP33.GE.3) THEN
        !Ebihara-Johnson or Drouhard-Roberge sodium model
        INa = PARAM(GNa)*Y(m)*Y(m)*Y(m)*Y(h)
      ENDIF !sodium model
      INa = INa*(Y(Vm)-PARAM(VNa))

      Vs = -82.3d0 - 13.0287d0*DLOG(Y(Cai)*0.001d0)
      Is = PARAM(Gs)*Y(d)*Y(f1)*(Y(Vm)-Vs)

C     1.d-2 correction factor in the next two equations is for the
C     change from cm^2 to mm^2
C SNHs 1/2/00 This correction was removed but has been put back in
C so that all currents are in uA/mm^2.

      Ix1 = 0.8d0*Y(x1)*(DEXP(0.04d0*(Y(Vm)+77.d0))-1.d0)/
     '  (DEXP(0.04d0*(Y(Vm)+35.d0)))*1.d-2
      IK1 = 0.35d0 * (4.d0*(DEXP(0.04d0*(Y(Vm)+85.d0))-1.d0) /
     '  (DEXP(0.08d0*(Y(Vm)+53.d0))+DEXP(0.04d0*(Y(Vm)+53.d0))) +
     '  0.2d0*(Y(Vm)+23.d0)/(1.d0-DEXP(-0.04d0*(Y(Vm)+23.d0))))*1.d-2

C The total current is now in uA/mm^2.
C *** DPN 15 March 2000 - fixing up units
C      DY(Vm)=(INa+Is+Ix1+IK1)*1.0d-3 !correction factor for [mA]
      DY(Vm)=INa+Is+Ix1+IK1
C SNHe
C      V=Y(Vm)+DY(Vm)*DT !prediction for V
      V=Y(Vm)

      IF(KTYP33.EQ.1) THEN !Beeler-Reuter sodium model
        Alm = -1.d0*(V+47.d0)/(DEXP(-0.1d0*(V+47.d0))-1.d0)
        Bm = 40.d0*DEXP(-0.056d0*(V+72.d0))
        Ah = 0.126d0*DEXP(-0.25d0*(V+77.d0))
        Bh = 1.7d0/(DEXP(-.082d0*(V+22.5d0))+1.d0)
        Aj = 0.055d0*DEXP(-0.25d0*(V+78.d0))/
     '    (DEXP(-0.2d0*(V+78.d0))+1.d0)
        Bj = 0.3d0/(DEXP(-0.1d0*(V+32.d0))+1.d0)
      ELSE IF(KTYP33.EQ.2) THEN !Ebihara-Johnson sodium model
        Alm = (0.32d0*(V+47.13d0))/(1.d0-DEXP(-V-47.13d0))
        Bm = 0.08d0*DEXP(-V/11.d0)
        IF(V.GE.-40.d0) THEN
          Ah = 0.d0
          Bh = 1.d0/(0.13d0*(DEXP((V+10.66d0)/(-11.1d0))+1.d0))
        ELSE
          Ah = 0.135d0*DEXP((-80.d0-V)/6.8d0)
          Bh = 3.56d0*DEXP(0.079d0*V)+3.1d5*DEXP(0.35d0*V)
        ENDIF
        Aj = 0.d0
        Bj = 0.d0
      ELSE IF(KTYP33.EQ.3) THEN !Drouhard-Roberge sodium model
        Alm = 0.9d0*(V+42.65d0)/(1.d0-DEXP(-0.22d0*(V+42.65d0)))
        Bm = 1.437d0*DEXP(-0.085d0*(V+39.75d0))
        Ah = 0.1d0*DEXP(-0.193d0*(V+79.65d0))
        Bh = 1.7d0/(1.d0+DEXP(-0.095d0*(V+20.5d0)))
        Aj = 0.d0
        Bj = 0.d0
      ELSE IF(KTYP33.EQ.4) THEN !Modified Drouhard-Roberge Na model
        !Alpha_m and Beta_m
        IF(V.LE.100.0d0) THEN
          Alm=0.9d0*((V+42.65d0)/(1.0d0-DEXP((-0.22d0*V)-9.383d0)))
        ELSE
          Alm=890.943789d0*((DEXP((0.0486479d0*V)-4.8647916d0))/
     '      (1.0d0+(5.93962526d0*DEXP((0.0486479d0*V)-4.8647916d0))))
        ENDIF
        IF(V.LE.-85.0d0) THEN
          Bm=100.0d0/(1.0d0+(0.4864082d0*DEXP((0.2597504d0*V)+
     '      22.0787804d0)))
        ELSE
          Bm=1.437d0*DEXP((-0.085d0*V)-3.37875d0)
        ENDIF

        !Alpha_h and Beta_h
        IF(V.LE.-90.0d0) THEN
          Ah=-12.0662845d0-(0.1422598d0*V)
        ELSE
          Ah=0.1d0*DEXP((-0.193d0*V)-15.37245d0)
        ENDIF
        Bh=1.7d0/(1.0d0+DEXP((-0.095d0*V)-1.9475d0))
        Aj = 0.d0
        Bj = 0.d0
      ENDIF !sodium model
      DY(m)=Alm*(1.d0-Y(m))-Bm*Y(m)
      DY(h)=Ah*(1.d0-Y(h))-Bh*Y(h)
      DY(j)=Aj*(1.d0-Y(j))-Bj*Y(j)

C DPN 15 March 2000 - remove repeated code
c      Vs = -82.3d0 - 13.0287d0*DLOG(Y(Cai))
c      Is = PARAM(Gs)*Y(d)*Y(f1)*(V-Vs)

C SNHs 1/2/00 The current must be converted to uA/cm^2 for this equation.
C      IF(KTYP33.EQ.4) THEN
C        !d[Ca]i/dt
C        IF(V.LE.200.0d0) THEN
C          !DY(Cai)=-1.d-7*Is + 0.07d0*(1.d-7-Y(Cai))
C          DY(Cai)=-1.d-5*Is + 0.07d0*(1.d-7-Y(Cai))
C        ELSE
C          DY(Cai)=0.0d0
C        ENDIF
C      ELSE
C        !DY(Cai)=-1.d-7*Is + 0.07d0*(1.d-7-Y(Cai))
C        DY(Cai)=-1.d-5*Is + 0.07d0*(1.d-7-Y(Cai))
C      ENDIF
C MLB 10August2000, concentration in mM not M
      IF(KTYP33.EQ.4) THEN
        !d[Ca]i/dt
        IF(V.LE.200.0d0) THEN
          DY(Cai)=-1.d-2*Is + 0.07d0*(1.d-4-Y(Cai))
        ELSE
          DY(Cai)=0.0d0
        ENDIF
      ELSE
        DY(Cai)=-1.d-2*Is + 0.07d0*(1.d-4-Y(Cai))
      ENDIF
C SNHe

      Ad = 0.095d0*DEXP(-0.01d0*(V-5.d0))/
     '  (1.d0+DEXP(-0.072d0*(V-5.d0)))
      Bd = 0.07d0*DEXP(-(V+44.d0)/59.d0)/
     '  (1.d0+DEXP(0.05d0*(V+44.d0)))
      Af = 0.012d0*DEXP(-0.008d0*(V+28.d0))/
     '  (1.d0+DEXP(0.15d0*(V+28.d0)))
      Bf = 0.0065d0*DEXP(-0.02d0*(V+30.d0))/
     '  (1.d0+DEXP(-0.2d0*(V+30.d0)))
C MLB it's meant to be divided by, not multiplied!
C      DY(d)=KTYP39R*Ad*(1.d0-Y(d))-Bd*Y(d)
C      DY(f1)=KTYP39R*Af*(1.d0-Y(f1))-Bf*Y(f1)
      DY(d)=(1.0d0/KTYP39R)*Ad*(1.d0-Y(d))-Bd*Y(d)
      DY(f1)=(1.0d0/KTYP39R)*Af*(1.d0-Y(f1))-Bf*Y(f1)

      IF(KTYP33.EQ.4) THEN
        !Alpha_x1 and Beta_x1
        IF(V.LE.400.0d0) THEN
          Ax1=0.0005d0*(DEXP((0.083d0*V)+4.150d0)/(DEXP((0.057d0*V)+
     '      2.85d0)+1.0d0))
        ELSE
          Ax1=151.7994692d0*(DEXP((0.0654679d0*V)-26.1871448d0)/
     '      (1.0d0+(1.5179947d0*DEXP((0.0654679d0*V)-26.1871448d0))))
        ENDIF
        Bx1=0.0013d0*(DEXP((-0.06d0*V)-1.2d0)/((DEXP((-0.04d0*V)
     '    -0.8d0))+1.0d0))
      ELSE
C AJPs 30/11/99 Ax1 has an error - no negative in the first exponent.
C Also put formulae to reflect original paper.
C        Ax1 = 5.d-4*DEXP(-(V+50.d0)/12.1d0)/(1.d0+DEXP((V+50.d0)/
C     '    17.5d0))
        Ax1 = 5.d-4*DEXP(0.083d0*(V+50.0d0))/
     '    (1.d0+DEXP(0.057d0*(V+50.d0)))
C AJPe
        Bx1 = 0.0013d0*DEXP(-0.06d0*(V+20.d0))/
     '    (1.d0+DEXP(-0.04d0*(V+20.d0)))
      ENDIF
      DY(x1)=Ax1*(1.d0-Y(x1))-Bx1*Y(x1)

C SNHs 1/2/00 The fiddle factors below are not correct.
!      DY(m)=DY(m)*1.d-1
!      DY(h)=DY(h)*1.d-1
!      DY(j)=DY(j)*1.d-1
!      DY(d)=DY(d)*1.d-1
!      DY(f1)=DY(f1)*1.d-1
!      DY(x1)=DY(x1)*1.d-1
!      DY(Cai)=DY(Cai)*1.d-1
C SNH e
C SNHs 10/1/00 Adding electroporation, essentially a time dependent
C increase in conductivity, using the variable G.

      IF(KTYP33.EQ.4.AND.KTYP39.EQ.1) THEN !defib and electroporation
C        alpha = 2.5d-3
C MLB incorrect units!
        alpha = 2.5d-5
        beta  = 2.5d-5
        gamma = 1.0d-9
        DY(G)=alpha * DEXP(beta*(Y(Vm)-PARAM(Vmrest))*
     '    (Y(Vm)-PARAM(Vmrest)))*(1-DEXP(-gamma*(Y(Vm)-
     '    PARAM(Vmrest))*(Y(Vm)-PARAM(Vmrest))))
        DY(Vm)=-(DY(Vm)/PARAM(Cm))+((PROTOCOL(PseudoIs)+
     '    PROTOCOL(Is1current))/(PARAM(Cm)*PARAM(Am)))
     '    -Y(G)*Y(Vm)/(PARAM(Cm)*PARAM(Am))
      ELSE
        DY(G)=0.0d0
        DY(Vm)=-(DY(Vm)/PARAM(Cm))+((PROTOCOL(PseudoIs)+
     '    PROTOCOL(Is1current))/(PARAM(Cm)*PARAM(Am)))
      ENDIF
C SNHe

      RETURN
      END


