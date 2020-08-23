      SUBROUTINE VCD(TIME,Y,DY,CONTROL,MODEL,SIZES,VARIANT,DERIVED,
     '  PARAM,PROTOCOL,AII,AIO,ARI,ARO,ERR)

C#### Subroutine: VCD
C###  Description:
C###    VCD returns the right hand side vector for the
C###    VanCapelle-Durrer electrical cellular model, with
C###    modifications from UCLA provided by Alan
C###    Garfinkel. Original reference in fax dated March 8, 1994.
C**** (re)Written by Martin Buist, 7 May 1999

C***  SIZES(1) is 5 equations (Vm,Recov,Cai,DVm,DRecov).
C***  SIZES(2) is not used.
C***  SIZES(3) is not used.
C***  SIZES(4) is not used.
C***  SIZES(5) is 6 parameters.
C***  PARAM(1) is the membrane capacitance Cm (corrected by 10^-6).
C***  PARAM(2) is the surface to volume ration Am.
C***  PARAM(3) is resting potential in mV.
C***  PARAM(4) is the time constant.
C***  PARAM(5) is the Ca2+ time constant.
C***  PARAM(6) is the repolarisation time constant.
C***  SIZES(6) is 4 protocols.
C***  PROTOCOL(1) is pseudo stimulus current.
C***  PROTOCOL(4) is applied stimulus current.
C***  SIZES(7) is not used.
C***  SIZES(8) is not used.
C***  SIZES(9) is not used.
C***  SIZES(10) is not used.
C***  SIZES(11) is not used.

      IMPLICIT NONE

      INCLUDE 'cell_vcd.inc'
      INCLUDE 'cell_reserved.inc'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'tol00.cmn'

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
      REAL*8 F1,I0,I1,NORM_PHI,PHI,T_CONST

CMB Units correction was done incorrectly, original in uA/cm^2
CMB changing to uA/mm^2
      IF(Y(Vm).LT.-70.0d0) THEN
CMB        I1=5.0d0+0.5d0*(Y(Vm)+70.0d0)
        I1=0.05d0+0.005d0*(Y(Vm)+70.0d0)
      ELSE IF(Y(Vm).GT.0.0d0) THEN
CMB        I1=6.0d0+0.425d0*Y(Vm)
        I1=0.06d0+0.00425d0*Y(Vm)
      ELSE
CMB        I1=5.0d0+((Y(Vm)+70.0d0)/70.0d0)
        I1=0.05d0+0.01d0*((Y(Vm)+70.0d0)/70.0d0)
      ENDIF

      IF(Y(Vm).LT.-74.3d0) THEN
CMB        F1=7.84d0+2.0d0*(Y(Vm)+74.3d0)
        F1=0.0784d0+0.02d0*(Y(Vm)+74.3d0)
      ELSE IF(Y(Vm).GT.-27.8d0) THEN
CMB        F1=-98.84d0+1.71d0*(Y(Vm)+27.8d0)
        F1=-0.9884d0+0.0171d0*(Y(Vm)+27.8d0)
      ELSE
CMB        F1=((3.837854d-3*Y(Vm)+0.584649d0)*Y(Vm)+25.31834d0)*Y(Vm)
CMB     '    +235.6256d0
        F1=((3.837854d-5*Y(Vm)+0.00584649d0)*Y(Vm)+0.2531834d0)*Y(Vm)
     '    +2.356256d0
      ENDIF

      IF(KTYP33.EQ.2) THEN  !Calif. mods to VCD
        F1=4.0d0*F1
        IF(KTYP34.EQ.2) THEN    !"Ischemic" APD (112ms)
          IF(Y(DVm).LT.0) I1=2.0d0*I1
        ENDIF
      ENDIF

      I0=I1+F1
CMB      DY(Vm)=(Y(Recov)*I1+(1.0d0-Y(Recov))*I0)*1.0d-3
      DY(Vm)=(Y(Recov)*I1+(1.0d0-Y(Recov))*I0)

! recovery
C      PHI=Y(Vm)+DY(Vm)*DT
      PHI=Y(Vm)
      IF(KTYP33.EQ.1) THEN !Original VCD
        T_CONST=PARAM(TConst)
      ELSE IF(KTYP33.EQ.2) THEN !VCDC mods
        ! From modifications provided by Alan Garfinkel
        IF(Y(DRecov).GE.-ZERO_TOL) THEN
          !recovery var Y is increasing
          IF(KTYP34.EQ.1) THEN !"Normal" APD (209ms)
            T_CONST=0.5d0
          ELSE IF(KTYP34.EQ.2) THEN !"Ischemic" APD (112ms)
            T_CONST=0.33d0
          ENDIF
        ELSE IF(Y(Recov).GT.0.85d0) THEN !Y decreasing & Y>0.85
          IF(KTYP34.EQ.1) THEN !"Normal" APD (209ms)
            T_CONST=0.1d0
          ELSE IF(KTYP34.EQ.2) THEN !"Ischemic" APD (112ms)
            T_CONST=0.066d0
          ENDIF
        ELSE                          !Y decreasing & Y<0.85
          IF(KTYP34.EQ.1) THEN !"Normal" APD (209ms)
            T_CONST=3.0d0
          ELSE IF(KTYP34.EQ.2) THEN !"Ischemic" APD (112ms)
            !PJH 5/7/98   T_CONST=3.31d0
            T_CONST=PARAM(RepTConst)
          ENDIF
        ENDIF
        ! Scale Factor for time constant
        T_CONST=T_CONST*PARAM(TConst)
      ENDIF

      IF(PHI.LT.-80.0d0) THEN
        DY(Recov)=-Y(Recov)/T_CONST
      ELSE IF(PHI.GT.-60.0d0) THEN
        DY(Recov)=(1.0d0-Y(Recov))/T_CONST
      ELSE
        DY(Recov)=((PHI+80.0d0)/20.0d0-Y(Recov))/T_CONST
      ENDIF
      IF(KTYP33.EQ.2) THEN !VCDC mods
        Y(DRecov)=DY(Recov)
      ENDIF

! calcium
      IF(DABS(PARAM(CaTConst)).GT.ZERO_TOL) THEN
        NORM_PHI=(PHI-PARAM(Vmrest))/100.0d0
        IF(NORM_PHI.LT.ZERO_TOL) NORM_PHI=0.0d0
        DY(Cai)=(NORM_PHI-Y(Cai))/PARAM(CaTConst)
      ELSE
        Y(Cai)=0.0d0
        DY(Cai)=0.0d0
      ENDIF

      DY(DVm)=0.0d0
      DY(DRecov)=0.0d0

      DY(Vm)=-(DY(Vm)/PARAM(Cm))+((PROTOCOL(PseudoIs)+
     '  PROTOCOL(Is1current))/(PARAM(Cm)*PARAM(Am)))

      RETURN
      END


