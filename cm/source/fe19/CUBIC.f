      SUBROUTINE CUBIC(TIME,Y,DY,CONTROL,MODEL,SIZES,VARIANT,DERIVED,
     '  PARAM,PROTOCOL,AII,AIO,ARI,ARO,ERR)

C#### Subroutine: CUBIC
C###  Description:
C###    CUBIC returns the right hand side vector for the cubic
C###    cellular model as given by Hunter,McNaughton & Noble,1975.
C**** (re)Written by Martin Buist, 6 May 1999

C***  SIZES(1) is 1 equation (Vm).
C***  SIZES(2) is not used.
C***  SIZES(3) is not used.
C***  SIZES(4) is not used.
C***  SIZES(5) is 6 parameters.
C***  PARAM(1) is the membrane capacitance Cm (corrected by 10^-6).
C***  PARAM(2) is the surface to volume ration Am.
C***  PARAM(3) is resting potential in mV.
C***  PARAM(4) is plateau potential in mV.
C***  PARAM(5) is threshold potential in mV.
C***  PARAM(6) is membrane conductance g (corrected by 10^-6).
C***  SIZES(6) is 4 protocols.
C***  PROTOCOL(1) is pseudo stimulus current.
C***  PROTOCOL(4) is applied stimulus current.
C***  SIZES(7) is not used.
C***  SIZES(8) is not used.
C***  SIZES(9) is not used.
C***  SIZES(10) is not used.
C***  SIZES(11) is not used.
C***  Result is in units of mA/mm^2.

      IMPLICIT NONE

      INCLUDE 'cell_cubic.inc'
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
      REAL*8 PHI,PLATEAU,THRESHOLD,VTH,VPL

      PLATEAU=PARAM(Vmplat)-PARAM(Vmrest)
      THRESHOLD=PARAM(Vmthres)-PARAM(Vmrest)

      PHI=Y(Vm)-PARAM(Vmrest)
      VTH=PHI/THRESHOLD
      VPL=PHI/PLATEAU
      IF(KTYP33.EQ.1) THEN      !Cubic
        DY(Vm)=PARAM(g)*PHI*(1.d0-VTH)*(1.d0-VPL)
      ELSE IF(KTYP33.EQ.2) THEN !Quintic
        DY(Vm)=PARAM(g)*PHI*(1.d0-VTH*VTH)*(1.d0-VPL*VPL)
      ELSE IF(KTYP33.EQ.3) THEN !Seventh-order
        DY(Vm)=PARAM(g)*PHI*(1.d0-VTH*VTH*VTH)*(1.d0-VPL*VPL*VPL)
      ENDIF
      DY(Vm)=-(DY(Vm)/PARAM(Cm))+((PROTOCOL(PseudoIs)+
     '  PROTOCOL(Is1current))/(PARAM(Cm)*PARAM(Am)))

      RETURN
      END


