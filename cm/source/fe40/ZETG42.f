      SUBROUTINE ZETG42(NBH,ng,CG,DS0,DS1,DSDXI,EIV1,EIV2,PG,PV,RL,XE,
     '  ZE,ERROR,*)

C#### Subroutine: ZETG42
C###  Description:
C###    ZETG42 is for batten elements only. Tension is in kN. Element
C###    basis must be cubic Hermite.

C**** CG(1)  is Young's modulus in GPa
C**** CG(2)  is Xsection area in cm^2
C**** CG(3)  is 2nd moment of area Izz in cm^4
C**** CG(4)  is axial force/unit length in kN/m
C**** CG(5)  is transverse force/unit length in kN/m
C**** DXDXI0 is deriv of undeformed X wrt Xi
C**** DXDXI1 is deriv of deformed   X wrt Xi
C**** DYDXI0 is deriv of undeformed Y wrt Xi
C**** DYDXI1 is deriv of deformed   Y wrt Xi
C**** DSDXI0 is rate of change of arclength with Xi in undeformed element
C**** DSDXI1 is rate of change of arclength with Xi in deformed   element
C**** DSDXI  is difference in rate of change of arc length with Xi
C**** EI     is Young's modulus * 2nd moment of area in kN.m^2
C**** DYDS0  is dY/dS undeformed, D2YDS20 is d^2Y/dS^2 undeformed
C**** DYDS1  is dY/dS deformed  , D2YDS21 is d^2Y/dS^2 deformed
C**** CURV0  is undeformed curvature
C**** CURV1  is deformed   curvature
C**** EIV1   is EIV*SQRT(1-dY/dS^2)
C**** EIV2   is EIV*d^2Y/dS^2*dY/dS/(1-dY/dS^2)**(3/2)
C**** PV     is CG(4)*dY/dS*SQRT(1-d^2Y/dS^2)
C**** RL     is change in arc length dS

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER NBH(NHM),ng
      REAL*8 CG(NMM),DS0,DS1,DSDXI,EIV1,EIV2,PG(NSM,NUM,NGM,NBM),PV,RL,
     '  XE(NSM,NJM),ZE(NSM,NHM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,ns,nx
      REAL*8 D2YDS21,DSDXI0,DSDXI1,DXDXI0,DXDXI1,DYDS1,DYDXI0,DYDXI1,EI

      CALL ENTERS('ZETG42',*9999)
      nx=1 !temporary
      nb=NBH(NH_LOC(1,nx))
      DXDXI0=0.0d0
      DXDXI1=0.0d0
      DO ns=1,NST(nb)
        DXDXI0=DXDXI0+PG(ns,2,ng,nb)*XE(ns,1)
        DXDXI1=DXDXI1+PG(ns,2,ng,nb)*ZE(ns,1)
      ENDDO
      nb=NBH(NH_LOC(2,nx))
      DYDXI0=0.0d0
      DYDXI1=0.0d0
      DO ns=1,NST(nb)
        DYDXI0=DYDXI0+PG(ns,2,ng,nb)*XE(ns,2)
        DYDXI1=DYDXI1+PG(ns,2,ng,nb)*ZE(ns,2)
      ENDDO
      DSDXI0=DSQRT(DXDXI0**2+DYDXI0**2)
      DSDXI1=DSQRT(DXDXI1**2+DYDXI1**2)
      DSDXI =DSDXI1-DSDXI0
       IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$       call mp_setlock()
         WRITE(OP_STRING,'('' dS/dXi_0='',E11.3)')DSDXI0
         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
         WRITE(OP_STRING,'('' dS/dXi_1='',E11.3)')DSDXI1
         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
         WRITE(OP_STRING,'('' dS/dXi='',E11.3)')DSDXI
         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      DS0=1.0d0/DSDXI0
      DS1=1.0d0/DSDXI1

      nb=NBH(NH_LOC(2,nx))
C      DYDS0=0.0d0
      DYDS1=0.0d0
      DO 100 ns=1,NST(nb)
C        DYDS0=DYDS0+PG(ns,2,ng,nb)*DS0*XE(ns,2)
        DYDS1=DYDS1+PG(ns,2,ng,nb)*DS1*ZE(ns,2)
 100  CONTINUE
      IF(DABS(DYDS1-1.d0).LT.1.0D-6) DYDS1=0.99999d0
C      IF(DOP) WRITE(IO4,'('' dY/dS undeformed='',E11.3)') DYDS0
      IF(DOP) THEN
        WRITE(OP_STRING,'('' dY/dS deformed  ='',E11.3)') DYDS1
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      EI=CG(1)*CG(3)*1.0D-2
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'(''  EI='',E11.3)') EI
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF
      nb=NBH(NH_LOC(2,nx))
C      D2YDS20=0.0d0
      D2YDS21=0.0d0
      DO 200 ns=1,NST(nb)
C        D2YDS20=D2YDS20+PG(ns,3,ng,nb)*DS0*DS0*XE(ns,2)
        D2YDS21=D2YDS21+PG(ns,3,ng,nb)*DS1*DS1*ZE(ns,2)
 200  CONTINUE
C      IF(DOP) WRITE(IO4,'('' d^2Y/dS^2 undeformed='',E11.3)') D2YDS20
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' d^2Y/dS^2 deformed  ='',E11.3)') D2YDS21
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

c Fried
      EIV1=EI*D2YDS21/(1.d0-DYDS1**2)
      EIV2=EI*DYDS1*D2YDS21**2/(1.d0-DYDS1**2)**2
      PV=CG(4)*DYDS1/DSQRT(DABS(1.d0-DYDS1**2))
c
c length constraint
       RL=2.d0*DSDXI*DXDXI1/DSDXI1

C
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' EIV1='',E11.3)') EIV1
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' EIV2='',E11.3)') EIV2
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' PV='',E11.3)') PV
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' RL='',E11.3)') RL
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('ZETG42')
      RETURN
 9999 CALL ERRORS('ZETG42',ERROR)
      CALL EXITS('ZETG42')
      RETURN 1
      END


