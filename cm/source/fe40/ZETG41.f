      SUBROUTINE ZETG41(CE,AXIAL_STRAIN,SZ,TENSION,XE,ZE,ERROR,*)

C#### Subroutine: ZETG41
C###  Description:
C###    ZETG41 is for truss elements only. Tension is in kN.

C**** Basis function is assumed to be 1D linear.
C**** CE(1) is Young's modulus in GPa
C**** CE(2) is truss Xsect.area in cm^2
C**** XE(ns,nj) has coordinates  in m
C**** ZE(ns,nj) has displacement in cm

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'

!     Parameter List
      REAL*8 CE(NMM),AXIAL_STRAIN,SZ,TENSION,XE(NSM,NJM),ZE(NSM,NHM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nj
      REAL*8 SX

      CALL ENTERS('ZETG41',*9999)
      SX=0.d0
      SZ=0.d0
      !add displacement in m onto initial coords
      DO nj=1,NJ_LOC(NJL_GEOM,0,1) !RGB temporary nr = 1
        SX=SX+ (XE(2,nj)-XE(1,nj))**2
        SZ=SZ+((XE(2,nj)+ZE(2,nj)*1.0D-2)-(XE(1,nj)+ZE(1,nj)*1.0D-2))**2
      ENDDO
      SX=DSQRT(SX)
      SZ=DSQRT(SZ)
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,*) ' SX=',SX,' SZ=',SZ
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF
      AXIAL_STRAIN=SZ/SX-1.d0
      TENSION=CE(1)*CE(2)*1.d2*AXIAL_STRAIN
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,*) ' Tension=',TENSION
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('ZETG41')
      RETURN
 9999 CALL ERRORS('ZETG41',ERROR)
      CALL EXITS('ZETG41')
      RETURN 1
      END


