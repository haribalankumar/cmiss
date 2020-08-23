      REAL*8 FUNCTION REFINE_GETXIPOS(i1,IBT,IDRN,nb,XII)

C#### Function: REFINE_GETXIPOS
C###  Type: REAL*8
C###  Description:
C###    REFINE_GETXIPOS returns the Xi location (in the IDRN) direction
C###    that corresponds to step i1 of the refine process.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER i1,IBT(3,NIM,NBFM),IDRN,nb
      REAL*8 XII
!     Local Variables
      REAL*8 XI

      IF((IBT(1,IDRN,nb).EQ.1.OR.IBT(1,IDRN,nb).EQ.5.OR.
     '  IBT(1,IDRN,nb).EQ.6).AND.IBT(2,IDRN,nb).EQ.2)
     '  THEN !quadratic Lagrange
        IF(XII.EQ.0.5d0) THEN
          IF(i1.EQ.1) THEN
            XI=0.25d0
          ELSE
            XI=0.75d0
          ENDIF
        ELSE IF(XII.GT.0.5d0) THEN
          IF(i1.EQ.1) THEN
            XI=XII
          ELSE
            XI=(XII+1.0d0)/2.0d0
          ENDIF
        ELSE
          IF(i1.EQ.1) THEN
            XI=XII/2.0d0
          ELSE
            XI=XII
          ENDIF
        ENDIF
      ELSE IF((IBT(1,IDRN,nb).EQ.1.OR.
     '    IBT(1,IDRN,nb).EQ.5.OR.IBT(1,IDRN,nb).EQ.6).AND.
     '    IBT(2,IDRN,nb).EQ.3) THEN !cubic Lagrange
        IF(XII.EQ.(1.0d0/3.0d0)) THEN
          IF(i1.EQ.1) THEN
            XI=1.0d0/9.0d0
          ELSE IF(i1.EQ.2) THEN
            XI=2.0d0/9.0d0
          ELSE
            XI=0.5d0
          ENDIF
        ELSE IF(XII.EQ.(2.0d0/3.0d0)) THEN
          IF(i1.EQ.1) THEN
            XI=0.5d0
          ELSE IF(i1.EQ.2) THEN
            XI=7.0d0/9.0d0
          ELSE
            XI=8.0d0/9.0d0
          ENDIF
        ELSE IF(XII.GT.(2.0d0/3.0d0)) THEN
          IF(i1.EQ.1) THEN
            XI=XII
          ELSE IF(i1.EQ.2) THEN
            XI=(2.0d0*XII+1.0d0)/3.0d0
          ELSE
            XI=(XII+2.0d0)/3.0d0
          ENDIF
        ELSE IF(XII.GT.(1.0d0/3.0d0)) THEN
          IF(i1.EQ.1) THEN
            XI=(3.0d0*XII+1.0d0)/6.0d0
          ELSE IF(i1.EQ.2) THEN
            XI=XII
          ELSE
            XI=(3.0d0*XII+2.0d0)/6.0d0
          ENDIF
        ELSE
          IF(i1.EQ.1) THEN
            XI=XII/3.0d0
          ELSE IF(i1.EQ.2) THEN
            XI=2.0d0*XII/3.0d0
          ELSE
            XI=XII
          ENDIF
        ENDIF
      ELSE
        XI=XII
      ENDIF

      REFINE_GETXIPOS=XI

      RETURN
      END


C FE05 Functions
C ==============

