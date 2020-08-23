      REAL*8 FUNCTION DACOS_MOD(X)

C#### Function: DACOS_MOD
C###  Type: REAL*8
C###  Description:
C###    DACOS_MOD returns a DACOS of X and if X is within TOLERANCE
C###    of 1.0d0 or -1.0d0 uses returns zero (ie DACOS(1.0d0))

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      REAL*8 X
!     Local Variables
      REAL*8 DACOS,ONE,TOLERANCE
      CHARACTER ERROR*10
      DATA ONE/1.0d0/, TOLERANCE/1.0d-10/

      IF(DABS(X).GT.(ONE+TOLERANCE)) THEN
        WRITE(OP_STRING,'('' >>ERROR: X is invalid in DACOS_MOD'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        DACOS_MOD=0.0d0
      ELSE IF(DABS(DABS(X)-ONE).LT.TOLERANCE) THEN
        DACOS_MOD=0.0d0
      ELSE
        DACOS_MOD=DACOS(X)
      ENDIF

 9999 RETURN
      END


