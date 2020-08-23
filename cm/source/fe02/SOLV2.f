      SUBROUTINE SOLV2(A11,A12,A21,A22,B1,B2,X1,X2,SINGLR,ERROR,*)

C#### Subroutine: SOLV2
C###  Description:
C###    SOLV2 solves a 2*2 system of linear equations.

      IMPLICIT NONE
!     Parameter List
      REAL*8 A11,A12,A21,A22,B1,B2,X1,X2
      CHARACTER ERROR*(*)
      LOGICAL SINGLR
!     Local Variables
      REAL*8 DET,TOL
      DATA TOL/1.D-07/

      CALL ENTERS('SOLV2',*9999)
      X1=0.0D0
      X2=0.0D0
      DET=A11*A22-A12*A21
      IF(DABS(DET).LT.TOL) THEN
        SINGLR=.TRUE.
        GO TO 9998
      ELSE
        SINGLR=.FALSE.
      ENDIF
      X1=(A22*B1-A12*B2)/DET
      X2=(A11*B2-A21*B1)/DET

 9998 CALL EXITS('SOLV2')
      RETURN
 9999 CALL ERRORS('SOLV2',ERROR)
      CALL EXITS('SOLV2')
      RETURN 1
      END

