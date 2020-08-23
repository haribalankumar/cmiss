      SUBROUTINE SCALE2(AMOUNT,TRANS)

C#### Subroutine: SCALE2
C###  Description:
C###    SCALE2 updates the transformation matrix TRANS with a scaling
C###    by AMOUNT along the Z(2) axis.

      IMPLICIT NONE
!     Parameter List
      REAL*8 AMOUNT,TRANS(3,4)
!     Local Variables
      INTEGER n

C     CALL ENTERS('SCALE2',*9999)
      IF(amount.NE.1.0D0) THEN
        DO n=1,4
          TRANS(2,n)= TRANS(2,n)*AMOUNT
        ENDDO
      ENDIF

C     CALL EXITS('SCALE2')
      RETURN
      END


