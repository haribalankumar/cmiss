      SUBROUTINE SCALE3(AMOUNT,TRANS)

C#### Subroutine: SCALE3
C###  Description:
C###    SCALE3 updates the transformation matrix TRANS with a scaling
C###    by AMOUNT along the Z(3) axis.

      IMPLICIT NONE
!     Parameter List
      REAL*8 AMOUNT,TRANS(3,4)
!     Local Variables
      INTEGER n

C     CALL ENTERS('SCALE3',*9999)
      IF(amount.NE.1.0D0) THEN
        DO n=1,4
          TRANS(3,n)= TRANS(3,n)*AMOUNT
        ENDDO
      ENDIF

C     CALL EXITS('SCALE3')
      RETURN
      END


