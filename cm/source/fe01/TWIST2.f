      SUBROUTINE TWIST2(AMOUNT,TRANS)

C#### Subroutine: TWIST2
C###  Description:
C###    TWIST2 updates the transformation matrix TRANS with a rotation
C###    of AMOUNT radians about the Z(2) axis.

      IMPLICIT NONE
!     Parameter List
      REAL*8 AMOUNT,TRANS(3,4)
!     Local Variables
      INTEGER n
      REAL*8 COSA,SINA,TRANS1,TRANS3

C     CALL ENTERS('TWIST2',*9999)
      IF(amount.NE.0.0D0) THEN
        SINA=DSIN(AMOUNT)
        COSA=DCOS(AMOUNT)
        DO n=1,4
          TRANS3=TRANS(3,n)
          TRANS1=TRANS(1,n)
          TRANS(3,n)=TRANS3*COSA-TRANS1*SINA
          TRANS(1,n)=TRANS3*SINA+TRANS1*COSA
        ENDDO
      ENDIF

C     CALL EXITS('TWIST2')
      RETURN
      END


