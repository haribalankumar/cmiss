      SUBROUTINE TWIST3(AMOUNT,TRANS)

C#### Subroutine: TWIST3
C###  Description:
C###    TWIST3 updates the transformation matrix TRANS with a rotation
C###    of AMOUNT radians about the Z(3) axis.

      IMPLICIT NONE
!     Parameter List
      REAL*8 AMOUNT,TRANS(3,4)
!     Local Variables
      INTEGER n
      REAL*8 COSA,SINA,TRANS1,TRANS2

C     CALL ENTERS('TWIST3',*9999)
      IF(amount.NE.0.0D0) THEN
        SINA=DSIN(AMOUNT)
        COSA=DCOS(AMOUNT)
        DO n=1,4
          TRANS1=TRANS(1,n)
          TRANS2=TRANS(2,n)
          TRANS(1,n)=TRANS1*COSA-TRANS2*SINA
          TRANS(2,n)=TRANS1*SINA+TRANS2*COSA
        ENDDO
      ENDIF

C     CALL EXITS('TWIST3')
      RETURN
      END


