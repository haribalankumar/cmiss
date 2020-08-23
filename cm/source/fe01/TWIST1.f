      SUBROUTINE TWIST1(AMOUNT,TRANS)

C#### Subroutine: TWIST1
C###  Description:
C###    TWIST1 updates the transformation matrix TRANS with a rotation
C###    of AMOUNT radians about the Z(1) axis.

      IMPLICIT NONE
!     Parameter List
      REAL*8 AMOUNT,TRANS(3,4)
!     Local Variables
      INTEGER n
      REAL*8 COSA,SINA,TRANS2,TRANS3

C     CALL ENTERS('TWIST1',*9999)
      IF(amount.NE.0.0D0) THEN
        SINA=DSIN(AMOUNT)
        COSA=DCOS(AMOUNT)
        DO n=1,4
          TRANS2=TRANS(2,N)
          TRANS3=TRANS(3,N)
          TRANS(2,N)=TRANS2*COSA-TRANS3*SINA
          TRANS(3,N)=TRANS2*SINA+TRANS3*COSA
        ENDDO
      ENDIF

C     CALL EXITS('TWIST1')
      RETURN
      END


