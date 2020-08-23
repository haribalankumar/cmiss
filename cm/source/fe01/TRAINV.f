      SUBROUTINE TRAINV(T1,T2)

C#### Subroutine: NAME
C###  Description:
C###    NAME inverts transformation matrix T1 into T2.

      IMPLICIT NONE
!     Parameter List
      REAL*8 T1(3,4),T2(3,4)
!     Local Variables
      INTEGER i1,i2
      REAL*8 A,Z(3),Z2(3)

C     CALL ENTERS('TRAINV',*9999)
      CALL INVERT(3,T1,T2,A)
      Z(1)=-T1(1,4)
      Z(2)=-T1(2,4)
      Z(3)=-T1(3,4)
      DO 2 i2=1,3
        Z2(i2)=0.0d0
        DO 1 i1=1,3
          Z2(i2)=Z2(i2)+Z(i1)*T2(i2,i1)
 1      CONTINUE
 2    CONTINUE
      DO 3 i2=1,3
        T2(i2,4)=Z2(i2)
 3    CONTINUE

C     CALL EXITS('TRAINV')
      RETURN
      END


