      SUBROUTINE ZZ(Z1,Z2,TRANS)

C#### Subroutine: ZZ
C###  Description:
C###    ZZ transforms rectangular cartesian coordinates Z1 to Z2
C###    using the transformation matrix TRANS.

      IMPLICIT NONE
!     Parameter List
C Note: Leave these arrays dimensioned dynamically. PJH 22Feb96
      REAL*8 TRANS(3,4),Z1(*),Z2(*)
!     Local Variables
      INTEGER i1,i2
      REAL*8 Z(3)

C     CALL ENTERS('ZZ',*9999)
      Z(1)=Z1(1)
      Z(2)=Z1(2)
      Z(3)=Z1(3)
      DO 2 i2=1,3
        Z2(i2)=0.0D0
        DO 1 i1=1,3
          Z2(i2)=Z2(i2)+Z(i1)*TRANS(i2,i1)
 1      CONTINUE
        Z2(i2)=Z2(i2)+TRANS(i2,4)
 2    CONTINUE

C     CALL EXITS('ZZ')
      RETURN
      END
