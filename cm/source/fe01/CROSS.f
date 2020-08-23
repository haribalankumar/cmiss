      SUBROUTINE CROSS(A,B,C)

C#### Subroutine: CROSS
C###  Description:
C###    CROSS returns the vector cross product of A*B in C.

      IMPLICIT NONE
!     Parameter List
      REAL*8 A(3),B(3),C(3)

C     CALL ENTERS('CROSS',*9999)
      C(1)=A(2)*B(3)-A(3)*B(2)
      C(2)=A(3)*B(1)-A(1)*B(3)
      C(3)=A(1)*B(2)-A(2)*B(1)

C     CALL EXITS('CROSS')
      RETURN
      END


