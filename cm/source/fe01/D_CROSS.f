      SUBROUTINE D_CROSS(NITB,A,B,C,D_A,D_B,D_C)

C#### Subroutine: D_CROSS
C###  Description:
C###    D_CROSS returns the vector cross product of A and B in C
C###    and the NITB derivatives D_C of C given the derivatives
C###    D_A and D_B of A and B.

      IMPLICIT NONE
!     Parameter List
      INTEGER NITB
      REAL*8 A(3),B(3),C(3),D_A(3,3),D_B(3,3),D_C(3,3)
!     Local Variables
      INTEGER ni
C     CALL ENTERS('CROSS',*9999)
      C(1)=A(2)*B(3)-A(3)*B(2)
      C(2)=A(3)*B(1)-A(1)*B(3)
      C(3)=A(1)*B(2)-A(2)*B(1)
      DO ni=1,NITB
        D_C(1,ni)=D_A(2,ni)*B(3)-D_A(3,ni)*B(2)
     '    +A(2)*D_B(3,ni)-A(3)*D_B(2,ni)
        D_C(2,ni)=D_A(3,ni)*B(1)-D_A(1,ni)*B(3)
     '    +A(3)*D_B(1,ni)-A(1)*D_B(3,ni)
        D_C(3,ni)=D_A(1,ni)*B(2)-D_A(2,ni)*B(1)
     '    +A(1)*D_B(2,ni)-A(2)*D_B(1,ni)
      ENDDO
C     CALL EXITS('CROSS')
      RETURN
      END


