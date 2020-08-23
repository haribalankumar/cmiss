      REAL*8 FUNCTION DET(A)

C#### Function: DET
C###  Type: REAL*8
C###  Description:
C###    DET evaluates determinant of 3*3 matrix A.

      IMPLICIT NONE
!     Parameter List
      REAL*8 A(3,3)

      DET=A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
     '   +A(1,2)*(A(2,3)*A(3,1)-A(3,3)*A(2,1))
     '   +A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))

      RETURN
      END


