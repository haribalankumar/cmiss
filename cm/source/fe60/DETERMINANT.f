      REAL*8 FUNCTION DETERMINANT(I,J,K,A,LDA)

C#### Function: DETERMINANT
C###  Description: Computes the determinant of a 3x3 matrix A with row
C###    indices I, J and K.
CC JMB 15-NOV-2001

      IMPLICIT NONE
!     Parameter List
      INTEGER I,J,K,LDA
      REAL*8 A(LDA,*)

      DETERMINANT=A(I,1)*A(J,2)*A(K,3)-A(I,1)*A(K,2)*A(J,3)+A(J,1)
     '  *A(K,2)*A(I,3)-A(J,1)*A(I,2)*A(K,3)+A(K,1)*A(I,2)*A(J,3)
     '  -A(K,1)*A(J,2)*A(I,3)

      RETURN
      END


