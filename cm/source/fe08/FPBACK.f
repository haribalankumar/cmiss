      SUBROUTINE FPBACK(A,Z,N,K,C,NEST)
C#### Subroutine: FPBACK
C###  Description:
C###    Auxililary routine for CURFIT.
c  subroutine fpback calculates the solution of the system of
c  equations a*c = z with a a n x n upper triangular matrix
c  of bandwidth k.
c  ..
      IMPLICIT NONE
c  ..scalar arguments..
      INTEGER N,K,NEST
c  ..array arguments..
      REAL*8 A(NEST,K),Z(N),C(N)
c  ..local scalars..
      REAL*8 STORE
      INTEGER I,I1,J,K1,L,M
c  ..
      K1 = K-1
      C(N) = Z(N)/A(N,1)
      I = N-1
      IF(I.EQ.0) GO TO 30
      DO 20 J=2,N
        STORE = Z(I)
        I1 = K1
        IF(J.LE.K1) I1 = J-1
        M = I
        DO 10 L=1,I1
          M = M+1
          STORE = STORE-C(M)*A(I,L+1)
  10    CONTINUE
        C(I) = STORE/A(I,1)
        I = I-1
  20  CONTINUE
  30  RETURN
      END

