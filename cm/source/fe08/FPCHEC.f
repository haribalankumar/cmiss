      SUBROUTINE FPCHEC(X,M,T,N,K,IER)
C#### Subroutine: FPCHEC
C###  Description:
C###    Auxililary routine for CURFIT.
c  subroutine fpchec verifies the number and the position of the knots
c  t(j),j=1,2,...,n of a spline of degree k, in relation to the number
c  and the position of the data points x(i),i=1,2,...,m. if all of the
c  following conditions are fulfilled, the error parameter ier is set
c  to zero. if one of the conditions is violated ier is set to ten.
c      1) k+1 <= n-k-1 <= m
c      2) t(1) <= t(2) <= ... <= t(k+1)
c         t(n-k) <= t(n-k+1) <= ... <= t(n)
c      3) t(k+1) < t(k+2) < ... < t(n-k)
c      4) t(k+1) <= x(i) <= t(n-k)
c      5) the conditions specified by schoenberg and whitney must hold
c         for at least one subset of data points, i.e. there must be a
c         subset of data points y(j) such that
c             t(j) < y(j) < t(j+k+1), j=1,2,...,n-k-1
c  ..
      IMPLICIT NONE
c  ..scalar arguments..
      INTEGER M,N,K,IER
c  ..array arguments..
      REAL*8 X(M),T(N)
c  ..local scalars..
      INTEGER I,J,K1,K2,L,NK1,NK2,NK3
      REAL*8 TJ,TL
c  ..
      K1 = K+1
      K2 = K1+1
      NK1 = N-K1
      NK2 = NK1+1
      IER = 10
c  check condition no 1
      IF(NK1.LT.K1 .OR. NK1.GT.M) GO TO 80
c  check condition no 2
      J = N
      DO 20 I=1,K
        IF(T(I).GT.T(I+1)) GO TO 80
        IF(T(J).LT.T(J-1)) GO TO 80
        J = J-1
  20  CONTINUE
c  check condition no 3
      DO 30 I=K2,NK2
        IF(T(I).LE.T(I-1)) GO TO 80
  30  CONTINUE
c  check condition no 4
      IF(X(1).LT.T(K1) .OR. X(M).GT.T(NK2)) GO TO 80
c  check condition no 5
      IF(X(1).GE.T(K2) .OR. X(M).LE.T(NK1)) GO TO 80
      I = 1
      L = K2
      NK3 = NK1-1
      IF(NK3.LT.2) GO TO 70
      DO 60 J=2,NK3
        TJ = T(J)
        L = L+1
        TL = T(L)
  40    I = I+1
        IF(I.GE.M) GO TO 80
        IF(X(I).LE.TJ) GO TO 40
        IF(X(I).GE.TL) GO TO 80
  60  CONTINUE
  70  IER = 0
  80  RETURN
      END

