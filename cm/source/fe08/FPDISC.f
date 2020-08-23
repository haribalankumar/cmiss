      SUBROUTINE FPDISC(T,N,K2,B,NEST)
C#### Subroutine: FPDISC
C###  Description:
C###    Auxililary routine for CURFIT.
c  subroutine fpdisc calculates the discontinuity jumps of the kth
c  derivative of the b-splines of degree k at the knots t(k+2)..t(n-k-1)
c  ..
      IMPLICIT NONE
c  ..scalar arguments..
      INTEGER N,K2,NEST
c  ..array arguments..
      REAL*8 T(N),B(NEST,K2)
c  ..local scalars..
      REAL*8 AN,FAC,PROD
      INTEGER I,IK,J,JK,K,K1,L,LJ,LK,LMK,LP,NK1,NRINT
c  ..local array..
      REAL*8 H(12)
c  ..
      K1 = K2-1
      K = K1-1
      NK1 = N-K1
      NRINT = NK1-K
      AN = NRINT
      FAC = AN/(T(NK1+1)-T(K1))
      DO 40 L=K2,NK1
        LMK = L-K1
        DO 10 J=1,K1
          IK = J+K1
          LJ = L+J
          LK = LJ-K2
          H(J) = T(L)-T(LK)
          H(IK) = T(L)-T(LJ)
  10    CONTINUE
        LP = LMK
        DO 30 J=1,K2
          JK = J
          PROD = H(J)
          DO 20 I=1,K
            JK = JK+1
            PROD = PROD*H(JK)*FAC
  20      CONTINUE
          LK = LP+K1
          B(LMK,J) = (T(LK)-T(LP))/PROD
          LP = LP+1
  30    CONTINUE
  40  CONTINUE
      RETURN
      END

