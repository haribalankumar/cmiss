      SUBROUTINE FPBSPL(T,N,K,X,L,H)
C#### Subroutine: FPBSPL
C###  Description:
C###    Auxililary routine for SPLDER.
c  subroutine fpbspl evaluates the (k+1) non-zero b-splines of
c  degree k at t(l) <= x < t(l+1) using the stable recurrence
c  relation of de boor and cox.
c  ..
      IMPLICIT NONE
c  ..scalar arguments..
      REAL*8 X
      INTEGER N,K,L
c  ..array arguments..
      REAL*8 T(N),H(6)
c  ..local scalars..
      REAL*8 F,ONE
      INTEGER I,J,LI,LJ
c  ..local arrays..
      REAL*8 HH(5)
c  ..
      ONE = 0.1D+01
      H(1) = ONE
      DO 20 J=1,K
        DO 10 I=1,J
          HH(I) = H(I)
  10    CONTINUE
        H(1) = 0.0d0
        DO 20 I=1,J
          LI = L+I
          LJ = LI-J
          F = HH(I)/(T(LI)-T(LJ))
          H(I) = H(I)+F*(T(LI)-X)
          H(I+1) = F*(X-T(LJ))
  20  CONTINUE
      RETURN
      END

