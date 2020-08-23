      REAL*8 FUNCTION FPRATI(P1,F1,P2,F2,P3,F3)
C#### Subroutine: FPRATI
C###  Description:
C###    Auxililary routine for CURFIT.
c  given three points (p1,f1),(p2,f2) and (p3,f3), function fprati
c  gives the value of p such that the rational interpolating function
c  of the form r(p) = (u*p+v)/(p+w) equals zero at p.
c  ..
      IMPLICIT NONE
c  ..scalar arguments..
      REAL*8 P1,F1,P2,F2,P3,F3
c  ..local scalars..
      REAL*8 H1,H2,H3,P
c  ..
      IF(P3.GT.0.0d0) GO TO 10
c  value of p in case p3 = infinity.
      P = (P1*(F1-F3)*F2-P2*(F2-F3)*F1)/((F1-F2)*F3)
      GO TO 20
c  value of p in case p3 ^= infinity.
  10  H1 = F1*(F2-F3)
      H2 = F2*(F3-F1)
      H3 = F3*(F1-F2)
      P = -(P1*P2*H3+P2*P3*H1+P3*P1*H2)/(P1*H1+P2*H2+P3*H3)
c  adjust the value of p1,f1,p3 and f3 such that f1 > 0 and f3 < 0.
  20  IF(F2.LT.0.0d0) GO TO 30
      P1 = P2
      F1 = F2
      GO TO 40
  30  P3 = P2
      F3 = F2
  40  FPRATI = P
      RETURN
      END

