      REAL*8 FUNCTION I0(X)

C#### Function: I0
C###  Type: REAL*8
C###  Description:
C###    I0 calculates the modified Bessel function of the first kind of
C###    order 0 using the approximation of Abromowitz and Stegun.

      IMPLICIT NONE
!     Parameter List
      REAL*8 X
!     Local Variables
      REAL*8 A1,A2,A3,A4,A5,A6,A7,T

      !Calculates I0(x) for x < 3.75
      T=X*X/(3.75D0*3.75D0)
      A1=1.0D0
      A2=3.5156229D0
      A3=3.0899424D0
      A4=1.2067492D0
      A5=0.2659732D0
      A6=0.0360768D0
      A7=0.0045813D0
      I0=A1+(A2+(A3+(A4+(A5+(A6+A7*T)*T)*T)*T)*T)*T
      RETURN
      END


