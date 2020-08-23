      REAL*8 FUNCTION K1(X)

C#### Function: K1
C###  Type: REAL*8
C###  Description:
C###    K1 calculates the modified Bessel function of the second kind of
C###    order 1 using the approximation of Abromowitz and Stegun.

      IMPLICIT NONE
!     Paramter List
      REAL*8 X
!     Local Variables
      REAL*8 I1
      REAL*8 A1,A2,A3,A4,A5,A6,A7,A8,T

      !Calculates K1(x)
      IF (X.LE.2.0D0) THEN
         T=(X/2.0D0)*(X/2.0D0)
         A1=DLOG(X/2.0D0)*I1(X)
         A2=1.0D0
         A3=0.15443144D0
         A4=-0.67278579D0
         A5=-0.18156897D0
         A6=-0.01919402D0
         A7=-0.00110404D0
         A8=-0.00004686D0
         K1=A1+A2/X+(A3+(A4+(A5+(A6+(A7+A8*T)*T)*T)*T)*T)*X/4
      ELSE
         IF (X .GT. 174.0D0) THEN
            K1=0.0D0
         ELSE
            T=2.0D0/X
            A1=1.25331414D0
            A2=0.23498619D0
            A3=-0.03655620D0
            A4=0.01504268D0
            A5=-0.00780355D0
            A6=0.00325614D0
            A7=-0.00068245D0
            K1=(A1+(A2+(A3+(A4+(A5+(A6+A7*T)*T)*T)*T)*T)*T)*
     '          DEXP(-X)/DSQRT(X)
         ENDIF
      ENDIF
      RETURN
      END


