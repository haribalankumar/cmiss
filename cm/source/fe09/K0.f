      REAL*8 FUNCTION K0(X)

C#### Function: K0
C###  Type: REAL*8
C###  Description:
C###    K0 calculates the modified Bessel function of the second kind of
C###    order 1 using the approximation of Abromowitz and Stegun.

      IMPLICIT NONE
!     Parameter Values
      REAL*8 X
!     Local Variables
      REAL*8 I0
      REAL*8 A1,A2,A3,A4,A5,A6,A7,T

      !Calculates K0(x)
      IF ( X .LE. 2.0D0) THEN
         T=X*X/4.0D0
         A1=-0.57721566D0
         A2=0.42278420D0
         A3=0.23069756D0
         A4=0.03488590D0
         A5=0.00262698D0
         A6=0.00010750D0
         A7=0.00000740D0
         K0=-DLOG(X/2.0D0)*I0(X)+
     '       (A1+(A2+(A3+(A4+(A5+(A6+A7*T)*T)*T)*T)*T)*T)
      ELSE
         IF (X .GT. 174.0D0) THEN
            K0=0.0D0
         ELSE
            T=2.0D0/X
            A1=1.25331414D0
            A2=-0.07832358D0
            A3=0.02189568D0
            A4=-0.01062446D0
            A5=0.00587872D0
            A6=-0.00251540D0
            A7=0.00053208D0
            K0=(A1+(A2+(A3+(A4+(A5+(A6+A7*T)*T)*T)*T)*T)*T)*
     '          DEXP(-X)/DSQRT(X)
         ENDIF
      ENDIF
      RETURN
      END


