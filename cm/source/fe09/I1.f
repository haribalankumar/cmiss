      REAL*8 FUNCTION I1(X)

C#### Function: I1
C###  Type: REAL*8
C###  Description:
C###    I1 calculates the modified Bessel function of the first kind of
C###    order 1 using the approximation of Abromowitz and Stegun.

      IMPLICIT NONE
!     Parameter List
      REAL*8 X
!     Local Variables
      REAL*8 A1,A2,A3,A4,A5,A6,A7,T

      !Calculates I1(x)
      T=(X/3.75D0)*(X/3.75D0)
      A1=0.5D0
      A2=0.87890594D0
      A3=0.51498869D0
      A4=0.15084934D0
      A5=0.02658733D0
      A6=0.00301532D0
      A7=0.00032411D0
      I1=(A1+(A2+(A3+(A4+(A5+(A6+A7*T)*T)*T)*T)*T)*T)*X
      RETURN
      END


