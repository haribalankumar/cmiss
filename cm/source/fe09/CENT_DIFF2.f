      REAL*8 FUNCTION CENT_DIFF2(x1,x2,xstar,delta1,delta2)

C#### Function: CENT_DIFF2
C###  Type: REAL*8
C###  Description:
C###    <HTML>
C###       Calculates second order central difference
C###       at the node xstar surrounded by nodes x1 and x2.
C###       <B>NOTE:</B> a non-uniform mesh is assumed.
C###    </HTML>

      REAL*8 x1,x2,xstar,delta1,delta2

      CENT_DIFF2=(x2- 2.D0*xstar + x1)*2.D0/((delta1**2)+(delta2**2))

      RETURN
      END




