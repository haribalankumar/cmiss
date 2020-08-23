      REAL*8 FUNCTION CONC_EVAL(T)

C#### Function: CONC_EVAL
C###  Description:
C###    CONC_EVAL evaluates the saturation concentration (g.mm-3)
C###    at temperature T, using the formula and coefficients
C###    given by Perry's Chemical Engineers' Handbook (7th ed)
C###  Eds RH Perry, DW Green, and JO Maloney, McGraw-Hill, 1997. pp 2
C###  -54, table 2-6 (Vapour pressure of inorganic and organic liquids).
C###    T is in K. CONC_EVAL is returned in mg/L.

      IMPLICIT NONE
!     Parameter List
      REAL*8 T
!     Local variables
      REAL*8 C(4)
      DATA C/73.649d0,-7258.2d0,-7.3037d0,4.1653d-6/

      CONC_EVAL=2.166417d-9/T
     '  *DEXP(C(1)+C(2)/T+C(3)*DLOG(T)+C(4)*T**2)
      CONC_EVAL=CONC_EVAL*1.d9 !using mg/L

      RETURN
      END


