      REAL*8 FUNCTION B_COEFF(X)

C#### Function: B_COEFF
C###  Description:
C###    B_COEFF is the coefficient (B) 0f Uxt in the PDE.
C###    This function must be altered whenever solving a new PDE.

      IMPLICIT NONE
!     Parameter List
      REAL*8 X

      B_COEFF=1.0d0-2.0d0*X

      RETURN
      END


