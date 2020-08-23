      REAL*8 FUNCTION A_COEFF()

C#### Function: A_COEFF
C###  Description:
C###    A_COEFF is the coefficient (A) 0f Uxx in the PDE.
C###    This function must be altered whenever solving a new PDE.

      IMPLICIT NONE
!     Parameter List

      A_COEFF=1.0d0

      RETURN
      END

