      REAL*8 FUNCTION C_COEFF(X)

C#### Function: C_COEFF
C###  Description:
C###    C_COEFF is the coefficient (C) 0f Utt in the PDE.
C###    This function must be altered whenever solving a new PDE.

      IMPLICIT NONE
!     Parameter List
      REAL*8 X

      C_COEFF=X*X-X-2.0d0

      RETURN
      END


