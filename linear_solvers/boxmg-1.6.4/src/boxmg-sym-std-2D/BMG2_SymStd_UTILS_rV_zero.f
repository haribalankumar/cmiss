      SUBROUTINE BMG2_SymStd_UTILS_rV_zero( v, Nx, Ny )

      IMPLICIT NONE

      INCLUDE 'BMG_constants.h'

      INTEGER Nx, Ny
      REAL*8  v(Nx*Ny)

      INTEGER i

      DO i=1, Nx*Ny
         v(i) = rZERO
      END DO

      RETURN
      END

