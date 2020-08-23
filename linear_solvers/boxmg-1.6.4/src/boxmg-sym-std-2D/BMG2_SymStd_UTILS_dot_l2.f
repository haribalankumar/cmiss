      SUBROUTINE BMG2_SymStd_UTILS_dot_l2( u, v, Nx, Ny, l2dot )

C ==========================================================================
C
C   BMG2_SymStd_dot_l2.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_SymStd_dot_l2.f computes the dot product of a grid function
C   u with grid function v.  It assumes that ghost points should be 
C   neglected.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2003/07/03 (JDM)
C
C
C ==================================================================
C   INPUT:
C ========================
C
C
C
C ==================================================================
C   OUTPUT:
C ===========================
C
C
C
C ==================================================================
C   LOCAL:
C ========================
C
C
C
C ==========================================================================

      IMPLICIT NONE

C -----------------------------
C     Includes
C
      INCLUDE 'BMG_constants.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER  Nx, Ny
      REAL*8   l2dot, u(Nx,Ny), v(Nx,Ny)

C ----------------------------
C     Local Declarations
C
      INTEGER  i, j

C =========================================================================

      l2dot=rZERO

      DO j = 2, Ny-1
         DO i = 2, Nx-1
            l2dot=l2dot+u(i,j)*v(i,j)
         END DO
      END DO


C =========================================================================

      RETURN
      END
