      SUBROUTINE BMG2_SymStd_UTILS_norm_l2( u, Nx, Ny, l2norm )

C ==========================================================================
C
C   BMG2_SymStd_norm_l2.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_SymStd_norm_l2.f computes the l2-norm of a grid function
C   (the vector u).  It assumes that ghost points should be neglected.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2000/02/17 (JDM)
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
      REAL*8   l2norm, u(Nx,Ny)

C ----------------------------
C     Local Declarations
C
      INTEGER  i, j, kk

C =========================================================================

      l2norm=rZERO

CC$OMP PARALLEL DO REDUCTION(+:l2norm) 
CC$OMP& DEFAULT(NONE)
CC$OMP& PRIVATE(i,j,kk)
CC$OMP& SHARED(Nx,Ny,u)
	DO kk = 1,(Nx-2)*(Ny-2)
	   i = mod(kk-1,Nx-2)+2
	   j = (kk-1)/(Nx-2)+2
           l2norm=l2norm+u(i,j)*u(i,j)
      END DO
CC$OMP END PARALLEL DO

      l2norm=SQRT(l2norm)

C =========================================================================

      RETURN
      END
