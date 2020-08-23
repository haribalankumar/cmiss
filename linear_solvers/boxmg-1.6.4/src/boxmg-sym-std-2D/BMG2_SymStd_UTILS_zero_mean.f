      SUBROUTINE BMG2_SymStd_UTILS_zero_mean( x, Nx, Ny )

C ==========================================================================
C
C   BMG2_SymStd_UTILS_zero_mean
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_SymStd_UTILS_zreo_mean.f replaces x with a*x.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2004/12/22 (TMA)
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
      REAL*8   x(Nx,Ny)

C ----------------------------
C     Local Declarations
C
      INTEGER  i, j, kk
      REAL*8   C

C =========================================================================

      C = 0.D0
C$OMP PARALLEL DO REDUCTION(+:C)
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,j,kk)
C$OMP& SHARED(Nx,Ny,x)
	DO kk = 1,(Nx-2)*(Ny-2)
	   i = mod(kk-1,Nx-2)+2
	   j = (kk-1)/(Nx-2)+2
           C=C+x(i,j)
        END DO
C$OMP END PARALLEL DO

      C = C /((Nx-2)*(Ny-2))

C$OMP PARALLEL DO 
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,j,kk)
C$OMP& SHARED(Nx,Ny,x,C)
	DO kk = 1,(Nx-2)*(Ny-2)
	   i = mod(kk-1,Nx-2)+2
	   j = (kk-1)/(Nx-2)+2
           x(i,j)=x(i,j)-C
      END DO
C$OMP END PARALLEL DO


C =========================================================================

      RETURN
      END
