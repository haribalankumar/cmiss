      SUBROUTINE BMG3_SymStd_UTILS_diag_scal( alpha, x, Nx, Ny, Nz)

C ==========================================================================
C
C   BMG3_SymStd_UTILS_dscal.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_UTILS_dscal.f replaces x with a*x.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2003/07/03 (TMA)
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
      INTEGER  Nx, Ny, Nz
      REAL*8   alpha(Nx,Ny,Nz), x(Nx,Ny,Nz)

C ----------------------------
C     Local Declarations
C
      INTEGER  i, j, k, kk, maxXY, maxXYZ, Nxm2, Nym2, Nzm2


C =========================================================================

      Nxm2 = Nx-2
      Nym2 = Ny-2
      Nzm2 = Nz-2

      maxXYZ = Nxm2*Nym2*Nzm2
      maxXY  = Nxm2*Nym2

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,j,k,kk)
C$OMP& SHARED(alpha,maxXY,maxXYZ,Nxm2,Nym2,Nzm2,x)
      DO kk = 0,maxXYZ-1
         i = mod(kk,Nxm2)+2
         j = mod(kk/Nxm2,Nym2)+2
         k = (kk/maxXY)+2
         x(i,j,k) = alpha(i,j,k)*x(i,j,k)
      ENDDO
C$OMP END PARALLEL DO

C =========================================================================

      RETURN
      END
