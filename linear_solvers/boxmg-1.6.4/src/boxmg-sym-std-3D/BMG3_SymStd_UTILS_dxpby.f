      SUBROUTINE BMG3_SymStd_UTILS_dxpby( beta, x, y, Nx, Ny, Nz)

C ==========================================================================
C
C   BMG3_SymStd_UTILS_dxpby.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_UTILS_dxpby.f replaces y with x + beta*y.
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
      REAL*8   beta, x(Nx,Ny,Nz), y(Nx,Ny,Nz)

C ----------------------------
C     Local Declarations
C
      INTEGER  i, j, k, kk

C =========================================================================

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,j,k,kk)
C$OMP& SHARED(beta,Nx,Ny,Nz,x,y)
      DO kk = 1, (Nx-2)*(Ny-2)*(Nz-2)
         i = mod((kk-1),(Nx-2))+2
         j = mod((kk-1)/(Nx-2),(Ny-2))+2
         k = (kk-1)/((Nx-2)*(Ny-2))+2
         y(i,j,k)=beta*y(i,j,k) + x(i,j,k)
      ENDDO
C$OMP END PARALLEL DO

C =========================================================================

      RETURN
      END
