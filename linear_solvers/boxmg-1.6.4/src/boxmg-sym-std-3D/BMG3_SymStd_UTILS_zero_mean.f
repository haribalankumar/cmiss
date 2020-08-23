      SUBROUTINE BMG3_SymStd_UTILS_zeRo_mean( x, Nx, Ny, Nz)

C ==========================================================================
C
C   BMG3_SymStd_UTILS_zero_mean
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_UTILS_zreo_mean.f replaces x with x + c, c is a constant,
C     such that the integral of x + c over domain D is zero.
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
      INTEGER  Nx, Ny, Nz
      REAL*8   x(Nx,Ny,Nz)

C ----------------------------
C     Local Declarations
C
      INTEGER  i, j, k, kk, maxXY, maxXYZ, Nxm2, Nym2, Nzm2
      REAL*8   C

C =========================================================================

      C = 0.D0

      Nxm2 = Nx-2
      Nym2 = Ny-2
      Nzm2 = Nz-2

      maxXYZ  = Nxm2*Nym2*Nzm2
      maxXY = Nxm2*Nym2

C$OMP PARALLEL DO REDUCTION(+:C)
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,j,k,kk)
C$OMP& SHARED(maxXY,maxXYZ,Nxm2,Nym2,Nzm2,x)
      DO kk = 0, maxXYZ-1
         i = mod(kk,Nxm2)+2
         j = mod(kk/Nxm2,Nym2)+2
         k = (kk/maxXY)+2
         C = C + x(i,j,k)
      ENDDO
C$OMP END PARALLEL DO

      C = C/MaxXYZ

C$OMP PARALLEL DO 
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,j,k,kk)
C$OMP& SHARED(maxXY,maxXYZ,Nxm2,Nym2,Nzm2,x,C)
      DO kk = 0, maxXYZ-1
         i = mod(kk,Nxm2)+2
         j = mod(kk/Nxm2,Nym2)+2
         k = (kk/maxXY)+2
         x(i,j,k) = x(i,j,k) - C
      ENDDO
C$OMP END PARALLEL DO

C =========================================================================

      RETURN
      END
