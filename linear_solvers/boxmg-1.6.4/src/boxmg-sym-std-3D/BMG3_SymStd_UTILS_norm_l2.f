      SUBROUTINE BMG3_SymStd_UTILS_norm_l2( 
     &                u, Nx, Ny, Nz, l2norm
     &                )

C ==========================================================================
C
C   BMG3_SymStd_UTILS_norm_l2.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_UTILS_norm_l2.f computes the l2-norm of a grid function
C   (the vector u).  It assumes that ghost points should be neglected.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   2000/03/06  - written (JDM)
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
      REAL*8   l2norm, u(Nx,Ny,Nz)

C ----------------------------
C     Local Declarations
C
      INTEGER  i, j, k, kk, maxXY, maxXYZ, Nxm2, Nym2, Nzm2

C =========================================================================

      l2norm=rZERO

      Nxm2 = Nx-2
      Nym2 = Ny-2
      Nzm2 = Nz-2

      maxXYZ = Nxm2*Nym2*Nzm2
      maxXY  = Nxm2*Nym2

C$OMP PARALLEL DO REDUCTION(+:l2norm)
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,j,k,kk)
C$OMP& SHARED(maxXY,maxXYZ,Nxm2,Nym2,Nzm2,u)
      DO kk = 0,maxXYZ-1
         i = mod(kk,Nxm2)+2
         j = mod(kk/Nxm2,Nym2)+2
         k = (kk/maxXY)+2
         l2norm=l2norm+u(i,j,k)*u(i,j,k)
      ENDDO
C$OMP END PARALLEL DO

      l2norm=SQRT(l2norm)

C =========================================================================

      RETURN
      END
