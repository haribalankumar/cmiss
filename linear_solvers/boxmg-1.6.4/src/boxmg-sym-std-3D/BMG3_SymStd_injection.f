      SUBROUTINE BMG3_SymStd_injection(
     &                KFG, KCG, 
     &                Q, QC, CI, Nx, Ny, Nz, Nxc, Nyc, Nzc 
     &                )

C ==========================================================================
C
C   BMG3_SymStd_restrict.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_restrict computes the restriction of a vector on the
C   fine grid, Q, to a vector on the coarse grid, QC.  The weights 
C   involve the transpose of the interpolation operator from the coarse 
C   grid to the fine grid.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C
C   2000/03/06  - written (JDM)
C               - the core computation was cut from mgrcl3.f
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
      INCLUDE 'BMG_stencils.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER  Nx, Nxc, Ny, Nyc, Nz, Nzc

      INTEGER  KCG, KFG
      REAL*8   CI(Nxc,Nyc,Nzc,26), Q(Nx,Ny,Nz), QC(Nxc,Nyc,Nzc)

C ----------------------------
C     Local Declarations
C
      INTEGER  i, ic, j, jc, k, kc, kl, maxXY, maxXYZ, Nxc2, Nyc2, Nzc2

C ==========================================================================

C ---------------------------------------------
C     Restrict the vector Q -> QC:
C ---------------------------------------------

      Nxc2 = Nxc-2
      Nyc2 = Nyc-2
      Nzc2 = Nzc-2

      maxXYZ = Nxc2*Nyc2*Nzc2
      maxXY  = Nxc2*Nyc2

      PRINT *,' ****** INJECTION ******* '

C$OMP PARALLEL DO IF(maxXYZ.GT.100)
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,ic,j,jc,k,kc)
C$OMP& SHARED(maxXY,maxXYZ,Nxc2,Nyc2,Nzc2,CI,Q,QC)
      DO kl = 0,maxXYZ-1
         
         ic = MOD(kl,Nxc2)+2
         i = 2*(ic-1)
         jc = MOD(kl/Nxc2,Nyc2)+2
         j = 2*(jc-1)
         kc = (kl/maxXY)+2
         k = 2*(kc-1)

         QC(ic,jc,kc) = Q(i,j,k)
              
      ENDDO
C$OMP END PARALLEL DO
         
C ==========================================================================

      return
      end






