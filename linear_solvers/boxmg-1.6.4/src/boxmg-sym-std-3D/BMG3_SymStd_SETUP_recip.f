      SUBROUTINE BMG3_SymStd_SETUP_recip( 
     &                SO, SOR, Nx, Ny, Nz, NStncl, NSORv
     &                )

C ==========================================================================
C
C   BMG3_SymStd_SETUP_recip.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_SymStd_SETUP_recip.f computes the reciprocal of the 
C   central stencil coefficient for use in GS relaxation.  
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2000/11/14 (JDM)
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

C ---------------------------
C    Argument Declarations:
C
      INTEGER   NSORv, NStncl, Nx, Ny, Nz
      REAL*8    SO(Nx,Ny,Nz,NStncl), SOR(Nx,Ny,Nz,NSORv)

C --------------------------
C     Local Declarations:
C
      INTEGER   i, j, k, kk

C ==========================================================================

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,j,k,kk)
C$OMP& SHARED(Nx,Ny,Nz,SO,SOR)
      DO kk = 1, (Nx-2)*(Ny-2)*(Nz-2)
         i = mod((kk-1),Nx-2)+2
         j = mod((kk-1)/(Nx-2),Ny-2)+2
         k = (kk-1)/((Nx-2)*(Ny-2))+2
         SOR(i,j,k,msor)=rONE/SO(i,j,k,kp)
      ENDDO
C$OMP END PARALLEL DO

C ==========================================================================

C =======================

         RETURN
         END
