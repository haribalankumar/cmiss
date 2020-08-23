      SUBROUTINE BMG2_SymStd_GET_pointers(
     &                k, IGRD, NOG, Nx, Ny, p_U, p_SO, p_SOR, p_CI
     &                )

C ===========================================================================
C
C   BMG2_SymStd_GET_pointers.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_SymStd_GET_pointers.f retrieves pointers to data at grid
C   level k within the internal pointer array IGRD.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2000/12/10 (JDM)
C
C ==================================================================
C
C
C ===========================================================================

      IMPLICIT NONE

C ---------------------------
C    Includes:
C
      INCLUDE 'BMG_workspace.h'

C ---------------------------
C    Argument Declarations:
C     
      INTEGER k, Nx, Ny, NOG

      INTEGER IGRD(NOG,9), p_CI, p_SO, p_SOR, p_U

C ===========================================================================

      Nx = IGRD(k,idL_BMG_Nx)
      Ny = IGRD(k,idL_BMG_Ny)

      p_U   = IGRD(k,ipL_BMG_U)
      p_SO  = IGRD(k,ipL_BMG_SO)
      p_SOR = IGRD(k,ipL_BMG_SOR)
      p_CI  = IGRD(k,ipL_BMG_CI)

C ===========================================================================

      RETURN
      END
