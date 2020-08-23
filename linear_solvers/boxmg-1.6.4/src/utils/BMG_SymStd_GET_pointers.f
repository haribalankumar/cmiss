      SUBROUTINE BMG_SymStd_GET_pointers(
     &                Ndim, k, BMG_iWORK, BMG_pWORK, NOG, Nx, Ny, Nz, 
     &                p_U, p_SO, p_SOR, p_CI
     &                )

C ===========================================================================
C
C   BMG_SymStd_GET_pointers.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG_SymStd_GET_pointers.f retrieves pointers to data at grid
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
      INTEGER k, Ndim, Nx, Ny, Nz, NOG, p_CI, p_SO, p_SOR, p_U

      INTEGER BMG_iWORK(*), BMG_pWORK(NBMG_pWORK)

C ===========================================================================


      
      IF( Ndim.EQ.2 ) THEN
         Nz = 1
         CALL BMG2_SymStd_GET_pointers(
     &                k, BMG_iWORK(BMG_pWORK(ip_iG)), NOG, 
     &                Nx, Ny, p_U, p_SO, p_SOR, p_CI
     &                )
      ELSEIF( Ndim.EQ.3 ) THEN
         CALL BMG3_SymStd_GET_pointers(
     &                k, BMG_iWORK(BMG_pWORK(ip_iG)), NOG, 
     &                p_U, p_SO, p_SOR, p_CI, Nx, Ny, Nz
     &                )
      ENDIF
         
C ===========================================================================

      RETURN
      END
