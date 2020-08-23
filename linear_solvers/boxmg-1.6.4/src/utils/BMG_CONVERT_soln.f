      SUBROUTINE BMG_CONVERT_soln(
     &              Type, Nd, Nx, Ny, Nz, 
     &              BMG_MAP, BMG_TO_CM, 
     &              BMG_rWORK, BMG_pWORK, X
     &              )

C =======================================================================
C 
C   BMG_CONVERT_rhs 
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG_CONVERT_rhs converts a standard RHS array into a BOXMG 
C   (padded) RHS array. 
C
C   Written:    2005/01/22 (TMA)
C
C =======================================================================
C   INPUT:
C ========================
C
C   ----------------------
C   Fine Grid Dimensions:
C   ---------------------- 
C
C   Nd       Dimension of problem
C   Nx       Number of points in x-direction (excluding ghost points)
C   Ny       Number of points in y-direction (excluding ghost points)
C   Nz       Number of points in z-direction (excluding ghost points)
C  Type      Tells if call is before BOXMG solve call or after
C
C
C   -----------------------
C   BOXMG Parameter Arrays:
C   -----------------------
C
C   BMG_rWORK      Arrays containing solution, rhs, matrix, residual, etc.
C   BMG_pWORK      Array containing pointers to locations in rWORK and iWORK.
C
C =======================================================================
C   LOCAL:
C ========================
C
C ==========================================================================
      
      IMPLICIT NONE
      INCLUDE 'BMG_workspace.h'

C ========================================================================== 

C ---------------- 
C     Arguments  
C
      INTEGER  Nd, Nx, Ny, Nz, BMG_MAP(*), BMG_TO_CM(*), Type
      INTEGER  BMG_pWORK(*)
      REAL*8   BMG_rWORK(*), X(*)

C ========================================================================== 


      IF( Nd.EQ.2 ) THEN

         CALL BMG2_CONVERT_soln( 
     &           Type, BMG_MAP, BMG_TO_CM, 
     &           BMG_rWORK(BMG_pWORK(ip_U)), X, Nx, Ny 
     &           )

      ELSEIF( Nd.EQ.3 ) THEN

         CALL BMG3_CONVERT_soln( 
     &           Type, BMG_MAP, BMG_TO_CM,
     &           BMG_rWORK(BMG_pWORK(ip_U)), X, Nx, Ny, Nz 
     &           )

      ENDIF


      RETURN
      END
