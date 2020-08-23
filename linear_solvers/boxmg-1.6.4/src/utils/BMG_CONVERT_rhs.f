      SUBROUTINE BMG_CONVERT_rhs( 
     &              Nd, Nx, Ny, Nz, 
     &              BMG_rWORK, BMG_pWORK, 
     &              BMG_iPARMS, BMG_BDY, BMG_MAP, BMG_TO_CM, 
     &              AA, IA, JA, B 
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
      INCLUDE 'BMG_parameters.h'
      
C ========================================================================== 

C ---------------- 
C     Arguments  
C
      INTEGER  Nd, Nx, Ny, Nz
      INTEGER  BMG_BC, BMG_BDY(*), BMG_MAP(*), BMG_pWORK(*),
     &         BMG_iPARMS(*), BMG_TO_CM(*), IA(*), JA(*) 
      REAL*8   AA(*), B(*), BMG_rWORK(*)

C ========================================================================== 

      IF(Nd.EQ.2) THEN
         BMG_BC = BMG_iPARMS(id_BMG2_BC) 
      ELSEIF(Nd.EQ.3) THEN
         BMG_BC = BMG_iPARMS(id_BMG3_BC) 
      ENDIF

      IF( Nd.EQ.2 ) THEN

         CALL BMG2_CONVERT_rhs( 
     &           BMG_rWORK(BMG_pWORK(ip_Q)),BMG_BC,BMG_BDY,BMG_MAP,
     &           BMG_TO_CM,AA,IA,JA,B,Nx,Ny 
     &           )

      ELSEIF( Nd.EQ.3 ) THEN

         CALL BMG3_CONVERT_rhs( 
     &           BMG_rWORK(BMG_pWORK(ip_Q)),BMG_BC,BMG_BDY,BMG_MAP,
     &           BMG_TO_CM,AA,IA,JA,B,Nx,Ny,Nz 
     &           )

      ENDIF


      RETURN
      END
