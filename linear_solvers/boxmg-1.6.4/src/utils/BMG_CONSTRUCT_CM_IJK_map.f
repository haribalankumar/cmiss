      SUBROUTINE BMG_CONSTRUCT_CM_IJK_map( 
     &              Nd, Nx, Ny, Nz, BMG_IJK_MAP
     &              )

C =======================================================================
C 
C   BMG_CONSTRUCT_IJK_MAP
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG_CONSTRUCT_IJK_MAP creates MAP from KK loop to (I,J,K).
C
C   Written:    2005/05/30 (TMA)
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
      INTEGER BMG_IJK_MAP(*), Nd, Nx, Ny, Nz

C ---------------------
C     Local Variables 
C

      INTEGER I, J, K, KK, COUNT

C ========================================================================== 

      IF(Nd.EQ.2) THEN
         COUNT = 0
         DO J=1,Ny
            DO I=1,Nx
               COUNT = COUNT + 1
               BMG_IJK_MAP(COUNT) = I
               COUNT = COUNT + 1
               BMG_IJK_MAP(COUNT) = J
            ENDDO
         ENDDO
      ELSEIF(Nd.EQ.3) THEN
         COUNT = 0
         DO K=1,Nz
            DO J=1,Ny
               DO I=1,Nx
                  COUNT = COUNT + 1
                  BMG_IJK_MAP(COUNT) = I
                  COUNT = COUNT + 1
                  BMG_IJK_MAP(COUNT) = J
                  COUNT = COUNT + 1
                  BMG_IJK_MAP(COUNT) = K
               ENDDO
            ENDDO
         ENDDO

      ENDIF

      RETURN
      END
