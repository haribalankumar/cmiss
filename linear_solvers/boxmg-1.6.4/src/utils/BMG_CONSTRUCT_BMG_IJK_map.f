      SUBROUTINE BMG_CONSTRUCT_BMG_IJK_map( 
     &              Nd, NOG, IGRD, NOGm, BMG_IJK_MAP
     &              )

C =======================================================================
C 
C   BMG_CONSTRUCT_IJK_map
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG_CONSTRUCT_IJK_map creates map from KK loop to (I,J,K).
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
      INTEGER  Nd, NOG, NOGm
      INTEGER  IGRD, BMG_IJK_MAP(*)

      INTEGER II, JJ, KK, I, J, K, COUNT, GRID
      INTEGER p_U, p_SO, p_SOR, p_CI, iSTART

C ========================================================================== 


      iSTART = 20

      IF(Nd.EQ.2) THEN

         DO GRID=1,20
            BMG_IJK_MAP(GRID) = 0
         ENDDO

         COUNT = 0
         DO GRID=NOG,1,-1

            CALL BMG2_SymStd_GET_pointers( 
     &           GRID, IGRD, NOGm, II, JJ, p_U, p_SO, p_SOR, p_CI
     &           )    

            BMG_IJK_MAP(GRID) = COUNT + 21

            DO J=1,JJ
               DO I=1,II
                  IF((I.GT.1.AND.I.LT.II).AND.
     &               (J.GT.1.AND.J.LT.JJ)) THEN
                     COUNT = COUNT + 1
                     BMG_IJK_MAP(iSTART+COUNT) = I
                     COUNT = COUNT + 1
                     BMG_IJK_MAP(iSTART+COUNT) = J
                  ENDIF
               ENDDO
            ENDDO
            
         ENDDO

    
      ELSEIF(Nd.EQ.3) THEN

         DO GRID=1,20
            BMG_IJK_MAP(GRID) = 0
         ENDDO

         COUNT = 0
         DO GRID=NOG,1,-1
            !
            CALL BMG3_SymStd_GET_pointers( 
     &           GRID, IGRD, NOGm, p_U, p_SO, p_SOR, 
     &           p_CI, II, JJ, KK 
     &           )    
            !
            BMG_IJK_MAP(GRID) = COUNT + 21
            !
            DO K=1,KK
               DO J=1,JJ
                  DO I=1,II
                     IF((I.GT.1.AND.I.LT.II).AND.
     &                  (J.GT.1.AND.J.LT.JJ).AND.
     &                  (K.GT.1.AND.K.LT.KK)) THEN
                        !
                        COUNT = COUNT + 1
                        BMG_IJK_MAP(iSTART+COUNT) = I
                        COUNT = COUNT + 1
                        BMG_IJK_MAP(iSTART+COUNT) = J
                        COUNT = COUNT + 1
                        BMG_IJK_MAP(iSTART+COUNT) = K
                        !
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO

         ENDDO
                        
      ENDIF

      RETURN
      END
