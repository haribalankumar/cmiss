      SUBROUTINE BMG2_CONVERT_soln( 
     &               Type, CM_IJK_MAP, BMG_TO_CM, Q, X, NX, NY 
     &               )

C =======================================================================
C 
C   BMG2_CONVERT_soln 
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_CONVERT_soln converts a standard 2D SOLN array into a BOXMG 
C   (padded) 2D SOLN array. 
C
C   Written:    2005/01/22 (TMA)
C
C   Author:  Austin, Travis
C   EMAIL:   t.austin@auckland.ac.nz
C
C =======================================================================
C   INPUT:
C ========================
C
C   ----------------------
C   Fine Grid Dimensions:
C   ---------------------- 
C
C   NX       Number of points in x-direction (excluding ghost points)
C   NY       Number of points in y-direction (excluding ghost points)
C
C   -----------------------
C   BOXMG Solution Array:
C   -----------------------
C
C   Q        Solution array from boxmg with extra padding.
C
C =======================================================================
C   OUTPUT:
C ========================
C
C   X        Standard CMISS solution array.
C
C ==========================================================================
      
      IMPLICIT NONE

C ========================================================================== 

C ---------------- 
C     Arguments  
C
      INTEGER  NX, NY, CM_IJK_MAP(*), BMG_TO_CM(*), Type
      REAL*8   Q(NX+2,NY+2), X(NX*NY) 

C ---------------- 
C     Arguments  
C
      INTEGER CMKK, I, J, KK

C ========================================================================== 

      IF( Type.EQ.0 ) THEN

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(CMKK,I,J,KK)
C$OMP& SHARED(BMG_TO_CM,CM_IJK_MAP,NX,NY,Q,X)
         DO KK = 1, NX*NY
            !
            I = CM_IJK_MAP(2*KK-1)
            J = CM_IJK_MAP(2*KK  )
            !
            CMKK = BMG_TO_CM(KK)
            !
            Q(I+1,J+1) = X(CMKK) 
            !
         ENDDO
C$OMP END PARALLEL DO

      ELSEIF( Type.EQ.1 ) THEN

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(CMKK,I,J,KK)
C$OMP& SHARED(BMG_TO_CM,CM_IJK_MAP,NX,NY,Q,X)
         DO KK = 1, NX*NY
            !
            I = CM_IJK_MAP(2*KK-1)
            J = CM_IJK_MAP(2*KK  )
            !
            CMKK = BMG_TO_CM(KK)
            !
            X(CMKK) = Q(I+1,J+1)
            !
         ENDDO
C$OMP END PARALLEL DO

      ENDIF

      RETURN
      END
