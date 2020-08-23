      SUBROUTINE BMG3_CONVERT_soln( 
     &               Type, CM_IJK_MAP, BMG_TO_CM, Q, X, NX, NY, NZ 
     &               )

C =======================================================================
C 
C   BMG3_CONVERT_soln 
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_CONVERT_soln converts a standard 3D SOLN array into a BOXMG 
C   (padded) 3D SOLN array. 
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
C   NZ       Number of points in z-direction (excluding ghost points)
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
      INTEGER  CM_IJK_MAP(*), BMG_TO_CM(*), NX, NY, NZ, Type
      REAL*8   Q(NX+2,NY+2,NZ+2), X(NX*NY*NZ) 

C ---------------- 
C     Arguments  
C
      INTEGER CMKK, I, J, K, KK, NXNY, NXNYNZ

C ========================================================================== 

      NXNYNZ = NX*NY*NZ
      NXNY   = NX*NY
      
      IF( Type.EQ.0 ) THEN

C$OMP PARALLEL DO 
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(CMKK,I,J,K,KK)
C$OMP& SHARED(CM_IJK_MAP,BMG_TO_CM,NXNYNZ,Q,X)
         DO KK = 1,NXNYNZ
            !
            I = CM_IJK_MAP(3*KK-2)
            J = CM_IJK_MAP(3*KK-1)
            K = CM_IJK_MAP(3*KK  )
            !
            CMKK = BMG_TO_CM(KK)
            !
            Q(I+1,J+1,K+1) = X(CMKK) 
            !
         ENDDO
C$OMP END PARALLEL DO

      ELSEIF( Type.EQ.1 ) THEN

C$OMP PARALLEL DO 
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(CMKK,I,J,K,KK)
C$OMP& SHARED(BMG_TO_CM,CM_IJK_MAP,NXNYNZ,Q,X)
         DO KK = 1,NXNYNZ
            !
            I = CM_IJK_MAP(3*KK-2)
            J = CM_IJK_MAP(3*KK-1)
            K = CM_IJK_MAP(3*KK  )
            !
            CMKK = BMG_TO_CM(KK)
            !
            X(CMKK) = Q(I+1,J+1,K+1) 
            !
         ENDDO
C$OMP END PARALLEL DO

      ENDIF

      RETURN
      END
