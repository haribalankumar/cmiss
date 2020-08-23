      SUBROUTINE BMG2_CONVERT_csr_matrix( 
     &               NX, NY,
     &               SO, BMG_BDY, BMG_iPARMS, BMG_TO_CM,
     &               AA, IA, JA, NStncl
     &               )

C ==========================================================================
C
C***BEGIN PROLOGUE  BMG2_CONVERT_csr_matrix
C
C***PURPOSE  
C 
C     This subroutine converts a matrix stored in CSR format to a matrix
C     stored in the Dendy-Cell-Symmetric format for 2D structured grids.  It
C     is assumed that the matrix has been generated from a lexicographical
C     construction so that the unknowns may be accessed via (ij) entries.
C
C***AUTHOR  AUSTIN, TRAVIS
C           E-MAIL: t.austin@auckland.ac.nz
C
C***DESCRIPTION
C
C***PARAMETERS
C
C***INPUT
C
C   NStncl    Determines what type of stencil we have.
C     NX      x-dimension of the grid, excluding fictitious points
C     NY      y-dimension of the grid, excluding fictitious points
C
C     AA      Nonzeros in matrix stored in CSR format.
C     IA      Pointers to rows in matrix A.
C     JA      Column number for corresponding entry in AA 
C
C   Note:  Size(AA) = Size(JA)
C
C***OUTPUT
C
C    SO       Matrix stored in Dendy-Cell-Symmetric format
C
C ==========================================================================

      IMPLICIT NONE

C -----------------------------
C     Includes
C
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_parameters.h'

C ========================================================================== 

C ---------------- 
C     Arguments  
C
      INTEGER  NX, NY, NStncl 
      INTEGER  BMG_BDY(*),BMG_iPARMS(*), BMG_TO_CM(*), IA(*), JA(*) 
      REAL*8   SO(NX+2,NY+2,NStncl), AA(*)

C ---------------- 
C     Local 
C 
      INTEGER  BC, IT, I, J, KL, IJ, IM1J, IJM1, IM1JM1, IM1JP1, N, SGN
      INTEGER  CMIJ, CMIM1J, CMIJM1, CMIM1JM1, CMIM1JP1, CMKL, NXNY
      REAL*8   FACT, SUM
C ========================================================================== 



C ==================================================
C
C  Zero out matrix.
C
C =================================================

      NXNY = NX*NY

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(I,J,N)
C$OMP& SHARED(NX,NY,NStncl,SO)
      DO KL=1,(NX+2)*(NY+2)
	 !
         I = MOD(KL-1,(NX+2)) + 1
         J = (KL-1)/(NX+2) + 1
         !
         DO N=1,NStncl
            SO(I,J,N) = 0.D0
         ENDDO
         !
      ENDDO
C$OMP END PARALLEL DO
         
C ==================================================
C
C  Determine whether or not matrix is indefinite.
C
C =================================================

      BC = BMG_iPARMS(id_BMG2_BC)

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(CMKL,J,SUM)
C$OMP& SHARED(AA,BC,BMG_BDY,BMG_iPARMS,BMG_TO_CM,IA,NXNY)
      DO KL=1,NXNY
         !
         SUM = 0.D0
         !
         CMKL = BMG_TO_CM(KL)  ! Get CM numbering for grid point
         !
         DO J=IA(CMKL),IA(CMKL+1)-1
            SUM = SUM + AA(J)
         ENDDO
         IF( DABS(SUM).GT.1E-10 ) THEN
             BMG_iPARMS(id_BMG2_BC) = ABS(BC)
         ENDIF
         !
         BMG_BDY(KL) = 0
         !
      ENDDO
C$OMP END PARALLEL DO

      !
      ! See boxmg-1.6.4/include/BMG_parameters.h for insight into 
      ! this line (we are setting appropriate definite bcs)
      !
      IF( BMG_iPARMS(id_BMG2_BC).EQ.5 ) BMG_iPARMS(id_BMG2_BC) = 0

C =================================================
C
C  Determine which rows correspond to Dirichlet
C  boundary points and also which have a negative
C  diagonal value.
C
C =================================================

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(CMKL,I,IJ,IJM1,IM1J,IM1JM1,IT,J,KL)
C$OMP& SHARED(AA,BMG_BDY,BMG_TO_CM,IA,JA,NX,NY,NXNY,SGN,SO)
      DO KL=1,NXNY
         !
         CMKL = BMG_TO_CM(KL)  ! Get CM numbering for grid point
         !
         DO IT=IA(CMKL),IA(CMKL+1)-1
            !
            ! Check sign of diagonal 
            !
            IF( JA(IT).EQ.CMKL ) THEN 
               IF( AA(IT).LE.0.D0 ) THEN
                  SGN = -1      ! Negative diagonal
               ELSE
                  SGN = 1       ! Positive diagonal
               ENDIF
            ENDIF
            !
         ENDDO
         !
         ! Check if all off-diagonals are zero
         !
         DO IT=IA(CMKL),IA(CMKL+1)-1
           !
            IF( ABS(AA(IT)).GE.1E-8 .AND. JA(IT).NE.CMKL ) THEN
               BMG_BDY(KL) = SGN ! Zero if dirichlet bc
            ENDIF 
            !
         ENDDO
         !
      ENDDO	
C$OMP END PARALLEL DO

C =================================================
C
C  Fill-up BOXMG matrix.
C
C =================================================
      
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(CMKL,CMIJ,CMIM1J,CMIJM1,CMIM1JM1,CMIM1JP1,  
C$OMP&         FACT,I,IJ,IJM1,IM1J,IM1JM1,IM1JP1,IT,J,KL)
C$OMP& SHARED(AA,BMG_BDY,BMG_TO_CM,IA,JA,NX,NY,NXNY,SO)
      DO KL=1,NXNY
	 !
         I = MOD(KL-1,NX) + 1
         J = (KL-1)/NX + 1
	 !
         IJ = (J-1)*NX + I
         IM1J   = IJ-1
         IJM1   = IJ-NX
         IM1JM1 = IJ-NX-1 
         IM1JP1 = IJ+NX-1
	 !
         ! Get corresponding CM numbering for grid points
         !
         CMIJ     = 0
         CMIM1J   = 0
         CMIJM1   = 0
         CMIM1JM1 = 0
         CMIM1JP1 = 0
         !
         IF( IJ.GT.0.AND.IJ.LT.NXNY+1 ) 
     &        CMIJ   = BMG_TO_CM(IJ)
         IF( IM1J.GT.0.AND.IM1J.LT.NXNY+1 ) 
     &        CMIM1J = BMG_TO_CM(IM1J)
         IF( IJM1.GT.0.AND.IJM1.LT.NXNY+1 ) 
     &        CMIJM1 = BMG_TO_CM(IJM1)
         IF( IM1JM1.GT.0.AND.IM1JM1.LT.NXNY+1 ) 
     &        CMIM1JM1 = BMG_TO_CM(IM1JM1)
         IF( IM1JP1.GT.0.AND.IM1JP1.LT.NXNY+1 ) 
     &        CMIM1JP1 = BMG_TO_CM(IM1JP1)
         
         !
         IF( BMG_BDY(KL).NE.0 ) THEN
            !
            !
            DO IT=IA(CMIJ),IA(CMIJ+1)-1
               !
               IF( JA(IT).EQ.CMIJ ) THEN
                  !
                  SO(I+1,J+1,ko) = BMG_BDY(IJ)*AA(IT)
                  !
               ELSEIF( I.NE.1 .AND. JA(IT).EQ.CMIM1J ) THEN
                  !
                  FACT=ABS(BMG_BDY(IM1J))*BMG_BDY(IJ)
                  SO(I+1,J+1,kw) = -FACT*AA(IT)
                  !
               ELSEIF( J.NE.1 .AND. JA(IT).EQ.CMIJM1 ) THEN
                  !
                  FACT=ABS(BMG_BDY(IJM1))*BMG_BDY(IJ)
                  SO(I+1,J+1,ks) = -FACT*AA(IT)
                  !
               ELSEIF( (I.NE.1 .OR. J.NE.1) 
     &                 .AND. JA(IT).EQ.CMIM1JM1 ) THEN
                  !
                  FACT=ABS(BMG_BDY(IM1JM1))*BMG_BDY(IJ)
                  SO(I+1,J+1,ksw) = -FACT*AA(IT)
                  !
               ELSEIF( (I.NE.1 .OR. J.NE.NY)
     &                 .AND. JA(IT).EQ.CMIM1JP1 ) THEN
                  !
                  FACT=ABS(BMG_BDY(IM1JP1))*BMG_BDY(IJ)
                  SO(I+1,J+2,knw) = -FACT*AA(IT)
                  !
               ENDIF
               !
            ENDDO
            !
         ELSE
            !
            SO(I+1,J+1,ko)  = 1.D0
            SO(I+1,J+1,kw)  = 0.D0
            SO(I+1,J+1,ks)  = 0.D0
            SO(I+1,J+1,ksw) = 0.D0
            SO(I+1,J+2,knw) = 0.D0
            !
         ENDIF

      ENDDO	
C$OMP END PARALLEL DO

      RETURN
      END 
