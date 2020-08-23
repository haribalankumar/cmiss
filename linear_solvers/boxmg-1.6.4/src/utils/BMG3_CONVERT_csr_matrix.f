      SUBROUTINE BMG3_CONVERT_csr_matrix(
     &               NX, NY, NZ,
     &               SO, BMG_BDY, BMG_iPARMS, BMG_TO_CM,
     &               AA, IA, JA, NStncl
     &               )

C ==========================================================================
C
C***BEGIN PROLOGUE  BMG3_CONVERT_csr_matrix
C
C***PURPOSE  
C 
C     This subroutine converts a matrix stored in CSR format to a matrix
C     stored in the Dendy-Cell-Symmetric format for 3D structured grids.  It
C     is assumed that the matrix has been generated from a lexicographical
C     construction so that the unknowns may be accessed via (ijk) entries.
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
C     NZ      z-dimension of the grid, excluding fictitious points
C
C     AA      Nonzeros in matrix stored in CSR format.
C     IA      Pointers to rows in matrix A.
C     JA      Column number for corresponding entry in AA 
C
C   Note:  Size(AA) = Size(JA)
C
C
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
      INTEGER  NX, NY, NZ, NStncl 
      INTEGER  BMG_BDY(*), BMG_iPARMS(*), BMG_TO_CM(*), IA(*), JA(*) 
      REAL*8   SO(NX+2,NY+2,NZ+2,NStncl), AA(*) 

C ---------------- 
C     Local 
C 
      INTEGER  BC, CMIJK, CMIM1JK, CMIJM1K, CMIJKM1, CMIJP1KM1,
     &         CMIM1JM1K, CMIM1JKM1, CMIJM1KM1, 
     &         CMIM1JM1KM1, CMIM1JP1K, CMIP1JKM1,
     &         CMIM1JP1KM1, CMIP1JP1KM1, CMIP1JM1KM1, CMKL,
     &         IFD, IT, I, J, K, KL, 
     &         IJK, IM1JK, IJM1K, IJKM1, IJP1KM1,
     &         IM1JM1K, IM1JKM1, IJM1KM1, 
     &         IM1JM1KM1, IM1JP1K, IP1JKM1,
     &         IM1JP1KM1, IP1JP1KM1, IP1JM1KM1, 
     &         N, NXNY, NXNYNZ, NXP2, NYP2, NZP2, 
     &         NXP2NYP2, SGN
      REAL*8   FACT,SUM

C ========================================================================== 

      IFD = BMG_iPARMS(id_BMG3_STENCIL)

C ==================================================
C
C  Zero out matrix.
C
C =================================================

      NXP2 = NX+2
      NYP2 = NY+2
      NZP2 = NZ+2

      NXP2NYP2 = NXP2*NYP2

      NXNY   = NX*NY
      NXNYNZ = NX*NY*NZ


C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(I,J,K,N)
C$OMP& SHARED(NXP2,NYP2,NZP2,NXP2NYP2,NStncl,SO)
      DO KL=1,NXP2*NYP2*NZP2
	 !
         I = MOD((KL-1),NXP2)+1
         J = MOD((KL-1)/NXP2,NYP2)+1
         K = (KL-1)/NXP2NYP2+1
         !
         DO N=1,NStncl
            SO(I,J,K,N) = 0.D0
         ENDDO
         !
      ENDDO
C$OMP END PARALLEL DO


C ==================================================
C
C  Determine whether or not matrix is indefinite.
C
C =================================================

      BC = BMG_iPARMS(id_BMG3_BC)

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(CMKL,N,SUM)
C$OMP& SHARED(AA,BC,BMG_BDY,BMG_iPARMS,BMG_TO_CM,IA,NXNYNZ)
      DO KL=1,NXNYNZ
         !
         SUM = 0.D0
         !
         CMKL = BMG_TO_CM(KL)  ! Get CM numbering for grid point
         !
         DO N=IA(CMKL),IA(CMKL+1)-1
            SUM = SUM + AA(N)
         ENDDO
         IF( DABS(SUM).GT.1E-10 ) THEN
            BMG_iPARMS(id_BMG3_BC) = ABS(BC)
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
      IF( BMG_iPARMS(id_BMG3_BC).EQ.5 ) BMG_iPARMS(id_BMG3_BC) = 0

C =================================================
C
C  Determine which rows correspond to Dirichlet
C  boundary points and also which have a negative
C  diagonal value.
C
C =================================================

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(CMKL,IT,KL,SGN)
C$OMP& SHARED(AA,BMG_BDY,BMG_TO_CM,IA,JA,NXNYNZ,SO)
      DO KL=1,NXNYNZ
         !
         CMKL = BMG_TO_CM(KL) 
         !
         DO IT=IA(CMKL),IA(CMKL+1)-1
            !
            IF( JA(IT).EQ.CMKL ) THEN ! Check diagonal
               IF( AA(IT).LE.0.D0 ) THEN
                  SGN = -1      ! Negative diagonal
               ELSE
                  SGN =  1      ! Positive diagonal
               ENDIF
            ENDIF
            !
         ENDDO
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
      
      IF( IFD.EQ.2 ) THEN

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(CMIJK,CMIM1JK,CMIJM1K,CMIJKM1,CMIJM1KM1,
C$OMP&         CMIM1JKM1,CMIM1JM1K,CMIM1JM1KM1,CMIM1JP1K,
C$OMP&         CMIP1JKM1,CMIJP1KM1,CMIM1JP1KM1,CMIP1JP1KM1,
C$OMP&         CMIP1JM1KM1,FACT,I,IJK,IM1JK,IJM1K,IJKM1,
C$OMP&         IJM1KM1,IM1JKM1,IM1JM1K,IM1JM1KM1,IM1JP1K,
C$OMP&         IP1JKM1,IJP1KM1,IM1JP1KM1,IP1JP1KM1,IP1JM1KM1,
C$OMP&         IT,J,K,KL)
C$OMP& SHARED(AA,BMG_BDY,BMG_TO_CM,IA,JA,NX,NY,NXNY,NXNYNZ,SO)
         DO KL = 1,NXNYNZ
	    !
            I = MOD((KL-1),NX)+1
            J = MOD((KL-1)/NX,NY)+1
            K = (KL-1)/NXNY+1
	    !
            IJK         = (K-1)*NXNY+(J-1)*NX+I
            IM1JK       = IJK - 1
            IJM1K       = IJK - NX
            IJKM1       = IJK - NXNY
            IJM1KM1     = IJK - NX - NXNY
            IM1JKM1     = IJK - NXNY - 1
            IM1JM1K     = IJK - NX - 1
            IM1JM1KM1   = IJK - NX - NXNY - 1
	    !
            IM1JP1K     = IJK + NX - 1
            IP1JKM1     = IJK - NXNY + 1
            IJP1KM1     = IJK - NXNY + NX
            IM1JP1KM1   = IJK - NXNY + NX - 1
            IP1JP1KM1   = IJK - NXNY + NX + 1
            IP1JM1KM1   = IJK - NXNY - NX + 1
            !
            CMIJK       = 0
            CMIM1JK     = 0
            CMIJM1K     = 0
            CMIJKM1     = 0
            CMIJM1KM1   = 0
            CMIM1JKM1   = 0
            CMIM1JM1K   = 0
            CMIM1JM1KM1 = 0
	    !
            CMIM1JP1K   = 0
            CMIP1JKM1   = 0
            CMIJP1KM1   = 0
            CMIM1JP1KM1 = 0
            CMIP1JP1KM1 = 0
            CMIP1JM1KM1 = 0
            !
            IF( IJK.GT.0.AND.IJK.LT.(NXNYNZ+1)) 
     &           CMIJK = BMG_TO_CM(IJK) 

            IF( IM1JK.GT.0.AND.IM1JK.LT.(NXNYNZ+1)) 
     &           CMIM1JK = BMG_TO_CM(IM1JK)

            IF( IJM1K.GT.0.AND.IJM1K.LT.(NXNYNZ+1)) 
     &           CMIJM1K = BMG_TO_CM(IJM1K)

            IF( IJKM1.GT.0.AND.IJKM1.LT.(NXNYNZ+1)) 
     &           CMIJKM1 = BMG_TO_CM(IJKM1)

            IF( IJM1KM1.GT.0.AND.IJM1KM1.LT.(NXNYNZ+1)) 
     &           CMIJM1KM1 = BMG_TO_CM(IJM1KM1)

            IF( IM1JKM1.GT.0.AND.IM1JKM1.LT.(NXNYNZ+1)) 
     &           CMIM1JKM1   = BMG_TO_CM(IM1JKM1) 

            IF( IM1JM1K.GT.0.AND.IM1JM1K.LT.(NXNYNZ+1)) 
     &            CMIM1JM1K   = BMG_TO_CM(IM1JM1K) 
            
            IF( IM1JM1KM1.GT.0.AND.IM1JM1KM1.LT.(NXNYNZ+1)) 
     &           CMIM1JM1KM1 = BMG_TO_CM(IM1JM1KM1)
 
            IF( IM1JP1K.GT.0.AND.IM1JP1K.LT.(NXNYNZ+1)) 
     &             CMIM1JP1K   = BMG_TO_CM(IM1JP1K)

            IF( IP1JKM1.GT.0.AND.IP1JKM1.LT.(NXNYNZ+1)) 
     &           CMIP1JKM1   = BMG_TO_CM(IP1JKM1)

            IF( IJP1KM1.GT.0.AND.IJP1KM1.LT.(NXNYNZ+1)) 
     &           CMIJP1KM1   = BMG_TO_CM(IJP1KM1)

            IF( IM1JP1KM1.GT.0.AND.IM1JP1KM1.LT.(NXNYNZ+1)) 
     &           CMIM1JP1KM1 = BMG_TO_CM(IM1JP1KM1) 

            IF( IP1JP1KM1.GT.0.AND.IP1JP1KM1.LT.(NXNYNZ+1)) 
     &           CMIP1JP1KM1 = BMG_TO_CM(IP1JP1KM1)
            
            IF( IP1JM1KM1.GT.0.AND.IP1JM1KM1.LT.(NXNYNZ+1)) 
     &           CMIP1JM1KM1 = BMG_TO_CM(IP1JM1KM1) 
            !
            IF( BMG_BDY(KL).NE.0 ) THEN
               !
               DO IT=IA(CMIJK),IA(CMIJK+1)-1
                  !
                  IF( JA(IT).EQ.CMIJK ) THEN
                     !
                     SO(I+1,J+1,K+1,kp) = BMG_BDY(KL)*AA(IT)
                     !
                  ELSEIF( JA(IT).EQ.CMIM1JK .AND. I.NE.1 ) THEN
                     !
                     FACT=ABS(BMG_BDY(IM1JK))*BMG_BDY(KL)
                     SO(I+1,J+1,K+1,kpw) = -FACT*AA(IT)
                     !
                  ELSEIF( JA(IT).EQ.CMIJM1K .AND. J.NE.1 ) THEN
                     !
                     FACT=ABS(BMG_BDY(IJM1K))*BMG_BDY(KL)
                     SO(I+1,J+1,K+1,kps) = -FACT*AA(IT)
                     !
                  ELSEIF( JA(IT).EQ.CMIJKM1 .AND. K.NE.1 ) THEN
                     !
                     FACT=ABS(BMG_BDY(IJKM1))*BMG_BDY(KL)
                     SO(I+1,J+1,K+1,kb) = -FACT*AA(IT)
                     !
                  ELSEIF( JA(IT).EQ.CMIM1JM1K .AND. 
     &                    (I.NE.1 .OR. J.NE.1) ) THEN
                     !
                     FACT=ABS(BMG_BDY(IM1JM1K))*BMG_BDY(KL)
                     SO(I+1,J+1,K+1,kpsw) = -FACT*AA(IT)
                     !
                  ELSEIF( JA(IT).EQ.CMIM1JKM1 .AND. 
     &                    (I.NE.1 .OR. K.NE.1) ) THEN
                     !
                     FACT=ABS(BMG_BDY(IM1JKM1))*BMG_BDY(KL)
                     SO(I+1,J+1,K+1,kbw) = -FACT*AA(IT)
                     !
                  ELSEIF( JA(IT).EQ.CMIJM1KM1 .AND. 
     &                    (J.NE.1 .OR. K.NE.1) ) THEN
                     !
                     FACT=ABS(BMG_BDY(IJM1KM1))*BMG_BDY(KL)
                     SO(I+1,J+1,K+1,kbs) = -FACT*AA(IT)
                     !
                  ELSEIF( JA(IT).EQ.CMIM1JM1KM1 .AND. 
     &                    (I.NE.1 .OR. J.NE.1 .OR. K.NE.1) ) THEN
                     !
                     FACT=ABS(BMG_BDY(IM1JM1KM1))*BMG_BDY(KL)
                     SO(I+1,J+1,K+1,kbsw) = -FACT*AA(IT)
                     !
                  ELSEIF( JA(IT).EQ.CMIM1JP1K .AND. 
     &                    (I.NE.1 .OR. J.NE.NY) )  THEN
                     !
                     FACT=ABS(BMG_BDY(IM1JP1K))*BMG_BDY(KL)
                     SO(I+1,J+2,K+1,kpnw) = -FACT*AA(IT)
                     !
                  ELSEIF( JA(IT).EQ.CMIM1JP1KM1 .AND. 
     &                    (I.NE.1 .OR. J.NE.NY .OR. K.NE.1) ) THEN
                     !
                     FACT=ABS(BMG_BDY(IM1JP1KM1))*BMG_BDY(KL)
                     SO(I+1,J+2,K+1,kbnw) = -FACT*AA(IT)
                     !
                  ELSEIF( JA(IT).EQ.CMIJP1KM1 .AND.
     &                    (J.NE.NY .OR. K.NE.1)  ) THEN
                     !
                     FACT=ABS(BMG_BDY(IJP1KM1))*BMG_BDY(KL)
                     SO(I+1,J+2,K+1,kbn ) = -FACT*AA(IT)
                     !
                  ELSEIF( JA(IT).EQ.CMIP1JP1KM1 .AND. 
     &                    (I.NE.NX .OR. J.NE.NY .OR. K.NE.1) ) THEN
                     !
                     FACT=ABS(BMG_BDY(IP1JP1KM1))*BMG_BDY(KL)
                     SO(I+2,J+2,K+1,kbne) = -FACT*AA(IT)
                     !
                  ELSEIF( JA(IT).EQ.CMIP1JKM1 .AND.
     &                    (I.NE.NX .OR. K.NE.1) ) THEN
                     !
                     FACT=ABS(BMG_BDY(IP1JKM1))*BMG_BDY(KL)
                     SO(I+2,J+1,K+1,kbe) = -FACT*AA(IT)
                     !
                  ELSEIF( JA(IT).EQ.CMIP1JM1KM1 .AND.
     &                    (I.NE.NX .OR. J.NE.1 .OR. K.NE.1) ) THEN
                     !
                     FACT=ABS(BMG_BDY(IP1JM1KM1))*BMG_BDY(KL)
                     SO(I+2,J+1,K+1,kbse) = -FACT*AA(IT)
                     !
                  ENDIF
	          !
               ENDDO
              !
            ELSE
               !
               SO(I+1,J+1,K+1,kp  ) = 1.D0
               SO(I+1,J+1,K+1,kpw ) = 0.D0 ! west
               SO(I+1,J+1,K+1,kps ) = 0.D0 ! south
               SO(I+1,J+1,K+1,kpsw) = 0.D0 ! south-west
               SO(I+1,J+1,K+1,kb  ) = 0.D0 ! below
               SO(I+1,J+1,K+1,kbw ) = 0.D0 ! below-west
               SO(I+1,J+1,K+1,kbs ) = 0.D0 ! below-south
               SO(I+1,J+1,K+1,kbsw) = 0.D0 ! below-south-west
               !
               SO(I+1,J+2,K+1,kpnw) = 0.D0 ! north-west
               SO(I+1,J+2,K+1,kbnw) = 0.D0 ! below-north-west
               SO(I+1,J+2,K+1,kbn ) = 0.D0 ! below-north
               SO(I+1,J+1,K+2,kbne) = 0.D0 ! below-north-east
               SO(I+2,J+1,K+1,kbe ) = 0.D0 ! below-east
               SO(I+2,J+1,K+1,kbse) = 0.D0 ! below-south-east
               !
            ENDIF
            
         ENDDO	
C$OMP END PARALLEL DO

      ELSE


C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(CMIJK,CMIM1JK,CMIJM1K,CMIJKM1,FACT,
C$OMP&         I,IJK,IM1JK,IJM1K,IJKM1,IT,J,K,KL)
C$OMP& SHARED(AA,BMG_BDY,BMG_TO_CM,IA,JA,NX,NY,NXNY,NXNYNZ,SO)
         DO KL = 1,NXNYNZ
	    !
            I = MOD((KL-1),NX)+1
            J = MOD((KL-1)/NX,NY)+1
            K = (KL-1)/NXNY+1
	    !
            IJK       = (K-1)*NXNY+(J-1)*NX+I
            IM1JK     = IJK - 1
            IJM1K     = IJK - NX
            IJKM1     = IJK - NXNY
            !
            CMIJK     = BMG_TO_CM(IJK)
            IF( I.NE.1 ) THEN
               CMIM1JK = BMG_TO_CM(IM1JK)
            ELSE
               CMIM1JK = -1
            ENDIF
            !
            IF( J.NE.1 ) THEN
               CMIJM1K = BMG_TO_CM(IJM1K)
            ELSE
               CMIJM1K = -1
            ENDIF
            !
            IF( K.NE.1 ) THEN
               CMIJKM1 = BMG_TO_CM(IJKM1)
            ELSE
               CMIJKM1 = -1
            ENDIF
            !
            IF( BMG_BDY(KL).NE.0 ) THEN
               !
               DO IT=IA(CMIJK),IA(CMIJK+1)-1
                  !
                  IF( JA(IT).EQ.CMIJK ) THEN
                     !
                     SO(I+1,J+1,K+1,kp) = BMG_BDY(KL)*AA(IT)
                     !
                  ELSEIF( JA(IT).EQ.CMIM1JK .AND. I.NE.1 ) THEN
                     !
                     FACT=ABS(BMG_BDY(IM1JK))*BMG_BDY(KL)
                     SO(I+1,J+1,K+1,kpw) = -FACT*AA(IT)
                     !
                  ELSEIF( JA(IT).EQ.CMIJM1K .AND. J.NE.1 ) THEN
                     !
                     FACT=ABS(BMG_BDY(IJM1K))*BMG_BDY(KL)
                     SO(I+1,J+1,K+1,kps) = -FACT*AA(IT)
                     !
                  ELSEIF( JA(IT).EQ.CMIJKM1 .AND. K.NE.1 ) THEN
                     !
                     FACT=ABS(BMG_BDY(IJKM1))*BMG_BDY(KL)
                     SO(I+1,J+1,K+1,kb) = -FACT*AA(IT)
                     !
                  ENDIF
	          !
               ENDDO
              !
            ELSE
               !
               SO(I+1,J+1,K+1,kp  ) = 1.D0
               SO(I+1,J+1,K+1,kpw ) = 0.D0 ! west
               SO(I+1,J+1,K+1,kps ) = 0.D0 ! south
               SO(I+1,J+1,K+1,kb  ) = 0.D0 ! below
               !
            ENDIF
            !
         ENDDO
         !
      ENDIF
      !
      ! Update ghost boundaries if periodic problem
      !
      BC = BMG_iPARMS(id_BMG3_BC)
      !
      IF( BC.EQ.BMG_BCs_def_per_x .OR. BC.EQ.BMG_BCs_indef_per_x ) THEN


      ELSEIF( BC.EQ.BMG_BCs_def_per_y .OR. BC.EQ.BMG_BCs_indef_per_y ) 
     &        THEN


      ELSEIF( BC.EQ.BMG_BCs_def_per_z .OR. BC.EQ.BMG_BCs_indef_per_z ) 
     &        THEN

         


      ENDIF


      RETURN
      END
