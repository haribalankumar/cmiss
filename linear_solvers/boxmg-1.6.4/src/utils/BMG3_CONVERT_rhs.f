      SUBROUTINE BMG3_CONVERT_rhs( 
     &              QF, BMG_BC, BMG_BDY, CM_IJK_MAP, BMG_TO_CM,
     &              AA, IA, JA, B, NX, NY, NZ 
     &              )

C =======================================================================
C 
C   BMG3_CONVERT_rhs 
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_CONVERT_rhs converts a standard 2D RHS array into a BOXMG 
C   (padded) 2D RHS array. 
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
C   CMISS RHS Array:
C   -----------------------
C
C   B        Standard CMISS right-hand side array.
C
C =======================================================================
C   OUTPUT:
C ========================
C
C   -----------------------
C   BOXMG Solution Array:
C   -----------------------
C
C   QF        2D BOXMG right-hand side array.
C
C ==========================================================================

      IMPLICIT NONE
      INCLUDE 'BMG_parameters.h'

C ========================================================================== 

C ---------------- 
C     Arguments  
C
      INTEGER  BMG_BC, BMG_BDY(*), BMG_TO_CM(*), CM_IJK_MAP(*), IA(*), 
     &         JA(*), NX, NY, NZ
      REAL*8   AA(*), B(NX*NY*NZ), QF(NX+2,NY+2,NZ+2)

C -------------------- 
C     Local Variables  
C
      INTEGER CMIP, CMJP, CMKL, CMKP, CMKLP, I, IP, J, JP, 
     &        K, KP, KL, KLP, KT, NXNY, NXNYNZ, MINIP, MAXIP, 
     &        MINJP, MAXJP, MINKP, MAXKP

C ========================================================================== 

      NXNYNZ = NX*NY*NZ
      NXNY   = NX*NY

      IF( BMG_BC.NE.BMG_BCs_indef_nonper ) THEN

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(CMKL,I,J,K,KL)
C$OMP& SHARED(B,BMG_BDY,BMG_TO_CM,CM_IJK_MAP,NXNYNZ,QF)
         DO KL = 1,NXNYNZ
            !
            I = CM_IJK_MAP(3*KL-2)
            J = CM_IJK_MAP(3*KL-1)
            K = CM_IJK_MAP(3*KL  )
            ! 
            CMKL = BMG_TO_CM(KL) 
            !
            IF( BMG_BDY(KL).NE.0 ) THEN
               QF(I+1,J+1,K+1) = BMG_BDY(KL)*B(CMKL)
            ELSE
               QF(I+1,J+1,K+1) = B(CMKL)
            ENDIF
            !
         ENDDO
C$OMP END PARALLEL DO

      ELSE

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(CMKL,I,J,K,KL)
C$OMP& SHARED(B,BMG_BDY,BMG_TO_CM,CM_IJK_MAP,NXNYNZ,QF)
         DO KL = 1,NXNYNZ
            !
            I = CM_IJK_MAP(3*KL-2)
            J = CM_IJK_MAP(3*KL-1)
            K = CM_IJK_MAP(3*KL  )
            ! 
            CMKL = BMG_TO_CM(KL) 
            !
            QF(I+1,J+1,K+1) = BMG_BDY(KL)*B(CMKL)
            !
         ENDDO
C$OMP END PARALLEL DO
      ENDIF         

      IF( BMG_BC.NE.BMG_BCs_indef_nonper ) THEN

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(CMIP,CMJP,CMKP,CMKL,CMKLP,I,IP,J,JP,K,KP,KL,KT,
C$OMP&         KLP,MAXIP,MAXJP,MAXKP,MINIP,MINJP,MINKP)
C$OMP& SHARED(AA,B,BMG_BDY,BMG_TO_CM,CM_IJK_MAP,
C$OMP&        IA,JA,NXNYNZ,NX,NY,NZ,NXNY,QF)
         DO KL = 1, NXNYNZ
            !
            I = CM_IJK_MAP(3*KL-2)
            J = CM_IJK_MAP(3*KL-1)
            K = CM_IJK_MAP(3*KL  )
            !
            CMKL = BMG_TO_CM(KL) ! Get CM Numbering
            !
            IF( BMG_BDY(KL).EQ.0 ) THEN
               !
               MINIP = -1
               MAXIP = -1
               !
               MINJP = -1
               MAXJP = 1
               !
               MINKP = -1
               MAXKP = 1
               !
               IF( I.EQ.1) THEN
                  MINIP = 0
                  MAXIP = 1
               ENDIF
               !
               IF( I.EQ.NX) THEN
                  MINIP = -1
                  MAXIP = 0
               ENDIF
               !
               IF( J.EQ.1 ) THEN
                  MINJP = 0
                  MAXJP = 1
               ENDIF
               !
               IF( J.EQ.NY ) THEN
                  MINJP = -1
                  MAXJP = 0
               ENDIF
               !
               IF( K.EQ.1 ) THEN
                  MINKP = 0
                  MAXKP = 1
               ENDIF
               !
               IF( K.EQ.NZ ) THEN
                  MINKP = -1
                  MAXKP = 0
               ENDIF
               !
               DO KP=MINKP,MAXKP,1
               DO JP=MINJP,MAXJP,1
               DO IP=MINIP,MAXIP,1
                  !
                  IF( .NOT.(IP.EQ.0 .AND. JP.EQ.0 .AND. KP.EQ.0) ) THEN
                     !
                     KLP = KL + KP*NXNY + JP*NX + IP
                     !
                     CMKLP = BMG_TO_CM(KLP)
                     !
                     CMIP = CM_IJK_MAP(3*KLP-2)
                     CMJP = CM_IJK_MAP(3*KLP-1)
                     CMKP = CM_IJK_MAP(3*KLP  )
                     !
                     IF( BMG_BDY(KLP).NE.0 ) THEN
                        !
                        DO KT=IA(CMKLP),IA(CMKLP+1)-1
                           !
                           IF( JA(KT).EQ.CMKL ) THEN
                              QF(CMIP+1,CMJP+1,CMKP+1) = 
     &                         QF(CMIP+1,CMJP+1,CMKP+1) - AA(KT)*B(CMKL)
                           ENDIF
                           !
                        ENDDO
                        !
                     ENDIF
                     !
                  ENDIF
                  !
               ENDDO
               ENDDO
               ENDDO
               !
            ENDIF
            !
         ENDDO
C$OMP END PARALLEL DO

      ENDIF

      IF( BMG_BC.EQ.BMG_BCs_indef_nonper ) THEN
         CALL BMG3_SymStd_UTILS_zero_mean(QF,NX,NY,NZ)
      ENDIF

      RETURN
      END
