      SUBROUTINE BMG2_CONVERT_rhs( 
     &              QF, BMG_BC, BMG_BDY, CM_IJK_MAP, BMG_TO_CM,
     &              AA, IA, JA, B, NX, NY 
     &              )

C =======================================================================
C 
C   BMG2_CONVERT_rhs 
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_CONVERT_rhs converts a standard 2D RHS array into a BOXMG 
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
      INTEGER  BMG_BC, BMG_BDY(*), CM_IJK_MAP(*), BMG_TO_CM(*), IA(*), 
     &         JA(*), NX, NY
      REAL*8   AA(*), B(NX*NY), QF(NX+2,NY+2)

C -------------------- 
C     Local Variables
C
      INTEGER CMIP, CMJP, CMKL, CMKLP,  
     &        I, IP, J, JP, KL, KLP, KT, NXNY,
     &        MINIP, MAXIP, MINJP, MAXJP

C ========================================================================== 

      NXNY = NX*NY

      IF( BMG_BC.NE.BMG_BCs_indef_nonper ) THEN

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(CMKL,I,J,KL)
C$OMP& SHARED(B,BMG_BDY,BMG_TO_CM,CM_IJK_MAP,NXNY,QF)
         DO KL = 1, NXNY
            !
            I = CM_IJK_MAP(2*KL-1)
            J = CM_IJK_MAP(2*KL  )
            ! 
            CMKL = BMG_TO_CM(KL)  ! Get CM Numbering
            !
            !PRINT *, CMKL, KL, I, J
            IF( BMG_BDY(KL).NE.0 ) THEN
               QF(I+1,J+1) = BMG_BDY(KL)*B(CMKL)
            ELSE
               QF(I+1,J+1) = B(CMKL)
            ENDIF
            !
         ENDDO
C$OMP END PARALLEL DO

      ELSE
         
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(CMKL,I,J,KL)
C$OMP& SHARED(B,BMG_BDY,BMG_TO_CM,CM_IJK_MAP,NXNY,QF)
         DO KL = 1, NXNY
            !
            I = CM_IJK_MAP(2*KL-1)
            J = CM_IJK_MAP(2*KL  )
            !
            CMKL = BMG_TO_CM(KL)  ! Get CM Numbering
            !
            QF(I+1,J+1) = BMG_BDY(KL)*B(CMKL)
            !
         ENDDO
C$OMP END PARALLEL DO

      ENDIF

      IF( BMG_BC.NE.BMG_BCs_indef_nonper ) THEN
C$OMP PARALLEL DO 
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(CMIP,CMJP,CMKL,CMKLP,I,IP,J,JP,KL,KLP,KT,
C$OMP&         MAXIP,MAXJP,MINIP,MINJP)
C$OMP& SHARED(AA,B,BMG_BDY,BMG_TO_CM,CM_IJK_MAP,IA,JA,NX,NY,
C$OMP&        NXNY,QF)
         DO KL = 1, NXNY                                    
            !
            I = CM_IJK_MAP(2*KL-1)
            J = CM_IJK_MAP(2*KL  )
            !
            CMKL = BMG_TO_CM(KL) 
            !
            IF( BMG_BDY(KL).EQ.0 ) THEN
              !
               MINIP = -1
               MAXIP = -1
               !
               MINJP = -1
               MAXJP = 1
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
               IF( J.EQ.NY) THEN
                  MINJP = -1
                  MAXJP = 0
               ENDIF
               !
               DO JP=MINJP,MAXJP,1
               DO IP=MINIP,MAXIP,1
                  !
                  IF( .NOT.(IP.EQ.0 .AND. JP.EQ.0) ) THEN
                     !
                     KLP = KL + JP*NX + IP
                     !
                     CMKLP = BMG_TO_CM(KLP)
                     !
                     CMIP = CM_IJK_MAP(2*KLP-1)
                     CMJP = CM_IJK_MAP(2*KLP  )
                     !
                     IF( BMG_BDY(KLP).NE.0 ) THEN
                        !
                        DO KT=IA(CMKLP),IA(CMKLP+1)-1
                           !
                           IF( JA(KT).EQ.CMKL ) THEN
                              QF(CMIP+1,CMJP+1) = 
     &                             QF(CMIP+1,CMJP+1) - AA(KT)*B(CMKL)
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
               !
            ENDIF
            !
         ENDDO
C$OMP END PARALLEL DO

      ENDIF
            
      RETURN
      END
