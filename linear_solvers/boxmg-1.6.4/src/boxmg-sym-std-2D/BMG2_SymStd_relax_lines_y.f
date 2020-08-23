      SUBROUTINE BMG2_SymStd_relax_lines_y( 
     &                K, SO, QF, Q, SOR, B, II, JJ, 
     &                KF, IFD, NStncl, IRELAX_SYM, UPDOWN 
     &                )

C ==========================================================================

      IMPLICIT NONE

C -----------------------------
C     Includes
C
      INCLUDE 'BMG_constants.h'
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_parameters.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER II, JJ, NStncl

      INTEGER IFD, IRELAX_SYM, K, KF, UPDOWN
      REAL*8  B(JJ), Q(II,JJ), QF(II,JJ), SO(II,JJ,NStncl), SOR(JJ,II,2)

C ----------------------------
C     Local Declarations
C
      INTEGER I, I1, IBEG, IEND, J, J1
      INTEGER IBEG_START, IBEG_END, IBEG_STRIDE
      INTEGER INFO

C =========================================================================

      J1=JJ-1
      I1=II-1
      
      IF (IRELAX_SYM.EQ.BMG_RELAX_NONSYM.OR.(UPDOWN.EQ.BMG_DOWN 
     &     .AND.IRELAX_SYM.EQ.BMG_RELAX_SYM)) THEN
         !
         ! Relax red lines, then black lines.
         !     
         IBEG_START = 3
         IBEG_END   = 2
         IBEG_STRIDE= -1
      ELSE
         !
         ! Relax black lines, then red lines.
         !     
         IBEG_START = 2
         IBEG_END   = 3
         IBEG_STRIDE= 1
      ENDIF

      
      IF ( K.LT.KF .OR. IFD.NE.1 ) THEN
         !     
         !  9 pt. operator
         !  Relax red lines, then black lines.
         !
         DO IBEG=IBEG_START,IBEG_END,IBEG_STRIDE
            IEND=2*((I1-IBEG)/2)+IBEG

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(B,I,J)
C$OMP& SHARED(IBEG,IEND,J1,INFO,Q,QF,SO,SOR)
            DO  I=IBEG,IEND,2
               DO  J=2,J1
                  B(J)= QF(I,J)+SO(I,J,KW)*Q(I-1,J)+SO(I+1,J,KW)
     &                 *Q(I+1,J)+SO(I,J,KSW)*Q(I-1,J-1)+SO(I+1,J,KNW)
     &                 *Q(I+1,J-1)+SO(I,J+1,KNW)*Q(I-1,J+1)
     &                 +SO(I+1,J+1,KSW)*Q(I+1,J+1)
               ENDDO

               CALL DPTTRS (J1-1, 1, SOR(2,I,1), SOR(3,I,2), 
     &              B(2), J1-1, INFO)

               DO j=2,J1
                  Q(I,J) = B(J)
               ENDDO
            ENDDO
C$OMP END PARALLEL DO

         ENDDO
         !
      ELSE
         !     
         ! 5 pt. operator
         ! Relax red lines, then black lines.
         !
         DO IBEG=IBEG_START,IBEG_END,IBEG_STRIDE
            IEND=2*((I1-IBEG)/2)+IBEG

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& FIRSTPRIVATE(J1,IBEG,IEND)
C$OMP& PRIVATE(B,I,J)
C$OMP& SHARED(INFO,Q,QF,SO,SOR)
            DO  I=IBEG,IEND,2
               DO  J=2,J1
                  B(J)=QF(I,J)+SO(I,J,KW)*Q(I-1,J)+SO(I+1,J,KW)
     &                 *Q(I+1,J)
               ENDDO

               CALL DPTTRS (J1-1, 1, SOR(2,I,1), SOR(3,I,2), 
     &              B(2), J1-1, INFO)
               DO j=2,J1
                  Q(I,J) = B(J)
               ENDDO
            ENDDO
C$OMP END PARALLEL DO
            
         ENDDO
	   !
      ENDIF
      
      RETURN
      END
