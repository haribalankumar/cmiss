      SUBROUTINE BMG2_SymStd_relax_lines_x ( 
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

      INTEGER IRELAX_SYM, IFD, K, KF, UPDOWN
      REAL*8  B(II), Q(II,JJ), QF(II,JJ), SO(II,JJ,NStncl), SOR(II,JJ,2)

C ----------------------------
C     Local Declarations
C
      INTEGER I, I1, J, J1, JBEG, JEND
      INTEGER JBEG_START, JBEG_END, JBEG_STRIDE
      INTEGER INFO

C =========================================================================

      J1=JJ-1
      I1=II-1
C
C     on the way down we relax red lines and then black lines
C     on the way up we relax in the opposite order
C

      IF (IRELAX_SYM.EQ.BMG_RELAX_NONSYM.OR.(UPDOWN.EQ.BMG_DOWN 
     &     .AND.IRELAX_SYM.EQ.BMG_RELAX_SYM)) THEN
         !
         ! Relax red lines, then black lines.
         !     
         JBEG_START = 3
         JBEG_END   = 2
         JBEG_STRIDE= -1
      ELSE
         !
         ! Relax black lines, then red lines.
         !     
         JBEG_START = 2
         JBEG_END   = 3
         JBEG_STRIDE= 1
      ENDIF

      IF ( K.LT.KF .OR. IFD.NE.1 ) THEN
         !
         !  9-point stencil
         !
         DO JBEG=JBEG_START, JBEG_END, JBEG_STRIDE
            JEND=2*((J1-JBEG)/2)+JBEG
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& FIRSTPRIVATE(I1,JBEG,JEND)
C$OMP& PRIVATE(I,J)
C$OMP& SHARED(INFO,Q,QF,SO,SOR)
            DO J=JBEG,JEND,2
               DO I=2,I1
                  Q(I,J)=QF(I,J)+SO(I,J,KS)*Q(I,J-1)+SO(I,J+1,KS)
     &                 *Q(I,J+1)+SO(I,J,KSW)*Q(I-1,J-1)+SO(I+1,J,KNW)
     &                 *Q(I+1,J-1)+SO(I,J+1,KNW)*Q(I-1,J+1)
     &                 +SO(I+1,J+1,KSW)*Q(I+1,J+1)
               ENDDO
               CALL DPTTRS (I1-1, 1, SOR(2,J,1), SOR(3,J,2), 
     &              Q(2,J), I1-1, INFO)
            ENDDO
C$OMP END PARALLEL DO

         ENDDO
         !
      ELSE
         !
         !  5-point stencil
         !
         DO JBEG=JBEG_START, JBEG_END, JBEG_STRIDE
            JEND=2*((J1-JBEG)/2)+JBEG
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(I,J)
C$OMP& SHARED(I1,INFO,JBEG,JEND,Q,QF,SO,SOR)
            DO J=JBEG,JEND,2
               DO I=2,I1
                  Q(I,J)=QF(I,J)+SO(I,J,KS)*Q(I,J-1)+SO(I,J+1,KS)
     &                 *Q(I,J+1)
               ENDDO

               CALL DPTTRS (I1-1, 1, SOR(2,J,1), SOR(3,J,2), 
     &              Q(2,J), I1-1, INFO)
            ENDDO
C$OMP END PARALLEL DO

         ENDDO
         !
      ENDIF

      RETURN
      END
