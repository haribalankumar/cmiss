      SUBROUTINE BMG2_SymStd_relax_GS ( 
     &                K, SO, QF, Q, SOR, II, JJ, 
     &                KF, IFD, NStncl, NSORv, IRELAX_SYM, UPDOWN 
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
      INTEGER II, JJ, NStncl, NSORv

      INTEGER IFD, IRELAX_SYM, K, KF, UPDOWN
      REAL*8  Q(II,JJ), QF(II,JJ), SO(II,JJ,NStncl), SOR(II,JJ,NSORv)

C ----------------------------
C     Local Declarations
C
      INTEGER I, I1, I2, IBEG, IEND, J, J1, J2, JBEG, JEND, JO
      INTEGER LSTART, LEND, LSTRIDE

C =========================================================================

      J1=JJ-1
      I1=II-1
      J2=JJ-2
      I2=II-2

      IF (IRELAX_SYM.EQ.BMG_RELAX_NONSYM.OR.(UPDOWN.EQ.BMG_DOWN 
     &     .AND.IRELAX_SYM.EQ.BMG_RELAX_SYM)) THEN
         LSTART = 2
         LEND   = 3
         LSTRIDE= 1
      ELSE
         LSTART = 3
         LEND   = 2
         LSTRIDE=-1
      ENDIF
      
      IF ( K.LT.KF .OR. IFD.NE.1 ) THEN
         !
         ! 9-point stencil
         !
         DO JBEG=LSTART,LEND,LSTRIDE
            JEND=2*((J1-JBEG)/2)+JBEG
C$OMP PARALLEL DO 
CC$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(I,IBEG,IEND,J)
C$OMP& SHARED(I1,JBEG,JEND,JO,J1,
C$OMP&        LSTART,LEND,LSTRIDE,
C$OMP&        Q,QF,SO,SOR)
            DO  J=JBEG,JEND,2
               DO  IBEG=LSTART,LEND,LSTRIDE
                  IEND=2*((I1-IBEG)/2)+IBEG
                  DO  I=IBEG,IEND,2
                     
                     Q(I,J) = ( QF(I,J) 
     &                         + SO(I,J,KW)*Q(I-1,J)
     &                         + SO(I+1,J,KW)*Q(I+1,J)
     &                         + SO(I,J,KS)*Q(I,J-1)
     &                         + SO(I,J+1,KS)*Q(I,J+1)
     &                         + SO(I,J,KSW)*Q(I-1,J-1)
     &                         + SO(I+1,J,KNW)*Q(I+1,J-1)
     &                         + SO(I,J+1,KNW)*Q(I-1,J+1)
     &                         + SO(I+1,J+1,KSW)*Q(I+1,J+1)
     &                        )*SOR(I,J,MSOR)
                  ENDDO
               ENDDO
            ENDDO
C$OMP END PARALLEL DO
         ENDDO
         !
      ELSE
         !
         ! 5-point stencil
         !
         DO JO=LSTART,LEND,LSTRIDE
C$OMP PARALLEL DO 
CC$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(I,J,IBEG,IEND)
C$OMP& SHARED(I1,JO,J1,Q,QF,SO,SOR)
            DO J=2,J1
               IBEG=MOD(J+JO,2)+2
               IEND=2*((I1-IBEG)/2)+IBEG
               DO I=IBEG,IEND,2
                  Q(I,J) = ( QF(I,J) 
     &                      + SO(I,J,KW)*Q(I-1,J)
     &                      + SO(I+1,J,KW)*Q(I+1,J)
     &                      + SO(I,J,KS)*Q(I,J-1)
     &                      + SO(I,J+1,KS)*Q(I,J+1)
     &                     )*SOR(I,J,MSOR)                 
                  SOR(I,J,MTOT)=RZERO

               ENDDO
            ENDDO
C$OMP END PARALLEL DO
         ENDDO
         !
      ENDIF

      RETURN
      END
