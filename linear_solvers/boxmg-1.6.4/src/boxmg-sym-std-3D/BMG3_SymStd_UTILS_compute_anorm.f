      SUBROUTINE BMG3_SymStd_UTILS_compute_anorm(
     &                KG, SO, II, JJ, KK, 
     &                KF, IFD, NStncl, ANORM
     &                )

C ==========================================================================
C
C   BMG3_SymStd_UTILS_compute_anorm
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_UTILS_compute_anorm.f is self-explanatory.
C
C   ---------------------
C   HISTORY:
C   ---------------------
C
C   Written:    2004/12/22 (TMA)
C
C
C ==================================================================
C   INPUT:
C ========================
C
C
C
C ==================================================================
C   OUTPUT:
C ===========================
C
C
C
C ==================================================================
C   LOCAL:
C ========================
C
C
C
C ==========================================================================

      IMPLICIT NONE

C -----------------------------
C     Includes
C
      INCLUDE 'BMG_constants.h'
      INCLUDE 'BMG_stencils.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER II, JJ, KK, KG, KF, NStncl, IFD

      REAL*8  SO(II,JJ,KK,NStncl), ANORM

C ----------------------------
C     Local Declarations
C
      INTEGER I, I1, J, J1, K, K1
      REAL*8  SUM

C =========================================================================

      K1=KK-1
      J1=JJ-1
      I1=II-1

C     -------------------------------------------------------------------

      IF ( KG.LT.KF .OR. IFD.NE.1 ) THEN
         !
         !  27-point stencil
         ! 
         ANORM = 0.D0
         DO K=3,K1-1
            DO J=3,J1-1
               DO I=3,I1-1
                  SUM = ABS(SO(I,J,K,kp)) 
     &              + ABS(SO(i,j,k,kpw))     + ABS(SO(i+1,j,k,kpw))
     &              + ABS(SO(i,j,k,kps))     + ABS(SO(i,j+1,k,kps))
     &              + ABS(SO(i,j,k,kpsw))    + ABS(SO(i+1,j+1,k,kpsw))
     &              + ABS(SO(i,j, k ,kb))    + ABS(SO(i,j,k+1,kb))
     &              + ABS(SO(i,j,k,kbs))     + ABS(SO(i,j+1,k+1,kbs))
     &              + ABS(SO(i,j,k,kbsw))    + ABS(SO(i+1,j+1,k+1,kbsw))
     &              + ABS(SO(i,j+1,k,kpnw))  + ABS(SO(i+1,j,k,kpnw))
     &              + ABS(SO(i,j+1,k,kbnw))  + ABS(SO(i+1,j,k+1,kbnw))
     &              + ABS(SO(i,j+1,k,kbn))   + ABS(SO(i, j ,k+1,kbn))
     &              + ABS(SO(i+1,j+1,k,kbne))+ ABS(SO(i,j,k+1,kbne))
     &              + ABS(SO(i+1,j,k,kbe))   + ABS(SO(i,j,k+1,kbe))
     &              + ABS(SO(i+1,j,k,kbse))  + ABS(SO(i,j+1,k+1,kbse))
     &              + ABS(SO(i+1,j,k+1,kbw)) + ABS(SO(i,j,k,kbw))
                  IF( SUM.GT.ANORM ) ANORM = SUM
               ENDDO
            ENDDO
         ENDDO
         !
      ELSE
         !
         !  7-point stencil
         ! 
         ANORM = 0.D0
         DO K=2,K1       
            DO J=2,J1
               DO I=2,I1
                  SUM  = ABS(SO(I,J,K,kp)) 
     &                 + ABS(SO(I,J,K,kpw)) + ABS(SO(I+1,J,K,kpw)) 
     &                 + ABS(SO(I,J,K,kps)) + ABS(SO(I,J+1,K,kps))
     &                 + ABS(SO(I,J,K,kb )) + ABS(SO(I,J,K+1,kb ))
                  IF( SUM.GT.ANORM ) ANORM = SUM
               ENDDO
            ENDDO
         ENDDO
         !
      ENDIF

C =========================================================================

      RETURN
      END

