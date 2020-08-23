      SUBROUTINE BMG2_SymStd_UTILS_compute_anorm(
     &                K, SO, II, JJ, 
     &                KF, IFD, NStncl, ANORM
     &                )

C ==========================================================================
C
C   BMG2_SymStd_UTILS_compute_anorm
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_SymStd_UTILS_compute_anorm.f is self-explanatory.
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
      INTEGER II, JJ, K, KF, NStncl, IFD

      REAL*8  SO(II,JJ,NStncl), ANORM

C ----------------------------
C     Local Declarations
C
      INTEGER I, I1, J, J1
      REAL*8  SUM

C =========================================================================

      J1=JJ-1
      I1=II-1

C     -------------------------------------------------------------------

      IF ( K.LT.KF .OR. IFD.NE.1 ) THEN
         !
         !  9-point stencil
         ! 
         ANORM = 10000.D0
         DO J=2,J1
            DO I=2,I1
               IF( ABS(SO(I,J,KSW)).LT.ANORM .AND. 
     &             ABS(SO(I,J,KSW)).GT.0.0D0       ) THEN
                  ANORM = ABS(SO(I,J,KSW))
               ENDIF
             ENDDO
         ENDDO
         !
      ELSE
         !
         !  5-point stencil
         !        
         ANORM = 0.0D0
         DO J=2,J1
            DO I=2,I1
               IF( ABS(SO(I,J,KSW)).GT.ANORM ) THEN
                  ANORM = ABS(SO(I,J,KS))
               ENDIF
            ENDDO
         ENDDO
         !
      ENDIF

C =========================================================================

      RETURN
      END

