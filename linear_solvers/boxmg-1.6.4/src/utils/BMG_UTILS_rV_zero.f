
      SUBROUTINE BMG_UTILS_rV_zero( Q, N )

C 
C   BMG_UTILS_rV_zero
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG_UTILS_rV_zero zeros out an arbitrary one-dimensional real array.
C
C   Written:    2005/01/22 (TMA)
C
C =======================================================================
C   INPUT:
C ========================
C
C   N    Length of array.
C
C =======================================================================
C   OUTPUT:
C ========================
C
C   Q    Array that is zeroed out.
C
C =======================================================================

      IMPLICIT NONE

C ---------------------------
C     Argument Declarations:
C
      INTEGER N
      REAL*8  Q(*)

C ---------------------------
C     Local Arguments:
C
      INTEGER I

C ========================================================================


C$OMP PARALLEL DO 
      DO I = 1, N
         Q(I)= 0.D0
      ENDDO
C$OMP END PARALLEL DO

      RETURN
      END
