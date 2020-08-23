      SUBROUTINE BMG2_SymStd_PRECON_diag( 
     &                SO, QF, Q, II, JJ, NStncl
     &                )

C ==========================================================================

C***BEGIN PROLOGUE BMG2_SymStd_PRECON_diag
C***SUBSIDIARY
C***PURPOSE  BMG_SymStd_PRECON_diag solves the system D*Q = QF
C            where D is the diagonal of the matrix in SO.
C
C***LIBRARY   SLATEC
C***AUTHOR  DENDY, J. E. JR.
C             LOS ALAMOS NATIONAL LABORATORY
C           VOYTKO, M. H.
C             LAWRENCE LIVERMORE NATIONAL LABORATORY
C
C          Copyright, 1988. The Regents of the University of California.
C          This software was produced under a U.S. Government contract
C          (W-7405-ENG-36) by Los Alamso National Laboratory, which is
C          operated by the University of California for the U.S. Depart-
C          ment of Energy. The U.S. Government is licensed to use,
C          reproduce, and distribute this software. Permission is
C          granted to the public to copy and use this software without
C          charge, provided that this Notice and any statement of
C          authorship are reproduced on all copies. Neither the
C          Government nor the University makes any warranty, expressed
C          or implied, or assumes any liability or responsibility for
C          the use of this software.
C
C***PARAMETERS
C***INPUT
C
C   SOR       Refer to BOXMG.
C   QF        Refer to BOXMG.
C   II        Number of grid points in x direction, including
C             two fictitious points.
C   JJ        Number of grid points in y direction, including
C             two fictitious points.
C
C***OUTPUT
C
C   Q         Refer to BOXMG.
C
C
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   030710  DATE WRITTEN (TMA)
C
C***END PROLOGUE  BMG2_SymStd_PRECON_diag

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
      INTEGER II, JJ, NStncl

      REAL*8  Q(II,JJ), QF(II,JJ), SO(II,JJ,NStncl)

C ----------------------------
C     Local Declarations
C
      INTEGER I, I1, J, J1, ERR
      REAL*8 tiny

      PARAMETER ( tiny = 1e-10 )

C =========================================================================

      J1=JJ-1
      I1=II-1

C -------------------------------------------------------------------

C     ----------------------------------------
C     PRE-CHECK DIAGONAL FOR ZEROS
C     ----------------------------------------

      DO J=2,J1
         DO I=2,I1
            IF ( ABS(SO(I,J,KO)).LT.tiny ) THEN  
               WRITE(*,*)'Error in BMG2_SymStd_PRECON_diag:  '
               WRITE(*,*)'divide by zero.  Set D = I'
               ERR = 1
            ENDIF
         ENDDO
      ENDDO

      IF(ERR.EQ.0) THEN

C       ----------------------------------------
C       SOLVE D*Q = QF where D is diagonal of A.
C       ----------------------------------------

         DO J=2,J1
            DO I=2,I1
               Q(I,J) = QF(I,J)/SO(I,J,KO)  
            ENDDO
	 ENDDO

      ELSE
         
C       ----------------------------------------
C       SET Q = QF by calling dcopy routine.
C       ----------------------------------------
        
         CALL BMG2_SymStd_UTILS_dcopy( 
     &             QF, Q, II, JJ 
     &             )

      ENDIF


C =========================================================================

      RETURN
      END

