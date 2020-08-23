      SUBROUTINE BMG3_SymStd_PRECON_diag( 
     &                SO, QF, Q, II, JJ, KK, NStncl
     &                )

C ==========================================================================

C***BEGIN PROLOGUE BMG3_SymStd_PRECON_diag
C***SUBSIDIARY
C***PURPOSE  BMG3_SymStd_PRECON_diag solves the system D*Q = QF
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
C***END PROLOGUE  BMG3_SymStd_PRECON_diag

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
      INTEGER II, JJ, KK, NStncl

      REAL*8  Q(II,JJ,KK), QF(II,JJ,KK), SO(II,JJ,KK,NStncl)

C ----------------------------
C     Local Declarations
C
      INTEGER I, I1, J, J1, K, K1, KL, ERR
      REAL*8 tiny

      PARAMETER ( tiny = 1e-10 )

C =========================================================================

      K1=KK-1
      J1=JJ-1
      I1=II-1

C -------------------------------------------------------------------

C     ----------------------------------------
C     PRE-CHECK DIAGONAL FOR ZEROS
C     ----------------------------------------

      DO K=2,K1
         DO J=2,J1
            DO I=2,I1
               IF ( ABS(SO(I,J,K,KO)).LT.tiny ) THEN  
                  WRITE(*,*)'Error in BMG3_SymStd_PRECON_diag:  '
                  WRITE(*,*)'divide by zero.  Set D = I'
                  ERR = 1
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      IF(ERR.EQ.0) THEN

C       ----------------------------------------
C       SOLVE D*Q = QF where D is diagonal of A.
C       ----------------------------------------

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(I,J,K,KL)
C$OMP& SHARED(I1,J1,K1,so,qf,q)
         DO kl=1,(I1-2)*(J1-2)*(K1-2)

            I = MOD((kl-1),(I1-2))+2
            J = MOD((kl-1)/(I1-2),(J1-2))+2
            K = (kl-1)/((I1-2)*(J1-2))+2

            Q(I,J,K) = QF(I,J,K)/SO(I,J,K,KO)  
            
         ENDDO
C$OMP END PARALLEL DO

      ELSE
         
C       ----------------------------------------
C       SET Q = QF by calling dcopy routine.
C       ----------------------------------------
        
         CALL BMG3_SymStd_UTILS_dcopy(
     &             QF, Q, II, JJ, KK 
     &             )

      ENDIF


C =========================================================================

      RETURN
      END

