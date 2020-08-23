      SUBROUTINE BMG2_SymStd_SOLVE_cg(
     &                       Q, QF, II, JJ, ABD, BBD, NABD1, NABD2,
     &                       BMG_IOFLAG, BMG_iPARMS
     &                       )

C ==========================================================================
C
C***BEGIN PROLOGUE  MGSADD
C***SUBSIDIARY
C***PURPOSE  mgsadd does a direct solve on the coarsest
C            grid. it uses the linpack routine spbsl.
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
C   QF        Refer to BOXMG.
     
C   II        Number of grid points in x direction, including
C             two fictitious points.
C   JJ        Number of grid points in y direction, including
C             two fictitious points.
C   ABD       Refer to BOXMG.
C
C   NABD1     Refer to BOXMG.
C***INPUT/OUTPUT
C   ABD       Refer to BOXMG.
C   BBD       Refer to BOXMG.
C   NABD1     Refer to BOXMG.
C***OUTPUT
C   Q         Refer to BOXMG
C
C***ROUTINES CALLED  SPBSL
C***REVISION HISTORY  (YYMMDD)
C   830925  DATE WRITTEN
C   900627  modified to conform to the 4/10/90 "Guide to the SLATEC
C           Common Mathematical Library" by Victor A. Bandy
C   991206  indented the code, made do-loops endo on enddo
C           by M. Berndt
C***END PROLOGUE  MGSADD
C
C ==========================================================================

      IMPLICIT NONE

      INCLUDE  'BMG_parameters.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER  II, JJ, NABD1, NABD2
      REAL*8   ABD(NABD1,NABD2), BBD(NABD2), Q(II,JJ), QF(II,JJ)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)
      INTEGER  BMG_iPARMS(NBMG_iPARMS)

C ----------------------------
C     Local Declarations
C
      INTEGER  I, I1, I2, J, J1, KK, N, INFO, BC

C =========================================================================

      BC = BMG_iPARMS(id_BMG2_BC)

C -------------------------------------------
C     Copy the RHS into LAPACK
C -------------------------------------------

      I1=II-1
      J1=JJ-1

      I2=I1-1
      N=I2*(J1-1)

      KK=0
      DO J=2,J1
         DO I=2,I1
            KK=KK+1
            BBD(KK)=QF(I,J)
         ENDDO
      ENDDO

C -------------------------------------------
C     Solve with LAPACK routine DPBTRS
C -------------------------------------------

      INFO = 0
      CALL DPBTRS ('U', KK, I1, 1, ABD, NABD1, BBD, NABD2, INFO) 

      IF (INFO .NE. 0) THEN

         IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'Coarse grid solve failed!'
            WRITE(*,510) 'INFO = ', INFO
         END IF

         BMG_iPARMS(id_BMG2_Ext_Err_Code) = INFO
         CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,20)
         PRINT *
         STOP

      ENDIF

C -------------------------------------------
C     Copy the solution back from LAPACK
C -------------------------------------------

      KK=0
      DO J=2,J1
         DO I=2,I1
            KK=KK+1
            Q(I,J)=BBD(KK)
         ENDDO
      ENDDO

C =========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG2_SymStd_SOLVE_cg.f',/,5X,A)
 510  FORMAT (5X,A,1X,I3)

C ===========================================

      RETURN
      END
