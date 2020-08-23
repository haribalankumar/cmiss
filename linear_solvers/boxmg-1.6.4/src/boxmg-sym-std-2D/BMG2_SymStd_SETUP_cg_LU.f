      SUBROUTINE BMG2_SymStd_SETUP_cg_LU(
     &                SO, II, JJ, NStncl, ABD, NABD1, NABD2,
     &                BMG_IOFLAG, BMG_iPARMS
     &                )

C ==========================================================================
C
C***BEGIN PROLOGUE  MGSSET
C***SUBSIDIARY
C***PURPOSE  mgsset sets up the matrix on the coarsest grid, and
C            using the linpack routine it forms the l-u decomposition
C            of the matrix.
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
C   SO        Refer to BOXMG.
C   II        Number of grid points in x direction, including
C             two fictitious points.
C   JJ        Number of grid points in y direction, including
C             two fictitious points.
C   ABD       Refer to BOXMG.
C
C***ROUTINES CALLED  SPBFA
C***REVISION HISTORY  (YYMMDD)
C   830925  DATE WRITTEN
C   900627  modified to conform to the 4/10/90 "Guide to the SLATEC
C           Common Mathematical Library" by Victor A. Bandy
C   991206  indented the code, made do-loops end on enddo
C           by M. Berndt
C***END PROLOGUE  MGSSET
C

C ==========================================================================

      IMPLICIT NONE

C -----------------------------
C     Includes
C
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_parameters.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER  II, JJ, NABD1, NABD2, NStncl
      REAL*8   ABD(NABD1,NABD2), SO(II,JJ,NStncl)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)
      INTEGER  BMG_iPARMS(NBMG_iPARMS)

C ----------------------------
C     Local Declarations
C
      INTEGER  I, INFO, I1, I2, J, J1, KK, N, BC

C ==========================================================================
      
      BC = BMG_iPARMS(id_BMG2_BC)

C -------------------------------------------------------
C     Copy the operator on the coarsest grid into ABD 
C -------------------------------------------------------

      I1=II-1
      J1=JJ-1

      I2=I1-1
      N=I2*(J1-1)

      IF ( NStncl.EQ.5 ) THEN
         KK=0
         DO J=2,J1
            DO I=2,I1
               KK=KK+1
               ABD(II,KK)=SO(I,J,KO)
               ABD(I1,KK)=-SO(I,J,KW)
               ABD(3,KK)=-SO(I+1,J,KNW)
               ABD(2,KK)=-SO(I,J,KS)
               ABD(1,KK)=-SO(I,J,KSW)
            ENDDO
         ENDDO
      ELSEIF ( NStncl.EQ.3 ) THEN
         KK=0
         DO J=2,J1
            DO I=2,I1
               KK=KK+1
               ABD(II,KK)=SO(I,J,KO)
               ABD(I1,KK)=-SO(I,J,KW)
               ABD(3,KK)=0
               ABD(2,KK)=-SO(I,J,KS)
               ABD(1,KK)=0
            ENDDO
         ENDDO
      ELSE
         IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 
            WRITE(*,510) 'NEED: NStncl = 3 or 5 '
            WRITE(*,510) 'HAVE: NStncl = ', NStncl
         END IF
         
         CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,18)
         RETURN

      ENDIF

      IF( BC.EQ.BMG_BCs_indef_nonper ) THEN
         ABD(II,1) = 2.D0*ABD(II,1)
      ENDIF         

C -------------------------------------------------------
C     Factor using the LAPACK routine DPBTRF
C -------------------------------------------------------

      CALL DPBTRF('U', KK, I1, ABD, NABD1, INFO) 
      
      IF (INFO .NE. 0) THEN

         IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'Coarse grid Cholesky decomposition failed!'
            WRITE(*,510) 'INFO = ', INFO
         END IF

         BMG_iPARMS(id_BMG2_Ext_Err_Code) = INFO
         CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,17)
         PRINT *
         STOP

      ENDIF



C ==========================================================================

 500    FORMAT (/,'FATAL ERROR: BMG2_SymStd_SETUP_cg_LU.f',/,5X,A)
 510    FORMAT (5X,A,1X,I3)

C ====================

      RETURN
      END
