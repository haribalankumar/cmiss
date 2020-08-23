      SUBROUTINE BMG2_PerSymStd_SOLVE_cg(
     &                Q, QF, II, JJ, ABD, BBD, NABD1, NABD2, IBC 
     &                )
C
C***BEGIN PROLOGUE  MGSADP
C***SUBSIDIARY
C***PURPOSE  mgsaddp does a direct solve on the coarsest
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
C   AIBC       Refer to BOXMG.
C***INPUT/OUTPUT
C   ABD       Refer to BOXMG.
C   BBD       Refer to BOXMG.
C   NABD1     Refer to BOXMG.
C***OUTPUT
C   Q         Refer to BOXMG
C
C***ROUTINES CALLED  SPBSL, SPOSL
C***REVISION HISTORY  (YYMMDD)
C   830925  DATE WRITTEN
C   900627  modified to conform to the 4/10/90 "Guide to the SLATEC
C           Common Mathematical Library" by Victor A. Bandy
C***END PROLOGUE  MGSADP
C
      IMPLICIT NONE

C     CALLING ARGUMENTS

      integer II, JJ, IBC, NABD1, NABD2
      real*8 ABD(NABD1,NABD2), BBD(NABD2), Q(II,JJ), QF(II,JJ)

C     LOCAL VARIABLES

      integer I, I1, I2, AIBC, J, J1, KK, N, INFO
      real*8 C, CINT, QINT, RZERO
C
C   direct solve on coarsest grid
C
C***FIRST EXECUTABLE STATEMENT  MGSADP
C

      RZERO = 0

      AIBC=IABS(IBC)
      I1=II-1
      J1=JJ-1
      I2=I1-1
      N=I2*(J1-1)
      KK=0

      DO 11 J=2,J1
         DO 10 I=2,I1
            KK=KK+1
            BBD(KK)=QF(I,J)
 10      CONTINUE
 11   CONTINUE

C
C   spbsl and sposl are linpack routines.
C
      IF( AIBC.NE.1 .AND. AIBC.NE.2 .AND. AIBC.NE.3 ) THEN
         CALL DPBTRS('U',KK,I1,1,ABD,NABD1,BBD,NABD2,INFO)
      ELSE
         CALL DPOTRS('U',KK,1,ABD,NABD1,BBD,NABD2,INFO)
      ENDIF
      
      KK=0
      DO 21 J=2,J1
         DO 20 I=2,I1
            KK=KK+1
            Q(I,J)=BBD(KK)
 20      CONTINUE
 21   CONTINUE

      IF ( IBC.GE.0 ) GO TO 60

      CINT=RZERO
      QINT=RZERO
      DO 71 J=2,J1
         DO 70 I=2,I1
            QINT=QINT+Q(I,J)
            CINT=CINT+1
 70      CONTINUE
 71   CONTINUE
      C=-QINT/CINT
      DO 81 J=2,J1
         DO 80 I=2,I1
            Q(I,J)=Q(I,J)+C
 80      CONTINUE
 81   CONTINUE

 60   CONTINUE

      IF( AIBC.NE.1 .AND. AIBC.NE.3 ) GO TO 40

      DO 30 I=2,I1
         Q(I,JJ)=Q(I,2)
         Q(I,1)=Q(I,J1)
 30   CONTINUE

 40   CONTINUE

      IF ( AIBC.NE.2 .AND. AIBC.NE.3 ) RETURN

      DO 50 J=2,J1
         Q(II,J)=Q(2,J)
         Q(1,J)=Q(I1,J)
 50   CONTINUE

      IF ( AIBC.NE.3 ) RETURN

      Q(1,1)   = Q(I1,J1)
      Q(II,1)  = Q(2,J1)
      Q(1,JJ)  = Q(I1,2)
      Q(II,JJ) = Q(2,2)

      RETURN
      END
