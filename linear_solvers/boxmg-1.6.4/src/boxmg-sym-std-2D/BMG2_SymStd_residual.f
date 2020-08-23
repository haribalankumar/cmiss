      SUBROUTINE BMG2_SymStd_residual( 
     &                K, SO, QF, Q, RES, II, JJ,
     &                KF, IFD, NStncl
     &                )

C ==========================================================================

C***BEGIN PROLOGUE BMG2_SymStd_resl2
C***SUBSIDIARY
C***PURPOSE  BMG_SymStd_resl2 calculates the l2 residual 
C            on the fine grid depending on the value of irelax.
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
C   K         K is the grid number.
C   SO        Refer to BOXMG.
C   QF        Refer to BOXMG.
C   RES       Refer to BOXMG.
C   II        Number of grid points in x direction, including
C             two fictitious points.
C   JJ        Number of grid points in y direction, including
C             two fictitious points.
C   KF        index of the finest grid
C   IFD       Refer to BOXMG
C   IRELAX    Refer to BOXMG.
C***INPUT/OUTPUT
C   Q         Refer to BOXMG.
C***OUTPUT
C   ERR       ERR is the l2 error of the resiuals.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   830925  DATE WRITTEN
C   900627  modified to conform to the 4/10/90 "Guide to the SLATEC
C           Common Mathematical Library" by Victor A. Bandy
C
C   991206  Changed all gotos to if statements, added indentation, and 
C           changed do-loops to end on enddo (M. Berndt)
C
C   000218  Separated the routine mgrlx into relaxation and 
C           residual calculation (M. Berndt)
C
C***END PROLOGUE  BMG2_SymStd_resl2

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

      INTEGER IFD, K, KF
      REAL*8  Q(II,JJ), QF(II,JJ), SO(II,JJ,NStncl), RES(II,JJ), SUM

C ----------------------------
C     Local Declarations
C
      INTEGER I, I1, J, J1, KK

C =========================================================================

      J1=JJ-1
      I1=II-1

C     -------------------------------------------------------------------

      IF ( K.LT.KF .OR. IFD.NE.1 ) THEN
         !
         !  9-point stencil
         ! 
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& FIRSTPRIVATE(I1,J1)
C$OMP& PRIVATE(I,J,KK)
C$OMP& SHARED(Q,QF,RES,SO)
	   DO KK=1,(I1-1)*(J1-1)

		I=MOD(KK-1,I1-1)+2
		J=(KK-1)/(I1-1)+2

		RES(I,J) = QF(I,J)
     &               + SO(I  ,J  ,KW )*Q(I-1,J)
     &               + SO(I+1,J  ,KW )*Q(I+1,J)
     &               + SO(I  ,J  ,KS )*Q(I  ,J-1)
     &               + SO(I  ,J+1,KS )*Q(I  ,J+1)
     &               + SO(I  ,J  ,KSW)*Q(I-1,J-1)
     &               + SO(I+1,J  ,KNW)*Q(I+1,J-1)
     &               + SO(I  ,J+1,KNW)*Q(I-1,J+1)
     &               + SO(I+1,J+1,KSW)*Q(I+1,J+1)
     &               - SO(I  ,J  ,KO )*Q(I  ,J)
	   ENDDO
C$OMP END PARALLEL DO
         !
      ELSE
         !
         !  5-point stencil
         !      
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(I,J,KK)
C$OMP& SHARED(I1,J1,Q,QF,RES,SO)
	   DO KK=1,(I1-1)*(J1-1)

		I=MOD(KK-1,I1-1)+2
		J=(KK-1)/(I1-1)+2

                RES(I,J) = QF(I,J)
     &               + SO(I  ,J  ,KW)*Q(I-1,J)
     &               + SO(I+1,J  ,KW)*Q(I+1,J)
     &               + SO(I  ,J  ,KS)*Q(I  ,J-1)
     &               + SO(I  ,J+1,KS)*Q(I  ,J+1)
     &               - SO(I  ,J  ,KO)*Q(I  ,J)
         ENDDO
C$OMP END PARALLEL DO
         !
      ENDIF
C =========================================================================

      RETURN
      END

