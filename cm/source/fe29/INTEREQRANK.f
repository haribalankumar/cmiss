      SUBROUTINE INTEREQRANK(MT,K,RANK,SB,IWORK,WORK,LWORK,ERROR,*)

C#### Function: INTEREQRANK
C###  Description:
C###    Calculates the inter-equation truncation rank.
C###  Note:
C###   LWORK >= 4*min(M,T) + min(M,T)*(K + 1) + NEST*(10 + 3*K)
C###   where NEST = min(M,T) + K + 1
CC JMB 13-OCT-2000

      IMPLICIT NONE
!     Parameter List
      INTEGER IWORK(*), K, LWORK, MT, RANK
      REAL*8 SB(*), WORK(*)
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 FP, ONE, S, TOL, XB, XE
      INTEGER i, IFAIL, IC, IT, IW, IWRK, IY, LIWRK, N, NEST
      DATA TOL / 1.D-3 /
      PARAMETER (ONE = 1.0d0, S = 0.2d0)

      CALL ENTERS('INTEREQRANK',*9999)

      ! Initialisation
      NEST = MT + K + 1
      IF( LWORK.LT.(4*MT + MT*(K + 1) + NEST*(10 + 3*K)) ) GOTO 9999
      IY = MT
      IW = IY + MT
      IT = IW + MT
      IC = IT + NEST
      IWRK = IC + NEST
      LIWRK = MT*(K + 1) + NEST*(7 + 3*K)

      DO i = 1,MT
        WORK(i) = DBLE(i)
        WORK(IY + i) = DLOG(SB(i))
        WORK(IW + i) = ONE
      ENDDO

      XB = ONE
      XE = DBLE(MT)
      N = 0

      ! NAG E02BEF equivalent
      CALL CURFIT(0, MT, WORK(1), WORK(IY + 1), WORK(IW + 1), XB, XE, K,
     '  S, NEST, N, WORK(IT + 1), WORK(IC + 1), FP, WORK(IWRK + 1),
     '  LIWRK, IWORK, IFAIL)
      IF( IFAIL.NE.0 ) GOTO 9999

      ! Evaluate the curvature at each singular value
      CALL SPLDER(WORK(IT + 1), N, WORK(IC + 1), K, 2, WORK(1),
     '  WORK(IY + 1), MT, WORK(IWRK + 1), IFAIL)
      IF( IFAIL.NE.0 ) GOTO 9999

      i = 1
      DO WHILE( DABS(WORK(IY + i)).GT.TOL )
        i = i + 1
      ENDDO
      RANK = i - 1

      CALL EXITS('INTEREQRANK')
      RETURN
 9999 CALL ERRORS('INTEREQRANK',ERROR)
      CALL EXITS('INTEREQRANK')
      RETURN 1
      END

C---------------------------------------------------------------------
