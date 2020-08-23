      REAL*8 FUNCTION REFUN(LAMBDA,M,N,BETA,S,WORK)

C#### Function: REFUN
C###  Type: REAL*8
C###  Description:
C###    Auxililary routine for REOPT. Computes the RE function.
CC JMB 21-FEB-2000

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'inver00.cmn'
!     Parameter List
      INTEGER M, N
      REAL*8 BETA(*), LAMBDA, S(*), WORK(*)
!     Local Variables
      INTEGER i, IX, IXREG, IVTB, IZETA, j, MN
      REAL*8 ONE, RE, SCALE, SSQ, ZERO
      PARAMETER (ONE = 1.0d0, ZERO = 0.0d0)
!     Functions
      REAL*8 DNRM2

      ! Initialisation
      MN = MIN(M,N)
      IX = MN*N
      IF( ICOUPLING.EQ.1 ) THEN
        IZETA = IX + N
      ELSE
        IVTB = IX + N*NTST
        IZETA = IVTB + NTST
      ENDIF
      IXREG = IZETA + MN

      DO i = 1,MN
        WORK(IZETA + i) = S(i)*BETA(i)/(S(i)**2 + LAMBDA**2)
      ENDDO
      CALL DGEMV('T', MN, N, ONE, WORK(1), MN, WORK(IZETA + 1), 1, ZERO,
     '  WORK(IXREG + 1), 1)

      IF( ICOUPLING.EQ.1 ) THEN
        CALL DAXPY(N, -ONE, WORK(IX + 1), 1, WORK(IXREG + 1), 1)
        RE = DNRM2(N, WORK(IXREG + 1), 1)
      ELSE
        SCALE = ZERO
        SSQ = ONE
        DO j = 1,NTST
          DO i = 1,N
            WORK(IXREG + N + i) = WORK(IVTB + j)*WORK(IXREG + i)
     '        - WORK(IX + N*(j - 1) + i)
          ENDDO
          CALL DLASSQ(N, WORK(IXREG + N + 1), 1, SCALE, SSQ)
        ENDDO
        RE = SCALE*DSQRT(SSQ)
      ENDIF

      REFUN = RE

      RETURN
      END

C---------------------------------------------------------------------
