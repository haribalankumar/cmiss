      SUBROUTINE TGSVD(CONSTR,M,N,U,LDU,SM,LDSM,VT,LDVT,X,B,K,WORK,
     '  LWORK,ERROR,*)

C#### Subroutine: TGSVD
C###  Description:
C###    Computes the truncated GSVD regularised solution.
C###  Reference:
C###   P. C. Hansen, "Regularization tools, a Matlab package for
C###  analysis and solution of discrete ill-posed problems," UNI.C, 1998.
C###  Note:
C###   LWORK >= min(M,N) + N
CC JMB 10-OCT-2000

      IMPLICIT NONE
!     Parameter List
      INTEGER K, LDSM, LDU, LDVT, M, N, LWORK
      REAL*8 B(*), SM(LDSM,*), U(LDU,*), VT(LDVT,*), WORK(*), X(*)
      CHARACTER CONSTR, ERROR*(*)
!     Local Variables
      INTEGER i, MN
      REAL*8 ONE, ZERO
      PARAMETER (ONE = 1.0d0, ZERO = 0.0d0)
!     Functions
      LOGICAL LSAME

      CALL ENTERS('TGSVD',*9999)

      ! Initialisation
      MN = MIN(M,N)
      IF( LWORK.LT.(MN + N) ) GOTO 9999

      CALL DGEMV('T', M, MN, ONE, U, LDU, B, 1, ZERO, WORK(1), 1)

      DO i = 1,MN
        WORK(i) = WORK(i)/SM(i,1)
      ENDDO

      IF( LSAME(CONSTR, 'I') ) THEN
        ! TSVD
        CALL DGEMV('T', K, N, ONE, VT, LDVT, WORK(1) , 1, ZERO, X, 1)
      ELSE
        ! TGSVD
        IF( MN.EQ.N ) THEN
          DO i = 1,N
            WORK(MN + i) = ZERO
          ENDDO
        ELSE
          CALL DGEMV('T', M, N - MN, ONE, U(1,MN + 1), LDU, B, 1, ZERO,
     '      WORK(MN + 1), 1)
          CALL DGEMV('T', N - MN, N, ONE, VT(1,MN + 1), LDVT,
     '      WORK(MN + 1), 1, ZERO, X, 1)
        ENDIF
        IF( K.NE.0 ) THEN
          CALL DGEMV('T', K, N, ONE, VT(1,MN - K + 1), LDVT,
     '      WORK(MN - K + 1), 1, ONE, X, 1)
        ENDIF
      ENDIF

      CALL EXITS('TGSVD')
      RETURN
 9999 CALL ERRORS('TGSVD',ERROR)
      CALL EXITS('TGSVD')
      RETURN 1
      END

C---------------------------------------------------------------------
