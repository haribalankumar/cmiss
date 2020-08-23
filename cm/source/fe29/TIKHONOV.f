      SUBROUTINE TIKHONOV(CONSTR,M,N,U,LDU,SM,LDSM,VT,LDVT,X,B,LAMBDA,
     '  WORK,LWORK,ERROR,*)

C#### Subroutine: TIKHONOV
C###  Description:
C###    Computes the Tikhonov regularised solution.
C###  Reference:
C###   A. N. Tikhonov & V. Y. Arsenin, "Solutions of Ill-Posed Problems",
C###  Wiley, 1977.
C###   P. C. Hansen, "Regularization tools, a Matlab package for
C###  analysis and solution of discrete ill-posed problems," UNI.C, 1998.
C###  Note:
C###   LWORK >= min(M,N) + N
CC JMB 10-OCT-2000

      IMPLICIT NONE
!     Parameter List
      INTEGER LDSM, LDU, LDVT, LWORK, M, N
      REAL*8 B(*), LAMBDA, SM(LDSM,*), U(LDU,*), VT(LDVT,*), WORK(*),
     '  X(*)
      CHARACTER CONSTR, ERROR*(*)
!     Local Variables
      INTEGER i, MN
      REAL*8 ONE,ZERO
      PARAMETER (ONE = 1.0d0,ZERO = 0.0d0)
!     Functions
      LOGICAL LSAME

      CALL ENTERS('TIKHONOV',*9999)

      ! Initialisation
      MN = MIN(M,N)
      IF( LWORK.LT.(MN + N) ) GOTO 9999

      CALL DGEMV('T', M, MN, ONE, U, LDU, B, 1, ZERO, WORK(1), 1)
      IF( LSAME(CONSTR, 'I') ) THEN
        ! Zero-order Tikhonov
        DO i = 1,MN
          WORK(i) = SM(i,1)*WORK(i)/(SM(i,1)**2 + LAMBDA**2)
        ENDDO
        CALL DGEMV('T', MN, N, ONE, VT, LDVT, WORK(1), 1, ZERO, X, 1)
      ELSE
        ! First or second order Tikhonov
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
        DO i = 1,MN
          WORK(i) = SM(i,1)*WORK(i)/(SM(i,1)**2 + LAMBDA**2*SM(i,2)**2)
        ENDDO
        CALL DGEMV('T', M, N, ONE, VT, LDVT, WORK(1), 1, ONE, X ,1)
      ENDIF

      CALL EXITS('TIKHONOV')
      RETURN
 9999 CALL ERRORS('TIKHONOV',ERROR)
      CALL EXITS('TIKHONOV')
      RETURN 1
      END

C---------------------------------------------------------------------
