      SUBROUTINE QUASIOPT(CONSTR,M,N,NRHS,U,LDU,SM,LDSM,B,LDB,
     '  REG_PARAMETER,WORK,LWORK,ERROR,*)

C#### Subroutine: QUASIOPT
C###  Description:
C###    Computes the regularisation parameters for the quasi-optimality
C###  criterion.
C###  Reference:
C###   P. C. Hansen, "Regularization tools, a Matlab package for
C###  analysis and solution of discrete ill-posed problems," UNI.C, 1998.
C###  Note:
C###    LWORK >= 3*min(M,N)
CC JMB 13-OCT-2000

      IMPLICIT NONE
      INCLUDE 'inver00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER LDB, LDSM, LDU, LWORK, M, N
      REAL*8 B(LDB,*), REG_PARAMETER(*), SM(LDSM,*), U(LDU,*), WORK(*)
      CHARACTER CONSTR, ERROR*(*)
!     Local Variables
      INTEGER i, j, MN, NRHS
      REAL*8 AX, BX, MINQ, ONE, Q
      PARAMETER (ONE = 1.0d0)
!     Functions
      REAL*8 FMIN, QUASIFUN
      EXTERNAL QUASIFUN

      CALL ENTERS('QUASIOPT',*9999)

      ! Initialisation
      MN = MIN(M,N)
      IF( LWORK.LT.3*MN ) GOTO 9999
      CALL GSVALUES(CONSTR, MN, SM, LDSM, WORK(1), ERROR, *9999)

      IF( ISTABILISE.EQ.2 ) THEN
        ! TGSVD
        DO j = 1,NRHS
          CALL FOURIERCOEFFS(CONSTR, M, N, U, LDU, B(1,j),
     '      WORK(MN + 1), ERROR, *9999)
          MINQ = 1.D+6
          REG_PARAMETER(j) = ONE
          DO i = 1,MN
            Q = DABS(WORK(MN + i)/WORK(i))
            IF( Q.LT.MINQ ) THEN
              MINQ = Q
              REG_PARAMETER(j) = DBLE(i)
            ENDIF
          ENDDO
        ENDDO
      ELSEIF( ISTABILISE.EQ.3 ) THEN
        ! Tikhonov
        AX = MAX(WORK(MN),DSQRT(LOOSE_TOL))
        BX = WORK(1)
        DO j = 1,NRHS
          CALL FOURIERCOEFFS(CONSTR, M, N, U, LDU, B(1,j),
     '      WORK(MN + 1), ERROR, *9999)
          REG_PARAMETER(j) = FMIN(AX, BX, QUASIFUN, M, N, WORK(MN + 1),
     '      WORK(1), WORK(2*MN + 1))
        ENDDO
      ENDIF

      CALL EXITS('QUASIOPT')
      RETURN
 9999 CALL ERRORS('QUASIOPT',ERROR)
      CALL EXITS('QUASIOPT')
      RETURN 1
      END

C---------------------------------------------------------------------
