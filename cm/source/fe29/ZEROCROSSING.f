      SUBROUTINE ZEROCROSSING(CONSTR,M,N,NRHS,U,LDU,SM,LDSM,B,LDB,
     '  REG_PARAMETER,WORK,LWORK,ERROR,*)

C#### Subroutine: ZEROCROSSING
C###  Description:
C###    Computes the regularisation parameters for the zero-crossing
C###  criterion.
C###  Reference:
C###    P. R. Johnston, and R. M. Gulrajani, "A New Method for
C###  Regularization Parameter Determination in the Inverse Problem of
C###  Electrocardiography," IEEE Trans. Biomed. Eng., vol. 44 pp. 19-39,
C###  1997.
C###  Note:
C###    LWORK >= 4*min(M,N) + 1
CC JMB 03-FEB-2000

      IMPLICIT NONE
      INCLUDE 'inver00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER LDB, LDSM, LDU, LWORK, M, N, NRHS
      REAL*8 B(LDB,*), REG_PARAMETER(*), SM(LDSM,*), U(LDU,*), WORK(*)
      CHARACTER CONSTR, ERROR*(*)
!     Local Variables
      INTEGER j, MN
      REAL*8 AX, BX
!     Functions
      REAL*8 FMIN, RESIDFUN, FZERO, ZEROFUN
      EXTERNAL ZEROFUN

      CALL ENTERS('ZEROCROSSING',*9999)

      ! Initialisation
      MN = MIN(M,N)
      IF( LWORK.LT.(4*MN + 1) ) GOTO 9999
      CALL GSVALUES(CONSTR, MN, SM, LDSM, WORK(1), ERROR, *9999)

      IF( ISTABILISE.EQ.3 ) THEN
        ! Tikhonov
        AX = MAX(WORK(MN),DSQRT(LOOSE_TOL))
        BX = WORK(1)
        DO j = 1,NRHS
          CALL FOURIERCOEFFS(CONSTR, M, MN, U, LDU, B(1,j),
     '      WORK(MN + 1), ERROR, *9999)
          WORK(2*MN + 1) = RESIDFUN(M, N, B(1,j), WORK(MN + 1))
          REG_PARAMETER(j) = FMIN(AX, BX, ZEROFUN, M, N, WORK(MN + 1),
     '      WORK(1), WORK(2*MN + 1))
          IF( DABS(REG_PARAMETER(j) - DSQRT(LOOSE_TOL)).LE.
     '      LOOSE_TOL ) THEN
            REG_PARAMETER(j) = WORK(1)
          ELSE
            REG_PARAMETER(j) = FZERO(DSQRT(LOOSE_TOL), REG_PARAMETER(j),
     '        ZEROFUN, M, N, WORK(MN + 1), WORK(1), WORK(2*MN + 1))
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('ZEROCROSSING')
      RETURN
 9999 CALL ERRORS('ZEROCROSSING',ERROR)
      CALL EXITS('ZEROCROSSING')
      RETURN 1
      END

C---------------------------------------------------------------------
