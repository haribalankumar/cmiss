      SUBROUTINE GCV(CONSTR,M,N,NRHS,U,LDU,SM,LDSM,B,LDB,REG_PARAMETER,
     '  WORK,LWORK,ERROR,*)

C#### Subroutine: GCV
C###  Description:
C###    Computes the regularisation parameters for the GCV criterion.
C###  Reference:
C###   P. C. Hansen, "Regularization tools, a Matlab package for
C###  analysis and solution of discrete ill-posed problems," UNI.C, 1998.
C###  Note:
C###    LWORK >= 4*min(M,N) + 1
CC JMB 13-OCT-2000

      IMPLICIT NONE
      INCLUDE 'inver00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER LDB, LDSM, LDU, LWORK, M, N, NRHS
      REAL*8 B(LDB,*), REG_PARAMETER(*), SM(LDSM,*), U(LDU,*), WORK(*)
      CHARACTER CONSTR, ERROR*(*)
!     Local Variables
      INTEGER i, j, MN
      REAL*8 AX, BX, DELTA, G, MING, ONE
      PARAMETER (ONE = 1.0d0)
!     Functions
      REAL*8 FMIN, RESIDFUN, GCVFUN
      EXTERNAL GCVFUN

      CALL ENTERS('GCV',*9999)

      ! Initilaisation
      MN = MIN(M,N)
      IF( LWORK.LT.(4*MN + 1) ) GOTO 9999
      CALL GSVALUES(CONSTR, MN, SM, LDSM, WORK(1), ERROR, *9999)

      IF( ISTABILISE.EQ.2 ) THEN
        ! TGSVD
        DO j = 1,NRHS
          CALL FOURIERCOEFFS(CONSTR, M, MN, U, LDU, B(1,j),
     '     WORK(MN + 1), ERROR, *9999)
          DELTA = RESIDFUN(M, N, B(1,j), WORK(MN + 1))
          WORK(3*MN - 1) = WORK(2*MN)**2 + DELTA
          DO i = MN - 2,1,-1
            WORK(2*MN + i) = WORK(2*MN + i + 1) + WORK(MN + i + 1)**2
          ENDDO
          MING = 1.D+6
          REG_PARAMETER(j) = ONE
          DO i = 1,MN - 1
            G = WORK(2*MN + i)/DBLE((M - i + (N - MN)))**2
            IF ( G.LT.MING ) THEN
              MING = G
              REG_PARAMETER(j) = DBLE(i)
            ENDIF
          ENDDO
        ENDDO
      ELSEIF( ISTABILISE.EQ.3 ) THEN
        ! Tikhonov
        AX = MAX(WORK(MN),DSQRT(LOOSE_TOL))
        BX = WORK(1)
        DO j = 1,NRHS
          CALL FOURIERCOEFFS(CONSTR, M, MN, U, LDU, B(1,j),
     '     WORK(MN + 1), ERROR, *9999)
          WORK(2*MN + 1) = RESIDFUN(M, N, B(1,j), WORK(MN + 1))
          REG_PARAMETER(j) = FMIN(AX, BX, GCVFUN, M, N, WORK(MN + 1),
     '      WORK(1), WORK(2*MN + 1))
        ENDDO
      ENDIF

      CALL EXITS('GCV')
      RETURN
 9999 CALL ERRORS('GCV',ERROR)
      CALL EXITS('GCV')
      RETURN 1
      END

C---------------------------------------------------------------------
