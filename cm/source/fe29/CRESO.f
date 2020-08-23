      SUBROUTINE CRESO(CONSTR,M,N,NRHS,U,LDU,SM,LDSM,B,LDB,
     '  REG_PARAMETER,WORK,LWORK,ERROR,*)

C#### Subroutine: CRESO
C###  Description
C###    Computes the regularisation parameters for the CRESO criterion.
C###  Reference:
C###    P. R. Johnston, and R. M. Gulrajani, "A New Method for
C###  Regularization Parameter Determination in the Inverse Problem of
C###  Electrocardiography," IEEE Trans. Biomed. Eng., vol. 44 pp. 19-39,
C###  1997.
C###  Note:
C###    LWORK >= 2*min(M,N)
CC JMB 13-OCT-2000

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
      REAL*8 FMIN, CRESOFUN
      EXTERNAL CRESOFUN

      CALL ENTERS('CRESO',*9999)

      ! Initilaisation
      MN = MIN(M,N)
      IF( LWORK.LT.2*MN ) GOTO 9999
      CALL GSVALUES(CONSTR, MN, SM, LDSM, WORK(1), ERROR, *9999)

      IF( ISTABILISE.EQ.3 ) THEN
        ! Tikhonov
        AX = MAX(WORK(MN),DSQRT(LOOSE_TOL))
        BX = WORK(1)
        DO j = 1,NRHS
          CALL FOURIERCOEFFS(CONSTR, M, MN, U, LDU, B(1,j),
     '     WORK(MN + 1), ERROR, *9999)
          REG_PARAMETER(j) = FMIN(AX, BX, CRESOFUN, M, N, WORK(MN + 1),
     '      WORK(1), %VAL(0))
        ENDDO
      ENDIF

      CALL EXITS('CRESO')
      RETURN
 9999 CALL ERRORS('CRESO',ERROR)
      CALL EXITS('CRESO')
      RETURN 1
      END

C---------------------------------------------------------------------
