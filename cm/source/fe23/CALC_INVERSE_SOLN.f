      SUBROUTINE CALC_INVERSE_SOLN(nBody,nHeart,NYNR,S_START,S_END,
     '  dPHIdTAU,dPHIdTAU2,d_PHIDIFF,d_TAU,LAMBDA,LAPL,LAPLSQR,LNORM,
     '  PHIDIFF,RNORM,YP,OPFILE,ERROR,*)

C#### Subroutine: CALC_INVERSE_SOLN
C###  Description:
C###    CALC_INVERSE_SOLN computes a new inverse solution for a given
C###    LAMBDA.  Computes the solution norm and the residual norm.
C###    Inputs:
C###      dPHIdTAU  is (dPHI/dtau)
C###      dPHIdTAU2 is (dPHI/dtau)^T.(dPHI/dtau)
C###      PHIDIFF   is (PHI-PHI^)
C###      d_PHIDIFF is (dPHI/dtau)^T.(PHI-PHI^)
C###      d_TAU     is (dPHI/dtau)^T.tau_k
C###      LAMBDA    is regularisation parameter (^2)
C###      YP(ny,5)  is L^T.L(tau_k) where L is Laplacian
C###      YP(ny,8)  is activation time of previous iteration
C###    Outputs:
C###      YP(ny,1)  is new activation time
C###      LNORM     is the solution norm
C###      RNORM     is the residual norm

C     Note: All indices except the ny have been dropped from NYNR
C           nx dropped from YP

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'time02.cmn'

!     Parameter List
      INTEGER nBody,nHeart,NYNR(0:NY_R_M),S_START,S_END
      REAL*8 dPHIdTAU(nBody,nHeart,NTSM),dPHIdTAU2(nHeart,nHeart),
     '  d_PHIDIFF(nHeart),d_TAU(nBody,NTSM),LAMBDA,
     '  LAPL(NY_TRANSFER_M,NY_TRANSFER_M),
     '  LAPLSQR(NY_TRANSFER_M,NY_TRANSFER_M),LNORM,
     '  PHIDIFF(nBody,NTSM),RNORM,YP(NYM,NIYM)
      CHARACTER ERROR*(*)
      LOGICAL OPFILE

!     Local Variables
      INTEGER IFAIL,iH,iHH,IPIV(nHeart),ny
      REAL*8 SUM,WK(nHeart),WK1D(nHeart),WK2D(nHeart,nHeart)

!     Functions
      REAL*8 CALC_LNORM,CALC_RNORM

      CALL ENTERS('CALC_INVERSE_SOLN',*9999)

C***  WK2D = dPHIdTAU2 + lambda^2 L^tL
C***  WK1D = d_PHIDIFF - lambda^2 L^tL(Tau)
      DO iH = 1,nHeart
        DO iHH = 1,nHeart
          WK2D(iH,iHH) = dPHIdTAU2(iH,iHH) +
     '      LAMBDA*LAMBDA*LAPLSQR(iH,iHH)
        ENDDO
        ny = NYNR(iH)
        WK1D(iH) = d_PHIDIFF(iH) - LAMBDA*LAMBDA*YP(ny,5)
      ENDDO

C***  and then invert WK2D (symmetric routines for speed)
      IFAIL = 1

C LKC 25-SEPT-2002 the work space was not the same size as the array.
C  Not sure what the optimal work space is though.
C     CALL DSYTRF('U',nHeart,WK2D,nHeart,IPIV,WK,NY_TRANSFER_M,IFAIL)
C
      CALL DSYTRF('U',nHeart,WK2D,nHeart,IPIV,WK,nHeart,IFAIL)
        IF (IFAIL.NE.0) THEN
          WRITE (ERROR,*) '>>IFAIL=',IFAIL,' in DGETRF'
          GOTO 9999
        ENDIF
      CALL DSYTRI('U',nHeart,WK2D,nHeart,IPIV,WK,IFAIL)
        IF (IFAIL.NE.0) THEN
          WRITE (ERROR,*) '>>IFAIL=',IFAIL,' in DGETRI'
          GOTO 9999
        ENDIF
C     Fill in symmetric inverse from calculated values
      DO iH = 1,nHeart
        DO iHH = iH+1,nHeart
          WK2D(iHH,iH) = WK2D(iH,iHH)
        ENDDO
      ENDDO

C***  Compute new activation times
      DO iH = 1,nHeart
        ny = NYNR(iH)
        SUM = 0.0d0
        DO iHH = 1,nHeart
          SUM = SUM + WK2D(iH,iHH)*WK1D(iHH)
        ENDDO
        YP(ny,1) = YP(ny,8) + SUM
      ENDDO

C***  Compute solution norm
      LNORM = CALC_LNORM(NYNR,LAPL,YP)
C***  Compute residual norm
      RNORM = CALC_RNORM(nBody,nHeart,NYNR,S_START,S_END,dPHIdTAU,d_TAU,
     '  PHIDIFF,YP)

      IF(OPFILE) THEN
        WRITE(IOFI,'(3E15.8)') LAMBDA,LNORM,RNORM
      ENDIF

      CALL EXITS('CALC_INVERSE_SOLN')
      RETURN
 9999 CALL ERRORS('CALC_INVERSE_SOLN',ERROR)
      CALL EXITS('CALC_INVERSE_SOLN')
      RETURN 1
      END


