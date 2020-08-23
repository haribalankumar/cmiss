      REAL*8 FUNCTION CALC_LNORM(NYNR,LAPL,YP)

C#### Function: CALC_LNORM
C###  Description:
C###    Compute the solution norm, which is the norm of L(tau), where
C###      L   is the Laplacian (or identity if no regularisation)
C###      tau is the current activation times
C###    Inputs:
C###      LAPL      is surface laplacian
C###      YP(ny,1)  is current activation time
C###    Returns L(tau) in YP(ny,6).

C     Note: All indices except the ny have been dropped from NYNR
C           nx dropped from YP

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'inver00.cmn'

!     Parameter List
      INTEGER NYNR(0:NY_R_M)
      REAL*8 LAPL(NY_TRANSFER_M,NY_TRANSFER_M),YP(NYM,NIYM)

!     Local Variables
      INTEGER no_nynr,no_nynr1,ny,ny1
      REAL*8 LNORM

C***  Compute solution norm
      LNORM = 0.0d0
      DO no_nynr = 1,NYNR(0)
        ny = NYNR(no_nynr)
        IF(IREGULARISE.EQ.1) THEN       !no regularisation
          YP(ny,6) = YP(ny,1)
        ELSE IF(IREGULARISE.EQ.2) THEN  !surface Laplacian
          YP(ny,6) = 0.0d0
          DO no_nynr1 = 1,NYNR(0)
            ny1 = NYNR(no_nynr1)
            YP(ny,6) = YP(ny,6) + LAPL(no_nynr,no_nynr1)*YP(ny1,1)
          ENDDO
        ENDIF
        LNORM = LNORM + YP(ny,6)**2
      ENDDO

      CALC_LNORM = SQRT(LNORM)
      RETURN
      END


