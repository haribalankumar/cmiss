      SUBROUTINE LCURVE(CONSTR,M,N,NRHS,U,LDU,SM,LDSM,B,LDB,
     '  REG_PARAMETER,WORK,LWORK,ERROR,*)

C#### Subroutine: LCURVE
C###  Description:
C###    Computes the regularisation parameters for the L-curve
C###  criterion.
C###  Reference:
C###   P. C. Hansen, "Regularization tools, a Matlab package for
C###  analysis and solution of discrete ill-posed problems," UNI.C, 1998.
C###  Note:
C###    LWORK >= 8*min(M,N) + K + 2
CC JMB 13-OCT-2000

      IMPLICIT NONE
      INCLUDE 'inver00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER LDB, LDSM, LDU, LWORK, M, N, NRHS
      REAL*8 B(LDB,*), REG_PARAMETER(*), SM(LDSM,*), U(LDU,*), WORK(*)
      CHARACTER CONSTR, ERROR*(*)
!     Local Varaibles
      INTEGER DEG, K, LTWORK, NB, Q
      PARAMETER (K = 3, Q = 2, DEG = 2, NB = 64,
     '  LTWORK = (2*Q + 1) + (2*Q + 1)*NB)
      INTEGER i, IFAIL, iKMAX, iL2MIN, j, l, MN,
     '  NSV, p, V(2*Q+1)
      REAL*8 A(2*Q+1,DEG+1), AX, BX, C(2*Q+1,2), CURV, DELTA, KAPPA,
     '  KAPPAMAX, K1, K2, L2MIN, L2, ONE, SPX(0:NPOINTS,0:K-1),
     '  SPY(0:NPOINTS,0:K-1), TWORK(LTWORK), X, XCORNER,
     '  YCORNER, ZERO
      PARAMETER (ONE = 1.0d0, ZERO = 0.0d0)
!     Functions
      REAL*8 BVALUE, FMIN, LCFUN, RESIDFUN
      EXTERNAL LCFUN

      CALL ENTERS('LCURVE',*9999)

      ! Initilaisation
      MN = MIN(M,N)
      IF( LWORK.LT.(8*MN + K + 2) ) GOTO 9999
      CALL GSVALUES(CONSTR, MN, SM, LDSM, WORK(1), ERROR, *9999)

      IF( ISTABILISE.EQ.2 ) THEN
        ! TGSVD
        DO l = 1,NRHS
          CALL FOURIERCOEFFS(CONSTR, M, MN, U, LDU, B(1,l),
     '      WORK(MN + 1), ERROR, *9999)
          DELTA = RESIDFUN(M, N, B(1,l), WORK(MN + 1))
          DO i = 1,MN
            WORK(2*MN + i) = WORK(MN + i)/WORK(i)
          ENDDO
          WORK(3*MN + 1) = WORK(2*MN + 1)**2
          DO i = 2,MN
            WORK(3*MN + i) = WORK(3*MN + i - 1) + WORK(2*MN + i)**2
          ENDDO

          WORK(5*MN) = DELTA
          DO i = MN - 1,1,-1
            WORK(4*MN + i) = WORK(4*MN + i + 1) + WORK(MN + i + 1)**2
          ENDDO

          DO i = 1,MN
            WORK(3*MN + i) = DSQRT(WORK(3*MN + i))
            WORK(4*MN + i) = DSQRT(WORK(4*MN + i))
          ENDDO

          ! Set threshold for skipping very small singular values
          NSV = MN
          DO WHILE( WORK(NSV).LE.ZERO_TOL .AND. NSV.GE.1 )
            NSV = NSV-1
          ENDDO
          DO i = 1,NSV
            WORK(3*MN + i) = DLOG(WORK(3*MN + i))
            WORK(4*MN + i) = DLOG(WORK(4*MN + i))
            WORK(5*MN + i) = WORK(3*MN + i)
            WORK(6*MN + i) = WORK(4*MN + i)
          ENDDO

          ! Perform local smoothing to the points j-Q,j+Q.
          DO i = -Q,Q
            V(i+Q+1) = i
          ENDDO
          DO j = Q + 1,NSV - Q - 1
            DO i = 1,2*Q + 1
              A(i,1) = ONE
              DO p = 2,DEG + 1
                A(i,p) = A(i,p-1)*V(i)
              ENDDO
            ENDDO
            DO i = 1,2*Q + 1
              C(i,1) = WORK(3*MN + j + V(i))
              C(i,2) = WORK(4*MN + j + V(i))
            ENDDO
            IFAIL = 1
            CALL DGELS('N', 2*Q + 1, DEG + 1, 2, A, 2*Q + 1, C,
     '        2*Q + 1, TWORK, LTWORK, IFAIL)
              WORK(5*MN + j) = C(1,1)
              WORK(6*MN + j) = C(1,2)
          ENDDO

          ! Fit a 2-D spline curve to the smoothed discrete L-curve
          DO i = 1,NSV + K + 1
            WORK(7*MN + i) = DBLE(i)
          ENDDO

          ! Evaluate spline between T(K+1)..T(N+1)
          DO i=0,NPOINTS
             X=DBLE(WORK(7*MN + K + 1)) + DBLE(i)*
     '         DBLE(WORK(7*MN + NSV + 1) - WORK(7*MN + K + 1))/
     '         DBLE(NPOINTS)
            DO j = 0,2 ! Value + first and second derivative
              SPX(i,j) = BVALUE(WORK(7*MN + 1), WORK(6*MN + 1), NSV,
     '          K + 1, X, j)
              SPY(i,j) = BVALUE(WORK(7*MN + 1), WORK(5*MN + 1), NSV,
     '          K + 1, X, j)
            ENDDO
          ENDDO

          ! Compute the corner of the discretised spline
          KAPPAMAX = ZERO
          DO i = 0,NPOINTS
            K1 = SPX(i,1)*SPY(i,2) - SPX(i,2)*SPY(i,1)
            K2 = (SPX(i,1)**2 + SPY(i,1)**2)**(1.5d0)
            IF( DABS(K2).GT.ZERO_TOL ) THEN
              KAPPA = -K1/K2
              IF( KAPPA.GT.KAPPAMAX ) THEN
                KAPPAMAX = KAPPA
                iKMAX = i
              ENDIF
            ENDIF
          ENDDO
          XCORNER = SPX(iKMAX,0)
          YCORNER = SPY(iKMAX,0)

          ! Locate corner of discrete L-curve
          IF( KAPPAMAX.LE.ZERO ) THEN
            REG_PARAMETER(l) = ONE
          ELSE
            i = 1
            iL2MIN = i
            L2MIN = (WORK(4*MN + i) - XCORNER)**2 + (WORK(3*MN + i) -
     '        YCORNER)**2
            DO WHILE( iL2MIN.EQ.i )
              i = i + 1
              L2 = (WORK(4*MN + i) - XCORNER)**2 + (WORK(3*MN + i) -
     '          YCORNER)**2
              IF( L2.LT.L2MIN ) THEN
                L2MIN = L2
                iL2MIN = i
              ENDIF
            ENDDO
            REG_PARAMETER(l) = DBLE(iL2MIN)
          ENDIF
        ENDDO
      ELSEIF(ISTABILISE.EQ.3) THEN
        ! Tikhonov
        IF(ICONSTRAINT.NE.4) THEN
          CURV = 0.5d0
        ELSE
          CURV = ZERO
        ENDIF
        AX = MAX(WORK(MN),DSQRT(LOOSE_TOL))
        BX = WORK(1)
        DO j = 1,NRHS
          CALL FOURIERCOEFFS(CONSTR, M, MN, U, LDU, B(1,j),
     '      WORK(MN + 1), ERROR, *9999)
          WORK(2*MN + 1) = RESIDFUN(M, N, B(1,j), WORK(MN + 1))
          REG_PARAMETER(j) = FMIN(AX, BX, LCFUN, M, N, WORK(MN + 1),
     '      WORK(1), WORK(2*MN + 1))
          KAPPAMAX = -ONE*LCFUN(REG_PARAMETER(j), M, N, WORK(MN + 1),
     '      WORK(1), WORK(2*MN + 1))
          IF(KAPPAMAX.LE.CURV) THEN
            REG_PARAMETER(j) = WORK(1)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('LCURVE')
      RETURN
 9999 CALL ERRORS('LCURVE',ERROR)
      CALL EXITS('LCURVE')
      RETURN 1
      END

C---------------------------------------------------------------------
