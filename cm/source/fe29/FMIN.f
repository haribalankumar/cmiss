      REAL*8 FUNCTION FMIN(AX,BX,FUNCT,M,N,BETA,S,WORK)

C#### Function: FMIN
C###  Type: REAL*8
C###  Description:
C###    Computes the local minimum of the function FUNCT(X).
C###  X = FMIN(AX,BX,FUNCT) attemps to return a value of X which is a
C###  local minimiser of FUNCT(X) in the interval AX < X < BX.
C###  Reference:
C###    Forsythe, Malcolm, and Moler, "Computuer Methods for
C###  Mathematical Computations," Prentice-Hall, 1977.
CC JMB 13-OCT-2000

      IMPLICIT NONE
      INCLUDE 'inver00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER M, N
      REAL*8 AX, BETA(*), BX, FUNCT, S(*), WORK(*)
      EXTERNAL FUNCT
!     Local Varaiables
      INTEGER i
      LOGICAL GS
      REAL*8 A, B, C, D, E, FU, FV, FW, FX, MINF, MINFX, ONE, P, Q, R,
     '  RATIO, TOL1, TOL2, U, V, W, X, XM, ZERO
      PARAMETER (ONE = 1.0d0, ZERO = 0.0d0)

      ! C is the squared inverse of the golden ratio
      C = 0.5d0*(3.0d0 - DSQRT(5.0d0))


      ! Find appromimate location of global minimum using NPOINTS
      ! logrithmically located points.
      X = AX
      MINF = FUNCT(X, M, N, BETA, S, WORK)
      MINFX = X
      RATIO = (BX/AX)**(ONE/(DBLE(NPOINTS) - ONE))
      DO i = (NPOINTS - 1),1,-1
         X = RATIO*X
         FX = FUNCT(X, M, N, BETA, S, WORK)
         IF( FX.LT.MINF ) THEN
          MINF = FX
          MINFX = X
        ENDIF
      ENDDO

      ! Identify minimum region
      A = MAX(AX,MINFX/RATIO)
      B = MIN(BX,MINFX*RATIO)
      V = A + C*(B - A)
      W = V
      X = V
      D = ZERO
      E = ZERO
      FX = FUNCT(X, M, N, BETA, S, WORK)
      FV = FX
      FW = FX

      ! Main loop starts here
      XM = 0.50d0*(A + B)
      TOL1 = LOOSE_TOL*ABS(X) + CONVERG_TOL/3.0d0
      TOL2 = 2.0d0*TOL1

      DO WHILE( ABS(X - XM).GT.(TOL2 - 0.50d0*(B - A)) )
        GS = .TRUE.
        IF( ABS(E).GT.TOL1 ) THEN
          GS = .FALSE.
          ! Fit parabola
          R = (X - W)*(FX - FV)
          Q = (X - V)*(FX - FW)
          P = (X - V)*Q - (X - W)*R
          Q = 2.0d0*(Q - R)
          IF( Q.GT.ZERO ) P = -P
          Q = ABS(Q)
          R = E
          E = D

          ! Is the parabola aceptable
          IF(( ABS(P).LT.ABS(0.50d0*Q*R)).AND.(P.GT.Q*(A - X)).AND.
     '      (P.LT.Q*(B - X)) ) THEN
            ! A parabolic interpolation step
            D = P/Q
            U = X + D

            ! F must not be evaluated too close to A or B
            IF( ((U - A).LT.TOL2).OR.((B - U).LT.TOL2) ) THEN
              D = SIGN(TOL1,XM - X)
            ENDIF
          ELSE
            ! Not acceptable, must do a golden section step
            GS = .TRUE.
          ENDIF
        ENDIF
        IF( GS ) THEN
          ! A golden-section step is required
          IF( X.GE.XM ) E = A - X
          IF( X.LT.XM ) E = B - X
          D = C*E
        ENDIF

        ! F must not be evaluated too close to X
        IF( ABS(D).GE.TOL1 ) U = X + D
        IF( ABS(D).LT.TOL1 ) U = X + SIGN(TOL1,D)
        FU = FUNCT(U, M, N, BETA, S, WORK)

        ! Update A, B, V, W and X
        IF( FU.LE.FX ) THEN
          IF( U.GE.X ) A = X
          IF( U.LT.X ) B = X
          V = W
          FV = FW
          W = X
          FW = FX
          X = U
          FX = FU
        ELSE
          IF( U.LT.X ) A = U
          IF( U.GE.X ) B = U
          IF( (FU.LE.FW).OR.(W.EQ.X) ) THEN
            V = W
            FV = FW
            W = U
            FW = FU
          ELSEIF( (FU.LE.FV).OR.(V.EQ.X).OR.(V.EQ.W) ) THEN
            V = U
            FV = FU
          ENDIF
        ENDIF
        XM = 0.50d0*(A + B)
        TOL1 = LOOSE_TOL*ABS(X) + CONVERG_TOL/3.0d0
        TOL2 = 2.0d0*TOL1
      ENDDO
      FMIN = X

      RETURN
      END

C---------------------------------------------------------------------
