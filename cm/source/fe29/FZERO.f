      REAL*8 FUNCTION FZERO(AX,BX,FUNCT,M,N,BETA,S,WORK)

C#### Function: FZERO
C###  Type: REAL*8
C###  Description:
C###    Computes the zero of the function F(X) in the interval AX,BX.
C###  X = FZERO(AX,BX,FUNCT) attemps to return a value of X which is a
C###  zero of FUNCT(X) in the interval AX < X < BX.
C###  Reference:
C###    Forsythe, Malcolm, and Moler, "Computuer Methods for
C###  Mathematical Computations," Prentice-Hall, 1977.
CC JMB-13-OCT-2000

      IMPLICIT NONE
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER M, N
      REAL*8 AX, BETA(*), BX, FUNCT, S(*), WORK(*)
!     Local Variables
      REAL*8 A, B, C, D, E, FA, FB, FC, ONE, P, Q, R, T, TOL1, XM, ZERO
      PARAMETER (ONE = 1.0d0, ZERO = 0.0d0)

      ! Initalisation
      A = AX
      B = BX
      FA = FUNCT(A, M, N, BETA, S, WORK)
      FB = FUNCT(B, M, N, BETA, S, WORK)
      FC = FB
      DO WHILE ( FB.NE.0 )
        IF( (FB*(FC/ABS(FC))).GT.ZERO ) THEN
          C = A
          FC = FA
          D = B - A
          E = D
        ENDIF
        IF( ABS(FC).LT.ABS(FB) ) THEN
          A = B
          B = C
          C = A
          FA = FB
          FB = FC
          FC = FA
        ENDIF

        ! Convergence test and possible exit
        TOL1 = 2.0d0*LOOSE_TOL*ABS(B) + 0.50d0*CONVERG_TOL
        XM = 0.50d0*(C - B)
        IF( (ABS(XM).LE.TOL1).OR.(FB.EQ.ZERO) ) GOTO 10

        ! Choose bisection or interpolation
        IF( (ABS(E).LT.TOL1).OR.(ABS(FA).LE.ABS(FB)) ) THEN
          D = XM
          E = D
        ELSE
          ! Interpolation
          T = FB/FA
          IF( A.EQ.C ) THEN
            ! Linear interpolation
            P = 2.0d0*XM*T
            Q = ONE - T
          ELSE
            ! Inverse quadratic interpolation
            Q = FA/FC
            R = FB/FC
            P = T*(2.0d0*XM*Q*(Q - R) - (B - A)*(R - ONE))
            Q = (Q - ONE)*(R - ONE)*(T - ONE)
          ENDIF
          IF( P.GT.ZERO ) THEN
            Q = -Q
          ELSE
            P = ABS(P)
          ENDIF

          ! Is the interpolated point acceptable
          IF( ((2.0d0*P).LT.(3.0d0*XM*Q - ABS(TOL1*Q)) ).OR.
     '      (P.LT.ABS(0.50d0*E*Q)) ) THEN
            E = D
            D = P/Q
          ELSE
            D = XM
            E = D
          ENDIF
        ENDIF

        ! Complete step
        A = B
        FA = FB
        IF( ABS(D).GT.TOL1 ) THEN
          B = B + D
        ELSE
          B = B + SIGN(TOL1,XM)
        ENDIF
        FB = FUNCT(B, M, N, BETA, S, WORK)
      ENDDO
   10 FZERO = B

      RETURN
      END

C---------------------------------------------------------------------
