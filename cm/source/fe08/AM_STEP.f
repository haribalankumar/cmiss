      SUBROUTINE AM_STEP(AII,AIO,CONTROL,ERR_CODE,K,KOLD,MAX_ORDER,
     '  MAX_RND,MODEL,NUMBER_EQN,NUM_VARS,NUMBER_STEPS,SIZES,VARIANT,
     '  ALPHA,ARI,ARO,BETA,DERIVED,DY,EPS,EPSILON,ERROR_WEIGHT,G,H,
     '  HOLD,MAX_STEP,PARAMETERS,PHI,PREDICTION,PROTOCOL,PSI,SIGMA,
     '  T,V,W,Y,CRASH,PHASE1,ROUND_CTRL,START,FUNC)

C#### Subroutine: AM_STEP
C###  Description:
C###    Subroutine to perfom the Adams-Moulton step from T to T+H

      IMPLICIT NONE

      INCLUDE 'cell_reserved.inc'

!     Parameter list
      INTEGER AII(*),AIO(*),CONTROL(*),ERR_CODE,K,KOLD,MAX_ORDER,
     '  MAX_RND,MODEL(*),NUMBER_EQN,NUM_VARS,NUMBER_STEPS,SIZES(*),
     '  VARIANT
      REAL*8 ALPHA(MAX_ORDER),ARI(*),ARO(*),BETA(MAX_ORDER),
     '  DERIVED(*),DY(NUM_VARS),EPS,EPSILON,ERROR_WEIGHT(NUMBER_EQN),
     '  G(MAX_ORDER+1),H,HOLD,MAX_STEP,PARAMETERS(*),
     '  PHI(NUMBER_EQN,MAX_ORDER+2+MAX_RND),PREDICTION(NUM_VARS),
     '  PROTOCOL(*),PSI(MAX_ORDER),SIGMA(MAX_ORDER+1),T,V(MAX_ORDER),
     '  W(MAX_ORDER),Y(NUM_VARS)
      LOGICAL CRASH,PHASE1,ROUND_CTRL,START
      EXTERNAL FUNC
!     Local Variables
      INTEGER i,IFAIL,iq,j,K_MINUS1,K_MINUS2,KNEW,K_PLUS1,K_PLUS2,l,
     '  NUMBER_STEPS_MINUS2,NUMBER_STEPS_PLUS1,NUMBER_STEPS_PLUS2
      REAL*8  ABSH,ERROR_K,ERROR_K_MINUS1,ERROR_K_MINUS2,
     '  ERROR_K_PLUS1,ERR,FOUR_EPSILON,GAMMASTAR(13),HALF_EPS,
     '  HNEW,R,RHO,ROUND,SUM,TAU,TEMP1,TEMP2,TOLD,POWERTWO(13),
     '  TWO_EPSILON
      LOGICAL GOODSTEP

C     Gamma*_(i) is defined to be:
C       Gamma*_(i)=1/i! \int{0}{1}[(s-1)(s)(s+1)...(s+i-2)]ds i=1,2,3,..
C     with Gamma*_(0)=1, but can be generated with the recursion
C     formula:
C       Gamma*_(m)+1/2.Gamma*_(m-1)+1/3.Gamma*_(m-2)+....+
C       1/(m+1).Gamma*_(0)=0
C     NOTE: For all the cases below Gamma* is negative however the
C     absolute value is stored as we are only interested in the
C     absolute error.
      REAL*8 GAMMASTAR_1,GAMMASTAR_2,GAMMASTAR_3,GAMMASTAR_4,
     '  GAMMASTAR_5,GAMMASTAR_6,GAMMASTAR_7,GAMMASTAR_8,GAMMASTAR_9,
     '  GAMMASTAR_10,GAMMASTAR_11,GAMMASTAR_12,GAMMASTAR_13
      PARAMETER(GAMMASTAR_1=1.0d0/2.0d0,GAMMASTAR_2=1.0d0/12.0d0,
     '  GAMMASTAR_3=1.0d0/24.0d0,GAMMASTAR_4=19.0d0/720.0d0,
     '  GAMMASTAR_5=27.0d0/1440.0d0,GAMMASTAR_6=863.0d0/60480.0d0,
     '  GAMMASTAR_7=275.0d0/24192.0d0,GAMMASTAR_8=33953.0d0/3628800.0d0,
     '  GAMMASTAR_9=8183.0d0/1036800.0d0,
     '  GAMMASTAR_10=3250433.0d0/479001600.0d0,
     '  GAMMASTAR_11=4671.0d0/788480.0d0,
     '  GAMMASTAR_12=5852897.0d0/1117670400.0d0,
     '  GAMMASTAR_13=78418523.0d0/16765056000.0d0)
      DATA GAMMASTAR / GAMMASTAR_1,GAMMASTAR_2,GAMMASTAR_3,GAMMASTAR_4,
     '  GAMMASTAR_5,GAMMASTAR_6,GAMMASTAR_7,GAMMASTAR_8,GAMMASTAR_9,
     '  GAMMASTAR_10,GAMMASTAR_11,GAMMASTAR_12,GAMMASTAR_13 /

      DATA POWERTWO /2.0d0,4.0d0,8.0d0,16.0d0,32.0d0,64.0d0,128.0d0,
     '  256.0d0,512.0d0,1024.0d0,2048.0d0,4096.0d0,8192.0d0/

      REAL*8 TIME(2)

      TWO_EPSILON=2.0d0*EPSILON
      FOUR_EPSILON=4.0d0*EPSILON

C     *** BLOCK 0 ***
C     Check if the step size or error tolerance is too small for the
C     machine precision. If it is the first step then initialise the
C     PHI array and estimate a starting step size.

      CRASH=.TRUE.
      IF(DABS(H).LT.FOUR_EPSILON*DABS(T)) THEN
C       If step size is too small the determine an acceptable one and
C       exit
        H=DSIGN(FOUR_EPSILON*DABS(T),H)
        IF(H.GT.MAX_STEP) H=MAX_STEP
      ELSE
        HALF_EPS=EPS/2.0d0
C
        ROUND=0.0d0
        DO l=1,NUMBER_EQN
          ROUND=ROUND+(Y(l)/ERROR_WEIGHT(l))**2
        ENDDO !l
        ROUND=TWO_EPSILON*DSQRT(ROUND)
        IF(HALF_EPS.LT.ROUND) THEN
C         If error tolerance is too small then increase it to an
C         acceptable value and exit
          EPS=2.0d0*ROUND*(1.0d0+FOUR_EPSILON)
        ELSE
          CRASH=.FALSE.
          G(1)=1.0d0
          G(2)=0.5d0
          SIGMA(1)=1.0d0

          IF(START) THEN
C           Initialise and compute appropriate step size for the first
C           step
            TIME(TCell) = T
            TIME(DTCell) = H
            CALL FUNC(TIME,Y,DY,CONTROL,MODEL,SIZES,VARIANT,DERIVED,
     '        PARAMETERS,PROTOCOL,AII,AIO,ARI,ARO,ERR_CODE)
            IF(ERR_CODE.EQ.0) THEN
              SUM=0.0d0
              DO l=1,NUMBER_EQN
                PHI(l,1)=DY(l)
                PHI(l,2)=0.0d0
                SUM=SUM+(DY(l)/ERROR_WEIGHT(l))**2
              ENDDO !l
              SUM=DSQRT(SUM)
              IF(EPS.LT.16.0d0*SUM*H*H) THEN
                ABSH=0.25d0*DSQRT(EPS/SUM)
              ELSE
                ABSH=DABS(H)
              ENDIF
              H=DSIGN(DMAX1(ABSH,FOUR_EPSILON*DABS(T)),H)
              IF(H.GT.MAX_STEP) H=MAX_STEP
              HOLD=0.0d0
              K=1
              KOLD=0
              START=.FALSE.
              PHASE1=.TRUE.
              IF(MAX_RND.GT.0) THEN
                ROUND_CTRL=.FALSE.
                IF(HALF_EPS.LE.100.0d0*ROUND) THEN
                  ROUND_CTRL=.TRUE.
                  DO l=1,NUMBER_EQN
                    PHI(l,MAX_ORDER+3)=0.0d0
                  ENDDO !l
                ENDIF
              ELSE
                ROUND_CTRL=.FALSE.
              ENDIF
            ENDIF
          ENDIF
C         *** End BLOCK 0 ***

          IF(ERR_CODE.EQ.0) THEN
            IFAIL=0

C           *** BLOCK 1 ***
C           Compute the coefficients of the formulas for this step.
C           Avoid computing those quantities not changed when the
C           step size is not changed.

            GOODSTEP=.FALSE.
            DO WHILE(.NOT.GOODSTEP.AND..NOT.CRASH.AND.ERR_CODE.EQ.0)

              K_PLUS1=K+1
              K_PLUS2=K+2
              K_MINUS1=K-1
              K_MINUS2=K-2

C             NUMBER_STEPS is the number of steps taken with step size
C             H. When the current order K is less than NUMBER_STEPS
C             then no coefficients change.

              IF(H.NE.HOLD) NUMBER_STEPS = 0
              IF(NUMBER_STEPS.LE.KOLD) NUMBER_STEPS=NUMBER_STEPS+1
              NUMBER_STEPS_PLUS1 = NUMBER_STEPS+1
              IF (K.GE.NUMBER_STEPS) THEN
C               Compute those components of ALPHA, BETA, PSI and SIGMA
C               which are changed
C               The formula are:
C                 Psi_i(n+1)=h_(n+1)+h_n+...+h_(n+2-i)      i>=1
C                 Alpha_i(n+1)=h_(n+1)/Psi_i(n+1)           i>=1
C                 Beta_1(n+1)=1.0                           i=1
C                 Beta_i(n+1)=Psi_1(n+1).Psi_2(n+1)...Psi_(i-1)(n+1)/
C                   Psi_1(n).Psi_2(n)...Psi_(i-1)(n)        i>1
C                 Sigma_1(n+1)=1.0                          i=1
C                 Sigma_i(n+1)=h.2h...(i-1)h/
C                   Psi_1(n+1).Psi_2(n+1)...Psi_(i-1)(n+1)  i=>1
                BETA(NUMBER_STEPS)=1.0d0
                ALPHA(NUMBER_STEPS)=1.0d0/DBLE(NUMBER_STEPS)
                TEMP1=H*DBLE(NUMBER_STEPS)
                SIGMA(NUMBER_STEPS_PLUS1)=1.0d0
                DO i=NUMBER_STEPS_PLUS1,K
                  TEMP2=PSI(i-1)
                  PSI(i-1)=TEMP1
                  BETA(i)=BETA(i-1)*PSI(i-1)/TEMP2
                  TEMP1=TEMP2+H
                  ALPHA(i)=H/TEMP1
                  SIGMA(i+1)=DBLE(i)*ALPHA(i)*SIGMA(i)
                ENDDO !i
                PSI(K)=TEMP1

C               Compute the coefficients G
                IF(NUMBER_STEPS.LE.1) THEN
C                 Initialise V and set W
                  DO iq=1,K
                    V(iq)=1.0d0/DBLE(iq*(iq+1))
                    W(iq)=V(iq)
                  ENDDO !iq
                ELSE
                  IF(K.GT.KOLD) THEN
C                   If the order was raised update the diagonal part
C                   of V
                    V(K)=1.0d0/DBLE(K*K_PLUS1)
                    NUMBER_STEPS_MINUS2=NUMBER_STEPS-2
                    DO j=1,NUMBER_STEPS_MINUS2
                      i=K-j
                      V(i)=V(i)-ALPHA(j+1)*V(i+1)
                    ENDDO !J
                  ENDIF
C                 Update V and set W
                  DO iq=1,K_PLUS1-NUMBER_STEPS
                    V(iq)=V(iq)-ALPHA(NUMBER_STEPS)*V(iq+1)
                    W(iq)=V(iq)
                  ENDDO !IQ
                  G(NUMBER_STEPS_PLUS1)=W(1)
C
                ENDIF
C               Compute the G in the work vector W
                NUMBER_STEPS_PLUS2 = NUMBER_STEPS + 2
                DO i=NUMBER_STEPS_PLUS2,K_PLUS1
                  DO iq=1,K_PLUS2-i
                    W(iq)=W(iq)-ALPHA(i-1)*W(iq+1)
                  ENDDO !iq
                  G(i)=W(1)
                ENDDO !i

              ENDIF
C             *** End BLOCK 1 ***

C             *** BLOCK 2 ***
C             Predicit a solution, PREDICTION, and evaluate the
C             derivatives, DY, using the predicited solution.
C             Estimate the local error at order K and the errors
C             at orders K, K-1 and K-2 as if a constant step size
C             were used.

C             Change PHI to PHI* i.e. PHI*_i(n)=BETA_i(n+1).PHI_i(n)
              DO i=NUMBER_STEPS_PLUS1,K
                DO l=1,NUMBER_EQN
                  PHI(l,i)=BETA(i)*PHI(l,i)
                ENDDO !l
              ENDDO !i

C             Predicit the solution and differences
              DO l=1,NUMBER_EQN
                PHI(l,K_PLUS2)=PHI(l,K_PLUS1)
                PHI(l,K_PLUS1)=0.0d0
                PREDICTION(l)=0.0d0
              ENDDO !l
              !initialise the non-ODE variables
              DO l=NUMBER_EQN+1,NUM_VARS
                PREDICTION(l)=Y(l)
              ENDDO !l
              DO j=1,K
                i=K_PLUS1-j
                DO l=1,NUMBER_EQN
                  PREDICTION(l)=PREDICTION(l)+G(i)*PHI(l,i)
                  PHI(l,i)=PHI(l,i)+PHI(l,i+1)
                ENDDO !l
              ENDDO !i
              IF(ROUND_CTRL) THEN
                DO l=1,NUMBER_EQN
                  TAU=H*PREDICTION(l)-PHI(l,MAX_ORDER+3)
                  PREDICTION(l)=Y(l)+TAU
                  PHI(l,MAX_ORDER+4)=(PREDICTION(l)-Y(l))-TAU
                ENDDO !l
              ELSE
                DO l=1,NUMBER_EQN
                  PREDICTION(l)=Y(l)+H*PREDICTION(l)
                ENDDO !l
              ENDIF
              TOLD=T
              T=T+H
              ABSH=DABS(H)
              TIME(TCell) = T
              TIME(DTCell) = H
              CALL FUNC(TIME,PREDICTION,DY,CONTROL,MODEL,SIZES,VARIANT,
     '          DERIVED,PARAMETERS,PROTOCOL,AII,AIO,ARI,ARO,ERR_CODE)
              IF(ERR_CODE.EQ.0) THEN
C               Estimate errors at orders K,K-1,K-2
C               For a constant step size, h, the local error at x_(n+1)
C               is given by:
C                 ERROR_K=|h.Gamma*_K.Sigma_(K+1)(n+1).Phi_(K+1)^P(n+1)|
                ERROR_K=0.0d0
                DO l=1,NUMBER_EQN
                  ERROR_K=ERROR_K+
     '              ((DY(l)-PHI(l,1))/ERROR_WEIGHT(l))**2
                ENDDO !l
                ERR=ABSH*DSQRT(ERROR_K)*(G(K)-G(K_PLUS1))
                ERROR_K=ABSH*DSQRT(ERROR_K)*SIGMA(K_PLUS1)*
     '            GAMMASTAR(K)
                IF(K_MINUS2.EQ.0) THEN
                  ERROR_K_MINUS1=0.0d0
                  DO l=1,NUMBER_EQN
                    ERROR_K_MINUS1=ERROR_K_MINUS1+
     '                ((DY(l)+PHI(L,K)-PHI(l,1))/ERROR_WEIGHT(l))**2
                  ENDDO !l
                  ERROR_K_MINUS1=ABSH*SIGMA(K)*GAMMASTAR(K_MINUS1)*
     '              DSQRT(ERROR_K_MINUS1)
                ELSE IF(K_MINUS2.GT.0) THEN
                  ERROR_K_MINUS2=0.0d0
                  ERROR_K_MINUS1=0.0d0
                  DO l=1,NUMBER_EQN
                    ERROR_K_MINUS2=ERROR_K_MINUS2+
     '                ((DY(l)+PHI(L,K_MINUS1)-PHI(l,1))/
     '                ERROR_WEIGHT(l))**2
                    ERROR_K_MINUS1=ERROR_K_MINUS1+
     '                ((DY(l)+PHI(L,K)-PHI(l,1))/ERROR_WEIGHT(l))**2
                  ENDDO !l
                  ERROR_K_MINUS2=ABSH*SIGMA(K_MINUS1)*
     '              GAMMASTAR(K_MINUS2)*DSQRT(ERROR_K_MINUS2)
                  ERROR_K_MINUS1=ABSH*SIGMA(K)*GAMMASTAR(K_MINUS1)*
     '              DSQRT(ERROR_K_MINUS1)
                ENDIF
                KNEW=K

C               Test if order should be lowered
                IF(K_MINUS2.EQ.0) THEN
                  IF(ERROR_K_MINUS1.LE.0.5D0*ERROR_K) KNEW = K_MINUS1
                ELSE IF(K_MINUS2.GT.0) THEN
                  IF(DMAX1(ERROR_K_MINUS1,ERROR_K_MINUS2).LE.ERROR_K)
     '              KNEW = K_MINUS1
                ENDIF
C               *** End BLOCK 2 ***

C               Test if the step was successful
                IF(ERR.GT.EPS) THEN

C                 *** BLOCK 3 ***
C                 The step is unsuccessful.

C                 Restore T, PHI and PSI
                  PHASE1=.FALSE.
                  T=TOLD
                  DO i=1,K
                    DO l=1,NUMBER_EQN
                      PHI(l,i)=(PHI(l,i)-PHI(l,i+1))/BETA(i)
                    ENDDO !l
                  ENDDO !i
                  IF(K.GE.2) THEN
                    DO i=2,K
                      PSI(i-1)=PSI(i)-H
                    ENDDO !i
                  ENDIF

                  IFAIL=IFAIL+1
C                 Double the step size. If the step fails three times
C                 set the order to 1, thereafter use optimal step size.
C                 This procedure will exit if the estimated step size
C                 is too small for the machine precision.

                  IF(IFAIL.LT.3) THEN
                    H=H/2.0d0
                  ELSE IF(IFAIL.EQ.3) THEN
                    KNEW = 1
                    H=H/2.0d0
                  ELSE
                    KNEW = 1
                    IF(HALF_EPS.LT.ERROR_K/4.0d0) THEN
                      H=DSQRT(HALF_EPS/ERROR_K)*H
                    ELSE
                      H=H/2.0d0
                    ENDIF
                  ENDIF
                  K = KNEW
                  IF(DABS(H).LT.FOUR_EPSILON*DABS(T)) THEN
                    CRASH =.TRUE.
                    H=DSIGN(FOUR_EPSILON*DABS(T),H)
                    IF(H.GT.MAX_STEP) H=MAX_STEP
                    EPS=EPS+EPS
                  ENDIF
C                 *** End BLOCK 3 ***

                ELSE
                  GOODSTEP=.TRUE.
                ENDIF
              ENDIF

            ENDDO

            IF(.NOT.CRASH.AND.ERR_CODE.EQ.0) THEN

C             *** BLOCK 4 ***
C             The step is successfull. Correct the predicited solution,
C             evaluate the derivatives using the correct solution and
C             update the differences. Determine best order and step
C             size for the next step.

              KOLD=K
              HOLD=H
C             Correct and Evaluate
              IF(ROUND_CTRL) THEN
                DO l=1,NUMBER_EQN
                  RHO=H*G(K_PLUS1)*(DY(l)-PHI(l,1))-PHI(l,MAX_ORDER+4)
                  Y(l)=PREDICTION(l)+RHO
                  PHI(l,MAX_ORDER+3)=(Y(l)-PREDICTION(l))-RHO
                ENDDO !l
              ELSE
                DO l=1,NUMBER_EQN
                  Y(l)=PREDICTION(l)+H*G(K_PLUS1)*(DY(l)-PHI(l,1))
                ENDDO !l
              ENDIF
              TIME(TCell) = T
              TIME(DTCell) = H
              CALL FUNC(TIME,Y,DY,CONTROL,MODEL,SIZES,VARIANT,DERIVED,
     '          PARAMETERS,PROTOCOL,AII,AIO,ARI,ARO,ERR_CODE)
              IF(ERR_CODE.EQ.0) THEN
C               Update differences for next step
                DO l=1,NUMBER_EQN
                  PHI(l,K_PLUS1)=DY(l)-PHI(l,1)
                  PHI(l,K_PLUS2)=PHI(l,K_PLUS1)-PHI(l,K_PLUS2)
                ENDDO !l
                DO i=1,K
                  DO l=1,NUMBER_EQN
                    PHI(l,i)=PHI(l,i)+PHI(l,K_PLUS1)
                  ENDDO !l
                ENDDO !k

C               Estimate the error at order K+1 unless: we are in
C               the first phase in which case always raise the order;
C               we have already decieded to lower the order; or the
C               step size is not constant so the error estimate is
C               unreliable.
                ERROR_K_PLUS1=0.0d0
                IF(KNEW.EQ.K_MINUS1.OR.K.EQ.MAX_ORDER) PHASE1=.FALSE.
                IF(PHASE1) THEN
C                 Raise the order
                  K=K_PLUS1
                  ERROR_K=ERROR_K_PLUS1
                ELSE
                  IF(KNEW.EQ.K_MINUS1) THEN
C                   Lower the order
                    K=K_MINUS1
                    ERROR_K=ERROR_K_MINUS1
                  ELSE IF(K_PLUS1.LE.NUMBER_STEPS) THEN
                    DO l=1,NUMBER_EQN
                      ERROR_K_PLUS1=ERROR_K_PLUS1+
     '                  (PHI(l,K_PLUS2)/ERROR_WEIGHT(l))**2
                    ENDDO !l
                    ERROR_K_PLUS1=ABSH*GAMMASTAR(K_PLUS1)*
     '                DSQRT(ERROR_K_PLUS1)
C                   Using the estimated error at order K+1 determine
C                   the appropriate order for the next step.
                    IF(K.LE.1) THEN
                      IF(ERROR_K_PLUS1.LT.ERROR_K/2.0d0.AND.
     '                  K.NE.MAX_ORDER) THEN
C                       Raise the order
                        K=K_PLUS1
                        ERROR_K=ERROR_K_PLUS1
                      ENDIF
                    ELSE
                      IF(ERROR_K_MINUS1.LE.
     '                  DMIN1(ERROR_K,ERROR_K_PLUS1)) THEN
C                       Lower the order
                        K=K_MINUS1
                        ERROR_K=ERROR_K_MINUS1
                      ELSE IF(ERROR_K_PLUS1.LT.ERROR_K.AND.
     '                    K.NE.MAX_ORDER) THEN
C                       Here the error at K+1 < the error at K < the
C                       maximum error of K-1 and K-2 otherwise the
C                       order would have been lowered in block 2, thus
C                       raise the order
                        K = K_PLUS1
                        ERROR_K = ERROR_K_PLUS1
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF

C               With the new order determine the appropriate step
C               size for the next step.
                HNEW=H+H
                IF(.NOT.PHASE1.AND.
     '            HALF_EPS.LT.ERROR_K*POWERTWO(K+1)) THEN
                  HNEW=H
                  IF(HALF_EPS.LT.ERROR_K) THEN
                    R=(HALF_EPS/ERROR_K)**(1.0d0/DBLE(K+1))
                    HNEW=ABSH*DMAX1(0.5d0,DMIN1(0.9d0,R))
                    HNEW=DSIGN(DMAX1(HNEW,FOUR_EPSILON*DABS(T)),H)
                  ENDIF
                ENDIF
                H=HNEW
                IF(H.GT.MAX_STEP) H=MAX_STEP
              ENDIF
C             *** End BLOCK 4 ***

            ENDIF
          ENDIF
        ENDIF
      ENDIF

      RETURN
      END

