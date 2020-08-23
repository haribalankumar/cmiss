      SUBROUTINE AM_DE(AII,AIO,CONTROL,ERROR_TYPE,IFAIL,ISNOLD,K,
     '  KOLD,MAX_ITERS,MAX_ORDER,MAX_RND,MODEL,NUMBER_EQN,NUM_VARS,
     '  NUMBER_STEPS,SIZES,VARIANT,ABS_ERR,ALPHA,ARI,ARO,BETA,DELSGN,
     '  DERIVED,DY,DYOUT,ERROR_WEIGHT,G,H,HOLD,MAX_STEP,PARAMETERS,
     '  PHI,PREDICTION,PROTOCOL,PSI,REL_ERR,SIGMA,T,TT,TOLD,TOUT,V,W,Y,
     '  YY,EXTEND_INTERVAL,PHASE1,ROUND_CTRL,START,STIFF_EQNS,FUNC)

C#### Subroutine: AM_DE
C###  Description:
C###    Adams-Moulton differential equation integrator. Called from
C###    the ADAMS buffer subroutine.

      IMPLICIT NONE

      INCLUDE 'cell_reserved.inc'

!     Parameter list
      INTEGER AII(*),AIO(*),CONTROL(*),ERROR_TYPE,IFAIL,ISNOLD,K,KOLD,
     '  MAX_ITERS,MAX_ORDER,MAX_RND,MODEL(*),NUMBER_EQN,NUM_VARS,
     '  NUMBER_STEPS,SIZES(*),VARIANT
      REAL*8 ABS_ERR,ALPHA(MAX_ORDER),ARI(*),ARO(*),BETA(MAX_ORDER),
     '  DELSGN,DERIVED(*),DY(NUM_VARS),DYOUT(NUM_VARS),
     '  ERROR_WEIGHT(NUM_VARS),G(MAX_ORDER+1),H,HOLD,MAX_STEP,
     '  PARAMETERS(*),PHI(NUM_VARS,MAX_ORDER+2+MAX_RND),
     '  PREDICTION(NUM_VARS),PROTOCOL(*),PSI(MAX_ORDER),REL_ERR,
     '  SIGMA(MAX_ORDER+1),T,TT,TOLD,TOUT,V(MAX_ORDER),W(MAX_ORDER),
     '  Y(NUM_VARS),YY(NUM_VARS)
      LOGICAL EXTEND_INTERVAL,PHASE1,ROUND_CTRL,START,STIFF_EQNS
      EXTERNAL FUNC
!     Local Variables
      INTEGER ERR_CODE,ISN,NUM_TIMES_LOW_ORDER,L,IT_NUM
      REAL*8 ABS_DEL,ABS_EPS,DEL,DLAMCH,EPS,EPSILON,FOUR_EPSILON,
     '  REL_EPS,TEND,TIME(2)
      LOGICAL CRASH,FINISHED

      EPSILON=DLAMCH('EPS')
      FOUR_EPSILON=4.0d0*EPSILON

C     Test for improper parameters
      EPS=DMAX1(REL_ERR,ABS_ERR)
      IFAIL=IABS(IFAIL)
      IF(NUMBER_EQN.LT.1.OR.
     '  T.EQ.TOUT.OR.
     '  REL_ERR.LT.0.0d0.OR.ABS_ERR.LT.0.0d0.OR.EPS.LE.0.0d0.OR.
     '  (IFAIL.LT.1.OR.IFAIL.GT.5).OR.
     '  (IFAIL.NE.1.AND.T.NE.TOLD).OR.
     '  (ERROR_TYPE.LT.1.OR.ERROR_TYPE.GT.4).OR.
     '  (MAX_ORDER.LT.1.OR.MAX_ORDER.GT.12)) THEN
        IFAIL = 6
      ELSE

C       On each call set interval of integration and counter for the
C       number of steps.
        DEL=TOUT-T
        ABS_DEL=DABS(DEL)
        TEND=T+10.0d0*DEL
        ISN=ISIGN(1,IFAIL)
        IF(ISN.LT.0) TEND = TOUT
        IT_NUM=0
        NUM_TIMES_LOW_ORDER=0
        STIFF_EQNS=.FALSE.
        REL_EPS=REL_ERR/EPS
        ABS_EPS=ABS_ERR/EPS
        IF(IFAIL.EQ.1.OR.ISNOLD.LT.0.OR.DELSGN*DEL.LE.0.0d0) THEN
C         On start and restart also set the work variables TT and YY,
C         store the direction of integration and initialise the step
C         size
          START=.TRUE.
          TT=T
          DO l=1,NUM_VARS
            YY(l)=Y(l)
          ENDDO !l
          DELSGN=DSIGN(1.0d0,DEL)
          H=DSIGN(DMAX1(DABS(TOUT-TT),FOUR_EPSILON*DABS(TT)),TOUT-TT)
          IF(H.GT.MAX_STEP) H=MAX_STEP
        ENDIF

C       Iterate on the solution until finished.
        FINISHED=.FALSE.
        DO WHILE(.NOT.FINISHED)

          IF(DABS(TT-T).GE.ABS_DEL) THEN
C           If already past the output point then interpolate and
C           finish
            CALL AM_INTERPOLATE(KOLD,MAX_ORDER,MAX_RND,
     '        NUMBER_EQN,NUM_VARS,DYOUT,PHI,PSI,TT,TOUT,YY,Y)
            IFAIL=2
            T=TOUT
            TOLD=T
            ISNOLD=ISN
            FINISHED=.TRUE.
          ELSE IF(.NOT.EXTEND_INTERVAL.AND.
     '        DABS(TOUT-TT).LT.FOUR_EPSILON*DABS(TT)) THEN
C           If you cannot extend the interval to go past the output
C           output point and you are sufficiently close then
C           extroplate and finish.
            H=TOUT-TT
            TIME(TCell) = TT
            TIME(DTCell) = H
            CALL FUNC(TIME,YY,DY,CONTROL,MODEL,SIZES,VARIANT,DERIVED,
     '        PARAMETERS,PROTOCOL,AII,AIO,ARI,ARO,ERR_CODE)
            IF(ERR_CODE.EQ.0) THEN
              DO l=1,NUMBER_EQN
                Y(l)=YY(l)+H*DY(l)
              ENDDO !l
              !don't need to update non-ODE variables
              IFAIL=2
              T=TOUT
              TOLD=T
              ISNOLD=ISN
              FINISHED=.TRUE.
            ELSE
              IFAIL=7
              FINISHED=.TRUE.
            ENDIF
          ELSE IF(IT_NUM.GE.MAX_ITERS) THEN
C           Test for too many steps
            IFAIL=ISN*4
            IF(STIFF_EQNS) IFAIL=ISN*5
            DO l=1,NUM_VARS
              Y(l)=YY(l)
            ENDDO !l
            T=TT
            TOLD=T
            ISNOLD=1
            FINISHED=.TRUE.
          ELSE
C           Limit step size, set the error weight vector and take
C           a step
            H=DSIGN(DMIN1(DABS(H),DABS(TEND-TT)),H)
            IF(H.GT.MAX_STEP) H=MAX_STEP
            ERR_CODE=0
            IF(ERROR_TYPE.EQ.1) THEN
              DO l=1,NUMBER_EQN
                ERROR_WEIGHT(l)=1.0d0
              ENDDO !l
            ELSE IF(ERROR_TYPE.EQ.2) THEN
              DO l=1,NUMBER_EQN
                ERROR_WEIGHT(l)=DABS(YY(l))
                IF(DABS(YY(l)).LE.EPSILON) ERR_CODE=1
              ENDDO !l
            ELSE IF(ERROR_TYPE.EQ.3) THEN
              DO l=1,NUMBER_EQN
                ERROR_WEIGHT(l)=DABS(DY(l))
                IF(DABS(DY(l)).LE.EPSILON) ERR_CODE=1
              ENDDO !l
            ELSE IF(ERROR_TYPE.EQ.4) THEN
              DO l=1,NUMBER_EQN
                ERROR_WEIGHT(l)=REL_EPS*DABS(YY(l))+ABS_EPS
              ENDDO !l
            ENDIF
            IF(ERR_CODE.EQ.0) THEN

              CALL AM_STEP(AII,AIO,CONTROL,ERR_CODE,K,KOLD,
     '          MAX_ORDER,MAX_RND,MODEL,NUMBER_EQN,NUM_VARS,
     '          NUMBER_STEPS,SIZES,VARIANT,ALPHA,ARI,ARO,BETA,DERIVED,
     '          DY,EPS,EPSILON,ERROR_WEIGHT,G,H,HOLD,MAX_STEP,
     '          PARAMETERS,PHI,PREDICTION,PROTOCOL,PSI,SIGMA,TT,V,W,YY,
     '          CRASH,PHASE1,ROUND_CTRL,START,FUNC)

              IF(ERR_CODE.EQ.0) THEN
                IF(CRASH) THEN
C                 Tolerances too small
                  IFAIL=ISN*3
                  REL_ERR=EPS*REL_EPS
                  ABS_ERR=EPS*ABS_EPS
                  DO l=1,NUM_VARS
                    Y(l)=YY(l)
                  ENDDO !l
                  T=TT
                  TOLD=T
                  ISNOLD=1
                  FINISHED=.TRUE.
                ELSE
C                 Adjust number of steps and test for stiffness
                  IT_NUM=IT_NUM+1
                  NUM_TIMES_LOW_ORDER=NUM_TIMES_LOW_ORDER+1
                  IF(KOLD.GT.4) NUM_TIMES_LOW_ORDER=0
                  IF(NUM_TIMES_LOW_ORDER.GE.50) STIFF_EQNS=.TRUE.
                ENDIF
              ELSE
C               Error from func
                IFAIL=7
                FINISHED=.TRUE.
              ENDIF
            ELSE
C             Error weight is zero
              IFAIL=8
              FINISHED=.TRUE.
            ENDIF
          ENDIF
        ENDDO
      ENDIF

      RETURN
      END


