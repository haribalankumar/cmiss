      SUBROUTINE ADAMS(AII,AIO,CONTROL,ERROR_TYPE,IFAIL,
     '  IWORK,L_IWORK,L_WORK,MAX_ITERS,MAX_ORDER,MODEL,NUMBER_EQN,
     '  NUM_VARS,SIZES,VARIANT,ABS_ERR,ARI,ARO,DERIVED,DY,MAX_STEP,
     '  PARAMETERS,PROTOCOL,REL_ERR,T,TOUT,WORK,Y,EXTEND_INTERVAL,
     '  STIFF_EQNS,USE_ROUND_CTRL,FUNC,ERROR)

C#### Subroutine: ADAMS
C###  Description:
C###    <HTML><PRE>
C###    This procedure integrates a system of NUMBER_EQN ordinary
C###    differential equations of the form
C###      dy(i)/dt=FUNC(t,y(1),y(2),...,y(NUMBER_EQN))
C###    from t=T to t=TOUT given y(i) at t=T using an Adams-Moulton
C###    integration scheme. The integrator used is adapted from the
C###    Adams-Moulton integrator of L.F. Shampine and M.K. Gordon
C###    and is detailed in their book "Computer Solution of Ordinary
C###    Differential Equations: The Initial Value Problem".
C###    The exact form of the rhs function is
C###      CALL FUNC(T,Y,DY,CONTROL,MODEL,SIZES,VARIANT,DERIVED,
C###        PARAMETERS,PROTOCOL,AII,AIO,ARI,ARO,ERR)
C###    where CONTROL, MODEL, SIZES ARI, ARO are integer vectors
C###    and DERIVED, PARAMETERS, PROTOCOL, Y, DY, ARI, ARO are
C###    double precision vectors.
C###    The integrator uses a modified divided difference form of the
C###    Adams PECE formulas and local extropolation. It iterates on
C###    the solution (up to MAX_ITERS) and adjusts the order (up to
C###    MAX_ORDER and no more than 12) and step size (up to MAX_STEP)
C###    to control the local error per unit step in a generalised
C###    sense.
C###    The error control is based on the L2 Norm of the weighted
C###    local error vector. The weighting is controlled by ERROR_TYPE.
C###    If ERROR_TYPE=1 then no weighting is used; if ERROR_TYPE=2
C###    then the error vector components are weighted by the most
C###    recent component of the solution vector; if ERROR_TYPE=3 then
C###    the error vector components are weighted by the most recent
C###    component of the residual (derivative) vector; if ERROR_TYPE=4
C###    then a mixed relative weighting is used. This weighting is
C###    calculated from (REL_ERR*y(i)+ABS_ERR)/MAX(REL_ERR,ABS_ERR).
C###    For reasons of efficiency the integrator integrates beyond
C###    TOUT internally (though never beyond T+10*(TOUT-T)) and
C###    interpolates the solution at TOUT. If it is not possible to
C###    integrate beyond TOUT then EXTEND_INTERVAL should be set to
C###    .FALSE.
C###    The integrator can perform propagative rounding error control
C###    by setting USE_ROUND_CONTROL to .TRUE.
C###    The integrator needs workspace allocated as follows:
C###      IWORK(L_IWORK) where L_IWORK=7 and
C###      WORK(L_WORK) where L_WORK=7+7*MAX_ORDER+
C###        NUMBER_EQN*(MAX_ORDER+6) in general, and
C###        L_WORK=7+7*MAX_ORDER+NUMBER_EQN*(MAX_ORDER+8) if rounding
C###        control is required.
C###    The workspace needs to be mainted between calls. Interesting
C###    values contained within the workspace are:
C###      IWORK(1) - the order, K, of the Adams polynomial to be used
C###        with the next call,
C###      IWORK(2) - the order of the Adams polynomial used in the
C###        previous (current) step,
C###      WORK(1) - the step size, H, to be used with the next call,
C###      WORK(2) - the step size used in the previous (current) step.
C###      WORK(3) - the finish time of the last step taken.
C###      WORK(4) - the finish time of the last step output.
C###    On entry the following parameters are required:
C###      NUMBER_EQN - the number of equations to be integrated,
C###      NUM_VARS - the total number of state variables
C###      Y - vector of initial conditions,
C###      T - starting point of the integration,
C###      TOUT - point at which the solution is desired,
C###      ABS_ERR, REL_ERR - absolute and relative local error
C###        tolerances,
C###      ERROR_TYPE - the type of error test to be used,
C###      IFAIL - error code indicator. Must be set to 1 for start-up,
C###      MAX_ITERS - the maximum number of Adams iterations (steps)
C###        allowed to reach TOUT,
C###      MAX_ORDER - the maximum order allowed for the Adams
C###        polynomial,
C###      MAX_STEP - the maximum step size allowed to reach TOUT,
C###      EXTEND_INTERVAL - set to .TRUE. if it is possible to
C###        integrate beyond TOUT.
C###      USE_ROUND_CTRL - set to .TRUE. if rounding control is to be
C###        used.
C###    On output the following parameters are set:
C###      Y - vector of solutions at TOUT,
C###      DY - vector of the derivates at TOUT,
C###      T - last point reach in the integration. Normal return has
C###        T=TOUT,
C###      ABS_ERR, REL_ERR - normal return has the tolerances unchanged.
C###        If IFAIL=3 then the tolerances are set to increased values.
C###      STIFF_EQNS - set to .TRUE. if the equations appear to be
C###        stiff,
C###      IFAIL - error flag. If
C###        IFAIL=2 - normal return, integration reached TOUT,
C###        IFAIL=3 - integration did not reach TOUT because the error
C###          are too small. ABS_ERR and REL_ERR have been increased
C###          to appropriate levels.
C###        IFAIL=4 - integration did not reach TOUT because more than
C###          MAX_ITERS steps were needed,
C###        IFAIL=5 - integration did not reach TOUT because equations
C###          appear to be stiff,
C###        IFAIL=6 - invalid input parameters,
C###        IFAIL=7 - error returned from FUNC,
C###        IFAIL=8 - if ERROR_TYPE is 2 or 3 then one of error weights
C###          is zero,
C###        IFAIL=9 - WORK array is too small,
C###        IFAIL=10 - IWORK array is too small.
C###    </PRE></HTML>

      IMPLICIT NONE

!     Parameter list
      INTEGER L_IWORK,L_WORK,AII(*),AIO(*),CONTROL(*),ERROR_TYPE,
     '  IFAIL,IWORK(L_IWORK),MAX_ITERS,MAX_ORDER,MODEL(*),
     '  NUMBER_EQN,NUM_VARS,SIZES(*),VARIANT
      REAL*8 ABS_ERR,ARI(*),ARO(*),DERIVED(*),DY(*),
     '  MAX_STEP,PARAMETERS(*),PROTOCOL(*),REL_ERR,T,TOUT,
     '  WORK(L_WORK),Y(*)
      LOGICAL EXTEND_INTERVAL,STIFF_EQNS,USE_ROUND_CTRL
      EXTERNAL FUNC
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER MAX_RND
      INTEGER ISNOLD_IDX,K_IDX,KOLD_IDX,NUMBER_STEPS_IDX,
     '  PHASE1_IDX,REQUIRED_IDX,ROUND_CTRL_IDX,START_IDX
      INTEGER H_IDX,HOLD_IDX,TT_IDX,TOLD_IDX,ALPHA_IDX,BETA_IDX,
     '  DELSGN_IDX,DY_IDX,ERROR_WEIGHT_IDX,G_IDX,PREDICTION_IDX,
     '  PSI_IDX,SIGMA_IDX,V_IDX,W_IDX,YY_IDX,PHI_IDX
      LOGICAL ROUND_CTRL,PHASE1,START

      DATA ISNOLD_IDX,K_IDX,KOLD_IDX,NUMBER_STEPS_IDX,PHASE1_IDX,
     '  ROUND_CTRL_IDX,START_IDX / 3,1,2,4,5,6,7 /
      DATA H_IDX,HOLD_IDX,TT_IDX,TOLD_IDX,DELSGN_IDX,
     '  ALPHA_IDX / 1,2,3,4,5,6 /


      IF(L_IWORK.GE.7) THEN
        REQUIRED_IDX=7+7*MAX_ORDER+NUM_VARS*(MAX_ORDER+6)
        IF(USE_ROUND_CTRL) REQUIRED_IDX=REQUIRED_IDX+2*NUM_VARS
        IF(L_WORK.GE.REQUIRED_IDX) THEN
          BETA_IDX=ALPHA_IDX+MAX_ORDER
          DY_IDX=BETA_IDX+MAX_ORDER
          ERROR_WEIGHT_IDX=DY_IDX+NUM_VARS
          G_IDX=ERROR_WEIGHT_IDX+NUM_VARS
          PSI_IDX=G_IDX+MAX_ORDER+1
          PREDICTION_IDX=PSI_IDX+MAX_ORDER
          SIGMA_IDX=PREDICTION_IDX+NUM_VARS
          V_IDX=SIGMA_IDX+MAX_ORDER+1
          W_IDX=V_IDX+MAX_ORDER
          YY_IDX=W_IDX+MAX_ORDER
          PHI_IDX=YY_IDX+NUM_VARS

          IF(USE_ROUND_CTRL) THEN
            MAX_RND=2
          ELSE
            MAX_RND=0
          ENDIF

          IF(ABS(IFAIL).NE.1) THEN
            START=IWORK(START_IDX).NE.0
            PHASE1=IWORK(PHASE1_IDX).NE.0
            IF(USE_ROUND_CTRL) THEN
              ROUND_CTRL=IWORK(ROUND_CTRL_IDX).NE.0
            ELSE
              ROUND_CTRL=.FALSE.
            ENDIF
          ENDIF

          CALL AM_DE(AII,AIO,CONTROL,ERROR_TYPE,IFAIL,
     +      IWORK(ISNOLD_IDX),IWORK(K_IDX),IWORK(KOLD_IDX),
     +      MAX_ITERS,MAX_ORDER,MAX_RND,MODEL,NUMBER_EQN,NUM_VARS,
     +      IWORK(NUMBER_STEPS_IDX),SIZES,VARIANT,ABS_ERR,
     +      WORK(ALPHA_IDX),ARI,ARO,WORK(BETA_IDX),
     +      WORK(DELSGN_IDX),DERIVED,WORK(DY_IDX),DY,
     +      WORK(ERROR_WEIGHT_IDX),WORK(G_IDX),WORK(H_IDX),
     +      WORK(HOLD_IDX),MAX_STEP,PARAMETERS,WORK(PHI_IDX),
     +      WORK(PREDICTION_IDX),PROTOCOL,WORK(PSI_IDX),REL_ERR,
     +      WORK(SIGMA_IDX),T,WORK(TT_IDX),WORK(TOLD_IDX),
     +      TOUT,WORK(V_IDX),WORK(W_IDX),Y,WORK(YY_IDX),
     +      EXTEND_INTERVAL,PHASE1,ROUND_CTRL,START,STIFF_EQNS,FUNC)

          IF(START) THEN
            IWORK(START_IDX)=1
          ELSE
            IWORK(START_IDX)=0
          ENDIF
          IF(PHASE1) THEN
            IWORK(PHASE1_IDX)=1
          ELSE
            IWORK(PHASE1_IDX)=0
          ENDIF
          IF(USE_ROUND_CTRL) THEN
            IF(ROUND_CTRL) THEN
              IWORK(ROUND_CTRL_IDX)=1
            ELSE
              IWORK(ROUND_CTRL_IDX)=0
            ENDIF
          ENDIF

          IF(ABS(IFAIL).EQ.3) THEN
            ERROR='>>Integration did not reach TOUT, Tolerances too '
     '        //'small'
          ELSE IF(ABS(IFAIL).EQ.4) THEN
            ERROR='>>Integration did not reach TOUT, Maximum '
     '        //'iterations exceeded'
          ELSE IF(ABS(IFAIL).EQ.5) THEN
            ERROR='>>Integration did not reach TOUT, Equations are '
     '        //'too stiff'
          ELSE IF(ABS(IFAIL).EQ.6) THEN
            ERROR='>>Input parameters invalid'
          ELSE IF(ABS(IFAIL).EQ.7) THEN
            ERROR='>>Error returned from FUNC'
          ELSE IF(ABS(IFAIL).EQ.8) THEN
            ERROR='>>Invalid zero error weight to divide'
          ENDIF
        ELSE
          WRITE(ERROR,'(''>>WORK array is too small, increase L_WORK '
     '      //'to '',I5)') REQUIRED_IDX
          IFAIL=9
        ENDIF
      ELSE
        ERROR='>>IWORK array is too small, increase L_IWORK to 7'
        IFAIL=10
      ENDIF

      RETURN
      END

