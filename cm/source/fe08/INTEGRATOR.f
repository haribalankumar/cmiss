      SUBROUTINE INTEGRATOR(INTEGRATOR_IWORK,AII,AIO,CONTROL,maqdt,
     '  MODEL,
     '  NUM_ODES,NUM_VARS,nr,nx,SIZES,VARIANT,INTEGRATOR_WORK,AQ,ARI,
     '  ARO,
     '  DERIVED,PARAMETERS,PROTOCOL,STIMULUS,T,YQS,RHSROUTINE,ERROR,*)

C#### Subroutine: INTEGRATOR
C###  Description:
C###    <html>
C###    The INTEGRATOR integrates a system of ODE's of the form
C###      <center>dy(i)/dt=RHSROUTINE(t,y(1),y(2),...,y(NUM_ODES))
C###      </center>
C###    where there are NUM_VARS state variables, of which the first
C###    NUM_ODES variables are ODE's.
C###    </html>
C***  Created by Martin Buist, June 1998

      IMPLICIT NONE

      INCLUDE 'adam00.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cell_reserved.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'integrator.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'lsoda00.cmn'
      INCLUDE 'nqloc00.inc'
      INCLUDE 'oxs005.cmn'
      INCLUDE 'tol00.cmn'
      INCLUDE 'integrator_reserved.inc'

      INTEGER INTEGRATOR_IWORK(INTEGRATOR_LIWORK),AII(*),AIO(*),
     '  CONTROL(*),maqdt,
     '  MODEL(*),NUM_ODES,NUM_VARS,nr,nx,SIZES(11),VARIANT
      REAL*8 INTEGRATOR_WORK(INTEGRATOR_LWORK),AQ(NMAQM),ARI(*),ARO(*),
     '  DERIVED(*),
     '  PARAMETERS(*),PROTOCOL(*),STIMULUS,T,YQS(*)
C MLB use (*) to fix linux problems
C      INTEGER SIZES(10)
C      INTEGER ADAMS_IWORK(ADAMS_LIWORK),AII(SIZES(NUM_AII)),
C     '  AIO(SIZES(NUM_AIO)),CONTROL(SIZES(NUM_CONTROL)),maqdt,
C     '  MODEL(SIZES(NUM_MODEL)),nr,nx,VARIANT
C      REAL*8 ADAMS_WORK(ADAMS_LWORK),AQ(NMAQM),ARI(SIZES(NUM_ARI)),
C     '  ARO(SIZES(NUM_ARO)),DERIVED(SIZES(NUM_DERIVED)),
C     '  PARAMETERS(SIZES(NUM_PARAM)),PROTOCOL(SIZES(NUM_PROTOCOL)),
C     '  STIMULUS,T,YQS(SIZES(NUM_EQN))
      CHARACTER ERROR*(*)
      EXTERNAL RHSROUTINE
!     Local Variables
      INTEGER MAX_NUM_VARS
      PARAMETER(MAX_NUM_VARS=100)
      INTEGER ERR,i,loctime,LOCTIMET
      REAL*8 REAL_WORK(5*MAX_NUM_VARS),locdt,TIME(2)
      REAL*8 DY1(MAX_NUM_VARS)
      REAL*8 STARTT,ENDT
      LOGICAL CHDT,STIFF_EQNS
      INTEGER INT_CONTROL(SIZE_INTEGER_CONTROL),jdum
      REAL*8 REAL_CONTROL(SIZE_REAL_CONTROL)
      integer liw,lrw
C      REAL*8 DY2(MAX_NUM_VARS),DY3(MAX_NUM_VARS),DY4(MAX_NUM_VARS),
C     '  YSTAR(MAX_NUM_VARS)
C      LOGICAL GATE

      CALL ENTERS('INTEGRATOR',*9999)


C      DT2=DT
C      T2=T
C      IF(ITYP3(nr,nx).EQ.8) THEN !N98
C        DT2=DT2/1000.0d0 !seconds
C        T2=T2/1000.0d0   !seconds
C      ENDIF

C *** If not solving any real cellular models shouldn't
C     do this - the stimulus can be used straight from ipcell
      IF(ITYP19(nr,nx).EQ.1.AND.ITYP3(nr,nx).LT.5)
     '  PROTOCOL(Is1current)=STIMULUS

      ERR=0
      IF(NUM_VARS.LE.MAX_NUM_VARS) THEN
        IF(KTYP37.EQ.1) THEN !Euler
          CHDT=.FALSE.

 500      IF(KTYP23.EQ.1) THEN
            LOCTIMET=1
            AQ(maqdt)=DT
          ELSE
            LOCTIMET=INT(DT/AQ(maqdt))
          ENDIF
          DO loctime=1,LOCTIMET
            locdt=AQ(maqdt)

            IF(.NOT.CHDT) THEN
              TIME(TCell) = T
              TIME(DTCell) = locdt
              CALL RHSROUTINE(TIME,YQS,DY1,CONTROL,MODEL,SIZES,VARIANT,
     '          DERIVED,PARAMETERS,PROTOCOL,AII,AIO,ARI,ARO,ERR)
              IF(ERR.NE.0) THEN
                ERROR='>>Error reported from RHSROUTINE'
                GOTO 9999
              ENDIF
            ENDIF

            CHDT=.FALSE.
            IF((loctime.EQ.1).AND.(KTYP23.EQ.2)) THEN
              IF(DABS(DY1(1)*locdt).GT.(1.0d0+ZERO_TOL)) THEN
                AQ(maqdt)=AQ(maqdt)/10.0d0
                CHDT=.TRUE.
                GOTO 500
              ELSEIF(DABS(DY1(1)*locdt).LT.(0.1d0-ZERO_TOL)) THEN
                IF(AQ(maqdt)*10.0d0.LE.DT) THEN
                  AQ(maqdt)=AQ(maqdt)*10.0d0
                  CHDT=.TRUE.
                  GOTO 500
                ENDIF
              ENDIF
            ENDIF

            DO i=1,NUM_ODES
              YQS(i)=YQS(i)+DY1(i)*locDT
            ENDDO
          ENDDO !time

C          IF(ITYP3(nr,nx).EQ.6) THEN !LR
C            IF(YQS(8).LT.ZERO_TOL) YQS(8)=0.0d0
C          ENDIF

        ELSE IF(KTYP37.EQ.2) THEN !Improved Euler

          ! Set-up the control arrays
          INT_CONTROL(ERROR_CODE) = 0
          INT_CONTROL(REAL_WORK_SIZE) = 5*MAX_NUM_VARS
          INT_CONTROL(GET_WORK_SIZES) = 0
          INT_CONTROL(SET_DEFAULT_VALUES) = 0
          CALL IMPROVED_EULER(AII,AIO,CONTROL,MODEL,NUM_ODES,NUM_VARS,
     '      SIZES,VARIANT,ARI,ARO,DERIVED,PARAMETERS,PROTOCOL,
     '      T,DT,TIME,YQS,RHSROUTINE,INT_CONTROL,REAL_CONTROL,
     '      REAL_WORK,ERROR)
          IF(INT_CONTROL(ERROR_CODE).NE.0) THEN
            GOTO 9999
          ENDIF

        ELSE IF(KTYP37.EQ.3) THEN !4th order Runge-Kutta

          ! Set-up the control arrays
          INT_CONTROL(ERROR_CODE) = 0
          INT_CONTROL(REAL_WORK_SIZE) = 5*MAX_NUM_VARS
          INT_CONTROL(GET_WORK_SIZES) = 0
          INT_CONTROL(SET_DEFAULT_VALUES) = 0
          CALL RUNGE_KUTTA(AII,AIO,CONTROL,MODEL,NUM_ODES,NUM_VARS,
     '      SIZES,VARIANT,ARI,ARO,DERIVED,PARAMETERS,PROTOCOL,
     '      T,DT,TIME,YQS,RHSROUTINE,INT_CONTROL,REAL_CONTROL,
     '      REAL_WORK,ERROR)
          IF(INT_CONTROL(ERROR_CODE).NE.0) THEN
            GOTO 9999
          ENDIF

        ELSE IF(KTYP37.EQ.4) THEN !not defined

          ERROR='>>Not implemented'
          GOTO 9999

        ELSE IF(KTYP37.EQ.5) THEN !Adams-Moulton

          IF(INTEGRATOR_IWORK(1).EQ.0) THEN
            ERR=1 !Startup
          ELSE
            IF(ITYP19(nr,nx).EQ.1.AND.ITYP3(nr,nx).LT.5) THEN
              IF(DABS(STIMULUS).GT.ZERO_TOL) THEN
                ERR=1 !Startup
              ELSE
                ERR=2 !Continue
              ENDIF
            ELSE
              ERR=2 !Continue
            ENDIF
          ENDIF
          ERR=1
          IF(ERR.EQ.1) THEN
            STARTT=T
            ENDT=T+DT
            INTEGRATOR_WORK(4) = 0.0d0 !Initialise Told for Adams integrator
            INTEGRATOR_WORK(5) = 0.0d0
            INTEGRATOR_WORK(6) = 0.0d0
          ELSE
            STARTT=INTEGRATOR_WORK(4)
            ENDT=T+DT
          ENDIF

          STIFF_EQNS=.FALSE.

          CALL ADAMS(AII,AIO,CONTROL,ADAMS_ERROR_TYPE,ERR,
     '      INTEGRATOR_IWORK,INTEGRATOR_LIWORK,INTEGRATOR_LWORK,
     '      ADAMS_MAX_ITERS,ADAMS_MAX_ORDER,MODEL,NUM_ODES,NUM_VARS,
     '      SIZES,VARIANT,ADAMS_ABS_ERR,ARI,ARO,DERIVED,DY1,
     '      ADAMS_MAX_STEP,PARAMETERS,PROTOCOL,ADAMS_REL_ERR,STARTT,
     '      ENDT,INTEGRATOR_WORK,YQS,.TRUE.,STIFF_EQNS,
     '      ADAMS_USE_ROUND_CTRL,
     '      RHSROUTINE,ERROR)

          IF(ERR.EQ.2) THEN
            ERR=0
          ELSE
            GOTO 9999
          ENDIF
          IF(STIFF_EQNS) THEN
            WRITE(OP_STRING,'('' >>WARNING: Equations appear to be '
     '        //'stiff'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF

        ELSEIF(KTYP37.EQ.6) THEN
           STARTT=T
           ENDT=T+DT

           liw = INTEGRATOR_LIWORK - LSODA_NUM_INT_COMMON
           lrw = INTEGRATOR_LWORK - LSODA_NUM_REAL_COMMON

           CALL LSODA_WRAPPER(RHSROUTINE,NUM_ODES,YQS,STARTT,ENDT,
     1       LSODA_ERROR_TYPE,LSODA_REL_ERR,LSODA_ABS_ERR,LSODA_ITASK,
     2       INTEGRATOR_IWORK(LSODA_ISTATE),LSODA_IOPT,
     3       INTEGRATOR_WORK(LSODA_NUM_REAL_COMMON+1),lrw,
     4       INTEGRATOR_IWORK(LSODA_NUM_INT_COMMON+1),liw,
     5       jdum,LSODA_JACOBIAN_TYPE,AII,AIO,CONTROL,MODEL,
     6       SIZES,VARIANT,ARI,ARO,DERIVED,PARAMETERS,PROTOCOL,ERR,
     7       INTEGRATOR_WORK,INTEGRATOR_IWORK,ERROR,*9999)

        ELSE
          ERROR='>>Unknown integration scheme'
          GOTO 9999

        ENDIF

      ELSE

        ERROR='>>Increase MAX_NUM_VARS in INTEGRATOR'
        GOTO 9999

      ENDIF

      CALL EXITS('INTEGRATOR')
      RETURN
 9999 CALL ERRORS('INTEGRATOR',ERROR)
      CALL EXITS('INTEGRATOR')
      RETURN 1
      END

