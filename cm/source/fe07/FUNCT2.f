      SUBROUTINE FUNCT2(MODE,NTRES,n,LDFJ,XC,FC,FJAC,NSTATE,IUSER,USER)

C#### Subroutine: FUNCT2
C###  Description:
C###    FUNCT2 evaluates a set of residuals, FC, for optimisation
C###    routine MINLSSQP, during the optimisation of material parameters.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'aero00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'coup00.cmn'
      INCLUDE 'chmesh0.cmn'
      INCLUDE 'data00.cmn'
      INCLUDE 'error0.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'gks000.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'ofst00.cmn'
      INCLUDE 'opti00.cmn'
      INCLUDE 'trsf00.cmn'
!     Parameter List
      INTEGER IUSER(*),LDFJ,MODE,n,NSTATE,NTRES
      REAL*8 FC(*),FJAC(LDFJ,*),USER(*),XC(*)

!     Local Variables
      INTEGER MXCOQU
      PARAMETER(MXCOQU=25)
      INTEGER IBEG,IEND,noopti,nores,nr,nx
      REAL*8 XCOLD(3)
C!!! CS 5/4/2000 changed length of ERROR to avoid interpreter bug
C!!! KAT to fix interpreter sometime
      CHARACTER BASELINE*(5),CHAR*(50),ERROR*(ERRSTRLEN),FREQUENCY*(8),
     '  OPTI_REGION*(2),SQUID_CONFIG*(1),STRING*(MXCH),TIME*(8)
      LOGICAL EV_FUN,EV_JAC,END

      CALL ENTERS('FUNCT2',*9999)

      nr=1 !temporary
      nx=1 !temporary

      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP   CRITICAL(FUNCT2_1)
        WRITE(OP_STRING,'('' MODE='',I3,'' NTRES='',I3,'' N='',I3,'
     '    //''' NSTATE='',I3,'' LDFJ='',I6)') MODE,NTRES,n,NSTATE,LDFJ
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP   END CRITICAL(FUNCT2_1)
      ENDIF

      IF(KTYP26.EQ.1.AND.KTYP27.EQ.2) THEN
C **    Residuals are reaction diff.s
        DO noco=1,MXCOQU
          NTCOQU(noco)=0
        ENDDO
        DO noopti=1,NTOPTI
          USER(OS_PAOPTI+noopti-1)=XC(noopti)
        ENDDO

        IF(ITYP2(nr,nx).EQ.2) THEN !Finite elasticity problem
          IF(KTYP51(nr).EQ.3.AND.KTYP52(nr).EQ.3) THEN
C ***       For incompressibilty + fluid movement in 3D case:
C ***       update the aux elem params, given the current set of
C ***       mat params, the given undeformed, deformed geometries
C ***       and the pressure bc's
            CO(1)='FEM'
            CO(2)='UPDATE'
            CO(3)='PRESSURE'
            CO(4)='OPTIMISATION_PARAMETERS'
            NTCO=4
            noco=1
            IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP         CRITICAL(FUNCT2_2)
              WRITE(OP_STRING,'('' fem update pressure optim.n..'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP         END CRITICAL(FUNCT2_2)
            ENDIF
            CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '        ERROR,*9999)
          ENDIF

          IF(KTYP28.EQ.0) THEN
            CO(1)='FEM'
            CO(2)='SOLVE'
            CO(3)='STEP'
            CO(4)='1'
            CO(5)='UPDATE'
            CO(6)='1'
            NTCO=6
            noco=1
            IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP         CRITICAL(FUNCT2_3)
              WRITE(OP_STRING,'('' fem solve step 1 update 1..'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP         END CRITICAL(FUNCT2_3)
            ENDIF
            CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '        ERROR,*9999)
          ENDIF

        ELSE IF(ITYP5(nr,nx).EQ.2.AND.ITYP19(nr,nx).EQ.1.AND.
     '      ITYP2(nr,nx).EQ.9) THEN !activation model

          CO(1)='FEM'
          CO(2)='UPDATE'
          CO(3)='OPTIMISATION'
          NTCO=3
          noco=1
          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP       CRITICAL(FUNCT2_4)
            WRITE(OP_STRING,'('' fem update optimisation..'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP       END CRITICAL(FUNCT2_4)
          ENDIF
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)

          CO(1)='FEM'
          CO(2)='UPDATE'
          CO(3)='GRID'
          CO(4)='MATERIAL'
          NTCO=4
          noco=1
          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP       CRITICAL(FUNCT2_5)
            WRITE(OP_STRING,'('' fem update grid material..'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP       END CRITICAL(FUNCT2_5)
          ENDIF
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)

          CO(1)='FEM'
          CO(2)='DEFINE'
          CO(3)='INITIAL'
          NTCO=3
          noco=1
          COQU(3,1)='R'
          COQU(3,2)='OPTIMISE'     !Name of input file
          NTCOQU(1)=0
          NTCOQU(2)=0
          NTCOQU(3)=2
          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP       CRITICAL(FUNCT2_6)
            WRITE(OP_STRING,'('' fem define initial;r;optimise..'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP       END CRITICAL(FUNCT2_6)
          ENDIF
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)

          CO(1)='FEM'
          CO(2)='DEFINE'
          CO(3)='SOLVE'
          NTCO=3
          noco=1
          COQU(3,1)='R'
          COQU(3,2)='OPTIMISE'     !Name of input file
          NTCOQU(1)=0
          NTCOQU(2)=0
          NTCOQU(3)=2
          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP       CRITICAL(FUNCT2_7)
            WRITE(OP_STRING,'('' fem define solve;r;optimise..'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP       END CRITICAL(FUNCT2_7)
          ENDIF
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)
          NTCOQU(1)=0
          NTCOQU(2)=0
          NTCOQU(3)=0

          CO(1)='FEM'
          CO(2)='SOLVE'
          NTCO=2
          noco=1
          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP       CRITICAL(FUNCT2_8)
            WRITE(OP_STRING,'('' fem solve..'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP       END CRITICAL(FUNCT2_8)
          ENDIF
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)

        ENDIF !ityp2(1)

        IF(MODE.EQ.0.OR.MODE.EQ.2) THEN
C ***     evaluate residuals
          CO(1)='FEM'
          CO(2)='EVALUATE'
          CO(3)='RESIDUALS'
          CO(4)='WRT'
          CO(5)='MAT_PARAMS'
          CO(6)='NOVIEW'
          NTCO=6
          noco=1
          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP       CRITICAL(FUNCT2_9)
            WRITE(OP_STRING,'('' fem evaluate residuals '
     '        //'wrt mat_params noview..'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP       END CRITICAL(FUNCT2_9)
          ENDIF
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '     ERROR,*9999)
          DO nores=1,NTRES
           FC(nores)=USER(OS_RESID+nores-1)
         ENDDO
        ENDIF

        IF(IPDLEV.GT.0.AND.(MODE.EQ.1.OR.MODE.EQ.2)) THEN
C ***     evaluate first derivatives of residuals wrt opt params
          CO(1)='FEM'
          CO(2)='EVALUATE'
          CO(3)='RESIDUALS'
          CO(4)='ANALYT_DERIVS'
          CO(5)='WRT'
          CO(6)='MAT_PARAMS'
          CO(7)='NOVIEW'
          NTCO=7
          noco=1
          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP       CRITICAL(FUNCT2_10)
            WRITE(OP_STRING,'('' fem evaluate residuals analyt_derivs'
     '        //' wrt mat_params noview..'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP       END CRITICAL(FUNCT2_10)
          ENDIF
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)
          DO nores=1,NTRES
            DO noopti=1,NTOPTI
              FJAC(nores,noopti)=USER(OS_RESJAC+nores-1+
     '          (noopti-1)*NREM)
            ENDDO
          ENDDO
        ENDIF

      ELSE IF(KTYP26.EQ.2.AND.KTYP27.EQ.1) THEN
C **    Objective function is minimum area of trapezoids
        DO noco=1,6
          NTCOQU(noco)=0
        ENDDO
        DO noopti=1,NTOPTI
          USER(OS_PAOPTI+noopti-1)=XC(noopti)
        ENDDO

        CO(1)='FEM'
        CO(2)='DEFINE'
        CO(3)='POLYLINE'
        COQU(3,1)='c'
        CO(4)='WITH'
        CO(5)='OPT'
        CO(6)='ENDPOINT'
        CO(7)='30'
        NTCO=7
        NTCOQU(3)=1
        noco=1
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)

        CO(1)='FEM'
        CO(2)='DEFINE'
        CO(3)='FIT'
        COQU(3,1)='r'
        CO(4)='GEOM'
        NTCO=4
        NTCOQU(3)=1
        noco=1
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)

        CO(1)='FEM'
        CO(2)='FIT'
        CO(3)='GEOM'
        NTCO=3
        noco=1
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)

        CO(1)='FEM'
        CO(2)='UPDATE'
        CO(3)='NODES'
        NTCO=3
        noco=1
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)

        CO(1)='FEM'
        CO(2)='LIST'
        CO(3)='DATA'
        CO(4)='ERROR'
        NTCO=4
        noco=1
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)

        FUNC=SQED
        FC(1)=FUNC

        IF(UPVUOP) THEN
          CO(1)='FEM'
          CO(2)='DEFINE'
          CO(3)='LINE'
          COQU(3,1)='s'
          NTCO=3
          NTCOQU(3)=1
          noco=1
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)

          CO(1)='FEM'
          CO(2)='DEFINE'
          CO(3)='DATA'
          COQU(3,1)='s'
          NTCO=3
          NTCOQU(3)=1
          noco=1
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)

          CO(1)='FEM'
          CO(2)='DEFINE'
          CO(3)='DATA'
          COQU(3,1)='s'
          CO(4)='PROJECTIONS'
          NTCO=4
          NTCOQU(3)=1
          noco=1
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)
        ENDIF

      ELSE IF(KTYP26.LE.2.AND.KTYP27.EQ.4) THEN
C **    Objective function is hydrostatic pressure condition
        DO noco=1,6
          NTCOQU(noco)=0
        ENDDO
        DO noopti=1,NTOPTI
          USER(OS_PAOPTI+noopti-1)=XC(noopti)
        ENDDO

        CO(1)='FEM'
        CO(2)='SOLVE'
        NTCO=2
        noco=1
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)

        CO(1)='FEM'
        CO(2)='EVALUATE'
        CO(3)='OBJECTIVE'
        NTCO=3
        noco=1
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)
        FC(1)=FUNC

        IF(UPVUOP) THEN
          CO(1)='FEM'
          CO(2)='UPDATE'
          CO(3)='MESH'
          NTCO=3
          noco=1
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)
        ENDIF

      ELSE IF(KTYP26.EQ.2.AND.KTYP27.EQ.5) THEN !data fitting
C **    Objective function is sum of squares of all the orthogonal
C **    projection distances of the data points

        CALL ASSERT(MODE.GE.0.AND.MODE.LE.2,'Invalid mode',ERROR,*9999)
        IF(KTYP29B.EQ.1) THEN !residuals are components
C         Only need to calculated Jacobian once
          EV_JAC=NSTATE.EQ.1
        ELSE
          EV_JAC=MODE.EQ.1.OR.MODE.EQ.2
        ENDIF
        EV_FUN=MODE.EQ.0.OR.MODE.EQ.2

        DO noco=1,6
          NTCOQU(noco)=0
        ENDDO
        DO noopti=1,NTOPTI
          USER(OS_PAOPTI+noopti-1)=XC(noopti)
        ENDDO

        IF(EV_JAC) THEN
          DO noopti=1,NTOPTI
            DO nores=1,NTRES
              USER(OS_RESJAC+nores-1+(noopti-1)*NREM)=FJAC(nores,noopti)
            ENDDO
          ENDDO
        ENDIF

        IF(EV_FUN.OR.EV_JAC) THEN
          CO(1)='FEM'
          CO(2)='EVALUATE'
          CO(3)='RESIDUALS'
          CO(4)='WRT'
          CO(5)='DATA_FITTING'
          CO(6)='NOVIEW'
          IF(EV_JAC.AND.EV_FUN) THEN
            CO(7)='BOTH'
          ELSE IF(EV_JAC) THEN
            CO(7)='JACOBIAN'
          ELSE IF(EV_FUN) THEN
            CO(7)='FUNCTION'
          ENDIF
          NTCO=7
          noco=1
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)
        ENDIF !EV_FUN.OR.EV_JAC

        IF(MODE.EQ.0.OR.MODE.EQ.2) THEN
          DO nores=1,NTRES
            FC(nores)=USER(OS_RESID+nores-1)
          ENDDO
        ENDIF

        IF(IPDLEV.GE.0.AND.(MODE.EQ.1.OR.MODE.EQ.2)) THEN
          DO noopti=1,NTOPTI
            DO nores=1,NTRES
              FJAC(nores,noopti)=USER(OS_RESJAC+nores-1+(noopti-1)*NREM)
            ENDDO
          ENDDO
        ENDIF

C CS 28/4/2000 KTYP26.EQ.3.AND.KTYP27.EQ.1 is not a valid option
C This code appears to be obsolete
C      ELSE IF(KTYP26.EQ.3.AND.KTYP27.EQ.1) THEN !Dave's stripe calcs.
C        xksia=xc(1)
C        xksib=xc(2)
C        if (xksia.lt.1.0d-6) xksia=1.0d-6
C        if (xksia.gt.0.999999d0) xksia=0.999999d0     !Ask Peter how
C        b(1)=(1.0d0-xksia)**3                        !to deal with this
C        b(2)=3.0d0*xksia*(1.0d0-xksia)**2               !underflow problem
C        b(3)=3.0d0*xksia*xksia*(1.0d0-xksia)
C        b(4)=xksia*xksia*xksia
C        xxa=0.0d0
C        yya=0.0d0
C        do i=1,4
C          xxa=xxa+b(i)*xa(i)
C          yya=yya+b(i)*ya(i)
C        end do   !end i

C        if (xksib.lt.1.0d-6) xksib=1.0d-6
C        if (xksib.gt.0.999999d0) xksib=0.999999d0
C        b(1)=(1.0d0-xksib)**3
C        b(2)=3.0d0*xksib*(1.0d0-xksib)**2
C        b(3)=3.0d0*xksib*xksib*(1.0d0-xksib)
C        b(4)=xksib*xksib*xksib
C        xxb=0.0d0
C        yyb=0.0d0
C        do i=1,4
C          xxb=xxb+b(i)*xb(i)
C          yyb=yyb+b(i)*yb(i)
C        end do    !end i

C        fc(1)=((yya-yyb)**2+(xxa-xxb)**2)**0.5d0

C CS 10/9/2001 This appear to be obsolete
C KTYP26.EQ.4 now refers to holmes constitutive law constants optimisation
C CS start
C      ELSE IF(KTYP26.EQ.4.AND.KTYP27.EQ.1) THEN !Dave's stripe calcs.
c ...   xksia passed in with common block
C        xksib=xc(1)
C        acdiff=0.0d0
C        acn=0.0d0
C        bcn=0.0d0
C        w1=0.0d0
C        xka=0.0d0
C        do kk=1,nsimp+1      !Simpson's Rule
C          xka=xksia+(xksib-xksia)*(kk-1.0d0)/dble(nsimp)
C          acn=xa(1)*(-3.0d0*(1.0d0-xka)**2)
C          acn=acn+xa(2)*(3.0d0*(1.0d0-xka)*(1.0d0-3.0d0*xka))
C          acn=acn+xa(3)*(3.0d0*xka*(2.0d0-3.0d0*xka))
C          acn=acn+xa(4)*3.0d0*xka*xka
C
C          bcn=ya(1)*(-3.0d0*(1.0d0-xka)**2)
C          bcn=bcn+ya(2)*(3.0d0*(1.0d0-xka)*(1.0d0-3.0d0*xka))
C          bcn=bcn+ya(3)*(3.0d0*xka*(2.0d0-3.0d0*xka))
C          bcn=bcn+ya(4)*3.0d0*xka*xka
C
C          w1=2.0d0                 !presume kk is odd
C          kk1=kk/2
C          kk1=kk1*2
C          if (kk1.eq.kk) w1=4.0d0   !kk is even
C          if (kk.eq.1.or.kk.eq.nsimp+1) w1=1.0d0    !kk is endpt.
C          acdiff=acdiff+((acn*acn+bcn*bcn)**0.5d0)*w1*(xksib-xksia)
C     /           /dble(nsimp)/3.0d0
C        end do !end kk
C
C        fc(1)=1.0d0+(slen-acdiff)**2     !exact length - estimated length
C
C
C      ELSE IF(KTYP26.EQ.4.AND.KTYP27.EQ.2) THEN !xi calculation in
C                                                !cubic-linear element
Cc ...   Find x_endo         !x coordinate of endocardial point
C                            !corresponding to trial value of xi1
C        XI1=XC(1)
C        S01=1.0d0-3.0d0*XI1*XI1+2.0d0*XI1*XI1*XI1  !evaluate hermite basis fn
C        S02=3.0d0*XI1*XI1-2.0d0*XI1*XI1*XI1
C        S11=XI1-2.0d0*XI1*XI1+XI1*XI1*XI1
C        S12=-(XI1*XI1-XI1*XI1*XI1)
C        F01=XNODVAL(1,1)                  !nodal values
C        F02=XNODVAL(3,1)
C        F11=XNODVAL(2,1)
C        F12=XNODVAL(4,1)
C        X_ENDO=S01*F01+S02*F02+S11*F11+S12*F12
C
CC ...   Find y_endo         !y coordinate of endocardial point
C                            !corresponding to trial value of xi1
C        F01=XNODVAL(1,2)
C        F02=XNODVAL(3,2)
C        F11=XNODVAL(2,2)
C        F12=XNODVAL(4,2)
C        Y_ENDO=S01*F01+S02*F02+S11*F11+S12*F12
C
CC ...   Find x_epi          !x coordinate of epicardial point
C                            !corresponding to trial value of xi1
C        F01=XNODVAL(5,1)
C        F02=XNODVAL(7,1)
C        F11=XNODVAL(6,1)
C        F12=XNODVAL(8,1)
C        X_EPI=S01*F01+S02*F02+S11*F11+S12*F12
C
CC ...   Find y_epi          !y coordinate of epicardial point
C                            !corresponding to trial value of xi1
C        F01=XNODVAL(5,2)
C        F02=XNODVAL(7,2)
C        F11=XNODVAL(6,2)
C        F12=XNODVAL(8,2)
C        Y_EPI=S01*F01+S02*F02+S11*F11+S12*F12
C
CC ...   Calculate (XNPB,YNPB) the nearest point to (XPXIF,YPXIF) on the
CC       line drawn between (X_ENDO,Y_ENDO) and (X_EPI,Y_EPI)
C        XDIFFERENCE=DABS(X_EPI-X_ENDO)
C        IF(XDIFFERENCE.LT.1.0D-8) THEN      !Avoid singular case
C          XNPB=X_EPI
C          YNPB=YPXIF
C        ELSE                                !Nonsingular case
C          S1=X_ENDO-X_EPI
C          S2=Y_ENDO-Y_EPI
C          XNPB=(XPXIF*S1*S1+X_ENDO*S2*S2-S1*S2*(Y_ENDO-YPXIF))
C     '          /(S1*S1+S2*S2)
C          YNPB=Y_ENDO+S2*(XNPB-X_ENDO)/S1
C        ENDIF
C
C        FC(1)=(XNPB-XPXIF)**2+(YNPB-YPXIF)**2   !value of obj. function
C
C
C      ELSE IF(KTYP26.EQ.4.AND.KTYP27.EQ.3) THEN !xi calculation in
C                                                !Lagrange quad element
C        XI1=XC(1)
C        XI2=XC(2)
C
C        VL20=2.0d0*(XI1-0.5d0)*(XI1-1.0d0)     !shape factors
C        VL21=-4.0d0*(XI1-1.0d0)*XI1
C        VL22=2.0d0*XI1*(XI1-0.5d0)
C        CL20=2.0d0*(XI2-0.5d0)*(XI2-1.0d0)
C        CL21=-4.0d0*(XI2-1.0d0)*XI2
C        CL22=2.0d0*XI2*(XI2-0.5d0)
C
Cc ...   Pass array XNODVAL containing nodal values via common
C
C        XNPB=VL20*CL20*XNODVAL(1,1)+VL21*CL20*XNODVAL(2,1)+
C     '       VL22*CL20*XNODVAL(3,1)+VL20*CL21*XNODVAL(4,1)+
C     '       VL21*CL21*XNODVAL(5,1)+VL22*CL21*XNODVAL(6,1)+
C     '       VL20*CL22*XNODVAL(7,1)+VL21*CL22*XNODVAL(8,1)+
C     '       VL22*CL22*XNODVAL(9,1)
C
C        YNPB=VL20*CL20*XNODVAL(1,2)+VL21*CL20*XNODVAL(2,2)+
C     '       VL22*CL20*XNODVAL(3,2)+VL20*CL21*XNODVAL(4,2)+
C     '       VL21*CL21*XNODVAL(5,2)+VL22*CL21*XNODVAL(6,2)+
C     '       VL20*CL22*XNODVAL(7,2)+VL21*CL22*XNODVAL(8,2)+
C     '       VL22*CL22*XNODVAL(9,2)
C
C        FC(1)=(XNPB-XPXIF)**2+(YNPB-YPXIF)**2 !value of obj. function

      ELSE IF(KTYP26.EQ.2.AND.KTYP27.EQ.6) THEN
C **    Optimise geom params by minimizing fluid interface residual
        DO noco=1,6
          NTCOQU(noco)=0
        ENDDO
        DO noopti=1,NTOPTI
          USER(OS_PAOPTI+noopti-1)=XC(noopti)
        ENDDO

        CO(1)='FEM'
        CO(2)='UPDATE'
        CO(3)='NODE'
        CO(4)='INTERFACE'
        CO(5)='FLUX'
        NTCO=5
        noco=1
        IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP     CRITICAL(FUNCT2_11)
          WRITE(OP_STRING,'('' fem update node interface flux..'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP     END CRITICAL(FUNCT2_11)
        ENDIF
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)

        CO(1)='FEM'
        CO(2)='SOLVE'
        CO(3)='FOR'
        CO(4)='1'
        NTCO=4
        noco=1
        IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP     CRITICAL(FUNCT2_12)
          WRITE(OP_STRING,'('' fem solve for 1..'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP     END CRITICAL(FUNCT2_12)
        ENDIF
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)

        CO(1)='FEM'
        CO(2)='SOLVE'
        CO(3)='FOR'
        CO(4)='2'
        NTCO=4
        noco=1
        IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP     CRITICAL(FUNCT2_13)
          WRITE(OP_STRING,'('' fem solve for 2..'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP     END CRITICAL(FUNCT2_13)
        ENDIF
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)

        CO(1)='FEM'
        CO(2)='EVALUATE'
        CO(3)='OBJECTIVE'
        NTCO=3
        noco=1
        IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP     CRITICAL(FUNCT2_14)
          WRITE(OP_STRING,'('' fem evaluate objective..'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP     END CRITICAL(FUNCT2_14)
        ENDIF
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)
        FC(1)=FUNC

        IF(UPVUOP) THEN
          CO(1)='FEM'
          CO(2)='UPDATE'
          CO(3)='NODE'
          CO(4)='INTERFACE'
          CO(5)='INCREMENT'
          NTCO=5
          noco=1
          WRITE(OP_STRING,'('' fem update node interface pos..'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)

          CO(1)='FEM'
          CO(2)='DEFINE'
          CO(3)='LINE'
          COQU(3,1)='s'
          NTCO=3
          noco=1
          NTCOQU(3)=1
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)
        ENDIF

      ELSE IF(KTYP26.EQ.2.AND.KTYP27.EQ.7) THEN
C **    Objective fn resids are aerofoil wake press diffs + stress
        DO noco=1,7
          NTCOQU(noco)=0
        ENDDO
        DO noopti=1,NTOPTI
          USER(OS_PAOPTI+noopti-1)=XC(noopti) !-->PAOPTI in FEM
c         IF(DOP) THEN
            WRITE(OP_STRING,'('' XC('',I2,'')= '',D20.12)')
     '        noopti,XC(noopti)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c         ENDIF
        ENDDO

C       put opti params into wake fe nodes for wake and sail positions
        CO(1)='FEM'
        CO(2)='UPDATE'
        CO(3)='OPTIMISATION'
        NTCO=3
        noco=1
        IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP     CRITICAL(FUNCT2_15)
          WRITE(OP_STRING,'('' fem update optimisation..  NTRES='','
     '      //'I3)') NTRES
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP     END CRITICAL(FUNCT2_15)
        ENDIF
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)

C       update phi on entry face and flux on exit face
        CO(1)='FEM'
        CO(2)='UPDATE'
        CO(3)='AEROFOIL'
        NTCO=3
        noco=1
        IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP     CRITICAL(FUNCT2_16)
          WRITE(OP_STRING,'('' fem update aerofoil..'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP     END CRITICAL(FUNCT2_16)
        ENDIF
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)

C       solve for phi
        CO(1)='FEM'
        CO(2)='SOLVE'
        NTCO=2
        noco=1
        IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP     CRITICAL(FUNCT2_17)
          WRITE(OP_STRING,'('' fem solve.. (for pot. field)'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP     END CRITICAL(FUNCT2_17)
        ENDIF
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)

C       calculate pressure on sail
        CO(1)='FEM'
        CO(2)='EVALUATE'
        CO(3)='AEROFOIL'
        NTCO=3
        noco=1
        IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP     CRITICAL(FUNCT2_18)
          WRITE(OP_STRING,'('' fem evaluate aerofoil..'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP     END CRITICAL(FUNCT2_18)
        ENDIF
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)

        IF(N_OPTI(2).GT.0) THEN !sail stress parameters included
C         Transfer sail pressure to array PF
          CO(1)='FEM'
          CO(2)='UPDATE'
          CO(3)='COUPLING'
          NTCO=3
          noco=1
          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP       CRITICAL(FUNCT2_19)
            WRITE(OP_STRING,'('' fem update coupling..'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP       END CRITICAL(FUNCT2_19)
          ENDIF
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)

C         Evaluate stress residuals
          CO(1)='FEM'
          CO(2)='EVALUATE'
          CO(3)='RESIDUALS'
          CO(4)='WRT'
          CO(5)='GEOM_PARAMS'
          CO(6)='FOR'
          CO(7)='2'
          CO(8)='NOVIEW'
          NTCO=8
          noco=1
          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP       CRITICAL(FUNCT2_20)
            WRITE(OP_STRING,'('' fem eval resid wrt geom for 2 '
     '        //'noview..'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP       END CRITICAL(FUNCT2_20)
          ENDIF
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)
        ENDIF

        DO nores=1,NTRES
          IF(nores.LE.N_OPTI(1)) THEN !wake parameters
            FC(nores)=dPHI_resid(nores)
c           IF(DOP) THEN
              WRITE(OP_STRING,'('' FC('',I2,'')= '',D20.12,'
     '          //''' (wake dPHI residuals)'')') nores,FC(nores)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c           ENDIF
          ELSE IF(nores.GT.N_OPTI(1)) THEN !sail stress parameters
            FC(nores)=USER(OS_RESID+nores-1)
c           IF(DOP) THEN
              WRITE(OP_STRING,'('' FC('',I2,'')= '',D20.12,'
     '          //''' (sail stress residual)'')')     nores,FC(nores)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c           ENDIF
          ENDIF
        ENDDO

        IF(MODE.EQ.0.OR.MODE.EQ.2) THEN      !evaluate residuals
        ELSE IF(MODE.EQ.1.OR.MODE.EQ.2) THEN !evaluate 1st derivs
        ENDIF

      ELSE IF(KTYP26.EQ.2.AND.KTYP27.EQ.8) THEN
C **    Objective fn residual is aerofoil lift
        DO noco=1,6
          NTCOQU(noco)=0
        ENDDO
        DO noopti=1,NTOPTI
          USER(OS_PAOPTI+noopti-1)=XC(noopti) !-->PAOPTI in FEM
          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP       CRITICAL(FUNCT2_21)
            WRITE(OP_STRING,'('' XC('',I2,'')= '',D20.12)')
     '        noopti,XC(noopti)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP       END CRITICAL(FUNCT2_21)
          ENDIF
        ENDDO

        CO(1)='FEM'
        CO(2)='UPDATE'
        CO(3)='OPTIMISATION'
        NTCO=3
        noco=1
        IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP     CRITICAL(FUNCT2_22)
          WRITE(OP_STRING,'('' fem update optimisation..  NTRES='','
     '      //'I3)') NTRES
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP     END CRITICAL(FUNCT2_22)
        ENDIF
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '    ERROR,*9999)

        CO(1)='FEM'
        CO(2)='UPDATE'
        CO(3)='AEROFOIL'
        NTCO=3
        noco=1
        IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP     CRITICAL(FUNCT2_23)
          WRITE(OP_STRING,'('' fem update aerofoil..'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP     END CRITICAL(FUNCT2_23)
        ENDIF
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)

        CO(1)='FEM'
        CO(2)='SOLVE'
        NTCO=2
        noco=1
        IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP     CRITICAL(FUNCT2_24)
          WRITE(OP_STRING,'('' fem solve.. (for potential field)'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP     END CRITICAL(FUNCT2_24)
        ENDIF
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)

        CO(1)='FEM'
        CO(2)='EVALUATE'
        CO(3)='AEROFOIL'
        NTCO=3
        noco=1
        IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP     CRITICAL(FUNCT2_25)
          WRITE(OP_STRING,'('' fem evaluate aerofoil..'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP     END CRITICAL(FUNCT2_25)
        ENDIF
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)

        DO nores=1,NTRES-1
          FC(nores)=PRESS_DIFF_WAKE(nores)
        ENDDO
        IF(UPPER_BOUND.GT.1.0D-6) THEN
          FC(NTRES)=UPPER_BOUND-DABS(TOT_LIFT)
        ELSE
          FC(NTRES)=DABS(TOT_LIFT)
        ENDIF
        IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP     CRITICAL(FUNCT2_26)
          DO nores=1,NTRES
            IF(nores.LT.NTRES.AND.NTRES.GT.1) THEN
              WRITE(OP_STRING,'('' FC('',I2,'')= '',D20.12,'
     '          //''' (wake pressure difference)'')')
     '          nores,FC(nores)
            ELSE
              WRITE(OP_STRING,'('' FC('',I2,'')= '',D20.12,'
     '          //''' (upper bound+/-abs(aerofoil lift))'')')
     '          nores,FC(nores)
            ENDIF
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO
CC$OMP     END CRITICAL(FUNCT2_26)
        ENDIF

        IF(MODE.EQ.0.OR.MODE.EQ.2) THEN !evaluate residuals
        ENDIF

      ELSE IF(KTYP26.EQ.2.AND.KTYP27.EQ.10) THEN !torso opti
        DO noco=1,7
          NTCOQU(noco)=0
        ENDDO
        DO noopti=1,NTOPTI
          USER(OS_PAOPTI+noopti-1)=XC(noopti)
        ENDDO

        CO(1)='FEM'
        CO(2)='UPDATE'
        CO(3)='GEOMETRY'
        CO(4)='FROM'
        CO(5)='FIELD'
        CO(6)='REGION'
        CO(7)='ALL'
        NTCO=7
        noco=1
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)

        IF(CUSTOMISATION_TYPE.EQ.1.OR.CUSTOMISATION_TYPE.EQ.3.OR.
     '    CUSTOMISATION_TYPE.EQ.4) THEN
          CO(1)='FEM'
          CO(2)='CHANGE'
          CO(3)='MESH'
          CO(4)='FROM'
          CO(5)='OPTIMISER'
          NTCO=5
          noco=1
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,
     '      *9999)

        ELSEIF(CUSTOMISATION_TYPE.EQ.2) THEN
          CO(1)='FEM'
          CO(2)='CHANGE'
          CO(3)='MESH'
          CO(4)='FROM'
          CO(5)='OPTIMISER'
          CO(6)='REGION'
          CO(7)='2'
          NTCO=7
          noco=1
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)

        ELSEIF(CUSTOMISATION_TYPE.EQ.3) THEN


        ENDIF


        CO(1)='FEM'
        CO(2)='EVALUATE'
        CO(3)='RESIDUALS'
        CO(4)='WRT'
        CO(5)='TORSO_CUSTOMISE'
        CO(6)='NOVIEW'
        NTCO=6
        noco=1
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)


        IF(MODE.EQ.0.OR.MODE.EQ.2) THEN
          DO nores=1,NTRES
            FC(nores)=USER(OS_RESID+nores-1)
          ENDDO
        ENDIF
C No jacobian values available
C        IF(MODE.EQ.1.OR.MODE.EQ.2) THEN
C         DO noopti=1,NTOPTI
C           DO nores=1,NTRES
C             FJAC(nores,noopti)=USER(OS_RESJAC+nores-1+
C     '          (noopti-1)*NREM)
C           ENDDO
C         ENDDO
C        ENDIF
      ELSE IF(KTYP26.EQ.2.AND.KTYP27.EQ.11) THEN !second moments
        DO noco=1,6
          NTCOQU(noco)=0
        ENDDO

        DO noopti=1,NTOPTI
          XCOLD(noopti)=USER(OS_PAOPTI+noopti-1)
        ENDDO

        WRITE(CHAR,'(D12.4,'','',D12.4,'','',D12.4)')-XCOLD(1),
     '    -XCOLD(2),-XCOLD(3)
        CO(1)='FEM'
        CO(2)='CHANGE'
        CO(3)='DATA'
        CO(4)='ROTATE'
        CO(5)='BY'
        CO(6)=CHAR
        NTCO=6
        noco=1
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)

        DO noopti=1,NTOPTI
          USER(OS_PAOPTI+noopti-1)=XC(noopti)
        ENDDO

        WRITE(CHAR,'(D12.4,'','',D12.4,'','',D12.4)')XC(1),XC(2),XC(3)
        CO(1)='FEM'
        CO(2)='CHANGE'
        CO(3)='DATA'
        CO(4)='ROTATE'
        CO(5)='BY'
        CO(6)=CHAR
        NTCO=6
        noco=1
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)

        CO(1)='FEM'
        CO(2)='EVALUATE'
        CO(3)='RESIDUALS'
        CO(4)='WRT'
        CO(5)='MOMENTS'
        CO(6)='NOVIEW'
        NTCO=6
        noco=1
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)

        IF(MODE.EQ.0.OR.MODE.EQ.2) THEN
          DO nores=1,NTRES
            FC(nores)=USER(OS_RESID+nores-1)
          ENDDO
        ENDIF
C No jacobian values available
C        IF(MODE.EQ.1.OR.MODE.EQ.2) THEN
C         DO noopti=1,NTOPTI
C           DO nores=1,NTRES
C             FJAC(nores,noopti)=USER(OS_RESJAC+nores-1+
C     '          (noopti-1)*NREM)
C           ENDDO
C         ENDDO
C        ENDIF







      ELSE IF(KTYP26.EQ.2.AND.KTYP27.EQ.12) THEN !activation times
C **    Objective function is the norm of the difference between
C **    measured and computed torso potentials.
        DO noco=1,13
          NTCOQU(noco)=0
        ENDDO
        DO noopti=1,NTOPTI
          USER(OS_PAOPTI+noopti-1)=XC(noopti)
        ENDDO

        CO(1)='FEM'
        CO(2)='UPDATE'
        CO(3)='OPTIMISATION'
        NTCO=3
        noco=1
C***    !Puts current activation sequence into YP
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)

        CO(1)='FEM'
        CO(2)='APPLY'
        CO(3)='TRANSFER'
        CO(4)='ACTIVATION'
        CO(5)='OUTARRAY'
        CO(6)='TSTART'
        WRITE(CO(7),'(D12.3)') ACTN_MIN(2)
        CO(8)='TEND'
        WRITE(CO(9),'(D12.3)') ACTN_MAX(2)
        NTCO=9

       noco=1
C***    !Puts estimate of torso potentials based on current activation
C***    !sequence into PHI_H array.
       CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)

        CO(1)='FEM'
        CO(2)='EVALUATE'
        CO(3)='RESIDUALS'
        CO(4)='WRT'
        CO(5)='POTENTIAL'
        CO(6)='NOVIEW'
        IF(MODE.EQ.0) THEN
          CO(7)='FUNCTION'
        ELSEIF(MODE.EQ.1) THEN
          CO(7)='JACOBIAN'
        ELSEIF(MODE.EQ.2) THEN
          CO(7)='BOTH'
        ELSE
          ERROR='Invalid mode'
          GOTO 9999
        ENDIF
        CO(8)='TSTART'
        WRITE(CO(9),'(D12.3)') ACTN_MIN(2)
        CO(10)='TEND'
        WRITE(CO(11),'(D12.3)') ACTN_MAX(2)
        NTCO=11
       noco=1
C***    !Computes the differences between potentials in the PHI and
C***    !PHI_H arrays.
       CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)

       IF(MODE.EQ.0.OR.MODE.EQ.2) THEN
          DO nores=1,NTRES
            FC(nores)=USER(OS_RESID+nores-1)
          ENDDO
        ENDIF
        IF(MODE.EQ.1.OR.MODE.EQ.2) THEN
         DO noopti=1,NTOPTI
           DO nores=1,NTRES
             FJAC(nores,noopti)=USER(OS_RESJAC+nores-1+
     '          (noopti-1)*NREM)
           ENDDO
         ENDDO
        ENDIF




      ELSE IF(KTYP26.EQ.2.AND.KTYP27.EQ.13) THEN !dipole
        IF(DOP) THEN
          WRITE(*,*) 'In FUNCT2...'
          WRITE(*,*) '  Current time is ...',USER(OS_WORK)
          WRITE(*,*) '  Dipole region   ...',IUSER(OS_IWORK)
          WRITE(*,*) '  Current freq is ...',TRSF_FREQUENCY
          WRITE(*,*)
        ENDIF
        
C*** Write the variables into strings for passing into FEM
        WRITE(OPTI_REGION,'(I2)') IUSER(OS_IWORK) !nr being optimised
C        WRITE(OPTI_REGION,'(I2)') IUSER(OS_IWORK) !coupled problem

C LKC 13-MAR-2003 New option for diff squid configurations
        WRITE(SQUID_CONFIG,'(I1)') IUSER(OS_IWORK+1) !squid config

        
        WRITE(TIME,'(F8.2)') USER(OS_WORK)
        WRITE(FREQUENCY,'(F8.2)') TRSF_FREQUENCY
        
C LKC 9-JUL-2003 Allow baseline calculations 
        WRITE(BASELINE,'(F5.1)') USER(OS_WORK+1) !Baseline

        DO noopti=1,NTOPTI
          USER(OS_PAOPTI+noopti-1)=XC(noopti) !optim vars
        ENDDO

        IF(DOP) THEN
          WRITE(*,'(''** FUNCT2: Dipole Cent : '',3F15.9)')
     &      XC(1),XC(2),XC(3)
          WRITE(*,'(''** FUNCT2: Dipole Dirn : '',3F15.9)')
     &      XC(4),XC(5),XC(6)
        ENDIF

C*** Puts optimised dipole parameters into cm arrays        
        CO(1)='FEM'
        CO(2)='UPDATE'
        CO(3)='OPTIMISATION'
        CO(4)='TSTART'
        CO(5)=TIME
        CO(6)='TEND'
        CO(7)=TIME
        CO(8)='REGION'
        CO(9)=OPTI_REGION ! this is the dipole region number ....
        NTCO=9
        noco=1

C        WRITE(*,*) '  *** FUNCT2 -- upopti',TIME
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)

C -- FEM update GD for new dipole source
C -- FEM solve for dipole

        CO(1)='FEM'
        CO(2)='UPDATE'
        CO(3)='SOURCE'
        CO(4)='REGION'
        CO(5)=OPTI_REGION
        CO(6)='TIME'
        CO(7)=TIME
        NTCO=7
        noco=1
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '    ERROR,*9999)

        IF(IS_COUPLED(nx)) THEN !note nx=1 is temporary
          CO(1)='FEM'
          CO(2)='SOLVE'
          CO(3)='COUPLED'
          CO(4)='AT'
          CO(5)=TIME
          NTCO=5
          noco=1
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)
        ELSE
          CO(1)='FEM'
          CO(2)='SOLVE'
          CO(3)='AT'
          CO(4)=TIME
          NTCO=4
          noco=1
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)
        ENDIF


C*** Computes estimate of MAGNETIC field based on current dipole
C*** parameters. The estimated solution is stored in nss 2

C LKC 10-DEC-2002 read in the positions of the magnetic sensors
        IF(KTYP27B.EQ.3) THEN

          CO(1)='FEM'
          CO(2)='DEFINE'
          CO(3)='DATA'
          CO(4)='NOVIEW'
          NTCO=4
          noco=1

          COQU(3,1)='R'
          COQU(3,2)=DATASET2     !Name of input file
          COQU(3,3)='EXAMPLE'     !Name of input file
          NTCOQU(1)=0
          NTCOQU(2)=0
          NTCOQU(3)=3

          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)

        ENDIF


        IF(KTYP27B.EQ.1.OR.KTYP27B.EQ.3) THEN
          CO(1)='FEM'
          CO(2)='EVALUATE'
          CO(3)='SOLUTION'
          CO(4)='MAGNETIC'
          CO(5)='FROM'
          CO(6)='MFI'
          CO(7)='NSS'
          CO(8)='2'
          CO(9)='NOVIEW'
          CO(10)='TIME'
          CO(11)=TIME
          CO(12)='FREQUENCY'
          CO(13)=FREQUENCY
          CO(14)='REGION'
          CO(15)='1' ! this is the potential field region number
          CO(16)='SQUID_CONFIG' ! squid configuration
          CO(17)=SQUID_CONFIG 
          CO(18)='BASELINE' ! squid configuration
          CO(19)=BASELINE 
          
          NTCO=19
          noco=1
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)
        ENDIF


C LKC 10-DEC-2002 read in the positions of the potential sensors
        IF(KTYP27B.EQ.3) THEN

          CO(1)='FEM'
          CO(2)='DEFINE'
          CO(3)='DATA'
          CO(4)='NOVIEW'
          NTCO=4
          noco=1
          COQU(3,1)='R'
          COQU(3,2)=DATASET1     !Name of input file
          COQU(3,3)='EXAMPLE'     !Name of input file
          NTCOQU(1)=0
          NTCOQU(2)=0
          NTCOQU(3)=3
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)


          CO(1)='FEM'
          CO(2)='DEFINE'
          CO(3)='XI'
          NTCO=3
          noco=1
          COQU(3,1)='R'
          COQU(3,2)=DATASET1     !Name of input file
          COQU(3,3)='EXAMPLE'
          NTCOQU(1)=0
          NTCOQU(2)=0
          NTCOQU(3)=3
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)

        ENDIF


C*** Computes estimate of POTENTIAL field based on current dipole
C*** parameters
        IF(KTYP27B.EQ.2.OR.KTYP27B.EQ.3) THEN

C -- FEM eval PHI directly from YP
C!!! Note this is defaulted to region 1 data points
          CO(1)='FEM'
          CO(2)='UP'
          CO(3)='PHI'
          CO(4)='TSTART'
          CO(5)=TIME
          CO(6)='TEND'
          CO(7)=TIME
          CO(8)='FROM'
          CO(9)='DATA'
          CO(10)='PHI_H'
          NTCO=10
          noco=1
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)

        ENDIF

C*** Computes the differences between magnetic/potential fields
C*** nss2 (current estimate) and nss1 (measured solution) in MFI arrays
        CO(1)='FEM'
        CO(2)='EVALUATE'
        CO(3)='RESIDUALS'
        CO(4)='WRT'
        CO(5)='MAGNETIC'
        CO(6)='LIST'
        CO(7)='0'
        CO(8)='TSTART'
        CO(9)=TIME
        CO(10)='TEND'
        CO(11)=TIME
        NTCO=11
        noco=1

C        WRITE(*,*) '  *** FUNCT2 evresid -- ',TIME
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)

        DO nores=1,NTRES
          FC(nores)=USER(OS_RESID+nores-1)
          IF(DOP) THEN
            WRITE(OP_STRING,'('' RESID '',I5,D12.3)') nores,FC(nores)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO ! nores

      ELSE IF(KTYP26.EQ.3) THEN ! Micro-structure
        DO noopti=1,NTOPTI
          USER(OS_PAOPTI+noopti-1)=XC(noopti) !optimisation vars
        ENDDO
        IF(ITYP2(nr,nx).EQ.2) THEN !Finite elasticity problem
          CO(1)='FEM'
          CO(2)='UPDATE'
          CO(3)='FIELD'
          CO(4)='OPTIMISE'
          NTCO=4
          noco=1
          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP       CRITICAL(FUNCT2_27)
            WRITE(OP_STRING,'('' fem update field optimise..'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP       END CRITICAL(FUNCT2_27)
          ENDIF
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)

          CALL STRING_TRIM(COM_FILE,IBEG,IEND)
          CO(1)='FEM'
          CO(2)='READ'
          CO(3)='COMMAND'
          NTCO=3
          noco=1
          COQU(3,1)=COM_FILE(IBEG:IEND)
          NTCOQU(1)=0
          NTCOQU(2)=0
          NTCOQU(3)=1
          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP       CRITICAL(FUNCT2_28)
            WRITE(OP_STRING,'('' fem read com;filename..'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP       END CRITICAL(FUNCT2_28)
          ENDIF
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)

        ELSE
          CALL ASSERT(.FALSE.,'Problem type not implemented',
     '      ERROR,*9999)
        ENDIF

        IF(MODE.EQ.0.OR.MODE.EQ.2) THEN
C ***     evaluate residuals
          IF(KTYP27.LT.4) THEN
            CO(1)='FEM'
            CO(2)='EVALUATE'
            CO(3)='RESIDUALS'
            CO(4)='WRT'
            CO(5)='STRESS'
            CO(6)='FIBRE'
            CO(7)='NOVIEW'
            NTCO=8
          ELSE
            CO(1)='FEM'
            CO(2)='EVALUATE'
            CO(3)='RESIDUALS'
            CO(4)='WRT'
            CO(5)='STRAIN'
            CO(6)='FIBRE'
            CO(7)='NOVIEW'
            NTCO=8
          ENDIF
          noco=1
          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP       CRITICAL(FUNCT2_29)
            WRITE(OP_STRING,'('' fem evaluate residuals ..'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP       END CRITICAL(FUNCT2_29)
          ENDIF
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)
          DO nores=1,NTRES
            FC(nores)=USER(OS_RESID+nores-1)
            IF(DOP) THEN
              WRITE(OP_STRING,'('' RESID '',I5,D12.3)')
     '          nores,FC(nores)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO ! nores
C LKC 20-MAR-2002 Unused?
C          OPTIMISING=.FALSE.
        ENDIF

        IF(IPDLEV.GT.0.AND.(MODE.EQ.1.OR.MODE.EQ.2)) THEN
C ***     evaluate first derivatives of residuals wrt opt params
          CALL ASSERT(.FALSE.,'Not Implemented',ERROR,*9999)
        ENDIF

      ELSE IF(KTYP26.EQ.4) THEN ! Holmes constants
        DO noopti=1,NTOPTI
          USER(OS_PAOPTI+noopti-1)=XC(noopti) !optimisation vars
        ENDDO
        IF(ITYP2(nr,nx).EQ.2) THEN !Finite elasticity problem
          CO(1)='FEM'
          CO(2)='UPDATE'
          CO(3)='MATERIAL'
          CO(4)='OPTIMISE'
          NTCO=4
          noco=1
          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP       CRITICAL(FUNCT2_27)
            WRITE(OP_STRING,'('' fem update field optimise..'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP       END CRITICAL(FUNCT2_27)
          ENDIF
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)

          CALL STRING_TRIM(COM_FILE,IBEG,IEND)
          CO(1)='FEM'
          CO(2)='READ'
          CO(3)='COMMAND'
          NTCO=3
          noco=1
          COQU(3,1)=COM_FILE(IBEG:IEND)
          NTCOQU(1)=0
          NTCOQU(2)=0
          NTCOQU(3)=1
          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP       CRITICAL(FUNCT2_28)
            WRITE(OP_STRING,'('' fem read com;filename..'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP       END CRITICAL(FUNCT2_28)
          ENDIF
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)

        ELSE
          CALL ASSERT(.FALSE.,'Problem type not implemented',
     '      ERROR,*9999)
        ENDIF

        IF(MODE.EQ.0.OR.MODE.EQ.2) THEN
C ***     evaluate residuals
          CO(1)='FEM'
          CO(2)='EVALUATE'
          CO(3)='RESIDUALS'
          CO(4)='WRT'
          CO(5)='PRESSURE'
          CO(7)='NOVIEW'
          NTCO=8
          noco=1
          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP       CRITICAL(FUNCT2_29)
            WRITE(OP_STRING,'('' fem evaluate residuals ..'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP       END CRITICAL(FUNCT2_29)
          ENDIF
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)
          DO nores=1,NTRES
            FC(nores)=USER(OS_RESID+nores-1)
            IF(DOP) THEN
              WRITE(OP_STRING,'('' RESID '',I5,D12.3)')
     '          nores,FC(nores)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO ! nores
C LKC 20-MAR-2002 Unused?
C          OPTIMISING=.FALSE.
        ENDIF

        IF(IPDLEV.GT.0.AND.(MODE.EQ.1.OR.MODE.EQ.2)) THEN
C ***     evaluate first derivatives of residuals wrt opt params
          CALL ASSERT(.FALSE.,'Not Implemented',ERROR,*9999)
        ENDIF
      ELSE IF(KTYP26.EQ.1.AND.KTYP27.EQ.5) THEN
C **    Mat'l estimation with residuals being data error
        DO noopti=1,NTOPTI
          USER(OS_PAOPTI+noopti-1)=XC(noopti)
C          WRITE(OP_STRING,'('' XC('',I5,'') = '',D12.3)')
C     '         noopti,XC(noopti)
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
        CO(1)='FEM'
        CO(2)='UPDATE'
        CO(3)='OPTIMISATION'
        NTCO=3
        noco=1
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '    ERROR,*9999)

        CALL STRING_TRIM(COM_FILE,IBEG,IEND)
        CO(1)='FEM'
        CO(2)='READ'
        CO(3)='COMMAND'
        NTCO=3
        noco=1
        COQU(3,1)=COM_FILE(IBEG:IEND)
        NTCOQU(1)=0
        NTCOQU(2)=0
        NTCOQU(3)=1
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '    ERROR,*9999)
        DO nores=1,NTRES
          FC(nores)=USER(OS_RESID+nores-1)
C          WRITE(OP_STRING,'(''FC('',I5,'')='',D10.4)') 
C     '      nores,FC(nores)
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO

      ELSE IF(KTYP26.EQ.1.AND.KTYP27.EQ.6) THEN
C **    Mat'l estimation with residuals being data error
        DO noopti=1,NTOPTI
          USER(OS_PAOPTI+noopti-1)=XC(noopti)
C          WRITE(OP_STRING,'('' XC('',I5,'') = '',D12.3)')
C     '         noopti,XC(noopti)
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
        CO(1)='FEM'
        CO(2)='UPDATE'
        CO(3)='OPTIMISATION'
        NTCO=3
        noco=1
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '    ERROR,*9999)

        CALL STRING_TRIM(COM_FILE,IBEG,IEND)
        CO(1)='FEM'
        CO(2)='READ'
        CO(3)='COMMAND'
        NTCO=3
        noco=1
        COQU(3,1)=COM_FILE(IBEG:IEND)
        NTCOQU(1)=0
        NTCOQU(2)=0
        NTCOQU(3)=1
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '    ERROR,*9999)

C ***   evaluate residuals
        CO(1)='FEM'
        CO(2)='EVALUATE'
        CO(3)='RESIDUALS'
        CO(4)='WRT'
        CO(5)='REACTION'
        CO(7)='NOVIEW'
        NTCO=8
        noco=1
        IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C            Critical section is not essential.
CC$OMP     CRITICAL(FUNCT2_29)
          WRITE(OP_STRING,'('' fem evaluate residuals ..'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP     END CRITICAL(FUNCT2_29)
        ENDIF
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '    ERROR,*9999)
        DO nores=1,NTRES
          FC(nores)=USER(OS_RESID+nores-1)
          IF(DOP) THEN
            WRITE(OP_STRING,'('' RESID '',I5,D12.3)')
     '        nores,FC(nores)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO ! nores
      ENDIF

      CALL EXITS('FUNCT2')
      RETURN
 9999 CALL ERRORS('FUNCT2',ERROR)
      CALL STRING_TRIM(ERROR,IBEG,IEND)
      WRITE(OP_STRING,'(A)') ' '//ERROR(IBEG:IEND)
      CALL WRITES(IOER,OP_STRING,ERROR,*9998)
      ERROR(1:)=' '
      MODE=-9999
 9998 CALL EXITS('FUNCT2')
      RETURN
      END


