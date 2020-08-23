      SUBROUTINE MARCH6(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,IBT,IDO,
     '  INP,ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,ISR_GKK,ISR_GQ,LGE,NBH,NBJ,
     '  NDET,NDIPOLES,NEELEM,NENP,NGAP,NHE,NHP,NHQ,NKH,NKHE,NKJE,
     '  NLL,NONY,NORD,NPB,NPF,NP_INTERFACE,NPNE,NPNODE,NPNY,NQNY,nr,
     '  NRE,NRLIST,NRLIST2,NVHE,NVHP,NVJE,NW,nx,NYNE,NYNO,NYNP,
     '  NYNR,NYQNR,CE,CGE,CONY,CP,CURVCORRECT,CYNO,DET,DIPOLE_CEN,
     '  DIPOLE_DIR,DL,DRDN,DRDNO,GD,
     '  GK,GKK,GQ,GR,GRR,PG,RAD,RD,RG,SE,WG,XA,XE,XG,XG1,XIG,XN,
     '  XN_GRAD,XO,XP,XR,XR_GRAD,YG,YP,YQ,YQS,ZA,ZP,FIX,ERROR,*)
C     SMAR009 18/01/99 removed ,NWP from list
C#### Subroutine: MARCH6
C###  Description:
C###    MARCH6 solves quasi-static problems.

C-**** Time stepping is governed by KTYP23=1,2 for fixed or automatic
C**** time steps. For KTYP23=2 the error is calculated & time step
C**** adjusted if KTYP23=2.
C**** On exit YP(ny,1) is solution vector.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'coup00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'marc00.cmn'
      INCLUDE 'quas00.cmn'
C     SMAR009 removed 16/12/98 INCLUDE 'cmiss$reference:solv00.cmn'
      INCLUDE 'time02.cmn'


!     Parameter List
      INTEGER DIPOLE_CEN_NTIME(NDIPOLEM,NRM,NXM),
     '  DIPOLE_DIR_NTIME(NDIPOLEM,NRM,NXM),IBT(3,NIM,NBFM),
     '  IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),ISC_GK(NISC_GKM),
     '  ISC_GKK(NISC_GKKM),ISC_GQ(NISC_GQM),ISR_GK(NISR_GKM),
     '  ISR_GKK(NISR_GKKM),ISR_GQ(NISR_GQM),LGE(NHM*NSM,NRCM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NDET(NBFM,0:NNM),
     '  NDIPOLES(NRM,NXM),NEELEM(0:NE_R_M,0:NRM),NENP(NPM,0:NEPM,0:NRM),
     '  NGAP(NIM,NBM),NHE(NEM),NHP(NPM,0:NRM),NHQ(NRM),
     '  NKH(NHM,NPM,NCM,0:NRM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),NONY(0:NOYM,NYM,NRCM,0:NRM),
     '  NORD(5,NE_R_M),NPB(0:NP_R_M,5),NPF(9,NFM),
     '  NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM),NQNY(2,NYQM,0:NRCM),
     '  nr,NRE(NEM),NRLIST(0:NRM),NRLIST2(0:NRM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),
     '  nx,NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNO(0:NYOM,NOOPM,NRCM,0:NRM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM),NYQNR(0:NYQM,0:NRCM,NCM,0:NRM)
C     SMAR009 18/01/99 removed NWP(NPM,2),
      REAL*8 CE(NMM,NEM),CGE(NMM,NGM,NEM),CONY(0:NOYM,NYM,NRCM,0:NRM),
     '  CP(NMM,NPM),CURVCORRECT(2,2,NNM,NEM),
     '  CYNO(0:NYOM,NOOPM,NRCM,0:NRM),
     '  DET(NBFM,0:NNM,NGM,6),DIPOLE_CEN(4,0:NDIPTIMM,NDIPOLEM,NRM,NXM),
     '  DIPOLE_DIR(4,0:NDIPTIMM,NDIPOLEM,NRM,NXM),
     '  DL(3,NLM),DRDN(NGM),DRDNO(NGM,NKM),
     '  GD(NZ_GD_M),GK(NZ_GK_M),GKK(NZ_GKK_M),
     '  GQ(NZ_GQ_M),GR(NYROWM),GRR(NOM),PG(NSM,NUM,NGM,NBM),
     '  RAD(NGM),RD(NGM),RG(NGM),SE(NSM,NBFM,NEM),
     '  WG(NGM,NBM),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XG(NJM,NUM),XG1(NJM,NUM,NGM),XIG(NIM,NGM,NBM),
     '  XN(NJM,NGM),XN_GRAD(NJM,NGM),XO(NOM),XP(NKM,NVM,NJM,NPM),
     '  XR(NJM,NGM),XR_GRAD(NJM,NGM),YG(NIYGM,NGM,NEM),YP(NYM,NIYM),
     '  YQ(NYQM,NIQM,NAM),YQS(NIQSM,NQM),ZA(NAM,NHM,NCM,NEM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM)
!     Local Variables
      INTEGER CERROR(50),ERR,nc,na,NIQLIST(0:1),NIQSLIST(0:1),
     '  NIYLIST(0:16),no_nrlist,no_nynr,NOTIMES,np,
     '  nrr,NSTEP,NUMTIMEDATA,
     '  NUMTIMEDATA1,ny
      REAL AVETIME,ELAPSED_TIME,TIME_START1(1),TIME_START2(1),
     '  TIME_STOP(1)
      REAL*8 QDT,T,TIME,YPMAX(16),YPMIN(16)
      CHARACTER FILEFORMAT*6
      LOGICAL ENDFILE,ISBINFILEOPEN,
     '  FIRST,FIRST_SOLVE,OUTPUT,
     '  UPDATE_GLOB_MATRIX,UPDATE_GLOB_SOURCE,UPDATE_GLOB_VECTOR,
     '  UPDATE_SOLU_MATRIX,UPDATE_SOLU_SOURCE,UPDATE_SOLU_VECTOR,
     '  YPDATA,YQDATA,YQSDATA
      SAVE NSTEP

      CALL ENTERS('MARCH6',*9999)

      CALL ASSERT(NIYM.GE.16,'>>Increase NIYM to at least 16',
     '  ERROR,*9999)
      CALL CPU_TIMER(CPU_USER,TIME_START1)
      CALL CPU_TIMER(CPU_USER,TIME_START2)

      QDT=TINCR

      IF(BINTIMEFILE.GT.0) THEN
        FILEFORMAT='BINARY'
      ELSE
        FILEFORMAT='ASCII'
      ENDIF

      YPDATA=.TRUE.
      YQDATA=.FALSE.
      YQSDATA=.FALSE.
      IF(.NOT.RESTART) THEN !perform initial tasks

        FIRST=.TRUE.
        FIRST_SOLVE=.TRUE.
        T=T0
        NSTEP=0

C*** Put initial conditions YP(ny,3) into current solution YP(ny,1) &
C*** previous solution YP(ny,8)

        DO nc=1,NCT(nr,nx)
          DO no_nynr=1,NYNR(0,0,nc,nr)
            ny=NYNR(no_nynr,0,nc,nr)
            IF(NPNY(0,ny,0).EQ.1) THEN
              np=NPNY(4,ny,0)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
            ENDIF
            IF(.NOT.FIX(ny,1)) YP(ny,1)=YP(ny,3)
            YP(ny,8)=YP(ny,1)
          ENDDO !no_nynr
        ENDDO !nc

        NIYLIST(0)=1
        NIYLIST(1)=1
        NIQLIST(0)=0
        NIQLIST(1)=0
        NIQSLIST(0)=0
C        na=1

        CALL IOHIST(IOFILE1,na,NHQ,NIQLIST,NIQSLIST,NIYLIST,NPNY,
     '    NQNY,NRLIST,NRLIST2,NUMTIMEDATA,nx,NYNR,NYQNR,TIME,YP,
     '    YPMAX,YPMIN,YQ,YQS,'WRITE',FILEFORMAT,QUASI_OUTHISTFILE,
     '    'OPEN',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,ERROR,*9999)

      ELSE !IF(RESTART) THEN
        T=T_RESTART(nx)
        FIRST=.FALSE.
        FIRST_SOLVE=.FALSE.

        CALL IOHIST(IOFILE1,na,NHQ,NIQLIST,NIQSLIST,NIYLIST,NPNY,
     '    NQNY,NRLIST,NRLIST2,NUMTIMEDATA,nx,NYNR,NYQNR,TIME,YP,
     '    YPMAX,YPMIN,YQ,YQS,'WRITE',FILEFORMAT,QUASI_OUTHISTFILE,
     '    'APPEND',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,ERROR,*9999)

      ENDIF !.NOT.restart

C***  Main time loop

      IF(QUASIPROB.EQ.1) THEN
        UPDATE_GLOB_SOURCE=.TRUE.
        UPDATE_GLOB_VECTOR=.FALSE.
        UPDATE_GLOB_MATRIX=.FALSE.
      ELSE IF(QUASIPROB.EQ.2) THEN
        UPDATE_GLOB_SOURCE=.TRUE.
        UPDATE_GLOB_VECTOR=.TRUE.
        UPDATE_GLOB_MATRIX=.FALSE.

C***    Open up the initial condition history file

        NIYLIST(0)=1
        NIYLIST(1)=1
        NRLIST(0)=1
        NRLIST(1)=nr
        NRLIST2(0)=1
        NRLIST2(1)=nr
        na=1

        CALL IOHIST(IOFILE2,na,NHQ,NIQLIST,NIQSLIST,NIYLIST,NPNY,
     '    NQNY,NRLIST,NRLIST2,NUMTIMEDATA,nx,NYNR,NYQNR,TIME,YP,
     '    YPMAX,YPMIN,YQ,YQS,'READ',FILEFORMAT,QUASI_INHISTFILE,
     '    'OPEN',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,ERROR,*9999)

      ELSE IF(QUASIPROB.EQ.3) THEN
        UPDATE_GLOB_SOURCE=.TRUE.
        UPDATE_GLOB_VECTOR=.TRUE.
        UPDATE_GLOB_MATRIX=.TRUE.
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
      AVETIME=0.0
      NOTIMES=0
      IF(OUTPUT_SOLTIMES) THEN
        WRITE(OP_STRING,'(/'' CPU time for setup ='',D12.4,'' s'')')
     '    ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

 10   CONTINUE !loop start

        CALL CPU_TIMER(CPU_USER,TIME_START2)

        NSTEP=NSTEP+1

        IF(IWRIT1(nr,nx).NE.0) THEN
          IF(MOD(NSTEP,IWRIT1(nr,nx)).EQ.0) THEN
            OUTPUT=.TRUE.
          ELSE
            OUTPUT=.FALSE.
          ENDIF
        ELSE
          OUTPUT=.FALSE.
        ENDIF

C*** The global matrices are assumed to be assembled already so solve
C*** the problem for the current time step.

        IF(FIRST) THEN
          UPDATE_SOLU_SOURCE=.TRUE.
          UPDATE_SOLU_VECTOR=.TRUE.
          UPDATE_SOLU_MATRIX=.TRUE.
          FIRST=.FALSE.
        ELSE
          IF(QUASIPROB.EQ.1) THEN
            UPDATE_SOLU_SOURCE=.TRUE.
            UPDATE_SOLU_VECTOR=.FALSE.
            UPDATE_SOLU_MATRIX=.FALSE.
          ELSE IF(QUASIPROB.EQ.2) THEN
            UPDATE_SOLU_SOURCE=.TRUE.
            UPDATE_SOLU_VECTOR=.TRUE.
            UPDATE_SOLU_MATRIX=.FALSE.
          ELSE IF(QUASIPROB.EQ.3) THEN
            UPDATE_SOLU_SOURCE=.TRUE.
            UPDATE_SOLU_VECTOR=.TRUE.
            UPDATE_SOLU_MATRIX=.TRUE.
          ENDIF
        ENDIF

        IF(QUASIPROB.EQ.2) THEN !Time dependent boundary values

C***      Read in the value of YP for the time step
          CALL IOHIST(IOFILE2,na,NHQ,NIQLIST,NIQSLIST,NIYLIST,NPNY,
     '      NQNY,NRLIST,NRLIST2,NUMTIMEDATA1,nx,NYNR,NYQNR,TIME,YP,
     '      YPMAX,YPMIN,YQ,YQS,'READ',FILEFORMAT,QUASI_INHISTFILE,
     '      'TIME_DATA',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,
     '      ERROR,*9999)

        ENDIF

        IF(nr.EQ.0) THEN !Coupled

          CALL SOLVE9(ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,ISR_GKK,ISR_GQ,
     '      LGE,NBH,NENP,NHE,NDIPOLES,NONY(0,1,1,nr),NP_INTERFACE,
     '      NPNE,NPNY,NRE,NVHE,nx,NYNE,NYNO(0,1,1,nr),NYNP,
     '      NYNR(0,0,1,nr),CONY(0,1,1,nr),CYNO(0,1,1,nr),GD,GK,GKK,
     '      GQ,GR,GRR,XO,YP,FIRST_SOLVE,FIX,UPDATE_SOLU_MATRIX,
     '      UPDATE_SOLU_SOURCE,UPDATE_SOLU_VECTOR,ERROR,*9999)

        ELSE IF(ITYP4(nr,nx).EQ.1) THEN !FEM

          CALL SOLVE1(ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,ISR_GKK,ISR_GQ,
     '      LGE,NBH,NENP,NHE,NONY(0,1,1,nr),NP_INTERFACE,NPNE,NPNY,
     '      nr,NRE,NVHE,nx,NYNE,NYNO(0,1,1,nr),NYNP,
     '      NYNR(0,0,1,nr),CONY(0,1,1,nr),CYNO(0,1,1,nr),GK,GKK,GQ,
     '      GR,GRR,XO,YP,FIRST_SOLVE,FIX,UPDATE_SOLU_MATRIX,
     '      UPDATE_SOLU_VECTOR,ERROR,*9999)

        ELSE IF(ITYP4(nr,nx).EQ.2) THEN !BEM

          CALL SOLVE2(ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,ISR_GKK,ISR_GQ,
     '      LGE,NBH,NENP,NHE,NDIPOLES,NONY(0,1,1,nr),NP_INTERFACE,
     '      NPNE,NPNY,nr,NRE,NVHE,nx,NYNE,NYNO(0,1,1,nr),NYNP,
     '      NYNR(0,0,1,nr),CONY(0,1,1,nr),CYNO(0,1,1,nr),GD,GK,GKK,
     '      GQ,GRR,XO,YP,FIRST_SOLVE,FIX,UPDATE_SOLU_MATRIX,
     '      UPDATE_SOLU_SOURCE,UPDATE_SOLU_VECTOR,ERROR,*9999)

        ENDIF

        IF(FIRST_SOLVE) FIRST_SOLVE=.FALSE.

C*** Adjust any time varying parameters.

C*** Write output to screen

        IF(OUTPUT) THEN
          WRITE(OP_STRING,'(/'' Solution at time T='',D11.4,'
     '      //''' with DT='',D11.4)') T,QDT
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          DO nc=1,NCT(nr,nx)
C LKC 15-APR-2002 doesn't work if entering nr=0 which is what is
C   passed in. Using the first region for now....
C            CALL YPZP(1,NBH,NEELEM,NHE,NHP(1,nr),NKH(1,1,1,nr),NPNODE,
C     '        nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
C            CALL ZPOP(4,NBH,nc,NEELEM,NHP(1,nr),NKH(1,1,1,nr),
C     '        NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,
C     '        YP(1,1),ZA,ZP,FIX(1,1),ERROR,*9999)
            CALL YPZP(1,NBH,NEELEM,NHE,NHP(1,nr),NKH(1,1,1,nr),NPNODE,
     '        NRLIST(1),NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZA,ZP,
     '        ERROR,*9999)
            CALL ZPOP(4,NBH,nc,NEELEM,NHP(1,nr),NKH(1,1,1,nr),
     '        NPNODE,NRLIST(1),NVHP(1,1,1,nr),nx,NYNE,NYNP,
     '        YP(1,1),ZA,ZP,FIX(1,1),ERROR,*9999)
          ENDDO !nc
        ENDIF

        TIME=T
        YPDATA=.TRUE.
        YQDATA=.FALSE.
        na=1
        CALL IOHIST(IOFILE1,na,NHQ,NIQLIST,NIQSLIST,NIYLIST,NPNY,
     '    NQNY,NRLIST,NRLIST2,NUMTIMEDATA,nx,NYNR,NYQNR,TIME,YP,
     '    YPMAX,YPMIN,YQ,YQS,'WRITE',FILEFORMAT,QUASI_OUTHISTFILE,
     '    'TIME_DATA',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,ERROR,*9999)

C*** Adjust time stepping parameters if required

C        IF(KTYP23.EQ.2) THEN
C          KFAC=2
C          ERR=0.0d0
C          DO no_nynr=1,NYNR(0,0,1,nr) !Loop over global variables of GK
C            ny=NYNR(no_nynr,0,1,nr) !global variable #
C            YP(ny,8)=(YP(ny,13)-YP(ny,14))*QDT**1/KFAC
C            ERR=ERR+YP(ny,8)**2
C          ENDDO
C          ERR=DSQRT(ERR)
C          SNORM=0.0D0
C          DO no_nynr=1,NYNR(0,0,1,nr) !Loop over global variables of GK
C            ny=NYNR(no_nynr,0,1,nr) !global variable #
C            SNORM=SNORM+YP(ny,1)**2
C          ENDDO
C          SNORM=DSQRT(SNORM)
C          TOL=SNORM
C          IF(ERR.GT.TOL) THEN
C            QDT=QDT/2.0d0
C            INCR=0
C          ELSE IF(ERR.LT.TOL/2.0d0) THEN
C            INCR=INCR+1
C            IF(INCR.GT.1) THEN
C              QDT=QDT*1.250d0
C              INCR=0
C            ENDIF
C          ENDIF
C        ENDIF

C*** Transfer current solution to previous solution and increment time

        T=T+QDT
        DO no_nynr=1,NYNR(0,0,1,nr) !loop over global variables of GK
          ny=NYNR(no_nynr,0,1,nr) !global variable #
          IF(NPNY(0,ny,0).EQ.1) THEN
            np=NPNY(4,ny,0)
C GMH 8/1/97 Update cmgui link
            CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
          ENDIF
          YP(ny,8)=YP(ny,1)  !is solution vector at time T
        ENDDO

        IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP     CRITICAL(MARCH6_1)
          WRITE(OP_STRING,'('' TIMES : CURRENT ='',D10.4,'
     '      //''' DELTA='',D10.4,'' INIT='',D10.4,'' FINAL='',D10.4)')
     '      T,QDT,T0,T1
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP     END CRITICAL(MARCH6_1)
        ENDIF

C*** Check for numeric keypad entry interrupt
C KAT 2Mar01: does nothing
C        CALL GETSTR2(ERROR,*9998)

C*** Check if all the current time is less than the finish time and
C*** repeat solution if so.

      IF(T.LE.T1) THEN

C*** Check if the stiffness matrices needs to be (re)calculated and if
C*** so (re)calculate them.

        DO no_nrlist=1,QUASIREGLIST(0)
          nrr=QUASIREGLIST(no_nrlist)

          IF(ITYP4(nrr,nx).EQ.1) THEN !FEM

            CALL ASSEMBLE1(IBT,IDO,INP,ISC_GK,ISC_GKK,ISC_GQ,
     '        ISR_GK,ISR_GKK,ISR_GQ,NBH,NBJ,NEELEM,NHE,NHP(1,nrr),
     '        NKJE,NONY(0,1,1,0),NORD,NPF,NP_INTERFACE,NPNE,
     '        NPNY,nrr,0,NRE,NVHE,NVJE,NW,nx,NYNE,NYNP,
     '        NYNR(0,0,1,nrr),CE,CGE,CONY(0,1,1,0),CP,CURVCORRECT,
     '        GK,GKK,GQ,GR,PG,SE,WG,XA,XP,YG,COUPLED_BEM(nx),
     '        UPDATE_GLOB_MATRIX,UPDATE_GLOB_VECTOR,ERROR,*9999)

          ELSE IF(ITYP4(nrr,nx).EQ.2) THEN !BEM

            CALL ASSEMBLE2(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,
     '        IBT,IDO,INP,ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,ISR_GKK,ISR_GQ,
     '        NBH,NBJ,NDET,NDIPOLES,NEELEM,NENP,NGAP,NHE,NHP(1,nrr),
     '        NKH(1,1,1,nrr),NKHE,NKJE,NLL,NONY(0,1,1,0),NPF,
     '        NP_INTERFACE,NPNE,NPNODE,NPNY,nrr,0,NRE,NVHP(1,1,1,nrr),
     '        NVJE,NW,nx,NYNE,NYNP,NYNR(0,0,1,nrr),CE,CONY(0,1,1,0),
     '        CURVCORRECT,DET,DIPOLE_CEN,DIPOLE_DIR,DL,
     '        GD,GK,GKK,GQ,PG,SE,T,WG,XA,XE,XG,
     '        XIG,XP,UPDATE_GLOB_MATRIX,
     '        UPDATE_GLOB_SOURCE,ERROR,*9999)

          ENDIF
        ENDDO !no_nrlist

        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
        AVETIME=AVETIME+ELAPSED_TIME
        NOTIMES=NOTIMES+1
        IF(OUTPUT_SOLTIMES) THEN
          WRITE(OP_STRING,'(/'' Step ='',I4,'', Time ='',D12.4)')
     '      NSTEP,T
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' CPU time for the time step ='',D12.4,'
     '      //''' s'')') ELAPSED_TIME
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

        GOTO 10
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START1(1)
      IF(OUTPUT_SOLTIMES) THEN
        AVETIME=AVETIME/NOTIMES
        WRITE(OP_STRING,'(/'' Total CPU time for the problem ='',D12.4,'
     '    //''' s'')') ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Total number of time steps = '',I4)')
     '    NOTIMES
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Average CPU time per time step ='',D12.4,'
     '    //''' s'')') AVETIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

C     still have YPDATA=true,YQDATA=false here.
      CALL IOHIST(IOFILE1,na,NHQ,NIQLIST,NIQSLIST,NIYLIST,NPNY,NQNY,
     '  NRLIST,NRLIST2,NUMTIMEDATA,nx,NYNR,NYQNR,TIME,YP,YPMAX,YPMIN,
     '  YQ,YQS,'CLOSE',FILEFORMAT,QUASI_OUTHISTFILE,' ',ENDFILE,.TRUE.,
     '  YPDATA,YQDATA,YQSDATA,ERROR,*9999)

      IF(QUASIPROB.EQ.2) THEN
        CALL IOHIST(IOFILE2,na,NHQ,NIQLIST,NIQSLIST,NIYLIST,NPNY,
     '    NQNY,NRLIST,NRLIST2,NUMTIMEDATA1,nx,NYNR,NYQNR,TIME,YP,
     '    YPMAX,YPMIN,YQ,YQS,'CLOSE',FILEFORMAT,QUASI_INHISTFILE,
     '    ' ',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,ERROR,*9999)
      ENDIF

      CALL EXITS('MARCH6')
      RETURN

C 9998 T_RESTART=T
C      IF(ISBINFILEOPEN(IOFILE1))
C     '  CALL BINARYCLOSEFILE(IOFILE1,ERR,CERROR)
C      IF(ISBINFILEOPEN(IOFILE2))
C     '  CALL BINARYCLOSEFILE(IOFILE2,ERR,CERROR)
C      CALL EXITS('MARCH6')
C      RETURN
 9999 CALL ERRORS('MARCH6',ERROR)
      IF(ISBINFILEOPEN(IOFILE1))
     '  CALL BINARYCLOSEFILE(IOFILE1,ERR,CERROR)
      IF(ISBINFILEOPEN(IOFILE2))
     '  CALL BINARYCLOSEFILE(IOFILE2,ERR,CERROR)
      CALL EXITS('MARCH6')
      RETURN 1
      END

