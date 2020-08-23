      SUBROUTINE MARCH11(IBT,IDO,INP,ISC_GD,ISC_GK,ISC_GKK,ISR_GD,
     &  ISR_GK,ISR_GKK,NBH,NBJ,NEELEM,NELIST,NENP,NHE,NHP,NKJE,
     &  NO_NE,NONY,
     &  NORD,NPF,NPLIST2,NPNE,NPNODE,NPNY,nr,NTIME_POINTS,NVJE,NVJP,nx,
     &  NXI,NYNE,NYNO,NYNP,NYNR,NY_OFST,NZ_ESED,ACINUS,BBM,CE,CG,CGE,
     &  CONY,CYNO,CP,ED,EM,ER,ES,GD,GK,GK2,GKK,GR,GRR,PG,RG,
     &  SE,STACK_ED,STACK_EM,STACK_ES,TIME_VALUES,WG,
     &  XA,XAB,XG,XO,XP,YG,YP,ZE,ZG,FIRST_NX,FIX,ERROR,*)
      
C#### Subroutine: MARCH11
C###  Description:
C###    MARCH11 performs time integration of pulmonary transport
C###    problems.  ITYP3=1 is an advective-diffusive transport equation
C###    for gas mixing in alveolated airways; ITYP3=2 is a coupled
C###    system of temperature and water vapour transfer equations;
C###    ITYP3=3 is a coupled system of pressure and flow equations for
C###    non-Newtonian fluid in small capillaries.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'loc00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lungex00.cmn'
      INCLUDE 'pulm00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'marc00.cmn'
      INCLUDE 'moti00.cmn'
      INCLUDE 'solv00.cmn'
      INCLUDE 'time02.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISC_GD(NISC_GDM),ISC_GK(NISC_GKM),ISC_GKK(NISC_GKKM),
     '  ISR_GD(NISR_GDM),ISR_GK(NISR_GKM),ISR_GKK(NISR_GKKM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M),NELIST(0:NEM),
     '  NENP(NPM,0:NEPM),NHE(NEM),NHP(NPM),
     '  NKJE(NKM,NNM,NJM,NEM),NO_NE(NEM),NONY(0:NOYM,NYM,NRCM),
     '  NORD(5,NE_R_M),NPF(9,NFM),NPLIST2(0:NPM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M),NPNY(0:6,NYM,0:NRCM),nr,
     '  NTIME_POINTS(NTIMEVARSM),NVJE(NNM,NBFM,NJM,NEM),
     &  NVJP(NJM,NPM),nx,
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNO(0:NYOM,NOOPM,NRCM),NY_OFST(NYM),
     &  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM),NZ_ESED(0:12,NE_R_M)
      REAL*8 ACINUS(4,NEM),BBM(2,NEM),CE(NMM,NEM),CG(NMM,NGM),
     '  CGE(NMM,NGM,NEM),CONY(0:NOYM,NYM,NRCM),CP(NMM,NPM),
     &  CYNO(0:NYOM,NOOPM,NRCM),GD(NZ_GD_M),GK(NZ_GK_M),GK2(NZ_GK_M),
     &  GKK(NZ_GKK_M),GR(NYROWM),GRR(NOM),PG(NSM,NUM,NGM,NBM),
     &  SE(NSM,NBFM,NEM),STACK_ED(12,NE_R_M),STACK_EM(12,NE_R_M),
     &  STACK_ES(12,NE_R_M),TIME_VALUES(2,0:NTIMEPOINTSM+1,NTIMEVARSM),
     &  WG(NGM,NBM),XA(NAM,NJM,NEM),XAB(NORM,NEM),XO(NOM),
     &  XP(NKM,NVM,NJM,NPM),YG(NIYGM,NGM,NEM),YP(NYM,NIYM),ZE(NSM,NHM),
     &  ZG(NHM,NUM),ED(NHM*NSM,NHM*NSM),EM(NHM*NSM,NHM*NSM),ER(NHM*NSM),
     &  ES(NHM*NSM,NHM*NSM),RG(NGM),XG(NJM,NUM)
      CHARACTER ERROR*(*)
      LOGICAL FIRST_NX,FIX(NYM,NIYFIXM)
!     Local Variables
      INTEGER N_TM,nb,nbreath,NHST,ntp,ntv,
     '  nxl,ny,SPARSE_P,SPARSE_S
      INTEGER*4 NELIST2_PTR
      REAL*8 gradient,MASS,MASS_EXP,SUM_N2,SUM_WATER,T,water_loss
      LOGICAL CONTINUE,DYNAM1,FIRST_A,GAS,INTERVAL,
     '  LINEAR,SCALE,SMOOTHING,UPDATE_MATRIX,UPDATE_VECTOR,
     '  UPDATE_VELOCITY,update_ny

      CALL ENTERS('MARCH11',*9999)

!       IF(ITYP3(nr,nx).EQ.7) GOTO 998 ! AJS redundant code
      
      IF(ITYP3(nr,nx).LT.3)THEN
        CALL MESH_VOLUME(NBJ,NEELEM,NORD,NPNE,NVJE,NXI,BBM,CE,XP,ERROR,
     &    *9999)
        CURRENT_VOLUME(nx)=CE(4,1)
        
        IF(.NOT.RESTART) INITIAL_VOLUME=CE(4,1) !volume below stem branch

        CALL CALCMASS(NBJ,NEELEM,NH_LOC(1,nx),NORD,NPNE,NVJE,NXI,
     &    NYNP,BBM,MASS,XAB,XP,YP(1,1),ERROR,*9999)
c        IF(.NOT.RESTART)THEN
          IDEAL_MASS(nx)=MASS !initial mass
c        ENDIF

! AJS: Causes mass errors for iterated solve steps (i.e. after the first solve)
!           IDEAL_MASS(nx)=CONC_INIT*INITIAL_VOLUME

      ENDIF

      UPDATE_MATRIX=.TRUE.
      NHST=NHP(NPNODE(1)) !find better place for this
      SMOOTHING=.FALSE.
C***  Set initial (or restarted) parameters for problem type
      IF(ITYP3(nr,nx).EQ.3.OR.ITYP3(nr,nx).EQ.6) THEN
C     Set initial parameters for pulmonary circulation
        CALL INIT_PULM_CIRC(ISC_GK,ISR_GK,ISC_GKK,ISR_GKK,nb,NBJ,NEELEM,
     &    NENP,NHP,NHST,NPNE,NPNODE,nr,NYNE,NYNP,NYNR,nx,NXI,nxl,
     &    SPARSE_P,SPARSE_S,T,DYNAM1,FIX,LINEAR,UPDATE_VECTOR,ERROR,
     &    *9999)
      ELSE
        CALL INIT_PULM(nb,NBJ,nbreath,NEELEM,NHP,NHST,NKJE,NORD,
     &    NPNE,NPNODE,nr,NVJE,NVJP,NYNP,nx,NXI,N_TM,SPARSE_P,SPARSE_S,
     &    ACINUS,BBM,CE,MASS,SUM_N2,SUM_WATER,T,water_loss,XAB,XP,YP,
     &    DYNAM1,FIRST_A,FIX,GAS,LINEAR,SCALE,SMOOTHING,UPDATE_NY,
     &    UPDATE_VECTOR,ERROR,*9999)
      ENDIF
      CONTINUE=.TRUE. !time-stepping continues while this is set
      DO WHILE(CONTINUE)
        T=T+DT !increment the time

        BreathIterations=BreathIterations+1
        N_SOLN=N_SOLN+1 !counts solution # for output at #th time
        IF(TFINISH+0.5d0*DT-T.GT.ZERO_TOL)THEN !continue
          IF(ITYP3(nr,nx).GE.3) THEN !capillaries
C... Defines time dependent pressure BC, currently only for CAP_INLET
C... currently only does linear interpolation between intervals
            IF(KTYP3_INIT(nx).EQ.4) THEN ! time variable through iptime
              DO ntv=1,NTIMEVARST
                INTERVAL=.FALSE.
                ntp=0
                DO WHILE(.NOT.INTERVAL) !determines which variable BC
                  ntp=ntp+1 !time interval we are in
                  IF(ntp.EQ.NTIME_POINTS(1)) THEN !final time interval
                    IF(T.EQ.TIME_VALUES(1,ntp,ntv)) INTERVAL=.TRUE.
                  ELSE
                    IF(T.GE.TIME_VALUES(1,ntp,ntv).AND.T.LT.
     '                TIME_VALUES(1,ntp+1,ntv)) INTERVAL=.TRUE.
                  ENDIF
                ENDDO !WHILE
C!!! NB/ CAP_INLET not defined anymore when multiple inlets
                IF(CAP_INLET.GT.0) THEN
                  ny=NYNP(1,1,1,NPNE(1,nb,CAP_INLET),0,1) !inlet ny #
                  gradient=(TIME_VALUES(2,ntp+1,ntv)-
     '              TIME_VALUES(2,ntp,ntv))/(TIME_VALUES(1,ntp+1,ntv)
     '              -TIME_VALUES(1,ntp,ntv))
                  YP(ny,1)=TIME_VALUES(2,ntp,ntv)+
     '              (T-TIME_VALUES(1,ntp,ntv))*gradient
                  YP(ny,3)=YP(ny,1)
                ENDIF !CAP_INLET
              ENDDO !ntv
            ENDIF !KTYP
            IF(CAP_INLET.GT.0) THEN
              ny=NYNP(1,1,1,NPNE(1,nb,CAP_INLET),0,1) !inlet ny #
              WRITE(OP_STRING,'('' Time='',F12.6,'' Pressure= '',
     &          F12.6)') t,YP(ny,1)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF !ITYP
C***      Call a pre-solution routine to calculate volume changes,
C***      flow distributions, and calculate boundary condition values.
          CALL PRE_SOL(nb,NBJ,NEELEM,NENP,NORD,NPNE,nr,
     &      NVJE,nx,NXI,NYNE,NYNP,BBM,CE,XAB,XP,YP(1,1),
     &      T,FIRST_A,GAS,SCALE,UPDATE_MATRIX,
     &      UPDATE_VECTOR,update_ny,ERROR,*9999)
C***      Assemble the system of matrices, reduce, and solve

          CALL CALCMASS(NBJ,NEELEM,NH_LOC(1,nx),NORD,NPNE,NVJE,NXI,
     &    NYNP,BBM,MASS,XAB,XP,YP(1,1),ERROR,*9999)
          IF(FIRST_NX) FIRST_A=.TRUE.
          CALL SOLVE11(IBT,IDO,INP,ISC_GD,ISR_GD,ISC_GK,ISC_GKK,ISR_GK,
     &      ISR_GKK,nb,NBH,NBJ,NEELEM,NELIST,NENP,NHE,NHST,NKJE,NO_NE,
     &      NONY,NORD,NPF,NPLIST2,NPNE,NPNODE,NPNY,nr,NVJE,NVJP,nx,NXI,
     &      NYNE,NYNO,NYNP,NYNR,NY_OFST,NZ_ESED,SPARSE_S,BBM,CE,CG,CGE,
     &      CONY,CYNO,CP,ED,EM,ER,ES,GD,GK,GKK,GK2,GR,GRR,PG,
     &      RG,SE,STACK_ED,STACK_EM,STACK_ES,WG,XA,
     &      XAB,XG,XO,XP,YG,YP,ZE,ZG,DYNAM1,FIRST_A,FIX,LINEAR,
     &      SMOOTHING,UPDATE_MATRIX,UPDATE_VECTOR,UPDATE_VELOCITY,
     &      FIRST_NX,ERROR,*9999)
          IF(ITYP3(nr,nx).LT.3)THEN
C AJS calculate gas exchange flux for each BBM acinus
            IF(ITYP3(nr,nx).EQ.1.AND.ITYP7(nr,nx).GT.1)THEN
              CALL GAS_EXCH_ACINUS(NBJ,NEELEM,NENP,NH_LOC(1,nx),NPNE,nr,
     &          NVJE,NXI,NYNP(1,1,1,1,0,1),ACINUS,BBM,CE,T,XAB,XP,YP,
     &          ERROR,*9999)
            ENDIF
            CALL CALCMASS(NBJ,NEELEM,NH_LOC(1,nx),NORD,NPNE,NVJE,NXI,
     &        NYNP,BBM,MASS_EXP,XAB,XP,YP(1,1),ERROR,*9999)
            CALL MESH_DEFORM(nb,NBJ,NEELEM,NELIST,NORD,NPNE,nr,NVJE,nx,
     &        NXI,NYNP,BBM,CE,XAB,XP,YP(1,1),FIX,GAS,ERROR,
     &        *9999)
            CALL CALCMASS(NBJ,NEELEM,NH_LOC(1,nx),NORD,NPNE,NVJE,NXI,
     &        NYNP,BBM,MASS,XAB,XP,YP(1,1),ERROR,*9999)
            MASS_EXP=MASS-MASS_EXP !the mass change due to acinus expansion
          ENDIF
          IF(ITYP3(nr,nx).EQ.4.AND.COUPLE_VIA_LPM.EQ.'Y') THEN
            !CALL PERFUSION_OUTPUTS to write out anything relating to perfusion models
            CALL PERFUSION_OUTPUTS(nb,NEELEM,NORD,NPNE,nx,NXI,NVJE,NYNE,
     &        NYNP,XAB,XP,YP,ERROR,*9999) 
          ELSEIF(ITYP3(nr,nx).EQ.2.OR.ITYP3(nr,nx).EQ.3.OR.
     &        ITYP3(nr,nx).EQ.5)THEN
            NELIST2_PTR=0
            CALL ALLOCATE_MEMORY(NE_R_M+1,1,INTTYPE,NELIST2_PTR,
     '        MEM_INIT,ERROR,*9999)
            CALL POST_SOL(nb,NBJ,NEELEM,NELIST,%VAL(NELIST2_PTR),NENP,
     &        NORD,NPNE,NPNODE,nr,NVJE,nx,NXI,NYNE,NYNP(1,1,1,1,0,1),
     '        ACINUS,BBM,CE,CP,T,XAB,XP,YP,GAS,UPDATE_MATRIX,
     &        ERROR,*9999)
            CALL FREE_MEMORY(NELIST2_PTR,ERROR,*9999)
          ENDIF
          FIRST_BREATH=.FALSE.
        ELSE
          CONTINUE=.FALSE.
        ENDIF !T
      ENDDO !CONTINUE
      T_LAST=T-DT
c       CLOSE(51)
      IF(OPENF_HISTORY)THEN !close output files (opened in INIT_PULM)
        IF(LUNG_OP.GE.1)CALL CLOSEF(IOFILE2,ERROR,*9999)
        IF(LUNG_OP.GE.2)CALL CLOSEF(IOFILE3,ERROR,*9999)
        IF(LUNG_OP.GE.3)CALL CLOSEF(IOFILE4,ERROR,*9999)
      ENDIF
      
C AJS redundant code
! 998   IF(ITYP3(nr,nx).EQ.7)CALL GAS_EXCHANGE(NPNODE,nr,nx,XP,
!      &   ERROR,*9999)!Gas exchange

      CALL EXITS('MARCH11')
      RETURN

 9999 CALL ERRORS('MARCH11',ERROR)
      CALL EXITS('MARCH11')
      RETURN 1
      END



