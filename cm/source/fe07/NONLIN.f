      SUBROUTINE NONLIN(IBT,IDO,INP,ISEG,ISELNO,ISFIBR,ISFIEL,
     '  ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,ISR_GKK,ISR_GQ,ISLINE,ISLINO,
     '  ISNONO,ISTATE,IUSER,IWK,LGE,MXI,NAN,NBH,NBHF,NBJ,NBJF,NEELEM,
     '  NENP,NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,NKHE,NKJE,
     '  NLATNE,NLL,NLLIST,NMNO,NNB,NNF,NONY,NPF,NP_INTERFACE,NPL,NPNE,
     '  NPNODE,NPNY,NQNE,NQNLAT,NQS,NQXI,nr_solve,NRE,NRLIST,NSB,
     '  NTIW,NVHE,NVHP,NVJE,NW,nx,NXI,NYNE,NYNO,NYNP,NYNR,Z_CONT_LIST,
     '  CE,CG,CGE,CONY,CP,CURVCORRECT,CYNO,DL,
     '  ERRMAX,FEXT,GK,GKK,GQ,GR,GRR,PAOPTI,PG,PMIN,PMAX,
     '  R,RE1,RESID,RESJAC,RG,SE,USER,WARM_START,WG,XA,
     '  XC,XE,XG,XO,XP,XQ,YG,YGF,YP,YQ,ZA,ZA1,Z_CONT,ZE,ZE1,ZP,
     '  ZP1,CSEG,FIX,ERROR,*)

C#### Subroutine: NONLIN
C###  Description:
C###    NONLIN controls the solution of non-linear,quasi-static problems
C###    using Newton-Raphson, Modified-Newton or Sequential Quadratic.

C**** Programming iterative schemes with line search.
C****  YP(ny,1) and ZP,ZA has current equilibrium solution
C****  YP(ny,2) has prescribed dep var/force increms set by FIX(ny,1)
C****  YP(ny,3) has prescribed initial equilibrium solution
C****  YP(ny,4) has current set of equilibrium equation residuals
C****  YP(ny,5) is temporary storage 1 used to store solution
C****           increments to add to YP(ny,1)
C****  YP(ny,10) is temporary storage 2 used to store current reference
C****           solution for FE50 cavity elements (set up in IPINI5/
C****           UPSOLU)
C 28/3/96 These iy's need to be checked below.
C****  YP(ny,8..16) is used for update information in Broyden
C****    schemes, or for previous time solutions in Newton &
C****    modified Newton schemes.
C**** ITER1 counts the number of equilibrium iterations.
C**** ITER2    "         "       iterations since last update.
C**** NWRIT    "         "       load steps since last output.
C**** NTLOAD is number of load steps
C**** NTITER is maximum number of iterations
C**** REITER is set .TRUE. if NTLOAD=0,indicating restart of iterations
C**** If KTYP71=1 pressure loads are read from FILE01.PRESSURE
C OLD C**** If KTYP57(nr)=3 pressure loads are read from FILE07.PRESSURE
C**** KTYP1A is 1,2 for series/parallel element stiffness matrix calcs
C**** KTYP1D determines the method for derivative calculation
C**** Form the constraint reduced matrix (non-zeros only) GKK from GK.
C**** The mapping coeffs NONY(noy,ny,nrc,nr),CONY(noy,ny,nrc,nr),
C**** noy=1,NONY(0,ny,nrc,nr)
C**** are calc.d to reduce system of equns NYNR(no_nynr,nrc,nc,nr),
C**** no_nynr=1,NYNR(0,nrc,nc,nr)
C**** to the system no=1,NOT(nrc,nc,nr,nx) by removing constraints.
C**** NONY(0,ny,nrc,nr) is 0,1 when FIX(ny,1) is .true.,.false., resp.

C**** JHC 24/06/05 Restructured to have while loop rather than GOTO's
C**** and to eliminate repetitions of codes 

C**** JHC 26/06/05 for contact problems, arrays & counters used are:
C****  YP(ny,1)  current equilibrium solution - T+dt
C****  YP(ny,2)  prescribed dep var/force increms set by FIX(ny,1)
C****  YP(ny,3)  prescribed initial equilibrium solution
C****  YP(ny,4)  current set of equilibrium equation residuals
C****  YP(ny,5)  temporary storage 1 used to store solution
C****            increments to add to YP(ny,1)
C****  YP(ny,6)  storage of contact residual vector
C****  YP(ny,8)  previous solution at iteration i-1 
C****  YP(ny,11) velocity at iteration i

C***   CONT_IT - Total number of Newton Iteration counter (all loads, all augumentations)
C***   LOAD_IT - Load step counter (including Newton step iterations for augmentations for a particular load)
C***   AUG_IT - Augmentation counter
C***   CONV_IT - Convergence check counter (for a particular load, for an augumentation)

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'acti01.cmn'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'disp00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'fsklib.inc'
      INCLUDE 'gen000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'gks000.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'host00.cmn'
      INCLUDE 'host00.inc'
      INCLUDE 'ipma50.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'load00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'ktyp70.cmn'
      INCLUDE 'nonl00.cmn'
      INCLUDE 'ofst00.cmn'
      INCLUDE 'press00.cmn'
      INCLUDE 'ptr00.cmn'
      INCLUDE 'solv00.cmn'
      INCLUDE 'time02.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),ISEG(*),ISELNO(NWM,NEM),
     '  ISFIBR(NWM,NEM,NGRSEGM),ISFIEL(NWM,NEM),ISC_GK(NISC_GKM),
     '  ISC_GKK(NISC_GKKM),ISC_GQ(NISC_GQM),ISR_GK(NISR_GKM),
     '  ISR_GKK(NISR_GKKM),ISR_GQ(NISR_GQM),ISLINE(NWM,2*NGRSEGM),
     '  ISLINO(NWM),ISNONO(NWM,NPM),ISTATE(*),IUSER(*),IWK(*),
     '  LGE(NHM*NSM,NRCM),MXI(2,NEM),NAN(NIM,NAM,NBFM),
     '  NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NEELEM(0:NE_R_M,0:NRM),NENP(NPM,0:NEPM,0:NRM),
     '  NFF(6,NEM),NFFACE(0:NF_R_M,NRM),NGAP(NIM,NBM),
     '  NHE(NEM),NHP(NPM,0:NRM),NKB(2,2,2,NNM,NBFM),
     '  NKEF(0:4,16,6,NBFM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NLATNE(NEQM+1),
     '  NLL(12,NEM),NLLIST(0:NLM),
     '  NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),NONY(0:NOYM,NYM,NRCM,0:NRM),
     '  NMNO(1:2,0:NOPM),NPF(9,NFM),NP_INTERFACE(0:NPM,0:3),
     '  NPL(5,0:3,NLM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM),
     '  NQNE(NEQM,NQEM),NQNLAT(NEQM*NQEM),NQS(NEQM),NQXI(0:NIM,NQSCM),
     '  nr_solve,NRE(NEM),NRLIST(0:NRM),NSB(NKM,NNM,NBFM),NTIW,
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx,
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNO(0:NYOM,NOOPM,NRCM,0:NRM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM),Z_CONT_LIST(NDM,2,7)
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),CGE(NMM,NGM,NEM),
     '  CONY(0:NOYM,NYM,NRCM,0:NRM),CP(NMM,NPM),
     '  CURVCORRECT(2,2,NNM,NEM),CYNO(0:NYOM,NOOPM,NRCM,0:NRM),
     '  DL(3,NLM),
     '  ERRMAX,FEXT(NIFEXTM,NGM,NEM),GK(NZ_GK_M),
     '  GKK(NZ_GKK_M),GQ(NZ_GQ_M),GR(NYROWM),GRR(NOM),
     '  PAOPTI(*),PG(NSM,NUM,NGM,NBM),PMIN(*),PMAX(*),R(NOPM,*),
     '  RE1(NSM,NHM),RESID(*),RESJAC(NREM,*),RG(NGM),
     '  SE(NSM,NBFM,NEM),USER(*),WG(NGM,NBM),
     '  XA(NAM,NJM,NEM),
     '  XC(*),XE(NSM,NJM),XG(NJM,NUM),XO(NOM),XP(NKM,NVM,NJM,NPM),
     '  XQ(NJM,NQM),YG(NIYGM,NGM,NEM),YGF(NIYGFM,NGFM,NFM),
     '  YP(NYM,NIYM),YQ(NYQM,NIQM,NAM,NXM),
     '  ZA(NAM,NHM,NCM,NEM),ZA1(NAM,NHM,NCM,NEM),
     &  Z_CONT(NDM,2,67),ZE(NSM,NHM),ZE1(NSM,NHM),
     '  ZP(NKM,NVM,NHM,NPM,NCM),ZP1(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM),WARM_START
!     Local Variables
      INTEGER ERR,i,IBEG,
     '  IE,IEND,il,IOSTAT,IP,
     '  IREC,ITER,IEVAL,ITER1,ITER2,j,METHOD,
     '  nc,nd,ne,ng,nh,nhx,nj,njj,nk,no,noelem,noload,
     '  nonrlist,no_nynr,noy,np,nr,nro,nr1,nv,NWRIT,ny,nyo,ny1,ny2
      REAL TIMEA_CPU,TIMEA_ELAPSED,TIMEB_CPU,
     '  ELAPSED_TIME,TIMEB_ELAPSED,TIME_START1(1),TIME_START1_REAL(1),
     '  TIME_STOP(1),TIMED_CPU,TIMED_ELAPSED
      REAL*8 ALPHA,CONT_PTS_CONT,CONT_PTS_TIED,GAP_CHECK_CONT,
     '  GAP_CHECK_TIED,GAP_SUM_CONT,GAP_SUM_TIED,GRR_TEMP(NOM),
     '  OLDRAT,OLDRS3,OPTION(56),
     '  RATIO(0:6),REF_ENERGY,RSUM_SOLINCR,
     '  ZERO_TOL,TOTAL_PRESS,LAG_MUL_CONV,LAG_MUL_RATIO,CURRENT_LAG_MUL

C MHT 26-11-10. Modification to inappropriately declared large arrays ZAA,ZPA
C These cause segmentation faults for mechanics, and do not appear to be
C used in NONLIN other than being passed to YPZP. 
C If my reduction of these arrays causes a problem, then whoever put 
C these in needs to include code to declare them correctly.

C !!! Whoop! Whoop! Large arrays only used for some problems.
c      REAL*8 ZAA(NAM,NHM,NCM,NEM),ZPA(NKM,NVM,NHM,NPM,NCM)
       REAL*8 ZAA(1,1,1,1),ZPA(1,1,1,1,1)

C ***  ZPA,ZAA store acceleration information 

!     REAL*8 LENGTH_SOL,LENGTH
      CHARACTER CHAR1*100
      LOGICAL ADD_GRAVITY,CONVERGED,CONTACT,ERROR_FLAG,FREE_VAR,OPENED,
     &  OUTPUT,REITER,SOLVE_SYSTEM,UPDATE_MATRIX,UPNODE
!     Functions
!     MINLSSQP Variables
      INTEGER LDFJ,LDR
      INTEGER IFAIL
      REAL*8 OBJF
      EXTERNAL FUNCT4

      DATA ZERO_TOL /1.0d-8/

      CALL ENTERS('NONLIN',*9999)
      IF(KTYP10.EQ.2) THEN
        CALL ASSERT(NIYM.GE.8,'>>Increase NIYM to 8 for line searches',
     '    ERROR,*9999)
      ELSE IF(KTYP26.EQ.1.AND.KTYP27.EQ.2.AND.KTYP28.EQ.0) THEN
        CALL ASSERT(NIYM.GE.8,'>>Increase NIYM to 6',
     '    ERROR,*9999)
      ELSE
        CALL ASSERT(NIYM.GE.5,
     '    '>>Increase NIYM to 5 for general NONLIN problems',
     '    ERROR,*9999)
      ENDIF
      IF(UPVU) THEN
        CALL ASSERT(NIYFIXM.GE.2,'>>Increase NIYFIXM to 2 for graphics',
     '    ERROR,*9999)
      ELSE
        CALL ASSERT(NIYFIXM.GE.1,
     '    '>>Increase NIYFIXM to 1 for general NONLIN problems',
     '    ERROR,*9999)
      ENDIF

C NEWS 30/06/05 JHC check for coupled or not coupled Contact Mechanics problems
      CONTACT=.FALSE.
      IF(nr_solve.EQ.0) THEN ! coupled problem
        CONTACT=.TRUE.
        DO nonrlist=1,NRLIST(0)
          nr=NRLIST(nonrlist)  
          IF(KTYP5G(nr).GT.0.AND.
     &      CONTACT) THEN!contact
            CONTACT=.TRUE. ! coupled contact
          ELSE
            CONTACT=.FALSE. ! it is coupled problem but not contact
          ENDIF
        ENDDO
      ELSE IF(KTYP5G(nr_solve).GT.0) THEN ! not coupled contact
        CONTACT=.TRUE.
      ENDIF
C NEWE JHC   

C *** NJM must be 67 or greater to store contact point field values.
C *** in WD(NJM,nd) and ZD(NJM,nd) updated in UPDATA (fe23)
      IF(CONTACT) THEN ! contact
        CALL ASSERT(NJM.GE.67,
     &    '>>Increase NJM to 67 for general CONTACT problems',
     &    ERROR,*9999)
C*** 22/04/08 JHC Check if .ipcont has been defined for contact problems
        CALL ASSERT(CALL_CONT,'>>Contact parameters not defined',
     &    ERROR,*9999)
      ENDIF
      
      DO nonrlist=1,NRLIST(0)
        nr=NRLIST(nonrlist)
        IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN
C         constant volume constraint (FE50 cavity problems)
          CALL ASSERT(NIYM.GE.10,'>>Increase NIYM in parameters.inc '
     '      //'(>=10 for Const vol constr problems)',ERROR,*9999)
        ENDIF
C KAT 2005-12-14: NTACTV never set
C         IF(ITYP5(nr,nx).EQ.3.AND.NIFEXTM.LT.4+NTACTV) THEN
C         WRITE(CHAR,'(I1)') IDIGITS(4+NTACTV)
C         WRITE(ERROR,'(''>>Increase NIFEXTM to '',I'//CHAR//')')
C      '    4+NTACTV
        IF(ITYP5(nr,nx).EQ.3) THEN
          ERROR='Modal analysis is not implemented'
          GO TO 9999
        ENDIF
        IF(ITYP1(nr,nx).EQ.5) THEN
          CALL ASSERT(USE_NONLIN.GT.0,'>>Cannot solve NONLIN problems '
     '      //'as USE_NONLIN is zero in parameters.inc',ERROR,*9999)
          CALL ASSERT(NIFEXTM.GE.7,
     '      '>>Increase NIFEXTM to 7 for finite deformation problems',
     '      ERROR,*9999)
        ENDIF
      ENDDO !nr

      LDFJ=NREM
      LDR=NOPM
      TIME=TSTART
      DTIME=DT
      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP   CRITICAL(NONLIN_1)
        WRITE(OP_STRING,'('' Time='',D11.4,'' DTime='',D11.4,'
     '    //''' NTLOAD='',I2)')TIME,DTIME,NTLOAD
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP   END CRITICAL(NONLIN_1)
      ENDIF

      nr1=NRLIST(1)
      METHOD=ITYP9(nr1,nx)
C *** Check that solution methods are the same in all regions
      DO nonrlist=2,NRLIST(0)
        nr=NRLIST(nonrlist)
        CALL ASSERT(ITYP9(nr,nx).EQ.METHOD,
     '    '>>Iteration method not same in all regions',ERROR,*9999)
        CALL ASSERT(ITYP5(nr,nx).EQ.ITYP5(nr1,nx),
     '    '>>Problem type not same in all regions',ERROR,*9999)
      ENDDO

C *** Loop over all regions in list
      FREE_VAR=.FALSE.
      DO nonrlist=1,NRLIST(0)
        nr=NRLIST(nonrlist)
        IF(KTYP5.EQ.1.AND.ITYP1(nr,nx).EQ.5) THEN
C new MPN 1Feb2000: OMP parallel proc directive
C old C$DOACROSS local (ne,noelem,ng) share (nr,nx)
C$OMP     PARALLEL DO
C$OMP&      PRIVATE(ne,noelem,ng),
C$OMP&      SHARED(nr,nx)
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            DO ng=1,NGT(NBH(NH_LOC(1,nx),1,ne))
              FEXT(2,ng,ne)=1.0D0
              FEXT(5,ng,ne)=0.0D0
              FEXT(6,ng,ne)=0.0D0
              FEXT(7,ng,ne)=0.0D0
            ENDDO !ng
          ENDDO !noelem (ne)
C$OMP     END PARALLEL DO
        ENDIF
C ***   Solve equations (if at least one dof is not fixed)
        DO no_nynr=1,NYNR(0,0,1,nr) !loop over global variables
          ny=NYNR(no_nynr,0,1,nr) !is global variable number
          IF(.NOT.FIX(ny,1)) FREE_VAR=.TRUE.
        ENDDO !no_nynr
      ENDDO !nonrlist (nr)

C NEW 24/06/05 JHC Parallel stuff is moved to *** in Archive directory

      IF(NTLOAD.EQ.0) THEN
        REITER=.TRUE.
        NTLOAD=1
      ELSEIF(ITYP5(nr1,nx).EQ.5) THEN
C       Wavefront continuation method.
C       FACTOR is used in ZERE30 as the continuation parameter
        CALL ASSERT(NTLOAD.EQ.1,
     '    '>>Can only perform one increment at a time',ERROR,*9999)
        REITER=.TRUE.
      ELSE
        REITER=.FALSE.
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_START1)
      CALL REAL_TIMER(REAL_TOTAL,TIME_START1_REAL)
      TIMEA_CPU=0.0
      TIMEA_ELAPSED=0.0

      OUTPUT=.FALSE.
      NWRIT=1
      DO noload=1,NTLOAD
        NOSTEP=noload !to pass via /NONL/ (nonl00.cmn) to other routines
        RATIO(0)=0.0d0
C KAT 5Nov98: I don't think this is what FIRSTS is for.
C        FIRSTS=.TRUE.
        IF(IWRIT1(nr_solve,nx).GT.0) THEN
          IF(noload.EQ.1.OR.IWRIT2(nr_solve,nx).EQ.2.
     '      OR.NWRIT.EQ.IWRIT1(nr_solve,nx)) THEN
            OUTPUT=.TRUE. !output initial residual
            NWRIT=1
          ELSE
            OUTPUT=.FALSE.
            NWRIT=NWRIT+1
          ENDIF
        ENDIF

C*** 19/03/08 JHC Set the flag ADD_GRAVITY to be true
C                 This is to allow gravity vector to be set up 
C                 only once for multi-region problems such as contact mechanics 
C                 And there is no reason to have different gravity vector for
C                 different region
        ADD_GRAVITY=.TRUE.

        DO nonrlist=1,NRLIST(0)
          nr=NRLIST(nonrlist)

          IF(.NOT.REITER) THEN
            IF(nonrlist.EQ.1) THEN
              WRITE(OP_STRING,'(/'' Load step'',I3/,1X,12(''=''))')
     '          noload
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF

C ***       Apply external pressure loads
            IF(ITYP1(nr,nx).EQ.4) THEN !FE40 problem
              IF(KTYP71.EQ.1) THEN
                CALL STRING_TRIM(FILE01,IBEG,IEND)
 201            CALL OPENF(IOFILE2,'DISK',
     '            FILE01(IBEG:IEND)//'.pressure',
     '            'OLD','DIRECT','FORMATTED',132,ERROR,*203)
 203            IF(ERROR(1:10).EQ.'Iostat= 30') THEN
                  WRITE(OP_STRING,'('' Pressure file is locked'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  GO TO 201
                ELSE IF(ERROR(2:2).ne.' ') THEN
                  GO TO 9999
                ENDIF
                IF(NW(ne,1).EQ.5) THEN
                  DO noelem=1,NEELEM(0,nr)
                    ne=NEELEM(noelem,nr)
                    INQUIRE(UNIT=IOFILE2,OPENED=OPENED,NEXTREC=IREC)
                    IF(OPENED) THEN
                      READ(IOFILE2,FMT='(9D12.4)',REC=IREC,
     '                  IOSTAT=IOSTAT) (PRESS(ng,ne),ng=1,9)
                      IF(IWRIT4(nr,nx).GE.5) THEN
                        WRITE(OP_STRING,
     '                    '('' Pressures in element '',I3,'' are '','
     '                    //'D12.4)') ne,(PRESS(ng,ne),ng=1,9)
                        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      ENDIF
                      IF(iostat.ne.0) THEN
                        IF(IOSTAT.EQ.36) THEN
                          ERROR=' End of pressure file'
                        ELSE
                          ERROR=' Real data read error in pressure file'
                        ENDIF
                        GO TO 9999
                      ENDIF
                    ELSE
                      ERROR=' Unit 14 is not open'
                      GO TO 9999
                    ENDIF
                  ENDDO !noelem (ne)
                ENDIF
                CALL CLOSEF(IOFILE2,ERROR,*9999)
              ENDIF
            ENDIF !ityp1 FE40problem type
C NEWS 26/06/05 JHC Add contact information
C ***       Setup previous solution and velocity information at t=0.0
            IF (KTYP5G(nr).GT.0) THEN ! contact
              IF(CONT_IT.EQ.1) THEN
                DO nc=1,NCT(nr,nx) !loop over RHS(displ) and LHS(force) vars
                  DO no_nynr=1,NYNR(0,0,nc,nr) !loop over global variables
                    ny=NYNR(no_nynr,0,nc,nr) !is global variable number
                    IF(NPNY(0,ny,0).EQ.1) THEN
                      np=NPNY(4,ny,0)
                      CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                    ENDIF
                    IF (KTYP5I(nr).EQ.1) THEN ! inertia
                      YP(ny,12)=0.0d0 ! Previous velocity=0.0 at time=0.0
                      YP(ny,8)=YP(ny,1) ! Previous solution at time=0.0
                    ENDIF 
                  ENDDO !no_nynr (ny)
                ENDDO !nc
              ENDIF ! CONT_IT

              IF(LOAD_IT.EQ.1) THEN
C ***           Check which points are in contact
C ***           Only do once at beginning of each load step
                DO nd=1,NDT
C ***             Check signed normal gap Z_CONT(nd,j,13) for slave (j=1)
                  IF (Z_CONT_LIST(nd,1,4).NE.2) THEN ! contact
                    IF (Z_CONT(nd,1,13).GE.GAP_TOL) THEN ! then in contact
                      Z_CONT_LIST(nd,1,3)=1
                      Z_CONT_LIST(nd,2,3)=1
                    ELSE ! not in contact
C ***                 Check if problem starts with active contraints
C*** 22/04/08 JHC If user defines CHECK_CONTACT to be true in .ipcont file
C                 then make all contact points to be active
                      IF (KTYP5J(nr).GE.1.OR.CHECK_CONTACT) THEN ! active
                        Z_CONT_LIST(nd,1,3)=1
                        Z_CONT_LIST(nd,2,3)=1
                      ELSE ! not active
                        Z_CONT_LIST(nd,1,3)=0
                        Z_CONT_LIST(nd,2,3)=0
                      ENDIF ! KTYP5J
                    ENDIF
                  ELSE ! tied contact -always in contact
                    Z_CONT_LIST(nd,1,3)=1
                    Z_CONT_LIST(nd,2,3)=1
                  ENDIF                  
                ENDDO !nd 
                KTYP5J(nr)=0 ! turn off initial active constraint
              ENDIF !LOAD_IT
            ENDIF ! KTYP5G
C NEWE JHC
C ***       Apply increments to displacement and force b.c.s
C NEWS 26/06/05 JHC for contact problem, only increment at the first load step
            IF((KTYP5G(nr).EQ.0).OR. ! non-contact
     &        ((KTYP5G(nr).GT.0).AND.(LOAD_IT.EQ.1))) THEN ! contact and 1st load step
              i=0
              DO nc=1,NCT(nr,nx) !loop over RHS(displ) and LHS(force) vars
                DO no_nynr=1,NYNR(0,0,nc,nr) !loop over global variables
                  ny=NYNR(no_nynr,0,nc,nr) !is global variable number
                  IF(NPNY(0,ny,0).EQ.1) THEN
                    np=NPNY(4,ny,0)
C GMH 8/1/97 Update cmgui link
                    CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                  ENDIF
C??? KAT 27Aug97:  Shouldn't this be FIX(ny,2)?
                  IF(FIX(ny,1)) THEN !ny has incremented ess bdry cond
C KAT 2005-12-14: TBDRY never set
C                     IF(ITYP5(nr,nx).EQ.3) THEN !modal analysis
C                        YP(ny,1)=TBDRY(noload)
C                     ELSE
                      IF(KTYP26.EQ.1.AND.KTYP27.EQ.2.
     &                AND.KTYP28.EQ.0) THEN

c!!! cpb 15/5/95 this is wrong iy=6 is not undeformed geometry

                      CALL ASSERT(.FALSE.,'>>Check code wrt new iy '
     '                  //'scheme',ERROR,*9999)

C                     Update nodal params from optimisation params
C                     (YP(ny,6) contains undeformed geometry) ??? CHECK
                      i=i+1
                      YP(ny,1)=YP(ny,6)+PAOPTI(NMNO(1,0)+i)
                    ELSE !increment boundary conditions
C                     YP(ny,2) contains dep. var. increments
C news VJ 18Nov2004: Adding option to solve for undeformed state given the deformed state
                      IF(NPNY(0,ny,0).EQ.1) THEN !node based ny
                        nh=NPNY(3,ny,0)
                      ELSE IF(NPNY(0,ny,0).EQ.2) THEN !element based ny
                        nh=NPNY(2,ny,0)
                      ENDIF
                      IF(KTYP5L.EQ.1.OR. !solving for deformed coordinates
     &                  nh.GT.NJ_LOC(NJL_GEOM,0,nr).OR. !nh not a geometric variable
     &                  NPNY(0,ny,0).EQ.2.OR. !ny is element based
     &                  nc.EQ.2) THEN !dealing with force components
                        YP(ny,1)=YP(ny,1)+YP(ny,2)*FACTOR
                      ELSE IF(KTYP5L.EQ.2) THEN !Solving for undeformed coordinates
C ***  RY/MPN  14/11/06     Applying BC to cavity reference state (YP(ny,10)) when solving for undeformed
                        UPNODE=.TRUE.
                        nk=NPNY(1,ny,0)
                        nv=NPNY(2,ny,0)
                        njj=nh !sensible default
C VJ check previous regions to see if np has another ny associated with it due to coupling
                        DO nro=1,nr
                          ny1=NYNP(nk,nv,nh,np,0,1,nro) !nrc=0 selects global variable
                          IF(ny1.LT.ny.AND.ny1.GT.0) THEN
                            UPNODE=.FALSE.
                          ENDIF
                        ENDDO
                        IF(UPNODE) THEN
                          DO nhx=1,NH_LOC(0,nx)
                            IF(nh.EQ.NH_LOC(nhx,nx)) THEN
                              njj=nhx
                            ENDIF
                          ENDDO
                          nj=NJ_LOC(NJL_GEOM,njj,nr)
                          XP(nk,nv,nj,np)=YP(ny,1)-YP(ny,2)*FACTOR
                        ENDIF
                        IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) 
     &                  THEN !const vol
                          YP(ny,10)=YP(ny,10)-YP(ny,2)*FACTOR
                        ENDIF
                      ENDIF  
                    ENDIF
                  ENDIF !FIX(ny,1)
C NEWS 26/06/05 JHC for contact,
                  IF(KTYP5G(nr).GT.0) THEN  
                    IF (KTYP5I(nr).EQ.1) THEN ! inertia
                      YP(ny,11)=(YP(ny,1)-YP(ny,8))/T_inc ! velocity 
                      YP(ny,14)=(YP(ny,11)-YP(ny,12))/T_inc ! acceleration 
                    ENDIF
                    YP(ny,8)=YP(ny,1) ! Update previous solution (for augmentation) 
                  ENDIF ! KTYP5G
C NEWE JHC
                ENDDO !no_nynr (ny)
              ENDDO !nc
C VJ 6Mar2005: Adding ability to increment body force vector
C b_inc is the incremental body force vector updated at each fem solve command
C b contains the force values entered in the ipinit file
C*** 19/03/08 JHC Prevent accumulation of gravity vector over multiple region loop 
              IF(ADD_GRAVITY) THEN
                b_inc(1)=b_inc(1)+b(1)*FACTOR
                b_inc(2)=b_inc(2)+b(2)*FACTOR
                b_inc(3)=b_inc(3)+b(3)*FACTOR
                ADD_GRAVITY=.FALSE. ! only increment gravity vector once
c                WRITE(OP_STRING,'(''Gravity: '',3(F8.2))') b_inc(1),
c     &            b_inc(2),b_inc(3)
c                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
C ***         Apply increment to material parameters
              IF(KTYP14.GT.0) THEN
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  IE=NW(ne,1)
                  DO il=1,ILT(IE,nr,nx)
                    CE(il,ne)=CE(il,ne)+CE(il+ILT(IE,nr,nx),ne)*FACTOR
                  ENDDO !il
                ENDDO !noelem (ne)
              ENDIF ! KTYP14
            ENDIF ! contact

            IF(nonrlist.EQ.1) TIME=TIME+DTIME
C NEWS 26/06/05 JHC Add contact information
            IF(KTYP5G(nr).GT.0) THEN ! for contact mechanics
              IF (CONT_IT.EQ.1) THEN
                IF (KTYP5I(nr).EQ.1) THEN ! inertia
C ***             Get ZPA,ZAA terms from YP(ny,11) : velocity 
                  CALL YPZP(11,NBH,NEELEM,NHE,NHP(1,nr),NKH(1,1,1,nr),
     '              NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZAA,ZPA,
     '              ERROR,*9999)
                ENDIF !KTYP5I
              ENDIF !CONT_IT
            ENDIF !KTYP5G
C NEWE JHC  
          ENDIF ! NOT.REITER
  
        ENDDO !nonrlist (nr)

C NEWS 23/06/05 JHC NEWTON STEP BEGINS HERE
        CONVERGED=.FALSE.

        IF(.NOT.FREE_VAR) THEN
          WRITE(OP_STRING,'(/'' no free variables'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          SOLVE_SYSTEM=.FALSE.
        ELSE
          IF(NTITER.EQ.0) THEN
            SOLVE_SYSTEM=.FALSE.
          ELSE
            SOLVE_SYSTEM=.TRUE.
          ENDIF
        ENDIF

        ITER1=0
        ITER2=0

C NEWS 01/07/05 JHC initialise OLDRAT and OLDRS so that before performing
C 1st Newton step, always assemble the matrix
        IF(ITER1.EQ.0) THEN
          OLDRAT=-1.0d0 ! ratio(0) will always be greater than OLDRAT
          OLDRS3=-1.0d0 ! sum of soln increments will always be greater than OLDRS3
        ENDIF
C NEWE JHC

        IF(WARM_START.AND.METHOD.EQ.2) THEN
          UPDATE_MATRIX=.FALSE.
        ELSE
          UPDATE_MATRIX=.TRUE.
        ENDIF
  
        DO WHILE (.NOT.CONVERGED.AND.
     &    SOLVE_SYSTEM)
C ***     Get residual
          DO nonrlist=1,NRLIST(0)
            nr=NRLIST(nonrlist)
            CALL YPZP(1,NBH,NEELEM,NHE,NHP(1,nr),NKH(1,1,1,nr),NPNODE,
     '        nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
C NEWS 28/06/05 JHC adding contact information           
            IF (KTYP5G(nr).GT.0.AND. ! contact
     &        KTYP5I(nr).EQ.1) THEN ! inertia
C ***         Get ZPA,ZAA terms from YP(ny,11) : velocity 
              CALL YPZP(11,NBH,NEELEM,NHE,NHP(1,nr),NKH(1,1,1,nr),
     '          NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZAA,ZPA,
     '          ERROR,*9999)
            ENDIF
C NEWE JHC
            IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN !const Vol
C             Put reference state for cavity from YP(ny,10) into
C             ZA1,ZP1 for ZERE55
              CALL YPZP(10,NBH,NEELEM,NHE,NHP(1,nr),NKH(1,1,1,nr),
     '          NPNODE,nr,
     '          NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZA1,ZP1,ERROR,*9999)
            ENDIF

            CALL ZPRP(IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,NEELEM,
     '        NFF,NFFACE,NGAP,NHE,NHP(1,nr),NKB,NKEF,NKH(1,1,1,nr),
     '        NKHE,NKJE,NNB,NNF,NPF,NPNE,NPNODE,NPNY,nr,NRE,NSB,
     '        NVHE,NVHP(1,1,1,nr),NVJE,NW,nx,NXI,NYNE,NYNP,
     '        NYNR(0,0,1,nr),Z_CONT_LIST,
     '        CE,CG,CGE,
     '        CP,CURVCORRECT,FEXT,PG,RE1,RG,SE,
     '        WG,XA,XG,XP,YG,YGF,YP,ZA,ZA1,%VAL(0),Z_CONT,ZE,ZE1,ZP,
     '        ZP1,%VAL(0),FIX,
     '        ERROR,*9999)

C 25/06/07 XSL Code for modifying residual vector with contact contributions
C is moved inside ZPRP

C           Output residual stats
            IF(OUTPUT.AND.IWRIT3(nr,nx).GT.0) THEN
              IF(REITER) THEN
                IP=1
              ELSE
                IP=2
              ENDIF
              CALL ZPOP(IP,NBH,1,NEELEM,NHP(1,nr),NKH(1,1,1,nr),NPNODE,
     '          nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,4),ZA,ZP,FIX(1,1),
     '          ERROR,*9999)
            ENDIF
          ENDDO ! nonrlist

C 28/07/08 XSL NEWS Normalise initial energy for the 1st Newton step
          IF((.NOT.CONTACT.AND.ITER1.EQ.0).OR. ! 1st iteration, not contact
     &      (CONTACT.AND.CONV_IT.EQ.1)) THEN ! 1st iteration, contact
            INITIAL_ENERGY=0.0D0
          ENDIF
C NEWE

C*** 01/08/07 XSL Calculate and output convergence ratio
          CALL CALC_CONV_RATIO(ITER1,IOOP,NBH,ne,NEELEM,
     &      NONY,NPNY,nr,nr_solve,nx,NYNO,NYNR,
     &      CONY,ERRMAX,GRR,RATIO,YG,YP,
     &      CONTACT,OUTPUT,
     &      ERROR,*9999)

C 27/11/08 XSL NEWS Calculate reference energy for energy norm using line search
C The valule for REF_ENERGY is exactly the same as the output from CALC_COV_RATIO above
C Only do this calculation if energy_norm is used
          REF_ENERGY=0.0D0
          IF((KTYP007.EQ.4).AND.(KTYP10.EQ.2)) THEN ! energy norm with line search
            IF((.NOT.CONTACT.AND.ITER1.GT.0).OR. ! >1st iteration, not contact
     &        (CONTACT.AND.CONV_IT.GT.1)) THEN ! >1st iteration, contact
              DO no=1,NOT(1,1,nr_solve,nx) !loop over global soln rows
                DO nyo=1,NYNO(0,no,1,nr_solve)
                  ny1=NYNO(nyo,no,1,nr_solve) !is row number
                  ny2=NYNO(nyo,no,2,nr_solve) !is global variable number
                  REF_ENERGY=REF_ENERGY+YP(ny2,5)*YP(ny1,4)
                ENDDO !nyo
              ENDDO !no
              REF_ENERGY=DABS(REF_ENERGY)/INITIAL_ENERGY
            ENDIF !not contact and contact
          ENDIF !ktyp
C 27/11/08 NEWE

C         Calc sum of absolute solution vector increments
C NEWS 04/07/05 JHC Before performing the 1st Newton step, YP(ny,5) is empty as nothing has been solved
C Hence when iter1=0, i.e. before performing 1st Newton step, do not loop over YP(ny,5) but just set it to zero
          RSUM_SOLINCR=0.0D0
          IF((.NOT.CONTACT.AND.ITER1.GT.0).OR. ! after 1st Newton step, non-contact
     &      (CONTACT.AND.CONV_IT.GT.1)) THEN  ! after 1st Newton step, contact
            DO no=1,NOT(2,1,nr_solve,nx) !loop over global soln variables
              IF(NYNO(0,no,2,nr_solve).GT.0) THEN
                ny=NYNO(1,no,2,nr_solve)   !is first global variable number
                RSUM_SOLINCR=RSUM_SOLINCR+DABS(YP(ny,5)) !coupled to no
              ENDIF
            ENDDO !no
          ENDIF
C NEWE JHC

C NEWS 01/08/07 XSL Output sum of soln inc
C          IF(OUTPUT.OR.IWRIT4(nr_solve,nx).GT.1) THEN !output residuals
C            IF((IWRIT4(nr_solve,nx).GT.1).AND.KTYP007.EQ.1) THEN !output the different components of residuals, ratios, and solution inc.
C              WRITE(OP_STRING,
C     &          '('' Sum of solution increments              '',
C     &            ''           ='',D11.4)')
C     &          RSUM_SOLINCR
C            ELSE IF(OUTPUT.AND.(KTYP007.EQ.1)) THEN !only output sum of ratios and sum of soln inc.
C              WRITE(OP_STRING,
C     &          '('' Sum of solution increments               ''
C     &            ''                ='',D11.4)')
C     &          RSUM_SOLINCR
C              ELSE IF(KTYP007.EQ.2) THEN 
C              WRITE(OP_STRING,
C     &          '('' Sum of solution increments          ''
C     &            '' ='',D11.4)')
C     &          RSUM_SOLINCR
C            ENDIF
C            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C          ENDIF
C NEWE XSL

C NEWS 27/11/08 XSL Output sum of soln inc
          IF(OUTPUT.OR.IWRIT4(nr_solve,nx).GT.1) THEN
            IF(KTYP007.EQ.1.OR.KTYP007.EQ.2.OR.KTYP007.EQ.3) THEN
              WRITE(OP_STRING,
     &          '('' Sum of solution increments              '',
     &            ''                 ='',D11.4)')
     &          RSUM_SOLINCR
            ENDIF
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
C NEWE XSL

C***      Check Convergence
          IF(UPVU) THEN
            CALL UPVIEW(IBT,IDO,INP,ISEG,ISELNO,ISFIBR,ISFIEL,ISLINE,
     '        ISLINO,ISNONO,IWK,MXI,NAN,NBH,NBJ,NEELEM,NGAP,NHE,
     '        NHP,NKH,NKHE,NKJE,NLATNE,NLL,NLLIST,NPF,NPL,NPNE,NPNODE,
     '        NQNE,NQNLAT,
     '        NQS,NQXI,NRE,NTIW,NVHE,NVHP,NVJE,NW,nx,NYNE,
     '        NYNP,CURVCORRECT,DL,SE,XA,XE,XG,XP,XQ,YG,YP,YQ,ZA,ZE,ZP,
     '        CSEG,FIX,ERROR,*9999)
          ENDIF
 
          IF(METHOD.EQ.4) THEN
            WRITE(OP_STRING,'(/'' Solution took:'',D11.4,'
     '        //''' s CPU time; '',D11.4,'' s elapsed (wall) time.'')')
     '        TIMEA_CPU,TIMEA_ELAPSED
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE !method is not 4
C 21/07/08 XSL NEWS Add seperate convergence check for saclar product delta.R
            IF(KTYP007.EQ.4) THEN
              IF(RATIO(0).GT.ERRMAX) THEN
                CONVERGED=.FALSE.
                IF(CONTACT) THEN ! contact
                  IF(CONV_IT-1.LT.NTITER) THEN
                    SOLVE_SYSTEM=.TRUE.
                  ELSE
                    SOLVE_SYSTEM=.FALSE.
C NEWS 29/06/05 JHC Contact convergence criteria. If Newton steps taken exceeds the max.
C number of iterations specified by the user, set $CONVERGED=1 to terminate the Newton loop
C in the comfile.
C                  CALL SET_USER_INTEGER('CONVERGED',1,ERR)
                    CONV_IT=CONT_IT+1 ! add 1 to CONV_IT so that CONV_IT-1 is greater than NTITER
C NEWE JHC
                  ENDIF
                ELSE
                  IF(ITER1.LT.NTITER) THEN
                    IF(METHOD.EQ.1) THEN
                      UPDATE_MATRIX=.TRUE.
                    ELSE IF(METHOD.EQ.2) THEN
                      IF(ITER1.GT.1.AND.RATIO(0).GT.0.75D0*OLDRAT) THEN
                        UPDATE_MATRIX=.TRUE.
                      ELSE
                        UPDATE_MATRIX=.FALSE.
                      ENDIF
                    ELSE IF(METHOD.EQ.3) THEN
                      IF((ITER1.GT.1.AND.RATIO(0).GT.0.75D0*OLDRAT
     '                  .AND.RSUM_SOLINCR.GT.0.75D0*OLDRS3)) THEN
                        UPDATE_MATRIX=.TRUE.
                      ELSE
                        UPDATE_MATRIX=.FALSE.
                      ENDIF
                    ENDIF
                    SOLVE_SYSTEM=.TRUE.
                  ELSE
                    SOLVE_SYSTEM=.FALSE.
                  ENDIF !ITER1
                ENDIF !CONTACT
              ELSE !ratio(0) is below err limit
                IF((CONTACT.AND.CONV_IT.EQ.1).OR.
     &           (.NOT.CONTACT.AND.ITER1.EQ.0)) THEN !1st iteration
                  CONVERGED=.FALSE.
                  SOLVE_SYSTEM=.TRUE.
                ELSE !not 1st iteration
                  IF(CONTACT) THEN
                    CALL SET_USER_INTEGER('CONVERGED',1,ERR)
                  ENDIF !contact
                  CONVERGED=.TRUE.
                  SOLVE_SYSTEM=.FALSE.
                ENDIF !if it's the 1st iteration                          
              ENDIF !ratio(0)
            ELSE !ktyp007 is not 4
C NEWE
C            IF((RATIO(0).GT.ERRMAX.AND.RATIO(0).LT.1.0d8.AND. !not sure why we need this
              IF((RATIO(0).GT.ERRMAX.AND.
     '          RSUM_SOLINCR.GT.ZERO_TOL).OR.
     &          RSUM_SOLINCR.GT.ERRMAX) THEN
                CONVERGED=.FALSE.
                IF(CONTACT) THEN ! contact
                  IF(CONV_IT-1.LT.NTITER) THEN
                    SOLVE_SYSTEM=.TRUE.
                  ELSE
                    SOLVE_SYSTEM=.FALSE.
C NEWS 29/06/05 JHC Contact convergence criteria. If Newton steps taken exceeds the max.
C number of iterations specified by the user, set $CONVERGED=1 to terminate the Newton loop
C in the comfile.
C                  CALL SET_USER_INTEGER('CONVERGED',1,ERR)
                    CONV_IT=CONT_IT+1 ! add 1 to CONV_IT so that CONV_IT-1 is greater than NTITER
C NEWE JHC
                  ENDIF
                ELSE
                  IF(ITER1.LT.NTITER) THEN
                    IF(METHOD.EQ.1) THEN
                      UPDATE_MATRIX=.TRUE.
                    ELSE IF(METHOD.EQ.2) THEN
                      IF(ITER1.GT.1.AND.RATIO(0).GT.0.75D0*OLDRAT) THEN
                        UPDATE_MATRIX=.TRUE.
                      ELSE
                        UPDATE_MATRIX=.FALSE.
                      ENDIF
                    ELSE IF(METHOD.EQ.3) THEN
                      IF((ITER1.GT.1.AND.RATIO(0).GT.0.75D0*OLDRAT.AND.
     '                  RSUM_SOLINCR.GT.0.75D0*OLDRS3)) THEN
                        UPDATE_MATRIX=.TRUE.
                      ELSE
                        UPDATE_MATRIX=.FALSE.
                      ENDIF
                    ENDIF
                    SOLVE_SYSTEM=.TRUE.
                  ELSE
                    SOLVE_SYSTEM=.FALSE.
                  ENDIF !ITER1
                ENDIF !CONTACT
              ELSE IF(RATIO(0).GT.ERRMAX.AND.
     &          RSUM_SOLINCR.LE.ZERO_TOL) THEN
                CONVERGED=.FALSE.
C NEWS 28/06/05 JHC need to check whether it is 1st NEWTON step. 
C If so, then continue on solving, if not then the solution is not converged  
                IF(CONTACT.AND.         ! contact
     &            CONV_IT.EQ.1) THEN    ! 1st iteration
                  SOLVE_SYSTEM=.TRUE.
                ELSE IF(ITER1.EQ.0) THEN ! 1st iteration and not contact
                  SOLVE_SYSTEM=.TRUE.
                ELSE
C NEWE JHC
                  SOLVE_SYSTEM=.FALSE.
                  WRITE(CHAR1,'(D7.1)') ZERO_TOL
                  WRITE(OP_STRING,'(/'' Exiting since sum of solution '
     '              //'increments < '//CHAR1(1:7)//''')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDIF    
              ELSE
C NEWS 29/06/05 JHC Contact convergence criteria. If the problem is converged,
C then set $CONVERGED=1 to terminate the Newton loop in the comfile.
                IF(CONTACT) THEN ! contact
                  CALL SET_USER_INTEGER('CONVERGED',1,ERR)
                ENDIF
C NEWE JHC
                CONVERGED=.TRUE.
                SOLVE_SYSTEM=.FALSE.
              ENDIF !RATIO(0)
            ENDIF !ktyp007=4
C new VJ 15Apr2004. Program will never get into the if block below
C            IF(RATIO(0).GT.ERRMAX.AND.
C     &        RSUM_CONSTRAINED(0).LT.ZERO_TOL) THEN
CC              CHAR1=CFROMR(ZERO_TOL,'(D7.1)')
C              WRITE(CHAR1,'(D7.1)') ZERO_TOL
C              WRITE(OP_STRING,'(/'' Exiting since sum of constrained '
C     '          //'residuals < '//CHAR1(1:7)//''')')
C              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C            ENDIF
          ENDIF !METHOD.EQ.4

C ****    If not converged, then solve
          IF (.NOT.CONVERGED.AND.
     &       SOLVE_SYSTEM) THEN
            OLDRAT=RATIO(0) 
            OLDRS3=RSUM_SOLINCR 
            ITER1=ITER1+1
            ITER2=ITER2+1
            OUTPUT=.FALSE.
            IF(IWRIT2(nr_solve,nx).EQ.2) THEN
              IF(NWRIT.EQ.IWRIT1(nr_solve,nx)) THEN
                OUTPUT=.TRUE.
                NWRIT=1
              ELSE
                NWRIT=NWRIT+1
              ENDIF
            ENDIF

            IF(METHOD.LE.3) THEN !Full, mod or quasi Newton methods.
              IF(UPDATE_MATRIX.OR.
     &          CONTACT) THEN ! for contact problem, always assemble 
C ***           Assemble global stiffness matrix GK for region nr_solve
                CALL ASSEMBLE5(IBT,IDO,INP,ISC_GK,ISR_GK,NAN,NBH,
     '            NBHF,NBJ,NBJF,NEELEM,NFF,NFFACE,NGAP,NHE,NHP,
     '            NKB,NKEF,NKH,NKHE,NKJE,NMNO,NNB,NNF,NPF,NPNE,NPNODE,
     '            NPNY,NRE,NRLIST,nr_solve,NSB,NVHE,NVHP,NVJE,NW,nx,NXI,
     '            NYNE,NYNP,NYNR,CE,CGE,CP,CURVCORRECT,
     '            FEXT,GK,PG,SE,WG,XA,XP,
     '            YG,YGF,YP,ZA,ZA1,%VAL(0),ZP,ZP1,%VAL(0),
     '            FIX,ERROR,*9999) 
                ITER2 = 0
 
              ENDIF
C NEWS 26/06/06 JHC for contact problems, modify GK
              IF(CONTACT) THEN ! contact
                ERROR_FLAG=.FALSE.
C ***           Modify global stiffness matrix GK with contact contribution
                CALL CPU_TIMER(CPU_USER,TIME_START1)
                CALL CONTACT_STIFFNESS(IBT,IDO,INP,ISC_GK,ISR_GK,NBH,
     '            NBHF,NFF,NHE,NKEF,NKHE,NNF,NPF,NPNE,NRE,NVHE,nx,
     '            NYNP,Z_CONT_LIST,GK,SE,Z_CONT,
     '            FIX,ERROR,*400)
                GO TO 401
C ***           This statement is designed to be skipped if no error
C ***           occurs.  However if a error occurs within a subroutine
C ***           the alternate return points to line 400 to set the flag

 400            CONTINUE
                ERROR_FLAG=.TRUE.
                WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
                CALL WRITES(IOER,OP_STRING,ERROR,*401)
 401            CONTINUE
                CALL ASSERT(.NOT.ERROR_FLAG,'>>An error occurred during'
     '            //' contact stiffness calculations',ERROR,*9999)
                CALL CPU_TIMER(CPU_USER,TIME_STOP)
                ELAPSED_TIME=TIME_STOP(1)-TIME_START1(1)
                IF(IWRIT4(nr,nx).GE.1) THEN
                  WRITE(OP_STRING,
     '              '(/'' For contact stiffness calcs: CPU time of 1 '
     &              //'thread:'','
     '              //'D11.4,'' s'')') ELAPSED_TIME
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
              
C ***           New matrix
                FIRSTS(nx)=.TRUE.
                UPDATE_MATRIX=.TRUE.
              ENDIF ! CONTACT
C NEWE JHC

C 18/07/08 NEWS XSL Store GRR into temp array for line search algorithms
C The old GRR (before updating in SOLVE5) should be used to evaluate
C the previous residual at alpha=0 when CALC_CONV_RATIO is called in line search
              IF(.NOT.CONTACT) THEN
                DO no=1,NOT(1,1,nr_solve,nx) !loop over global soln rows
                  GRR_TEMP(no)=GRR(no)
                ENDDO !no
              ENDIF ! NOT CONTACT
C NEWE XSL

C ***         Solve global system of equations
              CALL SOLVE5(ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,ISR_GKK,ISR_GQ,
     '          LGE,NBH,NENP,NHE,NONY(0,1,1,nr_solve),NP_INTERFACE,NPNE,
     '          NPNY,nr_solve,NRE,NVHE,nx,NYNE,NYNO(0,1,1,nr_solve),
     '          NYNP,NYNR(0,0,1,nr_solve),CONY(0,1,1,nr_solve),
     '          CYNO(0,1,1,nr_solve),GK,GKK,GQ,
     '          GR,GRR,XO,YP,FIRSTS(nx),UPDATE_MATRIX,ERROR,*9999 )

C 28/07/08 XSL NEWS Calculate initial energy for normalising
C the scale product of the residuals and the solution increments
C Does it after SOLVE5 where solution increment gets updated
              IF((.NOT.CONTACT.AND.ITER1.EQ.1).OR. ! 1st iteration, not contact
     &          (CONTACT.AND.CONV_IT.EQ.1)) THEN ! 1st iteration, contact
                DO no=1,NOT(1,1,nr_solve,nx) !loop over global soln rows
                  DO nyo=1,NYNO(0,no,1,nr_solve)
                    ny1=NYNO(nyo,no,1,nr_solve) !is row number
                    ny2=NYNO(nyo,no,2,nr_solve) !is global variable number
                    INITIAL_ENERGY=INITIAL_ENERGY+YP(ny2,5)*YP(ny1,4)
                  ENDDO !nyo
                ENDDO !no
                INITIAL_ENERGY=DABS(INITIAL_ENERGY)
              ENDIF
C NEWE

C 04/07/08 XSL NEWS Line search should be done outside the region loop

C              DO nonrlist=1,NRLIST(0)
C                nr=NRLIST(nonrlist)
C                IF(KTYP10.EQ.2) THEN !Do line search
C                  CALL LINMIN(ALPHA,CE,CG,CGE,CONY,CP,
C     '              CURVCORRECT,FEXT,
C     '              FIX,IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,NEELEM,
C     '              NFF,NFFACE,NGAP,NHE,NHP,NKB,NKH,NKEF,NKHE,NNB,NNF,
C     '              NKJE,NONY,NPF,NPNE,NPNODE,NPNY,nr,NRE,NSB,NVHE,NVHP,
C     '              NVJE,NW,nx,NXI,NYNE,NYNO,NYNP,NYNR,PG,RE1,RG,SE,
C     '              WG,
C     '              XA,XG,XP,YP,YG,YGF,ZA,ZA1,Z_CONT,Z_CONT_LIST,ZE,ZE1,
C     '              ZP,ZP1,ERROR,*9999)
C                ELSE !no line search
C                  ALPHA=1.0d0
C                ENDIF

              IF(KTYP10.EQ.2) THEN !Do line search
                IF(.NOT.CONTACT.AND.ITER1.GT.1) THEN ! not contact, after 1st Newton iteration
                  CALL LINMIN(ALPHA,CE,CG,CGE,CONY,CP,
     '              CURVCORRECT,FEXT,FIX,GRR_TEMP,IBT,
     '              IDO,INP,ITER1,LGE,
     &             NAN,NBH,NBHF,NBJ,NBJF,NEELEM,
     '              NFF,NFFACE,NGAP,NHE,NHP,NKB,NKH,NKEF,NKHE,NNB,NNF,
     '              NKJE,NONY,NPF,NPNE,NPNODE,NPNY,NRLIST,nr_solve,
     '              NRE,NSB,NVHE,NVHP,
     '              NVJE,NW,nx,NXI,NYNE,NYNO,NYNP,NYNR,PG,RE1,
     '              REF_ENERGY,RG,SE,WG,
     '              XA,XG,XP,YP,YG,YGF,ZA,ZA1,Z_CONT,Z_CONT_LIST,ZE,ZE1,
     '              ZP,ZP1,ERROR,*9999)
                ELSE IF(CONTACT.AND.CONV_IT.GT.1) THEN !contact, after 1st Newton iteration
                  CALL LINMIN(ALPHA,CE,CG,CGE,CONY,CP,
     '              CURVCORRECT,FEXT,FIX,GRR,IBT,
     '              IDO,INP,ITER1,LGE,
     &             NAN,NBH,NBHF,NBJ,NBJF,NEELEM,
     '              NFF,NFFACE,NGAP,NHE,NHP,NKB,NKH,NKEF,NKHE,NNB,NNF,
     '              NKJE,NONY,NPF,NPNE,NPNODE,NPNY,NRLIST,nr_solve,
     '              NRE,NSB,NVHE,NVHP,
     '              NVJE,NW,nx,NXI,NYNE,NYNO,NYNP,NYNR,PG,RE1,
     '              REF_ENERGY,RG,SE,WG,
     '              XA,XG,XP,YP,YG,YGF,ZA,ZA1,Z_CONT,Z_CONT_LIST,ZE,ZE1,
     '              ZP,ZP1,ERROR,*9999)
                ELSE !1st Newton iteration
                  ALPHA=1.0d0
C XSL For highly nonlinear and unstable problames use a small initial step
C                ALPHA=0.1d0
                ENDIF
              ELSE !no line search
                ALPHA=1.0d0
              ENDIF
C XSL NEWE
              DO nonrlist=1,NRLIST(0)
                nr=NRLIST(nonrlist)
C               Update current solution by adding increments
                DO no_nynr=1,NYNR(0,0,1,nr) !loop over global vars
                  ny=NYNR(no_nynr,0,1,nr) !is global variable number
                  IF(NPNY(0,ny,0).EQ.1) THEN
                    np=NPNY(4,ny,0)
                    CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                  ENDIF
                  YP(ny,5)=ALPHA*YP(ny,5)
                  IF(NPNY(0,ny,0).EQ.1) THEN !node based ny
                    nh=NPNY(3,ny,0)
                  ELSE IF(NPNY(0,ny,0).EQ.2) THEN !element based ny
                    nh=NPNY(2,ny,0)
                  ENDIF
                  IF(KTYP5L.EQ.1.OR. !solving for deformed coordinates
     &              nh.GT.NJ_LOC(NJL_GEOM,0,nr).OR. !nh not a geometric variable
     &              NPNY(0,ny,0).EQ.2) THEN !ny is element based
                    YP(ny,1)=YP(ny,1)+YP(ny,5)
C NEWS 26/06/05 JHC for conctact problems                    
                    IF(KTYP5G(nr).GT.0) THEN !contact
                      IF (KTYP5I(nr).EQ.1) THEN ! inertia
                        YP(ny,11)=(YP(ny,1)-YP(ny,8))/T_inc ! velocity
                        YP(ny,14)=(YP(ny,11)-YP(ny,12))/T_inc ! acceleration
                      ENDIF 
                    ENDIF !KTYP5G
C NEWE JHC
                  ELSE IF(KTYP5L.EQ.2) THEN !Solving for undeformed coordinates
                    UPNODE=.TRUE.
                    nk=NPNY(1,ny,0)
                    nv=NPNY(2,ny,0)
                    nh=NPNY(3,ny,0)
                    np=NPNY(4,ny,0)
                    njj=nh !sensible default
C VJ check previous regions to see if np has another ny associated with it due to coupling
                    DO nro=1,nr
                      ny1=NYNP(nk,nv,nh,np,0,1,nro) !nrc=0 selects global variable
                      IF(ny1.LT.ny.AND.ny1.GT.0) THEN
                        UPNODE=.FALSE.
                      ENDIF
                    ENDDO
                    IF(UPNODE) THEN
                      DO nhx=1,NH_LOC(0,nx)
                        IF(nh.EQ.NH_LOC(nhx,nx)) THEN
                          njj=nhx
                        ENDIF
                      ENDDO
                      nj=NJ_LOC(NJL_GEOM,njj,nr)
C need to only update the nodal attribute once for coupled problems which have multiple degrees of freedom for one 
C nodal degree of freedom - heart wall coupled to cavity for example             
                      XP(nk,nv,nj,np)=XP(nk,nv,nj,np)+YP(ny,5) 
                    ENDIF !UPNODE
C ***  RY/MPN  14/11/06     updating cavity reference state (YP(ny,10)) when solving for undeformed
                    IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN !const vol
                      YP(ny,10)=YP(ny,10)+YP(ny,5)
                    ENDIF !ITYP2
                  ENDIF
                ENDDO !no_nynr
              ENDDO !nonrlist (nr)
C NEWS 26/06/05 JHC for conctact problems                    
              IF(CONTACT) THEN ! contact
                IF (KTYP5I(nr).EQ.1) THEN ! inertia
C ***             Copy YP(ny,11) (current solution) into ZAA,ZPA
                  DO nonrlist=1,NRLIST(0)
                    nr=NRLIST(nonrlist)
                    CALL YPZP(11,NBH,NEELEM,NHE,NHP(1,nr),NKH(1,1,1,nr),
     '                NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZAA,ZPA,
     '                ERROR,*9999)
                  ENDDO
                ENDIF

                CONT_IT=CONT_IT+1
                LOAD_IT=LOAD_IT+1
                CONV_IT=CONV_IT+1
                
                SOLVE_SYSTEM=.FALSE. ! for contact, Newton step is controlled by comfile
              ENDIF ! CONTACT
C NEWE JHC
            ELSE IF(METHOD.EQ.4) THEN !Sequential Quad Programming (MINLSSQP)
              CALL ASSERT(USE_NPSOL.GT.0,'>>Cannot use Seq. quad. '
     '          //' prog. as USE_NPSOL is zero',ERROR,*9999)
              CALL ASSERT(KTYP5L.EQ.1,'>>Cannot use Seq. quad. prog. '
     '          //' as not implemented for solving for '
     &          //' undeformed state',ERROR,*9999)
              CALL ASSERT(nr_solve.EQ.1,'>>Seq. quad. prog. needs '
     '          //'generalising for problems not in region 1 (nr needs '
     '          //'to be passed to FUNCT4, and onto UPOPTI)',ERROR,
     &          *9999)
              CALL ASSERT(NREM.GE.NOT(1,1,nr_solve,nx),
     &          '>>NREM too small',ERROR,*9999)
              CALL ASSERT(NOPM.GE.NOT(2,1,nr_solve,nx),
     &          '>>NOPM too small',ERROR,*9999)
C             Fix lower/upper bounds and initial values for soln variables
              DO no_nynr=1,NYNR(0,0,1,nr_solve) !loop over global variables
                ny=NYNR(no_nynr,0,1,nr_solve) !is global variable number
                DO noy=1,NONY(0,ny,2,nr_solve)
                  no=NONY(noy,ny,2,nr_solve)
                  PMIN(no)=-RMAX
                  PMAX(no)= RMAX
                  XC(no)=YP(ny,1)
                ENDDO
              ENDDO
              CALL MINLSSQPOPTION('Default',OPTION)
              CALL MINLSSQPOPTION('Function Tolerance = 1.0D-12',OPTION)
              CALL MINLSSQPOPTIONR('Optimality Tolerance',ERRMAX,OPTION)
              CALL MINLSSQPOPTIONL('Verbose',(KTYP10.GE.10),OPTION)
              CALL MINLSSQPOPTIONL('QP Verb',(MOD(KTYP10,20).EQ.0),
     &          OPTION)
              CALL MINLSSQPOPTIONL('Line Verb',(MOD(KTYP10,30).EQ.0),
     &          OPTION)
              CALL MINLSSQPOPTIONI('Maximum Iterations',NTITER,OPTION)
              CALL MINLSSQPOPTIONI('Derivative level',KTYP1D,OPTION)
              IF(KTYP1D.GE.0) THEN !using (some) analytic derivatives
                CALL MINLSSQPOPTIONI('Verify level',KTYP15,OPTION)
              ENDIF
              IF(DIFF_INTERVAL.GT.0.D0) THEN !diff has been set by user
                CALL MINLSSQPOPTIONR('Difference DX',DIFF_INTERVAL,
     &            OPTION)
              ENDIF
C             There are NOT(1,1,nr_solve,nx) subfunctions (residuals)
C             and NOR(nr_solve) solution vars which are held in the
C             XC      vector with BU/BL as upper/lower bounds
              CALL MINLSSQP(OBJF,XC,RESID,RESJAC,LDFJ,
     '          NOT(1,1,nr_solve,nx),NOT(2,1,nr_solve,nx),ITER,IEVAL,
     '          R,LDR,PMIN,PMAX,ISTATE,FUNCT4,IUSER,USER,OPTION,IFAIL)
              DO no_nynr=1,NYNR(0,0,1,nr_solve) !loop over global variables
                ny=NYNR(no_nynr,0,1,nr_solve) !is global variable number
                IF(NPNY(0,ny,0).EQ.1) THEN
                  np=NPNY(4,ny,0)
C GMH 8/1/97 Update cmgui link
                  CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                ENDIF
                DO noy=1,NONY(0,ny,2,nr_solve)
                  no=NONY(noy,ny,2,nr_solve)
                  YP(ny,1)=XC(no)
                ENDDO !noy
              ENDDO !no_nynr
              ITER1=iter
              CONVERGED=.TRUE.!to simply print out the timing info
            ELSE IF(METHOD.EQ.5) THEN
C KFAKFA
C Nonlinear Conjugate Gradient Method

            ENDIF ! methods

C           report CPU and ELAPSED (real) times for solutions
            TIMEB_CPU=TIMEA_CPU
            TIMEB_ELAPSED=TIMEA_ELAPSED
            CALL CPU_TIMER(CPU_USER,TIME_STOP)
            TIMEA_CPU=TIME_STOP(1)-TIME_START1(1)
            TIMED_CPU=TIMEA_CPU-TIMEB_CPU
            CALL REAL_TIMER(REAL_TOTAL,TIME_STOP)
            TIMEA_ELAPSED=TIME_STOP(1)-TIME_START1_REAL(1)
            TIMED_ELAPSED=TIMEA_ELAPSED-TIMEB_ELAPSED
            IF((OUTPUT.OR.IWRIT4(nr_solve,nx).GT.0).
     &        AND.METHOD.NE.4) THEN
              IF(CONTACT) THEN ! contact
                WRITE(OP_STRING,'(/'' --''/'' Completed iteration '
     &            //'number '','
     &            //'I3,''  ('',I3,'' since last matrix update)'')')
     '            LOAD_IT-1, 0
              ELSE
                WRITE(OP_STRING,'(/'' --''/'' Completed iteration '
     &            //'number '','
     &            //'I3,''  ('',I3,'' since last matrix update)'')')
     '            ITER1,ITER2
              ENDIF
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              IF(IWRIT4(nr_solve,nx).GT.0) THEN
                WRITE(OP_STRING,'('' Iteration time: '','
     '            //'D11.4,'' s CPU;  '',D11.4,'' s elapsed'')')
     '            TIMED_CPU,TIMED_ELAPSED
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
            ENDIF
            IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP     CRITICAL(NONLIN_2)
              WRITE(OP_STRING,'(/'' **** Post step residuals ****''/)')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP     END CRITICAL(NONLIN_2)
            ENDIF
         
          ENDIF ! CONVERGED
        ENDDO !END OF NEWTON STEP
        
        IF(CONVERGED) THEN
          IF(CONTACT) THEN ! contact
C*** 19/03/08 JHC CONV_IT should be used for writing out convergence iterations 
C            WRITE(OP_STRING,'(/'' Convergence achieved after '','
C     '        //'I3,'' iterations'')') LOAD_IT-1
            WRITE(OP_STRING,'(/'' Convergence achieved after '','
     '        //'I3,'' iterations'')') CONV_IT-1
          ELSE
            WRITE(OP_STRING,'(/'' Convergence achieved after '','
     '        //'I3,'' iterations'')') ITER1
          ENDIF
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C *** DPN 07 May 2001 - Need to be able to test for convergence in the
C *** command file, and also check for bad divergence (i.e. very large
C *** convergence ratio)
          CALL SET_USER_INTEGER('CONVERGED',1,ERR)
          CALL SET_USER_DOUBLE('CONVERGENCE_RATIO',RATIO(0),ERR)
          CALL SET_USER_INTEGER('CONVERGENCE_ITERATIONS',ITER1,ERR)
          
C NEWS 26/06/05 JHC for contact mechanics
          IF(CONTACT) THEN ! contact
C ***       Check if gap convergence is achieved.
            GAP_SUM_TIED=0.0d0
            CONT_PTS_TIED=0.0d0
            GAP_SUM_CONT=0.0d0
            CONT_PTS_CONT=0.0d0
C*** 19/03/08 JHC TOTAL_PRESS is the sum of contact pressure, which should not be 
C                 very sensitive to the penalty value used
            TOTAL_PRESS=0.0d0
                     
            DO nd=1,NDT
C ***         Check signed normal gap Z_CONT(nd,j,13) for slave (j=1)
              IF (Z_CONT_LIST(nd,1,3).EQ.1) THEN ! contact point
                IF (Z_CONT_LIST(nd,1,4).NE.2) THEN ! contact (normal or friction)
C*** 25/03/08 JHC  should check for contact pressure > 0 as contact gap might be < 0 for augmentation problems
C                  IF (Z_CONT(nd,1,13).GE.0.0d0) THEN ! then in contact
                  IF (Z_CONT(nd,1,18).GT.0.0d0) THEN ! +ve contact pressure
                    GAP_SUM_CONT=GAP_SUM_CONT+(Z_CONT(nd,1,13)**2)
                    CONT_PTS_CONT=CONT_PTS_CONT+1.0d0
                    TOTAL_PRESS=TOTAL_PRESS+Z_CONT(nd,1,18)
                  ENDIF 
                ELSE ! tied contact point                 
                  GAP_SUM_TIED=GAP_SUM_TIED+(Z_CONT(nd,1,13)**2)
                  CONT_PTS_TIED=CONT_PTS_TIED+1.0d0
                ENDIF
              ENDIF
            ENDDO !nd

            IF (CONT_PTS_CONT.EQ.0.0d0) THEN
              GAP_CHECK_CONT=0.0d0
            ELSE
C XSL Average contact gap over the amount of contact points
              GAP_CHECK_CONT=DSQRT(GAP_SUM_CONT/CONT_PTS_CONT)
            ENDIF
            IF (CONT_PTS_TIED.EQ.0.0d0) THEN
              GAP_CHECK_TIED=0.0d0
            ELSE
              GAP_CHECK_TIED=DSQRT(GAP_SUM_TIED/CONT_PTS_TIED)
            ENDIF

            WRITE(OP_STRING,'(/'' Contact Gap '',D11.4,
     '        '''')') GAP_CHECK_CONT
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
  
            WRITE(OP_STRING,'(/'' Tied Gap '',D11.4,
     '        '''')') GAP_CHECK_TIED
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C*** JHC 19/03/08 JHC Added output for total contact pressure
C                     This value should be used for sensitivity analysis on penalty value                
            WRITE(OP_STRING,'(/'' Sum Of Contact Pressure = '',E12.5)')
     &        TOTAL_PRESS
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

            AUG_IT=AUG_IT+1
            IF (AUG_IT.LE.AUGMENT) THEN
          
C ***         Reset CONVERGED status
C ***         (although only reseting CONV_IT is needed)
              CONVERGED=.FALSE.
              SOLVE_SYSTEM=.TRUE.
              CONV_IT=1 ! convergence counter         
              CALL SET_USER_INTEGER('CONVERGED',0,ERR)

C ***         Perform augmentation
C*** 19/03/08 JHC Add convergence of Lag. multipliers
              LAG_MUL_RATIO=0.0d0
              DO nd=1,NDT
                IF (AUG_IT.GT.1) THEN
                  LAG_MUL_CONV=
     &              (Z_CONT(nd,1,18)-Z_CONT(nd,1,17))**2+
     &              (Z_CONT(nd,1,20)-Z_CONT(nd,1,19))**2+
     &              (Z_CONT(nd,1,22)-Z_CONT(nd,1,21))**2
                  CURRENT_LAG_MUL=
     &              (Z_CONT(nd,1,17)**2)+
     &              (Z_CONT(nd,1,19)**2)+
     &              (Z_CONT(nd,1,21)**2)
                  IF (CURRENT_LAG_MUL.LE.ZERO_TOL.AND.
     &              LAG_MUL_CONV.GE.ZERO_TOL) THEN
                    CONT_PTS_CONT=CONT_PTS_CONT-1.0d0
                  ELSE IF (LAG_MUL_CONV.GT.ZERO_TOL) THEN
                    LAG_MUL_RATIO=LAG_MUL_RATIO+
     &                DSQRT(LAG_MUL_CONV/CURRENT_LAG_MUL)
                  ENDIF
                ENDIF
                DO j=1,2
                  Z_CONT(nd,j,17)=Z_CONT(nd,j,18) ! Update normal contact force estimate
                  Z_CONT(nd,j,19)=Z_CONT(nd,j,20) ! Update tangential1 contact force estimate
                  Z_CONT(nd,j,21)=Z_CONT(nd,j,22) ! Update tangential2 contact force estimate
                ENDDO !j
              ENDDO !nd

C*** 19/03/08 JHC Write out convergence of Lag. multipliers
              IF (AUG_IT.GT.1) THEN
                WRITE(OP_STRING,'(/'' Lagrangian multiplier '
     &            //'convergence ='',E12.5)') 
     &            LAG_MUL_RATIO/CONT_PTS_CONT
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF

              WRITE(OP_STRING,'(/'' Augmentation '',I2,
     '          '''')') AUG_IT
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

C ***         Set YP(ny,1) back to previous solution
              DO nonrlist=1,NRLIST(0) ! update initial guess
                nr=NRLIST(nonrlist)
                DO no_nynr=1,NYNR(0,0,1,nr) !loop over global vars
                  ny=NYNR(no_nynr,0,1,nr) !is global variable number
                  IF(NPNY(0,ny,0).EQ.1) THEN
                    np=NPNY(4,ny,0)
                    CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                  ENDIF
                  YP(ny,1)=YP(ny,8) ! restore previous solution 
                  YP(ny,11)=YP(ny,12) ! restore previous velocity 
                  IF (KTYP5I(nr).EQ.1) THEN ! inertia
                    YP(ny,11)=(YP(ny,1)-YP(ny,8))/T_inc ! velocity
                    YP(ny,14)=(YP(ny,11)-YP(ny,12))/T_inc ! acceleration 
                  ENDIF            
                ENDDO !no_nynr
              ENDDO !nonrlist (nr)  

            ELSE ! no augmentation or end of augmentation
              IF (KTYP5I(nr).EQ.1) THEN ! inertia
                DO nonrlist=1,NRLIST(0) ! update initial guess
                  nr=NRLIST(nonrlist)
                  DO no_nynr=1,NYNR(0,0,1,nr) !loop over global vars
                    ny=NYNR(no_nynr,0,1,nr) !is global variable number
                    IF(NPNY(0,ny,0).EQ.1) THEN
                      np=NPNY(4,ny,0)
                      CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                    ENDIF
                    YP(ny,8)=YP(ny,1) ! store previous solution
                    YP(ny,12)=YP(ny,11) ! store previous velocity     
                  ENDDO !no_nynr
                ENDDO !nonrlist (nr)
              ENDIF ! inertia

              CONV_IT=1 ! convergence counter
              LOAD_IT=1 ! reset load step counter
              AUG_IT=0 ! reset augmentation count

C*** 19/03/08 JHC Add updates for tangent forces and Lagrangian multiplers for the next load step 
              LAG_MUL_RATIO=0.0d0
              DO nd=1,NDT
                IF (AUGMENT.GT.0) THEN ! End of augmentation
                  LAG_MUL_CONV=
     &              (Z_CONT(nd,1,18)-Z_CONT(nd,1,17))**2+
     &              (Z_CONT(nd,1,20)-Z_CONT(nd,1,19))**2+
     &              (Z_CONT(nd,1,22)-Z_CONT(nd,1,21))**2
                  CURRENT_LAG_MUL=
     &              (Z_CONT(nd,1,17)**2)+
     &              (Z_CONT(nd,1,19)**2)+
     &              (Z_CONT(nd,1,21)**2)
                  IF (CURRENT_LAG_MUL.LE.ZERO_TOL.AND.
     &              LAG_MUL_CONV.GE.ZERO_TOL) THEN
                    CONT_PTS_CONT=CONT_PTS_CONT-1.0d0
                  ELSE IF (LAG_MUL_CONV.GT.ZERO_TOL) THEN
                    LAG_MUL_RATIO=LAG_MUL_RATIO+
     &                DSQRT(LAG_MUL_CONV/CURRENT_LAG_MUL)
                  ENDIF
                ENDIF

                DO j=1,2
                  Z_CONT(nd,j,19)=Z_CONT(nd,j,20) ! Update lamda_t(1), the tangnet Lagrange multiplier
                  Z_CONT(nd,j,21)=Z_CONT(nd,j,22) ! Update lamda_t(2), the tangnet Lagrange multiplier
                  IF (AUGMENT.GT.0) THEN ! End of augmentation
                    Z_CONT(nd,j,17)=Z_CONT(nd,j,18) ! Update lamda_n, the normal Lagrange multiplier
                  ENDIF ! Augment
                ENDDO !j
              ENDDO ! nd

              IF (AUGMENT.GT.0) THEN ! End of augmentation
                WRITE(OP_STRING,'(/'' Lagrangian multiplier '
     &            //'convergence ='',E12.5)')
     &            LAG_MUL_RATIO/CONT_PTS_CONT
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(/'' End Of Augmentation'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
C NEWE JHC
            
            ENDIF !AUG_IT

          ENDIF ! CONTACT
C NEWE JHC
        ELSE ! not converged
          IF(CONTACT.AND. ! contact
     &      (.NOT.SOLVE_SYSTEM)) THEN ! not converged
C NEWS 29/06/05 JHC Needs to distinguish when contact mechanics is not converged:
C If the problem has not converged within maximum number of iterations specified by user
C then set $CONVERGED in the comfile to 1 so that it exits the NEWTON loop
C IF the problem is not converged but the number of Newton steps taken so far is less than
C the maximum iterations set by user then keep solving.
            IF(CONV_IT-1.GT.NTITER) THEN
              CALL SET_USER_INTEGER('CONVERGED',1,ERR)
              WRITE(OP_STRING,
     '          '(/'' Convergence has not been reached after '','
     '          //'I3,'' iterations'')') LOAD_IT-1
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C ***         resets contact counter            
              LOAD_IT=1
              CONV_IT=1
              AUG_IT=0
 
            ELSE ! not yet converged, keep solving
C NEWE JHC
              CALL SET_USER_INTEGER('CONVERGED',0,ERR) 
            ENDIF
          ELSE ! not contact
C *** DPN 07 May 2001 - Need to be able to test for convergence in the
C *** command file, and also check for bad divergence (i.e. very large
C *** convergence ratio)
            CALL SET_USER_INTEGER('CONVERGED',0,ERR)
            CALL SET_USER_DOUBLE('CONVERGENCE_RATIO',RATIO(0),ERR)
            CALL SET_USER_INTEGER('CONVERGENCE_ITERATIONS',ITER1,ERR)
            WRITE(OP_STRING,
     '        '(/'' Convergence has not been reached after '','
C     '        //'I3,'' iterations'')') NTITER
     '        //'I3,'' iterations'')') ITER1
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF ! contact
        ENDIF !CONVERGED
            
        IF(IWRIT4(nr_solve,nx).GT.0.AND.
     &    (.NOT.CONTACT)) THEN ! for contact problems, it is not yet possible to compute the total time as those variables are not global
          WRITE(OP_STRING,'('' Total time: '','
     '      //'D11.4,'' s CPU;  '',D11.4,'' s elapsed'')')
     '      TIMEA_CPU,TIMEA_ELAPSED
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

C KAT 2005-12-14: NTACTV, ALFA, DEL_T never set
C         DO nonrlist=1,NRLIST(0)
C           nr=NRLIST(nonrlist)
C           IF(ITYP5(nr,nx).EQ.3) THEN
C             DO noelem=1,NEELEM(0,nr)
C               ne=NEELEM(noelem,nr)
C               nb=NBH(NH_LOC(1,nx),1,ne)
C               DO ng=1,NGT(nb)
C                 DO na=1,NTACTV
C                   FEXT(na+4,ng,ne) = DEXP(-ALFA(na)*DEL_T)*FEXT(na+4,ng,
C      '              ne)+DEXP(-ALFA(na)*DEL_T/2.0D0)*(FEXT(1,ng,ne)
C      '              -FEXT(2,ng,ne))
C                 ENDDO !na
C                 FEXT(2,ng,ne)=FEXT(1,ng,ne)
C               ENDDO !ng
C             ENDDO !noelem
C           ENDIF
C         ENDDO !nonrlist (nr)

      ENDDO !noload (load step)

      CALL EXITS('NONLIN')
      RETURN

C NEW 24/06/05 JHC Parallel stuff is moved to *** in Archieve directory

 9999 CALL ERRORS('NONLIN',ERROR)
      CALL EXITS('NONLIN')
      RETURN 1
      END

