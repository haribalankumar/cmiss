      SUBROUTINE FITSIG(IBT,IDO,INP,ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,
     '  ISR_GKK,ISR_GQ,LD,LGE,LN,NBH,NBJ,
     '  NDDATA,NDDL,NDLT,NEELEM,NENP,
     '  NHE,NHP,NHQ,NKB,NKH,NKHE,NKJE,NLL,
     '  NNB,NNF,NNL,NONY,NPF,NP_INTERFACE,NPL,NPNE,NPNODE,
     '  NPNY,NQNY,nr,NRE,NRLIST,NRLIST2,NVHE,NVJE,NVHP,
     '  NW,NWP,nx,NXI,NYNE,
     '  NYNO,NYNP,NYNR,NYNY,NYQNR,CONY,CURVCORRECT,
     '  CYNO,CYNY,EDD,ER,ES,GK,GKK,GR,
     '  GRR,GQ,PG,RG,SE,SF,SP,TEND,TSTART,WD,WDL,WG,WU,XA,XE,XG,XID,
     '  XIDL,XO,XP,YP,YQ,YQS,ZA,ZD,ZDL,ZE,ZP,FILEFORMAT,FIRST_A,
     '  FIX,IN_PLANE,ERROR,*)

C#### Subroutine: FITSIG
C###  Description:
C###    FITSIG fits field variables to signal data.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'data00.cmn'
      INCLUDE 'fit000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'sign00.cmn'
      INCLUDE 'solv00.cmn'
      INCLUDE 'time02.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISC_GK(NISC_GKM),ISC_GKK(NISC_GKKM),ISC_GQ(NISC_GQM),
     '  ISR_GK(NISR_GKM),ISR_GKK(NISR_GKKM),ISR_GQ(NISR_GQM),
     '  LD(NDM),LGE(NHM*NSM,NRCM),LN(0:NEM),NBH(NHM,NCM,NEM),
     '  NBJ(NJM,NEM),
     '  NDDATA(0:NDM,0:NRM),NDDL(NEM,NDEM),NDLT(NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NENP(NPM,0:NEPM,0:NRM),NHE(NEM),
     '  NHP(NPM),NHQ(NRM),NKB(2,2,2,NNM,NBFM),NKH(NHM,NPM,NCM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),
     '  NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),
     '  NONY(0:NOYM,NYM,NRCM,0:NRM),NPF(9,NFM),NP_INTERFACE(0:NPM,0:3),
     '  NPL(5,0:3,NLM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM),NQNY(2,NYQM,0:NRCM),nr,
     '  NRE(NEM),NRLIST(0:NRM),NRLIST2(0:NRM),
     '  NVHE(NNM,NBFM,NHM,NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVHP(NHM,NPM,NCM),NW(NEM,3),nx,
     '  NWP(NPM,2),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNO(0:NYOM,NOOPM,NRCM,0:NRM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM),NYNY(0:NYYM,NYM,NRM),
     '  NYQNR(0:NYQM,0:NRCM,NCM,0:NRM)
      REAL*8 CONY(0:NOYM,NYM,NRCM,0:NRM),CURVCORRECT(2,2,NNM,NEM),
     '  CYNO(0:NYOM,NOOPM,NRCM,0:NRM),CYNY(0:NYYM,NYM,NRM),
     '  EDD(NDM),ER(NHM*NSM),ES(NHM*NSM,NHM*NSM),GK(NZ_GK_M),
     '  GKK(NZ_GKK_M),GR(NYROWM),GRR(NOM),GQ(NZ_GQ_M),
     '  PG(NSM,NUM,NGM,NBM),RG(NGM),
     '  SE(NSM,NBFM,NEM),SF(NSM,NBFM),
     '  SP(NKM,NBFM,NPM),TEND,TSTART,WD(NJM,NDM),
     '  WDL(NHM,NDEM),WG(NGM,NBM),WU(0:NUM+1,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),XID(NIM,NDM),
     '  XIDL(NIM,NDEM),XO(NOM),XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM),
     '  YQ(NYQM,NIQM,NAM),YQS(NIQSM,NQM),
     '  ZA(NAM,NHM,NCM,NEM),ZD(NJM,NDM),ZDL(NHM,NDEM),ZE(NSM,NHM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),FILEFORMAT*(*)
      LOGICAL FIRST_A,FIX(NYM,NIYFIXM),IN_PLANE
!     Local Variables
      INTEGER GETNYR,l,na,nb,nd,ndl,NTOT,ne,nh,nhj,nhs1,
     '  nhs2,NHST(2),nhx,NIQLIST(0:1),NIQSLIST(0:1),NIYLIST(0:1),
     '  nj,njj,nk,no1,no2,nonode,no_nynr1,no_nynr2,notime,no_times1,
     '  no_times2,noy1,noy2,np,NUMTIMEDATA,NUMTIMEDATA1,nv,ny1,ny2,
     '  ny3,nyo1,nz,nzz
      INTEGER*4 WORK_PTR
      REAL AVETIME1,AVETIME2,ELAPSED_TIME,TIME_START1(1),TIME_START2(1),
     '  TIME_START3(1),TIME_START4(1),TIME_START5(1),TIME_STOP(1)
      REAL*8 co1,co2,co3,PXI,SAED,SADED,
     '  SIGNALMAX(NSIGNALREGIONSMX),SIGNALMIN(NSIGNALREGIONSMX),
     '  SMED,SMDED,SQDED,TIME,YPMAX(1),YPMIN(1),Z(6)
      LOGICAL ENDFILE,FIRSTTIME,NEXTTIME,UPDATE_MATRIX,X_INIT,YPDATA,
     &  YQDATA,YQSDATA
      CHARACTER CERROR*50
C      CHARACTER CERROR(50)

C     31-OCT-97 INTEGER err, unused
C               LOGICAL ISBINFILEOPEN,


      CALL ENTERS('FITSIG',*9999)
C LKC 15-OCT-97 added start timers
      CALL CPU_TIMER(CPU_USER,TIME_START1)


C LKC 12-SEP-1999 This is not correct, the data in a signal
C file always starts (internally) from 1. The region
C numbers are ignored. The DATA_REGION variable attempts
C to map between the fitted region and the signal region.
C
CC LKC 24-MAY-1998
C      CALL ASSERT(nr.EQ.DATA_REGION,'>>Fit region NE Data region'
C     '  ,ERROR,*9999)



C LKC 14-JAN-98 Initialise to some value
      TIME=0.D0

C***  Open up signal file
      CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '  SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILENAME,
     '  'OPEN',ENDFILE,.TRUE.,ERROR,*9999)

      CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '  SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILENAME,
     '  'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)

C LKC 15-OCT-97 IF(ASCII) IOSIGN does not calculate the number of
C   time data points.

C***  Calculate the number of time data points
      IF(FILEFORMAT.EQ.'ASCII') THEN
        NUMTIMEDATA=0
        ENDFILE=.FALSE.
        DO WHILE(.NOT.ENDFILE)
          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILENAME,
     '      'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
          IF(.NOT.ENDFILE) NUMTIMEDATA=NUMTIMEDATA+1
        ENDDO
      ENDIF

      CALL ASSERT(NUMTIMEDATA.GT.0,'>>No signal data',ERROR,*9999)


      CALC_XI=SIGNAL_ELEMLOC(IOFILE1).EQ.1
      CALL ASSERT(CALC_XI,'>>Define xi positions first',ERROR,*9999)

C***  Calculate element data arrays
      CALL LDLDR(LD,NDDATA,NDDL,NDLT,ERROR,*9999)

C***  Open up the history file
      NRLIST(0)=1
      NRLIST(1)=nr
      NRLIST2(0)=1
      NRLIST2(1)=nr
      NIYLIST(0)=1
      NIYLIST(1)=1
      NIQLIST(0)=0
      NIQSLIST(0)=0
      NEXTTIME=.TRUE.
      YPDATA=.TRUE.
      YQDATA=.FALSE.
      YQSDATA=.FALSE.
      na=1
      CALL IOHIST(IOFILE2,na,NHQ,NIQLIST,NIQSLIST,NIYLIST,NPNY,NQNY,
     '  NRLIST,NRLIST2,NUMTIMEDATA1,nx,NYNR,NYQNR,TIME,YP,YPMAX,YPMIN,
     '  YQ,YQS,'WRITE',FILEFORMAT,OUTFILENAME,'OPEN',ENDFILE,NEXTTIME,
     '  YPDATA,YQDATA,YQSDATA,ERROR,*9999)

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START1(1)
      IF(IWRIT4(nr,nx).GE.1) THEN
        WRITE(OP_STRING,'(/'' CPU time for signal and history '
     '    //'initialisation: '',D11.4,'' s'')') ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

C***  Loop over the number of signals

      NO_TIMES1=0
      AVETIME1=0.0
      FIRSTTIME=.TRUE.

C LKC 15-OCT-97
C***    Reset both files
      CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '  SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILENAME,
     '  'RESET',ENDFILE,.TRUE.,ERROR,*9999)
      CALL ASSERT(NUMTIMEDATA.GT.0,
     '  '>>No signal data, NUMTIMEDATA = 0',ERROR,*9999)

C LKC 1-SEPT-1999 Moved from the inside the time loop
      IF(SIGNAL_SKIPSAMPLE.GT.1) THEN
        NUMTIMEDATA=NUMTIMEDATA-1
      ENDIF

      DO notime=1,NUMTIMEDATA
        CALL CPU_TIMER(CPU_USER,TIME_START2)

C***    Read in the signal data
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILENAME,
     '    'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)

        IF(SIGNAL_SKIPSAMPLE.GT.1.AND.notime.EQ.1) THEN
C*** fix so the times are still correct when using sampleskip
          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,INFILENAME,
     '      'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)

C LKC 1-SEP-1999 moved outside time loop
C          NUMTIMEDATA=NUMTIMEDATA-1

        ENDIF

        IF(TIME.GT.TEND) GOTO 9000
        IF( MOD(notime,SIGNAL_SKIPSAMPLE).EQ.0
     '    .AND.TIME.GE.TSTART
     '    .AND.TIME.LE.TEND) THEN

C LKC 18-OCT added "at"
          IF(SIGNAL_OUTPUT.GE.2) THEN
            WRITE(OP_STRING,'(/'' Fitting signal '',I5,'//
     '        ' '' at time ='',D12.4)') notime,TIME
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF

C cpb 8/11/95 Haven't put in the njj loop yet as it is not needed at
C moment
          njj=1

         CALL CPU_TIMER(CPU_USER,TIME_START3)

C cpb 5/7/96 Zero ZP each time

          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            DO nhj=1,NUM_FIT(njj)
              nj=NLH_FIT(nhj,1,njj)
              nhx=NLH_FIT(nhj,3,njj)
              nh=NH_LOC(nhx,nx)
              DO nv=1,NVHP(nh,np,1)
                DO nk=1,NKH(nh,np,1)
                  ny1=NYNP(nk,nv,nh,np,0,1,nr)
                  IF(FIX(ny1,1)) THEN
                    ZP(nk,nv,nh,np,1)=XP(nk,nv,nj,np)
                  ELSE
                    ZP(nk,nv,nh,np,1)=0.0d0
                  ENDIF
                ENDDO !nk
              ENDDO !nv
            ENDDO !nhj
          ENDDO !nonode (np)

          CALL ZPYP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,nr,NVHP,nx,
     '      NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)

          IF(FIRSTTIME) THEN !first time is being fitted
            UPDATE_MATRIX=.TRUE.
            FIRSTTIME=.FALSE.
          ELSE
            UPDATE_MATRIX=.FALSE.
          ENDIF

C***      Initialise global variables
          CALL CPU_TIMER(CPU_USER,TIME_START4)
          IF(UPDATE_MATRIX) THEN
            CALL INIT_SPARSE_MATRIX(NYNR(0,2,1,nr),NYNR(0,2,1,nr),
     '        ISC_GK,ISR_GK,0,0,1,NYT(1,1,nx),NZ_GK_M,NZT(1,nx),
     '        NYNR(0,1,1,nr),NYNR(0,1,1,nr),KTYP24,GK,ERROR,*9999)
          ENDIF
          DO no_nynr1=1,NYNR(0,1,1,nr)
            ny1=NYNR(no_nynr1,1,1,nr)
            GR(ny1)=0.d0
          ENDDO !no_nynr1 (ny1)

          CALL CPU_TIMER(CPU_USER,TIME_STOP)
          ELAPSED_TIME=TIME_STOP(1)-TIME_START4(1)
          IF(IWRIT4(nr,nx).GE.1) THEN
            WRITE(OP_STRING,'(/'' CPU time for global matrix setup and '
     '        //'initialisation: '',D11.4,'' s'')') ELAPSED_TIME
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF

          CALL CPU_TIMER(CPU_USER,TIME_START4)
          NO_TIMES2=0
          AVETIME2=0.0

          DO l=1,LN(0) !loop over elements in the fit

            CALL CPU_TIMER(CPU_USER,TIME_START5)

            ne=LN(l)

            CALL MELGEF(LGE,NBH(1,1,ne),ne,NHST,njj,NPNE(1,1,ne),
     '        nr,NVHE(1,1,1,ne),nx,NYNE,NYNP,ERROR,*9999)

            IF(IWRIT4(nr,nx).GE.5) THEN
              FORMAT='(/'' Element'',I5,'', Number of variables: '','
     '          //'''NHST(1)='',I3,'', NHST(2)='',I3)'
              WRITE(OP_STRING,FORMAT) ne,NHST(1),NHST(2)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              FORMAT='('' LGE(1..,1): '',14I5,:(/13X,14I5))'
              WRITE(OP_STRING,FORMAT) (LGE(nhs1,1),nhs1=1,NHST(1))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              FORMAT='('' LGE(1..,2): '',14I5,:(/13X,14I5))'
              WRITE(OP_STRING,FORMAT) (LGE(nhs1,2),nhs1=1,NHST(2))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            DO nhs1=1,NHST(1)
              ER(nhs1)=0.0d0
              IF(UPDATE_MATRIX) THEN
                DO nhs2=1,NHST(2)
                  ES(nhs1,nhs2)=0.0d0
                ENDDO !nhs2
              ENDIF
            ENDDO !nhs1

            CALL ZPZE(NBH(1,1,ne),1,NUM_FIT(njj),NKHE(1,1,1,ne),
     '        NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1),nx,
     '        CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,
     '        ERROR,*9999)

            CALL ZDER(IBT,IDO,INP,NBH(1,1,ne),NBJ(1,ne),NDDL,NDLT(ne),
     '        ne,njj,NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '        NRE(ne),NVJE(1,1,1,ne),nx,ER,PG,RG,
     '        SE(1,1,ne),SF,WD,WDL,WG,WU(0,ne),XA(1,1,ne),XE,XG,XID,
     '        XIDL,XP,ZD,ZDL,ZE,ZP,IN_PLANE,ERROR,*9999)

            IF(UPDATE_MATRIX) THEN
              CALL ZDES(IBT,IDO,INP,NBH(1,1,ne),NBJ(1,ne),NDLT(ne),njj,
     '          NKJE(1,1,1,ne),NPF,NPNE(1,1,ne),NRE(ne),
     '          NVJE(1,1,1,ne),nx,ES,PG,RG,
     '          SE(1,1,ne),SF,WDL,WG,WU(0,ne),XA(1,1,ne),XE,XG,XIDL,
     '          XP,IN_PLANE,ERROR,*9999)
            ENDIF

            IF(IWRIT4(nr,nx).GE.5) THEN
              WRITE(OP_STRING,'(/'' Element '',I4,'' rhs vector '
     '          //'and stiffness matrix:'')') ne
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              CALL OPESTFMAT(NHST,IOOP,ES,ER,'ES','ER',.TRUE.,.TRUE.,
     '          ERROR,*9999)
            ENDIF

C***        Assemble element stiffness matrix into global system.

            DO nhs1=1,NHST(1)
              ny1=IABS(LGE(nhs1,1))
              GR(ny1)=GR(ny1)+ER(nhs1)
              IF(UPDATE_MATRIX) THEN
                DO nhs2=1,NHST(2)
                  ny2=IABS(LGE(nhs2,2))
                  CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,NZ_GK_M,
     '              NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
                  GK(nz)=GK(nz)+ES(nhs1,nhs2)
                ENDDO !nhs2
              ENDIF
            ENDDO !nhs1

            CALL CPU_TIMER(CPU_USER,TIME_STOP)
            ELAPSED_TIME=TIME_STOP(1)-TIME_START5(1)
            AVETIME2=AVETIME2+ELAPSED_TIME
            NO_TIMES2=NO_TIMES2+1
            IF(IWRIT4(nr,nx).GE.1) THEN
              WRITE(OP_STRING,'('' CPU time for element '',I5,'
     '          //''' assembly: '',D11.4,'' s'')') ne,ELAPSED_TIME
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF

          ENDDO !l (ne)

          CALL CPU_TIMER(CPU_USER,TIME_STOP)
          ELAPSED_TIME=TIME_STOP(1)-TIME_START4(1)
          IF(IWRIT4(nr,nx).GE.1) THEN
            WRITE(OP_STRING,'(/'' CPU time for stiffness matrix '
     '        //'assembly : '',D11.4,'' s'')') ELAPSED_TIME
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' Average CPU time per element for '
     '        //'assembly:'',D11.4,'' s'')') AVETIME2/NO_TIMES2
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF

          IF(IWRIT4(nr,nx).GE.4) THEN
            CALL CPU_TIMER(CPU_USER,TIME_START4)
            IF(UPDATE_MATRIX) THEN
              WRITE(OP_STRING,
     '          '(/'' Global load vector GR & stiffness matrix GK:'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' NYNR(0,1,1,nr)='',I5,'
     '          //''', NYNR(0,2,1,nr)='',I5)') NYNR(0,1,1,nr),
     '          NYNR(0,2,1,nr)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              CALL OPSTFMAT(NYNR(0,2,1,nr),ISC_GK,ISR_GK,IOOP,
     '          NYT(1,1,nx),NYT(2,1,nx),NZT(1,nx),NYNR(0,1,1,nr),KTYP24,
     '          GK,GR,'GK ','GR ',.FALSE.,.TRUE.,.TRUE.,ERROR,*9999)
            ELSE
              WRITE(OP_STRING,
     '          '(/'' Global load vector GR:'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' NYNR(0,1,1,nr)='',I5)')
     '          NYNR(0,1,1,nr)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              CALL OPSTFMAT(NYNR(0,2,1,nr),ISC_GK,ISR_GK,IOOP,
     '          NYT(1,1,nx),NYT(2,1,nx),NZT(1,nx),NYNR(0,1,1,nr),KTYP24,
     '          GK,GR,'GK ','GR ',.FALSE.,.FALSE.,.TRUE.,
     '          ERROR,*9999)
            ENDIF

            CALL CPU_TIMER(CPU_USER,TIME_STOP)
            ELAPSED_TIME=TIME_STOP(1)-TIME_START4(1)
            WRITE(OP_STRING,'(/'' CPU time for stiffness matrices '
     '        //'output: '',D11.4,'' s'')') ELAPSED_TIME
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

          ENDIF

        CALL CPU_TIMER(CPU_USER,TIME_START4)

C*** Calculate solution mapping arrays for the current fit variable

          IF(UPDATE_MATRIX) THEN
            CALL GLOBALF(IBT,IDO,INP,NBH,NBJ,NENP,1,NKB,NKH,NKHE,NKJE,
     '        NLL,NNB,NNF,NNL,NONY,NPF,NPL,NPNE,
     '        NPNODE,NPNY,nr,NVHE,NVHP,NVJE,NWP,
     '        nx,NXI,NYNE,NYNO,NYNP,NYNR,NYNY,CONY,
     '        CYNO,CYNY,SE,SP,XA,XE,XP,FIX,ERROR,*9999)

            IF(NOT(2,1,nr,nx).EQ.0) THEN
              ERROR=' >>The number of unknowns is zero'
              GOTO 9999
            ENDIF

C*** Initialise solution variables

            IF(SPARSEGKK(nx).NE.0) THEN
              WORK_PTR=0
              CALL ALLOCATE_MEMORY(NOT(1,1,nr,nx)*NOT(2,1,nr,nx),1,
     '          CHARTYPE,WORK_PTR,MEM_INIT,ERROR,*9999)
            ENDIF
            CALL CALC_SPARSE_GKK(ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,ISR_GKK,
     '        ISR_GQ,LGE,NBH,NENP,NHE,NOT(1,1,nr,nx),NOT(2,1,nr,nx),
     '        NONY(0,1,1,nr),NP_INTERFACE,NPNE,NPNY,nr,NRE,NVHE,nx,
     '        NYNE,NYNP,NYNR(0,1,1,nr),GK,GQ,%VAL(WORK_PTR),
     '        .FALSE.,.TRUE.,ERROR,*9999)
            IF(SPARSEGKK(nx).NE.0) THEN
              CALL FREE_MEMORY(WORK_PTR,ERROR,*9999)
            ENDIF
            DO nzz=1,NZZT(1,nr,nx)
              GKK(nzz)=0.0d0
            ENDDO !nzz

          ENDIF
          DO no1=1,NOT(1,1,nr,nx)
            GRR(no1)=0.0d0
          ENDDO !no1

          DO no1=1,NOT(2,1,nr,nx)
            DO nyo1=1,NYNO(0,no1,2,nr)
              ny1=NYNO(nyo1,no1,2,nr)
              IF(NPNY(0,ny1,0).EQ.1) THEN
                np=NPNY(4,ny1,0)
C GMH 8/1/97 Update cmgui link
                CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
              ENDIF
              YP(ny1,1)=0.0d0
            ENDDO !nyo1
          ENDDO !no1

          CALL CPU_TIMER(CPU_USER,TIME_STOP)
          ELAPSED_TIME=TIME_STOP(1)-TIME_START4(1)
          IF(IWRIT4(nr,nx).GE.1) THEN
            WRITE(OP_STRING,'(/'' CPU time for solution matrix setup '
     '        //'and initialisation: '',D11.4,'' s'')') ELAPSED_TIME
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF

C*** Generate the reduced system of equations
          CALL CPU_TIMER(CPU_USER,TIME_START4)

c cpb 25/9/95 Putting fitting into standard form
          DO no_nynr1=1,NYNR(0,1,1,nr) !Loop over global rows of GK
            ny1=NYNR(no_nynr1,1,1,nr) !is row #
            DO noy1=1,NONY(0,ny1,1,nr) !loop over #no's attached to ny1
              no1=NONY(noy1,ny1,1,nr) !is no# attached to row ny1
              co1=CONY(noy1,ny1,1,nr) !is coupling coeff for row mappng
C                                    ie row_no1=a*row_ny1+b*row_ny2
              GRR(no1)=GRR(no1)+GR(ny1)*co1 !get reduced R.H.S.vector
              DO no_nynr2=1,NYNR(0,0,1,nr) !loop over the #cols of GK
                ny2=NYNR(no_nynr2,0,1,nr) !is global variable #
                ny3=GETNYR(1,NPNY,nr,2,0,ny2,NYNE,NYNP) !local GK var #
                CALL SPARSE(ny1,ny3,NYT(1,1,nx),nz,NZ_GK_M,
     '            NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
                IF(nz.NE.0) THEN
                  IF(UPDATE_MATRIX) THEN
                    DO noy2=1,NONY(0,ny2,2,nr) !loop over #no's for ny2
                      no2=NONY(noy2,ny2,2,nr) !no# attached to ny2
                      co2=CONY(noy2,ny2,2,nr) !coup coeff for col mapng
C                                     i.e. var_no1=a*var_ny1+b*var_ny2
                      CALL SPARSE(no1,no2,NOT(1,1,nr,nx),nzz,
     '                  NZ_GKK_M,NZZT(1,nr,nx),ISC_GKK,ISR_GKK,
     '                  SPARSEGKK(nx),ERROR,*9999)
                      IF(nzz.NE.0) GKK(nzz)=GKK(nzz)+GK(nz)*co1*co2
                    ENDDO !noy2
                  ENDIF
C!!! There is only one constant defined for each variable therefore
C!!! it is implied that this constant is applied to the last ny mapped
C!!! to the no. ie. no=a*ny1+b*ny2+(c*ny3+d) that is a*ny1=b*ny2 and
C!!! a*ny1=c*ny3+d and not a*ny1=b*ny2+d etc.
                  co3=CONY(0,ny2,2,nr) !is add constant applied to vars
                  GRR(no1)=GRR(no1)-GK(nz)*co3 !add const. in RHS vec
                ENDIF
              ENDDO !no_nynr2
            ENDDO !noy1
          ENDDO !no_nynr1

          CALL CPU_TIMER(CPU_USER,TIME_STOP)
          ELAPSED_TIME=TIME_STOP(1)-TIME_START4(1)
          IF(IWRIT4(nr,nx).GE.1) THEN
            WRITE(OP_STRING,'(/'' CPU time for solution matrix '
     '        //'assembly: '',D11.4,'' s'')') ELAPSED_TIME
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF


          IF(KTYP4.NE.0) THEN !Output global matrices
            CALL WRITE_SOL_MATRIX(ISC_GKK,ISR_GKK,nr,nx,GKK,GRR,
     '        ERROR,*9999)
          ENDIF

          IF(IWRIT4(nr,nx).GE.3) THEN
            CALL CPU_TIMER(CPU_USER,TIME_START4)
            WRITE(OP_STRING,
     '        '(/'' Global load vector GRR & stiffness matrix GKK:'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' NOT(1,1,nr,nx)='',I5,'
     '        //''', NOT(2,1,nr,nx)='',I5)') NOT(1,1,nr,nx),
     '        NOT(2,1,nr,nx)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            CALL OPSTFMAT(NYNR(0,1,1,nr),ISC_GKK,ISR_GKK,IOOP,
     '        NOT(1,1,nr,nx),NOT(2,1,nr,nx),NZZT(1,nr,nx),
     '        NYNR(0,2,1,nr),SPARSEGKK(nx),GKK,GRR,'GKK','GRR',
     '        .TRUE.,.TRUE.,.TRUE.,ERROR,*9999)
          ENDIF

          CALL CPU_TIMER(CPU_USER,TIME_STOP)
          ELAPSED_TIME=TIME_STOP(1)-TIME_START4(1)
          IF(IWRIT4(nr,nx).GE.1) THEN
            IF(IWRIT4(nr,nx).GE.4.OR.KTYP4.NE.0) THEN
              WRITE(OP_STRING,'(/'' CPU time for solution matrix '
     '          //'output: '',D11.4,'' s'')') ELAPSED_TIME
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF

C*** Solve reduced system

          ! Use zero as an initial solution guess for the first iteration.
          ! Otherwise use the solution from the previous time.
          ! A separate variable to FIRST_A is required because FIRST_A is
          ! modified in SOLVE_SYSTEM.
          X_INIT=FIRST_A

          CALL SOLVE_SYSTEM(ISC_GKK,ISR_GKK,NOT(1,1,nr,nx),
     '      NOT(1,1,nr,nx),NOT(2,1,nr,nx),NZZT(1,nr,nx),IWRIT4(nr,nx),
     '      PRECON_CODE(nx),SOLVEROPTION(nx),SPARSEGKK(nx),GKK,GRR,
     '      XO,FIRST_A,UPDATE_MATRIX,X_INIT,nx,ERROR,*9999)

          DO no1=1,NOT(2,1,nr,nx)
            DO nyo1=1,NYNO(0,no1,2,nr)
              ny1=NYNO(nyo1,no1,2,nr)
              IF(NPNY(0,ny1,0).EQ.1) THEN
                np=NPNY(4,ny1,0)
C GMH 8/1/97 Update cmgui link
                CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
              ENDIF
              co1=CYNO(nyo1,no1,2,nr)
              YP(ny1,1)=YP(ny1,1)+XO(no1)*co1
            ENDDO !nyo1
          ENDDO !no1

          CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,nr,NVHP,nx,
     '      NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)

C***      Write the history file
          NEXTTIME=.TRUE.
          YPDATA=.TRUE.
          YQDATA=.FALSE.
          na=1
          CALL IOHIST(IOFILE2,na,NHQ,NIQLIST,NIQSLIST,NIYLIST,NPNY,
     '      NQNY,NRLIST,NRLIST2,NUMTIMEDATA1,nx,NYNR,NYQNR,TIME,YP,
     '      YPMAX,YPMIN,YQ,YQS,'WRITE',FILEFORMAT,OUTFILENAME,
     '      'TIME_DATA',ENDFILE,NEXTTIME,YPDATA,YQDATA,YQSDATA,
     '      ERROR,*9999)

          CALL CPU_TIMER(CPU_USER,TIME_START4)

          SMED=0.0d0
          SAED=0.0d0
          SQED=0.0d0
          SMDED=0.0d0
          SADED=0.0d0
          SQDED=0.0d0
          NTOT=0
          DO l=1,LN(0)
            ne=LN(l)
            CALL ZPZE(NBH(1,1,ne),1,NUM_FIT(njj),NKHE(1,1,1,ne),
     '        NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1),
     '        nx,CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),
     '        ZE,ZP,ERROR,*9999)
            DO nhj=1,NUM_FIT(njj)
              nhx=NLH_FIT(nhj,3,njj)
              nh=NH_LOC(nhx,nx)
              nb=NBH(nh,1,ne)
              nj=NJ_FIT(nhj,njj)
              DO ndl=1,NDLT(ne)
                nd=NDDL(ne,ndl)
                Z(nh)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '            XID(1,nd),ZE(1,nhx))
                EDD(nd)=Z(nh)-ZD(nj,nd)
                SMED=SMED+EDD(nd)
                SAED=SAED+DABS(EDD(nd))
                SQED=SQED+EDD(nd)**2
                SMDED=SMDED+ZD(nj,nd)
                SADED=SADED+DABS(ZD(nj,nd))
                SQDED=SQDED+ZD(nj,nd)**2
                NTOT=NTOT+1
              ENDDO !ndl
            ENDDO !nhj
          ENDDO !l (ne)
C LEO 18-OCT-97 restricted output
          IF(NTOT.GT.1.AND.SIGNAL_OUTPUT.GE.3) THEN
            WRITE(OP_STRING,'(/'' Data: '','
     '        //'/''   Average value           : '',D12.6,'' +/- '','
     '        //'D12.6,'
     '        //'/''   Average absolute value  : '',D12.6,'' +/- '','
     '        //'D12.6,'
     '        //'/''   Root mean squared value : '',D12.6,'
     '        //'/'' Fit: '','
     '        //'/''   Average error           : '',D12.6,'' +/- '','
     '        //'D12.6,'
     '        //'/''   Average absolute error  : '',D12.6,'' +/- '','
     '        //'D12.6,'
     '        //'/''   Root mean squared error : '',D12.6)')
     '        SMDED/DBLE(NTOT),DSQRT((SQDED-SMDED**2/DBLE(NTOT))/
     '        DBLE(NTOT-1)),
     '        SADED/DBLE(NTOT),DSQRT((SQDED-SADED**2/DBLE(NTOT))/
     '        DBLE(NTOT-1)),
     '        DSQRT(SQDED/DBLE(NTOT)),
     '        SMED/DBLE(NTOT),DSQRT((SQED-SMED**2/DBLE(NTOT))/
     '        DBLE(NTOT-1)),
     '        SAED/DBLE(NTOT),DSQRT((SQED-SAED**2/DBLE(NTOT))/
     '        DBLE(NTOT-1)),
     '        DSQRT(SQED/DBLE(NTOT))
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF

          CALL CPU_TIMER(CPU_USER,TIME_STOP)
          ELAPSED_TIME=TIME_STOP(1)-TIME_START4(1)
          IF(IWRIT4(nr,nx).GE.1) THEN
            IF(IWRIT4(nr,nx).GE.4.OR.KTYP4.NE.0) THEN
              WRITE(OP_STRING,'(/'' CPU time for error estimate: '','
     '          //'D11.4,'' s'')') ELAPSED_TIME
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF

          CALL CPU_TIMER(CPU_USER,TIME_STOP)
          ELAPSED_TIME=TIME_STOP(1)-TIME_START3(1)
          IF(IWRIT4(nr,nx).GE.1) THEN
            IF(IWRIT4(nr,nx).GE.4.OR.KTYP4.NE.0) THEN
              WRITE(OP_STRING,'(/'' Total CPU time for fit variable '','
     '          //'I1,'': '',D11.4,'' s'')') njj,ELAPSED_TIME
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF
        ENDIF !time >= tstart & time <= tend

        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
        AVETIME1=AVETIME1+ELAPSED_TIME
        NO_TIMES1=NO_TIMES1+1

        IF( SIGNAL_OUTPUT.GE.1
     '   .AND. MOD(notime,SIGNAL_SKIPSAMPLE).EQ.0
     '   .AND. TIME.GE.TSTART
     '   .AND. TIME.LE.TEND) THEN
          WRITE(OP_STRING,'('' Signal '',I5,'//
     '      ''' completed ('',D11.4, ''). '//
     '      ' CPU time: '',D11.4,'' s'')')
     '      notime,TIME,ELAPSED_TIME
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

      ENDDO !notime

 9000 CONTINUE

C LKC 15-OCT-97 Replace with 'proper' close
C      CALL BINCLOSEFILE(IOFILE1,ERROR,*9999)
C      CALL BINCLOSEFILE(IOFILE2,ERROR,*9999)

C***  Close the history and signal file
      CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '  SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,INFILENAME,
     '  ' ',ENDFILE,.TRUE.,ERROR,*9998)
      NEXTTIME=.TRUE.
      YPDATA=.TRUE.
      YQDATA=.FALSE.
      YQSDATA=.FALSE.
      na=1
      CALL IOHIST(IOFILE2,na,NHQ,NIQLIST,NIQSLIST,NIYLIST,NPNY,NQNY,
     '  NRLIST,NRLIST2,NUMTIMEDATA1,nx,NYNR,NYQNR,TIME,YP,YPMAX,YPMIN,
     '  YQ,YQS,'CLOSE',FILEFORMAT,OUTFILENAME,' ',ENDFILE,NEXTTIME,
     '  YPDATA,YQDATA,YQSDATA,ERROR,*9999)

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START1(1)

C LKC 31-OCT-97 redirect output
C added the second CALL WRITS (but not WRITE(OP_STRING ....)
      IF(SIGNAL_OUTPUT.GE.1) THEN
        WRITE(OP_STRING,'(/'' Total CPU time for signal fitting '
     '    //'problem: '',D11.4,'' s'')') ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Average CPU time per signal: '',D11.4,'
     '    //''' s'')') AVETIME1/NO_TIMES1
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('FITSIG')
      RETURN
 9999 CALL ERRORS('FITSIG',ERROR)

C LKC 15-OCT-97 Close the files 'properly' below
C      IF(ISBINFILEOPEN(IOFILE1))
C     '  CALL BINARYCLOSEFILE(IOFILE1,ERR,CERROR)
C      IF(ISBINFILEOPEN(IOFILE2))
C     '  CALL BINARYCLOSEFILE(IOFILE2,ERR,CERROR)

C***  Close the history and signal file
      CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '  SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,INFILENAME,
     '  ' ',ENDFILE,.TRUE.,CERROR,*9998)
 9998 NEXTTIME=.TRUE.
      YPDATA=.TRUE.
      YQDATA=.FALSE.
      YQSDATA=.FALSE.
      na=1
      CALL IOHIST(IOFILE2,na,NHQ,NIQLIST,NIQSLIST,NIYLIST,NPNY,NQNY,
     '  NRLIST,NRLIST2,NUMTIMEDATA1,nx,NYNR,NYQNR,TIME,YP,YPMAX,
     '  YPMIN,YQ,YQS,'CLOSE',FILEFORMAT,OUTFILENAME,' ',ENDFILE,
     '  NEXTTIME,YPDATA,YQDATA,YQSDATA,CERROR,*9997)
 9997 CALL EXITS('FITSIG')
      RETURN 1
      END


