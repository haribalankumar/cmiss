      SUBROUTINE IPINIT(IBT,IDO,INP,ITHRES,NAN,NBH,NBHF,NBJ,NBJF,NEELEM,
     '  NEL,NENP,NENQ,NFF,NGAP,NHE,NHP,NHQ,NKEF,NKH,NKHE,NKJE,
     '  NLL,NNF,NNL,NODENVCB,NPF,NPL,NPLIST,NP_INTERFACE,NPNE,NPNF,
     '  NPNODE,NPNY,NQET,NQLIST,NQNE,NQS,nr,NTIME_NR,NVCB,NVCNODE,
     '  NVHE,NVHF,NVHP,NVJE,NVJF,NW,NWQ,nx,NXLIST,NXI,NYNE,NYNP,NYNQ,
     '  NYNR,TV_BC_SET,AQ,CE,CG,CGE,CP,CQ,CURVCORRECT,DF,DL,PG,
     '  RCQS,RG,SE,SF,THRES,WG,XA,XE,XG,XIG,XP,XQ,YG,YP,YQ,YQS,ZA,ZE,ZG,
     '  ZP,TIME_VARIABLE_NAMES,ACTIVATION,ALL_REGIONS,FIX,FIXQ,FIX_ZERO,
     '  GENER,NOFIX,UPDATE_REF,ERROR,*)

C#### Subroutine: IPINIT
C###  Description:
C###    IPINIT inputs initial conditions.

C#### Variable: KTYP5
C###  Type: INTEGER
C###  Set_up: IPINIT
C###  Description:
C###    KTYP5 is 1/2/3 for the initial solution being Zero/Read in/
C###    Restarted or 0 if there is no initial solution.

C#### Variable: KTYP3_init
C###  Type: INTEGER
C###  Set_up: IPINIT
C###  Description:
C###    KTYP3_init specifies whether boundary condition variables are
C###    time varying:
C###      1 - constant wrt time,
C###      2 - defined in subroutine USER_IPINIT,
C###      3 - read from File.IPINIT_time at each time step,
C###      4 - defined by a time variable.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'coup00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'fluid00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'time01.cmn'
      INCLUDE 'time_variable.cmn'
      INCLUDE 'lung00.cmn'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ITHRES(3,NGM,NEM),NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),
     '  NBHF(NHM,NCM,NFM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NEELEM(0:NE_R_M,0:NRM),NEL(0:NELM,NLM),
     '  NENP(NPM,0:NEPM,0:NRM),NENQ(0:8,NQM),NFF(6,NEM),NGAP(NIM,NBM),
     '  NHE(NEM),NHP(NPM),NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NKEF(0:4,16,6,NBFM),NKH(NHM,NPM,NCM),
     '  NLL(12,NEM),NNF(0:17,6,NBFM),
     '  NNL(0:4,12,NBFM),NODENVCB(NVCBM),NPF(9,NFM),NPL(5,0:3,NLM),
     '  NPLIST(0:NP_R_M),NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),
     '  NPNF(NNM,NBFM),NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM,NXM),
     '  NQET(NQSCM),NQLIST(0:NQM),NQNE(NEQM,NQEM),NQS(NEQM),nr,
     &  NTIME_NR(0:NTIMEVARSM,NRM),NVCB(-1:3,NVCBM),NVCNODE(2,NP_R_M),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHF(NNM,NBFM,NHM),NVHP(NHM,NPM,NCM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJF(NNM,NBFM,NJM),NW(NEM,3),
     '  NWQ(8,0:NQM,NAM),nx,NXLIST(0:NXM),NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM),TV_BC_SET(0:NIQSM,0:NQM)
      REAL*8 AQ(NMAQM,NQM),CE(NMM,NEM),CG(NMM,NGM),CGE(NMM,NGM,NEM),
     '  CP(NMM,NPM),CQ(NMM,NQM),CURVCORRECT(2,2,NNM,NEM),DF(NFM),
     '  DL(3,NLM),PG(NSM,NUM,NGM,NBM),RCQS(NQRM),RG(NGM),
     '  SE(NSM,NBFM,NEM),SF(NSM,NBFM),THRES(3,NGM,NEM),WG(NGM,NBM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),XIG(NIM,NGM,NBM),
     '  XP(NKM,NVM,NJM,NPM),XQ(NJM,NQM),YG(NIYGM,NGM,NEM),YP(NYM,NIYM),
     '  YQ(NYQM,NIQM,NAM,NXM),YQS(NIQSM,NQM),ZA(NAM,NHM,NCM,NEM),
     '  ZE(NSM,NHM),ZG(NHM,NUM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),TIME_VARIABLE_NAMES(NTIMEVARSM)*(*)
      INTEGER NHQ,NYNQ(NHM,NQM,0:NRCM,NXM)

      LOGICAL ACTIVATION,ALL_REGIONS,FIX(NYM,NIYFIXM),
     '  FIXQ(NYQM,NIYFIXM,NXM),FIX_ZERO,GENER,NOFIX,UPDATE_REF


! Local variables
      INTEGER i,IBEG,IEND,n,ICHAR,INFO,NGROUP,NOQUES
      LOGICAL FILEIP
      CHARACTER CHAR*1,CHAR3*3,CHAR11*12

C DPN 24 August 1999 - time variable stuff
      INTEGER CURRENT_VAR,j,nqq,nq,VARIABLE,IBEG1,IEND1
      CHARACTER VARIABLE_NAME*(MAX_TIME_VARIABLE_NAME)
      LOGICAL NAME_FOUND

      CALL ENTERS('IPINIT',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999



      
! Time-varying parameters
      IF((ITYP5(nr,nx).EQ.2).AND.(.NOT.GENER)) THEN !time integration
        FORMAT='($,'' Are any bdry conditions time-varying [N]? '',A)'
        IF(IOTYPE.EQ.3) ADATA(1)='N'
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(ADATA(1).EQ.'Y') THEN
          FORMAT='(/'' Specify whether these coefficients are [1]:'''//
     '      '/''   (1) Constant with respect to time'''//
     '      '/''   (2) Defined in subroutine USER_IPINIT'''//
     '      '/''   (3) Read from File.IPINIT_time at each time step'''//
     '      '/''   (4) Defined by a time variable'''//
     '      '/$,''    '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=KTYP3_init(nx)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,5,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) KTYP3_init(nx)=IDATA(1)
        ELSE
          KTYP3_init(nx)=1
        ENDIF
      ELSE IF(ITYP5(nr,nx).NE.2) THEN !not time integration
        KTYP3_init(nx)=1
      ENDIF !ityp5

      IF(KTYP3_init(nx).EQ.3) THEN !time-varying from file.IPINIT_time
        FORMAT='($,'' Enter the filename [current]: '',A)'
        CDEFLT(1)=FILE00
        IF(IOTYPE.EQ.3) CDATA(1)=FILE03
        CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,30,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)

C PM 26-JUL-01
        IF(IOTYPE.NE.3) THEN
          CALL STRING_TRIM(CDATA(1),IBEG,IEND)
          FILE03=CDATA(1)(IBEG:IEND)//'.ipinit_time'
          WRITE(OP_STRING,'(2X,''File name is '',A)') FILE03
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C          WRITE(*,'('' Filename is '',A,''.ipinit'')') file03
        ENDIF
      ENDIF !ktyp3_init=3
C KSB for fluid flow in elastic tubes: time-dependent boundary
C     conditions:
      IF ((ITYP4(nr,nx).EQ.3).AND.(ITYP3(nr,nx).EQ.1).AND.
     &  KTYP3_init(nx).EQ.4) THEN
C ***     Initialisation
        CURRENT_VAR=1
        WRITE(CHAR3,'(I1)') nr
        DO j=1,2 !For each pressure boundary condition
          WRITE(CHAR,'(I1)') NTIME_NR(0,nr)+1 !variable number
          FORMAT='($,'' Enter name of time variable '//CHAR//
     &      ' for region '//CHAR3//': '',A20)'
          IF(IOTYPE.EQ.3) THEN
            CALL STRING_TRIM(TIME_VARIABLE_NAMES(NTIME_NR(CURRENT_VAR,
     &        nr)),IBEG,IEND)       
            CDATA(1)=
     &        TIME_VARIABLE_NAMES(NTIME_NR(CURRENT_VAR,nr))(IBEG:IEND)
          ENDIF
          CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,
     '      MAX_TIME_VARIABLE_NAME,IDATA,IDEFLT,IMIN,IMAX,LDATA,
     '      LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
C ***           Need to get the variable number from the name
C           WRITE(VARIABLE_NAME,'(A20)') CDATA(1)
            CALL STRING_TRIM(CDATA(1),IBEG1,IEND1)
            CALL CUPPER(CDATA(1)(IBEG1:IEND1),VARIABLE_NAME)
            NAME_FOUND=.FALSE.
C           CALL STRING_TRIM(VARIABLE_NAME,IBEG1,IEND1)
            DO i=1,NTIMEVARST
              CALL STRING_TRIM(TIME_VARIABLE_NAMES(i),IBEG,IEND)
              IF((IEND1-IBEG1).EQ.(IEND-IBEG)) THEN
                IF(VARIABLE_NAME(IBEG1:IEND1).EQ.
     '            TIME_VARIABLE_NAMES(i)(IBEG:IEND)) THEN
                  VARIABLE=i
                  NAME_FOUND=.TRUE.
                ENDIF
              ENDIF
            ENDDO !i
            IF(.NOT.NAME_FOUND) THEN
              WRITE(ERROR,'('' No such time variable: '',A20)')
     '          VARIABLE_NAME
              GOTO 9999
            ELSE
              NTIME_NR(0,nr)=NTIME_NR(0,nr)+1
              NTIME_NR(NTIME_NR(0,nr),nr)=VARIABLE
            ENDIF
          ENDIF
        ENDDO
      ENDIF
C PM 26-JUL-01 : for fluid in elastic tube
      IF ((ITYP4(nr,nx).EQ.3).AND.(ITYP3(nr,nx).EQ.1)) GOTO 6410

      IF(KTYP3_init(nx).GT.1) THEN !time-varying from subroutine or file
        !Note: a better way of doing this would be to use another
        !loop in ipinix and ask for the time varying ny's --> FIX

C       ..fe60 problems - get systole and diastole outlet nodes
        IF(ITYP4(nr,nx).EQ.5.AND.ITYP2(nr,nx).EQ.5.AND.
     '    ITYP3(nr,nx).EQ.3) THEN

          NGROUP=1
          IDEFLT(1)=NGROUP
          IF(IOTYPE.EQ.3) IDATA(1)=NGROUP
          FORMAT='($,'' The number of node groups for the '//
     '      'systole outlet is [1]: '',I2)'
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,
     '      FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     '      IDATA,IDEFLT,1,100,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     '      RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NGROUP=IDATA(1)
          SYSTOLE_OUTLET(0)=0
          DO i=1,NGROUP
            FORMAT='($,'' Enter the node numbers of the '//
     '        'systole outlet: '',I5)'
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,
     '        FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     '        IDATA,IONE,1,NPM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     '        RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              DO n=1,IDATA(0)
                SYSTOLE_OUTLET(SYSTOLE_OUTLET(0)+n)=IDATA(n)
              ENDDO !n
              SYSTOLE_OUTLET(0)=SYSTOLE_OUTLET(0)+IDATA(0)
            ENDIF !iotype.ne.3
          ENDDO !ngroup

          NGROUP=1
          IDEFLT(1)=NGROUP
          IF(IOTYPE.EQ.3) IDATA(1)=NGROUP
          FORMAT='($,'' The number of node groups for the '//
     '      'diastole outlet is [1]: '',I2)'
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,
     '      FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     '      IDATA,IDEFLT,1,100,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     '      RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NGROUP=IDATA(1)
          DIASTOLE_OUTLET(0)=0
          DO i=1,NGROUP
            FORMAT='($,'' Enter the node numbers of the '//
     '        'diastole outlet: '',I5)'
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,
     '        FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     '        IDATA,IONE,1,NPM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     '        RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              DO n=1,IDATA(0)
                DIASTOLE_OUTLET(DIASTOLE_OUTLET(0)+n)=IDATA(n)
              ENDDO !n
              DIASTOLE_OUTLET(0)=DIASTOLE_OUTLET(0)+IDATA(0)
            ENDIF !iotype.ne.3
          ENDDO !ngroup

          RDEFLT(1)=1.d0
          WRITE(CHAR11,'(E12.4)') RDEFLT(1)
          CALL STRING_TRIM(CHAR11,IBEG,IEND)
          FORMAT='($,'' The cyclical (heart) beat period is '
     '      //'['//CHAR11(IBEG:IEND)//']: '',E12.4)'
          IF(IOTYPE.EQ.3) RDATA(1)=HEART_PERIOD
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.ne.3) THEN
            HEART_PERIOD=RDATA(1)
          ENDIF

          WRITE(OP_STRING,'(''>>Remember to code up USER_10'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

        ELSE IF(KTYP3_init(nx).EQ.4) THEN

C ***     Want to use a time variable to specify a b.c.

C          CALL ASSERT(TV_NUM_VARIABLES.GT.0,
C     '      '>>Need to define time variables first',ERROR,*9999)
          CALL ASSERT(CALL_TIME,'>>Must define time first',ERROR,*9999)
          CALL ASSERT(USE_CELL.EQ.1,'>>USE_CELL must be 1',ERROR,*9999)

C ***     Initialisation
          CURRENT_VAR=1
          DO nqq=1,NIQST
            DO nq=1,NQT
              TV_BC_SET(nqq,nq)=0
            ENDDO
          ENDDO
          IF(IOTYPE.NE.3) TV_BC_SET(0,0)=0
C ***     First get the variable number - loop until exit
 6301     FORMAT='($,'' Enter the variable number [EXIT]: '',I5)'
          IF(IOTYPE.EQ.3) THEN
            IF(TV_BC_SET(0,0).GE.CURRENT_VAR) THEN
              IDATA(1)=TV_BC_SET(CURRENT_VAR,0)
            ELSE
              IDATA(1)=0
            ENDIF
          ENDIF
C         DPN dodgy hack - using 0-100 until ability to control non
C         state variables is fully implemented...
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,
     '      0,NIQST,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '      ERROR,*9999)
C          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
C     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,
C     '      0,100,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
C     '      ERROR,*9999)
          IF(IDATA(1).NE.0) THEN !not default exit
            IF(IOTYPE.NE.3) THEN
              TV_BC_SET(0,0)=TV_BC_SET(0,0)+1
              TV_BC_SET(CURRENT_VAR,0)=IDATA(1)
            ENDIF

C ***       Next get the grid points
C DPN 26 July 2001 - Need to use NQLIST but not IDATA (its too small)
            WRITE(CHAR3,'(I3)') TV_BC_SET(CURRENT_VAR,0)
            nqq=NQR(1,nr)
 6302       FORMAT='($,'' Enter collocation point #s/name '
     '        //'[EXIT]: '',I5)'
            IF(IOTYPE.EQ.3) THEN
              IF(nqq.LE.NQR(2,nr)) THEN
                !IDATA(1)=nqq
                NQLIST(0)=1
                NQLIST(1)=nqq
              ELSE
                !IDATA(0)=0
                !IDATA(1)=0
                NQLIST(0)=0
                NQLIST(1)=0
              ENDIF
              nqq=nqq+1
            ENDIF
 6402       CDATA(1)='GRIDS' !for use with group input
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,NQLIST,IZERO,
     '        0,NQT,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '        ERROR,*9999)
            IF(NQLIST(1).NE.0) THEN !not default exit
              IF(IOTYPE.NE.3) THEN
                !NQLIST(0)=IDATA(0)
                !DO n=1,IDATA(0)
                DO n=1,NQLIST(0)
                  !NQLIST(n)=IDATA(n)
                  !nq=IDATA(n)
                  nq=NQLIST(n)
                  IF(nq.LT.NQR(1,nr).OR.nq.GT.NQR(2,nr)) THEN
                    WRITE(OP_STRING,'('' >>Grid point '',I7,'' is not '
     '                //'in the current region'')') nq
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 6402
                  ENDIF
                ENDDO !n
              ENDIF !IOTYPE.NE.3
C             Define variable for first grid point in group. Rest of
C             group filled in at the end.
              nq=NQLIST(1)
              CALL STRING_TRIM(TIME_VARIABLE_NAMES(1),IBEG,IEND)
              CDEFLT(1)=TIME_VARIABLE_NAMES(1)(IBEG:IEND)
              CALL STRING_TRIM(CHAR3,IBEG1,IEND1)
              FORMAT='($,'' The time variable to use for variable '
     '          //CHAR3(IBEG1:IEND1)//' is ['
     '          //TIME_VARIABLE_NAMES(1)(IBEG:IEND)//']: '',A20)'
              IF(IOTYPE.EQ.3) CDATA(1)=
     '          TIME_VARIABLE_NAMES(TV_BC_SET(CURRENT_VAR,nq))
              CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,
     '          MAX_TIME_VARIABLE_NAME,IDATA,IDEFLT,IMIN,IMAX,LDATA,
     '          LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
C ***           Need to get the variable number from the name
                WRITE(VARIABLE_NAME,'(A20)') CDATA(1)
                NAME_FOUND=.FALSE.
C                VARIABLE=1
C                DO WHILE(VARIABLE.LE.TV_NUM_VARIABLES.AND.
C     '            .NOT.NAME_FOUND)
C                  CALL STRING_TRIM(TV_NAMES(VARIABLE),IBEG,IEND)
C                  CALL CUPPER(TV_NAMES(VARIABLE),BUFFER)
C                  IF(CBBREV(VARIABLE_NAME,BUFFER,IEND-IBEG,IBEG,IEND,
C     '              N3CO)) THEN
C                    NAME_FOUND=.TRUE.
C                  ELSE
C                    VARIABLE=VARIABLE+1
C                  ENDIF
C                ENDDO
                CALL STRING_TRIM(VARIABLE_NAME,IBEG1,IEND1)
                DO i=1,NTIMEVARST
                  CALL STRING_TRIM(TIME_VARIABLE_NAMES(i),IBEG,IEND)
                  IF((IEND1-IBEG1).EQ.(IEND-IBEG)) THEN
                    IF(VARIABLE_NAME(IBEG1:IEND1).EQ.
     '                TIME_VARIABLE_NAMES(i)(IBEG:IEND)) THEN
                      VARIABLE=i
                      NAME_FOUND=.TRUE.
                    ENDIF
                  ENDIF
                ENDDO !i
                IF(.NOT.NAME_FOUND) THEN
                  WRITE(ERROR,'('' No such time variable: '',A20)')
     '              VARIABLE_NAME
                  GOTO 9999
                ELSE
                  TV_BC_SET(CURRENT_VAR,nq)=VARIABLE
                ENDIF
C               Apply to all elements in the group
                DO n=2,NQLIST(0)
                  nq=NQLIST(n)
                  TV_BC_SET(CURRENT_VAR,nq)=VARIABLE
                ENDDO !n
              ENDIF
              GOTO 6302 !For more grid points
            ENDIF
            CURRENT_VAR=CURRENT_VAR+1
            GOTO 6301 !for more variables
          ENDIF
        ELSE
          FORMAT='($,'' Enter the number of nodes (up to 99)[1]: '',I2)'
          IF(IOTYPE.EQ.3) IDATA(1)=NODE_BC_time(0)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,99,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NODE_BC_time(0)=IDATA(1)

          DO n=1,NODE_BC_time(0)
            WRITE(CHAR3,'(I3)') n
            FORMAT='($,'' Enter node# ['//CHAR3//']:'',I6)'
            IDEFLT(1)=n
            IF(IOTYPE.EQ.3) IDATA(1)=NODE_BC_time(n)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &        NPM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) NODE_BC_time(n)=IDATA(1)
          ENDDO !n

        ENDIF !ityp
      ENDIF !ktyp3_init

 6410 IF(ITYP1(nr,nx).EQ.3) THEN
        CALL IPINI3(IBT,IDO,INP,ITHRES,NBH,NBHF,NBJ,NBJF,NEELEM,
     '    NEL,NENP,NENQ,NFF,NHE,NHP,NHQ,NKEF,NKH,NKHE,NKJE,NNF,NPF,
     '    NPLIST,NLL,NP_INTERFACE,NNL,NPL,NPNE,NPNF,NPNODE,NQET,
     '    NQLIST,NQNE,NQS,nr,NVHE,NVHF,
     '    NVHP,NVJE,NVJF,NW,NWQ,nx,NXLIST,NXI,NYNE,NYNP,NYNQ,
     '    NYNR(0,0,1,nr),AQ,CE,CQ,CURVCORRECT,DL,PG,RCQS,RG,SE,
     '    SF,THRES,WG,XA,XE,XG,XIG,XP,XQ,YG,YP,YQ,YQS,ZA,ZE,ZP,
     '    TIME_VARIABLE_NAMES,ALL_REGIONS,FIX,FIXQ,FIX_ZERO,GENER,
     '    ERROR,*9999)
      ELSE IF(ITYP1(nr,nx).EQ.4) THEN
        CALL IPINI4(IDO,INP,NBH,NBHF,NBJ,NEELEM,NEL,NENP,NFF,
     '    NGAP,NHE,NHP,NKEF,NKH,NKHE,NLL,NNF,NNL,NPF,NPL,NPLIST,NPNE,
     '    NPNODE,nr,NVHE,NVHP,NW,nx,NYNE,NYNP,DF,DL,PG,SE,
     '    WG,XP,YP,ZA,ZP,FIX,NOFIX,ERROR,*9999)
      ELSE IF(ITYP1(nr,nx).EQ.5) THEN
C news VJ added YG array to IPINI5 for grid coupling gauss stress feature
        CALL IPINI5(IBT,IDO,INP,NAN,NBH,NBHF,NEELEM,NEL,NENP,NFF,NHE,
     '    NHP,NKEF,NKH,NKHE,NLL,NNF,NNL,NPF,NPL,NPLIST,NPNE,NPNODE,nr,
     '    NVHE,NVHP,NW,nx,NYNE,NYNP,CE,CG,CGE,CP,DF,DL,PG,SE,WG,XE,XG,
     '    XP,YG,YP,ZA,ZE,ZG,ZP,FIX,NOFIX,UPDATE_REF,ERROR,*9999)
C newe VJ
      ELSE IF(ITYP1(nr,nx).EQ.6) THEN
!new rgb 1/5/98
        CALL IPINI6(NHP,NKH,NODENVCB,NPNODE,nr,NVCB,NVCNODE,
     '    NVHP,nx,NYNP,YP,ALL_REGIONS,ERROR,*9999)
!     SMAR009 23/12/98 removed ,NODENVC from list
      ELSE IF(ITYP1(nr,nx).EQ.9) THEN
!news AJP 12/4/95
        IF(ITYP5(nr,nx).EQ.1.AND.ITYP2(nr,nx).EQ.1) THEN
          !Static linear elasticity
          CALL IPINI4(IDO,INP,NBH,NBHF,NBJ,NEELEM,NEL,NENP,NFF,
     '      NGAP,NHE,NHP,NKEF,NKH,NKHE,NLL,NNF,NNL,NPF,NPL,NPLIST,NPNE,
     '      NPNODE,nr,NVHE,NVHP,NW,nx,NYNE,NYNP,DF,DL,PG,SE,
     '      WG,XP,YP,ZA,ZP,FIX,NOFIX,ERROR,*9999)
        ELSE !rest of BEM
CC AJPs 191297
          CALL IPINI9(IBT,NBH,NEELEM,NHP,NKH,NPLIST,NP_INTERFACE,
     '      NPNODE,NPNY,nr,NVHP,NW,nx,NYNP,NYNR(0,0,1,nr),XP,YP,
     '      ACTIVATION,ALL_REGIONS,FIX,FIX_ZERO,GENER,NOFIX,ERROR,*9999)
CC AJPe 191297
        ENDIF
      ENDIF
!newe

      CALL EXITS('IPINIT')
      RETURN
 9999 CALL ERRORS('IPINIT',ERROR)
      CALL EXITS('IPINIT')
      RETURN 1
      END


