      SUBROUTINE DEEQUA(IBT,IDO,INP,LGE,NBH,NBHF,NBJ,NEELEM,NELIST,
     '  NENP,NFF,NHE,NHP,NHQ,NKH,NKHE,NKJE,NLL,NNF,NPF,
     '  NP_INTERFACE,NPL,NPNE,NPNODE,NPNY,NQNY,NRLIST,NVHE,NVHP,
     '  NVJE,NVJP,NW,NXI,NXLIST,NYNE,NYNP,NYNQ,NYNR,NYQNR,CE,
     '  CURVCORRECT,SE,XA,XAB,XE,XP,FIX,STRING,ERROR,*)

C#### Subroutine: DEEQUA
C###  Description:
C###    DEEQUA defines equations.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'ptr00.cmn'
      INCLUDE 'solv00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LGE(NHM*NSM,NRCM),NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NENP(NPM,0:NEPM,0:NRM),NFF(6,NEM),
     '  NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),
     '  NHQ(NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),
     '  NNF(0:17,6,NBFM),NPF(9,NFM),NP_INTERFACE(0:NPM,0:3),
     '  NPL(5,0:3,NLM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),NRLIST(0:NRM),
     '  NQNY(2,NYQM,0:NRCM,NXM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),NW(NEM,3,NXM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),NXLIST(0:NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNQ(NHM,NQM,0:NRCM,NXM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),
     '  NYQNR(0:NYQM,0:NRCM,NCM,0:NRM,NXM)
      REAL*8 CE(NMM,NEM,NXM),CURVCORRECT(2,2,NNM,NEM),SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XAB(NORM,NEM),XE(NSM,NJM),XP(NKM,NVM,NJM,NPM)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,
     '  IPFILE,N3CO,ne,nee,no_nrlist,nr,nx,nx2,nxc,nxx,NX_ACTION
      CHARACTER FILE*(MXCH),STATUS*3
      LOGICAL ALL_REGIONS,CALCU,CBBREV,FILIO,
     '  FIRST_TIME,GENER,MOUSE

      CALL ENTERS('DEEQUA',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define equation;l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Define the equation type, solution method, and, if necessary,
C###    dependent variable interpolation type. These are read from or
C###    written to the file FILENAME.ipequa in the directory PATH.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Parameter:      <(lock/nolock)[lock]>
C###    Specify whether or not to lock this class number (of solve
C###    type).
C###  Parameter:      <compliance (lambert/vessel/none)>
C###    Specify the type of compliance relationship in flow equations        
CC###  Parameter:      <element (#s/all)[all]>


        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
C        OP_STRING(2)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(4)=BLANK(1:15)//'<(lock/nolock)[lock]>'
        OP_STRING(4)=BLANK(1:15)//'<(compliance (lambert/vessel)[none]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C        OP_STRING(1)=STRING(1:IEND)//';m'
C        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
C        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEEQUA',ERROR,*9999)
      ELSE
        IPFILE=2 !is input file version number on 7-Dec-1995
C GMH 7/12/95 1->2 Method of asking for element types changed

        CALL PARSE_QUALIFIERS('DLPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
C KAT 6Nov98: Not used anywhere.
c        CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
c     '    ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL ASSERT(NXM.GE.NXLIST(0),'>>Increase NXM in ippara',
     '    ERROR,*9999)

        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

        IF(CBBREV(CO,'NOLOCK',2,noco+1,NTCO,N3CO)) THEN
          NX_ACTION=NX_ALLOCATE
        ELSE
          NX_ACTION=NX_ALLOCATE_AND_LOCK
        ENDIF

        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        IF(nx.EQ.0) THEN
          CALL NX_LOC(NX_ACTION,nxc,nx,NX_SOLVE,ERROR,*9999)
        ENDIF

        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          CALL ASSERT(nr.LE.NRT,'>>Region not defined',ERROR,
     '        *9999)
          CALL ASSERT(NPT(nr).GT.0,'>>No nodes defined',ERROR,
     '      *9999)
          CALL ASSERT(NET(nr).GT.0,'>>No elements defined',ERROR,
     '      *9999)
        ENDDO

        IF(FILIO) THEN
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'equa',STATUS,
     '        ERR,ERROR,*9999)

            IF(ISC_GK_PTR(nx).EQ.0) THEN
              CALL ALLOCATE_MEMORY(NISC_GKM*USE_SPARSE,0,INTTYPE,
     '        ISC_GK_PTR(nx),MEM_INIT,ERROR,*9999)
            ENDIF
            IF(ISR_GK_PTR(nx).EQ.0) THEN
              CALL ALLOCATE_MEMORY(NISR_GKM*USE_SPARSE,0,INTTYPE,
     '        ISR_GK_PTR(nx),MEM_INIT,ERROR,*9999)
            ENDIF
            IF(ISC_GQ_PTR(nx).EQ.0) THEN
              CALL ALLOCATE_MEMORY(NISC_GQM*USE_SPARSE,0,INTTYPE,
     '        ISC_GQ_PTR(nx),MEM_INIT,ERROR,*9999)
            ENDIF
            IF(ISR_GQ_PTR(nx).EQ.0) THEN
              CALL ALLOCATE_MEMORY(NISR_GQM*USE_SPARSE,0,INTTYPE,
     '        ISR_GQ_PTR(nx),MEM_INIT,ERROR,*9999)
            ENDIF

            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              IF(ALL_REGIONS) CALL PROMPT_REGION_ALL(nr,ERROR,*9999)
              CALL IPEQUA(IBT,IDO,INP,%VAL(ISC_GK_PTR(nx)),
     '          %VAL(ISC_GQ_PTR(nx)),%VAL(ISR_GK_PTR(nx)),
     &          %VAL(ISR_GQ_PTR(nx)),LGE,NBH,NBHF,NBJ,NEELEM,NELIST,
     &          NENP,NFF,NHE(1,nx),NHP(1,0,nx),NHQ(1,nx),NKH,NKHE,NKJE,
     &          NLL,NNF,NPF,NP_INTERFACE,NPL,NPNE,NPNODE,NPNY(0,1,0,nx),
     '          NQNY(1,1,0,nx),nr,NRLIST,NVHE,NVHP,NVJE,NVJP,NW(1,1,nx),
     '          nx,NXI,NYNE,NYNP,NYNQ(1,1,0,nx),NYNR(0,0,1,0,nx),
     '          NYQNR(0,0,1,0,nx),CE(1,1,nx),CURVCORRECT,SE,XA,XAB,XE,
     '          XP,FIX(1,1,nx),ERROR,*9999)
            ENDDO !no_nrlist
            CALL CLOSEF(IFILE,ERROR,*9999)
            IF(BACKUP) IOTYPE=2
          ENDDO !while first_time or backup
        ENDIF !filio

        IF(IOTYPE.NE.3) THEN
          DO nr=1,NRT
            ASSEMBLE_GLOBAL(nr,nx)=.FALSE.
          ENDDO !nr
        ENDIF

        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          IF(ITYP4(nr,nx).EQ.2) THEN !BE method
            DO nee=1,NEELEM(0,nr)
              ne=NEELEM(nee,nr)
              NW(ne,3,nx)=0
            ENDDO
          ENDIF
        ENDDO

        IF(KTYP32.EQ.2) THEN
          !Bidomain problem defined
          CALL ASSERT(NXLIST(0).GE.2,
     '      '>>Error - must define two classes for bidomain problems',
     '      ERROR,*9999)
          DO nxx=2,NXLIST(0)
            nxc=NXLIST(nxx)
            CALL NX_LOC(NX_INQUIRE,nxc,nx2,NX_SOLVE,ERROR,*9999)
            IF(nx2.EQ.0) THEN
              CALL NX_LOC(NX_ACTION,nxc,nx2,NX_SOLVE,ERROR,*9999)
            ENDIF
            ITYP5(NRLIST(1),nx2)=2 !time dep.
            ITYP4(NRLIST(1),nx2)=4 !collocation - no information regarding
                                   ! choice of method is known at this point
                                   ! defaults to collocation for bidomain MLT
            ITYP3(NRLIST(1),nx2)=ITYP3(NRLIST(1),nx) !model
            ITYP2(NRLIST(1),nx2)=9 !cellular based modelling
            ITYP19(NRLIST(1),nx2)=1 !electrical modelling
            ITYP16(NRLIST(1),nx2)=2 !explicit
            NH_LOC(1,nx2)=NH_LOC(1,nx)
            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              NHQ(nr,nx2)=NHQ(nr,nx)
              CALL CALC_NY_GRID_DEP(NHQ(1,nx2),NQNY(1,1,0,nx2),nr,nx2,
     '          NYNQ(1,1,0,nx2),NYQNR(0,0,1,0,nx2),ERROR,*9999)
C MLT 13Nov2002  initialising NHE and NBH for bidomain second class
              DO nee=1,NEELEM(0,nr)
                ne=NEELEM(nee,nr)
                NHE(ne,nx2) = NHE(ne,nx)
                NBH(NH_LOC(1,nx2),1,ne) = NBH(NH_LOC(1,nx),1,ne)
              ENDDO
            ENDDO
          ENDDO
C MLT 13Nov2002 Unsure what function of following code is
C SGM 26Oct2000 grid-based Finite element also
C MLT 29Nov02 Added grid Finite Volumes for completeness
        ELSE IF((ITYP4(NRLIST(1),nx).EQ.4.OR.
     '           ITYP4(NRLIST(1),nx).EQ.6.OR.
     '           ITYP4(NRLIST(1),nx).EQ.7)
     '      .AND.(NXLIST(0).GT.1)) THEN
          DO nxx=2,NXLIST(0)
            nxc=NXLIST(nxx)
            CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
            IF(nx.EQ.0) THEN
              CALL NX_LOC(NX_ACTION,nxc,nx,NX_SOLVE,ERROR,*9999)
            ENDIF
          ENDDO
        ENDIF

        CALL_EQUA=.TRUE.
      ENDIF

      CALL EXITS('DEEQUA')
      RETURN
 9999 CALL ERRORS('DEEQUA',ERROR)
      IF(nx.GT.0) THEN
        CALL NX_LOC(NX_FREE,nxc,nx,NX_SOLVE,ERROR,*9999)
      ENDIF
      CALL EXITS('DEEQUA')
      RETURN 1
      END


