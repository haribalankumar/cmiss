      SUBROUTINE LIOUTP(IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,NEELEM,
     &  NELIST,NENP,NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,NKHE,NKJE,NLL,
     &  NNB,NNF,NP_INTERFACE,NPF,NPNE,NPNODE,NPNY,NRE,NRLIST,NSB,
     '  NVHE,NVHP,NVJE,NVJP,NW,NXI,NXLIST,NYNE,NYNP,NYNR,Z_CONT_LIST,
     '  CE,CG,CGE,CP,CURVCORRECT,DET,DL,DRDN,FEXT,PG,RAD,
     '  RD,RE1,RG,SE,
     '  WG,XA,XE,XG,XG1,XIG,XN,XP,XR,YD,YG,YGF,YP,ZA,ZA1,
     '  Z_CONT,ZE,
     '  ZF,ZG,ZP,ZP1,STRING,FIX,ERROR,*)
      
C#### Subroutine: LIOUTP
C###  Description:
C###    LIOUTP outputs solutions.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'pulm00.cmn'
      INCLUDE 'ityp00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LGE(NHM*NSM),NAN(NIM,NAM,NBFM),
     '  NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),NENP(NPM,0:NEPM,0:NRM),
     &  NFF(6,NEM),NFFACE(0:NF_R_M,NRM),NGAP(NIM,NBM),NHE(NEM,NXM),
     &  NHP(NPM,0:NRM,NXM),NKB(2,2,2,NNM,NBFM),NKEF(0:4,16,6,NBFM),
     &  NKH(NHM,NPM,NCM,0:NRM),NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,
     &  NEM),NLL(12,NEM),NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),
     &  NP_INTERFACE(0:NPM,0:3),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     &  NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM,NXM),NRE(NEM),
     &  NRLIST(0:NRM),NSB(NKM,NNM,NBFM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),
     '  NW(NEM,3,NXM),NXI(-NIM:NIM,0:NEIM,0:NEM),NXLIST(0:NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),Z_CONT_LIST(NDM,2,7)
      REAL*8 CE(NMM,NEM,NXM),CG(NMM,NGM),CGE(NMM,NGM,NEM,NXM),
     '  CP(NMM,NPM,NXM),CURVCORRECT(2,2,NNM,NEM),
     '  DET(NBFM,0:NNM,NGM,6),DL(3,NLM),DRDN(NGM),
     '  FEXT(NIFEXTM,NGM,NEM),PG(NSM,NUM,NGM,NBM),RAD(NGM),
     '  RD(NGM),RE1(NSM,NHM),
     '  RG(NGM),SE(NSM,NBFM,NEM),WG(NGM,NBM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),XG1(NJM,NUM,NGM),
     '  XIG(NIM,NGM,NBM),XN(NJM,NGM),XP(NKM,NVM,NJM,NPM),
     '  XR(NJM,NGM),YD(NHM),
     '  YG(NIYGM,NGM,NEM),YGF(NIYGFM,NGFM,NFM),YP(NYM,NIYM,NXM),
     '  ZA(NAM,NHM,NCM,NEM),ZA1(NAM,NHM,NCM,NEM),Z_CONT(NDM,2,67),
     '  ZE(NSM,NHM),ZF(NSM,NHM),ZG(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM),ZP1(NKM,NVM,NHM,NPM,NCM)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,N3CO,no_nrlist,nr,nx,nxc
      REAL*8 ZE1(NSM,NHM)
      CHARACTER FILE*100
      LOGICAL ALL_REGIONS,CBBREV,OPFILE

      CALL ENTERS('LIOUTP',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list output<;FILENAME>
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Lists solution output to screen or file FILENAME.opoutp if
C###    qualifier present.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(5)=BLANK(1:15)//'<for pulmonary blood flow:>'
        OP_STRING(6)=BLANK(1:15)//'      <INLET #[NEELEM(1)]>'
        OP_STRING(7)=BLANK(1:15)//'      <OUTLET #[NEELEM(NEELEM(0))]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIOUTP',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opoutp','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,ERROR,
     &    *9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        nr=NRLIST(1)
        IF(ITYP1(nr,nx).EQ.3.AND.ITYP2(nr,nx).EQ.11) THEN
        !pulmonary transport problems (currently only blood flow!)
          IF(CBBREV(CO,'INLET',3,noco+1,NTCO,N3CO))THEN
            INLETS(1)=IFROMC(CO(N3CO+1)) !specified inlet vessel
            INLETS(0)=1
          ENDIF
          IF(CBBREV(CO,'OUTLET',3,noco+1,NTCO,N3CO))THEN
            OUTLETS(1)=IFROMC(CO(N3CO+1)) !specified inlet vessel
            OUTLETS(0)=1
          ENDIF
          IF(INLETS(0).EQ.0) THEN
            INLETS(0)=1
            INLETS(1)=NEELEM(1,nr) !default = 1st element
          ENDIF
          IF(OUTLETS(0).EQ.0) THEN
            OUTLETS(0)=1
            OUTLETS(1)=NEELEM(NEELEM(0,nr),nr) !default = last element
          ENDIF
        ENDIF
        
        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          CALL OPOUTP(IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,
     '      NEELEM,NELIST,NENP(1,0,nr),NFF,
     '      NFFACE,NGAP,NHE(1,nx),NHP(1,nr,nx),NKB,NKEF,
     '      NKH(1,1,1,nr),NKHE,NKJE,NLL,NNB,NNF,NP_INTERFACE,NPF,NPNE,
     '      NPNODE,NPNY(0,1,0,nx),nr,NRE,NSB,NVHE,NVHP(1,1,1,nr),NVJE,
     '      NVJP,NW(1,1,nx),nx,NXI,NYNE,NYNP,NYNR(0,0,1,nr,nx),
     '      Z_CONT_LIST,
     '      CE(1,1,nx),CG,CGE(1,1,1,nx),CP(1,1,nx),
     '      CURVCORRECT,DET,DL,DRDN,FEXT,PG,RAD,RD,RE1,RG,SE,
     '      WG,XA,XE,
     '      XG,XG1,XIG,XN,XP,XR,YD,YG,YGF,YP(1,1,nx),ZA,ZA1,Z_CONT,ZE,
     '      ZE1,ZF,ZG,ZP,ZP1,FIX(1,1,nx),ERROR,*9999)
        ENDDO !no_nrlist

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LIOUTP')
      RETURN
 9999 CALL ERRORS('LIOUTP',ERROR)
      CALL EXITS('LIOUTP')
      IF(OPFILE) THEN
        CALL CLOSEF(IOFI,ERROR,*9999)
        IOFI=IOOP
      ENDIF
      RETURN 1
      END


