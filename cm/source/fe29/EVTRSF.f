      SUBROUTINE EVTRSF(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,
     '  IBT,IDO,INP,IPIV,ISC_GKK,ISIZE_TBH,
     '  ISR_GKK,LD,LGE,NAN,NBH,NBJ,NDET,NDIPOLES,
     '  NEELEM,NELIST,NENP,NGAP,NHE,NHP,NKB,NKH,NKHE,NKJE,NLL,NNB,NNF,
     '  NNL,NONY,NORD,NPF,NP_INTERFACE,NPL,NPLIST3,NPLIST4,
     '  NPNE,NPNODE,NPNY,NRE,NVHE,NVHP,NVJE,NW,NWP,NXI,NXLIST,
     '  NYNE,NYNO,NYNP,NYNR,NYNY,NYQNR,CE,CGE,CONY,CP,CURVCORRECT,
     '  CYNO,CYNY,DET,DIPOLE_CEN,DIPOLE_DIR,DL,GD,
     '  GKK,GR,GRR,PG,SE,SP,T_BH,WG,WK1_INV,WK2_INV,
     '  WK3_INV,WK4_INV,XA,XE,XG,XID,XIG,XO,XP,
     '  YG,YP,ZA,ZE,ZP,STRING,FIX,ERROR,*)

C#### Subroutine: EVTRSF
C###  Description:
C###    EVTRSF is a buffer routine for EVTRSF_DYNAM

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'ptr00.cmn'

!     Parameter List
      INTEGER DIPOLE_CEN_NTIME(NDIPOLEM,NRM,NXM),
     '  DIPOLE_DIR_NTIME(NDIPOLEM,NRM,NXM),IBT(3,NIM,NBFM),
     '  IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),IPIV(NY_TRANSFER_M),
     '  ISC_GKK(NISC_GKKM,NXM),
     '  ISIZE_TBH(2),
     '  ISR_GKK(NISR_GKKM,NXM),LD(NDM),LGE(NHM*NSM,NRCM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NDET(NBFM,0:NNM),NDIPOLES(NRM,NXM),NEELEM(0:NE_R_M,0:NRM),
     '  NELIST(0:NEM),
     '  NENP(NPM,0:NEPM,0:NRM),NGAP(NIM,NBM),NHE(NEM,NXM),
     '  NHP(NPM,0:NRM,NXM),NKB(2,2,2,NNM,NBFM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),
     '  NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),
     '  NONY(0:NOYM,NYM,NRCM,0:NRM,NXM),NORD(5,NE_R_M),
     '  NPF(9,NFM),NP_INTERFACE(0:NPM,0:3),
     '  NPL(5,0:3,NLM),NPLIST3(0:NPM),NPLIST4(0:NPM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),NRE(NEM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM),NWP(NPM,2),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),NXLIST(0:NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),
     '  NYNY(0:NYYM,NYM,NRM,NXM),
     '  NYQNR(0:NYQM,0:NRCM,NCM,0:NRM,NXM)
      REAL*8 CE(NMM,NEM,NXM),CGE(NMM,NGM,NEM,NXM),
     '  CONY(0:NOYM,NYM,NRCM,0:NRM,NXM),
     '  CP(NMM,NPM,NXM),
     '  CURVCORRECT(2,2,NNM,NEM),CYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),
     '  CYNY(0:NYYM,NYM,NRM,NXM),
     '  DET(NBFM,0:NNM,NGM,6),DIPOLE_CEN(4,0:NDIPTIMM,NDIPOLEM,NRM,NXM),
     '  DIPOLE_DIR(4,0:NDIPTIMM,NDIPOLEM,NRM,NXM),
     '  DL(3,NLM),GD(NZ_GD_M),GKK(NZ_GKK_M,NXM),GR(NYROWM),
     '  GRR(NOM),PG(NSM,NUM,NGM,NBM),SE(NSM,NBFM,NEM),SP(NKM,NBFM,NPM),
     '  T_BH(NY_TRANSFER_M,NY_TRANSFER_M),
     '  WG(NGM,NBM),WK1_INV(NY_TRANSFER_M,NY_TRANSFER_M),
     '  WK2_INV(NY_TRANSFER_M,NY_TRANSFER_M),
     '  WK3_INV(NY_TRANSFER_M,NY_TRANSFER_M),WK4_INV(NY_TRANSFER_M),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XG(NJM,NUM),XID(NIM,NDM),XIG(NIM,NGM,NBM),
     '  XO(NOM,NXM),XP(NKM,NVM,NJM,NPM),
     '  YG(NIYGM,NGM,NEM),YP(NYM,NIYM,NXM),
     '  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER STRING*(MXCH),ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER IBEG,IEND,N3CO,nx,nxc
      CHARACTER FILE*(MXCH)
      LOGICAL AT_NODES,CBBREV,CHECK,OPFILE,SVD

      CALL ENTERS('EVTRSF',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM evaluate transfer<;FILENAME>
C###  Parameter:        <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Parameter:        <check>
C###    Provides additional output. Row sums of T_BH
C###    corresponding to potential should sum to 0 since
C###     a closed double layer should generate no external potential.
C###  Parameter:        <svd>
C###    Evaluates a singular value decomposition of T_BH.
C###    Finds the singular values and condition number of T_BH.
C###  Parameter:        <at (nodes/electrodes)[nodes]>
C###    By default, signals are at nodes, but may be specified to be
C###    at electrode (data) locations.
C###  Description:
C###    Evaluates a transfer matrix that maps the solution from one
C###    surface to those on another. Assumes the second surface is a
C###    no flux surface.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<check>'
        OP_STRING(3)=BLANK(1:15)//'<svd>'
        OP_STRING(4)=BLANK(1:15)//'<at (nodes/electrodes)[nodes]>'
        OP_STRING(5)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','EVTRSF_DYNAM',ERROR,*9999)
      ELSE
        CALL ASSERT(CALL_TRANSFER,
     '    '>>Define the transfer first',ERROR,*9999)
        CALL ASSERT(USE_TRANSFER.EQ.1,
     '    '>>Set USE_TRANSFER to 1 in parameters file',ERROR,*9999)
        IF(NTCOQU(noco).GT.0) THEN
          FILE=COQU(noco,1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          OPFILE=.TRUE.
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.optrsf','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
          IOFI=IOOP
        ENDIF

        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)
        CALL ASSERT(CALL_SOLV,'>>Need to define solve first',
     '    ERROR,*9999)

        IF(CBBREV(CO,'CHECK',2,noco+1,NTCO,N3CO)) THEN
          CHECK=.TRUE.
        ELSE
          CHECK=.FALSE.
        ENDIF

        IF(CBBREV(CO,'SVD',2,noco+1,NTCO,N3CO)) THEN
          SVD=.TRUE.
        ELSE
          SVD=.FALSE.
        ENDIF

        IF(CBBREV(CO,'AT',2,noco+1,NTCO,N3CO)) THEN
          IF(CBBREV(CO,'ELECTRODES',1,n3co+1,n3co+2,N3CO)) THEN
            CALL ASSERT(CALL_XI,'>>Calculate xi first',ERROR,*9999)
            AT_NODES=.FALSE.
          ELSE
            AT_NODES=.TRUE.
          ENDIF
        ELSE
          AT_NODES=.TRUE.
        ENDIF

        IF(ISC_GK_PTR(nx).EQ.0) THEN
          CALL ALLOCATE_MEMORY(NISC_GKM*USE_SPARSE,0,INTTYPE,
     '      ISC_GK_PTR(nx),MEM_INIT,ERROR,*9999)
        ENDIF
        IF(ISR_GK_PTR(nx).EQ.0) THEN
          CALL ALLOCATE_MEMORY(NISR_GKM*USE_SPARSE,0,INTTYPE,
     '      ISR_GK_PTR(nx),MEM_INIT,ERROR,*9999)
        ENDIF
        IF(GK_PTR(nx).EQ.0) THEN
          CALL ALLOCATE_MEMORY(NZ_GK_M,1,DPTYPE,GK_PTR(nx),MEM_INIT,
     '      ERROR,*9999)
        ENDIF
        IF(ISC_GQ_PTR(nx).EQ.0) THEN
          CALL ALLOCATE_MEMORY(NISC_GQM*USE_SPARSE,0,INTTYPE,
     '      ISC_GQ_PTR(nx),MEM_INIT,ERROR,*9999)
        ENDIF
        IF(ISR_GQ_PTR(nx).EQ.0) THEN
          CALL ALLOCATE_MEMORY(NISR_GQM*USE_SPARSE,0,INTTYPE,
     '      ISR_GQ_PTR(nx),MEM_INIT,ERROR,*9999)
        ENDIF
        IF(GQ_PTR(nx).EQ.0) THEN
          CALL ALLOCATE_MEMORY(NZ_GQ_M,1,DPTYPE,GQ_PTR(nx),MEM_INIT,
     '      ERROR,*9999)
        ENDIF

        CALL EVTRSF_DYNAM(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,
     '    IBT,IDO,INP,IPIV,%VAL(ISC_GK_PTR(nx)),ISC_GKK,
     '    %VAL(ISC_GQ_PTR(nx)),
     '    ISIZE_TBH,%VAL(ISR_GK_PTR(nx)),
     '    ISR_GKK,%VAL(ISR_GQ_PTR(nx)),LD,LGE,NAN,NBH,
     '    NBJ,NDET,NDIPOLES,
     '    NEELEM,NELIST,NENP,NGAP,NHE,NHP,NKB,NKH,NKHE,NKJE,NLL,NNB,NNF,
     '    NNL,NONY,NORD,NPF,NP_INTERFACE,NPL,NPLIST3,NPLIST4,
     '    NPNE,NPNODE,NPNY,NRE,NVHE,NVHP,NVJE,NW(1,1,nx),NWP,nx,NXI,
     '    NYNE,NYNO,NYNP,NYNR,NYNY,NYQNR,CE,CGE,CONY,CP,CURVCORRECT,
     '    CYNO,CYNY,DET,DIPOLE_CEN,DIPOLE_DIR,DL,GD,
     '    %VAL(GK_PTR(nx)),
     '    GKK,%VAL(GQ_PTR(nx)),GR,GRR,PG,SE,SP,
     '    T_BH,WG,WK1_INV,WK2_INV,WK3_INV,WK4_INV,
     '    XA,XE,XG,XID,XIG,XO,XP,
     '    YG,YP,ZA,ZE,ZP,AT_NODES,CHECK,FIX,OPFILE,SVD,ERROR,*9999)

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
        EVALUATE_TRANSFER=.TRUE.
      ENDIF

      CALL EXITS('EVTRSF')
      RETURN
 9999 CALL ERRORS('EVTRSF',ERROR)
      CALL EXITS('EVTRSF')
      RETURN 1
      END


