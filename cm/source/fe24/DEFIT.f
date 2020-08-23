      SUBROUTINE DEFIT(IBT,IDO,INP,LD,LGE,LN,NBH,NBJ,NBJF,
     '  NEELEM,NELIST,NFLIST,NFFACE,NHE,NHP,NKEF,NKH,NKHE,NKJE,
     '  NMNO,NNF,NPF,NP_INTERFACE,NPLIST,NPNE,NPNF,NPNODE,
     '  NPNY,NRLIST,NVHE,NVJE,NVJF,NVHP,NVJP,NXLIST,NYNE,NYNP,NYNR,
     '  SE,SF,WU,XA,XE,XID,XP,ZD,ZP,STRING,FIX,ERROR,*)

C#### Subroutine: DEFIT
C###  Description:
C###    DEFIT defines optimisation parameters for fitting geometry or
C###    field data points.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'data00.cmn'
      INCLUDE 'file00.cmn'
C LKC 24-OCT-97 unused block
      INCLUDE 'fit000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'ptr00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),LD(NDM),LGE(NHM*NSM,NRCM),LN(0:NEM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),NFLIST(0:NFM),
     '  NFFACE(0:NF_R_M,NRM),
     '  NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),NKEF(0:4,16,6,NBFM),
     '  NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NMNO(1:2,0:NOPM,NXM),NNF(0:17,6,NBFM),
     '  NPF(9,NFM),NP_INTERFACE(0:NPM,0:3),
     '  NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),
     '  NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM,NXM),
     '  NRLIST(0:NRM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVJE(NNM,NBFM,NJM,NEM),
     '  NVJF(NNM,NBFM,NJM),NVHP(NHM,NPM,NCM,0:NRM),NVJP(NJM,NPM),
     '  NXLIST(0:NXM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM)
      REAL*8 SE(NSM,NBFM,NEM),SF(NSM,NBFM),
     '  WU(0:NUM+1,NEM),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XID(NIM,NDM),XP(NKM,NVM,NJM,NPM),ZD(NJM,NDM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(*)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IBEG3,IBEG4,
     '  IEND,IEND1,IEND2,IEND3,IEND4,IPFILE,nf,
     '  N3CO,no_nrlist,nr,nx,nxc,NX_ACTION
      CHARACTER FILE*(MXCH),STATUS*3,DEFAULT1*30,DEFAULT2*30,
     '  DATA_TYPE*20
      LOGICAL ALL_REGIONS,CALCU,CBBREV,FILIO,GENER,MOUSE


      CALL ENTERS('DEFIT',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        IF(KTYP8.EQ.1) THEN
          DEFAULT1='geometry'
        ELSE IF(KTYP8.EQ.2) THEN
          DEFAULT1='fibre'
        ELSE IF(KTYP8.EQ.3) THEN
          DEFAULT1='field'
        ELSE IF(KTYP8.EQ.4) THEN
          DEFAULT1='signal'
        ELSE IF(KTYP8.EQ.5) THEN
          DEFAULT1='fourier'
        ELSE IF(KTYP8.EQ.6) THEN
          DEFAULT1='optimisation'
        ELSE IF(KTYP8.EQ.7) THEN
          DEFAULT1='material'
        ELSE
          DEFAULT1='field'
        ENDIF
        IF(KTYP6.EQ.1) THEN
          DEFAULT2='gauss'
        ELSE IF(KTYP6.EQ.2) THEN
          DEFAULT2='resid'
        ELSE IF(KTYP6.EQ.4) THEN
          DEFAULT2='grid'
        ELSE
          DEFAULT2='data'
        ENDIF
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)
        CALL STRING_TRIM(DEFAULT1,IBEG3,IEND3)
        CALL STRING_TRIM(DEFAULT2,IBEG4,IEND4)

C---------------------------------------------------------------------

C#### Command: FEM define fit;d/l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Define fit properties. Fit properties are read from or
C###    written to the file FILENAME (with extension .ipfit) in
C###    the directory specified by PATH. The fit can be of geometry
C###    data, fibre data, or motion data (default is geometry data).
C###    The fit can be linear or nonlinear, and smoothing constraints
C###    can be optionally specified.  This command allows the
C###    definition of which parameters and derivatives are to be
C###    included in the fit.  The fit takes the data points at their
C###    arbitrary points, and fits a field at the nodes which
C###    corresponds to that data information.
C###  Parameter:      <(fibre/field/geometry/material/patch/sheet/signal)[field]>
C###    Specify what is being fitted to.
C###  Parameter:      <(data/gauss/grid/residuals/node_group)[data]>
C###    Specify what information is to be fitted.
C###  Parameter:      <element (#s/calculated)[calculated]>
C###    Specify specific elements to be included in the fit or
C###    if data is being fitted, the elements containing data
C###    points can be calculated and those elements used in the
C###    fit.
C###  Parameter:      <faces (#s/calculated)[calculated]>
C###    Specify specific faces to be included in the fit or
C###    if data is being fitted, the faces containing data
C###    points can be calculated and those faces used in the
C###    fit.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions. Use only one region if using
C###    node_group.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of fit type) to use.
C###  Parameter:      <(lock/nolock)[lock]>
C###    Specify whether or not to lock this class number
C###    (of fit type).

        OP_STRING(1)=STRING(1:IEND)//';d/l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)
     '    //'<(fibre/field/geometry/material/patch/sheet/signal)['//
     '    DEFAULT1(IBEG3:IEND3)//']>'
        OP_STRING(3)=BLANK(1:15)//
     '    '<(data/gauss/grid/residuals/node_group)['
     '    //DEFAULT2(IBEG4:IEND4)//']>'
        OP_STRING(4)=BLANK(1:15)//
     '    '<element (#s/calculated)[calculated]>'
        OP_STRING(5)=BLANK(1:15)//
     '    '<faces (#s/calculated)[calculated]>'
        OP_STRING(6)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(7)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(8)=BLANK(1:15)//'<(lock/nolock)[lock]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define fit;d/l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Define fit properties. Fit properties are read from or
C###    written to the file FILENAME (with extension .ipfit) in
C###    the directory specified by PATH.
C###  Parameter:      optimisation
C###    Specify that you wish to use nonlinear geometric fitting. If
C###    this option is given the fit will be invoked by the
C###    optimise command instead of the fit command.
C###  Parameter:      <element (#s/calculated)[calculated]>
C###    Specify specific elements to be included in the fit or
C###    the elements containing data points can be calculated
C###    and those elements used in the fit.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of fit type) to use.
C###  Parameter:      <(lock/nolock)[lock]>
C###    Specify whether or not to lock this class number
C###    (of fit type).

        OP_STRING(1)=STRING(1:IEND)//';d/l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'optimisation'
        OP_STRING(3)=BLANK(1:15)//
     '    '<element (#s/calculated)[calculated]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(5)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(6)=BLANK(1:15)//'<(lock/nolock)[lock]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define fit;d/l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Define fit properties. Fit properties are read from or
C###    written to the file FILENAME (with extension .ipfit) in
C###    the directory specified by PATH.
C###  Parameter:      spline
C###    The spline option is used for signal fitting only. Splines
C###    are fitted in space for a specific time which is specified
C###    at the fit stage.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of fit type) to use.
C###  Parameter:      <(lock/nolock)[lock]>
C###    Specify whether or not to lock this class number
C###    (of fit type).

        OP_STRING(1)=STRING(1:IEND)//';d/l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'spline'
        OP_STRING(3)=BLANK(1:15)//
     '    '<element (#s/calculated)[calculated]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(5)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(6)=BLANK(1:15)//'<(lock/nolock)[lock]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define fit;d/l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]> tag
C###  Description:
C###    Define fit properties. Fit properties are read from or
C###    written to the file FILENAME (with extension .ipfit) in
C###    the directory specified by PATH. This invokes a tag fitting
C###    procedure where data weights represent normals to the tagging
C###    plane of the data point.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of fit type) to use.
C###  Parameter:      <(lock/nolock)[lock]>
C###    Specify whether or not to lock this class number
C###    (of fit type).

        OP_STRING(1)=STRING(1:IEND)//';d/l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']> tag'
        OP_STRING(3)=BLANK(1:15)//
     '    '<element (#s/calculated)[calculated]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(5)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(6)=BLANK(1:15)//'<(lock/nolock)[lock]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEFIT',ERROR,*9999)
      ELSE
        IPFILE=2 !is input file version number on 19-Dec-1994

        CALL PARSE_QUALIFIERS('DLPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)

        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
C KFA 17-DEC-2001
        TAG_FIT = .FALSE.
        IF(CBBREV(CO,'FIBRE',3,noco+1,NTCO,N3CO)) THEN
          CALL ASSERT(CALL_ELFB,'>>Define fibre elements first',
     '      ERROR,*9999)
          KTYP8=2
        ELSE IF(CBBREV(CO,'FIELD',3,noco+1,NTCO,N3CO)) THEN
          CALL ASSERT(CALL_ELFD,'>>Define field elements first',
     '      ERROR,*9999)
          KTYP8=3
        ELSE IF(CBBREV(CO,'FOURIER',3,noco+1,NTCO,N3CO)) THEN
          KTYP8=5
        ELSE IF(CBBREV(CO,'GEOMETRY',3,noco+1,NTCO,N3CO)) THEN
C LKC 24-OCT-97 Added assert below
          CALL ASSERT(CALL_DATA,'>>Define Data first',ERROR,*9999)
          CALL ASSERT(CALL_ELFD,'>>Define field elements first',
     '      ERROR,*9999)
          CALL ASSERT(CALC_XI,'>>Define Xi positions first',ERROR,*9999)
          KTYP8=1
C KFA 17-DEC-2001 Added new type it is the same as above except
C   the extra flag, which is declared in fit000.cmn
        ELSE IF(CBBREV(CO,'TAG',3,noco+1,NTCO,N3CO)) THEN
          CALL ASSERT(CALL_DATA,'>>Define Data first',ERROR,*9999)
          CALL ASSERT(CALL_ELFD,'>>Define field elements first',
     '      ERROR,*9999)
          CALL ASSERT(CALC_XI,'>>Define Xi positions first',ERROR,*9999)
          KTYP8=1
          TAG_FIT = .TRUE.
        ELSE IF(CBBREV(CO,'MATERIAL',3,noco+1,NTCO,N3CO)) THEN
          CALL ASSERT(CALL_ELFD,'>>Define field elements first',
     '      ERROR,*9999)
          KTYP8=7
        ELSE IF(CBBREV(CO,'OPTIMISATION',2,noco+1,NTCO,N3CO)) THEN
          CALL ASSERT(CALC_XI,'>>Define Xi positions first',ERROR,*9999)
          KTYP8=6
        ELSE IF(CBBREV(CO,'PATCH',3,noco+1,NTCO,N3CO)) THEN
          KTYP8=8
        ELSE IF(CBBREV(CO,'SIGNAL',2,noco+1,NTCO,N3CO)) THEN
          CALL ASSERT(CALL_ELFD,'>>Define field elements first',
     '      ERROR,*9999)
          KTYP8=4
        ELSE IF(CBBREV(CO,'SPLINE',2,noco+1,NTCO,N3CO)) THEN
          CALL ASSERT(CALL_ELFD,'>>Define field elements first',
     '      ERROR,*9999)
          KTYP8=9
        ELSE IF(CBBREV(CO,'SHEET',3,noco+1,NTCO,N3CO)) THEN
          CALL ASSERT(CALL_ELFB,'>>Define fibre elements first',ERROR,
     '      *9999)
          KTYP8=2
        ELSE
C MPN 17-Aug-95: set default only if current value is not set
          IF(KTYP8.LT.1.OR.KTYP8.GT.6) THEN
            CALL ASSERT(CALL_ELFD,'>>Define field elements first',
     '        ERROR,*9999)
            KTYP8=3
          ENDIF
        ENDIF

        IF(CBBREV(CO,'FACES',3,noco+1,NTCO,N3CO)) THEN
          CALL PARSE_FACES(NFFACE,NFLIST,noco-1,NRLIST,NTCO,CO,
     '      ERROR,*9999)
C PM 14Aug02: Put the faces into LN
C          NELIST(0)=0
          LN(0)=0
          DO nf=1,NFLIST(0)
C            NELIST(0)=NELIST(0)+1
C            NELIST(NELIST(0))=NPF(6,NFLIST(nf))
            LN(0)=LN(0)+1
            LN(LN(0))=NFLIST(nf)
          ENDDO !nf
          KTYP11=2
        ELSE
C          NFLIST(0)=0
C!!! Where is KTYP11 initialized???
          IF(KTYP11.NE.2) THEN  !not face fitting
            KTYP11=1
            NFLIST(0)=0
          ELSEIF(KTYP11.EQ.2) THEN ! face fitting and external faces
            LN(0)=0
C!!! Where is NFLIST initialized???
            DO nf=1,NFLIST(0)
              LN(0)=LN(0)+1
              LN(LN(0))=NFLIST(nf)
            ENDDO !nf
          ENDIF
        ENDIF

        IF(CBBREV(CO,'GAUSS',1,noco+1,NTCO,N3CO)) THEN
          DATA_TYPE='GAUSS'
          KTYP6=1
        ELSE IF(CBBREV(CO,'GRID',3,noco+1,NTCO,N3CO)) THEN
          DATA_TYPE='GRID'
          KTYP6=4
        ELSE IF(CBBREV(CO,'RESIDUALS',3,noco+1,NTCO,N3CO)) THEN
          DATA_TYPE='RESID'
          KTYP6=2
C TVK 14-04-2000 added 'node_group' data type for field fitting
        ELSE IF(CBBREV(CO,'NODE_GROUP',3,noco+1,NTCO,N3CO)) THEN
          DATA_TYPE='NODE_GROUP'
          DATATYPE_NG=DATA_TYPE
          KTYP6=0
        ELSE
C MPN 17-Aug-95: set default only if current value is not set
          IF(KTYP6.LT.0.OR.KTYP6.GT.1) KTYP6=0
          IF(KTYP6.NE.1) DATA_TYPE='DATA'
        ENDIF

C!!! MPN 8-Apr-96: LN is blocked by USE_DATA. Could perhaps use
C!!!               NELIST here instead if LN doesn't need to be kept
C!!!               around. Alternatively take out USE_DATA from memory
C!!!               allocation of LN in FEMA
        CALL ASSERT(USE_DATA.EQ.1,'>>USE_DATAMX must be 1',ERROR,*9999)
        IF(KTYP11.NE.2) THEN
          CALL PARSE_ELEMENTS(NEELEM,LN,noco,NRLIST,NTCO,CO,
     '      ERROR,*9999)
        ENDIF

        IF(CBBREV(CO,'NOLOCK',2,noco+1,NTCO,N3CO)) THEN
          NX_ACTION=NX_ALLOCATE
        ELSE
          NX_ACTION=NX_ALLOCATE_AND_LOCK
        ENDIF

        IF(KTYP8.EQ.6) THEN
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_OPTI,ERROR,*9999)
          IF(nx.EQ.0) THEN
            CALL NX_LOC(NX_ACTION,nxc,nx,NX_OPTI,ERROR,*9999)
          ENDIF
        ELSE
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_FIT,ERROR,*9999)
          IF(nx.EQ.0) THEN
            CALL NX_LOC(NX_ACTION,nxc,nx,NX_FIT,ERROR,*9999)
          ENDIF
       ENDIF

        IF(FILIO) THEN
          IF(iotype.NE.5) THEN
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'fit',STATUS,
     '        ERR,ERROR,*9999)
          ENDIF

          IF(ISC_GK_PTR(nx).EQ.0) THEN
            CALL ALLOCATE_MEMORY(NISC_GKM*USE_SPARSE,0,INTTYPE,
     '      ISC_GK_PTR(nx),MEM_INIT,ERROR,*9999)
          ENDIF
          IF(ISR_GK_PTR(nx).EQ.0) THEN
            CALL ALLOCATE_MEMORY(NISR_GKM*USE_SPARSE,0,INTTYPE,
     '      ISR_GK_PTR(nx),MEM_INIT,ERROR,*9999)
          ENDIF

c cpb 20/3/95 Need to put in region prompting
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            IF(KTYP8.EQ.1.OR.KTYP8.EQ.3.OR.KTYP8.EQ.4.OR.
     '        KTYP8.EQ.7.OR.KTYP8.EQ.8.) THEN
              CALL ASSERT(NJ_LOC(NJL_FIEL,0,nr).GT.0,
     '          '>>Define field first',ERROR,*9999)
            ELSE IF(KTYP8.EQ.2) THEN
              CALL ASSERT(NJ_LOC(NJL_FIBR,0,nr).GT.0,
     '          '>>Define field first',ERROR,*9999)
            ENDIF
            CALL IPFIT(IBT,IDO,INP,%VAL(ISC_GK_PTR(nx)),
     '        %VAL(ISR_GK_PTR(nx)),LD,LGE,LN,NBH,NBJ,NBJF,
     '        NEELEM,NELIST,NFFACE,NFLIST,
     '        NHE(1,nx),NHP(1,0,nx),NKEF,NKH,NKHE,NKJE,
     '        NMNO(1,0,nx),NNF,NPF,
     '        NP_INTERFACE,NPLIST,NPNE,NPNF,NPNODE,NPNY(0,1,0,nx),
     '        nr,NRLIST,NVHE,NVJE,NVJF,NVHP,
     '        NVJP,nx,NYNE,NYNP,NYNR(0,0,1,0,nx),SE,SF,WU,
     '        XA,XE,XID,XP,ZD,ZP,DATA_TYPE,FIX(1,1,nx),ERROR,*9999)
            ENDDO !no_nrlist
          CALL CLOSEF(IFILE,ERROR,*9999)

        ENDIF !filio

        CALL_FIT=.TRUE.
      ENDIF

      CALL EXITS('DEFIT')
      RETURN
 9999 CALL ERRORS('DEFIT',ERROR)
      IF(KTYP8.EQ.6) THEN
        IF(nx.GT.0) THEN
          CALL NX_LOC(NX_FREE,nxc,nx,NX_OPTI,ERROR,*9998)
        ENDIF
      ELSE
        IF(nx.GT.0) THEN
          CALL NX_LOC(NX_FREE,nxc,nx,NX_FIT,ERROR,*9998)
        ENDIF
      ENDIF
 9998 CALL EXITS('DEFIT')
      RETURN 1
      END


