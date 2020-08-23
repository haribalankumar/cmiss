      SUBROUTINE FIT(IBT,IDO,INP,IPIVOT,ISC_GKK,
     '  ISR_GKK,LD,LGE,LN,NAN,NBH,NBHF,NBJ,NBJF,NDDATA,NDDL,NDLT,
     '  NEELEM,NELIST,NENP,NENQ,NEP,NFF,NFFACE,NGAP,NHE,NHP,NHQ,
     '  NKB,NKEF,NKH,NKHE,NKJE,NLL,
     '  NMNO,NNB,NNF,NNL,NONY,NPF,
     '  NP_INTERFACE,NPL,NPLIST,NPLIST2,NPNE,NPNF,NPNODE,
     '  NPNY,NQNE,NQNY,NQS,NQXI,
     '  NRE,NRLIST,NRLIST2,NSB,NVHE,NVHP,NVJE,NVJF,NW,
     '  NWP,NXI,NXLIST,NYNE,NYNO,NYNP,NYNR,NYNY,NYQNR,Z_CONT_LIST,
     '  CE,CG,CGE,CONY,CYNY,CP,CURVCORRECT,CYNO,EDD,
     '  ER,ES,FEXT,GKK,GR,GRR,PG,RE1,RG,SE,SF,SP,WD,WDL,WG,
     '  WU,XA,XE,
     '  XG,XID,XIDL,XIG,XIP,XO,XP,YG,YGF,YP,YQ,YQS,ZA,ZA1,Z_CONT,ZD,
     '  ZDL,ZE,ZP,ZP1,STRING,FIX,ERROR,*)

C#### Subroutine: FIT
C###  Description:
C###    FIT fits geometry or field parameters to data points defined
C###    in ZD, YG, or YQS
C###    or fits material parameters to minimise eqtn residuals.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'fit001.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'ptr00.cmn'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),IPIVOT(NOM),
     '  ISC_GKK(NISC_GKKM,NXM),
     '  ISR_GKK(NISR_GKKM,NXM),
     '  LD(NDM),LGE(NHM*NSM,NRCM),LN(0:NEM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),
     '  NBJ(NJM,NEM),NBJF(NJM,NFM),NDDATA(0:NDM,0:NRM),NDDL(NEM,NDEM),
     '  NDLT(NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NENP(NPM,0:NEPM,0:NRM),NENQ(0:8,NQM),NEP(NPM),
     '  NFF(6,NEM),NFFACE(0:NF_R_M,NRM),
     '  NGAP(NIM,NBM),NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),
     '  NHQ(NRM,NXM),NKB(2,2,2,NNM,NBFM),
     '  NKEF(0:4,16,6,NBFM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),
     '  NMNO(1:2,0:NOPM,NXM),NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),
     '  NNL(0:4,12,NBFM),NONY(0:NOYM,NYM,NRCM,0:NRM,NXM),NPF(9,NFM),
     '  NP_INTERFACE(0:NPM,0:3),NPL(5,0:3,NLM),
     '  NPLIST(0:NPM),NPLIST2(0:NPM),
     '  NPNF(NNM,NBFM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM,NXM),
     '  NQNE(NEQM,NQEM),NQNY(2,NYQM,0:NRCM,NXM),NQS(NEQM),
     '  NQXI(0:NIM,NQSCM),NRE(NEM),NRLIST(0:NRM),NRLIST2(0:NRM),
     '  NSB(NKM,NNM,NBFM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NVJE(NNM,NBFM,NJM,NEM),
     '  NVJF(NNM,NBFM,NJM),NWP(NPM,2),NW(NEM,3,NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),NXLIST(0:NXM),
     '  NYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),
     '  NYNY(0:NYYM,NYM,NRM,NXM),
     '  NYQNR(0:NYQM,0:NRCM,NCM,0:NRM,NXM),
     '  Z_CONT_LIST(NDM,2,7)
      REAL*8 CE(NMM,NEM,NXM),CG(NMM,NGM),CGE(NMM,NGM,NEM,NXM),
     '  CONY(0:NOYM,NYM,NRCM,0:NRM,NXM),
     '  CP(NMM,NPM,NXM),CURVCORRECT(2,2,NNM,NEM),
     '  CYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),CYNY(0:NYYM,NYM,NRM,NXM),
     '  EDD(NDM),ER(NHM*NSM),ES(NHM*NSM,NHM*NSM),FEXT(NIFEXTM,NGM,NEM),
     '  GKK(NZ_GKK_M,NXM),GR(NYROWM),
     '  GRR(NOM),PG(NSM,NUM,NGM,NBM),RE1(NSM,NHM),RG(NGM),
     '  SE(NSM,NBFM,NEM),SF(NSM,NBFM),
     '  SP(NKM,NBFM,NPM),WD(NJM,NDM),WDL(NHM,NDEM),
     '  WG(NGM,NBM),
     '  WU(0:NUM+1,NEM),XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),
     '  XID(NIM,NDM), XIDL(NIM,NDEM),XIG(NIM,NGM,NBM),XIP(NIM,NPM),
     '  XO(NOM,NXM),XP(NKM,NVM,NJM,NPM),
     '  YG(NIYGM,NGM,NEM),YGF(NIYGFM,NGFM,NFM),YP(NYM,NIYM,NXM),
     '  YQ(NYQM,NIQM,NAM,NXM),YQS(NIQSM,NQM),ZA(NAM,NHM,NCM,NEM),
     '  ZA1(NAM,NHM,NCM,NEM),Z_CONT(NDM,2,67),
     '  ZD(NJM,NDM),ZDL(NHM,NDEM),
     '  ZE(NSM,NHM),
     '  ZP(NKM,NVM,NHM,NPM,NCM),ZP1(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER IBEG,IEND,IBEG1,IEND1,IFROMC,N3CO,nelec,no_nrlist,
     '  nr,nslice,nxFIT,nxSOLVE,nxcFIT,nxcSOLVE,START_NODE,SNODE_NUM
      REAL*8 RFROMC,TSTART,TEND,ZE1(NSM,NHM)
      CHARACTER FILEFORMAT*6
      LOGICAL ALL_REGIONS,CBBREV,IN_PLANE,SPLINE

      CALL ENTERS('FIT',*9999)

      write(*,*) 'entering femfit'

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM fit
C###  Parameter:   <(calculate_gsm/previous_gsm)[calculate_gsm]>
C###    Calculate the global stiffness matrix or reuse a
C###    previous one.
C###  Parameter:   <in_plane>
C###    Takes only into account the error between the target point
C###    and the point in the element that lies in a plane given by the
C###    normal-vector of the plane specified as the three weights in 
C###    the ipdata file. The error component along the normal-vector
C###    will be ignored.
C###  Parameter:   <region (#s/all)[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not
C###    specified will be skipped.
C###  Parameter:   <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Parameter:   <solve_class #[1]>
C###  Description:
C###    Fit nodal parameters to data using least squares.
C###    Geometry: fit nodal coordinates to give best fit of line or
C###    surface to data points.
C###    Fibre: fit nodal fibre parameters to give best fit to fibre
C###    data defined over a surface or volume.
C###    Field: fit nodal field parameters to give best fit to field
C###    data defined over 2D or 3D space.
C###    Fourier:.
C###    Gauss:.
C###    In all cases the Xi-coordinate positions of the data points
C###    must have been calculated with >fem de data;c xi (or read in)
C###    and the fitting parameters defined with >fem de fit etc.
C###    After a geometric fit the nodal coordinates of the mesh are
C###    updated by...

        OP_STRING(1)=STRING(1:IEND)
     '    //' <(calculate_gsm/previous_gsm)[calculate_gsm]>'
        OP_STRING(2)=BLANK(1:15)//'<in_plane>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(5)=BLANK(1:15)//'<solve_class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM fit
C###  Parameter:   <signal>
C###    Indicates a field will be fitted from signal data.
C###  Parameter:   <tstart #[beginning]>
C###    Sets the start time of the input signal to fit from
C###  Parameter:   <tend #[end]>
C###    Sets the end time of the input signal to fit to.
C###  Parameter:   <(ascii/binary)[ascii]>
C###    Specify whether the signal file is stored as a binary or
C###    ascii file (ipsign or binsig).
C###  Parameter:   <region #[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped.
C###  Parameter:   <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description: Fits a signal fit to create a history file.
C###

        OP_STRING(1)=STRING(1:IEND)//' <signal>'
        OP_STRING(2)=BLANK(1:15)//'<tstart #[beginning]>'
        OP_STRING(3)=BLANK(1:15)//'<tend #[end]>'
        OP_STRING(4)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(5)=BLANK(1:15)//'<region #[1]>'
        OP_STRING(6)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM fit
C###  Parameter:   <signal> <spline>
C###    Specifies that a signal is being fitted using a cubic spline
C###  Parameter:   <electrodes_per_row> #[0]
C###    The number of electrodes in each row of the grid
C###  Parameter:   <rows> #[0]
C###    Specify the number of rows of electrodes
C###  Parameter:   <offset_fit> #[0]
C###    Specify the offset of the nodes to fit to.
C###  Parameter:   <mesh_offset> #[0]
C###    Specify the offset of the mesh fitted from.
C###  Parameter:   <from (data/nodes)[data]>
C###    Specifies whether the fit is from data or nodes
C###  Parameter:   <tstart #[beginning]>
C###    Specifies the limits of the fit interval
C###  Parameter:   <tend #[end]>
C###    Specifies the limits of the fit interval
C###  Parameter:   <(ascii/binary)[ascii]>
C###    Specify whether the file is stored as a binary or
C###    ascii file (ipsign or binsig).
C###  Parameter:   <region #[1]>]
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped.
C###  Parameter:   <class #[1]>
C###    Specify the class number (of solve type) to use.

C###  Description: Fits a signal using cubic splines to create a
C###    history file
C###
        OP_STRING(1)=STRING(1:IEND)//' <signal> <spline>'
        OP_STRING(2)=BLANK(1:15)//'<electrodes_per_row> #[0]'
        OP_STRING(3)=BLANK(1:15)//'<rows> #[0]'
        OP_STRING(4)=BLANK(1:15)//'<offset_fit> #[0]'
        OP_STRING(5)=BLANK(1:15)//'<mesh_offset> #[0]'
        OP_STRING(6)=BLANK(1:15)//'<from (data/nodes)[data]>'
        OP_STRING(7)=BLANK(1:15)//'<tstart #[beginning]>'
        OP_STRING(8)=BLANK(1:15)//'<tend #[end]>'
        OP_STRING(9)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(10)=BLANK(1:15)//'<region #[1]>'
        OP_STRING(11)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
C#### Command: FEM fit
C###  Parameter:   <patch>
C###    Specify that patch fitting is the method to use.
C###  Parameter:   <nodes (#s/all)[all]>
C###    Specify the nodes for which fitted values are to be found.
C###  Parameter:   <region #[1]>
C###    Specify the region number.
C###  Parameter:   <class #[1]>
C###    Specify the class number (of solve type) to fit.
C###  Description:
C###    Fits a polynomial to Gauss point data in an element
C###    patch surrounding a node. The nodal value is then recovered
C###    by evaluating the polynomial at the nodes coordinates.
C###    Fields for the recovered values must already be set up
C###    using >FEM define field and >FEM define elements field
C###    and the fit variables set with >FEM define fit gauss patch
C###  See-Also: FIT_PATCH

        OP_STRING(1)=STRING(1:IEND)//' <patch>'
        OP_STRING(2)=BLANK(1:15)//'<nodes (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<region #[1]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','FIT',ERROR,*9999)
      ELSE
        CALL ASSERT(KTYP8.GT.0,'>>Define fit first',ERROR,*9999)

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_NODES(NPNODE,NPLIST,noco,NRLIST,NTCO,CO,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxcFIT=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxcFIT,nxFIT,NX_FIT,ERROR,*9999)
        CALL ASSERT(nxFIT.GT.0,
     '    '>>no nx defined for this fit class',ERROR,*9999)

        IF(ISC_GK_PTR(nxFIT).EQ.0) THEN
          CALL ALLOCATE_MEMORY(NISC_GKM*USE_SPARSE,0,INTTYPE,
     '      ISC_GK_PTR(nxFIT),MEM_INIT,ERROR,*9999)
        ENDIF
        IF(ISR_GK_PTR(nxFIT).EQ.0) THEN
          CALL ALLOCATE_MEMORY(NISR_GKM*USE_SPARSE,0,INTTYPE,
     '      ISR_GK_PTR(nxFIT),MEM_INIT,ERROR,*9999)
        ENDIF
        IF(GK_PTR(nxFIT).EQ.0) THEN
          CALL ALLOCATE_MEMORY(NZ_GK_M,1,DPTYPE,GK_PTR(nxFIT),MEM_INIT,
     '      ERROR,*9999)
        ENDIF
        IF(ISC_GQ_PTR(nxFIT).EQ.0) THEN
          CALL ALLOCATE_MEMORY(NISC_GQM*USE_SPARSE,0,INTTYPE,
     '      ISC_GQ_PTR(nxFIT),MEM_INIT,ERROR,*9999)
        ENDIF
        IF(ISR_GQ_PTR(nxFIT).EQ.0) THEN
          CALL ALLOCATE_MEMORY(NISR_GQM*USE_SPARSE,0,INTTYPE,
     '      ISR_GQ_PTR(nxFIT),MEM_INIT,ERROR,*9999)
        ENDIF
        IF(GQ_PTR(nxFIT).EQ.0) THEN
          CALL ALLOCATE_MEMORY(NZ_GQ_M,1,DPTYPE,GQ_PTR(nxFIT),MEM_INIT,
     '      ERROR,*9999)
        ENDIF

        IF(KTYP8.EQ.7) THEN !material parameter fitting
          IF(CBBREV(CO,'solve_class',7,noco+1,NTCO,n3co)) THEN
            nxcSOLVE=IFROMC(CO(n3co+1))
          ELSE
            nxcSOLVE=1
          ENDIF
          CALL NX_LOC(NX_INQUIRE,nxcSOLVE,nxSOLVE,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nxSOLVE.GT.0,
     '      '>>no nx defined for this solve class',ERROR,*9999)
        ELSE                !all other fitting
          nxSOLVE=1
        ENDIF !ktyp8

C*** Input for in_plane geometric fitting

        IF(CBBREV(CO,'IN_PLANE',3,noco+1,NTCO,N3CO)) THEN
          IN_PLANE=.TRUE.
        ELSE
          IN_PLANE=.FALSE.
        ENDIF

C*** Input for SPLINE SIGNAL FITTING

        IF(CBBREV(CO,'SIGNAL',2,noco+1,NTCO,N3CO)) THEN
          IF(CBBREV(CO,'SPLINE',2,noco+1,NTCO,N3CO)) THEN
             CALL ASSERT(KTYP8.EQ.9,'>>Define spline signal fit first',
     '        ERROR,*9999)
            SPLINE=.TRUE.

C nelec=ELECTRODES_PER_ROW
            IF(CBBREV(CO,'ELECTRODES_PER_ROW',3,noco+1,NTCO,N3CO)) THEN
              CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
              nelec=IFROMC(CO(N3CO+1)(IBEG1:IEND1))
            ELSE
              ERROR='Missing ELECTRODES_PER_ROW'
              GOTO 9999
            ENDIF

C ROWS=NELEC
            IF(CBBREV(CO,'ROWS',3,noco+1,NTCO,N3CO)) THEN
              CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
              nslice=IFROMC(CO(N3CO+1)(IBEG1:IEND1))
            ELSE
              ERROR='Missing ROW parameter'
              GOTO 9999
            ENDIF

C OFFSET_FIT=NSLICE
            IF(CBBREV(CO,'OFFSET_FIT',3,noco+1,NTCO,N3CO)) THEN
              CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
              START_NODE=IFROMC(CO(N3CO+1)(IBEG1:IEND1))
            ELSE
              ERROR='Missing OFFSET_FIT parameter'
              GOTO 9999
            ENDIF

C MESH_OFFSET=START NODENUM
            IF(CBBREV(CO,'MESH_OFFSET',3,noco+1,NTCO,N3CO)) THEN
              CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
              SNODE_NUM=IFROMC(CO(N3CO+1)(IBEG1:IEND1))
            ELSE
              ERROR='Missing MESH_OFFSET parameter'
              GOTO 9999
            ENDIF

          ELSE
            CALL ASSERT(KTYP8.EQ.4,'>>Define signal fit first',
     '        ERROR,*9999)
            START_NODE=0
            SPLINE=.FALSE.
          ENDIF
        ENDIF

        IF(CBBREV(CO,'PATCH',2,noco+1,NTCO,N3CO)) THEN
          CALL ASSERT(KTYP8.EQ.8,
     '      '>>define fit patch first',ERROR,*9999)
        ENDIF

        IF(KTYP6.NE.1) THEN
          IF(KTYP8.EQ.6) THEN
            ERROR='>>Use optimise to fit data by optimisation'
            GOTO 9999
          ELSE
            IF((KTYP8.LT.1).OR.(KTYP8.GT.9)
     '        .OR.(KTYP8.EQ.8)) THEN
              ERROR='>>Fit is not defined'
              GOTO 9999
            ENDIF
          ENDIF
        ENDIF !ktyp26.ne.1

        IF(KTYP8.EQ.4.OR.KTYP8.EQ.9) THEN
          IF(CBBREV(CO,'TSTART',2,noco+1,NTCO,N3CO)) THEN
            TSTART=RFROMC(CO(N3CO+1))
          ELSE
            TSTART=-RMAX
          ENDIF

          IF(CBBREV(CO,'TEND',2,noco+1,NTCO,N3CO)) THEN
            TEND=RFROMC(CO(N3CO+1))
          ELSE
            TEND=RMAX
          ENDIF

          IF(CBBREV(CO,'BINARY',2,noco+1,NTCO,N3CO)) THEN
            FILEFORMAT='BINARY'
          ELSE
            FILEFORMAT='ASCII'
          ENDIF
        ENDIF !ktyp8=4 or 9

        FIRSTA=.TRUE.  !for full LU decomposition first time round
        IOTYPE=1

        IF(KTYP8.EQ.4.OR.(KTYP8.EQ.9)) THEN
          nr=NRLIST(1)

          IF(.NOT.SPLINE) THEN
            CALL FITSIG(IBT,IDO,INP,%VAL(ISC_GK_PTR(nxFIT)),
     '        ISC_GKK(1,nxFIT),
     '        %VAL(ISC_GQ_PTR(nxFIT)),%VAL(ISR_GK_PTR(nxFIT)),
     '        ISR_GKK(1,nxFIT),
     '        %VAL(ISR_GQ_PTR(nxFIT)),LD,LGE,LN,NBH,NBJ,
     '        NDDATA,NDDL,NDLT,NEELEM,NENP,NHE(1,nxFIT),NHP(1,nr,nxFIT),
     '        NHQ(1,nxFIT),NKB,NKH(1,1,1,nr),NKHE,NKJE,NLL,NNB,NNF,NNL,
     '        NONY(0,1,1,0,nxFIT),NPF,NP_INTERFACE,NPL,NPNE,NPNODE,
     '        NPNY(0,1,0,nxFIT),NQNY(1,1,0,nxFIT),nr,NRE,NRLIST,NRLIST2,
     '        NVHE,NVJE,NVHP(1,1,1,nr),
     '        NW(1,1,nxFIT),NWP,nxFIT,NXI,NYNE,
     '        NYNO(0,1,1,0,nxFIT),NYNP,
     '        NYNR(0,0,1,0,nxFIT),NYNY(0,1,1,nxFIT),
     '        NYQNR(0,0,1,0,nxFIT),CONY(0,1,1,0,nxFIT),CURVCORRECT,
     '        CYNO(0,1,1,0,nxFIT),CYNY(0,1,1,nxFIT),EDD,ER,ES,
     '        %VAL(GK_PTR(nxFIT)),GKK(1,nxFIT),GR,GRR,
     '        %VAL(GQ_PTR(nxFIT)),PG,RG,SE,SF,SP,TEND,TSTART,
     '        WD,WDL,WG,WU,XA,XE,XG,XID,XIDL,XO(1,nxFIT),XP,
     '        YP(1,1,nxFIT),YQ(1,1,1,nxFIT),YQS,ZA,ZD,ZDL,ZE,ZP,
     '        FILEFORMAT,FIRSTA,FIX(1,1,nxFIT),IN_PLANE,ERROR,*9999)
           ELSE

             CALL FITSIG_SPLINE(IBT,IDO,INP,LD,LN,NBH,NBJ,NDDATA,
     '         NEELEM,nelec,NHQ(1,nxFIT),NKJE,NPF,NPNE,NPNODE,
     '         NPNY(0,1,0,nxFIT),NQNY(1,1,0,nxFIT),nslice,nr,NRE,
     '         NRLIST,NRLIST2,NVJE,nxFIT,NYNR(0,0,1,0,nxFIT),
     '         NYQNR(0,0,1,0,nxFIT),SE,START_NODE,SNODE_NUM,TEND,
     '         TSTART,WD,XID,XE,XP,YP(1,1,nxFIT),YQ(1,1,1,nxFIT),YQS,
     '         ZD,FILEFORMAT,ERROR,*9999)
           ENDIF

        ELSE !all other ktyp8s
          IF(KTYP8.NE.8) THEN
            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              CALL FITFLD(IBT,IDO,INP,
     '          %VAL(ISC_GK_PTR(nxFIT)),ISC_GKK(1,nxFIT),
     '          %VAL(ISC_GQ_PTR(nxFIT)),%VAL(ISR_GK_PTR(nxFIT)),
     '          ISR_GKK(1,nxFIT),
     '          %VAL(ISR_GQ_PTR(nxFIT)),LGE,LN,
     '          NAN,NBH,NBHF,NBJ,NBJF,NDDL,NDLT,NEELEM,NENP,NENQ,
     '          NEP,NFF,NFFACE,NGAP,NHE(1,nxSOLVE),
     '          NHP(1,nr,nxSOLVE),NKB,NKEF,NKH,NKHE,NKJE,NLL,
     '          NMNO(1,0,nxFIT),NNB,NNF,NNL,NONY(0,1,1,0,nxFIT),
     '          NPF,NP_INTERFACE,NPL,NPLIST,NPLIST2,
     '          NPNE,NPNF,NPNODE,NPNY,NQNE,NQS,NQXI,nr,NRE,NSB,NVHE,
     '          NVHP(1,1,1,nr),NVJE,NVJF,NW(1,1,nxFIT),NWP,
     '          nxFIT,nxSOLVE,NXI,NYNE,NYNO(0,1,1,0,nxFIT),
     '          NYNP,NYNR,NYNY(0,1,1,nxFIT),Z_CONT_LIST,
     '          CE,CG,CGE,CONY(0,1,1,0,nxFIT),
     '          CP,CURVCORRECT,
     '          CYNO(0,1,1,0,nxFIT),CYNY(0,1,1,nxFIT),EDD,ER,ES,FEXT,
     '          %VAL(GK_PTR(nxFIT)),GKK(1,nxFIT),GR,
     '          GRR,%VAL(GQ_PTR(nxFIT)),PG,RE1,RG,SE,SF,SP,
     '          WD,WDL,WG,WU,XA,XE,XG,XID,XIDL,XIG,XIP,XO(1,nxFIT),
     '          XP,YG,YGF,YP,YQS,ZA,ZA1,Z_CONT,ZD,ZDL,ZE,ZE1,
     '          ZP,ZP1,FIRSTA,FIX(1,1,nxFIT),IN_PLANE,ERROR,*9999)
            ENDDO !no_nrlist (nr)
          ELSE
            nr=NRLIST(1)
            CALL FIT_PATCH(IPIVOT,NELIST,NENP,NKJE,NKH,NPF,NPLIST,
     '        NPLIST2,NPNE,nr,NVHP,NVJE,nxFIT,%VAL(GK_PTR(nxFIT)),
     '        GR,PG,SE,XA,XE,XG,XP,YG,ZP,ERROR,*9999)
          ENDIF
        ENDIF
      ENDIF

      IF(KTYP12.EQ.3) THEN !Strain energy smoothing term

C cpb 8/3/95

        ERROR='>>Code needs updating'
        GOTO 9999

CC ...   Carry current solution to subsequent calls to ZDER via ZP
C        DO nonode=1,NPNODE(0,nr)
C          np=NPNODE(nonode,nr)
C          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
C            DO nv=1,NVJP(nj,np)
C              DO nk=1,NKJ(nj,np)
CC cpb 8/6/94 I think this needs to be generalised with nj_loc
C                ZP(nk,nv,nj,np,nc)=XP(nk,nv,nj+NJT,np)
C              ENDDO
C            ENDDO
C          ENDDO
C        ENDDO
      ENDIF

      CALL EXITS('FIT')
      RETURN
 9999 CALL ERRORS('FIT',ERROR)
      CALL EXITS('FIT')
      RETURN 1
      END


