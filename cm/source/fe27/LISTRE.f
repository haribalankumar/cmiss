      SUBROUTINE LISTRE(IBT,IDO,INP,NAN,NBH,NBJ,
     '  NDDL,NDLT,NEELEM,NELIST,NGLIST,
     '  NHE,NHP,NKH,NKHE,NKJE,NPF,NPNE,NPNODE,NRLIST,
     '  NVHE,NVHP,NVJE,NW,NXLIST,NYNE,NYNP,
     '  CE,CG,CGE,CP,CURVCORRECT,FEXT,PG,RG,SE,WG,XA,XE,XG,XID,XIG,XP,
     '  YG,YP,ZA,ZE,ZG,ZP,STRING,ERROR,*)

C#### Subroutine: LISTRE
C###  Description:
C###    LISTRE lists stress tensors and principal stresses at Gauss
C###    points or at equally spaced points along an arbitrary xi
C###    coordinate line.
C**** PRSTMAX is maximum principle stress
C**** PRSTMIN is minimum principle stress
C****   (absolute value for FE40, actual value for FE50).

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbst02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'stra00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NDDL(NEM,NDEM),NDLT(NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NGLIST(0:NGM),
     '  NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),
     '  NRLIST(0:NRM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NVJE(NNM,NBFM,NJM,NEM),
     '  NW(NEM,3,NXM),NXLIST(0:NXM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM,NXM),CG(NMM,NGM),CGE(NMM,NGM,NEM,NXM),
     '  CP(NMM,NPM,NXM),
     '  CURVCORRECT(2,2,NNM,NEM),FEXT(NIFEXTM,NGM,NEM),
     '  PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM,NEM),
     '  WG(NGM,NBM),XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XG(NJM,NUM),XID(NIM,NDM),
     '  XIG(NIM,NGM,NBM),XP(NKM,NVM,NJM,NPM),YG(NIYGM,NGM,NEM),
     '  YP(NYM,NIYM,NXM),ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZG(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,IPOINTTYP,IXI,N3CO,no_nglist,
     '  no_nrlist,nr,NTPOIN,nx,nxc
      REAL*8 RFROMC,XIPOS(3)
      CHARACTER COORDS*9,FILE*(MXCH),STRESSTYPE*17
      LOGICAL ABBREV,ALL_REGIONS,CBBREV,EXTREMA_ONLY,FULL,OPFILE,
     &  STRAINENERGY 

      CALL ENTERS('LISTRE',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list stress<;FILENAME><;(table/list)[list]>
C###  Parameter:    <at (all/GAUSS_PT#s/data)[all]>
C###    Specifies the gauss points for which the stress is to be listed
C###    for. The 'all' command prompts all gauss points to be listed.
C###  Parameter:    <(fibre/reference/wall)[fibre]>
C###    Specify the coordinate system used to list the stress tensors in
C###  Parameter:    <(total/passive/active/total_no_hydro)[total]>
C###    Specifies whether to output total stress or just passive or
C###    active stress components. Total_no_hydro lists total stress
C###    without the hydrostatic pressure component.
C###  Parameter:    <element (all/GROUP/#s)[all]>
C###    Specifies the element number or groups for which the stress
C###    is to be iutputted. The all command prompts for all currently
C###    defined elements to be included.
C###  Parameter:    <(all/extrema_only)[all]>
C###    Option to only list the extrema. (max and min stress)
C###  Parameter:    <(full/reduced)[reduced]>
C###    Option to list without the principle stress
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Lists stress components wrt fibre/reference/wall coordinates
C###    at specified Gauss points in selected elements to screen
C###    or file (FILENAME.opstre).
C###    TABLE outputs stresses in an abbreviated (tabular) format.
C###    EXTREMA_ONLY option prints out only the max and min principal
C###    stresses.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME><;(table/list)[list]>'
        OP_STRING(2)=BLANK(1:15)//'<at (all/GAUSS_PT#s/data)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<(fibre/reference/wall)[fibre]>'
        OP_STRING(4)=BLANK(1:15)//'<(total/passive/active/'
     &    //'total_no_hydro)[total]>'
        OP_STRING(5)=BLANK(1:15)//'<element (all/GROUP/#s)[all]>'
        OP_STRING(6)=BLANK(1:15)//'<(all/extrema_only)[all]>'
        OP_STRING(7)=BLANK(1:15)//'<(full/reduced)[reduced]>'
        OP_STRING(8)=BLANK(1:15)//'<strain_energy>'
        OP_STRING(9)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(10)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list stress<;FILENAME><;(table/list)[list]>
C###  Parameter:    xi XI_DIRECTION#{>=1,<=3}
C###    Specify the xi direction which the points will be placed
C###    for the stress to be listed.
C###  Parameter:    <xi_1 VALUE#[0.5]{>=0.0,<=1.0}>
C###    Define a location for the stress to be outputted in the
C###    Xi_1 direction for each element specified.
C###  Parameter:    <xi_2 VALUE#[0.5]{>=0.0,<=1.0}>
C###    Define a location for the stress to be outputted in the
C###    Xi_2 direction for each element specified.
C###  Parameter:    <xi_3 VALUE#[0.5]{>=0.0,<=1.0}>
C###    Define a location for the stress to be outputted in the
C###    Xi_3 direction for each element specified.
C###  Parameter:    <points #PTS[5]>
C###    Specifies the number of point in an xi direction
C###  Parameter:    <(fibre/reference/wall)[fibre]>
C###    Specify the coordinate system used to list the stress tensors in
C###  Parameter:    <(total/passive/active/total_no_hydro)[total]>
C###    Specifies whether to output total stress or just passive or
C###    active stress components. Total_no_hydro lists total stress 
C###    without the hydrostatic pressure component.
C###  Parameter:    <element (all/GROUP/#s)[all]>
C###    Specify the element groups. The "all" keyword will use
C###    all element groups.
C###  Parameter:    <(all/extrema_only)[all]>
C###    Specifies whether all values or only maximum and minimum are
C###    listed.
C###  Parameter:    <(full/reduced)[reduced]>
C###    Option to list without the principle stress
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Lists stress components wrt fibre/reference/wall coordinates
C###    at specified Xi locations in selected elements to screen
C###    or file (FILENAME.opstre).
C###    TABLE outputs stresses in an abbreviated (tabular) format.
C###    EXTREMA_ONLY option prints out only the max and min principal
C###    stresses.

        OP_STRING(1)= STRING(1:IEND)//'<;FILENAME><;(table/list)[list]>'
        OP_STRING(2)= BLANK(1:15)//'xi XI_DIRECTION#{>=1,<=3}'
        OP_STRING(3)= BLANK(1:15)//'<xi_1 VALUE#[0.5]{>=0.0,<=1.0}>'
        OP_STRING(4)= BLANK(1:15)//'<xi_2 VALUE#[0.5]{>=0.0,<=1.0}>'
        OP_STRING(5)= BLANK(1:15)//'<xi_3 VALUE#[0.5]{>=0.0,<=1.0}>'
        OP_STRING(6)= BLANK(1:15)//'<points #PTS[5]>'
        OP_STRING(7)= BLANK(1:15)//'<(fibre/reference/wall)[fibre]>'
        OP_STRING(8)= BLANK(1:15)//'<(total/passive/active/'
     &    //'total_no_hydro)[total]>'
        OP_STRING(9)= BLANK(1:15)//'<element (all/GROUP/#s)[all]>'
        OP_STRING(10)=BLANK(1:15)//'<(all/extrema_only)[all]>'
        OP_STRING(11)=BLANK(1:15)//'<(full/reduced)[reduced]>'
        OP_STRING(12)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(13)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LISTRE',ERROR,*9999)
      ELSE
        TABLE=.FALSE.
        IF(NTCOQU(noco).GT.0) THEN
          FILE=COQU(noco,1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          OPFILE=.TRUE.
          IOFI=IOFILE1
          IF(NTCOQU(noco).GT.1) THEN
            IF(ABBREV(COQU(noco,2),'TABLE',1)) TABLE=.TRUE.
          ENDIF
          IF(TABLE) THEN
            CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.ipdata','NEW',
     '       'SEQUEN','FORMATTED',132,ERROR,*9999)
          ELSE
            CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opstre','NEW',
     '        'SEQUEN','FORMATTED',132,ERROR,*9999)
          ENDIF
        ELSE
          OPFILE=.FALSE.
        ENDIF

        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '    ERROR,*9999)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        IF(CBBREV(CO,'EXTREMA_ONLY',2,noco+1,NTCO,N3CO)) THEN
          EXTREMA_ONLY=.TRUE.
        ELSE
          EXTREMA_ONLY=.FALSE.
        ENDIF

        IF(CBBREV(CO,'FULL',3,noco+1,NTCO,N3CO)) THEN
          FULL=.TRUE.
        ELSE
          FULL=.FALSE.
        ENDIF

        IF(CBBREV(CO,'AT',2,noco+1,NTCO,N3CO)) THEN
          IF(CBBREV(CO,'DATA',3,N3CO+1,NTCO,N3CO)) THEN
            IPOINTTYP=2 !data points
          ELSE
            IPOINTTYP=1 !Gauss Points
            CALL PARSIL(CO(N3CO+1),NGM,NGLIST(0),NGLIST(1),ERROR,*9999)
          ENDIF
        ELSE
          IPOINTTYP=1 !Gauss Points
          NGLIST(0)=NGT(NBH(NH_LOC(1,nx),1,NELIST(1)))
          DO no_nglist=1,NGLIST(0)
            NGLIST(no_nglist)=no_nglist
          ENDDO
        ENDIF

        IF(CBBREV(CO,'XI',1,noco+1,NTCO,N3CO)) THEN
          IXI=IFROMC(CO(N3CO+1))
        ELSE
          IXI=0
        ENDIF

        IF(IXI.EQ.1) THEN
          IF(CBBREV(CO,'XI_2',4,noco+1,NTCO,N3CO)) THEN
            XIPOS(2)=RFROMC(CO(N3CO+1))
          ELSE
            XIPOS(2)=0.5D0
          ENDIF
          IF(CBBREV(CO,'XI_3',4,noco+1,NTCO,N3CO)) THEN
            XIPOS(3)=RFROMC(CO(N3CO+1))
          ELSE
            XIPOS(3)=0.5D0
          ENDIF

        ELSE IF(IXI.EQ.2) THEN
          IF(CBBREV(CO,'XI_1',4,noco+1,NTCO,N3CO)) THEN
            XIPOS(1)=RFROMC(CO(N3CO+1))
          ELSE
            XIPOS(1)=0.5D0
          ENDIF
          IF(CBBREV(CO,'XI_3',4,noco+1,NTCO,N3CO)) THEN
            XIPOS(3)=RFROMC(CO(N3CO+1))
          ELSE
            XIPOS(3)=0.5D0
          ENDIF

        ELSE IF(IXI.EQ.3) THEN
          IF(CBBREV(CO,'XI_1',4,noco+1,NTCO,N3CO)) THEN
            XIPOS(1)=RFROMC(CO(N3CO+1))
          ELSE
            XIPOS(1)=0.5D0
          ENDIF
          IF(CBBREV(CO,'XI_2',4,noco+1,NTCO,N3CO)) THEN
            XIPOS(2)=RFROMC(CO(N3CO+1))
          ELSE
            XIPOS(2)=0.5D0
          ENDIF
        ENDIF

        IF(CBBREV(CO,'POINTS',1,noco+1,NTCO,N3CO)) THEN
          NTPOIN=IFROMC(CO(N3CO+1))
          CALL ASSERT(NTPOIN.GT.1,'>>Must be more than one point.',
     '      ERROR,*9999)
        ELSE
          NTPOIN=5
        ENDIF

        IF(CBBREV(CO,'REFERENCE',3,noco+1,NTCO,N3CO)) THEN
          COORDS='Reference'
        ELSE IF(CBBREV(CO,'WALL',3,noco+1,NTCO,N3CO)) THEN
          COORDS='Wall'
        ELSE
          COORDS='Fibre'
        ENDIF

C MPN 24May2000: option for total or just passive or active stress cmpts
        IF(CBBREV(CO,'TOTAL',3,noco+1,NTCO,N3CO)) THEN
          STRESSTYPE='Total'
        ELSE IF(CBBREV(CO,'TOTAL_NO_HYDRO',8,noco+1,NTCO,N3CO)) THEN
          STRESSTYPE='Total_no_hydro'
        ELSE IF(CBBREV(CO,'PASSIVE',3,noco+1,NTCO,N3CO)) THEN
          STRESSTYPE='Passive'
        ELSE IF(CBBREV(CO,'ACTIVE',3,noco+1,NTCO,N3CO)) THEN
          STRESSTYPE='Active'
        ELSE
          STRESSTYPE='Total'
        ENDIF
C newe
C GR 05August2008 news
        IF(CBBREV(CO,'STRAIN_ENERGY',8,noco+1,NTCO,N3CO)) THEN
          STRAINENERGY=.TRUE.
        ELSE
          STRAINENERGY=.FALSE.
        ENDIF
C GR newe
C       Initialise max and min stresses
        PRSTMAX=-1.0d6
        PRSTMIN= 1.0d6

        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          CALL ASSERT((EXTREMA_ONLY.AND.ITYP1(nr,nx).EQ.5).OR.
     '      .NOT.EXTREMA_ONLY,'>>ERROR: option EXTREMA_ONLY '
     '      //'only implemented for FE50 problems',ERROR,*9999)
          CALL OPSTRE(IBT,IDO,INP,IPOINTTYP,IXI,NAN,NBH,NBJ,
     '      NDDL,NDLT,NEELEM,NELIST,
     '      NGLIST,NHE(1,nx),NHP(1,nr,nx),NKH(1,1,1,nr),NKHE,NKJE,
     '      NPF,NPNE,NPNODE,nr,NTPOIN,NVHE,
     '      NVHP(1,1,1,nr),NVJE,NW(1,1,nx),nx,nxc,NYNE,NYNP,
     '      CE(1,1,nx),CG,CGE(1,1,1,nx),
     '      CP(1,1,nx),CURVCORRECT,
     '      FEXT,PG,RG,SE,WG,XA,XE,XG,XID,XIG,XIPOS,XP,
     '      YG,YP(1,1,nx),ZA,ZE,ZG,ZP,COORDS,STRESSTYPE,EXTREMA_ONLY,
     '      FULL,STRAINENERGY,ERROR,*9999)
        ENDDO !no_nrlist (nr)

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
          TABLE=.FALSE.
        ENDIF
      ENDIF

      CALL EXITS('LISTRE')
      RETURN
 9999 CALL ERRORS('LISTRE',ERROR)
      CALL EXITS('LISTRE')
      IF(OPFILE) THEN
        CALL CLOSEF(IOFI,ERROR,*9999)
        IOFI=IOOP
        TABLE=.FALSE.
      ENDIF
      RETURN 1
      END


