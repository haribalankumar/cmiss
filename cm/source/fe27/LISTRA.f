      SUBROUTINE LISTRA(IBT,IDO,INP,NAN,NBH,NBJ,NDDL,NDLT,NEELEM,NELIST,
     '  NGLIST,NHE,NHP,NKH,NKHE,NKJE,NPF,NPNE,NPNODE,
     '  NRLIST,NVHE,NVHP,NVJE,NW,NXLIST,NYNE,NYNP,
     '  CE,CG,CGE,CP,CURVCORRECT,PG,RG,SE,XA,XE,XG,XID,XIG,XP,YP,ZA,
     '  ZE,ZG,ZP,STRING,ERROR,*)

C#### Subroutine: LISTRA
C###  Description:
C###    LISTRA lists strain tensors, invariants and principal strains
C###    at Gauss points or at equally spaced points along an arbitrary
C###    xi coordinate line.
C**** PRSTMAX is maximum principle strain
C**** PRSTMIN is minimum principle strain

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbst02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'stra00.cmn'
      INCLUDE 'time00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NDDL(NEM,NDEM),NDLT(NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NGLIST(0:NGM),NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),
     '  NKH(NHM,NPM,NCM,0:NRM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NRLIST(0:NRM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NVJE(NNM,NBFM,NJM,NEM),
     '  NW(NEM,3,NXM),NXLIST(0:NXM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM,NXM),CG(NMM,NGM),CGE(NMM,NGM,NEM,NXM),
     '  CP(NMM,NPM,NXM),
     '  CURVCORRECT(2,2,NNM,NEM),PG(NSM,NUM,NGM,NBM),RG(NGM),
     '  SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),
     '  XID(NIM,NDM),XIG(NIM,NGM,NBM),XP(NKM,NVM,NJM,NPM),
     '  YP(NYM,NIYM,NXM),ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZG(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,INDEX,
     '  IPOINTTYP,IXI,N3CO,no_nglist,no_nrlist,
     '  nr,NTPOIN,nx,nxc
      REAL*8 RFROMC,XIPOS(3)
      CHARACTER COORDS*9,FILE*(MXCH),TYPE*3
      LOGICAL ABBREV,ALL_REGIONS,CBBREV,EXTREMA_ONLY,FULL,
     '  OPFILE,STATIC,HIGH_PRECISION

      CALL ENTERS('LISTRA',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list strain<;FILENAME><;(table/list)[list]>
C###  Description:
C###    Lists strain components wrt fibre/reference/wall coordinates
C###    at specified Gauss points or at data points in selected elements
C###    to screen or file (FILENAME.opstra) in the directory
C###    specified by PATH with $current specifing the current default file..
C###    TABLE outputs strains in an abbreviated (tabular) format.
C###    EXTREMA_ONLY option prints out only the max and min principal
C###    strains.
C###  Parameter:    <(gauss/data)[gauss]>
C###    Specify wiether strains are listed a data points or gauss points
C###  Parameter:    <at (GAUSS_POINT#s/all)[all]>
C###    Specify the numbers of the gauss points. The "all" keyword will use
C###    all gauss points.
C###  Parameter:    <(fibre/reference/wall)[fibre]>
C###    Specify the coordinate system used to list the strain tensors in
C###  Parameter:    <element (all/GROUP/#s)[all]>
C###    Specify the element groups. The "all" keyword will use
C###    all element groups.
C###  Parameter:    <when (TIME#/static)[static]>
C###    Specify the time in time dependant problems or static for
C###    qausi static problems
C###  Parameter:    <wrt TIME_REF#[0]>
C###    Specify the reference  time
C###  Parameter:    <(normal/high_precision)[normal]>
C###    Specifies whether higher precision output is required
C###  Parameter:    <(all/extrema_only)[all]>
C###    Specifies wiether all values or only maximum and minimum are listed
C###  Parameter:    <(full/reduced)[reduced]>
C###    Option to list without the principle stress
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.



        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME><;(table/list)[list]>'
        OP_STRING(2)=BLANK(1:15)//'<(gauss/data)[gauss]>'
        OP_STRING(3)=BLANK(1:15)//'<at (all/GAUSS_POINT#s)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<(fibre/reference/wall)[fibre]>'
        OP_STRING(5)=BLANK(1:15)//'<element (all/GROUP/#s)[all]>'
        OP_STRING(6)=BLANK(1:15)//'<when (static/TIME#)[static]>'
        OP_STRING(7)=BLANK(1:15)//'<wrt TIME_REF#[0]>'
        OP_STRING(8)=BLANK(1:15)//'<(normal/high_precision)[normal]>'
        OP_STRING(9)=BLANK(1:15)//'<(all/extrema_only)[all]>'
        OP_STRING(10)=BLANK(1:15)//'<(full/reduced)[reduced]>'
        OP_STRING(11)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(12)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list strain<;FILENAME><;(table/list)[list]>
C###  Parameter:    xi XI_DIRECTION#{>=1,<=3}
C###    Specifies the xi direction along which points will be
C###    positioned so the strain can be listed.
C###  Parameter:    <xi_1 VALUE#[0.5]{>=0.0,<=1.0}>
C###    Specifies the value of the xi 1 position
C###  Parameter:    <xi_2 VALUE#[0.5]{>=0.0,<=1.0}>
C###    Specifies the value of the xi 2 position
C###  Parameter:    <xi_3 VALUE#[0.5]{>=0.0,<=1.0}>
C###    Specifies the value of the xi 3 position
C###  Parameter:    <points #POINTS[5]>
C###    Specifies the number of point in an xi direction
C###  Parameter:    <(fibre/reference/wall)[fibre]>
C###    Specify the coordinate system used to list the strain tensors in
C###  Parameter:    <element (all/GROUP/#s)[all]>
C###    Specify the element groups. The "all" keyword will use
C###    all element groups.
C###  Parameter:    <(normal/high_precision)[normal]>
C###    Specifies whether higher precision output is required
C###  Parameter:    <(all/extrema_only)[all]>
C###    Specifies whether all values or only maximum and minimum are
C###    listed.
C###  Parameter:    <(full/reduced)[reduced]>
C###    Option to list without the principle stress
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Lists strain components wrt fibre/reference/wall coordinates
C###    at specified Xi locations in selected elements to screen
C###    or file (FILENAME.opstra).
C###    TABLE outputs strains in an abbreviated (tabular) format.
C###    EXTREMA_ONLY option prints out only the max and min principal
C###    strains.

        OP_STRING(1)= STRING(1:IEND)//'<;FILENAME><;(table/list)[list]>'
        OP_STRING(2)= BLANK(1:15)//'xi XI_DIRECTION#{>=1,<=3}'
        OP_STRING(3)= BLANK(1:15)//'<xi_1 VALUE#[0.5]{>=0.0,<=1.0}>'
        OP_STRING(4)= BLANK(1:15)//'<xi_2 VALUE#[0.5]{>=0.0,<=1.0}>'
        OP_STRING(5)= BLANK(1:15)//'<xi_3 VALUE#[0.5]{>=0.0,<=1.0}>'
        OP_STRING(6)= BLANK(1:15)//'<points #POINTS[5]>'
        OP_STRING(7)= BLANK(1:15)//'<(fibre/reference/wall)[fibre]>'
        OP_STRING(8)= BLANK(1:15)//'<element (all/GROUP/#s)[all]>'
        OP_STRING(9)=BLANK(1:15)//'<(normal/high_precision)[normal]>'
        OP_STRING(10)= BLANK(1:15)//'<(all/extrema_only)[all]>'
        OP_STRING(11)=BLANK(1:15)//'<(full/reduced)[reduced]>'
        OP_STRING(12)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(13)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
C#### Command: FEM list strain<;FILENAME><;(table/list)[list]>
C###  Parameter:    xi_point
C###    Specifies the strain is to listed for one point
C###  Parameter:    <xi_1 VALUE#[0.5]{>=0.0,<=1.0}>
C###    Specifies the value of the xi 1 position
C###  Parameter:    <xi_2 VALUE#[0.5]{>=0.0,<=1.0}>
C###    Specifies the value of the xi 2 position
C###  Parameter:    <xi_3 VALUE#[0.5]{>=0.0,<=1.0}>
C###    Specifies the value of the xi 3 position
C###  Parameter:    <(fibre/reference/wall)[fibre]>
C###    Specify the coordinate system used to list the strain tensors in
C###  Parameter:    <element (all/GROUP/#s)[all]>
C###    Specify the element groups. The "all" keyword will use
C###    all element groups.
C###  Parameter:    <(normal/high_precision)[normal]>
C###    Specifies whether higher precision output is required
C###  Parameter:    <(all/extrema_only)[all]>
C###    Specifies whether all values or only maximum and minimum are
C###    listed.
C###  Parameter:    <(full/reduced)[reduced]>
C###    Option to list without the principle stress
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Lists strain components wrt fibre/reference/wall coordinates
C###    at specified Xi locations in selected elements to screen
C###    or file (FILENAME.opstra).
C###    TABLE outputs strains in an abbreviated (tabular) format.
C###    EXTREMA_ONLY option prints out only the max and min principal
C###    strains.

        OP_STRING(1)= STRING(1:IEND)//'<;FILENAME><;(table/list)[list]>'
        OP_STRING(2)= BLANK(1:15)//'xi_point'
        OP_STRING(3)= BLANK(1:15)//'<xi_1 VALUE#[0.5]{>=0.0,<=1.0}>'
        OP_STRING(4)= BLANK(1:15)//'<xi_2 VALUE#[0.5]{>=0.0,<=1.0}>'
        OP_STRING(5)= BLANK(1:15)//'<xi_3 VALUE#[0.5]{>=0.0,<=1.0}>'
        OP_STRING(6)= BLANK(1:15)//'<(fibre/reference/wall)[fibre]>'
        OP_STRING(7)= BLANK(1:15)//'<element (all/GROUP/#s)[all]>'
        OP_STRING(8)=BLANK(1:15)//'<(normal/high_precision)[normal]>'
        OP_STRING(9)= BLANK(1:15)//'<(all/extrema_only)[all]>'
        OP_STRING(10)=BLANK(1:15)//'<(full/reduced)[reduced]>'
        OP_STRING(11)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(12)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LISTRA',ERROR,*9999)
      ELSE
        TABLE=.FALSE.
        IF(NTCOQU(noco).GT.0) THEN
          FILE=COQU(noco,1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          OPFILE=.TRUE.
          IOFI=IOFILE1
          IF(NTCOQU(noco).GT.1) THEN
            IF(ABBREV(COQU(noco,2),'TABLE',1)) THEN
              TABLE=.TRUE.
            ENDIF
          ENDIF
          IF(TABLE) THEN
            CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.ipdata','NEW',
     '       'SEQUEN','FORMATTED',132,ERROR,*9999)
          ELSE
            CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opstra','NEW',
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

C TVK 16Nov1999: Option for High_Precision output
        IF(CBBREV(CO,'HIGH_PRECISION',2,noco+1,NTCO,N3CO)) THEN
          HIGH_PRECISION=.TRUE.
        ELSE
          HIGH_PRECISION=.FALSE.
        ENDIF

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

        IXI=0
        IF(CBBREV(CO,'DATA',1,noco+1,NTCO,N3CO)) THEN
          IPOINTTYP=2 !data points
        ELSE IF(CBBREV(CO,'XI_POINT',4,noco+1,NTCO,N3CO)) THEN
          IPOINTTYP=3 !single xi point
          IXI=4
        ELSE IF(CBBREV(CO,'XI',1,noco+1,NTCO,N3CO)) THEN
          IPOINTTYP=0 !xi-coordinate line
          IXI=IFROMC(CO(N3CO+1))
        ELSE
          IPOINTTYP=1 !Gauss Points
        ENDIF

        IF(IPOINTTYP.EQ.3) THEN
          IF(CBBREV(CO,'E11',3,noco+1,NTCO,N3CO)) THEN
            TYPE='E11'
          ELSE IF(CBBREV(CO,'E22',3,noco+1,NTCO,N3CO)) THEN
            TYPE='E22'
          ELSE IF(CBBREV(CO,'E33',3,noco+1,NTCO,N3CO)) THEN
            TYPE='E33'
          ELSE IF(CBBREV(CO,'E12',3,noco+1,NTCO,N3CO)) THEN
            TYPE='E12'
          ELSE IF(CBBREV(CO,'E13',3,noco+1,NTCO,N3CO)) THEN
            TYPE='E13'
          ELSE IF(CBBREV(CO,'E23',3,noco+1,NTCO,N3CO)) THEN
            TYPE='E23'
          ELSE
            TYPE='E11'
          ENDIF
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
        ELSE IF(IPOINTTYP.EQ.3) THEN
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
          IF(CBBREV(CO,'XI_3',4,noco+1,NTCO,N3CO)) THEN
            XIPOS(3)=RFROMC(CO(N3CO+1))
          ELSE
            XIPOS(3)=0.5D0
          ENDIF
        ENDIF

        IF(CBBREV(CO,'POINTS',1,noco+1,NTCO,N3CO)) THEN
          NTPOIN=IFROMC(CO(N3CO+1))
        ELSE
          IF(IPOINTTYP.EQ.3) THEN
            NTPOIN=1
          ELSE
            NTPOIN=5
          ENDIF
        ENDIF

        IF(CBBREV(CO,'YP',1,noco+1,NTCO,N3CO)) THEN
          INDEX=IFROMC(CO(N3CO+1))
        ELSE
          INDEX=1
        ENDIF

        IF(CBBREV(CO,'REFERENCE',3,noco+1,NTCO,N3CO)) THEN
          COORDS='Reference'
        ELSE IF(CBBREV(CO,'WALL',3,noco+1,NTCO,N3CO)) THEN
          COORDS='Wall'
        ELSE
          COORDS='Fibre'
        ENDIF

        IF(CBBREV(CO,'WHEN',2,noco+1,NTCO,N3CO)) THEN
          TIME=RFROMC(CO(N3CO+1))
          STATIC=.FALSE.
        ELSE
          STATIC=.TRUE.
        ENDIF

        IF(CBBREV(CO,'WRT',2,noco+1,NTCO,N3CO)) THEN
          TIME_REF=RFROMC(CO(N3CO+1))
        ELSE
          TIME_REF=0.0D0
        ENDIF

        IF(CBBREV(CO,'AT',1,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NGM,NGLIST(0),NGLIST(1),ERROR,*9999)
        ELSE
          IF(STATIC)THEN
            NGLIST(0)=NGT(NBH(NH_LOC(1,nx),1,NELIST(1)))
            DO no_nglist=1,NGLIST(0)
              NGLIST(no_nglist)=no_nglist
            ENDDO
          ELSE IF(.NOT.STATIC)THEN !use geometric gauss points
            NGLIST(0)=NGT(NBJ(1,1))
            DO no_nglist=1,NGLIST(0)
              NGLIST(no_nglist)=no_nglist
            ENDDO
          ENDIF
        ENDIF

C       Initialise max and min strains
        PRSTMAX=-1.0D6
        PRSTMIN= 1.0D6

        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          CALL OPSTRA(IBT,IDO,INDEX,INP,
     '      IPOINTTYP,IXI,NAN,NBH,NBJ,NDDL,NDLT,
     '      NEELEM,NELIST,NGLIST,NHE(1,nx),NHP(1,nr,nx),
     '      NKH(1,1,1,nr),NKHE,NKJE,NPF,NPNE,NPNODE,nr,NTPOIN,
     '      NVHE,NVHP(1,1,1,nr),NVJE,NW(1,1,nx),nx,nxc,NYNE,NYNP,
     '      CE(1,1,nx),CG,CGE(1,1,1,nx),
     '      CP(1,1,nx),CURVCORRECT,PG,RG,SE,XA,XE,
     '      XG,XID,XIG,XIPOS,XP,YP,ZA,ZE,ZG,ZP,
     '      COORDS,STATIC,HIGH_PRECISION,EXTREMA_ONLY,FULL,
     '      TYPE,ERROR,*9999)
        ENDDO

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
          TABLE=.FALSE.
        ENDIF
      ENDIF

      CALL EXITS('LISTRA')
      RETURN
 9999 CALL ERRORS('LISTRA',ERROR)
      CALL EXITS('LISTRA')
      IF(OPFILE) THEN
        CALL CLOSEF(IOFI,ERROR,*9999)
        IOFI=IOOP
        TABLE=.FALSE.
      ENDIF
      RETURN 1
      END


