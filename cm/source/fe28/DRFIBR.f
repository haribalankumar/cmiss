      SUBROUTINE DRFIBR(IBT,IDO,INP,ISEG,ISFIBR,MXI,NAN,NBH,NBJ,
     '  NEELEM,NELIST,NENP,NHE,NHP,NKH,NKHE,NKJE,NNB,NPF,NPNE,
     '  NPNODE,NRLIST,NVHE,NVHP,NVJE,NW,NXI,NXLIST,NYNE,NYNP,
     '  CURVCORRECT,SE,XA,XE,XG,XP,YP,ZA,ZE,ZP,CSEG,STRING,ERROR,*)

C#### Subroutine: DRFIBR
C###  Description:
C###    DRFIBR draws element fibre orientations.  Constant vectors of
C###    Xi-coordinate length DXIF are drawn on plane Xi(3)=XIF at
C###    increments of DXI1 & DXI2 in Xi(1) and Xi(2) direction.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'fibr00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISEG(*),ISFIBR(NWM,NEM,NGRSEGM),MXI(2,NEM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),NENP(NPM,0:NEPM,0:NRM),
     '  NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NNB(4,4,4,NBFM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NRLIST(0:NRM),NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),NXLIST(0:NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM),ZA(NAM,NHM,NCM,NEM),
     '  ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG1,IEND,IEND1,INDEX,INDEX_POLYLINE,iw,
     '  IWK(6),N3CO,nb,ne,NEE,NEE_CURRENT,ni,noiw,nolist,no_nrlist,
     '  nr,NTDXI,NTIW,nx,nxc
      REAL*8 DXI(3),RFROMC,XI(3),XI_INITIAL(3)
      CHARACTER TYPE*5
      LOGICAL ALL_REGIONS,CBBREV,CONTINUE,DEFORMED

      CALL ENTERS('DRFIBR',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM draw fibres <region #>[1]
C###  Parameter:           <elements (all/#s/GROUP)[all]>
C###    Define what elements are to be included.
C###  Parameter:           <xi_3 VALUE#[0.0]{>=0.0,<=1.0}>
C###    Distance along Xi_3 direction where fibres will be drawn.
C###  Parameter:           <dxi VALUE#s[0.2,0.2,0.13]>
C###    Fibre separation in Xi(1,2 and 3) directions.
C###  Parameter:           <on (all/WS#s)[all]>
C###    Specify the worksation (GX window) to draw the
C###    points on.
C###  Parameter:           <rgb=RGB[blue]>
C###    Specify the colour of the data.  The colours
C###    allowed are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY
C###  Parameter:           <(undeformed/deformed)[undeformed]>
C###    Draw undeformed/deformed geometry
C###  Parameter:           <regions (all/#s/all)[all]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:           <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Draws element fibre orientations.  Constant vectors of
C###    Xi-coordinate length xi_3 are drawn on plane xi_3 at
C###    increments of dxi(1) & dxi(2) in Xi(1) and Xi(2) direction.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<region #[1]>'
        OP_STRING(3)=BLANK(1:15)//'<elements (all/#s/GROUP)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<xi_3 VALUE#[0.0]{>=0.0,<=1.0}>'
        OP_STRING(5)=BLANK(1:15)//'<dxi VALUE#s[0.2,0.2,0.13]>'
        OP_STRING(6)=BLANK(1:15)//'<on (all/WS#s)[all]>'
        OP_STRING(7)=BLANK(1:15)//'<rgb=RGB[blue]>'
        OP_STRING(8)=BLANK(1:15)//'<(undeformed/deformed)[undeformed]>'
        OP_STRING(9)=BLANK(1:15)//'<regions (#s/all)[all]>'
        OP_STRING(10)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM draw fibres sheet/helix
C###  Parameter:           <elements (all/#s/GROUP)[all]>
C###    Define what elements are to be included.
C###  Parameter:           <at XI_COORD#s[1.0,1.0,1.0]{>=0.0,<=1.0}>
C###    Specify the Xi coordinates to start at (not relivant for helix)
C###  Parameter:           <step XI#[0.01]>
C###    Specify the step in the XI# direction
C###  Parameter:           <on (all/WS#s)[all]>
C###    Specify the worksation (GX window) to draw the
C###    points on.
C###  Parameter:           <rgb=RGB[blue]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Parameter:           <(undeformed/deformed)[undeformed]>
C###  Parameter:           <regions (#s/all)[all]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:           <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Draws element fibre orientations. The sheet option specifies
C###    fibres to be draw in the plan of the sheet, the helix option
C###    draws a helix through the tangent field.

        OP_STRING(1)=STRING(1:IEND)//' sheet/helix'
        OP_STRING(2)=BLANK(1:15)//'<elements (all/#s/GROUP)[all]>'
        OP_STRING(3)=BLANK(1:15)
     '    //'<at XI_COORD#s[1.0,1.0,1.0]{>=0.0,<=1.0}>'
        OP_STRING(4)=BLANK(1:15)//'<step XI#[0.01]>'
        OP_STRING(5)=BLANK(1:15)//'<on (all/WS#s)[all]>'
        OP_STRING(6)=BLANK(1:15)//'<rgb=RGB[blue]>'
        OP_STRING(7)=BLANK(1:15)//'<(undeformed/deformed)[undeformed]>'
        OP_STRING(8)=BLANK(1:15)//'<regions (#s/all)[all]>'
        OP_STRING(9)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRFIBR',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        CALL ASSERT(CALL_ELFB,'>>No fibre elements defined',ERROR,*9999)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '    ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        IF(CBBREV(CO,'SHEET',2,noco+1,NTCO,N3CO)) THEN
          TYPE='SHEET'
        ELSE IF(CBBREV(CO,'HELIX',2,noco+1,NTCO,N3CO)) THEN
          TYPE='HELIX'
        ELSE
          TYPE='AXIS'
        ENDIF

        IF(TYPE(1:4).EQ.'AXIS') THEN
          IF(CBBREV(CO,'XI_3',1,noco+1,NTCO,N3CO)) THEN
            XIF=RFROMC(CO(N3CO+1))
          ELSE
            XIF=0.0D0
          ENDIF
          IF(CBBREV(CO,'DXI',2,noco+1,NTCO,N3CO)) THEN
            CALL PARSRL(CO(N3CO+1),3,NTDXI,DXI,ERROR,*9999)
            DXI1=DXI(1)
            DXI2=DXI(2)
            IF(NTDXI.GT.2) THEN
              DXIF=DXI(3)
            ELSE
              DXIF=0.13D0
            ENDIF
          ELSE
            DXI1=0.20D0
            DXI2=0.20D0
            DXIF=0.13D0
          ENDIF

        ELSE IF(TYPE(1:5).EQ.'SHEET') THEN
          IF(CBBREV(CO,'AT',1,noco+1,NTCO,N3CO)) THEN
            CALL PARSRL(CO(N3CO+1),3,NTDXI,DXI,ERROR,*9999)
            DXI1=DXI(1) !are the starting Xi coords
            DXI2=DXI(2)
          ELSE
            DXI1=0.0D0
            DXI2=0.0D0
          ENDIF

          IF(CBBREV(CO,'STEP',2,noco+1,NTCO,N3CO)) THEN
            DXIF=RFROMC(CO(N3CO+1)) !is Xi step length
          ELSE
            DXIF=0.01D0
          ENDIF

        ELSE IF(TYPE(1:5).EQ.'HELIX') THEN
          IF(CBBREV(CO,'AT',1,noco+1,NTCO,N3CO)) THEN
            CALL PARSRL(CO(N3CO+1),3,NTDXI,XI_INITIAL,ERROR,*9999)
          ELSE
            XI_INITIAL(1)=1.0D0 !are the starting Xi coords
            XI_INITIAL(2)=1.0D0
            XI_INITIAL(3)=1.0D0
          ENDIF

          IF(CBBREV(CO,'STEP',2,noco+1,NTCO,N3CO)) THEN
            DXIF=RFROMC(CO(N3CO+1)) !is Xi step length
          ELSE
            DXIF=0.01D0
          ENDIF
        ENDIF

        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
        ELSE
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLUE')
        ENDIF

        IF(CBBREV(CO,'DEFORMED',2,noco+1,NTCO,N3CO)) THEN
          DEFORMED=.TRUE.
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '      ERROR,*9999)
        ELSE
          nx=1 !temporary (for ityps)
          DEFORMED=.FALSE.
        ENDIF

        IF(ADD) THEN
          NTFIBR=NTFIBR+1
        ELSE IF(NTFIBR.EQ.0) THEN
          NTFIBR=1
        ENDIF
        CALL ASSERT(NTFIBR.LE.NRM,'>>NRM too small',ERROR,*9999)
        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          CALL ASSERT(NJ_LOC(NJL_FIBR,0,nr).GT.0,'>>No fibres defined',
     '      ERROR,*9999)
          IF(ITYP2(nr,nx).EQ.14.OR.ITYP2(nr,nx).EQ.15.OR.DEFORMED) THEN
            CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '        NKH(1,1,1,nr),NPNODE,
     '        nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,
     '        ERROR,*9999)
          ENDIF
          DO noiw=1,NTIW
            IW=IWK(noiw)
            DO nolist=1,NELIST(0)
              ne=NELIST(nolist)
              nb=NBJ(1,ne)
              IF(NIT(nb).GT.1) THEN
                IF(TYPE(1:4).EQ.'AXIS') THEN
                  CALL ACWK(iw,1,ERROR,*9999)
                  CALL SGFIBR(INDEX,IBT,IDO,INP,ISEG,
     '              ISFIBR(iw,ne,NTFIBR),iw,MXI(1,ne),NAN,NBJ(1,ne),ne,
     '              NKHE(1,1,1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '              nr,NVHE(1,1,1,ne),NVJE(1,1,1,ne),
     '              NW(ne,1,nx),nx,CSEG,
     '              CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '              TYPE,XA,XE,XG,XI,XP,ZA(1,1,1,ne),ZE,ZP,
     '              DEFORMED,ERROR,*9999)
                  CALL DAWK(iw,1,ERROR,*9999)
                ELSE IF(TYPE(1:5).EQ.'SHEET') THEN
                  CALL ACWK(iw,1,ERROR,*9999)
                  CALL SGFIBR(INDEX,IBT,IDO,INP,ISEG,
     '              ISFIBR(iw,ne,NTFIBR),iw,MXI(1,ne),NAN,NBJ(1,ne),ne,
     '              NKHE(1,1,1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '              nr,NVHE(1,1,1,ne),NVJE(1,1,1,ne),
     '              NW(ne,1,nx),nx,CSEG,CURVCORRECT(1,1,1,ne),
     '              SE(1,1,ne),TYPE,XA,XE,XG,XI,XP,ZA(1,1,1,ne),ZE,ZP,
     '              DEFORMED,ERROR,*9999)
                  CALL DAWK(iw,1,ERROR,*9999)
                ELSE IF(TYPE(1:5).EQ.'HELIX') THEN
                  CONTINUE=.TRUE.
                  XI(1)=XI_INITIAL(1) !is initial Xi_1 postion
                  XI(2)=XI_INITIAL(2) !is initial Xi_2 postion
                  XI(3)=XI_INITIAL(3) !is initial Xi_3 postion
                  DIRECTION=1.0D0
                  NEE=ne
                  WRITE(OP_STRING,'('' Start helix from element '',I4,'
     '              //''' Xi coords:'',3F7.3)') NEE,(XI(ni),ni=1,3)
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  DO WHILE (CONTINUE)
                    CALL ACWK(iw,1,ERROR,*9999)
                    CALL SGFIBR(INDEX,IBT,IDO,INP,ISEG,
     '                ISFIBR(iw,NEE,NTFIBR),iw,MXI(1,NEE),NAN,
     '                NBJ(1,NEE),NEE,NKHE(1,1,1,NEE),NKJE(1,1,1,NEE),
     '                NPF(1,1),NPNE(1,1,NEE),nr,NVHE(1,1,1,NEE),
     '                NVJE(1,1,1,NEE),NW(NEE,1,nx),nx,CSEG,
     '                CURVCORRECT(1,1,1,nee),SE(1,1,NEE),
     '                TYPE,XA,XE,XG,XI,XP,ZA(1,1,1,NEE),ZE,ZP,DEFORMED,
     '                ERROR,*9999)
                    CALL DAWK(iw,1,ERROR,*9999)
                    CALL NENXI(NBJ,NEELEM,NENP,NNB,NPNE,NXI,
     '                ERROR,*9999)
                    NEE_CURRENT=NEE
                    WRITE(OP_STRING,'('' Current element is '',I4)')
     '                NEE_CURRENT
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    IF(XI(2).LT.0.0D0) THEN
                      NEE=NXI(-2,1,NEE)
                      WRITE(OP_STRING,
     '                  '('' Adj. element in -Xi(2) dir: '',I4)') NEE
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      XI(2)=1.0D0
                    ELSE IF(XI(2).GT.1.0D0) THEN
                      NEE=NXI( 2,1,NEE)
                      WRITE(OP_STRING,
     '                  '('' Adj. element in  Xi(2) dir: '',I4)') NEE
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      XI(2)=0.0D0
                    ELSE IF(XI(1).LT.0.0D0) THEN
                      NEE=NXI(-1,1,NEE)
                      WRITE(OP_STRING,
     '                  '('' Adj. element in -Xi(1) dir: '',I4)') NEE
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      XI(1)=1.0D0
                    ELSE IF(XI(1).GT.1.0D0) THEN
                      NEE=NXI( 1,1,NEE)
                      WRITE(OP_STRING,
     '                  '('' Adj. element in  Xi(1) dir: '',I4)') NEE
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      XI(1)=0.0D0
                    ELSE IF(XI(3).LT.0.0D0) THEN
                      NEE=NXI(-3,1,NEE)
                      WRITE(OP_STRING,
     '                  '('' Adj. element in -Xi(3) dir: '',I4)') NEE
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      XI(3)=1.0D0
                    ELSE IF(XI(3).GT.1.0D0) THEN
                      NEE=NXI( 3,1,NEE)
                      WRITE(OP_STRING,
     '                  '('' Adj. element in  Xi(3) dir: '',I4)') NEE
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      XI(3)=0.0D0
                    ENDIF
                    IF(NEE.EQ.0) CONTINUE=.FALSE.
                  ENDDO
                ENDIF !TYPE
              ENDIF
            ENDDO !ne
          ENDDO !iw
        ENDDO !nr
      ENDIF

      CALL EXITS('DRFIBR')
      RETURN
 9999 CALL ERRORS('DRFIBR',ERROR)
      CALL EXITS('DRFIBR')
      RETURN 1
      END


