      SUBROUTINE DRDATA(IBT,IDO,INP,ISDANO,ISDAPR,ISDATA,ISDATR,ISEG,
     '  LD,LN,MXI,NBJ,NBH,NDDL,NDLT,NDP,NEELEM,NELIST,NHP,
     '  NKH,NKHE,NKJE,NPF,NPNE,NPNODE,NRE,NRLIST,NVHE,NVHP,
     '  NVJE,NW,NYNE,NYNP,
     '  CE,CG,CGE,CP,CURVCORRECT,PG,SE,WD,WDL,XA,XE,XG,XID,XIDL,XP,
     '  YP,ZA,ZD,ZDD,ZDL,ZE,ZP,CSEG,DATYPE,STRING,ERROR,*)

C#### Subroutine: DRDATA
C###  Description:
C###    DRDATA draws geometric data points.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'data00.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'hist00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'moti00.cmn'
      INCLUDE 'ntsg00.cmn'
      INCLUDE 'scal01.cmn'
      INCLUDE 'time00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISDANO(NWM,NEM),ISDAPR(NWM,NEM),
     '  ISDATA(NWM,NGRSEGM),ISDATR(NWM,NEM),
     '  ISEG(*),LD(NDM),LN(0:NEM),MXI(2,NEM),
     '  NBJ(NJM,NEM),NBH(NHM,NCM,NEM),NDDL(NEM,NDEM),NDLT(NEM),NDP(NDM),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),NHP(NPM,0:NRM,NXM),
     '  NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),NRE(NEM),
     '  NRLIST(0:NRM),NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM,NXM),CG(NMM,NGM),CGE(NMM,NGM,NEM,NXM),
     '  CP(NMM,NPM,NXM),
     '  CURVCORRECT(2,2,NNM,NEM),PG(NSM,NUM,NGM,NBM),SE(NSM,NBFM,NEM),
     '  WD(NJM,NDM),WDL(NHM,NDEM),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XG(NJM,NUM),XID(NIM,NDM),XIDL(NIM,NDEM),
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM),ZA(NAM,NHM,NCM,NEM),
     '  ZD(NJM,NDM),ZDD(NJM,NDM),ZDL(NHM,NDEM),ZE(NSM,NHM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),DATYPE(*)*(*),
     '  ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,INDEX,INDEX_POLYLINE,
     '  INDEX_POLYMARKER,iw,IW_DEFAULT,IWK(6),IWSLICE,N3CO,
     '  NCLUSTER(2),noiw,NO_NDHIST(2),
     '  NO_NJHIST(2),nonr,nr,NTIL,NTIW,nx
      REAL*8 NO_STEP(2),RANGE,RFROMC,STEP,VAL,ZDLMIN
      CHARACTER TITLE*80,TYPE*20
      LOGICAL ABBREV,ALL_REGIONS,CBBREV,STATIC

      CALL ENTERS('DRDATA',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM draw data
C###  Description:
C###    Draws geometric data points.
C###  Parameter:    <(geometry/fibre/field)[geometry]>
C###    Draw either the geometric positions of the data points or
C###    a line segement indicating the fibre angle at the data
C###    points or colour the data points according to a field value.
C###  Parameter:    <elements (all/#s/GROUP)[all]>
C###    Specify which elements to draw data points in. The default
C###    is to draw data points in all elements. An element list
C###    or group name may be supplied.
C###  Parameter:    <from XI_3#[0.0]{>=0.0,<=1.0}>
C###    When drawing 3d data sets in the 2d window this parameter
C###    specifies the minimum xi 3 value and only data points with
C###    a higher xi 3 value will be drawn.
C###  Parameter:    <to XI_3#[1.0]{>=0.0,<=1.0}>
C###    When drawing 3d data sets in the 2d window this parameter
C###    specifies the maximum xi 3 value and only data points with
C###    a lower xi 3 value will be drawn.
C###  Parameter:    <(undeformed/deformed)[undeformed]>
C###    Specify whether the data points are to be drawn in their
C###    original positions or in their updated positions after a
C###    mechanics solution or a fit etc.
C###  Parameter:    <scale FACTOR#[1]>
C###    This parameter scales the size of the data graphics.
C###  Parameter:    <when (static/TIME#)[static]>
C###    Specify if the problem is static and if it is not a time
C###    value can be given at which the data information will be
C###    drawn.
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify which window to draw the data on. The default
C###    is to draw the data on all windows.
C###  Parameter:    <rgb=(RGB/computed)[blue]>
C###    Specify the colour of the data.  The colours
C###    allowed are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###    If the computed option is specified the colour of the data
C###    points is calculated from the value stored at the data point.
C###  Parameter:    <regions (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.

        OP_STRING(1)=STRING(1:IEND)
     '    //' <(geometry/fibre/field)[geometry]>'
        OP_STRING(2)=BLANK(1:15)//'<elements (#s/GROUP)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<from XI_3#[0.0]{>=0.0.<=1.0}>'
        OP_STRING(4)=BLANK(1:15)//'<to XI_3#[1.0]{>=0.0,<=1.0}>'
        OP_STRING(5)=BLANK(1:15)//'<(undeformed/deformed)[undeformed]>'
        OP_STRING(6)=BLANK(1:15)//'<scale FACTOR#[1]>'
        OP_STRING(7)=BLANK(1:15)//'<when (static/TIME)[static]>'
        OP_STRING(8)=BLANK(1:15)//'<on (all/WS#s)[all]>'
        OP_STRING(9)=BLANK(1:15)//'<rgb=(RGB/computed)[blue]>'
        OP_STRING(10)=BLANK(1:15)//'<regions (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM draw data sheet
C###  Description:
C###    ??
C###  Parameter:    <theta=ANGLE#[0.0]{degrees}>
C###  Parameter:    <range D_THETA#[10.0]{degrees}>

        OP_STRING(1)=STRING(1:IEND)//' sheet'
        OP_STRING(2)=BLANK(1:15)//'<theta=ANGLE#[0.0]{degrees}>'
        OP_STRING(3)=BLANK(1:15)//'<range D_THETA#[10.0]{degrees}>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM draw data numbers/values/trace
C###  Description:
C###    Draw data. Numbers draws the data point numbers in their
C###    geometric position on the window. Values draws the value
C###    of each data point on the screen as text. Trace draws a
C###    small circle at the data point projection.
C###  Parameter:    <elements (all/#s/GROUP)[all]>
C###    Specify which elements to draw data points in. The default
C###    is to draw data points in all elements. An element list
C###    or group name may be supplied.
C###  Parameter:    <from XI_3#[0.0]{>=0.0,<=1.0}>
C###    When drawing 3d data sets in the 2d window this parameter
C###    specifies the minimum xi 3 value and only data points with
C###    a higher xi 3 value will be drawn.
C###  Parameter:    <to XI_3#[1.0]{>=0.0,<=1.0}>
C###    When drawing 3d data sets in the 2d window this parameter
C###    specifies the maximum xi 3 value and only data points with
C###    a lower xi 3 value will be drawn.
C###  Parameter:    <when (static/TIME)[static]>
C###    Specify if the problem is static and if it is not a time
C###    value can be given at which the data information will be
C###    drawn.
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify which window to draw the data on. The default
C###    is to draw the data on all windows.
C###  Parameter:    <rgb=RGB[cyan]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Parameter:    <regions (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.

        OP_STRING(1)=STRING(1:IEND)//' numbers/values/trace'
        OP_STRING(2)=BLANK(1:15)//'<elements (all/#s/GROUP)[all]>'
        OP_STRING(3)=BLANK(1:15)
     '    //'<from XI_3#[0.0]{>=0.0}> <to XI_3#[1.0]{<=1.0}>'
        OP_STRING(4)=BLANK(1:15)//'<when (static/TIME)[static]>'
        OP_STRING(5)=BLANK(1:15)//'<on (all/WS#s)[all]>'
        OP_STRING(6)=BLANK(1:15)//'<rgb=RGB[cyan]>'
        OP_STRING(7)=BLANK(1:15)//'<regions (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM draw data parameter NUMBER
C###  Description:
C###    ??
C###  Parameter:    <elements (all/#s/GROUP#)[all]>
C###    Specify which elements to draw data points in. The default
C###    is to draw data points in all elements. An element list
C###    or group name may be supplied.
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify which window to draw the data on. The default
C###    is to draw the data on all windows.
C###  Parameter:    <scale FACTOR#[1]>
C###    This parameter scales the size of the data graphics.
C###  Parameter:    <rgb=RGB[cyan]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Parameter:    <regions (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.

        OP_STRING(1)=STRING(1:IEND)//' parameter NUMBER'
        OP_STRING(2)=BLANK(1:15)//'<elements (all/#s/GROUP#)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<on (all/WS#s)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<scale FACTOR#[1]>'
        OP_STRING(5)=BLANK(1:15)//'<rgb=RGB[cyan]>'
        OP_STRING(6)=BLANK(1:15)//'<regions (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM draw data projections
C###  Description:
C###    This command draws the projections of data points onto
C###    elements.
C###  Parameter:    <elements (all/#s/GROUP#)[all]>
C###    Specify which elements to draw data projections in. The default
C###    is to draw data projections in all elements. An element list
C###    or group name may be supplied.
C###  Parameter:    <from XI_3#[0.0]{>=0.0}>
C###    When drawing 3d data sets in the 2d window this parameter
C###    specifies the minimum xi 3 value and only data points with
C###    a higher xi 3 value will be drawn.
C###  Parameter:    <to XI_3#[1.0]{<=1.0}>
C###    When drawing 3d data sets in the 2d window this parameter
C###    specifies the maximum xi 3 value and only data points with
C###    a lower xi 3 value will be drawn.
C###  Parameter:    <when (static/TIME)[static]>
C###    Specify if the problem is static and if it is not a time
C###    value can be given at which the data information will be
C###    drawn.
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify which window to draw the data on. The default
C###    is to draw the data on all windows.
C###  Parameter:    <rgb=RGB[cyan]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Parameter:    <greater VALUE#[0.0]>
C###    Specify a minimum length of data point projection. Only
C###    projections greater than this length are drawn.
C###  Parameter:    <regions (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.

        OP_STRING(1)=STRING(1:IEND)//' projections'
        OP_STRING(2)=BLANK(1:15)//'<elements (all/#s/GROUP#[all]>'
        OP_STRING(3)=BLANK(1:15)//'<from XI_3#[0.0]{>=0.0}>'
        OP_STRING(4)=BLANK(1:15)//'<to XI_3#[1.0]{<=1.0}>'
        OP_STRING(5)=BLANK(1:15)//'<when (static/TIME)[static]>'
        OP_STRING(6)=BLANK(1:15)//'<on (all/WS#s)[all]>'
        OP_STRING(7)=BLANK(1:15)//'<rgb=RGB[cyan]>'
        OP_STRING(8)=BLANK(1:15)//'<greater VALUE#[0.0]>'
        OP_STRING(9)=BLANK(1:15)//'<regions (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM draw data slice
C###  Description:
C###    Draw data points which lie in +/- RANGE around POSITION
C###  Parameter:    <x/y/z POSITION#>
C###    Specify the center value of either x or y or z about which
C###    data points are to be drawn.
C###  Parameter:    <range=VALUE#>
C###    Specify the width of the slice to draw data points in around
C###    the x or y or z position.
C###  Parameter:    <rgb=RGB[cyan]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.

        OP_STRING(1)=STRING(1:IEND)//' slice'
        OP_STRING(2)=BLANK(1:15)//'x/y/z POSITION#'
        OP_STRING(3)=BLANK(1:15)//'range=VALUE#'
        OP_STRING(4)=BLANK(1:15)//'<rgb=RGB[cyan]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM draw data history
C###  Description:
C###    Draw a time series history of a data point.
C###  Parameter:    <data NDHIST#[1]>
C###    Specify the data point number of which the time series will
C###    be drawn.
C###  Parameter:    <coord NJHIST#[1]>
C###    Specify the nj (global direction) number to be drawn in
C###    the history plot.
C###  Parameter:    <rgb=RGB[cyan]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.

        OP_STRING(1)=STRING(1:IEND)//' history'
        OP_STRING(2)=BLANK(1:15)//'<data NDHIST#[1]>'
        OP_STRING(3)=BLANK(1:15)//'<coord NJHIST#[1]>'
        OP_STRING(4)=BLANK(1:15)//'<rgb=RGB[cyan]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRDATA',ERROR,*9999)
      ELSE
        IF(CALL_FIT) THEN
          WRITE(OP_STRING,'('' >>WARNING: Fitting may have '
     '      //'corrupted data point coordinates'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        nx=1 !temporary cpb 22/11/94

C GMH 2/9/95 Unused        BEM   =.FALSE.
        IF(CBBREV(CO,'GEOMETRY',1,noco+1,NTCO,N3CO)) THEN
          TYPE='GEOMETRY'
          IF(CBBREV(CO,'BEM',1,noco+1,NTCO,N3CO)) THEN
C GMH 2/9/95 Unused            BEM=.TRUE.
          ENDIF
        ELSE IF(CBBREV(CO,'FIBRES',3,noco+1,NTCO,N3CO)) THEN
          TYPE='FIBRES'
          CALL ASSERT(CALC_XI,'>>Define xi positions first',ERROR,*9999)
        ELSE IF(CBBREV(CO,'FIELD',3,noco+1,NTCO,N3CO)) THEN
          TYPE='FIELD'
        ELSE IF(CBBREV(CO,'NUMBERS',1,noco+1,NTCO,N3CO)) THEN
          TYPE='NUMBERS'
        ELSE IF(CBBREV(CO,'VALUES',1,noco+1,NTCO,N3CO)) THEN
          TYPE='VALUES'
        ELSE IF(CBBREV(CO,'PROJECTIONS',2,noco+1,NTCO,N3CO)) THEN
          TYPE='PROJECTIONS'
        ELSE IF(CBBREV(CO,'SLICE',2,noco+1,NTCO,N3CO)) THEN
          TYPE='SLICE'
        ELSE IF(CBBREV(CO,'SHEETS',2,noco+1,NTCO,N3CO)) THEN
          TYPE='SHEETS'
        ELSE IF(CBBREV(CO,'SIGNAL',2,noco+1,NTCO,N3CO)) THEN
          TYPE='SIGNAL'
        ELSE IF(CBBREV(CO,'STRIPE',2,noco+1,NTCO,N3CO)) THEN
          TYPE='STRIPE'
        ELSE IF(CBBREV(CO,'TRACE',1,noco+1,NTCO,N3CO)) THEN
          TYPE='TRACE'
        ELSE IF(CBBREV(CO,'PARAMETER',2,noco+1,NTCO,N3CO)) THEN
          TYPE='PARAMETER'
        ELSE IF(CBBREV(CO,'WITH',2,noco+1,NTCO,N3CO)) THEN
          TYPE='SERIES'
        ELSE IF(CBBREV(CO,'XI',2,noco+1,NTCO,N3CO)) THEN
          TYPE='XIPOS'
          IF(CBBREV(CO,'LINEAR',1,noco+1,NTCO,N3CO)) THEN
C GMH 2/9/95 Unused            CUBIC=.FALSE.
C GMH 2/9/95 Unused            QUADRATIC=.FALSE.
C GMH 2/9/95 Unused            ORTHOG=.FALSE.
          ELSE IF(CBBREV(CO,'ORTHOGONAL',2,noco+1,NTCO,N3CO)) THEN
C GMH 2/9/95 Unused            CUBIC=.FALSE.
C GMH 2/9/95 Unused            QUADRATIC=.FALSE.
C GMH 2/9/95 Unused            ORTHOG=.TRUE.
          ELSE IF(CBBREV(CO,'CUBIC',2,noco+1,NTCO,N3CO)) THEN
C GMH 2/9/95 Unused            CUBIC=.TRUE.
C GMH 2/9/95 Unused            QUADRATIC=.FALSE.
C GMH 2/9/95 Unused            ORTHOG=.FALSE.
          ELSE IF(CBBREV(CO,'QUADRATIC',2,noco+1,NTCO,N3CO)) THEN
C GMH 2/9/95 Unused            QUADRATIC=.TRUE.
C GMH 2/9/95 Unused            CUBIC=.FALSE.
C GMH 2/9/95 Unused            ORTHOG=.FALSE.
          ELSE
C GMH 2/9/95 Unused            CUBIC=.FALSE.
C GMH 2/9/95 Unused            QUADRATIC=.FALSE.
C GMH 2/9/95 Unused            ORTHOG=.FALSE.
          ENDIF
          IF(CBBREV(CO,'NEW',1,noco+1,NTCO,N3CO)) THEN
C GMH 2/9/95 Unused            NEW=.TRUE.
          ELSE IF(CBBREV(CO,'OLD',2,noco+1,NTCO,N3CO)) THEN
C GMH 2/9/95 Unused            NEW=.FALSE.
          ELSE
C GMH 2/9/95 Unused            NEW=.TRUE.
          ENDIF
          IF(CBBREV(CO,'CENTROID',3,noco+1,NTCO,N3CO)) THEN
C GMH 2/9/95 Unused            CENTROID=.TRUE.
            IF(CBBREV(CO,'CLUSTER',3,noco+1,NTCO,N3CO)) THEN
              CALL PARSIL(CO(N3CO+1),1,NTIL,NCLUSTER,ERROR,*9999)
            ELSE
              NCLUSTER(1)=72
            ENDIF
          ELSE
C GMH 2/9/95 Unused            CENTROID=.FALSE.
          ENDIF
        ELSE IF(CBBREV(CO,'CAT',1,noco+1,NTCO,N3CO)) THEN
          TYPE='CATDAT'
          IF(CBBREV(CO,'STEP',2,noco+1,NTCO,N3CO)) THEN
            CALL PARSRL(CO(N3CO+1),1,NTIL,NO_STEP,ERROR,*9999)
            STEP=NO_STEP(1)
          ELSE
            STEP=5.0D0
          ENDIF
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' stepsize=',STEP
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ELSE IF(CBBREV(CO,'IMAGE',2,noco+1,NTCO,N3CO)) THEN
          TYPE='IMAGE'
        ELSE IF(CBBREV(CO,'BEADS',1,noco+1,NTCO,N3CO)) THEN
          TYPE='BEADS'
        ELSE IF(CBBREV(CO,'HISTORY',1,noco+1,NTCO,N3CO)) THEN
          TYPE='HISTORY'
          IF(CBBREV(CO,'DATA',2,noco+1,NTCO,N3CO)) THEN
            NO_NDHIST(1)=NDHIST
            CALL PARSIL(CO(N3CO+1),1,NTIL,NO_NDHIST,ERROR,*9999)
          ELSE
            NDHIST=1
          ENDIF
          IF(CBBREV(CO,'COORD',2,noco+1,NTCO,N3CO)) THEN
            NO_NJHIST(1)=NJHIST
            CALL PARSIL(CO(N3CO+1),1,NTIL,NO_NJHIST,ERROR,*9999)
          ELSE
            NJHIST=1
          ENDIF
          IF(DOP) THEN
            WRITE(OP_STRING,*)' history data at ndhist,njhist',
     '        ndhist,njhist
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ELSE
          TYPE='GEOMETRY'
        ENDIF

        IF(TYPE(1:5).EQ.'SLICE') THEN
          IF(CBBREV(CO,'X',1,noco+1,NTCO,N3CO)) THEN
C GMH 2/9/95 Unused            nj=1
            IWSLICE=2
          ELSE IF(CBBREV(CO,'Y',1,noco+1,NTCO,N3CO)) THEN
C GMH 2/9/95 Unused            nj=2
            IWSLICE=1
          ELSE IF(CBBREV(CO,'Z',1,noco+1,NTCO,N3CO)) THEN
C GMH 2/9/95 Unused            nj=3
            IWSLICE=3
          ENDIF
          VAL=RFROMC(CO(N3CO+1))
          IF(CBBREV(CO,'RANGE',1,noco+1,NTCO,N3CO)) THEN
            RANGE=RFROMC(CO(N3CO+1))
          ELSE
            RANGE=0.0D0
          ENDIF
          XI3MIN=VAL-RANGE
          XI3MAX=VAL+RANGE

        ELSE IF(TYPE(1:6).EQ.'STRIPE') THEN
          IF(CBBREV(CO,'GROUP',1,noco+1,NTCO,N3CO)) THEN
C GMH 2/9/95 Unused            GROUP_NAME=CO(N3CO+1)
          ELSE
            ERROR='>>>Group name needed'
            GO TO 9999
          ENDIF
          IF(CBBREV(CO,'TOTAL',1,noco+1,NTCO,N3CO)) THEN
C GMH 2/9/95 Unused            NT_STRIPE_DATA=IFROMC(CO(N3CO+1))
          ELSE
C GMH 2/9/95 Unused            NT_STRIPE_DATA=2000
          ENDIF
        ENDIF

        IF(TYPE(1:8).EQ.'GEOMETRY'.OR.TYPE(1:6).EQ.'FIELD'.OR
     '    .TYPE(1:7).EQ.'NUMBERS'.OR.TYPE(1:6).EQ.'VALUES'.OR
     '    .TYPE(1:11).EQ.'PROJECTIONS'.OR.TYPE(1:6).EQ.'FIBRES'.OR
     '    .TYPE(1:5).EQ.'TRACE'.OR.TYPE(1:9).EQ.'PARAMETER'.OR
     '    .TYPE(1:5).EQ.'XIPOS'.OR.TYPE(1:6).EQ.'CATDAT'.OR
     '    .TYPE(1:6).EQ.'SHEETS'.OR.TYPE(1:5).EQ.'BEADS') THEN
          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,
     '      ERROR,*9999)
          CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '      ERROR,*9999)
        ENDIF

        IF(TYPE(1:8).EQ.'GEOMETRY'.OR.TYPE(1:5).EQ.'FIELD') THEN
          CALL STRING_TRIM(TYPE,IBEG,IEND)
          IF(CBBREV(CO,'DEFORMED',3,noco+1,NTCO,N3CO)) THEN
            DO nonr=1,NRLIST(0)
              nr=NRLIST(nonr)
              CALL ASSERT(NJ_LOC(njl_fibr,0,nr).EQ.NJT,
     '          '>>No deformed coords',ERROR,*9999)
            ENDDO !nr
            TYPE=TYPE(IBEG:IEND)//'_DEFORMED'
          ENDIF
          IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
            IF(ABBREV(CO(N3CO+1),'COMPUTED',3)) THEN
              INDEX=0  ! Index calculated inside SGDATA
            ELSE
              INDEX=INDEX_POLYMARKER(0,'PLUS','SIZE1',CO(N3CO+1))
            ENDIF
          ELSE
            INDEX=INDEX_POLYMARKER(0,'PLUS','SIZE1','BLUE')
          ENDIF
        ELSE IF(TYPE(1:6).EQ.'FIBRES'.OR.TYPE(1:6).EQ.'SHEETS') THEN
          IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
            INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
          ELSE
            INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLUE')
          ENDIF
        ELSE IF(TYPE(1:5).EQ.'TRACE') THEN
          IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
            INDEX=INDEX_POLYMARKER(0,'CIRCLE','SIZE1',CO(N3CO+1))
          ELSE
            INDEX=INDEX_POLYMARKER(0,'CIRCLE','SIZE1','CYAN')
          ENDIF
        ELSE IF(TYPE(1:11).EQ.'PROJECTIONS') THEN
          IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
            INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
          ELSE
            INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','CYAN')
          ENDIF
          IF(CBBREV(CO,'GREATER',2,noco+1,NTCO,N3CO)) THEN
             ZDLMIN=RFROMC(CO(N3CO+1))
          ELSE
             ZDLMIN=0.0D0
          ENDIF
        ENDIF
        IF(TYPE(1:8).EQ.'GEOMETRY'.OR.TYPE(1:5).EQ.'FIELD'.OR
     '      .TYPE(1:5).EQ.'TRACE'.OR.TYPE(1:6).EQ.'FIBRES'
     '      .OR.TYPE(1:11).EQ.'PROJECTIONS') THEN
          IF(CBBREV(CO,'FROM',2,noco+1,NTCO,N3CO)) THEN
            XI3MIN=RFROMC(CO(N3CO+1))
          ELSE
            XI3MIN=0.0D0
          ENDIF
          IF(CBBREV(CO,'TO',2,noco+1,NTCO,N3CO)) THEN
            XI3MAX=RFROMC(CO(N3CO+1))
          ELSE
            XI3MAX=1.0D0
          ENDIF
          IF(CBBREV(CO,'WHEN',2,noco+1,NTCO,N3CO)) THEN
            TIME=RFROMC(CO(N3CO+1))
            STATIC=.FALSE.
          ELSE
            STATIC=.TRUE.
            TIME=0.0D0
          ENDIF
        ELSE IF(TYPE(1:6).EQ.'SHEETS') THEN
          IF(CBBREV(CO,'THETA',2,noco+1,NTCO,N3CO)) THEN
            SHEET_THETA=PI*RFROMC(CO(N3CO+1))/180.0D0
          ELSE
            SHEET_THETA=0.0D0
          ENDIF
          IF(CBBREV(CO,'RANGE',1,noco+1,NTCO,N3CO)) THEN
            SHEET_RANGE=RFROMC(CO(N3CO+1))*PI/180.0D0
          ELSE
C CPB 28/3/94 should this be = PI/180.0D0 ?????
            SHEET_RANGE=PI/18.0D0
          ENDIF
          STATIC=.TRUE.
        ENDIF

!Set default workstation numbers
        IF(TYPE(1:7).EQ.'HISTORY') THEN
          IW_DEFAULT=10
          IF(DOP) THEN
            WRITE(OP_STRING,*)' HISTORY ON WKST 10'
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ELSE IF(TYPE(1:5).EQ.'SLICE') THEN
          IW_DEFAULT=IWSLICE
        ELSE IF(TYPE(1:6).EQ.'SHEETS') THEN
          IW_DEFAULT=13
        ELSE
C MPN 3Nov97: draw on all available windows (not just 2 or 3)
C old          IW_DEFAULT=2*NJT-3+IMAP
          IW_DEFAULT=0
        ENDIF

        CALL WS_LIST(IWK,IW_DEFAULT,NTIW,noco,NTCO,CO,ERROR,*9999)
        IF(DOP) THEN
          WRITE(OP_STRING,*)' iwk(1..NTIW)=',(IWK(noiw),noiw=1,NTIW)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

        IF(ADD) THEN
          NTDATA=NTDATA+1
        ELSE IF(NTDATA.EQ.0) THEN
          NTDATA=1
        ENDIF
        CALL ASSERT(NTDATA.LE.NRM,'>>NRM too small',ERROR,*9999)
        DATYPE(NTDATA)=TYPE(1:11)
        nr=1                    !Needs fixing

        IF(.NOT.STATIC) THEN !motion
          CALL YPZP(MOTION_IY,NBH,NEELEM,NJ_LOC(NJL_GEOM,0,nr),
     '      NHP(1,nr,nx),NKH(1,1,1,nr),NPNODE,
     '      nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)
        ENDIF
        IF(TYPE(1:9).EQ.'PARAMETER') THEN
          TYPE=TYPE(1:9)//CO(noco+2)(1:2)
          IF(DOP) THEN
            WRITE(OP_STRING,*) 'DATYPE(ntdata)=',DATYPE(NTDATA)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
        IF(TYPE(1:9).EQ.'PARAMETER'.OR.TYPE(1:5).EQ.'FIBRE') THEN
          IF(CBBREV(CO,'SCALE',1,noco+1,NTCO,N3CO)) THEN
            SCALE=RFROMC(CO(N3CO+1))
          ELSE
            SCALE=1.0D0
          ENDIF
        ENDIF
        IF(TYPE(1:8).EQ.'GEOMETRY'.OR.TYPE(1:5).EQ.'FIELD'.OR
     '    .TYPE(1:7).EQ.'NUMBERS'.OR.TYPE(1:6).EQ.'VALUES'.OR
     '    .TYPE(1:11).EQ.'PROJECTIONS'.OR
     '    .TYPE(1:6).EQ.'FIBRES'.OR.TYPE(1:6).EQ.'SHEETS'.OR
     '    .TYPE(1:5).EQ.'TRACE'.OR.TYPE(1:9).EQ.'PARAMETER'.OR
     '    .TYPE(1:7).EQ.'HISTORY'.OR.TYPE(1:5).EQ.'SLICE') THEN
          DO noiw=1,NTIW
            IW=IWK(noiw)
            CALL ACWK(iw,1,ERROR,*9999)
            CALL SGDATA(INDEX,IBT,IDO,INP,ISDANO,ISDAPR,
     '        ISDATA(iw,NTDATA),ISDATR,ISEG,iw,LD,MXI,
     '        NBJ,NBH,NDDL,NDLT,NDP,
     '        NEELEM,NELIST,NKHE,NKJE,NPF,NPNE,NRE,NVHE,NVJE,
     '        NW(1,1,nx),CE(1,1,nx),CG,CGE(1,1,1,nx),
     '        CP(1,1,nx),CSEG,CURVCORRECT,PG,SE,STATIC,
     '        TITLE,TYPE,WD,WDL,
     '        XA,XE,XG,XID,XIDL,XP,ZA,ZD,ZDD,ZDL,ZDLMIN,ZE,ZP,ERROR,
     '        *9999)
            CALL DAWK(iw,1,ERROR,*9999)
          ENDDO
        ELSE IF(TYPE(1:6).EQ.'SERIES') THEN
          IW=2*NJT-3
          CALL ACWK(iw,1,ERROR,*9999)
          CALL SGDATA(INDEX,IBT,IDO,INP,ISDANO,ISDAPR,ISDATA(iw,NTDATA),
     '      ISDATR,ISEG,iw,LN,MXI,NBJ,NBH,NDDL,NDLT,NDP,
     '      NEELEM,NELIST,NKHE,NKJE,NPF,NPNE,NRE,NVHE,NVJE,NW(1,1,nx),
     '      CE(1,1,nx),CG,CGE(1,1,1,nx),
     '      CP(1,1,nx),CSEG,CURVCORRECT,PG,SE,STATIC,
     '      TITLE,TYPE,WD,WDL,
     '      XA,XE,XG,XID,XIDL,XP,ZA,ZD,ZDD,ZDL,ZDLMIN,ZE,ZP,ERROR,*9999)
          CALL DAWK(iw,1,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('DRDATA')
      RETURN
 9999 CALL ERRORS('DRDATA',ERROR)
      CALL EXITS('DRDATA')
      RETURN 1
      END


