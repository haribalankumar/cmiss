C     CMISS Module FE14: Routines which call GL graphics via C

C     Function INDEX_FILL_AREA      returns bundle index for fill-areas
C     Function INDEX_POLYLINE       returns bundle index for polylines
C     Function INDEX_POLYMARKER     returns bundle index for polymarkers
C     Function INDEX_TEXT           returns bundle index for text
C     Subroutine ACWK               activates workstation
C     Subroutine ANIMATE            set up or perform anamation
C     Subroutine ARCHIVE_GRAPHICS   archives PHIGS structure or GKS segments
C     Subroutine CELL_ARRAY         draws bit-image cell array
C     Subroutine CHOICE             sets up choice window
C     Subroutine CIRCLE             draws circle
C     Subroutine CLOSE_PRINT_FILE   deactivates print ws
C     Subroutine CLOSE_SEGMENT      closes graphics segment
C     Subroutine CLOSE_WS           closes graphics workstation
C     Subroutine CLWS               closes workstations
C     Subroutine COLOUR             returns colour index
C     Subroutine COWK               closes and reopens a workstation
C     Subroutine CREATE_SEGMENT     create graphics segment without ISEG,CSEG
C     Subroutine DAWK               deactivates a workstation
C     Subroutine DBOX               draws box
C     Subroutine DELETE_SEGMENT     delete graphics segment
C     Subroutine DETECT             change segment detectability
C     Subroutine DETRAN             Defines 3D transformations using Phigs
C     Subroutine DEVIEW             defines 3D views using Phigs
C     Subroutine DISPLAY            displays string in a window
C     Subroutine DISPLAY_FILE       finds files in current directory
C     Subroutine DOCUM              displays documentation window
C     Subroutine ELLIPSE            draws fill-area ellipse
C     Subroutine EVENT              handles event mode i/p
C     Subroutine FILL_AREA          draws fill-area
C     Subroutine FLUSH_DEVICE_EVENTS flush queues
C     Subroutine FREEWK             returns first available workstation number
C     Subroutine GKS_DRAW           handles GKS draw functions
C     Subroutine GKS_STRG           GKS string input
C     Subroutine GKSTEXT            text entered by GKS string input
C     Subroutine INPUT_MODE         set input device mode (request,event,sample)
C     Subroutine LINE3D             draws line on 3D viewport
C     Subroutine LICOLREP           inquires phigs colour representations
C     Subroutine LISTRU             inquires phigs data structures
C     Subroutine LOCATOR            calls GKS locator
C     Subroutine OPEN_PRINT_FILE    activates print ws and sets indices, etc
C     Subroutine OPEN_SEGMENT       opens  graphics segment
C     Subroutine PHIG               defines PHIGS initial transf.s and view mat.
C     Subroutine PICK               use pick input device
C     Subroutine POLYLINE           draws polyline
C     Subroutine POLYMARKER         draws polymarker
C     Subroutine PRECHOICE1         sets up choice window for special case
C     Subroutine PRECHOICE2         sets up choice window for special case
C     Subroutine PRINT_IMAGE_FILL_AREAS   not used?
C     Subroutine QUIT_GRAPHICS      close any workstations and gks or phigs
C     Subroutine RECALL_GRAPHICS    recalls PHIGS structure or GKS segments
C     Subroutine ROTATE             rotate a window
C     Subroutine SET_COLOUR_LUT     set colour lookup table
C     Subroutine SET_COLOUR_LUT_RANGE set colour lookup table for a given range
C     Subroutine SET_COLOUR_ONE     set RGB colour rep for one LUT index
C     Subroutine SET_COLOUR_REP     set RGB colour rep for all LUT indices
C     Subroutine SET_FILL_REP       resets fill area representation
C     Subroutine SET_FILL_AREA_REP  set polyline representation
C     Subroutine SET_GLOBAL_XFORM   Updates global object transformation
C     Subroutine SET_POLYLINE_REP   set polyline representation
C     Subroutine SET_POLYMARKER_REP set polyline representation
C     Subroutine SET_TEXT_REP       set polyline representation
C     Subroutine SET_VIEW           Updates viewing transform
C     Subroutine SETUP              performs setup operations for workstations
C     Subroutine STROKE             calls GKS stroke
C     Subroutine SURFACE            draws shaded surface
C     Subroutine TEXT               draws text
C     Subroutine VALUATOR           calls GKS valuator
C     Subroutine VISIB              change segment visibility
C     Subroutine WKST_WINDOW        change workstation window


      FUNCTION INDEX_FILL_AREA(ICOLOUR,PATTERN_FILL,STYLE_FILL,RGB_FILL)

C**** Returns fill-area bundle index:
C**** If ICOLOUR=0, for given pattern and style, as follows:
C****   Fill-area Index = 1
C**** Else if 1<=ICOLOUR<=233 for given ICOLOUR as follows:
C****   Fill-area Index = 17-249

      CHARACTER PATTERN_FILL*(*),STYLE_FILL*(*),RGB_FILL*(*),
     '  CUPPER*8,PATTERN*8,STYLE*8,RGB*8
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'


      IF(ICOLOUR.EQ.0) THEN
        PATTERN=CUPPER(PATTERN_FILL)
        STYLE  =CUPPER(STYLE_FILL)
        RGB    =CUPPER(RGB_FILL)
        IF(     RGB(1:5).EQ.'BLACK' .OR.RGB(1:3).EQ.'000') THEN !black
        ELSE IF(RGB(1:3).EQ.'RED'   .OR.RGB(1:3).EQ.'100') THEN !red
        ELSE IF(RGB(1:5).EQ.'GREEN' .OR.RGB(1:3).EQ.'010') THEN !green
        ELSE IF(RGB(1:6).EQ.'YELLOW'.OR.RGB(1:3).EQ.'110') THEN !yellow
        ELSE IF(RGB(1:4).EQ.'BLUE'  .OR.RGB(1:3).EQ.'001') THEN !blue
        ELSE IF(RGB(1:4).EQ.'CYAN'  .OR.RGB(1:3).EQ.'011') THEN !cyan
        ELSE IF(RGB(1:5).EQ.'WHITE' .OR.RGB(1:3).EQ.'111') THEN !white
          INDEX=0
        ENDIF
        IF(DOP) THEN
          WRITE(IOOP,'('' Index_Fill_Area: PATTERN='',A,'' STYLE='',A,
     '      '' RGB='',I3,'' INDEX='',I3)') PATTERN,STYLE,RGB,INDEX
        ENDIF

      ELSE IF(ICOLOUR.GE.1.AND.ICOLOUR.LE.233) THEN

        INDEX=ICOLOUR+16

        IF(DOP) THEN
          WRITE(IOOP,'('' Index_Fill_Area: ICOLOUR='',I3,'' INDEX='','
     '      //'I3)') ICOLOUR,INDEX
        ENDIF

      ENDIF

      INDEX_FILL_AREA=INDEX

      RETURN
      END


      FUNCTION INDEX_POLYLINE(ICOLOUR,TYPE_LINE,WIDTH_LINE,RGB_LINE)

C**** Returns polyline bundle index:
C**** If ICOLOUR=0, for given line type, width and RGB value, as follows:
C****   Polyline Index = 1       black solid    width1
C****                    2         "   dotted      "
C****                    3         "   dashed      "
C****                    4         "   dot-dash    "
C****                    5-8         as above   width2
C****
C****                    9 -16   red    as above
C****                    17-24   green      "
C****                    25-32   blue       "
C****                    33-40   cyan       "
C****
C**** Else if 1<=ICOLOUR<=216 for given ICOLOUR as follows:
C****   Polyline Index = 41-256  216 colours solid width2

      CHARACTER TYPE_LINE*(*),WIDTH_LINE*(*),RGB_LINE*(*),
     '  CUPPER*8,TYPE*8,WIDTH*8,RGB*8
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'


      IF(ICOLOUR.EQ.0) THEN
        TYPE =CUPPER(TYPE_LINE)
        WIDTH=CUPPER(WIDTH_LINE)
        RGB  =CUPPER(RGB_LINE)

        IF(     RGB(1:5).EQ.'BLACK'.OR.RGB(1:3).EQ.'000') THEN
          INDEX_OFFSET_1=0
        ELSE IF(RGB(1:3).EQ.'RED'  .OR.RGB(1:3).EQ.'100') THEN
          INDEX_OFFSET_1=8
        ELSE IF(RGB(1:5).EQ.'GREEN'.OR.RGB(1:3).EQ.'010') THEN
          INDEX_OFFSET_1=16
        ELSE IF(RGB(1:4).EQ.'BLUE' .OR.RGB(1:3).EQ.'001') THEN
          INDEX_OFFSET_1=24
        ELSE IF(RGB(1:4).EQ.'CYAN' .OR.RGB(1:3).EQ.'011') THEN
          INDEX_OFFSET_1=32
        ELSE
          WRITE(IOOP,'(''>>Colour not defined'')')
        ENDIF

        IF(WIDTH(1:6).EQ.'WIDTH1') THEN
          INDEX_OFFSET_2=0
        ELSE IF(WIDTH(1:6).EQ.'WIDTH2') THEN
          INDEX_OFFSET_2=4
        ENDIF
        IF(TYPE(1:5).EQ.'SOLID') THEN
          INDEX=1
        ELSE IF(TYPE(1:6).EQ.'DOTTED') THEN
          INDEX=2
        ELSE IF(TYPE(1:6).EQ.'DASHED') THEN
          INDEX=3
        ELSE IF(TYPE(1:8).EQ.'DOT-DASH') THEN
          INDEX=4
        ENDIF

        INDEX=INDEX+INDEX_OFFSET_1+INDEX_OFFSET_2

        IF(DOP) THEN
          WRITE(IOOP,'('' Index_Polyline: Type='',A,'' Width='',A,
     '    '' RGB='',A,'' INDEX='',I3)') TYPE,WIDTH,RGB,INDEX
        ENDIF

      ELSE IF(ICOLOUR.GE.1.AND.ICOLOUR.LE.216) THEN

        INDEX=ICOLOUR+40

        IF(DOP) THEN
          WRITE(IOOP,'('' Index_Polyline: ICOLOUR='',I3,'' INDEX='','
     '      //'I3)') ICOLOUR,INDEX
        ENDIF

      ENDIF

      INDEX_POLYLINE=INDEX

      RETURN
      END


      FUNCTION INDEX_POLYMARKER(ICOLOUR,TYPE_MARKER,SIZE_MARKER,
     '  RGB_MARKER)

C**** Returns polymarker bundle index for given marker type, size & RGB value.
C**** If ICOLOUR=0  for given marker type, size and RGB value, as follows:
C****   Polymarker Index = 1       black plus      size1
C****                      2         "   asterisk    "
C****                      3         "   circle      "
C****                      4         "   point       "
C****                      5-8         as above    size2
C****
C****                      9 -16   red    as above
C****                      17-24   green      "
C****                      25-32   blue       "
C****                      33-40   cyan       "
C****
C**** Else if 1<=ICOLOUR<=216 for given ICOLOUR as follows:
C****   Polymarker Index = 41-256  216 colours plus size2

      CHARACTER TYPE_MARKER*(*),SIZE_MARKER*(*),RGB_MARKER*(*),
     '  CUPPER*8,TYPE*8,SIZE*8,RGB*8
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'


      IF(ICOLOUR.EQ.0) THEN
        TYPE=CUPPER(TYPE_MARKER)
        SIZE=CUPPER(SIZE_MARKER)
        RGB =CUPPER(RGB_MARKER)

        IF(     RGB(1:5).EQ.'BLACK'.OR.RGB(1:3).EQ.'000') THEN
          INDEX_OFFSET_1=0
        ELSE IF(RGB(1:3).EQ.'RED'  .OR.RGB(1:3).EQ.'100') THEN
          INDEX_OFFSET_1=8
        ELSE IF(RGB(1:5).EQ.'GREEN'.OR.RGB(1:3).EQ.'010') THEN
          INDEX_OFFSET_1=16
        ELSE IF(RGB(1:4).EQ.'BLUE' .OR.RGB(1:3).EQ.'001') THEN
          INDEX_OFFSET_1=24
        ELSE IF(RGB(1:4).EQ.'CYAN' .OR.RGB(1:3).EQ.'011') THEN
          INDEX_OFFSET_1=32
        ELSE
          WRITE(IOOP,'(''>>Colour not defined'')')
        ENDIF

        IF(SIZE(1:5).EQ.'SIZE1') THEN
          INDEX_OFFSET_2=0
        ELSE IF(SIZE(1:5).EQ.'SIZE2') THEN
          INDEX_OFFSET_2=4
        ENDIF
        IF(TYPE(1:4).EQ.'PLUS') THEN
          INDEX=1
        ELSE IF(TYPE(1:8).EQ.'ASTERISK') THEN
          INDEX=2
        ELSE IF(TYPE(1:6).EQ.'CIRCLE') THEN
          INDEX=3
        ELSE IF(TYPE(1:5).EQ.'POINT') THEN
          INDEX=4
        ENDIF

        INDEX=INDEX+INDEX_OFFSET_1+INDEX_OFFSET_2

        IF(DOP) THEN
          WRITE(IOOP,'('' Index_Polymarker: Type='',A,'' Width='',A,
     '    '' RGB='',A,'' INDEX='',I3)') TYPE,WIDTH,RGB,INDEX
        ENDIF

      ELSE IF(ICOLOUR.GE.1.AND.ICOLOUR.LE.216) THEN

        INDEX=ICOLOUR+40

        IF(DOP) THEN
          WRITE(IOOP,'('' Index_Polymarker: ICOLOUR='',I3,'' INDEX='','
     '      //'I3)') ICOLOUR,INDEX
        ENDIF

      ENDIF

      INDEX_POLYMARKER=INDEX

      RETURN
      END


      FUNCTION INDEX_TEXT(ICOLOUR,WIDTH_TEXT,FONT_TEXT,RGB_TEXT)

C**** Returns text bundle index for given text size,font and RGB value.
C**** For given text size and RGB value, as follows:
C****         Text Index = 1       black font1    width1
C****                      2         "   font2      "
C****                      3         "   font3      "
C****                      4         "   font4      "
C****                      5-8         as above for width2
C****
C****                      9 -16   red    as above
C****                      17-24   green      "
C****                      25-32   blue       "
C****                      33-40   cyan       "

      CHARACTER WIDTH_TEXT*(*),RGB_TEXT*(*),FONT_TEXT*(*),CUPPER*8,
     '  FONT*8,WIDTH*8,RGB*8
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'


      IF(ICOLOUR.EQ.0) THEN
        FONT=CUPPER(FONT_TEXT)
        WIDTH=CUPPER(WIDTH_TEXT)
        RGB =CUPPER(RGB_TEXT)
        IF(     RGB(1:5).EQ.'BLACK'.OR.RGB(1:3).EQ.'000') THEN
          INDEX_OFFSET_1=0
        ELSE IF(RGB(1:3).EQ.'RED'  .OR.RGB(1:3).EQ.'100') THEN
          INDEX_OFFSET_1=8
        ELSE IF(RGB(1:5).EQ.'GREEN'.OR.RGB(1:3).EQ.'010') THEN
          INDEX_OFFSET_1=16
        ELSE IF(RGB(1:4).EQ.'BLUE' .OR.RGB(1:3).EQ.'001') THEN
          INDEX_OFFSET_1=24
        ELSE IF(RGB(1:4).EQ.'CYAN' .OR.RGB(1:3).EQ.'011') THEN
          INDEX_OFFSET_1=32
        ELSE
          WRITE(IOOP,'(''>>Colour not defined'')')
        ENDIF

        IF(WIDTH(1:6).EQ.'WIDTH1') THEN
          INDEX_OFFSET_2=0
        ELSE IF(WIDTH(1:6).EQ.'WIDTH2') THEN
          INDEX_OFFSET_2=4
        ENDIF
        IF(FONT(1:5).EQ.'FONT1') THEN
          INDEX=1
        ELSE IF(FONT(1:5).EQ.'FONT2') THEN
          INDEX=2
        ELSE IF(FONT(1:5).EQ.'FONT3') THEN
          INDEX=3
        ELSE IF(FONT(1:5).EQ.'FONT4') THEN
          INDEX=4
        ENDIF

        INDEX=INDEX+INDEX_OFFSET_1+INDEX_OFFSET_2

        IF(DOP) THEN
          WRITE(IOOP,'('' Index_Text: Width='',A,'' RGB='',A,'
     '      //''' INDEX='',I3)') WIDTH,RGB,INDEX
        ENDIF
      ENDIF
      INDEX_TEXT=INDEX

      RETURN
      END


      SUBROUTINE ACWK(IW,ID,ERROR,*)

C**** Activate workstation identified by IW.
C**** If ID=1 output is deferred.
C**** If IW=5 or 6 buffers are cleared.

      CHARACTER ERROR*(*),OPTION(25)*15,CO(1)*1,STRING*1
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:post00.cmn'

      CALL ENTERS('ACWK',*9999)
      IF(DOP) WRITE(IOOP,'('' IW='',I3,'' IWKS(iw)='',I2)') IW,IWKS(IW)
      IF(IWKS(IW).EQ.0) THEN
        CALL SETUP(IW,ERROR,*9999)

C ***   Display choice windows
        IF(NJT.LE.2) THEN
          IF(IW.EQ.1) THEN
            OPTION( 1)='Archive..'
            OPTION( 2)='Cancel'
            OPTION( 3)='Change..'
            OPTION( 4)='Define..'
            OPTION( 5)='Display..'
            OPTION( 6)='Hide/Show..'
            OPTION( 7)='Label..'
            OPTION( 8)='List..'
            OPTION( 9)='Pick..'
            OPTION(10)='Print..'
            OPTION(11)='Read..'
            OPTION(12)='Recall..'
            OPTION(13)='Refine..'
            OPTION(14)='Transform..'
            OPTION(15)='...'
            OPTION(16)='Exit'
            NTCH=16
            CALL CHOICE('DETRAN',1,1,INSTAT,91,'EVENT',1,NTCH,NOCH,
     '        NOCO,9,CO,OPTION,STRING,0.0,YDISP-0.50*XDISP,ERROR,*9999)
          ENDIF
        ELSE IF(NJT.EQ.3) THEN
          IF(IW.EQ.1) THEN
            OPTION( 1)='Archive..'
            OPTION( 2)='Cancel'
            OPTION( 3)='Change..'
            OPTION( 4)='Define..'
            OPTION( 5)='Display..'
            OPTION( 6)='Hide/Show..'
            OPTION( 7)='Label..'
            OPTION( 8)='List..'
            OPTION( 9)='Pick..'
            OPTION(10)='Print..'
            OPTION(11)='Read..'
            OPTION(12)='Recall..'
            OPTION(13)='Refine..'
            OPTION(14)='Transform..'
            OPTION(15)='...'
            OPTION(16)='Exit'
            NTCH=16
            CALL CHOICE('DETRAN',1,1,INSTAT,91,'EVENT',1,NTCH,NOCH,
     '        NOCO,9,CO,OPTION,STRING,0.0,0.50*DISP,ERROR,*9999)
          ELSE IF(IW.EQ.2) THEN
            OPTION( 1)='Archive..'
            OPTION( 2)='Cancel'
            OPTION( 3)='Change..'
            OPTION( 4)='Define..'
            OPTION( 5)='Hide/Show..'
            OPTION( 6)='Label..'
            OPTION( 7)='Pick..'
            OPTION( 8)='Print..'
            OPTION( 9)='Recall..'
            OPTION(10)='Transform..'
            OPTION(11)='Exit'
            NTCH=11
            CALL CHOICE('DETRAN',1,1,INSTAT,92,'EVENT',1,NTCH,NOCH,
     '        NOCO,9,CO,OPTION,STRING,0.99*DISP,0.99*DISP,ERROR,*9999)
          ELSE IF(IW.EQ.3) THEN
            OPTION( 1)='Archive..'
            OPTION( 2)='Back plane'
            OPTION( 3)='Cancel'
            OPTION( 4)='Define(s)..'
            OPTION( 5)='Front plane'
            OPTION( 6)='Hide/Show..'
            OPTION( 7)='Label..'
            OPTION( 8)='Pan'
            OPTION( 9)='Parallel'
            OPTION(10)='Perspective'
            OPTION(11)='Print..'
            OPTION(12)='Proj ref pt'
            OPTION(13)='Recall..'
            OPTION(14)='Rescale'
            OPTION(15)='Reset'
            OPTION(16)='Rotate data'
            OPTION(17)='Rotate view'
            OPTION(18)='Save view'
            OPTION(19)='Select view'
            OPTION(20)='View ref pt'
            OPTION(21)='View plane'
            OPTION(22)='View dist'
            OPTION(23)='View up'
            OPTION(24)='Zoom'
            OPTION(25)='Exit'
            NTCH=25
            CALL CHOICE('DETRAN',1,1,INSTAT,93,'EVENT',1,NTCH,NOCH,
     '        NOCO,9,CO,OPTION,STRING,0.99*DISP,0.48*DISP,ERROR,*9999)
          ELSE IF(IW.EQ.4) THEN
            OPTION( 1)='Archive..'
            OPTION( 2)='Cancel'
            OPTION( 3)='Define..'
            OPTION( 4)='Hide/Show..'
            OPTION( 5)='Label..'
            OPTION( 6)='Pick..'
            OPTION( 7)='Print..'
            OPTION( 8)='Recall..'
            OPTION( 9)='Transform..'
            OPTION(10)='Exit'
            NTCH=10
            CALL CHOICE('DETRAN',1,1,INSTAT,94,'EVENT',1,NTCH,NOCH,
     '        NOCO,9,CO,OPTION,STRING,0.10*DISP,0.25*DISP,ERROR,*9999)
          ENDIF
        ENDIF
        IF(IW.EQ.35) THEN !signal trace window
          OPTION( 1)='Accept'
          OPTION( 2)='Reject'
          OPTION( 3)='Modify'
          OPTION( 4)='Choose Signal'
          OPTION( 5)='Exit'
          NTCH=5
          CALL CHOICE('DITRAC',1,1,INSTAT,36,'EVENT',1,NTCH,NOCH,NOCO,
     '      10,CO,OPTION,STRING,0.75*XDISP,YDISP-0.51*XDISP,ERROR,*9999)
        ENDIF
      ENDIF

      IF(IWKS(IW).EQ.1) THEN
        IF(IW.NE.5.AND.IW.NE.6)THEN
          CALL GL_ACWK(IW)           !activate this window
          IWKS(IW)=2
        ENDIF
      ENDIF

      CALL EXITS('ACWK')
      RETURN
 9999 CALL ERRORS('ACWK',ERROR)
      CALL EXITS('ACWK')
      RETURN 1
      END


!news AAY 16 Jun 91 put animate here since it calls gl routine
      SUBROUTINE ANIMATE(ISEG,NOCO,NTCO,NTCOQU,CO,COQU,CSEG,STRING,
     '  ERROR,*)

C     rotate window

      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
      INTEGER ISEG(*),NTCOQU(*),TIMES(60)
      CHARACTER CO(*)*(*),COQU(NXCO,*)*(*),CSEG(*)*(*),ERROR*(*),
     '  STRING*(MXCH)
      LOGICAL CBBREV,ABBREV

      CALL ENTERS('ANIMATE',*9999)
      IF(CO(NOCO+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        WRITE(IOOP,'(X,A,/)') STRING(1:IEND)
     '    //' loop times TIMELIST[1..NFRAME] on IW[1] delay DELAY[1]',
     '    ' step times TIMELIST[1..NFRAME] on IW[1]',
     '    ' start/end',
     '    ' time=TIME',
     '    ' display TIME'
      ELSE
        IF(CBBREV(CO,'LOOP',2,NOCO+1,NTCO,N3CO).OR.
     '     CBBREV(CO,'STEP',2,NOCO+1,NTCO,N3CO)) THEN
    IF(CBBREV(CO,'LOOP',2,NOCO+1,NTCO,N3CO))THEN
      LOOP=1
          ELSE
      LOOP=0
          ENDIF
          IF(CBBREV(CO,'TIMES',2,NOCO+1,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),60,NTIMES,TIMES,ERROR,*9999)
          ELSE
            NTIMES=1
            TIMES(1)=1
          ENDIF
          IF(CBBREV(CO,'ON',2,NOCO+1,NTCO,N3CO)) THEN
            IW=IFROMC(CO(N3CO+1))
          ELSE
            IW=1
          ENDIF
          IF(CBBREV(CO,'DELAY',2,NOCO+1,NTCO,N3CO)) THEN
            IDELAY=IFROMC(CO(N3CO+1))
          ELSE
            IDELAY=1
          ENDIF
          CALL ACWK(IW,0,ERROR,*9999)
          CALL GL_ANIMATE(IW,NTIMES,TIMES,IDELAY,LOOP)
          CALL DAWK(IW,0,ERROR,*9999)
        ELSE IF(CBBREV(CO,'START',2,NOCO+1,NTCO,N3CO)) THEN
    CALL GL_SET_ANIMATION_FLAG(1)
        ELSE IF(CBBREV(CO,'END',2,NOCO+1,NTCO,N3CO)) THEN
    CALL GL_SET_ANIMATION_FLAG(0)
        ELSE IF(CBBREV(CO,'TIME',2,NOCO+1,NTCO,N3CO)) THEN
          ITIME=IFROMC(CO(N3CO+1))
    CALL GL_SET_ANIMATION_FRAME(ITIME)
        ELSE IF(CBBREV(CO,'DISPLAY',2,NOCO+1,NTCO,N3CO)) THEN
          ITIME=IFROMC(CO(N3CO+1))
          CALL ACWK(IW,0,ERROR,*9999)
    CALL GL_DRAW_STRUCT_AT_TIME(IW,ITIME)
          CALL DAWK(IW,0,ERROR,*9999)
        ENDIF
      ENDIF
      CALL EXITS('ANIMATE')
      RETURN
 9999 CALL ERRORS('ANIMATE',ERROR)
      CALL EXITS('ANIMATE')
      RETURN 1
      END
!newe


      SUBROUTINE ARCHIVE_GRAPHICS(IW,FILE_NAME,ERROR,*)

C**** Archives structure.

      CHARACTER FILE_NAME*(*),ERROR*(*)

      CALL ENTERS('ARCHIVE_GRAPHICS',*9999)
C     not supported yet

      CALL EXITS('ARCHIVE_GRAPHICS')
      RETURN
 9999 CALL ERRORS('ARCHIVE_GRAPHICS',ERROR)
      CALL EXITS('ARCHIVE_GRAPHICS')
      RETURN 1
      END


!news gl Bezier function not supported yet
      SUBROUTINE BEZIER(INDEX,INDEX_PLIN,ISBEZE,ISEG,IW,NTL,
     '  XPTS,YPTS,CSEG,ERROR,*)

C**** Creates Bezier control points on workstation IW.
C**** INDEX is the GKS line type index
C**** INDEX_PLIN is the polyline index
C**** XBEZ(I),YBEZ(I),I=1 & 4 are x,y coords of nodes.
C**** XBEZ(I),YBEZ(I),I=2 & 3 are x,y coords of slope control points.
C**** Note: This routine creates temporary segments which are deleted at
C**** the end and the total segment count NTSG is reduced again.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'cmiss$reference:ntsg00.cmn'
!     Parameter List
      INTEGER INDEX,INDEX_PLIN,ISBEZE,ISEG(*),IW,NTL
      REAL XPTS(*),YPTS(*)
      CHARACTER CSEG(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER NTLM
      PARAMETER(NTLM=20)
      INTEGER I,IFROMC,INDEX_OLD,INSTAT,IPICK,ISEGM,ISL2BE(NTLM),
     '  ISL3BE(NTLM),ISN2BE(NTLM),ISN3BE(NTLM),LD1,N2LI,NL
      REAL BEZLEN(NTLM),PL(3,21),PT(3,2),XBEZ(4,NTLM),XSLOPE,XWC,
     '  YBEZ(4,NTLM),YSLOPE,YWC
      CHARACTER BLANK*1,CFROMI*4
      LOGICAL AUX
      DATA LD1/1/,BLANK/' '/

      CALL ENTERS('BEZIER',*9999)
      CALL EXITS('BEZIER')
      RETURN
 9999 CALL ERRORS('BEZIER',ERROR)
      CALL EXITS('BEZIER')
      RETURN 1
      END
!newe


!news
      SUBROUTINE BEZIER_POINTS(PL,XBEZ,YBEZ,ERROR,*)

C**** Calculates 21 points along Bezier curve given by XBEZ(i),YBEZ(i),i=1,4.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
!     Parameter List
      REAL PL(3,*),XBEZ(*),YBEZ(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J,K
      REAL XI,XI2,XI3

      CALL ENTERS('BEZIER_POINTS',*9999)
      IF(DOP) THEN
        WRITE(*,'('' XBEZ: '',4E12.3)') (XBEZ(I),I=1,4)
        WRITE(*,'('' YBEZ: '',4E12.3)') (YBEZ(I),I=1,4)
      ENDIF
      DO J=0,20
        XI=REAL(J)/20.0
        XI2=XI*XI
        XI3=XI2*XI
        PL(1,J+1)=(1.-3.*XI+3.*XI2-XI3) * XBEZ(1) +
     '            (3.*XI-6.*XI2+3.*XI3) * XBEZ(2) +
     '            (3.*XI2-3.*XI3      ) * XBEZ(3) +
     '            (XI3                ) * XBEZ(4)
        PL(2,J+1)=(1.-3.*XI+3.*XI2-XI3) * YBEZ(1) +
     '            (3.*XI-6.*XI2+3.*XI3) * YBEZ(2) +
     '            (3.*XI2-3.*XI3      ) * YBEZ(3) +
     '            (XI3                ) * YBEZ(4)
        IF(DOP) THEN
          WRITE(IOOP,'('' J,XI= '',I3,F12.5)') J,XI
          WRITE(IOOP,'('' PL 1,2 = '',2F12.4)') (PL(K,J+1),K=1,2)
        ENDIF
      ENDDO

      CALL EXITS('BEZIER_POINTS')
      RETURN
 9999 CALL ERRORS('BEZIER_POINTS',ERROR)
      CALL EXITS('BEZIER_POINTS')
      RETURN 1
      END
!newe


!new  BUILD_XFORM_MATRIX3 left out


      SUBROUTINE CELL_ARRAY(IW,NDIM,XMIN_CA,XMAX_CA,YMIN_CA,YMAX_CA,
     '  ERROR,*)

C**** Draws a cell array given by array I2P(*,*,IW) (512x512) on w/s IW(1/2).
C**** If NDIM is not equal to 512, the program only draws every 512/NDIM
C**** pixel.  This requires an even division i.e. NDIM=256,128,....

      CHARACTER ERROR*(*)
      INCLUDE 'cmiss$reference:view00.cmn'
      INCLUDE 'cmiss$reference:pics00.cmn'

      CALL ENTERS('CELL_ARRAY',*9999)
!old  NSTEP=NINT(512.0/REAL(NDIM))
!old  DO I=1,NDIM
!old    DO J=1,NDIM
!old      I2P(I,J,IW+IOFFSET)=255-I2P(I*NSTEP,J*NSTEP,IW)
!old    ENDDO
!old  ENDDO
!old  IF(IW.EQ.1.OR.IW.EQ.2) THEN
        CALL GL_CELL_ARRAY(XMIN_CA,YMIN_CA,XMAX_CA,YMAX_CA,NXM,NXM,
     '    NDIM,NDIM,I2P(1,1,IW))
!old  ENDIF

      CALL EXITS('CELL_ARRAY')
      RETURN
 9999 CALL ERRORS('CELL_ARRAY',ERROR)
      CALL EXITS('CELL_ARRAY')
      RETURN 1
      END


      SUBROUTINE CHOICE(LABEL,IDEVICE,INIT,INSTAT,IW,MODE,INCH,NTCH,
     '  NOCH,NOCO,NTYPE,CO,OPTION,STRING,XREF,YREF,ERROR,*)

C**** Prompts the user with NTCH choices listed in the array of strings OPTION
C**** INIT specifies whether choice is to be initialised (1:yes,2:no)
C**** IDEVICE indicates input device: 1 for normal 4 for mouse buttons
C**** MODE can be 'REQUEST' or 'EVENT'
C**** The list of choices is headed with the string LABEL
C**** INCH=the initial choice
C**** NOCH=the returned choice
C**** Up to 55 (MAXCH) choices are permitted.
C**** XREF,YREF are coordinates of reference point
C**** NTYPE=0 is for command list when type ?
C****         (reference point is bottom RH corner of screen)
C**** NTYPE=1  is for reference point at bottom LH corner (screen coords)
C**** NTYPE=2  is  "      "       "    " bottom RH corner    "      "
C**** NTYPE=3  is  "      "       "    " top    LH corner (world  coords)
C**** NTYPE=4  is  "      "       "    " top    RH corner    "      "
C**** NTYPE=5  is  "      "       "    " bottom LH corner (screen coords)
C**** NTYPE=6  is  "      "       "    " bottom RH corner    "      "
C**** NTYPE=7  is  "      "       "    " top    LH corner (screen coords)
C**** NTYPE=8  is  "      "       "    " top    RH corner    "      "
C**** NTYPE=9  is  "      "       "    " top    LH corner (screen coords)
C**** NTYPE=10 is  "      "       "    " top    RH corner    "      "
C****   (for NTYPE=5,6,7,8 XREF & YREF are factors multiplying
C****    XDISP & YDISP so that XDISP,YDISP need not be precalculated)
C**** LNCHMX is maximum number of characters in option list
C**** CHSZ   is character height
C**** CHLN   is character width

      INTEGER LNCH(55)
      CHARACTER CLASS*6,CO(*)*(*),ERROR*(*),OPTION(*)*(*),MODE*(*),
     '  STRING*(*),DATA(55)*80,LABEL*(*)
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:disp00.cmn'
      INCLUDE 'cmiss$reference:echo00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:gks000.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'

      DATA LD1/1/,LNDATA/50/,MAXCH/55/

      CALL ENTERS('CHOICE',*9999)
      IF(INIT.EQ.1) THEN !initialise choice
        CALL SETUP(IW,ERROR,*9999)
        LNCHMX=1
        DO N1CH=1,NTCH
          CALL STRING_TRIM(OPTION(N1CH),IBEG,IEND)
          LNCH(N1CH)=IEND-IBEG+1
          IF(LNCH(N1CH).GT.LNCHMX) LNCHMX=LNCH(N1CH)
        ENDDO

        CHSZ=0.005
        CHLN=0.004
        IF(NTYPE.EQ.0) THEN
          ECAREA(2)=XDISP
          ECAREA(1)=ECAREA(2)-CHLN*LNCHMX
          ECAREA(3)=0.0
          ECAREA(4)=ECAREA(3)+CHSZ*NTCH
        ELSE IF(NTYPE.EQ.1) THEN
          ECAREA(1)=XREF
          ECAREA(2)=ECAREA(1)+CHLN*LNCHMX
          ECAREA(3)=YREF
          ECAREA(4)=ECAREA(3)+CHSZ*NTCH
        ELSE IF(NTYPE.EQ.2) THEN
          ECAREA(2)=XREF
          ECAREA(1)=ECAREA(2)-CHLN*LNCHMX
          ECAREA(3)=YREF
          ECAREA(4)=ECAREA(3)+CHSZ*NTCH
        ELSE IF(NTYPE.EQ.3) THEN
          ECAREA(1)=(XREF-XMIN)/(XMAX-XMIN)*0.49*XDISP
          ECAREA(2)=ECAREA(1)+CHLN*LNCHMX
          IF(NJT.EQ.2) THEN
            ECAREA(4)=YDISP-0.49*XDISP+
     '        (YREF-YMIN)/(YMAX-YMIN)*0.49*XDISP
          ELSE IF(NJT.EQ.3) THEN
            ECAREA(4)=YDISP-0.49*XDISP+
     '        (YREF-ZMIN)/(ZMAX-ZMIN)*0.49*XDISP
          ENDIF
          ECAREA(3)=ECAREA(4)-CHSZ*NTCH
        ELSE IF(NTYPE.EQ.4) THEN
          ECAREA(2)=(XREF-XMIN)/(XMAX-XMIN)*0.49*XDISP
          ECAREA(1)=ECAREA(2)-CHLN*LNCHMX
          IF(NJT.EQ.2) THEN
            ECAREA(4)=YDISP-0.49*XDISP+
     '        (YREF-YMIN)/(YMAX-YMIN)*0.49*XDISP
          ELSE IF(NJT.EQ.3) THEN
            ECAREA(4)=YDISP-0.49*XDISP+
     '        (YREF-ZMIN)/(ZMAX-ZMIN)*0.49*XDISP
          ENDIF
          ECAREA(3)=ECAREA(4)-CHSZ*NTCH
        ELSE IF(NTYPE.EQ.5) THEN
          ECAREA(1)=XREF*XDISP
          ECAREA(2)=ECAREA(1)+CHLN*LNCHMX
          ECAREA(3)=YREF*YDISP
          ECAREA(4)=ECAREA(3)+CHSZ*NTCH
        ELSE IF(NTYPE.EQ.6) THEN
          ECAREA(2)=XREF*XDISP
          ECAREA(1)=ECAREA(2)-CHLN*LNCHMX
          ECAREA(3)=YREF*YDISP
          ECAREA(4)=ECAREA(3)+CHSZ*NTCH
        ELSE IF(NTYPE.EQ.7) THEN
          ECAREA(1)=XREF*XDISP
          ECAREA(2)=ECAREA(1)+CHLN*LNCHMX
          ECAREA(4)=YREF*YDISP
          ECAREA(3)=ECAREA(4)-CHSZ*NTCH
        ELSE IF(NTYPE.EQ.8) THEN
          ECAREA(2)=XREF*XDISP
          ECAREA(1)=ECAREA(2)-CHLN*LNCHMX
          ECAREA(4)=YREF*YDISP
          ECAREA(3)=ECAREA(4)-CHSZ*NTCH
        ELSE IF(NTYPE.EQ.9) THEN
          ECAREA(1)=XREF
          ECAREA(2)=ECAREA(1)+CHLN*LNCHMX
          ECAREA(4)=YREF
          ECAREA(3)=ECAREA(4)-CHSZ*NTCH
        ELSE IF(NTYPE.EQ.10) THEN
          ECAREA(2)=XREF
          ECAREA(1)=ECAREA(2)-CHLN*LNCHMX
          ECAREA(4)=YREF
          ECAREA(3)=ECAREA(4)-CHSZ*NTCH
        ENDIF
        IF(DOP) THEN
          WRITE(IOOP,'('' XREF='',E13.5,'' YREF='',E13.5)') XREF,YREF
          WRITE(IOOP,'('' DISP='',E13.5,'' ECAREA(1..4):'',4E13.5)')
     '      DISP,(ECAREA(I),I=1,4)
        ENDIF
        CALL GL_INIT_CHOICE(IW,IDEVICE,ECAREA)
      ENDIF
      IF(DOP) WRITE(IOOP,*) 'Mode=',MODE
      IF(MODE.EQ.'REQUEST') THEN
C       get next event
        CALL EVENT(ID_WS,ID_DEVICE,INSTAT,CLASS,IDATA,RDATA,SDATA,
     '    ERROR,*9999)
C       if its a choice option in the correct window
        IF(INSTAT.EQ.1.AND.CLASS(1:6).EQ.'CHOICE'.AND.ID_WS.EQ.IW)THEN
C         return choice
          NOCH=IDATA
        ELSE
        ENDIF
C       delete menu
        CALL GL_DELETE_CHOICE(IW)
      ENDIF

      IF(MODE.EQ.'REQUEST'.AND.NTYPE.EQ.0) THEN
        CO(NOCO)=OPTION(NOCH)
        CALL STRING_TRIM(STRING,IBEG,IEND)
        STRING=STRING(1:IEND)//' '//CO(NOCO)
        IF(.NOT.NULL(CO(NOCO))) THEN
          CO(NOCO+1)='?'
        ENDIF
      ENDIF

      CALL EXITS('CHOICE')
      RETURN
 9999 CALL ERRORS('CHOICE',ERROR)
      CALL EXITS('CHOICE')
      RETURN 1
      END


      SUBROUTINE CIRCLE(TYPE,IW,RADIUS,Z,NLINES,ERROR,*)

C**** Draws polyline (TYPE='POLYLINE') or fill-area (TYPE='FILL AREA') circle
C**** on workstation IW with radius RADIUS about the point Z(nj)
C**** with NLINES lines (max 100).

      REAL Z(*),PTS(3,100)
      CHARACTER TYPE*(*)
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:trans00.cmn'

      CALL ENTERS('CIRCLE',*9999)
      DO I=1,NLINES+1
        THETA=2.0*PI*REAL(I-1)/REAL(NLINES)
        IF(NJT.EQ.2) THEN
          PTS(1,I)=Z(1)+RADIUS*COS(THETA)
          PTS(2,I)=Z(2)+RADIUS*SIN(THETA)
          PTS(3,I)=0.0
        ELSE IF(NJT.EQ.3) THEN
          IF(IW.EQ.1) THEN
            PTS(1,I)=Z(1)+RADIUS*COS(THETA)
            PTS(2,I)=0.0
            PTS(3,I)=Z(3)+RADIUS*SIN(THETA)
          ELSE IF(IW.EQ.2) THEN
            PTS(1,I)=0.0
            PTS(2,I)=Z(2)+RADIUS*COS(THETA)
            PTS(3,I)=Z(3)+RADIUS*SIN(THETA)
          ELSE IF(IW.EQ.3) THEN
            CALL ZZ(Z,Z,TRANS)
            PTS(1,I)=Z(1)+RADIUS*COS(THETA)
            PTS(2,I)=Z(2)+RADIUS*SIN(THETA)
            PTS(3,I)=0.0
          ENDIF
        ENDIF
      ENDDO

      IF(TYPE(1:8).EQ.'POLYLINE') THEN
        CALL POLYLINE(1,IW,NLINES+1,PTS,ERROR,*9999)
      ELSE IF(TYPE(1:9).EQ.'FILL AREA') THEN
        CALL FILL_AREA(1,IW,NLINES,PTS,ERROR,*9999)
      ENDIF

      CALL EXITS('CIRCLE')
      RETURN
 9999 CALL ERRORS('CIRCLE',ERROR)
      CALL EXITS('CIRCLE')
      RETURN 1
      END


      SUBROUTINE CLOSE_PRINT_FILE(TYPE,ERROR,*)

C**** Dectivates print workstation and creates another printfile with
C**** extension .POST which is easier to edit than the .GKS file.
C**** For IW=3 print wkst is 16
C**** TYPE can be 'POSTSCRIPT' (print wkst is 15 )
C****          or 'METAFILE'   (print wkst is 17 )

      CHARACTER ERROR*(*),TYPE*(*),CFROMI*5
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:gks000.cmn'
      INCLUDE 'cmiss$reference:gks001.cmn'


      CALL ENTERS('CLOSE_PRINT_FILE',*9999)

      IF(TYPE(1:10).EQ.'POSTSCRIPT') THEN
        IF(IW.EQ.3) THEN
          CALL DAWK(16,1,ERROR,*9999)
          CLOSE(UNIT=16)
        ELSE
          CALL DAWK(15,1,ERROR,*9999)
          CLOSE(UNIT=15)
        ENDIF

        !Convert linefeeds to returns in GKS or PHIGS postscript files
        CALL POST(IW,FILE00,ERROR,*9999)

      ELSE IF(TYPE(1:8).EQ.'METAFILE') THEN
        CALL DAWK(17,1,ERROR,*9999)
        CLOSE(UNIT=17)
      ENDIF

      CALL EXITS('CLOSE_PRINT_FILE')
      RETURN
 9999 CALL ERRORS('CLOSE_PRINT_FILE',ERROR)
      CALL EXITS('CLOSE_PRINT_FILE')
      RETURN 1
      END


      SUBROUTINE CLOSE_SEGMENT(ISEGNUM,IW,ERROR,*)

C**** Closes graphics segment ISEGNUM.

      CHARACTER ERROR*(*)

      CALL ENTERS('CLOSE_SEGMENT',*9999)

      CALL GL_CLOSE_SEGMENT(ISEGNUM)

      CALL EXITS('CLOSE_SEGMENT')
      RETURN
 9999 CALL ERRORS('CLOSE_SEGMENT',ERROR)
      CALL EXITS('CLOSE_SEGMENT')
      RETURN 1
      END


      SUBROUTINE CLOSE_WS(IW,ERROR,*)

C**** Closes graphics workstation IW.

      CHARACTER ERROR*(*)

      CALL ENTERS('CLOSE_WS',*9999)

      CALL GL_CLOSE_WS(IW)

      CALL EXITS('CLOSE_WS')
      RETURN
 9999 CALL ERRORS('CLOSE_WS',ERROR)
      CALL EXITS('CLOSE_WS')
      RETURN 1
      END


      SUBROUTINE CLWS(ERROR,*)

C**** Deactivates and closes workstations.

      CHARACTER ERROR*(*)
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'

      CALL ENTERS('CLWS',*9999)
      DO IW=1,99
        IF(IW.NE.5.AND.IW.NE.6)THEN
    IF(IWKS(IW).GT.0) CALL GL_CLOSE_WS(IW)
        ENDIF
        IWKS(IW)=0
      ENDDO

      CALL EXITS('CLWS')
      RETURN
 9999 CALL ERRORS('CLWS',ERROR)
      CALL EXITS('CLWS')
      RETURN 1
      END


!new  CO and STRING passed to this routine don't do anything
      SUBROUTINE COLOUR(COLOUR_INDEX,IW,CO,STRING,XREF,YREF,ERROR,*)

C**** Returns colour index

      INTEGER COLOUR_INDEX
      CHARACTER CLASS*8,ERROR(*)*(*),OPTION(3)*20
      CHARACTER CO(*)*(*),STRING*(*)
      LOGICAL ABBREV,CONTINUE,UPDATE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'

      DATA LD1/1/,YFACTOR/0.136/

      CALL ENTERS('COLOUR',*9999)
      OPTION(1)='Set_Valuators'
      OPTION(2)='Set_Current'
      OPTION(3)='Return'
      CALL CHOICE('COLOUR_INDEX',1,1,INSTAT,72,'EVENT',1,3,NOCH,NOCO,2,
     '  CO,OPTION,STRING,XREF,YREF+1.2*YFACTOR*YDISP,ERROR,*9999)
      CALL VALUATOR('RED'  ,81,'EVENT',2,0.,1.,0.5,
     '  VALUE,XREF,YREF+4.0*YFACTOR*YDISP,ERROR,*9999)
      CALL VALUATOR('GREEN',82,'EVENT',2,0.,1.,0.5,
     '  VALUE,XREF,YREF+3.0*YFACTOR*YDISP,ERROR,*9999)
      CALL VALUATOR('BLUE' ,83,'EVENT',2,0.,1.,0.5,
     '  VALUE,XREF,YREF+2.0*YFACTOR*YDISP,ERROR,*9999)
      CONTINUE=.TRUE.
      DO WHILE (CONTINUE)
        CALL EVENT(ID_WS,ID_Device,Input_Status,CLASS,NOCH,
     '    VALUE,SDATA,ERROR,*9999)
        IF(DOP) WRITE(*,*) ' Input_Class=',Class
        IF(CLASS(1:6).EQ.'CHOICE') THEN
          IF(NOCH.GT.0) THEN
            IF(ABBREV(OPTION(NOCH),'SET_CURRENT',8)) THEN
              CONTINUE=.FALSE.
              UPDATE=.TRUE.
            ELSE IF(ABBREV(OPTION(NOCH),'RETURN',6)) THEN
              CONTINUE=.FALSE.
            ENDIF
          ENDIF
        ELSE IF(CLASS(1:8).EQ.'VALUATOR') THEN
          IF(ID_WS.EQ.81) THEN
            RED_INTENS=VALUE
          ELSE IF(ID_WS.EQ.82) THEN
            GREEN_INTENS=VALUE
          ELSE IF(ID_WS.EQ.83) THEN
            BLUE_INTENS=VALUE
          ENDIF
        ENDIF
      ENDDO
      IF(UPDATE) THEN
        CALL GL_SET_COLOUR_REP(IW,COLOUR_INDEX,RED_INTENS,GREEN_INTENS,
     '    BLUE_INTENS)
      ENDIF
      CALL GL_DELETE_CHOICE(72)
      CALL GL_DELETE_VALUATOR(81)
      CALL GL_DELETE_VALUATOR(82)
      CALL GL_DELETE_VALUATOR(83)

      CALL EXITS('COLOUR')
      RETURN

 9999 CALL ERRORS('COLOUR',ERROR)
      CALL EXITS('COLOUR')
      RETURN 1
      END


      SUBROUTINE COWK(IW,ERROR,*)

C**** Closes and reopens workstation identified by IW (in order to remove
C**** workstation viewport from screen).

      INTEGER LINE_COLOUR_INDEX,LINE_INDEX,
     '  LINE_TYPE1,LINE_TYPE2,LINE_TYPE3,LINE_TYPE4,
     '  MARKER_COLOUR_INDEX,MARKER_INDEX,MARKER_SIZE,
     '  MARKER_TYPE1,MARKER_TYPE2,MARKER_TYPE3,MARKER_TYPE4,
     '  MARKER_TYPE5,
     '  AREA_COLOUR_INDEX,AREA_INDEX,AREA_INT_STYLE,AREA_STYLE_INDEX,
     '  TEXT_INDEX,TEXT_COLOUR_INDEX,
     '  TEXT_FONT1,TEXT_FONT2,TEXT_FONT3,
     '  TEXT_PRECISION1,TEXT_PRECISION2,TEXT_PRECISION3
      REAL LINE_WIDTH1,LINE_WIDTH2,LINE_WIDTH3,LINE_WIDTH4
      CHARACTER ERROR(*)*(*)
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:gks001.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'

      DATA LD1/1/

      CALL ENTERS('COWK',*9999)
      IF(IWKS(IW).EQ.2) THEN
        CALL GL_DAWK(IW)
        IWKS(IW)=1
      ENDIF
      IF(IWKS(IW).EQ.1) THEN
        CALL GL_CLOSE_WS(IW)
        IWKS(IW)=0
      ENDIF
      IF(IW.EQ.1) THEN      !close choice menu and valuators
        CALL GL_DELETE_CHOICE(91)
        CALL GL_DELETE_VALUATOR(83)
        CALL GL_DELETE_VALUATOR(84)
      ELSE IF(IW.EQ.2) THEN !close choice menu and valuators
        CALL GL_DELETE_CHOICE(92)
        CALL GL_DELETE_VALUATOR(85)
        CALL GL_DELETE_VALUATOR(86)
      ELSE IF(IW.EQ.3) THEN !close choice menu and valuators
        CALL GL_DELETE_CHOICE(93)
        CALL GL_DELETE_VALUATOR(87)
        CALL GL_DELETE_VALUATOR(88)
        CALL GL_DELETE_VALUATOR(89)
      ELSE IF(IW.EQ.4) THEN !close choice menu and valuators
        CALL GL_DELETE_CHOICE(94)
        CALL GL_DELETE_VALUATOR(81)
        CALL GL_DELETE_VALUATOR(82)
      ENDIF
      IF(IWKS(IW).EQ.0) THEN
        CALL SETUP(IW,ERROR)
      ENDIF

      CALL EXITS('COWK')
      RETURN

 9999 CALL ERRORS('COWK',ERROR)
      CALL EXITS('COWK')
      RETURN 1
      END


      SUBROUTINE CREATE_SEGMENT(ISEGNUM,IW,ERROR,*)

C**** create graphics segment ISEGNUM.

      CHARACTER ERROR*(*)

      CALL ENTERS('CREATE_SEGMENT',*9999)

      CALL GL_CREATE_SEGMENT(ISEGNUM,IW)

      CALL EXITS('CREATE_SEGMENT')
      RETURN
 9999 CALL ERRORS('CREATE_SEGMENT',ERROR)
      CALL EXITS('CREATE_SEGMENT')
      RETURN 1
      END


      SUBROUTINE DAWK(IW,ID,ERROR,*)

C**** Deactivate workstation identified by IW.
C**** If ID=1 deferred output is released.

      CHARACTER ERROR*(*)
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'


      CALL ENTERS('DAWK',*9999)
      IWKS(IW)=1
      IF(IW.NE.5.AND.IW.NE.6)THEN
        CALL GL_DAWK(IW)           !deactivate this window
      ENDIF

      CALL EXITS('DAWK')
      RETURN
 9999 CALL ERRORS('DAWK',ERROR)
      CALL EXITS('DAWK')
      RETURN 1
      END


      SUBROUTINE DBOX(IW,XDIM,YDIM,XWC,YWC,ZWC,ERROR,*)

C**** Draws box of dimensions 2*XDIM,2*YDIM around the point XWC,YWC,ZWC

      REAL Z(3,5)
      INCLUDE 'cmiss$reference:geom00.cmn'


      CALL ENTERS('DBOX',*9999)
      IF(IW.EQ.1) THEN
        IF(NJT.EQ.2) THEN !plot x against y
          Z(1,1)=XWC-XDIM
          Z(2,1)=YWC-YDIM
          Z(1,2)=XWC+XDIM
          Z(2,2)=Z(2,1)
          Z(1,3)=Z(1,2)
          Z(2,3)=YWC+YDIM
          Z(1,4)=Z(1,1)
          Z(2,4)=Z(2,3)
          Z(1,5)=Z(1,1)
          Z(2,5)=Z(2,1)
          CALL POLYLINE(1,IW,5,Z,ERROR,*9999)
        ELSE IF(NJT.EQ.3) THEN !plot x against z
          Z(1,1)=XWC-XDIM
          Z(3,1)=ZWC-YDIM
          Z(1,2)=XWC+XDIM
          Z(3,2)=Z(3,1)
          Z(1,3)=Z(1,2)
          Z(3,3)=ZWC+YDIM
          Z(1,4)=Z(1,1)
          Z(3,4)=Z(3,3)
          Z(1,5)=Z(1,1)
          Z(3,5)=Z(3,1)
          CALL POLYLINE(1,IW,5,Z,ERROR,*9999)
        ENDIF
      ELSE IF(IW.EQ.2) THEN !plot Y against Z
        Z(2,1)=YWC-XDIM
        Z(3,1)=ZWC-YDIM
        Z(2,2)=YWC+XDIM
        Z(3,2)=Z(3,1)
        Z(2,3)=Z(2,2)
        Z(3,3)=ZWC+YDIM
        Z(2,4)=Z(2,1)
        Z(3,4)=Z(3,3)
        Z(2,5)=Z(2,1)
        Z(3,5)=Z(3,1)
        CALL POLYLINE(1,IW,5,Z,ERROR,*9999)
      ELSE IF(IW.EQ.3) THEN !plot X against Y
        Z(1,1)=XWC-XDIM
        Z(2,1)=YWC-YDIM
        Z(1,2)=XWC+XDIM
        Z(2,2)=Z(2,1)
        Z(1,3)=Z(1,2)
        Z(2,3)=YWC+YDIM
        Z(1,4)=Z(1,1)
        Z(2,4)=Z(2,3)
        Z(1,5)=Z(1,1)
        Z(2,5)=Z(2,1)
        CALL POLYLINE(1,IW,5,Z,ERROR,*9999)
      ELSE IF(IW.EQ.4) THEN
      ENDIF

      CALL EXITS('DBOX')
      RETURN
 9999 CALL ERRORS('DBOX',ERROR)
      CALL EXITS('DBOX')
      RETURN 1
      END


      SUBROUTINE DELETE_SEGMENT(ISEGNUM,ISEG,IW,ERROR,*)

C**** Deletes graphics segment ISEGNUM.

      INTEGER ISEG(*)
      CHARACTER ERROR*(*)

      CALL ENTERS('DELETE_SEGMENT',*9999)

      CALL GL_DELETE_SEGMENT(ISEGNUM,IW)
      ISEG(ISEGNUM)=0
      ISEGNUM=0

      CALL EXITS('DELETE_SEGMENT')
      RETURN
 9999 CALL ERRORS('DELETE_SEGMENT',ERROR)
      CALL EXITS('DELETE_SEGMENT')
      RETURN 1
      END


      SUBROUTINE DETECT(IW,ISEG,ISEGNUM,CLASS,ERROR,*)

C**** Change ISEGNUM to CLASS='DETECTABLE' or 'UNDETECTABLE'
C**** Change ISEG(ISEGNUM) to 3 if detectable, 2 if not.

      INTEGER ISEG(*)
      CHARACTER CLASS*(*),ERROR*(*)

      CALL ENTERS('DETECT',*9999)
      IF(CLASS(1:10).EQ.'DETECTABLE') THEN
        ISEG(ISEGNUM)=3
      ELSE
        ISEG(ISEGNUM)=2
      ENDIF
      CALL GL_DETECT(ISEGNUM,ISEG(ISEGNUM),IW)

      CALL EXITS('DETECT')
      RETURN
 9999 CALL ERRORS('DETECT',ERROR)
      CALL EXITS('DETECT')
      RETURN 1
      END


      SUBROUTINE DETRAN(ISPLOT,NOCO,NTCO,NTCOQU,CO,COQU,STRING,ERROR,*)

C**** Defines 3D transformations using Phigs

      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:back00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'cmiss$reference:ntsg00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
      INCLUDE 'cmiss$reference:phig00.cmn'
      INCLUDE 'cmiss$reference:view00.cmn'
      INTEGER ISPLOT(NHM,0:NEM,*),NTCOQU(*)
      REAL FPT(3),SCL(3),SFT(3)
      LOGICAL ABBREV,CHANGE,CONTINUE,EVENT,FILIO
      CHARACTER CO(*)*(*),COQU(NXCO,*)*(*),CUPPER*(255),ERROR*(*),
     '  STRING*(MXCH),CHOOSE*30,FILE*50,
     '  OPTION(10)*80,STATUS*3

      DATA LD1/1/

      CALL ENTERS('DETRAN',*9999)
 1    IF(CO(NOCO+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        WRITE(IOOP,'(X,A)') STRING(1:IEND),';l/p/r/w'//
     '    '<;FILENAME>['//FILE00(IBEG1:IEND1)//']<;doc>'
      ELSE IF(CO(NOCO+1).EQ.'??') THEN
        CALL DOCUM('fe14','doc','DETRAN',ERROR,*9999)
      ELSE
        CALL CHECKQ('ELPRW',NOCO,1,CO,COQU,STRING,*1)
        FILIO=.FALSE.
        EVENT=.FALSE.
        IF(ABBREV(COQU(NOCO,1),'P',1)) THEN
          FILIO=.TRUE.
          IOTYPE=1
          STATUS='NEW'
        ELSE IF(ABBREV(COQU(NOCO,1),'R',1)) THEN
          FILIO=.TRUE.
          IOTYPE=2
          STATUS='OLD'
        ELSE IF(ABBREV(COQU(NOCO,1),'W',1)) THEN
          FILIO=.TRUE.
          IOTYPE=3
          STATUS='NEW'
        ELSE IF(ABBREV(COQU(NOCO,1),'L',1)) THEN
          FILIO=.TRUE.
          IOTYPE=4
          STATUS='OLD'
        ELSE IF(ABBREV(COQU(NOCO,1),'E',1)) THEN
          EVENT=.TRUE.
        ENDIF

        IF(FILIO) THEN
          IVDU=IOIP
          IFILE=1
          CALL CHECKF(2,NOCO,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.IPTRAN',STATUS,
     '      'DIRECT','FORMATTED',132,ERROR,*9999)
!old      CALL IPT1(ERROR,*9999)
          CLOSE(UNIT=IFILE)
          IF(IOTYPE.EQ.1.OR.IOTYPE.EQ.2)THEN
            CALL SET_GLOBAL_XFORM(3,ISVIEW,ANGLE(1),ANGLE(2),ANGLE(3),
     '        FIXED_PT,SCALE,SHIFT,ERROR,*9999)
            MODE_PROJ=0 !PHIGS$K_PARALLEL
            FPT(1)=0.0
            FPT(2)=0.0
            FPT(3)=0.0
            SFT(1)=0.0
            SFT(2)=0.0
            SFT(3)=0.0
            SCL(1)=1.0
            SCL(2)=1.0
            SCALE(3)=1.0
            CALL SET_VIEW(ISTATUS,3,MODE_PROJ,NPC_CLIP,
     '        0,0,0,
     '        BACK_PLANE_DIST_NEW,FRONT_PLANE_DIST_NEW,
     '        FPT,SCL,SFT,PROJ_REF_PT_NEW,VIEW_PLANE_DIST_NEW,
     '        VIEW_PLANE_NEW,VIEW_REF_PT_NEW,VIEW_UP_NEW,
     '        VIEWPORT,WINDOW_NEW,ERROR,*9999)
          ENDIF
        ENDIF

      ENDIF

      CALL EXITS('DETRAN')
      RETURN
 9999 CALL ERRORS('DETRAN',ERROR)
      CALL EXITS('DETRAN')
      RETURN 1
      END


!new  DEVIEW gone


      SUBROUTINE DISPLAY(ISEG,NOCO,NTCO,NTCOQU,CO,COQU,CSEG,END,STRING,
     '    ERROR,*)

C**** Displays string in a window.

      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'cmiss$reference:echo00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
      INTEGER NTCOQU(*)
      LOGICAL END,ABBREV,CONTINUE
      CHARACTER STRING*(MXCH),CO(*)*(*),COQU(NXCO,*)*(*),ERROR*(*),
     '  BLANK*1,MESSAGE*80

      DATA BLANK/' '/

      CALL ENTERS('DISPLAY',*9999)
      MESSAGE=CO(2)
      ECAREA(1)=0.5*XDISP
      ECAREA(2)=ECAREA(1)+0.5*XDISP
      ECAREA(3)=0.3*YDISP
      ECAREA(4)=ECAREA(3)+0.2*YDISP
      CALL STRING_TRIM(MESSAGE,IBEG,IEND)

      CALL EXITS('DISPLAY')
      RETURN
 9999 CALL ERRORS('DISPLAY',ERROR)
      CALL EXITS('DISPLAY')
      RETURN 1
      END


      SUBROUTINE DISPLAY_FILE(IW,NOCO,NTFILE,CO,FILE_EXT,FILE_NAME,
     '  STRING,ERROR,*)

C**** Displays files in current directory and returns chosen one.
C**** If none chosen FILE_NAME is 'Exit'.

      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
      LOGICAL MORE_FILES
      CHARACTER CO(*)*(*),ERROR*(*),FILE_EXT*(*),FILE_NAME*(*),
     '  STRING*(MXCH),CHOOSE*20,OPTION(22)*20

      CALL ENTERS('DISPLAY_FILE',*9999)
      NOFILE_START=0
 10   CALL FIND_FILE(NOFILE_START,NTFILE,FILE_EXT,OPTION,ERROR,*9999)
      IF(NTFILE.EQ.20) THEN
        MORE_FILES=.TRUE.
      ELSE
        MORE_FILES=.FALSE.
      ENDIF
      NOFILE_DEFAULT=0
      DO NOFILE=1,NTFILE !to display default in square brackets
        CALL STRING_TRIM(OPTION(NOFILE),IBEG1,IEND1)
        CALL STRING_TRIM(FILE00,IBEG2,IEND2)
        IF(OPTION(NOFILE)(IBEG1:IEND1).EQ.FILE00(IBEG2:IEND2)) THEN
          DO N1FILE=NOFILE,2,-1
            OPTION(N1FILE)=OPTION(N1FILE-1)
          ENDDO
          OPTION(1)='['//FILE00(IBEG2:IEND2)//']'
          NOFILE_DEFAULT=1
        ENDIF
      ENDDO
      IF(MORE_FILES) THEN
        OPTION(NTFILE+1)='...'
        NTCH=NTFILE+2
      ELSE
        NTCH=NTFILE+1
      ENDIF
      OPTION(NTCH)='Exit'
      CALL PRECHOICE1(2,IW,NOCH,NOCO,NTCH,CO,'REQUEST',OPTION,STRING,
     '  ERROR,*9999)
      CHOOSE=OPTION(NOCH)
      IF(CHOOSE(1:3).EQ.'...') THEN
        NOFILE_START=NOFILE_START+20
        GO TO 10 !to display more files
      ELSE IF(CHOOSE(1:4).EQ.'Exit') THEN
        FILE_NAME='Exit'
      ELSE       !return chosen filename
        IF(NOCH.EQ.NOFILE_DEFAULT) THEN
          FILE_NAME=FILE00
        ELSE
          FILE_NAME=OPTION(NOCH)
          FILE00=FILE_NAME
        ENDIF
      ENDIF

      CALL EXITS('DISPLAY_FILE')
      RETURN
 9999 CALL ERRORS('DISPLAY_FILE',ERROR)
      CALL EXITS('DISPLAY_FILE')
      RETURN 1
      END


      SUBROUTINE DOCUM(FILE,EXTEN,STRING,ERROR,*)

C**** Displays documentation
C**** FILE identifies file to be read (file ID=8 and workstation ID=8)
C**** EXTEN is name of file extension (eg .DOC)
C**** STRING identifies string in documentation file at beginning & end
C****   of text to be displayed
C**** Documentation has root plus number of branches (max 20)
C**** Root and branches can each have arbitrary number of windows
C**** NTWIND is number of windows for root or branch
C**** NTLINE is number of lines in root or branch (this is divided by 23
C****   to give the number of windows)

      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cmis00.cmn'
      INCLUDE 'cmiss$reference:comm00.cmn'
      INCLUDE 'cmiss$reference:docd00.cmn'
      INCLUDE 'cmiss$reference:docu00.cmn'
      INCLUDE 'cmiss$reference:gks000.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
      INCLUDE 'cmiss$reference:ntsg00.cmn'
      REAL P(3)
      LOGICAL ABBREV,FOUND,ROOT
      CHARACTER ERROR*(*),FILE*(*),EXTEN*(*),STRING*(MXCH),
     '  CFROMI*3,CUPPER*80

      DATA LD1/1/

      CALL ENTERS('DOCUM',*9999)
      CALL STRING_TRIM(FILE,IBEG1,IEND1)
      CALL STRING_TRIM(EXTEN,IBEG2,IEND2)
      CALL OPENF(8,'DISK','/usr/people/cmiss/vms/document/'
     '  //FILE(IBEG1:IEND1)//'.'//EXTEN(IBEG2:IEND2),
     '  'OLD','DIRECT','FORMATTED',132,ERROR,*9998)
      CALL ACWK(8,0,ERROR,*9999)
      ROOT=.TRUE.
      CDATA(2)=CUPPER(STRING)
 10   CALL STRING_TRIM(CDATA(2),IBEG3,IEND3)
      IF(DOP) WRITE(IOOP,'('' String: '',A)') CDATA(2)(IBEG3:IEND3)
      IF(ROOT) THEN
        CALL FIND(8,CDATA(2)(IBEG3:IEND3),1,2000,LINE1,FOUND,
     '    ERROR,*9999)
        if(dop) write(ioop,*) ' line1=',line1,' found=',found
        CALL FIND(8,CDATA(2)(IBEG3:IEND3),LINE1+1,LINE1+100,LINE2,FOUND,
     '    ERROR,*9999)
        if(dop) write(ioop,*) ' line2=',line2,' found=',found
        LINE0=LINE2
        ROOT=.FALSE.
      ELSE IF(.NOT.ROOT) THEN
        CALL FIND(8,CDATA(2)(IBEG3:IEND3),LINE0,2000,LINE1,FOUND,
     '    ERROR,*9999)
        if(dop) write(ioop,*) ' line1=',line1,' found=',found
        CALL FIND(8,CDATA(2)(IBEG3:IEND3),LINE1+1,LINE1+100,LINE2,FOUND,
     '    ERROR,*9999)
        if(dop) write(ioop,*) ' line2=',line2,' found=',found
      ENDIF
      NTLINE=LINE2-LINE1-1
      NTWIND=1+INT((NTLINE-1)/23)
      if(dop) write(ioop,*) ' ntwind=',ntwind
      LINE3=0
      DO NOWIND=1,NTWIND
        DO I=1,MIN(23,NTLINE-LINE3)
          READ(8,FMT='(A)',REC=LINE1+LINE3+I) CDATA(1)(1:80)
          if(dop) write(ioop,'(1x,a)') cdata(1)(1:80)
          P(1)=0.0
          P(2)=1.0-I*0.025
          CALL TEXT(1,2,8,CDATA(1)(1:80),P,ERROR,*9999)
        ENDDO
        IF(NOWIND.LT.NTWIND) THEN
          IF(COMMAND) THEN
            P(1)=0.0
            P(2)=0.375
            CALL TEXT(1,2,8,'Press RETURN for next window',P,
     '        ERROR,*9999)
            FORMAT='($,'' >>[next]: '')'
          ELSE
            CALL CREATE_SEGMENT(NTSG+1,8,ERROR,*9999)
            P(1)=0.9
            P(2)=0.375
            CALL TEXT(1,2,8,'<NEXT>',P,ERROR,*9999)
            CALL CLOSE_SEGMENT(NTSG+1,8,ERROR,*9999)

            CALL DETECT(NTSG+1,'DETECTABLE')
            CALL PICK(8,LD1,INSTAT,ISEGM,IPICK)
            IF(DOP) WRITE(*,*) 'INSTAT=',INSTAT
            CALL DELETE_SEGMENT(NTSG+1,ISEG,IW,ERROR,*9999)
          ENDIF
          LINE3=23*NOWIND
        ELSE IF(NOWIND.EQ.NTWIND) THEN
          IF(COMMAND) THEN
            P(1)=0.0
            P(2)=0.375
            CALL TEXT(1,2,8,'Press RETURN to close window',P,
     '        ERROR,*9999)
            FORMAT='($,'' >>[close]: '')'
          ELSE
            CALL CREATE_SEGMENT(NTSG+1,8,ERROR,*9999)
            P(1)=0.9
            P(2)=0.375
            CALL TEXT(1,2,8,'<CLOSE>',P,ERROR,*9999)
            CALL CLOSE_SEGMENT(NTSG+1,8,ERROR,*9999)
            CALL DETECT(NTSG+1,'DETECTABLE')
            CALL PICK(8,LD1,INSTAT,ISEGM,IPICK)
            IF(DOP) WRITE(*,*) 'INSTAT=',INSTAT
            CALL DELETE_SEGMENT(NTSG+1,ISEG,IW,ERROR,*9999)
          ENDIF
        ENDIF
        IF(COMMAND) THEN
          CALL GINOUT(IOTYPE,2,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CBLANK,80,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
c         CALL CINOUT(1,IOIP,0,FORMAT,1,CDATA,CBLANK,80,INFO,
c    '      ERROR,*9999)
        ELSE
          CDATA(1)=' '
        ENDIF
  CALL GL_CLEAR_WS(8)
        CDATA(1)=CUPPER(CDATA(1))
        IF(.NOT.NULL(CDATA(1))) THEN
          CALL STRING_TRIM(CDATA(1),IBEG4,IEND4)
          IF(ABBREV(CDATA(1),'EXIT',1)) THEN
            GO TO 160
          ELSE
            CALL STRING_TRIM(CDATA(2),IBEG3,IEND3)
            CDATA(2)=CDATA(2)(IBEG3:IEND3)//' '//CDATA(1)(IBEG4:IEND4)
            GO TO 10
          ENDIF
        ENDIF
      ENDDO
 160  CALL DAWK(8,0,ERROR,*9999)
      CALL COWK(8,ERROR,*9999)
      CLOSE(UNIT=8)

      CALL EXITS('DOCUM')
      RETURN
 9998 CALL ERRORS('DOCUM',ERROR)
      CALL EXITS('DOCUM')
      RETURN 1
 9999 CALL ERRORS('DOCUM',ERROR)
      CALL EXITS('DOCUM')
      CALL DAWK(8,0,ERROR,*9999)
      CALL COWK(8,ERROR,*9999)
      CLOSE(UNIT=8)
      RETURN 1
      END


      SUBROUTINE ELLIPSE(SEMIAXIS1,SEMIAXIS2,XWC,YWC,NLINES,ERROR,*)

C**** Draws fill-area ellipse about the point XWC,YWC with semi-axis
C**** lengths SEMIAXIS1 and SEMIAXIS2 with NLINES lines (max 100).

      REAL ELLPSE(3,100)
      CHARACTER ERROR*(*)
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'

      CALL ENTERS('ELLIPSE',*9999)
      IW=1 !needs to be passed to ellipse
      DO I=1,NLINES
        THETA=2.0*PI*REAL(I-1)/REAL(NLINES)
        ELLPSE(1,I)=XWC+SEMIAXIS1*COS(THETA)
        ELLPSE(2,I)=YWC+SEMIAXIS2*SIN(THETA)
      ENDDO
      CALL FILL_AREA(0,IW,NLINES,ELLPSE,ERROR,*9999)

      CALL EXITS('ELLIPSE')
      RETURN
 9999 CALL ERRORS('ELLIPSE',ERROR)
      CALL EXITS('ELLIPSE')
      RETURN 1
      END


      SUBROUTINE EVENT(ID_WS,ID_DEVICE,ID_STATUS,CLASS,IDATA,RDATA,
     '  SDATA,ERROR,*)

C***  Handles event mode input from physical graphics devices.
C***  Returns: ID_WS      workstation identifier
C***           CLASS      input class string, can be one of the following:
C***                      'NONE'
C***                      'LOCATOR'
C***                      'VALUATOR'
C***                      'CHOICE'
C***                      'PICK'
C***                      'STRING'
C***           ID_DEVICE  logical input device number
C***           ID_STATUS  1 for input ok, 0 for end of input.
C***           IDATA      data input by choice
C***                      or segment and pick number for pick
C***                      or view trans.number for locator
C***           RDATA      data input by valuator or locator
C***           SDATA       "     "    " string

      INTEGER IDATA(*)
      REAL RDATA(*)
      CHARACTER CLASS*(*),SDATA*(*),ERROR*(*)
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'

      CALL ENTERS('EVENT',*9999)
C     get next event from event queue if present
      CALL GL_EVENT(ID_WS,ICLASS,ID_DEVICE,ID_STATUS)
      IF(DOP) THEN
        WRITE(IOOP,'('' ID_WS='',I4)') ID_WS
        WRITE(IOOP,'('' ID_Device='',I4)') ID_Device
      ENDIF
      !a status of one signals end of input
      IF(ID_STATUS.EQ.1)THEN
        IF(ICLASS.EQ.1) THEN
          CLASS='CHOICE'
          CALL GL_GET_CHOICE(ID_WS,IDATA,ID_STATUS)
        ELSE IF(ICLASS.EQ.5) THEN
          CLASS='VALUATOR'
          CALL GL_GET_VALUATOR(ID_WS,RDATA,ID_STATUS)
        ELSE IF(ICLASS.EQ.2) THEN
          CLASS='LOCATOR'
          CALL GL_GET_LOCATOR(ID_WS,RDATA(1),RDATA(2),ID_STATUS)
        ELSE IF(ICLASS.EQ.3) THEN
          CLASS='PICK'
          CALL GL_GET_PICK(ID_WS,ID,ID_STATUS,IDATA(1),IDATA(2))
        ELSE IF(ICLASS.EQ.4) THEN
          CLASS='STRING'
          CALL GL_GET_STRING(IW_WS,SDATA,ISIZE,ID_STATUS) !this needs checking
        ELSE
          CLASS='NONE'
          IDATA(1)=0
        ENDIF
      ENDIF
      IF(DOP) WRITE(IOOP,'('' Input Class='',A)') CLASS

      CALL EXITS('EVENT')
      RETURN
 9999 CALL ERRORS('EVENT',ERROR)
      CALL EXITS('EVENT')
      RETURN 1
      END


      SUBROUTINE FILL_AREA(IBUNDLE,IW,NPOINTS,POINTS,ERROR,*)

C**** Draws a fill-area on IW with index IBUNDLE. POINTS(1..3,np) is an array
C**** containing the 3D coordinates of each point. If the IBUNDLE is 0 the
C**** primative will use the previously defined fill-area index.
C**** WARNING: IW=4 will change the values in the POINTS array!

      REAL POINTS(3,*)
      CHARACTER ERROR*(*)
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:curr00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:map000.cmn'

      CALL ENTERS('FILL_AREA',*9999)
      IF(NPOINTS.EQ.0) GOTO 9998
      IF(DOP) THEN
        WRITE(IOOP,'(/'' POINTS(1,1..):'',4E12.3)')
     '    (POINTS(1,NP),NP=1,4)
        WRITE(IOOP,'( '' POINTS(2,1..):'',4E12.3)')
     '    (POINTS(2,NP),NP=1,4)
        WRITE(IOOP,'( '' POINTS(3,1..):'',4E12.3)')
     '    (POINTS(3,NP),NP=1,4)
      ENDIF

      IF(IBUNDLE.NE.0) CALL GL_SET_FILL_AREA_INDEX(IW,IBUNDLE)

      IF(IW.LE.3)THEN
  CALL GL_FILL_AREA(NPOINTS,POINTS)

      ELSE IF(IW.EQ.4) THEN
        IF(PROJEC(1:11).EQ.'RECTANGULAR') THEN !points x and y
        ELSE IF(PROJEC(1:2).EQ.'XI') THEN !assume points are in Xi coords
          DO NP=1,NPOINTS
            POINTS(1,NP)=-1.0+2.*(REAL(MXI1-1)+POINTS(1,NP))/MAX_XI
            POINTS(2,NP)=-1.0+2.*(REAL(MXI2-1)+POINTS(2,NP))/MAX_XI
          ENDDO
        ELSE  !assume points are in polar coords
          CALL MAP4(NPOINTS,POINTS,ERROR,*9999)
        ENDIF
        CALL GL_FILL_AREA(NPOINTS,POINTS)
      ENDIF

 9998 CALL EXITS('FILL_AREA')
      RETURN
 9999 CALL ERRORS('FILL_AREA',ERROR)
      CALL EXITS('FILL_AREA')
      RETURN 1
      END


      SUBROUTINE FLUSH_DEVICE_EVENTS(IW,CLASS,IFLAG,ERROR,*)

C**** Flush queues
C**** CLASS can be 'PICK', 'CHOICE', 'LOCATOR', 'VALUATOR' or 'ALL'

      CHARACTER CLASS*(*),ERROR*(*)

      CALL ENTERS('FLUSH_DEVICE_EVENTS',*9999)
      CALL GL_FLUSH_DEVICE_EVENTS(CLASS,IFLAG)

      CALL EXITS('FLUSH_DEVICE_EVENTS')
      RETURN
 9999 CALL ERRORS('FLUSH_DEVICE_EVENTS',ERROR)
      CALL EXITS('FLUSH_DEVICE_EVENTS')
      RETURN 1
      END


      SUBROUTINE FREEWK(MNWK,MXWK,NOWK,ERROR,*)

C**** Finds the first free workstation NOWK between MNWK and MXWK.

      INTEGER MNWK,MXWK,NOWK,ERSTAT,CONID,WKTY
      CHARACTER ERROR*(*)

      CALL ENTERS('FREEWK',*9999)
      CALL ASSERT((MNWK.GE.0).AND.(MXWK.LE.99),
     '  'Workstation limits are invalid',ERROR,*9999)
      DO 1 N1WK=MNWK,MXWK
        CALL GL_INQ_WS(N1WK,ERSTAT)
        IF (ERSTAT.EQ.25) THEN
          NOWK=N1WK
          GOTO 2
        ELSE IF (ERSTAT.NE.0) THEN
          ERROR='Error occurred while searching for free workstation'
          GOTO 9999
        ENDIF
 1    CONTINUE
      ERROR='No free workstations exist between specified limits'
      GOTO 9999
 2    CALL EXITS('FREEWK')
      RETURN
 9999 CALL ERRORS('FREEWK',ERROR)
      CALL EXITS('FREEWK')
      RETURN 1
      END

      SUBROUTINE GCLRWK(A,*)
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      DIMENSION A(*)
      CHARACTER ERROR*10
      WRITE(OP_STRING,*) '>>Link with GKS'
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
9999  RETURN 1
      END


      SUBROUTINE GKS_DRAW(IW,ISEG,CSEG,STRING,ERROR,*)


      CALL EXITS('GKS_DRAW')
      RETURN
 9999 CALL ERRORS('GKS_DRAW',ERROR)
      CALL EXITS('GKS_DRAW')
      RETURN 1
      END


      SUBROUTINE GKS_STRG(IW,NCHAR,PROMPT_STRING,TEXT_STRING,ERROR,*)

C**** Prompts the user for text under GKS.
C**** Note: Workstation IW is activated and deactivated in this routine.


      CALL EXITS('GKS_STRG')
      RETURN
 9999 CALL ERRORS('GKS_STRG',ERROR)
      CALL EXITS('GKS_STRG')
      RETURN 1
      END


      SUBROUTINE GKSTEXT(IW,INSTAT,PROMPT_STRING,TEXT_STRING,XWC,YWC,
     '  ERROR,*)

C**** Prompts the user for text under GKS.
C**** XDC1,YDC1 are device coords of bottom left of IW viewport
C**** XDC2,YDC2 are device coords of top right   of IW viewport
C**** X_MIN,Y_MIN are world  coords of bottom left of IW viewport
C**** X_MAX,Y_MAX are world  coords of top right   of IW viewport
C**** XWC,YWC are returned world coordinates of reference point
C**** INSTAT is 1 if successful 0 if not

      CALL EXITS('GKSTEXT')
      RETURN
 9999 CALL ERRORS('GKSTEXT',ERROR)
      CALL EXITS('GKSTEXT')
      RETURN 1
      END


      SUBROUTINE INPUT_MODE(IW,ID,CLASS,MODE,ERROR,*)

C**** resets input mode for workstation IW device ID.
C***  CLASS  can be  'LOCATOR','CHOICE','PICK','VALUATOR'
C***  MODE   can be  'REQUEST','EVENT','SAMPLE'

      CHARACTER CLASS*(*),MODE*(*),ERROR*(*)

      CALL ENTERS('INPUT_MODE',*9999)

      CALL EXITS('INPUT_MODE')
      RETURN
 9999 CALL ERRORS('INPUT_MODE',ERROR)
      CALL EXITS('INPUT_MODE')
      RETURN 1
      END


!new  INVIS not called by anything?


      SUBROUTINE LINE3D(NTDX,XL,ERROR,*)

C**** Draws line on 3D viewport at rect.cart. world coords.

      REAL XL(20,*),Z1(3),Z2(3)
      CHARACTER ERROR*(*)
      INCLUDE 'cmiss$reference:fbgr00.cmn'
      INCLUDE 'cmiss$reference:trans00.cmn'

      CALL ENTERS('LINE3D',*9999)
      DO 200 NODX=1,NTDX
        Z1(1)=XL(NODX,1)
        Z1(2)=XL(NODX,2)
        Z1(3)=XL(NODX,3)
        IF(RHTRAN) THEN
          CALL FBIMGPT(5,Z1(1),Z1(2),Z1(3),Z2(1),Z2(2))
          XL(NODX,1)=Z2(1)
          XL(NODX,2)=Z2(2)
        ELSE IF(LFTRAN) THEN
          CALL FBIMGPT(6,Z1(1),Z1(2),Z1(3),Z2(1),Z2(2))
          XL(NODX,1)=Z2(1)
          XL(NODX,2)=Z2(2)
        ELSE
          CALL ZZ(Z1,Z2,TRANS)
          XL(NODX,1)=Z2(1)
          XL(NODX,2)=Z2(2)
          XL(NODX,3)=Z2(3)
        ENDIF
 200  CONTINUE

      CALL EXITS('LINE3D')
      RETURN
 9999 CALL ERRORS('LINE3D',ERROR)
      CALL EXITS('LINE3D')
      RETURN 1
      END


      SUBROUTINE LICOLREP(NOCO,NTCO,NTCOQU,CO,COQU,STRING,ERROR,*)

C***  Lists colour representation thru inquiry functions

      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:gks001.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
      REAL COL(3)
      LOGICAL CBBREV
      CHARACTER STRING*(MXCH),CO(*)*(*),COQU(NXCO,*)*(*),ERROR*(*)

      CALL ENTERS('LICOLREP',*9999)

      IF(CO(NOCO+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        WRITE(IOOP,'(X,A)') 'FEM List colours on WS[3]'
      ELSE IF(CO(NOCO+1).EQ.'??') THEN
        CALL DOCUM('fe14','doc','LICOLREP',ERROR,*9999)
      ELSE
        IF(CBBREV(CO,'ON',1,NOCO+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),1,NTIL,IW,ERROR,*9999)
        ELSE
          IW=3
        ENDIF
        IF(IW.EQ.3)THEN
          WRITE(*,*)' NINDICES=',NINDICES
          DO INDEX=0,NINDICES-1
            CALL GL_INQ_COLOUR_REP(IW,INDEX,IERR,COL)
            WRITE(*,*)INDEX,COL(1),COL(2),COL(3)
          ENDDO
        ENDIF
      ENDIF

      CALL EXITS('LICOLREP')
      RETURN
 9999 CALL ERRORS('LICOLREP',ERROR)
      CALL EXITS('LICOLREP')
      RETURN 1
      END


!new  added STRING doesn't work
      SUBROUTINE LISTRU(CO,STRING,ERROR,*)

C***  Lists graphical data structure thru inquiry functions

      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
      INCLUDE 'cmiss$reference:view00.cmn'
      INTEGER ELEM_NUM,ELEM_TYPE,ELEM_SIZE
      CHARACTER CO(*)*(*),STRING*(MXCH),ERROR*(*)

      CALL ENTERS('LISTRU',*9999)

 1    IF(CO(NOCO+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        WRITE(IOOP,'(X,A)') 'FEM List Structure'
      ELSE IF(CO(NOCO+1).EQ.'??') THEN
        CALL DOCUM('fe14','doc','LISTRU',ERROR,*9999)
      ELSE
        ELEM_NUM=0
        ELEM_TYPE=1
        DO WHILE (ELEM_NUM.LE.20)
          CALL GL_INQ_ELEM(ISVIEW,ELEM_NUM,IERR,ELEM_TYPE,
     '      ELEM_SIZE)
          WRITE(*,'('' Isview element'',I3,'' error,type,size'',3I5)')
     '      ELEM_NUM,IERR,ELEM_TYPE,ELEM_SIZE
          ELEM_NUM=ELEM_NUM+1
        ENDDO
      ENDIF

      CALL EXITS('LISTRU')
      RETURN
 9999 CALL ERRORS('LISTRU',ERROR)
      CALL EXITS('LISTRU')
      RETURN 1
      END


      SUBROUTINE LOCATOR(INIT,IW,INSTAT,MODE,NECHO,XREF,XWC,YREF,YWC,
     '  ERROR,*)

C**** Calls GKS locator
C**** INIT specifies whether locator is to be initialised (1:yes,2:no)
C**** IW specifies workstation number (which is also transformation no)
C**** INSTAT is returned as 1 if locate is successful, 0 otherwise
C**** MODE is input mode - 'REQUEST','SAMPLE' or 'EVENT'
C**** NECHO is 0 on entry for no echo                 (not supported)
C**** NECHO is 3 on entry for default   prompt/echo
C**** NECHO is 4 on entry for line      prompt/echo   (not supported)
C**** NECHO is 5 on entry for rectangle prompt/echo   (not supported)
C**** NECHO is 6 on entry for digital   prompt/echo   (not supported)
C**** XREF,YREF are initial coords of echo
C**** XWC,YWC are returned world coords of located point

      INTEGER IARRAY(2)
      REAL RDATA(2)
      CHARACTER ERROR*(*),MODE*(*),CLASS*8,RECORD(2)*80
      LOGICAL TRIGGERED
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'

      DATA LD1/1/

      CALL ENTERS('LOCATOR',*9999)
      IF(INIT.EQ.1) THEN !initialise locator device
        IF(NECHO.EQ.3) THEN
          CALL GL_INIT_LOCATOR(IW,LD1,XREF,YREF,0,XDISP,0,YDISP)
        ELSE IF(NECHO.EQ.4) THEN
        ELSE IF(NECHO.EQ.5) THEN
        ELSE IF(NECHO.EQ.6) THEN
        ENDIF
C       make sure the appropriate pixel to world transformation is used
  CALL GL_SELECT_XFORM(IW)
      ENDIF
      IF(MODE(1:7).EQ.'REQUEST') THEN
  TRIGGERED=.FALSE.
  DO WHILE(.NOT.TRIGGERED)
          CALL EVENT(ID_WS,ID_DEVICE,INSTAT,CLASS,IDATA,RDATA,SDATA,
     '      ERROR,*9999)
C         if its a locator
          IF(INSTAT.EQ.1.AND.CLASS(1:8).EQ.'LOCATOR'.AND.ID_WS.EQ.IW)
     '      THEN
C           return cursor location
            XWC=RDATA(1)
            YWC=RDATA(2)
            TRIGGERED=.TRUE.
          ELSE IF(INSTAT.EQ.0) THEN
            TRIGGERED=.TRUE.
          ENDIF
  ENDDO
C       delete menu
  CALL GL_DELETE_LOCATOR(IW)
      ELSE IF(MODE(1:6).EQ.'SAMPLE') THEN
  CALL GL_GET_LOCATOR(IW,RDATA(1),RDATA(2),INSTAT)
  XWC=RDATA(1)
        YWC=RDATA(2)
      ENDIF
      IF(DOP) THEN
        WRITE(IOOP,'('' XWC='',E11.3,'' YWC='',E11.3)') XWC,YWC
      ENDIF

      CALL EXITS('LOCATOR')
      RETURN

 9999 CALL ERRORS('LOCATOR',ERROR)
      CALL EXITS('LOCATOR')
      RETURN 1
      END


      SUBROUTINE OPEN_PRINT_FILE(IW,TYPE,ERROR,*)

C**** Activates print workstation. If IW is 3 the phigs viewing transformations
C**** are set up. Otherwise the GKS workstation bundles are copied from IW to
C**** the print workstation.
C**** For IW=3 print wkst is 16
C**** TYPE can be 'POSTSCRIPT' (print wkst is 15 )
C**** or 'METAFILE' (print wkst is 17 )

      INTEGER TEXT_HORIZ_ALIGN,TEXT_VERT_ALIGN,TEXT_COLOUR_INDEX,
     '  TEXT_FONT,TEXT_PATH,TEXT_PRECISION
      CHARACTER ERROR*(*),TYPE*(*),CFROMI*5
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:cbzm00.cmn'
      INCLUDE 'cmiss$reference:ntsg00.cmn'
      INCLUDE 'cmiss$reference:phig00.cmn'
      INCLUDE 'cmiss$reference:view00.cmn'

      CALL ENTERS('OPEN_PRINT_FILE',*9999)


      CALL EXITS('OPEN_PRINT_FILE')
      RETURN
 9999 CALL ERRORS('OPEN_PRINT_FILE',ERROR)
      CALL EXITS('OPEN_PRINT_FILE')
      RETURN 1
      END


      SUBROUTINE OPEN_SEGMENT(ISEGNUM,ISEG,IW,CLABEL,INDEX,INDEX_OLD,
     '  NLABEL,IVIS,CSEG,ERROR,*)

C**** Opens graphics segment ISEGNUM. If ISEGNUM does not already exist
C**** (has the value 0), then NTSG is incremented and a new segment is created.
C**** Otherwise the old segment structure is replaced by the new structure.
C**** All subsequent call to graphics primitives, etc will be part of this
C**** structure. Note that GKS segments cannot be nested, so a call to
C**** OPEN_SEGMENT must be followed by a call to CLOSE_SEGMENT before the next
C**** OPEN_SEGMENT.
C**** When created CSEG(ISEGNUM) stores a string comprising the first 48
C**** characters of CLABEL together with INDEX(4 chars) & NLABEL(5 chars)
C**** then a '/' and the IW number(2 chars).
C**** INDEX is the bundle table index for the current segment and if an old
C**** segment is being replaced, INDEX_OLD is read from the CSEG string and
C**** returned for possible use in the replacement segment creation.
C**** IVIS=1 creates a visible segment
C**** IVIS=2 creates an invisible segment

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:ntsg00.cmn'
      INCLUDE 'cmiss$reference:view00.cmn'
!     Parameter List
      INTEGER INDEX,INDEX_OLD,ISEG(*),ISEGNUM,IVIS,IW,NLABEL
      CHARACTER CSEG(*)*(*),CLABEL*(*),ERROR*(*)
!     Local Variables
      INTEGER IFROMC
      CHARACTER CFROMI*5,CHAR2*5,CHAR4*5,CHAR5*5,CLABEL2*(52)

      CALL ENTERS('OPEN_SEGMENT',*9999)

      IF(ISEGNUM.EQ.0) THEN !segment does not already exist
        NTSG=NTSG+1
        ISEGNUM=NTSG
      ELSE !recover INDEX, delete old segment and use the old ISEG,CSEG labels
        INDEX_OLD=IFROMC(CSEG(ISEGNUM)(49:52))
        IF(DOP) WRITE(IOOP,'('' INDEX='',I3)') INDEX
        CALL GL_DELETE_SEGMENT(ISEGNUM,IW)
      ENDIF

! MPN 6/10/93 For some unexplained reason CFROMI doesn't work here
      write(unit=CHAR2,fmt='(I2)') IW
      write(unit=CHAR4,fmt='(I4)') INDEX
      write(unit=CHAR5,fmt='(I5)') NLABEL
c      CHAR2=CFROMI(IW,'(I2)')
c      CHAR4=CFROMI(INDEX,'(I4)')
c      CHAR5=CFROMI(NLABEL,'(I5)')
      CLABEL2(1:)=CLABEL
      ISEG(ISEGNUM)=2
      IF(IVIS.EQ.2) ISEG(ISEGNUM)=1
      CALL GL_CREATE_SEGMENT(ISEGNUM,IW)
      CALL GL_VISIB(ISEGNUM,IVIS,IW)
      CSEG(ISEGNUM)=CLABEL2(1:48)//CHAR4(1:4)//CHAR5(1:5)//'/'
     '  //CHAR2(1:2)

      CALL EXITS('OPEN_SEGMENT')
      RETURN
 9999 CALL ERRORS('OPEN_SEGMENT',ERROR)
      CALL EXITS('OPEN_SEGMENT')
      RETURN 1
      END


      SUBROUTINE PHIG(IW,ERROR,*)

C**** Defines initial transformations and view matrices for PHIGS.
C**** A_TRANS is the 4*4 transformation matrix containing rotations
C**** shifts and scaling of an object in world coords
C**** A_ORIENT is a 4*4 matrix indicating the orientation of the
C**** world coords in the viewing reference coords (z-axis into the
C**** screen, y-axis vertical)
C**** A_MAP is a 4*4 matrix specifying parallel or perspective views
C**** together with clipping planes

      CHARACTER ERROR*(*)
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:gks000.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'cmiss$reference:phig00.cmn'
      INCLUDE 'cmiss$reference:trac00.cmn'
      INCLUDE 'cmiss$reference:trans00.cmn'
      INCLUDE 'cmiss$reference:gl_constants.f'

      CALL ENTERS('PHIG',*9999)

C *** Reference point, rotation angles, scaling & shift vector
C *** for world coord transformation
C     FIXED_PT(1)=0.5*(XMIN+XMAX)
C     FIXED_PT(2)=0.5*(YMIN+YMAX)
      MODE_PROJ=GL_PARALLEL
      FIXED_PT(1)=0.0
      FIXED_PT(2)=0.0
      FIXED_PT(3)=0.0
      IF(DOP) WRITE(IOOP,'('' FIXED_PT:'',3E12.3)')
     '  (FIXED_PT(NJ),NJ=1,3)
      DO J=1,3
        ANGLE(J)=0.0
        SCALE(J)=1.0
        SHIFT(J)=0.0
      ENDDO

C *** Reference point, view plane normal & view up vector
C *** in world coords to orient viewing reference coords
      VIEW_REF_PT(1)=0.5*(XMIN+XMAX)
      VIEW_REF_PT(2)=0.5*(YMIN+YMAX)
      VIEW_REF_PT(3)=0.5*(ZMIN+ZMAX)
      VIEW_PLANE(1)=0.0
      VIEW_PLANE(2)=0.0
      VIEW_PLANE(3)=1.0
      VIEW_UP(1)=0.0
      VIEW_UP(2)=1.0
      VIEW_UP(3)=0.0

C *** Window, projection ref pt and viewplane, backplane & frontplane
C *** positions in viewing reference coords and viewport in norm proj
C *** coords for mapping to normalized projection coord system
      PROJ_REF_PT(1)= VIEW_REF_PT(1)      !changed from 0 for gl sgi
      PROJ_REF_PT(2)= VIEW_REF_PT(2)      !changed from 0 for gl sgi
      PROJ_REF_PT(3)=3*DIAG+VIEW_REF_PT(3) !changed from 3.0*DIAG for gl sgi
      IF(MODE_PROJ.EQ.GL_PARALLEL)THEN
        WINDOW(1)=XMIN-VIEW_REF_PT(1)
        WINDOW(2)=XMAX-VIEW_REF_PT(1)
        WINDOW(3)=YMIN-VIEW_REF_PT(2)
        WINDOW(4)=YMAX-VIEW_REF_PT(2)
        FRONT_PLANE_DIST= 2.5*DIAG
        BACK_PLANE_DIST = 3.5*DIAG
      ELSE
        WINDOW(1)=XMIN-VIEW_REF_PT(1)
        WINDOW(2)=XMAX-VIEW_REF_PT(1)
        WINDOW(3)=YMIN-VIEW_REF_PT(2)
        WINDOW(4)=YMAX-VIEW_REF_PT(2)
        FRONT_PLANE_DIST= 2.5*DIAG
        BACK_PLANE_DIST = 3.5*DIAG
      ENDIF
      VIEW_PLANE_DIST = 0.0
      IF(DOP) WRITE(IOOP,'('' Frontplane at '',E12.3,
     '  '' Backplane at '',E12.3)') FRONT_PLANE_DIST,BACK_PLANE_DIST
      VIEWPORT(1)= 0.0
      VIEWPORT(2)= 1.0
      VIEWPORT(3)= 0.0
      VIEWPORT(4)= 1.0
      VIEWPORT(5)= 0.0   !7-feb-1990
      VIEWPORT(6)= 1.0   !7-feb-1990

C *** Clipping limits in normalized projection coords
      NPC_CLIP(1)= 0.0
      NPC_CLIP(2)= 1.0
      NPC_CLIP(3)= 0.0
      NPC_CLIP(4)= 1.0
      NPC_CLIP(5)= 0.0   !7-feb-1990
      NPC_CLIP(6)= 1.0   !7-feb-1990

      CALL SET_VIEW(ISTATUS,IW,MODE_PROJ,NPC_CLIP,
     '   ANGLE(1),ANGLE(2),ANGLE(3),
     '   BACK_PLANE_DIST,FRONT_PLANE_DIST,
     '   FIXED_PT,SCALE,SHIFT,PROJ_REF_PT,VIEW_PLANE_DIST,
     '   VIEW_PLANE,VIEW_REF_PT,VIEW_UP,
     '   VIEWPORT,WINDOW,ERROR,*9999)
      CALL SET_GLOBAL_XFORM(IW,ISVIEW,ANGLE(1),ANGLE(2),ANGLE(3),
     '   FIXED_PT,SCALE,SHIFT,ERROR,*9999)

C *** Initialize
      VIEW_PLANE_DIST_NEW =VIEW_PLANE_DIST
      BACK_PLANE_DIST_NEW =BACK_PLANE_DIST
      FRONT_PLANE_DIST_NEW=FRONT_PLANE_DIST
      DO NJ=1,3
        PROJ_REF_PT_NEW(NJ)=PROJ_REF_PT(NJ)
        VIEW_REF_PT_NEW(NJ)=VIEW_REF_PT(NJ)
        VIEW_PLANE_NEW(NJ)=VIEW_PLANE(NJ)
        VIEW_UP_NEW(NJ)=VIEW_UP(NJ)
      ENDDO
      WINDOW_NEW(1)= WINDOW(1)
      WINDOW_NEW(2)= WINDOW(2)
      WINDOW_NEW(3)= WINDOW(3)
      WINDOW_NEW(4)= WINDOW(4)

      CALL EXITS('PHIG')
      RETURN
 9999 CALL ERRORS('PHIG',ERROR)
      CALL EXITS('PHIG')
      RETURN 1
      END


      SUBROUTINE PICK(IW,ID,MODE,INSTAT,IPICKRET,IPICKID,ERROR,*)

C**** Calls GKS pick
C**** MODE can be 'REQUEST' or 'EVENT'
C**** IPICKRET is returned segment on exit if in request mode
C**** IPCKID   is returned pick identifier if in request mode
C**** INSTAT   ir returned 1 if successful input triggered or else 0

      LOGICAL TRIGGERED
      INTEGER IDATA(2)
      CHARACTER MODE*(*),ERROR*(*),CLASS*8
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:gks000.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'

      DATA LD1/1/,XFACTOR/0.15/,YFACTOR/0.1/,DUMMY/' '/

      CALL ENTERS('PICK',*9999)

      IF(MODE(1:7).EQ.'REQUEST') THEN
        CALL GL_INIT_PICK(IW,ID)
  TRIGGERED=.FALSE.
  DO WHILE(.NOT.TRIGGERED)
          CALL EVENT(ID_WS,ID_DEVICE,INSTAT,CLASS,IDATA,RDATA,SDATA,
     '      ERROR,*9999)
C         if its a pick
          IF(INSTAT.EQ.1.AND.CLASS(1:8).EQ.'PICK'.AND.ID_WS.EQ.IW
     '      .AND.IDATA(1).GT.0)THEN
C            return picked segment
      IPICKRET=IDATA(1)
            TRIGGERED=.TRUE.
          ELSE IF(INSTAT.EQ.0) THEN
            TRIGGERED=.TRUE.
          ENDIF
  ENDDO
        CALL GL_DELETE_PICK(IW)
      ELSE IF(MODE(1:5).EQ.'EVENT') THEN
      ENDIF

      CALL EXITS('PICK')
      RETURN

 9999 CALL ERRORS('PICK',ERROR)
      CALL EXITS('PICK')
      RETURN 1
      END


      SUBROUTINE POLYLINE(IBUNDLE,IW,NPOINTS,POINTS,ERROR,*)

C**** Draws a polyline on IW with index IBUNDLE. POINTS(1..3,np) is an array
C**** containing the 3D coordinates of each point. If the IBUNDLE is 0 the
C**** primative will use the previously defined polyline index.
C**** WARNING: IW=4 will change the values in the POINTS array!

      REAL POINTS(3,*)
      CHARACTER ERROR*(*)
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:curr00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:map000.cmn'


      CALL ENTERS('POLYLINE',*9999)
      IF(NPOINTS.EQ.0) GOTO 9998

      IF(NJT.EQ.2.OR.IW.EQ.31.OR.IW.EQ.32)THEN
  DO N=1,NPOINTS
    POINTS(3,N)=0.0
        ENDDO
      ENDIF
      IF(DOP) THEN
        WRITE(IOOP,'('' IW='',I2,'' NPOINTS='',I6)') IW,NPOINTS
        DO N=1,NPOINTS
          WRITE(IOOP,'('' POINTS(nj,'',I6,''): '',3E12.3)')
     '      N,(POINTS(NJ,N),NJ=1,NJT)
        ENDDO
      ENDIF

      IF(IBUNDLE.NE.0) CALL GL_SET_POLYLINE_INDEX(IW,IBUNDLE)
      IF(IW.LE.3) THEN
        CALL GL_POLYLINE(NPOINTS,POINTS)

      ELSE IF(IW.EQ.4) THEN !Map window
        IF(PROJEC(1:11).EQ.'RECTANGULAR') THEN !points x and y
          CALL GL_POLYLINE(NPOINTS,POINTS)
        ELSE IF(PROJEC(1:2).EQ.'XI') THEN !assume points are in Xi coords
          DO NP=1,NPOINTS
            POINTS(1,NP)=-1.0+2.*(REAL(MXI1-1)+POINTS(1,NP))/MAX_XI
            POINTS(2,NP)=-1.0+2.*(REAL(MXI2-1)+POINTS(2,NP))/MAX_XI
          ENDDO
          CALL GL_POLYLINE(NPOINTS,POINTS)
        ELSE  !assume points are in polar coords
Cold      DO NP=1,NPOINTS-1
Cold        CALL MAP4(2,POINTS(1,NP),ERROR,*9999)
Cold        IF(ABS(POINTS(1,NP)-POINTS(2,NP)).LT.0.5) THEN !ok to plot
Cold          CALL GL_POLYLINE(NPOINTS,POINTS(1,NP))
Cold        ENDIF
Cold      ENDDO
Cnews     AAY 4 June 91
    CALL MAP4(1,POINTS(1,1),ERROR,*9999)
    DO NP=2,NPOINTS
      CALL MAP4(1,POINTS(1,NP),ERROR,*9999)
      IF(ABS(POINTS(1,NP)-POINTS(1,NP-1)).LT.0.5) THEN !ok to plot
              CALL GL_POLYLINE(2,POINTS(1,NP-1))
            ENDIF
          ENDDO
Cnewe
        ENDIF

      ELSE IF(IW.EQ.5.OR.IW.EQ.6) THEN
        !frame grabber windows

      ELSE IF(IW.EQ.7) THEN
        !plot Y against z

      ELSE IF(IW.EQ.10) THEN
        ! Time history plots
        CALL MAP10(NPOINTS,POINTS)
        CALL GL_POLYLINE(NPOINTS,POINTS)

      ELSE IF(IW.EQ.11) THEN
        CALL MAP10(NPOINTS,POINTS)
        CALL GL_POLYLINE(NPOINTS,POINTS)

      ELSE IF(IW.EQ.21.OR.IW.EQ.22.OR.IW.EQ.23
     '  .OR.IW.EQ.30.OR.IW.EQ.31.OR.IW.EQ.32.OR.IW.EQ.34.OR.IW.EQ.35
     '  .OR.IW.EQ.40.OR.IW.EQ.41.OR.IW.EQ.42.OR.IW.EQ.43.OR.IW.EQ.44
     '  .OR.IW.EQ.45.OR.IW.EQ.46.OR.IW.EQ.47
     '  .OR.IW.EQ.50.OR.IW.EQ.51
     '  .OR.IW.EQ.60.OR.IW.EQ.61.OR.IW.EQ.62.OR.IW.EQ.63.OR.IW.EQ.64
     '  .OR.IW.EQ.65.OR.IW.EQ.66.OR.IW.EQ.68.OR.IW.EQ.69) THEN !plots
        !plot x against y
        CALL GL_POLYLINE(NPOINTS,POINTS)

      ENDIF


 9998 CALL EXITS('POLYLINE')
      RETURN
 9999 CALL ERRORS('POLYLINE',ERROR)
      CALL EXITS('POLYLINE')
      RETURN 1
      END


      SUBROUTINE POLYMARKER(IBUNDLE,IW,NPOINTS,POINTS,ERROR,*)

C**** Draws a polymarker on IW with index IBUNDLE. POINTS(1..3,np) is an array
C**** containing the 3D coordinates of each point. If the IBUNDLE is 0 the
C**** primitive will use the previously defined polymarker index.
C**** WARNING: IW=4 will change the values in the POINTS array!

      REAL POINTS(3,*),Z(3)
      CHARACTER ERROR*(*)
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'cmiss$reference:map000.cmn'

      CALL ENTERS('POLYMARKER',*9999)
      IF(NPOINTS.EQ.0) GOTO 9998

      IF(DOP) THEN
        WRITE(IOOP,'('' IW='',I2,'' NPOINTS='',I6)') IW,NPOINTS
        DO N=1,NPOINTS
          WRITE(IOOP,'('' POINTS(nj,'',I6,''): '',3E12.3)')
     '      N,(POINTS(NJ,N),NJ=1,NJT)
        ENDDO
      ENDIF

      IF(IBUNDLE.NE.0) CALL GL_SET_POLYMARKER_INDEX(IW,IBUNDLE)
      IF(IW.LE.3) THEN
        CALL GL_POLYMARKER(NPOINTS,POINTS)

      ELSE IF(IW.EQ.4) THEN !Map window
        IF(PROJEC(1:11).EQ.'RECTANGULAR') THEN !points x and y
          CALL GL_POLYMARKER(NPOINTS,POINTS)
        ELSE IF(PROJEC(1:2).EQ.'XI') THEN !assume points are in Xi coords
          DO NP=1,NPOINTS
            POINTS(1,NP)=-1.0+2.*(REAL(MXI1-1)+POINTS(1,NP))/MAX_XI
            POINTS(2,NP)=-1.0+2.*(REAL(MXI2-1)+POINTS(2,NP))/MAX_XI
          ENDDO
          CALL GL_POLYMARKER(NPOINTS,POINTS)
        ELSE  !assume points are in polar coords
Cold      DO NP=1,NPOINTS-1
Cold        CALL MAP4(2,POINTS(1,NP),ERROR,*9999)
Cold        IF(ABS(POINTS(1,NP)-POINTS(2,NP)).LT.0.5) THEN !ok to plot
Cold          CALL GL_POLYMARKER(NPOINTS,POINTS(1,NP))
Cold        ENDIF
Cold      ENDDO
Cnews     AAY 4 June 91
    CALL MAP4(NPOINTS,POINTS,ERROR,*9999)
    CALL GL_POLYMARKER(NPOINTS,POINTS)
Cnewe
        ENDIF

      ELSE IF(IW.EQ.5.OR.IW.EQ.6) THEN
        !frame grabber windows

      ELSE IF(IW.EQ.7) THEN
        !plot Y against z

      ELSE IF(IW.EQ.10) THEN
        ! Time history plots
        CALL MAP10(NPOINTS,POINTS)
        CALL GL_POLYMARKER(NPOINTS,POINTS)

      ELSE IF(IW.EQ.11) THEN
        CALL MAP10(NPOINTS,POINTS)
        CALL GL_POLYMARKER(NPOINTS,POINTS)

      ELSE IF(IW.EQ.21.OR.IW.EQ.22.OR.IW.EQ.23
     '  .OR.IW.EQ.30.OR.IW.EQ.31.OR.IW.EQ.32.OR.IW.EQ.34.OR.IW.EQ.35
     '  .OR.IW.EQ.40.OR.IW.EQ.41.OR.IW.EQ.42.OR.IW.EQ.43.OR.IW.EQ.44
     '  .OR.IW.EQ.45.OR.IW.EQ.46.OR.IW.EQ.47
     '  .OR.IW.EQ.50.OR.IW.EQ.51
     '  .OR.IW.EQ.60.OR.IW.EQ.61.OR.IW.EQ.62.OR.IW.EQ.63.OR.IW.EQ.64
     '  .OR.IW.EQ.65.OR.IW.EQ.66.OR.IW.EQ.68.OR.IW.EQ.69) THEN !plots
        !plot x against y
        CALL GL_POLYMARKER(NPOINTS,POINTS)

      ENDIF


 9998 CALL EXITS('POLYMARKER')
      RETURN
 9999 CALL ERRORS('POLYMARKER',ERROR)
      CALL EXITS('POLYMARKER')
      RETURN 1
      END


      SUBROUTINE PRECHOICE1(IP,IW,NOCH,NOCO,NTCH,CO,MODE,OPTION,
     '  STRING,ERROR,*)

C**** Sets up choice menu associated with finite element windows.
C**** IP=1 is standard finite element window choices.
C**** IP=2 is subsequent choice.

      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
      CHARACTER CO(*)*(*),ERROR*(*),MODE*(*),OPTION(*)*(*),STRING*(MXCH)

      CALL ENTERS('PRECHOICE1',*9999)
      IF(IP.EQ.1) THEN
        IF(NJT.LE.2) THEN
          IF(IW.EQ.1) THEN
            CALL CHOICE('DEFINE',1,1,INSTAT,7,MODE,1,NTCH,NOCH,NOCO,9,
     '        COD,OPTION,STRING,0.01*XDISP,YDISP-0.45*XDISP,ERROR,*9999)
          ENDIF
        ELSE IF(NJT.EQ.3) THEN
          IF(IW.EQ.1) THEN
            CALL CHOICE('DEFINE',1,1,INSTAT,7,MODE,1,NTCH,NOCH,NOCO,9,
     '        COD,OPTION,STRING,0.01*DISP,0.50*DISP,ERROR,*9999)
          ELSE IF(IW.EQ.2) THEN
            CALL CHOICE('DEFINE',1,1,INSTAT,7,MODE,1,NTCH,NOCH,NOCO,9,
     '        COD,OPTION,STRING,1.00*DISP,0.99*DISP,ERROR,*9999)
          ELSE IF(IW.EQ.3) THEN
            CALL CHOICE('DEFINE',1,1,INSTAT,7,MODE,1,NTCH,NOCH,NOCO,10,
     '        COD,OPTION,STRING,1.00*DISP,0.48*DISP,ERROR,*9999)
          ELSE IF(IW.EQ.4) THEN
            CALL CHOICE('DEFINE',1,1,INSTAT,7,MODE,1,NTCH,NOCH,NOCO,9,
     '        COD,OPTION,STRING,0.11*DISP,0.19*DISP,ERROR,*9999)
          ENDIF
        ENDIF

      ELSE IF(IP.EQ.2) THEN
        IF(NJT.LE.2) THEN
          IF(IW.EQ.1) THEN
            CALL CHOICE('DEFINE',1,1,INSTAT,7,MODE,1,NTCH,NOCH,NOCO,9,
     '        CO,OPTION,STRING,0.01*XDISP,YDISP-0.45*XDISP,ERROR,*9999)
          ENDIF
        ELSE IF(NJT.EQ.3) THEN
          IF(IW.EQ.1) THEN
            CALL CHOICE('DEFINE',1,1,INSTAT,7,MODE,1,NTCH,NOCH,NOCO,9,
     '        CO,OPTION,STRING,0.01*DISP,0.50*DISP,ERROR,*9999)
          ELSE IF(IW.EQ.2) THEN
            CALL CHOICE('DEFINE',1,1,INSTAT,7,MODE,1,NTCH,NOCH,NOCO,9,
     '        CO,OPTION,STRING,1.00*DISP,0.99*DISP,ERROR,*9999)
          ELSE IF(IW.EQ.3) THEN
            CALL CHOICE('DEFINE',1,1,INSTAT,7,MODE,1,NTCH,NOCH,NOCO,9,
     '        CO,OPTION,STRING,0.95*DISP,0.48*DISP,ERROR,*9999)
          ELSE IF(IW.EQ.4) THEN
            CALL CHOICE('DEFINE',1,1,INSTAT,7,MODE,1,NTCH,NOCH,NOCO,9,
     '        CO,OPTION,STRING,0.11*DISP,0.19*DISP,ERROR,*9999)
          ENDIF
        ENDIF
      ENDIF

      CALL EXITS('PRECHOICE1')
      RETURN
 9999 CALL ERRORS('PRECHOICE1',ERROR)
      CALL EXITS('PRECHOICE1')
      RETURN 1
      END


      SUBROUTINE PRECHOICE2(IW,NOCO,CO,FILE_EXT,FILE_NAME,PTYPE,
     '  Q,QUALIFIERS,STRING,ERROR,*)

C**** Sets up choice menu for type of qualifier.

      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
      CHARACTER CO(*)*(*),Q*(*),ERROR*(*),FILE_EXT*(*),FILE_NAME*(*),
     '  PTYPE*(*),QUALIFIERS*(*),STRING*(MXCH),
     '  CHOOSE*20,FULL_BRIEF*5,OPTION(22)*20
      LOGICAL CELEM,CHOOSE_FULL_BRIEF

      CALL ENTERS('PRECHOICE2',*9999)
      NOCH=0
      IF(CELEM('c',QUALIFIERS)) THEN
        NOCH=NOCH+1
        OPTION(NOCH)='calculate'
      ENDIF
      IF(CELEM('d',QUALIFIERS)) THEN
        NOCH=NOCH+1
        OPTION(NOCH)='default'
      ENDIF
      IF(CELEM('l',QUALIFIERS)) THEN
        NOCH=NOCH+1
        OPTION(NOCH)='list'
      ENDIF
      IF(CELEM('m',QUALIFIERS)) THEN
        NOCH=NOCH+1
        OPTION(NOCH)='mouse'
      ENDIF
      IF(CELEM('p',QUALIFIERS)) THEN
        NOCH=NOCH+1
        OPTION(NOCH)='prompt'
      ENDIF
      IF(CELEM('r',QUALIFIERS)) THEN
        NOCH=NOCH+1
        OPTION(NOCH)='read'
      ENDIF
      IF(CELEM('s',QUALIFIERS)) THEN
        NOCH=NOCH+1
        OPTION(NOCH)='segment'
      ENDIF
      IF(CELEM('w',QUALIFIERS)) THEN
        NOCH=NOCH+1
        OPTION(NOCH)='write'
      ENDIF
      OPTION(NOCH+1)='FULL'
      OPTION(NOCH+2)='Exit'
      NTCH=NOCH+2
      FULL_BRIEF='BRIEF'
      CHOOSE_FULL_BRIEF=.TRUE.
      DO WHILE(CHOOSE_FULL_BRIEF)
        CALL PRECHOICE1(2,IW,NOCH,NOCO,NTCH,CO,'REQUEST',OPTION,STRING,
     '    ERROR,*9999)
        CHOOSE=OPTION(NOCH)
        IF(CHOOSE(1:4).EQ.'FULL') THEN
          OPTION(NTCH-1)='BRIEF'
          FULL_BRIEF='FULL'
          CHOOSE_FULL_BRIEF=.TRUE.
        ELSE IF(CHOOSE(1:5).EQ.'BRIEF') THEN
          OPTION(NTCH-1)='FULL'
          FULL_BRIEF='BRIEF'
          CHOOSE_FULL_BRIEF=.TRUE.
        ELSE IF(CHOOSE(1:9).EQ.'calculate') THEN
          Q='c'
          CHOOSE_FULL_BRIEF=.FALSE.
        ELSE IF(CHOOSE(1:7).EQ.'default') THEN
          Q='d'
          CHOOSE_FULL_BRIEF=.FALSE.
        ELSE IF(CHOOSE(1:4).EQ.'list') THEN
          Q='l'
          CHOOSE_FULL_BRIEF=.FALSE.
        ELSE IF(CHOOSE(1:5).EQ.'mouse') THEN
          Q='m'
          CHOOSE_FULL_BRIEF=.FALSE.
        ELSE IF(CHOOSE(1:6).EQ.'prompt') THEN
          Q='p'
          CHOOSE_FULL_BRIEF=.FALSE.
        ELSE IF(CHOOSE(1:4).EQ.'read') THEN
          Q='r'
          CHOOSE_FULL_BRIEF=.FALSE.
        ELSE IF(CHOOSE(1:7).EQ.'segment') THEN
          Q='s'
          CHOOSE_FULL_BRIEF=.FALSE.
        ELSE IF(CHOOSE(1:5).EQ.'write') THEN
          Q='w'
          CHOOSE_FULL_BRIEF=.FALSE.
        ELSE IF(CHOOSE(1:4).EQ.'Exit') THEN
          Q='e'
          CHOOSE_FULL_BRIEF=.FALSE.
        ENDIF
      ENDDO

      PTYPE(1:)=' '
      IF(FILE_EXT(1:6).EQ.'IPDATA'.AND.Q.EQ.'s') THEN
        OPTION( 1)='geometry'
        OPTION( 2)='fibre'
        OPTION( 3)='field'
        OPTION( 4)='numbers'
        OPTION( 5)='projections'
        OPTION( 6)='values'
        OPTION( 7)='trace'
        OPTION( 8)='Exit'
        NTCH=8
        CALL PRECHOICE1(2,IW,NOCH,NOCO,NTCH,CO,'REQUEST',OPTION,STRING,
     '    ERROR,*9999)
        CHOOSE=OPTION(NOCH)
        IF(CHOOSE(1:4).NE.'Exit') THEN
          PTYPE=CHOOSE
        ELSE
          Q='e'
        ENDIF

      ELSE IF(FILE_EXT(1:5).EQ.'IPFIT'.AND.(Q.EQ.'d'.OR.Q.EQ.'r')) THEN
        OPTION( 1)='geometry'
        OPTION( 2)='fibre'
        OPTION( 3)='field'
        OPTION( 4)='Fourier'
        OPTION( 5)='Exit'
        NTCH=5
        CALL PRECHOICE1(2,IW,NOCH,NOCO,NTCH,CO,'REQUEST',OPTION,STRING,
     '    ERROR,*9999)
        CHOOSE=OPTION(NOCH)
        IF(CHOOSE(1:4).NE.'Exit') THEN
          PTYPE(1:)=CHOOSE(1:)
        ELSE
          Q='e'
        ENDIF
      ENDIF

      IF(Q.EQ.'r') THEN !find & display valid files in current directory
        CALL DISPLAY_FILE(IW,NOCO,NTFILE,CO,FILE_EXT,FILE_NAME,
     '    STRING,ERROR,*9999)
        IF(FILE_NAME(1:4).EQ.'Exit') Q='e'
      ENDIF

      CALL EXITS('PRECHOICE2')
      RETURN
 9999 CALL ERRORS('PRECHOICE2',ERROR)
      CALL EXITS('PRECHOICE2')
      RETURN 1
      END


      SUBROUTINE PRINT_IMAGE_FILL_AREAS(NOCO,NTCO,CO,NTCOQU,COQU,
     '  ERROR,*)

C**** Print fill area representation of image

      REAL XPTS(5),YPTS(5)
      INTEGER NTCOQU(*)
      CHARACTER CO(*)*(*),COQU(16,*)*(*),ERROR*(*)
      LOGICAL ABBREV,CBBREV
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:pics00.cmn'

C     CALL ENTERS('PRINT_IMAGE_FILL_AREAS',*9999)


 9998 CALL EXITS('PRINT_IMAGE_FILL_AREAS')
      RETURN
 9999 CALL ERRORS('PRINT_IMAGE_FILL_AREAS',ERROR)
      CALL EXITS('PRINT_IMAGE_FILL_AREAS')
      RETURN 1
      END


      SUBROUTINE QUIT_GRAPHICS(ERROR,*)

C**** Close any workstations and gks or phigs

      CHARACTER ERROR*(*)
      INCLUDE 'cmiss$reference:gks000.cmn'
      INCLUDE 'cmiss$reference:phig00.cmn'

      CALL ENTERS('QUIT_GRAPHICS',*9999)
      IF(GKS.OR.PHIGS) THEN
        CALL CLWS(ERROR,*9999)
      ENDIF
      IF(GKS) THEN
        GKS=.FALSE.
      ENDIF
      IF(PHIGS) THEN
        PHIGS=.FALSE.
      ENDIF

      CALL EXITS('QUIT_GRAPHICS')
      RETURN
 9999 CALL ERRORS('QUIT_GRAPHICS',ERROR)
      CALL EXITS('QUIT_GRAPHICS')
      RETURN 1
      END


      SUBROUTINE RECALL_GRAPHICS(IW,FILE_NAME,ERROR,*)

C**** Archives PHIGS structure or GKS segments.

      CHARACTER FILE_NAME*(*),ERROR*(*)

      CALL ENTERS('RECALL_GRAPHICS',*9999)

      CALL EXITS('RECALL_GRAPHICS')
      RETURN
 9999 CALL ERRORS('RECALL_GRAPHICS',ERROR)
      CALL EXITS('RECALL_GRAPHICS')
      RETURN 1
      END


!news AAY 16 Jun 91 put rotate here since it calls gl routine
      SUBROUTINE ROTATE(ISEG,NOCO,NTCO,NTCOQU,CO,COQU,CSEG,STRING,
     '  ERROR,*)

C     rotate window

      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
      INTEGER ISEG(*),NTCOQU(*)
      CHARACTER CO(*)*(*),COQU(NXCO,*)*(*),CSEG(*)*(*),ERROR*(*),
     '  STRING*(MXCH)
      LOGICAL CBBREV,ABBREV

      CALL ENTERS('ROTATE',*9999)
      IF(CO(NOCO+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        WRITE(IOOP,'(X,A)') STRING(1:IEND)
     '    //';m on IW[1]'
      ELSE
        IF(CBBREV(CO,'ON',1,NOCO+1,NTCO,N3CO)) THEN
          IW=IFROMC(CO(N3CO+1))
        ELSE
          IW=1
        ENDIF
        IF(ABBREV(COQU(NOCO,1),'M',1)) THEN
          CALL ACWK(IW,0,ERROR,*9999)
          CALL GL_ROTATE(IW)
          CALL DAWK(IW,0,ERROR,*9999)
        ENDIF
      ENDIF
      CALL EXITS('ROTATE')
      RETURN
 9999 CALL ERRORS('ROTATE',ERROR)
      CALL EXITS('ROTATE')
      RETURN 1
      END
!newe


      SUBROUTINE SET_COLOUR_LUT(COLOUR_LUT,ERROR,*)

C**** Define colour lookup table.

      REAL COLOUR_LUT(3,0:*)
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'


      CALL ENTERS('SET_COLOUR_LUT',*9999)

C **  Set up colours for colour indices 0(grey) and 1(black)
      COLOUR_LUT(1,0)=0.15
      COLOUR_LUT(2,0)=0.20
      COLOUR_LUT(3,0)=0.25
      COLOUR_LUT(1,1)=0.0
      COLOUR_LUT(2,1)=0.0
      COLOUR_LUT(3,1)=0.0

C **  Set up colours for colour indices 2(red) to 249(blue)
      DO I=2,249
        THETA=REAL(I-2)*6.0/247.0                 !THETA range
        IF(THETA.LT.2.0) THEN                     !0 - 2
          COLOUR_LUT(1,I)=1.0
          COLOUR_LUT(3,I)=0.0
          IF(THETA.LT.1.0) THEN                   !  0 - 1
            COLOUR_LUT(2,I)=THETA*0.75
          ELSE                                    !  1 - 2
            COLOUR_LUT(2,I)=0.75+(THETA-1.0)/4.0
          ENDIF
        ELSE IF(THETA.LT.4.0) THEN                !2 - 4
          COLOUR_LUT(1,I)=(4.0-THETA)/2.0
          COLOUR_LUT(2,I)=1.0
          COLOUR_LUT(3,I)=(THETA-2.0)/2.0
        ELSE                                      !4 - 6
          COLOUR_LUT(1,I)=0.0
          COLOUR_LUT(3,I)=1.0
          IF(THETA.LT.5.0) THEN                   !  4 - 5
            COLOUR_LUT(2,I)=1.0-(THETA-4.0)/4.0
          ELSE                                    !  5 - 6
            COLOUR_LUT(2,I)=0.75-(THETA-5.0)*0.75
          ENDIF
        ENDIF
      ENDDO

      CALL EXITS('SET_COLOUR_LUT')
      RETURN
 9999 CALL ERRORS('SET_COLOUR_LUT',ERROR)
      CALL EXITS('SET_COLOUR_LUT')
      RETURN 1
      END


      SUBROUTINE SET_COLOUR_LUT_RANGE(INDEX_START,INDEX_FINISH,
     '  COLOUR_LUT,RGB_START,RGB_FINISH,ERROR,*)

C**** Define colour lookup table for a given range INDEX_START to INDEX_FINISH
C**** as percentages (0-100), and for a given RGB range.

      REAL COLOUR_LUT(3,0:*),RGBSCOL(3),RGBFCOL(3),RGBSTEP(3)
      CHARACTER RGB_START*(*),RGB_FINISH*(*),CUPPER*8,RGBS*8,RGBF*8
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'

      CALL ENTERS('SET_COLOUR_LUT_RANGE',*9999)

      CALL ASSERT(INDEX_START.GE.0,'>> Invalid Range (Start < 0)',
     '  ERROR,*9999)
      CALL ASSERT(INDEX_FINISH.LE.249,'>> Invalid Range (Finish > 100)',
     '  ERROR,*9999)
      CALL ASSERT(INDEX_FINISH.GE.INDEX_START,'>> Invalid Range',
     '  ERROR,*9999)

      RGBS=CUPPER(RGB_START)
      RGBF=CUPPER(RGB_FINISH)

      IF (RGBS(1:5).EQ.'BLACK') THEN
        RGBS='000'
      ELSE IF(RGBS(1:5).EQ.'WHITE') THEN
        RGBS='111'
      ELSE IF(RGBS(1:3).EQ.'RED') THEN
        RGBS='100'
      ELSE IF(RGBS(1:5).EQ.'GREEN') THEN
        RGBS='010'
      ELSE IF(RGBS(1:4).EQ.'BLUE') THEN
        RGBS='001'
      ELSE IF(RGBS(1:4).EQ.'CYAN') THEN
        RGBS='011'
      ELSE IF(RGBS(1:6).EQ.'YELLOW') THEN
        RGBS='110'
      ELSE IF(RGBS(1:6).EQ.'PURPLE') THEN
        RGBS='101'
      ENDIF

      IF (RGBF(1:5).EQ.'BLACK') THEN
        RGBF='000'
      ELSE IF(RGBF(1:5).EQ.'WHITE') THEN
        RGBF='111'
      ELSE IF(RGBF(1:3).EQ.'RED') THEN
        RGBF='100'
      ELSE IF(RGBF(1:5).EQ.'GREEN') THEN
        RGBF='010'
      ELSE IF(RGBF(1:4).EQ.'BLUE') THEN
        RGBF='001'
      ELSE IF(RGBF(1:4).EQ.'CYAN') THEN
        RGBF='011'
      ELSE IF(RGBS(1:6).EQ.'YELLOW') THEN
        RGBS='110'
      ELSE IF(RGBS(1:6).EQ.'PURPLE') THEN
        RGBS='101'
      ENDIF

      DO I=1,3
        RGBSCOL(I)=RFROMC(RGBS(I:I))
        RGBFCOL(I)=RFROMC(RGBF(I:I))
        RGBSTEP(I)=RGBFCOL(I)-RGBSCOL(I)
      ENDDO

      ISTART=NINT(REAL(INDEX_START)*2.47+2.)
      IFINISH=NINT(REAL(INDEX_FINISH)*2.47+2.)
      IF (ISTART.EQ.IFINISH) THEN
        DO I=1,3
          COLOUR_LUT(I,ISTART)=RGBSCOL(I)
        ENDDO
      ELSE
        DO J=ISTART,IFINISH
          THETA=REAL(J-ISTART)/REAL(IFINISH-ISTART)
          DO I=1,3
            COLOUR_LUT(I,J)=RGBSCOL(I)+THETA*RGBSTEP(I)
          ENDDO
        ENDDO
      ENDIF

      CALL EXITS('SET_COLOUR_LUT_RANGE')
      RETURN
 9999 CALL ERRORS('SET_COLOUR_LUT_RANGE',ERROR)
      CALL EXITS('SET_COLOUR_LUT_RANGE')
      RETURN 1
      END


      SUBROUTINE SET_COLOUR_ONE(IW,ILUT,COLOUR,ERROR,*)

C**** Sets colour representation for index ICOL on workstation IW.

      REAL COLOUR(*)
      CHARACTER ERROR*(*)

      CALL ENTERS('SET_COLOUR_ONE',*9999)
      CALL GL_SET_COLOUR_REP(IW,ILUT,COLOUR(1),COLOUR(2),COLOUR(3))

      CALL EXITS('SET_COLOUR_ONE')
      RETURN
 9999 CALL ERRORS('SET_COLOUR_ONE',ERROR)
      CALL EXITS('SET_COLOUR_ONE')
      RETURN 1
      END


!new  took COLOUR_LUT out
      SUBROUTINE SET_COLOUR_REP(IW,ILUT_MIN,ILUT_MAX,ERROR,*)

C**** Sets colour representation for colour indices on workstation IW.
C**** ILUT=0   is grey
C**** ILUT=1   is black
C**** ILUT=ILUT_MIN..ILUT_MAX are red .. blue,
C**** ILUT_MIN is 2 and ILUT_MAX is 249 unless redefined by 'change colour'

      CHARACTER ERROR*(*)
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'

      CALL ENTERS('SET_COLOUR_REP',*9999)

      IF(DOP) THEN
        WRITE(IO4,'('' ILUT_MIN='',I4,'' ILUT_MAX='',I4)')
     '    ILUT_MIN,ILUT_MAX
      ENDIF

      DO ILUT=0,1
        ICOL=ILUT
        CALL GL_SET_COLOUR_REP(IW,ILUT,COLOUR_LUT(1,ICOL),
     '    COLOUR_LUT(2,ICOL),COLOUR_LUT(3,ICOL))
      ENDDO

      IF(COLOUR_WS) THEN
        DO ILUT=2,ILUT_MIN-1 !puts all indices up to ILUT_MIN as red
          ICOL=2
          CALL GL_SET_COLOUR_REP(IW,ILUT,COLOUR_LUT(1,ICOL),
     '      COLOUR_LUT(2,ICOL),COLOUR_LUT(3,ICOL))
        ENDDO
        DO ILUT=ILUT_MIN,ILUT_MAX
          ICOL=2+INT(REAL(ILUT-ILUT_MIN)/REAL(ILUT_MAX-ILUT_MIN)*247.)
          CALL GL_SET_COLOUR_REP(IW,ILUT,COLOUR_LUT(1,ICOL),
     '      COLOUR_LUT(2,ICOL),COLOUR_LUT(3,ICOL))
        ENDDO
        DO ILUT=ILUT_MAX+1,249 !puts all indices above ILUT_MAX as blue
          ICOL=249
          CALL GL_SET_COLOUR_REP(IW,ILUT,COLOUR_LUT(1,ICOL),
     '      COLOUR_LUT(2,ICOL),COLOUR_LUT(3,ICOL))
        ENDDO
      ENDIF

      CALL EXITS('SET_COLOUR_REP')
      RETURN
 9999 CALL ERRORS('SET_COLOUR_REP',ERROR)
      CALL EXITS('SET_COLOUR_REP')
      RETURN 1
      END


      SUBROUTINE SET_FILL_REP(IW,INDEX,STYLE,ISTYLE,ICOLOUR,ERROR,*)

C**** Resets fill area representation

      INCLUDE 'cmiss$reference:gl_constants.f'
      CHARACTER STYLE*(*),ERROR*(*)

      CALL ENTERS('SET_FILL_REP',*9999)
      IF(IW.LE.2) THEN
        IF(STYLE(1:7).EQ.'PATTERN') THEN
          CALL GL_SET_FILL_REP(IW,INDEX,GL_FILL_PATTERN,ISTYLE,ICOLOUR)
        ELSE IF(STYLE(1:5).EQ.'SOLID') THEN
          CALL GL_SET_FILL_REP(IW,INDEX,GL_FILL_SOLID,ISTYLE,ICOLOUR)
        ENDIF
      ELSE
      ENDIF

      CALL EXITS('SET_FILL_REP')
      RETURN
 9999 CALL ERRORS('SET_FILL_REP',ERROR)
      CALL EXITS('SET_FILL_REP')
      RETURN 1
      END


      SUBROUTINE SET_FILL_AREA_REP(IW,ERROR,*)

C**** Set FILL_AREA representation on workstation IW.

      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'
      INCLUDE 'cmiss$reference:gl_constants.f'

      CALL ENTERS('SET_FILL_AREA_REP',*9999)

      !Set fill indices 1 & 16 to be black and white
      CALL GL_SET_FILL_REP(IW, 1,GL_FILL_SOLID,14,1) !index 1 (black)
      CALL GL_SET_FILL_REP(IW,16,GL_FILL_SOLID,14,0) !index 16(white)
      !Set fill indices 2..15 to be patterns giving 16 shades of grey
      DO IPAT=2,15
        CALL GL_SET_FILL_REP(IW,IPAT,GL_FILL_PATTERN,12+(17-IPAT),1)
      ENDDO
      IF(COLOUR_WS) THEN !for colour use predefined colours
        DO ICOLOUR=1,233
          ILUT=INT(2.+REAL(ICOLOUR-1)*247./232.)
          INDEX=INDEX_FILL_AREA(ICOLOUR,'SOLID',' ',' ')
          CALL GL_SET_FILL_REP(IW,INDEX,GL_FILL_SOLID,1,ILUT)
        ENDDO
      ENDIF

      CALL EXITS('SET_FILL_AREA_REP')
      RETURN
 9999 CALL ERRORS('SET_FILL_AREA_REP',ERROR)
      CALL EXITS('SET_FILL_AREA_REP')
      RETURN 1
      END


      SUBROUTINE SET_GLOBAL_XFORM(IW,ISVIEW,ANGLE1,ANGLE2,ANGLE3,
     '  FIXED_PT,SCALE,SHIFT,ERROR,*)

C**** Updates global object transformation. At present one transformation is
C**** applied to all graphical segments.

      REAL A_TRANS(4,4),FIXED_PT(*),SCALE(*),SHIFT(*)
      CHARACTER ERROR*(*)

      CALL ENTERS('SET_GLOBAL_XFORM',*9999)
      CALL GL_SET_XFORM(IW,FIXED_PT,SHIFT,
     '  ANGLE1,ANGLE2,ANGLE3,SCALE,ISTATUS)

      CALL EXITS('SET_GLOBAL_XFORM')
      RETURN
 9999 CALL ERRORS('SET_GLOBAL_XFORM',ERROR)
      CALL EXITS('SET_GLOBAL_XFORM')
      RETURN 1
      END


      SUBROUTINE SET_POLYLINE_REP(IW,ERROR,*)

C**** Set POLYLINE representation on workstation IW.

      INTEGER LUT_COLOUR_INDEX(4)
      REAL LINE_WIDTH1,LINE_WIDTH2
      CHARACTER LINE_COLOUR_NAME(4)*6
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'
      INCLUDE 'cmiss$reference:gl_constants.f'

      DATA LINE_COLOUR_NAME /'RED','GREEN','BLUE','CYAN'/
      DATA LUT_COLOUR_INDEX /    2,    125,   249,   165/

      CALL ENTERS('SET_POLYLINE_REP',*9999)
      LINE_TYPE1=GL_LINE_SOLID
      LINE_TYPE2=GL_LINE_DOT
      LINE_TYPE3=GL_LINE_DASH
      LINE_TYPE4=GL_LINE_DASHDOT
      LINE_WIDTH1=1.0
      LINE_WIDTH2=4.0

      INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLACK')
      CALL GL_SET_POLYLINE_REP(IW,INDEX,LINE_TYPE1,LINE_WIDTH1,1)
      INDEX=INDEX_POLYLINE(0,'DOTTED','WIDTH1','BLACK')
      CALL GL_SET_POLYLINE_REP(IW,INDEX,LINE_TYPE2,LINE_WIDTH1,1)
      INDEX=INDEX_POLYLINE(0,'DASHED','WIDTH1','BLACK')
      CALL GL_SET_POLYLINE_REP(IW,INDEX,LINE_TYPE3,LINE_WIDTH1,1)
      INDEX=INDEX_POLYLINE(0,'DOT-DASH','WIDTH1','BLACK')
      CALL GL_SET_POLYLINE_REP(IW,INDEX,LINE_TYPE4,LINE_WIDTH1,1)

      INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH2','BLACK')
      CALL GL_SET_POLYLINE_REP(IW,INDEX,LINE_TYPE1,LINE_WIDTH2,1)
      INDEX=INDEX_POLYLINE(0,'DOTTED','WIDTH2','BLACK')
      CALL GL_SET_POLYLINE_REP(IW,INDEX,LINE_TYPE2,LINE_WIDTH2,1)
      INDEX=INDEX_POLYLINE(0,'DASHED','WIDTH2','BLACK')
      CALL GL_SET_POLYLINE_REP(IW,INDEX,LINE_TYPE3,LINE_WIDTH2,1)
      INDEX=INDEX_POLYLINE(0,'DOT-DASH','WIDTH2','BLACK')
      CALL GL_SET_POLYLINE_REP(IW,INDEX,LINE_TYPE4,LINE_WIDTH2,1)

      IF(COLOUR_WS) THEN
        DO I=1,4
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',LINE_COLOUR_NAME(I))
          CALL GL_SET_POLYLINE_REP(IW,INDEX,LINE_TYPE1,LINE_WIDTH1,
     '      LUT_COLOUR_INDEX(I))
          INDEX=INDEX_POLYLINE(0,'DOTTED','WIDTH1',LINE_COLOUR_NAME(I))
          CALL GL_SET_POLYLINE_REP(IW,INDEX,LINE_TYPE2,LINE_WIDTH1,
     '      LUT_COLOUR_INDEX(I))
          INDEX=INDEX_POLYLINE(0,'DASHED','WIDTH1',LINE_COLOUR_NAME(I))
          CALL GL_SET_POLYLINE_REP(IW,INDEX,LINE_TYPE3,LINE_WIDTH1,
     '      LUT_COLOUR_INDEX(I))
          INDEX=INDEX_POLYLINE(0,'DOT-DASH','WIDTH1',
     '      LINE_COLOUR_NAME(I))
          CALL GL_SET_POLYLINE_REP(IW,INDEX,LINE_TYPE4,LINE_WIDTH1,
     '      LUT_COLOUR_INDEX(I))

          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH2',LINE_COLOUR_NAME(I))
          CALL GL_SET_POLYLINE_REP(IW,INDEX,LINE_TYPE1,LINE_WIDTH2,
     '      LUT_COLOUR_INDEX(I))
          INDEX=INDEX_POLYLINE(0,'DOTTED','WIDTH2',LINE_COLOUR_NAME(I))
          CALL GL_SET_POLYLINE_REP(IW,INDEX,LINE_TYPE2,LINE_WIDTH2,
     '      LUT_COLOUR_INDEX(I))
          INDEX=INDEX_POLYLINE(0,'DASHED','WIDTH2',LINE_COLOUR_NAME(I))
          CALL GL_SET_POLYLINE_REP(IW,INDEX,LINE_TYPE3,LINE_WIDTH2,
     '      LUT_COLOUR_INDEX(I))
          INDEX=INDEX_POLYLINE(0,'DOT-DASH','WIDTH2',
     '      LINE_COLOUR_NAME(I))
          CALL GL_SET_POLYLINE_REP(IW,INDEX,LINE_TYPE4,LINE_WIDTH2,
     '      LUT_COLOUR_INDEX(I))
        ENDDO

        DO ICOLOUR=1,216
          ILUT=INT(2.+REAL(ICOLOUR-1)*247./215.)
          INDEX=INDEX_POLYLINE(ICOLOUR,'SOLID','WIDTH2',' ')
          CALL GL_SET_POLYLINE_REP(IW,INDEX,LINE_TYPE1,LINE_WIDTH2,ILUT)
        ENDDO
      ENDIF

      CALL EXITS('SET_POLYLINE_REP')
      RETURN
 9999 CALL ERRORS('SET_POLYLINE_REP',ERROR)
      CALL EXITS('SET_POLYLINE_REP')
      RETURN 1
      END


      SUBROUTINE SET_POLYMARKER_REP(IW,ERROR,*)

C**** Set POLYMARKER representation on workstation IW.

      INTEGER LUT_COLOUR_INDEX(4)
      REAL MARKER_SIZE
      CHARACTER MARKER_COLOUR_NAME(4)*6
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'
      INCLUDE 'cmiss$reference:gl_constants.f'

      DATA MARKER_COLOUR_NAME /'RED','GREEN','BLUE','CYAN'/
      DATA LUT_COLOUR_INDEX /    2,    125,   249,   165/

      CALL ENTERS('SET_POLYMARKER_REP',*9999)

      MARKER_SIZE1=1.0
      MARKER_SIZE2=2.0
      MARKER_TYPE1=GL_MARKER_PLUS
      MARKER_TYPE2=GL_MARKER_ASTERIX
      MARKER_TYPE3=GL_MARKER_CIRCLE
      MARKER_TYPE4=GL_MARKER_POINT
      INDEX=INDEX_POLYMARKER(0,'PLUS','SIZE1','BLACK')
      CALL GL_SET_POLYMARKER_REP(IW,INDEX,MARKER_TYPE1,MARKER_SIZE1,1)
      INDEX=INDEX_POLYMARKER(0,'ASTERISK','SIZE1','BLACK')
      CALL GL_SET_POLYMARKER_REP(IW,INDEX,MARKER_TYPE2,MARKER_SIZE1,1)
      INDEX=INDEX_POLYMARKER(0,'CIRCLE','SIZE1','BLACK')
      CALL GL_SET_POLYMARKER_REP(IW,INDEX,MARKER_TYPE3,MARKER_SIZE1,1)
      INDEX=INDEX_POLYMARKER(0,'POINT','SIZE1','BLACK')
      CALL GL_SET_POLYMARKER_REP(IW,INDEX,MARKER_TYPE4,MARKER_SIZE1,1)

      INDEX=INDEX_POLYMARKER(0,'PLUS','SIZE2','BLACK')
      CALL GL_SET_POLYMARKER_REP(IW,INDEX,MARKER_TYPE1,MARKER_SIZE2,1)
      INDEX=INDEX_POLYMARKER(0,'ASTERISK','SIZE2','BLACK')
      CALL GL_SET_POLYMARKER_REP(IW,INDEX,MARKER_TYPE2,MARKER_SIZE2,1)
      INDEX=INDEX_POLYMARKER(0,'CIRCLE','SIZE2','BLACK')
      CALL GL_SET_POLYMARKER_REP(IW,INDEX,MARKER_TYPE3,MARKER_SIZE2,1)
      INDEX=INDEX_POLYMARKER(0,'POINT','SIZE2','BLACK')
      CALL GL_SET_POLYMARKER_REP(IW,INDEX,MARKER_TYPE4,MARKER_SIZE2,1)

      IF(COLOUR_WS) THEN
        DO I=1,4
          INDEX=INDEX_POLYMARKER(0,'PLUS','SIZE1',MARKER_COLOUR_NAME(I))
          CALL GL_SET_POLYMARKER_REP(IW,INDEX,MARKER_TYPE1,
     '      MARKER_SIZE1,LUT_COLOUR_INDEX(I))
          INDEX=INDEX_POLYMARKER(0,'ASTERISK','SIZE1',
     '      MARKER_COLOUR_NAME(I))
          CALL GL_SET_POLYMARKER_REP(IW,INDEX,MARKER_TYPE2,
     '      MARKER_SIZE1,LUT_COLOUR_INDEX(I))
          INDEX=INDEX_POLYMARKER(0,'CIRCLE','SIZE1',
     '      MARKER_COLOUR_NAME(I))
          CALL GL_SET_POLYMARKER_REP(IW,INDEX,MARKER_TYPE3,
     '      MARKER_SIZE1,LUT_COLOUR_INDEX(I))
          INDEX=INDEX_POLYMARKER(0,'POINT','SIZE1',
     '      MARKER_COLOUR_NAME(I))
          CALL GL_SET_POLYMARKER_REP(IW,INDEX,MARKER_TYPE4,
     '      MARKER_SIZE1,LUT_COLOUR_INDEX(I))

          INDEX=INDEX_POLYMARKER(0,'PLUS','SIZE2',MARKER_COLOUR_NAME(I))
          CALL GL_SET_POLYMARKER_REP(IW,INDEX,MARKER_TYPE1,
     '      MARKER_SIZE2,LUT_COLOUR_INDEX(I))
          INDEX=INDEX_POLYMARKER(0,'ASTERISK','SIZE2',
     '      MARKER_COLOUR_NAME(I))
          CALL GL_SET_POLYMARKER_REP(IW,INDEX,MARKER_TYPE2,
     '      MARKER_SIZE2,LUT_COLOUR_INDEX(I))
          INDEX=INDEX_POLYMARKER(0,'CIRCLE','SIZE2',
     '      MARKER_COLOUR_NAME(I))
          CALL GL_SET_POLYMARKER_REP(IW,INDEX,MARKER_TYPE3,
     '      MARKER_SIZE2,LUT_COLOUR_INDEX(I))
          INDEX=INDEX_POLYMARKER(0,'POINT','SIZE2',
     '      MARKER_COLOUR_NAME(I))
          CALL GL_SET_POLYMARKER_REP(IW,INDEX,MARKER_TYPE4,
     '      MARKER_SIZE2,LUT_COLOUR_INDEX(I))
        ENDDO

        DO ICOLOUR=1,216
          ILUT=INT(2.+REAL(ICOLOUR-1)*247./215.)
          INDEX=INDEX_POLYMARKER(ICOLOUR,'PLUS','SIZE2',' ')
          CALL GL_SET_POLYMARKER_REP(IW,INDEX,MARKER_TYPE1,
     '      MARKER_SIZE2,ILUT)
        ENDDO
      ENDIF

      CALL EXITS('SET_POLYMARKER_REP')
      RETURN
 9999 CALL ERRORS('SET_POLYMARKER_REP',ERROR)
      CALL EXITS('SET_POLYMARKER_REP')
      RETURN 1
      END


      SUBROUTINE SET_TEXT_REP(IW,ERROR,*)

C**** Set TEXT representation on workstation IW.

      INTEGER LUT_COLOUR_INDEX(4)
      INTEGER TEXT_FONT1,TEXT_FONT2,TEXT_FONT3,TEXT_FONT4,
     '  TEXT_PRECISION1,TEXT_PRECISION2,TEXT_PRECISION3
      CHARACTER MARKER_COLOUR_NAME(4)*6
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'
      INCLUDE 'cmiss$reference:gl_constants.f'

      DATA MARKER_COLOUR_NAME /'RED','GREEN','BLUE','CYAN'/
      DATA LUT_COLOUR_INDEX /    2,    125,   249,   165/

      CALL ENTERS('SET_TEXT_REP',*9999)

      TEXT_WIDTH1=0.5
      TEXT_WIDTH2=1.0
      TEXT_FONT1=1      !Standard font in string precision
      TEXT_FONT2=-9     !Serif        font in stroke precision
      TEXT_FONT3=-11    !Serif italic   "   "    "       "
      TEXT_FONT4=-12    !Sans Serif     "   "    "       "
      TEXT_PRECISION1=GL_TEXT_STRING    !=0
      TEXT_PRECISION2=GL_TEXT_CHAR      !=1
      TEXT_PRECISION3=GL_TEXT_STROKE    !=2
      TEXT_SPACING=0.0

        INDEX=INDEX_TEXT(0,'WIDTH1','FONT1','BLACK')
        CALL GL_SET_TEXT_REP(IW,INDEX,TEXT_FONT1,TEXT_PRECISION1,
     '    TEXT_WIDTH1,TEXT_SPACING,1)
        INDEX=INDEX_TEXT(0,'WIDTH1','FONT2','BLACK')
        CALL GL_SET_TEXT_REP(IW,INDEX,TEXT_FONT2,TEXT_PRECISION3,
     '    TEXT_WIDTH1,TEXT_SPACING,1)
        INDEX=INDEX_TEXT(0,'WIDTH1','FONT3','BLACK')
        CALL GL_SET_TEXT_REP(IW,INDEX,TEXT_FONT3,TEXT_PRECISION3,
     '    TEXT_WIDTH1,TEXT_SPACING,1)
        INDEX=INDEX_TEXT(0,'WIDTH1','FONT4','BLACK')
        CALL GL_SET_TEXT_REP(IW,INDEX,TEXT_FONT4,TEXT_PRECISION3,
     '    TEXT_WIDTH1,TEXT_SPACING,1)

        INDEX=INDEX_TEXT(0,'WIDTH2','FONT1','BLACK')
        CALL GL_SET_TEXT_REP(IW,INDEX,TEXT_FONT1,TEXT_PRECISION1,
     '    TEXT_WIDTH2,TEXT_SPACING,1)
        INDEX=INDEX_TEXT(0,'WIDTH2','FONT2','BLACK')
        CALL GL_SET_TEXT_REP(IW,INDEX,TEXT_FONT2,TEXT_PRECISION3,
     '    TEXT_WIDTH2,TEXT_SPACING,1)
        INDEX=INDEX_TEXT(0,'WIDTH2','FONT3','BLACK')
        CALL GL_SET_TEXT_REP(IW,INDEX,TEXT_FONT3,TEXT_PRECISION3,
     '    TEXT_WIDTH2,TEXT_SPACING,1)
        INDEX=INDEX_TEXT(0,'WIDTH2','FONT4','BLACK')
        CALL GL_SET_TEXT_REP(IW,INDEX,TEXT_FONT4,TEXT_PRECISION3,
     '    TEXT_WIDTH2,TEXT_SPACING,1)

        IF(COLOUR_WS) THEN
          DO I=1,4
            INDEX=INDEX_TEXT(0,'WIDTH1','FONT1',MARKER_COLOUR_NAME(I))
            CALL GL_SET_TEXT_REP(IW,INDEX,TEXT_FONT1,TEXT_PRECISION1,
     '        TEXT_WIDTH1,TEXT_SPACING,LUT_COLOUR_INDEX(I))
            INDEX=INDEX_TEXT(0,'WIDTH1','FONT2',MARKER_COLOUR_NAME(I))
            CALL GL_SET_TEXT_REP(IW,INDEX,TEXT_FONT2,TEXT_PRECISION3,
     '        TEXT_WIDTH1,TEXT_SPACING,LUT_COLOUR_INDEX(I))
            INDEX=INDEX_TEXT(0,'WIDTH1','FONT3',MARKER_COLOUR_NAME(I))
            CALL GL_SET_TEXT_REP(IW,INDEX,TEXT_FONT3,TEXT_PRECISION3,
     '        TEXT_WIDTH1,TEXT_SPACING,LUT_COLOUR_INDEX(I))
            INDEX=INDEX_TEXT(0,'WIDTH1','FONT4',MARKER_COLOUR_NAME(I))
            CALL GL_SET_TEXT_REP(IW,INDEX,TEXT_FONT4,TEXT_PRECISION3,
     '        TEXT_WIDTH1,TEXT_SPACING,LUT_COLOUR_INDEX(I))

            INDEX=INDEX_TEXT(0,'WIDTH2','FONT1',MARKER_COLOUR_NAME(I))
            CALL GL_SET_TEXT_REP(IW,INDEX,TEXT_FONT1,TEXT_PRECISION1,
     '        TEXT_WIDTH2,TEXT_SPACING,LUT_COLOUR_INDEX(I))
            INDEX=INDEX_TEXT(0,'WIDTH2','FONT2',MARKER_COLOUR_NAME(I))
            CALL GL_SET_TEXT_REP(IW,INDEX,TEXT_FONT2,TEXT_PRECISION3,
     '        TEXT_WIDTH2,TEXT_SPACING,LUT_COLOUR_INDEX(I))
            INDEX=INDEX_TEXT(0,'WIDTH2','FONT3',MARKER_COLOUR_NAME(I))
            CALL GL_SET_TEXT_REP(IW,INDEX,TEXT_FONT3,TEXT_PRECISION3,
     '        TEXT_WIDTH2,TEXT_SPACING,LUT_COLOUR_INDEX(I))
            INDEX=INDEX_TEXT(0,'WIDTH2','FONT4',MARKER_COLOUR_NAME(I))
            CALL GL_SET_TEXT_REP(IW,INDEX,TEXT_FONT4,TEXT_PRECISION3,
     '        TEXT_WIDTH2,TEXT_SPACING,LUT_COLOUR_INDEX(I))
          ENDDO

          DO ICOLOUR=1,216
            ILUT=INT(2.+REAL(ICOLOUR-1)*247./215.)
            INDEX=INDEX_TEXT(ICOLOUR,'WIDTH1','FONT1',' ')
            CALL GL_SET_TEXT_REP(IW,INDEX,TEXT_FONT1,TEXT_PRECISION1,
     '        TEXT_WIDTH1,TEXT_SPACING,ILUT)
          ENDDO
        ENDIF


      CALL EXITS('SET_TEXT_REP')
      RETURN
 9999 CALL ERRORS('SET_TEXT_REP',ERROR)
      CALL EXITS('SET_TEXT_REP')
      RETURN 1
      END


      SUBROUTINE SETUP(IW,ERROR,*)

C**** Performs setup operations for workstation IW.
C**** Workstation_id=1 for z,x viewport in 3D or x,y in 2D
C****      "         2  "  z,y    "      "  "
C****      "         3  "  x,y    "      "  "
C****      "         4  "  map viewport (3D only)
C****      "         5  "  DT2651 frame-buffer 1
C****      "         6  "  DT2651 frame-buffer 2
C****      "         7  "  choice viewport
C****      "         8  "  help documentation viewport
C****      "         9  "  PHIGS 3D plot viewport
C****      "        10  "  nodal time history viewport
C****      "        11  "  sections viewport
C****      "        12  "  fibre angle distribution (3D only)
C****      "        15  "  GKS postscript file
C****      "        16  "  Phigs postscript file & .POST file
C****      "        17  "  Metafile
C****      "        18  "  x,y viewport in 3D for chords
C****      "        19  "  z,y viewport in 3D for chords
C****      "        20  "  record output
C****      "        21  "  ode output
C****      "        22  "  ode output
C****      "        23  "  ode output
C****      "        24  "  ode output
C****      "        31  "  stress/strain profile viewport
C****      "        32  "  stress/strain profile viewport
C****      "        33  "  fibre angle profile viewport
C****      "        34  "  signals
C****      "        35  "  signal trace
C****      "        36  "  trace choice viewport
C****      "        40  "  s-plane plot
C****      "        41  "  Bode plot
C****      "        42  "  Nyquist plot
C****      "        43  "  transient response plot
C****      "        44  "  ionic current model - voltage
C****      "        45  "  ionic current model - gates
C****      "        46  "  ionic current model - currents
C****      "        47  "  ionic current model - concentrations
C****      "        50  "  Matrix decomposition window
C****      "        51  "  Matrix decomposition window
C****      "        55  "  String input window
C****      "        56  "  String input window
C****      "        60  "  Spreadsheet window
C****      "        61  "  Spreadsheet GKS-plot window
C****      "        62  "  Spreadsheet GKS-plot window
C****      "        63  "  Spreadsheet GKS-plot window
C****      "        64  "  Spreadsheet GKS-plot window
C****      "        65  "  Spreadsheet GKS-plot window
C****      "        66  "  Spreadsheet GKS-plot window
C****      "        67  "  Spreadsheet PHIGS-plot window
C****      "        68  "  Spreadsheet PHIGS-plot window
C****      "        69  "  Spreadsheet GKS-plot window
C****      "        70  "  choice viewport (main spreadsheet menu)
C****      "        71  "  choice viewport (3D rotation etc)
C****      "        72  "  choice viewport (page scroll/col menu)
C****      "        73  "  choice viewport (parameter display)
C****      "        74  "  choice viewport (various)
C****      "        75  "  choice viewport (ODE choice menu)
C****      "        76  "  choice viewport (s-plane choice menu)
C****      "        77  "  choice viewport
C****      "        78  "  choice viewport
C****      "        79  "  choice viewport
C****      "        81  "  valuator
C****      "        82  "  valuator
C****      "        83  "  valuator for 1st window
C****      "        84  "  valuator  "   "     "
C****      "        85  "  valuator  "  2nd    "
C****      "        86  "  valuator  "   "     "
C****      "        87  "  valuator  "  3rd    "
C****      "        88  "  valuator  "   "     "
C****      "        89  "  valuator  "   "     "
C****      "        91  "  choice viewport (window 1 menu; row bar plot)
C****      "        92  "  choice viewport (window 2 menu; histogram)
C****      "        93  "  choice viewport (window 3 menu; line/area plot)
C****      "        94  "  choice viewport (window 4 menu; phase plane)
C****      "        95  "  choice viewport (point/bar plot)
C****      "        96  "  choice viewport (2D scatter plot)
C****      "        97  "  choice viewport (GKS_DRAW)
C****      "        98  "  choice viewport (GKS_DRAW)
C****      "        99  "  choice viewport (GKS_DRAW)

      INTEGER IFLAGS(13),
     '  TEXT_FONT1,TEXT_FONT2,TEXT_FONT3,TEXT_FONT4,
     '  TEXT_PRECISION1,TEXT_PRECISION2,TEXT_PRECISION3
      REAL AXIS(3),COLOUR(3),
     '  PLANES(2),SCALES(2),XDC(6),XNPC(6)
      CHARACTER ERROR*(*)
      STRUCTURE /lightdatatype/   !AAY 8/4/90
        INTEGER*4 MODEL
        REAL*4 COLOUR(3)
        REAL*4 VECTOR(3)
      END STRUCTURE
      RECORD/lightdatatype/LIGHTDATA
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:cbzm00.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'
      INCLUDE 'cmiss$reference:disp00.cmn'
      INCLUDE 'cmiss$reference:echo00.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:gks000.cmn'
      INCLUDE 'cmiss$reference:gks001.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'cmiss$reference:ntsg00.cmn'
      INCLUDE 'cmiss$reference:ode000.cmn'
      INCLUDE 'cmiss$reference:trac00.cmn'
      INCLUDE 'cmiss$reference:paper00.cmn'
      INCLUDE 'cmiss$reference:phig00.cmn'
      INCLUDE 'cmiss$reference:trans00.cmn'
      INCLUDE 'cmiss$reference:view00.cmn'
      INCLUDE 'cmiss$reference:gl_constants.f'

      CALL ENTERS('SETUP',*9999)

      IF(.NOT.GKS.AND..NOT.PHIGS) THEN
        CALL SET_COLOUR_LUT(COLOUR_LUT,ERROR,*9999)
        DO I=1,4
          IZOOM(I)=1
          XNDC(1,1,I)=0.
          XNDC(1,2,I)=1.
          XNDC(1,3,I)=0.
          XNDC(1,4,I)=1.
        ENDDO
C       put graphics processes in the foreground and not the background
        CALL FOREGROUND()
C       open graphics
        CALL GL_OPEN(IOER)
C       set graphics flags
        PHIGS=.TRUE.
        GKS=.TRUE.
C       Find screen size
        CALL GL_INQ_DISP(IERST,XDISP,YDISP)
        DISP=MIN(XDISP,YDISP)
        IF(DOP) WRITE(*,*) ' xdisp, ydisp =',XDISP,YDISP
        RATIO=YDISP/XDISP
C       Inquire whether colour or monochrome workstation
        CALL GL_INQ_COLOUR(IERST,NCOLOURS,ICOLOUR_FLAG,NINDICES)
        IF(DOP) THEN
          WRITE(IOOP,'('' No colours='',I8,'' ICOLOUR_FLAG='',I2,
     '      '' No predefined indices='',I4)')
     '      NCOLOURS,ICOLOUR_FLAG,NINDICES
        ENDIF
        IF(NINDICES.GE.250) THEN
          COLOUR_WS=.TRUE.
        ELSE
          COLOUR_WS=.FALSE.
        ENDIF
      ENDIF

      IF(IW.EQ.1.AND.IWKS(1).EQ.0) THEN
        CALL WKST_WINDOW(IW,0.,1.,0.,1.,ERROR,*9999)
        IF(NJT.LE.2) THEN
          CALL WKST_VIEWPORT(IW,0.01*DISP,0.01*DISP+511,0.40*DISP,
     '      0.40*DISP+511)
        ELSE
          CALL WKST_VIEWPORT(IW,0.,0.48*DISP,0.51*DISP,0.99*DISP)
        ENDIF
        CALL GL_WINDOW(IW,XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX)
        CALL SET_COLOUR_REP(IW,2,249,COLOUR_LUT,ERROR,*9999)
        CALL GL_OPEN_WKST(1)
        IWKS(IW)=1
        IWKG(IW)=1
        CALL PHIG(IW,ERROR,*9999)
        CALL SET_FILL_AREA_REP(IW,ERROR,*9999)
        CALL SET_POLYLINE_REP(IW,ERROR,*9999)
        CALL SET_POLYMARKER_REP(IW,ERROR,*9999)
        CALL SET_TEXT_REP(IW,ERROR,*9999)

      ELSE IF(IW.EQ.2.AND.IWKS(2).EQ.0) THEN
        CALL WKST_WINDOW(IW,0.,1.,0.,1.,ERROR,*9999)
        CALL WKST_VIEWPORT(IW,0.50*DISP,0.98*DISP,0.51*DISP,0.99*DISP)
        CALL GL_WINDOW(IW,XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX)
        CALL SET_COLOUR_REP(IW,2,249,COLOUR_LUT,ERROR,*9999)
        CALL GL_OPEN_WKST(2)
        IWKS(IW)=1
        IWKG(IW)=1
        CALL PHIG(IW,ERROR,*9999)
        CALL SET_FILL_AREA_REP(IW,ERROR,*9999)
        CALL SET_POLYLINE_REP(IW,ERROR,*9999)
        CALL SET_POLYMARKER_REP(IW,ERROR,*9999)
        CALL SET_TEXT_REP(IW,ERROR,*9999)

      ELSE IF(IW.EQ.3.AND.IWKS(3).EQ.0) THEN
        CALL WKST_WINDOW(IW,0.,1.,0.,1.,ERROR,*9999)
        CALL WKST_VIEWPORT(IW,0.40*DISP,0.98*DISP,0.4*DISP,0.98*DISP)
        CALL GL_WINDOW(IW,XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX)
        CALL SET_COLOUR_REP(IW,2,249,COLOUR_LUT,ERROR,*9999)
        CALL GL_OPEN_WKST(3)
        IWKS(IW)=1
        IWKG(IW)=1
        NTSG=NTSG+1
        ISVIEW=NTSG
C ***   Define transformation and view matrices
        CALL PHIG(IW,ERROR,*9999)
        CALL SET_POLYLINE_REP(IW,ERROR,*9999)
        CALL SET_POLYMARKER_REP(IW,ERROR,*9999)
        CALL SET_TEXT_REP(IW,ERROR,*9999)


C ***   set up two lights for surface shading  !AAY 8/4/90
C       first is an ambient light
        LIGHTDATA.MODEL=GL_RGB
        LIGHTDATA.COLOUR(1)=0.7
        LIGHTDATA.COLOUR(2)=0.7
        LIGHTDATA.COLOUR(3)=0.7
        CALL GL_SET_LIGHT_SRC_REP(IW,1,GL_LIGHT_TYPE_AMBIENT,
     '    4*4,LIGHTDATA)
C       second is a vector light
        LIGHTDATA.MODEL=GL_RGB
        LIGHTDATA.COLOUR(1)=0.7
        LIGHTDATA.COLOUR(2)=0.7
        LIGHTDATA.COLOUR(3)=0.7
        LIGHTDATA.VECTOR(1)=-1.0
        LIGHTDATA.VECTOR(2)=-1.0
        LIGHTDATA.VECTOR(3)=-1.0
        CALL GL_SET_LIGHT_SRC_REP(IW,2,GL_LIGHT_TYPE_INFINITE,
     '    4*7,LIGHTDATA)

C ***   set depth cue representation   !AAY 8/4/90
        PLANES(1)=0.0
        PLANES(2)=1.0
        SCALES(1)=0.0
        SCALES(2)=1.0
C       this works so long as the colour is white or black (or grey)?
        CALL GL_SET_DEPTH_CUE_REP(IW,1,GL_DEPTH_CUE_SUPPRESSED,
     '    PLANES,SCALES,GL_RGB,WHITE)
        CALL GL_SET_DEPTH_CUE_REP(IW,2,GL_DEPTH_CUE_ALLOWED,
     '    PLANES,SCALES,GL_RGB,WHITE)

        CALL GL_SET_EDGE_REP(IW,1,GL_EDGE_FLAG_ON,
     '    GL_EDGE_SOLID,1,1)
        CALL GL_SET_EDGE_REP(IW,2,GL_EDGE_FLAG_OFF,
     '    GL_EDGE_SOLID,1,1)

C ***   Set 1st extended interior representation : black hatch no lightsource
        CALL GL_SET_SURFACE_REP(IW,1,   !workstation id and index
     '    GL_FILL_HATCH,GL_FILL_HATCH, !front and back styles
     '    -1,-1,                           !front and back style indices
     '    GL_RGB,BLACK,               !rgb front colour
     '    GL_RGB,BLACK,               !rgb back colour
C         these next two don't seem to do anything?
     '    GL_SHADE_NONE,GL_SHADE_NONE, !front & back shading method
     '    GL_LIGHTING_NONE,GL_LIGHTING_NONE, !front & back lighting
     '    1.0,1.0,1.0, !ambient,diffuse and specular reflection coeffs
     '    GL_RGB,WHITE,100.0,0.0,!rgb specular colour,exponent,transparency
     '    1.0,1.0,1.0, !ambient,diffuse and specular reflection back faces
     '    GL_RGB,WHITE,100.0,0.0,!rgb specular values for back faces
     '    1,4, !surface approx type and value
     '    1,4) !trim approx type and value

C ***   2nd interior representation : white solid no lightsource
        CALL GL_SET_SURFACE_REP(IW,2,   !workstation id and index
     '    GL_FILL_SOLID,GL_FILL_SOLID, !front and back styles
     '    1,1,                             !front and back style indices
     '    GL_INDEXED_COLOUR,0,               !rgb front colour
     '    GL_INDEXED_COLOUR,0,               !rgb back colour
C         these next two don't seem to do anything?
     '    GL_SHADE_NONE,GL_SHADE_NONE, !front & back shading method
     '    GL_LIGHTING_NONE,
     '    GL_LIGHTING_NONE, !front & back lighting
     '    1.0,1.0,1.0, !ambient,diffuse and specular reflection coeffs
     '    GL_INDEXED_COLOUR,0,100.0,1.0,!rgb specular colour,exponent,transparency
     '    1.0,1.0,1.0, !ambient,diffuse and specular reflection back faces
     '    GL_INDEXED_COLOUR,0,100.0,1.0,!rgb specular values for back faces
     '    1,4, !surface approx type and value
     '    1,4) !trim approx type and value

C ***   3rd interior representation : white hollow no lightsource
        CALL GL_SET_SURFACE_REP(IW,3,   !workstation id and index
     '    GL_FILL_HOLLOW,
     '    GL_FILL_HOLLOW,         !front and back styles
     '    1,1,                             !front and back style indices
     '    GL_RGB,BLACK,               !rgb front colour
     '    GL_RGB,BLACK,               !rgb back colour
C         these next two don't seem to do anything?
     '    GL_SHADE_NONE,GL_SHADE_NONE, !front & back shading method
     '    GL_LIGHTING_NONE,
     '    GL_LIGHTING_NONE, !front & back lighting
     '    0.5,0.5,0.9, !ambient,diffuse and specular reflection coeffs
     '    GL_RGB,WHITE,100.0,1.0,!rgb specular colour,exponent,transparency
     '    0.5,0.5,0.9, !ambient,diffuse and specular reflection back faces
     '    GL_RGB,WHITE,0.9,0.4,!rgb specular values for back faces
     '    1,4, !surface approx type and value
     '    1,4) !trim approx type and value

C ***   4th interior representation : white solid with lightsource
        CALL GL_SET_SURFACE_REP(IW,4,   !workstation id and index
     '    GL_FILL_SOLID,
     '    GL_FILL_SOLID,          !front and back styles
     '    1,1,                             !front and back style indices
     '    GL_INDEXED_COLOUR,2,        !index front colour
     '    GL_INDEXED_COLOUR,2,        !index back colour
C         these next two don't seem to do anything?
     '    GL_SHADE_NONE,GL_SHADE_NONE, !front & back shading method
     '    GL_LIGHTING_DIFFUSE,
     '    GL_LIGHTING_DIFFUSE, !front & back lighting
     '    0.3,0.5,0.9, !ambient,diffuse and specular reflection coeffs
     '    GL_INDEXED_COLOUR,2,100.0,1.0,!rgb specular colour,exponent,transparency
     '    1.0,1.0,1.0, !ambient,diffuse and specular reflection back faces
     '    GL_INDEXED_COLOUR,2,100.0,1.0,!rgb specular values for back faces
     '    1,4, !surface approx type and value
     '    1,4) !trim approx type and value

C ***   5th interior representation : green/blue solid shaded
        CALL GL_SET_SURFACE_REP(IW,5,   !workstation id and index
     '    GL_FILL_SOLID,
     '    GL_FILL_SOLID,          !front and back styles
     '    1,1,                             !front and back style indices
     '    GL_RGB,GREEN,               !rgb front colour
     '    GL_RGB,BLUE,               !rgb back colour
C         these next two don't seem to do anything?
     '    GL_SHADE_NONE,GL_SHADE_NONE, !front & back shading method
     '    GL_LIGHTING_DIFFUSE,
     '    GL_LIGHTING_SPECULAR, !front & back lighting
     '    0.5,0.5,0.9, !ambient,diffuse and specular reflection coeffs
     '    GL_RGB,WHITE,100.0,1.0,!rgb specular colour,exponent,transparency
     '    1.0,0.5,0.9, !ambient,diffuse and specular reflection back faces
     '    GL_RGB,BLUE,0.9,0.4,!rgb specular values for back faces
     '    1,4, !surface approx type and value
     '    1,4) !trim approx type and value



      ELSE IF(IW.EQ.4.AND.IWKS(4).EQ.0) THEN
        CALL SET_COLOUR_REP(IW,2,249,COLOUR_LUT,ERROR,*9999)
        CALL SET_FILL_AREA_REP(IW,ERROR,*9999)
        CALL SET_POLYLINE_REP(IW,ERROR,*9999)
        CALL SET_POLYMARKER_REP(IW,ERROR,*9999)
        CALL SET_TEXT_REP(IW,ERROR,*9999)
        CALL WKST_WINDOW(IW,0.,1.,0.,1.,ERROR,*9999)
        CALL GL_WINDOW(IW,-1.0,1.0,-1.0,1.0,-1.0,1.0)
        CALL WKST_VIEWPORT(IW,0.1*DISP,0.9*DISP,0.2*DISP,1.0*DISP)
        CALL GL_OPEN_WKST(IW)
  CALL GL_INIT_VIEW(IW)
        IWKS(IW)=1

      ELSE IF(IW.EQ.5.AND.IWKS(5).EQ.0) THEN
C       At some time this window on the frame buffer should be opened thru GKS
C       open(unit=105,file='fem.w5',status='unknown',dispose='delete')
C       CALL GOPWK(5,105,GWSDEF)
C       CALL GSWN(IW,0.0,1.0,0.0,1.0)
        IWKS(IW)=1

      ELSE IF(IW.EQ.6.AND.IWKS(6).EQ.0) THEN
C       At some time this window on the frame buffer should be opened thru GKS
C       open(unit=106,file='fem.w6',status='unknown',dispose='delete')
C       CALL GOPWK(6,106,GWSDEF)
C       CALL GSWN(IW,0.0,1.0,0.0,1.0)
        IWKS(IW)=1

      ELSE IF(IW.EQ.7.AND.IWKS(7).EQ.0) THEN
        CALL SET_COLOUR_REP(IW,2,249,COLOUR_LUT,ERROR,*9999)
        CALL SET_FILL_AREA_REP(IW,ERROR,*9999)
        CALL SET_POLYLINE_REP(IW,ERROR,*9999)
        CALL SET_POLYMARKER_REP(IW,ERROR,*9999)
        CALL SET_TEXT_REP(IW,ERROR,*9999)
        CALL WKST_WINDOW(IW,0.,1.,0.,1.,ERROR,*9999)
        CALL GL_WINDOW(IW,-1.0,1.0,-1.0,1.0,-1.0,1.0)
        CALL WKST_VIEWPORT(IW,0.,XDISP,0.,YDISP)
        CALL GL_OPEN_WKST(IW)
        IWKS(IW)=1

      ELSE IF(IW.EQ.8.AND.IWKS(8).EQ.0) THEN
        CALL SET_COLOUR_REP(IW,2,249,COLOUR_LUT,ERROR,*9999)
        CALL SET_FILL_AREA_REP(IW,ERROR,*9999)
        CALL SET_POLYLINE_REP(IW,ERROR,*9999)
        CALL SET_POLYMARKER_REP(IW,ERROR,*9999)
        CALL SET_TEXT_REP(IW,ERROR,*9999)
        CALL WKST_WINDOW(IW,0.,1.,0.3688,1.,ERROR,*9999)
        CALL GL_WINDOW(IW,-1.0,1.0,-1.0,1.0,-1.0,1.0)
        CALL WKST_VIEWPORT(IW,0.2*DISP,1.0*DISP,0.295*DISP,0.8*DISP)
        CALL GL_OPEN_WKST(IW)
        IWKS(IW)=1

      ELSE IF(IW.EQ.9.AND.IWKS(9).EQ.0) THEN

      ELSE IF(IW.EQ.10.AND.IWKS(10).EQ.0) THEN
        CALL SET_COLOUR_REP(IW,2,249,COLOUR_LUT,ERROR,*9999)
        CALL SET_FILL_AREA_REP(IW,ERROR,*9999)
        CALL SET_POLYLINE_REP(IW,ERROR,*9999)
        CALL SET_POLYMARKER_REP(IW,ERROR,*9999)
        CALL SET_TEXT_REP(IW,ERROR,*9999)
        CALL WKST_WINDOW(IW,0.,1.,0.5-0.23/0.49/2.,0.5+0.23/0.49/2.,
     '    ERROR,*9999)
        CALL GL_WINDOW(IW,-0.14,1.1,-2.3,2.3,-1.0,1.0)
        CALL WKST_VIEWPORT(IW,0.5*XDISP,0.99*XDISP,YDISP-0.23*XDISP,
     '    YDISP)
        CALL GL_OPEN_WKST(IW)
        IWKS(IW)=1

      ELSE IF(IW.EQ.11.AND.IWKS(11).EQ.0) THEN
        CALL SET_COLOUR_REP(IW,2,249,COLOUR_LUT,ERROR,*9999)
        CALL SET_FILL_AREA_REP(IW,ERROR,*9999)
        CALL SET_POLYLINE_REP(IW,ERROR,*9999)
        CALL SET_POLYMARKER_REP(IW,ERROR,*9999)
        CALL SET_TEXT_REP(IW,ERROR,*9999)
        CALL WKST_WINDOW(IW,0.,1.,0.5-0.23/0.49/2.,0.5+0.23/0.49/2.,
     '    ERROR,*9999)
        CALL GL_WINDOW(IW,-0.14,1.1,-2.3,2.3,-1.0,1.0)
        CALL WKST_VIEWPORT(IW,0.5*XDISP,0.99*XDISP,YDISP-0.23*XDISP,
     '    YDISP-0.26*XDISP)
        CALL GL_OPEN_WKST(IW)
        IWKS(IW)=1

      ELSE IF(IW.EQ.12.AND.IWKS(12).EQ.0) THEN
        CALL SET_COLOUR_REP(IW,2,249,COLOUR_LUT,ERROR,*9999)
        CALL SET_FILL_AREA_REP(IW,ERROR,*9999)
        CALL SET_POLYLINE_REP(IW,ERROR,*9999)
        CALL SET_POLYMARKER_REP(IW,ERROR,*9999)
        CALL SET_TEXT_REP(IW,ERROR,*9999)
        CALL WKST_WINDOW(IW,0.,1.,0.,1.,ERROR,*9999)
        CALL GL_WINDOW(IW,-0.20,1.1,-1.1*PI/2.,1.1*PI/2.,-1.0,1.0)
        CALL WKST_VIEWPORT(IW,0.96*XDISP,0.99*XDISP,YDISP-0.49*XDISP,
     '    YDISP-0.26*XDISP)
        CALL GL_OPEN_WKST(IW)
        IWKS(IW)=1

      ELSE IF(IW.EQ.15.AND.IWKS(15).EQ.0) THEN


      ELSE IF(IW.EQ.16.AND.IWKS(16).EQ.0) THEN

      ELSE IF(IW.EQ.17.AND.IWKS(17).EQ.0) THEN

      ELSE IF(IW.EQ.18.AND.IWKS(18).EQ.0) THEN

      ELSE IF(IW.EQ.19.AND.IWKS(19).EQ.0) THEN

      ELSE IF(IW.EQ.20.AND.IWKS(20).EQ.0) THEN

      ELSE IF(IW.EQ.21.AND.IWKS(21).EQ.0) THEN

      ELSE IF(IW.EQ.22.AND.IWKS(22).EQ.0) THEN

      ELSE IF(IW.EQ.23.AND.IWKS(23).EQ.0) THEN

      ELSE IF(IW.EQ.31.AND.IWKS(31).EQ.0) THEN
        CALL WKST_WINDOW(IW,0.,1.,0.,1.,ERROR,*9999)
        CALL WKST_VIEWPORT(IW,0.5*XDISP,0.99*XDISP,YDISP-0.23*XDISP,
     '    YDISP)
        CALL GL_WINDOW(IW,-0.1,1.2,-1.3,1.3,-1.0,1.0)
        CALL SET_COLOUR_REP(IW,2,249,COLOUR_LUT,ERROR,*9999)
        CALL GL_OPEN_WKST(IW)
        IWKS(IW)=1
Cold    CALL PHIG(IW,ERROR,*9999)
!news   eventually replace all PHIG calls with GL_INIT_VIEW AAY 28 May 91
  CALL GL_INIT_VIEW(IW)
!newe
        CALL SET_FILL_AREA_REP(IW,ERROR,*9999)
        CALL SET_POLYLINE_REP(IW,ERROR,*9999)
        CALL SET_POLYMARKER_REP(IW,ERROR,*9999)
        CALL SET_TEXT_REP(IW,ERROR,*9999)

      ELSE IF(IW.EQ.32.AND.IWKS(32).EQ.0) THEN
        CALL WKST_WINDOW(IW,0.,1.,0.,1.,ERROR,*9999)
        CALL WKST_VIEWPORT(IW,0.5*XDISP,0.99*XDISP,YDISP-0.49*XDISP,
     '     YDISP-0.26*XDISP)
        CALL GL_WINDOW(IW,-0.1,1.2,-1.3,1.3,-1.0,1.0)
        CALL SET_COLOUR_REP(IW,2,249,COLOUR_LUT,ERROR,*9999)
        CALL GL_OPEN_WKST(IW)
        IWKS(IW)=1
Cold    CALL PHIG(IW,ERROR,*9999)
!news   eventually replace all PHIG calls with GL_INIT_VIEW AAY 28 May 91
  CALL GL_INIT_VIEW(IW)
!newe
        CALL SET_FILL_AREA_REP(IW,ERROR,*9999)
        CALL SET_POLYLINE_REP(IW,ERROR,*9999)
        CALL SET_POLYMARKER_REP(IW,ERROR,*9999)
        CALL SET_TEXT_REP(IW,ERROR,*9999)

      ELSE IF(IW.EQ.33.AND.IWKS(33).EQ.0) THEN

      ELSE IF(IW.EQ.34.AND.IWKS(34).EQ.0) THEN

      ELSE IF(IW.EQ.35.AND.IWKS(35).EQ.0) THEN

      ELSE IF(IW.EQ.36.AND.IWKS(36).EQ.0) THEN

      ELSE IF(IW.EQ.40.AND.IWKS(40).EQ.0) THEN

      ELSE IF(IW.EQ.41.AND.IWKS(41).EQ.0) THEN

      ELSE IF(IW.EQ.42.AND.IWKS(42).EQ.0) THEN

      ELSE IF(IW.EQ.43.AND.IWKS(43).EQ.0) THEN

      ELSE IF(IW.EQ.44.AND.IWKS(44).EQ.0) THEN

      ELSE IF(IW.EQ.45.AND.IWKS(45).EQ.0) THEN

      ELSE IF(IW.EQ.46.AND.IWKS(46).EQ.0) THEN

      ELSE IF(IW.EQ.47.AND.IWKS(47).EQ.0) THEN

      ELSE IF(IW.EQ.51.AND.IWKS(51).EQ.0) THEN

      ELSE IF(IW.EQ.55.AND.IWKS(55).EQ.0) THEN

      ELSE IF(IW.EQ.56.AND.IWKS(56).EQ.0) THEN

      ELSE IF(IW.EQ.60.AND.IWKS(60).EQ.0) THEN

      ELSE IF(IW.EQ.61.AND.IWKS(61).EQ.0) THEN

      ELSE IF(IW.EQ.62.AND.IWKS(62).EQ.0) THEN

      ELSE IF(IW.EQ.63.AND.IWKS(63).EQ.0) THEN

      ELSE IF(IW.EQ.64.AND.IWKS(64).EQ.0) THEN

      ELSE IF(IW.EQ.65.AND.IWKS(65).EQ.0) THEN

      ELSE IF(IW.EQ.66.AND.IWKS(66).EQ.0) THEN

      ELSE IF(IW.EQ.67.AND.IWKS(67).EQ.0) THEN

      ELSE IF(IW.EQ.68.AND.IWKS(68).EQ.0) THEN

      ELSE IF(IW.EQ.69.AND.IWKS(69).EQ.0) THEN

      ELSE IF(IW.EQ.70.AND.IWKS(70).EQ.0) THEN

      ELSE IF(IW.EQ.71.AND.IWKS(71).EQ.0) THEN

      ELSE IF(IW.EQ.72.AND.IWKS(72).EQ.0) THEN

      ELSE IF(IW.EQ.73.AND.IWKS(73).EQ.0) THEN

      ELSE IF(IW.EQ.74.AND.IWKS(74).EQ.0) THEN

      ELSE IF(IW.EQ.75.AND.IWKS(75).EQ.0) THEN

      ELSE IF(IW.EQ.76.AND.IWKS(76).EQ.0) THEN

      ELSE IF(IW.EQ.77.AND.IWKS(77).EQ.0) THEN

      ELSE IF(IW.EQ.78.AND.IWKS(78).EQ.0) THEN

      ELSE IF(IW.EQ.79.AND.IWKS(79).EQ.0) THEN

      ELSE IF(IW.EQ.81.AND.IWKS(81).EQ.0) THEN
        IWKS(IW)=1

      ELSE IF(IW.EQ.82.AND.IWKS(82).EQ.0) THEN
        IWKS(IW)=1

      ELSE IF(IW.EQ.83.AND.IWKS(83).EQ.0) THEN
        IWKS(IW)=1

      ELSE IF(IW.EQ.84.AND.IWKS(84).EQ.0) THEN
        IWKS(IW)=1

      ELSE IF(IW.EQ.85.AND.IWKS(85).EQ.0) THEN
        IWKS(IW)=1

      ELSE IF(IW.EQ.86.AND.IWKS(86).EQ.0) THEN
        IWKS(IW)=1

      ELSE IF(IW.EQ.87.AND.IWKS(87).EQ.0) THEN
        IWKS(IW)=1

      ELSE IF(IW.EQ.88.AND.IWKS(88).EQ.0) THEN
        IWKS(IW)=1

      ELSE IF(IW.EQ.89.AND.IWKS(89).EQ.0) THEN
        IWKS(IW)=1

      ELSE IF(IW.EQ.91.AND.IWKS(91).EQ.0) THEN
        IWKS(IW)=1

      ELSE IF(IW.EQ.92.AND.IWKS(92).EQ.0) THEN
        IWKS(IW)=1

      ELSE IF(IW.EQ.93.AND.IWKS(93).EQ.0) THEN
        IWKS(IW)=1

      ELSE IF(IW.EQ.94.AND.IWKS(94).EQ.0) THEN
        IWKS(IW)=1

      ELSE IF(IW.EQ.95.AND.IWKS(95).EQ.0) THEN

      ELSE IF(IW.EQ.96.AND.IWKS(96).EQ.0) THEN

      ELSE IF(IW.EQ.97.AND.IWKS(97).EQ.0) THEN

      ELSE IF(IW.EQ.98.AND.IWKS(98).EQ.0) THEN

      ELSE IF(IW.EQ.99.AND.IWKS(99).EQ.0) THEN

      ENDIF

      IF(DOP) THEN
        WRITE(OP_STRING,'('' Reset to:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' IW='',I3,'' IWKS(iw)='',I2,'
     '    //''' IWKT(iw)='',I2,'' IWKG(iw)='',I2)')
     '    IW,IWKS(IW),IWKT(IW),IWKG(IW)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      NOIW=0
      DO IIW=1,99      !to record defined workstations in IWKDEF
        IF(IWKS(IIW).GT.0) THEN !IIW is defined
          NOIW=NOIW+1
          IWKDEF(NOIW)=IIW
        ENDIF
      ENDDO
      IWKDEF(0)=NOIW

      CALL EXITS('SETUP')
      RETURN
 9999 CALL ERRORS('SETUP',ERROR)
      CALL EXITS('SETUP')
      RETURN 1
      END


      SUBROUTINE SET_VIEW(ISTATUS,IW,MODE_PROJ,NPC_CLIP,
     '  ANGLE1,ANGLE2,ANGLE3,
     '  BACK_PLANE_DIST,FRONT_PLANE_DIST,
     '  FPT,FSCALE,FSHFT,PROJ_REF_PT,VIEW_PLANE_DIST,
     '  VIEW_PLANE,VIEW_REF_PT,VIEW_UP,
     '  VIEWPORT,WINDOW,ERROR,*)

C**** Updates viewing transform, region of view and clipping planes.
C     MODE_PROJ is 0 for orthographic projection
C                  1 for perspective projection.
C     NPC_CLIP defines clipping planes (usually not changed).
C     VIEW_REF_PT is the point in the object at which you look.
C     PROJ_REF_PT is the point from which you look (usually 0,0,3*DIAG)
C     VIEW_PLANE defines the normal to the plane at which you look
C       (usually 0,0,1).
C     VIEW_UP defines the orientation of the viewer ie which way is up
C       (usually 0,1,0).
C     VIEW_PLANE_DIST I forget what this is (usually 0).
C     BACK_PLANE_DIST distance from PROJ_REF_PT to back clipping plane
C     FRONT_PLANE_DIST distance from PROJ_REF_PT to front clipping plane
C     VIEWPORT area of workstation window in which the object will appear.
C       (usually not changed from default)
C     WINDOW area of the view plane which gets mapped onto the viewport.
C       (usually from XMIN-VIEW_REF_PT(1),YMIN-VIEW_REF_PT(2) to
C       XMAX-VIEW_REF_PT(1),YMAX-VIEW_REF_PT(2) )
C     ANGLE1,ANGLE2,ANGLE3,FPT,FSCALE,FSHFT define a transformation which will
C       operate on the view. This is the way rotations, translations etc
C       of the view are achieved, rather than by changing PROJ_REF_PT,
C       VIEW_PLANE,VIEW_UP, etc separately.

      INTEGER NPC_CLIP(*)
      REAL A_MAP(4,4),A_ORIENT(4,4),B_ORIENT(4,4),A_TRANS(4,4),
     '  B_MATRIX(4,4),FPT(*),FSCALE(*),FSHFT(*),
     '  PROJ_REF_PT(*),VIEW_PLANE(*),VIEW_REF_PT(*),VIEW_UP(*),
     '  VIEWPORT(*),WINDOW(*)
      CHARACTER ERROR*(*)

C     there must be a better way of doing this AAY
      CALL ENTERS('SET_VIEW',*9999)

      CALL GL_SET_VIEW(ISTATUS,IW,MODE_PROJ,NPC_CLIP,
     '  ANGLE1,ANGLE2,ANGLE3,
     '  BACK_PLANE_DIST,FRONT_PLANE_DIST,
     '  FPT,FSCALE,FSHFT,PROJ_REF_PT,VIEW_PLANE_DIST,
     '  VIEW_PLANE,VIEW_REF_PT,VIEW_UP,
     '  VIEWPORT,WINDOW)

Cnew  I don't like any of the BUILD_XFORM or EVAL_VIEW or COMPOSE_MATRIX stuff
Cnew  let above call handle it all. AAY May 1991
Csun  CALL GL_BUILD_XFORM_MATRIX3(FPT,FSHFT,
Csun '  ANGLE1,ANGLE2,ANGLE3,FSCALE,ISTATUS,B_MATRIX)
Csun  CALL GL_EVAL_VIEW_ORIEN_MATRIX3(VIEW_REF_PT,
Csun    VIEW_PLANE,VIEW_UP,ISTATUS,A_ORIENT)
Csun  CALL GL_COMPOSE_MATRIX3(B_MATRIX,A_ORIENT,ISTATUS,
Csun '  B_ORIENT)
Csun  IF(ISTATUS.NE.0)THEN
Csun    WRITE(*,*)' ---Error from eval_view_orient_matrix=',ISTATUS
Csun  ELSE
Csun    CALL GL_EVAL_VIEW_MAP_MATRIX3(WINDOW,VIEWPORT,
Csun '    MODE_PROJ,PROJ_REF_PT,VIEW_PLANE_DIST,
Csun '    BACK_PLANE_DIST,FRONT_PLANE_DIST,ISTATUS,A_MAP)
Csun    IF(ISTATUS.NE.0)THEN
Csun      WRITE(*,*)' ---Error from eval_view_map_matrix=',ISTATUS
Csun    ELSE
Csun      CALL GL_SET_VIEW_REP3(IW,1,B_ORIENT,A_MAP,NPC_CLIP)
Csun    ENDIF
Csun  ENDIF


      CALL EXITS('SET_VIEW')
      RETURN
 9999 CALL ERRORS('SET_VIEW',ERROR)
      CALL EXITS('SET_VIEW')
      RETURN 1
      END


      SUBROUTINE STROKE(INIT,IW,INSTAT,MODE,NECHO,NTMPTS,NTPTS,XPTS,
     '  YPTS,ERROR,*)

C**** Calls GKS stroke
C**** INIT specifies whether stroke is to be initialised (1:yes,2:no)
C**** IW specifies workstation number (which is also transformation no)
C**** INSTAT is returned as 1 if locate is successful, 0 otherwise
C**** MODE is input mode - 'REQUEST','SAMPLE' or 'EVENT'
C**** NECHO is 0 on entry for no echo
C**** NECHO is 3 on entry for default   prompt/echo
C**** NECHO is 4 on entry for line      prompt/echo
C**** NECHO is 5 on entry for rectangle prompt/echo
C**** NECHO is 6 on entry for digital   prompt/echo
C**** XPTS,YPTS are returned world coords of stroke

      INTEGER IARRAY(2)
      REAL XPTS(*),YPTS(*),RDATA(2,100)
      CHARACTER ERROR*(*),MODE*(*),CLASS*8,RECORD(2)*80
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'

      DATA LD1/1/

      CALL ENTERS('STROKE',*9999)
      IF(INIT.EQ.1) THEN !initialise stroke device
        CALL GL_INIT_STROKE(IW,LD1)
      ENDIF
      CALL GL_SELECT_XFORM(IW)
      IF(MODE(1:7).EQ.'REQUEST') THEN
        CALL EVENT(ID_WS,ID_DEVICE,INSTAT,CLASS,IDATA,RDATA,SDATA,
     '    ERROR,*9999)
C       if its a stroke
        IF(INSTAT.EQ.1.AND.CLASS(1:8).EQ.'STROKE'.AND.ID_WS.EQ.IW)THEN
C         return cursor location
          NTPTS=IDATA
          DO I=1,NTPTS
            XPTS(I)=RDATA(1,I)
            YPTS(I)=RDATA(2,I)
          ENDDO
        ELSE
        ENDIF
C       delete stroke
  CALL GL_DELETE_STROKE(IW)
      ENDIF

      CALL EXITS('STROKE')
      RETURN

 9999 CALL ERRORS('STROKE',ERROR)
      CALL EXITS('STROKE')
      RETURN 1
      END


      SUBROUTINE SURFACE(IBUNDLE,IW,NVERTICES,VERTICES,FACET_COLOURS,
     '  FACET_NORMAL,VERTEX_COLOURS,VERTEX_NORMAL,NFACETS,NVPF,
     '  IVERTICES,ERROR,*)

C**** draws a surface on IW with index IBUNDLE. VERTICES(1..3,np) is an array
C**** containing the 3D coordinates of each vertex. If the IBUNDLE is 0 the
C**** primative will use the previously defined surface index.
C**** FACET_COLOURS  array of colour indices for each facet
C**** FACET_NORMAL   array of facet normals
C**** VERTEX_COLOURS array of vertex colour indices for each facet
C**** VERTEX_NORMAL  array of vertex normals
C**** NFACETS        number of facets
C**** NVPF           array containing #vertices for each facet
C**** IVERTICES      array of vertex indices for each facet

      INTEGER ACTIVE(2),DEACT(2),IVERTICES(3,*),FACET_COLOURS(*),
     '  NVPF(*),VERTEX_COLOURS(*),IEDGEVIS(600)
      REAL FACET_NORMAL(3,*),VERTICES(3,*),
     '  VERTEX_NORMAL(3,*),Z(3)
      CHARACTER ERROR*(*)
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:surf00.cmn'
      INCLUDE 'cmiss$reference:gl_constants.f'

      CALL ENTERS('SURFACE',*9999)
      IF(NVERTICES.EQ.0)GOTO 9998

C ***   set depth cue representation
        IF(DEPTHCUE(1:2).EQ.'ON') THEN
          CALL GL_SET_DEPTH_CUE_INDEX(2)
        ELSE
          CALL GL_SET_DEPTH_CUE_INDEX(1)
        ENDIF
C ***   set edge representation
        IF(EDGES(1:2).EQ.'ON') THEN
          IEDGEVIS_FLAG=GL_EDGE_VISIBILE
          CALL GL_SET_EDGE_INDEX(1)
        ELSE
          IEDGEVIS_FLAG=GL_EDGE_INVISIBILE
          CALL GL_SET_EDGE_INDEX(2)
        ENDIF
C ***   face culling
        IF(CULL(1:4).EQ.'NONE') THEN
          CALL GL_SET_FACE_DIST_MODE(GL_DIST_NO)
          CALL GL_SET_FACE_CULL_MODE(GL_CULL_NONE)
        ELSE IF(CULL(1:4).EQ.'BACK') THEN
          CALL GL_SET_FACE_DIST_MODE(GL_DIST_YES)
          CALL GL_SET_FACE_CULL_MODE(GL_CULL_BACK)
        ELSE IF(CULL(1:5).EQ.'FRONT') THEN
          CALL GL_SET_FACE_DIST_MODE(GL_DIST_YES)
          CALL GL_SET_FACE_CULL_MODE(GL_CULL_FRONT)
        ENDIF
C ***   activate the two lights
        IF(INDEX_SURF.EQ.4.OR.INDEX_SURF.EQ.5)THEN
          ACTIVE(1)=1
          ACTIVE(2)=2
          CALL GL_SET_LIGHT_SRC_STATE(2,ACTIVE,0,DEACT)
        ENDIF
C ***   colour indexing
        IF(VARIABLE_TYPE(1:8).EQ.'GEOMETRY') THEN
          INDEX_COLOUR=GL_PFA_NONE
          ICOLOUR_TYPE=GL_RGB
          IF(INDEX_SURF.EQ.2)THEN !white fill area
            INDEX_COLOUR=GL_PFA_COLOUR
            ICOLOUR_TYPE=GL_INDEXED_COLOUR
            DO NF=1,NFACETS
              FACET_COLOURS(NF)=0
            ENDDO
          ENDIF
        ELSE IF(VARIABLE_TYPE(1:5).EQ.'FIELD') THEN
          INDEX_COLOUR=GL_PFA_COLOUR
          ICOLOUR_TYPE=GL_INDEXED_COLOUR
        ELSE IF(VARIABLE_TYPE(1:6).EQ.'STRAIN') THEN
          INDEX_COLOUR=GL_PV_COLOUR
          ICOLOUR_TYPE=GL_INDEXED_COLOUR
        ENDIF

        CALL GL_SURFACE(
     '    IW,
     '    GL_SHAPE_CONVEX,       !convex facets
     '    GL_PFA_NONE,          !per facet data flag
     '    INDEX_COLOUR,          ! per vertex data flag
     '    IEDGEVIS_FLAG,         !per edge data flag
     '    IEDGEVIS,               !array of edge visibility flags
     '    ICOLOUR_TYPE,          !colour type
     '    NVERTICES,             !number of vertices
     '    FACET_NORMAL,          !array of facet normals
     '    FACET_COLOURS,         !facet colour values
     '    VERTICES,              !array of vertices
     '    VERTEX_COLOURS,        !array of vertex colours
     '    VERTEX_NORMAL,         !array of vertex normals
     '    NFACETS,               !number of facets
     '    NVPF,                  !list of #vertices in each facet
     '    IVERTICES)             !list of vertex indices for each facet



 9998 CALL EXITS('SURFACE')
      RETURN
 9999 CALL ERRORS('SURFACE',ERROR)
      CALL EXITS('SURFACE')
      RETURN 1
      END


      SUBROUTINE TEXT(IBUNDLE,IGEOM,IW,STRING,POINT,ERROR,*)

C**** Draws text on IW with non-geometric attributes bundled in IBUNDLE, and
C**** geometric attributes according to IGEOM.
C**** STRING contains the text to be shown at POINT(nj)
C**** If IBUNDLE or IGEOM is 0 the primitive will use the previously
C**** defined text index.

      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'cmiss$reference:map000.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
      REAL POINT(*),TVECS(2,3),X(3)
      CHARACTER STRING*(MXCH)

      CALL ENTERS('TEXT',*9999)
      IF(IBUNDLE.NE.0) CALL GL_SET_TEXT_INDEX(IW,IBUNDLE)
      IF(IGEOM.NE.0) CALL GL_SET_TEXT_GEOM(IGEOM)
      IF(IW.LE.3) THEN
        CALL GL_TEXT(POINT,STRING,LEN(STRING))

      ELSE IF(IW.EQ.4) THEN
        IF(PROJEC(1:11).EQ.'RECTANGULAR') THEN !points x and y
        ELSE IF(PROJEC(1:2).EQ.'XI') THEN !assume points are in Xi coords
          POINT(1)=-1.0+2.*(REAL(MXI1-1)+POINT(1))/MAX_XI
          POINT(2)=-1.0+2.*(REAL(MXI2-1)+POINT(2))/MAX_XI
        ELSE !plot Y against Z - assume has been transformed to polar coords
          CALL MAP4(1,POINT,ERROR,*9999)
        ENDIF
        CALL GL_TEXT(POINT,STRING,LEN(STRING))

      ELSE
        CALL GL_TEXT(POINT,STRING,LEN(STRING))
      ENDIF


      CALL EXITS('TEXT')
      RETURN
 9999 CALL ERRORS('TEXT',ERROR)
      CALL EXITS('TEXT')
      RETURN 1
      END


!news
      SUBROUTINE TRANSFORM_SEGMENT(ISEGM_LIST,IW,TYPE,
     '  XCENTRE,YCENTRE,VALUE,ERROR,*)

C**** Transforms graphics segments in ISEGM_LIST(isegm),isegm=1,ISEGM_LIST(0)
C**** on IW with VALUE for:
C**** TYPE is 'rotate','scale','x-translate' or 'y-translate'.
C**** Rotate and scale operate about XCENTRE,YCENTRE.

      INTEGER ISEGM_LIST(0:*),IW
      REAL XCENTRE,YCENTRE,VALUE
      CHARACTER ERROR*(*),TYPE*(*)
      CALL ENTERS('TRANSFORM_SEGMENT',*9999)
      CALL EXITS('TRANSFORM_SEGMENT')
      RETURN
 9999 CALL ERRORS('TRANSFORM_SEGMENT',ERROR)
      CALL EXITS('TRANSFORM_SEGMENT')
      RETURN 1
      END
!newe


      SUBROUTINE VALUATOR(LABEL,IW,MODE,NTYPE,VALMIN,VALMAX,VALINI,
     '  VALUE,XREF,YREF,ERROR,*)

C**** Calls GKS valuator
C**** MODE can be 'REQUEST' or 'EVENT' or 'SAMPLE'
C**** VALINI is initial value on entry
C**** VALUE  is returned value on exit if in request mode
C**** XREF,YREF are coordinates of reference point
C**** NTYPE=1  is for reference point at bottom LH corner (screen coords)
C**** NTYPE=2  is  "      "       "    " bottom RH corner    "      "
C**** NTYPE=3  is  "      "       "    " top    LH corner (world  coords)
C**** NTYPE=4  is  "      "       "    " top    RH corner    "      "
C**** NTYPE=5  is  "      "       "    " top    LH corner (x:0..1,y:-1..1)
C**** NTYPE=6  is  "      "       "    " top    RH corner    "      "
C**** NTYPE=9  is  "      "       "    " top    LH corner (screen coords)
C**** NTYPE=10 is  "      "       "    " top    RH corner    "      "

      CHARACTER ERROR*(*),MODE*(*),CLASS*8,DUMMY,LABEL*20
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:disp00.cmn'
      INCLUDE 'cmiss$reference:echo00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:gks000.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'

      DATA LD1/1/,XFACTOR/0.15/,YFACTOR/0.1/,DUMMY/' '/

      CALL ENTERS('VALUATOR',*9999)
      IF(INIT.EQ.1)THEN !initialise valuator
        CALL SETUP(IW,ERROR,*9999)

        IF(NTYPE.EQ.1) THEN
          ECAREA(1)=XREF
          ECAREA(2)=ECAREA(1)+XFACTOR*XDISP
          ECAREA(3)=YREF
          ECAREA(4)=ECAREA(3)+YFACTOR*YDISP
        ELSE IF(NTYPE.EQ.2) THEN
          ECAREA(2)=XREF
          ECAREA(1)=ECAREA(2)-XFACTOR*XDISP
          ECAREA(3)=YREF
          ECAREA(4)=ECAREA(3)+YFACTOR*YDISP
        ELSE IF(NTYPE.EQ.3) THEN
          ECAREA(1)=                 (XREF-XMIN)/(XMAX-XMIN)*0.49*XDISP
          ECAREA(2)=ECAREA(1)+XFACTOR*XDISP
          ECAREA(4)=YDISP-0.49*XDISP+(YREF-YMIN)/(YMAX-YMIN)*0.49*XDISP
          ECAREA(3)=ECAREA(4)-YFACTOR*YDISP
        ELSE IF(NTYPE.EQ.4) THEN
          ECAREA(2)=                 (XREF-XMIN)/(XMAX-XMIN)*0.49*XDISP
          ECAREA(1)=ECAREA(2)-XFACTOR*XDISP
          ECAREA(4)=YDISP-0.49*XDISP+(YREF-YMIN)/(YMAX-YMIN)*0.49*XDISP
          ECAREA(3)=ECAREA(4)-YFACTOR*YDISP
        ELSE IF(NTYPE.EQ.5) THEN
          ECAREA(1)=0.5*XDISP+XREF*0.49*XDISP
          ECAREA(2)=ECAREA(1)+XFACTOR*XDISP
          ECAREA(4)=YDISP-0.375*XDISP+YREF*0.23/2.0*XDISP
          ECAREA(3)=ECAREA(4)-YFACTOR*YDISP
        ELSE IF(NTYPE.EQ.6) THEN
          ECAREA(2)=          XREF*0.49*XDISP
          ECAREA(1)=ECAREA(2)-XFACTOR*XDISP
          ECAREA(4)=YDISP-0.375*XDISP+YREF*0.23/2.0*XDISP
          ECAREA(3)=ECAREA(4)-YFACTOR*YDISP
        ELSE IF(NTYPE.EQ.9) THEN
          ECAREA(1)=XREF
          ECAREA(2)=ECAREA(1)+XFACTOR*XDISP
          ECAREA(4)=YREF
          ECAREA(3)=ECAREA(4)-YFACTOR*YDISP
        ELSE IF(NTYPE.EQ.10) THEN
          ECAREA(2)=XREF
          ECAREA(1)=ECAREA(2)-XFACTOR*XDISP
          ECAREA(4)=YREF
          ECAREA(3)=ECAREA(4)-YFACTOR*YDISP
        ENDIF
        IF(DOP) THEN
          WRITE(IOOP,'('' XREF='',E13.5,'' YREF='',E13.5)') XREF,YREF
          WRITE(IOOP,'('' DISP='',E13.5,'' ECAREA(1..4):'',4E13.5)')
     '      DISP,(ECAREA(I),I=1,4)
        ENDIF
        CALL GL_INIT_VALUATOR(IW,LD1,ECAREA,VALMIN,VALMAX)
      ENDIF

      IF(MODE.EQ.'REQUEST') THEN
C       get next event
        CALL EVENT(ID_WS,ID_DEVICE,INSTAT,CLASS,IDATA,RDATA,SDATA,
     '    ERROR,*9999)
C       if its a choice option in the correct window
        IF(INSTAT.EQ.1.AND.CLASS(1:6).EQ.'VALUATOR'.AND.ID_WS.EQ.IW)THEN
C         return choice
          VALUE=RDATA
        ELSE
        ENDIF
C       delete valuator
        CALL GL_DELETE_VALUATOR(IW)
      ENDIF

      CALL EXITS('VALUATOR')
      RETURN

 9999 CALL ERRORS('VALUATOR',ERROR)
      CALL EXITS('VALUATOR')
      RETURN 1
      END


      SUBROUTINE VISIB(IW,ISEG,ISEGNUM,CLASS,ERROR,*)

C**** change ISEGNUM to CLASS='VISIBLE' or 'INVISIBLE'
C**** change ISEG(ISEGNUM) to 2 if visible, 1 if not.

      INTEGER ISEG(*)
      CHARACTER CLASS*(*),ERROR*(*)

      CALL ENTERS('VISIB',*9999)
      IF(CLASS(1:7).EQ.'VISIBLE') THEN
        ISEG(ISEGNUM)=2
      ELSE
        ISEG(ISEGNUM)=1
      ENDIF
      CALL GL_INVIS(IW,ISEGNUM,ISEG(ISEGNUM))

      CALL EXITS('VISIB')
      RETURN
 9999 CALL ERRORS('VISIB',ERROR)
      CALL EXITS('VISIB')
      RETURN 1
      END


      SUBROUTINE WKST_VIEWPORT(IW,XMIN,XMAX,YMIN,YMAX)

C**** Change workstation window.

      CALL GL_WKST_VIEWPORT(IW,XMIN,XMAX,YMIN,YMAX)

      RETURN
      END


      SUBROUTINE WKST_WINDOW(IW,XNDC1,XNDC2,YNDC1,YNDC2,ERROR,*)

C**** Change workstation window.

      CHARACTER ERROR*(*)

      CALL ENTERS('WKST_WINDOW',*9999)
      CALL GL_WKST_WINDOW(IW,XNDC1,XNDC2,YNDC1,YNDC2)

      CALL EXITS('WKST_WINDOW')
      RETURN
 9999 CALL ERRORS('WKST_WINDOW',ERROR)
      CALL EXITS('WKST_WINDOW')
      RETURN 1
      END
