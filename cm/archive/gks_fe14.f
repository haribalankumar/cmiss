C**** CMISS Module FE14: Routines which call graphics buffer routines

C!!!! Note: Replace directory string "PRODUCT_1:[PRODUCT.CMISS.DOCUMENT]"
C           if necessary

C     Function INDEX_FILL_AREA      returns bundle index for fill-areas
C     Function INDEX_POLYLINE       returns bundle index for polylines
C     Function INDEX_POLYMARKER     returns bundle index for polymarkers
C     Function INDEX_TEXT           returns bundle index for text
C     Subroutine ACWK               activates workstation
C     Subroutine ARCHIVE_GRAPHICS   archives PHIGS structure or GKS segments
C     Subroutine BEZIER             creates Bezier curve
C     Subroutine BEZIER_POINTS      calculates Bezier curve points
C     Subroutine BUILD_XFORM_MATRIX3 redefine phigs transformation matrix
C     Subroutine CELL_ARRAY         draws bit-image cell array
C     Subroutine CHECK_GKS_OPEN     opens GKS if not already open
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
C     Subroutine DELETE_SEGMENT     deletes graphic segment
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
C     Subroutine INVIS              updates workstation invisibility filter
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
C     Subroutine PRINT_IMAGE_FILL_AREAS *** To FE_ARCHIVE ***
C     Subroutine QUIT_GRAPHICS      close any workstations and gks or phigs
C     Subroutine RECALL_GRAPHICS    recalls PHIGS structure or GKS segments
C     Subroutine SET_COLOUR_LUT     set colour lookup table
C     Subroutine SET_COLOUR_LUT_RANGE set colour lookup table for a given range
C     Subroutine SET_COLOUR_ONE     set RGB colour rep for one LUT index
C     Subroutine SET_COLOUR_REP     set RGB colour rep for all LUT indices
C     Subroutine SET_FILL_REP       resets fill area representation
C     Subroutine SET_FILL_AREA_REP  set polyline representation
C     Subroutine SET_POLYLINE_REP   set polyline representation
C     Subroutine SET_POLYMARKER_REP set polyline representation
C     Subroutine SET_TEXT_REP       set polyline representation
C     Subroutine SETUP              performs setup operations for workstations
C     Subroutine STROKE             calls GKS stroke
C     Subroutine SURFACE            draws shaded surface
C     Subroutine TEXT               draws text
C     Subroutine TRANSFORM_SEGMENT  transform graphics segment
C     Subroutine VALUATOR           calls GKS valuator
C     Subroutine VISIB              change segment visibility
C     Subroutine WKST_WINDOW        change workstation window


      INTEGER FUNCTION INDEX_FILL_AREA(icolour,PATTERN_FILL,STYLE_FILL,
     '  RGB_FILL)

C**** Returns fill-area bundle index:
C**** If icolour=0, for given pattern and style, as follows:
C****   Fill-area Index = 1
C**** Else if 1<=icolour<=233 for given icolour as follows:
C****   Fill-area Index = 17-249

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'
!     Parameter List
      INTEGER icolour
      CHARACTER PATTERN_FILL*(*),RGB_FILL*(*),STYLE_FILL*(*)
!     Local Variables
      INTEGER index
      CHARACTER ERROR*10,CUPPER*8,PATTERN*8,RGB*8,STYLE*8

      IF(icolour.EQ.0) THEN
        PATTERN=CUPPER(PATTERN_FILL)
        STYLE  =CUPPER(STYLE_FILL)
        RGB    =CUPPER(RGB_FILL)
        IF(.NOT.COLOUR_WS) RGB='BLACK'

        IF(     RGB(1:5).EQ.'BLACK' .OR.RGB(1:3).EQ.'000') THEN !black
        ELSE IF(RGB(1:3).EQ.'RED'   .OR.RGB(1:3).EQ.'100') THEN !red
        ELSE IF(RGB(1:5).EQ.'GREEN' .OR.RGB(1:3).EQ.'010') THEN !green
        ELSE IF(RGB(1:6).EQ.'YELLOW'.OR.RGB(1:3).EQ.'110') THEN !yellow
        ELSE IF(RGB(1:4).EQ.'BLUE'  .OR.RGB(1:3).EQ.'001') THEN !blue
        ELSE IF(RGB(1:4).EQ.'CYAN'  .OR.RGB(1:3).EQ.'011') THEN !cyan
        ELSE IF(RGB(1:5).EQ.'WHITE' .OR.RGB(1:3).EQ.'111') THEN !white
          index=0
        ENDIF
        IF(DOP) THEN
          WRITE(OP_STRING,
     '      '('' Index_Fill_Area: PATTERN='',A,'' STYLE='',A,'
     '      //''' RGB='',I3,'' INDEX='',I3)') PATTERN,STYLE,RGB,index
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

      ELSE IF(icolour.GE.1.AND.icolour.LE.(MAXCOLOURS-16)) THEN

        index=icolour+16

        IF(DOP) THEN
          WRITE(OP_STRING,'('' Index_Fill_Area: icolour='',I3,'//
     '      ''' INDEX='',I3)') icolour,index
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

      ENDIF

      INDEX_FILL_AREA=index

 9999 RETURN
      END


      INTEGER FUNCTION INDEX_POLYLINE(icolour,TYPE_LINE,WIDTH_LINE,
     ' RGB_LINE)

C**** Returns polyline bundle index:
C**** If icolour=0, for given line type, width and RGB value, as follows:
C****   Polyline Index = 1       black solid    width1
C****                    2         "   dotted      "
C****                    3         "   dashed      "
C****                    4         "   dot-dash    "
C****                    5-8         as above   width2
C****
C****                    9 -12   red    as above width1
C****                    13-16   green      "
C****                    17-20   blue       "
C****                    21-24   cyan       "
C****                    25-28   yellow     "
C****                    29-32   white      "
C****                    33-36   light blue "
C****                    37-40   grey       "
C****
C**** Else if 1<=icolour<=216 for given icolour as follows:
C****   Polyline Index = 41-256  216 colours solid width2

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'
!     Parameter List
      INTEGER icolour
      CHARACTER RGB_LINE*(*),TYPE_LINE*(*),WIDTH_LINE*(*)
!     Local Variables
      INTEGER index,INDEX_OFFSET_1,INDEX_OFFSET_2
      CHARACTER ERROR*10,CUPPER*8,TYPE*8,WIDTH*8,RGB*8

      IF(icolour.EQ.0) THEN
        TYPE =CUPPER(TYPE_LINE)
        WIDTH=CUPPER(WIDTH_LINE)
        RGB  =CUPPER(RGB_LINE)
        IF(.NOT.COLOUR_WS) RGB='BLACK'

        IF(     RGB(1:5).EQ.'BLACK' .OR.RGB(1:3).EQ.'000') THEN
          INDEX_OFFSET_1=0
        ELSE IF(RGB(1:3).EQ.'RED'   .OR.RGB(1:3).EQ.'100') THEN
          INDEX_OFFSET_1=8
        ELSE IF(RGB(1:5).EQ.'GREEN' .OR.RGB(1:3).EQ.'010') THEN
          INDEX_OFFSET_1=12
        ELSE IF(RGB(1:4).EQ.'BLUE'  .OR.RGB(1:3).EQ.'001') THEN
          INDEX_OFFSET_1=16
        ELSE IF(RGB(1:4).EQ.'CYAN'  .OR.RGB(1:3).EQ.'011') THEN
          INDEX_OFFSET_1=20
        ELSE IF(RGB(1:6).EQ.'YELLOW'.OR.RGB(1:3).EQ.'110') THEN
          INDEX_OFFSET_1=24
        ELSE IF(RGB(1:5).EQ.'WHITE' .OR.RGB(1:3).EQ.'111') THEN
          INDEX_OFFSET_1=28
        ELSE IF(RGB(1:6).EQ.'LTBLUE'.OR.RGB(1:3).EQ.'012') THEN
          INDEX_OFFSET_1=32
        ELSE IF(RGB(1:4).EQ.'GREY'  .OR.RGB(1:3).EQ.'210') THEN
          INDEX_OFFSET_1=36
        ELSE
          WRITE(OP_STRING,'(''>>Colour not defined'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

        IF(WIDTH(1:6).EQ.'WIDTH1') THEN
          INDEX_OFFSET_2=0
        ELSE IF(WIDTH(1:6).EQ.'WIDTH2') THEN
          INDEX_OFFSET_2=4
        ENDIF
        IF(TYPE(1:5).EQ.'SOLID') THEN
          index=1
        ELSE IF(TYPE(1:6).EQ.'DOTTED') THEN
          index=2
        ELSE IF(TYPE(1:6).EQ.'DASHED') THEN
          index=3
        ELSE IF(TYPE(1:8).EQ.'DOT-DASH') THEN
          index=4
        ENDIF

        index=index+INDEX_OFFSET_1+INDEX_OFFSET_2

        IF(DOP) THEN
          WRITE(OP_STRING,'('' Index_Polyline: Type='',A,'
     '      //''' Width='',A,'' RGB='',A,'' INDEX='',I3)')
     '      TYPE,WIDTH,RGB,index
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

      ELSE IF(icolour.GE.1.AND.icolour.LE.(MAXCOLOURS-40)) THEN

        index=icolour+40

        IF(DOP) THEN
          WRITE(OP_STRING,'('' Index_Polyline: icolour='',I3,'//
     '    '''INDEX='',I3)') icolour,index
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

      ENDIF

      INDEX_POLYLINE=index

 9999 RETURN
      END


      INTEGER FUNCTION INDEX_POLYMARKER(icolour,TYPE_MARKER,SIZE_MARKER,
     '  RGB_MARKER)

C**** Returns polymarker bundle index for given marker type, size & RGB value.
C**** If icolour=0  for given marker type, size and RGB value, as follows:
C****   Polymarker Index = 1       black plus      size1
C****                      2         "   asterisk    "
C****                      3         "   circle      "
C****                      4         "   point       "
C****                      5-8         as above    size2
C****
C****                      9 -12   red    as above  size1
C****                      13-16   green      "
C****                      17-20   blue       "
C****                      21-24   cyan       "
C****                      25-28   yellow     "
C****                      29-32   white      "
C****                      33-36   light blue "
C****                      37-40   grey       "
C****
C**** Else if 1<=icolour<=216 for given icolour as follows:
C****   Polymarker Index = 41-256  216 colours plus size2

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'
!     Parameter List
      INTEGER icolour
      CHARACTER RGB_MARKER*(*),SIZE_MARKER*(*),TYPE_MARKER*(*)
!     Local Variables
      INTEGER index,INDEX_OFFSET_1,INDEX_OFFSET_2
      REAL*8 WIDTH
      CHARACTER ERROR*10,CUPPER*8,RGB*8,SIZE*8,TYPE*8

      IF(icolour.EQ.0) THEN
        TYPE=CUPPER(TYPE_MARKER)
        SIZE=CUPPER(SIZE_MARKER)
        RGB =CUPPER(RGB_MARKER)
        IF(.NOT.COLOUR_WS) RGB='BLACK'

        IF(     RGB(1:5).EQ.'BLACK'.OR.RGB(1:3).EQ.'000') THEN
          INDEX_OFFSET_1=0
        ELSE IF(RGB(1:3).EQ.'RED'   .OR.RGB(1:3).EQ.'100') THEN
          INDEX_OFFSET_1=8
        ELSE IF(RGB(1:5).EQ.'GREEN' .OR.RGB(1:3).EQ.'010') THEN
          INDEX_OFFSET_1=12
        ELSE IF(RGB(1:4).EQ.'BLUE'  .OR.RGB(1:3).EQ.'001') THEN
          INDEX_OFFSET_1=16
        ELSE IF(RGB(1:4).EQ.'CYAN'  .OR.RGB(1:3).EQ.'011') THEN
          INDEX_OFFSET_1=20
        ELSE IF(RGB(1:6).EQ.'YELLOW'.OR.RGB(1:3).EQ.'110') THEN
          INDEX_OFFSET_1=24
        ELSE IF(RGB(1:5).EQ.'WHITE' .OR.RGB(1:3).EQ.'111') THEN
          INDEX_OFFSET_1=28
        ELSE IF(RGB(1:6).EQ.'LTBLUE'.OR.RGB(1:3).EQ.'012') THEN
          INDEX_OFFSET_1=32
        ELSE IF(RGB(1:4).EQ.'GREY'  .OR.RGB(1:3).EQ.'210') THEN
          INDEX_OFFSET_1=36
        ELSE
          WRITE(OP_STRING,'(''>>Colour not defined'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

        IF(SIZE(1:5).EQ.'SIZE1') THEN
          INDEX_OFFSET_2=0
        ELSE IF(SIZE(1:5).EQ.'SIZE2') THEN
          INDEX_OFFSET_2=4
        ENDIF
        IF(TYPE(1:4).EQ.'PLUS') THEN
          index=1
        ELSE IF(TYPE(1:8).EQ.'ASTERISK') THEN
          index=2
        ELSE IF(TYPE(1:6).EQ.'CIRCLE') THEN
          index=3
        ELSE IF(TYPE(1:5).EQ.'POINT') THEN
          index=4
        ENDIF

        index=index+INDEX_OFFSET_1+INDEX_OFFSET_2

        IF(DOP) THEN
          WRITE(OP_STRING,'('' Index_Polymarker: Type='',A,'
     '    //''' Width='',A,'' RGB='',A,'' INDEX='',I3)')
     '    TYPE,WIDTH,RGB,index
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

      ELSE IF(icolour.GE.1.AND.icolour.LE.(MAXCOLOURS-40)) THEN

        index=icolour+40

        IF(DOP) THEN
          WRITE(OP_STRING,'('' Index_Polymarker: icolour='',I3,'//
     '      ''' INDEX='',I3)')icolour,index
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

      ENDIF

      INDEX_POLYMARKER=index

 9999 RETURN
      END


      INTEGER FUNCTION INDEX_TEXT(icolour,WIDTH_TEXT,FONT_TEXT,RGB_TEXT)

C**** Returns text bundle index for given text size,font and RGB value.
C**** For given text size and RGB value, as follows:
C****         Text Index = 1       black font1    width1
C****                      2         "   font2      "
C****                      3         "   font3      "
C****                      4         "   font4      "
C****                      5-8         as above for width2
C****
C****                      9 -10   red    font1,2 width1
C****                      11-12    "     font1,2 width2
C****                      13-16   green  as above
C****                      17-20   blue       "
C****                      21-24   cyan       "
C****                      25-28   yellow     "
C****                      29-32   white      "
C****                      33-36   light blue "
C****                      37-40   grey       "
C****

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'
!     Parameter List
      INTEGER icolour
      CHARACTER FONT_TEXT*(*),RGB_TEXT*(*),WIDTH_TEXT*(*)
!     Local Variables
      INTEGER index,INDEX_OFFSET_1,INDEX_OFFSET_2
      CHARACTER ERROR*10,CUPPER*8,FONT*8,RGB*8,WIDTH*8

      IF(icolour.EQ.0) THEN
        FONT=CUPPER(FONT_TEXT)
        WIDTH=CUPPER(WIDTH_TEXT)
        RGB =CUPPER(RGB_TEXT)
        IF(.NOT.COLOUR_WS) RGB='BLACK'

        IF(     RGB(1:5).EQ.'BLACK'.OR.RGB(1:3).EQ.'000') THEN
          INDEX_OFFSET_1=0
        ELSE IF(RGB(1:3).EQ.'RED'   .OR.RGB(1:3).EQ.'100') THEN
          INDEX_OFFSET_1=8
        ELSE IF(RGB(1:5).EQ.'GREEN' .OR.RGB(1:3).EQ.'010') THEN
          INDEX_OFFSET_1=12
        ELSE IF(RGB(1:4).EQ.'BLUE'  .OR.RGB(1:3).EQ.'001') THEN
          INDEX_OFFSET_1=16
        ELSE IF(RGB(1:4).EQ.'CYAN'  .OR.RGB(1:3).EQ.'011') THEN
          INDEX_OFFSET_1=20
        ELSE IF(RGB(1:6).EQ.'YELLOW'.OR.RGB(1:3).EQ.'110') THEN
          INDEX_OFFSET_1=24
        ELSE IF(RGB(1:5).EQ.'WHITE' .OR.RGB(1:3).EQ.'111') THEN
          INDEX_OFFSET_1=28
        ELSE IF(RGB(1:6).EQ.'LTBLUE'.OR.RGB(1:3).EQ.'012') THEN
          INDEX_OFFSET_1=32
        ELSE IF(RGB(1:4).EQ.'GREY'  .OR.RGB(1:3).EQ.'210') THEN
          INDEX_OFFSET_1=36
        ELSE
          WRITE(OP_STRING,'(''>>Colour not defined'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

        IF(WIDTH(1:6).EQ.'WIDTH1') THEN
          INDEX_OFFSET_2=0
        ELSE IF(WIDTH(1:6).EQ.'WIDTH2') THEN
          IF(RGB(1:5).EQ.'BLACK'.OR.RGB(1:3).EQ.'000') THEN
            INDEX_OFFSET_2=4
          ELSE
            INDEX_OFFSET_2=2
          ENDIF
        ENDIF
        IF(FONT(1:5).EQ.'FONT1') THEN
          index=1
        ELSE IF(FONT(1:5).EQ.'FONT2') THEN
          index=2
        ELSE IF(FONT(1:5).EQ.'FONT3') THEN
          index=3
        ELSE IF(FONT(1:5).EQ.'FONT4') THEN
          index=4
        ENDIF

        index=index+INDEX_OFFSET_1+INDEX_OFFSET_2

        IF(DOP) THEN
          WRITE(OP_STRING,'('' Index_Text: Width='',A,'' RGB='',A,'//
     '      ''' INDEX='',I3)') WIDTH,RGB,index
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

      ELSE IF(icolour.GE.1.AND.icolour.LE.(MAXCOLOURS-40)) THEN

        index=icolour+40

        IF(DOP) THEN
          WRITE(OP_STRING,'('' Index_Text: icolour='',I3,'//
     '      ''' INDEX='',I3)')icolour,index
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

      ENDIF
      INDEX_TEXT=index

 9999 RETURN
      END


      SUBROUTINE ACWK(iw,ID,ERROR,*)

C#### Subroutine: ACWK
C###  Description:
C**** Activate workstation identified by iw.
C**** If ID=1 output is deferred.
C**** If iw=5 or 6 buffers are cleared.
C**** If POSTSCRIPT is .true. (as set by call to WS_LIST)
C**** then workstation 15 (GKS) or 16 (Phigs) is activated also.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:disp00.cmn'
      INCLUDE 'cmiss$reference:post00.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.FOR'
      INCLUDE 'SYS$LIBRARY:PHIGS$DEFS.FOR'
!     Parameter List
      INTEGER iw,ID
      CHARACTER ERROR*(*)
!     Local Variables
      CHARACTER CFROMI*2,CHAR2*2

      CALL ENTERS('ACWK',*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,'('' iw='',I3,'' IWKS(iw)='',I2,'//
     '    ''' IWKT(iw)='',I2)') iw,IWKS(iw),IWKT(iw)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(iw.EQ.8) THEN !help window
        IF(IWKS(iw).EQ.0) CALL SETUP(iw,ERROR,*9999)
      ELSE IF(iw.EQ.10.OR.iw.EQ.11) THEN !history & section windows
        IF(IWKS(iw).EQ.0) CALL SETUP(iw,ERROR,*9999)
      ELSE IF(iw.EQ.12.OR.iw.EQ.13) THEN !fibre or sheet angle windows
        IF(IWKS(iw).EQ.0) CALL SETUP(iw,ERROR,*9999)
      ELSE IF(iw.EQ.15.OR.iw.EQ.16) THEN !postscript windows
        IF(IWKS(iw).EQ.0) CALL SETUP(iw,ERROR,*9999)
      ELSE IF(iw.EQ.50.OR.iw.EQ.51) THEN !gen strain windows
        IF(IWKS(iw).EQ.0) CALL SETUP(iw,ERROR,*9999)
      ELSE IF(iw.GE.60.AND.iw.LE.69) THEN !spreadsheet windows
        IF(IWKS(iw).EQ.0) CALL SETUP(iw,ERROR,*9999)
      ENDIF

      IF(IWKS(iw).EQ.1) THEN !workstation iw is defined
        IF(IWKT(iw).EQ.1) THEN      !GKS window
          IF(DOP) THEN
            WRITE(OP_STRING,'('' Activate workstation iw='',i3)') iw
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CALL GKS_ACWK(iw,ERROR,*9999)
          IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            CALL GKS_SELNT(iw,ERROR,*9999)
            IF(ID.EQ.0) THEN      !perform screen write immediately
              CALL GKS_SDS(iw,GASAP,GALLOW,ERROR,*9999)
            ELSE IF(ID.EQ.1) THEN !defer screen write until update workstation
              CALL GKS_SDS(iw,GASTI,GSUPPD,ERROR,*9999)
            ENDIF
          ENDIF
        ELSE IF(IWKT(iw).EQ.2) THEN !Phigs window
          IF(ID.EQ.0) THEN      !use 'quick update' mode
            CALL PHIGS$SET_DISPLAY_UPDATE_STATE(iw,PHIGS$K_ASAP,
     '        PHIGS$K_UQUM)
          ELSE IF(ID.EQ.1) THEN !use 'no immediate visual effect' mode
            CALL PHIGS$SET_DISPLAY_UPDATE_STATE(iw,PHIGS$K_WAIT,
     '        PHIGS$K_NIVE)
          ENDIF
        ELSE IF(IWKT(iw).EQ.3) THEN !Frame grabber display screen
        ENDIF
        IWKS(iw)=2

      ELSE
        CHAR2=CFROMI(iw,'(I2)')
        ERROR=' >>Workstation '//CHAR2//' is not defined '
        GO TO 9999
      ENDIF

      IF(POSTSCRIPT) THEN
        IF(IWKT(iw).EQ.1) THEN      !iw is GKS window
          IF(DOP) THEN
            WRITE(OP_STRING,'('' Activate workstation iw=15'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CALL GKS_ACWK(15,ERROR,*9999)
        ELSE IF(IWKT(iw).EQ.2) THEN !iw is Phigs window
          CALL PHIGS$SET_DISPLAY_UPDATE_STATE(16,PHIGS$K_ASAP,
     '      PHIGS$K_UQUM)
        ENDIF
      ENDIF

      CALL EXITS('ACWK')
      RETURN
 9999 CALL ERRORS('ACWK',ERROR)
      CALL EXITS('ACWK')
      RETURN 1
      END


      SUBROUTINE ANIMATE(A,*)

C#### Subroutine: ANIMATE
C###  Description:

      INCLUDE 'cmiss$reference:cbdi02.cmn'

      DIMENSION A(*)

      CHARACTER ERROR*10
      WRITE(OP_STRING,*) '>>Cannot do under Phigs'
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
 9999 RETURN 1
      END


      SUBROUTINE ARCHIVE_GRAPHICS(FILE_NAME,ERROR,*)

C#### Subroutine: ARCHIVE_GRAPHICS
C###  Description:
C**** Archives PHIGS structure or GKS segments.

      IMPLICIT NONE
      INCLUDE 'SYS$LIBRARY:PHIGS$DEFS.FOR'
!     Parameter List
      CHARACTER ERROR*(*),FILE_NAME*(*)
!     Local Variables
      INTEGER IBEG,IEND

      CALL ENTERS('ARCHIVE_GRAPHICS',*9999)
      CALL PHIGS$SET_CONFLICT_RESOLUTION(PHIGS$K_UPDATE,PHIGS$K_UPDATE)
      CALL STRING_TRIM(FILE_NAME,IBEG,IEND)
      CALL PHIGS$OPEN_ARCHIVE(1,FILE_NAME(IBEG:IEND)//'.ARCHIVE')
      CALL PHIGS$ARCHIVE_ALL_STRUCT(1)
      CALL PHIGS$CLOSE_ARCHIVE(1)

      CALL EXITS('ARCHIVE_GRAPHICS')
      RETURN
 9999 CALL ERRORS('ARCHIVE_GRAPHICS',ERROR)
      CALL EXITS('ARCHIVE_GRAPHICS')
      RETURN 1
      END


      SUBROUTINE BEZIER(index,INDEX_PLIN,ISBEZE,ISEG,iw,NTL,
     '  XPTS,YPTS,CSEG,ERROR,*)

C#### Subroutine: BEZIER
C###  Description:
C**** Creates Bezier control points on workstation iw.
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
      INTEGER index,INDEX_PLIN,ISBEZE,ISEG(*),iw,NTL
      REAL*8 XPTS(*),YPTS(*)
      CHARACTER CSEG(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER NTLM
      PARAMETER(NTLM=20)
      INTEGER i,IFROMC,INDEX_OLD,INSTAT,IPICK,isegm,ISL2BE(NTLM),
     '  ISL3BE(NTLM),ISN2BE(NTLM),ISN3BE(NTLM),LD1,N2LI,nl
      REAL*8 D_XWC,D_YWC
      REAL*8 BEZLEN(NTLM),PL(3,21),PT(3,2),XBEZ(4,NTLM),XSLOPE,
     '  YBEZ(4,NTLM),YSLOPE
      DATA LD1/1/,BLANK/' '/

      CALL ENTERS('BEZIER',*9999)
      CALL ACWK(iw,0,ERROR,*9999)

      DO nl=1,NTL
        IF(DOP) THEN
          WRITE(OP_STRING,*) 'nl=',nl
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
C ***   Draw polymarkers at Bezier control points and tangent lines
        DO i=1,4
          XBEZ(I,nl)=XPTS(3*(nl-1)+i)
          YBEZ(I,nl)=YPTS(3*(nl-1)+i)
        ENDDO
        BEZLEN(nl)=DSQRT((XBEZ(4,nl)-XBEZ(1,nl))**2+
     '    (YBEZ(4,nl)-YBEZ(1,nl))**2)
        IF(DOP) THEN
          WRITE(OP_STRING,'('' XBEZ: '',4E12.3)') (XBEZ(i,nl),i=1,4)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' YBEZ: '',4E12.3)') (YBEZ(i,nl),i=1,4)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

        CALL OPEN_SEGMENT(ISN2BE(nl),ISEG,iw,'CPT1',
     '    index,INDEX_OLD,nl,1,CSEG,ERROR,*9999)
        PT(1,2)=XBEZ(2,nl)
        PT(2,2)=YBEZ(2,nl)
        CALL POLYMARKER(1,iw,1,PT(1,2),ERROR,*9999)
        CALL CLOSE_SEGMENT(ISN2BE(nl),iw,ERROR,*9999)
        CALL DETECT(iw,ISEG,ISN2BE(nl),'DETECTABLE',ERROR,*9999)

        CALL OPEN_SEGMENT(ISL2BE(nl),ISEG,iw,'TAN1',
     '    index,INDEX_OLD,nl,1,CSEG,ERROR,*9999)
        PT(1,1)=XBEZ(1,nl)
        PT(2,1)=YBEZ(1,nl)
        CALL POLYLINE(3,iw,2,PT(1,1),ERROR,*9999)
        CALL CLOSE_SEGMENT(ISL2BE(nl),iw,ERROR,*9999)

        CALL OPEN_SEGMENT(ISN3BE(nl),ISEG,iw,'CPT2',
     '    index,INDEX_OLD,nl,1,CSEG,ERROR,*9999)
        PT(1,1)=XBEZ(3,nl)
        PT(2,1)=YBEZ(3,nl)
        CALL POLYMARKER(1,iw,1,PT(1,1),ERROR,*9999)
        CALL CLOSE_SEGMENT(ISN3BE(nl),iw,ERROR,*9999)
        CALL DETECT(iw,ISEG,ISN3BE(nl),'DETECTABLE',ERROR,*9999)

        CALL OPEN_SEGMENT(ISL3BE(nl),ISEG,iw,'TAN2',
     '    index,INDEX_OLD,nl,1,CSEG,ERROR,*9999)
        PT(1,2)=XBEZ(4,nl)
        PT(2,2)=YBEZ(4,nl)
        CALL POLYLINE(3,iw,2,PT(1,1),ERROR,*9999)
        CALL CLOSE_SEGMENT(ISL3BE(nl),iw,ERROR,*9999)
      ENDDO

      INSTAT=1
      DO WHILE(INSTAT.EQ.1)
C ***   Pick control point segments
        WRITE(OP_STRING,
     '    '('' >>Pick control point on '',I1,'':'')') iw
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        CALL PICK(iw,LD1,'REQUEST',INSTAT,isegm,IPICK,ERROR,*9999)
        IF(INSTAT.EQ.1) THEN
          IF(CSEG(isegm)(1:4).EQ.'CPT1') THEN
            nl=IFROMC(CSEG(isegm)(53:57))
C ***       Locate new control point
            WRITE(OP_STRING,
     '        '('' >>Relocate 1st control point for line '','
     '        //'I4,'' on '',I1,'':'')') nl,iw
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            CALL LOCATOR(1,iw,INSTAT,'REQUEST',3,0.0D0,D_XWC,0.0D0,
     '        D_YWC,ERROR,*9999)
            XBEZ(2,nl)=D_XWC
            YBEZ(2,nl)=D_YWC

            CALL OPEN_SEGMENT(ISN2BE(nl),ISEG,iw,'CPT1',
     '        index,INDEX_OLD,nl,1,CSEG,ERROR,*9999)
            PT(1,2)=XBEZ(2,nl)
            PT(2,2)=YBEZ(2,nl)
            CALL POLYMARKER(1,iw,1,PT(1,2),ERROR,*9999)
            CALL CLOSE_SEGMENT(ISN2BE(nl),iw,ERROR,*9999)
            CALL DETECT(iw,ISEG,ISN2BE(nl),'DETECTABLE',ERROR,*9999)

            CALL OPEN_SEGMENT(ISL2BE(nl),ISEG,iw,'TAN1',
     '        index,INDEX_OLD,nl,1,CSEG,ERROR,*9999)
            PT(1,1)=XBEZ(1,nl)
            PT(2,1)=YBEZ(1,nl)
            CALL POLYLINE(3,iw,2,PT(1,1),ERROR,*9999)
            CALL CLOSE_SEGMENT(ISL2BE(nl),iw,ERROR,*9999)

C ***       Find adjacent control point if colinear
            IF(nl.ne.1) THEN
              N2LI=nl-1
              XSLOPE=(XBEZ(2,nl)-XBEZ(1,nl))/BEZLEN(nl)
              YSLOPE=(YBEZ(2,nl)-YBEZ(1,nl))/BEZLEN(nl)
              XBEZ(3,N2LI)=XBEZ(1,nl)-XSLOPE*BEZLEN(N2LI)
              YBEZ(3,N2LI)=YBEZ(1,nl)-YSLOPE*BEZLEN(N2LI)

              CALL OPEN_SEGMENT(ISN3BE(N2LI),ISEG,iw,'CPT2',
     '          index,INDEX_OLD,N2LI,1,CSEG,ERROR,*9999)
              PT(1,1)=XBEZ(3,N2LI)
              PT(2,1)=YBEZ(3,N2LI)
              CALL POLYMARKER(1,iw,1,PT(1,1),ERROR,*9999)
              CALL CLOSE_SEGMENT(ISN3BE(N2LI),iw,ERROR,*9999)
              CALL DETECT(iw,ISEG,ISN3BE(N2LI),'DETECTABLE',ERROR,*9999)

              CALL OPEN_SEGMENT(ISL3BE(N2LI),ISEG,iw,'TAN2',
     '          index,INDEX_OLD,N2LI,1,CSEG,ERROR,*9999)
              PT(1,2)=XBEZ(4,N2LI)
              PT(2,2)=YBEZ(4,N2LI)
              CALL POLYLINE(3,iw,2,PT(1,1),ERROR,*9999)
              CALL CLOSE_SEGMENT(ISL3BE(N2LI),iw,ERROR,*9999)
            ENDIF

          ELSE IF(CSEG(isegm)(1:4).EQ.'CPT2') THEN
            nl=IFROMC(CSEG(isegm)(53:57))
C ***       Locate new control point
            WRITE(OP_STRING,
     '        '('' >>Relocate 2nd control point for line '','
     '        //'I4,'' on '',I1,'':'')') nl,iw
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            CALL LOCATOR(1,iw,INSTAT,'REQUEST',3,0.0D0,D_XWC,0.0D0,
     '        D_YWC,ERROR,*9999)
            XBEZ(3,nl)=REAL(D_XWC)
            YBEZ(3,nl)=REAL(D_YWC)

            CALL OPEN_SEGMENT(ISN3BE(nl),ISEG,iw,'CPT2',
     '        index,INDEX_OLD,nl,1,CSEG,ERROR,*9999)
            PT(1,1)=XBEZ(3,nl)
            PT(2,1)=YBEZ(3,nl)
            CALL POLYMARKER(1,iw,1,PT(1,1),ERROR,*9999)
            CALL CLOSE_SEGMENT(ISN3BE(nl),iw,ERROR,*9999)
            CALL DETECT(iw,ISEG,ISN3BE(nl),'DETECTABLE',ERROR,*9999)

            CALL OPEN_SEGMENT(ISL3BE(nl),ISEG,iw,'TAN2',
     '        index,INDEX_OLD,nl,1,CSEG,ERROR,*9999)
            PT(1,2)=XBEZ(4,nl)
            PT(2,2)=YBEZ(4,nl)
            CALL POLYLINE(3,iw,2,PT(1,1),ERROR,*9999)
            CALL CLOSE_SEGMENT(ISL3BE(nl),iw,ERROR,*9999)

C ***       Find adjacent control point if colinear
            IF(nl.ne.NTL) THEN
              N2LI=nl+1
              XSLOPE=(XBEZ(4,nl)-XBEZ(3,nl))/BEZLEN(nl)
              YSLOPE=(YBEZ(4,nl)-YBEZ(3,nl))/BEZLEN(nl)
              XBEZ(2,N2LI)=XBEZ(4,nl)+XSLOPE*BEZLEN(N2LI)
              YBEZ(2,N2LI)=YBEZ(4,nl)+YSLOPE*BEZLEN(N2LI)

              CALL OPEN_SEGMENT(ISN2BE(N2LI),ISEG,iw,'CPT1',
     '          index,INDEX_OLD,N2LI,1,CSEG,ERROR,*9999)
              PT(1,2)=XBEZ(2,N2LI)
              PT(2,2)=YBEZ(2,N2LI)
              CALL POLYMARKER(1,iw,1,PT(1,2),ERROR,*9999)
              CALL CLOSE_SEGMENT(ISN2BE(N2LI),iw,ERROR,*9999)
              CALL DETECT(iw,ISEG,ISN2BE(N2LI),'DETECTABLE',ERROR,*9999)

              CALL OPEN_SEGMENT(ISL2BE(N2LI),ISEG,iw,'TAN1',
     '          index,INDEX_OLD,N2LI,1,CSEG,ERROR,*9999)
              PT(1,1)=XBEZ(1,N2LI)
              PT(2,1)=YBEZ(1,N2LI)
              CALL POLYLINE(3,iw,2,PT(1,1),ERROR,*9999)
              CALL CLOSE_SEGMENT(ISL2BE(N2LI),iw,ERROR,*9999)
            ENDIF
          ENDIF

C ***     Calculate line length and redraw line segments
          CALL OPEN_SEGMENT(ISBEZE,ISEG,iw,'polyline_bezier',index,
     '      INDEX_OLD,INDEX_PLIN,1,CSEG,ERROR,*9999)
          DO nl=1,NTL
            IF(DOP) THEN
              WRITE(OP_STRING,'(/'' nl= '',I3)') nl
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            CALL BEZIER_POINTS(PL,XBEZ(1,nl),YBEZ(1,nl),ERROR,*9999)
            CALL POLYLINE(1,iw,21,PL,ERROR,*9999)
          ENDDO
          CALL CLOSE_SEGMENT(ISBEZE,iw,ERROR,*9999)
        ENDIF
      ENDDO
      CALL DAWK(iw,0,ERROR,*9999)

      CALL ACWK(iw,1,ERROR,*9999)
      DO nl=NTL,1,-1
        IF(ISEG(ISL3BE(nl)).GT.0) THEN
          CALL DELETE_SEGMENT(ISL3BE(nl),ISEG,iw,ERROR,*9999)
          NTSG=NTSG-1 !reduce total segment count by 1
        ENDIF
        IF(ISEG(ISN3BE(nl)).GT.0) THEN
          CALL DELETE_SEGMENT(ISN3BE(nl),ISEG,iw,ERROR,*9999)
          NTSG=NTSG-1 !reduce total segment count by 1
        ENDIF
        IF(ISEG(ISL2BE(nl)).GT.0) THEN
          CALL DELETE_SEGMENT(ISL2BE(nl),ISEG,iw,ERROR,*9999)
          NTSG=NTSG-1 !reduce total segment count by 1
        ENDIF
        IF(ISEG(ISN2BE(nl)).GT.0) THEN
          CALL DELETE_SEGMENT(ISN2BE(nl),ISEG,iw,ERROR,*9999)
          NTSG=NTSG-1 !reduce total segment count by 1
        ENDIF
        DO i=1,4
          XPTS(3*(nl-1)+i)=XBEZ(i,nl)
          YPTS(3*(nl-1)+i)=YBEZ(i,nl)
        ENDDO
      ENDDO
      CALL DAWK(iw,1,ERROR,*9999)

      CALL EXITS('BEZIER')
      RETURN
 9999 CALL ERRORS('BEZIER',ERROR)
      CALL EXITS('BEZIER')
      RETURN 1
      END


      SUBROUTINE BEZIER_POINTS(PL,XBEZ,YBEZ,ERROR,*)

C#### Subroutine: BEZIER_POINTS
C###  Description:
C**** Calculates 21 points along Bezier curve given by XBEZ(i),YBEZ(i),i=1,4.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
!     Parameter List
      REAL*8 PL(3,*),XBEZ(*),YBEZ(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,k
      REAL*8 XI,XI2,XI3

      CALL ENTERS('BEZIER_POINTS',*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,'('' XBEZ: '',4E12.3)') (XBEZ(i),i=1,4)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' YBEZ: '',4E12.3)') (YBEZ(i),i=1,4)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      DO j=0,20
        XI=DBLE(j)/20.0D0
        XI2=XI*XI
        XI3=XI2*XI
        PL(1,j+1)=(1.0D0-3.0D0*XI+3.0D0*XI2-XI3) * XBEZ(1) +
     '            (3.0D0*XI-6.0D0*XI2+3.0D0*XI3) * XBEZ(2) +
     '            (3.0D0*XI2-3.0D0*XI3         ) * XBEZ(3) +
     '            (XI3                         ) * XBEZ(4)
        PL(2,J+1)=(1.0D0-3.0D0*XI+3.0D0*XI2-XI3) * YBEZ(1) +
     '            (3.0D0*XI-6.0D0*XI2+3.0D0*XI3) * YBEZ(2) +
     '            (3.0D0*XI2-3.0D0*XI3         ) * YBEZ(3) +
     '            (XI3                         ) * YBEZ(4)
        IF(DOP) THEN
          WRITE(OP_STRING,'('' J,XI= '',I3,F12.5)') j,XI
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' PL 1,2 = '',2F12.4)')
     '      (PL(k,j+1),k=1,2)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDDO

      CALL EXITS('BEZIER_POINTS')
      RETURN
 9999 CALL ERRORS('BEZIER_POINTS',ERROR)
      CALL EXITS('BEZIER_POINTS')
      RETURN 1
      END


      SUBROUTINE BUILD_XFORM_MATRIX3(FIXED_PT,SHIFT,
     '  ANGLE1,ANGLE2,ANGLE3,SCALE,ISTATUS,A_TRANS,ERROR,*)

C#### Subroutine: BUILD_XFORM_MATRIX3
C###  Description:
C**** Redefines transformation matrix in root structure (for phigs only)

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:view00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'SYS$LIBRARY:PHIGS$DEFS.FOR'
!     Parameter List
      INTEGER ISTATUS
      REAL ANGLE1,ANGLE2,ANGLE3,A_TRANS(4,4),FIXED_PT(*),SCALE,SHIFT(*)
      CHARACTER ERROR*(*)

      CALL ENTERS('BUILD_XFORM_MATRIX3',*9999)
      CALL PHIGS$BUILD_XFORM_MATRIX3(FIXED_PT,SHIFT,
     '  ANGLE1,ANGLE2,ANGLE3,SCALE,ISTATUS,A_TRANS)
      IF(istatus.ne.0) THEN
        WRITE(OP_STRING,*)' ---Error from BUILD_XFORM_MATRIX3=',
     '    ISTATUS
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ELSE
        CALL PHIGS$SET_EDIT_MODE(PHIGS$K_EDIT_REPLACE)
        CALL PHIGS$OPEN_STRUCT(ISVIEW)
        CALL PHIGS$SET_ELEM_POINTER(1)
        CALL PHIGS$SET_GLOBAL_XFORM3(A_TRANS)
        CALL PHIGS$CLOSE_STRUCT()
      ENDIF

      CALL EXITS('BUILD_XFORM_MATRIX3')
      RETURN
 9999 CALL ERRORS('BUILD_XFORM_MATRIX3',ERROR)
      CALL EXITS('BUILD_XFORM_MATRIX3')
      RETURN 1
      END


      SUBROUTINE CELL_ARRAY(iw,NDIMX,NDIMY,XMIN_CA,XMAX_CA,YMIN_CA,
     '  YMAX_CA,ERROR,*)

C#### Subroutine: CELL_ARRAY
C###  Description:
C**** Draws a cell array given by array I2P(512*512,iw) on w/s IW(1,2).
C**** If NDIMX,Y are not equal to 512, the program only draws every
C**** 512/NDIMX,Yth pixel.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'
      INCLUDE 'cmiss$reference:disp00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'gx$path:gx.inc'
      INCLUDE 'cmiss$reference:pics00.cmn'
      INCLUDE 'cmiss$reference:pics01.cmn'
!     Parameter List
      INTEGER iw,NDIMX,NDIMY
      REAL XMAX_CA,XMIN_CA,YMAX_CA,YMIN_CA
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ERR,i,j,NSTEPX,NSTEPY
      REAL BNDRECT(4)

      CALL ENTERS('CELL_ARRAY',*9999)
      NSTEPX=NINT(DBLE(IMGX)/DBLE(NDIMX))
      NSTEPY=NINT(DBLE(IMGY)/DBLE(NDIMY))
      IF(DOP) THEN
        WRITE(OP_STRING,'('' Stepping by '',2I4)') NSTEPX,NSTEPY
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      IF(iw.EQ.1.OR.iw.EQ.2) THEN
        IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
          DO j=1,NDIMY
            DO i=1,NDIMX
              I2P(i,j,iw+2)=INT(DBLE(I2P(i*NSTEPX,j*NSTEPY,iw))
     '          *(gxNSPECT)/256.0)+gxSPOFF
            ENDDO
          ENDDO
          BNDRECT(1)=XMIN_CA
          BNDRECT(2)=YMIN_CA
          BNDRECT(3)=XMAX_CA
          BNDRECT(4)=YMAX_CA
          CALL DRWIMG(BNDRECT,NDIMX,NDIMY,I2P(1,1,iw+2),ERR)
          IF(ERR.GT.0) THEN
            ERROR='>>Error from GX'
            GOTO 9999
          ENDIF
        ELSE IF(WINDOW_TYPE(1:5).EQ.'VWS') THEN
          DO j=1,NDIMY
            DO i=1,NDIMX
              I2P(i,j,iw+2)=MAXCOLOURS-INT(DBLE(I2P(i*NSTEPX,j*NSTEPY,
     '          iw))*(MAXCOLOURS-3)/256.0)
            ENDDO
          ENDDO
          CALL GKS_CA(XMIN_CA,YMIN_CA,XMAX_CA,YMAX_CA,NXM,NXM,1,1,
     '      NDIMX,NDIMY,I2P(1,1,iw+2),ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('CELL_ARRAY')
      RETURN
 9999 CALL ERRORS('CELL_ARRAY',ERROR)
      CALL EXITS('CELL_ARRAY')
      RETURN 1
      END


      SUBROUTINE CHECK_GKS_OPEN(iw,ERROR,*)

C#### Subroutine: CHECK_GKS_OPEN
C###  Description:
C**** Opens GKS if not already open.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'
      INCLUDE 'cmiss$reference:disp00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:gks000.cmn'
      INCLUDE 'cmiss$reference:gks001.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'cmiss$reference:phig00.cmn'
      INCLUDE 'cmiss$reference:post00.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.FOR'
!     Parameter List
      INTEGER iw
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,ICOLOUR_FLAG,IERST,IDUMMY,IFLAGS(13),NCOLOURS,WSTYPE

      CALL ENTERS('CHECK_GKS_OPEN',*9999)

      IF(.NOT.GKS) THEN
        WSTYPE=GWSDEF
        IF(DOP) THEN
          WRITE(OP_STRING,*)' Setup opening GKS'
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        CALL GKS_OPKS(IOER,ERROR,*9999)
        GKS=.TRUE.
        IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
          COLOUR_WS=.TRUE.
          MAXCOLOURS=240
          CALL GKS_QCF(WSTYPE,IERST,NCOLOURS,ICOLOUR_FLAG,NINDICES,
     '      ERROR,*9999)
        ELSE IF(WINDOW_TYPE(1:5).EQ.'VWS') THEN
          DO i=1,13
            IFLAGS(i)=GBUNDL
          ENDDO
          CALL GKS_SASF(IFLAGS,ERROR,*9999)
!         Find screen size
          CALL GKS_QDSP(WSTYPE,IERST,IDUMMY,XDISP,YDISP,
     '      IDUMMY,IDUMMY,ERROR,*9999)
          DISP=MIN(XDISP,YDISP)
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' xdisp, ydisp =',XDISP,YDISP
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          RATIO=YDISP/XDISP
!         Inquire whether colour or monochrome workstation
          IF(iw.EQ.15) THEN
            CALL GKS_QCF(PSTYPE,IERST,NCOLOURS,ICOLOUR_FLAG,NINDICES,
     '        ERROR,*9999)
          ELSE
            CALL GKS_QCF(WSTYPE,IERST,NCOLOURS,ICOLOUR_FLAG,NINDICES,
     '        ERROR,*9999)
          ENDIF
          IF(DOP) THEN
            WRITE(OP_STRING,
     '        '('' no colours='',I8,'' ICOLOUR_FLAG='',I2,'
     '        //''' no predefined indices='',I4)')
     '        NCOLOURS,ICOLOUR_FLAG,NINDICES
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(NINDICES.GT.2) THEN
            COLOUR_WS=.TRUE.
            MAXCOLOURS=249
          ELSE
            COLOUR_WS=.FALSE.
            MAXCOLOURS=2
          ENDIF
        ENDIF
        IF(.NOT.PHIGS) CALL SET_COLOUR_LUT(COLOUR_LUT,ERROR,*9999)
      ENDIF

      CALL EXITS('CHECK_GKS_OPEN')
      RETURN
 9999 CALL ERRORS('CHECK_GKS_OPEN',ERROR)
      CALL EXITS('CHECK_GKS_OPEN')
      RETURN 1
      END


      SUBROUTINE CHOICE(LABEL,IDEVICE,IPAGE,INSTAT,iw,MODE,INCH,NTCH,
     '  NOCH,NOCO,NTYPE,CO,OPTION,STRING,XREF,YREF,ERROR,*)

C#### Subroutine: CHOICE
C###  Description:
C**** Prompts the user with NTCH choices listed in the array of
C**** strings OPTION
C**** IPAGE : page number to use
C**** IDEVICE indicates input device: 1 for normal 4 for mouse buttons
C**** MODE can be 'REQUEST' or 'EVENT'
C**** The list of choices is headed with the string LABEL
C**** INCH=the no of options displayed in one list
C**** Up to 100 (MAXCH) choices are permitted.
C**** NOCH=the returned choice
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
C**** (for NTYPE=5,6,7,8 XREF & YREF are factors multiplying
C****  XDISP & YDISP so that XDISP,YDISP need not be precalculated)
C**** LNCHMX is maximum number of characters in choice list

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:disp00.cmn'
      INCLUDE 'cmiss$reference:echo00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:gks000.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
!     Parameter List
      INTEGER IDEVICE,INCH,INSTAT,IPAGE,iw,NOCH,NOCO,NTCH,NTYPE
      REAL XREF,YREF
      CHARACTER CO(*)*(*),ERROR*(*),LABEL*(*),MODE*(*),OPTION(*)*(*),
     '  STRING*(*)
!     Local Variables
      INTEGER i,IBEG,ICOUNT,IEND,IMODE,LNCH(100),LNCHMX,
     '  n1ch,NSTART,NUM_CHOICES,NTPAGE,NULL,TOT_CHOICES
      REAL CHLN,CHSZ
      CHARACTER CHOICES(100)*(60)
      LOGICAL BREAK,CONTINUE

      CALL ENTERS('CHOICE',*9999)

      CALL CHECK_GKS_OPEN(iw,ERROR,*9999)
      NTPAGE=INT(NTCH/(INCH+1))+1
      IF(DOP) THEN
        WRITE(OP_STRING,'('' Menu is '',I3,'' pages of'')') NTPAGE
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(I3,'' options per page totalling '',I3,'//
     '    ''' options.'')') INCH,NTCH
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' NTYPE='',I2)') NTYPE
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      CONTINUE=.TRUE.
      DO WHILE(CONTINUE)
        CONTINUE=.FALSE.
        IF(IPAGE.LT.1) IPAGE=1
        IF(IPAGE.GT.NTPAGE) IPAGE=NTPAGE

        NSTART=(IPAGE-1)*INCH+1
        IF(NTCH.LE.INCH) THEN
          DO n1ch=1,NTCH
            CHOICES(n1ch)=OPTION(n1ch)
          ENDDO
          NUM_CHOICES=NTCH
        ELSE
          ICOUNT=1
          BREAK=.FALSE.
          DO n1ch=1,INCH
            IF(BREAK.OR.n1ch+NSTART.GT.NTCH) THEN
              CHOICES(ICOUNT)=' '
              BREAK=.TRUE.
            ELSE
              CHOICES(ICOUNT)=OPTION(n1ch+NSTART-1)
            ENDIF
            ICOUNT=ICOUNT+1
          ENDDO
          IF(NTCH.GT.INCH+NSTART-1) THEN
            CHOICES(ICOUNT)='>>>'
          ELSE
            CHOICES(ICOUNT)=' '
          ENDIF
          ICOUNT=ICOUNT+1
          IF(NSTART.GT.1) THEN
            CHOICES(ICOUNT)='<<<'
          ELSE
            CHOICES(ICOUNT)=' '
          ENDIF
          ICOUNT=ICOUNT+1
          CHOICES(ICOUNT)=OPTION(NTCH)
          NUM_CHOICES=ICOUNT
          IF(DOP) THEN
            WRITE(OP_STRING,'('' Page # '',I2)') IPAGE
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF

        CHSZ=0.005 !character height
        CHLN=0.004 !character width

        IF(IDEVICE.EQ.4) THEN !mouse buttons
          TOT_CHOICES=3
        ELSE                  !choice menu
          TOT_CHOICES=100
        ENDIF

!       Find choice string lengths LNCH(n1ch) and max length LNCHMX
        LNCHMX=1
        DO n1ch=1,NUM_CHOICES
          CALL STRING_TRIM(CHOICES(n1ch),IBEG,IEND)
          LNCH(n1ch)=IEND-IBEG+1
          IF(LNCH(n1ch).GT.LNCHMX) LNCHMX=LNCH(n1ch)
        ENDDO
        IF(DOP) THEN
          WRITE(OP_STRING,'('' LNCHMX = '',I3)') LNCHMX
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

!       Define echo area
        IF(NTYPE.EQ.0) THEN
          ECAREA(2)=XDISP
          ECAREA(1)=ECAREA(2)-CHLN*LNCHMX
          ECAREA(3)=0.0
          ECAREA(4)=ECAREA(3)+CHSZ*REAL(NUM_CHOICES)
        ELSE IF(NTYPE.EQ.1) THEN
          ECAREA(1)=XREF
          ECAREA(2)=ECAREA(1)+CHLN*LNCHMX
          ECAREA(3)=YREF
          ECAREA(4)=ECAREA(3)+CHSZ*REAL(NUM_CHOICES)
        ELSE IF(NTYPE.EQ.2) THEN
          ECAREA(2)=XREF
          ECAREA(1)=ECAREA(2)-CHLN*LNCHMX
          ECAREA(3)=YREF
          ECAREA(4)=ECAREA(3)+CHSZ*REAL(NUM_CHOICES)
        ELSE IF(NTYPE.EQ.3) THEN
          ECAREA(1)=                   (XREF-XMIN)/(XMAX-XMIN)*
     '      0.49*XDISP
          ECAREA(2)=ECAREA(1)+CHLN*LNCHMX
          IF(NJT.EQ.2) THEN
            ECAREA(4)=YDISP-0.49*XDISP+(YREF-YMIN)/(YMAX-YMIN)*
     '        0.49*XDISP
          ELSE IF(NJT.EQ.3) THEN
            ECAREA(4)=YDISP-0.49*XDISP+(YREF-ZMIN)/(ZMAX-ZMIN)*
     '        0.49*XDISP
          ENDIF
          ECAREA(3)=ECAREA(4)-CHSZ*REAL(NUM_CHOICES)
        ELSE IF(NTYPE.EQ.4) THEN
          ECAREA(2)=                   (XREF-XMIN)/(XMAX-XMIN)*
     '      0.49*XDISP
          ECAREA(1)=ECAREA(2)-CHLN*LNCHMX
          IF(NJT.EQ.2) THEN
            ECAREA(4)=YDISP-0.49*XDISP+(YREF-YMIN)/(YMAX-YMIN)*
     '        0.49*XDISP
          ELSE IF(NJT.EQ.3) THEN
            ECAREA(4)=YDISP-0.49*XDISP+(YREF-ZMIN)/(ZMAX-ZMIN)*
     '        0.49*XDISP
          ENDIF
          ECAREA(3)=ECAREA(4)-CHSZ*REAL(NUM_CHOICES)
        ELSE IF(NTYPE.EQ.5) THEN
          ECAREA(1)=XREF*XDISP
          ECAREA(2)=ECAREA(1)+CHLN*LNCHMX
          ECAREA(3)=YREF*YDISP
          ECAREA(4)=ECAREA(3)+CHSZ*REAL(NUM_CHOICES)
        ELSE IF(NTYPE.EQ.6) THEN
          ECAREA(2)=XREF*XDISP
          ECAREA(1)=ECAREA(2)-CHLN*LNCHMX
          ECAREA(3)=YREF*YDISP
          ECAREA(4)=ECAREA(3)+CHSZ*REAL(NUM_CHOICES)
        ELSE IF(NTYPE.EQ.7) THEN
          ECAREA(1)=XREF*XDISP
          ECAREA(2)=ECAREA(1)+CHLN*LNCHMX
          ECAREA(4)=YREF*YDISP
          ECAREA(3)=ECAREA(4)-CHSZ*REAL(NUM_CHOICES)
        ELSE IF(NTYPE.EQ.8) THEN
          ECAREA(2)=XREF*XDISP
          ECAREA(1)=ECAREA(2)-CHLN*LNCHMX
          ECAREA(4)=YREF*YDISP
          ECAREA(3)=ECAREA(4)-CHSZ*REAL(NUM_CHOICES)
        ELSE IF(NTYPE.EQ.9) THEN
          ECAREA(1)=XREF
          ECAREA(2)=ECAREA(1)+CHLN*LNCHMX
          ECAREA(4)=YREF
          ECAREA(3)=ECAREA(4)-CHSZ*REAL(NUM_CHOICES)
        ELSE IF(NTYPE.EQ.10) THEN
          ECAREA(2)=XREF
          ECAREA(1)=ECAREA(2)-CHLN*LNCHMX
          ECAREA(4)=YREF
          ECAREA(3)=ECAREA(4)-CHSZ*REAL(NUM_CHOICES)
        ENDIF

        IF(ECAREA(1).LT.0.0)   ECAREA(1)=0.0
        IF(ECAREA(2).GT.XDISP) ECAREA(2)=XDISP
        IF(ECAREA(3).LT.0.0)   ECAREA(3)=0.0
        IF(ECAREA(4).GT.YDISP) ECAREA(4)=YDISP

        IF(DOP) THEN
          WRITE(OP_STRING,
     '      '('' XREF='',E13.5,'' YREF='',E13.5)') XREF,YREF
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' DISP='',E13.5,'' ECAREA(1..4):'','
     '      //'4E13.5)') DISP,(ECAREA(i),i=1,4)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

        CALL SETUP(iw,ERROR,*9999)
        CALL GKS_SCHM(iw,IDEVICE,GREQU,GECHO,ERROR,*9999)
        CALL GKS_SFACI(3,ERROR,*9999)
        CALL GKS_SFAIS(GHATCH,ERROR,*9999)
        CALL GKS_INCH(iw,IDEVICE,GOK,NUM_CHOICES,3,CHOICES,
     '    NUM_CHOICES,LNCH,TOT_CHOICES,ERROR,*9999)
        IF(MODE(1:5).EQ.'EVENT') THEN
          IMODE=GEVENT
        ELSE IF(MODE(1:7).EQ.'REQUEST') THEN
          IMODE=GREQU
        ENDIF
        IF(IDEVICE.EQ.4.AND.iw.LE.3) THEN
          CALL GKS_SCHM(iw,IDEVICE,IMODE,GNECHO,ERROR,*9999)
        ELSE
          CALL GKS_SCHM(iw,IDEVICE,IMODE,GECHO,ERROR,*9999)
        ENDIF

        IF(DOP) THEN
          WRITE(OP_STRING,*) 'Mode=',MODE
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

        IF(IMODE.EQ.GREQU) THEN
          CALL GKS_RQCH(iw,IDEVICE,INSTAT,NOCH,ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' Call to GKS_RQCH completed'
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' Returned choice # '',I3)') NOCH
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(NOCH.EQ.NUM_CHOICES-1.and.ipage.ne.1) THEN
            CONTINUE=.TRUE.
            IPAGE=IPAGE-1
          ELSE IF(NOCH.EQ.NUM_CHOICES-2.and.ipage.ne.NTPAGE) THEN
            CONTINUE=.TRUE.
            IPAGE=IPAGE+1
          ELSE IF(NOCH.EQ.NUM_CHOICES) THEN
            NOCH=NTCH
          ELSE
            NOCH=(IPAGE-1)*INCH+NOCH
          ENDIF
        ENDIF
      ENDDO
      IF(DOP) THEN
        WRITE(OP_STRING,'('' Original choice # '',I3)') NOCH
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(IMODE.EQ.GREQU.AND.NTYPE.EQ.0) THEN
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


      SUBROUTINE CIRCLE(TYPE,iw,RADIUS,Z,NLINES,ERROR,*)

C#### Subroutine: CIRCLE
C###  Description:
C**** Draws polyline (TYPE='POLYLINE') or fill-area (TYPE='FILL AREA')
C**** circle on workstation iw with radius RADIUS about the point Z(nj)
C**** with NLINES lines (max 100).

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:trans00.cmn'
!     Parameter List
      INTEGER iw,NLINES
      REAL*8 RADIUS,Z(*)
      CHARACTER ERROR*(*),TYPE*(*)
!     Local Variables
      INTEGER i,nj
      REAL*8 PTS(3,100),THETA

      CALL ENTERS('CIRCLE',*9999)

      IF(DOP) THEN
        WRITE(OP_STRING,'('' Radius='',E12.4,'
     '    //''' Position='',3(E12.4,X))') RADIUS,(Z(nj),nj=1,NJT)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      DO i=1,NLINES+1
        THETA=2.0D0*PI*DBLE(i-1)/DBLE(NLINES)
        IF(NJT.EQ.2) THEN
          PTS(1,i)=Z(1)+RADIUS*DCOS(THETA)
          PTS(2,i)=Z(2)+RADIUS*DSIN(THETA)
          PTS(3,i)=0.0D0
        ELSE IF(NJT.EQ.3) THEN
          IF(iw.EQ.1) THEN
            PTS(1,i)=Z(1)+RADIUS*DCOS(THETA)
            PTS(2,i)=0.0D0
            PTS(3,i)=Z(3)+RADIUS*DSIN(THETA)
          ELSE IF(iw.EQ.2) THEN
            PTS(1,i)=0.0D0
            PTS(2,i)=Z(2)+RADIUS*DCOS(THETA)
            PTS(3,i)=Z(3)+RADIUS*DSIN(THETA)
          ELSE IF(iw.EQ.3) THEN
            CALL ZZ(Z,Z,TRANS)
            PTS(1,i)=Z(1)+RADIUS*DCOS(THETA)
            PTS(2,i)=Z(2)+RADIUS*DSIN(THETA)
            PTS(3,i)=0.0D0
          ENDIF
        ENDIF
      ENDDO

      IF(TYPE(1:8).EQ.'POLYLINE') THEN
        CALL POLYLINE(1,iw,NLINES+1,PTS,ERROR,*9999)
      ELSE IF(TYPE(1:9).EQ.'FILL AREA') THEN
        CALL FILL_AREA(1,iw,NLINES+1,PTS,ERROR,*9999)
      ENDIF

      CALL EXITS('CIRCLE')
      RETURN
 9999 CALL ERRORS('CIRCLE',ERROR)
      CALL EXITS('CIRCLE')
      RETURN 1
      END


      SUBROUTINE CLOSE_PRINT_FILE(TYPE,ERROR,*)

C#### Subroutine: CLOSE_PRINT_FILE
C###  Description:
C**** Dectivates print workstation and creates another printfile with
C**** extension .POST which is easier to edit than the .GKS file.
C**** (POST converts linefeeds to returns in GKS or PHIGS postscript files)
C**** TYPE can be 'POSTSCRIPT' , print wkst is 15(GKS) or 16(PHIGS)
C****          or 'METAFILE'   , print wkst is 17

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:disp00.cmn'
      INCLUDE 'cmiss$reference:gks000.cmn'
      INCLUDE 'cmiss$reference:gks001.cmn'
      INCLUDE 'cmiss$reference:post00.cmn'
!     Parameter List
      CHARACTER ERROR*(*),TYPE*(*)
!     Local Variables
      INTEGER noiw

      CALL ENTERS('CLOSE_PRINT_FILE',*9999)

      IF(TYPE(1:8).EQ.'METAFILE') THEN
        CALL DAWK(17,1,ERROR,*9999)
        CALL CLOSE_WS(17,ERROR,*9999)

      ELSE IF(TYPE(1:10).EQ.'POSTSCRIPT') THEN
        IF(IWKS(15).GT.0) THEN      !(GKS) 15 is open and possibly active
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL DAWK(15,1,ERROR,*9999)
C MPN 25-Mar-94: not implemented yet
C            CALL CLPRNT()
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            CALL DAWK(15,1,ERROR,*9999)
            CALL CLOSE_WS(15,ERROR,*9999)
            CALL POST(15,POSTFILE,ERROR,*9999)
          ENDIF
        ENDIF
        IF(IWKS(16).GT.0) THEN      !(PHIGS) 16 is open and possibly active
          CALL DAWK(16,1,ERROR,*9999)
          CALL CLOSE_WS(16,ERROR,*9999)
          CALL POST(16,POSTFILE,ERROR,*9999)
        ENDIF
        POSTSCRIPT=.FALSE.
        EPS=.FALSE.
      ENDIF

      GKS_WS_OPEN=.FALSE.        !Check for at least one open workstation
      DO noiw=1,99
        IF(IWKS(noiw).GT.0) GKS_WS_OPEN=.TRUE.
      ENDDO

      CALL EXITS('CLOSE_PRINT_FILE')
      RETURN
 9999 CALL ERRORS('CLOSE_PRINT_FILE',ERROR)
      CALL EXITS('CLOSE_PRINT_FILE')
      RETURN 1
      END


      SUBROUTINE CLOSE_SEGMENT(ISEGNUM,iw,ERROR,*)

C#### Subroutine: CLOSE_SEGMENT
C###  Description:
C**** Closes graphics segment ISEGNUM.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:view00.cmn'
!     Parameter List
      INTEGER ISEGNUM,iw
      CHARACTER ERROR*(*)

      CALL ENTERS('CLOSE_SEGMENT',*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,
     '    '('' iw='',I3,'' IWKS(iw)='',I2,'' IWKT(iw)='',I2,'
     '    //''' IWKG(iw)='',I2)') iw,IWKS(iw),IWKT(iw),IWKG(iw)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' ISEGNUM='',I5)') ISEGNUM
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      IF(IWKT(iw).EQ.1) THEN !GKS
        CALL GKS_CLSG(ERROR,*9999)
        IF(DOP) THEN
          WRITE(OP_STRING,
     '      '('' Close segment isegnum='',I5,'' on iw='',I2)')
     '      ISEGNUM,iw
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
      ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS
        CALL PHIGS$REMOVE_NAMES_FROM_SET(1,ISEGNUM)
        CALL PHIGS$CLOSE_STRUCT()
        CALL PHIGS$OPEN_STRUCT(ISVIEW)
        CALL PHIGS$EXECUTE_STRUCT(ISEGNUM)
        CALL PHIGS$CLOSE_STRUCT()
      ENDIF

      CALL EXITS('CLOSE_SEGMENT')
      RETURN
 9999 CALL ERRORS('CLOSE_SEGMENT',ERROR)
      CALL EXITS('CLOSE_SEGMENT')
      RETURN 1
      END


      SUBROUTINE CLOSE_WS(iw,ERROR,*)

C#### Subroutine: CLOSE_WS
C###  Description:
C**** Closes graphics workstation iw.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:post00.cmn'
      INCLUDE 'cmiss$reference:view00.cmn'
      INCLUDE 'SYS$LIBRARY:PHIGS$DEFS.FOR'
!     Parameter List
      INTEGER iw
      CHARACTER ERROR*(*)
!     Local Variables

      CALL ENTERS('CLOSE_WS',*9999)
      IF(IWKT(iw).EQ.1) THEN      !GKS window
        CALL GKS_CLWK(iw,ERROR,*9999)
        CALL CLOSEF(iw,ERROR,*9999)
      ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS window
C MPN 11/9/92 - If postscript output then view has not been
C posted yet, so post now (to prevent unwanted page breaks).
        IF(DEWIND_POSTSCRIPT) THEN
          CALL PHIGS$POST_STRUCT(iw,ISVIEW,1.0)
          CALL PHIGS$UPDATE_WS(iw,phigs$k_perform_flag)
          DEWIND_POSTSCRIPT=.FALSE.
        ENDIF
        CALL PHIGS$CLOSE_WS(iw)
C cpb 9/8/93 no need to close the file ???
      ENDIF
      IF(DOP) THEN
        WRITE(OP_STRING,
     '    '('' Unit='',I2,'' has been closed.'')') iw
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      IWKS(iw)=0

      CALL EXITS('CLOSE_WS')
      RETURN
 9999 CALL ERRORS('CLOSE_WS',ERROR)
      CALL EXITS('CLOSE_WS')
      RETURN 1
      END


      SUBROUTINE CLWS(ERROR,*)

C#### Subroutine: CLWS
C###  Description:
C**** Deactivates and closes all workstations.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER iw,noiw

      CALL ENTERS('CLWS',*9999)
      DO noiw=1,IWKDEF(0)
        iw=IWKDEF(noiw)
        IF(IWKT(iw).EQ.1) THEN      !GKS
          IF(IWKS(iw).EQ.2) THEN
            CALL GKS_DAWK(iw,ERROR,*9999)
            IWKS(iw)=1
          ENDIF
          IF(IWKS(iw).EQ.1) THEN
            CALL CLOSE_WS(iw,ERROR,*9999)
          ENDIF
        ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS
          IF(IWKS(iw).GT.0) THEN
            CALL PHIGS$CLOSE_WS(iw)
            IWKS(iw)=0
          ENDIF
        ENDIF
      ENDDO
      IWKDEF(0)=0

      CALL EXITS('CLWS')
      RETURN
 9999 CALL ERRORS('CLWS',ERROR)
      CALL EXITS('CLWS')
      RETURN 1
      END


      SUBROUTINE COLOUR(COLOUR_INDEX,iw,CO,STRING,XREF,YREF,ERROR,*)

C#### Subroutine: COLOUR
C###  Description:
C**** Returns colour index

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.FOR'
!     Parameter List
      INTEGER COLOUR_INDEX,iw
      REAL XREF,YREF
      CHARACTER CO(*)*(*),STRING*(*),ERROR*(*)
!     Local Variables
      INTEGER ID_DEVICE,ID_WS,INPUT_CLASS,INPUT_STATUS,INSTAT,LD1,NOCH,
     '  NOCO
      REAL BLUE_INTENS,GREEN_INTENS,RED_INTENS,YFACTOR,VALUE
      CHARACTER OPTION(3)*20
      LOGICAL ABBREV,CONTINUE,UPDATE

      DATA LD1/1/,YFACTOR/0.136/

      CALL ENTERS('COLOUR',*9999)
      OPTION(1)='Set_Valuators'
      OPTION(2)='Set_Current'
      OPTION(3)='Return'
      CALL CHOICE('COLOUR_INDEX',1,1,INSTAT,72,'EVENT',3,3,NOCH,NOCO,2,
     '  CO,OPTION,STRING,XREF,YREF+1.2*YFACTOR*YDISP,ERROR,*9999)
      CALL VALUATOR('RED'  ,81,'EVENT',2,0.,1.,0.5,
     '  VALUE,XREF,YREF+4.0*YFACTOR*YDISP,ERROR,*9999)
      CALL VALUATOR('GREEN',82,'EVENT',2,0.,1.,0.5,
     '  VALUE,XREF,YREF+3.0*YFACTOR*YDISP,ERROR,*9999)
      CALL VALUATOR('BLUE' ,83,'EVENT',2,0.,1.,0.5,
     '  VALUE,XREF,YREF+2.0*YFACTOR*YDISP,ERROR,*9999)
      CONTINUE=.TRUE.
      DO WHILE (CONTINUE)
        CALL GKS_WAIT(0.0,ID_WS,Input_Class,ID_Device,ERROR,*9999)
        IF(DOP) THEN
          WRITE(OP_STRING,*) ' Input_Class=',Input_Class
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,*) ' ID_WS=',ID_WS
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        IF(Input_Class.eq.GKS$K_Input_Class_Choice) THEN
          CALL GKS_GET_CHOICE(Input_Status,NOCH,ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' Input_Choice=',NOCH
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(NOCH.GT.0) THEN
            IF(ABBREV(OPTION(NOCH),'SET_CURRENT',8)) THEN
              CONTINUE=.FALSE.
              UPDATE=.TRUE.
            ELSE IF(ABBREV(OPTION(NOCH),'RETURN',6)) THEN
              CONTINUE=.FALSE.
            ENDIF
          ENDIF
        ELSE IF(Input_Class.eq.GKS$K_Input_Class_Valuator) THEN
          CALL GKS_GET_VALUATOR(VALUE,ERROR,*9999)
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
        CALL GKS_SCR(iw,COLOUR_INDEX,RED_INTENS,GREEN_INTENS,
     '    BLUE_INTENS,ERROR,*9999)
      ENDIF
      CALL GKS_SCHM(72,LD1,GREQU,GECHO,ERROR,*9999)
      CALL GKS_SVLM(81,LD1,GREQU,GECHO,ERROR,*9999)
      CALL GKS_SVLM(82,LD1,GREQU,GECHO,ERROR,*9999)
      CALL GKS_SVLM(83,LD1,GREQU,GECHO,ERROR,*9999)

      CALL EXITS('COLOUR')
      RETURN

 9999 CALL ERRORS('COLOUR',ERROR)
      CALL EXITS('COLOUR')
      RETURN 1
      END


      SUBROUTINE COWK(iw,ERROR,*)

C#### Subroutine: COWK
C###  Description:
C**** Closes and reopens workstation identified by iw (in order to remove
C**** workstation viewport from screen).

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:gks001.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.FOR'
      INCLUDE 'SYS$LIBRARY:PHIGS$DEFS.FOR'
!     Parameter List
      INTEGER iw
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ipat,LD1,LINE_COLOUR_INDEX,LINE_INDEX,
     '  LINE_TYPE1,LINE_TYPE2,LINE_TYPE3,LINE_TYPE4,
     '  MARKER_COLOUR_INDEX,MARKER_INDEX,
     '  MARKER_TYPE1,MARKER_TYPE2,MARKER_TYPE3,MARKER_TYPE4,
     '  MARKER_TYPE5,
     '  TEXT_COLOUR_INDEX,TEXT_FONT1,TEXT_FONT2,TEXT_FONT3,TEXT_INDEX,
     '  TEXT_PRECISION1,TEXT_PRECISION2,TEXT_PRECISION3
      REAL LINE_WIDTH1,LINE_WIDTH2,LINE_WIDTH3,LINE_WIDTH4,MARKER_SIZE,
     '  TEXT_EXPFAC,TEXT_SPACING
      DATA LD1/1/

      CALL ENTERS('COWK',*9999)
      LINE_COLOUR_INDEX=1
      LINE_INDEX=1
      LINE_TYPE1=GLSOLI
      LINE_TYPE2=GLDASH
      LINE_TYPE3=GLDOT
      LINE_TYPE4=GLSOLI
      LINE_WIDTH1=1.0
      LINE_WIDTH2=1.0
      LINE_WIDTH3=1.0
      LINE_WIDTH4=4.0
      MARKER_COLOUR_INDEX=1
      MARKER_INDEX=1
      MARKER_SIZE=1.0
      MARKER_TYPE1=GPLUS
      MARKER_TYPE2=GAST
      MARKER_TYPE3=GOMARK
      MARKER_TYPE4=GXMARK
      MARKER_TYPE5=GPOINT
      TEXT_COLOUR_INDEX=1
      TEXT_EXPFAC=0.5
      TEXT_FONT1=1
      TEXT_FONT2=-105
      TEXT_FONT3=1
      TEXT_INDEX=1
      TEXT_PRECISION1=GSTRP
      TEXT_PRECISION2=GCHARP
      TEXT_PRECISION3=GSTRKP
      TEXT_SPACING=0.0
      IF(IWKS(iw).EQ.2) THEN
        CALL GKS_DAWK(iw,ERROR,*9999)
        IWKS(iw)=1
      ENDIF
      IF(IWKS(iw).EQ.1) THEN
        CALL CLOSE_WS(iw,ERROR,*9999)
      ENDIF
      IF(iw.EQ.1) THEN      !close choice menu and valuators
        CALL GKS_SCHM(91,LD1,GREQU,GECHO,ERROR,*9999)
        CALL GKS_SVLM(83,LD1,GREQU,GECHO,ERROR,*9999)
        CALL GKS_SVLM(84,LD1,GREQU,GECHO,ERROR,*9999)
      ELSE IF(iw.EQ.2) THEN !close choice menu and valuators
        CALL GKS_SCHM(92,LD1,GREQU,GECHO,ERROR,*9999)
        CALL GKS_SVLM(85,LD1,GREQU,GECHO,ERROR,*9999)
        CALL GKS_SVLM(86,LD1,GREQU,GECHO,ERROR,*9999)
      ELSE IF(iw.EQ.3) THEN !close choice menu and valuators
        CALL GKS_SCHM(93,LD1,GREQU,GECHO,ERROR,*9999)
        CALL GKS_SVLM(87,LD1,GREQU,GECHO,ERROR,*9999)
        CALL GKS_SVLM(88,LD1,GREQU,GECHO,ERROR,*9999)
        CALL GKS_SVLM(89,LD1,GREQU,GECHO,ERROR,*9999)
      ELSE IF(iw.EQ.4) THEN !close choice menu and valuators
        CALL GKS_SCHM(94,LD1,GREQU,GECHO,ERROR,*9999)
        CALL GKS_SVLM(81,LD1,GREQU,GECHO,ERROR,*9999)
        CALL GKS_SVLM(82,LD1,GREQU,GECHO,ERROR,*9999)
      ENDIF
      IF(IWKS(iw).EQ.0) THEN
        CALL GKS_OPWK(iw,GCONID,GWSDEF,ERROR,*9999)
        IWKS(iw)=1
        IF(iw.EQ.4) THEN
          CALL GKS_SWKWN(4,0.0,1.0,0.0,1.0,ERROR,*9999)
          CALL GKS_SWKVP(4,0.2*DISP,1.0*DISP,0.0,0.8*DISP,ERROR,*9999)
          CALL GKS_SPLR(4,1,LINE_TYPE1,LINE_WIDTH1,
     '      LINE_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SPLR(4,2,LINE_TYPE2,LINE_WIDTH2,
     '      LINE_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SPLR(4,3,LINE_TYPE3,LINE_WIDTH3,
     '      LINE_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SPMR(4,1,MARKER_TYPE1,MARKER_SIZE,
     '      MARKER_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SPMR(4,2,MARKER_TYPE2,MARKER_SIZE,
     '      MARKER_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SPMR(4,3,MARKER_TYPE3,MARKER_SIZE,
     '      MARKER_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_STXR(4,1,TEXT_FONT1,TEXT_PRECISION1,
     '      TEXT_EXPFAC,TEXT_SPACING,TEXT_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SWN(4,-1.0,1.0,-1.0,1.0,ERROR,*9999)
        ELSE IF(iw.EQ.8) THEN
          CALL GKS_SWKWN(8,0.0,1.0,0.3688,1.0,ERROR,*9999)
          CALL GKS_SWKVP(8,0.2*DISP,1.0*DISP,0.295*DISP,0.8*DISP,
     '      ERROR,*9999)
          CALL GKS_SPLR(8,1,LINE_TYPE1,LINE_WIDTH1,
     '      LINE_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SPLR(8,2,LINE_TYPE2,LINE_WIDTH2,
     '      LINE_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SPLR(8,3,LINE_TYPE3,LINE_WIDTH3,
     '      LINE_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SPMR(8,1,MARKER_TYPE1,MARKER_SIZE,
     '      MARKER_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SPMR(8,2,MARKER_TYPE2,MARKER_SIZE,
     '      MARKER_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SPMR(8,3,MARKER_TYPE3,MARKER_SIZE,
     '      MARKER_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_STXR(8,1,TEXT_FONT1,TEXT_PRECISION1,
     '      TEXT_EXPFAC,TEXT_SPACING,TEXT_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SWN(8,0.0,1.0,0.0,1.0,ERROR,*9999)
        ELSE IF(iw.EQ.10) THEN
          CALL GKS_SWKWN(10,0.0,1.0,0.5-0.23/0.49/2.0,0.5+0.23/0.49/2.0,
     '      ERROR,*9999)
          CALL GKS_SWKVP(10,0.5*XDISP,0.99*XDISP,YDISP-0.23*XDISP,YDISP,
     '      ERROR,*9999)
          CALL GKS_SPLR(10,1,LINE_TYPE1,LINE_WIDTH1,
     '      LINE_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SPLR(10,2,LINE_TYPE2,LINE_WIDTH2,
     '      LINE_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SPLR(10,3,LINE_TYPE3,LINE_WIDTH3,
     '      LINE_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SPMR(10,1,MARKER_TYPE1,MARKER_SIZE,
     '      MARKER_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SPMR(10,2,MARKER_TYPE2,MARKER_SIZE,
     '      MARKER_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SPMR(10,3,MARKER_TYPE3,MARKER_SIZE,
     '      MARKER_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_STXR(10,1,TEXT_FONT1,TEXT_PRECISION1,
     '      TEXT_EXPFAC,TEXT_SPACING,TEXT_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SWN(10,-0.14,1.1,-2.3,2.3,ERROR,*9999)
        ELSE IF(iw.EQ.11) THEN
          CALL GKS_SWKWN(11,0.0,1.0,0.5-0.23/0.49/2.0,0.5+0.23/0.49/2.0,
     '      ERROR,*9999)
          CALL GKS_SWKVP(11,0.5*XDISP,0.99*XDISP,YDISP-0.49*XDISP,
     '      YDISP-0.26*XDISP,ERROR,*9999)
          CALL GKS_SPLR(11,1,LINE_TYPE1,LINE_WIDTH1,
     '      LINE_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SPLR(11,2,LINE_TYPE2,LINE_WIDTH2,
     '      LINE_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SPLR(11,3,LINE_TYPE3,LINE_WIDTH3,
     '      LINE_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SPMR(11,1,MARKER_TYPE1,MARKER_SIZE,
     '      MARKER_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SPMR(11,2,MARKER_TYPE2,MARKER_SIZE,
     '      MARKER_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SPMR(11,3,MARKER_TYPE3,MARKER_SIZE,
     '      MARKER_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_STXR(11,1,TEXT_FONT1,TEXT_PRECISION1,
     '      TEXT_EXPFAC,TEXT_SPACING,TEXT_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SWN(11,-0.14,1.1,-2.3,2.3,ERROR,*9999)
        ELSE IF(iw.EQ.61) THEN
          CALL GKS_STXR(61,1,TEXT_FONT1,TEXT_PRECISION1,
     '      TEXT_EXPFAC,TEXT_SPACING,TEXT_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SWN(61,-0.1,1.1,-0.1,1.1,ERROR,*9999)
          CALL GKS_SWKWN(61,0.0,1.0,0.0,1.0,ERROR,*9999)
          CALL GKS_SWKVP(61,0.03*XDISP,0.52*XDISP,
     '       YDISP-0.52*XDISP,YDISP-0.03*XDISP,ERROR,*9999)
        ELSE IF(iw.EQ.62) THEN
          CALL GKS_STXR(62,1,TEXT_FONT1,TEXT_PRECISION1,
     '      TEXT_EXPFAC,TEXT_SPACING,TEXT_COLOUR_INDEX,ERROR,*9999)
C ***     Set up fill area patterns to be 16 shades of grey
          IF(NINDICES.GE.16) THEN
            CALL GKS_SFAR(62,1,GSOLID,14,1,ERROR,*9999)
            CALL GKS_SFAR(62,16,GSOLID,14,0,ERROR,*9999)
            DO ipat=2,15
              CALL GKS_SFAR(62,ipat,GPATTR,12+(17-ipat),1,ERROR,*9999)
            ENDDO
          ENDIF
          CALL GKS_SWN(62,-0.1,1.1,-0.1,1.1,ERROR,*9999)
          CALL GKS_SWKWN(62,0.0,1.0,0.0,1.0,ERROR,*9999)
          CALL GKS_SWKVP(62,0.06*XDISP,0.55*XDISP,
     '       YDISP-0.55*XDISP,YDISP-0.06*XDISP,ERROR,*9999)
        ELSE IF(iw.EQ.63) THEN
          CALL GKS_STXR(63,1,TEXT_FONT1,TEXT_PRECISION1,
     '      TEXT_EXPFAC,TEXT_SPACING,TEXT_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SWN(63,-0.1,1.1,-0.1,1.1,ERROR,*9999)
          CALL GKS_SWKWN(63,0.0,1.0,0.0,1.0,ERROR,*9999)
          CALL GKS_SWKVP(63,0.09*XDISP,0.58*XDISP,
     '       YDISP-0.58*XDISP,YDISP-0.09*XDISP,ERROR,*9999)
        ELSE IF(iw.EQ.64) THEN
          CALL GKS_SPMR(64,1,MARKER_TYPE1,MARKER_SIZE,
     '      MARKER_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SPMR(64,2,MARKER_TYPE2,MARKER_SIZE,
     '      MARKER_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SPMR(64,3,MARKER_TYPE3,MARKER_SIZE,
     '      MARKER_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SPMR(64,4,MARKER_TYPE4,MARKER_SIZE,
     '      MARKER_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SPMR(64,5,MARKER_TYPE5,MARKER_SIZE,
     '      MARKER_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_STXR(64,1,TEXT_FONT1,TEXT_PRECISION1,
     '      TEXT_EXPFAC,TEXT_SPACING,TEXT_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SWN(64,-5.5,5.5,-5.5,5.5,ERROR,*9999)
          CALL GKS_SWKWN(64,0.0,1.0,0.0,1.0,ERROR,*9999)
          CALL GKS_SWKVP(64,0.12*XDISP,0.61*XDISP,
     '       YDISP-0.61*XDISP,YDISP-0.12*XDISP,ERROR,*9999)
        ELSE IF(iw.EQ.65) THEN
          CALL GKS_STXR(65,1,TEXT_FONT1,TEXT_PRECISION1,
     '      TEXT_EXPFAC,TEXT_SPACING,TEXT_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SWN(65,-0.1,1.1,-0.1,1.1,ERROR,*9999)
          CALL GKS_SWKWN(65,0.0,1.0,0.0,1.0,ERROR,*9999)
          CALL GKS_SWKVP(65,0.15*XDISP,0.64*XDISP,
     '       YDISP-0.64*XDISP,YDISP-0.15*XDISP,ERROR,*9999)
        ELSE IF(iw.EQ.66) THEN
          CALL GKS_STXR(66,1,TEXT_FONT1,TEXT_PRECISION1,
     '      TEXT_EXPFAC,TEXT_SPACING,TEXT_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SWN(66,-0.1,1.1,-0.1,1.1,ERROR,*9999)
          CALL GKS_SWKWN(66,0.0,1.0,0.0,1.0,ERROR,*9999)
          CALL GKS_SWKVP(66,0.18*XDISP,0.67*XDISP,
     '       YDISP-0.67*XDISP,YDISP-0.18*XDISP,ERROR,*9999)
        ELSE IF(iw.EQ.67) THEN
        ELSE IF(iw.EQ.68) THEN
        ELSE IF(iw.EQ.69) THEN
          CALL GKS_SPMR(69,1,MARKER_TYPE1,MARKER_SIZE,
     '      MARKER_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SPMR(69,2,MARKER_TYPE2,MARKER_SIZE,
     '      MARKER_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SPMR(69,3,MARKER_TYPE3,MARKER_SIZE,
     '      MARKER_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SPMR(69,4,MARKER_TYPE4,MARKER_SIZE,
     '      MARKER_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SPMR(69,5,MARKER_TYPE5,MARKER_SIZE,
     '      MARKER_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_STXR(69,1,TEXT_FONT1,TEXT_PRECISION1,
     '      TEXT_EXPFAC,TEXT_SPACING,TEXT_COLOUR_INDEX,ERROR,*9999)
          CALL GKS_SWN(69,-0.1,1.1,-0.1,1.1,ERROR,*9999)
          CALL GKS_SWKWN(69,0.0,1.0,0.0,1.0,ERROR,*9999)
          CALL GKS_SWKVP(69,0.12*XDISP,0.61*XDISP,
     '       YDISP-0.61*XDISP,YDISP-0.12*XDISP,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('COWK')
      RETURN

 9999 CALL ERRORS('COWK',ERROR)
      CALL EXITS('COWK')
      RETURN 1
      END


      SUBROUTINE CREATE_SEGMENT(ISEGNUM,iw,ERROR,*)

C#### Subroutine: CREATE_SEGMENT
C###  Description:
C**** create graphics segment ISEGNUM.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:view00.cmn'
      INCLUDE 'SYS$LIBRARY:PHIGS$DEFS.FOR'
!     Parameter List
      INTEGER ISEGNUM,iw
      CHARACTER ERROR*(*)

      CALL ENTERS('CREATE_SEGMENT',*9999)
      IF(IWKT(iw).EQ.1) THEN      !GKS
        CALL GKS_CRSG(ISEGNUM,ERROR,*9999)
      ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS
        CALL PHIGS$SET_EDIT_MODE(PHIGS$K_EDIT_INSERT)
        CALL PHIGS$OPEN_STRUCT(ISEGNUM)
        CALL PHIGS$ADD_NAMES_TO_SET(1,ISEGNUM)
      ENDIF

      CALL EXITS('CREATE_SEGMENT')
      RETURN
 9999 CALL ERRORS('CREATE_SEGMENT',ERROR)
      CALL EXITS('CREATE_SEGMENT')
      RETURN 1
      END


      SUBROUTINE DAWK(iw,ID,ERROR,*)

C#### Subroutine: DAWK
C###  Description:
C**** Deactivate workstation identified by iw.
C**** If ID=1 deferred output is released.
C**** If POSTSCRIPT is .true. (as set by call to WS_LIST)
C**** then workstation 15 (GKS) or 16 (Phigs) is deactivated also,
C**** and POSTSCRIPT is set to .false.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:disp00.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:post00.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.FOR'
      INCLUDE 'SYS$LIBRARY:PHIGS$DEFS.FOR'
!     Parameter List
      INTEGER ID,iw
      CHARACTER ERROR*(*)

      CALL ENTERS('DAWK',*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,
     '    '('' iw='',I3,'' IWKS(iw)='',I2,'//
     '    ''' IWKT(iw)='',I2)') iw,IWKS(iw),IWKT(iw)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(IWKT(iw).EQ.1) THEN      !GKS window
        IF(IWKS(iw).EQ.2) THEN
          IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            IF(ID.EQ.0) THEN
            ELSE IF(ID.EQ.1) THEN !update workstation
              CALL GKS_UWK(iw,GPERFO,ERROR,*9999)
              CALL GKS_SDS(iw,GASAP,GALLOW,ERROR,*9999)
            ENDIF
          ENDIF
          CALL GKS_DAWK(iw,ERROR,*9999)
          IWKS(iw)=1
        ENDIF

      ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS window
        IF(IWKS(iw).EQ.2) THEN
          IF(ID.EQ.0) THEN      !use 'quick update' mode
          ELSE IF(ID.EQ.1) THEN !use 'no immediate visual effect' mode
            IF(iw.ne.16) THEN ! Not a postscript workstation
              CALL PHIGS$UPDATE_WS(iw,phigs$k_perform_flag)
            ENDIF
          ENDIF
          IWKS(iw)=1
        ENDIF

      ELSE IF(IWKT(iw).EQ.3) THEN !Frame grabber display screen
        IF(IWKS(iw).EQ.2) THEN
          IWKS(iw)=1
        ENDIF
      ENDIF

C CPB 8/7/93 Is it necessary to update the PHIGS workstation here as it gets
C done in CLOSE_WS - could be why a blank page is printed????
      IF(POSTSCRIPT.OR.EPS) THEN
        IF(IWKT(iw).EQ.1) THEN      !iw is GKS window
          CALL GKS_DAWK(15,ERROR,*9999)
        ELSE IF(IWKT(iw).EQ.2) THEN !iw is Phigs window
C          CALL PHIGS$UPDATE_WS(16,phigs$k_perform_flag)
        ENDIF
      ENDIF

      CALL EXITS('DAWK')
      RETURN
 9999 CALL ERRORS('DAWK',ERROR)
      CALL EXITS('DAWK')
      RETURN 1
      END


      SUBROUTINE DBOX(iw,XDIM,YDIM,XWC,YWC,ZWC,ERROR,*)

C#### Subroutine: DBOX
C###  Description:
C**** Draws box of dimensions 2*XDIM,2*YDIM around the point XWC,YWC,ZWC

      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER iw
      REAL*8 XDIM,XWC,YDIM,YWC,ZWC
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 Z(3,5)

      CALL ENTERS('DBOX',*9999)
      IF(iw.EQ.1) THEN
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
          CALL POLYLINE(1,iw,5,Z,ERROR,*9999)
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
          CALL POLYLINE(1,iw,5,Z,ERROR,*9999)
        ENDIF
      ELSE IF(iw.EQ.2) THEN !plot Y against Z
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
        CALL POLYLINE(1,iw,5,Z,ERROR,*9999)
      ELSE IF(iw.EQ.3) THEN !plot X against Y
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
        CALL POLYLINE(1,iw,5,Z,ERROR,*9999)
      ELSE IF(iw.EQ.4) THEN
      ENDIF

      CALL EXITS('DBOX')
      RETURN
 9999 CALL ERRORS('DBOX',ERROR)
      CALL EXITS('DBOX')
      RETURN 1
      END


      SUBROUTINE DELETE_SEGMENT(ISEGNUM,ISEG,iw,ERROR,*)

C#### Subroutine: DELETE_SEGMENT
C###  Description:
C**** Deletes graphics segment ISEGNUM.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:view00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISEGNUM,iw
      CHARACTER ERROR*(*)

      CALL ENTERS('DELETE_SEGMENT',*9999)
      IF(IWKT(iw).EQ.1) THEN      !GKS
        CALL GKS_DSG(ISEGNUM,ERROR,*9999)
      ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS
        CALL PHIGS$DEL_STRUCT(ISEGNUM)
c       CALL PHIGS$POST_STRUCT(iw,ISVIEW,1.0) !why is this needed
      ENDIF
      ISEG(ISEGNUM)=0
      ISEGNUM=0

      CALL EXITS('DELETE_SEGMENT')
      RETURN
 9999 CALL ERRORS('DELETE_SEGMENT',ERROR)
      CALL EXITS('DELETE_SEGMENT')
      RETURN 1
      END


      SUBROUTINE DETECT(iw,ISEG,ISEGNUM,CLASS,ERROR,*)

C#### Subroutine: DETECT
C###  Description:
C**** Change ISEGNUM to CLASS='DETECTABLE' or 'UNDETECTABLE'
C**** Change ISEG(ISEGNUM) to 3 if detectable, 2 if not.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:disp00.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
!     Parameter List
      INTEGER ISEG(*),ISEGNUM,iw
      CHARACTER CLASS*(*),ERROR*(*)
!     Local Variables
      INTEGER ERR

      CALL ENTERS('DETECT',*9999)
      IF(iw.ne.3) THEN !gks
        IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
          CALL SETWIN(iw,ERR)
          IF(err.ne.0) THEN
            ERROR='>>Error from GX'
            GOTO 9999
          ENDIF
        ENDIF
        IF(CLASS(1:10).EQ.'DETECTABLE') THEN
          CALL GKS_SDTEC(ISEGNUM,GDETEC,ERROR,*9999)
          ISEG(ISEGNUM)=3
        ELSE
          CALL GKS_SDTEC(ISEGNUM,GUNDET,ERROR,*9999)
          ISEG(ISEGNUM)=2
        ENDIF
      ELSE IF(iw.EQ.3) THEN
      ENDIF

      CALL EXITS('DETECT')
      RETURN
 9999 CALL ERRORS('DETECT',ERROR)
      CALL EXITS('DETECT')
      RETURN 1
      END


      SUBROUTINE DETRAN(NOCO,NTCO,NTCOQU,CO,COQU,STRING,ERROR,*)

C#### Subroutine: DETRAN
C###  Description:
C**** Defines transformations and view matrices for PHIGS.
C**** A_TRANS is the 4*4 transformation matrix containing rotations
C**** shifts and scaling of an object in world coords
C**** A_ORIENT is a 4*4 matrix indicating the orientation of the
C**** world coords in the viewing reference coords (z-axis into the
C**** screen, y-axis vertical)
C**** A_MAP is a 4*4 matrix specifying parallel or perspective views
C**** together with clipping planes

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:back00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'cmiss$reference:ntsg00.cmn'
      INCLUDE 'cmiss$reference:phig00.cmn'
      INCLUDE 'cmiss$reference:view00.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.FOR'
      INCLUDE 'SYS$LIBRARY:PHIGS$DEFS.FOR'
!     Parameter List
      INTEGER NOCO,NTCO,NTCOQU(*)
      CHARACTER CO(*)*(*),COQU(NXCO,*)*(*),ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER IBEG,IBEG1,IEND,IEND1,IFROMC,ISTATUS,iw,MODE_PROJ,N3CO
      CHARACTER FILE*100,STATUS*3
      LOGICAL ABBREV,CBBREV,EVENT,FILIO

      CALL ENTERS('DETRAN',*9999)
 1    IF(CO(NOCO+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        WRITE(OP_STRING,'(X,A)') STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME>['//FILE00(IBEG1:IEND1)//']<;doc>'
     '    //' <on WS_ID>[3]'
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
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

        IF(CBBREV(CO,'ON',1,NOCO+1,NTCO,N3CO)) THEN
          IF(ABBREV(CO(N3CO+1),'POSTSCRIPT',1).OR.
     1       ABBREV(CO(N3CO+1),'ENCAPSULATED',1)) THEN
            IF(NJT.LE.2) THEN
              iw=15
            ELSE
              iw=16
            ENDIF
          ELSE
            iw=IFROMC(CO(N3CO+1))
          ENDIF
        ELSE
          iw=3
        ENDIF

        IF(FILIO) THEN
          IVDU=IOIP
          IFILE=1
          CALL CHECKF(2,NOCO,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.IPTRAN',STATUS,
     '      'DIRECT','FORMATTED',132,ERROR,*9999)
          CALL IPTRAN(ERROR,*9999)
          CALL CLOSEF(IFILE,ERROR,*9999)
          IF(IOTYPE.EQ.1.OR.IOTYPE.EQ.2) THEN
            CALL PHIGS$BUILD_XFORM_MATRIX3(FIXED_PT,SHIFT,
     '        ANGLE(1),ANGLE(2),ANGLE(3),SCALE,ISTATUS,A_TRANS)
            IF(istatus.ne.0) THEN
              WRITE(OP_STRING,*)
     '          ' >>Error from BUILD_XFORM_MATRIX3=',ISTATUS
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ELSE
              CALL PHIGS$SET_EDIT_MODE(PHIGS$K_EDIT_REPLACE)
              CALL PHIGS$OPEN_STRUCT(ISVIEW)
              CALL PHIGS$SET_ELEM_POINTER(1)
              CALL PHIGS$SET_GLOBAL_XFORM3(A_TRANS)
              CALL PHIGS$CLOSE_STRUCT()
            ENDIF
            MODE_PROJ=PHIGS$K_PARALLEL
            CALL PHIGS$EVAL_VIEW_ORIEN_MATRIX3(VIEW_REF_PT_NEW,
     '        VIEW_PLANE_NEW,VIEW_UP_NEW,ISTATUS,A_ORIENT)
            IF(istatus.ne.0) THEN
              WRITE(OP_STRING,*)
     '          ' >>Error from eval_view_orient_matrix=',ISTATUS
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ELSE
              CALL PHIGS$EVAL_VIEW_MAP_MATRIX3(WINDOW_NEW,VIEWPORT,
     '          MODE_PROJ,PROJ_REF_PT_NEW,VIEW_PLANE_DIST_NEW,
     '          BACK_PLANE_DIST_NEW,FRONT_PLANE_DIST_NEW,ISTATUS,A_MAP)
              IF(istatus.ne.0) THEN
                WRITE(OP_STRING,*)
     '            ' >>Error from eval_view_map_matrix=',ISTATUS
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ELSE
                CALL PHIGS$SET_VIEW_REP3(iw,1,A_ORIENT,A_MAP,NPC_CLIP,
     '            PHIGS$K_NOCLIP,PHIGS$K_NOCLIP,PHIGS$K_NOCLIP)
              ENDIF
            ENDIF
          ENDIF
        ENDIF

      ENDIF

      CALL EXITS('DETRAN')
      RETURN
 9999 CALL ERRORS('DETRAN',ERROR)
      CALL EXITS('DETRAN')
      RETURN 1
      END


      SUBROUTINE DISPLAY(CO,ERROR,*)

C#### Subroutine: DISPLAY
C###  Description:
C**** Displays string in a window.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:echo00.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.FOR'
!     Parameter List
      CHARACTER CO(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER IBEG,IEND,INSTAT,NCHAR
      CHARACTER MESSAGE*80
      DATA BLANK/' '/

      CALL ENTERS('DISPLAY',*9999)
      MESSAGE=CO(2)
      ECAREA(1)=0.5*XDISP
      ECAREA(2)=ECAREA(1)+0.5*XDISP
      ECAREA(3)=0.3*YDISP
      ECAREA(4)=ECAREA(3)+0.2*YDISP
      CALL STRING_TRIM(MESSAGE,IBEG,IEND)
      CALL GKS_INST(1,1,1,MESSAGE(IBEG:IEND),1,40,1,0,BLANK,ERROR,*9999)
      CALL GKS_RQST(1,1,INSTAT,NCHAR,MESSAGE,ERROR,*9999)

      CALL EXITS('DISPLAY')
      RETURN
 9999 CALL ERRORS('DISPLAY',ERROR)
      CALL EXITS('DISPLAY')
      RETURN 1
      END


      SUBROUTINE DISPLAY_FILE(iw,NOCO,NTFILE,CO,FILE_EXT,FILE_NAME,
     '  STRING,ERROR,*)

C#### Subroutine: DISPLAY_FILE
C###  Description:
C**** Displays files in current directory and returns chosen one.
C**** If none chosen FILE_NAME is 'Exit'.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:file00.cmn'
!     Parameter List
      INTEGER iw,NOCO,NTFILE
      CHARACTER CO(*)*(*),ERROR*(*),FILE_EXT*(*),FILE_NAME*(*),
     '  STRING*(*)
!     Local Variables
      INTEGER IBEG1,IBEG2,IEND1,IEND2,n1file,NOCH,nofile,NOFILE_DEFAULT,
     '  NOFILE_START,NTCH
      CHARACTER CHOOSE*20,OPTION(22)*20
      LOGICAL MORE_FILES

      CALL ENTERS('DISPLAY_FILE',*9999)
      NOFILE_START=0
 10   CALL FIND_FILE(NOFILE_START,NTFILE,FILE_EXT,OPTION,ERROR,*9999)
      IF(NTFILE.EQ.20) THEN
        MORE_FILES=.TRUE.
      ELSE
        MORE_FILES=.FALSE.
      ENDIF
      NOFILE_DEFAULT=0
      DO nofile=1,NTFILE !to display default in square brackets
        CALL STRING_TRIM(OPTION(nofile),IBEG1,IEND1)
        CALL STRING_TRIM(FILE00,IBEG2,IEND2)
        IF(OPTION(nofile)(IBEG1:IEND1).EQ.FILE00(IBEG2:IEND2)) THEN
          DO n1file=nofile,2,-1
            OPTION(n1file)=OPTION(n1file-1)
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
      CALL PRECHOICE1(2,iw,NOCH,NOCO,NTCH,CO,'REQUEST',OPTION,STRING,
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

C#### Subroutine: DOCUM
C###  Description:
C**** Displays documentation
C**** FILE identifies file to be read (file ID=8 and workstation ID=8)
C**** EXTEN is name of file extension (eg .DOC)
C**** STRING identifies string in documentation file at beginning & end
C****   of text to be displayed
C**** Documentation has root plus number of branches (max 20)
C**** Root and branches can each have arbitrary number of windows
C**** NT_WIND is number of windows for root or branch
C**** NT_LINE is number of lines in root or branch (this is divided by 23
C****   to give the number of windows)

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cmis00.cmn'
      INCLUDE 'cmiss$reference:comm00.cmn'
      INCLUDE 'cmiss$reference:disp00.cmn'
      INCLUDE 'cmiss$reference:docd00.cmn'
      INCLUDE 'cmiss$reference:docu00.cmn'
      INCLUDE 'cmiss$reference:gks000.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
!     Parameter List
      CHARACTER ERROR*(*),FILE*(*),EXTEN*(*),STRING*(*)
!     Local Variables
      INTEGER i,IBEG1,IBEG2,IBEG4,IEND1,IEND2,IEND4,
     '  KEY,LINE0,LINE1,LINE2,LINE3,LINE_MAX,NOWIND,NT_LINE,NT_WIND
      REAL*8 P(3)
      CHARACTER COPYFILE*12,CUPPER*80,STRG*80,
     '  TERM1*1,TERM2*1,TERM3*1,TERM4*1,TERM5*1,TOPIC*20
      LOGICAL CONTINUE,EXAMPLES,
     '  FIND_COMMAND,FIND_COMMON_BLOCK,FIND_COMMON_VARIABLE,
     '  FIND_INCLUDE_FILE,FIND_SUBROUTINE,FIND_VARIABLE,
     '  FOUND,
     '  KEY3,KEY4,KEY5,KEY32,
     '  OPTION1,OPTION2,OPTION3,OPTION4,OPTION5,
     '  READ_COM_FILE,ROOT

      CALL ENTERS('DOCUM',*9999)
      DOCDIR=.TRUE.
      DO_EXAMPLE=.FALSE.
      SELECT_EXAMPLE=.FALSE.

      KEY3 =.FALSE.
      KEY4 =.FALSE.
      KEY5 =.FALSE.
      KEY32=.FALSE.

      TOPIC=FILE
      IF(TOPIC(1:8).EQ.'EXAMPLES') THEN
        EXAMPLES=.TRUE.
        LINE_MAX=1000
      ELSE IF(TOPIC(1:11).EQ.'SUBROUTINES') THEN
        EXAMPLES=.FALSE.
        LINE_MAX=10000
      ELSE IF(TOPIC(1:7).EQ.'MODULES') THEN
        EXAMPLES=.FALSE.
        LINE_MAX=2000
      ELSE IF(TOPIC(1:13).EQ.'COMMON_BLOCKS') THEN
        EXAMPLES=.FALSE.
        LINE_MAX=1500
      ELSE
        EXAMPLES=.FALSE.
        LINE_MAX=1000
      ENDIF

      CALL STRING_TRIM(TOPIC,IBEG1,IEND1)
      CALL STRING_TRIM(EXTEN,IBEG2,IEND2)
      CALL OPENF(8,'DISK','PRODUCT_1:[PRODUCT.CMISS.DOCUMENT]'
     '  //TOPIC(IBEG1:IEND1)//'.'
     '  //EXTEN(IBEG2:IEND2),'OLD','DIRECT','FORMATTED',132,
     '  ERROR,*9999)
      WRITE(OP_STRING,'('' [Prev Screen] scrolls up one page'')')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' [Next Screen] scrolls down one page'')')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      IF(EXAMPLES) THEN
        WRITE(OP_STRING,'(''          [Do] exits from help and '
     '    //'reads current example file'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''      [Select] exits from help and '
     '    //'opens current example file'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''      [Insert] inserts current example '
     '    //'into local directory'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        KEY3 =.TRUE.
        KEY4 =.TRUE.
        KEY5 =.TRUE.
      ENDIF
      WRITE(OP_STRING,'(''      [Delete] returns to top page'')')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''      [Remove] exits from help'')')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ROOT=.TRUE.

      FIND_COMMAND     =.FALSE.
      FIND_COMMON_BLOCK=.FALSE.
      FIND_COMMON_VARIABLE=.FALSE.
      FIND_INCLUDE_FILE=.FALSE.
      FIND_SUBROUTINE  =.FALSE.
      FIND_VARIABLE    =.FALSE.
      IF(TOPIC(1:8).EQ.'COMMANDS'.
     '  AND.STRING(1:7).EQ.'Command') THEN
        FIND_COMMAND     =.TRUE.      !to display a particular command
      ELSE IF(TOPIC(1:13).EQ.'COMMON_BLOCKS') THEN
        IF(STRING(1:7).EQ.'Include') THEN
          FIND_INCLUDE_FILE=.TRUE.    !to display a particular include file
        ELSE IF(STRING(1:8).EQ.'COMMON /') THEN
          FIND_COMMON_BLOCK=.TRUE.    !to display a particular commonblock
        ELSE IF(STRING(1:8).EQ.'VARIABLE') THEN
          FIND_COMMON_VARIABLE=.TRUE. !to display a commonblock variable
        ENDIF
      ELSE IF(TOPIC(1:11).EQ.'SUBROUTINES'.
     '  AND.STRING(1:10).EQ.'Subroutine') THEN
        FIND_SUBROUTINE  =.TRUE.      !to display a particular subroutine
      ELSE IF(TOPIC(1:9).EQ.'VARIABLES'.
     '  AND.STRING(1:8).EQ.'Variable') THEN
        FIND_VARIABLE    =.TRUE.      !to display a particular variable
      ENDIF
      STRG=CUPPER(STRING)
      CDATA(2)=STRG

      OPTION1=.FALSE.
      OPTION2=.FALSE.
      OPTION3=.FALSE.
      OPTION4=.FALSE.
      OPTION5=.FALSE.
      TERM1=' '
      TERM2=' '
      TERM3=' '
      TERM4=' '
      TERM5=' '
      IF(.NOT.USE_SOCKET) THEN
        CALL ACWK(8,0,ERROR,*9999)
      ENDIF

 10   CALL STRING_TRIM(CDATA(2),IBEG2,IEND2)
      IF(DOP) THEN
        WRITE(OP_STRING,'('' String: '',A)') CDATA(2)(IBEG2:IEND2)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(EXAMPLES) THEN !enable numeric key press
        KEY32=.TRUE.
      ENDIF
      IF(ROOT) THEN
        IF(FIND_COMMAND) THEN              !Find command name
          CALL FIND(8,'Command    '//CDATA(2)(12:IEND2),1,
     '      LINE_MAX,LINE1,FOUND,ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' LINE1=',LINE1,' FOUND=',FOUND
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          endif
          CALL FIND(8,'Command',LINE1+1,LINE1+50,LINE2,FOUND,
     '      ERROR,*9999)
          LINE1=LINE1-1 !to ensure that command name is displayed
        ELSE IF(FIND_COMMON_BLOCK) THEN    !Find common-block name
          CALL FIND(8,'COMMON '//CDATA(2)(8:IEND2),1,
     '      LINE_MAX,LINE1,FOUND,ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' LINE1=',LINE1,' FOUND=',FOUND
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CALL FIND(8,'Include',LINE1,LINE1-20,LINE1,FOUND,ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' LINE1=',LINE1,' FOUND=',FOUND
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          LINE1=LINE1-1
          LINE2=LINE1+22
        ELSE IF(FIND_COMMON_VARIABLE) THEN !Find common-block variable
          CALL FIND(8,CDATA(2)(10:IEND2),1,
     '      LINE_MAX,LINE1,FOUND,ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' LINE1=',LINE1,' FOUND=',FOUND
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CALL FIND(8,'Include',LINE1,LINE1-20,LINE1,FOUND,ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' LINE1=',LINE1,' FOUND=',FOUND
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          LINE1=LINE1-1
          LINE2=LINE1+22
        ELSE IF(FIND_INCLUDE_FILE) THEN    !Find include-file name
          CALL FIND(8,'Include '//CDATA(2)(9:IEND2)//'.cmn',1,
     '      LINE_MAX,LINE1,FOUND,ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' LINE1=',LINE1,' FOUND=',FOUND
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          LINE1=LINE1-1
          LINE2=LINE1+22
        ELSE IF(FIND_SUBROUTINE) THEN      !Find subroutine name
          CALL FIND(8,'Subroutine    '//CDATA(2)(12:IEND2),1,
     '      LINE_MAX,LINE1,FOUND,ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' LINE1=',LINE1,' FOUND=',FOUND
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CALL FIND(8,'Subroutine',LINE1+1,LINE1+50,LINE2,FOUND,
     '      ERROR,*9999)
          LINE1=LINE1-1 !to ensure that subroutine name is displayed
        ELSE IF(FIND_VARIABLE) THEN        !Find variable name
          CALL FIND(8,CDATA(2)(10:IEND2),1,
     '      LINE_MAX,LINE1,FOUND,ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' LINE1=',LINE1,' FOUND=',FOUND
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          LINE1=LINE1-1 !to ensure that subroutine name is displayed
          LINE2=LINE1+22
        ELSE
          CALL FIND(8,CDATA(2)(IBEG2:IEND2),1,LINE_MAX,LINE1,FOUND,
     '      ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' LINE1=',LINE1,' FOUND=',FOUND
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CALL FIND(8,CDATA(2)(IBEG2:IEND2),LINE1+1,LINE1+LINE_MAX,
     '      LINE2,FOUND,ERROR,*9999)
          IF(TOPIC(1:7).EQ.'MODULES') LINE1=LINE1-1 !to ensure name is displayed
        ENDIF
        IF(DOP) THEN
          WRITE(OP_STRING,*) ' LINE2=',LINE2,' FOUND=',FOUND
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        LINE0=LINE2
        ROOT=.FALSE.

      ELSE IF(.NOT.ROOT) THEN
        CALL FIND(8,CDATA(2)(IBEG2:IEND2),LINE0,LINE_MAX,LINE1,FOUND,
     '    ERROR,*9999)
        IF(DOP) THEN
          WRITE(OP_STRING,*) ' LINE1=',LINE1,' FOUND=',FOUND
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        CALL FIND(8,CDATA(2)(IBEG2:IEND2),LINE1+1,LINE1+LINE_MAX,LINE2,
     '    FOUND,ERROR,*9999)
        IF(DOP) THEN
          WRITE(OP_STRING,*) ' LINE2=',LINE2,' FOUND=',FOUND
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF

      READ(8,FMT='(A)',REC=LINE1+1) CDATA(1)(1:80)
      IF(DOP) THEN
        WRITE(OP_STRING,'(1X,A)') CDATA(1)(1:80)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      IF(CDATA(1)(2:13).EQ.'FILE=Example') THEN !read com file
        READ_COM_FILE=.TRUE.
        CALL STRING_TRIM(CDATA(1),IBEG1,IEND1)
        CALL STRING_TRIM(EXAMPLES_DIR,IBEG2,IEND2)
        CALL OPENF(9,'DISK',EXAMPLES_DIR(IBEG2:IEND2)//CDATA(1)(7:25),
     '    'OLD','DIRECT','FORMATTED',132,ERROR,*9999)
        LINE1=1
        CALL FIND(9,'fem quit',1,LINE_MAX,LINE2,FOUND,ERROR,*9999)
        IF(DOP) THEN
          WRITE(OP_STRING,*) ' LINE2=',LINE2,' FOUND=',FOUND
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        NT_LINE=LINE2
        KEY32=.FALSE. !to disable numeric keys
      ELSE !continue reading current file
        READ_COM_FILE=.FALSE.
        NT_LINE=LINE2-LINE1-1
      ENDIF
      NT_WIND=1+INT((NT_LINE-1)/23)
      IF(DOP) THEN
        WRITE(OP_STRING,*) ' NT_WIND=',NT_WIND
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      LINE3=0
      NOWIND=1
      CONTINUE=.TRUE.
      DO WHILE(CONTINUE)
        DO i=1,MIN(23,NT_LINE-LINE3)
          IF(READ_COM_FILE) THEN
            READ(9,FMT='(A)',REC=LINE3+i) CDATA(1)(1:80)
          ELSE
            READ(8,FMT='(A)',REC=LINE1+LINE3+i) CDATA(1)(1:80)
          ENDIF
          IF(DOP) THEN
            WRITE(OP_STRING,'(1X,A)') CDATA(1)(1:80)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(USE_SOCKET) THEN
            CALL WRITES(IOH2,CDATA(1)(1:80),ERROR,*9999)
          ELSE
            P(1)=0.0
            P(2)=1.0-i*0.025
            CALL TEXT(1,2,8,CDATA(1)(1:80),P,ERROR,*9999)
          ENDIF
        ENDDO
        CALL GETSTR3(KEY,ERROR,*9999)
        IF(KEY.EQ.1.AND.NOWIND.GT.1) THEN           !Previous screen key
          NOWIND=NOWIND-1
        ELSE IF(KEY.EQ.2.AND.NOWIND.LT.NT_WIND) THEN !Next screen key
          NOWIND=NOWIND+1
        ELSE IF(KEY3.AND.KEY.EQ.3) THEN             !Do key
          DO_EXAMPLE=.TRUE.
          CONTINUE=.FALSE.
        ELSE IF(KEY4.AND.KEY.EQ.4) THEN             !Select key
          SELECT_EXAMPLE=.TRUE.
          CONTINUE=.FALSE.
        ELSE IF(KEY5.AND.KEY.EQ.5) THEN             !Insert key
          COPYFILE='EXAMPLE_'//TERM1//TERM2//TERM3
          CALL COPY_FILE(COPYFILE)
          WRITE(OP_STRING,
     '      '(1X,A,'' copied to local directory'')') COPYFILE
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ELSE IF(KEY.EQ.6) THEN                      !Delete key
          CDATA(2)=STRG
          IF(.NOT.USE_SOCKET) THEN
            CALL GKS_CLRWK(8,GALWAY,ERROR,*9999)
          ENDIF
          ROOT=.TRUE.
          OPTION1=.FALSE.
          OPTION2=.FALSE.
          OPTION3=.FALSE.
          OPTION4=.FALSE.
          OPTION5=.FALSE.
          TERM1=' '
          TERM2=' '
          TERM3=' '
          TERM4=' '
          TERM5=' '
          GO TO 10
        ELSE IF(KEY.EQ.7) THEN                      !Remove key
          CONTINUE=.FALSE.
        ELSE IF(KEY32.AND.KEY.GE.32) THEN           !Numeric key
          CALL STRING_TRIM(STRG,IBEG4,IEND4)
          IF(.NOT.OPTION1) THEN      !this is the 1st numeric key press
            OPTION1=.TRUE.
            TERM1=CHAR(KEY)
            CDATA(2)=STRG(IBEG4:IEND4)//' '//TERM1
            EXAMPLE_NAME='EXAMPLE_'//TERM1
          ELSE IF(.NOT.OPTION2) THEN !this is the 2nd numeric key press
            OPTION2=.TRUE.
            TERM2=CHAR(KEY)
            CDATA(2)=STRG(IBEG4:IEND4)//' '//TERM1//TERM2
            EXAMPLE_NAME='EXAMPLE_'//TERM1//TERM2
          ELSE IF(.NOT.OPTION3) THEN !this is the 3rd numeric key press
            OPTION3=.TRUE.
            TERM3=CHAR(KEY)
            CDATA(2)=STRG(IBEG4:IEND4)//' '//TERM1//TERM2//TERM3
            EXAMPLE_NAME='EXAMPLE_'//TERM1//TERM2//TERM3
          ELSE IF(.NOT.OPTION4) THEN !this is the 4th numeric key press
            OPTION4=.TRUE.
            TERM4=CHAR(KEY)
            CDATA(2)=STRG(IBEG4:IEND4)//' '//TERM1//TERM2//TERM3//TERM4
            EXAMPLE_NAME='EXAMPLE_'//TERM1//TERM2//TERM3//TERM4
          ELSE                       !this is the 5th numeric key press
            OPTION5=.TRUE.
            TERM5=CHAR(KEY)
            CDATA(2)=STRG(IBEG4:IEND4)//' '//TERM1//TERM2//TERM3//
     '        TERM4//TERM5
            EXAMPLE_NAME='EXAMPLE_'//TERM1//TERM2//TERM3//TERM4//TERM5
          ENDIF
          IF(.NOT.USE_SOCKET) THEN
            CALL GKS_CLRWK(8,GALWAY,ERROR,*9999)
          ENDIF
          GO TO 10
        ENDIF
        LINE3=23*(NOWIND-1)
        IF(.NOT.USE_SOCKET) THEN
          CALL GKS_CLRWK(8,GALWAY,ERROR,*9999)
        ENDIF
      ENDDO
      IF(.NOT.USE_SOCKET) THEN
        CALL DAWK(8,0,ERROR,*9999)
        CALL CLOSE_WS(8,ERROR,*9999)     !new
! old   CALL COWK(8,ERROR,*9999)     !GBS 3/8/92
      ENDIF
      CALL CLOSEF(8,ERROR,*9999)
      CALL CLOSEF(9,ERROR,*9999)
      DOCDIR=.FALSE.

      CALL EXITS('DOCUM')
      RETURN
 9999 CALL ERRORS('DOCUM',ERROR)
      CALL EXITS('DOCUM')
      CALL DAWK(8,0,ERROR,*9999)
      CALL CLOSE_WS(8,ERROR,*9999)     !new
! old CALL COWK(8,ERROR,*9999)     !GBS 3/8/92
      CLOSE(UNIT=8)
      CLOSE(UNIT=9)
      RETURN 1
      END


      SUBROUTINE ELLIPSE(SEMIAXIS1,SEMIAXIS2,XWC,YWC,NLINES,ERROR,*)

C#### Subroutine: ELLIPSE
C###  Description:
C**** Draws fill-area ellipse about the point XWC,YWC with semi-axis
C**** lengths SEMIAXIS1 and SEMIAXIS2 with NLINES lines (max 100).

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER NLINES
      REAL*8 SEMIAXIS1,SEMIAXIS2,XWC,YWC
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i
      REAL*8 D_ELLIPSE(3,100),THETA

      CALL ENTERS('ELLIPSE',*9999)
      DO i=1,NLINES
        THETA=2.0D0*PI*DBLE(i-1)/DBLE(NLINES)
        D_ELLIPSE(1,i)=XWC+SEMIAXIS1*DCOS(THETA)
        D_ELLIPSE(2,i)=YWC+SEMIAXIS2*DSIN(THETA)
      ENDDO
      CALL GKS_FA(1,2,NLINES,D_ELLIPSE,ERROR,*9999)

      CALL EXITS('ELLIPSE')
      RETURN
 9999 CALL ERRORS('ELLIPSE',ERROR)
      CALL EXITS('ELLIPSE')
      RETURN 1
      END


      SUBROUTINE EVENT(ID_WS,ID_DEVICE,ID_STATUS,CLASS,IDATA,R4DATA,
     '  SDATA,ERROR,*)

C#### Subroutine: EVENT
C###  Description:
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
C***           ID_STATUS  0 for input ok, 1 for error.
C***           IDATA      data input by choice
C***                      or segment and pick number for pick
C***                      or view trans.number for locator
C***           R4DATA     data input by valuator or locator
C***           SDATA       "     "    " string

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.FOR'
!     Parameter List
      INTEGER IDATA(*),ID_DEVICE,ID_STATUS,ID_WS,INPUT_CLASS,ISIZE
      REAL R4DATA(*)
      CHARACTER CLASS*(*),ERROR*(*),SDATA*(*)

      CALL ENTERS('EVENT',*9999)
C     get next event from event queue if present
      CALL GKS_WAIT(0.0,ID_WS,Input_Class,ID_Device,ERROR,*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,'('' ID_WS='',I4,'' ID_Device='',I4)')
     '    ID_WS,ID_Device
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      IF(INPUT_CLASS.EQ.GKS$K_INPUT_CLASS_CHOICE) THEN
        CLASS='CHOICE'
        CALL GKS_GET_CHOICE(ID_STATUS,IDATA,ERROR,*9999)
        ID_STATUS=1
      ELSE IF(INPUT_CLASS.EQ.GKS$K_INPUT_CLASS_VALUATOR) THEN
        CLASS='VALUATOR'
        CALL GKS_GET_VALUATOR(R4DATA,ERROR,*9999)
        ID_STATUS=1
      ELSE IF(INPUT_CLASS.EQ.GKS$K_INPUT_CLASS_LOCATOR) THEN
        CLASS='LOCATOR'
        CALL GKS_GET_LOCATOR(IDATA,R4DATA(1),R4DATA(2),ERROR,*9999)
        ID_STATUS=1
      ELSE IF(INPUT_CLASS.EQ.GKS$K_INPUT_CLASS_PICK) THEN
        CLASS='PICK'
        CALL GKS_GET_PICK(ID_STATUS,IDATA(1),IDATA(2),ERROR,*9999)
      ELSE IF(INPUT_CLASS.EQ.GKS$K_INPUT_CLASS_STRING) THEN
        CLASS='STRING'
        CALL GKS_GET_STRING(ID_STATUS,SDATA,ISIZE,ERROR,*9999) !this needs checking
      ELSE
        CLASS='NONE'
        IDATA(1)=0
      ENDIF
      IF(DOP) THEN
        WRITE(OP_STRING,'('' Input Class='',A)') CLASS
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('EVENT')
      RETURN
 9999 CALL ERRORS('EVENT',ERROR)
      CALL EXITS('EVENT')
      RETURN 1
      END


      SUBROUTINE FILL_AREA(IBUNDLE,iw,NTPTS,PTS,ERROR,*)

C#### Subroutine: FILL_AREA
C###  Description:
C**** Draws a fill-area on iw with index IBUNDLE.
C**** PTS(1..3,nopts) is a real*8 array containing the 3D coords of
C**** each point.
C**** If the IBUNDLE is 0 the primitive will use the previously
C**** defined fill-area index.
C**** NOTE: If iw is 15 or 16 (postscript) IBUNDLE is reset to be black.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:curr00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:map000.cmn'
!     Parameter List
      INTEGER IBUNDLE,iw,NTPTS
      REAL*8 PTS(3,*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nopts

      CALL ENTERS('FILL_AREA',*9999)
      IF(NTPTS.EQ.0) GOTO 9998

      IF(DOP) THEN
        WRITE(OP_STRING,'(/'' PTS(1,1..):'',4E12.3)')
     '    (PTS(1,nopts),nopts=1,4)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'( '' PTS(2,1..):'',4E12.3)')
     '    (PTS(2,nopts),nopts=1,4)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'( '' PTS(3,1..):'',4E12.3)')
     '    (PTS(3,nopts),nopts=1,4)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(iw.EQ.1) THEN !plot x against y
        IF(ibundle.ne.0) CALL GKS_SFAI(IBUNDLE,ERROR,*9999)
!       IF(NJT.EQ.2) THEN !plot y against x
          CALL GKS_FA(1,2,NTPTS,PTS,ERROR,*9999)
!       ELSE IF(NJT.EQ.3) THEN !plot z against x
!         CALL GKS_FA(1,3,NTPTS,PTS,ERROR,*9999)
!       ENDIF

      ELSE IF(iw.EQ.2) THEN !plot y against z
        IF(ibundle.ne.0) CALL GKS_SFAI(IBUNDLE,ERROR,*9999)
        CALL GKS_FA(2,3,NTPTS,PTS,ERROR,*9999)

      ELSE IF(iw.EQ.3) THEN !plot x against z
        IF(ibundle.ne.0) CALL GKS_SFAI(IBUNDLE,ERROR,*9999)
        CALL GKS_FA(1,3,NTPTS,PTS,ERROR,*9999)

!     ELSE IF(iw.EQ.3) THEN
C       IF(ibundle.ne.0) CALL PHIGS$SET_FAREA_INDEX(IBUNDLE)
!       CALL PHIGS$FILL_AREA_SET3(1,NTPTS,PTS)

      ELSE IF(iw.EQ.4) THEN
        IF(ibundle.ne.0) CALL GKS_SFAI(IBUNDLE,ERROR,*9999)
        IF(PROJEC(1:11).EQ.'RECTANGULAR') THEN !points x and y
        ELSE IF(PROJEC(1:2).EQ.'XI') THEN !points in Xi coords
          DO nopts=1,NTPTS
            PTS(1,nopts)=-1.0+2.*(REAL(MXI1-1)+PTS(1,nopts))/MAX_XI
            PTS(2,nopts)=-1.0+2.*(REAL(MXI2-1)+PTS(2,nopts))/MAX_XI
          ENDDO
        ELSE  !assume points are in polar coords
          CALL MAP4(NTPTS,PTS,PTS,ERROR,*9999)
        ENDIF
        CALL GKS_FA(1,2,NTPTS,PTS,ERROR,*9999)

      ELSE IF(iw.EQ.15) THEN !GKS postscript
        IF(ibundle.ne.0) THEN
C CPB 1/11/92 - Using Colour postscript so no need to change bundle ???
C          IF(IBUNDLE.GT.8) THEN !reset to black for postscript
C            IBUNDLE=IBUNDLE-8
C            IF(DOP) WRITE(IOOP,'('' Reset IBUNDLE='',I2)') IBUNDLE
C          ENDIF
C CPB 1/11/92 Changed call to GKS_SPLI (??? Lines in Fill Area) to
C GKS_SFAI
          CALL GKS_SFAI(IBUNDLE,ERROR,*9999)
        ENDIF
        IF(NJT.EQ.2) THEN !plot y against x
          CALL GKS_FA(1,2,NTPTS,PTS,ERROR,*9999)
        ELSE IF(NJT.EQ.3) THEN !plot z against x
          CALL GKS_FA(1,3,NTPTS,PTS,ERROR,*9999)
        ENDIF

      ELSE IF(iw.EQ.51) THEN !homogeneous strain
        IF(ibundle.ne.0) CALL GKS_SFAI(IBUNDLE,ERROR,*9999)
        CALL GKS_FA(1,2,NTPTS,PTS,ERROR,*9999)

      ENDIF

 9998 CALL EXITS('FILL_AREA')
      RETURN
 9999 CALL ERRORS('FILL_AREA',ERROR)
      CALL EXITS('FILL_AREA')
      RETURN 1
      END


      SUBROUTINE FLUSH_DEVICE_EVENTS(iw,CLASS,ERROR,*)

C#### Subroutine: FLUSH_DEVICE_EVENTS
C###  Description:
C**** Flush queues
C**** CLASS can be 'PICK', 'CHOICE', 'LOCATOR', 'VALUATOR' or 'ALL'

      IMPLICIT NONE
      INCLUDE 'SYS$LIBRARY:GKSDEFS.FOR'
!     Parameter List
      INTEGER iw
      CHARACTER CLASS*(*),ERROR*(*)

      CALL ENTERS('FLUSH_DEVICE_EVENTS',*9999)
      IF(iw.ne.3) THEN !gks
        IF(CLASS(1:3).EQ.'ALL') THEN
          CALL GKS_FLUSH(iw,GKS$K_Input_Class_Pick,1,ERROR,*9999)
          CALL GKS_FLUSH(iw,GKS$K_Input_Class_Choice,1,ERROR,*9999)
          CALL GKS_FLUSH(iw,GKS$K_Input_Class_Locator,1,ERROR,*9999)
          CALL GKS_FLUSH(iw,GKS$K_Input_Class_Valuator,1,ERROR,*9999)
        ELSE IF(CLASS(1:4).EQ.'PICK') THEN
          CALL GKS_FLUSH(iw,GKS$K_Input_Class_Pick,1,ERROR,*9999)
        ELSE IF(CLASS(1:6).EQ.'CHOICE') THEN
          CALL GKS_FLUSH(iw,GKS$K_Input_Class_Choice,1,ERROR,*9999)
        ELSE IF(CLASS(1:7).EQ.'LOCATOR') THEN
          CALL GKS_FLUSH(iw,GKS$K_Input_Class_Locator,1,ERROR,*9999)
        ELSE IF(CLASS(1:8).EQ.'VALUATOR') THEN
          CALL GKS_FLUSH(iw,GKS$K_Input_Class_Valuator,1,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('FLUSH_DEVICE_EVENTS')
      RETURN
 9999 CALL ERRORS('FLUSH_DEVICE_EVENTS',ERROR)
      CALL EXITS('FLUSH_DEVICE_EVENTS')
      RETURN 1
      END


      SUBROUTINE FREEWK(MNWK,MXWK,NOWK,ERROR,*)

C#### Subroutine: FREEWK
C###  Description:
C**** Finds the first free workstation NOWK between MNWK and MXWK.

      IMPLICIT NONE
!     Parameter List
      INTEGER MNWK,MXWK,NOWK
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER CONID,ERSTAT,n1wk,NO1WK,WKTY

      CALL ENTERS('FREEWK',*9999)
      CALL ASSERT((MNWK.GE.0).AND.(MXWK.LE.99),
     '  'Workstation limits are invalid',ERROR,*9999)
      DO 1 n1wk=MNWK,MXWK
        CALL GKS_QWKC(NO1WK,ERSTAT,CONID,WKTY,ERROR,*9999)
        IF(ERSTAT.EQ.25) THEN
          NOWK=n1wk
          GOTO 2
        ELSE IF(erstat.ne.0) THEN
          ERROR='Error occurred while searching for free workstation'
          GOTO 9999
        ENDIF
 1    CONTINUE
      ERROR='no free workstations exist between specified limits'
      GOTO 9999
 2    CALL EXITS('FREEWK')
      RETURN
 9999 CALL ERRORS('FREEWK',ERROR)
      CALL EXITS('FREEWK')
      RETURN 1
      END


      SUBROUTINE GKS_DRAW(iw,ISEG,CO,CSEG,STRING,ERROR,*)

C#### Subroutine: GKS_DRAW
C###  Description:
C**** Handles GKS draw functions with event mode input for iw.
C**** TRANSF(i,j,noseg) is segment transformation matrix
C**** LTRANS(noseg) is .TRUE. if a transf. has been defined previously
C**** Choice menus are:   97 = GKS_DRAW (event)
C****                     98 = GKS_DRAW (event)
C****                     99 = GKS_DRAW (request)
C**** Valuator is 81

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:draw00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'cmiss$reference:ntsg00.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.FOR'
!     Parameter List
      INTEGER ISEG(*),iw
      CHARACTER CO(*)*(*),CSEG(*)*(*),ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER AREA_COLOUR_INDEX,AREA_INDEX,AREA_INT_STYLE,
     '  AREA_STYLE_INDEX,i,IBEG,ID_DEVICE,ID_WS,IEND,
     '  index,INDEX_OLD,INPUT_CHOICE,INPUT_STATUS,INS,
     '  INSTAT,IPICK,ISDRAW(50),isegm,ITOT,LD1,
     '  LINE_COLOUR_INDEX,LINE_INDEX,LINE_TYPE,
     '  MARKER_COLOUR_INDEX,MARKER_INDEX,MARKER_TYPE,
     '  n1ch,NCHAR,NOCH,NOCO,nodraw,noseg,nosg,NTCH,NTDRAW,NTPTS,
     '  TEXT_HORIZ_ALIGN,TEXT_INDEX,TEXT_VERT_ALIGN,
     '  TEXT_COLOUR_INDEX,TEXT_FONT,TEXT_PATH,TEXT_PRECISION
      REAL CHSZ,ECAREA(4),LINE_WIDTH,MARKER_SIZE,
     '  R4DATA(10),
     '  TEXT_EXPFAC,TEXT_HEIGHT,TEXT_SPACING,VALUE,
     '  XREFA,XREFB,XVECT,
     '  YREFA,YREFB,YVECT
      REAL*8 A,B,DIST,D_PTS(2,500),D_XWC,D_YWC,D_XWC1,D_XWC2,D_YWC1,
     '  D_YWC2,RADIUS,RFROMC,THETA,THETA1,THETA2,X0,XSIDE,Y0,YSIDE
      CHARACTER CFROMI*4,CHOOSE*30,CIW*2,
     '  CLASS*8,OPTION0(50)*80,OPTION2(50)*80,
     '  OPTION3(50)*80,SDATA*(10),
     '  STRG*80,TEXT_STRING*80
      LOGICAL ABBREV,CONTINUEA,CONTINUEB,END,FIRST
      DATA LD1/1/

      CALL ENTERS('GKS_DRAW',*9999)
      OPTION0( 1)='Cell_Array_Attributes'
      OPTION0( 2)='Cell_Array'
      OPTION0( 3)='Fill_Area_Attributes'
      OPTION0( 4)='Fill_Area_Box'
      OPTION0( 5)='Fill_Area_Circle'
      OPTION0( 6)='Fill_Area_Ellipse'
      OPTION0( 7)='Polyline_Attributes'
      OPTION0( 8)='Polyline_Bezier'
      OPTION0( 9)='Polyline_Locate'
      OPTION0(10)='Polyline_Stroke'
      OPTION0(11)='Polyline_Arrow'
      OPTION0(12)='Polyline_Box'
      OPTION0(13)='Polyline_Box..'
      OPTION0(14)='Polyline_Circle'
      OPTION0(15)='Polyline_Circle..'
      OPTION0(16)='Polyline_Semicircle'
      OPTION0(17)='Polyline_Semicircle..'
      OPTION0(18)='Polyline_Ellipse'
      OPTION0(19)='Polyline_Semiellipse'
      OPTION0(20)='Polyline_Semiellipse..'
      OPTION0(21)='Polymarker_Attributes'
      OPTION0(22)='Polymarker_Locate'
      OPTION0(23)='Text_Attributes'
      OPTION0(24)='Text_Locate'
      OPTION0(25)='Segment_Remove'
      OPTION0(26)='Segment_Reposition'
      OPTION0(27)='Segment_Reposition..'
      OPTION0(28)='Segment_Resize'
      OPTION0(29)='Segment_Rotate'
      OPTION0(30)='Set_WC_Viewport'
      OPTION0(31)='Set_WC_Window'
      OPTION0(32)='Set_WS_Viewport'
      OPTION0(33)='Set_WS_Window'
      OPTION0(34)='Inquire'
      OPTION0(35)='Clear_WS'
      OPTION0(36)='Exit'
      NTCH=36
      IF(NJT.EQ.2) THEN
        CALL CHOICE('GKS',1,1,INS,97,'EVENT',NTCH,NTCH,NOCH,NOCO,9,
     '    CO,OPTION0,STRING,0.50*XDISP,YDISP,ERROR,*9999)
      ELSE IF(NJT.EQ.3) THEN
        CALL CHOICE('GKS',1,1,INS,97,'EVENT',NTCH,NTCH,NOCH,NOCO,2,
     '    CO,OPTION0,STRING,XDISP,0.,ERROR,*9999)
      ENDIF

      CHSZ=0.005
      AREA_COLOUR_INDEX=1
      AREA_INDEX=1
      AREA_INT_STYLE=GHOLLO
      AREA_STYLE_INDEX=1
      LINE_COLOUR_INDEX=1
      LINE_INDEX=1
      LINE_TYPE=GLSOLI
      LINE_WIDTH=1.0
      MARKER_INDEX=1
      MARKER_SIZE=1.0
      MARKER_TYPE=GPLUS
      MARKER_COLOUR_INDEX=1
      TEXT_COLOUR_INDEX=1
      TEXT_EXPFAC=1.0
      TEXT_FONT=1
      IF(DIAG.EQ.0.0) THEN
        TEXT_HEIGHT=0.02
      ELSE IF(DIAG.GT.0.0) THEN
        TEXT_HEIGHT=0.01*DIAG
      ENDIF
      TEXT_HORIZ_ALIGN=GAHNOR
      TEXT_VERT_ALIGN=GAVNOR
      TEXT_INDEX=1
      TEXT_PATH=GRIGHT
      TEXT_PRECISION=GSTRP
      TEXT_SPACING=0.0
      DO noseg=1,100
        LTRANS(noseg)=.FALSE.
      ENDDO
      NTDRAW=0 !initializes segment# for segments created in current session
      DO nodraw=1,50
        ISDRAW(nodraw)=0 !all segments will be new segments
      ENDDO

      WRITE(OP_STRING,'('' >>Choose from menu'')')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      CONTINUEA=.TRUE.
      DO WHILE (CONTINUEA)
        CALL EVENT(ID_WS,ID_Device,Input_Status,CLASS,Input_Choice,
     '    R4DATA,SDATA,ERROR,*9999)
        IF(DOP) THEN
          WRITE(OP_STRING,*) ' Input_Class=',Class
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        IF(CLASS(1:6).EQ.'CHOICE'.AND.ID_WS.EQ.97) THEN
          IF(INPUT_CHOICE.GT.0) THEN
            CHOOSE=OPTION0(Input_Choice)
            XREFA=0.8*XDISP
            YREFA=REAL(23-Input_Choice)*CHSZ
            XREFB=0.7*XDISP
            YREFB=YREFA
            IF(DOP) THEN
              WRITE(OP_STRING,*) ' CHOOSE=',CHOOSE
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            IF(CHOOSE(1:7).EQ.'Segment') THEN
              CALL GKS_SPKM(iw,LD1,GREQU,GECHO,ERROR,*9999)
              DO nosg=1,NTSG
                IF(ISEG(nosg).GT.0) THEN
                  CALL GKS_SDTEC(nosg,GDETEC,ERROR,*9999)
                ENDIF
              ENDDO
            ENDIF

            IF(ABBREV(CHOOSE,'CELL_ARRAY_ATTRIBUTES',21)) THEN
              OPTION2(1)='Colour_index'
              OPTION2(2)='Line_type'
              OPTION2(3)='Line_width'
              OPTION2(4)='Set_index'
              OPTION2(5)='Pick_segment'
              OPTION2(6)='Set_current'
              CALL CHOICE('CELL_ARRAY',1,1,INS,98,'EVENT',6,6,NOCH,
     '          NOCO,2,CO,OPTION2,STRING,XREFA,YREFA,ERROR,*9999)
              CONTINUEB=.TRUE.
              DO WHILE (CONTINUEB)
                CALL EVENT(ID_WS,ID_Device,Input_Status,CLASS,NOCH,
     '            R4DATA,SDATA,ERROR,*9999)
                IF(DOP) THEN
                  WRITE(OP_STRING,*)
     '              '        Input_Class=',Class
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                IF(CLASS(1:6).EQ.'CHOICE') THEN
                  IF(DOP) THEN
                    WRITE(OP_STRING,*)
     '                '        Input_Choice=',NOCH
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                  IF(NOCH.GT.0) THEN
                    IF(ABBREV(OPTION2(NOCH),'COLOUR_INDEX',8)) THEN
                      CALL COLOUR(LINE_COLOUR_INDEX,iw,CO,STRING,XREFA,
     '                  YREFA,ERROR,*9999)
                    ELSE IF(ABBREV(OPTION2(NOCH),'LINE_TYPE',8)) THEN
                      OPTION3(1)='Solid'
                      OPTION3(2)='Dashed'
                      OPTION3(3)='Dotted'
                      OPTION3(4)='Dash_dot'
                      CALL CHOICE('PL',1,1,INS,99,'REQUEST',4,4,NOCH,
     '                  NOCO,2,CO,OPTION3,STRING,0.7*XDISP,YREFA,
     '                  ERROR,*9999)
                      IF(ABBREV(OPTION3(NOCH),'SOLID',        5)) THEN
                        LINE_TYPE=GLSOLI
                      ELSE IF(ABBREV(OPTION3(NOCH),'DASHED',  5)) THEN
                        LINE_TYPE=GLDASH
                      ELSE IF(ABBREV(OPTION3(NOCH),'DOTTED',  5)) THEN
                        LINE_TYPE=GLDOT
                      ELSE IF(ABBREV(OPTION3(NOCH),'DASH_DOT',5)) THEN
                        LINE_TYPE=GLDASD
                      ENDIF
                      IF(DOP) THEN
                        WRITE(OP_STRING,
     '                    '('' LINE_TYPE='',I1)') LINE_TYPE
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      ENDIF
                    ELSE IF(ABBREV(OPTION2(NOCH),'LINE_WIDTH',8)) THEN
                      CALL VALUATOR('LINE_WIDTH',81,'REQUEST',2,0.,5.,
     '                  1.,VALUE,XREFA,YREFA,ERROR,*9999)
                    ELSE IF(ABBREV(OPTION2(NOCH),'PICK_SEGMENT',8)) THEN
                      DO nosg=1,NTSG
                        IF(ISEG(nosg).GT.0)
     '                    CALL DETECT(iw,ISEG,nosg,'DETECTABLE',
     '                      ERROR,*9999)
                      ENDDO
                      WRITE(OP_STRING,
     '                  '('' >>Pick segment to change '//
     '                  '(use 2nd button to end)'')')
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      CALL PICK(iw,LD1,'REQUEST',INSTAT,isegm,IPICK,
     '                  ERROR,*9999)
                      DO WHILE (INSTAT.EQ.1)
                        WRITE(OP_STRING,'('' Segment number='',I3,'
     '                    //''' Type='',A10)') isegm,CSEG(isegm)
                        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                        CALL GKS_SPLR(iw,isegm,LINE_TYPE,LINE_WIDTH,
     '                    LINE_COLOUR_INDEX,ERROR,*9999)
                        CALL GKS_UWK(iw,GPERFO,ERROR,*9999)
                        CALL PICK(iw,LD1,'REQUEST',INSTAT,isegm,IPICK,
     '                    ERROR,*9999)
                      ENDDO
                      CONTINUEB=.FALSE.
                    ELSE IF(ABBREV(OPTION2(NOCH),'SET_CURRENT',8)) THEN
                      CONTINUEB=.FALSE.
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO
              CALL INPUT_MODE(98,1,'CHOICE','REQUEST',ERROR,*9999)

            ELSE IF(ABBREV(CHOOSE,'CELL_ARRAY',10)) THEN

            ELSE IF(ABBREV(CHOOSE,'FILL_AREA_ATTRIBUTES',20)) THEN
              OPTION2(1)='Colour_index'
              OPTION2(2)='Int_Style'
              OPTION2(3)='Reference_pt'
              OPTION2(4)='Pattern_size'
              OPTION2(5)='Set_index'
              OPTION2(6)='Pick_segment'
              OPTION2(7)='Set_current'
              CALL CHOICE('FILL_AREA',1,1,INS,98,'EVENT',7,7,NOCH,NOCO,
     '          2,CO,OPTION2,STRING,XREFA,YREFA,ERROR,*9999)
              CONTINUEB=.TRUE.
              DO WHILE (CONTINUEB)
                CALL EVENT(ID_WS,ID_Device,Input_Status,CLASS,NOCH,
     '            R4DATA,SDATA,ERROR,*9999)
                IF(DOP) THEN
                  WRITE(OP_STRING,*)
     '              '        Input_Class=',Class
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                IF(CLASS(1:6).EQ.'CHOICE') THEN
                  IF(DOP) THEN
                    WRITE(OP_STRING,*)
     '                '        Input_Choice=',NOCH
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                  IF(NOCH.GT.0) THEN
                    IF(ABBREV(OPTION2(NOCH),'COLOUR_INDEX',8)) THEN
                      CALL COLOUR(AREA_COLOUR_INDEX,iw,CO,STRING,XREFA,
     '                  YREFA,ERROR,*9999)
                    ELSE IF(ABBREV(OPTION2(NOCH),'INT_STYLE',8)) THEN
                      OPTION3(1)='Hollow'
                      OPTION3(2)='Solid'
                      OPTION3(3)='Pattern'
                      OPTION3(4)='Hatch'
                      CALL CHOICE('INT_STYLE',1,1,INS,99,'REQUEST',4,4,
     '                  NOCH,NOCO,2,
     '                  CO,OPTION3,STRING,0.7*XDISP,YREFA,ERROR,*9999)
                      IF(ABBREV(OPTION3(NOCH),'HOLLOW',5)) THEN
                        AREA_INT_STYLE=GHOLLO
                      ELSE IF(ABBREV(OPTION3(NOCH),'SOLID'  ,5)) THEN
                        AREA_INT_STYLE=GSOLID
                      ELSE IF(ABBREV(OPTION3(NOCH),'PATTERN',5)) THEN
                        AREA_INT_STYLE=GPATTR
                      ELSE IF(ABBREV(OPTION3(NOCH),'HATCH'  ,5)) THEN
                        AREA_INT_STYLE=GHATCH
                      ENDIF
                    ELSE IF(ABBREV(OPTION2(NOCH),'REFERENCE_PT',8)) THEN
                    ELSE IF(ABBREV(OPTION2(NOCH),'PATTERN_SIZE',8)) THEN
                      CALL VALUATOR('PATTERN_SIZE',81,'REQUEST',2,0.,
     '                  5.,1.,VALUE,XREFB,YREFB,ERROR,*9999)
                    ELSE IF(ABBREV(OPTION2(NOCH),'SET_INDEX',8)) THEN
                      OPTION3(1)='Index_1'
                      OPTION3(2)='Index_2'
                      OPTION3(3)='Index_3'
                      OPTION3(4)='Index_4'
                      CALL CHOICE('SET_INDEX',1,1,INS,99,'REQUEST',4,4,
     '                  NOCH,NOCO,2,
     '                  CO,OPTION3,STRING,0.7*XDISP,YREFA,ERROR,*9999)
                      IF(ABBREV(OPTION3(NOCH),'INDEX_1',7)) THEN
                        AREA_INDEX=1
                      ELSE IF(ABBREV(OPTION3(NOCH),'INDEX_2',7)) THEN
                        AREA_INDEX=2
                      ELSE IF(ABBREV(OPTION3(NOCH),'INDEX_3',7)) THEN
                        AREA_INDEX=3
                      ELSE IF(ABBREV(OPTION3(NOCH),'INDEX_4',7)) THEN
                        AREA_INDEX=4
                      ENDIF
                    ELSE IF(ABBREV(OPTION2(NOCH),'PICK_SEGMENT',8)) THEN
                      DO nosg=1,NTSG
                        IF(ISEG(nosg).GT.0.AND.CSEG(nosg)(5:8).EQ.
     '                    'area')
     '                    CALL DETECT(iw,ISEG,nosg,'DETECTABLE',
     '                      ERROR,*9999)
                      ENDDO
                      WRITE(OP_STRING,'('' >>Pick segment '
     '                  //'to change (use 2nd button to end)'')')
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      CALL PICK(iw,LD1,'REQUEST',INSTAT,isegm,IPICK,
     '                  ERROR,*9999)
                      DO WHILE (INSTAT.EQ.1)
                        WRITE(OP_STRING,'('' Segment number='',I3,'
     '                    //''' Type='',A10)') isegm,CSEG(isegm)
                        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                        CALL GKS_SFAR(iw,isegm,AREA_INT_STYLE,
     '                    AREA_STYLE_INDEX,AREA_COLOUR_INDEX,
     '                    ERROR,*9999)
                        CALL GKS_UWK(iw,GPERFO,ERROR,*9999)
                        CALL PICK(iw,LD1,'REQUEST',INSTAT,isegm,IPICK,
     '                    ERROR,*9999)
                      ENDDO
                      DO nosg=1,NTSG
                        IF(ISEG(nosg).GT.0.AND.CSEG(nosg)(5:8).EQ.
     '                    'area')
     '                    CALL DETECT(iw,ISEG,nosg,'UNDETECTABLE',
     '                      ERROR,*9999)
                      ENDDO
                      CONTINUEB=.FALSE.
                    ELSE IF(ABBREV(OPTION2(NOCH),'SET_CURRENT',8)) THEN
                      CALL GKS_SFAI(AREA_INDEX,ERROR,*9999)
                      CALL GKS_SFAR(iw,AREA_INDEX,AREA_INT_STYLE,
     '                  AREA_STYLE_INDEX,AREA_COLOUR_INDEX,ERROR,*9999)
                      CONTINUEB=.FALSE.
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO
              CALL INPUT_MODE(98,1,'CHOICE','REQUEST',ERROR,*9999)

            ELSE IF(ABBREV(CHOOSE,'FILL_AREA_BOX',13)) THEN
              CALL ACWK(iw,0,ERROR,*9999)
              NTDRAW=NTDRAW+1
              CALL OPEN_SEGMENT(ISDRAW(NTDRAW),ISEG,iw,
     '          'DRAW fill_area_box',index,INDEX_OLD,0,1,
     '          CSEG,ERROR,*9999)
              CALL GKS_SFAI(AREA_INDEX,ERROR,*9999)
              CALL GKS_SFAR(iw,AREA_INDEX,AREA_INT_STYLE,
     '          AREA_STYLE_INDEX,AREA_COLOUR_INDEX,ERROR,*9999)
              INSTAT=1
              DO WHILE (INSTAT.EQ.1)
                WRITE(OP_STRING,'('' >>Locate first corner of box '//
     '            '(use 2nd button to end)'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                CALL LOCATOR(1,iw,INSTAT,'REQUEST',3,
     '            0.0D0,D_XWC1,0.0D0,D_YWC1,ERROR,*9999)
                IF(INSTAT.EQ.1) THEN
                  WRITE(OP_STRING,
     '              '('' >>Locate second corner of box'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  CALL LOCATOR(1,iw,INSTAT,'REQUEST',5,
     '              D_XWC1,D_XWC2,D_YWC1,D_YWC2,ERROR,*9999)
                  D_XWC=0.5D0*(D_XWC1+D_XWC2)
                  D_YWC=0.5D0*(D_YWC1+D_YWC2)
                  XSIDE=DABS(D_XWC2-D_XWC1)/2.0D0
                  YSIDE=DABS(D_YWC2-D_YWC1)/2.0D0
                  D_PTS(1,1)=D_XWC-XSIDE
                  D_PTS(2,1)=D_YWC-YSIDE
                  D_PTS(1,2)=D_XWC+XSIDE
                  D_PTS(2,2)=D_PTS(2,1)
                  D_PTS(1,3)=D_PTS(1,2)
                  D_PTS(2,3)=D_YWC+YSIDE
                  D_PTS(1,4)=D_PTS(1,1)
                  D_PTS(2,4)=D_PTS(2,3)
                  D_PTS(1,5)=D_PTS(1,1)
                  D_PTS(2,5)=D_PTS(2,1)
                  CALL GKS_FA(1,2,5,D_PTS,ERROR,*9999)
                  XSEGMENT_DATA(1,NTSG)=REAL(D_XWC)
                  YSEGMENT_DATA(1,NTSG)=REAL(D_YWC)
                  XSEGMENT_DATA(2,NTSG)=REAL(XSIDE)
                  YSEGMENT_DATA(2,NTSG)=REAL(YSIDE)
                ENDIF
              ENDDO
              CALL CLOSE_SEGMENT(ISDRAW(NTDRAW),iw,ERROR,*9999)
              CALL DAWK(iw,0,ERROR,*9999)

            ELSE IF(ABBREV(CHOOSE,'FILL_AREA_CIRCLE',16)) THEN
              CALL ACWK(iw,0,ERROR,*9999)
              NTDRAW=NTDRAW+1
              CALL OPEN_SEGMENT(ISDRAW(NTDRAW),ISEG,iw,
     '          'DRAW fill_area_circle',index,INDEX_OLD,0,1,
     '          CSEG,ERROR,*9999)
              CALL GKS_SFAI(AREA_INDEX,ERROR,*9999)
              CALL GKS_SFAR(iw,AREA_INDEX,AREA_INT_STYLE,
     '          AREA_STYLE_INDEX,AREA_COLOUR_INDEX,ERROR,*9999)
              INSTAT=1
              DO WHILE (INSTAT.EQ.1)
                WRITE(OP_STRING,'('' >>Locate centre of circle '//
     '            '(use 2nd button to end)'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                CALL LOCATOR(1,iw,INSTAT,'REQUEST',3,
     '            0.0D0,D_XWC1,0.0D0,D_YWC1,ERROR,*9999)
                IF(INSTAT.EQ.1) THEN
                  WRITE(OP_STRING,
     '              '('' >>Locate radius of circle'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  CALL LOCATOR(1,iw,INSTAT,'REQUEST',4,
     '              D_XWC1,D_XWC2,D_YWC1,D_YWC2,ERROR,*9999)
                  RADIUS=DSQRT((D_XWC2-D_XWC1)**2+(D_YWC2-D_YWC1)**2)
                  DO i=1,100
                    THETA=2.0D0*PI*DBLE(i-1)/100.0D0
                    D_PTS(1,i)=D_XWC1+RADIUS*DCOS(THETA)
                    D_PTS(2,i)=D_YWC1+RADIUS*DSIN(THETA)
                  ENDDO
                  CALL GKS_FA(1,2,100,D_PTS,ERROR,*9999)
                  XSEGMENT_DATA(1,NTSG)=REAL(D_XWC1)
                  YSEGMENT_DATA(1,NTSG)=REAL(D_YWC1)
                  XSEGMENT_DATA(2,NTSG)=REAL(RADIUS)
                  XSEGMENT_DATA(3,NTSG)=0.0
                  YSEGMENT_DATA(3,NTSG)=REAL(2.0D0*PI)
                ENDIF
              ENDDO
              CALL CLOSE_SEGMENT(ISDRAW(NTDRAW),iw,ERROR,*9999)
              CALL DAWK(iw,0,ERROR,*9999)

            ELSE IF(ABBREV(CHOOSE,'FILL_AREA_ELLIPSE',16)) THEN
              CALL ACWK(iw,0,ERROR,*9999)
              NTDRAW=NTDRAW+1
              CALL OPEN_SEGMENT(ISDRAW(NTDRAW),ISEG,iw,
     '          'DRAW fill_area_ellipse',index,INDEX_OLD,0,1,
     '          CSEG,ERROR,*9999)
              CALL GKS_SFAI(AREA_INDEX,ERROR,*9999)
              CALL GKS_SFAR(iw,AREA_INDEX,AREA_INT_STYLE,
     '          AREA_STYLE_INDEX,AREA_COLOUR_INDEX,ERROR,*9999)
              INSTAT=1
              DO WHILE (INSTAT.EQ.1)
                WRITE(OP_STRING,
     '            '('' >>Locate first corner of enclosing'//
     '            ' box (use 2nd button to end)'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                CALL LOCATOR(1,iw,INSTAT,'REQUEST',3,
     '            0.0D0,D_XWC1,0.0D0,D_YWC1,ERROR,*9999)
                IF(INSTAT.EQ.1) THEN
                  WRITE(OP_STRING,
     '              '('' >>Locate second corner of enclosing'//
     '              ' box'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  CALL LOCATOR(1,iw,INSTAT,'REQUEST',5,
     '              D_XWC1,D_XWC2,D_YWC1,D_YWC2,ERROR,*9999)
                  X0=0.5D0*(D_XWC1+D_XWC2)
                  Y0=0.5D0*(D_YWC1+D_YWC2)
                  A=DABS(D_XWC2-D_XWC1)/2.0D0
                  B=DABS(D_YWC2-D_YWC1)/2.0D0
                  DO I=1,100
                    THETA=2.0D0*PI*DBLE(I-1)/99.0D0
                    D_PTS(1,I)=X0+A*DCOS(THETA)
                    D_PTS(2,I)=Y0+B*DSIN(THETA)
                  ENDDO
                  CALL GKS_FA(1,2,100,D_PTS,ERROR,*9999)
                  XSEGMENT_DATA(1,NTSG)=REAL(X0)
                  YSEGMENT_DATA(1,NTSG)=REAL(Y0)
                  XSEGMENT_DATA(2,NTSG)=0.0
                  XSEGMENT_DATA(3,NTSG)=0.0
                  YSEGMENT_DATA(3,NTSG)=REAL(2.0D0*PI)
                ENDIF
              ENDDO
              CALL CLOSE_SEGMENT(ISDRAW(NTDRAW),iw,ERROR,*9999)
              CALL DAWK(iw,0,ERROR,*9999)

            ELSE IF(ABBREV(CHOOSE,'POLYLINE_ATTRIBUTES',19)) THEN
              OPTION2(1)='Colour_index'
              OPTION2(2)='Line_type'
              OPTION2(3)='Line_width'
              OPTION2(4)='Set_index'
              OPTION2(5)='Pick_segment'
              OPTION2(6)='Set_current'
              CALL CHOICE('POLYLINE',1,1,INS,98,'EVENT',6,6,n1ch,NOCO,2,
     '          CO,OPTION2,STRING,XREFA,YREFA,ERROR,*9999)
              CONTINUEB=.TRUE.
              DO WHILE (CONTINUEB)
                CALL EVENT(ID_WS,ID_Device,Input_Status,CLASS,NOCH,
     '            R4DATA,SDATA,ERROR,*9999)
                IF(DOP) THEN
                  WRITE(OP_STRING,*)
     '              '        Input_Class=',Class
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                IF(CLASS(1:6).EQ.'CHOICE') THEN
                  IF(DOP) THEN
                    WRITE(OP_STRING,*)
     '                '        Input_Choice=',NOCH
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                  IF(NOCH.GT.0) THEN
                    IF(ABBREV(OPTION2(NOCH),'COLOUR_INDEX',8)) THEN
                      CALL COLOUR(LINE_COLOUR_INDEX,iw,CO,STRING,XREFA,
     '                  YREFA,ERROR,*9999)
                    ELSE IF(ABBREV(OPTION2(NOCH),'LINE_TYPE',8)) THEN
                      OPTION3(1)='Solid'
                      OPTION3(2)='Dashed'
                      OPTION3(3)='Dotted'
                      OPTION3(4)='Dash_dot'
                      CALL CHOICE('LINE_TYPE',1,1,INS,99,'REQUEST',4,4,
     '                  NOCH,NOCO,2,CO,OPTION3,STRING,0.7*XDISP,YREFA,
     '                  ERROR,*9999)
                      IF(ABBREV(OPTION3(NOCH),'SOLID',        5)) THEN
                        LINE_TYPE=GLSOLI
                      ELSE IF(ABBREV(OPTION3(NOCH),'DASHED',  5)) THEN
                        LINE_TYPE=GLDASH
                      ELSE IF(ABBREV(OPTION3(NOCH),'DOTTED',  5)) THEN
                        LINE_TYPE=GLDOT
                      ELSE IF(ABBREV(OPTION3(NOCH),'DASH_DOT',5)) THEN
                        LINE_TYPE=GLDASD
                      ENDIF
                      IF(DOP) THEN
                        WRITE(OP_STRING,
     '                    '('' LINE_TYPE='',I1)') LINE_TYPE
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      ENDIF
                    ELSE IF(ABBREV(OPTION2(NOCH),'LINE_WIDTH',8)) THEN
                      CALL VALUATOR('LINE_WIDTH',81,'REQUEST',2,0.,5.,
     '                  1.,LINE_WIDTH,XREFB,YREFB,ERROR,*9999)
                    ELSE IF(ABBREV(OPTION2(NOCH),'SET_INDEX',8)) THEN
                      OPTION3(1)='Index_1'
                      OPTION3(2)='Index_2'
                      OPTION3(3)='Index_3'
                      OPTION3(4)='Index_4'
                      CALL CHOICE('SET_INDEX',1,1,INS,99,'REQUEST',4,4,
     '                  NOCH,NOCO,2,CO,OPTION3,STRING,0.7*XDISP,YREFA,
     '                  ERROR,*9999)
                      IF(ABBREV(OPTION3(NOCH),'INDEX_1',7)) THEN
                        LINE_INDEX=1
                      ELSE IF(ABBREV(OPTION3(NOCH),'INDEX_2',7)) THEN
                        LINE_INDEX=2
                      ELSE IF(ABBREV(OPTION3(NOCH),'INDEX_3',7)) THEN
                        LINE_INDEX=3
                      ELSE IF(ABBREV(OPTION3(NOCH),'INDEX_4',7)) THEN
                        LINE_INDEX=4
                      ENDIF
                    ELSE IF(ABBREV(OPTION2(NOCH),'PICK_SEGMENT',8)) THEN
                      DO nosg=1,NTSG
                        IF(ISEG(nosg).GT.0.AND.CSEG(nosg)(5:8).EQ.
     '                    'line')
     '                    CALL DETECT(iw,ISEG,nosg,'DETECTABLE',
     '                      ERROR,*9999)
                      ENDDO
                      WRITE(OP_STRING,
     '                  '('' >>Pick segment to change '//
     '                  '(use 2nd button to end)'')')
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      CALL PICK(iw,LD1,'REQUEST',INSTAT,isegm,IPICK,
     '                  ERROR,*9999)
                      DO WHILE (INSTAT.EQ.1)
                        WRITE(OP_STRING,'('' Segment number='',I3,'
     '                    //''' Type='',A10)') isegm,CSEG(isegm)
                        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                        CALL GKS_SPLR(iw,isegm,LINE_TYPE,LINE_WIDTH,
     '                    LINE_COLOUR_INDEX,ERROR,*9999)
                        CALL GKS_UWK(iw,GPERFO,ERROR,*9999)
                        CALL PICK(iw,LD1,'REQUEST',INSTAT,isegm,IPICK,
     '                    ERROR,*9999)
                      ENDDO
                      DO nosg=1,NTSG
                        IF(ISEG(nosg).GT.0.AND.CSEG(nosg)(5:8).EQ.
     '                    'line')
     '                    CALL DETECT(iw,ISEG,nosg,'UNDETECTABLE',
     '                      ERROR,*9999)
                      ENDDO
                      CONTINUEB=.FALSE.
                    ELSE IF(ABBREV(OPTION2(NOCH),'SET_CURRENT',8)) THEN
                      CALL GKS_SPLI(LINE_INDEX,ERROR,*9999)
                      CALL GKS_SPLR(iw,LINE_INDEX,LINE_TYPE,LINE_WIDTH,
     '                  LINE_COLOUR_INDEX,ERROR,*9999)
                      CONTINUEB=.FALSE.
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO
              CALL INPUT_MODE(98,1,'CHOICE','REQUEST',ERROR,*9999)

            ELSE IF(ABBREV(CHOOSE,'POLYLINE_BEZIER',15)) THEN
              CIW=CFROMI(iw,'(I1)')
              STRG='FEM define polyline;m bezier on '//CIW
              CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*9999)

            ELSE IF(ABBREV(CHOOSE,'POLYLINE_LOCATE',15)) THEN
              CALL ACWK(iw,0,ERROR,*9999)
              NTDRAW=NTDRAW+1
              CALL OPEN_SEGMENT(ISDRAW(NTDRAW),ISEG,iw,
     '          'DRAW polyline_locate',index,INDEX_OLD,0,1,
     '          CSEG,ERROR,*9999)
              CALL GKS_SPLI(LINE_INDEX,ERROR,*9999)
              CALL GKS_SPLR(iw,LINE_INDEX,LINE_TYPE,LINE_WIDTH,
     '          LINE_COLOUR_INDEX,ERROR,*9999)
              WRITE(OP_STRING,
     '          '('' >>Locate start of polyline '')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              CALL LOCATOR(1,iw,INSTAT,'REQUEST',3,
     '          0.0D0,D_PTS(1,1),0.0D0,D_PTS(2,1),ERROR,*9999)
              XSEGMENT_DATA(1,NTSG)=REAL(D_PTS(1,1))
              YSEGMENT_DATA(1,NTSG)=REAL(D_PTS(2,1))
              WRITE(OP_STRING,
     '          '('' >>Locate points on polyline '//
     '          '(use 2nd button to end)'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              NTPTS=0
              DO WHILE (INSTAT.EQ.1.AND.NTPTS.LT.20)
                CALL LOCATOR(1,iw,INSTAT,'REQUEST',4,
     '            D_PTS(1,1),D_PTS(1,2),D_PTS(2,1),D_PTS(2,2),
     '            ERROR,*9999)
                IF(INSTAT.EQ.1) THEN
                  NTPTS=NTPTS+1
                  CALL GKS_PL(1,2,2,D_PTS,ERROR,*9999)
                  XSEGMENT_DATA(NTPTS+1,NTSG)=REAL(D_PTS(1,2))
                  YSEGMENT_DATA(NTPTS+1,NTSG)=REAL(D_PTS(2,2))
                  D_PTS(1,1)=D_PTS(1,2)
                  D_PTS(2,1)=D_PTS(2,2)
                ENDIF
              ENDDO
              ISEGMENT_DATA(1,NTSG)=NTPTS
              CALL CLOSE_SEGMENT(ISDRAW(NTDRAW),iw,ERROR,*9999)
              CALL DAWK(iw,0,ERROR,*9999)

            ELSE IF(ABBREV(CHOOSE,'POLYLINE_STROKE',15)) THEN
              CALL ACWK(iw,0,ERROR,*9999)
              NTDRAW=NTDRAW+1
              CALL OPEN_SEGMENT(ISDRAW(NTDRAW),ISEG,iw,
     '          'DRAW polyline_stroke',index,INDEX_OLD,0,1,
     '          CSEG,ERROR,*9999)
              CALL GKS_SPLI(LINE_INDEX,ERROR,*9999)
              CALL GKS_SPLR(iw,LINE_INDEX,LINE_TYPE,LINE_WIDTH,
     '          LINE_COLOUR_INDEX,ERROR,*9999)
              ECAREA(1)=0.0
              ECAREA(2)=XDISP
              ECAREA(3)=0.0
              ECAREA(4)=YDISP
              CALL GKS_INLC(iw,LD1,iw,D_XWC,D_YWC,0,ERROR,*9999)
              INSTAT=1
              DO WHILE (INSTAT.EQ.1)
                WRITE(OP_STRING,'('' >>Locate initial '//
     '            'stroke point (use 2nd button to end)'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                CALL LOCATOR(1,iw,INSTAT,'REQUEST',3,0.0D0,D_XWC,
     '            0.0D0,D_YWC,ERROR,*9999)
                IF(INSTAT.EQ.1) THEN
                  WRITE(OP_STRING,
     '              '('' >>Use 1st button to end stroke '//
     '              'input'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

                  CALL STROKE(1,iw,INSTAT,'REQUEST',NTPTS,
     '              ERROR,*9999)
                  WRITE(OP_STRING,*) ' Number of stroke points = ',
     '              NTPTS
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  CALL GKS_PL(1,2,NTPTS,D_PTS,ERROR,*9999)
                  XSEGMENT_DATA(1,NTSG)=REAL(D_XWC)
                  YSEGMENT_DATA(1,NTSG)=REAL(D_YWC)
                ENDIF
              ENDDO
              CALL CLOSE_SEGMENT(ISDRAW(NTDRAW),iw,ERROR,*9999)
              CALL DAWK(iw,0,ERROR,*9999)

            ELSE IF(ABBREV(CHOOSE,'POLYLINE_ARROW',14)) THEN
              CALL ACWK(iw,0,ERROR,*9999)
              NTDRAW=NTDRAW+1
              CALL OPEN_SEGMENT(ISDRAW(NTDRAW),ISEG,iw,
     '          'DRAW polyline_arrow',index,INDEX_OLD,0,1,
     '          CSEG,ERROR,*9999)
              CALL GKS_SPLI(LINE_INDEX,ERROR,*9999)
              CALL GKS_SPLR(iw,LINE_INDEX,LINE_TYPE,LINE_WIDTH,
     '          LINE_COLOUR_INDEX,ERROR,*9999)
              AREA_INT_STYLE=GSOLID
              CALL GKS_SFAI(AREA_INDEX,ERROR,*9999)
              CALL GKS_SFAR(iw,AREA_INDEX,AREA_INT_STYLE,
     '          AREA_STYLE_INDEX,AREA_COLOUR_INDEX,ERROR,*9999)
              INSTAT=1
              DO WHILE (INSTAT.EQ.1)
                WRITE(OP_STRING,'('' >>Locate arrow head '//
     '            '(use 2nd button to end)'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                CALL LOCATOR(1,iw,INSTAT,'REQUEST',3,0.0D0,D_PTS(1,1),
     '            0.0D0,D_PTS(2,1),ERROR,*9999)
                IF(INSTAT.EQ.1) THEN
                  WRITE(OP_STRING,'('' >>Locate arrow tail'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  CALL LOCATOR(1,iw,INSTAT,'REQUEST',4,
     '              D_PTS(1,1),D_PTS(1,2),D_PTS(2,1),D_PTS(2,2),
     '              ERROR,*9999)
                  CALL GKS_PL(1,2,2,D_PTS,ERROR,*9999)
                  THETA=DATAN2(D_PTS(2,2)-D_PTS(2,1),
     '                        D_PTS(1,2)-D_PTS(1,1))
                  DIST =DSQRT((D_PTS(1,2)-D_PTS(1,1))**2
     '              +(D_PTS(2,2)-D_PTS(2,1))**2)
                  D_PTS(1,2)=D_PTS(1,1)+DIST/3.0D0*DCOS(THETA-PI/10.0D0)
                  D_PTS(2,2)=D_PTS(2,1)+DIST/3.0D0*DSIN(THETA-PI/10.0D0)
                  D_PTS(1,3)=D_PTS(1,1)+DIST/3.0D0*DCOS(THETA+PI/10.0D0)
                  D_PTS(2,3)=D_PTS(2,1)+DIST/3.0D0*DSIN(THETA+PI/10.0D0)
                  CALL GKS_FA(1,2,3,D_PTS,ERROR,*9999)
                ENDIF
              ENDDO
              XSEGMENT_DATA(1,NTSG)=REAL(D_PTS(1,1))
              YSEGMENT_DATA(1,NTSG)=REAL(D_PTS(2,1))
              CALL CLOSE_SEGMENT(ISDRAW(NTDRAW),iw,ERROR,*9999)
              CALL DAWK(iw,0,ERROR,*9999)

            ELSE IF(ABBREV(CHOOSE,'POLYLINE_BOX',12)) THEN
              CALL ACWK(iw,0,ERROR,*9999)
              NTDRAW=NTDRAW+1
              CALL OPEN_SEGMENT(ISDRAW(NTDRAW),ISEG,iw,
     '          'DRAW polyline_box',index,INDEX_OLD,0,1,
     '          CSEG,ERROR,*9999)
              CALL GKS_SPLI(LINE_INDEX,ERROR,*9999)
              CALL GKS_SPLR(iw,LINE_INDEX,LINE_TYPE,LINE_WIDTH,
     '          LINE_COLOUR_INDEX,ERROR,*9999)
              INSTAT=1
              DO WHILE (INSTAT.EQ.1)
                WRITE(OP_STRING,
     '            '('' >>Locate first corner of box '//
     '            '(use 2nd button to end)'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                CALL LOCATOR(1,iw,INSTAT,'REQUEST',3,
     '            0.0D0,D_XWC1,0.0D0,D_YWC1,ERROR,*9999)
                IF(INSTAT.EQ.1) THEN
                  WRITE(OP_STRING,
     '              '('' >>Locate second corner of box'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  CALL LOCATOR(1,iw,INSTAT,'REQUEST',5,
     '              D_XWC1,D_XWC2,D_YWC1,D_YWC2,ERROR,*9999)
                  D_XWC=0.5D0*(D_XWC1+D_XWC2)
                  D_YWC=0.5D0*(D_YWC1+D_YWC2)
                  XSIDE=DABS(D_XWC2-D_XWC1)/2.0D0
                  YSIDE=DABS(D_YWC2-D_YWC1)/2.0D0
                  D_PTS(1,1)=D_XWC-XSIDE
                  D_PTS(2,1)=D_YWC-YSIDE
                  D_PTS(1,2)=D_XWC+XSIDE
                  D_PTS(2,2)=D_PTS(2,1)
                  D_PTS(1,3)=D_PTS(1,2)
                  D_PTS(2,3)=D_YWC+YSIDE
                  D_PTS(1,4)=D_PTS(1,1)
                  D_PTS(2,4)=D_PTS(2,3)
                  D_PTS(1,5)=D_PTS(1,1)
                  D_PTS(2,5)=D_PTS(2,1)
                  CALL GKS_PL(1,2,5,D_PTS,ERROR,*9999)
                  XSEGMENT_DATA(1,NTSG)=REAL(D_XWC)
                  YSEGMENT_DATA(1,NTSG)=REAL(D_YWC)
                  XSEGMENT_DATA(2,NTSG)=REAL(XSIDE)
                  YSEGMENT_DATA(2,NTSG)=REAL(YSIDE)
                ENDIF
              ENDDO
              CALL CLOSE_SEGMENT(ISDRAW(NTDRAW),iw,ERROR,*9999)
              CALL DAWK(iw,0,ERROR,*9999)

            ELSE IF(ABBREV(CHOOSE,'POLYLINE_BOX..',14)) THEN
              CALL GKS_STRG(56,NCHAR,'Enter X-coord of any corner [0]:',
     '          TEXT_STRING,ERROR,*9999)
              IF(NCHAR.GT.0) THEN
                D_XWC1=RFROMC(TEXT_STRING)
              ELSE
                D_XWC1=0.0D0
              ENDIF
              CALL GKS_STRG(56,NCHAR,
     '          'Enter Y-coord of that corner [0]:',
     '          TEXT_STRING,ERROR,*9999)
              IF(NCHAR.GT.0) THEN
                D_YWC1=RFROMC(TEXT_STRING)
              ELSE
                D_YWC1=0.0D0
              ENDIF
              CALL GKS_STRG(56,NCHAR,
     '          'Enter X-coord of diag. opp. corner [0]:',
     '          TEXT_STRING,ERROR,*9999)
              IF(NCHAR.GT.0) THEN
                D_XWC2=RFROMC(TEXT_STRING)
              ELSE
                D_XWC2=0.0D0
              ENDIF
              CALL GKS_STRG(56,NCHAR,
     '          'Enter Y-coord of that corner [0]:',
     '          TEXT_STRING,ERROR,*9999)
              IF(NCHAR.GT.0) THEN
                D_YWC2=RFROMC(TEXT_STRING)
              ELSE
                D_YWC2=0.0D0
              ENDIF
              CALL ACWK(iw,0,ERROR,*9999)
              NTDRAW=NTDRAW+1
              CALL OPEN_SEGMENT(ISDRAW(NTDRAW),ISEG,iw,
     '          'DRAW polyline_box',index,INDEX_OLD,0,1,
     '          CSEG,ERROR,*9999)
              CALL GKS_SPLI(LINE_INDEX,ERROR,*9999)
              CALL GKS_SPLR(iw,LINE_INDEX,LINE_TYPE,LINE_WIDTH,
     '          LINE_COLOUR_INDEX,ERROR,*9999)
              INSTAT=1
              D_XWC=0.5D0*(D_XWC1+D_XWC2)
              D_YWC=0.5D0*(D_YWC1+D_YWC2)
              XSIDE=DABS(D_XWC2-D_XWC1)/2.0D0
              YSIDE=DABS(D_YWC2-D_YWC1)/2.0D0
              D_PTS(1,1)=D_XWC-XSIDE
              D_PTS(2,1)=D_YWC-YSIDE
              D_PTS(1,2)=D_XWC+XSIDE
              D_PTS(2,2)=D_PTS(2,1)
              D_PTS(1,3)=D_PTS(1,2)
              D_PTS(2,3)=D_YWC+YSIDE
              D_PTS(1,4)=D_PTS(1,1)
              D_PTS(2,4)=D_PTS(2,3)
              D_PTS(1,5)=D_PTS(1,1)
              D_PTS(2,5)=D_PTS(2,1)
              CALL GKS_PL(1,2,5,D_PTS,ERROR,*9999)
              XSEGMENT_DATA(1,NTSG)=REAL(D_XWC)
              YSEGMENT_DATA(1,NTSG)=REAL(D_YWC)
              XSEGMENT_DATA(2,NTSG)=REAL(XSIDE)
              YSEGMENT_DATA(2,NTSG)=REAL(YSIDE)
              CALL CLOSE_SEGMENT(ISDRAW(NTDRAW),iw,ERROR,*9999)
              CALL DAWK(iw,0,ERROR,*9999)

            ELSE IF(ABBREV(CHOOSE,'POLYLINE_CIRCLE',15)) THEN
              CALL ACWK(iw,0,ERROR,*9999)
              NTDRAW=NTDRAW+1
              CALL OPEN_SEGMENT(ISDRAW(NTDRAW),ISEG,iw,
     '          'DRAW polyline_circle',index,INDEX_OLD,0,1,
     '          CSEG,ERROR,*9999)
              CALL GKS_SPLI(LINE_INDEX,ERROR,*9999)
              CALL GKS_SPLR(iw,LINE_INDEX,LINE_TYPE,LINE_WIDTH,
     '          LINE_COLOUR_INDEX,ERROR,*9999)
              INSTAT=1
              DO WHILE (INSTAT.EQ.1)
                WRITE(OP_STRING,
     '            '('' >>Locate centre of circle '//
     '            '(use 2nd button to end)'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                CALL LOCATOR(1,iw,INSTAT,'REQUEST',3,
     '            0.0D0,D_XWC1,0.0D0,D_YWC1,ERROR,*9999)
                IF(INSTAT.EQ.1) THEN
                  WRITE(OP_STRING,'('' >>Locate radius of circle'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  CALL LOCATOR(1,iw,INSTAT,'REQUEST',4,
     '              D_XWC1,D_XWC2,D_YWC1,D_YWC2,ERROR,*9999)
                  RADIUS=DSQRT((D_XWC2-D_XWC1)**2+(D_YWC2-D_YWC1)**2)
                  DO I=1,100
                    THETA=2.0D0*PI*DBLE(I-1)/99.0D0
                    D_PTS(1,I)=D_XWC1+RADIUS*DCOS(THETA)
                    D_PTS(2,I)=D_YWC1+RADIUS*DSIN(THETA)
                  ENDDO
                  CALL GKS_PL(1,2,100,D_PTS,ERROR,*9999)
                  XSEGMENT_DATA(1,NTSG)=REAL(D_XWC1)
                  YSEGMENT_DATA(1,NTSG)=REAL(D_YWC1)
                  XSEGMENT_DATA(2,NTSG)=REAL(RADIUS)
                  XSEGMENT_DATA(3,NTSG)=REAL(0.0D0)
                  YSEGMENT_DATA(3,NTSG)=REAL(2.0D0*PI)
                ENDIF
              ENDDO
              CALL CLOSE_SEGMENT(ISDRAW(NTDRAW),iw,ERROR,*9999)
              CALL DAWK(iw,0,ERROR,*9999)

            ELSE IF(ABBREV(CHOOSE,'POLYLINE_CIRCLE..',17)) THEN
              CALL GKS_STRG(56,NCHAR,'Enter X-coord of centre [0]:',
     '          TEXT_STRING,ERROR,*9999)
              IF(NCHAR.GT.0) THEN
                D_XWC1=RFROMC(TEXT_STRING)
              ELSE
                D_XWC1=0.0D0
              ENDIF
              CALL GKS_STRG(56,NCHAR,'Enter Y-coord of centre [0]:',
     '          TEXT_STRING,ERROR,*9999)
              IF(NCHAR.GT.0) THEN
                D_YWC1=RFROMC(TEXT_STRING)
              ELSE
                D_YWC1=0.0D0
              ENDIF
              CALL GKS_STRG(56,NCHAR,'Enter radius [0]:',
     '          TEXT_STRING,ERROR,*9999)
              IF(NCHAR.GT.0) THEN
                RADIUS=RFROMC(TEXT_STRING)
              ELSE
                RADIUS=0.0D0
              ENDIF
              DO i=1,100
                THETA=2.0D0*PI*DBLE(I-1)/99.0D0
                D_PTS(1,i)=D_XWC1+RADIUS*DCOS(THETA)
                D_PTS(2,i)=D_YWC1+RADIUS*DSIN(THETA)
              ENDDO
              CALL ACWK(iw,0,ERROR,*9999)
              NTDRAW=NTDRAW+1
              CALL OPEN_SEGMENT(ISDRAW(NTDRAW),ISEG,iw,
     '          'DRAW polyline_circle..',index,INDEX_OLD,0,1,
     '          CSEG,ERROR,*9999)
              CALL GKS_SPLI(LINE_INDEX,ERROR,*9999)
              CALL GKS_SPLR(iw,LINE_INDEX,LINE_TYPE,LINE_WIDTH,
     '          LINE_COLOUR_INDEX,ERROR,*9999)
              CALL GKS_PL(1,2,100,D_PTS,ERROR,*9999)
              XSEGMENT_DATA(1,NTSG)=REAL(D_XWC1)
              YSEGMENT_DATA(1,NTSG)=REAL(D_YWC1)
              XSEGMENT_DATA(2,NTSG)=REAL(RADIUS)
              XSEGMENT_DATA(3,NTSG)=0.0
              YSEGMENT_DATA(3,NTSG)=REAL(2.0D0*PI)
              CALL CLOSE_SEGMENT(ISDRAW(NTDRAW),iw,ERROR,*9999)
              CALL DAWK(iw,0,ERROR,*9999)

            ELSE IF(ABBREV(CHOOSE,'POLYLINE_SEMICIRCLE',19)) THEN
              CALL ACWK(iw,0,ERROR,*9999)
              NTDRAW=NTDRAW+1
              CALL OPEN_SEGMENT(ISDRAW(NTDRAW),ISEG,iw,
     '          'DRAW polyline_semicircle',index,INDEX_OLD,0,1,
     '          CSEG,ERROR,*9999)
              CALL GKS_SPLI(LINE_INDEX,ERROR,*9999)
              CALL GKS_SPLR(iw,LINE_INDEX,LINE_TYPE,LINE_WIDTH,
     '          LINE_COLOUR_INDEX,ERROR,*9999)
              INSTAT=1
              DO WHILE (INSTAT.EQ.1)
                WRITE(OP_STRING,
     '            '('' >>Locate centre of circle '//
     '            '(use 2nd button to end)'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                CALL LOCATOR(1,iw,INSTAT,'REQUEST',3,
     '            0.0D0,D_XWC1,0.0D0,D_YWC1,ERROR,*9999)
                IF(INSTAT.EQ.1) THEN
                  WRITE(OP_STRING,
     '              '('' >>Locate radius and start of semicircle'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  CALL LOCATOR(1,iw,INSTAT,'REQUEST',4,
     '              D_XWC1,D_XWC2,D_YWC1,D_YWC2,ERROR,*9999)
                  RADIUS=DSQRT((D_XWC2-D_XWC1)**2+(D_YWC2-D_YWC1)**2)
                  THETA1=DATAN2(D_YWC2-D_YWC1,D_XWC2-D_XWC1)
                  WRITE(OP_STRING,
     '              '('' >>Locate second angle of semicircle'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  CALL LOCATOR(1,iw,INSTAT,'REQUEST',4,
     '              D_XWC1,D_XWC2,D_YWC1,D_YWC2,ERROR,*9999)
                  THETA2=DATAN2(D_YWC2-D_YWC1,D_XWC2-D_XWC1)
                  ITOT=DINT(100.0D0*(THETA2-THETA1)/(2.0D0*PI))
                  DO i=1,ITOT
                    THETA=THETA1+DBLE(i-1)/DBLE(ITOT-1)*(THETA2-THETA1)
                    D_PTS(1,i)=D_XWC1+RADIUS*DCOS(THETA)
                    D_PTS(2,i)=D_YWC1+RADIUS*DSIN(THETA)
                  ENDDO
                  CALL GKS_PL(1,2,ITOT,D_PTS,ERROR,*9999)
                  XSEGMENT_DATA(1,NTSG)=REAL(D_XWC1)
                  YSEGMENT_DATA(1,NTSG)=REAL(D_YWC1)
                  XSEGMENT_DATA(2,NTSG)=REAL(RADIUS)
                  XSEGMENT_DATA(3,NTSG)=REAL(THETA1)
                  YSEGMENT_DATA(3,NTSG)=REAL(THETA2)
                ENDIF
              ENDDO
              CALL CLOSE_SEGMENT(ISDRAW(NTDRAW),iw,ERROR,*9999)
              CALL DAWK(iw,0,ERROR,*9999)

            ELSE IF(ABBREV(CHOOSE,'POLYLINE_SEMICIRCLE..',21)) THEN
              CALL GKS_STRG(56,NCHAR,'Enter X-coord of centre [0]:',
     '          TEXT_STRING,ERROR,*9999)
              IF(NCHAR.GT.0) THEN
                D_XWC1=RFROMC(TEXT_STRING)
              ELSE
                D_XWC1=0.0D0
              ENDIF
              CALL GKS_STRG(56,NCHAR,'Enter Y-coord of centre [0]:',
     '          TEXT_STRING,ERROR,*9999)
              IF(NCHAR.GT.0) THEN
                D_YWC1=RFROMC(TEXT_STRING)
              ELSE
                D_YWC1=0.0D0
              ENDIF
              CALL GKS_STRG(56,NCHAR,'Enter radius [0]:',
     '          TEXT_STRING,ERROR,*9999)
              IF(NCHAR.GT.0) THEN
                RADIUS=RFROMC(TEXT_STRING)
              ELSE
                RADIUS=0.0D0
              ENDIF
              CALL GKS_STRG(56,NCHAR,'Enter 1st angle [0deg]:',
     '          TEXT_STRING,ERROR,*9999)
              IF(NCHAR.GT.0) THEN
                THETA1=RFROMC(TEXT_STRING)
              ELSE
                THETA1=0.0D0
              ENDIF
              CALL GKS_STRG(56,NCHAR,'Enter 2nd angle [0deg]:',
     '          TEXT_STRING,ERROR,*9999)
              IF(NCHAR.GT.0) THEN
                THETA2=RFROMC(TEXT_STRING)
              ELSE
                THETA2=0.0D0
              ENDIF
              ITOT=DINT(100.0D0*(THETA2-THETA1)/360.0D0)
              DO i=1,ITOT
                THETA=THETA1+DBLE(I-1)/DBLE(ITOT-1)*(THETA2-THETA1)
                THETA=THETA*PI/180.0D0
                D_PTS(1,i)=D_XWC1+RADIUS*DCOS(THETA)
                D_PTS(2,i)=D_YWC1+RADIUS*DSIN(THETA)
              ENDDO
              CALL ACWK(iw,0,ERROR,*9999)
              NTDRAW=NTDRAW+1
              CALL OPEN_SEGMENT(ISDRAW(NTDRAW),ISEG,iw,
     '          'DRAW polyline_semicircle..',index,INDEX_OLD,0,1,
     '          CSEG,ERROR,*9999)
              CALL GKS_SPLI(LINE_INDEX,ERROR,*9999)
              CALL GKS_SPLR(iw,LINE_INDEX,LINE_TYPE,LINE_WIDTH,
     '          LINE_COLOUR_INDEX,ERROR,*9999)
              CALL GKS_PL(1,2,ITOT,D_PTS,ERROR,*9999)
              XSEGMENT_DATA(1,NTSG)=REAL(D_XWC1)
              YSEGMENT_DATA(1,NTSG)=REAL(D_YWC1)
              XSEGMENT_DATA(2,NTSG)=REAL(RADIUS)
              XSEGMENT_DATA(3,NTSG)=REAL(THETA1*PI/180.0D0)
              YSEGMENT_DATA(3,NTSG)=REAL(THETA2*PI/180.0D0)
              CALL CLOSE_SEGMENT(ISDRAW(NTDRAW),iw,ERROR,*9999)
              CALL DAWK(iw,0,ERROR,*9999)

            ELSE IF(ABBREV(CHOOSE,'POLYLINE_ELLIPSE',14)) THEN
              CALL ACWK(iw,0,ERROR,*9999)
              NTDRAW=NTDRAW+1
              CALL OPEN_SEGMENT(ISDRAW(NTDRAW),ISEG,iw,
     '          'DRAW polyline_ellipse',index,INDEX_OLD,0,1,
     '          CSEG,ERROR,*9999)
              CALL GKS_SPLI(LINE_INDEX,ERROR,*9999)
              CALL GKS_SPLR(iw,LINE_INDEX,LINE_TYPE,LINE_WIDTH,
     '          LINE_COLOUR_INDEX,ERROR,*9999)
              INSTAT=1
              DO WHILE (INSTAT.EQ.1)
                WRITE(OP_STRING,
     '            '('' >>Locate first corner of enclosing'//
     '            'box (use 2nd button to end)'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                CALL LOCATOR(1,iw,INSTAT,'REQUEST',3,
     '            0.0D0,D_XWC1,0.0D0,D_YWC1,ERROR,*9999)
                IF(INSTAT.EQ.1) THEN
                  WRITE(OP_STRING,
     '              '('' >>Locate second corner of enclosing'//
     '              ' box'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  CALL LOCATOR(1,iw,INSTAT,'REQUEST',5,
     '              D_XWC1,D_XWC2,D_YWC1,D_YWC2,ERROR,*9999)
                  X0=0.5D0*(D_XWC1+D_XWC2)
                  Y0=0.5D0*(D_YWC1+D_YWC2)
                  A=DABS(D_XWC2-D_XWC1)/2.0D0
                  B=DABS(D_YWC2-D_YWC1)/2.0D0
                  DO i=1,100
                    THETA=2.0D0*PI*DBLE(i-1)/99.0D0
                    D_PTS(1,i)=X0+A*DCOS(THETA)
                    D_PTS(2,i)=Y0+B*DSIN(THETA)
                  ENDDO
                  CALL GKS_PL(1,2,100,D_PTS,ERROR,*9999)
                  XSEGMENT_DATA(1,NTSG)=REAL(X0)
                  YSEGMENT_DATA(1,NTSG)=REAL(Y0)
                ENDIF
              ENDDO
              CALL CLOSE_SEGMENT(ISDRAW(NTDRAW),iw,ERROR,*9999)
              CALL DAWK(iw,0,ERROR,*9999)

            ELSE IF(ABBREV(CHOOSE,'POLYLINE_SEMIELLIPSE',20)) THEN

            ELSE IF(ABBREV(CHOOSE,'POLYLINE_SEMIELLIPSE..',22)) THEN

            ELSE IF(ABBREV(CHOOSE,'POLYMARKER_ATTRIBUTES',21)) THEN
              OPTION2(1)='Colour_index'
              OPTION2(2)='Marker_type'
              OPTION2(3)='Marker_size'
              OPTION2(4)='Set_index'
              OPTION2(5)='Pick_segment'
              OPTION2(6)='Set_current'
              CALL CHOICE('PM',1,1,INS,98,'EVENET',6,6,n1ch,NOCO,2,
     '          CO,OPTION2,STRING,XREFA,YREFA,ERROR,*9999)
              CONTINUEB=.TRUE.
              DO WHILE (CONTINUEB)
                CALL EVENT(ID_WS,ID_Device,Input_Status,CLASS,NOCH,
     '            R4DATA,SDATA,ERROR,*9999)
                IF(DOP) THEN
                  WRITE(OP_STRING,*)
     '              '        Input_Class=',Class
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                IF(CLASS(1:6).EQ.'CHOICE') THEN
                  IF(DOP) THEN
                     WRITE(OP_STRING,*)
     '                 '        Input_Choice=',NOCH
                     CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                  IF(NOCH.GT.0) THEN
                    IF(ABBREV(OPTION2(NOCH),'COLOUR_INDEX',8)) THEN
                      CALL COLOUR(MARKER_COLOUR_INDEX,iw,CO,STRING,
     '                  XREFA,YREFA,ERROR,*9999)
                    ELSE IF(ABBREV(OPTION2(NOCH),'MARKER_TYPE',8)) THEN
                      OPTION3(1)='Asterisk'
                      OPTION3(2)='Point'
                      OPTION3(3)='Plus_sign'
                      OPTION3(4)='Circle'
                      OPTION3(5)='Cross'
                      CALL CHOICE('MARKER_TYPE',1,1,INS,99,'REQUEST',
     '                  5,5,NOCH,NOCO,
     '                  2,CO,OPTION3,STRING,0.7*XDISP,YREFA,ERROR,*9999)
                      IF(ABBREV(OPTION3(NOCH),'ASTERISK',5)) THEN
                        MARKER_TYPE=GAST
                      ELSE IF(ABBREV(OPTION3(NOCH),'POINT'    ,5)) THEN
                        MARKER_TYPE=GPOINT
                      ELSE IF(ABBREV(OPTION3(NOCH),'PLUS_SIGN',5)) THEN
                        MARKER_TYPE=GPLUS
                      ELSE IF(ABBREV(OPTION3(NOCH),'CIRCLE'   ,5)) THEN
                        MARKER_TYPE=GOMARK
                      ELSE IF(ABBREV(OPTION3(NOCH),'CROSS'    ,5)) THEN
                        MARKER_TYPE=GXMARK
                      ENDIF
                    ELSE IF(ABBREV(OPTION2(NOCH),'MARKER_SIZE',8)) THEN
                      CALL VALUATOR('MARKER_SIZE',81,'REQUEST',2,0.,5.,
     '                  1.,MARKER_SIZE,XREFB,YREFB,ERROR,*9999)
                    ELSE IF(ABBREV(OPTION2(NOCH),'SET_INDEX',8)) THEN
                      OPTION3(1)='Index_1'
                      OPTION3(2)='Index_2'
                      OPTION3(3)='Index_3'
                      OPTION3(4)='Index_4'
                      CALL CHOICE('SET_INDEX',1,1,INS,99,'REQUEST',4,4,
     '                  NOCH,NOCO,2,
     '                  CO,OPTION3,STRING,0.7*XDISP,YREFA,ERROR,*9999)
                      IF(ABBREV(OPTION3(NOCH),'INDEX_1',7)) THEN
                        MARKER_INDEX=1
                      ELSE IF(ABBREV(OPTION3(NOCH),'INDEX_2',7)) THEN
                        MARKER_INDEX=2
                      ELSE IF(ABBREV(OPTION3(NOCH),'INDEX_3',7)) THEN
                        MARKER_INDEX=3
                      ELSE IF(ABBREV(OPTION3(NOCH),'INDEX_4',7)) THEN
                        MARKER_INDEX=4
                      ENDIF
                    ELSE IF(ABBREV(OPTION2(NOCH),'PICK_SEGMENT',8)) THEN
                      DO nosg=1,NTSG
                        IF(ISEG(nosg).GT.0.AND.CSEG(nosg)(5:8).EQ.
     '                    'mark')
     '                    CALL DETECT(iw,ISEG,nosg,'DETECTABLE',
     '                      ERROR,*9999)
                      ENDDO
                      WRITE(OP_STRING,
     '                  '('' >>Pick segment to change '//
     '                  '(use 2nd button to end)'')')
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      CALL PICK(iw,LD1,'REQUEST',INSTAT,isegm,IPICK,
     '                  ERROR,*9999)
                      DO WHILE (INSTAT.EQ.1  )
                        WRITE(OP_STRING,
     '                    '('' Segment number='',I3,'
     '                    //''' Type='',A10)') isegm,CSEG(isegm)
                        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                        CALL GKS_SPMR(iw,isegm,MARKER_TYPE,MARKER_SIZE,
     '                    MARKER_COLOUR_INDEX,ERROR,*9999)
                        CALL GKS_UWK(iw,GPERFO,ERROR,*9999)
                        CALL PICK(iw,LD1,'REQUEST',INSTAT,isegm,IPICK,
     '                    ERROR,*9999)
                      ENDDO
                      DO nosg=1,NTSG
                        IF(ISEG(nosg).GT.0.AND.CSEG(nosg)(5:8).EQ.
     '                    'mark')
     '                    CALL DETECT(iw,ISEG,nosg,'UNDETECTABLE',
     '                      ERROR,*9999)
                      ENDDO
                      CONTINUEB=.FALSE.
                    ELSE IF(ABBREV(OPTION2(NOCH),'SET_CURRENT',8)) THEN
                      CALL GKS_SPMI(MARKER_INDEX,ERROR,*9999)
                      CALL GKS_SPMR(iw,MARKER_INDEX,MARKER_TYPE,
     '                  MARKER_SIZE,MARKER_COLOUR_INDEX,ERROR,*9999)
                      CONTINUEB=.FALSE.
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO
              CALL INPUT_MODE(98,1,'CHOICE','REQUEST',ERROR,*9999)

            ELSE IF(ABBREV(CHOOSE,'POLYMARKER_LOCATE',17)) THEN
              CALL ACWK(iw,0,ERROR,*9999)
              NTDRAW=NTDRAW+1
              CALL OPEN_SEGMENT(ISDRAW(NTDRAW),ISEG,iw,
     '          'DRAW polymarker_locate',index,INDEX_OLD,0,1,
     '          CSEG,ERROR,*9999)
              CALL GKS_SPMI(MARKER_INDEX,ERROR,*9999)
              CALL GKS_SPMR(iw,MARKER_INDEX,MARKER_TYPE,MARKER_SIZE,
     '          MARKER_COLOUR_INDEX,ERROR,*9999)
              WRITE(OP_STRING,'('' >>Locate polymarker position '//
     '          '(use 2nd button to end)'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              CALL LOCATOR(1,iw,INSTAT,'REQUEST',3,
     '          0.0D0,D_PTS(1,1),0.0D0,D_PTS(2,1),ERROR,*9999)
              DO WHILE (INSTAT.EQ.1)
                CALL GKS_PM(1,2,1,D_PTS,ERROR,*9999)
                XSEGMENT_DATA(1,NTSG)=REAL(D_PTS(1,1))
                YSEGMENT_DATA(1,NTSG)=REAL(D_PTS(2,1))
                CALL LOCATOR(1,iw,INSTAT,'REQUEST',3,
     '            0.0D0,D_PTS(1,1),0.0D0,D_PTS(2,1),ERROR,*9999)
              ENDDO
              CALL CLOSE_SEGMENT(ISDRAW(NTDRAW),iw,ERROR,*9999)
              CALL DAWK(iw,0,ERROR,*9999)

            ELSE IF(ABBREV(CHOOSE,'TEXT_ATTRIBUTES',15)) THEN
              OPTION2( 1)='Alignment'
              OPTION2( 2)='Colour_index'
              OPTION2( 3)='Expansion_factor'
              OPTION2( 4)='Font'
              OPTION2( 5)='Height'
              OPTION2( 6)='Path'
              OPTION2( 7)='Precision'
              OPTION2( 8)='Spacing'
              OPTION2( 9)='Upvector'
              OPTION2(10)='Set_index'
              OPTION2(11)='Pick_segment'
              OPTION2(12)='Set_current'
              CALL CHOICE('TEXT',1,1,INS,98,'EVENT',12,12,n1ch,NOCO,2,
     '          CO,OPTION2,STRING,XREFA,YREFA,ERROR,*9999)
              CONTINUEB=.TRUE.
              DO WHILE (CONTINUEB)
                CALL EVENT(ID_WS,ID_Device,Input_Status,CLASS,NOCH,
     '            R4DATA,SDATA,ERROR,*9999)
                IF(DOP) THEN
                  WRITE(OP_STRING,*) '        Input_Class=',Class
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                IF(CLASS(1:6).EQ.'CHOICE') THEN
                  IF(DOP) THEN
                    WRITE(OP_STRING,*)
     '                '        Input_Choice=',NOCH
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                  IF(NOCH.GT.0) THEN
                    IF(ABBREV(OPTION2(NOCH),'ALIGNMENT',4)) THEN
                      OPTION3(1)='Horizontal'
                      OPTION3(2)='Vertical'
                      CALL CHOICE('ALIGNMENT',1,1,INS,99,'REQUEST',2,2,
     '                  NOCH,NOCO,2,
     '                  CO,OPTION3,STRING,0.7*XDISP,YREFA,ERROR,*9999)
                      IF(     ABBREV(OPTION3(NOCH),'HORIZONTAL',2)) THEN
                        OPTION3(1)='Normal'
                        OPTION3(2)='Left'
                        OPTION3(3)='Centre'
                        OPTION3(4)='Right'
                        CALL CHOICE('HOROZONTAL',1,1,INS,99,'REQUEST',
     '                    4,4,NOCH,NOCO,2,CO,OPTION3,STRING,0.7*XDISP,
     '                    YREFA,ERROR,*9999)
                        IF(     ABBREV(OPTION3(NOCH),'NORMAL',2)) THEN
                          TEXT_HORIZ_ALIGN=GAHNOR
                        ELSE IF(ABBREV(OPTION3(NOCH),'LEFT'  ,2)) THEN
                          TEXT_HORIZ_ALIGN=GALEFT
                        ELSE IF(ABBREV(OPTION3(NOCH),'CENTRE',2)) THEN
                          TEXT_HORIZ_ALIGN=GACENT
                        ELSE IF(ABBREV(OPTION3(NOCH),'RIGHT' ,2)) THEN
                          TEXT_HORIZ_ALIGN=GARITE
                        ENDIF
                      ELSE IF(ABBREV(OPTION3(NOCH),'VERTICAL' , 2)) THEN
                        OPTION3(1)='Normal'
                        OPTION3(2)='Top'
                        OPTION3(3)='Cap'
                        OPTION3(4)='Half'
                        OPTION3(5)='Base'
                        OPTION3(6)='Bottom'
                        CALL CHOICE('VERTICAL',1,1,INS,99,'REQUEST',6,
     '                    6,NOCH,NOCO,2,
     '                    CO,OPTION3,STRING,0.7*XDISP,YREFA,ERROR,*9999)
                        IF(     ABBREV(OPTION3(NOCH),'NORMAL',2)) THEN
                          TEXT_VERT_ALIGN=GAVNOR
                        ELSE IF(ABBREV(OPTION3(NOCH),'TOP'   ,2)) THEN
                          TEXT_VERT_ALIGN=GATOP
                        ELSE IF(ABBREV(OPTION3(NOCH),'CAP'   ,2)) THEN
                          TEXT_VERT_ALIGN=GACAP
                        ELSE IF(ABBREV(OPTION3(NOCH),'HALF'  ,2)) THEN
                          TEXT_VERT_ALIGN=GAHALF
                        ELSE IF(ABBREV(OPTION3(NOCH),'BASE'  ,2)) THEN
                          TEXT_VERT_ALIGN=GABASE
                        ELSE IF(ABBREV(OPTION3(NOCH),'BOTTOM',2)) THEN
                          TEXT_VERT_ALIGN=GABOTT
                        ENDIF
                      ENDIF
                      CALL GKS_STXAL(TEXT_HORIZ_ALIGN,TEXT_VERT_ALIGN,
     '                  ERROR,*9999)
                    ELSE IF(ABBREV(OPTION2(NOCH),'COLOUR_INDEX',4)) THEN
                      CALL COLOUR(TEXT_COLOUR_INDEX,iw,CO,STRING,XREFA,
     '                  YREFA,ERROR,*9999)
                    ELSE IF(ABBREV(OPTION2(NOCH),'EXPANSION_FACTOR',4))
     '                THEN
                      CALL VALUATOR('EXPANSION',81,'REQUEST',2,0.,5.,1.,
     '                  TEXT_EXPFAC,XREFB,YREFB,ERROR,*9999)
                    ELSE IF(ABBREV(OPTION2(NOCH),'FONT',4)) THEN
                      OPTION3( 1)='+1'
                      OPTION3( 2)='-1'
                      OPTION3( 3)='-2'
                      OPTION3( 4)='-3'
                      OPTION3( 5)='-4'
                      OPTION3( 6)='-5'
                      OPTION3( 7)='-6'
                      OPTION3( 8)='-7'
                      OPTION3( 9)='-8'
                      OPTION3(10)='-9'
                      OPTION3(11)='-10'
                      OPTION3(12)='-11'
                      OPTION3(13)='-12'
                      OPTION3(14)='-13'
                      OPTION3(15)='-14'
                      OPTION3(16)='-15'
                      OPTION3(17)='-16'
                      OPTION3(18)='-17'
                      OPTION3(19)='-18'
                      OPTION3(20)='-19'
                      OPTION3(21)='-20'
                      OPTION3(22)='-21'
                      OPTION3(23)='-22'
                      OPTION3(24)='-23'
                      CALL CHOICE('FONT',1,1,INS,99,'REQUEST',24,24,
     '                  NOCH,NOCO,2,CO,OPTION3,STRING,0.7*XDISP,YREFA,
     '                  ERROR,*9999)
                      IF(NOCH.EQ.1) THEN
                        TEXT_FONT=1
                      ELSE
                        TEXT_FONT=1-NOCH
                      ENDIF
                    ELSE IF(ABBREV(OPTION2(NOCH),'HEIGHT',4)) THEN
                      CALL VALUATOR('TEXT_HEIGHT',81,'REQUEST',2,0.,5.,
     '                  1.,TEXT_HEIGHT,XREFB,YREFB,ERROR,*9999)
                    ELSE IF(ABBREV(OPTION2(NOCH),'PATH',4)) THEN
                      OPTION3(1)='Right'
                      OPTION3(2)='Left'
                      OPTION3(3)='Up'
                      OPTION3(4)='Down'
                      CALL CHOICE('PATH',1,1,INS,99,'REQUEST',4,4,NOCH,
     '                  NOCO,2,CO,OPTION3,STRING,0.7*XDISP,YREFA,ERROR,
     '                  *9999)
                      IF(     ABBREV(OPTION3(NOCH),'RIGHT',2)) THEN
                        TEXT_PATH=GRIGHT
                      ELSE IF(ABBREV(OPTION3(NOCH),'LEFT', 2)) THEN
                        TEXT_PATH=GLEFT
                      ELSE IF(ABBREV(OPTION3(NOCH),'UP'  , 2)) THEN
                        TEXT_PATH=GDOWN
                      ELSE IF(ABBREV(OPTION3(NOCH),'DOWN', 2)) THEN
                        TEXT_PATH=GDOWN
                      ENDIF
                    ELSE IF(ABBREV(OPTION2(NOCH),'PRECISION',4)) THEN
                      OPTION3(1)='String    (fonts 1,-1.. -3)'
                      OPTION3(2)='Character (fonts 1,-1.. -3)'
                      OPTION3(3)='Stroke    (fonts 1,-2..-23)'
                      CALL CHOICE('PRECISION',1,1,INS,99,'REQUEST',3,3,
     '                  NOCH,NOCO,2,
     '                  CO,OPTION3,STRING,0.7*XDISP,YREFA,ERROR,*9999)
                      IF(     OPTION3(NOCH)(1:6).eq.'String') THEN
                        TEXT_PRECISION=GSTRP
                      ELSE IF(OPTION3(NOCH)(1:9).eq.'Character') THEN
                        TEXT_PRECISION=GCHARP
                      ELSE IF(OPTION3(NOCH)(1:6).eq.'Stroke') THEN
                        TEXT_PRECISION=GSTRKP
                      ENDIF
                    ELSE IF(ABBREV(OPTION2(NOCH),'SPACING',4)) THEN
                      CALL VALUATOR('TEXT_SPACING',81,'REQUEST',2,0.,
     '                  5.,1.,TEXT_SPACING,XREFB,YREFB,ERROR,*9999)
                    ELSE IF(ABBREV(OPTION2(NOCH),'UPVECTOR',4)) THEN
                    ELSE IF(ABBREV(OPTION2(NOCH),'SET_INDEX',8)) THEN
                      OPTION3(1)='Index_1'
                      OPTION3(2)='Index_2'
                      OPTION3(3)='Index_3'
                      OPTION3(4)='Index_4'
                      CALL CHOICE('SET_INDEX',1,1,INS,99,'REQUEST',4,4,
     '                  NOCH,NOCO,2,
     '                  CO,OPTION3,STRING,0.7*XDISP,YREFA,ERROR,*9999)
                      IF(ABBREV(OPTION3(NOCH),'INDEX_1',7)) THEN
                        TEXT_INDEX=1
                      ELSE IF(ABBREV(OPTION3(NOCH),'INDEX_2',7)) THEN
                        TEXT_INDEX=2
                      ELSE IF(ABBREV(OPTION3(NOCH),'INDEX_3',7)) THEN
                        TEXT_INDEX=3
                      ELSE IF(ABBREV(OPTION3(NOCH),'INDEX_4',7)) THEN
                        TEXT_INDEX=4
                      ENDIF
                    ELSE IF(ABBREV(OPTION2(NOCH),'PICK_SEGMENT',4)) THEN
                      DO nosg=1,NTSG
                        IF(ISEG(nosg).GT.0.AND.CSEG(nosg)(5:8).EQ.
     '                    'text')
     '                    CALL DETECT(iw,ISEG,nosg,'DETECTABLE',ERROR,
     '                      *9999)
                      ENDDO
                      WRITE(OP_STRING,'('' >>Pick segment '//
     '                  'to change (use 2nd button to end)'')')
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      CALL PICK(iw,LD1,'REQUEST',INSTAT,isegm,IPICK,
     '                  ERROR,*9999)
                      DO WHILE (INSTAT.EQ.1  )
                        WRITE(OP_STRING,'('' Segment number='',I3,'
     '                    //''' Type='',A10)') isegm,CSEG(isegm)
                        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                        CALL GKS_STXR(iw,NTSG,TEXT_FONT,TEXT_PRECISION,
     '                    TEXT_EXPFAC,TEXT_SPACING,TEXT_COLOUR_INDEX,
     '                    ERROR,*9999)
                        CALL GKS_SCHH(TEXT_HEIGHT,ERROR,*9999)
                        CALL GKS_UWK(iw,GPERFO,ERROR,*9999)
                        CALL PICK(iw,LD1,'REQUEST',INSTAT,isegm,IPICK,
     '                    ERROR,*9999)
                      ENDDO
                      DO nosg=1,NTSG
                        IF(ISEG(nosg).GT.0.AND.CSEG(nosg)(5:8).EQ.
     '                    'text')
     '                    CALL DETECT(iw,ISEG,nosg,'UNDETECTABLE',
     '                      ERROR,*9999)
                      ENDDO
                      CONTINUEB=.FALSE.
                    ELSE IF(ABBREV(OPTION2(NOCH),'SET_CURRENT',4)) THEN
                      CALL GKS_STXI(TEXT_INDEX,ERROR,*9999)
                      CALL GKS_STXR(iw,TEXT_INDEX,TEXT_FONT,
     '                  TEXT_PRECISION,TEXT_EXPFAC,TEXT_SPACING,
     '                  TEXT_COLOUR_INDEX,ERROR,*9999)
                      CALL GKS_SCHH(TEXT_HEIGHT,ERROR,*9999)
                      CONTINUEB=.FALSE.
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO
              CALL INPUT_MODE(98,1,'CHOICE','REQUEST',ERROR,*9999)

            ELSE IF(ABBREV(CHOOSE,'TEXT_LOCATE',11)) THEN
C             CALL GKS_STXI(TEXT_INDEX,ERROR,*9999)
C             CALL GKS_STXR(iw,TEXT_INDEX,TEXT_FONT,TEXT_PRECISION,
C    '          TEXT_EXPFAC,TEXT_SPACING,TEXT_COLOUR_INDEX,ERROR,*9999)
C             CALL GKS_SCHH(TEXT_HEIGHT,ERROR,*9999)
C             CALL GKS_STXP(TEXT_PATH,ERROR,*9999)
C             CALL GKS_STXAL(TEXT_HORIZ_ALIGN,TEXT_VERT_ALIGN,
C    '          ERROR,*9999)
C             TEXT_STRING(1:80)=' '
C             WRITE(IOOP,'($,'' >>Enter text:'')')
C             READ(IOIP,*) TEXT_STRING
              CALL ACWK(iw,0,ERROR,*9999)
              NTDRAW=NTDRAW+1
              CALL OPEN_SEGMENT(ISDRAW(NTDRAW),ISEG,iw,
     '          'DRAW text_locate',index,INDEX_OLD,0,1,
     '          CSEG,ERROR,*9999)
              INSTAT=GOK
              DO WHILE (INSTAT.EQ.GOK)
                WRITE(OP_STRING,'('' >>Locate text start point '//
     '            '(use 2nd button to end)'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                CALL GKSTEXT(iw,INSTAT,'>',TEXT_STRING,D_XWC,D_YWC,
     '            ERROR,*9999)
                IF(INSTAT.EQ.GOK) THEN
                  CALL STRING_TRIM(TEXT_STRING,IBEG,IEND)
                  CALL GKS_STXAL(GALEFT,GABASE,ERROR,*9999)
                  CALL GKS_TX(D_XWC,D_YWC,TEXT_STRING(IBEG:IEND),
     '              ERROR,*9999)
                  XSEGMENT_DATA(1,NTSG)=D_XWC
                  YSEGMENT_DATA(1,NTSG)=D_YWC
                ENDIF
              ENDDO
              CALL CLOSE_SEGMENT(ISDRAW(NTDRAW),iw,ERROR,*9999)
              CALL DAWK(iw,0,ERROR,*9999)

            ELSE IF(ABBREV(CHOOSE,'SEGMENT_REMOVE',14)) THEN
              CALL ACWK(iw,0,ERROR,*9999)
              WRITE(OP_STRING,'('' >>Pick segment to remove '//
     '          '(use 2nd button to end)'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              CALL PICK(iw,LD1,'REQUEST',INSTAT,isegm,IPICK,
     '          ERROR,*9999)
              DO WHILE (INSTAT.EQ.1)
                  WRITE(OP_STRING,'('' Segment number='',I3,'
     '              //''' Type='',A10)') isegm,CSEG(isegm)
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  IF(ISEG(isegm).GT.0) THEN
                    CALL GKS_DSG(isegm,ERROR,*9999)
                    ISEG(isegm)=0
                  ENDIF
                  CALL GKS_UWK(iw,GPERFO,ERROR,*9999)
                  CALL PICK(iw,LD1,'REQUEST',INSTAT,isegm,IPICK,
     '              ERROR,*9999)
              ENDDO
              CALL DAWK(iw,0,ERROR,*9999)

            ELSE IF(ABBREV(CHOOSE,'SEGMENT_REPOSITION',18)) THEN
              CALL ACWK(iw,0,ERROR,*9999)
              WRITE(OP_STRING,'('' >>Pick segment to reposition '//
     '          '(use 2nd button to end)'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              CALL PICK(iw,LD1,'REQUEST',INSTAT,isegm,IPICK,
     '          ERROR,*9999)
              DO WHILE (INSTAT.EQ.1)
                WRITE(OP_STRING,'('' Segment number='',I3,'
     '            //''' Type='',A10)') isegm,CSEG(isegm)
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                IF(ISEG(isegm).GT.0) THEN
                  WRITE(OP_STRING,'('' >>Locate new position'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  CALL LOCATOR(1,iw,INSTAT,'REQUEST',3,0.0D0,D_XWC,
     '              0.0D0,D_YWC,ERROR,*9999)
                  XVECT=REAL(D_XWC)-XSEGMENT_DATA(1,isegm)
                  YVECT=REAL(D_YWC)-YSEGMENT_DATA(1,isegm)
                  IF(.NOT.LTRANS(isegm)) THEN
                    LTRANS(isegm)=.TRUE.
                    CALL GKS_EVTM(0.0,0.0,XVECT,YVECT,0.0,1.0,1.0,GWC,
     '                TRANSF(1,1,isegm),ERROR,*9999)
                  ELSE IF(LTRANS(isegm)) THEN
                    CALL GKS_ACTM(TRANSF(1,1,isegm),0.0,0.0,XVECT,YVECT,
     '                0.0,1.0,1.0,GWC,TRANSF(1,1,isegm),ERROR,*9999)
                  ENDIF
                  CALL GKS_SSGT(isegm,TRANSF(1,1,isegm),ERROR,*9999)
                  CALL GKS_UWK(iw,GPERFO,ERROR,*9999)
                  XSEGMENT_DATA(1,isegm)=REAL(D_XWC)
                  YSEGMENT_DATA(1,isegm)=REAL(D_YWC)
                ENDIF
                CALL PICK(iw,LD1,'REQUEST',INSTAT,isegm,IPICK,
     '            ERROR,*9999)
              ENDDO
              CALL DAWK(iw,0,ERROR,*9999)

            ELSE IF(ABBREV(CHOOSE,'SEGMENT_RESIZE',14)) THEN
              CALL ACWK(iw,0,ERROR,*9999)
              WRITE(OP_STRING,'('' >>Pick segment to resize '//
     '          '(use 2nd button to end)'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              CALL PICK(iw,LD1,'REQUEST',INSTAT,isegm,IPICK,
     '          ERROR,*9999)
              DO WHILE (INSTAT.EQ.1)
                WRITE(OP_STRING,'('' Segment number='',I3,'
     '            //''' Type='',A10)') isegm,CSEG(isegm)
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                IF(ISEG(isegm).GT.0) THEN
                  CALL VALUATOR('SEG_RESIZE',81,'REQUEST',3,0.,5.,1.,
     '              VALUE,XSEGMENT_DATA(1,isegm),YSEGMENT_DATA(1,isegm),
     '              ERROR,*9999)
                  WRITE(OP_STRING,'('' Value='',E11.4)') VALUE
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  IF(.NOT.LTRANS(isegm)) THEN
                    LTRANS(isegm)=.TRUE.
                    CALL GKS_EVTM(XSEGMENT_DATA(1,isegm),
     '                YSEGMENT_DATA(1,isegm),0.,0.,0.,VALUE,VALUE,GWC,
     '                TRANSF(1,1,isegm),ERROR,*9999)
                  ELSE IF(LTRANS(isegm)) THEN
                    CALL GKS_ACTM(TRANSF(1,1,isegm),
     '                XSEGMENT_DATA(1,isegm),YSEGMENT_DATA(1,isegm),0.,
     '                0.,0.,VALUE,VALUE,GWC,TRANSF(1,1,isegm),
     '                ERROR,*9999)
                  ENDIF
                  CALL GKS_SSGT(isegm,TRANSF(1,1,isegm),ERROR,*9999)
                  CALL GKS_UWK(iw,GPERFO,ERROR,*9999)
                ENDIF
                CALL PICK(iw,LD1,'REQUEST',INSTAT,isegm,IPICK,
     '            ERROR,*9999)
              ENDDO
              CALL DAWK(iw,0,ERROR,*9999)

            ELSE IF(ABBREV(CHOOSE,'SEGMENT_ROTATE',14)) THEN
              CALL ACWK(iw,0,ERROR,*9999)
              WRITE(OP_STRING,'('' >>Pick segment to rotate '//
     '          '(use 2nd button to end)'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              CALL PICK(iw,LD1,'REQUEST',INSTAT,isegm,IPICK,
     '          ERROR,*9999)
              DO WHILE (INSTAT.EQ.1)
                WRITE(OP_STRING,'('' Segment number='',I3,'
     '            //''' Type='',A10)') isegm,CSEG(isegm)
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                IF(ISEG(isegm).GT.0) THEN
                  CALL VALUATOR('SEG_ROTATE',81,'REQUEST',3,0.,5.,1.,
     '              VALUE,XSEGMENT_DATA(1,isegm),YSEGMENT_DATA(1,isegm),
     '              ERROR,*9999)
                  WRITE(OP_STRING,'('' Value='',E11.4)') VALUE
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  IF(.NOT.LTRANS(isegm)) THEN
                    LTRANS(isegm)=.TRUE.
                    CALL GKS_EVTM(XSEGMENT_DATA(1,isegm),
     '                YSEGMENT_DATA(1,isegm),0.,0.,VALUE,1.,1.,GWC,
     '                TRANSF(1,1,isegm),ERROR,*9999)
                  ELSE IF(LTRANS(isegm)) THEN
                    CALL GKS_ACTM(TRANSF(1,1,isegm),
     '                XSEGMENT_DATA(1,isegm),YSEGMENT_DATA(1,isegm),0.,
     '                0.,VALUE,1.,1.,GWC,TRANSF(1,1,isegm),ERROR,*9999)
                  ENDIF
                  CALL GKS_SSGT(isegm,TRANSF(1,1,isegm),ERROR,*9999)
                  CALL GKS_UWK(iw,GPERFO,ERROR,*9999)
                ENDIF
                CALL PICK(iw,LD1,'REQUEST',INSTAT,isegm,IPICK,
     '            ERROR,*9999)
              ENDDO
              CALL DAWK(iw,0,ERROR,*9999)

            ELSE IF(ABBREV(CHOOSE,'SET_WC_VIEWPORT',15)) THEN

            ELSE IF(ABBREV(CHOOSE,'SET_WC_WINDOW',13)) THEN

            ELSE IF(ABBREV(CHOOSE,'SET_WS_VIEWPORT',15)) THEN

            ELSE IF(ABBREV(CHOOSE,'SET_WS_WINDOW',13)) THEN

            ELSE IF(ABBREV(CHOOSE,'INQUIRE',7)) THEN

            ELSE IF(ABBREV(CHOOSE,'CLEAR_WS',8)) THEN
              CALL GKS_CLRWK(iw,GALWAY,ERROR,*9999)

            ELSE IF(ABBREV(CHOOSE,'EXIT',4)) THEN
              CONTINUEA=.FALSE.
            ENDIF

c            IF(   ABBREV(CHOOSE,'CELL_ARRAY',            10)
c     '        .OR.ABBREV(CHOOSE,'FILL_AREA_BOX',         13)
c     '        .OR.ABBREV(CHOOSE,'FILL_AREA_CIRCLE',      16)
c     '        .OR.ABBREV(CHOOSE,'FILL_AREA_ELLIPSE',     17)
c     '        .OR.ABBREV(CHOOSE,'POLYLINE_LOCATE',       15)
c     '        .OR.ABBREV(CHOOSE,'POLYLINE_STROKE',       15)
c     '        .OR.ABBREV(CHOOSE,'POLYLINE_ARROW',        14)
c     '        .OR.ABBREV(CHOOSE,'POLYLINE_BOX',          12)
c     '        .OR.ABBREV(CHOOSE,'POLYLINE_BOX..',        14)
c     '        .OR.ABBREV(CHOOSE,'POLYLINE_CIRCLE',       15)
c     '        .OR.ABBREV(CHOOSE,'POLYLINE_SEMICIRCLE',   19)
c     '        .OR.ABBREV(CHOOSE,'POLYLINE_SEMICIRCLE..', 21)
c     '        .OR.ABBREV(CHOOSE,'POLYLINE_ELLIPSE',      16)
c     '        .OR.ABBREV(CHOOSE,'POLYLINE_SEMIELLIPSE',  20)
c     '        .OR.ABBREV(CHOOSE,'POLYLINE_SEMIELLIPSE..',22)
c     '        .OR.ABBREV(CHOOSE,'POLYMARKER_LOCATE',     17)
c     '        .OR.ABBREV(CHOOSE,'POLYLINE_SEMIELLIPSE',  20)
c     '        .OR.ABBREV(CHOOSE,'TEXT_LOCATE',           11)) THEN
c              CALL GKS_CLSG(ERROR,*9999)
            IF(CHOOSE(1:7).EQ.'Segment') THEN
              DO nosg=1,NTSG
                IF(ISEG(nosg).GT.0)
     '            CALL DETECT(iw,ISEG,nosg,'UNDETECTABLE',ERROR,*9999)
              ENDDO
              WRITE(OP_STRING,
     '          '('' >>Pick segment on '',I2,'//
     '          ''' (use 2nd button to end)'')') iw
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            WRITE(OP_STRING,'('' >>Choose from menu'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

          ENDIF
        ELSE IF(CLASS(1:6).EQ.'CHOICE') THEN   !Choice from another menu
          CONTINUEA=.FALSE.
        ENDIF
      ENDDO

      CALL INPUT_MODE(97,1,'CHOICE','REQUEST',ERROR,*9999)

      CALL EXITS('GKS_DRAW')
      RETURN
 9999 CALL ERRORS('GKS_DRAW',ERROR)
      CALL EXITS('GKS_DRAW')
      RETURN 1
      END


      SUBROUTINE GKS_STRG(iw,NCHAR,PROMPT_STRING,TEXT_STRING,ERROR,*)

C#### Subroutine: GKS_STRG
C###  Description:
C**** Prompts the user for text under GKS.
C**** Note: Workstation iw is activated and deactivated in this routine.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:echo00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
!     Parameter List
      INTEGER iw,NCHAR
      CHARACTER ERROR*(*),PROMPT_STRING*(*),TEXT_STRING*(*)
!     Local Variables
      INTEGER IBEG,IEND,INSTAT,NCHAR_PROMPT
      REAL XDC1,XDC2,YDC1,YDC2

      CALL ENTERS('GKS_STRG',*9999)
      CALL SETUP(iw,ERROR,*9999)
      TEXT_STRING=' '
      CALL STRING_TRIM(PROMPT_STRING,IBEG,IEND)
      NCHAR_PROMPT=IEND-IBEG+1
      IF(DOP) THEN
        WRITE(OP_STRING,
     '    '('' no characters in prompt string = '',I3)')
     '    NCHAR_PROMPT
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      IF(iw.EQ.55) THEN      !window on RHS of screen
        XDC1=0.50*XDISP
        XDC2=1.00*XDISP
        YDC1=0.60*YDISP
        YDC2=0.75*YDISP
      ELSE IF(iw.EQ.56) THEN !window in centre of screen
        XDC1=0.15*XDISP
        XDC2=0.85*XDISP
        YDC1=0.40*YDISP
        YDC2=0.60*YDISP
      ENDIF
      ECAREA(1)=XDC1
      ECAREA(2)=XDC2
      ECAREA(3)=YDC1
      ECAREA(4)=YDC2
      CALL GKS_INST(iw,1,NCHAR_PROMPT,PROMPT_STRING(IBEG:IEND),1,
     '  80,NCHAR_PROMPT+1,0,BLANK,ERROR,*9999)
      CALL GKS_RQST(iw,1,INSTAT,NCHAR,TEXT_STRING,ERROR,*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,
     '    '('' no characters in string = '',I3)') NCHAR
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      IF(NCHAR.EQ.NCHAR_PROMPT) THEN
        TEXT_STRING=' '
        NCHAR=0
      ELSE
        TEXT_STRING=TEXT_STRING(NCHAR_PROMPT+1:)
      ENDIF

      CALL EXITS('GKS_STRG')
      RETURN
 9999 CALL ERRORS('GKS_STRG',ERROR)
      CALL EXITS('GKS_STRG')
      RETURN 1
      END


      SUBROUTINE GKSTEXT(iw,INSTAT,PROMPT_STRING,TEXT_STRING,XWC,YWC,
     '  ERROR,*)

C#### Subroutine: GKSTEXT
C###  Description:
C**** Prompts the user for text under GKS.
C**** XDC1,YDC1 are device coords of bottom left of iw viewport
C**** XDC2,YDC2 are device coords of top right   of iw viewport
C**** X_MIN,Y_MIN are world  coords of bottom left of iw viewport
C**** X_MAX,Y_MAX are world  coords of top right   of iw viewport
C**** XWC,YWC are returned world coordinates of reference point
C**** INSTAT is 1 if successful 0 if not

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:echo00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
!     Parameter List
      INTEGER INSTAT,iw
      REAL XWC,YWC
      CHARACTER ERROR*(*),PROMPT_STRING*(*),TEXT_STRING*(*)
!     Local Variables
      INTEGER IBEG,IEND,INSTAT1,IWCNDC,LD1,LEN_TEXT,NCHAR
      REAL XDC1,XDC2,X_MAX,X_MIN,
     '     YDC1,YDC2,Y_MAX,Y_MIN
      DATA LD1/1/

      CALL ENTERS('GKSTEXT',*9999)
      LEN_TEXT=LEN(TEXT_STRING)
c     write(*,'('' no bytes in text string='',I4)') len_text
      TEXT_STRING=' '
      CALL STRING_TRIM(PROMPT_STRING,IBEG,IEND)
      IF(iw.EQ.1) THEN !1st window in FEM environment
        IF(NJT.EQ.2) THEN
          XDC1=0.0
          XDC2=0.49*XDISP
          YDC1=YDISP-0.49*XDISP
          YDC2=YDISP
          X_MIN=XMIN
          X_MAX=XMAX
          Y_MIN=YMIN
          Y_MAX=YMAX
        ELSE IF(NJT.EQ.3) THEN
          XDC1=0.0
          XDC2=0.48*DISP
          YDC1=0.51*DISP
          YDC2=0.99*DISP
          X_MIN=XMIN
          X_MAX=XMAX
          Y_MIN=ZMIN
          Y_MAX=ZMAX
        ENDIF
      ELSE IF(iw.EQ.2) THEN !2nd window in FEM environment
        XDC1=0.50*DISP
        XDC2=0.98*DISP
        YDC1=0.51*DISP
        YDC2=0.99*DISP
        X_MIN=YMIN
        X_MAX=YMAX
        Y_MIN=ZMIN
        Y_MAX=ZMAX
      ELSE IF(iw.EQ.3) THEN !3rd window in FEM environment
        XDC1=0.50*DISP
        XDC2=0.98*DISP
        YDC1=0.0
        YDC2=0.48*DISP
        X_MIN=XMIN
        X_MAX=XMAX
        Y_MIN=YMIN
        Y_MAX=YMAX
      ELSE IF(iw.EQ.4) THEN !4th window in FEM environment
        XDC1=0.10*DISP
        XDC2=0.90*DISP
        YDC1=0.20*DISP
        YDC2=1.00*DISP
        X_MIN=-1.0
        X_MAX= 1.0
        Y_MIN=-1.0
        Y_MAX= 1.0
      ELSE IF(iw.EQ.20) THEN !angiography o/p from ANG environment
        XDC1=0.50*XDISP
        XDC2=0.99*XDISP
        YDC1=YDISP-0.69*XDISP
        YDC2=YDISP
        X_MIN=-0.14
        X_MAX= 1.1
        Y_MIN=-3.2
        Y_MAX= 3.2
      ELSE IF(iw.EQ.62) THEN !histogram plot from GEN environment
        XDC1=0.06*XDISP
        XDC2=0.55*XDISP
        YDC1=YDISP-0.55*XDISP
        YDC2=YDISP-0.06*XDISP
        X_MIN=-0.1
        X_MAX= 1.1
        Y_MIN=-0.1
        Y_MAX= 1.1
      ELSE IF(iw.EQ.63) THEN !line/area plot from GEN environment
        XDC1=0.09*XDISP
        XDC2=0.58*XDISP
        YDC1=YDISP-0.58*XDISP
        YDC2=YDISP-0.09*XDISP
        X_MIN=-0.1
        X_MAX= 1.1
        Y_MIN=-0.1
        Y_MAX= 1.1
      ELSE IF(iw.EQ.65) THEN !point/bar plot from GEN environment
        XDC1=0.15*XDISP
        XDC2=0.64*XDISP
        YDC1=YDISP-0.64*XDISP
        YDC2=YDISP-0.15*XDISP
        X_MIN=-0.1
        X_MAX= 1.1
        Y_MIN=-0.1
        Y_MAX= 1.1
      ELSE IF(iw.EQ.66) THEN !2D scatter plot from GEN environment
        XDC1=0.18*XDISP
        XDC2=0.67*XDISP
        YDC1=YDISP-0.67*XDISP
        YDC2=YDISP-0.18*XDISP
        X_MIN=-0.1
        X_MAX= 1.1
        Y_MIN=-0.1
        Y_MAX= 1.1
      ENDIF
      CALL GKS_SVPIP(iw,0,GHIGHR,ERROR,*9999)
      ECAREA(1)=0.0
      ECAREA(2)=XDISP
      ECAREA(3)=0.0
      ECAREA(4)=YDISP
      CALL GKS_INLC(iw,LD1,iw,XWC,YWC,0,ERROR,*9999)
      CALL GKS_RQLC(iw,LD1,INSTAT,IWCNDC,XWC,YWC,ERROR,*9999)
      IF(INSTAT.EQ.GOK) THEN
        INSTAT=1
        ECAREA(1)=XDC1+(XWC-X_MIN)/(X_MAX-X_MIN)*(XDC2-XDC1)
        ECAREA(2)=ECAREA(1)+0.3*XDISP
        ECAREA(3)=YDC1+(YWC-Y_MIN)/(Y_MAX-Y_MIN)*(YDC2-YDC1)
        ECAREA(4)=ECAREA(3)+0.05*YDISP
        CALL GKS_INST(iw,1,IEND-IBEG+1,PROMPT_STRING(IBEG:IEND),1,
     '    80,IEND,0,BLANK,ERROR,*9999)
        CALL GKS_RQST(iw,1,INSTAT1,NCHAR,TEXT_STRING,ERROR,*9999)
        TEXT_STRING=TEXT_STRING(IEND+1:)
      ELSE
        INSTAT=0
      ENDIF

      CALL EXITS('GKSTEXT')
      RETURN
 9999 CALL ERRORS('GKSTEXT',ERROR)
      CALL EXITS('GKSTEXT')
      RETURN 1
      END


      SUBROUTINE INPUT_MODE(iw,ID,CLASS,MODE,ERROR,*)

C#### Subroutine: INPUT_MODE
C###  Description:
C**** resets input mode for workstation iw device ID.
C***  CLASS  can be  'LOCATOR','CHOICE','PICK','VALUATOR'
C***  MODE   can be  'REQUEST','EVENT','SAMPLE'

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:disp00.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
!     Parameter List
      INTEGER ID,iw
      CHARACTER CLASS*(*),ERROR*(*),MODE*(*)
!     Local Variables
      INTEGER IMODE

      CALL ENTERS('INPUT_MODE',*9999)
      IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
        IF(MODE(1:5).EQ.'EVENT') THEN
          IMODE=GEVENT
        ELSE IF(MODE(1:6).EQ.'SAMPLE') THEN
          IMODE=GSAMPL
        ELSE IF(MODE(1:7).EQ.'REQUEST') THEN
          IMODE=GREQU
        ENDIF
        IF(CLASS(1:4).EQ.'PICK') THEN
          CALL GKS_SPKM(iw,ID,IMODE,GECHO,ERROR,*9999)
        ELSE IF(CLASS(1:6).EQ.'CHOICE') THEN
          CALL GKS_SCHM(iw,ID,IMODE,GECHO,ERROR,*9999)
        ELSE IF(CLASS(1:7).EQ.'LOCATOR') THEN
          CALL GKS_SLCM(iw,ID,IMODE,GECHO,ERROR,*9999)
        ELSE IF(CLASS(1:8).EQ.'VALUATOR') THEN
          CALL GKS_SVLM(iw,ID,IMODE,GECHO,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('INPUT_MODE')
      RETURN
 9999 CALL ERRORS('INPUT_MODE',ERROR)
      CALL EXITS('INPUT_MODE')
      RETURN 1
      END


      SUBROUTINE INVIS(iw,ISEG,ERROR,*)

C#### Subroutine: INVIS
C###  Description:
C**** Updates workstation iw's invisibility filter (for phigs wkst only)

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:ntsg00.cmn'
      INCLUDE 'SYS$LIBRARY:PHIGS$DEFS.FOR'
!     Parameter List
      INTEGER ISEG(*),iw
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IEXCL,IINCL(255),NOEXCL,NOINCL,noseg

      CALL ENTERS('INVIS',*9999)
      NOINCL=0
      NOEXCL=0
      DO noseg=1,NTSG
        IF(ISEG(noseg).EQ.1) THEN
          NOINCL=NOINCL+1
          IINCL(NOINCL)=noseg
        ENDIF
      ENDDO
      CALL PHIGS$SET_INVISIBILITY_FILTER(iw,NOINCL,IINCL,NOEXCL,IEXCL)
      CALL PHIGS$REDRAW_ALL_STRUCT(iw,PHIGS$K_CLEAR_ALWAYS)

      CALL EXITS('INVIS')
      RETURN
 9999 CALL ERRORS('INVIS',ERROR)
      CALL EXITS('INVIS')
      RETURN 1
      END


      SUBROUTINE LINE3D(NTDX,XL,ERROR,*)

C#### Subroutine: LINE3D
C###  Description:
C**** Draws line on 3D viewport at rect.cart. world coords.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:fbgr00.cmn'
      INCLUDE 'cmiss$reference:trans00.cmn'
!     Parameter List
      INTEGER NTDX
      REAL*8 XL(20,*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nodx
      REAL Z1(3),Z2(3)

      CALL ENTERS('LINE3D',*9999)
      DO 200 nodx=1,NTDX
        Z1(1)=REAL(XL(nodx,1))
        Z1(2)=REAL(XL(nodx,2))
        Z1(3)=REAL(XL(nodx,3))
        IF(RHTRAN) THEN
          CALL FBIMGPT(5,Z1(1),Z1(2),Z1(3),Z2(1),Z2(2))
          XL(nodx,1)=DBLE(Z2(1))
          XL(nodx,2)=DBLE(Z2(2))
        ELSE IF(LFTRAN) THEN
          CALL FBIMGPT(6,Z1(1),Z1(2),Z1(3),Z2(1),Z2(2))
          XL(nodx,1)=DBLE(Z2(1))
          XL(nodx,2)=DBLE(Z2(2))
        ELSE
          CALL ZZ(Z1,Z2,TRANS)
          XL(nodx,1)=DBLE(Z2(1))
          XL(nodx,2)=DBLE(Z2(2))
          XL(nodx,3)=DBLE(Z2(3))
        ENDIF
 200  CONTINUE

      CALL EXITS('LINE3D')
      RETURN
 9999 CALL ERRORS('LINE3D',ERROR)
      CALL EXITS('LINE3D')
      RETURN 1
      END


      SUBROUTINE LICOLREP(iw,ERROR,*)

C#### Subroutine: LICOLREP
C###  Description:
C***  Lists colour representation thru inquiry functions

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:gks001.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.FOR'
      INCLUDE 'SYS$LIBRARY:PHIGS$DEFS.FOR'
!     Parameter List
      INTEGER iw
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER N,IERR,index,NUM_INDEXES,
     '  LIST_INDEXES(300),NTH_INDEX
      REAL COL(3)

      CALL ENTERS('LICOLREP',*9999)

      WRITE(OP_STRING,'('' iw = '',I2)') iw
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' NINDICES = '',I4)') NINDICES
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      IF(IWKT(iw).EQ.1) THEN      !GKS
        CALL GKS_QECI(iw,1,IERR,NUM_INDEXES,NTH_INDEX,ERROR,*9999)
        DO N=1,NUM_INDEXES
          CALL GKS_QECI(iw,N,IERR,NUM_INDEXES,NTH_INDEX,ERROR,*9999)
          LIST_INDEXES(N)=NTH_INDEX
        ENDDO
        WRITE(OP_STRING,'('' Defined indices:'',/(20I4))')
     '    (LIST_INDEXES(N),N=1,NUM_INDEXES)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF
      DO index=0,NINDICES-1
        IF(IWKT(iw).EQ.1) THEN      !GKS
          CALL GKS_QCR(iw,index,0,IERR,COL,ERROR,*9999)
          CALL GKS_QCR(iw,index,1,IERR,COL,ERROR,*9999)
        ELSE IF(IWKT(iw).EQ.2) THEN !Phigs
          CALL PHIGS$INQ_COLOUR_REP(iw,index,PHIGS$K_VALUE_REALIZED,
     '      IERR,COL)
        ENDIF
        WRITE(OP_STRING,'('' Index='',I4,'' Colours: '',3E12.3)')
     '    index,COL(1),COL(2),COL(3)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDDO

      CALL EXITS('LICOLREP')
      RETURN
 9999 CALL ERRORS('LICOLREP',ERROR)
      CALL EXITS('LICOLREP')
      RETURN 1
      END


      SUBROUTINE LISTRU(NOCO,CO,STRING,ERROR,*)

C#### Subroutine: LISTRU
C###  Description:
C***  Lists phigs data structure thru inquiry functions

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:view00.cmn'
!     Parameter List
      INTEGER NOCO
      CHARACTER CO(*)*(*),STRING*(*),ERROR*(*)
!     Local Variables
      INTEGER ELEM_NUM,ELEM_SIZE,ELEM_TYPE,IBEG,IEND,IERR

      CALL ENTERS('LISTRU',*9999)

      IF(CO(NOCO+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        WRITE(OP_STRING,'(X,A)') 'FEM List Structure'
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ELSE IF(CO(NOCO+1).EQ.'??') THEN
        CALL DOCUM('fe14','doc','LISTRU',ERROR,*9999)
      ELSE
        ELEM_NUM=0
        ELEM_TYPE=1
        DO WHILE (ELEM_NUM.LE.20)
          CALL PHIGS$INQ_ELEM_TYPE_SIZE(ISVIEW,ELEM_NUM,IERR,ELEM_TYPE,
     '      ELEM_SIZE)
          WRITE(OP_STRING,
     '      '('' Isview element'',I3,'' error,type,size'',3I5)')
     '      ELEM_NUM,IERR,ELEM_TYPE,ELEM_SIZE
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELEM_NUM=ELEM_NUM+1
        ENDDO
      ENDIF

      CALL EXITS('LISTRU')
      RETURN
 9999 CALL ERRORS('LISTRU',ERROR)
      CALL EXITS('LISTRU')
      RETURN 1
      END


      SUBROUTINE LOCATOR(INIT,iw,INSTAT,MODE,NECHO,
     '  D_XREF,D_XWC,D_YREF,D_YWC,ERROR,*)

C#### Subroutine: LOCATOR
C###  Description:
C**** Calls GKS locator
C**** INIT specifies whether locator is to be initialised (1:yes,2:no)
C**** iw specifies workstation number (which is also transformation no)
C**** INSTAT is returned as 1 if locate is successful, 0 otherwise
C**** MODE is input mode - 'REQUEST','SAMPLE' or 'EVENT'
C**** NECHO is 0 on entry for no echo
C**** NECHO is 3 on entry for default   prompt/echo
C**** NECHO is 4 on entry for line      prompt/echo
C**** NECHO is 5 on entry for rectangle prompt/echo
C**** NECHO is 6 on entry for digital   prompt/echo
C**** XREF,YREF are initial coords of echo
C**** XWC,YWC are returned world coords of located point

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:disp00.cmn'
      INCLUDE 'cmiss$reference:echo00.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'gx$path:gx.inc'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
!     Parameter List
      INTEGER INIT,INSTAT,iw,NECHO
      REAL*8 D_XREF,D_XWC,D_YREF,D_YWC
      CHARACTER ERROR*(*),MODE*(*)
!     Local Variables
      INTEGER CANCEL,ERR,IERR,IWCNDC,
     '  LD1,NTT,NTRN
      REAL BNDRCT(4),XREF,XWC,YREF,YWC
      DATA LD1/1/

      CALL ENTERS('LOCATOR',*9999)
      XREF=REAL(D_XREF)
      YREF=REAL(D_YREF)
      IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
        BNDRCT(1)=0.0
        BNDRCT(2)=0.0
        BNDRCT(3)=0.0
        BNDRCT(4)=0.0
        CALL LOCATR(gxLOC_POINT,BNDRCT,XREF,YREF,XWC,YWC,CANCEL,ERR)
        IF(err.ne.0) THEN
          ERROR='>>Error from GX'
          GOTO 9999
        ENDIF
        IF(CANCEL.GT.0) THEN
          INSTAT=0
        ELSE
          INSTAT=1
        ENDIF
      ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
        IF(INIT.EQ.1) THEN !initialise locator device
!         must be in request mode to initialize
          CALL GKS_SLCM(iw,LD1,GREQU,GECHO,ERROR,*9999)
          ECAREA(1)=0.0
          ECAREA(2)=XDISP
          ECAREA(3)=0.0
          ECAREA(4)=YDISP
          CALL GKS_INLC(iw,LD1,iw,XREF,YREF,NECHO,ERROR,*9999)
          IF(NECHO.EQ.0) THEN !no echo
            IF(MODE(1:5).EQ.'EVENT') THEN
              CALL GKS_SLCM(iw,1,GEVENT,GNECHO,ERROR,*9999)
            ELSE IF(MODE(1:6).EQ.'SAMPLE') THEN
              CALL GKS_SLCM(iw,1,GSAMPL,GNECHO,ERROR,*9999)
            ELSE IF(MODE(1:7).EQ.'REQUEST') THEN
              CALL GKS_SLCM(iw,1,GREQU,GNECHO,ERROR,*9999)
            ENDIF
          ELSE IF(NECHO.GT.0) THEN !echo
            IF(MODE(1:5).EQ.'EVENT') THEN
              CALL GKS_SLCM(iw,1,GEVENT,GECHO,ERROR,*9999)
            ELSE IF(MODE(1:6).EQ.'SAMPLE') THEN
              CALL GKS_SLCM(iw,1,GSAMPL,GECHO,ERROR,*9999)
            ELSE IF(MODE(1:7).EQ.'REQUEST') THEN
              CALL GKS_SLCM(iw,1,GREQU,GECHO,ERROR,*9999)
            ENDIF
          ENDIF
          CALL GKS_QENTN(1,IERR,NTT,NTRN,ERROR,*9999)
          IF(iw.ne.NTRN) CALL GKS_SVPIP(iw,NTRN,GHIGHR,ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,*)' iw=',iw,
     '        ' call to gqentn returns ierr,ntt,ntrn:',IERR,NTT,NTRN
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
        IF(MODE(1:6).EQ.'SAMPLE') THEN
          CALL GKS_SMLC(iw,LD1,iw,XWC,YWC,ERROR,*9999)
        ELSE IF(MODE(1:7).EQ.'REQUEST') THEN
          CALL GKS_RQLC(iw,LD1,INSTAT,IWCNDC,XWC,YWC,ERROR,*9999)
        ENDIF
        IF(INSTAT.EQ.GOK) THEN
          INSTAT=1
        ELSE
          INSTAT=0
        ENDIF
      ENDIF
      IF(DOP) THEN
        WRITE(OP_STRING,'('' IWCNDC='',I2)') IWCNDC
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '('' XWC='',E11.3,'' YWC='',E11.3)') XWC,YWC
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      D_XWC=DBLE(XWC)
      D_YWC=DBLE(YWC)

      CALL EXITS('LOCATOR')
      RETURN

 9999 CALL ERRORS('LOCATOR',ERROR)
      CALL EXITS('LOCATOR')
      RETURN 1
      END


      SUBROUTINE OPEN_PRINT_FILE(iw,TYPE,ERROR,*)

C#### Subroutine: OPEN_PRINT_FILE
C###  Description:
C**** Activates print workstation for GKS or PHIGS output.
C**** If Phigs the phigs viewing transformations are set up.
C**** TYPE can be 'POSTSCRIPT' (print wkst is 15 for GKS, 16 for Phigs)
C****          or 'METAFILE'   (print wkst is 17 )

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:cbzm00.cmn'
      INCLUDE 'cmiss$reference:ntsg00.cmn'
      INCLUDE 'cmiss$reference:phig00.cmn'
      INCLUDE 'cmiss$reference:view00.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.FOR'
      INCLUDE 'SYS$LIBRARY:PHIGS$DEFS.FOR'
!     Parameter List
      INTEGER iw
      CHARACTER ERROR*(*),TYPE*(*)
!     Local Variables
      INTEGER i,ISTATUS,IZ,j,MODE_PROJ
      CHARACTER GRAPHICS*5
      RECORD/phigs$typ_colour_val/CLR_VAL,CLR_BACK_VAL,CLR_SPEC_VAL,
     '  CBR_SPEC_VAL

      CALL ENTERS('OPEN_PRINT_FILE',*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,'('' iw='',I3,'' IWKS(iw)='',I2,'//
     '    ''' IWKT(iw)='',I2)') iw,IWKS(iw),IWKT(iw)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(IWKT(iw).EQ.1) THEN      !GKS
        GRAPHICS='GKS'
      ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS
        GRAPHICS='PHIGS'
      ENDIF

      IF(GRAPHICS(1:3).EQ.'GKS') THEN
        IF(TYPE(1:8).EQ.'METAFILE') THEN
          CALL ACWK(17,1,ERROR,*9999)

        ELSE IF(TYPE(1:10).EQ.'POSTSCRIPT') THEN
          CALL ACWK(15,1,ERROR,*9999)

C         copy current indices 1-8 from iw to 15
Commented out by PJH on 1-Jun-91 since INDEX now defined separately for each
c         segment
c         DO index=1,8
c           CALL GKS_QPLR(iw,index,GSET,IERR,LINE_TYPE,
c    '        LINE_WIDTH,LINE_COLOUR_index,ERROR,*9999)
c           CALL GKS_QPMR(iw,index,GSET,IERR,MARKER_TYPE,
c    '        MARKER_SIZE,MARKER_COLOUR_INDEX,ERROR,*9999)
c           CALL GKS_QTXR(iw,index,GSET,IERR,TEXT_FONT,TEXT_PRECISION,
c    '        TEXT_EXPFAC,TEXT_SPACING,TEXT_COLOUR_INDEX,ERROR,*9999)
c           TEXT_FONT=-105
c           IF(ierr.ne.0) THEN
c             ERROR='Error return from gq**r, ierr='//
C          CFROMI(IERR,'(I5)')
c             WRITE(IOOP,'(1X,A)') ERROR(1:40)
c           ENDIF
c           CALL GKS_SPLR(15,index,LINE_TYPE,LINE_WIDTH,1,ERROR,*9999)
c           CALL GKS_SPMR(15,index,MARKER_TYPE,MARKER_SIZE,1,
C    '        ERROR,*9999)
c           CALL GKS_STXR(15,index,TEXT_FONT,TEXT_PRECISION,
c    '        TEXT_EXPFAC,TEXT_SPACING,1,ERROR,*9999)
c         ENDDO
          IF(iw.LE.3) THEN
            IZ=IZOOM(iw)
            IF(DOP) THEN
              WRITE(OP_STRING,
     '          '('' iw='',I2,'' IZOOM(iw)='',I2)') iw,IZ
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,
     '          '('' XNDC:'',4E12.3)') (XNDC(IZ,i,iw),i=1,4)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            CALL GKS_SWKWN(15,XNDC(IZ,1,iw),XNDC(IZ,2,iw),
     '        XNDC(IZ,3,iw),XNDC(IZ,4,iw),ERROR,*9999)
          ELSE IF(iw.EQ.10) THEN
C           set up appropriate window,viewport and workstation viewport
C           CALL GKS_SWKWN(15,0.,1.,0.5-0.23/0.49/2.,0.5+0.23/0.49/2.,
C    '        ERROR,*9999)
          ENDIF
        ENDIF

      ELSE IF(GRAPHICS(1:5).EQ.'PHIGS') THEN
        CALL ACWK(16,1,ERROR,*9999)
        MODE_PROJ=PHIGS$K_PARALLEL
        CALL PHIGS$EVAL_VIEW_ORIEN_MATRIX3(VIEW_REF_PT_NEW,
     '     VIEW_PLANE_NEW,VIEW_UP_NEW,ISTATUS,A_ORIENT)
        IF(DOP) THEN
          DO i=1,4
            WRITE(OP_STRING,'('' A_ORIENT('',I1,'',j): '',4E12.3)')
     '        i,(A_ORIENT(i,j),j=1,4)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDIF
        IF(istatus.ne.0) THEN
          WRITE(OP_STRING,*)
     '      ' ---Error from eval_view_orient_matrix=',ISTATUS
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ELSE
          CALL PHIGS$EVAL_VIEW_MAP_MATRIX3(WINDOW_NEW,VIEWPORT,
     '      MODE_PROJ,PROJ_REF_PT_NEW,VIEW_PLANE_DIST_NEW,
     '      BACK_PLANE_DIST_NEW,FRONT_PLANE_DIST_NEW,ISTATUS,A_MAP)
          IF(DOP) THEN
            DO i=1,4
              WRITE(OP_STRING,'(''    A_MAP('',I1,'',j): '',4E12.3)')
     '          i,(A_MAP(i,j),j=1,4)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDIF
          IF(istatus.ne.0) THEN
            WRITE(OP_STRING,*)
     '        ' ---Error from eval_view_map_matrix=',ISTATUS
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE
            CALL PHIGS$SET_VIEW_REP3(16,1,A_ORIENT,A_MAP,NPC_CLIP,
     '        PHIGS$K_NOCLIP,PHIGS$K_NOCLIP,PHIGS$K_NOCLIP)
          ENDIF
        ENDIF

        IF(TYPE(1:8).EQ.'METAFILE') THEN
        ELSE IF(TYPE(1:10).EQ.'POSTSCRIPT') THEN
C         copy current indices 1-8 from iw to 16
c         DO index=1,8
c           CALL PHIGS$INQ_PLINE_REP(iw,index,PHIGS$K_VALUE_SET,IERR,
c    '        LINE_TYPE,LINE_WIDTH,LINE_COLOUR_INDEX)
c           CALL PHIGS$INQ_PMARKER_REP(iw,index,PHIGS$K_VALUE_SET,IERR,
c    '        MARKER_TYPE,MARKER_SIZE,MARKER_COLOUR_INDEX)
c           CALL PHIGS$INQ_TEXT_REP(iw,index,PHIGS$K_VALUE_SET,IERR,
c    '        TEXT_FONT,TEXT_PRECISION,
c    '        TEXT_EXPFAC,TEXT_SPACING,TEXT_COLOUR_INDEX)
c           TEXT_FONT=-105
c           CALL PHIGS$INQ_EXT_INT_REP(iw,index,PHIGS$K_VALUE_SET,IERR,
c    '        INTERIOR_STYLE,INTERIOR_BACK_STYLE,
c    '        INT_STYLE_INDEX,INT_BACK_STYLE,
c    '        ICOLOUR_TYPE,CLR_VAL,ICOLOUR_BACK_TYPE, CLR_BACK_VAL,
c    '        ISHAD,IBACKSHAD,LIGHT,LIGHT_BACK,
c    '        AMB,DIFF,SPEC,ICOLOUR_SPEC_TYPE,CLR_SPEC_VAL,SPEC_EXPO,
C    '        TRANSPAR,
c    '        ABB,DBFF,SBEC,IBOLOUR_SPEC_TYPE,CBR_SPEC_VAL,SBEC_EXPO,
C    '        TBANSPAR,
c    '        IAPPROX_TYPE,IAPPROX_VAL,ITRIM_TYPE,ITRIM_VAL)
c           IF(ierr.ne.0) THEN
c             ERROR='Error return from inquire rep, ierr='//
C          CFROMI(IERR,'(I5)')
c             GOTO 9999
c           ENDIF
c           CALL PHIGS$SET_PLINE_REP(16,index,LINE_TYPE,
c    '        LINE_WIDTH,1)
c           CALL PHIGS$SET_PMARKER_REP(16,index,MARKER_TYPE,
c    '        MARKER_SIZE,1)
c           CALL PHIGS$SET_TEXT_REP(16,index,TEXT_FONT,TEXT_PRECISION,
c    '        TEXT_EXPFAC,TEXT_SPACING,1)
c           CALL PHIGS$SET_EXT_INT_REP(16,index,
c    '        INTERIOR_STYLE,INTERIOR_BACK_STYLE,
c    '        INT_STYLE_index,INT_BACK_STYLE,
c    '        ICOLOUR_TYPE,CLR_VAL,ICOLOUR_BACK_TYPE,CLR_BACK_VAL,
c    '        ISHAD,IBACKSHAD,LIGHT,LIGHT_BACK,
c    '        AMB,DIFF,SPEC,ICOLOUR_SPEC_TYPE,CLR_SPEC_VAL,SPEC_EXPO,
C    '        TRANSPAR,
c    '        ABB,DBFF,SBEC,IBOLOUR_SPEC_TYPE,CBR_SPEC_VAL,SBEC_EXPO,
C    '        TBANSPAR,
c    '        IAPPROX_TYPE,IAPPROX_VAL,ITRIM_TYPE,ITRIM_VAL)
c         ENDDO
        ENDIF

      ENDIF

      CALL EXITS('OPEN_PRINT_FILE')
      RETURN
 9999 CALL ERRORS('OPEN_PRINT_FILE',ERROR)
      CALL EXITS('OPEN_PRINT_FILE')
      RETURN 1
      END


      SUBROUTINE OPEN_SEGMENT(ISEGNUM,ISEG,iw,CLABEL,index,INDEX_OLD,
     '  NLABEL,IVIS,CSEG,ERROR,*)

C#### Subroutine: OPEN_SEGMENT
C###  Description:
C**** Opens graphics segment ISEGNUM. If ISEGNUM does not already exist
C**** (has the value 0), then NTSG is incremented and a new segment is created.
C**** Otherwise the old segment structure is replaced by the new structure.
C**** All subsequent call to graphics primitives, etc will be part of this
C**** structure. Note that GKS segments cannot be nested, so a call to
C**** OPEN_SEGMENT must be followed by a call to CLOSE_SEGMENT before the next
C**** OPEN_SEGMENT.
C**** When created CSEG(ISEGNUM) stores a string comprising the first 48
C**** characters of CLABEL together with index(4 chars) & NLABEL(5 chars)
C**** then a '/' and the iw number(2 chars).
C**** index is the bundle table index for the current segment and if an old
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
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.FOR'
      INCLUDE 'SYS$LIBRARY:PHIGS$DEFS.FOR'
!     Parameter List
      INTEGER index,INDEX_OLD,ISEG(*),ISEGNUM,IVIS,iw,NLABEL
      CHARACTER CLABEL*(*),CSEG(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER IFROMC
      CHARACTER CFROMI*5,CHAR2*2,CHAR4*4,CHAR5*5,CLABEL2*(52)

      CALL ENTERS('OPEN_SEGMENT',*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,
     '    '('' iw='',I3,'' IWKS(iw)='',I2,'' IWKT(iw)='',I2,'
     '    //''' IWKG(iw)='',I2)') iw,IWKS(iw),IWKT(iw),IWKG(iw)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' ISEGNUM='',I5)') ISEGNUM
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(ISEGNUM.EQ.0) THEN !segment does not already exist
        NTSG=NTSG+1
        ISEGNUM=NTSG
      ELSE !recover index, delete old segment and use the old ISEG,CSEG labels
        INDEX_OLD=IFROMC(CSEG(ISEGNUM)(49:52))
        IF(DOP) THEN
          WRITE(OP_STRING,'('' INDEX='',I3)') index
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        IF(IWKT(iw).EQ.1) THEN      !GKS
          IF(DOP) THEN
            WRITE(OP_STRING,
     '        '('' Delete segment number '',I3)') ISEGNUM
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CALL GKS_DSG(ISEGNUM,ERROR,*9999)
        ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS
          CALL PHIGS$DEL_STRUCT(ISEGNUM)
c         CALL PHIGS$POST_STRUCT(iw,ISVIEW,1.0)  !why is this necessary?
        ENDIF
      ENDIF

      ISEG(ISEGNUM)=2
      IF(IVIS.EQ.2) ISEG(ISEGNUM)=1

C     open structure for subsequent graphic, set insert mode
      IF(IWKT(iw).EQ.1) THEN      !GKS
        IF(DOP) THEN
          WRITE(OP_STRING,
     '      '('' Create segment isegnum='',I5,'' on iw='',I2)')
     '      ISEGNUM,iw
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        CALL GKS_CRSG(ISEGNUM,ERROR,*9999)
        IF(ISEG(ISEGNUM).EQ.1) THEN
          CALL GKS_SVIS(ISEGNUM,GINVIS,ERROR,*9999)
        ELSE IF(ISEG(ISEGNUM).EQ.2) THEN
          CALL GKS_SVIS(ISEGNUM,GVISI,ERROR,*9999)
        ENDIF
        CALL GKS_SDTEC(ISEGNUM,GUNDET,ERROR,*9999)
      ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS
        CALL PHIGS$SET_EDIT_MODE(PHIGS$K_EDIT_INSERT)
        CALL PHIGS$OPEN_STRUCT(ISEGNUM)
        CALL PHIGS$ADD_NAMES_TO_SET(1,ISEGNUM)
      ENDIF

      CHAR2=CFROMI(iw,'(I2)')
      CHAR4=CFROMI(index,'(I4)')
      CHAR5=CFROMI(NLABEL,'(I5)')
      CLABEL2(1:)=CLABEL
      CSEG(ISEGNUM)=CLABEL2(1:48)//CHAR4(1:4)//CHAR5(1:5)//'/'//
     '  CHAR2(1:2)

      CALL EXITS('OPEN_SEGMENT')
      RETURN
 9999 CALL ERRORS('OPEN_SEGMENT',ERROR)
      CALL EXITS('OPEN_SEGMENT')
      RETURN 1
      END


      SUBROUTINE PHIG(iw,ERROR,*)

C#### Subroutine: PHIG
C###  Description:
C**** Defines initial transformations and view matrices for PHIGS.
C**** A_TRANS is the 4*4 transformation matrix containing rotations
C**** shifts and scaling of an object in world coords
C**** A_ORIENT is a 4*4 matrix indicating the orientation of the
C**** world coords in the viewing reference coords (z-axis into the
C**** screen, y-axis vertical)
C**** A_MAP is a 4*4 matrix specifying parallel or perspective views
C**** together with clipping planes

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:gks000.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'cmiss$reference:phig00.cmn'
      INCLUDE 'cmiss$reference:trac00.cmn'
      INCLUDE 'cmiss$reference:trans00.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
      INCLUDE 'SYS$LIBRARY:PHIGS$DEFS.FOR'
!     Parameter List
      INTEGER iw
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,ISTATUS,j,nj

      CALL ENTERS('PHIG',*9999)

C *** Reference point, rotation angles, scaling & shift vector
C *** for world coord transformation
C     FIXED_PT(1)=0.5*(XMIN+XMAX)
C     FIXED_PT(2)=0.5*(YMIN+YMAX)
      FIXED_PT(1)=0.0
      FIXED_PT(2)=0.0
      FIXED_PT(3)=0.0
      IF(DOP) THEN
        WRITE(OP_STRING,'('' FIXED_PT:'',3E12.3)')
     '    (FIXED_PT(nj),nj=1,3)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      DO j=1,3
        ANGLE(j)=0.0
        SCALE(j)=1.0
        SHIFT(j)=0.0
      ENDDO
      CALL PHIGS$BUILD_XFORM_MATRIX3(FIXED_PT,SHIFT,
     '  ANGLE(1),ANGLE(2),ANGLE(3),SCALE,ISTATUS,A_TRANS)
      IF(DOP) THEN
        WRITE(OP_STRING,'('' Transformation matrix:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO i=1,4
          WRITE(OP_STRING,'(4E12.3)') (A_TRANS(i,j),j=1,4)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

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
      CALL PHIGS$EVAL_VIEW_ORIEN_MATRIX3(VIEW_REF_PT,VIEW_PLANE,
     '  VIEW_UP,ISTATUS,A_ORIENT)
      CALL PHIGS$EVAL_VIEW_ORIEN_MATRIX3(VIEW_REF_PT,VIEW_PLANE,
     '  VIEW_UP,ISTATUS,A_ORIENT_NEW)
      IF(DOP) THEN
        WRITE(OP_STRING,'('' View orientation matrix:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO i=1,4
          WRITE(OP_STRING,'(4E12.3)') (A_ORIENT(i,j),j=1,4)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

C *** Window, projection ref pt and viewplane, backplane & frontplane
C *** positions in viewing reference coords and viewport in norm proj
C *** coords for mapping to normalized projection coord system
      WINDOW(1)=XMIN-VIEW_REF_PT(1)
      WINDOW(2)=XMAX-VIEW_REF_PT(1)
      WINDOW(3)=YMIN-VIEW_REF_PT(2)
      WINDOW(4)=YMAX-VIEW_REF_PT(2)
      PROJ_REF_PT(1)= 0.0
      PROJ_REF_PT(2)= 0.0
      PROJ_REF_PT(3)=3.0*DIAG
      BACK_PLANE_DIST =-DIAG
      FRONT_PLANE_DIST=DIAG
      VIEW_PLANE_DIST = 0.0
c     Z_MIN= 1.E6
c     Z_MAX=-1.E6
c     DO np=1,NPT
c       IF(ITYP1.EQ.1.OR.KTYP1.EQ.11) THEN
c         Z=XP(1,1,3,np)
c       ELSE
c         Z=ZP(1,1,nh,np)
c       ENDIF
c       IF(Z.LT.Z_MIN) Z_MIN=Z
c       IF(Z.GT.Z_MAX) Z_MAX=Z
c     ENDDO
c     IF(Z_MAX.LE.Z_MIN) THEN
c       Z_MAX= 1.0
c       Z_MIN=-1.0
c     ENDIF
c     IF(Z_MIN.GT.0.0) Z_MIN=0.0
c     BACK_PLANE_DIST =Z_MIN
c     FRONT_PLANE_DIST=Z_MAX
c     VIEW_PLANE_DIST =0.5*(Z_MIN+Z_MAX)
      IF(DOP) THEN
        WRITE(OP_STRING,'('' Frontplane at '',E12.3,'
     '    //''' Backplane at '',E12.3)')
     '    FRONT_PLANE_DIST,BACK_PLANE_DIST
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      VIEWPORT(1)= 0.0
      VIEWPORT(2)= 1.0
      VIEWPORT(3)= 0.0
      VIEWPORT(4)= 1.0
C     VIEWPORT(5)=-1.0   !old phigs version
C     VIEWPORT(6)= 0.0   !old phigs version
      VIEWPORT(5)= 0.0   !7-feb-1990
      VIEWPORT(6)= 1.0   !7-feb-1990
      CALL PHIGS$EVAL_VIEW_MAP_MATRIX3(WINDOW,VIEWPORT,
     '  PHIGS$K_PARALLEL,
     '  PROJ_REF_PT,VIEW_PLANE_DIST,BACK_PLANE_DIST,
     '  FRONT_PLANE_DIST,ISTATUS,A_MAP)
      IF(DOP) THEN
        WRITE(OP_STRING,'('' Mapping matrix:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO i=1,4
          WRITE(OP_STRING,'(4E12.3)') (A_MAP(i,j),j=1,4)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

C *** Clipping limits in normalized projection coords
      NPC_CLIP(1)= 0.0
      NPC_CLIP(2)= 1.0
      NPC_CLIP(3)= 0.0
      NPC_CLIP(4)= 1.0
C     NPC_CLIP(5)=-1.0   !old phigs version
C     NPC_CLIP(6)= 0.0   !old phigs version
      NPC_CLIP(5)= 0.0   !7-feb-1990
      NPC_CLIP(6)= 1.0   !7-feb-1990
      CALL PHIGS$SET_VIEW_REP3(iw,1,A_ORIENT,A_MAP,NPC_CLIP,
     '  PHIGS$K_NOCLIP,PHIGS$K_NOCLIP,PHIGS$K_NOCLIP)
C *** Initialize
      VIEW_PLANE_DIST_NEW =VIEW_PLANE_DIST
      BACK_PLANE_DIST_NEW =BACK_PLANE_DIST
      FRONT_PLANE_DIST_NEW=FRONT_PLANE_DIST
      DO nj=1,3
        PROJ_REF_PT_NEW(nj)=PROJ_REF_PT(nj)
        VIEW_REF_PT_NEW(nj)=VIEW_REF_PT(nj)
        VIEW_PLANE_NEW(nj)=VIEW_PLANE(nj)
        VIEW_UP_NEW(nj)=VIEW_UP(nj)
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


      SUBROUTINE PICK(iw,ID,MODE,INSTAT,IPICKRET,IPICKID,ERROR,*)

C#### Subroutine: PICK
C###  Description:
C**** Calls GKS pick
C**** MODE can be 'REQUEST' or 'EVENT'
C**** IPICKRET is returned segment on exit if in request mode
C**** IPCKID   is returned pick identifier if in request mode
C**** INSTAT   ir returned 1 if successful input triggered

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:disp00.cmn'
      INCLUDE 'cmiss$reference:echo00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:gks000.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
!     Parameter List
      INTEGER ID,INSTAT,IPICKID,IPICKRET,iw
      CHARACTER ERROR*(*),MODE*(*)
!     Local variables
      INTEGER ERR

      CALL ENTERS('PICK',*9999)

      CALL INPUT_MODE(iw,ID,'PICK',MODE,ERROR,*9999)
      IF(iw.ne.3) THEN !gks
!       Set pick area to full screen

        IF(MODE(1:5).EQ.'EVENT') THEN
        ELSE IF(MODE(1:7).EQ.'REQUEST') THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            ECAREA(1)=0.0
            ECAREA(2)=0.0
            ECAREA(3)=0.0
            ECAREA(4)=0.0
            CALL FPICK(ECAREA,IPICKRET,INSTAT,ERR)
            IF(err.ne.0) THEN
              ERROR='>>Error from GX'
              GOTO 9999
            ENDIF
            IF(INSTAT.EQ.0) THEN
              INSTAT=1
            ELSE
              INSTAT=0
            ENDIF
            IPICKID=IPICKRET
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            ECAREA(1)=0.0
            ECAREA(2)=XDISP
            ECAREA(3)=0.0
            ECAREA(4)=YDISP
            CALL GKS_INPK(iw,ID,GNPICK,1,1,1,0,' ',ERROR,*9999)
            CALL GKS_RQPK(iw,ID,INSTAT,IPICKRET,IPICKID,ERROR,*9999)
            IF(INSTAT.EQ.GOK) THEN
              INSTAT=1
            ELSE
              INSTAT=0
            ENDIF
          ENDIF
        ENDIF
      ENDIF

      CALL EXITS('PICK')
      RETURN

 9999 CALL ERRORS('PICK',ERROR)
      CALL EXITS('PICK')
      RETURN 1
      END


      SUBROUTINE POLYLINE(IBUNDLE,iw,NTPTS,D_PTS,ERROR,*)

C#### Subroutine: POLYLINE
C###  Description:
C**** Draws a polyline on iw with index IBUNDLE.
C**** D_PTS(1..3,nopts) contains the REAL*8 3D coords of each point.
C**** IF iw=4 coords are curvilinear coords, else rect. cart.
C**** If IBUNDLE is 0 the primitive will use the previously
C**** defined polyline index.
C**** NOTE: If iw is 15 or 16 (postscript) IBUNDLE is reset to be black.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:cmis00.cmn'
      INCLUDE 'cmiss$reference:curr00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:map000.cmn'
!     Parameter List
      INTEGER IBUNDLE,iw,NTPTS
      REAL*8 D_PTS(3,*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nj,nopts

      CALL ENTERS('POLYLINE',*9999)
      IF(NTPTS.EQ.0) GOTO 9998

      IF(DOP) THEN
        WRITE(OP_STRING,
     '    '('' iw='',I2,'' NTPTS='',I6)') iw,NTPTS
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nopts=1,NTPTS
          WRITE(OP_STRING,
     '      '('' D_PTS(nj,'',I6,''): '',3E12.3)')
     '      nopts,(D_PTS(nj,nopts),nj=1,NJT)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

      IF(iw.EQ.1) THEN !plot x against y
        IF(ibundle.ne.0) CALL GKS_SPLI(IBUNDLE,ERROR,*9999)
!       IF(NJT.EQ.2) THEN !plot x against y
          CALL GKS_PL(1,2,NTPTS,D_PTS,ERROR,*9999)
!       ELSE IF(NJT.EQ.3) THEN !plot x against z
!         CALL GKS_PL(1,3,NTPTS,D_PTS,ERROR,*9999)
!       ENDIF

      ELSE IF(iw.EQ.2) THEN !plot y against z
        IF(ibundle.ne.0) CALL GKS_SPLI(IBUNDLE,ERROR,*9999)
        CALL GKS_PL(2,3,NTPTS,D_PTS,ERROR,*9999)

      ELSE IF(iw.EQ.3) THEN !plot x against z
        IF(ibundle.ne.0) CALL GKS_SPLI(IBUNDLE,ERROR,*9999)
        CALL GKS_PL(1,3,NTPTS,D_PTS,ERROR,*9999)

!     ELSE IF(iw.EQ.3) THEN
!       IF(ibundle.ne.0) CALL PHIGS$SET_PLINE_INDEX(IBUNDLE)
!       CALL PHIGS_POLYLINE3(NTPTS,D_PTS,ERROR,*9999)

      ELSE IF(iw.EQ.4) THEN !Map window
        IF(ibundle.ne.0) CALL GKS_SPLI(IBUNDLE,ERROR,*9999)
        IF(PROJEC(1:11).EQ.'RECTANGULAR') THEN !points x and y
          CALL GKS_PL(1,2,NTPTS,D_PTS,ERROR,*9999)
        ELSE IF(PROJEC(1:2).EQ.'XI') THEN !assume pts are in Xi coords
          DO nopts=1,NTPTS
            D_PTS(1,nopts)=-1.0D0+2.0D0*(DBLE(MXI1-1)+
     '        (D_PTS(1,nopts))/MAX_XI)
            D_PTS(2,nopts)=-1.0D0+2.0D0*(DBLE(MXI2-1)+(D_PTS(2,nopts))
     '        /MAX_XI)
          ENDDO
          CALL GKS_PL(1,2,NTPTS,D_PTS,ERROR,*9999)
        ELSE  !assume points are in polar coords
           CALL MAP4(NTPTS,D_PTS,D_PTS,ERROR,*9999)
           CALL GKS_PL(1,2,NTPTS,D_PTS,ERROR,*9999)
c          DO nopts=1,NTPTS-1
c            DO nj=1,NJT
c              TWO_PTS(nj,1)=D_PTS(nj,nopts  )
c              TWO_PTS(nj,2)=D_PTS(nj,nopts+1)
c            ENDDO
!           Map4 converts mu,theta coords to map coords (-1<.<1)
c            CALL MAP4(2,TWO_PTS,TWO_PTS,ERROR,*9999)
! PJH 20aug92 IF(DABS(TWO_PTS(1,1)-TWO_PTS(2,2)).LT.0.5D0) THEN !ok to plot
c            CALL GKS_PL(1,2,2,TWO_PTS,ERROR,*9999)
! PJH 20aug92 ENDIF
c          ENDDO
        ENDIF

      ELSE IF(iw.EQ.5.OR.iw.EQ.6) THEN
!       frame grabber windows
c PJH 22-Jul-92 need to sort out real*8 here
c       IF(PROJEC(1:6).EQ.'HAMMER') THEN
c         IF(JTYP3.EQ.1) THEN !points x and y are xid coords
c           DO nopts=1,NTPTS
c             X(nopts)=REAL(D_PTS(1,nopts))
c             Y(nopts)=REAL(D_PTS(2,nopts))
c           ENDDO
c         ELSE !has been transformed to polar coords
c           CALL MAP4(NTPTS,D_PTS,D_PTS,ERROR,*9999)
c         ENDIF
c         COLOUR=254
c ***     Scale to pixel window size
c         DO N=1,NTPTS
c           XX=(X(N)+1.0)/2.0
c           YY=(Y(N)+1.0)/2.0
c           X(N)=(1.0-XX)*WINDOW(1)+XX*WINDOW(3)
c           Y(N)=(1.0-YY)*WINDOW(2)+YY*WINDOW(4)
c         ENDDO
c       ELSE
c         DO nopts=1,NTPTS
c           this call removed on 26/4/91 by Steve to make scale correct
c           CALL FBIMGPT(iw,D_PTS(1,nopts),D_PTS(2,nopts),
c             D_PTS(3,nopts),X(nopts),Y(nopts))
c           X(nopts)=REAL(D_PTS(1,nopts))
c           Y(nopts)=REAL(D_PTS(2,nopts))
c         ENDDO
c       ENDIF
c       CALL FBPL(iw,NTPTS,X,Y)

      ELSE IF(iw.EQ.7) THEN
        IF(ibundle.ne.0) CALL GKS_SPLI(IBUNDLE,ERROR,*9999)
!       plot Y against z
        CALL GKS_PL(2,3,NTPTS,D_PTS,ERROR,*9999)

      ELSE IF(iw.EQ.10) THEN
        IF(ibundle.ne.0) CALL GKS_SPLI(IBUNDLE,ERROR,*9999)
!       Time history plots
        DO nopts=1,NTPTS
          CALL MAP10(1,D_PTS(2,nopts))
        ENDDO
        CALL GKS_PL(1,2,NTPTS,D_PTS,ERROR,*9999)

      ELSE IF(iw.EQ.11) THEN
        IF(ibundle.ne.0) CALL GKS_SPLI(IBUNDLE,ERROR,*9999)
        ! Section plots
        DO nopts=1,NTPTS
          CALL MAP10(1,D_PTS(2,nopts))
        ENDDO
        CALL GKS_PL(1,2,NTPTS,D_PTS,ERROR,*9999)

      ELSE IF(iw.EQ.12) THEN !fibre plot
        IF(ibundle.ne.0) CALL GKS_SPLI(IBUNDLE,ERROR,*9999)
        CALL GKS_PL(1,2,NTPTS,D_PTS,ERROR,*9999)

      ELSE IF(iw.EQ.13) THEN !sheet plot
        IF(ibundle.ne.0) CALL GKS_SPLI(IBUNDLE,ERROR,*9999)
        CALL GKS_PL(1,2,NTPTS,D_PTS,ERROR,*9999)

      ELSE IF(iw.EQ.15) THEN !GKS postscript
        IF(ibundle.ne.0) THEN
C CPB 1/11/92 removed as using colour postscript
C          IF(IBUNDLE.GT.8) THEN !reset to black for postscript
C            IBUNDLE=IBUNDLE-8
C            IF(DOP) WRITE(IOOP,'('' Reset IBUNDLE='',I2)') IBUNDLE
C          ENDIF
          CALL GKS_SPLI(IBUNDLE,ERROR,*9999)
        ENDIF
        IF(NJT.EQ.2) THEN !plot x against y
          CALL GKS_PL(1,2,NTPTS,D_PTS,ERROR,*9999)
        ELSE IF(NJT.EQ.3) THEN !plot x against z
          CALL GKS_PL(1,3,NTPTS,D_PTS,ERROR,*9999)
        ENDIF

      ELSE IF(iw.EQ.16) THEN !PHIGS postscript
        IF(ibundle.ne.0) THEN
C CPB 1/11/92 removed as using colour postscript
C          IF(IBUNDLE.GT.8) THEN !reset to black for postscript
C            IBUNDLE=IBUNDLE-8
C            IF(DOP) WRITE(IOOP,'('' Reset IBUNDLE='',I2)') IBUNDLE
C          ENDIF
          CALL PHIGS$SET_PLINE_INDEX(IBUNDLE)
        ENDIF
        CALL PHIGS_POLYLINE3(NTPTS,D_PTS,ERROR,*9999)

      ELSE IF(iw.EQ.21.OR.iw.EQ.22.OR.iw.EQ.23
     '  .OR.iw.EQ.30.OR.iw.EQ.31.OR.iw.EQ.32.OR.iw.EQ.34.OR.iw.EQ.35
     '  .OR.iw.EQ.40.OR.iw.EQ.41.OR.iw.EQ.42.OR.iw.EQ.43.OR.iw.EQ.44
     '  .OR.iw.EQ.45.OR.iw.EQ.46.OR.iw.EQ.47
     '  .OR.iw.EQ.50.OR.iw.EQ.51
     '  .OR.iw.EQ.60.OR.iw.EQ.61.OR.iw.EQ.62.OR.iw.EQ.63.OR.iw.EQ.64
     '  .OR.iw.EQ.65.OR.iw.EQ.66.OR.iw.EQ.68.OR.iw.EQ.69) THEN !plots
        IF(ibundle.ne.0) CALL GKS_SPLI(IBUNDLE,ERROR,*9999)
        !plot x against y
        CALL GKS_PL(1,2,NTPTS,D_PTS,ERROR,*9999)

      ELSE IF(iw.EQ.67) THEN
        IF(ibundle.ne.0) CALL PHIGS$SET_PLINE_INDEX(IBUNDLE)
        CALL PHIGS_POLYLINE3(NTPTS,D_PTS,ERROR,*9999)

      ENDIF

 9998 CALL EXITS('POLYLINE')
      RETURN
 9999 CALL ERRORS('POLYLINE',ERROR)
      CALL EXITS('POLYLINE')
      RETURN 1
      END


      SUBROUTINE POLYMARKER(IBUNDLE,iw,NTPTS,D_PTS,ERROR,*)

C#### Subroutine: POLYMARKER
C###  Description:
C**** Draws a polymarker on iw with index IBUNDLE.
C**** D_PTS(1..3,nopts) contains the REAL*8 3D coords of each point.
C**** IF iw=4 coords are curvilinear coords, else rect. cart.
C**** If IBUNDLE is 0 the primitive will use the previously
C**** defined polymarker index.
C**** NOTE: If iw is 15 or 16 (postscript) IBUNDLE is reset to be black.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'cmiss$reference:map000.cmn'
!     Parameter List
      INTEGER IBUNDLE,iw,NTPTS
      REAL*8 D_PTS(3,*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nj,nopts

      CALL ENTERS('POLYMARKER',*9999)
      IF(NTPTS.EQ.0) GOTO 9998

      IF(DOP) THEN
        WRITE(OP_STRING,
     '    '('' iw='',I2,'' NTPTS='',I6)') iw,NTPTS
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nopts=1,NTPTS
          WRITE(OP_STRING,'('' D_PTS(nj,'',I6,''): '',3E12.3)')
     '      nopts,(D_PTS(nj,nopts),nj=1,NJT)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

      IF(iw.EQ.1) THEN !plot x against y
        IF(ibundle.ne.0) CALL GKS_SPMI(IBUNDLE,ERROR,*9999)
!       IF(NJT.EQ.2) THEN !plot x against y
          CALL GKS_PM(1,2,NTPTS,D_PTS,ERROR,*9999)
!       ELSE IF(NJT.EQ.3) THEN !plot x against z
!         CALL GKS_PM(1,3,NTPTS,D_PTS,ERROR,*9999)
!       ENDIF

      ELSE IF(iw.EQ.2) THEN !plot y against z
        IF(ibundle.ne.0) CALL GKS_SPMI(IBUNDLE,ERROR,*9999)
        CALL GKS_PM(2,3,NTPTS,D_PTS,ERROR,*9999)

      ELSE IF(iw.EQ.3) THEN !plot x against z
        IF(ibundle.ne.0) CALL GKS_SPMI(IBUNDLE,ERROR,*9999)
        CALL GKS_PM(1,3,NTPTS,D_PTS,ERROR,*9999)

!     ELSE IF(iw.EQ.3.OR.iw.EQ.67) THEN
!       IF(ibundle.ne.0) CALL PHIGS$SET_PMARKER_INDEX(IBUNDLE)
!       CALL PHIGS_POLYMARKER3(NTPTS,D_PTS,ERROR,*9999)

      ELSE IF(iw.EQ.4) THEN
        IF(ibundle.ne.0) CALL GKS_SPMI(IBUNDLE,ERROR,*9999)
        IF(PROJEC(1:11).EQ.'RECTANGULAR') THEN !points x and y
        ELSE IF(PROJEC(1:2).EQ.'XI') THEN !assume pts are in Xi coords
          DO nopts=1,NTPTS
            D_PTS(1,nopts)=-1.0D0+2.0D0*(DBLE(MXI1-1)+(D_PTS(1,nopts))
     '        /MAX_XI)
            D_PTS(2,nopts)=-1.0D0+2.0D0*(DBLE(MXI2-1)+(D_PTS(2,nopts))
     '        /MAX_XI)
          ENDDO
        ELSE  !assume points are in polar coords
!         Map4 converts mu,theta coords to map coords (-1<.<1)
          CALL MAP4(NTPTS,D_PTS,D_PTS,ERROR,*9999)
        ENDIF
        CALL GKS_PM(1,2,NTPTS,D_PTS,ERROR,*9999)

      ELSE IF(iw.EQ.5.OR.iw.EQ.6) THEN
!       frame grabber windows
c PJH 22-Jul-92 need to sort out real*8 here
c       DO nopts=1,NTPTS
c         IF(PROJEC(1:5).EQ.'IMAGE') THEN !transform to image coordinates
c           CALL FBIMGPT(iw,REAL(D_PTS(1,nopts)),REAL(D_PTS(2,nopts)),
c    '        REAL(D_PTS(3,nopts)),X(nopts),Y(nopts))
c         ELSE
c           X(nopts)=REAL(D_PTS(1,nopts))
c           Y(nopts)=REAL(D_PTS(2,nopts))
c         ENDIF
c       ENDDO
c       CALL FBPM(iw,NTPTS,X,Y)

      ELSE IF(iw.EQ.10) THEN
        IF(ibundle.ne.0) CALL GKS_SPMI(IBUNDLE,ERROR,*9999)
!       Time history plots
        DO nopts=1,NTPTS
          CALL MAP10(1,D_PTS(2,nopts))
        ENDDO
        CALL GKS_PM(1,2,NTPTS,D_PTS,ERROR,*9999)

      ELSE IF(iw.EQ.11) THEN
        IF(ibundle.ne.0) CALL GKS_SPMI(IBUNDLE,ERROR,*9999)
!       Section plots
        DO nopts=1,NTPTS
          CALL MAP10(1,D_PTS(2,nopts))
        ENDDO
        CALL GKS_PM(1,2,NTPTS,D_PTS,ERROR,*9999)

      ELSE IF(iw.EQ.12) THEN !fibre plot
        IF(ibundle.ne.0) CALL GKS_SPLI(IBUNDLE,ERROR,*9999)
        CALL GKS_PM(1,2,NTPTS,D_PTS,ERROR,*9999)

      ELSE IF(iw.EQ.13) THEN !sheet plot
        IF(ibundle.ne.0) CALL GKS_SPLI(IBUNDLE,ERROR,*9999)
        DO nopts=1,NTPTS
          D_PTS(1,nopts)=DBLE(XMIN+XMAX)-D_PTS(1,nopts)
        ENDDO
        CALL GKS_PM(2,1,NTPTS,D_PTS,ERROR,*9999) !note reversal

      ELSE IF(iw.EQ.15) THEN !GKS postscript
        IF(ibundle.ne.0) THEN
C CPB 1/11/92 removed as using colour postscript
C          IF(IBUNDLE.GT.8) THEN !reset to black for postscript
C            IBUNDLE=IBUNDLE-8
C            IF(DOP) WRITE(IOOP,'('' Reset IBUNDLE='',I2)') IBUNDLE
C          ENDIF
C CPB 1/11/12 Changed from GKS_SPLI to GKS_SPMI
          CALL GKS_SPMI(IBUNDLE,ERROR,*9999)
        ENDIF
        IF(NJT.EQ.2) THEN !plot x against y
          CALL GKS_PM(1,2,NTPTS,D_PTS,ERROR,*9999)
        ELSE IF(NJT.EQ.3) THEN !plot x against z
          CALL GKS_PM(1,3,NTPTS,D_PTS,ERROR,*9999)
        ENDIF

      ELSE IF(iw.EQ.16) THEN !PHIGS postscript
        IF(ibundle.ne.0) THEN
C cpb 17/9/92 removed as the default wstype is now colour postscript
C          IF(IBUNDLE.GT.8) THEN !reset to black for postscript
C            IBUNDLE=IBUNDLE-8
C            IF(DOP) WRITE(IOOP,'('' Reset IBUNDLE='',I2)') IBUNDLE
C          ENDIF
          CALL PHIGS$SET_PMARKER_INDEX(IBUNDLE)
        ENDIF
        CALL PHIGS_POLYMARKER3(NTPTS,D_PTS,ERROR,*9999)

      ELSE IF(iw.EQ.21.OR.iw.EQ.22.OR.iw.EQ.23.OR.iw.EQ.34.OR.iw.EQ.35
     '  .OR.iw.EQ.40.OR.iw.EQ.41.OR.iw.EQ.42.OR.iw.EQ.43.OR.iw.EQ.44
     '  .OR.iw.EQ.45.OR.iw.EQ.46.OR.iw.EQ.47
     '  .OR.iw.EQ.50.OR.iw.EQ.51
     '  .OR.iw.EQ.60.OR.iw.EQ.61.OR.iw.EQ.62.OR.iw.EQ.63.OR.iw.EQ.64
     '  .OR.iw.EQ.65.OR.iw.EQ.66.OR.iw.EQ.68.OR.iw.EQ.69) THEN !plots
        IF(ibundle.ne.0) CALL GKS_SPMI(IBUNDLE,ERROR,*9999)
        CALL GKS_PM(1,2,NTPTS,D_PTS,ERROR,*9999)

      ENDIF

 9998 CALL EXITS('POLYMARKER')
      RETURN
 9999 CALL ERRORS('POLYMARKER',ERROR)
      CALL EXITS('POLYMARKER')
      RETURN 1
      END


      SUBROUTINE PRECHOICE1(IP,iw,NOCH,NOCO,NTCH,CO,MODE,OPTION,STRING,
     '  ERROR,*)

C#### Subroutine: PRECHOICE1
C###  Description:
C**** Sets up choice menu associated with finite element windows.
C**** IP=1 is standard finite element window choices.
C**** IP=2 is subsequent choice.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
!     Parameter List
      INTEGER IP,iw,NOCH,NOCO,NTCH
      CHARACTER CO(*)*(*),ERROR*(*),MODE*(*),OPTION(*)*(*),STRING*(*)
!     Local Variables
      INTEGER INCH,INSTAT

      CALL ENTERS('PRECHOICE1',*9999)
      INCH=NTCH
      IF(IP.EQ.1) THEN
        IF(NJT.LE.2) THEN
          IF(iw.EQ.1) THEN
            CALL CHOICE('DEFINE',1,1,INSTAT,7,MODE,INCH,NTCH,NOCH,NOCO,
     '        9,CO,OPTION,STRING,0.01*XDISP,YDISP-0.4*XDISP,ERROR,*9999)
          ENDIF
        ELSE IF(NJT.EQ.3) THEN
          IF(iw.EQ.1) THEN
            CALL CHOICE('DEFINE',1,1,INSTAT,7,MODE,INCH,NTCH,NOCH,NOCO,
     '        9,CO,OPTION,STRING,0.01*DISP,0.50*DISP,ERROR,*9999)
          ELSE IF(iw.EQ.2) THEN
            CALL CHOICE('DEFINE',1,1,INSTAT,7,MODE,INCH,NTCH,NOCH,NOCO,
     '        9,CO,OPTION,STRING,1.00*DISP,0.99*DISP,ERROR,*9999)
          ELSE IF(iw.EQ.3) THEN
            CALL CHOICE('DEFINE',1,1,INSTAT,7,MODE,INCH,NTCH,NOCH,NOCO,
     '        10,CO,OPTION,STRING,1.00*DISP,0.48*DISP,ERROR,*9999)
          ELSE IF(iw.EQ.4) THEN
            CALL CHOICE('DEFINE',1,1,INSTAT,7,MODE,INCH,NTCH,NOCH,NOCO,
     '        9,CO,OPTION,STRING,0.11*DISP,0.19*DISP,ERROR,*9999)
          ENDIF
        ENDIF

      ELSE IF(IP.EQ.2) THEN
        IF(NJT.LE.2) THEN
          IF(iw.EQ.1) THEN
            CALL CHOICE('DEFINE',1,1,INSTAT,7,MODE,INCH,NTCH,NOCH,NOCO,
     '        9,CO,OPTION,STRING,0.01*XDISP,YDISP-0.45*XDISP,ERROR,
     '        *9999)
          ENDIF
        ELSE IF(NJT.EQ.3) THEN
          IF(iw.EQ.1) THEN
            CALL CHOICE('DEFINE',1,1,INSTAT,7,MODE,INCH,NTCH,NOCH,NOCO,
     '        9,CO,OPTION,STRING,0.01*DISP,0.50*DISP,ERROR,*9999)
          ELSE IF(iw.EQ.2) THEN
            CALL CHOICE('DEFINE',1,1,INSTAT,7,MODE,INCH,NTCH,NOCH,NOCO,
     '        9,CO,OPTION,STRING,1.00*DISP,0.99*DISP,ERROR,*9999)
          ELSE IF(iw.EQ.3) THEN
            CALL CHOICE('DEFINE',1,1,INSTAT,7,MODE,INCH,NTCH,NOCH,NOCO,
     '        9,CO,OPTION,STRING,0.95*DISP,0.48*DISP,ERROR,*9999)
          ELSE IF(iw.EQ.4) THEN
            CALL CHOICE('DEFINE',1,1,INSTAT,7,MODE,INCH,NTCH,NOCH,NOCO,
     '        9,CO,OPTION,STRING,0.11*DISP,0.19*DISP,ERROR,*9999)
          ENDIF
        ENDIF
      ENDIF

      CALL EXITS('PRECHOICE1')
      RETURN
 9999 CALL ERRORS('PRECHOICE1',ERROR)
      CALL EXITS('PRECHOICE1')
      RETURN 1
      END


      SUBROUTINE PRECHOICE2(iw,NOCO,CO,FILE_EXT,FILE_NAME,PARAM_TYPE,
     '  Q,QUALIFIERS,STRING,ERROR,*)

C#### Subroutine: PRECHOICE2
C###  Description:
C**** Sets up choice menu for type of qualifier.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
!     Parameter List
      INTEGER iw,NOCO
      CHARACTER CO(*)*(*),ERROR*(*),FILE_EXT*(*),FILE_NAME*(*),
     '  PARAM_TYPE*(*),Q*(*),QUALIFIERS*(*),STRING*(*)
!     Local Variables
      INTEGER IBEG,IEND,NCHAR,NOCH,NTCH,NTFILE
      CHARACTER CHOOSE*20,FULL_BRIEF*5,OPTION(22)*20,TEXT_STRING*50
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
      IF(CELEM('w',QUALIFIERS)) THEN
        NOCH=NOCH+1
        OPTION(NOCH)='write'
      ENDIF
      IF(NOCH.EQ.0) THEN
        Q='r'
      ELSE
        OPTION(NOCH+1)='FULL'
        OPTION(NOCH+2)='Exit'
        NTCH=NOCH+2
        FULL_BRIEF='BRIEF'
        CHOOSE_FULL_BRIEF=.TRUE.
        DO WHILE(CHOOSE_FULL_BRIEF)
          CALL PRECHOICE1(2,iw,NOCH,NOCO,NTCH,CO,'REQUEST',OPTION,
     '      STRING,ERROR,*9999)
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
          ELSE IF(CHOOSE(1:5).EQ.'write') THEN
            Q='w'
            CHOOSE_FULL_BRIEF=.FALSE.
          ELSE IF(CHOOSE(1:4).EQ.'Exit') THEN
            Q='e'
            CHOOSE_FULL_BRIEF=.FALSE.
          ENDIF
        ENDDO
      ENDIF

      IF(FILE_EXT(1:5).EQ.'IPFIT'.AND.(Q.EQ.'d'.OR.Q.EQ.'r')) THEN
        OPTION( 1)='geometry'
        OPTION( 2)='fibre'
        OPTION( 3)='field'
        OPTION( 4)='Fourier'
        OPTION( 5)='Exit'
        NTCH=5
        CALL PRECHOICE1(2,iw,NOCH,NOCO,NTCH,CO,'REQUEST',OPTION,STRING,
     '    ERROR,*9999)
        CHOOSE=OPTION(NOCH)
        IF(CHOOSE(1:4).ne.'Exit') THEN
          PARAM_TYPE=CHOOSE
        ELSE
          Q='e'
        ENDIF
      ENDIF

      IF(Q.EQ.'r') THEN !find & display valid files in current directory
        CALL DISPLAY_FILE(iw,NOCO,NTFILE,CO,FILE_EXT,FILE_NAME,
     '    STRING,ERROR,*9999)
        IF(FILE_NAME(1:4).EQ.'Exit') Q='e'
      ELSE IF(Q.EQ.'w') THEN !prompt for file name
        CALL STRING_TRIM(FILE00,IBEG,IEND)
        CALL GKS_STRG(56,NCHAR,'Enter file name ['
     '    //FILE00(IBEG:IEND)//']',TEXT_STRING,ERROR,*9999)
        IF(NCHAR.GT.0) THEN
          FILE_NAME=TEXT_STRING
          FILE00=FILE_NAME
        ELSE
          FILE_NAME=FILE00
        ENDIF
      ENDIF

      CALL EXITS('PRECHOICE2')
      RETURN
 9999 CALL ERRORS('PRECHOICE2',ERROR)
      CALL EXITS('PRECHOICE2')
      RETURN 1
      END


      SUBROUTINE QUIT_GRAPHICS(ERROR,*)

C#### Subroutine: QUIT_GRAPHICS
C###  Description:
C**** Close any workstations and gks or phigs

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:gks000.cmn'
      INCLUDE 'cmiss$reference:phig00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)

      CALL ENTERS('QUIT_GRAPHICS',*9999)
      IF(GKS.OR.PHIGS) THEN
        CALL CLWS(ERROR,*9999)
      ENDIF
      IF(GKS) THEN
        CALL GKS_CLKS(ERROR,*9999)
        GKS=.FALSE.
      ENDIF
      IF(PHIGS) THEN
        CALL PHIGS$CLOSE_PHIGS()
        PHIGS=.FALSE.
      ENDIF

      CALL EXITS('QUIT_GRAPHICS')
      RETURN
 9999 CALL ERRORS('QUIT_GRAPHICS',ERROR)
      CALL EXITS('QUIT_GRAPHICS')
      RETURN 1
      END


      SUBROUTINE RECALL_GRAPHICS(iw,FILE_NAME,ERROR,*)

C#### Subroutine: RECALL_GRAPHICS
C###  Description:
C**** Archives PHIGS structure or GKS segments.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:view00.cmn'
      INCLUDE 'SYS$LIBRARY:PHIGS$DEFS.FOR'
!     Parameter List
      INTEGER IBEG,IEND,iw
      CHARACTER ERROR*(*),FILE_NAME*(*)

      CALL ENTERS('RECALL_GRAPHICS',*9999)
      CALL PHIGS$SET_CONFLICT_RESOLUTION(PHIGS$K_UPDATE,
     '  PHIGS$K_UPDATE)
      CALL STRING_TRIM(FILE_NAME,IBEG,IEND)
      CALL PHIGS$OPEN_ARCHIVE(1,FILE_NAME(IBEG:IEND)//'.ARCHIVE')
      CALL PHIGS$RETRIEVE_ALL_STRUCT(1)
      CALL PHIGS$CLOSE_ARCHIVE(1)
      CALL PHIGS$POST_STRUCT(iw,ISVIEW,1.0)

      CALL EXITS('RECALL_GRAPHICS')
      RETURN
 9999 CALL ERRORS('RECALL_GRAPHICS',ERROR)
      CALL EXITS('RECALL_GRAPHICS')
      RETURN 1
      END


      SUBROUTINE ROTATE(A,*)
      DIMENSION A(*)
      CHARACTER ERROR*10
      WRITE(OP_STRING,*) '>>Cannot do under Phigs'
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
 9999 RETURN 1
      END


      SUBROUTINE SET_COLOUR_LUT(LUT,ERROR,*)

C#### Subroutine: ROTATE
C###  Description:
C**** Define colour lookup table.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'
!     Parameter List
      REAL LUT(3,0:*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j
      REAL THETA

      CALL ENTERS('SET_COLOUR_LUT',*9999)

C **  Set up colours for colour indices 0(grey), 1(black) and 2(white)
      DO j=1,3
        LUT(j,0)=0.65 !grey
        LUT(j,1)=0.0  !black
        LUT(j,2)=1.0  !white
      ENDDO

C **  Set up colours for colour indices 3 (red) to MAXCOLOURS (blue)
      DO i=3,MAXCOLOURS
        THETA=REAL(i-3)*6.0/REAL(MAXCOLOURS-3)    !THETA range
        IF(THETA.LT.2.0) THEN                     !0 - 2
          LUT(1,i)=1.0
          LUT(3,i)=0.0
          IF(THETA.LT.1.0) THEN                   !  0 - 1
            LUT(2,i)=THETA*0.75
          ELSE                                    !  1 - 2
            LUT(2,i)=0.75+(THETA-1.0)/4.0
          ENDIF
        ELSE IF(THETA.LT.4.0) THEN                !2 - 4
          LUT(1,i)=(4.0-THETA)/2.0
          LUT(2,i)=1.0
          LUT(3,i)=(THETA-2.0)/2.0
        ELSE                                      !4 - 6
          LUT(1,i)=0.0
          LUT(3,i)=1.0
          IF(THETA.LT.5.0) THEN                   !  4 - 5
            LUT(2,i)=1.0-(THETA-4.0)/4.0
          ELSE                                    !  5 - 6
            LUT(2,i)=0.75-(THETA-5.0)*0.75
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
     '  LUT,RGB_START,RGB_FINISH,ERROR,*)

C#### Subroutine: SET_COLOUR_LUT_RANGE
C###  Description:
C**** Define colour lookup table for a given range INDEX_START to INDEX_FINISH
C**** as percentages (0-100), and for a given RGB range.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'
!     Parameter List
      INTEGER INDEX_START,INDEX_FINISH
      REAL LUT(3,0:*)
      CHARACTER ERROR*(*),RGB_FINISH*(*),RGB_START*(*)
!     Local Variables
      INTEGER i,IFINISH,ISTART,j
      REAL RGBFCOL(3),RGBSCOL(3),RGBSTEP(3),THETA
      REAL*8 RFROMC
      CHARACTER CUPPER*10,RGBF*10,RGBS*10

      CALL ENTERS('SET_COLOUR_LUT_RANGE',*9999)

      CALL ASSERT(INDEX_START.GE.0,'>> Invalid Range (Start < 0)',
     '  ERROR,*9999)
      CALL ASSERT(INDEX_FINISH.LE.100,'>> Invalid Range (Finish > 100)',
     '  ERROR,*9999)
      CALL ASSERT(INDEX_FINISH.GE.INDEX_START,'>> Invalid Range',
     '  ERROR,*9999)

      RGBS=CUPPER(RGB_START)
      RGBF=CUPPER(RGB_FINISH)

      IF(RGBS(1:5).EQ.'BLACK'.OR.RGBS(1:6).EQ.'000000') THEN
        RGBSCOL(1)=0.0
        RGBSCOL(2)=0.0
        RGBSCOL(3)=0.0
      ELSE IF(RGBS(1:5).EQ.'WHITE'.OR.RGBS(1:6).EQ.'999999') THEN
        RGBSCOL(1)=1.0
        RGBSCOL(2)=1.0
        RGBSCOL(3)=1.0
      ELSE IF(RGBS(1:3).EQ.'RED'.OR.RGBS(1:6).EQ.'990000') THEN
        RGBSCOL(1)=1.0
        RGBSCOL(2)=0.0
        RGBSCOL(3)=0.0
      ELSE IF(RGBS(1:5).EQ.'GREEN'.OR.RGBS(1:6).EQ.'009900') THEN
        RGBSCOL(1)=0.0
        RGBSCOL(2)=1.0
        RGBSCOL(3)=0.0
      ELSE IF(RGBS(1:4).EQ.'BLUE'.OR.RGBS(1:6).EQ.'000099') THEN
        RGBSCOL(1)=0.0
        RGBSCOL(2)=0.0
        RGBSCOL(3)=1.0
      ELSE IF(RGBS(1:4).EQ.'CYAN'.OR.RGBS(1:6).EQ.'009999') THEN
        RGBSCOL(1)=0.0
        RGBSCOL(2)=1.0
        RGBSCOL(3)=1.0
      ELSE IF(RGBS(1:6).EQ.'YELLOW'.OR.RGBS(1:6).EQ.'999900') THEN
        RGBSCOL(1)=1.0
        RGBSCOL(2)=1.0
        RGBSCOL(3)=0.0
      ELSE IF(RGBS(1:6).EQ.'PURPLE'.OR.RGBS(1:6).EQ.'990099') THEN
        RGBSCOL(1)=1.0
        RGBSCOL(2)=0.0
        RGBSCOL(3)=1.0
      ELSE
        DO I=1,5,2
          RGBSCOL((I+1)/2)=REAL(RFROMC(RGBS(I:I+1))/99.0)
        ENDDO
      ENDIF

      IF(RGBF(1:5).EQ.'BLACK'.OR.RGBF(1:6).EQ.'000000') THEN
        RGBFCOL(1)=0.0
        RGBFCOL(2)=0.0
        RGBFCOL(3)=0.0
      ELSE IF(RGBF(1:5).EQ.'WHITE'.OR.RGBF(1:6).EQ.'999999') THEN
        RGBFCOL(1)=1.0
        RGBFCOL(2)=1.0
        RGBFCOL(3)=1.0
      ELSE IF(RGBF(1:3).EQ.'RED'.OR.RGBF(1:6).EQ.'990000') THEN
        RGBFCOL(1)=1.0
        RGBFCOL(2)=0.0
        RGBFCOL(3)=0.0
      ELSE IF(RGBF(1:5).EQ.'GREEN'.OR.RGBF(1:6).EQ.'009900') THEN
        RGBFCOL(1)=0.0
        RGBFCOL(2)=1.0
        RGBFCOL(3)=0.0
      ELSE IF(RGBF(1:4).EQ.'BLUE'.OR.RGBF(1:6).EQ.'000099') THEN
        RGBFCOL(1)=0.0
        RGBFCOL(2)=0.0
        RGBFCOL(3)=1.0
      ELSE IF(RGBF(1:4).EQ.'CYAN'.OR.RGBF(1:6).EQ.'009999') THEN
        RGBFCOL(1)=0.0
        RGBFCOL(2)=1.0
        RGBFCOL(3)=1.0
      ELSE IF(RGBF(1:6).EQ.'YELLOW'.OR.RGBF(1:6).EQ.'999900') THEN
        RGBFCOL(1)=1.0
        RGBFCOL(2)=1.0
        RGBFCOL(3)=0.0
      ELSE IF(RGBF(1:6).EQ.'PURPLE'.OR.RGBF(1:6).EQ.'990099') THEN
        RGBFCOL(1)=1.0
        RGBFCOL(2)=0.0
        RGBFCOL(3)=1.0
      ELSE
        DO I=1,5,2
          RGBFCOL((I+1)/2)=REAL(RFROMC(RGBF(I:I+1))/99.0)
        ENDDO
      ENDIF

      DO I=1,3
        RGBSTEP(I)=RGBFCOL(I)-RGBSCOL(I)
      ENDDO

      ISTART=NINT(DBLE(INDEX_START)/100.0*REAL(MAXCOLOURS-3)+3.0)
      IFINISH=NINT(DBLE(INDEX_FINISH)/100.0*REAL(MAXCOLOURS-3)+3.0)
      IF(ISTART.EQ.IFINISH) THEN
        DO I=1,3
          LUT(I,ISTART)=RGBSCOL(I)
        ENDDO
      ELSE
        DO J=ISTART,IFINISH
          THETA=REAL(J-ISTART)/REAL(IFINISH-ISTART)
          DO I=1,3
            LUT(I,J)=RGBSCOL(I)+THETA*RGBSTEP(I)
          ENDDO
        ENDDO
      ENDIF

      CALL EXITS('SET_COLOUR_LUT_RANGE')
      RETURN
 9999 CALL ERRORS('SET_COLOUR_LUT_RANGE',ERROR)
      CALL EXITS('SET_COLOUR_LUT_RANGE')
      RETURN 1
      END


      SUBROUTINE SET_COLOUR_ONE(iw,ilut,COLOUR,ERROR,*)

C#### Subroutine: SET_COLOUR_ONE
C###  Description:
C**** Sets colour representation for index ICOL on workstation iw.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
!     Parameter List
      INTEGER ilut,iw
      REAL COLOUR(*)
      CHARACTER ERROR*(*)

      CALL ENTERS('SET_COLOUR_ONE',*9999)
      IF(IWKT(iw).EQ.1) THEN      !GKS
        CALL GKS_SCR(iw,ilut,COLOUR(1),COLOUR(2),COLOUR(3),ERROR,*9999)
      ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS
        CALL PHIGS$SET_COLOUR_REP(iw,ilut,COLOUR)
      ENDIF

      CALL EXITS('SET_COLOUR_ONE')
      RETURN
 9999 CALL ERRORS('SET_COLOUR_ONE',ERROR)
      CALL EXITS('SET_COLOUR_ONE')
      RETURN 1
      END


      SUBROUTINE SET_COLOUR_REP(iw,ILUT_MIN,ILUT_MAX,ERROR,*)

C#### Subroutine: SET_COLOUR_REP
C###  Description:
C**** Sets colour representation for colour indices on workstation iw.
C**** ilut=0   is grey
C**** ilut=1   is black
C**** ilut=2   is white
C**** ilut=ILUT_MIN..ILUT_MAX are red .. blue,
C**** ILUT_MIN is 3 and ILUT_MAX is MAXCOLOURS unless
C****   redefined by 'change colour'

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
!     Parameter List
      INTEGER ILUT_MAX,ILUT_MIN,iw
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,ICOL,ilut

      CALL ENTERS('SET_COLOUR_REP',*9999)

      IF(DOP) THEN
        WRITE(OP_STRING,
     '    '('' ILUT_MIN='',I4,'' ILUT_MAX='',I4)') ILUT_MIN,
     '    ILUT_MAX
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      DO ilut=0,1            !Set grey and black on all w/s
        ICOL=ilut
        IF(DOP) THEN
          WRITE(OP_STRING,
     '      '('' ILUT='',I4,'' ICOL='',I4,'' COLOUR_LUT:''3F6.3)')
     '      ilut,ICOL,(COLOUR_LUT(i,ICOL),i=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        IF(IWKT(iw).EQ.1) THEN      !GKS
          CALL GKS_SCR(iw,ilut,COLOUR_LUT(1,ICOL),COLOUR_LUT(2,ICOL),
     '      COLOUR_LUT(3,ICOL),ERROR,*9999)
        ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS
          CALL PHIGS$SET_COLOUR_REP(iw,ilut,COLOUR_LUT(1,ICOL))
        ENDIF
      ENDDO

      IF(COLOUR_WS) THEN
        ilut=2                !Set white only on colour w/s
        ICOL=ilut
        IF(DOP) THEN
          WRITE(OP_STRING,'('' ILUT='',I4,'' ICOL='',I4,'
     '      //''' COLOUR_LUT:''3F6.3)')
     '      ilut,ICOL,(COLOUR_LUT(I,ICOL),I=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        IF(IWKT(iw).EQ.1) THEN      !GKS
          CALL GKS_SCR(iw,ilut,COLOUR_LUT(1,ICOL),COLOUR_LUT(2,ICOL),
     '      COLOUR_LUT(3,ICOL),ERROR,*9999)
        ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS
          CALL PHIGS$SET_COLOUR_REP(iw,ilut,COLOUR_LUT(1,ICOL))
        ENDIF
        DO ilut=3,ILUT_MIN-1 !puts all indices up to ILUT_MIN as red
          ICOL=3
          IF(DOP) THEN
            WRITE(OP_STRING,
     '        '('' ILUT='',I4,'' ICOL='',I4,'' COLOUR_LUT:''3F6.3)')
     '        ilut,ICOL,(COLOUR_LUT(i,ICOL),i=1,3)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(IWKT(iw).EQ.1) THEN      !GKS
            CALL GKS_SCR(iw,ilut,COLOUR_LUT(1,ICOL),COLOUR_LUT(2,ICOL),
     '        COLOUR_LUT(3,ICOL),ERROR,*9999)
          ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS
            CALL PHIGS$SET_COLOUR_REP(iw,ilut,COLOUR_LUT(1,ICOL))
          ENDIF
        ENDDO
        DO ilut=ILUT_MIN,ILUT_MAX
          ICOL=3+INT(DBLE(ilut-ILUT_MIN)/DBLE(ILUT_MAX-ILUT_MIN)*
     '      (MAXCOLOURS-3))
          IF(DOP) THEN
            WRITE(OP_STRING,
     '        '('' ILUT='',I4,'' ICOL='',I4,'' COLOUR_LUT:''3F6.3)')
     '        ilut,ICOL,(COLOUR_LUT(i,ICOL),i=1,3)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(IWKT(iw).EQ.1) THEN      !GKS
            CALL GKS_SCR(iw,ilut,COLOUR_LUT(1,ICOL),COLOUR_LUT(2,ICOL),
     '        COLOUR_LUT(3,ICOL),ERROR,*9999)
          ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS
            CALL PHIGS$SET_COLOUR_REP(iw,ilut,COLOUR_LUT(1,ICOL))
          ENDIF
        ENDDO
        DO ilut=ILUT_MAX+1,MAXCOLOURS
                    !puts all indices above ILUT_MAX as blue
          ICOL=MAXCOLOURS
          IF(DOP) THEN
            WRITE(OP_STRING,
     '        '('' ILUT='',I4,'' ICOL='',I4,'' COLOUR_LUT:''3F6.3)')
     '        ilut,ICOL,(COLOUR_LUT(i,ICOL),i=1,3)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(IWKT(iw).EQ.1) THEN      !GKS
            CALL GKS_SCR(iw,ilut,COLOUR_LUT(1,ICOL),COLOUR_LUT(2,ICOL),
     '        COLOUR_LUT(3,ICOL),ERROR,*9999)
          ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS
            CALL PHIGS$SET_COLOUR_REP(iw,ilut,COLOUR_LUT(1,ICOL))
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('SET_COLOUR_REP')
      RETURN
 9999 CALL ERRORS('SET_COLOUR_REP',ERROR)
      CALL EXITS('SET_COLOUR_REP')
      RETURN 1
      END


      SUBROUTINE SET_FILL_REP(iw,INDEX,STYLE,ISTYLE,icolour,ERROR,*)

C#### Subroutine: SET_FILL_REP
C###  Description:
C**** Resets fill area representation

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
!     Parameter List
      INTEGER icolour,index,ISTYLE,iw
      CHARACTER STYLE*(*),ERROR*(*)

      CALL ENTERS('SET_FILL_REP',*9999)
      IF(IWKT(iw).EQ.1) THEN      !GKS
        IF(STYLE(1:7).EQ.'PATTERN') THEN
          CALL GKS_SFAR(iw,index,GPATTR,ISTYLE,icolour,ERROR,*9999)
        ELSE IF(STYLE(1:5).EQ.'SOLID') THEN
          CALL GKS_SFAR(iw,index,GSOLID,ISTYLE,icolour,ERROR,*9999)
        ENDIF
      ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS
      ENDIF

      CALL EXITS('SET_FILL_REP')
      RETURN
 9999 CALL ERRORS('SET_FILL_REP',ERROR)
      CALL EXITS('SET_FILL_REP')
      RETURN 1
      END


      SUBROUTINE SET_FILL_AREA_REP(iw,ERROR,*)

C#### Subroutine: SET_FILL_AREA_REP
C###  Description:
C**** Set FILL_AREA representation on workstation iw.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.FOR'
      INCLUDE 'SYS$LIBRARY:PHIGS$DEFS.FOR'
!     Parameter List
      INTEGER iw
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER icolour,ilut,index,INDEX_FILL_AREA,ipat

      CALL ENTERS('SET_FILL_AREA_REP',*9999)

!     Set fill indices 1 & 16 to be black and white
      CALL GKS_SFAR(iw, 1,GSOLID,14,1,ERROR,*9999) !index 1 (black)
      CALL GKS_SFAR(iw,16,GSOLID,14,0,ERROR,*9999) !index 16(white)
!     Set fill indices 2..15 to be patterns giving 16 shades of grey
      IF(iw.ne.15) THEN
        DO ipat=2,15
          CALL GKS_SFAR(iw,ipat,GPATTR,12+(17-ipat),1,ERROR,*9999)
        ENDDO
      ENDIF
      IF(COLOUR_WS) THEN !for colour use predefined colours
        DO icolour=1,MAXCOLOURS-16
          ilut=INT(3.+DBLE(icolour-1)*(MAXCOLOURS-3)/(MAXCOLOURS-17))
          index=INDEX_FILL_AREA(icolour,'SOLID',' ',' ')
          CALL GKS_SFAR(iw,index,GSOLID,1,ilut,ERROR,*9999)
        ENDDO
      ENDIF

      CALL EXITS('SET_FILL_AREA_REP')
      RETURN
 9999 CALL ERRORS('SET_FILL_AREA_REP',ERROR)
      CALL EXITS('SET_FILL_AREA_REP')
      RETURN 1
      END


      SUBROUTINE SET_POLYLINE_REP(iw,ERROR,*)

C#### Subroutine: SET_POLYLINE_REP
C###  Description:
C**** Set POLYLINE representation on workstation iw.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.FOR'
      INCLUDE 'SYS$LIBRARY:PHIGS$DEFS.FOR'
!     Parameter List
      INTEGER iw
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,icolour,ilut,index,INDEX_POLYLINE,LINE_TYPE1,LINE_TYPE2,
     '  LINE_TYPE3,LINE_TYPE4,LUT_COLOUR_INDEX(8)
      REAL LINE_WIDTH1,LINE_WIDTH2
      CHARACTER LINE_COLOUR_NAME(8)*6

      DATA LINE_COLOUR_NAME /'RED','GREEN','BLUE','CYAN','YELLOW',
     '  'WHITE','LTBLUE','GREY'/
      LUT_COLOUR_INDEX(1) = 3
      LUT_COLOUR_INDEX(2) = (3+(MAXCOLOURS-3)/2)
      LUT_COLOUR_INDEX(3) = MAXCOLOURS
      LUT_COLOUR_INDEX(4) = 3+2*(MAXCOLOURS-3)/3
      LUT_COLOUR_INDEX(5) = 3+(MAXCOLOURS-3)/3
      LUT_COLOUR_INDEX(6) = 2
      LUT_COLOUR_INDEX(7) = 3+(MAXCOLOURS-3)*3/4
      LUT_COLOUR_INDEX(8) = 0

      CALL ENTERS('SET_POLYLINE_REP',*9999)
      LINE_TYPE1=GLSOLI
      LINE_TYPE2=GLDOT
      LINE_TYPE3=GLDASH
      LINE_TYPE4=GLDASD
      LINE_WIDTH1=1.0
      LINE_WIDTH2=4.0

      IF(IWKT(iw).EQ.1) THEN      !GKS
        index=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLACK')
        CALL GKS_SPLR(iw,index,LINE_TYPE1,LINE_WIDTH1,1,ERROR,*9999)
        index=INDEX_POLYLINE(0,'DOTTED','WIDTH1','BLACK')
        CALL GKS_SPLR(iw,index,LINE_TYPE2,LINE_WIDTH1,1,ERROR,*9999)
        index=INDEX_POLYLINE(0,'DASHED','WIDTH1','BLACK')
        CALL GKS_SPLR(iw,index,LINE_TYPE3,LINE_WIDTH1,1,ERROR,*9999)
        index=INDEX_POLYLINE(0,'DOT-DASH','WIDTH1','BLACK')
        CALL GKS_SPLR(iw,index,LINE_TYPE4,LINE_WIDTH1,1,ERROR,*9999)

        index=INDEX_POLYLINE(0,'SOLID','WIDTH2','BLACK')
        CALL GKS_SPLR(iw,index,LINE_TYPE1,LINE_WIDTH2,1,ERROR,*9999)
        index=INDEX_POLYLINE(0,'DOTTED','WIDTH2','BLACK')
        CALL GKS_SPLR(iw,index,LINE_TYPE2,LINE_WIDTH2,1,ERROR,*9999)
        index=INDEX_POLYLINE(0,'DASHED','WIDTH2','BLACK')
        CALL GKS_SPLR(iw,index,LINE_TYPE3,LINE_WIDTH2,1,ERROR,*9999)
        index=INDEX_POLYLINE(0,'DOT-DASH','WIDTH2','BLACK')
        CALL GKS_SPLR(iw,index,LINE_TYPE4,LINE_WIDTH2,1,ERROR,*9999)

        IF(COLOUR_WS) THEN
          DO i=1,8
            index=INDEX_POLYLINE(0,'SOLID','WIDTH1',LINE_COLOUR_NAME(i))
            CALL GKS_SPLR(iw,index,LINE_TYPE1,LINE_WIDTH1,
     '        LUT_COLOUR_INDEX(i),ERROR,*9999)
            index=INDEX_POLYLINE(0,'DOTTED','WIDTH1',
     '        LINE_COLOUR_NAME(i))
            CALL GKS_SPLR(iw,index,LINE_TYPE2,LINE_WIDTH1,
     '        LUT_COLOUR_INDEX(i),ERROR,*9999)
            index=INDEX_POLYLINE(0,'DASHED','WIDTH1',
     '        LINE_COLOUR_NAME(i))
            CALL GKS_SPLR(iw,index,LINE_TYPE3,LINE_WIDTH1,
     '        LUT_COLOUR_INDEX(i),ERROR,*9999)
            index=INDEX_POLYLINE(0,'DOT-DASH','WIDTH1',
     '        LINE_COLOUR_NAME(i))
            CALL GKS_SPLR(iw,index,LINE_TYPE4,LINE_WIDTH1,
     '        LUT_COLOUR_INDEX(i),ERROR,*9999)
          ENDDO

          DO icolour=1,MAXCOLOURS-40
            ilut=INT(3.+DBLE(icolour-1)*(MAXCOLOURS-3)/(MAXCOLOURS-41))
            index=INDEX_POLYLINE(icolour,'SOLID','WIDTH2',' ')
            CALL GKS_SPLR(iw,index,LINE_TYPE1,LINE_WIDTH2,ilut,
     '        ERROR,*9999)
          ENDDO
        ENDIF

      ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS
        index=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLACK')
        CALL PHIGS$SET_PLINE_REP(iw,index,LINE_TYPE1,LINE_WIDTH1,1)
        index=INDEX_POLYLINE(0,'DOTTED','WIDTH1','BLACK')
        CALL PHIGS$SET_PLINE_REP(iw,index,LINE_TYPE2,LINE_WIDTH1,1)
        index=INDEX_POLYLINE(0,'DASHED','WIDTH1','BLACK')
        CALL PHIGS$SET_PLINE_REP(iw,index,LINE_TYPE3,LINE_WIDTH1,1)
        index=INDEX_POLYLINE(0,'DOT-DASH','WIDTH1','BLACK')
        CALL PHIGS$SET_PLINE_REP(iw,index,LINE_TYPE4,LINE_WIDTH1,1)

        index=INDEX_POLYLINE(0,'SOLID','WIDTH2','BLACK')
        CALL PHIGS$SET_PLINE_REP(iw,index,LINE_TYPE1,LINE_WIDTH2,1)
        index=INDEX_POLYLINE(0,'DOTTED','WIDTH2','BLACK')
        CALL PHIGS$SET_PLINE_REP(iw,index,LINE_TYPE2,LINE_WIDTH2,1)
        index=INDEX_POLYLINE(0,'DASHED','WIDTH2','BLACK')
        CALL PHIGS$SET_PLINE_REP(iw,index,LINE_TYPE3,LINE_WIDTH2,1)
        index=INDEX_POLYLINE(0,'DOT-DASH','WIDTH2','BLACK')
        CALL PHIGS$SET_PLINE_REP(iw,index,LINE_TYPE4,LINE_WIDTH2,1)

        IF(COLOUR_WS) THEN
          DO i=1,8
            index=INDEX_POLYLINE(0,'SOLID','WIDTH1',LINE_COLOUR_NAME(i))
            CALL PHIGS$SET_PLINE_REP(iw,index,LINE_TYPE1,LINE_WIDTH1,
     '        LUT_COLOUR_INDEX(i))
            index=INDEX_POLYLINE(0,'DOTTED','WIDTH1',
     '        LINE_COLOUR_NAME(i))
            CALL PHIGS$SET_PLINE_REP(iw,index,LINE_TYPE2,LINE_WIDTH1,
     '        LUT_COLOUR_INDEX(i))
            index=INDEX_POLYLINE(0,'DASHED','WIDTH1',
     '        LINE_COLOUR_NAME(i))
            CALL PHIGS$SET_PLINE_REP(iw,index,LINE_TYPE3,LINE_WIDTH1,
     '        LUT_COLOUR_INDEX(i))
            index=INDEX_POLYLINE(0,'DOT-DASH','WIDTH1',
     '        LINE_COLOUR_NAME(i))
            CALL PHIGS$SET_PLINE_REP(iw,index,LINE_TYPE4,LINE_WIDTH1,
     '        LUT_COLOUR_INDEX(i))
          ENDDO

          DO icolour=1,MAXCOLOURS-40
            ilut=INT(3.+DBLE(icolour-1)*(MAXCOLOURS-3)/(MAXCOLOURS-41))
            index=INDEX_POLYLINE(icolour,'SOLID','WIDTH2',' ')
            CALL PHIGS$SET_PLINE_REP(iw,index,LINE_TYPE1,LINE_WIDTH2,
     '        ilut)
          ENDDO
        ENDIF
      ENDIF

      CALL EXITS('SET_POLYLINE_REP')
      RETURN
 9999 CALL ERRORS('SET_POLYLINE_REP',ERROR)
      CALL EXITS('SET_POLYLINE_REP')
      RETURN 1
      END


      SUBROUTINE SET_POLYMARKER_REP(iw,ERROR,*)

C#### Subroutine: SET_POLYMARKER_REP
C###  Description:
C**** Set POLYMARKER representation on workstation iw.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.FOR'
      INCLUDE 'SYS$LIBRARY:PHIGS$DEFS.FOR'
!     Parameter List
      INTEGER iw
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,icolour,ilut,index,INDEX_POLYMARKER,LUT_COLOUR_INDEX(8),
     '  MARKER_TYPE1,MARKER_TYPE2,MARKER_TYPE3,MARKER_TYPE4
      REAL MARKER_SIZE1,MARKER_SIZE2
      CHARACTER MARKER_COLOUR_NAME(8)*6

      DATA MARKER_COLOUR_NAME /'RED','GREEN','BLUE','CYAN','YELLOW',
     '  'WHITE','LTBLUE','GREY'/
      LUT_COLOUR_INDEX(1) = 3
      LUT_COLOUR_INDEX(2) = (3+(MAXCOLOURS-3)/2)
      LUT_COLOUR_INDEX(3) = MAXCOLOURS
      LUT_COLOUR_INDEX(4) = 3+2*(MAXCOLOURS-3)/3
      LUT_COLOUR_INDEX(5) = 3+(MAXCOLOURS-3)/3
      LUT_COLOUR_INDEX(6) = 2
      LUT_COLOUR_INDEX(7) = 3+(MAXCOLOURS-3)*3/4
      LUT_COLOUR_INDEX(8) = 0

      CALL ENTERS('SET_POLYMARKER_REP',*9999)

      MARKER_SIZE1=1.0
      MARKER_SIZE2=2.0
      MARKER_TYPE1=GPLUS
      MARKER_TYPE2=GAST
      MARKER_TYPE3=GOMARK
      MARKER_TYPE4=GPOINT

      IF(IWKT(iw).EQ.1) THEN      !GKS
        index=INDEX_POLYMARKER(0,'PLUS','SIZE1','BLACK')
        CALL GKS_SPMR(iw,index,MARKER_TYPE1,MARKER_SIZE1,1,ERROR,*9999)
        index=INDEX_POLYMARKER(0,'ASTERISK','SIZE1','BLACK')
        CALL GKS_SPMR(iw,index,MARKER_TYPE2,MARKER_SIZE1,1,ERROR,*9999)
        index=INDEX_POLYMARKER(0,'CIRCLE','SIZE1','BLACK')
        CALL GKS_SPMR(iw,index,MARKER_TYPE3,MARKER_SIZE1,1,ERROR,*9999)
        index=INDEX_POLYMARKER(0,'POINT','SIZE1','BLACK')
        CALL GKS_SPMR(iw,index,MARKER_TYPE4,MARKER_SIZE1,1,ERROR,*9999)

        index=INDEX_POLYMARKER(0,'PLUS','SIZE2','BLACK')
        CALL GKS_SPMR(iw,index,MARKER_TYPE1,MARKER_SIZE2,1,ERROR,*9999)
        index=INDEX_POLYMARKER(0,'ASTERISK','SIZE2','BLACK')
        CALL GKS_SPMR(iw,index,MARKER_TYPE2,MARKER_SIZE2,1,ERROR,*9999)
        index=INDEX_POLYMARKER(0,'CIRCLE','SIZE2','BLACK')
        CALL GKS_SPMR(iw,index,MARKER_TYPE3,MARKER_SIZE2,1,ERROR,*9999)
        index=INDEX_POLYMARKER(0,'POINT','SIZE2','BLACK')
        CALL GKS_SPMR(iw,index,MARKER_TYPE4,MARKER_SIZE2,1,ERROR,*9999)

        IF(COLOUR_WS) THEN
          DO i=1,8
            index=INDEX_POLYMARKER(0,'PLUS','SIZE1',
     '        MARKER_COLOUR_NAME(i))
            CALL GKS_SPMR(iw,index,MARKER_TYPE1,MARKER_SIZE1,
     '        LUT_COLOUR_INDEX(i),ERROR,*9999)
            index=INDEX_POLYMARKER(0,'ASTERISK','SIZE1',
     '        MARKER_COLOUR_NAME(i))
            CALL GKS_SPMR(iw,index,MARKER_TYPE2,MARKER_SIZE1,
     '        LUT_COLOUR_INDEX(i),ERROR,*9999)
            index=INDEX_POLYMARKER(0,'CIRCLE','SIZE1',
     '        MARKER_COLOUR_NAME(i))
            CALL GKS_SPMR(iw,index,MARKER_TYPE3,MARKER_SIZE1,
     '        LUT_COLOUR_INDEX(i),ERROR,*9999)
            index=INDEX_POLYMARKER(0,'POINT','SIZE1',
     '        MARKER_COLOUR_NAME(i))
            CALL GKS_SPMR(iw,index,MARKER_TYPE4,MARKER_SIZE1,
     '        LUT_COLOUR_INDEX(i),ERROR,*9999)
          ENDDO

          DO icolour=1,MAXCOLOURS-40
            ilut=INT(3.+DBLE(icolour-1)*(MAXCOLOURS-3)/(MAXCOLOURS-41))
            index=INDEX_POLYMARKER(icolour,'PLUS','SIZE2',' ')
            CALL GKS_SPMR(iw,index,MARKER_TYPE1,MARKER_SIZE2,ilut,
     '        ERROR,*9999)
          ENDDO
        ENDIF

      ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS
        index=INDEX_POLYMARKER(0,'PLUS','SIZE1','BLACK')
        CALL PHIGS$SET_PMARKER_REP(iw,index,MARKER_TYPE1,MARKER_SIZE1,1)
        index=INDEX_POLYMARKER(0,'ASTERISK','SIZE1','BLACK')
        CALL PHIGS$SET_PMARKER_REP(iw,index,MARKER_TYPE2,MARKER_SIZE1,1)
        index=INDEX_POLYMARKER(0,'CIRCLE','SIZE1','BLACK')
        CALL PHIGS$SET_PMARKER_REP(iw,index,MARKER_TYPE3,MARKER_SIZE1,1)
        index=INDEX_POLYMARKER(0,'POINT','SIZE1','BLACK')
        CALL PHIGS$SET_PMARKER_REP(iw,index,MARKER_TYPE4,MARKER_SIZE1,1)

        index=INDEX_POLYMARKER(0,'PLUS','SIZE2','BLACK')
        CALL PHIGS$SET_PMARKER_REP(iw,index,MARKER_TYPE1,MARKER_SIZE2,1)
        index=INDEX_POLYMARKER(0,'ASTERISK','SIZE2','BLACK')
        CALL PHIGS$SET_PMARKER_REP(iw,index,MARKER_TYPE2,MARKER_SIZE2,1)
        index=INDEX_POLYMARKER(0,'CIRCLE','SIZE2','BLACK')
        CALL PHIGS$SET_PMARKER_REP(iw,index,MARKER_TYPE3,MARKER_SIZE2,1)
        index=INDEX_POLYMARKER(0,'POINT','SIZE2','BLACK')
        CALL PHIGS$SET_PMARKER_REP(iw,index,MARKER_TYPE4,MARKER_SIZE2,1)

        IF(COLOUR_WS) THEN
          DO i=1,8
            index=INDEX_POLYMARKER(0,'PLUS','SIZE1',
     '        MARKER_COLOUR_NAME(i))
            CALL PHIGS$SET_PMARKER_REP(iw,index,MARKER_TYPE1,
     '        MARKER_SIZE1,LUT_COLOUR_INDEX(i))
            index=INDEX_POLYMARKER(0,'ASTERISK','SIZE1',
     '        MARKER_COLOUR_NAME(i))
            CALL PHIGS$SET_PMARKER_REP(iw,index,MARKER_TYPE2,
     '        MARKER_SIZE1,LUT_COLOUR_INDEX(I))
            index=INDEX_POLYMARKER(0,'CIRCLE','SIZE1',
     '        MARKER_COLOUR_NAME(i))
            CALL PHIGS$SET_PMARKER_REP(iw,index,MARKER_TYPE3,
     '        MARKER_SIZE1,LUT_COLOUR_INDEX(I))
            index=INDEX_POLYMARKER(0,'POINT','SIZE1',
     '        MARKER_COLOUR_NAME(i))
            CALL PHIGS$SET_PMARKER_REP(iw,index,MARKER_TYPE4,
     '        MARKER_SIZE1,LUT_COLOUR_INDEX(i))
          ENDDO

          DO icolour=1,216
            ilut=INT(3.+DBLE(icolour-1)*246./215.)
            index=INDEX_POLYMARKER(icolour,'PLUS','SIZE2',' ')
            CALL PHIGS$SET_PMARKER_REP(iw,index,MARKER_TYPE1,
     '        MARKER_SIZE2,ilut)
          ENDDO
        ENDIF
      ENDIF

      CALL EXITS('SET_POLYMARKER_REP')
      RETURN
 9999 CALL ERRORS('SET_POLYMARKER_REP',ERROR)
      CALL EXITS('SET_POLYMARKER_REP')
      RETURN 1
      END


      SUBROUTINE SET_TEXT_REP(iw,ERROR,*)

C#### Subroutine: SET_TEXT_REP
C###  Description:
C**** Set TEXT representation on workstation iw.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.FOR'
      INCLUDE 'SYS$LIBRARY:PHIGS$DEFS.FOR'
!     Parameter List
      INTEGER iw
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,icolour,ilut,index,INDEX_TEXT,LUT_COLOUR_INDEX(8),
     '  TEXT_FONT1,TEXT_FONT2,TEXT_FONT3,TEXT_FONT4,
     '  TEXT_PRECISION1,TEXT_PRECISION2,TEXT_PRECISION3
      REAL TEXT_SPACING,TEXT_WIDTH1,TEXT_WIDTH2
      CHARACTER TEXT_COLOUR_NAME(8)*6

      DATA TEXT_COLOUR_NAME /'RED','GREEN','BLUE','CYAN','YELLOW',
     '  'WHITE','LTBLUE','GREY'/
      LUT_COLOUR_INDEX(1) = 3
      LUT_COLOUR_INDEX(2) = (3+(MAXCOLOURS-3)/2)
      LUT_COLOUR_INDEX(3) = MAXCOLOURS
      LUT_COLOUR_INDEX(4) = 3+2*(MAXCOLOURS-3)/3
      LUT_COLOUR_INDEX(5) = 3+(MAXCOLOURS-3)/3
      LUT_COLOUR_INDEX(6) = 2
      LUT_COLOUR_INDEX(7) = 3+(MAXCOLOURS-3)*3/4
      LUT_COLOUR_INDEX(8) = 0

      CALL ENTERS('SET_TEXT_REP',*9999)

      TEXT_WIDTH1=0.5
      TEXT_WIDTH2=1.0
      TEXT_FONT1=1      !Standard font in string precision
      TEXT_FONT2=-9     !Serif        font in stroke precision
      TEXT_FONT3=-11    !Serif italic   "   "    "       "
      TEXT_FONT4=-12    !Sans Serif     "   "    "       "
      TEXT_PRECISION1=GSTRP    !=0
      TEXT_PRECISION2=GCHARP   !=1
      TEXT_PRECISION3=GSTRKP   !=2
      TEXT_SPACING=0.0

      IF(IWKT(iw).EQ.1) THEN      !GKS
        index=INDEX_TEXT(0,'WIDTH1','FONT1','BLACK')
        CALL GKS_STXR(iw,index,TEXT_FONT1,TEXT_PRECISION1,TEXT_WIDTH1,
     '    TEXT_SPACING,1,ERROR,*9999)
        index=INDEX_TEXT(0,'WIDTH1','FONT2','BLACK')
        CALL GKS_STXR(iw,index,TEXT_FONT2,TEXT_PRECISION3,TEXT_WIDTH1,
     '    TEXT_SPACING,1,ERROR,*9999)
        index=INDEX_TEXT(0,'WIDTH1','FONT3','BLACK')
        CALL GKS_STXR(iw,index,TEXT_FONT3,TEXT_PRECISION3,TEXT_WIDTH1,
     '    TEXT_SPACING,1,ERROR,*9999)
        index=INDEX_TEXT(0,'WIDTH1','FONT4','BLACK')
        CALL GKS_STXR(iw,index,TEXT_FONT4,TEXT_PRECISION3,TEXT_WIDTH1,
     '    TEXT_SPACING,1,ERROR,*9999)

        index=INDEX_TEXT(0,'WIDTH2','FONT1','BLACK')
        CALL GKS_STXR(iw,index,TEXT_FONT1,TEXT_PRECISION1,TEXT_WIDTH2,
     '    TEXT_SPACING,1,ERROR,*9999)
        index=INDEX_TEXT(0,'WIDTH2','FONT2','BLACK')
        CALL GKS_STXR(iw,index,TEXT_FONT2,TEXT_PRECISION3,TEXT_WIDTH2,
     '    TEXT_SPACING,1,ERROR,*9999)
        index=INDEX_TEXT(0,'WIDTH2','FONT3','BLACK')
        CALL GKS_STXR(iw,index,TEXT_FONT3,TEXT_PRECISION3,TEXT_WIDTH2,
     '    TEXT_SPACING,1,ERROR,*9999)
        index=INDEX_TEXT(0,'WIDTH2','FONT4','BLACK')
        CALL GKS_STXR(iw,index,TEXT_FONT4,TEXT_PRECISION3,TEXT_WIDTH2,
     '    TEXT_SPACING,1,ERROR,*9999)

        IF(COLOUR_WS) THEN
          DO i=1,8
            index=INDEX_TEXT(0,'WIDTH1','FONT1',TEXT_COLOUR_NAME(i))
            CALL GKS_STXR(iw,index,TEXT_FONT1,TEXT_PRECISION1,
     '        TEXT_WIDTH1,TEXT_SPACING,LUT_COLOUR_INDEX(i),ERROR,*9999)
            index=INDEX_TEXT(0,'WIDTH1','FONT2',TEXT_COLOUR_NAME(i))
            CALL GKS_STXR(iw,index,TEXT_FONT2,TEXT_PRECISION3,
     '        TEXT_WIDTH1,TEXT_SPACING,LUT_COLOUR_INDEX(i),ERROR,*9999)

            index=INDEX_TEXT(0,'WIDTH2','FONT1',TEXT_COLOUR_NAME(i))
            CALL GKS_STXR(iw,index,TEXT_FONT1,TEXT_PRECISION1,
     '        TEXT_WIDTH2,TEXT_SPACING,LUT_COLOUR_INDEX(i),ERROR,*9999)
            index=INDEX_TEXT(0,'WIDTH2','FONT2',TEXT_COLOUR_NAME(i))
            CALL GKS_STXR(iw,index,TEXT_FONT2,TEXT_PRECISION3,
     '        TEXT_WIDTH2,TEXT_SPACING,LUT_COLOUR_INDEX(i),ERROR,*9999)
          ENDDO

          DO icolour=1,MAXCOLOURS-40
            ilut=INT(3.+DBLE(icolour-1)*(MAXCOLOURS-3)/(MAXCOLOURS-41))
            index=INDEX_TEXT(icolour,'WIDTH1','FONT1',' ')
            CALL GKS_STXR(iw,index,TEXT_FONT1,TEXT_PRECISION1,
     '        TEXT_WIDTH1,TEXT_SPACING,ilut,ERROR,*9999)
          ENDDO
        ENDIF

      ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS
        index=INDEX_TEXT(0,'WIDTH1','FONT1','BLACK')
        CALL PHIGS$SET_TEXT_REP(iw,index,TEXT_FONT1,TEXT_PRECISION1,
     '    TEXT_WIDTH1,TEXT_SPACING,1)
        index=INDEX_TEXT(0,'WIDTH1','FONT2','BLACK')
        CALL PHIGS$SET_TEXT_REP(iw,index,TEXT_FONT2,TEXT_PRECISION3,
     '    TEXT_WIDTH1,TEXT_SPACING,1)
        index=INDEX_TEXT(0,'WIDTH1','FONT3','BLACK')
        CALL PHIGS$SET_TEXT_REP(iw,index,TEXT_FONT3,TEXT_PRECISION3,
     '    TEXT_WIDTH1,TEXT_SPACING,1)
        index=INDEX_TEXT(0,'WIDTH1','FONT4','BLACK')
        CALL PHIGS$SET_TEXT_REP(iw,index,TEXT_FONT4,TEXT_PRECISION3,
     '    TEXT_WIDTH1,TEXT_SPACING,1)

        index=INDEX_TEXT(0,'WIDTH2','FONT1','BLACK')
        CALL PHIGS$SET_TEXT_REP(iw,index,TEXT_FONT1,TEXT_PRECISION1,
     '    TEXT_WIDTH2,TEXT_SPACING,1)
        index=INDEX_TEXT(0,'WIDTH2','FONT2','BLACK')
        CALL PHIGS$SET_TEXT_REP(iw,index,TEXT_FONT2,TEXT_PRECISION3,
     '    TEXT_WIDTH2,TEXT_SPACING,1)
        index=INDEX_TEXT(0,'WIDTH2','FONT3','BLACK')
        CALL PHIGS$SET_TEXT_REP(iw,index,TEXT_FONT3,TEXT_PRECISION3,
     '    TEXT_WIDTH2,TEXT_SPACING,1)
        index=INDEX_TEXT(0,'WIDTH2','FONT4','BLACK')
        CALL PHIGS$SET_TEXT_REP(iw,index,TEXT_FONT4,TEXT_PRECISION3,
     '    TEXT_WIDTH2,TEXT_SPACING,1)

        IF(COLOUR_WS) THEN
          DO i=1,4
            index=INDEX_TEXT(0,'WIDTH1','FONT1',TEXT_COLOUR_NAME(i))
            CALL PHIGS$SET_TEXT_REP(iw,index,TEXT_FONT1,TEXT_PRECISION1,
     '        TEXT_WIDTH1,TEXT_SPACING,LUT_COLOUR_INDEX(i))
            index=INDEX_TEXT(0,'WIDTH1','FONT2',TEXT_COLOUR_NAME(i))
            CALL PHIGS$SET_TEXT_REP(iw,index,TEXT_FONT2,TEXT_PRECISION3,
     '        TEXT_WIDTH1,TEXT_SPACING,LUT_COLOUR_INDEX(i))

            index=INDEX_TEXT(0,'WIDTH2','FONT1',TEXT_COLOUR_NAME(i))
            CALL PHIGS$SET_TEXT_REP(iw,index,TEXT_FONT1,TEXT_PRECISION1,
     '        TEXT_WIDTH2,TEXT_SPACING,LUT_COLOUR_INDEX(i))
            index=INDEX_TEXT(0,'WIDTH2','FONT2',TEXT_COLOUR_NAME(i))
            CALL PHIGS$SET_TEXT_REP(iw,index,TEXT_FONT2,TEXT_PRECISION3,
     '        TEXT_WIDTH2,TEXT_SPACING,LUT_COLOUR_INDEX(i))
          ENDDO

          DO icolour=1,MAXCOLOURS-40
            ilut=INT(3.+DBLE(icolour-1)*(MAXCOLOURS-3)/(MAXCOLOURS-41))
            index=INDEX_TEXT(icolour,'WIDTH1','FONT1',' ')
            CALL PHIGS$SET_TEXT_REP(iw,index,TEXT_FONT1,TEXT_PRECISION1,
     '        TEXT_WIDTH1,TEXT_SPACING,ilut)
          ENDDO
        ENDIF
      ENDIF

      CALL EXITS('SET_TEXT_REP')
      RETURN
 9999 CALL ERRORS('SET_TEXT_REP',ERROR)
      CALL EXITS('SET_TEXT_REP')
      RETURN 1
      END


      SUBROUTINE SETUP(iw,ERROR,*)

C#### Subroutine: SETUP
C###  Description:
C**** Performs setup operations for workstation iw.
C**** Note that default workstation type is used (GWSDEF=41 on
C**** microvax). If using microvax from another workstation use
C**** DCL command Define GKS$WSTYPE to set appropriate value.
C****
C**** IWKDEF(0) is number of defined workstations
C**** IWKDEF(noiw),noiw=1,IWKDEF(0) is list of defined workstations
C**** IWKG(iw) is 0 for nongraphics window (eg menu)
C****    "        1 for graphics output window
C**** IWKS(iw) is 0 when workstation iw is not defined
C****    "        1  "      "        iw is defined but not active
C****    "        2  "      "        iw is defined and active
C**** IWKT(iw) is 1 for GX or GKS workstation
C****    "        2 for PHIGS workstation
C****    "        3 for Frame-grabber display screen
C****
C**** Workstation_id=1 for x,y viewport in 2D or 3D
C****      "         2  "  y,z    "      "  "
C****      "         3  "  x,z    "      "  "
C****      "         4  "  map viewport (3D only)
C****      "         5  "  DT2651 frame-buffer 1
C****      "         6  "  DT2651 frame-buffer 2
C****      "         7  "  choice viewport
C****      "         8  "  help documentation viewport
C****      "         9  "  PHIGS 3D plot viewport
C****      "        10  "  nodal time history viewport
C****      "        11  "  sections viewport
C****      "        12  "  fibre angle distribution (3D only)
C****      "        13  "  sheet angle distribution (3D only)
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
C****      "        51  "  Homogeneous strain window
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
C****      "        95  "  choice viewport (window 5 menu; point/bar plot)
C****      "        96  "  choice viewport (2D scatter plot)
C****      "        97  "  choice viewport (GKS_DRAW)
C****      "        98  "  choice viewport (GKS_DRAW)
C****      "        99  "  choice viewport (GKS_DRAW)

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:cbzm00.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'
      INCLUDE 'cmiss$reference:cmis00.cmn'
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
      INCLUDE 'cmiss$reference:phig00.cmn'
      INCLUDE 'cmiss$reference:post00.cmn'
      INCLUDE 'cmiss$reference:trans00.cmn'
      INCLUDE 'cmiss$reference:view00.cmn'
      INCLUDE 'gx$path:gx.inc'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.FOR'
      INCLUDE 'SYS$LIBRARY:PHIGS$DEFS.FOR'
!     Parameter List
      INTEGER iw
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,IBEG,ICOLOUR_FLAG,IEND,IERST,
     '  iiw,ipat,NCOLOURS,noiw,
     '  TEXT_COLOUR_INDEX,TEXT_FONT1,TEXT_FONT2,TEXT_FONT3,TEXT_FONT4,
     '  TEXT_INDEX,TEXT_PRECISION1,TEXT_PRECISION2,TEXT_PRECISION3,
     '  WSTYPE
      REAL BACKGROUND_COLOUR(3),BLACK,BLUE,GREEN,PLANES(2),SCALES(2),
     '  TEXT_EXPFAC,TEXT_SPACING,WHITE,XDC(6),XDIFF,XNPC(6)
      CHARACTER NODE*10
      STRUCTURE /lightdatatype/   !AAY 8/4/90
        INTEGER*4 MODEL
        REAL*4 COLOUR(3)
        REAL*4 VECTOR(3)
      END STRUCTURE
      RECORD/lightdatatype/LIGHTDATA

      INTEGER IOSTAT,IREC
      CHARACTER ACCESS*10,FILE*100
      LOGICAL EXIST,OPENED

      CALL ENTERS('SETUP',*9999)

      IF(.NOT.GKS.AND..NOT.PHIGS) THEN
        DO i=1,4
          IZOOM(i)=1
          XNDC(1,1,i)=0.0
          XNDC(1,2,i)=1.0
          XNDC(1,3,i)=0.0
          XNDC(1,4,i)=1.0
        ENDDO
      ENDIF

      IF(DOP) THEN
        WRITE(OP_STRING,
     '    '('' iw='',I3,'' IWKS(iw)='',I2,'' IWKT(iw)='',I2,'
     '    //''' IWKG(iw)='',I2)') iw,IWKS(iw),IWKT(iw),IWKG(iw)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(iw.EQ.9.OR.iw.EQ.16.OR.iw.EQ.67.OR.iw.EQ.68) THEN !DPB. 18feb95 got
                                                           ! rid of iw.eq.3
        IWKT(iw)=2 !indicates iw is PHIGS workstation
      ELSE IF(iw.EQ.5.OR.iw.EQ.6) THEN
        IWKT(iw)=3 !indicates iw is Frame-grabber screen
      ELSE
        IWKT(iw)=1 !indicates iw is GKS   workstation
      ENDIF

      IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
        WSTYPE=32506067   !(01F000D3) try for 240 colours
      ELSE
        WSTYPE=GWSDEF
      ENDIF

      IF(IWKT(iw).EQ.1) THEN      !GKS
        CALL CHECK_GKS_OPEN(iw,ERROR,*9999)
        GKS_WS_OPEN=.TRUE. !at least one GKS workstation is open

      ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS
        IF(.NOT.PHIGS) THEN
          IF(DOP) THEN
            WRITE(OP_STRING,*)' Setup opening PHIGS'
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CALL PHIGS$OPEN_PHIGS('sys$output',0)
          PHIGS=.TRUE.
          TEXT_COLOUR_INDEX=1
          TEXT_EXPFAC=0.5
          TEXT_FONT1=1
          TEXT_FONT2=-105  !For PostScript w/s (iw = 15/16) only
          TEXT_FONT3=1
          TEXT_FONT4=1
          TEXT_INDEX=1
          TEXT_PRECISION1=GSTRP
          TEXT_PRECISION2=GCHARP
          TEXT_PRECISION3=GSTRKP
          TEXT_SPACING=0.0
          IF(.NOT.GKS) THEN
            XDISP=0.3315748
            YDISP=0.2735246
            DISP=MIN(XDISP,YDISP)
            IF(DOP) THEN
              WRITE(OP_STRING,*) ' xdisp, ydisp =',XDISP,YDISP
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            RATIO=YDISP/XDISP
!           Inquire whether colour or monochrome workstation
            IF(iw.EQ.16) THEN
              CALL PHIGS$INQ_COLOUR_FAC(PSTYPE,IERST,NCOLOURS,
     '          ICOLOUR_FLAG,NINDICES)
            ELSE
              CALL PHIGS$INQ_COLOUR_FAC(PHIGS$K_WSTYPE_DEFAULT,IERST,
     '          NCOLOURS,ICOLOUR_FLAG,NINDICES)
            ENDIF
            IF(DOP) THEN
              WRITE(OP_STRING,
     '          '('' no colours='',I8,'' ICOLOUR_FLAG='',I2,'
     '          //''' no predefined indices='',I4)') NCOLOURS,
     '          ICOLOUR_FLAG,NINDICES
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            IF(NINDICES.GT.2) THEN
              COLOUR_WS=.TRUE.
              IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          MAXCOLOURS=249
              ELSE
          MAXCOLOURS=240
              ENDIF
            ELSE
              COLOUR_WS=.FALSE.
              MAXCOLOURS=2
            ENDIF
            CALL SET_COLOUR_LUT(COLOUR_LUT,ERROR,*9999)
          ENDIF
        ENDIF
      ENDIF

! iw=1
      IF(iw.EQ.1.AND.IWKS(1).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
C PJH 18Feb95. Don't use until GX in CMgui
C       IF(USE_SOCKET) THEN
C         CALL X_OPWK(iw,0.,0.48,0.51,0.99,ERROR,*9999)
C       ELSE IF(.NOT.USE_SOCKET) THEN
        IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
          CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
          NODE=DISPLAY_NAME(IBEG:IEND)//'::0AA'
          CALL STRING_TRIM(NODE,IBEG,IEND)
          CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
        ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          open(unit=101,file='fem.w1',status='unknown',
     '      dispose='delete')
          CALL GKS_OPWK(1,101,GWSDEF,ERROR,*9999)
          IF(NJT.LE.2) THEN
            CALL GKS_SWKVP(iw,0.,SCREEN_WIDTH/100.*XDISP,
     '        YDISP-SCREEN_WIDTH/100.*XDISP,YDISP,ERROR,*9999)
          ELSE
            CALL GKS_SWKVP(iw,0.,0.48*DISP,0.51*DISP,0.99*DISP,
     '        ERROR,*9999)
          ENDIF
        ENDIF
C       ENDIF
        IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
          CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
          CALL SET_POLYLINE_REP(iw,ERROR,*9999)
          CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
          CALL SET_TEXT_REP(iw,ERROR,*9999)
          CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)
        ENDIF
        IF(NJT.LE.2) THEN
          CALL GKS_SWN(iw,XMIN,XMAX,YMIN,YMAX,ERROR,*9999)
        ELSE IF(NJT.EQ.3) THEN
          CALL GKS_SWN(iw,XMIN,XMAX,ZMIN,ZMAX,ERROR,*9999)
        ENDIF

! iw=2
      ELSE IF(iw.EQ.2.AND.IWKS(2).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
c       IF(USE_SOCKET) THEN
c         CALL X_OPWK(iw,0.50,0.98,0.51,0.99,ERROR,*9999)
c       ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0AB'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=102,file='fem.w2',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(2,102,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.50*DISP,0.98*DISP,0.51*DISP,0.99*DISP,
     '        ERROR,*9999)
          ENDIF
c       ENDIF
        IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
          CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
          CALL SET_POLYLINE_REP(iw,ERROR,*9999)
          CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
          CALL SET_TEXT_REP(iw,ERROR,*9999)
          CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)
        ENDIF
        CALL GKS_SWN(iw,YMIN,YMAX,ZMIN,ZMAX,ERROR,*9999)

! iw=3
      ELSE IF(iw.EQ.3.AND.IWKS(3).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
c       IF(USE_SOCKET) THEN
c         CALL X_OPWK(iw,0.50,0.98,0.51,0.99,ERROR,*9999)
c       ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0AB'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=103,file='fem.w2',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(3,103,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.50*DISP,0.98*DISP,0.51*DISP,0.99*DISP,
     '        ERROR,*9999)
          ENDIF
c       ENDIF
        IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
          CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
          CALL SET_POLYLINE_REP(iw,ERROR,*9999)
          CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
          CALL SET_TEXT_REP(iw,ERROR,*9999)
          CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)
        ENDIF
        CALL GKS_SWN(iw,YMIN,YMAX,ZMIN,ZMAX,ERROR,*9999)

! iw=3 PHIGS version - not currently used (15Feb95)
      ELSE IF(iw.EQ.3.AND.IWKS(3).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=2 !indicates PHIGS workstation
        IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
          CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
          NODE=DISPLAY_NAME(IBEG:IEND)//'::0AC'
          CALL STRING_TRIM(NODE,IBEG,IEND)
          CALL PHIGS$OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE)
        ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          CALL PHIGS$OPEN_WS(iw,'FEM.W3',41)
          CALL PHIGS$SET_COLOUR_MODEL(iw,PHIGS$K_RGB)
          CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
          CALL SET_POLYLINE_REP(iw,ERROR,*9999)
          CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
          CALL SET_TEXT_REP(iw,ERROR,*9999)
          XNPC(1)= 0.0
          XNPC(2)= 1.0
          XNPC(3)= 0.0
          XNPC(4)= 1.0
          XNPC(5)= 0.0
          XNPC(6)= 1.0
          CALL PHIGS$SET_WS_WINDOW3(iw,XNPC)
          XDC(1)= 0.50*DISP
          XDC(2)= 0.98*DISP
          XDC(3)= 0.0
          XDC(4)= 0.48*DISP
          XDC(5)= 0.0
          XDC(6)= 0.0
          CALL PHIGS$SET_WS_VIEWPORT(iw,XDC)

C ***     set up two lights for surface shading  !AAY 8/4/90
C         first is an ambient light
          LIGHTDATA.MODEL=PHIGS$K_RGB
          LIGHTDATA.COLOUR(1)=0.7
          LIGHTDATA.COLOUR(2)=0.7
          LIGHTDATA.COLOUR(3)=0.7
          CALL PHIGS$SET_LIGHT_SRC_REP(iw,1,PHIGS$K_LIGHT_TYPE_AMBIENT,
     '      4*4,LIGHTDATA)
C         second is a vector light
          LIGHTDATA.MODEL=PHIGS$K_RGB
          LIGHTDATA.COLOUR(1)=0.7
          LIGHTDATA.COLOUR(2)=0.7
          LIGHTDATA.COLOUR(3)=0.7
          LIGHTDATA.VECTOR(1)=-1.0
          LIGHTDATA.VECTOR(2)= 1.0
          LIGHTDATA.VECTOR(3)=-1.0
          CALL PHIGS$SET_LIGHT_SRC_REP(iw,2,PHIGS$K_LIGHT_TYPE_INFINITE,
     '      4*7,LIGHTDATA)

C ***     set depth cue representation   !AAY 8/4/90
          PLANES(1)=0.0
          PLANES(2)=1.0
          SCALES(1)=0.0
          SCALES(2)=1.0
C         this works so long as the colour is white or black (or grey)?
!PJH 25-feb-1992 next 2 statements commented since giving error
          CALL PHIGS$SET_DEPTH_CUE_REP(iw,1,
     '      PHIGS$K_DEPTH_CUE_SUPPRESSED,PLANES,SCALES,PHIGS$K_RGB,
     '      WHITE)
          CALL PHIGS$SET_DEPTH_CUE_REP(iw,2,PHIGS$K_DEPTH_CUE_ALLOWED,
     '      PLANES,SCALES,PHIGS$K_RGB,WHITE)

C ***     Define transformation and view matrices
          CALL PHIG(iw,ERROR,*9999)
          CALL PHIGS$SET_VIEW_REP3(iw,1,A_ORIENT,A_MAP,NPC_CLIP,
     '      PHIGS$K_NOCLIP,PHIGS$K_NOCLIP,PHIGS$K_NOCLIP)
          NTSG=NTSG+1
          ISVIEW=NTSG
C ***     isview is the root structure for all the phigs output
C ***     Insert global transformation matrix into structure ISVIEW
          CALL PHIGS$SET_EDIT_MODE(PHIGS$K_EDIT_INSERT)
          CALL PHIGS$OPEN_STRUCT(ISVIEW)

          CALL PHIGS$SET_GLOBAL_XFORM3(A_TRANS)
          CALL PHIGS$SET_VIEW_INDEX(1)
          CALL PHIGS$SET_CHAR_HEIGHT(0.05*DIAG/SQRT(12.))
C ***     use bundled representation for all attributes AAY 5/4/90
          DO i=0,20
            CALL PHIGS$SET_INDIV_ASF(i,PHIGS$K_ASF_BUNDLED)
          ENDDO

C ***     Set edge representation
          CALL PHIGS$SET_EDGE_REP(iw,1,PHIGS$K_EDGE_FLAG_ON,
     '      PHIGS$K_EDGE_SOLID,1,1)
          CALL PHIGS$SET_EDGE_REP(iw,2,PHIGS$K_EDGE_FLAG_OFF,
     '      PHIGS$K_EDGE_SOLID,1,1)

C ***     Set 1st extended interior representation : black hatch no lightsource
          CALL PHIGS$SET_EXT_INT_REP(iw,1,   !workstation id and index
     '      PHIGS$K_INTSTYLE_HATCH,PHIGS$K_INTSTYLE_HATCH, !front and back styles
     '      -1,-1,                           !front and back style indices
     '      PHIGS$K_RGB,BLACK,               !rgb front colour
     '      PHIGS$K_RGB,BLACK,               !rgb back colour
C           these next two don't seem to do anything?
     '      PHIGS$K_SHADE_NONE,PHIGS$K_SHADE_NONE, !front & back shading method
     '      PHIGS$K_LIGHTING_NONE,PHIGS$K_LIGHTING_NONE, !front & back lighting
     '      1.0,1.0,1.0, !ambient,diffuse and specular reflection coeffs
     '      PHIGS$K_RGB,WHITE,100.0,0.0,!rgb specular colour,exponent,transparency
     '      1.0,1.0,1.0, !ambient,diffuse and specular reflection back faces
     '      PHIGS$K_RGB,WHITE,100.0,0.0,!rgb specular values for back faces
     '      1,4, !surface approx type and value
     '      1,4) !trim approx type and value

C ***     2nd interior representation : white solid no lightsource
          CALL PHIGS$SET_EXT_INT_REP(iw,2,   !workstation id and index
     '      PHIGS$K_INTSTYLE_SOLID,PHIGS$K_INTSTYLE_SOLID, !front and back styles
     '      1,1,                             !front and back style indices
     '      PHIGS$K_INDEXED_COLOUR,0,               !rgb front colour
     '      PHIGS$K_INDEXED_COLOUR,0,               !rgb back colour
C           these next two don't seem to do anything?
     '      PHIGS$K_SHADE_NONE,PHIGS$K_SHADE_NONE, !front & back shading method
     '      PHIGS$K_LIGHTING_NONE,
     '      PHIGS$K_LIGHTING_NONE, !front & back lighting
     '      1.0,1.0,1.0, !ambient,diffuse and specular reflection coeffs
     '      PHIGS$K_INDEXED_COLOUR,0,100.0,1.0,!rgb specular colour,exponent,transparency
     '      1.0,1.0,1.0, !ambient,diffuse and specular reflection back faces
     '      PHIGS$K_INDEXED_COLOUR,0,100.0,1.0,!rgb specular values for back faces
     '      1,4, !surface approx type and value
     '      1,4) !trim approx type and value

C ***     3rd interior representation : white hollow no lightsource
          CALL PHIGS$SET_EXT_INT_REP(iw,3,   !workstation id and index
     '      PHIGS$K_INTSTYLE_HOLLOW,
     '      PHIGS$K_INTSTYLE_HOLLOW,         !front and back styles
     '      1,1,                             !front and back style indices
     '      PHIGS$K_RGB,BLACK,               !rgb front colour
     '      PHIGS$K_RGB,BLACK,               !rgb back colour
C           these next two don't seem to do anything?
     '      PHIGS$K_SHADE_NONE,PHIGS$K_SHADE_NONE, !front & back shading method
     '      PHIGS$K_LIGHTING_NONE,
     '      PHIGS$K_LIGHTING_NONE, !front & back lighting
     '      0.5,0.5,0.9, !ambient,diffuse and specular reflection coeffs
     '      PHIGS$K_RGB,WHITE,100.0,1.0,!rgb specular colour,exponent,transparency
     '      0.5,0.5,0.9, !ambient,diffuse and specular reflection back faces
     '      PHIGS$K_RGB,WHITE,0.9,0.4,!rgb specular values for back faces
     '      1,4, !surface approx type and value
     '      1,4) !trim approx type and value

C ***     4th interior representation : white solid with lightsource
          CALL PHIGS$SET_EXT_INT_REP(iw,4,   !workstation id and index
     '      PHIGS$K_INTSTYLE_SOLID,
     '      PHIGS$K_INTSTYLE_SOLID,          !front and back styles
     '      1,1,                             !front and back style indices
     '      PHIGS$K_INDEXED_COLOUR,2,        !index front colour
     '      PHIGS$K_INDEXED_COLOUR,2,        !index back colour
C           these next two don't seem to do anything?
     '      PHIGS$K_SHADE_NONE,PHIGS$K_SHADE_NONE, !front & back shading method
     '      PHIGS$K_LIGHTING_DIFFUSE,
     '      PHIGS$K_LIGHTING_DIFFUSE, !front & back lighting
     '      0.3,0.5,0.9, !ambient,diffuse and specular reflection coeffs
     '      PHIGS$K_INDEXED_COLOUR,2,100.0,1.0,!rgb specular colour,exponent,transparency
     '      1.0,1.0,1.0, !ambient,diffuse and specular reflection back faces
     '      PHIGS$K_INDEXED_COLOUR,2,100.0,1.0,!rgb specular values for back faces
     '      1,4, !surface approx type and value
     '      1,4) !trim approx type and value

C ***     5th interior representation : green/blue solid shaded
          CALL PHIGS$SET_EXT_INT_REP(iw,5,   !workstation id and index
     '      PHIGS$K_INTSTYLE_SOLID,
     '      PHIGS$K_INTSTYLE_SOLID,          !front and back styles
     '      1,1,                             !front and back style indices
     '      PHIGS$K_RGB,GREEN,               !rgb front colour
     '      PHIGS$K_RGB,BLUE,               !rgb back colour
C           these next two don't seem to do anything?
     '      PHIGS$K_SHADE_NONE,PHIGS$K_SHADE_NONE, !front & back shading method
     '      PHIGS$K_LIGHTING_DIFFUSE,
     '      PHIGS$K_LIGHTING_SPECULAR, !front & back lighting
     '      0.5,0.5,0.9, !ambient,diffuse and specular reflection coeffs
     '      PHIGS$K_RGB,WHITE,100.0,1.0,!rgb specular colour,exponent,transparency
     '      1.0,0.5,0.9, !ambient,diffuse and specular reflection back faces
     '      PHIGS$K_RGB,BLUE,0.9,0.4,!rgb specular values for back faces
     '      1,4, !surface approx type and value
     '      1,4) !trim approx type and value

          CALL PHIGS$SET_GEOM_NORMAL_CAL_MODE(PHIGS$K_FORWARD)

C ***     turn on painters hidden line/surface removal              !AAY 5/4/90
          CALL PHIGS$SET_HLHSR_ID(1)

          CALL PHIGS$CLOSE_STRUCT()
          CALL PHIGS$POST_STRUCT(iw,ISVIEW,1.0)
        ENDIF

! iw=4
      ELSE IF(iw.EQ.4.AND.IWKS(4).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
c       IF(USE_SOCKET) THEN
c         CALL X_OPWK(iw,0.2,1.0,0.2,1.0,ERROR,*9999)
c       ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0AD'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=104,file='fem.w4',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(4,104,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.2*DISP,1.0*DISP,0.2*DISP,1.0*DISP,
     '      ERROR,*9999)
c       ENDIF
        IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
          CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
          CALL SET_POLYLINE_REP(iw,ERROR,*9999)
          CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
          CALL SET_TEXT_REP(iw,ERROR,*9999)
          CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)
        ENDIF
        CALL GKS_SWN(iw,-1.0,1.0,-1.0,1.0,ERROR,*9999)

! iw=5
      ELSE IF(iw.EQ.5.AND.IWKS(5).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=3 !indicates Frame-grabber screen
C       At some time this window on the frame buffer should be opened thru GKS
C       open(unit=105,file='fem.w5',status='unknown',dispose='delete')
C       CALL GKS_OPWK(5,105,GWSDEF,ERROR,*9999)
C       CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)

! iw=6
      ELSE IF(iw.EQ.6.AND.IWKS(6).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=3 !indicates Frame-grabber screen
C       At some time this window on the frame buffer should be opened thru GKS
C       open(unit=106,file='fem.w6',status='unknown',dispose='delete')
C       CALL GKS_OPWK(6,106,GWSDEF,ERROR,*9999)
C       CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)

! iw=7
      ELSE IF(iw.EQ.7.AND.IWKS(7).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
c       IF(USE_SOCKET) THEN
c         CALL X_OPWK(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
c    '        ERROR,*9999)
c       ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0AG'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
C            CALL GKS_SWKVP(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
C     '        ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=107,file='fem.w7',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(7,107,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
          ENDIF
c       ENDIF
        IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
          CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
          CALL SET_POLYLINE_REP(iw,ERROR,*9999)
          CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
          CALL SET_TEXT_REP(iw,ERROR,*9999)
          CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)
        ENDIF
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)

! iw=8
      ELSE IF(iw.EQ.8.AND.IWKS(8).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
c       IF(USE_SOCKET) THEN
c         CALL X_OPWK(iw,0.25,0.75,0.25,0.75,
c    '        ERROR,*9999)
c       ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0AH'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
C            CALL GKS_SWKVP(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
C     '        ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=108,file='fem.w8',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(8,108,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.2*DISP,1.0*DISP,0.295*DISP,0.8*DISP,
     '        ERROR,*9999)
          ENDIF
c       ENDIF
        IF(WINDOW_TYPE.EQ.'VWS') THEN
          CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
          CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
          CALL SET_POLYLINE_REP(iw,ERROR,*9999)
          CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
          CALL SET_TEXT_REP(iw,ERROR,*9999)
          CALL GKS_SWKWN(iw,0.,1.,0.3688,1.,ERROR,*9999)
        ENDIF
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)

! iw=9
      ELSE IF(iw.EQ.9.AND.IWKS(9).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=2 !indicates PHIGS workstation
        IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
          CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
          NODE=DISPLAY_NAME(IBEG:IEND)//'::0AI'
          CALL STRING_TRIM(NODE,IBEG,IEND)
          CALL PHIGS$OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE)
        ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          CALL PHIGS$OPEN_WS(9,'FEM.W29',41)
        ENDIF
        XNPC(1)= 0.0
        XNPC(2)= 1.0
        XNPC(3)= 0.0
        XNPC(4)= 1.0
C       XNPC(5)=-1.0   !old phigs version
C       XNPC(6)= 0.0   !old phigs version
        XNPC(5)= 0.0   !7-feb-1990
        XNPC(6)= 1.0   !7-feb-1990
        CALL PHIGS$SET_WS_WINDOW3(iw,XNPC)
        XDC(1)= 0.0
        XDC(2)= 0.48*XDISP
        XDC(3)= 0.51*DISP
        XDC(4)= 0.99*DISP
        CALL PHIGS$SET_WS_VIEWPORT(iw,XDC)
C ***   Define transformation and view matrices
        CALL PHIG(iw,ERROR,*9999)
        CALL PHIGS$SET_VIEW_REP3(iw,1,A_ORIENT,A_MAP,NPC_CLIP,
     '    PHIGS$K_NOCLIP,PHIGS$K_NOCLIP,PHIGS$K_NOCLIP)

! iw=10
      ELSE IF(iw.EQ.10.AND.IWKS(10).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.5,0.99,0.77,0.99,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0AJ'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=110,file='fem.w10',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(10,110,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.5*XDISP,0.99*XDISP,YDISP-0.23*XDISP,
     '        YDISP,ERROR,*9999)
          ENDIF
        ENDIF
        IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
          CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
          CALL SET_POLYLINE_REP(iw,ERROR,*9999)
          CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
          CALL SET_TEXT_REP(iw,ERROR,*9999)
          CALL GKS_SWKWN(iw,0.,1.,0.5-0.23/0.49/2.,0.5+0.23/0.49/2.,
     '      ERROR,*9999)
        ENDIF
        CALL GKS_SWN(iw,-0.14,1.1,-2.3,2.3,ERROR,*9999)

! iw=11
      ELSE IF(iw.EQ.11.AND.IWKS(11).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.5,0.99,0.51,0.74,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0AK'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=11,file='fem.w11',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(11,11,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.5*XDISP,0.99*XDISP,YDISP-0.49*XDISP,
     '        YDISP-0.26*XDISP,ERROR,*9999)
          ENDIF
        ENDIF
        IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
          CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
          CALL SET_POLYLINE_REP(iw,ERROR,*9999)
          CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
          CALL SET_TEXT_REP(iw,ERROR,*9999)
          CALL GKS_SWKWN(iw,0.,1.,0.5-0.23/0.49/2.,0.5+0.23/0.49/2.,
     '      ERROR,*9999)
        ENDIF
        CALL GKS_SWN(iw,-0.14,1.1,-2.3,2.3,ERROR,*9999)

! iw=12
      ELSE IF(iw.EQ.12.AND.IWKS(12).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.76,0.99,0.51,0.74,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0AL'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=12,file='fem.w12',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(12,12,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.76*XDISP,0.99*XDISP,YDISP-0.49*XDISP,
     '        YDISP-0.26*XDISP,ERROR,*9999)
          ENDIF
        ENDIF
        IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
          CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
          CALL SET_POLYLINE_REP(iw,ERROR,*9999)
          CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
          CALL SET_TEXT_REP(iw,ERROR,*9999)
          CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)
        ENDIF
        CALL GKS_SWN(iw,-0.20,1.1,-1.1*PI/2.,1.1*PI/2.,ERROR,*9999)

! iw=13
      ELSE IF(iw.EQ.13.AND.IWKS(13).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.30,0.70,0.20,0.60,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0AM'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=13,file='fem.w13',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(13,13,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.30*XDISP,0.70*XDISP,YDISP-0.80*XDISP,
     '        YDISP-0.40*XDISP,ERROR,*9999)
          ENDIF
        ENDIF
        IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
          CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
          CALL SET_POLYLINE_REP(iw,ERROR,*9999)
          CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
          CALL SET_TEXT_REP(iw,ERROR,*9999)
          CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)
          CALL GKS_SWKVP(iw,0.30*XDISP,0.70*XDISP,YDISP-0.80*XDISP,
     '      YDISP-0.40*XDISP,ERROR,*9999)
        ENDIF
        XDIFF=XMAX-XMIN
        CALL GKS_SWN(iw,XMIN,XMAX,-0.2*XDIFF,0.8*XDIFF,ERROR,*9999)

! iw=15
      ELSE IF(iw.EQ.15.AND.IWKS(15).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        CALL STRING_TRIM(FILE00,IBEG,IEND)
        POSTFILE=FILE00(IBEG:IEND)
        IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
          POSTFILE=FILE00(IBEG:IEND)//'.ps'
          IF(PORTRAIT) THEN
C MPN 25-Mar-94: not implemented yet
C            CALL OPPRNT(POSTFILE,MPRTRT,0.0,COLOURPS)
          ELSE
C MPN 25-Mar-94: not implemented yet
C            CALL OPPRNT(POSTFILE,MLNSCP,0.0,COLOURPS)
          ENDIF
        ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          IF(DOP) THEN
            WRITE(OP_STRING,*)
     '        ' opening ',POSTFILE(IBEG:IEND),'.gks'
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CALL OPENF(15,'DISK',POSTFILE(IBEG:IEND)//'.gks','NEW',
     '      'SEQUEN','FORMATTED',256,ERROR,*9999)
C CPB   5/8/93 All the workstation types are now defined in dewind
C          IF(PORTRT) THEN !sets portrait bitmask
CC*            CALL GKS_OPWK(15,15,GPTSC.OR.GPRTRT.OR.GSIZA,ERROR,*9999)
CC* Co  lour, A4, portrait
C            CALL GKS_OPWK(15,15,'1050003E'X,ERROR,*9999)
C            IF(DOP) THEN
C              WRITE(OP_STRING,*) 'portrait'
C              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C            ENDIF
C          ELSE            !sets landscape bitmask
CC*            CALL GKS_OPWK(15,15,GPTSC.OR.GLDSCP.OR.GSIZA,ERROR,*9999)
CC* Co  lour, A4, landscape
C            CALL GKS_OPWK(15,15,'0050003E'X,ERROR,*9999)
C            IF(DOP) THEN
C              WRITE(OP_STRING,*) 'landscape'
C              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C            ENDIF
C          ENDIF
          CALL GKS_OPWK(15,15,PSTYPE,ERROR,*9999)
!         define COLOUR_LUT array
          CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
!         define background colour (index 0)
          BACKGROUND_COLOUR(1)=1.0
          BACKGROUND_COLOUR(2)=1.0
          BACKGROUND_COLOUR(3)=1.0
          CALL SET_COLOUR_ONE(iw,0,BACKGROUND_COLOUR,ERROR,*9999)
          CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
          CALL SET_POLYLINE_REP(iw,ERROR,*9999)
          CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
          CALL SET_TEXT_REP(iw,ERROR,*9999)
c any   font adjustment for the postscript ws sb done in set_text_rep? PJH
c         CALL GKS_STXR(15,1,TEXT_FONT2,TEXT_PRECISION2,
c    '      0.4,TEXT_SPACING,TEXT_COLOUR_INDEX,ERROR,*9999)
c         CALL GKS_STXR(15,2,TEXT_FONT2,TEXT_PRECISION2,
c    '      0.4,TEXT_SPACING,TEXT_COLOUR_INDEX,ERROR,*9999)
C ***     Set up fill area patterns to be 16 shades of grey
          IF(COLOUR_WS) THEN
            CALL GKS_SFAR(iw,1,GSOLID,14,1,ERROR,*9999)
            CALL GKS_SFAR(iw,16,GSOLID,14,0,ERROR,*9999)
            DO ipat=2,15
              CALL GKS_SFAR(iw,ipat,GPATTR,12+(19-ipat),1,ERROR,*9999)
            ENDDO
          ENDIF
          CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)
        ENDIF
        CALL GKS_SWN(iw,XMIN,XMAX,YMIN,YMAX,ERROR,*9999)

! iw=16
      ELSE IF(iw.EQ.16.AND.IWKS(16).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=2 !indicates PHIGS workstation
        CALL STRING_TRIM(FILE00,IBEG,IEND)
        POSTFILE=FILE00(IBEG:IEND)
        IF(DOP) THEN
          WRITE(OP_STRING,*)' opening ',POSTFILE(IBEG:IEND),'.phigs'
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
C CPB 5/8/93 PUT iwstype INTO DEWIND CALLED PSTYPE
C ***   1050003e hex (colour postscript portrait with A4 paper)
C ***   is 293678398 decimal
C        iwstype=273678398
        CALL PHIGS$OPEN_WS(16,POSTFILE(IBEG:IEND)//'.phigs',PSTYPE)
        CALL PHIGS$SET_COLOUR_MODEL(iw,PHIGS$K_RGB)
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
!       define background colour (index 0)
        BACKGROUND_COLOUR(1)=1.0
        BACKGROUND_COLOUR(2)=1.0
        BACKGROUND_COLOUR(3)=1.0
        CALL SET_COLOUR_ONE(iw,0,BACKGROUND_COLOUR,ERROR,*9999)
        CALL SET_POLYLINE_REP(iw,ERROR,*9999)
        CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
c       CALL PHIGS$SET_PMARKER_REP(iw,1,PHIGS$K_MARKERTYPE_DOT,1.,1)
c       CALL PHIGS$SET_PMARKER_REP(iw,2,PHIGS$K_MARKERTYPE_PLUS,1.,1)
c       CALL PHIGS$SET_PMARKER_REP(iw,3,PHIGS$K_MARKERTYPE_ASTERISK,1.,1)
c       CALL PHIGS$SET_PMARKER_REP(iw,4,PHIGS$K_MARKERTYPE_CIRCLE,1.,1)
c       CALL PHIGS$SET_PMARKER_REP(iw,5,PHIGS$K_MARKERTYPE_DIAG_CROSS,1.,1)
c       CALL PHIGS$SET_PLINE_REP(iw,1,PHIGS$K_LINETYPE_SOLID,1.,1)
c       CALL PHIGS$SET_PLINE_REP(iw,2,PHIGS$K_LINETYPE_DASHED,1.,1)
c       CALL PHIGS$SET_PLINE_REP(iw,3,PHIGS$K_LINETYPE_DOTTED,1.,1)
c       CALL PHIGS$SET_PLINE_REP(iw,4,PHIGS$K_LINETYPE_DASHED_DOTTED,1.,1)
c       CALL PHIGS$SET_TEXT_REP(iw,1,TEXT_FONT2,TEXT_PRECISION1,TEXT_EXPFAC,
c    '    TEXT_SPACING,TEXT_COLOUR_INDEX)
c       CALL PHIGS$SET_TEXT_REP(iw,2,TEXT_FONT2,TEXT_PRECISION1,TEXT_EXPFAC,
c    '    TEXT_SPACING,TEXT_COLOUR_INDEX)
        XNPC(1)= 0.0
        XNPC(2)= 1.0
        XNPC(3)= 0.0
        XNPC(4)= 1.0
        XNPC(5)= 0.0
        XNPC(6)= 1.0
        CALL PHIGS$SET_WS_WINDOW3(iw,XNPC)
C ***   Define transformation and view matrices
C CPB 2/9/92 - IF THE IW=3 WORKSTATION IS SETUP USE ITS VIEW AND ORIENTATION
C ETC. MATRICES
        IF(IWKS(3).EQ.0) THEN
          CALL PHIG(iw,ERROR,*9999)
          CALL PHIGS$SET_VIEW_REP3(iw,1,A_ORIENT,A_MAP,NPC_CLIP,
     '      PHIGS$K_NOCLIP,PHIGS$K_NOCLIP,PHIGS$K_NOCLIP)
        ELSE
          CALL PHIGS$SET_VIEW_REP3(iw,1,A_ORIENT_NEW,A_MAP,NPC_CLIP,
     '      PHIGS$K_NOCLIP,PHIGS$K_NOCLIP,PHIGS$K_NOCLIP)
        ENDIF
C CPB 29/3/92 added below code from iw.eq.3 case above
C ***   set up two lights for surface shading  !AAY 8/4/90
C       first is an ambient light
        LIGHTDATA.MODEL=PHIGS$K_RGB
        LIGHTDATA.COLOUR(1)=0.7
        LIGHTDATA.COLOUR(2)=0.7
        LIGHTDATA.COLOUR(3)=0.7
        CALL PHIGS$SET_LIGHT_SRC_REP(iw,1,PHIGS$K_LIGHT_TYPE_AMBIENT,
     '    4*4,LIGHTDATA)
C       second is a vector light
        LIGHTDATA.MODEL=PHIGS$K_RGB
        LIGHTDATA.COLOUR(1)=0.7
        LIGHTDATA.COLOUR(2)=0.7
        LIGHTDATA.COLOUR(3)=0.7
        LIGHTDATA.VECTOR(1)=-1.0
        LIGHTDATA.VECTOR(2)= 1.0
        LIGHTDATA.VECTOR(3)=-1.0
        CALL PHIGS$SET_LIGHT_SRC_REP(iw,2,PHIGS$K_LIGHT_TYPE_INFINITE,
     '    4*7,LIGHTDATA)

C ***   set depth cue representation   !AAY 8/4/90
        PLANES(1)=0.0
        PLANES(2)=1.0
        SCALES(1)=0.0
        SCALES(2)=1.0
C       this works so long as the colour is white or black (or grey)?
!PJH 25-feb-1992 next 2 statements commented since giving error
        CALL PHIGS$SET_DEPTH_CUE_REP(iw,1,PHIGS$K_DEPTH_CUE_SUPPRESSED,
     '    PLANES,SCALES,PHIGS$K_RGB,WHITE)
        CALL PHIGS$SET_DEPTH_CUE_REP(iw,2,PHIGS$K_DEPTH_CUE_ALLOWED,
     '    PLANES,SCALES,PHIGS$K_RGB,WHITE)

        NTSG=NTSG+1
        ISVIEW=NTSG
C ***   isview is the root structure for all the phigs output
C ***   Insert global transformation matrix into structure ISVIEW
        CALL PHIGS$SET_EDIT_MODE(PHIGS$K_EDIT_INSERT)
        CALL PHIGS$OPEN_STRUCT(ISVIEW)

        CALL PHIGS$SET_GLOBAL_XFORM3(A_TRANS)
        CALL PHIGS$SET_VIEW_INDEX(1)
        CALL PHIGS$SET_CHAR_HEIGHT(0.05*DIAG/SQRT(12.))
C ***   use bundled representation for all attributes AAY 5/4/90
        DO i=0,20
          CALL PHIGS$SET_INDIV_ASF(I,PHIGS$K_ASF_BUNDLED)
        ENDDO
C ***   Set edge representation
        CALL PHIGS$SET_EDGE_REP(iw,1,PHIGS$K_EDGE_FLAG_ON,
     '    PHIGS$K_EDGE_SOLID,1,1)
        CALL PHIGS$SET_EDGE_REP(iw,2,PHIGS$K_EDGE_FLAG_OFF,
     '    PHIGS$K_EDGE_SOLID,1,1)

C ***   Set 1st extended interior representation : black hatch no lightsource
        CALL PHIGS$SET_EXT_INT_REP(iw,1,   !workstation id and index
     '    PHIGS$K_INTSTYLE_HATCH,PHIGS$K_INTSTYLE_HATCH, !front and back styles
     '    -1,-1,                           !front and back style indices
     '    PHIGS$K_RGB,BLACK,               !rgb front colour
     '    PHIGS$K_RGB,BLACK,               !rgb back colour
C         these next two don't seem to do anything?
     '    PHIGS$K_SHADE_NONE,PHIGS$K_SHADE_NONE, !front & back shading method
     '    PHIGS$K_LIGHTING_NONE,PHIGS$K_LIGHTING_NONE, !front & back lighting
     '    1.0,1.0,1.0, !ambient,diffuse and specular reflection coeffs
     '    PHIGS$K_RGB,WHITE,100.0,0.0,!rgb specular colour,exponent,transparency
     '    1.0,1.0,1.0, !ambient,diffuse and specular reflection back faces
     '    PHIGS$K_RGB,WHITE,100.0,0.0,!rgb specular values for back faces
     '    1,4, !surface approx type and value
     '    1,4) !trim approx type and value

C ***   2nd interior representation : white solid no lightsource
        CALL PHIGS$SET_EXT_INT_REP(iw,2,   !workstation id and index
     '    PHIGS$K_INTSTYLE_SOLID,PHIGS$K_INTSTYLE_SOLID, !front and back styles
     '    1,1,                             !front and back style indices
     '    PHIGS$K_INDEXED_COLOUR,0,               !rgb front colour
     '    PHIGS$K_INDEXED_COLOUR,0,               !rgb back colour
C         these next two don't seem to do anything?
     '    PHIGS$K_SHADE_NONE,PHIGS$K_SHADE_NONE, !front & back shading method
     '    PHIGS$K_LIGHTING_NONE,
     '    PHIGS$K_LIGHTING_NONE, !front & back lighting
     '    1.0,1.0,1.0, !ambient,diffuse and specular reflection coeffs
     '    PHIGS$K_INDEXED_COLOUR,0,100.0,1.0,!rgb specular colour,exponent,transparency
     '    1.0,1.0,1.0, !ambient,diffuse and specular reflection back faces
     '    PHIGS$K_INDEXED_COLOUR,0,100.0,1.0,!rgb specular values for back faces
     '    1,4, !surface approx type and value
     '    1,4) !trim approx type and value

C ***   3rd interior representation : white hollow no lightsource
        CALL PHIGS$SET_EXT_INT_REP(iw,3,   !workstation id and index
     '    PHIGS$K_INTSTYLE_HOLLOW,
     '    PHIGS$K_INTSTYLE_HOLLOW,         !front and back styles
     '    1,1,                             !front and back style indices
     '    PHIGS$K_RGB,BLACK,               !rgb front colour
     '    PHIGS$K_RGB,BLACK,               !rgb back colour
C         these next two don't seem to do anything?
     '    PHIGS$K_SHADE_NONE,PHIGS$K_SHADE_NONE, !front & back shading method
     '    PHIGS$K_LIGHTING_NONE,
     '    PHIGS$K_LIGHTING_NONE, !front & back lighting
     '    0.5,0.5,0.9, !ambient,diffuse and specular reflection coeffs
     '    PHIGS$K_RGB,WHITE,100.0,1.0,!rgb specular colour,exponent,transparency
     '    0.5,0.5,0.9, !ambient,diffuse and specular reflection back faces
     '    PHIGS$K_RGB,WHITE,0.9,0.4,!rgb specular values for back faces
     '    1,4, !surface approx type and value
     '    1,4) !trim approx type and value

C ***   4th interior representation : white solid with lightsource
        CALL PHIGS$SET_EXT_INT_REP(iw,4,   !workstation id and index
     '    PHIGS$K_INTSTYLE_SOLID,
     '    PHIGS$K_INTSTYLE_SOLID,          !front and back styles
     '    1,1,                             !front and back style indices
     '    PHIGS$K_INDEXED_COLOUR,2,        !index front colour
     '    PHIGS$K_INDEXED_COLOUR,2,        !index back colour
C         these next two don't seem to do anything?
     '    PHIGS$K_SHADE_NONE,PHIGS$K_SHADE_NONE, !front & back shading method
     '    PHIGS$K_LIGHTING_DIFFUSE,
     '    PHIGS$K_LIGHTING_DIFFUSE, !front & back lighting
     '    0.3,0.5,0.9, !ambient,diffuse and specular reflection coeffs
     '    PHIGS$K_INDEXED_COLOUR,2,100.0,1.0,!rgb specular colour,exponent,transparency
     '    1.0,1.0,1.0, !ambient,diffuse and specular reflection back faces
     '    PHIGS$K_INDEXED_COLOUR,2,100.0,1.0,!rgb specular values for back faces
     '    1,4, !surface approx type and value
     '    1,4) !trim approx type and value

C ***   5th interior representation : green/blue solid shaded
        CALL PHIGS$SET_EXT_INT_REP(iw,5,   !workstation id and index
     '    PHIGS$K_INTSTYLE_SOLID,
     '    PHIGS$K_INTSTYLE_SOLID,          !front and back styles
     '    1,1,                             !front and back style indices
     '    PHIGS$K_RGB,GREEN,               !rgb front colour
     '    PHIGS$K_RGB,BLUE,               !rgb back colour
C         these next two don't seem to do anything?
     '    PHIGS$K_SHADE_NONE,PHIGS$K_SHADE_NONE, !front & back shading method
     '    PHIGS$K_LIGHTING_DIFFUSE,
     '    PHIGS$K_LIGHTING_SPECULAR, !front & back lighting
     '    0.5,0.5,0.9, !ambient,diffuse and specular reflection coeffs
     '    PHIGS$K_RGB,WHITE,100.0,1.0,!rgb specular colour,exponent,transparency
     '    1.0,0.5,0.9, !ambient,diffuse and specular reflection back faces
     '    PHIGS$K_RGB,BLUE,0.9,0.4,!rgb specular values for back faces
     '    1,4, !surface approx type and value
     '    1,4) !trim approx type and value

        CALL PHIGS$SET_GEOM_NORMAL_CAL_MODE(PHIGS$K_FORWARD)

C ***   turn on painters hidden line/surface removal              !AAY 5/4/90
        CALL PHIGS$SET_HLHSR_ID(1)
        CALL PHIGS$CLOSE_STRUCT()
C MPN 11/9/92 - put into CLOSE_WS
C         CALL PHIGS$POST_STRUCT(iw,ISVIEW,1.0)

! iw=17
      ELSE IF(iw.EQ.17.AND.IWKS(17).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        CALL STRING_TRIM(FILE00,IBEG,IEND)
        IF(DOP) THEN
          WRITE(OP_STRING,*)' opening',FILE00(IBEG:IEND),'.meta'
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        CALL GKS_OPEN_WS(17,FILE00(IBEG:IEND)//'.meta',
     '    GKS$K_GKSM_OUTPUT,ERROR,*9999)

! iw=18
      ELSE IF(iw.EQ.18.AND.IWKS(18).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.,0.48,0.288,0.48,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0AR'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=18,file='fem.w18',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(18,18,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.,0.48*DISP,0.288*DISP,0.48*DISP,ERROR,
     '      *9999)
        ENDIF
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.3,0.7,ERROR,*9999)

! iw=19
      ELSE IF(iw.EQ.19.AND.IWKS(19).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.,0.48,0.288,0.48,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0AS'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=19,file='fem.w19',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(19,19,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,DISP,1.192*DISP,0.51*DISP,0.99*DISP,
     '        ERROR,*9999)
          ENDIF
        ENDIF
        IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          CALL GKS_SWKWN(iw,0.3,0.7,0.,1.,ERROR,*9999)
        ENDIF
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)

! iw=20
      ELSE IF(iw.EQ.20.AND.IWKS(20).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.5,0.99,0.31,1.0,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0AT'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=20,file='Diameter',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(20,20,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.5*XDISP,0.99*XDISP,YDISP-0.69*XDISP,
     '        YDISP,ERROR,*9999)
          ENDIF
        ENDIF
        IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
          CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
          CALL SET_POLYLINE_REP(iw,ERROR,*9999)
          CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
          CALL SET_TEXT_REP(iw,ERROR,*9999)
          CALL GKS_SVP(iw,0.5-0.49/0.69/2.,0.5+0.49/0.69/2.,0.,1.,
     '      ERROR,*9999)
          CALL GKS_SWKWN(iw,0.5-0.49/0.69/2.,0.5+0.49/0.69/2.,0.,1.,
     '      ERROR,*9999)
        ENDIF
        CALL GKS_SWN(iw,-0.14,1.1,-3.2,3.2,ERROR,*9999)

! iw=21
      ELSE IF(iw.EQ.21.AND.IWKS(21).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.5,0.99,0.77,1.0,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0AU'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=21,file='fem.w21',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(21,21,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.5*XDISP,0.99*XDISP,YDISP-0.23*XDISP,
     '        YDISP,ERROR,*9999)
          ENDIF
        ENDIF
        IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
          CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
          CALL SET_POLYLINE_REP(iw,ERROR,*9999)
          CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
          CALL SET_TEXT_REP(iw,ERROR,*9999)
          CALL GKS_SWKWN(iw,0.,1.,0.5-0.23/0.49/2.,0.5+0.23/0.49/2.,
     '      ERROR,*9999)
        ENDIF
        CALL GKS_SWN(iw,-0.14,1.1,-2.3,2.3,ERROR,*9999)

! iw=22
      ELSE IF(iw.EQ.22.AND.IWKS(22).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.5,0.99,0.51,0.74,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0AV'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=22,file='fem.w22',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(22,22,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.5*XDISP,0.99*XDISP,YDISP-0.49*XDISP,
     '        YDISP-0.26*XDISP,ERROR,*9999)
          ENDIF
        ENDIF
        IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
          CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
          CALL SET_POLYLINE_REP(iw,ERROR,*9999)
          CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
          CALL SET_TEXT_REP(iw,ERROR,*9999)
          CALL GKS_SWKWN(iw,0.,1.,0.5-0.23/0.49/2.,0.5+0.23/0.49/2.,
     '      ERROR,*9999)
        ENDIF
        CALL GKS_SWN(iw,-0.14,1.1,-2.3,2.3,ERROR,*9999)

! iw=23
      ELSE IF(iw.EQ.23.AND.IWKS(23).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.5,0.99,0.25,0.48,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0AW'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=23,file='fem.w23',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(23,23,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.5*XDISP,0.99*XDISP,YDISP-0.75*XDISP,
     '        YDISP-0.52*XDISP,ERROR,*9999)
          ENDIF
        ENDIF
        IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
          CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
          CALL SET_POLYLINE_REP(iw,ERROR,*9999)
          CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
          CALL SET_TEXT_REP(iw,ERROR,*9999)
          CALL GKS_SWKWN(iw,0.,1.,0.5-0.23/0.49/2.,0.5+0.23/0.49/2.,
     '      ERROR,*9999)
        ENDIF
        CALL GKS_SWN(iw,-0.14,1.1,-2.3,2.3,ERROR,*9999)

! iw=31
      ELSE IF(iw.EQ.31.AND.IWKS(31).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.5,0.99,0.75,1.0,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0BE'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=31,file='fem.w31',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(31,31,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.5*XDISP,0.99*XDISP,YDISP-0.23*XDISP,
     '        YDISP,ERROR,*9999)
          ENDIF
        ENDIF
        IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
          CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
          CALL SET_POLYLINE_REP(iw,ERROR,*9999)
          CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
          CALL SET_TEXT_REP(iw,ERROR,*9999)
          CALL GKS_SWKWN(iw,0.,1.,0.5-0.23/0.49/2.,0.5+0.23/0.49/2.,
     '      ERROR,*9999)
        ENDIF
        CALL GKS_SWN(iw,-0.2,1.1,-2.3,2.3,ERROR,*9999)

! iw=32
      ELSE IF(iw.EQ.32.AND.IWKS(32).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.5,0.99,0.51,0.74,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0BF'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=32,file='fem.w32',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(32,32,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.5*XDISP,0.99*XDISP,YDISP-0.49*XDISP,
     '        YDISP-0.26*XDISP,ERROR,*9999)
          ENDIF
        ENDIF
        IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
          CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
          CALL SET_POLYLINE_REP(iw,ERROR,*9999)
          CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
          CALL SET_TEXT_REP(iw,ERROR,*9999)
          CALL GKS_SWKWN(iw,0.,1.,0.5-0.23/0.49/2.,0.5+0.23/0.49/2.,
     '      ERROR,*9999)
        ENDIF
        CALL GKS_SWN(iw,-0.2,1.1,-2.3,2.3,ERROR,*9999)

! iw=33
      ELSE IF(iw.EQ.33.AND.IWKS(33).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.5,0.99,0.77,1.0,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0BG'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=33,file='fem.w33',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(33,33,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.5*XDISP,0.99*XDISP,YDISP-0.23*XDISP,
     '        YDISP,ERROR,*9999)
          ENDIF
        ENDIF
        IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
          CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
          CALL SET_POLYLINE_REP(iw,ERROR,*9999)
          CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
          CALL SET_TEXT_REP(iw,ERROR,*9999)
          CALL GKS_SWKWN(iw,0.,1.,0.5-0.23/0.49/2.,0.5+0.23/0.49/2.,
     '      ERROR,*9999)
        ENDIF
        CALL GKS_SWN(iw,-0.2,1.1,-2.3,2.3,ERROR,*9999)

! iw=34
      ELSE IF(iw.EQ.34.AND.IWKS(34).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.4,1.0,0.0,1.0,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0BH'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=34,file='fem.w34',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(34,34,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.4*XDISP,XDISP,0.0,YDISP,ERROR,*9999)
          ENDIF
        ENDIF
        IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
          CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
          CALL SET_POLYLINE_REP(iw,ERROR,*9999)
          CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
          CALL SET_TEXT_REP(iw,ERROR,*9999)
          CALL GKS_SVP  (iw,0.,0.6*XDISP/YDISP,0.,1.,ERROR,*9999)
          CALL GKS_SWKWN(iw,0.,0.6*XDISP/YDISP,0.,1.,ERROR,*9999)
          CALL GKS_SWKVP(iw,0.4*XDISP,XDISP,0.0,YDISP,ERROR,*9999)
        ENDIF
        CALL GKS_SWN(iw,0.,30.,0.,430.,ERROR,*9999)

! iw=35
      ELSE IF(iw.EQ.35.AND.IWKS(35).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.0,0.39,0.61,1.0,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0BI'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=35,file='fem.w35',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(35,35,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.0,0.39*XDISP,YDISP-0.39*XDISP,YDISP,
     '        ERROR,*9999)
          ENDIF
        ENDIF
        IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
          CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
          CALL SET_POLYLINE_REP(iw,ERROR,*9999)
          CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
          CALL SET_TEXT_REP(iw,ERROR,*9999)
          CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)
        ENDIF
        CALL GKS_SWN(iw,0.,10.,0.,10.,ERROR,*9999)

! iw=36
      ELSE IF(iw.EQ.36.AND.IWKS(36).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0BJ'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
            CALL GKS_SWKVP(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=36,file='fem.w36',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(36,36,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
          ENDIF
        ENDIF
        IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
          CALL SET_TEXT_REP(iw,ERROR,*9999)
          CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)
        ENDIF
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)

! iw=40
      ELSE IF(iw.EQ.40.AND.IWKS(40).EQ.0) THEN
        IWKS(iw)=1 !indicates graphics output window
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.10,0.50,0.40,0.80,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0BN'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=40,file='fem.w40',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(40,40,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.10*XDISP,0.50*XDISP,
     '         YDISP-0.60*XDISP,YDISP-0.20*XDISP,ERROR,*9999)
          ENDIF
        ENDIF
        IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
          CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
          CALL SET_POLYLINE_REP(iw,ERROR,*9999)
          CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
          CALL SET_TEXT_REP(iw,ERROR,*9999)
          CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)
        ENDIF
        CALL GKS_SWN(iw,-1.1,1.1,-1.1,1.1,ERROR,*9999)

! iw=41
      ELSE IF(iw.EQ.41.AND.IWKS(41).EQ.0) THEN
        IWKS(iw)=1 !indicates graphics output window
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.51,0.91,0.40,0.80,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0BO'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=41,file='fem.w41',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(41,41,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.51*XDISP,0.91*XDISP,
     '         YDISP-0.60*XDISP,YDISP-0.20*XDISP,ERROR,*9999)
          ENDIF
        ENDIF
        IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
          CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
          CALL SET_POLYLINE_REP(iw,ERROR,*9999)
          CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
          CALL SET_TEXT_REP(iw,ERROR,*9999)
          CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)
        ENDIF
        CALL GKS_SWN(iw,-1.1,1.1,-1.1,1.1,ERROR,*9999)

! iw=42
      ELSE IF(iw.EQ.42.AND.IWKS(42).EQ.0) THEN
        IWKS(iw)=1 !indicates graphics output window
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.52,0.92,0.39,0.79,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0BP'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=42,file='fem.w42',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(42,42,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.52*XDISP,0.92*XDISP,
     '       YDISP-0.61*XDISP,YDISP-0.21*XDISP,ERROR,*9999)
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
        CALL SET_POLYLINE_REP(iw,ERROR,*9999)
        CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,-1.1,1.1,-1.1,1.1,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=43
      ELSE IF(iw.EQ.43.AND.IWKS(43).EQ.0) THEN
        IWKS(iw)=1 !indicates graphics output window
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.53,0.93,0.38,0.78,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0BQ'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=43,file='fem.w43',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(43,43,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.53*XDISP,0.93*XDISP,
     '       YDISP-0.62*XDISP,YDISP-0.22*XDISP,ERROR,*9999)
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
        CALL SET_POLYLINE_REP(iw,ERROR,*9999)
        CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,-1.1,1.1,-1.1,1.1,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=44
      ELSE IF(iw.EQ.44.AND.IWKS(44).EQ.0) THEN
        IWKS(iw)=1 !indicates graphics output window
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.20,0.69,0.44,0.67,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0BR'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=44,file='fem.w44',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(44,44,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.20*XDISP,0.69*XDISP,YDISP-0.56*XDISP,
     '      YDISP-0.33*XDISP,ERROR,*9999)
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
        CALL SET_POLYLINE_REP(iw,ERROR,*9999)
        CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN  (iw,-0.1,1.15,-0.1,1.1,ERROR,*9999)
        CALL GKS_SVP  (iw,0.,1.,0.5-0.23/0.49/2.,0.5+0.23/0.49/2.,
     '    ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.5-0.23/0.49/2.,0.5+0.23/0.49/2.,
     '    ERROR,*9999)

! iw=45
      ELSE IF(iw.EQ.45.AND.IWKS(45).EQ.0) THEN
        IWKS(iw)=1 !indicates graphics output window
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.19,0.68,0.54,0.77,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0BS'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=45,file='fem.w45',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(45,45,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.19*XDISP,0.68*XDISP,YDISP-0.46*XDISP,
     '      YDISP-0.23*XDISP,ERROR,*9999)
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
        CALL SET_POLYLINE_REP(iw,ERROR,*9999)
        CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN  (iw,-0.1,1.15,-0.1,1.1,ERROR,*9999)
        CALL GKS_SVP  (iw,0.,1.,0.5-0.23/0.49/2.,0.5+0.23/0.49/2.,
     '    ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.5-0.23/0.49/2.,0.5+0.23/0.49/2.,
     '    ERROR,*9999)

! iw=46
      ELSE IF(iw.EQ.46.AND.IWKS(46).EQ.0) THEN
        IWKS(iw)=1 !indicates graphics output window
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.18,0.67,0.64,0.87,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0BT'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=46,file='fem.w46',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(46,46,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.18*XDISP,0.67*XDISP,YDISP-0.36*XDISP,
     '      YDISP-0.13*XDISP,ERROR,*9999)
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
        CALL SET_POLYLINE_REP(iw,ERROR,*9999)
        CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN  (iw,-0.1,1.15,-0.1,1.1,ERROR,*9999)
        CALL GKS_SVP  (iw,0.,1.,0.5-0.23/0.49/2.,0.5+0.23/0.49/2.,
     '    ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.5-0.23/0.49/2.,0.5+0.23/0.49/2.,
     '    ERROR,*9999)

! iw=47
      ELSE IF(iw.EQ.47.AND.IWKS(47).EQ.0) THEN
        IWKS(iw)=1 !indicates graphics output window
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.17,0.66,0.74,0.97,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0BU'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=47,file='fem.w47',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(47,47,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.17*XDISP,0.66*XDISP,YDISP-0.26*XDISP,
     '      YDISP-0.03*XDISP,ERROR,*9999)
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
        CALL SET_POLYLINE_REP(iw,ERROR,*9999)
        CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN  (iw,-0.1,1.15,-0.1,1.1,ERROR,*9999)
        CALL GKS_SVP  (iw,0.,1.,0.5-0.23/0.49/2.,0.5+0.23/0.49/2.,
     '    ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.5-0.23/0.49/2.,0.5+0.23/0.49/2.,
     '    ERROR,*9999)

! iw=50
      ELSE IF(iw.EQ.50.AND.IWKS(50).EQ.0) THEN
        IWKS(iw)=1 !indicates graphics output window
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.,0.49,0.35,1.0,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0BX'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=50,file='fem.w50',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(50,50,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.,0.49*XDISP,YDISP-0.49*(4./3.)*XDISP,
     '      YDISP,ERROR,*9999)
        ENDIF
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN  (iw,0.0,0.75,-1.0,1.0,ERROR,*9999)
        CALL GKS_SVP  (iw,0.125,0.875,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.125,0.875,0.0,1.0,ERROR,*9999)

! iw=51
      ELSE IF(iw.EQ.51.AND.IWKS(51).EQ.0) THEN
        IWKS(iw)=1 !indicates graphics output window
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.50,0.99,0.51,1.0,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0BY'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=51,file='fem.w51',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(51,51,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.50*XDISP,0.99*XDISP,
     '       YDISP-0.49*XDISP,YDISP-0.00*XDISP,ERROR,*9999)
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
        CALL SET_POLYLINE_REP(iw,ERROR,*9999)
        CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,-0.1,1.1,-0.1,1.1,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)

! iw=55
      ELSE IF(iw.EQ.55.AND.IWKS(55).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0CC'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=55,file='fem.w55',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(55,55,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
        ENDIF
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=56
      ELSE IF(iw.EQ.56.AND.IWKS(56).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0CD'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=56,file='fem.w56',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(56,56,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
        ENDIF
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=60
      ELSE IF(iw.EQ.60.AND.IWKS(60).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.0,0.63,0.51,1.0,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0CH'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=60,file='fem.w60',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(60,60,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.,0.63*XDISP,YDISP-0.49*XDISP,YDISP,
     '      ERROR,*9999)
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
        CALL SET_POLYLINE_REP(iw,ERROR,*9999)
        CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=61
      ELSE IF(iw.EQ.61.AND.IWKS(61).EQ.0) THEN
        IWKS(iw)=1 !indicates graphics output window
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.03,0.52,0.48,0.97,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0CI'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=61,file='fem.w61',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(61,61,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.03*XDISP,0.52*XDISP,
     '       YDISP-0.52*XDISP,YDISP-0.03*XDISP,ERROR,*9999)
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
        CALL SET_POLYLINE_REP(iw,ERROR,*9999)
        CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,-0.1,1.1,-0.1,1.1,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)

! iw=62
      ELSE IF(iw.EQ.62.AND.IWKS(62).EQ.0) THEN
        IWKS(iw)=1 !indicates graphics output window
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.06,0.55,0.45,0.94,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0CJ'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=62,file='fem.w62',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(62,62,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.06*XDISP,0.55*XDISP,
     '       YDISP-0.55*XDISP,YDISP-0.06*XDISP,ERROR,*9999)
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
        CALL SET_POLYLINE_REP(iw,ERROR,*9999)
        CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
C ***   Set up fill area patterns to be 16 shades of grey
        IF(COLOUR_WS) THEN
          CALL GKS_SFAR(iw,1,GSOLID,14,1,ERROR,*9999)
          CALL GKS_SFAR(iw,16,GSOLID,14,0,ERROR,*9999)
          DO ipat=2,15
            CALL GKS_SFAR(iw,ipat,GPATTR,12+(17-ipat),1,ERROR,*9999)
          ENDDO
        ENDIF
        CALL GKS_SWN(iw,-0.1,1.1,-0.1,1.1,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=63
      ELSE IF(iw.EQ.63.AND.IWKS(63).EQ.0) THEN
        IWKS(iw)=1 !indicates graphics output window
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.09,0.58,0.42,0.91,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0CK'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=63,file='fem.w63',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(63,63,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.09*XDISP,0.58*XDISP,
     '       YDISP-0.58*XDISP,YDISP-0.09*XDISP,ERROR,*9999)
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
        CALL SET_POLYLINE_REP(iw,ERROR,*9999)
        CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,-0.1,1.1,-0.1,1.1,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=64
      ELSE IF(iw.EQ.64.AND.IWKS(64).EQ.0) THEN
        IWKS(iw)=1 !indicates graphics output window
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.12,0.61,0.39,0.88,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0CL'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=64,file='fem.w64',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(64,64,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.12*XDISP,0.61*XDISP,
     '       YDISP-0.61*XDISP,YDISP-0.12*XDISP,ERROR,*9999)
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
        CALL SET_POLYLINE_REP(iw,ERROR,*9999)
        CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,-5.5,5.5,-5.5,5.5,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=65
      ELSE IF(iw.EQ.65.AND.IWKS(65).EQ.0) THEN
        IWKS(iw)=1 !indicates graphics output window
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.15,0.64,0.36,0.85,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0CM'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=65,file='fem.w65',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(65,65,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.15*XDISP,0.64*XDISP,
     '       YDISP-0.64*XDISP,YDISP-0.15*XDISP,ERROR,*9999)
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
        CALL SET_POLYLINE_REP(iw,ERROR,*9999)
        CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,-0.1,1.1,-0.1,1.1,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=66
      ELSE IF(iw.EQ.66.AND.IWKS(66).EQ.0) THEN
        IWKS(iw)=1 !indicates graphics output window
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.18,0.67,0.33,0.82,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0CN'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=66,file='fem.w66',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(66,66,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.18*XDISP,0.67*XDISP,
     '       YDISP-0.67*XDISP,YDISP-0.18*XDISP,ERROR,*9999)
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
        CALL SET_POLYLINE_REP(iw,ERROR,*9999)
        CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,-0.1,1.1,-0.1,1.1,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=67
      ELSE IF(iw.EQ.67.AND.IWKS(67).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=2 !indicates PHIGS workstation
        IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
          CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
          NODE=DISPLAY_NAME(IBEG:IEND)//'::0CI'
          CALL STRING_TRIM(NODE,IBEG,IEND)
          CALL PHIGS$OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE)
        ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          CALL PHIGS$OPEN_WS(67,'FEM.W67',41)
        ENDIF
        XNPC(1)= 0.0
        XNPC(2)= 1.0
        XNPC(3)= 0.0
        XNPC(4)= 1.0
C       XNPC(5)=-1.0   !old phigs version
C       XNPC(6)= 0.0   !old phigs version
        XNPC(5)= 0.0   !7-feb-1990
        XNPC(6)= 1.0   !7-feb-1990
        CALL PHIGS$SET_WS_WINDOW3(67,XNPC)
        XDC(1)= 0.21*XDISP
        XDC(2)= 0.70*XDISP
        XDC(3)= YDISP-0.70*XDISP
        XDC(4)= YDISP-0.21*XDISP
        CALL PHIGS$SET_WS_VIEWPORT(67,XDC)
C ***   Define transformation and view matrices
        CALL PHIG(iw,ERROR,*9999)
        CALL PHIGS$SET_VIEW_REP3(iw,1,A_ORIENT,A_MAP,NPC_CLIP,
     '    PHIGS$K_NOCLIP,PHIGS$K_NOCLIP,PHIGS$K_NOCLIP)
        NTSG=NTSG+1
        ISVIEW=NTSG
C ***   Insert global transformation matrix into structure ISVIEW
        CALL PHIGS$SET_EDIT_MODE(PHIGS$K_EDIT_INSERT)
        CALL PHIGS$OPEN_STRUCT(ISVIEW)
        CALL PHIGS$SET_GLOBAL_XFORM3(A_TRANS)
        CALL PHIGS$SET_VIEW_INDEX(1)
        CALL PHIGS$SET_CHAR_HEIGHT(0.05*DIAG/SQRT(12.))
        CALL PHIGS$CLOSE_STRUCT()
        CALL PHIGS$POST_STRUCT(iw,ISVIEW,1.0)

! iw=68
      ELSE IF(iw.EQ.68.AND.IWKS(68).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=2 !indicates PHIGS workstation
        IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
          CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
          NODE=DISPLAY_NAME(IBEG:IEND)//'::0CP'
          CALL STRING_TRIM(NODE,IBEG,IEND)
          CALL PHIGS$OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE)
        ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          CALL PHIGS$OPEN_WS(68,'FEM.W68',41)
        ENDIF
        XNPC(1)= 0.0
        XNPC(2)= 1.0
        XNPC(3)= 0.0
        XNPC(4)= 1.0
C       XNPC(5)=-1.0   !old phigs version
C       XNPC(6)= 0.0   !old phigs version
        XNPC(5)= 0.0   !7-feb-1990
        XNPC(6)= 1.0   !7-feb-1990
        CALL PHIGS$SET_WS_WINDOW3(68,XNPC)
        XDC(1)= 0.12*XDISP
        XDC(2)= 0.61*XDISP
        XDC(3)= YDISP-0.61*XDISP
        XDC(4)= YDISP-0.12*XDISP
        CALL PHIGS$SET_WS_VIEWPORT(68,XDC)
C ***   Define transformation and view matrices
        CALL PHIG(iw,ERROR,*9999)
        CALL PHIGS$SET_VIEW_REP3(iw,1,A_ORIENT,A_MAP,NPC_CLIP,
     '    PHIGS$K_NOCLIP,PHIGS$K_NOCLIP,PHIGS$K_NOCLIP)
        NTSG=NTSG+1
        ISVIEW=NTSG
C ***   Insert global transformation matrix into structure ISVIEW
        CALL PHIGS$SET_EDIT_MODE(PHIGS$K_EDIT_INSERT)
        CALL PHIGS$OPEN_STRUCT(ISVIEW)
        CALL PHIGS$SET_GLOBAL_XFORM3(A_TRANS)
        CALL PHIGS$SET_VIEW_INDEX(1)
        CALL PHIGS$SET_CHAR_HEIGHT(0.05*DIAG/SQRT(12.0))
        CALL PHIGS$CLOSE_STRUCT()
        CALL PHIGS$POST_STRUCT(iw,ISVIEW,1.0)

! iw=69
      ELSE IF(iw.EQ.69.AND.IWKS(69).EQ.0) THEN
        IWKS(iw)=1 !indicates graphics output window
        IWKG(iw)=1 !indicates graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.12,0.61,0.39,0.88,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0CQ'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=69,file='fem.w69',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(69,69,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.12*XDISP,0.61*XDISP,
     '       YDISP-0.61*XDISP,YDISP-0.12*XDISP,ERROR,*9999)
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
        CALL SET_POLYLINE_REP(iw,ERROR,*9999)
        CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,-0.1,1.1,-0.1,1.1,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=70
      ELSE IF(iw.EQ.70.AND.IWKS(70).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0CR'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=70,file='fem.w70',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(70,70,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=71
      ELSE IF(iw.EQ.71.AND.IWKS(71).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0CS'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=71,file='fem.w71',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(71,71,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=72
      ELSE IF(iw.EQ.72.AND.IWKS(72).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0CT'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=72,file='fem.w72',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(72,72,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=73
      ELSE IF(iw.EQ.73.AND.IWKS(73).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0CU'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=73,file='fem.w73',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(73,73,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=74
      ELSE IF(iw.EQ.74.AND.IWKS(74).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0CV'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=74,file='fem.w74',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(74,74,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=75
      ELSE IF(iw.EQ.75.AND.IWKS(75).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0CW'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=75,file='fem.w75',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(75,75,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=76
      ELSE IF(iw.EQ.76.AND.IWKS(76).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0CX'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=76,file='fem.w76',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(76,76,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=77
      ELSE IF(iw.EQ.77.AND.IWKS(77).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0CY'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=77,file='fem.w77',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(77,77,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=78
      ELSE IF(iw.EQ.78.AND.IWKS(78).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0CZ'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=78,file='fem.w78',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(78,78,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=79
      ELSE IF(iw.EQ.79.AND.IWKS(79).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0DA'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=79,file='fem.w79',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(79,79,GWSDEF,ERROR,*9999)
          ENDIF
          CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=81
      ELSE IF(iw.EQ.81.AND.IWKS(81).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0DC'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
            CALL GKS_SWKVP(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=81,file='fem.w81',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(81,81,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
          ENDIF
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=82
      ELSE IF(iw.EQ.82.AND.IWKS(82).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0DD'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
            CALL GKS_SWKVP(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=82,file='fem.w82',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(82,82,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
          ENDIF
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=83
      ELSE IF(iw.EQ.83.AND.IWKS(83).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0DE'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
            CALL GKS_SWKVP(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=83,file='fem.w83',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(83,83,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
          ENDIF
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=84
      ELSE IF(iw.EQ.84.AND.IWKS(84).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0DF'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
            CALL GKS_SWKVP(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=84,file='fem.w84',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(84,84,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
          ENDIF
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=85
      ELSE IF(iw.EQ.85.AND.IWKS(85).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0DG'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
            CALL GKS_SWKVP(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=85,file='fem.w85',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(85,85,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
          ENDIF
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=86
      ELSE IF(iw.EQ.86.AND.IWKS(86).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0DH'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
            CALL GKS_SWKVP(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=86,file='fem.w86',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(86,86,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
          ENDIF
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=87
      ELSE IF(iw.EQ.87.AND.IWKS(87).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0DI'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
            CALL GKS_SWKVP(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=87,file='fem.w87',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(87,87,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
          ENDIF
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=88
      ELSE IF(iw.EQ.88.AND.IWKS(88).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0DJ'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
            CALL GKS_SWKVP(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=88,file='fem.w88',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(88,88,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
          ENDIF
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=89
      ELSE IF(iw.EQ.89.AND.IWKS(89).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0DK'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
            CALL GKS_SWKVP(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=89,file='fem.w89',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(89,89,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
          ENDIF
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=91
      ELSE IF(iw.EQ.91.AND.IWKS(91).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0DM'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
C            CALL GKS_SWKVP(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
C     '        ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=91,file='fem.w91',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(91,91,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
          ENDIF
        ENDIF
        IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
          CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
          CALL SET_FILL_AREA_REP(iw,ERROR,*9999)
          CALL SET_POLYLINE_REP(iw,ERROR,*9999)
          CALL SET_POLYMARKER_REP(iw,ERROR,*9999)
          CALL SET_TEXT_REP(iw,ERROR,*9999)
          CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
          CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)
        ENDIF

! iw=92
      ELSE IF(iw.EQ.92.AND.IWKS(92).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0DN'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
            CALL GKS_SWKVP(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=92,file='fem.w92',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(92,92,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
          ENDIF
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=93
      ELSE IF(iw.EQ.93.AND.IWKS(93).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0DO'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
            CALL GKS_SWKVP(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=93,file='fem.w93',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(93,93,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
          ENDIF
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=94
      ELSE IF(iw.EQ.94.AND.IWKS(94).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0DP'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
            CALL GKS_SWKVP(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=94,file='fem.w94',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(94,94,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
          ENDIF
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=95
      ELSE IF(iw.EQ.95.AND.IWKS(95).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0DQ'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
            CALL GKS_SWKVP(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=95,file='fem.w95',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(95,95,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
          ENDIF
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=96
      ELSE IF(iw.EQ.96.AND.IWKS(96).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0DR'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
            CALL GKS_SWKVP(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=96,file='fem.w96',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(96,96,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
          ENDIF
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=97
      ELSE IF(iw.EQ.97.AND.IWKS(97).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0DS'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
            CALL GKS_SWKVP(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=97,file='fem.w97',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(97,97,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
          ENDIF
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=98
      ELSE IF(iw.EQ.98.AND.IWKS(98).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0DT'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
            CALL GKS_SWKVP(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=98,file='fem.w98',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(98,98,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
          ENDIF
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

! iw=99
      ELSE IF(iw.EQ.99.AND.IWKS(99).EQ.0) THEN
        IWKS(iw)=1 !indicates window open
        IWKG(iw)=0 !indicates non-graphics output window
        IWKT(iw)=1 !indicates GKS workstation
        IF(USE_SOCKET) THEN
          CALL X_OPWK(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
        ELSE IF(.NOT.USE_SOCKET) THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL STRING_TRIM(DISPLAY_NAME,IBEG,IEND)
            NODE=DISPLAY_NAME(IBEG:IEND)//'::0DU'
            CALL STRING_TRIM(NODE,IBEG,IEND)
            CALL GKS_OPEN_WS(iw,NODE(IBEG:IEND),WSTYPE,ERROR,*9999)
            CALL GKS_SWKVP(iw,ECAREA(1),ECAREA(2),ECAREA(3),ECAREA(4),
     '        ERROR,*9999)
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            open(unit=99,file='fem.w99',status='unknown',
     '        dispose='delete')
            CALL GKS_OPWK(99,99,GWSDEF,ERROR,*9999)
            CALL GKS_SWKVP(iw,0.,XDISP,0.,YDISP,ERROR,*9999)
          ENDIF
        ENDIF
        CALL SET_COLOUR_REP(iw,3,MAXCOLOURS,ERROR,*9999)
        CALL SET_TEXT_REP(iw,ERROR,*9999)
        CALL GKS_SWN(iw,0.0,1.0,0.0,1.0,ERROR,*9999)
        CALL GKS_SWKWN(iw,0.,1.,0.,1.,ERROR,*9999)

      ENDIF

      IF(DOP) THEN
        WRITE(OP_STRING,'('' Reset to:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' IW='',I3,'' IWKS(iw)='',I2,'
     '    //''' IWKT(iw)='',I2,'' IWKG(iw)='',I2)')
     '    IW,IWKS(iw),IWKT(iw),IWKG(iw)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      noiw=0
      DO iiw=1,99      !to record defined workstations in IWKDEF
        IF(IWKS(iiw).GT.0) THEN !iiw is defined
          noiw=noiw+1
          IWKDEF(noiw)=iiw
        ENDIF
      ENDDO
      IWKDEF(0)=noiw

      CALL EXITS('SETUP')
      RETURN
 9999 CALL ERRORS('SETUP',ERROR)
      CALL EXITS('SETUP')
      RETURN 1
      END


      SUBROUTINE STROKE(INIT,iw,INSTAT,MODE,NTPTS,ERROR,*)

C#### Subroutine: STROKE
C###  Description:
C**** Calls GKS stroke
C**** INIT specifies whether stroke is to be initialised (1:yes,2:no)
C**** iw specifies workstation number (which is also transformation no)
C**** INSTAT is returned as 1 if locate is successful, 0 otherwise
C**** MODE is input mode - 'REQUEST','SAMPLE' or 'EVENT'
C**** NECHO is 0 on entry for no echo
C**** NECHO is 3 on entry for default   prompt/echo
C**** NECHO is 4 on entry for line      prompt/echo
C**** NECHO is 5 on entry for rectangle prompt/echo
C**** NECHO is 6 on entry for digital   prompt/echo
C**** XPTS,YPTS are returned world coords of stroke

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
!     Parameter List
      INTEGER INIT,INSTAT,iw,NTPTS
      CHARACTER ERROR*(*),MODE*(*)
!     Local Variables
      INTEGER IERR,IWCNDC,LD1,NTRN,NTT
      DATA LD1/1/

      CALL ENTERS('STROKE',*9999)
      IF(INIT.EQ.1) THEN !initialise stroke device
        CALL GKS_SSKM(
     '   iw,LD1,GREQU,GECHO,ERROR,*9999) !must be in request mode to initialize
      ENDIF
      CALL GKS_QENTN(1,IERR,NTT,NTRN,ERROR,*9999)
      IF(iw.ne.NTRN) CALL GKS_SVPIP(iw,NTRN,GHIGHR,ERROR,*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,*)' iw=',iw,
     ''   call to gqentn returns ierr,ntt,ntrn:',IERR,NTT,NTRN
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      IF(MODE(1:5).EQ.'EVENT') THEN
      ELSE IF(MODE(1:6).EQ.'SAMPLE') THEN
      ELSE IF(MODE(1:7).EQ.'REQUEST') THEN
        CALL GKS_RQSK(iw,LD1,500,INSTAT,IWCNDC,NTPTS,ERROR,*9999)
      ENDIF
      IF(INSTAT.EQ.GOK) THEN
        INSTAT=1
      ELSE
        INSTAT=0
      ENDIF
      IF(DOP) THEN
        WRITE(OP_STRING,'('' IWCNDC='',I2)') IWCNDC
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('STROKE')
      RETURN

 9999 CALL ERRORS('STROKE',ERROR)
      CALL EXITS('STROKE')
      RETURN 1
      END


      SUBROUTINE SURFACE(IBUNDLE,iw,NVERTICES,ZZP,FACET_COLOURS,
     '  FACET_NORMAL,VERTEX_COLOURS,VERTEX_NORMAL,NFACETS,NVPF,
     '  IVERTICES,NVERTMX,ERROR,*)

C#### Subroutine: SURFACE
C###  Description:
C**** draws a surface on iw with index IBUNDLE. VERTICES(1..3,np) is an array
C**** containing the 3D coordinates of each vertex. If the IBUNDLE is 0 the
C**** primative will use the previously defined surface index.
C**** Only works for iw=3 at present.
C**** FACET_COLOURS  array of colour indices for each facet
C**** FACET_NORMAL   array of facet normals
C**** VERTEX_COLOURS array of vertex colour indices for each facet
C**** VERTEX_NORMAL  array of vertex normals
C**** NFACETS        number of facets
C**** NVPF           array containing #vertices for each facet
C**** IVERTICES      array of vertex indices for each facet

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:surf00.cmn'
      INCLUDE 'SYS$LIBRARY:PHIGS$DEFS.FOR'
!     Parameter List
      INTEGER IBUNDLE,IVERTICES(3,*),iw,FACET_COLOURS(*),NFACETS,
     '  NVERTICES,NVERTMX,NVPF(*),VERTEX_COLOURS(*)
      REAL FACET_NORMAL(3,*),VERTEX_NORMAL(3,*)
      REAL*8 ZZP(3,*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ACTIVE(2),DEACT(2),ICOLOUR_TYPE,IEDGEVIS(1000),
     '  IEDGEVIS_FLAG,INDEX_COLOUR,nf,nj,nvert,nvertMAX
      PARAMETER(nvertMAX=1000)
      REAL VERTICES(3,nvertMAX)  !must be the same dimension as ZZP in subr SGSURF
      DATA IEDGEVIS/1000*1/

      CALL ENTERS('SURFACE',*9999)

      CALL ASSERT(NVERTMX.EQ.nvertMAX,
     '  '>>Redimension VERTICES in subr. surface',ERROR,*9999)
      DO nj=1,3
        DO nvert=1,nvertMAX
          VERTICES(nj,nvert)=REAL(ZZP(nj,nvert))
        ENDDO
      ENDDO

      IF(NVERTICES.EQ.0)GOTO 9998

      IF(iw.EQ.1) THEN

      ELSE IF(iw.EQ.2) THEN

      ELSE IF(iw.EQ.3) THEN
C ***   set rendering model
        CALL PHIGS$SET_RENDER_COLOUR_MODEL(PHIGS$K_RGB)
C ***   set surface properties
        IF(ibundle.ne.0) CALL PHIGS$SET_INT_INDEX(IBUNDLE)
C ***   set depth cue representation
        IF(DEPTHCUE(1:2).EQ.'ON') THEN
          CALL PHIGS$SET_DEPTH_CUE_INDEX(2)
        ELSE
          CALL PHIGS$SET_DEPTH_CUE_INDEX(1)
        ENDIF
C ***   set edge representation
        IF(EDGES(1:2).EQ.'ON') THEN
          IEDGEVIS_FLAG=PHIGS$K_PE_VISIBILITY
          CALL PHIGS$SET_EDGE_INDEX(1)
          CALL PHIGS$SET_EDGE_COLOUR_INDEX(1)
        ELSE
          IEDGEVIS_FLAG=PHIGS$K_PE_NONE
          CALL PHIGS$SET_EDGE_INDEX(2)
        ENDIF
C ***   face culling
        IF(CULL(1:4).EQ.'NONE') THEN
          CALL PHIGS$SET_FACE_DIST_MODE(PHIGS$K_DIST_NO)
          CALL PHIGS$SET_FACE_CULL_MODE(PHIGS$K_CULL_NONE)
        ELSE IF(CULL(1:4).EQ.'BACK') THEN
          CALL PHIGS$SET_FACE_DIST_MODE(PHIGS$K_DIST_YES)
          CALL PHIGS$SET_FACE_CULL_MODE(PHIGS$K_CULL_BACK)
        ELSE IF(CULL(1:5).EQ.'FRONT') THEN
          CALL PHIGS$SET_FACE_DIST_MODE(PHIGS$K_DIST_YES)
          CALL PHIGS$SET_FACE_CULL_MODE(PHIGS$K_CULL_FRONT)
        ENDIF
C ***   activate the two lights
        IF(INDEX_SURF.EQ.4.OR.INDEX_SURF.EQ.5) THEN
          ACTIVE(1)=1
          ACTIVE(2)=2
          CALL PHIGS$SET_LIGHT_SRC_STATE(2,ACTIVE,0,DEACT)
        ENDIF
C ***   colour indexing
        IF(VARIABLE_TYPE(1:5).EQ.'FIELD') THEN
          INDEX_COLOUR=PHIGS$K_PFA_COLOUR
          ICOLOUR_TYPE=PHIGS$K_INDEXED_COLOUR
        ELSE IF(VARIABLE_TYPE(1:6).EQ.'STRAIN') THEN
          INDEX_COLOUR=PHIGS$K_PFA_COLOUR
          ICOLOUR_TYPE=PHIGS$K_INDEXED_COLOUR
        ELSE IF(VARIABLE_TYPE(1:8).EQ.'GEOMETRY') THEN
          INDEX_COLOUR=PHIGS$K_PFA_NONE
          ICOLOUR_TYPE=PHIGS$K_RGB
          IF(INDEX_SURF.EQ.2) THEN !white fill area
            INDEX_COLOUR=PHIGS$K_PFA_COLOUR
            ICOLOUR_TYPE=PHIGS$K_INDEXED_COLOUR
            DO nf=1,NFACETS
              FACET_COLOURS(nf)=0
            ENDDO
          ENDIF
        ENDIF

        CALL PHIGS$INDEX_POLYGONS_WITH_DATA(
     '    PHIGS$K_SHAPE_CONVEX,  !convex facets
     '    INDEX_COLOUR,          !per facet data flag
     '    PHIGS$K_PV_NONE,       !per vertex data flag
     '    IEDGEVIS_FLAG,         !per edge data flag
     '    IEDGEVIS,              !array of edge visibility flags
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
C
C CPB 29/3/92 adding surface postscript routine
C
      ELSE IF(iw.EQ.16) THEN
C ***   set rendering model
        CALL PHIGS$SET_RENDER_COLOUR_MODEL(PHIGS$K_RGB)
C ***   set surface properties
        IF(ibundle.ne.0) CALL PHIGS$SET_INT_INDEX(IBUNDLE)
C ***   set depth cue representation
        IF(DEPTHCUE(1:2).EQ.'ON') THEN
          CALL PHIGS$SET_DEPTH_CUE_INDEX(2)
        ELSE
          CALL PHIGS$SET_DEPTH_CUE_INDEX(1)
        ENDIF
C ***   set edge representation
        IF(EDGES(1:2).EQ.'ON') THEN
          IEDGEVIS_FLAG=PHIGS$K_PE_VISIBILITY
          CALL PHIGS$SET_EDGE_INDEX(1)
          CALL PHIGS$SET_EDGE_COLOUR_INDEX(1)
        ELSE
          IEDGEVIS_FLAG=PHIGS$K_PE_NONE
          CALL PHIGS$SET_EDGE_INDEX(2)
        ENDIF
C ***   face culling
        IF(CULL(1:4).EQ.'NONE') THEN
          CALL PHIGS$SET_FACE_DIST_MODE(PHIGS$K_DIST_NO)
          CALL PHIGS$SET_FACE_CULL_MODE(PHIGS$K_CULL_NONE)
        ELSE IF(CULL(1:4).EQ.'BACK') THEN
          CALL PHIGS$SET_FACE_DIST_MODE(PHIGS$K_DIST_YES)
          CALL PHIGS$SET_FACE_CULL_MODE(PHIGS$K_CULL_BACK)
        ELSE IF(CULL(1:5).EQ.'FRONT') THEN
          CALL PHIGS$SET_FACE_DIST_MODE(PHIGS$K_DIST_YES)
          CALL PHIGS$SET_FACE_CULL_MODE(PHIGS$K_CULL_FRONT)
        ENDIF
C ***   activate the two lights
        IF(INDEX_SURF.EQ.4.OR.INDEX_SURF.EQ.5) THEN
          ACTIVE(1)=1
          ACTIVE(2)=2
          CALL PHIGS$SET_LIGHT_SRC_STATE(2,ACTIVE,0,DEACT)
        ENDIF
C ***   colour indexing
        IF(VARIABLE_TYPE(1:5).EQ.'FIELD') THEN
          INDEX_COLOUR=PHIGS$K_PFA_COLOUR
          ICOLOUR_TYPE=PHIGS$K_INDEXED_COLOUR
        ELSE IF(VARIABLE_TYPE(1:6).EQ.'STRAIN') THEN
          INDEX_COLOUR=PHIGS$K_PFA_COLOUR
          ICOLOUR_TYPE=PHIGS$K_INDEXED_COLOUR
        ELSE IF(VARIABLE_TYPE(1:8).EQ.'GEOMETRY') THEN
          INDEX_COLOUR=PHIGS$K_PFA_NONE
          ICOLOUR_TYPE=PHIGS$K_RGB
          IF(INDEX_SURF.EQ.2) THEN !white fill area
            INDEX_COLOUR=PHIGS$K_PFA_COLOUR
            ICOLOUR_TYPE=PHIGS$K_INDEXED_COLOUR
            DO nf=1,NFACETS
              FACET_COLOURS(nf)=0
            ENDDO
          ENDIF
        ENDIF

        CALL PHIGS$INDEX_POLYGONS_WITH_DATA(
     '    PHIGS$K_SHAPE_CONVEX,  !convex facets
     '    INDEX_COLOUR,          !per facet data flag
     '    PHIGS$K_PV_NONE,       !per vertex data flag
     '    IEDGEVIS_FLAG,         !per edge data flag
     '    IEDGEVIS,              !array of edge visibility flags
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

      ENDIF

 9998 CALL EXITS('SURFACE')
      RETURN
 9999 CALL ERRORS('SURFACE',ERROR)
      CALL EXITS('SURFACE')
      RETURN 1
      END


      SUBROUTINE TEXT(IBUNDLE,IGEOM,iw,STRING,D_PT,ERROR,*)

C#### Subroutine: TEXT
C###  Description:
C**** Draws text on iw with non-geometric attributes bundled in
C**** IBUNDLE, and geometric attributes according to IGEOM.
C**** D_PT(1..3) contains the REAL*8 3D coords of each point.
C**** IF iw=4 coords are curvilinear coords, else rect. cart.
C**** STRING contains the text to be shown at D_PT(nj)
C**** If IBUNDLE or IGEOM is 0 the primitive will use the previously
C**** defined text index.
C**** NOTE: If iw is 15 or 16 (postscript) IBUNDLE is reset to be black.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:disp00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'cmiss$reference:map000.cmn'
      INCLUDE 'gx$path:gx.inc'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
      INCLUDE 'SYS$LIBRARY:PHIGS$DEFS.FOR'
!     Parameter List
      INTEGER IBUNDLE,IGEOM,iw
      REAL*8 D_PT(*)
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER ERR,nj
      REAL TVECS(2,3),X(3)

      CALL ENTERS('TEXT',*9999)

      IF(DOP) THEN
        WRITE(OP_STRING,'('' D_PT(nj): '',3E12.3)')
     '    (D_PT(nj),nj=1,NJT)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(iw.EQ.1) THEN
        IF(ibundle.ne.0) CALL GKS_STXI(IBUNDLE,ERROR,*9999)
        IF(igeom.ne.0) THEN
C         set up geometric attributes
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            IF(IGEOM.EQ.1) THEN
              CALL STXALN(gxCNTRE,gxCNTRE,ERR)
            ELSE IF(IGEOM.EQ.2) THEN
              CALL STXALN(gxCNTRE,gxCNTRE,ERR)
            ENDIF
            IF(err.ne.0) THEN
              ERROR='>>Error from GX'
              GOTO 9999
            ENDIF
            CALL STXHGT(0.01*DIAG,ERR)
            IF(err.ne.0) THEN
              ERROR='>>Error from GX'
              GOTO 9999
            ENDIF
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            IF(IGEOM.EQ.1) THEN
              CALL GKS_STXAL(GACENT,GAHALF,ERROR,*9999)
            ELSE IF(IGEOM.EQ.2) THEN
              CALL GKS_STXAL(GAHNOR,GAVNOR,ERROR,*9999)
            ENDIF
            CALL GKS_SCHH(0.01*DIAG,ERROR,*9999)
            CALL GKS_STXP(GRIGHT,ERROR,*9999)
            CALL GKS_SCHUP(0.0,1.0,ERROR,*9999)
          ENDIF
        ENDIF
!       IF(NJT.EQ.2) THEN !plot x against y
          CALL GKS_TX(D_PT(1),D_PT(2),STRING,ERROR,*9999)
!       ELSE IF(NJT.EQ.3) THEN !plot x against z
!         CALL GKS_TX(D_PT(1),D_PT(3),STRING,ERROR,*9999)
!       ENDIF

      ELSE IF(iw.EQ.2) THEN
        IF(ibundle.ne.0) CALL GKS_STXI(IBUNDLE,ERROR,*9999)
        IF(igeom.ne.0) THEN
C         set up geometric attributes
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            IF(IGEOM.EQ.1) THEN
              CALL STXALN(gxCNTRE,gxCNTRE,ERR)
            ELSE IF(IGEOM.EQ.2) THEN
              CALL STXALN(gxCNTRE,gxCNTRE,ERR)
            ENDIF
            IF(err.ne.0) THEN
              ERROR='>>Error from GX'
              GOTO 9999
            ENDIF
            CALL STXHGT(0.01*DIAG,ERR)
            IF(err.ne.0) THEN
              ERROR='>>Error from GX'
              GOTO 9999
            ENDIF
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            IF(IGEOM.EQ.1) THEN
              CALL GKS_STXAL(GACENT,GAHALF,ERROR,*9999)
            ELSE IF(IGEOM.EQ.2) THEN
              CALL GKS_STXAL(GAHNOR,GAVNOR,ERROR,*9999)
            ENDIF
            CALL GKS_SCHH(0.01*DIAG,ERROR,*9999)
            CALL GKS_STXP(GRIGHT,ERROR,*9999)
            CALL GKS_SCHUP(0.0,1.0,ERROR,*9999)
          ENDIF
        ENDIF
!       plot y against z
        CALL GKS_TX(D_PT(2),D_PT(3),STRING,ERROR,*9999)

      ELSE IF(iw.EQ.3) THEN
        IF(ibundle.ne.0) CALL GKS_STXI(IBUNDLE,ERROR,*9999)
        IF(igeom.ne.0) THEN
C         set up geometric attributes
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            IF(IGEOM.EQ.1) THEN
              CALL STXALN(gxCNTRE,gxCNTRE,ERR)
            ELSE IF(IGEOM.EQ.2) THEN
              CALL STXALN(gxCNTRE,gxCNTRE,ERR)
            ENDIF
            IF(err.ne.0) THEN
              ERROR='>>Error from GX'
              GOTO 9999
            ENDIF
            CALL STXHGT(0.01*DIAG,ERR)
            IF(err.ne.0) THEN
              ERROR='>>Error from GX'
              GOTO 9999
            ENDIF
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            IF(IGEOM.EQ.1) THEN
              CALL GKS_STXAL(GACENT,GAHALF,ERROR,*9999)
            ELSE IF(IGEOM.EQ.2) THEN
              CALL GKS_STXAL(GAHNOR,GAVNOR,ERROR,*9999)
            ENDIF
            CALL GKS_SCHH(0.01*DIAG,ERROR,*9999)
            CALL GKS_STXP(GRIGHT,ERROR,*9999)
            CALL GKS_SCHUP(0.0,1.0,ERROR,*9999)
          ENDIF
        ENDIF
!       plot x against z
        CALL GKS_TX(D_PT(1),D_PT(3),STRING,ERROR,*9999)

      ELSE IF(iw.EQ.3.OR.iw.EQ.67.OR.iw.EQ.68) THEN
        IF(ibundle.ne.0) CALL PHIGS$SET_PLINE_INDEX(IBUNDLE)
        IF(igeom.ne.0) THEN
C         set up geometric attributes
          IF(IGEOM.EQ.1) THEN
            TVECS(1,1)=1.0D0
            TVECS(1,2)=0.0D0
            TVECS(1,3)=0.0D0
            TVECS(2,1)=0.0D0
            TVECS(2,2)=1.0D0
            TVECS(2,3)=0.0D0
C           set character height etc.
          ENDIF
        ENDIF
        CALL PHIGS_TEXT3(D_PT,TVECS,STRING,ERROR,*9999)

      ELSE IF(iw.EQ.4) THEN
        IF(ibundle.ne.0) CALL GKS_STXI(IBUNDLE,ERROR,*9999)
        IF(igeom.ne.0) THEN
C         set up geometric attributes
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            IF(IGEOM.EQ.1) THEN
              CALL STXALN(gxCNTRE,MCNTRE,ERR)
            ELSE IF(IGEOM.EQ.2) THEN
              CALL STXALN(gxCNTRE,MCNTRE,ERR)
            ELSE IF(IGEOM.EQ.3) THEN
              CALL STXALN(gxLEFT,MCNTRE,ERR)
            ENDIF
            IF(err.ne.0) THEN
              ERROR='>>Error from GX'
              GOTO 9999
            ENDIF
            IF(IGEOM.EQ.1) THEN
              CALL STXHGT(0.03,ERR)
            ELSE IF(IGEOM.EQ.2) THEN
              CALL STXHGT(0.02,ERR)
            ELSE IF(IGEOM.EQ.3) THEN
              CALL STXHGT(0.03,ERR)
            ENDIF
            IF(err.ne.0) THEN
              ERROR='>>Error from GX'
              GOTO 9999
            ENDIF
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            IF(IGEOM.EQ.1) THEN
              CALL GKS_STXAL(GACENT,GAHALF,ERROR,*9999)
              CALL GKS_SCHH(0.03,ERROR,*9999)
            ELSE IF(IGEOM.EQ.2) THEN
              CALL GKS_STXAL(GAHNOR,GAVNOR,ERROR,*9999)
              CALL GKS_SCHH(0.02,ERROR,*9999)
            ELSE IF(IGEOM.EQ.3) THEN
              CALL GKS_STXAL(GALEFT,GAHALF,ERROR,*9999)
              CALL GKS_SCHH(0.03,ERROR,*9999)
            ENDIF
            CALL GKS_STXP(GRIGHT,ERROR,*9999)
            CALL GKS_SCHUP(0.0,1.0,ERROR,*9999)
          ENDIF
        ENDIF
        IF(PROJEC(1:11).EQ.'RECTANGULAR') THEN !points x and y
        ELSE IF(PROJEC(1:2).EQ.'XI') THEN !assume pts are in Xi coords
          D_PT(1)=-1.0D0+2.0D0*(DBLE(MXI1-1)+D_PT(1))/MAX_XI
          D_PT(2)=-1.0D0+2.0D0*(DBLE(MXI2-1)+D_PT(2))/MAX_XI
        ELSE !plot Y against Z - assume been transf.d to polar coords
          CALL MAP4(1,D_PT,D_PT,ERROR,*9999)
        ENDIF
        CALL GKS_TX(D_PT(1),D_PT(2),STRING,ERROR,*9999)

      ELSE IF(iw.EQ.8) THEN  !help window
        IF(ibundle.ne.0) CALL GKS_STXI(IBUNDLE,ERROR,*9999)
        IF(igeom.ne.0) THEN
C         set up geometric attributes
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            IF(IGEOM.EQ.1) THEN
              CALL STXALN(gxCNTRE,gxCNTRE,ERR)
            ELSE IF(IGEOM.EQ.2) THEN
              CALL STXALN(gxCNTRE,gxCNTRE,ERR)
            ENDIF
            IF(err.ne.0) THEN
              ERROR='>>Error from GX'
              GOTO 9999
            ENDIF
            CALL STXHGT(0.01*DIAG,ERR)
            IF(err.ne.0) THEN
              ERROR='>>Error from GX'
              GOTO 9999
            ENDIF
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            IF(IGEOM.EQ.1) THEN
              CALL GKS_STXAL(GACENT,GAHALF,ERROR,*9999)
            ELSE IF(IGEOM.EQ.2) THEN
              CALL GKS_STXAL(GAHNOR,GAVNOR,ERROR,*9999)
            ENDIF
            CALL GKS_SCHH(0.01,ERROR,*9999)
            CALL GKS_STXP(GRIGHT,ERROR,*9999)
            CALL GKS_SCHUP(0.0,1.0,ERROR,*9999)
          ENDIF
        ENDIF
!       plot x against y
        CALL GKS_TX(D_PT(1),D_PT(2),STRING,ERROR,*9999)

      ELSE IF(iw.EQ.10) THEN !history profiles
        IF(ibundle.ne.0) CALL GKS_STXI(IBUNDLE,ERROR,*9999)
        IF(igeom.ne.0) THEN
C         set up geometric attributes
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            IF(IGEOM.EQ.1) THEN
              CALL STXALN(gxCNTRE,gxCNTRE,ERR)
            ELSE IF(IGEOM.EQ.2) THEN
              CALL STXALN(gxCNTRE,gxCNTRE,ERR)
            ENDIF
            IF(err.ne.0) THEN
              ERROR='>>Error from GX'
              GOTO 9999
            ENDIF
            CALL STXHGT(0.01*DIAG,ERR)
            IF(err.ne.0) THEN
              ERROR='>>Error from GX'
              GOTO 9999
            ENDIF
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            IF(IGEOM.EQ.1) THEN
              CALL GKS_STXAL(GACENT,GAHALF,ERROR,*9999)
            ELSE IF(IGEOM.EQ.2) THEN
              CALL GKS_STXAL(GACENT,GATOP,ERROR,*9999)
            ENDIF
            CALL GKS_SCHH(0.01,ERROR,*9999)
            CALL GKS_STXP(GRIGHT,ERROR,*9999)
            CALL GKS_SCHUP(0.0,1.0,ERROR,*9999)
          ENDIF
        ENDIF
!       plot x against y
        CALL GKS_TX(D_PT(1),D_PT(2),STRING,ERROR,*9999)

      ELSE IF(iw.EQ.11) THEN !pressure section
        IF(ibundle.ne.0) CALL GKS_STXI(IBUNDLE,ERROR,*9999)
        IF(igeom.ne.0) THEN
C         set up geometric attributes
          IF(IGEOM.EQ.1) THEN
            CALL GKS_STXAL(GACENT,GAHALF,ERROR,*9999)
          ELSE IF(IGEOM.EQ.2) THEN
            CALL GKS_STXAL(GARITE,GAHALF,ERROR,*9999)
          ENDIF
          CALL GKS_SCHH(0.01,ERROR,*9999)
          IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            CALL GKS_STXP(GRIGHT,ERROR,*9999)
            CALL GKS_SCHUP(0.0,1.0,ERROR,*9999)
          ENDIF
        ENDIF
!       plot x against y
        CALL GKS_TX(D_PT(1),D_PT(2),STRING,ERROR,*9999)

      ELSE IF(iw.EQ.12) THEN !fibre section
        IF(ibundle.ne.0) CALL GKS_STXI(IBUNDLE,ERROR,*9999)
        IF(igeom.ne.0) THEN
C         set up geometric attributes
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            IF(IGEOM.EQ.1) THEN
              CALL STXALN(gxCNTRE,gxCNTRE,ERR)
            ELSE IF(IGEOM.EQ.2) THEN
              CALL STXALN(gxCNTRE,gxCNTRE,ERR)
            ELSE IF(IGEOM.EQ.3) THEN
              CALL STXALN(gxCNTRE,gxCNTRE,ERR)
            ENDIF
            IF(err.ne.0) THEN
              ERROR='>>Error from GX'
              GOTO 9999
            ENDIF
            CALL STXHGT(0.01*DIAG,ERR)
            IF(err.ne.0) THEN
              ERROR='>>Error from GX'
              GOTO 9999
            ENDIF
          ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            IF(IGEOM.EQ.1) THEN
              CALL GKS_STXAL(GACENT,GAHALF,ERROR,*9999)
            ELSE IF(IGEOM.EQ.2) THEN
              CALL GKS_STXAL(GARITE,GAHALF,ERROR,*9999)
            ELSE IF(IGEOM.EQ.3) THEN
              CALL GKS_STXAL(GACENT,GATOP,ERROR,*9999)
            ENDIF
            CALL GKS_SCHH(0.1,ERROR,*9999)
            CALL GKS_STXP(GRIGHT,ERROR,*9999)
            CALL GKS_SCHUP(0.0,1.0,ERROR,*9999)
          ENDIF
        ENDIF
!       plot x against y
        CALL GKS_TX(D_PT(1),D_PT(2),STRING,ERROR,*9999)

      ELSE IF(iw.EQ.13) THEN !sheet plot
        IF(ibundle.ne.0) CALL GKS_SPLI(IBUNDLE,ERROR,*9999)
        X(1)=D_PT(1)
        X(2)=D_PT(2)
        CALL GKS_TX(X(1),X(2),STRING,ERROR,*9999)

      ELSE IF(iw.EQ.15) THEN !GKS postscript
        IF(ibundle.ne.0) THEN
C CPB 1/11/92 removed as using colour postscript
C          IF(IBUNDLE.GT.8) THEN !reset to black for postscript
C            IBUNDLE=IBUNDLE-8
C            IF(DOP) WRITE(IOOP,'('' Reset IBUNDLE='',I2)') IBUNDLE
C          ENDIF
            CALL GKS_SCHUP(0.0,1.0,ERROR,*9999)
            CALL GKS_SCHH(0.01*DIAG,ERROR,*9999)
C CPB 1/11/92 Changed call from GKS_SPLI to GKS_STXI
          CALL GKS_STXI(IBUNDLE,ERROR,*9999)
        ENDIF
        IF(igeom.ne.0) THEN
C         set up geometric attributes
          IF(IGEOM.EQ.1) THEN
            CALL GKS_STXAL(GACENT,GAHALF,ERROR,*9999)
          ELSE IF(IGEOM.EQ.2) THEN
            CALL GKS_STXAL(GAHNOR,GAVNOR,ERROR,*9999)
          ENDIF
          CALL GKS_SCHH(0.01*DIAG,ERROR,*9999)
          IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            CALL GKS_STXP(GRIGHT,ERROR,*9999)
            CALL GKS_SCHUP(0.0,1.0,ERROR,*9999)
          ENDIF
        ENDIF
        IF(NJT.EQ.2) THEN !plot x against y
          CALL GKS_TX(D_PT(1),D_PT(2),STRING,ERROR,*9999)
        ELSE IF(NJT.EQ.3) THEN !plot x against z
          CALL GKS_TX(D_PT(1),D_PT(3),STRING,ERROR,*9999)
        ENDIF

      ELSE IF(iw.EQ.16) THEN !PHIGS postscript
        IF(ibundle.ne.0) THEN
          IF(IBUNDLE.GT.8) THEN !reset to black for postscript
            IBUNDLE=IBUNDLE-8
            IF(DOP) THEN
              WRITE(OP_STRING,'('' Reset IBUNDLE='',I2)') IBUNDLE
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF
          CALL PHIGS$SET_PLINE_INDEX(IBUNDLE)
        ENDIF
        IF(igeom.ne.0) THEN
C         set up geometric attributes
          IF(IGEOM.EQ.1) THEN
            TVECS(1,1)=1.0
            TVECS(1,2)=0.0
            TVECS(1,3)=0.0
            TVECS(2,1)=0.0
            TVECS(2,2)=1.0
            TVECS(2,3)=0.0
C           set character height etc.
          ENDIF
        ENDIF
        CALL PHIGS_TEXT3(D_PT,TVECS,STRING,ERROR,*9999)

      ELSE IF(iw.EQ.21.OR.iw.EQ.22.OR.iw.EQ.23 !oxsoft plots
     '  .OR.iw.EQ.31.OR.iw.EQ.32) THEN !stress plots
        IF(ibundle.ne.0) CALL GKS_STXI(IBUNDLE,ERROR,*9999)
        IF(igeom.ne.0) THEN
C         set up geometric attributes
          IF(IGEOM.EQ.1) THEN
            CALL GKS_STXAL(GACENT,GAHALF,ERROR,*9999)
          ELSE IF(IGEOM.EQ.2) THEN
            CALL GKS_STXAL(GACENT,GATOP,ERROR,*9999)
          ELSE IF(IGEOM.EQ.3) THEN
            CALL GKS_STXAL(GARITE,GAHALF,ERROR,*9999)
          ENDIF
          CALL GKS_SCHH(0.01,ERROR,*9999)
          IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            CALL GKS_STXP(GRIGHT,ERROR,*9999)
            CALL GKS_SCHUP(0.0,1.0,ERROR,*9999)
          ENDIF
        ENDIF
!       plot x against y
        CALL GKS_TX(D_PT(1),D_PT(2),STRING,ERROR,*9999)

      ELSE IF(iw.EQ.33) THEN !fibre profiles
        IF(ibundle.ne.0) CALL GKS_STXI(IBUNDLE,ERROR,*9999)
        IF(igeom.ne.0) THEN
C         set up geometric attributes
          IF(IGEOM.EQ.1) THEN
            CALL GKS_STXAL(GACENT,GAHALF,ERROR,*9999)
          ELSE IF(IGEOM.EQ.2) THEN
            CALL GKS_STXAL(GARITE,GAHALF,ERROR,*9999)
          ENDIF
          CALL GKS_SCHH(0.01,ERROR,*9999)
          IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            CALL GKS_STXP(GRIGHT,ERROR,*9999)
            CALL GKS_SCHUP(0.0,1.0,ERROR,*9999)
          ENDIF
        ENDIF
!       plot x against y
        CALL GKS_TX(D_PT(1),D_PT(2),STRING,ERROR,*9999)

      ELSE IF(iw.EQ.34) THEN !signal display
        IF(ibundle.ne.0) CALL GKS_STXI(IBUNDLE,ERROR,*9999)
        IF(igeom.ne.0) THEN
C         set up geometric attributes
          IF(IGEOM.EQ.1) THEN
            CALL GKS_STXAL(GARITE,GAHALF,ERROR,*9999)
            CALL GKS_SCHH(0.8,ERROR,*9999)
          ELSE IF(IGEOM.EQ.2) THEN
            CALL GKS_STXAL(GACENT,GAHALF,ERROR,*9999)
            CALL GKS_SCHH(1.0,ERROR,*9999)
          ENDIF
          IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            CALL GKS_STXP(GRIGHT,ERROR,*9999)
            CALL GKS_SCHUP(0.0,1.0,ERROR,*9999)
          ENDIF
        ENDIF
!       plot x against y
        CALL GKS_TX(D_PT(1),D_PT(2),STRING,ERROR,*9999)

      ELSE IF(iw.EQ.35) THEN !trace display
        IF(ibundle.ne.0) CALL GKS_STXI(IBUNDLE,ERROR,*9999)
        IF(igeom.ne.0) THEN
C         set up geometric attributes
          IF(IGEOM.EQ.1) THEN
            CALL GKS_STXAL(GARITE,GAHALF,ERROR,*9999)
            CALL GKS_SCHH(0.1,ERROR,*9999)
          ELSE IF(IGEOM.EQ.2) THEN
            CALL GKS_STXAL(GACENT,GAHALF,ERROR,*9999)
            CALL GKS_SCHH(0.01,ERROR,*9999)
          ENDIF
          IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
            CALL GKS_STXP(GRIGHT,ERROR,*9999)
            CALL GKS_SCHUP(0.0,1.0,ERROR,*9999)
          ENDIF
        ENDIF
!       plot x against y
        CALL GKS_TX(D_PT(1),D_PT(2),STRING,ERROR,*9999)

      ELSE IF(iw.EQ.40.OR.iw.EQ.41.OR.iw.EQ.42.OR.iw.EQ.43.OR.iw.EQ.44
     '  .OR.iw.EQ.45.OR.iw.EQ.46.OR.iw.EQ.47
     '  .OR.iw.EQ.50.OR.iw.EQ.51
     '  .OR.iw.EQ.60.OR.iw.EQ.61.OR.iw.EQ.62.OR.iw.EQ.63.OR.iw.EQ.64
     '  .OR.iw.EQ.65.OR.iw.EQ.66.OR.iw.EQ.68.OR.iw.EQ.69) THEN !plots
        IF(ibundle.ne.0) CALL GKS_STXI(IBUNDLE,ERROR,*9999)
        IF(igeom.ne.0) THEN
C         set up geometric attributes
          IF(IGEOM.EQ.1) THEN
            CALL GKS_STXAL(GALEFT,GAHALF,ERROR,*9999)
            CALL GKS_SCHH(0.02,ERROR,*9999)
            IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
              CALL GKS_STXP(GRIGHT,ERROR,*9999)
              CALL GKS_SCHUP(0.0,1.0,ERROR,*9999)
            ENDIF
          ELSE IF(IGEOM.EQ.2) THEN
            CALL GKS_STXAL(GACENT,GAHALF,ERROR,*9999)
            CALL GKS_SCHH(0.02,ERROR,*9999)
            IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
              CALL GKS_STXP(GRIGHT,ERROR,*9999)
              CALL GKS_SCHUP(0.0,1.0,ERROR,*9999)
            ENDIF
          ELSE IF(IGEOM.EQ.3) THEN
            CALL GKS_STXAL(GACENT,GAHALF,ERROR,*9999)
            CALL GKS_SCHH(0.02,ERROR,*9999)
            IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
              CALL GKS_STXP(GDOWN,ERROR,*9999)
              CALL GKS_SCHUP(0.0,1.0,ERROR,*9999)
            ENDIF
          ELSE IF(IGEOM.EQ.4) THEN
            CALL GKS_STXAL(GARITE,GAHALF,ERROR,*9999)
            CALL GKS_SCHH(0.02,ERROR,*9999)
            IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
              CALL GKS_STXP(GRIGHT,ERROR,*9999)
              CALL GKS_SCHUP(0.0,1.0,ERROR,*9999)
            ENDIF
          ENDIF
        ENDIF
!       plot x against y
        CALL GKS_TX(D_PT(1),D_PT(2),STRING,ERROR,*9999)

      ENDIF

      CALL EXITS('TEXT')
      RETURN
 9999 CALL ERRORS('TEXT',ERROR)
      CALL EXITS('TEXT')
      RETURN 1
      END


      SUBROUTINE TRANSFORM_SEGMENT(ISEGM_LIST,iw,TYPE,XCENTRE,YCENTRE,
     '  VALUE,ERROR,*)

C#### Subroutine: TRANSFORM_SEGMENT
C###  Description:
C**** Transforms graphics segments in ISEGM_LIST(isegm),isegm=1,
C**** ISEGM_LIST(0) on iw with VALUE for:
C**** TYPE is 'rotate','scale','x-translate' or 'y-translate'.
C**** Rotate and scale operate about XCENTRE,YCENTRE.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:draw00.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.FOR'
!     Parameter List
      INTEGER ISEGM_LIST(0:*),iw
      REAL XCENTRE,YCENTRE,VALUE
      CHARACTER ERROR*(*),TYPE*(*)
!     Local Variables
      INTEGER isegm,nosg
      CHARACTER LOCAL_TYPE*11

      CALL ENTERS('TRANSFORM_SEGMENT',*9999)
      LOCAL_TYPE=TYPE
      IF(LOCAL_TYPE(1:5).EQ.'scale') THEN
        IF(VALUE.GT.0.0) THEN
          VALUE= 1.0+VALUE
        ELSE
          VALUE=1.0/(1.0-VALUE)
        ENDIF
      ENDIF
      IF(DOP) THEN
        WRITE(OP_STRING,'('' Type='',A,'' value='',E12.3)')
     '    LOCAL_TYPE,VALUE
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      DO isegm=1,ISEGM_LIST(0)
        nosg=ISEGM_LIST(isegm)

        IF(LOCAL_TYPE(1:6).EQ.'rotate') THEN
          IF(.NOT.LTRANS(nosg)) THEN
            LTRANS(nosg)=.TRUE.
            CALL GKS_EVTM(XCENTRE,YCENTRE,0.0,0.0,VALUE,1.0,1.0,GWC,
     '        TRANSF(1,1,nosg),ERROR,*9999)
          ELSE IF(LTRANS(nosg)) THEN
            CALL GKS_ACTM(TRANSF(1,1,nosg),XCENTRE,YCENTRE,0.0,0.0,
     '        VALUE,1.0,1.0,GWC,TRANSF(1,1,nosg),ERROR,*9999)
          ENDIF

        ELSE IF(LOCAL_TYPE(1:5).EQ.'scale') THEN
          IF(.NOT.LTRANS(nosg)) THEN
            LTRANS(nosg)=.TRUE.
            CALL GKS_EVTM(XCENTRE,YCENTRE,0.0,0.0,0.0,VALUE,VALUE,GWC,
     '        TRANSF(1,1,nosg),ERROR,*9999)
          ELSE IF(LTRANS(nosg)) THEN
            CALL GKS_ACTM(TRANSF(1,1,nosg),XCENTRE,YCENTRE,0.0,0.0,0.0,
     '        VALUE,VALUE,GWC,TRANSF(1,1,nosg),ERROR,*9999)
          ENDIF

        ELSE IF(LOCAL_TYPE(1:11).EQ.'x-translate') THEN
          IF(.NOT.LTRANS(nosg)) THEN
            LTRANS(nosg)=.TRUE.
            CALL GKS_EVTM(0.0,0.0,VALUE,0.0,0.0,1.0,1.0,GWC,
     '        TRANSF(1,1,nosg),ERROR,*9999)
          ELSE IF(LTRANS(nosg)) THEN
            CALL GKS_ACTM(TRANSF(1,1,nosg),0.0,0.0,VALUE,0.0,0.0,1.0,
     '        1.0,GWC,TRANSF(1,1,nosg),ERROR,*9999)
          ENDIF
          XSEGMENT_DATA(1,nosg)=XSEGMENT_DATA(1,nosg)+VALUE

        ELSE IF(LOCAL_TYPE(1:11).EQ.'y-translate') THEN
          IF(.NOT.LTRANS(nosg)) THEN
            LTRANS(nosg)=.TRUE.
            CALL GKS_EVTM(0.,0.,0.,VALUE,0.,1.,1.,GWC,TRANSF(1,1,nosg),
     '        ERROR,*9999)
          ELSE IF(LTRANS(nosg)) THEN
            CALL GKS_ACTM(TRANSF(1,1,nosg),0.,0.,0.,VALUE,0.,1.,1.,GWC,
     '        TRANSF(1,1,nosg),ERROR,*9999)
          ENDIF
          YSEGMENT_DATA(1,nosg)=YSEGMENT_DATA(1,nosg)+VALUE
        ENDIF

        CALL GKS_SSGT(nosg,TRANSF(1,1,nosg),ERROR,*9999)
        CALL GKS_UWK(iw,GPERFO,ERROR,*9999)
      ENDDO

      CALL EXITS('TRANSFORM_SEGMENT')
      RETURN
 9999 CALL ERRORS('TRANSFORM_SEGMENT',ERROR)
      CALL EXITS('TRANSFORM_SEGMENT')
      RETURN 1
      END


      SUBROUTINE VALUATOR(LABEL,iw,MODE,NTYPE,VALMIN,VALMAX,VALINI,
     '  VALUE,XREF,YREF,ERROR,*)

C#### Subroutine: VALUATOR
C###  Description:
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

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:disp00.cmn'
      INCLUDE 'cmiss$reference:echo00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:gks000.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
!     Parameter List
      INTEGER iw,NTYPE
      REAL VALINI,VALMAX,VALMIN,VALUE,XREF,YREF
      CHARACTER ERROR*(*),LABEL*(*),MODE*(*)
!     Local Variables
      INTEGER i,IMODE,INSTAT,LD1
      REAL XFACTOR,YFACTOR
      CHARACTER DUMMY

      DATA LD1/1/,XFACTOR/0.15/,YFACTOR/0.1/,DUMMY/' '/

      CALL ENTERS('VALUATOR',*9999)

      IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
        XFACTOR=0.15
        YFACTOR=0.1
      ELSE IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
        XFACTOR=0.15
        YFACTOR=0.125
      ENDIF

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
        WRITE(OP_STRING,'('' XREF='',E13.5,'' YREF='',E13.5)')
     '    XREF,YREF
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' DISP='',E13.5,'' ECAREA(1..4):'','
     '    //'4E13.5)')
     '    DISP,(ECAREA(i),i=1,4)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      CALL SETUP(iw,ERROR,*9999)
      CALL GKS_SVLM(iw,LD1,GREQU,GECHO,ERROR,*9999)
      CALL GKS_INVL(iw,LD1,VALINI,1,VALMIN,VALMAX,0,DUMMY,ERROR,*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,*) ' Call to GINVL completed'
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      IF(MODE(1:5).EQ.'EVENT') THEN
        IMODE=GEVENT
      ELSE IF(MODE(1:6).EQ.'SAMPLE') THEN
        IMODE=GSAMPL
      ELSE IF(MODE(1:7).EQ.'REQUEST') THEN
        IMODE=GREQU
      ENDIF
      CALL GKS_SVLM(iw,LD1,IMODE,GECHO,ERROR,*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,*) 'Mode=',MODE
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      IF(IMODE.EQ.GREQU) THEN
        CALL GKS_RQVL(iw,LD1,INSTAT,VALUE,ERROR,*9999)
        IF(instat.ne.GOK) THEN
          VALUE=VALINI
          INSTAT=0
        ELSE
          INSTAT=1
        ENDIF
      ENDIF

      CALL EXITS('VALUATOR')
      RETURN

 9999 CALL ERRORS('VALUATOR',ERROR)
      CALL EXITS('VALUATOR')
      RETURN 1
      END


      SUBROUTINE VISIB(iw,ISEG,ISEGNUM,CLASS,ERROR,*)

C#### Subroutine: VISIB
C###  Description:
C**** change ISEGNUM to CLASS='VISIBLE' or 'INVISIBLE'
C**** change ISEG(ISEGNUM) to 2 if visible, 1 if not.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
!     Parameter List
      INTEGER ISEG(*),ISEGNUM,iw
      CHARACTER CLASS*(*),ERROR*(*)

      CALL ENTERS('VISIB',*9999)
      IF(IWKT(iw).EQ.1) THEN      !GKS window
        IF(CLASS(1:7).EQ.'VISIBLE') THEN
          CALL GKS_SVIS(ISEGNUM,GVISI,ERROR,*9999)
          ISEG(ISEGNUM)=2
        ELSE
          CALL GKS_SVIS(ISEGNUM,GINVIS,ERROR,*9999)
          ISEG(ISEGNUM)=1
        ENDIF
      ELSE IF(IWKT(iw).EQ.2) THEN !Phigs window
        IF(CLASS(1:7).EQ.'VISIBLE') THEN
          ISEG(ISEGNUM)=2
        ELSE
          ISEG(ISEGNUM)=1
        ENDIF
        CALL INVIS(iw,ISEG,ERROR,*9999)
      ENDIF

      CALL EXITS('VISIB')
      RETURN
 9999 CALL ERRORS('VISIB',ERROR)
      CALL EXITS('VISIB')
      RETURN 1
      END


      SUBROUTINE WKST_WINDOW(iw,XNDC1,XNDC2,YNDC1,YNDC2,ERROR,*)

C#### Subroutine: WKST_WINDOW
C###  Description:
C**** Change workstation window.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbwk01.cmn'
!     Parameter List
      INTEGER iw
      REAL XNDC1,XNDC2,YNDC1,YNDC2
      CHARACTER ERROR*(*)

      CALL ENTERS('WKST_WINDOW',*9999)
      IF(IWKT(iw).EQ.1) THEN      !GKS window
        CALL GKS_SWKWN(iw,XNDC1,XNDC2,YNDC1,YNDC2,ERROR,*9999)
      ELSE IF(IWKT(iw).EQ.2) THEN !Phigs window
      ENDIF

      CALL EXITS('WKST_WINDOW')
      RETURN
 9999 CALL ERRORS('WKST_WINDOW',ERROR)
      CALL EXITS('WKST_WINDOW')
      RETURN 1
      END
