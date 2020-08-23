      INTEGER FUNCTION INDEX_POLYLINE(icolour,TYPE_LINE,WIDTH_LINE,
     ' RGB_LINE)

C#### Function: INDEX_POLYLINE
C###  Type: INTEGER
C###  Description:
C###    <HTML> <PRE>
C###    INDEX_POLYLINE returns polyline bundle index:
C###    If icolour=0, for given line type, width and RGB value,
C###    as follows:
C###    Polyline Index = 1       black solid    width1
C###                     2         "   dotted      "
C###                     3         "   dashed      "
C###                     4         "   dot-dash    "
C###                     5-8         as above   width2
C###
C###                     9 -12   red    as above width1
C###                     13-16   green      "
C###                     17-20   blue       "
C###                     21-24   cyan       "
C###                     25-28   yellow     "
C###                     29-32   white      "
C###                     33-36   light blue "
C###                     37-40   grey       "
C###
C###    Else if 1<=icolour<=216 for given icolour as follows:
C###    Polyline Index = 41-256  216 colours solid width2.
C###    </PRE> </HTML>

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'colo00.cmn'
!     Parameter List
      INTEGER icolour
      CHARACTER RGB_LINE*(*),TYPE_LINE*(*),WIDTH_LINE*(*)
!     Local Variables
      INTEGER index,INDEX_OFFSET_1,INDEX_OFFSET_2
      CHARACTER ERROR*10,TYPE*8,WIDTH*8,RGB*255

      IF(icolour.EQ.0) THEN
C        TYPE =CUPPER(TYPE_LINE)
C        WIDTH=CUPPER(WIDTH_LINE)
C        RGB  =CUPPER(RGB_LINE)
        CALL CUPPER(TYPE_LINE,TYPE)
        CALL CUPPER(WIDTH_LINE,WIDTH)
        CALL CUPPER(RGB_LINE,RGB)

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


