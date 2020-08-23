      INTEGER FUNCTION INDEX_TEXT(icolour,WIDTH_TEXT,FONT_TEXT,RGB_TEXT)

C#### Function: INDEX_TEXT
C###  Type: INTEGER
C###  Description:
C###    <HTML> <PRE>
C###    INDEX_TEXT returns text bundle index for given text size, font
C###    and RGB value. For given text size and RGB value, as follows:
C###    Text Index = 1       black font1    width1
C###                 2         "   font2      "
C###                 3         "   font3      "
C###                 4         "   font4      "
C###                 5-8         as above for width2
C###
C###                 9 -10   red    font1,2 width1
C###                 11-12    "     font1,2 width2
C###                 13-16   green  as above
C###                 17-20   blue       "
C###                 21-24   cyan       "
C###                 25-28   yellow     "
C###                 29-32   white      "
C###                 33-36   light blue "
C###                 37-40   grey       "
C###    </PRE> </HTML>

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'colo00.cmn'
!     Parameter List
      INTEGER icolour
      CHARACTER FONT_TEXT*(*),RGB_TEXT*(*),WIDTH_TEXT*(*)
!     Local Variables
      INTEGER index,INDEX_OFFSET_1,INDEX_OFFSET_2
      CHARACTER ERROR*10,FONT*8,RGB*8,WIDTH*8

      IF(icolour.EQ.0) THEN
C        FONT=CUPPER(FONT_TEXT)
C        WIDTH=CUPPER(WIDTH_TEXT)
C        RGB =CUPPER(RGB_TEXT)
        CALL CUPPER(FONT_TEXT,FONT)
        CALL CUPPER(WIDTH_TEXT,WIDTH)
        CALL CUPPER(RGB_TEXT,RGB)

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


