      INTEGER FUNCTION INDEX_POLYMARKER(icolour,TYPE_MARKER,SIZE_MARKER,
     '  RGB_MARKER)

C#### Function: INDEX_POLYMARKER
C###  Type: INTEGER
C###  Description:
C###    <HTML>
C###    INDEX_POLYMARKER returns polymarker bundle index for given
C###    marker type, size and RGB value.
C###    <PRE>
C###    If icolour=0  for given marker type, size and RGB value,
C###    as follows:
C###    Polymarker Index = 1       black plus      size1
C###                       2         "   asterisk    "
C###                       3         "   circle      "
C###                       4         "   point       "
C###                       5-8         as above    size2
C###
C###                       9 -12   red    as above  size1
C###                       13-16   green      "
C###                       17-20   blue       "
C###                       21-24   cyan       "
C###                       25-28   yellow     "
C###                       29-32   white      "
C###                       33-36   light blue "
C###                       37-40   grey       "
C###
C###    Else if 1<=icolour<=216 for given icolour as follows:
C###    Polymarker Index = 41-256  216 colours plus size2.
C###    </PRE> </HTML>

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'colo00.cmn'
!     Parameter List
      INTEGER icolour
      CHARACTER RGB_MARKER*(*),SIZE_MARKER*(*),TYPE_MARKER*(*)
!     Local Variables
      INTEGER IBEG,IEND,index,INDEX_OFFSET_1,INDEX_OFFSET_2
      CHARACTER ERROR*10,RGB*8,SIZE*8,TYPE*8

      IF(icolour.EQ.0) THEN
C        TYPE=CUPPER(TYPE_MARKER)
C        SIZE=CUPPER(SIZE_MARKER)
C        RGB =CUPPER(RGB_MARKER)
CC AJPs 28/4/98
        CALL STRING_TRIM(TYPE_MARKER,IBEG,IEND)
        CALL CUPPER(TYPE_MARKER(IBEG:IEND),TYPE)
        CALL STRING_TRIM(SIZE_MARKER,IBEG,IEND)
        CALL CUPPER(SIZE_MARKER(IBEG:IEND),SIZE)
        CALL STRING_TRIM(RGB_MARKER,IBEG,IEND)
        CALL CUPPER(RGB_MARKER(IBEG:IEND),RGB)
c        CALL CUPPER(TYPE_MARKER,TYPE)
c        CALL CUPPER(SIZE_MARKER,SIZE)
c        CALL CUPPER(RGB_MARKER,RGB)
CC AJPe 28/4/98

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
     '    //''' RGB='',A,'' INDEX='',I3)') TYPE,RGB,index
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


