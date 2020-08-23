      SUBROUTINE PARSE_QUALIFIERS(QLIST,noco,nocoqu,CO,COQU,
     '  CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*)

C#### Subroutine: PARSE_QUALIFIERS
C###  Description:
C###    PARSE_QUALIFIERS parses command qualifier COQU(noco,nocoqu)
C###    to see if member of QLIST and sets:
C###    logical variables CALCU,FILIO,GENER,MOUSE
C###    character  "      STATUS
C###    integer    "      IOTYPE (held in include file b01.cmn)

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
!     Parameter List
      INTEGER noco,nocoqu
      LOGICAL CALCU,FILIO,GENER,MOUSE
      CHARACTER CO(*)*(*),COQU(25,*)*(*),ERROR*(*),QLIST*(*),
     '  STATUS*3
!     Local Variables
      LOGICAL ABBREV,CELEM

      CALL ENTERS('PARSE_QUALIFIERS',*9999)
      IF(CELEM(COQU(noco,nocoqu)(1:1),QLIST).
     '  AND.COQU(noco,nocoqu)(2:2).EQ.' ') THEN
        CALCU=.FALSE.
        FILIO=.FALSE.
        GENER=.FALSE.
        MOUSE=.FALSE.
        IF(ABBREV(COQU(noco,1),'C',1)) THEN
          CALCU=.TRUE.
        ELSE IF(ABBREV(COQU(noco,1),'D',1)) THEN
          FILIO=.TRUE.
          IOTYPE=5
        ELSE IF(ABBREV(COQU(noco,1),'G',1)) THEN
          FILIO=.TRUE.
          GENER=.TRUE.
          IOTYPE=1
          STATUS='NEW'
        ELSE IF(ABBREV(COQU(noco,1),'L',1)) THEN
          FILIO=.TRUE.
          IOTYPE=4
          STATUS='OLD'
        ELSE IF(ABBREV(COQU(noco,1),'M',1)) THEN
          MOUSE=.TRUE.
          IOTYPE=3
          STATUS='NEW'
        ELSE IF(ABBREV(COQU(noco,1),'P',1)) THEN
          FILIO=.TRUE.
          IOTYPE=1
          STATUS='NEW'
        ELSE IF(ABBREV(COQU(noco,1),'R',1)) THEN
          FILIO=.TRUE.
          IOTYPE=2
          STATUS='OLD'
        ELSE IF(ABBREV(COQU(noco,1),'W',1)) THEN
          FILIO=.TRUE.
          IOTYPE=3
          STATUS='NEW'
        ENDIF

      ELSE
        CO(noco+1)='?'
        ERROR='>>Reenter: '//CO(noco)
        GO TO 9999
      ENDIF

      CALL EXITS('PARSE_QUALIFIERS')
      RETURN
 9999 CALL EXITS('PARSE_QUALIFIERS')
      RETURN 1
      END


