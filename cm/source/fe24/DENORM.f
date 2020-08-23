      SUBROUTINE DENORM(NELIST,NW,STRING,ERROR,*)

C#### Subroutine: DENORM
C###  Description:
C###    DENORM defines which elements have their normals reversed
C###    for boundary element problems. A value of NW(ne,3)=1
C###    means the normal for element ne will be reversed. Any
C###    other value will not change the normal direction. The
C###    default value is zero.
C***  Created by Martin Buist April 1998

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'

!     Parameter List
      INTEGER NELIST(0:NEM),NW(NEM,3,NXM)
      CHARACTER ERROR*(*),STRING*(*)

!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,IPFILE
      LOGICAL CALCU,FILIO,FIRST_TIME,GENER,MOUSE
      CHARACTER FILE*(MXCH),STATUS*3

      CALL ENTERS('DENORM',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C--------------------------------------------------------------------

C#### Command: FEM define normal;l;p;r;w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Defines which elements have their normals reversed
C###    for boundary element problems.

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C--------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DENORM',ERROR,*9999)
      ELSE
        IPFILE=1
        CALL ASSERT(CALL_ELEM,'>>Must define elements first',
     '    ERROR,*9999)
        CALL ASSERT(CALL_EQUA,'>>Must define equation first',
     '    ERROR,*9999)
        CALL ASSERT(USE_BEM.EQ.1,'>>Set USE_BEM to 1 in ippara',
     '    ERROR,*9999)

C        nx=1 ! MPN 14Jun2000  may need generalising

        CALL PARSE_QUALIFIERS('LPRW',noco,1,CO,COQU,CALCU,FILIO,
     '    GENER,MOUSE,STATUS,ERROR,*1)
        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

        IF(FILIO) THEN
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            CALL OPEN_FILE(.FALSE.,IPFILE,FILE,'norm',STATUS,
     '        ERR,ERROR,*9999)
            CALL IPNORM(NELIST,NW,ERROR,*9999)
            CALL CLOSEF(IFILE,ERROR,*9999)
            IF(BACKUP) IOTYPE=2
          ENDDO
        ENDIF
      ENDIF

      CALL EXITS('DENORM')
      RETURN
 9999 CALL ERRORS('DENORM',ERROR)
      CALL EXITS('DENORM')
      RETURN 1
      END


