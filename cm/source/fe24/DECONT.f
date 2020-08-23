      SUBROUTINE DECONT(STRING,ERROR,*)

C#### Subroutine: DECONT
C###  Description:
C###    DECONT defines parameters for contact mechanics problems
C###    with prompted input or from filename.ipcont.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'

!     Parameter List
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,IPFILE
      CHARACTER FILE*(MXCH),STATUS*3
      LOGICAL ALL_REGIONS,CALCU,FILIO,FIRST_TIME,GENER,MOUSE

      CALL ENTERS('DECONT',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define contact;l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Define contact mechanics parameters. Parameter values are
C###    read from or written to the file FILENAME.ipcont in the directory
C###    specified by PATH with $current specifing the current default file.

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     &    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     &    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DECONT',ERROR,*9999)
      ELSE

        CALL PARSE_QUALIFIERS('DLPRW',noco,1,CO,COQU,
     &    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)

        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

        CALL ASSERT(CALL_EQUA,'>>Equation not defined',ERROR,*9999)

        IF(FILIO) THEN
          IPFILE=1 ! file version number
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          ALL_REGIONS=.FALSE.
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'cont',STATUS,
     &        ERR,ERROR,*9999)
            CALL IPCONT(ERROR,*9999)
            CALL CLOSEF(IFILE,ERROR,*9999)
            IF(BACKUP) IOTYPE=2
          ENDDO
        ENDIF !filio
        CALL_CONT=.TRUE.
      ENDIF

      CALL EXITS('DECONT')
      RETURN
 9999 CALL ERRORS('DECONT',ERROR)
      CALL EXITS('DECONT')
      RETURN 1
      END


