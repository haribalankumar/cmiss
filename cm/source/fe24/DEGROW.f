      SUBROUTINE DEGROW(STRING,ERROR,*)

C#### Subroutine: DEGROW
C###  Description:
C###    DEGROW defines growth law.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2
      CHARACTER FILE*(MXCH),STATUS*3
      LOGICAL CALCU,FILIO,GENER,MOUSE

      CALL ENTERS('DEGROW',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define growth;d/l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Defines growth law parameters. Parameters can be read from
C###    or written to the file FILENAME.ipgrow, with $current specifing
C###    the current default file.

        OP_STRING(1)=STRING(1:IEND)//';d/l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEGROW',ERROR,*9999)
      ELSE
        CALL PARSE_QUALIFIERS('DLPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

        IF(FILIO) THEN
          CALL STRING_TRIM(FILE,IBEG,IEND)
          CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.ipgrow',STATUS,
     '      'DIRECT','FORMATTED',132,ERROR,*9999)
          CALL IPGROW(ERROR,*9999)
          CALL CLOSEF(IFILE,ERROR,*9999)
        ENDIF
        CALL_GROW=.TRUE.
      ENDIF

      CALL EXITS('DEGROW')
      RETURN
 9999 CALL ERRORS('DEGROW',ERROR)
      CALL EXITS('DEGROW')
      RETURN 1
      END


