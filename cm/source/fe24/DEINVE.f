      SUBROUTINE DEINVE(STRING,ERROR,*)

C#### Subroutine: DEINVE
C###  Description:
C###    DEINVE defines regularisation options for the inverse of the
C###    transfer matrix.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,IPFILE
      LOGICAL ALL_REGIONS,CALCU,FILIO,FIRST_TIME,GENER,MOUSE
      CHARACTER FILE*(MXCH),STATUS*3

      CALL ENTERS('DEINVE',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define inverse;l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Defines inverse parameters. Parameters can be read from and
C###    written to the file FILENAME.ipinve, with $current specifing
C###    current default file.

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)


C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEINVE',ERROR,*9999)
      ELSE
        CALL PARSE_QUALIFIERS(' DLMPRW',NOCO,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        IF(FILIO) CALL CHECKF(2,NOCO,NTCOQU,CO,COQU,FILE,STRING,*1)

        IF(FILIO) THEN
          IPFILE=2
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          ALL_REGIONS=.FALSE. !when parse_regions not called
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'inve',STATUS,
     '        ERR,ERROR,*9999)
            CALL IPINVE(ERROR,*9999)
            CALL CLOSEF(IFILE,ERROR,*9999)
            IF(BACKUP) IOTYPE=2
          ENDDO
        ENDIF !filio
      ENDIF

      CALL EXITS('DEINVE')
      RETURN
 9999 CALL ERRORS('DEINVE',ERROR)
      CALL EXITS('DEINVE')
      RETURN 1
      END


