      SUBROUTINE DEITER(STRING,ERROR,*)

C#### Subroutine: DEITER
C###  Description:
C###    DEITER defines iterations

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,IPFILE
      CHARACTER FILE*(MXCH),STATUS*3
      LOGICAL ALL_REGIONS,CALCU,FILIO,GENER,MOUSE

      CALL ENTERS('DEITER',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: define iteration;d/l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Defines iteration parameters.

        OP_STRING(1)=STRING(1:IEND)//';d/l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEITER',ERROR,*9999)
      ELSE
        IPFILE=1 !is input file version number - CPB 8/2/94
        CALL PARSE_QUALIFIERS('DLPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

        IF(FILIO) THEN
          ALL_REGIONS=.FALSE. !when parse_regions not called
          CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'iter',STATUS,
     '      ERR,ERROR,*9999)
          CALL IPITER(ERROR,*9999)
          CALL CLOSEF(IFILE,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('DEITER')
      RETURN
 9999 CALL ERRORS('DEITER',ERROR)
      CALL EXITS('DEITER')
      RETURN 1
      END


