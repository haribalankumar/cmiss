      SUBROUTINE DELEAD(STRING,ERROR,*)

C#### Subroutine: DELEAD
C###  Description:
C###    DELEAD defines electrocardiographic leads.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,IPFILE
      CHARACTER FILE*(MXCH),STATUS*3
      LOGICAL ALL_REGIONS,CALCU,FILEIO,GENER,MOUSE

      CALL ENTERS('DELEAD',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define leads;l/m/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Defines electrocardiographic lead parameters. Parameters can be
C###    read from or written to the file FILENAME.iplead, with $current
C###    specifing the current default file.

        OP_STRING(1)=STRING(1:IEND)//';l/m/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DELEAD',ERROR,*9999)
      ELSE

        IPFILE=1 !is input file version number on 28/5/95

        CALL PARSE_QUALIFIERS('LPRW',noco,1,CO,COQU,
     '    CALCU,FILEIO,GENER,MOUSE,STATUS,ERROR,*1)

        IF(FILEIO) THEN
          ALL_REGIONS=.FALSE. !when parse_regions not called
          CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'lead',STATUS,
     '      ERR,ERROR,*9999)
          CALL IPLEAD(ERROR,*9999)
          CALL CLOSEF(IFILE,ERROR,*9999)

C LKC 22-JUL-98 add logical
          CALL_LEAD=.TRUE.
        ENDIF

      ENDIF

      CALL EXITS('DELEAD')
      RETURN
 9999 CALL ERRORS('DELEAD',ERROR)
      CALL EXITS('DELEAD')
      RETURN 1
      END


