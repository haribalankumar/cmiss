      SUBROUTINE DEREGI(STRING,ERROR,*)

C#### Subroutine: DEREGI
C###  Description:
C###    DEREGI defines number of regions (NRT).

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,IPFILE
      CHARACTER FILE*(MXCH),STATUS*3
      LOGICAL ALL_REGIONS,CALCU,FILIO,FIRST_TIME,GENER,MOUSE

      CALL ENTERS('DEREGI',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define regions;l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Defines the number of regions. Parameters can be
C###    read from or written to the file FILENAME.ipregi in the directory
C###    specified by PATH with $current specifing the current
C###    default file.

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEREGI',ERROR,*9999)
      ELSE
        IPFILE=1 !is input file version number on 22-May-93
        CALL PARSE_QUALIFIERS('DLPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

        IF(FILIO) THEN
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          ALL_REGIONS=.FALSE. !when parse_regions not called
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'regi',STATUS,
     '        ERR,ERROR,*9999)
            CALL IPREGI(ERROR,*9999)
            CALL CLOSEF(IFILE,ERROR,*9999)
            IF(BACKUP) IOTYPE=2
          ENDDO
        ENDIF
        CALL ASSERT(NRT.LE.NRM,'>>NRM too small',ERROR,*9999)
      ENDIF

      CALL EXITS('DEREGI')
      RETURN
 9999 CALL ERRORS('DEREGI',ERROR)
      CALL EXITS('DEREGI')
      RETURN 1
      END


