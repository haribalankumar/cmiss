      SUBROUTINE DEVALU(STRING,ERROR,*)

C#### Subroutine: DEVALU
C###  Description:
C###    DEVALU defines parameter values in an input file. Currently set 
C###    up for pulmonary problems only. Specific parameters depend on the 
C###    problem type. This is for non-spatially varying parameters with
C###    variable names that are contained in common block files.
C###  Created by AJS Nov 2010

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lungex00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'valu00.cmn'
!     Parameter List
!       INTEGER 
!       REAL*8 
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,IPFILE,N3CO
      CHARACTER STATUS*3,FILE*(MXCH)
      LOGICAL ALL_REGIONS,BACKUP,CALCU,CBBREV,FILIO,FIRST_TIME,
     '  GENER,MOUSE,VENTILATION

      CALL ENTERS('DEVALU',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define values;l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Defines parameters values in an input file. Parameters can be read
C###    from or written to the file FILENAME.ipvalu, with $current
C###    specifing the current default file.
C###  Parameter:      <ventilation>
C###    Pulmonary ventilation problem type.

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<ventilation>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEVALU',ERROR,*9999)
      ELSE
        IPFILE=1 !is input file version number on 24-Jan-1990
        CALL PARSE_QUALIFIERS('LPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)

        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
        IF(CBBREV(CO,'VENTILATION',3,noco+1,NTCO,N3CO)) THEN
          VENTILATION=.TRUE.
        ELSE
        VENTILATION=.FALSE.
        ENDIF

        IF(FILIO) THEN
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          ALL_REGIONS=.FALSE. !when parse_regions not called
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'valu',
     '        STATUS,ERR,ERROR,*9999)
            CALL IPVALU(VENTILATION,ERROR,*9999)
            CALL CLOSEF(IFILE,ERROR,*9999)
            IF(BACKUP) IOTYPE=2
          ENDDO
        ENDIF !filio

        PARAMETERS_DEFINED=.TRUE.

      ENDIF

      CALL EXITS('DEVALU')
      RETURN
 9999 CALL ERRORS('DEVALU',ERROR)
      CALL EXITS('DEVALU')
      RETURN 1
      END



