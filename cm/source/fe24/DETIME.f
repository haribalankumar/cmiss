      SUBROUTINE DETIME(NTIME_INTERP,NTIME_POINTS,TIME_VALUES,
     '  STRING,TIME_VARIABLE_NAMES,ERROR,*)

C#### Subroutine: DETIME
C###  Description:
C###    Reads in time variable(s) from a .iptime file. This allows
C###    the user to set specific time points and the value of the
C###    time variables at those points. This is currently used to
C###    input time dependent boundary conditions for some grid
C###    problems.
C###    Created by Martin Buist, 15 November 1999.

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'

!     Parameter List
      INTEGER NTIME_INTERP(NTIMEVARSM),NTIME_POINTS(NTIMEVARSM)
      REAL*8 TIME_VALUES(2,0:NTIMEPOINTSM+1,NTIMEVARSM)
      CHARACTER ERROR*(*),TIME_VARIABLE_NAMES(NTIMEVARSM)*(*),
     '  STRING*(*)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,IPFILE
      LOGICAL ALL_REGIONS,CALCU,FILIO,GENER,MOUSE
      CHARACTER FILE*(MXCH),STATUS*3

      CALL ENTERS('DETIME',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C--------------------------------------------------------------------

C#### Command: FEM define time_variable;l/r/p/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Reads in time variable(s) from FILENAME.iptime. This allows
C###    the user to set specific time points and the value of the
C###    time variables at those points. This is currently used to
C###    input time dependent boundary conditions for some grid
C###    problems.

        OP_STRING(1)=STRING(1:IEND)//';l/r/p/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C--------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DETIME',ERROR,*9999)
      ELSE

        CALL PARSE_QUALIFIERS('LRPW',noco,1,CO,COQU,CALCU,FILIO,GENER,
     '    MOUSE,STATUS,ERROR,*1)
        IPFILE=3
        ALL_REGIONS=.FALSE.

        IF(FILIO) THEN
          CALL CHECKF(2,NOCO,NTCOQU,CO,COQU,FILE,STRING,*1)

          CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'time',
     '      STATUS,ERR,ERROR,*9999)

          CALL IPTIME(NTIME_INTERP,NTIME_POINTS,TIME_VALUES,
     '      TIME_VARIABLE_NAMES,ERROR,*9999)

          CALL CLOSEF(IFILE,ERROR,*9999)
        ENDIF
        CALL_TIME=.TRUE.
      ENDIF

      CALL EXITS('DETIME')
      RETURN
 9999 CALL ERRORS('DETIME',ERROR)
      CALL EXITS('DETIME')
      RETURN 1
      END


CC AJPs 11-11-97 - 191297 - RGB
