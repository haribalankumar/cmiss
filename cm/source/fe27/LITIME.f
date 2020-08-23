      SUBROUTINE LITIME(NTIME_INTERP,NTIME_POINTS,TIME_VALUES,STRING,
     '  TIME_VARIABLE_NAMES,ERROR,*)

C#### Subroutine: LITIME
C###  Description:
C###    LITIME lists time variable information
C**** Created by Martin Buist, 15 November 1999.

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'

!     Parameter List
      INTEGER NTIME_INTERP(NTIMEVARSM),NTIME_POINTS(NTIMEVARSM)
      REAL*8 TIME_VALUES(2,0:NTIMEPOINTSM+1,NTIMEVARSM)
      CHARACTER ERROR*(*),TIME_VARIABLE_NAMES(NTIMEVARSM)*(*),
     '  STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND
      CHARACTER FILE*100
      LOGICAL OPFILE

      CALL ENTERS('LITIME',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list time_variable<;FILENAME>
C###  Description:
C###    Lists time variable(s) to the screen or to FILENAME.optime

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LITIME',ERROR,*9999)
      ELSE
        CALL ASSERT(CALL_TIME,'>>Must define time first',ERROR,*9999)
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.optime','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
          IOFI=IOOP
        ENDIF

        CALL OPTIME(NTIME_INTERP,NTIME_POINTS,TIME_VALUES,
     '    TIME_VARIABLE_NAMES,ERROR,*9999)

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LITIME')
      RETURN
 9999 CALL ERRORS('LITIME',ERROR)
      CALL EXITS('LITIME')
      RETURN 1
      END


C LKC 11-NOV-1998 *** archived ***
C SUBROUTINE LITRAN(STRING,ERROR,*)


