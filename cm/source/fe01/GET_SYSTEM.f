      SUBROUTINE GET_SYSTEM(ERROR,*)

C#### Subroutine: GET_SYSTEM
C###  Description:
C###    GET_SYSTEM returns current system nodename, architecture and
C###    operating system.

      IMPLICIT NONE
      INCLUDE 'disp00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER CSTRSIZE
      PARAMETER (CSTRSIZE=31)
      CHARACTER C_SYSNAME(CSTRSIZE),C_WINDTYPE(CSTRSIZE)
!     Functions
      INTEGER C_STRLEN

      CALL ENTERS('GET_SYSTEM',*9999)

c cpb 16/10/95 Old Way
CC
CC Get the imagename of the CMISS executable being run
CC
C      CALL GETARG(0,EXENAME)
C      CALL STRING_TRIM(EXENAME,BEG,END)
CC Store the name of the executable image (used in fe07)
C      IMAGE_NAME=EXENAME(BEG:END)
CC
CC Find out the file status and information for that executable
CC
C      RETURNCODE=STAT(EXENAME(BEG:END),STAB)
CC
CC Convert the modification time into an ascii string
CC
C      IF(RETURNCODE.EQ.0) THEN
C        MODTIME=CTIME(STAB(10))
C        TYPE *,'CMISS revision time ',MODTIME
C      ENDIF

C cpb 26/10/95 Old way
C!      I=GETHOSTNAME(NODE_NAME,LENGTH)

      OS_TYPE='UNIX'

      CALL GETSYSTEMINFO(%VAL(CSTRSIZE),
     '  %REF(C_SYSNAME(1)),%REF(C_WINDTYPE(1)))
      CALL C2FSTRING(%REF(C_SYSNAME(1)),C_STRLEN(%REF(C_SYSNAME(1))),
     '  SYSNAME)
      CALL C2FSTRING(%REF(C_WINDTYPE(1)),C_STRLEN(%REF(C_WINDTYPE(1))),
     '  WINDOW_TYPE)

      CALL EXITS('GET_SYSTEM')
      RETURN
 9999 CALL ERRORS('GET_SYSTEM',ERROR)
      CALL EXITS('GET_SYSTEM')
      RETURN 1
      END


