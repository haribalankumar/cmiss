      LOGICAL FUNCTION ABBREV(SHORT,LONG,MNCH)

C#### Function: ABBREV
C###  Type: LOGICAL
C###  Description:
C###    ABBREV returns .TRUE. if the character string SHORT is an
C###    abbreviation of LONG.  SHORT must be at least MNCH characters
C###    long.  As a side effect, if ABBREV is .TRUE., SHORT is
C###    overwritten with LONG.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
!     Parameter List
      INTEGER MNCH
      CHARACTER LONG*(*),SHORT*(*)
!     Local Variables
      CHARACTER BUFFER*(MXCH)
      LOGICAL ABBRV

C      BUFFER=CUPPER(SHORT)
      CALL CUPPER(SHORT,BUFFER)
      ABBREV=ABBRV(BUFFER,LONG,MNCH)
      IF(ABBREV) THEN
        SHORT=LONG
      ENDIF
      RETURN
C GMH 12/12/96 This does not make sense
C 9999 WRITE(OP_STRING,'('' >>ERROR: RETURN FROM CUPPER'',A)') ERROR
C      CALL WRITES(IOER,OP_STRING,ERROR,*9998)
C 9998 RETURN
      END


