      SUBROUTINE PARSRE(STRING,RE,ERROR,*)

C#### Subroutine: PARSRE
C###  Description:
C###    PARSRE converts the character string STRING into a REAL*8
C###    variable RE.

      IMPLICIT NONE
!     Parameter List
      REAL*8 RE
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER IBEG,IEND,IOSTAT

      CALL ENTERS('PARSRE',*9999)

C LKC 31-MAY-2022 - Move lower so that it is not alwaysed called.
C   Probably not really necessary as the strings are already trimmed.
C      CALL STRING_TRIM(STRING,IBEG,IEND)

      READ(UNIT=STRING,FMT=*,IOSTAT=IOSTAT) RE

      IF(IOSTAT.LT.0) THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        ERROR=' no items were found in STRING="'
     '    //STRING(IBEG:IEND)//'"'
        GOTO 9999
      ELSE IF(IOSTAT.GT.0) THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        ERROR=' Invalid real number in STRING="'
     '    //STRING(IBEG:IEND)//'"'
        GOTO 9999
      ENDIF

      CALL EXITS('PARSRE')
      RETURN
 9999 CALL ERRORS('PARSRE',ERROR)
      CALL EXITS('PARSRE')
      RETURN 1
      END


