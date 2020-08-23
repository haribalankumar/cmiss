      SUBROUTINE APPENDI(IEND,IDATA,STRING)

C#### Subroutine: APPENDI
C###  Description:
C###    APPENDI writes an integer after a specified character in a string
C###    and calculates the end of the integer in the string.

      IMPLICIT NONE
!     Parameter List
      INTEGER IEND,IDATA
      CHARACTER STRING*(*)
!     Local Variables
      INTEGER IBEG,LENGTH1,LENGTH2
      CHARACTER FORMAT*5
!     Functions
      INTEGER IDIGITS

      DATA FORMAT/'(I**)'/

C      CALL ENTERS('APPENDI',*9999)

      LENGTH2=LEN(STRING)
      IF(IEND.LT.LENGTH2) THEN
        LENGTH1=IDIGITS(IDATA)
        WRITE(FORMAT(3:4),'(I2.2)') LENGTH1
        IBEG=IEND+1
        IEND=MIN(IEND+LENGTH1,LENGTH2)
        WRITE(STRING(IBEG:IEND),FORMAT) IDATA
      ENDIF

C      CALL EXITS('APPENDI')
      RETURN
      END


