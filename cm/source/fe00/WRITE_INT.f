      SUBROUTINE WRITE_INT(IUNIT,INT,ERR)

C#### Subroutine: WRITE_INT
C###  Description:
C###    Writes INT to output identified by IUNIT using the minimum
C###    number of characters.
C###  See-Also: WRITE_STRING

      IMPLICIT NONE
!     Parameter List
      INTEGER IUNIT,INT
      INTEGER ERR
!     Local Variables
      INTEGER LENGTH,BUFSIZE
      PARAMETER(BUFSIZE=99) ! will an integer require >99 chars
      CHARACTER FORMAT*5,BUFFER*(BUFSIZE)
!     Functions
      INTEGER IDIGITS

      LENGTH=IDIGITS(INT)
      IF(LENGTH.GT.BUFSIZE) THEN
        CALL FLAG_ERROR(-1,'BUFSIZE too small in WRITE_INT')
        CALL WRITE_STRING(IUNIT,1,'*',ERR)
        ERR=1
        RETURN
      ENDIF
      WRITE(FORMAT,'(''(I'',I2.2,'')'')') LENGTH
      WRITE(BUFFER(:LENGTH),FORMAT) INT
      CALL WRITE_STRING(IUNIT,LENGTH,BUFFER,ERR)
      RETURN
      END


