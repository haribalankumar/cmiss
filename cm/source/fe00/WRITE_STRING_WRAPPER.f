      SUBROUTINE WRITE_STRING_WRAPPER(IUNIT,STRLEN,STRING,ERR)

C#### Subroutine: WRITE_STRING_WRAPPER
C###  Description:
C###    Merely calls WRITE_STRING with the same arguments.  g77 does not
C###    allow a routine to call itself.  This works around that.
C###  See-Also: WRITE_STRING

      IMPLICIT NONE
!     Parameter List
      INTEGER IUNIT,STRLEN
      CHARACTER STRING(STRLEN)
      INTEGER ERR

      CALL WRITE_STRING(IUNIT,STRLEN,STRING,ERR)
        
      END


