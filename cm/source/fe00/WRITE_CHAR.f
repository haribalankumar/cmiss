      SUBROUTINE WRITE_CHAR(IUNIT,LINE,ERR)

C#### Subroutine: WRITE_CHAR
C###  Description:
C###    Simple interface to WRITE_STRING that allows writing of a fortran
C###    character variable without supplying a length.
C###  See-Also: WRITE_STRING, WRITE_LINE

      IMPLICIT NONE
!     Parameter List
      INTEGER IUNIT
      CHARACTER LINE*(*)
      INTEGER ERR

C     Under the Fortran 90 Standard, we can't pass a scalar to an array,
C     but this seems to work.  If necessary we could have another
C     intermediate routine to copy the character to an automatic
C     character array.
      CALL WRITE_STRING(IUNIT,LEN(LINE),LINE,ERR)
      RETURN
      END


