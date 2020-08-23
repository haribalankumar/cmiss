      INTEGER FUNCTION LEN_TRIM1(STRING)

C#### Function: LEN_TRIM1
C###  Type: INTEGER
C###  Description:
C###    A FORTRAN 77 string must be at least one character long.  This
C###    function therefore calculates the minimum length of a string
C###    without the trailing blanks.  i.e. Returns the position of the
C###    last non-blank character in STRING, or 1 if there are no
C###    non-blank characters.
C###    KAT 20/3/00.

!     Parameter List
      CHARACTER*(*) STRING
!     Functions
      INTEGER LEN_TRIM

      LEN_TRIM1=MAX(1,LEN_TRIM(STRING))

      RETURN
      END


