      INTEGER FUNCTION FSTRINGLEN(STRING)

C#### Function: FSTRINGLEN
C###  Description:
C###    Finds the length of a left justified, blank padded fortran
C###    string.

      IMPLICIT NONE
!     Parameter List
      CHARACTER STRING*(*)
!     Local Variables
      INTEGER I

      DO I=LEN(STRING),1,-1
        IF(STRING(I:I).NE.' ') THEN
          FSTRINGLEN=I
          RETURN
        ENDIF
      ENDDO !I

      FSTRINGLEN=0

      RETURN
      END


