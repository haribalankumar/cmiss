      INTEGER FUNCTION LEN_TRIM(STRING)

C#### Function: LEN_TRIM
C###  Type: INTEGER
C###  Description:
C###    Simulates Fortran 90 function of the same name.
C###    Returns the position of the last non-blank character in STRING.
C###    This is zero if there are no non-blank characters.
C###    KAT 20/3/00.

!     Parameter List
      CHARACTER*(*) STRING
!     Local Variables
      INTEGER I

      I=LEN(STRING)
      IF(I.GT.0) THEN
        DO WHILE(STRING(I:I).EQ.' ')
          I=I-1
          IF(I.EQ.0) GO TO 1
        ENDDO
      ENDIF
 1    LEN_TRIM=I

      RETURN
      END


