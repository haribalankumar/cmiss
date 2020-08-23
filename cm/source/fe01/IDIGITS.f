      INTEGER FUNCTION IDIGITS(I)

C#### Function: IDIGITS
C###  Type: INTEGER
C###  Description:
C###    Counts the number of digits in the integer I.

!     Parameter List
      INTEGER I
!     Local Variables
      INTEGER ITEMP

      IF(I.GT.0) THEN
        ITEMP=I
        IDIGITS=0
      ELSE !zero needs one digit, negative numbers need a minus sign
        ITEMP=-I
        IDIGITS=1
      ENDIF
      DO WHILE(ITEMP.GT.0)
        ITEMP=ITEMP/10
        IDIGITS=IDIGITS+1
      ENDDO

      RETURN
      END


