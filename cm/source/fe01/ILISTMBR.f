      INTEGER FUNCTION ILISTMBR(ILIST,I)

C#### Function: ILISTMBR
C###  Type: INTEGER
C###  Description:
C###    ILISTMBR returns an integer in a list.
C###    The routine is useful in FORTRAN77 when only a pointer to
C###    the list is available

!     Parameter List
      INTEGER ILIST(*),I

      ILISTMBR=ILIST(I)

      RETURN
      END


