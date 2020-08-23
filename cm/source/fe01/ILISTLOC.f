      INTEGER*4 FUNCTION ILISTLOC(ILIST,I)

C#### Function: ILISTLOC
C###  Type: INTEGER
C###  Description:
C###    ILISTLOC returns a pointer to an integer in a list.
C###    The routine is useful in FORTRAN77 when only a pointer to
C###    the list is available

!     Parameter List
      INTEGER ILIST(*),I

C***  Fortran 90 has LOC instead of %LOC
      ILISTLOC=%LOC(ILIST(I))

      RETURN
      END


