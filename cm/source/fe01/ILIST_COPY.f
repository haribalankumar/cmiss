      SUBROUTINE ILIST_COPY(LEN,LIST_FROM,LIST_TO)

C#### Subroutine: ILIST_COPY
C###  Description:
C###    ILIST_COPY copies a list of integers to another array.
C###    The routine is useful in FORTRAN77 when only a pointer to one of
C###    the arrays is available

      IMPLICIT NONE
!     Parameter List
      INTEGER LEN,LIST_FROM(*),LIST_TO(*)
!     Local Variables
      INTEGER i
C     CALL ENTERS('ILIST_COPY',*9999)
      DO i=1,LEN
        LIST_TO(i)=LIST_FROM(i)
      ENDDO !i
C     CALL EXITS('ILIST_COPY')
      RETURN
      END


