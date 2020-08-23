      SUBROUTINE ISORT(n,IDATA)

C#### Subroutine: ISORT
C###  Description:
C###    ISORT sorts N integer IDATA values into a non-decreasing
C###    sequence using ISHELLSORT if N<50 or IHEAPSORT if N>50 as
C###    recommended by numerical recipes.

      IMPLICIT NONE

!     Parameter List
      INTEGER IDATA(*),N
!     Local Variables

      IF(N.LE.50) THEN
        CALL ISHELLSORT(N,IDATA)
      ELSE
        CALL IHEAPSORT(N,IDATA)
      ENDIF

      RETURN
      END


