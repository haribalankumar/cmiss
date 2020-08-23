      SUBROUTINE GET_EQUATION_TOTAL(TOTAL,ERROR,*)

C#### Subroutine: GET_EQUATION_TOTAL
C###  Description:
C###    Gets the total number of equations.

      IMPLICIT NONE
      INCLUDE 'equation00.cmn'

!     Parameter List
      INTEGER TOTAL
      CHARACTER ERROR*(*)

      CALL ENTERS('GET_EQUATION_TOTAL',*9999)

      TOTAL=TOTAL_EQUATION_LABELS
      
      CALL EXITS('GET_EQUATION_TOTAL')
      RETURN
 9999 CALL ERRORS('GET_EQUATION_TOTAL',ERROR)
      CALL EXITS('GET_EQUATION_TOTAL')
      RETURN 1
      END


