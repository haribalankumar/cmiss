      SUBROUTINE GET_CONSTANT_TOTAL(TOTAL,ERROR,*)

C#### Subroutine: GET_CONSTANT_TOTAL
C###  Description:
C###    Gets the total number of constant labels.

      IMPLICIT NONE
      INCLUDE 'constant00.cmn'

!     Parameter List
      INTEGER TOTAL
      CHARACTER ERROR*(*)

      CALL ENTERS('GET_CONSTANT_TOTAL',*9999)

      TOTAL=TOTAL_CONSTANT_LABELS
      
      CALL EXITS('GET_CONSTANT_TOTAL')
      RETURN
 9999 CALL ERRORS('GET_CONSTANT_TOTAL',ERROR)
      CALL EXITS('GET_CONSTANT_TOTAL')
      RETURN 1
      END


