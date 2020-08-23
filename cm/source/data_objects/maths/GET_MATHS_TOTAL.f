      SUBROUTINE GET_MATHS_TOTAL(TOTAL,ERROR,*)

C#### Subroutine: GET_MATHS_TOTAL
C###  Description:
C###    Gets the total number of maths objects.

      IMPLICIT NONE
      INCLUDE 'maths00.cmn'

!     Parameter List
      INTEGER TOTAL
      CHARACTER ERROR*(*)

      CALL ENTERS('GET_MATHS_TOTAL',*9999)

      TOTAL=TOTAL_MATHS_LABELS
      
      
      CALL EXITS('GET_MATHS_TOTAL')
      RETURN
 9999 CALL ERRORS('GET_MATHS_TOTAL',ERROR)
      CALL EXITS('GET_MATHS_TOTAL')
      RETURN 1
      END


