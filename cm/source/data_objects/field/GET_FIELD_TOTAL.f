      SUBROUTINE GET_FIELD_TOTAL(TOTAL,ERROR,*)

C#### Subroutine: GET_FIELD_TOTAL
C###  Description:
C###    Gets the total number of fields.

      IMPLICIT NONE
      INCLUDE 'field00.cmn'

!     Parameter List
      INTEGER TOTAL
      CHARACTER ERROR*(*)

      CALL ENTERS('GET_FIELD_TOTAL',*9999)

      TOTAL=TOTAL_FIELD_LABELS
      
      
      CALL EXITS('GET_FIELD_TOTAL')
      RETURN
 9999 CALL ERRORS('GET_FIELD_TOTAL',ERROR)
      CALL EXITS('GET_FIELD_TOTAL')
      RETURN 1
      END


