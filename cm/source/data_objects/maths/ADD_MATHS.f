      SUBROUTINE ADD_MATHS(LABEL,ERROR,*)

C#### Subroutine: ADD_MATHS
C###  Description:
C###    Adds a maths object.

      IMPLICIT NONE
      INCLUDE 'maths00.cmn'

!     Parameter List
      CHARACTER LABEL*(*),ERROR*(*)

      CALL ENTERS('ADD_MATHS',*9999)

C DMAL 21 JULY 2004: These routine have been written but not used.
C This block is to satisfy thefortran check errors.      
      IF(.FALSE.)THEN
        CALL PRINT_MATHS_VALUES(%VAL(0),ERROR,*9999)
      ENDIF

      CALL SET_LABEL(MATHS_LABELS_PTR,TOTAL_MATHS_LABELS,
     &  ALLOCATED_MATHS_LABELS,MATHS_LABELS_LEN,LABEL,ERROR,*9999)

      CALL EXITS('ADD_MATHS')
      RETURN
 9999 CALL ERRORS('ADD_MATHS',ERROR)
      CALL EXITS('ADD_MATHS')
      RETURN 1
      END


