      SUBROUTINE GET_LABEL_CORE(LABELS_PTR,ALLOCATED_LABELS,
     &  INDEX,LABEL,ERROR,*)

C#### Subroutine: GET_LABEL_CORE
C###  Description:
C###    The core routing for getting a label from an array given an index.

      IMPLICIT NONE

!     Parameter List
      INTEGER ALLOCATED_LABELS,INDEX
      CHARACTER LABEL*(*),LABELS_PTR(ALLOCATED_LABELS)*(*),ERROR*(*)

      CALL ENTERS('GET_LABEL_CORE',*9999)

      LABEL=LABELS_PTR(INDEX)

      CALL EXITS('GET_LABEL_CORE')
      RETURN
 9999 CALL ERRORS('GET_LABEL_CORE',ERROR)
      CALL EXITS('GET_LABEL_CORE')
      RETURN 1
      END


