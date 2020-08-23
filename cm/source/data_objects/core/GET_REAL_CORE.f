      SUBROUTINE GET_REAL_CORE(REALS_PTR,ALLOCATED_REALS,
     &  INDEX,VALUE,ERROR,*)

C#### Subroutine: GET_REAL_CORE
C###  Description:
C###    The core routine for getting a real value from an array.

      IMPLICIT NONE

!     Parameter List
      INTEGER ALLOCATED_REALS,INDEX
      REAL*8 VALUE
      REAL*8 REALS_PTR(ALLOCATED_REALS)
      CHARACTER ERROR*(*)

      CALL ENTERS('GET_REAL_CORE',*9999)

      VALUE=REALS_PTR(INDEX)

      CALL EXITS('GET_REAL_CORE')
      RETURN
 9999 CALL ERRORS('GET_REAL_CORE',ERROR)
      CALL EXITS('GET_REAL_CORE')
      RETURN 1
      END


