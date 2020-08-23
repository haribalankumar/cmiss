      SUBROUTINE CRASH(DUMMY,ERROR,*)

C#### Subroutine: CRASH
C###  Description:
C###    CRASH causes crash by assigning 3 parameters in call to two
C###    parameters in subroutine. Used for deliberately crashing from
C###    CMISS.

      IMPLICIT NONE
!     Parameter List
      REAL*8 DUMMY
      CHARACTER ERROR*(*)

      CALL ENTERS('CRASH',*9999)

      ERROR='crash'

      CALL EXITS('CRASH')
      RETURN
9999  CALL ERRORS('CRASH',ERROR)
      CALL EXITS('CRASH')
      RETURN 1
      END


