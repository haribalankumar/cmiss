      SUBROUTINE ASSEMBLE4(ERROR,*)

C#### Subroutine: ASSEMBLE4
C###  Description:
C###    ASSEMBLE4 assembles the global unreduced matrices GK,
C###    GM, etc. for Time dependent BEM problems.

      IMPLICIT NONE
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables

      CALL ENTERS('ASSEMBLE4',*9999)

      CALL ASSERT(.FALSE.,'Not implemented yet',ERROR,*9999)

      CALL EXITS('ASSEMBLE4')
      RETURN
 9999 CALL ERRORS('ASSEMBLE4',ERROR)
      CALL EXITS('ASSEMBLE4')
      RETURN 1
      END


