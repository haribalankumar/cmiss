      SUBROUTINE SLEEPER(IDELAY)

C#### Subroutine: SLEEPER
C###  Description:
C###    SLEEPER pauses for an integer period of seconds.

      IMPLICIT NONE
!     Parameter List
      INTEGER IDELAY

      CALL SLEEP(IDELAY)

      RETURN
      END


