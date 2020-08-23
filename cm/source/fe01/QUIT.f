      SUBROUTINE QUIT(END,ERROR)

C#### Subroutine: QUIT
C###  Description:
C###    QUIT performs a graceful end to the plotting session.

      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='QUIT')
      INCLUDE 'gks000.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
      LOGICAL END
!     Local Variables
!     INTEGER IERR,MBCLS

      CALL ENTERS(ROUTINENAME,*99)
      END=.TRUE.
      IF(GKS) CALL QUIT_GRAPHICS(ERROR,*99)

      GOTO 100
 99   CALL ERRORS(ROUTINENAME,ERROR)

 100  CALL EXITS(ROUTINENAME)
      RETURN
      END


