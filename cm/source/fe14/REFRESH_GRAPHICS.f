      SUBROUTINE REFRESH_GRAPHICS(DELAY,ERROR,*)

C#### Subroutine: REFRESH_GRAPHICS
C###  Description:
C###    REFRESH_GRAPHICS is a buffered call to GXWAIT. DELAY is the
C###    time to wait (note that DELAY is a REAL variable).

      IMPLICIT NONE

      INCLUDE 'cbwk01.cmn'

!     Parameter List
      REAL DELAY
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ERR,IW,noiw
      LOGICAL WINDOWISOPEN

      CALL ENTERS('REFRESH_GRAPHICS',*9999)

C     MLB 20-April-2000 Must check to see that a window has been
C     opened before calling GXWAIT.
      WINDOWISOPEN=.FALSE.
      DO noiw=1,IWKDEF(0)
        IW=IWKDEF(noiw)
        IF(IWKS(IW).GT.0) WINDOWISOPEN=.TRUE.
      ENDDO
      IF(WINDOWISOPEN) THEN
        CALL GXWAIT(DELAY,ERR)
        IF(ERR.NE.0) THEN
          ERROR='>>Error from GX'
          GOTO 9999
        ENDIF
      ENDIF

      CALL EXITS('REFRESH_GRAPHICS')
      RETURN
 9999 CALL ERRORS('REFRESH_GRAPHICS',ERROR)
      CALL EXITS('REFRESH_GRAPHICS')
      RETURN 1
      END


