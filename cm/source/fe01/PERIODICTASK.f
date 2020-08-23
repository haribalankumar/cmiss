      SUBROUTINE PERIODICTASK

C#### Subroutine: PERIODICTASK
C###  Description:
C###    Handles periodic tasks in order to keep cmiss running properly

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
!     Local Variables
      INTEGER ERR

C GMH 15/11/95 making call conditional
      IF(USE_GRAPHICS.EQ.1) THEN
        CALL GXWAIT(0.0,ERR) !Update graphics
      ENDIF

      RETURN
      END


