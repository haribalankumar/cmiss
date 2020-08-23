      SUBROUTINE QUIT_GRAPHICS(ERROR,*)

C#### Subroutine: QUIT_GRAPHICS
C###  Description:
C###    QUIT_GRAPHICS closes any workstations and gks or phigs.

      IMPLICIT NONE
      INCLUDE 'disp00.cmn'
      INCLUDE 'gks000.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER ERR

      CALL ENTERS('QUIT_GRAPHICS',*9999)
      IF(GKS) THEN
        CALL CLWS(ERROR,*9999)
      ENDIF
      IF(GKS) THEN
        IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
          CALL CLGRPH(ERR)
          IF(err.ne.0) THEN
            ERROR='>>Error from GX'
            GOTO 9999
          ENDIF
        ENDIF
        GKS=.FALSE.
      ENDIF

      CALL EXITS('QUIT_GRAPHICS')
      RETURN
 9999 CALL ERRORS('QUIT_GRAPHICS',ERROR)
      CALL EXITS('QUIT_GRAPHICS')
      RETURN 1
      END


