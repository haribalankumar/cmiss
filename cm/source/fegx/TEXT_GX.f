      SUBROUTINE TEXT_GX(INDEX,STRING,HEIGHT,R_PT,ERROR,*)

C#### Subroutine: TEXT_GX
C###  Description:
C###    Draws text.

      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='TEXT_GX')
      INCLUDE 'gx.inc'
!     Parameter List
      INTEGER INDEX
      CHARACTER STRING*(*)
      REAL HEIGHT,R_PT(2)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ERR

      CALL ENTERS(ROUTINENAME,*9999)

      IF(INDEX.LE.8) THEN                       ! Black
        CALL STXCLR(gxBLACK,ERR)
      ELSE IF(INDEX.GE.9.AND.INDEX.LE.12) THEN  ! Red
        CALL STXCLR(gxRED,ERR)
      ELSE IF(INDEX.GE.13.AND.INDEX.LE.16) THEN ! Green
        CALL STXCLR(gxGREEN,ERR)
      ELSE IF(INDEX.GE.17.AND.INDEX.LE.20) THEN ! Blue
        CALL STXCLR(gxBLUE,ERR)
      ELSE IF(INDEX.GE.21.AND.INDEX.LE.24) THEN ! Cyan
        CALL STXCLR(gxCYAN,ERR)
      ELSE IF(INDEX.GE.25.AND.INDEX.LE.28) THEN ! Yellow
        CALL STXCLR(gxYELLW,ERR)
      ELSE IF(INDEX.GE.29.AND.INDEX.LE.32) THEN ! White
        CALL STXCLR(gxWHITE,ERR)
      ELSE IF(INDEX.GE.33.AND.INDEX.LE.36) THEN ! Light Blue
        CALL STXCLR(gxLBLUE,ERR)
      ELSE IF(INDEX.GE.37.AND.INDEX.LE.40) THEN ! Grey
        CALL STXCLR(gxGREY,ERR)
      ENDIF
      IF(err.ne.0) THEN
        ERROR='>>Error from GX'
        GOTO 9999
      ENDIF
      CALL STXALN(gxCNTRE,gxCNTRE,ERR)
      IF(err.ne.0) THEN
        ERROR='>>Error from GX'
        GOTO 9999
      ENDIF
      CALL STXHGT(HEIGHT,ERR)
      IF(err.ne.0) THEN
        ERROR='>>Error from GX'
        GOTO 9999
      ENDIF
      CALL WRTTXT(R_PT(1),R_PT(2),STRING,ERR)
      IF(err.ne.0) THEN
        ERROR='>>Error from GX'
        GOTO 9999
      ENDIF

      CALL EXITS(ROUTINENAME)
      RETURN
 9999 CALL ERRORS(ROUTINENAME,ERROR)
      CALL EXITS(ROUTINENAME)
      RETURN 1
      END


