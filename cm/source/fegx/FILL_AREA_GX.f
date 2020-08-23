      SUBROUTINE FILL_AREA_GX(INDEX,NTPTS,R_PTS,ERROR,*)

C#### Subroutine: FILL_AREA_GX
C###  Description:
C###    Draws a fill-area.

C**** If the IBUNDLE is 0 the primitive will use the previously
C**** defined fill-area index.
C**** NOTE: If iw is 15 or 16 (postscript) IBUNDLE is reset to be black.

      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='FILL_AREA_GX')
      INCLUDE 'gx.inc'
!     Parameter List
      INTEGER INDEX,NTPTS
      REAL R_PTS(NTPTS,2)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ERR,FLCOL

      CALL ENTERS(ROUTINENAME,*9999)

      IF(INDEX.EQ.1) THEN
        CALL SFLCLR(gxBLACK,ERR)
        IF(err.ne.0) THEN
          ERROR='>>Error from GX'
          GOTO 9999
        ENDIF
      ELSE
        FLCOL=gxSPOFF+INT((INDEX-17.0)/(249.0-17.0)*gxNSPECT)
        CALL SFLCLR(FLCOL,ERR)
        IF(err.ne.0) THEN
          ERROR='>>Error from GX'
          GOTO 9999
        ENDIF
      ENDIF

      CALL FLAREA(NTPTS,R_PTS(1,1),R_PTS(1,2),ERR)
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


