      SUBROUTINE FILL_LINE_GX(INDEX,NTPTS,R_PTS,ERROR,*)

C#### Subroutine: FILL_LINE_GX
C###  Description:
C###    Draws a polyline.
C###    It calculates the colours better when drawing 1-d fields.

      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='FILL_LINE_GX')
      INCLUDE 'gx.inc'
!     Parameter List
      INTEGER INDEX,NTPTS
      REAL R_PTS(NTPTS,2)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ERR,LNCOL

      CALL ENTERS(ROUTINENAME,*9999)

      CALL SLNTYP(gxSOLID,ERR)
      IF(err.ne.0) THEN
        ERROR='>>Error from GX'
        GOTO 9999
      ENDIF
      CALL SLNWDH(gxTHICK,ERR)
      IF(err.ne.0) THEN
        ERROR='>>Error from GX'
        GOTO 9999
      ENDIF
      LNCOL=gxSPOFF+INT((INDEX-17.0)/(249.0-17.0)*gxNSPECT)
      CALL SLNCLR(LNCOL,ERR)
      IF(err.ne.0) THEN
        ERROR='>>Error from GX'
        GOTO 9999
      ENDIF

      CALL DRWPLN(NTPTS,R_PTS(1,1),R_PTS(1,2),ERR)
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


