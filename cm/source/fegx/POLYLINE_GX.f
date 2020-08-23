      SUBROUTINE POLYLINE_GX(INDEX,NTPTS,R_XPTS,R_YPTS,ERROR,*)

C#### Subroutine: POLYLINE_GX
C###  Description:
C###    Draws a polyline.

      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='POLYLINE_GX')
      INCLUDE 'gx.inc'
!     Parameter List
      INTEGER INDEX,NTPTS
      REAL R_XPTS(NTPTS),R_YPTS(NTPTS)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ERR,LNCOL,LNTYPE

      CALL ENTERS(ROUTINENAME,*9999)

        IF(INDEX.LT.41) THEN
          LNTYPE=MOD(INDEX,4)
          IF(LNTYPE.EQ.1) THEN
            CALL SLNTYP(gxSOLID,ERR)
          ELSE IF(LNTYPE.EQ.2) THEN
            CALL SLNTYP(gxDOT,ERR)
          ELSE IF(LNTYPE.EQ.3) THEN
            CALL SLNTYP(gxDASH,ERR)
          ELSE
            CALL SLNTYP(gxDASDT,ERR)
          ENDIF
          IF(err.ne.0) THEN
            ERROR='>>Error from GX'
            GOTO 9999
          ENDIF
          IF(INDEX.GE.5.AND.INDEX.LE.8) THEN
            CALL SLNWDH(gxTHICK,ERR)
          ELSE
            CALL SLNWDH(gxTHIN,ERR)
          ENDIF
          IF(err.ne.0) THEN
            ERROR='>>Error from GX'
            GOTO 9999
          ENDIF
          IF(INDEX.LE.8) THEN                       ! Black
            CALL SLNCLR(gxBLACK,ERR)
          ELSE IF(INDEX.GE.9.AND.INDEX.LE.12) THEN  ! Red
            CALL SLNCLR(gxRED,ERR)
          ELSE IF(INDEX.GE.13.AND.INDEX.LE.16) THEN ! Green
            CALL SLNCLR(gxGREEN,ERR)
          ELSE IF(INDEX.GE.17.AND.INDEX.LE.20) THEN ! Blue
            CALL SLNCLR(gxBLUE,ERR)
          ELSE IF(INDEX.GE.21.AND.INDEX.LE.24) THEN ! Cyan
            CALL SLNCLR(gxCYAN,ERR)
          ELSE IF(INDEX.GE.25.AND.INDEX.LE.28) THEN ! Yellow
            CALL SLNCLR(gxYELLW,ERR)
          ELSE IF(INDEX.GE.29.AND.INDEX.LE.32) THEN ! White
            CALL SLNCLR(gxWHITE,ERR)
          ELSE IF(INDEX.GE.33.AND.INDEX.LE.36) THEN ! Light Blue
            CALL SLNCLR(gxLBLUE,ERR)
          ELSE IF(INDEX.GE.37.AND.INDEX.LE.40) THEN ! Grey
            CALL SLNCLR(gxGREY,ERR)
          ENDIF
          IF(err.ne.0) THEN
            ERROR='>>Error from GX'
            GOTO 9999
          ENDIF
        ELSE
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
          LNCOL=gxSPOFF+INT((INDEX-40.0)/(256.0-40.0)*gxNSPECT)
          CALL SLNCLR(LNCOL,ERR)
          IF(err.ne.0) THEN
            ERROR='>>Error from GX'
            GOTO 9999
          ENDIF
        ENDIF

        CALL DRWPLN(NTPTS,R_XPTS,R_YPTS,ERR)
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


