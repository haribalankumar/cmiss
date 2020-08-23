      SUBROUTINE POLYMARKER_GX(INDEX,NTPTS,R_XPTS,R_YPTS,ERROR,*)

C#### Subroutine: POLYMARKER_GX
C###  Description:
C###    Draws a polymarker.

      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='POLYMARKER_GX')
      INCLUDE 'gx.inc'
!     Parameter List
      INTEGER INDEX,NTPTS
      REAL R_XPTS(NTPTS),R_YPTS(NTPTS)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ERR,PMCOL,PMTYPE

      CALL ENTERS(ROUTINENAME,*9999)

        IF(INDEX.LT.41) THEN
          PMTYPE=MOD(INDEX,4)
          IF(PMTYPE.EQ.1) THEN
            CALL SMKTYP(gxPLUS,ERR)
          ELSE IF(PMTYPE.EQ.2) THEN
            CALL SMKTYP(gxATRSK,ERR)
          ELSE IF(PMTYPE.EQ.3) THEN
            CALL SMKTYP(gxCRCLE,ERR)
          ELSE
            CALL SMKTYP(gxPOINT,ERR)
          ENDIF
          IF(err.ne.0) THEN
            ERROR='>>Error from GX'
            GOTO 9999
          ENDIF
          IF(INDEX.LE.8) THEN                       ! Black
            CALL SMKCLR(gxBLACK,ERR)
          ELSE IF(INDEX.GE.9.AND.INDEX.LE.12) THEN  ! Red
            CALL SMKCLR(gxRED,ERR)
          ELSE IF(INDEX.GE.13.AND.INDEX.LE.16) THEN ! Green
            CALL SMKCLR(gxGREEN,ERR)
          ELSE IF(INDEX.GE.17.AND.INDEX.LE.20) THEN ! Blue
            CALL SMKCLR(gxBLUE,ERR)
          ELSE IF(INDEX.GE.21.AND.INDEX.LE.24) THEN ! Cyan
            CALL SMKCLR(gxCYAN,ERR)
          ELSE IF(INDEX.GE.25.AND.INDEX.LE.28) THEN ! Yellow
            CALL SMKCLR(gxYELLW,ERR)
          ELSE IF(INDEX.GE.29.AND.INDEX.LE.32) THEN ! White
            CALL SMKCLR(gxWHITE,ERR)
          ELSE IF(INDEX.GE.33.AND.INDEX.LE.36) THEN ! Light Blue
            CALL SMKCLR(gxLBLUE,ERR)
          ELSE IF(INDEX.GE.37.AND.INDEX.LE.40) THEN ! Grey
            CALL SMKCLR(gxGREY,ERR)
          ENDIF
          IF(err.ne.0) THEN
            ERROR='>>Error from GX'
            GOTO 9999
          ENDIF
        ELSE
          CALL SMKTYP(gxPLUS,ERR)
          IF(err.ne.0) THEN
            ERROR='>>Error from GX'
            GOTO 9999
          ENDIF
          PMCOL=gxSPOFF+INT((INDEX-40.0)/(256.0-40.0)*gxNSPECT)
          CALL SMKCLR(PMCOL,ERR)
          IF(err.ne.0) THEN
            ERROR='>>Error from GX'
            GOTO 9999
          ENDIF
        ENDIF

        CALL DRWPMK(NTPTS,R_XPTS,R_YPTS,ERR)
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


