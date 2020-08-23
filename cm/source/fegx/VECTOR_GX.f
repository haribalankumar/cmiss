      SUBROUTINE VECTOR_GX(INDEX,CENTRE,DIRECTION,ERROR,*)

C#### Subroutine: POLYLINE_GX
C###  Description:
C###    Draws a polyline.

      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='VECTOR_GX')
      INCLUDE 'b00.cmn'
      INCLUDE 'gx.inc'
!     Parameter List
      INTEGER INDEX
      REAL CENTRE(2),DIRECTION(2)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ERR,LNTYPE,LNCOL
      REAL ANGLE
      REAL LENGTH,HEADLENGTH,R_PTS(2,2)

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

C*** Draw main body

      R_PTS(1,1)=CENTRE(1)
      R_PTS(1,2)=CENTRE(2)
      R_PTS(2,1)=CENTRE(1)+DIRECTION(1)
      R_PTS(2,2)=CENTRE(2)+DIRECTION(2)
      CALL DRWPLN(2,R_PTS(1,1),R_PTS(1,2),ERR)
      IF(err.ne.0) THEN
        ERROR='>>Error from GX'
        GOTO 9999
      ENDIF

C*** Draw head

      LENGTH=SQRT(DIRECTION(1)**2+DIRECTION(2)**2)
      IF(ABS(LENGTH).GE.1.0e-6) THEN
        IF(DIRECTION(1).EQ.0.0) THEN
          IF(DIRECTION(2).EQ.0.0) THEN
            ANGLE=0.0
          ELSE IF(DIRECTION(2).GT.0.0) THEN
            ANGLE=REAL(PI)
          ELSE
            ANGLE=REAL(-PI)
          ENDIF
        ELSE
          ANGLE=ATAN2(DIRECTION(2),DIRECTION(1))
        ENDIF
        HEADLENGTH=0.10*LENGTH
        ANGLE=ANGLE-REAL(PI-PI/12.0d0)
        R_PTS(1,1)=R_PTS(2,1)+HEADLENGTH*COS(ANGLE)
        R_PTS(1,2)=R_PTS(2,2)+HEADLENGTH*SIN(ANGLE)
        CALL DRWPLN(2,R_PTS(1,1),R_PTS(1,2),ERR)
        IF(err.ne.0) THEN
          ERROR='>>Error from GX'
          GOTO 9999
        ENDIF
        ANGLE=ANGLE+REAL(PI/6.0d0)
        R_PTS(1,1)=R_PTS(2,1)+HEADLENGTH*COS(ANGLE)
        R_PTS(1,2)=R_PTS(2,2)+HEADLENGTH*SIN(ANGLE)
        CALL DRWPLN(2,R_PTS(1,1),R_PTS(1,2),ERR)
        IF(err.ne.0) THEN
          ERROR='>>Error from GX'
          GOTO 9999
        ENDIF
      ENDIF

      CALL EXITS(ROUTINENAME)
      RETURN
 9999 CALL ERRORS(ROUTINENAME,ERROR)
      CALL EXITS(ROUTINENAME)
      RETURN 1
      END


