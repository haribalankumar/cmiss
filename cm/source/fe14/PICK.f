      SUBROUTINE PICK(iw,MODE,INSTAT,IPICKRET,IPICKID,ERROR,*)

C#### Subroutine: PICK
C###  Description:
C###    PICK calls GKS pick.

C**** MODE can be 'REQUEST' or 'EVENT'
C**** IPICKRET is returned segment on exit if in request mode
C**** IPCKID   is returned pick identifier if in request mode
C**** INSTAT   ir returned 1 if successful input triggered

      IMPLICIT NONE
      INCLUDE 'disp00.cmn'
      INCLUDE 'echo00.cmn'
!     Parameter List
      INTEGER INSTAT,IPICKID,IPICKRET,iw
      CHARACTER ERROR*(*),MODE*(*)
!     Local variables
      INTEGER ERR

      CALL ENTERS('PICK',*9999)

c     CALL INPUT_MODE(iw,ID,'PICK',MODE,ERROR,*9999)
      IF(iw.ne.3) THEN !gks
!       Set pick area to full screen

        IF(MODE(1:5).EQ.'EVENT') THEN
        ELSE IF(MODE(1:7).EQ.'REQUEST') THEN
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            ECAREA(1)=0.0
            ECAREA(2)=0.0
            ECAREA(3)=0.0
            ECAREA(4)=0.0
            CALL FPICK(ECAREA,IPICKRET,INSTAT,ERR)
            IF(err.ne.0) THEN
              ERROR='>>Error from GX'
              GOTO 9999
            ENDIF
            IF(INSTAT.EQ.0) THEN
              INSTAT=1
            ELSE
              INSTAT=0
            ENDIF
            IPICKID=IPICKRET
c         ELSE IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
c           ECAREA(1)=0.0
c           ECAREA(2)=XDISP
c           ECAREA(3)=0.0
c           ECAREA(4)=YDISP
c           CALL GKS_INPK(iw,ID,GNPICK,1,1,1,0,' ',ERROR,*9999)
c           CALL GKS_RQPK(iw,ID,INSTAT,IPICKRET,IPICKID,ERROR,*9999)
c           IF(INSTAT.EQ.GOK) THEN
c             INSTAT=1
c           ELSE
c             INSTAT=0
c           ENDIF
          ENDIF
        ENDIF
      ENDIF

      CALL EXITS('PICK')
      RETURN

 9999 CALL ERRORS('PICK',ERROR)
      CALL EXITS('PICK')
      RETURN 1
      END


