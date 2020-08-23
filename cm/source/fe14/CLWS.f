      SUBROUTINE CLWS(ERROR,*)

C#### Subroutine: CLWS
C###  Description:
C###    CLWS deactivates and closes all workstations.

      IMPLICIT NONE
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'disp00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ERR,iw,noiw

      CALL ENTERS('CLWS',*9999)
      DO noiw=1,IWKDEF(0)
        iw=IWKDEF(noiw)
        IF(IWKT(iw).EQ.1) THEN      !GKS
          IF(IWKS(iw).EQ.2) THEN
            IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
              CALL UNSWIN(ERR)
              IF(err.ne.0) THEN
                ERROR='>>Error from GX'
                GOTO 9999
              ENDIF
            ENDIF
            IWKS(iw)=1
          ENDIF
          IF(IWKS(iw).EQ.1) THEN
            CALL CLOSE_WS(iw,ERROR,*9999)
          ENDIF
        ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS
        ENDIF
      ENDDO
      IWKDEF(0)=0

      CALL EXITS('CLWS')
      RETURN
 9999 CALL ERRORS('CLWS',ERROR)
      CALL EXITS('CLWS')
      RETURN 1
      END


