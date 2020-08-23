      SUBROUTINE CLOSE_WS(iw,ERROR,*)

C#### Subroutine: CLOSE_WS
C###  Description:
C###    CLOSE_WS closes graphics workstation iw.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'disp00.cmn'
!     Parameter List
      INTEGER iw
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ERR

      CALL ENTERS('CLOSE_WS',*9999)

      IF(IWKT(iw).EQ.1) THEN      !GKS window
        IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
          CALL CLWIN(IW,ERR)
          IF(err.ne.0) THEN
            ERROR='>>Error from GX'
            GOTO 9999
          ENDIF
        ENDIF
      ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS window
      ENDIF

      IF(DOP) THEN
        WRITE(OP_STRING,
     '    '('' Unit='',I2,'' has been closed.'')') iw
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      IWKS(iw)=0

      CALL EXITS('CLOSE_WS')
      RETURN
 9999 CALL ERRORS('CLOSE_WS',ERROR)
      CALL EXITS('CLOSE_WS')
      RETURN 1
      END


