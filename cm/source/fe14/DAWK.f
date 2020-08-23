      SUBROUTINE DAWK(iw,ID,ERROR,*)

C#### Subroutine: DAWK
C###  Description:
C###    DAWK deactivates workstation identified by iw.

C**** If ID=1 deferred output is released.
C**** If POSTSCRIPT is .true. (as set by call to WS_LIST)
C**** then workstation 15 (GKS) or 16 (Phigs) is deactivated also,
C**** and POSTSCRIPT is set to .false.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'disp00.cmn'
      INCLUDE 'cbwk01.cmn'
!     Parameter List
      INTEGER ID,iw
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER ERR

      CALL ENTERS('DAWK',*9999)

      IF(DOP) THEN
        WRITE(OP_STRING,
     '    '('' iw='',I3,'' IWKS(iw)='',I2,'//
     '    ''' IWKT(iw)='',I2)') iw,IWKS(iw),IWKT(iw)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(IWKT(iw).EQ.1) THEN      !GKS window
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
      ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS window
      ELSE IF(IWKT(iw).EQ.3) THEN !Frame grabber display screen
      ENDIF

      CALL EXITS('DAWK')
      RETURN
 9999 CALL ERRORS('DAWK',ERROR)
      CALL EXITS('DAWK')
      RETURN 1
      END


