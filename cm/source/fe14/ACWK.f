      SUBROUTINE ACWK(iw,ID,ERROR,*)

C#### Subroutine: ACWK
C###  Description:
C###    ACWK activates workstation identified by iw.
C**** If ID=1 output is deferred.
C**** If iw=5 or 6 buffers are cleared.
C**** If POSTSCRIPT is .true. (as set by call to WS_LIST)
C**** then workstation 15 (GKS) or 16 (Phigs) is activated also.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'disp00.cmn'
      INTEGER iw,ID
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ERR

      CALL ENTERS('ACWK',*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,'('' iw='',I3,'' IWKS(iw)='',I2,'//
     '    ''' IWKT(iw)='',I2)') iw,IWKS(iw),IWKT(iw)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(iw.EQ.9) THEN !electrocardiographic lead
        IF(IWKS(iw).EQ.0) CALL SETUP(iw,ERROR,*9999)
      ELSE IF(iw.EQ.10.OR.iw.EQ.11) THEN !history & section windows
        IF(IWKS(iw).EQ.0) CALL SETUP(iw,ERROR,*9999)
      ELSE IF(iw.EQ.12.OR.iw.EQ.13) THEN !fibre or sheet angle windows
        IF(IWKS(iw).EQ.0) CALL SETUP(iw,ERROR,*9999)
      ENDIF

      IF(IWKS(iw).EQ.1) THEN !workstation iw is defined
        IF(IWKT(iw).EQ.1) THEN      !GKS window
          IF(DOP) THEN
            WRITE(OP_STRING,'('' Activate workstation iw='',i3)') iw
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL SETWIN(IW,ERR)
            IF(err.ne.0) THEN
              ERROR='>>Error from GX'
              GOTO 9999
            ENDIF
          ENDIF
        ELSE IF(IWKT(iw).EQ.2) THEN !Phigs window
        ENDIF
      ENDIF

      CALL EXITS('ACWK')
      RETURN
 9999 CALL ERRORS('ACWK',ERROR)
      CALL EXITS('ACWK')
      RETURN 1
      END


