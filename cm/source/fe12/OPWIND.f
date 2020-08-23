      SUBROUTINE OPWIND(FULL,ERROR,*)
      
C#### Subroutine: OPWIND
C###  Description:
C###    OPWIND lists window characteristics.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'graf00.cmn'
      INCLUDE 'inout00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
      LOGICAL FULL
!     Local Variables
      INTEGER IBEG,IEND,IW,noiw

      CALL ENTERS('OPWIND',*9999)

      DO noiw=1,IWKDEF(0)
        FORMAT=' '
        IW=IWKDEF(noiw)
        IF(IWKT(IW).EQ.1) THEN      !GKS workstation
          IF(IWKS(IW).EQ.0) THEN
            FORMAT='('' Workstation   (GKS) '',I2,'' is not open.'')'
          ELSE IF(IWKS(IW).EQ.1) THEN
            FORMAT='('' Workstation   (GKS) '',I2,'' is open.'')'
          ELSE IF(IWKS(IW).EQ.2) THEN
            FORMAT='('' Workstation   (GKS) '',I2,'
     '        //''' is open & active.'')'
          ENDIF
        ELSE IF(IWKT(IW).EQ.2) THEN !PHIGS workstation
          IF(IWKS(IW).EQ.0) THEN
            FORMAT='('' Workstation (PHIGS) '',I2,'' is not open.'')'
          ELSE IF(IWKS(IW).EQ.1) THEN
            FORMAT='('' Workstation (PHIGS) '',I2,'' is open.'')'
          ELSE IF(IWKS(IW).EQ.2) THEN
            FORMAT='('' Workstation (PHIGS) '',I2,'
     '        //''' is open & active.'')'
          ENDIF
        ENDIF

        CALL STRING_TRIM(FORMAT,IBEG,IEND)
        IF(IWKG(IW).EQ.0) THEN !iw is non-graphics window
          FORMAT=FORMAT(IBEG:IEND-1)//','' Non-graphics output'')'
        ELSE IF(IWKG(IW).EQ.1) THEN !iw is graphics output window
          FORMAT=FORMAT(IBEG:IEND-1)//','' Graphics output'')'
        ENDIF
        WRITE(OP_STRING,FORMAT) IW
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDDO

      WRITE(OP_STRING,'(/'' Xmin='',E13.4,'' Xmax='',E13.4)')
     '  XMIN,XMAX
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'( '' Ymin='',E13.4,'' Ymax='',E13.4)')
     '  YMIN,YMAX
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'( '' Zmin='',E13.4,'' Zmax='',E13.4)')
     '  ZMIN,ZMAX
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,
     '  '(/'' Window width='',F5.1,''% of screen'')') SCREEN_WIDTH
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      IF(FULL) THEN
        WRITE(OP_STRING,'(/'' # open workstations IWKDEF(0)='',I2)')
     '    IWKDEF(0)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO noiw=1,IWKDEF(0)
          IW=IWKDEF(noiw)
          WRITE(OP_STRING,'('' Window # IWKDEF('',I2,'')='',I2)')
     '      noiw,IWKDEF(noiw)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''   Window type   IWKT('',I2,'
     '      //''')='',I1)') IW,IWKT(IW)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''   Open status   IWKS('',I2,'
     '      //''')='',I1)') IW,IWKS(IW)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''   Graphics type IWKG('',I2,'
     '      //''')='',I1)') IW,IWKG(IW)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

      CALL EXITS('OPWIND')
      RETURN
 9999 CALL ERRORS('OPWIND',ERROR)
      CALL EXITS('OPWIND')
      RETURN 1
      END


