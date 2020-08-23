      SUBROUTINE CLOSE_SEGMENT(ISEGNUM,iw,ERROR,*)

C#### Subroutine: CLOSE_SEGMENT
C###  Description:
C###    CLOSE_SEGMENT closes graphics segment ISEGNUM.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'disp00.cmn'
!     Parameter List
      INTEGER ISEGNUM,iw,ERR
      CHARACTER ERROR*(*)

      CALL ENTERS('CLOSE_SEGMENT',*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,
     '    '('' iw='',I3,'' IWKS(iw)='',I2,'' IWKT(iw)='',I2,'
     '    //''' IWKG(iw)='',I2)') iw,IWKS(iw),IWKT(iw),IWKG(iw)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' ISEGNUM='',I5)') ISEGNUM
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      IF(IWKT(iw).EQ.1) THEN !GKS
        IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
          CALL CLOBJT(ERR)
          IF(err.ne.0) THEN
            ERROR='>>Error from GX'
            GOTO 9999
          ENDIF
        ENDIF
        IF(DOP) THEN
          WRITE(OP_STRING,
     '      '('' Close segment isegnum='',I5,'' on iw='',I2)')
     '      ISEGNUM,iw
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
      ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS
      ENDIF

      CALL EXITS('CLOSE_SEGMENT')
      RETURN
 9999 CALL ERRORS('CLOSE_SEGMENT',ERROR)
      CALL EXITS('CLOSE_SEGMENT')
      RETURN 1
      END


