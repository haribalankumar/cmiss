      SUBROUTINE FIND(IUNIT,STRING,LINE0,LINE1,line,FOUND,ERROR,*)

C#### Subroutine: FIND
C###  Description:
C###    FIND locates the first occurrence of a character string in the
C###    direct access file connected to unit number IUNIT.
C**** If the string is found between record numbers LINE0 and LINE1,
C**** inclusive, the logical variable FOUND is returned with the value
C**** .TRUE., LINE is the record number at which the string was found,
C**** and the file pointer is left at record LINE. If the string is not
C**** found or a read beyond the end of the file is attempted, then
C**** FOUND = .FALSE. and LINE = the last record read.
C**** Programme execution is terminated if IUNIT is not opened.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER line,IUNIT,LINE0,LINE1
      CHARACTER ERROR*(*),STRING*(*)
      LOGICAL FOUND
!     Local Variables
      INTEGER IBEG,IEND,INCR,IOSTAT,RECL
      LOGICAL EXIST,OPENED
      CHARACTER ACCESS*10,CHAR*2,CLINE*500,FORM*11

      CALL ENTERS('FIND',*9999)
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' String: '',A)') STRING
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF
      INQUIRE(UNIT=IUNIT,EXIST=EXIST,OPENED=OPENED)
      IF(EXIST) THEN
        IF(OPENED) THEN
          INQUIRE(UNIT=IUNIT,ACCESS=ACCESS,FORM=FORM)
          IF(ACCESS.EQ.'DIRECT') THEN
            IF(FORM.EQ.'FORMATTED') THEN
              INQUIRE(UNIT=IUNIT,RECL=RECL)
              IF(LINE1.GE.LINE0) THEN
                INCR= 1
              ELSE
                INCR=-1
              ENDIF
              DO line=LINE0,LINE1,INCR
                READ(UNIT=IUNIT,FMT='(A)',REC=line,IOSTAT=IOSTAT,
     '            ERR=9999) CLINE(1:RECL)
                IF(IOSTAT.EQ.0) THEN
                  IF(INDEX(CLINE(1:RECL),STRING).NE.0) THEN
                    FOUND=.TRUE.
                    ERROR=' '
                    GOTO 9998
                  ENDIF
                ELSE
                  FOUND=.FALSE.
                  ERROR=' '
                  GOTO 9998
                ENDIF
              ENDDO
              FOUND=.FALSE.
              GOTO 9998
            ELSE
              WRITE(CHAR,'(I2)') IUNIT
              CALL STRING_TRIM(CHAR,IBEG,IEND)
              ERROR=' File '//CHAR(IBEG:IEND)//
     '          ' is not connected for formatted i/o'
              GOTO 9999
            ENDIF
          ELSE
            WRITE(CHAR,'(I2)') IUNIT
            CALL STRING_TRIM(CHAR,IBEG,IEND)
            ERROR=' File '//CHAR(IBEG:IEND)//
     '        ' is not connected for direct access'
            GOTO 9999
          ENDIF
        ELSE
          WRITE(CHAR,'(I2)') IUNIT
          CALL STRING_TRIM(CHAR,IBEG,IEND)
          ERROR=' File '//CHAR(IBEG:IEND)//' is not open'
          GOTO 9999
        ENDIF
      ELSE
        WRITE(CHAR,'(I2)') IUNIT
        CALL STRING_TRIM(CHAR,IBEG,IEND)
        ERROR=' File '//CHAR(IBEG:IEND)//' does not exist'
        GOTO 9999
      ENDIF

 9998 CALL EXITS('FIND')
      RETURN
 9999 CALL ERRORS('FIND',ERROR)
      CALL EXITS('FIND')
      RETURN 1
      END


