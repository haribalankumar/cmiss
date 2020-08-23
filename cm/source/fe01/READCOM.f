      SUBROUTINE READCOM(FILENAME,QUIT,ERROR,*)

C#### Subroutine: READCOM
C###  Description:
C###    <HTML>
C###    READCOM reads command (*.com) files.
C###    <PRE>
C###    PAUSE on    pause before every exec command
C###    PAUSE off   do not pause
C###    </PRE></HTML>

C***  Command file statements:
C***    !  or # indicates comment to be printed
C***    !! or ##   "        "    "  "  ignored
C***    DO          DO loop
C***    ENDDO       end of DO loop
C***    IF    <-|   simple conditional IF statement
C***    ELSEIF  |   eg: IF I < 100
C***    ELSE    |   or: IF I=1 [NB must use spaces for all ops except =]
C***    ENDIF <-|
C***    PAUSE on    pause before every exec command
C***    PAUSE off   do not pause
C***    </PRE></HTML>
C****  Note: When DO is encountered LOOP_NAME & LIST of indices are
C****        defined
C****        At top of loop LOOP_NAME is assigned to next member of LIST
C****        (This is signified by START which is reset to .true. when
C****        an ENDDO is encountered - provided LIST has not been
C****        exhausted)

      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='READCOM')
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
C      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'gtstr00.cmn'
      INCLUDE 'inout00.cmn'
!     Parameter List
      CHARACTER FILENAME*(*),ERROR*(*)
      LOGICAL QUIT
!     Local Variables
      CHARACTER COMMENT
      PARAMETER (COMMENT='#')
      INTEGER COMLINE_BEG,COMLINE_END,ERR,NOQUES,ICHAR,
     '  INFO,IOSTAT,POSITION
      CHARACTER CHAR,COMMAND_LINE*(MXCH)
      LOGICAL ABBREV,CONTINUE,FILEIP,PAUSE

      CALL ENTERS(ROUTINENAME,*9999)

C KAT 25/2/00: Opening and closing COM_UNIT in the same routine.
      COM_UNIT=COM_UNIT+1
C      IEND=LEN(FILENAME)
C      IF(IEND.GT.4.AND.FILENAME(IEND-3:IEND).EQ.'.com') IEND=IEND-4
C KAT 28/2/00: Doesn't work with g77
C      CALL OPENF(COM_UNIT,'DISK',FILENAME(:IEND)//'.com',
C     '  'OLD','DIRECT','FORMATTED',132,ERROR,*9998)
      CALL OPENF(COM_UNIT,'DISK',FILENAME,
     '  'OLD','DIRECT','FORMATTED',MXCH,ERROR,*9997)
C KAT 2/3/00: First comfile read stays open in OPEN_COM_UNIT.
      IF(NEST.EQ.0.AND.IREC_COMFILE(0).LT.0.AND.FILENAME.NE.'cmiss.com')
     '  THEN
C       Must change OPEN_COM_UNIT in case `open com' is executed in comfile.
        OPEN_COM_UNIT=COM_UNIT
        IREC_COMFILE(0)=1
      ENDIF
      NEST=NEST+1

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999
      PAUSE=.FALSE.
      IREC_COMFILE(NEST)=1
      CONTINUE=.TRUE.
      DO WHILE(CONTINUE)

        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' irec='',I3)') IREC_COMFILE(NEST)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF

C CS 26/2/99 Initialise string incase the last line in the com
C file is empty
        COMMAND_LINE=' '

        READ(COM_UNIT,FMT='(A)',REC=IREC_COMFILE(NEST),IOSTAT=IOSTAT)
     '    COMMAND_LINE

        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,
     '      '('' iostat='',I4,'' COMMAND_LINE(1:MXCH)='',A)') iostat,
     '      COMMAND_LINE(1:MXCH)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF

C KAT 22Dec99: Not executing command if there is an error condition
C       At end of file, there may be a command without a carriage
C       return.  Execute it anyway.

C Portland linux comp uses iostat 253 for reading past end of records
C        IF(IOSTAT.EQ.253) IOSTAT=-1


C DPN 21-May-2002 - From GINOUT()....
C -------------------------------------------
C KAT 18Aug00:
C Note 9.23 in X3J3/96-007 draft standard says "An end-of-file condition
C may occur only during execution of a sequential input file".  We
C should therefore get a +ve (error) IOSTAT if the record didn't exist.
C However SGI's f77 gives a -1 to indicate end-of-file.  It would be
C nice to check that we don't have another error condition but
C error codes differ between compilers.  Therefore we assume that
C a non-zero IOSTAT means end-of-file.
C --------------------------------------------
C
C The IBM compiler will only set IOSTAT to -1 to indicate an end-of-file
C if the END= specifier is also set on the above read, but the Fortran
C standard states that it is illegal to specify both END= and REC=
C specifiers, which breaks the Linux compilation. Therefore, simply check
C for the known IBM error here and manually set to -1.

C The Intel compiler sets IOSTAT to 158 for record number out of range.
C Intel compiler 8.0 uses iostat 36 for accessing non-existant records.

C We created this file ourselves so we could calculate how many records
C there are when we create it.

C         IF(IOSTAT.EQ.36) IOSTAT=-1

C         IF(IOSTAT.EQ.1.OR.IOSTAT.EQ.158) IOSTAT=-1

        IF(IOSTAT.GT.0) IOSTAT=-1

        IF(IOSTAT.LE.0) THEN
CC!!! CS 25/2/99 Changed to work under linux too
CC        IF(IOSTAT.EQ.0) THEN
C        IF(IOSTAT.GE.-1) THEN
C          IF(IOSTAT.EQ.-1) CONTINUE=.FALSE.

C KAT 2/3/00: Setting IREC_COMFILE(0) here in case there is a jump to
C         the fatal signal handler.
          IF(IOSTAT.EQ.0) IREC_COMFILE(NEST)=IREC_COMFILE(NEST)+1
          IF(COM_UNIT.EQ.OPEN_COM_UNIT) IREC_COMFILE(0)=IREC_COMFILE(1)

          CALL STRING_TRIM(COMMAND_LINE,COMLINE_BEG,COMLINE_END)
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' COMMAND_LINE='',A)')
     '        COMMAND_LINE(COMLINE_BEG:COMLINE_END)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          POSITION=INDEX(COMMAND_LINE,COMMENT)
          IF(POSITION.NE.0.AND.POSITION.LT.COMLINE_END) THEN
            IF(COMMAND_LINE(POSITION+1:POSITION+1).EQ.COMMENT)
     '        COMLINE_END=POSITION-1
          ENDIF

          IF(COMMAND_LINE(COMLINE_BEG:COMLINE_END).NE.' ') THEN

            IF(PAUSE) THEN
              FORMAT='($,'' Continue [yes]?  '',A)'
              IOTYPE=0
              CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,AYES,CDATA,CDEFLT,
     '          ICHAR,IDATA,IDEFLT,0,1,LDATA,LDEFLT,
     '          RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'(1X,A)') ' ADATA(1)=',ADATA(1)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF
              IF(ADATA(1).EQ.'N') THEN
                CONTINUE=.FALSE.
              ENDIF
            ENDIF
            IF(CONTINUE) THEN
              !Echo the command
              IF(ECHO_RAW_COM) THEN
                OP_STRING(1)='> '//COMMAND_LINE(:COMLINE_END)
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C               This is so that the command is written to log files
C               before a fatal signal is received.  This function is not
C               in the Fortran standard so perhaps it may need to be
C               replaced with a call to a C function to flush STDOUT.
                CALL FLUSH(IOOP)
              ENDIF

              IF(ABBREV(COMMAND_LINE(COMLINE_BEG:COMLINE_BEG+4),
     '          'PAUSE',5)) THEN
                IF(ABBREV(COMMAND_LINE(COMLINE_BEG+6:
     '            COMLINE_BEG+8),'OFF',2)) THEN
                  PAUSE=.FALSE.
                ELSE !IF(ABBREV(COMMAND_LINE(COMLINE_BEG+6:COMLINE_BEG+7),
!     '            'ON',2)) THEN
                  PAUSE=.TRUE.
                ENDIF
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                  WRITE(OP_STRING,'(1X,A)') ' PAUSE=',PAUSE
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                ENDIF

              ELSE
                IF(NEST.EQ.1) THEN !comfile run from cm argument
                  CALL STORE_COMMAND(COMMAND_LINE(:COMLINE_END),ERR,
     '              ERROR,*9999)
                ENDIF
C ***           Parse command
                CALL INTERPRET_COMMAND_LINE(QUIT,
     '            COMMAND_LINE(:COMLINE_END),ERR)
                IF(ERR.NE.0) THEN
                  ERROR=' '
                  GO TO 9999
                ENDIF !ERR
                IF(QUIT) CONTINUE=.FALSE.
              ENDIF !PAUSE
            ENDIF !CONTINUE

          ENDIF !blank line

          IF(IOSTAT.LT.0) CONTINUE=.FALSE. !end of file
C KAT 22Dec99: I think this is less platform dependent.
        ELSE
          NEST=NEST-1
          CALL FLAG_ERROR(0,' ')
          CALL WRITE_CHAR(IOER,'Command file read failed: IOSTAT=',ERR)
          CALL WRITE_INT(IOER,IOSTAT,ERR)
          CALL WRITE_CHAR(IOER,NEWLINE,ERR)
          ERROR=' '
          CALL CLOSEF(COM_UNIT,CHAR,*9997)
          GO TO 9997
C        ELSE IF(IOSTAT.EQ.36) THEN
C          CALL CLOSEF(COM_UNIT,ERROR,*9999)
C          COM_UNIT=COM_UNIT-1
C          NEST=NEST-1
        ENDIF
      ENDDO

C KAT 1Feb00: First comfile read stays open in OPEN_COM_UNIT
      NEST=NEST-1
      IF(COM_UNIT.NE.OPEN_COM_UNIT) THEN
        CALL CLOSEF(COM_UNIT,ERROR,*9997)
        COM_UNIT=COM_UNIT-1
      ENDIF

      CALL EXITS(ROUTINENAME)
      RETURN

C KAT 1Feb00: First comfile read stays open in OPEN_COM_UNIT
 9999 NEST=NEST-1
      IF(COM_UNIT.NE.OPEN_COM_UNIT) THEN
        CALL CLOSEF(COM_UNIT,CHAR,*9997)
        COM_UNIT=COM_UNIT-1
      ENDIF
      GOTO 9990
C KAT 6/6/00: Not trying to open file if it does not exist.
C 9998 IF(ERROR(1:14).EQ.'File not found'
C     '  .AND.FILENAME.EQ.'./cmiss.com') THEN
C        WRITE(OP_STRING,'('' >>no cmiss.com file defined'')')
C        CALL WRITES(IODI,OP_STRING,ERROR,*9997)
C        COM_UNIT=COM_UNIT-1
C        RETURN
C      ENDIF
 9997 COM_UNIT=COM_UNIT-1
 9990 CALL ERRORS(ROUTINENAME,ERROR)
      CALL EXITS(ROUTINENAME)
      RETURN 1
      END


