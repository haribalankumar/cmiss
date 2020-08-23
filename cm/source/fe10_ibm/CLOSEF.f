      SUBROUTINE CLOSEF(IUNIT,ERROR,*)

C#### Subroutine: CLOSEF
C###  Description:
C###    CLOSEF closes a file using the FORTRAN CLOSE command.
C###    IUNIT is the name by which the file is referred to in the
C###    program.  Direct access files are truncated at their current
C###    record.  ERROR gives diagnostics in the event of failure.
C###    If an error is detected control is returned to the statement
C###    number of the star.  ERROR is a character string.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'file00.cmn'
!     Parameter List
      INTEGER IUNIT
      CHARACTER ERROR*(*)
!     Local Variables
c      INTEGER I,IBEG,ICLEN,IEND,INTSTR(1024),IOSTAT,IREC
      INTEGER I,IBEG,IEND,IOSTAT,IREC
      CHARACTER ACCESS*10,FILE*100,FILENAME*(MXCH),
     '  FORMAT*11,ISTAT*5,LINE*132
      LOGICAL EXIST,OPENED

      CALL ENTERS('CLOSEF',*9999)
C DB.  Make sure that OPENED and EXISTs before getting more information
C      INQUIRE(UNIT=IUNIT,NAME=FILE,IOSTAT=IOSTAT,EXIST=EXIST,
C     '  OPENED=OPENED,FORM=FORMAT,ACCESS=ACCESS,NEXTREC=IREC)
C      IF(DOP) THEN
C        WRITE(OP_STRING,'('' UNIT='',I4,'' FILE='',A)') IUNIT,FILE
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' EXIST='',L1,'' OPENED='',L1)') EXIST,OPENED
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' ACCESS='',A,'' FORMAT='',A)') ACCESS,
C     '    FORMAT
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' NEXTREC='',I5)') IREC
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C      ENDIF
      INQUIRE(UNIT=IUNIT,EXIST=EXIST,OPENED=OPENED)
      IF(DOP) THEN
        WRITE(OP_STRING,'('' UNIT='',I4)') IUNIT
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' EXIST='',L1,'' OPENED='',L1)') EXIST,OPENED
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      IF(EXIST.AND.OPENED) THEN
        INQUIRE(UNIT=IUNIT,NAME=FILE,IOSTAT=IOSTAT,FORM=FORMAT,
     '    ACCESS=ACCESS,NEXTREC=IREC)
        IF(DOP) THEN
          WRITE(OP_STRING,'('' FILE='',A)') FILE
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' ACCESS='',A,'' FORMAT='',A)') ACCESS,
     '      FORMAT
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NEXTREC='',I5)') IREC
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
C CPB 21/4/94 Moving all files to sequential access stream_lf files
C If a file has a direct access it needs to be written to a
C old start
C        CLOSE(UNIT=IUNIT,IOSTAT=IOSTAT)
C        IF(IOSTAT.EQ.0) THEN
C          CALL STRING_TRIM(FILE,IBEG,IEND)
C          ICLEN=IEND-IBEG+1
C          CALL FSKF2C(FILE,ICLEN,INTSTR)
C          CALL CONVERT_FILE(INTSTR)  !in file "convert_file.c"
C        ELSE
C          CIUNIT=CFROMI(IUNIT,'(I10)')
C          CALL STRING_TRIM(CIUNIT,IBEG,IEND)
C          ERROR=' File close error. Unit = '//CIUNIT(IBEG:IEND)
C          GOTO 9999
C        ENDIF
C      ENDIF
C old end

        IF(ACCESS(1:6).EQ.'DIRECT') THEN
C CPB 21/4/94 Copy the direct access scratch file to a sequential
C access stream_lf file.
          IF(FSTATUS(IUNIT).EQ.'N') THEN ! file needs to be copied i.e. was a new file
            CALL STRING_TRIM(FILES(IUNIT),IBEG,IEND)
            FILENAME=FILES(IUNIT)(IBEG:IEND)
            IF(DOP) THEN
              WRITE(OP_STRING,'('' Copying scratch file to : '',A)')
     '          FILENAME
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            OPEN(UNIT=99,FILE=FILENAME,STATUS='unknown',
     '        ACCESS='SEQUENTIAL',FORM='FORMATTED',IOSTAT=IOSTAT)
            IF(IOSTAT.LE.0) THEN
              DO I=1,IREC-1
                READ(UNIT=IUNIT,REC=I,FMT='(A)',IOSTAT=IOSTAT) LINE
                IF(IOSTAT.LE.0) THEN
                  CALL STRING_TRIM(LINE,IBEG,IEND)
                  WRITE(UNIT=99,FMT=*,IOSTAT=IOSTAT) LINE(2:IEND)
                  IF(DOP) THEN
                    WRITE(OP_STRING,'('' Input line> '',A)')LINE(1:IEND)
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                  IF(IOSTAT.NE.0) GOTO 9998
                ELSE
                  GOTO 9998
                ENDIF
              ENDDO
            ELSE
              GOTO 9998
            ENDIF
            CLOSE(UNIT=99,IOSTAT=IOSTAT)
            IF(IOSTAT.GT.0) GOTO 9998

          ELSE IF((FSTATUS(IUNIT).EQ.'O').AND.MODIFYFILE(IUNIT)) THEN
            CALL STRING_TRIM(FILES(IUNIT),IBEG,IEND)
            FILENAME=FILES(IUNIT)(IBEG:IEND)
            IF(DOP) THEN
              WRITE(OP_STRING,'('' Copying scratch file to : '',A)')
     '          FILENAME
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            OPEN(UNIT=99,FILE=FILENAME,STATUS='OLD',
     '        ACCESS='SEQUENTIAL',FORM='FORMATTED',IOSTAT=IOSTAT)
            IF(IOSTAT.LE.0) THEN
              DO i=1,IREC-1
                READ(UNIT=IUNIT,REC=I,FMT='(A)',IOSTAT=IOSTAT) LINE
                IF(IOSTAT.LE.0) THEN
                  CALL STRING_TRIM(LINE,IBEG,IEND)
                  WRITE(UNIT=99,FMT='(A)',IOSTAT=IOSTAT) LINE(1:IEND)
                  IF(iostat.ne.0) GOTO 9998
                ELSE
                  GOTO 9998
                ENDIF !iostat
              ENDDO !i
            ELSE
              GOTO 9998
            ENDIF !iostat
            CLOSE(UNIT=99,IOSTAT=IOSTAT)
            IF(iostat.GT.0) GOTO 9998
          ENDIF !iold
C         We will get a close error when the fortran library tries to
C         unlink the scratch file (which we have already unlinked).
C         Therefore ignore close errors.
          CLOSE(UNIT=IUNIT,IOSTAT=IOSTAT)
        ELSE !not direct access
          CLOSE(UNIT=IUNIT,IOSTAT=IOSTAT)
          IF(IOSTAT.NE.0) GOTO 9998
        ENDIF
      ENDIF

      ERROR=' '
      CALL EXITS('CLOSEF')
      RETURN
C 9998 ISTAT=CFROMI(IOSTAT,'(I3)')
 9998 WRITE(ISTAT,'(I3)') IOSTAT
      ERROR='>>Iostat = '//ISTAT(1:3)
 9999 CALL ERRORS('CLOSEF',ERROR)
      CALL EXITS('CLOSEF')
      RETURN 1
      END


C unused!
C      SUBROUTINE FIND_FILE(NOFILE_START,NTFILE,FILE_EXT,FILE_LIST,
C     '  ERROR,*)
C
CC**** Finds files in current directory.
C
C      IMPLICIT NONE
C!     Parameter List
C      INTEGER NOFILE_START,NTFILE
C      CHARACTER ERROR*(*),FILE_EXT*(*),FILE_LIST(*)*(*)
C!     Local Variables
C      INTEGER CLOCAT,CONTEXT,IBEG,IBRA,IDOT,IEND,INDEX,LIB$FIND_FILE,
C     '  LIB$FIND_FILE_END,NOFILE
C      CHARACTER RESULT*100
C
C      CALL ENTERS('FIND_FILE',*9999)
C
C      CALL EXITS('FIND_FILE')
C      RETURN
C 9999 CALL ERRORS('FIND_FILE',ERROR)
C      CALL EXITS('FIND_FILE')
C      RETURN 1
C      END


