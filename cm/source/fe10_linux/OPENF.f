      SUBROUTINE OPENF(IUNIT,DEVICE,FILE,STATUS,ACCESS,FORM,IRECL,
     '  ERROR,*)

C#### Subroutine: OPENF
C###  Description:
C###    OPENF opens a file.

C**** IUNIT is the name by which the file is referred to in the program.
C**** Valid unit identifiers are non negative integers.
C**** DEVICE is the device on which the file resides.
C**** Valid devices are: 'DISK' 'TERM'
C**** FILE is the filename [filetype [filemode]]
C**** STATUS is 'NEW' or 'OLD' or 'SCRATCH'.  'NEW' is equivalent to
C**** Fortran 'UNKNOWN'.  i.e. it overwrites existing files.
C**** ACCESS specifies if the file is 'DIRECT' or 'SEQUEN' or ' APPEND'.
C**** Direct access files have a length of 2000 records.
C**** FORM specifies whether the file is 'FORMATTED' or 'UNFORMATTED'
C**** Direct access files have a length of 2000 records.
C**** IRECL is the logical record length.
C**** ERROR gives diagnostics in the event of failure.
C**** If an error is detected control is returned to the statement
C**** number of the star.
C**** DEVICE, FILE, ACCESS and ERROR are all character strings.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'disp00.cmn'
      INCLUDE 'docd00.cmn'
      INCLUDE 'file00.cmn'
!     Parameter List
      INTEGER IRECL,IUNIT
      CHARACTER ACCESS*(*),DEVICE*(*),ERROR*(*),FILE*(*),FORM*(*),
     '  STATUS*(*)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IEND,IEND1,
     '  IOSTAT,IREC
      CHARACTER LINE*256,
     &  FILENAME*(MXCH),
       ! max path length on AIX 5.1 and IRIX 6.5 is 1023; +1 for null
     &  SCRATCHNAME*(1024),
     &  STATUS_UNIX*10,UNIT*100
!     Functions
      INTEGER LEN_TRIM

      CALL ENTERS('OPENF',*9999)
      IF(DOP) THEN
        WRITE(OP_STRING(1),'(''  IUNIT = '',I2)') IUNIT
        WRITE(OP_STRING(2),'('' DEVICE = '', A)') DEVICE
        WRITE(OP_STRING(3),'(''   FILE = '', A)') FILE
        WRITE(OP_STRING(4),'('' STATUS = '', A)') STATUS
        WRITE(OP_STRING(5),'('' ACCESS = '', A)') ACCESS
        WRITE(OP_STRING(6),'(''   FORM = '', A)') FORM
        WRITE(OP_STRING(7),'(''  IRECL = '',I4)') IRECL
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      MODIFYFILE(IUNIT)=.FALSE.

      IF(DOCDIR) THEN !read from documentation directory

C!!! cpb 2/10/98 Set example path with a set directory command
C
C        IF(Open_COM_file) THEN
C          Backend_COM_file=.TRUE.
C
CC Locate the root location for the CMISS examples
C          CALL STRING_TRIM(EXAMPLE_PATH,IBEG1,IEND1)
C          IF(DOP) THEN
C            WRITE(OP_STRING,'('' example_path='',A)')
C     '        EXAMPLE_PATH(IBEG1:IEND1)
C            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C          ENDIF
C
C          IF(FILE(1:8).EQ."example_") THEN
CC Find the example subdirectory from the com file name
C            IBEG2=8 !e.g. example_31   ibeg2=8
C            IEND2=CLOCAT('.',FILE)-1 !iend2=10
Cc           write(*,'('' ibeg2='',I4,''iend2='',I4)')   ibeg2,iend2
C            SUB_DIR=FILE(IBEG2+1:IBEG2+1) !sub_dir='3'
C            i_tot=IEND2-IBEG2 !i_tot=2
C            DO i=2,i_tot !do i=2,2
C              CALL STRING_TRIM(SUB_DIR,IBEG3,IEND3)
Cc             write(*,'('' ibeg3='',I4,''iend3='',I4)') ibeg3,iend3
C              SUB_DIR=SUB_DIR(1:IEND3)//'/'//FILE(IBEG2+1:IBEG2+i)
C              if(dop) write(*,'('' sub_dir='',A)') sub_dir(1:60)
C            ENDDO !sub_dir='3/31'
C
CC Append the example subdirectory to the cmiss$examples logical
C           if(dop) write(*,'('' ibeg1='',I4,''iend1='',I4)') ibeg1,iend1
C           CALL STRING_TRIM(SUB_DIR,IBEG3,IEND3)
C           if(dop) write(*,'('' ibeg3='',I4,''iend3='',I4)') ibeg3,iend3
C           EXAMPLES_DIR=EXAMPLE_PATH(IBEG1:IEND1)//'/'
C     '       //SUB_DIR(IBEG3:IEND3)//'/'
C           CALL STRING_TRIM(EXAMPLES_DIR,IBEG1,IEND1)
C           if(dop) then
C             write(op_string,'('' examples_dir='',a)')
C     '         examples_dir(ibeg1:iend1)
C             call writes(iodi,op_string,error,*9999)
C           endif
C
CC Retain the example directory string for later use
C           Example_directory=EXAMPLES_DIR(IBEG1:IEND1)
C           if(dop) then
C             write(op_string,'('' example_directory='',a)')
C     '         example_directory(ibeg1:iend1)
C             call writes(iodi,op_string,error,*9999)
C           endif
C         ENDIF
C
C        ELSE !define commands from an already opened COM file
C          IF(Backend_COM_file) THEN
C            EXAMPLES_DIR=Example_directory !defined when OPEN file
C          ELSE !an examples subdirectory has been set
C            CALL STRING_TRIM(EXAMPLE_PATH,IBEG1,IEND1)
CC Append the example subdirectory to the cmiss$examples logical
C            CALL STRING_TRIM(Example_subdir,IBEG2,IEND2)
C            IF(IBEG2.EQ.1.AND.IEND2.EQ.1) THEN !zero length subdir str
C              EXAMPLES_DIR=EXAMPLE_PATH(IBEG1:IEND1)//'/'
C            ELSE
C              IF(Example_subdir(IEND2:IEND2).EQ.'/') IEND2=IEND2-1
C              EXAMPLES_DIR=EXAMPLE_PATH(IBEG1:IEND1)
C     '          //Example_subdir(IBEG2:IEND2)//'/'
C            ENDIF
C          ENDIF
C
C        ENDIF !Open_COM_file

C Append filename to examples pathname

        CALL STRING_TRIM(EXAMPLES_DIR,IBEG1,IEND1)
C KAT 27/4/00: I_end_of_path not used
C        I_end_of_path=IEND1-IBEG1+1
        CALL STRING_TRIM(FILE,IBEG,IEND)
        FILENAME=EXAMPLES_DIR(IBEG1:IEND1)//FILE(IBEG:IEND)
        IF(DOP) THEN
          WRITE(OP_STRING,'('' EXAMPLES_DIR='',A,/,'' FILENAME='',A)')
     '      EXAMPLES_DIR(IBEG1:IEND1),
     '      FILENAME(1:IEND1-IBEG1+IEND-IBEG+2)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

        DOCDIR=.FALSE.

      ELSE IF(.NOT.DOCDIR) THEN
C KAT 27/4/00: I_end_of_path not used
CC       Locate position in FILE of end of directory path string
C        I=0
C        I_end_of_path=0
C        LAST=.FALSE.
C        DO WHILE (.NOT.LAST)
C          I=CLOCAT('/',FILE(I_end_of_path+1:))
C          if(dop) write(*,'('' i='',I2)') i
C          IF(I.EQ.0) THEN
C            LAST=.TRUE.
C          ELSE
C            I_end_of_path=I_end_of_path+I
C            if(dop) write(*,'('' I_end_of_path='',I2)') I_end_of_path
C          ENDIF
C        ENDDO !while .not.last

CC KAT 1/5/00: prepending current directory if not absolute path
C        IF(FILE(:1).EQ.'/') THEN
C          FILENAME=FILE
C        ELSE
C          FILENAME=PATH00(:LEN_TRIM(PATH00))//FILE
C        ENDIF
        FILENAME=FILE
        if(dop) then
          write(op_string(1),'(''   filename = '', a)') filename
          call writes(iodi,op_string,error,*9999)
        endif

      ENDIF !docdir

      CALL STRING_TRIM(FILENAME,IBEG,IEND)

      IF(DEVICE.EQ.'TERM') THEN
        OPEN(UNIT=IUNIT,FILE=FILENAME(IBEG:IEND),STATUS=STATUS,
     '    RECL=IRECL,IOSTAT=IOSTAT)

      ELSE IF(DEVICE.EQ.'DISK') THEN
        IF(ACCESS.EQ.'SEQUEN') THEN
C          STATUS_UNIX=CUPPER(STATUS)
          CALL CUPPER(STATUS,STATUS_UNIX)
          IF(STATUS_UNIX(1:3).EQ.'NEW') STATUS_UNIX='UNKNOWN'
          OPEN(UNIT=IUNIT,FILE=FILENAME(IBEG:IEND),STATUS=STATUS_UNIX,
     '      ACCESS='SEQUENTIAL',FORM=FORM,IOSTAT=IOSTAT)

        ELSE IF(ACCESS.EQ.'APPEND') THEN
C***      Fortran 90 has ACCESS='SEQUENTIAL',POSITION='APPEND'
          OPEN(UNIT=IUNIT,FILE=FILENAME(IBEG:IEND),STATUS=STATUS,
     '      ACCESS='APPEND',FORM=FORM,IOSTAT=IOSTAT)

        ELSE IF(ACCESS.EQ.'DIRECT') THEN
          IF(STATUS.EQ.'OLD') THEN
            OPEN(UNIT=99,FILE=FILENAME(IBEG:IEND),STATUS=STATUS,
     '        ACCESS='SEQUENTIAL',FORM=FORM,IOSTAT=IOSTAT)

            IF(IOSTAT.EQ.0) THEN
              OPEN(UNIT=IUNIT,STATUS='SCRATCH',
     '          ACCESS='DIRECT',FORM=FORM,IOSTAT=IOSTAT,
     '          RECL=IRECL)
              IREC=0
C KAT 22Dec99: Checking that file is created
              IF(IOSTAT.NE.0) THEN
                CALL FLAG_ERROR(IOSTAT,'Cannot create scratch file')
                GOTO 9998
              ENDIF
C!!! Intel fortran version 8.0 and probably 9.1 will start seg faulting
C!!! when opening sequential files if we unlink the scratch file then
C!!! close it.

C             Unlink the file so that it does not hang around if we exit
C             abnormally.
              INQUIRE(UNIT=IUNIT,NAME=SCRATCHNAME,IOSTAT=IOSTAT)
              IF(IOSTAT.NE.0) THEN
                CALL FLAG_ERROR(IOSTAT,'Inquire on scratch file failed')
                GOTO 9998
              ENDIF
              IEND=LEN_TRIM(SCRATCHNAME)
              IF(IEND.EQ.0) THEN
                CALL WRITE_CHAR(IOER,
     &            'WARNING: Can''t unlink scratch file: name is blank'//
     &            NEWLINE,ERR)
              ELSEIF(IEND.EQ.LEN(SCRATCHNAME)) THEN
C               The name didn't fit into the buffer.
C               Don't accidentally delete another file with the
C               truncated name.
                CALL WRITE_CHAR(IOER,
     &            'WARNING: Can''t unlink scratch file: '//
     &            'name too long for buffer'//NEWLINE,ERR)
              ELSE
                IEND=IEND+1
                SCRATCHNAME(IEND:IEND)=CHAR(0)
                CALL UNLINK_FILE(SCRATCHNAME,ERR)
                ! Any error message is already output
              ENDIF
              DO WHILE(IOSTAT.EQ.0)
                READ(99,'(A)',IOSTAT=IOSTAT) LINE
                IF(IOSTAT.EQ.0) THEN
                  IREC=IREC+1
                  WRITE(UNIT=IUNIT,FMT='(A)',REC=IREC,IOSTAT=IOSTAT)
     '              LINE(1:IRECL)
                  IF(IOSTAT.NE.0) THEN
                    CALL FLAG_ERROR(0,' ')
                    CALL WRITE_CHAR(IOER,
     '                'can''t write to scratch file: IOSTAT=',ERR)
                    CALL WRITE_INT(IOER,IOSTAT,ERR)
                    CALL WRITE_CHAR(IOER,NEWLINE,ERR)
                  ENDIF
                ELSEIF(IOSTAT.GT.0) THEN ! file error
                  CALL FLAG_ERROR(0,' ')
                  CALL WRITE_CHAR(IOER,'can''t read file ',ERR)
                  CALL WRITE_STRING(IOER,LEN_TRIM(FILENAME),FILENAME,
     '              ERR)
                  CALL WRITE_CHAR(IOER,': IOSTAT=',ERR)
                  CALL WRITE_INT(IOER,IOSTAT,ERR)
                  CALL WRITE_CHAR(IOER,NEWLINE,ERR)
                ENDIF
              ENDDO !IOSTAT.EQ.0
              CLOSE(UNIT=99)
C 10           READ(99,'(A)',IOSTAT=IOSTAT,END=20) LINE
CC cpb 20/4/94 This inquire doesn't seem to pick up the next record
CC properly
CC             INQUIRE(UNIT=IUNIT,NEXTREC=IREC)
C              IF(IOSTAT.EQ.0) THEN
C                IREC=IREC+1
C                WRITE(UNIT=IUNIT,FMT='(A)',REC=IREC,IOSTAT=IOSTAT)
C     '            LINE(1:IRECL)
C                GOTO 10
C              ELSE ! file error
C                GOTO 20
C              ENDIF
C 20           CLOSE(UNIT=99)
              IF(IOSTAT.LT.0) THEN ! EOF
                IOSTAT=0
              ELSEIF(IOSTAT.GT.0) THEN ! file error
                CLOSE(UNIT=IUNIT,IOSTAT=ERR) ! ignore failure to unlink
                GO TO 9999
              ENDIF
C KAT 21Jan00: There is no toslf command
C            ELSE IF(IOSTAT.EQ.44) THEN ! Temporary - need to inform
C               ! user to convert the file.
C              WRITE(OP_STRING,'('' >>Convert DIRECT access file to '',
C     '          ''a SEQUENTIAL access STREAM_LF file using the TOSLF '',
C     '          ''command'')')
C              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF

          ELSE IF(STATUS.EQ.'NEW') THEN
            OPEN(UNIT=IUNIT,STATUS='SCRATCH',ACCESS='DIRECT',
     '        FORM=FORM,RECL=IRECL,IOSTAT=IOSTAT)
          ENDIF

        ELSE
          ERROR=' ACCESS='//ACCESS//' is invalid'
          GOTO 9999
        ENDIF

      ELSE
        ERROR=' DEVICE='//DEVICE//' is invalid'
        GOTO 9999
      ENDIF

      IF(IOSTAT.NE.0) THEN
C        ERROR=CFROMI(IOSTAT,'(I3)')
        WRITE(ERROR,'(I3)') IOSTAT
C CPB 14/6/94 Not sure about these iostats for the sgi
C        UNIT=CFROMI(IUNIT,'(I4)')
        WRITE(UNIT,'(I4)') IUNIT
        ERROR=' Iostat='//ERROR(1:3)//' error in OPENF(UNIT='
     '    //UNIT(1:4)//')'
Cnew    for sgi the iostat for file not found is 2 or 157?AAY 21 Feb 92
C       for linux port f90 file not found is 209
        IF(IOSTAT.EQ.29.OR.IOSTAT.EQ.2.OR.IOSTAT.EQ.157
     '       .OR.IOSTAT.EQ.209) THEN
          CALL STRING_TRIM(ERROR,IBEG,IEND)
          ERROR=' File not found: '//ERROR(IBEG:IEND)
        ELSE IF(IOSTAT.EQ.13) THEN
          CALL STRING_TRIM(ERROR,IBEG,IEND)
          ERROR=' Permission denied: '//ERROR(IBEG:IEND)
        ELSE IF(IOSTAT.EQ.43) THEN
          CALL STRING_TRIM(ERROR,IBEG,IEND)
          ERROR=' Filename invalid: '//ERROR(IBEG:IEND)
        ELSE IF(IOSTAT.EQ.44) THEN
          CALL STRING_TRIM(ERROR,IBEG,IEND)
          ERROR=' Inconsistent record type: '//ERROR(IBEG:IEND)
        ENDIF
        WRITE(OP_STRING,'('' Filename: '',A)') FILENAME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        GOTO 9999
      ENDIF

      FILES(IUNIT)=FILENAME
      FSTATUS(IUNIT)=STATUS(1:1)
      ERROR='  '
      CALL EXITS('OPENF')
      RETURN
 9998 ERROR=' '
 9999 CALL ERRORS('OPENF',ERROR)
      CALL EXITS('OPENF')
      RETURN 1
      END


C      SUBROUTINE POST_FILE(FILE_NAME)
C
CC**** Posts file to printer.
C
C      IMPLICIT NONE
C      INCLUDE 'cmiss$reference:cbdi02.cmn'
C!     Parameter List
C      CHARACTER FILE_NAME*(*)
C!     Local Variables
C      INTEGER IBEG,IEND,ISTATUS,LIB$SPAWN
C
C      CALL STRING_TRIM(FILE_NAME,IBEG,IEND)
C
C      RETURN
C      END
C
C
C      SUBROUTINE PURGE_FILE(FILE_NAME)
C
CC**** Purges file versions.
C
C      IMPLICIT NONE
C      INCLUDE 'cmiss$reference:cbdi02.cmn'
C!     Parameter List
C      CHARACTER FILE_NAME*(*)
C!     Local Variables
C      INTEGER IBEG,IEND,ISTATUS,LIB$SPAWN
C
C      CALL STRING_TRIM(FILE_NAME,IBEG,IEND)
C
C      RETURN
C      END


