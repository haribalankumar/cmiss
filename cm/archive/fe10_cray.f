C#### Module: FE10_CRAY
C###  Description:
C###    System specific routines (CRAY)

C!!!! Note: If necessary replace directory string
C!!!!       "/product/cmiss/vms/document/examples/"

C###  Routine: NEW_HANDLER     (fn) condition handling function
C###  Routine: COPY_FILE     copies file between directories
C###  Routine: CLOSEF        closes file
C###  Routine: DIALOG        controls command parsing environment
C###  Routine: DO_COMMAND    executes an operating system command
C###  Routine: ENTERS        traces entry to a subprogram
C###  Routine: EXECUTE_COMFILE Executes an operating system com file
C###  Routine: FIND_FILE     finds files in current directory
C###  Routine: GET_COMMAND_LINE  gets the command line arguments
C###  Routine: GET_DATE_TIME uses rtl call to return date and time
C###  Routine: GET_DISPLAY   uses rtl call to return display name
C###  Routine: GET_REV_TIME  return CMISS revision time
C###  Routine: GET_SYSTEM    uses rtl call to return system parameters
C###  Routine: GETSTR1       rets str typed by user from keybd mapping
C###  Routine: GETSTR2       returns interrupt key number
C###  Routine: GETSTR3       returns key for DOCUM
C###  Routine: OPENF         opens file
C###  Routine: POST_FILE     posts file to printer
C###  Routine: PURGE_FILE    purges file versions
C###  Routine: SETCMISSCOMMANDPARAMS sets command parameters
C###  Routine: SETUPCMISS   setup initial CMISS constants
C###  Routine: SLEEPER       sleeps for an integer period of seconds


      SUBROUTINE COPY_FILE(FILE_NAME)

C**** Copies file between directories.

      IMPLICIT NONE
!     Parameter List
      CHARACTER FILE_NAME*(*)
!     Local Variables
      INTEGER IBEG,IEND,ISTATUS


      RETURN
      END


      SUBROUTINE CLOSEF(IUNIT,ERROR,*)

C**** Closes a file using the FORTRAN CLOSE command.
C**** IUNIT is the name by which the file is referred to in the program.
C**** Direct access files are truncated at their current record.
C**** ERROR gives diagnostics in the event of failure.
C**** If an error is detected control is returned to the statement
C**** number of the star.
C**** ERROR is a character string.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:mxch.inc'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
!     Parameter List
      INTEGER IUNIT
      CHARACTER ERROR*(*)
!     Local Variables
c      INTEGER I,IBEG,ICLEN,IEND,INTSTR(1024),IOSTAT,IREC
      INTEGER I,IBEG,IEND,IOSTAT,IREC
      CHARACTER ACCESS*10,FILE*100,FILENAME*100,
     '  FORMAT*11,ISTAT*5,LINE*256
      LOGICAL EXIST,OPENED

      CALL ENTERS('CLOSEF',*9999)
      INQUIRE(UNIT=IUNIT,NAME=FILE,IOSTAT=IOSTAT,EXIST=EXIST,
     '  OPENED=OPENED,FORM=FORMAT,ACCESS=ACCESS,NEXTREC=IREC)
      IF(DOP) THEN
        WRITE(OP_STRING,'('' UNIT='',I4,'' FILE='',A)') IUNIT,FILE
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' EXIST='',L1,'' OPENED='',L1)') EXIST,OPENED
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' ACCESS='',A,'' FORMAT='',A)') ACCESS,
     '    FORMAT
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' NEXTREC='',I5)') IREC
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      IF(EXIST.AND.OPENED) THEN

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
     '        ACCESS='SEQUENTIAL',FORM='FORMATTED',IOSTAT=IOSTAT,
     '        RECORDTYPE='STREAM_LF')
            IF(IOSTAT.EQ.0) THEN
              DO I=1,IREC-1
                READ(UNIT=IUNIT,REC=I,FMT='(A)',IOSTAT=IOSTAT) LINE
                IF(IOSTAT.EQ.0) THEN
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
            IF(IOSTAT.NE.0) GOTO 9998
          ENDIF
        ENDIF
        CLOSE(UNIT=IUNIT,IOSTAT=IOSTAT)
        IF(IOSTAT.NE.0) GOTO 9998
      ENDIF

      ERROR=' '
      CALL EXITS('CLOSEF')
      RETURN
 9998 WRITE(ISTAT,'(I3)') IOSTAT
      ERROR='>>Iostat = '//ISTAT(1:3)
 9999 CALL ERRORS('CLOSEF',ERROR)
      CALL EXITS('CLOSEF')
      RETURN 1
      END


      SUBROUTINE DIALOG

C**** The dialog interpreter

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:fsklib.inc'
      INCLUDE 'cmiss$reference:mxch.inc'
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:binf00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:cbfe01.cmn'
      INCLUDE 'cmiss$reference:cbpr00.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:cmis00.cmn'
C      INCLUDE 'cmiss$reference:cspi00.cmn'
C      INCLUDE 'cmiss$reference:ctrl00.cmn'
      INCLUDE 'cmiss$reference:diag00.cmn'
      INCLUDE 'cmiss$reference:dial00.cmn'
      INCLUDE 'cmiss$reference:disp00.cmn'
      INCLUDE 'cmiss$reference:docu00.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:gen000.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:gks000.cmn'
      INCLUDE 'cmiss$reference:gks001.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'cmiss$reference:gtstr00.cmn'
      INCLUDE 'cmiss$reference:head00.cmn'
      INCLUDE 'cmiss$reference:host00.cmn'
      INCLUDE 'cmiss$reference:host00.inc'
      INCLUDE 'cmiss$reference:ktyp00.cmn'
      INCLUDE 'cmiss$reference:mach00.cmn'
      INCLUDE 'cmiss$reference:mach00.inc'
      INCLUDE 'cmiss$reference:map000.cmn'
      INCLUDE 'cmiss$reference:ntsg00.cmn'
      INCLUDE 'cmiss$reference:trac00.cmn'
!     Local Variables
      INTEGER MXCO,MXCOQU,MXSG                !MXCH in 'mxch.inc'
      PARAMETER (MXCO=25,MXCOQU=25,MXSG=2000)
      INTEGER I,IBEG,ICODE,ID,ID_DEVICE,ID_STATUS,ID_WS,IEND,IFROMC,
     '  INDEX_PLIN,INLIST,INPUT_CHAN,INPUT_CHOICE,INTSTR(MXCH),IRET,
     '  ISEG(MXSG),ISEGM,ISEGM_LIST(0:100),ISTATUS,IW,J,LEN,
     '  MODE_PROJ,N,N1LIST,NCENTRE,NCHAR,NEW_HANDLER,NJ,nhost,NOCH,NOCO,
     '  NO_COL,NO_GROUP_TRANSFORM,NOGRPL,NOMACRO,NOPTS,NOSG,NSUB,NTCH,
     '  NTCOD,NTCOQUD(8),NTFILE,NT_VIEW,ERR,RETURNVAL
      REAL DX,DY,FANGLE(3),FPT(3),FSCALE(3),FSHFT(3),
     '  SDATA,XNDC1,XNDC2,XNDC3,XNDC4,VALUE,VMATRIX(4,4),X0,XCENTRE,
     '  Y0,YCENTRE
      LOGICAL CHANGE,CONTINUE,DEFINE,END,FILEIP,FIRST_ZOOM,GRID,
     '  QUALIFY,REFINE_ACTIVE,TRANSFORM_ACTIVE,UPDATE,UPDATE_ZOOM,
     '  DISPLAY_VALUATOR_81,DISPLAY_VALUATOR_82,DISPLAY_VALUATOR_83,
     '  DISPLAY_VALUATOR_84,DISPLAY_VALUATOR_85,DISPLAY_VALUATOR_86,
     '  DISPLAY_VALUATOR_87,DISPLAY_VALUATOR_88,DISPLAY_VALUATOR_89
      CHARACTER CHAR1*1,CHAR3*3,CHAR5*5,CHOOSE*40,
     '  CIW*1,CLASS*8,CO(MXCO)*30,COD(MXCO)*90,COQUD(MXCO,MXCOQU)*30,
     '  CSEG(MXSG)*60,ERROR*(MXCH),FILE_NAME*50,OPTION(30)*40,
     '  OPTION1(18)*11,OPTION2(11)*11,OPTION3(25)*11,OPTION4(10)*11,
     '  OPTION5(9)*20,OPTION10(10)*20,OPTION11(10)*20,
     '  PARAMETER_TYPE*20,Q*1,STATSTR*80,STRG*(MXCH),STRING*(MXCH),
     '  TEXT_STRING*80,TRANSFORM_TYPE*7

      DATA FPT/0.0,0.0,0.0/,FSHFT/0.0,0.0,0.0/,FSCALE/1.0,1.0,1.0/,
     '  FANGLE/0.0,0.0,0.0/

      DO I=1,100
        DO J=1,MXCH
          OP_STRING(I)(J:J)=' '
        ENDDO
      ENDDO
      GKS=.FALSE.
C      MAPOPN=.FALSE.
      SMG_READ=.TRUE.
c cpb 22/10/95 Set machine constants (don't change)
      CHARSIZE=1
      INTSIZE=8
      SINTSIZE=2
      SPSIZE=4
      DPSIZE=8
      LOGSIZE=4
      ENDIANTYPE=CHAR(MACH_UNKNOWN) !?
      SPFORMTYPE=CHAR(MACH_UNKNOWN) !?
      DPFORMTYPE=CHAR(MACH_UNKNOWN) !?
      MACHTYPE=CHAR(MACH_CRAY) !CRAY
      OSTYPE=CHAR(MACH_UNKNOWN) !?
c
      TR01=.FALSE.
      TR02=.FALSE.
      TR03=.FALSE.
      TR04=.FALSE.
      FIRSTF=.TRUE.
      FIRST_ZOOM =.TRUE.
      UPDATE_ZOOM=.FALSE.
      BUFFER_COUNT=0
      IREC_COMFILE=0
      CALL ENTERS('DIALOG',*9999)
      CALL GET_REV_TIME(ERROR,*9999)
      CALL GET_SYSTEM(ERROR,*9999)
      CALL GET_DISPLAY(ERROR,*9999)
      EXAMPLES_DIR='/product/cmiss/vms/document/examples/'
      IOIP=1 !input  file
      IOOP=2 !usual output
      IODI=3 !diagnostics output
      IOTR=4 !trace output
      IOER=5 !error output
      IOH1=6 !help 1 output
      IOH2=7 !help 2 output
      IOOUT=9 !output FILE
      CALL OPENF(1,'TERM','SYS$INPUT' ,'UNKNOWN',' ',' ',132,
     '  ERROR,*9999)
      CALL OPENF(2,'TERM','SYS$OUTPUT','UNKNOWN',' ',' ',255,
     '  ERROR,*9999)
      CALL OPENF(3,'TERM','SYS$OUTPUT','UNKNOWN',' ',' ',255,
     '  ERROR,*9999)
      CALL OPENF(4,'TERM','SYS$OUTPUT','UNKNOWN',' ',' ',255,
     '  ERROR,*9999)
      CALL OPENF(5,'TERM','SYS$OUTPUT','UNKNOWN',' ',' ',255,
     '  ERROR,*9999)
      CALL OPENF(6,'TERM','SYS$OUTPUT','UNKNOWN',' ',' ',255,
     '  ERROR,*9999)
      CALL OPENF(7,'TERM','SYS$OUTPUT','UNKNOWN',' ',' ',255,
     '  ERROR,*9999)
      CALL OPENF(8,'TERM','SYS$OUTPUT','UNKNOWN',' ',' ',255,
     '  ERROR,*9999)
      IOFI=IOOP !default for listing to file commands
      IFILE=10 !ifile is the file for input files i.e. ip* files
      IOFILE1=11 !output file i.e. op* files
      IOFILE2=12 ! iofile2 -> iofile6 are general purpose output files
      IOFILE3=13
      IOFILE4=14
      IOFILE5=17
      IOFILE6=18
C CPB 30/11/92 the IO# units below should (?) be redundant.
      IO1=1
      IO2=2
      IO3=IOIP
      IO4=IOOP
      IVDU=IOIP
      IMAP=0
      DO NUM_STACK=1,5                !<-|
        DOP_STACK(NUM_STACK)=.FALSE.  !  |
      ENDDO                           !  |
      NUM_STACK=5 !5 levels below FEM !  |- used for diagnostics
      NT_SUB=0                        !  |
      DO NSUB=1,MXSUB                 !  |
        SUBNAM(NSUB)=' '              !  |
      ENDDO                           !  |
      DIAGNO=.FALSE.                  !  |
      ALLSUB=.FALSE.                  !  |
      FROMSUB=.FALSE.                 !  |
      DOP=.FALSE.                     !<-|
      NTSG=0
      NT_KEY=0
      NT_MACRO(0)=0
      MACRO_DEFINE=.FALSE.
      MACRO_COMMAND_EXECUTE=.FALSE.
      MACRO_KEY_EXECUTE=.FALSE.
      DISPLAY_VALUATOR_81=.FALSE.
      DISPLAY_VALUATOR_82=.FALSE.
      DISPLAY_VALUATOR_83=.FALSE.
      DISPLAY_VALUATOR_84=.FALSE.
      DISPLAY_VALUATOR_85=.FALSE.
      DISPLAY_VALUATOR_86=.FALSE.
      DISPLAY_VALUATOR_87=.FALSE.
      DISPLAY_VALUATOR_88=.FALSE.
      DISPLAY_VALUATOR_89=.FALSE.
      REFINE_ACTIVE=.FALSE.
      TRANSFORM_ACTIVE=.FALSE.
      DO_EXAMPLE=.FALSE.
      SELECT_EXAMPLE=.FALSE.
      NT_VIEW=0
      HEADING=' '
C      DO NO_COL=1,9
C        NAME_COL(NO_COL)='UNDEFINED'
C      ENDDO
      DO I=1,99
        IWKS(I)=0
      ENDDO
c     IWKDEF(0)=1
c     IWKDEF(1)=1
      DO NOSG=1,MXSG
        ISEG(NOSG)=0
        CSEG(NOSG)=' '
      ENDDO
      XMIN=-1.0
      XMAX= 1.0
      YMIN=-1.0
      YMAX= 1.0
      ZMIN=-1.0
      ZMAX= 1.0
      DIAG=SQRT(12.0)
      NJT=2
      FILE00='file'
      MENU=.FALSE.
      COD(2)='PI'
      WRITE(COD(3),'(E11.5)') PI
      CALL ASSIGN(STRING,1,COD,COQUD,ERROR,*9999)

      IF(USE_SOCKET.EQ..FALSE.) THEN

C ***   Check for cmiss.com file
        END=.FALSE.
C news MPN 6-Jul-95: changing "read file;com" to "read com;file" etc
        STRG='read com;cmiss'
C old        STRG='read cmiss;com'
        CALL PARSE(0,ISEG,CSEG,STRG,END,ERROR,*9999)

        IF(END) THEN !quit issued from cmiss.com
          CONTINUE=.FALSE.
        ELSE
          CONTINUE=.TRUE.
        ENDIF
        DO WHILE(CONTINUE) !is main program loop

C GMH 15/11/95 making call conditional
          IF(USE_GRAPHICS.EQ.1) THEN
            CALL GXWAIT(0.0,ERR) !Update graphics
          ENDIF

          CALL SETSTR(STRING,ERROR,*150)
          CALL GETSTR(STRING(LNPR-2:),PR(:LNPR),END,ERROR,*150)

          CALL PARSE(0,ISEG,CSEG,STRING,END,ERROR,*150)
          STRING=' ' !Test 4-3-1991

          IF(END) THEN       !The QUIT command has been used to quit
            CONTINUE=.FALSE.
          ENDIF
          GOTO 200

C ***     Handle error condition
 150      CALL STRING_TRIM(ERROR,IBEG,IEND)
          WRITE(OP_STRING,'(A)') ' '//ERROR(IBEG:IEND)//'>DIALOG'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ERROR(1:)=' '

 200      CONTINUE
        ENDDO !end of main program loop

      ELSE   !USE_SOCKET = .TRUE.

        CONNID1=0
        CONNID2=1

        IF (FSKLISTEN(CONNID1,PORT1) .EQ. -1) GOTO 9998
        IF (FSKLISTEN(CONNID2,PORT2) .EQ. -1) GOTO 9998

        CONTINUE=.TRUE.
        DO WHILE (CONTINUE)
C CPB 6/4/95 Adding timeout for sockets to enable graphics updates
          RETURNVAL=0
          DO WHILE(RETURNVAL.EQ.0)
            RETURNVAL=FSKSELECT(CONNID1,200) ! timeout after 200ms
            IF(RETURNVAL.EQ.-1) GOTO 9998
C GMH 15/11/95 making call conditional
            IF(USE_GRAPHICS.EQ.1) THEN
              CALL GXWAIT(0.0,ERR) !Update graphics
            ENDIF
          ENDDO
          IF (FSKREAD(LEN,SK_LONG_INT,1,CONNID1) .EQ. -1) GOTO 9998
          IF (FSKREAD(INTSTR,SK_CHAR,LEN+1,CONNID1) .EQ. -1) GOTO 9998
          CALL FSKC2F(STRING,LEN,INTSTR)

          CALL STRING_TRIM(STRING,IBEG,IEND)
          WRITE(OP_STRING,'(A)') STRING(IBEG:IEND)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9998)
          IF(STRING(IBEG:IBEG+3).NE.'QUIT') THEN
            CALL PARSE(0,ISEG,CSEG,STRING,END,ERROR,*160)
            WRITE(OP_STRING,'(/1X,A)') '> '
            CALL WRITES(IOOP,OP_STRING,ERROR,*9998)
          ENDIF

          GOTO 210

C ***   Handle error condition
 160      CALL STRING_TRIM(ERROR,IBEG,IEND)
          WRITE(OP_STRING,'(A)') ' '//ERROR(IBEG:IEND)//'>DIALOG'
          CALL WRITES(IOER,OP_STRING,ERROR,*9998)
          ERROR(1:)=' '

C ***   Tell front end that CMISS is back in the main loop (used for
C ***   popping down the prompt window)
 210      IF(FSKWRITE(-1,SK_LONG_INT,1,CONNID2).EQ.-1) GOTO 9998

          IF (STRING(1:4) .EQ. 'QUIT') THEN
C ***   First tell front-end to close down sockets, then close down own
C ***   sockets.
            IF (FSKWRITE(-2,SK_LONG_INT,1,CONNID2).EQ.-1) GOTO 9998
            IRET=FSKCLOSE(CONNID1)
            IRET=FSKCLOSE(CONNID2)
            CONTINUE=.FALSE.
          ENDIF

        ENDDO

      ENDIF

C ??? MPN 27-oct-94 Is there a better place this?
!news MPN 13-Sep-94: Parallel element stiffness matrix computations
      IF(KTYP1A.EQ.2) THEN  !parallel element stiffness matrix calcs
        DO nhost=1,NUMHOSTS_USED
          IF(SOCKET_OPEN(nhost)) THEN
            !Signal to slave processes to stop
            IF(FSKWRITE(QUIT_PROCESS,SK_LONG_INT,1,ICONNID(nhost))
     '        .EQ.-1) GOTO 9999
            IRET=FSKCLOSE(ICONNID(nhost)) !close socket
          ENDIF
        ENDDO
      ENDIF
!newe

      CALL EXITS('DIALOG')

      RETURN
 9999 CALL ERRORS('DIALOG',ERROR)
      WRITE(OP_STRING,'(132A)') ' Fatal error propagated to DIALOG:'
      CALL WRITES(IOER,OP_STRING,ERROR,*9998)
      CALL STRING_TRIM(ERROR,IBEG,IEND)
      WRITE(OP_STRING,'(132A)') ' '//ERROR(IBEG:IEND)//'>DIALOG'
      CALL WRITES(IOER,OP_STRING,ERROR,*9998)
      CLOSE(UNIT=IOER)
      CALL EXITS('DIALOG')
 9998 IF(USE_SOCKET) THEN
        IRET=FSKCLOSE(CONNID1)
        IRET=FSKCLOSE(CONNID2)
      ENDIF
      RETURN
      END


      SUBROUTINE DO_COMMAND(COMMAND,ERROR,*)

C#### Subroutine: DO_COMMAND
C###  Description:
C###    DO_COMMAND executes an operating system command.

      IMPLICIT NONE
!     Parameter List
      CHARACTER COMMAND*(*),ERROR*(*)
!     Local Variables
      INTEGER C_COMMAND(50),C_ERROR(50),CSTRLEN,ERR,IBEG,IEND

      CALL ENTERS('DO_COMMAND',*9999)

      CALL STRING_TRIM(COMMAND,IBEG,IEND)
      COMMAND=COMMAND(IBEG:IEND)//' &'
      CALL F2CSTRING(C_COMMAND,COMMAND)
      CALL SYSTEMCOMMAND(C_COMMAND,ERR,C_ERROR)
      IF(ERR.NE.0) THEN
        CALL CSTRINGLEN(CSTRLEN,C_ERROR)
        CALL C2FSTRING(C_ERROR,CSTRLEN,ERROR)
        GOTO 9999
      ENDIF

      CALL EXITS('DO_COMMAND')
      RETURN
 9999 CALL ERRORS('DO_COMMAND',ERROR)
      CALL EXITS('DO_COMMAND')
      RETURN 1
      END


      SUBROUTINE ENTERS(NAME,*)

C**** Traces entry to a subprogram recording the level of nesting,
C**** the entry time, and writing the subprogram name to a trace file.
C**** Diagnostic o/p is turned on if DIAGNO=.TRUE. and ALLSUB=.TRUE. or
C**** ALLSUB=.FALSE. & NAME=SUBNAM (a subroutine name).
C**** TRSB         is subprogram name which turns trace on
C**** IOTR         is trace file number (diagnostic o/p file)
C**** TR01         is true if basic tracing on
C**** TR02         is true if full tracing on
C**** NOLV         is current level number
C**** NTLV         is current total number of levels called
C**** NXLV         is maximum number of levels which can be traced
C**** NOSB         is integer number for current subroutine
C**** NTSB         is current total number of subroutines called
C**** NXSB         is maximum no. of subroutines which can be traced
C**** NOSBLV(nolv) is NOSB number at level NOLV
C**** NOSBSM(nosb) is no. of times subroutine no. NOSB has been called
C**** SB(nosb)     is subroutine name for subroutine no. NOSB
C**** TM           is current time
C**** TMST         is trace start time for entry to level 1
C**** TMEL(nolv)   is elapsed time at level NOLV (ie =0 @ entry to NOLV)
C**** TMELSM(nosb) is sum of elapsed times spent in subroutine no. NOSB
C**** TMEN(nolv)   is elapsed time (from start of trace) to level NOLV
C**** TMTLSM(nosb) is total time spent in NOSB & all subs below it

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
C      INCLUDE 'cmiss$reference:ctrl00.cmn'
      INCLUDE 'cmiss$reference:diag00.cmn'
      INCLUDE 'cmiss$reference:trac00.cmn'
!     Parameter List
      CHARACTER NAME*(*)
!     Local Variables
      INTEGER IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,N1SB,NSUB
      REAL TM,VTIME
      CHARACTER COLV*10,ERROR*10
c      EXTERNAL CTRLC_AST

      IF(DIAGNO) THEN
        CALL STRING_TRIM(NAME,IBEG1,IEND1)
        IF(ALLSUB) THEN !turn diagnostics on in all subroutines
          WRITE(OP_STRING,'('' *** Enter '',A)') NAME(IBEG1:IEND1)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ELSE IF(.NOT.ALLSUB) THEN !diags on in selected subrs
          NUM_STACK=NUM_STACK+1
          DOP_STACK(NUM_STACK)=.FALSE.
          DO NSUB=1,NT_SUB
            CALL STRING_TRIM(SUBNAM(NSUB),IBEG2,IEND2)
          IF(NAME(IBEG1:IEND1).EQ.SUBNAM(NSUB)(IBEG2:IEND2))
     '        DOP_STACK(NUM_STACK)=.TRUE.
          ENDDO
          IF(FROMSUB) THEN
            IF(DOP_STACK(NUM_STACK-1)) DOP_STACK(NUM_STACK)=.TRUE.
          ENDIF
          IF(DOP_STACK(NUM_STACK)) THEN
            WRITE(OP_STRING,'('' *** Enter '',A)') NAME(IBEG1:IEND1)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE IF(DOP_STACK(NUM_STACK-1)) THEN
            WRITE(OP_STRING,'('' *** Calls '',A)') NAME(IBEG1:IEND1)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
          DOP=DOP_STACK(NUM_STACK)
        ENDIF
      ENDIF

      IF(TR01) THEN
        IF(NAME.EQ.TRSB) THEN
          TR02=.TRUE.
        ENDIF
        TR01=.FALSE.
        NOLV=NOLV+1
        IF(NOLV.GT.NTLV) THEN
          NTLV=NOLV
        ENDIF
        IF((NOLV.GT.0).AND.(NOLV.LE.NXLV)) THEN
          TM=VTIME()
          IF(NOLV.EQ.1) THEN
            TMST=TM
            IF(TR02) THEN
              WRITE(OP_STRING,'(/''      Time:    Calls:    Level:'//
     '          ' >Subprogram entered'''//
     '          '/''      Time:    Total:   Actual:'//
     '          ' <Subprogram exited'')')
              CALL WRITES(IOTR,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF
          TMEN(NOLV)=TM-TMST
          TMEL(NOLV)=0.0D0
          DO 1 N1SB=1,NTSB
            IF(SB(N1SB).EQ.NAME) THEN
              NOSB=N1SB
              GOTO 2
            ENDIF
    1     CONTINUE
          IF(NTSB.LT.NXSB) THEN
            NTSB=NTSB+1
            SB(NTSB)=NAME
            NOSB=NTSB
          ELSE
            NTSB=NXSB+1
            NOSB=NXSB+1
          ENDIF
    2     NOSBLV(NOLV)=NOSB
          NOSBSM(NOSB)=NOSBSM(NOSB)+1
          IF(TR02) THEN
            WRITE(COLV,'(I10)') NOLV
            CALL STRING_TRIM(COLV,IBEG,IEND)
            WRITE(OP_STRING,'(1X,F10.3,I10,I10,'//COLV(IBEG:IEND)//
     '        '('' >''),A)')
     '        TMEN(NOLV),NOSBSM(NOSB),NOLV,NAME
            CALL WRITES(IOTR,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
        TR01=.TRUE.
      ENDIF

      RETURN
 9999 RETURN 1
      END


      SUBROUTINE EXECUTE_COMFILE(FILE_NAME)

C**** Executes an operating system command file or shell script file

      IMPLICIT NONE
!     Parameter List
      CHARACTER FILE_NAME*(*)
!     Local Variables
      INTEGER IBEG,IEND

      CALL STRING_TRIM(FILE_NAME,IBEG,IEND)
      CALL SYSTEM(FILE_NAME(IBEG:IEND))

      RETURN
      END


      SUBROUTINE FIND_FILE(NOFILE_START,NTFILE,FILE_EXT,FILE_LIST,
     '  ERROR,*)

C**** Finds files in current directory.

      IMPLICIT NONE
!     Parameter List
      INTEGER NOFILE_START,NTFILE
      CHARACTER ERROR*(*),FILE_EXT*(*),FILE_LIST(*)*(*)
!     Local Variables
      INTEGER CLOCAT,CONTEXT,IBEG,IBRA,IDOT,IEND,INDEX,LIB$FIND_FILE,
     '  LIB$FIND_FILE_END,NOFILE
      CHARACTER RESULT*100

      CALL ENTERS('FIND_FILE',*9999)

      CALL EXITS('FIND_FILE')
      RETURN
 9999 CALL ERRORS('FIND_FILE',ERROR)
      CALL EXITS('FIND_FILE')
      RETURN 1
      END


      SUBROUTINE GET_COMMAND_LINE(ARGS,NUMARGS)

C**** Gets the command line arguments.

      IMPLICIT NONE
!     Parameter List
      INTEGER NUMARGS
      CHARACTER ARGS(10)*80
!     Local Variables
      INTEGER I,IARGC

      NUMARGS = IARGC() !SGI Fortran routine
      DO I=1,NUMARGS
        CALL GETARG(I,ARGS(I)) !SGI Fortran routine
      ENDDO

      RETURN
      END


      SUBROUTINE GET_DATE_TIME(DATE,IDATEFMT,ERROR,*)

C**** Uses runtime library call to return date and time
C**** With IDATEFMT equal to:
C****     1  -  returns date in DEC format "DD-MMM-YYYY HH:MM:SS.HH"
C****           i.e 3:45 in the afternoon on 27th of November 1990
C****           becomes "27-NOV-1990 15:45:00.00". (Default).
C****     2  -  returns date as all numbers in IGES format "YYMMDD.HHMMSS"
C****           i.e. the same date is returned as "901127.154500"
C****     3  -  returns a 7 or 8 digit random number based on the time
C****           as "HHSSMMHH" working from hundredths back to hours

      IMPLICIT NONE
!     Parameter List
      INTEGER IDATEFMT
      CHARACTER DATE*(*),ERROR*(*)
!     Local Variables
1      INTEGER MONTHNUM
      CHARACTER MONTH*(3),TEMPDATE*(30)

      CALL ENTERS('GET_DATE_TIME',*9999)

      TEMPDATE='01-JAN-1992 00:00:00.00'
      IF (IDATEFMT.EQ.2) THEN
        MONTH=TEMPDATE(4:6)
        IF (MONTH.EQ.'JAN') THEN
          MONTH='01'
        ELSE IF (MONTH.EQ.'FEB') THEN
          MONTH='02'
        ELSE IF (MONTH.EQ.'MAR') THEN
          MONTH='03'
        ELSE IF (MONTH.EQ.'APR') THEN
          MONTH='04'
        ELSE IF (MONTH.EQ.'MAY') THEN
          MONTH='05'
        ELSE IF (MONTH.EQ.'JUN') THEN
          MONTH='06'
        ELSE IF (MONTH.EQ.'JUL') THEN
          MONTH='07'
        ELSE IF (MONTH.EQ.'AUG') THEN
          MONTH='08'
        ELSE IF (MONTH.EQ.'SEP') THEN
          MONTH='09'
        ELSE IF (MONTH.EQ.'OCT') THEN
          MONTH='10'
        ELSE IF (MONTH.EQ.'NOV') THEN
          MONTH='11'
        ELSE IF (MONTH.EQ.'DEC') THEN
          MONTH='12'
        ENDIF
        DATE=TEMPDATE(10:11)//MONTH(1:2)//TEMPDATE(1:2)//'.'//
     '       TEMPDATE(13:14)//TEMPDATE(16:17)//TEMPDATE(19:20)
      ELSE IF (IDATEFMT.EQ.3) THEN
        DATE=TEMPDATE(22:23)//TEMPDATE(19:20)//TEMPDATE(16:17)
     '    //TEMPDATE(13:14)
      ELSE
        DATE=TEMPDATE
      ENDIF

      CALL EXITS('GET_DATE_TIME')
      RETURN
 9999 CALL ERRORS('GET_DATE_TIME',ERROR)
      CALL EXITS('GET_DATE_TIME')
      RETURN 1
      END


      SUBROUTINE GET_DISPLAY(ERROR,*)

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:disp00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,CLOCAT
      CHARACTER NAME*30

      DATA NAME/' '/

      CALL ENTERS('GET_DISPLAY',*9999)

      !Get display name
      CALL GETENV('DISPLAY',NAME)
      IF(NAME(1:1).EQ.' ') THEN
        DISPLAY_NAME='UNKNOWN'
      ELSE IF(NAME(1:1).EQ.':') THEN
        DISPLAY_NAME=NODE_NAME
      ELSE
        I=CLOCAT(':',NAME)
        DISPLAY_NAME=NAME(1:I-1)
      ENDIF

      TYPE *,'Display name     is ',DISPLAY_NAME

      CALL EXITS('GET_DISPLAY')
      RETURN
 9999 CALL ERRORS('GET_DISPLAY',ERROR)
      CALL EXITS('GET_DISPLAY')
      RETURN 1
      END


      SUBROUTINE GET_REV_TIME(ERROR,*)

C#### Subroutine: GET_REV_TIME
C###  Description:
C###    Returns the current revision time for the CMISS executable.

      IMPLICIT NONE

      INCLUDE 'cmiss$reference:disp00.cmn'

!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER BEG,END,RETURNCODE,STAB(12),STAT
      CHARACTER CTIME*24,EXENAME*50,MODTIME*24


      CALL ENTERS('GET_REV_TIME',*9999)
C
C Get the imagename of the CMISS executable being run
C
      CALL GETARG(0,EXENAME)
      CALL STRING_TRIM(EXENAME,BEG,END)
C Store the name of the executable image (used in fe07)
      IMAGE_NAME=EXENAME(BEG:END)
C
C Find out the file status and information for that executable
C
      RETURNCODE=STAT(EXENAME(BEG:END),STAB)
C
C Convert the modification time into an ascii string
C
      IF(RETURNCODE.EQ.0) THEN
        MODTIME=CTIME(STAB(10))
        TYPE *,'CMISS revision time ',MODTIME
      ENDIF

      CALL EXITS('GET_REV_TIME')
      RETURN
 9999 CALL ERRORS('GET_REV_TIME',ERROR)
      CALL EXITS('GET_REV_TIME')
      RETURN 1
      END


      SUBROUTINE GET_SYSTEM(ERROR,*)

C**** Returns current system nodename & window system.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:disp00.cmn'
      INCLUDE 'cmiss$reference:cmis00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)

      CALL ENTERS('GET_SYSTEM',*9999)

      OS_TYPE='UNIX'
!new
      CALL GETENV('HOST',NODE_NAME)
!old - MPN 30/6/93 - Doesn't seem to work
!      I=GETHOSTNAME(NODE_NAME,LENGTH)
      IF(USE_SOCKET)THEN
C cpb 6/4/95 Need Motif for GX
C        WINDOW_TYPE='TERMINAL'
        WINDOW_TYPE='MOTIF'
      ELSE
        WINDOW_TYPE='MOTIF'
      ENDIF

      TYPE *,'System nodename  is ',NODE_NAME
      TYPE *,'Window system    is ',WINDOW_TYPE
      TYPE *,'Operating system is ',OS_TYPE

      CALL EXITS('GET_SYSTEM')
      RETURN
 9999 CALL ERRORS('GET_SYSTEM',ERROR)
      CALL EXITS('GET_SYSTEM')
      RETURN 1
      END


      SUBROUTINE GETSTR2(A,*)

      IMPLICIT NONE
      CHARACTER A(*)

      RETURN
      END


      SUBROUTINE OPENF(IUNIT,DEVICE,FILE,STATUS,ACCESS,FORM,IRECL,
     '  ERROR,*)

C**** Opens a file.
C**** IUNIT is the name by which the file is referred to in the program.
C**** Valid unit identifiers are non negative integers.
C**** DEVICE is the device on which the file resides.
C**** Valid devices are: 'DISK' 'TERM'
C**** FILE is the filename [filetype [filemode]]
C**** STATUS is 'NEW' or 'OLD' or 'SCRATCH'
C**** ACCESS specifies if the file is 'DIRECT' or 'SEQUEN' or 'APPEND'
C**** Direct access files have a length of 2000 records.
C**** FORM specifies whether the file is 'FORMATTED' or 'UNFORMATTED'
C**** Direct access files have a length of 2000 records.
C**** IRECL is the logical record length
C**** ERROR gives diagnostics in the event of failure.
C**** If an error is detected control is returned to the statement
C**** number of the star.
C**** DEVICE, FILE, ACCESS and ERROR are all character strings.
C**** NOTE: The keyword /RECORDTYPE in OPEN is nonstandard!

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:mxch.inc'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
!     Parameter List
      INTEGER IRECL,IUNIT
      CHARACTER ACCESS*(*),DEVICE*(*),ERROR*(*),FILE*(*),FORM*(*),
     '  STATUS*(*)
!     Local Variables
      INTEGER IBEG,IEND,IOSTAT,IREC
      CHARACTER COMMAND*80,DUMMY*256,
     '  RECL*5,STATUS_UNIX*10,UNIT*4
      DATA RECFM/'F'/,XTENT/'2000'/

      CALL ENTERS('OPENF',*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,'('' Open file '',A)') FILE
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      WRITE(UNIT,'(I4)') IUNIT
      WRITE(RECL,'(I5)') IRECL

      IF(DEVICE.EQ.'TERM') THEN
        OPEN(UNIT=IUNIT,FILE=FILE,STATUS=STATUS,RECL=IRECL,
     '    IOSTAT=IOSTAT)

      ELSE IF(DEVICE.EQ.'DISK') THEN
        IF(ACCESS.EQ.'SEQUEN') THEN
C          STATUS_UNIX=CUPPER(STATUS)
          CALL CUPPER(STATUS,STATUS_UNIX)
          IF(STATUS_UNIX(1:3).EQ.'NEW') STATUS_UNIX='UNKNOWN'
          OPEN(UNIT=IUNIT,FILE=FILE,STATUS=STATUS_UNIX,
     '      ACCESS='SEQUENTIAL',FORM=FORM,IOSTAT=IOSTAT,
     '      RECORDTYPE='STREAM_LF')

        ELSE IF(ACCESS.EQ.'APPEND') THEN
          OPEN(UNIT=IUNIT,FILE=FILE,STATUS=STATUS,
     '      ACCESS='APPEND',FORM=FORM,IOSTAT=IOSTAT,
     '      RECORDTYPE='STREAM_LF')

        ELSE IF(ACCESS.EQ.'DIRECT') THEN
          IF(STATUS.EQ.'OLD') THEN
            OPEN(UNIT=99,FILE=FILE,STATUS=STATUS,
     '        ACCESS='SEQUENTIAL',FORM=FORM,IOSTAT=IOSTAT,
     '        RECORDTYPE='STREAM_LF',
     '        READONLY)
            IF(IOSTAT.EQ.0) THEN
              OPEN(UNIT=IUNIT,STATUS='SCRATCH',
     '    ACCESS='DIRECT',FORM=FORM,IOSTAT=IOSTAT,
     '          RECL=IRECL)
              IREC=0
 10           READ(99,'(A)',IOSTAT=IOSTAT,END=20) DUMMY
              IF(IOSTAT.EQ.0) THEN
C cpb 20/4/94 This inquire doesn't seem to pick up the next record
C properly
C          INQUIRE(UNIT=IUNIT,NEXTREC=IREC)
                IREC=IREC+1
                WRITE(UNIT=IUNIT,FMT='(A)',REC=IREC,IOSTAT=IOSTAT)
     '            DUMMY(1:IRECL)
                GOTO 10
              ELSE ! file error
                GOTO 20
              ENDIF
 20           CLOSE(UNIT=99)
              IF(IOSTAT.EQ.-1) THEN ! EOF
                IOSTAT=0
              ENDIF
            ELSE IF(IOSTAT.EQ.44) THEN ! Temporary - need to inform
               ! user to convert the file.
              WRITE(OP_STRING,'('' >>Convert DIRECT access file to '',
     '          ''a SEQUENTIAL access STREAM_LF file using the TOSLF '',
     '          ''command'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF

          ELSE IF(STATUS.EQ.'NEW') THEN
            OPEN(UNIT=IUNIT,STATUS='SCRATCH',ACCESS='DIRECT',
     '        FORM=FORM,RECL=IRECL,IOSTAT=IOSTAT,RECORDTYPE='FIXED')
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
        WRITE(ERROR,'(I3)') IOSTAT
C CPB 14/6/94 Not sure about these iostats for the sgi
        ERROR=' Iostat='//ERROR(1:3)//' error in OPENF(UNIT='
     '    //UNIT//')'
        IF(IOSTAT.EQ.29.OR.IOSTAT.EQ.2.OR.IOSTAT.EQ.157) THEN
          CALL STRING_TRIM(ERROR,IBEG,IEND)
          ERROR=' File not found: '//ERROR(IBEG:IEND)
        ELSE IF(IOSTAT.EQ.43) THEN
          CALL STRING_TRIM(ERROR,IBEG,IEND)
          ERROR=' Filename invalid: '//ERROR(IBEG:IEND)
        ELSE IF(IOSTAT.EQ.44) THEN
          CALL STRING_TRIM(ERROR,IBEG,IEND)
          ERROR=' Inconsistent record type: '//ERROR(IBEG:IEND)
        ENDIF
        WRITE(OP_STRING,'('' Filename: '',A)') FILE
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        GOTO 9999
      ENDIF

      FILES(IUNIT)=FILE
      FSTATUS(IUNIT)=STATUS(1:1)
      ERROR='  '
      CALL EXITS('OPENF')
      RETURN
 9999 CALL ERRORS('OPENF',ERROR)
      CALL EXITS('OPENF')
      RETURN 1
      END


      SUBROUTINE POST_FILE(FILE_NAME)

C**** Posts file to printer.

      IMPLICIT NONE
!     Parameter List
      CHARACTER FILE_NAME*(*)
!     Local Variables

      RETURN
      END


      SUBROUTINE PURGE_FILE(FILE_NAME)

C**** Purges file versions.

      IMPLICIT NONE
!     Parameter List
      CHARACTER FILE_NAME*(*)
!     Local Variables

      RETURN
      END


      SUBROUTINE SETCMISSCOMMANDPARAMS(PRT1,PRT2,USESOCK,CMVERSION,
     '  CIMAGENAME)

C#### Subroutine: SETCMISSCOMMANDPARAMS
C###  Description:
C###    Sets up the cmiss command parameters obtained from the
C###    command line.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cmis00.cmn'
      INCLUDE 'cmiss$reference:disp00.cmn'
!     Parameter List
      INTEGER PRT1,PRT2,USESOCK,CIMAGENAME(*)
      REAL CMVERSION
!     Local Variables
      INTEGER CSTRLEN,IBEG,IEND
      CHARACTER CMSTRING*5,IMAGENAME*255

      WRITE(CMSTRING,'(F4.2)') CMVERSION
      CALL STRING_TRIM(CMSTRING,IBEG,IEND)
      CMISS=CMSTRING(IBEG:IEND)
      IF(USESOCK.EQ.0) THEN
        USE_SOCKET=.FALSE.
      ELSE
        USE_SOCKET=.TRUE.
        PORT1=PRT1
        PORT2=PRT2
      ENDIF
      CALL CSTRINGLEN(CSTRLEN,CIMAGENAME)
      CALL C2FSTRING(CIMAGENAME,CSTRLEN,IMAGENAME)
      CALL STRING_TRIM(IMAGENAME,IBEG,IEND)
      IMAGE_NAME=IMAGENAME(IBEG:IEND)

      RETURN
      END


      SUBROUTINE SETUPCMISS(ERR)

C#### Subroutine: SETUPCMISS
C###  Description:
C###    Setups up the initial variables for CMISS.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:fsklib.inc'
      INCLUDE 'cmiss$reference:mxch.inc'
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:binf00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:cbfe01.cmn'
      INCLUDE 'cmiss$reference:cbpr00.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:cmis00.cmn'
C      INCLUDE 'cmiss$reference:cspi00.cmn'
C      INCLUDE 'cmiss$reference:ctrl00.cmn'
      INCLUDE 'cmiss$reference:diag00.cmn'
      INCLUDE 'cmiss$reference:dial00.cmn'
      INCLUDE 'cmiss$reference:disp00.cmn'
      INCLUDE 'cmiss$reference:docu00.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:gen000.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:gks000.cmn'
      INCLUDE 'cmiss$reference:gks001.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'cmiss$reference:gtstr00.cmn'
      INCLUDE 'cmiss$reference:head00.cmn'
      INCLUDE 'cmiss$reference:host00.cmn'
      INCLUDE 'cmiss$reference:host00.inc'
      INCLUDE 'cmiss$reference:ktyp00.cmn'
      INCLUDE 'cmiss$reference:mach00.cmn'
      INCLUDE 'cmiss$reference:mach00.inc'
      INCLUDE 'cmiss$reference:map000.cmn'
      INCLUDE 'cmiss$reference:ntsg00.cmn'
      INCLUDE 'cmiss$reference:trac00.cmn'
!     Local Variables
      INTEGER MXCO,MXCOQU,MXSG                !MXCH in 'mxch.inc'
      PARAMETER (MXCO=25,MXCOQU=25,MXSG=2000)
      INTEGER I,ISEG(MXSG),J,NO_COL,NOSG,NSUB,NT_VIEW,ERR
      LOGICAL FIRST_ZOOM,REFINE_ACTIVE,TRANSFORM_ACTIVE,
     '  UPDATE_ZOOM,
     '  DISPLAY_VALUATOR_81,DISPLAY_VALUATOR_82,DISPLAY_VALUATOR_83,
     '  DISPLAY_VALUATOR_84,DISPLAY_VALUATOR_85,DISPLAY_VALUATOR_86,
     '  DISPLAY_VALUATOR_87,DISPLAY_VALUATOR_88,DISPLAY_VALUATOR_89
      CHARACTER CSEG(MXSG)*60,ERROR*(MXCH)

      DATA FPT/0.0,0.0,0.0/,FSHFT/0.0,0.0,0.0/,FSCALE/1.0,1.0,1.0/,
     '  FANGLE/0.0,0.0,0.0/

      DO I=1,100
        DO J=1,MXCH
          OP_STRING(I)(J:J)=' '
        ENDDO
      ENDDO
      GKS=.FALSE.
C      MAPOPN=.FALSE.
      SMG_READ=.TRUE.
c cpb 22/10/95 Set machine constants (don't change)
      CHARSIZE=1
      INTSIZE=8
      SINTSIZE=2
      SPSIZE=4
      DPSIZE=8
      LOGSIZE=4
      ENDIANTYPE=CHAR(MACH_UNKNOWN) !?
      SPFORMTYPE=CHAR(MACH_UNKNOWN) !?
      DPFORMTYPE=CHAR(MACH_UNKNOWN) !?
      MACHTYPE=CHAR(MACH_CRAY) !CRAY
      OSTYPE=CHAR(MACH_UNKNOWN) !?
c
      TR01=.FALSE.
      TR02=.FALSE.
      TR03=.FALSE.
      TR04=.FALSE.
      FIRSTF=.TRUE.
      FIRST_ZOOM =.TRUE.
      UPDATE_ZOOM=.FALSE.
      BUFFER_COUNT=0
      IREC_COMFILE=0
      CALL ENTERS('SETUPCMISS',*9999)
      CALL GET_REV_TIME(ERROR,*9999)
      CALL GET_SYSTEM(ERROR,*9999)
      CALL GET_DISPLAY(ERROR,*9999)
      EXAMPLES_DIR='/product/cmiss/vms/document/examples/'
      IOIP=1 !input  file
      IOOP=2 !usual output
      IODI=3 !diagnostics output
      IOTR=4 !trace output
      IOER=5 !error output
      IOH1=6 !help 1 output
      IOH2=7 !help 2 output
      IOOUT=9 !output FILE
      CALL OPENF(1,'TERM','SYS$INPUT' ,'UNKNOWN',' ',' ',132,
     '  ERROR,*9999)
      CALL OPENF(2,'TERM','SYS$OUTPUT','UNKNOWN',' ',' ',255,
     '  ERROR,*9999)
      CALL OPENF(3,'TERM','SYS$OUTPUT','UNKNOWN',' ',' ',255,
     '  ERROR,*9999)
      CALL OPENF(4,'TERM','SYS$OUTPUT','UNKNOWN',' ',' ',255,
     '  ERROR,*9999)
      CALL OPENF(5,'TERM','SYS$OUTPUT','UNKNOWN',' ',' ',255,
     '  ERROR,*9999)
      CALL OPENF(6,'TERM','SYS$OUTPUT','UNKNOWN',' ',' ',255,
     '  ERROR,*9999)
      CALL OPENF(7,'TERM','SYS$OUTPUT','UNKNOWN',' ',' ',255,
     '  ERROR,*9999)
      CALL OPENF(8,'TERM','SYS$OUTPUT','UNKNOWN',' ',' ',255,
     '  ERROR,*9999)
      IOFI=IOOP !default for listing to file commands
      IFILE=10 !ifile is the file for input files i.e. ip* files
      IOFILE1=11 !output file i.e. op* files
      IOFILE2=12 ! iofile2 -> iofile6 are general purpose output files
      IOFILE3=13
      IOFILE4=14
      IOFILE5=17
      IOFILE6=18
C CPB 30/11/92 the IO# units below should (?) be redundant.
      IO1=1
      IO2=2
      IO3=IOIP
      IO4=IOOP
      IVDU=IOIP
      IMAP=0
      DO NUM_STACK=1,5                !<-|
        DOP_STACK(NUM_STACK)=.FALSE.  !  |
      ENDDO                           !  |
      NUM_STACK=5 !5 levels below FEM !  |- used for diagnostics
      NT_SUB=0                        !  |
      DO NSUB=1,MXSUB                 !  |
        SUBNAM(NSUB)=' '              !  |
      ENDDO                           !  |
      DIAGNO=.FALSE.                  !  |
      ALLSUB=.FALSE.                  !  |
      FROMSUB=.FALSE.                 !  |
      DOP=.FALSE.                     !<-|
      NTSG=0
      NT_KEY=0
      NT_MACRO(0)=0
      MACRO_DEFINE=.FALSE.
      MACRO_COMMAND_EXECUTE=.FALSE.
      MACRO_KEY_EXECUTE=.FALSE.
      DISPLAY_VALUATOR_81=.FALSE.
      DISPLAY_VALUATOR_82=.FALSE.
      DISPLAY_VALUATOR_83=.FALSE.
      DISPLAY_VALUATOR_84=.FALSE.
      DISPLAY_VALUATOR_85=.FALSE.
      DISPLAY_VALUATOR_86=.FALSE.
      DISPLAY_VALUATOR_87=.FALSE.
      DISPLAY_VALUATOR_88=.FALSE.
      DISPLAY_VALUATOR_89=.FALSE.
      REFINE_ACTIVE=.FALSE.
      TRANSFORM_ACTIVE=.FALSE.
      DO_EXAMPLE=.FALSE.
      SELECT_EXAMPLE=.FALSE.
      NT_VIEW=0
      HEADING=' '
C      DO NO_COL=1,9
C        NAME_COL(NO_COL)='UNDEFINED'
C      ENDDO
      DO I=1,99
        IWKS(I)=0
      ENDDO
c     IWKDEF(0)=1
c     IWKDEF(1)=1
      DO NOSG=1,MXSG
        ISEG(NOSG)=0
        CSEG(NOSG)=' '
      ENDDO
      XMIN=-1.0
      XMAX= 1.0
      YMIN=-1.0
      YMAX= 1.0
      ZMIN=-1.0
      ZMAX= 1.0
      DIAG=SQRT(12.0)
      NJT=2
      FILE00='file'
      MENU=.FALSE.

C      RETURN
C 9999 WRITE(*,'('' Fatal Error: '',A)') ERROR
C      ERR=1
C      RETURN
C      END

      CALL EXITS('SETUPCMISS')
      RETURN
 9999 CALL ERRORS('SETUPCMISS',ERROR)
      ERR=1
      CALL EXITS('SETUPCMISS')
      RETURN 1
      END

      SUBROUTINE SLEEPER(IDELAY)

C**** Sleeps for an integer period of seconds

      IMPLICIT NONE
!     Parameter List
      INTEGER IDELAY

      CALL SLEEP(IDELAY)

      RETURN
      END







