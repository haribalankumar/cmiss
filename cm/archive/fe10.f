C#### Module: FE10
C###  Description:
C###    System specific routines (Vax-VMS)

C###  Routine: NEW_HANDLER   (fn) condition handling function
C###  Routine: ADDSTRTOBUFF  Adds a command string to the buffer
C###  Routine: COPY_FILE     copy file between directories
C###  Routine: CLOSEF        close file
C###  Routine: CMINITIALISE Initialises global CM variables
C###  Routine: CTRLC_AST     set CTRLC=true when asynch CTRL_C issued
C###  Routine: DIALOG        control command parsing environment
C###  Routine: DO_COMMAND    executes an operating system command
C###  Routine: ENTERS        trace entry to a subprogram
C###  Routine: EXECUTE_COMFILE execute an op system command file
C###  Routine: FIND_FILE     *** ARCHIVED ***
C###  Routine: GET_CMISS_EXAMPLES return CMISS examples directory
C###  Routine: GET_COMMAND_LINE  get the command line arguments
C###  Routine: GET_DATE_TIME return date and time
C###  Routine: GET_DISPLAY   return DecWindows display
C###  Routine: GET_REV_TIME  return CMISS revision time
C###  Routine: GET_SYSTEM    return system parameters
C###  Routine: GETSTR1       return user string from keyboard mapping
C###  Routine: GETSTR2       return interrupt key number
C###  Routine: GETSTR3       return key for DOCUM
C###  Routine: MAINCMLOOP    Main CMISS Loop.
C###  Routine: OPENF         open file
C###  Routine: POST_FILE     *** ARCHIVED ***
C###  Routine: PURGE_FILE    *** ARCHIVED ***
C###  Routine: RESET_HANDLER reset old error handler
C###  Routine: SET_CTRLC_AST set the asynchronous CTRL/C handler
C###  Routine: SET_HANDLER   sets new error handler
C###  Routine: SETUPCMISS    setup initial CMISS constants
C###  Routine: SLEEPER       sleeps for an integer period of seconds

      INTEGER*4 FUNCTION NEW_HANDLER(SIGARGS,MECHARGS)

C#### Function: NEW_HANDLER
C###  Type: INTEGER*4
C###  Description:
C###   NEW_HANDLER is the condition handler (written by Nick Burke).
C###   Note: Using the stack frame depth in mechrgs(3) as the first
C###   argument of SYS$UNWIND causes execution to resume at the
C###   location specified by the program counter in the call frame of
C###   the procedure that receives control after the unwind has
C###   completed.  Setting the first arg of sys$unwind to 0 will cause
C###   execution to resume at the location specified by the program
C###   counter in the call frame of the procedure that called the
C###   procedure that established the condition handler.

C**** CPB 23/10/93 Changing function, MECHARGS and SIGARGS types from
C**** INTEGER*2 to INTEGER*4 for use on OpenVMS AXP systems

      IMPLICIT NONE
      INCLUDE '($CHFDEF)'
      INCLUDE '($SSDEF)'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
!     Parameter List
      INTEGER*4 SIGARGS(*)
C CPB 23/10/93 changing mechargs def so system constants can be used on
C both the VMS and OpenVMS systems.
      RECORD /CHFDEF2/ MECHARGS
!     Local Variables
      CHARACTER ERROR*10
      INTEGER INDEX,LIB$MATCH_COND,SEVERITY
      LOGICAL*4 SYS$UNWIND,STATUS

C     Check to see if a stack unwind is already in progress
C     If it is, do nothing
      INDEX = LIB$MATCH_COND(SIGARGS(2),SS$_UNWIND)

      IF(index.ne.1) THEN
C       Mask off all but the 4 ls bits (ie the severity code)
        SEVERITY=JIAND(SIGARGS(2),7)
        WRITE(OP_STRING,'('' >>Error severity:'',I5)') SEVERITY
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' >>Sub call depth:'',I5)')
     '    MECHARGS.CHF$IS_MCH_DEPTH
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)

C       Unwind the stack if it is a severe error
        IF(SEVERITY.EQ.4) THEN
          WRITE(OP_STRING,'('' >>Severe error:'')')
          CALL WRITES(IOER,OP_STRING,ERROR,*9999)
C         Display the error message
          SIGARGS(1) = SIGARGS(1)-2    !subtract last two elements
          CALL SYS$PUTMSG(SIGARGS,,,)
          SIGARGS(1) = SIGARGS(1) + 2  !restore last two elements
C         Then unwind the stack
C CPB 23/10/93 changing mechargs def so system constants can be used
C on both the VMS and OpenVMS systems.
          STATUS=SYS$UNWIND(MECHARGS.CHF$IS_MCH_DEPTH,)
          IF(status.ne.SS$_NORMAL) THEN
            WRITE(OP_STRING,'('' >>Fatal error, unwind failed'')')
            CALL WRITES(IOER,OP_STRING,ERROR,*9999)
            IF(STATUS .EQ. SS$_ACCVIO) THEN
              WRITE(OP_STRING,'('' >>STATUS: access violation'')')
              CALL WRITES(IOER,OP_STRING,ERROR,*9999)
            ELSE IF(STATUS .EQ. SS$_INSFRAME) THEN
              WRITE(OP_STRING,'('' >>STATUS: insufficient call '
     '          //'frames'')')
              CALL WRITES(IOER,OP_STRING,ERROR,*9999)
            ELSE IF(STATUS .EQ. SS$_NOSIGNAL) THEN
              WRITE(OP_STRING,'('' >>STATUS: no signal active'')')
              CALL WRITES(IOER,OP_STRING,ERROR,*9999)
            ELSE IF(STATUS .EQ. SS$_UNWINDING) THEN
              WRITE(OP_STRING,'('' >>STATUS: unwind already in '
     '          //'progress'')')
              CALL WRITES(IOER,OP_STRING,ERROR,*9999)
            ELSE
              WRITE(OP_STRING,'('' >>STATUS: unknown error:'',I5)')
     '          status
              CALL WRITES(IOER,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF

        ELSE
          NEW_HANDLER = SS$_RESIGNAL
        ENDIF
      ENDIF

 9999 RETURN
      END


      SUBROUTINE ADDSTRTOBUFF(STRING,ERROR,*)

C#### Subroutine: ADDSTRTOBUFF
C###  Description:
C###    ADDSTRTOBUFF adds the command string STRING to the buffer of
C###    command strings

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:gtstr00.cmn'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER IBEG,IEND,no_buffer

      CALL ENTERS('ADDSTRTOBUFF',*9999)

      BUFFER_COUNT=BUFFER_COUNT+1
      IF(BUFFER_COUNT.GT.99) BUFFER_COUNT=99
      DO no_buffer=BUFFER_COUNT,2,-1
        BUFFER(no_buffer)=BUFFER(no_buffer-1)
        BUFFER_LENGTH(no_buffer)=BUFFER_LENGTH(no_buffer-1)
      ENDDO
      BUFFER(1)=' >'//STRING
      CALL STRING_TRIM(BUFFER(1),IBEG,IEND)
      BUFFER_LENGTH(1)=IEND

      CALL EXITS('ADDSTRTOBUFF')
      RETURN
 9999 CALL ERRORS('ADDSTRTOBUFF',ERROR)
      CALL EXITS('ADDSTRTOBUFF')
      RETURN 1
      END


      SUBROUTINE COPY_FILE(FILE_NAME)

C#### Subroutine: COPY_FILE
C###  Description:
C###    COPY_FILE copies file between directories.

      IMPLICIT NONE
!     Parameter List
      CHARACTER FILE_NAME*(*)
!     Local Variables
      INTEGER IBEG,IEND,ISTATUS,LIB$SPAWN

      CALL STRING_TRIM(FILE_NAME,IBEG,IEND)
      ISTATUS=LIB$SPAWN('copy esv1$dkb200:[user.cmiss.doc]'
     '  //FILE_NAME(IBEG:IEND)//'.COM *.*',,,,,,,,,,,)

      RETURN
      END


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
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:gino00.cmn'
!     Parameter List
      INTEGER IUNIT
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER CLOCAT,i,IBEG,IEND,iNEW,iOLD,IOSTAT,IREC
      CHARACTER ACCESS*10,FILE*200,FILENAME*200,
     '  FORM*11,ISTAT*5,LINE*256
      LOGICAL EXIST,OPENED

      CALL ENTERS('CLOSEF',*9999)
      INQUIRE(UNIT=IUNIT,NAME=FILE,IOSTAT=IOSTAT,EXIST=EXIST,
     '  OPENED=OPENED,ACCESS=ACCESS,FORM=FORM,NEXTREC=IREC)
      IF(DOP) THEN
        WRITE(OP_STRING,'('' UNIT='',I4,'' FILE='',A)') IUNIT,FILE
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' EXIST='',L1,'' OPENED='',L1)') EXIST,OPENED
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' ACCESS='',A,'' FORMAT='',A)') ACCESS,FORM
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' NEXTREC='',I5)') IREC
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(EXIST.AND.OPENED) THEN

        IF(ACCESS(1:6).EQ.'DIRECT') THEN
C        Find end of file
c         IOSTAT=0
c         DO WHILE(IOSTAT.EQ.0)
c           irec=irec+1
c           READ(UNIT=IUNIT,REC=irec,FMT='(A)',IOSTAT=IOSTAT) LINE(1:80)
c         ENDDO
C        If file is new or modified, copy from direct access
C        scratch file to sequential stream_lf
          iNEW=CLOCAT('NEW_',FILE)
          iOLD=CLOCAT('OLD_',FILE)
          IF(iNEW.ne.0) THEN !copy new file
            CALL STRING_TRIM(FILE,IBEG,IEND)
            FILENAME=FILE(iNEW+4:IEND)
            IF(DOP) THEN
              WRITE(OP_STRING,'('' Copying scratch file to : '',A)')
     '          FILENAME
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            OPEN(UNIT=99,FILE=FILENAME,STATUS='NEW',
     '        ACCESS='SEQUENTIAL',FORM='FORMATTED',IOSTAT=IOSTAT,
     '        CARRIAGECONTROL='NONE',RECORDTYPE='STREAM_LF')
            IF(IOSTAT.EQ.0) THEN
              DO i=1,IREC-1
                READ(UNIT=IUNIT,REC=I,FMT='(A)',IOSTAT=IOSTAT) LINE
                IF(IOSTAT.EQ.0) THEN
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
            IF(iostat.ne.0) GOTO 9998

          ELSE IF(iOLD.ne.0.AND.Modify_File) THEN !copy old file
            CALL STRING_TRIM(FILE,IBEG,IEND)
            FILENAME=FILE(iOLD+4:IEND)
            IF(DOP) THEN
              WRITE(OP_STRING,'('' Copying scratch file to : '',A)')
     '          FILENAME
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            OPEN(UNIT=99,FILE=FILENAME,STATUS='NEW',
     '        ACCESS='SEQUENTIAL',FORM='FORMATTED',IOSTAT=IOSTAT,
     '        CARRIAGECONTROL='NONE',RECORDTYPE='STREAM_LF')
            IF(IOSTAT.EQ.0) THEN
              DO i=1,IREC-1
                READ(UNIT=IUNIT,REC=I,FMT='(A)',IOSTAT=IOSTAT) LINE
                IF(IOSTAT.EQ.0) THEN
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
            IF(iostat.ne.0) GOTO 9998
          ENDIF !iold

        ENDIF !direct access
        CLOSE(UNIT=IUNIT,IOSTAT=IOSTAT)
        IF(iostat.ne.0) GOTO 9998
      ENDIF !exist & opened

      ERROR=' '
      CALL EXITS('CLOSEF')
      RETURN
 9998 WRITE(ISTAT,'(I3)') IOSTAT
      ERROR='>>Error: Iostat = '//ISTAT(1:3)
 9999 CALL ERRORS('CLOSEF',ERROR)
      CALL EXITS('CLOSEF')
      RETURN 1
      END


      SUBROUTINE CMINITIALISE(ERR)

C#### Subroutine: CMINITIALISE
C###  Description:
C###    Sets up the initial global variables for CMISS.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:mxch.inc'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbfe01.cmn'
      INCLUDE 'cmiss$reference:cmis00.cmn'
C      INCLUDE 'cmiss$reference:cspi00.cmn'
      INCLUDE 'cmiss$reference:diag00.cmn'
      INCLUDE 'cmiss$reference:docu00.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:gks000.cmn'
      INCLUDE 'cmiss$reference:gtstr00.cmn'
      INCLUDE 'cmiss$reference:head00.cmn'
      INCLUDE 'cmiss$reference:map000.cmn'
C$    INCLUDE 'cmiss$reference:mp00.cmn'
      INCLUDE 'cmiss$reference:ntsg00.cmn'

!     SMAR009 removed 16 Dec 1998
!     INCLUDE 'cmiss$reference:tol00.cmn'

      INCLUDE 'cmiss$reference:trac00.cmn'
!     Parameter List
      INTEGER ERR
!     Local Variables
      INTEGER macro_command,no_col,
     '  nsub
      REAL*8 DLAMCH
      CHARACTER ERROR*(MXCH)

C     Initialise any variables used by enter/exit
      TR01=.FALSE.
      TR02=.FALSE.
      TR03=.FALSE.
      TR04=.FALSE.

      FILE00='file'
      PATH00=' '

C GMH 18/2/97 Initialise BLANK to spaces (not NULL's)
      DO no_col=1,80
        BLANK(no_col:no_col)=' '
      ENDDO !no_col
      GKS=.FALSE.
C      MAPOPN=.FALSE.
      NUM_LIBRARY=0

      FIRST_SYNTAX=.TRUE.
      FEM_ARRAYS_MALLOCED=.FALSE.
      BUFFER_COUNT=0
      DO com=1,20
        IREC_COMFILE(COM_UNIT)=-1
      ENDDO
      COM_UNIT=75
      NEST=0
      CALL GET_REV_TIME(ERROR,*9999)
      CALL GET_SYSTEM(ERROR,*9999)
      CALL GET_DISPLAY(ERROR,*9999)
      IF(USE_SOCKET) THEN
        IOIP=1 !input  file
        IOOP=2 !usual output
        IODI=3 !diagnostics output
        IOTR=4 !trace output
        IOER=5 !error output
        IOH1=6 !help 1 output
        IOH2=7 !help 2 output
      ELSE
        IOIP=5 !input  file
        IOOP=6 !usual output
        IODI=6 !diagnostics output
        IOTR=6 !trace output
        IOER=6 !error output
        IOH1=6 !help 1 output
        IOH2=6 !help 2 output
        CALL OPENF(5,'TERM','SYS$INPUT' ,'UNKNOWN',' ',' ',132,
     '    ERROR,*9999)
        CALL OPENF(6,'TERM','SYS$OUTPUT','UNKNOWN',' ',' ',255,
     '    ERROR,*9999)
      ENDIF
      IOOUT=9 !output FILE

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
      ECHO_OUTPUT=.FALSE.
      DO num_stack=1,5                !<-|
        DOP_STACK(num_stack)=.FALSE.  !  |
      ENDDO                           !  |
      num_stack=5 !5 levels below FEM !  |- used for diagnostics
      NT_SUB=0                        !  |
      DO nsub=1,MXSUB                 !  |
        SUBNAM(nsub)=' '              !  |
      ENDDO                           !  |
      DIAGNO=.FALSE.                  !  |
      ALLSUB=.FALSE.                  !  |
      FROMSUB=.FALSE.                 !  |
      DOP=.FALSE.                     !  |
C$    THREAD_NUM=-1                   !<-|
      NTSG=0
      NT_KEY=0
      NT_MACRO(0)=0
      NT_MACRO_names(0)=0
      DO macro_command=1,9
        NT_MACRO_names(macro_command)=0
        MACRO_names(macro_command)=' '
      ENDDO
      MACRO_DEFINE=.FALSE.
      MACRO_COMMAND_EXECUTE=.FALSE.
      MACRO_KEY_EXECUTE=.FALSE.
      DO_EXAMPLE=.FALSE.
      SELECT_EXAMPLE=.FALSE.
      HEADING=' '
C      DO no_col=1,9
C        NAME_COL(no_col)='UNDEFINED'
C      ENDDO

C     Initialise any platform-dependent variables
      CALL SETUPCMISS(ERROR,*9999)
C     Initialise the cm-cmgui link
      CALL CMGUI_LINK_INITIALISE(ERROR,*9999)
C     Initialise any general variables

      RETURN
 9999 WRITE(*,'('' Fatal Error: '',A)') ERROR
      ERR=1
      RETURN
      END


       SUBROUTINE CTRLC_AST()

C#### Subroutine: CTRLC_AST
C###  Description:
C###    CTRLC_AST sets CTRLC to .TRUE. to indicate that CTRL/C has
C###    been pressed and then reestablishes the asynchrouns ctrlc
C###    handler with a call to SET_CTRLC_AST.

C**** CPB 23/10/93 Moved this subroutine from FE01

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:ctrl00.cmn'

      CTRLC=.TRUE.

C CPB 23/10/93 Resetting the CTRLC_AST routine here to avoid generating
C CTRL/Y the second time CTRL/C is pressed. The resetting has to be done
C in a separate subroutine to avoid recursion.

      CALL SET_CTRLC_AST

      RETURN
      END


C LC not used ? 26/2/97

C      SUBROUTINE DIALOG
C
CC#### Subroutine: DIALOG
CC###  Description:
CC###    <HTML>
CC###    DIALOG is the dialog interpreter.
CC###    Note: pi and PI are assigned as user-defined names
CC###    <PRE>
CC###    Condition Handler References:
CC###       ch2 Introduction to System Routines
CC###       ch4 LIB$ Manual -An Overview of VAX Condition Handling
CC###       ch9 VAX Fortran User Manual
CC###       ch8 VAX C Run time library reference manual
CC###    section 8.4 VAX Debugger Manual
CC###       ch9 Guide to Programming Resources ***THIS IS THE BEST***
CC###       ch10 Introduction to VMS System Services
CC###    </PRE>
CC###    ISTATUS is >0 on successful completion of RTL call
CC###    (ISTATUS=.true.).
CC###    </HTML>
C
C      IMPLICIT NONE
C      INCLUDE '($IODEF)'
C      INCLUDE '($JPIDEF)'
C      INCLUDE '($SMGDEF)'
C      INCLUDE '($TRMDEF)'
C      INCLUDE 'cmiss$reference:fsklib.inc'
C      INCLUDE 'cmiss$reference:mxch.inc'
C      INCLUDE 'cmiss$reference:b00.cmn'
C      INCLUDE 'cmiss$reference:b01.cmn'
C      INCLUDE 'cmiss$reference:binf00.cmn'
C      INCLUDE 'cmiss$reference:cbdi02.cmn'
C      INCLUDE 'cmiss$reference:cbdi10.cmn'
C      INCLUDE 'cmiss$reference:cbfe01.cmn'
C      INCLUDE 'cmiss$reference:cbpr00.cmn'
C      INCLUDE 'cmiss$reference:cbwk01.cmn'
C      INCLUDE 'cmiss$reference:cmis00.cmn'
C      INCLUDE 'cmiss$reference:cspi00.cmn'
C      INCLUDE 'cmiss$reference:ctrl00.cmn'
C      INCLUDE 'cmiss$reference:diag00.cmn'
C      INCLUDE 'cmiss$reference:dial00.cmn'
C      INCLUDE 'cmiss$reference:disp00.cmn'
C      INCLUDE 'cmiss$reference:docu00.cmn'
C      INCLUDE 'cmiss$reference:file00.cmn'
C      INCLUDE 'cmiss$reference:gen000.cmn'
C      INCLUDE 'cmiss$reference:geom00.cmn'
C      INCLUDE 'cmiss$reference:gks000.cmn'
C      INCLUDE 'cmiss$reference:gks001.cmn'
C      INCLUDE 'cmiss$reference:graf00.cmn'
C      INCLUDE 'cmiss$reference:gtstr00.cmn'
C      INCLUDE 'cmiss$reference:head00.cmn'
C      INCLUDE 'cmiss$reference:host00.cmn'
C      INCLUDE 'cmiss$reference:iter00.cmn'
C      INCLUDE 'cmiss$reference:ktyp00.cmn'
C      INCLUDE 'cmiss$reference:mach00.cmn'
C      INCLUDE 'cmiss$reference:map000.cmn'
C      INCLUDE 'cmiss$reference:ntsg00.cmn'
C      INCLUDE 'cmiss$reference:trac00.cmn'
C!     Local Variables
C      INTEGER MXCO,MXCOQU,MXSG
C      PARAMETER (MXCO=25,MXCOQU=25,MXSG=10000)
C      INTEGER ERR,i,IBEG,IBEG4,ID,ID_DEVICE,ID_STATUS,ID_WS,IEND,IEND4,
C     '  IFROMC,INDEX_PLIN,INLIST,INPUT_CHOICE,INTSTR(MXCH),IRET,
C     '  ISEG(MXSG),ISEGM,ISEGM_LIST(0:100),ISTATUS,iter,IW,LEN,
C     '  LIB$ESTABLISH,LIB$GETJPI,LIB$SPAWN,macro_command,
C     '  MODE_PROJ,n,N1LIST,NCENTRE,NCHAR,NEW_HANDLER,nj,nhost,NOCH,
C     '  no_col,NO_GROUP_TRANSFORM,nogrpl,nomacro,NOPTS,
C     '  nosg,nsub,NTCH,NTFILE,NT_VIEW,RETURNVAL
C      REAL DX,DY,FANGLE(3),FPT(3),FSCALE(3),FSHFT(3),
C     '  R4DATA(2),XNDC1,XNDC2,XNDC3,XNDC4,VALUE,VMATRIX(4,4),
C     '  X0,XCENTRE,Y0,YCENTRE
C      LOGICAL CBBREV,CHANGE,CONTINUE,DEFINE,END,FIRST_ZOOM,
C     '  REFINE_ACTIVE,TRANSFORM_ACTIVE,UPDATE,UPDATE_ZOOM,
C     '  DISPLAY_VALUATOR_81,DISPLAY_VALUATOR_82,DISPLAY_VALUATOR_83,
C     '  DISPLAY_VALUATOR_84,DISPLAY_VALUATOR_85,DISPLAY_VALUATOR_86,
C     '  DISPLAY_VALUATOR_87,DISPLAY_VALUATOR_88,DISPLAY_VALUATOR_89
C      CHARACTER CFROMI*5,CFROMR*11,CHAR3*3,CHAR4*4,CHAR5*5,CHOOSE*40,
C     '  CIW*1,CLASS*8,COD(MXCO)*90,COQUD(MXCO,MXCOQU)*30,
C     '  CSEG(MXSG)*60,ERROR*(MXCH),FILE_NAME*50,
C     '  OPTION(30)*40,
C     '  OPTION1(19)*11,OPTION2(12)*11,OPTION3(26)*11,OPTION4(11)*11,
C     '  OPTION5(9)*20,OPTION10(10)*20,OPTION11(10)*20,
C     '  PARAMETER_TYPE*20,Q*1,SDATA*10,STATSTR*80,STRG*(MXCH),
C     '  STRING*(MXCH),TEXT_STRING*80,TRANSFORM_TYPE*7
C
C      EXTERNAL CTRLC_AST,NEW_HANDLER
C      DATA FPT/0.0,0.0,0.0/,FSHFT/0.0,0.0,0.0/,FSCALE/1.0,1.0,1.0/,
C     '  FANGLE/0.0,0.0,0.0/
C
C      CALL ENTERS('DIALOG',*9999)
C      DO i=1,100
C        OP_STRING(i)(1:1)=CHAR(0)
C      ENDDO
C      GKS=.FALSE.
C      MAPOPN=.FALSE.
C      SMG_READ=.TRUE.
C      INPUT_CHAN=0
C      NUM_LIBRARY=0
Cc cpb 22/10/95 Set machine constants (don't change)
C      CHARSIZE=1
C      INTSIZE=4
C      SINTSIZE=2
C      SPSIZE=4
C      DPSIZE=8
C      LOGSIZE=4
C      ENDIANTYPE=CHAR(MACH_LITTLEENDIAN) !Little endian native
C      SPFORMTYPE=CHAR(MACH_SPIEEE) !IEEE single precision format
C      DPFORMTYPE=CHAR(MACH_DPIEEE) !IEEE double precision format
C      MACHTYPE=CHAR(MACH_DECALPHA) !Alpha
C      OSTYPE=CHAR(MACH_VMS) !VMS
Cc
C      TR01=.FALSE.
C      TR02=.FALSE.
C      TR03=.FALSE.
C      TR04=.FALSE.
C      FIRST_SYNTAX=.TRUE.
C      FIRST_ZOOM =.TRUE.
C      UPDATE_ZOOM=.FALSE.
C      BUFFER_COUNT=0
C      IREC_COMFILE=0
C      CALL GET_REV_TIME(ERROR,*9999)
C      CALL GET_SYSTEM(ERROR,*9999)
C      CALL GET_DISPLAY(ERROR,*9999)
C      IOIP=1 !input  file
C      IOOP=2 !usual output
C      IODI=3 !diagnostics output
C      IOTR=4 !trace output
C      IOER=5 !error output
C      IOH1=6 !help 1 output
C      IOH2=7 !help 2 output
C      IOOUT=9 !output FILE
C      CALL OPENF(1,'TERM','SYS$INPUT' ,'NEW',' ',' ',132,ERROR,*9999)
C      CALL OPENF(2,'TERM','SYS$OUTPUT','NEW',' ',' ',255,ERROR,*9999)
C      CALL OPENF(3,'TERM','SYS$OUTPUT','NEW',' ',' ',255,ERROR,*9999)
C      CALL OPENF(4,'TERM','SYS$OUTPUT','NEW',' ',' ',255,ERROR,*9999)
C      CALL OPENF(5,'TERM','SYS$OUTPUT','NEW',' ',' ',255,ERROR,*9999)
C      CALL OPENF(6,'TERM','SYS$OUTPUT','NEW',' ',' ',255,ERROR,*9999)
C      CALL OPENF(7,'TERM','SYS$OUTPUT','NEW',' ',' ',255,ERROR,*9999)
C      CALL OPENF(8,'TERM','SYS$OUTPUT','NEW',' ',' ',255,ERROR,*9999)
C      IOFI=IOOP !default for listing to file commands
C      IFILE=10 !ifile is the file for input files i.e. ip* files
C      IOFILE1=11 !output file i.e. op* files
C      IOFILE2=12 ! iofile2 -> iofile6 are general purpose output files
C      IOFILE3=13
C      IOFILE4=14
C      IOFILE5=17
C      IOFILE6=18
CC CPB 30/11/92 the IO# units below should (?) be redundant.
C      IO1=1
C      IO2=2
C      IO3=IOIP
C      IO4=IOOP
C      IVDU=IOIP
C      IMAP=0
C      ECHO_OUTPUT=.FALSE.
C      DO num_stack=1,5                !<-|
C        DOP_STACK(num_stack)=.FALSE.  !  |
C      ENDDO                           !  |
C      num_stack=5 !5 levels below FEM !  |- used for diagnostics
C      NT_SUB=0                        !  |
C      DO nsub=1,MXSUB                 !  |
C        SUBNAM(nsub)=' '              !  |
C      ENDDO                           !  |
C      DIAGNO=.FALSE.                  !  |
C      ALLSUB=.FALSE.                  !  |
C      FROMSUB=.FALSE.                 !  |
C      DOP=.FALSE.                     !<-|
C      NTSG=0
C      NT_KEY=0
C      NT_MACRO(0)=0
C      NT_MACRO_names(0)=0
C      DO macro_command=1,9
C        NT_MACRO_names(macro_command)=0
C        MACRO_names(macro_command)=' '
C      ENDDO
C      MACRO_DEFINE=.FALSE.
C      MACRO_COMMAND_EXECUTE=.FALSE.
C      MACRO_KEY_EXECUTE=.FALSE.
C      REFINE_ACTIVE=.FALSE.
C      TRANSFORM_ACTIVE=.FALSE.
C      DO_EXAMPLE=.FALSE.
C      SELECT_EXAMPLE=.FALSE.
C      NT_VIEW=0
C      HEADING=' '
C      DO no_col=1,9
C        NAME_COL(no_col)='UNDEFINED'
C      ENDDO
C      DO i=1,99
C        IWKS(i)=0
C      ENDDO
Cc     IWKDEF(0)=1
Cc     IWKDEF(1)=1
C      DO nosg=1,MXSG
C        ISEG(nosg)=0
C        CSEG(nosg)=' '
C      ENDDO
C      XMIN=-1.0
C      XMAX= 1.0
C      YMIN=-1.0
C      YMAX= 1.0
C      ZMIN=-1.0
C      ZMAX= 1.0
C      DIAG=SQRT(12.0)
C      NJT=2
C      FILE00='file'
C      PATH00=' '
C      MENU=.FALSE.
C      COD(2)='PI'
C      COD(3)=CFROMR(PI,'(E11.5)')
C      CALL ASSIGN(STRING,1,COD,COQUD,ERROR,*9999)
C!     COD(2)='pi'  !interferes with pick node command (PJH 10-feb-92)
C!     COD(3)=CFROMR(PI,'(E11.5)')
C!     CALL ASSIGN(STRING,1,COD,COQUD,ERROR,*9999)
C
C      IF(USE_SOCKET.EQ..FALSE.) THEN
C
CC ***   Query mode
C        ISTATUS = LIB$GETJPI(JPI$_MODE,,,,STATSTR,LEN)
C        IF(.NOT.ISTATUS) CALL LIB$SIGNAL(%VAL(ISTATUS))
C
CC ***   Check whether batch mode
C        IF(STATSTR(1:1).EQ.'B') THEN !Cannot use SMG read routines
C          SMG_READ=.FALSE.
C        ELSE                        !Set up control-c trap
C          CALL SET_CTRLC_AST
C        ENDIF
C
CC CPB 26/10/95 Condition handler now changed in syntax
CC ***   Establish condition handler
C
C        FATAL_HANDLER=.TRUE.
C        CHANGE_HANDLER=.FALSE.
C
CCC CPB 6/10/93 Don't store old handler as OpenVMS cannot handle
CCC the typing of the procedure. Use LIB$REVERT instead.
C        OLD_HANDLER=LIB$ESTABLISH(NEW_HANDLER) !address of deflt handler
CC
C
CC ***   Check for cmiss.com file
C        END=.FALSE.
CC news MPN 6-Jul-95: changing "read file;com" to "read com;file" etc
CC old        STRG='read cmiss;com'
C!PJH 5dec95        STRG='read com;cmiss'
C!PJH 5dec95        CALL PARSE(0,ISEG,CSEG,STRG,END,ERROR,*9999)
C
C        IF(END) THEN !quit issued from cmiss.com
C          CONTINUE=.FALSE.
C        ELSE
C          CONTINUE=.TRUE.
C        ENDIF
C        DO WHILE(CONTINUE) !is main program loop
C
CC GMH 15/11/95 making call conditional
C          IF(USE_GRAPHICS.EQ.1) THEN
C            CALL GXWAIT(0.0,ERR) !Update graphics
C          ENDIF
C
C          IF(MACRO_COMMAND_EXECUTE) THEN !parse command macros
C            MACRO_COMMAND_EXECUTE=.FALSE.
C            DO nomacro=1,NT_MACRO_names(MACRO_command_ID)
C              STRG=MACRO_COMMAND_buffer(nomacro,MACRO_command_ID)
C              CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*150)
C            ENDDO
C          ELSE IF(MACRO_KEY_EXECUTE) THEN !parse key macros
C            MACRO_KEY_EXECUTE=.FALSE.
C            DO nomacro=1,NT_MACRO(macro_key)
C              STRG=MACRO_KEY_buffer(nomacro,macro_key)
C              CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*150)
C            ENDDO
C          ELSE IF(ITERATE_COMMAND_EXECUTE) THEN !parse command macros
C            CALL STRING_TRIM(FILE00,IBEG,IEND)
C            DO iter_counter=1,NT_ITERATE
C              !Note that iter_counter is kept in iter00.cmn
C              CHAR4=CFROMI(iter_counter,'(I4)')
C              CALL STRING_TRIM(CHAR4,IBEG4,IEND4)
C              FILE00=FILE00(IBEG:IEND)//'_'//CHAR4(IBEG4:IEND4)
C              DO nomacro=1,NT_MACRO_names(MACRO_command_ID)
C                STRG=MACRO_COMMAND_buffer(nomacro,MACRO_command_ID)
C                CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*150)
C              ENDDO !nomacro
C            ENDDO !iter_counter
C            ITERATE_COMMAND_EXECUTE=.FALSE.
C          ELSE IF(DO_EXAMPLE) THEN !read selected example file
C            DO_EXAMPLE=.FALSE.
C            CALL STRING_TRIM(EXAMPLE_NAME,IBEG,IEND)
C            STRING=' > read '//EXAMPLE_NAME(IBEG:IEND)//';com;doc'
C            WRITE(OP_STRING,'(A)') STRING(1:30)
C            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C          ELSE IF(SELECT_EXAMPLE) THEN !open selected example file
C            SELECT_EXAMPLE=.FALSE.
C            CALL STRING_TRIM(EXAMPLE_NAME,IBEG,IEND)
C            STRING=' > open '//EXAMPLE_NAME(IBEG:IEND)//';com;doc'
C            WRITE(OP_STRING,'(A)') STRING(1:30)
C            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C          ELSE !set prompt and get string from command line
C            CALL SETSTR(STRING,ERROR,*150)
C            IF(SMG_READ) THEN
C              CALL GETSTR1(STRING,PR(:LNPR),LNPR,END,ERROR,*150)
C            ELSE
C              CALL GETSTR(STRING(LNPR-2:),PR(:LNPR),END,ERROR,*150)
C            ENDIF
C          ENDIF
C
C          IF(END) THEN       !CTRL-Z has been used to quit
C            CALL QUIT(END,ERROR,*150)
C            CONTINUE=.FALSE.
C            GOTO 200
C          ENDIF
C
CC ***     Parse string
C          IF(SMG_READ)THEN
C            STRING=STRING(3:)
C            CALL PARSE(0,ISEG,CSEG,STRING,END,ERROR,*150)
C          ELSE
C            CALL PARSE(0,ISEG,CSEG,STRING,END,ERROR,*150)
C          ENDIF
C          STRING=' ' !Test 4-3-1991
C
CC ***     Check for changes to condition handler
C          IF(CHANGE_HANDLER) THEN
C            IF(FATAL_HANDLER) THEN
C              CALL LIB$ESTABLISH(NEW_HANDLER) !establishs new condition handler
C            ELSE IF(.NOT.FATAL_HANDLER) THEN
CC
CC CPB 13/7/93 - OpenVMS has changed the type of the argument to
CC LIB$ESTABLISH. As a result changing back to the old condition handler
CC can not be handled in this way (OLD_HANDLER must be a procedure
CC value). Removing this call to LIB$ESTABLISH.
CC
CC              CALL LIB$ESTABLISH(OLD_HANDLER) !establishs VAX condition handler
CC
CC CPB 6/10/93 - Changing the way the old handler is restored. Use a LIB$REVERT
CC rather than storing the address of the old handler. This should now be
CC compatible with OpenVMS
CC
CC              WRITE(OP_STRING,'(''>>Can not reset condition '',
CC     '          ''handler in OpenVMS'')')
CC              CALL WRITES(IOER,OP_STRING,ERROR,*9999)
CC
C              CALL LIB$REVERT() ! Reverts to old condition handler
C            ENDIF
C          ENDIF
C
C          IF(END) THEN       !The QUIT command has been used to quit
C            CONTINUE=.FALSE.
C          ENDIF
C          GOTO 200
C
CC ***     Handle error condition
C 150      CALL STRING_TRIM(ERROR,IBEG,IEND)
C          WRITE(OP_STRING,'(A)') ' '//ERROR(IBEG:IEND)//'>DIALOG'
C          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C          IF(CTRLC) THEN
C            CTRLC=.FALSE.
CC CPB 23/10/93 This is done when CTRLC_AST is called so no need to
CC do it here
CC            !if mode is not batch
CC            IF(STATSTR(1:1).ne.'B')THEN
CC              ISTATUS=SYS$QIOW(,%VAL(INPUT_CHAN),%VAL(ICODE),,,,
CC     '          CTRLC_AST,,,,,)
CC              IF(.NOT.ISTATUS) CALL LIB$SIGNAL(%VAL(ISTATUS))
CC            ENDIF
C          ENDIF
C          ERROR(1:)=' '
C
C 200      CONTINUE
C        ENDDO !end of main program loop
C
C      ELSE   !USE_SOCKET = .TRUE.
C
C        CONNID1=0
C        CONNID2=1
C
C        IF (FSKLISTEN(CONNID1,PORT1) .EQ. -1) GOTO 9998
C        IF (FSKLISTEN(CONNID2,PORT2) .EQ. -1) GOTO 9998
C
CC ***   Establish condition handler
C        FATAL_HANDLER=.TRUE.
C        CHANGE_HANDLER=.FALSE.
C        CALL LIB$ESTABLISH(NEW_HANDLER) !estab new condit handler
C
C        CONTINUE=.TRUE.
C        DO WHILE (CONTINUE)
CC CPB 14/3/94 Adding timeout for sockets to enable graphics updates
C          RETURNVAL=0
C          DO WHILE(RETURNVAL.EQ.0)
C            RETURNVAL=FSKSELECT(CONNID1,200) ! timeout after 200ms
C            IF(RETURNVAL.EQ.-1) GOTO 9998
CC GMH 15/11/95 making call conditional
C          IF(USE_GRAPHICS.EQ.1) THEN
C            CALL GXWAIT(0.0,ERR) !Update graphics
C          ENDIF
C    ENDDO
C          IF (FSKREAD(LEN,SK_LONG_INT,1,CONNID1) .EQ. -1) GOTO 9998
C          IF (FSKREAD(INTSTR,SK_CHAR,LEN+1,CONNID1) .EQ. -1) GOTO 9998
C          CALL FSKC2F(STRING,LEN,INTSTR)
C
C          CALL STRING_TRIM(STRING,IBEG,IEND)
C          WRITE(OP_STRING,'(A)') STRING(IBEG:IEND)
C          CALL WRITES(IOOP,OP_STRING,ERROR,*9998)
C          IF(STRING(IBEG:IBEG+3).ne.'QUIT') THEN
C            CALL PARSE(0,ISEG,CSEG,STRING,END,ERROR,*160)
C            WRITE(OP_STRING,'(/1X,A)') '> '
C            CALL WRITES(IOOP,OP_STRING,ERROR,*9998)
C          ENDIF
C
C          GOTO 210
C
CC ***   Handle error condition
C 160      CALL STRING_TRIM(ERROR,IBEG,IEND)
C          WRITE(OP_STRING,'(A)') ' '//ERROR(IBEG:IEND)//'>DIALOG'
C          CALL WRITES(IOER,OP_STRING,ERROR,*9998)
C          ERROR(1:)=' '
C
CC ***   Tell front end that CMISS is back in the main loop (used for
CC ***   popping down the prompt window)
C 210      IF(FSKWRITE(-1,SK_LONG_INT,1,CONNID2).EQ.-1) GOTO 9998
C
C          IF (STRING(1:4) .EQ. 'QUIT') THEN
CC ***   First tell front-end to close down sockets, then close down own
CC ***   sockets.
C            IF (FSKWRITE(-2,SK_LONG_INT,1,CONNID2).EQ.-1) GOTO 9998
C            IRET=FSKCLOSE(CONNID1)
C            IRET=FSKCLOSE(CONNID2)
C            CALL QUIT(END,ERROR,*9999)
C            CONTINUE=.FALSE.
C          ENDIF
C
CC ***     Check for changes to condition handler
C          IF(CHANGE_HANDLER) THEN
C            IF(FATAL_HANDLER) THEN
C              CALL LIB$ESTABLISH(NEW_HANDLER) !estab new condit handler
C            ELSE
C              CALL LIB$REVERT() ! Reverts to old condition handler
C            ENDIF
C          ENDIF
C
C        ENDDO
C
C      ENDIF
C
CC ??? MPN 27-oct-94 Is there a better place this?
C!news MPN 13-Sep-94: Parallel element stiffness matrix computations
C      IF(KTYP1A.EQ.2) THEN  !parallel element stiffness matrix calcs
C        DO nhost=1,NUMHOSTS_USED
C          IF(SOCKET_OPEN(nhost)) THEN
C            !Signal to slave processes to stop
C            IF(FSKWRITE(QUIT_PROCESS,SK_LONG_INT,1,ICONNID(nhost))
C     '        .EQ.-1) GOTO 9999
C            IRET=FSKCLOSE(ICONNID(nhost)) !close socket
C          ENDIF
C        ENDDO
C      ENDIF
C!newe
C
C      CALL EXITS('DIALOG')
C
C      RETURN
C 9999 CALL ERRORS('DIALOG',ERROR)
C      WRITE(OP_STRING,'(132A)') ' Fatal error propagated to DIALOG:'
C      CALL WRITES(IOER,OP_STRING,ERROR,*9998)
C      CALL STRING_TRIM(ERROR,IBEG,IEND)
C      WRITE(OP_STRING,'(132A)') ' '//ERROR(IBEG:IEND)//'>DIALOG'
C      CALL WRITES(IOER,OP_STRING,ERROR,*9998)
C      CLOSE(UNIT=IOER)
C      CALL EXITS('DIALOG')
C 9998 IF(USE_SOCKET) THEN
C        IRET=FSKCLOSE(CONNID1)
C        IRET=FSKCLOSE(CONNID2)
C      ENDIF
C      RETURN
C      END


      SUBROUTINE DO_COMMAND(COMMAND,ERROR,*)

C#### Subroutine: DO_COMMAND
C###  Description:
C###    DO_COMMAND executes an operating system command.

      IMPLICIT NONE
!     Parameter List
      CHARACTER COMMAND*(*),ERROR*(*)
!     Local Variables
      INTEGER IBEG,IEND,ISTATUS,LIB$SPAWN,NOWAIT

      CALL ENTERS('DO_COMMAND',*9999)

      NOWAIT=1 !don't wait for process to finish
      CALL STRING_TRIM(COMMAND,IBEG,IEND)
      ISTATUS=LIB$SPAWN(COMMAND(IBEG:IEND),
     '  'cmiss$reference:null_do_not_delete.',,NOWAIT,,,,,,,,,)

      CALL EXITS('DO_COMMAND')
      RETURN
 9999 CALL ERRORS('DO_COMMAND',ERROR)
      CALL EXITS('DO_COMMAND')
      RETURN 1
      END


      SUBROUTINE ENTERS(NAME,*)

C#### Subroutine: ENTERS
C###  Description:
C###    ENTERS traces entry to a subprogram recording the level of
C###    nesting, the entry time, and writing the subprogram name to a
C###    trace file.  Diagnostic o/p is turned on if DIAGNO=.TRUE. and
C###    ALLSUB=.TRUE. or ALLSUB=.FALSE. and NAME=SUBNAM (a subroutine
C###    name).
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
      INCLUDE '($IODEF)'
      INCLUDE '($JPIDEF)'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:ctrl00.cmn'
      INCLUDE 'cmiss$reference:diag00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
      INCLUDE 'cmiss$reference:trac00.cmn'
!     Parameter List
      CHARACTER NAME*(*)
!     Local Variables
      INTEGER IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,n1sb,nsub
      REAL*8 TM,VTIME
      CHARACTER COLV*10,C1*(MXCH),ERROR*10
C      EXTERNAL CTRLC_AST

      IF(CTRLC) THEN
        CTRLC=.FALSE.
C CPB 23/10/93 This is done when CTRLC_AST is called so no need to do it here
C        ISTATUS=SYS$QIOW(,%VAL(INPUT_CHAN),%VAL(ICODE),,,,CTRLC_AST,
C     '    ,,,,)
C        IF(.NOT.ISTATUS) CALL LIB$SIGNAL(%VAL(ISTATUS))
        GO TO 9999
      ENDIF

      IF(DIAGNO) THEN
        CALL STRING_TRIM(NAME,IBEG1,IEND1)
        IF(ALLSUB) THEN !turn diagnostics on in all subroutines
          WRITE(OP_STRING,'('' *** Enter '',A)') NAME(IBEG1:IEND1)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ELSE IF(.NOT.ALLSUB) THEN !diags on in selected subrs
          num_stack=num_stack+1
          DOP_STACK(num_stack)=.FALSE.
          DO nsub=1,NT_SUB
            CALL STRING_TRIM(SUBNAM(nsub),IBEG2,IEND2)
            CALL CUPPER(NAME(IBEG1:IEND1),C1)
C          IF(CUPPER(NAME(IBEG1:IEND1)).EQ.SUBNAM(nsub)(IBEG2:IEND2))
            IF(C1(IBEG1:IEND1).EQ.SUBNAM(nsub)(IBEG2:IEND2))
     '        DOP_STACK(num_stack)=.TRUE.
          ENDDO
          IF(FROMSUB) THEN
            IF(DOP_STACK(num_stack-1)) DOP_STACK(num_stack)=.TRUE.
          ENDIF
          IF(DOP_STACK(num_stack)) THEN
            WRITE(OP_STRING,'('' *** Enter '',A)') NAME(IBEG1:IEND1)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE IF(DOP_STACK(num_stack-1)) THEN
            WRITE(OP_STRING,'('' *** Calls '',A)') NAME(IBEG1:IEND1)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
          DOP=DOP_STACK(num_stack)
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
          DO n1sb=1,NTSB
            IF(SB(n1sb).EQ.NAME) THEN
              NOSB=n1sb
              GOTO 2
            ENDIF
          ENDDO
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

C#### Subroutine: EXECUTE_COMFILE
C###  Description:
C###    EXECUTE_COMFILE executes an operating system command file or
C###    shell script file.

      IMPLICIT NONE
!     Parameter List
      CHARACTER FILE_NAME*(*)
!     Local Variables
      INTEGER IBEG,IEND,ISTATUS,LIB$SPAWN

      CALL STRING_TRIM(FILE_NAME,IBEG,IEND)
      ISTATUS=LIB$SPAWN('@'//FILE_NAME(IBEG:IEND)//'.com',,,,,,,,,,,)

      RETURN
      END


      SUBROUTINE GET_CMISS_EXAMPLES(CMISS_EXAMPLES,ERROR,*)

C#### Subroutine: GET_CMISS_EXAMPLES
C###  Description:
C###    GET_CMISS_EXAMPLES returns CMISS examples directory as
C###    string in CMISS_EXAMPLES.

      IMPLICIT NONE
      INCLUDE '($DVIDEF)'
      INCLUDE '($IODEF)'
      INCLUDE '($LNMDEF)'
      INCLUDE '($SSDEF)'
!     Parameter List
      CHARACTER CMISS_EXAMPLES*255,ERROR*(*)
!     Local Variables
      INTEGER*4 CMISS_EXAMPLES_LEN,STATUS
      INTEGER SYS$TRNLNM
C Item list common block

      STRUCTURE /ITMLST/
        UNION
          MAP
            INTEGER*2 BUFLEN
            INTEGER*2 ITMCOD
            INTEGER*4 BUFADR
            INTEGER*4 RETADR
          END MAP
          MAP
            INTEGER*4 END_LIST
          END MAP
        END UNION
      END STRUCTURE

      RECORD /ITMLST/ ITEMLIST(2)

C GMH 11/12/96 Trace, etc uninitialised for this call
C      CALL ENTERS('GET_CMISS_EXAMPLES',*9999)

C Translate the CMISS$EXAMPLES logical

      ITEMLIST(1).ITMCOD=LNM$_STRING
      ITEMLIST(1).BUFLEN=LEN(CMISS_EXAMPLES)
      ITEMLIST(1).BUFADR=%LOC(CMISS_EXAMPLES)
      ITEMLIST(1).RETADR=%LOC(CMISS_EXAMPLES_LEN)
      ITEMLIST(2).ITMCOD=0
      ITEMLIST(2).BUFLEN=0

      STATUS=SYS$TRNLNM(,'LNM$JOB','CMISS$EXAMPLES',,ITEMLIST)
      IF(STATUS.EQ.SS$_NOLOGNAM) THEN
        STATUS=SYS$TRNLNM(,'LNM$PROCESS_TABLE','CMISS$EXAMPLES',,
     '    ITEMLIST)
      ENDIF

C GMH 11/12/96 Trace, etc uninitialised for this call
C      CALL EXITS('GET_CMISS_EXAMPLES')
C      RETURN
C 9999 CALL ERRORS('GET_CMISS_EXAMPLES',ERROR)
C      CALL EXITS('GET_CMISS_EXAMPLES')
C      RETURN 1
      RETURN
      END


      SUBROUTINE GET_COMMAND_LINE(ARGS,NUMARGS)

C#### Subroutine: GET_COMMAND_LINE
C###  Description:
C###    GET_COMMAND_LINE gets the command line arguments.

      IMPLICIT NONE
!     Parameter List
      INTEGER NUMARGS
      CHARACTER ARGS(10)*80
!     Local Variables
      INTEGER CLOCAT,i,IBEG,IEND,LENGTH,OLDPOSN,POSN
      CHARACTER COMMAND_LINE*80

      CALL LIB$GET_FOREIGN(COMMAND_LINE,,LENGTH,)
      CALL STRING_TRIM(COMMAND_LINE,IBEG,IEND)
      i=0
      OLDPOSN=IBEG
      POSN=CLOCAT(' ',COMMAND_LINE(OLDPOSN:))
      DO WHILE(posn.ne.1)
        i=i+1
        ARGS(I)=COMMAND_LINE(OLDPOSN:OLDPOSN+POSN-2)
        OLDPOSN=OLDPOSN+POSN
        POSN=CLOCAT(' ',COMMAND_LINE(OLDPOSN:))
        IF(i.EQ.10) POSN=1 !no more than 10 args at present
      ENDDO
      NUMARGS=I

      RETURN
      END


      SUBROUTINE GET_DATE_TIME(DATE,IDATEFMT,ERROR,*)

C#### Subroutine: GET_DATE_TIME
C###  Description:
C###    <HTML> <PRE>
C###    GET_DATE_TIME uses runtime library call to return date and time
C###    with IDATEFMT equal to:
C###    1  -  returns date in DEC format "DD-MMM-YYYY HH:MM:SS.HH"
C###          ie 3:45 in the afternoon on 27th of November 1990
C###          becomes "27-NOV-1990 15:45:00.00". (Default).
C###    2  -  returns date as all numbers in IGES format "YYMMDD.HHMMSS"
C###          ie the same date is returned as "901127.154500"
C###    3  -  returns a 7 or 8 digit random number based on the time
C###          as "HHSSMMHH" working from hundredths back to hours
C###    </PRE> </HTML>

      IMPLICIT NONE
!     Parameter List
      INTEGER IDATEFMT
      CHARACTER DATE*(*),ERROR*(*)
!     Local Variables
1      INTEGER MONTHNUM
      CHARACTER MONTH*(3),TEMPDATE*(30)

      CALL ENTERS('GET_DATE_TIME',*9999)

      CALL LIB$DATE_TIME(TEMPDATE)
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

C#### Subroutine: GET_DISPLAY
C###  Description:
C###    GET_DISPLAY returns the Xwindows display node (ie node on which
C###    graphics is displayed).  If user has defined 'display' as
C###    remote nodename with command
C###       setenv display nodename::0.0 (via Decnet)
C###    or setenv display nodename:0.0  (via TCP/IP)
C###  then NAME is remote display name,
C###  else NAME is local  display name.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:disp00.cmn'
      INCLUDE '($DVIDEF)'
      INCLUDE '($IODEF)'
      INCLUDE '($LNMDEF)'
      INCLUDE '($SSDEF)'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER DEVICE_CHANNEL,DISP_LEN
      INTEGER DECW$C_WS_DSP_NODE,IO$M_WS_DISPLAY
      INTEGER SYS$ASSIGN,SYS$QIO,SYS$TRNLNM
      INTEGER*4 STATUS,FUNC_CODE,DEV_STR_LEN
      CHARACTER DISPLAY_DEVICE*15,DISP_NAME*15
      LOGICAL UNKNOWN_DISP
C Item list common block
c     COMMON /ITLIST/ NAME_LEN,NAME_CODE,NAME_ADDR,RET_ADDR,END_LIST
C IO status common block
      INTEGER*2 IOSTAT
c     COMMON /IOBLOCK/ IOSTAT,MSG_LEN,READER_PID

      DATA DECW$C_WS_DSP_NODE/1/
      DATA IO$M_WS_DISPLAY/64/

      INTEGER*2 IWINDOW_SYSTEM_LEN
      INTEGER IWINDOW_SYSTEM,SYS$GETSYIW
      STRUCTURE /ITMLST/
        UNION
          MAP
            INTEGER*2 BUFLEN
            INTEGER*2 ITMCOD
            INTEGER*4 BUFADR
            INTEGER*4 RETADR
          END MAP
          MAP
            INTEGER*4 END_LIST
          END MAP
        END UNION
      END STRUCTURE
      EXTERNAL SYI$_WINDOW_SYSTEM

      RECORD /ITMLST/ ITEMLIST(2)

      CALL ENTERS('GET_DISPLAY',*9999)

      !Get display name

C CPB 8/3/94 Changing the display name calculation to pick up the
C display from the DECW$DISPLAY logical

C Translate the DECW$DISPLAY logical

      ITEMLIST(1).ITMCOD=LNM$_STRING
      ITEMLIST(1).BUFLEN=LEN(DISPLAY_DEVICE)
      ITEMLIST(1).BUFADR=%LOC(DISPLAY_DEVICE)
      ITEMLIST(1).RETADR=%LOC(DEV_STR_LEN)
      ITEMLIST(2).ITMCOD=0
      ITEMLIST(2).BUFLEN=0

      UNKNOWN_DISP=.FALSE.
      STATUS=SYS$TRNLNM(,'LNM$JOB','DECW$DISPLAY',,ITEMLIST)
      IF(STATUS.EQ.SS$_NOLOGNAM) THEN
        STATUS=SYS$TRNLNM(,'LNM$PROCESS_TABLE','DECW$DISPLAY',,ITEMLIST)
      ENDIF
      IF(STATUS.EQ.SS$_NOLOGNAM) THEN

C Initialise item list for the window_system

        ITEMLIST(1).ITMCOD=%LOC(SYI$_WINDOW_SYSTEM)
        ITEMLIST(1).BUFLEN=4 !bytes for window_system output
        ITEMLIST(1).BUFADR=%LOC(IWINDOW_SYSTEM)
        ITEMLIST(1).RETADR=%LOC(IWINDOW_SYSTEM_LEN)
        ITEMLIST(2).ITMCOD=0
        ITEMLIST(2).BUFLEN=0

        STATUS=SYS$GETSYIW(,,,ITEMLIST,IOSTAT,,)
        IF(IWINDOW_SYSTEM.EQ.0) THEN
          WINDOW_TYPE='TERMINAL'
        ELSE IF(IWINDOW_SYSTEM.EQ.1) THEN
          WINDOW_TYPE='MOTIF'
          UNKNOWN_DISP=.TRUE.  ! cannot translate DECW$DISPLAY and the
             ! window system is MOTIF hence error
        ELSE IF(IWINDOW_SYSTEM.EQ.2) THEN
          WINDOW_TYPE='VWS'
        ELSE
          ERROR='Unknown window type'
          GO TO 9999
        ENDIF

      ELSE
        WINDOW_TYPE='MOTIF'
      ENDIF

C Assign a channel to the WSA device

      STATUS=SYS$ASSIGN(DISPLAY_DEVICE,DEVICE_CHANNEL,,)
      IF(STATUS.EQ.SS$_NORMAL) THEN

C Get the node string

        FUNC_CODE=IO$_SENSEMODE.OR.IO$M_WS_DISPLAY
        DISP_LEN=15
        STATUS=SYS$QIO(,%VAL(DEVICE_CHANNEL),%VAL(FUNC_CODE),IOSTAT,,,
     '    %REF(DISP_NAME),%VAL(DISP_LEN),%VAL(DECW$C_WS_DSP_NODE),
     '    0,0,0)
      ENDIF

      DISPLAY_NAME=DISP_NAME
      IF(DISPLAY_NAME(1:1).EQ.'0') THEN
        DISPLAY_NAME=NODE_NAME
      ENDIF

      IF(DISPLAY_NAME(1:5).EQ.NODE_NAME(1:5)) THEN
        IF(UNKNOWN_DISP) THEN
          TYPE *,'Local display name  is UNKNOWN'
          TYPE *,'ERROR: Set DECW$DISPLAY with SET DISPLAY/CREATE',
     '      ' command'
        ELSE
          TYPE *,'Local display name  is ',DISPLAY_NAME
        ENDIF
      ELSE
        IF(UNKNOWN_DISP) THEN
          TYPE *,'Remote display name is UNKNOWN'
          TYPE *,'ERROR: Set DECW$DISPLAY with SET DISPLAY/CREATE/',
     '      'NODE=???? command'
        ELSE
          TYPE *,'Remote display name is ',DISPLAY_NAME
        ENDIF
      ENDIF
      TYPE *,'Window system       is ',WINDOW_TYPE

      CALL EXITS('GET_DISPLAY')
      RETURN
 9999 CALL ERRORS('GET_DISPLAY',ERROR)
      CALL EXITS('GET_DISPLAY')
      RETURN 1
      END


      SUBROUTINE GET_REV_TIME(ERROR,*)

C#### Subroutine: GET_REV_TIME
C###  Description:
C###    GET_REV_TIME returns the current revision time for the CMISS
C###    executable.

      IMPLICIT NONE

      INCLUDE '($JPIDEF)'
      INCLUDE '($SYSSRVNAM)'
      INCLUDE '($RMSDEF)'
      INCLUDE '($FABDEF)'
      INCLUDE '($XABDATDEF)'
      INCLUDE '($XABDEF)'
      INCLUDE 'cmiss$reference:disp00.cmn'

!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      RECORD /FABDEF/ MYFAB
      RECORD /XABDEF/ MYXAB
      CHARACTER*100 RESSTRING, FILENAME
      CHARACTER*50 REVTIMESTRING
      INTEGER*4 PROCESSID,STATUS
      INTEGER RESVALUE,RESLENGTH

      CALL ENTERS('GET_REV_TIME',*9999)

C
C Get the imagename of the current process i.e. the name of the current
C CMISS executable being run
C
      PROCESSID=0
      CALL LIB$GETJPI(JPI$_IMAGNAME,PROCESSID,,RESVALUE,
     '  RESSTRING,RESLENGTH)
      !store the name of the executable image (used in fe07)
      IMAGE_NAME=RESSTRING(1:RESLENGTH)
C
C Initialise the File Access Block (FAB) fields and set the filename
C to be the name of the CMISS executable
C
      MYFAB.FAB$B_BID=FAB$C_BID
      MYFAB.FAB$B_BLN=FAB$C_BLN
      MYFAB.FAB$L_FNA=%LOC(RESSTRING)
      MYFAB.FAB$B_FNS=RESLENGTH
C
C Initialise the eXtended Access Block (XAB) fields
C
      MYXAB.XAB$B_COD=XAB$C_DAT
      MYXAB.XAB$B_BLN=XAB$C_DATLEN
      MYFAB.FAB$L_XAB=%LOC(MYXAB)
C
C Use RMS to open the file and fill in the missing FAB and XAB fields
C
      STATUS=SYS$OPEN(MYFAB)
      IF(STATUS.EQ.RMS$_NORMAL) THEN
C
C Convert the 64-bit Revision time to an ASCII date string
C
        CALL SYS$ASCTIM(,REVTIMESTRING,MYXAB.XAB$Q_RDT,)
        TYPE *,'CMISS revision time ',REVTIMESTRING
      ENDIF

      CALL EXITS('GET_REV_TIME')
      RETURN
 9999 CALL ERRORS('GET_REV_TIME',ERROR)
      CALL EXITS('GET_REV_TIME')
      RETURN 1
      END


      SUBROUTINE GET_SYSTEM(ERROR,*)

C#### Subroutine: GET_SYSTEM
C###  Description:
C###    GET_SYSTEM returns current system nodename, architecture and
C###    operating system.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:disp00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER*2 ARCHNAME_LEN,HWNAME_LEN,NODENAME_LEN,VERSION_LEN
      INTEGER SYS$GETSYIW,STATUS
      CHARACTER ARCHNAME*15,HWNAME*50,NODENAME*15,VERSION*8
      STRUCTURE /ITMLST/
        UNION
          MAP
            INTEGER*2 BUFLEN
            INTEGER*2 ITMCOD
            INTEGER*4 BUFADR
            INTEGER*4 RETADR
          END MAP
          MAP
            INTEGER*4 END_LIST
          END MAP
        END UNION
      END STRUCTURE
      STRUCTURE /IOSBLK/
        INTEGER*4 STS,RESERVED
      END STRUCTURE
      EXTERNAL SYI$_ARCH_NAME,SYI$_HW_NAME,SYI$_NODENAME,SYI$_VERSION

      RECORD /ITMLST/ GETSYI_LIST(5)
      RECORD /IOSBLK/ IOSB

      CALL ENTERS('GET_SYSTEM',*9999)


C CPB 12/10/94 Adding architecture type and VMS version number

      !Initialise item list for node_name,architecture and OS version
      GETSYI_LIST(1).BUFLEN   =15 !bytes for node_name output
      GETSYI_LIST(1).ITMCOD   =%LOC(SYI$_NODENAME)
      GETSYI_LIST(1).BUFADR   =%LOC(NODENAME)
      GETSYI_LIST(1).RETADR   =%LOC(NODENAME_LEN)
      GETSYI_LIST(2).BUFLEN   =15 !bytes for architecture output
      GETSYI_LIST(2).ITMCOD   =%LOC(SYI$_ARCH_NAME)
      GETSYI_LIST(2).BUFADR   =%LOC(ARCHNAME)
      GETSYI_LIST(2).RETADR   =%LOC(ARCHNAME_LEN)
      GETSYI_LIST(3).BUFLEN   =50 !bytes for architecture output
      GETSYI_LIST(3).ITMCOD   =%LOC(SYI$_HW_NAME)
      GETSYI_LIST(3).BUFADR   =%LOC(HWNAME)
      GETSYI_LIST(3).RETADR   =%LOC(HWNAME_LEN)
      GETSYI_LIST(4).BUFLEN   =8 !bytes for version output
      GETSYI_LIST(4).ITMCOD   =%LOC(SYI$_VERSION)
      GETSYI_LIST(4).BUFADR   =%LOC(VERSION)
      GETSYI_LIST(4).RETADR   =%LOC(VERSION_LEN)
      GETSYI_LIST(5).END_LIST =0 !terminator
      STATUS=SYS$GETSYIW(,,,GETSYI_LIST,IOSB,,)
      NODE_NAME=NODENAME
      ARCH_TYPE=ARCHNAME
      OS_TYPE='VMS'
      OS_VERSION='VMS '//VERSION

      TYPE *,'System nodename     is ',NODE_NAME(1:NODENAME_LEN)
      TYPE *,'System architecture is ',ARCH_TYPE(1:ARCHNAME_LEN)
      TYPE *,'System type         is ',HWNAME(1:HWNAME_LEN)
      TYPE *,'Operating system    is ',OS_TYPE(1:VERSION_LEN+4)

      CALL EXITS('GET_SYSTEM')
      RETURN
 9999 CALL ERRORS('GET_SYSTEM',ERROR)
      CALL EXITS('GET_SYSTEM')
      RETURN 1
      END


      SUBROUTINE GETSTR1(STRING,PR,LNPR,ENDCOM,ERROR,*)

C#### Subroutine: GETSTR1
C###  Description:
C###    GETSTR1 inserts a prompt, PROMPT, in the input file.  It
C###    returns the string typed by the user.  If LEARN is true the
C###    string is written to FILE00.  The F10 key indicates an eof -
C###    END=100 in read sets ENDCOM=.TRUE..  PROMPT is prompt for SMG
C###    routine.  Pressing the enter key causes a sequential read from
C###    unit 76 (com file).

C**** MACRO_KEY_buffer(nomacr,nokey),nomacr=1,NT_MACRO(nokey)
C**** are command lines associated with macro key number NOKEY.
C**** MACRO_DEFINE is .true. when macro is being defined ('Do' key to
C**** start and finish).
C**** MACRO_KEY_EXECUTE is .true. when macro is executing.
C**** NOTE: $ format is NONSTANDARD FORTRAN 77.

      IMPLICIT NONE
      INCLUDE '($SMGDEF)'
      INCLUDE '($TRMDEF)'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:gtstr00.cmn'
      INCLUDE 'cmiss$reference:learn00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER LNPR
      CHARACTER ERROR*(*),PR*(*),STRING*(MXCH)
      LOGICAL ENDCOM
!     Local Variables
      INTEGER*2 TERMINATOR,TERMINATOR2
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,IOSTAT,
     '  IREC_COMMENT,ISTATUS,
     '  KEYBOARD,LN_PROMPT,LN_PROMPT_OLD,LN_STRING,nomacro,
     '  SMG$CREATE_VIRTUAL_KEYBOARD,SMG$READ_KEYSTROKE,
     '  UP_COUNT
      CHARACTER CHAR*1,CHAR1*1,COMMAND*80,COMMENT*80,
     '  C1*(MXCH),LAST_COMMENT*80,KEY_BUFFER(4)*132,PROMPT*132
      LOGICAL COMMENTS,LEFT_ARROW

      CALL ENTERS('GETSTR1',*9999)
      WRITE(OP_STRING,*)    !gives blank line for SMG routine
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ISTATUS=SMG$CREATE_VIRTUAL_KEYBOARD(KEYBOARD,'SYS$INPUT')
      IF(.NOT.ISTATUS) CALL LIB$SIGNAL(%VAL(ISTATUS))

      UP_COUNT=0
      LEFT_ARROW=.FALSE.
      PROMPT=PR
      LN_PROMPT=LNPR
      STRING=PROMPT(1:LN_PROMPT)
      LN_STRING=LNPR

      TERMINATOR=0
      DO WHILE(terminator.ne.13)
        ISTATUS=SMG$READ_KEYSTROKE(KEYBOARD,TERMINATOR,
     '    PROMPT(2:LN_PROMPT),1,)
        PROMPT=PROMPT(1:LN_PROMPT)//CHAR(TERMINATOR)

        IF(TERMINATOR.EQ.SMG$K_TRM_UP) THEN !Up_arrow
C CPB 28/4/94 Deleting these lines because they corrupt the type-ahead
C terminal buffer
C          ISTATUS=SMG$READ_KEYSTROKE(KEYBOARD,TERMINATOR2,
C     '      BLANK(2:LN_PROMPT),0,) !blanks out current line & times out
          IF(UP_COUNT.GE.0.AND.UP_COUNT.LT.BUFFER_COUNT) THEN
            UP_COUNT=UP_COUNT+1
            PROMPT=BUFFER(UP_COUNT)
            LN_PROMPT=BUFFER_LENGTH(UP_COUNT)
          ELSE
            PROMPT=PR
            LN_PROMPT=LNPR
          ENDIF
          STRING=PROMPT(1:LN_PROMPT)
          LN_STRING=LN_PROMPT

        ELSE IF(TERMINATOR.EQ.SMG$K_TRM_DOWN) THEN !Down_arrow
C CPB 28/4/94 Deleting these lines because they corrupt the type-ahead
C terminal buffer
C          ISTATUS=SMG$READ_KEYSTROKE(KEYBOARD,TERMINATOR2,
C     '      BLANK(2:LN_PROMPT),0,) !blanks out current line & times out
          IF(UP_COUNT.GT.0.AND.UP_COUNT.LE.BUFFER_COUNT) THEN
            PROMPT=BUFFER(UP_COUNT)
            LN_PROMPT=BUFFER_LENGTH(UP_COUNT)
            UP_COUNT=UP_COUNT-1
          ELSE
            PROMPT=PR
            LN_PROMPT=LNPR
          ENDIF
          STRING=PROMPT(1:LN_PROMPT)
          LN_STRING=LN_PROMPT

        ELSE IF(TERMINATOR.EQ.SMG$K_TRM_LEFT) THEN !Left_arrow
          IF(.NOT.LEFT_ARROW) THEN
            LEFT_ARROW=.TRUE.
            STRING=PROMPT(1:LN_PROMPT)
            LN_STRING=LN_PROMPT
          ENDIF
          LN_PROMPT=LN_PROMPT-1
          IF(LN_PROMPT.LT.2) LN_PROMPT=2

        ELSE IF(TERMINATOR.EQ.SMG$K_TRM_RIGHT) THEN !Right_arrow
          LN_PROMPT=LN_PROMPT+1
          IF(LN_PROMPT.GT.132) LN_PROMPT=132
          PROMPT=STRING(1:LN_PROMPT)

        ELSE IF(TERMINATOR.EQ.127) THEN !Delete a character
C CPB 28/4/94 Deleting these lines because they corrupt the type-ahead
C terminal buffer
C          ISTATUS=SMG$READ_KEYSTROKE(KEYBOARD,TERMINATOR2,
C     '      BLANK(2:LN_STRING),0,) !blanks out current line & times out
          LN_PROMPT=LN_PROMPT-1
          IF(LN_PROMPT.LT.2) LN_PROMPT=2
          IF(LN_STRING.LE.LN_PROMPT+1) THEN
            STRING=STRING(1:LN_PROMPT)
          ELSE IF(LN_STRING.GT.LN_PROMPT+1) THEN
            STRING=STRING(1:LN_PROMPT)//STRING(LN_PROMPT+2:LN_STRING)
          ENDIF
          ISTATUS=SMG$READ_KEYSTROKE(KEYBOARD,TERMINATOR2,
     '      STRING(2:LN_STRING),0,)
          LN_STRING=LN_STRING-1
          IF(LN_STRING.LT.2) LN_STRING=2

        ELSE IF(TERMINATOR.EQ.256) THEN !PF1 key
C CPB 28/4/94 Deleting these lines because they corrupt the type-ahead
C terminal buffer
          ISTATUS=SMG$READ_KEYSTROKE(KEYBOARD,TERMINATOR2,
     '      BLANK(2:LN_PROMPT),0,) !blanks out current line & times out
          PROMPT=KEY_BUFFER(1)(1:132)
          CALL STRING_TRIM(PROMPT,IBEG,IEND)
          STRING=PROMPT(1:IEND)
          LN_PROMPT=IEND

        ELSE IF(TERMINATOR.EQ.257) THEN !PF2 key
C CPB 28/4/94 Deleting these lines because they corrupt the type-ahead
C terminal buffer
C          ISTATUS=SMG$READ_KEYSTROKE(KEYBOARD,TERMINATOR2,
C     '      BLANK(2:LN_PROMPT),0,) !blanks out current line & times out
          PROMPT=KEY_BUFFER(2)(1:132)
          CALL STRING_TRIM(PROMPT,IBEG,IEND)
          STRING=PROMPT(1:IEND)
          LN_PROMPT=IEND

        ELSE IF(TERMINATOR.EQ.258) THEN !PF3 key
C CPB 28/4/94 Deleting these lines because they corrupt the type-ahead
C terminal buffer
C          ISTATUS=SMG$READ_KEYSTROKE(KEYBOARD,TERMINATOR2,
C     '      BLANK(2:LN_PROMPT),0,) !blanks out current line & times out
          PROMPT=KEY_BUFFER(3)(1:132)
          CALL STRING_TRIM(PROMPT,IBEG,IEND)
          STRING=PROMPT(1:IEND)
          LN_PROMPT=IEND

        ELSE IF(TERMINATOR.EQ.259) THEN !PF4 key
C CPB 28/4/94 Deleting these lines because they corrupt the type-ahead
C terminal buffer
C          ISTATUS=SMG$READ_KEYSTROKE(KEYBOARD,TERMINATOR2,
C     '      BLANK(2:LN_PROMPT),0,) !blanks out current line & times out
          PROMPT=KEY_BUFFER(4)(1:132)
          CALL STRING_TRIM(PROMPT,IBEG,IEND)
          STRING=PROMPT(1:IEND)
          LN_PROMPT=IEND

        ELSE IF(TERMINATOR.EQ.260) THEN !RH 0 key
          WRITE(OP_STRING,*) '>>RH 0 key not defined'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

        ELSE IF(TERMINATOR.GE.261.AND.TERMINATOR.LE.269) THEN
          MACRO_KEY_EXECUTE=.TRUE.
          MACRO_KEY=TERMINATOR-260
          WRITE(CHAR1,'(I1)') MACRO_KEY
C PJH 23Jan96
C         STRING(4:)='READ MACRO KEY '//CHAR1
          STRING(4:)=MACRO_names(MACRO_KEY)
          IF(DOP) THEN
            WRITE(OP_STRING,'('' >>Key '',I1,'' pressed'')') MACRO_KEY
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' >>Macro is '',A)')
     '        MACRO_names(MACRO_KEY)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CALL STRING_TRIM(STRING,IBEG,IEND)
          LN_STRING=IEND
          TERMINATOR=13

        ELSE IF(TERMINATOR.EQ.270) THEN !Enter key
          WRITE(OP_STRING,*) '>>Enter key not defined'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

        ELSE IF(TERMINATOR.EQ.295) THEN !Help key
          CALL STRING_TRIM(STRING,IBEG,IEND)
          STRING=STRING(1:IEND)//' ?'
          LN_STRING=IEND+2
          TERMINATOR=13

        ELSE IF(TERMINATOR.EQ.296) THEN !Do key
          IF(.NOT.MACRO_DEFINE) THEN
            WRITE(OP_STRING,'('' >>Use Insert key to put '
     '        //' commands into macro'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            MACRO_DEFINE=.TRUE.
            LN_PROMPT_OLD=LN_PROMPT
          ELSE IF(MACRO_DEFINE) THEN
            WRITE(OP_STRING,'('' >>Press key to define macro '
     '        //'(1..9 on RH keypad)'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
 20         ISTATUS=SMG$READ_KEYSTROKE(KEYBOARD,TERMINATOR2,
     '        PROMPT(2:LN_PROMPT),,)
            IF(TERMINATOR2.GE.261.AND.TERMINATOR2.LE.269) THEN
              MACRO_KEY=TERMINATOR2-260
            ELSE
              WRITE(OP_STRING,'('' >>Invalid key: use 1..9 on RH '','
     '          //'''keypad'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              GO TO 20
            ENDIF
            WRITE(OP_STRING,'('' >>Macro defined as key '',I1)')
     '        MACRO_KEY
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            NT_MACRO(MACRO_KEY)=NT_MACRO(0)   !defines #lines in macro
            DO nomacro=1,NT_MACRO(MACRO_KEY) !defines lines in macro
              MACRO_KEY_buffer(nomacro,MACRO_KEY)
     '          =MACRO_KEY_buffer(nomacro,0)
            ENDDO
            MACRO_DEFINE=.FALSE.
            NT_MACRO(0)=0
            LN_PROMPT=LN_PROMPT_OLD
            STRING=PROMPT(1:LN_PROMPT)
            LN_STRING=LN_PROMPT
          ENDIF

        ELSE IF(TERMINATOR.EQ.311) THEN !Find key
          WRITE(OP_STRING,*) '>>Find key not defined'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

        ELSE IF(TERMINATOR.EQ.312) THEN !Insert key
          NT_MACRO(0)=NT_MACRO(0)+1
          MACRO_KEY_buffer(NT_MACRO(0),0)=PROMPT(3:LN_PROMPT)
          WRITE(OP_STRING,'('' Macro line '',I2)') NT_MACRO(0)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C CPB 28/4/94 Deleting these lines because they corrupt the type-ahead
C terminal buffer
C          ISTATUS=SMG$READ_KEYSTROKE(KEYBOARD,TERMINATOR2,
C     '      BLANK(2:LN_STRING),0,) !blanks out current line & times out
          LN_PROMPT=LN_PROMPT_OLD

        ELSE IF(TERMINATOR.EQ.313) THEN !Remove key
          WRITE(OP_STRING,*) '>>Remove key not defined'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

        ELSE IF(TERMINATOR.EQ.314) THEN !Select key
          NT_KEY=MOD(NT_KEY,4)+1
          WRITE(CHAR1,'(I1)') NT_KEY
          KEY_BUFFER(NT_KEY)=PROMPT(1:LN_PROMPT)
          WRITE(OP_STRING,*) ' PF'//CHAR1(1:1)//' key defined'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

        ELSE IF(TERMINATOR.EQ.315) THEN !Previous screen key
C CPB 28/4/94 Deleting these lines because they corrupt the type-ahead
C terminal buffer
C          ISTATUS=SMG$READ_KEYSTROKE(KEYBOARD,TERMINATOR2,
C     '      BLANK(2:LN_PROMPT),0,) !blanks out current line & times out
          IREC_COMFILE=IREC_COMFILE-1
          PROMPT(3:)=' '
          READ(76,FMT='(A)',REC=IREC_COMFILE,IOSTAT=IOSTAT)
     '      PROMPT(4:100)
          IF(IOSTAT.EQ.0) THEN
            CALL STRING_TRIM(PROMPT,IBEG,IEND)
            STRING=PROMPT(1:IEND)
            LN_PROMPT=IEND
          ELSE IF(IOSTAT.EQ.25.OR.IOSTAT.EQ.39) THEN
            WRITE(OP_STRING,'(''  Top of file'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            LN_PROMPT=2
          ELSE
            WRITE(OP_STRING,*) 'iostat=',iostat
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
          LN_STRING=LN_PROMPT

        ELSE IF(TERMINATOR.EQ.316) THEN !Next screen key
C CPB 28/4/94 Deleting these lines because they corrupt the type-ahead
C terminal buffer
C          ISTATUS=SMG$READ_KEYSTROKE(KEYBOARD,TERMINATOR2,
C     '      BLANK(2:LN_PROMPT),0,) !blanks out current line & times out
          IREC_COMFILE=IREC_COMFILE+1
          PROMPT(3:)=' '

          IF(DOP) write(*,'('' IREC_COMFILE='',I4)') IREC_COMFILE
          READ(76,FMT='(A)',REC=IREC_COMFILE,IOSTAT=IOSTAT) COMMAND
          IF(DOP) write(*,'('' iostat='',I4,'' command='',A)')
     '            iostat,command

          IF(IOSTAT.EQ.0) THEN
            COMMENTS=.TRUE.
            IREC_COMMENT=0
            DO WHILE (COMMENTS)
              IREC_COMMENT=IREC_COMMENT+1
              READ(76,FMT='(A)',REC=IREC_COMFILE+IREC_COMMENT,
     '          IOSTAT=IOSTAT) COMMENT
              IF(COMMENT(1:1).EQ.'#') THEN !write comment line
                COMMENTS=.TRUE.
                IF(IREC_COMMENT.EQ.1) THEN
                  WRITE(OP_STRING,'(A)') PROMPT(1:3)//COMMAND
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDIF
                WRITE(OP_STRING,'(4X,A)') COMMENT(2:)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                LAST_COMMENT=' '//COMMENT(2:)
                COMMENT=' '
              ELSE
                COMMENTS=.FALSE.
              ENDIF
            ENDDO
            IF(IREC_COMMENT.GT.1) THEN !additional comment lines found
              IREC_COMFILE=IREC_COMFILE+IREC_COMMENT-1
              PROMPT(4:100)=LAST_COMMENT
            ELSE                       !no addditional lines
              PROMPT(4:100)=COMMAND
            ENDIF
            CALL STRING_TRIM(PROMPT,IBEG,IEND)
            LN_PROMPT=IEND
            CALL STRING_TRIM(COMMAND,IBEG,IEND)
            STRING=PROMPT(1:3)//COMMAND(1:IEND)
          ELSE IF(IOSTAT.EQ.36) THEN
            WRITE(OP_STRING,'(''  End of file'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            LN_PROMPT=2
          ENDIF
          LN_STRING=LN_PROMPT

        ELSE IF(TERMINATOR.EQ.SMG$K_TRM_CTRLC) THEN ! (ctrl C)
          GO TO 9999

        ELSE IF(TERMINATOR.EQ.290) THEN !F10 key (ctrl Z)
          GO TO 9998

        ELSE IF(TERMINATOR.EQ.13) THEN  !return key
C CPB 28/4/94 Deleting these lines because they corrupt the type-ahead
C terminal buffer
C          ISTATUS=SMG$READ_KEYSTROKE(KEYBOARD,TERMINATOR2,
C     '      BLANK(2:132),0,) !blanks out current line & times out

C CPB 22/10/93 Adding time out escape to enable correct type ahead and
C graphics updating
        ELSE IF(TERMINATOR.EQ.SMG$K_TRM_TIMEOUT) THEN  !time out
C GMH 15/11/95 making call conditional
          IF(USE_GRAPHICS.EQ.1) THEN
            CALL GXWAIT(0.0,ERR) !Update graphics
          ENDIF
        ELSE  !add a character
C CPB 22/10/93 Deleting these lines because they corrupt the type-ahead
C terminal buffer
C          ISTATUS=SMG$READ_KEYSTROKE(KEYBOARD,TERMINATOR2,
C     '      BLANK(2:LN_PROMPT),0,) !blanks out current line & times out
          LN_PROMPT=LN_PROMPT+1
          IF(LN_STRING.LE.LN_PROMPT) THEN
            STRING=PROMPT(1:LN_PROMPT)
            LN_STRING=LN_PROMPT
          ELSE IF(LN_STRING.GT.LN_PROMPT) THEN
            STRING=PROMPT(1:LN_PROMPT)//STRING(LN_PROMPT+1:LN_STRING)
          ENDIF

        ENDIF
      ENDDO
      WRITE(OP_STRING,*) STRING(2:LN_STRING)
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      IF(LEARN) THEN
        CALL STRING_TRIM(PR,IBEG1,IEND1)
        CALL STRING_TRIM(STRING,IBEG2,IEND2)
        CALL CUPPER(STRING(IBEG2:IBEG2+1),C1)
C        IF(CUPPER(STRING(IBEG2:IBEG2+1)).ne.'LE') THEN
        IF(C1.ne.'LE') THEN
          WRITE(76,'(A)') PR(IBEG1+2:IEND1+1)//STRING(IBEG2:IEND2)
        ENDIF
      ENDIF
C**** CALL CSTACK('PUSH',STRING,ERROR)

      CALL EXITS('GETSTR1')
      RETURN
 9998 ENDCOM=.TRUE.      !If CTRL-Z used to exit
      RETURN

 9999 CALL ERRORS('GETSTR1',ERROR)
      CALL EXITS('GETSTR1')
      RETURN 1
      END


      SUBROUTINE GETSTR2(ERROR,*)

C#### Subroutine: GETSTR2
C###  Description:
C###    GETSTR2 returns interrupt key number from mapped keyboard.
C###    Called from MARCH1 in module FE07.

      IMPLICIT NONE
      INCLUDE '($SMGDEF)'
      INCLUDE '($TRMDEF)'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:gtstr00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Parameter List
      INTEGER*2 TERMINATOR
      INTEGER ISTATUS,KEYBOARD,SMG$CREATE_VIRTUAL_KEYBOARD,
     '  SMG$READ_KEYSTROKE
      CHARACTER PROMPT*1
      DATA PROMPT/' '/
      SAVE KEYBOARD

      CALL ENTERS('GETSTR2',*9999)

      IF(KEYBOARD.EQ.0) THEN
        ISTATUS=SMG$CREATE_VIRTUAL_KEYBOARD(KEYBOARD,'SYS$INPUT')
        IF(.NOT.ISTATUS) CALL LIB$SIGNAL(%VAL(ISTATUS))
      ENDIF

      ISTATUS=SMG$READ_KEYSTROKE(KEYBOARD,TERMINATOR,PROMPT(1:1),0,)
      IF(TERMINATOR.GE.261.AND.TERMINATOR.LE.269) THEN
        MACRO_KEY_EXECUTE=.TRUE.
        MACRO_KEY=TERMINATOR-260
        WRITE(OP_STRING,'('' >>Key '',I1,'' pressed'')') MACRO_KEY
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        GO TO 9999
      ELSE
        MACRO_KEY_EXECUTE=.FALSE.
      ENDIF

      CALL EXITS('GETSTR2')
      RETURN

 9999 CALL ERRORS('GETSTR2',ERROR)
      CALL EXITS('GETSTR2')
      RETURN 1
      END


      SUBROUTINE GETSTR3(KEY,ERROR,*)

C#### Subroutine: GETSTR3
C###  Description:
C###    <HTML> <PRE>
C###    GETSTR3 returns the key typed by the user for DOCUM window.
C###    Prev_Screen  =  1
C###    Next_Screen  =  2
C###    Do           =  3
C###    Select       =  4
C###    Insert       =  5
C###    Delete       =  6
C###    Remove       =  7
C###    Alphanumeric keys return their ASCII character code.
C###    </PRE> </HTML>

      IMPLICIT NONE
      INCLUDE '($SMGDEF)'
      INCLUDE '($TRMDEF)'
!     Parameter List
      INTEGER KEY
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER*2 TERMINATOR
      INTEGER ISTATUS,KEYBOARD,SMG$CREATE_VIRTUAL_KEYBOARD,
     '  SMG$READ_KEYSTROKE
      LOGICAL CONTINUE

      CALL ENTERS('GETSTR3',*9999)

      ISTATUS=SMG$CREATE_VIRTUAL_KEYBOARD(KEYBOARD,'SYS$INPUT')
      IF(.NOT.ISTATUS) CALL LIB$SIGNAL(%VAL(ISTATUS))

      TERMINATOR=0
      CONTINUE=.TRUE.
      DO WHILE(CONTINUE)
        CONTINUE=.FALSE.
        ISTATUS=SMG$READ_KEYSTROKE(KEYBOARD,TERMINATOR,,,)
        IF(TERMINATOR.EQ.SMG$K_TRM_PREV_SCREEN) THEN
          KEY=1
        ELSE IF(TERMINATOR.EQ.SMG$K_TRM_NEXT_SCREEN) THEN
          KEY=2
        ELSE IF(TERMINATOR.EQ.SMG$K_TRM_DO) THEN
          KEY=3
        ELSE IF(TERMINATOR.EQ.SMG$K_TRM_SELECT) THEN
          KEY=4
C CPB 21/10/93 Wrong SMG constant name
C        ELSE IF(TERMINATOR.EQ.SMG$K_TRM_INSERT) THEN    !CHECK THIS ONE
        ELSE IF(TERMINATOR.EQ.SMG$K_TRM_INSERT_HERE) THEN
          KEY=5
        ELSE IF(TERMINATOR.EQ.SMG$K_TRM_DELETE) THEN
          KEY=6
        ELSE IF(TERMINATOR.EQ.SMG$K_TRM_REMOVE) THEN
          KEY=7
        ELSE IF(TERMINATOR.GE.32.AND.TERMINATOR.LE.127) THEN
          KEY=TERMINATOR
        ELSE
          CONTINUE=.TRUE.
        ENDIF
      ENDDO

      CALL EXITS('GETSTR3')
      RETURN
 9999 CALL ERRORS('GETSTR3',ERROR)
      CALL EXITS('GETSTR3')
      RETURN 1
      END


      SUBROUTINE MAINCMLOOP(EXAMPLE,EXNUM,COMFNAME,PARAMETERS,
     '  PARAMETERSFILENAME,FIRST)

C#### Subroutine: MAINCMLOOP
C###  Description:
C###    Main CMISS Loop

      IMPLICIT NONE
      INCLUDE '($IODEF)'
      INCLUDE '($JPIDEF)'
      INCLUDE '($SMGDEF)'
      INCLUDE '($TRMDEF)'
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
      INCLUDE 'cmiss$reference:cmgui00.cmn'
      INCLUDE 'cmiss$reference:cmis00.cmn'
      INCLUDE 'cmiss$reference:ctrl00.cmn'
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
      INCLUDE 'cmiss$reference:iter00.cmn'
      INCLUDE 'cmiss$reference:ktyp00.cmn'
      INCLUDE 'cmiss$reference:mach00.cmn'
      INCLUDE 'cmiss$reference:mach00.inc'
      INCLUDE 'cmiss$reference:map000.cmn'
      INCLUDE 'cmiss$reference:ntsg00.cmn'
      INCLUDE 'cmiss$reference:ptr00.cmn'
      INCLUDE 'cmiss$reference:trac00.cmn'
!     Parameter List
      INTEGER EXAMPLE,EXNUM(*),COMFNAME(*),PARAMETERS,
     '  PARAMETERSFILENAME(*),FIRST
!     Local Variables
      INTEGER MXCO,MXCOQU,MXSG
      PARAMETER (MXCO=25,MXCOQU=25,MXSG=10000)
      INTEGER CODE,CSTRLEN,CURRENT_TIME,ERR,IBEG,IBEG1,IBEG4,
     '  IEND,IEND1,IEND4,INTSTR(MXCH),
     '  IRET,ISEG(MXSG),ISTATUS,LEN,
     '  LIB$ESTABLISH,LIB$GETJPI,
     '  NEW_HANDLER,nhost,
     '  nomacro,
     '  nosg,PRIMARY,RETURNVAL,SECONDARY,SLEEP_DELAY
      LOGICAL CONTINUE,DEFINE,END,UPDATE
      CHARACTER CHAR4*4,COD(MXCO)*90,
     '  COMFILENAME*255,COQUD(MXCO,MXCOQU)*30,
     '  CSEG(MXSG)*60,EXAMPLENUM*255,ERROR*(MXCH),STATSTR*80,
     '  STRG*(MXCH),STRING*(MXCH)

      EXTERNAL CTRLC_AST,NEW_HANDLER

      CONTINUE=.TRUE.
      END=.FALSE.

      IF(.NOT.USE_SOCKET) THEN
        IF(CMGUI_LINK) THEN
C GMH 12/2/97 Move condition handler from the 'first' block
        OLD_HANDLER=LIB$ESTABLISH(NEW_HANDLER) !address of deflt handler

        IF(FIRST.EQ.0) THEN
          FIRST=1
          DO NOSG=1,MXSG
            ISEG(NOSG)=0
            CSEG(NOSG)=' '
          ENDDO
          COD(2)='PI'
          WRITE(COD(3),'(E11.5)') PI
          CALL ASSIGN(STRING,1,COD,COQUD,ERROR,*9999)

          FATAL_HANDLER=.TRUE.
          CHANGE_HANDLER=.FALSE.

          USEPARAMFILE=PARAMETERS.EQ.1
          IF(USEPARAMFILE) THEN
            CALL CSTRINGLEN(CSTRLEN,PARAMETERSFILENAME)
            CALL C2FSTRING(PARAMETERSFILENAME,CSTRLEN,PARAMETERSFILE)
            PARAMETERSFILE(CSTRLEN+1:)=' '
            CALL STRING_TRIM(PARAMETERSFILE,IBEG,IEND)
            STRG='fem define parameters;r;'//PARAMETERSFILE(IBEG:IEND)
C           temporarily set FIRST params to .FALSE. so that arrays won't
C           be allocated in FEM until after params read in.
            FIRST_SYNTAX=.FALSE.
            CALL PARSE(0,ISEG,CSEG,STRG,END,ERROR,*9999)
            FIRST_SYNTAX=.TRUE.
          ENDIF
        ENDIF !FIRST

        CONTINUE=.TRUE.
        DO WHILE(CONTINUE) !is main program loop
C         Perform any periodic tasks
          CALL PERIODICTASK
C         Check for available messages (rest in SYNTAX)
          CALL WH_INPUT_F_UPDATE(CMGUI_PROMPT_I,CODE)
          IF(CODE.EQ.0) THEN
            ERROR='Could not update prompt input'
            GOTO 9999
          ENDIF
          CALL WH_INPUT_F_UPDATE(CMGUI_COMMAND_I,CODE)
          IF(CODE.EQ.0) THEN
            ERROR='Could not update command input'
            GOTO 9999
          ENDIF
          CALL WH_OUTPUT_F_UPDATE(CMGUI_COMMAND_O,CODE)
          IF(CODE.EQ.0) THEN
            ERROR='Could not update command output'
            GOTO 9999
          ENDIF
          CALL WH_OUTPUT_F_CAN_OPEN(CMGUI_COMMAND_O,
     '      CODE)
          IF(CODE.NE.0) THEN
            CMGUI_MESSAGE_TIME=0
C           We have a command to process
            CALL WH_OUTPUT_F_OPEN_MESSAGE(CMGUI_COMMAND_O,
     '        PRIMARY,SECONDARY,CODE)
            IF(CODE.EQ.0) THEN
              ERROR='Could not open command message'
              GOTO 9999
            ENDIF
            IF(SECONDARY.EQ.1) THEN !idle message - do nothing
            ELSEIF(SECONDARY.EQ.2) THEN !command
              CALL WH_OUTPUT_F_NUM_ITEMS(CMGUI_COMMAND_O,
     '          LEN)
              IF(LEN.LE.MXCH) THEN
C               Zero the string
                STRING=' '
                CALL WH_OUTPUT_F_GET_CHAR(CMGUI_COMMAND_O,
     '          LEN,STRING,CODE)
                IF(CODE.EQ.0) THEN
                  ERROR='Could not get command'
                  GOTO 9999
                ENDIF
C               Send acknowledgement of command receipt
C               CALL WH_INPUT_F_OPEN_MESSAGE(CMGUI_COMMAND_O,
C    '            0,4,CODE)
C               IF(CODE.EQ.0) THEN
C                 ERROR='Could not open command receipt message'
C                 GOTO 9999
C               ENDIF
C               CALL WH_INPUT_F_CLOSE_MESSAGE(CMGUI_COMMAND_O,
C    '            CODE)
C               IF(CODE.EQ.0) THEN
C                 ERROR='Could not close command receipt message'
C                 GOTO 9999
C               ENDIF
                CALL PARSE(0,ISEG,CSEG,STRING,END,ERROR,*170)
                GOTO 220
C               Handle error condition
 170            CALL STRING_TRIM(ERROR,IBEG,IEND)
                WRITE(OP_STRING,'(A)') ' '//ERROR(IBEG:IEND)//
     '            '>MAIN_LOOP'
                CALL WRITES(IOER,OP_STRING,ERROR,*9998)
                ERROR(1:)=' '

C ***     Check for changes to condition handler
 220            IF(CHANGE_HANDLER) THEN
C GMH 12/2/97 Shouldn't this be reset?
                  CHANGE_HANDLER=.FALSE.
                  IF(FATAL_HANDLER) THEN
                    CALL LIB$ESTABLISH(NEW_HANDLER) !estab new handler
                  ELSE
                    CALL LIB$REVERT() ! Reverts to old condition handler
                  ENDIF
                ENDIF
C               Send command completion
C               CALL WH_INPUT_F_OPEN_MESSAGE(CMGUI_COMMAND_O,
C    '            0,5,CODE)
C               IF(CODE.EQ.0) THEN
C                 ERROR='Could not open command completion message'
C                 GOTO 9999
C               ENDIF
C               CALL WH_INPUT_F_CLOSE_MESSAGE(CMGUI_COMMAND_O,
C    '            CODE)
C               IF(CODE.EQ.0) THEN
C                 ERROR='Could not close command completion message'
C                 GOTO 9999
C               ENDIF
C               Tell CMGUI to pop down the prompt window if necessary
                CALL WH_INPUT_F_OPEN_MESSAGE(CMGUI_PROMPT_I,
     '            0,3,CODE)
                IF(CODE.EQ.0) THEN
                  ERROR='Could not open prompt message'
                  GOTO 9999
                ENDIF
                CALL WH_INPUT_F_CLOSE_MESSAGE(CMGUI_PROMPT_I,
     '            CODE)
                IF(CODE.EQ.0) THEN
                  ERROR='Could not close prompt message'
                  GOTO 9999
                ENDIF
C               Check for quitting
                IF (END) THEN !time to quit
                  CONTINUE=.FALSE.
                ENDIF
              ELSE
                WRITE(OP_STRING,'(''Command is longer than'',I5)') MXCH
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                CALL WH_OUTPUT_F_GET_REMAINDER(CMGUI_COMMAND_O,
     '            CODE)
                IF(CODE.EQ.0) THEN
                  ERROR='Could not use up command'
                  GOTO 9999
                ENDIF
              ENDIF
            ELSE
              ERROR='Invalid secondary type'
              GOTO 9999
            ENDIF !primary
            CALL WH_OUTPUT_F_CLOSE_MESSAGE(CMGUI_COMMAND_O,
     '        CODE)
            IF(CODE.EQ.0) THEN
              ERROR='Could not close command message'
              GOTO 9999
            ENDIF
          ELSE !no command
C           CMGUI_MESSAGE_TIME is when we entered the idle loop
            IF(CMGUI_MESSAGE_TIME.EQ.0) THEN
              CALL GET_SECONDS(CMGUI_MESSAGE_TIME)
              IDLE_SENT=.FALSE.
            ELSE
C GMH 29/1/97 Ensure that we dont have cmiss running with no cmgui
              CALL GET_SECONDS(CURRENT_TIME)
              IF(CURRENT_TIME-CMGUI_MESSAGE_TIME.GT.
     '          CMGUI_IDLE_TIME) THEN
                IF(IDLE_SENT) THEN
                  IF(CURRENT_TIME-CMGUI_MESSAGE_TIME
     '              .GT.2*CMGUI_IDLE_TIME) THEN !time to quit
                    CONTINUE=.FALSE.
                  ENDIF
                ELSE
                  CALL WH_INPUT_F_OPEN_MESSAGE(CMGUI_COMMAND_I,
     '              0,1,CODE)
                  IF(CODE.EQ.0) THEN
                    ERROR='Could not open idle message'
                    GOTO 9999
                  ENDIF
                  CALL WH_INPUT_F_CLOSE_MESSAGE(CMGUI_COMMAND_I,
     '              CODE)
                  IF(CODE.EQ.0) THEN
                    ERROR='Could not close idle message'
                    GOTO 9999
                  ENDIF
                  IDLE_SENT=.TRUE.
                ENDIF
              ENDIF !past idle
C             Try to send any information we may have in idle time
              CALL WH_INPUT_F_UPDATE(CMGUI_DATA_I,CODE)
              IF(CODE.EQ.0) THEN
                ERROR='Could not update data input'
                GOTO 9999
              ENDIF
              CALL WH_OUTPUT_F_UPDATE(CMGUI_DATA_O,CODE)
              IF(CODE.EQ.0) THEN
                ERROR='Could not update data output'
                GOTO 9999
              ENDIF
              SLEEP_DELAY=1 !seconds
              CALL SLEEPER(SLEEP_DELAY)
            ENDIF
          ENDIF

        ENDDO !main loop
C       Notify the front end by sending an empty message
C       Secondary ID=1
        CALL WH_INPUT_F_OPEN_MESSAGE(CMGUI_COMMAND_I,
     '    0,3,CODE)
        IF(CODE.EQ.0) THEN
          ERROR='Could not open quit message'
          GOTO 9999
        ENDIF
        CALL WH_INPUT_F_CLOSE_MESSAGE(CMGUI_COMMAND_I,
     '    CODE)
        IF(CODE.EQ.0) THEN
          ERROR='Could not close quit message'
          GOTO 9999
        ENDIF

        ELSE !cmgui_link
        IF(FIRST.EQ.0) THEN
          DO NOSG=1,MXSG
            ISEG(NOSG)=0
            CSEG(NOSG)=' '
          ENDDO
          COD(2)='PI'
          WRITE(COD(3),'(E11.5)') PI
          CALL ASSIGN(STRING,1,COD,COQUD,ERROR,*9999)

C***      Query mode
          ISTATUS = LIB$GETJPI(JPI$_MODE,,,,STATSTR,LEN)
          IF(.NOT.ISTATUS) CALL LIB$SIGNAL(%VAL(ISTATUS))

C***      Check whether batch mode
          IF(STATSTR(1:1).EQ.'B') THEN !Cannot use SMG read routines
            SMG_READ=.FALSE.
          ELSE !Set up control-c trap
            CALL SET_CTRLC_AST
          ENDIF

C CPB 26/10/95 Condition handler now changed in syntax
C ***     Establish condition handler

          FATAL_HANDLER=.TRUE.
          CHANGE_HANDLER=.FALSE.

CC CPB 6/10/93 Don't store old handler as OpenVMS cannot handle
CC the typing of the procedure. Use LIB$REVERT instead.
          OLD_HANDLER=LIB$ESTABLISH(NEW_HANDLER) !address of deflt handler

          FIRST=1

          USEPARAMFILE=PARAMETERS.EQ.1
          IF(USEPARAMFILE) THEN
            CALL CSTRINGLEN(CSTRLEN,PARAMETERSFILENAME)
            CALL C2FSTRING(PARAMETERSFILENAME,CSTRLEN,PARAMETERSFILE)
            PARAMETERSFILE(CSTRLEN+1:)=' '
            CALL STRING_TRIM(PARAMETERSFILE,IBEG,IEND)
            STRG='fem define parameters;r;'//PARAMETERSFILE(IBEG:IEND)
C           temporarily set FIRST params to .FALSE. so that arrays won't
C           be allocated in FEM until after params read in.
            FIRST_SYNTAX=.FALSE.
            CALL PARSE(0,ISEG,CSEG,STRG,END,ERROR,*9999)
            FIRST_SYNTAX=.TRUE.
          ENDIF

          CALL CSTRINGLEN(CSTRLEN,COMFNAME)
          CALL C2FSTRING(COMFNAME,CSTRLEN,COMFILENAME)
          COMFILENAME(CSTRLEN+1:)=' '
          CALL STRING_TRIM(COMFILENAME,IBEG,IEND)
          IF(EXAMPLE.EQ.1) THEN

C CPB 2/10/98 Use set directory command to set example path
            CALL CSTRINGLEN(CSTRLEN,EXNUM)
            CALL C2FSTRING(EXNUM,CSTRLEN,EXAMPLENUM)
            CALL STRING_TRIM(EXAMPLENUM,IBEG1,IEND1)
            CALL SETEXAMPLEDIR(EXAMPLENUM,ERROR,*9999)

            STRG='read command;'//COMFILENAME(IBEG:IEND)//';example'
          ELSE
            STRG='read command;'//COMFILENAME(IBEG:IEND)
          ENDIF
          CALL PARSE(0,ISEG,CSEG,STRG,END,ERROR,*9999)

          IF(END) THEN !quit issued from command file
            CONTINUE=.FALSE.
          ELSE
            CONTINUE=.TRUE.
          ENDIF

        ENDIF

        DO WHILE(CONTINUE) !is main program loop

C GMH 15/11/95 making call conditional
          IF(USE_GRAPHICS.EQ.1) THEN
            CALL GXWAIT(0.0,ERR) !Update graphics
          ENDIF

          IF(MACRO_COMMAND_EXECUTE) THEN !parse command macros
            MACRO_COMMAND_EXECUTE=.FALSE.
            DO nomacro=1,NT_MACRO_names(MACRO_command_ID)
              STRG=MACRO_COMMAND_buffer(nomacro,MACRO_command_ID)
              CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*150)
            ENDDO
          ELSE IF(MACRO_KEY_EXECUTE) THEN !parse key macros
            MACRO_KEY_EXECUTE=.FALSE.
            DO nomacro=1,NT_MACRO(macro_key)
              STRG=MACRO_KEY_buffer(nomacro,macro_key)
              CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*150)
            ENDDO
          ELSE IF(ITERATE_COMMAND_EXECUTE) THEN !parse command macros
            CALL STRING_TRIM(FILE00,IBEG,IEND)
            DO iter_counter=1,NT_ITERATE
C               Note that iter_counter is kept in iter00.cmn
              WRITE(CHAR4,'(I4)') iter_counter
              CALL STRING_TRIM(CHAR4,IBEG4,IEND4)
              FILE00=FILE00(IBEG:IEND)//'_'//CHAR4(IBEG4:IEND4)
              DO nomacro=1,NT_MACRO_names(MACRO_command_ID)
                STRG=MACRO_COMMAND_buffer(nomacro,MACRO_command_ID)
                CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*150)
              ENDDO !nomacro
            ENDDO !iter_counter
            ITERATE_COMMAND_EXECUTE=.FALSE.
          ELSE IF(DO_EXAMPLE) THEN !read selected example file
            DO_EXAMPLE=.FALSE.
            CALL STRING_TRIM(EXAMPLE_NAME,IBEG,IEND)
            STRING=' > read '//EXAMPLE_NAME(IBEG:IEND)//';com;example'
            WRITE(OP_STRING,'(A)') STRING(1:30)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE IF(SELECT_EXAMPLE) THEN !open selected example file
            SELECT_EXAMPLE=.FALSE.
            CALL STRING_TRIM(EXAMPLE_NAME,IBEG,IEND)
            STRING=' > open '//EXAMPLE_NAME(IBEG:IEND)//';com;example'
            WRITE(OP_STRING,'(A)') STRING(1:30)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE !set prompt and get string from command line
            CALL SETSTR(STRING,ERROR,*150)
            IF(SMG_READ) THEN
              CALL GETSTR1(STRING,PR(:LNPR),LNPR,END,ERROR,*150)
            ELSE
              CALL GETSTR(STRING(LNPR-2:),PR(:LNPR),END,ERROR,*150)
            ENDIF
          ENDIF

          IF(END) THEN !CTRL-Z has been used to quit
            CALL QUIT(END,ERROR,*150)
            CONTINUE=.FALSE.
            GOTO 200
          ENDIF

C***      Parse string
          IF(SMG_READ)THEN
            STRING=STRING(3:)
            CALL PARSE(0,ISEG,CSEG,STRING,END,ERROR,*150)
          ELSE
            CALL PARSE(0,ISEG,CSEG,STRING,END,ERROR,*150)
          ENDIF
          STRING=' ' !Test 4-3-1991

C ***     Check for changes to condition handler
          IF(CHANGE_HANDLER) THEN
            IF(FATAL_HANDLER) THEN
              CALL LIB$ESTABLISH(NEW_HANDLER) !establishs new condition handler
            ELSE IF(.NOT.FATAL_HANDLER) THEN
C
C CPB 13/7/93 - OpenVMS has changed the type of the argument to
C LIB$ESTABLISH. As a result changing back to the old condition handler
C can not be handled in this way (OLD_HANDLER must be a procedure
C value). Removing this call to LIB$ESTABLISH.
C
C              CALL LIB$ESTABLISH(OLD_HANDLER) !establishs VAX condition handler
C
C CPB 6/10/93 - Changing the way the old handler is restored. Use a LIB$REVERT
C rather than storing the address of the old handler. This should now be
C compatible with OpenVMS
C
C              WRITE(OP_STRING,'(''>>Can not reset condition '',
C     '          ''handler in OpenVMS'')')
C              CALL WRITES(IOER,OP_STRING,ERROR,*9999)
C
              CALL LIB$REVERT() ! Reverts to old condition handler
            ENDIF
          ENDIF

          IF(END) THEN !The QUIT command has been used to quit
            CONTINUE=.FALSE.
          ENDIF
          GOTO 200

C***        Handle error condition
 150      CALL STRING_TRIM(ERROR,IBEG,IEND)
          WRITE(OP_STRING,'(A)') ' '//ERROR(IBEG:IEND)//'>DIALOG'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          IF(CTRLC) THEN
            CTRLC=.FALSE.
C CPB 23/10/93 This is done when CTRLC_AST is called so no need to
C do it here
C            !if mode is not batch
C            IF(STATSTR(1:1).ne.'B')THEN
C              ISTATUS=SYS$QIOW(,%VAL(INPUT_CHAN),%VAL(ICODE),,,,
C     '          CTRLC_AST,,,,,)
C              IF(.NOT.ISTATUS) CALL LIB$SIGNAL(%VAL(ISTATUS))
C            ENDIF
          ENDIF
          ERROR(1:)=' '

 200      CONTINUE
        ENDDO !end of main program loop
        ENDIF !cmgui_link
      ELSE !USE_SOCKET = .TRUE.

        CONNID1=0
        CONNID2=1

        IF (FSKLISTEN(CONNID1,PORT1) .EQ. -1) GOTO 9998
        IF (FSKLISTEN(CONNID2,PORT2) .EQ. -1) GOTO 9998

C ***   Establish condition handler
        FATAL_HANDLER=.TRUE.
        CHANGE_HANDLER=.FALSE.
        CALL LIB$ESTABLISH(NEW_HANDLER) !estab new condit handler

        FIRST=1

        USEPARAMFILE=PARAMETERS.EQ.1
        IF(USEPARAMFILE) THEN
          CALL CSTRINGLEN(CSTRLEN,PARAMETERSFILENAME)
          CALL C2FSTRING(PARAMETERSFILENAME,CSTRLEN,PARAMETERSFILE)
          PARAMETERSFILE(CSTRLEN+1:)=' '
          CALL STRING_TRIM(PARAMETERSFILE,IBEG,IEND)
          STRG='fem define parameters;r;'//PARAMETERSFILE(IBEG:IEND)
C         temporarily set FIRST params to .FALSE. so that arrays won't
C         be allocated in FEM until after params read in.
          FIRST_SYNTAX=.FALSE.
          CALL PARSE(0,ISEG,CSEG,STRG,END,ERROR,*9999)
          FIRST_SYNTAX=.TRUE.
        ENDIF

        CONTINUE=.TRUE.
        DO WHILE (CONTINUE)
C CPB 14/3/94 Adding timeout for sockets to enable graphics updates
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
          IF(STRING(IBEG:IBEG+3).ne.'QUIT') THEN
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
            CALL QUIT(END,ERROR,*9999)
            CONTINUE=.FALSE.
          ENDIF

C ***     Check for changes to condition handler
          IF(CHANGE_HANDLER) THEN
            IF(FATAL_HANDLER) THEN
              CALL LIB$ESTABLISH(NEW_HANDLER) !estab new condit handler
            ELSE
              CALL LIB$REVERT() ! Reverts to old condition handler
            ENDIF
          ENDIF

        ENDDO

      ENDIF

C ??? MPN 27-oct-94 There must be a better place for this?
C news MPN 13-Sep-94: Parallel element stiffness matrix computations
      IF(KTYP1A.EQ.2) THEN !parallel element stiffness matrix calcs
        DO nhost=1,NUMHOSTS_USED
          IF(SOCKET_OPEN(nhost)) THEN
C           Signal to slave processes to stop
            IF(FSKWRITE(QUIT_PROCESS,SK_LONG_INT,1,ICONNID(nhost))
     '        .EQ.-1) GOTO 9999
            IRET=FSKCLOSE(ICONNID(nhost)) !close socket
          ENDIF
        ENDDO
      ENDIF
C newe

      RETURN
 9999 WRITE(*,'('' >>Fatal Error in Main Loop'')')
      WRITE(*,'('' >>Error : '',A)') ERROR
      GOTO 9997
 9998 WRITE(*,'('' >>Fatal Socket Error in Main Loop'')')
 9997 IF(USE_SOCKET) THEN
        IRET=FSKCLOSE(CONNID1)
        IRET=FSKCLOSE(CONNID2)
      ENDIF
      RETURN
      END


      SUBROUTINE OPENF(IUNIT,DEVICE,FILE,STATUS,ACCESS,FORM,
     '  IRECL,ERROR,*)

C#### Subroutine: OPENF
C###  Description:
C###    OPENF opens a file.

C**** IUNIT is the name by which the file is referred to in the
C**** program.  Valid unit identifiers are non negative integers.
C**** DEVICE is the device on which the file resides.
C**** Valid devices are: 'DISK' 'TERM'
C**** FILE is the filename [filetype [filemode]]
C**** STATUS is 'NEW' or 'OLD' or 'SCRATCH'
C**** ACCESS specifies whether file is 'DIRECT' or 'SEQUEN' or 'APPEND'
C**** Direct access files have a length of 2000 records.
C**** FORM specifies whether the file is 'FORMATTED' or 'UNFORMATTED'
C**** IRECL is the logical record length
C**** ERROR gives diagnostics in the event of failure.
C**** If an error is detected control is returned to the statement
C**** number of the star.
C**** DEVICE, FILE, ACCESS and ERROR are all character strings.
C**** NOTE: The keyword /RECORDTYPE in OPEN is nonstandard!

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cmis00.cmn'
      INCLUDE 'cmiss$reference:disp00.cmn'
      INCLUDE 'cmiss$reference:docd00.cmn'
      INCLUDE 'cmiss$reference:gino00.cmn'
!     Parameter List
      INTEGER IRECL,IUNIT
      CHARACTER ACCESS*(*),DEVICE*(*),ERROR*(*),FILE*(*),
     '  STATUS*(*),FORM*(*)
!     Local Variables
      INTEGER CLOCAT,i,i_tot,
     '  IBEG,IBEG1,IBEG2,IBEG3,
     '  IEND,IEND1,IEND2,IEND3,
     '  IOSTAT,IREC,I_end_of_path
      CHARACTER DUMMY*256,FILENAME*200,SFILE*100,RECL*5,
     '  RECORDTYPE*10,SUB_DIR*100,UNIT*4
      LOGICAL LAST

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
      WRITE(UNIT,'(I4)') IUNIT
      WRITE(RECL,'(I5)') IRECL
      Modify_File=.FALSE.

      IF(DOCDIR) THEN !read from documentation directory

C!!! cpb 2/10/98 Set example path with a set directory command
C        IF(Open_COM_file) THEN !open com file from backend CMISS
C          Backend_COM_file=.TRUE.
C
C          IBEG1=1
C          IEND1=CLOCAT(']',EXAMPLE_PATH)
C          IF(DOP) THEN
C            WRITE(OP_STRING,'('' example_path='',A)')
C     '        EXAMPLE_PATH(IBEG1:IEND1)
C            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C          ENDIF
C
CC Find the example subdirectory from the com file name
C          IBEG2=8                   !e.g. example_31   ibeg2=8
C          IEND2=CLOCAT('.',FILE)-1                    !iend2=10
Cc         write(*,'('' ibeg2='',I4,''iend2='',I4)')   ibeg2,iend2
C          SUB_DIR=FILE(IBEG2+1:IBEG2+1)               !sub_dir='3'
C          i_tot=IEND2-IBEG2                           !i_tot=2
C          DO i=2,i_tot                                !do i=2,2
C            CALL STRING_TRIM(SUB_DIR,IBEG3,IEND3)
Cc           write(*,'('' ibeg3='',I4,''iend3='',I4)') ibeg3,iend3
C            SUB_DIR=SUB_DIR(1:IEND3)//'.'//FILE(IBEG2+1:IBEG2+i)
C            if(dop) write(*,'('' sub_dir='',A)') sub_dir(1:60)
C          ENDDO                                       !sub_dir='3.31'
C
CC Append the example subdirectory to the example_path
C          if(dop) write(*,'('' ibeg1='',I4,''iend1='',I4)') ibeg1,iend1
C          CALL STRING_TRIM(SUB_DIR,IBEG3,IEND3)
C          if(dop) write(*,'('' ibeg3='',I4,''iend3='',I4)') ibeg3,iend3
C          EXAMPLES_DIR=EXAMPLE_PATH(IBEG1:IEND1-1)
C     '      //'.'//SUB_DIR(IBEG3:IEND3)//']'
C          CALL STRING_TRIM(EXAMPLES_DIR,IBEG1,IEND1)
C          if(dop) then
C            write(op_string,'('' examples_dir='',a)')
C     '        examples_dir(ibeg1:iend1)
C            call writes(iodi,op_string,error,*9999)
C          endif
C
CC Retain the example directory string for later use
C          Example_directory=EXAMPLES_DIR(IBEG1:IEND1)
C          if(dop) then
C            write(op_string,'('' example_directory='',a)')
C     '        example_directory(ibeg1:iend1)
C            call writes(iodi,op_string,error,*9999)
C          endif
C
C        ELSE !define commands from an already opened COM file
C          IF(Backend_COM_file) THEN
C            EXAMPLES_DIR=Example_directory !defined when OPEN file
C          ELSE !an examples subdirectory has been set
C            IBEG1=1
C            IEND1=CLOCAT(']',EXAMPLE_PATH)
CC Append the example subdirectory to the example path
C            CALL STRING_TRIM(Example_subdir,IBEG2,IEND2)
C            IF(Example_subdir(IEND2:IEND2).EQ.'/') IEND2=IEND2-1
CC         Replace Unix '/' with VMS '.'
C            CALL INSERT(Example_subdir,'/','.')
C            EXAMPLES_DIR=EXAMPLE_PATH(IBEG1:IEND1-1)
C     '        //Example_subdir(IBEG2:IEND2)//']'
C          ENDIF
C
C        ENDIF !Open_COM_file

C Append filename to examples pathname

        CALL STRING_TRIM(EXAMPLES_DIR,IBEG1,IEND1)
        I_end_of_path=IEND1-IBEG1+1
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
C     Locate position in FILE of end of directory path string
        I=0
        I_end_of_path=0
        LAST=.FALSE.
        DO WHILE (.NOT.LAST)
          I=CLOCAT('/',FILE(I_end_of_path+1:))
          if(dop) write(*,'('' i='',I2)') i
          IF(I.EQ.0) THEN
            LAST=.TRUE.
          ELSE
            I_end_of_path=I_end_of_path+I
            if(dop) write(*,'('' I_end_of_path='',I2)') I_end_of_path
          ENDIF
        ENDDO !while .not.last

C     Replace Unix pathname with VMS pathname
        IF(I_end_of_path.GT.0) THEN !filename contains path info
          IF(FILE(1:1).EQ.'/') THEN !absolute filepath
C cpb 18/4/95 Changing this so that if an absolute path is specified
C the first 'directory' is a VMS diskname
C            FILENAME(1:1)='['
C            FILENAME(2:I_end_of_path-1)=FILE(2:I_end_of_path-1)
C            FILENAME(I_end_of_path:I_end_of_path)=']'
C            FILENAME(I_end_of_path+1:)=FILE(I_end_of_path+1:)
            I=CLOCAT('/',FILE(2:))
            IF(I.EQ.0) THEN
              ERROR='>>Invalid filename'
              GOTO 9999
            ENDIF
            FILENAME=FILE(2:I)
            FILENAME(I:I+1)=':['
            FILE(I_end_of_path:I_end_of_path)=']'
            FILENAME(I+2:)=FILE(I+2:)
          ELSE                      !relative filepath
            FILENAME(1:2)='[.'
            FILENAME(3:I_end_of_path+1)=FILE(1:I_end_of_path-1)
            FILENAME(I_end_of_path+2:I_end_of_path+2)=']'
            FILENAME(I_end_of_path+3:)=FILE(I_end_of_path+1:)
          ENDIF
          CALL INSERT(FILENAME,'/','.')
        ELSE
          FILENAME(1:)=FILE(1:)
        ENDIF
        if(dop) then
          write(op_string(1),'(''   filename = '', a)') filename
          call writes(iodi,op_string,error,*9999)
        endif

      ENDIF !docdir

      IF(DEVICE.EQ.'TERM') THEN
        OPEN(UNIT=IUNIT,FILE=FILENAME,STATUS=STATUS,RECL=IRECL,
     '    IOSTAT=IOSTAT)

      ELSE IF(DEVICE.EQ.'DISK') THEN
        IF(ACCESS.EQ.'SEQUEN') THEN
          OPEN(UNIT=IUNIT,FILE=FILENAME,STATUS=STATUS,
     '      ACCESS='SEQUENTIAL',FORM=FORM,IOSTAT=IOSTAT,
     '      CARRIAGECONTROL='NONE',RECORDTYPE='STREAM_LF')

        ELSE IF(ACCESS.EQ.'APPEND') THEN
          OPEN(UNIT=IUNIT,FILE=FILENAME,STATUS=STATUS,
     '      ACCESS='APPEND',FORM=FORM,IOSTAT=IOSTAT,
     '      CARRIAGECONTROL='NONE',RECORDTYPE='STREAM_LF')

        ELSE IF(ACCESS.EQ.'DIRECT') THEN
          IF(STATUS.EQ.'OLD') THEN
C           open sequential file in specified directory
            OPEN(UNIT=99,FILE=FILENAME,STATUS=STATUS,
     '        ACCESS='SEQUENTIAL',FORM=FORM,IOSTAT=IOSTAT,
     '        CARRIAGECONTROL='NONE',RECORDTYPE='STREAM_LF',
     '        READONLY)

            IF(IOSTAT.EQ.0) THEN
C             open direct access scratch file in local directory
              SFILE='OLD_'//FILENAME(I_end_of_path+1:)
              IF(DOP) WRITE(*,'('' Scratch file='',A)') SFILE
              OPEN(UNIT=IUNIT,FILE=SFILE,STATUS='SCRATCH',
     '          ACCESS='DIRECT',FORM=FORM,IOSTAT=IOSTAT,RECL=IRECL)
              IREC=0
 10           READ(99,'(A)',IOSTAT=IOSTAT,END=20) DUMMY
              IF(IOSTAT.EQ.0) THEN
C cpb 20/4/94 This inquire doesn't seem to pick up the next record
C properly
C          INQUIRE(UNIT=IUNIT,NEXTREC=IREC)
                IREC=IREC+1
                WRITE(UNIT=IUNIT,FMT='(A)',REC=IREC,IOSTAT=IOSTAT)
     '            DUMMY(1:IRECL)
                IF(DOP) WRITE(*,'('' DUMMY='',A)') DUMMY(1:60)
                GOTO 10
              ELSE !file error
                GOTO 20
              ENDIF
 20           CLOSE(UNIT=99)
              IF(IOSTAT.EQ.-1) THEN !EOF
                IOSTAT=0
              ENDIF

            ELSE IF(IOSTAT.EQ.44) THEN !Tell user to convert the file.
              WRITE(OP_STRING,'('' >>Convert DIRECT access file to '','
     '          //'''a SEQUENTIAL access STREAM_LF file using the '
     '          //'TOSLF command'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF

          ELSE IF(STATUS.EQ.'NEW') THEN
            SFILE='NEW_'//FILENAME
            OPEN(UNIT=IUNIT,FILE=SFILE,STATUS='SCRATCH',ACCESS='DIRECT',
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

      IF(iostat.ne.0) THEN
        WRITE(ERROR,'(I3)') IOSTAT
        ERROR=' Iostat='//ERROR(1:3)//' error in OPENF(UNIT='
     '    //UNIT//')'
        IF(IOSTAT.EQ.29) THEN
          CALL STRING_TRIM(ERROR,IBEG,IEND)
          ERROR=' File not found: '//ERROR(IBEG:IEND)
        ELSE IF(IOSTAT.EQ.25) THEN
          CALL STRING_TRIM(ERROR,IBEG,IEND)
          ERROR=' Record number outside range: '//ERROR(IBEG:IEND)
        ELSE IF(IOSTAT.EQ.30) THEN
          CALL STRING_TRIM(ERROR,IBEG,IEND)
          ERROR=' Open failure: '//ERROR(IBEG:IEND)
        ELSE IF(IOSTAT.EQ.32) THEN
          CALL STRING_TRIM(ERROR,IBEG,IEND)
          ERROR=' Invalid unit number: '//ERROR(IBEG:IEND)
        ELSE IF(IOSTAT.EQ.34) THEN
          CALL STRING_TRIM(ERROR,IBEG,IEND)
          ERROR=' Unit already open: '//ERROR(IBEG:IEND)
        ELSE IF(IOSTAT.EQ.37) THEN
          CALL STRING_TRIM(ERROR,IBEG,IEND)
          ERROR=' Inconsistent record length: '//ERROR(IBEG:IEND)
        ELSE IF(IOSTAT.EQ.43) THEN
          CALL STRING_TRIM(ERROR,IBEG,IEND)
          ERROR=' Filename invalid: '//ERROR(IBEG:IEND)
        ELSE IF(IOSTAT.EQ.44) THEN
          CALL STRING_TRIM(ERROR,IBEG,IEND)
          ERROR=' Inconsistent record type: '//ERROR(IBEG:IEND)
        ELSE IF(IOSTAT.EQ.51) THEN
          CALL STRING_TRIM(ERROR,IBEG,IEND)
          ERROR=' Inconsistent file organisation: '
     '      //ERROR(IBEG:IEND)
        ENDIF
        WRITE(OP_STRING,'('' Filename: '',A)') FILENAME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        GOTO 9999
      ELSE
        ERROR='  '
      ENDIF

      CALL EXITS('OPENF')
      RETURN
 9999 CALL ERRORS('OPENF',ERROR)
      CALL EXITS('OPENF')
      RETURN 1
      END


      SUBROUTINE RESET_HANDLER(ERROR,*)

C#### Subroutine: RESET_HANDLER
C###  Description:
C###    Resets the error handler to the old (default) error handling
C###    routines.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:dial00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables

      CALL ENTERS('RESET_HANDLER',*9999)

      CHANGE_HANDLER=.TRUE.
      FATAL_HANDLER=.FALSE.

      CALL EXITS('RESET_HANDLER')
      RETURN
 9999 CALL ERRORS('RESET_HANDLER',ERROR)
      CALL EXITS('RESET_HANDLER')
      RETURN 1
      END


      SUBROUTINE SET_CTRLC_AST()

C#### Subroutine: SET_CTRLC_AST
C###  Description:
C###    SET_CTRLC_AST establishes the asynchrouns CTRL/C handler so
C###    that CTRLC_AST is called when CTRL/C is pressed.

      IMPLICIT NONE
      INCLUDE '($SYSSRVNAM)'
      INCLUDE '($IODEF)'
      INCLUDE 'cmiss$reference:ctrl00.cmn'
!     Parameter List
!     Local Variables
      INTEGER ISTATUS,IOSB(2)

      EXTERNAL CTRLC_AST

C Assign channel if not already assigned
      IF(INPUT_CHAN.EQ.0) THEN
        ISTATUS=SYS$ASSIGN('SYS$INPUT',INPUT_CHAN,,)
        IF(.NOT.ISTATUS) CALL LIB$SIGNAL(%VAL(ISTATUS))
      ENDIF
C Enable AST so that CTRLC_AST is called when CTRL/C is pressed
      ICODE=IO$_SETMODE.OR.IO$M_CTRLCAST
      ISTATUS=SYS$QIOW(,%VAL(INPUT_CHAN),%VAL(ICODE),IOSB,,,
     '  CTRLC_AST,,,,,)
      IF(.NOT.ISTATUS) CALL LIB$SIGNAL(%VAL(ISTATUS))

      RETURN
      END


      SUBROUTINE SET_HANDLER(ERROR,*)

C#### Subroutine: SET_HANDLER
C###  Description:
C###    Sets the error handling routines from the old (default) routines

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:dial00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables

      CALL ENTERS('SET_HANDLER',*9999)

      CHANGE_HANDLER=.TRUE.
      FATAL_HANDLER=.TRUE.

      CALL EXITS('SET_HANDLER')
      RETURN
 9999 CALL ERRORS('SET_HANDLER',ERROR)
      CALL EXITS('SET_HANDLER')
      RETURN 1
      END


      SUBROUTINE SETUPCMISS(ERROR,*)

C#### Subroutine: SETUPCMISS
C###  Description:
C###    Sets up up the initial platform dependent variables for CMISS.

      IMPLICIT NONE
      INCLUDE '($IODEF)'
      INCLUDE '($JPIDEF)'
      INCLUDE '($SMGDEF)'
      INCLUDE '($TRMDEF)'
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
      INCLUDE 'cmiss$reference:ctrl00.cmn'
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
      INCLUDE 'cmiss$reference:iter00.cmn'
      INCLUDE 'cmiss$reference:ktyp00.cmn'
      INCLUDE 'cmiss$reference:mach00.cmn'
      INCLUDE 'cmiss$reference:mach00.inc'
      INCLUDE 'cmiss$reference:map000.cmn'
      INCLUDE 'cmiss$reference:ntsg00.cmn'
      INCLUDE 'cmiss$reference:tol00.cmn'
      INCLUDE 'cmiss$reference:trac00.cmn'
!     Parameter List
      CHARACTER ERROR*(MXCH)
!     Local Variables
      INTEGER MXCO,MXCOQU,MXSG
      PARAMETER (MXCO=25,MXCOQU=25,MXSG=10000)
      INTEGER i,ISEG(MXSG),macro_command,no_col,
     '  nosg,nsub,NT_VIEW
      REAL*8 DLAMCH
      LOGICAL FIRST_ZOOM,
     '  REFINE_ACTIVE,TRANSFORM_ACTIVE,UPDATE_ZOOM
      CHARACTER CSEG(MXSG)*60

      EXTERNAL CTRLC_AST,NEW_HANDLER

      CALL ENTERS('SETUPCMISS',*9999)

      DO i=1,100
        OP_STRING(i)(1:1)=CHAR(0)
      ENDDO
      SMG_READ=.TRUE.
      INPUT_CHAN=0
c cpb 22/10/95 Set machine constants (don't change)
      CHARSIZE=1
      INTSIZE=4
      SINTSIZE=2
      LINTSIZE=8
      SPSIZE=4
      DPSIZE=8
      LOGSIZE=4
      SPCSIZE=8
      DPCSIZE=16
      ENDIANTYPE=CHAR(MACH_LITTLEENDIAN) !Little endian native
      SPFORMTYPE=CHAR(MACH_SPIEEE) !IEEE single precision format
      DPFORMTYPE=CHAR(MACH_DPIEEE) !IEEE double precision format
      CHARFORMTYPE=CHAR(MACH_CHARASCII) !ASCII character format
      INTFORMTYPE=CHAR(MACH_TWOSCOMP) !Two's complement format
      MACHTYPE=CHAR(MACH_DECALPHA) !Alpha
      OSTYPE=CHAR(MACH_VMS) !VMS
C      NUM_LIBRARY=0
C      GKS=.FALSE.
C      MAPOPN=.FALSE.
C      ZERO_TOL=1.0d-10
C      CONVERG_TOL=DLAMCH('EPS')
C      LOOSE_TOL=DSQRT(CONVERG_TOL)
c
C      TR01=.FALSE.
C      TR02=.FALSE.
C      TR03=.FALSE.
C      TR04=.FALSE.
C      FIRST_SYNTAX=.TRUE.
C      FEM_ARRAYS_MALLOCED=.FALSE.
C      FIRST_ZOOM =.TRUE.
C      UPDATE_ZOOM=.FALSE.
C      BUFFER_COUNT=0
C      IREC_COMFILE=0
C      CALL GET_REV_TIME(ERROR,*9999)
C      CALL GET_SYSTEM(ERROR,*9999)
C      CALL GET_DISPLAY(ERROR,*9999)
C      IOIP=1 !input  file
C      IOOP=2 !usual output
C      IODI=3 !diagnostics output
C      IOTR=4 !trace output
C      IOER=5 !error output
C      IOH1=6 !help 1 output
C      IOH2=7 !help 2 output
C      IOOUT=9 !output FILE
C      CALL OPENF(1,'TERM','SYS$INPUT' ,'NEW',' ',' ',132,ERROR,*9999)
C      CALL OPENF(2,'TERM','SYS$OUTPUT','NEW',' ',' ',255,ERROR,*9999)
C      CALL OPENF(3,'TERM','SYS$OUTPUT','NEW',' ',' ',255,ERROR,*9999)
C      CALL OPENF(4,'TERM','SYS$OUTPUT','NEW',' ',' ',255,ERROR,*9999)
C      CALL OPENF(5,'TERM','SYS$OUTPUT','NEW',' ',' ',255,ERROR,*9999)
C      CALL OPENF(6,'TERM','SYS$OUTPUT','NEW',' ',' ',255,ERROR,*9999)
C      CALL OPENF(7,'TERM','SYS$OUTPUT','NEW',' ',' ',255,ERROR,*9999)
C      CALL OPENF(8,'TERM','SYS$OUTPUT','NEW',' ',' ',255,ERROR,*9999)
C      IOFI=IOOP !default for listing to file commands
C      IFILE=10 !ifile is the file for input files i.e. ip* files
C      IOFILE1=11 !output file i.e. op* files
C      IOFILE2=12 ! iofile2 -> iofile6 are general purpose output files
C      IOFILE3=13
C      IOFILE4=14
C      IOFILE5=17
C      IOFILE6=18
CC CPB 30/11/92 the IO# units below should (?) be redundant.
C      IO1=1
C      IO2=2
C      IO3=IOIP
C      IO4=IOOP
C      IVDU=IOIP
C      IMAP=0
C      ECHO_OUTPUT=.FALSE.
C      DO num_stack=1,5                !<-|
C        DOP_STACK(num_stack)=.FALSE.  !  |
C      ENDDO                           !  |
C      num_stack=5 !5 levels below FEM !  |- used for diagnostics
C      NT_SUB=0                        !  |
C      DO nsub=1,MXSUB                 !  |
C        SUBNAM(nsub)=' '              !  |
C      ENDDO                           !  |
C      DIAGNO=.FALSE.                  !  |
C      ALLSUB=.FALSE.                  !  |
C      FROMSUB=.FALSE.                 !  |
C      DOP=.FALSE.                     !<-|
C      NTSG=0
C      NT_KEY=0
C      NT_MACRO(0)=0
C      NT_MACRO_names(0)=0
C      DO macro_command=1,9
C        NT_MACRO_names(macro_command)=0
C        MACRO_names(macro_command)=' '
C      ENDDO
C      MACRO_DEFINE=.FALSE.
C      MACRO_COMMAND_EXECUTE=.FALSE.
C      MACRO_KEY_EXECUTE=.FALSE.
C      REFINE_ACTIVE=.FALSE.
C      TRANSFORM_ACTIVE=.FALSE.
C      DO_EXAMPLE=.FALSE.
C      SELECT_EXAMPLE=.FALSE.
C      NT_VIEW=0
C      HEADING=' '
C      DO no_col=1,9
C        NAME_COL(no_col)='UNDEFINED'
C      ENDDO
C      DO i=1,99
C        IWKS(i)=0
C      ENDDO
Cc     IWKDEF(0)=1
Cc     IWKDEF(1)=1
C      DO nosg=1,MXSG
C        ISEG(nosg)=0
C        CSEG(nosg)=' '
C      ENDDO
C      XMIN=-1.0
C      XMAX= 1.0
C      YMIN=-1.0
C      YMAX= 1.0
C      ZMIN=-1.0
C      ZMAX= 1.0
C      DIAG=SQRT(12.0)
C      NJT=2
C      FILE00='file'
C      PATH00=' '
C      MENU=.FALSE.

      CALL EXITS('SETUPCMISS')
      RETURN
 9999 CALL ERRORS('SETUPCMISS',ERROR)
      CALL EXITS('SETUPCMISS')
      RETURN 1
      END


      SUBROUTINE SLEEPER(IDELAY)

C#### Subroutine: SLEEPER
C###  Description:
C###    SLEEPER pauses for an integer period of seconds.

      IMPLICIT NONE
!     Parameter List
      INTEGER IDELAY

      CALL SLEEP(%VAL(IDELAY))  !C run-time library routine

      RETURN
      END
