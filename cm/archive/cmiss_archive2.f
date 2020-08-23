C#### Module: CMISS_ARCHIVE2
C###  Description:
C###    Contains archived code from modules FE10 -> FE19

CFE10 Function TIMER        call RTL to return system times
C     Subroutine DIALOG
C     Subroutine FIND_FILE     find files in current directory
C###  Routine: GET_CMISS_EXAMPLES return CMISS examples directory
C     Subroutine OPENF
C     Subroutine POST_FILE     post file to printer
C     Subroutine PURGE_FILE    purge file versions
                
CFE11 Subroutine IPFIT      input opt.n params for geometry or field fit
C###  Routine: IPFOUR   input Fourier Transform analysis parameters
C     Subroutine IPSHEE     input sheet direction field
C     Subroutine IPGRID     input grid parameters    
C     Subroutine IPTIME     input time variables
              
CFE12 Subroutine OPCOLO     output workstation index colours
C     Subroutine OPFIBR     output fibre direction field
C     Subroutine OPMATR     output matrix
C     Subroutine OPST80     output strains at Gauss points.
C###  Routine: OPTEXT   output text
C     Subroutine OPTIME     output time variables
C     SUBROUTINE OPTRAN     outputs Phigs transformation
C###  Routine: OPVSAE   output VSaero parameters

CFE13 Subroutine IPMESH2    old version of mesh routine for lungs
C###  Routine: IPMESH9_DYNAM  input arbitrary 2D domain mesh
C###  Routine: IPMESH9  dynamic memory wrapper for IPMESH9_DYNAM

CFE14
C###  Routine: CELL_ARRAY       draws bit-image cell array
C###  Routine: CREATE_SEGMENT   create graphics seg without ISEG,CSEG
C###  Routine: ELLIPSE          draws fill-area ellipse
C     Subroutine PRINT_IMAGE_FILL_AREAS   not used?

CFE19 Subroutine ADAMS_MOULTON ode solver
C     Subroutine AM_DE         used by adams_moulton
C     Subroutine AM_STEP         "       "       "
C     Subroutine AM_INTERPOLATE  "       "       "
C     Subroutine ADAMS         older version of above ode solver
C     Subroutine LR            old version of LR2 equations
C     Subroutine LR_CURRENTS     "       "       "
C     Subroutine L_R_RATES       "       "       "
C     Subroutine GATES      computes time constants and activation vars
C     Subroutine NOBLE98       old version of Noble 98
C     Subroutine DEFINED_NOBLE98 "       "       "
C     Subroutine NOBLE98_RATES   "       "       "
C     Subroutine NOBLE98_CHANGE  "       "       "
C     Subroutine NOBLE98_CURRENTS"       "       "
C     Function   FN_CA_CELL    HMT Cai
C     Function   FN_TN_CELL    HMT troponin kinetics
C     Function   FN_TM_CELL    HMT tropomyosin kinetics
C     Function   FN_TO_CELL    HMT isometric tension
C     Function   FN_Q_CELL     HMT fading memory model

Module FE10
=========== 

      
      REAL*8 FUNCTION TIMER(IHANDLE,T_FLAG)

C#### Function: TIMER
C###  Type: REAL*8
C###  Description:
C###    TIMER initialises timer if STATUS=0, otherwise returns elapsed 
C###    time (if T_FLAG=T_ELAPSED) or CPU time (if T_FLAG=T_CPU) in 
C###    seconds since initialisation.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:time02.cmn'
!     Parameter List
      INTEGER IHANDLE,T_FLAG
!     Local Variables
      INTEGER IRET_CPU,ISTATUS,LIB$INIT_TIMER,LIB$STAT_TIMER
      INTEGER*8 IRET_ELAPSED

      IF(T_FLAG.EQ.T_INITIALISE) THEN
        ISTATUS=LIB$INIT_TIMER(IHANDLE)
        TIMER=0.0D0
      ELSE IF(T_FLAG.EQ.T_CPU) THEN
        ISTATUS=LIB$STAT_TIMER(T_CPU,IRET_CPU,IHANDLE)
        TIMER=DBLE(IRET_CPU)/100.0D0
      ELSE IF(T_FLAG.EQ.T_ELAPSED) THEN
        ISTATUS=LIB$STAT_TIMER(T_ELAPSED,IRET_ELAPSED,IHANDLE)
C       Elapsed timed is returned in 100 nanosecond units and is
C       negative since a 'delta' time is returned (check this!)
        TIMER=-DBLE(IRET_ELAPSED)*100.0d-9
      ENDIF

      RETURN
      END


      SUBROUTINE DIALOG

C#### Subroutine: DIALOG
C###  Description:
C**** The dialog interpreter
C**** Note: pi and PI are assigned as user-defined names
C**** Condition Handler References:
C****   ch2 Introduction to System Routines
C****   ch4 LIB$ Manual -An Overview of VAX Condition Handling 
C****   ch9 VAX Fortran User Manual
C****   ch8 VAX C Run time library reference manual
C****   section 8.4 VAX Debugger Manual
C****   ch9 Guide to Programming Resources ***THIS IS THE BEST***
C****   ch10 Introduction to VMS System Services
C**** ISTATUS is >0 on successful complet.n of RTL call (ISTATUS=.true.)

      IMPLICIT NONE
      INCLUDE '($IODEF)'
      INCLUDE '($JPIDEF)'
      INCLUDE '($SMGDEF)'
      INCLUDE '($TRMDEF)'
      INCLUDE 'cmiss$reference:fsklib.inc'
      INCLUDE 'cmiss$reference:mxch.inc'
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:call00.cmn'
      INCLUDE 'cmiss$reference:cbdi00.cmn'
      INCLUDE 'cmiss$reference:cbdi01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:cbfe01.cmn'
      INCLUDE 'cmiss$reference:cbpr00.cmn'
      INCLUDE 'cmiss$reference:cbsy01.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:cbzm00.cmn'
      INCLUDE 'cmiss$reference:cmis00.cmn'
      INCLUDE 'cmiss$reference:cspi00.cmn'
      INCLUDE 'cmiss$reference:ctrl00.cmn'
      INCLUDE 'cmiss$reference:diag00.cmn'
      INCLUDE 'cmiss$reference:dial00.cmn'
      INCLUDE 'cmiss$reference:docu00.cmn'
      INCLUDE 'cmiss$reference:draw00.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:gen000.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:gks000.cmn'
      INCLUDE 'cmiss$reference:gks001.cmn'
      INCLUDE 'cmiss$reference:gks002.cmn'
      INCLUDE 'cmiss$reference:gks003.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'cmiss$reference:grou00.cmn'
      INCLUDE 'cmiss$reference:gtstr00.cmn'
      INCLUDE 'cmiss$reference:head00.cmn'
      INCLUDE 'cmiss$reference:host00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:ktyp00.cmn'
      INCLUDE 'cmiss$reference:map000.cmn'
      INCLUDE 'cmiss$reference:ntsg00.cmn'
      INCLUDE 'cmiss$reference:obje00.cmn'
      INCLUDE 'cmiss$reference:oxs000.cmn'
      INCLUDE 'cmiss$reference:phig00.cmn'
      INCLUDE 'cmiss$reference:plin00.cmn'
      INCLUDE 'cmiss$reference:trac00.cmn'
      INCLUDE 'cmiss$reference:view00.cmn'
!     Local Variables
      INTEGER MXCO,MXCOQU,MXSG
      PARAMETER (MXCO=16,MXCOQU=16,MXSG=2000)
      INTEGER ERR,i,IBEG,ID,ID_DEVICE,ID_STATUS,ID_WS,IEND,IFROMC,
     '  INDEX_PLIN,INLIST,INPUT_CHOICE,INTSTR(MXCH),IRET,
     '  ISEG(MXSG),ISEGM,ISEGM_LIST(0:100),ISTATUS,IW,LEN,
     '  LIB$ESTABLISH,LIB$GETJPI,LIB$SPAWN,macro_command,
     '  MODE_PROJ,n,N1LIST,NCENTRE,NCHAR,NEW_HANDLER,nj,nhost,NOCH,NOCO,
     '  no_col,NO_GROUP_TRANSFORM,nogrpl,nomacro,NOPTS,
     '  nosg,nsub,NTCH,NTFILE,NT_VIEW,RETURNVAL
      REAL DX,DY,FANGLE(3),FPT(3),FSCALE(3),FSHFT(3),
     '  R4DATA(2),XNDC1,XNDC2,XNDC3,XNDC4,VALUE,VMATRIX(4,4),
     '  X0,XCENTRE,Y0,YCENTRE
      LOGICAL CHANGE,CONTINUE,DEFINE,END,FIRST_ZOOM,
     '  REFINE_ACTIVE,TRANSFORM_ACTIVE,UPDATE,UPDATE_ZOOM,
     '  DISPLAY_VALUATOR_81,DISPLAY_VALUATOR_82,DISPLAY_VALUATOR_83,
     '  DISPLAY_VALUATOR_84,DISPLAY_VALUATOR_85,DISPLAY_VALUATOR_86,
     '  DISPLAY_VALUATOR_87,DISPLAY_VALUATOR_88,DISPLAY_VALUATOR_89
      CHARACTER CFROMI*5,CFROMR*11,CHAR3*3,CHAR5*5,CHOOSE*40,
     '  CIW*1,CLASS*8,CO(MXCO)*30,COD(MXCO)*90,COQUD(MXCO,MXCOQU)*30,
     '  CSEG(MXSG)*60,ERROR*(MXCH),FILE_NAME*50,OPTION(30)*40,
     '  OPTION1(19)*11,OPTION2(12)*11,OPTION3(26)*11,OPTION4(11)*11,
     '  OPTION5(9)*20,OPTION10(10)*20,OPTION11(10)*20,
     '  PARAMETER_TYPE*20,Q*1,SDATA*10,STATSTR*80,STRG*(MXCH),
     '  STRING*(MXCH),TEXT_STRING*80,TRANSFORM_TYPE*7

      DATA NXCH/MXCH/
      DATA MXIN/2147483647/
      DATA RXRE/1.70141178D+38/
      DATA RNRE/2.93873588D-39/
      DATA E/2.71828182846D0/

      EXTERNAL CTRLC_AST,NEW_HANDLER
      DATA OPTION1/'Archive..',
     '             'Cancel..',
     '             'Change..',
     '             'Check..',
     '             'Define..',
     '             'Display..',
     '             'Draw..',
     '             'Hide/Show..',
     '             'Label..',
     '             'List..',
     '             'Pick..',
     '             'Print..',
     '             'Read..',
     '             'Recall..',
     '             'Refine..',
     '             'Transform..',
     '             'Write..',
     '             '...',
     '             'Exit'/
      DATA OPTION2/'Archive..',
     '             'Cancel..',
     '             'Change..',
     '             'Define..',
     '             'Draw..',
     '             'Hide/Show..',
     '             'Label..',
     '             'Pick..',
     '             'Print..',
     '             'Recall..',
     '             'Transform..',
     '             'Exit'/
      DATA OPTION3/'Archive..',
     '             'Back plane',
     '             'Cancel..',
     '             'Define..',
     '             'Draw..',
     '             'Front plane',
     '             'Hide/Show..',
     '             'Label..',
     '             'Pan',
     '             'Parallel',
     '             'Perspective',
     '             'Print..',
     '             'Proj ref pt',
     '             'Recall..',
     '             'Rescale',
     '             'Reset',
     '             'Rotate data',
     '             'Rotate view',
     '             'Save view',
     '             'Select view',
     '             'View ref pt',
     '             'View plane',
     '             'View dist',
     '             'View up',
     '             'Zoom',
     '             'Exit'/
      DATA OPTION4/'Archive..',
     '             'Cancel..',
     '             'Define..',
     '             'Draw..',
     '             'Hide/Show..',
     '             'Label..',
     '             'Pick..',
     '             'Print..',
     '             'Recall..',
     '             'Transform..',
     '             'Exit'/
      DATA OPTION5/'Read signals..',
     '             'Setup',
     '             'Display signals',
     '             'Write signals',
     '             'Display field',
     '             'Pick electrode',
     '             'Cancel..',
     '             'Change colour..',
     '             'Exit'/
      DATA FPT/0.0,0.0,0.0/,FSHFT/0.0,0.0,0.0/,FSCALE/1.0,1.0,1.0/,
     '  FANGLE/0.0,0.0,0.0/

      DO i=1,100
        OP_STRING(i)(1:1)=CHAR(0)
      ENDDO
      NXCO=MXCO
      NXCOQU=MXCOQU
      GKS=.FALSE.
      GKS_WS_OPEN=.FALSE.
      MAPOPN=.FALSE.
      SMG_READ=.TRUE.
      NUM_LIBRARY=0
      TR01=.FALSE.
      TR02=.FALSE.
      TR03=.FALSE.
      TR04=.FALSE.
      FIRST =.TRUE.
      FIRSTF=.TRUE.
      FIRST_DRAW =.TRUE.
      FIRST_PLOT =.TRUE.
      FIRST_OXS  =.TRUE.
      FIRST_ZOOM =.TRUE.
      UPDATE_ZOOM=.FALSE.
      BUFFER_COUNT=0
      IREC_COMFILE=0
      CALL ENTERS('DIALOG',*9999)
      CALL GET_REV_TIME(ERROR,*9999)
      CALL GET_SYSTEM(ERROR,*9999)
      CALL GET_DISPLAY(ERROR,*9999)
      IOIP=1 !input  file
      IOOP=2 !usual output
      IODI=3 !diagnostics output
      IOTR=4 !trace output
      IOER=5 !error output
      IOH1=6 !help 1 output
      IOH2=7 !help 2 output
      IOH3=8 !help 3 output
      IOOUT=9 !output FILE
      CALL OPENF(1,'TERM','SYS$INPUT' ,'NEW',' ',' ',132,ERROR,*9999)
      CALL OPENF(2,'TERM','SYS$OUTPUT','NEW',' ',' ',255,ERROR,*9999)
      CALL OPENF(3,'TERM','SYS$OUTPUT','NEW',' ',' ',255,ERROR,*9999)
      CALL OPENF(4,'TERM','SYS$OUTPUT','NEW',' ',' ',255,ERROR,*9999)
      CALL OPENF(5,'TERM','SYS$OUTPUT','NEW',' ',' ',255,ERROR,*9999)
      CALL OPENF(6,'TERM','SYS$OUTPUT','NEW',' ',' ',255,ERROR,*9999)
      CALL OPENF(7,'TERM','SYS$OUTPUT','NEW',' ',' ',255,ERROR,*9999)
      CALL OPENF(8,'TERM','SYS$OUTPUT','NEW',' ',' ',255,ERROR,*9999)
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
      DOP=.FALSE.                     !<-|
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
      DO no_col=1,9
        NAME_COL(no_col)='UNDEFINED'
      ENDDO
      DO i=1,99
        IWKS(i)=0
      ENDDO
c     IWKDEF(0)=1
c     IWKDEF(1)=1
      DO nosg=1,MXSG
        ISEG(nosg)=0
        CSEG(nosg)=' '
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
      COD(3)=CFROMR(PI,'(E11.5)')
      CALL ASSIGN(STRING,1,COD,COQUD,ERROR,*9999)
!     COD(2)='pi'  !interferes with pick node command (PJH 10-feb-92)
!     COD(3)=CFROMR(PI,'(E11.5)')
!     CALL ASSIGN(STRING,1,COD,COQUD,ERROR,*9999)

      IF(USE_SOCKET.EQ..FALSE.) THEN

C *** Query mode
        ISTATUS = LIB$GETJPI(JPI$_MODE,,,,STATSTR,LEN)
        IF(.NOT.ISTATUS) CALL LIB$SIGNAL(%VAL(ISTATUS))

C ***   Check whether batch mode
        IF(STATSTR(1:1).EQ.'B') THEN !Cannot use SMG read routines
          SMG_READ=.FALSE.
        ELSE                        !Set up control-c trap
          CALL SET_CTRLC_AST
        ENDIF

C ***   Establish condition handler
        FATAL_HANDLER=.TRUE.
        CHANGE_HANDLER=.FALSE.
C CPB 6/10/93 Don't store old handler as OpenVMS cannot handle the typing of
C the procedure. Use LIB$REVERT instead.
        OLD_HANDLER=LIB$ESTABLISH(NEW_HANDLER) !keeps address of default handler

C ***   Check for cmiss.com file
        STRG='read cmiss;com'
        CALL PARSE(0,ISEG,CSEG,STRG,END,ERROR,*9999)
        
        CONTINUE=.TRUE.
        DO WHILE(CONTINUE) !is main program loop
        
          CALL GXWAIT(0.0,ERR)   ! Update graphics
			
C ***     Check for changes to first 4 windows
          IF(GKS.AND.GKS_WS_OPEN) THEN
            UPDATE=.TRUE.
            IF(MENU) THEN !go into mouse mode until 'Exit' is chosen
              CHOOSE=' '
              MENU=.FALSE.
            ELSE IF(.NOT.MENU) THEN !immed exit from mouse mode unless event
              CHOOSE='Exit'
            ENDIF
            DO WHILE(UPDATE)
              CHANGE=.FALSE.
              CALL EVENT(ID_WS,ID_DEVICE,ID_STATUS,CLASS,IDATA,R4DATA,
     '          SDATA,ERROR,*9999)
              IF(ID_WS.EQ.91) THEN
                IW=1
                CHAR3='1st'
              ELSE IF(ID_WS.EQ.92) THEN
                IW=2
                CHAR3='2nd'
              ELSE IF(ID_WS.EQ.93) THEN
                IW=3
                CHAR3='3rd'
              ELSE IF(ID_WS.EQ.94) THEN
                IW=4
                CHAR3='4th'
              ELSE IF(ID_WS.EQ.95) THEN
                IW=1
                CHAR3=' '
              ENDIF
              CIW=CFROMI(IW,'(I1)')
              IF(((ID_WS.GE.91.AND.ID_WS.LE.95).OR.ID_WS.EQ.7).AND.
     '          CLASS(1:6).EQ.'CHOICE') THEN
                Input_Choice=IDATA(1)
                IF(INPUT_CHOICE.GT.0) THEN
                  IF(ID_WS.EQ.91) THEN
                    CHOOSE=OPTION1(Input_Choice)
                  ELSE IF(ID_WS.EQ.92) THEN
                    CHOOSE=OPTION2(Input_Choice)
                  ELSE IF(ID_WS.EQ.93) THEN
                    CHOOSE=OPTION3(Input_Choice)
                  ELSE IF(ID_WS.EQ.94) THEN
                    CHOOSE=OPTION4(Input_Choice)
                  ELSE IF(ID_WS.EQ.95) THEN
                    CHOOSE=OPTION5(Input_Choice)
                  ELSE IF(ID_WS.EQ.7) THEN
                    IF(REFINE_ACTIVE) THEN
                      CHOOSE=OPTION10(Input_Choice)
                    ELSE IF(TRANSFORM_ACTIVE) THEN
                      CHOOSE=OPTION11(Input_Choice)
                    ENDIF
                  ENDIF
                  CALL TRIM(CHOOSE,IBEG,IEND)
                  IF((ID_WS.GE.91.AND.ID_WS.LE.94).OR.ID_WS.EQ.7) THEN
                    FORMAT=' >>'//CHOOSE(IBEG:IEND)//' the '//CHAR3
     '                //' window'
                    WRITE(OP_STRING,'(A)') FORMAT(1:30)
      		    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    IF(FIRST_ZOOM) THEN
                      UPDATE_ZOOM=.TRUE.
                      FIRST_ZOOM =.FALSE.
                    ENDIF
                    IF(CHOOSE(1:7).EQ.'Archive'.OR.CHOOSE(1:6).EQ.
     '                'Recall') THEN
                      CALL TRIM(FILE00,IBEG,IEND)
                      CALL GKS_STRG(56,NCHAR,'Enter file name ['
     '                  //FILE00(IBEG:IEND)//']',TEXT_STRING,ERROR,*100)
                      IF(NCHAR.GT.0) THEN
                        FILE00=TEXT_STRING
                      ENDIF
                    ENDIF
        
                    IF(CHOOSE(1:7).EQ.'Archive') THEN
                      IF(IW.EQ.3) THEN !Phigs
                        CALL ARCHIVE_GRAPHICS(FILE00,ERROR,*9999)
                      ELSE !GKS
                        STRG='FEM print;'//FILE00(IBEG:IEND)
     '                    //' metafile on '//CIW
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      ENDIF
                      WRITE(OP_STRING,'('' ...archiving complete'')')
      		      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                    ELSE IF(CHOOSE(1:6).EQ.'Cancel') THEN
                      NTCH=0
                      IF(CALL_BASE) THEN
                        OPTION(NTCH+1)='bases'
                        NTCH=NTCH+1
                      ENDIF
                      IF(CALL_DATA) THEN
                        OPTION(NTCH+1)='data'
                        NTCH=NTCH+1
                      ENDIF
                      IF(CALL_ELEM) THEN
                        OPTION(NTCH+1)='elements'
                        NTCH=NTCH+1
                      ENDIF
                      IF(CALL_NODE) THEN
                        OPTION(NTCH+1)='nodes'
                        NTCH=NTCH+1
                      ENDIF
                      OPTION(NTCH+1)='Exit'
                      NTCH=NTCH+1
c                     CALL COWK(IW,ERROR,*100)
c                     IF(NJT.EQ.2.AND.IW.EQ.1) UPDATE=.FALSE.
                      CALL PRECHOICE1(1,IW,NOCH,NOCO,NTCH,COD,
     '                  'REQUEST',OPTION,STRING,ERROR,*100)
                      CHOOSE=OPTION(NOCH)
        
                      IF(CHOOSE(1:5).EQ.'bases') THEN
                        STRG='FEM cancel bases'
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      ELSE IF(CHOOSE(1:4).EQ.'data') THEN
                        STRG='FEM cancel data'
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      ELSE IF(CHOOSE(1:8).EQ.'elements') THEN
                        STRG='FEM cancel elements'
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      ELSE IF(CHOOSE(1:5).EQ.'nodes') THEN
                        STRG='FEM cancel nodes'
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      ENDIF
        
                    ELSE IF(CHOOSE(1:6).EQ.'Change') THEN
                      OPTION(1)='colour'
                      OPTION(2)='line'
                      OPTION(3)='node'
                      OPTION(4)='exit'
                      NTCH=4
                      CALL PRECHOICE1(1,IW,NOCH,NOCO,NTCH,COD,
     '                  'REQUEST',OPTION,STRING,ERROR,*100)
                      CHOOSE=OPTION(NOCH)
                      CALL TRIM(CHOOSE,IBEG,IEND)
                      IF(CHOOSE(1:4).ne.'exit') THEN
c PJH   12AUG91         IF(CHOOSE(1:4).EQ.'line') THEN
c                         STRG='FEM update node deriv'
c                         CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
c                       ENDIF
                        STRG='FEM change '//CHOOSE(IBEG:IEND)
     '                    //';m on '//CIW
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      ELSE IF(CHOOSE(1:4).EQ.'Exit') THEN
                        CHOOSE=' '
                      ENDIF
        
                    ELSE IF(CHOOSE(1:7).EQ.'Check..') THEN
                      IF(CALL_NODE) THEN
                        OPTION(1)='nodes     *'
                      ELSE
                        OPTION(1)='nodes'
                      ENDIF
                      IF(CALL_BASE) THEN
                        OPTION(2)='bases     *'
                      ELSE
                        OPTION(2)='bases'
                      ENDIF
                      IF(CALL_ELEM) THEN
                        OPTION(3)='elements  *'
                      ELSE
                        OPTION(3)='elements'
                      ENDIF
                      IF(CALL_EQUA) THEN
                        OPTION(4)='equations *'
                      ELSE
                        OPTION(4)='equations'
                      ENDIF
                      IF(CALL_MATE) THEN
                        OPTION(5)='materials *'
                      ELSE
                        OPTION(5)='materials'
                      ENDIF
                      IF(CALL_INIT) THEN
                        OPTION(6)='initial   *'
                      ELSE
                        OPTION(6)='initial'
                      ENDIF
                      IF(CALL_SOLV) THEN
                        OPTION(7)='solve     *'
                      ELSE
                        OPTION(7)='solve'
                      ENDIF
                      OPTION(8)='exit'
                      NTCH=8
                      CALL PRECHOICE1(1,IW,NOCH,NOCO,NTCH,COD,
     '                  'REQUEST',OPTION,STRING,ERROR,*100)
                      CHOOSE=' '
        
                    ELSE IF(CHOOSE(1:6).EQ.'Define') THEN
                      NTCH=0
                      OPTION(NTCH+1)='bases'
                      OPTION(NTCH+2)='coordinates'
                      OPTION(NTCH+3)='data'
                      NTCH=NTCH+3
                      IF(CALL_BASE) THEN
                        OPTION(NTCH+1)='elements'
                        NTCH=NTCH+1
                      ENDIF
                      IF(CALL_NODE.AND.CALL_ELEM) THEN
                        OPTION(NTCH+1)='equation'
                        OPTION(NTCH+2)='field'
                        OPTION(NTCH+3)='fit'
                        OPTION(NTCH+4)='gauss'
                        NTCH=NTCH+4
                      ENDIF
                      IF(CALL_EQUA) THEN
                        OPTION(NTCH+1)='growth'
                        NTCH=NTCH+1
                      ENDIF
                      IF(CALL_EQUA) THEN
                        OPTION(NTCH+1)='initial'
                        NTCH=NTCH+1
                      ENDIF
                      IF(CALL_NODE.AND.CALL_ELEM) THEN
                        OPTION(NTCH+1)='lines'
                        NTCH=NTCH+1
                      ENDIF
                      IF(CALL_EQUA) THEN
                        OPTION(NTCH+1)='materials'
                        NTCH=NTCH+1
                      ENDIF
                      IF(CALL_NODE.AND.CALL_ELEM) THEN
                        OPTION(NTCH+1)='motion'
                        NTCH=NTCH+1
                      ENDIF
                      OPTION(NTCH+1)='nodes'
                      OPTION(NTCH+2)='object'
                      NTCH=NTCH+2
                      IF(CALL_EQUA) THEN
                        OPTION(NTCH+1)='solve'
                        NTCH=NTCH+1
                      ENDIF
                      OPTION(NTCH+1)='Exit'
                      NTCH=NTCH+1
                      CALL PRECHOICE1(1,IW,NOCH,NOCO,NTCH,COD,
     '                  'REQUEST',OPTION,STRING,ERROR,*100)
                      CHOOSE=OPTION(NOCH)
        
                      IF(CHOOSE(1:5).EQ.'bases') THEN
                        CALL PRECHOICE2(IW,NOCO,COD,'IPBASE',FILE_NAME,
     '                    PARAMETER_TYPE,Q,'mprw',STRING,ERROR,*100)
                        IF(Q.EQ.'m') THEN
                          STRG='FEM define bases;'//Q
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE IF(Q.EQ.'p'.OR.Q.EQ.'r') THEN
                          STRG='FEM define bases;'//Q//';'//FILE_NAME
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE IF(Q.EQ.'w') THEN
                          STRG='FEM define bases;w;'//FILE_NAME
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ENDIF
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:11).EQ.'coordinates') THEN
                        CALL PRECHOICE2(IW,NOCO,COD,'IPCOOR',FILE_NAME,
     '                    PARAMETER_TYPE,Q,'prw',STRING,ERROR,*100)
                        IF(Q.EQ.'p'.OR.Q.EQ.'r') THEN
                          STRG='FEM define coordinates;'//Q//';'
     '                      //FILE_NAME
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE IF(Q.EQ.'w') THEN
                          STRG='FEM define coordinates;w;'//FILE_NAME
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ENDIF
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:4).EQ.'data') THEN
                        CALL PRECHOICE2(IW,NOCO,COD,'IPDATA',FILE_NAME,
     '                    PARAMETER_TYPE,Q,'cmprw',STRING,ERROR,*100)
                        IF(Q.EQ.'c') THEN
                          STRG='FEM define data;c xi'
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE IF(Q.EQ.'m') THEN
                          STRG='FEM define data;m on '//CIW
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE IF(Q.EQ.'p'.OR.Q.EQ.'r') THEN
                          STRG='FEM define data;'//Q//';'//FILE_NAME
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          STRG='FEM draw data on '//CIW
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          STRG='FEM define data;c xi'
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE IF(Q.EQ.'w') THEN
                          STRG='FEM define data;w;'//FILE_NAME
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ENDIF
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:8).EQ.'elements') THEN
                        CALL PRECHOICE2(IW,NOCO,COD,'IPELEM',FILE_NAME,
     '                    PARAMETER_TYPE,Q,'mprw',STRING,ERROR,*100)
                        IF(Q.EQ.'m') THEN
                          STRG='FEM define elements;m on '//CIW
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          STRG='FEM draw lines on '//CIW
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE IF(Q.EQ.'p'.OR.Q.EQ.'r') THEN
                          STRG='FEM define elements;'//Q//';'//FILE_NAME
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          STRG='FEM draw elements on '//CIW
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE IF(Q.EQ.'w') THEN
                          STRG='FEM define elements;w;'//FILE_NAME
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ENDIF
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:8).EQ.'equation') THEN
                        CALL PRECHOICE2(IW,NOCO,COD,'IPEQUA',FILE_NAME,
     '                    PARAMETER_TYPE,Q,'prw',STRING,ERROR,*100)
                        IF(Q.EQ.'p'.OR.Q.EQ.'r') THEN
                          STRG='FEM define equation;'//Q//';'//FILE_NAME
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE IF(Q.EQ.'w') THEN
                          STRG='FEM define equation;w;'//FILE_NAME
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ENDIF
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:5).EQ.'field') THEN
                        CALL PRECHOICE2(IW,NOCO,COD,'IPFIEL',FILE_NAME,
     '                    PARAMETER_TYPE,Q,'dprw',STRING,ERROR,*100)
                        IF(Q.EQ.'d') THEN
                          STRG='FEM define field;d'
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE IF(Q.EQ.'p'.OR.Q.EQ.'r') THEN
                          STRG='FEM define field;'//Q//';'//FILE_NAME
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          STRG='FEM draw field on '//CIW
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE IF(Q.EQ.'w') THEN
                          STRG='FEM define field;w;'//FILE_NAME
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ENDIF
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:3).EQ.'fit') THEN
                        CALL PRECHOICE2(IW,NOCO,COD,'IPFIT',FILE_NAME,
     '                    PARAMETER_TYPE,Q,'dprw',STRING,ERROR,*100)
                        IF(Q.EQ.'d') THEN
                          STRG='FEM define fit;d '
     '                      //PARAMETER_TYPE(IBEG:IEND)
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE IF(Q.EQ.'p'.OR.Q.EQ.'r') THEN
                          STRG='FEM define fit;'//Q//';'//FILE_NAME
     '                      //PARAMETER_TYPE(IBEG:IEND)
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE IF(Q.EQ.'w') THEN
                          STRG='FEM define fit;w;'//FILE_NAME
     '                      //PARAMETER_TYPE(IBEG:IEND)
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ENDIF
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:6).EQ.'growth') THEN
                        CALL PRECHOICE2(IW,NOCO,COD,'IPGROW',FILE_NAME,
     '                    PARAMETER_TYPE,Q,'dprw',STRING,ERROR,*100)
                        IF(Q.EQ.'d') THEN
                          STRG='FEM define growth;d'
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE IF(Q.EQ.'p'.OR.Q.EQ.'r') THEN
                          STRG='FEM define growth;'//Q//';'//FILE_NAME
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE IF(Q.EQ.'w') THEN
                          STRG='FEM define growth;w;'//FILE_NAME
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ENDIF
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:7).EQ.'initial') THEN
                        CALL PRECHOICE2(IW,NOCO,COD,'IPINIT',FILE_NAME,
     '                    PARAMETER_TYPE,Q,'dprw',STRING,ERROR,*100)
                        IF(Q.EQ.'d') THEN
                          STRG='FEM define initial;d'
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          STRG='FEM draw nodes on '//CIW
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          STRG='FEM draw reaction on '//CIW
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE IF(Q.EQ.'p'.OR.Q.EQ.'r') THEN
                          STRG='FEM define initial;'//Q//';'//FILE_NAME
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          STRG='FEM draw nodes on '//CIW
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          STRG='FEM draw reaction on '//CIW
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE IF(Q.EQ.'w') THEN
                          STRG='FEM define initial;w;'//FILE_NAME
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ENDIF
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:5).EQ.'lines') THEN
                        CALL PRECHOICE2(IW,NOCO,COD,'IPLINE',FILE_NAME,
     '                    PARAMETER_TYPE,Q,'cmprw',STRING,ERROR,*100)
                        IF(Q.EQ.'c') THEN
                          STRG='FEM define lines;c'
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE IF(Q.EQ.'m') THEN
                          STRG='FEM define lines;m on '//CIW
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE IF(Q.EQ.'p'.OR.Q.EQ.'r') THEN
                          STRG='FEM define lines;'//Q//';'//FILE_NAME
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          STRG='FEM draw lines on '//CIW
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE IF(Q.EQ.'w') THEN
                          STRG='FEM define lines;w;'//FILE_NAME
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ENDIF
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:9).EQ.'materials') THEN
                        CALL PRECHOICE2(IW,NOCO,COD,'IPMATE',FILE_NAME,
     '                    PARAMETER_TYPE,Q,'dprw',STRING,ERROR,*100)
                        IF(Q.EQ.'d') THEN
                          STRG='FEM define materials;d'
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE IF(Q.EQ.'p'.OR.Q.EQ.'r') THEN
                          STRG='FEM define materials;'//Q//';'
     '                      //FILE_NAME
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE IF(Q.EQ.'w') THEN
                          STRG='FEM define materials;w;'//FILE_NAME
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ENDIF
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:6).EQ.'motion') THEN
                        CALL PRECHOICE2(IW,NOCO,COD,'IPMOTI',FILE_NAME,
     '                    PARAMETER_TYPE,Q,'dprw',STRING,ERROR,*100)
                        IF(Q.EQ.'d') THEN
                          STRG='FEM define motion;d '
     '                      //PARAMETER_TYPE(IBEG:IEND)
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE IF(Q.EQ.'p'.OR.Q.EQ.'r') THEN
                          STRG='FEM define motion;'//Q//';'//FILE_NAME
     '                      //PARAMETER_TYPE(IBEG:IEND)
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE IF(Q.EQ.'w') THEN
                          STRG='FEM define motion;w;'//FILE_NAME
     '                      //PARAMETER_TYPE(IBEG:IEND)
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ENDIF
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:5).EQ.'nodes') THEN
                        CALL PRECHOICE2(IW,NOCO,COD,'IPNODE',FILE_NAME,
     '                    PARAMETER_TYPE,Q,'mprw',STRING,ERROR,*100)
                        IF(Q.EQ.'m') THEN
                          STRG='FEM define nodes;m on '//CIW
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE IF(Q.EQ.'p'.OR.Q.EQ.'r') THEN
                          STRG='FEM define nodes;'//Q//';'//FILE_NAME
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          STRG='FEM define window;c'
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          STRG='FEM draw nodes on '//CIW
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE IF(Q.EQ.'w') THEN
                          STRG='FEM define nodes;w;'//FILE_NAME
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ENDIF
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:6).EQ.'object') THEN
                        STRG='FEM define object'
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        IF(CALL_NODE.AND.CALL_ELEM) THEN
                          STRG='FEM define data;c xi'
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ENDIF
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:5).EQ.'solve') THEN
                        CALL PRECHOICE2(IW,NOCO,COD,'IPSOLV',FILE_NAME,
     '                    PARAMETER_TYPE,Q,'dprw',STRING,ERROR,*100)
                        IF(Q.EQ.'d') THEN
                          STRG='FEM define solve;d'
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE IF(Q.EQ.'p'.OR.Q.EQ.'r') THEN
                          STRG='FEM define solve;'//Q//';'//FILE_NAME
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE IF(Q.EQ.'w') THEN
                          STRG='FEM define solve;w;'//FILE_NAME
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ENDIF
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:4).EQ.'Exit') THEN
                        CHOOSE=' '
                      ENDIF
        
                    ELSE IF(CHOOSE(1:7).EQ.'Display') THEN
                      OPTION(1)='history'
                      OPTION(2)='profile'
                      OPTION(3)='section'
                      OPTION(4)='Exit'
                      NTCH=4
                      CALL PRECHOICE1(1,IW,NOCH,NOCO,NTCH,COD,
     '                  'REQUEST',OPTION,STRING,ERROR,*100)
                      CHOOSE=OPTION(NOCH)
        
                      IF(CHOOSE(1:7).EQ.'history') THEN
                        STRG='FEM display history'
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:7).EQ.'profile') THEN
                        NTCH=0
                        IF(CALL_SOLV) THEN
                          OPTION(NTCH+1)='dependent'
                          NTCH=NTCH+1
                        ENDIF
                        IF(CALL_FIEL) THEN
                          OPTION(NTCH+1)='field'
                          NTCH=NTCH+1
                        ENDIF
                        IF(CALL_IMAG) THEN
                          OPTION(NTCH+1)='image'
                          NTCH=NTCH+1
                        ENDIF
                        IF(CALL_SOLV) THEN
                          OPTION(NTCH+1)='strain'
                          OPTION(NTCH+2)='cauchy'
                          OPTION(NTCH+3)='nominal'
                          OPTION(NTCH+4)='piola'
                          NTCH=NTCH+4
                        ENDIF
                        OPTION(NTCH+1)='Exit'
                        NTCH=NTCH+1
                        CALL PRECHOICE1(2,IW,NOCH,NOCO,NTCH,COD,
     '                    'REQUEST',OPTION,STRING,ERROR,*100)
                        CHOOSE=OPTION(NOCH)
                        IF(CHOOSE(1:5).EQ.'image') THEN
                          CALL GKS_DRAW(IW,ISEG,COD,CSEG,STRING,
     '                      ERROR,*100)
                          CHAR5=CFROMI(NTSG,'(I5)')
                          STRG='FEM define object with '//CHAR5
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          STRG='FEM display profile image'
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE IF(CHOOSE(1:4).EQ.'Exit') THEN
                          CHOOSE=' '
                        ELSE
                          CALL GKS_DRAW(IW,ISEG,COD,CSEG,STRING,
     '                      ERROR,*100)
                          CHAR5=CFROMI(NTSG,'(I5)')
                          STRG='FEM define object with '//CHAR5
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          STRG='FEM define data;c xi'
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          CALL TRIM(CHOOSE,IBEG,IEND)
                          STRG='FEM display profile '
     '                      //CHOOSE(IBEG:IEND)//' 1'
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ENDIF
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:7).EQ.'section') THEN
                        STRG='FEM display section'
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

                      ELSE IF(CHOOSE(1:4).EQ.'Exit') THEN
                        CHOOSE=' '
                      ENDIF
        
                    ELSE IF(CHOOSE(1:4).EQ.'Draw') THEN
                      NTCH=0
                      OPTION(NTCH+1)='alignment'
                      OPTION(NTCH+2)='axes'
                      NTCH=NTCH+2
                      IF(CALL_FIEL.OR.CALL_SOLV) THEN
                        OPTION(NTCH+1)='contour'
                        NTCH=NTCH+1
                      ENDIF
                      OPTION(NTCH+1)='data'
                      NTCH=NTCH+1
                      IF(CALL_BASE) THEN
                        OPTION(NTCH+1)='elements'
                        NTCH=NTCH+1
                      ENDIF
                      IF(CALL_NODE.AND.CALL_ELEM) THEN
                        OPTION(NTCH+1)='field'
                        OPTION(NTCH+2)='gauss'
                        NTCH=NTCH+2
                      ENDIF
                      IF(CALL_SOLV) THEN
                        OPTION(NTCH+1)='gradient'
                        NTCH=NTCH+1
                      ENDIF
                      IF(CALL_NODE.AND.CALL_ELEM) THEN
                        OPTION(NTCH+1)='lines'
                        NTCH=NTCH+1
                      ENDIF
                      OPTION(NTCH+1)='nodes'
                      NTCH=NTCH+1
                      IF(CALL_SOLV) THEN
                        OPTION(NTCH+1)='pressure'
                        NTCH=NTCH+1
                      ENDIF
                      IF(CALL_SOLV) THEN
                        OPTION(NTCH+1)='reaction'
                        NTCH=NTCH+1
                      ENDIF
                      OPTION(NTCH+1)='rule'
                      NTCH=NTCH+1
                      IF(CALL_SOLV) THEN
                        OPTION(NTCH+1)='strain'
                        OPTION(NTCH+2)='stress'
                        NTCH=NTCH+2
                      ENDIF
                      OPTION(NTCH+1)='Exit'
                      NTCH=NTCH+1
                      CALL PRECHOICE1(1,IW,NOCH,NOCO,NTCH,COD,
     '                  'REQUEST',OPTION,STRING,ERROR,*100)
                      CHOOSE=OPTION(NOCH)
        
                      IF(CHOOSE(1:9).EQ.'alignment') THEN
                        STRG='FEM draw alignment on '//CIW
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:4).EQ.'axes') THEN
                        STRG='FEM draw axes on '//CIW
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:7).EQ.'contour') THEN
                        STRG='FEM draw contour on '//CIW
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:4).EQ.'data') THEN
                        OPTION( 1)='geometry'
                        OPTION( 2)='fibre'
                        OPTION( 3)='field'
                        OPTION( 4)='numbers'
                        OPTION( 5)='projections'
                        OPTION( 6)='values'
                        OPTION( 7)='trace'
                        OPTION( 8)='Exit'
                        NTCH=8
                        CALL PRECHOICE1(2,IW,NOCH,NOCO,NTCH,COD,
     '                    'REQUEST',OPTION,STRING,ERROR,*100)
                        CHOOSE=OPTION(NOCH)
                        CALL TRIM(CHOOSE,IBEG,IEND)
                        STRG='FEM draw data '//CHOOSE(IBEG:IEND)
     '                    //' on '//CIW
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:8).EQ.'elements') THEN
                        STRG='FEM draw elements on '//CIW
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:5).EQ.'field') THEN
                        STRG='FEM draw field on '//CIW
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:5).EQ.'gauss') THEN
                        STRG='FEM draw gauss on '//CIW
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:8).EQ.'gradient') THEN
                        STRG='FEM draw gradient on '//CIW
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:5).EQ.'lines') THEN
                        STRG='FEM draw lines on '//CIW
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:5).EQ.'nodes') THEN
                        STRG='FEM draw nodes on '//CIW
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:8).EQ.'pressure') THEN
                        STRG='FEM update Gauss pressure from 1'
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        STRG='FEM define field;d'
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        STRG='FEM define fit;d Gauss'
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        STRG='FEM fit Gauss'
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        STRG='FEM draw field on '//CIW
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:8).EQ.'reaction') THEN
                        STRG='FEM draw reaction on '//CIW
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:4).EQ.'rule') THEN
                        STRG='FEM draw rule on '//CIW
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:6).EQ.'strain') THEN
                        OPTION(1)='field'
                        OPTION(2)='vectors'
                        OPTION(3)='Exit'
                        NTCH=3
                        CALL PRECHOICE1(2,IW,NOCH,NOCO,NTCH,COD,
     '                    'REQUEST',OPTION,STRING,ERROR,*100)
                        CHOOSE=OPTION(NOCH)
                        IF(CHOOSE(1:5).eq.'field') THEN
                          STRG='FEM update Gauss strain'
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          STRG='FEM define field;d'
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          STRG='FEM define fit;d Gauss'
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          STRG='FEM fit Gauss'
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          STRG='FEM draw field on '//CIW
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          WRITE(OP_STRING,'('' ...completed'')')
      			  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                        ELSE IF(CHOOSE(1:7).eq.'vectors') THEN
                          STRG='FEM draw strain on '//CIW
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          WRITE(OP_STRING,'('' ...completed'')')
      			  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                        ELSE IF(CHOOSE(1:4).EQ.'Exit') THEN
                          CHOOSE=' '
                        ENDIF
        
                      ELSE IF(CHOOSE(1:8).EQ.'stress') THEN
                        OPTION(1)='field'
                        OPTION(2)='vectors'
                        OPTION(3)='Exit'
                        NTCH=3
                        CALL PRECHOICE1(2,IW,NOCH,NOCO,NTCH,COD,
     '                    'REQUEST',OPTION,STRING,ERROR,*100)
                        CHOOSE=OPTION(NOCH)
                        IF(CHOOSE(1:5).eq.'field') THEN
                          STRG='FEM update Gauss stress'
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          STRG='FEM define field;d'
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          STRG='FEM define fit;d Gauss'
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          STRG='FEM fit Gauss'
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          STRG='FEM draw field on '//CIW
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          WRITE(OP_STRING,'('' ...completed'')')
      			  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                        ELSE IF(CHOOSE(1:7).eq.'vectors') THEN
                          STRG='FEM draw stress on '//CIW
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          WRITE(OP_STRING,'('' ...completed'')')
      			  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                        ELSE IF(CHOOSE(1:4).EQ.'Exit') THEN
                          CHOOSE=' '
                        ENDIF
        
                      ELSE IF(CHOOSE(1:4).EQ.'Exit') THEN
                        CHOOSE=' '
                      ENDIF

                    ELSE IF(CHOOSE(1:9).EQ.'Hide/Show') THEN
                      SHOW_HIDE=.TRUE.
                      IWINDOW=IW
                      DO WHILE(SHOW_HIDE) !i.e. until Exit is chosen
                        COD(1)='FEM'
                        NOCO=1
C CPB 20/7/93 Commenting out fem call because cannot parse work arrays
C through from syntax.
C                        CALL FEM(ISEG,NOCO,NTCOD,NTCOQUD,COD,COQUD,
C     '                    CSEG,END,STRING,ERROR,*100)
                      ENDDO
        
                    ELSE IF(CHOOSE(1:5).EQ.'Label') THEN
                      CALL GKS_DRAW(IW,ISEG,COD,CSEG,STRING,ERROR,*100)
        
                    ELSE IF(CHOOSE(1:4).EQ.'List') THEN
                      NTCH=0
                      IF(CALL_BASE) THEN
                        OPTION(NTCH+1)='bases'
                        NTCH=NTCH+1
                      ENDIF
                      OPTION(NTCH+1)='coordinates'
                      NTCH=NTCH+1
                      IF(CALL_DATA) THEN
                        OPTION(NTCH+3)='data'
                        NTCH=NTCH+1
                      ENDIF
                      IF(CALL_ELEM) THEN
                        OPTION(NTCH+1)='elements'
                        NTCH=NTCH+1
                      ENDIF
                      IF(CALL_EQUA) THEN
                        OPTION(NTCH+1)='equation'
                        NTCH=NTCH+1
                      ENDIF
                      IF(CALL_FIT) THEN
                        OPTION(NTCH+1)='fit'
                        NTCH=NTCH+1
                      ENDIF
                      IF(CALL_ELEM) THEN
                        OPTION(NTCH+1)='gauss'
                        NTCH=NTCH+1
                      ENDIF
                      IF(CALL_GROW) THEN
                        OPTION(NTCH+1)='growth'
                        NTCH=NTCH+1
                      ENDIF
                      IF(CALL_INIT) THEN
                        OPTION(NTCH+1)='initial'
                        NTCH=NTCH+1
                      ENDIF
                      IF(CALL_LINE) THEN
                        OPTION(NTCH+1)='lines'
                        NTCH=NTCH+1
                      ENDIF
                      IF(CALL_MATE) THEN
                        OPTION(NTCH+1)='materials'
                        NTCH=NTCH+1
                      ENDIF
                      IF(CALL_MOTI) THEN
                        OPTION(NTCH+1)='motion'
                        NTCH=NTCH+1
                      ENDIF
                      IF(CALL_NODE) THEN
                        OPTION(NTCH+1)='nodes'
                        NTCH=NTCH+1
                      ENDIF
                      IF(CALL_OBJE) THEN
                        OPTION(NTCH+1)='object'
                        NTCH=NTCH+1
                      ENDIF
                      IF(CALL_SOLV) THEN
                        OPTION(NTCH+1)='solve'
                        OPTION(NTCH+2)='strain'
                        OPTION(NTCH+3)='stress'
                        NTCH=NTCH+3
                      ENDIF
                      OPTION(NTCH+1)='Exit'
                      NTCH=NTCH+1
                      CALL PRECHOICE1(1,IW,NOCH,NOCO,NTCH,COD,
     '                    'REQUEST',OPTION,STRING,ERROR,*100)
                      CHOOSE=OPTION(NOCH)
                      CALL TRIM(CHOOSE,IBEG,IEND)
                      IF(CHOOSE(1:4).ne.'Exit') THEN
                        IF(CALL_SOLV.AND.CHOOSE(1:5).EQ.'nodes') THEN
                          OPTION(1)='cartesian'
                          OPTION(2)='displacement'
                          OPTION(3)='flux'
                          OPTION(4)='reaction'
                          OPTION(5)='solution'
                          OPTION(6)='Exit'
                          NTCH=6
                          CALL PRECHOICE1(2,IW,NOCH,NOCO,NTCH,COD,
     '                      'REQUEST',OPTION,STRING,ERROR,*100)
                          STRG='FEM list '//CHOOSE(IBEG:IEND)//' '
     '                      //OPTION(NOCH)
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ELSE
                          STRG='FEM list '//CHOOSE(IBEG:IEND)
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ENDIF
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      ELSE IF(CHOOSE(1:4).EQ.'Exit') THEN
                        CHOOSE=' '
                      ENDIF
        
                    ELSE IF(CHOOSE(1:3).EQ.'Pan'.OR
     '                     .CHOOSE(1:9).EQ.'Translate'.OR
     '                     .CHOOSE(1:11).EQ.'View ref pt'.OR
     '                     .CHOOSE(1:10).EQ.'View plane'.OR
     '                     .CHOOSE(1:7).EQ.'View up') THEN
                      TRANSFORM_ACTIVE=.TRUE.
                      TRANSFORM_TYPE='WINDOW'
                      IF(NJT.LE.2) THEN
                        IF(IW.EQ.1) THEN
                          IF(.NOT.DISPLAY_VALUATOR_83) THEN
                            DISPLAY_VALUATOR_83=.TRUE.
                            CALL VALUATOR('1',83,'EVENT',1,-1.,1.,0.,
     '                        VALUE,0.18*XDISP,
     '                        YDISP-0.49*XDISP-0.110*YDISP,ERROR,*100)
                          ENDIF
                          IF(.NOT.DISPLAY_VALUATOR_84) THEN
                            DISPLAY_VALUATOR_84=.TRUE.
                            CALL VALUATOR('2',84,'EVENT',2,-1.,1.,0.,
     '                        VALUE,
     '                        0.49*XDISP,YDISP-0.49*XDISP-0.110*YDISP,
     '                        ERROR,*100)
                          ENDIF
                        ENDIF
                      ELSE IF(NJT.EQ.3) THEN
                        IF(IW.EQ.1) THEN
                          IF(.NOT.DISPLAY_VALUATOR_84) THEN
                            DISPLAY_VALUATOR_84=.TRUE.
                            CALL VALUATOR('2',84,'EVENT',9,-1.,1.,0.,
     '                        VALUE,0.01*DISP,0.47*DISP,ERROR,*100)
                          ENDIF
                          IF(.NOT.DISPLAY_VALUATOR_83) THEN
                            DISPLAY_VALUATOR_83=.TRUE.
                            CALL VALUATOR('1',83,'EVENT',9,-1.,1.,0.,
     '                        VALUE,0.0,0.48*DISP,ERROR,*100)
                          ENDIF
                        ELSE IF(IW.EQ.2) THEN
                          IF(.NOT.DISPLAY_VALUATOR_86) THEN
                            DISPLAY_VALUATOR_86=.TRUE.
                            CALL VALUATOR('2',86,'EVENT',9,-1.,1.,0.,
     '                        VALUE,1.00*DISP,0.97*DISP,ERROR,*100)
                          ENDIF
                          IF(.NOT.DISPLAY_VALUATOR_85) THEN
                            DISPLAY_VALUATOR_85=.TRUE.
                            CALL VALUATOR('1',85,'EVENT',9,-1.,1.,0.,
     '                        VALUE,0.99*DISP,0.98*DISP,ERROR,*100)
                          ENDIF
                        ELSE IF(IW.EQ.3) THEN
                          IF(.NOT.DISPLAY_VALUATOR_89) THEN
                            DISPLAY_VALUATOR_89=.TRUE.
                            CALL VALUATOR('3',89,'EVENT',9,-1.,1.,0.,
     '                        VALUE,1.01*DISP,0.45*DISP,ERROR,*100)
                          ENDIF
                          IF(.NOT.DISPLAY_VALUATOR_88) THEN
                            DISPLAY_VALUATOR_88=.TRUE.
                            CALL VALUATOR('2',88,'EVENT',9,-1.,1.,0.,
     '                        VALUE,1.00*DISP,0.46*DISP,ERROR,*100)
                          ENDIF
                          IF(.NOT.DISPLAY_VALUATOR_87) THEN
                            DISPLAY_VALUATOR_87=.TRUE.
                            CALL VALUATOR('1',87,'EVENT',9,-1.,1.,0.,
     '                        VALUE,0.99*DISP,0.47*DISP,ERROR,*100)
                          ENDIF
                        ELSE IF(IW.EQ.4) THEN
                          IF(.NOT.DISPLAY_VALUATOR_82) THEN
                            DISPLAY_VALUATOR_82=.TRUE.
                            CALL VALUATOR('2',82,'EVENT',9,-1.,1.,0.,
     '                        VALUE,0.11*DISP,0.16*DISP,ERROR,*100)
                          ENDIF
                          IF(.NOT.DISPLAY_VALUATOR_81) THEN
                            DISPLAY_VALUATOR_81=.TRUE.
                            CALL VALUATOR('1',81,'EVENT',9,-1.,1.,0.,
     '                        VALUE,0.10*DISP,0.17*DISP,ERROR,*100)
                          ENDIF
                        ENDIF
                      ENDIF
        
                    ELSE IF(CHOOSE(1:8).EQ.'Parallel') THEN
                      MODE_PROJ=0  !PHIGS$K_PARALLEL
                      CHANGE=.TRUE.
        
                    ELSE IF(CHOOSE(1:11).EQ.'Perspective') THEN
                      MODE_PROJ=1  !PHIGS$K_PERSPECTIVE
                      CHANGE=.TRUE.
        
                    ELSE IF(CHOOSE(1:4).EQ.'Pick') THEN
                      OPTION( 1)='data'
                      OPTION( 2)='elements'
                      OPTION( 3)='fibres'
                      OPTION( 4)='lines'
                      OPTION( 5)='nodes'
                      OPTION( 6)='exit'
                      NTCH=6
                      CALL PRECHOICE1(1,IW,NOCH,NOCO,NTCH,COD,
     '                  'REQUEST',OPTION,STRING,ERROR,*100)
                      CHOOSE=OPTION(NOCH)
                      CALL TRIM(CHOOSE,IBEG,IEND)
                      IF(CHOOSE(1:4).ne.'exit') THEN
                        STRG='FEM pick '//CHOOSE(IBEG:IEND)//' on '
     '                    //CFROMI(IW,'(I1)')
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      ELSE IF(CHOOSE(1:4).EQ.'Exit') THEN
                        CHOOSE=' '
                      ENDIF
        
                    ELSE IF(CHOOSE(1:5).EQ.'Print') THEN
                      IWINDOW=IW
                      OPTION(1)='Print to file..'
                      OPTION(2)='Print now'
                      OPTION(3)='Exit'
                      NTCH=3
                      CALL PRECHOICE1(1,IW,NOCH,NOCO,NTCH,COD,
     '                  'REQUEST',OPTION,STRING,ERROR,*100)
                      CHOOSE=OPTION(NOCH)
        
                      IF(CHOOSE(1:13).EQ.'Print to file') THEN
                        CALL TRIM(FILE00,IBEG,IEND)
                        CALL GKS_STRG(56,NCHAR,'Enter file name ['
     '                    //FILE00(IBEG:IEND)//']',TEXT_STRING,ERROR,
     '                    *100)
                        IF(NCHAR.GT.0) THEN
                          FILE00=TEXT_STRING
                        ENDIF
                        IF(IWKT(IW).EQ.1) THEN      !GKS
                          STRG='FEM define window on 15'
                        ELSE IF(IWKT(IW).EQ.2) THEN !PHIGS
                          STRG='FEM define window on 16'
                        ENDIF
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        CALL TRIM(FILE00,IBEG,IEND)
                        STRG='FEM print;'//FILE00(IBEG:IEND)
     '                    //' postscript on '//CFROMI(IW,'(I1)')
                        WRITE(OP_STRING,'('' >'',A)') STRG(1:80)
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        STRG='close postscript'
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        WRITE(OP_STRING,'('' ..printing complete'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:9).EQ.'Print now') THEN
                        IF(IWKT(IW).EQ.1) THEN      !GKS
                          STRG='FEM define window on 15'
                        ELSE IF(IWKT(IW).EQ.2) THEN !PHIGS
                          STRG='FEM define window on 16'
                        ENDIF
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        CALL TRIM(FILE00,IBEG,IEND)
                        STRG='FEM print;'//FILE00(IBEG:IEND)
     '                    //' postscript on '//CFROMI(IW,'(I1)')
                        WRITE(OP_STRING,'('' >'',A)') STRG(1:80)
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        STRG='close postscript'
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        ISTATUS=LIB$SPAWN('post '//FILE00(IBEG:IEND)
     '                    //'.ps',,,,,,,,,,,)
                        ISTATUS=LIB$SPAWN('purge '//FILE00(IBEG:IEND)
     '                    //'.ps',,,,,,,,,,,)
                        WRITE(OP_STRING,'('' ...'',A,'
     '                    //'''.ps submitted to print queue'')')
     '                    FILE00(IBEG:IEND)
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:4).EQ.'Exit') THEN
                        CHOOSE=' '
                      ENDIF
        
                    ELSE IF(CHOOSE(1:6).EQ.'Read..') THEN
                      OPTION(1)='com file'
                      OPTION(2)='iod file'
                      OPTION(3)='geometry'
                      OPTION(4)='problem'
                      OPTION(5)='Exit'
                      NTCH=5
                      CALL PRECHOICE1(1,IW,NOCH,NOCO,NTCH,COD,
     '                  'REQUEST',OPTION,STRING,ERROR,*100)
                      CHOOSE=OPTION(NOCH)
        
                      IF(CHOOSE(1:8).EQ.'com file') THEN
                        CALL DISPLAY_FILE(IW,NOCO,NTFILE,CO,'com',
     '                    FILE_NAME,STRING,ERROR,*100)
                        IF(FILE_NAME(1:4).ne.'Exit') THEN
                          CALL TRIM(FILE_NAME,IBEG,IEND)
                          STRG='FEM read '//FILE_NAME(IBEG:IEND)//';com'
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          WRITE(OP_STRING,'('' ...completed'')')
      			  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                        ELSE IF(FILE_NAME(1:4).EQ.'Exit') THEN
                          CHOOSE=' '
                        ENDIF
        
                      ELSE IF(CHOOSE(1:8).EQ.'iod file') THEN
                        CALL DISPLAY_FILE(IW,NOCO,NTFILE,CO,'iod',
     '                    FILE_NAME,STRING,ERROR,*100)
                        IF(FILE_NAME(1:4).ne.'Exit') THEN
                          CALL TRIM(FILE_NAME,IBEG,IEND)
                          STRG='FEM read '//FILE_NAME(IBEG:IEND)//';iod'
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          IF(NJT.EQ.3) CIW='3'
                          STRG='FEM define window;c'
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          STRG='FEM draw axes on '//CIW
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          STRG='FEM draw lines on '//CIW
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          IF(NJT.EQ.2) THEN
                            STRG='FEM draw nodes on '//CIW
                            CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                            STRG='FEM draw elements on '//CIW
                            CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          ENDIF
                          WRITE(OP_STRING,'('' ...completed'')')
      			  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                        ELSE IF(FILE_NAME(1:4).EQ.'Exit') THEN
                          CHOOSE=' '
                        ENDIF
        
                      ELSE IF(CHOOSE(1:8).EQ.'geometry') THEN
                        OPTION(1)='cmiss file'
                        OPTION(2)='user file'
                        OPTION(3)='Exit'
                        NTCH=3
                        CALL PRECHOICE1(2,IW,NOCH,NOCO,NTCH,COD,
     '                    'REQUEST',OPTION,STRING,ERROR,*100)
                        CHOOSE=OPTION(NOCH)
                        IF(CHOOSE(1:10).EQ.'cmiss file') THEN
                          OPTION(1)='1D/rect.cart./linear'
                          OPTION(2)='2D/rect.cart./bilinear'
                          OPTION(3)='2D/cyl.polar /bilinear'
                          OPTION(4)='3D/rect.cart./trilinear'
                          OPTION(5)='3D/cyl.polar /trilinear'
                          OPTION(6)='3D/prolate   /trilinear'
                          OPTION(7)='3D/Heart'
                          OPTION(8)='Exit'
                          NTCH=8
                          CALL PRECHOICE1(2,IW,NOCH,NOCO,NTCH,COD,
     '                      'REQUEST',OPTION,STRING,ERROR,*100)
                          CHOOSE=OPTION(NOCH)
                IF(CHOOSE(1:4).ne.'Exit') THEN
                  IF(CHOOSE(1:20).EQ.'1D/rect.cart./linear') THEN
                    FILE_NAME='1D_RC_Linear'
                  ELSE IF(CHOOSE(1:22).EQ.'2D/rect.cart./bilinear') THEN
                    FILE_NAME='2D_RC_Bilinear'
                  ELSE IF(CHOOSE(1:22).EQ.'2D/cyl.polar /bilinear')
     '              THEN
                    FILE_NAME='2D_CP_Bilinear'
                  ELSE IF(CHOOSE(1:23).EQ.'3D/rect.cart./trilinear')
     '              THEN
                    FILE_NAME='3D_RC_Trilinear'
                  ELSE IF(CHOOSE(1:23).EQ.'3D/cyl.polar /trilinear')
     '              THEN
                    FILE_NAME='3D_CP_Trilinear'
                  ELSE IF(CHOOSE(1:23).EQ.'3D/prolate   /trilinear')
     '              THEN
                    FILE_NAME='3D_PS_Trilinear'
                  ELSE IF(CHOOSE(1:23).EQ.'3D/Heart') THEN
                    FILE_NAME='H27_F_S0'
                  ENDIF
                  CALL TRIM(FILE_NAME,IBEG,IEND)
                  STRG='FEM read '//FILE_NAME(IBEG:IEND)//';iod;doc'
                  CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                ELSE IF(CHOOSE(1:4).EQ.'Exit') THEN
                  FILE_NAME='Exit'
                ENDIF
        
                        ELSE IF(CHOOSE(1:9).EQ.'user file') THEN
                          CALL DISPLAY_FILE(IW,NOCO,NTFILE,CO,'ipbase',
     '                      FILE_NAME,STRING,ERROR,*100)
                          IF(FILE_NAME(1:4).ne.'Exit') THEN
                            STRG='FEM define nodes;r;'//FILE_NAME
                            CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                            STRG='FEM define bases;r;'//FILE_NAME
                            CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                            STRG='FEM define elements;r;'//FILE_NAME
                            CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          ENDIF
                        ENDIF
        
                        IF(FILE_NAME(1:4).ne.'Exit') THEN
                          IF(NJT.EQ.3) CIW='3'
                          STRG='FEM define window;c'
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          STRG='FEM draw axes on '//CIW
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          STRG='FEM draw lines on '//CIW
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          IF(NJT.EQ.2) THEN
                            STRG='FEM draw nodes on '//CIW
                            CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                            STRG='FEM draw elements on '//CIW
                            CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          ENDIF
                          WRITE(OP_STRING,'('' ...completed'')')
      			  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                        ELSE IF(FILE_NAME(1:4).EQ.'Exit') THEN
                          CHOOSE=' '
                        ENDIF
        
                      ELSE IF(CHOOSE(1:7).EQ.'problem') THEN
                        OPTION(1)='cmiss file'
                        OPTION(2)='user file'
                        OPTION(3)='Exit'
                        NTCH=3
                        CALL PRECHOICE1(2,IW,NOCH,NOCO,NTCH,COD,
     '                    'REQUEST',OPTION,STRING,ERROR,*100)
                        CHOOSE=OPTION(NOCH)
                        IF(CHOOSE(1:10).EQ.'cmiss file') THEN
                          OPTION(1)='1D/rc/Transient heatflow/linear'
                          OPTION(2)='2D/rc/Heat flow/bilinear'
                          OPTION(3)=' '
                          OPTION(4)=' '
                          OPTION(5)=' '
                          OPTION(6)=' '
                          OPTION(7)='Exit'
                          NTCH=7
                          CALL PRECHOICE1(2,IW,NOCH,NOCO,NTCH,COD,
     '                      'REQUEST',OPTION,STRING,ERROR,*100)
                          CHOOSE=OPTION(NOCH)
                IF(CHOOSE(1:4).ne.'Exit') THEN
                  IF(CHOOSE(1:31).EQ.'1D/rc/Transient heatflow/linear')
     '              THEN
                    FILE_NAME='1D_RC_Linear_Transient_Heat'
                  ELSE IF(CHOOSE(1:24).EQ.'2D/rc/bilinear/Heat-flow')
     '              THEN
                    FILE_NAME='2D_RC_Bilinear_Heat-flow'
                  ENDIF
                  CALL TRIM(FILE_NAME,IBEG,IEND)
                  STRG='FEM read '//FILE_NAME(IBEG:IEND)//';iod;doc'
                  CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                ELSE IF(CHOOSE(1:4).EQ.'Exit') THEN
                  FILE_NAME='Exit'
                ENDIF
        
                        ELSE IF(CHOOSE(1:9).EQ.'user file') THEN
                          CALL DISPLAY_FILE(IW,NOCO,NTFILE,CO,
     '                      'ipequa',FILE_NAME,STRING,ERROR,*100)
                          IF(FILE_NAME(1:4).ne.'Exit') THEN
                            STRG='FEM define equation;r;'//FILE_NAME
                            CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                            STRG='FEM define materials;r;'//FILE_NAME
                            CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                            STRG='FEM define initial;r;'//FILE_NAME
                            CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                            STRG='FEM define solve;r;'//FILE_NAME
                            CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          ENDIF
                        ENDIF
        
                        IF(FILE_NAME(1:4).ne.'Exit') THEN
                          STRG='FEM draw nodes on '//CIW
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          STRG='FEM draw reaction on '//CIW
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          WRITE(OP_STRING,'('' ...completed'')')
      			  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                        ELSE IF(FILE_NAME(1:4).EQ.'Exit') THEN
                          CHOOSE=' '
                        ENDIF
        
                      ELSE IF(CHOOSE(1:4).EQ.'Exit') THEN
                        CHOOSE=' '
                      ENDIF
        
                    ELSE IF(CHOOSE(1:8).EQ.'Recall..') THEN
                      IF(IW.EQ.3) THEN !Phigs
                        CALL RECALL_GRAPHICS(IW,FILE_NAME,ERROR,*9999)
                      ELSE !GKS
                      ENDIF
                      WRITE(OP_STRING,'('' ...recalling complete'')')
      		      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                    ELSE IF(CHOOSE(1:8).EQ.'Refine..') THEN
                      OPTION10(1)='Refine Xi 1'
                      OPTION10(2)='Refine Xi 2'
                      OPTION10(3)='Refine Xi 3'
                      OPTION10(4)='Exit'
                      NTCH=4
                      CALL PRECHOICE1(1,IW,NOCH,NOCO,NTCH,COD,'EVENT',
     '                  OPTION10,STRING,ERROR,*100)
                      REFINE_ACTIVE=.TRUE.
        
                    ELSE IF(CHOOSE(1:9).EQ.'Refine Xi') THEN
                      STRG='FEM refine Xi '//CHOOSE(11:11)
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      STRG='FEM update mesh'
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      WRITE(OP_STRING,'('' ...completed'')')
      		      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                    ELSE IF(CHOOSE(1:11).EQ.'Rotate data'.OR
     '                     .CHOOSE(1:11).EQ.'Proj ref pt'.OR
     '                     .CHOOSE(1:11).EQ.'Rotate view') THEN
                      TRANSFORM_ACTIVE=.TRUE.
                      TRANSFORM_TYPE='WINDOW'
                      IF(IW.EQ.3) THEN
                        IF(.NOT.DISPLAY_VALUATOR_89) THEN
                          DISPLAY_VALUATOR_89=.TRUE.
                          CALL VALUATOR('3',89,'EVENT',9,-1.,1.,0.,
     '                      VALUE,1.01*DISP,0.45*DISP,ERROR,*100)
                        ENDIF
                        IF(.NOT.DISPLAY_VALUATOR_88) THEN
                          DISPLAY_VALUATOR_88=.TRUE.
                          CALL VALUATOR('2',88,'EVENT',9,-1.,1.,0.,
     '                      VALUE,1.00*DISP,0.46*DISP,ERROR,*100)
                        ENDIF
                        IF(.NOT.DISPLAY_VALUATOR_87) THEN
                          DISPLAY_VALUATOR_87=.TRUE.
                          CALL VALUATOR('1',87,'EVENT',9,-1.,1.,0.,
     '                      VALUE,0.99*DISP,0.47*DISP,ERROR,*100)
                        ENDIF
                      ENDIF
        
                    ELSE IF(CHOOSE(1:7).EQ.'Rescale'.OR
     '                     .CHOOSE(1:6).EQ.'Rotate'.OR
     '                     .CHOOSE(1:5).EQ.'Scale'.OR
     '                     .CHOOSE(1:12).EQ.'Zoom (scale)'.OR
     '                     .CHOOSE(1:4).EQ.'Zoom'.AND.IW.EQ.3.OR
     '                     .CHOOSE(1:10).EQ.'View plane'.OR
     '                     .CHOOSE(1:11).EQ.'Front plane') THEN
                      TRANSFORM_ACTIVE=.TRUE.
                      TRANSFORM_TYPE='WINDOW'
                      IF(NJT.LE.2) THEN
                        IF(IW.EQ.1) THEN
                          IF(.NOT.DISPLAY_VALUATOR_83) THEN
                            DISPLAY_VALUATOR_83=.TRUE.
                            CALL VALUATOR('1',83,'EVENT',1,-1.,1.,0.,
     '                        VALUE,0.18*XDISP,
     '                        YDISP-0.49*XDISP-0.110*YDISP,ERROR,*100)
                          ENDIF
                        ENDIF
                      ELSE IF(NJT.EQ.3) THEN
                        IF(IW.EQ.1) THEN
                          IF(.NOT.DISPLAY_VALUATOR_83) THEN
                            DISPLAY_VALUATOR_83=.TRUE.
                            CALL VALUATOR('1',83,'EVENT',9,-1.,1.,0.,
     '                        VALUE,0.0,0.48*DISP,ERROR,*100)
                          ENDIF
                        ELSE IF(IW.EQ.2) THEN
                          IF(.NOT.DISPLAY_VALUATOR_85) THEN
                            DISPLAY_VALUATOR_85=.TRUE.
                            CALL VALUATOR('1',85,'EVENT',9,-1.,1.,0.,
     '                        VALUE,0.99*DISP,0.98*DISP,ERROR,*100)
                          ENDIF
                        ELSE IF(IW.EQ.3) THEN
                          IF(.NOT.DISPLAY_VALUATOR_87) THEN
                            DISPLAY_VALUATOR_87=.TRUE.
                            CALL VALUATOR('1',87,'EVENT',9,-1.,1.,0.,
     '                        VALUE,0.99*DISP,0.47*DISP,ERROR,*100)
                          ENDIF
                        ELSE IF(IW.EQ.4) THEN
                          IF(.NOT.DISPLAY_VALUATOR_81) THEN
                            DISPLAY_VALUATOR_81=.TRUE.
                            CALL VALUATOR('1',81,'EVENT',9,-1.,1.,0.,
     '                        VALUE,0.10*DISP,0.17*DISP,ERROR,*100)
                          ENDIF
                        ENDIF
                      ENDIF
        
                    ELSE IF(CHOOSE(1:5).EQ.'Reset') THEN
                      CALL ACWK(IW,1,ERROR,*100)
                      CALL WKST_WINDOW(IW,XNDC(1,1,IW),XNDC(1,2,IW),
     '                  XNDC(1,3,IW),XNDC(1,4,IW),ERROR,*100)
                      CALL DAWK(IW,1,ERROR,*100)
                      IZOOM(IW)=1
                      FIRST_ZOOM=.TRUE.
                      CHANGE=.TRUE.
        
                      IF(IW.EQ.3) THEN !is this needed AAY?
                        CHANGE=.TRUE.
                        VIEW_PLANE_DIST_NEW =VIEW_PLANE_DIST
                        BACK_PLANE_DIST_NEW =BACK_PLANE_DIST
                        FRONT_PLANE_DIST_NEW=FRONT_PLANE_DIST
                        MODE_PROJ=0  !PHIGS$K_PARALLEL
                        DO nj=1,3
                          ANGLE(nj)=0.0
                          SCALE(nj)=1.0
                          SHIFT(nj)=0.0
                          PROJ_REF_PT_NEW(nj)=PROJ_REF_PT(nj)
                          VIEW_REF_PT_NEW(nj)=VIEW_REF_PT(nj)
                          VIEW_PLANE_NEW(nj)=VIEW_PLANE(nj)
                          VIEW_UP_NEW(nj)=VIEW_UP(nj)
                        ENDDO
                        WINDOW_NEW(1)= WINDOW(1)
                        WINDOW_NEW(2)= WINDOW(2)
                        WINDOW_NEW(3)= WINDOW(3)
                        WINDOW_NEW(4)= WINDOW(4)
                        XNDC1=WINDOW_NEW(1)
                        XNDC2=WINDOW_NEW(2)
                        XNDC3=WINDOW_NEW(3)
                        XNDC4=WINDOW_NEW(4)
                        IZOOM(IW)=1
                        FIRST_ZOOM=.TRUE.
                      ENDIF
        
                    ELSE IF(CHOOSE(1:9).EQ.'Save view') THEN
                      NT_VIEW=NT_VIEW+1
                      CHAR3=CFROMI(NT_VIEW,'(I3)')
                      CALL TRIM(CHAR3,IBEG,IEND)
                      FILE_NAME='VIEW_'//CHAR3(IBEG:IEND)
                      CALL TRIM(FILE_NAME,IBEG,IEND)
                      CALL GKS_STRG(56,NCHAR,'Enter file name ['
     '                  //FILE_NAME(IBEG:IEND)//']',TEXT_STRING,
     '                  ERROR,*100)
                      IF(NCHAR.GT.0) THEN
                        FILE_NAME=TEXT_STRING
                      ENDIF
                      STRG='FEM define transform;w;'//FILE_NAME
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
        
                    ELSE IF(CHOOSE(1:11).EQ.'Select view') THEN
                      CALL DISPLAY_FILE(IW,NOCO,NTFILE,CO,'iptran',
     '                  FILE_NAME,STRING,ERROR,*100)
                      STRG='FEM define transform;r;'//FILE_NAME
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      STRG='FEM update mesh on 3'
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
        
                    ELSE IF(CHOOSE(1:9).EQ.'Transform') THEN
                      OPTION(1)='group'
                      OPTION(2)='object'
                      OPTION(3)='segment'
                      OPTION(4)='window'
                      OPTION(5)='Exit'
                      NTCH=5
                      CALL PRECHOICE1(1,IW,NOCH,NOCO,NTCH,COD,
     '                  'REQUEST',OPTION,STRING,ERROR,*100)
                      CHOOSE=OPTION(NOCH)
                      IF(CHOOSE(1:5).EQ.'group') THEN
                        DO nogrpl=1,NTGRPL
                          OPTION(nogrpl)=LAGRPL(nogrpl)
                        ENDDO
                        OPTION(NTGRPL+1)='Exit'
                        NTCH=NTGRPL+1
                        CALL PRECHOICE1(1,IW,NOCH,NOCO,NTCH,COD,
     '                    'REQUEST',OPTION,STRING,ERROR,*100)
                        CHOOSE=OPTION(NOCH)
                        IF(CHOOSE(1:4).ne.'Exit') THEN
                          OPTION11(1)='Rotate'
                          OPTION11(2)='Scale'
                          OPTION11(3)='Translate'
                          OPTION11(4)='Exit'
                          NTCH=4
                          CALL PRECHOICE1(1,IW,NOCH,NOCO,NTCH,COD,
     '                      'EVENT',OPTION11,STRING,ERROR,*100)
                          TRANSFORM_ACTIVE=.TRUE.
                          TRANSFORM_TYPE='GROUP'
                          NO_GROUP_TRANSFORM=NOCH
                          !Calc centre of all nodes of all polylines in group
                          ISEGM=0
                          NCENTRE=0
                          XCENTRE=0.0
                          YCENTRE=0.0
                          DO nosg=1,NTSG
      			    IF(CSEG(nosg)(1:8).EQ.'polyline') THEN
      			      INDEX_PLIN=IFROMC(CSEG(nosg)(53:57))
      			      IF(DOP) THEN
      				WRITE(OP_STRING,'('' polyline '
     '                            //'index = '',i2)') INDEX_PLIN
      				CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      			      ENDIF      
      			      IF(INLIST(INDEX_PLIN,
     '                          LIGRPL(1,NO_GROUP_TRANSFORM),
     '                          LIGRPL(0,NO_GROUP_TRANSFORM),N1LIST))
     '                          THEN
      				IF(DOP) THEN
      				  WRITE(OP_STRING,'('' ..is in '','
     '                              //'''group'')')
      				  CALL WRITES(IODI,OP_STRING,ERROR,
     '                              *9999)
      				  WRITE(OP_STRING,'('' Seg no='',I3,'
     '                              //''' CSEG(52:60) is '',A)') 
     '                              nosg,CSEG(nosg)(52:60)
      				  CALL WRITES(IODI,OP_STRING,ERROR,
     '                              *9999)
      				ENDIF
      				DO NOPTS=1,
     '                            NT_PLIN_SECTIONS(INDEX_PLIN)+1
      				  NCENTRE=NCENTRE+1
      				  XCENTRE=XCENTRE+
     '                              PLIN_DATA(1,NOPTS,INDEX_PLIN)
      				  YCENTRE=YCENTRE+
     '                              PLIN_DATA(2,NOPTS,INDEX_PLIN)
      				  WRITE(OP_STRING,'('' x='',E12.3,'
     '                              //''' y='',E12.3)') XCENTRE,YCENTRE
      				  CALL WRITES(IOOP,OP_STRING,ERROR,
     '                              *9999)
      				ENDDO
      				ISEGM=ISEGM+1
      				ISEGM_LIST(ISEGM)=nosg
      			      ENDIF
      			    ENDIF
                          ENDDO
                          ISEGM_LIST(0)=ISEGM
                          IF(NCENTRE.GT.0) THEN
                            XCENTRE=XCENTRE/REAL(NCENTRE)
                            YCENTRE=YCENTRE/REAL(NCENTRE)
                          ENDIF
                          IF(DOP) THEN
                            WRITE(OP_STRING,'('' isegm_list:'','
     '                        //'20I3)')
     '                        (ISEGM_LIST(ISEGM),ISEGM=1,ISEGM_LIST(0))
      			    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                          ENDIF
                        ENDIF
        
                      ELSE IF(CHOOSE(1:6).EQ.'object') THEN
                        TRANSFORM_TYPE='OBJECT'
        
                      ELSE IF(CHOOSE(1:7).EQ.'segment') THEN
                        TRANSFORM_TYPE='SEGMENT'
        
                      ELSE IF(CHOOSE(1:6).EQ.'window') THEN
                        OPTION11(1)='Pan'
                        OPTION11(2)='Rescale'
                        OPTION11(3)='Reset'
                        OPTION11(4)='Zoom (scale)'
                        OPTION11(5)='Zoom (box)'
                        OPTION11(6)='Zoom out'
                        OPTION11(7)='Exit'
                        NTCH=7
                        CALL PRECHOICE1(1,IW,NOCH,NOCO,NTCH,COD,
     '                    'EVENT',OPTION11,STRING,ERROR,*100)
                        TRANSFORM_ACTIVE=.TRUE.
                        TRANSFORM_TYPE='WINDOW'
                      ENDIF
        
                    ELSE IF(CHOOSE(1:7).EQ.'Write..') THEN
                      OPTION(1)='buffer'
                      OPTION(2)='iod file'
                      OPTION(3)='geometry'
                      OPTION(4)='problem'
                      OPTION(5)='Exit'
                      NTCH=5
                      CALL PRECHOICE1(1,IW,NOCH,NOCO,NTCH,COD,
     '                  'REQUEST',OPTION,STRING,ERROR,*100)
                      CHOOSE=OPTION(NOCH)
                      IF(CHOOSE(1:4).ne.'Exit') THEN
                        CALL TRIM(FILE00,IBEG,IEND)
                        CALL GKS_STRG(56,NCHAR,'Enter file name ['
     '                    //FILE00(IBEG:IEND)//']',TEXT_STRING,
     '                    ERROR,*100)
                        IF(NCHAR.GT.0) THEN
                          FILE_NAME=TEXT_STRING
                          FILE00=FILE_NAME
                        ELSE
                          FILE_NAME=FILE00
                        ENDIF
                      ENDIF
        
                      IF(CHOOSE(1:6).EQ.'buffer') THEN
                        CALL TRIM(FILE_NAME,IBEG,IEND)
                        STRG='write buffer'
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:8).EQ.'iod file') THEN
                        CALL TRIM(FILE_NAME,IBEG,IEND)
                        STRG='FEM write '//FILE_NAME(IBEG:IEND)//';iod'
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:8).EQ.'geometry') THEN
                        STRG='FEM define nodes;w;'//FILE_NAME
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        STRG='FEM define bases;w;'//FILE_NAME
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        STRG='FEM define elements;w;'//FILE_NAME
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:7).EQ.'problem') THEN
                        STRG='FEM define equation;w;'//FILE_NAME
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        STRG='FEM define materials;w;'//FILE_NAME
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        STRG='FEM define initial;w;'//FILE_NAME
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        STRG='FEM define solve;w;'//FILE_NAME
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
                      ELSE IF(CHOOSE(1:4).EQ.'Exit') THEN
                        CHOOSE=' '
                      ENDIF
        
                    ELSE IF(CHOOSE(1:10).EQ.'Zoom (box)') THEN
                      STRG='FEM zoom in on '//CFROMI(IW,'(I1)')
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      CHANGE=.FALSE.
        
                    ELSE IF(CHOOSE(1:8).EQ.'Zoom out') THEN
                      STRG='FEM zoom out on '//CFROMI(IW,'(I1)')
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      CHANGE=.FALSE.
        
                    ELSE IF(CHOOSE(1:3).EQ.'...') THEN
                      OPTION(1)='Display'
                      OPTION(2)='Fit'
                      OPTION(3)='Solve'
                      OPTION(4)='Update'
                      OPTION(5)='Exit'
                      NTCH=5
                      CALL PRECHOICE1(1,IW,NOCH,NOCO,NTCH,COD,
     '                  'REQUEST',OPTION,STRING,ERROR,*100)
                      CHOOSE=OPTION(NOCH)
                      IF(CHOOSE(1:7).EQ.'Display') THEN
                      ELSE IF(CHOOSE(1:3).EQ.'Fit') THEN
                        OPTION(1)='geometry'
                        OPTION(2)='fibre'
                        OPTION(3)='field'
                        OPTION(4)='Fourier'
                        OPTION(5)='Gauss'
                        OPTION(6)='Exit'
                        NTCH=6
                        CALL PRECHOICE1(1,IW,NOCH,NOCO,NTCH,COD,
     '                    'REQUEST',OPTION,STRING,ERROR,*100)
                        CHOOSE=OPTION(NOCH)
                        IF(CHOOSE(1:4).ne.'Exit') THEN
                          STRG='FEM fit '//CHOOSE
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          IF(CHOOSE(1:8).EQ.'geometry') THEN
                            STRG='FEM update node'
                            CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                            STRG='FEM update mesh'
                            CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          ENDIF
                          WRITE(OP_STRING,'('' ...completed'')')
      			  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                        ELSE IF(CHOOSE(1:4).EQ.'Exit') THEN
                          CHOOSE=' '
                        ENDIF
                      ELSE IF(CHOOSE(1:5).EQ.'Solve') THEN
                        STRG='FEM solve'
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      ELSE IF(CHOOSE(1:6).EQ.'Update') THEN
                        OPTION(1)='Gauss'
                        OPTION(2)='growth'
                        OPTION(3)='mesh'
                        OPTION(4)='nodes'
                        OPTION(5)='residuals'
                        OPTION(6)='xi'
                        OPTION(7)='Exit'
                        NTCH=7
                        CALL PRECHOICE1(1,IW,NOCH,NOCO,NTCH,COD,
     '                    'REQUEST',OPTION,STRING,ERROR,*100)
                        CHOOSE=OPTION(NOCH)
                        IF(CHOOSE(1:4).ne.'Exit') THEN
                          STRG='FEM update '//CHOOSE
                          CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                          WRITE(OP_STRING,'('' ...completed'')')
      			  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                        ELSE IF(CHOOSE(1:4).EQ.'Exit') THEN
                          CHOOSE=' '
                        ENDIF
                      ELSE IF(CHOOSE(1:4).EQ.'Exit') THEN
                        CHOOSE=' '
                      ENDIF
        
                    ELSE IF(CHOOSE(1:4).EQ.'Exit') THEN
                      IF(DISPLAY_VALUATOR_81) THEN
                        DISPLAY_VALUATOR_81=.FALSE.
                        CALL INPUT_MODE(81,1,'CHOICE','REQUEST',
     '                    ERROR,*100)
                      ENDIF
                      IF(DISPLAY_VALUATOR_82) THEN
                        DISPLAY_VALUATOR_82=.FALSE.
                        CALL INPUT_MODE(82,1,'CHOICE','REQUEST',
     '                    ERROR,*100)
                      ENDIF
                      IF(DISPLAY_VALUATOR_83) THEN
                        DISPLAY_VALUATOR_83=.FALSE.
                        CALL INPUT_MODE(83,1,'VALUATOR','REQUEST',
     '                    ERROR,*100)
                      ENDIF
                      IF(DISPLAY_VALUATOR_84) THEN
                        DISPLAY_VALUATOR_84=.FALSE.
                        CALL INPUT_MODE(84,1,'VALUATOR','REQUEST',
     '                    ERROR,*100)
                      ENDIF
                      IF(DISPLAY_VALUATOR_85) THEN
                        DISPLAY_VALUATOR_85=.FALSE.
                        CALL INPUT_MODE(85,1,'VALUATOR','REQUEST',
     '                    ERROR,*100)
                      ENDIF
                      IF(DISPLAY_VALUATOR_86) THEN
                        DISPLAY_VALUATOR_86=.FALSE.
                        CALL INPUT_MODE(86,1,'VALUATOR','REQUEST',
     '                    ERROR,*100)
                      ENDIF
                      IF(DISPLAY_VALUATOR_87) THEN
                        DISPLAY_VALUATOR_87=.FALSE.
                        CALL INPUT_MODE(87,1,'VALUATOR','REQUEST',
     '                    ERROR,*100)
                      ENDIF
                      IF(DISPLAY_VALUATOR_88) THEN
                        DISPLAY_VALUATOR_88=.FALSE.
                        CALL INPUT_MODE(88,1,'VALUATOR','REQUEST',
     '                    ERROR,*100)
                      ENDIF
                      IF(DISPLAY_VALUATOR_89) THEN
                        DISPLAY_VALUATOR_89=.FALSE.
                        CALL INPUT_MODE(89,1,'VALUATOR','REQUEST',
     '                    ERROR,*100)
                      ENDIF
                      IF(REFINE_ACTIVE) THEN
                        CALL INPUT_MODE(7,1,'CHOICE','REQUEST',
     '                    ERROR,*100)
                        REFINE_ACTIVE=.FALSE.
                      ELSE IF(TRANSFORM_ACTIVE) THEN
                        CALL INPUT_MODE(7,1,'CHOICE','REQUEST',
     '                    ERROR,*100)
                        TRANSFORM_ACTIVE=.FALSE.
                      ELSE
                        UPDATE=.FALSE.
                      ENDIF
                      CHANGE=.FALSE.
                      CHOOSE=' '
                    ENDIF
        
                  ELSE IF(ID_WS.EQ.95) THEN
                    IF(CHOOSE(1:12).EQ.'Read signals') THEN
                      CALL PRECHOICE2(IW,NOCO,COD,'SIGNAL',FILE_NAME,
     '                  PARAMETER_TYPE,Q,' ',STRING,ERROR,*100)
                      IF(FILE_NAME(1:4).ne.'Exit') THEN
                        CALL TRIM(FILE_NAME,IBEG,IEND)
                        FILE02=FILE_NAME
                        FILE00=FILE_NAME
                        STRG='FEM read '//FILE_NAME(IBEG:IEND)
     '                    //';signal'
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      ELSE IF(FILE_NAME(1:4).EQ.'Exit') THEN
                        CHOOSE=' '
                      ENDIF
                    ELSE IF(CHOOSE(1:5).EQ.'Setup') THEN
                      STRG='FEM define coordinates;r;h00_ep'
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      STRG='FEM define bases;r;h00_ep'
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      STRG='FEM define elements;r;h00_ep'
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      STRG='FEM define nodes;r;h00_ep'
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      STRG='FEM define field;r;h00_ep'
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      STRG='FEM define fit;r;h00_ep field'
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      STRG='FEM define data;r;h00_ep geometry'
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      STRG='FEM define data;c xi'
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      WRITE(OP_STRING,'('' ...completed'')')
      		      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    ELSE IF(CHOOSE(1:15).EQ.'Display signals') THEN
                      STRG='FEM display trace'
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      WRITE(OP_STRING,'('' ...completed'')')
      		      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    ELSE IF(CHOOSE(1:13).EQ.'Write signals') THEN
                      STRG='FEM define data;c signal'
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      CALL TRIM(FILE02,IBEG,IEND)
                      STRG='FEM write '//FILE02(IBEG:IEND)//';signal'
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      WRITE(OP_STRING,'('' ...completed'')')
      		      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    ELSE IF(CHOOSE(1:13).EQ.'Display field') THEN
                      STRG='FEM define map hammer label'
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      STRG='FEM cancel field;s on 4'
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      STRG='FEM cancel contour;s on 4'
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      STRG='FEM cancel data;s on 4'
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      STRG='FEM fit field'
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      STRG='FEM draw field on 4'
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      STRG='FEM draw contour on 4'
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      STRG='FEM hide contour on 4'
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      STRG='FEM draw data on 4'
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      STRG='FEM draw line on 4 rgb=grey'
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      STRG='FEM draw electrode on 4'
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      WRITE(OP_STRING,'('' ...completed'')')
      		      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    ELSE IF(CHOOSE(1:14).EQ.'Pick electrode') THEN
                      STRG='FEM pick electrode on 4'
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      STRG='FEM display trace'
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      WRITE(OP_STRING,'('' ...completed'')')
      		      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    ELSE IF(CHOOSE(1:8).EQ.'Cancel..') THEN
                      OPTION(1)='Map'
                      OPTION(2)='Signals'
                      OPTION(3)='Electrodes'
                      OPTION(4)='Exit'
                      NTCH=3
                      CALL PRECHOICE1(1,IW,NOCH,NOCO,NTCH,COD,
     '                  'REQUEST',OPTION,STRING,ERROR,*100)
                      CHOOSE=OPTION(NOCH)
                      IF(CHOOSE(1:3).EQ.'Map') THEN
                        STRG='FEM cancel map;s'
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      ELSE IF(CHOOSE(1:7).EQ.'Signals') THEN
                        STRG='FEM cancel trace;s'
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        WRITE(OP_STRING,'('' ...completed'')')
     			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      ELSE IF(CHOOSE(1:10).EQ.'Electrodes') THEN
                        STRG='FEM cancel electrodes;s'
                        CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                        WRITE(OP_STRING,'('' ...completed'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      ELSE IF(CHOOSE(1:4).EQ.'Exit') THEN
                        CHOOSE=' '
                      ENDIF
                    ELSE IF(CHOOSE(1:15).EQ.'Change colour..') THEN
                      STRG='FEM change colour;m on 4'
                      CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*100)
                      WRITE(OP_STRING,'('' ...completed'')')
      		      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    ELSE IF(CHOOSE(1:4).EQ.'Exit') THEN
                      UPDATE=.FALSE.
                      CHANGE=.FALSE.
                      CHOOSE=' '
                    ENDIF
                  ENDIF
                ENDIF
        
              ELSE IF(CLASS(1:8).EQ.'VALUATOR'
     '          .AND..NOT.FIRST_ZOOM.AND.ID_WS.GT.80) THEN
                VALUE=R4DATA(1)
                ID=ID_WS-80
                IF(TRANSFORM_ACTIVE.AND.
     '            TRANSFORM_TYPE(1:6).EQ.'WINDOW') THEN
                  !Apply transformation to window
                  CHANGE=.TRUE.
                  IF(IW.LE.2) THEN
                    IF(DOP) THEN
                      WRITE(OP_STRING,'('' ID='',I2,'' Value='','
     '                  //'E12.3)') ID,VALUE
      		      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      WRITE(OP_STRING,'('' IZOOM(iw)='',I1,'
     '                  //''' XNDC:'',4F6.2)') 
     '                  IZOOM(IW),(XNDC(IZOOM(IW),i,IW),i=1,4)
      		      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                    X0=(XNDC(IZOOM(IW),1,IW)+XNDC(IZOOM(IW),2,IW))/2.0 !x centre
                    Y0=(XNDC(IZOOM(IW),3,IW)+XNDC(IZOOM(IW),4,IW))/2.0 !y centre
                    DX= XNDC(IZOOM(IW),2,IW)-XNDC(IZOOM(IW),1,IW)      !x width
                    DY= XNDC(IZOOM(IW),4,IW)-XNDC(IZOOM(IW),3,IW)      !y width
                  ELSE IF(IW.EQ.3) THEN
                    X0=(WINDOW_NEW(1)+WINDOW_NEW(2))/2.0
                    Y0=(WINDOW_NEW(3)+WINDOW_NEW(4))/2.0
                    DX= WINDOW_NEW(2)-WINDOW_NEW(1)
                    DY= WINDOW_NEW(4)-WINDOW_NEW(3)
                    XNDC1=WINDOW_NEW(1)
                    XNDC2=WINDOW_NEW(2)
                    XNDC3=WINDOW_NEW(3)
                    XNDC4=WINDOW_NEW(4)
                  ENDIF
                  IF(CHOOSE(1:3).EQ.'Pan') THEN
                    IF(ID.EQ.1.OR.ID.EQ.3.OR.ID.EQ.5) THEN
                      XNDC1=XNDC(IZOOM(IW),1,IW)-VALUE*DX
                      XNDC2=XNDC(IZOOM(IW),2,IW)-VALUE*DX
                      XNDC3=XNDC(IZOOM(IW),3,IW)
                      XNDC4=XNDC(IZOOM(IW),4,IW)
                    ELSE IF(ID.EQ.2.OR.ID.EQ.4.OR.ID.EQ.6) THEN
                      XNDC1=XNDC(IZOOM(IW),1,IW)
                      XNDC2=XNDC(IZOOM(IW),2,IW)
                      XNDC3=XNDC(IZOOM(IW),3,IW)-VALUE*DY
                      XNDC4=XNDC(IZOOM(IW),4,IW)-VALUE*DY
                    ELSE IF(ID.EQ.7.OR.ID.EQ.8.OR.ID.EQ.9) THEN
                      SHIFT(ID-6)=VALUE*DIAG
                    ENDIF
                  ELSE IF(CHOOSE(1:7).EQ.'Rescale') THEN
                    IF(ID.EQ.1.OR.ID.EQ.3.OR.ID.EQ.5) THEN
                      XNDC1=X0-(1.0+VALUE)*DX/2.0
                      XNDC2=X0+(1.0+VALUE)*DX/2.0
                      XNDC3=XNDC(IZOOM(IW),3,IW)
                      XNDC4=XNDC(IZOOM(IW),4,IW)
                    ELSE IF(ID.EQ.2.OR.ID.EQ.4.OR.ID.EQ.6) THEN
                      XNDC1=XNDC(IZOOM(IW),1,IW)
                      XNDC2=XNDC(IZOOM(IW),2,IW)
                      XNDC3=Y0-(1.0+VALUE)*DY/2.0
                      XNDC4=Y0+(1.0+VALUE)*DY/2.0
                    ELSE IF(ID.EQ.7.OR.ID.EQ.8.OR.ID.EQ.9) THEN
                      IF(VALUE.GE.0.0) THEN
                        SCALE(ID-6)=VALUE
                      ELSE IF(VALUE.LT.-0.0001) THEN
                        SCALE(ID-6)=-1.0/VALUE
                      ENDIF
                    ENDIF
                  ELSE IF(CHOOSE(1:11).EQ.'Rotate data') THEN
                    IF(ID.EQ.7) THEN
                      FANGLE(1)=VALUE*PI
                    ELSE IF(ID.EQ.8) THEN
                      FANGLE(2)=VALUE*PI
                    ELSE IF(ID.EQ.9) THEN
                      FANGLE(3)=VALUE*PI
                    ENDIF
                  ELSE IF(CHOOSE(1:11).EQ.'Rotate view') THEN
                    IF(ID.EQ.7) THEN
                      ANGLE(1)=VALUE*PI
                    ELSE IF(ID.EQ.8) THEN
                      ANGLE(2)=VALUE*PI
                    ELSE IF(ID.EQ.9) THEN
                      ANGLE(3)=VALUE*PI
                    ENDIF
                  ELSE IF(CHOOSE(1:11).EQ.'Proj ref pt') THEN
                    IF(VALUE*DIAG*10.LE.FRONT_PLANE_DIST_NEW.AND.
     '                IW.EQ.3.AND.ID.EQ.9)THEN
                      WRITE(OP_STRING,*)' ---Error: Proj_ref_point ',
     '                  'will be behind the front clipping plane'
      		      CALL WRITES(IOER,OP_STRING,ERROR,*9999)
                      CHANGE=.FALSE.
                    ELSE
                      PROJ_REF_PT_NEW(ID-6)=VALUE*DIAG*10
                      IF(DOP) THEN
                        WRITE(OP_STRING,*) ' Proj_ref_pt(',ID-6,')=',
     '                    PROJ_REF_PT_NEW(ID-6)
      			CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      ENDIF
                    ENDIF
                  ELSE IF(CHOOSE(1:10).EQ.'View plane') THEN
                    VIEW_PLANE_DIST_NEW=VALUE*DIAG
                    IF(DOP) THEN
                      WRITE(OP_STRING,*) ' View plane dist=',
     '                  VIEW_PLANE_DIST_NEW
      		      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    ENDIF
                  ELSE IF(CHOOSE(1:10).EQ.'Back plane') THEN
                    IF(VALUE*DIAG.GE.FRONT_PLANE_DIST_NEW)THEN
                      WRITE(OP_STRING,*)' ---Error: Back ',
     '                  'plane will be in front of the front plane'
      		      CALL WRITES(IOER,OP_STRING,ERROR,*9999)
                      CHANGE=.FALSE.
                    ELSE
                      BACK_PLANE_DIST_NEW=VALUE*DIAG
                      WRITE(OP_STRING,*) ' Back plane dist=',
     '                  BACK_PLANE_DIST_NEW
      		      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    ENDIF
                  ELSE IF(CHOOSE(1:11).EQ.'Front plane') THEN
                    IF(VALUE*DIAG.LE.BACK_PLANE_DIST_NEW)THEN
                      WRITE(OP_STRING,*)' ---Error: Front ',
     '                  'plane will be behind the front plane'
      		      CALL WRITES(IOER,OP_STRING,ERROR,*9999)
                      CHANGE=.FALSE.
                    ELSE
                      FRONT_PLANE_DIST_NEW=VALUE*DIAG
                      WRITE(OP_STRING,*) ' Front plane dist=',
     '                  FRONT_PLANE_DIST_NEW
      		      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    ENDIF
                  ELSE IF(CHOOSE(1:11).EQ.'View ref pt') THEN
                    VIEW_REF_PT_NEW(ID-6)=VALUE*DIAG
                    IF(DOP) THEN
                      WRITE(OP_STRING,*) ' View ref pt(',ID-6,')=',
     '                  VIEW_REF_PT_NEW(ID-6)
      		      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                  ELSE IF(CHOOSE(1:10).EQ.'View plane') THEN
                    VIEW_PLANE_NEW(ID-6)=VALUE
                    IF(DOP) THEN
                      WRITE(OP_STRING,*) ' View plane(',ID-6,')=',
     '                  VIEW_PLANE_NEW(ID-6)
      		      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                  ELSE IF(CHOOSE(1:7).EQ.'View up') THEN
                    VIEW_UP_NEW(ID-6)=VALUE
                    IF(DOP) THEN
                      WRITE(OP_STRING,*) ' View_up(',ID-6,')=',
     '                  VIEW_UP_NEW(ID-6)
      		      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                  ELSE IF(CHOOSE(1:4).EQ.'Zoom') THEN
                    IF(IW.LT.3)THEN
                      XNDC1=X0-(1.0+VALUE)*DX/2.0
                      XNDC2=X0+(1.0+VALUE)*DX/2.0
                      XNDC3=Y0-(1.0+VALUE)*DY/2.0
                      XNDC4=Y0+(1.0+VALUE)*DY/2.0
                    ELSE IF(IW.EQ.3)THEN
                      XNDC1=X0-(1.0+VALUE)*DIAG
                      XNDC2=X0+(1.0+VALUE)*DIAG
                      XNDC3=Y0-(1.0+VALUE)*DIAG
                      XNDC4=Y0+(1.0+VALUE)*DIAG
                    ENDIF
                  ENDIF
        
                ELSE IF(TRANSFORM_ACTIVE.AND.TRANSFORM_TYPE(1:5).EQ.
     '            'GROUP') THEN
                  !Apply transformation to specified group
                  IF(DOP.AND.NO_GROUP_TRANSFORM.GT.0) THEN
                    WRITE(OP_STRING,'('' ..Transform Polyline '
     '                //'group'',I2,'' Label ='',A)') 
     '                NO_GROUP_TRANSFORM,LAGRPL(NO_GROUP_TRANSFORM)
      		    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    WRITE(OP_STRING,'('' ..Polyline indices:'','
     '                //'30I3)')
     '                (LIGRPL(N,NO_GROUP_TRANSFORM),
     '                N=1,LIGRPL(0,NO_GROUP_TRANSFORM))
      		    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                  IF(CHOOSE(1:6).EQ.'Rotate') THEN
                    CALL TRANSFORM_SEGMENT(ISEGM_LIST,IW,'rotate',
     '                XCENTRE,YCENTRE,VALUE,ERROR,*9999)
                  ELSE IF(CHOOSE(1:7).EQ.'Scale') THEN
                    CALL TRANSFORM_SEGMENT(ISEGM_LIST,IW,'scale',
     '                XCENTRE,YCENTRE,VALUE,ERROR,*9999)
                  ELSE IF(CHOOSE(1:9).EQ.'Translate') THEN
                    IF(ID.EQ.3) THEN
                      CALL TRANSFORM_SEGMENT(ISEGM_LIST,IW,
     '                  'x-translate',
     '                  XCENTRE,YCENTRE,VALUE,ERROR,*9999)
                    ELSE IF(ID.EQ.4) THEN
                      CALL TRANSFORM_SEGMENT(ISEGM_LIST,IW,
     '                  'y-translate',
     '                  XCENTRE,YCENTRE,VALUE,ERROR,*9999)
                    ENDIF
                  ENDIF
        
                ELSE IF(TRANSFORM_ACTIVE.AND.TRANSFORM_TYPE(1:6).EQ.
     '            'OBJECT') THEN
                  !Apply transformation to specified object
        
                ELSE IF(TRANSFORM_ACTIVE.AND.TRANSFORM_TYPE(1:7).EQ.
     '            'SEGMENT') THEN
                  !Apply transformation to specified segment
                ENDIF
        
              ELSE IF(CHOOSE(1:4).EQ.'Exit') THEN !exit when no event
                CHANGE=.FALSE.
                UPDATE=.FALSE.
              ENDIF
        
              IF(CHANGE)THEN
                IF(IW.LE.2) THEN
                  IF(XNDC1.GE.0.0.AND.XNDC2.LE.1.0.AND
     '              .XNDC3.GE.0.0.AND.XNDC4.LE.1.0) THEN
                    IF(UPDATE_ZOOM) THEN
                      IZOOM(IW)=IZOOM(IW)+1
                      UPDATE_ZOOM=.FALSE.
                    ENDIF
                    CALL ACWK(IW,1,ERROR,*100)
                    XNDC(IZOOM(IW),1,IW)=XNDC1
                    XNDC(IZOOM(IW),2,IW)=XNDC2
                    XNDC(IZOOM(IW),3,IW)=XNDC3
                    XNDC(IZOOM(IW),4,IW)=XNDC4
                    CALL WKST_WINDOW(IW,XNDC(IZOOM(IW),1,IW),
     '                XNDC(IZOOM(IW),2,IW),XNDC(IZOOM(IW),3,IW),
     '                XNDC(IZOOM(IW),4,IW),ERROR,*100)
                    CALL DAWK(IW,1,ERROR,*100)
                    IF(DOP) THEN
                      WRITE(OP_STRING,'('' IZOOM(iw)='',I1,'
     '                  //''' XNDC:'',4F6.2)')
     '                  IZOOM(IW),(XNDC(IZOOM(IW),i,IW),i=1,4)
      		      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                  ENDIF
                ELSE IF(IW.EQ.3) THEN
!old MPN 7-Apr-95: deleting all PHIGS calls
c                  WINDOW_NEW(1)= XNDC1
c                  WINDOW_NEW(2)= XNDC2
c                  WINDOW_NEW(3)= XNDC3
c                  WINDOW_NEW(4)= XNDC4
c                  CALL PHIGS_BUILD_XFORM_MATRIX3(FPT,FSHFT,
c     '              ANGLE,FSCALE,ISTATUS,VMATRIX,ERROR,*9999)
c                  CALL PHIGS_EVAL_VIEW_ORIEN_MATRIX3(VIEW_REF_PT_NEW,
c     '              VIEW_PLANE_NEW,VIEW_UP_NEW,ISTATUS,A_ORIENT,
c     '              ERROR,*9999)
c                  CALL PHIGS_COMPOSE_MATRIX3(VMATRIX,A_ORIENT,ISTATUS,
c     '              A_ORIENT_NEW,ERROR,*9999)
c                  IF(istatus.ne.0)THEN
c                    WRITE(OP_STRING,*)' ---Error from ',
c     '                'eval_view_orient_matrix=',ISTATUS
c      		    CALL WRITES(IOER,OP_STRING,ERROR,*9999)
c                  ELSE
c                    CALL PHIGS_EVAL_VIEW_MAP_MATRIX3(WINDOW_NEW,
c     '                VIEWPORT,MODE_PROJ,PROJ_REF_PT_NEW,
c     '                VIEW_PLANE_DIST_NEW,BACK_PLANE_DIST_NEW,
c     '                FRONT_PLANE_DIST_NEW,ISTATUS,A_MAP,ERROR,*9999)
c                    IF(istatus.ne.0)THEN
c                      WRITE(OP_STRING,*) ' ---Error from ',
c     '                  'eval_view_map_matrix=',ISTATUS
c      		      CALL WRITES(IOER,OP_STRING,ERROR,*9999)
c                    ELSE
c                      CALL PHIGS_SET_VIEW_REP3(IW,A_ORIENT_NEW,A_MAP,
c     '                  NPC_CLIP,ERROR,*9999)
c                    ENDIF
c                  ENDIF
c                  CALL PHIGS_BUILD_XFORM_MATRIX3(FIXED_PT,SHIFT,
c     '              FANGLE,SCALE,ISTATUS,A_TRANS,ERROR,*9999)
c                  IF(istatus.ne.0)THEN
c                    WRITE(OP_STRING,*)' ---Error from ',
c     '                'BUILD_XFORM_MATRIX3=',ISTATUS
c      		    CALL WRITES(IOER,OP_STRING,ERROR,*9999)
c                  ELSE
c                    CALL ACWK(3,1,ERROR,*100)
c                    CALL PHIGS_SET_EDIT_MODE(ERROR,*9999)
c                    CALL PHIGS_OPEN_STRUCT(ISVIEW,ERROR,*9999)
c                    CALL PHIGS_SET_ELEM_POINTER(ERROR,*9999)
c                    CALL PHIGS_SET_GLOBAL_XFORM3(A_TRANS,ERROR,*9999)
c                    CALL PHIGS_CLOSE_STRUCT(ERROR,*9999)
c                    CALL DAWK(3,1,ERROR,*100)
c                  ENDIF
                ENDIF
              ENDIF
              GO TO 101
        
C ***         Handle error condition
 100          CALL TRIM(ERROR,IBEG,IEND)
              WRITE(OP_STRING,'(A)') ' '//ERROR(IBEG:IEND)//'>DIALOG'
      	      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
 101          CONTINUE
            ENDDO
          ENDIF !gks & gks_ws_open
        
          IF(MACRO_KEY_EXECUTE) THEN !parse key macros
            MACRO_KEY_EXECUTE=.FALSE.
            DO nomacro=1,NT_MACRO(macro_key)
              STRG=MACRO_KEY_buffer(nomacro,macro_key)
              CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*150)
            ENDDO
          ELSE IF(MACRO_COMMAND_EXECUTE) THEN !parse command macros
            MACRO_COMMAND_EXECUTE=.FALSE.
            macro_command=1 !temporary
            DO nomacro=1,NT_MACRO_names(macro_command)
              STRG=MACRO_COMMAND_buffer(nomacro,macro_command)
              CALL PARSE(1,ISEG,CSEG,STRG,END,ERROR,*150)
            ENDDO
          ELSE IF(DO_EXAMPLE) THEN !read selected example file
            DO_EXAMPLE=.FALSE.
            CALL TRIM(EXAMPLE_NAME,IBEG,IEND)
            STRING=' > read '//EXAMPLE_NAME(IBEG:IEND)//';com;doc'
            WRITE(OP_STRING,'(A)') STRING(1:30)
      	    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE IF(SELECT_EXAMPLE) THEN !open selected example file
            SELECT_EXAMPLE=.FALSE.
            CALL TRIM(EXAMPLE_NAME,IBEG,IEND)
            STRING=' > open '//EXAMPLE_NAME(IBEG:IEND)//';com;doc'
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
        
          IF(END) THEN       !CTRL-Z has been used to quit
            CALL QUIT(END,ERROR,*150)
            CONTINUE=.FALSE.
            GOTO 200
          ENDIF
        
C ***     Parse string
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
        
          IF(END) THEN       !The QUIT command has been used to quit
            CONTINUE=.FALSE.
          ENDIF
          GOTO 200
        
C ***     Handle error condition
 150      CALL TRIM(ERROR,IBEG,IEND)
          WRITE(OP_STRING,'(A)') ' '//ERROR(IBEG:IEND)//'>DIALOG'
      	  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          IF(CTRLC) THEN
            CTRLC=.FALSE.
C CPB 23/10/93 This is done when CTRLC_AST is called so no need to do it here
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

      ELSE   !USE_SOCKET = .TRUE.

        CONNID1=0
        CONNID2=1

        IF (FSKLISTEN(CONNID1,PORT1) .EQ. -1) GOTO 9998
        IF (FSKLISTEN(CONNID2,PORT2) .EQ. -1) GOTO 9998
        
        CONTINUE=.TRUE.
        DO WHILE (CONTINUE)
C CPB 14/3/94 Adding timeout for sockets to enable graphics updates
          RETURNVAL=0
          DO WHILE(RETURNVAL.EQ.0)
            RETURNVAL=FSKSELECT(CONNID1,200) ! timeout after 200ms
            IF(RETURNVAL.EQ.-1) GOTO 9998
            CALL GXWAIT(0.0,ERR)              ! Update graphics
	  ENDDO				  
          IF (FSKREAD(LEN,SK_LONG_INT,1,CONNID1) .EQ. -1) GOTO 9998
          IF (FSKREAD(INTSTR,SK_CHAR,LEN+1,CONNID1) .EQ. -1) GOTO 9998
          CALL FSKC2F(STRING,LEN,INTSTR)

          CALL TRIM(STRING,IBEG,IEND)
      	  WRITE(OP_STRING,'(A)') STRING(IBEG:IEND)
      	  CALL WRITES(IOOP,OP_STRING,ERROR,*9998)
          IF(STRING(IBEG:IBEG+3).ne.'QUIT') THEN
      	    CALL PARSE(0,ISEG,CSEG,STRING,END,ERROR,*160)
      	    WRITE(OP_STRING,'(/1X,A)') '> '
      	    CALL WRITES(IOOP,OP_STRING,ERROR,*9998)
          ENDIF

          GOTO 210

C ***   Handle error condition
 160      CALL TRIM(ERROR,IBEG,IEND)
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
      CALL TRIM(ERROR,IBEG,IEND)
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
      

      SUBROUTINE FIND_FILE(NOFILE_START,NTFILE,FILE_EXT,FILE_LIST,
     '  ERROR,*)

C#### Subroutine: FIND_FILE
C###  Description:
C###    FIND_FILE finds files in current directory.

      IMPLICIT NONE
      INCLUDE '($RMSDEF)'
!     Parameter List
      INTEGER NOFILE_START,NTFILE
      CHARACTER ERROR*(*),FILE_EXT*(*),FILE_LIST(*)*(*)
!     Local Variables
      INTEGER CLOCAT,CONTEXT,IBEG,IBRA,IDOT,IEND,INDEX,LIB$FIND_FILE,
     '  LIB$FIND_FILE_END,NOFILE
      CHARACTER RESULT*100

      CALL ENTERS('FIND_FILE',*9999)
      CONTEXT=0
      NOFILE=0
      DO WHILE(LIB$FIND_FILE('*.'//FILE_EXT,RESULT,CONTEXT,,,,)
     '  .EQ.RMS$_NORMAL.AND.NOFILE.LE.NOFILE_START+19)
        NOFILE=NOFILE+1
        IF(NOFILE.GT.NOFILE_START) THEN
          CALL TRIM(RESULT,IBEG,IEND)
          IBRA=CLOCAT(']',RESULT(IBEG:IEND))
          IDOT=CLOCAT('.',RESULT(IBRA+1:IEND))
          FILE_LIST(NOFILE-NOFILE_START)=RESULT(IBRA+1:IBRA+IDOT-1)
        ENDIF
      ENDDO
      INDEX=LIB$FIND_FILE_END(CONTEXT)
      NTFILE=NOFILE-NOFILE_START

      CALL EXITS('FIND_FILE')
      RETURN
 9999 CALL ERRORS('FIND_FILE',ERROR)
      CALL EXITS('FIND_FILE')
      RETURN 1
      END


      SUBROUTINE GET_CMISS_EXAMPLES(CMISS_EXAMPLES,*)

C#### Subroutine: GET_CMISS_EXAMPLES
C###  Description:
C###    GET_CMISS_EXAMPLES returns CMISS examples directory as
C###    string in CMISS_EXAMPLES.

      IMPLICIT NONE
!     Parameter List
      CHARACTER CMISS_EXAMPLES*255
!     Local Variables

C GMH 11/12/96 Trace, etc uninitialised for this call
C      CALL ENTERS('GET_CMISS_EXAMPLES',*9999)

      CMISS_EXAMPLES='cmiss$examples:'

C GMH 11/12/96 Trace, etc uninitialised for this call
C      CALL EXITS('GET_CMISS_EXAMPLES')
C      RETURN
C 9999 CALL ERRORS('GET_CMISS_EXAMPLES',ERROR)
C      CALL EXITS('GET_CMISS_EXAMPLES')
C      RETURN 1
      RETURN
      END


      SUBROUTINE OPENF - archive old file handling

C CPB 20/4/94 Changing over all external cmiss formatted files to
C Sequential stream_lf files. Direct access files are handled by
C copying the external file to a dummy scratch interal direct access
C file.
C old start
C      IF(DEVICE.EQ.'TERM') THEN
C        OPEN(UNIT=IUNIT,FILE=FILE,STATUS=STATUS,RECL=IRECL,
C     '    IOSTAT=IOSTAT)
C        IF(iostat.ne.0) THEN
C          ERROR=CFROMI(IOSTAT,'(I3)')
C          ERROR=' Iostat='//ERROR(1:3)//
C     '      ' error occurred in OPENF(UNIT='//UNIT//')'
C          GOTO 9999
C        ENDIF
C
C      ELSE IF(DEVICE.EQ.'DISK') THEN
C        IF(ACCESS.EQ.'SEQUEN') THEN
C!news - GBS 12-08-92  -  write new files as stream_lf format
C! CPB 3/9/92 added carriage control for new files
C      	  IF(STATUS.EQ.'NEW') THEN
C      	    OPEN(UNIT=IUNIT,FILE=FILE,STATUS=STATUS,
C     '	      ACCESS='SEQUENTIAL',FORM=FORM,RECL=IRECL,IOSTAT=IOSTAT,
C     '	      CARRIAGECONTROL='NONE',RECORDTYPE='STREAM_LF')
C          ELSE
C      	    OPEN(UNIT=IUNIT,FILE=FILE,STATUS=STATUS,
C     '	      ACCESS='SEQUENTIAL',FORM=FORM,RECL=IRECL,IOSTAT=IOSTAT)
C          ENDIF
C!newe
C          IF(iostat.ne.0) THEN
C            ERROR=' Iostat='//CFROMI(IOSTAT,'(I3)')
C            IF(IOSTAT.EQ.29) THEN
C              CALL TRIM(FILE,IBEG,IEND)
C              ERROR=ERROR(1:11)
C     '          //' (file '//FILE(IBEG:IEND)//' not found)'
C     '          //' in OPENF(UNIT='//UNIT//')'
C            ELSE
C              ERROR=ERROR(1:11)
C     '          //' in OPENF(UNIT='//UNIT//')'
C            ENDIF
C            GOTO 9999
C          ENDIF
C
C        ELSE IF(ACCESS.EQ.'APPEND') THEN
C          OPEN(UNIT=IUNIT,FILE=FILE,STATUS=STATUS,
C     '         ACCESS='APPEND',FORM=FORM,RECL=IRECL,IOSTAT=IOSTAT)
C          IF(iostat.ne.0) THEN
C            ERROR=' Iostat='//CFROMI(IOSTAT,'(I3)')
C            IF(IOSTAT.EQ.29) THEN
C              CALL TRIM(FILE,IBEG,IEND)
C              ERROR=ERROR(1:11)
C     '          //' (file '//FILE(IBEG:IEND)//' not found)'
C     '          //' in OPEN(UNIT='//UNIT//')'
C            ELSE
C              ERROR=ERROR(1:11)
C     '          //' in OPEN(UNIT='//UNIT//')'
C            ENDIF
C            GOTO 9999
C          ENDIF
C
C        ELSE IF(ACCESS.EQ.'DIRECT') THEN
C          IF(STATUS.EQ.'OLD') THEN !old, direct access, disk file.
C            !Open input file as sequential access in order to copy it to
C            !a direct access file. This copes with the problem of files
C            !being of fixed or variable record length (files are left as
C            !sequential, fixed length records when opened as direct access
C            !but EVE leaves files as sequential variable record length)
C            OPEN(UNIT=99,FILE=FILE,STATUS='OLD',ACCESS='SEQUENTIAL',
C     '        FORM=FORM,IOSTAT=IOSTAT,RECORDTYPE='VARIABLE',READONLY)
C            IF(IOSTAT.EQ.0) THEN !Sequential,variable record-length file
C              OPEN(UNIT=IUNIT,FILE=FILE,STATUS='NEW',ACCESS='DIRECT',
C     '          FORM=FORM,RECL=IRECL,IOSTAT=IOSTAT,RECORDTYPE='FIXED')
C	      IREC=0 !new
C 10           READ(99,'(A)',IOSTAT=IOSTAT) DUMMY
C              IF(IOSTAT.EQ.0) THEN
C                INQUIRE(UNIT=IUNIT,NEXTREC=IREC)
C                WRITE(IUNIT,'(A)',REC=IREC,IOSTAT=IOSTAT) DUMMY
C                GO TO 10
C              ELSE !finished reading file
C                CALL CLOSEF(99,ERROR,*9999)
C                CALL CLOSEF(IUNIT,ERROR,*9999)
C                OPEN(UNIT=IUNIT,FILE=FILE,STATUS='OLD',ACCESS='DIRECT',
C     '            FORM=FORM,RECL=IRECL,IOSTAT=IOSTAT,RECORDTYPE='FIXED',
C     '            READONLY)
C                IF(DOP) THEN
C                  WRITE(OP_STRING,'('' Sequential,variable'',
C     '              '' record-length file'',
C     '              '' copied to readonly,direct,fixed'',
C     '              '' record-length file'')')
C      		  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                ENDIF
C              ENDIF
C            ELSE IF(IOSTAT.EQ.44) THEN !Sequential,fixed record-length file
C              OPEN(UNIT=IUNIT,FILE=FILE,STATUS='OLD',ACCESS='DIRECT',
C     '          FORM=FORM,RECL=IRECL,IOSTAT=IOSTAT,RECORDTYPE='FIXED',
C     '          READONLY)
C              IF(iostat.ne.44)THEN
C                IF(DOP) THEN
C                  WRITE(OP_STRING,'('' Sequential,fixed '
C     '              //'record-length file re-opened as readonly,'
C     '              //'direct,fixed, record-length file'')')
C      		  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                ENDIF
C              ELSEIF(IOSTAT.EQ.44)THEN !For ftp files 11-2-93 AJP,CPB
C                WRITE(OP_STRING,'('' File needs to be converted '',
C     '              ''to a fixed-length record format file to '',
C     '              ''enable direct access of records'')')
C      		  CALL WRITES(IOER,OP_STRING,ERROR,*9999)
C                  GOTO 9999
C              ENDIF
C            ENDIF
C
C          ELSE IF(STATUS.EQ.'NEW') THEN
C            OPEN(UNIT=IUNIT,FILE=FILE,STATUS='NEW',ACCESS='DIRECT',
C     '        FORM=FORM,RECL=IRECL,IOSTAT=IOSTAT,RECORDTYPE='FIXED')
C          ENDIF
C
C          IF(iostat.ne.0) THEN
C            ERROR=CFROMI(IOSTAT,'(I3)')
C            ERROR=' Iostat='//ERROR(1:3)//' error in OPENF(UNIT='
C     '        //UNIT//')'
C            IF(IOSTAT.EQ.29) THEN
C              CALL TRIM(ERROR,IBEG,IEND)
C              ERROR=' File not found: '//ERROR(IBEG:IEND)
C            ELSE IF(IOSTAT.EQ.43) THEN
C              CALL TRIM(ERROR,IBEG,IEND)
C              ERROR=' Filename invalid: '//ERROR(IBEG:IEND)
C            ELSE IF(IOSTAT.EQ.44) THEN
C              CALL TRIM(ERROR,IBEG,IEND)
C              ERROR=' Inconsistent record type: '//ERROR(IBEG:IEND)
C            ENDIF
C            WRITE(OP_STRING,'('' Filename: '',A)') FILE
C      	    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C            GOTO 9999
C          ENDIF
C        ELSE
C          ERROR=' ACCESS='//ACCESS//' is invalid'
C          GOTO 9999
C        ENDIF
C
C      ELSE
C        ERROR=' DEVICE='//DEVICE//' is invalid'
C        GOTO 9999
C      ENDIF
C old end



      SUBROUTINE POST_FILE(FILE_NAME)

C#### Subroutine: POST_FILE
C###  Description:
C###    POST_FILE posts file to printer.

      IMPLICIT NONE
!     Parameter List
      CHARACTER FILE_NAME*(*)
!     Local Variables
      INTEGER IBEG,IEND,ISTATUS,LIB$SPAWN

      CALL TRIM(FILE_NAME,IBEG,IEND)
      ISTATUS=LIB$SPAWN('post '//FILE_NAME(IBEG:IEND)//'.ps',
     '  ,,,,,,,,,,)

      RETURN
      END


      SUBROUTINE PURGE_FILE(FILE_NAME)

C#### Subroutine: PURGE_FILE
C###  Description:
C###    PURGE_FILE purges file versions.

      IMPLICIT NONE
!     Parameter List
      CHARACTER FILE_NAME*(*)
!     Local Variables
      INTEGER IBEG,IEND,ISTATUS,LIB$SPAWN

      CALL TRIM(FILE_NAME,IBEG,IEND)
      ISTATUS=LIB$SPAWN('purge '//FILE_NAME(IBEG:IEND)//'.ps',
     '  ,,,,,,,,,,)

      RETURN
      END
    
 

Module FE11
===========

C KAT 2002-07-16: 2D Triangle stuff archived from IPREFI

        IF(IBT(1,1,nb).EQ.3) THEN !Simplex elements
          TRIANGLE=.TRUE. !use Triangle lib.
C         Flags follow the sequence prqa_aAcevngBPNEIOXzoYYSiFlsCQVVVhf
          DO i=1,40
            FLAGS(i)=0
          ENDDO
          FLAGS(2)=1 !r
          FLAGS(30)=1 !Q
          FLAGS(3)=1 !q
          IF(INPUT_TYPE.EQ.1) FLAGS(1)=1 !p
          IF(NNT(nb).EQ.6) FLAGS(21)=1 !o2
          
          FORMAT='(/'' Enter method for defining refinement [1]:'''//
     '    '/''  (1) By defining maximum area for specified elements'''//
     '      '/''  (2) By field description of maximum element area'''//
     '      '/''  (3) By automatic adaption'''//
     '      '/$,''    '',I1)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=1
          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) REFINE_TYPE=IDATA(1)
          
          IF(REFINE_TYPE.EQ.1) THEN
            FLAGS(6)=1 !a
            
            noelem=0
 6720       FORMAT='($,'' Enter element #s/name [EXIT]: '',I5)'
            IF(IOTYPE.EQ.3) THEN
              noelem=noelem+1
              IF(noelem.LE.NEELEM(0,nr)) THEN
                ne=NEELEM(noelem,nr)
                IDATA(1)=ne
                NELIST(0)=1
                NELIST(1)=ne
              ELSE
                IDATA(0)=0
              ENDIF
            ENDIF
 6760       CDATA(1)='ELEMENTS' !for use with group input
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,
     '        0,NET(nr),LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '        ERROR,*9999)
            IF(IDATA(1).NE.0) THEN !not default exit
              NELIST(0)=IDATA(0)
              DO n=1,IDATA(0)
                NELIST(n)=IDATA(n)
                ne=IDATA(n)
                IF(.NOT.INLIST(ne,NEELEM(1,nr),
     '            NEELEM(0,nr),N1)) THEN
                  WRITE(OP_STRING,'('' Element '',I5,'' is not '
     '              //'in the current region'')') ne
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  GOTO 6760
                ENDIF
              ENDDO !n
              
C             Get area for first element in group
              ne=NELIST(1) !rest of group filled at end of loop
              AREA=0.0d0
              CALL XPXE(NBJ,NKJE(1,1,1,ne),NPF(1,1),
     '          NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '          SE(1,1,ne),XA,XE,XP,ERROR,*9999)
              DO ng=1,NGT(nb)
                CALL XEXG(NBJ(1,ne),ng,nr,PG,XE,XG,ERROR,*9999)
                CALL XGMG(NJ_LOC(njl_fibr,0,nr),NIT(nb),nb,nr,DXIX,GL,
     '            GU,RG(ng),XG,ERROR,*9999)
                RWG=RG(ng)*WG(ng,nb)
                AREA=AREA+RWG
              ENDDO
              RDEFLT(1)=AREA
              WRITE(CHAR2,'(F12.4)') RDEFLT(1)
              CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
              FORMAT='($,'' Element area criteria is '
     '          //'['//CHAR2(IBEG2:IEND2)//']: '',F12.4)'
              IF(IOTYPE.EQ.3) RDATA(1)=RDEFLT(1)
              CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '          IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '          INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) AREA=RDATA(1)
C             Apply area to all of elements group
              DO n=1,NELIST(0)
                ne=NELIST(n)
                PLSG_ELEMENT_AREA(ne)=AREA 
              ENDDO !n
              
              GO TO 6720 !for more elements
            ENDIF !idata(1).ne.0
            
          ELSE IF(REFINE_TYPE.EQ.2) THEN
            FLAGS(34)=1 !field
            FORMAT='('' Enter the field number [1]:'''//
     '        '/$,''    '',I1)'
            IDEFLT(1)=1
            IF(IOTYPE.EQ.3) IDATA(1)=IDEFLT(1)
            CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,6,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) AREA_FIELD=IDATA(1)
            AREA_FIELD=NJ_LOC(NJL_FIEL,AREA_FIELD,nr)
          ELSE
            FLAGS(6)=1 !a
            
            FORMAT='('' Enter the desired percentage error [5]:'''//
     '        '/$,''    '',I1)'
            RDEFLT(1)=5.0d0
            IF(IOTYPE.EQ.3) RDATA(1)=RDEFLT(1)
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) MAX_PERCENTAGE_ERR=RDATA(1)/100.0d0
            
            FORMAT='('' Enter singularity strength parameter '
     '        //'[NO_SINGULARITY]:''/$,''    '',I1)'
            RDEFLT(1)=NBSC(2,nb)
            IF(IOTYPE.EQ.3) RDATA(1)=RDEFLT(1)
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) SINGULARITY_FACTOR=RDATA(1)
            
            MAX_ELEMENT_ERR_NORM=MAX_PERCENTAGE_ERR*
     '        DSQRT((ENERGY_NORM_TOTAL**2+ERR_NORM_TOTAL**2)/
     '        NEELEM(0,nr))
            
            IF(DOP) THEN 
              WRITE(OP_STRING,'(/''  MAX_ELEMENT_ERR_NORM = '',D12.5)') 
     '          MAX_ELEMENT_ERR_NORM 
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
            
            DO no_nelist=1,NELIST(0)
              ne=NELIST(no_nelist)
              
              REFINEMENT_FACTOR=NEERR(ne,2)/MAX_ELEMENT_ERR_NORM
              IF(DOP) THEN
                WRITE(OP_STRING,'(/''  Refinement factor = '',D12.5)') 
     '            REFINEMENT_FACTOR
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
              
              AREA=0.0d0
              CALL XPXE(NBJ,NKJE(1,1,1,ne),NPF(1,1),
     '          NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '          SE(1,1,ne),XA,XE,XP,ERROR,*9999)
              DO ng=1,NGT(nb)
                CALL XEXG(NBJ(1,ne),ng,nr,PG,XE,XG,ERROR,*9999)
                CALL XGMG(NJ_LOC(njl_fibr,0,nr),NIT(nb),nb,nr,DXIX,GL,
     '          GU,RG(ng),XG,ERROR,*9999)
                RWG=RG(ng)*WG(ng,nb)
                AREA=AREA+RWG
              ENDDO
              
              PLSG_ELEMENT_AREA(ne)=(REFINEMENT_FACTOR**(-1.0d0/
     '          SINGULARITY_FACTOR))*AREA 
              
              IF (DOP) THEN
                WRITE(OP_STRING,'(/''  Current area = '',D12.5)') AREA
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(/''  Required max. area = '',D12.5)') 
     '            PLSG_ELEMENT_AREA(ne)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO
          ENDIF
          
          RDEFLT(1)=RMAX
          FORMAT='($,'' Enter the maximum triangle area  '
     '      //'[NO_MAXIMUM]: '',F12.4)'
          IF(IOTYPE.EQ.3) RDATA(1)=RDEFLT(1)
          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) MAX_AREA=RDATA(1)
          IF(MAX_AREA.LT.RMAX) THEN     
            FLAGS(4)=1 !a
            FLAGS(5)=1 !fixed
          ENDIF
          
          RDEFLT(1)=MIN_ANGLE
          WRITE(CHAR2,'(F12.4)') RDEFLT(1)
          CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
          FORMAT='($,'' Enter the minimum angle size in '
     '      //'triangulation ['//CHAR2(IBEG2:IEND2)//']: '',F12.4)'
          IF(IOTYPE.EQ.3) RDATA(1)=RDEFLT(1)
          CALL GINOUT(IOTYPE,5,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            MIN_ANGLE=RDATA(1)
            FLAGS(3)=1 !q
          ENDIF
          
          IF(INPUT_TYPE.EQ.1) THEN
            FORMAT='($,'' Allow boundary segment splitting [Y]? '',A)'
            IF(IOTYPE.EQ.3) ADATA(1)='Y'
            CALL GINOUT(IOTYPE,1,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,AYES,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              IF(ADATA(1).EQ.'N') FLAGS(22)=1 !Y
            ENDIF
            
            FORMAT='($,'' Allow internal segment splitting [Y]? '',A)'
            IF(IOTYPE.EQ.3) ADATA(1)='Y'
            CALL GINOUT(IOTYPE,1,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,AYES,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              IF(ADATA(1).EQ.'N') FLAGS(22)=1 !Y
            ENDIF
          ENDIF
          
          FORMAT='($,'' Enter the maximum number of Steiner '
     '      //'points [INF]? '',I10)'
          IDEFLT(1)=IMAX
          IF(IOTYPE.EQ.3) IDATA(1)=IDEFLT(1)
          CALL GINOUT(IOTYPE,3,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IF(IDATA(1).LT.IMAX) THEN
              NUM_STEINER_POINTS=IDATA(1)
              FLAGS(23)=1 !S
            ENDIF
          ENDIF
          
          FORMAT='($,'' Set up Voronoi information [N]? '',A)'
          IF(IOTYPE.EQ.3) ADATA(1)='N'
          CALL GINOUT(IOTYPE,1,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IF(ADATA(1).EQ.'Y') THEN
              FLAGS(10)=1 !v
              FLAGS(11)=1 !n
            ENDIF
          ENDIF
          
          FORMAT='('' Enter Delauney Triangulation method [1]:'''//
     '      '/''   (1) Divide-and-Conquer'''//
     '      '/''   (2) Divide-and-Conquer with vertical cuts only'''//
     '      '/''   (3) Incremental'''//
     '      '/''   (4) Sweepline'''//
     '      '/$,''    '',I1)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=1
          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) METHOD=IDATA(1)
          IF(METHOD.EQ.2) THEN
            FLAGS(27)=1 !l
          ENDIF
          IF(METHOD.EQ.3) THEN
            FLAGS(25)=1 !i
          ENDIF
          IF(METHOD.EQ.4) THEN
            FLAGS(26)=1 !F
          ENDIF

        ELSE


C MPN 17July2000: archived from IPACTI for HMT

      IF(KTYP59(nr).EQ.1.OR.      !SS tension-length-Ca relation
     '   KTYP59(nr).EQ.2) THEN    !Steady State HMT
        IF(KTYP59(nr).EQ.1) THEN         !SS tension-length-Ca relation
          RDEFLT(1)=1.0d2
          FORMAT='($,'' Enter max isometric tension at ext.ratio=1 '
     '      //'(Tref) [100kPa]: '',D11.4)'
        ELSE IF(KTYP59(nr).EQ.2) THEN    !Steady State HMT
          RDEFLT(1)=1.25d2
          FORMAT='($,'' Enter max isometric tension at ext.ratio=1 '
     '      //'(Tref) [125kPa]: '',D11.4)'
        ENDIF
        IF(IOTYPE.EQ.3) RDATA(1)=Tref
        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) Tref=RDATA(1)

        RDEFLT(1)=1.45d0
        IF(KTYP59(nr).EQ.1) THEN         !SS tension-length-Ca relation
          FORMAT='($,'' Enter non-dimensional slope parameter (beta)'
     '      //' [1.45]: '',D11.4)'
        ELSE IF(KTYP59(nr).EQ.2) THEN    !Steady State HMT
          FORMAT='($,'' Enter non-dim. slope parameter for T0 (beta0)'
     '      //' [1.45]: '',D11.4)'
        ENDIF
        IF(IOTYPE.EQ.3) RDATA(1)=T0_beta
        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) T0_beta=RDATA(1)

        IF(KTYP59(nr).EQ.1) THEN         !SS tension-length-Ca relation
          RDEFLT(1)=0.5d0
          FORMAT='($,'' Enter c50 for [Ca]i saturation curve (0<c<1)'
     '      //' [0.5]: '',D11.4)'
          IF(IOTYPE.EQ.3) RDATA(1)=Ca_c50
          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) Ca_c50=RDATA(1)

          RDEFLT(1)=3.0d0
          FORMAT='($,'' Enter Hill coeff. for [Ca]i saturation curve'
     '      //' (h) [3.0]: '',D11.4)'
          IF(IOTYPE.EQ.3) RDATA(1)=Ca_h
          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) Ca_h=RDATA(1)

          RDEFLT(1)=1.0d0
          FORMAT='($,'' Enter max [Ca]i (Ca_max) [1]: '',D11.4)'
          IF(IOTYPE.EQ.3) RDATA(1)=Ca_max
          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) Ca_max=RDATA(1)

          RDEFLT(1)=0.0d0
          FORMAT='($,'' Enter initial calcium level [Ca]i'
     '      //' [0]: '',D11.4)'
          IF(IOTYPE.EQ.3) RDATA(1)=FEXT(4,1,1) !1st Gauss pt in 1st elem
          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              DO ng=1,NGT(NBH(NH_LOC(1,nx),1,ne))
                FEXT(4,ng,ne)=RDATA(1)
              ENDDO !ng
            ENDDO !noelem (ne)
          ENDIF

        ELSE IF(KTYP59(nr).EQ.2) THEN    !Steady State HMT
          RDEFLT(1)=4.25d0
          FORMAT='($,'' Enter ref. Hill coeff. for [Ca]i saturation'
     '      //' curve (n_ref) [4.25]: '',D11.4)'
          IF(IOTYPE.EQ.3) RDATA(1)=HMT_n_ref
          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) HMT_n_ref=RDATA(1)

          RDEFLT(1)=1.95d0
          FORMAT='($,'' Enter non-dim. slope parameter for'
     '      //' Hill coeff. (beta1) [1.95]: '',D11.4)'
          IF(IOTYPE.EQ.3) RDATA(1)=HMT_n_beta
          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) HMT_n_beta=RDATA(1)

          RDEFLT(1)=5.33d0
          FORMAT='($,'' Enter ref. pC50 for [Ca]i saturation'
     '      //' curve (pC50_ref) [5.33]: '',D11.4)'
          IF(IOTYPE.EQ.3) RDATA(1)=HMT_pC50_ref
          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) HMT_pC50_ref=RDATA(1)

          RDEFLT(1)=0.31d0
          FORMAT='($,'' Enter non-dim. slope parameter for'
     '      //' pC50 (beta2) [0.31]: '',D11.4)'
          IF(IOTYPE.EQ.3) RDATA(1)=HMT_pC50_beta
          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) HMT_pC50_beta=RDATA(1)

          RDEFLT(1)=0.0d0
          FORMAT='($,'' Enter initial calcium level [Ca]i'
     '      //' [0.0mM]: '',D11.4)'
          IF(IOTYPE.EQ.3) RDATA(1)=FEXT(4,1,1) !1st Gauss pt in 1st elem
          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              DO ng=1,NGT(NBH(NH_LOC(1,nx),1,ne))
                FEXT(4,ng,ne)=RDATA(1)
              ENDDO !ng
            ENDDO !noelem (ne)
          ENDIF

        ENDIF

C MPN 17July2000: archived from IPACTI for fading memory model

        RDEFLT(1)=1.0D-2
        FORMAT='(/$,'' Enter the time step [0.01 s]: '',D11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=DEL_T
        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) DEL_T=RDATA(1)

        RDEFLT(1)=2.0D0
        FORMAT='(/$,'' Enter the slope of the stress/velocity relation'
     '    //' in stretching [2.0 kPa]: '',D11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=TV_SLO
        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) TV_SLO=RDATA(1)
                           
        RDEFLT(1)=1.8D0
        FORMAT='($,'' Enter the ratio of the yield tension to the'//
     '         ' isometric tension [1.8]: '',D11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=YIELDR
        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) YIELDR=RDATA(1)

        RDEFLT(1)=2.0D0
        FORMAT='($,'' Enter the static non-linearity parameter, "a" '//
     '         '[2.0]: '',D11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=SNLPA
        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) SNLPA=RDATA(1)

        IDEFLT(1)=3
        FORMAT='(/$,'' Enter the number of dynamic terms in the'//
     '    ' material response function [3]: '',I2)'
        IF(IOTYPE.EQ.3) IDATA(1)=NTACTV
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NTACTV=IDATA(1)

        DO nactv=1,NTACTV
          IF(nactv.EQ.1) THEN
            CHAR1='1st'
            RDEFLT(1)=17.0D0
            RDEFLT(2)=32.0D0
          ELSE IF(nactv.EQ.2) THEN
            CHAR1='2nd'
            RDEFLT(1)=26.0D0
            RDEFLT(2)=1000.0D0
          ELSE IF(nactv.EQ.3) THEN
            CHAR1='3rd'
            RDEFLT(1)=42.0D0
            RDEFLT(2)=5000.0D0
          ENDIF
          CALL TRIM(CHAR1,IBEG1,IEND1)
          WRITE(CHAR2,'(F5.0)') RDEFLT(1)
          WRITE(CHAR3,'(F6.0)') RDEFLT(2)
          CALL TRIM(CHAR2,IBEG2,IEND2)
          CALL TRIM(CHAR3,IBEG3,IEND3)
          CHAR2=CHAR2(IBEG2:IEND2)//','//CHAR3(IBEG3:IEND3)//' per sec'
          CALL TRIM(CHAR2,IBEG2,IEND2)
          FORMAT='($,'' Enter the coefficient and rate constant '
     '      //'for the '//CHAR1(IBEG1:IEND1)//' term ['
     '      //CHAR2(IBEG2:IEND2)//']: '',2D11.4)'
          IF(IOTYPE.EQ.3) THEN
            RDATA(1)=ACOEFF(nactv)
            RDATA(2)=ALFA(nactv)
          ENDIF
          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,2,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            ACOEFF(nactv)=RDATA(1)
            ALFA(nactv)=RDATA(2)
          ENDIF
        ENDDO



C 18 August 1997
      SUBROUTINE IPGRID(IBT,IDO,INP,NAN,NBH,NBJ,NEELEM,NGAP,NHE,NJE,
     '  NKE,NPF,NPNE,NQE,NQGE,nr,NRE,NVHE,NVJE,NVQ,NW,NWQ,nx_d,NXI,NXQ,
     '  DNUDXQ,DXDXIQ,GCHQ,GUQ,PG,SE,XA,XE,XG,XGRC,XIG,XP,XQ,
     '  YQ,ZA,ZE,ZG,ZP,DEFORMED,NOINITIAL,ACTION,ERROR,*)

C#### Subroutine: IPGRID
C###  Description:
C###    IPGRID defines collocation grids na=1,NMGT (NMGT>1 for 
C###    multigrid if ITYP4(nr,nx)=4).  Steps through all elements, then 
C###    through all gauss points in each element, defining a grid point 
C###    for each unique gauss point (as gauss points on the boundary 
C###    should only have one grid point, though they are shared between 
C###    two elements.  After this, defines grid points to use for 
C###    calculating no-flux boundary conditions at grid points on the 
C###    external boundary of the mesh.  Now also stores values of 
C###    dx/dxi and dnu/dx for each gp in DXQ to be used later in the 
C###    solution process.  Coordinate position (in rectangular 
C###    cartesian) is now stored in XQ, including fibre angle and
C###    sheet angle.
C###
C###    Now also reads and writes YQ information for saving current 
C###    grid state.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b12.cmn'
      INCLUDE 'cmiss$reference:call00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:defn00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:grid00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
      INCLUDE 'cmiss$reference:loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NGAP(NIM,NBM),NHE(NEM,NXM),NJE(NEM),
     '  NKE(NKM,NNM,NBFM,NEFM),NPF(15,NFM),NPNE(NNM,NBFM,NEFM),
     '  NQE(NSM,NBFM,NEFM),NQGE(NGM,NEFM,NBM),nr,NRE(NEM),
     '  NVHE(NNM,NBFM,NHM,NEFM),NVJE(NNM,NBFM,NJM,NEFM),NVQ(NQM,NAM),
     '  NW(NEM,2),NWQ(6,0:NQM,NAM),nx_d,
     '  NXI(-NIM:NIM,0:4,0:NEM),NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      REAL*8 DNUDXQ(3,3,NQM),DXDXIQ(3,3,NQM),GCHQ(3,NQM),GUQ(3,3,NQM),
     '  PG(NSM,NUM,NGM,NBM),SE(NSM,NBFM,NEFM),
     '  XA(NAM,NJM,NQM),XE(NSM,NJM),XG(NJM,NUM),XGRC(NJM,NUM),
     '  XIG(NIM,NGM,NBM),XP(NKM,NVM,NJM,NPM),XQ(NJM,NQM),
     '  YQ(NQM,NIQM,NAM),ZA(NAM,NHM,NCM,NEM),
     '  ZE(NSM,NHM),ZG(NHM,NUM),ZP(NKM,NVM,NHM,NPM,NCM)
      LOGICAL DEFORMED,NOINITIAL
      CHARACTER ACTION*(*),ERROR*(*)
!     Local Variables
      INTEGER IEDGE,IJ,IK,ISUM,IUPDATE,
     '  ME,me_adjac,MG,MGI,mq,mq1,mq2,mq3,mqq,na,nb,nbb,nb_extended,
     '  nc,ndir,ne,ne_adjac,ng,ng_temp,ng1,ng2,NGI(3),ngA,ngFIRST,
     '  ngLAST,ngSTEP,NGTB,nh,nhx,ni,ni1,ni2,ni3,
     '  nii,nij,nik,niq,NITB,nj,nk,noelem,NOFFSET,NONZERO,npos,
     '  nq,nqq,nstep,nu,NUM_ADJ,NU1(3),nxielem,NXN(3),SAME_XI_1D
      REAL*8 A_VECTOR(3),B_VECTOR(3),C_VECTOR(3),
     '  CHTOFF(3,3,3),DBM(3,3,3),DET,
     '  DXIDXI(9,3),DXIX(3,3),GL(3,3),GU(3,3),SUM,
     '  X3G(4,3),XGDIFF,XI(3)
      LOGICAL COLLAPSED,INTERNAL,MU
      CHARACTER ROW*6
      DATA NU1/2,4,7/

      CALL ENTERS('IPGRID',*9999)
      nc=1

      CALL ASSERT(CALL_EQUA,
     '  'Equation must be defined first',ERROR,*9999)

      CALL ASSERT(NIM.EQ.3,'>>NIM must =3',ERROR,*9999)
      CALL ASSERT(NAM.GE.3,'>>NAM must be >= 3',ERROR,*9999)

      DO na=1,NAM
        DO niq=1,NIQM
          DO nq=1,NQM
            YQ(nq,niq,na)=0.0d0
          ENDDO
        ENDDO
      ENDDO

      IF(ACTION(1:4).EQ.'FILE') THEN !read/write grid info

C       Solution time
        IF(IOTYPE.EQ.2) THEN
          READ(IFILE,'(D16.8)') T_RESTART
        ELSE IF(IOTYPE.EQ.3) THEN
          WRITE(IFILE,'(D16.8)') T_RESTART
        ENDIF
C       Number of grid points
        IF(IOTYPE.EQ.2) THEN
          READ(IFILE,'(I6)') nq
          CALL ASSERT(nq.EQ.NQT,'>>NQT does not match',ERROR,*9999)
        ELSE IF(IOTYPE.EQ.3) THEN
          WRITE(IFILE,'(I6)') NQT
        ENDIF
C       Number of levels
        IF(IOTYPE.EQ.2) THEN
          READ(IFILE,'(I6)') na
          CALL ASSERT(na.EQ.NAM,'>>NAM does not match',ERROR,*9999)
        ELSE IF(IOTYPE.EQ.3) THEN
          WRITE(IFILE,'(I6)') NAM
        ENDIF
C       Number of grid point indices
        IF(IOTYPE.EQ.2) THEN
          READ(IFILE,'(I6)') niq
          CALL ASSERT(niq.EQ.NIQM,'>>NIQM does not match',ERROR,*9999)
        ELSE IF(IOTYPE.EQ.3) THEN
          WRITE(IFILE,'(I6)') NIQM
        ENDIF
C       Grid point information
        DO nq=1,NQT
          DO na=1,NAM
            IF(IOTYPE.EQ.2) THEN
              IF(na.EQ.2.AND.NOINITIAL) THEN
                ! Don't read in old initial conditions
                READ(IFILE,'(10D16.8)') (SUM,niq=1,NIQM)
              ELSE
                READ(IFILE,'(10D16.8)') (YQ(nq,niq,na),niq=1,NIQM)
              ENDIF
            ELSE IF(IOTYPE.EQ.3) THEN
              WRITE(IFILE,'(10D16.8)') (YQ(nq,niq,na),niq=1,NIQM)
            ENDIF
          ENDDO !na
        ENDDO !nq

      ELSE IF(ACTION(1:4).EQ.'CALC') THEN !calculate grid points
!Initialize real arrays setup here and auxiliary arrays

        IF(.NOT.ADD) THEN !Meshing first region

C$DOACROSS local(nq,ni,nj),
C$&        share(GUQ,DXDXIQ,DNUDXQ,GCHQ,NQM)
          DO nq=1,NQM
            DO ni=1,3
              DO nj=1,3
                GUQ(nj,ni,nq)=0.d0
                DXDXIQ(nj,ni,nq)=0.d0
                DNUDXQ(nj,ni,nq)=0.d0
              ENDDO
              GCHQ(ni,nq)=0.d0
              XQ(ni,nq)=0.0d0
            ENDDO
          ENDDO
        
          DO ni=1,3
            DO nj=1,3
              GU(nj,ni)=0.d0
              GL(nj,ni)=0.d0
              DO nk=1,3
                CHTOFF(nk,nj,ni)=0.d0
                DBM(nk,nj,ni)=0.d0
              ENDDO
            ENDDO
            DO nj=1,4
              X3G(nj,ni)=0.d0
            ENDDO
          ENDDO

!Initialize NVQ,NWQ & NXQ
C$DOACROSS local(na,nq,ni),
C$&        share(NXQ,NVQ,NMGT,NQM)
          DO na=1,NMGT
            DO nq=0,NQM
              IF(nq.GT.0) NVQ(nq,na)=0 !not dimensioned from zero
              DO ni=-3,3
                NXQ(ni,0,nq,na)=1 
                NXQ(ni,1,nq,na)=0
                NXQ(ni,2,nq,na)=0
                NXQ(ni,3,nq,na)=0
                NXQ(ni,4,nq,na)=0
              ENDDO
            ENDDO !nq
          ENDDO !na
        
C$DOACROSS local(nq,ni),
C$&        share(NWQ,NQM)
          DO nq=1,NQM
!Don't erase everything because of possibility of deforming mesh
!Need to keep at least NWQ(4:5,nq,1) for solution
            DO ni=1,3
              NWQ(ni,nq,1)=0
            ENDDO
            DO ni=1,6
              NWQ(ni,nq,2)=0
            ENDDO
          ENDDO !nq
        ENDIF !not add
      
        nxielem=NIT(NBJ(1,NEELEM(1,nr)))
!Find basis function using Extended Lagrange 
        nb_extended=1
        DO nbb=1,NBT
          IF((nxielem.NE.NIT(nb_extended)).OR.(NBC(nb_extended).NE.7))
     '      THEN
            nb_extended=nb_extended+1
          ENDIF
        ENDDO

C MLB swapped for above loop with modified conditions
C        DO WHILE (nb_extended.LE.NBT.AND.NBC(nb_extended).NE.7)
C          nb_extended=nb_extended+1
C        ENDDO

        CALL ASSERT((NBC(nb_extended).EQ.7),
     '    'Extended basis function not defined',ERROR,*9999)
        WRITE(OP_STRING,'('' Extended basis is #'',I2)') nb_extended
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

!Find element adjacency
        CALL NENXI(IBT,INP,NBJ,NEELEM,NPNE,NXI,ERROR,*9999)

!     Calculate dxi(local)/dxi(global) at each gauss pt (assuming all
!     elements have the same basis)  i.e. relation of local quadratic
!     element to global geometry element.
C  Doesn't work for zero length elements (e.g. at apex) because the zero
C  distance is set to 1.0
        NITB=NIT(nb_extended)       
        DO ni=1,NITB
          ni1=0
          ni2=0
          ni3=0
          IF(ni.GE.1) ni1=1
          IF(ni.GE.2) ni2=1
          IF(ni.GE.3) ni3=1
          ngFIRST=1
          nb=nb_extended
          ngLAST=MAX(NGAP(1,nb),ni2*NGAP(2,nb)*NGAP(1,nb),
     '      ni3*NGAP(3,nb)*NGAP(2,nb)*NGAP(1,nb))
          ngSTEP=MAX(ni1,ni2*NGAP(1,nb),ni3*NGAP(2,nb)*NGAP(1,nb))
          
!     e.g. for 7x7 2d:  ni=1: do ng=1,7,1,  ni=2: do ng=1,49,7
          DO ng=ngFIRST,ngLAST,ngSTEP
!     This is the gauss point number in this direction (1..ngap(ni,nb))
            ngA=((ng-1)/ngSTEP)+1
!     Step forward and back one in this direction
            ng1=ng+ngSTEP
            IF(ng1.GT.ngLAST) ng1=ngFIRST+ngSTEP
            ng2=ng-ngSTEP
            IF(ng2.LT.ngFIRST) ng2=ngLAST-ngSTEP
!     dxi/dxi is xi distance between the two adjacent points in ni dirn
            DXIDXI(NGA,ni)=XIG(ni,ng1,nb)-XIG(ni,ng2,nb)
            IF(DXIDXI(ngA,ni).LT.0.0d-6)
     '        DXIDXI(ngA,ni)=DXIDXI(ngA,ni)+1.0d0
          ENDDO !ng
        ENDDO !ni
        IF(DOP) THEN
          WRITE(OP_STRING,'('' dxi/dxi:'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          DO ni=1,NITB
            WRITE(OP_STRING,'(9F12.5)')
     '        (DXIDXI(ng,ni),ng=1,NGAP(ni,nb))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDIF
      
C$DOACROSS local(noelem,ne,nb,ng),
C$&        share(NEELEM,NGT,NQGE,nb_extended,nr)
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          DO nb=1,NBT
            DO ng=1,NGT(nb)
              NQGE(ng,ne,nb)=0
            ENDDO !nb
          ENDDO !ng
        ENDDO !noelem

        NQR(1,nr)=0
        NQR(2,nr)=0

        IF(ADD) THEN
          nq=NQT
          NQR(1,nr)=NQT+1
        ELSE
          nq=0
          NQR(1,nr)=1
        ENDIF

!     Do for all elements : Define grid points
!!!   CANNOT actually parallelize this loop due to dependencies - but
!!!   this loop is not significant timewise in the subroutine.
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          nb=nb_extended
          NITB=NIT(nb)
          NGTB=NGT(nb)
        
          CALL XPXE(NBJ(1,ne),NKE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '      NQE(1,1,ne),nr,NVJE(1,1,1,ne),SE(1,1,ne),
     '      XA,XE,XP,ERROR,*9999)

          IF(DEFORMED) THEN
            CALL ZPZE(NBH(1,1,ne),nc,NHE(ne,nx_d),NKE(1,1,1,ne),
     '        NPF(1,1),NPNE(1,1,ne),NRE(ne),NVHE(1,1,1,ne),
     '        NW(ne,1),nx_d,SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
          ENDIF

!     Do for all gauss pts within element
          DO ng=1,NGTB
            IF(DOP) THEN
              WRITE(OP_STRING,
     '          '('' Element '',I4,'' - Gauss pt '',I4)') ne,ng
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
!     NGAP(i,nb) is number of g.p.s there are along Xi(i) in basis nb
!     NGI(i) is this g.p. number in dirn Xi(i) in element
!      - between (1..ngap(i,nb))
            NGI(1)=MOD((ng-1),NGAP(1,nb))+1  
            IF(NITB.GE.2) THEN
              NGI(2)=MOD(INT((ng-1)/NGAP(1,nb)),NGAP(2,nb))+1
              IF(NITB.EQ.3) THEN
                NGI(3)=MOD(INT((ng-1)/(NGAP(1,nb)*NGAP(2,nb))),
     '            NGAP(3,nb))+1
              ELSE
                NGI(3)=0
              ENDIF
            ELSE
              NGI(2)=0
              NGI(3)=0
            ENDIF !nitb
!     Find if g.p. is on edge of elem
            IEDGE=0
            DO ni=1,3
!     NXN contains edge that gp is on
              IF(ni.LE.NITB) THEN   
                IF(NGI(ni).EQ.1) THEN
                  NXN(ni)=-ni
                  IEDGE=1
                ELSE IF(NGI(ni).EQ.NGAP(ni,nb)) THEN
                  NXN(ni)=ni
                  IEDGE=1
                ELSE
                  NXN(ni)=0
                ENDIF
              ELSE
                NXN(ni)=0
              ENDIF
            ENDDO !ni
            IUPDATE=0
            IF(IEDGE.EQ.1) THEN
!     If on an edge, then check whether it shares that edge
!     with an element for whom g.p. have already been
!     defined. (ie ME less than NE)
              DO nik=MIN(NXN(3),0),MAX(NXN(3),0),3
                DO nij=MIN(NXN(2),0),MAX(NXN(2),0),2
                  DO nii=MIN(NXN(1),0),MAX(NXN(1),0)
                    IF(NITB.eq.1) THEN
                      NUM_ADJ=NXI(nii,0,ne)
                    ELSE IF(NITB.eq.2) THEN
                      NUM_ADJ=1 
C NOTE: if elements branch in 2D then consistent Xi directions are
C assumed with in the 3D host mesh, if this is not the case this
C code will have to be extended NPS 14/2/97
                    ELSE IF(NITB.eq.3) THEN
                      NUM_ADJ=1
                    ENDIF
                    DO ne_adjac=1,NUM_ADJ
                      ME=NXI(nik,1,NXI(nij,1,NXI(nii,ne_adjac,ne))) 
C Neighbour elem
 !     If exists and g.p. already defined in that element
                      IF(ME.LT.NE.AND.ME.GT.0) THEN 
                        SAME_XI_1D=1
                        IF(NITB.eq.1) THEN
                          DO me_adjac=1,NXI(nii,0,ME)
                            IF(ne.EQ.NXI(nii,me_adjac,ME)) THEN
                              SAME_XI_1D=0
                            ENDIF
                          ENDDO
                        ENDIF
C This loop is added to test wiether their is consistant Xi directions
C between elements in 1D (it may at some stage need to be extended to
C 2D & 3D. If the Xi directions are opposite in 1D then the local gauss
C point number MG=nq
                        MG=ng - nii*(NGAP(1,nb)-1)*SAME_XI_1D
     '                    - nij/2*(NGAP(2,nb)-1)*NGAP(1,nb)
     '                    - nik/3*(NGAP(3,nb)-1)*NGAP(2,nb)*NGAP(1,nb)
                        NQGE(ng,ne,nb)=NQGE(MG,ME,nb) !Use global gp
                        IUPDATE=1 !already defined
                      ENDIF
                    ENDDO
                  ENDDO !nii
                ENDDO !nij
              ENDDO !nik
            ENDIF !iedge=1
            DO ni=1,NITB
              XI(ni)=XIG(ni,ng,nb_extended)
            ENDDO
            IF(DEFORMED) THEN
C           Interpolate dependent var.s ZG and derivs wrt Xi (JP=0)
              CALL ZEZW(0,IBT,IDO,INP,NAN,NBH(1,1,ne),NHE(ne,nx_d),
     '          nr,nx_d,DXIX,ZE,ZG,XI,ERROR,*9999)
              DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
                nj=NJ_LOC(NJL_GEOM,nhx,nr)
                nh=NH_LOC(nhx,nx_d)
                DO nu=1,NUT(NBH(nh,nc,ne))
                  XG(nj,nu)=ZG(nhx,nu)
                ENDDO !nu
              ENDDO !nh/nj
            ELSE
C           Interpolate geometric vars XG and derivs wrt Xi
              CALL XEXW(IBT,IDO,INP,NAN,NBJ(1,ne),nr,XE,XG,XI,
     '          ERROR,*9999)
            ENDIF
C 25-APR-1996 GBS Check if g.p. is on a collapsed edge
C    Only works for r.c. I think
            COLLAPSED=.FALSE.
            DO ni=1,NITB
              NONZERO=0
              DO nj=1,NJT
                IF(DABS(XG(nj,NU1(ni))).GT.1.d-10) NONZERO=NONZERO+1
              ENDDO
              COLLAPSED=COLLAPSED.OR.(NONZERO.EQ.0)
            ENDDO
C 26-APR-1996 GBS If mu=0 for prolate, then do not compute grid point
            IF(ITYP10(nr).EQ.4.AND.DABS(XG(2,1)).LT.1.d-10) THEN
              MU=.TRUE. !prolate AND mu=0
            ELSE
              MU=.FALSE. !not prolate or mu<>0
            ENDIF

!     If the grid point has not already been defined, 
!     and not on a collapsed edge or mu=0 for prolate, then define a
!     new grid point nq
            IF(IUPDATE.EQ.0.AND. !gp not already defined
     '        .NOT.COLLAPSED.AND. !not on a collapsed edge
     '        .NOT.MU) THEN !mu<>0
              nq=nq+1
              IF(DOP) THEN
                WRITE(OP_STRING,'('' nq = '',I5)') nq
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              CALL ASSERT(nq.LE.NQM,'Too many grid points:'
     '          //' increase NQT',ERROR,*9999)
              NQGE(ng,ne,nb)=nq
              IF(DOP) THEN
                WRITE(OP_STRING,'('' Coordinate position:'')')
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(3F12.3)')(XG(nj,1),nj=1,NJT)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              
!     Store position of nq (convert to r.c.)
              CALL XZ_DERIV(ITYP10(nr),1,XG,XGRC)
              DO nj=1,NJT
                XQ(nj,nq)=XGRC(nj,1)
              ENDDO
              IF(DOP) THEN
                WRITE(OP_STRING,'('' Coordinate position (r.c.):'')')
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(3F12.3)')(xq(nj,nq),nj=1,njt)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF

!     Store fibre and sheet angles      
C GBS 26-7-96  Removed because these are not used
c            DO nj=1,NJ_LOC(NJL_FIBR,0,nr)
c              IF(DEFORMED) THEN
C Must have previously called "update gauss deformed_fibres collocation"
c                XQ(njt+nj,nq)=YG(ng,nj,ne)
c              ELSE
c                XQ(njt+nj,nq)=XG(NJ_LOC(NJL_FIBR,nj,nr),1)
c              ENDIF
c            ENDDO
c      	    IF(DOP) THEN
c              WRITE(OP_STRING,'('' Fibre and sheet angles:'')')
c              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c              WRITE(OP_STRING,'(3F12.3)')(xq(nj,nq),nj=njt+1,njt+3)
c              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c            ENDIF

!     Compute and store dx/dxi      
!     Multiply by dxi/dxi in order to adjust for local quadratic element
              DO ni=1,NITB
                CALL XZ_DERIV(ITYP10(nr),NU1(ni),XG,XGRC)
                DO nj=1,NJT
                  DXDXIQ(nj,ni,nq)=XGRC(nj,NU1(ni))*DXIDXI(NGI(ni),ni)
                ENDDO
              ENDDO
              IF(DOP) THEN
                WRITE(OP_STRING,'('' dx/dxi:'')')
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                DO ni=1,NJT
                  WRITE(OP_STRING,'(3F12.3)')
     '              (DXDXIQ(nj,ni,nq),nj=1,njt)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDDO
              ENDIF

!     Compute and store Gij (upper)
!     GL (=g(lower ij)) = dx(k)/dxi(i)*dx(k)/dxi(j)      
              DO ni2=1,NITB
                DO ni1=1,NITB
                  GL(ni1,ni2)=0.0D0
                  DO nj=1,NJT
                    GL(ni1,ni2)=GL(ni1,ni2)+
     '                DXDXIQ(nj,ni1,nq)*DXDXIQ(nj,ni2,nq)
                  ENDDO
                ENDDO
              ENDDO
              IF(DOP) THEN
                WRITE(OP_STRING,'('' GL'')')
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                DO nj=1,NITB
                  WRITE(OP_STRING,'(3F12.5)') (gl(ni,nj),ni=1,nitb)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDDO
              ENDIF
              
!     	    GU = inv(GL)
              IF(NITB.GE.2) THEN
                CALL INVERT(NITB,GL,GU,DET)
              ELSE IF(GL(1,1).NE.0.0D0) THEN
                GU(1,1)=1.0D0/GL(1,1)
              ELSE
                WRITE(OP_STRING,*) 
     '            ' >>>Zero value for GL - cannot invert'
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
              IF(DOP) THEN
                WRITE(OP_STRING,'('' GU'')')
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                DO nj=1,NITB
                  WRITE(OP_STRING,'(3F12.5)') (GU(ni,nj),ni=1,NITB)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDDO
              ENDIF
              DO ni=1,3
                DO nj=1,3
C Temporary hack GBS 2-10-96
                  if(gu(ni,nj).gt.5000.0d0) then
                    write (*,'(5i7,f20.12)') ne,ng,nq,ni,nj,gu(ni,nj)
                    gu(ni,nj)=5000.0d0
                  endif
C End of temporary hack
                  GUQ(ni,nj,nq)=GU(ni,nj)
                ENDDO
              ENDDO
       
!     Compute Christoffels CHTOFF and product with GU (GCHQ)
              CALL TOFFEL(nb,NJE(ne),nr,CHTOFF,DBM,GU,XG,X3G,.FALSE.,
     '          ERROR,*9999)
              DO nik=1,NITB
                SUM=0.D0
                DO nii=1,NITB
                  DO nij=1,NITB
                    SUM=SUM+CHTOFF(nik,nii,nij)*GU(nii,nij)
     '                *DXIDXI(NGI(nii),nii)*DXIDXI(NGI(nij),nij)
                  ENDDO
                ENDDO
                GCHQ(nik,nq)=SUM*DXIDXI(NGI(nik),nik)
              ENDDO
              IF(DOP) THEN
                WRITE(OP_STRING,'('' GCHQ(nik=1..,nq)'',3E12.5)') 
     '            (GCHQ(nik,nq),nik=1,NITB)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF

!     Compute and store dnu/dx (dirn cosines of material coords)
              IF(DEFORMED) THEN
                ng_temp=0 ! to compute at xi position (computed above)
                CALL MAT_VEC_DEF(IBT,IDO,INP,NAN,NBH,NBJ(1,ne),ng_temp,
     '            NHE(ne,nx_d),nr,nx_d,A_VECTOR,B_VECTOR,C_VECTOR,
     '            PG,XE,XG,XI,ZE,ZG,ERROR,*9999)
              ELSE
                CALL MAT_VEC_XI(IBT,IDO,INP,NAN,NBJ(1,ne),nr,
     '            A_VECTOR,B_VECTOR,C_VECTOR,XE,XG,XI,ERROR,*9999)
              ENDIF
              DO nj=1,3
                DNUDXQ(1,nj,nq)=A_VECTOR(nj)          
                DNUDXQ(2,nj,nq)=B_VECTOR(nj)
                DNUDXQ(3,nj,nq)=C_VECTOR(nj)
              ENDDO
              IF(DOP) THEN
                WRITE(OP_STRING,'('' dnu/dx:'')')
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                DO ni=1,3
                  WRITE(OP_STRING,'(3F12.3)')(DNUDXQ(ni,nj,nq),nj=1,3)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDDO
              ENDIF
                
            ENDIF !iupdate (grid pt not defined)
          ENDDO !ng
        ENDDO !noelem
        NQT=nq
        NQR(2,nr)=NQT

!   Compute connectivity of grid points
!   NWQ(ni,nq,2) used as temporary storage of neighbouring points
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          nb=nb_extended
          NITB=NIT(nb)
          NGTB=NGT(nb)
          DO ng=1,NGTB
            nq=NQGE(ng,ne,nb)
            IF(DOP) THEN
              WRITE(OP_STRING,'(/'' ne,ng,nq '',3I8)') ne,ng,nq
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
              
            IF(nq.GT.0) THEN
! NXQ stores position of neighbouring grid point.
! Loop through Gauss pts ng (->nq) of current element & establish
! mapping to neighbours mq both ways (ie NXQ(nq)=mq & NXQ(mq)=nq). 
! If the 2-way mapping has already been established then do nothing
! but if find a 1-way mapping then nq must have multiple neigbours
! on one side (i.e. the element Xi line splits in two). 
                                      
              NXQ(0,1,nq,1)=nq
              NGI(1)=MOD((NG-1),NGAP(1,nb))+1  
              IF(NITB.GE.2) THEN
                NGI(2)=MOD(INT((NG-1)/NGAP(1,nb)),NGAP(2,nb))+1
                IF(NITB.EQ.3) THEN
                  NGI(3)=MOD(INT((NG-1)/(NGAP(1,nb)*NGAP(2,nb))),
     '              NGAP(3,nb))+1
                ELSE
                  NGI(3)=0
                ENDIF
              ELSE
                NGI(2)=0
                NGI(3)=0
              ENDIF
              NOFFSET=1
              DO ni=1,NITB
                DO ndir=-1,1,2
                  MGI=NGI(ni)+ndir
                  IF(MGI.GT.0.AND.MGI.LE.NGAP(ni,nb)) THEN
                    MG=NG+ndir*NOFFSET
! Neighbouring grid pt mq
                    mq=NQGE(MG,ne,nb)
                    IF(DOP) THEN
                      WRITE(OP_STRING,'('' ni,ndir,mq '',3I8)')
     '                  ni,ndir,mq
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF

                    IF(mq.GT.0) THEN
                      mqq=NXQ( ndir*ni,1,nq,1) !forward mapping
                      nqq=NXQ(-ndir*ni,1,mq,1) !reverse mapping

                      IF((mqq.EQ.0.AND.nqq.NE.0).OR
     '                  .(nqq.EQ.0.AND.mqq.NE.0)) THEN !1-way mapping 
                        IF(DOP) THEN
                          WRITE(OP_STRING,'('' found 1-way mapping'
     '                      //' at nq='',I6,'' mq='',I6)') nq,mq
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                        ENDIF
C MLB 1 October 1996
C Changed to use multiple branching (>2)
C                        IF(mqq.EQ.0) THEN !1st slot for nq & 2nd for mq 
C                          NXQ( ndir*ni,1,nq,1)=mq
C                          NWQ(ni,nq,2)=NWQ(ni,nq,2)+ndir
C                          NXQ(-ndir*ni,0,mq,1)=2
C                          NXQ(-ndir*ni,2,mq,1)=nq
C                        ELSE IF(nqq.EQ.0) THEN !1st for nq & 2nd for mq 
C                          NXQ(-ndir*ni,1,mq,1)=nq
C                          NWQ(ni,mq,2)=NWQ(ni,mq,2)-ndir
C                          NXQ( ndir*ni,0,nq,1)=2
C                          NXQ( ndir*ni,2,nq,1)=mq
C                        ENDIF                      
                    
                        IF(mqq.EQ.0) THEN !1st slot for nq & 2nd for mq 
                          NXQ( ndir*ni,1,nq,1)=mq
                          NWQ(ni,nq,2)=NWQ(ni,nq,2)+ndir
                          NXQ(-ndir*ni,0,mq,1)=NXQ(-ndir*ni,0,mq,1)+1
                          NXQ(-ndir*ni,NXQ(-ndir*ni,0,mq,1),mq,1)=nq
                        ELSE IF(nqq.EQ.0) THEN !1st for nq & 2nd for mq 
                          NXQ(-ndir*ni,1,mq,1)=nq
                          NWQ(ni,mq,2)=NWQ(ni,mq,2)-ndir
                          NXQ( ndir*ni,0,nq,1)=NXQ( ndir*ni,0,nq,1)+1
                          NXQ( ndir*ni,NXQ( ndir*ni,0,nq,1),nq,1)=mq
                        ENDIF                                          

                      ELSE !no mapping or 2-way mapping
                        IF(mqq.EQ.0) THEN
                          NXQ( ndir*ni,1,nq,1)=mq
                          NWQ(ni,nq,2)=NWQ(ni,nq,2)+ndir
                        ENDIF
                        IF(nqq.EQ.0) THEN
                          NXQ(-ndir*ni,1,mq,1)=nq
                          NWQ(ni,mq,2)=NWQ(ni,mq,2)-ndir
                        ENDIF                                          
                      ENDIF !mqq,nqq

c                 ELSE
c                   WRITE(OP_STRING,'('' Error: mq < 0 '')')
c                   CALL WRITES(IOER,OP_STRING,ERROR,*9999)
                    ENDIF
                  ENDIF !mgi
                ENDDO !ndir
                NOFFSET=NOFFSET*NGAP(ni,nb)
              ENDDO !ni
            ENDIF !nq
          ENDDO !ng
        ENDDO !ne

C branching of 3d elements may not work correctly
C news AJP 25/7/96.  It seems that nowhere is NXQ(ni,0,nq,na) put to
C zero if NXQ(ni,1,nq,na) is zero.  This is the case for a 2d heart
C cross section.
        DO ni=-3,3
          DO nq=NQR(1,nr),NQR(2,nr)
            DO na=1,NMGT
              IF(NXQ(ni,1,nq,na).EQ.0) NXQ(ni,0,nq,na)=0
!           do I need to check NXQ(ni,2..4,nq,na) as well??
            ENDDO !na
          ENDDO !nq
        ENDDO !ni
C newe AJP 25/7/96

! Compute NQGE for non-extended basis functions of the same order
! Stores closest grid point to the gauss point for that basis
        DO nb=1,NBT
          IF(nb.NE.nb_extended.AND.NIT(nb).EQ.NIT(nb_extended)) THEN
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
C loop over gauss pts for mechanics basis
              DO ng=1,NGT(nb) 
                ng1=1
C find gauss pt in extended basis to lower left of ng in xi-space
                DO ni=NIT(nb),1,-1 
                  npos=INT(XIG(ni,ng,nb)*(NGAP(ni,nb_extended)-1))+1
                  ng1=(ng1-1)*NGAP(ni,nb_extended)+npos
                ENDDO
C associated grid point is nq
                nq=NQGE(ng1,ne,nb_extended)
C check to see if ng is closer to any of the neighbouring gauss pts
C upwards in xi-space
                DO ni=1,NIT(nb)
                  XGDIFF=(XIG(ni,ng,nb)-XIG(ni,ng1,nb_extended))*
     '              DBLE(NGAP(ni,nb_extended)-1)
                  IF(XGDIFF.GE.0.5d0) nq=NXQ(ni,1,nq,1)
                ENDDO !ni
                NQGE(ng,ne,nb)=nq
              ENDDO !ng
            ENDDO !ne
          ENDIF !non-extended basis
        ENDDO !nb

!     Store g.p.s to use for calculating noflux b.c.
        IF(DOP) THEN
          WRITE(OP_STRING,'(/'' Boundary point checking:'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        DO nq=NQR(1,nr),NQR(2,nr)
!       Check if on external bdy by checking all neighbouring g.p.
          INTERNAL=.TRUE.
          ISUM=0
          IK=MAX(0,NITB-2) !zero for 1,2-D, one for 3-D
          IJ=MIN(NITB-1,1) !zero for 1-D, one for 2,3-D
          DO nik=-IK,IK
            DO nij=-IJ,IJ
              DO nii=-1,1
                mq=NXQ(nii,1,NXQ(nij*2,1,NXQ(nik*3,1,nq,1),1),1)
                IF(mq.EQ.0) THEN
                  INTERNAL=.FALSE.
                  ISUM=ISUM+1
                ENDIF
              ENDDO !nii
            ENDDO !nij
          ENDDO !nik
          IF(INTERNAL) THEN !internal g.p.
            DO ni=1,2
              NWQ(ni,nq,1)=0
            ENDDO
          ELSE !on external bdy
            IF(NITB.EQ.2) THEN
              IF(ISUM.EQ.1) THEN
!             Special case for "270 degree" corner (internal corner)
!             Find which g.p. is missing, and head away from that.
                DO nij=-1,1
                  DO nii=-1,1
                    mq=NXQ(nii,1,NXQ(nij*2,1,nq,1),1)
                    IF(mq.EQ.0) THEN
                      NWQ(1,nq,2)=-nii
                      NWQ(2,nq,2)=-nij
                    ENDIF
                  ENDDO
                ENDDO
              ENDIF
            ELSE IF(NITB.EQ.3) THEN
            ENDIF
!         Now find g.p.s to use for no-flux b.c.s
            mq=nq
            DO nstep=1,2
 !           Step twice in dirn indicated from NWQ(ni,nq,2)
              DO ni=1,NITB
                mq=NXQ(ni*NWQ(ni,nq,2),1,mq,1)
              ENDDO
              NWQ(nstep,nq,1)=mq
            ENDDO
          ENDIF !int/ext 
        ENDDO !nq
        DO nq=NQR(1,nr),NQR(2,nr)
          IF(DOP) THEN
            WRITE(OP_STRING,'(3I6)') nq,(NWQ(ni,nq,1),ni=1,2)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
!     	Set temporary array to zero
          DO ni=1,3
            NWQ(ni,nq,2)=0
          ENDDO
        ENDDO !nq

! Calculate NVQ(nq,na) array
C     nq=1
C     DO WHILE(NXQ(1,1,1,nq).NE.0) !calc #nq's in 1st row
C       nq=nq+1
C     ENDDO
C     NQ_ROW1=nq
C     IF(DOP) THEN 
C       WRITE(OP_STRING,'('' #nq in 1st row = '',I4)') NQ_ROW1
C       CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C     ENDIF
C     DO na=1,NMGT
C       DO nq=1,NQT
C         NVQ(na,nq)=-1 !initialize to not in grid 
C       ENDDO
C       IF(DOP) THEN
C         WRITE(OP_STRING,'(/'' Multigrid level na='',I5)') na
C         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C       ENDIF
C       nq=1-2**(na-1)
c       DO WHILE(nq.LE.0.OR.NXQ(1,1,3,nq+1).NE.0)
C         DO WHILE(nq.LE.0.OR.NXQ(1,1,2,nq).NE.0.OR.NXQ(1,1,1,nq).NE.0)
C           nq=nq+2**(na-1) !step nq in increments of 2**(na-1)
C           IF(DOP) THEN
C             WRITE(OP_STRING,'('' nq='',I5,' is 1st nq in row'')')nq
C             CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C           ENDIF
C           IF(MOD(nq,2**na).EQ.1) THEN !row has g.p.s at level na+1
C             IF(DOP) THEN
C               WRITE(OP_STRING,'('' row has g.p.s at level na+1'')')
C               CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C             ENDIF
C             NVQ(na,nq)=0
C             IF(DOP) THEN
C               WRITE(OP_STRING,'('' NVQ(na,'',I4,'')='',I1)')
C    '            nq,NVQ(na,nq)
C               CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C             ENDIF
C             Row1_tot=1
C             DO WHILE(NXQ(1,1,1,nq).NE.0)  !loop over Xi(1)
C               nq=nq+2**(na-1)
C               Row1_tot=Row1_tot+1
C               IF(MOD(nq,2**na).EQ.1) THEN
C                 NVQ(na,nq)=0 !nq coincides with coarse g.p.
C               ELSE
C                 NVQ(na,nq)=1 !nq lies between 2 coarse g.p.s in Xi(1)
C               ENDIF
C               IF(DOP) THEN
C                 WRITE(OP_STRING,'('' NVQ(na,'',I4,'')='',I1)') 
C    '              nq,NVQ(na,nq)
C       	  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C               ENDIF
C             ENDDO
C             nq=nq+(NQ_ROW1-1)*(2**(na-1)-1) !skip rows
C             IF(DOP) THEN
C               WRITE(OP_STRING,'('' Row1_tot='',I3,'' skip'',I4,
C    '            '' Current nq='',I5)') 
C    '            Row1_tot,(NQ_ROW1-1)*(2**(na-1)-1),nq
C       	CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C             ENDIF
C       
C             nq_start=nq+2**(na-1) !is 1st nq in new row
C             IF(nq_start.LT.NQT) THEN
C               !loop over row with g.p.s midway betw pts at level na
C               DO Row1=1,Row1_tot
C                 nq=nq+2**(na-1)
C                 IF(MOD(nq-nq_start,2**na).EQ.0) THEN
C                   NVQ(na,nq)=2 !nq lies betw 2 coarse g.p.s in Xi(2)
C                 ELSE
C                   NVQ(na,nq)=4 !nq lies betw 4 coarse g.p.s in Xi1,2
C                 ENDIF
C                 IF(DOP) THEN
C                   WRITE(OP_STRING,'('' NVQ(na,'',I4,'')='',I1)') 
C    '                nq,NVQ(na,nq)
C                   CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                 ENDIF
C               ENDDO !Row1
C             ENDIF !nq<NQT
C             nq=nq+(NQ_ROW1-1)*(2**(na-1)-1) !skip rows
C             IF(DOP) THEN
C              WRITE(OP_STRING,'(14X,''skip'',I4,'' Curr nq='',I5)') 
C    '            (NQ_ROW1-1)*(2**(na-1)-1),nq
C               CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C             ENDIF
C           ENDIF
C         
C         ENDDO !while
c       ENDDO !while
C     ENDDO !na

! NWQ is used in CONSTRUCTNVQ as a temporary array
        CALL CONSTRUCTNVQ(NITB,NVQ,NXQ,NWQ(1,1,3),ERROR,*9999)

! Calculate NXQ(-3:3,1,nq,na+1) array
        DO na=1,NMGT-1
          DO nq=NQR(1,nr),NQR(2,nr)
            IF(NVQ(nq,na+1).GE.0) THEN !nq is in grid na+1
              NXQ(0,1,nq,na+1)=nq
              DO ni=-1,1
                mq3=NXQ(3*ni,1,nq,na) !Xi3 neighbour at grid level na
                IF(mq3.ne.0) NXQ(3*ni,1,nq,na+1)=NXQ(3*ni,1,mq3,na)
                mq2=NXQ(2*ni,1,nq,na) !Xi2 neighbour at grid level na
                IF(mq2.ne.0) NXQ(2*ni,1,nq,na+1)=NXQ(2*ni,1,mq2,na)
                mq1=NXQ(  ni,1,nq,na) !Xi1 neighbour at grid level na
                IF(mq1.ne.0) NXQ(  ni,1,nq,na+1)=NXQ(  ni,1,mq1,na)
              ENDDO !ni
            ENDIF !nq is in grid na+1
          ENDDO !nq
        ENDDO !na

! Calculate NWQ(ni,nq,na=2..) array
        DO na=2,NMGT
          IF(DOP) THEN
            WRITE(OP_STRING,'('' Calculate NWQ for grid na='',I2)') na
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF                        
          DO nq=NQR(1,nr),NQR(2,nr)
            DO ni=1,6
              NWQ(ni,nq,na)=0
            ENDDO
            IF(NVQ(nq,na).GE.0) THEN !nq is in grid na
              IF(NXQ(-2,1,nq,na).eq.0) THEN !bottom row
                ROW='Lower'
                mq1=NXQ(2,1,nq,na) !is one row  above  
                mq2=NXQ(2,1,mq1,na) !is two rows above  
              ELSE IF( NXQ(2,1,nq,na).eq.0) THEN !top row
                ROW='Upper'
                mq1=NXQ(-2,1,nq,na) !is one row  below
                mq2=NXQ(-2,1,mq1,na) !is two rows below  
              ELSE
                ROW='Centre'
              ENDIF
              IF(ROW(1:5).EQ.'Lower'.OR.ROW(1:5).EQ.'Upper') THEN 
                IF(NXQ(-1,1,nq,na).EQ.0) THEN !LH corner
                  NWQ(1,nq,na)=NXQ(1,1,mq1,na)
                  NWQ(2,nq,na)=NXQ(1,1,NXQ(1,1,mq2,na),na)
                ELSE IF(NXQ(1,1,nq,na).EQ.0) THEN !RH corner
                  NWQ(1,nq,na)=NXQ(-1,1,mq1,na)
                  NWQ(2,nq,na)=NXQ(-1,1,NXQ(-1,1,mq2,na),na)
                ELSE !midside
                  NWQ(1,nq,na)=mq1
                  NWQ(2,nq,na)=mq2
                ENDIF
              ELSE IF(ROW(1:6).EQ.'Centre') THEN
                IF(NXQ(-1,1,nq,na).eq.0) THEN !LH edge
                  NWQ(1,nq,na)=NXQ(1,1,nq,na)
                  NWQ(2,nq,na)=NXQ(1,1,NXQ(1,1,nq,na),na)
                ELSE IF(NXQ(1,1,nq,na).eq.0) THEN !RH edge
                  NWQ(1,nq,na)=NXQ(-1,1,nq,na)
                  NWQ(2,nq,na)=NXQ(-1,1,NXQ(-1,1,nq,na),na)
                ENDIF
              ENDIF
              IF(DOP) THEN
                WRITE(OP_STRING,'(''***** nq='',I6,'' ROW='',A,'
     '            //''' NWQ='',2I6)') nq,ROW,NWQ(1,nq,na),NWQ(2,nq,na)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              
            ENDIF !nq is in grid na
          ENDDO !nq
        ENDDO !na
      ENDIF !r/w/c
        
      CALL EXITS('IPGRID')
      RETURN
 9999 CALL ERRORS('IPGRID',ERROR)
      CALL EXITS('IPGRID')
      RETURN 1
      END 


C 24/2/97 LC removed from :

C#### Subroutine: IPEQUA
C###  Description: 
C###    IPEQUA defines equation for region nr and problem nx.

C old GMH 5/12/95 Change method of asking for element type
C        DO noelem=1,NEELEM(0,nr)
C          ne=NEELEM(noelem,nr)
C          IF(ne.EQ.1) THEN
C            FORMAT='('' Enter element type:'''//
C     '        '/''   (1) Truss or cable      (5) Membrane           '//
C     '        ' ( 9) 3-dimensional'''//
C     '        '/''   (2) Batten              (6) Plate (Kirchhoff)  '//
C     '        ' (10) Tank Bottom '''//
C     '        '/''   (3) Beam  (Kirchhoff)   (7) Shell              '//
C     '        ' (11) Plane stress '''//
C     '        '/''   (4) Link                (8) Shell/Fluid I-face '//
C     '        ' (12) Plane strain '''//
C     '        '/'' '')'
C            CALL GINOUT(IOTYPE,0,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
C     '        ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
C     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C          ENDIF
C          WRITE(CHAR1,'(I4)') ne
C          CALL TRIM(CHAR1,IBEG1,IEND1)
C          IDEFLT(1)=0
C          IF(ne.NE.1) IDEFLT(1)=NW(NE-1,1)
C          WRITE(CHAR2,'(I3)') IDEFLT(1)
C          CALL TRIM(CHAR2,IBEG2,IEND2)
C          FORMAT='($,'' Element '//CHAR1(IBEG1:IEND1)//
C     '      ' ['//CHAR2(IBEG2:IEND2)//']: '',I3)'
C          IF(IOTYPE.EQ.3) IDATA(1)=NW(ne,1)
C          CALL GINOUT(IOTYPE,3,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,12,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C          IF(IOTYPE.NE.3) NW(ne,1)=IDATA(1)
C          IF(NW(ne,1).EQ.4) THEN
CC ***       Link elements: length approx zero so set SE=1 for y-deriv
C            nb=NBJ(2,ne)
C            SE(2,nb,ne)=1.d0
C            SE(4,nb,ne)=1.d0
C          ELSE IF(NW(ne,1).EQ.6) THEN !Plate element s.b. bilin geom
C            nb=NBJ(1,ne)
C            CALL ASSERT(NNT(nb).EQ.4,
C     '        '>>Plate must have bilinear geometry',ERROR,*9999)
C            IF(NJT.EQ.3) THEN
C              CALL ASSERT(JTYP9.GT.0,'>>3D Plate must have fibre field',
C     '          ERROR,*9999)
C              nc=1 !Temporary AJP 17-12-91
C              CALL ASSERT(NHM*NSM.GE.12*NKT(0,NBH(3,nc,ne)),
C     '          '>>NVM is too small',ERROR,*9999)
C            ENDIF
C          ELSE IF(NW(ne,1).EQ.9) THEN
C            CALL ASSERT(NJT.EQ.3,'>>Need 3D coordinates',ERROR,*9999)
C            nb=NBJ(1,ne)
C            CALL ASSERT(NIT(nb).EQ.3,'>>Need 3D basis',ERROR,*9999)
C          ENDIF
C        ENDDO
C endold GMH 5/12/95


C 24/2/97 LC  removed from :

C#### Subroutine: IPFIT
C###  Description:
C###    IPFIT inputs parameters for geometry or field fit for region nr.

C*** See archive version for non-linear and Fourier fitting problems


C GMH 13/3/96 Add #/name to fit stuff
C        CONTINUE=.TRUE.
C        nonode=0
C        DO WHILE(CONTINUE)
C          IDEFLT(1)=0
C          IF(ENTERCOUPLING) THEN
C            FORMAT='(/$,'' Enter node to fix or to specify coupling'
C     '        //' [EXIT]: '',I4)'
C          ELSE
C            FORMAT='(/$,'' Enter node number to fix [EXIT]: '',I4)'
C          ENDIF
C          IF(IOTYPE.EQ.3) THEN
C            nonode=nonode+1
C            IF(nonode.LE.NPNODE(0,nr)) THEN
C              np=NPNODE(nonode,nr)
C              IDATA(1)=np
C            ELSE
C              IDATA(1)=0
C              CONTINUE=.FALSE.
C            ENDIF
C          ENDIF
C 6500     CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,1,NPM,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C          IF(IOTYPE.NE.3) THEN
C            IF(IDATA(1).EQ.0) THEN
C              CONTINUE=.FALSE.
C            ELSE IF(IDATA(1).NE.0) THEN
C              CONTINUE=.TRUE.
C              np=IDATA(1)
C              IF(.NOT.INLIST(np,NPNODE(1,nr),NPNODE(0,nr),N1)) THEN
C                WRITE(OP_STRING,'(A)')
C     '            ' >>This node does not belong to the current region'
C                CALL WRITES(IOER,OP_STRING,ERROR,*9999)
C                GOTO 6500
C              ENDIF
C            ENDIF
C          ENDIF
C          IF(CONTINUE) THEN
C            DO njj=1,NUM_FIT(0)
C              WRITE(CHAR1,'(I10)') njj
C              CALL TRIM(CHAR1,IBEG1,IEND1)
C              DO nhj=1,NUM_FIT(njj)
C                WRITE(CHAR2,'(I1)') nhj
C                nhx=NLH_FIT(nhj,3,njj)
C                nh=NH_LOC(nhx,nx)
C                ADEFLT(1)='N'
C                IF(NVHP(nh,np,1,nr).EQ.1.AND.NKH(nh,np,1,nr).EQ.1) THEN
C                  FORMAT='($,'' Is variable '//CHAR2(1:1)//' of fit '
C     '              //'variable '//CHAR1(IBEG1:IEND1)//' fixed '
C     '              //'[N]?: '',A)'
C                ELSE
C                  FORMAT='($,'' Are any variables for variable '
C     '              //CHAR2(1:1)//' of fit variable '
C     '              //CHAR1(IBEG1:IEND1)//' fixed [N]?: '',A)'
C                ENDIF
C                IF(IOTYPE.EQ.3) THEN
C                  ISFIXED=.FALSE.
C                  DO nv=1,NVHP(nh,np,1,nr)
C                    DO nk=1,NKH(nh,np,1,nr)
C                      ny=NYNP(nk,nv,nh,np,0,1,nr)
C                      IF(FIX(ny,1)) ISFIXED=.TRUE.
C                    ENDDO !nk
C                  ENDDO !nv
C                  IF(ISFIXED) THEN
C                    ADATA(1)='Y'
C                  ELSE
C                    ADATA(1)='N'
C                  ENDIF
C                ENDIF      
C                CALL GINOUT(IOTYPE,1,IVDU,IFILE,0,0,NOQUES,FILEIP,
C     '            FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
C     '            0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
C     '            ERROR,*9999)
C                IF(IOTYPE.NE.3) THEN
C                  IF(ADATA(1).EQ.'Y'.OR.ADATA(1).EQ.'y') THEN
C                    ISFIXED=.TRUE.
C                  ELSE
C                    ISFIXED=.FALSE.
C                    DO nv=1,NVHP(nh,np,1,nr)
C                      DO nk=1,NKH(nh,np,1,nr)
C                        ny=NYNP(nk,nv,nh,np,0,1,nr)
C                        FIX(ny,1)=.FALSE.
C                      ENDDO !nk
C                    ENDDO !nv
C                  ENDIF
C                ENDIF
C                IF(ISFIXED) THEN
C                  IF(NVHP(nh,np,1,nr).EQ.1.AND.NKH(nh,np,1,nr).EQ.1) 
C     '              THEN
CC*** Don't need to ask a question as it has already been answered
C                    ny=NYNP(1,1,nh,np,0,1,nr)
C                    FIX(ny,1)=.TRUE.
C                  ELSE
C                    DO nv=1,NVHP(nh,np,1,nr) !loop over versions
C                      IF(NVHP(nh,np,1,nr).GT.1) THEN
C                        WRITE(CHAR1,'(I2)') nv
C                        FORMAT='('' For version number '//CHAR1(1:2)
C     '                    //':'')'
C                        CALL GINOUT(IOTYPE,0,IVDU,IFILE,1,0,NOQUES,
C     '                    FILEIP,FORMAT,1,ADATA,ANO,CDATA,CDEFLT,
C     '                    ICHAR,IDATA,IDEFLT,0,IMAX,LDATA,LDEFLT,
C     '                    RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C                      ENDIF
C                      DO nk=1,NKH(nh,np,1,nr) !loop over spatial derivs
C                        ny=NYNP(nk,nv,nh,np,0,1,nr)
C                        ADEFLT(1)='N'
C                        IF(nk.EQ.1) THEN
C                          FORMAT='($,'' Is variable '//CHAR2(1:1)
C     '                      //' of fit variable '//CHAR1(IBEG1:IEND1)
C     '                      //' fixed [N]?: '',A)'
C                        ELSE IF(nk.GT.1) THEN
C                          WRITE(CHAR3,'(I1)') nk
C                          FORMAT='($,'' Is variable '//CHAR2(1:1)
C     '                      //' of fit variable '
C     '                      //CHAR1(IBEG1:IEND1)//' derivative '
C     '                      //CHAR3(1:1)//' fixed [N]?: '',A)'
C                        ENDIF
C                        IF(IOTYPE.EQ.3) THEN
C                          IF(FIX(ny,1)) THEN
C                            ADATA(1)='Y'
C                          ELSE
C                            ADATA(1)='N'
C                          ENDIF
C                        ENDIF
C                        CALL GINOUT(IOTYPE,1,IVDU,IFILE,0,0,NOQUES,
C     '                    FILEIP,FORMAT,1,ADATA,ANO,CDATA,CDEFLT,
C     '                    ICHAR,IDATA,IDEFLT,0,IMAX,LDATA,LDEFLT,
C     '                    RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C                        IF(IOTYPE.NE.3) THEN
C                          IF(ADATA(1).EQ.'Y'.OR.ADATA(1).EQ.'y') THEN
C                            FIX(ny,1)=.TRUE.
C                          ELSE
C                            FIX(ny,1)=.FALSE.
C                          ENDIF
C                        ENDIF
C                        IF(FIX(ny,1).AND.ENTERCOUPLING) THEN
CC*** variable is fixed - find out if it is coupled to another var.
C                          ERROR='>> Not implemented'
C                          GOTO 9999
C                        ENDIF
C                      ENDDO !nk
C                    ENDDO !nv
C                  ENDIF
C                ENDIF
C              ENDDO !nhj
C            ENDDO !njj
C          ENDIF
C        ENDDO !continue


 

      SUBROUTINE IPFIT(IBT,IDO,INP,LD,LN,NBH,NBJ,NEELEM,NJE,NKE,NKH,NKJ,
     '  NONY,NPF,NPNE,NPNODE,NPO,NQE,nr,NRE,NVNE,NVNP,nx,
     '  CONY,SCALE,SE,TYPE,WD,WU,XA,XE,XID,XP,ZD,XO,FIX_FIT,
     '  ERROR,*)

C#### Optimisation parameters for geometry or field fit for region nr.
C###  KTYP8=1 ITYP6(nr) =1 is linear    geometric fitting problem
C###  KTYP8=1 ITYP6(nr) =2 is nonlinear geometric fitting problem
C###  KTYP8=2 JTYP9 =1..3 is fibre/sheet fitting problem (ITYP6(nr)=1)
C###  KTYP8=3 JTYP11>0 is field fitting problem (ITYP6(nr)=1) (& Gauss)
C###  KTYP8=4 is a potential field fitting problem
C###  KTYP8=5 is motion fitting problem with Fourier basis
C###  KTYP8=6 is a fitting problem with optimisation
C###  KTYP8=7 is a fitting problem on gauss points
C###  KTYP12=0 for no smoothing
C###  KTYP12=1 for Sobelov smoothing  (weights defined in WU(i,ne))
C###  KTYP12=2 for Strain energy smoothing
C###  NJG is geometric variable number in linear field fitting.
C###  NJO is field variable number in linear field fitting (eg =NJT+1).
C###  NHO(njj) is Gauss variable to be fitted in Gauss point fitting
C###    for fit variable njj
C**** Note: XO(no) is in real*8 ..corrections needed in nonlinear part

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b21.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cmis00.cmn'
      INCLUDE 'cmiss$reference:fit000.cmn'
      INCLUDE 'cmiss$reference:four00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:head00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:ktyp00.cmn'
      INCLUDE 'cmiss$reference:nsurf00.cmn'
!     Parameter List
      INTEGER IBT(2,NIM,*),IDO(NKM,0:NIM,*),INP(NNM,NIM,*),
     '  LD(*),LN(0:*),NBH(NHM,NCM,*),NBJ(NJM,*),
     '  NEELEM(0:NEM,0:*),NJE(*),NKE(NKM,NNM,NBFM,*),NKH(NHM,NPM,*),
     '  NKJ(NJM,*),NONY(0:NOYM,NYM,*),NPF(12,*),NPNE(NNM,NBFM,*),
     '  NPNODE(0:NPM,0:*),NPO(0:*),NQE(NSM,NBFM,*),nr,NRE(*),
     '  NVNE(NNM,NBFM,NJM,*),NVNP(NJM,NPM,*),nx
      REAL*8 CONY(0:NOYM,NYM,*),SCALE(*),SE(NSM,NBFM,*),WD(NJM,*),
     '  WU(0:10,*),XA(NAM,NJM,*),XE(NSM,*),XID(NJM,*),XP(NKM,NVM,NJM,*),
     '  ZD(NJM,*)
      REAL*8 XO(*)
      CHARACTER ERROR*(*),TYPE*(*)
      LOGICAL FIX_FIT(NYM,*)
!     Local Variables
      INTEGER IBEG,IBEG1,IBEG2,IBEG3,IBEG4,ICHAR,IEND,IEND1,IEND2,
     '  IEND3,IEND4,INFO,ILIM,IPFILE,L,N,n1node,na,NAYT,nb,nd,ne,ni,
     '  nj,njj,NJO,nk,nk2,NKTOT,NLIST,nn,no,noelem,nonode,NOQUES,noy,np,
     '  nrc,nu,ny,NYPJK,NNODE,nv
      REAL*8 PXI,SSCALE,THETA1,THETA2,X(4),XI(3)
      CHARACTER CFROMI*100,CFROMR*100,CHAR1*10,CHAR2*10,CHAR3*10,
     '  CHAR4*10,CIOT*1,CJET*1,CJOT*1,CNJ*4,CNK*4,CNO*4,CNP*4,
     '  CNY*4,CSC*12,CXO*12
      LOGICAL CONTINUE,DATA,FILEIP,FIRST,INLIST,NODE

      CALL ENTERS('IPFIT',*9999)

      nrc=2 !temporary
      nv=1 ! temporary cpb 23/11/94
      
      IPFILE=1 !is input file version number on 24-Jan-1990
      FILEIP=.FALSE.
      NOQUES=0

C CPB 19/4/94 This subroutine needs to be fixed to allow .ipfit files
C to be written (i.e. check on iotypes for ginouts)

C cpb 2/2/94 moved to OPEN_FILE call in DEFIT
C      IPFILE=1 !is input file version number on 24-Jan-1990
C      IF(IOTYPE.EQ.1.OR.IOTYPE.EQ.3) THEN
C        WRITE(UNIT=IFILE,REC=1,FMT='(A,I2)') 'CMISS Version '//CMISS
C     '    //' ipfit File Version ',IPFILE
C        WRITE(UNIT=IFILE,REC=2,FMT='(A)') 'Heading: '//HEADING
C        WRITE(UNIT=IFILE,REC=3,FMT='(1X)')
C      ELSE IF(IOTYPE.EQ.2.OR.IOTYPE.EQ.4) THEN
C        READ(UNIT=IFILE,REC=1,FMT='(39X,I2)') IPFILE
C        IF(DOP) THEN
C          WRITE(OP_STRING,'('' File version number is '',I2)') IPFILE
C      	  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C        ENDIF
C        READ(UNIT=IFILE,REC=2,FMT='(9X,A)') HEADING
C        WRITE(OP_STRING,'('' File heading: '',A)') HEADING
C      	CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C        READ(UNIT=IFILE,REC=3,FMT='(1X)')
C      ENDIF

      IF(KTYP8.NE.6) THEN
      	FORMAT='($,'' Specify whether problem is (1)linear or '//
     '	  '(2)nonlinear [1]: '',I1)'
      	CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IONE,1,2,INFO,
     '	  ERROR,*9999)
      	ITYP6(nr)=IDATA(1)
      ELSE
        ITYP6(nr)=1
      ENDIF

C CPB 19/4/94 Adding NJ_FIT
      IF(KTYP8.EQ.1) THEN         ! geometric fitting problem
        IF(ITYP6(nr).EQ.1) THEN   ! ...is linear
          NJOT=JTYP11             ! is number of field variables defined
C          NJO0=NJT+1
C          NJO1=NJT+NJOT
          NJ_FIT(0)=NJOT
          DO njj=1,NJ_LOC(NJL_FIEL,0)
            NJ_FIT(njj)=NJ_LOC(NJL_FIEL,njj)
          ENDDO
          IF(JTYP11.EQ.1) THEN    ! one geom coord being fitted
            FORMAT='($,'' Specify the coordinate number to be fitted'//
     '        ' [1]: '',I4)'
            CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IONE,1,NJT,
     '        INFO,ERROR,*9999)
            NJG=IDATA(1)
            NJ_FIT(1)=NJ_LOC(NJL_FIEL,NJG)
            NBO=NBJ(NJG,NEELEM(1,1))
            DO nd=1,NDT
              CALL ZX(ITYP10(nr),ZD(1,nd),X)
C              ZD(NJT+1,nd)=X(NJG)
C              WD(NJT+1,nd)=WD(NJG,nd)
              ZD(NJ_FIT(1),nd)=X(NJG)
              WD(NJ_FIT(1),nd)=WD(NJG,nd)
            ENDDO
            DO noelem=1,NEELEM(0,1)
              ne=NEELEM(noelem,1)
C              NBJ(NJT+1,ne)=NBJ(NJG,ne)
              NBJ(NJ_FIT(1),ne)=NBJ(NJG,ne)
            ENDDO
            DO nonode=1,NPNODE(0,1)
              np=NPNODE(nonode,1)
C              NKJ(NJT+1,np)=NKJ(NJG,np)
              NKJ(NJ_FIT(1),np)=NKJ(NJG,np)
C              DO nk=1,NKJ(NJT+1,np)
C                XP(nk,nv,NJT+1,np)=XP(nk,nv,NJG,np)
              DO nk=1,NKJ(NJ_FIT(1),np)
                XP(nk,nv,NJ_FIT(1),np)=XP(nk,nv,NJG,np)
              ENDDO
            ENDDO
C CPB 26/3/93 - Generalising the # of geom coords being fitted from
C 1,njt to njo0,nj01
C CPB 19/4/94 Adding NJ_FIT
          ELSE IF(JTYP11.GT.1) THEN ! > 1 geom coords being fitted
            NBO=NBJ(1,NEELEM(1,1))
            DO nd=1,NDT
              CALL ZX(ITYP10(nr),ZD(1,nd),X)
C              DO nj=NJO0,NJO1
C                ZD(nj,nd)=X(nj-NJT)
C                WD(nj,nd)=WD(nj-NJT,nd)
              DO njj=1,NJ_FIT(0)
                nj=NJ_FIT(njj)
                ZD(nj,nd)=X(NJ_LOC(NJL_GEOM,njj))
                WD(nj,nd)=WD(NJ_LOC(NJL_GEOM,njj),nd)
              ENDDO
            ENDDO
            DO noelem=1,NEELEM(0,1)
              ne=NEELEM(noelem,1)
C              DO nj=NJO0,NJO1
C                NBJ(nj,ne)=NBJ(nj-NJT,ne)
              DO njj=1,NJ_FIT(0)
                NBJ(nj,ne)=NBJ(NJ_LOC(NJL_GEOM,njj),ne)
              ENDDO
            ENDDO
            DO nonode=1,NPNODE(0,1)
              np=NPNODE(nonode,1)
C              DO nj=NJO0,NJO1
C                NKJ(nj,np)=NKJ(nj-NJT,np)
C                DO nk=1,NKJ(nj,np)
C                  XP(nk,nv,nj,np)=XP(nk,nv,NJ-NJT,np)
              DO njj=1,NJ_FIT(0)
                nj=NJ_FIT(njj)
                NKJ(nj,np)=NKJ(NJ_LOC(NJL_GEOM,njj),np)
                DO nk=1,NKJ(nj,np)
                  XP(nk,nv,nj,np)=XP(nk,nv,NJ_LOC(NJL_GEOM,njj),np)
                ENDDO
              ENDDO
            ENDDO
          ELSE
            ERROR='>>No. of field variables nonsensical'
            GO TO 9999
          ENDIF

          IF(ITYP10(nr).GT.1.AND.((JTYP11.EQ.1.AND.NJG.NE.1).OR.JTYP11
     '      .EQ.NJT)) THEN
            !polar coords and not radius
            !Correct angles for 2*PI differences
c cpb 7/4/94 this needs to be generalised for nj_loc
            DO nj=1,NJT
              IF((JTYP11.EQ.1.AND.NJ.EQ.1).OR      !one coord in fit
     '          .(JTYP11.EQ.NJT.AND.NJ.GT.1)) THEN !all coords/pick theta
                DO nd=1,NDT
C                  THETA1=ZD(NJT+NJ,nd) !is theta data value
                  THETA1=ZD(NJ_LOC(NJL_FIEL,nj),nd) !is theta data value
                  ne=LD(nd)
                  IF(ne.GT.0) THEN !ignore this data point
C                    nb=NBJ(NJT+nj,ne)
                    nb=NBJ(NJ_LOC(NJL_FIEL,nj),ne)
                    DO ni=1,NIT(nb)
                      XI(ni)=XID(ni,nd)
                    ENDDO
                    CALL XPXE(NBJ(1,ne),NJE(ne),NKE(1,1,1,ne),
     '                NPF(1,1),NPNE(1,1,ne),NQE(1,1,ne),nr,
     '                NVNE(1,1,1,ne),NVNP,
     '                SE(1,1,ne),XA,XE,XP,ERROR,*9999)
                    IF(nj.EQ.1) THEN       !one coord in fit is NJG
                      nb=NBJ(NJG,ne)
                      THETA2=PXI(IBT(1,1,nb),IDO(1,0,nb),INP(1,1,nb),
     '                  nb,1,XI,XE(1,NJG)) !is interpolated theta
                    ELSE IF(nj.GT.1) THEN !all coords
                      nb=NBJ(nj,ne)
                      THETA2=PXI(IBT(1,1,nb),IDO(1,0,nb),INP(1,1,nb),
     '                  nb,1,XI,XE(1,nj)) !is interpolated theta
                    ENDIF
                    IF(DOP) THEN
                      WRITE(OP_STRING,'('' Original data angle ='','
     '                  //'E12.3,'' deg. Interpolated angle ='',E12.3,'
     '                  //''' deg.'','' Difference = '',E12.3,'
     '                  //''' deg.'')')
     '                  THETA1*180.0D0/PI,THETA2*180.0D0/PI,
     '                  (THETA1-THETA2)*180.0D0/PI
      		      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                    IF(DABS(THETA1-THETA2).GT.PI) THEN !needs 2n*pi correction
                      IF(DOP) THEN
                        WRITE(OP_STRING,'('' Correction needed..'')')
      			CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      ENDIF
                      N=-3
                      CONTINUE=.TRUE.
                      DO WHILE(CONTINUE)
                        N=N+1
                        IF(DABS(THETA1-THETA2+2.0D0*DBLE(N)*PI).LT.PI)
     '                    THEN
                          THETA1=THETA1+2.0D0*DBLE(N)*PI
                          CONTINUE=.FALSE.
                        ELSE IF(N.GT.5) THEN
                          WRITE(OP_STRING,
     '                      '('' Warning: Data correction failed for'
     '                      //' nd='',I6)') ND
      			  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                          CONTINUE=.FALSE.
                        ENDIF
                      ENDDO
C                      ZD(NJT+nj,nd)=THETA1 !is new value
                      ZD(NJ_LOC(NJL_FIEL,nj),nd)=THETA1 !is new value
                      IF(DOP) THEN
                        WRITE(OP_STRING,'('' New ZD('',I1,'') for nd'
     '                    //' = '',I3,'' is '',E12.3)') 
     '                    NJT+nj,nd,THETA1
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      ENDIF
                    ENDIF
                  ENDIF
                ENDDO
              ENDIF
            ENDDO
          ENDIF

        ELSE IF(ITYP6(nr).EQ.2) THEN ! nonlinear geometry fit
          NJOT=NJT
C CPB 19/4/94 Adding NJ_FIT
C          NJO0=1
C          NJO1=NJT
          NJ_FIT(0)=NJT
          DO nj=1,NJT
            NJ_FIT(nj)=nj
          ENDDO
          NBO=NBJ(1,1)
        ENDIF
C       NIOT=NJT-1         !!!!!Check whether this is needed
        NIOT=NIT(NBO)      !try this AAY 17/4/90

      ELSE IF(KTYP8.EQ.2) THEN       ! fibre/sheet fitting problem
        IF(JTYP9.LE.2) THEN          !...is fibre type
          NJOT=1
C CPB 19/4/94 Adding NJ_FIT
C          NJO0=NJT+JTYP9
C          NJO1=NJO0
          NJ_FIT(0)=1
          NJ_FIT(1)=NJ_LOC(NJL_FIBR,1)
        ELSE IF(JTYP9.EQ.3) THEN     !...is sheet type
          NJOT=1
C CPB 19/4/94 Adding NJ_FIT
C          NJO0=NJT+3
C          NJO1=NJO0
          NJ_FIT(0)=1
          NJ_FIT(1)=NJ_LOC(NJL_FIBR,3)
        ENDIF      
C        NBO=NBJ(NJO0,NEELEM(1,nr))
        NBO=NBJ(NJ_FIT(1),NEELEM(1,nr))
        NIOT=NIT(NBO)
        IF(TYPE(1:5).EQ.'GAUSS') THEN
          FORMAT='($,'' Specify the Gauss point variable to be fitted'
     '      //' [1]: '',I4)'
          CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IONE,1,NJM,INFO,
     '      ERROR,*9999)
          NHO(1)=IDATA(1)
          L=0
          DO noelem=1,NEELEM(0,1)
            ne=NEELEM(noelem,1)
            L=L+1
            LN(L)=ne
          ENDDO
          LN(0)=L
        ENDIF

      ELSE IF(KTYP8.EQ.3.OR.KTYP8.EQ.7) THEN ! field or Gauss fitting problem
        FORMAT='($,'' Specify the number of field variables'//
     '    ' to be fitted [1]: '',I4)'
        CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IONE,1,JTYP11,
     '    INFO,ERROR,*9999)
        NJOT=IDATA(1)
C CPB 19/4/94 Adding NJ_FIT
C        NJO0=NJT+1
C        NJO1=NJT+NJOT
        NJ_FIT(0)=NJOT
        DO nj=1,NJ_FIT(0)
          IDEFLT(1)=nj
          CHAR1=CFROMI(IDEFLT(1),'(I1)')
          FORMAT='($,'' Specify the field variable # to be fitted'
     '      //' for fit variable '//CHAR1(1:1)//' ['//CHAR1(1:1)
     '      //']: '',I4)'
          CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IDEFLT(1),1,
     '      JTYP11,INFO,ERROR,*9999)
          NJ_FIT(nj)=NJ_LOC(NJL_FIEL,IDATA(1))
        ENDDO
C CPB 19/4/94 Not sure why this statement is in
C        JTYP11=NJOT
C        NBO=NBJ(NJO0,NEELEM(1,nr))
        NBO=NBJ(NJ_FIT(1),NEELEM(1,nr))
        NIOT=NIT(NBO)
        IF(TYPE(1:5).EQ.'GAUSS') THEN
          DO nj=1,NJ_FIT(0)
            CHAR1=CFROMI(nj,'(I1)')
            FORMAT='($,'' Specify the Gauss point variable to be '
     '        //'fitted for fit variable '//CHAR1(1:1)//' ['
     '        //CHAR1(1:1)//']: '',I2)'
            CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,nj,1,NJM,INFO,
     '        ERROR,*9999)
            NHO(nj)=IDATA(1)
          ENDDO
          L=0
          DO noelem=1,NEELEM(0,1)
            ne=NEELEM(noelem,1)
            L=L+1
            LN(L)=ne
          ENDDO
          LN(0)=L
        ENDIF

      ELSE IF(KTYP8.EQ.4) THEN ! potential field fitting problem
        FORMAT='($,'' Specify the field variable to fit'//
     '    ' [1]: '',I4)'
        CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IONE,1,JTYP11,
     '    INFO,ERROR,*9999)
        NJF=IDATA(1)
        NJOT=1
C CPB 19/4/94 Adding NJ_FIT
C        NJO0=NJT+1
C        NJO1=NJT+NJOT
        NJ_FIT(0)=NJOT
        NJ_FIT(1)=NJ_LOC(NJL_FIEL,NJF)
C CPB 19/4/94 Not sure why this statement is in
C        JTYP11=NJOT
C        NBO=NBJ(NJO0,NEELEM(1,nr))
        NBO=NBJ(NJ_FIT(1),NEELEM(1,nr))
        NIOT=NIT(NBO)
        IF(TYPE(1:5).EQ.'GAUSS') THEN
          FORMAT='($,'' Specify the Gauss point variable to be fitted'
     '      //' [1]: '',I4)'
          CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IONE,1,NJM,INFO,
     '      ERROR,*9999)
          NHO(1)=IDATA(1)
          L=0
          DO noelem=1,NEELEM(0,1)
            ne=NEELEM(noelem,1)
            L=L+1
            LN(L)=ne
          ENDDO
          LN(0)=L
        ENDIF

      ELSE IF(KTYP8.EQ.5) THEN ! motion fitting with Fourier basis
        CALL ASSERT(ITYP6(nr).EQ.1,'>>Linear fits only',ERROR,*9999)
        WRITE(OP_STRING,'(''>>Check code for correctness'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        NJOT=1
C CPB 19/4/94 Adding NJ_FIT
C        NJO0=NJT+1
C        NJO1=NJT+NJOT
        NJ_FIT(0)=1
        NJ_FIT(1)=NJ_LOC(NJL_FIEL,1)
        FORMAT='($,'' Specify the coordinate number to be fitted'//
     '    ' [1]: '',I4)'
        CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IONE,1,NJT,INFO,
     '    ERROR,*9999)
        NJG=IDATA(1)
C       basis type number assumed to be the same as the
C       (static) geometric variable
        NBO=NBJ(NJG,NEELEM(1,1))
        DO nd=1,NDT
C CPB 19/4/94 Adding NJ_FIT
C          DO nj=NJO0,NJO1
          DO njj=1,NJ_FIT(0)
            nj=NJ_FIT(njj)
            ZD(nj,nd)=ZD(NJG,nd)
            WD(nj,nd)=WD(NJG,nd)
          ENDDO
        ENDDO
        DO noelem=1,NEELEM(0,1)
          ne=NEELEM(noelem,1)
C          NBJ(NJO0,ne)=NBH(NJG,1,ne)
          NBJ(NJ_FIT(1),ne)=NBH(NJG,1,ne)
          CALL ASSERT(NBH(NJG,1,ne).GT.0,'>>Define motion first',ERROR,
     '      *9999)
        ENDDO
        DO nonode=1,NPNODE(0,1)
          np=NPNODE(nonode,1)
C          NKJ(NJO0,np)=NKH(NJG,np,1)
          NKJ(NJ_FIT(1),np)=NKH(NJG,np,1)
          CALL ASSERT(NKH(NJG,np,1).GT.0,'>>Define motion first',ERROR,
     '      *9999)
        ENDDO
        NIOT=NIT(NBO)
      ELSE IF(KTYP8.EQ.6) THEN   ! data fitting with optimisation
      	NJOT=NJT
C CPB 19/4/94 Adding NJ_FIT
C      	NJO0=1
C      	NJO1=NJT
        NJ_FIT(0)=NJT
        DO njj=1,NJT
          NJ_FIT(njj)=njj
        ENDDO
      	NBO=NBJ(1,1)
        NIOT=NIT(NBO)  
      ENDIF

      CJOT=CFROMI(NJOT,'(I1)')
      CIOT=CFROMI(NIOT,'(I1)')
      CJET=CFROMI(NJT,'(I1)')

      IF(KTYP8.EQ.6) THEN  ! Fitting by optimisaton
      	FORMAT='('' Enter smoothing type [0]:'''
     '	  //'/''   (0) None'''
     '	  //'/''   (1) Sobelov'''
     '	  //'/$,4X,I1)'
        ILIM=1
      ELSE
      	FORMAT='('' Enter smoothing type [0]:'''
     '	  //'/''   (0) None'''
     '	  //'/''   (1) Sobelov'''
     '	  //'/''   (2) Strain energy'''
     '	  //'/$,4X,I1)'
        ILIM=2
      ENDIF
      IF(IOTYPE.EQ.3) THEN
        IDATA(1)=KTYP12
      ENDIF
      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,ILIM,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) KTYP12=IDATA(1)

 2000 IF(KTYP12.EQ.0.AND.LN(0).EQ.0) THEN
        L=0
        DO noelem=1,NEELEM(0,1)
          ne=NEELEM(noelem,1)
          DATA=.FALSE.
          DO nd=1,NDT
            IF(LD(nd).EQ.ne) THEN
              DATA=.TRUE.
              GO TO 2002
            ENDIF
          ENDDO
 2002     IF(DATA) THEN
            L=L+1
            LN(L)=ne
          ENDIF
        ENDDO
        LN(0)=L
      ELSE IF(LN(0).EQ.0)THEN
        L=0
        DO noelem=1,NEELEM(0,1)
          ne=NEELEM(noelem,1)
          IF(NIT(NBJ(1,ne)).EQ.NIOT)THEN
            L=L+1
            LN(L)=ne
          ENDIF
        ENDDO
        LN(0)=L
      ENDIF
      IF(DOP) THEN
        CHAR1=CFROMI(LN(0),'(I10)')
        CALL TRIM(CHAR1,IBEG1,IEND1)
        FORMAT='(/'' The '//CHAR1(IBEG1:IEND1)//
     '    ' elements involved in the optimisation are:'')'
        WRITE(OP_STRING,FORMAT)
      	CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(20I4)') (LN(L),L=1,LN(0))
      	CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

C *** Calculate list of nodes in fit NPO(n),n=1,NPO(0)
      NPO(0)=0
      DO nonode=1,NPNODE(0,1)
        np=NPNODE(nonode,1)
        NODE=.FALSE.
        DO L=1,LN(0)
          ne=LN(L)
          DO nn=1,NNT(NBO)
            IF(np.EQ.NPNE(nn,NBO,ne)) NODE=.TRUE.
          ENDDO
        ENDDO
        IF(NODE) THEN
          NPO(0)=NPO(0)+1
          NPO(NPO(0))=np
        ENDIF
      ENDDO

C *** Check that NYM is big enough
      NYT(nrc,1,nx)=0
      DO nonode=1,NPNODE(0,nr)
        np=NPNODE(nonode,nr)
C cpb 5/9/93 use all nj's to calculate nyt for optimisation problems
        IF(KTYP8.EQ.6) THEN
          DO nj=1,NJT
            DO nk=1,NKJ(nj,np)
              NYT(nrc,1,nx)=NYT(nrc,1,nx)+1
            ENDDO
          ENDDO
        ELSE
C CPB 14/2/93 Most fits just fit in one variable so only calculate NYT 
C for one NJ
C          nj=NJO0 
          nj=NJ_FIT(1) 
          DO nk=1,NKJ(nj,np)
            NYT(nrc,1,nx)=NYT(nrc,1,nx)+1
          ENDDO
        ENDIF
      ENDDO
      CALL ASSERT(NYT(nrc,1,nx).LE.NYM,'>>NYM too small',ERROR,*9999)
C     next line commented due to space problems AAY
C     IF(KTYP8.EQ.5) CALL ASSERT(NYT(1,nx)*NJG.LE.NYM,'>>NYM too small',
C    '  ERROR,*9999)


C cpb 12/10/94 Commenting this out as it does not appear to be used
C and corrupts NXI

CC CPB 7/10/94 this call should be replaced by NENXI
C*** LNLXI to archive
C      CALL LNLXI(INP,LN,NPNE,NXI)
C
C      IF(DOP) THEN
C        FORMAT='('' '')'
C        WRITE(OP_STRING,FORMAT)
C      	CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        DO L=1,LN(0)
C          CDATA(1)=' '
C          DO ni=-NIOT,NIOT
C            CALL TRIM(CDATA(1),IBEG1,IEND1)
C            CDATA(1)=CDATA(1)(1:IEND1)//CFROMI(LN(NXI(ni,L)),'(I4)')
C          ENDDO
C          CALL TRIM(CDATA(1),IBEG1,IEND1)
C          CHAR2=CFROMI(LN(L),'(I10)')
C          CALL TRIM(CHAR2,IBEG2,IEND2)
C          FORMAT='('' For element number '//CHAR2(IBEG2:IEND2)//
C     '      ' the neighbouring elements are: '//
C     '      CDATA(1)(1:IEND1)//''')'
C          WRITE(OP_STRING,FORMAT)
C      	  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        ENDDO
C      ENDIF

C *** Specify smoothing constraints on each element
      IF(KTYP12.EQ.1) THEN !Sobelov smoothing
        DO nu=1,9
          RDEFLT(nu)=0.0D0
        ENDDO
        DO L=1,LN(0)
          ne=LN(L)
          WU(0,ne)=1.0D0 ! the scaling factor for the Sobelov weights
          CHAR4=CFROMI(ne,'(I4)')
          CALL TRIM(CHAR4,IBEG,IEND)
          IF(NIT(NBJ(1,NEELEM(1,1))).EQ.1) THEN
            NKTOT=2
            FORMAT='($,'' For element '//CHAR4(IBEG:IEND)//':''/,'
     '        //''' The 2 weights on derivs wrt Xi_1/_11 '//
     '        'are [prev]:'',2E9.2)'
          ELSE IF(NIT(NBJ(1,NEELEM(1,1))).EQ.2) THEN
            NKTOT=5
            FORMAT='($,'' For element '//CHAR4(IBEG:IEND)//':''/,'
     '        //''' The 5 weights on derivs wrt Xi_1/_11/_2/_22/'
     '        //'_12 are [prev]:'',5E9.2)'
          ELSE IF(NIT(NBJ(1,NEELEM(1,1))).EQ.3) THEN
            NKTOT=9
            FORMAT='($,'' For element '//CHAR4(IBEG:IEND)//':''/,'
     '        //''' The 9 weights on derivs wrt Xi_1/_11/_2/_22/'
     '        //'_12/_3/_33/_23/_31 are [prev]:'',/,9E9.2)'
          ENDIF
          CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,NKTOT,RDATA,
     '      RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          DO nu=1,NKTOT
            WU(nu,ne)=RDATA(nu)
            RDEFLT(nu)=RDATA(nu)
          ENDDO
          IF(DOP) THEN
            WRITE(OP_STRING,'('' WU(0..9,ne)='',/,9E9.2)') 
     '        (WU(nu,ne),NU=0,9)
      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO

      ELSE IF(KTYP12.EQ.2) THEN !Strain energy smoothing
        DO L=1,LN(0)
          ne=LN(L)
          CHAR4=CFROMI(ne,'(I4)')
          CALL TRIM(CHAR4,IBEG,IEND)
          FORMAT='($,'' Enter weight for element '//CHAR4(IBEG:IEND)
     '      //' [0]: '',E12.3)'
          IF(IOTYPE.EQ.3) RDATA(1)=WU(1,ne)
          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            WU(0,ne)=1.0D0    !total weight
            WU(1,ne)=RDATA(1) !Xi(1) weight
            WU(3,ne)=RDATA(1) !Xi(2) weight
            WU(6,ne)=RDATA(1) !Xi(3) weight
            WU(2,ne)=0.0D0
            WU(4,ne)=0.0D0
            WU(5,ne)=0.0D0
            WU(7,ne)=0.0D0
            WU(8,ne)=0.0D0
            WU(9,ne)=0.0D0
          ENDIF
        ENDDO
      ENDIF

C *** Enter constraints on fitting
      IF(KTYP8.NE.5) THEN !geom/field/fibre/sheet/optimisation/potential
        DO ny=1,NYT(nrc,1,nx)
          FIX_FIT(ny,1)=.FALSE.
          FIX_FIT(ny,2)=.FALSE.
          NONY(0,ny,2)=1
        ENDDO
        CONTINUE=.TRUE.
        NNODE=0
        DO WHILE(CONTINUE)
          FORMAT='(/$,'' Enter node to fix or to specify coupling'
     '      //' [EXIT]: '',I3)'
          IF(IOTYPE.EQ.3) THEN
            IF(np.GT.NPNODE(0,nr)) THEN
              CONTINUE=.FALSE.
              IDATA(1)=0
            ELSE
              NNODE=NNODE+1
              np=NPNODE(NNODE,nr)
              IDATA(1)=np
            ENDIF
          ENDIF
          CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IZERO,1,NPM,
     '      INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IF(IDATA(1).EQ.0) THEN
              CONTINUE=.FALSE.
            ELSE IF(IDATA(1).NE.0) THEN
              CONTINUE=.TRUE.
              np=IDATA(1)
            ENDIF
          ENDIF
          nonode=0
          n1node=0
          DO WHILE(n1node.EQ.0.AND.nonode.LT.NPNODE(0,nr))
            nonode=nonode+1
            IF(NPNODE(nonode,nr).EQ.NP) THEN
              n1node=nonode
            ENDIF
          ENDDO
          ny=0
          DO nonode=1,n1node-1
C CPB 19/4/94 Adding NJ_FIT
C            DO nj=NJO0,NJO1
            DO njj=1,NJ_FIT(0)
              nj=NJ_FIT(njj)
              DO nk=1,NKJ(nj,NPNODE(nonode,nr))
                ny=ny+1
              ENDDO
            ENDDO
          ENDDO
          IF(CONTINUE) THEN
C CPB 19/4/94 Adding NJ_FIT
C            DO nj=NJO0,NJO1
c             nj=NJO0
            DO njj=1,NJ_FIT(0)
              nj=NJ_FIT(njj)
C              IF(NJO0.EQ.NJO1) THEN !only one fitting variable
              IF(NJ_FIT(0).EQ.1) THEN !only one fitting variable
                ADATA(2)='Y'
              ELSE       !more than one fitting variable
                IF(KTYP8.EQ.6) THEN
      		  CHAR2=CFROMI(nj,'(I10)')
                ELSE
C      		  CHAR2=CFROMI(nj-NJO0+1,'(I10)')
      		  CHAR2=CFROMI(njj,'(I10)')
                ENDIF
                CALL TRIM(CHAR2,IBEG2,IEND2)
                FORMAT='($,'' Are there any optimisation variables'//
     '            ' coupled to mesh variable '//CHAR2(IBEG2:IEND2)//
     '            ' [N]? '',A)'
                CALL AINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,ADATA(2),ANO,
     '            INFO,ERROR,*9999)
              ENDIF
              DO nk=1,NKJ(nj,np) !loop over spatial derivatives
                ny=ny+1
                IF(ADATA(2).EQ.'Y') THEN
                  CHAR3=CFROMI(nk-1,'(I10)')
                  CALL TRIM(CHAR3,IBEG3,IEND3)
                  FORMAT='($,'' Are there any optimisation variables'//
     '              ' coupled to derivative '//CHAR3(IBEG3:IEND3)//
     '              ' [N]? '',A)'
                  CALL AINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,ADATA(3),ANO,
     '              INFO,ERROR,*9999)
                ELSE
                  ADATA(3)='N'
                ENDIF
                IF(ADATA(3).EQ.'Y') THEN
                  FIX_FIT(ny,2)=.TRUE.
      		  IF(KTYP8.EQ.6) THEN
                    NONY(0,ny,2)=1
                    FIX_FIT(ny,1)=.FALSE.
                    FIX_FIT(ny,2)=.FALSE.
      		  ELSE
C      		    CHAR2=CFROMI(nj-NJO0+1,'(I10)')
      		    CHAR2=CFROMI(njj,'(I10)')
                    CALL TRIM(CHAR2,IBEG2,IEND2)
                    CHAR3=CFROMI(nk-1,'(I10)')
                    CALL TRIM(CHAR3,IBEG3,IEND3)
                    FORMAT='($,'' Mesh variable '//CHAR2(IBEG2:IEND2)//
     '                ' derivative '//CHAR3(IBEG3:IEND3)//
     '                ' Number of optimisation variables is [1]:'',I4)'
                    CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IONE,1,
     '                NYT(2,1,nx),INFO,ERROR,*9999)
                    NONY(0,ny,2)=IDATA(1)
                    IF(NONY(0,ny,2).EQ.0) THEN
                      FIX_FIT(ny,1)=.TRUE.
                    ELSE
                      FIX_FIT(ny,1)=.FALSE.
                      DO noy=1,NONY(0,ny,2)
                        IDEFLT(noy)=-1
                      ENDDO
                      CHAR4=CFROMI(NONY(0,ny,2),'(I10)')
                      CALL TRIM(CHAR4,IBEG4,IEND4)
                      FORMAT='($,'' Specify the '//CHAR4(IBEG4:IEND4)//
     '                  ' optimisation degrees of freedom'//
     '                  ' [1 to 1]: '',I4,:,/(20I4))'
                      CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,NONY(0,ny,2),
     '                  IDATA,IDEFLT,0,NOM,INFO,ERROR,*9999)
                      DO noy=1,NONY(0,ny,2)
                        IF(IDATA(noy).GT.0) THEN
                          NONY(noy,ny,2)=IDATA(noy)
                        ELSE
                          FIX_FIT(ny,2)=.FALSE.
                        ENDIF
                      ENDDO
                      FORMAT='($,'' Specify the '//CHAR4(IBEG4:IEND4)//
     '                  ' associated coupling coefficients'//
     '                  ' [1.0]: '',E12.6,:, /1X,(5(E12.6,4X)))'
                      CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,NONY(0,ny,2),
     '                  RDATA,RONE,-RMAX,RMAX,INFO,ERROR,*9999)
                      DO noy=1,NONY(0,ny,2)
                        CONY(noy,ny,2)=RDATA(noy)
                      ENDDO
      		    ENDIF
                  ENDIF
                  IF(DOP) THEN
                    CNP=CFROMI(np,'(I4)')
      		    IF(KTYP8.EQ.6) THEN
      		      CNJ=CFROMI(nj,'(I10)')
      		    ELSE
C      		      CNJ=CFROMI(nj-NJO0+1,'(I10)')
      		      CNJ=CFROMI(njj,'(I10)')
      		    ENDIF
                    CNK=CFROMI(nk-1,'(I4)')
                    CNY=CFROMI(ny,'(I4)')
                    CNO=CFROMI(NONY(0,ny,2),'(I4)')
                    FORMAT='('' Node'//CNP//
     '                ' Variable'//CNJ//
     '                ' Derivative'//CNK//
     '                ' Mesh freedom'//CNY//
     '                ' Number of optimisation freedoms'//CNO//''')'
                    WRITE(OP_STRING,FORMAT)
      		    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    WRITE(OP_STRING,'('' Freedom number  '',8I12)')
     '                (NONY(noy,ny,2),noy=1,NONY(0,ny,2))
      		    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    WRITE(OP_STRING,'('' Freedom coeff   '',8E12.5)')
     '                (CONY(noy,ny,2),noy=1,NONY(0,ny,2))
      		    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                ELSE IF(ADATA(3).EQ.'N') THEN
                  NONY(0,ny,2)=0
                  FIX_FIT(ny,1)=.TRUE.
                  IF(DOP) THEN
                    CNP=CFROMI(np,'(I4)')
      		    IF(KTYP8.EQ.6) THEN
      		      CNJ=CFROMI(nj,'(I10)')
      		    ELSE
C      		      CNJ=CFROMI(nj-NJO0+1,'(I10)')
      		      CNJ=CFROMI(njj,'(I10)')
      		    ENDIF
                    CNK=CFROMI(nk-1,'(I4)')
                    CNY=CFROMI(ny,'(I4)')
                    CNO=CFROMI(NONY(0,ny,2),'(I4)')
                    FORMAT='('' Node'//CNP//
     '                ' Variable'//CNJ//
     '                ' Derivative'//CNK//
     '                ' Mesh freedom'//CNY//
     '                ' Number of optimisation freedoms'//CNO//''')'
                    WRITE(OP_STRING,FORMAT)
      		    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDIF
        ENDDO

      ELSE IF(KTYP8.EQ.5) THEN !motion fitting with Fourier basis
        DO ny=1,NYT(nrc,1,nx)
          FIX_FIT(ny,1)=.FALSE.
          FIX_FIT(ny,2)=.FALSE.
          NONY(0,ny,2)=1
        ENDDO
        CONTINUE=.TRUE.
        DO WHILE(CONTINUE)
          FORMAT='(/$,'' Enter node to fix or to specify coupling'
     '      //' [EXIT]: '',I3)'
          IF(IOTYPE.EQ.3) THEN
            CONTINUE=.FALSE.
            IDATA(1)=0
          ENDIF
          CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IZERO,1,NPM,
     '      INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IF(IDATA(1).EQ.0) THEN
              CONTINUE=.FALSE.
            ELSE IF(IDATA(1).NE.0) THEN
              CONTINUE=.TRUE.
              np=IDATA(1)
              CONTINUE=INLIST(np,NPO(1),NPO(0),NLIST)
            ENDIF
          ENDIF
          IF(CONTINUE) THEN
C            DO nj=NJO0,NJO1  !loop over geometric variables (coords)
            DO njj=1,NJ_FIT(0)  !loop over geometric variables (coords)
              nj=NJ_FIT(njj)
C              IF(nj.EQ.NJO0) THEN
              IF(njj.EQ.1) THEN
                ADATA(2)='Y'
              ELSE
C                CHAR2=CFROMI(nj-NJO0+1,'(I10)')
                CHAR2=CFROMI(njj,'(I10)')
                CALL TRIM(CHAR2,IBEG2,IEND2)
                FORMAT='($,'' Are there any optimisation variables'//
     '            ' coupled to mesh variable '//CHAR2(IBEG2:IEND2)//
     '            ' [N]? '',A)'
                CALL AINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,ADATA(2),ANO,
     '            INFO,ERROR,*9999)
              ENDIF
              DO nk=1,NKJ(NJG,np) !loop over spatial derivatives
                IF(ADATA(2).EQ.'Y') THEN
                  CHAR3=CFROMI(nk-1,'(I10)')
                  CALL TRIM(CHAR3,IBEG3,IEND3)
                  FORMAT='($,'' Are there any optimisation variables'//
     '              ' coupled to spatial derivative '
     '              //CHAR3(IBEG3:IEND3)//' [N]? '',A)'
                  CALL AINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,ADATA(3),ANO,
     '              INFO,ERROR,*9999)
                ELSE
                  ADATA(3)='N'
                ENDIF
                DO nk2=NK,NKJ(nj,np)-NKJ(NJG,np)+nk,NKJ(NJG,np)
C                 loop over ny corresp. to np,nk2
C                 NJO not used in function NYPJK    GBS 24/1/92
                  ny=NYPJK(NJO,nk2,NKJ,np,NPO)
                  IF(ADATA(3).EQ.'Y') THEN
                    FIX_FIT(ny,2)=.TRUE.
C                    CHAR2=CFROMI(nj-NJO0+1,'(I10)')
                    CHAR2=CFROMI(njj,'(I10)')
                    CALL TRIM(CHAR2,IBEG2,IEND2)
                    CHAR3=CFROMI(nk-1,'(I10)')
                    CALL TRIM(CHAR3,IBEG3,IEND3)
                    CHAR4=CFROMI(na,'(I10)')
                    CALL TRIM(CHAR4,IBEG4,IEND4)
                    FORMAT='($,'' Mesh variable '//CHAR2(IBEG2:IEND2)//
     '                ' derivative '//CHAR3(IBEG3:IEND3)//
     '                ' fourier coeff '//CHAR4(IBEG4:IEND4)//
     '                ' Number of optimisation variables is [1]:'',I4)'
                    CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IONE,
     '                1,NAYT,INFO,ERROR,*9999)
                    NONY(0,ny,2)=IDATA(1)
                    IF(NONY(0,ny,2).EQ.0) THEN
                      FIX_FIT(ny,1)=.TRUE.
                    ELSE
                      FIX_FIT(ny,1)=.FALSE.
                      DO noy=1,NONY(0,ny,2)
                        IDEFLT(noy)=-1
                      ENDDO
                      CHAR4=CFROMI(NONY(0,ny,2),'(I10)')
                      CALL TRIM(CHAR4,IBEG4,IEND4)
                      FORMAT='($,'' Specify the '//CHAR4(IBEG4:IEND4)//
     '                  ' optimisation degrees of freedom'//
     '                  ' [1 to 1]: '',I4,:,/(20I4))'
                      CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,NONY(0,ny,2),
     '                  IDATA,IDEFLT,0,NOM,INFO,ERROR,*9999)
                      DO noy=1,NONY(0,ny,2)
                        IF(IDATA(noy).GT.0) THEN
                          NONY(noy,ny,2)=IDATA(noy)
                        else
                          fix_fit(ny,2)=.false.
                        ENDIF
                      ENDDO
                      FORMAT='($,'' Specify the '//CHAR4(IBEG4:IEND4)//
     '                  ' associated coupling coefficients'//
     '                  ' [1.0]: '',E12.6,:, /1X,(5(E12.6,4X)))'
                      CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,NONY(0,ny,2),
     '                  RDATA,RONE,-RMAX,RMAX,INFO,ERROR,*9999)
                      DO noy=1,NONY(0,ny,2)
                        CONY(noy,ny,2)=RDATA(noy)
                      ENDDO
                    ENDIF
                    IF(DOP) THEN
                      CNP=CFROMI(np,'(I4)')
C                      CNJ=CFROMI(nj-NJO0+1,'(I4)')
                      CNJ=CFROMI(njj,'(I4)')
                      CNK=CFROMI(nk-1,'(I4)')
                      CNY=CFROMI(ny,'(I4)')
                      CNO=CFROMI(NONY(0,ny,2),'(I4)')
                      FORMAT='('' Node'//CNP
     '                  //' Variable'//CNJ
     '                  //' Derivative'//CNK
     '                  //' Mesh freedom'//CNY
     '                  //' Number of optimisation freedoms'//CNO//''')'
                      WRITE(OP_STRING,FORMAT)
      		      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      WRITE(OP_STRING,'('' Freedom number     '','
     '                  //'8I12)') (NONY(noy,ny,2),noy=1,NONY(0,ny,2))
      		      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      WRITE(OP_STRING,'('' Freedom coefficient'','
     '                  //'8E12.5)') (CONY(noy,ny,2),noy=1,NONY(0,ny,2))
      		      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                  ELSE IF(ADATA(3).EQ.'N') THEN
                    NONY(0,ny,2)=0
                    FIX_FIT(ny,1)=.TRUE.
                    IF(DOP) THEN
                      CNP=CFROMI(np,'(I4)')
C                      CNJ=CFROMI(nj-NJO0+1,'(I4)')
                      CNJ=CFROMI(njj,'(I4)')
                      CNK=CFROMI(nk-1,'(I4)')
                      CNY=CFROMI(ny,'(I4)')
                      CNO=CFROMI(NONY(0,ny,2),'(I4)')
                      FORMAT='('' Node'//CNP//
     '                  ' Variable'//CNJ//
     '                  ' Derivative'//CNK//
     '                  ' Mesh freedom'//CNY//
     '                  ' Number of optimisation freedoms'//CNO//''')'
                      WRITE(OP_STRING,FORMAT)
      		      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
        NOT(nrc,nr,nx)=0
        DO ny=1,NYT(nrc,1,nx)
          IF(.NOT.FIX_FIT(ny,1)) THEN !mesh dof ny is in fit
            IF(FIX_FIT(ny,2)) THEN 
              !mesh dof ny has NONY(0,ny,2) specified opt vars.
            ELSE IF(.NOT.FIX_FIT(ny,2)) THEN 
              !mesh dof ny has 1 unspecified opt var.
              NOT(nrc,nr,nx)=NOT(nrc,nr,nx)+1
              NONY(1,ny,2)=NOT(nrc,nr,nx)
              CONY(1,ny,2)=1.0D0
            ENDIF
          ENDIF
        ENDDO
      ENDIF

 5000 IF(ITYP6(nr).EQ.2) THEN !nonlinear fitting
C       IO5=25
C       IF(EXIST(FOPTM//' IOHESS '//OPMODE,ERROR)) THEN
C         CALL OPENF(IO5,'DISK',FOPTM//' IOHESS '//OPMODE,'SEQUENTIAL',
C    '      'FORMATTED',132,ERROR,*9999)
C         NOT1=NOT(nr,nx)
C         CALL IOHESS('READ',IO5,NOT1,XO,WK1,ERROR,*9999)
C         IF(NOT.EQ.NOT1) THEN
C           FORMAT='(/'' Do wish to read initial values for the'//
C    '        ' Hessian matrix (Default is NO)? '',A)'
C           CALL AINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,ADATA,ANO,
C    '        INFO,ERROR,*9999)
C           IF(ADATA(1).EQ.'Y') THEN
C             MODE=3
C           ELSE
C             MODE=1
C           ENDIF
C         ENDIF
C       ELSE
C         CALL OPENF(IO5,'DISK',FOPTM//' IOHESS '//OPMODE,'SEQUENTIAL',
C    '      'FORMATTED',132,ERROR,*9999)
C         MODE=1
C       ENDIF
C       IF(MODE.EQ.1) THEN
C****     CALL XPXO(IWORK,NJP,NKJ,NOT(nr,nx),NONY,NPO,CONY,XA,XO,XP,WK1)
C       ENDIF
        FORMAT='(/$,'' Do you wish to set initial values for'//
     '    ' optimisation variables [N]? '',A)'
        CALL AINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,ADATA,ANO,
     '    INFO,ERROR,*9999)
        IF(ADATA(1).EQ.'Y') THEN
          DO no=1,NOT(nrc,nr,nx)
            CHAR1=CFROMI(no,'(I10)')
            CALL TRIM(CHAR1,IBEG1,IEND1)
            CHAR2=CFROMR(XO(no),'(E12.6)')
            CALL TRIM(CHAR2,IBEG2,IEND2)
            FORMAT='($,'' Specify optimisation variable '//
     '        CHAR1(IBEG1:IEND1)//' ['//CHAR2(IBEG2:IEND2)//
     '        ' ]: '',E12.6)'
            CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,XO(no),
     '        -RMAX,RMAX,INFO,ERROR,*9999)
            XO(no)=RDATA(1)
          ENDDO
        ENDIF
        SSCALE=0.0D0
        DO no=1,NOT(nrc,nr,nx)
          SSCALE=SSCALE+DABS(XO(no))
        ENDDO
        IF(SSCALE.LE.RMIN) THEN
          SSCALE=1.0D0
        ELSE
          SSCALE=SSCALE/(NOT(nrc,nr,nx)*10.0D0)
        ENDIF
        FORMAT='($,'' Do you wish to specify the scale associated'//
     '    ' with each optimisation variable [N]? '',A)'
        CALL AINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,ADATA,ANO,
     '    INFO,ERROR,*9999)
        DO no=1,NOT(nrc,nr,nx)
          IF(ADATA(1).EQ.'Y') THEN
            CHAR1=CFROMI(no,'(I10)')
            CALL TRIM(CHAR1,IBEG1,IEND1)
            CHAR2=CFROMR(SSCALE,'(E12.6)')
            CALL TRIM(CHAR2,IBEG2,IEND2)
            FORMAT='($,'' Specify the scale associated with'//
     '        ' degree of freedom '//CHAR1(IBEG1:IEND1)//
     '        ' ['//CHAR2(IBEG2:IEND2)//' ]: '',E12.6)'
            CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,SSCALE,
     '        RMIN,RMAX,INFO,ERROR,*9999)
            SCALE(no)=RDATA(1)
          ELSE
            SCALE(no)=SSCALE
          ENDIF
          CNO=CFROMI(no,'(I4)')
          CXO=CFROMR(XO(no),'(E12.5)')
          CSC=CFROMR(SCALE(no),'(E12.5)')
          FORMAT='('' Optimisation freedom'//CNO//
     '      ' has value'//CXO//
     '      ' and scale'//CSC//''')'
          WRITE(OP_STRING,FORMAT)
      	  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

      CALL EXITS('IPFIT')
      RETURN
 9999 CALL ERRORS('IPFIT',ERROR)
      CALL EXITS('IPFIT')
      RETURN 1
      END

      SUBROUTINE IPFOUR(ERROR,*)

C#### Subroutine: IPFOUR
C###  Description:
C###    IPFOUR does input for Fourier Transform analysis parameters.

C**** DF is the logarithmic frequency increment in LPT.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:mxch.inc'
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b15.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:ktyp00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IBEG,ICHAR,IEND,INFO,NOQUES
      LOGICAL FILEIP

      CALL ENTERS('IPFOUR',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      FORMAT='('' Specify lower & upper frequency limits (Hz) &'''//
     '  '/$,'' log frequency increment [1,100,1]: '',3D11.4)'
      IF(IOTYPE.EQ.3) THEN
        RDATA(1)=F0
        RDATA(2)=F1
        RDATA(3)=dFreq
      ENDIF
      RDEFLT(1)=1.0D0
      RDEFLT(2)=100.0D0
      RDEFLT(3)=1.0D0
      CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,3,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        F0=RDATA(1)
        F1=RDATA(2)
        dFreq=RDATA(3)
      ENDIF
      FORMAT='('' Specify type of driving function [1]: '''//
     '  '/''   (1) Impulse function '''//
     '  '/''   (2) Step function '''//
     '  '/''   (3) Sine wave '''//
     '  '/''   (4) White noise '''//
     '  '/$,''    '',I1)'
      IF(IOTYPE.EQ.3) IDATA(1)=KTYP25
      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) KTYP25=IDATA(1)

      IF(KTYP25.EQ.1) THEN !Impulse function

      ELSE IF(KTYP25.EQ.2) THEN

      ELSE IF(KTYP25.EQ.3) THEN
        FORMAT='($,'' Specify frequency,amplitude & phase'//
     '    ' [1.0, 1.0, 1.0]: '',3D11.4)'
        IF(IOTYPE.EQ.3) THEN
          RDATA(1)=FREQ
          RDATA(2)=AMPL
          RDATA(3)=PHASE
        ENDIF
        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RONE,-RMAX,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          FREQ=RDATA(1)
          AMPL=RDATA(2)
          PHASE=RDATA(3)
        ENDIF

      ELSE IF(KTYP25.EQ.4) THEN
        FORMAT='($,'' Specify amplitude [1.0]: '',D11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=AMPL
        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RONE,-RMAX,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) AMPL=RDATA(1)

      ELSE IF(KTYP25.EQ.5) THEN

      ELSE IF(KTYP25.EQ.6) THEN

      ENDIF

      FORMAT='($,'' Specify mass and stiffness damping '
     '  //'coeffs [0,0]: '',2D11.4)'
      IF(IOTYPE.EQ.3) THEN
        RDATA(1)=DAMPING_FACTOR1
        RDATA(2)=DAMPING_FACTOR2
      ENDIF
      CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,2,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        DAMPING_FACTOR1=RDATA(1)
        DAMPING_FACTOR2=RDATA(2)
      ENDIF

      FORMAT='($,'' Specify ang. freq. for damping [0]: '',D11.4)'
      IF(IOTYPE.EQ.3) RDATA(1)=DAMPING_FREQUENCY
      CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) DAMPING_FREQUENCY=RDATA(1)

      CALL STRING_TRIM(FILE00,IBEG,IEND)
      FORMAT='($,'' Enter output file name for nodal var.s ['
     '  //FILE00(IBEG:IEND)//'] (file ext is .FOURIER): '',A)'
      CDEFLT(1)=FILE00
      IF(IOTYPE.EQ.3) CDATA(1)=FILE02
      CALL GINOUT(IOTYPE,2,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,30,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) FILE02=CDATA(1)(1:100)

      CALL EXITS('IPFOUR')
      RETURN
 9999 CALL ERRORS('IPFOUR',ERROR)
      CALL EXITS('IPFOUR')
      RETURN 1
      END


      SUBROUTINE IPSHEE(IDO,NBJ,NEELEM,NKJ,NPNODE,nr,NVJE,NVJP,
     '  XA,XP,ERROR,*)

C#### Subroutine: IPSHEE
C###  Description:
C###    IPSHEE inputs sheet direction field for region nr.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:fibr01.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
      INCLUDE 'cmiss$reference:loc00.inc'
!     Parameter List
      INTEGER IDO(NKM,NNM,0:NIM,NBFM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NKJ(NJM,NPM),NPNODE(0:NP_R_M,0:NRM),nr,
     '  NVJE(NNM,NBFM,NJM,NEFM),NVJP(NJM,NPM)
      REAL*8 XA(NAM,NJM,NQM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER FREELIST(NJ_LOC_MX),IBEG1,ICHAR,IEND1,INFO,
     '  ne,nj,nk,NKJP(6),nn,noelem,nonode,
     '  NOQUES,np,NTNODE,nu,NUK(8),numf,NUMFREE,nv
      CHARACTER CHAR1*5,CHAR2*1,CHAR3*1,CHAR4*12
      LOGICAL BYELEMSHEE,ELEM,FILEIP,PROMPT_NV

      DATA NUK/1,2,4,6,7,9,10,11/

      CALL ENTERS('IPSHEE',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      CALL ASSERT(NJ_LOC(NJL_FIBR,0,nr).GE.1,
     '  '>>Define fibre AND imbrication angles first',ERROR,*9999)

      IF(NET(nr).GT.0) THEN
        ELEM=.TRUE.
      ELSE
        ELEM=.FALSE.
      ENDIF

      IF(IOTYPE.NE.3) THEN

C ***   Check to see if there is enough room to store the sheets 
C ***   (& imbrication angles if not already defined)
        nj=0
        NUMFREE=0
C ***   Will only check for room if JTYP9 does not already equal 3
        DO WHILE((NUMFREE.NE.3-JTYP9).AND.(nj.LE.NJM))
          nj=nj+1
          IF(NJ_TYPE(nj,1).EQ.0.OR.(NJ_TYPE(nj,1).EQ.NJL_FIBR.AND.
     '      NJ_TYPE(nj,2).EQ.3)) THEN  
C ***       Empty space or sheet info in that nj location
            NUMFREE=NUMFREE+1
            FREELIST(NUMFREE)=NJ
          ENDIF
      	  CALL ASSERT(nj.LE.NJM,' >>Increase NJMX',ERROR,*9999)
        ENDDO
        
C ***   Clear any existing sheets unless doing an add
c        nj=NJ_LOC(NJL_FIBR,3)
c        IF(nj.GT.0) THEN
c          NJ_TYPE(nj,1)=0
c          NJ_TYPE(nj,2)=0
c          NJ_LOC(NJL_FIBR,3)=0
c          NJ_LOC(NJL_FIBR,0)=2
c          IF(nj.EQ.NJ_LOC(0,0)) THEN
c            NJ_LOC(0,0)=0
c            DO njj1=1,3
c              DO njj2=1,NJ_LOC(njj1,0)
c                nj=NJ_LOC(njj1,njj2)
c                IF(nj.GT.NJ_LOC(0,0)) NJ_LOC(0,0)=nj
c              ENDDO
c            ENDDO
c          ENDIF
c        ENDIF
C ***   Store sheets (and imbric angles if not defined) in free space
C ***   Will not do loop if jtyp9 already equals 3
        DO numf=1,NUMFREE
          nj=FREELIST(numf)
          NJ_LOC(NJL_FIBR,numf+JTYP9)=nj
          NJ_TYPE(nj,1)=NJL_FIBR
          NJ_TYPE(nj,2)=numf+JTYP9
          IF(nj.GT.NJ_LOC(0,0)) NJ_LOC(0,0)=nj  
        ENDDO
        NJ_LOC(NJL_FIBR,0)=3
        JTYP9=3
      ENDIF

      FORMAT='('' Specify whether angles entered in [1]: '''//
     '  '/''   (1) degrees'''//
     '  '/''   (2) radians'''// 
     '  '/$,''    '',I1)'
      IF(IOTYPE.EQ.3) IDATA(1)=JTYP13
      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,2,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) JTYP13=IDATA(1)

C cpb 10/3/97 fibre element information now set up witht a define
C element fibre

CCcc cpb 16/9/95 Adding element dependent sheet basis types
C      FORMAT='($,'' Is the basis function for the sheet angle '
C     '    //'element dependent [N]? '',A)'
C      CDEFLT(1)='N'
C      IF(IOTYPE.EQ.3) THEN
C        ADATA(1)='N'
C        nb1=NBJ(NJ_LOC(NJL_FIBR,3),NEELEM(1,nr))
C        DO noelem=2,NEELEM(0,nr)
C          ne=NEELEM(noelem,nr)
C          nb2=NBJ(NJ_LOC(NJL_FIBR,3),ne)
C          IF(nb2.NE.nb1) ADATA(1)='Y'
C        ENDDO !ne
C      ENDIF
C      CALL GINOUT(IOTYPE,1,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '  ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
C     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C      IF(ADATA(1).EQ.'N') THEN
C        BYELEMSHEE=.FALSE.
C        FORMAT='($,'' The basis function type number for the '
C     '    //'sheet angle is [1]: '',I1)'
C        IF(IOTYPE.EQ.3) THEN
C          NB_SHEET=NBJ(NJ_LOC(NJL_FIBR,3),NEELEM(1,nr))
C          IDATA(1)=NB_SHEET
C        ENDIF
C        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NBT,
C     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C        IF(IOTYPE.NE.3) THEN
C          NB_SHEET=IDATA(1)
C          DO noelem=1,NEELEM(0,nr)
C            ne=NEELEM(noelem,nr)
C            NB_GEOM=NBJ(1,ne)
C            IF(NIT(NB_GEOM).NE.NIT(NB_SHEET)) THEN
C              NB_SHEET=NB_GEOM
C            ENDIF
C            NBJ(NJ_LOC(NJL_FIBR,3),ne)=NB_SHEET
C          ENDDO !noelem (ne)
C          DO nonode=1,NPNODE(0,nr)
C            np=NPNODE(nonode,nr)
C            NKJ(NJ_LOC(NJL_FIBR,3),np)=NKT(0,NB_SHEET)
C          ENDDO !nonode (np)
C        ENDIF
C      ELSE
C        CALL ASSERT(CALL_ELEM,'>>Define elements first',ERROR,*9999)
C        BYELEMSHEE=.TRUE.
C        NKTMAX=1
C        DO noelem=1,NEELEM(0,nr)
C          ne=NEELEM(noelem,nr)
C          IDEFLT(1)=NBJ(1,ne)
C          WRITE(CHAR2,'(I1)') IDEFLT(1)
C          WRITE(CHAR6,'(I5)') ne
C          FORMAT='($,'' The basis function type number for the sheet '
C     '      //'angle in element '//CHAR6(1:5)//' is ['//CHAR2//']: '','
C     '      //'I1)'
C          IF(IOTYPE.EQ.3) THEN
C            NB_SHEET=NBJ(NJ_LOC(NJL_FIBR,3),ne)
C            IDATA(1)=NB_SHEET
C          ENDIF
C          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NBT,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C          IF(IOTYPE.NE.3) THEN
C            NB_SHEET=IDATA(1)
C            NB_GEOM=NBJ(1,ne)
C            IF(NIT(NB_GEOM).NE.NIT(NB_SHEET)) THEN
C              NB_SHEET=NB_GEOM
C            ENDIF
C            NBJ(NJ_LOC(NJL_FIBR,3),ne)=NB_SHEET
C          ENDIF !iotype
C          IF(NKT(0,NB_SHEET).GT.NKTMAX) NKTMAX=NKT(0,NB_SHEET)
C        ENDDO !noelem
C        DO nonode=1,NPNODE(0,nr)
C          np=NPNODE(nonode,nr)
C          NKJ(NJ_LOC(NJL_FIBR,3),np)=NKTMAX
C        ENDDO !nonode (np)
C      ENDIF
C
C      IF(NNT(NB_SHEET).GT.0) THEN

        IDEFLT(1)=NPNODE(0,nr)
        WRITE(CHAR1,'(I5)') IDEFLT(1)
        IF(IOTYPE.EQ.3) THEN
          NTNODE=NPNODE(0,nr)
          IDATA(1)=NTNODE
        ENDIF
        FORMAT='($,'' The number of nodes is ['//CHAR1(1:5)//']: '',I5)'
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,NPM,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NTNODE=IDATA(1)

!news MPN 16-Nov-94
        !ask for version prompting
        FORMAT='($,'' Do you want prompting for different versions '
     '    //'of the sheet field [N]? '',A)'
        IF(IOTYPE.EQ.3) THEN
          PROMPT_NV=.FALSE.
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            IF(NVJP(NJ_LOC(NJL_FIBR,3),np).GT.1) PROMPT_NV=.TRUE.
          ENDDO
          IF(PROMPT_NV) THEN
            ADATA(1)='Y'
          ELSE
            ADATA(1)='N'
          ENDIF
        ENDIF
        CALL GINOUT(IOTYPE,1,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          IF(ADATA(1).EQ.'Y') THEN
            PROMPT_NV=.TRUE.
          ELSE
            PROMPT_NV=.FALSE.
          ENDIF
        ENDIF
!newe

        IDEFLT(1)=NKJ(NJ_LOC(NJL_FIBR,3),NPNODE(1,nr))-1
        WRITE(CHAR1,'(I1)') IDEFLT(1)
        IF(IOTYPE.EQ.3) THEN
          NKJP(NJ_LOC(NJL_FIBR,3))=NKJ(NJ_LOC(NJL_FIBR,3),NPNODE(1,nr))
          IDATA(1)=NKJP(NJ_LOC(NJL_FIBR,3))-1
        ENDIF
        FORMAT='($,'' The number of derivatives for the sheet field'
     '    //' is ['//CHAR1(1:1)//']: '',I1)'
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,NKM-1,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NKJP(NJ_LOC(NJL_FIBR,3))=IDATA(1)+1

C      ELSE
C        NTNODE=0
C      ENDIF

      IF(NTNODE.GT.0) THEN
        DO nonode=1,NTNODE
          IDEFLT(1)=NPNODE(nonode,nr)
          WRITE(CHAR1,'(I5)') IDEFLT(1)
          IF(IOTYPE.EQ.3) THEN
            np=NPNODE(nonode,nr)
            IDATA(1)=np
          ENDIF
          FORMAT='(/$,'' Node number ['//CHAR1(1:5)//']: '',I5)'
          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NPM,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) np=IDATA(1)

          WRITE(CHAR1,'(I1)') NJ_LOC(NJL_FIBR,3)
!news MPN 16-Nov-94
          IF(PROMPT_NV) THEN  !prompt for diff versions of sheets
            FORMAT='($,'' The number of versions for the sheet'
     '        //' field is [1]: '',I2)'
            IF(IOTYPE.EQ.3) IDATA(1)=NVJP(NJ_LOC(NJL_FIBR,3),np)
            CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '        1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NVM,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) NVJP(NJ_LOC(NJL_FIBR,3),np)=IDATA(1)
          ELSE
            NVJP(NJ_LOC(NJL_FIBR,3),np)=1
          ENDIF
                    
          DO nv=1,NVJP(NJ_LOC(NJL_FIBR,3),np)
            IF(NVJP(NJ_LOC(NJL_FIBR,3),np).GT.1) THEN 
              !ask for diff nj versions
              WRITE(CHAR1,'(I2)') nv
              FORMAT='('' For version number'//CHAR1(1:2)//':'')'
              CALL GINOUT(IOTYPE,0,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '          0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '          ERROR,*9999)
c             CALL INOUT(IOTYPE,IVDU,IFILE,FORMAT,ERROR,*9999)
            ENDIF
!newe
            DO nk=1,NKJP(NJ_LOC(NJL_FIBR,3))
              IF(np.EQ.1) THEN
                RDEFLT(1)=0.0D0
              ELSE IF(np.GT.1) THEN
                IF(JTYP13.EQ.1) THEN      !in degrees
                  RDEFLT(1)=XP(nk,nv,NJ_LOC(NJL_FIBR,3),NP-1)*180.0D0/PI
                ELSE IF(JTYP13.EQ.2) THEN !in radians
                  RDEFLT(1)=XP(nk,nv,NJ_LOC(NJL_FIBR,3),NP-1)
                ENDIF
              ENDIF
              WRITE(CHAR4,'(E12.5)') RDEFLT(1)
c CPB 17/9/95 Only use non-standard mapping if there is one basis
C!news MPN 13-Dec-94: always want to use IDO here
C                nu=IDO(nk,1,0,NB_SHEET)
C!old
C!              IF(ELEM) THEN
C!                nu=IDO(nk,1,0,NB_SHEET)
C!              ELSE
C!                nu=NUK(nk)
C!              ENDIF

C MLB 16/4/97 BYELEMSHEE is not set.
C              
              CALL ASSERT(.FALSE.,
     '          'ERROR, code has variable not set correctly',
     '          ERROR,*9999)
              BYELEMSHEE=.FALSE.
C 
              IF(BYELEMSHEE) THEN
                nu=NUK(nk)
              ELSE
                nu=IDO(nk,1,0,NB_SHEET)
              ENDIF
              IF(nu.EQ.1) THEN
                FORMAT='($,'' The Xj('//CHAR1(1:1)//') coordinate'//
     '            ' is ['//CHAR4(1:12)//']: '',E12.5)'
              ELSE IF(nu.EQ.2.OR.nu.EQ.4.OR.nu.EQ.7) THEN
                IF(nu.EQ.2) THEN
                  CHAR2='1'
                ELSE IF(nu.EQ.4) THEN
                  CHAR2='2'
                ELSE IF(nu.EQ.7) THEN
                  CHAR2='3'
                ENDIF
                FORMAT='($,'' The Xj('//CHAR1(1:1)//') derivative'//
     '            ' wrt s('//CHAR2(1:1)//') is ['//
     '            CHAR4(1:12)//']: '',E12.5)'
              ELSE IF(nu.EQ.6.OR.nu.EQ.9.OR.nu.EQ.10) THEN
                IF(nu.EQ.6) THEN
                  CHAR2='1'
                  CHAR3='2'
                ELSE IF(nu.EQ.9) THEN
                  CHAR2='2'
                  CHAR3='3'
                ELSE IF(nu.EQ.10) THEN
                  CHAR2='3'
                  CHAR3='1'
                ENDIF
                FORMAT='($,'' The Xj('//CHAR1(1:1)//') derivative'//
     '            ' wrt s('//CHAR2(1:1)//') & s('//CHAR3(1:1)//
     '            ') is ['//CHAR4(1:12)//']: '',E12.5)'
              ELSE IF(nu.EQ.11) THEN
                FORMAT='($,'' The Xj('//CHAR1(1:1)//') derivative'//
     '            ' wrt s(1), s(2) & s(3) is ['//
     '            CHAR4(1:12)//']: '',E12.5)'
              ENDIF
              IF(IOTYPE.EQ.3) THEN
                IF(JTYP13.EQ.1) THEN
                  RDATA(1)=XP(nk,nv,NJ_LOC(NJL_FIBR,3),np)*180.0D0/PI
                ELSE
                  RDATA(1)=XP(nk,nv,NJ_LOC(NJL_FIBR,3),np)
                ENDIF
              ENDIF
              CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,
     '          ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '          LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
                IF(JTYP13.EQ.1) THEN      !in degrees
                  XP(nk,nv,NJ_LOC(NJL_FIBR,3),np)=RDATA(1)*PI/180.0D0
                ELSE IF(JTYP13.EQ.2) THEN !in radians
                  XP(nk,nv,NJ_LOC(NJL_FIBR,3),np)=RDATA(1)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
!news MPN 16-Nov-94
        IF(IOTYPE.NE.3) THEN
          IF(ELEM) THEN
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              DO nn=1,NNT(NB_SHEET)
                NVJE(nn,NB_SHEET,NJ_LOC(NJL_FIBR,3),ne)=1 !deflt vers. #
              ENDDO  !nn
            ENDDO !noelem
            IF(PROMPT_NV) THEN
              WRITE(OP_STRING,'(''>>Redefine elements to use correct '
     '          //'versions for sheet field'')')
      	      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF
        ENDIF
!newe
      ELSE IF(NAT(NB_SHEET).GT.0) THEN
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          IF(NIT(NBJ(1,ne)).GT.1) THEN
            WRITE(CHAR1,'(I4)') ne
            CALL TRIM(CHAR1,IBEG1,IEND1)
            FORMAT='('' element number '//CHAR1(IBEG1:IEND1)//': '')'
            CALL GINOUT(IOTYPE,0,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C CPB 29/3/94 Adding NJ_LOC
C            nj=NJE(ne)+JTYP9
            nj=NJ_LOC(NJL_FIBR,3)
            WRITE(CHAR1,'(I1)') nj
            IF(ne.EQ.1) THEN
              RDEFLT(1)=0.0d0
            ELSE IF(ne.GT.1) THEN
              RDEFLT(1)=XA(1,nj,ne-1)
            ENDIF
            WRITE(CHAR4,'(E12.5)') RDEFLT(1)
            FORMAT='($,'' the Xj('//CHAR1(1:1)//
     '        ') coordinate is ['//CHAR4(1:12)//']: '',E12.5)'
            IF(IOTYPE.EQ.3) THEN
              IF(JTYP13.EQ.1) THEN
                RDATA(1)=XA(1,nj,ne)*180.0d0/PI
              ELSE IF(JTYP13.EQ.2) THEN
                RDATA(1)=XA(1,nj,ne)
              ENDIF
            ENDIF
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              IF(JTYP13.EQ.1) THEN
                XA(1,nj,ne)=RDATA(1)*PI/180.0d0
              ELSE IF(JTYP13.EQ.2) THEN
                XA(1,nj,ne)=RDATA(1)
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('IPSHEE')
      RETURN
 9999 CALL ERRORS('IPSHEE',ERROR)
      CALL EXITS('IPSHEE')
      RETURN 1
      END

C18-nov-1999
      SUBROUTINE IPTIME(IBT,IDO,INP,NAN,NGAP,VARIABLE,END_TIME,
     '  START_TIME,END_TIME_SPECIFIED,START_TIME_SPECIFIED,ERROR,*)

C#### Subroutine: IPTIME
C###  Description:
C###    Inputs time variable information
C *** Created: David Nickerson, 10 August 1999

      IMPLICIT NONE

      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:time_variable.cmn'
      INCLUDE 'cmiss$reference:tol00.cmn'

      !Parameter list
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NAN(NIM,NAM,NBFM),NGAP(NIM,NBM),VARIABLE
      REAL*8 END_TIME,START_TIME
      LOGICAL END_TIME_SPECIFIED,START_TIME_SPECIFIED
      CHARACTER ERROR*(*)
      !Local variables
      INTEGER END_VAR,MAX_OPT,START_VAR,var,i,input_var,IBEG,IEND,
     '  ICHAR,INFO,NOQUES,node,elem,current_var,nb_found,nb,
     '  POSITION(4),ni,nn,IDO1,IDO2,IDO3,nk
      REAL*8 SCALE_FACTOR
      CHARACTER CHAR*1,CHAR2*4,CHAR3*20
      LOGICAL FILEIP,BASIS_FOUND
      
      CALL ENTERS('IPTIME',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      INFO=0
      ICHAR=999

C *** When reading, START_VAR and END_VAR will be the same so only do
C *** the var loop once, but loop through prompting for time variables 
C *** until user exit's increasing TV_NUM_VARIABLES for each variable.
C ***
C *** When writing the variable(s) the var loop is used to loop through
C *** all required variables, and writes each one out.


      IF(IOTYPE.NE.3) THEN !reading
        START_VAR=TV_NUM_VARIABLES+1
        END_VAR=START_VAR
      ELSE !writing
        CALL ASSERT(TV_NUM_VARIABLES.GT.0,'>>no time variables defined',
     '    ERROR,*9999)
        IF(VARIABLE.GT.0) THEN
          START_VAR=VARIABLE
          END_VAR=VARIABLE
        ELSE
          START_VAR=1
          END_VAR=TV_NUM_VARIABLES
        ENDIF
      ENDIF

      CALL ASSERT(START_VAR.LE.TV_NUM_VARIABLES+1.AND.
     '  END_VAR.LE.TV_NUM_VARIABLES+1,'Invalid variable number',
     '  ERROR,*9999)

      var=START_VAR

      CDEFLT(1)='EXIT'
 6100 FORMAT='($,'' Enter the variable name [EXIT]: '',A20)'
      IF(IOTYPE.EQ.3) CDATA(1)=TV_NAMES(var)
      CALL GINOUT(IOTYPE,2,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,MAX_TIME_VARIABLE_NAME,IDATA,
     '  IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '  ERROR,*9999)
      IF(CDATA(1).NE.'EXIT'.AND.var.LE.END_VAR) THEN !not default exit 
                                                     !or more variables
                                                     !to write out.
        IF(IOTYPE.NE.3) THEN
          TV_NUM_VARIABLES=TV_NUM_VARIABLES+1
          IF(TV_NUM_VARIABLES.GT.MAX_TIME_VARIABLES) THEN
            ERROR='Exceeded maximum number of time variables allowed'
            TV_NUM_VARIABLES=TV_NUM_VARIABLES-1
            GOTO 9999
          ELSE
            input_var=TV_NUM_VARIABLES
          ENDIF
          DO i=1,TV_NUM_VARIABLES
            IF(input_var.NE.i.AND.
     '        CDATA(1)(1:MAX_TIME_VARIABLE_NAME).EQ.TV_NAMES(i)) THEN
              TV_NUM_VARIABLES=TV_NUM_VARIABLES-1
              WRITE(ERROR,'(''Name already exists: '',A20)') 
     '          CDATA(1)(1:MAX_TIME_VARIABLE_NAME)
              GOTO 9999
            ENDIF
          ENDDO
          TV_NAMES(input_var)=CDATA(1)(1:MAX_TIME_VARIABLE_NAME)
        ENDIF
        
C ***   Get the basis type
        IDEFLT(1)=1
        WRITE(CHAR,'(I1)') IDEFLT(1)
        FORMAT='('' Basis Type ['//CHAR//']:'''//
     '    '/''   (1) Linear Lagrange    '''//
     '    '/''   (2) Quadratic Lagrange '''//
     '    '/''   (3) Cubic Lagrange     '''//
     '    '/''  *(4) Cubic Hermite      '''//
     '    '/$,''    '',I1)'
        MAX_OPT=3
        IF(IOTYPE.EQ.3) IDATA(1)=TV_BASIS_TYPES(var)
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,MAX_OPT,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) TV_BASIS_TYPES(input_var)=IDATA(1)

C ***   Check if the requested basis exists and add it if it doesn't
        
        TV_BASIS(input_var)=-1 !initialise
        BASIS_FOUND=.FALSE.
        !loop through all existing bases
        DO nb=1,NBT
          !check the number of local Xi-coords for basis nb - only 
          !want 1-D
          IF(NIT(nb).EQ.1) THEN
            !check for required basis
            IF(IBT(1,1,nb).EQ.1) THEN !lagrange
              IF(IBT(2,1,nb).EQ.1) THEN !linear lagrange
                IF(TV_BASIS_TYPES(input_var).EQ.TV_LINEAR_LAGRANGE) THEN
                  nb_found=nb
                  BASIS_FOUND=.TRUE.
                ENDIF
              ELSEIF(IBT(2,1,nb).EQ.2) THEN !quadratic lagrange
                IF(TV_BASIS_TYPES(input_var).EQ.TV_QUADRATIC_LAGRANGE) 
     '            THEN
                  nb_found=nb
                  BASIS_FOUND=.TRUE.
                ENDIF
              ELSEIF(IBT(2,1,nb).EQ.3) THEN !cubic lagrange
                IF(TV_BASIS_TYPES(input_var).EQ.TV_CUBIC_LAGRANGE) THEN
                  nb_found=nb
                  BASIS_FOUND=.TRUE.
                ENDIF
              ENDIF !lagrange basis type
            ENDIF !lagrange basis
          ENDIF !same # xi coords
        ENDDO !nb

        IF(BASIS_FOUND) THEN !basis already defined
          TV_BASIS(input_var)=nb_found
        ELSE !create the required basis

          IF(NBT+1.GT.NBM) THEN
            ERROR=' >>Increase NBM'
            TV_NUM_VARIABLES=TV_NUM_VARIABLES-1
            GOTO 9999
          ENDIF
          IF(NBFT+1.GT.NBFM) THEN
            ERROR=' >>Increase NBFM'
            TV_NUM_VARIABLES=TV_NUM_VARIABLES-1
            GOTO 9999
          ENDIF

          IF(TV_BASIS_TYPES(input_var).EQ.TV_LINEAR_LAGRANGE) THEN
            WRITE(OP_STRING,'('' Adding a linear lagrange basis'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSEIF(TV_BASIS_TYPES(input_var).EQ.TV_QUADRATIC_LAGRANGE) 
     '        THEN
            WRITE(OP_STRING,'('' Adding a quadratic lagrange basis'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSEIF(TV_BASIS_TYPES(input_var).EQ.TV_CUBIC_LAGRANGE) THEN
            WRITE(OP_STRING,'('' Adding a cubic lagrange basis'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
              
          nb=NBT+1
          NBT=NBT+1
          NBFT=NBFT+1

          NBC(nb)=1
          NBI(nb)=1
          NBASEF(nb,0)=1
          NBASEF(nb,1)=nb
          NBASEF(nb,2)=nb
          NFBASE(1,nb)=nb
          NFBASE(2,nb)=nb

          NAT(nb)=0 !aux. bases
          NAN(1,1,nb)=0
          NIT(nb)=1 !xi dirns
          NUT(nb)=NIT(nb)*NIT(nb)+2
          IF(NUT(nb).GT.NUM) THEN
            ERROR=' >>Increase NUM'
            TV_NUM_VARIABLES=TV_NUM_VARIABLES-1
            GOTO 9999
          ENDIF
          DO ni=1,NIT(nb)
            IBT(1,ni,nb)=1 !lagrange
            IF(TV_BASIS_TYPES(input_var).EQ.TV_LINEAR_LAGRANGE) THEN
              IBT(2,ni,nb)=1 !linear
            ELSEIF(TV_BASIS_TYPES(input_var).EQ.TV_QUADRATIC_LAGRANGE) 
     '          THEN
              IBT(2,ni,nb)=2 !quadratic
            ELSEIF(TV_BASIS_TYPES(input_var).EQ.TV_CUBIC_LAGRANGE) THEN
              IBT(2,ni,nb)=3 !cubic
            ENDIF
            NGAP(ni,nb)=0 !gauss pts
          ENDDO !ni
          NGT(nb)=0
          !set the number of nodes in each element
          NNT(nb)=TV_NODES_PER_ELEMENT(TV_BASIS_TYPES(input_var))
          IF(NNT(nb).GT.NNM) THEN
            ERROR=' >>Increase NNM'
            TV_NUM_VARIABLES=TV_NUM_VARIABLES-1
            GOTO 9999
          ENDIF

          !INP
          DO ni=1,4
            POSITION(ni)=1
          ENDDO !ni
          DO nn=1,NNT(nb)
            ni=1
            DO WHILE(POSITION(ni).GT.NNT(nb))
              POSITION(ni)=1
              ni=ni+1
              POSITION(ni)=POSITION(ni)+1
            ENDDO !ni
            DO ni=1,NIT(nb)
              INP(nn,ni,nb)=POSITION(ni)
            ENDDO !ni
            POSITION(1)=POSITION(1)+1
          ENDDO !nn

          !IDO
          NST(nb)=0
          NKT(0,nb)=1 !the maximum number of nodal derivatives
          DO nn=1,NNT(nb)
            NKT(nn,nb)=1
            DO ni=1,NIT(nb)
              DO nk=1,NKT(nn,nb)
                IDO(nk,nn,ni,nb)=1
              ENDDO !nk
            ENDDO !ni
            IF(NKT(nn,nb).GT.NKT(0,nb)) NKT(0,nb)=NKT(nn,nb)
            NST(nb)=NST(nb)+NKT(nn,nb)
          ENDDO !nn
          IF(NKT(0,nb).GT.NKM) THEN
            ERROR=' >>Increase NKM'
            TV_NUM_VARIABLES=TV_NUM_VARIABLES-1
            GOTO 9999
          ENDIF
          IF(NST(nb).GT.NSM) THEN
            ERROR=' >>Increase NSM'
            TV_NUM_VARIABLES=TV_NUM_VARIABLES-1
            GOTO 9999
          ENDIF
          DO nn=1,NNT(nb)
            DO nk=1,NKT(nn,nb)
              IDO1=IDO(nk,nn,1,nb)
              IDO2=1
              IDO3=1
              IF(NIT(nb).GE.2) IDO2=IDO(nk,nn,2,nb)
              IF(NIT(nb).EQ.3) IDO3=IDO(nk,nn,3,nb)
              IF(IDO1.EQ.1.AND.IDO2.EQ.1.AND.IDO3.EQ.1)
     '          IDO(nk,nn,0,nb)=1
              IF(IDO1.EQ.2.AND.IDO2.EQ.1.AND.IDO3.EQ.1)
     '          IDO(nk,nn,0,nb)=2
              IF(IDO1.EQ.1.AND.IDO2.EQ.2.AND.IDO3.EQ.1)
     '          IDO(nk,nn,0,nb)=4
              IF(IDO1.EQ.2.AND.IDO2.EQ.2.AND.IDO3.EQ.1)
     '          IDO(nk,nn,0,nb)=6
              IF(IDO1.EQ.1.AND.IDO2.EQ.1.AND.IDO3.EQ.2)
     '          IDO(nk,nn,0,nb)=7
              IF(IDO1.EQ.2.AND.IDO2.EQ.1.AND.IDO3.EQ.2)
     '          IDO(nk,nn,0,nb)=9
              IF(IDO1.EQ.1.AND.IDO2.EQ.2.AND.IDO3.EQ.2)
     '          IDO(nk,nn,0,nb)=10
              IF(IDO1.EQ.2.AND.IDO2.EQ.2.AND.IDO3.EQ.2)
     '          IDO(nk,nn,0,nb)=11
            ENDDO !nk
          ENDDO !nn
 
C ***     Assign the basis number
          TV_BASIS(input_var)=nb

        ENDIF !quad basis found/created



C ***   Get the elements
        IDEFLT(1)=1
        WRITE(CHAR,'(I1)') IDEFLT(1)
        FORMAT='($,'' The number of elements is ['//CHAR//']: '',I3)'
        IF(IOTYPE.NE.3) THEN
          MAX_OPT=MAX_TIME_VARIABLE_NODES*
     '      TV_NODES_PER_ELEMENT(TV_BASIS_TYPES(input_var))
        ELSE
          MAX_OPT=MAX_TIME_VARIABLE_NODES*
     '      TV_NODES_PER_ELEMENT(TV_BASIS_TYPES(var))
        ENDIF
        IF(IOTYPE.EQ.3) IDATA(1)=TV_NUM_ELEMS(var)
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,MAX_OPT,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) TV_NUM_ELEMS(input_var)=IDATA(1)

        IF(IOTYPE.NE.3) THEN
          current_var=input_var
        ELSE
          current_var=var
        ENDIF

        DO elem=1,TV_NUM_ELEMS(current_var)
          WRITE(CHAR2,'(I3)') elem
          CALL TRIM(CHAR2,IBEG,IEND)
          FORMAT='(''   Element: '//CHAR2(IBEG:IEND)//''')'
          CALL GINOUT(IOTYPE,0,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,MAX_OPT,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          DO node=1,TV_NODES_PER_ELEMENT(TV_BASIS_TYPES(current_var))
            IF(elem.EQ.1.OR.node.GT.1) THEN
              WRITE(CHAR3,'(I3 I3)') elem,node
              CALL TRIM(CHAR3,IBEG,IEND)
              FORMAT='(''     Node: '//CHAR3(IBEG:IEND)//''')'
              CALL GINOUT(IOTYPE,0,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '          1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     '          MAX_OPT,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '          ERROR,*9999)
C ***         Get the time of the node (if it is the first node in
C ***         the element, then the time is the same as the last node
C ***         in the previous element)
C ***         Also need to check that time is always increasing
              IF(elem.GT.1.OR.node.GT.1) THEN
                RDEFLT(1)=TV_NODE_TIMES(current_var,(elem-1)*
     '            TV_NODES_PER_ELEMENT(TV_BASIS_TYPES(current_var))+
     '            node-1)
              ELSE
                RDEFLT(1)=0.0d0
              ENDIF
              WRITE(CHAR3,'(D12.5)') RDEFLT(1)
              CALL TRIM(CHAR3,IBEG,IEND)
              FORMAT=
     '          '($,''       Time ['//CHAR3(IBEG:IEND)//']: '',D12.5)'
              IF(IOTYPE.EQ.3) RDATA(1)=TV_NODE_TIMES(current_var,
     '          (elem-1)*
     '          TV_NODES_PER_ELEMENT(TV_BASIS_TYPES(current_var))+node)
              CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '          IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '          RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) TV_NODE_TIMES(current_var,(elem-1)*
     '          TV_NODES_PER_ELEMENT(TV_BASIS_TYPES(current_var))+node)=
     '          RDATA(1)
              IF(node.GT.1) THEN
                IF(TV_NODE_TIMES(current_var,(elem-1)*
     '            TV_NODES_PER_ELEMENT(TV_BASIS_TYPES(current_var))+
     '            node).LT.TV_NODE_TIMES(current_var,(elem-1)*
     '            TV_NODES_PER_ELEMENT(TV_BASIS_TYPES(current_var))+
     '            node-1)) THEN
                  ERROR='Time must always increase'
                  TV_NUM_VARIABLES=TV_NUM_VARIABLES-1
                  GOTO 9999
                ENDIF
              ENDIF
C ***         Get the coordinates of the node (if it is the first node 
C ***         in the element, then the coordinate is the same as the 
C ***         last node in the previous element)
              IF(elem.GT.1.OR.node.GT.1) THEN
                RDEFLT(1)=TV_NODE_COORDS(current_var,(elem-1)*
     '            TV_NODES_PER_ELEMENT(TV_BASIS_TYPES(current_var))+
     '            node-1)
              ELSE
                RDEFLT(1)=0.0d0
              ENDIF
              WRITE(CHAR3,'(D12.5)') RDEFLT(1)
              CALL TRIM(CHAR3,IBEG,IEND)
              FORMAT=
     '          '($,''       Coordinate ['//CHAR3(IBEG:IEND)
     '          //']: '',D12.5)'
              IF(IOTYPE.EQ.3) RDATA(1)=TV_NODE_COORDS(current_var,
     '          (elem-1)*
     '          TV_NODES_PER_ELEMENT(TV_BASIS_TYPES(current_var))+node)
              CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '          IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '          RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) TV_NODE_COORDS(current_var,(elem-1)*
     '          TV_NODES_PER_ELEMENT(TV_BASIS_TYPES(current_var))+node)=
     '          RDATA(1)
            ELSE
C ***         Set the coord and time to be the same as the last node
C ***         in the previous element
              TV_NODE_TIMES(current_var,(elem-1)*
     '          TV_NODES_PER_ELEMENT(TV_BASIS_TYPES(current_var))+node)=
     '          TV_NODE_TIMES(current_var,(elem-1)*
     '          TV_NODES_PER_ELEMENT(TV_BASIS_TYPES(current_var))+
     '          node-1)
              TV_NODE_COORDS(current_var,(elem-1)*
     '          TV_NODES_PER_ELEMENT(TV_BASIS_TYPES(current_var))+node)=
     '          TV_NODE_COORDS(current_var,(elem-1)*
     '          TV_NODES_PER_ELEMENT(TV_BASIS_TYPES(current_var))+
     '          node-1)
            ENDIF !elem.OR.node
          ENDDO !node
        ENDDO !elem
        
C ***   Set the start and end times
        IF(IOTYPE.NE.3) THEN
          TV_START_TIME(input_var)=TV_NODE_TIMES(input_var,1)
          TV_END_TIME(input_var)=TV_NODE_TIMES(input_var,
     '      TV_NUM_ELEMS(input_var)*
     '      TV_NODES_PER_ELEMENT(TV_BASIS_TYPES(input_var)))
        ENDIF

        IF(IOTYPE.EQ.3) var=var+1

        GOTO 6100 !for more variables
        
      ENDIF !CDATA(1).NE.'EXIT'

      IF(IOTYPE.NE.3) THEN !not writing
        IF(START_TIME_SPECIFIED) THEN
          IF(DABS(START_TIME-TV_START_TIME(input_var)).GT.ZERO_TOL) THEN
C ***       The specified start time is different from that read in,
C ***       so need to scale the times of the given nodes.
            SCALE_FACTOR=(TV_END_TIME(input_var)-
     '        TV_START_TIME(input_var))/(TV_END_TIME(input_var)-
     '        START_TIME)
            TV_NODE_TIMES(input_var,1)=START_TIME
            TV_START_TIME(input_var)=START_TIME
            DO node=2,TV_NUM_ELEMS(input_var)*
     '        TV_NODES_PER_ELEMENT(TV_BASIS_TYPES(input_var))
              TV_NODE_TIMES(input_var,node)=START_TIME+
     '          (TV_NODE_TIMES(input_var,node)-START_TIME)*
     '          SCALE_FACTOR
            ENDDO !node
            TV_END_TIME(input_var)=TV_NODE_TIMES(input_var,
     '        TV_NUM_ELEMS(input_var)*
     '        TV_NODES_PER_ELEMENT(TV_BASIS_TYPES(input_var)))
          ENDIF
        ENDIF !start time specified

        IF(END_TIME_SPECIFIED) THEN
          IF(DABS(END_TIME-TV_END_TIME(input_var)).GT.ZERO_TOL) THEN
C ***       The specified end time is different from that read in,
C ***       so need to scale the times of the given nodes.
            SCALE_FACTOR=(TV_END_TIME(input_var)-
     '        TV_START_TIME(input_var))/(END_TIME-
     '        TV_START_TIME(input_var))
            TV_NODE_TIMES(input_var,TV_NUM_ELEMS(input_var)*
     '        TV_NODES_PER_ELEMENT(TV_BASIS_TYPES(input_var)))=
     '        END_TIME
            TV_END_TIME(input_var)=END_TIME
            DO node=1,TV_NUM_ELEMS(input_var)*
     '        TV_NODES_PER_ELEMENT(TV_BASIS_TYPES(input_var))-1
              TV_NODE_TIMES(input_var,node)=TV_END_TIME(input_var)-
     '          (END_TIME-TV_NODE_TIMES(input_var,node))*
     '          SCALE_FACTOR
            ENDDO !node
            TV_START_TIME(input_var)=TV_NODE_TIMES(input_var,1)
          ENDIF
        ENDIF !end time specified
      ENDIF !iotype

      CALL EXITS('IPTIME')
      RETURN
 9999 CALL ERRORS('IPTIME',ERROR)
      CALL EXITS('IPTIME')
      RETURN 1
      END


Module FE12
=========== 

C mpn 17July2000:  From OPACTI for HMT/Fading memory
C news MPN 23May2000: adding Steady State HMT (several changes below)
      ELSE IF(KTYP59(nr).EQ.2) THEN    !Steady State HMT
        WRITE(OP_STRING,'('' Steady State HMT'','
     '    //'/''  Max isometric tension at ext.ratio=1 (Tref)     = '','
     '    //'D11.4,'' kPa'','
     '    //'/''  Non-dim. slope parameter for T0 (beta0)         = '','
     '    //'D11.4,'
     '    //'/''  Ref. Hill coeff for [Ca]i sat. curve (n_ref)    = '','
     '    //'D11.4,'
     '    //'/''  Non-dim. slope parameter for Hill coeff (beta1) = '','
     '    //'D11.4,'
     '    //'/''  Ref. pC50 for [Ca]i saturation curve (pC50_ref) = '','
     '    //'D11.4,'
     '    //'/''  Non-dim. slope parameter for pC50 (beta2)       = '','
     '    //'D11.4,'
     '    //'/''  Current [Ca]i at Gauss pt 1 in element 1        = '','
     '    //'D11.4,'' mM'')') 
     '      Tref,T0_beta,HMT_n_ref,HMT_n_beta,
     '      HMT_pC50_ref,HMT_pC50_beta,FEXT(4,1,1)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      ELSE IF(KTYP59(nr).EQ.3) THEN  !Fading Memory
        WRITE(OP_STRING,'(''   Fading Memory Model.'''
     '    //'/''   Time step = '',E11.4,'' s'''
     '    //'/''   Description of Stress/Velocity curve:'''
     '    //'/''   Slope of stress/veloc relation in stretching = '','
     '    //'E11.4,'' kPa.s'''
     '    //'/''   Yield tension to Isometric tension ratio     = '','
     '    //'E11.4,'' kPa'''
     '    //'/''   Static non-linearity parameter "a"           = '','
     '    //'E11.4,'
     '    //'//''   Number of dynamic terms in material response '','
     '    //'''function:'',I3)') DEL_T,TV_SLO,YIELDR,SNLPA,NTACTV  
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(ntactv.ne.0) THEN
          DO nactv=1,NTACTV
            WRITE(OP_STRING,'(''   Coefficient'',I2,'' ='',E11.4,'
     '        //'''     Rate Constant'',I2,'' ='',E11.4,'' per sec'')')
     '        nactv,ACOEFF(nactv),nactv,ALFA(nactv) 
      	    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDIF

C 24/2/97 LC removed from :

C#### Subroutine: OPMAP
C###  Description:
C###    Outputs mapping arrays for region nr.

c cpb 12/3/95 old
C                    ENDDO !nk
C                  ENDDO !nv
C                ENDDO !nh
C              ENDDO !nonode (np)
C              DO noelem=1,NEELEM(0,nr)
C                ne=NEELEM(noelem,nr)
C                DO nh=1,NHE(ne)
C                  IF(nrc.EQ.1) THEN
CC!!! Use the L.H.S. (nc=1) basis to determine the # of equations
C                    nb=NBH(nh,1,ne)
C                  ELSE
C                    nb=NBH(nh,nc,ne)
C                  ENDIF
C                  DO na=1,NAT(nb)
C                    ny=NYNE(na,nh,nrcc,nc,ne)
C                    IF(NONY(0,ny,nrc,nr).EQ.0) THEN
C                      WRITE(OP_STRING,'('' NYNE='',I5,'
C     '                  //''' NONY(0,ny)  ='',I6)')
C     '                  ny,NONY(0,ny,nrc,nr)
C                      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C                    ELSE IF(NONY(0,ny,nrc,nr).GT.0) THEN
C                      WRITE(OP_STRING,'('' NYNE='',I5,'
C     '                  //''' NONY(0,ny)  ='',I6,'
C     '                  //''' CONY(0,ny)='',F10.4)')
C     '                  ny,NONY(0,ny,nrc,nr),CONY(0,ny,nrc,nr)
C                      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C                      WRITE(OP_STRING,'(12X,''NONY(1..,ny)='','
C     '                  //'9(X,I5))') (NONY(noy,ny,nrc,nr),noy=1,
C     '                  NONY(0,ny,nrc,nr))
C                      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C                      WRITE(OP_STRING,'(12X,''CONY(1..,ny)='','
C     '                  //'9(X,F10.4))') (CONY(noy,ny,nrc,nr),noy=1,
C     '                  NONY(0,ny,nrc,nr))
C                      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C                    ENDIF
C                  ENDDO !na
C                ENDDO !nh
C              ENDDO !noelem (ne) 


      SUBROUTINE OPCOLO(IW,ERROR,*)

C#### Subroutine: OPCOLO
C###  Description:
C**** Output workstation index colours.

      IMPLICIT NONE
!     Parameter List
      INTEGER IW
      CHARACTER ERROR*(*)
!     Local Variables

      CALL ENTERS('OPCOLO',*9999)
      CALL LICOLREP(IW,ERROR,*9999)

      CALL EXITS('OPCOLO')
      RETURN
 9999 CALL ERRORS('OPCOLO',ERROR)
      CALL EXITS('OPCOLO')
      RETURN 1
      END


      SUBROUTINE OPFIBR(IDO,NBJ,NEELEM,NKJ,NPNODE,nr,XA,XP,
     '  ERROR,*)

C#### Subroutine: OPFIBR
C###  Description:
C###    OPFIBR outputs fibre direction field for region nr.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
      INCLUDE 'cmiss$reference:loc00.inc'
!     Parameter List
      INTEGER IDO(NKM,NNM,0:NIM,NBFM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NKJ(NJM,NPM),
     '  NPNODE(0:NP_R_M,0:NRM),nr
      REAL*8 XA(NAM,NJM,NQM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IBEG,IBEG1,IEND,IEND1,nb,ne,nj,nk,
     '  NKJP(5),noelem,nonode,np,NTNODE,nu,NUK(8)
      CHARACTER CHAR1*4,CHAR2*1,
     '  CHAR3*1,CHAR4*12
      LOGICAL ELEM

      DATA NUK/1,2,4,6,7,9,10,11/

C!!!  This doesn't work correctly for non-zero numbers of derivatives.
C!!!  This needs to be updated for versions. MPN 6-Jan-95

      CALL ENTERS('OPFIBR',*9999)
      CALL ASSERT(NJ_LOC(njl_fibr,0,nr).GT.0,
     '  '>>Fibre field not defined',ERROR,*9999)
      IF(NET(nr).GT.0) THEN
        ELEM=.TRUE.
      ELSE
        ELEM=.FALSE.
      ENDIF

      IF(NJ_LOC(njl_fibr,0,nr).EQ.1) THEN
        CHAR1='1-2'
        WRITE(CHAR2,'(I1)') JTYP12
      ELSE IF(NJ_LOC(njl_fibr,0,nr).EQ.2) THEN
        CHAR1='2-3'
        WRITE(CHAR2,'(I1)') JTYP12+1
      ENDIF
      WRITE(OP_STRING,
     '  '(/$,''     Fibre field is in '//CHAR1(1:3)//' plane.'')')
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,
     '  '( $,''     Fibre angle is defined wrt Xi('//CHAR2(1:1)//
     '  ') coordinate.'')')
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      IF(JTYP13.EQ.1) THEN
        CHAR4='degrees'
      ELSE IF(JTYP13.EQ.2) THEN
        CHAR4='radians'
      ENDIF
      WRITE(OP_STRING,'( $,''     Angles are entered in '//CHAR4(1:7)
     '  //'.'')')
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C CPB 29/3/94 Adding NJ_LOC
      nb=NBJ(NJ_LOC(NJL_FIBR,1,nr),NEELEM(1,nr))
      WRITE(OP_STRING,
     '  '(/$,''     The basis function type number for the fibre'//
     '  ' field is: '',I2)') nb
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      IF(NNT(nb).GT.0) THEN
        NTNODE=NPNODE(0,nr)
        WRITE(OP_STRING,'( $,''     The number of nodes is '',I4)')
     '    NTNODE
      	CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        NKJP(NJ_LOC(NJL_FIBR,1,nr))=NKJ(NJ_LOC(NJL_FIBR,1,nr),
     '    NPNODE(1,nr))
        WRITE(OP_STRING,
     '    '( $,''     The number of derivatives for the fibre field'
     '    //' is: '',I1)') NKJP(NJ_LOC(NJL_FIBR,1,nr))-1
      	CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSE
        NTNODE=0
      ENDIF

      IF(NTNODE.GT.0) THEN
        DO nonode=1,NTNODE
          np=NPNODE(nonode,nr)
          WRITE(OP_STRING,'(/$,''   Node number: '',I4)') np
      	  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(CHAR1,'(I2)') NJ_LOC(NJL_FIBR,1,nr)
          CALL TRIM(CHAR1,IBEG,IEND)
          DO nk=1,NKJP(NJ_LOC(NJL_FIBR,1,nr))
            IF(ELEM) THEN
              nu=IDO(nk,1,0,nb)
            ELSE
              nu=NUK(nk)
            ENDIF
            IF(nu.EQ.1) THEN
              FORMAT=
     '          '($,''   The Xj('//CHAR1(IBEG:IEND)
     '          //') coordinate is: '',D12.5)'
            ELSE IF(nu.eq.2.or.nu.eq.4.or.nu.EQ.7) THEN
              IF(nu.EQ.2) THEN
                CHAR2='1'
              ELSE IF(nu.EQ.4) THEN
                CHAR2='2'
              ELSE IF(nu.EQ.7) THEN
                CHAR2='3'
              ENDIF
              FORMAT='($,''   The Xj('//CHAR1(IBEG:IEND)
     '          //') derivative'//' wrt s('//CHAR2(1:1)
     '          //') is: '',D12.5)'
            ELSE IF(nu.eq.6.or.nu.eq.9.or.nu.EQ.10) THEN
              IF(nu.EQ.6) THEN
                CHAR2='1'                                    
                CHAR3='2'
              ELSE IF(nu.EQ.9) THEN
                CHAR2='2'
                CHAR3='3'
              ELSE IF(nu.EQ.10) THEN
                CHAR2='3'
                CHAR3='1'
              ENDIF
              FORMAT='($,''   The Xj('//CHAR1(IBEG:IEND)
     '          //') derivative'//' wrt s('//CHAR2(1:1)//') & s('
     '          //CHAR3(1:1)//') is: '',D12.5)'
            ELSE IF(nu.eq.11) THEN
              FORMAT='($,''   The Xj('//CHAR1(IBEG:IEND)
     '          //') derivative'//' wrt s(1), s(2) & s(3)'
     '          //' is: '',D12.5)'
            ENDIF
            IF(JTYP13.EQ.1) THEN
              WRITE(OP_STRING,FORMAT) 
     '          XP(nk,1,NJ_LOC(NJL_FIBR,1,nr),np)*
     '          180.0D0/PI
      	      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSE
              WRITE(OP_STRING,FORMAT) XP(nk,1,NJ_LOC(NJL_FIBR,1,nr),np)
      	      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO
        ENDDO

      ELSE IF(NAT(nb).GT.0) THEN
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          IF(NIT(NBJ(1,ne)).GT.1) THEN
            WRITE(CHAR1,'(I4)') ne
            CALL TRIM(CHAR1,IBEG1,IEND1)
            FORMAT='('' element number '//CHAR1(IBEG1:IEND1)//': '')'
C CPB 29/3/94 Adding NJ_LOC
C            nj=NJE(ne)+JTYP9
            nj=NJ_LOC(NJL_FIBR,1,nr)
            WRITE(CHAR1,'(I2)') nj
            CALL TRIM(CHAR1,IBEG,IEND)
            FORMAT='($,'' the Xj('//CHAR1(IBEG:IEND)
     '        //') coordinate is: '',D12.5)'
            IF(JTYP13.EQ.1) THEN
              WRITE(OP_STRING,FORMAT) XA(1,nj,ne)*180.0d0/PI
      	      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSE IF(JTYP13.EQ.2) THEN
              WRITE(OP_STRING,FORMAT) XA(1,nj,ne)
      	      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('OPFIBR')
      RETURN
 9999 CALL ERRORS('OPFIBR',ERROR)
      CALL EXITS('OPFIBR')
      RETURN 1
      END


      SUBROUTINE OPMATR(ERROR,*)

C#### Subroutine: OPMATR
C###  Description: Output matrix.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:coef00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER NCMATR,nrmatr

      CALL ENTERS('OPMATR',*9999)

      DO nrmatr=1,NT_MATR
        WRITE(OP_STRING,'(10E12.4)')
     '    (AMATR(nrmatr,NCMATR),NCMATR=1,NT_MATR)
      	CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDDO

      CALL EXITS('OPMATR')
      RETURN
 9999 CALL ERRORS('OPMATR',ERROR)
      CALL EXITS('OPMATR')
      RETURN 1
      END


      SUBROUTINE OPST80(NBJ,ne,NJE,PG,VE,XE,XG,ZE,ZG,ERROR,*)

C**** Output stresses and strains at Gauss points. 
C**** Only called for Fourier time strain calculation. 
C**** It's days are numbered.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbst02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:press00.cmn'
      INCLUDE 'cmiss$reference:trans00.cmn'
!     Parameter List
      INTEGER NBJ(*),ne,NJE
      REAL*8 PG(NSM,NUM,NGM,*),VE(NSM,*),
     '  XE(NSM,*),XG(NJM,*),ZE(NSM,*),ZG(NHM,*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J,K,ng,NGTB,nj,nu
      REAL*8 AXL(3,3),AZL(3,3),DUJNU(3,3),DXINU(3,3),
     '  DXJNU(3,3),
     '  EG(3,3),ETA1,ETA2,EVAL(3),EVEC(3),GL(3,3),GU(3,3),PHI,
     '  RGX,VG(4,10),XI(3)

      CALL ENTERS('OPST80',*9999)
C ***   Fourier motion:
        NGTB=NGT(NBJ(1))

        WRITE(OP_STRING,'(/'' Element '',I3,'' (Fourier):'')') NE
      	CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO 300 NG=1,NGTB

C ***     Calculate derivs of Xi wrt Nu (fibre) coords (DXINU)
          CALL XEXG(NBJ,ng,NJE,PG,VE,XE,XG,ERROR,*9999)
          CALL XGMG(1,0,NBJ(1),NJE,DXINU,GL,GU,RGX,XG,ERROR,*9999)
          IF(JTYP12.LE.1) THEN
            ETA1=XG(NJE+1,1)
            WRITE(OP_STRING,'(/'' Gauss point '',I2,'
     '        //''' Fibre angle='',F6.1,'' degrees wrt Xi(1) coord'')')
     '        NG,ETA1*180.0D0/PI
      	    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF(JTYP12.EQ.2) THEN
            ETA2=XG(NJE+1,1)
            WRITE(OP_STRING,'(/'' Gauss point '',I2,'
     '        //''' Fibre angle='',F6.1,'' degrees wrt Xi(2) coord'')')
     '        NG,ETA2*180.0D0/PI
      	    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF

C ***     undeformed r.c. coord derivs wrt Xi (XG)
          DO nu=1,4
            CALL XZ_DERIV(ITYP10(1),nu,XG,ZG)
          ENDDO
          DO nj=1,NJT
            DO nu=1,4
              XG(NJ,NU)=ZG(NJ,NU)
            ENDDO
          ENDDO
C ***     add displacements at time=0 - these are given by VE
          CALL XEXG(NBJ,ng,NJE,PG,VE,VE,VG,ERROR,*9999)
          DO J=1,3
            DO I=1,4
              VG(J,I)=XG(J,I)+VG(J,I)
            ENDDO
          ENDDO
C ***     add displacements at time=T - these are given by ZE
          CALL XEXG(NBJ,ng,NJE,PG,VE,ZE,ZG,ERROR,*9999)
          DO J=1,3
            DO I=1,4
              ZG(J,I)=XG(J,I)+ZG(J,I)
            ENDDO
          ENDDO
C ***     initial r.c. coord derivs wrt Nu (DXJNU)
          DO J=1,3
            DO I=1,2
              DXJNU(J,I)=VG(J,2)*DXINU(1,I)+VG(J,4)*DXINU(2,I)
            ENDDO
          ENDDO
C ***     derivs of displacement wrt Nu (DUJNU)
          DO J=1,3
            DO I=1,2
              DUJNU(J,I)=ZG(J,2)*DXINU(1,I)+ZG(J,4)*DXINU(2,I)
            ENDDO
          ENDDO
C ***     calculate undeformed metric (AXL) and deformed metric (AZL)
          DO J=1,2
            DO I=1,2
              AXL(I,J)=0.0D0
              AZL(I,J)=0.0D0
              DO K=1,3
                AXL(I,J)=AXL(I,J)+DXJNU(K,I)*DXJNU(K,J)
                AZL(I,J)=AZL(I,J)+DUJNU(K,I)*DUJNU(K,J)
              ENDDO
              EG(I,J)=0.5D0*(AZL(I,J)-AXL(I,J))
            ENDDO
          ENDDO
          CALL EVALUE(2,EG,EVAL,ERROR,*9999)
          CALL EVECTR(2,EG,EVAL(1),EVEC,ERROR,*9999)
          WRITE(OP_STRING,'('' EG(1,1)='',E12.4,'' EG(2,2)='',E12.4,
     '      '' EG(1,2)='',E12.4)') EG(1,1),EG(2,2),EG(1,2)
      	  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          PHI=DATAN2(EVEC(2),EVEC(1))
          WRITE(OP_STRING,'('' Princ. strains:'',2E12.4,''  PHI = '',
     '      F6.1,'' degrees'')') EVAL(1),EVAL(2),PHI*180.0D0/PI
      	  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(EVAL(1).GT.PRSTMAX) PRSTMAX=EVAL(1)
          IF(DABS(EVAL(1)).LT.PRSTMIN) PRSTMIN=DABS(EVAL(1))
 300    CONTINUE

 9999 CALL EXITS('OPST80')
      RETURN
      END

C KAT 2001-12-14
      SUBROUTINE OPTEXT(ERROR,*)

C#### Subroutine: OPTEXT
C###  Description:
C###    OPTEXT outputs text.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:ntsg00.cmn'
      INCLUDE 'cmiss$reference:text00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER notext

      CALL ENTERS('OPTEXT',*9999)

      DO notext=1,NTTEXT
        WRITE(OP_STRING,'(1X,I2,'') '',A)') notext,TEXT(notext)
      	CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDDO

      CALL EXITS('OPTEXT')
      RETURN
 9999 CALL ERRORS('OPTEXT',ERROR)
      CALL EXITS('OPTEXT')
      RETURN 1
      END


C18-nov-1999
      SUBROUTINE OPTIME(VARIABLE,NODAL_INFO,ERROR,*)

C#### Subroutine: OPTIME
C###  Description:
C###    OPTIME outputs time variables
C *** Created 10 August 1999 - David Nickerson

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:time_variable.cmn'
!     Parameter List
      INTEGER VARIABLE
      LOGICAL NODAL_INFO
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER var,START_VAR,END_VAR,node,elem
      CHARACTER BASIS_TYPES(4)*20

      DATA BASIS_TYPES /'Linear Lagrange',
     '                  'Quadratic Lagrange',
     '                  'Cubic Lagrange',
     '                  'Cubic Hermite'/

      CALL ENTERS('OPTIME',*9999)

      CALL ASSERT(TV_NUM_VARIABLES.GT.0,'>>no time variables defined',
     '  ERROR,*9999)

      IF(VARIABLE.GT.0) THEN
        START_VAR=VARIABLE
        END_VAR=VARIABLE
      ELSE
        START_VAR=1
        END_VAR=TV_NUM_VARIABLES
      ENDIF
      CALL ASSERT(START_VAR.LE.TV_NUM_VARIABLES.AND.
     '  END_VAR.LE.TV_NUM_VARIABLES,'Invalid variable number',
     '  ERROR,*9999)

      DO var=START_VAR,END_VAR

        WRITE(OP_STRING,'()')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C ***   Write out the variable's name
        WRITE(OP_STRING,'(''Variable name: '',A)') TV_NAMES(var)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C ***   Write out the basis type
        WRITE(OP_STRING,'(''Basis type: '',A)') 
     '    BASIS_TYPES(TV_BASIS_TYPES(var))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C ***   Write out the start and end times
        WRITE(OP_STRING,'(''Start time: '',D12.5)') TV_START_TIME(var)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''End time: '',D12.5)') TV_END_TIME(var)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C ***   Write out the elements
        WRITE(OP_STRING,'(''Number of elements: '',I5)') 
     '    TV_NUM_ELEMS(var)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        IF(NODAL_INFO) THEN
          DO elem=1,TV_NUM_ELEMS(var)
            WRITE(OP_STRING,'(''Element: '',I3)') elem
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            DO node=1,TV_NODES_PER_ELEMENT(TV_BASIS_TYPES(var))
              WRITE(OP_STRING,'(2X,''Node: '',I3,'' '',I3)') elem,node
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'(4X,''Time: '',D12.5)') 
     '          TV_NODE_TIMES(var,
     '          (elem-1)*TV_NODES_PER_ELEMENT(TV_BASIS_TYPES(var))+node)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'(4X,''Coordinate: '',D12.5)') 
     '          TV_NODE_COORDS(var,
     '          (elem-1)*TV_NODES_PER_ELEMENT(TV_BASIS_TYPES(var))+node)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO !node
          ENDDO !elem
        ENDIF!nodal_info

      ENDDO !var

      WRITE(OP_STRING,'()')
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      CALL EXITS('OPTIME')
      RETURN
 9999 CALL ERRORS('OPTIME',ERROR)
      CALL EXITS('OPTIME')
      RETURN 1
      END


      SUBROUTINE OPTRAN(ERROR,*)

C#### Subroutine: OPTRAN
C###  Description:
C###    OPTRAN outputs transformation parameters for Phigs 3D window.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:phig00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i

      CALL ENTERS('OPTRAN',*9999)

C *** Reference point, rotation angles, scaling & shift vector
C *** for world coord transformation
      WRITE(OP_STRING,'('' Rotation angles   :'',3E12.3)')
     '  (ANGLE(i),i=1,3)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Scale coefficients:'',3E12.3)')
     '  (SCALE(i),i=1,3)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Shift vector      :'',3E12.3)')
     '  (SHIFT(i),i=1,3)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,
     '  '('' Reference point   :'',3E12.3)') (VIEW_REF_PT_NEW(i),i=1,3)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(/'' Transformation matrix:'')')
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C *** Reference point, view plane normal & view up vector
C *** in world coords to orient viewing reference coords
      WRITE(OP_STRING,
     '  '(/'' View plane vector :'',3E12.3)') (VIEW_PLANE_NEW(i),i=1,3)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,
     '  '( '' View up vector    :'',3E12.3)') (VIEW_UP_NEW(i),i=1,3)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,
     '  '( '' Projection ref pt :'',3E12.3)') (PROJ_REF_PT_NEW(i),i=1,3)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(/'' Initial view orientation matrix:'')')
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(/'' Current view orientation matrix:'')')
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C *** Window, projection ref pt and viewplane, backplane & frontplane
C *** positions in viewing reference coords and viewport in norm proj
C *** coords for mapping to normalized projection coord system
      WRITE(OP_STRING,
     '  '(/'' Window limits     :'',6E12.3)') (WINDOW_NEW(i),i=1,6)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,
     '   '( '' Front plane distance = '',E12.3)') FRONT_PLANE_DIST_NEW
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,
     '  '( '' Back  plane distance = '',E12.3)') BACK_PLANE_DIST_NEW
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(/'' Mapping matrix:'')')
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      CALL EXITS('OPTRAN')
      RETURN
 9999 CALL ERRORS('OPTRAN',ERROR)
      CALL EXITS('OPTRAN')
      RETURN 1
      END


      SUBROUTINE OPVSAE(TYPE,ERROR,*)

C#### Subroutine: OPVSAE
C###  Description:
C###    OPVSAE outputs VSaero parameters.

      IMPLICIT NONE
!     Parameter List
      CHARACTER TYPE*(*),ERROR*(*)
!     Local Variables

      CALL ENTERS('OPVSAE',*9999)

      IF(TYPE(1:11).EQ.'BASIC_INPUT') THEN
        CALL OPVSA(ERROR,*9999)
      ELSE IF(TYPE(1:14).EQ.'PATCH_GEOMETRY') THEN
        CALL OPVSB(ERROR,*9999)
      ELSE IF(TYPE(1:10).EQ.'WAKE_INPUT') THEN
        CALL OPVSC(ERROR,*9999)
      ELSE IF(TYPE(1:19).EQ.'SURFACE_STREAMLINES') THEN
        CALL OPVSD(ERROR,*9999)
      ELSE IF(TYPE(1:20).EQ.'BOUNDARY_LAYER_INPUT') THEN
        CALL OPVSE(ERROR,*9999)
      ELSE IF(TYPE(1:22).EQ.'OFF_BODY_VELOCITY_SCAN') THEN
        CALL OPVSF(ERROR,*9999)
      ELSE IF(TYPE(1:20).EQ.'OFF_BODY_STREAMLINES') THEN
        CALL OPVSG(ERROR,*9999)
      ENDIF

      CALL EXITS('OPVSAE')
      RETURN
 9999 CALL ERRORS('OPVSAE',ERROR)
      CALL EXITS('OPVSAE')
      RETURN 1
      END


Module FE13
=========== 

      SUBROUTINE IPMESH2(IBT,IDO,INP,LD,LDR,NBJ,NCO,NDP,NE_TERM,
     '  NEELEM,NENP,NE_OLD,NE_TEMP,NJE,NKE,NKJ,N_SPACE,NP_INTERFACE,
     '  NPC_MIN,NPF,NPNE,NPNODE,nr,NRE,NVJE,NVJP,NW,NXI,PN_TEMP,
     '  BBM,CE,MIN_DIST,SE,WD,XA,XE,XP,ZD,ERROR,*)

C#### Subroutine: IPMESH2
C###  Description:
C###    IPMESH2 defines mesh parameters for fractal trees.

C****   CE(1,ne)   Generation number
C****   CE(2,ne)   Branch length
C****   CE(3,ne)   Branch diameter
C****   CE(4,ne)   Branch volume
C****   CE(5,ne)   Branch subtended volume
C****   CE(6,ne)   Flow
C****   CE(7,ne)   Diffusion Coefficient
C****   CE(8,ne)   Oxygen loss coefficient
C****   CE(9,ne)   Branch diameter ratio
C****   CE(10,ne)  Branch diameter,d (excluding sacs and ducts)
C****   CE(11,ne)  Proportion each branch is of total volume
C****   CE(12,ne)  Initial branch length
C****   CE(13,ne)  ****UNUSED*************
C****   CE(14,ne)  Initial branch diameter

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:mesh00.cmn'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LD(NDM),LDR(0:NDM),NBJ(NJM,NEM),NCO(NEM),
     '  NDP(NDM),NE_TERM(NE_R_M),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:60,0:NRM),NJE(NEM),NKE(NKM,NNM,NBFM,NEM),
     '  NKJ(NJM,NPM),NPC_MIN(NP_R_M),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),N_SPACE(0:NDM),NP_INTERFACE(0:NPM,0:3),
     '  NE_OLD(NORM),NE_TEMP(NORM),nr,NRE(NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),NW(NEM,2),
     '  NXI(-NIM:NIM,0:4,0:NEM),PN_TEMP(NE_R_M)
      REAL*8 BBM(NMM,NORM,NRM),CE(NMM,NEM),MIN_DIST(NE_R_M),
     '  SE(NSM,NBFM,NEM),WD(NJM,NDM),XA(NAM,NJM,NQM),XE(NSM,NJM),
     '  XP(NKM,NVM,NJM,NPM),ZD(NJM,NDM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER cnt,B_DIV(20),gen,G_LIM,IBEG,ICHAR,IEND,
     '  INDEX_IP1,INDEX_IP2,INDEX_IP3,INFO,INUM_BRANCH,ISEED1,M,N,nb,
     '  nbl,nba,nbh,nb2,nc,nd,nd1,nd2,nd3,NDT_ELEM,NDT2,ne,ne1,ne2,
     '  necond,ne_host,NGEN(5),nj,nk,nn,noelem,noelem2,nogen,nonode,
     '  NOQUES,np,NPOLD,np1,np2,np3,np4,npcond,nreg,NRB,nr_high,ns,nsp,
     '  NT_BNS,nu,N_ELM,N_GEN,N_TRM,NUM_NODES,NUM_SPLIT,N_DIRN,
     '  N_OLD_BRANCH(50),N_PARORD,NTMP,N_TEMP_BRANCH(50),N_TEMP_ORD,ord,
     '  PE_TEMP,RANDC_METHOD,TRAC
      REAL*8 A_LIM,ANG,ANG_XY,ANG_XY_MN,ANG_Y,ANG_Y_MN,BAS_FAC,BBM_VOL,
     '  BRANCH_D,BRANCH_L,B_DIAM,B_DIAM_SD,B_LEN,B_LEN_SD,COFM(3),
     '  DIAMETER,DIST,DIVLENGTH,FRAC(50),LEN,LENGTH,Lgen,LIMIT_ND,
     '  L_LIM,NORM_RAND_NUM,PANG1,PANG2,PXI,RANDOM_NUMBER,
     '  RATIO_ANG,sum,S_RATIO(99),temp,TEMPRAND,X(3),XP1(3),XI(3)
      CHARACTER CHAR1*1,CHAR2*2,CHAR3*17,OPT2(6)*17
      LOGICAL BRANCH,FAIL,FILEIP,PBRANCH,SS(2),TCHECK

      DATA OPT2(1)/'Right Upper Lobe '/,
     '     OPT2(2)/'Right Middle Lobe'/,
     '     OPT2(3)/'Right Lower Lobe'/,
     '     OPT2(4)/'Left Upper Lobe'/,
     '     OPT2(5)/'Left Lower Lobe'/,
     '     OPT2(6)/'Central Airways'/

      CALL ENTERS('IPMESH2',*9999)

      CALL ASSERT(NKM.GE.2,'>>NKM must be 2 or higher',ERROR,*9999)
      CALL ASSERT(NMM.GE.14,'>>NMM must be 14 or higher',ERROR,*9999)
      FILEIP=.FALSE.
      NOQUES=0
      IST_POINT=0
      DO gen=1,20
        NDIFF_DAUGHTER(gen)=0
        DO N=0,9
          NO_BRANCH_PATTERN(gen,N)=0
        ENDDO !N
      ENDDO !gen

      FORMAT='('' Enter model type [1]:'''//
     '  '/''   (1) Full lung model'''//
     '  '/''   (2) Conducting airway model'''//
     '  '/''   (3) Respiratory airway model'''//
     '  '/$,''    '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) IDATA(1)=1
      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) MODEL_TYPE=IDATA(1)
      IF((MODEL_TYPE.EQ.1).OR.(MODEL_TYPE.EQ.2))THEN !Conducting airways
        FORMAT='('' Enter conducting airway model type [3]:'''//
     '    '/''   (1) Connectivity delta model'''//
     '    '/''   (2) Volume-filling mesh'''//
     '    '/''   (3) Tube representation of symmetric airways'''//
     '    '/$,''    '',I1)'
        IDEFLT(1)=3
        IF(IOTYPE.EQ.3) IDATA(1)=3
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) CONDUCTING_TYPE=IDATA(1)
      ENDIF
      IF((MODEL_TYPE.EQ.1).OR.(MODEL_TYPE.EQ.3))THEN !Resp. airways
        FORMAT='('' Enter respiratory airway model type [1]:'''//
     '    '/''   (1) Multi-branching acinar mesh'''//
     '    '/''   (2) Black-Box model'''//
     '    '/$,''    '',I1)'
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=1
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) RESPIRATORY_TYPE=IDATA(1)
      ENDIF
      IF((MODEL_TYPE.EQ.1).OR.(MODEL_TYPE.EQ.2))THEN !Conducting airways
        IF(CONDUCTING_TYPE.EQ.1)THEN
C*** This code was written a long time ago, and has not been used for
C*** a number of years.  Not guaranteed to work.
          !growing a conducting airway model based on Horsfield's
          !delta model.  ie. generates conductivity based on delta
          !(the difference in Horsfield order of daughter branches)
          FORMAT='($,'' Enter basis function for airways [1]: '',I1)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=1
          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NBM,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) nba=IDATA(1)
          FORMAT='('' Enter starting point of model [2]:'''//
     '      '/''   (1) Trachea'''//
     '      '/''   (2) Beyond lobar bronchi'''//
     '      '/$,''    '',I1)'
          IDEFLT(1)=2
          IF(IOTYPE.EQ.3) IDATA(1)=2
          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3)IST_POINT=IDATA(1)
          IF(IST_POINT.EQ.1) THEN !Starting at Trachea
            N1_ORDER=27 !The order for the lobar bronchi
            FORMAT='($,'' Enter order of last branch [1]: '',I2)'
            IDEFLT(1)=1
            IF(IOTYPE.EQ.3) IDATA(1)=2
            CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,N1_ORDER-1,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) N2_ORDER=IDATA(1)
            NT_GEN=N1_ORDER-N2_ORDER+5
            B_ANGLE_Y(1)=0.0d0
            FORMAT='($,'' Enter angle to y-axis of 2nd branch '
     '        //'[70deg]: '',D12.4)'
            RDEFLT(1)=70.0d0
            IF(IOTYPE.EQ.3) RDATA(1)=70.0d0
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) B_ANGLE_Y(2)=RDATA(1)*PI/180.0d0
            DO N=1,9 !Enter information for trachea
              WRITE(CHAR1,'(I1)') N
              FORMAT='($,'' Enter length of branch '//CHAR1
     '          //': '',D12.4)'
              RDEFLT(1)=0.0d0
              IF(IOTYPE.EQ.3) RDATA(1)=0.0d0
              CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '          1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '          LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) B_LENGTH(N+90)=RDATA(1)
              FORMAT='($,'' Enter diameter : '',D12.4)'
              RDEFLT(1)=0.0d0
              IF(IOTYPE.EQ.3) RDATA(1)=0.0d0
              CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '          1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '          LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) B_DIAMETER(N+90)=RDATA(1)
              FORMAT='($,'' Enter angle for branch '//CHAR1
     '          //'[70deg]: '',D12.4)'
              RDEFLT(1)=70.0d0
              IF(IOTYPE.EQ.3) RDATA(1)=70.0d0
              CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '          1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '          LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) B_ANGLE_Y(N+90)=RDATA(1)*PI/180.0d0
            ENDDO !N
            FORMAT='($,'' Enter angle ratio for lobes [1.05]: '',D12.4)'
            RDEFLT(1)=1.05d0
            IF(IOTYPE.EQ.3) RDATA(1)=1.05d0
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) RATIO_ANGLE_2=RDATA(1)
            RATIO_ANGLE_1=RATIO_ANGLE_2
          ELSE IF(IST_POINT.EQ.2)THEN
            FORMAT='($,'' Enter order of first branch [1]: '',I2)'
            IDEFLT(1)=1
            IF(IOTYPE.EQ.3) IDATA(1)=1
            CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) N1_ORDER=IDATA(1)
            FORMAT='($,'' Enter order of last branch [1]: '',I2)'
            IDEFLT(1)=1
            IF(IOTYPE.EQ.3) IDATA(1)=1
            CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,N1_ORDER-1,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) N2_ORDER=IDATA(1)
            NT_GEN=N1_ORDER-N2_ORDER+1
            B_ANGLE_Y(1)=0.0d0
            FORMAT='($,'' Enter angle to y-axis of 2nd branch '
     '        //'[70deg]: '',D12.4)'
            RDEFLT(1)=70.0d0
            IF(IOTYPE.EQ.3) RDATA(1)=70.0d0
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) B_ANGLE_Y(2)=RDATA(1)*PI/180.0d0
            FORMAT='($,'' Enter angle ratio for gens 1-3 [2.0]: '','
     '        //'D12.4)'
            RDEFLT(1)=2.0d0
            IF(IOTYPE.EQ.3) RDATA(1)=2.0d0
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) RATIO_ANGLE_1=RDATA(1)
            FORMAT='($,'' Enter angle ratio for gens 4.. [1.2]: '','
     '        //'D12.4)'
            RDEFLT(1)=1.20d0
            IF(IOTYPE.EQ.3) RDATA(1)=1.20d0
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) RATIO_ANGLE_2=RDATA(1)
          ENDIF !IST_POINT
        ELSE IF(CONDUCTING_TYPE.EQ.2)THEN
          !volume-filling mesh
          IF(CONDUCTING_TYPE.EQ.2)THEN
            FORMAT='('' Enter volume-filling model type [2]:'''//
     '        '/''   (1) Full lung model (5 lobes)'''//
     '        '/''   (2) Single lobe model'''//
     '        '/$,''    '',I1)'
            IDEFLT(1)=2
            IF(IOTYPE.EQ.3) IDATA(1)=2
            CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3)IST_POINT=IDATA(1)
            IF(IST_POINT.EQ.1) THEN !a full lung model (all 5 lobes)
              nr_high=6
              FORMAT='('' Enter method of model entry [2]:'''//
     '          '/''   (1) Defined by .ipmesh file'''//
     '          '/''   (2) Read in from .ipelem and .ipnode files'''//
     '          '/$,''    '',I1)'
              IDEFLT(1)=2
              IF(IOTYPE.EQ.3) IDATA(1)=2
              CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '          1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '          LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3)INDEX_IP2=IDATA(1)
              FORMAT='('' Enter method for solution [2]:'''//
     '          '/''   (1) Solve whole model simultaneously'''//
     '          '/''   (2) Solve lobes seperately for inspiration'''//
     '          '/$,''    '',I1)'
              IDEFLT(1)=2
              IF(IOTYPE.EQ.3) IDATA(1)=2
              CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '          1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '          LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3)SOLVE_TYPE=IDATA(1)
            ENDIF
            IF(INDEX_IP2.EQ.2)THEN !read in
              CALL ASSERT(NPNODE(0,nr).GT.0,'>>Define nodes first',
     '          ERROR,*9999)
              FORMAT='($,'' Enter element # of trachea: '',I2)'
              IF(IOTYPE.EQ.3) IDATA(1)=1
              CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '          1,99,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     '          *9999)
              IF(IOTYPE.NE.3) TRAC=IDATA(1)
            ENDIF
            DO nreg=1,nr_high
              WRITE(CHAR3,'(A)') OPT2(nreg)
              IF(INDEX_IP2.EQ.1)THEN !grown
                FORMAT='($,'' Enter region # of '//CHAR3(1:17)//': '','
     '            //'I2)'
                IF(IOTYPE.EQ.3) IDATA(1)=1
                CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '            1,99,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     '            *9999)
                IF(IOTYPE.NE.3) INDEX_REG(nreg)=IDATA(1)
              ELSE IF(INDEX_IP2.EQ.2)THEN !read in
                WRITE(OP_STRING,'($,'' '//CHAR3(1:17)//': '')') 
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
              IF(nreg.NE.6)THEN
                FORMAT='($,'' Enter element # of lobe bronchus '
     '            //'[1]:'',I2)'
                IDEFLT(1)=1
                IF(IOTYPE.EQ.3) IDATA(1)=1
                CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '            1,NE_R_M,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '            ERROR,*9999)
                IF(IOTYPE.NE.3) BRN_STORE(nreg)=IDATA(1)
              ENDIF
            ENDDO !nreg
            IF(INDEX_IP2.EQ.2) nr_high=1
          ELSE IF(IST_POINT.EQ.2)THEN !single lobe model
            INDEX_IP2=1
            FORMAT='('' Enter lobe to be modelled [1]:'''//
     '        '/''   (1) Right upper lobe'''//
     '        '/''   (2) Right middle lobe'''//
     '        '/''   (3) Right lower lobe'''//
     '        '/''   (4) Left upper lobe'''//
     '        '/''   (5) Left lower lobe'''//
     '        '/$,''    '',I1)'
            IDEFLT(1)=1
            IF(IOTYPE.EQ.3) IDATA(1)=1
            CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,5,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3)INDEX_REG(1)=IDATA(1)
            INDEX_IP2=1
            FORMAT='('' Enter option for element/node #s [2]:'''//
     '        '/''   (1) Specify start global node and element #s'''//
     '        '/''   (2) Default from highest host #s'''//
     '        '/$,''    '',I1)'
            IDEFLT(1)=2
            IF(IOTYPE.EQ.3) IDATA(1)=2
            CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3)INDEX_IP3=IDATA(1)
            IF(INDEX_IP3.EQ.1)THEN !prompt for start node and elem #s
              FORMAT='($,'' Start global node # for mesh: '',I5)'
              IF(IOTYPE.EQ.3) IDATA(1)=1
              CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '          1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '          NPNODE(0,0),100000,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '          INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) NPT(0)=IDATA(1)
              FORMAT='($,'' Start global element # for mesh: '',I5)'
              IF(IOTYPE.EQ.3) IDATA(1)=1
              CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '          1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '          NEELEM(0,0),100000,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '          INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) NET(0)=IDATA(1)
            ENDIF !INDEX_IP3
            nr_high=2
            WRITE(CHAR3,'(A)') OPT2(INDEX_REG(1))        
            FORMAT='($,'' Enter region # of '//CHAR3(1:17)//': '',I2)'
            IF(IOTYPE.EQ.3) IDATA(1)=1
            CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) INDEX_REG(1)=IDATA(1)
            FORMAT='($,'' Enter element # of lobe bronchus [1]:'',I2)'
            IDEFLT(1)=1
            IF(IOTYPE.EQ.3) IDATA(1)=1
            CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '        1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NE_R_M,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) BRN_STORE(1)=IDATA(1)
            WRITE(CHAR3,'(A)') OPT2(6) !region num for central airways
            FORMAT='($,'' Enter region # of '//CHAR3(1:17)//': '',I2)'
            IF(IOTYPE.EQ.3) IDATA(1)=1
            CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) INDEX_REG(2)=IDATA(1)
          ENDIF !IST_POINT
          IF(INDEX_IP2.EQ.1)THEN
            FORMAT='($,'' Enter basis function for lobes [1]: '',I1)'
            IDEFLT(1)=1
            IF(IOTYPE.EQ.3) IDATA(1)=1
            CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NBM,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) nbl=IDATA(1)
            FORMAT='($,'' Enter basis function for airways [1]: '',I1)'
            IDEFLT(1)=1
            IF(IOTYPE.EQ.3) IDATA(1)=1
            CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NBM,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) nba=IDATA(1)
            FORMAT='($,'' Enter the angle limit [90deg]: '',D12.4)'
            RDEFLT(1)=90.0d0
            IF(IOTYPE.EQ.3) RDATA(1)=1.0d0
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) A_LIM=RDATA(1)*PI/180.0d0
            FORMAT='($,'' Enter the length limit [0.1]: '',D12.4)'
            RDEFLT(1)=0.1d0
            IF(IOTYPE.EQ.3) RDATA(1)=1.0d0
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) L_LIM=RDATA(1)
            FORMAT='($,'' Enter highest gen for reass [1]: '',I4)'
            IDEFLT(1)=1
            IF(IOTYPE.EQ.3) IDATA(1)=1
            CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) G_LIM=IDATA(1)
            FORMAT='('' Enter method for random points [2]:'''//
     '        '/''   (1) Read in'''//
     '        '/''   (2) Generated by RNS'''//
     '        '/$,''    '',I1)'
            IDEFLT(1)=2
            IF(IOTYPE.EQ.3) IDATA(1)=2
            CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) RANDC_METHOD=IDATA(1)
            IF(RANDC_METHOD.EQ.1)THEN
              CALL ASSERT(NDT.GE.1,'>>Read data pnts first',ERROR,*9999)
            ELSE IF(RANDC_METHOD.EQ.2)THEN
              FORMAT='($,'' Enter # of pnts per element [5000]: '',I4)'
              IDEFLT(1)=5000
              IF(IOTYPE.EQ.3) IDATA(1)=5000
              CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '          1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NDM,
     '          LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) NDT_ELEM=IDATA(1)
              FORMAT='($,'' Enter the spacing limit [0.05]: '',D12.4)'
              RDEFLT(1)=0.05d0
              IF(IOTYPE.EQ.3) RDATA(1)=0.05d0
              CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '          1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,
     '          IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     '          *9999)
              IF(IOTYPE.NE.3) LIMIT_ND=RDATA(1)
              FORMAT='($,'' Enter RNS (0-9) for dimensions [0]: '',I1)'
              IDEFLT(1)=0
              IF(IOTYPE.EQ.3) IDATA(1)=0
              CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '          1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,9,
     '          LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) ISEED_GEOM=IDATA(1)
            ENDIF !RANDC_METHOD
            DO gen=4,20
              WRITE(CHAR2,'(I2)') gen        
              FORMAT='($,'' Enter fractal for generation '//
     '          CHAR2(1:2)//' [0.5]: '',D12.4)'
              RDEFLT(1)=0.5d0
              IF(IOTYPE.EQ.3) RDATA(1)=0.5d0
              CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '          1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '          LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) FRAC(gen)=RDATA(1)
            ENDDO !gen
            FORMAT='($,'' Enter fractal for generation 21... [0.5]: '','
     '        //'D12.4)'
            RDEFLT(1)=0.5d0
            IF(IOTYPE.EQ.3) RDATA(1)=0.5d0 
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) FRAC(21)=RDATA(1)
            DO gen=22,35
              FRAC(gen)=FRAC(21)
            ENDDO !gen
          ELSE IF(INDEX_IP2.EQ.2)THEN
            FORMAT='($,'' Enter basis function for airways [1]: '',I1)'
            IDEFLT(1)=1
            IF(IOTYPE.EQ.3) IDATA(1)=1
            CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NBM,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) nba=IDATA(1)
          ENDIF
        ELSE IF(CONDUCTING_TYPE.EQ.3)THEN !tube for symmetric airways
          FORMAT='($,'' Enter basis function for airways [1]: '',I1)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=1
          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NBM,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) nba=IDATA(1)
          FORMAT='($,'' Enter number of tube generations [1]: '',I2)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=1
          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NTT_GEN=IDATA(1)
          DO gen=1,NTT_GEN
            WRITE(CHAR2,'(I2)') gen
            FORMAT='($,'' Enter length of tube branch '//CHAR2
     '        //' [1.0]: '',D12.4)'
            RDEFLT(1)=1.0d0
            IF(IOTYPE.EQ.3) RDATA(1)=1.0d0
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '        1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) B_LENGTH(gen)=RDATA(1)
            FORMAT='($,'' Enter tube diameter [0.1]: '',D12.4)'
            RDEFLT(1)=0.1d0
            IF(IOTYPE.EQ.3) RDATA(1)=0.1d0
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '        1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) B_DIAMETER(gen)=RDATA(1)
            FORMAT='($,'' Enter number of divisions [1]: '',I4)'
            IDEFLT(1)=1
            IF(IOTYPE.EQ.3) IDATA(1)=1
            CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NE_R_M,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) B_DIV(gen)=IDATA(1)
          ENDDO !gen
        ENDIF !CONDUCTING_TYPE
      ENDIF !MODEL_TYPE
      IF((MODEL_TYPE.EQ.1).OR.(MODEL_TYPE.EQ.3))THEN !Resp. airways
        IF(MODEL_TYPE.EQ.3)THEN !isolated acinar model
          FORMAT='($,'' Enter basis function for airways [1]: '',I1)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=1
          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NBM,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) nba=IDATA(1)
        ENDIF
        IF(RESPIRATORY_TYPE.EQ.1)THEN
          FORMAT='($,'' Enter number of acinar gens [1]: '',I2)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=1
          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NT_GEN=IDATA(1)
          FORMAT='('' Enter method of data entry [2]:'''//
     '      '/''   (1) Enter data for first branch only'''//
     '      '/''   (2) Enter lengths, diameters and angles'''//
     '      '/$,''    '',I1)'
          IDEFLT(1)=2
          IF(IOTYPE.EQ.3) IDATA(1)=2
          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) INDEX_IP1=IDATA(1)
          IF(INDEX_IP1.EQ.2) THEN !Lengths and diameters from input
            DO nogen=1,NT_GEN
              IF(CONDUCTING_TYPE.EQ.3)THEN
                gen=NTT_GEN+nogen
              ELSE IF(MODEL_TYPE.EQ.3)THEN
                gen=nogen
              ENDIF 
              WRITE(CHAR2,'(I2)') gen
              FORMAT='($,'' Enter length of branch '//CHAR2
     '          //' [1.0]: '',D12.4)'
              RDEFLT(1)=1.0d0
              IF(IOTYPE.EQ.3) RDATA(1)=1.0d0
              CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '          1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '          LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) B_LENGTH(gen)=RDATA(1)
              FORMAT='($,'' Enter diameter [0.1]: '',D12.4)'
              RDEFLT(1)=0.1d0
              IF(IOTYPE.EQ.3) RDATA(1)=0.1d0
              CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '          1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '          LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) B_DIAMETER(gen)=RDATA(1)
              FORMAT='($,'' Enter ratio (s/S) [1.0]: '',D12.4)'
              RDEFLT(1)=1.0d0
              IF(IOTYPE.EQ.3) RDATA(1)=1.0d0
              CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '          1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '          LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) S_RATIO(gen)=RDATA(1)
              FORMAT='($,'' Enter angle [70.0]: '',D12.4)'
              RDEFLT(1)=70.0d0
              IF(IOTYPE.EQ.3) RDATA(1)=70.0d0
              CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '          1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '          LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) W_ANGLE_Y(gen)=RDATA(1)*PI/180.d0
            ENDDO !nogen (gen)
            FORMAT='($,'' Enter scale factor base/apex [1.0]: '',D12.4)'
            RDEFLT(1)=1.0d0
            IF(IOTYPE.EQ.3) RDATA(1)=1.d0
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '        1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) BAS_FAC=RDATA(1)
          ELSE IF(INDEX_IP1.EQ.1) THEN !Input 1st values, then scale
            FORMAT='($,'' Enter length of branch 1 [1.0]: '',D12.4)'
            RDEFLT(1)=1.0d0
            IF(IOTYPE.EQ.3) RDATA(1)=1.0d0
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '        1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) B_LENGTH(1)=RDATA(1)
            FORMAT='($,'' Enter diameter [0.1]: '',D12.4)'
            RDEFLT(1)=0.1d0
            IF(IOTYPE.EQ.3) RDATA(1)=0.1d0
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '        1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) B_DIAMETER(1)=RDATA(1)
            FORMAT='($,'' Enter angle 1st gen [0.0]: '',D12.4)'
            RDEFLT(1)=0.0d0
            IF(IOTYPE.EQ.3) RDATA(1)=0.0d0
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '        1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) W_ANGLE_Y(1)=RDATA(1)*PI/180.0d0
            FORMAT='($,'' Enter angle 2nd gen [70.0]: '',D12.4)'
            RDEFLT(1)=70.0d0
            IF(IOTYPE.EQ.3) RDATA(1)=70.0d0
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '        1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) W_ANGLE_Y(2)=RDATA(1)*PI/180.0d0
            FORMAT='($,'' Enter length ratio 2nd gen [1.5]: '',D12.4)'
            RDEFLT(1)=1.5d0
            IF(IOTYPE.EQ.3) RDATA(1)=1.5d0
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '        1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) RATIO_LENGTH=RDATA(1)
            FORMAT='($,'' Enter diameter ratio 2nd gen [1.4]: '',D12.4)'
            RDEFLT(1)=1.4d0
            IF(IOTYPE.EQ.3) RDATA(1)=1.4d0
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '        1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) RATIO_DIAMETER=RDATA(1)
            FORMAT='($,'' Enter length ratio gen 3..[1.5]: '',D12.4)'
            RDEFLT(1)=1.5d0
            IF(IOTYPE.EQ.3) RDATA(1)=1.5d0
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '        1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) RATIO2_LENGTH=RDATA(1)
            FORMAT='($,'' Enter diameter ratio gen 3..[1.4]: '',D12.4)'
            RDEFLT(1)=1.4d0
            IF(IOTYPE.EQ.3) RDATA(1)=1.4d0
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '        1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) RATIO2_DIAMETER=RDATA(1)
            FORMAT='($,'' Enter angle ratio [1.2]: '',D12.4)'
            RDEFLT(1)=1.2d0
            IF(IOTYPE.EQ.3) RDATA(1)=1.2d0
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '        1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) RATIO_ANG=RDATA(1)
            B_LENGTH(2)=B_LENGTH(1)/RATIO_LENGTH
            B_DIAMETER(2)=B_DIAMETER(1)/RATIO_DIAMETER
            DO gen=3,NT_GEN
              B_LENGTH(gen)=B_LENGTH(gen-1)/RATIO2_LENGTH
              B_DIAMETER(gen)=B_DIAMETER(gen-1)/RATIO2_DIAMETER
              W_ANGLE_Y(gen)=W_ANGLE_Y(gen-1)/RATIO_ANG
            ENDDO !gen
          ENDIF !INDEX_IP1
          DO nogen=1,NT_GEN !Loop over each generation
            IF(CONDUCTING_TYPE.EQ.3)THEN
              gen=nogen+NTT_GEN
            ELSE IF(MODEL_TYPE.EQ.3)THEN
              gen=nogen
            ENDIF
            WRITE(CHAR2,'(I2)') gen
            FORMAT='($,'' Enter mean branching ratio (0-3) for gen '
     '        //CHAR2//' [2.0]: '',D12.4)'
            RDEFLT(1)=2.0d0
            IF(IOTYPE.EQ.3) RDATA(1)=2.0d0
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '        1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,3.0d0,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) MNB(gen)=RDATA(1)
            FORMAT='($,'' Enter std dev of branching ratio for gen '
     '        //CHAR2//' [1.0]: '',D12.4)'
            RDEFLT(1)=1.0d0
            IF(IOTYPE.EQ.3) RDATA(1)=1.0d0
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '        1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) SDNB(gen)=RDATA(1)
          ENDDO !nogen (gen)
        ELSE IF(RESPIRATORY_TYPE.EQ.2)THEN
          FORMAT='('' Enter acinar volume distribution [1]:'''//
     '      '/''   (1) Same volume throughout lung model'''//
     '      '/''   (2) Distributed by some function...'''//
     '      '/$,''    '',I1)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=1
          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) VOLUME_TYPE=IDATA(1)
          IF(VOLUME_TYPE.EQ.1)THEN
            FORMAT='($,'' Enter the BBM volume [181.4mm^3]: '',D12.4)'
            RDEFLT(1)=181.4d0
            IF(IOTYPE.EQ.3) RDATA(1)=181.4d0
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3)BBM_VOL=RDATA(1)
          ELSE IF(VOLUME_TYPE.EQ.2)THEN
            WRITE(OP_STRING,'('' Sorry, not implemented yet'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF !VOLUME_TYPE
        ENDIF !RESPIRATORY_TYPE
      ENDIF !MODEL_TYPE
      FORMAT='('' Enter ventilation type [1]:'''//
     '  '/''   (1) Uniform ventilation'''//
     '  '/''   (2) Hydrodynamic calculation'''//
     '  '/$,''    '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) IDATA(1)=1
      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,ADATA,
     '  ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,LDATA,LDEFLT,RDATA,
     '  RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) INDEX_VENTILATION=IDATA(1)
      FORMAT='('' Enter compliance relationship [1]:'''//
     '  '/''   (1) Cube-root of volume change'''//
     '  '/''   (2) Pressure volume relation'''//
     '  '/$,''    '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) IDATA(1)=1
      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,ADATA,
     '  ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,LDATA,LDEFLT,RDATA,
     '  RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) INDEX_COMPLIANCE=IDATA(1)
      IF((MODEL_TYPE.EQ.1).OR.(MODEL_TYPE.EQ.2))THEN !Conducting airways
        IF(CONDUCTING_TYPE.EQ.1)THEN
          IF(NJT.EQ.3) THEN !3D tree
            B_ANGLE_XY(1)=0.0d0
            FORMAT='($,'' Enter XYangle 2nd branch [90.0]: '',D12.4)'
            RDEFLT(1)=90.0d0
            IF(IOTYPE.EQ.3) RDATA(1)=90.0d0
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) B_ANGLE_XY(2)=RDATA(1)*PI/180.0d0
          ENDIF !NJT
          IF(JTYP14.EQ.3) THEN !Regular fractal tree
            B_ANGLE_CV=0.0d0
            FORMAT='($,'' Enter COV for angle [0.1]: '',D12.4)'
            RDEFLT(1)=0.1d0
            IF(IOTYPE.EQ.3) RDATA(1)=0.1d0
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) B_ANGLE_CV=RDATA(1)
          ELSE IF(JTYP14.EQ.4) THEN !Stochastic fractal tree
            FORMAT='('' Enter type of probability dist moment [1]:'''//
     '        '/''   (1) Moments scaled by fractal scaling law'''//
     '        '/''   (2) Moments entered for each generation'''//
     '        '/$,''    '',I1)'
            IDEFLT(1)=1
            IF(IOTYPE.EQ.3) IDATA(1)=1
            CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,LDATA,
     '        LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) INDEX_MOMENTS=IDATA(1)
            IF(INDEX_MOMENTS.EQ.1) THEN !Moment scaled by fract scal law
              FORMAT='($,'' Enter COV for length [0.1]: '',D12.4)'
              RDEFLT(1)=0.1d0
              IF(IOTYPE.EQ.3) RDATA(1)=0.1d0
              CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '          1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '          LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) B_LENGTH_CV=RDATA(1)
              FORMAT='($,'' Enter COV for diameter [0.1]: '',D12.4)'
              RDEFLT(1)=0.1d0
              IF(IOTYPE.EQ.3) RDATA(1)=0.1d0
              CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '          1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '          LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) B_DIAMETER_CV=RDATA(1)
            ENDIF !INDEX_MOMENTS
            IF(IST_POINT.EQ.1) THEN !starting at Trachea
              DO N=91,99 !Calculate standard deviations
                B_LENGTH_SD(N)=B_LENGTH_CV*B_LENGTH(N)
                B_DIAMETER_SD(N)=B_DIAMETER_CV*B_DIAMETER(N)
              ENDDO !N
              B_LENGTH_SD(1)=B_LENGTH_SD(N)
              B_DIAMETER_SD(1)=B_DIAMETER_SD(N)
            ENDIF !IST_POINT
            DO ord=N1_ORDER,N2_ORDER,-1 !Calc std devs
              B_LENGTH_SD(ord)=B_LENGTH_CV*B_LENGTH(ord)
              B_DIAMETER_SD(ord)=B_DIAMETER_CV*B_DIAMETER(ord)
            ENDDO !ord
          ENDIF !INDEX_BRANCH
          FORMAT='($,'' Enter RNS (0-9) for # branches [0]: '',I1)'
          IDEFLT(1)=0
          IF(IOTYPE.EQ.3) IDATA(1)=0
          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,9,LDATA,
     '      LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) ISEED_BRANCH=IDATA(1)
          FORMAT='($,'' Enter RNS (0-9) for dimensions [0]: '',I1)'
          IDEFLT(1)=0
          IF(IOTYPE.EQ.3) IDATA(1)=0
          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,9,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) ISEED_GEOM=IDATA(1)
          FORMAT='($,'' Enter RNS (0-9) for angles [0]: '',I1)'
          IDEFLT(1)=0
          IF(IOTYPE.EQ.3) IDATA(1)=0
          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,9,LDATA,
     '      LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) ISEED1=IDATA(1)
          B_ANGLE_Y(N1_ORDER)=B_ANGLE_Y(2)
          IF(NJT.EQ.3) B_ANGLE_XY(N1_ORDER)=B_ANGLE_XY(2)
          B_ANGLE_SD(N1_ORDER)=B_ANGLE_CV*B_ANGLE_Y(2)
          IF(IST_POINT.EQ.1)THEN !Starting at Trachea
            IF(NJT.EQ.3) B_ANGLE_XY(91)=B_ANGLE_XY(91)
            B_ANGLE_SD(91)=B_ANGLE_CV*B_ANGLE_Y(91)
            DO ord=92,99 !Calculate std dev & initial angles
              B_ANGLE_SD(ord)=B_ANGLE_CV*B_ANGLE_Y(ord)
              IF(NJT.EQ.3) THEN !3D tree
                B_ANGLE_XY(ord)=B_ANGLE_XY(ord-1)/RATIO_ANGLE_2
              ENDIF !NJT
            ENDDO !ord
          ENDIF !IST_POINT
          DO ord=N1_ORDER-1,N2_ORDER,-1
            IF(ord.GT.N1_ORDER-3) THEN
              B_ANGLE_Y(ord)=B_ANGLE_Y(ord+1)/RATIO_ANGLE_1
              B_ANGLE_SD(ord)=B_ANGLE_CV*B_ANGLE_Y(ord)
              IF(NJT.EQ.3) B_ANGLE_XY(ord)=B_ANGLE_XY(ord+1)/
     '          RATIO_ANGLE_1
            ELSE IF(ord.LE.N1_ORDER-3) THEN !Calc std dev & angles
              B_ANGLE_Y(ord)=B_ANGLE_Y(ord+1)/RATIO_ANGLE_2
              B_ANGLE_SD(ord)=B_ANGLE_CV*B_ANGLE_Y(ord)
              IF(NJT.EQ.3) B_ANGLE_XY(ord)=B_ANGLE_XY(ord+1)/
     '          RATIO_ANGLE_2
            ENDIF !ord
          ENDDO !ord
        ENDIF
      ENDIF 
      IF((MODEL_TYPE.EQ.1).OR.(MODEL_TYPE.EQ.3))THEN
        IF(RESPIRATORY_TYPE.EQ.1)THEN
          IF(NJT.EQ.3) THEN !3D tree
            B_ANGLE_XY(1)=0.0d0
            FORMAT='($,'' Enter XYangle 2nd bnch [90deg]: '',D12.4)'
            RDEFLT(1)=90.0d0
            IF(IOTYPE.EQ.3) RDATA(1)=90.0d0
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) B_ANGLE_XY(2)=RDATA(1)*PI/180.0d0
          ENDIF !NJT
          IF(JTYP14.EQ.3) THEN !Regular fractal tree
            B_ANGLE_CV=0.0d0
            FORMAT='($,'' Enter COV for angle [0.1]: '',D12.4)'
            RDEFLT(1)=0.1d0
            IF(IOTYPE.EQ.3) RDATA(1)=0.1d0
            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) B_ANGLE_CV=RDATA(1)
          ELSE IF(JTYP14.EQ.4) THEN !Stochastic fractal tree
            FORMAT='('' Enter type of probability dist moment [1]:'''//
     '        '/''   (1) Moments scaled by fractal scaling law'''//
     '        '/''   (2) Moments entered for each generation'''//
     '        '/$,''    '',I1)'
            IDEFLT(1)=1
            IF(IOTYPE.EQ.3) IDATA(1)=1
            CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) INDEX_MOMENTS=IDATA(1)
            IF(INDEX_MOMENTS.EQ.1) THEN !Moment scaled by fract scal law
              FORMAT='($,'' Enter COV for length [0.1]: '',D12.4)'
              RDEFLT(1)=0.1d0
              IF(IOTYPE.EQ.3) RDATA(1)=0.1d0
              CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '          1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '          LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) B_LENGTH_CV=RDATA(1)
              FORMAT='($,'' Enter COV for diameter [0.1]: '',D12.4)'
              RDEFLT(1)=0.1d0
              IF(IOTYPE.EQ.3) RDATA(1)=0.1d0
              CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '          1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '          LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) B_DIAMETER_CV=RDATA(1)
            ENDIF
          ENDIF
          DO nogen=1,NT_GEN
            IF(CONDUCTING_TYPE.EQ.3)THEN
              gen=nogen+NTT_GEN
            ELSE IF(MODEL_TYPE.EQ.3)THEN
              gen=nogen
            ENDIF 
            IF(INDEX_IP1.EQ.2) THEN
              B_LENGTH(gen)=B_LENGTH(gen)*BAS_FAC
              B_DIAMETER(gen)=B_DIAMETER(gen)*BAS_FAC
            ENDIF !INDEX_IP1
            B_LENGTH_SD(gen)=B_LENGTH_CV*B_LENGTH(gen)
            B_DIAMETER_SD(gen)=B_DIAMETER_CV*B_DIAMETER(gen)
          ENDDO !nogen (gen)
          FORMAT='($,'' Enter RNS (0-9) for # branches [0]: '',I1)'
          IDEFLT(1)=0
          IF(IOTYPE.EQ.3) IDATA(1)=0
          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,9,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) ISEED_BRANCH=IDATA(1)
          FORMAT='($,'' Enter RNS (0-9) for dimensions [0]: '',I1)'
          IDEFLT(1)=0
          IF(IOTYPE.EQ.3) IDATA(1)=0
          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,9,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) ISEED_GEOM=IDATA(1)
          FORMAT='($,'' Enter RNS (0-9) for angle [0]: '',I1)'
          IDEFLT(1)=0
          IF(IOTYPE.EQ.3) IDATA(1)=0
          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,9,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) ISEED1=IDATA(1)
          W_ANGLE_SD(1)=0.0d0
          DO nogen=2,NT_GEN !Calc std devs and initial angles
            IF(CONDUCTING_TYPE.EQ.3)THEN
              gen=nogen+NTT_GEN
            ELSE IF(MODEL_TYPE.EQ.3)THEN
              gen=nogen
            ENDIF 
            W_ANGLE_SD(gen)=B_ANGLE_CV*W_ANGLE_Y(gen)
            IF(NJT.EQ.3) B_ANGLE_XY(gen)=B_ANGLE_XY(gen-1)
          ENDDO !nogen (gen)
        ENDIF !CONDUCTING_TYPE
      ENDIF !MODEL_TYPE
      IF((MODEL_TYPE.EQ.1).OR.(MODEL_TYPE.EQ.2))THEN
        IF(CONDUCTING_TYPE.EQ.1)THEN
C*** MHT 10-09-97 
C*** I have tried to decipher and comment the following section as best 
C*** as possible, but it is not guaranteed to be right.
C GMH 8/1/97 Update cmgui link
          CALL NODE_CHANGE(1,.FALSE.,ERROR,*9999)
          XP(1,1,1,1)=0.0d0 !node 1 x is at centre
          XP(1,1,2,1)=0.0d0 !node 1 y is at bottom
          XP(1,1,3,1)=0.0d0 !node 1 z is at bottom
          XP(1,1,4,1)=0.0d0 !node 1 branch angle is zero
          XP(2,1,4,1)=0.0d0 !node 1 branch angle (3D) is zero
          np=1
          ne=0
          NET(nr)=ne
          NPT(nr)=np
          nb=nba
          NBT=1
          NFBASE(1,1)=1
          NBASEF(1,0)=1
          NBASEF(1,1)=1
          NBFT=1
          NE_OLD(1)=np !NE_OLD(n) stores previous generation node #s
          NE_OLD(1)=N1_ORDER
          N_OLD_BRANCH(1)=1
          NUM_NODES=1
          NUM_SPLIT=1
          DO gen=1,NT_GEN !for each generation in the conducting airways
            NT_BNS=0
            DO M=1,NUM_NODES
              DO N=1,NUM_SPLIT
                IF(gen.LE.5) THEN !only needed for 1st 5 gens
                  IF(gen.EQ.1) THEN !Set up pattern for 1st 9 gens
                    IF(IST_POINT.EQ.1) THEN !starts at Trachea
                      N_TEMP_ORD=31 !order of the trachea
                      TCHECK=.TRUE. !Checks if branch in Trachea section
                      PBRANCH=.TRUE. !Checks if parent has branch # GT 0
                    ELSE !starts at a lobar bronchus
                      N_TEMP_ORD=N1_ORDER !highest order in lobe
                      TCHECK=.FALSE. !not in trachea section
                      PBRANCH=.FALSE.
                    ENDIF !IST_POINT
                    NTMP=2-IST_POINT !1 if start Trachea else 0
                  ELSE IF(N_OLD_BRANCH(M).EQ.1) THEN
                    IF(N.EQ.1) THEN
                      NTMP=5
                      N_TEMP_ORD=30
                    ELSE IF(N.EQ.2) THEN
                      NTMP=2
                      N_TEMP_ORD=29
                    ENDIF !N
                    TCHECK=.TRUE.
                    PBRANCH=.TRUE.
                  ELSE IF(N_OLD_BRANCH(M).EQ.2) THEN
                    IF(N.EQ.1) THEN
                      NTMP=4
                      N_TEMP_ORD=28
                    ELSE IF(N.EQ.2) THEN
                      NTMP=3
                      N_TEMP_ORD=27
                    ENDIF !N
                    TCHECK=.TRUE.
                    PBRANCH=.TRUE.
                  ELSE IF(N_OLD_BRANCH(M).EQ.3) THEN
                    IF(N.EQ.1) THEN
                      NTMP=0
                      N_TEMP_ORD=26
                    ELSE IF(N.EQ.2) THEN
                      NTMP=0
                      N_TEMP_ORD=23
                    ENDIF !N
                    TCHECK=.FALSE.
                    PBRANCH=.TRUE.
                  ELSE IF(N_OLD_BRANCH(M).EQ.4) THEN
                    IF(N.EQ.1) THEN
                      NTMP=0
                      N_TEMP_ORD=27
                    ELSE IF(N.EQ.2) THEN
                      NTMP=0
                      N_TEMP_ORD=24
                    ENDIF !N
                    TCHECK=.FALSE.
                    PBRANCH=.TRUE.
                  ELSE IF(N_OLD_BRANCH(M).EQ.5) THEN
                    IF(N.EQ.1) THEN
                      NTMP=6
                      N_TEMP_ORD=27
                    ELSE IF(N.EQ.2) THEN
                      NTMP=7
                      N_TEMP_ORD=29
                    ENDIF !N
                    TCHECK=.TRUE.
                    PBRANCH=.TRUE.
                  ELSE IF(N_OLD_BRANCH(M).EQ.6) THEN
                    IF(N.EQ.1) THEN
                      NTMP=0
                      N_TEMP_ORD=23
                    ELSE IF(N.EQ.2) THEN
                      NTMP=0
                      N_TEMP_ORD=26
                    ENDIF !N
                    TCHECK=.FALSE.
                    PBRANCH=.TRUE.
                  ELSE IF(N_OLD_BRANCH(M).EQ.7) THEN
                    IF(N.EQ.1) THEN
                      NTMP=9
                      N_TEMP_ORD=28
                    ELSE IF(N.EQ.2) THEN
                      NTMP=8
                      N_TEMP_ORD=25
                    ENDIF !N
                    TCHECK=.TRUE.
                    PBRANCH=.TRUE.
                  ELSE IF(N_OLD_BRANCH(M).EQ.8) THEN
                    IF(N.EQ.1) THEN
                      NTMP=0
                      N_TEMP_ORD=21
                    ELSE IF(N.EQ.2) THEN
                      NTMP=0
                      N_TEMP_ORD=24
                    ENDIF !N
                    TCHECK=.FALSE.
                    PBRANCH=.TRUE.
                  ELSE IF(N_OLD_BRANCH(M).EQ.9) THEN
                    IF(N.EQ.1) THEN
                      NTMP=0
                      N_TEMP_ORD=24
                    ELSE IF(N.EQ.2) THEN
                      NTMP=0
                      N_TEMP_ORD=27
                    ENDIF !N
                    TCHECK=.FALSE.
                    PBRANCH=.TRUE.
                  ELSE IF(N.EQ.1) THEN !N is odd
                    N_PARORD =NE_OLD(M)
                    ord=N_PARORD-1
                    N_TEMP_ORD=ord
                    NTMP=0
                    TCHECK=.FALSE.
                    PBRANCH=.FALSE.
                  ELSE IF(N.EQ.2) THEN !N is even
                    N_PARORD =NE_OLD(M)
                    ord=N_PARORD-1-NDIFF_DAUGHTER(N_PARORD)
                    N_TEMP_ORD=ord
                    NTMP=0
                    TCHECK=.FALSE.
                    PBRANCH=.FALSE.
                  ENDIF !gen
                  ord=N_TEMP_ORD
                ELSE !gen greater than 5
                  N_PARORD=NE_OLD(M)
                  IF(N.EQ.1) THEN !N is odd
                    ord=N_PARORD-1
                  ELSE IF(N.EQ.2) THEN !N is even
                    ord=N_PARORD-1-NDIFF_DAUGHTER(N_PARORD)
                  ENDIF !N
                  N_TEMP_ORD=ord
                  TCHECK=.FALSE.
                ENDIF !gen
                IF(ord.GE.N2_ORDER) THEN !branch exists
                  BRANCH=.TRUE.
                ELSE
                  BRANCH=.FALSE.
                ENDIF !ord
                IF(BRANCH) THEN
                  IF(PBRANCH) THEN !Parent is a trachea branch
                    IF(TCHECK) THEN !Branch is a trachea branch
                      IF(gen.EQ.1) THEN !Calculate means for angles
                        ANG_Y_MN=0.0d0
                      ELSE
                        ANG_Y_MN=XP(1,1,4,NE_OLD(M))+B_ANGLE_Y(NTMP+90)
                      ENDIF !gen
                    ELSE IF(.NOT.TCHECK) THEN !parent is trach branch
                      IF(N.EQ.1) THEN !N is odd
                        N_DIRN=NINT(RANDOM_NUMBER(ISEED1))
                        IF(N_DIRN.EQ.0) THEN
                          ANG_Y_MN=XP(1,1,4,NE_OLD(M))-B_ANGLE_Y(ord)
                        ELSE
                          ANG_Y_MN=XP(1,1,4,NE_OLD(M))+B_ANGLE_Y(ord)
                        ENDIF !N_DIRN
                      ELSE IF(N.EQ.2) THEN !N is even
                        IF(N_DIRN.EQ.0) THEN
                          ANG_Y_MN=XP(1,1,4,NE_OLD(M))+B_ANGLE_Y(ord)
                        ELSE
                          ANG_Y_MN=XP(1,1,4,NE_OLD(M))-B_ANGLE_Y(ord)
                        ENDIF !N_DIRN
                      ENDIF !N
                      B_ANGLE_SD(ord)=B_ANGLE_CV*B_ANGLE_Y(ord)
                    ENDIF !TCHECK
                  ELSE
                    IF(gen.EQ.1) THEN !Calculate means for angles
                      ANG_Y_MN=0.0d0
                    ELSE IF(N.EQ.1) THEN !N is odd
                      N_DIRN=NINT(RANDOM_NUMBER(ISEED1))
                      IF(N_DIRN.EQ.0) THEN
                        ANG_Y_MN=XP(1,1,4,NE_OLD(M))-B_ANGLE_Y(ord)
                      ELSE
                        ANG_Y_MN=XP(1,1,4,NE_OLD(M))+B_ANGLE_Y(ord)
                      ENDIF !N_DIRN
                    ELSE IF(N.EQ.2) THEN !N is even
                      IF(N_DIRN.EQ.0) THEN
                        ANG_Y_MN=XP(1,1,4,NE_OLD(M))+B_ANGLE_Y(ord)
                      ELSE
                        ANG_Y_MN=XP(1,1,4,NE_OLD(M))-B_ANGLE_Y(ord)
                      ENDIF !N_DIRN
                    ENDIF !gen
                  ENDIF !PBRANCH
                  ANG_Y=NORM_RAND_NUM(ISEED_GEOM,ANG_Y_MN,
     '              B_ANGLE_SD(ord))
                  IF(NJT.EQ.3) THEN !3-dimensional mesh
                    IF(gen.EQ.1) THEN
                      ANG_XY_MN=0.0d0
                    ELSE IF(N.EQ.1) THEN !N is odd
                      N_DIRN=NINT(RANDOM_NUMBER(ISEED1))
                      IF(N_DIRN.EQ.0) THEN
                        ANG_XY_MN=XP(2,1,4,NE_OLD(M))-B_ANGLE_XY(ord)
                      ELSE
                        ANG_XY_MN=XP(2,1,4,NE_OLD(M))+B_ANGLE_XY(ord)
                      ENDIF !N_DIRN
                    ELSE IF(N.EQ.2) THEN !N is even
                      IF(N_DIRN.EQ.1) THEN
                        ANG_XY_MN=XP(2,1,4,NE_OLD(M))-B_ANGLE_XY(ord)
                      ELSE
                        ANG_XY_MN=XP(2,1,4,NE_OLD(M))+B_ANGLE_XY(ord)
                      ENDIF !N_DIRN
                    ENDIF !gen
                    ANG_XY=NORM_RAND_NUM(ISEED_GEOM,ANG_XY_MN,
     '                B_ANGLE_SD(ord))
                  ENDIF !NJT
                  IF(TCHECK) THEN !Takes length from branch #
                    B_LEN=B_LENGTH(NTMP+90)
                    B_LEN_SD=B_LENGTH_SD(NTMP+90)
                    B_DIAM=B_DIAMETER(NTMP+90)
                    B_DIAM_SD=B_DIAMETER_SD(NTMP+90)
                  ELSE !Takes length from order #
                    B_LEN=B_LENGTH(ord)
                    B_LEN_SD=B_LENGTH_SD(ord)
                    B_DIAM=B_DIAMETER(ord)
                    B_DIAM_SD=B_DIAMETER_SD(ord)
                  ENDIF !TCHECK
                  IF(JTYP14.EQ.3) THEN !regular tree
                    BRANCH_L=B_LEN
                    BRANCH_D=B_DIAM
                  ELSE IF(JTYP14.EQ.4) THEN !stochastic tree
                    BRANCH_L=NORM_RAND_NUM(ISEED_GEOM,B_LEN,B_LEN_SD)
                    BRANCH_D=NORM_RAND_NUM(ISEED_GEOM,B_DIAM,B_DIAM_SD)
                  ENDIF !JTYP14
                  NT_BNS=NT_BNS+1 !Update variables
                  np=np+1
                  ne=NE+1
                  CE(1,ne)=gen
                  CE(2,ne)=BRANCH_L
                  CE(3,ne)=BRANCH_D
                  CE(4,ne)=PI*CE(3,ne)**2*0.25d0*CE(2,ne) !Branch vol
                  NW(ne,1)=ord
                  B_VOLUME(ord)=B_VOLUME(ord)+CE(4,ne)
                  NPNE(1,nb,ne)=NE_OLD(M) !node # at start of branch
                  NPNE(2,nb,ne)=np !is node # at end of branch
C GMH 8/1/97 Update cmgui link
                  CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                  XP(1,1,1,np)=XP(1,1,1,NE_OLD(M))+BRANCH_L*DSIN(ANG_Y)
     '              *DCOS(ANG_XY)
                  XP(1,1,2,np)=XP(1,1,2,NE_OLD(M))-BRANCH_L*DCOS(ANG_Y)
                  XP(1,1,3,np)=XP(1,1,3,NE_OLD(M))
     '              -BRANCH_L*DSIN(ANG_Y)*DSIN(ANG_XY)
                  XP(1,1,4,np)=ANG_Y
                  IF(NJT.EQ.3) XP(2,1,4,np)=ANG_XY
                  NE_TEMP(NT_BNS)=np !records node #s at this gen
                  NE_TEMP(NT_BNS)=N_TEMP_ORD
                  N_TEMP_BRANCH(NT_BNS)=NTMP
                ENDIF !BRANCH
              ENDDO !N
              NUM_SPLIT=2
            ENDDO !M
            NO_BRANCHES(gen)=NT_BNS
            NUM_NODES=NT_BNS
            DO N=1,NUM_NODES
              NE_OLD(N)=NE_TEMP(N)
              NE_OLD(N)=NE_TEMP(N) !updates parent order#s
              IF(gen.LE.5) THEN
                N_OLD_BRANCH(N)=N_TEMP_BRANCH(N) !upd parent branch #s
              ENDIF !gen
            ENDDO !N
          ENDDO !gen
          NPT(nr)=np !highest node# in region nr
          NET(nr)=ne !highest element# in region nr
          necond=ne !number of elements in model
          npcond=np !number of nodes in model
        ELSE IF (CONDUCTING_TYPE.EQ.2)THEN
 !calculate the generations for the central airways
          IF(INDEX_IP2.EQ.1)THEN
            nr=INDEX_REG(nr_high) !region number for central airways
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              np1=NPNE(1,nba,ne) !start node of element ne
              np2=NPNE(2,nba,ne) !end node of element ne
              IF(noelem.EQ.1)THEN
 !assume the first branch is the trachea
                CE(1,ne)=1.0d0
              ELSE
                ne1=NXI(-1,1,ne) !parent global element #
                CE(1,ne)=CE(1,ne1)+1.0d0 !increment generation by one
              ENDIF !noelem
              sum=0.0d0
              DO nj=1,NJT !calculate the length of the branch
                sum=sum+(XP(1,1,nj,np2)-XP(1,1,nj,np1))**2.d0
              ENDDO !nj
              CE(2,ne)=DSQRT(sum) !length of branch element
            ENDDO !noelem (ne)
            DO nreg=1,nr_high-1
              nr=INDEX_REG(nreg) !for each lobe
              nb=nbl
              IF(RANDC_METHOD.EQ.1)THEN
C             check that the random numbers have been read in
                CALL ASSERT(NDT.GE.0,'>>Define data points',ERROR,*9999)
                DO nd1=1,NDT !can only read in rps
                  NDP(nd1)=NDT !when generating a single
                  LD(nd1)=BRN_STORE(nreg) !lobe mesh
                  DO nj=1,NJT
                    WD(nj,NDT)=1.0d0
                  ENDDO !nj
                ENDDO !nd1
              ELSE IF(RANDC_METHOD.EQ.2)THEN !generate random numbers
                DO nb2=1,NBFT
                  IBT(1,1,nb2)=1
                  IBT(2,1,nb2)=1
                ENDDO !nb2
                NDT=0
                DO noelem=1,NEELEM(0,nr) !for each element in the lobe
                  ne=NEELEM(noelem,nr) !global host element #
                  CALL XPXE(NBJ(1,ne),NKE(1,1,1,ne),NPF(1,1),
     '              NPNE(1,1,ne),nr,NVJE(1,1,1,ne),SE(1,1,ne),XA,XE,XP,
     '              ERROR,*9999)
                  NDT2=0
                  DO nd1=1,NDT_ELEM !generate NDT_ELEM random points
                    nu=1
                    DO nj=1,NJT !random point in element coords
                      XI(nj)=RANDOM_NUMBER(ISEED_GEOM)
                    ENDDO !nj
                    DO nj=1,NJT !random point in world coordinates
                      X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '                  nb,nu,XI,XE(1,nj)) !in world coordinates
                    ENDDO !nj
                    FAIL=.FALSE.
                    nd2=0
                    DO WHILE((nd2.LT.NDT2).AND.(.NOT.FAIL))
                      nd2=nd2+1
                      nd3=NDT-nd2+1
                      DIST=0.0d0 
                      DO nj=1,NJT
                        DIST=DIST+(ZD(nj,nd3)-X(nj))**2.d0
                      ENDDO !nj
                      DIST=DSQRT(DIST)
                      IF(DIST.LT.LIMIT_ND) FAIL=.TRUE.
                    ENDDO
                    IF(.NOT.FAIL)THEN !add to random points
                      NDT=NDT+1    !counts random points for whole lobe
                      NDT2=NDT2+1  !counts random points for element
                      NDP(NDT)=NDT !assigns points to lobar bronchus
                      LD(NDT)=BRN_STORE(nreg)
                      DO nj=1,NJT
                        ZD(nj,NDT)=X(nj) !stores random point coords
                        WD(nj,NDT)=1.0d0
                      ENDDO !nj
                    ENDIF !DIST
                  ENDDO !nd1
                ENDDO !noelem (ne)
              ENDIF !RANDC_METHOD
              WRITE(OP_STRING,'('' number of random points: '',I7)') NDT
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C             generate branches for the current lobe
              ne=NET(0) !highest element number
              np=NPT(0) !highest node number
              noelem=NEELEM(0,nr) !number of elements in region nr
              nonode=NPNODE(0,nr) !number of nodes in region nr
              ne1=BRN_STORE(nreg) !global element # of lobar bronchus
              nb=nba
              Lgen=NINT(CE(1,ne1)) !generation of lobar bronchus
              gen=Lgen
              N_GEN=1
              N_ELM=1
              NE_OLD(1)=BRN_STORE(nreg)
              WRITE(OP_STRING,'(''   gen  #brn  #TBs'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              DO WHILE(N_GEN.NE.0)
                N_GEN=0
                gen=gen+1
                NT_BNS=N_ELM
                N_ELM=0
                N_TRM=0
                DO nd=1,NDT
                  LDR(nd)=0
                ENDDO !nd
                IF(nr.EQ.6)THEN
                  DO M=1,NT_BNS
                    N_SPACE(NE_OLD(M))=NE_OLD(M)
                  ENDDO
                ENDIF
                DO M=1,NT_BNS
                  ne1=NE_OLD(M) !parent global element #
                  ne2=NXI(-1,1,ne1) !grandparent global element #
                  np1=NPNE(2,nba,ne1) !parent global end node #
                  np2=NPNE(1,nba,ne1) !parent global start node #
                  np3=NPNE(1,nba,ne2) !grandparent global start node #
                  IF(NPNE(1,nba,ne1+1).EQ.np2)THEN
                    np4=NPNE(2,nba,ne1+1) !aunty end node number
                  ELSE
                    np4=NPNE(2,nba,ne1-1)
                  ENDIF !NPNE
                  CALL SPLIT_SPACE(LD,ne1,ne,np1,np2,np3,XP,ZD,SS,ERROR,
     '              *9999)
                  DO N=1,2
                    noelem=noelem+1
                    nonode=nonode+1
                    ne=ne+1
                    np=np+1
                    CALL CALC_COFM(LD,ne,COFM,ZD,ERROR,*9999)
                    FAIL=.TRUE.
                    noelem2=0
                    DO nj=1,NJT
                      X(nj)=COFM(nj)
                    ENDDO !nj
                    DO WHILE((FAIL).AND.noelem2.LE.NEELEM(0,nr))
                      noelem2=noelem2+1
                      ne_host=NEELEM(noelem2,nr)
                      nbh=NBJ(nj,ne_host)
                      IF(nbh.EQ.nbl)THEN
                        CALL XPXE(NBJ(1,ne_host),NKE(1,1,1,ne_host),
     '                    NPF(1,1),NPNE(1,1,ne_host),nr,
     '                    NVJE(1,1,1,ne_host),SE(1,1,ne_host),XA,XE,XP,
     '                    ERROR,*9999)
                        CALL FINDXI(IBT,IDO,INP,NBJ(1,ne_host),X,XI,XE,
     '                    FAIL,ERROR,*9999)
                      ENDIF
                    ENDDO !WHILE
                    CALL BRANCH_COFM(np1,np2,ANG,COFM,FRAC(gen),A_LIM,
     '                L_LIM,XP,XP1,BRANCH,FAIL,ERROR,*9999)
                    CALL ASSERT(np.LE.NPM,'>>NPM too small',ERROR,*9999)
                    CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                    DO nj=1,NJT
                      XP(1,1,nj,np)=XP1(nj)
                      XP(2,1,nj,np)=XP1(nj)
                    ENDDO !nj
                    IF(BRANCH.AND.SS(N))THEN
                      N_GEN=N_GEN+1 !short branches are TB
                      N_ELM=N_ELM+1
                      NE_TEMP(N_ELM)=ne
                    ENDIF !BRANCH
                    IF(gen.EQ.Lgen)THEN
                      NXI(N,1,ne1)=noelem
                      NXI(1,0,ne1)=2
                    ENDIF
                    NEELEM(noelem,nr)=ne
                    NPNODE(nonode,nr)=np
                    NPNE(1,nba,ne)=np1
                    NPNE(2,nba,ne)=np
                    NXI(-1,1,ne)=ne1
                    NXI(1,1,ne)=noelem
                    CE(1,ne)=gen
                    LEN=0.0d0
                    DO nj=1,NJT
                      LEN=LEN+(XP(1,1,nj,np)-XP(1,1,nj,np1))**2.0d0
                    ENDDO !nj
                    CE(2,ne)=DSQRT(LEN) !length of branch
                    CE(16,ne)=ANG !branch angle to parent (radians)
                    N_SPACE(N)=N_SPACE(ne1)
                    IF((.NOT.BRANCH).OR.(.NOT.SS(N)))THEN
                      CE(5,ne)=1
                      N_TRM=N_TRM+1
                      NE_TERM(N_TRM)=ne !element # of TB
                      DO nd=1,NDT
                        nsp=LD(nd)
                        IF(nsp.EQ.ne)THEN !point to reassign
                          LDR(nd)=ne !store space #
                          LD(nd)=0
                        ENDIF !nsp
                      ENDDO !nd
                    ENDIF !.NOT.BRANCH
                  ENDDO !N
                  CALL CALC_PLANE_ANGLE(nonode,np,np1,np2,np4,PANG1,
     '              PANG2,XP,ERROR,*9999)
                  CE(13,ne)=PANG1 !plane angle
                  CE(13,ne-1)=PANG2 !plane angle
                ENDDO !M
                WRITE(OP_STRING,'(3(I6))')gen,N_ELM,N_TRM
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                DO N=1,N_ELM
                  NE_OLD(N)=NE_TEMP(N)
                ENDDO !N
                IF(N_TRM.NE.0)THEN
                  IF(DOP)THEN
                    WRITE(OP_STRING,'('' #TBs   gen   g_lim'',3(I6))')
     '                N_TRM,gen,G_LIM
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF                    
                  IF(gen.LE.G_LIM)THEN
                    CALL REASSTB(LD,LDR,nb,NE_OLD,NE_TERM,NPC_MIN,
     '                NPNE,N_ELM,N_TRM,MIN_DIST,XP,ERROR,*9999)
                  ENDIF !gen
                ENDIF !N_TRM
              ENDDO !WHILE N_GEN
              NET(nr)=ne
              NPT(nr)=np
              IF(ne.GT.NET(0)) NET(0)=ne
              IF(np.GT.NPT(0)) NPT(0)=np
              NEELEM(0,nr)=noelem !# of elements in region nr
              NPNODE(0,nr)=nonode !# of nodes in region nr
              NEELEM(0,0)=NEELEM(0,0)+noelem !# of elements in mesh
              NPNODE(0,0)=NPNODE(0,0)+nonode !# of nodes in mesh
            ENDDO !nreg (nr)
            necond=NEELEM(0,0)
            npcond=NPNODE(0,0)
            DO nreg=1,nr_high
              nr=INDEX_REG(nreg)
              nb=nba
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                IF(CE(1,ne).NE.0.0d0)THEN !branching element
                  NRE(ne)=nr
                  NJE(ne)=NJT
                  DO nj=1,NJE(ne)
                    NBJ(nj,ne)=nb
                  ENDDO !nj
                  NCO(ne)=1
                  DO ns=1,NST(nb)+NAT(nb)
                    SE(ns,nb,ne)=1.0d0
                  ENDDO !ns
                  DO nn=1,NNT(nb)
                    DO nk=1,NKT(nn,nb)
                      NKE(nk,nn,nb,ne)=1
                    ENDDO !nk
                    DO nj=1,NJE(ne)
                      NVJE(nn,nb,nj,ne)=1
                    ENDDO !nj
                  ENDDO !nn
                ENDIF !CE
              ENDDO !noelem(ne)
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                IF(CE(1,ne).NE.0.0d0)THEN !branching element
                  DO nj=1,NJT
                    NKJ(nj,np)=1
                    DO nc=1,NCM
                      NVJP(nj,np)=1
                    ENDDO !nc
                  ENDDO !nj
                ENDIF !CE
              ENDDO !noelem(ne)
            ENDDO !nreg (nr)
            CALL NENXI(IBT,INP,NBJ,NEELEM,NPNE,NXI,ERROR,*9999)
            DO nreg=1,nr_high-1 !for the branches across regions
              nr=INDEX_REG(nreg) !lobe region
              nb=nba
              ne1=BRN_STORE(nreg) !lobar bronchus element #
              np1=NPNE(2,nb,ne1) !end node # of lobar bronchus
              NXI(1,0,ne1)=0
              DO noelem2=1,NEELEM(0,nr) !check each element
                ne2=NEELEM(noelem2,nr)
                IF(CE(1,ne2).NE.0.0d0)THEN !a branching element
                  np2=NPNE(1,nb,ne2) !start node # of branch
                  IF(np1.EQ.np2)THEN
                    NXI(1,0,ne1)=NXI(1,0,ne1)+1 !incr # of daughters
                    NTMP=NXI(1,0,ne1)
                    IF(DOP)THEN
                      WRITE(OP_STRING,'('' ne1,np1,ne2,np2,#d: '','
     '                  //'5(I4))') ne1,np1,ne2,np2,NTMP
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF                    
                    NXI(1,NTMP,ne1)=ne2 !the daughter element #
                  ENDIF
                ENDIF
              ENDDO !noelem2 (ne2)
            ENDDO !nreg (nr)
          ELSE IF(INDEX_IP2.EQ.2)THEN
            nb=nba
            ne=TRAC !element # of trachea
            CE(1,ne)=1.d0 !first generation
            np1=NPNE(1,nb,ne)
            np2=NPNE(2,nb,ne)
            sum=0.d0
            DO nj=1,NJT
              sum=sum+(XP(1,1,nj,np2)-XP(1,1,nj,np1))**2.d0
            ENDDO !nj
            CE(2,ne)=DSQRT(sum)
            gen=1
            NT_BNS=1
            NE_OLD(1)=ne
            DO WHILE(NT_BNS.NE.0)
              NUM_NODES=NT_BNS
              NT_BNS=0
              gen=gen+1
              DO M=1,NUM_NODES
                ne=NE_OLD(M)
                IF(NXI(1,0,ne).NE.0)THEN
                  DO N=1,2
                    NT_BNS=NT_BNS+1
                    ne1=NXI(1,N,ne) !Nth daughter branch
                    CE(1,ne1)=DBLE(gen)
                    np1=NPNE(1,nb,ne1)
                    np2=NPNE(2,nb,ne1)
                    sum=0.d0
                    DO nj=1,NJT
                      sum=sum+(XP(1,1,nj,np2)-XP(1,1,nj,np1))**2.d0
                    ENDDO !nj
                    CE(2,ne1)=DSQRT(sum)
                    NE_TEMP(NT_BNS)=ne1 !records elements in this gen
                  ENDDO !N
                ENDIF
              ENDDO !M
              DO N=1,NT_BNS
                NE_OLD(N)=NE_TEMP(N) !updates list of prev gen elem#s
              ENDDO !N
            ENDDO !WHILE
            NTT_GEN=gen
            DO nreg=1,6
              nr=INDEX_REG(nreg)
              nb=nba
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                NRE(ne)=nr
                NJE(ne)=NJT
                DO nj=1,NJE(ne)
                  NBJ(nj,ne)=nb
                ENDDO !nj
                NCO(ne)=1
                DO ns=1,NST(nb)+NAT(nb)
                  SE(ns,nb,ne)=1.0d0
                ENDDO !ns
                DO nn=1,NNT(nb)
                  DO nk=1,NKT(nn,nb)
                    NKE(nk,nn,nb,ne)=1
                  ENDDO !nk
                  DO nj=1,NJE(ne)
                    NVJE(nn,nb,nj,ne)=1
                  ENDDO !nj
                ENDDO !nn
              ENDDO !noelem(ne)
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                DO nj=1,NJT
                  NKJ(nj,np)=1
                  DO nc=1,NCM
                    NVJP(nj,np)=1
                  ENDDO !nc
                ENDDO !nj
              ENDDO !noelem(ne)
            ENDDO !nreg (nr)
          ENDIF !INDEX_IP2
          DO nreg=1,nr_high !calculate the branch orders
            nr=INDEX_REG(nreg)
            DO noelem=NEELEM(0,nr),1,-1
              ne=NEELEM(noelem,nr)
              IF(CE(1,ne).NE.0.0d0)THEN
                ord=0
                IF(NXI(1,0,ne).NE.0)THEN !branch has daughters
                  DO noelem2=1,2 !for both daughters
                    ne2=NXI(1,noelem2,ne) !global element # of daughter
                    temp=NINT(CE(10,ne2)) !order of daughter element
                    IF(temp.GT.ord) ord=temp
                  ENDDO !noelem2 (ne2)
                ENDIF !NXI
                ord=ord+1 !Horsfield ordering
                CE(10,ne)=DBLE(ord) !store the branch order
              ENDIF !CE
            ENDDO !noelem (ne)
          ENDDO !nreg (nr)
C       Calculate and apply the diameters etc.
          DO nreg=1,nr_high
            nr=INDEX_REG(nreg)
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              ord=NINT(CE(10,ne)) !Horsfield order of branch
              BRANCH_L=CE(2,ne) !branch length
              BRANCH_D=BRANCH_L**(1.d0/3.d0) !calulation for diameter
              CE(3,ne)=(PI*BRANCH_D**2.d0)/4.d0 !xsec area
              CE(4,ne)=CE(2,ne)*CE(3,ne) !branch volume
              CE(5,ne)=CE(4,ne) !initialise for volumes below
              CE(9,ne)=1.d0 !ratio of duct to alveolar+duct xsec area
              CE(12,ne)=CE(2,ne) !initial branch length
              CE(13,ne)=BRANCH_D !initial branch diameter
              CE(14,ne)=CE(4,ne) !initial branch volume
            ENDDO !noelem (ne)
          ENDDO !nreg (nr)
        ELSE IF(CONDUCTING_TYPE.EQ.3)THEN
          nb=nba
          ne=0
          np=1
          necond=0
          npcond=1
          NPNODE(np,nr)=np
          XP(1,1,1,np)=0.0d0
          XP(1,1,2,np)=0.0d0
          XP(1,1,3,np)=0.0d0
          XP(1,1,4,np)=0.0d0
          XP(2,1,4,np)=0.0d0
          NXI(-1,1,1)=0
          DO gen=1,NTT_GEN
            DIVLENGTH=B_LENGTH(gen)/B_DIV(gen)
            DO N=1,B_DIV(gen)
              np=np+1
              ne=ne+1
              CE(1,ne)=gen !branch generation
              CE(2,ne)=DIVLENGTH !element length
              CE(13,ne)=B_DIAMETER(gen) !element diameter
              CE(3,ne)=(PI*CE(13,ne)**2.0d0)/4.0d0 !xsec area
              CE(4,ne)=CE(3,ne)*CE(2,ne) !volume
              CE(5,ne)=CE(4,ne) !initial volume below
              CE(9,ne)=1.0d0 !no alveolar area
              CE(10,ne)=CE(3,ne) !initial xsec area
              CE(12,ne)=CE(2,ne) !initial branch length
              CE(14,ne)=CE(4,ne) !initial branch volume
              NPNE(1,1,ne)=np-1 !start node number
              NPNE(2,1,ne)=np !end node number
              NEELEM(ne,nr)=ne !global element number
              NPNODE(np,nr)=np !global node number
              XP(1,1,1,np)=0.0d0 !x coord
              XP(1,1,2,np)=0.0d0 !y coord
              XP(1,1,3,np)=XP(1,1,3,np-1)-CE(2,ne) !z coord
              XP(1,1,4,np)=0.0d0
              XP(2,1,4,np)=0.0d0
              NXI(-1,0,ne)=1 !one adjacent element in -Xi(1) direction
              NXI(-1,1,ne)=ne-1 !global element number of -Xi(1) element
              NXI(1,0,ne-1)=0
              IF(ne.NE.1) NXI(1,0,ne-1)=1
              NXI(1,1,ne)=ne+1
              necond=necond+1
              npcond=npcond+1
              NRE(ne)=nr
            ENDDO !N
            CE(15,ne)=1.0d0 !identifies the end of a true branch
          ENDDO !gen
          CE(15,ne)=0.0d0 !end of terminal bronchiole
          NPT(nr)=np
          NET(nr)=ne
          NPT(0)=NPT(nr)
          NET(0)=NET(nr)
          NEELEM(0,nr)=NET(nr)
          NPNODE(0,nr)=NPT(nr)
        ENDIF !CONDUCTING_TYPE
      ENDIF !MODEL_TYPE
      IF((MODEL_TYPE.EQ.1).OR.(MODEL_TYPE.EQ.3))THEN
        IF(RESPIRATORY_TYPE.EQ.1)THEN
          IF(MODEL_TYPE.EQ.3)THEN !Isolated acinar model
C GMH 8/1/97 Update cmgui link
            CALL NODE_CHANGE(1,.FALSE.,ERROR,*9999)
            XP(1,1,1,1)=0.0d0 !node 1 x is at centre
            XP(1,1,2,1)=0.0d0 !node 1 y is at bottom
            XP(1,1,3,1)=0.0d0 !node 1 z is at bottom
            XP(1,1,4,1)=0.0d0 !node 1 branch angle is zero
            XP(2,1,4,1)=0.0d0 !node 1 branch angle (3D) is zero
            nb=nba
            np=1
            ne=0
            NET(nr)=ne
            NPT(nr)=np
            PE_TEMP=0
            necond=0
            npcond=0
          ELSE IF(MODEL_TYPE.EQ.1)THEN !Attached to cond airway tube
            np=NPT(nr)
            ne=NET(nr)
            PE_TEMP=ne
          ENDIF !MODEL_TYPE
          NBT=1
          NFBASE(1,1)=1
          NBASEF(1,0)=1
          NBASEF(1,1)=1
          NBFT=1
          NE_OLD(1)=np !NE_OLD(n) stores previous generation node #s
          NT_BNS=1 !Total number of branches in generation
          DO nogen=1,NT_GEN
            IF(CONDUCTING_TYPE.EQ.3)THEN
              gen=nogen+NTT_GEN
              TEMPRAND=NORM_RAND_NUM(ISEED_BRANCH,MNB(NTT_GEN),
     '          SDNB(NTT_GEN))
            ELSE IF(MODEL_TYPE.EQ.3)THEN
              gen=nogen
            ENDIF
            NUM_NODES=NT_BNS
            NT_BNS=0
            DO M=1,NUM_NODES
              TEMPRAND=NORM_RAND_NUM(ISEED_BRANCH,MNB(gen),SDNB(gen))
              INUM_BRANCH=NINT(TEMPRAND)
              IF(INUM_BRANCH.LE.1.AND.gen.GT.NTT_GEN+1) INUM_BRANCH=0
              NO_BRANCH_PATTERN(gen,INUM_BRANCH)=
     '          NO_BRANCH_PATTERN(gen,INUM_BRANCH)+1
              DO N=1,INUM_BRANCH
                NT_BNS=NT_BNS+1
                np=np+1
                ne=ne+1
                NET(nr)=NET(nr)+1
                NPT(nr)=NPT(nr)+1
                CE(1,ne)=gen
                CE(2,ne)=N !!Temporarily
                NW(ne,1)=gen
                NPNE(1,nb,ne)=NE_OLD(M) !node# at start of branch
                NPNE(2,nb,ne)=np !node# at end of current branch
                NE_TEMP(NT_BNS)=np !records node#s at this gen
                PN_TEMP(ne)=PE_TEMP
                NXI(-1,1,ne)=PE_TEMP
              ENDDO !N
              PE_TEMP=PE_TEMP+1
            ENDDO !M
            NO_BRANCHES(gen)=NT_BNS
            DO N=1,NT_BNS
              NE_OLD(N)=NE_TEMP(N) !updates list of prev gen node#s
            ENDDO !N
          ENDDO !nogen
          DO ne1=necond+1,NET(nr)
            N=INT(CE(2,ne1))
            gen=INT(CE(1,ne1))
            NPOLD=NPNE(1,nb,ne1)
            np1=NPNE(2,nb,ne1)
            ne2=PN_TEMP(ne1) !parent element number
            IF(N.EQ.1)THEN !Calculates first branch
              N_DIRN=NINT(RANDOM_NUMBER(ISEED1))
              IF(N_DIRN.EQ.0) THEN
                ANG_Y_MN=XP(1,1,4,NPOLD)+W_ANGLE_Y(gen)
              ELSE
                ANG_Y_MN=XP(1,1,4,NPOLD)-W_ANGLE_Y(gen)
              ENDIF !N_DIRN
            ELSE IF(N.EQ.2) THEN !Calculates second branch
              IF(N_DIRN.EQ.0) THEN
                ANG_Y_MN=XP(1,1,4,NPOLD)-W_ANGLE_Y(gen)
              ELSE
                ANG_Y_MN=XP(1,1,4,NPOLD)+W_ANGLE_Y(gen)
              ENDIF !N_DIRN
            ELSE !Trichotomy or more
              ANG_Y_MN=XP(1,1,4,NPOLD)
            ENDIF !N
            ANG_Y=NORM_RAND_NUM(ISEED_GEOM,ANG_Y_MN,W_ANGLE_SD(gen))
            IF(NJT.EQ.3)THEN !3D Tree
              IF(N.EQ.1)THEN !Calculates first branch
                ANG_XY_MN=XP(2,1,4,NPOLD)-B_ANGLE_XY(gen)
              ELSE IF(N.EQ.2)THEN !Calculates second branch
                ANG_XY_MN=XP(2,1,4,NPOLD)+B_ANGLE_XY(gen)
              ELSE IF(N.EQ.3)THEN !Trichotomy
                ANG_XY_MN=0.0d0
              ENDIF !N
              ANG_XY=NORM_RAND_NUM(ISEED_GEOM,ANG_XY_MN,B_ANGLE_SD(gen))
            ELSE
              ANG_XY=0.0d0
            ENDIF !NJT
            IF(JTYP14.EQ.3) THEN !regular tree
              BRANCH_L=B_LENGTH(gen)
              BRANCH_D=B_DIAMETER(gen)
              BRANCH=.TRUE.
            ELSE IF(JTYP14.EQ.4) THEN !stochastic tree
              BRANCH_L=0.0d0
              DO WHILE ((BRANCH_L.LE.0.0d0).OR.(BRANCH_D.LE.0.0d0))
                BRANCH_L=NORM_RAND_NUM(ISEED_GEOM,B_LENGTH(gen),
     '            B_LENGTH_SD(gen))
                BRANCH_D=NORM_RAND_NUM(ISEED_GEOM,B_DIAMETER(gen),
     '            B_DIAMETER_SD(gen))
              ENDDO !BRANCH_L
            ENDIF !JTYP14
            CE(1,ne1)=gen
            CE(2,ne1)=BRANCH_L
            CE(13,ne1)=BRANCH_D !branch diameter
            CE(9,ne1)=S_RATIO(gen) !branch diameter ratio
            CE(3,ne1)=(PI*BRANCH_D**2.0d0)/(4.0d0*CE(9,ne1))
            CE(10,ne1)=CE(3,ne1)
            CE(4,ne1)=CE(2,ne1)*CE(3,ne1) !vol by cylinder         
            CE(12,ne1)=BRANCH_L
            CE(14,ne1)=CE(4,ne1) !initial volumes
            CE(5,ne1)=CE(14,ne1) !initialise volume below 
            B_VOLUME(gen)=B_VOLUME(gen)+CE(14,ne1)
C GMH 8/1/97 Update cmgui link
            CALL NODE_CHANGE(NP1,.FALSE.,ERROR,*9999)
            XP(1,1,1,NP1)=XP(1,1,1,NPOLD)+BRANCH_L*DSIN(ANG_Y)*
     '        DCOS(ANG_XY)
            XP(1,1,2,NP1)=XP(1,1,2,NPOLD)-BRANCH_L*DCOS(ANG_Y)
            XP(1,1,3,NP1)=XP(1,1,3,NPOLD)-BRANCH_L*DSIN(ANG_Y)*
     '        DSIN(ANG_XY)
            XP(1,1,4,NP1)=ANG_Y
            IF(NJT.EQ.3) THEN
              XP(2,1,4,NP1)=ANG_XY
            ELSE
              XP(2,1,4,NP1)=0.0d0
            ENDIF !NJT
          ENDDO !ne1
        ELSE IF(RESPIRATORY_TYPE.EQ.2)THEN
          DO nreg=1,nr_high
            nr=INDEX_REG(nreg)
            cnt=0
            necond=0
            npcond=0
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              IF(NXI(1,0,ne).EQ.0)THEN
                cnt=cnt+1
                BBM(1,cnt,nr)=ne
                BBM(2,cnt,nr)=BBM_VOL
                BBM(3,cnt,nr)=BBM_VOL
                BBM(4,cnt,nr)=BBM_VOL
                BBM(5,cnt,nr)=0.0d0
                BBM(6,cnt,nr)=0.0d0
                BBM(7,cnt,nr)=0.0d0
                BBM(8,cnt,nr)=0.0d0
                BBM(9,cnt,nr)=0.0d0
                BBM(10,cnt,nr)=0.0d0
                CE(5,ne)=CE(5,ne)+BBM_VOL
                CE(15,ne)=0.0d0 
              ENDIF !nxi
            ENDDO !noelem (ne)
            NTB=cnt
          ENDDO !nreg (nr)
        ENDIF !RESPIRATORY_TYPE
      ENDIF !MODEL_TYPE
      IF((MODEL_TYPE.EQ.1).OR.(MODEL_TYPE.EQ.2))THEN
        IF(CONDUCTING_TYPE.EQ.1)THEN
          CALL ASSERT(NPT(nr).LE.NPM,'>>NPM too small',ERROR,*9999)
          CALL ASSERT(NET(nr).LE.NEM,'>>NEM too small',ERROR,*9999)
          NPT(0)=NPT(nr) !highest node# in all regions
          NET(0)=NET(nr) !highest element# in all regions
          NPNODE(0,nr)=NPT(nr) !number nodes in region nr
          NPNODE(0, 0)=NPT(0) !number nodes in all regions
          NEELEM(0,nr)=NET(nr) !number elements in region nr
          NEELEM(0, 0)=NET(0) !number elements in all regions
          DO ne=1,NEELEM(0,nr)
            NEELEM(ne,nr)=ne
            NRE(ne)=nr
            NP1=NPNE(1,nb,ne) !is node # at start of current branch
            NP2=NPNE(2,nb,ne) !is node # at end of current branch
C GMH 8/1/97 Update cmgui link
            CALL NODE_CHANGE(NP2,.FALSE.,ERROR,*9999)
            XP(2,1,1,NP2)=XP(1,1,1,NP2)
            XP(2,1,2,NP2)=XP(1,1,2,NP2)
            XP(2,1,3,NP2)=XP(1,1,3,NP2)
          ENDDO !ne
          DO ne=NEELEM(0,nr),1,-1
            ne1=NEELEM(ne,nr) !element number
            ne2=NXI(-1,1,ne1) !parent element number
          ENDDO !ne
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            NJE(ne)=NJT
            DO nj=1,NJE(ne)
              NBJ(nj,ne)=1
            ENDDO !nj
            IF(DABS(CE(5,1)).GT.1.0D-6) THEN
              CE(11,ne)=CE(4,ne)/CE(5,1)
            ELSE
              CE(11,ne)=CE(4,ne)
              WRITE(OP_STRING,'('' WARNING... TOT_INIT_VOL=0'')') 
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF !CE(5,1)
          ENDDO !noelem
          DO nonode=1,NPNODE(0,nr)
            np=nonode
            NPNODE(nonode,nr)=NP
            DO nj=1,NJT
              NKJ(nj,np)=1
              DO nc=1,NCM
                NVJP(nj,np)=1
              ENDDO !nc
            ENDDO !nj
          ENDDO !nonode (np)
          DO nb=1,NBFT
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              NCO(ne)=1
              DO ns=1,NST(nb)+NAT(nb)
                SE(ns,nb,ne)=1.0d0
              ENDDO !ns
              DO nn=1,NNT(nb)
                DO nk=1,NKT(nn,nb)
                  NKE(nk,nn,nb,ne)=NK
                ENDDO !nk
                DO nj=1,NJE(ne)
                  NVJE(nn,nb,nj,ne)=1 !version one of nn,nj in elem ne
                ENDDO !nj
              ENDDO !nn
            ENDDO !noelem (ne)
          ENDDO !nb
        ENDIF !CONDUCTING_TYPE
      ENDIF !MODEL_TYPE
      IF((MODEL_TYPE.EQ.1).OR.(MODEL_TYPE.EQ.3))THEN
        IF(CONDUCTING_TYPE.NE.2)THEN
          CALL ASSERT(NPT(nr).LE.NPM,'>>NPM too small',ERROR,*9999)
          CALL ASSERT(NET(nr).LE.NEM,'>>NEM too small',ERROR,*9999)
          NPT(0)=NPT(nr) !highest node# in all regions
          NET(0)=NET(nr) !highest element# in all regions
          NPNODE(0,nr)=NPT(nr) !number nodes in region nr
          NPNODE(0, 0)=NPT(0) !number nodes in all regions
          NEELEM(0,nr)=NET(nr) !number elements in region nr
          NEELEM(0, 0)=NET(0) !number elements in all regions
          nb=nba
          DO ne=necond+1,NEELEM(0,nr)
            NEELEM(ne,nr)=ne
            NRE(ne)=nr
            NP1=NPNE(1,nb,ne) !is node # at start of current branch
            NP2=NPNE(2,nb,ne) !is node # at end of current branch
C GMH 8/1/97 Update cmgui link
            CALL NODE_CHANGE(NP2,.FALSE.,ERROR,*9999)
            XP(2,1,1,NP2)=XP(1,1,1,NP2)
            XP(2,1,2,NP2)=XP(1,1,2,NP2)
            XP(2,1,3,NP2)=XP(1,1,3,NP2)
          ENDDO !ne
          DO ne=NEELEM(0,nr),necond+1,-1
            ne1=NEELEM(ne,nr) !element number
            ne2=NXI(-1,1,ne1) !parent element number
            IF(MODEL_TYPE.EQ.3)THEN
              IF(ne2.NE.0) CE(5,ne2)=CE(5,ne2)+CE(5,ne1)
            ELSE IF(CONDUCTING_TYPE.EQ.3)THEN
              IF(RESPIRATORY_TYPE.EQ.1)THEN
                IF(ne2.GT.necond-1)THEN !beyond the cond airway tube
                  CE(5,ne2)=CE(5,ne2)+CE(5,ne1)
                ELSE !part of the conducting airway tube
                  IF(ne2.NE.0)THEN !not top of trachea
                    IF(CE(15,ne2).EQ.1.0d0)THEN !bottom of a tube branch
                      CE(5,ne2)=CE(5,ne2)+2.0d0*CE(5,ne1)
                    ELSE !within a tube branch
                      CE(5,ne2)=CE(5,ne2)+CE(5,ne1)
                    ENDIF
                  ENDIF
                ENDIF !ne2
              ELSE IF(RESPIRATORY_TYPE.EQ.2)THEN
                IF(ne2.NE.0)THEN !not top of trachea
                  IF(CE(15,ne2).EQ.1.0d0)THEN !bottom of a tube branch
                    CE(5,ne2)=CE(5,ne2)+2.0d0*CE(5,ne1)
                  ELSE !within a tube branch
                    CE(5,ne2)=CE(5,ne2)+CE(5,ne1)
                  ENDIF
                ENDIF !ne2
              ENDIF
            ENDIF
          ENDDO !ne
          IF(CONDUCTING_TYPE.EQ.3)THEN
            DO ne=necond,1,-1
              ne1=NEELEM(ne,nr)
              ne2=NXI(-1,1,ne1)
              IF(ne2.NE.0)THEN !not top of trachea
                IF(CE(15,ne2).EQ.1.0d0)THEN !bottom of a tube branch
                  CE(5,ne2)=CE(5,ne2)+2.0d0*CE(5,ne1)
                ELSE !within a tube branch
                  CE(5,ne2)=CE(5,ne2)+CE(5,ne1)
                ENDIF
              ENDIF
            ENDDO !noelem
          ENDIF
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            NJE(ne)=NJT
            DO nj=1,NJE(ne)
              NBJ(nj,ne)=1
            ENDDO !nj
            IF(DABS(CE(5,1)).GT.1.0D-6) THEN
              CE(11,ne)=CE(4,ne)/CE(5,1)
            ELSE
              CE(11,ne)=CE(4,ne)
              WRITE(OP_STRING,'('' WARNING... TOT_INIT_VOL=0'')') 
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF !CE(5,1)
          ENDDO !noelem
          DO nonode=1,NPNODE(0,nr)
            np=nonode
            NPNODE(nonode,nr)=np
            DO nj=1,NJT
              NKJ(nj,np)=1
              DO nc=1,NCM
                NVJP(nj,np)=1
              ENDDO !nc
            ENDDO !nj
          ENDDO !nonode (np)
          DO nb=1,NBFT
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              NCO(ne)=1
              DO ns=1,NST(nb)+NAT(nb)
                SE(ns,nb,ne)=1.0d0
              ENDDO !ns
              DO nn=1,NNT(nb)
                DO nk=1,NKT(nn,nb)
                  NKE(nk,nn,nb,ne)=NK
                ENDDO !nk
                DO nj=1,NJE(ne)
                  NVJE(nn,nb,nj,ne)=1 !version one of nn,nj in elem ne
                ENDDO !nj
              ENDDO !nn
            ENDDO !noelem (ne)
          ENDDO !nb
          IF(RESPIRATORY_TYPE.EQ.2)THEN
            DO cnt=1,NTB
              BBM(11,cnt,nr)=BBM_VOL/CE(5,1)
            ENDDO !cnt
          ENDIF !respiratory_type
        ELSE IF(CONDUCTING_TYPE.EQ.2)THEN !volume-filling CA
          IF(INDEX_IP2.EQ.1)THEN
            DO nreg=1,nr_high !calculate the volumes below each branch
              nr=INDEX_REG(nreg)
              DO noelem=NEELEM(0,nr),1,-1
                ne=NEELEM(noelem,nr)
                ne1=NXI(-1,1,ne) !parent element #
                IF(ne1.NE.0) CE(5,ne1)=CE(5,ne1)+CE(5,ne)
              ENDDO !noelem (ne)
            ENDDO !nreg (nr)
            CALL ASSERT(CE(5,1).GT.0.0d0,'>>Total inital volume is'
     '        //' zero',ERROR,*9999)
            DO nreg=1,nr_high !calculate the proportion of volume
              nr=INDEX_REG(nreg)
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                IF(CE(5,1).NE.0.0d0) CE(11,ne)=CE(4,ne)/CE(5,1)
              ENDDO !noelem (ne)
            ENDDO !nreg (nr)
          ELSE IF(INDEX_IP2.EQ.2)THEN
            DO nreg=1,nr_high
              nr=INDEX_REG(nreg) !lobe region
              DO noelem=NEELEM(0,nr),1,-1
                ne=NEELEM(noelem,nr)
                ne1=NXI(-1,1,ne) !parent element #
                IF(ne1.NE.0) CE(5,ne1)=CE(5,ne1)+CE(5,ne)
              ENDDO !noelem (ne)
            ENDDO
            CALL ASSERT(CE(5,TRAC).GT.0.0d0,'>>Total inital volume is'
     '        //' zero',ERROR,*9999)
            DO nreg=1,6
              nr=INDEX_REG(nreg) !lobe region
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                CE(11,ne)=CE(4,ne)/CE(5,TRAC)
              ENDDO !noelem (ne)
            ENDDO !nreg
          ENDIF
          IF(RESPIRATORY_TYPE.EQ.2)THEN !calculate BBM proportions
            DO cnt=1,NTB
              BBM(11,cnt,nr)=BBM_VOL/CE(5,1)
            ENDDO !cnt
          ENDIF !RESPIRATORY_TYPE
        ENDIF !conducting_type
      ENDIF !model_type
      IF(INDEX_IP2.NE.2)THEN
        CALL INTERFACE(NP_INTERFACE,NPNODE,ERROR,*9999)
        CALL CALC_NENP(NBJ,NEELEM,NENP,NPNE,NPNODE,ERROR,*9999)
      ENDIF

      CALL EXITS('IPMESH2')
      RETURN
 9999 CALL ERRORS('IPMESH2',ERROR)
      CALL EXITS('IPMESH2')
      RETURN 1
      END


      SUBROUTINE IPMESH9(IBT,IDO,INP,NBJ,NEELEM,NENP,
     '  NKJE,NKJ,NPF,NPNE,NPNODE,nr,NRE,NVJE,NVJP,NXI,
     '  PLSG_EDGE_LIST,PLSG_EDGE_MARKER_LIST,PLSG_POINT_MARKER_LIST,
     '  PLSG_SEGMENT_LIST,PLSG_SEGMENT_MARKER_LIST,
     '  PLSG_ELEMENT_AREA,PLSG_HOLE_LIST,PLSG_PARTITION_LIST,
     ,  SE,SP,XA,XE,XP,ZA,ERROR,*)

C#### Subroutine: IPMESH9
C###  Description:
C###    IPMESH9 wraps IPMESH9_DYNAM dynamically allocating and freeing
C###    the arrays required in IPMESH9_DYNAM

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:loc00.inc'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),NKJE(NKM,NNM,NJM,NEM),
     '  NKJ(NJM,NPM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,NRE(NEM),NVJE(NNM,NBFM,NJM,NEM),
     '  NVJP(NJM,NPM),NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  PLSG_EDGE_LIST(NLM*2),
     '  PLSG_EDGE_MARKER_LIST(NLM),PLSG_POINT_MARKER_LIST(NPM),
     '  PLSG_SEGMENT_LIST(NLM*2),PLSG_SEGMENT_MARKER_LIST(NLM)
      REAL*8 PLSG_ELEMENT_AREA(NEM),
     '  PLSG_HOLE_LIST(NRM*2),PLSG_PARTITION_LIST(NEM*2),
     '  SE(NSM,NBFM,NEM),SP(NKM,NBFM,NPM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XP(NKM,NVM,NJM,NPM),ZA(NAM,NHM,NCM,NEM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER*4 NUM_NODE_VERSIONS_PTR,POINT_ATTRIBUTE_LIST_PTR,
     '  TRIANGLE_ATTRIBUTE_LIST_PTR

      CALL ENTERS('IPMESH9',*9999)

      NUM_NODE_VERSIONS_PTR=0
      POINT_ATTRIBUTE_LIST_PTR=0
      TRIANGLE_ATTRIBUTE_LIST_PTR=0

      CALL ALLOCATE_MEMORY(NPM,1,1,NUM_NODE_VERSIONS_PTR,.TRUE.,
     '  ERROR,*9999)
      CALL ALLOCATE_MEMORY(NPM*NAM,1,3,POINT_ATTRIBUTE_LIST_PTR,.TRUE.,
     '  ERROR,*9999)
      CALL ALLOCATE_MEMORY(NEM*NAM,1,3,TRIANGLE_ATTRIBUTE_LIST_PTR,
     '  .TRUE.,ERROR,*9999)

      CALL IPMESH9_DYNAM(IBT,IDO,INP,NBJ,NEELEM,NENP,
     '  NKJE,NKJ,NPF,NPNE,NPNODE,nr,NRE,
     '  %VAL(NUM_NODE_VERSIONS_PTR),NVJE,
     '  NVJP,NXI,PLSG_EDGE_LIST,PLSG_EDGE_MARKER_LIST,
     '  PLSG_POINT_MARKER_LIST,PLSG_SEGMENT_LIST,
     '  PLSG_SEGMENT_MARKER_LIST,%VAL(POINT_ATTRIBUTE_LIST_PTR),
     '  PLSG_ELEMENT_AREA,PLSG_HOLE_LIST,PLSG_PARTITION_LIST,
     '  %VAL(TRIANGLE_ATTRIBUTE_LIST_PTR),
     '  SE,SP,XA,XE,XP,ZA,ERROR,*9999)

      CALL FREE_MEMORY(NUM_NODE_VERSIONS_PTR,ERROR,*9999)
      CALL FREE_MEMORY(POINT_ATTRIBUTE_LIST_PTR,ERROR,*9999)
      CALL FREE_MEMORY(TRIANGLE_ATTRIBUTE_LIST_PTR,ERROR,*9999)

      CALL EXITS('IPMESH9')
      RETURN
 9999 CALL ERRORS('IPMESH9',ERROR)
      CALL EXITS('IPMESH9')
      RETURN 1
      END

      SUBROUTINE IPMESH9_DYNAM(IBT,IDO,INP,NBJ,NEELEM,NENP,
     '  NKJE,NKJ,NPF,NPNE,NPNODE,nr,NRE,NUM_NODE_VERSIONS,NVJE,
     '  NVJP,NXI,PLSG_EDGE_LIST,PLSG_EDGE_MARKER_LIST,
     '  PLSG_POINT_MARKER_LIST,PLSG_SEGMENT_LIST,
     '  PLSG_SEGMENT_MARKER_LIST,POINT_ATTRIBUTE_LIST,
     '  PLSG_ELEMENT_AREA,PLSG_HOLE_LIST,PLSG_PARTITION_LIST,
     '  TRIANGLE_ATTRIBUTE_LIST,
     '  SE,SP,XA,XE,XP,ZA,ERROR,*)

C#### Subroutine: IPMESH9_DYNAM
C###  Description:
C###    IPMESH9_DYNAM defines mesh parameters for a mesh of 
C###    triangular elements, of an arbitrary 2D domain. 

C#### Variable: PLSG_EDGE_LIST()
C###  Type: INTEGER
C###  Set_up: IPMESH9_DYNAM
C###  Description: 
C###    PLSG_EDGE_LIST is an array of edges in a Delaunay 
C###    Triangulation. The first edges end points are at indices (1)
C###    and (2), followed by the remaining edges.

C#### Variable: PLSG_EDGE_MARKER_LIST()
C###  Type: INTEGER
C###  Set_up: IPMESH9_DYNAM
C###  Description: 
C###    PLSG_EDGE_MARKER_LIST is an array of edge markers.

C#### Variable: PLSG_POINT_MARKER_LIST(np)
C###  Type: INTEGER
C###  Set_up: IPMESH9_DYNAM
C###  Description: 
C###    PLSG_POINT_MARKER_LIST is an array of point/node markers.

C#### Variable: PLSG_SEGMENT_LIST()
C###  Type: INTEGER
C###  Set_up: IPMESH9_DYNAM
C###  Description: 
C###    PLSG_SEGMENT_LIST is an array of segment endpoints in a 
C###    Planar Straight Line Graph. The first segment's endpoints are
C###    at indices (1) and (2) followed by the remaining segments.

C#### Variable: PLSG_SEGMENT_MARKER_LIST()
C###  Type: INTEGER
C###  Set_up: IPMESH9_DYNAM
C###  Description: 
C###    PLSG_SEGMENT_MARKER_LIST is an array of segment markers.
C###    For use in Planar Straight Line Graph definition.

C#### Variable: PLSG_HOLE_LIST()
C###  Type: REAL*8
C###  Set_up: IPMESH9_DYNAM
C###  Description: 
C###    PLSG_HOLE_LIST is an array of holes in a Planar Straight Line
C###    Graph. The first hole's x and y coordinates are
C###    at indices (1) and (2) followed by the remaining holes.

C#### Variable: PLSG_PARTITION_LIST() 
C###  Type: REAL*8
C###  Set_up: IPMESH9_DYNAM
C###  Description: 
C###    PLSG_PARTITION_LIST an array of partitions in a Planar 
C###    Straight Line Graph. Refered to as 'regions' in Triangle V1.3
C###    The first partion's x and y coordinates are at indices (1) and
C###    (2) followed by an attribute at index (3) and the maximum
C###    triangle area at index (4), followed by the remaining 
C###    partitions. 

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:call00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
      INCLUDE 'cmiss$reference:loc00.inc'
      INCLUDE 'cmiss$reference:mesh00.cmn'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM), 
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),NKJE(NKM,NNM,NJM,NEM),
     '  NKJ(NJM,NPM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,
     '  NRE(NEM),NUM_NODE_VERSIONS(NPM),NVJE(NNM,NBFM,NJM,NEM),
     '  NVJP(NJM,NPM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),PLSG_EDGE_LIST(NLM*2),
     '  PLSG_EDGE_MARKER_LIST(NLM),PLSG_POINT_MARKER_LIST(NPM),
     '  PLSG_SEGMENT_LIST(NLM*2),PLSG_SEGMENT_MARKER_LIST(NLM)
      REAL*8 PLSG_ELEMENT_AREA(NEM),PLSG_HOLE_LIST(NRM*2),
     '  PLSG_PARTITION_LIST(NEM*2),POINT_ATTRIBUTE_LIST(NPM*NAM),
     '  SE(NSM,NBFM,NEM),SP(NKM,NBFM,NPM),
     '  TRIANGLE_ATTRIBUTE_LIST(NEM*NAM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XP(NKM,NVM,NJM,NPM),
     '  ZA(NAM,NHM,NCM,NEM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ERROR_FLAG,ERROR_STRINGC(500),
      !SMAR009 22/12/98 GENER_TYPE,
     '  i,IBEG1,IBEG2,IEND1,IEND2,ICHAR,INFO,
     '  j,LENGTH,METHOD,nb,nb1,ne,nj,nk,nn,ns,np,noelem,
     '  nonode,NOQUES,NUM_POINTS,
     '  NUM_TRIANGLE_ATTRIBUTES,TAKEOFF
      CHARACTER CHAR1*20,CHAR2*5,ERROR_STRINGF*100
      LOGICAL FILEIP

C CS 23/7/98 removing RB's obsolete vars
C      INTEGER GENER_TYPE,IBEG3,IEND3,LOOP,NUM_CRNRS,
C      CHARACTER CHAR3*5,
C      REAL*8 PERTURB

      CALL ENTERS('IPMESH9_DYNAM',*9999)

      CALL ASSERT(USE_TRIANGLE.EQ.1,
     '  '>>USE_TRIANGLE must be set to 1',ERROR,*9999)
      CALL ASSERT(NJT.EQ.2,
     '  '>>The number of global coordinates must be 2',ERROR,*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

C     Flags follow the sequence prqa_aAcevngBPNEIOXzoYSiFlsCQVVVhf
      DO i=1,40
        FLAGS(i)=0
      ENDDO
      FLAGS(30)=1  !Q
  
      IDEFLT(1)=1
      FORMAT='($,'' Enter basis function # for mesh [1]: '',I1)'
      IF(IOTYPE.EQ.3) IDATA(1)=IDEFLT(1)
      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NBT,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) nb=IDATA(1)

      FORMAT='('' Enter triangulation input type [1]:'''//
     '  '/''   (1) Planar Straight Line Graph'''//
     '  '/''   (2) Set of nodes'''//
     '  '/$,''    '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) IDATA(1)=IDEFLT(1)
      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) INPUT_TYPE=IDATA(1)

      FORMAT='($,'' Do you wish to generate boundary/internal '
     '  //'boundary nodes [N]? '',A)'      
      IF(IOTYPE.EQ.3) ADATA(1)='N'
      CALL GINOUT(IOTYPE,1,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        IF(ADATA(1).EQ.'Y') THEN
C rgb 25/6/98 obsolete
C          FLAGS(36)=1 ! I - ignore type 1 of in triangulation
C          INPUT_TYPE=2
C          FORMAT='(''  Enter how B/IB nodes are to be made up[1]:'''//
C     '      '/''    (1) From elements'''//
C     '      '/''    (2) From a set of nodes'''//
C     '      '/$,''     '',I1)'      
C          IDEFLT(1)=1
C          IF(IOTYPE.EQ.3) IDATA(1)=IDEFLT(1)
C          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C          IF(IOTYPE.NE.3) GENER_TYPE=IDATA(1)
C          RDEFLT(1)=0.01d0
C          FORMAT='($,''  Enter the distance between the B/IB node'
C     '      //'pairs [0.01]: '',F12.4)'
C          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
C          IF(IOTYPE.NE.3) PERTURB=RDATA(1)
C          IDEFLT(1)=0
C          FORMAT='($,''  Enter the number of corners in the mesh'
C     '      //' [0]: '',I3)'
C          IF(IOTYPE.EQ.3) IDATA(1)=IDEFLT(1)
C          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C          IF(IOTYPE.NE.3) NUM_CRNRS=IDATA(1)
C          DO i=1,NUM_CRNRS
C            IDEFLT(1)=i
C            WRITE(CHAR2,'(I3)') i
C            CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
C            FORMAT='($,''   The corner node number is ['//CHAR2(IBEG2:
C     '        IEND2)//']: '',I3)'
C            IF(IOTYPE.EQ.3) IDATA(1)=IDEFLT(1)
C            CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
C     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C            IF(IOTYPE.NE.3) CRNRS(0,i)=IDATA(1)  
C            IDEFLT(1)=1
C            IF(NJT.EQ.2) THEN
C              WRITE(CHAR2,'(I3)') i
C              CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
C              FORMAT='(''   Corner node '//CHAR2(IBEG2:IEND2)//
C     '          ' is [1]:'''//
C     '          '/''    (1) Convex  '''//
C     '          '/''    (2) Concave'''//
C     '          '/$,''     '',I1)'
C              IF(IOTYPE.EQ.3) IDATA(1)=IDEFLT(1)
C              CALL GINOUT(IOTYPE,3,IVDU,IFILE,1,0,NOQUES,FILEIP,
C     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
C     '          1,6,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
C     '          *9999)
C              IF(IOTYPE.NE.3) CRNRS(-1,i)=IDATA(1)
C              WRITE(CHAR2,'(I1)') NJT
C              WRITE(CHAR3,'(I3)') i
C              CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
C              CALL STRING_TRIM(CHAR3,IBEG3,IEND3)
C              FORMAT='($,''   Enter the other '//CHAR2(IBEG2:IEND2)//
C     '          ' numbers of corner node  '//CHAR3(IBEG3:IEND3)//
C     '          ' : '',18(1X,I5))'
C              IF(IOTYPE.NE.3) THEN
C                CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
C     '            FORMAT,NJT,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
C     '            IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
C     '            INFO,ERROR,*9999)
C                DO j=1,NJT
C                  CRNRS(j,i)=IDATA(j)
C                ENDDO
C              ENDIF        
C            ELSE IF(NJT.EQ.3) THEN
C              WRITE(CHAR2,'(I3)') i
C              CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
C              IDEFLT(1)=1
C              FORMAT='(''   Corner node '//CHAR2(IBEG2:IEND2)//
C     '          ' is [1]:'''//
C     '          '/''     (1) Convex edge   - made up of three nodes'''//
C     '          '/''     (2) Concave edge  - made up of three nodes'''//
C     '          '/''     (3) Convex corner  - made up of four nodes'''//
C     '          '/''     (4) Concave corner - made up of four nodes'''//
C     '          '/$,''     '',I1)'
C              IF(IOTYPE.EQ.3) IDATA(1)=IDEFLT(1)
C              CALL GINOUT(IOTYPE,3,IVDU,IFILE,1,0,NOQUES,FILEIP,
C     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
C     '          1,6,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
C     '          *9999)
C              IF(IOTYPE.NE.3) CRNRS(0,i)=IDATA(1)
C              IF((CRNRS(0,i).EQ.1).OR.(CRNRS(0,i).EQ.2)) LOOP=2
C              IF((CRNRS(0,i).EQ.3).OR.(CRNRS(0,i).EQ.4)) LOOP=3
C              WRITE(CHAR2,'(I1)') LOOP
C              WRITE(CHAR3,'(I3)') i
C              CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
C              CALL STRING_TRIM(CHAR3,IBEG3,IEND3)
C              FORMAT='($,''  Enter the other '//CHAR2(IBEG2:IEND2)//
C     '          ' numbers of corner node in counter clockwise order
C     '          '//CHAR3(IBEG3:IEND3)//
C     '          ' : '',18(1X,I5))'
C              IF(IOTYPE.NE.3) THEN
C                CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
C     '            FORMAT,NJT,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
C     '            IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
C     '            INFO,ERROR,*9999)
C                DO j=1,LOOP
C                  CRNRS(j,i)=IDATA(j)
C                ENDDO
C              ENDIF
C            ENDIF
C          ENDDO ! i
C          DO np=1,NPNODE(0,nr)
C            NUM_POINTS=NUM_POINTS+(NUM_NODE_VERSIONS(np)-1)
C          ENDDO
        ELSE
          DO i=1,NPNODE(0,nr)
            NUM_NODE_VERSIONS(i) = 1
          ENDDO
          NUM_POINTS=NPNODE(0,nr)
        ENDIF
      ENDIF 

      IF(INPUT_TYPE.EQ.1) THEN 
        FLAGS(1)=1   !p
        IDEFLT(1)=NPNODE(0,nr)
        WRITE(CHAR1,'(I3)') IDEFLT(1)
        CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
        FORMAT='($,'' Enter the number of segments ['//CHAR1(IBEG1
     '      :IEND1)//']: '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=IDEFLT(1)
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NUM_SEGMENTS=IDATA(1)

        j=1
        DO i=1,NUM_SEGMENTS
          WRITE(CHAR1,'(I3)') i
          CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
          IF(i.NE.NUM_SEGMENTS) THEN
            WRITE(CHAR2,'(I3)') i+1
          IDEFLT(2)=i+1
          ELSE
            WRITE(CHAR2,'(I3)') 1
          IDEFLT(2)=1
          ENDIF
          IDEFLT(1)=i
          CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
          FORMAT='($,'' Enter the end nodes of segment '//CHAR1(IBEG1
     '      :IEND1)//' ['//CHAR1(IBEG1:IEND1)//','//CHAR2(IBEG2:
     '      IEND2)//']: '',2I3)'
          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,2,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
           PLSG_SEGMENT_LIST(j)=IDATA(1)
           PLSG_SEGMENT_LIST(j+1)=IDATA(2)
           j=j+2
          ENDIF
          PLSG_SEGMENT_MARKER_LIST(i)=0
        ENDDO !i        

        IDEFLT(1)=0
        FORMAT='($,'' Enter the number of holes [0]: '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=IDEFLT(1)
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NUM_HOLES=IDATA(1)

        j=1
        DO i=1,NUM_HOLES
          WRITE(CHAR1,'(I3)') i
          CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
          RDEFLT(1)=0.0d0
          RDEFLT(2)=0.0d0
          FORMAT='($,'' Enter the coordinates of a point in hole '//
     '      CHAR1(IBEG1:IEND1)//' [0,0]: '',2F12.4)'
          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,2,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            PLSG_HOLE_LIST(j)=RDATA(1)
            PLSG_HOLE_LIST(j+1)=RDATA(2)
            j=j+2
          ENDIF
        ENDDO !i

        IDEFLT(1)=0
        FORMAT='($,'' Enter the number of local area constraints'
     '    //' [0]: '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=IDEFLT(1)
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NUM_PARTITIONS=IDATA(1)
        IF(NUM_PARTITIONS.NE.0) THEN     
          FLAGS(6)=1 !a (partitian area constraints)
        ENDIF
 
        j=1
        DO i=1,NUM_PARTITIONS
          WRITE(CHAR1,'(I3)') i
          CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
          RDEFLT(1)=0.0d0
          FORMAT='($,'' Enter the x coordinate of a point in '
     '      //'partitian '//CHAR1(IBEG1:IEND1)//' [0.0]: '',F12.4)'
          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) PLSG_PARTITION_LIST(j)=RDATA(1)
          RDEFLT(1)=0.0d0
          FORMAT='($,'' Enter the y coordinate of a point in '
     '      //'partitian '//CHAR1(IBEG1:IEND1)//' [0.0]: '',F12.4)'
          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) PLSG_PARTITION_LIST(j+1)=RDATA(1)
          RDEFLT(1)=1.0d0
          FORMAT='($,'' Enter the maximum triangle area in partitian '//
     '      CHAR1(IBEG1:IEND1)//'[1]: '',F12.4)'
          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
           PLSG_PARTITION_LIST(j+2)=0.0d0
           PLSG_PARTITION_LIST(j+3)=RDATA(1)
           j=j+4
          ENDIF
        ENDDO !i

        FORMAT='($,'' Allow boundary segment splitting [Y]? '',A)'
        IF(IOTYPE.EQ.3) ADATA(1)='Y'
        CALL GINOUT(IOTYPE,1,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,AYES,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          IF(ADATA(1).EQ.'N') THEN
            FLAGS(22)=1 !Y
          ENDIF
        ENDIF

        FORMAT='($,'' Allow internal segment splitting [Y]? '',A)'
        IF(IOTYPE.EQ.3) ADATA(1)='Y'
        CALL GINOUT(IOTYPE,1,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,AYES,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          IF(ADATA(1).EQ.'N') THEN
            FLAGS(22)=1 !Y
          ENDIF
        ENDIF

      ELSE
        NUM_HOLES=0
        NUM_PARTITIONS=0
      ENDIF

      RDEFLT(1)=RMAX
      FORMAT='($,'' Enter the maximum triangle area  [NO_MAXIMUM]: '',
     '  F12.4)'
      IF(IOTYPE.EQ.3) RDATA(1)=RDEFLT(1)
      CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) MAX_AREA=RDATA(1)
      IF(MAX_AREA.LT.RMAX) THEN     
        FLAGS(4)=1 !a
        FLAGS(5)=1 !fixed
      ENDIF

      FORMAT='($,'' Enter the minimum angle size in '
     '  //'triangulation [20.0]: '',F12.4)'
      RDEFLT(1)=20.0d0
      IF(IOTYPE.EQ.3) RDATA(1)=RDEFLT(1)
      CALL GINOUT(IOTYPE,5,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        MIN_ANGLE=RDATA(1)
        FLAGS(3)=1 !q
      ENDIF

      FORMAT='($,'' Enter the maximum number of Steiner '
     '  //'points [INF]? '',I10)'
      IDEFLT(1)=IMAX
      IF(IOTYPE.EQ.3) IDATA(1)=IDEFLT(1)
      CALL GINOUT(IOTYPE,3,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        IF(IDATA(1).LT.IMAX) THEN
          NUM_STEINER_POINTS=IDATA(1)
          FLAGS(23)=1 !S
        ENDIF
      ENDIF

      FORMAT='($,'' Set up Voronoi information [N]? '',A)'
      IF(IOTYPE.EQ.3) ADATA(1)='N'
      CALL GINOUT(IOTYPE,1,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        IF(ADATA(1).EQ.'Y') THEN
          FLAGS(10)=1 !v
          FLAGS(11)=1 !n
        ENDIF
      ENDIF

      FORMAT='('' Enter Delauney Triangulation method [1]:'''//
     '  '/''   (1) Divide-and-Conquer'''//
     '  '/''   (2) Divide-and-Conquer with vertical cuts only'''//
     '  '/''   (3) Incremental'''//
     '  '/''   (4) Sweepline'''//
     '  '/$,''    '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) IDATA(1)=1
      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) METHOD=IDATA(1)
      IF(METHOD.EQ.2) THEN
        FLAGS(27)=1  !l
      ENDIF
      IF(METHOD.EQ.3) THEN
        FLAGS(25)=1  !i
      ENDIF
      IF(METHOD.EQ.4) THEN
        FLAGS(26)=1  !F
      ENDIF
 
      IF(NNT(nb).EQ.6) FLAGS(21)=1  !o2
      ERROR_FLAG=0
      NUM_TRIANGLE_ATTRIBUTES=0
      NUM_EDGES=0

      CALL triangulate(AREA_FIELD,ERROR_FLAG,FLAGS,IBT,
     '  IDO,INP,NAM,NBFM,NBJ,NCM,
     '  NEELEM,NEM,NHM,NIM,NJM,NKJE,NKM,NLM,NNM,NPF,NPM,NPNE,nr,
     '  NNT(nb),NUM_EDGES,NUM_HOLES,NPNODE(0,nr),NUM_NODE_VERSIONS,
     '  NUM_PARTITIONS,
     '  NUM_POINTS,0,NUM_SEGMENTS,
     '  NUM_STEINER_POINTS,NEELEM(0,nr),
     '  NUM_TRIANGLE_ATTRIBUTES,NVJE,NVJP,NVM,NXI,PLSG_EDGE_LIST,
     '  PLSG_EDGE_MARKER_LIST,PLSG_POINT_MARKER_LIST,PLSG_SEGMENT_LIST,
     '  PLSG_SEGMENT_MARKER_LIST,MAX_AREA,MIN_ANGLE,
     '  PLSG_ELEMENT_AREA,PLSG_HOLE_LIST,
     '  PLSG_PARTITION_LIST,POINT_ATTRIBUTE_LIST,SE,
     '  TRIANGLE_ATTRIBUTE_LIST,XA,XE,XP,ZA,ERROR_STRINGC)

      IF(ERROR_FLAG.NE.0) THEN
        CALL CSTRINGLEN(LENGTH,ERROR_STRINGC)
        CALL C2FSTRING(ERROR_STRINGC,LENGTH,ERROR_STRINGF)
        CALL ASSERT(.FALSE.,ERROR_STRINGF,ERROR,*9999)
      ENDIF

!     This gives the original number of nodes plus the 
!     number of new points created
      IF(FLAGS(36).EQ.1) THEN
        TAKEOFF=2
      ELSE
        TAKEOFF=1
      ENDIF
      DO np=1,NPNODE(0,nr)
        NUM_POINTS=NUM_POINTS-(NUM_NODE_VERSIONS(np)-TAKEOFF)
      ENDDO

      NET(nr)=NEELEM(0,nr)       !highest element# in region nr
      NET(0) =NEELEM(0,nr)       !highest element# in all regions
      NEELEM(0,0) =NEELEM(0,nr)  !number elements in all regions
      NPT(nr)=NUM_POINTS         !highest node# in region nr
      NPT(0) =NUM_POINTS         !highest node# in all regions
      CALL ASSERT(NPT(nr).LE.NP_R_M,'>>Increase NP_R_M',ERROR,*9999)
      CALL ASSERT(NET(nr).LE.NE_R_M,'>>Increase NE_R_M',ERROR,*9999)

      IF(NBT.GT.1) THEN
        DO nb1=2,NBT
          DO ne=1,NET(nr)
            DO nn=1,NNT(nb1)
              NPNE(nn,nb1,ne)=NPNE(nn,1,ne)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      DO noelem=1,NEELEM(0,nr)
        ne=noelem
        NEELEM(noelem,nr)=NE
        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
          NBJ(nj,ne)=nb
          DO nn=1,NNT(nb)
            DO nk=1,NKT(nn,nb)
              NKJE(nk,nn,nj,ne)=nk
            ENDDO !nk
          ENDDO !nn
        ENDDO
        NRE(ne)=nr
      ENDDO

      DO nonode=1,NUM_POINTS
        np=nonode
        NPNODE(nonode,nr)=np
        DO nj=1,NJT
          NKJ(nj,np)=NKT(0,nb)
          NVJP(nj,np)=1 !One version for all
        ENDDO !nj                 !extra nodes created
      ENDDO !nonode (np)
          
      NPNODE(0,nr)=NUM_POINTS
      NPNODE(0,0) =NUM_POINTS

      DO nb1=1,NBFT
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          DO ns=1,NST(nb1)+NAT(nb1)
            SE(ns,nb1,ne)=1.0D0
          ENDDO !ns
C KAT 23Feb01: now handled by NKJE above
C          DO nn=1,NNT(nb1)
C            DO nk=1,NKT(nn,nb1)
C              NKE(nk,nn,nb1,ne)=NK
C            ENDDO !nk
C          ENDDO !nn
        ENDDO !noelem (ne)
      ENDDO !nb1

      CALL CALC_NENP(NBJ,NEELEM,NENP,NPNE,ERROR,*9999)

      CALL_TRIANGLE=.TRUE.

      CALL EXITS('IPMESH9_DYNAM')
      RETURN
 9999 CALL ERRORS('IPMESH9_DYNAM',ERROR)
      CALL EXITS('IPMESH9_DYNAM')
      RETURN 1
      END


Module FE14
=========== 

      SUBROUTINE CELL_ARRAY(iw,NDIMX,NDIMY,XMIN_CA,XMAX_CA,YMIN_CA,
     '  YMAX_CA,ERROR,*)

C#### Subroutine: CELL_ARRAY
C###  Description:
C###    CELL_ARRAY draws a cell array given by array I2P(512*512,iw) 
C###    on w/s IW(1,2).  If NDIMX,Y are not equal to 512, the program 
C###    only draws every 512/NDIMX,Yth pixel.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:disp00.cmn'
      INCLUDE 'gx$path:gx.inc'
      INCLUDE 'cmiss$reference:pics00.cmn'
      INCLUDE 'cmiss$reference:pics01.cmn'
!     Parameter List
      INTEGER iw,NDIMX,NDIMY
      REAL XMAX_CA,XMIN_CA,YMAX_CA,YMIN_CA
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ERR,i,j,NSTEPX,NSTEPY
      REAL BNDRECT(4)

      CALL ENTERS('CELL_ARRAY',*9999)
      NSTEPX=NINT(DBLE(IMGX)/DBLE(NDIMX))
      NSTEPY=NINT(DBLE(IMGY)/DBLE(NDIMY))
      IF(DOP) THEN
        WRITE(OP_STRING,'('' Stepping by '',2I4)') NSTEPX,NSTEPY
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      IF(iw.EQ.1.OR.iw.EQ.2) THEN
        IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
          DO j=1,NDIMY
            DO i=1,NDIMX
              I2P(i,j,iw+2)=INT(DBLE(I2P(i*NSTEPX,j*NSTEPY,iw))
     '          *(gxNSPECT)/256.d0)+gxSPOFF
            ENDDO
          ENDDO
          BNDRECT(1)=XMIN_CA
          BNDRECT(2)=YMIN_CA
          BNDRECT(3)=XMAX_CA
          BNDRECT(4)=YMAX_CA
          CALL DRWIMG(BNDRECT,NDIMX,NDIMY,I2P(1,1,iw+2),ERR)
          IF(ERR.GT.0) THEN
            ERROR='>>Error from GX'
            GOTO 9999
          ENDIF
        ENDIF
      ENDIF

      CALL EXITS('CELL_ARRAY')
      RETURN
 9999 CALL ERRORS('CELL_ARRAY',ERROR)
      CALL EXITS('CELL_ARRAY')
      RETURN 1
      END


      SUBROUTINE CREATE_SEGMENT(ISEGNUM,iw,ERROR,*)

C#### Subroutine: CREATE_SEGMENT
C###  Description:
C###    CREATE_SEGMENT creates graphics segment ISEGNUM.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:disp00.cmn'
!     Parameter List
      INTEGER ISEGNUM,iw
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER ERR

      CALL ENTERS('CREATE_SEGMENT',*9999)
      IF(IWKT(iw).EQ.1) THEN      !GKS
        IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
          CALL OPOBJT(ISEGNUM,ERR)
          IF(err.ne.0) THEN
            ERROR='>>Error from GX'
            GOTO 9999
          ENDIF
        ENDIF
      ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS
      ENDIF

      CALL EXITS('CREATE_SEGMENT')
      RETURN
 9999 CALL ERRORS('CREATE_SEGMENT',ERROR)
      CALL EXITS('CREATE_SEGMENT')
      RETURN 1
      END


      SUBROUTINE ELLIPSE(SEMIAXIS1,SEMIAXIS2,XWC,YWC,NLINES,ERROR,*)

C#### Subroutine: ELLIPSE
C###  Description:
C###    ELLIPSE draws fill-area ellipse about the point XWC,YWC with
C###    semi-axis lengths SEMIAXIS1 and SEMIAXIS2 with NLINES lines
C###    (max 100).

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:disp00.cmn'
!     Parameter List
      INTEGER NLINES
      REAL*8 SEMIAXIS1,SEMIAXIS2,XWC,YWC
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,ERR
      REAL*8 D_ELLIPSE(3,100),THETA
      REAL R_PTS(100,3)

      CALL ENTERS('ELLIPSE',*9999)
      DO i=1,NLINES
        THETA=2.0D0*PI*DBLE(i-1)/DBLE(NLINES)
        D_ELLIPSE(1,i)=XWC+SEMIAXIS1*DCOS(THETA)
        D_ELLIPSE(2,i)=YWC+SEMIAXIS2*DSIN(THETA)
        R_PTS(i,1)=REAL(D_ELLIPSE(1,i))   !is real x
        R_PTS(i,2)=REAL(D_ELLIPSE(2,i))   !is real y
      ENDDO
      IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
        CALL FLAREA(NLINES,R_PTS(1,1),R_PTS(1,2),ERR)
        IF(err.ne.0) THEN
          ERROR='>>Error from GX'
          GOTO 9999
        ENDIF
      ENDIF

      CALL EXITS('ELLIPSE')
      RETURN
 9999 CALL ERRORS('ELLIPSE',ERROR)
      CALL EXITS('ELLIPSE')
      RETURN 1
      END


      SUBROUTINE PRINT_IMAGE_FILL_AREAS(NOCO,NTCO,CO,NTCOQU,COQU,
     '  ERROR,*)

C**** Print fill area representation of image

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:pics00.cmn'
      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
!     Parameter List
      INTEGER NOCO,NTCO,NTCOQU(*)
      CHARACTER CO(*)*(*),COQU(16,*)*(*),ERROR*(*)
!     Local Variables
      INTEGER IBEG,IDUMMY,IEND,IERST,IFROMC,I_GREY,IP,IW,IX,IY,N3CO,NTIW
      REAL DISP,DISPX,DISPY,SIZEX,SIZEX_MIN,SIZEY,SIZEY_MIN,
     '  XMAX,XMIN,YMAX,YMIN
      REAL*8 PTS(3,5)
      LOGICAL CBBREV

C     CALL ENTERS('PRINT_IMAGE_FILL_AREAS',*9999)
      IF(CO(NOCO).EQ.'?') THEN
        WRITE(OP_STRING,*) 'print I2P* <on WS>[on 15]'
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ELSE IF(CO(NOCO).EQ.'??') THEN
        CALL DOCUM('imp','doc','PRINT_IMAGE',ERROR,*9999)
      ELSE
        CALL TRIM(CO(NOCO),IBEG,IEND)
        IP=IFROMC(CO(NOCO)(IEND:IEND))
        IF(CBBREV(CO,'ON',1,NOCO+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),1,NTIW,IW,ERROR,*9999)
        ELSE
          IW=15
        ENDIF
        CALL ACWK(IW,0,ERROR,*9999)
        IF(IW.EQ.1) THEN
          XMIN=1.0
          YMIN=1.0
          XMAX=65.0
          YMAX=65.0
        ELSE IF(IW.EQ.15) THEN
          XMIN=1.0
          YMIN=1.0
          XMAX=129.0
          YMAX=129.0
        ENDIF
        CALL WKST_WINDOW(IW,XMIN,XMAX,YMIN,YMAX,ERROR,*9999)
        CALL GKS_SELNT(IW,ERROR,*9999)
        CALL GKS_QDSP(61,IERST,IDUMMY,DISPX,DISPY,IDUMMY,IDUMMY,
     '    ERROR,*9999)
        DISP=MIN(DISPX,DISPY)
        SIZEX_MIN=0.0254*8.0/300.0*(XMAX-XMIN)/DISP
        SIZEY_MIN=0.0254*8.0/300.0*(YMIN-YMAX)/DISP
        SIZEX=1.0
        SIZEY=1.0
        IF(SIZEX.LT.SIZEX_MIN)SIZEX=SIZEX_MIN
        IF(SIZEY.LT.SIZEY_MIN)SIZEY=SIZEY_MIN

        CALL GKS_SPARF(1.0,1.0,ERROR,*9999)
        CALL GKS_SPA(SIZEX,SIZEY,ERROR,*9999)
        IF(DOP) THEN
          WRITE(OP_STRING,*)' SIZEX,SIZEY,SIZEX_MIN,SIZEY_MIN=',
     '    SIZEX,SIZEY,SIZEX_MIN,SIZEY_MIN
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

        DO IX=XMIN,XMAX-1
          DO IY=YMIN,YMAX-1
            I_GREY=I2P(IX,IY,IP)/16 + 1
            PTS(1,1)=DBLE(IX)
            PTS(1,2)=DBLE(IX+1)
            PTS(1,3)=DBLE(IX+1)
            PTS(1,4)=DBLE(IX)
            PTS(1,5)=DBLE(IX)
            PTS(2,1)=DBLE(IY)
            PTS(2,2)=DBLE(IY)
            PTS(2,3)=DBLE(IY+1)
            PTS(2,4)=DBLE(IY+1)
            PTS(2,5)=DBLE(IY)
            CALL GKS_SFAI(I_GREY,ERROR,*9999)
C            IF(I_GREY.EQ.0) THEN
C              CALL GKS_SFAIS(GHOLLO,ERROR,*9999)
C            ELSE IF(I_GREY.EQ.16) THEN
C              CALL GKS_SFAIS(GSOLID,ERROR,*9999)
C            ELSE
C              CALL GKS_SFAIS(GPATTR,ERROR,*9999)
C              CALL GKS_SFASI(13+I_GREY,ERROR,*9999)
C            ENDIF
            CALL GKS_FA(1,2,5,PTS,ERROR,*9999)
          ENDDO
        ENDDO
        CALL DAWK(IW,0,ERROR,*9999)

      ENDIF

 9998 CALL EXITS('PRINT_IMAGE_FILL_AREAS')
      RETURN
 9999 CALL ERRORS('PRINT_IMAGE_FILL_AREAS',ERROR)
      CALL EXITS('PRINT_IMAGE_FILL_AREAS')
      RETURN 1
      END




Module FE19
=========== 


      SUBROUTINE ADAMS_MOULTON(AII,AIO,CONTROL,ERROR_TYPE,IFAIL,IWORK,
     '  L_AII,L_AIO,L_ARI,L_ARO,L_CONTROL,L_IWORK,L_MODEL,
     '  L_PARAMETERS,L_WORK,MAX_ITERS,MAX_ORDER,MODEL,NUMBER_EQN,
     '  ABS_ERR,ARI,ARO,DY,MAX_STEP,PARAMETERS,REL_ERR,T,TOUT,WORK,Y,
     '  EXTEND_INTERVAL,STIFF_EQNS,USE_ROUND_CTRL,FUNC,ERROR)

C#### Subroutine: ADAMS_MOULTON
C###  Description:
C###    <html><pre>
C###    This procedure integrates a system of NUMBER_EQN ordinary
C###    differential equations of the form 
C###      dy(i)/dt=FUNC(t,y(1),y(2),...,y(NUMBER_EQN)) 
C###    from t=T to t=TOUT given y(i) at t=T using an Adams-Moulton
C###    integration scheme. The integrator used is adapted from the
C###    Adams-Moulton integrator of L.F. Shampine and M.K. Gordon 
C###    and is detailed in their book "Computer Solution of Ordinary
C###    Differential Equations: The Initial Value Problem".
C###    The exact form of the rhs function is
C###      CALL FUNC(T,Y,DY,AII,AIO,CONTROL,L_AII,L_AIO,L_ARI,L_ARO,
C###        L_CONTROL,L_MODEL,L_PARAMETERS,MODEL,NUMBER_EQN,
C###        ARI,ARO,PARAMETERS,ERR)
C###    where CONTROL, MODEL, ARI, ARO are integer vectors of sizes
C###    L_CONTROL, L_MODEL, L_ARI and L_ARO respectively and PARAMETERS,
C###    Y, DY, ARI, ARO are double precision vectors of sizes 
C###    L_PARAMETERS, NUMBER_EQN, NUMBER_EQN, L_ARI and L_ARO 
C###    respectively.
C###    The integrator uses a modified divided difference form of the
C###    Adams PECE formulas and local extropolation. It iterates on 
C###    the solution (up to MAX_ITERS) and adjusts the order (up to 
C###    MAX_ORDER and no more than 12) and step size (up to MAX_STEP)
C###    to control the local error per unit step in a generalised 
C###    sense. 
C###    The error control is based on the L2 Norm of the weighted 
C###    local error vector. The weighting is controlled by ERROR_TYPE.
C###    If ERROR_TYPE=1 then no weighting is used; if ERROR_TYPE=2 
C###    then the error vector components are weighted by the most 
C###    recent component of the solution vector; if ERROR_TYPE=3 then
C###    the error vector components are weighted by the most recent 
C###    component of the residual (derivative) vector; if ERROR_TYPE=4
C###    then a mixed relative weighting is used. This weighting is 
C###    calculated from (REL_ERR*y(i)+ABS_ERR)/MAX(REL_ERR,ABS_ERR).
C###    For reasons of efficiency the integrator integrates beyond
C###    TOUT internally (though never beyond T+10*(TOUT-T)) and 
C###    interpolates the solution at TOUT. If it is not possible to
C###    integrate beyond TOUT then EXTEND_INTERVAL should be set to
C###    .FALSE.
C###    The integrator can perform propagative rounding error control
C###    by setting USE_ROUND_CONTROL to .TRUE.
C###    The integrator needs workspace allocated as follows:
C###      IWORK(L_IWORK) where L_IWORK=7 and
C###      WORK(L_WORK) where L_WORK=7+7*MAX_ORDER+
C###        NUMBER_EQN*(MAX_ORDER+6) in general, and 
C###        L_WORK=7+7*MAX_ORDER+NUMBER_EQN*(MAX_ORDER+8) if rounding
C###        control is required.
C###    The workspace needs to be mainted between calls. Interesting
C###    values contain within the workspace are: 
C###      IWORK(1) - the order, K, of the Adams polynomial to be used 
C###        with the next call,
C###      IWORK(2) - the order of the Adams polynomial used in the 
C###        previous (current) step,
C###      WORK(1) - the step size, H, to be used with the next call,
C###      WORK(2) - the step size used in the previous (current) step.
C###      WORK(3) - the finish time of the last step taken.
C###      WORK(4) - the finish time of the last step output.
C###    On entry the following parameters are required:
C###      NUMBER_EQN - the number of equations to be integrated,
C###      Y - vector of initial conditions,
C###      T - starting point of the integration,
C###      TOUT - point at which the solution is desired,
C###      ABS_ERR, REL_ERR - absolute and relative local error 
C###        tolerances,
C###      ERROR_TYPE - the type of error test to be used,
C###      IFAIL - error code indicator. Must be set to 1 for start-up,
C###      MAX_ITERS - the maximum number of Adams iterations (steps)
C###        allowed to reach TOUT,
C###      MAX_ORDER - the maximum order allowed for the Adams 
C###        polynomial,
C###      MAX_STEP - the maximum step size allowed to reach TOUT,
C###      EXTEND_INTERVAL - set to .TRUE. if it is possible to 
C###        integrate beyond TOUT.
C###      USE_ROUND_CTRL - set to .TRUE. if roudning control is to be
C###        used.
C###    On output the following parameters are set:
C###      Y - vector of solutions at TOUT,
C###      DY - vector of the derivates at TOUT,
C###      T - last point reach in the integration. Normal return has
C###        T=TOUT,
C###      ABS_ERR, REL_ERR - normal return has the tolerances unchanged.
C###        If IFAIL=3 then the tolerances are set to increased values.
C###      STIFF_EQNS - set to .TRUE. if the equations appear to be
C###        stiff,
C###      IFAIL - error flag. If 
C###        IFAIL=2 - normal return, integration reached TOUT,
C###        IFAIL=3 - integration did not reach TOUT because the error
C###          are too small. ABS_ERR and REL_ERR have been increased
C###          to appropriate levels.
C###        IFAIL=4 - integration did not reach TOUT because more than
C###          MAX_ITERS steps were needed,
C###        IFAIL=5 - integration did not reach TOUT because equations
C###          appear to be stiff,
C###        IFAIL=6 - invalid input parameters,
C###        IFAIL=7 - error returned from FUNC,
C###        IFAIL=8 - if ERROR_TYPE is 2 or 3 then one of error weights
C###          is zero,
C###        IFAIL=9 - WORK array is too small,
C###        IFAIL=10 - IWORK array is too small.
C###    </pre></html>

      IMPLICIT NONE
      
!     Parameter list
      INTEGER L_AII,L_AIO,L_ARI,L_ARO,L_CONTROL,L_IWORK,L_MODEL,
     '  L_PARAMETERS,L_WORK
      INTEGER AII(L_AII),AIO(L_AIO),CONTROL(L_CONTROL),ERROR_TYPE,
     '  IFAIL,IWORK(L_IWORK),MAX_ITERS,MAX_ORDER,MODEL(L_MODEL),
     '  NUMBER_EQN
      REAL*8 ABS_ERR,ARI(L_ARI),ARO(L_ARO),DY(NUMBER_EQN),
     '  MAX_STEP,PARAMETERS(L_PARAMETERS),REL_ERR,T,TOUT,WORK(L_WORK),
     '  Y(NUMBER_EQN)
      LOGICAL EXTEND_INTERVAL,STIFF_EQNS,USE_ROUND_CTRL
      EXTERNAL FUNC
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER MAX_RND
      INTEGER ISNOLD_IDX,K_IDX,KOLD_IDX,NUMBER_STEPS_IDX,
     '  PHASE1_IDX,REQUIRED_IDX,ROUND_CTRL_IDX,START_IDX
      INTEGER H_IDX,HOLD_IDX,TT_IDX,TOLD_IDX,ALPHA_IDX,BETA_IDX,
     '  DELSGN_IDX,DY_IDX,ERROR_WEIGHT_IDX,G_IDX,PREDICTION_IDX,
     '  PSI_IDX,SIGMA_IDX,V_IDX,W_IDX,YY_IDX,PHI_IDX
      LOGICAL ROUND_CTRL,PHASE1,START

      DATA ISNOLD_IDX,K_IDX,KOLD_IDX,NUMBER_STEPS_IDX,PHASE1_IDX,
     '  ROUND_CTRL_IDX,START_IDX / 3,1,2,4,5,6,7 /
      DATA H_IDX,HOLD_IDX,TT_IDX,TOLD_IDX,DELSGN_IDX,
     '  ALPHA_IDX / 1,2,3,4,5,6 /


      IF(L_IWORK.GE.7) THEN
        REQUIRED_IDX=7+7*MAX_ORDER+NUMBER_EQN*(MAX_ORDER+6)
        IF(USE_ROUND_CTRL) REQUIRED_IDX=REQUIRED_IDX+2*NUMBER_EQN
        IF(L_WORK.GE.REQUIRED_IDX) THEN
          BETA_IDX=ALPHA_IDX+MAX_ORDER
          DY_IDX=BETA_IDX+MAX_ORDER
          ERROR_WEIGHT_IDX=DY_IDX+NUMBER_EQN
          G_IDX=ERROR_WEIGHT_IDX+NUMBER_EQN
          PSI_IDX=G_IDX+MAX_ORDER+1
          PREDICTION_IDX=PSI_IDX+MAX_ORDER
          SIGMA_IDX=PREDICTION_IDX+NUMBER_EQN
          V_IDX=SIGMA_IDX+MAX_ORDER+1
          W_IDX=V_IDX+MAX_ORDER
          YY_IDX=W_IDX+MAX_ORDER 
          PHI_IDX=YY_IDX+NUMBER_EQN
          
          IF(USE_ROUND_CTRL) THEN
            MAX_RND=2
          ELSE
            MAX_RND=0
          ENDIF
          
          IF(ABS(IFAIL).NE.1) THEN
            START=IWORK(START_IDX).NE.0
            PHASE1=IWORK(PHASE1_IDX).NE.0
            IF(USE_ROUND_CTRL) THEN
              ROUND_CTRL=IWORK(ROUND_CTRL_IDX).NE.0
            ELSE
              ROUND_CTRL=.FALSE.
            ENDIF
          ENDIF
          
          CALL AM_DE(AII,AIO,CONTROL,ERROR_TYPE,IFAIL,
     +      IWORK(ISNOLD_IDX),IWORK(K_IDX),IWORK(KOLD_IDX),
     +      L_AII,L_AIO,L_ARI,L_ARO,L_CONTROL,L_MODEL,
     +      L_PARAMETERS,MAX_ITERS,MAX_ORDER,MAX_RND,MODEL,
     +      NUMBER_EQN,IWORK(NUMBER_STEPS_IDX),ABS_ERR,WORK(ALPHA_IDX),
     +      ARI,ARO,WORK(BETA_IDX),WORK(DELSGN_IDX),WORK(DY_IDX),DY,
     +      WORK(ERROR_WEIGHT_IDX),WORK(G_IDX),WORK(H_IDX),
     +      WORK(HOLD_IDX),MAX_STEP,PARAMETERS,WORK(PHI_IDX),
     +      WORK(PREDICTION_IDX),WORK(PSI_IDX),REL_ERR,
     +      WORK(SIGMA_IDX),T,WORK(TT_IDX),WORK(TOLD_IDX),
     +      TOUT,WORK(V_IDX),WORK(W_IDX),Y,WORK(YY_IDX),
     +      EXTEND_INTERVAL,PHASE1,ROUND_CTRL,START,STIFF_EQNS,FUNC)
          
          IF(START) THEN
            IWORK(START_IDX)=1
          ELSE
            IWORK(START_IDX)=0
          ENDIF
          IF(PHASE1) THEN
            IWORK(PHASE1_IDX)=1
          ELSE
            IWORK(PHASE1_IDX)=0
          ENDIF
          IF(USE_ROUND_CTRL) THEN
            IF(ROUND_CTRL) THEN
              IWORK(ROUND_CTRL_IDX)=1
            ELSE
              IWORK(ROUND_CTRL_IDX)=0
            ENDIF
          ENDIF
          
          IF(ABS(IFAIL).EQ.3) THEN
            ERROR='>>Integration did not reach TOUT, Tolerances too '
     '        //'small'
          ELSE IF(ABS(IFAIL).EQ.4) THEN
            ERROR='>>Integration did not reach TOUT, Maximum '
     '        //'iterations exceeded'
          ELSE IF(ABS(IFAIL).EQ.5) THEN
            ERROR='>>Integration did not reach TOUT, Equations are '
     '        //'too stiff'
          ELSE IF(ABS(IFAIL).EQ.6) THEN
            ERROR='>>Input parameters invalid'
          ENDIF
        ELSE
          WRITE(ERROR,'(''>>WORK array is too small, increase L_WORK '
     '      //'to '',I5)') REQUIRED_IDX
          IFAIL=9
        ENDIF
      ELSE
        ERROR='>>IWORK array is too small, increase L_IWORK to 7'
        IFAIL=10
      ENDIF

      RETURN
      END


      SUBROUTINE AM_DE(AII,AIO,CONTROL,ERROR_TYPE,IFAIL,ISNOLD,K,
     '  KOLD,L_AII,L_AIO,L_ARI,L_ARO,L_CONTROL,L_MODEL,L_PARAMETERS,
     '  MAX_ITERS,MAX_ORDER,MAX_RND,MODEL,NUMBER_EQN,NUMBER_STEPS,
     '  ABS_ERR,ALPHA,ARI,ARO,BETA,DELSGN,DY,DYOUT,ERROR_WEIGHT,G,
     '  H,HOLD,MAX_STEP,PARAMETERS,PHI,PREDICTION,PSI,REL_ERR,SIGMA,
     '  T,TT,TOLD,TOUT,V,W,Y,YY,EXTEND_INTERVAL,PHASE1,ROUND_CTRL,
     '  START,STIFF_EQNS,FUNC)

C#### Subroutine: AM_DE
C###  Description:
C###    Adams-Moulton differential equation integrator. Called from
C###    the ADAMS_MOULTON buffer subroutine.
      
      IMPLICIT NONE
      
!     Parameter list
      INTEGER L_AII,L_AIO,L_ARI,L_ARO,L_CONTROL,L_MODEL,L_PARAMETERS
      INTEGER AII(L_AII),AIO(L_AIO),CONTROL(L_CONTROL),ERROR_TYPE,
     '  IFAIL,ISNOLD,K,KOLD,MAX_ITERS,MAX_ORDER,MAX_RND,MODEL(L_MODEL),
     '  NUMBER_EQN,NUMBER_STEPS
      REAL*8 ABS_ERR,ALPHA(MAX_ORDER),ARI(L_ARI),ARO(L_ARO),
     '  BETA(MAX_ORDER),DELSGN,DY(NUMBER_EQN),DYOUT(NUMBER_EQN),
     '  ERROR_WEIGHT(NUMBER_EQN),G(MAX_ORDER+1),H,HOLD,MAX_STEP,
     '  PARAMETERS(L_PARAMETERS),PHI(NUMBER_EQN,MAX_ORDER+2+MAX_RND),
     '  PREDICTION(NUMBER_EQN),PSI(MAX_ORDER),REL_ERR,
     '  SIGMA(MAX_ORDER+1),T,TT,TOLD,TOUT,V(MAX_ORDER),W(MAX_ORDER),
     '  Y(NUMBER_EQN),YY(NUMBER_EQN)
      LOGICAL EXTEND_INTERVAL,PHASE1,ROUND_CTRL,START,STIFF_EQNS
      EXTERNAL FUNC
!     Local Variables
      INTEGER ERR_CODE,ISN,NUM_TIMES_LOW_ORDER,L,IT_NUM
      REAL*8 ABS_DEL,ABS_EPS,DEL,DLAMCH,EPS,EPSILON,FOUR_EPSILON,
     '  REL_EPS,TEND
      LOGICAL CRASH,FINISHED

      EPSILON=DLAMCH('EPS')
      FOUR_EPSILON=4.0d0*EPSILON

C     Test for improper parameters
      EPS=DMAX1(REL_ERR,ABS_ERR)
      IFAIL=IABS(IFAIL)
      IF(NUMBER_EQN.LT.1.OR.
     '  T.EQ.TOUT.OR.
     '  REL_ERR.LT.0.0d0.OR.ABS_ERR.LT.0.0d0.OR.EPS.LE.0.0d0.OR.
     '  (IFAIL.LT.1.OR.IFAIL.GT.5).OR.
     '  (IFAIL.NE.1.AND.T.NE.TOLD).OR.
     '  (ERROR_TYPE.LT.1.OR.ERROR_TYPE.GT.4).OR.
     '  (MAX_ORDER.LT.1.OR.MAX_ORDER.GT.12)) THEN
        IFAIL = 6
      ELSE

C       On each call set interval of integration and counter for the
C       number of steps.
        DEL=TOUT-T
        ABS_DEL=DABS(DEL)
        TEND=T+10.0d0*DEL
        ISN=ISIGN(1,IFAIL)
        IF(ISN.LT.0) TEND = TOUT
        IT_NUM=0
        NUM_TIMES_LOW_ORDER=0
        STIFF_EQNS=.FALSE.
        REL_EPS=REL_ERR/EPS
        ABS_EPS=ABS_ERR/EPS
        IF(IFAIL.EQ.1.OR.ISNOLD.LT.0.OR.DELSGN*DEL.LE.0.0d0) THEN
C         On start and restart also set the work variables TT and YY,
C         store the direction of integration and initialise the step
C         size
          START=.TRUE.
          TT=T
          DO l=1,NUMBER_EQN
            YY(l)=Y(l)
          ENDDO !l
          DELSGN=DSIGN(1.0d0,DEL)
          H=DSIGN(DMAX1(DABS(TOUT-TT),FOUR_EPSILON*DABS(TT)),TOUT-TT)
          IF(H.GT.MAX_STEP) H=MAX_STEP
        ENDIF

C       Iterate on the solution until finished.
        FINISHED=.FALSE.
        DO WHILE(.NOT.FINISHED)

          IF(DABS(TT-T).GE.ABS_DEL) THEN
C           If already past the output point then interpolate and 
C           finish
            CALL AM_INTERPOLATE(KOLD,MAX_ORDER,MAX_RND,
     '        NUMBER_EQN,DYOUT,PHI,PSI,TT,TOUT,YY,Y)
            IFAIL=2
            T=TOUT
            TOLD=T
            ISNOLD=ISN
            FINISHED=.TRUE.
          ELSE IF(.NOT.EXTEND_INTERVAL.AND.
     '        DABS(TOUT-TT).LT.FOUR_EPSILON*DABS(TT)) THEN
C           If you cannot extend the interval to go past the output
C           output point and you are sufficiently close then
C           extroplate and finish.
            H=TOUT-TT
            CALL FUNC(TT,YY,DY,AII,AIO,CONTROL,L_AII,L_AIO,
     '        L_ARI,L_ARO,L_CONTROL,L_MODEL,L_PARAMETERS,MODEL,
     '        NUMBER_EQN,ARI,ARO,PARAMETERS,ERR_CODE)
            IF(ERR_CODE.EQ.0) THEN
              DO l=1,NUMBER_EQN
                Y(l)=YY(l)+H*DY(l)
              ENDDO !l
              IFAIL=2
              T=TOUT
              TOLD=T
              ISNOLD=ISN
              FINISHED=.TRUE.
            ELSE
              IFAIL=7
              FINISHED=.TRUE.
            ENDIF
          ELSE IF(IT_NUM.GE.MAX_ITERS) THEN
C           Test for too many steps            
            IFAIL=ISN*4
            IF(STIFF_EQNS) IFAIL=ISN*5
            DO l=1,NUMBER_EQN
              Y(l)=YY(l)
            ENDDO !l
            T=TT
            TOLD=T
            ISNOLD=1
            FINISHED=.TRUE.
          ELSE
C           Limit step size, set the error weight vector and take 
C           a step
            H=DSIGN(DMIN1(DABS(H),DABS(TEND-TT)),H)
            IF(H.GT.MAX_STEP) H=MAX_STEP
            ERR_CODE=0
            IF(ERROR_TYPE.EQ.1) THEN
              DO l=1,NUMBER_EQN
                ERROR_WEIGHT(l)=1.0d0
              ENDDO !l
            ELSE IF(ERROR_TYPE.EQ.2) THEN              
              DO l=1,NUMBER_EQN
                ERROR_WEIGHT(l)=DABS(YY(l))
                IF(DABS(YY(l)).LE.EPSILON) ERR_CODE=1
              ENDDO !l
            ELSE IF(ERROR_TYPE.EQ.3) THEN
              DO l=1,NUMBER_EQN
                ERROR_WEIGHT(l)=DABS(DY(l))
                IF(DABS(DY(l)).LE.EPSILON) ERR_CODE=1
              ENDDO !l
            ELSE IF(ERROR_TYPE.EQ.4) THEN
              DO l=1,NUMBER_EQN
                ERROR_WEIGHT(l)=REL_EPS*DABS(YY(l))+ABS_EPS
              ENDDO !l
            ENDIF
            IF(ERR_CODE.EQ.0) THEN
            
              CALL AM_STEP(AII,AIO,CONTROL,ERR_CODE,K,KOLD,L_AII,
     '          L_AIO,L_ARI,L_ARO,L_CONTROL,L_MODEL,L_PARAMETERS,
     '          MAX_ORDER,MAX_RND,MODEL,NUMBER_EQN,NUMBER_STEPS,
     '          ALPHA,ARI,ARO,BETA,DY,EPS,EPSILON,ERROR_WEIGHT,G,H,
     '          HOLD,MAX_STEP,PARAMETERS,PHI,PREDICTION,PSI,SIGMA,
     '          TT,V,W,YY,CRASH,PHASE1,ROUND_CTRL,START,FUNC)

              IF(ERR_CODE.EQ.0) THEN
                IF(CRASH) THEN
C                 Tolerances too small
                  IFAIL=ISN*3
                  REL_ERR=EPS*REL_EPS
                  ABS_ERR=EPS*ABS_EPS
                  DO l=1,NUMBER_EQN
                    Y(l)=YY(l)
                  ENDDO !l
                  T=TT
                  TOLD=T
                  ISNOLD=1
                  FINISHED=.TRUE.
                ELSE 
C                 Adjust number of steps and test for stiffness
                  IT_NUM=IT_NUM+1
                  NUM_TIMES_LOW_ORDER=NUM_TIMES_LOW_ORDER+1
                  IF(KOLD.GT.4) NUM_TIMES_LOW_ORDER=0
                  IF(NUM_TIMES_LOW_ORDER.GE.50) STIFF_EQNS=.TRUE.
                ENDIF
              ELSE
C               Error from func
                IFAIL=7
                FINISHED=.TRUE.
              ENDIF
            ELSE
C             Error weight is zero
              IFAIL=8
              FINISHED=.TRUE.
            ENDIF
          ENDIF          
        ENDDO
      ENDIF

      RETURN
      END


      SUBROUTINE AM_STEP(AII,AIO,CONTROL,ERR_CODE,K,KOLD,L_AII,
     '  L_AIO,L_ARI,L_ARO,L_CONTROL,L_MODEL,L_PARAMETERS,MAX_ORDER,
     '  MAX_RND,MODEL,NUMBER_EQN,NUMBER_STEPS,ALPHA,ARI,ARO,BETA,
     '  DY,EPS,EPSILON,ERROR_WEIGHT,G,H,HOLD,MAX_STEP,PARAMETERS,PHI,
     '  PREDICTION,PSI,SIGMA,T,V,W,Y,CRASH,PHASE1,ROUND_CTRL,START,
     '  FUNC)

C#### Subroutine: AM_STEP
C###  Description:
C###    Subroutine to perfom the Adams-Moulton step from T to T+H
      
      IMPLICIT NONE
!     Parameter list
      INTEGER L_AII,L_AIO,L_ARI,L_ARO,L_CONTROL,L_MODEL,L_PARAMETERS
      INTEGER AII(L_AII),AIO(L_AIO),CONTROL(L_CONTROL),ERR_CODE,
     '  K,KOLD,MAX_ORDER,MAX_RND,MODEL(L_MODEL),NUMBER_EQN,NUMBER_STEPS
      REAL*8 ALPHA(MAX_ORDER),ARI(L_ARI),ARO(L_ARO),BETA(MAX_ORDER),
     '  DY(NUMBER_EQN),EPS,EPSILON,ERROR_WEIGHT(NUMBER_EQN),
     '  G(MAX_ORDER+1),H,HOLD,PREDICTION(NUMBER_EQN),MAX_STEP,
     '  PARAMETERS(L_PARAMETERS),PHI(NUMBER_EQN,MAX_ORDER+2+MAX_RND),
     '  PSI(MAX_ORDER),SIGMA(MAX_ORDER+1),T,V(MAX_ORDER),W(MAX_ORDER),
     '  Y(NUMBER_EQN)
      LOGICAL CRASH,PHASE1,ROUND_CTRL,START
      EXTERNAL FUNC
!     Local Variables
      INTEGER i,IFAIL,iq,j,K_MINUS1,K_MINUS2,KNEW,K_PLUS1,K_PLUS2,l,
     '  NUMBER_STEPS_MINUS2,NUMBER_STEPS_PLUS1,NUMBER_STEPS_PLUS2
      REAL*8  ABSH,ERROR_K,ERROR_K_MINUS1,ERROR_K_MINUS2,
     '  ERROR_K_PLUS1,ERR,FOUR_EPSILON,GAMMASTAR(13),HALF_EPS,
     '  HNEW,R,RHO,ROUND,SUM,TAU,TEMP1,TEMP2,TOLD,POWERTWO(13),
     '  TWO_EPSILON
      LOGICAL GOODSTEP

C     Gamma*_(i) is defined to be:
C       Gamma*_(i)=1/i! \int{0}{1}[(s-1)(s)(s+1)...(s+i-2)]ds i=1,2,3,..
C     with Gamma*_(0)=1, but can be generated with the recursion 
C     formula:
C       Gamma*_(m)+1/2.Gamma*_(m-1)+1/3.Gamma*_(m-2)+....+
C       1/(m+1).Gamma*_(0)=0
C     NOTE: For all the cases below Gamma* is negative however the
C     absolute value is stored as we are only interested in the 
C     absolute error.
      REAL*8 GAMMASTAR_1,GAMMASTAR_2,GAMMASTAR_3,GAMMASTAR_4,
     '  GAMMASTAR_5,GAMMASTAR_6,GAMMASTAR_7,GAMMASTAR_8,GAMMASTAR_9,
     '  GAMMASTAR_10,GAMMASTAR_11,GAMMASTAR_12,GAMMASTAR_13
      PARAMETER(GAMMASTAR_1=1.0d0/2.0d0,GAMMASTAR_2=1.0d0/12.0d0,
     '  GAMMASTAR_3=1.0d0/24.0d0,GAMMASTAR_4=19.0d0/720.0d0,
     '  GAMMASTAR_5=27.0d0/1440.0d0,GAMMASTAR_6=863.0d0/60480.0d0,
     '  GAMMASTAR_7=275.0d0/24192.0d0,GAMMASTAR_8=33953.0d0/3628800.0d0,
     '  GAMMASTAR_9=8183.0d0/1036800.0d0,
     '  GAMMASTAR_10=3250433.0d0/479001600.0d0,
     '  GAMMASTAR_11=4671.0d0/788480.0d0,
     '  GAMMASTAR_12=5852897.0d0/1117670400.0d0,
     '  GAMMASTAR_13=78418523.0d0/16765056000.0d0)
      DATA GAMMASTAR / GAMMASTAR_1,GAMMASTAR_2,GAMMASTAR_3,GAMMASTAR_4,
     '  GAMMASTAR_5,GAMMASTAR_6,GAMMASTAR_7,GAMMASTAR_8,GAMMASTAR_9,
     '  GAMMASTAR_10,GAMMASTAR_11,GAMMASTAR_12,GAMMASTAR_13 /

      DATA POWERTWO /2.0d0,4.0d0,8.0d0,16.0d0,32.0d0,64.0d0,128.0d0,
     '  256.0d0,512.0d0,1024.0d0,2048.0d0,4096.0d0,8192.0d0/


      TWO_EPSILON=2.0d0*EPSILON
      FOUR_EPSILON=4.0d0*EPSILON

C     *** BLOCK 0 ***
C     Check if the step size or error tolerance is too small for the
C     machine precision. If it is the first step then initialise the
C     PHI array and estimate a starting step size.

      CRASH=.TRUE.
      IF(DABS(H).LT.FOUR_EPSILON*DABS(T)) THEN
C       If step size is too small the determine an acceptable one and
C       exit
        H=DSIGN(FOUR_EPSILON*DABS(T),H)
        IF(H.GT.MAX_STEP) H=MAX_STEP
      ELSE
        HALF_EPS=EPS/2.0d0
C
        ROUND=0.0d0
        DO l=1,NUMBER_EQN
          ROUND=ROUND+(Y(l)/ERROR_WEIGHT(l))**2
        ENDDO !l
        ROUND=TWO_EPSILON*DSQRT(ROUND)
        IF(HALF_EPS.LT.ROUND) THEN
C         If error tolerance is too small then increase it to an
C         acceptable value and exit
          EPS=2.0d0*ROUND*(1.0d0+FOUR_EPSILON)
        ELSE
          CRASH=.FALSE.
          G(1)=1.0d0
          G(2)=0.5d0
          SIGMA(1)=1.0d0
          
          IF(START) THEN
C           Initialise and compute appropriate step size for the first
C           step
            CALL FUNC(T,Y,DY,AII,AIO,CONTROL,L_AII,L_AIO,L_ARI,
     '        L_ARO,L_CONTROL,L_MODEL,L_PARAMETERS,MODEL,NUMBER_EQN,
     '        ARI,ARO,PARAMETERS,ERR_CODE)
            IF(ERR_CODE.EQ.0) THEN
              SUM=0.0d0
              DO l=1,NUMBER_EQN
                PHI(l,1)=DY(l)
                PHI(l,2)=0.0d0
                SUM=SUM+(DY(l)/ERROR_WEIGHT(l))**2
              ENDDO !l
              SUM=DSQRT(SUM)
              IF(EPS.LT.16.0d0*SUM*H*H) THEN
                ABSH=0.25d0*DSQRT(EPS/SUM)
              ELSE
                ABSH=DABS(H)
              ENDIF
              H=DSIGN(DMAX1(ABSH,FOUR_EPSILON*DABS(T)),H)
              IF(H.GT.MAX_STEP) H=MAX_STEP
              HOLD=0.0d0
              K=1
              KOLD=0
              START=.FALSE.
              PHASE1=.TRUE.
              IF(MAX_RND.GT.0) THEN
                ROUND_CTRL=.FALSE.
                IF(HALF_EPS.LE.100.0d0*ROUND) THEN
                  ROUND_CTRL=.TRUE.
                  DO l=1,NUMBER_EQN
                    PHI(l,MAX_ORDER+3)=0.0d0
                  ENDDO !l
                ENDIF
              ELSE
                ROUND_CTRL=.FALSE.
              ENDIF
            ENDIF
          ENDIF
C         *** End BLOCK 0 ***

          IF(ERR_CODE.EQ.0) THEN
            IFAIL=0
          
C           *** BLOCK 1 ***
C           Compute the coefficients of the formulas for this step.
C           Avoid computing those quantities not changed when the 
C           step size is not changed.   

            GOODSTEP=.FALSE.
            DO WHILE(.NOT.GOODSTEP.AND..NOT.CRASH.AND.ERR_CODE.EQ.0)
            
              K_PLUS1=K+1
              K_PLUS2=K+2
              K_MINUS1=K-1
              K_MINUS2=K-2

C             NUMBER_STEPS is the number of steps taken with step size
C             H. When the current order K is less than NUMBER_STEPS
C             then no coefficients change.

              IF(H.NE.HOLD) NUMBER_STEPS = 0
              IF(NUMBER_STEPS.LE.KOLD) NUMBER_STEPS=NUMBER_STEPS+1
              NUMBER_STEPS_PLUS1 = NUMBER_STEPS+1
              IF (K.GE.NUMBER_STEPS) THEN
C               Compute those components of ALPHA, BETA, PSI and SIGMA
C               which are changed
C               The formula are:
C                 Psi_i(n+1)=h_(n+1)+h_n+...+h_(n+2-i)      i>=1
C                 Alpha_i(n+1)=h_(n+1)/Psi_i(n+1)           i>=1
C                 Beta_1(n+1)=1.0                           i=1
C                 Beta_i(n+1)=Psi_1(n+1).Psi_2(n+1)...Psi_(i-1)(n+1)/
C                   Psi_1(n).Psi_2(n)...Psi_(i-1)(n)        i>1
C                 Sigma_1(n+1)=1.0                          i=1
C                 Sigma_i(n+1)=h.2h...(i-1)h/
C                   Psi_1(n+1).Psi_2(n+1)...Psi_(i-1)(n+1)  i=>1
                BETA(NUMBER_STEPS)=1.0d0
                ALPHA(NUMBER_STEPS)=1.0d0/DBLE(NUMBER_STEPS)
                TEMP1=H*DBLE(NUMBER_STEPS)
                SIGMA(NUMBER_STEPS_PLUS1)=1.0d0
                DO i=NUMBER_STEPS_PLUS1,K
                  TEMP2=PSI(i-1)
                  PSI(i-1)=TEMP1
                  BETA(i)=BETA(i-1)*PSI(i-1)/TEMP2
                  TEMP1=TEMP2+H
                  ALPHA(i)=H/TEMP1
                  SIGMA(i+1)=DBLE(i)*ALPHA(i)*SIGMA(i)
                ENDDO !i
                PSI(K)=TEMP1

C               Compute the coefficients G
                IF(NUMBER_STEPS.LE.1) THEN
C                 Initialise V and set W
                  DO iq=1,K
                    V(iq)=1.0d0/DBLE(iq*(iq+1))
                    W(iq)=V(iq)
                  ENDDO !iq
                ELSE 
                  IF(K.GT.KOLD) THEN
C                   If the order was raised update the diagonal part
C                   of V
                    V(K)=1.0d0/DBLE(K*K_PLUS1)
                    NUMBER_STEPS_MINUS2=NUMBER_STEPS-2
                    DO j=1,NUMBER_STEPS_MINUS2
                      i=K-j
                      V(i)=V(i)-ALPHA(j+1)*V(i+1)
                    ENDDO !J
                  ENDIF
C                 Update V and set W
                  DO iq=1,K_PLUS1-NUMBER_STEPS
                    V(iq)=V(iq)-ALPHA(NUMBER_STEPS)*V(iq+1)
                    W(iq)=V(iq)
                  ENDDO !IQ
                  G(NUMBER_STEPS_PLUS1)=W(1)
C
                ENDIF
C               Compute the G in the work vector W
                NUMBER_STEPS_PLUS2 = NUMBER_STEPS + 2
                DO i=NUMBER_STEPS_PLUS2,K_PLUS1
                  DO iq=1,K_PLUS2-i
                    W(iq)=W(iq)-ALPHA(i-1)*W(iq+1)
                  ENDDO !iq
                  G(i)=W(1)
                ENDDO !i
                
              ENDIF
C             *** End BLOCK 1 ***
     
C             *** BLOCK 2 ***
C             Predicit a solution, PREDICTION, and evaluate the 
C             derivatives, DY, using the predicited solution. 
C             Estimate the local error at order K and the errors
C             at orders K, K-1 and K-2 as if a constant step size
C             were used.

C             Change PHI to PHI* i.e. PHI*_i(n)=BETA_i(n+1).PHI_i(n)
              DO i=NUMBER_STEPS_PLUS1,K
                DO l=1,NUMBER_EQN
                  PHI(l,i)=BETA(i)*PHI(l,i)
                ENDDO !l
              ENDDO !i
              
C             Predicit the solution and differences
              DO l=1,NUMBER_EQN
                PHI(l,K_PLUS2)=PHI(l,K_PLUS1)
                PHI(l,K_PLUS1)=0.0d0
                PREDICTION(l)=0.0d0
              ENDDO !l
              DO j=1,K
                i=K_PLUS1-j
                DO l=1,NUMBER_EQN
                  PREDICTION(l)=PREDICTION(l)+G(i)*PHI(l,i)
                  PHI(l,i)=PHI(l,i)+PHI(l,i+1)
                ENDDO !l
              ENDDO !i
              IF(ROUND_CTRL) THEN
                DO l=1,NUMBER_EQN
                  TAU=H*PREDICTION(l)-PHI(l,MAX_ORDER+3)
                  PREDICTION(l)=Y(l)+TAU
                  PHI(l,MAX_ORDER+4)=(PREDICTION(l)-Y(l))-TAU
                ENDDO !l
              ELSE
                DO l=1,NUMBER_EQN
                  PREDICTION(l)=Y(l)+H*PREDICTION(l)
                ENDDO !l
              ENDIF
              TOLD=T
              T=T+H
              ABSH=DABS(H)
              CALL FUNC(T,PREDICTION,DY,AII,AIO,CONTROL,L_AII,
     '          L_AIO,L_ARI,L_ARO,L_CONTROL,L_MODEL,L_PARAMETERS,
     '          MODEL,NUMBER_EQN,ARI,ARO,PARAMETERS,ERR_CODE)
              IF(ERR_CODE.EQ.0) THEN
C               Estimate errors at orders K,K-1,K-2
C               For a constant step size, h, the local error at x_(n+1)
C               is given by:
C                 ERROR_K=|h.Gamma*_K.Sigma_(K+1)(n+1).Phi_(K+1)^P(n+1)|
                ERROR_K=0.0d0
                DO l=1,NUMBER_EQN
                  ERROR_K=ERROR_K+
     '              ((DY(l)-PHI(l,1))/ERROR_WEIGHT(l))**2
                ENDDO !l
                ERR=ABSH*DSQRT(ERROR_K)*(G(K)-G(K_PLUS1))
                ERROR_K=ABSH*DSQRT(ERROR_K)*SIGMA(K_PLUS1)*
     '            GAMMASTAR(K)
                IF(K_MINUS2.EQ.0) THEN
                  ERROR_K_MINUS1=0.0d0
                  DO l=1,NUMBER_EQN
                    ERROR_K_MINUS1=ERROR_K_MINUS1+
     '                ((DY(l)+PHI(L,K)-PHI(l,1))/ERROR_WEIGHT(l))**2
                  ENDDO !l
                  ERROR_K_MINUS1=ABSH*SIGMA(K)*GAMMASTAR(K_MINUS1)*
     '              DSQRT(ERROR_K_MINUS1)
                ELSE IF(K_MINUS2.GT.0) THEN
                  ERROR_K_MINUS2=0.0d0
                  ERROR_K_MINUS1=0.0d0
                  DO l=1,NUMBER_EQN
                    ERROR_K_MINUS2=ERROR_K_MINUS2+
     '                ((DY(l)+PHI(L,K_MINUS1)-PHI(l,1))/
     '                ERROR_WEIGHT(l))**2
                    ERROR_K_MINUS1=ERROR_K_MINUS1+ 
     '                ((DY(l)+PHI(L,K)-PHI(l,1))/ERROR_WEIGHT(l))**2
                  ENDDO !l
                  ERROR_K_MINUS2=ABSH*SIGMA(K_MINUS1)*
     '              GAMMASTAR(K_MINUS2)*DSQRT(ERROR_K_MINUS2)
                  ERROR_K_MINUS1=ABSH*SIGMA(K)*GAMMASTAR(K_MINUS1)*
     '              DSQRT(ERROR_K_MINUS1)
                ENDIF
                KNEW=K

C               Test if order should be lowered                
                IF(K_MINUS2.EQ.0) THEN
                  IF(ERROR_K_MINUS1.LE.0.5D0*ERROR_K) KNEW = K_MINUS1
                ELSE IF(K_MINUS2.GT.0) THEN
                  IF(DMAX1(ERROR_K_MINUS1,ERROR_K_MINUS2).LE.ERROR_K) 
     '              KNEW = K_MINUS1
                ENDIF
C               *** End BLOCK 2 ***
        
C               Test if the step was successful
                IF(ERR.GT.EPS) THEN

C                 *** BLOCK 3 ***
C                 The step is unsuccessful. 

C                 Restore T, PHI and PSI
                  PHASE1=.FALSE.
                  T=TOLD
                  DO i=1,K
                    DO l=1,NUMBER_EQN
                      PHI(l,i)=(PHI(l,i)-PHI(l,i+1))/BETA(i)
                    ENDDO !l
                  ENDDO !i
                  IF(K.GE.2) THEN
                    DO i=2,K
                      PSI(i-1)=PSI(i)-H
                    ENDDO !i
                  ENDIF

                  IFAIL=IFAIL+1
C                 Double the step size. If the step fails three times
C                 set the order to 1, thereafter use optimal step size.
C                 This procedure will exit if the estimated step size
C                 is too small for the machine precision.

                  IF(IFAIL.LT.3) THEN
                    H=H/2.0d0
                  ELSE IF(IFAIL.EQ.3) THEN
                    KNEW = 1
                    H=H/2.0d0
                  ELSE 
                    KNEW = 1
                    IF(HALF_EPS.LT.ERROR_K/4.0d0) THEN
                      H=DSQRT(HALF_EPS/ERROR_K)*H
                    ELSE
                      H=H/2.0d0
                    ENDIF
                  ENDIF
                  K = KNEW
                  IF(DABS(H).LT.FOUR_EPSILON*DABS(T)) THEN
                    CRASH =.TRUE.
                    H=DSIGN(FOUR_EPSILON*DABS(T),H)
                    IF(H.GT.MAX_STEP) H=MAX_STEP
                    EPS=EPS+EPS                
                  ENDIF
C                 *** End BLOCK 3 ***
              
                ELSE
                  GOODSTEP=.TRUE.
                ENDIF
              ENDIF
              
            ENDDO
            
            IF(.NOT.CRASH.AND.ERR_CODE.EQ.0) THEN
      
C             *** BLOCK 4 ***
C             The step is successfull. Correct the predicited solution,
C             evaluate the derivatives using the correct solution and
C             update the differences. Determine best order and step
C             size for the next step.

              KOLD=K
              HOLD=H
C             Correct and Evaluate
              IF(ROUND_CTRL) THEN
                DO l=1,NUMBER_EQN
                  RHO=H*G(K_PLUS1)*(DY(l)-PHI(l,1))-PHI(l,MAX_ORDER+4)
                  Y(l)=PREDICTION(l)+RHO
                  PHI(l,MAX_ORDER+3)=(Y(l)-PREDICTION(l))-RHO
                ENDDO !l
              ELSE
                DO l=1,NUMBER_EQN
                  Y(l)=PREDICTION(l)+H*G(K_PLUS1)*(DY(l)-PHI(l,1))
                ENDDO !l
              ENDIF
              CALL FUNC(T,Y,DY,AII,AIO,CONTROL,L_AII,L_AIO,L_ARI,
     '          L_ARO,L_CONTROL,L_MODEL,L_PARAMETERS,MODEL,
     '          NUMBER_EQN,ARI,ARO,PARAMETERS,ERR_CODE)
              IF(ERR_CODE.EQ.0) THEN
C               Update differences for next step
                DO l=1,NUMBER_EQN
                  PHI(l,K_PLUS1)=DY(l)-PHI(l,1)
                  PHI(l,K_PLUS2)=PHI(l,K_PLUS1)-PHI(l,K_PLUS2)
                ENDDO !l
                DO i=1,K
                  DO l=1,NUMBER_EQN
                    PHI(l,i)=PHI(l,i)+PHI(l,K_PLUS1)
                  ENDDO !l
                ENDDO !k

C               Estimate the error at order K+1 unless: we are in
C               the first phase in which case always raise the order;
C               we have already decieded to lower the order; or the
C               step size is not constant so the error estimate is
C               unreliable.
                ERROR_K_PLUS1=0.0d0
                IF(KNEW.EQ.K_MINUS1.OR.K.EQ.MAX_ORDER) PHASE1=.FALSE.
                IF(PHASE1) THEN 
C                 Raise the order
                  K=K_PLUS1
                  ERROR_K=ERROR_K_PLUS1
                ELSE              
                  IF(KNEW.EQ.K_MINUS1) THEN
C                   Lower the order
                    K=K_MINUS1
                    ERROR_K=ERROR_K_MINUS1                   
                  ELSE IF(K_PLUS1.LE.NUMBER_STEPS) THEN
                    DO l=1,NUMBER_EQN
                      ERROR_K_PLUS1=ERROR_K_PLUS1+
     '                  (PHI(l,K_PLUS2)/ERROR_WEIGHT(l))**2
                    ENDDO !l
                    ERROR_K_PLUS1=ABSH*GAMMASTAR(K_PLUS1)*
     '                DSQRT(ERROR_K_PLUS1)
C                   Using the estimated error at order K+1 determine
C                   the appropriate order for the next step.
                    IF(K.LE.1) THEN
                      IF(ERROR_K_PLUS1.LT.ERROR_K/2.0d0.AND.
     '                  K.NE.MAX_ORDER) THEN
C                       Raise the order
                        K=K_PLUS1
                        ERROR_K=ERROR_K_PLUS1
                      ENDIF
                    ELSE
                      IF(ERROR_K_MINUS1.LE.
     '                  DMIN1(ERROR_K,ERROR_K_PLUS1)) THEN
C                       Lower the order
                        K=K_MINUS1
                        ERROR_K=ERROR_K_MINUS1
                      ELSE IF(ERROR_K_PLUS1.LT.ERROR_K.AND.
     '                    K.NE.MAX_ORDER) THEN
C                       Here the error at K+1 < the error at K < the
C                       maximum error of K-1 and K-2 otherwise the
C                       order would have been lowered in block 2, thus
C                       raise the order
                        K = K_PLUS1
                        ERROR_K = ERROR_K_PLUS1
                      ENDIF
                    ENDIF                
                  ENDIF
                ENDIF

C               With the new order determine the appropriate step
C               size for the next step.
                HNEW=H+H
                IF(.NOT.PHASE1.AND.
     '            HALF_EPS.LT.ERROR_K*POWERTWO(K+1)) THEN
                  HNEW=H
                  IF(HALF_EPS.LT.ERROR_K) THEN
                    R=(HALF_EPS/ERROR_K)**(1.0d0/DBLE(K+1))
                    HNEW=ABSH*DMAX1(0.5d0,DMIN1(0.9d0,R))
                    HNEW=DSIGN(DMAX1(HNEW,FOUR_EPSILON*DABS(T)),H)
                  ENDIF
                ENDIF
                H=HNEW
                IF(H.GT.MAX_STEP) H=MAX_STEP
              ENDIF
C             *** End BLOCK 4 ***
            
            ENDIF
          ENDIF
        ENDIF
      ENDIF

      RETURN
      END


      SUBROUTINE AM_INTERPOLATE(K,MAX_ORDER,MAX_RND,
     '  NUMBER_EQN,DYOUT,PHI,PSI,T,TOUT,Y,YOUT)

C#### Subroutine: AM_INTERPOLATE
C###  Description:
C###    Interpolates the solution at TOUT using the Kth order Adams
C###    polynomial calculated near T in AM_STEP. On output YOUT
C###    contains the solution at TOUT and DYOUT contains the 
C###    derivative of the solution at TOUT.

      IMPLICIT NONE
!     Parameter list
      INTEGER K,MAX_ORDER,MAX_RND,NUMBER_EQN
      REAL*8 PHI(NUMBER_EQN,MAX_ORDER+2+MAX_RND),PSI(MAX_ORDER),T,
     '  TOUT,Y(NUMBER_EQN),YOUT(NUMBER_EQN),DYOUT(NUMBER_EQN)
!     Local Variables
      INTEGER i,j,l
      REAL*8 ETA,G(13),GAMMA,H_I,RHO(13),TERM,W(13)

      DATA G(1) /1.0d0/
      DATA RHO(1) /1.0d0/

C     The formula for interpolation is
C       y_out=y_(n+1)+h_I\sum{i=1}{K+1}[g_i,1^I.Phi_i(n+1)]
C       dy_out=\sum{i=1}{K+1}[rho_i^I.Phi_i(n+1)]
C     where
C       h_I=t_out-t_(n+1)
C       g_i,q^I=1/q                                       i=1
C       g_i,q^I=GAMMA_(I-1)(1).g_(i-1),q^I-
C         eta_(i-1).g_(i-1),(q+1)^I                       i>=2
C       eta_i=h_I/Psi_i(n+1)
C       GAMMA_i(s)=s.h_I/Psi_1(n+1)                       i=1
C       GAMMA_i(s)=(s.h_I+Psi_(i-1)(n+1))/Psi_i(n+1)      i>=2
C       rho_1=1.0                                         i=1
C       rho_i=rho_(i-1).GAMMA_(i-1)(1)                    i=2,..,K+1

      H_I=TOUT-T

C     Initialise W for computing G
      DO i=1,K+1
        W(i)=1.0d0/DBLE(i)
      ENDDO !i
      TERM=0.0d0

C     Compute G
      DO j=2,K+1
        GAMMA=(H_I+TERM)/PSI(j-1)
        ETA=H_I/PSI(j-1)
        DO i=1,K+2-j
          W(i)=GAMMA*W(i)-ETA*W(i+1)
        ENDDO !i
        G(j)=W(1)
        RHO(j)=GAMMA*RHO(j-1)
        TERM=PSI(j-1)
      ENDDO !j

C     Interpolate
      DO l=1,NUMBER_EQN
        YOUT(l)=0.0d0
        DYOUT(l)=0.0d0
      ENDDO !l
      DO j=1,K+1
        i=K+2-j
        DO l=1,NUMBER_EQN
          YOUT(l)=YOUT(l)+G(i)*PHI(l,i)
          DYOUT(l)=DYOUT(l)+RHO(i)*PHI(l,i)
        ENDDO !l
      ENDDO !j
      DO l=1,NUMBER_EQN
        YOUT(l)=Y(l)+H_I*YOUT(l)
      ENDDO !l
      
      RETURN
      END


      SUBROUTINE ADAMS(SUB_NAME,NT_EQN,TBEG,HMAX,ALPH,H,HOLD,
     '  T,THIRD,H202,H303,Y,F,D1,DDIF,JFLAG,NFE,NSTEP,NBUMP,IFLAG)

C#### Subroutine: ADAMS
C###  Description:
C###    This procedure integrates the solution y from t to t+h using a 
C###    3rd order variable stepsize Adams p2ec3e formula. 
C###    On first call Adams also sets coefficients and does a pe1ce2 
C###    step. After each step the estimated local error is tested. 
C###    If the error is too large, the step is rejected and repeated 
C###    with smaller h. The variables are defined as they are used.

      IMPLICIT NONE
c      INCLUDE 'cmiss$reference:oxs005.cmn'
!     Parameter List
      INTEGER IFLAG,JFLAG,NBUMP,NFE,NSTEP,NT_EQN
      REAL*8 ALPH,BET,D1(*),DDIF(*),F(*),H,H202,H303,HMAX,HOLD,T,TBEG,
     '  THIRD,Y(*)
      CHARACTER SUB_NAME*(*)
!     Local Variables
      INTEGER J,M,NREDO!,I
      REAL*8 D2(99),EPS,EPSP1,ERR,ERRD,ERREPS,ERRMAX,ERTEST,FP(99),
     '  H1,HDEMX,HMIN,P(99),R,RMAX,TTST,Z,membrane_current
      LOGICAL DISCONT

      VOLATILE ERRMAX,HMIN,ALPH,BET

!      write(*,'('' >>>call adams'')')
      DISCONT=.FALSE.
 1    IF(IFLAG.EQ.0) THEN
        IF(DISCONT) GOTO 11

C ***   Compute machine epsilon.
        EPS=1.0d0
        EPSP1=2.0d0
        DO WHILE(EPSP1.GT.1.0d0)
          EPS=0.5d0*EPS
          EPSP1=EPS+1.0d0
        ENDDO
        EPS=2.0d0*EPS
        
!!!        WRITE(*,*) 'eps',EPS,EPSP1

C ***   hmin is smallest allowed step size and starting stepsize.
        HMIN=100.0d0*EPS

C ***   errmax is largest acceptable relative local error.
        ERRMAX=0.005d0
!!!        WRITE(*,*) 'ERRMAX = ',ERRMAX
C ***   erreps is inserted into denom of error to prevent / by zero.
        ERREPS=0.001d0*ERRMAX
!!!        WRITE(*,*) 'ERREPS = ',ERREPS

C ***   After each step the value of next stepsize is estimated. This
C ***   value is multiplied by alph to be on safe side. If a step is
C ***   rejected, h is multiplied by bet before trying again.
        ALPH=0.5d0
        BET=0.7d0
        THIRD=1.0d0/3.0d0
        NBUMP=0
        T=TBEG
!!!        WRITE(*,*) ALPH,BET,THIRD,NBUMP,T

C Resetting initial H to avoid underflow with optimisation
C 11     H=HMIN
 11     H=0.0001d0
        JFLAG=1

C ***   First step uses euler predictor and trapezoidal corrector.
        H202=0.0d0
        H303=0.5d0*H*H
        HOLD=1.0d0-H
!!!        WRITE(*,*) 'h202,h303,hold',h202,h303,hold

        IF(SUB_NAME(1:2).EQ.'HH') THEN       !Hodgkin-Huxley eqns
          CALL HH(T,Y,F)
        ELSE IF(SUB_NAME(1:2).EQ.'BR') THEN  !Beeler-Reuter eqns
!         CALL BR(T,Y,F)
        ELSE IF(SUB_NAME(1:3).EQ.'RBR') THEN !Reduced Beeler-Reuter eqns
!         CALL RBR(T,Y,F)
        ELSE IF(SUB_NAME(1:4).EQ.'DN') THEN  !diFrancesco-Noble eqns
          CALL DN(T,Y,F)
        ELSE IF(SUB_NAME(1:2).EQ.'LR') THEN  !Luo-Rudy eqns
          CALL LR(T,Y,F)
!         IF(DOP) THEN
!            WRITE(99,'('' T='',E11.3,'' IFLAG='',I3)') T,IFLAG
!            WRITE(99,'('' Y( 1..10):'',10E11.3)') (Y(I),I= 1,10)
!            WRITE(99,'('' Y(11..20):'',10E11.3)') (Y(I),I=11,20)
!            WRITE(99,'('' Y(21..30):'',10E11.3)') (Y(I),I=21,30)
c            WRITE(99,'('' T='',E11.3,'' IFLAG='',I3)') T,IFLAG
c            WRITE(99,'('' Y( 1.. 7):'',10E11.3)') (Y(I),I= 1, 7)
c            WRITE(99,'('' Y( 8..14):'',10E11.3)') (Y(I),I= 8,14)
!         ENDIF
        ELSE IF(SUB_NAME(1:3).EQ.'JRW') THEN  !Jafri-Rice-Winslow eqns
          CALL JRW(T,Y,F)
        ELSE IF(SUB_NAME(1:2).EQ.'DM') THEN  !Distribution moment model
          CALL FCN_DM(T,Y,F)
        ELSE IF(SUB_NAME(1:3).EQ.'N98') THEN  !Noble98 model
          CALL NOBLE98(T,Y,F)
        ELSE IF(SUB_NAME(1:3).EQ.'ATR') THEN  !atrial model
          CALL ATR(T,Y,F)
        ENDIF
        NFE=NFE+1
        DO J=1,NT_EQN
          D1(J)=0.0d0
!!!          WRITE(*,*) 'y,f',y(j),f(j)
        ENDDO
      ENDIF

!!!      WRITE(*,*) 'h',H
!!!      WRITE(*,*) 'nfe',nfe
!!!      WRITE(*,*) 'jflag',jflag
!!!      WRITE(*,*) 'hmin',hmin

      IF(IFLAG.LT.2) IFLAG=IFLAG-1
      DO WHILE(IFLAG.LT.1)
        IFLAG=IFLAG+1

C ***   rmax is maximum relative stepsize increase 10 usually,
C ***   but set to 2 after rejected step

        EPS=1.0d0
        EPSP1=2.0d0
        DO WHILE(EPSP1.GT.1.0d0)
          EPS=0.5d0*EPS
          EPSP1=EPS+1.0d0
        ENDDO
        EPS=2.0d0*EPS
        ALPH=0.5d0
        BET=0.7d0
        HMIN=100.0d0*EPS
        RMAX=10.0d0
        IF(JFLAG.EQ.2) RMAX=2.0d0
        IF(IFLAG.EQ.0) RMAX=1.0d5
        NREDO=0
        ERRMAX=0.005d0
        ERREPS=0.001d0*ERRMAX
        JFLAG=1
        M=0

!!!        WRITE(*,*) 'IFLAG,RMAX,NREDO,JFLAG,M',IFLAG,RMAX,NREDO,JFLAG,M

        DO WHILE((IFLAG.LT.2).AND.(M.LT.2))

C ***     Compute predicted values and store in array p
          DO J=1,NT_EQN
            P(J)=Y(J)+H*F(J)+H202*D1(J)
          ENDDO

          M=0
          TTST=T+H
          ERR=0.0d0

!!!          WRITE(*,*) 'ERR,TTST,M',ERR,TTST,M

          DO WHILE((ERR.LT.1.0d0).AND.(M.LT.2))

C ***       Compute derivatives using values in array p
            IF(DABS(P(1)).GT.200.0d0) IFLAG=5
            IF(SUB_NAME(1:2).EQ.'HH') THEN       !Hodgkin-Huxley eqns
              CALL HH(TTST,P,FP)
            ELSE IF(SUB_NAME(1:2).EQ.'BR') THEN  !Beeler-Reuter eqns
!             CALL BR(TTST,P,FP)
            ELSE IF(SUB_NAME(1:3).EQ.'RBR') THEN !Reduced Beeler-Reuter eqns
!             CALL RBR(TTST,P,FP)
            ELSE IF(SUB_NAME(1:4).EQ.'JRWP') THEN  !JRW (Princeton)
              CALL FCN_JRWP(1,2.0d0,2,P(1),P,FP,membrane_current)       
            ELSE IF(SUB_NAME(1:2).EQ.'LR') THEN  !Luo-Rudy eqns
              CALL LR(TTST,P,FP)
!              WRITE(*,'('' T='',E11.3,'' IFLAG='',I3)') T,IFLAG
!              WRITE(*,'('' Y( 1..10):'',10E11.3)') (Y(I),I= 1,10)
!              WRITE(*,'('' Y(11..20):'',10E11.3)') (Y(I),I=11,20)
!              WRITE(*,'('' Y(21..30):'',10E11.3)') (Y(I),I=21,30)
c              WRITE(99,'('' T='',E11.3,'' IFLAG='',I3)') T,IFLAG
c              WRITE(99,'('' Y( 1.. 7):'',10E11.3)') (Y(I),I= 1, 7)
c              WRITE(99,'('' Y( 8..14):'',10E11.3)') (Y(I),I= 8,14)
            ELSE IF(SUB_NAME(1:3).EQ.'JRW') THEN  !JRW eqns
              CALL JRW(TTST,P,FP)
            ELSE IF(SUB_NAME(1:2).EQ.'DM') THEN  !Distribution moment
              CALL FCN_DM(TTST,P,FP)
            ELSE IF(SUB_NAME(1:3).EQ.'N98') THEN !Noble98 model
              CALL NOBLE98(TTST,P,FP)
            ELSE IF(SUB_NAME(1:3).EQ.'ATR') THEN !atrial model
              CALL ATR(TTST,P,FP)
            ENDIF
            NFE=NFE+1

C ***       Compute first and second divided differences
            DO J=1,NT_EQN
              D2(J)=(FP(J)-F(J))/H
              DDIF(J)=(D2(J)-D1(J))/(H+HOLD)
!!!              WRITE(*,*) 'd2,dif',d2(j),ddif(j)
            ENDDO
!!!            WRITE(*,*) 'h,hold',H,HOLD
            IF(M.EQ.0) THEN
              H1=H303+HOLD*H202

C ***         Compute corrected values & store pending error test.
              DO J=1,NT_EQN
                P(J)=P(J)+H1*DDIF(J)
              ENDDO

C ***         Test estimated local error using l(infinity) norm.
              ERR=0.0d0
              HDEMX=0.5d0*H303/ERRMAX
              IF(IFLAG.EQ.0) HDEMX=0.5d0*HDEMX
!!!              WRITE(*,*) 'hdemx,h303,errmax',hdemx,h303,errmax
              DO J=1,NT_EQN
                ERRD=DABS(P(J))
                IF(ERRD.LT.ERREPS) ERRD=ERREPS
                ERTEST=HDEMX*DABS(DDIF(J))/ERRD
                IF(ERR.LT.ERTEST) ERR=ERTEST
!!!                WRITE(*,*) 'ERRD,ERREPS,ERTEST,DDIF,ERR'
!!!                WRITE(*,*) ERRD,ERREPS,ERTEST,DDIF(j),ERR
              ENDDO
            ENDIF
            M=M+1
          ENDDO

C ***     Error test fails. Reduce h and try again if h .GT. hmin.
C ***     If 6 straight reductions occur, assume discontinuity and start over.
          IF(M.LT.2) THEN
            JFLAG=2
            NREDO=NREDO+1
            IF(NREDO.GE.6) THEN
              IFLAG=0
              DISCONT=.TRUE.
!!!              WRITE(*,*) 'going to 1'
              GOTO 1
            ENDIF
            H=BET*H
            H202=0.5d0*H*H
            H303=THIRD*H*H*H
            IF(H.LE.HMIN) IFLAG=2
!!!            WRITE(*,*) 'H,h202,h303',h,h202,h303
!!!            WRITE(*,*) 'hmin,bet,iflag,jflag,nredo',hmin,bet,iflag,
!!!     '        jflag,nredo
          ENDIF
        ENDDO

        IF(IFLAG.LT.2) THEN
C ***     Step accepted. Update y array and prepare for next step.
          T=TTST
          NSTEP=NSTEP+1
          DO J=1,NT_EQN
            D1(J)=D2(J)
            F(J)=FP(J)
            Y(J)=P(J)
          ENDDO
          IF(SUB_NAME(1:2).EQ.'DN') THEN  !diFrancesco-Noble eqns
!           IF(DOP) THEN
!              WRITE(*,'('' T='',E11.3,'' IFLAG='',I3)') T,IFLAG
!              WRITE(*,'('' Y( 1..10):'',10E11.3)')
!     '    (Y(I),I= 1,10)
!              WRITE(*,'('' Y(11..20):'',10E11.3)')
!     '    (Y(I),I=11,20)
!              WRITE(*,'('' Y(21..30):'',10E11.3)')
!     '    (Y(I),I=21,30)
!           ENDIF
          ENDIF

C ***     Compute stepsize for next step.
          IF(ERR.GT.0) THEN
            Z=THIRD
            IF(IFLAG.EQ.0) Z=0.5d0
!!!            WRITE(*,*) 'ERR',ERR
            R=(ALPH/ERR)**Z
            IF(R.GT.RMAX) R=RMAX
          ELSE
            R=RMAX
          ENDIF
          HOLD=H
          H=R*H
          IF(H.LT.HMIN) H=HMIN
          IF(H.GE.HMAX) THEN
            H=HMAX
C ***       If h bumps against hmax too often hmax is increased.
            NBUMP=NBUMP+1
            IF(NBUMP.GE.20) THEN
              HMAX=2.0d0*HMAX
              NBUMP=0
            ENDIF
          ENDIF
          H202=0.5d0*H*H
          H303=THIRD*H*H*H
        ENDIF
      ENDDO

      RETURN
      END


      SUBROUTINE LR(t,Y,Func)

C#### Subroutine: LR
C###  Description:
C###    Computes RHS of Luo-Rudy equations.
C###    Voltage is in mV & time in ms.

      IMPLICIT NONE
c      INCLUDE 'cmiss$reference:b12.cmn'
      INCLUDE 'cmiss$reference:deoxs00.cmn'
      INCLUDE 'cmiss$reference:lr001.cmn'
      INCLUDE 'cmiss$reference:oxs005.cmn'
!     Parameter List
      REAL*8 t,Y(*),Func(*)
!     Common blocks
      REAL*8        I_Na,I_Ca,I_K,F_NSR,F_JSR
      COMMON /isum/ I_Na,I_Ca,I_K,F_NSR,F_JSR
cp      REAL*8        TRPN,CMDN,CSQN
cp      COMMON /buff/ TRPN,CMDN,CSQN
      REAL*8        I(19),ICa !DPN 3/10/97 - extending I()
      COMMON /curr/ I,ICa
C      REAL*8        TIMESCALE,TSTART,TEND,DT,TABT,TPS,TP
C      COMMON /times/ TIMESCALE,TSTART,TEND,DT,TABT,TPS,TP
      !REAL*8        Nai,Ki,Cai,CaNSR,CaJSR,Nao,Ko,Cao
      !COMMON /conc/ Nai,Ki,Cai,CaNSR,CaJSR,Nao,Ko,Cao
!     Local Variables
      !REAL*8 ACap,VolMyo,VolNSR,VolJSR,Farad
      REAL*8 VRatio_MJ, VRatio_MN, flux
      REAL*8 ALFA(7),BETA(7),I_stim,m,h,j,d,f,X !,Cm,V
      REAL*8 Iv_temp,RTONF
      REAL*8 LR_I(19),LR_ICa
      REAL*8 ISUMS(5),LR_I_Na,LR_I_Ca,LR_I_K,LR_F_NSR,LR_F_JSR
      INTEGER count

      RTONF = (Temp+273.d0)*86.16d-3 !mV

! Assign names to variables
C dpn 17/02/98      V = Y( 1)
      !?? DPN 3/10/97 - wrong indexing ??
      m = Y(8)     !Y( 2)
      h = Y(7)     !Y( 3)
      j = Y(9)     !Y( 4)
      d = Y(5)     !Y( 5)
      f = Y(6)     !Y( 6)
      X = Y(4)     !Y( 7)

      VRatio_MJ = VolJSR/VolMyo
      VRatio_MN = VolNSR/VolMyo

! Setting the stimulus current between t=5.0 and t=5.5ms
C??? DPN 7/10/97 - remove hard-coded times
C      IF ((t.GE.Tstim).AND.(t.LE.Tstim+0.5d0)) THEN
C        I_stim=-1.d2
      IF ((t.GE.TPS).AND.(t.LE.TPS+TP)) THEN
        I_stim=Istim
      ELSE
        I_stim=0.d0
      ENDIF
       
! Compute rate constants
      CALL L_R_RATES(Y,ALFA,BETA,RTONF)

! Compute ionic currents
      CALL LR_CURRENTS(t,Y,ALFA,BETA,RTONF,LR_I,ISUMS,LR_ICa)
      LR_I_Na = ISUMS(1)
      LR_I_Ca = ISUMS(2)
      LR_I_K = ISUMS(3)
      LR_F_NSR = ISUMS(4)
      LR_F_JSR = ISUMS(5)
      IF (SINGLE_CELL) THEN
        DO count=1,19
          I(count) = LR_I(count)
        ENDDO
        ICa = LR_ICa
        I_Na = ISUMS(1)
        I_Ca = ISUMS(2)
        I_K = ISUMS(3)
        F_NSR = ISUMS(4)
        F_JSR = ISUMS(5)
      ENDIF

      flux = LR_I(13)*VRatio_MJ+(-LR_I(14)+LR_I(15))*VRatio_MN      

! Compute o.d.e. RHS 
C      Func( 1) =  -(I_stim+I(1)+I(2)+I(3)+I(6)+I(8)+I(12))/Cm
C      Func( 2) =  ALFA(1)*(1.d0-m) - BETA(1)*m
C      Func( 3) =  ALFA(2)*(1.d0-h) - BETA(2)*h
C      Func( 4) =  ALFA(3)*(1.d0-j) - BETA(3)*j
C      Func( 5) =  ALFA(4)*(1.d0-d) - BETA(4)*d
C      Func( 6) =  ALFA(5)*(1.d0-f) - BETA(5)*f
C      Func( 7) =  ALFA(6)*(1.d0-X) - BETA(6)*X
C      FUNC( 8) =  -(I_Na*ACap)/(VolMyo*Farad)
C      FUNC( 9) =  -(I_Ca*ACap)/(VolMyo*2.d0*Farad)+flux
C      FUNC(10) =  -(I_K*ACap) /(VolMyo*Farad)
C      FUNC(11) =  -F_JSR
C      FUNC(12) =  -F_NSR
C      FUNC(13) =  -((ICa+I(9)+I(10))*1.d-3)/2.d0
C      FUNC(14) = 0 ![Ca]o constant

C *** DPN 15/04/98 - adding switch to each current in the total
      Iv_temp = ISWTCH(14)*LR_I(4)+ISWTCH(21)*LR_I(5)+ISWTCH(18)*LR_I(7)
     '  +ISWTCH(26)*LR_I(9)+ISWTCH(2)*LR_I(10)+ISWTCH(4)*LR_I(11)
      Iv_temp = ISWTCH(24)*Iv_temp ! turns all time indep. currents off
      Func( 1) =  -(I_stim+ISWTCH(16)*LR_I(1)+ISWTCH(5)*LR_I(2)
     '  +ISWTCH(13)*LR_I(3)+ISWTCH(17)*LR_I(6)+ISWTCH(22)*LR_I(8)
     '  +Iv_temp)/Cm
      Func( 8) =  ALFA(1)*(1.d0-m) - BETA(1)*m
      Func( 7) =  ALFA(2)*(1.d0-h) - BETA(2)*h
      Func( 9) =  ALFA(3)*(1.d0-j) - BETA(3)*j
      Func( 5) =  ALFA(4)*(1.d0-d) - BETA(4)*d
      Func( 6) =  ALFA(5)*(1.d0-f) - BETA(5)*f
      Func( 4) =  ALFA(6)*(1.d0-X) - BETA(6)*X
      FUNC( 2) =  -(LR_I_Na*ACap)/(VolMyo*Farad)
      FUNC( 3) =  -(LR_I_Ca*ACap)/(VolMyo*2.d0*Farad)+flux
      FUNC(10) =  -(LR_I_K*ACap) /(VolMyo*Farad)
      FUNC(12) =  -LR_F_JSR
      FUNC(11) =  -LR_F_NSR
      FUNC(13) =  -((LR_ICa+LR_I(9)+LR_I(10))*1.d-3)/2.d0
      FUNC(14) = 0.0d0 ![Ca]o constant

      RETURN
      END 


      SUBROUTINE LR_CURRENTS(t,Y,ALFA,BETA,RTONF,I,ISUMS,ICa)

C#### Subroutine: LR_CURRENTS
C###  Description:
C###    Computes ionic currents for Luo-Rudy model

C??? DPN 3/10/97 - wrong indexing of Y array ???

C**** Potentials, gating variables, internal concentrations
C**** Y( 1) = V                                             Y( 1)
C**** Y( 2) = m   (Na V- & t- dep. activation gate)         Y( 8)
C**** Y( 3) = h   (Na fast V- & t- dep. inactivation gate)  Y( 7)
C**** Y( 4) = j   (Na slow V- & t- dep. inactivation gate)  Y( 9)
C**** Y( 5) = d   (Ca V- & t- dep. activation gate)         Y( 5)
C**** Y( 6) = f   (Ca V- & t- dep. inactivation gate)       Y( 6)
C****         fCa (Ca Ca-dep. inactivation gate)    
C**** Y( 7) = X   (K  V- & t- dep. activation gate)         Y( 4)
C****         Xi  (K  V-dep. inactivation gate)
C**** Y( 8) = internal Na concentration                     Y( 2)
C**** Y( 9) = internal Ca concentration                     Y( 3)
C**** Y(10) = internal K concentration                      Y(10)
C**** Y(11) = CaJSR concentration                           Y(12)
C**** Y(12) = CaNSR concentration                           Y(11)
C**** Y(13) = internal Ca concentration at T tubules due to 
C****         L-type channels                               Y(13)
C**** Y(14) = extracellular [Ca] (constant???)

C**** Currents (uA/uF)
C**** I( 1) = INa
C**** I( 2) = ICa_t
C**** I( 3) = IK
C**** I( 4) = IK1
C**** I( 5) = IKp
C**** I( 6) = INaCa
C**** I( 7) = INaK
C**** I( 8) = InsCa
C**** I( 9) = IpCa
C**** I(10) = ICab
C**** I(11) = INab
C**** I(12) = Iv
C**** I(13) = Irel (flux)
C**** I(14) = Iup (flux)
C**** I(15) = Ileak (flux)
C**** I(16) = Itr (flux)
C??? DPN 3/10/97 - extend I() to include ICa, ICaK, and ICaNa (cpts of 
C                  I(Ca(L)) = ICa_t)
C**** I(17) = ICa
C**** I(18) = ICaK
C**** I(19) = ICaNa

C**** total current for each ions
C**** I_Na  = sum of currents containing Na ions
C**** I_Ca  = sum of currents containing Ca ions
C**** I_K   = sum of currents containing K ions
C**** F_JSR = sum of fluxes moving into JSR
C**** F_NSR = sum of fluxes moving into NSR

      IMPLICIT NONE
c      INCLUDE 'cmiss$reference:b12.cmn'
      INCLUDE 'cmiss$reference:deoxs00.cmn'
      INCLUDE 'cmiss$reference:oxs005.cmn'
      INCLUDE 'cmiss$reference:lr001.cmn'
!     Parameter List
      REAL*8 t,Y(*),ALFA(*),BETA(*),RTONF,I(19),ISUMS(5),ICa
!     Common blocks
cp      REAL*8        I_Na,I_Ca,I_K,F_NSR,F_JSR
cp      COMMON /isum/ I_Na,I_Ca,I_K,F_NSR,F_JSR
cp      REAL*8        I(19),ICa !DPN 3/10/97 - extending I()
cp      COMMON /curr/ I,ICa
cp      REAL*8        CaNSR,CaJSR,Ko !,Nai,Ki,Cai,Nao,Cao
cp      COMMON /conc/ CaNSR,CaJSR,Ko !,Nai,Ki,Cai,Nao,Cao
C      REAL*8        TIMESCALE,TSTART,TEND,DT,TABT,TPS,TP
C      COMMON /times/ TIMESCALE,TSTART,TEND,DT,TABT,TPS,TP
      REAL*8        Cai_on,Cai2
      COMMON /Cai_/ Cai_on,Cai2

      REAL*8                   IICa,IICaK,IICaNa,IINa,IIK,IIK1,fca
      COMMON /LR_MAX_CURRENTS/ IICa,IICaK,IICaNa,IINa,IIK,IIK1,fca

!     Local Variables
      REAL*8 V,m,h,j,d,f,X
      REAL*8 ENa!,GGNa
      REAL*8 ICaK,ICaNa
      REAL*8 EK,GGK,Xi!,PNaK
      REAL*8 EK1,GGK1,K10
      REAL*8 EKp,Kp
      REAL*8 fNaK,sigma
      REAL*8 EnsCa,Vns,InsNa,IInsNa,InsK,IInsK
      REAL*8 ECaN!, GGCab
      REAL*8 ENaN!, GGNab
      REAL*8 CaJSRbuffer !!SMAR009 23/12/98 Cabuffer,
      REAL*8 Grel,GGrel
c      REAL*8 Cai_on,Cai2,Caith, Kmrel, t_CICR
C DPN 3/11/97 - move Cai_on and Cai2 to a common block to allow them to
C               be reset between runs, without exiting program.
      REAL*8 t_CICR
      REAL*8 Kleak
      REAL*8 VRatio_JN
      REAL*8 LR_Cai,LR_Nai,LR_Ki
      REAL*8       CSQN !SMAR009 21/12/98 ,TRPN,CMDN
cp      COMMON /buff/ CSQN,TRPN,CMDN
      REAL*8        LR_CaNSR,LR_CaJSR,LR_Ko
      REAL*8        I_Na,I_Ca,I_K,F_NSR,F_JSR

! Initialising the volume parameters
      !VolJSR = 0.182d-6    !uL
      !VolNSR = 2.098d-6    !uL
      VRatio_JN = VolNSR/VolJSR       !volume ratio between JSR and NSR

! Initialising the thermodynamic parameters
      V     = Y(1)                    !membrane potential (mV)
      !Temp  = 37.d0                   !temperature (deg C)
      !RTONF = (Temp+273.d0)*86.16d-3  !RT/F (mV)
      !Farad = 96.5d3                  !C/mol

! Standard ionic concentrations in mmol/L
      LR_Nai   = Y(2)   !Y(8)       ![Na]i
      !Nao   = 140.d0     ![Na]o
      LR_Cai   = Y(3)   !Y(9)       ![Ca]i
      !Cao   = 1.8d0      ![Ca]o
      LR_Ki    = Y(10)      ![K]i
      LR_Ko    = Kb !DPN - fix for multiprocessor
      !Ko    = 5.4d0      ![K]o
 

! Ca buffers in the myoplasm
      !??? DPN - move to common block ???
      !TRPNTRPN = 70.d-3               !mmol/L
      !CMDNCMDN = 50.d-3               !mmol/L
      !KmTRPN   = 0.5d-3               !mmol/L
      !KmCMDN   = 2.38d-3              !mmol/L
c      TRPN     = TRPNTRPN*(LR_Cai/(LR_Cai+KmTRPN))
c      CMDN     = CMDNCMDN*(LR_Cai/(LR_Cai+KmCMDN))


! Calcium concentration after taking the buffer TRPN and CMDN into 
! account
C dpn 03/06/98 - uncomment
!      Cabuffer = Cai-TRPN-CMDN
!      IF(t.GE.Tstim)THEN
!        IF(Cabuffer.GT.1.d-6)THEN
!          Cai = Cabuffer
!        ELSE 
!          Cai = 1.d-6
!        ENDIF
!      ENDIF
c      Cabuffer = Cai-TRPN-CMDN
c      IF(t.GE.TPS)THEN
c        IF(Cabuffer.GT.1.d-6)THEN
c          Cai = Cabuffer
c        ELSE
c          Cai = 1.d-6
c        ENDIF
c      ENDIF
C dpn 03/06/98

! Ca fluxes in the sarcoplasmic reticulum
      LR_CaJSR = Y(12)  !Y(11)                   !mmol/L
      LR_CaNSR = Y(11)  !Y(12)                   !mmol/L


! Ca buffer in JSR and CSQN
      !??? DPN - move to common block ???
      !KmCSQN    = 0.8d0               !mmol/L
      !CSQNCSQN  = 10.d0               !mmol/L
      CSQN      = CSQNCSQN*(LR_CaJSR/(LR_CaJSR+KmCSQN))


! Ca flux in the SR after taking CSQN buffering into account 
C dpn 03/06/98 - uncomment
!      CaJSRbuffer = CaJSR-CSQN
!      IF(t.GE.Tstim)THEN
!        IF(CaJSRbuffer.GT.1.d-3)THEN
!          CaJSR = CaJSRbuffer
!        ELSE
!          CaJSR = 1.d-3
!        ENDIF
!      ENDIF
      CaJSRbuffer = LR_CaJSR-CSQN
      IF(t.GE.TPS)THEN
        IF(CaJSRbuffer.GT.1.d-3)THEN
          LR_CaJSR = CaJSRbuffer
        ELSE
          LR_CaJSR = 1.d-3
        ENDIF
      ENDIF
C dpn 03/06/98

! Ionic currents in the sarcolemma
! Fast Na current INa
      m = Y(8)  !Y(2)
      h = Y(7)  !Y(3)
      j = Y(9)  !Y(4)

      !GGNa = 16.d0                    !Na conductance (millisiemens/uF)
      ENa  = RTONF*dlog(Nao/LR_Nai)      !Na reversal potential (mV)  
C *** DPN 29 April 1998 - grab max INa
c      I(1) = GGNa*m*m*m*h*j*(V-ENa)
      IINa = GGNa * (V-ENa)
      I(1) = m*m*m * h * j * IINa

! L-type Ca current ICa_t
      d = Y(5)
      f = Y(6)

      KmCa = KmCa_L !0.6d-3                   !mmol/L
      fCa  = 1.d0/(1.d0+(LR_Cai/KmCa)**2)!Ca-dep. inact.n gate of L-type Ca

      !PCa  = 5.4d-4                   !permeablity of membrane to Ca (cm/s)
      !gCai = 1.d0                     !activity coefficent of Cai
      !gCao = 0.341d0                  !activity coefficent of Cao
      IICa = PCa*4.d0*(V/RTONF)*Farad*(gCai*LR_Cai*dexp(2.d0*V/RTONF)
     '       -gCao*Cao)/(dexp(2.d0*V/RTONF)-1.d0)

      !PNa    = 6.75d-7                !permeablity of membrane to Na (cm/s)
      !gNai   = 0.75d0                 !activity coefficent of Nai
      !gNao   = 0.75d0                 !activity coefficent of Nao
      IICaNa = PNa*(V/RTONF)*Farad*(gNai*LR_Nai*dexp(V/RTONF)-gNao*Nao)
     '                           /(dexp(V/RTONF)-1.d0)

      !PK     = 1.93d-7                !permeablity of membrane to K (cm/s)
      !gKi    = 0.75d0                 !activity coefficent of Ki
      !gKo    = 0.75d0                 !activity coefficent of Ko
      IICaK  = PK*(V/RTONF)*Farad*(gKi*LR_Ki*dexp(V/RTONF)-gKo*LR_Ko)
     '                           /(dexp(V/RTONF)-1.d0)

      ICa   = d*f*fCa*IICa
      ICaK  = d*f*fCa*IICaK
      ICaNa = d*f*fCa*IICaNa
      
      ! DPN 3/10/97 - include ICa,ICaK,ICaNa in I()
      I(17) = ICa
      I(18) = ICaK
      I(19) = ICaNa

      I(2)   = ICa + ICaK + ICaNa

! Time-dep. K current IK
      X    = Y(4)  !Y(7)                     !V- & t- dep. K-activation gate
      Xi   = 1.d0/(1.d0+dexp((V-56.26d0)/32.1d0)) 
     '                                !V-dep. K-inactivation gate
      GGK  = 0.282d0*dsqrt(LR_Ko/5.4d0)  !millisiemens/uF
      !PNaK = 0.01833d0                !permeablity ratio of Na:K
      EK   = RTONF*dlog((LR_Ko+PNaK*Nao)/(LR_Ki+PNaK*LR_Nai)) 
                                      !K reversal potential
C *** DPN 29 April 1998 - grab the max I(K)
c      I(3) = GGK*Xi*X*X*(V-EK)
      IIK = GGK*(V-EK)
      I(3) = IIK * Xi * X * X

! Time-indep. K current IK1
      GGK1 = 0.75d0*dsqrt(LR_Ko/5.4d0)   !millisiemens/uF
      EK1  = RTONF*dlog(LR_Ko/LR_Ki)        !K1 reversal potential  
      K10  = ALFA(7)/(ALFA(7)+BETA(7))
C *** DPN 29 April 1998 - grab max I(K1)
c      I(4) = GGK1*K10*(V-EK1)
      IIK1 = GGK1 * (V - EK1)
      I(4) = IIK1 * K10

! Plateau K current IKp
      !GGKp = 0.0183d0                !millisiemens/uF
      EKp  = EK1                      !Kp reversal potential  
      Kp   = 1.d0/(1.d0+dexp((7.488d0-V)/5.98d0))
      I(5) = GGKp*Kp*(V-EKp)

! Na-Ca exchanger current INaCa
      !DPN 7/10/97 - moved to common block
      !kNaCa = 2000.d0   !scaling factor of INaCa (uA/uF)
      !KmNa  = 87.5d0    !half sat.n conc of Na channel (mmol/L)
      KmCa  = KmCa_NaCa !1.38d0    !half sat.n conc of Ca channel (mmol/L)
      !eta   = 0.35d0    !position of the energy barrier for V-dep. of INaCa
      !ksat  = 0.1d0     !saturation factor of INaCa at very -ve potentials
      I(6)  = kNaCa*(dexp(eta*V/RTONF)*LR_Nai**3*Cao
     '             -dexp((eta-1.d0)*V/RTONF)*Nao**3*LR_Cai)
     '             /((KmNa**3+Nao**3)*(KmCa+Cao)
     '             *(1.d0+ksat*dexp((eta-1.d0)*V/RTONF)))

! Na-K pump INaK
      !IINaK = 1.5d0                   !uA/uF
      !KmNai = 10.d0                   !mmol/L
      !KmKo  = 1.5d0                   !mmol/L
      sigma = (1.d0/7.d0)*(dexp(Nao/67.3d0)-1.d0) 
     '                                ![Na]o-dependence factor of fNak
      fNaK  = 1.d0/(1.d0+0.1245d0*dexp(-0.1d0*V/RTONF)
     '        +0.0365d0*sigma*dexp(-V/RTONF))
      I(7)  = IINaK*fNaK*LR_Ko/((LR_Ko+KmKo)*(1.d0+
     '  dsqrt((KmNai/LR_Nai)**2)))

! Non-specific Ca-activated current InsCa
      !KmnsCa = 1.2d-3                 !mmol/L
      !PnsCa  = 1.75d-7                !cm/s
      EnsCa  = RTONF*dlog((LR_Ko+Nao)/(LR_Ki+LR_Nai))
      Vns    = V-EnsCa
      IInsK  = PnsCa*(Vns/RTONF)*Farad*(gNai*LR_Nai*dexp(Vns/RTONF)
     '                          -gNao*Nao)/(dexp(Vns/RTONF)-1.d0)
      InsK   = IInsK /(1.d0+(KmnsCa/LR_Cai)**3) 
      IInsNa = PnsCa*(Vns/RTONF)*Farad*(gKi*LR_Ki*dexp(Vns/RTONF)
     '                          -gKo*LR_Ko)/(dexp(Vns/RTONF)-1.d0)
      InsNa  = IInsNa/(1.d0+(KmnsCa/LR_Cai)**3)
      I(8)   = InsK + InsNa

! Sarcolemmal Ca pump IpCa
      !KmpCa  = 0.5d-3                 !mmol/L
      !IIpCa  = 1.15d0                 !uA/uF
      I(9)   = IIpCa*(LR_Cai/(KmpCa+LR_Cai))

! Ca background current ICab
      ECaN   = 0.5d0*RTONF*dlog(Cao/LR_Cai) !Nernst potential of Ca
      !GGCab  = 0.003016d0                !millisiemens/uF
      I(10)  = GGCab*(V-ECaN)

! Na background current INab
      ENaN   = ENa                    !Nernst potential of Na
      !GGNab  = 0.00141d0              !millisiemens/uF
      I(11)  = GGNab*(V-ENaN)

! Total time independent current
      I(12)  = I(4) + I(5) + I(7) + I(9) + I(10) + I(11)


! Compute Cai change at 2 ms after onset of stimulus
      !IF((t.GE.TPS).AND.(t.LT.0.1d0+TPS))THEN
      IF((t.GE.TPS).AND.(t.LT.TABT+TPS))THEN
        Cai_on = Y(13)
        IF(DEBUG) write (*,*) 'Cai_on = ', Cai_on
      ENDIF

      IF((t.GE.2.d0+TPS).AND.(t.LT.2.d0+TABT+TPS))THEN
        Cai2 = Y(13)-Cai_on
        IF(DEBUG) write (*,*) 'Cai2 = ', Cai2
      ENDIF

! Ca induced Ca release of JSR
      !Kmrel    = 0.8d-3               !mmol/L
      !Tau_on   = 2.d0                 !ms
      !Tau_off  = 2.d0                 !ms
      !Caith    = 0.18d-3              !mmol/L
      IF(Cai2.GT.Caith) THEN
c        GGrel = 60.d0                 !ms-1 for action potential simulations
        GGrel = GGrel_
      ELSE
        GGrel = 0.d0
      ENDIF
      !time of Ca-induced Ca-release
c dpn 03/06/98      t_CICR  = t-TPS
      t_CICR  = t-(TPS+2.0d0)
      Grel  = GGrel*((Cai2-Caith)/(Kmrel+Cai2-Caith))
     '       *(1.d0-dexp(-t_CICR/Tau_on))*dexp(-t_CICR/Tau_off)
      I(13) = Grel*(LR_CaJSR-LR_Cai)

! Ca release of JSR under Ca-overload conditions
!      CSQNth  = 0.7d0
!      IF(CSQN.GE.CSQNth) THEN
!         GGrel = 4.d0                !ms-1
!      ELSE
!         GGrel = 0.d0                !ms-1
!      ENDIF
!      Grel    = GGrel*(1.d0-dexp(-t/Tau_on))*dexp(-t/Tau_off)
!      I(13)   = Grel*(CaJSR-Cai)


! Ca uptake and leakage of NSR Iup and Ileak
      !Kmup       = 0.92d-3                !mmol/L
      !CaNSRCaNSR = 15.d0                  !mmol/L
      !IIup       = 0.005d0                !mmol/L/ms
      Kleak      = IIup/CaNSRCaNSR        !ms-1
      I(15)      = Kleak*LR_CaNSR            !mmol/L/ms
      I(14)      = IIup*(LR_Cai/(LR_Cai+Kmup))  !mmol/L/ms


! Translocation of Ca ions from NSR to JSR Itr
      !Tau_tr = 180.d0                     !ms
      I(16)  = (LR_CaNSR-LR_CaJSR)/Tau_tr       !mmol/L/ms

! Computes current carrying the five ions
c      I_Na  = I(1)+ICaNa+I(6)+I(7)+InsNa+I(11)
      I_Na  = I(1)+ICaNa+3.0d0*I(6)+3.0d0*I(7)+InsNa+I(11)
      I_Ca  = ICa-I(6)+I(9)+I(10)
c      I_K   = ICaK+I(3)+I(4)+I(5)+I(7)+InsK
      I_K   = ICaK+I(3)+I(4)+I(5)-2.0d0*I(7)+InsK
      F_NSR = -I(14)+I(15)+I(16)
      F_JSR = I(13)-I(16)*VRatio_JN

      ISUMS(1) = I_Na
      ISUMS(2) = I_Ca
      ISUMS(3) = I_K
      ISUMS(4) = F_NSR
      ISUMS(5) = F_JSR

      RETURN
      END 


      SUBROUTINE L_R_RATES(Y,ALFA,BETA,RTONF)

C#### Subroutine: L_R_RATES
C###  Description:
C###    <html><pre>
C###    Computes the rate coefficients for the gating variables.
C###    ALFA (1) = m
C###      "   2  = h
C###      "   3  = j
C###      "   4  = d
C###      "   5  = f
C###      "   6  = X
C###      "   7  = K1 
C###    BETA (1) = m
C###      "   2  = h
C###      "   3  = j
C###      "   4  = d
C###      "   5  = f
C###      "   6  = X
C###      "   7  = K1
C###    </pre></html>

      IMPLICIT NONE
c      INCLUDE 'cmiss$reference:b12.cmn'
c      INCLUDE 'cmiss$reference:deoxs00.cmn'
      INCLUDE 'cmiss$reference:lr001.cmn'
!     Parameter List
      REAL*8 ALFA(*),BETA(*),Y(*),RTONF
!     Common blocks
      REAL*8        CaNSR,CaJSR,Ko !,Nai,Ki,Cai,Nao,Cao
      COMMON /conc/ CaNSR,CaJSR,Ko !,Nai,Ki,Cai,Nao,Cao
!     Local Variables
      REAL*8 d0,f0,Tau_d,Tau_f
      REAL*8 V,EK1,Ki ! SMAR009 removed ->Ko !,RTONF,TEMP,Ko,Ki
      
      V = Y(1)
      !Temp  = 37.d0                   !temperature (deg C)
      !RTONF = (Temp+273.d0)*86.16d-3  !RT/F

! Standard ionic concentrations
      Ki   = Y(10)                    ![K]i
      !Ko   = 5.4d0                    ![K]o

! Compute the K1 reversal potential (mV)
      EK1 = RTONF*dlog(Ko/Ki)         !K1 reversal potential (mV)

! Fast Na channel gates m,h,j
      ALFA(1) = 0.32d0*(V+47.13d0)/(1.d0-dexp(-0.1d0*(V+47.13d0)))   !m
      BETA(1) = 0.08d0*dexp(-V/11.d0)                                  !m
      IF(V.GE.-40.d0) THEN       
        ALFA(2) = 0.d0                                                 !h
        BETA(2) = 1.d0/(0.13d0*(1.d0+dexp((V+10.66d0)/(-11.1d0))))     !h
        ALFA(3) = 0.d0                                                 !j
        BETA(3) = 0.3d0*dexp(-2.535d-7*V)/(1.d0+dexp(-0.1d0*(V+32.d0)))!j
      ELSE IF(V.LT.-40.d0) THEN
        ALFA(2) = 0.135d0*dexp((80.d0+V)/(-6.8d0))                     !h
        BETA(2) = 3.56d0*dexp(0.079d0*V)+(3.1d5*dexp(0.35d0*V))        !h
        ALFA(3) = (-1.2714d5*dexp(0.2444d0*V)-3.474d-5*
     '            dexp(-0.04391d0*V))*(V+37.78d0)/(1.d0+dexp(0.311d0*
     '            (V+79.23d0)))                                        !j
        BETA(3) = 0.1212d0*dexp(-0.01052d0*V)/(1.d0+dexp(-0.1378d0
     '            *(V+40.14d0)))                                       !j
      ENDIF

! L-type Ca channel gates d,f
      d0    = 1.d0/(1.d0+dexp(-(V+10.d0)/6.24d0))
      Tau_d = d0*(1.d0-dexp(-(V+10.d0)/6.24d0))/(0.035d0*(V+10.d0))
      f0    = (1.d0/(1.d0+dexp((V+35.06d0)/8.6d0)))+(0.6d0/
     '        (1.d0+dexp((50.d0-V)/20.d0)))
      Tau_f = 1.d0/(0.0197d0*dexp(-((0.0337d0*(V+10.d0))**2))+0.02d0)

      ALFA(4) = d0/Tau_d              !d
      BETA(4) = (1.d0-d0)/Tau_d       !d
      ALFA(5) = f0/Tau_f              !f
      BETA(5) = (1.d0-f0)/Tau_f       !f

! Time-dep. K channel gate X
      ALFA(6) = 7.19d-5*(V+30.d0)/(1.d0-dexp(-0.148d0*(V+30.d0)))     ! X
      BETA(6) = 1.31d-4*(V+30.d0)/(-1.d0+dexp(0.0687d0*(V+30.d0)))    ! X

! Time-indep. K channel gate K1
      ALFA(7) = 1.02d0/(1.d0+dexp(0.2385d0*(V-EK1-59.215d0)))
      BETA(7) = (0.49124d0*dexp(0.08032d0*(V-EK1+5.476d0))
     '          +dexp(0.06175d0*(V-EK1-594.31d0)))
     '          /(1.d0+dexp(-0.5143d0*(V-EK1+4.753d0)))

      RETURN
      END 


      SUBROUTINE GATES(Y,F)

C**** Computes time constants and activation variables.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:deoxs00.cmn'
      INCLUDE 'cmiss$reference:oxs001.cmn'
      INCLUDE 'cmiss$reference:oxs003.cmn'
      INCLUDE 'cmiss$reference:oxs004.cmn'
!     Parameter List
      REAL*8 F(*),Y(*)
!     Local Variables
      REAL*8 EM

      PMODE=3
      PLOTT=.TRUE.
      EM=-100.0
      DO WHILE(EM.LT.40.0)
        Y(1)=EM
        CALL RATES(Y,F)
        Z1=Y(5)/(SPEED(5)*ALPHA(5))
C    ETC ETC
      ENDDO

      RETURN
      END

      SUBROUTINE NOBLE98(t,Y,F)

C#### Subroutine: NOBLE98
C###  Description:
C###    Calculates the RHS of the Noble '98 set of ODE's

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:deoxs00.cmn'
      INCLUDE 'cmiss$reference:oxs005.cmn'
      !passed variables
      REAL*8 t,Y(*),F(*)
      !common blocks
      REAL*8 INa98,IpNa98,IbNa98,IK198,IKr98,IKs98,IbK98,IKATP98,
     '  IKNa98,IK98,ICaLCa98,ICaLK98,ICaLNa98,ICaLCads98,ICaLKds98,
     '  ICaLNads98,IbCa98,IKACH98,INaK98,INaCa98,INaCads98,Iup98,
     '  Itr98,Irel98,XSRrel98,Ito98
      COMMON /OH5_CURRENTS/
     '  INa98,IpNa98,IbNa98,IK198,IKr98,IKs98,IbK98,IKATP98,
     '  IKNa98,IK98,ICaLCa98,ICaLK98,ICaLNa98,ICaLCads98,ICaLKds98,
     '  ICaLNads98,IbCa98,IKACH98,INaK98,INaCa98,INaCads98,Iup98,
     '  Itr98,Irel98,XSRrel98,Ito98
      REAL*8 OH5_Vcell,OH5_Vi,OH5_Vds,OH5_VSRup
      COMMON /OH5_VOLUMES/
     '  OH5_Vcell,OH5_Vi,OH5_Vds,OH5_VSRup
      REAL*8 OH5_Istim
      COMMON /OH5_STIMULUS/
     '  OH5_Istim
      REAL*8    OH5_kcachoff,OH5_kmk1,OH5_kmK,
     '          OH5_kmNa,OH5_kNaCa,OH5_kmCa,
     '          OH5_k1,OH5_k2,OH5_k3,OH5_k4,
     '          OH5_kATP,OH5_KD,OH5_kkNa,
     '          OH5_kdsoff,OH5_kmCads,OH5_kdecay,
     '          OH5_Rdecay,OH5_INaKmax
      COMMON /OH5_PARAM1/
     '  OH5_kcachoff,OH5_kmk1,OH5_kmK,
     '  OH5_kmNa,OH5_kNaCa,OH5_kmCa,
     '  OH5_k1,OH5_k2,OH5_k3,OH5_k4,
     '  OH5_kATP,OH5_KD,OH5_kkNa,
     '  OH5_kdsoff,OH5_kmCads,OH5_kdecay,
     '  OH5_Rdecay,OH5_INaKmax
      REAL*8    OH5_GNa,OH5_GbK,OH5_GK1,OH5_GbNa,
     '          OH5_GpNa,OH5_GbCa,
     '          OH5_GKr1,OH5_GKr2,OH5_GKs,OH5_GKATP,OH5_GKACh,
     '          OH5_PCa,OH5_PCaK,OH5_PCaNa,OH5_PKNa
      COMMON /OH5_PARAM2/
     '  OH5_GNa,OH5_GbK,OH5_GK1,OH5_GbNa,
     '  OH5_GpNa,OH5_GbCa,
     '  OH5_GKr1,OH5_GKr2,OH5_GKs,OH5_GKATP,OH5_GKACh,
     '  OH5_PCa,OH5_PCaK,OH5_PCaNa,OH5_PKNa
      REAL*8    OH5_Cab,OH5_Ko,OH5_Nao,
     '          OH5_FractICaL,OH5_FractINaCa,
     '          OH5_DIFFCa,OH5_dNaCa,OH5_gama,OH5_nNaCa,
     '          OH5_kcyCa,OH5_kxcs,OH5_ksrCa,OH5_Cm,
     '          OH5_ALPHA12,OH5_BETA12
      COMMON /OH5_PARAM3/
     '  OH5_Cab,OH5_Ko,OH5_Nao,
     '  OH5_FractICaL,OH5_FractINaCa,
     '  OH5_DIFFCa,OH5_dNaCa,OH5_gama,OH5_nNaCa,
     '  OH5_kcyCa,OH5_kxcs,OH5_ksrCa,OH5_Cm,
     '  OH5_ALPHA12,OH5_BETA12
      REAL*8    OH5_F,OH5_R,OH5_T,
     '          OH5_Vecs,OH5_radius,OH5_length,
     '          OH5_Vup,OH5_Vrel,OH5_kmCa2,
     '          OH5_Mtrop,OH5_Ctrop,OH5_alfatrop,
     '          OH5_betatrop,OH5_gamatropSL,
     '          OH5_gamaSRSL,OH5_sacSL
      COMMON /OH5_PARAM4/
     '  OH5_F,OH5_R,OH5_T,
     '  OH5_Vecs,OH5_radius,OH5_length,
     '  OH5_Vup,OH5_Vrel,OH5_kmCa2,
     '  OH5_Mtrop,OH5_Ctrop,OH5_alfatrop,
     '  OH5_betatrop,OH5_gamatropSL,
     '  OH5_gamaSRSL,OH5_sacSL
      REAL*8    OH5_G_Ca_stretch,OH5_G_Na_stretch,OH5_G_K_stretch,
     '          OH5_G_Ns_stretch,OH5_G_An_stretch,
     '          OH5_E_Ns_stretch,OH5_E_An_stretch
      COMMON /OH5_PARAM5/
     '  OH5_G_Ca_stretch,OH5_G_Na_stretch,OH5_G_K_stretch,
     '  OH5_G_Ns_stretch,OH5_G_An_stretch,
     '  OH5_E_Ns_stretch,OH5_E_An_stretch
      REAL*8    OH5_alfaS,OH5_alfaV,
     '          OH5_alfaACh,OH5_betaACh,
     '          OH5_SLhst,OH5_IThst
      COMMON /OH5_PARAM6/
     '  OH5_alfaS,OH5_alfaV,
     '  OH5_alfaACh,OH5_betaACh,
     '  OH5_SLhst,OH5_IThst
      REAL*8 OH5_kdsdecay,OH5_Jsrleak
      COMMON /OH5_PARAM7/
     '  OH5_kdsdecay,OH5_Jsrleak
      !local variables
      REAL*8 V,mNa,hNa,dCa,fCa,f2Ca,xr1,xr2,xs,xACH,Cads,
     '  CaCALMOD,CaTROP,factivator,fproduct,rto,qto !SMAR009 21/12/98 Caup,Carel,
C *** stored in deoxs00.cmn
      REAL*8 OH5_Nai,OH5_Ki,OH5_Cao,OH5_Cai
      REAL*8 RTF,ENa,EK,EKs,ECa,Emh,POTS(5)
      REAL*8 alfa_mNa,beta_mNa,alfa_hNa,beta_hNa,alfa_dCa,beta_dCa,
     '  alfa_fCa,beta_fCa,alfa_xr1,beta_xr1,alfa_xr2,beta_xr2,
     '  alfa_xs,beta_xs,alfa_xACH,beta_xACH,alfa_rto,beta_rto
      REAL*8 ALPHA(9),BETA(9)
      REAL*8 INa,IpNa,IbNa,IK1,IbK,IKATP, !SMAR009 21/12/98 ,IKr,IKs,
     '  IKNa,IK,ICaLCa,ICaLK,ICaLNa,ICaLCads,ICaLKds,
     '  ICaLNads,IbCa,IKACH,INaK,INaCa,INaCads,Iup,
     '  Itr,Irel,XSRrel,Ito,I(26)

      RTF = OH5_R*OH5_T/OH5_F

! State variables
C ** DPN 04/11/98 - add f2Ca, slow ICaL inactivation
C ** DPN 05/11/98 - add Ito - r & q variables
      V          = y( 1) !membrane potential
      mNa        = y( 7) !Na channel m gate
      hNa        = y( 6) !Na channel h gate
      dCa        = y( 4) !Ca channel d gate
      fCa        = y( 5) !Ca channel f gate
      f2Ca       = y(11) !Ca channel slow inactivation gate
      xr1        = y(20) !fast K channel, fast rate
      xr2        = y(21) !fast K channel, slow rate
      xs         = y(22) !slow K channel
      xACh       = y(23) !??
      rto        = y(12) !inactivation variable for ito
      qto        = y(13) !activation variable for ito
      Oh5_Nai        = y( 2) !intracellular [Na] 
      Oh5_Ki         = y( 8) !intracellular [K]
      Oh5_Cao        = y(16) !extracellular [Ca]
      Oh5_Cai        = y( 3) !intracellular [Ca]
      Cads       = y(19) !sub-space [Ca]
c      Caup       = y( 9) ![Ca]NSR
c      Carel      = y(10) ![Ca]JSR
      CaCALMOD   = y(14) ![Ca] bound to calmodulin
      CaTROP     = y(15) ![Ca] bound to troponin
      factivator = y(17) !??? Ca SR release rate ???
      fproduct   = y(18) !???

! Equilibrium potentials
c      ENa  = RTF*DLOG(OH5_Nao/Oh5_Nai)
                        !Na channel
      IF(OH5_Nao/Oh5_Nai.GT.0.0d0) THEN
        ENa  = RTF*DLOG(OH5_Nao/Oh5_Nai)                 !Na channel
        POTS(1) = ENa
      ELSE
        WRITE(*,*) ' >>Attempted to log a -ve (Na)'
      ENDIF
      IF(OH5_Ko/Oh5_Ki.GT.0.0d0) THEN
        EK   = RTF*DLOG(OH5_Ko/Oh5_Ki)
        POTS(2) = EK
      ELSE
        WRITE(*,*) ' >>Attempted to log a -ve (K)'
      ENDIF 
                         !K  channel
c      EKs  = RTF*DLOG((OH5_Ko+OH5_PkNa*OH5_Nao)/(Oh5_Ki+OH5_PkNa*Oh5_Nai)) 
      IF((OH5_Ko+OH5_PkNa*OH5_Nao)/(Oh5_Ki+Oh5_Nai).GT.0.0d0) THEN
        EKs  = RTF*DLOG((OH5_Ko+OH5_PkNa*OH5_Nao)/(Oh5_Ki+Oh5_Nai)) 
        POTS(3) = EKs
      ELSE
        WRITE(*,*) ' >>Attempted to log a -ve (K2)'
      ENDIF

!Ks channel
c      ECa  = RTF*DLOG(Oh5_Cao/Oh5_Cai)                           !Ca
      IF(Oh5_Cai/Oh5_Cao.GT.0.0d0) THEN
        ECa  = -0.5d0*RTF*DLOG(Oh5_Cai/Oh5_Cao)                     !Ca
        POTS(4) = ECa
      ELSE
        WRITE(*,*) ' >>Attempted to log a -ve (Ca)'
        WRITE(*,'(''   [Ca]i = '',F12.6,'' [Ca]o = '',F12.6)') Oh5_Cai,
     '    Oh5_Cao
      ENDIF
      IF((OH5_Nao+0.12d0*OH5_Ko)/(Oh5_Nai+0.12d0*Oh5_Ki).GT.0.0d0) THEN
        Emh  = RTF*DLOG((OH5_Nao+0.12d0*OH5_Ko)/(Oh5_Nai+0.12d0*Oh5_Ki))
        POTS(5) = Emh
      ELSE
        WRITE(*,*) ' >>Attempted to log a -ve (NaK)'
      ENDIF

      CALL NOBLE98_CHANGE(t)

      CALL NOBLE98_RATES(V,Oh5_Cai,ALPHA,BETA)
      alfa_mNa = ALPHA(1)
      beta_mNa = BETA(1)
      alfa_hNa = ALPHA(2)
      beta_hNa = BETA(2)
      alfa_dCa = ALPHA(3)
      beta_dCa = BETA(3)
      alfa_fCa = ALPHA(4)
      beta_fCa = BETA(4)
      alfa_xr1 = ALPHA(5)
      beta_xr1 = BETA(5)
      alfa_xr2 = ALPHA(6)
      beta_xr2 = BETA(6)
      alfa_xs = ALPHA(7)
      beta_xs = BETA(7)
      alfa_xACH = ALPHA(8)
      beta_xACH = BETA(8)
      alfa_rto = ALPHA(9)
      beta_rto = BETA(9)

      CALL NOBLE98_CURRENTS(Y,RTF,POTS,I)  ! calculate ionic currents
      INa = I(1)
      IpNa = I(2)
      IbNa = I(3)
      IK1 = I(4)
c      IKr = I(5)
c      IKs = I(6)
      IbK = I(7)
      IKATP = I(8)
      IKNa = I(9)
      IK = I(10)
      ICaLCa = I(11)
      ICaLK = I(12)
      ICaLNa = I(13)
      ICaLCads = I(14)
      ICaLKds = I(15)
      ICaLNads = I(16)
      IbCa = I(17)
      IKACH = I(18)
      INaK = I(19)
      INaCa = I(20)
      INaCads = I(21)
      Iup = I(22)
      Itr = I(23)
      Irel = I(24)
      XSRrel = I(25)
      Ito = I(26)

      IF (SINGLE_CELL) THEN
        !need this to get currents back to UNEMAP for plotting
        !"fast fix, temp hack" - Chris Bradley, 02 December 1998 !!!!!
        INa98 = I(1)
        IpNa98 = I(2)
        IbNa98 = I(3)
        IK198 = I(4)
        IKr98 = I(5)
        IKs98 = I(6)
        IbK98 = I(7)
        IKATP98 = I(8)
        IKNa98 = I(9)
        IK98 = I(10)
        ICaLCa98 = I(11)
        ICaLK98 = I(12)
        ICaLNa98 = I(13)
        ICaLCads98 = I(14)
        ICaLKds98 = I(15)
        ICaLNads98 = I(16)
        IbCa98 = I(17)
        IKACH98 = I(18)
        INaK98 = I(19)
        INaCa98 = I(20)
        INaCads98 = I(21)
        Iup98 = I(22)
        Itr98 = I(23)
        Irel98 = I(24)
        XSRrel98 = I(25)
        Ito98 = I(26)
      ENDIF

C      ISWTCH(1) = 0.0d0 !INa
C      ISWTCH(2) = 1.0d0 !IbNa
C      ISWTCH(3) = 1.0d0 !IpNa
C      ISWTCH(4) = 1.0d0 !INaK
C      ISWTCH(5) = 0.0d0 !INaCa
C      ISWTCH(6) = 1.0d0 !INaCads
C      ISWTCH(7) = 1.0d0 !IK1
C      ISWTCH(8) = 1.0d0 !IK
C      ISWTCH(9) = 1.0d0 !IKATP
C      ISWTCH(10) = 1.0d0 !IKACh
C      ISWTCH(11) = 1.0d0 !IbK
C      ISWTCH(12) = 1.0d0 !ICaLNa
C      ISWTCH(13) = 1.0d0 !ICaLK
C      ISWTCH(14) = 1.0d0 !ICaLCa
C      ISWTCH(15) = 1.0d0 !ICaLNads
C      ISWTCH(16) = 1.0d0 !ICaLKds
C      ISWTCH(17) = 1.0d0 !ICaLCads
C      ISWTCH(18) = 1.0d0 !IbCa

!d(V)/dt     
      F( 1) = -(ISWTCH(16)*INa+ISWTCH(4)*IbNa+ISWTCH(23)*IpNa+
     '  ISWTCH(18)*INaK+ISWTCH(17)*INaCa+ISWTCH(19)*INaCads !eqn 1
     '  +ISWTCH(14)*IK1+ISWTCH(13)*IK+ISWTCH(15)*IKATP+
     '  ISWTCH(1)*IKACh+ISWTCH(3)*IbK 
     '  +ISWTCH(10)*ICaLNa+ISWTCH(9)*ICaLK+ISWTCH(8)*ICaLCa+
     '  ISWTCH(7)*ICaLNads+ISWTCH(6)*ICaLKds+ISWTCH(5)*ICaLCads
     '  +ISWTCH(2)*IbCa+ISWTCH(20)*Ito+ISWTCH(21)*IKNa+
     '  OH5_Istim) / OH5_Cm
c      F( 1) = -(INa+IbNa+IpNa+INaK+INaCa+INaCads                !eqn 1
c     '  +IK1+IK+IKATP+IKACh+IbK 
c     '  +ICaLNa+ICaLK+ICaLCa+ICaLNads+ICaLKds+ICaLCads
c     '  +IbCa+OH5_Istim) / OH5_Cm
!d(mNa)/dt                                                          
      F( 7) = alfa_mNa *(1.d0-mNa ) - beta_mNa *mNa            !eqn 2
!d(hNa)/dt                                                         
      F( 6) = alfa_hNa *(1.d0-hNa ) - beta_hNa *hNa            !eqn 3 
!d(dCa)/dt                                                          
      F( 4) = alfa_dCa *(1.d0-dCa ) - beta_dCa *dCa            !eqn 4
!d(fCa)/dt                                                          
      F( 5) = alfa_fCa *(1.d0-fCa ) - beta_fCa *fCa            !eqn 5
C ** DPN 04/11/98 - add f2Ca, slow ICaL inactivation
!d(f2Ca)/dt - slow ICaL inactivation - from PASCAL code
      F(11) = OH5_Rdecay*(1.0d0-Cads/(Cads+OH5_Kdsoff)-f2Ca)
!d(xr1)/dt                                                          
      F(20) = alfa_xr1 *(1.d0-xr1 ) - beta_xr1 *xr1            !eqn 7
!d(xr2)/dt                                                          
      F(21) = alfa_xr2 *(1.d0-xr2 ) - beta_xr2 *xr2            !eqn 8
!d(xs)/dt
      F(22) = alfa_xs  *(1.d0-xs  ) - beta_xs  *xs             !eqn 9
!d(xACh)/dt
      F(23) = alfa_xACh*(1.d0-xACh) - beta_xACh*xACh           !eqn 10
C ** DPN 05/11/98 - add Ito - r & q variables
!d(rto)/dt
      F(12) = alfa_rto*(1.0d0-rto) - beta_rto*rto
c      F(12) = 0.0d0
!d(qto)/dt SPEED[19] = 1.0 ???????
      F(13) = 333.0d0*(1.0d0/(1.0d0+DEXP(-(V+4.0d0)/5.0d0))-qto)
c      F(13) = 0.0d0
!d(Nai)/dt
C ** DPN 04/11/98 - new from HEART code
C      F( 2) = -(INa+IpNa+IbNa+3.d0*(INaK+INaCa+      !eqn 11
C     '  INaCads)+ICaLNa)/(OH5_Vi*OH5_F)
      F( 2) = -(INa+IpNa+(IbNa*OH5_Nao/140.0d0)+
     '  3.d0*(INaK+INaCa)+ICaLNa)/(OH5_Vi*OH5_F)
!d(Ki)/dt
C ** DPN 04/11/98 - new from HEART code
C      F( 8) = -(IK1+ICaLK+IbK+IK-2.d0*INaK)          !eqn 12
C     '        /(OH5_Vi*OH5_F)
      F( 8) = -(IK1+ICaLK+IbK+IK-2.d0*INaK+ICaLKds
     '  +IKATP+IKNa+Ito)/(OH5_Vi*OH5_F)
!d(Cao)/dt
c      F(16) =  (ICaLCa+ICaLCads+IbCa-2.d0*INaCa)       !eqn 13
c     '        /(2.d0*OH5_Vcell*OH5_Vecs*OH5_F)
c     '        - OH5_DIFFCa*(Oh5_Cao-OH5_Cab)         
c      F(16) =  (ICaLCa+IbCa-2.d0*INaCa)       !eqn 13
c     '        /(2.d0*OH5_Vcell*OH5_Vecs*OH5_F)
c     '        - OH5_DIFFCa*(Oh5_Cao-OH5_Cab)         
      F(16) = 0.0d0
C ** DPN 05/11/98 - need F(14) and F(15) before can calculate F(3)
!d(CaCALMOD)/dt
      F(14) = 1.d5*(OH5_Mtrop-CaCALMOD)*Oh5_Cai - 50.d0*CaCALMOD   !eqn 18
!d(CaTROP)/dt
C ** DPN 04/11/98 - new from HEART code
C      F(13) = OH5_alfatrop*DEXP(OH5_gamatropSL)*Oh5_Cai*(OH5_Ctrop-CaTROP)
C     '  -OH5_betatrop*CaTROP                                   !eqn 19
      F(15) = OH5_alfatrop*Oh5_Cai*(OH5_Ctrop-CaTROP)
     '  -OH5_betatrop*CaTROP                                   !eqn 19
!d(Cai)/dt   
C ** DPN 04/11/98 - new from HEART code
C      F( 3) =  -(ICaLCa+IbCa-2.d0*INaCa)/(2.d0*OH5_Vi*OH5_F) !eqn 14
C     '        - Iup + Irel*OH5_VSRup*OH5_Vrel/(OH5_Vi*OH5_Vup)
C     '        - CaCALMOD - CaTROP 
C     '        + Cads/(Cads+OH5_kdecay)*OH5_Rdecay
c      F( 3) =  -(ICaLCa+IbCa-2.d0*INaCa)/(2.d0*OH5_Vi*OH5_F) !eqn 14
c     '        - Iup + Irel*OH5_VSRup*OH5_Vrel/(OH5_Vi*OH5_Vup)
c     '        - F(12) - F(13) + Cads*OH5_Kdsdecay*OH5_Vds
      F( 3) =  -(ICaLCa+IbCa-2.d0*INaCa)/(2.d0*OH5_Vi*OH5_F) !eqn 14
     '        - Iup + Irel*OH5_VSRup*OH5_Vrel/(OH5_Vi*OH5_Vup)
     '        - F(14) - F(15) + Cads*OH5_Kdsdecay*OH5_Vds/OH5_Vi
!d(Cads)/dt
C ** DPN 04/11/98 - new from HEART code
C      F(17) = -(ICaLCads-2.d0*INaCads)                     !eqn 15
C     '        /(2.d0*OH5_Vds*OH5_F)
C     '        - Cads/(Cads+OH5_kdecay)*OH5_Rdecay
c      F(17) = -(ICaLCads)/(2.d0*OH5_Vds*OH5_Vi*OH5_F)
c     '        - Cads*OH5_Kdsdecay
      F(19) = -(ICaLCads)/(2.d0*OH5_Vds*OH5_F)
     '        - Cads*OH5_Kdsdecay
!d(Caup)/dt
      F( 9) = OH5_Vi/OH5_VSRup*Iup - Itr                   !eqn 16
!d(Carel)/dt
      F(10) = OH5_Vup/OH5_Vrel*Itr - Irel                  !eqn 17
!d(factivator)dt
C ** DPN 04/11/98 - new from HEART code
C      F(15) = (1.0d0-factivator-fproduct)*(5.d2*XSRrel**2+6.d2* !eqn20
C     '  DEXP((V-4.d1)*8.d-2))-factivator*(5.d2*XSRrel**2+6.d1)
      F(17) = (1.0d0-factivator-fproduct)*(5.d2*XSRrel**2+0.0d0*
     '  DEXP((V-4.d1)*8.d-2))-factivator*(5.d2*XSRrel**2+6.d1)
!d(fproduct)/dt
      F(18) = factivator*(5.d2*XSRrel**2+6.d1)-fproduct      !eqn 21
      IF (V.LT.-50.0d0) THEN
        F(17) = F(17) * 5.0d0
        F(18) = F(18) * 5.0d0
      ENDIF
      RETURN
      END

      SUBROUTINE DEFINE_NOBLE98(PARAM)

C#### Subroutine: DEFINE_NOBLE98
C###  Description:
C###    Initialise Noble '98 variables/parameters

      IMPLICIT NONE
      !passed variables
      REAL*8 PARAM(*)
      !common blocks
      REAL*8    OH5_kcachoff,OH5_kmk1,OH5_kmK,
     '          OH5_kmNa,OH5_kNaCa,OH5_kmCa,
     '          OH5_k1,OH5_k2,OH5_k3,OH5_k4,
     '          OH5_kATP,OH5_KD,OH5_kkNa,
     '          OH5_kdsoff,OH5_kmCads,OH5_kdecay,
     '          OH5_Rdecay,OH5_INaKmax
      COMMON /OH5_PARAM1/
     '  OH5_kcachoff,OH5_kmk1,OH5_kmK,
     '  OH5_kmNa,OH5_kNaCa,OH5_kmCa,
     '  OH5_k1,OH5_k2,OH5_k3,OH5_k4,
     '  OH5_kATP,OH5_KD,OH5_kkNa,
     '  OH5_kdsoff,OH5_kmCads,OH5_kdecay,
     '  OH5_Rdecay,OH5_INaKmax
      REAL*8    OH5_GNa,OH5_GbK,OH5_GK1,OH5_GbNa,
     '          OH5_GpNa,OH5_GbCa,
     '          OH5_GKr1,OH5_GKr2,OH5_GKs,OH5_GKATP,OH5_GKACh,
     '          OH5_PCa,OH5_PCaK,OH5_PCaNa,OH5_PKNa
      COMMON /OH5_PARAM2/
     '  OH5_GNa,OH5_GbK,OH5_GK1,OH5_GbNa,
     '  OH5_GpNa,OH5_GbCa,
     '  OH5_GKr1,OH5_GKr2,OH5_GKs,OH5_GKATP,OH5_GKACh,
     '  OH5_PCa,OH5_PCaK,OH5_PCaNa,OH5_PKNa
      REAL*8    OH5_Cab,OH5_Ko,OH5_Nao,
     '          OH5_FractICaL,OH5_FractINaCa,
     '          OH5_DIFFCa,OH5_dNaCa,OH5_gama,OH5_nNaCa,
     '          OH5_kcyCa,OH5_kxcs,OH5_ksrCa,OH5_Cm,
     '          OH5_ALPHA12,OH5_BETA12
      COMMON /OH5_PARAM3/
     '  OH5_Cab,OH5_Ko,OH5_Nao,
     '  OH5_FractICaL,OH5_FractINaCa,
     '  OH5_DIFFCa,OH5_dNaCa,OH5_gama,OH5_nNaCa,
     '  OH5_kcyCa,OH5_kxcs,OH5_ksrCa,OH5_Cm,
     '  OH5_ALPHA12,OH5_BETA12
      REAL*8    OH5_F,OH5_R,OH5_T,
     '          OH5_Vecs,OH5_radius,OH5_length,
     '          OH5_Vup,OH5_Vrel,OH5_kmCa2,
     '          OH5_Mtrop,OH5_Ctrop,OH5_alfatrop,
     '          OH5_betatrop,OH5_gamatropSL,
     '          OH5_gamaSRSL,OH5_sacSL
      COMMON /OH5_PARAM4/
     '  OH5_F,OH5_R,OH5_T,
     '  OH5_Vecs,OH5_radius,OH5_length,
     '  OH5_Vup,OH5_Vrel,OH5_kmCa2,
     '  OH5_Mtrop,OH5_Ctrop,OH5_alfatrop,
     '  OH5_betatrop,OH5_gamatropSL,
     '  OH5_gamaSRSL,OH5_sacSL
      REAL*8    OH5_G_Ca_stretch,OH5_G_Na_stretch,OH5_G_K_stretch,
     '          OH5_G_Ns_stretch,OH5_G_An_stretch,
     '          OH5_E_Ns_stretch,OH5_E_An_stretch
      COMMON /OH5_PARAM5/
     '  OH5_G_Ca_stretch,OH5_G_Na_stretch,OH5_G_K_stretch,
     '  OH5_G_Ns_stretch,OH5_G_An_stretch,
     '  OH5_E_Ns_stretch,OH5_E_An_stretch
      REAL*8    OH5_alfaS,OH5_alfaV,
     '          OH5_alfaACh,OH5_betaACh,
     '          OH5_SLhst,OH5_IThst
      COMMON /OH5_PARAM6/
     '  OH5_alfaS,OH5_alfaV,
     '  OH5_alfaACh,OH5_betaACh,
     '  OH5_SLhst,OH5_IThst
      REAL*8 OH5_kdsdecay,OH5_Jsrleak
      COMMON /OH5_PARAM7/
     '  OH5_kdsdecay,OH5_Jsrleak
      REAL*8 OH5_Vcell,OH5_Vi,OH5_Vds,OH5_VSRup
      COMMON /OH5_VOLUMES/
     '  OH5_Vcell,OH5_Vi,OH5_Vds,OH5_VSRup
      REAL*8 OH5_STIMSIZE,OH5_FREQ
      COMMON /OH5_STIM/ 
     '  OH5_STIMSIZE,OH5_FREQ
      !local variables
      REAL*8 PI

      OH5_STIMSIZE   = PARAM( 1)
      OH5_Cm         = PARAM( 2)
      OH5_F          = PARAM( 3)
      OH5_R          = PARAM( 4)
      OH5_T          = PARAM( 5)
      OH5_Rdecay     = PARAM( 6)
      OH5_kdsoff     = PARAM( 7)
      OH5_Nao        = PARAM( 8)
      OH5_Ko         = PARAM( 9)
      OH5_Vecs       = PARAM(10)
      OH5_Vup        = PARAM(11)
      OH5_Vrel       = PARAM(12)
      OH5_kdsdecay   = PARAM(13)
      OH5_Mtrop      = PARAM(14)
      OH5_alfatrop   = PARAM(15)
      OH5_betatrop   = PARAM(16)
      OH5_Ctrop      = PARAM(17)
      OH5_gK1        = PARAM(18)
      OH5_kmK1       = PARAM(19)
      OH5_gKr1       = PARAM(20)
      OH5_gKr2       = PARAM(21)
      OH5_gKs        = PARAM(22)
      OH5_gNa        = PARAM(23)
      OH5_gpNa       = PARAM(24)
      OH5_gbNa       = PARAM(25)
      OH5_PCa        = PARAM(26)
      OH5_PCaK       = PARAM(27)
      OH5_PCaNa      = PARAM(28)
      OH5_gbCa       = PARAM(29)
      OH5_INaKmax    = PARAM(30)
      OH5_kmK        = PARAM(31)
      OH5_kmNa       = PARAM(32)
      OH5_FractINaCa = PARAM(33)
      OH5_kNaCa      = PARAM(34)
      OH5_gama       = PARAM(35)
      OH5_dNaCa      = PARAM(36)
      OH5_kcyCa      = PARAM(37)
      OH5_kxcs       = PARAM(38)
      OH5_ksrCa      = PARAM(39)
      OH5_kmCa2      = PARAM(40)
      OH5_Jsrleak    = PARAM(41)
      OH5_kmCa       = PARAM(42)
      OH5_kmCads     = PARAM(43)
      OH5_PKNa       = PARAM(44)
      OH5_FREQ       = PARAM(45)
      OH5_radius     = PARAM(46)
      OH5_length     = PARAM(47)
      OH5_ALPHA12    = PARAM(48)
      OH5_BETA12     = PARAM(49)

      OH5_FractICaL  = 1.d0   ! Fraction of ICaL -> sub-space

! Cell space volumes
      PI = 4.0d0*ATAN(1.0d0)
      OH5_Vcell = 1.d-9*PI*OH5_radius**2*OH5_length          !uL
      OH5_Vi    = (1.d0-OH5_Vecs-OH5_Vup-OH5_Vrel)*OH5_Vcell !uL
      OH5_Vds   = 0.1d0*OH5_Vi                               !uL
      OH5_VSRup = OH5_Vcell*OH5_Vup                          !uL

C      IF (DABS(PARAM(62)).LT.1.0d-3) THEN
C        OH5_STIMSIZE = -3.0d0
C      ELSE
C        OH5_STIMSIZE = PARAM(62)
C      ENDIF
C      OH5_FREQ = PARAM(63)
C
C      OH5_Cm       = 95.0d-6 !uF
C
C      OH5_kcachoff = 1.d-3   !mM
C      OH5_kmk1     = 10.d0   !mM
C      OH5_kmK      = 1.d0    !mM
C      OH5_kmNa     = 40.d0   !mM
C      OH5_kNaCa    = 5.0d-4 !1.d-4   !nA
C      OH5_kmCa     = 5.0d-4 !1.d-3   !mM
C      OH5_k1       = 1.2d4   !mM
C      OH5_k2       = 100.d0  !mM
C      OH5_k3       = 60.d0   !mM
C      OH5_k4       = 25.d0   !mM
C      OH5_kATP     = 0.1d0   !mM
C      OH5_KD       = 0.13d-3 !mM ACh !uM ACh
C      OH5_kkNa     = 20.d0   !mM
C      OH5_kdsoff   = 1.d-3   !mM
C      OH5_kmCads   = 1.d-2   !mM
C      OH5_kdecay   = 100.d-3 ! s ???? !ms
c      OH5_Rdecay   = 20.d0   !mM/s
C      OH5_INaKmax  = 0.7d0 !0.14d0  !nA
C      OH5_GNa      = 2.5d0  !0.5d0   !uS
C      OH5_GbK      = 0.0d0  !1.7d-3  !uS
C      OH5_GK1      = 0.5d0  !1.7d-2  !uS
C      OH5_GbNa     = 6.0d-4 !1.2d-4  !uS
C      OH5_GpNa     = 4.d-3   !uS
C      OH5_GbCa     = 2.5d-4 !5.d-5   !uS
C      OH5_GKr1     = 2.8d-3  !uS
C      OH5_GKr2     = 1.7d-3  !uS
C      OH5_GKs      = 3.2d-3  !uS
C      OH5_GKATP    = 0.0d0   !uS - value not given
C      OH5_GKACh    = 5.d-5   !uS
C      OH5_PCa      = 0.1 !5.d-2   !nA/mM
C      OH5_PCaK     = 2.d-3   !nA/mM
C      OH5_PCaNa    = 2.d-3   !nA/mM
C      OH5_PKNa     = 3.d-2   !nA/mM
C      OH5_Cab        = 2.d0  !mM
c      OH5_Ko         = 4.d0    !mM
c      OH5_Nao        = 140.d0  !mM
c      OH5_FractICaL  = 1.d0   ! Fraction of ICaL -> sub-space
c      OH5_FractINaCa = 1.d-3  !  "   "    " INaCa -> "  "   "
c      OH5_DIFFCa     = 5.d-4   !NA ?
c      OH5_dNaCa      = 0.0d0 !1.d-4   !
C      OH5_gama       = 0.5d0   !
C      OH5_nNaCa      = 3.d0    !
C      OH5_kcyCa      = 3.d-4   !mM
C      OH5_kxcs       = 0.4d0   !mM
C      OH5_ksrCa      = 0.5d0   !mM
c      OH5_F         = 96485.0d0  !Coulombs/mole
c      OH5_R         = 8314.41d0  !mJoules/mole/degK
c      OH5_T            = 310.d0  !degK
c      OH5_Vecs         = 0.4d0   !
c      OH5_radius       = 10.d0 !um
c      OH5_length       = 80.d0 !um
c      OH5_Vup          = 1.d-2   !
c      OH5_Vrel         = 1.d-1   !
c      OH5_kmCa2        = 250.d0  !nA/mM
c      OH5_Mtrop        = 0.02d0  !mM
C      OH5_Ctrop        = 0.15d0  !mM
C      OH5_alfatrop     = 1.0d5 ! mM/s   !5.d3    !/s
C      OH5_betatrop     = 2.d2    !/s
C      OH5_gamatropSL   = 2.5d0   !
C      OH5_gamaSRSL         = 1.5d0   !
C      OH5_sacSL        = 2.5d0   !
C      !OH5_GCa_stretch = 0.d0  !uS
C      !OH5_GNa_stretch = 0.d0   !uS
C      !OH5_GK_stretch  = 0.d0   !uS
C      !OH5_GNs_stretch = 0.d0   !uS
C      !OH5_GAn_stretch = 0.d0   !uS
C      !OH5_ENs_stretch = 0.d0   !mV
C      !OH5_EAn_stretch = 0.d0   !mV
C      OH5_alfaS     = 2.d0  !
C      OH5_alfaV     = 2.d0    !
C      OH5_alfaACh   = 0.5d0   !/s
C      OH5_betaACh   = 0.5d0   !/s
C      OH5_SLhst     = 2.d0    !um
C      OH5_IThst     = 1.d0    !au ?

      RETURN
      END

      SUBROUTINE NOBLE98_RATES(V,Cai,ALPHA,BETA)

C#### Subroutine: NOBLE98_RATES
C###  Description:
C###    Calculates the gating rate constants for the Noble '98
C###    rat ventricular myocyte model.

      IMPLICIT NONE
      !passed variables
      REAL*8 V,Cai,ALPHA(9),BETA(9)
      !common blocks
      REAL*8    OH5_kcachoff,OH5_kmk1,OH5_kmK,
     '          OH5_kmNa,OH5_kNaCa,OH5_kmCa,
     '          OH5_k1,OH5_k2,OH5_k3,OH5_k4,
     '          OH5_kATP,OH5_KD,OH5_kkNa,
     '          OH5_kdsoff,OH5_kmCads,OH5_kdecay,
     '          OH5_Rdecay,OH5_INaKmax
      COMMON /OH5_PARAM1/
     '  OH5_kcachoff,OH5_kmk1,OH5_kmK,
     '  OH5_kmNa,OH5_kNaCa,OH5_kmCa,
     '  OH5_k1,OH5_k2,OH5_k3,OH5_k4,
     '  OH5_kATP,OH5_KD,OH5_kkNa,
     '  OH5_kdsoff,OH5_kmCads,OH5_kdecay,
     '  OH5_Rdecay,OH5_INaKmax
      REAL*8    OH5_GNa,OH5_GbK,OH5_GK1,OH5_GbNa,
     '          OH5_GpNa,OH5_GbCa,
     '          OH5_GKr1,OH5_GKr2,OH5_GKs,OH5_GKATP,OH5_GKACh,
     '          OH5_PCa,OH5_PCaK,OH5_PCaNa,OH5_PKNa
      COMMON /OH5_PARAM2/
     '  OH5_GNa,OH5_GbK,OH5_GK1,OH5_GbNa,
     '  OH5_GpNa,OH5_GbCa,
     '  OH5_GKr1,OH5_GKr2,OH5_GKs,OH5_GKATP,OH5_GKACh,
     '  OH5_PCa,OH5_PCaK,OH5_PCaNa,OH5_PKNa
      REAL*8    OH5_Cab,OH5_Ko,OH5_Nao,
     '          OH5_FractICaL,OH5_FractINaCa,
     '          OH5_DIFFCa,OH5_dNaCa,OH5_gama,OH5_nNaCa,
     '          OH5_kcyCa,OH5_kxcs,OH5_ksrCa,OH5_Cm,
     '          OH5_ALPHA12,OH5_BETA12
      COMMON /OH5_PARAM3/
     '  OH5_Cab,OH5_Ko,OH5_Nao,
     '  OH5_FractICaL,OH5_FractINaCa,
     '  OH5_DIFFCa,OH5_dNaCa,OH5_gama,OH5_nNaCa,
     '  OH5_kcyCa,OH5_kxcs,OH5_ksrCa,OH5_Cm,
     '  OH5_ALPHA12,OH5_BETA12
      REAL*8    OH5_F,OH5_R,OH5_T,
     '          OH5_Vecs,OH5_radius,OH5_length,
     '          OH5_Vup,OH5_Vrel,OH5_kmCa2,
     '          OH5_Mtrop,OH5_Ctrop,OH5_alfatrop,
     '          OH5_betatrop,OH5_gamatropSL,
     '          OH5_gamaSRSL,OH5_sacSL
      COMMON /OH5_PARAM4/
     '  OH5_F,OH5_R,OH5_T,
     '  OH5_Vecs,OH5_radius,OH5_length,
     '  OH5_Vup,OH5_Vrel,OH5_kmCa2,
     '  OH5_Mtrop,OH5_Ctrop,OH5_alfatrop,
     '  OH5_betatrop,OH5_gamatropSL,
     '  OH5_gamaSRSL,OH5_sacSL
      REAL*8    OH5_G_Ca_stretch,OH5_G_Na_stretch,OH5_G_K_stretch,
     '          OH5_G_Ns_stretch,OH5_G_An_stretch,
     '          OH5_E_Ns_stretch,OH5_E_An_stretch
      COMMON /OH5_PARAM5/
     '  OH5_G_Ca_stretch,OH5_G_Na_stretch,OH5_G_K_stretch,
     '  OH5_G_Ns_stretch,OH5_G_An_stretch,
     '  OH5_E_Ns_stretch,OH5_E_An_stretch
      REAL*8    OH5_alfaS,OH5_alfaV,
     '          OH5_alfaACh,OH5_betaACh,
     '          OH5_SLhst,OH5_IThst
      COMMON /OH5_PARAM6/
     '  OH5_alfaS,OH5_alfaV,
     '  OH5_alfaACh,OH5_betaACh,
     '  OH5_SLhst,OH5_IThst
      !local variables
!     SMAR009 23/12/98      REAL*8 alfa_f
      REAL*8 alfa_mNa,beta_mNa,alfa_hNa,beta_hNa,alfa_dCa,beta_dCa,
     '  alfa_fCa,beta_fCa,alfa_xr1,beta_xr1,alfa_xr2,beta_xr2,
     '  alfa_xs,beta_xs,alfa_xACH,beta_xACH,alfa_rto,beta_rto

! Na channel gating rate constants
      alfa_mNa = 2.d2*(V+41.d0)/(1.d0-DEXP(-0.1d0*(V+41.d0)))
      beta_mNa = 8.d3*DEXP(-0.056d0*(V+66.d0))
      alfa_hNa = 2.d1*DEXP(-0.125d0*(V+75.d0))
      beta_hNa = 2.d3/(1.d0+3.2d2*DEXP(-0.1d0*(V+75.d0)))

! Ca channel gating rate constants
C ** DPN 04/11/98 - new from HEART code
C      alfa_dCa = 90.d0*(V+19.d0)/(1.d0-DEXP(-0.25d0*(V+19.d0))) !eqn 4
C      beta_dCa = 36.d0*(V+19.d0)/(DEXP(0.1d0*(V+19.d0))-1.d0)   !eqn 4
C      alfa_f   = 12.d0/(1.d0+DEXP(-0.25d0*(V+34.d0)))           !eqn 5
C      alfa_fCa = alfa_f*(1.19d2*Cai/(OH5_Kcachoff+Cai)+1.d0)    !eqn 5
C      beta_fCa = 6.25d0*(V+34.d0)/(DEXP(0.25d0*(V+34.d0))-1.d0) !eqn 5
      alfa_dCa = 30.d0*(V+19.d0)/(1.d0-DEXP(-0.25d0*(V+19.d0))) !eqn 4
c      alfa_dCa = 30.d0*(V+24.d0)/(1.d0-DEXP(-0.25d0*(V+24.d0))) !eqn 4
      beta_dCa = 12.d0*(V+19.d0)/(DEXP(0.1d0*(V+19.d0))-1.d0)   !eqn 4
c      beta_dCa = 12.d0*(V+24.d0)/(DEXP(0.1d0*(V+24.d0))-1.d0)   !eqn 4
      alfa_fCa = 6.25d0*(V+34.d0)/(DEXP(0.25d0*(V+34.d0))-1.d0) !eqn 5
      beta_fCa = 12.d0/(1.d0+DEXP(-0.25d0*(V+34.d0)))           !eqn 5

! K  channel gating rate constants
      alfa_xr1 = 50.d0/(1.d0 +DEXP(-(V-5.d0)/9.d0))              !eqn 7
      beta_xr1 = 0.05d0*DEXP(-(V-20.d0)/15.d0)                  !eqn 7
      alfa_xr2 = alfa_xr1                                       !eqn 8
      beta_xr2 = 0.4d0*DEXP(-((V+30.d0)/30.d0)**3)              !eqn 8
      alfa_xs  = 14.d0/(1.d0+DEXP(-(V-40.d0)/9.d0))             !eqn 9
      beta_xs  = DEXP(-V/45.d0)                                 !eqn 9
      alfa_xACh= OH5_alfaACh                                    !eqn 10
      beta_xACh= OH5_betaACh                                    !eqn 10

! transient outward channel
      alfa_rto = 0.033d0*DEXP(-V/17.0d0)
      beta_rto = 33.0d0/(1.0d0+DEXP(-(V+10)/8.0d0))

      ALPHA(1) = alfa_mNa
      BETA(1)  = beta_mNa
      ALPHA(2) = alfa_hNa
      BETA(2)  = beta_hNa
      ALPHA(3) = alfa_dCa
      BETA(3)  = beta_dCa
      ALPHA(4) = alfa_fCa
      BETA(4)  = beta_fCa
      ALPHA(5) = alfa_xr1
      BETA(5)  = beta_xr1
      ALPHA(6) = alfa_xr2
      BETA(6)  = beta_xr2
      ALPHA(7) = alfa_xs
      BETA(7)  = beta_xs
      ALPHA(8) = alfa_xACH
      BETA(8)  = beta_xACH
      ALPHA(9) = alfa_rto
      BETA(9)  = beta_rto

      RETURN
      END

      SUBROUTINE NOBLE98_CHANGE(TIME)

C#### Subroutine: NOBLE98_CHANGE
C###  Description:
C###    Sets the stimulation current for the Noble '98 model based on
C###    the current time.
C###
C###    ??? DPN - needs upgrading to include different pacing schemes

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:deoxs00.cmn'
      !passed variables
      REAL*8 TIME
      !common blocks
      REAL*8 OH5_Istim
      COMMON /OH5_STIMULUS/
     '  OH5_Istim
      REAL*8 OH5_STIMSIZE,OH5_FREQ
      COMMON /OH5_STIM/ 
     '  OH5_STIMSIZE,OH5_FREQ
      !local variables
      REAL*8 PERIOD !SMAR009 23/12/98 ,stimulation_time

C *** Set the stimulus curent
      IF (time.GE.TPS.AND.time.LE.TPS+TP) THEN
        OH5_Istim = OH5_STIMSIZE
      ELSE
        OH5_Istim = 0.0d0
      ENDIF

C *** Evaluate the new stimulus time
      IF(OH5_FREQ.GT.0.0d0) THEN
        IF (time.GT.TPS+TP) THEN
          PERIOD = 1.0d0 / OH5_FREQ !freq in Hz
          TPS = TPS + PERIOD
        ENDIF
      ENDIF

      RETURN
      END

      SUBROUTINE NOBLE98_CURRENTS(Y,RTF,POTS,I)

C#### Subroutine: NOBLE98_CURRENTS
C###  Description:
C###    Calculates the ionic currents for the Noble '98 model.

      IMPLICIT NONE
!     SMAR009  INCLUDE 'cmiss$reference:cell00.cmn'
      !passed variables
      REAL*8 Y(*),RTF,POTS(5),I(26)
      !common blocks
      REAL*8    OH5_kcachoff,OH5_kmk1,OH5_kmK,
     '          OH5_kmNa,OH5_kNaCa,OH5_kmCa,
     '          OH5_k1,OH5_k2,OH5_k3,OH5_k4,
     '          OH5_kATP,OH5_KD,OH5_kkNa,
     '          OH5_kdsoff,OH5_kmCads,OH5_kdecay,
     '          OH5_Rdecay,OH5_INaKmax
      COMMON /OH5_PARAM1/
     '  OH5_kcachoff,OH5_kmk1,OH5_kmK,
     '  OH5_kmNa,OH5_kNaCa,OH5_kmCa,
     '  OH5_k1,OH5_k2,OH5_k3,OH5_k4,
     '  OH5_kATP,OH5_KD,OH5_kkNa,
     '  OH5_kdsoff,OH5_kmCads,OH5_kdecay,
     '  OH5_Rdecay,OH5_INaKmax
      REAL*8    OH5_GNa,OH5_GbK,OH5_GK1,OH5_GbNa,
     '          OH5_GpNa,OH5_GbCa,
     '          OH5_GKr1,OH5_GKr2,OH5_GKs,OH5_GKATP,OH5_GKACh,
     '          OH5_PCa,OH5_PCaK,OH5_PCaNa,OH5_PKNa
      COMMON /OH5_PARAM2/
     '  OH5_GNa,OH5_GbK,OH5_GK1,OH5_GbNa,
     '  OH5_GpNa,OH5_GbCa,
     '  OH5_GKr1,OH5_GKr2,OH5_GKs,OH5_GKATP,OH5_GKACh,
     '  OH5_PCa,OH5_PCaK,OH5_PCaNa,OH5_PKNa
      REAL*8    OH5_Cab,OH5_Ko,OH5_Nao,
     '          OH5_FractICaL,OH5_FractINaCa,
     '          OH5_DIFFCa,OH5_dNaCa,OH5_gama,OH5_nNaCa,
     '          OH5_kcyCa,OH5_kxcs,OH5_ksrCa,OH5_Cm,
     '          OH5_ALPHA12,OH5_BETA12
      COMMON /OH5_PARAM3/
     '  OH5_Cab,OH5_Ko,OH5_Nao,
     '  OH5_FractICaL,OH5_FractINaCa,
     '  OH5_DIFFCa,OH5_dNaCa,OH5_gama,OH5_nNaCa,
     '  OH5_kcyCa,OH5_kxcs,OH5_ksrCa,OH5_Cm,
     '  OH5_ALPHA12,OH5_BETA12
      REAL*8    OH5_F,OH5_R,OH5_T,
     '          OH5_Vecs,OH5_radius,OH5_length,
     '          OH5_Vup,OH5_Vrel,OH5_kmCa2,
     '          OH5_Mtrop,OH5_Ctrop,OH5_alfatrop,
     '          OH5_betatrop,OH5_gamatropSL,
     '          OH5_gamaSRSL,OH5_sacSL
      COMMON /OH5_PARAM4/
     '  OH5_F,OH5_R,OH5_T,
     '  OH5_Vecs,OH5_radius,OH5_length,
     '  OH5_Vup,OH5_Vrel,OH5_kmCa2,
     '  OH5_Mtrop,OH5_Ctrop,OH5_alfatrop,
     '  OH5_betatrop,OH5_gamatropSL,
     '  OH5_gamaSRSL,OH5_sacSL
      REAL*8    OH5_G_Ca_stretch,OH5_G_Na_stretch,OH5_G_K_stretch,
     '          OH5_G_Ns_stretch,OH5_G_An_stretch,
     '          OH5_E_Ns_stretch,OH5_E_An_stretch
      COMMON /OH5_PARAM5/
     '  OH5_G_Ca_stretch,OH5_G_Na_stretch,OH5_G_K_stretch,
     '  OH5_G_Ns_stretch,OH5_G_An_stretch,
     '  OH5_E_Ns_stretch,OH5_E_An_stretch
      REAL*8    OH5_alfaS,OH5_alfaV,
     '          OH5_alfaACh,OH5_betaACh,
     '          OH5_SLhst,OH5_IThst
      COMMON /OH5_PARAM6/
     '  OH5_alfaS,OH5_alfaV,
     '  OH5_alfaACh,OH5_betaACh,
     '  OH5_SLhst,OH5_IThst
      REAL*8 OH5_kdsdecay,OH5_Jsrleak
      COMMON /OH5_PARAM7/
     '  OH5_kdsdecay,OH5_Jsrleak
      !local variables
      REAL*8 VRTF,c1,c2,z1,z2 !SMAR009 21/12/98 ,z3,a1,b1,a2,a3,b3,a4,b4,b5,b6,
      REAL*8 V,mNa,hNa,dCa,fCa,f2Ca,xr1,xr2,xs,xACH,Nai,Ki,Cao,Cai,Cads,
     '  Caup,Carel,factivator,rto,qto !SMAR009 21/12/98 ,CaCALMOD,CaTROP,fproduct
      REAL*8 ATP,OH5_GKNa,ACh
      REAL*8 ENa,EK,EKs,ECa,Emh
      REAL*8 INa,IpNa,IbNa,IK1,IKr,IKs,IbK,IKATP,
     '  IKNa,IK,ICaLCa,ICaLK,ICaLNa,ICaLCads,ICaLKds,
     '  ICaLNads,IbCa,IKACH,INaK,INaCa,INaCads,Iup,
     '  Itr,Irel,XSRrel,Ito

      ENa = POTS(1)
      EK = POTS(2)
      EKs = POTS(3)
      ECa = POTS(4)
      Emh = POTS(5)

      ATP = 5.0d0
      OH5_GKNa = 0.0d0
      ACh = 0.0d0

! State variables
C ** DPN 04/11/98 - add f2Ca, slow ICaL inactivation
C ** DPN 05/11/98 - add Ito - r & q variables
      V          = y( 1) !membrane potential
      mNa        = y( 7) !Na channel m gate
      hNa        = y( 6) !Na channel h gate
      dCa        = y( 4) !Ca channel d gate
      fCa        = y( 5) !Ca channel f gate
      f2Ca       = y(11) !Ca channel slow inactivation gate
      xr1        = y(20) !fast K channel, fast rate
      xr2        = y(21) !fast K channel, slow rate
      xs         = y(22) !slow K channel
      xACh       = y(23) !??
      rto        = y(12) !inactivation variable for ito
      qto        = y(13) !activation variable for ito
      Nai        = y( 2) !intracellular [Na] 
      Ki         = y( 8) !intracellular [K]
      Cao        = y(16) !extracellular [Ca]
      Cai        = y( 3) !intracellular [Ca]
      Cads       = y(19) !sub-space [Ca]
      Caup       = y( 9) ![Ca]NSR
      Carel      = y(10) ![Ca]JSR
c      CaCALMOD   = y(14) ![Ca] bound to calmodulin
c      CaTROP     = y(15) ![Ca] bound to troponin
      factivator = y(17) !??? Ca SR release rate ???
c      fproduct   = y(18) !???

C      V          = y( 1) !membrane potential
C      mNa        = y( 7) !Na channel m gate
C      hNa        = y( 6) !Na channel h gate
C      dCa        = y( 4) !Ca channel d gate
C      fCa        = y( 5) !Ca channel f gate
C      xr1        = y(17) !fast K channel, fast rate
C      xr2        = y(18) !fast K channel, slow rate
C      xs         = y(19) !slow K channel
C      xACh       = y(20) !??
C      Nai        = y( 2) !intracellular [Na] 
C      Ki         = y( 8) !intracellular [K]
C      Cao        = y(13) !extracellular [Ca]
C      Cai        = y( 3) !intracellular [Ca]
C      Cads       = y(16) !sub-space [Ca]
C      Caup       = y( 9) ![Ca]NSR
C      Carel      = y(10) ![Ca]JSR
C      CaCALMOD   = y(11) ![Ca] bound to calmodulin
C      CaTROP     = y(12) ![Ca] bound to troponin
C      factivator = y(14) !??? Ca SR release rate ???
C      fproduct   = y(15) !???

C      V          = y( 1) !membrane potential
C      mNa        = y( 2) !Na channel m gate
C      hNa        = y( 3) !Na channel h gate
C      dCa        = y( 4) !Ca channel d gate
C      fCa        = y( 5) !Ca channel f gate
C      xr1        = y( 6) !fast K channel, fast rate
C      xr2        = y( 7) !fast K channel, slow rate
C      xs         = y( 8) !slow K channel
C      xACh       = y( 9) !??
C      Nai        = y(10) !intracellular [Na] 
C      Ki         = y(11) !intracellular [K]
C      Cao        = y(12) !intracellular [Ca]
C      Cai        = y(13) !extracellular [Ca]
C      Cads       = y(14) !sub-space [Ca]
C      Caup       = y(15) ![Ca]NSR
C      Carel      = y(16) ![Ca]JSR
C      CaCALMOD   = y(17) ![Ca] bound to calmodulin
C      CaTROP     = y(18) ![Ca] bound to troponin
C      factivator = y(19) !??? Ca SR release rate ???
C      fproduct   = y(20) !???

! Na currents
      INa   = OH5_GNa*mNa*mNa*mNa*hNa*(V-Emh)                  !eqn 36
C ** DPN 04/11/98 - new from HEART code
C      IpNa98  = OH5_GpNa/(1.d0+DEXP(0.125d0*(V+52.d0)))*(V-ENa)  !eqn 37
      IpNa  = OH5_GpNa/(1.d0+DEXP(-0.125d0*(V+52.d0)))*(V-ENa)  !eqn 37
C ** DPN 04/11/98 - new from HEART code
C      IbNa  = OH5_GbNa*Nai/140.d0*(V-ENa)                      !eqn 38
      IbNa  = OH5_GbNa*(V-ENa)

! K currents
c      IK1   = OH5_GK1*OH5_Ko/(OH5_Ko+OH5_kmK1)*(V-EK)          !eqn 29
c     '        /(1.d0+DEXP((V-EK-10.d0)/(RTF/2.0d0)))
c      IK1   = OH5_GK1*OH5_Ko/(OH5_Ko+OH5_kmK1)*(V-EK)          !eqn 29
c     '        /(1.d0+DEXP((V-EK+10.d0)/(RTF/2.0d0)))
      IK1   = OH5_GK1*OH5_Ko/(OH5_Ko+OH5_kmK1)*(V-EK)          !eqn 29
     '        /(1.d0+DEXP(1.25d0*(V-EK-10.d0)/RTF))
      IKr   = (OH5_GKr1*xr1+OH5_GKr2*xr2)*(V-EK)               !eqn 31
     '        /(1.d0+DEXP((V+9.d0)/22.4d0))         
      IKs   = OH5_GKs*xs*xs*(V-EKs)                            !eqn 32
      IbK   = OH5_GbK*(V-EK)                                   !eqn 33
      IKATP = OH5_GKATP*(V+80.d0)/(1.d0+(ATP/OH5_kATP)**2)     !eqn 34
      IKNa  = OH5_GKNa*(V-EK)*Nai/(Nai+OH5_kkNa)               !eqn 35
c      IK    = IKr+IKs+IKNa                               !eqn 30 
      IK    = IKr+IKs

! Transient outward current
      Ito = 0.005d0*rto*qto*Cai/(Cai+0.0005d0)*(V-EK)
c      Ito = 0.0d0
      
! Ca currents
C ** DPN 04/11/98 - new from HEART code
C      a1 = OH5_PCa*dCa*(1.d0-fCa)
C     '    *OH5_kcachoff/(OH5_kcachoff+Cai)
C      b1 = OH5_PCa*dCa*(1.d0-fCa)
C     '    *OH5_kdsoff/(OH5_kdsoff+Cads)
C      a2 = (V-50.d0)/RTF
C      a3 = DEXP(-a2)
C      b3 = DEXP(-a2/2.d0)
C      a4 = a2/(1.d0-a3)
C      b4 = a2/(1.d0-b3)
C      b5 = DEXP(1.d2/RTF)
C      b6 = DEXP(50.0d0/RTF)
C
C      ICaLCa   = 4.d0     *a1*b4*(Cai *b5-Cao*b3)             !eqn 39
C      ICaLK    = OH5_PCaK *a1*a4*(Ki  *b6-OH5_Ko *a3)         !eqn 40
Cc      ICaLK    = 1.0d0/oh5_pca *a1*b4*(Ki  *b6-OH5_Ko *a3)         !eqn 40
C      ICaLNa   = OH5_PCaNa*a1*a4*(Nai *b6-OH5_Nao*a3)         !eqn 41
Ccc      ICaLCads = 4.d0     *b1*b4*(Cads*b5-Cao*b3)             !eqn 42
Cc      ICaLKds  = OH5_PCaK *b1*a4*(Ki  *b6-OH5_Ko *a3)         !eqn 43
Ccc      ICaLKds  = OH5_PCaK *b1*b4*(Ki  *b6-OH5_Ko *a3)         !eqn 43
Ccc      ICaLNads = OH5_PCaNa*b1*a4*(Nai *b6-OH5_Nao*a3)         !eqn 44
C      
C      z3 = OH5_kcachoff/(OH5_kcachoff+Cai)
C
C      z1 = b4
C      ICaLK = dCa*(1.0d0-fCa)*z1*(Ki*b6-OH5_Ko*a3)*z3
C
C      z1 = b4
C      z2 = z1*(Cai*b5-Cao*b3)
C      ICaLCa = 4.0d0*OH5_PCa*dCa*(1.0d0-fCa)*z2*z3
C
C      z1 = OH5_PCaNa*OH5_PCa*a4
C      ICaLNa = dCa*(1.0d0-fCa)*z1*(Nai*b6-OH5_Nao*a3)*z3
      
      z1 = OH5_PCaK*OH5_PCa*(V-50.0d0)
     '  /(RTF*(1.0d0-DEXP(-(V-50.0d0)/RTF)))
      ICaLK = dCa*fCa*z1*(Ki*DEXP(50.0d0/RTF)-
     '  OH5_Ko*DEXP(-(V-50.0d0)/RTF))
      
      z1 = OH5_PCaNa*OH5_PCa*(V-50.0d0)
     '  /(RTF*(1.0d0-DEXP(-(V-50.0d0)/RTF)))
      ICaLNa = dCa*fCa*z1*(Nai*DEXP(50.0d0/RTF)-
     '  OH5_Nao*DEXP(-(V-50.0d0)/RTF))

      z1 = (V-50.0d0)/(RTF*(1.0d0-DEXP(-(V-50.0d0)/(0.5d0*RTF))))
      z2 = z1*(Cai*DEXP(50.0d0/(0.5d0*RTF))-
     '  Cao*DEXP(-(V-50.0d0)/(0.5d0*RTF)))
c      ICaLCa = 4.0d0*OH5_PCa*dCa*fCa*z2
c      z2 = z1*(Cads*DEXP(50.0d0/(0.5d0*RTF))-
c     '  Cao*DEXP(-(V-50.0d0)/(0.5d0*RTF)))
      ICaLCa = 4.0d0*OH5_PCa*dCa*fCa*z2

      ICaLCa = ICaLCa*f2Ca
      ICaLK  = ICaLK *f2Ca
      ICaLNa = ICaLNa*f2Ca
c      z1 = OH5_kcachoff/(OH5_kcachoff+Cai)
c      ICaLCa = ICaLCa*z1
c      ICaLK  = ICaLK *z1
c      ICaLNa = ICaLNa*z1

      IF (OH5_FractICaL.EQ.1.0d0) THEN
        ICaLCads = ICaLCa
        ICaLNads = ICaLNa
        ICaLKds  = ICaLK
        ICaLCa = 0.0d0
        ICaLNa = 0.0d0
        ICaLK  = 0.0d0
      ELSE
        write(*,*) 'BOB BOB'
c        ICaLCads = OH5_FractICaL * ICaLCa
c        ICaLNads = OH5_FractICaL * ICaLNa
c        ICaLKds  = OH5_FractICaL * ICaLK
c        ICaLCa = 0.0d0
c        ICaLNa = 0.0d0
c        ICaLK  = 0.0d0
      ENDIF

      IbCa     = OH5_GbCa*(V-ECa)                             !eqn 45
C MLB divide by zero
C      IKACh    = OH5_GKACh*OH5_Ko/(OH5_Ko+OH5_kmk1)*xACh      !eqn 46
C     '  /(1.d0+(OH5_kD/ACh)**2)
C     '  *(V-EK)/(1.d0+DEXP((V-EK-1.d1)*2.d0/RTF))
      IKAch = 0.0d0

! Pumps and exchanger currents
      INaK     = OH5_INaKmax*OH5_Ko /(OH5_Ko +OH5_kmK)
     '                      *Nai/(Nai+OH5_kmNa)

      VRTF = V/RTF 
      c1   = DEXP(OH5_gama*VRTF)
      c2   = DEXP((OH5_gama-1.d0)*VRTF)
C ** DPN 04/11/98 - new from HEART code
C      INaCa =   OH5_kNaCa*(c1*Nai**3*Cao - c2*OH5_Nao**3*Cai )  !eqn 48
C     '  /(1.d0+OH5_dNaCa*(   OH5_Nao**3*Cai +    Nai**3*Cao ))
C      INaCads = OH5_kNaCa*(c1*Nai**3*Cao - c2*OH5_Nao**3*Cads)  !eqn 49
C     '  / (1.d0+OH5_dNaCa*(   OH5_Nao**3*Cads+    Nai**3*Cao ))

      z1 = 1.0d0+OH5_dNaCa*(Cai*OH5_Nao**3+Cao*Nai**3)
      z2 = Nai**3*Cao*c1-OH5_Nao**3*Cai*c2
      INaCa = (1.0d0-OH5_FractINaCa)*OH5_kNaCa*z2/z1
      INaCa = INaCa/(1.0d0+Cai/0.0069d0)

      z1 = 1.0d0+OH5_dNaCa*(Cads*OH5_Nao**3+Cao*Nai**3)
      z2 = Nai**3*Cao*c1-OH5_Nao**3*Cads*c2
      INaCads = OH5_FractINaCa*OH5_kNaCa*z2/z1
      INaCads = INaCads/(1.0d0+Cads/0.0069d0)

! Ca sequestration
C ** DPN 04/11/98 - new from HEART code
C      Iup    = 3.d0*Cai-0.23d0*Caup*OH5_kcyCa*OH5_kxcs/OH5_ksrCa !eqn 50
C     '            /(Cai+       Caup*OH5_kcyCa*OH5_kxcs/OH5_ksrCa
C     '                             +OH5_kcyCa*OH5_kxcs+OH5_kcyCa)
c      Iup    = (0.4d0*Cai-0.03d0*Caup*OH5_kcyCa*OH5_kxcs/OH5_ksrCa) !eqn 50
c     '            /(Cai+       Caup*OH5_kcyCa*OH5_kxcs/OH5_ksrCa
c     '                             +OH5_kcyCa*OH5_kxcs+OH5_kcyCa)
      Iup    = (OH5_ALPHA12*Cai-
     '  OH5_BETA12*Caup*OH5_kcyCa*OH5_kxcs/OH5_ksrCa)
     '  /(Cai+       Caup*OH5_kcyCa*OH5_kxcs/OH5_ksrCa
     '  +OH5_kcyCa*OH5_kxcs+OH5_kcyCa)

      Itr     = 5.d1*(Caup-Carel)                                !eqn 51
c      Irel    = (factivator/(factivator+0.25d0))**2              !eqn 52
c     '  *OH5_kmCa2*Carel+OH5_Jsrleak*DEXP(OH5_gamaSRSL)
      Irel    = (factivator/(factivator+0.25d0))**2              !eqn 52
     '  *OH5_kmCa2*Carel+OH5_Jsrleak*Carel
      XSRrel  = Cai/(Cai+OH5_kmCa)                               !eqn 53
     '  +(1.d0-Cai/(Cai+OH5_kmCa))*Cads/(Cads+OH5_kmCads)

      I( 1) = INa
      I( 2) = IpNa
      I( 3) = IbNa
      I( 4) = IK1
      I( 5) = IKr
      I( 6) = IKs
      I( 7) = IbK
      I( 8) = IKATP
      I( 9) = IKNa
      I(10) = IK
      I(11) = ICaLCa
      I(12) = ICaLK
      I(13) = ICaLNa
      I(14) = ICaLCads
      I(15) = ICaLKds
      I(16) = ICaLNads
      I(17) = IbCa
      I(18) = IKACH
      I(19) = INaK
      I(20) = INaCa
      I(21) = INaCads
      I(22) = Iup
      I(23) = Itr
      I(24) = Irel
      I(25) = XSRrel
      I(26) = Ito

      RETURN
      END


      REAL*8 FUNCTION FN_CA_CELL(PARAM,t)

C#### Function: FN_CA_CELL
C###  Type: REAL*8
C###  Description:
C###    Calcium twitch function for use with the HMT model
      
      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cell_hmt.inc'

!     Parameter List
      REAL*8 PARAM(*),t

      FN_CA_CELL=PARAM(Ca_max)*t/PARAM(Ca_tau)*
     '  DEXP(1.d0-t/PARAM(Ca_tau))

      RETURN
      END


      REAL*8 FUNCTION FN_TN_CELL(PARAM,Ca_conc,boundCa,To,T)

C#### Function: FN_TN_CELL
C###  Type: REAL*8
C###  Description:
C###    Troponin kinetics function for HMT model

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cell_hmt.inc'
      INCLUDE 'cmiss$reference:tol00.cmn'

!     Parameter List
      REAL*8 PARAM(*),Ca_conc,boundCa,To,T
!     Local Variables
      REAL*8 TTo

      IF (DABS(To).LT.ZERO_TOL) THEN
        TTo = 0.0d0
      ELSE
        TTo = T/(PARAM(gamma)*To)
      ENDIF

      FN_TN_CELL=PARAM(Rho0)*Ca_conc*(PARAM(Cab_max)-boundCa)-
     '  PARAM(Rho1)*(1.0d0-TTo)*boundCa

      RETURN
      END


      REAL*8 FUNCTION FN_TM_CELL(PARAM,boundCa,sites,length)

C#### Function: FN_TM_CELL
C###  Type: REAL*8
C###  Description:
C###    Tropomyosin kinetics for HMT model (dz/dt)

      IMPLICIT NONE

      INCLUDE 'cmiss$reference:cell_hmt.inc'

!     Parameter List
      REAL*8 PARAM(*),boundCa,sites,length
!     Local Variables 
      REAL*8 C50,Tm_n,Tm_p50,Cb_norm

      !length dependence for n
      Tm_n   =PARAM(Tm_n_0)*(1.d0+PARAM(beta1)*(length-1.d0)) 
      !length dependence for p50
      Tm_p50 =PARAM(Tm_p50_0)*(1.d0+PARAM(beta2)*
     '  (length-1.d0))
C ***      C50=10**(6.d0-Tm_p50) !uM
      C50=10**(3.d0-Tm_p50) !mM

C DPN 05/08/98 - scale dz/dt by normalised [Ca]b
      Cb_norm =boundCa/PARAM(Cab_max)

      FN_TM_CELL=PARAM(alfa0)*(((boundCa/C50)*Cb_norm)**Tm_n*
     '  (1.d0-sites) - sites)

      RETURN
      END


      REAL*8 FUNCTION FN_TO_CELL(PARAM,length,sites)

C#### Function: FnTo
C###  Type: REAL*8
C###  Description:
C###    Tension-length-pCa relation for HMT model

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cell_hmt.inc'

!     Parameter List
      REAL*8 PARAM(*),length,sites
         
      FN_TO_CELL=PARAM(Tref)*(1.d0+PARAM(beta0)*
     '  (length-1.d0))*sites

      RETURN
      END


      REAL*8 FUNCTION FN_Q_CELL(DT,LENGTH,PARAM,ALPHA,AA,PHI)

C#### Function: FN_Q_CELL
C###  Type: REAL*8
C###  Description:
C###    Returns updated Q & PHI for current timestep dt, from fading   
C###    memory model.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cell_hmt.inc'
!     Parameter List
      REAL*8 DT,LENGTH,PARAM(*),ALPHA(3),AA(3),PHI(3)
!     Local Variables
      REAL*8 dPHI,PREVIOUS_LENGTH,LENGTH_DIFF
      REAL*8 FM_INTEGRAND1_CELL

      DATA PREVIOUS_LENGTH/1.0d0/
      SAVE PREVIOUS_LENGTH

      LENGTH_DIFF = LENGTH-PREVIOUS_LENGTH

      ! 1st Fading Memory term
      dPHI=0.5d0*(FM_INTEGRAND1_CELL(1,DT,0.d0,LENGTH_DIFF,PARAM,
     '  ALPHA)+FM_INTEGRAND1_CELL(1,DT,DT,LENGTH_DIFF,PARAM,
     '  ALPHA))*DT
      !PARAM(PHI1)=DEXP(-ALPHA(1)*dt)*PARAM(PHI1)+dPHI
      PHI(1)=DEXP(-ALPHA(1)*dt)*PHI(1)+dPHI

      ! 2nd Fading Memory term
      dPHI=0.5d0*(FM_INTEGRAND1_CELL(2,DT,0.d0,LENGTH_DIFF,PARAM,
     '  ALPHA)+FM_INTEGRAND1_CELL(2,DT,DT,LENGTH_DIFF,PARAM,
     '  ALPHA))*DT
      !PARAM(PHI2)=EXP(-ALPHA(2)*dt)*PARAM(PHI2)+dPHI
      PHI(2)=EXP(-ALPHA(2)*dt)*PHI(2)+dPHI

      ! 3rd Fading Memory term
      dPHI=0.5d0*(FM_INTEGRAND1_CELL(3,DT,0.d0,LENGTH_DIFF,PARAM,
     '  ALPHA)+FM_INTEGRAND1_CELL(3,DT,DT,LENGTH_DIFF,PARAM,
     '  ALPHA))*dt
      !PARAM(PHI3)=EXP(-ALPHA(3)*dt)*PARAM(PHI3)+dPHI
      PHI(3)=EXP(-ALPHA(3)*dt)*PHI(3)+dPHI


c      FN_Q_CELL = AA(1)*PARAM(PHI1)+AA(2)*PARAM(PHI2)+
c     '  AA(3)*PARAM(PHI3)
      FN_Q_CELL = AA(1)*PHI(1)+AA(2)*PHI(2)+AA(3)*PHI(3)

      PREVIOUS_LENGTH = LENGTH

      RETURN
      END

