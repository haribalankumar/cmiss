C**** CMISS Module FE10: Routines which involve Vax-VMS run-time library

C     Subroutine DIALOG     controls command parsing environment
C     Subroutine ENTERS     traces entry to a subprogram
C     Subroutine SYNTAX     executes valid commands
C     Function CM_RANDOM_NUMBER returns random number
C     Function TIMER        calls RTL to return system times
C     Subroutine COPY_FILE  copies file between directories DUMMY
C     Subroutine GETSTR1    returns string typed by user from keyboard mapping
C                                                           DUMMY
C     Subroutine GETSTR2    returns interrupt key number DUMMY
C     Subroutine GETSTR3    returns key for DOCUM DUMMY
C     Subroutine POST_FILE  posts file to printer DUMMY
C     Subroutine PURGE_FILE purges file versions  DUMMY



      SUBROUTINE DIALOG

C**** The dialog interpreter
C**** Note: pi and PI are assigned as user-defined names
C**** Condition Handler References:
C****   ch2 Introduction to System Routines
C****   ch4 LIB$ Manual -An Overview of the VAX Condition Handling Facility
C****   ch9 VAX Fortran User Manual
C****   ch8 VAX C Run time library reference manual
C****   section 8.4 VAX Debugger Manual
C****   ch9 Guide to Programming Resources ***THIS IS THE BEST***
C****   ch10 Introduction to VMS System Services

      PARAMETER (MXCH=255,MXSG=2000)
      PARAMETER (MXCO=16)
      PARAMETER (MXCOQU=16)
      INTEGER ISEG(MXSG),NTCOQUD(8)
      REAL A_TRANS(4,4),ANGLE(3),FIXED_PT(3),
     '  FPT(3),FSHFT(3),FSCALE(3),FANGLE(3),
     '  SCALE(3),SHIFT(3),
     '  PROJ_REF_PT_NEW(3),VIEW_PLANE_NEW(3),VIEW_REF_PT_NEW(3),
     '  VIEW_UP_NEW(3),
     '  A_MAP(4,4),A_ORIENT(4,4),A_ORIENT_NEW(4,4),PROJ_REF_PT(3),
     '  VIEW_REF_PT(3),
     '  VIEW_PLANE(3),VIEWPORT(6),VIEW_UP(3),WINDOW(4),WINDOW_NEW(4),
     '  NPC_CLIP(6),VMATRIX(4,4)
      LOGICAL CHANGE_HANDLER,CONTINUE,CTRLC,CHANGE,DEFINE,DIAGNO,
     '  DOP,END,FATAL_HANDLER,FILEIP,FIRST,FIRSTF,
     '  FIRST_DRAW,FIRST_PLOT,FIRST_OXS,FIRST_ZOOM,GKS,GRID,MENU,
     '  PHIGS,QUALIFY,SHOW_HIDE,SMG_READ,UPDATE,UPDATE_ZOOM,
     '  DISPLAY_VALUATOR_81,DISPLAY_VALUATOR_82,DISPLAY_VALUATOR_83,
     '  DISPLAY_VALUATOR_84,DISPLAY_VALUATOR_85,DISPLAY_VALUATOR_86,
     '  DISPLAY_VALUATOR_87,DISPLAY_VALUATOR_88,DISPLAY_VALUATOR_89,
     '  REFINE_ACTIVE,TRANSFORM_ACTIVE
      CHARACTER CSEG(MXSG)*60,STRING*(MXCH),PR*(MXCH),ERROR*(MXCH),
     '  CFROMI*5,CFROMR*11,CIW*1,CLASS*8,CHAR1*1,CHAR3*3,CHAR5*5,
     '  CHOOSE*16,
     '  COD(MXCO)*30,COQUD(MXCO,MXCOQU)*30,
     '  FILE00*50,FILE_NAME*20,FORMAT*200,
     '  HEADING*80,NAME_COL(0:9)*40,
     '  OPTION(30)*15,
     '  OPTION1(16)*11,OPTION2(11)*11,OPTION3(25)*11,OPTION4(10)*11,
     '  OPTION5(10)*20,OPTION6(10)*20,
     '  Q*1,STATSTR*80,STRG*132,TEXT_STRING*80,PARAMETER_TYPE*20

      DATA NXCH/MXCH/
      DATA MXIN/2147483647/
      DATA RXRE/1.70141178E+38/
      DATA RNRE/2.93873588E-39/
      DATA PI/3.14159265359/
      DATA E/2.71828182846/
      PARAMETER (IOCMX=4, IOAMX=4, IOLMX=4, IOIMX=99, IORMX=16)
      CHARACTER CDATA(IOCMX)*500,CDEFLT(IOCMX)*500,CBLANK(IOCMX)*500
      CHARACTER ADATA(IOAMX),ADEFLT(IOAMX),AYES(IOAMX),ANO(IOAMX)
      LOGICAL   LDATA(IOLMX),LDEFLT(IOLMX),LTRUE(IOLMX),LFALSE(IOLMX)
      INTEGER   IDATA(IOIMX),IDEFLT(IOIMX),IZERO(IOIMX),IONE(IOIMX)
      REAL      RDATA(IORMX),RDEFLT(IORMX),RZERO(IORMX),RONE(IORMX)
      COMMON /CIO/ CDATA,CDEFLT,CBLANK
      COMMON /AIO/ ADATA,ADEFLT,AYES,ANO
      COMMON /LIO/ LDATA,LDEFLT,LTRUE,LFALSE
      COMMON /IIO/ IDATA,IDEFLT,IZERO,IONE
      COMMON /RIO/ RDATA,RDEFLT,RZERO,RONE
      LOGICAL         CALL_BASE,CALL_ELEM,CALL_LINE,CALL_NODE,CALL_FIEL
      COMMON /CALL00/ CALL_BASE,CALL_ELEM,CALL_LINE,CALL_NODE,CALL_FIEL
      LOGICAL         CALL_EQUA,CALL_INIT,CALL_MATE,CALL_SOLV
      COMMON /CALL01/ CALL_EQUA,CALL_INIT,CALL_MATE,CALL_SOLV
      LOGICAL         CALL_DATA,CALL_FIT, CALL_GROW,CALL_MOTI
      COMMON /CALL02/ CALL_DATA,CALL_FIT, CALL_GROW,CALL_MOTI
      COMMON /CBSY01/ MXIN,RXRE,RNRE,PI,E

      COMMON /CBDI00/ FIRST
      COMMON /CBDI01/ NXCH
      COMMON /CBDI02/ IOIP,IOOP,IOER,IOHE
      COMMON /CBDI10/ NXCO,NXCOQU
      COMMON /CBPR01/ PR
      COMMON /CBPR02/ LNPR
      COMMON /CBFE01/ FIRSTF

      COMMON /GKS/ GKS
      COMMON /GKS001/ NINDICES
      COMMON /GKSDRAW/ FIRST_DRAW
      COMMON /GKSPLOT/ FIRST_PLOT
      COMMON /PHIG00/ PHIGS
      COMMON /PHIG01/ ANGLE,FIXED_PT,SCALE,SHIFT
      COMMON /PHIG03/ VIEW_REF_PT,VIEW_PLANE,VIEW_UP
      COMMON /PHIG04/ BACK_PLANE_DIST,FRONT_PLANE_DIST,PROJ_REF_PT,
     '  VIEW_PLANE_DIST,VIEWPORT,WINDOW,NPC_CLIP
      COMMON /PHIG05/ VIEW_REF_PT_NEW,VIEW_PLANE_NEW,VIEW_UP_NEW
      COMMON /PHIG06/ BACK_PLANE_DIST_NEW,FRONT_PLANE_DIST_NEW,
     '  PROJ_REF_PT_NEW,VIEW_PLANE_DIST_NEW,WINDOW_NEW
      COMMON /VIEW00/ ISVIEW

      COMMON /OXS000/ FIRST_OXS
      COMMON /CBSG00/ NTSG
      COMMON /DIAG00/ DIAGNO
      COMMON /DIAL00/ SMG_READ
      COMMON /DIAL01/ IWINDOW,SHOW_HIDE
      COMMON /DIAL02/ FATAL_HANDLER,CHANGE_HANDLER
      COMMON /DIAL03/ OLD_HANDLER
      COMMON /DIAL04/ MENU
      COMMON /CBWK01/ IWKS(99),IWKDEF(0:4)
      COMMON /CBWK02/ XDISP,YDISP,DISP
      COMMON /CBWC00/ XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,DIAG
      COMMON /FILE00/ FILE00
      COMMON /GEN004/ NAME_COL

      COMMON /OBJE00/ NTOBJE
      COMMON /OBJE01/ NSOBJE(20,50)
      CHARACTER OBJECT_NAME(20)*20
      COMMON /OBJE02/ OBJECT_NAME

      COMMON /B01/ IVDU,IFILE,IO1,IO2,IO3,IO4,IOTYPE,DOP
      COMMON /B03/ NBT,NDT,NET(9),NFT,NJT,NLT,NOT,NPT(9),NQT,NRT,NYT,NZT
      COMMON /B05/ ITYP1(9),ITYP2(9),ITYP3(9)

      LOGICAL TR01,TR02,TR03,TR04
      COMMON /CBTR01/ IOTR,TR01,TR02,TR03,TR04
      COMMON /CBZM00/ IZOOM(4),XNDC(20,4,4)
      COMMON /MAP00/ NTMAP,IMAP
      COMMON /CTRLC/ CTRLC
      COMMON /HEAD00/ HEADING

      INTEGER BUFFER_COUNT,BUFFER_LENGTH(99)
      COMMON /GETSTR00/ BUFFER_COUNT,BUFFER_LENGTH
      COMMON /GETSTR02/ IREC_COMFILE
      COMMON /GETSTR03/ NTKEY
      LOGICAL MACRO_DEFINE,MACRO_EXECUTE
      CHARACTER MACRO_BUFFER(50,0:9)*132
      COMMON /GETSTR04/ NTMACRO(0:9)
      COMMON /GETSTR05/ MACRO_DEFINE,MACRO_EXECUTE
      COMMON /GETSTR06/ MACRO_BUFFER
      COMMON /GETSTR07/ MACRO_KEY

      LOGICAL MAPOPN
      COMMON /CSPI/ MPLUN,MAPOPN

      DATA FPT/0.0,0.0,0.0/,FSHFT/0.0,0.0,0.0/,FSCALE/1.0,1.0,1.0/,
     '  FANGLE/0.0,0.0,0.0/

      NXCO=MXCO
      NXCOQU=MXCOQU
      GKS=.FALSE.
      MAPOPN=.FALSE.
      SMG_READ=.FALSE.
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
      IOIP=5
      IOOP=6
      IOER=6
      IOHE=6
      IOTR=6
      IO1=1
      IO2=2
      IO3=IOIP
      IO4=IOOP
      IVDU=IOIP
      IFILE=1
      IMAP=0
      NTSG=0
      NTKEY=0
      NTMACRO(0)=0
      MACRO_DEFINE=.FALSE.
      MACRO_EXECUTE=.FALSE.
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
      NT_VIEW=0
      HEADING=' '
      DO NO_COL=1,9
        NAME_COL(NO_COL)='UNDEFINED'
      ENDDO
      DO I=1,99
        IWKS(I)=0
      ENDDO
      IWKDEF(0)=1
      IWKDEF(1)=1
      DO NOSG=1,MXSG
        ISEG(NOSG)=0
        CSEG(NOSG)=' '
      ENDDO
      XMIN=-1.
      XMAX= 1.
      YMIN=-1.
      YMAX= 1.
      ZMIN=-1.
      ZMAX= 1.
      DIAG=SQRT(12.0)
      NJT=2
      FILE00='FILE'
      MENU=.FALSE.
      COD(2)='PI'
      COD(3)=CFROMR(PI,'(E11.5)')
      CALL ASSIGN(STRING,1,3,COD,NTCOQUD,COQUD,ERROR,*9999)
      COD(2)='pi'
      COD(3)=CFROMR(PI,'(E11.5)')
      CALL ASSIGN(STRING,1,3,COD,NTCOQUD,COQUD,ERROR,*9999)

C *** Establish condition handler
      FATAL_HANDLER=.FALSE.
      CHANGE_HANDLER=.FALSE.

C *** Check for cmiss.com file
      STRG='read cmiss;com'
      CALL PARSE(0,ISEG,CSEG,STRG,END,ERROR,*9999)

      CONTINUE=.TRUE.
      DO WHILE(CONTINUE) !is main program loop

          CALL SETSTR(STRING,ERROR,*150)
          IF(SMG_READ) THEN
            CALL GETSTR1(STRING,PR(:LNPR),LNPR,END,ERROR,*150)
          ELSE
            CALL GETSTR(STRING(LNPR-2:),PR(:LNPR),END,ERROR,*150)
          ENDIF

        IF(END) THEN       !CTRL-Z has been used to quit
          CALL QUIT(END,ERROR,*150)
          CONTINUE=.FALSE.
          GOTO 200
        ENDIF

C ***   Parse string
        IF(SMG_READ)THEN
          STRING=STRING(3:)
          CALL PARSE(0,ISEG,CSEG,STRING,END,ERROR,*150)
        ELSE
          CALL PARSE(0,ISEG,CSEG,STRING,END,ERROR,*150)
        ENDIF

C ***   Check for changes to condition handler
        IF(CHANGE_HANDLER) THEN
          IF(FATAL_HANDLER) THEN
          ENDIF
        ENDIF

        IF(END) THEN       !The QUIT command has been used to quit
          CONTINUE=.FALSE.
        ENDIF
        GOTO 200

C ***   Handle error condition
 150    CALL STRING_TRIM(ERROR,IBEG,IEND)
        WRITE(IOOP,'(A)') ' '//ERROR(IBEG:IEND)//'>DIALOG'
        IF(CTRLC) THEN
          CTRLC=.FALSE.
C         if mode is not batch
          IF(STATSTR(1:1).NE.'B')THEN
          ENDIF
        ENDIF
        ERROR(1:)=' '

 200    CONTINUE
      ENDDO !end of main program loop
      CALL EXITS('DIALOG')

      RETURN
 9999 CALL ERRORS('DIALOG',ERROR)
      WRITE(IOER,'(132A)') ' Fatal error propagated to DIALOG:'
      CALL STRING_TRIM(ERROR,IBEG,IEND)
      WRITE(IOER,'(132A)') ' '//ERROR(IBEG:IEND)//'>DIALOG'
      CLOSE(UNIT=IOER)
      CALL EXITS('DIALOG')
      RETURN
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

      IMPLICIT REAL (A-H,O-Z)
      CHARACTER NAME*(*),CFROMI*10,COLV*10,SUBNAM*60,TRSB*10
      LOGICAL ALLSUB,CTRLC,DIAGNO,DOP,TR01,TR02,TR03,TR04
      CHARACTER SB*60
      COMMON /CBTR00/ TRSB
      COMMON /CBTR01/ IOTR,TR01,TR02,TR03,TR04
      COMMON /CBTR02/ NOLV,NOSB,NTLV,NTSB,NXLV,NXSB,TMST
      COMMON /CBTR03/ NOSBLV(1)
      COMMON /CBTR04/ NOSBSM(1)
      COMMON /CBTR05/ TMEL(1)
      COMMON /CBTR06/ TMELSM(1)
      COMMON /CBTR07/ TMEN(1)
      COMMON /CBTR08/ TMTLSM(1)
      COMMON /CBTR09/ SB(1)
      COMMON /CBDI02/ IOIP,IOOP,IOER,IOHE
      COMMON /DIAG00/ DIAGNO,ALLSUB
      COMMON /DIAG01/ SUBNAM
      COMMON /CTRLC/ CTRLC
      COMMON /B01/ IVDU,IFILE,IO1,IO2,IO3,IO4,IOTYPE,DOP

      IF(DIAGNO) THEN
        CALL STRING_TRIM(NAME,IBEG1,IEND1)
        WRITE(IOOP,'('' *** Enter '',A)') NAME(IBEG1:IEND1)
        IF(ALLSUB) THEN
          DOP=.TRUE.
        ELSE IF(.NOT.ALLSUB) THEN
          CALL STRING_TRIM(SUBNAM,IBEG2,IEND2)
          IF(NAME(IBEG1:IEND1).EQ.SUBNAM(IBEG2:IEND2)) DOP=.TRUE.
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
              WRITE(IOTR,'(/''      Time:    Calls:    Level:'//
     '          ' >Subprogram entered'''//
     '          '/''      Time:    Total:   Actual:'//
     '          ' <Subprogram exited'')')
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
            COLV=CFROMI(NOLV,'(I10)')
            CALL STRING_TRIM(COLV,IBEG,IEND)
            WRITE(IOTR,'(1X,F10.3,I10,I10,'//COLV(IBEG:IEND)//
     '        '('' >''),A)')
     '        TMEN(NOLV),NOSBSM(NOSB),NOLV,NAME
          ENDIF
        ENDIF
        TR01=.TRUE.
      ENDIF

      RETURN
 9999 RETURN 1
      END

      SUBROUTINE SYNTAX(ISEG,CSEG,STRING,NTCO,CO,NTCOQU,COQU,END,
     '  ERROR,*)

C**** Determines whether CO is a valid command
C**** If CO is known the command qualifiers are checked.
C**** When the qualifiers are understood the command is executed.

      PARAMETER (MXCH=255)
      INTEGER ISEG(*),NTCOQU(*),BUFFER_COUNT,BUFFER_LENGTH(99),
     '  IWK(6)
      REAL COLOUR(3)
      LOGICAL ABBREV,ALLSUB,CBBREV,CHANGE_HANDLER,
     '  DIAGNO,DOP,END,GKS,FATAL_HANDLER,LEARN,MENU,SMG_READ,
     '  UPVU,UPVUOP
      CHARACTER STRING*(*),CO(*)*(*),COQU(NXCO,*)*(*),CSEG(*)*(*),
     '  ERROR*(*),CUPPER*(MXCH),OPTION(30)*80,BLANK*80,BUFFER(99)*132,
     '  FILE*50,FILE00*50,MACRO_BUFFER(50,0:9)*132,NODE*5,SUBNAM*60
      COMMON /CBDI02/ IOIP,IOOP,IOER,IOHE
      COMMON /CBDI10/ NXCO,NXCOQU
      COMMON /GKS/ GKS,UPVU,UPVUOP
      COMMON /MAP00/ NTMAP,IMAP
      COMMON /LEARN/ LEARN
      COMMON /FILE00/ FILE00
      COMMON /DIAG00/ DIAGNO,ALLSUB
      COMMON /DIAG01/ SUBNAM
      COMMON /DIAL00/ SMG_READ
      COMMON /DIAL02/ FATAL_HANDLER,CHANGE_HANDLER
      COMMON /DIAL04/ MENU
      COMMON /GETSTR00/ BUFFER_COUNT,BUFFER_LENGTH
      COMMON /GETSTR01/ BUFFER
      COMMON /GETSTR04/ NTMACRO(0:9)
      COMMON /GETSTR06/ MACRO_BUFFER

      COMMON /B01/ IVDU,IFILE,IO1,IO2,IO3,IO4,IOTYPE,DOP
      COMMON /B03/ NBT,NDT,NET(9),NFT,NJT,NLT,NOT,NPT(9),NQT,NRT,NYT,NZT
      COMMON /K25/ KTYP25,KTYP26,KTYP27,KTYP28,KTYP29
      PARAMETER (TOL=1.E-4)

C *** Define the variables needed to open the MAP-4000
      LOGICAL INF     /.TRUE./      ! Display message mode for MBOPN on
      LOGICAL LOAD    /.TRUE./      ! Load map mode for MBOPN on
      LOGICAL GETALT  /.TRUE./      ! Get alternate MAP if map_wanted is busy
      INTEGER MPLUN                 ! is MAP logical unit number
      INTEGER WT      /0/           ! is wait interval for busy map
      INTEGER TIMEOUT /86400/       ! seconds to wait for timeout
      INTEGER LFREE,MFREE           ! are bytes available in local,main memory
      INTEGER LMSIZ,MMSIZ           ! are size of local,main memory in bytes
      LOGICAL MAPOPN,PUTPICS,PUTMAPIO,PUTMASKMAP
      COMMON /CSPI/ MPLUN,MAPOPN,PUTPICS,PUTMAPIO,PUTMASKMAP
      DATA BLANK/' '/

      CALL ENTERS('SYNTAX',*9999)
      STRING=' '
      DO NOCO=1,NTCO-1
        CALL STRING_TRIM(STRING,IBEG,IEND)
        STRING=STRING(1:IEND)//' '//CO(NOCO)
      ENDDO
      NOCO=1

      IF(ABBREV(CO(NOCO),'FEM',1)) THEN
        CALL FEM(ISEG,NOCO,NTCO,NTCOQU,CO,COQU,CSEG,END,STRING,
     '    ERROR,*9999)
      ELSE
        CALL STAND(ISEG,NOCO,NTCO,NTCOQU,CO,COQU,CSEG,END,STRING,
     '    ERROR,*9999)
      ENDIF

      CALL EXITS('SYNTAX')
      RETURN
 9999 CALL ERRORS('SYNTAX',ERROR)
      CALL EXITS('SYNTAX')
      RETURN 1
      END



      INTEGER FUNCTION TIMER(ICODE,IVALUE,IHANDLE)

C**** Calls RTL to return system times

      IMPLICIT NONE
!     Parameter List
      INTEGER ICODE,IHANDLE,IVALUE
      TIMER= 0
      RETURN
      END



      SUBROUTINE GETSTR1(STRING,PR,LNPR,ENDCOM,ERROR,*)

C**** Inserts a prompt, PROMPT, in the input file.
C**** Returns the string typed by the user
C**** If LEARN is true the string is written to FILE00
C**** The F10 key indicates an eof - END=100 in read sets ENDCOM=.TRUE.
C**** PROMPT is prompt for SMG routine
C**** Pressing the enter key causes a sequential read from unit 76 (com file)
C**** MACRO_BUFFER(nomacr,nokey),nomacr=1,NT_MACRO(nokey) are command lines
C**** associated with macro key number NOKEY.
C**** MACRO_DEFINE is .true. when macro is being defined ('Do' key to start
C**** and finish).
C**** MACRO_EXECUTE is .true. when macro is executing.
C**** NOTE: $ format is NONSTANDARD FORTRAN 77.

      IMPLICIT NONE
      INCLUDE 'cmiss$common:cbdi02.cmn'
      INCLUDE 'cmiss$common:gtstr00.cmn'
      INCLUDE 'cmiss$common:learn00.cmn'
      INCLUDE 'cmiss$common:trac00.cmn'
!     Parameter List
      INTEGER LNPR
      CHARACTER ERROR*(*),PR*(*),STRING*(*)
      LOGICAL ENDCOM

      RETURN

      END


      SUBROUTINE GETSTR2(ERROR,*)

C**** Returns interrupt key number from mapped keyboard.
C**** Called from MARCH1 in module FE07.

      IMPLICIT NONE
      INCLUDE 'cmiss$common:gtstr00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
      RETURN
      END


      SUBROUTINE GETSTR3(KEY,ERROR,*)

C**** Returns the key typed by the user for DOCUM window.
C**** Prev_Screen  =  1
C**** Next_Screen  =  2
C**** Do           =  3
C**** Select       =  4
C**** Insert       =  5
C**** Delete       =  6
C**** Remove       =  7
C**** Alphanumeric keys return their ASCII character code.

      IMPLICIT NONE
!     Parameter List
      INTEGER KEY
      CHARACTER ERROR*(*)
      RETURN
      END


      SUBROUTINE POST_FILE(FILE_NAME)

C**** Posts file to printer.

!     Parameter List
      CHARACTER FILE_NAME*(*)
      RETURN
      END


      SUBROUTINE PURGE_FILE(FILE_NAME)

C**** Purges file versions.

!     Parameter List
      CHARACTER FILE_NAME*(*)
      RETURN
      END


      REAL FUNCTION CM_RANDOM_NUMBER(ISEED)

C**** Returns real random number between 0.0 and 1.0

      IMPLICIT NONE
!     Parameter List
      INTEGER ISEED
!     Local Variables
      REAL RAN

      CALL RANDOM(RAN)
      CM_RANDOM_NUMBER=RAN

      RETURN
      END
