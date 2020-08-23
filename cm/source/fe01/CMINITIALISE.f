      SUBROUTINE CMINITIALISE(BATCH_MODE,PRT1,PRT2,
     '  USESOCK,ERR)

C#### Subroutine: CMINITIALISE
C###  Description:
C###    Sets up the initial global variables for CMISS.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi02.inc'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'cmis00.cmn'
C      INCLUDE 'cmiss$reference:cspi00.cmn'
      INCLUDE 'diag00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'gks000.cmn'
      INCLUDE 'gtstr00.cmn'
      INCLUDE 'head00.cmn'
      INCLUDE 'learn00.cmn'
      INCLUDE 'map000.cmn'
C$    INCLUDE 'mp00.cmn'
      INCLUDE 'ntsg00.cmn'
      INCLUDE 'time02.cmn'
c      INCLUDE 'cmiss$reference:tol00.cmn'
      INCLUDE 'trac00.cmn'
!     Parameter List
      INTEGER BATCH_MODE,PRT1,PRT2,USESOCK,ERR
!     Local Variables
      INTEGER no_col,nsub
      REAL REALTIME
c      REAL*8 DLAMCH
      CHARACTER ERROR*(MXCH)

C SEN 5/1/03 Initialise the wallclock timer, and the cpu timer.
C     (We have to initialise the CPU timer, since a bug under
C     AIX/mp results in an erronious initial CPU time).
C     Also initialise the delta timers.
      CALL CPU_TIMER(CPU_TOTAL,REALTIME)
      CALL REAL_TIMER(REAL_TOTAL,REALTIME)
      CALL DCPU_TIMER(REALTIME)
      CALL DREAL_TIMER(REALTIME)

C     Initialise any variables used by enter/exit
      TR01=.FALSE.
      TR02=.FALSE.
      TR03=.FALSE.
      TR04=.FALSE.

      FILE00='file'
      PATH00='./'

C GMH 18/2/97 Initialise BLANK to spaces (not NULL's)
      DO no_col=1,80
        BLANK(no_col:no_col)=' '
      ENDDO !no_col
      GKS=.FALSE.
C      MAPOPN=.FALSE.
      NUM_LIBRARY=0

      FIRST_SYNTAX=.TRUE.
      LEARN=.FALSE.
C KAT 22/2/00: Only IREC_COMFILE(0) is used before set
      IREC_COMFILE(0)=-1
C      DO com=0,20
C        IREC_COMFILE(com)=-1
C      ENDDO
C      NEST=0
C KAT 22/2/00: First comfile read goes into FIRST_COM_UNIT+1 to leave
C     room for `open'ed comfiles.
      COM_UNIT=FIRST_COM_UNIT
      CALL GET_SYSTEM(ERROR,*9999)
      IF(USESOCK.EQ.0) THEN
        USE_SOCKET=.FALSE.
      ELSE
        USE_SOCKET=.TRUE.
        PORT1=PRT1
        PORT2=PRT2
      ENDIF
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
        IOER=0 !error output (stderr could be 7 on some systems)
        IOH1=6 !help 1 output
        IOH2=6 !help 2 output
C KAT 1Feb00: This seems unnecessary
C        CALL OPENF(5,'TERM','SYS$INPUT' ,'UNKNOWN',' ',' ',132,
C     '    ERROR,*9999)
C        CALL OPENF(6,'TERM','SYS$OUTPUT','UNKNOWN',' ',' ',255,
C     '    ERROR,*9999)
      ENDIF
      IF(BATCH_MODE.NE.0) THEN
        IOIP=-1 !no input. this should cause an error if used
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

      NEWLINE=ACHAR(10)

C KAT 28/3/00: Level of command echoing
      ECHO_RAW_COM=.TRUE.
      ECHO_INTERP_COM=.FALSE.

C CPB 30/11/92 the IO# units below should (?) be redundant.
      IO1=1
      IO2=2
      IO3=IOIP
      IO4=IOOP
      IVDU=IOIP
      IMAP=0
      ECHO_OUTPUT=.FALSE.
      NUM_STACK=1                     !<-|
      DOP_STACK(NUM_STACK)=.FALSE.    !  |- used for diagnostics
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
C KAT 7Nov00: Not used
C      DO_EXAMPLE=.FALSE.
C      SELECT_EXAMPLE=.FALSE.
      HEADING=' '
C      DO no_col=1,9
C        NAME_COL(no_col)='UNDEFINED'
C      ENDDO

C     Initialise any platform-dependent variables
      CALL SETUPCMISS(ERROR,*9999)
C     Initialise the cm-cmgui link
C KAT 26/4/00: Moved to maincmloop as GETCMISSCOMMANDPARAMS has not yet
C     been called.
C      CALL CMGUI_LINK_INITIALISE(ERROR,*9999)
CC     Initialise any general variables

      RETURN
 9999 WRITE(*,'('' Fatal Error: '',A)') ERROR
      ERR=1
      RETURN
      END


