      SUBROUTINE SETUPCMISS(ERROR,*)

C#### Subroutine: SETUPCMISS
C###  Description:
C###    Setups up the initial platform dependent variables for CMISS.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
c      INCLUDE 'cmiss$reference:b01.cmn'
c      INCLUDE 'cmiss$reference:binf00.cmn'
c      INCLUDE 'cmiss$reference:cbfe01.cmn'
      INCLUDE 'cbdi02.cmn'
c      INCLUDE 'cmiss$reference:cbdi10.cmn'
c      INCLUDE 'cmiss$reference:cbwk01.cmn'
c      INCLUDE 'cmiss$reference:cmis00.cmn'
c      INCLUDE 'cmiss$reference:cspi00.cmn'
c      INCLUDE 'cmiss$reference:diag00.cmn'
      INCLUDE 'disp00.cmn'
c      INCLUDE 'cmiss$reference:docu00.cmn'
c      INCLUDE 'cmiss$reference:file00.cmn'
c      INCLUDE 'cmiss$reference:gen000.cmn'
c      INCLUDE 'cmiss$reference:geom00.cmn'
c      INCLUDE 'cmiss$reference:gks000.cmn'
c      INCLUDE 'cmiss$reference:graf00.cmn'
c      INCLUDE 'cmiss$reference:gtstr00.cmn'
c      INCLUDE 'cmiss$reference:head00.cmn'
      INCLUDE 'mach00.cmn'
      INCLUDE 'mach00.inc'
c      INCLUDE 'cmiss$reference:map000.cmn'
c      INCLUDE 'cmiss$reference:ntsg00.cmn'
c      INCLUDE 'cmiss$reference:tol00.cmn'
c      INCLUDE 'cmiss$reference:trac00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
c     SMAR009   INTEGER MXCO,MXCOQU,MXSG        !MXCH in 'mxch.inc'
c     22/12/98  PARAMETER (MXCO=25,MXCOQU=25,MXSG=10000)
      INTEGER I,J !,IBEG,IEND,NO_COL,NOSG,NSUB
c      CHARACTER STRING*(MXCH)

      CALL ENTERS('SETUPCMISS',*9999)

      DO I=1,100
        DO J=1,MXCH
          OP_STRING(I)(J:J)=' '
        ENDDO
      ENDDO
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
      ENDIANTYPE=CHAR(MACH_BIGENDIAN) !Big endian native
      CHARFORMTYPE=CHAR(MACH_CHARASCII) !ASCII character format
      INTFORMTYPE=CHAR(MACH_TWOSCOMP) !Two's complement format
      SPFORMTYPE=CHAR(MACH_SPIEEE) !IEEE single precision format
      DPFORMTYPE=CHAR(MACH_DPIEEE) !IEEE double precision format
      MACHTYPE=CHAR(MACH_SGI) !SGI
      OSTYPE=CHAR(MACH_UNIX) !UNIX

      EXAMPLES_DIR=' '

C      GKS=.FALSE.
C      MAPOPN=.FALSE.
C      ZERO_TOL=1.0E-10
C      CONVERG_TOL=DLAMCH('E')
C      LOOSE_TOL=DSQRT(CONVERG_TOL)
c
C      TR01=.FALSE.
C      TR02=.FALSE.
C      TR03=.FALSE.
C      TR04=.FALSE.
C      FIRST_SYNTAX=.TRUE.
C      FEM_ARRAYS_MALLOCED=.FALSE.
C      BUFFER_COUNT=0
C      IREC_COMFILE=0

C      CALL GET_REV_TIME(ERROR,*9999)
C      CALL GET_SYSTEM(ERROR,*9999)
C      CALL GET_DISPLAY(ERROR,*9999)


C CPB 11/10/95 Changing over unit number for terminal mode

C      IF(USE_SOCKET) THEN
C        IOIP=1 !input  file
C        IOOP=2 !usual output
C        IODI=3 !diagnostics output
C        IOTR=4 !trace output
C        IOER=5 !error output
C        IOH1=6 !help 1 output
C        IOH2=7 !help 2 output
C      ELSE
C        IOIP=5 !input  file
C        IOOP=6 !usual output
C        IODI=6 !diagnostics output
C        IOTR=6 !trace output
C        IOER=6 !error output
C        IOH1=6 !help 1 output
C        IOH2=6 !help 2 output
C        CALL OPENF(5,'TERM','SYS$INPUT' ,'UNKNOWN',' ',' ',132,
C     '    ERROR,*9999)
C        CALL OPENF(6,'TERM','SYS$OUTPUT','UNKNOWN',' ',' ',255,
C     '    ERROR,*9999)
C      ENDIF
C      IOOUT=9 !output FILE
C      IOFI=IOOP !default for listing to file commands
C      IFILE=10 !ifile is the file for input files i.e. ip* files
C      IOFILE1=11 !output file i.e. op* files
C      IOFILE2=12 ! iofile2 -> iofile6 are general purpose output files
C      IOFILE3=13
C      IOFILE4=14
C      IOFILE5=17
C      IOFILE6=18
C CPB 30/11/92 the IO# units below should (?) be redundant.
C      IO1=1
C      IO2=2
C      IO3=IOIP
C      IO4=IOOP
C      IVDU=IOIP
C      IMAP=0
C      DO NUM_STACK=1,5                !<-|
C        DOP_STACK(NUM_STACK)=.FALSE.  !  |
C      ENDDO                           !  |
C      NUM_STACK=5 !5 levels below FEM !  |- used for diagnostics
C      NT_SUB=0                        !  |
C      DO NSUB=1,MXSUB                 !  |
C        SUBNAM(NSUB)=' '              !  |
C      ENDDO                           !  |
C      DIAGNO=.FALSE.                  !  |
C      ALLSUB=.FALSE.                  !  |
C      FROMSUB=.FALSE.                 !  |
C      DOP=.FALSE.                     !<-|
C      NTSG=0
C      NT_KEY=0
C      NT_MACRO(0)=0
C      MACRO_DEFINE=.FALSE.
C      MACRO_COMMAND_EXECUTE=.FALSE.
C      MACRO_KEY_EXECUTE=.FALSE.
C      DO_EXAMPLE=.FALSE.
C      SELECT_EXAMPLE=.FALSE.
C      HEADING=' '
C      DO NO_COL=1,9
C        NAME_COL(NO_COL)='UNDEFINED'
C      ENDDO
C      DO I=1,99
C        IWKS(I)=0
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

      CALL EXITS('SETUPCMISS')
      RETURN
 9999 CALL ERRORS('SETUPCMISS',ERROR)
      CALL EXITS('SETUPCMISS')
      RETURN 1
      END

