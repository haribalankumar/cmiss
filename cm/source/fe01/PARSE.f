      SUBROUTINE PARSE(QUIT,LEN_COMMAND,COMMAND_LINE,*)

C#### Subroutine: PARSE
C###  Description:
C###    PARSE parses the string COMMAND_LINE to extract and execute command(s).

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'error0.inc'
      INCLUDE 'trac00.cmn'
      INCLUDE 'mxch.inc'
!     Parameter List
      LOGICAL QUIT
      INTEGER LEN_COMMAND
      CHARACTER COMMAND_LINE(LEN_COMMAND)
!     Local Variables
      CHARACTER ERROR*(ERRSTRLEN)
C KAT 28Jan00: This is not an appropriate place for ISEG and CSEG to be set
C     up.  SYNTAX is a better place, but I think they are only used for
C     graphics and so should be dynamically allocated.
      INTEGER MXSG
      PARAMETER (MXSG=10000)
      INTEGER i,IEND,ISEG(MXSG)
      CHARACTER CSEG(MXSG)*60
      DATA ISEG/MXSG*0/,CSEG/MXSG*' '/

      CALL ENTERS('PARSE',*9999)

      IEND=LEN_COMMAND
      DO WHILE(COMMAND_LINE(IEND).EQ.' '.AND.IEND.GT.1)
        IEND=IEND-1
      ENDDO !i
      IF(COMMAND_LINE(IEND).NE.' ') THEN
C KAT 28/3/00: Output interpreted command if requested
        IF(ECHO_INTERP_COM) THEN
          OP_STRING(1)='|'
          DO i=1,IEND
            OP_STRING(1)(1+i:1+i)=COMMAND_LINE(i)
          ENDDO !i
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
CC**** If IP=1 the command_line is written to IOOP.
CC ***   Write command_line to IOOP if IP=1
C        IF(ip.EQ.1) THEN
C          CALL STRING_TRIM(COMMAND_LINE,IBEG,IEND)
C          WRITE(OP_STRING,'('' > '',A)') COMMAND_LINE(IBEG:IEND)
C          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C        ENDIF
CC ***   Store command string in buffer
C        CALL ADDSTRTOBUFF(COMMAND_LINE,ERROR,*9999)

C ***   Split string into commands, qualifiers and parameters
        CALL SPLIT(IEND,COMMAND_LINE,ERROR,*9999)
        IF(TR04) THEN
          CALL SHOWCO(ERROR,*9999)
        ENDIF

C ***   Carry out commands
        CALL SYNTAX(ISEG,CSEG,QUIT,ERROR,*9999)
      ENDIF

      CALL EXITS('PARSE')
      RETURN
 9999 CALL ERRORS('PARSE',ERROR)
      CALL EXITS('PARSE')
      RETURN 1
      END


C old MLB 2 Sept 1997
C      SUBROUTINE PARSE_CLASS(noco,NTCO,nxc,CO,ERROR,*)
C
CC#### Subroutine: PARSE_CLASS
CC###  Description:
CC###    PARSE_CLASS parses command CO for keyword 'class'
CC###    and returns nxc.
CC**** Default is 1. Max number is 9.
C
C      IMPLICIT NONE
C!     Parameter List
C      INTEGER noco,NTCO,nxc
C      CHARACTER CO(*)*(*),ERROR*(*)
C!     Local Variables
C      INTEGER IFROMC,n3co
C      LOGICAL CBBREV
C
C      CALL ENTERS('PARSE_CLASS',*9999)
C      IF(CBBREV(CO,'CLASS',1,noco+1,NTCO,n3co)) THEN
C        nxc=IFROMC(CO(n3co+1))
C        CALL ASSERT(nxc.LE.9,' >>Too many classes',ERROR,*9999)
C      ELSE
C        nxc=1
C      ENDIF
C
C      CALL EXITS('PARSE_CLASS')
C      RETURN
C 9999 CALL ERRORS('PARSE_CLASS',ERROR)
C      CALL EXITS('PARSE_CLASS')
C      RETURN 1
C      END


