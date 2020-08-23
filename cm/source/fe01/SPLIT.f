      SUBROUTINE SPLIT(LEN_COMMAND,COMMAND_LINE,ERROR,*)

C#### Subroutine: SPLIT
C###  Description:
C###    SPLIT calls SUBSTR and splits the string COMMAND_LINE into NTCO
C###    commands CO, each with NOCOQU command qualifiers COQU.
C**** COMMAND_LINE is intent(in) and should not be modified.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'trac00.cmn'
!     Parameter List
      INTEGER LEN_COMMAND
      CHARACTER COMMAND_LINE(LEN_COMMAND)
      CHARACTER*(*) ERROR
!     Local Variables
      INTEGER i,IBEG,IBEG1,IEND,IEND1,nocoqu
!     Functions
      INTEGER LEN_TRIM

      CALL ENTERS('SPLIT',*9999)

C KAT 20Nov00: user defined strings are now obselete.
C      CALL STRING_TRIM(STRING,IBEG,IEND)
C      IPOS_space=CLOCAT(' ',STRING(IBEG:IEND))
C      IF(IPOS_space.EQ.0) IPOS_space=IEND+1
CC      CALL CUPPER(STRING(IBEG:IPOS_space-1),BUFFER)
C      CALL CUPPER(STRING(IBEG:IBEG+IPOS_space-1),BUFFER)
C      CALL STRING_TRIM(BUFFER,IBEG,IEND)
C      IF(.NOT.ABBRV(BUFFER(IBEG:IEND),'ASSIGN',2).AND.
C     '   .NOT.ABBRV(BUFFER(IBEG:IEND),'DEASSIGN',3)) THEN
CC       Substitute user defined strings into STRING for all commands
CC       except the ASSIGN and DEASSIGN commands
C        CALL USER(STRING,' ;,.:=+-*/()[]^',ERROR,*9999)
C      ENDIF
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' COMMAND_LINE is:'',9999A)')
     '    (COMMAND_LINE(i),i=1,LEN_COMMAND)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF
C     Put substrings in SB & number in NTSB
C     Separators are 'space,equals,semicolon'
      CALL SUBSTR(LEN_COMMAND,COMMAND_LINE,' =;',ERROR,*9999)
C KAT 20Nov00: user defined strings are now obselete.
C      IF(.NOT.ABBRV(BUFFER(IBEG:IEND),'ASSIGN',2).AND.
C     '   .NOT.ABBRV(BUFFER(IBEG:IEND),'DEASSIGN',3)) THEN
CC       Concatenate any strings separated by // unless the command is
CC       assign or deassign, in which case concatenation is done within
CC       ASSIGN and DEASSI, resp.
C        DO nosb=1,NTSB
C          END=.FALSE.
C          DO WHILE (.NOT.END)
C            IPOS=CLOCAT('//',SB(nosb))
C            IF(IPOS.GT.0) THEN
C              CALL STRING_TRIM(SB(nosb),IBEG,IEND)
C              SB(nosb)=SB(nosb)(1:IPOS-1)//SB(nosb)(IPOS+2:IEND)
C            ELSE
C              END=.TRUE.
C            ENDIF
C          ENDDO
Cc         IF(DOP) THEN
Cc      WRITE(OP_STRING,'('' ipos='',I2,'' SB('',I2,'')='',A)')
Cc    '        IPOS,nosb,SB(nosb)(1:20)
Cc            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
Cc         ENDIF
C        ENDDO
C      ENDIF
      DO noco=1,25
        CO(noco)=' '
        DO nocoqu=1,25
          COQU(noco,nocoqu)=' '
        ENDDO
      ENDDO

      IF(TR02) THEN
        WRITE(OP_STRING,*) ' NTSB=',NTSB
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,*) '   SB=',(SB(nosb),nosb=1,NTSB)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF
      nosb=1
      DO 2 noco=1,25
        IF(TR02) THEN
          WRITE(OP_STRING,*) ' NOCO=',noco
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,*) ' nosb=',nosb
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        IF(nosb.LE.NTSB) THEN
          IF(SB(nosb)(1:1).EQ.' '.OR.SB(nosb)(1:1).EQ.'='.OR.
     '      (SB(nosb)(1:1).EQ.'"')) THEN
            IF(SB(nosb)(2:2).EQ.'!'.OR.SB(nosb)(2:2).EQ.'#') THEN
              NTCO=noco-1
              GO TO 3
            ENDIF
            CO(noco)=SB(nosb)(2:)
            nosb=nosb+1
            DO nocoqu=1,25
              IF(TR02) THEN
                WRITE(OP_STRING,*) ' nocoqu=',nocoqu
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,*) '   nosb=',nosb
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
              IF(nosb.LE.NTSB) THEN
                IF(SB(nosb)(1:1).EQ.';') THEN
C KAT 24/3/00: checking truncation
                  IEND=LEN_TRIM(SB(nosb)(2:))
                  IF(IEND.LE.LEN_COQU) THEN
                    COQU(noco,nocoqu)(:IEND)=SB(nosb)(2:IEND+1)
                    nosb=nosb+1
                  ELSE
                    ERROR='Command qualifier `'//SB(nosb)(2:IEND+1)//
     '                ''' truncated'
                    GOTO 9999
                  ENDIF
                ELSE
                  NTCOQU(noco)=nocoqu-1
                  GOTO 2
                ENDIF
              ELSE
                NTCO=noco
                NTCOQU(noco)=nocoqu-1
                GOTO 3
              ENDIF
            ENDDO
          ELSE
            ERROR=' Command could not be parsed'
            GOTO 9999
          ENDIF
        ELSE
          NTCO=noco-1
          GOTO 3
        ENDIF
 2    CONTINUE

C     Trim off leading/trailing blanks/tabs/null characters
 3    DO noco=1,NTCO
        CALL STRING_TRIM(CO(noco),IBEG,IEND)
        CO(noco)=CO(noco)(IBEG:IEND)
        DO nocoqu=1,NTCOQU(noco)
          CALL STRING_TRIM(COQU(noco,nocoqu),IBEG1,IEND1)
          COQU(noco,nocoqu)=COQU(noco,nocoqu)(IBEG1:IEND1)
        ENDDO
      ENDDO

      CALL EXITS('SPLIT')
      RETURN
 9999 CALL ERRORS('SPLIT',ERROR)
      CALL EXITS('SPLIT')
      RETURN 1
      END



