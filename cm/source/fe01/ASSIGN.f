      SUBROUTINE ASSIGN(STRING,noco,CO,COQU,ERROR,*)

C#### Subroutine: ASSIGN
C###  Description:
C###    ASSIGN defines a character string and places it within the list
C###    of user variables.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'mxch.inc'
      INCLUDE 'user00.cmn'
!     Parameter List
      INTEGER noco
      CHARACTER CO(*)*(*),COQU(25,*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER CLOCAT,IBEG,IBEG1,IBEG2,IBEG3,IEND,IEND1,IEND2,IEND3,
     '  IEXPR,IPOS,N1US,nous,NTITLV(20),NTLV
      REAL*8 RA(100),RFROMC
      CHARACTER ASSIGN_STR*200,STRG*40
      LOGICAL ABBREV,END,EXPR,PROMPT

      CALL ENTERS('ASSIGN',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
C news MPN 27Jul97: new style help
C---------------------------------------------------------------------

C#### Command: assign<;p>' STRING_NAME=(STRING/`EXPRESSION`)'
C###  Parameter:     <truncate>
C###  Description:
C###    Assigns the user defined variable STRING_NAME to be STRING
C###    or an arithmetic REAL*8 EXPRESSION delimited by ` `.
C###    The 'truncate' option may be used to truncate the
C###    REAL*8 EXPRESSION to an integer value before is is assigned
C###    to STRING_NAME.

        OP_STRING(1)=STRING(1:IEND)//'<;p>'
     '    //' STRING_NAME=(STRING/`EXPRESSION`)'
        OP_STRING(2)=BLANK(1:IEND)//'<truncate>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
C new end
C old        WRITE(OP_STRING,'(1X,A)') STRING(1:IEND)//'<;p>'
C old     '    //' STRING_NAME=STRING'
C old        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ELSE
        IF(ABBREV(COQU(noco,1),'P',1)) THEN
          PROMPT=.TRUE.
        ELSE
          PROMPT=.FALSE.
        ENDIF
        CALL STRING_TRIM(CO(noco+1),IBEG1,IEND1)
        IF(PROMPT) THEN
          WRITE(OP_STRING,'($,'' Enter '',A,'': '')')
     '      CO(noco+1)(IBEG1:IEND1)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          READ(IOIP,'(A)') STRG
          CO(noco+2)=STRG
        ENDIF
        CALL STRING_TRIM(CO(noco+2),IBEG2,IEND2)
C news MPN 28Jul97
C       If assign string starts and ends with a ` then it is an
C       expression, so convert it.
        IF(CO(noco+2)(IBEG2:IBEG2).EQ.'`'.AND.
     '    CO(noco+2)(IEND2:IEND2).EQ.'`') THEN
          ASSIGN_STR=CO(noco+2)(IBEG2+1:IEND2-1)
          EXPR=.TRUE.
        ELSE
          ASSIGN_STR=CO(noco+2)(IBEG2:IEND2)
          EXPR=.FALSE.
        ENDIF
C       Substitute previous user defined strings into current
C       string definition
        CALL USER(ASSIGN_STR,' ;,.:=+-*/()[]^',ERROR,*9999)
C       Concatenate any strings separated by //
        END=.FALSE.
        DO WHILE (.NOT.END)
          IPOS=CLOCAT('//',ASSIGN_STR)
          IF(IPOS.GT.0) THEN
            CALL STRING_TRIM(ASSIGN_STR,IBEG,IEND)
            ASSIGN_STR=ASSIGN_STR(1:IPOS-1)//ASSIGN_STR(IPOS+2:IEND)
          ELSE
            END=.TRUE.
          ENDIF
        ENDDO !while
        IF(EXPR) THEN
C         Evaluate the expression
          CALL PARSTR(ASSIGN_STR,20,NTLV,NTITLV,100,RA,ERROR,*9999)
          IF(ABBREV(CO(noco+3),'TRUNCATE',1)) THEN
C           truncate the expression to a simple integer value
            IEXPR=INT(RFROMC(ASSIGN_STR))
            WRITE(ASSIGN_STR,'(I20)') IEXPR
          ENDIF
        ENDIF
        CALL STRING_TRIM(ASSIGN_STR,IBEG2,IEND2)
C end new
        N1US=0
        DO nous=1,NTUS
C new MPN 29Jan98: make sure user variable ID's exactly match
C                  (ie. don't allow abbreviations)
          CALL STRING_TRIM(USID(nous),IBEG3,IEND3)
          IF(USID(nous)(IBEG3:IEND3).EQ.
     '      CO(noco+1)(IBEG1:IEND1)) THEN
C old
C old !old      IF(USID(nous)(1:IEND1-IBEG1+1).EQ.
C old !old '      CUPPER(CO(noco+1)(IBEG1:IEND1))) THEN
C old !news     AAY don't convert to upper case
C old           IF(USID(nous)(1:IEND1-IBEG1+1).EQ.
C old     '      CO(noco+1)(IBEG1:IEND1)) THEN
            N1US=nous
            GO TO 100
          ENDIF
        ENDDO
 100    IF(N1US.EQ.0) THEN
          CAll ASSERT(NTUS.LT.NXUS,
     '      '>>>ERROR: Too many user-defined variables. Increase NXUS',
     '      ERROR,*9999)
          NTUS=NTUS+1
          N1US=NTUS
        ENDIF
!old    USID(N1US)=CUPPER(CO(noco+1)(IBEG1:IEND1))
        USID(N1US)=CO(noco+1)(IBEG1:IEND1)
        LNUSID(N1US)=IEND1-IBEG1+1
!old    US(N1US)=CUPPER(CO(noco+2)(IBEG2:IEND2))
!old    US(N1US)=CO(noco+2)(IBEG2:IEND2)
        US(N1US)=ASSIGN_STR(IBEG2:IEND2)
        LNUS(N1US)=IEND2-IBEG2+1
        ACUS(N1US)=.TRUE.
      ENDIF

      CALL EXITS('ASSIGN')
      RETURN
 9999 CALL ERRORS('ASSIGN',ERROR)
      CALL EXITS('ASSIGN')
      RETURN 1
      END


