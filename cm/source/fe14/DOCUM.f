      SUBROUTINE DOCUM(FILE,EXTEN,STRING,ERROR,*)

C#### Subroutine: DOCUM
C###  Description:
C###    DOCUM displays documentation.

C**** FILE identifies file to be read (file ID=8 and workstation ID=8)
C**** EXTEN is name of file extension (eg .DOC)
C**** STRING identifies string in documentation file at beginning & end
C****   of text to be displayed
C**** Documentation has root plus number of branches (max 20)
C**** Root and branches can each have arbitrary number of windows
C**** NT_WIND is number of windows for root or branch
C**** NT_LINE is number of lines in root or branch (this is divided by 23
C****   to give the number of windows)

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cmis00.cmn'
      INCLUDE 'disp00.cmn'
      INCLUDE 'docd00.cmn'
      INCLUDE 'inout00.cmn'
!     Parameter List
      CHARACTER ERROR*(*),FILE*(*),EXTEN*(*),STRING*(*)
!     Local Variables
      INTEGER i,IBEG1,IBEG2,IEND1,IEND2,
     '  LINE0,LINE1,LINE2,LINE3,LINE_MAX,NOWIND,NT_LINE,NT_WIND
      REAL*8 P(3)
      CHARACTER STRG*255,TOPIC*20
      LOGICAL CONTINUE,EXAMPLES,
     '  FIND_COMMAND,FIND_COMMON_BLOCK,FIND_COMMON_VARIABLE,
     '  FIND_INCLUDE_FILE,FIND_SUBROUTINE,FIND_VARIABLE,
     '  FOUND,READ_COM_FILE,ROOT
C     INTEGER IBEG4,IEND4,KEY
C     LOGICAL KEY5

      CALL ENTERS('DOCUM',*9999)
      DOCDIR=.TRUE.
C KAT 7Nov00: Not used
C      DO_EXAMPLE=.FALSE.
C      SELECT_EXAMPLE=.FALSE.

C      KEY3 =.FALSE.
C      KEY4 =.FALSE.
C      KEY5 =.FALSE.
C      KEY32=.FALSE.

      TOPIC=FILE
      IF(TOPIC(1:8).EQ.'EXAMPLES') THEN
        EXAMPLES=.TRUE.
        LINE_MAX=1000
      ELSE IF(TOPIC(1:11).EQ.'SUBROUTINES') THEN
        EXAMPLES=.FALSE.
        LINE_MAX=10000
      ELSE IF(TOPIC(1:7).EQ.'MODULES') THEN
        EXAMPLES=.FALSE.
        LINE_MAX=2000
      ELSE IF(TOPIC(1:13).EQ.'COMMON_BLOCKS') THEN
        EXAMPLES=.FALSE.
        LINE_MAX=1500
      ELSE
        EXAMPLES=.FALSE.
        LINE_MAX=1000
      ENDIF

      CALL STRING_TRIM(TOPIC,IBEG1,IEND1)
      CALL STRING_TRIM(EXTEN,IBEG2,IEND2)
C      CALL OPENF(8,'DISK','PRODUCT_1:[PRODUCT.CMISS.DOCUMENT]'
C     '  //TOPIC(IBEG1:IEND1)//'.'
C     '  //EXTEN(IBEG2:IEND2),'OLD','DIRECT','FORMATTED',132,
C     '  ERROR,*9999)
      WRITE(OP_STRING,'('' [Prev Screen] scrolls up one page'')')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' [Next Screen] scrolls down one page'')')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      IF(EXAMPLES) THEN
        WRITE(OP_STRING,'(''          [Do] exits from help and '
     '    //'reads current example file'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''      [Select] exits from help and '
     '    //'opens current example file'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''      [Insert] inserts current example '
     '    //'into local directory'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C KAT 2001-04-10: Not used
C        KEY3 =.TRUE.
C        KEY4 =.TRUE.
C        KEY5 =.TRUE.
      ENDIF
      WRITE(OP_STRING,'(''      [Delete] returns to top page'')')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''      [Remove] exits from help'')')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ROOT=.TRUE.

      FIND_COMMAND     =.FALSE.
      FIND_COMMON_BLOCK=.FALSE.
      FIND_COMMON_VARIABLE=.FALSE.
      FIND_INCLUDE_FILE=.FALSE.
      FIND_SUBROUTINE  =.FALSE.
      FIND_VARIABLE    =.FALSE.
      IF(TOPIC(1:8).EQ.'COMMANDS'.
     '  AND.STRING(1:7).EQ.'Command') THEN
        FIND_COMMAND     =.TRUE.      !to display a particular command
      ELSE IF(TOPIC(1:13).EQ.'COMMON_BLOCKS') THEN
        IF(STRING(1:7).EQ.'Include') THEN
          FIND_INCLUDE_FILE=.TRUE.    !to display a particular include file
        ELSE IF(STRING(1:8).EQ.'COMMON /') THEN
          FIND_COMMON_BLOCK=.TRUE.    !to display a particular commonblock
        ELSE IF(STRING(1:8).EQ.'VARIABLE') THEN
          FIND_COMMON_VARIABLE=.TRUE. !to display a commonblock variable
        ENDIF
      ELSE IF(TOPIC(1:11).EQ.'SUBROUTINES'.
     '  AND.STRING(1:10).EQ.'Subroutine') THEN
        FIND_SUBROUTINE  =.TRUE.      !to display a particular subroutine
      ELSE IF(TOPIC(1:9).EQ.'VARIABLES'.
     '  AND.STRING(1:8).EQ.'Variable') THEN
        FIND_VARIABLE    =.TRUE.      !to display a particular variable
      ENDIF
C      STRG=CUPPER(STRING)
      CALL CUPPER(STRING,STRG)
      CDATA(2)=STRG

C KAT 2001-04-10: Not used
C      OPTION1=.FALSE.
C      OPTION2=.FALSE.
C      OPTION3=.FALSE.
C      OPTION4=.FALSE.
C      TERM1=' '
C      TERM2=' '
C      TERM3=' '
C      TERM4=' '
C      TERM5=' '
      IF(.NOT.USE_SOCKET) THEN
        CALL ACWK(8,0,ERROR,*9999)
      ENDIF

C 10   CONTINUE
      CALL STRING_TRIM(CDATA(2),IBEG2,IEND2)
      IF(DOP) THEN
        WRITE(OP_STRING,'('' String: '',A)') CDATA(2)(IBEG2:IEND2)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

C KAT 2001-04-10: Not used
C      IF(EXAMPLES) THEN !enable numeric key press
C        KEY32=.TRUE.
C      ENDIF
      IF(ROOT) THEN
        IF(FIND_COMMAND) THEN              !Find command name
          CALL FIND(8,'Command    '//CDATA(2)(12:IEND2),1,
     '      LINE_MAX,LINE1,FOUND,ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' LINE1=',LINE1,' FOUND=',FOUND
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          endif
          CALL FIND(8,'Command',LINE1+1,LINE1+50,LINE2,FOUND,
     '      ERROR,*9999)
          LINE1=LINE1-1 !to ensure that command name is displayed
        ELSE IF(FIND_COMMON_BLOCK) THEN    !Find common-block name
          CALL FIND(8,'COMMON '//CDATA(2)(8:IEND2),1,
     '      LINE_MAX,LINE1,FOUND,ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' LINE1=',LINE1,' FOUND=',FOUND
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CALL FIND(8,'Include',LINE1,LINE1-20,LINE1,FOUND,ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' LINE1=',LINE1,' FOUND=',FOUND
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          LINE1=LINE1-1
          LINE2=LINE1+22
        ELSE IF(FIND_COMMON_VARIABLE) THEN !Find common-block variable
          CALL FIND(8,CDATA(2)(10:IEND2),1,
     '      LINE_MAX,LINE1,FOUND,ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' LINE1=',LINE1,' FOUND=',FOUND
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CALL FIND(8,'Include',LINE1,LINE1-20,LINE1,FOUND,ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' LINE1=',LINE1,' FOUND=',FOUND
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          LINE1=LINE1-1
          LINE2=LINE1+22
        ELSE IF(FIND_INCLUDE_FILE) THEN    !Find include-file name
          CALL FIND(8,'Include '//CDATA(2)(9:IEND2)//'.cmn',1,
     '      LINE_MAX,LINE1,FOUND,ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' LINE1=',LINE1,' FOUND=',FOUND
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          LINE1=LINE1-1
          LINE2=LINE1+22
        ELSE IF(FIND_SUBROUTINE) THEN      !Find subroutine name
          CALL FIND(8,'Subroutine    '//CDATA(2)(12:IEND2),1,
     '      LINE_MAX,LINE1,FOUND,ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' LINE1=',LINE1,' FOUND=',FOUND
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CALL FIND(8,'Subroutine',LINE1+1,LINE1+50,LINE2,FOUND,
     '      ERROR,*9999)
          LINE1=LINE1-1 !to ensure that subroutine name is displayed
        ELSE IF(FIND_VARIABLE) THEN        !Find variable name
          CALL FIND(8,CDATA(2)(10:IEND2),1,
     '      LINE_MAX,LINE1,FOUND,ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' LINE1=',LINE1,' FOUND=',FOUND
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          LINE1=LINE1-1 !to ensure that subroutine name is displayed
          LINE2=LINE1+22
        ELSE
          CALL FIND(8,CDATA(2)(IBEG2:IEND2),1,LINE_MAX,LINE1,FOUND,
     '      ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' LINE1=',LINE1,' FOUND=',FOUND
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CALL FIND(8,CDATA(2)(IBEG2:IEND2),LINE1+1,LINE1+LINE_MAX,
     '      LINE2,FOUND,ERROR,*9999)
          IF(TOPIC(1:7).EQ.'MODULES') LINE1=LINE1-1 !to ensure name is displayed
        ENDIF
        IF(DOP) THEN
          WRITE(OP_STRING,*) ' LINE2=',LINE2,' FOUND=',FOUND
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        LINE0=LINE2
        ROOT=.FALSE.

      ELSE IF(.NOT.ROOT) THEN
        CALL FIND(8,CDATA(2)(IBEG2:IEND2),LINE0,LINE_MAX,LINE1,FOUND,
     '    ERROR,*9999)
        IF(DOP) THEN
          WRITE(OP_STRING,*) ' LINE1=',LINE1,' FOUND=',FOUND
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        CALL FIND(8,CDATA(2)(IBEG2:IEND2),LINE1+1,LINE1+LINE_MAX,LINE2,
     '    FOUND,ERROR,*9999)
        IF(DOP) THEN
          WRITE(OP_STRING,*) ' LINE2=',LINE2,' FOUND=',FOUND
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF

      READ(8,FMT='(A)',REC=LINE1+1) CDATA(1)(1:80)
      IF(DOP) THEN
        WRITE(OP_STRING,'(1X,A)') CDATA(1)(1:80)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      IF(CDATA(1)(2:13).EQ.'FILE=Example') THEN !read com file
        READ_COM_FILE=.TRUE.
        CALL STRING_TRIM(CDATA(1),IBEG1,IEND1)
        CALL STRING_TRIM(EXAMPLES_DIR,IBEG2,IEND2)
C!!! same unit as IOOUT for `set output'!
        CALL OPENF(9,'DISK',EXAMPLES_DIR(IBEG2:IEND2)//CDATA(1)(7:25),
     '    'OLD','DIRECT','FORMATTED',132,ERROR,*9999)
        LINE1=1
        CALL FIND(9,'fem quit',1,LINE_MAX,LINE2,FOUND,ERROR,*9999)
        IF(DOP) THEN
          WRITE(OP_STRING,*) ' LINE2=',LINE2,' FOUND=',FOUND
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        NT_LINE=LINE2
C KAT 2001-04-10: Not used
C        KEY32=.FALSE. !to disable numeric keys
      ELSE !continue reading current file
        READ_COM_FILE=.FALSE.
        NT_LINE=LINE2-LINE1-1
      ENDIF
      NT_WIND=1+INT((NT_LINE-1)/23)
      IF(DOP) THEN
        WRITE(OP_STRING,*) ' NT_WIND=',NT_WIND
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      LINE3=0
      NOWIND=1
      CONTINUE=.TRUE.
      DO WHILE(CONTINUE)
        DO i=1,MIN(23,NT_LINE-LINE3)
          IF(READ_COM_FILE) THEN
            READ(9,FMT='(A)',REC=LINE3+i) CDATA(1)(1:80)
          ELSE
            READ(8,FMT='(A)',REC=LINE1+LINE3+i) CDATA(1)(1:80)
          ENDIF
          IF(DOP) THEN
            WRITE(OP_STRING,'(1X,A)') CDATA(1)(1:80)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(USE_SOCKET) THEN
            CALL WRITES(IOH2,CDATA(1)(1:80),ERROR,*9999)
          ELSE
            P(1)=0.d0
            P(2)=1.d0-DBLE(i)*0.025d0
            CALL TEXT(1,8,CDATA(1)(1:80),P,ERROR,*9999)
          ENDIF
        ENDDO
C KAT 2Mar01: does nothing.  KEY not set.
C        CALL GETSTR3(KEY,ERROR,*9999)
C        IF(KEY.EQ.1.AND.NOWIND.GT.1) THEN           !Previous screen key
C          NOWIND=NOWIND-1
C        ELSE IF(KEY.EQ.2.AND.NOWIND.LT.NT_WIND) THEN !Next screen key
C          NOWIND=NOWIND+1
C        ELSE IF(KEY3.AND.KEY.EQ.3) THEN             !Do key
CC KAT 7Nov00: Not used
CC          DO_EXAMPLE=.TRUE.
C          CONTINUE=.FALSE.
C        ELSE IF(KEY4.AND.KEY.EQ.4) THEN             !Select key
CC KAT 7Nov00: Not used
CC          SELECT_EXAMPLE=.TRUE.
C          CONTINUE=.FALSE.
CC KAT 2Mar01: COPY_FILE does nothing
CC        ELSE IF(KEY5.AND.KEY.EQ.5) THEN             !Insert key
CC          COPYFILE='EXAMPLE_'//TERM1//TERM2//TERM3
CC          CALL COPY_FILE(COPYFILE)
CC          WRITE(OP_STRING,
CC     '      '(1X,A,'' copied to local directory'')') COPYFILE
CC          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C        ELSE IF(KEY.EQ.6) THEN                      !Delete key
C          CDATA(2)=STRG
C          IF(.NOT.USE_SOCKET) THEN
Cc            CALL GKS_CLRWK(8,GALWAY,ERROR,*9999)
C          ENDIF
C          ROOT=.TRUE.
C          OPTION1=.FALSE.
C          OPTION2=.FALSE.
C          OPTION3=.FALSE.
C          OPTION4=.FALSE.
C          TERM1=' '
C          TERM2=' '
C          TERM3=' '
C          TERM4=' '
C          TERM5=' '
C          GO TO 10
C        ELSE IF(KEY.EQ.7) THEN                      !Remove key
C          CONTINUE=.FALSE.
C        ELSE IF(KEY32.AND.KEY.GE.32) THEN           !Numeric key
C          CALL STRING_TRIM(STRG,IBEG4,IEND4)
C          IF(.NOT.OPTION1) THEN      !this is the 1st numeric key press
C            OPTION1=.TRUE.
C            TERM1=CHAR(KEY)
C            CDATA(2)=STRG(IBEG4:IEND4)//' '//TERM1
CC KAT 7Nov00: Not used
CC            EXAMPLE_NAME='EXAMPLE_'//TERM1
C          ELSE IF(.NOT.OPTION2) THEN !this is the 2nd numeric key press
C            OPTION2=.TRUE.
C            TERM2=CHAR(KEY)
C            CDATA(2)=STRG(IBEG4:IEND4)//' '//TERM1//TERM2
CC KAT 7Nov00: Not used
CC            EXAMPLE_NAME='EXAMPLE_'//TERM1//TERM2
C          ELSE IF(.NOT.OPTION3) THEN !this is the 3rd numeric key press
C            OPTION3=.TRUE.
C            TERM3=CHAR(KEY)
C            CDATA(2)=STRG(IBEG4:IEND4)//' '//TERM1//TERM2//TERM3
CC KAT 7Nov00: Not used
CC            EXAMPLE_NAME='EXAMPLE_'//TERM1//TERM2//TERM3
C          ELSE IF(.NOT.OPTION4) THEN !this is the 4th numeric key press
C            OPTION4=.TRUE.
C            TERM4=CHAR(KEY)
C            CDATA(2)=STRG(IBEG4:IEND4)//' '//TERM1//TERM2//TERM3//TERM4
CC KAT 7Nov00: Not used
CC            EXAMPLE_NAME='EXAMPLE_'//TERM1//TERM2//TERM3//TERM4
C          ELSE                       !this is the 5th numeric key press
C            TERM5=CHAR(KEY)
C            CDATA(2)=STRG(IBEG4:IEND4)//' '//TERM1//TERM2//TERM3//
C     '        TERM4//TERM5
CC KAT 7Nov00: Not used
CC            EXAMPLE_NAME='EXAMPLE_'//TERM1//TERM2//TERM3//TERM4//TERM5
C          ENDIF
C          IF(.NOT.USE_SOCKET) THEN
Cc            CALL GKS_CLRWK(8,GALWAY,ERROR,*9999)
C          ENDIF
C          GO TO 10
C        ENDIF
        LINE3=23*(NOWIND-1)
C        IF(.NOT.USE_SOCKET) THEN
c          CALL GKS_CLRWK(8,GALWAY,ERROR,*9999)
C        ENDIF
      ENDDO
      IF(.NOT.USE_SOCKET) THEN
        CALL DAWK(8,0,ERROR,*9999)
        CALL CLOSE_WS(8,ERROR,*9999)     !new
! old   CALL COWK(8,ERROR,*9999)     !GBS 3/8/92
      ENDIF
      CALL CLOSEF(8,ERROR,*9999)
      CALL CLOSEF(9,ERROR,*9999)
      DOCDIR=.FALSE.

      CALL EXITS('DOCUM')
      RETURN
 9999 CALL ERRORS('DOCUM',ERROR)
      CALL EXITS('DOCUM')
      CALL DAWK(8,0,ERROR,*9999)
      CALL CLOSE_WS(8,ERROR,*9999)     !new
! old CALL COWK(8,ERROR,*9999)     !GBS 3/8/92
      CLOSE(UNIT=8)
      CLOSE(UNIT=9)
      RETURN 1
      END


