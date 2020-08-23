      SUBROUTINE GINOUT(IO_TYPE,IPTYPE,IFREE,IFIXED,N1dummy,N2dummy,
     '  NOQUES,FILEIP,G_FORMAT,NDATA,A_DATA,A_DEFLT,C_DATA,C_DEFLT,
     '  C_LEN,I_DATA,I_DEFLT,IMIN,IMAX,L_DATA,L_DEFLT,R_DATA,R_DEFLT,
     '  RMIN,RMAX,INFO,ERROR,*)

C#### Subroutine: GINOUT
C###  Description:
C###    <HTML>
C###    GINOUT handles i/o of character,integer,logical and
C###    real*8 values.
C###    <PRE>
C###    If IO_TYPE=0 prompts written to IFREE
C###     "    "    1   "     written to IFREE and (with data) to IFIXED
C###     "    "    2   "     read with data from IFIXED
C###     "    "    3   "     written with data to IFIXED
C###     "    "    4   "     read with data from IFIXED & written to IFREE
C###     "    "    5   default values are used
C###    If IPTYPE=IPMESS a message is printed to the terminal only
C###     "    "   IPANSW input is answer(yes/no)
C###     "    "   IPCHAR   "    " character
C###     "    "   IPINTE   "    " integer
C###     "    "   IPLOGI   "    " logical
C###     "    "   IPREAL   "    " real*8
C###    NOQUES is #questions since start of command (sb zero initially)
C###    If IFIXED is not open no IO is performed on that file.
C###    Note: data may include lists & group name. eg 1,4,2,5..9:2,epi,17.
C###    If group name is used, the calling routine must set C_DATA to
C###    'NODES','ELEMENTS','FACES','LINES' or 'GRIDS'.
C###    N1dummy and N2dummy are unused and unaltered.
C###   
C###    IPTYPE=IPANSW: (answers)
C###      Answers separated by commas or blanks are read into
C###      A_DATA(NDAT) NDAT=1,NDATA.
C###    IPTYPE=IPCHAR: (characters)
C###      Strings separated by commas read into C_DATA(NDAT) NDAT=1,NDATA
C###      Strings are truncated after C_LEN characters.
C###    IPTYPE=IPINTE: (integers)
C###      Integer values separated by commas or blanks are read into
C###      I_DATA(NDAT) NDAT=1,NDATA.
C###      The limiting values for the read are defined by IMIN and IMAX.
C###    IPTYPE=IPLOGI: (logicals)
C###      Logical values separated by commas or blanks are read into
C###      L_DATA(NDAT) NDAT=1,NDATA.
C###    IPTYPE=IPREAL: (real*8s)
C###      Real values separated by commas or blanks are read into
C###      R_DATA(NDAT) NDAT=1,NDATA.
C###   
C###    If an end of file or carriage return is detected while attempting
C###    to read *DATA(NDAT) it is assigned the value *DEFLT(NDAT).
C###    FILEIP is .true. if response indicates file input.
C###    If an error is detected a diagnostic message is returned in ERROR
C###    and control is returned to the statement number of the asterisk.
C###   </PRE>
C###   </HTML>
      
C****  LINE_BUFFER contains line #s for previous FORMATS
C****  NDATA_BUFFER          "   NDATA  #s        "        "

      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='GINOUT')
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'disp00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
!     Parameter List
      INTEGER C_LEN,IFIXED,IFREE,IMAX,IMIN,INFO,IO_TYPE,IPTYPE,
     '  I_DATA(0:*),I_DEFLT(*),NDATA,NOQUES,N1dummy,N2dummy
      REAL*8 RMAX,RMIN,R_DATA(*),R_DEFLT(*)
      CHARACTER A_DATA(*)*(*),A_DEFLT(*)*(*),C_DATA(*)*(*),
     '  C_DEFLT(*)*(*),ERROR*(*),G_FORMAT*(*)
      LOGICAL FILEIP,L_DATA(*),L_DEFLT(*)
!     Local Variables
      INTEGER ERR,FORM_CHECK_BEG,FORM_CHECK_END,FORM_END,
     '  IBEG,IBEG1,IDATA_TEMP,
     '  IEND,IEND1,IOIM,IOSTAT,
     '  irec,irec_at_first_mismatch,
     '  irec_at_insert_point,irec_at_line_1_MCQ,irec_last_line_MCQ,
     '  irec_shift,irec_at_EOF,j,
     '  len_default_FORM,len_default_LINE,LINE_BUFFER,
     '  LINE_CHECK_BEG,LINE_CHECK_END,LINE_QUES_END,
     '  m,MAXTRY,n,nd,NDAT,NDATA_BUFFER,NDATA_RET,
     '  No_lines,No_lines_in_MCQ,No_lines_in_format,
     '  N1char,N2char,POS_dollar,POS_leftsqb1,
     '  POS_quote,POS_quote1,POS_quote2,POS_search,POS_2quote
      CHARACTER CIUNIT*10,CLINE*(MXCH),LINE*132,OUTFORMAT*800
      LOGICAL BLANKLINE,FILEIP_FIRST,FOUND_BLANK,FIRST_LOOP,
     '  FIRST_MISMATCH,LAST_LINE,
     '  MISMATCH_MESSAGE,MULTICHOICE,
     '  OPENED,REVERT1,REVERT2,SKIP_LINE
      ! Functions
      INTEGER CLOCAT,LEN_TRIM
      LOGICAL LFROMC

      DATA MAXTRY/8/

      CALL ENTERS(ROUTINENAME,*9999)

C GMH 13/2/97 Initialise variables
      irec_at_EOF=0
      irec_at_first_mismatch=0
      irec_at_line_1_MCQ=0

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'(/'' NDATA='',I3)') NDATA
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      FIRST_LOOP=.TRUE.

      IF(IPTYPE.EQ.IPANSW.AND.NDATA.GT.IOAMX) THEN
        ERROR=' NDATA exceeds the length of the array A_DATA'
        GOTO 9999
      ELSE IF(IPTYPE.EQ.IPCHAR.AND.NDATA.GT.IOCMX) THEN
        ERROR=' NDATA exceeds the length of the array C_DATA'
        GOTO 9999
      ELSE IF(IPTYPE.EQ.IPINTE.AND.NDATA.GT.IOIMX) THEN
        ERROR=' NDATA exceeds the length of the array I_DATA'
        GOTO 9999
      ELSE IF(IPTYPE.EQ.IPLOGI.AND.NDATA.GT.IOLMX) THEN
        ERROR=' NDATA exceeds the length of the array L_DATA'
        GOTO 9999
      ELSE IF(IPTYPE.EQ.IPREAL.AND.NDATA.GT.IORMX) THEN
        ERROR=' NDATA exceeds the length of the array R_DATA'
        GOTO 9999
      ENDIF
      REVERT1=.FALSE.
      REVERT2=.FALSE.

C     Prompted input
 5    IF(IO_TYPE.LE.1.OR.REVERT1.OR.REVERT2) THEN !prompted input
        NOQUES=NOQUES+1
        NDATA_BUFFER=NDATA
C KAT 2002-06-04: backup doesn't work anyway: removed buffers
C        BUFFER(1)=G_FORMAT
C        IPTYPE_BUFFER(1)=IPTYPE
C        IF(IPTYPE_BUFFER(1).EQ.1) THEN
C          DO nd=1,NDATA
C            A_DEFLT_BUFFER(nd,1)=A_DEFLT(nd)
C          ENDDO
C        ELSE IF(IPTYPE_BUFFER(1).EQ.2) THEN
C          ICHAR_BUFFER(1)=ICHAR
C          DO nd=1,NDATA
C            C_DEFLT_BUFFER(nd,1)=C_DEFLT(nd)
C          ENDDO
C        ELSE IF(IPTYPE_BUFFER(1).EQ.3) THEN
C          IMIN_BUFFER(1)=IMIN
C          IMAX_BUFFER(1)=IMAX
C          DO nd=1,NDATA
C            I_DEFLT_BUFFER(nd,1)=I_DEFLT(nd)
C          ENDDO
C        ELSE IF(IPTYPE_BUFFER(1).EQ.4) THEN
C          DO nd=1,NDATA
C            L_DEFLT_BUFFER(nd,1)=L_DEFLT(nd)
C          ENDDO
C        ELSE IF(IPTYPE_BUFFER(1).EQ.5) THEN
C          RMIN_BUFFER(1)=RMIN
C          RMAX_BUFFER(1)=RMAX
C          DO nd=1,NDATA
C            R_DEFLT_BUFFER(nd,1)=R_DEFLT(nd)
C          ENDDO
C        ENDIF !iptype_buffer
C        IBUFFER=1
        IF(IO_TYPE.EQ.1) THEN
          INQUIRE(UNIT=IFIXED,OPENED=OPENED,NEXTREC=irec)
          LINE_BUFFER=irec
        ENDIF
      ENDIF !io_type.le.1.or.revert1.or.revert2

C     Locate key characters in format string
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' G_FORMAT='',A)') G_FORMAT(1:80)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF
      N1char=3
      POS_search=2 !to skip leading (
      SKIP_LINE=G_FORMAT(POS_search:POS_search).EQ.'/' !Look for / at start
      IF(SKIP_LINE) THEN
        No_lines_in_format=1
        POS_search=POS_search+1
      ELSE
        No_lines_in_format=0
      ENDIF
      IF(G_FORMAT(POS_search:POS_search+1).EQ.'$,') THEN
        POS_dollar=POS_search
        POS_search=POS_search+2
      ELSE
        POS_dollar=0
      ENDIF
      POS_quote1=POS_search !first quote
      POS_2quote=POS_search !might be second to last quote
      DO WHILE(G_FORMAT(POS_search:POS_search).EQ.'''')
        POS_quote=CLOCAT('''',G_FORMAT(POS_search+1:))
        POS_quote2=POS_search+POS_quote !second quote
        POS_search=POS_quote2+1
      ENDDO !check that this isn't an apostrophe in the line
C     Record #lines in format prior to data line
      DO WHILE(G_FORMAT(POS_search:POS_search).EQ.'/')
        No_lines_in_format=No_lines_in_format+1
        POS_search=POS_search+1
        IF(G_FORMAT(POS_search:POS_search+1).EQ.'$,') THEN
          POS_dollar=POS_search
          POS_search=POS_search+2
        ENDIF
        POS_2quote=POS_search !might be first quote on last line
        DO WHILE(G_FORMAT(POS_search:POS_search).EQ.'''')
          POS_search=POS_search+1
          POS_quote=CLOCAT('''',G_FORMAT(POS_search:))
          POS_search=POS_search+POS_quote
        ENDDO !check that this isn't an apostrophy in the line
      ENDDO !check for more newlines
      IF(POS_dollar.EQ.0) POS_dollar=POS_2quote-2
      FORM_END=POS_search-2
      N2char=POS_search-POS_2quote !length of last line + 2
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' POS_dollar='',I4)') POS_dollar
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' N1&2char  ='',2I4)') N1char,N2char
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

C     Prompted input
 10   IF(IO_TYPE.LE.1.OR.REVERT1.OR.REVERT2) THEN !prompted input
        DO m=1,MAXTRY
          NDAT=NDATA_BUFFER
C KAT 2002-06-04: backup doesn't work anyway: removed buffers
C          G_FORMAT=BUFFER(IBUFFER)
C          IPTYP=IPTYPE_BUFFER(IBUFFER)
C          IF(IPTYP.EQ.1) THEN                      !yes/no data
C            DO nd=1,NDAT
C              A_DEFLT(nd)=A_DEFLT_BUFFER(nd,IBUFFER)
C            ENDDO
C          ELSE IF(IPTYP.EQ.2) THEN                 !character data
C            ICHARAC=ICHAR_BUFFER(IBUFFER)
C            DO nd=1,NDAT
C              C_DEFLT(nd)=C_DEFLT_BUFFER(nd,IBUFFER)
C            ENDDO
C          ELSE IF(IPTYP.EQ.3) THEN                 !integer data
C            IMINIM=IMIN_BUFFER(IBUFFER)
C            IMAXIM=IMAX_BUFFER(IBUFFER)
C            DO nd=1,NDAT
C              I_DEFLT(nd)=I_DEFLT_BUFFER(nd,IBUFFER)
C            ENDDO
C          ELSE IF(IPTYP.EQ.4) THEN                 !logical data
C            DO nd=1,NDAT
C              L_DEFLT(nd)=L_DEFLT_BUFFER(nd,IBUFFER)
C            ENDDO
C          ELSE IF(IPTYP.EQ.5) THEN                 !real data
C            RMINIM=RMIN_BUFFER(IBUFFER)
C            RMAXIM=RMAX_BUFFER(IBUFFER)
C            DO nd=1,NDAT
C              R_DEFLT(nd)=R_DEFLT_BUFFER(nd,IBUFFER)
C            ENDDO
C          ENDIF

C        Write out the prompt string
          CALL DISPLAY_PROMPT(G_FORMAT,ERROR,*9999)

          IF(IPTYPE.NE.IPMESS) THEN
            IF(IFREE.LT.0) THEN
              CALL DISPLAY_PROMPT('()',ERROR,*9999) !move to new line
              CALL FLAG_ERROR(0,'No interactive input available')
              ERROR=' '
              GOTO 9999
            ENDIF

C        Read in the response
            CALL GINPUT(IPTYPE,IFREE,NDAT,A_DATA,A_DEFLT,C_DATA,C_DEFLT,
     '        CLINE,I_DATA,I_DEFLT,IMIN,IMAX,L_DATA,L_DEFLT,
     '        R_DATA,R_DEFLT,RMIN,RMAX,FILEIP,FILEIP_FIRST,INFO,
     '        ERROR,*9999,*9995)
          ELSE IF(IPTYPE.EQ.IPMESS) THEN
            ERROR=' '
          ENDIF
          IF(ERROR.EQ.' ') THEN
            GOTO 2
          ELSE IF(ERROR(1:2).NE.' ?') THEN
            OP_STRING(1)=ERROR
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO
        ERROR=' More than 8 invalid input attempts made'
        GOTO 9999

 9995   BACKUP=.TRUE.
C KAT 2002-06-04: backup doesn't work anyway: removed buffers
C        IF(IBUFFER.LT.NOQUES.AND.IBUFFER.LT.NT_BUFFER) THEN
C          IBUFFER=IBUFFER+1
C          GO TO 10
C        ELSE
          WRITE(OP_STRING,'('' >>No further questions in buffer''/)')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          GO TO 10
C        ENDIF

    2   IF(IO_TYPE.EQ.1.OR.REVERT2) THEN !write prompted input to file
          irec=LINE_BUFFER
          IF(OPENED.AND.IFIXED.NE.0) THEN

            IF(REVERT2) THEN
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'('' ###Write to file: ifixed='',I2)')
     '            ifixed
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' G_FORMAT='',A)') G_FORMAT(1:60)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF !dop
C             Shift lines to make room for new question
              IF(FIRST_LOOP) THEN
                irec_at_insert_point=irec_at_first_mismatch
              ELSE
                irec_at_insert_point=irec_at_line_1_MCQ
              ENDIF
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'('' irec_at_insert_point='',I4)')
     '            irec_at_insert_point
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF
              DO irec=irec_at_EOF-1,irec_at_insert_point,-1
                READ(UNIT=IFIXED,FMT='(A)',REC=irec,IOSTAT=IOSTAT) LINE
                WRITE(UNIT=IFIXED,FMT='(A)',
     '            REC=irec+No_lines_in_format+1,IOSTAT=IOSTAT)  LINE
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                  WRITE(OP_STRING,'('' read from irec='',I3,'
     '              //''' write to irec+#lines+1='',I3,'
     '              //''' LINE: '',A)')
     '              irec,irec+No_lines_in_format+1,LINE(1:15)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                ENDIF
              ENDDO !irec
              irec_at_EOF=irec_at_EOF+No_lines_in_format+1
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'('' irec_at_EOF='',I4)') irec_at_EOF
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF
C             Set file record pointer to insert point of new question
              irec=irec_at_insert_point
            ENDIF !REVERT2

            IF(IPTYPE.EQ.IPMESS) THEN                  !output string
              WRITE(UNIT=IFIXED,FMT=G_FORMAT,REC=irec,IOSTAT=IOSTAT)

            ELSE IF(FILEIP) THEN
              G_FORMAT=G_FORMAT(1:N2char)//'FILE'
              WRITE(UNIT=IFIXED,FMT='(A)',REC=irec,IOSTAT=IOSTAT)
     '          G_FORMAT(N1char:N2char+4)
              IF(FILEIP_FIRST) THEN !prompt for filename
                WRITE(UNIT=IFIXED,FMT='('' Enter FILENAME: '',A)',
     '            REC=irec+1,IOSTAT=IOSTAT) FILE07
              ENDIF

            ELSE IF(.NOT.FILEIP) THEN !write to IPfile
              IF(N1char.EQ.0) THEN      !write integer values
                IF(IPTYPE.EQ.IPANSW) THEN                  !yes/no data
                  WRITE(UNIT=IFIXED,FMT=G_FORMAT,REC=irec,IOSTAT=IOSTAT)
     '              (A_DATA(n),n=1,NDAT)
                ELSE IF(IPTYPE.EQ.IPCHAR) THEN             !character data
                  WRITE(UNIT=IFIXED,FMT=G_FORMAT,REC=irec,IOSTAT=IOSTAT)
     '            (C_DATA(n)(1:C_LEN),n=1,NDATA)
                ELSE IF(IPTYPE.EQ.IPINTE) THEN             !integer data
                  WRITE(UNIT=IFIXED,FMT=G_FORMAT,REC=irec,IOSTAT=IOSTAT)
     '              (I_DATA(n),n=1,NDAT)
                ELSE IF(IPTYPE.EQ.IPLOGI) THEN             !logical data
                  WRITE(UNIT=IFIXED,FMT=G_FORMAT,REC=irec,IOSTAT=IOSTAT)
     '              (L_DATA(n),n=1,NDAT)
                ELSE IF(IPTYPE.EQ.IPREAL) THEN             !real data
                  WRITE(UNIT=IFIXED,FMT=G_FORMAT,REC=irec,IOSTAT=IOSTAT)
     '              (R_DATA(n),n=1,NDAT)
                ENDIF !iptype

              ELSE IF(N1char.GT.0) THEN !write input string
                CALL STRING_TRIM(CLINE,IBEG,IEND)
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                  WRITE(OP_STRING,'('' CLINE='',A)') CLINE(IBEG:IEND)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                ENDIF
                IF(POS_dollar+N2char+4.GT.LEN(OUTFORMAT)) THEN
                  CALL FLAG_ERROR(-1,
     &              'OUTFORMAT not long enough for format string')
                  ERROR=' '
                  GOTO 9999
                ENDIF
                OUTFORMAT=G_FORMAT(1:POS_dollar+N2char)//''',A)'
C KAT 25Jul99: Extending over two lines if longer than one.
                IEND1=132-MAX(N2char-2,0) !max length of CLINE
                IF(IEND.LE.IEND1) THEN
                  WRITE(UNIT=IFIXED,FMT=OUTFORMAT,
     '              REC=irec,IOSTAT=IOSTAT) CLINE(IBEG:IEND)
                ELSEIF(REVERT2) THEN
                  ! Making space for a broken line not implemented.
                  CALL WRITE_LINE(IOER,
     &              'WARNING: Input too long to update file correctly.'
     &              //'  Needs manual correction.',ERROR,*9999)
                  WRITE(UNIT=IFIXED,FMT=OUTFORMAT,
     '              REC=irec,IOSTAT=IOSTAT) CLINE(IBEG:IEND)
                ELSE !look for a place to break the line
                  IEND1=IEND1+1
                  DO WHILE(CLINE(IEND1:IEND1).NE.' '
     '              .AND.CLINE(IEND1:IEND1).NE.','.AND.IEND1.GT.0)
                    IEND1=IEND1-1
                  ENDDO
                  IF(IEND1.GT.0) THEN !write some on this line
                    WRITE(UNIT=IFIXED,FMT=OUTFORMAT,
     '                REC=irec,IOSTAT=IOSTAT) CLINE(IBEG:IEND1-1)
                  ELSE !start new line
                    WRITE(UNIT=IFIXED,FMT=OUTFORMAT,
     '                REC=irec,IOSTAT=IOSTAT) ' '
                  ENDIF
                  ! The record may have advanced by more than one line.
                  ! INQUIRE seems more reliable than No_lines_in_format.
                  INQUIRE(UNIT=IFIXED,OPENED=OPENED,NEXTREC=irec)
                  WRITE(UNIT=IFIXED,FMT='(1X,A)',
     '              REC=irec,IOSTAT=IOSTAT) CLINE(IEND1+1:IEND)
                ENDIF

                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                  call mp_setlock()
                  WRITE(OP_STRING,'('' IOSTAT='',I3)') IOSTAT
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                  call mp_unsetlock()
                ENDIF
C                ENDIF !revert2

              ENDIF !N1char
            ENDIF !fileip

            IF(REVERT2) THEN
              MODIFYFILE(IFIXED)=.TRUE.
              WRITE(OP_STRING,'('' >>Input file updated'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF !REVERT2

            IF(IOSTAT.EQ.0) THEN
              ERROR=' '
            ELSE
              ERROR=' Integer data write error'
              GOTO 9999
            ENDIF
          ELSE
            ERROR=' '
          ENDIF !opened
C KAT 2002-06-04: backup doesn't work anyway: removed buffers
C          IF(IBUFFER.GT.1) THEN
C            IBUFFER=IBUFFER-1
C            GO TO 10
C          ENDIF
        ENDIF !write prompted input to file

C     Read input from file
      ELSE IF(IO_TYPE.EQ.2.OR.IO_TYPE.EQ.4) THEN !read ip from file
        I_DATA(0)=0 !initialise
        I_DATA(1)=0 !initialise

        INQUIRE(UNIT=IFIXED,OPENED=OPENED,NEXTREC=irec)
        IF(OPENED) THEN
          FIRST_MISMATCH=.TRUE.
          FIRST_LOOP=.TRUE.
          FOUND_BLANK=.FALSE.
          MISMATCH_MESSAGE=.FALSE.
          IF(SKIP_LINE) THEN
            MULTICHOICE=No_lines_in_format.GT.1
          ELSE
            MULTICHOICE=No_lines_in_format.GT.0
          ENDIF
          LINE_CHECK_BEG=2
          FORM_CHECK_BEG=POS_quote1+LINE_CHECK_BEG
C KAT 9Dec99: Trimming spaces off expected input.
          DO WHILE(FORM_END.GT.FORM_CHECK_BEG
     '      .AND.G_FORMAT(FORM_END:FORM_END).EQ.' ')
            FORM_END=FORM_END-1
          ENDDO
C         Locate left square brackets in G_FORMAT & LINE
          POS_leftsqb1=
     '      CLOCAT('[',G_FORMAT(FORM_CHECK_BEG:POS_quote2-2))
          IF(POS_leftsqb1.NE.0) THEN
            LINE_CHECK_END=LINE_CHECK_BEG-1+POS_leftsqb1
            FORM_CHECK_END=FORM_CHECK_BEG-1+POS_leftsqb1
          ELSE
            IF(MULTICHOICE) THEN
              LINE_CHECK_END=10
              FORM_CHECK_END=FORM_CHECK_BEG
     '          +LINE_CHECK_END-LINE_CHECK_BEG
            ELSE
              FORM_CHECK_END=FORM_END
              LINE_CHECK_END=LINE_CHECK_BEG
     '          +FORM_CHECK_END-FORM_CHECK_BEG
            ENDIF
          ENDIF
 3        READ(UNIT=IFIXED,FMT='(A)',REC=irec,IOSTAT=IOSTAT) LINE
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' $$$Read line irec='',I4,'' LINE:'',A)')
     '        irec,LINE(1:20)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF

          IF(IOSTAT.EQ.0) THEN

C KAT 21May99: moved outside loop.
C            LINE_CHECK_BEG=2

            IF(SKIP_LINE.AND..NOT.FOUND_BLANK) THEN !skip line
              FOUND_BLANK=
     '          LINE(LINE_CHECK_BEG:10).EQ.BLANK(LINE_CHECK_BEG:10)
              IF(FOUND_BLANK) THEN !blank line has been found
                irec=irec+1
                GO TO 3 !to read next line
C             Else check question anyway.
              ENDIF
            ENDIF

            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' N1char='',I2,'' N2char='',I4,'
     '          //''' irec='',I6)') N1char,N2char,irec
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,
     '          '('' LINE_CHECK_BEG='',I4,'' LINE_CHECK_END='',I4,'
     '          //'/'' FORM_CHECK_BEG='',I4,'' FORM_CHECK_END='',I4)')
     '          LINE_CHECK_BEG,LINE_CHECK_END,
     '          FORM_CHECK_BEG,FORM_CHECK_END
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' ######LINE='',A,''##'')')
     '          LINE(LINE_CHECK_BEG-1:LINE_CHECK_END)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' ##G_FORMAT='',A,''##'')')
     '          G_FORMAT(FORM_CHECK_BEG-1:FORM_CHECK_END)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
            ENDIF !dop

C           Check that format string matches file string
            IF(   LINE(LINE_CHECK_BEG:LINE_CHECK_END).NE.
     '        G_FORMAT(FORM_CHECK_BEG:FORM_CHECK_END)) THEN
              IF(.NOT.MULTICHOICE.OR.FIRST_LOOP) THEN
C                IF(FIRST_MISMATCH) THEN
                IF(.NOT.MISMATCH_MESSAGE) THEN !no mismatch message yet
                  BLANKLINE=.TRUE.
C KAT 17Oct98: Check comment character is the first on the line and
C              check the first two characters as well.
C                  COMMENTLINE=.FALSE.
C                  DO j=LINE_CHECK_BEG,LINE_CHECK_END
C                    IF(LINE(j:j).NE.' ') BLANKLINE=.FALSE.
C                    IF(LINE(j:j).EQ.'!') COMMENTLINE=.TRUE.
C                  ENDDO
C                  IF((.NOT.BLANKLINE).AND.(.NOT.COMMENTLINE)) THEN
                  j=0
                  DO WHILE(BLANKLINE.AND.j.LT.LINE_CHECK_END)
                    j=j+1
                    BLANKLINE=LINE(j:j).EQ.' '
                  ENDDO
                  IF(.NOT.BLANKLINE) THEN
                    IF(LINE(j:j).NE.'!') THEN !not a comment line
C                      FIRST_MISMATCH=.FALSE.
CC$                    call mp_setlock()
C                      WRITE(OP_STRING,'('' >>Format out of place:'
C     '                  //' you may want to rewrite the file'')')
C                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
CC$                    call mp_unsetlock()
                      MISMATCH_MESSAGE=.TRUE.
                    ENDIF !not comment
                  ENDIF !not blankline
C KAT 17Oct98: Setting the first mismatch position even for blank lines
C              and comment lines ensures that any new questions are
C              inserted before the blank or comment which probably is
C              related to the next question.
                  IF(FIRST_MISMATCH) THEN
                    FIRST_MISMATCH=.FALSE.
                    IF(FOUND_BLANK) THEN
                      irec_at_first_mismatch=irec-1
                    ELSE
                      irec_at_first_mismatch=irec
                    ENDIF
                    IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                    call mp_setlock()
                      WRITE(OP_STRING,'('' irec_at_first_mismatch='','
     '                  //'I4)') irec_at_first_mismatch
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                    call mp_unsetlock()
                    ENDIF
                  ENDIF
                ENDIF !no mismatch message
                IF(SKIP_LINE) FOUND_BLANK=
     '            LINE(LINE_CHECK_BEG:10).EQ.BLANK(LINE_CHECK_BEG:10)
                irec=irec+1
                GO TO 3 !to read next line in search for match

              ELSE !IF(.NOT.FIRST_LOOP)  THEN
                IF(FIRST_MISMATCH) THEN
                  FIRST_MISMATCH=.FALSE.
                  irec_at_first_mismatch=irec
                  IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                    call mp_setlock()
                    WRITE(OP_STRING,'('' irec_at_first_mismatch='','
     '                //'I4)') irec_at_first_mismatch
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                    call mp_unsetlock()
                  ENDIF
                ENDIF
                WRITE(OP_STRING,'('' >>No match found: revert to '
     '            //'prompt'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                REVERT2=.TRUE.
                GO TO 5 !to use prompted input for current question
              ENDIF !not MULTICHOICE & FIRST_LOOP

            ELSE !check whether first character in LINE is ! or *
              CALL CHECK_IP_CHAR(irec,N1char,N2char,POS_dollar,
     '          G_FORMAT,LINE,REVERT1,ERROR,*3,*5,*9999)
            ENDIF !line.NE.g_format

C KAT17Oct98: Moved this here to avoid printing the message if the file
C             is updated correctly anyway.  The mesage is now only
C             printed if the format is present but out of place.  The
C             user will be prompted if the format is not present.
            IF(MISMATCH_MESSAGE.AND.FIRST_LOOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' >>Format out of place:'
     '          //' you may want to rewrite the file'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF

            IF(MULTICHOICE.AND.FIRST_LOOP) THEN !check lines of MCQ
              FIRST_LOOP=.FALSE.
C             Record position of start of question
C KAT 13Oct98: new SKIP_LINE adjustment
C              irec_at_line_1_MCQ=irec
              IF(FOUND_BLANK) THEN !skipped blank line
                irec_at_line_1_MCQ=irec-1 !include blank line
              ELSE
                irec_at_line_1_MCQ=irec
              ENDIF
              No_lines=1
              LAST_LINE=.FALSE.
              DO WHILE(.NOT.LAST_LINE) !to loop thru lines of MCQ
                irec=irec+1
                READ(UNIT=IFIXED,FMT='(A)',REC=irec,IOSTAT=IOSTAT)
     '            LINE
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                  call mp_setlock()
                  WRITE(OP_STRING,'('' LINE:'',A)') LINE(1:6)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                  call mp_unsetlock()
                ENDIF
                IF(LINE(1:4).EQ.'   ('.OR.
     '            LINE(1:3).EQ.'  ('.OR.
     '            LINE(1:4).EQ.'  *('.OR.
     '            LINE(1:3).EQ.' *(') THEN
                  No_lines=No_lines+1
                ELSE
                  LAST_LINE=.TRUE.
                  No_lines_in_MCQ=No_lines
                  IF(SKIP_LINE) No_lines_in_MCQ=No_lines_in_MCQ+1
                  irec_last_line_MCQ=irec !end position in file
                ENDIF
              ENDDO !while.not.last_line
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                WRITE(OP_STRING(1),'('' No_lines_in_MCQ   ='',I3)')
     '            No_lines_in_MCQ
                WRITE(OP_STRING(2),'('' No_lines_in_format='',I3)')
     '            No_lines_in_format
                WRITE(OP_STRING(3),'('' irec_at_line_1_MCQ='',I3)')
     '            irec_at_line_1_MCQ
                WRITE(OP_STRING(4),'('' irec_last_line_MCQ='',I3)')
     '            irec_last_line_MCQ
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
              ENDIF
              IF(No_lines_in_MCQ.NE.No_lines_in_format) THEN
                WRITE(OP_STRING,'('' >>Incorrect #options in: '',A)')
     '            G_FORMAT(FORM_CHECK_BEG:FORM_CHECK_END)
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C                 Find end of file
                IOSTAT=0
                DO WHILE(IOSTAT.EQ.0)
                  irec=irec+1
                  READ(UNIT=IFIXED,REC=irec,FMT='(A)',IOSTAT=IOSTAT)
     '              LINE
                ENDDO
                irec_at_EOF=irec
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                  call mp_setlock()
                  WRITE(OP_STRING,'('' irec_at_EOF='',I4)')irec_at_EOF
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                  call mp_unsetlock()
                ENDIF
C               Shift lines down to remove incomplete question
                irec_shift=irec_last_line_MCQ-irec_at_line_1_MCQ+1
                irec_at_EOF=irec_at_EOF-irec_shift
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                  call mp_setlock()
                  WRITE(OP_STRING(1),
     '              '('' $$$Shift lines up to remove MCQ'')')
                  WRITE(OP_STRING(2),'('' irec_shift ='',I4)')
     '              irec_shift
                  WRITE(OP_STRING(3),'('' irec_at_EOF='',I4)')
     '              irec_at_EOF
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                  call mp_unsetlock()
                ENDIF
                DO irec=irec_at_line_1_MCQ,irec_at_EOF
                  READ(UNIT=IFIXED,FMT='(A)',REC=irec+irec_shift,
     '              IOSTAT=IOSTAT) LINE
                  WRITE(UNIT=IFIXED,FMT='(A)',REC=irec,IOSTAT=IOSTAT)
     '              LINE
                  IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                  call mp_setlock()
                    WRITE(OP_STRING,'('' read from irec+irec_shift='''
     '                //',I3,'' write to irec='',I3,'' LINE: '',A)')
     '                irec+irec_shift,irec,LINE(1:15)
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                  call mp_unsetlock()
                  ENDIF
                ENDDO !irec
                REVERT2=.TRUE.
                GO TO 5 !to use prompted input for current question
              ELSE
                CALL CHECK_IP_CHAR(irec_last_line_MCQ,N1char,N2char,
     '            POS_dollar,G_FORMAT,LINE,REVERT1,ERROR,*3,*5,*9999)
              ENDIF !lines don't match
            ENDIF

C KAT 18Aug00:
C Note 9.23 in X3J3/96-007 draft standard says "An end-of-file condition
C may occur only during execution of a sequential input file".  We
C should therefore get a +ve (error) IOSTAT if the record didn't exist.
C However SGI's f77 gives a -1 to indicate end-of-file.  It would be
C nice to check that we don't have another error condition but
C error codes differ between compilers.  Therefore we assume that
C a non-zero IOSTAT means end-of-file.
          ELSE ! IOSTAT /= 0
C          ELSE IF(((OS_TYPE(1:3).EQ.'VMS').AND.(IOSTAT.EQ.36)).OR.
C     '        ((OS_TYPE(1:4).EQ.'UNIX').AND.(IOSTAT.EQ.-1))) THEN
            IF(FIRST_MISMATCH) THEN !hit EOF before file search
              FIRST_MISMATCH=.FALSE.
              IF(FOUND_BLANK) THEN
                irec_at_first_mismatch=irec-1
              ELSE
                irec_at_first_mismatch=irec
              ENDIF
            ENDIF
            irec_at_EOF=irec
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' irec_at_EOF='',I3)') irec_at_EOF
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
            WRITE(OP_STRING,'('' >>No match found: revert to '
     '        //'prompt'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            REVERT2=.TRUE.
            GO TO 5 !to use prompted input for current question
C          ELSE
C            CALL FLAG_ERROR(IOSTAT,'File read error')
C            ERROR=' '
C            GOTO 9999
          ENDIF !iostat

          IF(MULTICHOICE.AND.FIRST_LOOP) THEN
            ERROR=' Error with multichoice Q ?'
            GO TO 9999
          ELSE
            LINE_QUES_END=FORM_END-POS_2quote
C new MPN 22Apr97: allow for varying length defaults in single line
C                  questions by adjusting N2char
            IF(.NOT.MULTICHOICE.AND.POS_leftsqb1.NE.0) THEN
!Determine the length of the default value in LINE
              len_default_LINE=CLOCAT(']',LINE(LINE_CHECK_END:))-2
!Locate index of left square bracket in G_FORMAT
C                POS_leftsqb2=CLOCAT('[',G_FORMAT(1:))
!Determine the length of the default value in G_FORMAT
C                len_default_FORM=CLOCAT(']',G_FORMAT(POS_leftsqb2:))-2
              len_default_FORM=CLOCAT(']',G_FORMAT(FORM_CHECK_END:))-2
C KAT 9Dec99: Trimming spaces off expected input.
              LINE_QUES_END=LINE_QUES_END
     '          +len_default_LINE-len_default_FORM
C              N2char=N2char+len_default_LINE-len_default_FORM
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                WRITE(OP_STRING,
     '            '( '' POS_leftsqb1='',I4,'' len_default_LINE='',I4,'
     '            //''' POS_leftsqb2='',I4,'' len_default_FORM='',I4)')
     '            POS_leftsqb1,len_default_LINE,
     '            FORM_CHECK_END,len_default_FORM
C     '              POS_leftsqb2,len_default_FORM
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
              ENDIF !DOP
            ENDIF !not MULTICHOICE and POS_leftsqb1.NE.0
C end new
C KAT 9Dec99: Trimming spaces off expected input.
C            LINE_QUES_END=N2char-N1char+1
            CALL STRING_TRIM(LINE(LINE_QUES_END+1:),IBEG1,IEND1)
            IEND1=LINE_QUES_END+IEND1
            IBEG1=LINE_QUES_END+IBEG1
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' CHAR='',A)') LINE(IBEG1:IEND1)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF

            IF(LINE(IBEG1:IEND1).EQ.'FILE'
     '        .OR.LINE(IBEG1:IEND1).EQ.'file') THEN
C              CHAR1(1)(IBEG1:IEND1)=LINE(IBEG1:IEND1)
C              IF(ABBRV(CHAR1(1),'FILE',2)) THEN
C              IF(CBBREV(CHAR1(1),'FILE',2,1,1,N3CO)) THEN
              IF(.NOT.FILEIP) THEN
                FILEIP=.TRUE.
                READ(UNIT=IFIXED,FMT='(17X,A)',REC=irec+1,
     '            IOSTAT=IOSTAT) FILE07
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                  WRITE(OP_STRING,*) ' FILE07=',FILE07
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                ENDIF
                CALL STRING_TRIM(FILE07,IBEG,IEND)
                CALL OPENF(7,'DISK',FILE07(IBEG:IEND)//'.data','OLD',
     '            'SEQUEN','FORMATTED',132,ERROR,*9999)
              ENDIF
              IF(IPTYPE.EQ.IPANSW) THEN !yes/no data
                READ(7,*) (A_DATA(n),n=1,NDATA)
              ELSE IF(IPTYPE.EQ.IPCHAR) THEN !character data
                READ(7,*) (C_DATA(n),n=1,NDATA)
              ELSE IF(IPTYPE.EQ.IPINTE) THEN !integer data
                READ(7,*) (I_DATA(n),n=1,NDATA)
              ELSE IF(IPTYPE.EQ.IPLOGI) THEN !logical data
                READ(7,*) (L_DATA(n),n=1,NDATA)
              ELSE IF(IPTYPE.EQ.IPREAL) THEN !real data
                READ(7,*) (R_DATA(n),n=1,NDATA)
              ENDIF

            ELSE !no additional file input. This is normal file input
              NDATA_RET=NDATA
              IF(IPTYPE.EQ.IPANSW) THEN !yes/no data
C                  A_DATA(1)=CUPPER(LINE(IBEG1:IEND1))
                CALL CUPPER(LINE(IBEG1:IEND1),A_DATA(1))
              ELSE IF(IPTYPE.EQ.IPCHAR) THEN !character data
                CALL PARSSL(LINE(IBEG1:IEND1),NDATA,NDATA_RET,C_DATA,
     '            ERROR,*9999)
                DO nd=NDATA_RET+1,NDATA
                  C_DATA(nd)=C_DEFLT(nd)
                ENDDO !nd
              ELSE IF(IPTYPE.EQ.IPINTE) THEN !integer data
C KAT 12May99:  To avoid large unadjustable static arrays, NQLIST should
C               be passed to this routine instead of IDATA.
                IF(C_DATA(1)(1:5).EQ.'GRIDS') THEN
C                 Maximum number of data is size of NQLIST
                  IOIM=NQM
                ELSE
C                 Maximum number of data is size of IDATA
                  IOIM=IOIMX
                ENDIF
                CALL PARSILG(I_DATA,IOIM,C_DATA(1),LINE(IBEG1:IEND1),
     '            ERROR,*9999)
C KAT 25Jul99:  If some but not enough data is supplied, try the next
C               line for the rest of the data values.
C               Error handling in this needs to be improved.
                IF(I_DATA(0).NE.0.AND.I_DATA(0).LT.NDATA) THEN
                  irec=irec+1
                  READ(UNIT=IFIXED,FMT='(A)',REC=irec,IOSTAT=IOSTAT)
     '              LINE
                  NDATA_RET=I_DATA(0)
                  IDATA_TEMP=I_DATA(NDATA_RET)
                  CALL PARSILG(I_DATA(NDATA_RET),IOIM-NDATA_RET,
     '              C_DATA(1),LINE,ERROR,*9999)
                  I_DATA(0)=NDATA_RET+I_DATA(NDATA_RET)
                  I_DATA(NDATA_RET)=IDATA_TEMP
                ENDIF
                ! DPN 14/09/98 - check bounds
                DO n=1,I_DATA(0)
C KAT 20Jan00: Output range min if violated
                  IF(I_DATA(n).LT.IMIN) THEN
                    IEND=0
                    CALL APPENDC(IEND,'Data l.t. min. ',ERROR)
                    CALL APPENDI(IEND,IMIN,ERROR)
                    GOTO 9999
                  ELSE IF(I_DATA(n).GT.IMAX) THEN
                    IEND=0
                    CALL APPENDC(IEND,'Data g.t. max. ',ERROR)
                    CALL APPENDI(IEND,IMAX,ERROR)
                    GOTO 9999
C                  IF(I_DATA(n).LT.IMIN.OR.I_DATA(n).GT.IMAX) THEN
CC LKC 21-DEC-1998 output max and mins
CC                    ERROR='Integer data range'
C                    WRITE(ERROR,'(''Integer data range (IMAX: ''
C     '                ,I8,'')'')') IMAX
C                    GOTO 9999
                  ENDIF
                ENDDO
                DO nd=I_DATA(0)+1,NDATA !use defaults for unsupplied data
                  I_DATA(nd)=I_DEFLT(nd)
                ENDDO !nd
              ELSE IF(IPTYPE.EQ.IPLOGI) THEN !logical data
                L_DATA(1)=LFROMC(LINE(IBEG1:IEND1))
              ELSE IF(IPTYPE.EQ.IPREAL) THEN !real data
                CALL PARSRL(LINE(IBEG1:IEND1),IORMX,NDATA_RET,R_DATA,
     '            ERROR,*9999)
                ! DPN 14/09/98 - check bounds
                DO n=1,NDATA_RET
                  IF(R_DATA(n).LT.RMIN.OR.R_DATA(n).GT.RMAX) THEN
                    ERROR='Real data range'
                    GOTO 9999
                  ENDIF
                ENDDO
                DO nd=NDATA_RET+1,NDATA
                  R_DATA(nd)=R_DEFLT(nd)
                ENDDO !nd
              ENDIF
              IF(NDATA_RET.GT.NDATA) THEN
                ERROR='Too many data items'
                GOTO 9999
              ENDIF
            ENDIF !file input
          ENDIF !MULTICHOICE and FIRST_LOOP

C          ENDIF !n1char
          IF(IOSTAT.EQ.0) THEN
            ERROR=' '
          ELSE IF((OS_TYPE(1:3).EQ.'VMS'.AND.IOSTAT.EQ.36).OR.
     '        (OS_TYPE(1:4).EQ.'UNIX'.AND.IOSTAT.EQ.-1)) THEN
            ERROR='End of file'
            GO TO 9999
          ELSE
            CALL FLAG_ERROR(IOSTAT,'File read error')
            ERROR=' '
            GOTO 9999
          ENDIF !iostat

        ELSE !not opened
          WRITE(UNIT=CIUNIT,FMT='(I10)') IFIXED
          CALL STRING_TRIM(CIUNIT,IBEG,IEND)
          ERROR=' Unit '//CIUNIT(IBEG:IEND)//' is not open'
          GOTO 9999
        ENDIF !opened

        IF(IO_TYPE.EQ.4) THEN !list out to screen as read from file
          IF(IPTYPE.EQ.IPMESS) THEN
            WRITE(OP_STRING,G_FORMAT)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE IF(IPTYPE.EQ.IPANSW) THEN
            WRITE(OP_STRING,G_FORMAT) (A_DATA(n),n=1,NDATA)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE IF(IPTYPE.EQ.IPCHAR) THEN
            CALL STRING_TRIM(C_DATA(1),IBEG,IEND)
            WRITE(OP_STRING,G_FORMAT) (C_DATA(n)(IBEG:IEND),n=1,
     '        NDATA)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE IF(IPTYPE.EQ.IPINTE) THEN
            WRITE(OP_STRING,G_FORMAT) (I_DATA(n),n=1,NDATA)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE IF(IPTYPE.EQ.IPLOGI) THEN
            WRITE(OP_STRING,G_FORMAT) (L_DATA(n),n=1,NDATA)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE IF(IPTYPE.EQ.IPREAL) THEN
            WRITE(OP_STRING,G_FORMAT) (R_DATA(n),n=1,NDATA)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF !io_type=4 (listed read)

C     Write data to file
      ELSE IF(IO_TYPE.EQ.3) THEN !write data to file
        INQUIRE(UNIT=IFIXED,OPENED=OPENED,NEXTREC=irec)
        IF(OPENED.AND.IFIXED.NE.0) THEN
          IF(IPTYPE.EQ.IPMESS) THEN
            WRITE(UNIT=IFIXED,FMT=G_FORMAT,REC=irec,IOSTAT=IOSTAT)
          ELSE IF(IPTYPE.EQ.IPANSW) THEN
            WRITE(UNIT=IFIXED,FMT=G_FORMAT,REC=irec,IOSTAT=IOSTAT)
     '        (A_DATA(n),n=1,NDATA)
          ELSE IF(IPTYPE.EQ.IPCHAR) THEN
            WRITE(UNIT=IFIXED,FMT=G_FORMAT,REC=irec,IOSTAT=IOSTAT)
     '        (C_DATA(n)(1:C_LEN),n=1,NDATA)
          ELSE IF(IPTYPE.EQ.IPINTE) THEN
            WRITE(UNIT=IFIXED,FMT=G_FORMAT,REC=irec,IOSTAT=IOSTAT)
     '        (I_DATA(n),n=1,NDATA)
          ELSE IF(IPTYPE.EQ.IPLOGI) THEN
            WRITE(UNIT=IFIXED,FMT=G_FORMAT,REC=irec,IOSTAT=IOSTAT)
     '        (L_DATA(n),n=1,NDATA)
          ELSE IF(IPTYPE.EQ.IPREAL) THEN
            WRITE(UNIT=IFIXED,FMT=G_FORMAT,REC=irec,IOSTAT=IOSTAT)
     '        (R_DATA(n),n=1,NDATA)
          ENDIF
          IF(IOSTAT.EQ.0) THEN
            ERROR=' '
          ELSE
            ERROR=' Integer data write error'
            GOTO 9999
          ENDIF
        ELSE
          ERROR=' '
        ENDIF

!     Set data to defaults
      ELSE IF(IO_TYPE.EQ.5) THEN !set to defaults
        DO nd=1,NDATA
          IF(IPTYPE.EQ.IPANSW) THEN
            A_DATA(nd)=A_DEFLT(nd)
          ELSE IF(IPTYPE.EQ.IPCHAR) THEN
            C_DATA(nd)=C_DEFLT(nd)
          ELSE IF(IPTYPE.EQ.IPINTE) THEN
            I_DATA(nd)=I_DEFLT(nd)
          ELSE IF(IPTYPE.EQ.IPLOGI) THEN
            L_DATA(nd)=L_DEFLT(nd)
          ELSE IF(IPTYPE.EQ.IPREAL) THEN
            R_DATA(nd)=R_DEFLT(nd)
          ENDIF
        ENDDO
      ENDIF !io_type

      CALL EXITS(ROUTINENAME)
      RETURN

 9999 CALL STRING_TRIM(ERROR,IBEG,IEND)
      IF(ERROR(1:IEND).EQ.'End of file') THEN
        IO_TYPE=1
        WRITE(OP_STRING,'('' End of file: revert to interactive '','
     '    //'''input'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' G_FORMAT='',A)') G_FORMAT(1:80)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        GO TO 10

C KAT 21Jan99: This is not implemented
C      ELSE IF(ERROR(1:5).EQ.'Apply') THEN
C        IO_TYPE=5
C        WRITE(OP_STRING,'('' Apply previous value to all '
C     '    //'subsequent'')')
C        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C        APPLY=.TRUE.
C        GO TO 10

      ELSE IF(ERROR(1:IEND).EQ.'Default') THEN
        IO_TYPE=5
        WRITE(OP_STRING,'('' Revert to default input'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        GOTO 10

      ELSE IF(ERROR(1:IEND).EQ.'Edit') THEN
        ERROR=' '
        CALL CLOSEF(IFIXED,ERROR,*9999)
        ERROR='Edit>GINPUT>GINOUT'

      ELSE IF(ERROR(1:IEND).EQ.'Restart') THEN
        ERROR(IEND+1:)='>GINPUT>GINOUT'
        CLOSE(UNIT=IFIXED,STATUS='DELETE')

      ELSE IF(ERROR(1:14).EQ.'File not found') THEN
        ERROR='Help file not found>GINOUT'

      ELSE
        IF(ERROR.NE.' ') THEN!FLAG_ERROR has not been called
          CALL FLAG_ERROR(0,ERROR(:LEN_TRIM(ERROR)))
          ERROR=' '
        ENDIF

C 21-DEC-1998 This doesn't actually give any useful information
C  and ends up pushing usefull info off the end of the string
C  Removing the IOSTAT and irec and increase the format string
C
C        WRITE(UNIT=CIOSTA,FMT='(I5)') IOSTAT
C        WRITE(UNIT=CIREC ,FMT='(I5)') irec
C        ERROR=ERROR(1:IEND)//': IOSTAT='//CIOSTA(1:5)//' irec='//
C     '    CIREC(1:5)//' G_FORMAT='//G_FORMAT(1:10)//'>GINOUT'

        CALL WRITE_CHAR(IOER,'  G_FORMAT: ',ERR)
        CALL WRITE_STRING(IOER,LEN_TRIM(G_FORMAT),G_FORMAT,ERR)
        CALL WRITE_STRING(IOER,1,NEWLINE,ERR)

        CALL ERRORIN(ROUTINENAME)
        CLOSE(UNIT=IFIXED)
      ENDIF

      CALL EXITS(ROUTINENAME)
      RETURN 1
      END


