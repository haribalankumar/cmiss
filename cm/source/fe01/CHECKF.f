      SUBROUTINE CHECKF(Qual_ID,noco,NTCOQU,CO,COQU,FILE,STRING,*)

C#### Subroutine: CHECKF
C###  Description:
C###    <HTML><PRE>
C###    If Qual_ID=1:
C###    Checks whether filename defined in command list CO. This is for:
C###         open/read com           or export nodes/elem
C###      or open/read com;file      or export nodes/elem;file
C###      or open/read com;file;path or export nodes/elem;file;path
C###      or open/read com;file;example  or export nodes/elem;file;example
C###    where path follows Unix convention
C###           e.g. path=subdir1/subdir2
C###             or path=/disk/usr/people/hunter/subdir..
C###    The path is interpreted for specific os in OPENF in FE10.
C###
C###    If Qual_ID=2 and NTCOQU(noco).ge.1:
C###    Checks whether filename defined at position 2 in command
C###    qualifier list COQU and, if so, puts this name into FILE and
C###    defines FILE00 also as this name. Otherwise FILE is obtained
C###    from FILE00, unless FILE00 empty, in which case ->error.
C###    This is define command with file (& maybe path) specified:
C###           e.g. define node;r
C###                define node;r;myfile
C###                define node;r;myfile;path
C###                define node;r;myfile;example
C###
C###    If Qual_ID=3 and NTCOQU(noco)=1:
C###    Checks whether filename defined at position 1 in command
C###    qualifier list COQU and, if so, puts this name into FILE.
C###    This is usual list/update/evaluate command with file specified.
C###    If the iterate command has been issued the filename
C###    is appended with the iteration number.
C###    </PRE></HTML>

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'docd00.cmn'
      INCLUDE 'file00.cmn'
!     Parameter List
      INTEGER noco,Qual_ID,NTCOQU(*)
      CHARACTER CO(*)*(*),COQU(25,*)*(*),FILE*(*),STRING*(MXCH)
!     Local Variables
      CHARACTER ERRCHAR
C DBs 25/3/01.  F77 does not allow functions in parameter statements
C      PARAMETER(ERRCHAR=ACHAR(7))
C DBe
      INTEGER IBEGF,IBEGS,IENDF,IENDS
C DBe
      CHARACTER ERROR*10,FILENAME*(MXCH),SUBPATH*(MXCH)
      LOGICAL PATHGIVEN
!     Functions
      INTEGER LEN_TRIM
      LOGICAL ABBREV

      CALL ENTERS('CHECKF',*9999)

C DBs 25/3/01.  F77 does not allow functions in parameter statements
      ERRCHAR=ACHAR(7)
      FILENAME(1:1)=ERRCHAR
      FILE(:1)=ERRCHAR
      PATHGIVEN=.FALSE.
C     GMH 8/3/96 Set to false by default
      DOCDIR=.FALSE.
C KAT 7/6/00: Changed so that Example and Doc must be in the directory
C     not the filename position.
      IF(Qual_ID.EQ.1) THEN
C ***   This is for open/read com;file and export node/elem etc

        IF(NTCOQU(noco).EQ.0) THEN
C         use default file in default directory
C         eg open com
          CALL STRING_TRIM(FILE00,IBEGF,IENDF)
          CALL ASSERT(IENDF.GE.IBEGF,'>>Default filename invalid',
     '      ERROR,*9999)
          FILENAME=FILE00(IBEGF:IENDF)

        ELSE IF(NTCOQU(noco).EQ.1) THEN
C        ELSE IF(NTCOQU(noco).EQ.1.AND.
C     '    .NOT.ABBREV(COQU(noco,1),'DOC',1).AND.
C     '    .NOT.ABBREV(COQU(noco,1),'EXAMPLE',1)) THEN
C         use specified file in default directory
C         e.g. open com;file
          FILE=COQU(noco,1)
          CALL STRING_TRIM(FILE,IBEGF,IENDF)
          IF(FILE(IBEGF:IENDF).NE.'cmiss') FILE00=FILE(IBEGF:IENDF)
          FILENAME=FILE(IBEGF:IENDF)

C        ELSE IF(NTCOQU(noco).EQ.1.AND.
C     '    (ABBREV(COQU(noco,1),'DOC',1).OR.
C     '    ABBREV(COQU(noco,1),'EXAMPLE',1))) THEN
C          IF(ABBREV(COQU(noco,1),'DOC',1)) THEN
C            WRITE(OP_STRING,'('' >>WARNING: doc is obsolescent. '
C     '        //'Change doc to example'')')
C            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C          ENDIF
CC         use default filename in documentation directory
CC         eg open com;;example
C          CALL STRING_TRIM(FILE00,IBEGF,IENDF)
C          CALL ASSERT(IENDF.GE.IBEGF,'>>Default filename invalid',
C     '      ERROR,*9999)
C          FILE=FILE00
C          DOCDIR=.TRUE.
C          IF(ABBREV(CO(noco),'COMMANDS',2)) THEN
C            Open_COM_file=.TRUE.
C          ELSE
C            Open_COM_file=.FALSE.
C          ENDIF

C KAT 1/5/00: SUBPATH stuff is unnecessary
        ELSE IF(NTCOQU(noco).EQ.2.AND.
     '    .NOT.ABBREV(COQU(noco,2),'DOC',3).AND.
     '    .NOT.ABBREV(COQU(noco,2),'EXAMPLE',3)) THEN
C         use specified filename in specified directory
C         e.g. open com;file;subdir1/subdir2/
          FILE=COQU(noco,1)
          IF(FILE.NE.'cmiss') FILE00=FILE
          CALL STRING_TRIM(FILE,IBEGF,IENDF)
          FILENAME=FILE(IBEGF:IENDF)
          CALL STRING_TRIM(COQU(noco,2),IBEGS,IENDS)
          SUBPATH=COQU(noco,2)(IBEGS:IENDS)
          PATHGIVEN=.TRUE.

        ELSE IF(NTCOQU(noco).EQ.2.AND.
     '    (ABBREV(COQU(noco,2),'DOC',3).OR.
     '     ABBREV(COQU(noco,2),'EXAMPLE',3))) THEN
          IF(ABBREV(COQU(noco,2),'DOC',1)) THEN
            WRITE(OP_STRING,'('' >>WARNING: doc is obsolescent. '
     '        //'Change doc to example'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
C         use specified filename in documentation directory
C         e.g. open com;file;example
          IF(COQU(noco,1).NE.' ') THEN
            FILE=COQU(noco,1)
            IF(FILE.NE.'cmiss') FILE00=FILE
          ELSE
            FILE=FILE00
          ENDIF
          CALL STRING_TRIM(FILE,IBEGF,IENDF)
          DOCDIR=.TRUE.
        ENDIF !ntcoqu/doc

      ELSE IF(Qual_ID.EQ.2.AND.NTCOQU(noco).GE.1) THEN
C ***   This is define command

        IF(NTCOQU(noco).EQ.1) THEN
C         use default file in default directory
C         eg define nodes;r
          CALL STRING_TRIM(FILE00,IBEGF,IENDF)
          CALL ASSERT(IENDF.GE.IBEGF,'>>Default filename invalid',
     '      ERROR,*9999)
          FILENAME=FILE00(IBEGF:IENDF)

        ELSE IF(NTCOQU(noco).EQ.2) THEN
C        ELSE IF(NTCOQU(noco).EQ.2.AND.
C     '    .NOT.ABBREV(COQU(noco,2),'DOC',1).AND.
C     '    .NOT.ABBREV(COQU(noco,2),'EXAMPLE',1)) THEN
C         use specified file in default directory
C         eg define nodes;r;file
          FILE=COQU(noco,2)
          CALL STRING_TRIM(FILE,IBEGF,IENDF)
          FILE00=FILE(IBEGF:IENDF) !sets default filename
          FILENAME=FILE(IBEGF:IENDF)

C        ELSE IF(NTCOQU(noco).EQ.2.AND.
C     '    (ABBREV(COQU(noco,2),'DOC',1).OR.
C     '     ABBREV(COQU(noco,2),'EXAMPLE',1))) THEN
C          IF(ABBREV(COQU(noco,2),'DOC',1)) THEN
C            WRITE(OP_STRING,'('' >>WARNING: doc is obsolescent. '
C     '        //'Change doc to example'')')
C            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C          ENDIF
CC         use default filename in documentation directory
CC         eg define nodes;r;;example
C          CALL STRING_TRIM(FILE00,IBEGF,IENDF)
C          CALL ASSERT(IENDF.GE.IBEGF,'>>Default filename invalid',
C     '      ERROR,*9999)
C          FILE=FILE00
C          DOCDIR=.TRUE.
C          Open_COM_file=.FALSE.

C KAT 1/5/00: SUBPATH stuff is unnecessary
        ELSE IF(NTCOQU(noco).EQ.3.AND.
     '    .NOT.ABBREV(COQU(noco,3),'DOC',3).AND.
     '    .NOT.ABBREV(COQU(noco,3),'EXAMPLE',3)) THEN
C         use specified filename in specified directory
C         eg define nodes;r;file;subdir1/subdir2/..
C         or define nodes;r;file;/dir/subdir/..
          FILE=COQU(noco,2)
          CALL STRING_TRIM(FILE,IBEGF,IENDF)
          FILE00=FILE(IBEGF:IENDF) !sets default filename
          FILENAME=FILE(IBEGF:IENDF)
          CALL STRING_TRIM(COQU(noco,3),IBEGS,IENDS)
          SUBPATH=COQU(noco,3)(IBEGS:IENDS)
          PATHGIVEN=.TRUE.

        ELSE IF(NTCOQU(noco).EQ.3.AND.
     '    (ABBREV(COQU(noco,3),'DOC',1).OR.
     '     ABBREV(COQU(noco,3),'EXAMPLE',1))) THEN
          IF(ABBREV(COQU(noco,3),'DOC',1)) THEN
            WRITE(OP_STRING,'('' >>WARNING: doc is obsolescent. '
     '        //'Change doc to example'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
C         use specified filename in documentation directory
C         eg define nodes;r;file;example
          IF(COQU(noco,2).NE.' ') THEN
            FILE=COQU(noco,2)
            IF(FILE.NE.'cmiss') FILE00=FILE !sets default filename
          ELSE
            FILE=FILE00
          ENDIF
          CALL STRING_TRIM(FILE,IBEGF,IENDF)
          DOCDIR=.TRUE.

        ENDIF !ntcoqu

      ELSE IF(Qual_ID.EQ.3) THEN
C ***   This is usual list/update/evaluate command with file specified

        IF(NTCOQU(noco).EQ.0) THEN
C         list to output terminal
C         eg list nodes

        ELSE IF(NTCOQU(noco).EQ.1) THEN
C         list to specified file
C         eg list nodes;file
          FILE=COQU(noco,1)
C new MPN 15Jul97: this may be useful at some stage
C          CALL STRING_TRIM(FILE,IBEG,IEND)
C          CALL STRING_TRIM(PATH00,IBEG1,IEND1)
C          CALL ASSERT(IEND1.GE.IBEG1,'>>Default directory invalid',
C     '      ERROR,*9999)
C          TEMPFILE=PATH00(IBEG1:IEND1)//FILE(IBEG:IEND)
C          FILE=TEMPFILE
C end new
C KAT 26Oct00: Not used
C          IF(ITERATE_COMMAND_EXECUTE) THEN !append iteration#
C            CALL STRING_TRIM(FILE,IBEGF,IENDF)
C             WRITE(CHAR4,'(I4)') iter_counter
C            CALL STRING_TRIM(CHAR4,IBEG4,IEND4)
C            FILENAME=FILE(IBEGF:IENDF)//'_'//CHAR4(IBEG4:IEND4)
C            FILE=FILENAME
C          ENDIF !iterate_command_execute
C new MPN 15Jul97: this may be useful at some stage
C        ELSE IF(NTCOQU(noco).EQ.2) THEN
CC         list to specified file in another directory
CC         eg list nodes;file;another_dir
C          CALL STRING_TRIM(COQU(noco,1),IBEG,IEND)
C          CALL STRING_TRIM(COQU(noco,2),IBEG2,IEND2)
C          TEMPFILE=COQU(noco,2)(IBEG2:IEND2)//COQU(noco,1)(IBEG:IEND)
C          FILE=TEMPFILE
C        ELSE
C          CO(noco+1)='?'
C          STRING='>>Reenter: '//CO(noco)
C          GO TO 9999
C end new
        ENDIF !ntcoqu

C      ELSE
C        CO(noco+1)='?'
C        STRING='>>Reenter: '//CO(noco)
C        GO TO 9999
      ENDIF !Qual_ID

      IF(.NOT.DOCDIR.AND.FILENAME(1:1).NE.ERRCHAR) THEN
C KAT 1/5/00: SUBPATH stuff is unnecessary
CC KAT 27Mar98: Allowing paths without ending '/',and subdirectories.
CC     Build up FILE from PATH00 and SUBPATH and FILENAME
C      IF(.NOT.DOCDIR.AND.Qual_ID.NE.3) THEN
C        CALL STRING_TRIM(FILENAME,IBEGF,IENDF) !is this necessary?
C        IF(PATHGIVEN) THEN
C          CALL STRING_TRIM(SUBPATH,IBEGS,IENDS)
C          ABSOLUTE=SUBPATH(IBEGS:IBEGS).EQ.'/' !absolute path
C        ELSE
C          ABSOLUTE=.FALSE.
C        ENDIF
C        IF(ABSOLUTE) THEN !absolute path is supplied
C          FILE=SUBPATH(IBEGS:IENDS)
C        ELSE              !use PATH00
C          CALL STRING_TRIM(PATH00,IBEGP,IENDP)
C          CALL ASSERT(IENDP.GE.IBEGP,'>>Default directory invalid',
C     '      ERROR,*9999)
C          FILE=PATH00(IBEGP:IENDP)
C          IF(PATHGIVEN) THEN
C            CALL STRING_TRIM(FILE,IBEGP,IENDP)
C            IF(FILE(IENDP:IENDP).EQ.' ') THEN
C              FILE=SUBPATH(IBEGS:IENDS)
C            ELSE IF(FILE(IENDP:IENDP).EQ.'/') THEN
C              FILE=FILE(:IENDP)//SUBPATH(IBEGS:IENDS)
C            ELSE
C              FILE=FILE(:IENDP)//'/'//SUBPATH(IBEGS:IENDS)
C            ENDIF
C          ENDIF
C        ENDIF
C        CALL STRING_TRIM(FILE,IBEGP,IENDP)
C        IF(FILE(IENDP:IENDP).EQ.' ') THEN
C          FILE=FILENAME(IBEGF:IENDF)
C        ELSE IF(FILE(IENDP:IENDP).EQ.'/') THEN
C          FILE=FILE(:IENDP)//FILENAME(IBEGF:IENDF)
C        ELSE
C          FILE=FILE(:IENDP)//'/'//FILENAME(IBEGF:IENDF)
C        ENDIF
        FILE=PATH00
C KAT 1/5/00: SUBPATH stuff is unnecessary
        IF(PATHGIVEN.AND.SUBPATH.NE.' ') THEN
          IF(SUBPATH(1:1).EQ.'/') THEN
            IBEGS=1
          ELSE
            IBEGS=LEN_TRIM(FILE)+1
          ENDIF
          FILE(IBEGS:)=SUBPATH
        ENDIF
        IF(FILENAME(1:1).EQ.'/') THEN
          IBEGF=1
        ELSE
          IBEGF=LEN_TRIM(FILE)+1
        ENDIF
        FILE(IBEGF:)=FILENAME
      ENDIF

      IF(FILE(:1).EQ.ERRCHAR) THEN
        GO TO 9999
      ENDIF

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        CALL STRING_TRIM(FILE,IBEGF,IENDF)
        WRITE(OP_STRING,'('' FILE='',A)') FILE(IBEGF:IENDF)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' DOCDIR='',L1)') DOCDIR
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('CHECKF')
      RETURN

 9999 CO(noco+1)='?'
      STRING='>>Reenter: '//CO(noco)
      CALL EXITS('CHECKF')
      RETURN 1
      END


