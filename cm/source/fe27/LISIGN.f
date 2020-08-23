      SUBROUTINE LISIGN(LD,NBJ,NDDATA,WD,XID,ZD,STRING,ERROR,*)

C#### Subroutine: LISIGN
C###  Description:
C###    LISIGN lists information about a signal file

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(MXCH)
      INTEGER LD(NDM),NBJ(NJM,NEM),NDDATA(0:NDM,0:NRM)
      REAL*8 WD(NJM,NDM),XID(NIM,NDM),ZD(NJM,NDM)
!     Local Variables
      INTEGER IBEG,IEND,N3CO
      CHARACTER FILE*100,SIGFNAME*100,FILEFORMAT*6
      LOGICAL CBBREV,OPFILE

      CALL ENTERS('LISIGN',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list signal<;FILENAME>
C###  Parameter:           <signal FILENAME[$current]>
C###    Specifies the signal filename to be listed to screen
C###  Parameter:           <ascii/binary [ascii]>
C###    Specify whether the file is stored as a binary or
C###    ascii file.
C###  Parameter:           <region #[1]>
C###    Specify the element file region numbers to be defined.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped.
C###  Description:
C###   Give general information about a signal file
C###

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<signal FILNAME[$current]>'
        OP_STRING(3)=BLANK(1:15)//'<ascii/binary [ascii]>'
        OP_STRING(4)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LISIGN',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opsign','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        IF(CBBREV(CO,'SIGNAL',3,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          SIGFNAME=CO(N3CO+1)(IBEG:IEND)
        ELSE
          CALL STRING_TRIM(FILE00,IBEG,IEND)
          SIGFNAME=FILE00(IBEG:IEND)
        ENDIF

        IF(CBBREV(CO,'BINARY',3,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          FILEFORMAT='BINARY'
        ELSE
          FILEFORMAT='ASCII'
        ENDIF

        CALL OPSIGN(LD,NBJ,NDDATA,WD,XID,ZD,
     '    FILEFORMAT,SIGFNAME,ERROR,*9999)

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LISIGN')
      RETURN
 9999 CALL ERRORS('LISIGN',ERROR)
      CALL EXITS('LISIGN')
      RETURN 1
      END


