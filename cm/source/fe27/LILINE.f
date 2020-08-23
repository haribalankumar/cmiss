      SUBROUTINE LILINE(NEL,NLLINE,NPL,DL,STRING,ERROR,*)

C#### Subroutine: LILINE
C###  Description:
C###    LILINE lists line parameters.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NEL(0:NELM,NLM),NLLINE(0:NL_R_M,0:NRM),NPL(5,0:3,NLM)
      REAL*8 DL(3,NLM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,N3CO
      CHARACTER FILE*100
      LOGICAL CBBREV,OPFILE

      CALL ENTERS('LILINE',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list lines<;FILENAME>
C###  Description:
C###    Lists line parameters to screen or file FILENAME.opline if
C###    qualifier present.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C#### Command: FEM list lines<;FILENAME> groups
C###  Description:
C###    Lists line groups to screen or to FILENAME.opline if qualifier
C###    present.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> groups'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LILINE',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opline','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        IF(CBBREV(CO,'GROUP',1,noco+1,NTCO,N3CO)) THEN
          CALL OPLINEG(ERROR,*9999)
        ELSE
          CALL OPLINE(NEL,NLLINE,NPL,DL,ERROR,*9999)
        ENDIF

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LILINE')
      RETURN
 9999 CALL ERRORS('LILINE',ERROR)
      CALL EXITS('LILINE')
      RETURN 1
      END


