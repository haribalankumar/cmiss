      SUBROUTINE LIPLIN(STRING,ERROR,*)

C#### Subroutine: LIPLIN
C###  Description:
C###    LIPLIN lists polyline parameters.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,N3CO
      CHARACTER FILE*100
      LOGICAL CBBREV,OPFILE

      CALL ENTERS('LIPLIN',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list polyline<;FILENAME> <groups>
C###  Description:
C###    List polyline parameters to screen or FILENAME.opplin if
C###    qualifier present.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> <groups>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIPLIN',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opplin','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        IF(CBBREV(CO,'GROUP',1,noco+1,NTCO,N3CO)) THEN
          CALL OPPLING(ERROR,*9999)
        ELSE
          CALL OPPLIN(ERROR,*9999)
        ENDIF

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LIPLIN')
      RETURN
 9999 CALL ERRORS('LIPLIN',ERROR)
      CALL EXITS('LIPLIN')
      RETURN 1
      END


