      SUBROUTINE LILEAD(STRING,ERROR,*)

C#### Subroutine: LILEAD
C###  Description:
C###    LILEAD lists electrocardiographic leads.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND
      CHARACTER FILE*100
      LOGICAL OPFILE

      CALL ENTERS('LILEAD',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list lead<;FILENAME>
C###  Description:
C###    Lists electrocardiographic leads to screen or to
C###    FILENAME with extension *.oplead if filename specified.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LILEAD',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.oplead','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        CALL OPLEAD(ERROR,*9999)

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF

      ENDIF

      CALL EXITS('LILEAD')
      RETURN
 9999 CALL ERRORS('LILEAD',ERROR)
      CALL EXITS('LILEAD')
      RETURN 1
      END


