      SUBROUTINE LIINVE(STRING,ERROR,*)

C#### Subroutine: LIINVE
C###  Description:
C###    LIINVE lists regularisation options for the inverse of the
C###    transfer matrix.

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

      CALL ENTERS('LIINVE',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list inverse<;FILENAME>
C###  Description:
C###    Lists regularisation options for the inverse of the transfer
C###    matrix to screen or to FILENAME with extension *.opinve
C###    if filename is specified.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIINVE',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.optran','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
          IOFI=IOOP
        ENDIF

        CALL OPINVE(ERROR,*9999)

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LIINVE')
      RETURN
 9999 CALL ERRORS('LIINVE',ERROR)
      CALL EXITS('LIINVE')
      RETURN 1
      END


