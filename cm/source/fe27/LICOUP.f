      SUBROUTINE LICOUP(NEELEM,NP_INTERFACE,NPNODE,STRING,ERROR,*)

C#### Subroutine: LICOUP
C###  Description:
C###    LICOUP lists coupling between different regions.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NEELEM(0:NE_R_M,0:NRM),NP_INTERFACE(0:NPM,0:3),
     '  NPNODE(0:NP_R_M,0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND
      CHARACTER FILE*100
      LOGICAL OPFILE

      CALL ENTERS('LICOUP',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list coupling<;FILENAME>
C###  Description:
C###    Lists coupling between different regions. Lists to the screen
C###    or to FILENAME.opcoup if the qualifier is present.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LICOUP',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opcoup','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        CALL OPCOUP(NEELEM,NP_INTERFACE,NPNODE,ERROR,*9999)

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LICOUP')
      RETURN
 9999 CALL ERRORS('LICOUP',ERROR)
      CALL EXITS('LICOUP')
      RETURN 1
      END


