      SUBROUTINE LICORN(STRING,ERROR,*)

C#### Subroutine: LICORN
C###  Description:
C###    LICORN lists corner nodes (for BE problems).

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,N3CO,nr
      CHARACTER FILE*100
      LOGICAL CBBREV,OPFILE

      CALL ENTERS('LICORN',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list corners<;FILENAME>
C###  Description:
C###    Lists corner nodes. Lists to the screen or to FILENAME.opnods
C###    if the qualifier is present.
C###  Parameter:    <region #[1]>
C###    Specify the region number to list.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LICORN',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opcorn','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        IF(CBBREV(CO,'REGION',2,noco+1,NTCO,N3CO)) THEN
          nr=IFROMC(CO(N3CO+1))
          IF(nr.GT.NRT) NRT=nr
        ELSE
          nr=1
        ENDIF

! Needs rewriting
!        nc=1 !Temporary
!        CALL OPCORN(NPNODE,nr,NVHP,ERROR,*9999)

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LICORN')
      RETURN
 9999 CALL ERRORS('LICORN',ERROR)
      CALL EXITS('LICORN')
      RETURN 1
      END


