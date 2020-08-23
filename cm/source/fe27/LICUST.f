      SUBROUTINE LICUST(NRLIST,STRING,ERROR,*)

C#### Subroutine: LICUST
C###  Description:
C###    LICUST lists customisation coefficients for deforming a mesh.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NRLIST(0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,nolist
C SMAR009 19/01/99 removed ,nr
      CHARACTER FILE*100
      LOGICAL ALL_REGIONS,OPFILE

      CALL ENTERS('LICUST',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list customisation<;FILENAME>
C###  Description:
C###    Lists customisation coefficients for deforming a mesh
C###    to screen or file FILENAME.opcust if
C###    qualifier present in the directory
C###    specified by PATH with $current specifing the current default file..
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LICUST',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opcust','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
C SMAR009 does this do loop need to be here ?
        DO nolist=1,NRLIST(0)
C SMAR009 18/01/99 not used nr=NRLIST(nolist)
          CALL OPCUST(ERROR,*9999)
C SMAR009 18/01/99 removed nr,
        ENDDO

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LICUST')
      RETURN
 9999 CALL ERRORS('LICUST',ERROR)
      CALL EXITS('LICUST')
      RETURN 1
      END


