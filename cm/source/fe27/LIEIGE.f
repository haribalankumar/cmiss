      SUBROUTINE LIEIGE(EIGV,YP,STRING,FIX,ERROR,*)

C#### Subroutine: LIEIGE
C###  Description:
C###    LEIGE lists eigenvalues and eigenvectors.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
!     Parameter List
      REAL*8 EIGV(NOM,NTM,2),YP(NYM,NIYM,NXM)
      CHARACTER ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER IBEG,IEND,nr,nx
      CHARACTER FILE*100
      LOGICAL OPFILE

      CALL ENTERS('LIEIGE',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list eigenvalues<;FILENAME>
C###  Description:
C###    Lists eigenvalues and eigenvectors to the screen or to
C###    FILENAME.opeige if qualifier present.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIEIGE',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opeige','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        nx=1 !Temporary
        nr=1 !Temporary
        CALL ASSERT(ITYP2(nr,nx).GT.0,
     '    '>>no problem type defined',ERROR,*9999)

        CALL OPEIGE(EIGV,FIX(1,1,nx),YP,ERROR,*9999)

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LIEIGE')
      RETURN
 9999 CALL ERRORS('LIEIGE',ERROR)
      CALL EXITS('LIEIGE')
      RETURN 1
      END


