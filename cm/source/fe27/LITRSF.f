      SUBROUTINE LITRSF(NPLIST3,NPLIST4,NPLIST5,STRING,ERROR,*)

C#### Subroutine: LITRSF
C###  Description:
C###    LITRSF lists parameters for transfer matrix from
C###    epicardium to torso surface.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NPLIST3(0:NPM),NPLIST4(0:NPM),NPLIST5(0:NPM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND
      CHARACTER FILE*100
      LOGICAL OPFILE

      CALL ENTERS('LITRSF',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list transfer<;FILENAME>
C###  Description:
C###    Lists transfer parameters for transfer matrix from one surface
C###    to another to the screen or to FILENAME with extension
C###    *.optrst if filename specified.

C---------------------------------------------------------------------

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LITRSF',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.optrsf','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
          IOFI=IOOP
        ENDIF

CC AJPs 11-11-97 - 191297 - rgb
        CALL OPTRSF(NPLIST3,NPLIST4,NPLIST5,ERROR,*9999)
CC AJPe

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LITRSF')
      RETURN
 9999 CALL ERRORS('LITRSF',ERROR)
      CALL EXITS('LITRSF')
      RETURN 1
      END


