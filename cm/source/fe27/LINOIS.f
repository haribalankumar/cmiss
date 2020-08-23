      SUBROUTINE LINOIS(STRING,ERROR,*)

C#### Subroutine: LINOIS
C###  Description:
C###    LINOIS lists noise parameters
C***  Created By   Leo Cheng 23-OCT-97

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'

!     Parameter List
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND

      CALL ENTERS('LINOIS',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C-------------------------------------------------------------------

C#### Command: FEM list noise
C###  Parameter:
C###  Description:
C###    Lists to noise parameters to screen.
C###

        OP_STRING(1)=STRING(1:IEND)
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C-------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LICOOR',ERROR,*9999)
      ELSE

C***    LKC 23-OCT-97 There is no file output as OPNOIS files are
C***                  created by the command FEM apply noise;file
C***                  may be useful latter
C        IF(NTCOQU(noco).GT.0) THEN !file output
C          OPFILE=.TRUE.
C          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
C          CALL STRING_TRIM(FILE,IBEG,IEND)
C          IOFI=IOFILE1
C          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opcoor','NEW',
C     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
C        ELSE
C          OPFILE=.FALSE.
C        ENDIF
C        OPFILE=.FALSE.

        CALL  OPNOIS(ERROR,*9999)
C        IF(OPFILE) THEN
C          CALL CLOSEF(IOFI,ERROR,*9999)
C          IOFI=IOOP
C        ENDIF
      ENDIF

      CALL EXITS('LINOIS')
      RETURN
 9999 CALL ERRORS('LINOIS',ERROR)
      CALL EXITS('LINOIS')
      RETURN 1
      END


