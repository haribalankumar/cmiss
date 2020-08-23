      SUBROUTINE LICOLO(STRING,ERROR,*)

C#### Subroutine: LICOLO
C###  Description:
C###    LICOLO lists workstation index colours.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw(2),N3CO,NTIL
      LOGICAL CBBREV

      CALL ENTERS('LICOLO',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list colours<;FILENAME>
C###  Parameter:    <on WS#[1]>
C###    Specify the workstation (GX window) to draw the
C###    points on.
C###  Description:
C###    Lists the index colours of the specified workstation. Lists
C###    to the screen or to FILENAME.opcolo if qualifier present.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<on WS#[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe14','doc','LICOLO',ERROR,*9999)
      ELSE
        IF(CBBREV(CO,'ON',1,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),1,NTIL,iw,ERROR,*9999)
        ELSE
          iw(1)=1
        ENDIF

c       CALL OPCOLO(iw,ERROR,*9999)
      ENDIF

      CALL EXITS('LICOLO')
      RETURN
 9999 CALL ERRORS('LICOLO',ERROR)
      CALL EXITS('LICOLO')
      RETURN 1
      END


