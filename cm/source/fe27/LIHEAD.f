      SUBROUTINE LIHEAD(STRING,ERROR,*)

C#### Subroutine: LIHEAD
C###  Description:

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND

      CALL ENTERS('LIHEAD',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list heading
C###  Description:
C###    List file heading to screen.

        OP_STRING(1)=STRING(1:IEND)
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIHEAD',ERROR,*9999)
      ELSE
        CALL OPHEAD(ERROR,*9999)
      ENDIF

      CALL EXITS('LIHEAD')
      RETURN
 9999 CALL ERRORS('LIHEAD',ERROR)
      CALL EXITS('LIHEAD')
      RETURN 1
      END


