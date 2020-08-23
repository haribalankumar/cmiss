      SUBROUTINE LIWIND(STRING,ERROR,*)

C#### Subroutine: LIWIND
C###  Description:
C###    LIWIND lists window parameters.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND
      LOGICAL ABBREV,FULL

      CALL ENTERS('LIWIND',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list window
C###  Description:
C###    List workstation window coordinate range to the screen.
C###  Parameter:    <full>
C###    Specify the addition listing of Window type,Open status
C###    and Graphics type.

        OP_STRING(1)=STRING(1:IEND)//' <full>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIWIND',ERROR,*9999)
      ELSE
        IF(ABBREV(CO(noco+1),'FULL',1)) THEN
          FULL=.TRUE.
        ELSE
          FULL=.FALSE.
        ENDIF

        CALL OPWIND(FULL,ERROR,*9999)

      ENDIF

      CALL EXITS('LIWIND')
      RETURN
 9999 CALL ERRORS('LIWIND',ERROR)
      CALL EXITS('LIWIND')
      RETURN 1
      END


