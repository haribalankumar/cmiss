      SUBROUTINE DEHEAD(STRING,ERROR,*)

C#### Subroutine: DEHEAD
C###  Description:
C###    DEHEAD defines heading.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'head00.cmn'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER IBEG,IEND

      CALL ENTERS('DEHEAD',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM define heading
C###  Description:
C###    Defines a heading string for input files.
C###    This string is inserted into all subsequent input files.
C###  Parameter:      <HEADING>
C###    Specify the heading string.

        OP_STRING(1)=STRING(1:IEND)//' <HEADING>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEHEAD',ERROR,*9999)
      ELSE
        CALL STRING_TRIM(CO(noco+1),IBEG,IEND)
        IF(IEND.GT.IBEG) THEN
          HEADING=CO(noco+1)(IBEG:IEND)
        ELSE
          HEADING=' '
        ENDIF
      ENDIF

      CALL EXITS('DEHEAD')
      RETURN
 9999 CALL ERRORS('DEHEAD',ERROR)
      CALL EXITS('DEHEAD')
      RETURN 1
      END


