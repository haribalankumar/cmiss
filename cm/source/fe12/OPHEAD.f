      SUBROUTINE OPHEAD(ERROR,*)

C#### Subroutine: OPHEAD
C###  Description:
C###    OPHEAD outputs heading.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'head00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IBEG,IEND

      CALL ENTERS('OPHEAD',*9999)
      CALL STRING_TRIM(HEADING,IBEG,IEND)
      WRITE(OP_STRING,'('' Heading: '',A)') HEADING(IBEG:IEND)
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      CALL EXITS('OPHEAD')
      RETURN
 9999 CALL ERRORS('OPHEAD',ERROR)
      CALL EXITS('OPHEAD')
      RETURN 1
      END


