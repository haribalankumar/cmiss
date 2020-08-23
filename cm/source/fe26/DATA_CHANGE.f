      SUBROUTINE DATA_CHANGE(nd,CHANGE_STRUCTURE,ERROR,*)

C#### Subroutine: DATA_CHANGE
C###  Description:
C###    Currently used to notify CMGUI when the values
C###    or structure of a data point change.  Structure is
C###    defined as the number of fields, number of
C###    components in each field, etc.

      IMPLICIT NONE
      INCLUDE 'cmgui00.cmn'
      INCLUDE 'mach00.inc'

!     Parameter List
      INTEGER nd
      LOGICAL CHANGE_STRUCTURE
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER*4 NEW_PTR

      CALL ENTERS('DATA_CHANGE',*9999)

      IF(nd.GT.CMGUI_DATA_LEN) THEN
C       Increase the size of the array
        NEW_PTR=0
        CALL ALLOCATE_MEMORY(2*nd,2*nd,INTTYPE,NEW_PTR,.FALSE.,
     '    ERROR,*9999)
C       Copy the existing data and zero the new data
        CALL COPY_AND_INIT_INT(2*nd,CMGUI_DATA_LEN,
     '    %VAL(NEW_PTR),%VAL(CMGUI_DATA_PTR),
     '    0,ERROR,*9999)
        IF(CMGUI_DATA_LEN.NE.0) THEN
          CALL FREE_MEMORY(CMGUI_DATA_PTR,ERROR,*9999)
        ENDIF
        CMGUI_DATA_LEN=2*nd
        CMGUI_DATA_PTR=NEW_PTR
      ENDIF
      CALL DATA_CHANGE_DYNAM(%VAL(CMGUI_DATA_PTR),
     '  nd,CHANGE_STRUCTURE,ERROR,*9999)

      CALL EXITS('DATA_CHANGE')
      RETURN
 9999 CALL ERRORS('DATA_CHANGE',ERROR)
      CALL EXITS('DATA_CHANGE')
      RETURN 1
      END


