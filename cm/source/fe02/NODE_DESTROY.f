      SUBROUTINE NODE_DESTROY(np,ERROR,*)

C#### Subroutine: NODE_DESTROY
C###  Description:
C###    Currently used to notify CMGUI when a node
C###    has been created.

      IMPLICIT NONE
      INCLUDE 'cmgui00.cmn'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER np
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER*4 NEW_PTR

      CALL ENTERS('NODE_DESTROY',*9999)

      IF(np.GT.CMGUI_NODE_LEN) THEN
C       Increase the size of the array
        NEW_PTR=0
        CALL ALLOCATE_MEMORY(2*np,2*np,INTTYPE,NEW_PTR,.FALSE.,
     '    ERROR,*9999)
C       Copy the existing data and zero the new data
        CALL COPY_AND_INIT_INT(2*np,CMGUI_NODE_LEN,
     '    %VAL(NEW_PTR),%VAL(CMGUI_NODE_PTR),
     '    0,ERROR,*9999)
        IF(CMGUI_NODE_LEN.NE.0) THEN
          CALL FREE_MEMORY(CMGUI_NODE_PTR,ERROR,*9999)
        ENDIF
        CMGUI_NODE_LEN=2*np
        CMGUI_NODE_PTR=NEW_PTR
      ENDIF
      CALL NODE_DESTROY_DYNAM(%VAL(CMGUI_NODE_PTR),np,ERROR,*9999)

      CALL EXITS('NODE_DESTROY')
      RETURN
 9999 CALL ERRORS('NODE_DESTROY',ERROR)
      CALL EXITS('NODE_DESTROY')
      RETURN 1
      END


