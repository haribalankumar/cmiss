      SUBROUTINE REGION_DESTROY(nr,ERROR,*)

C#### Subroutine: REGION_DESTROY
C###  Description:
C###    Currently used to notify CMGUI when a region
C###    has been destroyed.

      IMPLICIT NONE
      INCLUDE 'cmgui00.cmn'
      INCLUDE 'mach00.inc'

!     Parameter List
      INTEGER nr
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER*4 NEW_PTR

      CALL ENTERS('REGION_DESTROY',*9999)

      IF(nr.GT.CMGUI_REGION_LEN) THEN
C       Increase the size of the array
        NEW_PTR=0
        CALL ALLOCATE_MEMORY(2*nr,2*nr,INTTYPE,NEW_PTR,.FALSE.,
     '    ERROR,*9999)
C       Copy the existing data and zero the new data
        CALL COPY_AND_INIT_INT(2*nr,CMGUI_REGION_LEN,
     '    %VAL(NEW_PTR),%VAL(CMGUI_REGION_PTR),
     '    0,ERROR,*9999)
        IF(CMGUI_REGION_LEN.NE.0) THEN
          CALL FREE_MEMORY(CMGUI_REGION_PTR,ERROR,*9999)
        ENDIF
        CMGUI_REGION_LEN=2*nr
        CMGUI_REGION_PTR=NEW_PTR
      ENDIF
      CALL REGION_DESTROY_DYNAM(%VAL(CMGUI_REGION_PTR),nr,ERROR,*9999)

      CALL EXITS('REGION_DESTROY')
      RETURN
 9999 CALL ERRORS('REGION_DESTROY',ERROR)
      CALL EXITS('REGION_DESTROY')
      RETURN 1
      END


