      SUBROUTINE NODE_DESTROY_DYNAM(CMGUI_NODE,np,ERROR,*)

C#### Subroutine: NODE_DESTROY_DYNAM
C###  Description:
C###    Currently used to notify CMGUI when a node
C###    has been created (called by NODE_DESTROY).

      IMPLICIT NONE
      INCLUDE 'cmgui00.cmn'
!     Parameter List
      INTEGER CMGUI_NODE(CMGUI_NODE_LEN),np
      CHARACTER ERROR*(*)
!     Local Variables

      CALL ENTERS('NODE_DESTROY_DYNAM',*9999)

C     Set the flag in the array
      NODES_DELETED=.TRUE.
      IF(CMGUI_NODE(np).EQ.0) THEN
        CMGUI_NODE(np)=2
      ELSE IF(CMGUI_NODE(np).EQ.1) THEN
        CMGUI_NODE(np)=3
      ENDIF

      CALL EXITS('NODE_DESTROY_DYNAM')
      RETURN
 9999 CALL ERRORS('NODE_DESTROY_DYNAM',ERROR)
      CALL EXITS('NODE_DESTROY_DYNAM')
      RETURN 1
      END



