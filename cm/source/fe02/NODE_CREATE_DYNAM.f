      SUBROUTINE NODE_CREATE_DYNAM(CMGUI_NODE,np,ERROR,*)

C#### Subroutine: NODE_CREATE_DYNAM
C###  Description:
C###    Currently used to notify CMGUI when a node
C###    has been created (called by NODE_CREATE).

      IMPLICIT NONE
      INCLUDE 'cmgui00.cmn'
!     Parameter List
      INTEGER CMGUI_NODE(CMGUI_NODE_LEN),np
      CHARACTER ERROR*(*)
!     Local Variables

      CALL ENTERS('NODE_CREATE_DYNAM',*9999)

C     Set the flag in the array
      NODES_CHANGED=.TRUE.
      IF(CMGUI_NODE(np).EQ.0) THEN
        CMGUI_NODE(np)=1
      ELSE IF(CMGUI_NODE(np).EQ.2) THEN
        CMGUI_NODE(np)=3
      ENDIF

      CALL EXITS('NODE_CREATE_DYNAM')
      RETURN
 9999 CALL ERRORS('NODE_CREATE_DYNAM',ERROR)
      CALL EXITS('NODE_CREATE_DYNAM')
      RETURN 1
      END


