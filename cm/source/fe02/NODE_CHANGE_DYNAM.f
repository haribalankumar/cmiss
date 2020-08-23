      SUBROUTINE NODE_CHANGE_DYNAM(CMGUI_NODE,
     '  np,CHANGE_STRUCTURE,ERROR,*)

C#### Subroutine: NODE_CHANGE_DYNAM
C###  Description:
C###    Currently used to notify CMGUI when the values
C###    or structure of a node change (called by
C###    NODE_CHANGE).

      IMPLICIT NONE
      INCLUDE 'cmgui00.cmn'
!     Parameter List
      INTEGER CMGUI_NODE(CMGUI_NODE_LEN),np
      LOGICAL CHANGE_STRUCTURE
      CHARACTER ERROR*(*)
!     Local Variables

      CALL ENTERS('NODE_CHANGE_DYNAM',*9999)

C     Set the flag in the array
      NODES_CHANGED=.TRUE.
      IF(CMGUI_NODE(np).EQ.0) THEN
        CMGUI_NODE(np)=1
      ELSE IF(CMGUI_NODE(np).EQ.2) THEN
        CMGUI_NODE(np)=3
      ENDIF

      CALL EXITS('NODE_CHANGE_DYNAM')
      RETURN
 9999 CALL ERRORS('NODE_CHANGE_DYNAM',ERROR)
      CALL EXITS('NODE_CHANGE_DYNAM')
      RETURN 1
      END


