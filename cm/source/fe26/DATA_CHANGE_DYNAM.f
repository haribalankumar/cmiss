      SUBROUTINE DATA_CHANGE_DYNAM(CMGUI_DATA,
     '  nd,CHANGE_STRUCTURE,ERROR,*)

C#### Subroutine: DATA_CHANGE_DYNAM
C###  Description:
C###    Currently used to notify CMGUI when the values
C###    or structure of a data point change (called by
C###    DATA_CHANGE).

      IMPLICIT NONE
      INCLUDE 'cmgui00.cmn'
!     Parameter List
      INTEGER CMGUI_DATA(CMGUI_DATA_LEN),nd
      LOGICAL CHANGE_STRUCTURE
      CHARACTER ERROR*(*)
!     Local Variables

      CALL ENTERS('DATA_CHANGE_DYNAM',*9999)

C     Set the flag in the array
      DATA_CHANGED=.TRUE.
      IF(CMGUI_DATA(nd).EQ.0) THEN
        CMGUI_DATA(nd)=1
      ELSEIF(CMGUI_DATA(nd).EQ.2) THEN
        CMGUI_DATA(nd)=3
      ENDIF

      CALL EXITS('DATA_CHANGE_DYNAM')
      RETURN
 9999 CALL ERRORS('DATA_CHANGE_DYNAM',ERROR)
      CALL EXITS('DATA_CHANGE_DYNAM')
      RETURN 1
      END


