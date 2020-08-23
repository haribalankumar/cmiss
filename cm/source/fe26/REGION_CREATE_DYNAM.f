      SUBROUTINE REGION_CREATE_DYNAM(CMGUI_REGION,nr,ERROR,*)

C#### Subroutine: REGION_CREATE_DYNAM
C###  Description:
C###    Currently used to notify CMGUI when a region
C###    has been created (called by REGION_CREATE).

      IMPLICIT NONE
      INCLUDE 'cmgui00.cmn'
!     Parameter List
      INTEGER CMGUI_REGION(CMGUI_REGION_LEN),nr
      CHARACTER ERROR*(*)
!     Local Variables

      CALL ENTERS('REGION_CREATE_DYNAM',*9999)

C     Set the flag in the array
      REGIONS_CHANGED=.TRUE.
      IF(CMGUI_REGION(nr).EQ.0) THEN
        CMGUI_REGION(nr)=1
      ELSEIF(CMGUI_REGION(nr).EQ.2) THEN
        CMGUI_REGION(nr)=3
      ENDIF

      CALL EXITS('REGION_CREATE_DYNAM')
      RETURN
 9999 CALL ERRORS('REGION_CREATE_DYNAM',ERROR)
      CALL EXITS('REGION_CREATE_DYNAM')
      RETURN 1
      END


