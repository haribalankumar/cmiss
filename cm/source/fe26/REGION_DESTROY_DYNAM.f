      SUBROUTINE REGION_DESTROY_DYNAM(CMGUI_REGION,nr,ERROR,*)

C#### Subroutine: REGION_DESTROY_DYNAM
C###  Description:
C###    Currently used to notify CMGUI when a region
C###    has been created (called by REGION_DESTROY).

      IMPLICIT NONE
      INCLUDE 'cmgui00.cmn'
!     Parameter List
      INTEGER CMGUI_REGION(CMGUI_REGION_LEN),nr
      CHARACTER ERROR*(*)
!     Local Variables

      CALL ENTERS('REGION_DESTROY_DYNAM',*9999)

C     Set the flag in the array
      REGIONS_DELETED=.TRUE.
      IF(CMGUI_REGION(nr).EQ.0) THEN
        CMGUI_REGION(nr)=2
      ELSEIF(CMGUI_REGION(nr).EQ.1) THEN
        CMGUI_REGION(nr)=3
      ENDIF

      CALL EXITS('REGION_DESTROY_DYNAM')
      RETURN
 9999 CALL ERRORS('REGION_DESTROY_DYNAM',ERROR)
      CALL EXITS('REGION_DESTROY_DYNAM')
      RETURN 1
      END
