      SUBROUTINE OBJ_ELLIPSE(ISEGM,ERROR,*)

C#### Subroutine: OBJ_ELLIPSE
C###  Description:
C###    OBJ_ELLIPSE defines object from GKS ellipse created by label
C###    command.

      IMPLICIT NONE
      INCLUDE 'obje00.cmn'
!     Parameter List
      INTEGER ISEGM
      CHARACTER ERROR*(*)
!     Local Variables

      CALL ENTERS('OBJ_ELLIPSE',*9999)

      NTOBJE=NTOBJE+1
      OBJECT_NAME(NTOBJE)='ELLIPSE'
      NSOBJE(1,NTOBJE)=ISEGM !is segment number of object
      NSOBJE(2,NTOBJE)=1     !is number of parts to object

      CALL EXITS('OBJ_ELLIPSE')
      RETURN
 9999 CALL ERRORS('OBJ_ELLIPSE',ERROR)
      CALL EXITS('OBJ_ELLIPSE')
      RETURN 1
      END


