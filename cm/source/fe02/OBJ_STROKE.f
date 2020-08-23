      SUBROUTINE OBJ_STROKE(ERROR,*)

C#### Subroutine: OBJ_STROKE
C###  Description:
C###    OBJ_STROKE defines object from GKS stroke created by label
C###    command.

      IMPLICIT NONE
      INCLUDE 'obje00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables

      CALL ENTERS('OBJ_STROKE',*9999)

      NTOBJE=NTOBJE+1
      OBJECT_NAME(NTOBJE)='STROKE'
      NSOBJE(2,NTOBJE)=1 !is number of parts to object

      CALL EXITS('OBJ_STROKE')
      RETURN
 9999 CALL ERRORS('OBJ_STROKE',ERROR)
      CALL EXITS('OBJ_STROKE')
      RETURN 1
      END


