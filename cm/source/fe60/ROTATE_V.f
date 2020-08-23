      SUBROUTINE ROTATE_V(ANGLE,V,V_ROTATE,ERROR,*)

C#### Subroutine: ROTATE_V 
C###  Description:
C###    Outputs the transform matrix to rotate points about a vector V.


      IMPLICIT NONE
!     Parameter list
      REAL*8 ANGLE,V(3),V_ROTATE(3,3)
      CHARACTER ERROR*(*)

      CALL ENTERS('ROTATE_V',*9999)

      V_ROTATE(1,1)=DCOS(ANGLE)+(1.d0-DCOS(ANGLE))*V(1)**2.d0
      V_ROTATE(1,2)=(1.d0-DCOS(ANGLE))*V(1)*V(2)-(V(3)*DSIN(ANGLE))
      V_ROTATE(1,3)=(1.d0-DCOS(ANGLE))*V(1)*V(3)+V(2)*DSIN(ANGLE)
      V_ROTATE(2,1)=(1.d0-DCOS(ANGLE))*V(1)*V(2)+V(3)*DSIN(ANGLE)
      V_ROTATE(2,2)=DCOS(ANGLE)+(1.d0-DCOS(ANGLE))*V(2)**2.d0
      V_ROTATE(2,3)=(1.d0-DCOS(ANGLE))*V(2)*V(3)-V(1)*DSIN(ANGLE)
      V_ROTATE(3,1)=(1.d0-DCOS(ANGLE))*V(1)*V(3)-V(2)*DSIN(ANGLE)
      V_ROTATE(3,2)=(1.d0-DCOS(ANGLE))*V(2)*V(3)+V(1)*DSIN(ANGLE)
      V_ROTATE(3,3)=DCOS(ANGLE)+(1.d0-DCOS(ANGLE))*V(3)**2.d0

      CALL EXITS('ROTATE_V')
      RETURN
 9999 CALL ERRORS('ROTATE_V',ERROR)
      CALL EXITS('ROTATE_V')
      RETURN 1
      END

