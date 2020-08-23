      SUBROUTINE MESH_ANGLE(ANGLE,XP1,XP2,XP3,ERROR,*)

C#### Subroutine: MESH_ANGLE
C###  Description:
C###    MESH_ANGLE calculates the angle between vectors XP1-XP2 and XP2-XP3

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'
      REAL*8 ANGLE,XP1(3),XP2(3),XP3(3)
      CHARACTER ERROR*(*)
      !Local variables
      INTEGER nj
      REAL*8 SCALAR,U(3),V(3)


      CALL ENTERS('MESH_ANGLE',*9999)

      DO nj=1,3
        U(nj)=XP2(nj)-XP1(nj) !end - start
        V(nj)=XP3(nj)-XP2(nj) !end - start
      ENDDO !nj
      CALL NORMALISE(3,U,ERROR,*9999)
      CALL NORMALISE(3,V,ERROR,*9999)
      ANGLE=SCALAR(3,U,V)
      ANGLE=MAX(-1.d0,ANGLE)
      ANGLE=MIN(1.d0,ANGLE)
      ANGLE=DACOS(ANGLE)
      
      RETURN
 9999 CALL ERRORS('MESH_ANGLE',ERROR)
      CALL EXITS('MESH_ANGLE')
      RETURN 1
      END



