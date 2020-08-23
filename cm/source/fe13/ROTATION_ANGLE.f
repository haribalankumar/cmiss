      SUBROUTINE ROTATION_ANGLE(np1,np2,np3,np4,np5,ANGLE,XP,ERROR,
     &  *)

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter values
      INTEGER np1,np2,np3,np4,np5
      REAL*8 ANGLE,XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER nj
      REAL*8 CALC_ANGLE,XPOINT(3,5),norm_1(4),norm_2(4),temp
      
      CALL ENTERS('ROTATION_ANGLE',*9999)
      
      DO nj=1,3
        XPOINT(nj,1)=XP(1,1,nj,np1)
        XPOINT(nj,2)=XP(1,1,nj,np2)
        XPOINT(nj,3)=XP(1,1,nj,np3)
        XPOINT(nj,4)=XP(1,1,nj,np4)
        XPOINT(nj,5)=XP(1,1,nj,np5)
      ENDDO !nj
      CALL PLANE_FROM_3_PTS(norm_1,2,XPOINT(1,1),XPOINT(1,2),
     '  XPOINT(1,3),ERROR,*9999)
      CALL NORMALISE2(3,norm_1,temp,ERROR,*9999) !unit vector
      CALL PLANE_FROM_3_PTS(norm_2,2,XPOINT(1,2),XPOINT(1,4),
     '  XPOINT(1,5),ERROR,*9999)
      CALL NORMALISE2(3,norm_2,temp,ERROR,*9999) !unit vector
      ANGLE=CALC_ANGLE(norm_1,norm_2)

      CALL EXITS('ROTATION_ANGLE')
      RETURN
 9999 CALL ERRORS('ROTATION_ANGLE',ERROR)
      CALL EXITS('ROTATION_ANGLE')
      RETURN 1
      END
