      SUBROUTINE PLANE_FROM_3_PTS(NORML,NORMALTYPE,
     '  POINT1,POINT2,POINT3,ERROR,*)

C#### Subroutine: PLANE_FROM_3_PTS
C###  Description:
C###    PLANE_FROM_3_PTS finds the equation of a plane in three
C###    dimensions and a vector normal to the plane from three
C###    non-colinear points.
C###    NORMALTYPE=1 for raw normal and plane equation
C###    NORMALTYPE=2 for unit normal and plane equation
C###    The coefficients represent aX + bY + cZ + d = 0
C###    NORML(1)=a,NORML(2)=b,NORML(3)=c,NORML(4)=d
C***  Created by Martin Buist, Jan 1997

      IMPLICIT NONE

      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter list
      INTEGER NORMALTYPE
      REAL*8 POINT1(3),POINT2(3),POINT3(3),NORML(4)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER nj
      REAL*8 DIFF1(3),DIFF2(3),NORMSIZE
      LOGICAL COLINEAR

      CALL ENTERS('PLANE_FROM_3_PTS',*9999)

      CALL ASSERT(NJT.EQ.3,'>>ERROR, must have 3d for plane+normal'
     '  ,ERROR,*9999)

C Check for colinearity
      COLINEAR=.FALSE.
      CALL CHECK_COLINEAR(POINT1,POINT2,POINT3,COLINEAR,ERROR,*9999)
      IF(.NOT.COLINEAR) THEN
        DO nj=1,NJT
          DIFF1(nj)=POINT2(nj)-POINT1(nj)
          DIFF2(nj)=POINT2(nj)-POINT3(nj)
        ENDDO !nj

        NORML(1)=(DIFF1(2)*DIFF2(3))-(DIFF1(3)*DIFF2(2))
        NORML(2)=(DIFF1(3)*DIFF2(1))-(DIFF1(1)*DIFF2(3))
        NORML(3)=(DIFF1(1)*DIFF2(2))-(DIFF1(2)*DIFF2(1))

        IF(NORMALTYPE.EQ.2) THEN
          NORMSIZE=0.0d0
          DO nj=1,NJT
            NORMSIZE=NORMSIZE+(NORML(nj)**2.0d0)
          ENDDO !nj
          NORMSIZE=DSQRT(NORMSIZE)
          DO nj=1,NJT
            NORML(nj)=NORML(nj)/NORMSIZE
          ENDDO !nj
        ENDIF

        NORML(4)=0.0d0
        DO nj=1,NJT
          NORML(4)=NORML(4)-(NORML(nj)*POINT1(nj))
        ENDDO !nj

      ELSE !Colinear

        WRITE(OP_STRING,'('' COLINEAR points in PLANE_FROM_3_PTS '')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' No plane generated '')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nj=1,4
          NORML(nj)=0.0d0
        ENDDO !nj
      ENDIF

      CALL EXITS('PLANE_FROM_3_PTS')
      RETURN
 9999 CALL ERRORS('PLANE_FROM_3_PTS',ERROR)
      CALL EXITS('PLANE_FROM_3_PTS')
      RETURN 1
      END


