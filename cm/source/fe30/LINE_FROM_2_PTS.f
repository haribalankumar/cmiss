      SUBROUTINE LINE_FROM_2_PTS(COEF,POINT1,POINT2,LINETYPE,ERROR,*)

C#### Subroutine: LINE_FROM_2_PTS
C###  Description:
C###    LINE_FROM_2_PTS finds the equation of a line in two
C###    dimensions from two points
C###    LINETYPE=1 if the equation is returned as y=mx+c
C###    LINETYPE=2 if the equation is returned as x=c
C###    COEF(1)=m,COEF(2)=c (or COEF(1)=c;LINETYPE=2)
C***  Created by Martin Buist, Jan 1997

      IMPLICIT NONE

      INCLUDE 'tol00.cmn'

!     Parameter list
      INTEGER LINETYPE
      REAL*8 POINT1(2),POINT2(2),COEF(2)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER nj

      CALL ENTERS('LINE_FROM_2_PTS',*9999)

C      CALL ASSERT(NJT.EQ.2,'>>ERROR, must have 2d for line eqtn'
C     '  ,ERROR,*9999)

      DO nj=1,2
        COEF(nj)=0.0d0
      ENDDO !nj

      IF(DABS(POINT2(1)-POINT1(1)).GT.ZERO_TOL) THEN
        LINETYPE=1
        COEF(1)=(POINT2(2)-POINT1(2))/(POINT2(1)-POINT1(1))
        COEF(2)=POINT1(2)-(COEF(1)*POINT1(1))
      ELSE
        LINETYPE=2
        COEF(1)=POINT1(1)
      ENDIF

      CALL EXITS('LINE_FROM_2_PTS')
      RETURN
 9999 CALL ERRORS('LINE_FROM_2_PTS',ERROR)
      CALL EXITS('LINE_FROM_2_PTS')
      RETURN 1
      END



