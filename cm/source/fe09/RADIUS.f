      REAL*8 FUNCTION RADIUS(ORDER,BRANCH)

C#### FUNCTION: RADIUS
C#### Type: REAL*8
C###  Description:
C###    RADIUS calculates the radius of a segment from
C###    it's order

      IMPLICIT NONE

!     Parameter List
      INTEGER ORDER,BRANCH
!     Local Variables

C NB orders are reversed at the momnet and there is no stscastic
C component, or recognition of detintion between LAD,CX,RCA

c      SCALAR=1.5d0

      IF (BRANCH.EQ.1) THEN  !RCA
        IF (ORDER.EQ.11)THEN
          RADIUS=1.61d0
        ELSE IF(ORDER.EQ.10)THEN
          RADIUS=0.65d0
        ELSE IF(ORDER.EQ.9)THEN
          RADIUS=0.35d0
        ELSE IF(ORDER.EQ.8)THEN
          RADIUS=0.21d0
        ELSE IF(ORDER.EQ.7)THEN
          RADIUS=0.13d0
        ELSE IF(ORDER.EQ.6)THEN
          RADIUS=0.07d0
        ENDIF
      ELSEIF (BRANCH.EQ.2) THEN !LAD
         IF (ORDER.EQ.11)THEN
          RADIUS=1.61d0
        ELSE IF(ORDER.EQ.10)THEN
          RADIUS=0.65d0
        ELSE IF(ORDER.EQ.9)THEN
          RADIUS=0.35d0
        ELSE IF(ORDER.EQ.8)THEN
          RADIUS=0.21d0
        ELSE IF(ORDER.EQ.7)THEN
          RADIUS=0.13d0
        ELSE IF(ORDER.EQ.6)THEN
          RADIUS=0.07d0
        ENDIF
      ELSEIF (BRANCH.EQ.3) THEN !CX
        IF (ORDER.EQ.11)THEN
          RADIUS=0.0d0
        ELSE IF(ORDER.EQ.10)THEN
          RADIUS=1.26d0
        ELSE IF(ORDER.EQ.9)THEN
          RADIUS=0.48d0
        ELSE IF(ORDER.EQ.8)THEN
          RADIUS=0.23d0
        ELSE IF(ORDER.EQ.7)THEN
          RADIUS=0.14d0
        ELSE IF(ORDER.EQ.6)THEN
          RADIUS=0.07d0
        ENDIF
      ENDIF
      RETURN
      END


