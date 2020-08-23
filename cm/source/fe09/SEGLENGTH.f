      REAL*8 FUNCTION SEGLENGTH(ORDER,BRANCH)

C#### FUNCTION: SEGLENGTH
C#### Type: REAL*8
C###  Description:
C###    SEGLENGTH calculates the length of a segment from
C###    it's order

      IMPLICIT NONE

!     Parameter List
      INTEGER ORDER,BRANCH
      REAL*8 SCALAR
!     Local Variables

C NB orders are reversed at the momnet and there is no stscastic
C component, or recognition of detintion between LAD,CX,RCA

      SCALAR=1.4d0
      IF (BRANCH.EQ.1) THEN !RCA
        IF (ORDER.EQ.11)THEN
          SEGLENGTH=3.12d0*SCALAR
        ELSE IF(ORDER.EQ.10)THEN
          SEGLENGTH=1.89d0*SCALAR
        ELSE IF(ORDER.EQ.9)THEN
          SEGLENGTH=1.62d0*SCALAR
        ELSE IF(ORDER.EQ.8)THEN
          SEGLENGTH=1.26d0*SCALAR
        ELSE IF(ORDER.EQ.7)THEN
          SEGLENGTH=0.990d0*SCALAR
        ELSE IF(ORDER.EQ.6)THEN
          SEGLENGTH=1.5d0*SCALAR
        ENDIF
      ELSE IF (BRANCH.EQ.2) THEN !LAD
        IF (ORDER.EQ.11)THEN
          SEGLENGTH=3.0d0*SCALAR
        ELSE IF(ORDER.EQ.10)THEN
          SEGLENGTH=1.8d0*SCALAR
        ELSE IF(ORDER.EQ.9)THEN
          SEGLENGTH=1.6d0*SCALAR
        ELSE IF(ORDER.EQ.8)THEN
          SEGLENGTH=1.2d0*SCALAR
        ELSE IF(ORDER.EQ.7)THEN
          SEGLENGTH=1.0d0*SCALAR
        ELSE IF(ORDER.EQ.6)THEN
          SEGLENGTH=1.4d0*SCALAR
        ENDIF
      ELSE IF (BRANCH.EQ.3) THEN !CX
        IF (ORDER.EQ.11)THEN
          SEGLENGTH=0.0d0*SCALAR
        ELSE IF(ORDER.EQ.10)THEN
          SEGLENGTH=3.43d0*SCALAR
        ELSE IF(ORDER.EQ.9)THEN
          SEGLENGTH=3.18d0*SCALAR
        ELSE IF(ORDER.EQ.8)THEN
          SEGLENGTH=1.78d0*SCALAR
        ELSE IF(ORDER.EQ.7)THEN
          SEGLENGTH=1.60d0*SCALAR
        ELSE IF(ORDER.EQ.6)THEN
          SEGLENGTH=2.2d0*SCALAR
        ENDIF
      ENDIF

      RETURN
      END

C FE40 Functions
C ==============

