      INTEGER FUNCTION SEGNUM(ORDER,BRANCH)

C#### FUNCTION: SEGNUM
C#### Type: INTEGER
C###  Description:
C###   SEGNUM calculates the number of segments a vessel
C###    will have for a given order

      IMPLICIT NONE
!     Parameter List
      INTEGER ORDER,BRANCH
C      REAL*8 NORM_RAND_NUM


      IF(BRANCH.EQ.1) THEN !RCA
        IF(ORDER.EQ.11)THEN
          SEGNUM=26
        ELSE IF(ORDER.EQ.10)THEN
          SEGNUM=9
        ELSE IF(ORDER.EQ.9)THEN
          SEGNUM=3
        ELSE IF(ORDER.EQ.8)THEN
          SEGNUM=3
        ELSE IF(ORDER.EQ.7)THEN
          SEGNUM=2
        ELSE IF(ORDER.EQ.6)THEN
          SEGNUM=1
        ENDIF
      ELSE IF (BRANCH.EQ.2) THEN !LAD
        IF (ORDER.EQ.11)THEN
          SEGNUM=17
        ELSE IF(ORDER.EQ.10)THEN
          SEGNUM=9
        ELSE IF(ORDER.EQ.9)THEN
          SEGNUM=5
        ELSE IF(ORDER.EQ.8)THEN
          SEGNUM=3
        ELSE IF(ORDER.EQ.7)THEN
          SEGNUM=2
        ELSE IF(ORDER.EQ.6)THEN
          SEGNUM=1
        ENDIF
      ELSE IF (BRANCH.EQ.3) THEN !CX
        IF (ORDER.EQ.11)THEN
          SEGNUM=0
        ELSE IF(ORDER.EQ.10)THEN
          SEGNUM=14
        ELSE IF(ORDER.EQ.9)THEN
          SEGNUM=7
        ELSE IF(ORDER.EQ.8)THEN
          SEGNUM=3
        ELSE IF(ORDER.EQ.7)THEN
          SEGNUM=2
        ELSE IF(ORDER.EQ.6)THEN
          SEGNUM=1
        ENDIF
      ENDIF

      RETURN
      END


