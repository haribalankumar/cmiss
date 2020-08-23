      SUBROUTINE SE2CHANGE(SE2,ERROR,*)

C#### Subroutine: SE2CHANGE
C###  Description:
C###    SE2CHANGE matches the scale-factors in the SE2 array and changes
C###    different/conflicing factors to the largest of the two.
C###    Called from UPLINE.

      IMPLICIT NONE
!     Parameter List
      REAL*8 SE2(18,2)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER INDX(2,2),indx1,indx2

      CALL ENTERS('SE2CHANGE',*9999)
C  First derivatives
      IF(SE2(15,1).EQ.3) THEN
        INDX(1,1)=4
        INDX(2,1)=11
      ELSE
        INDX(1,1)=INT(SE2(15,1))
        INDX(2,1)=INT(SE2(15,1))+7
      ENDIF
      IF(SE2(16,1).EQ.SE2(16,2)) THEN !Same node order?
        indx1=1
        indx2=2
      ELSE
        indx1=2
        indx2=1
      ENDIF
      IF(SE2(15,2).EQ.3) THEN
        INDX(indx1,2)=4
        INDX(indx2,2)=11
      ELSE
        INDX(indx1,2)=INT(SE2(15,2))
        INDX(indx2,2)=INT(SE2(15,2))+7
      ENDIF
C  Node 1
      IF(SE2(INDX(1,1),1).GT.SE2(INDX(1,2),2)) THEN
        SE2(INDX(1,2),2)=SE2(INDX(1,1),1)
      ELSEIF(SE2(INDX(1,1),1).LT.SE2(INDX(1,2),2)) THEN
        SE2(INDX(1,1),1)=SE2(INDX(1,2),2)
      ENDIF
C  Node 2
      IF(SE2(INDX(2,1),1).GT.SE2(INDX(2,2),2)) THEN
        SE2(INDX(2,2),2)=SE2(INDX(2,1),1)
      ELSEIF(SE2(INDX(2,1),1).LT.SE2(INDX(2,2),2)) THEN
        SE2(INDX(2,1),1)=SE2(INDX(2,2),2)
      ENDIF
C  Cross derivatives
      IF(SE2(18,1).NE.0) THEN
        IF(SE2(15,1).EQ.1.AND.SE2(18,1).EQ.2.OR.
     '     SE2(15,1).EQ.2.AND.SE2(18,1).EQ.1) THEN
          INDX(1,1)=3
          INDX(2,1)=10
        ELSEIF(SE2(15,1).EQ.1.AND.SE2(18,1).EQ.3.OR.
     '         SE2(15,1).EQ.3.AND.SE2(18,1).EQ.1) THEN
          INDX(1,1)=5
          INDX(2,1)=12
        ELSEIF(SE2(15,1).EQ.2.AND.SE2(18,1).EQ.3.OR.
     '         SE2(15,1).EQ.3.AND.SE2(18,1).EQ.2) THEN
          INDX(1,1)=6
          INDX(2,1)=13
        ENDIF
        IF(SE2(15,2).EQ.1.AND.SE2(18,2).EQ.2.OR.
     '     SE2(15,2).EQ.2.AND.SE2(18,2).EQ.1) THEN
          INDX(indx1,2)=3
          INDX(indx2,2)=10
        ELSEIF(SE2(15,2).EQ.1.AND.SE2(18,2).EQ.3.OR.
     '         SE2(15,2).EQ.3.AND.SE2(18,2).EQ.1) THEN
          INDX(indx1,2)=5
          INDX(indx2,2)=12
        ELSEIF(SE2(15,2).EQ.2.AND.SE2(18,2).EQ.3.OR.
     '         SE2(15,2).EQ.3.AND.SE2(18,2).EQ.2) THEN
          INDX(indx1,2)=6
          INDX(indx2,2)=13
        ENDIF
C  Node 1
        IF(SE2(INDX(1,1),1).GT.SE2(INDX(1,2),2)) THEN
          SE2(INDX(1,2),2)=SE2(INDX(1,1),1)
        ELSEIF(SE2(INDX(1,1),1).LT.SE2(INDX(1,2),2)) THEN
          SE2(INDX(1,1),1)=SE2(INDX(1,2),2)
        ENDIF
C  Node 2
        IF(SE2(INDX(2,1),1).GT.SE2(INDX(2,2),2)) THEN
          SE2(INDX(2,2),2)=SE2(INDX(2,1),1)
        ELSEIF(SE2(INDX(2,1),1).LT.SE2(INDX(2,2),2)) THEN
          SE2(INDX(2,1),1)=SE2(INDX(2,2),2)
        ENDIF
      ENDIF

      CALL EXITS('SE2CHANGE')
      RETURN
 9999 CALL ERRORS('SE2CHANGE',ERROR)
      CALL EXITS('SE2CHANGE')
      RETURN 1
      END



