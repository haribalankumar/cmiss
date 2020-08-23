      SUBROUTINE SHAPEC(ERROR,*)

C#### Subroutine: SHAPEC
C###  Description:
C###    SHAPEC shapes individual chord.

      IMPLICIT NONE
      INCLUDE 'graf00.cmn'
      INCLUDE 'mxch.inc'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 P(3,2)

      CALL ENTERS('SHAPEC',*9999)
      CALL ACWK(6,0,ERROR,*9999)
C      CALL CRHERM(ISEG,ISLINE(1,1),ISNONO,6,NLCHOR(1),
C    '  NLCHOR(NLCHOR(0)),NPL,CSEG,XP,ZP)
C      CALL GCLRWK(6)  !is this needed? AAY
      CALL DAWK(6,0,ERROR,*9999)
      CALL ACWK(7,0,ERROR,*9999)
      P(2,1)=0.d0
      P(2,2)=0.d0
      P(3,1)=DBLE(ZMIN)
      P(3,2)=DBLE(ZMAX)
      CALL POLYLINE(1,7,2,P,ERROR,*9999)
      CALL DAWK(7,0,ERROR,*9999)

      CALL EXITS('SHAPEC')
      RETURN
 9999 CALL ERRORS('SHAPEC',ERROR)
      CALL EXITS('SHAPEC')
      RETURN 1
      END


