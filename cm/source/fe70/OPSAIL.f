      SUBROUTINE OPSAIL(XP,ERROR,*)

C#### Subroutine: OPSAIL
C###  Description:
C###    OPSAIL performs sail parameter output.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      REAL*8 XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER j,nj,nosect,np1,np2,np3,nv
      REAL*8 CLECH,CLUFF,COS,DRAFT,DXLECH,DXLUFF,SIN,XLECH,XLUFF,
     '  XYCHO1(2,2,2),
     '  XYLEC1(2,2,2),XYLEC2(2,2,2),XYLUF1(2,2,2),XYLUF2(2,2,2),
     '  XZLEC1(2,2,2),XZLEC2(2,2,2),XZLUF1(2,2,2),XZLUF2(2,2,2),
     '  ZCLEW(3),ZHEADL(3),ZHEADT(3),ZTACK(3)

      CALL ENTERS('OPSAIL',*9999)
      ERROR='out of date code'
      GO TO 9999

C LC 25/2/97 LC archived section :  Temporary MPN 12-Nov-94
      CALL EXITS('OPSAIL')
      RETURN
 9999 CALL ERRORS('OPSAIL',ERROR)
      CALL EXITS('OPSAIL')
      RETURN 1
      END


