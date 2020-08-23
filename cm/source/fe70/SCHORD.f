      SUBROUTINE SCHORD(X1,X2,XSCALE,ERROR,*)

C#### Subroutine: SCHORD
C###  Description:
C###    SCHORD shows shape of individual chord.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'graf00.cmn'
      INCLUDE 'ntsg00.cmn'
c     INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
!     Parameter List
      REAL*8 XSCALE,X1(2,2,*),X2(2,2,*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER n
      REAL X(2),XX(40,2),Y(2),Z(2)
      REAL*8 XH3,XI

      CALL ENTERS('SCHORD',*9999)
      CALL ACWK(18,0,ERROR,*9999)
      CALL GCLRWK(18,1)
      DO n=1,21
        XI=DBLE(n-1)/20.D0
        XX(n,1)=XH3(1,X1,XI)*XSCALE
        XX(n,2)=XH3(2,X1,XI)*XSCALE
      ENDDO
      DO n=21,40
        XI=DBLE(n-20)/20.D0
        XX(n,1)=XH3(1,X2,XI)*XSCALE
        XX(n,2)=XH3(2,X2,XI)*XSCALE
      ENDDO
      X(1)=XX( 1,1)
      X(2)=XX(40,1)
      Y(1)=0.
      Y(2)=0.
c     CALL GSLN(GLDOT)
      CALL GPL(2,X,Y)
c     CALL GSLN(GLSOLI)
      CALL GPL(40,XX(1,1),XX(1,2))
      CALL DAWK(18,0,ERROR,*9999)
      CALL ACWK(19,0,ERROR,*9999)
      Y(1)=0.D0
      Y(2)=0.D0
      Z(1)=ZMIN
      Z(2)=ZMAX
      CALL GPL(2,Y,Z)
      CALL DAWK(19,0,ERROR,*9999)

 9998 CALL EXITS('SCHORD')
      RETURN
 9999 CALL ERRORS('SCHORD',ERROR)
      CALL EXITS('SCHORD')
      RETURN 1
      END


