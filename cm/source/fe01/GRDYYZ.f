      SUBROUTINE GRDYYZ(INDEX,ERROR,*)

C#### Subroutine: GRDYYZ
C###  Description:
C###    GRDYYZ draws grid lines of constant y in y,z coordinates.

      IMPLICIT NONE
      INCLUDE 'graf00.cmn'
!     Parameter List
      INTEGER INDEX
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 DIFF,P(3,3),RFROMC,YINCR,YINCRT,YY
      CHARACTER CHAR*9

      CALL ENTERS('GRDYYZ',*9999)
      YINCR=DBLE(YMAX-YMIN)/10.0D0
C      CHAR=CFROMR(YINCR,'(E8.1)')
      WRITE(CHAR,'(E8.1)') YINCR
      YINCRT=RFROMC(CHAR)
      DIFF=DBLE(ZMAX-ZMIN)
      P(3,2)=DBLE(ZMAX)-0.05D0*DIFF
      DO YY=0.0D0,DBLE(YMAX)-YINCRT,YINCRT
        P(3,1)=DBLE(ZMIN)+0.05D0*DIFF
        P(2,1)=YY
        P(2,2)=YY
        CALL POLYLINE(INDEX,2,2,P,ERROR,*9999)
C        CHAR=CFROMR(YY,'(E9.2)')
        WRITE(CHAR,'(E9.2)') YY
        P(3,1)=P(3,1)-0.025D0*DIFF
        CALL TEXT(1,2,CHAR(3:5),P,ERROR,*9999)
      ENDDO
      DO YY=-YINCRT,DBLE(YMIN)+YINCRT,-YINCRT
        P(3,1)=DBLE(ZMIN)+0.05D0*DIFF
        P(2,1)=YY
        P(2,2)=YY
        CALL POLYLINE(INDEX,2,2,P,ERROR,*9999)
C        CHAR=CFROMR(YY,'(E9.2)')
        WRITE(CHAR,'(E9.2)') YY
        P(3,1)=P(3,1)-0.025D0*DIFF
        CALL TEXT(1,2,'-'//CHAR(3:5),P,ERROR,*9999)
      ENDDO

      CALL EXITS('GRDYYZ')
      RETURN
 9999 CALL ERRORS('GRDYYZ',ERROR)
      CALL EXITS('GRDYYZ')
      RETURN 1
      END


