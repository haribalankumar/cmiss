      SUBROUTINE GRDYXY(INDEX,ERROR,*)

C#### Subroutine: GRDYXY
C###  Description:
C###    GRDYXY draws grid lines of constant y in x,y coordinates.

      IMPLICIT NONE
      INCLUDE 'graf00.cmn'
!     Parameter List
      INTEGER INDEX
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 DIFF,P(3,3),RFROMC,YINCR,YINCRT,YY
      CHARACTER CHAR*9

      CALL ENTERS('GRDYXY',*9999)
      YINCR=DBLE(YMAX-YMIN)/10.0D0
C      CHAR=CFROMR(YINCR,'(E8.1)')
      WRITE(CHAR,'(E8.1)') YINCR
      YINCRT=RFROMC(CHAR)
      DIFF=DBLE(XMAX-XMIN)
      P(1,2)=DBLE(XMAX)-0.05D0*DIFF
      DO YY=0.0D0,DBLE(YMAX)-YINCRT,YINCRT
        P(1,1)=DBLE(XMIN)+0.05D0*DIFF
        P(2,1)=YY
        P(2,2)=YY
        CALL POLYLINE(INDEX,1,2,P,ERROR,*9999)
C        CHAR=CFROMR(YY,'(E9.2)')
        WRITE(CHAR,'(E9.2)') YY
        P(1,1)=P(1,1)-0.025D0*DIFF
        CALL TEXT(1,1,CHAR(3:5),P,ERROR,*9999)
      ENDDO
C      CHAR=CFROMR(10.0D0*YINCR,'(E8.1)')
      WRITE(CHAR,'(E8.1)') 10.0D0*YINCR
      P(1,1)=DBLE(XMIN)+0.02D0*DIFF
      P(2,1)=P(2,1)+0.5D0*YINCRT
      CALL TEXT(1,1,'*'//CHAR(5:8),P,ERROR,*9999)
      DO YY=-YINCRT,DBLE(YMIN)+YINCRT,-YINCRT
        P(1,1)=DBLE(XMIN)+0.05D0*DIFF
        P(2,1)=YY
        P(2,2)=YY
        CALL POLYLINE(INDEX,1,2,P,ERROR,*9999)
C        CHAR=CFROMR(YY,'(E9.2)')
        WRITE(CHAR,'(E9.2)') YY
        P(1,1)=P(1,1)-0.025D0*DIFF
        CALL TEXT(1,1,'-'//CHAR(3:5),P,ERROR,*9999)
      ENDDO

      CALL EXITS('GRDYXY')
      RETURN
 9999 CALL ERRORS('GRDYXY',ERROR)
      CALL EXITS('GRDYXY')
      RETURN 1
      END


