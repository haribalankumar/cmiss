      SUBROUTINE GRDXXY(INDEX,ERROR,*)

C#### Subroutine: GRDXXY
C###  Description:
C###    GRDXXY draws grid lines of constant x in x,y coordinates.

      IMPLICIT NONE
      INCLUDE 'graf00.cmn'
!     Parameter List
      INTEGER INDEX
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 DIFF,P(3,3),RFROMC,XINCR,XINCRT,XX
      CHARACTER CHAR*9

      CALL ENTERS('GRDXXY',*9999)
      XINCR=(DBLE(XMAX-XMIN))/10.0D0
C      CHAR=CFROMR(XINCR,'(E8.1)')
      WRITE(CHAR,'(E8.1)') XINCR
      XINCRT=RFROMC(CHAR)
      DIFF=DBLE(YMAX-YMIN)
      P(2,2)=DBLE(YMAX)-0.05D0*DIFF
      DO XX=0.0D0,DBLE(XMAX)-XINCRT,XINCRT
        P(2,1)=DBLE(YMIN)+0.05D0*DIFF
        P(1,1)=XX
        P(1,2)=XX
        CALL POLYLINE(INDEX,1,2,P,ERROR,*9999)
C        CHAR=CFROMR(XX,'(E9.2)')
        WRITE(CHAR,'(E9.2)') XX
        P(2,1)=P(2,1)-0.025d0*DIFF
        CALL TEXT(1,1,CHAR(3:5),P,ERROR,*9999)
      ENDDO
      DO XX=-XINCRT,DBLE(XMIN)+XINCRT,-XINCRT
        P(2,1)=DBLE(YMIN)+0.05D0*DIFF
        P(1,1)=XX
        P(1,2)=XX
        CALL POLYLINE(INDEX,1,2,P,ERROR,*9999)
C        CHAR=CFROMR(XX,'(E9.2)')
        WRITE(CHAR,'(E9.2)') XX
        P(2,1)=P(2,1)-0.025D0*DIFF
        CALL TEXT(1,1,'-'//CHAR(3:5),P,ERROR,*9999)
      ENDDO

      CALL EXITS('GRDXXY')
      RETURN
 9999 CALL ERRORS('GRDXXY',ERROR)
      CALL EXITS('GRDXXY')
      RETURN 1
      END


