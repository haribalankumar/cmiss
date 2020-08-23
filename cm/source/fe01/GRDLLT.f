      SUBROUTINE GRDLLT(INDEX,ERROR,*)

C#### Subroutine: GRDLLT
C###  Description:
C###    GRDLLT draws grid lines of constant lamda in lamda, theta
C###    coordinates.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'graf00.cmn'
!     Parameter List
      INTEGER INDEX
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,iRR
      REAL*8 P(3,50),RINCR,RR,R_MAX,THETA

      CALL ENTERS('GRDLLT',*9999)
      R_MAX=DMIN1(0.95D0*DBLE(YMAX),0.95D0*DBLE(ZMAX))
      RINCR=R_MAX/4.0D0
C     DO RR=RINCR,R_MAX,RINCR
      DO iRR=1,4
        RR=R_MAX*RINCR
        DO i=1,50
          THETA=2.0D0*PI*DBLE(i-1)/49.0D0
          P(2,i)=RR*DCOS(THETA)
          P(3,i)=RR*DSIN(THETA)
        ENDDO
        CALL POLYLINE(INDEX,2,50,P,ERROR,*9999)
      ENDDO

      CALL EXITS('GRDLLT')
      RETURN
 9999 CALL ERRORS('GRDLLT',ERROR)
      CALL EXITS('GRDLLT')
      RETURN 1
      END


