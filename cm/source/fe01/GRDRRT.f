      SUBROUTINE GRDRRT(INDEX,ERROR,*)

C#### Subroutine: GRDRRT
C###  Description:
C###    GRDRRT draws grid lines of constant radius in r, theta
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

      CALL ENTERS('GRDRRT',*9999)
      R_MAX=DMIN1(0.95D0*DBLE(YMAX),0.95D0*DBLE(ZMAX))
      RINCR=R_MAX/10.0D0
C     DO RR=RINCR,R_MAX,RINCR
      DO iRR=1,10
        RR=DBLE(iRR)*RINCR
        DO i=1,50
          THETA=2.0D0*PI*DBLE(i-1)/49.0D0
          P(1,i)=RR*DCOS(THETA)
          P(2,i)=RR*DSIN(THETA)
        ENDDO
        CALL POLYLINE(INDEX,1,50,P,ERROR,*9999)
      ENDDO

      CALL EXITS('GRDRRT')
      RETURN
 9999 CALL ERRORS('GRDRRT',ERROR)
      CALL EXITS('GRDRRT')
      RETURN 1
      END


