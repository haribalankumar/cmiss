      SUBROUTINE GRDTRT(INDEX,iw,ERROR,*)

C#### Subroutine: GRDTRT
C###  Description:
C###    GRDTRT draws grid lines of constant theta in r, theta
C###    coordinates.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'graf00.cmn'
!     Parameter List
      INTEGER INDEX
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IBEG,IEND,iTHETA,iw,nj1,nj2
      REAL*8 P(3,2),R_MAX,THETA
      CHARACTER CHAR*10

      CALL ENTERS('GRDTRT',*9999)
      IF(iw.EQ.1)THEN
        nj1=1
        nj2=2
      ELSE IF(iw.EQ.2)THEN
        nj1=2
        nj2=3
      ENDIF
      R_MAX=DMIN1(0.95D0*DBLE(YMAX),0.95D0*DBLE(ZMAX))
      DO iTHETA=10,360,10
        THETA=DBLE(iTHETA)
        P(nj1,1)=0.0D0
        P(nj2,1)=0.0D0
        P(nj1,2)=R_MAX*DCOS(THETA*PI/180.0D0)
        P(nj2,2)=R_MAX*DSIN(THETA*PI/180.0D0)
        CALL POLYLINE(INDEX,iw,2,P,ERROR,*9999)
        WRITE(CHAR,'(I3)') IFIX(REAL(THETA))
        CALL STRING_TRIM(CHAR,IBEG,IEND)
        CALL TEXT(1,iw,CHAR(IBEG:IEND),P(1,2),ERROR,*9999)
      ENDDO

      CALL EXITS('GRDTRT')
      RETURN
 9999 CALL ERRORS('GRDTRT',ERROR)
      CALL EXITS('GRDTRT')
      RETURN 1
      END


