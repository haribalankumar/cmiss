      SUBROUTINE GRDMLM(INDEX,ERROR,*)

C#### Subroutine: GRDMLM
C###  Description:
C###    GRDMLM draws grid lines of constant mu in lamda, mu coordinates.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'graf00.cmn'
!     Parameter List
      INTEGER INDEX
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,IBEG,IEND,iRMU,nj1,nj2
      REAL*8 P(3,50),RL,RLMAX,RMU,XM
      CHARACTER CHAR*10

      CALL ENTERS('GRDMLM',*9999)
      IF(NJT.EQ.2)THEN
        nj1=1
        nj2=2
      ELSE IF(NJT.EQ.3)THEN
        nj1=1
        nj2=3
      ENDIF
      XM=DBLE(ZMAX)/FOCUS
      RLMAX=DLOG(XM+DSQRT(XM*XM-1.0D0))
C     DO RMU=0.0D0,180.0D0,10.0D0
      DO iRMU=0,180,10
        RMU=DBLE(iRMU)
        DO i=1,50
          RL=RLMAX*DBLE(i-1)/49.0D0
          P(nj1,i)= FOCUS*DCOSH(RL)*DCOS(RMU*PI/180.0D0)
          P(nj2,i)=-FOCUS*DSINH(RL)*DSIN(RMU*PI/180.0D0)
        ENDDO
        CALL POLYLINE(INDEX,1,50,P,ERROR,*9999)
        DO i=1,50
          RL=RLMAX*DBLE(i-1)/49.0D0
          P(nj1,i)= FOCUS*DCOSH(RL)*DCOS(RMU*PI/180.0D0)
          P(nj2,i)= FOCUS*DSINH(RL)*DSIN(RMU*PI/180.0D0)
        ENDDO
        CALL POLYLINE(INDEX,1,50,P,ERROR,*9999)
        WRITE(CHAR,'(I3)') IFIX(REAL(RMU))
        CALL STRING_TRIM(CHAR,IBEG,IEND)
        CALL TEXT(1,1,CHAR(IBEG:IEND),P(1,50),ERROR,*9999)
      ENDDO

      CALL EXITS('GRDMLM')
      RETURN
 9999 CALL ERRORS('GRDMLM',ERROR)
      CALL EXITS('GRDMLM')
      RETURN 1
      END


