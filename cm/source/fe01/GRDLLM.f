      SUBROUTINE GRDLLM(INDEX,ERROR,*)

C#### Subroutine: GRDLLM
C###  Description:
C###    GRDLLM draws grid lines of constant lamda in lamda, mu
C###    coordinates.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'graf00.cmn'
!     Parameter List
      INTEGER INDEX
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,iRL,nj1,nj2
      REAL*8 P(3,100),RL,RLINCR,RLMAX,RMU,XM

      CALL ENTERS('GRDLLM',*9999)
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' zmax='',E12.3,'' focus='',E12.3)')
     '    ZMAX,FOCUS
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF
      XM=DBLE(ZMAX)/FOCUS
      IF(NJT.EQ.2)THEN
        nj1=1
        nj2=2
      ELSE IF(NJT.EQ.3)THEN
        nj1=1
        nj2=3
      ENDIF
      IF(XM.GT.1.0D0) THEN
        RLMAX=DLOG(XM+DSQRT(XM*XM-1.0D0))
      ELSE
        RLMAX=0.0D0
      ENDIF
      RLINCR=RLMAX/10.0D0
C     DO RL=RLINCR,RLMAX,RLINCR
      DO iRL=1,10
        RL=RLMAX*RLINCR
        DO i=1,100
          RMU=2.0D0*PI*DBLE(i-1)/99.0D0
          P(nj1,i)=FOCUS*DCOSH(RL)*DCOS(RMU)
          P(nj2,i)=FOCUS*DSINH(RL)*DSIN(RMU)
        ENDDO
        CALL POLYLINE(INDEX,1,100,P,ERROR,*9999)
      ENDDO

      CALL EXITS('GRDLLM')
      RETURN
 9999 CALL ERRORS('GRDLLM',ERROR)
      CALL EXITS('GRDLLM')
      RETURN 1
      END


