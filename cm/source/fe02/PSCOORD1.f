      SUBROUTINE PSCOORD1(XD,XWC,YWC,ERROR,*)

C#### Subroutine: PSCOORD1
C###  Description:
C###    PSCOORD1 finds prolate spheroidal coords XD(nj) of a point
C###    located on Hammer projection map at XWC,YWC.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER isub
      REAL*8 DET,DF1DX1,DF1DX2,DF2DX1,DF2DX2,
     '  DI1DX1,DI1DX2,DI2DX1,DI2DX2,DIF,DIFMAX,DK,DKDX1,DKDX2,
     '  F1,F2,THETA,TOL,RK,RMP,RMU,TP,XD(3),XWC,YWC
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER iter
      REAL*8 DELTAX(2),XNEW(2),XOLD(2)
      DATA TOL/1.0D-6/

      CALL ENTERS('PSCOORD1',*9999)
      WRITE(OP_STRING,'('' Map coordinates: '',2E11.4)') XWC,YWC
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      XOLD(1) = PI/2.2D0
      XOLD(2) = PI/2.3D0
      DO 150 iter = 1, 10000
        RMU   = XOLD(1)
        THETA = XOLD(2)
        RMP   = RMU-PI/2.0D0
        TP    = (THETA-PI)/2.0D0
        RK    = (1.0D0 + DCOS(RMP)*DCOS(TP))**(-0.5D0)
        F1    = XWC + RK*DCOS(RMP)*DSIN(TP)
        F2    = YWC - RK*DSIN(RMP)
        DK    =-0.5D0*(1.0D0+DCOS(RMP)*DCOS(TP))**(-1.5D0)
        DKDX1 =-DK*DSIN(RMP)*DCOS(TP)
        DKDX2 =-DK*0.5D0*DCOS(RMP)*DSIN(TP)
        DF1DX1=-RK*DSIN(RMP)*DSIN(TP) + DCOS(RMP)*DSIN(TP)*DKDX1
        DF1DX2= 0.5D0*RK*DCOS(RMP)*DCOS(TP) + DCOS(RMP)*DSIN(TP)*DKDX2
        DF2DX1=-RK*DCOS(RMP)- DSIN(RMP)*DKDX1
        DF2DX2=-DSIN(RMP)*DKDX2
        DET   = DF1DX1*DF2DX2 - DF1DX2*DF2DX1
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,
     '      '('' Iteration '',i4,3x,''rmu='',e10.4,'' theta='','
     '      //'E10.4,/,18x,''rmp='',e10.4,'' tp='',e10.4,'
     '      //''' rk='',e10.4,'' det='',e10.4)')
     '      iter,rmu,theta,rmp,tp,rk,det
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        DI1DX1 = DF2DX2/DET
        DI1DX2 =-DF1DX2/DET
        DI2DX1 =-DF2DX1/DET
        DI2DX2 = DF1DX1/DET
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'(10x,'' dfdu'',20x,''didu'',20x,''f'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        DELTAX(1) = -DI1DX1*F1 -DI1DX2*F2
        DELTAX(2) = -DI2DX1*F1 -DI2DX2*F2
        XNEW(1) = XOLD(1) + DELTAX(1)
        XNEW(2) = XOLD(2) + DELTAX(2)
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,
     '      '(10X,'' DELTAX='',2(1X,E10.4),'' U ='',2(1X,E10.4))')
     '      DELTAX,XNEW
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        DIFMAX = 0.0D0
        DO isub = 1, 2
          DIF = DABS(XNEW(isub) - XOLD(isub))
          IF(DIF.GT.DIFMAX) DIFMAX = DIF
        ENDDO
        IF(DIFMAX.LT.TOL) GOTO 200
        XOLD(1) = XNEW(1)
        XOLD(2) = XNEW(2)
 150  CONTINUE
 200  CONTINUE
      XD(1)=1.0D0
      XD(2)=XNEW(1)
      XD(3)=XNEW(2)
      WRITE(OP_STRING,'('' Prolate coords:   mu='',F5.1,'
     '  //''' theta='',F5.1)') XD(2)*180.0D0/PI,XD(3)*180.0D0/PI
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      CALL EXITS('PSCOORD1')
      RETURN
 9999 CALL ERRORS('PSCOORD1',ERROR)
      CALL EXITS('PSCOORD1')
      RETURN 1
      END


