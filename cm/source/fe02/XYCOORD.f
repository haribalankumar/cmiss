      SUBROUTINE XYCOORD(X,Y,XE,XI1,XI2,ERROR,*)

C#### Subroutine: XYCOORD
C###  Description:
C###    XYCOORD finds Xi_1,Xi_2 coordinates of the point X,Y in a
C###    bilinear element with nodes XE(nn,nj) by Newton iteration.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      REAL*8 X,XE(NSM,NJM),XI1,XI2,Y
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER KOUNT
      REAL*8 DELTA,DXI1,DXI2,F1,F11,F12,F2,F21,F22,FN1,FN2
      LOGICAL NOTFOUND

      FN1(XI1,XI2)=(1.0D0-XI1)*((1.0D0-XI2)*XE(1,1)+XI2*XE(3,1))
     '                  + XI1 *((1.0D0-XI2)*XE(2,1)+XI2*XE(4,1))
      FN2(XI1,XI2)=(1.0D0-XI1)*((1.0D0-XI2)*XE(1,2)+XI2*XE(3,2))
     '                  + XI1 *((1.0D0-XI2)*XE(2,2)+XI2*XE(4,2))

      CALL ENTERS('XYCOORD',*9999)
      NOTFOUND=.TRUE.
      KOUNT=0
      DO WHILE(NOTFOUND)
        KOUNT=KOUNT+1
        F11=(1.0D0-XI2)*(XE(2,1)-XE(1,1))+XI2*(XE(4,1)-XE(3,1))
        F12=(1.0D0-XI1)*(XE(3,1)-XE(1,1))+XI1*(XE(4,1)-XE(2,1))
        F21=(1.0D0-XI2)*(XE(2,2)-XE(1,2))+XI2*(XE(4,2)-XE(3,2))
        F22=(1.0D0-XI1)*(XE(3,2)-XE(1,2))+XI1*(XE(4,2)-XE(2,2))

        DELTA=F11*F22-F12*F21
        F1=FN1(XI1,XI2) !is x interpolated value
        F2=FN2(XI1,XI2) !is y interpolated value
        IF(ITYP10(1).EQ.4) THEN !prolate coords
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' Input theta='',E12.3)') X
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,
     '        '('' Interp theta (before correction)='',E12.3)') F1
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          IF(F1.GT.2.0D0*PI) THEN
             F1=F1-2.0D0*PI
          ELSE IF(F1.LT.0.0D0) THEN
             F1=F1+2.0D0*PI
          ENDIF
C NOTE!!!!  code below sb better than above but isn't
c          IF(F1.GT.2.0D0*PI.OR.F1.LT.0.0D0) THEN
c            FMULT=(F1-AMOD(F1,2.0D0*PI))/(2.0D0*PI)
c            F1=F1-2.0D0*PI*FMULT !F1 is theta coord
c          ENDIF
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,
     '        '('' Interp theta  (after correction)='',E12.3)') F1
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
        ENDIF
        F1=F1-X !is difference between interpolated and given x
        F2=F2-Y !is difference between interpolated and given y
        DXI1=(F12*F2-F22*F1)/DELTA !is increment in Xi_1
        DXI2=(F21*F1-F11*F2)/DELTA !is increment in Xi_2

        XI1=XI1+DXI1
        XI2=XI2+DXI2
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,
     '      '('' Xi_1='',E12.3,'' Xi_2='',E12.3)') XI1,XI2
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        IF((DXI1*DXI1+DXI2*DXI2).LT.1.0D-6) NOTFOUND=.FALSE.
        IF(KOUNT.GT.20) THEN
          WRITE(OP_STRING,'('' More than 20 iterations'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          NOTFOUND=.FALSE.
        ENDIF
      ENDDO

      CALL EXITS('XYCOORD')
      RETURN
 9999 CALL ERRORS('XYCOORD',ERROR)
      CALL EXITS('XYCOORD')
      RETURN 1
      END


