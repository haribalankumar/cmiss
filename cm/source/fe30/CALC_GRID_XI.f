      SUBROUTINE CALC_GRID_XI(NITB,NWQ,nq,NXQ,XI,ERROR,*)

C#### Subroutine: CALC_GRID_XI
C###  Description:
C###    CALC_GRID_XI calculates the xi position of a grid point
C###    nq within a local quadratic element in X space.
C**** Written by Martin Buist, 21 June 1999

      IMPLICIT NONE

      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'

!     Parameter List
      INTEGER NITB,NWQ(8,0:NQM),nq,NXQ(-NIM:NIM,0:4,0:NQM)
      REAL*8 XI(3)
      CHARACTER ERROR*(*)
!     Local Variables

      CALL ENTERS('CALC_GRID_XI',*9999)

      IF(NWQ(1,nq).EQ.0) THEN !internal
        IF(NITB.EQ.1) THEN
          XI(1)=0.5d0
          XI(2)=0.0d0
          XI(3)=0.0d0
        ELSE IF(NITB.EQ.2) THEN
          XI(1)=0.5d0
          XI(2)=0.5d0
          XI(3)=0.0d0
        ELSE IF(NITB.EQ.3) THEN
          XI(1)=0.5d0
          XI(2)=0.5d0
          XI(3)=0.5d0
        ENDIF !xi dirns
      ELSE !external
        IF(NITB.EQ.1) THEN !3-1=2 options
          IF(NXQ(-1,1,nq).EQ.NWQ(1,nq)) THEN
            XI(1)=1.0d0
            XI(2)=0.0d0
            XI(3)=0.0d0
          ELSE
            XI(1)=0.0d0
            XI(2)=0.0d0
            XI(3)=0.0d0
          ENDIF !options
        ELSE IF(NITB.EQ.2) THEN !9-1=8 options
          IF(NXQ(-1,1,nq).EQ.NWQ(1,nq)) THEN
            XI(1)=1.0d0
            XI(2)=0.5d0
            XI(3)=0.0d0
          ELSE IF(NXQ(1,1,nq).EQ.NWQ(1,nq)) THEN
            XI(1)=0.0d0
            XI(2)=0.5d0
            XI(3)=0.0d0
          ELSE IF(NXQ(-2,1,nq).EQ.NWQ(1,nq)) THEN
            XI(1)=0.5d0
            XI(2)=1.0d0
            XI(3)=0.0d0
          ELSE IF(NXQ(2,1,nq).EQ.NWQ(1,nq)) THEN
            XI(1)=0.5d0
            XI(2)=0.0d0
            XI(3)=0.0d0
          ELSE IF((NXQ(1,1,NXQ(2,1,nq)).EQ.NWQ(1,nq)).OR.
     '      (NXQ(2,1,NXQ(1,1,nq)).EQ.NWQ(1,nq))) THEN
            XI(1)=0.0d0
            XI(2)=0.0d0
            XI(3)=0.0d0
          ELSE IF((NXQ(-1,1,NXQ(-2,1,nq)).EQ.NWQ(1,nq)).OR.
     '      (NXQ(-2,1,NXQ(-1,1,nq)).EQ.NWQ(1,nq))) THEN
            XI(1)=1.0d0
            XI(2)=1.0d0
            XI(3)=0.0d0
          ELSE IF((NXQ(2,1,NXQ(-1,1,nq)).EQ.NWQ(1,nq)).OR.
     '      (NXQ(-1,1,NXQ(2,1,nq)).EQ.NWQ(1,nq))) THEN
            XI(1)=1.0d0
            XI(2)=0.0d0
            XI(3)=0.0d0
          ELSE IF((NXQ(1,1,NXQ(-2,1,nq)).EQ.NWQ(1,nq)).OR.
     '      (NXQ(-2,1,NXQ(1,1,nq)).EQ.NWQ(1,nq))) THEN
            XI(1)=0.0d0
            XI(2)=1.0d0
            XI(3)=0.0d0
          ELSE
            XI(1)=0.0d0
            XI(2)=0.0d0
            XI(3)=0.0d0
          ENDIF !options
        ELSE IF(NITB.EQ.3) THEN !27-1=26 options
          IF(NXQ(-1,1,nq).EQ.NWQ(1,nq)) THEN
            XI(1)=1.0d0
            XI(2)=0.5d0
            XI(3)=0.5d0
          ELSE IF(NXQ(1,1,nq).EQ.NWQ(1,nq)) THEN
            XI(1)=0.0d0
            XI(2)=0.5d0
            XI(3)=0.5d0
          ELSE IF(NXQ(-2,1,nq).EQ.NWQ(1,nq)) THEN
            XI(1)=0.5d0
            XI(2)=1.0d0
            XI(3)=0.5d0
          ELSE IF(NXQ(2,1,nq).EQ.NWQ(1,nq)) THEN
            XI(1)=0.5d0
            XI(2)=0.0d0
            XI(3)=0.5d0
          ELSE IF(NXQ(-3,1,nq).EQ.NWQ(1,nq)) THEN
            XI(1)=0.5d0
            XI(2)=0.5d0
            XI(3)=1.0d0
          ELSE IF(NXQ(3,1,nq).EQ.NWQ(1,nq)) THEN
            XI(1)=0.5d0
            XI(2)=0.5d0
            XI(3)=0.0d0
          ELSE IF((NXQ(1,1,NXQ(2,1,nq)).EQ.NWQ(1,nq)).OR.
     '      (NXQ(2,1,NXQ(1,1,nq)).EQ.NWQ(1,nq))) THEN
            XI(1)=0.0d0
            XI(2)=0.0d0
            XI(3)=0.5d0
          ELSE IF((NXQ(-1,1,NXQ(-2,1,nq)).EQ.NWQ(1,nq)).OR.
     '      (NXQ(-2,1,NXQ(-1,1,nq)).EQ.NWQ(1,nq))) THEN
            XI(1)=1.0d0
            XI(2)=1.0d0
            XI(3)=0.5d0
          ELSE IF((NXQ(-1,1,NXQ(2,1,nq)).EQ.NWQ(1,nq)).OR.
     '      (NXQ(2,1,NXQ(-1,1,nq)).EQ.NWQ(1,nq))) THEN
            XI(1)=1.0d0
            XI(2)=0.0d0
            XI(3)=0.5d0
          ELSE IF((NXQ(1,1,NXQ(-2,1,nq)).EQ.NWQ(1,nq)).OR.
     '      (NXQ(-2,1,NXQ(1,1,nq)).EQ.NWQ(1,nq))) THEN
            XI(1)=0.0d0
            XI(2)=1.0d0
            XI(3)=0.5d0
          ELSE IF((NXQ(3,1,NXQ(2,1,nq)).EQ.NWQ(1,nq)).OR.
     '      (NXQ(2,1,NXQ(3,1,nq)).EQ.NWQ(1,nq))) THEN
            XI(1)=0.5d0
            XI(2)=0.0d0
            XI(3)=0.0d0
          ELSE IF((NXQ(3,1,NXQ(-2,1,nq)).EQ.NWQ(1,nq)).OR.
     '      (NXQ(-2,1,NXQ(3,1,nq)).EQ.NWQ(1,nq))) THEN
            XI(1)=0.5d0
            XI(2)=1.0d0
            XI(3)=0.0d0
          ELSE IF((NXQ(1,1,NXQ(3,1,nq)).EQ.NWQ(1,nq)).OR.
     '      (NXQ(3,1,NXQ(1,1,nq)).EQ.NWQ(1,nq))) THEN
            XI(1)=0.0d0
            XI(2)=0.5d0
            XI(3)=0.0d0
          ELSE IF((NXQ(3,1,NXQ(-1,1,nq)).EQ.NWQ(1,nq)).OR.
     '      (NXQ(-1,1,NXQ(3,1,nq)).EQ.NWQ(1,nq))) THEN
            XI(1)=1.0d0
            XI(2)=0.5d0
            XI(3)=0.0d0
          ELSE IF((NXQ(-3,1,NXQ(2,1,nq)).EQ.NWQ(1,nq)).OR.
     '      (NXQ(2,1,NXQ(-3,1,nq)).EQ.NWQ(1,nq))) THEN
            XI(1)=0.5d0
            XI(2)=0.0d0
            XI(3)=1.0d0
          ELSE IF((NXQ(-3,1,NXQ(-2,1,nq)).EQ.NWQ(1,nq)).OR.
     '      (NXQ(-2,1,NXQ(-3,1,nq)).EQ.NWQ(1,nq))) THEN
            XI(1)=0.5d0
            XI(2)=1.0d0
            XI(3)=1.0d0
          ELSE IF((NXQ(1,1,NXQ(-3,1,nq)).EQ.NWQ(1,nq)).OR.
     '      (NXQ(-3,1,NXQ(1,1,nq)).EQ.NWQ(1,nq))) THEN
            XI(1)=0.0d0
            XI(2)=0.5d0
            XI(3)=1.0d0
          ELSE IF((NXQ(-3,1,NXQ(-1,1,nq)).EQ.NWQ(1,nq)).OR.
     '      (NXQ(-1,1,NXQ(-3,1,nq)).EQ.NWQ(1,nq))) THEN
            XI(1)=1.0d0
            XI(2)=0.5d0
            XI(3)=1.0d0
          ELSE IF(NXQ(3,1,NXQ(2,1,NXQ(1,1,nq))).EQ.NWQ(1,nq)) THEN
            XI(1)=0.0d0
            XI(2)=0.0d0
            XI(3)=0.0d0
          ELSE IF(NXQ(3,1,NXQ(2,1,NXQ(-1,1,nq))).EQ.NWQ(1,nq)) THEN
            XI(1)=1.0d0
            XI(2)=0.0d0
            XI(3)=0.0d0
          ELSE IF(NXQ(3,1,NXQ(-2,1,NXQ(1,1,nq))).EQ.NWQ(1,nq)) THEN
            XI(1)=0.0d0
            XI(2)=1.0d0
            XI(3)=0.0d0
          ELSE IF(NXQ(3,1,NXQ(-2,1,NXQ(-1,1,nq))).EQ.NWQ(1,nq)) THEN
            XI(1)=1.0d0
            XI(2)=1.0d0
            XI(3)=0.0d0
          ELSE IF(NXQ(-3,1,NXQ(2,1,NXQ(1,1,nq))).EQ.NWQ(1,nq)) THEN
            XI(1)=0.0d0
            XI(2)=0.0d0
            XI(3)=1.0d0
          ELSE IF(NXQ(-3,1,NXQ(2,1,NXQ(-1,1,nq))).EQ.NWQ(1,nq)) THEN
            XI(1)=1.0d0
            XI(2)=0.0d0
            XI(3)=1.0d0
          ELSE IF(NXQ(-3,1,NXQ(-2,1,NXQ(1,1,nq))).EQ.NWQ(1,nq)) THEN
            XI(1)=0.0d0
            XI(2)=1.0d0
            XI(3)=1.0d0
          ELSE IF(NXQ(-3,1,NXQ(-2,1,NXQ(-1,1,nq))).EQ.NWQ(1,nq)) THEN
            XI(1)=1.0d0
            XI(2)=1.0d0
            XI(3)=1.0d0
          ELSE
            XI(1)=0.0d0
            XI(2)=0.0d0
            XI(3)=0.0d0
          ENDIF !options
        ENDIF !xi dirns
      ENDIF !internal/external

      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(CALC_GRID_XI_1)
        WRITE(OP_STRING,'('' nq,xi(1..3) '',I8,3F8.4)')
     '    nq,XI(1),XI(2),XI(3)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(CALC_GRID_XI_1)
      ENDIF

      CALL EXITS('CALC_GRID_XI')
      RETURN
 9999 CALL ERRORS('CALC_GRID_XI',ERROR)
      CALL EXITS('CALC_GRID_XI')
      RETURN 1
      END



