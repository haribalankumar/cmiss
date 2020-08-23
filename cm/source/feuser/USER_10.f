      SUBROUTINE USER_10(NPLIST,NPNODE,nr,NVCB,NVCNODE,SIMTIME,
     '  VALVES_SWITCHED,ERROR,*)

C#### Subroutine: USER_10
C###  Description:
C###    Operates valves for fe60.f problems

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'fluid00.cmn'
      INCLUDE 'voro00.inc'
!     Parameter List
      INTEGER NPLIST(0:NPM),NPNODE(0:NP_R_M,0:NRM),nr,NVCB(-1:3,NVCBM),
     '  NVCNODE(2,NP_R_M)
      REAL*8 SIMTIME
      LOGICAL VALVES_SWITCHED
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 VALVE_FLAG
      INTEGER np,nonode,nvc,n

      CALL ENTERS('USER_10',*9999)

      DO nonode=1,NPNODE(0,nr)
        np=NPNODE(nonode,nr)
        NPLIST(np)=nonode
      ENDDO

      VALVE_FLAG=DSIN(2.d0*PI*SIMTIME/HEART_PERIOD)
      IF(VALVE_FLAG.GT.0.d0) THEN
        WRITE(*,'(/,'' Diastole '',D13.6)') VALVE_FLAG
        IF(SYSTOLE) THEN
          VALVES_SWITCHED=.TRUE.
        ENDIF
        SYSTOLE=.FALSE.
      ELSE
        WRITE(*,'(/,'' Systole '',D13.6)') VALVE_FLAG
        IF(.NOT.SYSTOLE) THEN
          VALVES_SWITCHED=.TRUE.
        ENDIF
        SYSTOLE=.TRUE.
      ENDIF

      IF(SYSTOLE) THEN
        NUM_OUTLETS=1
        OUTLET_FIXDNODES(1)=2700
      ELSE !diastole
        NUM_OUTLETS=1
        OUTLET_FIXDNODES(1)=3330
      ENDIF

      IF(VALVES_SWITCHED) THEN
        IF(SYSTOLE) THEN
          DO n=1,DIASTOLE_OUTLET(0)
            np=DIASTOLE_OUTLET(n)
            nonode=NPLIST(np)
            nvc=NVCNODE(MAP,nonode)
            CALL ASSERT(NVCNODE(TYPE,nonode).EQ.BOUNDARY,'>>Non '//
     '        'boundary node found in diastole outlet',ERROR,*9999)
            CALL ASSERT(NVCB(BCTYPE,nvc).EQ.OUTLET,'>> Diastole is '//
     '        'not in outlet mode',ERROR,*9999)
            CALL ASSERT(NVCB(0,nvc).EQ.1,'>> Diastole has '//
     '        'more than one adjoining node',ERROR,*9999)
            NVCB(BCTYPE,nvc)=WALL
          ENDDO
          DO n=1,SYSTOLE_OUTLET(0)
            np=SYSTOLE_OUTLET(n)
            nonode=NPLIST(np)
            nvc=NVCNODE(MAP,nonode)
            CALL ASSERT(NVCNODE(TYPE,nonode).EQ.BOUNDARY,'>>Non '//
     '        'boundary node found in systole outlet',ERROR,*9999)
            CALL ASSERT(NVCB(BCTYPE,nvc).EQ.WALL,'>> systole is '//
     '        'not in wall mode',ERROR,*9999)
            CALL ASSERT(NVCB(0,nvc).EQ.1,'>> Systole has '//
     '        'more than one adjoining node',ERROR,*9999)
            NVCB(BCTYPE,nvc)=OUTLET
          ENDDO
        ELSE !Diastole
          DO n=1,SYSTOLE_OUTLET(0)
            np=SYSTOLE_OUTLET(n)
            nonode=NPLIST(np)
            nvc=NVCNODE(MAP,nonode)
            CALL ASSERT(NVCNODE(TYPE,nonode).EQ.BOUNDARY,'>>Non '//
     '        'boundary node found in systole outlet',ERROR,*9999)
            CALL ASSERT(NVCB(BCTYPE,nvc).EQ.OUTLET,'>> systole is '//
     '        'not in outlet mode',ERROR,*9999)
            CALL ASSERT(NVCB(0,nvc).EQ.1,'>> Systole has '//
     '        'more than one adjoining node',ERROR,*9999)
            NVCB(BCTYPE,nvc)=WALL
          ENDDO
          DO n=1,DIASTOLE_OUTLET(0)
            np=DIASTOLE_OUTLET(n)
            nonode=NPLIST(np)
            nvc=NVCNODE(MAP,nonode)
            CALL ASSERT(NVCNODE(TYPE,nonode).EQ.BOUNDARY,'>>Non '//
     '        'boundary node found in diastole outlet',ERROR,*9999)
            CALL ASSERT(NVCB(BCTYPE,nvc).EQ.WALL,'>> Diastole is '//
     '        'not in wall mode',ERROR,*9999)
            CALL ASSERT(NVCB(0,nvc).EQ.1,'>> Diastole has '//
     '        'more than one adjoining node',ERROR,*9999)
            NVCB(BCTYPE,nvc)=OUTLET
          ENDDO
        ENDIF
      ENDIF

      CALL EXITS('USER_10')
      RETURN
 9999 CALL ERRORS('USER_10',ERROR)
      CALL EXITS('USER_10')
      RETURN 1
      END

