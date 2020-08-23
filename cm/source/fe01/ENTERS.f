      SUBROUTINE ENTERS(NAME,*)

C#### Subroutine: ENTERS
C###  Description:
C###    ENTERS traces entry to a subprogram recording the level of
C###    nesting, the entry time, and writing the subprogram name to a
C###    trace file.  Diagnostic o/p is turned on if DIAGNO=.TRUE. and
C###    ALLSUB=.TRUE. or ALLSUB=.FALSE. and NAME=SUBNAM (a subroutine
C###    name).
C**** TRSB         is subprogram name which turns trace on
C**** IOTR         is trace file number (diagnostic o/p file)
C**** TR01         is true if basic tracing on
C**** TR02         is true if full tracing on
C**** NOLV         is current level number
C**** NTLV         is current total number of levels called
C**** NXLV         is maximum number of levels which can be traced
C**** NOSB         is integer number for current subroutine
C**** NTSB         is current total number of subroutines called
C**** NXSB         is maximum no. of subroutines which can be traced
C**** NOSBLV(nolv) is NOSB number at level NOLV
C**** NOSBSM(nosb) is no. of times subroutine no. NOSB has been called
C**** SB(nosb)     is subroutine name for subroutine no. NOSB
C**** TM           is current time
C**** TMST         is trace start time for entry to level 1
C**** TMEL(nolv)   is elapsed time at level NOLV (ie =0 @ entry to NOLV)
C**** TMELSM(nosb) is sum of elapsed times spent in subroutine no. NOSB
C**** TMEN(nolv)   is elapsed time (from start of trace) to level NOLV
C**** TMTLSM(nosb) is total time spent in NOSB & all subs below it

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
C      INCLUDE 'cmiss$reference:ctrl00.cmn'
      INCLUDE 'diag00.cmn'
      INCLUDE 'mxch.inc'
      INCLUDE 'trac00.cmn'
!     Parameter List
      CHARACTER NAME*(*)
!     Local Variables
      INTEGER IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,
     '  N1SB,NSUB
      REAL*8 TM,VTIME
      CHARACTER COLV*10,C1*(MXCH),ERROR*10
c      EXTERNAL CTRLC_AST

      IF(DIAGNO) THEN
        CALL STRING_TRIM(NAME,IBEG1,IEND1)
        IF(ALLSUB) THEN !turn diagnostics on in all subroutines
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' *** Enter '',A)') NAME(IBEG1:IEND1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ELSE !IF(.NOT.ALLSUB) THEN !diags on in selected subrs
          NUM_STACK=NUM_STACK+1
          DOP_STACK(NUM_STACK)=.FALSE.
          DO NSUB=1,NT_SUB
            CALL STRING_TRIM(SUBNAM(NSUB),IBEG2,IEND2)
            CALL CUPPER(NAME(IBEG1:IEND1),C1)
C          IF(CUPPER(NAME(IBEG1:IEND1)).EQ.SUBNAM(NSUB)(IBEG2:IEND2))
            IF(C1(IBEG1:IEND1).EQ.SUBNAM(NSUB)(IBEG2:IEND2))
     '        DOP_STACK(NUM_STACK)=.TRUE.
          ENDDO
          IF(FROMSUB) THEN
            IF(DOP_STACK(NUM_STACK-1)) DOP_STACK(NUM_STACK)=.TRUE.
          ENDIF
          IF(DOP_STACK(NUM_STACK)) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' *** Enter '',A)') NAME(IBEG1:IEND1)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ELSE IF(DOP_STACK(NUM_STACK-1)) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' *** Calls '',A)') NAME(IBEG1:IEND1)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          DOP=DOP_STACK(NUM_STACK)
        ENDIF
      ENDIF

      IF(TR01) THEN
        IF(NAME.EQ.TRSB) THEN
          TR02=.TRUE.
        ENDIF
        TR01=.FALSE.
        NOLV=NOLV+1
        IF(NOLV.GT.NTLV) THEN
          NTLV=NOLV
        ENDIF
        IF((NOLV.GT.0).AND.(NOLV.LE.NXLV)) THEN
          TM=VTIME()
          IF(NOLV.EQ.1) THEN
            TMST=TM
            IF(TR02) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'(/''      Time:    Calls:    Level:'//
     '          ' >Subprogram entered'''//
     '          '/''      Time:    Total:   Actual:'//
     '          ' <Subprogram exited'')')
              CALL WRITES(IOTR,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
          ENDIF
          TMEN(NOLV)=TM-TMST
          TMEL(NOLV)=0.0D0
          DO N1SB=1,NTSB
            IF(SB(N1SB).EQ.NAME) THEN
              NOSB=N1SB
              GOTO 2
            ENDIF
          ENDDO
          IF(NTSB.LT.NXSB) THEN
            NTSB=NTSB+1
            SB(NTSB)=NAME
            NOSB=NTSB
          ELSE
            NTSB=NXSB+1
            NOSB=NXSB+1
          ENDIF
    2     NOSBLV(NOLV)=NOSB
          NOSBSM(NOSB)=NOSBSM(NOSB)+1
          IF(TR02) THEN
C            COLV=CFROMI(NOLV,'(I10)')
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(COLV,'(I10)') NOLV
            CALL STRING_TRIM(COLV,IBEG,IEND)
            WRITE(OP_STRING,'(1X,F10.3,I10,I10,'//COLV(IBEG:IEND)//
     '        '('' >''),A)')
     '        TMEN(NOLV),NOSBSM(NOSB),NOLV,NAME
            CALL WRITES(IOTR,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
        ENDIF
        TR01=.TRUE.
      ENDIF

      RETURN
 9999 RETURN 1
      END



