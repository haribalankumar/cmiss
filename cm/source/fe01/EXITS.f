      SUBROUTINE EXITS(NAME)

C#### Subroutine: EXITS
C###  Description:
C###    EXITS traces exit from a subprogram recording the level of
C###    nesting, the exit time, and writing the subprogram name to a
C###    trace file.  If diagnostics are on (DIAGNO=.TRUE.) DOP is
C###    switched off if .NOT.ALLSUB & NAME=SUBNAM.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'diag00.cmn'
      INCLUDE 'trac00.cmn'
!     Parameter List
      CHARACTER NAME*(*)
!     Local Variables
      INTEGER IBEG,IBEG1,IEND,IEND1,n1lv,n1sb
      REAL*8 TM,TMEX,TMTL,VTIME
      CHARACTER COLV*10,ERROR*10

      IF(DIAGNO) THEN
        CALL STRING_TRIM(NAME,IBEG1,IEND1)
        IF(ALLSUB) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' *** Exit  '',A)') NAME(IBEG1:IEND1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ELSE !IF(.NOT.ALLSUB) THEN
          IF(DOP_STACK(NUM_STACK)) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' *** Exit  '',A)') NAME(IBEG1:IEND1)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          IF(NUM_STACK.GT.1) NUM_STACK=NUM_STACK-1
          DOP=DOP_STACK(NUM_STACK)
        ENDIF
      ENDIF

      IF(TR01) THEN
        IF(NAME.EQ.TRSB) THEN
          TR02=.FALSE.
        ENDIF
        TR01=.FALSE.
        IF((nolv.GT.0).AND.(nolv.LE.NXLV)) THEN
          TM=VTIME()
          TMEX=TM-TMST
          TMTL=TMEX-TMEN(nolv)
          TMEL(nolv)=TMTL-TMEL(nolv)
          nosb=NOSBLV(nolv)
          IF(nosb.LE.NXSB) THEN
            IF(SB(nosb).NE.NAME) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'(T32,''!Error in EXIT: '''//
     '          ',A,'' entered, '',A,'' exited.'')')
     '          SB(nosb),NAME
              CALL WRITES(IOTR,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
              RETURN
            ENDIF
          ENDIF
          TMELSM(nosb)=TMELSM(nosb)+TMEL(nolv)
          TMTLSM(nosb)=TMTLSM(nosb)+TMTL
          IF(TR02) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(COLV,'(I10)') nolv
            CALL STRING_TRIM(COLV,IBEG,IEND)
            WRITE(OP_STRING,'(1X,3(F10.3),'//COLV(IBEG:IEND)//
     '        '('' <''),A)')
     '        TMEX,TMTL,TMEL(nolv),NAME
            CALL WRITES(IOTR,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          DO n1lv=1,nolv-1
            TMEL(n1lv)=TMEL(n1lv)+TMEL(nolv)
          ENDDO
        ENDIF
        nolv=nolv-1
        IF(nolv.EQ.0) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'(/''                 Average times:''/'//
     '      '''     Calls:    Total:   Actual:  Subprogram:'')')
          CALL WRITES(IOTR,OP_STRING,ERROR,*9999)
          DO n1sb=1,NTSB
            WRITE(OP_STRING,'(1X,I10,2(F10.3),2X,A)')
     '        NOSBSM(n1sb),TMTLSM(n1sb)/NOSBSM(n1sb),
     '        TMELSM(n1sb)/NOSBSM(n1sb),SB(n1sb)
            CALL WRITES(IOTR,OP_STRING,ERROR,*9999)
          ENDDO
CC$        call mp_unsetlock()
        ENDIF
        IF(NAME.EQ.TRSB) THEN
          TR02=.FALSE.
        ENDIF
        TR01=.TRUE.
      ENDIF

 9999 RETURN
      END


