      SUBROUTINE TRACE(STRING,noco,NTCO,CO,ERROR,*)

C#### Subroutine: TRACE
C###  Description:
C###    TRACE turns the trace on or off and diagnostic o/p on
C###    (DOP=.true.) at entry to subroutine TRSB.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'mxch.inc'
      INCLUDE 'trac00.cmn'
!     Parameter List
      INTEGER noco,NTCO
      CHARACTER CO(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND
      LOGICAL ABBREV
      CHARACTER COUP*(MXCH)

      CALL ENTERS('TRACE',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
C#### Command:  set trace commands on
C###  Parameter:          subroutines on <SUBROUTINE_NAME>
C###  Description:
C###
        IF(TR03) THEN
          WRITE(OP_STRING,'(1X,A)') BLANK(1:IEND)
     '      //' commands off'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ELSE
          WRITE(OP_STRING,'(1X,A)') STRING(1:IEND)
     '      //' commands on'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        IF(TR01) THEN
          WRITE(OP_STRING,'(1X,A)') BLANK(1:IEND)
     '      //' subroutines off <SUBROUTINE_NAME>'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ELSE
          WRITE(OP_STRING,'(1X,A)') BLANK(1:IEND)
     '      //' subroutines on <SUBROUTINE_NAME>'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
      ELSE
        IF(ABBREV(CO(noco+1),'COMMANDS',1)) THEN
          IF(ABBREV(CO(noco+2),'ON',2)) THEN
            TR03=.TRUE.
            TR04=.TRUE.
          ELSE IF(ABBREV(CO(noco+2),'OFF',2)) THEN
            TR03=.FALSE.
            TR04=.FALSE.
          ENDIF
        ELSE IF(ABBREV(CO(noco+1),'SUBROUTINES',1)) THEN
          IF(ABBREV(CO(noco+2),'ON',2)) THEN
            TR01=.TRUE.
            TR02=.TRUE.
          ELSE IF(ABBREV(CO(noco+2),'OFF',2)) THEN
            TR01=.FALSE.
            TR02=.FALSE.
          ENDIF
          IF(NTCO.EQ.noco+3) THEN
C            COUP=CUPPER(CO(noco+3))
      CALL CUPPER(CO(noco+3),COUP)
            IF(TR01) THEN
              TRSB=COUP(1:10)
            ELSE
              DOP=.FALSE.
            ENDIF
          ELSE
            DOP=.FALSE.
          ENDIF
        ENDIF
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,*) ' TR01=',TR01,' TRO2=',TR02,
     '                  ' TR03=',TR03,' TRO4=',TR04
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
      ENDIF

      CALL EXITS('TRACE')
      RETURN
 9999 CALL ERRORS('TRACE',ERROR)
      CALL EXITS('TRACE')
      RETURN 1
      END


