      SUBROUTINE CALCUL(STRING,SUBST,ERROR,*)

C#### Subroutine: CALCUL
C###  Description:
C###    CALCUL calculates arithmetic expressions by first performing
C###    ^ then *,/ and then +,-. Note that the first substring of a
C###    pair is replaced by the result and the subsequent strings
C###    moved back one.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'trac00.cmn'
!     Parameter List
      CHARACTER STRING*(*),SUBST*(*),ERROR*(*)
!     Local Variables
      INTEGER IBEG,IBEG1,IEND,IEND1,n1sb

      CALL ENTERS('CALCUL',*9999)
      CALL SUBSTR(LEN(STRING),STRING,'(+-*/^)',ERROR,*9999)
C *** delete start operator if not negative sign (AAY 6 JUL 90)
      IF(SB(1)(1:1).NE.'-') SB(1)(1:1)=' '
C *** reconcatenate substrings ending in E (ie E+00 should not be split)
      nosb=1
      DO WHILE (nosb.LT.NTSB)
        CALL STRING_TRIM(SB(nosb),IBEG,IEND)
        DO WHILE(SB(nosb)(IEND:IEND).EQ.'E')
          SB(nosb)=SB(nosb)(1:IEND)//SB(nosb+1)
          CALL STRING_TRIM(SB(nosb),IBEG,IEND)
          DO n1sb=nosb+1,NTSB-1
            SB(n1sb)=SB(n1sb+1)
          ENDDO
          NTSB=NTSB-1
        ENDDO
        nosb=nosb+1
      ENDDO
C *** reconcatenate substrings ending in , (ie 3,-2 should not be split)
      nosb=1
      DO WHILE(nosb.LT.NTSB)
        CALL STRING_TRIM(SB(nosb),IBEG,IEND)
        DO WHILE(SB(nosb)(IEND:IEND).EQ.',')
          SB(nosb)=SB(nosb)(1:IEND)//SB(nosb+1)
          CALL STRING_TRIM(SB(nosb),IBEG,IEND)
          DO n1sb=nosb+1,NTSB-1
            SB(n1sb)=SB(n1sb+1)
          ENDDO
          NTSB=NTSB-1
        ENDDO
        nosb=nosb+1
      ENDDO
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        DO nosb=1,NTSB
          WRITE(OP_STRING,'('' SB('',I2,''):'',A70,:,(/A80))')
     '      nosb,SB(nosb)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDDO
CC$      call mp_unsetlock()
      ENDIF

      nosb=2
      DO WHILE(nosb.LE.NTSB)
        IF(SB(nosb)(1:1).EQ.'^') THEN
          CALL EXPON(SB(nosb-1)(2:),SB(nosb)(2:),ERROR,*9999)
          DO n1sb=nosb,NTSB-1
            SB(n1sb)=SB(n1sb+1)
          ENDDO
          NTSB=NTSB-1
        ELSE
          nosb=nosb+1
        ENDIF
      ENDDO

      nosb=2
      DO WHILE(nosb.LE.NTSB)
        IF(SB(nosb)(1:1).EQ.'*') THEN
          IF(SB(nosb-1)(1:1).EQ.'-') THEN !-'ve is sign not operator
            CALL STRING_TRIM(SB(nosb-1),IBEG1,IEND1)
            SB(nosb-1)(1:1)=' '
            SB(nosb-1)(2:)='-'//SB(nosb-1)(2:IEND1)
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              OP_STRING(1)=' Correct for -ve before multiply:'
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' SB('',I2,''):'',A70)')
     '          nosb-1,SB(nosb-1)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
          ENDIF
          CALL MULTLY(SB(nosb-1)(2:),SB(nosb)(2:),ERROR,*9999)
          DO n1sb=nosb,NTSB-1
            SB(n1sb)=SB(n1sb+1)
          ENDDO
          NTSB=NTSB-1
        ELSE IF(SB(nosb)(1:1).EQ.'/') THEN
          IF(SB(nosb-1)(1:1).EQ.'-') THEN !treat -ve as sign not operator
            CALL STRING_TRIM(SB(nosb-1),IBEG1,IEND1)
            SB(nosb-1)(1:1)=' '
            SB(nosb-1)(2:)='-'//SB(nosb-1)(2:IEND1)
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              OP_STRING(1)=' Correct for -ve before divide:'
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' SB('',I2,''):'',A70)')
     '          nosb-1,SB(nosb-1)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
          ENDIF
          CALL DIVIDE(SB(nosb-1)(2:),SB(nosb)(2:),ERROR,*9999)
          DO n1sb=nosb,NTSB-1
            SB(n1sb)=SB(n1sb+1)
          ENDDO
          NTSB=NTSB-1
        ELSE
          nosb=nosb+1
        ENDIF
      ENDDO

      nosb=2
      DO WHILE(nosb.LE.NTSB)
        IF(SB(nosb)(1:1).EQ.'+') THEN
          IF(SB(nosb-1)(1:1).EQ.'-') THEN !treat -ve as sign not operator
            CALL STRING_TRIM(SB(nosb-1),IBEG1,IEND1)
            SB(nosb-1)(1:1)=' '
            SB(nosb-1)(2:)='-'//SB(nosb-1)(2:IEND1)
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              OP_STRING(1)=' Correct for -ve before add:'
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' SB('',I2,''):'',A70)')
     '          nosb-1,SB(nosb-1)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
          ENDIF
          CALL ADD(SB(nosb-1)(2:),SB(nosb)(2:),ERROR,*9999)
          DO n1sb=nosb,NTSB-1
            SB(n1sb)=SB(n1sb+1)
          ENDDO
          NTSB=NTSB-1
        ELSE IF(SB(nosb)(1:1).EQ.'-') THEN
          IF(SB(nosb-1)(1:1).EQ.'-') THEN !treat -ve as sign not operator
            CALL STRING_TRIM(SB(nosb-1),IBEG1,IEND1)
            SB(nosb-1)(1:1)=' '
            SB(nosb-1)(2:)='-'//SB(nosb-1)(2:IEND1)
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              OP_STRING(1)=' Correct for -ve before subtract:'
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' SB('',I2,''):'',A70)')
     '          nosb-1,SB(nosb-1)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
          ENDIF
          CALL SUBTR(SB(nosb-1)(2:),SB(nosb)(2:),ERROR,*9999)
          DO n1sb=nosb,NTSB-1
            SB(n1sb)=SB(n1sb+1)
          ENDDO
          NTSB=NTSB-1
        ENDIF
      ENDDO

      SUBST(1:)=SB(1)(2:)

      CALL EXITS('CALCUL')
      RETURN
 9999 CALL ERRORS('CALCUL',ERROR)
      CALL EXITS('CALCUL')
      RETURN 1
      END


