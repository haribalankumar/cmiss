      SUBROUTINE MULTLY(STR1,STR2,ERROR,*)

C#### Subroutine: MULTLY
C###  Description:
C###    MULTLY replaces STR1 with STR1*STR2 where * implies inner
C###    product.
C**** RL1(norl) contains NTRL1 REAL*8 numbers from STR1
C**** RL2(norl)    "     NTRL2  "      "      "  STR2
C**** NTLV1 is number of nested levels in  STR1
C**** NTLV2  "    "    "   "      "     "  STR2
C**** NTILV1(nolv) is number of REAL*8 nos at nesting level nolv in STR1
C**** NTILV2(nolv) "    "    "   "    "   "    "      "     "       STR2

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      CHARACTER ERROR*(*),STR1*(*),STR2*(*)
!     Local Variables
      INTEGER i,IEND,j,k,nolv,norl,NTILV1(10),NTILV2(10),NTLV1,NTLV2,
     '  NTRL1,NTRL2
      REAL*8 RL1(100),RL2(100),SUM
      CHARACTER CHARTEMP*20
      LOGICAL SCALR1,SCALR2

      CALL ENTERS('MULTLY',*9999)
      IF(STR1(1:1).EQ.'[') THEN
        SCALR1=.FALSE.
        CALL PARSRA(STR1,10,NTLV1,NTILV1,100,RL1,ERROR,*9999)
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' Numbers in each array: '',10I4)')
     '      (NTILV1(nolv),nolv=1,NTLV1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        DO 100 nolv=NTLV1+1,10
          NTILV1(nolv)=1
 100    CONTINUE
        NTRL1=0
        DO 110 nolv=1,NTLV1
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' RL1: '',10E11.4)')
     '        (RL1(NTRL1+norl),norl=1,NTILV1(nolv))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          NTRL1=NTRL1+NTILV1(nolv)
 110    CONTINUE
      ELSE
        SCALR1=.TRUE.
        NTLV1=1
        CALL PARSRL(STR1,100,NTRL1,RL1,ERROR,*9999)
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' RL1:'',8E10.3)')
     '      (RL1(norl),norl=1,NTRL1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
      ENDIF

      IF(STR2(1:1).EQ.'[') THEN
        SCALR2=.FALSE.
        CALL PARSRA(STR2,10,NTLV2,NTILV2,100,RL2,ERROR,*9999)
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' Numbers in each array: '',10I4)')
     '      (NTILV2(nolv),nolv=1,NTLV2)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        DO 200 nolv=NTLV2+1,10
          NTILV2(nolv)=1
 200    CONTINUE
        NTRL2=0
        DO 210 nolv=1,NTLV2
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' RL2: '',10E11.4)')
     '        (RL2(NTRL2+norl),norl=1,NTILV2(nolv))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          NTRL2=NTRL2+NTILV2(nolv)
 210    CONTINUE
      ELSE
        SCALR2=.TRUE.
        NTLV2=1
        CALL PARSRL(STR2,100,NTRL2,RL2,ERROR,*9999)
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' RL2:'',8E10.3)')
     '      (RL2(norl),norl=1,NTRL2)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
      ENDIF

      IF(NTLV1.EQ.1.AND.NTLV2.EQ.1.AND.
     '  ((SCALR1.AND.SCALR2).OR.(.NOT.SCALR1.AND..NOT.SCALR2))) THEN
        IF(NTRL1.EQ.NTRL2) THEN
          SUM=0.d0
          DO norl=1,NTRL1
            SUM=SUM+RL1(norl)*RL2(norl)
          ENDDO
C          STR1=CFROMR(SUM,'E15.8')
          WRITE(STR1,'(E15.8)') SUM
        ELSE
          OP_STRING(1)='Incompatible strings in multiply'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
      ELSE IF(SCALR1) THEN
        STR1(1:1)='['
        IEND=1
        DO i=1,NTILV2(1)
          IF(NTLV2.GT.1) THEN
            STR1=STR1(1:IEND)//'['
            IEND=IEND+1
          ENDIF
          DO j=1,NTILV2(2)
            SUM=RL1(1)*RL2((i-1)*NTILV2(2)+j)
            WRITE(CHARTEMP,'(E15.8)') SUM
            STR1=STR1(1:IEND)//CHARTEMP//','
C news MPN 28Jul97: CHARTEMP is now 20 chars long (not 15)
            IEND=IEND+21
C old            IEND=IEND+16
          ENDDO
          IF(NTILV2(2).GT.1) THEN
            STR1=STR1(1:IEND-1)//'],'
            IEND=IEND+1
          ENDIF
        ENDDO
        IF(NTILV2(1).EQ.1) THEN
          STR1(IEND:IEND)=' '
        ELSE
          STR1(IEND:IEND)=']'
        ENDIF
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,*) 'string1 scalar ',STR1(1:IEND)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
      ELSE IF(SCALR2) THEN
        STR1(1:1)='['
        IEND=1
        DO i=1,NTILV1(1)
          IF(NTLV1.GT.1) THEN
            STR1=STR1(1:IEND)//'['
            IEND=IEND+1
          ENDIF
          DO j=1,NTILV1(2)
            SUM=RL2(1)*RL1((i-1)*NTILV1(2)+j)
            WRITE(CHARTEMP,'(E15.8)') SUM
            STR1=STR1(1:IEND)//CHARTEMP//','
C news MPN 28Jul97: CHARTEMP is now 20 chars long (not 15)
            IEND=IEND+21
C old            IEND=IEND+16
          ENDDO
          STR1=STR1(1:IEND-1)//'],'
          IEND=IEND+1
        ENDDO
        IF(NTILV1(1).EQ.1) THEN
          STR1(IEND:IEND)=' '
        ELSE
          STR1(IEND:IEND)=']'
        ENDIF
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,*) 'string2 scalar ',STR1(1:IEND)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
      ELSE IF(NTILV1(2).EQ.NTILV2(1)) THEN
        STR1(1:1)='['
        IEND=1
        DO i=1,NTILV1(1)
          IF(NTLV1.GT.1) THEN
            STR1=STR1(1:IEND)//'['
            IEND=IEND+1
          ENDIF
          DO j=1,NTILV2(2)
            SUM=0.d0
            DO k=1,NTILV1(2)
              SUM=SUM+RL1((i-1)*NTILV1(2)+k)*RL2((k-1)*NTILV2(2)+j)
            ENDDO
            WRITE(CHARTEMP,'(E15.8)') SUM
            STR1=STR1(1:IEND)//CHARTEMP//','
C news MPN 28Jul97: CHARTEMP is now 20 chars long (not 15)
            IEND=IEND+21
C old            IEND=IEND+16
          ENDDO
          STR1=STR1(1:IEND-1)//'],'
          IEND=IEND+1
        ENDDO
        IF(NTILV1(1).EQ.1) THEN
          STR1(IEND:IEND)=' '
        ELSE
          STR1(IEND:IEND)=']'
        ENDIF
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          OP_STRING(1)=STR1(1:IEND)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
      ELSE
        OP_STRING(1)='Incompatible strings in multiply'
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('MULTLY')
      RETURN
 9999 CALL ERRORS('MULTLY',ERROR)
      CALL EXITS('MULTLY')
      RETURN 1
      END


