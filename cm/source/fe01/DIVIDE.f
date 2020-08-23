      SUBROUTINE DIVIDE(STR1,STR2,ERROR,*)

C#### Subroutine: DIVIDE
C###  Description:
C###    DIVIDE replaces STR1 with STR1/STR2.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      CHARACTER ERROR*(*),STR1*(*),STR2*(*)
!     Local Variables
      INTEGER i,IEND,j,nolv,norl,NTILV1(10),NTLV1,NTRL1,NTRL2
      REAL*8 RL1(100),RL2(100),SUM
      CHARACTER CHARTEMP*20
      LOGICAL SCALR1

      CALL ENTERS('DIVIDE',*9999)
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
        DO nolv=NTLV1+1,10
          NTILV1(nolv)=1
        ENDDO
        NTRL1=0
        DO nolv=1,NTLV1
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
        ENDDO
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
        OP_STRING(1)='Incompatible strings in divide'
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        GO TO 9999
      ELSE
        CALL PARSRL(STR2,100,NTRL2,RL2,ERROR,*9999)
      ENDIF

      IF(SCALR1) THEN
        SUM=RL1(1)/RL2(1)
C        STR1=CFROMR(SUM,'E15.8')
        WRITE(STR1,'(E15.8)') SUM
      ELSE
        STR1(1:1)='['
        IEND=1
        DO i=1,NTILV1(1)
          IF(NTLV1.GT.1) THEN
            STR1=STR1(1:IEND)//'['
            IEND=IEND+1
          ENDIF
          DO j=1,NTILV1(2)
            SUM=RL1((i-1)*NTILV1(2)+j)/RL2(1)
            WRITE(CHARTEMP,'(E15.8)') SUM
            STR1=STR1(1:IEND)//CHARTEMP//','
C news MPN 28Jul97: CHARTEMP is now 20 chars long (not 15)
            IEND=IEND+21
C old       IEND=IEND+16
          ENDDO
          IF(NTILV1(2).GT.1) THEN
            STR1=STR1(1:IEND-1)//'],'
            IEND=IEND+1
          ENDIF
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
      ENDIF

      CALL EXITS('DIVIDE')
      RETURN
 9999 CALL ERRORS('DIVIDE',ERROR)
      CALL EXITS('DIVIDE')
      RETURN 1
      END


