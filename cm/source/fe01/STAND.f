      SUBROUTINE STAND(END,STRING,ERROR,*)

C#### Subroutine: STAND
C###  Description:
C###    STAND checks the command string for the functions that may be
C###    invoked within any environment.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'user00.cmn'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(MXCH)
      LOGICAL END
!     Local Variables
      INTEGER i,IBEG,IBEG2,IBEG3,IBEG4,IEND,IEND2,IEND3,IEND4,IFROMC,
     '  nous,NTITLV(20),NTLV
      REAL*8 R,RA(100),RFROMC
      LOGICAL ABBREV

      CALL ENTERS('STAND',*9999)
! PJH 21Jul96 Now handled at root level
!      IF(ABBREV(CO(noco),'ASSIGN',2)) THEN
!        CALL ASSIGN(STRING,noco,CO,COQU,ERROR,*9999)
!      ELSE IF(ABBREV(CO(noco),'DEASSIGN',2)) THEN
!        CALL DEASSI(STRING,noco,CO,ERROR,*9999)

C MPN 26Mar2002: I'm pretty sure that the 'DISPLAY' and 'EVALUATE'
C                blocks below are no longer used and can be archived.
      IF(ABBREV(CO(noco),'DISPLAY',2)) THEN
        IF(CO(noco+1).EQ.'?') THEN
          CALL STRING_TRIM(STRING,IBEG,IEND)
          WRITE(OP_STRING,'(1X,A)') STRING(1:IEND)//
     '      ' FILENAME STRING'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ELSE
c         CALL DOCUM(CO(noco+1),'example',CO(noco+2),ERROR,*9999)
        ENDIF

      ELSE IF(ABBREV(CO(noco),'EVALUATE',2)) THEN
        IF(CO(noco+1).EQ.'?') THEN
          CALL STRING_TRIM(STRING,IBEG,IEND)
          WRITE(OP_STRING,'(1X,A)') STRING(1:IEND)//' STRING'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ELSE
          CO(2)=CO(noco+1)
          CO(3)=CO(noco+1)
          CALL PARSTR(CO(3),20,NTLV,NTITLV,100,RA,ERROR,*9999)
          R=RFROMC(CO(3))
          i=IFROMC(CO(3))
          IF(DABS(R-DBLE(i)).LT.1.0D-5) THEN
            WRITE(CO(3),'(I5)') i
          ENDIF
          CALL STRING_TRIM(CO(2),IBEG2,IEND2)
          CALL STRING_TRIM(CO(3),IBEG3,IEND3)
          DO nous=1,NTUS
            IF(US(nous)(1:IEND2-IBEG2+1).EQ.CO(2)(IBEG2:IEND2)) THEN
              US(nous)=CO(3)(IBEG3:IEND3)
              CALL STRING_TRIM(US(nous),IBEG4,IEND4)
              LNUS(nous)=IEND4-IBEG4+1
            ENDIF
          ENDDO
        ENDIF

      ELSE IF(ABBREV(CO(noco),'QUIT',1)) THEN
C#### Command: FEM quit
C###  Description:
C###    End CMISS session
C#### Command: quit
C###  Description:
C###    End CMISS session
        CALL QUIT(END,ERROR)

      ELSE IF(ABBREV(CO(noco),'WINDOW',1)) THEN
        IF(ABBREV(CO(noco-1),'LIST',1)) THEN
          CALL LIWIND(STRING,ERROR,*9999)
        ENDIF

      ELSE IF(ABBREV(CO(noco),'>',1)) THEN
C#### Command: >
C###  Description:
C###    Sets the prompt forward one level.
        CALL SETPR(noco-1,CO,NTCOQU,COQU,ERROR,*9999)

      ELSE IF(ABBREV(CO(noco),'<',1)) THEN
C#### Command: <
C###  Description:
C###    Sets the prompt back one level.
        CALL SETPR(noco-2,CO,NTCOQU,COQU,ERROR,*9999)

      ELSE IF(ABBREV(CO(noco),'<<',2)) THEN
C#### Command: <<
C###  Description:
C###    Sets the prompt to original level.
        CALL SETPR(0,CO,NTCOQU,COQU,ERROR,*9999)

      ELSE IF((INDEX(CO(noco),'+')+INDEX(CO(noco),'-')
     '        +INDEX(CO(noco),'*')+INDEX(CO(noco),'/')
     '        +INDEX(CO(noco),'(')+INDEX(CO(noco),'^')).GT.0) THEN
        STRING=CO(noco)
        CALL PARSTR(STRING,20,NTLV,NTITLV,100,RA,ERROR,*9999)
        CALL STRING_TRIM(STRING,IBEG,IEND)
        WRITE(OP_STRING,'('' ='',A)') STRING(IBEG:IEND)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      ELSE IF(ABBREV(CO(noco),'RETURN',6).OR.noco.GT.NTCO) THEN
C#### Command: return
C###  Description:
C###    Sets the prompt forward one level.
        CALL SETPR(noco-1,CO,NTCOQU,COQU,ERROR,*9999)

      ELSE IF(ABBREV(CO(noco),'??',2)) THEN

      ELSE IF(ABBREV(CO(noco),'!',1)) THEN

      ELSE IF(ABBREV(CO(noco),'#',1)) THEN

      ELSE
        CALL UNKNOW(noco,CO,ERROR,*9999)
      ENDIF

      CALL EXITS('STAND')
      RETURN
 9999 CALL ERRORS('STAND',ERROR)
      CALL EXITS('STAND')
      RETURN 1
      END


