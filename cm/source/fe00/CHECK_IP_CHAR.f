      SUBROUTINE CHECK_IP_CHAR(irec,N1char,N2char,POS_dollar,
     '  G_FORMAT,LINE,REVERT1,ERROR,*,*,*)

C#### Subroutine: CHECK_IP_CHAR
C###  Description:
C###    CHECK_IP_CHAR checks for * or ! in ip files.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER irec,N1char,N2char,POS_dollar
      CHARACTER G_FORMAT*(*),LINE*(*),ERROR*(*)
      LOGICAL REVERT1
!     Local Variables

      CALL ENTERS('CHECK_IP_CHAR',*9999)
      IF(LINE(1:1).EQ.'!') THEN   !This line is a comment
        IF(LINE(2:2).EQ.'!') THEN !Do not print
        ELSE                      !Print comment line
          OP_STRING(1)=LINE
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        irec=irec+1
        GO TO 9993 !go to 3 to read next line

      ELSE IF(LINE(1:1).EQ.'*') THEN
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' Found *:Revert to prompted input'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' File positioned at irec='',I4)')irec
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        REVERT1=.TRUE.
C      Replace $ in format for prompted read
        G_FORMAT(POS_dollar+2:)=G_FORMAT(POS_dollar:)
        G_FORMAT(POS_dollar:POS_dollar+1)='$,'
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' G_FORMAT: '',A)') G_FORMAT(1:80)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        IF(N1char.GT.0) N1char=N1char+2
        IF(N2char.GT.0) N2char=N2char+2
        GO TO 9995 !go to 5 to use prompted input
      ENDIF !line(1:1)=! or *

      CALL EXITS('CHECK_IP_CHAR')
      RETURN

 9993 CALL EXITS('CHECK_IP_CHAR')
      RETURN 1

 9995 CALL EXITS('CHECK_IP_CHAR')
      RETURN 2

 9999 CALL EXITS('CHECK_IP_CHAR')
      RETURN 3
      END


