      SUBROUTINE DEASSI(STRING,noco,CO,ERROR,*)

C#### Subroutine: DEASSI
C###  Description:
C###    DEASSI deassigns a previously assigned user variable.

C Note: May want ALL option on this - then set NTUS=0

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'mxch.inc'
      INCLUDE 'user00.cmn'
!     Parameter List
      INTEGER noco
      CHARACTER CO(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,nous

      CALL ENTERS('DEASSI',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
C news MPN 27Jul97: new style help
C---------------------------------------------------------------------

C#### Command: deassign' STRING_NAME'
C###  Description:
C###    Deassigns the previously assigned user defined variable
C###    STRING_NAME
        OP_STRING(1)=STRING(1:IEND)//' STRING_NAME'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
C new end
C old        WRITE(OP_STRING,'(1X,A)') STRING(1:IEND)
C old     '    //' STRING_NAME'
C old        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      ELSE
        CALL STRING_TRIM(CO(noco+1),IBEG,IEND)
        DO nous=1,NTUS
          IF(CO(noco+1)(IBEG:IEND).EQ.USID(nous)(:LNUSID(nous))) THEN
            ACUS(nous)=.FALSE.
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('DEASSI')
      RETURN
 9999 CALL ERRORS('DEASSI',ERROR)
      CALL EXITS('DEASSI')
      RETURN 1
      END


