      SUBROUTINE LIASSI(STRING,noco,CO,ERROR,*)

C#### Subroutine: LIASSI
C###  Description:
C###    LIASSI lists the currently assigned user variables, their
C###    values,and whether they are currently active.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'mxch.inc'
      INCLUDE 'user00.cmn'
!     Parameter List
      INTEGER noco
      CHARACTER CO(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,nous

      CALL ENTERS('LIASSI',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
C#### Command: list assign user
C###  Description:
C###    lists the currently assigned user variables, their
C###    values,and whether they are currently active.

        WRITE(OP_STRING,'(1X,A)') STRING(1:IEND)//' user'
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ELSE
        DO nous=1,NTUS
          IF(ACUS(nous)) THEN !user assigned string is active
            WRITE(OP_STRING,'('' User assigned string '',A,'' = '',A,'
     '        //''' is active'')')
     '        USID(nous)(:LNUSID(nous)),US(nous)(:LNUS(nous))
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE                !user assigned string is inactive
            WRITE(OP_STRING,'('' User assigned string '',A,'' = '',A,'
     '        //''' is inactive'')')
     '        USID(nous)(:LNUSID(nous)),US(nous)(:LNUS(nous))
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO !nous
      ENDIF

      CALL EXITS('LIASSI')
      RETURN
 9999 CALL ERRORS('LIASSI',ERROR)
      CALL EXITS('LIASSI')
      RETURN 1
      END


