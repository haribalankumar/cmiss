      SUBROUTINE SHOWCO(ERROR,*)

C#### Subroutine: SHOWCO
C###  Description:
C###    SHOWCO outputs the parsed commands, parameters,
C###    and their associated qualifiers.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'trac00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IBEG,IEND,nocoqu

      CALL ENTERS('SHOWCO',*9999)
      IF(TR04) THEN
        WRITE(OP_STRING,'(A32,I4)') ' Number of commands: ',NTCO
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF
      DO noco=1,NTCO
        CALL STRING_TRIM(CO(noco),IBEG,IEND)
        IF(TR04) THEN
          WRITE(OP_STRING,'(A32,A,A1)') ' Command: "',
     '      CO(noco)(IBEG:IEND),'"'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(A32,I4)')
     '      ' Number of command qualifiers: ',NTCOQU(noco)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        DO nocoqu=1,NTCOQU(noco)
          CALL STRING_TRIM(COQU(noco,nocoqu),IBEG,IEND)
          IF(TR04) THEN
            WRITE(OP_STRING,'(A40,A,A1)') ' Command qualifier: "',
     '        COQU(noco,nocoqu)(IBEG:IEND),'"'
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO
      ENDDO

      CALL EXITS('SHOWCO')
      RETURN
 9999 CALL ERRORS('SHOWCO',ERROR)
      CALL EXITS('SHOWCO')
      RETURN 1
      END


