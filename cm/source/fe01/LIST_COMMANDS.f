      SUBROUTINE LIST_COMMANDS(IP,NTCH,OPTION,ERROR,*)

C#### Subroutine: LIST_COMMANDS
C###  Description:
C###    LIST_COMMANDS lists CMISS commands when using help.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'

!     Parameter List
      INTEGER IP,NTCH
      CHARACTER OPTION(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER noch

      CALL ENTERS('LIST_COMMANDS',*9999)

      DO noch=1,NTCH
        IF(IP.EQ.1) THEN
          WRITE(OP_STRING,'(1X,A)') OPTION(noch)
        ELSE IF(IP.EQ.2) THEN
          WRITE(OP_STRING,'(3X,A)') OPTION(noch)
        ELSE IF(IP.EQ.3) THEN
          WRITE(OP_STRING,'(5X,A)') OPTION(noch)
        ENDIF
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDDO

      CALL EXITS('LIST_COMMANDS')
      RETURN
 9999 CALL ERRORS('LIST_COMMANDS',ERROR)
      CALL EXITS('LIST_COMMANDS')
      RETURN 1
      END


