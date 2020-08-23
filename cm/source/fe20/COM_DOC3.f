      SUBROUTINE COM_DOC3(CO,OPTION,STRING,ERROR,*)

C#### Subroutine: COM_DOC3
C###  Description:
C###    Used when documenting commands. Called from FE20.

!     Parameter List
      CHARACTER CO(*)*132,ERROR*(*),OPTION*(*),STRING*(*)
!     Local Variables
      INTEGER IBEG1,IEND1,IBEG2,IEND2,IBEG3,IEND3

      CALL ENTERS('COM_DOC3',*9999)

      CO(3)=OPTION
      CO(4)='?'
      CALL STRING_TRIM(CO(1),IBEG1,IEND1)
      CALL STRING_TRIM(CO(2),IBEG2,IEND2)
      CALL STRING_TRIM(CO(3),IBEG3,IEND3)
      STRING=' '//CO(1)(IBEG1:IEND1)
     '     //' '//CO(2)(IBEG2:IEND2)
     '     //' '//CO(3)(IBEG3:IEND3)

      CALL EXITS('COM_DOC3')
      RETURN
 9999 CALL ERRORS('COM_DOC3',ERROR)
      CALL EXITS('COM_DOC3')
      RETURN 1
      END



