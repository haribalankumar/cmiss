      SUBROUTINE CHFOCU(STRING,ERROR,*)

C#### Subroutine: CHFOCU
C###  Description:
C###    CHFOCU allows user to redefine position of focus.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b14.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND
      REAL*8 RFROMC

      CALL ENTERS('CHFOCU',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM change focus NEW_FOCUS_VALUE
C###  Description:
C###    Change the focus of the prolate spheroidal coordinate system.

        OP_STRING(1)=STRING(1:IEND)//' NEW_FOCUS_VALUE'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','CHFOCU',ERROR,*9999)
      ELSE
        FOCUS=RFROMC(CO(noco+1))
      ENDIF

      CALL EXITS('CHFOCU')
      RETURN
 9999 CALL ERRORS('CHFOCU',ERROR)
      CALL EXITS('CHFOCU')
      RETURN 1
      END


