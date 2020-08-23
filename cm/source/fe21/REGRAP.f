      SUBROUTINE REGRAP(STRING,ERROR,*)

C#### Subroutine: REGRAP
C###  Description:
C###    REGRAP updates graphics display.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,N3CO
      REAL DELAY
      REAL*8 RFROMC
      LOGICAL CBBREV

      CALL ENTERS('REGRAP',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: refresh graphics
C###  Parameter:      <delay #[0.0]>
C###  Description:
C###    Updates graphics display.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<delay #[0.0]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','REGRAP',ERROR,*9999)
      ELSE
        IF(CBBREV(CO,'DELAY',1,noco+1,NTCO,N3CO)) THEN
          DELAY=REAL(RFROMC(CO(N3CO+1)))
        ELSE
          DELAY=0.0
        ENDIF

        CALL REFRESH_GRAPHICS(DELAY,ERROR,*9999) !Update graphics
      ENDIF

      CALL EXITS('REGRAP')
      RETURN
 9999 CALL ERRORS('REGRAP',ERROR)
      CALL EXITS('REGRAP')
      RETURN 1
      END


