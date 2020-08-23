      SUBROUTINE SETPR(NTCO,CO,NTCOQU,COQU,ERROR,*)

C#### Subroutine: SETPR
C###  Description:
C###    SETPR sets the prompt string PR.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbpr00.cmn'
!     Parameter List
      INTEGER NTCO,NTCOQU(*)
      CHARACTER CO(*)*(*),COQU(25,*)*(*),ERROR*(*)
!     Local Variables
      INTEGER IBEG,IEND,noco,nocoqu,prend

      CALL ENTERS('SETPR',*9999)

      prend=0
      DO 2 noco=1,NTCO
        CALL STRING_TRIM(CO(noco),IBEG,IEND)
        CALL APPENDC(prend,CO(noco)(IBEG:IEND),PR)
        DO 1 nocoqu=1,NTCOQU(noco)
          CALL APPENDC(prend,';',PR)
          CALL STRING_TRIM(COQU(noco,nocoqu),IBEG,IEND)
          CALL APPENDC(prend,COQU(noco,nocoqu)(IBEG:IEND),PR)
 1      CONTINUE
        CALL APPENDC(prend,' ',PR)
 2    CONTINUE
      IF(prend.EQ.0) prend=1
      PR(prend:prend)=CHAR(0)

      CALL EXITS('SETPR')
      RETURN
 9999 CALL ERRORS('SETPR',ERROR)
      CALL EXITS('SETPR')
      RETURN 1
      END
