      SUBROUTINE CHECKQ(QLIST,noco,nocoqu,CO,COQU,STRING,*)

C#### Subroutine: CHECKQ
C###  Description:
C###    CHECKQ checks whether command qualifier COQU(noco,nocoqu)
C###    is a member of QLIST.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'mxch.inc'
!     Parameter List
      INTEGER noco,nocoqu
      LOGICAL CELEM
      CHARACTER CO(*)*(*),COQU(25,*)*(*),QLIST*(*),STRING*(MXCH)
!     Local Variables
      CHARACTER ERROR*10

      CALL ENTERS('CHECKQ',*9999)
      IF(.NOT.CELEM(COQU(noco,nocoqu)(1:1),QLIST)) THEN
        IF(COQU(noco,nocoqu)(1:1).EQ.'s'.OR
     '    .COQU(noco,nocoqu)(1:1).EQ.'S') THEN
          CO(noco+1)='?'
          WRITE(OP_STRING,'(//'' Use draw command''/)')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          GO TO 9999
        ELSE
          CO(noco+1)='?'
          STRING='>>Reenter: '//CO(noco)
          GO TO 9999
        ENDIF
      ENDIF

      CALL EXITS('CHECKQ')
      RETURN
 9999 CALL EXITS('CHECKQ')
      RETURN 1
      END


