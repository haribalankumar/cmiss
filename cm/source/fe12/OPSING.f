      SUBROUTINE OPSING(ERROR,*)

C#### Subroutine: OPSING
C###  Description:
C###    OPSING outputs corner node data (for BE problems).

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'sing00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i

      CALL ENTERS('OPSING',*9999)
      WRITE(OP_STRING,
     '  '(''  The number of physical singularity locations is '',I3)')
     '  NPSING(0)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      DO i=1,NPSING(0)
        WRITE(OP_STRING,'(''   Physical singularity '',I3)') i
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '(''   Located at node '',I3,'' and in element '',I3)')
     '    NPSING(i),NESING(i)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDDO
      CALL EXITS('OPSING')

      RETURN
 9999 CALL ERRORS('OPSING',ERROR)
      CALL EXITS('OPSING')
      RETURN 1
      END


