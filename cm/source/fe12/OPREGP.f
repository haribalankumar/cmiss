      SUBROUTINE OPREGP(REG_PARAMETER,ERROR,*)

C#### Subroutine: OPREGP
C###  Description:
C###    OPREGI outputs regularisation parameters.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'mach00.inc'
!     Parameter List
!     Parameter List
      REAL*8 REG_PARAMETER(0:NTSM)
      CHARACTER ERROR*(*)

!     Local Variables
      INTEGER nt

      CALL ENTERS('OPREGP',*9999)

      WRITE(OP_STRING,'(''#'',A12,A20)')
     '  'Time','Regul. Parameter'
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      DO nt=1,INT(REG_PARAMETER(0))
        WRITE(OP_STRING,'(I12,E16.5)')
     '    nt,REG_PARAMETER(nt)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      ENDDO

      CALL EXITS('OPREGP')
      RETURN
 9999 CALL ERRORS('OPREGP',ERROR)
      CALL EXITS('OPREGP')
      RETURN 1
      END


