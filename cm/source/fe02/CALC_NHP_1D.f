      SUBROUTINE CALC_NHP_1D(NHST,NHP,NPNODE,nr,ERROR,*)

C#### Subroutine: CALC_NHP_1D
C###  Description:
C###    CALC_NHP_1D calculates the array NHP.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NHP(NPM,0:NRM),NHST,NPNODE(0:NP_R_M,0:NRM),nr
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nonode,np

      CALL ENTERS('CALC_NHP_1D',*9999)

C     Initialise NHP for current region
      DO np=1,NPM
        NHP(np,nr)=0
      ENDDO !np

      DO nonode=1,NPNODE(0,nr)
        np=NPNODE(nonode,nr)
        NHP(np,nr)=NHST
        NHP(np,0)=NHST
      ENDDO !nonode (np)

      CALL EXITS('CALC_NHP_1D')
      RETURN
 9999 CALL ERRORS('CALC_NHP_1D',ERROR)
      CALL EXITS('CALC_NHP_1D')
      RETURN 1
      END


