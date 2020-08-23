      SUBROUTINE USER52(ERROR,*)

C#### Subroutine: USER52
C###  Description:
C###    USER52 returns the derivatives of the user-defined strain energy
C###    function of principal extension ratios L1, L2, and L3 to
C###    subroutine ENERGY at current Gauss point.

      IMPLICIT NONE
!     Parameter List
      CHARACTER ERROR*(*)

      CALL ENTERS('USER52',*9999)

C     DW(1)=
C     DW(2)=
C     DW(3)=

      CALL EXITS('USER52')
      RETURN
 9999 CALL ERRORS('USER52',ERROR)
      CALL EXITS('USER52')
      RETURN 1
      END


