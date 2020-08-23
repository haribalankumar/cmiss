      SUBROUTINE USER53(CG,DW,E11,E22,E33,E12,E13,E23,ERROR,*)

C#### Subroutine: USER53
C###  Description:
C###    USER53 returns the derivatives of the user-defined strain energy
C###    function of fibre and transverse physical strains E11...E23
C###    to subroutine ENERGY at current Gauss point.
C###    Double exponential law.
C###    Exponential-dependence on fibre/transverse extension ratios

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      REAL*8 CG(NMM),DW(6),E11,E12,E13,E22,E23,E33
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 CEXPQ1,CEXPQ2,TOL

      CALL ENTERS('USER53',*9999)

      TOL=1.0D-08
      IF((DABS(E11).LT.TOL).AND.(DABS(E22).LT.TOL).AND.
     '  (DABS(E33).LT.TOL).AND.(DABS(E12).LT.TOL).AND.
     '  (DABS(E13).LT.TOL).AND.(DABS(E23).LT.TOL)) THEN
        CEXPQ1=CG(1)
        CEXPQ2=CG(4)
      ELSE
        CEXPQ1=CG(1)*DEXP(CG(2)*E11*E11+CG(3)*(E12*E12+E13*E13))
        CEXPQ2=CG(4)*DEXP(CG(5)*(E22+E33)**2+CG(6)*(E22*E33-E23*E23))
      ENDIF
      DW(1)=CEXPQ1*2.0d0*CG(2)*E11
      DW(2)=CEXPQ2*(2.0d0*CG(5)*(E22+E33)+CG(6)*E33)
      DW(3)=CEXPQ2*(2.0d0*CG(5)*(E22+E33)+CG(6)*E22)
      DW(4)= CEXPQ1*CG(3)*2.0d0*E12
      DW(5)= CEXPQ1*CG(3)*2.0d0*E13
      DW(6)=-CEXPQ2*CG(6)*2.0d0*E23

      CALL EXITS('USER53')
      RETURN
 9999 CALL ERRORS('USER53',ERROR)
      CALL EXITS('USER53')
      RETURN 1
      END


