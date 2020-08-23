      SUBROUTINE USER51(CG,DW,RK1,ERROR,*)

C#### Subroutine: USER51
C###  Description:
C###    USER51 returns the derivatives of the user-defined strain energy
C###    function of principal strain invariants I1,I2,I3,K1,K2 to
C###    subroutine ENERGY at current Gauss point.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      REAL*8 CG(NMM),DW(6),RK1
      CHARACTER ERROR*(*)

      CALL ENTERS('USER51',*9999)
C     DW(1)=CG(1)*DEXP(CG(2)*(RI1-3.0d0))
C     DW(2)=0.0d0
C     DW(3)=0.0d0
C     DW(4)=CG(3)*(DEXP(CG(4)*RK1)-1.0d0)
C     DW(5)=CG(5)*(DEXP(CG(6)*RK2))
      DW(1)=CG(1)
      DW(2)=CG(2)
      DW(3)=0.0d0
      DW(4)=2.0d0*CG(3)*RK1
      DW(5)=CG(4)
!news 28-MAY-1991 Function proposed by Humphrey/Yin.  JSW.
!     RK12=DSQRT(2.0d0*RK1+1.0d0)
!     DW(1)=CG(3)+CG(4)*(RK12-1.0d0)+2.0d0*CG(5)*(RI1-3.0d0)
!     DW(2)=0.0d0
!     DW(3)=0.0d0
!     DW(4)=2.0d0*CG(1)*(RK12-1.0d0)+3.0d0*CG(2)*(RK12-1.0d0)**2.0d0
!    '  +CG(4)*(RI1-3.0d0)
!     DW(4)=DW(4)/RK12
!     DW(5)=0.0d0
!newe

      CALL EXITS('USER51')
      RETURN
 9999 CALL ERRORS('USER51',ERROR)
      CALL EXITS('USER51')
      RETURN 1
      END


