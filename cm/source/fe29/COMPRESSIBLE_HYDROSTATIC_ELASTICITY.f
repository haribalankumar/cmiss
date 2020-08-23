      SUBROUTINE COMPRESSIBLE_HYDROSTATIC_ELASTICITY(dP_dV_V0,CG,TC,ZD,
     &  ERROR,*)
      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      REAL*8 CG(NMM),dP_dV_V0,TC(3,3),ZD
      LOGICAL UNIFORM
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 exp_term,lambda,lambda15,min_E

      CALL ENTERS('COMPRESSIBLE_HYDROSTATIC_ELASTICITY',*9999)

c      ZD=MAX(1.d0,ZD) !so that E cannot increase in lowest range
      lambda=ZD**(1.d0/3.d0)
      exp_term=DEXP(0.75d0*(3.d0*CG(2)+CG(3))*(lambda**2.d0-1.d0)**2.d0)

C MHT 10.2.11 Identified that the elastance calculation will give
C  increasing E values for 'ratio'==ZD <1.5. Have put in a linear
C  approximation to E in the region where ZD<1.5.
      IF(ZD.GE.1.5d0)THEN
        dP_dV_V0 = CG(1)/6.d0*(3.d0*CG(2)+CG(3))*exp_term*
     &    (3.d0*(3.d0*CG(2)+CG(3))*(lambda**2-1.d0)**2/lambda**2
     &    +(lambda**2+1.d0)/lambda**4)
      ELSE
        lambda15 = 1.5d0**(1.d0/3.d0)
        min_E = CG(1)/6.d0*(3.d0*CG(2)+CG(3))*exp_term*
     &    (3.d0*(3.d0*CG(2)+CG(3))*(lambda15**2-1.d0)**2/lambda15**2
     &    +(lambda15**2+1.d0)/lambda15**4)
        dP_dV_V0 = CG(1)*0.15d0+2.d0*(ZD-1.d0)*(min_E-CG(1)*0.15d0)
      ENDIF
      
      CALL EXITS('COMPRESSIBLE_HYDROSTATIC_ELASTICITY')
      RETURN
 9999 CALL ERRORS('COMPRESSIBLE_HYDROSTATIC_ELASTICITY',ERROR)
      CALL EXITS('COMPRESSIBLE_HYDROSTATIC_ELASTICITY')
      RETURN 1
      END
