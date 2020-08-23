
      SUBROUTINE COMPRESSIBLE_HYDROSTATIC_PRESSURE(pressure,CG,TC,ZD,
     &  UNIFORM,ERROR,*)
      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      REAL*8 CG(NMM),pressure,TC(3,3),ZD
      LOGICAL UNIFORM
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 exp_term,lambda

      CALL ENTERS('COMPRESSIBLE_HYDROSTATIC_PRESSURE',*9999)

      IF(UNIFORM)THEN
        ZD=MAX(1.d0,ZD)
        lambda=ZD**(1.d0/3.d0)
        exp_term=0.75d0*(3.d0*CG(2)+CG(3))*(lambda**2.d0-1.d0)**2.d0
        pressure=CG(1)/2.d0*(3.d0*CG(2)+CG(3))*(lambda**2.d0-1.d0)
     &    *DEXP(exp_term)/lambda
      ELSE
        pressure=(TC(1,1)+TC(2,2)+TC(3,3))/3.d0 !PRESSURE
      ENDIF
      
      CALL EXITS('COMPRESSIBLE_HYDROSTATIC_PRESSURE')
      RETURN
 9999 CALL ERRORS('COMPRESSIBLE_HYDROSTATIC_PRESSURE',ERROR)
      CALL EXITS('COMPRESSIBLE_HYDROSTATIC_PRESSURE')
      RETURN 1
      END
