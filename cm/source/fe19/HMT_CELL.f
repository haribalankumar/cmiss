      SUBROUTINE HMT_CELL(TIME,Y,DY,CONTROL,MODEL,SIZES,VARIANT,DERIVED,
     '  PARAM,PROTOCOL,AII,AIO,ARI,ARO,ERR_CODE)

C#### Subroutine: HMT_CELL
C###  Description:
C###    Solve isometric twitch using HMT cardiac mechanics model.

      IMPLICIT NONE

      INCLUDE 'cell_hmt.inc'
      INCLUDE 'cell_reserved.inc'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER SIZES(11)
      INTEGER CONTROL(*),MODEL(*),VARIANT,AII(*),AIO(*),ERR_CODE
      REAL*8 TIME(*),Y(*),PARAM(*),PROTOCOL(*),ARI(*),DY(*),
     '  DERIVED(*),ARO(*)
!     Local variables
      REAL*8 Cb_norm,C50,dPHI1,dPHI2,dPHI3,LENGTH_DIFF,Q,
     '  TTo,Tm_n,Tm_p50
!MHT 5/7/5 variables declared and not set, or set but not used
!     REAL*8 pCa      
!     LOGICAL DEBUG

      ERR_CODE=0

!MHT 5/7/5 variable DEBUG set but not used
C      DEBUG=.FALSE.

C *** DPN 05 July 2000 - need to check that the extension ratio has not
C *** exceeded the physiological limit, and that the proportion of actin
C *** sites does not fall outside the range 0 <= z < 1
C *** DPN 27 September 2000 - Still solve the model even if the values
C     are out of range, simply set them to the limits - this is to allow
C     distributed models to exceed the limits while trying to converge.
      IF ((Y(ExtensionRatio)-PARAM(MaxExtRatio)).GT.ZERO_TOL) THEN
        WRITE(*,
     '    '(''WARNING: Extension ratio above upper limit '',E12.5)')
     '    Y(ExtensionRatio)
        Y(ExtensionRatio) = PARAM(MaxExtRatio)
      ELSEIF ((PARAM(MinExtRatio)-Y(ExtensionRatio)).GT.ZERO_TOL) THEN
        WRITE(*,
     '    '(''WARNING: Extension ratio below lower limit '',E12.5)')
     '    Y(ExtensionRatio)
        Y(ExtensionRatio) = PARAM(MinExtRatio)
      ENDIF
      IF (Y(z).LT.0.0d0) THEN
        WRITE(*,'(''WARNING: Z below lower limit '',E12.5)') Y(z)
        Y(z) = 0.0d0
      ELSEIF ((Y(z)-1.0d0).GT.ZERO_TOL) THEN
        WRITE(*,'(''WARNING: Z above upper limit '',E12.5)') Y(z)
        Y(z) = 1.0d0
      ENDIF

C *** DPN 16 March 2001 - If the cell is lengthening, then solve the
C     steady state isometric tension equations, otherwise, when the cell
C     is contracting solve the full dynamic active contraction
C     equations.

      LENGTH_DIFF = Y(ExtensionRatio)-Y(ExtensionRatio_prev)
C *** DPN 05 July 2000 - Allow user to solve steady state equations
CXXXXXXXXXXXXXXX
C      IF(CONTROL(STEADY_STATE).EQ.1.OR.LENGTH_DIFF.GT.ZERO_TOL) THEN
      IF(CONTROL(STEADY_STATE).EQ.1) THEN
        !Solve steady state equations
        IF(CONTROL(ODE).EQ.1) THEN !evalulating DY for ODEs
C ***     Vm is always used!!
          DY(Vm)=0.0d0
          ! These need to be set so the integrator doesn't mess up the
          ! values of the variables
          DY(Cab)=0.0d0
          DY(z) = 0.0d0
        ELSE !Non-ODE variables
C ***     Intracellular Ca++
C         Only calculate Cai if the user wants to, i.e. they might be
C         setting the value via a time variable.
          IF(MODEL(CALCULATE_Cai).NE.0)
     '      Y(Cai)=PARAM(Ca_min)+PARAM(Ca_max)*
     '      TIME(TCell)/PARAM(Ca_tau)*DEXP(1.d0-TIME(TCell)/
     '      PARAM(Ca_tau))
          Y(Cab) = Y(Cai)*PARAM(Cab_max)/(Y(Cai)+
     '      (PARAM(Rho1)/PARAM(Rho0))*(1.0d0-1.0d0/PARAM(gamma)))
          Tm_n = PARAM(Tm_n_0)*(1.d0+PARAM(beta1)*
     '      (Y(ExtensionRatio)-1.d0))
          Tm_p50 = PARAM(Tm_p50_0)*(1.d0+PARAM(beta2)*
     '      (Y(ExtensionRatio)-1.d0))
          C50 = 10**(3.d0-Tm_p50) !mM
C *** DPN 24 March 2001 - scale z_SS by normalised [Ca]b
          Cb_norm = Y(Cab)/PARAM(Cab_max)
          Y(z) = Y(Cab)**Tm_n/(Y(Cab)**Tm_n+(C50/Cb_norm)**Tm_n)
          !Y(z) = Y(Cab)**Tm_n/(Y(Cab)**Tm_n+C50**Tm_n)
          !Y(z) = Y(Cai)**Tm_n/(Y(Cai)**Tm_n+C50**Tm_n)
          Y(IsometricTension) = PARAM(Tref)*(1.0d0+PARAM(beta0)*
     '      (Y(ExtensionRatio)-1.0d0))*Y(z)
          Y(Tension) = Y(IsometricTension)
          Y(ExtensionRatio_prev) = Y(ExtensionRatio)
          Y(phi1) = 0.0d0
          Y(phi2) = 0.0d0
          Y(phi3) = 0.0d0
        ENDIF !ODE
      ELSE ! solve dynamic equations

        IF(CONTROL(ODE).EQ.1) THEN !evalulating DY for ODEs

C ***     Vm is always used!!
          DY(Vm)=0.0d0

C ***     Troponin Kinetics
          IF (DABS(Y(IsometricTension)).LT.ZERO_TOL) THEN
            TTo = 0.0d0
          ELSE
            TTo = Y(Tension)/(PARAM(gamma)*Y(IsometricTension))
          ENDIF
          DY(Cab)=PARAM(Rho0)*Y(Cai)*(PARAM(Cab_max)-Y(Cab))-
     '      PARAM(Rho1)*(1.0d0-TTo)*Y(Cab)

C ***     Tropomyosin Kinetics
C         length dependence for n
          Tm_n   = PARAM(Tm_n_0)*(1.d0+PARAM(beta1)*
     '      (Y(ExtensionRatio)-1.d0))
C         length dependence for p50
          Tm_p50 = PARAM(Tm_p50_0)*(1.d0+PARAM(beta2)*
     '      (Y(ExtensionRatio)-1.d0))
          C50 = 10**(3.d0-Tm_p50) !mM

C DPN 05/08/98 - scale dz/dt by normalised [Ca]b
          Cb_norm = Y(Cab)/PARAM(Cab_max)

          DY(z)=PARAM(alfa0)*(((Y(Cab)/C50)*Cb_norm)**Tm_n*
     '      (1.d0-Y(z)) - Y(z))

        ELSE !evaluating new values for non-ODE state variables

C ***     Intracellular Ca++
C         Only calculate Cai if the user wants to, i.e. they might be
C         setting the value via a time variable.
          IF(MODEL(CALCULATE_Cai).NE.0)
     '      Y(Cai)=PARAM(Ca_min)+PARAM(Ca_max)*
     '      TIME(TCell)/PARAM(Ca_tau)*DEXP(1.d0-TIME(TCell)/
     '      PARAM(Ca_tau))

C ***     Tension-length-pCa - Isometric Tension
          Y(IsometricTension) = PARAM(Tref)*(1.d0+PARAM(beta0)*
     '      (Y(ExtensionRatio)-1.d0))*Y(z)

C ***     X-bridge kinetics - Active tension
C         Fading memory model...

C         1st Fading Memory term
          dPHI1=0.5d0*LENGTH_DIFF*(1.d0+DEXP(-PARAM(alpha1)*
     '      TIME(DTCell)))
          Y(phi1)=DEXP(-PARAM(alpha1)*TIME(DTCell))*Y(phi1)+dPHI1

C         2nd Fading Memory term
          dPHI2=0.5d0*LENGTH_DIFF*(1.d0+DEXP(-PARAM(alpha2)*
     '      TIME(DTCell)))
          Y(phi2)=DEXP(-PARAM(alpha2)*TIME(DTCell))*Y(phi2)+dPHI2

C         3rd Fading Memory term
          dPHI3=0.5d0*LENGTH_DIFF*(1.d0+DEXP(-PARAM(alpha3)*
     '      TIME(DTCell)))
          Y(phi3)=DEXP(-PARAM(alpha3)*TIME(DTCell))*Y(phi3)+dPHI3

          Q = PARAM(A1)*Y(phi1)+PARAM(A2)*Y(phi2)+PARAM(A3)*Y(phi3)

C *** DPN 16 March 2001 - This should now be handled by switching
C         between steady state and dynamic equations

C *** DPN 20 July 2000 - Need to fix HMT for lengthening experiments
c          IF(Q.LT.0.00d0) THEN
c            Y(Tension)=Y(IsometricTension)*(1.d0+PARAM(a)*Q)/(1.d0-Q)
c          ELSE
c            Y(Tension) = Y(IsometricTension)
c          ENDIF
CXXXXXXXXXXXXXXXXXXXXXXXXX

          IF(Q.GT.0.10d0) Q=0.10d0

          Y(Tension)=Y(IsometricTension)*(1.d0+PARAM(a)*Q)/(1.d0-Q)

          IF(Y(Tension).LT.ZERO_TOL) Y(Tension) = 0.0d0

C         Update previous extension ratio
          Y(ExtensionRatio_prev) = Y(ExtensionRatio)

        ENDIF
      ENDIF !steady state/dynamic

C pCa is never set...
C      IF(DEBUG) WRITE(*,'('' t='',F5.3,'' Cai='',E12.3,''(uM)'''
C     '  //''' pCa='',E12.3,'
C     '  //' '' Cb='',E12.3,''(uM)'''
C     '  //'  '' z='',E12.4)')
C     '  TIME(TCell),Y(Cai),pCa,Y(Cab),Y(z)

      RETURN
      END

