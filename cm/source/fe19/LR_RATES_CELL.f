      SUBROUTINE LR_RATES_CELL(TIME,Y,DY,CONTROL,MODEL,SIZES,VARIANT,
     '  CURRENTS,PARAM,PROTOCOL,AII,AIO,ARI,RATES,ERR_CODE)

C#### Subroutine: LR_RATES_CELL
C###  Description:
C###    Computes the rate coefficients for the gating variables.

      IMPLICIT NONE

      INCLUDE 'cell_lr.inc'
      INCLUDE 'cell_reserved.inc'

!     Parameter List
      INTEGER SIZES(11)
      INTEGER CONTROL(*),MODEL(*),VARIANT,AII(*),AIO(*),ERR_CODE
      REAL*8 TIME(*),Y(*),PARAM(*),PROTOCOL(*),ARI(*),DY(*),
     '  CURRENTS(27),RATES(21)
!     Local Variables
      REAL*8 d0,f0,Tau_d,Tau_f
      REAL*8 EK1

!     Compute the K1 reversal potential (mV)
      EK1 = PARAM(RTONF)*DLOG(PARAM(Ko)/Y(Ki))

! Fast Na channel gates m,h,j
      RATES(alpha_m) = 0.32d0*(Y(Vm)+47.13d0)/
     '  (1.d0-dexp(-0.1d0*(Y(Vm)+47.13d0)))
      RATES(beta_m)  = 0.08d0*dexp(-Y(Vm)/11.d0)
      IF(Y(Vm).GE.-40.d0) THEN
        RATES(alpha_h) = 0.d0
        RATES(beta_h)  = 1.d0/(0.13d0*(1.d0+dexp((Y(Vm)+10.66d0)/
     '    (-11.1d0))))
        RATES(alpha_j) = 0.d0
        RATES(beta_j)  = 0.3d0*dexp(-2.535d-7*Y(Vm))/
     '    (1.d0+dexp(-0.1d0*(Y(Vm)+32.d0)))
      ELSE IF(Y(Vm).LT.-40.d0) THEN
        RATES(alpha_h) = 0.135d0*dexp((80.d0+Y(Vm))/(-6.8d0))
        RATES(beta_h)  = 3.56d0*dexp(0.079d0*Y(Vm))+
     '    (3.1d5*dexp(0.35d0*Y(Vm)))
        RATES(alpha_j) = (-1.2714d5*dexp(0.2444d0*Y(Vm))-3.474d-5*
     '    dexp(-0.04391d0*Y(Vm)))*(Y(Vm)+37.78d0)/(1.d0+dexp(0.311d0*
     '    (Y(Vm)+79.23d0)))
        RATES(beta_j)  = 0.1212d0*dexp(-0.01052d0*Y(Vm))/
     '    (1.d0+dexp(-0.1378d0*(Y(Vm)+40.14d0)))
      ENDIF

! L-type Ca channel gates d,f
      d0    = 1.d0/(1.d0+dexp(-(Y(Vm)+10.d0)/6.24d0))
      Tau_d = d0*(1.d0-dexp(-(Y(Vm)+10.d0)/6.24d0))/
     '  (0.035d0*(Y(Vm)+10.d0))
      f0    = (1.d0/(1.d0+dexp((Y(Vm)+35.06d0)/8.6d0)))+(0.6d0/
     '  (1.d0+dexp((50.d0-Y(Vm))/20.d0)))
      Tau_f = 1.d0/(0.0197d0*dexp(-((0.0337d0*(Y(Vm)+10.d0))**2))+
     '  0.02d0)

      RATES(alpha_d) = d0/Tau_d
      RATES(beta_d)  = (1.d0-d0)/Tau_d
      RATES(alpha_f) = f0/Tau_f
      RATES(beta_f)  = (1.d0-f0)/Tau_f

! Time-dep. K channel gate x
      RATES(alpha_x) = 7.19d-5*(Y(Vm)+30.d0)/
     '  (1.d0-dexp(-0.148d0*(Y(Vm)+30.d0)))
      RATES(beta_x)  = 1.31d-4*(Y(Vm)+30.d0)/
     '  (-1.d0+dexp(0.0687d0*(Y(Vm)+30.d0)))

! Time-indep. K channel gate K1
      RATES(alpha_K1) = 1.02d0/(1.d0+dexp(0.2385d0*(Y(Vm)-
     '  EK1-59.215d0)))
      RATES(beta_K1)  = (0.49124d0*dexp(0.08032d0*(Y(Vm)-EK1+5.476d0))
     '  +dexp(0.06175d0*(Y(Vm)-EK1-594.31d0)))/(1.d0+dexp(-0.5143d0*
     '  (Y(Vm)-EK1+4.753d0)))

C *** No error
      ERR_CODE = 0

      RETURN
      END

