      SUBROUTINE HH_CELL(TIME,Y,DY,CONTROL,MODEL,SIZES,VARIANT,DERIVED,
     '  PARAM,PROTOCOL,AII,AIO,ARI,ARO,ERR_CODE)

C#### Subroutine: HH_CELL
C###  Description:
C###    Computes RHS of Hodgkin-Huxley equations.
C###    Voltage is in mV & time in ms.

      IMPLICIT NONE

      INCLUDE 'cell_hh.inc'
      INCLUDE 'cell_reserved.inc'

!     Parameters variables
      INTEGER SIZES(11)
      INTEGER CONTROL(*),MODEL(*),VARIANT,AII(*),AIO(*),ERR_CODE
      REAL*8 TIME(*),Y(*),PARAM(*),PROTOCOL(*),ARI(*),DY(*),DERIVED(*),
     '  ARO(*)

      !Local variables
      REAL*8 V,IStim,INa,IK,Il,a_n,a_m,a_h,b_n,b_m,b_h,ENa,EK,El

      !Equilibrium potentials (mV)
      ENa = PARAM(Er)+115.0d0
      EK = PARAM(Er)-12.0d0
      El = PARAM(Er)+10.613d0

      IStim=0.0d0
      IStim=PROTOCOL(PseudoIs)
      IF(TIME(TCell).GE.PROTOCOL(Is1start).AND.
     '  TIME(TCell).LT.PROTOCOL(Is1stop))
     '  IStim=IStim+PROTOCOL(Is1current)
      IF(TIME(TCell).GE.PROTOCOL(Is2start).AND.
     '  TIME(TCell).LT.PROTOCOL(Is2stop))
     '  IStim=IStim+PROTOCOL(Is2current)
      IF(PROTOCOL(IsFreqPeriod).GT.1.0d-6) THEN
        IF(DMOD(TIME(TCell),PROTOCOL(IsFreqPeriod)).LT.
     '    PROTOCOL(IsFreqDuration))
     '    IStim=IStim+PROTOCOL(IsFreqMag)
      ENDIF
      !convert from volume to area current (uA.mm^-2)
      Istim=Istim/PARAM(Am)

      !Compute rate constants
      V=Y(Vm)-PARAM(Er)
      a_h= 0.07d0*DEXP(-V/20.0d0)
      a_m=0.10d0*(25.0d0-V)/(DEXP(0.1d0*(25.0d0-V))-1.0d0)
      a_n=0.01d0*(10.0d0-V)/(DEXP(0.1d0*(10.0d0-V))-1.0d0)

      b_h = 1.0d0/(DEXP(0.1d0*(30.0d0-V)) + 1.0d0)
      b_m = 4.0d0*DEXP(-V/18.0d0)
      b_n = 0.125d0*DEXP(-V/80.0d0)

      !Compute o.d.e. RHSs
      INa=PARAM(gNa_max)*Y(m)**3*Y(h)*(Y(Vm)-ENa)
      IK=PARAM(gK_max)*Y(n)**4*(Y(Vm)-EK)
      Il=PARAM(gl_max)*(Y(Vm)-El)

      DY(Vm)=(IStim-INa-IK-Il)/PARAM(Cm)
      DY(h) = a_h*(1.d0-Y(h)) - b_h*Y(h)
      DY(m) = a_m*(1.d0-Y(m)) - b_m*Y(m)
      DY(n) = a_n*(1.d0-Y(n)) - b_n*Y(n)

      IF (CONTROL(RETURN_CURRENTS).NE.0) THEN
        DERIVED(DIStim)=IStim
        DERIVED(DINa)=INa
        DERIVED(DIK)=IK
        DERIVED(DIl)=Il
      ENDIF

C *** No error
      ERR_CODE=0

      RETURN
      END

C      SUBROUTINE HH_INIT_CELL(ICQS,LD,NDDATA,CELL_INIT,
C     '  PARAMETERS,RCQS,WD,XID,YQS,ZD,ERROR,*)
C
CC#### Subroutine: HH_INIT_CELL
CC###  Description:
CC###    Initialise the Hodgkin-Huxley model
CC***  Created by David Nickerson, May 1999
C
C      IMPLICIT NONE
C
C      INCLUDE 'cmiss$reference:cell00.cmn'
C      INCLUDE 'cmiss$reference:cell02.cmn'
C      INCLUDE 'cmiss$reference:cell_hh.inc'
C      INCLUDE 'cmiss$reference:cell_reserved.inc'
C      INCLUDE 'cmiss$reference:geom00.cmn'
C
C      !parameter list
C      INTEGER ICQS(NQIM),LD(NDM),NDDATA(0:NDM,0:NRM)
C      REAL*8 CELL_INIT(NUM_INITS),PARAM(NUM_PARAMS),RCQS(NQRM),
C     '  WD(NJM,NDM),XID(NIM,NDM),YQS(NIQSM,NQM),ZD(NJM,NDM)
C      CHARACTER ERROR*(*)
C      !local variables
C      INTEGER nd
C      REAL*8 F,R,T
C
CC *** ICQS
C      CELL_NUM_AII=1
C      CELL_NUM_AIO=1
C      CELL_NUM_CONTROL=1
C      CELL_NUM_MODEL=1
C      CELL_NUM_VARIANTS=1
C
C      CELL_AII_OFFSET=1
C      CELL_AIO_OFFSET=CELL_AII_OFFSET+CELL_NUM_AII
C      CELL_CONTROL_OFFSET=CELL_AIO_OFFSET+CELL_NUM_AIO
C      CELL_MODEL_OFFSET=CELL_CONTROL_OFFSET+CELL_NUM_CONTROL
C      CELL_VARIANT_OFFSET=CELL_MODEL_OFFSET+CELL_NUM_MODEL
C
C      NQIT=CELL_NUM_AII+CELL_NUM_AIO+CELL_NUM_CONTROL+CELL_NUM_MODEL+
C     '  CELL_NUM_VARIANTS
C      CALL ASSERT(NQIM.GE.NQIT,'>>Increase NQIM',ERROR,*9999)
C
CC *** RCQS
C      CELL_NUM_ARI=3
C      CELL_NUM_ARO=1
C      CELL_NUM_PARAMETERS=NUM_PARAMS
C      CELL_NUM_PROTOCOL=1
C
C      CELL_ARI_OFFSET=1
C      CELL_ARO_OFFSET=CELL_ARI_OFFSET+CELL_NUM_ARI
C      CELL_PARAMETERS_OFFSET=CELL_ARO_OFFSET+CELL_NUM_ARO
C      CELL_PROTOCOL_OFFSET=CELL_PARAMETERS_OFFSET+CELL_NUM_PARAMETERS
C
C      NQRT=CELL_NUM_ARI+CELL_NUM_ARO+CELL_NUM_PARAMETERS+
C     '  CELL_NUM_PROTOCOL
C      CALL ASSERT(NQRM.GE.NQRT,'>>Increase NQRM',ERROR,*9999)
C
CC *** YQS
C      CELL_NUM_DERIVED=10
C      CELL_NUM_STATE=NUM_INITS
C
C      CELL_DERIVED_OFFSET=1
C      CELL_STATE_OFFSET=CELL_DERIVED_OFFSET+CELL_NUM_DERIVED
C
C      NIQST=CELL_NUM_DERIVED+CELL_NUM_STATE
C      CALL ASSERT(NIQSM.GE.NIQST,'>>Increase NIQSM',ERROR,*9999)
C
CC *** Assign the parameters
C      DO nd=0,CELL_NUM_PARAMETERS-1
C        RCQS(CELL_PARAMETERS_OFFSET+nd)=PARAM(nd+1)
C      ENDDO
C
CC *** Set the number of "signals" required
C      IF(OUTPUT_DERIVED) THEN
C        NDDATA(0,0)=CELL_NUM_STATE+CELL_NUM_DERIVED
C        NDDATA(0,1)=NDDATA(0,0)
C      ELSE
C        NDDATA(0,0)=CELL_NUM_STATE
C        NDDATA(0,1)=CELL_NUM_STATE
C      ENDIF
CC *** Check that the maximum number of data points is greater than the
CC     number of "signals" required
C      CALL ASSERT(NDDATA(0,0).LE.NDM,'>>Increase NDM',ERROR,*9999)
C      DO nd=1,NDDATA(0,0)
C        WD(1,nd)=1.0d0
C        WD(2,nd)=1.0d0
C        WD(3,nd)=1.0d0
C        ZD(1,nd)=0.0d0
C        ZD(2,nd)=0.0d0
C        XID(1,nd)=0.0d0
C        XID(2,nd)=0.0d0
C        LD(nd)=1
C        NDDATA(nd,1)=nd
C      ENDDO
C      !initialise stuff
C      RCQS(CELL_VARIANT_OFFSET)=1!type of preparation(unused yet)?????
C      R=8.314d0
C      T=RCQS(CELL_PARAMETERS_OFFSET-1+Temp)+273.15d0
C      F=96487.d0
C      IF(DABS(RCQS(CELL_PARAMETERS_OFFSET-1+ic_Na)).GT.1.d-12) THEN
C        RCQS(CELL_ARI_OFFSET-1+VNa) = 1000.d0*R*T/F *
C     '    DLOG(RCQS(CELL_PARAMETERS_OFFSET-1+ec_Na)/
C     '    RCQS(CELL_PARAMETERS_OFFSET-1+ic_Na)) !in mV
C      ELSE
C        RCQS(CELL_ARI_OFFSET-1+VNa) = 0.d0
C      ENDIF
C      IF(DABS(RCQS(CELL_PARAMETERS_OFFSET-1+ic_K)).GT.1.d-12) THEN
C        RCQS(CELL_ARI_OFFSET-1+VK) = 1000.d0*R*T/F *
C     '    DLOG(RCQS(CELL_PARAMETERS_OFFSET-1+ec_K)/
C     '    RCQS(CELL_PARAMETERS_OFFSET-1+ic_K)) !in mV
C      ELSE
C        RCQS(CELL_ARI_OFFSET-1+VK) = 0.d0
C      ENDIF
C !      Sea Hare (Aplysia)
CC      Vrest = 1000.d0*R*T/F*DLOG((ec_K+0.12d0*ec_Na+1.44d0*51.d0)/
CC     '  (ic_K+0.12d0*ic_Na+1.44d0*485.d0))
C !      Squid
C      RCQS(CELL_ARI_OFFSET-1+Vrest) = 1000.d0*R*T/F*
C     '  DLOG((RCQS(CELL_PARAMETERS_OFFSET-1+ec_K)+
C     '  0.04d0*RCQS(CELL_PARAMETERS_OFFSET-1+ec_Na)+
C     '  0.45d0*52.d0)/(RCQS(CELL_PARAMETERS_OFFSET-1+ic_K)+0.04d0*
C     '  RCQS(CELL_PARAMETERS_OFFSET-1+ic_Na)+0.45d0*560.d0))
C      YQS(CELL_DERIVED_OFFSET-1+Stimulus,1) = 0.0d0
C      YQS(CELL_DERIVED_OFFSET-1+INa,1)      = 0.0d0
C      YQS(CELL_DERIVED_OFFSET-1+IK,1)       = 0.0d0
C      YQS(CELL_DERIVED_OFFSET-1+Ileak,1)    = 0.0d0
C      !Use computed resting potential as initial
C      YQS(CELL_STATE_OFFSET-1+Vm,1) = RCQS(CELL_ARI_OFFSET-1+Vrest)
C      YQS(CELL_DERIVED_OFFSET-1+alpha_h,1)=0.07d0
C      YQS(CELL_DERIVED_OFFSET-1+beta_h,1)=1.d0/(DEXP(3.d0)+1.d0)
C      YQS(CELL_DERIVED_OFFSET-1+alpha_m,1)=2.5d0/(DEXP(2.5d0)-1.d0)
C      YQS(CELL_DERIVED_OFFSET-1+beta_m,1)=4.d0
C      YQS(CELL_DERIVED_OFFSET-1+alpha_n,1)=0.1d0/(DEXP(1.d0)-1)
C      YQS(CELL_DERIVED_OFFSET-1+beta_n,1)=0.125d0
C      YQS(CELL_STATE_OFFSET-1+h,1)=YQS(CELL_DERIVED_OFFSET-1+alpha_h,1)/
C     '  (YQS(CELL_DERIVED_OFFSET-1+alpha_h,1)+
C     '  YQS(CELL_DERIVED_OFFSET-1+beta_h,1))
C      YQS(CELL_STATE_OFFSET-1+m,1)=YQS(CELL_DERIVED_OFFSET-1+alpha_m,1)/
C     '  (YQS(CELL_DERIVED_OFFSET-1+alpha_m,1)+
C     '  YQS(CELL_DERIVED_OFFSET-1+beta_m,1))
C      YQS(CELL_STATE_OFFSET-1+n,1)=YQS(CELL_DERIVED_OFFSET-1+alpha_n,1)/
C     '  (YQS(CELL_DERIVED_OFFSET-1+alpha_n,1)+
C     '  YQS(CELL_DERIVED_OFFSET-1+beta_n,1))
CC *** Put initial values into ZD to be written to the signal file
C      DO nd=1,CELL_NUM_STATE
C        ZD(3,nd)=YQS(CELL_STATE_OFFSET-1+nd,1)
C      ENDDO
CC *** If required, put the derived values into ZD
C      IF (OUTPUT_DERIVED) THEN
C        DO nd=CELL_NUM_STATE+1,NDDATA(0,0)
C          ZD(3,nd)=YQS(CELL_DERIVED_OFFSET-1+nd-CELL_NUM_STATE,1)
C        ENDDO
C      ENDIF
C
C      CALL EXITS('HH_INIT_CELL')
C      RETURN
C 9999 CALL ERRORS('HH_INIT_CELL',ERROR)
C      CALL EXITS('HH_INIT_CELL')
C      RETURN 1
C      END


