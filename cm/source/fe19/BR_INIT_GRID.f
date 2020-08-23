      SUBROUTINE BR_INIT_GRID(NQLIST,CQ,RCQS,YQ,YQS,ERROR,*)

C#### Subroutine: BR_INIT_GRID
C###  Description:
C###    BR_INIT_GRID initialises arrays for the BR
C###    ionic current model.
C**** Written by Martin Buist, 17 May 1999

      IMPLICIT NONE

      INCLUDE 'cell02.cmn'
      INCLUDE 'cellml.cmn'
      INCLUDE 'cell_br.inc'
      INCLUDE 'cell_reserved.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'nqloc00.inc'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER NQLIST(0:NQM)
      REAL*8 CQ(NMM,NQM),RCQS(NQRM),YQ(NYQM,NIQM,NAM),
     '  YQS(NIQSM,NQM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nq
      REAL*8 Ad,Af,Ah,Aj,Alm,Ax1,Bd,Bf,Bh,Bj,Bm,Bx1,V

      CALL ENTERS('BR_INIT_GRID',*9999)

C *** DPN 18 June 1999
C *** Simplifed models not going through ipcell can't have spatially
C     varying material parameters, unless they are set-up here.
      CELL_SPATIALLY_VARYING=.FALSE.

!ICQS
      CELL_NUM_AII(1)=1
      CELL_NUM_AIO(1)=1
      CELL_NUM_CONTROL(1)=1
      CELL_NUM_MODEL(1)=1
      CELL_NUM_VARIANTS=1

      CELL_AII_OFFSET(1)=CELL_VARIANT_OFFSET+1
      CELL_AIO_OFFSET(1)=CELL_AII_OFFSET(1)+CELL_NUM_AII(1)
      CELL_CONTROL_OFFSET(1)=CELL_AIO_OFFSET(1)+CELL_NUM_AIO(1)
      CELL_MODEL_OFFSET(1)=CELL_CONTROL_OFFSET(1)+CELL_NUM_CONTROL(1)
      ! This is now a parameter
      !CELL_VARIANT_OFFSET=CELL_MODEL_OFFSET+CELL_NUM_MODEL

      NQIT=CELL_MODEL_OFFSET(1)+CELL_NUM_MODEL(1)-1
      CALL ASSERT(NQIM.GE.NQIT,' >>Increase NQIM (>=NQIT)',ERROR,*9999)

!RCQS
      CELL_NUM_ARI(1)=1
      CELL_NUM_ARO(1)=1
      CELL_NUM_PARAMETERS(1)=7
      CELL_NUM_PROTOCOL(1)=4

      CELL_ARI_OFFSET(1)=1
      CELL_ARO_OFFSET(1)=CELL_ARI_OFFSET(1)+CELL_NUM_ARI(1)
      CELL_PARAMETERS_OFFSET(1)=CELL_ARO_OFFSET(1)+CELL_NUM_ARO(1)
      CELL_PROTOCOL_OFFSET(1)=CELL_PARAMETERS_OFFSET(1)
     '  +CELL_NUM_PARAMETERS(1)

      NQRT=CELL_PROTOCOL_OFFSET(1)+CELL_NUM_PROTOCOL(1)-1
      CALL ASSERT(NQRM.GE.NQRT,' >>Increase NQRM (>=NQRT)',ERROR,*9999)

!YQS
      CELL_NUM_DERIVED(1)=1
      CELL_NUM_STATE(1)=9
      CELL_NUM_ODE(1)=CELL_NUM_STATE(1)

      CELL_STATE_OFFSET(1)=1
      CELL_DERIVED_OFFSET(1)=CELL_STATE_OFFSET(1)+CELL_NUM_STATE(1)

      NIQST=CELL_DERIVED_OFFSET(1)+CELL_NUM_DERIVED(1)-1
      CALL ASSERT(NIQSM.GE.NIQST,' >>Increase NIQSM (>=NIQST)',
     '  ERROR,*9999)

!YQ
      CALL ASSERT(NIQM.GE.3,'>>Increase NIQM (>=3)',ERROR,*9999)
      CALL ASSERT(NQM.GE.3,'>>Increase NQM (>=3)',ERROR,*9999)
C KAT 2001-04-12: unused
C      NQLIST(0)=3

      CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,NQLIST(1),NIQ_V,
     '  ERROR,*9999)
      IF(NQLIST(1).EQ.0) THEN
        CALL NIQ_LOC(NIQ_ALLOCATE_AND_LOCK,NIQ_ION,NQLIST(1),
     '    NIQ_V,ERROR,*9999)
      ENDIF

      CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,NQLIST(2),NIQ_BNDRY,
     '  ERROR,*9999)
      IF(NQLIST(2).EQ.0) THEN
        CALL NIQ_LOC(NIQ_ALLOCATE_AND_LOCK,NIQ_ION,NQLIST(2),
     '    NIQ_BNDRY,ERROR,*9999)
      ENDIF

C KAT 2001-04-12: moved to after user specified initial conditions
C      CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,NQLIST(3),NIQ_OLDSOLN,
C     '  ERROR,*9999)
C      IF(NQLIST(3).EQ.0) THEN
C        CALL NIQ_LOC(NIQ_ALLOCATE_AND_LOCK,NIQ_ION,NQLIST(3),
C     '    NIQ_OLDSOLN,ERROR,*9999)
C      ENDIF

      DO nq=1,NQT
!STATE
        V=CQ(9,nq)
        IF(KTYP33.EQ.1) THEN !Beeler-Reuter sodium model
          Alm=-1.d0*(V+47.d0)/(DEXP(-0.1d0*(V+47.d0))-1.d0)
          Bm=40.d0*DEXP(-0.056d0*(V+72.d0))
          Ah=0.126d0*DEXP(-0.25d0*(V+77.d0))
          Bh=1.7d0/(DEXP(-.082d0*(V+22.5d0))+1.d0)
          Aj=0.055d0*DEXP(-0.25d0*(V+78.d0))/
     '      (DEXP(-0.2d0*(V+78.d0))+1.d0)
          Bj=0.3d0/(DEXP(-0.1d0*(V+32.d0))+1.d0)

        ELSE IF(KTYP33.EQ.2) THEN !Ebihara-Johnson sodium model
          Alm=(0.32d0*(V+47.13d0))/(1.d0-DEXP(-V-47.13d0))
          Bm=0.08d0*DEXP(-V/11.d0)
          IF(V.GE.-40.d0) THEN
            Ah=0.d0
            Bh=1.d0/(0.13d0*(DEXP((V+10.66d0)/(-11.1d0))+1.d0))
          ELSE
            Ah=0.135d0*DEXP((-80.d0-V)/6.8d0)
            Bh=3.56d0*DEXP(0.079d0*V)+3.1d5*DEXP(0.35d0*V)
          ENDIF
          Aj=0.d0
          Bj=0.d0

        ELSE IF(KTYP33.EQ.3) THEN !Drouhard-Roberge sodium model
          Alm=0.9d0*(V+42.65d0)/(1.d0-DEXP(-0.22d0*(V+42.65d0)))
          Bm=1.437d0*DEXP(-0.085d0*(V+39.75d0))
          Ah=0.1d0*DEXP(-0.193d0*(V+79.65d0))
          Bh=1.7d0/(1.d0+DEXP(-0.095d0*(V+20.5d0)))
          Aj=0.d0
          Bj=0.d0
        ELSE IF(KTYP33.EQ.4) THEN
          !Modified Drouhard-Roberge Na model
          !Alpha_m and Beta_m
          IF(V.LE.100.0d0) THEN
            Alm=0.9d0*((V+42.65d0)/(1.0d0-DEXP((-0.22d0*V)-
     '        9.383d0)))
          ELSE
            Alm=890.943789d0*((DEXP((0.0486479d0*V)-4.8647916d0))/
     '        (1.0d0+(5.93962526d0*DEXP((0.0486479d0*V)
     '        -4.8647916d0))))
          ENDIF
          IF(V.LE.-85.0d0) THEN
            Bm=100.0d0/(1.0d0+(0.4864082d0*DEXP((0.2597504d0*V)+
     '        22.0787804d0)))
          ELSE
            Bm=1.437d0*DEXP((-0.085d0*V)-3.37875d0)
          ENDIF

          !Alpha_h and Beta_h
          IF(V.LE.-90.0d0) THEN
            Ah=-12.0662845d0-(0.1422598d0*V)
          ELSE
            Ah=0.1d0*DEXP((-0.193d0*V)-15.37245d0)
          ENDIF
          Bh=1.7d0/(1.0d0+DEXP((-0.095d0*V)-1.9475d0))
          Aj = 0.d0
          Bj = 0.d0
        ENDIF !sodium model

        Ad=0.095d0*DEXP(-0.01d0*(V-5.d0))/
     '    (1.d0+DEXP(-0.072d0*(V-5.d0)))
        Bd=0.07d0*DEXP(-(V+44.d0)/59.d0)/
     '    (1.d0+DEXP(0.05d0*(V+44.d0)))
        Af=0.012d0*DEXP(-0.008d0*(V+28.d0))/
     '    (1.d0+DEXP(0.15d0*(V+28.d0)))
        Bf=0.0065d0*DEXP(-0.02d0*(V+30.d0))/
     '    (1.d0+DEXP(-0.2d0*(V+30.d0)))
        Ax1 = 5.d-4*DEXP(0.083d0*(V+50.0d0))/
     '    (1.d0+DEXP(0.057d0*(V+50.d0)))
        Bx1=0.0013d0*DEXP(-0.06d0*(V+20.d0))/
     '    (1.d0+DEXP(-0.04d0*(V+20.d0)))

        YQ(nq,NQLIST(1),1)=CQ(9,nq)
        YQ(nq,NQLIST(2),1)=0.0d0
C        YQ(nq,NQLIST(3),1)=CQ(9,nq)
        YQS(CELL_STATE_OFFSET(1)+Vm-1,nq)=CQ(9,nq)
        IF(DABS(Alm+Bm).GT.ZERO_TOL) THEN
          YQS(CELL_STATE_OFFSET(1)+m-1,nq)=Alm/(Alm+Bm) !m
        ELSE
          YQS(CELL_STATE_OFFSET(1)+m-1,nq)=0.0d0
        ENDIF
        IF(DABS(Ah+Bh).GT.ZERO_TOL) THEN
          YQS(CELL_STATE_OFFSET(1)+h-1,nq)=Ah/(Ah+Bh) !h
        ELSE
          YQS(CELL_STATE_OFFSET(1)+h-1,nq)=0.0d0
        ENDIF
        IF(DABS(Aj+Bj).GT.ZERO_TOL) THEN
          YQS(CELL_STATE_OFFSET(1)+j-1,nq)=Aj/(Aj+Bj) !j
        ELSE
          YQS(CELL_STATE_OFFSET(1)+j-1,nq)=0.0d0
        ENDIF
        IF(DABS(Ad+Bd).GT.ZERO_TOL) THEN
          YQS(CELL_STATE_OFFSET(1)+d-1,nq)=Ad/(Ad+Bd) !d
        ELSE
          YQS(CELL_STATE_OFFSET(1)+d-1,nq)=0.0d0
        ENDIF
        IF(DABS(Af+Bf).GT.ZERO_TOL) THEN
          YQS(CELL_STATE_OFFSET(1)+f1-1,nq)=Af/(Af+Bf) !f1
        ELSE
          YQS(CELL_STATE_OFFSET(1)+f1-1,nq)=0.0d0
        ENDIF
        IF(DABS(Ax1+Bx1).GT.ZERO_TOL) THEN
          YQS(CELL_STATE_OFFSET(1)+x1-1,nq)=Ax1/(Ax1+Bx1) !x1
        ELSE
          YQS(CELL_STATE_OFFSET(1)+x1-1,nq)=0.0d0
        ENDIF
        YQS(CELL_STATE_OFFSET(1)+Cai-1,nq)=1.0d-4 !Cai (in nmol/mm^3)
        YQS(CELL_STATE_OFFSET(1)+G-1,nq)=0.0d0 !G
!PARAMETERS
C *** DPN 14 March 2000 - fixing up units
c        RCQS(CELL_PARAMETERS_OFFSET+Cm-1)=CQ(1,nq)*1.0d-6
        RCQS(CELL_PARAMETERS_OFFSET(1)+Cm-1)=CQ(1,nq)
        RCQS(CELL_PARAMETERS_OFFSET(1)+Am-1)=CQ(2,nq)
        RCQS(CELL_PARAMETERS_OFFSET(1)+Vmrest-1)=CQ(9,nq)
        RCQS(CELL_PARAMETERS_OFFSET(1)+VNa-1)=CQ(10,nq)
        RCQS(CELL_PARAMETERS_OFFSET(1)+GNa-1)=CQ(11,nq)
        RCQS(CELL_PARAMETERS_OFFSET(1)+GNaC-1)=CQ(12,nq)
        RCQS(CELL_PARAMETERS_OFFSET(1)+Gs-1)=CQ(13,nq)
!PROTOCOL
        RCQS(CELL_PROTOCOL_OFFSET(1)+PseudoIs-1)=0.0d0
        RCQS(CELL_PROTOCOL_OFFSET(1)+Is1current-1)=0.0d0
      ENDDO

      CALL EXITS('BR_INIT_GRID')
      RETURN
 9999 CALL ERRORS('BR_INIT_GRID',ERROR)
      CALL EXITS('BR_INIT_GRID')
      RETURN 1
      END


C      SUBROUTINE LR_INIT_GRID(NQLIST,CQ,RCQS,YQ,YQS,ERROR,*)
C
CC#### Subroutine: LR_INIT_GRID
CC###  Description:
CC###    LR_INIT_GRID initialises arrays for the Luo-Rudy 2
CC###    ionic current model.
CC**** Written by Martin Buist, 17 May 1999
C
C      IMPLICIT NONE
C
C      INCLUDE 'cmiss$reference:call00.cmn'
C      INCLUDE 'cmiss$reference:cell02.cmn'
C      INCLUDE 'cmiss$reference:cell_lr.inc'
C      INCLUDE 'cmiss$reference:cell_reserved.inc'
C      INCLUDE 'cmiss$reference:geom00.cmn'
C      INCLUDE 'cmiss$reference:nqloc00.inc'
C
C!     Parameter List
C      INTEGER NQLIST(0:NQM)
C      REAL*8 CQ(NMM,NQM),RCQS(NQRM),YQ(NYQM,NIQM,NAM),
C     '  YQS(NIQSM,NQM)
C      CHARACTER ERROR*(*)
C!     Local Variables
C      INTEGER nq
C
C      CALL ENTERS('LR_INIT_GRID',*9999)
C
C      CALL ASSERT(CALL_CELL,' >>Define cell first',ERROR,*9999)
C
C!ICQS
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
C      NQIT=CELL_VARIANT_OFFSET+CELL_NUM_VARIANTS-1
C      CALL ASSERT(NQIM.GE.NQIT,' >>Increase NQIM (>=NQIT)',ERROR,*9999)
C
C
C!RCQS
C      CELL_NUM_ARI=1
C      CELL_NUM_ARO=23
C      CELL_NUM_PARAMETERS=60
C      CELL_NUM_PROTOCOL=4
C
C      CELL_ARI_OFFSET=1
C      CELL_ARO_OFFSET=CELL_ARI_OFFSET+CELL_NUM_ARI
C      CELL_PARAMETERS_OFFSET=CELL_ARO_OFFSET+CELL_NUM_ARO
C      CELL_PROTOCOL_OFFSET=CELL_PARAMETERS_OFFSET+CELL_NUM_PARAMETERS
C
C      NQRT=CELL_PROTOCOL_OFFSET+CELL_NUM_PROTOCOL-1
C      CALL ASSERT(NQRM.GE.NQRT,' >>Increase NQRM (>=NQRT)',ERROR,*9999)
C
C!YQS
C      CELL_NUM_DERIVED=27
C      CELL_NUM_STATE=14
C
C      CELL_DERIVED_OFFSET=1
C      CELL_STATE_OFFSET=CELL_DERIVED_OFFSET+CELL_NUM_DERIVED
C
C      NIQST=CELL_STATE_OFFSET+CELL_NUM_STATE-1
C      CALL ASSERT(NIQSM.GE.NQISM,' >>Increase NIQSM (>=NQISM)',
C     '  ERROR,*9999)
C
C!YQ
C      CALL ASSERT(NIQM.GE.3,'>>Increase NIQM (>=3)',ERROR,*9999)
C      CALL ASSERT(NQM.GE.3,'>>Increase NQM (>=3)',ERROR,*9999)
C      NQLIST(0)=3
C
C      CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,NQLIST(1),NIQ_V,
C     '  ERROR,*9999)
C      IF(NQLIST(1).EQ.0) THEN
C        CALL NIQ_LOC(NIQ_ALLOCATE_AND_LOCK,NIQ_ION,NQLIST(1),
C     '    NIQ_V,ERROR,*9999)
C      ENDIF
C
C      CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,NQLIST(2),NIQ_BNDRY,
C     '  ERROR,*9999)
C      IF(NQLIST(2).EQ.0) THEN
C        CALL NIQ_LOC(NIQ_ALLOCATE_AND_LOCK,NIQ_ION,NQLIST(2),
C     '    NIQ_BNDRY,ERROR,*9999)
C      ENDIF
C
C      CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,NQLIST(3),NIQ_OLDSOLN,
C     '  ERROR,*9999)
C      IF(NQLIST(3).EQ.0) THEN
C        CALL NIQ_LOC(NIQ_ALLOCATE_AND_LOCK,NIQ_ION,NQLIST(3),
C     '    NIQ_OLDSOLN,ERROR,*9999)
C      ENDIF
C
C      DO nq=1,NQT
C!STATE
C        YQ(nq,NQLIST(1),1)=CELL_INIT(1)
C        YQ(nq,NQLIST(2),1)=0.0d0
C        YQ(nq,NQLIST(3),1)=CELL_INIT(1)
C        YQS(CELL_STATE_OFFSET+Vm-1,nq)=CELL_INIT(1)
C        YQS(CELL_STATE_OFFSET+Nai-1)=CELL_INIT(2)
C        YQS(CELL_STATE_OFFSET+Cai-1)=CELL_INIT(3)
C        YQS(CELL_STATE_OFFSET+x-1)=CELL_INIT(4)
C        YQS(CELL_STATE_OFFSET+d-1)=CELL_INIT(5)
C        YQS(CELL_STATE_OFFSET+f-1)=CELL_INIT(6)
C        YQS(CELL_STATE_OFFSET+h-1)=CELL_INIT(7)
C        YQS(CELL_STATE_OFFSET+m-1)=CELL_INIT(8)
C        YQS(CELL_STATE_OFFSET+j-1)=CELL_INIT(9)
C        YQS(CELL_STATE_OFFSET+Ki-1)=CELL_INIT(10)
C        YQS(CELL_STATE_OFFSET+CaNSR-1)=CELL_INIT(11)
C        YQS(CELL_STATE_OFFSET+CaJSR-1)=CELL_INIT(12)
C        YQS(CELL_STATE_OFFSET+Cab-1)=CELL_INIT(13)
C        YQS(CELL_STATE_OFFSET+Cao-1)=CELL_INIT(14)
C!DERIVED
C        YQS(CELL_DERIVED_OFFSET+Stimulus-1)=0.0d0
C        YQS(CELL_DERIVED_OFFSET+INa-1)=0.0d0
C        YQS(CELL_DERIVED_OFFSET+ICaLCa-1)=0.0d0
C        YQS(CELL_DERIVED_OFFSET+ICaLK-1)=0.0d0
C        YQS(CELL_DERIVED_OFFSET+ICaLNa-1)=0.0d0
C        YQS(CELL_DERIVED_OFFSET+ICaL-1)=0.0d0
C        YQS(CELL_DERIVED_OFFSET+IK-1)=0.0d0
C        YQS(CELL_DERIVED_OFFSET+IK1-1)=0.0d0
C        YQS(CELL_DERIVED_OFFSET+IKp-1)=0.0d0
C        YQS(CELL_DERIVED_OFFSET+INaCa-1)=0.0d0
C        YQS(CELL_DERIVED_OFFSET+INaK-1)=0.0d0
C        YQS(CELL_DERIVED_OFFSET+InsK-1)=0.0d0
C        YQS(CELL_DERIVED_OFFSET+InsNa-1)=0.0d0
C        YQS(CELL_DERIVED_OFFSET+InsCa-1)=0.0d0
C        YQS(CELL_DERIVED_OFFSET+IpCa-1)=0.0d0
C        YQS(CELL_DERIVED_OFFSET+ICab-1)=0.0d0
C        YQS(CELL_DERIVED_OFFSET+INab-1)=0.0d0
C        YQS(CELL_DERIVED_OFFSET+Iv-1)=0.0d0
C        YQS(CELL_DERIVED_OFFSET+Irel-1)=0.0d0
C        YQS(CELL_DERIVED_OFFSET+Ileak-1)=0.0d0
C        YQS(CELL_DERIVED_OFFSET+Iup-1)=0.0d0
C        YQS(CELL_DERIVED_OFFSET+I_Na-1)=0.0d0
C        YQS(CELL_DERIVED_OFFSET+I_Ca-1)=0.0d0
C        YQS(CELL_DERIVED_OFFSET+I_K-1)=0.0d0
C        YQS(CELL_DERIVED_OFFSET+F_NSR-1)=0.0d0
C        YQS(CELL_DERIVED_OFFSET+F_JSR-1)=0.0d0
C        YQS(CELL_DERIVED_OFFSET+Itr-1)=0.0d0
C
C!PARAMETERS
C        RCQS(CELL_PARAMETERS_OFFSET+Cm-1)=
C        RCQS(CELL_PARAMETERS_OFFSET+Am-1)=
C!PROTOCOL
C        RCQS(CELL_PROTOCOL_OFFSET+PseudoIs-1)=0.0d0
C        RCQS(CELL_PROTOCOL_OFFSET+Is1current-1)=0.0d0
C
C      ENDDO
C
C      CALL EXITS('LR_INIT_GRID')
C      RETURN
C 9999 CALL ERRORS('LR_INIT_GRID',ERROR)
C      CALL EXITS('LR_INIT_GRID')
C      RETURN 1
C      END


C      REAL*8 FUNCTION CUBIC_ION(PHIM,IONIC_PARAMS)
C
CC#### Function: CUBIC_ION
CC###  Type: REAL*8
CC###  Description:
CC###    CUBIC_ION is the cubic ionic equation for bidomain solution.
C
C      IMPLICIT NONE
C      INCLUDE 'cmiss$reference:ktyp30.cmn'
C!     Parameter List
C      REAL*8 IONIC_PARAMS(12),PHIM
C!     Local Variables
C      REAL*8 MEMCOND,PHI,PLATEAU,REST,THRESHOLD,VTH,VPL
C
CC *** Ionic current cubic equation as given by Colli Franzone
CC *** Result in mA/mm^2
CC *** IONIC_PARAMS(1) is resting potential in mV
CC *** IONIC_PARAMS(2) is plateau potential in mV
CC *** IONIC_PARAMS(3) is threshold potential in mV
CC *** IONIC_PARAMS(4) is membrane conductance in uMho/mm^2
C
C      REST=IONIC_PARAMS(1)
C      PLATEAU=IONIC_PARAMS(2)-REST
C      THRESHOLD=IONIC_PARAMS(3)-REST
CC ** 10^-6 correction for the membrane conductance
C      MEMCOND=IONIC_PARAMS(4)*1.d-6
C
C      IF(PHIM.LT.REST.OR.PHIM.GT.PLATEAU) THEN
C        CUBIC_ION=0.0d0
C      ELSE
C        PHI=PHIM-REST
C        VTH=PHI/THRESHOLD
C        VPL=PHI/PLATEAU
C        IF(KTYP33.EQ.1) THEN      !Cubic
C          CUBIC_ION=MEMCOND*PHI*(1.d0-VTH)*(1.d0-VPL)
C        ELSE IF(KTYP33.EQ.2) THEN !Quintic
C          CUBIC_ION=MEMCOND*PHI*(1.d0-VTH*VTH)*(1.d0-VPL*VPL)
C        ELSE IF(KTYP33.EQ.3) THEN !Seventh-order
C          CUBIC_ION=MEMCOND*PHI*(1.d0-VTH*VTH*VTH)*(1.d0-VPL*VPL*VPL)
C        ENDIF
C      ENDIF
C
C      RETURN
C      END


C      SUBROUTINE FHN(Y,IONIC_PARAMS,F)
C
CC#### Subroutine: FHN
CC###  Description:
CC###    FHN is the FitzHugh-Nagumo ionic equation for
CC###    bidomain solution.
C
C      IMPLICIT NONE
C
C      INCLUDE 'cmiss$reference:b01.cmn'
C      INCLUDE 'cmiss$reference:b12.cmn'
C      INCLUDE 'cmiss$reference:cbdi02.cmn'
C      INCLUDE 'cmiss$reference:ktyp30.cmn'
C      INCLUDE 'cmiss$reference:tol00.cmn'
C
C!     Parameter List
C      REAL*8 F(*),IONIC_PARAMS(12),Y(*)
C!     Local Variables
C      REAL*8 DECAY,DIFF,EPS,FHN1,NORM_PHI,PHI,PHIM,PLATEAU,RATE,RECOV,
C     '  REST,THRESHOLD,VSTAR
C      CHARACTER ERROR*20
C
CC *** IONIC_PARAMS(1) is rest potential
CC *** IONIC_PARAMS(2) is plateau potential
CC *** IONIC_PARAMS(3) is threshold potential
CC *** IONIC_PARAMS(4) is excitation rate constant
CC *** IONIC_PARAMS(5) is excitation decay constant
CC ***  (6...) are recovery constants
C
C      REST=IONIC_PARAMS(1)
C      PLATEAU=IONIC_PARAMS(2)
C      THRESHOLD=IONIC_PARAMS(3)
C      RATE=IONIC_PARAMS(4)
C      DECAY=IONIC_PARAMS(5)
C
C      PHIM=Y(1)
C      RECOV=Y(2)
C
C      IF(DOP) THEN
C        WRITE(OP_STRING,'('' REST='',D12.4)') REST
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' PLATEAU='',D12.4)') PLATEAU
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' THRESHOLD='',D12.4)') THRESHOLD
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' RATE='',D12.4)') RATE
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' DECAY='',D12.4)') DECAY
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C      ENDIF
C
CC     Normalise wrt resting and plateau potentials
C      DIFF=PLATEAU-REST
C
C      PHI=(PHIM-REST)/DIFF
C      THRESHOLD=(THRESHOLD-REST)/DIFF
C
C      FHN1=RATE*PHI*(PHI-THRESHOLD)*(PHI-1.d0)
C
C      IF(DOP) THEN
C        WRITE(OP_STRING,'('' DIFF='',D12.4)') DIFF
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' THRESHOLD='',D12.4)') THRESHOLD
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' FHN='',D12.4)') FHN1
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C      ENDIF
C
C      IF(KTYP33.EQ.1) THEN !Standard FHN
C        FHN1=FHN1+DECAY*RECOV
C      ELSE !Roger's FHN or Panfilov FHN
C        FHN1=FHN1+DECAY*PHI*RECOV
C      ENDIF
C
C      F(1)=FHN1
C
C      IF(DOP) THEN
C        WRITE(OP_STRING,'('' FHN_ION='',D12.4)') FHN1
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C      ENDIF
C
C! recovery
C      PHI=Y(1)+F(1)*DT
C      DIFF=IONIC_PARAMS(2)-IONIC_PARAMS(1)
C      NORM_PHI=(PHI-IONIC_PARAMS(1))/DIFF
C      IF(KTYP33.LT.3) THEN !Standard or Roger's FHN
C        F(2)=IONIC_PARAMS(6)*(NORM_PHI-IONIC_PARAMS(7)*RECOV)
C      ELSE !Panfilov FHN
C        EPS=IONIC_PARAMS(6)+IONIC_PARAMS(7)*RECOV/(NORM_PHI+
C     '    IONIC_PARAMS(8))
C        RATE=IONIC_PARAMS(4)/DIFF/DIFF
C        VSTAR=PHI+IONIC_PARAMS(1)-IONIC_PARAMS(2)-IONIC_PARAMS(3)
C        F(2)=EPS*(-RECOV-RATE*VSTAR*(PHI-IONIC_PARAMS(1)))
C      ENDIF
C
C! calcium
C      IF(DABS(IONIC_PARAMS(9)).GT.ZERO_TOL) THEN
C        IF(NORM_PHI.LT.ZERO_TOL) NORM_PHI=0.0d0
C        F(3)=(NORM_PHI-Y(3))/IONIC_PARAMS(9)
C      ELSE
C        Y(3)=0.0d0
C        F(3)=0.0d0
C      ENDIF
C
C 9999 RETURN
C      END


C      SUBROUTINE VCD(Y,CQ,F)
C
CC#### Subroutine: VCD
CC###  Description:
CC###    VCD is the VanCapelle-Durrer ionic equation for bidomain
CC###    solution, with modifications from UCLA provided by Alan
CC###    Garfinkel. Original reference in fax dated March 8, 1994.
C
C      IMPLICIT NONE
C
C      INCLUDE 'cmiss$reference:b12.cmn'
C      INCLUDE 'cmiss$reference:geom00.cmn'
C      INCLUDE 'cmiss$reference:ktyp30.cmn'
C      INCLUDE 'cmiss$reference:tol00.cmn'
C
C!     Parameter List
C      REAL*8 CQ(NMM),F(*),Y(*)
C!     Local Variables
C      REAL*8 D_PHIM,F1,I0,I1,NORM_PHI,PHI,PHIM,RECOV,T_CONST
C
C      PHIM=Y(1)
C      RECOV=Y(2)
C      D_PHIM=Y(4)
C
C      IF(PHIM.LT.-70.0d0) THEN
C        I1=5.0d0+0.5d0*(PHIM+70.0d0)
C      ELSE IF(PHIM.GT.0.0d0) THEN
C        I1=6.0d0+0.425d0*PHIM
C      ELSE
C        I1=5.0d0+((PHIM+70.0d0)/70.0d0)
C      ENDIF
C
C      IF(PHIM.LT.-74.3d0) THEN
C        F1=7.84d0+2.0d0*(PHIM+74.3d0)
C      ELSE IF(PHIM.GT.-27.8d0) THEN
C        F1=-98.84d0+1.71d0*(PHIM+27.8d0)
C      ELSE
C        F1=((3.837854d-3*PHIM +0.584649d0)*PHIM +25.31834d0)*PHIM
C     '    +235.6256d0
C      ENDIF
C
C      IF(KTYP33.EQ.2) THEN  !Calif. mods to VCD
C        F1=4.0d0*F1
C        IF(KTYP34.EQ.2) THEN    !"Ischemic" APD (112ms)
C          IF(D_PHIM.LT.0) I1=2.0d0*I1
C        ENDIF
C      ENDIF
C
C      I0=I1+F1
C      F(1)=(RECOV*I1+(1.0d0-RECOV)*I0)*1.0d-6
C
C! recovery
C      PHI=PHIM+F(1)*DT
C      IF(KTYP33.EQ.1) THEN !Original VCD
C        T_CONST=CQ(10)
C      ELSE IF(KTYP33.EQ.2) THEN !VCDC mods
C        ! From modifications provided by Alan Garfinkel
C        IF(Y(5).GE.-ZERO_TOL) THEN
C          !recovery var Y is increasing
C          IF(KTYP34.EQ.1) THEN !"Normal" APD (209ms)
C            T_CONST=0.5d0
C          ELSE IF(KTYP34.EQ.2) THEN !"Ischemic" APD (112ms)
C            T_CONST=0.33d0
C          ENDIF
C        ELSE IF(RECOV.GT.0.85d0) THEN !Y decreasing & Y>0.85
C          IF(KTYP34.EQ.1) THEN !"Normal" APD (209ms)
C            T_CONST=0.1d0
C          ELSE IF(KTYP34.EQ.2) THEN !"Ischemic" APD (112ms)
C            T_CONST=0.066d0
C          ENDIF
C        ELSE                          !Y decreasing & Y<0.85
C          IF(KTYP34.EQ.1) THEN !"Normal" APD (209ms)
C            T_CONST=3.0d0
C          ELSE IF(KTYP34.EQ.2) THEN !"Ischemic" APD (112ms)
C            !PJH 5/7/98   T_CONST=3.31d0
C            T_CONST=CQ(14)
C          ENDIF
C        ENDIF
C        ! Scale Factor for time constant
C        T_CONST=T_CONST*CQ(10) !CQ(10) is time constant
C      ENDIF
C      IF(PHI.LT.-80.0d0) THEN
C        F(2)=-RECOV/T_CONST
C      ELSE IF(PHI.GT.-60.0d0) THEN
C        F(2)=(1.0d0-RECOV)/T_CONST
C      ELSE
C        F(2)=((PHI+80.0d0)/20.0d0-RECOV)/T_CONST
C      ENDIF
C      IF(KTYP33.EQ.2) THEN !VCDC mods
C        Y(5)=F(2)
C      ENDIF
C
C! calcium
C      IF(DABS(CQ(12)).GT.ZERO_TOL) THEN
C        NORM_PHI=(PHI-CQ(9))/100.0d0
C        IF(NORM_PHI.LT.ZERO_TOL) NORM_PHI=0.0d0
C        F(3)=(NORM_PHI-Y(3))/CQ(12)
C      ELSE
C        Y(3)=0.0d0
C        F(3)=0.0d0
C      ENDIF
C
C      RETURN
C      END


C      SUBROUTINE BR(Y,IONIC_PARAMS,F)
C
CC#### Subroutine: BR
CC###  Description:
CC###    BR is the Beeler-Reuter model, with the option of sodium
CC###    kinetics determined by Ebihara-Johnson or Drouhard-Roberge
CC###    according to KTYP33
C
C      IMPLICIT NONE
C
C      INCLUDE 'cmiss$reference:b12.cmn'
C      INCLUDE 'cmiss$reference:ktyp30.cmn'
C
C!     Parameter List
C      REAL*8 F(*),IONIC_PARAMS(12),Y(*)
C!     Local Variables
C      REAL*8 Ad,Af,Ah,Aj,Am,Ax1,Bd,Bf,Bh,Bj,Bm,Bx1,Cai,d,DCai,Dd,Df,Dh,
C     '  Dj,Dm,Dx1,f1,GNa,GNaC,Gs,h,IK1,INa,Is,Ix1,j,m,PHIM,V,Vs,VNa,x1
C
CC     IONIC_PARAMS(1) !resting potential
CC     IONIC_PARAMS(2) !sodium equilibrium (reversal) potential
CC     IONIC_PARAMS(3) !sodium conductance
CC     IONIC_PARAMS(4) !steady-state sodium conductance
CC     IONIC_PARAMS(5) !slow current conductance
C
C      PHIM=Y(1)
C      m=Y(2)
C      h=Y(3)
C      j=Y(4)
C      d=Y(5)
C      f1=Y(6)
C      x1=Y(7)
C      Cai=Y(8)
C
C      VNa=IONIC_PARAMS(2)
C      GNa=IONIC_PARAMS(3)
C      GNaC=IONIC_PARAMS(4)
C      Gs=IONIC_PARAMS(5)
C
C      INa=0.d0
C      IF(KTYP33.EQ.1) THEN
C        !Beeler-Reuter sodium model
C        INa = GNa*m*m*m*h*j + GNaC
C      ELSE IF(KTYP33.EQ.2.OR.KTYP33.GE.3) THEN
C        !Ebihara-Johnson or Drouhard-Roberge sodium model
C        INa = GNa*m*m*m*h
C      ENDIF !sodium model
C      INa = INa*(PHIM-VNa)
C
C      Vs = -82.3d0 - 13.0287d0*DLOG(Cai)
C      Is = Gs*d*f1*(PHIM-Vs)
C
CC     1.d-2 correction factor in the next two equations is for the
CC     change from cm^2 to mm^2
C
C      Ix1 = 0.8d0*x1*(DEXP(0.04d0*(PHIM+77.d0))-1.d0)/
C     '  (DEXP(0.04d0*(PHIM+35.d0)))!*1.d-2
C      IK1 = 0.35d0 * (4.d0*(DEXP(0.04d0*(PHIM+85.d0))-1.d0) /
C     '  (DEXP(0.08d0*(PHIM+53.d0))+DEXP(0.04d0*(PHIM+53.d0))) +
C     '  0.2d0*(PHIM+23.d0)/(1.d0-DEXP(-0.04d0*(PHIM+23.d0))))!*1.d-2
C
C      F(1)=(INa+Is+Ix1+IK1)*1.0d-6 !correction factor for [mV] & [mA]
C
C      V=PHIM+F(1)*DT !prediction for V
C
C      IF(KTYP33.EQ.1) THEN !Beeler-Reuter sodium model
C        Am = -1.d0*(V+47.d0)/(DEXP(-0.1d0*(V+47.d0))-1.d0)
C        Bm = 40.d0*DEXP(-0.056d0*(V+72.d0))
C        Ah = 0.126d0*DEXP(-0.25d0*(V+77.d0))
C        Bh = 1.7d0/(DEXP(-.082d0*(V+22.5d0))+1.d0)
C        Aj = 0.055d0*DEXP(-0.25d0*(V+78.d0))/
C     '    (DEXP(-0.2d0*(V+78.d0))+1.d0)
C        Bj = 0.3d0/(DEXP(-0.1d0*(V+32.d0))+1.d0)
C      ELSE IF(KTYP33.EQ.2) THEN !Ebihara-Johnson sodium model
C        Am = (0.32d0*(V+47.13d0))/(1.d0-DEXP(-V-47.13d0))
C        Bm = 0.08d0*DEXP(-V/11.d0)
C        IF(V.GE.-40.d0) THEN
C          Ah = 0.d0
C          Bh = 1.d0/(0.13d0*(DEXP((V+10.66d0)/(-11.1d0))+1.d0))
C        ELSE
C          Ah = 0.135d0*DEXP((-80.d0-V)/6.8d0)
C          Bh = 3.56d0*DEXP(0.079d0*V)+3.1d5*DEXP(0.35d0*V)
C        ENDIF
C        Aj = 0.d0
C        Bj = 0.d0
C      ELSE IF(KTYP33.EQ.3) THEN !Drouhard-Roberge sodium model
C        Am = 0.9d0*(V+42.65d0)/(1.d0-DEXP(-0.22d0*(V+42.65d0)))
C        Bm = 1.437d0*DEXP(-0.085d0*(V+39.75d0))
C        Ah = 0.1d0*DEXP(-0.193d0*(V+79.65d0))
C        Bh = 1.7d0/(1.d0+DEXP(-0.095d0*(V+20.5d0)))
C        Aj = 0.d0
C        Bj = 0.d0
C      ELSE IF(KTYP33.EQ.4) THEN !Modified Drouhard-Roberge Na model
C        !Alpha_m and Beta_m
C        IF(V.LE.100.0d0) THEN
C          Am=0.9d0*((V+42.65d0)/(1.0d0-DEXP((-0.22d0*V)-9.383d0)))
C        ELSE
C          Am=890.943789d0*((DEXP((0.0486479d0*V)-4.8647916d0))/
C     '      (1.0d0+(5.93962526d0*DEXP((0.0486479d0*V)-4.8647916d0))))
C        ENDIF
C        IF(V.LE.-85.0d0) THEN
C          Bm=100.0d0/(1.0d0+(0.4864082d0*DEXP((0.2597504d0*V)+
C     '      22.0787804d0)))
C        ELSE
C          Bm=1.437d0*DEXP((-0.085d0*V)-3.37875d0)
C        ENDIF
C
C        !Alpha_h and Beta_h
C        IF(V.LE.-90.0d0) THEN
C          Ah=-12.0662845d0-(0.1422598d0*V)
C        ELSE
C          Ah=0.1d0*DEXP((-0.193d0*V)-15.37245d0)
C        ENDIF
C        Bh=1.7d0/(1.0d0+DEXP((-0.095d0*V)-1.9475d0))
C        Aj = 0.d0
C        Bj = 0.d0
C      ENDIF !sodium model
C      Dm=Am*(1.d0-m)-Bm*m
C      Dh=Ah*(1.d0-h)-Bh*h
C      Dj=Aj*(1.d0-j)-Bj*j
C
C      Vs = -82.3d0 - 13.0287d0*DLOG(Cai)
C      Is = Gs*d*f1*(V-Vs)
C      IF(KTYP33.EQ.4) THEN
C        !d[Ca]i/dt
C        IF(V.LE.200.0d0) THEN
C          DCai=-1.d-7*Is + 0.07d0*(1.d-7-Cai)
C        ELSE
C          DCai=0.0d0
C        ENDIF
C      ELSE
C        DCai=-1.d-7*Is + 0.07d0*(1.d-7-Cai)
C      ENDIF
C
C      Ad = 0.095d0*DEXP(-0.01d0*(V-5.d0))/
C     '  (1.d0+DEXP(-0.072d0*(V-5.d0)))
C      Bd = 0.07d0*DEXP(-(V+44.d0)/59.d0)/
C     '  (1.d0+DEXP(0.05d0*(V+44.d0)))
C      Af = 0.012d0*DEXP(-0.008d0*(V+28.d0))/
C     '  (1.d0+DEXP(0.15d0*(V+28.d0)))
C      Bf = 0.0065d0*DEXP(-0.02d0*(V+30.d0))/
C     '  (1.d0+DEXP(-0.2d0*(V+30.d0)))
C      Dd=Ad*(1.d0-d)-Bd*d
C      Df=Af*(1.d0-f1)-Bf*f1
C
C      IF(KTYP33.EQ.4) THEN
C        !Alpha_x1 and Beta_x1
C        IF(V.LE.400.0d0) THEN
C          Ax1=0.0005d0*(DEXP((0.083d0*V)+4.150d0)/(DEXP((0.057d0*V)+
C     '      2.85d0)+1.0d0))
C        ELSE
C          Ax1=151.7994692d0*(DEXP((0.0654679d0*V)-26.1871448d0)/
C     '      (1.0d0+(1.5179947d0*DEXP((0.0654679d0*V)-26.1871448d0))))
C        ENDIF
C        Bx1=0.0013d0*(DEXP((-0.06d0*V)-1.2d0)/((DEXP((-0.04d0*V)
C     '    -0.8d0))+1.0d0))
C      ELSE
C        Ax1 = 5.d-4*DEXP(-(V+50.d0)/12.1d0)/(1.d0+DEXP((V+50.d0)/
C     '    17.5d0))
C        Bx1 = 0.0013d0*DEXP(-0.06d0*(V+20.d0))/
C     '    (1.d0+DEXP(-0.04d0*(V+20.d0)))
C      ENDIF
C      Dx1=Ax1*(1.d0-x1)-Bx1*x1
C
C      F(2)=Dm*1.d-1
C      F(3)=Dh*1.d-1
C      F(4)=Dj*1.d-1
C      F(5)=Dd*1.d-1
C      F(6)=Df*1.d-1
C      F(7)=Dx1*1.d-1
C      F(8)=DCai*1.d-1
C
C      RETURN
C      END


! ??? DPN 19/02/98 CAREV unused ???
c      SUBROUTINE CAREV(Y)

C#### Subroutine: CAREV
C###  Description:
C###    Finds ICa reversal potential to within 0.1mv

c      IMPLICIT NONE
c      INCLUDE 'cmiss$reference:deoxs00.cmn'
c      INCLUDE 'cmiss$reference:oxs003.cmn'
!     Parameter List
c      REAL*8 Y(*)
!     Local Variables
c      REAL*8 Z2,Z4,Z5 !,EDF,Z3

!      write(*,'('' >>>call carev'')')

c      IF(PCA.NE.0) THEN
c        Z4 = Y(4)*DEXP(2.0d0*VSURFCA/RTONF)
c        Z5 = Y(11)*DEXP(VSURFCA/RTONF)
c        Z1 = PCAK*PCA*Z5
c        Z2 = PCAK*PCA*Y(3)
C dpn 17/02/98        Z3 = (Z2-Z1)*(Z2-Z1)+4.0d0*(Z1+4.0d0*PCA*Z4)*(Z2+4.0d0*PCA*CAO)
c        EDF= RTONF*DLOG((Z2-Z1+SQRT(Z3))/(2.0d0*(Z1+4.0d0*PCA*Z4)))
c     '    +VSURFCA
c      ELSE
c        EDF= 0.0d0
c      ENDIF

c      RETURN
c      END


