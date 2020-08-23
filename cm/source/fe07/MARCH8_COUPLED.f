      SUBROUTINE MARCH8_COUPLED(INTEGRATOR_IWORK,CELL_ICQS_VALUE,IBT,
     &  ICQS,ICQS_SPATIAL,IDO,IICQS_SPATIAL,IRCQS_SPATIAL,INP,ISC_GKK,
     &  ISR_GKK,NBH,NBJ,NEELEM,NENQ,NLATNE,NLATNQ,NLATPNQ,
     &  NQNLAT,NHE,NHP,NKHE,NKH,NLQ,NPF,NPNE,NP_INTERFACE,NPLIST,NPNODE,
     &  NQGP,NQGP_PIVOT,NQLIST,NQNP,NQS,NQSCNB,NQXI,NRLIST,NTIME_INTERP,
     &  NTIME_POINTS,NVHE,NVHP,NW,NWQ,NXLIST,NXQ,NYNE,NYNP,
     &  INTEGRATOR_WORK,ALPHA,AQ,CELL_RCQS_VALUE,CQ,CURVCORRECT,DNUDXQ,
     &  DXDXIQ,DXDXIQ2,GCHQ,GK,GK2,GK3,GKK,GM,GQ,GQ2,GQ3,GUQ,
     &  NQGW,PG,PROPQ,RCQS,RCQS_SPATIAL,RHS,SE,TIME_VALUES,WG,XIQ,XQ,YP,
     &  YQ,YQS,ZA,ZE,ZP,FIXQ,RHSROUTINE,ERROR,*)

C#### Subroutine: MARCH8_COUPLED
C###  Description:
C###    MARCH8_COUPLED is another time stepping routine for more
C###     specific grid ionic current activation problems.
C***  Created by Martin Buist, February 2000

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'adam00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'cell02.cmn'
      INCLUDE 'cellml.cmn'
      INCLUDE 'cell_reserved.inc'
      INCLUDE 'cell_vcd.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'integrator.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'loc00.cmn'
      INCLUDE 'lsoda00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'maqloc00.inc'
      INCLUDE 'marc00.cmn'
      INCLUDE 'nqloc00.inc'
      INCLUDE 'solv00.cmn'
      INCLUDE 'time02.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter list
      INTEGER INTEGRATOR_IWORK(INTEGRATOR_LIWORK,NQM),
     &  CELL_ICQS_VALUE(NQIM,NQVM),IBT(3,NIM,NBFM),ICQS(NQIM),
     &  ICQS_SPATIAL(NQISVM,NQM),IDO(NKM,NNM,0:NIM,NBFM),
     &  IICQS_SPATIAL(0:NQISVM,NQVM),IRCQS_SPATIAL(0:NQRSVM,NQVM),
     &  INP(NNM,NIM,NBFM),ISC_GKK(NISC_GKKM,NXM),ISR_GKK(NISR_GKKM,NXM),
     &  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     &  NEELEM(0:NE_R_M,0:NRM),NENQ(0:8,NQM),NLATNE(NEQM+1),
     &  NLATNQ(NEQM*NQEM),NLATPNQ(NQM),NQNLAT(NEQM*NQEM),NHE(NEM,NXM),
     &  NHP(NPM,0:NRM,NXM),NKHE(NKM,NNM,NHM,NEM),NKH(NHM,NPM,NCM,0:NRM),
     &  NLQ(NQM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),NP_INTERFACE(0:NPM,0:3),
     &  NPLIST(0:NPM),NPNODE(0:NP_R_M,0:NRM),NQGP(0:NQGM,NQM),
     &  NQGP_PIVOT(NQGM,NQM),NQLIST(0:NQM),NQNP(NPM),NQS(NEQM),
     &  NQSCNB(NQSCM),NQXI(0:NIM,NQSCM),NRLIST(0:NRM),
     &  NTIME_INTERP(NTIMEVARSM),NTIME_POINTS(NTIMEVARSM),
     &  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),NW(NEM,3),
     &  NWQ(8,0:NQM,NAM),NXLIST(0:NXM),NXQ(-NIM:NIM,0:4,0:NQM,NAM),
     &  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     &  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 INTEGRATOR_WORK(INTEGRATOR_LWORK,NQM),ALPHA,AQ(NMAQM,NQM),
     &  CELL_RCQS_VALUE(NQRM,NQVM),CQ(NMM,NQM),CURVCORRECT(2,2,NNM,NEM),
     &  DNUDXQ(3,3,NQM),DXDXIQ(3,3,NQM),DXDXIQ2(3,3,NQM),GCHQ(3,NQM),
     &  GK(NZ_GK_M),GK2(*),GK3(*),GKK(NZ_GKK_M,NXM),GM(NZ_GM_M),
     &  GQ(NZ_GQ_M),GQ2(*),GQ3(*),GUQ(3,3,NQM),
     &  NQGW(NQGM,NQM),PG(NSM,NUM,NGM,NBM),PROPQ(3,3,4,2,NQM),
     &  RCQS(NQRM),RCQS_SPATIAL(NQRSVM,NQM),RHS(NQM),SE(NSM,NBFM,NEM),
     &  TIME_VALUES(2,0:NTIMEPOINTSM+1,NTIMEVARSM),WG(NGM,NBM),
     &  XIQ(NIM,NQM),XQ(NJM,NQM),YP(NYM,NIYM,NXM),YQ(NYQM,NIQM,NAM,NXM),
     &  YQS(NIQSM,NQM),ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),
     &  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
      LOGICAL FIXQ(NYQM,NIYFIXM,NXM)
!     Local variables
      INTEGER atime,dpot,LOCAL_ICQS(NQIM),maqp1i,maqp1t0,maqp1t1,maqp2i,
     &  maqp2t0,maqp2t1,maqdt,nc,nh,niq_old,niqoldDV,niqPSTIM,niqV,
     &  niqBNDRY,nmq,no2,np,nq,nqq,nr,nr2,nr3,nr4,NRTEMP(0:2),NSIZE,
     &  NSTEP,nv,nx,nx2,nx3,nx4,nxc,nx_ext,NX_TORSO_LIST(9),nxx,ny,
     &  SIZES(11)
      INTEGER*4 NQLIST2_PTR,NQLIST3_PTR,NQLIST4_PTR,NQNP_LOCAL_PTR,
     &  NQNP_LOCAL_PIV_PTR,D_MATRIX_PTR,GQ_D_MATRIX_PTR,GKGK2_PTR,
     &  GQGQ2_PTR,VARIANT
      REAL*8 CMAMDT,GDPHIDN,LOCAL_RCQS(NQRM),NEWDVM,P1,PSTIMULUS,
     &  STIMULUS,SUM1,SUM2,T
      REAL ELAPSED_TIME,TIME_START(1),TIME_START2(1),TIME_STOP(1),
     &  TIME_STOP2(1)
      LOGICAL CONTINU,ERROR_FLAG,FIRST_A,FIRST_EXT_A,UPDATE_EXT_MATRIX,
     &  UPDATE_MATRIX,USESALU,X_INIT

      EXTERNAL RHSROUTINE

      SAVE maqdt,maqp1t0,maqp1t1,maqp1i,maqp2t0,maqp2t1,maqp2i,niq_old,
     &  niqV,niqBNDRY,niqPSTIM,niqoldDV,NSIZE,NSTEP,SIZES

      CALL ENTERS('MARCH8_COUPLED',*9999)

C**** CPTYPE=1 is iteration on epi & endo
C**** CPTYPE=2 is coupled grid-endo, no iteration
C**** CPTYPE=3 is coupled grid-endo, epi iteration
C**** CPTYPE=4 is coupled grid-2 endo, no iteration
C**** CPTYPE=5 is coupled grid-2 endo, epi iteration
C**** CPTYPE=6 is coupled epi-grid-endo, no iteration
C**** CPTYPE=7 is coupled epi-grid-2 endo, no iteration

      DT=TINCR
      nr=NRLIST(1)
      CONTINU=.TRUE.
      USESALU=.TRUE.

      !Get class information
      CALL ASSERT(NXLIST(0).GE.2,'>>Error - bidomain needs 2 classes',
     &  ERROR,*9999)

      nxc=NXLIST(1)
      CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
      CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     &  ERROR,*9999)

      nxc=NXLIST(2)
      CALL NX_LOC(NX_INQUIRE,nxc,nx_ext,NX_SOLVE,ERROR,*9999)
      CALL ASSERT(nx_ext.GT.0,'>>No nx defined for this solve class',
     &  ERROR,*9999)

      LOCAL_ICQS(1)=1 !variant always equals 1?
      DO nmq=2,NQIT
        LOCAL_ICQS(nmq)=ICQS(nmq)
      ENDDO
      DO nmq=NQIT+1,NQIM
        LOCAL_ICQS(nmq)=0
      ENDDO
      DO nmq=1,NQRT
        LOCAL_RCQS(nmq)=RCQS(nmq)
      ENDDO
      DO nmq=NQRT+1,NQRM
        LOCAL_RCQS(nmq)=0.0d0
      ENDDO

      !Check if first time through
      IF(RESTART) THEN
        T=T_RESTART(nx)
        FIRST_A=.FALSE.
        FIRST_EXT_A=.FALSE.
        UPDATE_EXT_MATRIX=.FALSE.
        UPDATE_MATRIX=.FALSE.

        IF((CPTYPE.EQ.2).OR.(CPTYPE.EQ.3)) THEN
          nr2=NRLIST(2)
          nxc=NXLIST(3)
          CALL NX_LOC(NX_INQUIRE,nxc,nx2,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx2.GT.0,'>>No nx defined for this solve class',
     &      ERROR,*9999)
        ENDIF
        IF((CPTYPE.EQ.4).OR.(CPTYPE.EQ.5)) THEN
          nr2=NRLIST(2)
          nxc=NXLIST(3)
          CALL NX_LOC(NX_INQUIRE,nxc,nx2,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx2.GT.0,'>>No nx defined for this solve class',
     &      ERROR,*9999)
          nr3=NRLIST(3)
          nxc=NXLIST(4)
          CALL NX_LOC(NX_INQUIRE,nxc,nx3,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx3.GT.0,'>>No nx defined for this solve class',
     &      ERROR,*9999)
        ENDIF
        IF(CPTYPE.EQ.6) THEN
          nr2=NRLIST(2)
          nxc=NXLIST(3)
          CALL NX_LOC(NX_INQUIRE,nxc,nx2,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx2.GT.0,'>>No nx defined for this solve class',
     &      ERROR,*9999)
          nr3=NRLIST(3)
          nxc=NXLIST(4)
          CALL NX_LOC(NX_INQUIRE,nxc,nx3,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx3.GT.0,'>>No nx defined for this solve class',
     &      ERROR,*9999)
        ENDIF
        IF(CPTYPE.EQ.7) THEN
          nr3=NRLIST(3)
          nxc=NXLIST(4)
          CALL NX_LOC(NX_INQUIRE,nxc,nx3,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx3.GT.0,'>>No nx defined for this solve class',
     &      ERROR,*9999)
          nr4=NRLIST(4)
          nxc=NXLIST(5)
          CALL NX_LOC(NX_INQUIRE,nxc,nx4,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx4.GT.0,'>>No nx defined for this solve class',
     &      ERROR,*9999)
        ENDIF
      ELSE
        T=TSTART
        X_INIT=.TRUE.
        NSTEP=0

        !Init time variable NWQ to 0
        DO nq=1,NQT
          NWQ(8,nq,1)=0
        ENDDO

        !Get indices for AQ auxiliary indices
        CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE1,maqp1t0,MAQ_START,
     &    ERROR,*9999)
        CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE1,maqp1t1,MAQ_STOP,
     &    ERROR,*9999)
        CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE1,maqp1i,MAQ_CURRENT,
     &    ERROR,*9999)
        CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE2,maqp2t0,MAQ_START,
     &    ERROR,*9999)
        CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE2,maqp2t1,MAQ_STOP,
     &    ERROR,*9999)
        CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE2,maqp2i,MAQ_CURRENT,
     &    ERROR,*9999)
        CALL MAQ_LOC(MAQ_INQUIRE,MAQ_TIME,maqdt,MAQ_DT,
     &    ERROR,*9999)

        !Get indices for NIQ solution indices
        CALL ASSERT(NIQM.GE.6,'>>Increase NIQM >=6',ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niq_old,NIQ_OLDSOLN,
     &    ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqV,NIQ_V,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqBNDRY,NIQ_BNDRY,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqPSTIM,NIQ_OLDSOLN2,
     &    ERROR,*9999)
        IF(niqPSTIM.EQ.0) THEN
          CALL NIQ_LOC(NIQ_ALLOCATE_AND_LOCK,NIQ_ION,niqPSTIM,
     &      NIQ_OLDSOLN2,ERROR,*9999)
        ENDIF
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqoldDV,NIQ_X,
     &    ERROR,*9999)
        IF(niqoldDV.EQ.0) THEN
          CALL NIQ_LOC(NIQ_ALLOCATE_AND_LOCK,NIQ_ION,niqoldDV,
     &      NIQ_X,ERROR,*9999)
        ENDIF

        !Get activation time indices and initialize array
        CALL MAQ_LOC(MAQ_INQUIRE,MAQ_TIME,atime,MAQ_ACTIV_TIME,
     &      ERROR,*9999)
        CALL MAQ_LOC(MAQ_INQUIRE,MAQ_TIME,dpot,MAQ_M_DPOT,
     &      ERROR,*9999)
        DO nq=1,NQT
          AQ(atime,nq)=0.0d0
          AQ(dpot,nq)=0.0d0
        ENDDO

C ***   Perform any first call initialisation required by the various
C       integrators
        IF(KTYP37.EQ.5) THEN !Adams-Moulton integrator
C$OMP     PARALLEL DO
C$OMP&    PRIVATE(nq),
C$OMP&    SHARED(INTEGRATOR_IWORK,INTEGRATOR_WORK,DT)
          DO nq=1,NQT
            INTEGRATOR_IWORK(1,nq)=0
            INTEGRATOR_WORK(1,nq)=DT
          ENDDO
C$OMP     END PARALLEL DO
        ELSEIF(KTYP37.EQ.6) THEN ! LSODA integrator
C$OMP     PARALLEL
C$OMP&      PRIVATE(nq),
C$OMP&      SHARED(INTEGRATOR_IWORK,INTEGRATOR_WORK)
C$OMP     DO
          DO nq=1,NQT
            INTEGRATOR_IWORK(LSODA_ISTATE,nq)=1
          ENDDO
C$OMP     END PARALLEL
        ENDIF

        !Assemble matrices first time through
        IF(CPTYPE.EQ.1) THEN
          NRTEMP(0)=1
          NRTEMP(1)=NRLIST(1)
C SGM 13 Nov 2000 call ASSEMBLE10_FE for grid-based FE
C MLT 29Nov02 Same for grid FV
          IF(ITYP4(NRLIST(1),nx).EQ.4) THEN !collocation
            CALL ASSEMBLE10(ISC_GKK,ISR_GKK,NEELEM,NENQ,NLATNE,NLATNQ,
     &        NLATPNQ,NLQ,NQGP,NQGP_PIVOT,NQNLAT,NQS,NQXI,NRTEMP,
     &        NWQ(1,0,1),nx_ext,nx,0,NXQ,AQ,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,
     &        GCHQ,GUQ,GKK,NQGW,PROPQ,XQ,.TRUE.,.FALSE.,FIRST_A,FIXQ,
     &        .TRUE.,.FALSE.,.FALSE.,UPDATE_MATRIX,.FALSE.,ERROR,*9999)
          ELSEIF(ITYP4(NRLIST(1),nx).EQ.6.OR.
     &           ITYP4(NRLIST(1),nx).EQ.7) THEN !grid-based FE/grid FV
            CALL ASSEMBLE10_FE(ISC_GKK,ISR_GKK,NEELEM,NLATNE,NQGP,
     &        NQGP_PIVOT,NQNLAT,NQS,NQSCNB,NQXI,NRTEMP,nx_ext,nx,0,NXQ,
     &        CQ,GKK,GM,NQGW,PG,WG,XQ,PROPQ,.TRUE.,.FALSE.,
     &        FIRST_A,FIXQ,.TRUE.,.FALSE.,.FALSE.,UPDATE_MATRIX,.FALSE.,
     &        ERROR,*9999)
          ENDIF
          FIRST_A=.TRUE.
          FIRST_EXT_A=.TRUE.
          UPDATE_EXT_MATRIX=.TRUE.
          UPDATE_MATRIX=.TRUE.
        ELSE IF((CPTYPE.EQ.2).OR.(CPTYPE.EQ.3).OR.(CPTYPE.EQ.4)
     &    .OR.(CPTYPE.EQ.5).OR.(CPTYPE.EQ.6).OR.(CPTYPE.EQ.7)) THEN
          NRTEMP(0)=1
          NRTEMP(1)=NRLIST(1)
C SGM 13 Nov 2000 call ASSEMBLE10_FE for grid-based FE
C MLT 29Nov02 Same for grid FV
          IF(ITYP4(NRLIST(1),nx).EQ.4) THEN !collocation
            CALL ASSEMBLE10(ISC_GKK,ISR_GKK,NEELEM,NENQ,NLATNE,NLATNQ,
     &        NLATPNQ,NLQ,NQGP,NQGP_PIVOT,NQNLAT,NQS,NQXI,NRTEMP,
     &        NWQ(1,0,1),nx_ext,nx,0,NXQ,AQ,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,
     &        GCHQ,GUQ,GKK,NQGW,PROPQ,XQ,.FALSE.,.FALSE.,FIRST_A,
     &        FIXQ,.TRUE.,.FALSE.,.FALSE.,UPDATE_MATRIX,.FALSE.,
     &        ERROR,*9999)
          ELSEIF(ITYP4(NRLIST(1),nx).EQ.6.OR.
     &           ITYP4(NRLIST(1),nx).EQ.7) THEN !grid-based FE/grid FV
            CALL ASSEMBLE10_FE(ISC_GKK,ISR_GKK,NEELEM,NLATNE,NQGP,
     &        NQGP_PIVOT,NQNLAT,NQS,NQSCNB,NQXI,NRTEMP,nx_ext,nx,0,NXQ,
     &        CQ,GKK,GM,NQGW,PG,WG,XQ,PROPQ,.FALSE.,.FALSE.,
     &        FIRST_A,FIXQ,.TRUE.,.FALSE.,.FALSE.,UPDATE_MATRIX,.FALSE.,
     &        ERROR,*9999)
          ENDIF
          FIRST_A=.TRUE.
          FIRST_EXT_A=.TRUE.
          UPDATE_EXT_MATRIX=.TRUE.
          UPDATE_MATRIX=.TRUE.

          nr2=NRLIST(2)
          nxc=NXLIST(3)
          CALL NX_LOC(NX_INQUIRE,nxc,nx2,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx2.GT.0,'>>No nx defined for this solve class',
     &      ERROR,*9999)
          IF((CPTYPE.EQ.4).OR.(CPTYPE.EQ.5)) THEN
            nr3=NRLIST(3)
            nxc=NXLIST(4)
            CALL NX_LOC(NX_INQUIRE,nxc,nx3,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx3.GT.0,'>>No nx defined for this solve class',
     &        ERROR,*9999)
          ELSE IF(CPTYPE.EQ.6) THEN
            nr3=NRLIST(3)
            nxc=NXLIST(4)
            CALL NX_LOC(NX_INQUIRE,nxc,nx3,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx3.GT.0,'>>No nx defined for this solve class',
     &        ERROR,*9999)
          ELSE IF(CPTYPE.EQ.7) THEN
            nr3=NRLIST(3)
            nxc=NXLIST(4)
            CALL NX_LOC(NX_INQUIRE,nxc,nx3,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx3.GT.0,'>>No nx defined for this solve class',
     &        ERROR,*9999)
            nr4=NRLIST(4)
            nxc=NXLIST(5)
            CALL NX_LOC(NX_INQUIRE,nxc,nx4,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx4.GT.0,'>>No nx defined for this solve class',
     &      ERROR,*9999)
          ENDIF

          !Get interface grid points
          CALL PARSE_GRID(NQLIST,noco,NTCO,CO,ERROR,*9999)

          NQLIST2_PTR=0
          NQLIST3_PTR=0
          NQLIST4_PTR=0
          NQNP_LOCAL_PTR=0
          NQNP_LOCAL_PIV_PTR=0
          D_MATRIX_PTR=0
          GQ_D_MATRIX_PTR=0
          GKGK2_PTR=0
          GQGQ2_PTR=0

          CALL ALLOCATE_MEMORY(NQM+1,1,INTTYPE,NQLIST2_PTR,
     &      MEM_INIT,ERROR,*9999)
          CALL ALLOCATE_MEMORY(NQM+1,1,INTTYPE,NQLIST3_PTR,
     &      MEM_INIT,ERROR,*9999)
          CALL ALLOCATE_MEMORY(NQM+1,1,INTTYPE,NQLIST4_PTR,
     &      MEM_INIT,ERROR,*9999)
          CALL ALLOCATE_MEMORY(NPM,1,INTTYPE,NQNP_LOCAL_PTR,
     &      MEM_INIT,ERROR,*9999)
          CALL ALLOCATE_MEMORY(NPM,1,INTTYPE,NQNP_LOCAL_PIV_PTR,
     &      MEM_INIT,ERROR,*9999)
          CALL ALLOCATE_MEMORY(NPM*NPM*3,1,DPTYPE,D_MATRIX_PTR,
     &      MEM_INIT,ERROR,*9999)
          CALL ALLOCATE_MEMORY(NPM*NPM*3,1,DPTYPE,GQ_D_MATRIX_PTR,
     &      MEM_INIT,ERROR,*9999)
          CALL ALLOCATE_MEMORY(NPM*NPM,1,DPTYPE,GKGK2_PTR,
     &      MEM_INIT,ERROR,*9999)
          CALL ALLOCATE_MEMORY(NPM*NPM,1,DPTYPE,GQGQ2_PTR,
     &      MEM_INIT,ERROR,*9999)

          IF((CPTYPE.EQ.2).OR.(CPTYPE.EQ.3)) THEN
            CALL ASSEMBLE11_2_3(ISC_GKK,ISR_GKK,NEELEM,NENQ,NPNODE,
     &        NQGP,NQGP_PIVOT,NQLIST,%VAL(NQLIST2_PTR),
     &        %VAL(NQLIST3_PTR),NQNP,%VAL(NQNP_LOCAL_PTR),
     &        %VAL(NQNP_LOCAL_PIV_PTR),NQS,NQXI,nr,nr2,NWQ,nx_ext,
     &        nx2,NXQ,CQ,%VAL(D_MATRIX_PTR),DNUDXQ,DXDXIQ,DXDXIQ2,
     &        GCHQ,GK,GKK,GQ,%VAL(GQ_D_MATRIX_PTR),
     &        %VAL(GKGK2_PTR),%VAL(GQGQ2_PTR),GUQ,NQGW,PROPQ,FIXQ,
     &        .FALSE.,ERROR,*9999)
          ELSE IF((CPTYPE.EQ.4).OR.(CPTYPE.EQ.5)) THEN
            CALL ASSEMBLE11_4_5(ISC_GKK,ISR_GKK,NEELEM,NENQ,NPNODE,
     &        NQGP,NQGP_PIVOT,NQLIST,%VAL(NQLIST2_PTR),
     &        %VAL(NQLIST3_PTR),%VAL(NQLIST4_PTR),
     &        NQNP,%VAL(NQNP_LOCAL_PTR),
     &        %VAL(NQNP_LOCAL_PIV_PTR),NQS,NQXI,nr,nr2,nr3,NWQ,nx_ext,
     &        nx2,nx3,NXQ,CQ,%VAL(D_MATRIX_PTR),DNUDXQ,DXDXIQ,DXDXIQ2,
     &        GCHQ,GK,GK2,GKK,GQ,GQ2,%VAL(GQ_D_MATRIX_PTR),
     &        %VAL(GKGK2_PTR),%VAL(GQGQ2_PTR),GUQ,NQGW,PROPQ,FIXQ,
     &        .FALSE.,ERROR,*9999)
          ELSE IF(CPTYPE.EQ.6) THEN
            CALL ASSEMBLE11_6(ISC_GKK,ISR_GKK,NEELEM,NENQ,
     &        NP_INTERFACE,NPLIST,NPNODE,NQGP,NQGP_PIVOT,NQLIST,
     &        %VAL(NQLIST2_PTR),%VAL(NQLIST3_PTR),NQNP,
     &        %VAL(NQNP_LOCAL_PTR),%VAL(NQNP_LOCAL_PIV_PTR),NQS,NQXI,
     &        nr,nr2,nr3,NSIZE,NWQ,nx_ext,nx2,nx3,NXQ,CQ,
     &        %VAL(D_MATRIX_PTR),DNUDXQ,DXDXIQ,DXDXIQ2,GCHQ,GK,GK2,
     &        GKK,GQ,GQ2,%VAL(GQ_D_MATRIX_PTR),%VAL(GKGK2_PTR),
     &        %VAL(GQGQ2_PTR),GUQ,NQGW,PROPQ,YQ,FIXQ,.FALSE.,
     &        ERROR,*9999)
          ELSE IF(CPTYPE.EQ.7) THEN
            CALL ASSEMBLE11_7(ISC_GKK,ISR_GKK,NEELEM,NENQ,
     &        NP_INTERFACE,NPLIST,NPNODE,NQGP,NQGP_PIVOT,NQLIST,
     &        %VAL(NQLIST2_PTR),%VAL(NQLIST3_PTR),%VAL(NQLIST4_PTR),
     &        NQNP,%VAL(NQNP_LOCAL_PTR),%VAL(NQNP_LOCAL_PIV_PTR),NQS,
     &        NQXI,nr,nr2,nr3,nr4,NSIZE,NWQ,nx_ext,nx2,nx3,nx4,NXQ,CQ,
     &        %VAL(D_MATRIX_PTR),DNUDXQ,DXDXIQ,DXDXIQ2,GCHQ,GK,GK2,GK3,
     &        GKK,GQ,GQ2,GQ3,%VAL(GQ_D_MATRIX_PTR),%VAL(GKGK2_PTR),
     &        %VAL(GQGQ2_PTR),GUQ,NQGW,PROPQ,YQ,FIXQ,.FALSE.,
     &        ERROR,*9999)
          ENDIF

          CALL FREE_MEMORY(NQLIST2_PTR,ERROR,*9999)
          CALL FREE_MEMORY(NQLIST3_PTR,ERROR,*9999)
          CALL FREE_MEMORY(NQLIST4_PTR,ERROR,*9999)
          CALL FREE_MEMORY(NQNP_LOCAL_PTR,ERROR,*9999)
          CALL FREE_MEMORY(NQNP_LOCAL_PIV_PTR,ERROR,*9999)
          CALL FREE_MEMORY(D_MATRIX_PTR,ERROR,*9999)
          CALL FREE_MEMORY(GQ_D_MATRIX_PTR,ERROR,*9999)
          CALL FREE_MEMORY(GKGK2_PTR,ERROR,*9999)
          CALL FREE_MEMORY(GQGQ2_PTR,ERROR,*9999)

          IF(USESALU) THEN
            IF(CPTYPE.EQ.6) THEN
              DO nq=1,NSIZE
                RHS(nq)=1.0d0
              ENDDO
              CALL SOLVE_SYSTEM(ISC_GKK(1,nx_ext),ISR_GKK(1,nx_ext),
     &          NSIZE,NSIZE,NSIZE,NZZT(1,NRLIST(1),nx_ext),
     &          IWRIT4(NRLIST(1),nx_ext),PRECON_CODE(nx_ext),
     &          SOLVEROPTION(nx_ext),SPARSEGKK(nx_ext),GKK(1,nx_ext),
     &          RHS,YQ(1,niqV,1,nx3),FIRST_EXT_A,UPDATE_EXT_MATRIX,
     &          X_INIT,nx_ext,ERROR,*9999)
              UPDATE_EXT_MATRIX=.FALSE.
            ENDIF
            IF(CPTYPE.EQ.7) THEN
              DO nq=1,NSIZE
                RHS(nq)=1.0d0
              ENDDO
              CALL SOLVE_SYSTEM(ISC_GKK(1,nx_ext),ISR_GKK(1,nx_ext),
     &          NSIZE,NSIZE,NSIZE,NZZT(1,NRLIST(1),nx_ext),
     &          IWRIT4(NRLIST(1),nx_ext),PRECON_CODE(nx_ext),
     &          SOLVEROPTION(nx_ext),SPARSEGKK(nx_ext),GKK(1,nx_ext),
     &          RHS,YQ(1,niqV,1,nx4),FIRST_EXT_A,UPDATE_EXT_MATRIX,
     &          X_INIT,nx_ext,ERROR,*9999)
              UPDATE_EXT_MATRIX=.FALSE.
            ENDIF
          ENDIF
        ENDIF
      ENDIF

C$OMP PARALLEL DO
C$OMP&PRIVATE(nq),
C$OMP&SHARED(niq_OLD,niqV,nr,nx,YQ)
      DO nq=NQR(1,nr),NQR(2,nr)
        YQ(nq,niqV,1,nx)=YQ(nq,niq_OLD,1,nx)
      ENDDO
C$OMP END PARALLEL DO

      CALL CPU_TIMER(CPU_USER,TIME_START)
      CALL CPU_TIMER(CPU_USER,TIME_START2)

      !Main time integration loop
      DO WHILE(CONTINU)
        IF(.NOT.STATIC) THEN
          NSTEP=NSTEP+1
C$OMP     PARALLEL DO
C$OMP&    PRIVATE(nq),
C$OMP&    SHARED(niqV,nr,nx,YQ,YQS)
          DO nq=NQR(1,nr),NQR(2,nr)
            YQS(Vm,nq)=YQ(nq,niqV,1,nx)
          ENDDO
C$OMP     END PARALLEL DO
        ENDIF !not static

        !Transmembrane potential equations
        IF(TRANSMEMBRANE) THEN
          IF(STATIC) THEN
            ERROR_FLAG=.FALSE.
C$OMP       PARALLEL DO
C$OMP&      PRIVATE(CMAMDT,NEWDVM,nq,PSTIMULUS,STIMULUS),
C$OMP&      SHARED(DT,maqp1t0,maqp1t1,maqp1i,maqp2t0,maqp2t1,maqp2i,
C$OMP&        niqoldDV,niqPSTIM,niqV,NQGP_PIVOT,nr,nx_ext,nx,AQ,CQ,
C$OMP&        NQGW,T,YQ)
            DO nq=NQR(1,nr),NQR(2,nr)
              IF(.NOT.ERROR_FLAG) THEN
                !Calculate pseudo stimulus
                CALL CALC_STIMULUS2(maqp1t0,maqp1t1,maqp1i,maqp2t0,
     &            maqp2t1,maqp2i,niqV,NQGP(0,nq),NQGP_PIVOT(1,nq),
     &            nx_ext,AQ(1,nq),NQGW(1,nq),
     &            PSTIMULUS,STIMULUS,T,YQ,.FALSE.,
     &            ERROR,*100)

                !Update YQ's from new pseudo stimulus
                CMAMDT=DT/(CQ(1,nq)*CQ(2,nq))
                NEWDVM=YQ(nq,niqoldDV,1,nx)-
     &            (YQ(nq,niqPSTIM,1,nx)*CMAMDT)+
     &            (PSTIMULUS*CMAMDT)

                YQ(nq,niqPSTIM,1,nx)=PSTIMULUS
                YQ(nq,niqoldDV,1,nx)=NEWDVM
                YQ(nq,niqV,1,nx)=YQ(nq,niqV,1,nx)+NEWDVM

                GOTO 102
 100            CONTINUE
C$OMP           CRITICAL(MARCH8_COUPLED_1)
                ERROR_FLAG=.TRUE.
                WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
                CALL WRITES(IOER,OP_STRING,ERROR,*101)
                WRITE(OP_STRING,'(/'' >>An error occurred - '
     &            //'results may be unreliable'')')
                CALL WRITES(IODI,OP_STRING,ERROR,*101)
 101            CONTINUE
C$OMP           END CRITICAL(MARCH8_COUPLED_1)
 102            CONTINUE
              ENDIF !.NOT.ERROR_FLAG
            ENDDO !nq
C$OMP       END PARALLEL DO
          ELSE
            ERROR_FLAG=.FALSE.
C$OMP       PARALLEL DO
C$OMP&      PRIVATE(nq,LOCAL_ICQS,LOCAL_RCQS,STIMULUS,VARIANT),
C$OMP&      SHARED(INTEGRATOR_IWORK,maqdt,maqp1t0,maqp1t1,maqp1i,
C$OMP&        maqp2t0,
C$OMP&        maqp2t1,maqp2i,niqPSTIM,niqoldDV,niqV,NQGP,NQGP_PIVOT,nr,
C$OMP&        nx_ext,nx,INTEGRATOR_WORK,AQ,CQ,NQGW,SIZES,T,YQ,YQS,
C$OMP&        ICQS,RCQS)
            DO nq=NQR(1,nr),NQR(2,nr)
              IF(.NOT.ERROR_FLAG) THEN

                !Commands only used in parallel
C$              LOCAL_ICQS(1)=1 !variant always equals 1?
C$              DO nmq=2,NQIT
C$                LOCAL_ICQS(nmq)=ICQS(nmq)
C$              ENDDO
C$              DO nmq=1,NQRT
C$                LOCAL_RCQS(nmq)=RCQS(nmq)
C$              ENDDO

                IF(CELL_SPATIALLY_VARYING)
     &            CALL CELL_ASSIGN_SPATIAL(CELL_ICQS_VALUE,
     &              LOCAL_ICQS,ICQS_SPATIAL,
     &              IICQS_SPATIAL,IRCQS_SPATIAL,nq,CELL_RCQS_VALUE,
     &              LOCAL_RCQS,
     &              RCQS_SPATIAL,ERROR,*200)

                VARIANT = LOCAL_ICQS(CELL_VARIANT_OFFSET)
                !Fill in the sizes array
                SIZES(NUM_EQN)=CELL_NUM_STATE(VARIANT)
                SIZES(NUM_CONTROL)=CELL_NUM_CONTROL(VARIANT)
                SIZES(NUM_MODEL)=CELL_NUM_MODEL(VARIANT)
                SIZES(NUM_DERIVED)=CELL_NUM_DERIVED(VARIANT)
                SIZES(NUM_PARAM)=CELL_NUM_PARAMETERS(VARIANT)
                SIZES(NUM_PROTOCOL)=CELL_NUM_PROTOCOL(VARIANT)
                SIZES(NUM_AII)=CELL_NUM_AII(VARIANT)
                SIZES(NUM_AIO)=CELL_NUM_AIO(VARIANT)
                SIZES(NUM_ARI)=CELL_NUM_ARI(VARIANT)
                SIZES(NUM_ARO)=CELL_NUM_ARO(VARIANT)
                SIZES(NUM_TIME)=2
                !Calculate stimulus and pseudo stimulus
                CALL CALC_STIMULUS2(maqp1t0,maqp1t1,maqp1i,maqp2t0,
     &            maqp2t1,maqp2i,niqV,NQGP(0,nq),NQGP_PIVOT(1,nq),
     &            nx_ext,AQ(1,nq),NQGW(1,nq),
     &            LOCAL_RCQS(CELL_PROTOCOL_OFFSET(VARIANT)),STIMULUS,T,
     &            YQ,.TRUE.,ERROR,*200)
                !Calculate ionic currents
                CALL INTEGRATOR(INTEGRATOR_IWORK(1,nq),
     &            LOCAL_ICQS(CELL_AII_OFFSET(VARIANT)),
     &            LOCAL_ICQS(CELL_AIO_OFFSET(VARIANT)),
     &            LOCAL_ICQS(CELL_CONTROL_OFFSET(VARIANT)),maqdt,
     &            LOCAL_ICQS(CELL_MODEL_OFFSET(VARIANT)),
     &            CELL_NUM_ODE(VARIANT),CELL_NUM_STATE(VARIANT),nr,nx,
     &            SIZES,LOCAL_ICQS(CELL_VARIANT_OFFSET),
     &            INTEGRATOR_WORK(1,nq),AQ(1,nq),
     &            LOCAL_RCQS(CELL_ARI_OFFSET(VARIANT)),
     &            LOCAL_RCQS(CELL_ARO_OFFSET(VARIANT)),
     &            YQS(CELL_DERIVED_OFFSET(VARIANT),nq),
     &            LOCAL_RCQS(CELL_PARAMETERS_OFFSET(VARIANT)),
     &            LOCAL_RCQS(CELL_PROTOCOL_OFFSET(VARIANT)),STIMULUS,T,
     &            YQS(CELL_STATE_OFFSET(VARIANT),nq),RHSROUTINE,ERROR,
     &            *200)

                YQ(nq,niqPSTIM,1,nx)=
     &            LOCAL_RCQS(CELL_PROTOCOL_OFFSET(VARIANT))
                YQ(nq,niqoldDV,1,nx)=YQS(CELL_STATE_OFFSET(VARIANT),nq)-
     &            YQ(nq,niqV,1,nx)
                YQ(nq,niqV,1,nx)=YQS(CELL_STATE_OFFSET(VARIANT),nq)

                GOTO 202
 200            CONTINUE
C$OMP           CRITICAL(MARCH8_COUPLED_2)
                ERROR_FLAG=.TRUE.
                WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
                CALL WRITES(IOER,OP_STRING,ERROR,*201)
                WRITE(OP_STRING,'(/'' >>An error occurred - '
     &            //'results may be unreliable'')')
                CALL WRITES(IODI,OP_STRING,ERROR,*201)
 201            CONTINUE
C$OMP           END CRITICAL(MARCH8_COUPLED_2)
 202            CONTINUE
              ENDIF !.NOT.ERROR_FLAG
            ENDDO !nq
C$OMP       END PARALLEL DO
          ENDIF !static

          !Calculate transmembrane RHS vector
          CALL CALC_TRANSMEMBRANE_RHS(NENQ,niqV,NQS,NQXI,nr,NWQ,nx,
     &      nx_ext,NXQ,AQ,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,RHS,YQ,ERROR,*9999)

          !Solve transmembrane potential equations
          CALL SOLVE_SYSTEM(ISC_GKK(1,nx),ISR_GKK(1,nx),NQT,NQT,
     &      NQT,NZZT(1,NRLIST(1),nx),IWRIT4(NRLIST(1),nx),
     &      PRECON_CODE(nx),SOLVEROPTION(nx),SPARSEGKK(nx),GKK(1,nx),
     &      RHS,YQ(1,niqV,1,nx),FIRST_A,UPDATE_MATRIX,X_INIT,nx,
     &      ERROR,*9999)
          UPDATE_MATRIX=.FALSE.

          !DVm is needed with modified VCD
          IF(ITYP3(nr,nx).EQ.3.AND.KTYP33.EQ.2) THEN !VCDC
C$OMP       PARALLEL DO
C$OMP&      PRIVATE(nq),
C$OMP&      SHARED(DT,niq_old,niqV,nx,YQ,YQS)
            DO nq=1,NQT
              YQS(CELL_STATE_OFFSET(1)+DVm-1,nq)=(YQ(nq,niqV,1,nx)
     &          -YQ(nq,niq_old,1,nx))/DT
            ENDDO
C$OMP       END PARALLEL DO
          ENDIF

          IF(CALC_ACTIV_TIMES) THEN
            CALL MAQ_LOC(MAQ_INQUIRE,MAQ_TIME,atime,MAQ_ACTIV_TIME,
     &        ERROR,*9999)
            CALL MAQ_LOC(MAQ_INQUIRE,MAQ_TIME,dpot,MAQ_M_DPOT,
     &        ERROR,*9999)
            DO nq=NQR(1,nr),NQR(2,nr)
              IF((YQ(nq,niqV,1,nx)-YQ(nq,niq_OLD,1,nx)).GT.
     &          AQ(dpot,nq)) THEN
                AQ(dpot,nq)=YQ(nq,niqV,1,nx)-YQ(nq,niq_OLD,1,nx)
                AQ(atime,nq)=T
              ENDIF
            ENDDO
          ENDIF
        ENDIF

        !Extracellular potential equations
        IF(EXTRACELLULAR) THEN
          !Grid region is coupled to BE regions
          IF(COUPLEDGRIDBEM) THEN
            IF(CPTYPE.EQ.1) THEN
              !Define the required classes
              DO nxx=1,9
                NX_TORSO_LIST(nxx)=0
              ENDDO
              DO nxx=3,NXLIST(0)
                nxc=NXLIST(nxx)
                CALL NX_LOC(NX_INQUIRE,nxc,NX_TORSO_LIST(nxx-1),
     &            NX_SOLVE,
     &          ERROR,*9999)
              ENDDO !nxx

              !Generate RHS vector from YP
              CALL GEN_GRID_POTE_RHS2(IBT,IDO,INP,NBH,NBJ,NEELEM,NENQ,
     &          NHE,NHP,niqV,NKHE,NKH,NPF,NPNE,NPNODE,NQGP,NQGP_PIVOT,
     &          NQS,NQXI,nr,NRLIST,NVHE,NVHP,NW,NWQ(1,0,1),
     &          nx_ext,nx,NX_TORSO_LIST,NXQ(-NIM,0,0,1),NYNE,
     &          NYNP,ALPHA,AQ,CQ,CURVCORRECT,DNUDXQ,DXDXIQ,DXDXIQ2,
     &          NQGW,RHS,SE,XIQ,YP,YQ,ZA,ZE,ZP,FIXQ,ERROR,*9999)

            ELSE IF((CPTYPE.EQ.2).OR.(CPTYPE.EQ.4)) THEN
              !Generate RHS vector for extracellular potential solution
              NRTEMP(0)=1
              NRTEMP(1)=NRLIST(1)
              CALL GEN_EXT_RHS(niqV,niqBNDRY,NQGP,NQGP_PIVOT,NQGW,GM,
     &          NRTEMP,NTIME_INTERP,NTIME_POINTS,NWQ(1,0,1),nx_ext,nx,
     &          PROPQ,RHS,T,TIME_VALUES,YQ,CQ,FIXQ(1,1,nx_ext),
     &          ERROR,*9999)
              DO nqq=1,NQLIST(0)
                nq=NQLIST(nqq)
                RHS(nq)=0.0d0
              ENDDO !nq

            ELSE IF((CPTYPE.EQ.3).OR.(CPTYPE.EQ.5)) THEN
              !Define the required classes
              DO nxx=1,9
                NX_TORSO_LIST(nxx)=0
              ENDDO
              nxc=NXLIST(NXLIST(0))
              CALL NX_LOC(NX_INQUIRE,nxc,NX_TORSO_LIST(2),NX_SOLVE,
     &          ERROR,*9999)
              NRTEMP(0)=2
              NRTEMP(1)=NRLIST(1)
              NRTEMP(2)=NRLIST(NRLIST(0))

              !Generate RHS vector from YP
              CALL GEN_GRID_POTE_RHS2(IBT,IDO,INP,NBH,NBJ,NEELEM,NENQ,
     &          NHE,NHP,niqV,NKHE,NKH,NPF,NPNE,NPNODE,NQGP,NQGP_PIVOT,
     &          NQS,NQXI,nr,NRTEMP,NVHE,NVHP,NW,NWQ(1,0,1),
     &          nx_ext,nx,NX_TORSO_LIST,NXQ(-NIM,0,0,1),NYNE,
     &          NYNP,ALPHA,AQ,CQ,CURVCORRECT,DNUDXQ,DXDXIQ,DXDXIQ2,
     &          NQGW,RHS,SE,XIQ,YP,YQ,ZA,ZE,ZP,FIXQ,ERROR,*9999)
            ELSE IF(CPTYPE.EQ.6) THEN
              !Generate RHS vector for extracellular potential solution
              NRTEMP(0)=1
              NRTEMP(1)=NRLIST(1)
              CALL GEN_EXT_RHS(niqV,niqBNDRY,NQGP,NQGP_PIVOT,NQGW,GM,
     &          NRTEMP,NTIME_INTERP,NTIME_POINTS,NWQ(1,0,1),nx_ext,nx,
     &          PROPQ,RHS,T,TIME_VALUES,YQ,CQ,FIXQ(1,1,nx_ext),
     &          ERROR,*9999)
C$OMP         PARALLEL DO
C$OMP&        PRIVATE(nq),
C$OMP&        SHARED(NSIZE,NWQ,RHS)
              DO nq=1,NSIZE
                IF(NWQ(1,nq,1).GT.0) RHS(nq)=0.0d0
                IF(nq.GT.NQT) RHS(nq)=0.0d0
              ENDDO
C$OMP         END PARALLEL DO
            ELSE IF(CPTYPE.EQ.7) THEN
              !Generate RHS vector for extracellular potential solution
              NRTEMP(0)=1
              NRTEMP(1)=NRLIST(1)
              CALL GEN_EXT_RHS(niqV,niqBNDRY,NQGP,NQGP_PIVOT,NQGW,GM,
     &          NRTEMP,NTIME_INTERP,NTIME_POINTS,NWQ(1,0,1),nx_ext,nx,
     &          PROPQ,RHS,T,TIME_VALUES,YQ,CQ,FIXQ(1,1,nx_ext),
     &          ERROR,*9999)
C$OMP         PARALLEL DO
C$OMP&        PRIVATE(nq),
C$OMP&        SHARED(NSIZE,NWQ,RHS)
              DO nq=1,NSIZE
                IF(NWQ(1,nq,1).GT.0) RHS(nq)=0.0d0
                IF(nq.GT.NQT) RHS(nq)=0.0d0
              ENDDO
C$OMP         END PARALLEL DO
            ENDIF
          ELSE
            !Generate RHS vector for extracellular potential solution
            NRTEMP(0)=1
            NRTEMP(1)=NRLIST(1)
            CALL GEN_EXT_RHS(niqV,niqBNDRY,NQGP,NQGP_PIVOT,NQGW,GM,
     &        NRTEMP,
     &        NTIME_INTERP,NTIME_POINTS,NWQ(1,0,1),nx_ext,nx,PROPQ,RHS,
     &        T,TIME_VALUES,YQ,CQ,FIXQ(1,1,nx_ext),ERROR,*9999)
          ENDIF

          IF(CPTYPE.GE.6) THEN
            CALL SOLVE_SYSTEM(ISC_GKK(1,nx_ext),ISR_GKK(1,nx_ext),
     &        NSIZE,NSIZE,NSIZE,NZZT(1,NRLIST(1),nx_ext),
     &        IWRIT4(NRLIST(1),nx_ext),PRECON_CODE(nx_ext),
     &        SOLVEROPTION(nx_ext),SPARSEGKK(nx_ext),GKK(1,nx_ext),
     &        RHS,YQ(1,niqV,1,nx_ext),FIRST_EXT_A,UPDATE_EXT_MATRIX,
     &        X_INIT,nx_ext,ERROR,*9999)
            UPDATE_EXT_MATRIX=.FALSE.

            IF(CPTYPE.EQ.6) THEN
              nv=1
              nh=NH_LOC(1,nx3)
              nc=1
              SUM1=0.0d0
              SUM2=0.0d0
              DO nq=1,NPNODE(0,nr3)
                np=NPNODE(nq,nr3)
                IF(np.EQ.SOL_ACT_FIX_NODE) no2=nq
              ENDDO

              IF(USESALU) THEN
                DO nq=1,NSIZE+1
                  IF(nq.EQ.NQNP(SOL_ACT_FIX_NODE)) THEN
                    !do nothing
                  ELSE IF(nq.GT.NQNP(SOL_ACT_FIX_NODE)) THEN
                    SUM1=SUM1+(YQ(nq,1,1,nx2)*YQ(nq-1,niqV,1,nx_ext))
                    SUM2=SUM2+(YQ(nq,1,1,nx2)*YQ(nq-1,niqv,1,nx3))
                  ELSE
                    SUM1=SUM1+(YQ(nq,1,1,nx2)*YQ(nq,niqV,1,nx_ext))
                    SUM2=SUM2+(YQ(nq,1,1,nx2)*YQ(nq,niqv,1,nx3))
                  ENDIF
                ENDDO
                ALPHA=-SUM1/SUM2
CMLB please leave
C                ALPHA=SUM1/(1.0d0-SUM2)
C                WRITE(OP_STRING,'('' ALPHA = '',D12.6)') ALPHA
C                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

                DO nq=1,NSIZE
                  YQ(nq,niqV,1,nx_ext)=YQ(nq,niqV,1,nx_ext)+
     &              (ALPHA*YQ(nq,niqV,1,nx3))
                ENDDO
              ENDIF

              DO nq=NSIZE+1,NQNP(SOL_ACT_FIX_NODE)+1,-1
                YQ(nq,niqV,1,nx_ext)=YQ(nq-1,niqV,1,nx_ext)
              ENDDO
              YQ(NQNP(SOL_ACT_FIX_NODE),niqV,1,nx_ext)=0.0d0

              !write node solution back into YP
              DO nq=1,NPNODE(0,nr3)
                np=NPNODE(nq,nr3)
                ny=NYNP(1,nv,nh,np,0,nc,nr3)
                YP(ny,1,nx3)=YQ(NQNP(np),niqV,1,nx_ext)
              ENDDO

CMLB please leave
C              !Write flux back into YP
C              nc=2
C              DO nq=1,NPNODE(0,nr3)
C                np=NPNODE(nq,nr3)
C                ny=NYNP(1,nv,nh,np,0,nc,nr3)
C                IF(NQNP(np).GT.NQT) THEN
C                  YP(ny,1,nx3)=0.0d0
C                ELSE
C                  CALL GGRADPHIQDN(NENQ,niqV,NQNP(np),NQS,NQXI,
C     '              NXQ(-NIM,0,0,1),AQ,CQ(6,NQNP(np)),DNUDXQ,DXDXIQ,
C     '              DXDXIQ2,GDPHIDN,YQ(1,1,1,nx_ext),ERROR,*9999)
C                  YP(ny,1,nx3)=-GDPHIDN
C                ENDIF
C              ENDDO
C              nc=1

              IF(DOP) THEN
                !Check out the error in overwritten equation
                SUM1=0.0d0
                DO nq=1,NPNODE(0,nr3)
                  np=NPNODE(nq,nr3)
                  ny=NYNP(1,nv,nh,np,0,nc,nr3)
                  SUM1=SUM1+GK2(((nq-1)*NPNODE(0,nr3))+
     &              no2)*YP(ny,1,nx3)
                ENDDO
                WRITE(OP_STRING,'('' New Error :'',F12.6)') SUM1
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF

            ELSE IF(CPTYPE.EQ.7) THEN
              nv=1
              nh=NH_LOC(1,nx4)
              nc=1
              P1=0.0d0
              SUM1=0.0d0
              SUM2=0.0d0
              DO nq=1,NPNODE(0,nr4)
                np=NPNODE(nq,nr4)
                IF(np.EQ.SOL_ACT_FIX_NODE) no2=nq
              ENDDO

              IF(USESALU) THEN
                DO nq=1,NSIZE+1
                  IF(nq.EQ.NQNP(SOL_ACT_FIX_NODE)) THEN
                    !do nothing
                  ELSE IF(nq.GT.NQNP(SOL_ACT_FIX_NODE)) THEN
                    SUM1=SUM1+(YQ(nq,1,1,nx3)*YQ(nq-1,niqV,1,nx_ext))
                    SUM2=SUM2+(YQ(nq,1,1,nx3)*YQ(nq-1,niqv,1,nx4))
                  ELSE
                    SUM1=SUM1+(YQ(nq,1,1,nx3)*YQ(nq,niqV,1,nx_ext))
                    SUM2=SUM2+(YQ(nq,1,1,nx3)*YQ(nq,niqv,1,nx4))
                  ENDIF
                ENDDO

                DO nq=1,NPNODE(0,nr4)
                  np=NPNODE(nq,nr4)
                  IF(NQNP(np).GT.NQT) THEN
                    !do nothing
                  ELSE
                    CALL GGRADPHIQDN(NENQ,niqV,NQNP(np),NQS,NQXI,
     &                NXQ(-NIM,0,0,1),AQ,CQ(6,NQNP(np)),DNUDXQ,DXDXIQ,
     &                DXDXIQ2,GDPHIDN,YQ(1,1,1,nx_ext),ERROR,*9999)
                    P1=P1+GK3(((nq-1)*NPNODE(0,nr4))+no2)*GDPHIDN
                  ENDIF
                ENDDO

                ALPHA=(SUM1+P1)/(1.0d0-SUM2)

CMLB please leave
C                ALPHA=-SUM1/SUM2
C                ALPHA=(SUM1-P1)/(1.0d0-SUM2)
C                ALPHA=SUM1/(1.0d0-SUM2)
C                WRITE(OP_STRING,'('' ALPHA = '',D12.6)') ALPHA
C                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

                DO nq=1,NSIZE
                  YQ(nq,niqV,1,nx_ext)=YQ(nq,niqV,1,nx_ext)+
     &              (ALPHA*YQ(nq,niqV,1,nx4))
                ENDDO

C                DO nq=1,NSIZE+1
C                  IF(nq.EQ.NQNP(SOL_ACT_FIX_NODE)) THEN
C                    !do nothing
C                  ELSE IF(nq.GT.NQNP(SOL_ACT_FIX_NODE)) THEN
C                    SUM1=SUM1+(YQ(nq,1,1,nx3)*YQ(nq-1,niqV,1,nx_ext))
C                    SUM2=SUM2+(YQ(nq,1,1,nx3)*YQ(nq-1,niqv,1,nx4))
C                  ELSE
C                    SUM1=SUM1+(YQ(nq,1,1,nx3)*YQ(nq,niqV,1,nx_ext))
C                    SUM2=SUM2+(YQ(nq,1,1,nx3)*YQ(nq,niqv,1,nx4))
C                  ENDIF
C                ENDDO
C                ALPHA=-SUM1/SUM2
CCMLB please leave
CC                ALPHA=SUM1/(1.0d0-SUM2)
CC                WRITE(OP_STRING,'('' ALPHA = '',D12.6)') ALPHA
CC                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C
C                DO nq=1,NSIZE
C                  YQ(nq,niqV,1,nx_ext)=YQ(nq,niqV,1,nx_ext)+
C     '              (ALPHA*YQ(nq,niqV,1,nx4))
C                ENDDO
              ENDIF

              DO nq=NSIZE+1,NQNP(SOL_ACT_FIX_NODE)+1,-1
                YQ(nq,niqV,1,nx_ext)=YQ(nq-1,niqV,1,nx_ext)
              ENDDO
              YQ(NQNP(SOL_ACT_FIX_NODE),niqV,1,nx_ext)=0.0d0

              DO nq=1,NPNODE(0,nr4)
                np=NPNODE(nq,nr4)
                ny=NYNP(1,nv,nh,np,0,nc,nr4)
                YP(ny,1,nx4)=YQ(NQNP(np),niqV,1,nx_ext)
              ENDDO

CMLB please leave
C              !Write flux back into YP
C              nc=2
C              DO nq=1,NPNODE(0,nr4)
C                np=NPNODE(nq,nr4)
C                ny=NYNP(1,nv,nh,np,0,nc,nr4)
C                IF(NQNP(np).GT.NQT) THEN
C                  YP(ny,1,nx4)=0.0d0
C                ELSE
C                  CALL GGRADPHIQDN(NENQ,niqV,NQNP(np),NQS,NQXI,
C     '              NXQ(-NIM,0,0,1),AQ,CQ(6,NQNP(np)),DNUDXQ,DXDXIQ,
C     '              DXDXIQ2,GDPHIDN,YQ(1,1,1,nx_ext),ERROR,*9999)
C                  YP(ny,1,nx4)=-GDPHIDN
C                ENDIF
C              ENDDO
C              nc=1

              IF(DOP) THEN
                !Check out the error in overwritten equation
                SUM1=0.0d0
                DO nq=1,NPNODE(0,nr4)
                  np=NPNODE(nq,nr4)
                  ny=NYNP(1,nv,nh,np,0,nc,nr4)
                  SUM1=SUM1+GK3(((nq-1)*NPNODE(0,nr4))+
     &              no2)*YP(ny,1,nx4)
                ENDDO
                WRITE(OP_STRING,'('' New Error :'',F12.6)') SUM1
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF

            ENDIF
          ELSE
          !Solve extracellular potential equations
            CALL SOLVE_SYSTEM(ISC_GKK(1,nx_ext),ISR_GKK(1,nx_ext),
     &        NQT,NQT,NQT,NZZT(1,NRLIST(1),nx_ext),
     &        IWRIT4(NRLIST(1),nx_ext),PRECON_CODE(nx_ext),
     &        SOLVEROPTION(nx_ext),SPARSEGKK(nx_ext),GKK(1,nx_ext),RHS,
     &        YQ(1,niqV,1,nx_ext),FIRST_EXT_A,UPDATE_EXT_MATRIX,X_INIT,
     &        nx_ext,ERROR,*9999)
            UPDATE_EXT_MATRIX=.FALSE.
          ENDIF

          IF((CPTYPE.EQ.2).OR.(CPTYPE.EQ.3)) THEN
            nv=1
            nh=NH_LOC(1,nx2)
            nc=1

            !Write potential back into YP
            DO nq=1,NPNODE(0,nr2)
              np=NPNODE(nq,nr2)
              ny=NYNP(1,nv,nh,np,0,nc,nr2)
              YP(ny,1,nx2)=YQ(NQNP(np),niqV,1,nx_ext)
            ENDDO

            !Write flux back into YP
            nc=2
            DO nq=1,NPNODE(0,nr2)
              np=NPNODE(nq,nr2)
              ny=NYNP(1,nv,nh,np,0,nc,nr2)
              IF(NQNP(np).GT.NQT) THEN
                YP(ny,1,nx2)=0.0d0
              ELSE
                CALL GGRADPHIQDN(NENQ,niqV,NQNP(np),NQS,NQXI,
     &            NXQ(-NIM,0,0,1),AQ,CQ(6,NQNP(np)),DNUDXQ,DXDXIQ,
     &            DXDXIQ2,GDPHIDN,YQ(1,1,1,nx_ext),ERROR,*9999)
                YP(ny,1,nx2)=-GDPHIDN
              ENDIF
            ENDDO
            nc=1
          ENDIF

          IF((CPTYPE.EQ.4).OR.(CPTYPE.EQ.5)) THEN
            nv=1
            nh=NH_LOC(1,nx2)
            nc=1

            !Write potential back into YP (class 2)
            DO nq=1,NPNODE(0,nr2)
              np=NPNODE(nq,nr2)
              ny=NYNP(1,nv,nh,np,0,nc,nr2)
              YP(ny,1,nx2)=YQ(NQNP(np),niqV,1,nx_ext)
            ENDDO

            !Write potential back into YP (class 3)
            nh=NH_LOC(1,nx3)
            DO nq=1,NPNODE(0,nr3)
              np=NPNODE(nq,nr3)
              ny=NYNP(1,nv,nh,np,0,nc,nr3)
              YP(ny,1,nx3)=YQ(NQNP(np),niqV,1,nx_ext)
            ENDDO
          ENDIF

        ENDIF

        IF(STATIC) THEN
          CONTINU=.FALSE.
        ELSE
          T=T+DT

C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP     CRITICAL(MARCH8_COUPLED_3)
          IF(IWRIT5(NRLIST(1),nx).EQ.2) THEN
            IF(MOD(NSTEP,IWRIT1(NRLIST(1),nx)).EQ.0) THEN
              CALL CPU_TIMER(CPU_USER,TIME_STOP2)
              ELAPSED_TIME=TIME_STOP2(1)-TIME_START2(1)
              WRITE(OP_STRING,'(1X,''Time: '',E12.4,'', DT: '',E12.4,'
     &          //''' Iters: '',I8,'' Active: '',I6,'' Cpu: '',F6.2)')
     &          T,DT,NSTEP,NWQ(4,0,1),ELAPSED_TIME
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              CALL CPU_TIMER(CPU_USER,TIME_START2)
            ENDIF
          ENDIF
CC$OMP     END CRITICAL(MARCH8_COUPLED_3)

          IF(T.GE.TFINISH) CONTINU=.FALSE.
        ENDIF !static
      ENDDO !time

      IF(.NOT.STATIC) THEN
        T_RESTART(nx)=T

C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP   CRITICAL(MARCH8_COUPLED_4)
        IF(IWRIT5(NRLIST(1),nx).GE.1) THEN
          CALL CPU_TIMER(CPU_USER,TIME_STOP)
          ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)
          WRITE(OP_STRING,'(1X,I6,'' iterations : '',F8.2,'
     &      //'''s cpu'')') NSTEP,ELAPSED_TIME
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
CC$OMP   END CRITICAL(MARCH8_COUPLED_4)
      ENDIF

      !Overwrite the stored Vm solution with the current
      IF(UPDATE_VM) THEN
C$OMP   PARALLEL DO
C$OMP&  PRIVATE(nq),
C$OMP&  SHARED(niq_OLD,niqV,nr,nx,YQ)
        DO nq=NQR(1,nr),NQR(2,nr)
          YQ(nq,niq_OLD,1,nx)=YQ(nq,niqV,1,nx)
        ENDDO
C$OMP   END PARALLEL DO
      ENDIF

      CALL EXITS('MARCH8_COUPLED')
      RETURN
 9999 CALL ERRORS('MARCH8_COUPLED',ERROR)
      CALL EXITS('MARCH8_COUPLED')
      RETURN 1
      END
