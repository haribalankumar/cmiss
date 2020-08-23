      SUBROUTINE MARCH8(INTEGRATOR_IWORK,CELL_ICQS_VALUE,IBT,ICQS,
     &  ICQS_SPATIAL,IDO,IICQS_SPATIAL,IRCQS_SPATIAL,INP,ISC_GKK,
     &  ISR_GKK,NAQ,NBH,NBJ,NEELEM,NENQ,  NXI,NQNE,  NLATNE,NLATNQ,
     &  NLATPNQ,NQNLAT,NHE,NHP,NKHE,NKH,NLQ,NPF,NPNE,NPNODE,NQGP,
     &  NQGP_PIVOT,NQGW,NQS,NQSCNB,NQXI,NRLIST,NTIME_INTERP,
     &  NTIME_POINTS,NVHE,NVHP,NW,NWQ,NXLIST,NXQ,NYNE,NYNP,TV_BC_SET,
     &  INTEGRATOR_WORK,ALPHA,AQ,CELL_RCQS_VALUE,CQ,CURVCORRECT,DNUDXQ,
     &  DXDXIQ,DXDXIQ2,GCHQ,GKK,GM,GUQ,PG,PROPQ,RCQS,
     &  RCQS_SPATIAL,RHS,SE,TIME_VALUES,WG,XIQ,XQ,YP,YQ,YQS,ZA,ZE,ZP,
     &  FIXQ,ITER8,RHSROUTINE,RET_ERROR,*)

C#### Subroutine: MARCH8
C###  Description:
C###    MARCH8 is the time stepping routine for the general grid
C###    ionic current activation problem.
C***  Created by Martin Buist, June 1998

      IMPLICIT NONE

      INCLUDE 'adam00.cmn'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
C      INCLUDE 'cmiss$reference:cbfe01.cmn'
      INCLUDE 'cell02.cmn'
      INCLUDE 'cell_vcd.inc'
      INCLUDE 'cell_reserved.inc'
      INCLUDE 'error0.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'integrator.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'solver.inc'
      INCLUDE 'lsoda00.cmn'
      INCLUDE 'maqloc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'marc00.cmn'
      INCLUDE 'mxch.inc'
      INCLUDE 'nqloc00.inc'
      INCLUDE 'solv00.cmn'
      INCLUDE 'time02.cmn'
      INCLUDE 'time_variable.cmn'
      INCLUDE 'tol00.cmn'
      INCLUDE 'cellml.cmn'

!     Parameter list
      INTEGER INTEGRATOR_IWORK(INTEGRATOR_LIWORK,NQM),
     &  CELL_ICQS_VALUE(NQIM,NQVM),IBT(3,NIM,NBFM),ICQS(NQIM),
     &  ICQS_SPATIAL(NQISVM,NQM),IDO(NKM,NNM,0:NIM,NBFM),
     &  IICQS_SPATIAL(0:NQISVM,NQVM),IRCQS_SPATIAL(0:NQRSVM,NQVM),
     &  INP(NNM,NIM,NBFM),ISC_GKK(NISC_GKKM,NXM),ISR_GKK(NISR_GKKM,NXM),
     &  NAQ(NQM,NAM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     &  NEELEM(0:NE_R_M,0:NRM),NENQ(0:8,NQM),NLATNE(NEQM+1),
     &  NLATNQ(NEQM*NQEM),NLATPNQ(NQM),NQNLAT(NEQM*NQEM),NHE(NEM,NXM),
     &  NHP(NPM,0:NRM,NXM),NKHE(NKM,NNM,NHM,NEM),NKH(NHM,NPM,NCM,0:NRM),
     &  NLQ(NQM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     &  NQGP(0:NQGM,NQM),NQGP_PIVOT(NQGM,NQM),NQNE(NEQM,NQEM),NQS(NEQM),
     &  NQSCNB(NQSCM),NQXI(0:NIM,NQSCM),NRLIST(0:NRM),
     &  NTIME_INTERP(NTIMEVARSM),NTIME_POINTS(NTIMEVARSM),
     &  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     &  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     &  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),NW(NEM,3),NWQ(8,0:NQM,NAM),
     &  NXLIST(0:NXM),NXQ(-NIM:NIM,0:4,0:NQM,NAM),
     &  NXI(-NIM:NIM,0:NEIM,0:NEM),TV_BC_SET(0:NIQSM,0:NQM)
      REAL*8 INTEGRATOR_WORK(INTEGRATOR_LWORK,NQM),ALPHA,AQ(NMAQM,NQM),
     &  CELL_RCQS_VALUE(NQRM,NQVM),CQ(NMM,NQM),CURVCORRECT(2,2,NNM,NEM),
     &  dNudXQ(3,3,NQM),dXdXiQ(3,3,NQM),DXDXIQ2(3,3,NQM),GCHQ(3,NQM),
     &  GKK(NZ_GKK_M,NXM),GM(NZ_GM_M),GUQ(3,3,NQM),
     &  NQGW(NQGM,NQM),PG(NSM,NUM,NGM,NBM),PROPQ(3,3,4,2,NQM),
     &  RCQS(NQRM),RCQS_SPATIAL(NQRSVM,NQM),RHS(NQM),SE(NSM,NBFM,NEM),
     &  TIME_VALUES(2,0:NTIMEPOINTSM+1,NTIMEVARSM),WG(NGM,NBM),
     &  XIQ(NIM,NQM),XQ(NJM,NQM),YP(NYM,NIYM,NXM),YQ(NYQM,NIQM,NAM,NXM),
     &  YQS(NIQSM,NQM),ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),
     &  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER RET_ERROR*(*)
      LOGICAL FIXQ(NYQM,NIYFIXM,NXM),ITER8
C MHT removed 21-11-00
C     CHARACTER STRING*(MXCH),TIME_VARIABLE_NAMES(NTIMEVARSM)*(*)
!     Local variables
      INTEGER atime,dpot,ERR,i,IFROMC,maqdt,maqp1t0,maqp1t1,maqp1i,
     &  maqp2t0,maqp2t1,maqp2i,na,niqBNDRY,niq_old,niqV,NITB,nnq,
     &  NNQ_MIN,nonrlist,nq,nr,nr_torso,nr_blood,nr_grid,NRTEMP(0:1),
     &  NSOL,NSTEP,nxc,nx_ext,nx,nx_torso,nx_blood,nx_upd,SIZES(11),
     &  niqs,nzero,time_variable,VARIANT
      REAL*8 CHANGE,COEFFSEXT(NQGM),DT1MTHETA,DY(500),OLDDT,PHI1,PHI2,
     &  STIMULUS,T,TIME(2),EVTIME_FCN
      REAL ELAPSED_TIME,TIME_START(1),TIME_START2(1),TIME_STOP(1),
     &  TIME_STOP2(1)
      LOGICAL ADAPTIVE,BIDOMAIN,CHMTRIX,CONTINU,FIRSTITER,FIRST_A,
     &  IMPLICIT,UPDATE_MATRIX,UP_GRID_DELTAT,X_INIT,ERROR_FLAG
      CHARACTER CFROMI*2,CHAR2*2
      CHARACTER ERROR*(ERRSTRLEN)

!     Functions
      INTEGER LEN_TRIM

C      INTEGER IJ,mq,ni,nj,
C      REAL*8 DENOM,DET,dXidX(3,3),dXiN(3),NUMER,VALUE,
C     INTEGER L_AII,L_AIO,L_ARI,L_ARO,L_CONTROL,L_MODEL,L_PARAMETERS
C     PARAMETER(L_AII=1,L_AIO=1,L_ARI=1,L_ARO=1,L_CONTROL=1,
C    '  L_MODEL=2,L_PARAMETERS=99)

C *** DPN 22 June 1999 - need local arrays to store parameters to avoid
C     multiprocessing problems
      INTEGER nmq
      INTEGER LOCAL_ICQS(NQIM)
      REAL*8  LOCAL_RCQS(NQRM)

      EXTERNAL RHSROUTINE

      SAVE FIRSTITER,NNQ_MIN,NSTEP,OLDDT
C    &  SIZES,TIME_START

      CALL ENTERS('MARCH8',*9999)

      DT=TINCR
      X_INIT=.FALSE.
      CHMTRIX=.FALSE.
      IF(NMGT.GT.1) THEN !more than one grid level
        ADAPTIVE=.TRUE.
      ELSE               !fine grid level only
        ADAPTIVE=.FALSE.
      ENDIF

! Get indices for AQ auxiliary parameters
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE1,maqp1t0,MAQ_START,
     &  ERROR,*9999)
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE1,maqp1t1,MAQ_STOP,
     &  ERROR,*9999)
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE1,maqp1i,MAQ_CURRENT,
     &  ERROR,*9999)
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE2,maqp2t0,MAQ_START,
     &  ERROR,*9999)
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE2,maqp2t1,MAQ_STOP,
     &  ERROR,*9999)
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE2,maqp2i,MAQ_CURRENT,
     &  ERROR,*9999)
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_TIME,maqdt,MAQ_DT,
     &  ERROR,*9999)
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_TIME,atime,MAQ_ACTIV_TIME,
     &  ERROR,*9999)
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_TIME,dpot,MAQ_M_DPOT,
     &  ERROR,*9999)

! Get class information
      nxc=NXLIST(1)
      CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
      CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     &  ERROR,*9999)

! Check for Bidomain solution
      IF(ITYP19(NRLIST(1),nx).EQ.1.OR.ITYP19(NRLIST(1),nx).EQ.7) THEN
!electrical or coupled
        IF(KTYP32.EQ.1) THEN !monodomain
          BIDOMAIN=.FALSE.
        ELSE IF(KTYP32.EQ.2) THEN !bidomain
          BIDOMAIN=.TRUE.
        ELSE
          CALL ASSERT(.FALSE.,
     &      '>>Invalid type in KTYP32 (mono/bidomain)',ERROR,*9999)
        ENDIF
      ELSE
        BIDOMAIN=.FALSE.
      ENDIF

      IF(BIDOMAIN) THEN
        CALL ASSERT(NXLIST(0).GE.2,'>>Error - bidomain needs 2 classes',
     &    ERROR,*9999)
        nxc=NXLIST(2)
        CALL NX_LOC(NX_INQUIRE,nxc,nx_ext,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx_ext.GT.0,'>>No nx defined for this solve class',
     &    ERROR,*9999)
      ELSE
        nx_ext=0
      ENDIF
      IF(NXLIST(0).GE.3) THEN
        nxc=NXLIST(3)
        CALL NX_LOC(NX_INQUIRE,nxc,nx_upd,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx_upd.GT.0,'>>No nx defined for this solve class',
     &    ERROR,*9999)
      ELSE
        nx_upd=0
      ENDIF

! Check for implicit/explicit solution method
C SGM 22 November 2000 set IMPLICIT to true so ASSEMBLE10_FE will be
C     called, even for explicit cases
C MLT 29Nov02 Added grid FV
      IF(THETA(1).GT.ZERO_TOL.OR.(ITYP4(NRLIST(1),nx).EQ.6.OR.
     &                            ITYP4(NRLIST(1),nx).EQ.7)) THEN
        IMPLICIT=.TRUE.
      ELSE
        IMPLICIT=.FALSE.
      ENDIF

      DT1MTHETA = DT*(1.0d0-THETA(1))

      IF(RESTART) THEN
        T=T_RESTART(nx)
      ELSE
        IF(.NOT.ITER8) THEN
          T=TSTART
          NSTEP=0
          FIRSTITER=.TRUE.
          NNQ_MIN=1

C ***     Perform any intialisation required by the integrators
          IF(KTYP37.EQ.5) THEN ! Adams-Moulton integrator
            DO nq=1,NQT
              INTEGRATOR_IWORK(1,nq)=0
              INTEGRATOR_WORK(1,nq)=DT
            ENDDO
          ELSEIF(KTYP37.EQ.6) THEN ! LSODA integrator
C$OMP       PARALLEL DO
C$OMP&        DEFAULT(none)
C$OMP&        PRIVATE(i,nq),
C$OMP&        SHARED(DT,INTEGRATOR_IWORK,INTEGRATOR_WORK,
C$OMP&          LSODA_MAX_ITERS,LSODA_MAX_STEP,NQT)
            DO nq=1,NQT
              DO i=1,20 ! Only need to initialise optional inputs
                INTEGRATOR_IWORK(i+LSODA_NUM_INT_COMMON,nq)=0
              ENDDO
              DO i=1,22 ! Only need to initialise optional inputs
                INTEGRATOR_WORK(i+LSODA_NUM_REAL_COMMON,nq)=0.d0
              ENDDO
              INTEGRATOR_WORK(5+LSODA_NUM_REAL_COMMON,nq)= DT
              INTEGRATOR_WORK(6+LSODA_NUM_REAL_COMMON,nq)=
     &          LSODA_MAX_STEP
              INTEGRATOR_WORK(7+LSODA_NUM_REAL_COMMON,nq)=0.d0
              INTEGRATOR_IWORK(LSODA_ISTATE,nq)=1
              INTEGRATOR_IWORK(6+LSODA_NUM_INT_COMMON,nq)=
     &          LSODA_MAX_ITERS
            ENDDO
C$OMP       END PARALLEL DO
          ENDIF

          IF(ITYP19(NRLIST(1),nx).EQ.1.OR.ITYP19(NRLIST(1),nx).EQ.7)
     &      THEN
            DO nq=1,NQT
              AQ(atime,nq)=0.0d0
              AQ(dpot,nq)=0.0d0
            ENDDO
          ENDIF

! Initialise dynamic tracking
          IF(KTYP36.EQ.1) THEN !using DTAR
C$OMP       PARALLEL DO
C$OMP&        DEFAULT(none)
C$OMP&        PRIVATE(nq)
C$OMP&        SHARED(NWQ,NQT)
            DO nq=1,NQT
              NWQ(5,nq,1)=0
              NWQ(4,nq,1)=0
            ENDDO
C$OMP       END PARALLEL DO
            NWQ(4,0,1)=0
            NWQ(5,0,1)=0
          ELSE                 !not using DTAR
C$OMP       PARALLEL DO
C$OMP&        DEFAULT(none)
C$OMP&        PRIVATE(nq)
C$OMP&        SHARED(NWQ,NQT)
            DO nq=1,NQT
              NWQ(5,nq,1)=nq
              NWQ(4,nq,1)=2
            ENDDO
C$OMP       END PARALLEL DO
            NWQ(4,0,1)=NQT
            NWQ(5,0,1)=NQT
          ENDIF !DTAR
        ENDIF !.NOT.ITER8
      ENDIF !RESTART

C!!! Shouldn't this only be IF(.NOT.RESTART)
      OLDDT=DT
      UP_GRID_DELTAT=DT.NE.OLDDT

! Create/update matrices if implicit
      IF(IMPLICIT.OR.BIDOMAIN) THEN
        IF(RESTART) THEN
          UPDATE_MATRIX=.FALSE.
          FIRST_A=.FALSE.
          IF(UP_GRID_MATERIAL.OR.UP_GRID_TENSOR.OR.UP_GRID_DELTAT) THEN
C SGM 13 Nov 2000 call ASSEMBLE10_FE for grid-based FE
C MLT 29Nov02 call ASSEMBLE10_FE for grid FV also
            IF(ITYP4(NRLIST(1),nx).EQ.4) THEN !collocation
              CALL ASSEMBLE10(ISC_GKK,ISR_GKK,NEELEM,NENQ,NLATNE,NLATNQ,
     &          NLATPNQ,NLQ,NQGP,NQGP_PIVOT,NQNLAT,NQS,NQXI,NRLIST,
     &          NWQ(1,0,1),nx_ext,nx,nx_upd,NXQ,AQ,CQ,DNUDXQ,DXDXIQ,
     &          DXDXIQ2,GCHQ,GUQ,GKK,NQGW,PROPQ,XQ,BIDOMAIN,.FALSE.,
     &          FIRST_A,FIXQ,IMPLICIT,.TRUE.,UP_GRID_DELTAT,
     &          UPDATE_MATRIX,.FALSE.,ERROR,*9999)
            ELSEIF(ITYP4(NRLIST(1),nx).EQ.6.OR.
     &             ITYP4(NRLIST(1),nx).EQ.7) THEN !grid-based FE/grid FV
              CALL ASSEMBLE10_FE(ISC_GKK,ISR_GKK,NEELEM,NLATNE,
     &          NQGP,NQGP_PIVOT,NQNLAT,NQS,NQSCNB,NQXI,NRLIST,
     &          nx_ext,nx,nx_upd,NXQ,CQ,GKK,GM,NQGW,PG,WG,XQ,PROPQ,
     &          BIDOMAIN,.FALSE.,FIRST_A,FIXQ,IMPLICIT,.TRUE.,
     &          UP_GRID_DELTAT,UPDATE_MATRIX,.FALSE.,ERROR,*9999)
            ENDIF
            CHMTRIX=.TRUE.
            FIRST_A=.FALSE.
          ENDIF
        ELSE !not restart
          IF(.NOT.ITER8) THEN
C SGM 13 Nov 2000 call ASSEMBLE10_FE for grid-based FE
C MLT 29Nov02 call ASSEMBLE10_FE for grid FV also
            IF(ITYP4(NRLIST(1),nx).EQ.4) THEN !collocation
              CALL ASSEMBLE10(ISC_GKK,ISR_GKK,NEELEM,NENQ,NLATNE,NLATNQ,
     &          NLATPNQ,NLQ,NQGP,NQGP_PIVOT,NQNLAT,NQS,NQXI,NRLIST,
     &          NWQ(1,0,1),nx_ext,nx,nx_upd,NXQ,AQ,CQ,DNUDXQ,DXDXIQ,
     &          DXDXIQ2,GCHQ,GUQ,GKK,NQGW,PROPQ,XQ,BIDOMAIN,.FALSE.,
     &          FIRST_A,FIXQ,IMPLICIT,.FALSE.,.FALSE.,UPDATE_MATRIX,
     &          .FALSE.,ERROR,*9999)
            ELSEIF(ITYP4(NRLIST(1),nx).EQ.6.OR.
     &             ITYP4(NRLIST(1),nx).EQ.7) THEN !grid-based FE/grid FV
              CALL ASSEMBLE10_FE(ISC_GKK,ISR_GKK,NEELEM,NLATNE,
     &          NQGP,NQGP_PIVOT,NQNLAT,NQS,NQSCNB,NQXI,NRLIST,
     &          nx_ext,nx,nx_upd,NXQ,CQ,GKK,GM,NQGW,PG,WG,XQ,PROPQ,
     &          BIDOMAIN,.FALSE.,FIRST_A,FIXQ,IMPLICIT,.FALSE.,.FALSE.,
     &          UPDATE_MATRIX,.FALSE.,ERROR,*9999)
            ENDIF
            CHMTRIX=.TRUE.
          ENDIF
        ENDIF !restart

      ELSE !explicit
        IF(.NOT.RESTART) THEN
C         Calculate NQGP and NQGW
          CALL GET_FD_POINTS(NEELEM,NENQ,NLQ,NQGP,NQGP_PIVOT,NQS,
     &      NQXI,NRLIST,NWQ(1,0,1),NXQ,ERROR,*9999)
          DO nonrlist=1,NRLIST(0)
            nr=NRLIST(nonrlist)
            NITB=NQXI(0,NQS(NEELEM(1,nr)))
            ERROR_FLAG=.FALSE.
C$OMP       PARALLEL DO
CC$OMP&        DEFAULT(none)
C$OMP&        PRIVATE(COEFFSEXT,ERROR,na,nq,OP_STRING)
C$OMP&        SHARED(ADAPTIVE,CQ,DNUDXQ,DOP,DXDXIQ,DXDXIQ2,ERROR_FLAG,
C$OMP&        FIXQ,GCHQ,GUQ,IMPLICIT,IODI,ITYP4,NENQ,NIM,NITB,NLQ,NMGT,
C$OMP&        NQGP,NQGW,NQR,NQS,NQXI,nr,NWQ,nx,NXQ,nx_ext,PROPQ)
            DO nq=NQR(1,nr),NQR(2,nr)
              IF(.NOT.ERROR_FLAG
     &          .AND.(.NOT.ADAPTIVE
     &          .OR.ADAPTIVE.AND.NLQ(nq).GT.0
     &          .AND.NLQ(nq).LE.NMGT)) THEN !need residual coeffs
                IF(ADAPTIVE) THEN !adaptive grid
                  na=NLQ(nq)
                  IF(DOP) THEN
                    WRITE(OP_STRING,'('' nq='',I6,'' na='',I2)') nq,na
                    CALL WRITES(IODI,OP_STRING,ERROR,*110)
                  ENDIF
                ELSE !not adaptive
                  na=1
                ENDIF !adaptive
C SGM 18-10-00 added check for Grid-based Finite Element 
C MLT 10Dec02 added check for grid FV also
                IF(ITYP4(nr,nx).NE.6.AND.
     &             ITYP4(nr,nx).NE.7) THEN ! Not Grid-based FE
                  CALL CALC_GRID_COEF(NENQ,NITB,nq,NQGP,NQS,NQXI,
     &              NWQ(1,nq,na),NXQ(-NIM,0,0,na),nx_ext,nx,COEFFSEXT,
     &              CQ,DNUDXQ,DXDXIQ,DXDXIQ2,GCHQ(1,nq),GUQ(1,1,nq),
     &              NQGW(1,nq),PROPQ(1,1,1,1,nq),.FALSE.,FIXQ,IMPLICIT,
     &              .FALSE.,ERROR,*110)
                ENDIF
                GOTO 112
 110              CONTINUE
                  ERROR_FLAG=.TRUE.
                  IF(ERROR.NE.' ') THEN
                    CALL FLAG_ERROR(0,ERROR(:LEN_TRIM(ERROR)))
                  ENDIF
 112            CONTINUE
              ENDIF !.NOT.ERROR_FLAG
            ENDDO !nq
C$OMP       END PARALLEL DO
C KAT 2001-04-11: Aborting on error
            IF(ERROR_FLAG) THEN
              ERROR=' '
              GOTO 9999
            ENDIF
          ENDDO !nr
          UP_GRID_MATERIAL=.FALSE.
          UP_GRID_TENSOR=.FALSE.

        ELSE !restart
          IF(UP_GRID_MATERIAL.OR.UP_GRID_TENSOR) THEN
            DO nonrlist=1,NRLIST(0)
              nr=NRLIST(nonrlist)
              NITB=NQXI(0,NQS(NEELEM(1,nr)))
              ERROR_FLAG=.FALSE.

C$OMP         PARALLEL DO 
CC$OMP&          DEFAULT(none)
C$OMP&          PRIVATE(ERROR,na,nq,COEFFSEXT)
C$OMP&          SHARED(ADAPTIVE,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,ERROR_FLAG,
C$OMP&            FIXQ,GCHQ,GUQ,IMPLICIT,NENQ,NIM,NITB,NLQ,NQGP,NQGW,
C$OMP&            NQR,NQS,NQXI,nr,NWQ,nx,NXQ,nx_ext,PROPQ)
              DO nq=NQR(1,nr),NQR(2,nr)
                IF(.NOT.ERROR_FLAG) THEN
                  IF(ADAPTIVE) THEN !adaptive grid
                    na=NLQ(nq)
                  ELSE !not adaptive
                    na=1
                  ENDIF
C KAT 23Dec00: .not.IMPLICIT implies ITYP4(nr,nx).NE.6
C SGM 18-10-00 added check for Grid-based Finite Element
C                  IF(ITYP4(nr,nx).NE.6) THEN ! Not Grid-based FE
                  CALL CALC_GRID_COEF(NENQ,NITB,nq,NQGP,NQS,NQXI,
     &              NWQ(1,nq,na),NXQ(-NIM,0,0,na),nx_ext,nx,COEFFSEXT,
     &              CQ,DNUDXQ,DXDXIQ,DXDXIQ2,GCHQ(1,nq),GUQ(1,1,nq),
     &              NQGW(1,nq),PROPQ(1,1,1,1,nq),.FALSE.,FIXQ,IMPLICIT,
     &              .FALSE.,ERROR,*120)
C                  ENDIF
                  GOTO 122
 120                CONTINUE
                    ERROR_FLAG=.TRUE.
                    IF(ERROR.NE.' ') THEN
                      CALL FLAG_ERROR(0,ERROR(:LEN_TRIM(ERROR)))
                    ENDIF
 122              CONTINUE
                ENDIF !.NOT.ERROR_FLAG
              ENDDO !nq
C$OMP         END PARALLEL DO
C KAT 2001-04-11: Aborting on error
              IF(ERROR_FLAG) THEN
                ERROR=' '
                GOTO 9999
              ENDIF
            ENDDO !nr
            UP_GRID_MATERIAL=.FALSE.
            UP_GRID_TENSOR=.FALSE.
          ENDIF !upgrid
        ENDIF !restart
      ENDIF !implicit/explicit

      CONTINU=.TRUE.
      CALL CPU_TIMER(CPU_USER,TIME_START)
      CALL CPU_TIMER(CPU_USER,TIME_START2)
      CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niq_old,NIQ_OLDSOLN,ERROR,*9999)
      CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqV,NIQ_V,ERROR,*9999)
C DPN 27 July 2001 - need to check for the case when a mechanics model
C                    is used and NIQ_LOC isn't set-up
      IF(niqV.EQ.0) THEN
        ! Assume that Vm is stored in the first index of YQ
        niqV = 1
      ENDIF
      CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqBNDRY,NIQ_BNDRY,ERROR,*9999)


C!!! CPB DO NOT COPY ICQS AND RCQS EACH TIME STEP.

CC *** DPN 22 June 1999 - Initialise the local copies of ICQS and RCQS.
CC     This should be the only time the arrays are copied for non-multi
CC     processor
CC ??? probably only works for simple (non-ipcell) models now ????
      IF(.NOT.CELL_SPATIALLY_VARYING) THEN
        LOCAL_ICQS(1)=1!variant always equals 1?
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
      ELSE
        DO nmq=1,NQIM
          LOCAL_ICQS(nmq)=0
        ENDDO
        DO nmq=1,NQRM
          LOCAL_RCQS(nmq)=0.0d0
        ENDDO
      ENDIF

C--------------------------- mono & bidomain ---------------------------
      IF(.NOT.ITER8) THEN !mono & bidomain problem
!     Main time integration loop
        DO WHILE(CONTINU)
          NSTEP=NSTEP+1
          NSOL=0

C!!! CPB UPDATE IN PARALLEL
          IF(ITYP19(NRLIST(1),nx).EQ.1.OR.ITYP19(NRLIST(1),nx).EQ.7)
     &      THEN !electrical or coupled model
C$OMP       PARALLEL DO 
C$OMP&        DEFAULT(none)
C$OMP&        PRIVATE(nq)
C$OMP&        SHARED(niq_old,niqV,NQT,nx,YQS,YQ)
            DO nq=1,NQT
C SGM 18Dec2000 store YQS old in YQ old temporarily
              ! Vm is always first
              YQ(nq,niq_old,1,nx) = YQS(Vm,nq)
              YQS(Vm,nq)=YQ(nq,niqV,1,nx)
            ENDDO
C$OMP       END PARALLEL DO
          ENDIF

          UP_GRID_DELTAT=DT.NE.OLDDT

          IF((IMPLICIT.OR.BIDOMAIN).AND.UP_GRID_DELTAT) THEN
C SGM 13 Nov 2000 call ASSEMBLE10_FE for grid-based FE
C MLT 29Nov02 call ASSEMBLE10_FE for grid FV also
            IF(ITYP4(NRLIST(1),nx).EQ.4) THEN !collocation
              CALL ASSEMBLE10(ISC_GKK,ISR_GKK,NEELEM,NENQ,NLATNE,NLATNQ,
     &          NLATPNQ,NLQ,NQGP,NQGP_PIVOT,NQNLAT,NQS,NQXI,NRLIST,
     &          NWQ(1,0,1),nx_ext,nx,nx_upd,NXQ(-NIM,0,0,1),AQ,CQ,
     &          DNUDXQ,DXDXIQ,DXDXIQ2,GCHQ,GUQ,GKK,NQGW,PROPQ,XQ,
     &          BIDOMAIN,.FALSE.,FIRST_A,FIXQ,IMPLICIT,.TRUE.,
     &          UP_GRID_DELTAT,UPDATE_MATRIX,.FALSE.,ERROR,*9999)
            ELSEIF(ITYP4(NRLIST(1),nx).EQ.6.OR.
     &             ITYP4(NRLIST(1),nx).EQ.7) THEN !grid-based FE/grid FV
              CALL ASSEMBLE10_FE(ISC_GKK,ISR_GKK,NEELEM,NLATNE,NQGP,
     &          NQGP_PIVOT,NQNLAT,NQS,NQSCNB,NQXI,NRLIST,nx_ext,nx,
     &          nx_upd,NXQ(-NIM,0,0,1),CQ,GKK,GM,NQGW,PG,WG,XQ,PROPQ,
     &          BIDOMAIN,.FALSE.,FIRST_A,FIXQ,IMPLICIT,.TRUE.,
     &          UP_GRID_DELTAT,UPDATE_MATRIX,.FALSE.,ERROR,*9999)
            ENDIF
            UPDATE_MATRIX=.TRUE.
            CHMTRIX=.TRUE.
          ENDIF

          IF(KTYP36.EQ.1) THEN !DTAR
            !add points
            CALL CALC_DTAR(maqp1t0,maqp1t1,maqp1i,maqp2t0,maqp2t1,
     &        maqp2i,niq_old,niqV,nnq_min,NQXI,NRLIST,NSOL,NWQ(1,0,1),
     &        nx,NXQ(-NIM,0,0,1),AQ,CQ,T,YQ(1,1,1,nx),.TRUE.,
     &        ERROR,*9999)
          ENDIF !DTAR

!       Loop over regions
          DO nonrlist=1,NRLIST(0)
            nr=NRLIST(nonrlist)
            ERR=0
            ERROR_FLAG=.FALSE.

            DO na=max(NMGT-1,1),1,-1 !loop over grid levels
              IF(DOP) THEN
                WRITE(OP_STRING,'(/'' Grid level na='',I2)') na
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF !DOP

!           Loop over grid points

C MLT 12March03 this list has not been alphabetised yet
C$OMP         PARALLEL DO 
CC$OMP&          DEFAULT(none)
C$OMP&          PRIVATE(DY,ERR,ERROR,LOCAL_ICQS,niqs,nmq,nnq,nq,
C$OMP&            LOCAL_RCQS,SIZES,STIMULUS,time_variable,VARIANT)
C$OMP&          SHARED(BIDOMAIN,CELL_ICQS_VALUE,CELL_NUM_ODE,
C$OMP&            CELL_NUM_STATE,CELL_PROTOCOL_OFFSET,CELL_RCQS_VALUE,
C$OMP&            CELL_SPATIALLY_VARYING,DOP,ICQS,ICQS_SPATIAL,
C$OMP&            IICQS_SPATIAL,IRCQS_SPATIAL,RCQS_SPATIAL,
C$OMP&            IMPLICIT,ITYP4,ITYP19,maqdt,maqp1t0,maqp1t1,maqp1i,
C$OMP&            maqp2t0,maqp2t1,maqp2i,na,NENQ,niqV,NLQ,nnq_min,NQIT,
C$OMP&            nr,NRLIST,NQGP,NQGP_PIVOT,NQGW,NQRT,NQS,NQXI,NWQ,nx,
C$OMP&            nx_ext,NXQ,AQ,CQ,GCHQ,GUQ,PROPQ,RCQS,RHS,T,YQ,
C$OMP&            ERROR_FLAG,OP_STRING,IOER,IODI,CELL_AII_OFFSET,
C$OMP&            INTEGRATOR_IWORK,CELL_MODEL_OFFSET,CELL_AIO_OFFSET,
C$OMP&            CELL_CONTROL_OFFSET,TIME_VALUES,KTYP33,ITYP3,
C$OMP&            INTEGRATOR_WORK,TIME,CELLML_ROUTINES,DT,
C$OMP&            CELL_ARO_OFFSET,CELL_ARI_OFFSET,
C$OMP&            CELL_PARAMETERS_OFFSET,CELL_DERIVED_OFFSET,
C$OMP&            CELL_NUM_PROTOCOL,CELL_NUM_PARAMETERS,CELL_NUM_AIO,
C$OMP&            CELL_NUM_AII,CELL_NUM_CONTROL,CELL_NUM_DERIVED,
C$OMP&            CELL_NUM_MODEL,CELL_NUM_ARI,CELL_STATE_OFFSET,YQS,
C$OMP&            NTIME_POINTS,NTIME_INTERP,THETA,CELL_NUM_ARO,
C$OMP&            KTYP3_init,TV_BC_SET)
C$OMP&          REDUCTION(+:NSOL)

C *** Loop over grid points in the active list - NWQ(5,nq,na)
              DO nnq=nnq_min,NWQ(5,0,1)
                IF(.NOT.ERROR_FLAG) THEN
                  nq=NWQ(5,nnq,1)
                  NSOL=NSOL+1
                  IF(NLQ(nq).EQ.na) THEN !nq belongs to level na

C SGM 22 November 2000 no difference for internal/external with grid-based FE
C MLT 29Nov02 Added grid Finite Volume
                    IF(NWQ(1,nq,1).EQ.0.OR.(ITYP4(NRLIST(1),nx).EQ.6.OR.
     &                                      ITYP4(NRLIST(1),nx).EQ.7))
     &                THEN !internal or grid-based FE method or grid FV
                      IF(DOP) THEN
                        WRITE(OP_STRING,'('' nq='',I6,'' na='',I2)')
     &                    nq,na
                        CALL WRITES(IODI,OP_STRING,ERROR,*100)
                      ENDIF !DOP

C ***                 Handle spatially varying parameters

C ***                 First, if multi-processing, copy all parameter
C                     values to the local arrays

C ???                 DPN 22 June 1999 - is there a better way to do this ??
C ???                 Again, probably only works for simple models now ???
C$                    IF(.NOT.CELL_SPATIALLY_VARYING) THEN
C$                      LOCAL_ICQS(1)=1 !variant always equals 1?
C$                      DO nmq=2,NQIT
C$                        LOCAL_ICQS(nmq)=ICQS(nmq)
C$                      ENDDO
C$                      DO nmq=1,NQRT
C$                        LOCAL_RCQS(nmq)=RCQS(nmq)
C$                      ENDDO
C$                    ENDIF

C ***                 Now overwite any spatial variance, this means that
C                     you're probably solving a "real" model, including
C                     variants.
                      IF(CELL_SPATIALLY_VARYING)
     &                  CALL CELL_ASSIGN_SPATIAL(CELL_ICQS_VALUE,
     &                    LOCAL_ICQS,ICQS_SPATIAL,IICQS_SPATIAL,
     &                    IRCQS_SPATIAL,nq,CELL_RCQS_VALUE,LOCAL_RCQS,
     &                    RCQS_SPATIAL,ERROR,*100)

                      VARIANT = LOCAL_ICQS(CELL_VARIANT_OFFSET)
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

C ***                 Calc Istim from diffusive flux, SACs &
C                     applied current
                      IF(ITYP19(nr,nx).EQ.1.OR.ITYP19(nr,nx).EQ.7) THEN
!electrical or coupled model
                        CALL CALC_STIMULUS(maqp1t0,maqp1t1,maqp1i,
     &                    maqp2t0,maqp2t1,maqp2i,%VAL(0),niqV,nq,
     &                    NQGP(0,nq),NQGP_PIVOT(1,nq),nr,nx_ext,nx,
     &                    AQ(1,nq),CQ(1,nq),NQGW(1,nq),
     &                    LOCAL_RCQS(CELL_PROTOCOL_OFFSET(VARIANT)),
     &                    STIMULUS,T,THETA(1),YQ,BIDOMAIN,.TRUE.,ERROR,
     &                    *100)

!PJH 25Oct99 Temporary code for SAC expts
C                        IF(DABS(XQ(2,nq)-1.d0).LT.1.d-3) THEN !pickup y=1 pts
C                          write(*,'(I7,2E15.4)') nq,XQ(1,nq),
C     '                      CELL_RCQS_VALUE(CELL_PROTOCOL_OFFSET,
C     '                    ICQS_SPATIAL(1,nq))
C                        ENDIF
                      ENDIF

                      !apply b.c.'s
                      IF(KTYP3_init(nx).EQ.4.AND.TV_BC_SET(0,0).GT.0)
     &                  THEN
C                        !dynamically allocate memory for XE
C                        XE_PTR=0
C                        CALL ALLOCATE_MEMORY(NSM*NJM,1,DPTYPE,XE_PTR,
C     '                    MEM_INIT,ERROR,*100)
                        DO niqs=1,TV_BC_SET(0,0)
                          time_variable=TV_BC_SET(niqs,nq)
                          IF(time_variable.GT.0) THEN
C                            VALUE=GET_TV_VALUE_AT_TIME(IBT,IDO,INP,
C     '                        time_variable,T,%VAL(XE_PTR))
c                            CALL EVTIME(NTIME_INTERP,NTIME_POINTS,
c     '                        TIME_VARIABLE,T,TIME_VALUES,VALUE,STRING,
c     '                        TIME_VARIABLE_NAMES,.FALSE.,ERROR,*100)
C                            WRITE(*,*) time_variable,VALUE
c                            YQS(CELL_STATE_OFFSET+
c     '                        TV_BC_SET(niqs,0)-1,nq)=VALUE
                            YQS(CELL_STATE_OFFSET(VARIANT)+
     &                        TV_BC_SET(niqs,0)-1,nq)=EVTIME_FCN(
     &                        NTIME_INTERP,NTIME_POINTS,time_variable,
     &                        T,TIME_VALUES)
                          ENDIF
                        ENDDO !niqs
C                        !free dynamically allocated memory
C                        CALL FREE_MEMORY(XE_PTR,ERROR,*100)
                      ENDIF

!                     Integrate cell ODEs for curent grid pt
                      IF (CELL_NUM_ODE(VARIANT).GT.0) THEN
                        IF (CELL_NUM_ODE(VARIANT).EQ.
     &                    CELL_NUM_STATE(VARIANT)) THEN
                          !probably safe to assume that CONTROL(ODE)
                          !is not used
                        ELSE
                          CALL ASSERT(CELL_NUM_CONTROL(VARIANT).GT.0,
     &                      'Must have at least 1 CONTROL variable '
     &                      //'(ODE)',ERROR,*100)
                          LOCAL_ICQS(CELL_CONTROL_OFFSET(VARIANT)-1
     &                      +ODE)=1

                        ENDIF

                        IF (ITYP19(nr,nx).EQ.1.AND.ITYP3(nr,
     &                    nx).EQ.10.AND.KTYP33.EQ.6) THEN ! cellml
                          CALL INTEGRATOR(INTEGRATOR_IWORK(1,nq),
     &                      LOCAL_ICQS(CELL_AII_OFFSET(VARIANT)),
     &                      LOCAL_ICQS(CELL_AIO_OFFSET(VARIANT)),
     &                      LOCAL_ICQS(CELL_CONTROL_OFFSET(VARIANT)),
     &                      maqdt,
     &                      LOCAL_ICQS(CELL_MODEL_OFFSET(VARIANT)),
     &                      CELL_NUM_ODE(VARIANT),
     &                      CELL_NUM_STATE(VARIANT),nr,nx,SIZES,
     &                      LOCAL_ICQS(CELL_VARIANT_OFFSET),
     &                      INTEGRATOR_WORK(1,nq),AQ(1,nq),
     &                      LOCAL_RCQS(CELL_ARI_OFFSET(VARIANT)),
     &                      LOCAL_RCQS(CELL_ARO_OFFSET(VARIANT)),
     &                      YQS(CELL_DERIVED_OFFSET(VARIANT),nq),
     &                      LOCAL_RCQS(CELL_PARAMETERS_OFFSET(VARIANT)),
     &                      LOCAL_RCQS(CELL_PROTOCOL_OFFSET(VARIANT)),
     &                      STIMULUS,T,
     &                      YQS(CELL_STATE_OFFSET(VARIANT),nq),
     &                      %VAL(CELLML_ROUTINES(VARIANT)),ERROR,*100)
                        ELSE
                          CALL INTEGRATOR(INTEGRATOR_IWORK(1,nq),
     &                      LOCAL_ICQS(CELL_AII_OFFSET(VARIANT)),
     &                      LOCAL_ICQS(CELL_AIO_OFFSET(VARIANT)),
     &                      LOCAL_ICQS(CELL_CONTROL_OFFSET(VARIANT)),
     &                      maqdt,
     &                      LOCAL_ICQS(CELL_MODEL_OFFSET(VARIANT)),
     &                      CELL_NUM_ODE(VARIANT),
     &                      CELL_NUM_STATE(VARIANT),nr,nx,SIZES,
     &                      LOCAL_ICQS(CELL_VARIANT_OFFSET),
     &                      INTEGRATOR_WORK(1,nq),AQ(1,nq),
     &                      LOCAL_RCQS(CELL_ARI_OFFSET(VARIANT)),
     &                      LOCAL_RCQS(CELL_ARO_OFFSET(VARIANT)),
     &                      YQS(CELL_DERIVED_OFFSET(VARIANT),nq),
     &                      LOCAL_RCQS(CELL_PARAMETERS_OFFSET(VARIANT)),
     &                      LOCAL_RCQS(CELL_PROTOCOL_OFFSET(VARIANT)),
     &                      STIMULUS,T,
     &                      YQS(CELL_STATE_OFFSET(VARIANT),nq),
     &                      RHSROUTINE,ERROR,*100)
                        ENDIF

                      ENDIF

!                     Evaluate non-ODEs for curent grid pt - only works
!                     for real models!!
                      IF (CELL_NUM_ODE(VARIANT).LT.
     &                  CELL_NUM_STATE(VARIANT).AND.
     &                  CELL_SPATIALLY_VARYING) THEN

                        CALL ASSERT(CELL_NUM_CONTROL(VARIANT).GT.0,
     &                    'Must have at least 1 CONTROL variable '
     &                    //'(ODE)',ERROR,*100)

c                        CELL_ICQS_VALUE(CELL_CONTROL_OFFSET-1+ODE,
c     '                    ICQS_SPATIAL(1,nq))=0

                        LOCAL_ICQS(CELL_CONTROL_OFFSET(VARIANT)-1+ODE)=0

C *** DPN 05 July 2000 - Now passing more time information into the
C ***   cellular models.
                        TIME(TCell) = T
                        TIME(DTCell) = DT

                        IF (ITYP19(nr,nx).EQ.1.AND.ITYP3(nr,
     &                    nx).EQ.10.AND.KTYP33.EQ.6) THEN ! cellml
                          ! do nothing
                        ELSE
                          CALL RHSROUTINE(TIME,
     &                      YQS(CELL_STATE_OFFSET(VARIANT),nq),DY,
     &                      LOCAL_ICQS(CELL_CONTROL_OFFSET(VARIANT)),
     &                      LOCAL_ICQS(CELL_MODEL_OFFSET(VARIANT)),
     &                      SIZES,LOCAL_ICQS(CELL_VARIANT_OFFSET),
     &                      YQS(CELL_DERIVED_OFFSET(VARIANT),nq),
     &                      LOCAL_RCQS(CELL_PARAMETERS_OFFSET(VARIANT)),
     &                      LOCAL_RCQS(CELL_PROTOCOL_OFFSET(VARIANT)),
     &                      LOCAL_ICQS(CELL_AII_OFFSET(VARIANT)),
     &                      LOCAL_ICQS(CELL_AIO_OFFSET(VARIANT)),
     &                      LOCAL_RCQS(CELL_ARI_OFFSET(VARIANT)),
     &                      LOCAL_RCQS(CELL_ARO_OFFSET(VARIANT)),ERR)

                          IF(ERR.NE.0) THEN
                            ERROR='Error in RHS Routine'
                            GOTO 100
                          ENDIF
                        ENDIF

                      ENDIF

                    ENDIF !internal

                    GOTO 102
 100                  CONTINUE
                      ERROR_FLAG=.TRUE.
                      IF(ERROR.NE.' ') THEN
                        CALL FLAG_ERROR(0,ERROR(:LEN_TRIM(ERROR)))
                      ENDIF
 102                CONTINUE
                  ENDIF !not adaptive or adaptive and NLQ(nq)=na
                ENDIF !.NOT.ERROR_FLAG
              ENDDO !grid points
C$OMP         END PARALLEL DO

C *** DPN 03 July 2000 - need to break out of here if an error has been
C ***   reported by the cell model.
C ??? Do we need to look for errors that can be ignored for the next
C ??? time step, or always break out for all errors ???
              CALL ASSERT(.NOT.ERROR_FLAG,
     &          'Error reported from cell model',ERROR,*9999)

              NITB=NQXI(0,NQS(NEELEM(1,nr)))

C MLT 8March03 call secondary flux calculation for Grid FV
              IF((ITYP4(NRLIST(1),nx).EQ.7)) THEN !Grid FV
                CALL CALC_FV_GRID_SECFLUX(NITB,nr,NXQ(-NIM,0,0,1),
     &            YQ(1,niqV,na,nx),YQ(1,niqV,na,nx_ext),PROPQ,BIDOMAIN,
     &            ERROR,*9999)
              ENDIF

C MLT 12March03 this list has not been alphabetised yet
C !!! KAT 28May03: DY,ERR,TIME,SIZES,STIMULUS set in each loop.
C$OMP         PARALLEL DO
CC$OMP&          DEFAULT(none)
C$OMP&          PRIVATE(nq,nzero,i,DY,ERR,ERROR,PHI1,PHI2,
C$OMP&            OP_STRING,LOCAL_ICQS,LOCAL_RCQS,VARIANT,niqs,
C$OMP&            SIZES,STIMULUS,TIME,time_variable)
C$OMP&          SHARED(NQR,nr,RHS,NLQ,na,ITYP4,NRLIST,ERROR_FLAG,
C$OMP&            nx,ITYP19,DT1MTHETA,BIDOMAIN,NQGP,NQGW,NQGP_PIVOT,
C$OMP&            YQ,niqV,GM,NQGM,YQS,CELL_STATE_OFFSET,
C$OMP&            nx_ext,PROPQ,ZERO_TOL,NWQ,IMPLICIT,niq_old,KTYP38,
C$OMP&            NENQ,NQS,NQXI,NXQ,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,DOP,
C$OMP&            IODI,CELL_ICQS_VALUE,ICQS_SPATIAL,IICQS_SPATIAL,
C$OMP&            IRCQS_SPATIAL,
C$OMP&            CELL_RCQS_VALUE,RCQS_SPATIAL,TV_BC_SET,KTYP3_init,
C$OMP&            NTIME_INTERP,NTIME_POINTS,T,TIME_VALUES,
C$OMP&            CELL_NUM_STATE,
C$OMP&            CELL_NUM_CONTROL,KTYP33,INTEGRATOR_IWORK,
C$OMP&            CELL_NUM_ODE,
C$OMP&            AQ,CELL_ARI_OFFSET,CELL_ARO_OFFSET,
C$OMP&            CELL_DERIVED_OFFSET,CELL_PARAMETERS_OFFSET,
C$OMP&            CELL_PROTOCOL_OFFSET,CELLML_ROUTINES,
C$OMP&            CELL_AII_OFFSET,CELL_AIO_OFFSET,CELL_CONTROL_OFFSET,
C$OMP&            CELL_NUM_AIO,CELL_NUM_ARI,CELL_NUM_ARO,CELL_NUM_AII,
C$OMP&            CELL_NUM_DERIVED,CELL_NUM_PARAMETERS,
C$OMP&            CELL_NUM_PROTOCOL,
C$OMP&            ITYP3,CELL_MODEL_OFFSET,CELL_NUM_MODEL,
C$OMP&            NIM,CELL_SPATIALLY_VARYING,maqdt,INTEGRATOR_WORK,
C$OMP&            DT,NITB,
C$OMP&            GKK,ISC_GKK,ISR_GKK,USE_LAT,SPARSEGKK)      
              DO nq=NQR(1,nr),NQR(2,nr)

C MLT 11March03 new initialisation of RHS for better parallelisation
                RHS(nq) = 0.0d0

                IF(.NOT.ERROR_FLAG) THEN
                  IF(NLQ(nq).EQ.na) THEN !nq belongs to level na

C SGM 22 November 2000 no difference for internal/external with grid-based FE
C MLT 29Nov02 Added Finite Volume
                    IF(ITYP4(NRLIST(1),nx).EQ.6.OR.
     &                 ITYP4(NRLIST(1),nx).EQ.7) THEN !Grid-based FE/grid FV
                      IF(ITYP19(nr,nx).EQ.1.OR.ITYP19(nr,nx).EQ.7)
     &                  THEN !electrical or coupled model

C loop over supporting grid points
                          DO nzero=1,NQGP(0,nq)

C MLT 11March03 new RHS construction for better parallelisation  
C MLT 11March03 removed reference to CELL_STATE_OFFSET(VARIANT) after
C talking to DPN
C Revision 1.1115 did the DT and NQGM multiplications and BIDOMAIN
C conditional 1/NQGP(0,nq) as many times.
                             RHS(nq) = RHS(nq) + 
     &                         NQGW(NQGP_PIVOT(nzero,nq),nq) *
     &                         DT1MTHETA * 
     &                         YQ(NQGP(nzero,nq),niqV,na,nx) +
     &                         GM(NQGP_PIVOT(nzero,nq)+NQGM*(nq-1))* 
     &                         YQS(Vm,NQGP(nzero,nq)) 

C add extracellular potential contribution for bidomain
                            IF (BIDOMAIN) THEN
                              RHS(nq) = RHS(nq) + 
     &                         NQGW(NQGP_PIVOT(nzero,nq),nq) *
     &                         DT * 
     &                         YQ(NQGP(nzero,nq),niqV,na,nx_ext)
                            ENDIF
                          ENDDO

C MLT 8March03 add in secondary diffusion contribution for grid FV
                        IF((ITYP4(NRLIST(1),nx).EQ.7)) THEN !Grid FV
                          DO i=1,NITB
                            RHS(nq) = RHS(nq) - DT*PROPQ(1,i,3,2,nq)- 
     &                                          DT*PROPQ(1,i,4,2,nq)
                          ENDDO
                        ENDIF ! Grid FV or not
                      ENDIF ! electrical or coupled model

                    ELSE !Not Grid-based FE
                      IF(NWQ(1,nq,1).GT.0) THEN
                        !external boundary grid point
                        IF(IMPLICIT) THEN
C KAT 23Dec2000 store YQ at old timestep
                          IF(ITYP19(nr,nx).EQ.1.OR.ITYP19(nr,nx).EQ.7)
     &                      THEN !electrical or coupled model
                            YQ(nq,niq_old,na,nx) = YQ(nq,niqV,na,nx)
                          ENDIF
                          IF(BIDOMAIN) THEN
                            IF(KTYP38.EQ.1) THEN !no Vm flux
                              RHS(nq)=0.0d0
                            ELSE IF(KTYP38.EQ.2) THEN !no PHIi flux
                              IF(USE_LAT.EQ.0) THEN
                                CALL GGRADPHIQDN(NENQ,niqV,nq,NQS,NQXI,
     &                            NXQ(-NIM,0,0,1),AQ,CQ(3,nq),DNUDXQ,
     &                            DXDXIQ,DXDXIQ2,RHS(nq),
     &                            YQ(1,1,1,nx_ext),ERROR,*200)
                              ELSE !for lattice method we look at the weights
                                   !in GKK and apply them to phi instead.
                                CALL ASSERT(SPARSEGKK(nx).EQ.1,
     &                            '>>Only implemented for SPARSENESS 1',
     &                            ERROR,*200)
                                DO i=ISR_GKK(nq,nx_ext),
     &                            (ISR_GKK(nq+1,nx_ext)-1)
                                  RHS(nq)=RHS(nq) + GKK(i,nx_ext)*
     &                              YQ(ISC_GKK(i,nx_ext),1,na,nx_ext)
                                ENDDO
                                RHS(nq)=-1.0d0*RHS(nq)
                              ENDIF
                              RHS(nq)=RHS(nq)/1000.0d0 !mS - S
C                              IF(nq.eq.20) THEN
C                                write(*,*) RHS(nq)
C                              ENDIF
                              IF(DABS(RHS(nq)).LT.1.0d-6) RHS(nq)=0.0d0
C                             IF(RHS(nq).LT.1.0d-6) RHS(nq)=0.0d0
                              !Needed to be DABS()!
                              !This is needed!
                              !Zero intracellular flux creates a feedback
                              !loop between Vm and Phie which can lead
                              !to an accumulation of numerical error,
                              !enough to corrupt solutions. This is
                              !especially a problem where meshes have
                              !lots of curvature.
                            ENDIF
                          ELSE !monodomain
                            RHS(nq)=0.0d0
                          ENDIF !bidomain/monodomain
                        ELSE !explicit
C KAT 26Oct00: UP_GRID_CONNECTIVITY is never set.
CC                        PJH 30Aug99
C                        IF(UP_GRID_CONNECTIVITY) THEN
CC                     Note: the zero flux bc is handled differently
CC                     if the fem update grid connectivity is used.
CC                     XQ(NJT+nj,nq) has coeffs of normal to cleavage plane.
C                          IF(DOP) THEN
C                            WRITE(OP_STRING,'('' nq='',I8)') nq
C                            CALL WRITES(IODI,OP_STRING,ERROR,*200)
C                          ENDIF !DOP
C                          CALL INVERT(NJT,dXdXiQ(1,1,nq),dXidX,DET)
C                          DO ni=1,NITB !calc dXidXj*unit_normal to bdry
C                            dXiN(ni)=0.d0
C                            DO nj=1,NJT
C                              dXiN(ni)=dXiN(ni)+dXidX(ni,nj)*
C     '                          XQ(nj+NJT,nq)
C                            ENDDO !nj
C                          ENDDO !ni
C                          NUMER=0.d0 !numerator   of expression obtained from GradPHI.n=0
C                          DENOM=0.d0 !denominator of expression obtained from GradPHI.n=0
C                          DO ni=1,NITB !loop over Xi direction
C                            DO IJ=-ni,ni,2*ni !..on either side of nq
C                              mq=NXQ(IJ,1,nq,na) !is neighbour in IJ direction
C                              IF(DOP) THEN
C                                WRITE(OP_STRING,'(''    ni='',I1,'
C     '                            //''' IJ='',I2,'' mq='',I8)') ni,IJ,mq
C                                CALL WRITES(IODI,OP_STRING,ERROR,*200)
C                              ENDIF !DOP
C                              IF(mq.GT.0) THEN !neighbour exists
C                                NUMER=NUMER+dXiN(ni)
C     '                            *YQS(CELL_STATE_OFFSET,mq)
C                                DENOM=DENOM+dXiN(ni)
C                                IF(DOP) THEN
C                                  WRITE(OP_STRING,'('' mq='',I8,'
C     '                              //''' NUMER='',E11.3,'' DENOM='','
C     '                              //'E11.3)') mq,NUMER,DENOM
C                                  CALL WRITES(IODI,OP_STRING,ERROR,
C     '                              *200)
C                                ENDIF !DOP
C                              ENDIF !mq>0
C                            ENDDO !IJ
C                          ENDDO !ni
C                          CALL ASSERT(DABS(DENOM).GT.ZERO_TOL,
C     '                      '>>DENOM=0 in MARCH8',ERROR,*200)
C                          YQ(nq,niqV,na,nx)=NUMER/DENOM
CC PJH 30Aug99
C                        ELSE !update grid connectivity not called
                          IF(ITYP19(nr,nx).EQ.1.OR.ITYP19(nr,nx).EQ.7)
     &                      THEN !electrical or coupled model
                            PHI1=YQS(Vm,NWQ(1,nq,1))
                            PHI2=YQS(Vm,NWQ(2,nq,1))
C                           Assumed no flux boundary condition
C SGM 21Dec2000 store YQ at old timestep
                            YQ(nq,niq_old,na,nx) = YQ(nq,niqV,na,nx)
                            YQ(nq,niqV,na,nx)=(4.0d0*PHI1-PHI2)/3.0d0
                          ENDIF
                          IF(DOP) THEN
                            WRITE(OP_STRING,'('' Boundary pt'//
     &                ' YQ('',I6,'','',I1,'','',I1,''): '',E12.3)')
     &                        nq,niqV,na,YQ(nq,niqV,na,nx)
                            CALL WRITES(IODI,OP_STRING,ERROR,*200)
                          ENDIF !DOP
                        ENDIF !implicit/explicit

C DPN 12 June 2000
C If solving a mechanics only cellular model, need to evaluate RHS
C for external points as well!!
C *** DPN 26 September 2000 - also for coupled models
                        IF((ITYP19(nr,nx).EQ.2).OR.(ITYP19(nr,nx).EQ.7))
     &                    THEN !mechanical or coupled model

                          IF(CELL_SPATIALLY_VARYING)
     &                      CALL CELL_ASSIGN_SPATIAL(CELL_ICQS_VALUE,
     &                        LOCAL_ICQS,ICQS_SPATIAL,IICQS_SPATIAL,
     &                        IRCQS_SPATIAL,nq,CELL_RCQS_VALUE,
     &                        LOCAL_RCQS,RCQS_SPATIAL,ERROR,*200)

                          VARIANT = LOCAL_ICQS(CELL_VARIANT_OFFSET)
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

                          !apply b.c.'s
                          IF(KTYP3_init(nx).EQ.4
     &                      .AND.TV_BC_SET(0,0).GT.0) THEN
                            DO niqs=1,TV_BC_SET(0,0)
                              time_variable=TV_BC_SET(niqs,nq)
                              IF(time_variable.GT.0) THEN
                                YQS(CELL_STATE_OFFSET(VARIANT)+
     &                            TV_BC_SET(niqs,0)-1,nq)=EVTIME_FCN(
     &                            NTIME_INTERP,NTIME_POINTS,
     &                            time_variable,T,TIME_VALUES)
                              ENDIF
                            ENDDO !niqs
                          ENDIF

!                         Integrate cell ODEs for curent grid pt
                          IF (CELL_NUM_ODE(VARIANT).GT.0) THEN
                            IF (CELL_NUM_ODE(VARIANT).EQ.
     &                        CELL_NUM_STATE(VARIANT)) THEN
                              !probably safe to assume that CONTROL(ODE)
                              !is not used
                            ELSE
                              CALL
     &                          ASSERT(CELL_NUM_CONTROL(VARIANT).GT.0,
     &                          'Must have at least 1 CONTROL variable '
     &                          //'(ODE)',ERROR,*200)

                              LOCAL_ICQS(CELL_CONTROL_OFFSET(VARIANT)-1
     &                          +ODE)=1
                            ENDIF
                            IF (ITYP19(nr,nx).EQ.1.AND.ITYP3(nr,
     &                        nx).EQ.10.AND.KTYP33.EQ.6) THEN !cellml
                              CALL INTEGRATOR(INTEGRATOR_IWORK(1,nq),
     &                          LOCAL_ICQS(CELL_AII_OFFSET(VARIANT)),
     &                          LOCAL_ICQS(CELL_AIO_OFFSET(VARIANT)),
     &                          LOCAL_ICQS(CELL_CONTROL_OFFSET(
     &                          VARIANT)),maqdt,
     &                          LOCAL_ICQS(CELL_MODEL_OFFSET(VARIANT)),
     &                          CELL_NUM_ODE(VARIANT),
     &                          CELL_NUM_STATE(VARIANT),nr,nx,SIZES,
     &                          LOCAL_ICQS(CELL_VARIANT_OFFSET),
     &                          INTEGRATOR_WORK(1,nq),AQ(1,nq),
     &                          LOCAL_RCQS(CELL_ARI_OFFSET(VARIANT)),
     &                          LOCAL_RCQS(CELL_ARO_OFFSET(VARIANT)),
     &                          YQS(CELL_DERIVED_OFFSET(VARIANT),nq),
     &                          LOCAL_RCQS(CELL_PARAMETERS_OFFSET(
     &                          VARIANT)),
     &                          LOCAL_RCQS(CELL_PROTOCOL_OFFSET(
     &                          VARIANT)),
     &                          %VAL(0),T,
     &                          YQS(CELL_STATE_OFFSET(VARIANT),nq),
     &                          %VAL(CELLML_ROUTINES(VARIANT)),ERROR,
     &                          *200)
                            ELSE
                              CALL INTEGRATOR(INTEGRATOR_IWORK(1,nq),
     &                          LOCAL_ICQS(CELL_AII_OFFSET(VARIANT)),
     &                          LOCAL_ICQS(CELL_AIO_OFFSET(VARIANT)),
     &                          LOCAL_ICQS(CELL_CONTROL_OFFSET(
     &                          VARIANT)),maqdt,
     &                          LOCAL_ICQS(CELL_MODEL_OFFSET(VARIANT)),
     &                          CELL_NUM_ODE(VARIANT),
     &                          CELL_NUM_STATE(VARIANT),nr,nx,SIZES,
     &                          LOCAL_ICQS(CELL_VARIANT_OFFSET),
     &                          INTEGRATOR_WORK(1,nq),AQ(1,nq),
     &                          LOCAL_RCQS(CELL_ARI_OFFSET(VARIANT)),
     &                          LOCAL_RCQS(CELL_ARO_OFFSET(VARIANT)),
     &                          YQS(CELL_DERIVED_OFFSET(VARIANT),nq),
     &                          LOCAL_RCQS(CELL_PARAMETERS_OFFSET(
     &                          VARIANT)),
     &                          LOCAL_RCQS(CELL_PROTOCOL_OFFSET(
     &                          VARIANT)),
     &                          %VAL(0),T,
     &                          YQS(CELL_STATE_OFFSET(VARIANT),nq),
     &                          RHSROUTINE,ERROR,*200)
                            ENDIF
                          ENDIF

!                         Evaluate non-ODEs for curent grid pt - only
!                         works for real models!!
                          IF (CELL_NUM_ODE(VARIANT).LT.
     &                      CELL_NUM_STATE(VARIANT).AND.
     &                      CELL_SPATIALLY_VARYING) THEN

                            CALL ASSERT(CELL_NUM_CONTROL(VARIANT).GT.0,
     &                        'Must have at least 1 CONTROL variable '
     &                        //'(ODE)',ERROR,*200)

C *** DPN 05 July 2000 - Now passing more time information into the
C ***   cellular models.
                            TIME(TCell) = T
                            TIME(DTCell) = DT

                            LOCAL_ICQS(CELL_CONTROL_OFFSET(VARIANT)-1
     &                        +ODE)=0
                            IF (ITYP19(nr,nx).EQ.1.AND.ITYP3(nr,
     &                        nx).EQ.10.AND.KTYP33.EQ.6) THEN !cellml
                              ! do nothing
                            ELSE
                              CALL RHSROUTINE(TIME,
     &                          YQS(CELL_STATE_OFFSET(VARIANT),nq),
     &                          DY,
     &                          LOCAL_ICQS(CELL_CONTROL_OFFSET(
     &                          VARIANT)),
     &                          LOCAL_ICQS(CELL_MODEL_OFFSET(VARIANT)),
     &                          SIZES,LOCAL_ICQS(CELL_VARIANT_OFFSET),
     &                          YQS(CELL_DERIVED_OFFSET(VARIANT),nq),
     &                          LOCAL_RCQS(CELL_PARAMETERS_OFFSET(
     &                          VARIANT)),
     &                          LOCAL_RCQS(CELL_PROTOCOL_OFFSET(
     &                          VARIANT)),
     &                          LOCAL_ICQS(CELL_AII_OFFSET(VARIANT)),
     &                          LOCAL_ICQS(CELL_AIO_OFFSET(VARIANT)),
     &                          LOCAL_RCQS(CELL_ARI_OFFSET(VARIANT)),
     &                          LOCAL_RCQS(CELL_ARO_OFFSET(VARIANT)),
     &                          ERR)

                              IF(ERR.NE.0) THEN
cC$OMP                         CRITICAL(MARCH8_ANDRE3)
                                ERROR='Error in RHS Routine'
cC$OMP                         END CRITICAL(MARCH8_ANDRE3)
                                GOTO 200
                              ENDIF
                            ENDIF

                          ENDIF

                        ENDIF !ITYP19(nr,nx).EQ.2 - mechanical model

                      ELSE !interior grid point
                        IF(ITYP19(nr,nx).EQ.1.OR.ITYP19(nr,nx).EQ.7)
     &                    THEN !electrical or coupled model
C SGM 21Dec2000 store YQ at old timestep
                          YQ(nq,niq_old,na,nx) = YQ(nq,niqV,na,nx)
                          YQ(nq,niqV,na,nx)=YQS(Vm,nq)
                          IF(IMPLICIT) RHS(nq)=YQS(Vm,nq)
                          IF(DOP) THEN
                            WRITE(OP_STRING,'('' Interior pt'
     &                //' YQ('',I6,'','',I1,'','',I1,''): '',E12.3)')
     &                        nq,niqV,na,YQ(nq,niqV,na,nx)
                            CALL WRITES(IODI,OP_STRING,ERROR,*200)
                          ENDIF !DOP
                        ENDIF !electrical
                      ENDIF !boundary/internal

                    ENDIF !grid-based FE/not

                  ENDIF !adaptive
                  GOTO 202
 200              CONTINUE
                  ERROR_FLAG=.TRUE.
                  IF(ERROR.NE.' ') THEN
                    CALL FLAG_ERROR(0,ERROR(:LEN_TRIM(ERROR)))
                  ENDIF
 202              CONTINUE
                ENDIF !.NOT.ERROR_FLAG
              ENDDO !grid points
C$OMP         END PARALLEL DO

C MLT 6March2003 placed the YQ update here so that it would not interfere
C with any attempt to use the values explicitly on the right hand side -
C note there is no reason not to use the same YQ update with collocation, 
C except that the example solutions might change slightly.
C SGM 18Dec2000 use YQS+(YQ-YQS_old) as initial guess of YQ at next timestep
C                        !YQS_old stored in YQ(,niq_old,,)
              IF((ITYP4(NRLIST(1),nx).EQ.6).OR.
     &           (ITYP4(NRLIST(1),nx).EQ.7)) THEN !Grid-based FE or Grid FV
C$OMP           PARALLEL DO
C$OMP&            DEFAULT(none)
C$OMP&            PRIVATE(CHANGE,nq)
C$OMP&            SHARED(na,niq_old,niqV,NQR,nr,nx,YQ,YQS)
                DO nq=NQR(1,nr),NQR(2,nr)
                  CHANGE=YQ(nq,niqV,na,nx)-YQ(nq,niq_old,na,nx)
C SGM 21Dec2000 store YQ at old timestep
                  YQ(nq,niq_old,na,nx) = YQ(nq,niqV,na,nx)
                  YQ(nq,niqV,na,nx)= YQS(Vm,nq) + CHANGE
                ENDDO !grid points
C$OMP           END PARALLEL DO
              ENDIF !Grid-based FE and FV

C KAT 2001-04-11: Aborting on error
              IF(ERROR_FLAG) THEN
                ERROR=' '
                GOTO 9999
              ENDIF

!           Interpolate some grid points if using adaptive grid
              IF(ADAPTIVE) CALL MG_INTERPOL(na,1,NAQ,NLQ,NXQ,
     &                                      YQ,.TRUE.,ERROR,*9999)
            ENDDO !grid level na

!         Copy from coarse grid level na=NLQ to fine level na=1
C$OMP       PARALLEL DO
C$OMP&        DEFAULT(none)
C$OMP&        PRIVATE(CHAR2,na,nq)
C$OMP&        SHARED(ITYP19,niqV,NLQ,NQT,nr,nx,YQ)
            DO nq=1,NQT
              IF(NLQ(nq).GT.0) THEN
                IF(NLQ(nq).LE.9) THEN
                  na=NLQ(nq)
                ELSE
                  CHAR2=CFROMI(NLQ(nq),'(I2)')
                  na=IFROMC(CHAR2(1:1))
                ENDIF
                IF(ITYP19(nr,nx).EQ.1.OR.ITYP19(nr,nx).EQ.7)
     &            YQ(nq,niqV,1,nx)=YQ(nq,niqV,na,nx)
                   !electrical or coupled model
              ENDIF !NLQ>0
            ENDDO !nq
C$OMP       END PARALLEL DO

          ENDDO !region

          IF(IMPLICIT) THEN
!           Solve the system for the transmembrane potential
            IF(SOLVEROPTION(nx).NE.SOLV_BMG .AND. 
     &          SOLVEROPTION(nx).NE.SOLV_PCG_BMG ) THEN
               ! Use Standard Solver library
               CALL SOLVE_SYSTEM(ISC_GKK(1,nx),ISR_GKK(1,nx),NQT,NQT,
     &              NQT,NZZT(1,NRLIST(1),nx),IWRIT4(NRLIST(1),nx),
     &              PRECON_CODE(nx),SOLVEROPTION(nx),SPARSEGKK(nx),
     &              GKK(1,nx),RHS,YQ(1,niqV,1,nx),FIRST_A,UPDATE_MATRIX,
     &              .FALSE.,nx,ERROR,*9999)
               UPDATE_MATRIX=.FALSE.
               !
             ELSE
                ! Use Boxmg package to solve for tensor product mesh           
                CALL BOXMG_SOLVE(
     &               ITYP4(NRLIST(1),nx),
     &               ISC_GKK(1,nx),ISR_GKK(1,nx),NQT,NQT,NQT,
     &               NZZT(1,NRLIST(1),nx),IWRIT4(NRLIST(1),nx),
     &               SOLVEROPTION(nx),SPARSEGKK(nx),GKK(1,nx),RHS,
     &               YQ(1,niqV,1,nx),FIRST_A,UPDATE_MATRIX,X_INIT,
     &               NEELEM,NQNE,NQS,NQXI,NXI,nr,nx,ERROR,*9999)
             ENDIF
             !
          ENDIF !implicit

          IF(BIDOMAIN) THEN
            IF(CHMTRIX) THEN
              IF(.NOT.RESTART) FIRST_A=.TRUE.
              UPDATE_MATRIX=.TRUE.
              CHMTRIX=.FALSE.
            ENDIF

C MLT 8March03 call secondary flux calculation for Grid FV -  this 
C second call uses the just computed Vm value to update the secondary
C flux. This updating is probably not necessary given the 
C timestep restrictions from elsewhere in the discretisation. The solution
C is not varying much from time to time and this call adds extra 
C computational work. The call is added here for completeness.
C            IF((ITYP4(NRLIST(1),nx).EQ.7)) THEN !Grid FV
C              CALL CALC_FV_GRID_SECFLUX(NITB,nr,
C     '          NXQ(-NIM,0,0,1),YQ(1,niqV,na,nx),
C     '          YQ(1,niqV,na,nx_ext),PROPQ,BIDOMAIN,ERROR,*9999)
C            ENDIF

!         Generate RHS vector for extracellular potential solution
            CALL GEN_EXT_RHS(niqV,niqBNDRY,NQGP,NQGP_PIVOT,NQGW,GM,
     &        NRLIST,NTIME_INTERP,NTIME_POINTS,NWQ(1,0,1),nx_ext,nx,
     &        PROPQ,RHS,T,TIME_VALUES,YQ,CQ,
     &        FIXQ(1,1,nx_ext),ERROR,*9999)

C KAT 31Mar01 use YQ +(YQ-YQ_old) as initial guess of YQ at next
C           timestep.  There is no reason why this can't be used for
C           collocation also except that example output may change
C           slightly.  See IPINI3 for initialization.
C MLT 11December02 added in grid FV
            IF((ITYP4(NRLIST(1),nx).EQ.6).OR.
     &         (ITYP4(NRLIST(1),nx).EQ.7)) THEN !Grid-based FE or Grid FV
C$OMP         PARALLEL DO
C$OMP&          DEFAULT(none)
C$OMP&          PRIVATE(CHANGE,nq)
C$OMP&          SHARED(na,niq_old,niqV,NQR,nr,nx,nx_ext,YQ,YQS)
              DO nq=NQR(1,nr),NQR(2,nr)
                CHANGE=YQ(nq,niqV,1,nx_ext)-YQ(nq,niq_old,1,nx_ext)
CC KAT 31Mar01 store YQ at old timestep
                YQ(nq,niq_old,1,nx_ext)=YQ(nq,niqV,1,nx_ext)
                YQ(nq,niqV,1,nx_ext)=YQ(nq,niqV,1,nx_ext)+CHANGE
              ENDDO !grid points
C$OMP         END PARALLEL DO
            ENDIF !Grid-based FE and FV

!           Solve system for the extracellular potential
            IF( SOLVEROPTION(nx_ext).NE.SOLV_BMG .AND. 
     &          SOLVEROPTION(nx_ext).NE.SOLV_PCG_BMG ) THEN
               ! Use Standard Solver library 
               CALL SOLVE_SYSTEM(ISC_GKK(1,nx_ext),ISR_GKK(1,nx_ext),
     &              NQT,NQT,NQT,NZZT(1,NRLIST(1),nx_ext),
     &              IWRIT4(NRLIST(1),nx_ext),PRECON_CODE(nx_ext),
     &              SOLVEROPTION(nx_ext),SPARSEGKK(nx_ext),
     &              GKK(1,nx_ext),RHS,YQ(1,niqV,1,nx_ext),FIRST_A,
     &              UPDATE_MATRIX,.FALSE.,nx_ext,ERROR,*9999)
               UPDATE_MATRIX=.FALSE.
               !
            ELSE
               ! Use Boxmg package to solve for tensor product mesh
               CALL BOXMG_SOLVE(
     &               ITYP4(NRLIST(1),nx),
     &               ISC_GKK(1,nx_ext),ISR_GKK(1,nx_ext),NQT,NQT,NQT,
     &               NZZT(1,NRLIST(1),nx_ext),IWRIT4(NRLIST(1),nx_ext),
     &               SOLVEROPTION(nx_ext),SPARSEGKK(nx_ext),
     &               GKK(1,nx_ext),RHS,YQ(1,niqV,1,nx_ext),FIRST_A,
     &               UPDATE_MATRIX,X_INIT,NEELEM,NQNE,NQS,NQXI,NXI,NR,
     &               nx_ext,ERROR,*9999)
               UPDATE_MATRIX=.FALSE.
               !
            ENDIF
            !
          ENDIF !bidomain

          
          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP       CRITICAL(MARCH8_4)
            WRITE(OP_STRING,'('' Time: '',F12.6)') T
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            DO nq=1,NQT
              WRITE(OP_STRING,'('' nq,Vm'',I6,F12.6)') nq,
     &          YQ(nq,niqV,1,nx)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
CC$OMP       END CRITICAL(MARCH8_4)
          ENDIF !DOP

          IF(KTYP36.EQ.1) THEN
            ! remove points
            CALL CALC_DTAR(maqp1t0,maqp1t1,maqp1i,maqp2t0,maqp2t1,
     &        maqp2i,niq_old,niqV,nnq_min,NQXI,NRLIST,NSOL,NWQ(1,0,1),
     &        nx,NXQ(-NIM,0,0,1),AQ,CQ,T,YQ(1,1,1,nx),.FALSE.,
     &        ERROR,*9999)
          ENDIF !DTAR

          IF(ITYP19(NRLIST(1),nx).EQ.1.AND.ITYP3(NRLIST(1),nx).EQ.3
     &      .AND.KTYP33.EQ.2) THEN !VCDC
C$OMP       PARALLEL DO
C$OMP&        DEFAULT(none)
C$OMP&        PRIVATE(nq)
C$OMP&        SHARED(CELL_STATE_OFFSET,DT,niq_old,niqV,NQT,nx,YQ,YQS)
            DO nq=1,NQT
              YQS(CELL_STATE_OFFSET(1)+DVm-1,nq)=(YQ(nq,niqV,1,nx)
     &          -YQ(nq,niq_old,1,nx))/DT
            ENDDO
C$OMP       END PARALLEL DO
          ENDIF

!        Store solution as basis for next time step
          IF(ITYP19(NRLIST(1),nx).EQ.1.OR.ITYP19(NRLIST(1),nx).EQ.7)
     &      THEN !electrical or coupled model
            IF(CALC_ACTIV_TIMES) THEN
              DO nq=NQR(1,nr),NQR(2,nr)
                IF((YQ(nq,niqV,1,nx)-YQ(nq,niq_OLD,1,nx)).GT.
     &            AQ(dpot,nq)) THEN
                  AQ(dpot,nq)=YQ(nq,niqV,1,nx)-YQ(nq,niq_OLD,1,nx)
                  AQ(atime,nq)=T
                ENDIF
              ENDDO
            ENDIF
          ENDIF

          T=T+DT

C Adaptive time step update for explicit adams
          IF((KTYP37.EQ.5).AND.(KTYP23.EQ.2)) THEN !Adams & adaptive dt
            CALL ASSERT(THETA(1).LT.ZERO_TOL,' >>Must use theta=0.0',
     &        ERROR,*9999)
C           Find the new smallest time step
            OLDDT=DT
            IF(NWQ(5,0,1).GT.0) THEN
              DT=RMAX
C KAT: not much point parallelizing this unless we use a min reduction
C operator as most of the loop is in a critical section.
C C$OMP         PARALLEL DO
C C$OMP&          DEFAULT(none)
C C$OMP&          PRIVATE(nq,nnq)
C C$OMP&          SHARED(INTEGRATOR_WORK,NWQ,DT,nx)
              DO nnq=1,NWQ(5,0,1)
                nq=NWQ(5,nnq,1)
                IF(NWQ(1,nq,1).EQ.0) THEN
C C$OMP             CRITICAL(MARCH8_NEW_4)
                  IF(INTEGRATOR_WORK(1,nq).LT.DT) THEN
                    DT=INTEGRATOR_WORK(1,nq)
                  ENDIF
C C$OMP             END CRITICAL(MARCH8_NEW_4)
                ENDIF
              ENDDO
C C$OMP         END PARALLEL DO
            ENDIF
            TINCR=DT
          ENDIF

          IF(IWRIT5(NRLIST(1),nx).EQ.2) THEN
            IF(MOD(NSTEP,IWRIT1(NRLIST(1),nx)).EQ.0) THEN
              CALL CPU_TIMER(CPU_USER,TIME_STOP2)
              ELAPSED_TIME=TIME_STOP2(1)-TIME_START2(1)
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP         CRITICAL(MARCH8_5)
              WRITE(OP_STRING,'(1X,''Time: '',E12.4,'', DT: '',E12.4,'
     &          //''' Iters: '',I8,'' Active: '',I6,'' Cpu: '',F6.2)')
     &          T,DT,NSTEP,NWQ(4,0,1),ELAPSED_TIME
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
CC$OMP         END CRITICAL(MARCH8_5)
              CALL CPU_TIMER(CPU_USER,TIME_START2)
            ENDIF
          ENDIF

          IF(T.GE.TFINISH) CONTINU=.FALSE.
        ENDDO !time
! End on main time integration loop

        T_RESTART(nx)=T

        IF(IWRIT5(NRLIST(1),nx).GE.1) THEN
          CALL CPU_TIMER(CPU_USER,TIME_STOP)
          ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP     CRITICAL(MARCH8_6)
          WRITE(OP_STRING,'(1X,I6,'' iterations : '',F8.2,'
     &      //'''s cpu'')') NSTEP,ELAPSED_TIME
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
CC$OMP     END CRITICAL(MARCH8_6)
        ENDIF

C------------------- coupled bidomain/bdry elements ------------------
      ELSE !iterate
C#### Comment: BIDOMAIN ITERATION
C###  Description:
C###    <HTML> <PRE>
C###    You must specify 3 regions and 5 classes for coupled
C###    bidomain-torso problem solutions.
C###
C###    REGIONS: The first region is the region where grid points
C###               are defined. (nr_grid)
C###             The second region is the region which encloses the
C###               grid region. (nr_torso)
C###             The third region is the region enclosed by the
C###               grid region. (nr_blood)
C###
C###    CLASSES: The first class is the transmembrane potential
C###               solution class.
C###             The second class is the extracellular potential
C###               soulution class.
C###             The third class is the extracellular iteration
C###               class.
C###             The fourth class is the class belonging to
C###               nr_torso and to classes outside this region.
C###             The fifth class is the class belonging to
C###               nr_blood and to classes inside this region.
C###
C###    MODELLING ASSUMPTIONS:
C###      There are several things which are assumed about the
C###      mesh on which coupled bidomain iterations are performed.
C###
C###      1: The grid (heart) regions always returns outward normals
C###         and the blood,torso cavity regions also return outward
C###         normals so there is a sign reversal in the code when
C###         moving between the grid region and the surrounding
C###         boundary element regions.
C###
C###      2: The arc length derivative on the boundary between the
C###         heart and the torso cavity is consistent. There is
C###         no arc length direction reversal between the heart and
C###         the torso cavity.
C###
C###      3: The blood region inside the heart is always a lower
C###         region number than the heart muscle. This allows the
C###         blood problem to be solved independently. If this is
C###         not the case then the blood region normals and arc
C###         length derivatives would be incorrect.
C###   </PRE> </HTML>
C**** Martin Buist 30-Oct-1998

C#### Comment: EXTRACELLULAR STIMULATION
C###  Description:
C###    <HTML> <PRE>
C###    You must define a time variable using the
C###
C###    FEM DEFINE TIME ...
C###
C###    command before defining any extracellular stimuli in
C###    the initial conditions file. The time variable will
C###    control the time course of the stimulus.
C###
C###    Internal points:
C###    ----------------
C###      All internal points are current injections.
C###
C###    Boundary points:
C###    ----------------
C###      What happens at a boundary point is determined by the
C###      type of boundary condition which has been applied.
C###      Current injection occurs at FLUX boundary points only.
C###      Setting a time variable at a potential boundary point
C###      sets that point to follow the specified time course
C###      where the boundary potential is equal to the time
C###      variable value.
C###   </PRE> </HTML>
C**** Martin Buist 22-Nov-1999

        CALL ASSERT(NXLIST(0).GE.5,
     &    '>>You need 5 classes for bidomain iterations',ERROR,*9999)
        CALL ASSERT(NRLIST(0).GE.3,
     &    '>>You need 3 regions for bidomain iterations',ERROR,*9999)

        !Define the required classes
        nxc=NXLIST(4)
        CALL NX_LOC(NX_INQUIRE,nxc,nx_torso,NX_SOLVE,ERROR,*9999)
C        CALL ASSERT(nx_torso.GT.0,
C     '    '>>No nx defined for this solve class',ERROR,*9999)

        nxc=NXLIST(5)
        CALL NX_LOC(NX_INQUIRE,nxc,nx_blood,NX_SOLVE,ERROR,*9999)
C        CALL ASSERT(nx_blood.GT.0,
C     '    '>>No nx defined for this solve class',ERROR,*9999)

        IF((nx_blood.EQ.0).AND.(nx_torso.EQ.0)) THEN
          CALL ASSERT(.FALSE.,' >>A BEM region class must be specified',
     &      ERROR,*9999)
        ENDIF

        !Define the required regions
        CALL ASSERT(NRLIST(0).GE.3,
     &    '>>You need 3 regions for solve8 iterations',ERROR,*9999)
        nr_grid=NRLIST(1)
        nr_torso=NRLIST(2)
        nr_blood=NRLIST(3)

        !Generate the extracellular update matrix
        IF(FIRSTITER) THEN
          FIRSTITER=.FALSE.
          NRTEMP(0)=1
          NRTEMP(1)=NRLIST(1)
          IF(SPARSEGKK(nx_upd).EQ.4) SPARSEGKK(nx_upd)=2

C SGM 13 Nov 2000 call ASSEMBLE10_FE for grid-based FE
C MLT 29Nov02 call ASSEMBLE10_FE for grid FV also
          IF(ITYP4(NRLIST(1),nx).EQ.4) THEN !collocation
            CALL ASSEMBLE10(ISC_GKK,ISR_GKK,NEELEM,NENQ,NLATNE,NLATNQ,
     &        NLATPNQ,NLQ,NQGP,NQGP_PIVOT,NQNLAT,NQS,NQXI,NRTEMP,
     &        NWQ(1,0,1),nx_ext,nx,nx_upd,NXQ(-NIM,0,0,1),AQ,CQ,DNUDXQ,
     &        DXDXIQ,DXDXIQ2,GCHQ,GUQ,GKK,NQGW,PROPQ,XQ,.FALSE.,
     &        .TRUE.,FIRST_A,FIXQ,.FALSE.,.FALSE.,.FALSE.,UPDATE_MATRIX,
     &        .FALSE.,ERROR,*9999)
          ELSEIF(ITYP4(NRLIST(1),nx).EQ.6.OR.
     &           ITYP4(NRLIST(1),nx).EQ.7) THEN !grid-based FE/grid FV
            CALL ASSEMBLE10_FE(ISC_GKK,ISR_GKK,NEELEM,NLATNE,NQGP,
     &        NQGP_PIVOT,NQNLAT,NQS,NQSCNB,NQXI,NRLIST,nx_ext,nx,nx_upd,
     &        NXQ(-NIM,0,0,1),CQ,GKK,GM,NQGW,PG,WG,XQ,PROPQ,.FALSE.,
     &        .TRUE.,FIRST_A,FIXQ,.FALSE.,.FALSE.,.FALSE.,UPDATE_MATRIX,
     &        .FALSE.,ERROR,*9999)
          ENDIF
          FIRST_A=.TRUE.
          UPDATE_MATRIX=.TRUE.
          X_INIT=.TRUE.
        ELSE
          UPDATE_MATRIX=.FALSE.
          FIRST_A=.FALSE.
          X_INIT=.FALSE.
        ENDIF

        !Generate the RHS vector from YP
        CALL GEN_GRID_POTE_RHS(IBT,IDO,INP,NBH,NBJ,NEELEM,NENQ,NHE,NHP,
     &    niqV,NKHE,NKH,NPF,NPNE,NPNODE,NQGP,NQGP_PIVOT,NQS,NQXI,
     &    nr_blood,nr_grid,nr_torso,NVHE,NVHP,NW,NWQ(1,0,1),nx_blood,
     &    nx_ext,nx_torso,nx,nx_upd,NXQ(-NIM,0,0,1),NYNE,NYNP,
     &    ALPHA,AQ,CQ,CURVCORRECT,DNUDXQ,DXDXIQ,DXDXIQ2,NQGW,RHS,SE,XIQ,
     &    YP,YQ,ZA,ZE,ZP,FIXQ,ERROR,*9999)

!       Solve system of equations
        IF(SOLVEROPTION(nx_ext).NE.SOLV_BMG .AND. 
     &     SOLVEROPTION(nx_ext).NE.SOLV_PCG_BMG ) THEN
           ! Use Standard Solver library
           CALL SOLVE_SYSTEM(ISC_GKK(1,nx_upd),ISR_GKK(1,nx_upd),NQT,
     &          NQT,NQT,NZZT(1,NRLIST(1),nx_upd),
     &          IWRIT4(NRLIST(1),nx_upd),PRECON_CODE(nx_upd),
     &          SOLVEROPTION(nx_upd),SPARSEGKK(nx_upd),GKK(1,nx_upd),
     &          RHS,YQ(1,niqV,1,nx_ext),FIRST_A,UPDATE_MATRIX,X_INIT,
     &          nx_upd,ERROR,*9999)
           !
        ELSE
           ! Use Boxmg package to solve for tensor product mesh
           CALL BOXMG_SOLVE(
     &          ITYP4(NRLIST(1),nx),
     &          ISC_GKK(1,nx_upd),ISR_GKK(1,nx_upd),NQT,NQT,NQT,
     &          NZZT(1,NRLIST(1),nx_upd),IWRIT4(NRLIST(1),nx_upd),
     &          SOLVEROPTION(nx_upd),SPARSEGKK(nx_upd),
     &          GKK(1,nx_upd),RHS,YQ(1,niqV,1,nx_ext),FIRST_A,
     &          UPDATE_MATRIX,X_INIT,NEELEM,NQNE,NQS,NQXI,NXI,nr,
     &          nx_upd,ERROR,*9999)
           !
        ENDIF
        !
      ENDIF !monodomain/bidomain

      CALL EXITS('MARCH8')
      RETURN
 9999 CALL ERRORS('MARCH8',ERROR)
      RET_ERROR=ERROR
      CALL EXITS('MARCH8')
      RETURN 1
      END
