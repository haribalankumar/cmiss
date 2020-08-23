      SUBROUTINE SOLVE8(IBT,IDO,INP,ISC_GKK,ISR_GKK,NBH,NBJ,NEELEM,NENQ,
     &  NLATNE,NLATNQ,NLATPNQ,NQNLAT,NHE,NHP,NKHE,NKH,NLQ,NPF,NPNE,
     &  NPNODE,NQGP,NQGP_PIVOT,NQS,NQSCNB,NQXI,NRLIST,NVHE,NVHP,NW,NWQ,
     &  nx,NXLIST,NXQ,NYNE,NYNP,ALPHA,AQ,CQ,CURVCORRECT,DNUDXQ,DXDXIQ,
     &  DXDXIQ2,GCHQ,GKK,GM,GRR,GUQ,NQGW,PG,PROPQ,RHS,SE,WG,XIQ,XQ,
     &  YP,YQ,ZA,ZE,ZP,FIXQ,ITER8,ERROR,*)

C#### Subroutine: SOLVE8
C###  Description:
C###    <HTML>
C###    <P>
C###    SOLVE8 solves Laplaces equation (standard)
C###    using collocation. Currently used to check against an
C###    analytic solution.
C###    </P> <P>
C###      <PRE>
C###      1d: u=x
C###      2d: u=x^2 - y^2
C###      3d: u=x^2 + y^2 - 2z^2
C###
C###      The solution is stored in iy=1 in YQ(nq,iy,na,nx).
C###      The absolute error is stored in iy=2.
C###      The relative error is stored in iy=3.
C###      </PRE>
C###    </P> <P>
C###    If the ITERATE option is used on the solve command then
C###    four classes and three regions are required as input in list
C###    form. Note that a standard solve must be done before any
C###    solve iterate commands are used.
C###    </P> <P>
C###      <PRE>
C###      class A,B,C,D
C###        where:
C###        A is the grid solution class
C###        B is the grid update class
C###        C is the torso solution class
C###        D is the ventricle solution class
C###
C###      region E,F,G
C###        where:
C###        E is the grid region
C###        F is the torso cavity region
C###        G is the ventricle region
C###      </PRE>
C###    </P>
C###    </HTML>
C***  Created by Martin Buist, June 1999

      IMPLICIT NONE

      INCLUDE 'b12.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'maqloc00.inc'
      INCLUDE 'solv00.cmn'
      INCLUDE 'time02.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter list
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     &  ISC_GKK(NISC_GKKM,NXM),ISR_GKK(NISR_GKKM,NXM),NBH(NHM,NCM,NEM),
     &  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NENQ(0:8,NQM),
     &  NLATNE(NEQM+1),NLATNQ(NEQM*NQEM),NLATPNQ(NQM),
     &  NQNLAT(NEQM*NQEM),NHE(NEM,NXM),
     &  NHP(NPM,0:NRM,NXM),NKHE(NKM,NNM,NHM,NEM),NKH(NHM,NPM,NCM,0:NRM),
     &  NLQ(NQM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     &  NQGP(0:NQGM,NQM),NQGP_PIVOT(NQGM,NQM),NQS(NEQM),NQSCNB(NQSCM),
     &  NQXI(0:NIM,NQSCM),NRLIST(0:NRM),NVHE(NNM,NBFM,NHM,NEM),
     &  NVHP(NHM,NPM,NCM,0:NRM),NW(NEM,3),NWQ(8,0:NQM),nx,NXLIST(0:NXM),
     &  NXQ(-NIM:NIM,0:4,0:NQM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     &  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 ALPHA,AQ(NMAQM,NQM),CQ(NMM,NQM),CURVCORRECT(2,2,NNM,NEM),
     &  DNUDXQ(3,3,NQM),DXDXIQ(3,3,NQM),DXDXIQ2(3,3,NQM),GCHQ(3,NQM),
     &  GKK(NZ_GKK_M,NXM),GM(NZ_GM_M),GRR(NOM),GUQ(3,3,NQM),
     &  NQGW(NQGM,NQM),PG(NSM,NUM,NGM,NBM),PROPQ(3,3,4,2,NQM),RHS(NQM),
     &  SE(NSM,NBFM,NEM),WG(NGM,NBM),XIQ(NIM,NQM),XQ(NJM,NQM),
     &  YP(NYM,NIYM,NXM),YQ(NYQM,NIQM,NAM,NXM),ZA(NAM,NHM,NCM,NEM),
     &  ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
      LOGICAL FIXQ(NYQM,NIYFIXM,NXM),ITER8
!     Local variables
      INTEGER COUNT,maq,niqV,nq,nr,nr_blood,nr_grid,NRTEMP(0:1),
     &  nr_torso,nxc,nx_ext,nx_blood,nx_torso,nx_upd
      REAL*8 AVE_ABS_ERR,AVE_REL_ERR,TOT_ABS_ERR,TOT_REL_ERR,XNLOCAL(3)
      LOGICAL FIRST_A,FIRSTITER,UPDATE_MATRIX,X_INIT

      SAVE FIRSTITER

      CALL ENTERS('SOLVE8',*9999)

      IF(.NOT.ITER8) THEN
        FIRSTITER=.TRUE.
!       Check the problem can be solved
        nr=NRLIST(1)
        CALL ASSERT(CALL_GRID,' >>Define grid first',ERROR,*9999)

C DAH 24-JUL-2002 Lattice grid scheme does not need these
        IF(ITYP4(nr,nx).NE.6.AND.USE_LAT.EQ.0) THEN
           CALL ASSERT(UP_GRID_TENSOR,' >>Update grid metrics first',
     &          ERROR,*9999)
           CALL ASSERT(UP_GRID_MATERIAL,' >>Update grid material first',
     &          ERROR,*9999)
        ENDIF

        CALL ASSERT(CALL_INIT,' >>Define initial first',ERROR,*9999)
        CALL ASSERT(CALL_SOLV,' >>Define solve first',ERROR,*9999)
        CALL ASSERT(ITYP9(nr,nx).EQ.1,' >>Must choose direct solve',
     &    ERROR,*9999)
        CALL ASSERT(ITYP5(nr,nx).EQ.1,' >>Must choose static equation',
     &    ERROR,*9999)
        CALL ASSERT(ITYP2(nr,nx).EQ.3,
     &    ' >>Must choose Laplaces equation',ERROR,*9999)
C SGM 26 Oct 2000 grid-based Finite element also
C MLT 10Oct02 grid FV also
        CALL ASSERT(ITYP4(nr,nx).EQ.4.OR.ITYP4(nr,nx).EQ.6.OR.
     &              ITYP4(nr,nx).EQ.7,
     &    ' >>Must choose collocation, grid-based FE or grid FV soln'
     &    ,ERROR,*9999)
        CALL ASSERT(NIQM.GE.3,' >>NIQM must be at least 3',ERROR,*9999)
        CALL ASSERT(NOM.GE.NQM,' >>NOM must be >= NQM',ERROR,*9999)

C SGM 1 Nov 2000 FIXQ not initialised as yet
        DO nq=1,NQT
          IF(NWQ(1,nq).NE.0) THEN
            FIXQ(nq,1,nx)=.TRUE.
          ELSE
            FIXQ(nq,1,nx)=.FALSE.
          END IF
          FIXQ(nq,3,nx)=.FALSE. !not a monodomain approx b.c.
        ENDDO
        IF(SPARSEGKK(nx).EQ.4) SPARSEGKK(nx)=2

!       Create matrix
C SGM 13 Nov 2000 call ASSEMBLE10_FE for grid-based FE
C MLT 29Nov02 call ASSEMBLE10_FE for grid FV also
        IF(ITYP4(nr,nx).EQ.4) THEN !collocation
          CALL ASSEMBLE10(ISC_GKK,ISR_GKK,NEELEM,NENQ,NLATNE,NLATNQ,
     &      NLATPNQ,NLQ,NQGP,NQGP_PIVOT,NQNLAT,NQS,NQXI,NRLIST,NWQ,nx,
     &      nx,nx,NXQ,AQ,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,GCHQ,GUQ,GKK,NQGW,
     &      PROPQ,XQ,.TRUE.,.FALSE.,FIRST_A,FIXQ,.FALSE.,.FALSE.,
     &      .FALSE.,UPDATE_MATRIX,.TRUE.,ERROR,*9999)
        ELSEIF(ITYP4(nr,nx).EQ.6.OR.ITYP4(nr,nx).EQ.7) THEN !grid-based FE
          CALL ASSEMBLE10_FE(ISC_GKK,ISR_GKK,NEELEM,NLATNE,NQGP,
     &      NQGP_PIVOT,NQNLAT,NQS,NQSCNB,NQXI,NRLIST,nx,0,0,NXQ,CQ,GKK,
     &      GM,NQGW,PG,WG,XQ,PROPQ,.TRUE.,.FALSE.,FIRST_A,FIXQ,.FALSE.,
     &      .FALSE.,.FALSE.,UPDATE_MATRIX,.TRUE.,ERROR,*9999)
        ENDIF

!       Create RHS vector
        IF(NIT(NBJ(1,NEELEM(1,nr))).EQ.1) THEN
          DO nq=1,NQT
C!!! DAH 25-JUL-2002 this boundary stuff in only for testing
C!!!     This is very suspect, as NWQ seems poorly calculated.
C!!!     check thoroughly.            
            IF (USE_LAT.EQ.0) THEN
              IF(NWQ(1,nq).EQ.0) THEN
                GRR(nq)=0.0d0
              ELSE
                GRR(nq)=XQ(1,nq)
              ENDIF
            ELSE
              IF(NWQ(1,nq).EQ.0) THEN
                GRR(nq)=0.0d0
              ELSE
                GRR(nq)=XQ(1,nq)
              ENDIF
            ENDIF
          ENDDO
        ELSE IF(NIT(NBJ(1,NEELEM(1,nr))).EQ.2) THEN
          IF(SOLVE8_FLUXBC) THEN
            GRR(1)=(XQ(1,1)**2.0d0)-(XQ(2,1)**2.0d0)
            DO nq=2,NQT
              IF(NWQ(1,nq).EQ.0) THEN
                GRR(nq)=0.0d0
              ELSE
                IF(USE_LAT.EQ.0) THEN
                  CALL NORM31(NJT,nq,NXQ,DXDXIQ,DXDXIQ2,
     &              XNLOCAL,ERROR,*9999)
                ELSE
C                  CALL NORM_LATTICE(NJT,nq,NWQ(1,nq),DXDXIQ,DXDXIQ2,
C     '              XNLOCAL,ERROR,*9999)
                  CALL MAQ_LOC(MAQ_INQUIRE,MAQ_COORD,maq,MAQ_NORMAL_X,
     &              ERROR,*9999)
                  XNLOCAL(1)=AQ(maq,nq)
                  CALL MAQ_LOC(MAQ_INQUIRE,MAQ_COORD,maq,MAQ_NORMAL_Y,
     &              ERROR,*9999)
                  XNLOCAL(2)=AQ(maq,nq)
                  CALL MAQ_LOC(MAQ_INQUIRE,MAQ_COORD,maq,MAQ_NORMAL_Z,
     &              ERROR,*9999)
                  XNLOCAL(3)=AQ(maq,nq)
                ENDIF
                GRR(nq)=(XQ(1,nq)*2.0d0*XNLOCAL(1))-
     &            (XQ(2,nq)*2.0d0*XNLOCAL(2))
              ENDIF
            ENDDO
          ELSE
            DO nq=1,NQT
              IF(NWQ(1,nq).EQ.0) THEN
                GRR(nq)=0.0d0
              ELSE
                GRR(nq)=(XQ(1,nq)**2.0d0)-(XQ(2,nq)**2.0d0)
              ENDIF
            ENDDO
          ENDIF
        ELSE IF(NIT(NBJ(1,NEELEM(1,nr))).EQ.3) THEN
          DO nq=1,NQT
            IF(NWQ(1,nq).EQ.0) THEN
              GRR(nq)=0.0d0
            ELSE
              GRR(nq)=(XQ(1,nq)**2.0d0)+(XQ(2,nq)**2.0d0)-
     &          (2.0d0*XQ(3,nq)**2.0d0)
            ENDIF
          ENDDO
        ENDIF

!       Solve the system
        FIRST_A=.TRUE.
        UPDATE_MATRIX=.TRUE.
        X_INIT=.TRUE.

        CALL SOLVE_SYSTEM(ISC_GKK(1,nx),ISR_GKK(1,nx),NQT,NQT,NQT,
     &    NZZT(1,NRLIST(1),nx),IWRIT4(NRLIST(1),nx),PRECON_CODE(nx),
     &    SOLVEROPTION(nx),SPARSEGKK(nx),GKK(1,nx),GRR,YQ(1,1,1,nx),
     &    FIRST_A,UPDATE_MATRIX,X_INIT,nx,ERROR,*9999)

!       Check the solution
        COUNT=NQT
        TOT_ABS_ERR=0.0d0
        TOT_REL_ERR=0.0d0
        IF(NIT(NBJ(1,NEELEM(1,nr))).EQ.1) THEN
          WRITE(OP_STRING,'(/'' Analytic solution is u = x'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          DO nq=1,NQT
            YQ(nq,2,1,nx)=DABS(YQ(nq,1,1,nx)-XQ(1,nq))
            TOT_ABS_ERR=TOT_ABS_ERR+YQ(nq,2,1,nx)
            IF(DABS(XQ(1,nq)).GT.LOOSE_TOL) THEN
              YQ(nq,3,1,nx)=DABS((YQ(nq,1,1,nx)-XQ(1,nq))/XQ(1,nq))
              TOT_REL_ERR=TOT_REL_ERR+YQ(nq,3,1,nx)
            ELSE
              COUNT=COUNT-1
              YQ(nq,3,1,nx)=0.0d0
            ENDIF
          ENDDO
        ELSE IF(NIT(NBJ(1,NEELEM(1,nr))).EQ.2) THEN
          WRITE(OP_STRING,'(/'' Analytic solution is u = x^2 - y^2'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          DO nq=1,NQT
            YQ(nq,2,1,nx)=
     &        DABS(YQ(nq,1,1,nx)-((XQ(1,nq)**2.0d0)-(XQ(2,nq)**2.0d0)))
            TOT_ABS_ERR=TOT_ABS_ERR+YQ(nq,2,1,nx)
            IF(DABS((XQ(1,nq)**2.0d0)-(XQ(2,nq)**2.0d0)).GT.LOOSE_TOL)
     &       THEN
             YQ(nq,3,1,nx)=
     &          DABS((YQ(nq,1,1,nx)-((XQ(1,nq)**2.0d0)-
     &          (XQ(2,nq)**2.0d0)))/((XQ(1,nq)**2.0d0)-
     &          (XQ(2,nq)**2.0d0)))
              TOT_REL_ERR=TOT_REL_ERR+YQ(nq,3,1,nx)
            ELSE
              COUNT=COUNT-1
              YQ(nq,3,1,nx)=0.0d0
            ENDIF
          ENDDO
        ELSE IF(NIT(NBJ(1,NEELEM(1,nr))).EQ.3) THEN
          WRITE(OP_STRING,'(/'' Analytic solution is u ='
     &      //' x^2 + y^2 - 2z^2'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          DO nq=1,NQT
            YQ(nq,2,1,nx)=
     &        DABS(YQ(nq,1,1,nx)-((XQ(1,nq)**2.0d0)+(XQ(2,nq)**2.0d0)-
     &        (2.0d0*XQ(3,nq)**2.0d0)))
            TOT_ABS_ERR=TOT_ABS_ERR+YQ(nq,2,1,nx)
            IF(DABS((XQ(1,nq)**2.0d0)+(XQ(2,nq)**2.0d0)-
     &        (2.0d0*XQ(3,nq)**2.0d0)).GT.LOOSE_TOL) THEN
              YQ(nq,3,1,nx)=
     &          DABS((YQ(nq,1,1,nx)-((XQ(1,nq)**2.0d0)+
     &          (XQ(2,nq)**2.0d0)-
     &          (2.0d0*XQ(3,nq)**2.0d0)))/((XQ(1,nq)**2.0d0)+
     &          (XQ(2,nq)**2.0d0)-(2.0d0*XQ(3,nq)**2.0d0)))
              TOT_REL_ERR=TOT_REL_ERR+YQ(nq,3,1,nx)
            ELSE
              COUNT=COUNT-1
              YQ(nq,3,1,nx)=0.0d0
            ENDIF
          ENDDO
        ENDIF

        AVE_ABS_ERR=TOT_ABS_ERR/DBLE(NQT)
        AVE_REL_ERR=TOT_REL_ERR/DBLE(COUNT)

        WRITE(OP_STRING,'(/'' Average absolute error '',E12.6)')
     &    AVE_ABS_ERR
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Average relative error '',E12.6)')
     &    AVE_REL_ERR
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      ELSE !iterate

!       Define the required classes
        CALL ASSERT(NXLIST(0).GE.4,
     &    '>>You need 4 classes for solve8 iterations',ERROR,*9999)

        nxc=NXLIST(1) !potential
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     &    ERROR,*9999)
        nx_ext=nx

        nxc=NXLIST(2) !potential update
        CALL NX_LOC(NX_INQUIRE,nxc,nx_upd,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx_upd.GT.0,'>>No nx defined for this solve class',
     &    ERROR,*9999)

        nxc=NXLIST(3) !torso
        CALL NX_LOC(NX_INQUIRE,nxc,nx_torso,NX_SOLVE,ERROR,*9999)
C        CALL ASSERT(nx_torso.GT.0,
C     '    '>>No nx defined for this solve class',ERROR,*9999)

        nxc=NXLIST(4) !blood
        CALL NX_LOC(NX_INQUIRE,nxc,nx_blood,NX_SOLVE,ERROR,*9999)
C        CALL ASSERT(nx_blood.GT.0,
C     '    '>>No nx defined for this solve class',ERROR,*9999)

!       Define the required regions
        CALL ASSERT(NRLIST(0).GE.3,
     &    '>>You need 3 regions for solve8 iterations',ERROR,*9999)
        nr_grid=NRLIST(1)
        nr_torso=NRLIST(2)
        nr_blood=NRLIST(3)
        nr=nr_grid

!       Create matrix
        IF(FIRSTITER) THEN
          FIRSTITER=.FALSE.
          NRTEMP(0)=1
          NRTEMP(1)=NRLIST(1)
          IF(SPARSEGKK(nx_upd).EQ.4) SPARSEGKK(nx_upd)=2

C SGM 13 Nov 2000 call ASSEMBLE10_FE for grid-based FE
C MLT 29Nov02 call ASSEMBLE10_FE for grid FV also
          IF(ITYP4(nr,nx).EQ.4) THEN !collocation
            CALL ASSEMBLE10(ISC_GKK,ISR_GKK,NEELEM,NENQ,NLATNE,NLATNQ,
     &        NLATPNQ,NLQ,NQGP,NQGP_PIVOT,NQNLAT,NQS,NQXI,NRTEMP,NWQ,
     &        nx_ext,nx,nx_upd,NXQ,AQ,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,GCHQ,GUQ,
     &        GKK,NQGW,PROPQ,XQ,.FALSE.,.TRUE.,FIRST_A,FIXQ,.FALSE.,
     &        .FALSE.,.FALSE.,UPDATE_MATRIX,.FALSE.,ERROR,*9999)
          ELSEIF(ITYP4(nr,nx).EQ.6.OR.ITYP4(nr,nx).EQ.7) THEN !grid-based FE/FV
            CALL ASSEMBLE10_FE(ISC_GKK,ISR_GKK,NEELEM,NLATNE,NQGP,
     &        NQGP_PIVOT,NQNLAT,NQS,NQSCNB,NQXI,NRTEMP,nx_ext,nx,nx_upd,
     &        NXQ,CQ,GKK,GM,NQGW,PG,WG,XQ,PROPQ,.FALSE.,.TRUE.,FIRST_A,
     &        FIXQ,.FALSE.,.FALSE.,.FALSE.,UPDATE_MATRIX,.FALSE.,
     &        ERROR,*9999)
          ENDIF

          FIRST_A=.TRUE.
          UPDATE_MATRIX=.TRUE.
          X_INIT=.TRUE.
        ELSE
          FIRST_A=.FALSE.
          UPDATE_MATRIX=.FALSE.
          X_INIT=.FALSE.
        ENDIF

!       Generate the RHS vector from YP
        niqV=1
        CALL GEN_GRID_POTE_RHS(IBT,IDO,INP,NBH,NBJ,NEELEM,NENQ,NHE,NHP,
     &    niqV,NKHE,NKH,NPF,NPNE,NPNODE,NQGP,NQGP_PIVOT,NQS,NQXI,
     &    nr_blood,nr_grid,nr_torso,NVHE,NVHP,NW,NWQ,nx_blood,nx,
     &    nx_torso,nx,nx_upd,NXQ,NYNE,NYNP,ALPHA,AQ,CQ,CURVCORRECT,
     &    DNUDXQ,DXDXIQ,DXDXIQ2,NQGW,RHS,SE,XIQ,YP,YQ,ZA,ZE,ZP,FIXQ,
     &    ERROR,*9999)

!       Solve system of equations
        CALL SOLVE_SYSTEM(ISC_GKK(1,nx_upd),ISR_GKK(1,nx_upd),NQT,
     &    NQT,NQT,NZZT(1,NRLIST(1),nx_upd),IWRIT4(NRLIST(1),nx_upd),
     &    PRECON_CODE(nx_upd),SOLVEROPTION(nx_upd),SPARSEGKK(nx_upd),
     &    GKK(1,nx_upd),RHS,YQ(1,niqV,1,nx_ext),FIRST_A,UPDATE_MATRIX,
     &    X_INIT,nx_upd,ERROR,*9999)

!       Check the solution
        COUNT=NQT
        TOT_ABS_ERR=0.0d0
        TOT_REL_ERR=0.0d0
        IF(NIT(NBJ(1,NEELEM(1,nr))).EQ.1) THEN
          WRITE(OP_STRING,'(/'' Analytic solution is u = x'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          DO nq=1,NQT
            YQ(nq,2,1,nx)=DABS(YQ(nq,1,1,nx)-XQ(1,nq))
            TOT_ABS_ERR=TOT_ABS_ERR+YQ(nq,2,1,nx)
            IF(DABS(XQ(1,nq)).GT.LOOSE_TOL) THEN
              YQ(nq,3,1,nx)=DABS((YQ(nq,1,1,nx)-XQ(1,nq))/XQ(1,nq))
              TOT_REL_ERR=TOT_REL_ERR+YQ(nq,3,1,nx)
            ELSE
              COUNT=COUNT-1
              YQ(nq,3,1,nx)=0.0d0
            ENDIF
          ENDDO
        ELSE IF(NIT(NBJ(1,NEELEM(1,nr))).EQ.2) THEN
          WRITE(OP_STRING,'(/'' Analytic solution is u = x^2 - y^2'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          DO nq=1,NQT
            YQ(nq,2,1,nx)=
     &        DABS(YQ(nq,1,1,nx)-((XQ(1,nq)**2.0d0)-(XQ(2,nq)**2.0d0)))
            TOT_ABS_ERR=TOT_ABS_ERR+YQ(nq,2,1,nx)
            IF(DABS((XQ(1,nq)**2.0d0)-(XQ(2,nq)**2.0d0)).GT.LOOSE_TOL)
     &       THEN
             YQ(nq,3,1,nx)=
     &          DABS((YQ(nq,1,1,nx)-((XQ(1,nq)**2.0d0)-
     &          (XQ(2,nq)**2.0d0)))/((XQ(1,nq)**2.0d0)-
     &          (XQ(2,nq)**2.0d0)))
              TOT_REL_ERR=TOT_REL_ERR+YQ(nq,3,1,nx)
            ELSE
              COUNT=COUNT-1
              YQ(nq,3,1,nx)=0.0d0
            ENDIF
          ENDDO
        ELSE IF(NIT(NBJ(1,NEELEM(1,nr))).EQ.3) THEN
          WRITE(OP_STRING,'(/'' Analytic solution is u ='
     &      //' x^2 + y^2 - 2z^2'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          DO nq=1,NQT
            YQ(nq,2,1,nx)=
     &        DABS(YQ(nq,1,1,nx)-((XQ(1,nq)**2.0d0)+(XQ(2,nq)**2.0d0)-
     &        (2.0d0*XQ(3,nq)**2.0d0)))
            TOT_ABS_ERR=TOT_ABS_ERR+YQ(nq,2,1,nx)
            IF(DABS((XQ(1,nq)**2.0d0)+(XQ(2,nq)**2.0d0)-
     &        (2.0d0*XQ(3,nq)**2.0d0)).GT.LOOSE_TOL) THEN
              YQ(nq,3,1,nx)=
     &          DABS((YQ(nq,1,1,nx)-((XQ(1,nq)**2.0d0)+
     &          (XQ(2,nq)**2.0d0)-
     &          (2.0d0*XQ(3,nq)**2.0d0)))/((XQ(1,nq)**2.0d0)+
     &          (XQ(2,nq)**2.0d0)-(2.0d0*XQ(3,nq)**2.0d0)))
              TOT_REL_ERR=TOT_REL_ERR+YQ(nq,3,1,nx)
            ELSE
              COUNT=COUNT-1
              YQ(nq,3,1,nx)=0.0d0
            ENDIF
          ENDDO
        ENDIF

        AVE_ABS_ERR=TOT_ABS_ERR/DBLE(NQT)
        AVE_REL_ERR=TOT_REL_ERR/DBLE(COUNT)

        WRITE(OP_STRING,'(/'' Average absolute error '',E12.6)')
     &    AVE_ABS_ERR
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Average relative error '',E12.6)')
     &    AVE_REL_ERR
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF !not iterate

      CALL EXITS('SOLVE8')
      RETURN
 9999 CALL ERRORS('SOLVE8',ERROR)
      CALL EXITS('SOLVE8')
      RETURN 1
      END


