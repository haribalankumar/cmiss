      SUBROUTINE SOLVE8_COUPLED(ISC_GKK,ISR_GKK,NEELEM,NENQ,NLQ,
     '  NP_INTERFACE,NPLIST,NPNODE,NQGP,NQGP_PIVOT,NQLIST,NQNP,NQS,
     '  NQXI,NRLIST,NWQ,NXLIST,NXQ,NYNP,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,GCHQ,
     '  GK,GK2,GK3,GKK,GQ,GQ2,GQ3,GRR,GUQ,NQGW,PROPQ,XQ,YP,YQ,FIXQ,
     '  ERROR,*)

C#### Subroutine: SOLVE8_COUPLED
C###  Description:
C***  Created by Martin Buist, April 2000

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'solv00.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter list
      INTEGER ISC_GKK(NISC_GKKM,NXM),ISR_GKK(NISR_GKKM,NXM),
     '  NEELEM(0:NE_R_M,0:NRM),NENQ(0:8,NQM),NLQ(NQM),
     '  NP_INTERFACE(0:NPM,0:3),NPLIST(0:NPM),NPNODE(0:NP_R_M,0:NRM),
     '  NQGP(0:NQGM,NQM),NQGP_PIVOT(NQGM,NQM),NQLIST(0:NQM),NQNP(NPM),
     '  NQS(NEQM),NQXI(0:NIM,NQSCM),NRLIST(0:NRM),NWQ(8,0:NQM),
     '  NXLIST(0:NXM),NXQ(-NIM:NIM,0:4,0:NQM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CQ(NMM,NQM),DNUDXQ(3,3,NQM),DXDXIQ(3,3,NQM),
     '  DXDXIQ2(3,3,NQM),GCHQ(3,NQM),GK(NZ_GK_M),GK2(*),GK3(*),
     '  GKK(NZ_GKK_M,NXM),GQ(NZ_GQ_M),GQ2(*),GQ3(*),GRR(NOM),
     '  GUQ(3,3,NQM),NQGW(NQGM,NQM),PROPQ(3,3,4,2,NQM),XQ(NJM,NQM),
     '  YP(NYM,NIYM,NXM),YQ(NYQM,NIQM,NAM,NXM)
      CHARACTER ERROR*(*)
      LOGICAL FIXQ(NYQM,NIYFIXM,NXM)
!     Local variables
      INTEGER COUNT,nc,nh,np,nq,nqq,nr,nr2,nr3,nr4,NRTEMP(0:1),NSIZE,
     '  nv,nx,nx2,nx3,nx4,nxc,ny
      INTEGER*4 NQLIST2_PTR,NQLIST3_PTR,NQLIST4_PTR,
     '  NQNP_LOCAL_PTR,NQNP_LOCAL_PIV_PTR,D_MATRIX_PTR,
     '  GQ_D_MATRIX_PTR,GKGK2_PTR,GQGQ2_PTR
      REAL*8 AVE_ABS_ERR,AVE_REL_ERR,TOT_ABS_ERR,TOT_REL_ERR
      LOGICAL FIRST_A,UPDATE_MATRIX,X_INIT

      CALL ENTERS('SOLVE8_COUPLED',*9999)

C**** CPTYPE=1 is iteration on epi & endo
C**** CPTYPE=2 is coupled grid-endo, no iteration
C**** CPTYPE=3 is coupled grid-endo, epi iteration
C**** CPTYPE=4 is coupled grid-2 endo, no iteration
C**** CPTYPE=5 is coupled grid-2 endo, epi iteration
C**** CPTYPE=6 is coupled epi-grid-endo, no iteration
C**** CPTYPE=7 is coupled epi-grid-2 endo, no iteration

      nr=NRLIST(1)
      nr2=NRLIST(2)
      nxc=NXLIST(1)
      CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
      CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '  ERROR,*9999)
      nxc=NXLIST(2)
      CALL NX_LOC(NX_INQUIRE,nxc,nx2,NX_SOLVE,ERROR,*9999)
      CALL ASSERT(nx2.GT.0,'>>No nx defined for this solve class',
     '  ERROR,*9999)
      IF(CPTYPE.EQ.4) THEN
        nr3=NRLIST(3)
        nxc=NXLIST(3)
        CALL NX_LOC(NX_INQUIRE,nxc,nx3,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx3.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)
      ELSE IF(CPTYPE.EQ.6) THEN
        nr3=NRLIST(3)
        nxc=NXLIST(3)
        CALL NX_LOC(NX_INQUIRE,nxc,nx3,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx3.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)
      ELSE IF(CPTYPE.EQ.7) THEN
        nr3=NRLIST(3)
        nxc=NXLIST(3)
        CALL NX_LOC(NX_INQUIRE,nxc,nx3,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx3.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)
        nr4=NRLIST(4)
        nxc=NXLIST(4)
        CALL NX_LOC(NX_INQUIRE,nxc,nx4,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx4.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)
      ENDIF

C     Get interface grid points
      CALL PARSE_GRID(NQLIST,noco,NTCO,CO,ERROR,*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,'('' Number of grid points selected '',I8)')
     '    NQLIST(0)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

C     Calc grid point quadratics
      NRTEMP(0)=1
      NRTEMP(1)=nr
      CALL GET_FD_POINTS(NEELEM,NENQ,NLQ,NQGP,NQGP_PIVOT,NQS,NQXI,
     '  NRTEMP,NWQ,NXQ,ERROR,*9999)

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
     '  MEM_INIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(NQM+1,1,INTTYPE,NQLIST3_PTR,
     '  MEM_INIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(NQM+1,1,INTTYPE,NQLIST4_PTR,
     '  MEM_INIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(NPM,1,INTTYPE,NQNP_LOCAL_PTR,
     '  MEM_INIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(NPM,1,INTTYPE,NQNP_LOCAL_PIV_PTR,
     '  MEM_INIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(NPM*NPM*3,1,DPTYPE,D_MATRIX_PTR,
     '  MEM_INIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(NPM*NPM*3,1,DPTYPE,GQ_D_MATRIX_PTR,
     '  MEM_INIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(NPM*NPM,1,DPTYPE,GKGK2_PTR,
     '  MEM_INIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(NPM*NPM,1,DPTYPE,GQGQ2_PTR,
     '  MEM_INIT,ERROR,*9999)

      IF((CPTYPE.EQ.2).OR.(CPTYPE.EQ.3)) THEN
        CALL ASSEMBLE11_2_3(ISC_GKK,ISR_GKK,NEELEM,NENQ,NPNODE,NQGP,
     '    NQGP_PIVOT,NQLIST,%VAL(NQLIST2_PTR),%VAL(NQLIST3_PTR),NQNP,
     '    %VAL(NQNP_LOCAL_PTR),%VAL(NQNP_LOCAL_PIV_PTR),NQS,NQXI,nr,
     '    nr2,NWQ,nx,nx2,NXQ,CQ,%VAL(D_MATRIX_PTR),DNUDXQ,
     '    DXDXIQ,DXDXIQ2,GCHQ,GK,GKK,GQ,%VAL(GQ_D_MATRIX_PTR),
     '    %VAL(GKGK2_PTR),%VAL(GQGQ2_PTR),GUQ,NQGW,PROPQ,FIXQ,
     '    .TRUE.,ERROR,*9999)
      ELSE IF((CPTYPE.EQ.4).OR.(CPTYPE.EQ.5)) THEN
        CALL ASSEMBLE11_4_5(ISC_GKK,ISR_GKK,NEELEM,NENQ,NPNODE,NQGP,
     '    NQGP_PIVOT,NQLIST,%VAL(NQLIST2_PTR),%VAL(NQLIST3_PTR),
     '    %VAL(NQLIST4_PTR),NQNP,
     '    %VAL(NQNP_LOCAL_PTR),%VAL(NQNP_LOCAL_PIV_PTR),NQS,NQXI,nr,
     '    nr2,nr3,NWQ,nx,nx2,nx3,NXQ,CQ,%VAL(D_MATRIX_PTR),DNUDXQ,
     '    DXDXIQ,DXDXIQ2,GCHQ,GK,GK2,GKK,GQ,GQ2,%VAL(GQ_D_MATRIX_PTR),
     '    %VAL(GKGK2_PTR),%VAL(GQGQ2_PTR),GUQ,NQGW,PROPQ,FIXQ,
     '    .TRUE.,ERROR,*9999)
      ELSE IF(CPTYPE.EQ.6) THEN
        CALL ASSEMBLE11_6(ISC_GKK,ISR_GKK,NEELEM,NENQ,NP_INTERFACE,
     '    NPLIST,NPNODE,NQGP,NQGP_PIVOT,NQLIST,%VAL(NQLIST2_PTR),
     '    %VAL(NQLIST3_PTR),NQNP,%VAL(NQNP_LOCAL_PTR),
     '    %VAL(NQNP_LOCAL_PIV_PTR),NQS,NQXI,nr,nr2,nr3,NSIZE,NWQ,
     '    nx,nx2,nx3,NXQ,CQ,%VAL(D_MATRIX_PTR),DNUDXQ,DXDXIQ,
     '    DXDXIQ2,GCHQ,GK,GK2,GKK,GQ,GQ2,%VAL(GQ_D_MATRIX_PTR),
     '    %VAL(GKGK2_PTR),%VAL(GQGQ2_PTR),GUQ,NQGW,PROPQ,YQ,FIXQ,
     '    .TRUE.,ERROR,*9999)
      ELSE IF(CPTYPE.EQ.7) THEN
        CALL ASSEMBLE11_7(ISC_GKK,ISR_GKK,NEELEM,NENQ,NP_INTERFACE,
     '    NPLIST,NPNODE,NQGP,NQGP_PIVOT,NQLIST,%VAL(NQLIST2_PTR),
     '    %VAL(NQLIST3_PTR),%VAL(NQLIST4_PTR),NQNP,%VAL(NQNP_LOCAL_PTR),
     '    %VAL(NQNP_LOCAL_PIV_PTR),NQS,NQXI,nr,nr2,nr3,nr4,NSIZE,NWQ,
     '    nx,nx2,nx3,nx4,NXQ,CQ,%VAL(D_MATRIX_PTR),DNUDXQ,DXDXIQ,
     '    DXDXIQ2,GCHQ,GK,GK2,GK3,GKK,GQ,GQ2,GQ3,%VAL(GQ_D_MATRIX_PTR),
     '    %VAL(GKGK2_PTR),%VAL(GQGQ2_PTR),GUQ,NQGW,PROPQ,YQ,FIXQ,
     '    .TRUE.,ERROR,*9999)
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

      DO nq=1,NQT
        IF(NWQ(1,nq).EQ.0) THEN
          GRR(nq)=0.0d0
        ELSE
          GRR(nq)=(XQ(1,nq)**2.0d0)-(XQ(2,nq)**2.0d0)
        ENDIF
      ENDDO
      DO nqq=1,NQLIST(0)
        nq=NQLIST(nqq)
        GRR(nq)=0.0d0
      ENDDO !nq

      IF(DOP) THEN
        DO nq=1,NQT
          WRITE(OP_STRING,'(/I8,F12.6)') nq,GRR(nq)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

      FIRST_A=.TRUE.
      UPDATE_MATRIX=.TRUE.
      X_INIT=.TRUE.

      IF(CPTYPE.GE.6) THEN
        DO nq=1,NSIZE
          GRR(nq)=0.0d0
          IF(SOL_ACT_FIX_NODE.GT.0) THEN
            IF(NQNP(SOL_ACT_FIX_NODE).EQ.nq) GRR(nq)=5.0d0
          ENDIF
        ENDDO
        CALL SOLVE_SYSTEM(ISC_GKK(1,nx),ISR_GKK(1,nx),NSIZE,NSIZE,NSIZE,
     '    NZZT(1,nr,nx),IWRIT4(nr,nx),PRECON_CODE(nx),SOLVEROPTION(nx),
     '    SPARSEGKK(nx),GKK(1,nx),GRR,YQ(1,1,1,nx),FIRST_A,
     '    UPDATE_MATRIX,X_INIT,nx,ERROR,*9999)
        nqq=0
        !write node solution back into YP
        IF(CPTYPE.EQ.6) THEN
          DO nq=1,NPNODE(0,nr3)
            np=NPNODE(nq,nr3)
            nv=1
            nh=NH_LOC(1,nx3)
            nc=1
            ny=NYNP(1,nv,nh,np,0,nc,nr3)
            YP(ny,1,nx3)=YQ(NQNP(np),1,1,nx)
          ENDDO
        ELSE IF(CPTYPE.EQ.7) THEN
          DO nq=1,NPNODE(0,nr4)
            np=NPNODE(nq,nr4)
            nv=1
            nh=NH_LOC(1,nx4)
            nc=1
            ny=NYNP(1,nv,nh,np,0,nc,nr4)
            YP(ny,1,nx4)=YQ(NQNP(np),1,1,nx)
          ENDDO
        ENDIF
      ELSE
        CALL SOLVE_SYSTEM(ISC_GKK(1,nx),ISR_GKK(1,nx),NQT,NQT,NQT,
     '    NZZT(1,nr,nx),IWRIT4(nr,nx),PRECON_CODE(nx),SOLVEROPTION(nx),
     '    SPARSEGKK(nx),GKK(1,nx),GRR,YQ(1,1,1,nx),FIRST_A,
     '    UPDATE_MATRIX,X_INIT,nx,ERROR,*9999)
      ENDIF

C     Check solution vs analytic
      IF(CPTYPE.LT.6) THEN
        COUNT=NQT
        TOT_ABS_ERR=0.0d0
        TOT_REL_ERR=0.0d0
        WRITE(OP_STRING,'(/'' Analytic solution is u = x^2 - y^2'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        DO nq=1,NQT
          YQ(nq,2,1,nx)=
     '      DABS(YQ(nq,1,1,nx)-((XQ(1,nq)**2.0d0)-(XQ(2,nq)**2.0d0)))
          TOT_ABS_ERR=TOT_ABS_ERR+YQ(nq,2,1,nx)
          IF(DABS((XQ(1,nq)**2.0d0)-(XQ(2,nq)**2.0d0)).GT.LOOSE_TOL)
     '      THEN
            YQ(nq,3,1,nx)=
     '        DABS((YQ(nq,1,1,nx)-((XQ(1,nq)**2.0d0)-
     '        (XQ(2,nq)**2.0d0)))/((XQ(1,nq)**2.0d0)-
     '        (XQ(2,nq)**2.0d0)))
            TOT_REL_ERR=TOT_REL_ERR+YQ(nq,3,1,nx)
          ELSE
            COUNT=COUNT-1
            YQ(nq,3,1,nx)=0.0d0
          ENDIF
        ENDDO

        AVE_ABS_ERR=TOT_ABS_ERR/DBLE(NQT)
        AVE_REL_ERR=TOT_REL_ERR/DBLE(COUNT)

        WRITE(OP_STRING,'(/'' Average absolute error '',E12.6)')
     '    AVE_ABS_ERR
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Average relative error '',E12.6)')
     '    AVE_REL_ERR
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      IF(CPTYPE.EQ.2) THEN
        nv=1
        nh=NH_LOC(1,nx2)
        nc=1

        !Write potential back into YP
        DO nq=1,NPNODE(0,nr2)
          np=NPNODE(nq,nr2)
          ny=NYNP(1,nv,nh,np,0,nc,nr2)
          YP(ny,1,nx2)=YQ(NQNP(np),1,1,nx)
        ENDDO
      ENDIF
      IF(CPTYPE.EQ.4) THEN
        nv=1
        nh=NH_LOC(1,nx2)
        nc=1

        !Write potential back into YP
        DO nq=1,NPNODE(0,nr2)
          np=NPNODE(nq,nr2)
          ny=NYNP(1,nv,nh,np,0,nc,nr2)
          YP(ny,1,nx2)=YQ(NQNP(np),1,1,nx)
        ENDDO

        nv=1
        nh=NH_LOC(1,nx3)
        nc=1

        !Write potential back into YP
        DO nq=1,NPNODE(0,nr3)
          np=NPNODE(nq,nr3)
          ny=NYNP(1,nv,nh,np,0,nc,nr3)
          YP(ny,1,nx3)=YQ(NQNP(np),1,1,nx)
        ENDDO
      ENDIF

      CALL EXITS('SOLVE8_COUPLED')
      RETURN
 9999 CALL ERRORS('SOLVE8_COUPLED',ERROR)
      CALL EXITS('SOLVE8_COUPLED')
      RETURN 1
      END


C      SUBROUTINE SOLVE8_COUPLED(ISC_GKK,ISR_GKK,NEELEM,NENQ,NLQ,NQGP,
C     '  NQGP_PIVOT,NQS,NQXI,NRLIST,NWQ,NXLIST,NXQ,CQ,DNUDXQ,
C     '  DXDXIQ,DXDXIQ2,GCHQ,GK,GK2,GKK,GQ,GQ2,GUQ,NQGW,PROPQ,FIXQ,
C     '  ERROR,*)
C
CC#### Subroutine: SOLVE8_COUPLED
CC###  Description:
CC***  Created by Martin Buist, April 2000
C
C      IMPLICIT NONE
C
C      INCLUDE 'cmiss$reference:b01.cmn'
C      INCLUDE 'cmiss$reference:cbdi02.cmn'
C      INCLUDE 'cmiss$reference:geom00.cmn'
C      INCLUDE 'cmiss$reference:grid00.cmn'
C      INCLUDE 'cmiss$reference:loc00.cmn'
C      INCLUDE 'cmiss$reference:loc00.inc'
C      INCLUDE 'cmiss$reference:solv00.cmn'
C
C!     Parameter list
C      INTEGER ISC_GKK(NISC_GKKM,NXM),ISR_GKK(NISR_GKKM,NXM),
C     '  NEELEM(0:NE_R_M,0:NRM),NENQ(0:8,NQM),NLQ(NQM),
C     '  NQGP(0:22,NQM),NQGP_PIVOT(22,NQM),NQS(NEQM),NQXI(0:NIM,NQSCM),
C     '  NRLIST(0:NRM),NWQ(8,0:NQM),NXLIST(0:NXM),NXQ(-NIM:NIM,0:4,0:NQM)
C      REAL*8 CQ(NMM,NQM),DNUDXQ(3,3,NQM),DXDXIQ(3,3,NQM),
C     '  DXDXIQ2(3,3,NQM),GCHQ(3,NQM),GK(NZ_GK_M),GK2(*),
C     '  GKK(NZ_GKK_M,NXM),GQ(NZ_GQ_M),GQ2(*),GUQ(3,3,NQM),
C     '  NQGW(22,NQM),PROPQ(3,3,4,2,NQM)
C      CHARACTER ERROR*(*)
C      LOGICAL FIXQ(NYQM,NIYFIXM,NXM)
C!     Local variables
C      INTEGER BCINNBND_NQ,BCOUTBND_NQ,DUMMY_LIST(0:1),INNBOUND_NQ,
C     '  INTERNAL_NQ,NITB,nq,nqq,nq2,
C     '  nr,nrb1,nrb2,NRTEMP(0:1),ntemp,NUM_BE_INT,NUM_BE_EXT,NUM_BE_T,
C     '  nx,nxb1,nxb2,nxc,nzero,nzz,
C     '  OUTBOUND_NQ,SYSTEMSIZE
C      REAL*8 COEFFSEXT(19)
C
C      CALL ENTERS('SOLVE8_COUPLED',*9999)
C
C
CCCCCCCCCCCCCCCCCCC Note: need NQNP so need to have UP_NQNP false!
CCCCCCCCCCCCCCCCCCC Note: need to make NQGW filled for boundary gp's!
C
CC     Get 3 regions
C      CALL ASSERT(NRLIST(0).EQ.3,'>> Must specify 3 classes',
C     '  ERROR,*9999)
C      nr=NRLIST(1)
C      nrb1=NRLIST(2)
C      nrb2=NRLIST(3)
C      NITB=NQXI(0,NQS(NEELEM(1,nr)))
C
CC     Get 3 classes
C      CALL ASSERT(NXLIST(0).EQ.3,'>> Must specify 3 classes',
C     '  ERROR,*9999)
C      nxc=NXLIST(1)
C      CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
C      CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
C     '  ERROR,*9999)
C      nxc=NXLIST(2)
C      CALL NX_LOC(NX_INQUIRE,nxc,nxb1,NX_SOLVE,ERROR,*9999)
C      CALL ASSERT(nxb1.GT.0,'>>No nx defined for this solve class',
C     '  ERROR,*9999)
C      nxc=NXLIST(3)
C      CALL NX_LOC(NX_INQUIRE,nxc,nxb2,NX_SOLVE,ERROR,*9999)
C      CALL ASSERT(nxb2.GT.0,'>>No nx defined for this solve class',
C     '  ERROR,*9999)
C
CC     Calc grid point quadratics
C      NRTEMP(0)=1
C      NRTEMP(1)=NRLIST(1)
C      CALL GET_FD_POINTS(NEELEM,NENQ,NLQ,NQGP,NQGP_PIVOT,NQS,NQXI,
C     '  NRTEMP,NWQ,NXQ,ERROR,*9999)
C
CC     Calc system sizes
C      DO nq=1,NQT
C        NWQ(3,nq)=0
C      ENDDO !nq
C      DO nq=1,NQT
C        IF(NWQ(1,nq).NE.0) THEN
C          IF(NWQ(1,nq).GT.nq) THEN
C            DO nqq=1,NQGP(0,nq)
C              nq2=NQGP(nqq,nq)
C              NWQ(3,nq2)=2
C            ENDDO !nqq
C            NWQ(3,nq)=1
C          ELSE
C            DO nqq=1,NQGP(0,nq)
C              nq2=NQGP(nqq,nq)
C              NWQ(3,nq2)=3
C            ENDDO !nqq
C            NWQ(3,nq)=4
C          ENDIF
C        ENDIF
C      ENDDO !nq
C      DO nq=1,NQT
C        IF(NWQ(1,nq).NE.0) THEN
C          IF(NWQ(1,nq).GT.nq) THEN
C           NWQ(3,nq)=1
C          ELSE
C            NWQ(3,nq)=4
C          ENDIF
C        ENDIF
C      ENDDO !nq
C
CC     Calculate grid array sizes
C      INTERNAL_NQ=0
C      INNBOUND_NQ=0
C      BCINNBND_NQ=0
C      BCOUTBND_NQ=0
C      OUTBOUND_NQ=0
C      DO nq=1,NQT
C        IF(NWQ(3,nq).EQ.0) THEN
C          INTERNAL_NQ=INTERNAL_NQ+1
C        ELSE IF(NWQ(3,nq).EQ.1) THEN
C          INNBOUND_NQ=INNBOUND_NQ+1
C        ELSE IF(NWQ(3,nq).EQ.2) THEN
C          BCINNBND_NQ=BCINNBND_NQ+1
C        ELSE IF(NWQ(3,nq).EQ.3) THEN
C          BCOUTBND_NQ=BCOUTBND_NQ+1
C        ELSE IF(NWQ(3,nq).EQ.4) THEN
C          OUTBOUND_NQ=OUTBOUND_NQ+1
C        ENDIF
C      ENDDO !nq
C      NUM_BE_INT=NOT(1,1,nrb1,nxb1)
C      NUM_BE_T=NOT(1,1,nrb2,nxb2)-OUTBOUND_NQ
C      NUM_BE_EXT=NOT(1,1,nrb2,nxb2)-NUM_BE_T
C      CALL ASSERT(NUM_BE_INT.EQ.INNBOUND_NQ,
C     '  '>>Number of boundary grid points and BE nodes differ',
C     '  ERROR,*9999)
C      CALL ASSERT(NUM_BE_EXT.EQ.OUTBOUND_NQ,
C     '  '>>Number of boundary grid points and BE nodes differ',
C     '  ERROR,*9999)
C      SYSTEMSIZE=NQT+NUM_BE_T
C
CC     Diagnostics
C      IF(DOP) THEN
C        WRITE(OP_STRING,'('' Total nq            = '',I8)') NQT
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' Internal nq         = '',I8)') INTERNAL_NQ
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' Inner Boundary nq   = '',I8)') INNBOUND_NQ
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' Inner Bound conn nq = '',I8)') BCINNBND_NQ
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' Outer Boundary nq   = '',I8)') OUTBOUND_NQ
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' Outer Bound conn nq = '',I8)') BCOUTBND_NQ
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'(''  '')')
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' Total BE            = '',I8)')
C     '    NOT(1,1,nrb1,nxb1)+NOT(1,1,nrb2,nxb2)
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' Epicardial BE       = '',I8)') NUM_BE_EXT
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' Endocardial BE      = '',I8)') NUM_BE_INT
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' Torso BE            = '',I8)') NUM_BE_T
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'(''  '')')
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' Matrix system size  = '',I8)') SYSTEMSIZE
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C      ENDIF
C
CC     Check array sizes
C      NOT(1,1,nr,nx)=SYSTEMSIZE
C      NOT(2,1,nr,nx)=SYSTEMSIZE
C      nzz=SYSTEMSIZE*SYSTEMSIZE
C      NZZT(1,nr,nx)=nzz
C      IF(NZ_GKK_M.LT.nzz) THEN
C        WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
C        GOTO 9999
C      ENDIF
C
CC     Initialising
C      DO nq=1,nzz
C        GKK(nq,nx)=0.0d0
C      ENDDO !nq
C
CC     Assemble grid coeffs
CC      ntemp=NX_LIST(0)
CC      NX_LIST(0)=1
CC      DO nq=NQR(1,nr),NQR(2,nr)
CC        IF(NWQ(1,nq).EQ.0) THEN
CC          CALL CALC_GRID_COEF(NENQ,NITB,nq,NQGP,NQS,NQXI,NWQ(1,nq),
CC     '      NXQ,nx,nx,COEFFSEXT,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,
CC     '      GCHQ(1,nq),GUQ(1,1,nq),NQGW(1,nq),PROPQ(1,1,1,1,nq),.TRUE.,
CC     '      FIXQ,.FALSE.,.TRUE.,ERROR,*9999)
CC          DO nzero=1,NQGP(0,nq)
CC            GKK(((NQGP(nzero,nq)-1)*SYSTEMSIZE)+nq,nx)=
CC     '        COEFFSEXT(NQGP_PIVOT(nzero,nq))
CC          ENDDO !nz
CC        ENDIF
CC      ENDDO !nq
CC      NX_LIST(0)=ntemp
C
CC     Boundary element coeffs
C
C
CC     Diagnostics
C      IF(DOP) THEN
C        WRITE(OP_STRING,
C     '    '(/'' Global stiffness matrix GKK - external:'')')
C        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' NOT(1,1,nr,nx)='',I5,'
C     '    //''', NOT(2,1,nr,nx)='',I5)') NOT(1,1,nr,nx),
C     '    NOT(2,1,nr,nx)
C        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C        CALL OPSTFMAT(DUMMY_LIST,ISC_GKK(1,nx),ISR_GKK(1,nx),IOOP,
C     '    NOT(1,1,nr,nx),NOT(2,1,nr,nx),
C     '    NZZT(1,nr,nx),DUMMY_LIST,SPARSEGKK(nx),
C     '    GKK(1,nx),GKK(1,nx),'GKK','GKK',.TRUE.,.TRUE.,.FALSE.,
C     '    ERROR,*9999)
C      ENDIF
C
C      CALL EXITS('SOLVE8_COUPLED')
C      RETURN
C 9999 CALL ERRORS('SOLVE8_COUPLED',ERROR)
C      CALL EXITS('SOLVE8_COUPLED')
C      RETURN 1
C      END


