      SUBROUTINE GEN_GRID_POTE_RHS2(IBT,IDO,INP,NBH,NBJ,NEELEM,NENQ,NHE,
     '  NHP,niqV,NKHE,NKH,NPF,NPNE,NPNODE,NQGP,NQGP_PIVOT,NQS,NQXI,
     '  nr_grid,NRLIST,NVHE,NVHP,NW,NWQ,nx_ext,nx_trans,
     '  NX_TORSO_LIST,NXQ,NYNE,NYNP,ALPHA,AQ,CQ,CURVCORRECT,
     '  DNUDXQ,DXDXIQ,DXDXIQ2,
     '  NQGW,RHS,SE,XIQ,YP,YQ,ZA,ZE,ZP,FIXQ,ERROR,*)

C#### Subroutine: GEN_GRID_POTE_RHS2
C###  Description:
C###    GEN_GRID_POTE_RHS2 is used for coupled bidomain/torso problems
C###    to generate a RHS vector for the extracellular potential
C###    matrix.
C***  Created by Martin Buist, October 1998

      IMPLICIT NONE

C      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'time02.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter list
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NENQ(0:8,NQM),NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),niqV,
     '  NKHE(NKM,NNM,NHM,NEM),NKH(NHM,NPM,NCM,0:NRM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),nr_grid,
     '  NRLIST(0:NRM),NQGP(0:NQGM,NQM),NQGP_PIVOT(NQGM,NQM),NQS(NEQM),
     '  NQXI(0:NIM,NQSCM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),NW(NEM,3),NWQ(8,0:NQM),
     '  nx_ext,nx_trans,NX_TORSO_LIST(9),NXQ(-NIM:NIM,0:4,0:NQM)
      REAL*8 ALPHA,AQ(NMAQM,NQM),CQ(NMM,NQM),CURVCORRECT(2,2,NNM,NEM),
     &  DNUDXQ(3,3,NQM),DXDXIQ(3,3,NQM),DXDXIQ2(3,3,NQM),
     &  RHS(NQM),NQGW(NQGM,NQM),SE(NSM,NBFM,NEM),XIQ(NIM,NQM),
     &  YP(NYM,NIYM,NXM),YQ(NYQM,NIQM,NAM,NXM),
     &  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
      LOGICAL FIXQ(NYQM,NIYFIXM,NXM)
!     Local variables
      INTEGER i,ne,ne1,ne2,nh,ni,nq,NITB,nc,nb,nr,nr1,nr2,
     '  nx,nhx
      REAL*8 CMAMDT,FLUX,ONEMINUSALPHA,PXI,XI(3)
      REAL ELAPSED_TIME,TIME_START(1),TIME_STOP(1)
      LOGICAL ERROR_FLAG

      CALL ENTERS('GEN_GRID_POTE_RHS2',*9999)

      CALL ASSERT(.NOT.UP_NENQ,' >>Update element grid first',
     '  ERROR,*9999)

      CALL ASSERT(ALPHA.GT.-ZERO_TOL,'>>Alpha must be [0-1]',
     '  ERROR,*9999)
      ONEMINUSALPHA=1.0d0-ALPHA
      CALL ASSERT(ONEMINUSALPHA.GT.-ZERO_TOL,'>>Alpha must be [0-1]',
     '  ERROR,*9999)

      !set the grid boundary conditions
      ERROR_FLAG=.FALSE.
      CALL CPU_TIMER(CPU_USER,TIME_START)

C$OMP PARALLEL DO
C$OMP&PRIVATE(CMAMDT,FLUX,i,nb,nc,ne,ne1,ne2,nh,nhx,ni,NITB,nq,nr,
C$OMP&  nr1,nr2,nx,XI,ZA,ZE,ZP)
C$OMP&SHARED(ALPHA,CQ,CURVCORRECT,DNUDXQ,DXDXIQ,DXDXIQ2,
C$OMP&  FIXQ,IBT,IDO,INP,NBH,
C$OMP&  NBJ,NEELEM,NENQ,NHE,NHP,NKHE,NKH,niqV,NPF,NPNE,NPNODE,NQGP,
C$OMP&  NQGP_PIVOT,NQGW,NQS,NQXI,nr_grid,NRLIST,NVHE,NVHP,
C$OMP&  NW,NWQ,NXQ,nx_ext,nx_trans,NX_TORSO_LIST,NYNE,NYNP,
C$OMP&  ONEMINUSALPHA,RHS,SE,XIQ,YP,YQ)
      DO nq=NQR(1,nr_grid),NQR(2,nr_grid)
        IF(.NOT.ERROR_FLAG) THEN
          IF(NWQ(1,nq).NE.0) THEN     !boundary grid point

            !initialise from previous info
            ne=NENQ(NENQ(0,nq),nq)
            IF(ne.GT.0) THEN
              nb=NBJ(1,ne)
              NITB=NIT(nb)
              DO ni=1,NITB
                XI(ni)=XIQ(ni,nq)
              ENDDO
            ENDIF

            DO nr1=2,NRLIST(0)
              nr2=NRLIST(nr1)
              DO ne1=1,NEELEM(0,nr2)
                ne2=NEELEM(ne1,nr2)
                IF(ne2.EQ.ne) THEN
                  nr=nr2
                  nx=NX_TORSO_LIST(nr1)
                ENDIF
              ENDDO !ne
            ENDDO !nr

            IF((nx.GT.0).AND.(nr.GT.0).AND.(ne.GT.0)) THEN
              IF(FIXQ(nq,1,nx_ext)) nc=1
              IF(FIXQ(nq,2,nx_ext)) nc=2

              nhx=1
              nh=NH_LOC(nhx,nx)
              nb=NBH(nh,nc,ne)
              IF(nb.EQ.0) nb=NBJ(1,ne)

              !calculate the b.c. value
              CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '          NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,
     '          YP(1,1,nx),ZA,ZP,ERROR,*200)

              CALL ZPZE(NBH(1,1,ne),nc,NHE(ne,nx),NKHE(1,1,1,ne),
     '          NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1),nx,
     '          CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,
     '          ERROR,*200)

              !use PXI to interpolate nodal values
              RHS(nq)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '          XI,ZE(1,nhx))
              IF(nc.EQ.2) RHS(nq)=-RHS(nq)

              IF(nc.EQ.1) THEN
                RHS(nq)=(ALPHA*YQ(nq,niqV,1,nx_ext))+
     '            (ONEMINUSALPHA*RHS(nq))
              ELSE IF(nc.EQ.2) THEN
                CALL GGRADPHIQDN(NENQ,niqV,nq,NQS,NQXI,NXQ,AQ,
     '            CQ(6,nq),DNUDXQ,DXDXIQ,DXDXIQ2,FLUX,
     '            YQ(1,1,1,nx_ext),ERROR,*200)
                RHS(nq)=(ALPHA*FLUX)+(ONEMINUSALPHA*RHS(nq))
C                RHS(nq)=RHS(nq)/1000.0d0 !S -> mS
              ENDIF

            ELSE !not in iteration
C MLB change for march8_coupled 7/4/00 to flux bc.
C              RHS(nq)=YQ(nq,niqV,1,nx_ext)
              RHS(nq)=0.0d0
            ENDIF

          ELSE !internal point
            CMAMDT=1.0d0
            RHS(nq)=0.0d0
            DO i=1,NQGP(0,nq)
              IF(NQGP(i,nq).GT.0) RHS(nq)=RHS(nq)-
     '          (NQGW(NQGP_PIVOT(i,nq),nq)*
     '          YQ(NQGP(i,nq),niqV,1,nx_trans)*CMAMDT)
            ENDDO !i
          ENDIF !external
          GO TO 202
          !this statement is designed to be skipped if no error
          !occurs. However if a error occurs within a subroutine
          !the alternate return jumps to line 100 to set the flag
 200      CONTINUE
C$OMP CRITICAL(GEN_GRID_POTE_RHS2_1)
          ERROR_FLAG=.TRUE.
          WRITE(OP_STRING,'(/'' >>ERROR: An error occurred - '
     '      //'results may be unreliable'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*201)
 201      CONTINUE
C$OMP END CRITICAL(GEN_GRID_POTE_RHS2_1)
 202      CONTINUE
        ENDIF !not error_flag
      ENDDO
C$OMP END PARALLEL DO

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      IF(IWRIT5(nr_grid,nx_trans).EQ.9) THEN
        ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(GEN_GRID_POTE_RHS2_2)
        WRITE(OP_STRING,'(1X,''Time to create RHS'',F6.2,''s cpu'')')
     '    ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(GEN_GRID_POTE_RHS2_2)
      ENDIF

      CALL EXITS('GEN_GRID_POTE_RHS2')
      RETURN
 9999 CALL ERRORS('GEN_GRID_POTE_RHS2',ERROR)
      CALL EXITS('GEN_GRID_POTE_RHS2')
      RETURN 1
      END
