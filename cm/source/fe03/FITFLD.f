      SUBROUTINE FITFLD(IBT,IDO,INP,ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,
     '  ISR_GKK,ISR_GQ,LGE,LN,NAN,NBH,NBHF,NBJ,NBJF,NDDL,NDLT,NEELEM,
     '  NENP,NENQ,NEP,NFF,NFFACE,NGAP,NHE,NHP,
     '  NKB,NKEF,NKH,NKHE,NKJE,NLL,
     '  NMNO,NNB,NNF,NNL,NONY,NPF,NP_INTERFACE,NPL,NPLIST,NPLIST2,NPNE,
     '  NPNF,NPNODE,NPNY,NQNE,NQS,NQXI,nr,NRE,NSB,NVHE,NVHP,
     '  NVJE,NVJF,NW,NWP,nxFIT,nxSOLVE,NXI,NYNE,
     '  NYNO,NYNP,NYNR,NYNY,Z_CONT_LIST,CE,CG,CGE,CONY,CP,
     '  CURVCORRECT,CYNO,CYNY,EDD,ER,ES,FEXT,GK,GKK,GR,GRR,GQ,PG,RE1,RG,
     '  SE,SF,SP,WD,WDL,WG,WU,XA,XE,XG,XID,
     '  XIDL,XIG,XIP,XO,XP,YG,YGF,YP,YQS,ZA,
     '  ZA1,Z_CONT,ZD,ZDL,ZE,ZE1,ZP,ZP1,FIRST_A,FIX,IN_PLANE,ERROR,*)

C#### Subroutine: FITFLD
C###  Description:
C###    FITFLD fits field variables to either data or gauss points.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'fit000.cmn'
C      INCLUDE 'cmiss$reference:fit001.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'grou00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'solv00.cmn'
      INCLUDE 'time02.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISC_GK(NISC_GKM),ISC_GKK(NISC_GKKM),ISC_GQ(NISC_GQM),
     '  ISR_GK(NISR_GKM),ISR_GKK(NISR_GKKM),ISR_GQ(NISR_GQM),
     '  LGE(NHM*NSM,NRCM),LN(0:NEM),NAN(NIM,NAM,NBFM),
     '  NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NDDL(NEM,NDEM),NDLT(NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),NENQ(0:8,NQM),
     '  NEP(NPM),NFF(6,NEM),NFFACE(0:NF_R_M,NRM),
     '  NGAP(NIM,NBM),NHE(NEM),NHP(NPM),
     '  NKB(2,2,2,NNM,NBFM),NKEF(0:4,16,6,NBFM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NLL(12,NEM),NMNO(1:2,0:NOPM),
     '  NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),
     '  NNL(0:4,12,NBFM),NONY(0:NOYM,NYM,NRCM,0:NRM),
     '  NPF(9,NFM),NP_INTERFACE(0:NPM,0:3),NPL(5,0:3,NLM),
     '  NPLIST(0:NPM),NPLIST2(0:NPM),
     '  NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),NQNE(NEQM,NQEM),NQS(NEQM),
     '  NQXI(0:NIM,NQSCM),nr,NRE(NEM),NSB(NKM,NNM,NBFM),
     '  NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM),NVJE(NNM,NBFM,NJM,NEM),NVJF(NNM,NBFM,NJM),
     '  NWP(NPM,2),NW(NEM,3),nxFIT,nxSOLVE,NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),NYNO(0:NYOM,NOOPM,NRCM,0:NRM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),NYNY(0:NYYM,NYM,NRM),
     '  Z_CONT_LIST(NDM,2,7)
      REAL*8 CE(NMM,NEM,NXM),CG(NMM,NGM),CGE(NNM,NGM,NEM,NXM),
     '  CONY(0:NOYM,NYM,NRCM,0:NRM),CYNY(0:NYYM,NYM,NRM),
     '  CP(NMM,NPM,NXM),CURVCORRECT(2,2,NNM,NEM),
     '  CYNO(0:NYOM,NOOPM,NRCM,0:NRM),
     '  EDD(NDM),ER(NHM*NSM),ES(NHM*NSM,NHM*NSM),
     '  FEXT(NIFEXTM,NGM,NEM),GK(NZ_GK_M),
     '  GKK(NZ_GKK_M),GR(NYROWM),GRR(NOM),GQ(NZ_GQ_M),
     '  PG(NSM,NUM,NGM,NBM),RE1(NSM,NHM),RG(NGM),
     '  SE(NSM,NBFM,NEM),SF(NSM,NBFM),
     '  WD(NJM,NDM),WDL(NHM,NDEM),
     '  WG(NGM,NBM),WU(0:NUM+1,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),XID(NIM,NDM),
     '  XIDL(NIM,NDEM),XIG(NIM,NGM,NBM),XIP(NIM,NPM),XO(NOM),
     '  XP(NKM,NVM,NJM,NPM),
     '  YG(NIYGM,NGM,NEM),YGF(NIYGFM,NGFM,NFM),YP(NYM,NIYM,NXM),
     '  YQS(NIQSM,NQM),ZA(NAM,NHM,NCM,NEM),ZA1(NAM,NHM,NCM,NEM),
     '  Z_CONT(NDM,2,67),
     '  ZD(NJM,NDM),ZDL(NHM,NDEM),
     '  ZE(NSM,NHM),ZE1(NSM,NHM),SP(NKM,NBFM,NPM),
     '  ZP(NKM,NVM,NHM,NPM,NCM),ZP1(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
      LOGICAL FIRST_A,FIX(NYM,NIYFIXM),IN_PLANE
!     Local Variables
      INTEGER elem2,GETNYR,IEND,l,LGE_FACE(48),
     '  matrix_njj,matrix_nh,matrix_ny,nb,nd,ndl,nef,nn,no_resid_set,
     '  NTOT,ne,nf,ng,nh,nhj,nhs1,nhs2,NHST(2),NHST_FACE(NRCM),nhx,
     '  nj,njj,nk,NKHF(NKM,NNM,NHM),
     '  NKJF(NKM,NNM,NJM),nm,no,no1,no2,NOTT,nonode,nogroup,no_nynr,
     '  no_nynr1,no_nynr2,no_times,noy1,noy2,np,npe,ns1,ns2,ns3,ns4,
     '  NUMRESID,nv,NVHF(NNM,NBFM,NHM),ny,ny1,ny2,ny3,nyo1,nz,
     '  nzz,vorogrno
      INTEGER*4 WORK_PTR
      REAL AVETIME,ELAPSED_TIME,TIME_START1(1),TIME_START2(1),
     '  TIME_START3(1),TIME_START4(1),TIME_STOP(1)
      REAL*8 co1,co2,DELTA,EDG,PXI,RESIDUAL,
     '  SAED,SMED,SQED,SUM,VALUE,Z(6)
      REAL*8  Esmooth1(4,4),Esmooth2(4,4),Esmooth3(4,4)
C KAT 14Jan00: using GQ (instead of E) for dynamic memory
C     '  E(1500,1500) !temporary until allocate memory (PJH)
      LOGICAL UPDATE_MATRIX
      PARAMETER (DELTA=1.0D-4)
C      INTEGER elem
      DATA Esmooth1/2.d0,-2.d0, 1.d0,-1.d0,
     '             -2.d0, 2.d0,-1.d0, 1.d0,
     '              1.d0,-1.d0, 2.d0,-2.d0,
     '             -1.d0, 1.d0,-2.d0, 2.d0/
      DATA Esmooth2/2.d0, 1.d0,-2.d0,-1.d0,
     '              1.d0, 2.d0,-1.d0,-2.d0,
     '             -2.d0,-1.d0, 2.d0, 1.d0,
     '             -1.d0,-2.d0, 1.d0, 2.d0/
      DATA Esmooth3/1.d0,-1.d0,-1.d0, 1.d0,
     '             -1.d0, 1.d0, 1.d0,-1.d0,
     '             -1.d0, 1.d0, 1.d0,-1.d0,
     '              1.d0,-1.d0,-1.d0, 1.d0/

      CALL ENTERS('FITFLD',*9999)

      CALL CPU_TIMER(CPU_USER,TIME_START1)

       write(*,*) 'entering fitfld' 

C KAT 14Jan00:
      IF(KTYP8.EQ.7) THEN !material fit
        NUMRESID=NYNR(0,1,1,nr,nxSOLVE)
        NOTT=NUM_FIT(1)*NPNODE(0,nr)
        IF(NUMRESID*NOTT.GT.NZ_GQ_M) THEN
          IEND=0
          CALL APPENDC(IEND,' >>Increase NZ_GQ_M to ',ERROR)
          CALL APPENDI(IEND,NUMRESID*NOTT,ERROR)
          GOTO 9999
        ENDIF
      ENDIF

      matrix_njj=0
      DO njj=1,NUM_FIT(0) !Loop of the number of fit problems

        CALL CPU_TIMER(CPU_USER,TIME_START2)

        WRITE(OP_STRING,'(/'' Fitting problem #'',I1,'':'')') njj
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

        CALL CPU_TIMER(CPU_USER,TIME_START3)

        IF(njj.EQ.1) THEN !first fit variable is being fitted
          UPDATE_MATRIX=.TRUE.
          IF(NUM_FIT(0).GT.1) FIRST_A=.TRUE. !may be new sparsity pattern
        ELSE
          UPDATE_MATRIX=.FALSE.
C***      Check if the GK matrix needs to be updated
C         Check number of variables
          IF(NUM_FIT(njj).NE.NUM_FIT(matrix_njj)) THEN
            UPDATE_MATRIX=.TRUE.
            FIRST_A=.TRUE. !matrix is also a different shape
          ELSE
            DO nhj=1,NUM_FIT(njj)
              nhx=NLH_FIT(nhj,3,njj)
              nh=NH_LOC(nhx,nxFIT)
              nhx=NLH_FIT(nhj,3,matrix_njj)
              matrix_nh=NH_LOC(nhx,nxFIT)
C             Check bases and versions
              l=0
              DO WHILE(l.LT.LN(0).AND..NOT.UPDATE_MATRIX)
                l=l+1
                ne=LN(l)
                nb=NBH(nh,1,ne)
                IF(nb.NE.NBH(matrix_nh,1,ne)) THEN
                  UPDATE_MATRIX=.TRUE.
                  FIRST_A=.TRUE. !matrix is also a different shape
                ELSE
                  DO nn=1,NNT(nb)
                    IF(NVHE(nn,nb,nh,ne).NE.NVHE(nn,nb,matrix_nh,ne))
     '                THEN
                      UPDATE_MATRIX=.TRUE.
                      FIRST_A=.TRUE. !matrix is also a different shape
                    ENDIF
                  ENDDO !nn
                ENDIF !nb
              ENDDO !l (ne)
C             Check dofs to see if sparsity pattern may have changed
              nonode=0
              DO WHILE(nonode.LT.NPNODE(0,nr).AND..NOT.FIRST_A)
                nonode=nonode+1
                np=NPNODE(nonode,nr)
                IF(NVHP(nh,np,1).NE.NVHP(nh,np,1).OR.
     '            NKH(nh,np,1,nr).NE.NKH(nh,np,1,nr)) THEN
                  FIRST_A=.TRUE.
                ELSE
                  DO nv=1,NVHP(nh,np,1)
                    DO nk=1,NKH(nh,np,1,nr)
                      ny=NYNP(nk,nv,nh,np,0,1,nr)
                      matrix_ny=NYNP(nk,nv,matrix_nh,np,0,1,nr)
                      IF(FIX(ny,1).NEQV.FIX(matrix_ny,1)) FIRST_A=.TRUE.
                    ENDDO !nk
                  ENDDO !nv
                ENDIF !NVHP,NKH
              ENDDO !nonode
            ENDDO !nhj
          ENDIF !num_fit
        ENDIF
C KAT 17Nov98: Quick hack.
C       GK should not need recalculating unless UPDATE_MATRIX is .TRUE.,
C       but the entries need to be copied to different ny locations for
C       the no-ny maps in the GKK assembly to work.
C       Until this is written GK is recalculated whenever GKK changes.
        UPDATE_MATRIX=FIRST_A
        IF(UPDATE_MATRIX) THEN
          matrix_njj=njj
        ENDIF

C       ..initialise global variables
C KAT 17Nov98: Whole matrix is initiallized for all njj so only do this
C       once.  (There are separate nys for each njj.)
        IF(njj.EQ.1) THEN
C        IF(UPDATE_MATRIX) THEN
          CALL INIT_SPARSE_MATRIX(NYNR(0,2,1,nr,nxFIT),
     '      NYNR(0,2,1,nr,nxFIT),ISC_GK,ISR_GK,0,0,1,NYT(1,1,nxFIT),
     '      NZ_GK_M,NZT(1,nxFIT),NYNR(0,1,1,nr,nxFIT),
     '      NYNR(0,1,1,nr,nxFIT),KTYP24,GK,ERROR,*9999)
C        ENDIF

          DO no_nynr1=1,NYNR(0,1,1,nr,nxFIT)
            ny1=NYNR(no_nynr1,1,1,nr,nxFIT)
            GR(ny1)=0.d0
          ENDDO !no_nynr1 (ny1)
        ENDIF !njj=1

        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START3(1)
        IF(IWRIT4(nr,nxFIT).GE.1) THEN
          WRITE(OP_STRING,'(/'' CPU time for global matrix setup and '
     '      //'initialisation: '',D11.4,'' s'')') ELAPSED_TIME
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

        NO_TIMES=0
        AVETIME=0.0
        CALL CPU_TIMER(CPU_USER,TIME_START3)


C------------------ begin material parameter fitting -------------------

        IF(KTYP8.EQ.7) THEN !calc global matrix directly (mat.l fits)

          DO no_resid_set=1,KTYP28 !loop over residual sets
            IF(IWRIT4(nr,nxFIT).GE.1.OR.DOP) THEN
              WRITE(OP_STRING,'(/'' Residual set#'',I2)') no_resid_set
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF !dop

C           ..put coords & resids for no_resid_set into YP(ny,1,nxSOLVE)
            DO no_nynr=1,NUMRESID      !loop over rows
              ny = NYNR(no_nynr,1,1,nr,nxSOLVE)      !is row number
              YP(ny ,1,nxSOLVE)=YP(ny ,10+no_resid_set,nxSOLVE)
              ny1=GETNYR(2,NPNY,nr,0,1,ny,NYNE,NYNP) !is corr force row#
              YP(ny1,1,nxSOLVE)=YP(ny1,10+no_resid_set,nxSOLVE)
            ENDDO !no_nynr

C           ..calc residuals at ref state (stored in YP(ny,4,nxSOLVE))
            CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH(1,1,1,nr),
     '        NPNODE,nr,NVHP,nxSOLVE,NYNE,NYNP,
     '        YP(1,1,nxSOLVE),ZA,ZP,ERROR,*9999)
            CALL ZPRP(IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,
     '        NEELEM,NFF,
     '        NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH(1,1,1,nr),NKHE,NKJE,
     '        NNB,NNF,NPF,NPNE,NPNODE,NPNY(0,1,0,nxSOLVE),nr,NRE,NSB,
     '        NVHE,NVHP,NVJE,NW,nxSOLVE,NXI,NYNE,NYNP,
     '        NYNR(0,0,1,nr,nxSOLVE),Z_CONT_LIST,CE(1,1,nxSOLVE),CG,
     '        CGE(1,1,1,nxSOLVE),CP(1,1,nxSOLVE),
     '        CURVCORRECT,FEXT,PG,
     '        RE1,RG,SE,WG,XA,XG,XP,YG,YGF,
     '        YP(1,1,nxSOLVE),ZA,ZA1,
     '        %VAL(0),Z_CONT,ZE,ZE1,ZP,ZP1,%VAL(0),FIX,ERROR,*9999)

C           ..temporarily store residuals in YP(ny,4,nxFIT)
            DO no_nynr=1,NUMRESID !loop over rows
              ny=NYNR(no_nynr,1,1,nr,nxSOLVE)   !is row number
              YP(ny,4,nxFIT)=YP(ny,4,nxSOLVE)
            ENDDO !no_nynr

            no=0 !initialise column variable
            DO nhj=1,NUM_FIT(1) !loop over material params
              nm=NMNO(1,nhj)
              IF(IWRIT4(nr,nxFIT).GE.1.OR.DOP) THEN
                WRITE(OP_STRING,'(/'' Material parameter#'',I2)') nm
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF !dop
C KAT 23Dec99: Added assert
              CALL ASSERT(ILP(nm,1,nr,nxSOLVE).EQ.3,
     '          'Material parameters must be piecewise linear',
     '          ERROR,*9999)
              DO nonode=1,NPNODE(0,nr) !loop over nodes
                np=NPNODE(nonode,nr)
                no=no+1 !increment column variable
                IF(DOP) THEN
                  WRITE(OP_STRING,'(/'' ..perturb matl param nm='',I3,'
     '              //''' at node '',I5,'': no='',I5)') nm,np,no
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF !dop

C               ..perturb material parameter by delta
                CP(nm,np,nxSOLVE)=CP(nm,np,nxSOLVE)+DELTA

C               ..calc residuals at perturbed state (stored in YP(ny,4))
                CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH(1,1,1,nr),NPNODE,
     '            nr,NVHP,nxSOLVE,NYNE,NYNP,
     '            YP(1,1,nxSOLVE),ZA,ZP,ERROR,*9999)
                CALL ZPRP(IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,
     '            NEELEM,NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,
     '            NKH(1,1,1,nr),NKHE,NKJE,NNB,NNF,NPF,NPNE,NPNODE,
     '            NPNY(0,1,0,nxSOLVE),nr,NRE,NSB,NVHE,NVHP,NVJE,
     '            NW,nxSOLVE,NXI,NYNE,NYNP,NYNR(0,0,1,nr,nxSOLVE),
     '            Z_CONT_LIST,
     '            CE(1,1,nxSOLVE),CG,CGE(1,1,1,nxSOLVE),
     '            CP(1,1,nxSOLVE),CURVCORRECT,
     '            FEXT,PG,RE1,RG,SE,WG,XA,XG,XP,YG,YGF,
     '            YP(1,1,nxSOLVE),
     '            ZA,ZA1,%VAL(0),Z_CONT,ZE,ZE1,ZP,ZP1,%VAL(0),FIX,
     '            ERROR,*9999)

C               ..reset material parameter
                CP(nm,np,nxSOLVE)=CP(nm,np,nxSOLVE) -DELTA

                DO no_nynr=1,NUMRESID !loop over resids
                  ny=NYNR(no_nynr,1,1,nr,nxSOLVE)   !is row number
C KAT 14Jan00: using GQ (instead of E) for dynamic memory
                  GQ(no_nynr+NUMRESID*(no-1))=
     '              (YP(ny,4,nxSOLVE)-YP(ny,4,nxFIT))/DELTA
C                  E(ny,no)=(YP(ny,4,nxSOLVE)-YP(ny,4,nxFIT))/DELTA
                  IF(DOP) THEN
C KAT 14Jan00: using GQ (instead of E) for dynamic memory
                    WRITE(OP_STRING,'('' ny='',I5,'' no='',I5,'
     '                //''' E(ny,no)='',E13.5)')
     '                ny,no,GQ(no_nynr+NUMRESID*(no-1))
C     '                //''' E(ny,no)='',E13.5)') ny,no,E(ny,no)
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF !dop
                ENDDO !no_nynr

              ENDDO !nonode
            ENDDO !nhj
C KAT 14Jan00: Need this earlier for checking array size
C            NOTT=no !temporary PJH

C           ..calc cpt of resids left when fitted matl params are zero
            DO nhj=1,NUM_FIT(1) !loop over material params
              nm=NMNO(1,nhj)
              DO nonode=1,NPNODE(0,nr) !loop over nodes
                np=NPNODE(nonode,nr)
                CP(nm,np,nxFIT)=CP(nm,np,nxSOLVE) !temporary storage
                CP(nm,np,nxSOLVE)=0.d0
              ENDDO !nonode
            ENDDO !nhj

            CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH(1,1,1,nr),
     '        NPNODE,nr,NVHP,nxSOLVE,NYNE,NYNP,
     '        YP(1,1,nxSOLVE),ZA,ZP,ERROR,*9999)
            CALL ZPRP(IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,
     '        NEELEM,NFF,
     '        NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH(1,1,1,nr),NKHE,NKJE,
     '        NNB,NNF,NPF,NPNE,NPNODE,NPNY(0,1,0,nxSOLVE),nr,NRE,NSB,
     '        NVHE,NVHP,NVJE,NW,nxSOLVE,NXI,NYNE,NYNP,
     '        NYNR(0,0,1,nr,nxSOLVE),Z_CONT_LIST,CE(1,1,nxSOLVE),CG,
     '        CGE(1,1,1,nxSOLVE),
     '        CP(1,1,nxSOLVE),CURVCORRECT,
     '        FEXT,PG,RE1,RG,SE,WG,XA,XG,XP,YG,
     '        YGF,YP(1,1,nxSOLVE),
     '        ZA,ZA1,%VAL(0),Z_CONT,ZE,ZE1,ZP,ZP1,%VAL(0),FIX,
     '        ERROR,*9999)

C           ..restore material parameters in nxSOLVE
            DO nhj=1,NUM_FIT(1) !loop over material params
              nm=NMNO(1,nhj)
              DO nonode=1,NPNODE(0,nr) !loop over nodes
                np=NPNODE(nonode,nr)
                CP(nm,np,nxSOLVE)=CP(nm,np,nxFIT)
              ENDDO !nonode
            ENDDO !nhj

C           ..temporarily store these residuals in YP(ny,4,nxFIT)
            DO no_nynr=1,NUMRESID !loop over rows
              ny=NYNR(no_nynr,1,1,nr,nxSOLVE)   !is row number
              YP(ny,4,nxFIT)=-YP(ny,4,nxSOLVE)
              IF(DOP) THEN
                WRITE(OP_STRING,'('' YP('',I5,'',4,nxFIT)='',E13.5)')
     '            ny,YP(ny,4,nxFIT)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF !dop
            ENDDO !no_nynr

            DO no1=1,NOTT !loop over material params

C             ..compute RHS vector
              SUM=0.d0
              DO no_nynr=1,NUMRESID !loop over resids
                ny=NYNR(no_nynr,1,1,nr,nxSOLVE)   !is row number
C KAT 14Jan00: using GQ (instead of E) for dynamic memory
                SUM=SUM+GQ(no_nynr+NUMRESID*(no1-1))*YP(ny,4,nxFIT)
C                SUM=SUM+E(ny,no1)*YP(ny,4,nxFIT)
              ENDDO !no_nynr
              GR(no1)=GR(no1)+SUM
              IF(DOP) THEN
                WRITE(OP_STRING,'(/'' GR('',I5,'')='',E13.5)')
     '            no1,GR(no1)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF !dop

C             ..compute stiffness matrix
              DO no2=1,NOTT !loop over material params
                SUM=0.d0
                DO no_nynr=1,NUMRESID !loop over resids
                  ny=NYNR(no_nynr,1,1,nr,nxSOLVE) !is row number
C KAT 14Jan00: using GQ (instead of E) for dynamic memory
                  SUM=SUM+GQ(no_nynr+NUMRESID*(no1-1))
     '              *GQ(no_nynr+NUMRESID*(no2-1))
C                  SUM=SUM+E(ny,no1)*E(ny,no2)
                ENDDO !no_nynr
                CALL SPARSE(no1,no2,NOTT,nz,NZ_GK_M,
     '            NZT(1,nxFIT),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
                GK(nz)=GK(nz)+SUM
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' no1='',I5,'' no2='',I5,'
     '              //''' nz='',I5,'' GK='',E13.5)') no1,no2,nz,GK(nz)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF !dop
              ENDDO !no2
            ENDDO !no1

          ENDDO !no_resid_set

C         .. add element smoothing terms into global system
          IF(IWRIT4(nr,nxFIT).GE.1.OR.DOP) THEN
            OP_STRING(1)=' Smoothing...'
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF !dop
          DO l=1,LN(0) !loop over elements in the fit
            ne=LN(l)
            CALL MELGEF(LGE,NBH(1,1,ne),ne,NHST,njj,NPNE(1,1,ne),nr,
     '        NVHE(1,1,1,ne),nxFIT,NYNE,NYNP,ERROR,*9999)
            DO nhs1=1,4
              no1=IABS(LGE(nhs1,1))
              DO nhs2=1,4
                no2=IABS(LGE(nhs2,2))
C             ..smoothing terms on derivs wrt Xi coords
                SUM=0.1666667d0*WU(2,ne)*Esmooth1(nhs1,nhs2) !1st Xi deriv
     '             +0.1666667d0*WU(3,ne)*Esmooth2(nhs1,nhs2) !2nd Xi deriv
     '                         +WU(4,ne)*Esmooth3(nhs1,nhs2) !cross deriv
                CALL SPARSE(no1,no2,NOTT,nz,NZ_GK_M,
     '            NZT(1,nxFIT),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
                GK(nz)=GK(nz)+SUM
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' nhs1='',I5,'' nhs2='',I5,'
     '              //''' SUM='',E13.5)') nhs1,nhs2,SUM
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'(''  no1='',I5,''  no2='',I5,'
     '              //''' nz='',I5,'' GK='',E13.5)') no1,no2,nz,GK(nz)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF !dop
              ENDDO !nhs2
            ENDDO !nhs1
          ENDDO !l (ne)

! May need this in future
!          DO no_nynr=1,NYNR(0,1,1,nr,nx_sol) !loop over rows
!            ny=NYNR(no_nynr,1,1,nr,nx_sol) !is row number
!            DO noy=1,NONY(0,ny,2,nr,nx_sol)
!              no=NONY(noy,ny,2,nr,nx_sol)
!              DO noopti=1,NMNO(1,0,nx_opt)
!                RESJAC(no,noopti)=D_RP(ny,noopti)
!              ENDDO !noopti
!            ENDDO !noy
!          ENDDO !no_nynr

C--------------------- all other fitting ------------------------------

        ELSE !calc global matrix from elem matrices (non-mat.l fitting)

          DO l=1,LN(0) !loop over elements/faces in the fit

            CALL CPU_TIMER(CPU_USER,TIME_START4)

            IF(KTYP11.EQ.1) THEN !element fitting
              ne=LN(l)
            ELSEIF(KTYP11.EQ.2) THEN !face fitting
              nf=LN(l)
              ne=NPF(6,nf)
              nef=NPF(8,nf)
              WRITE(OP_STRING,'(/'//
     &          ''' Element No / Local Face No / Global Face No : '','//
     &          '3(2x,I5))') ne,nef,nf
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF

            IF(KTYP11.EQ.2) THEN !face fitting
              CALL CALC_FACE_INFORMATION_IND(NBJ(1,ne),NBJF(1,nf),
     '          nef,NKJE(1,1,1,ne),NKEF,NKJF,NNF,NPNE(1,1,ne),NPNF,
     '          nr,NVJE(1,1,1,ne),NVJF,SE(1,1,ne),SF,ERROR,*9999)

              DO nhj=1,NUM_FIT(njj)
                nj=NLH_FIT(nhj,1,njj)
                nhx=NLH_FIT(nhj,3,njj)
                nh=NH_LOC(nhx,nxFIT)
                nb=NBJF(nj,nf)
                NBHF(nh,1,nf)=nb
                DO nn=1,NNT(nb)
                  NVHF(nn,nb,nh)=NVJF(nn,nb,nj)
                  DO nk=1,NKT(nn,nb)
                    NKHF(nk,nn,nh)=NKJF(nk,nn,nj)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF

            IF(KTYP11.EQ.1) THEN !element fitting
              CALL MELGEF(LGE,NBH(1,1,ne),ne,NHST,njj,NPNE(1,1,ne),nr,
     '          NVHE(1,1,1,ne),nxFIT,NYNE,NYNP,ERROR,*9999)
            ELSEIF(KTYP11.EQ.2) THEN !face fitting
              CALL MELGEF_FACE(LGE_FACE,NBHF(1,1,nf),
     '          NHST_FACE,njj,NKH,NKJF,NPNF,NPNODE,nr,NVHF,NVHP,
     '          nxFIT,NYNP,ERROR,*9999)
            ENDIF

            IF(IWRIT4(nr,nxFIT).GE.5) THEN
              WRITE(OP_STRING,
     &          '(/'' Element'',I5,'', Number of variables: '','//
     &          '''NHST(1)='',I3,'', NHST(2)='',I3)') ne,NHST(1),NHST(2)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' LGE(1..,1): '',14I5,:(/13X,14I5))')
     &          (LGE(nhs1,1),nhs1=1,NHST(1))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' LGE(1..,2): '',14I5,:(/13X,14I5))')
     &          (LGE(nhs1,2),nhs1=1,NHST(2))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF

            IF(KTYP11.EQ.1) THEN !element fitting
              DO nhs1=1,NHST(1)
                ER(nhs1)=0.0d0
                IF(UPDATE_MATRIX) THEN
                  DO nhs2=1,NHST(2)
                    ES(nhs1,nhs2)=0.0d0
                  ENDDO !nhs2
                ENDIF
              ENDDO !nhs1
            ELSEIF(KTYP11.EQ.2) THEN !face fitting
              DO nhs1=1,NHST_FACE(1)
                ER(nhs1)=0.0d0
                IF(UPDATE_MATRIX) THEN
                  DO nhs2=1,NHST_FACE(2)
                    ES(nhs1,nhs2)=0.0d0
                  ENDDO !nhs2
                ENDIF
              ENDDO !nhs1
            ENDIF

            IF(KTYP11.EQ.1) THEN !element fitting
              CALL ZPZE(NBH(1,1,ne),1,NH_LOC(0,nxFIT),NKHE(1,1,1,ne),
     '          NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1),nxFIT,
     '          CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,
     '          ERROR,*9999)
            ELSEIF(KTYP11.EQ.2) THEN !face fitting
              CALL ZPZE_FACE(nf,NBHF,1,NKHF,NPNF,NVHF,nxFIT,SF,ZE,
     '          ZP,ERROR,*9999)
            ENDIF

            IF(KTYP6.NE.1.AND.KTYP8.EQ.1.AND.ITYP6(nr,nxFIT).EQ.1) THEN
C             linear geometric fit
              IF(ITYP10(nr).EQ.2) THEN !cylindrical coords
                DO nhj=1,NUM_FIT(njj)
                  IF(nhj.EQ.2) THEN !theta
C                   correct the initial element field pos theta angles
                    nhx=NLH_FIT(nhj,3,njj)
                    nh=NH_LOC(nhx,nxFIT)
                    nb=NBH(nh,1,ne)
                    ns1=1
                    ns2=1+NKT(1,nb)
                    ns3=1+NKT(1,nb)+NKT(2,nb)
                    ns4=1+NKT(1,nb)+NKT(2,nb)+NKT(3,nb)
                    IF(ZE(ns1,nhx).GE.ZE(ns2,nhx))
     '                 ZE(ns2,nhx)=ZE(ns2,nhx)+2.d0*PI
                    IF(ZE(ns3,nhx).GE.ZE(ns4,nhx))
     '                 ZE(ns4,nhx)=ZE(ns4,nhx)+2.d0*PI
                    IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                    call mp_setlock()
                      WRITE(OP_STRING,'(/'' ZE(ns,'',I1,''): '','
     '                  //'5(1X,D11.4),/:(12X,5(1X,D11.4)))')
     '                  nhx,(ZE(ns1,nhx),ns1=1,NST(nb))
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                    call mp_unsetlock()
                    ENDIF
                  ENDIF
                ENDDO !nhj
              ELSE IF(ITYP10(nr).EQ.3) THEN !spherical coords
              ELSE IF(ITYP10(nr).EQ.4) THEN !prolate coords
              ELSE IF(ITYP10(nr).EQ.5) THEN !oblate coords
              ENDIF
            ENDIF !linear geometric fit

            IF(KTYP6.EQ.1) THEN !Gauss point fitting
              CALL YGER(NBH(1,1,ne),njj,nxFIT,ER,PG,SE(1,1,ne),
     '          YG(1,1,ne),ERROR,*9999)

            ELSE IF(KTYP6.EQ.4) THEN !Grid point fitting
              CALL YQSER(IBT,IDO,INP,NBH(1,1,ne),ne,NENQ,njj,
     '          NQNE,NQS,NQXI,
     '          nxFIT,ER,SE(1,1,ne),YQS,ERROR,*9999)


            ELSEIF(DATATYPE_NG.EQ.'NODE_GROUP') THEN
              DO nogroup=1,NTGRNO
                IF (LAGRNO(nogroup).EQ.GRNAME_NG) THEN
                  vorogrno=nogroup
                ENDIF
              ENDDO
              CALL ILIST_COPY(NLIGRNO(vorogrno),
     '          %VAL(LIGRNO_PTR(vorogrno)),NPLIST(1))
              nplist(0)=NLIGRNO(vorogrno)
              nplist2(0)=0 !nplist2 is np's of nodes in element ne
              DO nonode=1,nplist(0)
                np=nplist(nonode)
                IF (NEP(np).EQ.ne) THEN
                  nplist2(0)=nplist2(0)+1
                  nplist2(nplist2(0))=np
                ENDIF
              ENDDO
              CALL XPER(IBT,IDO,INP,NBH(1,1,ne),NBJ(1,ne),
     '          njj,NKJE(1,1,1,ne),NPF(1,1),nplist2,NPNE(1,1,ne),
     '          NRE(ne),NVJE(1,1,1,ne),nxFIT,ER,PG,RG,
     '          SE(1,1,ne),WDL,WG,WU(0,ne),XA(1,1,ne),XE,XG,XIP,
     '          XIDL,XP,ZDL,ZE,ZP,ERROR,*9999)
            ELSE  !Data fitting
              IF(KTYP11.EQ.1) THEN !element fitting
                CALL ZDER(IBT,IDO,INP,NBH(1,1,ne),NBJ(1,ne),NDDL,
     '            NDLT(ne),ne,njj,NKJE(1,1,1,ne),NPF(1,1),
     '            NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),nxFIT,ER,PG,
     '            RG,SE(1,1,ne),SF,WD,WDL,WG,WU(0,ne),XA(1,1,ne),XE,XG,
     '            XID,XIDL,XP,ZD,ZDL,ZE,ZP,IN_PLANE,ERROR,*9999)
              ELSE !face fitting
C PM 14Aug02 :
C               nef=NPF(8,ne)
C               CALL CALC_FACE_INFORMATION_IND(NBJ(1,elem),
C               '            NBJF(1,ne),nef,
C               '            NKJE(1,1,1,elem),NKEF,NKJF,NNF,
C               '            NPNE(1,1,elem),NPNF,nr,
C               '            NVJE(1,1,1,elem),NVJF,SE(1,1,elem),
C               '            SF,ERROR,*9999)
                CALL ZDER(IBT,IDO,INP,NBHF(1,1,nf),NBJ(1,ne),NDDL,
     '            NDLT(nf),nf,njj,NKJF,NPF,NPNE(1,1,ne),NRE(ne),NVJF,
     '            nxFIT,ER,PG,RG,SE(1,1,ne),SF,WD,WDL,WG,WU(0,nf),
     '            XA(1,1,ne),XE,XG,XID,XIDL,XP,ZD,ZDL,ZE,ZP,IN_PLANE,
     '            ERROR,*9999)
              ENDIF
              
            ENDIF !ktyp6
            
            IF(UPDATE_MATRIX) THEN
              IF(KTYP6.EQ.1) THEN !Gauss point fitting
                CALL YGES(NBH(1,1,ne),njj,nxFIT,ES,PG,SE(1,1,ne),WG,
     '            WU(0,ne),ERROR,*9999)


              ELSE IF(KTYP6.EQ.4) THEN !Grid point fitting
                  CALL YQSES(IBT,IDO,INP,NBH(1,1,ne),NBJ(1,ne),
     '              NDLT(ne),ne,NENQ,njj,NKJE(1,1,1,ne),NPF(1,1),
     '              NPNE(1,1,ne),NQNE,NQS,NQXI,
     '              NRE(ne),NVJE(1,1,1,ne),nxFIT,ES,PG,
     '              RG,SE(1,1,ne),SF,WDL,WG,WU(0,ne),XA(1,1,ne),XE,XG,
     '              XIDL,XP,IN_PLANE,ERROR,*9999)

              ELSEIF(DATATYPE_NG.EQ.'NODE_GROUP') THEN !node_group
                CALL XPES(IBT,IDO,INP,NBH(1,1,ne),NBJ(1,ne),
     '            njj,NKJE(1,1,1,ne),NPF(1,1),nplist2,NPNE(1,1,ne),
     '            NRE(ne),
     '            NVJE(1,1,1,ne),nxFIT,ES,PG,RG,SE(1,1,ne),WDL,WG,
     '            WU(0,ne),XA(1,1,ne),XE,XG,XIDL,XP,ERROR,*9999)
              ELSE                !Data fitting
                IF(KTYP11.EQ.1) THEN !element fitting
                  CALL ZDES(IBT,IDO,INP,NBH(1,1,ne),NBJ(1,ne),
     '              NDLT(ne),njj,NKJE(1,1,1,ne),NPF(1,1),
     '              NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),nxFIT,ES,PG,
     '              RG,SE(1,1,ne),SF,WDL,WG,WU(0,ne),XA(1,1,ne),XE,XG,
     '              XIDL,XP,IN_PLANE,ERROR,*9999)
                ELSEIF(KTYP11.EQ.2) THEN !face fitting
                  CALL ZDES(IBT,IDO,INP,NBHF(1,1,nf),NBJF(1,nf),
     '              NDLT(nf),njj,NKJF,NPF(1,1),NPNF,
     '              NRE(ne),NVJE(1,1,1,ne),nxFIT,ES,PG,RG,SF,
     '              SF,WDL,WG,WU(0,nf),XA(1,1,ne),XE,XG,XIDL,XP,
     '              IN_PLANE,ERROR,*9999)
                ENDIF
              ENDIF !ktyp6
            ENDIF

            IF(IWRIT4(nr,nxFIT).GE.5) THEN
              WRITE(OP_STRING,'(/'' Element '',I4,'' rhs vector '
     '          //'and stiffness matrix:'')') ne
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
c cpb 25/9/95 Adding generic stiffness matrix ouput
              CALL OPESTFMAT(NHST,IOOP,ES,ER,'ES','ER',.TRUE.,.TRUE.,
     '          ERROR,*9999)
            ENDIF

C*** Assemble element stiffness matrix into global system.

C PM 14Aug02 : Added for face fitting
            IF(KTYP11.EQ.1) THEN !element fitting
              DO nhs1=1,NHST(1)
                ny1=IABS(LGE(nhs1,1))
                GR(ny1)=GR(ny1)+ER(nhs1)
                IF(UPDATE_MATRIX) THEN
                  DO nhs2=1,NHST(2)
                    ny2=IABS(LGE(nhs2,2))
C CPB 25/9/95 Adding sparsity
C                    GK(ny1,ny2)=GK(ny1,ny2)+ES(nhs1,nhs2)
                    CALL SPARSE(ny1,ny2,NYT(1,1,nxFIT),nz,NZ_GK_M,
     '                NZT(1,nxFIT),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
                    GK(nz)=GK(nz)+ES(nhs1,nhs2)
                  ENDDO !nhs2
                ENDIF
              ENDDO !nhs1
            ELSEIF (KTYP11.EQ.2) THEN !face fitting
              DO nhs1=1,NHST_FACE(1)
                ny1=LGE_FACE(nhs1)
                GR(ny1)=GR(ny1)+ER(nhs1)
                IF(UPDATE_MATRIX) THEN
                  DO nhs2=1,NHST_FACE(2)
                    ny2=LGE_FACE(nhs2)
                    CALL SPARSE(ny1,ny2,NYT(1,1,nxFIT),nz,NZ_GK_M,
     '                NZT(1,nxFIT),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
                    GK(nz)=GK(nz)+ES(nhs1,nhs2)
                  ENDDO !nhs2
                ENDIF
              ENDDO !nhs1
            ENDIF

            CALL CPU_TIMER(CPU_USER,TIME_STOP)
            ELAPSED_TIME=TIME_STOP(1)-TIME_START4(1)
            AVETIME=AVETIME+ELAPSED_TIME
            NO_TIMES=NO_TIMES+1
            IF(IWRIT4(nr,nxFIT).GE.1) THEN
              WRITE(OP_STRING,'(/'' CPU time for element '',I5,'
     '          //''' assembly: '',D11.4,'' s'')') ne,ELAPSED_TIME
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF

          ENDDO !l (ne)

        ENDIF !global/element assembly


C------------------ end material and other fitting ---------------------

        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START3(1)
        IF(IWRIT4(nr,nxFIT).GE.1) THEN
          WRITE(OP_STRING,'(/'' CPU time for stiffness matrix assembly '
     '      //': '',D11.4,'' s'')') ELAPSED_TIME
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          IF(NO_TIMES.NE.0) THEN
            WRITE(OP_STRING,'('' Average CPU time per element for '
     '        //'assembly:'',D11.4,'' s'')') AVETIME/NO_TIMES
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF

        IF(UPDATE_MATRIX) THEN
          WRITE(OP_STRING,'(/'' Assembly of global matrix and RHS '','
     '      //'''complete'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ELSE
          WRITE(OP_STRING,'(/'' Assembly of RHS vector complete'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF !update_matrix

        IF(IWRIT4(nr,nxFIT).GE.4) THEN
          CALL CPU_TIMER(CPU_USER,TIME_START3)
          WRITE(OP_STRING,
     '      '(/'' Global load vector GR & stiffness matrix GK:'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NYNR(0,1,1,nr,nxFIT)='',I5,'
     '      //''', NYNR(0,2,1,nr,nxFIT)='',I5)')
     '      NYNR(0,1,1,nr,nxFIT),NYNR(0,2,1,nr,nxFIT)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)

C         ..generic stiffness matrix ouput
          CALL OPSTFMAT(NYNR(0,2,1,nr,nxFIT),ISC_GK,ISR_GK,IOOP,
     '      NYT(1,1,nxFIT),NYT(2,1,nxFIT),NZT(1,nxFIT),
     '      NYNR(0,1,1,nr,nxFIT),
     '      KTYP24,GK,GR,'GK ','GR ',.FALSE.,.TRUE.,.TRUE.,ERROR,*9999)

          CALL CPU_TIMER(CPU_USER,TIME_STOP)
          ELAPSED_TIME=TIME_STOP(1)-TIME_START3(1)
          WRITE(OP_STRING,'(/'' CPU time for stiffness matrices '
     '      //'output: '',D11.4,'' s'')') ELAPSED_TIME
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

        ENDIF !iwrit4

        CALL CPU_TIMER(CPU_USER,TIME_START3)

C       ..calc solution mapping arrays for the current fit variable

        CALL GLOBALF(IBT,IDO,INP,NBH,NBJ,NENP,njj,NKB,NKH(1,1,1,nr),
     '    NKHE,NKJE,NLL,NNB,NNF,NNL,NONY,NPF,NPL,NPNE,NPNODE,NPNY,
     '    nr,NVHE,NVHP,NVJE,NWP,nxFIT,NXI,NYNE,NYNO,NYNP,
     '    NYNR(0,0,1,0,nxFIT),NYNY,CONY,CYNO,CYNY,SE,SP,XA,XE,XP,
     '    FIX,ERROR,*9999)

        IF(NOT(2,1,nr,nxFIT).EQ.0) THEN
          ERROR=' >>The number of unknowns is zero'
          GOTO 9999
        ENDIF


C----------------------- generate reduced system -----------------------

C       ...allocate memory & initialise solution variables
        IF(FIRST_A) THEN

C         ..temporary work array allocation
          IF(SPARSEGKK(nxFIT).NE.0) THEN
            WORK_PTR=0
            CALL ALLOCATE_MEMORY(NOT(1,1,nr,nxFIT)*NOT(2,1,nr,nxFIT),1,
     '        CHARTYPE,WORK_PTR,MEM_INIT,ERROR,*9999)
          ENDIF
          CALL CALC_SPARSE_GKK(ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,ISR_GKK,
     '      ISR_GQ,LGE,NBH,NENP,NHE,NOT(1,1,nr,nxFIT),
     '      NOT(2,1,nr,nxFIT),NONY(0,1,1,nr),NP_INTERFACE,NPNE,
     '      NPNY(0,1,0,nxFIT),nr,NRE,NVHE,nxFIT,NYNE,NYNP,
     '      NYNR(0,0,1,nr,nxFIT),GK,GQ,%VAL(WORK_PTR),.FALSE.,.TRUE.,
     '      ERROR,*9999)
          IF(SPARSEGKK(nxFIT).NE.0)
     '      CALL FREE_MEMORY(WORK_PTR,ERROR,*9999)
          DO nzz=1,NZZT(1,nr,nxFIT)
            GKK(nzz)=0.0d0
          ENDDO !nzz

        ENDIF !update_matrix

        DO no1=1,NOT(1,1,nr,nxFIT)
          GRR(no1)=0.0d0
        ENDDO !no1

        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START3(1)
        IF(IWRIT4(nr,nxFIT).GE.1) THEN
          WRITE(OP_STRING,'(/'' CPU time for solution matrix setup and '
     '      //'initialisation: '',D11.4,'' s'')') ELAPSED_TIME
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

C       ..generate the reduced system of equations
        CALL CPU_TIMER(CPU_USER,TIME_START3)

        IF(KTYP8.EQ.7) THEN !material fitting

          DO no=1,NOTT !loop over material params
            GRR(no)=GR(no)
          ENDDO !no
          DO nz=1,NZZT(1,nr,nxFIT)
            GKK(nz)=GK(nz)
          ENDDO !nz

        ELSE !all other fitting

C         ..put fitting into standard form

          DO no_nynr1=1,NYNR(0,1,1,nr,nxFIT) !loop global rows of GK
            ny1=NYNR(no_nynr1,1,1,nr,nxFIT) !is row #
            IF(KTYP24.EQ.1) THEN !Compressed Row
              IF(UPDATE_MATRIX) THEN
                DO nz=ISR_GK(ny1),ISR_GK(ny1+1)-1
                  ny3=ISC_GK(nz) !local column #
C!!! CS 15/4/2003 dropping index on NPNY correctly using class
C!!! This seems like a severe bug that has gone unnoticed
C!!! Wonder what the consequences are....
                  ny2=GETNYR(1,NPNY(0,1,0,nxFIT),
     '              nr,0,2,ny3,NYNE,NYNP) !global column #
                  DO noy1=1,NONY(0,ny1,1,nr)
                    no1=NONY(noy1,ny1,1,nr) !solution row # attached to ny1
                    co1=CONY(noy1,ny1,1,nr) !row coupling coefficient
                    DO noy2=1,NONY(0,ny2,2,nr)
                      no2=NONY(noy2,ny2,2,nr) !solution var # attached to ny2
                      co2=CONY(noy2,ny2,2,nr) !column coupling coefficient
                      CALL SPARSE(no1,no2,NOT(1,1,nr,nxFIT),nzz,
     '                  NZ_GKK_M,NZZT(1,nr,nxFIT),ISC_GKK,ISR_GKK,
     '                  SPARSEGKK(nxFIT),ERROR,*9999)
                      IF(nzz.NE.0) GKK(nzz)=GKK(nzz)+GK(nz)*co1*co2
                    ENDDO !noy2
                  ENDDO !noy1
                ENDDO !nz
              ENDIF !UPDATE_MATRIX
              ny1=NYNR(no_nynr1,1,1,nr,nxFIT) !is row number
              DO noy1=1,NONY(0,ny1,1,nr)
                no1=NONY(noy1,ny1,1,nr) !solution row # attached to ny1
                co1=CONY(noy1,ny1,1,nr) !row coupling coefficient
                GRR(no1)=GRR(no1)+GR(ny1)*co1
              ENDDO !noy1
            ELSE !Compressed Row
              DO noy1=1,NONY(0,ny1,1,nr) !loop over #no's attached to ny1
                no1=NONY(noy1,ny1,1,nr) !no# attached to row ny1
                co1=CONY(noy1,ny1,1,nr) !coupling coeff for row mapping
C                                    ie row_no1=a*row_ny1+b*row_ny2
                GRR(no1)=GRR(no1)+GR(ny1)*co1 !get reduced R.H.S.vector
C             SAB 21Mar01  Shifted UPDATE_MATRIX up to here as the CONY
C             doesn't seem to be used.
                IF(UPDATE_MATRIX) THEN
C               SAB 21Mar01  Changed to use the NYNR columns of the local matrix
C                 rather than the global variable
                  DO no_nynr2=1,NYNR(0,2,1,nr,nxFIT) !loop over #cols of GK
                    ny2=NYNR(no_nynr2,0,1,nr,nxFIT) !is global variable #
                    ny3=GETNYR(1,NPNY(0,1,0,nxFIT),nr,2,0,ny2,NYNE,NYNP)
C                                         !local GK var #
                    CALL SPARSE(ny1,ny3,NYT(1,1,nxFIT),nz,NZ_GK_M,
     '                NZT(1,nxFIT),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
                    IF(nz.NE.0) THEN
                      DO noy2=1,NONY(0,ny2,2,nr) !loop over #no's for ny2
                        no2=NONY(noy2,ny2,2,nr) !no# attached to ny2
                        co2=CONY(noy2,ny2,2,nr) !coup coeff col mapping
C                                       i.e. var_no1=a*var_ny1+b*var_ny2
                        CALL SPARSE(no1,no2,NOT(1,1,nr,nxFIT),nzz,
     '                    NZ_GKK_M,NZZT(1,nr,nxFIT),ISC_GKK,ISR_GKK,
     '                    SPARSEGKK(nxFIT),ERROR,*9999)
                        IF(nzz.NE.0) GKK(nzz)=GKK(nzz)+GK(nz)*co1*co2
                      ENDDO !noy2
C!!! There is only one constant defined for each variable therefore
C!!! it is implied that this constant is applied to the last ny mapped
C!!! to the no. ie. no=a*ny1+b*ny2+(c*ny3+d) that is a*ny1=b*ny2 and
C!!! a*ny1=c*ny3+d and not a*ny1=b*ny2+d etc.
C                    SAB 21Mar01  The CONY(0,*) entries always seem to be
C                    zero so this isn't doing anything and we can then shift
C                    the UPDATE_MATRIX flag up
C                    co3=CONY(0,ny2,2,nr) !add.ve constant applied to vars
C                    GRR(no1)=GRR(no1)-GK(nz)*co3 !add.ve const in RHS vec
                    ENDIF
                  ENDDO !no_nynr2
                ENDIF
              ENDDO !noy1
            ENDIF !Compressed Row
          ENDDO !no_nynr1
        ENDIF !ktyp8

C--------------------------- timing and output -------------------------

        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START3(1)
        IF(IWRIT4(nr,nxFIT).GE.1) THEN
          WRITE(OP_STRING,'(/'' CPU time for solution matrix '
     '      //'assembly: '',D11.4,'' s'')') ELAPSED_TIME
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

        IF(KTYP4.NE.0) THEN !Output global matrices
          CALL WRITE_SOL_MATRIX(ISC_GKK,ISR_GKK,nr,nxFIT,GKK,GRR,
     '      ERROR,*9999)
        ENDIF

        IF(IWRIT4(nr,nxFIT).GE.3) THEN
          CALL CPU_TIMER(CPU_USER,TIME_START3)
          WRITE(OP_STRING,
     '      '(/'' Global load vector GRR & stiffness matrix GKK:'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NOT(1,1,nr,nxFIT)='',I5,'
     '      //''', NOT(2,1,nr,nxFIT)='',I5)')
     '      NOT(1,1,nr,nxFIT),NOT(2,1,nr,nxFIT)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
c cpb 24/9/95 Adding generic stiffness matrix output
          CALL OPSTFMAT(NYNR(0,1,1,nr,nxFIT),ISC_GKK,ISR_GKK,IOOP,
     '      NOT(1,1,nr,nxFIT),NOT(2,1,nr,nxFIT),NZZT(1,nr,nxFIT),
     '      NYNR(0,2,1,nr,nxFIT),SPARSEGKK(nxFIT),GKK,GRR,
     '      'GKK','GRR',.TRUE.,.TRUE.,.TRUE.,ERROR,*9999)
        ENDIF

        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START3(1)
        IF(IWRIT4(nr,nxFIT).GE.1) THEN
          IF(IWRIT4(nr,nxFIT).GE.4.OR.KTYP4.NE.0) THEN
            WRITE(OP_STRING,'(/'' CPU time for solution matrix '
     '        //'output: '',D11.4,'' s'')') ELAPSED_TIME
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF

C-------------- solve reduced system of linear equations ---------------

        IF(DOP) THEN
          WRITE(*,'(/''        NZT(1,nxFIT)='',I5,'
     '      //'     /''    NZZT(1,nr,nxFIT)='',I5,'
     '      //'     /''   NOT(1,1,nr,nxFIT)='',I5,'
     '      //'     /''   NOT(2,1,nr,nxFIT)='',I5,'
     '      //'     /'' SOLVEROPTION(nxFIT)='',I5)')
     '      NZT(1,nxFIT),NZZT(1,nr,nxFIT),
     '      NOT(1,1,nr,nxFIT),NOT(2,1,nr,nxFIT),
     '      SOLVEROPTION(nxFIT)
          DO no1=1,NOTT
            write(*,'('' no1='',I4,'' GKK:'',8E11.3,/(14X,8E11.3))')
     '        no1,(GKK((no1-1)*NOTT+no2),no2=1,NOTT)
          ENDDO
        ENDIF !dop

        CALL SOLVE_SYSTEM(ISC_GKK,ISR_GKK, NOT(1,1,nr,nxFIT),
     '    NOT(1,1,nr,nxFIT),NOT(2,1,nr,nxFIT),NZZT(1,nr,nxFIT),
     '    IWRIT4(nr,nxFIT),PRECON_CODE(nxFIT),SOLVEROPTION(nxFIT),
     '    SPARSEGKK(nxFIT),GKK,GRR,XO,FIRST_A,UPDATE_MATRIX,.TRUE.,
     '    nxFIT,ERROR,*9999)

        CALL CPU_TIMER(CPU_USER,TIME_START3)

        IF(DOP.AND.KTYP8.EQ.7) THEN
          WRITE(OP_STRING,'(/'' Residuals (Ax-b):'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          DO no1=1,NOTT
            RESIDUAL=0.d0
            DO no2=1,NOTT
              CALL SPARSE(no1,no2,NOTT,nz,NZ_GK_M,
     '          NZT(1,nxFIT),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
              RESIDUAL=RESIDUAL+GK(nz)*XO(no2)
            ENDDO !no2
            RESIDUAL=RESIDUAL-GR(no1)
            WRITE(OP_STRING,'('' XO('',I5,'')='',E13.4,'
     '        //''' Resid='',E15.6)') no1,XO(no1),RESIDUAL
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO !no1
        ENDIF !dop & ktyp8=7

        IF(KTYP8.EQ.7) THEN !material parameter fitting
          no=0 !initialise column variable
          DO nhj=1,NUM_FIT(1) !loop over material params
            nm=NMNO(1,nhj)
            DO nonode=1,NPNODE(0,nr) !loop over nodes
              np=NPNODE(nonode,nr)
              no=no+1 !increment column variable
              CP(nm,np,nxFIT)=XO(no)
            ENDDO !nonode
          ENDDO !nhj

        ELSE !all other fitting
          DO no1=1,NOT(2,1,nr,nxFIT)
            DO nyo1=1,NYNO(0,no1,2,nr)
              ny1=NYNO(nyo1,no1,2,nr)
              co1=CYNO(nyo1,no1,2,nr)
              IF(NPNY(0,ny1,0,nxFIT).EQ.1) THEN
                nk=NPNY(1,ny1,0,nxFIT)
                nv=NPNY(2,ny1,0,nxFIT)
                nh=NPNY(3,ny1,0,nxFIT)
                np=NPNY(4,ny1,0,nxFIT)
              ELSE
                ERROR='>> Element dofs not implemented'
                GOTO 9999
              ENDIF
              IF(KTYP6.EQ.1) THEN !Gauss fitting
                ZP(nk,nv,nh,np,1)=ZP(nk,nv,nh,np,1)+XO(no1)*co1 !no increment
              ELSE !Data fitting
                ZP(nk,nv,nh,np,1)=ZP(nk,nv,nh,np,1)+XO(no1)*co1!+incr
              ENDIF
            ENDDO !nyo1
          ENDDO !no1
        ENDIF !ktyp8


        IF(KTYP8.EQ.7) THEN !material fitting

        ELSE !all other fitting
          SMED=0.0d0
          SAED=0.0d0
          SQED=0.0d0
          NTOT=0
          nplist2(0)=0
          DO nd=1,NDT
            EDD(nd)=0.0d0
          ENDDO

          DO l=1,LN(0)

            IF(KTYP11.EQ.1) THEN !element fitting
              ne=LN(l)
            ELSEIF (KTYP11.EQ.2) THEN !face fitting
              nf=LN(l)
              ne=NPF(6,nf)
            ENDIF

C PM 14Aug02 :
C            IF(KTYP11.EQ.1) THEN !element fitting
C              CALL ZPZE(NBH(1,1,ne),1,NH_LOC(0,nxFIT),NKHE(1,1,1,ne),
C     '          NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1),
C     '          nxFIT,CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),
C     '          ZE,ZP,ERROR,*9999)
C            ELSE !face fitting
C              elem=NPF(6,ne)
C              CALL ZPZE_FACE(elem,ne,NBH(1,1,ne),NBJ,NBJF,1,
C     '          NH_LOC(0,nxFIT),NKEF,NKHE(1,1,1,ne),NKJE,NNF,
C     '          NPF,NPNE(1,1,elem),NPNF,
C     '          nr,NVHE(1,1,1,ne),NVJE,NVJF,NW(ne,1),nxFIT,
C     '          CURVCORRECT(1,1,1,elem),SE(1,1,elem),SF,
C     '          ZA(1,1,1,elem),
C     '          ZE,ZP,ERROR,*9999)
C            ENDIF

            CALL ZPZE(NBH(1,1,ne),1,NH_LOC(0,nxFIT),NKHE(1,1,1,ne),
     '        NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1),
     '        nxFIT,CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),
     '        ZE,ZP,ERROR,*9999)

            DO nhj=1,NUM_FIT(njj)
              nhx=NLH_FIT(nhj,3,njj)
              nh=NH_LOC(nhx,nxFIT)
              nb=NBH(nh,1,ne)
              IF(KTYP6.EQ.1) THEN !Gauss fitting
                DO ng=1,NGT(nb)
                  VALUE=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '              XIG(1,ng,nb),ZE(1,nhx))
                  EDG=YG(NG_FIT(nhj,njj),ng,ne)-VALUE
                  SMED=SMED+EDG
                  SAED=SAED+DABS(EDG)
                  SQED=SQED+EDG**2
                  NTOT=NTOT+1
                ENDDO !ng

              ELSE IF(KTYP6.EQ.4) THEN !Grid fitting
C                DO ng=1,NGT(nb)
C                  VALUE=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
C     '              XIG(1,ng,nb),ZE(1,nhx))
C                  EDG=YG(NG_FIT(nhj,njj),ng,ne)-VALUE
C                  SMED=SMED+EDG
C                  SAED=SAED+DABS(EDG)
C                  SQED=SQED+EDG**2
C                  NTOT=NTOT+1
C                ENDDO !ng

              ELSEIF(DATATYPE_NG.EQ.'NODE_GROUP') THEN !node group fit
                DO nonode=1,nplist(0)
                  np=nplist(nonode)
                  IF (NEP(np).EQ.ne) THEN
                    nplist2(0)=nplist2(0)+1
                    nplist2(nplist2(0))=np
                  ENDIF
                ENDDO
                nj=NJ_FIT(nhj,njj)
                DO npe=1,nplist2(0)
                  np=nplist2(npe)
                  Z(nh)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '              XIP(1,np),ZE(1,nhx))
                  EDD(np)=Z(nh)-XP(1,1,nj,np)
                  SMED=SMED+EDD(np)
                  SAED=SAED+DABS(EDD(np))
                  SQED=SQED+EDD(np)**2
                  NTOT=NTOT+1
                ENDDO !npe
              ELSE !Data fitting
                nj=NJ_FIT(nhj,njj)
                IF(KTYP11.EQ.1) THEN !element fitting
                  elem2=ne
                ELSEIF (KTYP11.EQ.2) THEN !face fitting
                  elem2=nf
                ENDIF
                DO ndl=1,NDLT(elem2)
                  nd=NDDL(elem2,ndl)
                  Z(nh)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '              XID(1,nd),ZE(1,nhx))

                  IF(TAG_FIT) THEN
                    EDD(nd)=EDD(nd)+WD(nj,nd)*(Z(nh)-ZD(nj,nd))
                  ELSE
                    EDD(nd)=Z(nh)-ZD(nj,nd)
                  ENDIF
C***            Note: this should check fibre is actually being fitted
                  IF(KTYP8.EQ.2.AND.(NJ_LOC(NJL_FIBR,0,nr).EQ.1.OR.
     '              NJ_LOC(NJL_FIBR,0,nr).EQ.2)) THEN
C                 fibre fit
                    IF(EDD(nd).GT.PI/2.0d0) THEN
                      EDD(nd)=EDD(nd)-PI
                      WRITE(OP_STRING,'('' >>Warning!! Fibre error has '
     '                  //'pi subtracted'')')
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    ELSE IF(EDD(nd).LT.-PI/2.d0) THEN
                      EDD(nd)=EDD(nd)+PI
                      WRITE(OP_STRING,'('' >>Warning!! Fibre error has '
     '                  //'pi added'')')
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    ENDIF
                  ELSE IF(KTYP8.EQ.2.AND.NJ_LOC(NJL_FIBR,0,nr).EQ.3)
     '                THEN !sheet fit
c                 IF(EDD(nd).GT.PI/2.d0) THEN
c                   EDD(nd)=EDD(nd)-PI
c                   WRITE(OP_STRING,'('' >>Warning!!! Sheet error has '
c     '               //'pi subtracted'')')
c                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
c                 ELSE IF(EDD(nd).LT.-PI/2.d0) THEN
c                   EDD(nd)=EDD(nd)+PI
c                   WRITE(OP_STRING,'('' >>Warning!!! Fibre error has '
c     '               //'pi added'')')
c                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
c                 ENDIF
                  ENDIF
                  SMED=SMED+EDD(nd)
                  SAED=SAED+DABS(EDD(nd))
                  SQED=SQED+EDD(nd)**2
                  NTOT=NTOT+1
                ENDDO !ndl
              ENDIF
            ENDDO !nhj
          ENDDO !l (ne)
        ENDIF !matl/other fitting
        IF(TAG_FIT) THEN
          SMED=0.0d0
          SAED=0.0d0
          SQED=0.0d0
          NTOT=0
          DO nd=1,NDT
            SMED=SMED+EDD(nd)
            SAED=SAED+DABS(EDD(nd))
            SQED=SQED+EDD(nd)**2
            NTOT=NTOT+1
          ENDDO
        ENDIF

        IF(NTOT.GT.1) THEN
          WRITE(OP_STRING,
     '       '(/'' Average error           : '',D12.6,'' +/- '',D12.6,'
     '      //'/'' Average absolute error  : '',D12.6,'' +/- '',D12.6,'
     '      //'/'' Root mean squared error : '',D12.6)')
     '      SMED/DBLE(NTOT),
     '      DSQRT((SQED-SMED**2/DBLE(NTOT))/DBLE(NTOT-1)),
     '      SAED/DBLE(NTOT),
     '      DSQRT((SQED-SAED**2/DBLE(NTOT))/DBLE(NTOT-1)),
     '      DSQRT(SQED/DBLE(NTOT))
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

C***    Store the fitted solution back into XP
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
          DO nhj=1,NUM_FIT(njj)
            nj=NLH_FIT(nhj,1,njj)
            nhx=NLH_FIT(nhj,3,njj)
            nh=NH_LOC(nhx,nxFIT)
            DO nv=1,NVHP(nh,np,1)
              DO nk=1,NKH(nh,np,1,nr)
                XP(nk,nv,nj,np)=ZP(nk,nv,nh,np,1)
              ENDDO !nk
            ENDDO !nv
          ENDDO !nhj
        ENDDO !nonode (np)

        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START3(1)
        IF(IWRIT4(nr,nxFIT).GE.1) THEN
          IF(IWRIT4(nr,nxFIT).GE.4.OR.KTYP4.NE.0) THEN
            WRITE(OP_STRING,'(/'' CPU time for solution storage '
     '        //'and error estimate: '',D11.4,'' s'')') ELAPSED_TIME
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF

        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
        IF(IWRIT4(nr,nxFIT).GE.1) THEN
          IF(IWRIT4(nr,nxFIT).GE.4.OR.KTYP4.NE.0) THEN
            WRITE(OP_STRING,'(/'' Total CPU time for fit variable '','
     '        //'I1,'': '',D11.4,'' s'')') njj,ELAPSED_TIME
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF

      ENDDO !njj

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START1(1)
      IF(IWRIT4(nr,nxFIT).GE.1) THEN
        IF(IWRIT4(nr,nxFIT).GE.4.OR.KTYP4.NE.0) THEN
          WRITE(OP_STRING,'(/'' Total CPU time for fitting '
     '      //'problem: '',D11.4,'' s'')') ELAPSED_TIME
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('FITFLD')
      RETURN
 9999 CALL ERRORS('FITFLD',ERROR)
      CALL EXITS('FITFLD')
      RETURN 1
      END



