      SUBROUTINE EVTRSF_DYNAM(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,
     '  IBT,IDO,INP,IPIV,ISC_GK,ISC_GKK,ISC_GQ,ISIZE_TBH,
     '  ISR_GK,ISR_GKK,ISR_GQ,LD,LGE,NAN,NBH,NBJ,NDET,NDIPOLES,
     '  NEELEM,NELIST,NENP,NGAP,NHE,NHP,NKB,NKH,NKHE,NKJE,NLL,NNB,NNF,
     '  NNL,NONY,NORD,NPF,NP_INTERFACE,NPL,NPLIST3,NPLIST4,
     '  NPNE,NPNODE,NPNY,NRE,NVHE,NVHP,NVJE,NW,NWP,nx,NXI,
     '  NYNE,NYNO,NYNP,NYNR,NYNY,NYQNR,CE,CGE,CONY,CP,CURVCORRECT,
     '  CYNO,CYNY,DET,DIPOLE_CEN,DIPOLE_DIR,DL,GD,GK,
     '  GKK,GQ,GR,GRR,PG,SE,SP,T_BH,WG,WK1_INV,WK2_INV,
     '  WK3_INV,WK4_INV,XA,XE,XG,XID,XIG,XO,XP,
     '  YG,YP,ZA,ZE,ZP,AT_NODES,CHECK,FIX,OPFILE,SVD,ERROR,*)

C#### Subroutine: EVTRSF_DYNAM
C###  Description:
C###    EVTRSF_DYNAM evaluates the transfer matrix T_BH for epi to
C###    torso mappings. In iptrsf the decision between a single vs
C###    double layer transfer matrix has been made.  The single layer
C###    transfer matrix is a mapping from epicardial to body surface
C###    potentials.This type of transfer matrix can be calculated either
C###    algebraicly or by multiple forward solutions (direct).  The
C###    algebraic approach has only been implemented for the homogenous
C###    torso, but the direct approach is general (i.e. should work
C###    for homogenouse and full torso models).
C###    The double layer transfer matrix is a mapping from transmembrane
C###    potentials to body surface potentials. It can only be
C###    constructed here by mutliple forward solutions.  A Poisson
C###    equation with a special source term must be solved inside
C###    the heart.  The construction of the right hand side is performed
C###    explicitly here and the no-flux boundary condition of the
C###    transmembrane potential across the heart surface is imposed.
C###
C**** Expression to evaluate for the algebraic homogeneous evaluation
C**** of the transfer matrix T_BH is
C**** inv(GK_bb-GQ_bh*inv[GQ_hh]*GK_hb)*(GQ_bh*inv[GQ_hh]*GK_hh-GK_bh)
C!!!! Currently T_BH matrix is set up square.  This is not needed.

      IMPLICIT NONE
      INCLUDE 'anal00.cmn'
      INCLUDE 'b00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'coup00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mxch.inc'
      INCLUDE 'solv00.cmn'
      INCLUDE 'tol00.cmn'
      INCLUDE 'trsf00.cmn'

!     Parameter List
      INTEGER DIPOLE_CEN_NTIME(NDIPOLEM,NRM,NXM),
     '  DIPOLE_DIR_NTIME(NDIPOLEM,NRM,NXM),IBT(3,NIM,NBFM),
     '  IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),IPIV(NY_TRANSFER_M),ISC_GK(NISC_GKM),
     '  ISC_GKK(NISC_GKKM,NXM),ISC_GQ(NISC_GQM),
     '  ISIZE_TBH(2),ISR_GK(NISR_GKM),
     '  ISR_GKK(NISR_GKKM,NXM),ISR_GQ(NISR_GQM),LD(NDM),
     '  LGE(NHM*NSM,NRCM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NDET(NBFM,0:NNM),NDIPOLES(NRM,NXM),NEELEM(0:NE_R_M,0:NRM),
     '  NELIST(0:NEM),
     '  NENP(NPM,0:NEPM,0:NRM),NGAP(NIM,NBM),NHE(NEM,NXM),
     '  NHP(NPM,0:NRM,NXM),NKB(2,2,2,NNM,NBFM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),
     '  NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),
     '  NONY(0:NOYM,NYM,NRCM,0:NRM,NXM),NORD(5,NE_R_M),
     '  NPF(9,NFM),
     '  NP_INTERFACE(0:NPM,0:3),NPL(5,0:3,NLM),
     '  NPLIST3(0:NPM),NPLIST4(0:NPM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),NRE(NEM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),NWP(NPM,2),nx,
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),
     '  NYNY(0:NYYM,NYM,NRM,NXM),
     '  NYQNR(0:NYQM,0:NRCM,NCM,0:NRM,NXM)
      REAL*8 CE(NMM,NEM,NXM),CGE(NMM,NGM,NEM,NXM),
     '  CONY(0:NOYM,NYM,NRCM,0:NRM,NXM),CP(NMM,NPM,NXM),
     '  CURVCORRECT(2,2,NNM,NEM),CYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),
     '  CYNY(0:NYYM,NYM,NRM,NXM),
     '  DET(NBFM,0:NNM,NGM,6),DIPOLE_CEN(4,0:NDIPTIMM,NDIPOLEM,NRM,NXM),
     '  DIPOLE_DIR(4,0:NDIPTIMM,NDIPOLEM,NRM,NXM),
     '  DL(3,NLM),
     '  GD(NZ_GD_M),GK(NZ_GK_M),
     '  GKK(NZ_GKK_M,NXM),GQ(NZ_GQ_M),GR(NYROWM),
     '  GRR(NOM),PG(NSM,NUM,NGM,NBM),
     '  SE(NSM,NBFM,NEM),SP(NKM,NBFM,NPM),
     '  T_BH(NY_TRANSFER_M,NY_TRANSFER_M),
     '  WG(NGM,NBM),WK1_INV(NY_TRANSFER_M,NY_TRANSFER_M),
     '  WK2_INV(NY_TRANSFER_M,NY_TRANSFER_M),
     '  WK3_INV(NY_TRANSFER_M,NY_TRANSFER_M),WK4_INV(NY_TRANSFER_M),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XG(NJM,NUM),XID(NIM,NDM),XIG(NIM,NGM,NBM),
     '  XO(NOM,NXM),
     '  XP(NKM,NVM,NJM,NPM),
     '  YG(NIYGM,NGM,NEM),YP(NYM,NIYM,NXM),
     '  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
      LOGICAL AT_NODES,CHECK,FIX(NYM,NIYFIXM,NXM),OPFILE,SVD

!     Local Variables
      INTEGER GETNYR,icol,icol2,IFAIL,irow,isize,
     '  ISIZE_GQ_HH,ISIZE_GQ_BB,LWORK,
     '  nb,nc,nd,ne,nh1,nh2,nhx1,nhx2,nj,nk1,nk2,nolist1,nolist2,nonr,
     '  no_nynr,no_nynr1,no_nynr2,
     '  no1,noy1,np,np1,np2,nr,nr_gkk,nr0,nr1,nv1,nv2,ny,
     '  ny1,ny2,ny3,ny4,nz,nz1
      REAL*8 co1,CONDITION_NUMBER,CORRECTION_DIPOLE(3),FACTOR,P(1),
     '  Q(1),SUM,SUM2,SUM4,SUMR,XDR(3)
      LOGICAL ASSEMBLE,FIRST_SOLVE,FULL_CHECK,
     '  GQ_ASSEM,UPDATE_MATRIX,UPDATE_VECTOR

!     Functions
      REAL*8 PXI

      CALL ENTERS('EVTRSF_DYNAM',*9999)

      IF(KTYP94.EQ.2.AND.NRT.GT.1) THEN
        WRITE(OP_STRING,'('' Warning: Algebraic Transfer matrix only '
     '  //'implemented for the homogenous case (and NRT>1 here?)'')')
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
      ENDIF

      TBH_OK=.FALSE.

C***    Check that the number of rows in GK and GQ are the same for
C***    both the heart (nplist3) and torso (nplist4) nodes.  Can only
C***    calculate the transfer matrix if this is satisfied.
      ISIZE_GK_HH=0
      ISIZE_GK_BB=0
      ISIZE_GQ_HH=0
      ISIZE_GQ_BB=0

      nr0=TRSF_NR_FIRST !region number of inner/heart surface
      nr1=TRSF_NR_SECOND !region number of second/body surface

      DO nolist1=1,NPLIST3(0) !list of heart nodes
        np1=NPLIST3(nolist1)
        DO nhx1=1,NHP(np1,nr0,nx)
          nh1=NH_LOC(nhx1,nx)
          DO nv1=1,NVHP(nh1,np1,1,nr0)
            DO nk1=1,MAX(NKH(nh1,np1,1,nr0)-KTYP93(1,nr0),1)
              ISIZE_GK_HH=ISIZE_GK_HH+1
            ENDDO !nk1
          ENDDO !nv1
          DO nv1=1,NVHP(nh1,np1,2,nr0)
            DO nk1=1,MAX(NKH(nh1,np1,2,nr0)-KTYP93(2,nr0),1)
              ISIZE_GQ_HH=ISIZE_GQ_HH+1
            ENDDO !nk1
          ENDDO !nv1
        ENDDO !nh1
      ENDDO !np1
      DO nolist1=1,NPLIST4(0) !list of torso nodes
        np1=NPLIST4(nolist1)
        DO nhx1=1,NHP(np1,nr1,nx)
          nh1=NH_LOC(nhx1,nx)
          DO nv1=1,NVHP(nh1,np1,1,nr1)
            DO nk1=1,MAX(NKH(nh1,np1,1,nr1)-KTYP93(1,nr1),1)
              ISIZE_GK_BB=ISIZE_GK_BB+1
            ENDDO !nk1
          ENDDO !nv1
          DO nv1=1,NVHP(nh1,np1,2,nr1)
            DO nk1=1,MAX(NKH(nh1,np1,2,nr1)-KTYP93(2,nr1),1)
              ISIZE_GQ_BB=ISIZE_GQ_BB+1
            ENDDO !nk1
          ENDDO !nv1
        ENDDO !nh1
      ENDDO !np1

      CALL ASSERT(ISIZE_GK_HH.EQ.ISIZE_GQ_HH,
     '  '>>GK and GQ heart submatrices are different sizes',
     '  ERROR,*9999)
      CALL ASSERT(ISIZE_GK_BB.EQ.ISIZE_GQ_BB,
     '  '>>GK and GQ torso submatrices are different sizes',
     '  ERROR,*9999)
      CALL ASSERT(ISIZE_GK_HH.LE.NY_TRANSFER_M,
     '  '>>ny_transfer_m too small',ERROR,*9999)
      CALL ASSERT(ISIZE_GK_BB.LE.NY_TRANSFER_M,
     '  '>>ny_transfer_m too small',ERROR,*9999)

      ISIZE_TBH(1)=ISIZE_GK_BB
      ISIZE_TBH(2)=ISIZE_GK_HH


C***  Initialise T_BH
      DO irow=1,NY_TRANSFER_M
        DO icol=1,NY_TRANSFER_M
          T_BH(irow,icol)=0.0d0
        ENDDO !icol
      ENDDO !irow

      IF(KTYP94.EQ.1) THEN !direct approach
C       Call globalh to set up NOT etc based on new fix values.
        DO nonr=1,TRSF_NRLIST(0)
          nr=TRSF_NRLIST(nonr)
          COUP_NRLIST(nonr,nx)=nr
          CALL GLOBALH(IBT,IDO,INP,NAN,NBH,NBJ,NELIST,NENP,
     '      NHE(1,nx),NKB,NKHE,NKJE,NLL,NNB,NNF,NNL,NONY(0,1,1,0,nx),
     '      NP_INTERFACE,NPF,NPL,NPNE,NPNY(0,1,0,nx),nr,NRE,
     '      NVHE,NVHP,NVJE,NWP,nx,NXI,NYNE,NYNO(0,1,1,0,nx),
     '      NYNP,NYNR(0,0,1,0,nx),NYNY(0,1,1,nx),NYQNR(0,0,1,0,nx),
     '      CONY(0,1,1,0,nx),CYNO(0,1,1,0,nx),CYNY(0,1,1,nx),
     '      SE,SP,XA,XE,XP,YP,FIX(1,1,nx),ERROR,*9999)
          COUP_NRLIST(nonr,nx)=nr
        ENDDO
        COUP_NRLIST(0,nx)=TRSF_NRLIST(0)
C       Call globalc to set up coupling coefficients for the list of
C       regions associated with the transfer matrix.
        IF(COUP_NRLIST(0,nx).GT.1) THEN
          CALL GLOBALC(NAN,NBH,NHE(1,nx),NONY(0,1,1,0,nx),
     '      NP_INTERFACE,NPNODE,NPNY(0,1,0,nx),NRE,NW,
     '      nx,NXI,NYNE,NYNO(0,1,1,0,nx),NYNP,NYNR(0,0,1,0,nx),
     '      CONY(0,1,1,0,nx),CYNO(0,1,1,0,nx),FIX,ERROR,*9999)
        ENDIF

        CALL_SOLV=.FALSE. !Have reset coupling coeff
        FIRST_SOLVE=.TRUE.

C       Assemble GK,GQ matrices (if it has not already been done).
        DO nonr=1,TRSF_NRLIST(0)
          nr=TRSF_NRLIST(nonr)
          IF(.NOT.ASSEMBLE_GLOBAL(nr,nx)) ASSEMBLE=.TRUE.
        ENDDO !nonr
        IF(ASSEMBLE) THEN
C         GK,GQ have not been assembled
          GQ_ASSEM=.FALSE.
          IF(TRSF_NRLIST(0).GT.1) THEN
            GQ_ASSEM=COUPLED_BEM(nx) !coupled_bem set up in DESOLV
          ENDIF
          IF(IS_COUPLED(nx)) THEN
            CALL INIT_SPARSE_MATRIX(NYNR(0,2,1,0,nx),
     '        NYNR(0,2,1,0,nx),ISC_GK,ISR_GK,0,0,0,NYT(1,1,nx),
     '        NZ_GK_M,NZT(1,nx),NYNR(0,1,1,0,nx),NYNR(0,1,1,0,nx),
     '        KTYP24,GK,ERROR,*9999)
            IF(COUPLED_BEM(nx)) THEN
              CALL INIT_SPARSE_MATRIX(NYNR(0,2,2,0,nx),
     '          NYNR(0,2,2,0,nx),ISC_GK,ISR_GQ,0,0,0,NYT(1,2,nx),
     '          NZ_GQ_M,NZT(2,nx),NYNR(0,1,2,0,nx),NYNR(0,1,2,0,nx),
     '          KTYP24,GQ,ERROR,*9999)
            ENDIF
          ENDIF
          DO nonr=1,TRSF_NRLIST(0)
            nr=TRSF_NRLIST(nonr)
            IF(IS_COUPLED(nx)) THEN
              nr_gkk=0
            ELSE
              nr_gkk=nr
            ENDIF
            IF(ITYP6(nr,nx).EQ.1) THEN !problem is linear
              IF(ITYP4(nr,nx).EQ.1) THEN !problem is f.e.s only
                CALL ASSEMBLE1(IBT,IDO,INP,ISC_GK,ISC_GKK(1,nx),
     '            ISC_GQ,ISR_GK,ISR_GKK(1,nx),ISR_GQ,NBH,NBJ,
     '            NEELEM,NHE(1,nx),NHP(1,nr,nx),NKJE,
     '            NONY(0,1,1,nr_gkk,nx),NORD,NPF,NP_INTERFACE,NPNE,
     '            NPNY(0,1,0,nx),nr,nr_gkk,NRE,NVHE,NVJE,NW,nx,NYNE,
     '            NYNP,NYNR(0,0,1,nr,nx),CE(1,1,nx),CGE(1,1,1,nx),
     '            CONY(0,1,1,nr_gkk,nx),CP(1,1,nx),CURVCORRECT,
     '            GK,GKK(1,nx),GQ,GR,PG,SE,WG,XA,XP,YG,
     '            GQ_ASSEM,.TRUE.,.TRUE.,ERROR,*9999)
              ELSE IF(ITYP4(nr,nx).EQ.2.OR.ITYP4(nr,nx).EQ.3) THEN
C               BEM
                CALL ASSEMBLE2(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,
     '            IBT,IDO,INP,ISC_GK,ISC_GKK(1,nx),ISC_GQ,ISR_GK,
     '            ISR_GKK(1,nx),ISR_GQ,NBH,NBJ,NDET,NDIPOLES,NEELEM,
     '            NENP,NGAP,NHE(1,nx),NHP(1,nr,nx),
     '            NKH(1,1,1,nr),NKHE,NKJE,NLL,NONY(0,1,1,nr_gkk,nx),
     '            NPF,NP_INTERFACE,NPNE,NPNODE,NPNY(0,1,0,nx),nr,nr_gkk,
     '            NRE,NVHP(1,1,1,nr),NVJE,NW,nx,NYNE,NYNP,
     '            NYNR(0,0,1,nr,nx),CE(1,1,nx),CONY(0,1,1,nr_gkk,nx),
     '            CURVCORRECT,DET,DIPOLE_CEN,
     '            DIPOLE_DIR,DL,GD,GK,GKK(1,nx),GQ,PG,
     '            SE,0.0d0,WG,XA,XE,XG,XIG,
     '            XP,.TRUE.,.TRUE.,ERROR,*9999)
              ENDIF
            ENDIF
            ASSEMBLE_GLOBAL(nr,nx)=.TRUE.
          ENDDO !nr
        ENDIF !assemble
C Need to solve even if a solution has been performed elsewhere since
C origanal solve may have had different bcs etc involved.
        ASSEMBLE_SOLUTION(0,nx)=.FALSE.

C For all Double layer transfer matrices
        IF(KTYP95.EQ.2.OR.KTYP95.EQ.3) THEN

C  We want the FACTOR=sigma_i/sigma_e :
C   sigma_i=const.Ge and sigma_e=(1+const)Ge thus FACTOR=const/(1+const)
C   Note: the number const is written up as k in the
C   electical imaging theory paper - see LKC/AJP.
C
C   After simplifying we get FACTOR=(s/k) where s and k are the values
C   entered in the material file - s=G_i and k=G_i+G_e
C
C  Note - factor2=1/(G_i+G_e). This is taken account for in the routines
C   XEPGKGQ and XEPGKGQ_3DL, XEPGKGQHYP_3DL and XEPGKGQHYP.
C   Should really be done here for speedup. This will produce a
C    modified GQ matrix.

          ne=NEELEM(1,TRSF_NR_FIRST)
          IF(DABS(CE(1,ne,nx)-CE(2,ne,nx)).LT.ZERO_TOL) THEN
            WRITE(OP_STRING,
     '        '('' >> WARNING: Negative extracellular conductivity'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
          FACTOR=CE(2,ne,nx)/CE(1,ne,nx) ! this should be (s/k)
        ELSE
          FACTOR=1.0d0
        ENDIF

C*** Do the standard transfer matrix

        IF(.NOT.CALL_ANAL_TRSF) THEN
C         Loop over columns of T_BH (inner surface dofs)
          icol=0
          DO nolist1=1,NPLIST3(0) !list of heart nodes
            np1=NPLIST3(nolist1)
            DO nhx1=1,NHP(np1,nr0,nx)
              nh1=NH_LOC(nhx1,nx)
              DO nv1=1,NVHP(nh1,np1,1,nr0)
                DO nk1=1,MAX(NKH(nh1,np1,1,nr0)-KTYP93(1,nr0),1)
                  icol=icol+1
                  ny4=NYNP(nk1,nv1,nh1,np1,0,1,nr0)
C                 global ny # for inner surface/heart dof
                  YP(ny4,1,nx)=1.0d0
C                 calculate GRR here since it will be quicker than
C                 doing it in the SOLVE routines. GRR is just the
C                 column of GK corresponding to ny4 (since YP(ny4, )=1
C                 and there should be no coupling on this inner
C                 surface [at least in the single layer case]).
                  UPDATE_VECTOR=.FALSE.
                  nr=nr0
                  IF(TRSF_NRLIST(0).GT.1) nr=0
                  DO no_nynr1=1,NYNR(0,1,1,nr,nx) !Loop over the # of
C                                                 rows of GK
                    ny1=NYNR(no_nynr1,1,1,nr,nx) !is the row #
                    DO noy1=1,NONY(0,ny1,1,nr,nx) !loop over the # of
C                                 no's attached to the row ny1
                      no1=NONY(noy1,ny1,1,nr,nx) !row number attached
C                                                to this ny
                      co1=CONY(noy1,ny1,1,nr,nx) !coupling coeff for
C                                               the row mapping
                      ny3=GETNYR(1,NPNY(0,1,0,nx),nr,2,0,ny4,NYNE,NYNP)
C                   the local GK variable # corresponding to ny4
                      CALL SPARSE(ny1,ny3,NYT(1,1,nx),nz1,
     '                  NZ_GK_M,NZT(1,nx),ISC_GK,ISR_GK,KTYP24,
     '                  ERROR,*9999)
                      IF(nz1.NE.0) THEN
                        GRR(no1)=-GK(nz1)*co1*FACTOR
                      ELSE
                        GRR(no1)=0.0d0
                      ENDIF
                    ENDDO !noy1
                  ENDDO !no_nynr1
                  UPDATE_MATRIX=.NOT.ASSEMBLE_SOLUTION(0,nx)
                  IF(TRSF_NRLIST(0).GT.1) THEN !Coupled solve
                    CALL SOLVE9(ISC_GK,ISC_GKK(1,nx),ISC_GQ,ISR_GK,
     '                ISR_GKK(1,nx),ISR_GQ,LGE,NBH,NENP,NHE(1,nx),
     '                NDIPOLES,NONY(0,1,1,0,nx),NP_INTERFACE,NPNE,
     '                NPNY(0,1,0,nx),NRE,NVHE,nx,NYNE,NYNO(0,1,1,0,nx),
     '                NYNP,NYNR(0,0,1,0,nx),CONY(0,1,1,0,nx),
     '                CYNO(0,1,1,0,nx),GD,GK,GKK(1,nx),GQ,GR,GRR,
     '                XO(1,nx),YP(1,1,nx),FIRST_SOLVE,FIX(1,1,nx),
     '                UPDATE_MATRIX,.FALSE.,UPDATE_VECTOR,ERROR,*9999)
                  ELSE !Only one solve region
                    CALL SOLVE2(ISC_GK,ISC_GKK(1,nx),ISC_GQ,ISR_GK,
     '                ISR_GKK(1,nx),ISR_GQ,LGE,NBH,NENP,NHE(1,nx),
     '                NDIPOLES,NONY(0,1,1,nr0,nx),NP_INTERFACE,NPNE,
     '                NPNY(0,1,0,nx),nr0,NRE,NVHE,nx,NYNE,
     '                NYNO(0,1,1,nr0,nx),NYNP,NYNR(0,0,1,nr0,nx),
     '                CONY(0,1,1,nr0,nx),CYNO(0,1,1,nr0,nx),GD,GK,
     '                GKK(1,nx),GQ,GRR,XO(1,nx),YP(1,1,nx),
     '                FIRST_SOLVE,FIX(1,1,nx),UPDATE_MATRIX,
     '                .FALSE.,UPDATE_VECTOR,ERROR,*9999)
                  ENDIF !TRSF_NRLIST(0)
                  ASSEMBLE_SOLUTION(0,nx)=.TRUE.
                  FIRST_SOLVE=.FALSE.
C                 Transfer solution to appropriate column of T_BH
                  irow=0
                  IF (AT_NODES) THEN
                    DO nolist2=1,NPLIST4(0) !list of body nodes
                      np2=NPLIST4(nolist2)
                      DO nhx2=1,NHP(np2,nr1,nx)
                        nh2=NH_LOC(nhx2,nx)
                        DO nv2=1,NVHP(nh2,np2,1,nr1)
                          DO nk2=1,
     '                      MAX(NKH(nh2,np2,1,nr1)-KTYP93(1,nr1),1)
                            ny3=NYNP(nk2,nv2,nh2,np2,0,1,nr1)
C                           global ny# for this var
                            irow=irow+1
                            T_BH(irow,icol)=YP(ny3,1,nx)
                            TBH_OK=(TBH_OK.OR.
     '                        T_BH(irow,icol).GT.ZERO_TOL)
                          ENDDO !nk2
                        ENDDO !nv2
                      ENDDO !nh2
                    ENDDO !np2
                  ELSE    ! at electrodes
                    CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr1,nx),
     '                NKH(1,1,1,nr1),NPNODE,nr1,NVHP(1,1,1,nr1),nx,
     '                NYNE,NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)
                    DO nolist2=1,NPLIST4(0) !list of body nodes
                      nd=NPLIST4(nolist2)
                      ne=LD(nd)
                      nh2=NH_LOC(1,nx)
                      nb=NBH(nh2,1,ne)
                      CALL ZPZE(NBH(1,1,ne),1,NHE(ne,nx),NKHE(1,1,1,ne),
     '                  NPF(1,1),NPNE(1,1,ne),NRE(ne),NVHE(1,1,1,ne),
     '                  NW(ne,1),nx,CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '                  ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
                      irow=irow+1
                      T_BH(irow,icol)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                  INP(1,1,nb),nb,1,XID(1,nd),ZE(1,nh2))
                      TBH_OK=(TBH_OK.OR.T_BH(irow,icol).GT.ZERO_TOL)
                    ENDDO !nd
                  ENDIF
                  YP(ny4,1,nx)=0.0d0 !reset the heart node to 0

                ENDDO !nk1
              ENDDO !nv1
            ENDDO !nh1
          ENDDO !nolist1


C*** Solve for an analytic solution

        ELSE

C cpb 15/1/01 Adding correction dipole for analytic solutions
          IF(NJT.EQ.2) THEN
            IF(ANAL_CHOICE(nr0).EQ.1) THEN
              CORRECTION_DIPOLE(1)=-1.6d0
              CORRECTION_DIPOLE(2)=-1.6d0
            ELSE IF(ANAL_CHOICE(nr0).EQ.2) THEN
              CORRECTION_DIPOLE(1)=-8.0d0
              CORRECTION_DIPOLE(2)=8.0d0
            ELSE IF(ANAL_CHOICE(nr0).EQ.3) THEN
              CORRECTION_DIPOLE(1)=-ANAL_A_1*ANAL_D_1
              CORRECTION_DIPOLE(2)=-ANAL_B_1*ANAL_D_1
            ELSE
              ERROR='>>Correction dipole not implemented'
              GOTO 9999
            ENDIF
          ELSE IF(NJT.EQ.3) THEN
            IF(ANAL_CHOICE(nr0).EQ.1) THEN !Case 1b
              CORRECTION_DIPOLE(1)=1.8525d0
              CORRECTION_DIPOLE(2)=2.0d0*1.8525d0
              CORRECTION_DIPOLE(3)=0.0d0
            ELSEIF(ANAL_CHOICE(nr0).EQ.2) THEN !Case 3a
              CORRECTION_DIPOLE(1)=1.6087d0
              CORRECTION_DIPOLE(2)=2.0d0*1.6087d0
              CORRECTION_DIPOLE(3)=1.9259d0*2.7115d0
            ELSEIF(ANAL_CHOICE(nr0).EQ.3) THEN !Case 2b
              CORRECTION_DIPOLE(1)=0.0d0
              CORRECTION_DIPOLE(2)=0.0d0
              CORRECTION_DIPOLE(3)=4.7407d0*3.3281d0
            ELSEIF(ANAL_CHOICE(nr0).EQ.4) THEN !Case 1a
              CORRECTION_DIPOLE(1)=1.6087d0
              CORRECTION_DIPOLE(2)=2.0d0*1.6087d0
              CORRECTION_DIPOLE(3)=0.0d0
            ELSEIF(ANAL_CHOICE(nr0).EQ.5) THEN !Case 2a
              CORRECTION_DIPOLE(1)=0.0d0
              CORRECTION_DIPOLE(2)=0.0d0
              CORRECTION_DIPOLE(3)=1.9259d0*2.7115d0
            ELSEIF(ANAL_CHOICE(nr0).EQ.6) THEN !Case 3b
              CORRECTION_DIPOLE(1)=1.8525d0
              CORRECTION_DIPOLE(2)=2.0d0*1.8525d0
              CORRECTION_DIPOLE(3)=4.7407d0*3.3281d0
            ELSE IF(ANAL_CHOICE(nr0).EQ.7) THEN !General solution
              CORRECTION_DIPOLE(1)=ANAL_A_11*ANAL_D_11
              CORRECTION_DIPOLE(2)=ANAL_B_11*ANAL_D_11
              CORRECTION_DIPOLE(3)=ANAL_D_1
            ELSE
              ERROR='>>Correction dipole not implemented'
              GOTO 9999
            ENDIF
          ENDIF

          DO no1=1,MAX(NOT(1,1,0,nx),NOT(2,1,0,nx))
            GRR(no1)=0.0d0
          ENDDO !no1

          icol=0
          np1=NPLIST3(1)
          DO nhx1=1,NHP(np1,nr0,nx)
            nh1=NH_LOC(nhx1,nx)
            DO nv1=1,NVHP(nh1,np1,1,nr0)
              DO nk1=1,MAX(NKH(nh1,np1,1,nr0)-KTYP93(1,nr0),1)
                icol=icol+1
                ny4=NYNP(nk1,nv1,nh1,np1,0,1,nr0)
C               global ny # for inner surface/heart dof

                UPDATE_VECTOR=.FALSE.
                nr=nr0
                IF(TRSF_NRLIST(0).GT.1) nr=0
                DO no_nynr1=1,NYNR(0,1,1,nr,nx) !Loop over the # of
C                                                 rows of GK
                  ny1=NYNR(no_nynr1,1,1,nr,nx) !is the row #
                  DO noy1=1,NONY(0,ny1,1,nr,nx) !loop over the # of
C                                                 no's attached to
C                                                 the row ny1
                    no1=NONY(noy1,ny1,1,nr,nx) !row number attached
C                                                to this ny
                    co1=CONY(noy1,ny1,1,nr,nx) !coupling coeff for
C                                               the row mapping

                    IF(ny1.LE.NPLIST3(0)) THEN !only for first block
                      nc=1 !move to where more appropriate
                      DO no_nynr2=1,NYNR(0,0,nc,nr,nx) !loop over GK col
                        ny2=NYNR(no_nynr2,0,nc,nr,nx) !global GK var #
                        ny3=GETNYR(nc,NPNY,nr,2,0,ny2,NYNE,NYNP)
C                     is local GK variable #
                        CALL SPARSE(ny1,ny3,NYT(1,nc,nx),nz1,NZ_GK_M,
     '                    NZT(nc,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C Only want to do correction factor for the inner sphere
                        IF(ny2.LE.NPLIST3(0)) THEN ! assumes linear?
                          GRR(no1)=GRR(no1)-GK(nz1)*YP(ny2,1,nx)
     '                      *co1 !*co2??
     '                      *FACTOR
                        ENDIF !ny2
                      ENDDO !no_nynr2

C cpb 15/1/01 Adding correction dipole for analytic solutions
                      np=NPNY(4,ny1,1,nx)
                      SUM=0.0d0
                      SUMR=0.0d0
                      DO nj=1,NJT
                        XDR(nj)=-XP(1,1,nj,np) !r is from source to field point
                        SUM=SUM+XDR(nj)*CORRECTION_DIPOLE(nj) !r.p
                        SUMR=SUMR+XDR(nj)*XDR(nj) !r*r
                      ENDDO
                      SUMR=DSQRT(SUMR)
                      IF(NJT.EQ.2) THEN
                        SUM=-SUM/(SUMR*SUMR)
                      ELSE IF(NJT.EQ.3) THEN
                        SUM=SUM/(SUMR*SUMR*SUMR)
                      ENDIF
                      GRR(no1)=GRR(no1)-SUM

                    ENDIF !ny1
                  ENDDO !nony1
                ENDDO !no_nynr1

                UPDATE_MATRIX=.NOT.ASSEMBLE_SOLUTION(0,nx)
                IF(TRSF_NRLIST(0).GT.1) THEN !Coupled solve
                  CALL SOLVE9(ISC_GK,ISC_GKK(1,nx),ISC_GQ,ISR_GK,
     '              ISR_GKK(1,nx),ISR_GQ,LGE,NBH,NENP,NHE(1,nx),
     '              NDIPOLES,NONY(0,1,1,0,nx),NP_INTERFACE,NPNE,
     '              NPNY(0,1,0,nx),NRE,NVHE,nx,NYNE,NYNO(0,1,1,0,nx),
     '              NYNP,NYNR(0,0,1,0,nx),CONY(0,1,1,0,nx),
     '              CYNO(0,1,1,0,nx),GD,GK,GKK(1,nx),GQ,GR,GRR,
     '              XO(1,nx),YP(1,1,nx),FIRST_SOLVE,FIX(1,1,nx),
     '              UPDATE_MATRIX,.FALSE.,UPDATE_VECTOR,ERROR,*9999)
                ELSE !Only one solve region
                  CALL SOLVE2(ISC_GK,ISC_GKK(1,nx),ISC_GQ,ISR_GK,
     '              ISR_GKK(1,nx),ISR_GQ,LGE,NBH,NENP,NHE(1,nx),
     '              NDIPOLES,NONY(0,1,1,nr0,nx),NP_INTERFACE,NPNE,
     '              NPNY(0,1,0,nx),nr0,NRE,NVHE,nx,NYNE,
     '              NYNO(0,1,1,nr0,nx),NYNP,NYNR(0,0,1,nr0,nx),
     '              CONY(0,1,1,nr0,nx),CYNO(0,1,1,nr0,nx),GD,GK,
     '              GKK(1,nx),GQ,GRR,XO(1,nx),YP(1,1,nx),
     '              FIRST_SOLVE,FIX(1,1,nx),UPDATE_MATRIX,
     '              .FALSE.,UPDATE_VECTOR,ERROR,*9999)
                ENDIF !TRSF_NRLIST(0)
                ASSEMBLE_SOLUTION(0,nx)=.TRUE.
                FIRST_SOLVE=.FALSE.

              ENDDO !nk1
            ENDDO !nv1
          ENDDO !nh1
        ENDIF !analytic solve or calculate transfer matrix


CC Put back solution from YP( ,6,nx)
CC  This will only make sense if Salu has not been used, since Salu
CC  overwrites this .
        IF(.NOT.SALU_CONSISTENCY(nx).AND..NOT.CALL_ANAL_TRSF) THEN
          DO nr=1,NRT
            DO nc=1,NCT(nr,nx)
              DO no_nynr=1,NYNR(0,0,nc,nr,nx)
                ny=NYNR(no_nynr,0,nc,nr,nx)
                YP(ny,1,nx)=YP(ny,6,nx)
                FIX(ny,1,nx)=FIX(ny,5,nx)
              ENDDO !no_nynr
            ENDDO !nc
          ENDDO !nr
        ENDIF !not Salu


C*** Check the transfer matrix is valid
        IF(CHECK) THEN

          IF(KTYP95.GE.2) THEN !double layer
            !Row sums of T_BH corresponding to potential should sum
            !to 0 since a closed double layer
            !should generate no external potential.
            WRITE(OP_STRING,'(/'' Potential row sums of T_BH '')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            IF(OPFILE) THEN
              WRITE(OP_STRING,'(/'' Potential row sums of T_BH '')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF

            DO irow=1,ISIZE_GK_BB
              SUM=0.0d0
              icol=0
              DO nolist1=1,NPLIST3(0) !list of heart nodes
                np1=NPLIST3(nolist1)
C                CALL NODE_CHANGE(np1,.FALSE.,ERROR,*9999)
                DO nhx1=1,NHP(np1,nr0,nx)
                  nh1=NH_LOC(nhx1,nx)
                  DO nv1=1,NVHP(nh1,np1,1,nr0)
                    DO nk1=1,MAX(NKH(nh1,np1,1,nr0)-KTYP93(1,nr0),1)
                      icol=icol+1
                      IF(nk1.EQ.1) SUM=SUM+T_BH(irow,icol)
                    ENDDO !nk1
                  ENDDO !nv1
                ENDDO !nh1
              ENDDO !np1
              WRITE(OP_STRING,'('' Row sum '',I4,'' = '',D12.4)')
     '          irow,SUM
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              IF(OPFILE) THEN
                WRITE(OP_STRING,'('' Row sum '',I4,'' = '',D12.4)')
     '            irow,SUM
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO !irow
          ENDIF !KTYP95
        ENDIF !CHECK


      ELSE IF(KTYP94.EQ.2) THEN !algebraic approach
        nr=TRSF_NRLIST(1)
C       Homogeneous at the moment, so inner and outer surfaces
C       are the same.
C***    Initialise work arrays
        DO irow=1,NY_TRANSFER_M
          WK4_INV(irow)=0.0d0
          DO icol=1,NY_TRANSFER_M
            WK1_INV(irow,icol)=0.0d0
            WK2_INV(irow,icol)=0.0d0
            WK3_INV(irow,icol)=0.0d0
          ENDDO !icol
        ENDDO !irow

        IF(CHECK) THEN
C***      Find norms of components of original equations to check on
C***      error accumulation in transfer matrix.
C***      Checking GQhh*dphihdn
          SUM=0.0d0
          SUM2=0.0d0
          FULL_CHECK=.TRUE.
          IF(NY_TRANSFER_M.LE.NYT(1,1,nx).OR.
     '      NY_TRANSFER_M.LE.NYT(1,2,nx))THEN
            WRITE(OP_STRING,'('' Cannot perform full check '','
     '        //'''since NY_TRANSFER_M not big enough'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            FULL_CHECK=.FALSE.
          ENDIF
          DO nolist1=1,NPLIST3(0) !list of heart nodes
            np1=NPLIST3(nolist1)
            DO nhx1=1,NHP(np1,nr0,nx)
              nh1=NH_LOC(nhx1,nx)
              DO nv1=1,NVHP(nh1,np1,2,nr0)
                DO nk1=1,MAX(NKH(nh1,np1,2,nr0)-KTYP93(2,nr0),1)
                  ny1=NYNP(nk1,nv1,nh1,np1,1,2,nr0) !heart row # of GQ
                  ny4=NYNP(nk1,nv1,nh1,np1,0,2,nr0) !global ny # for this var
                  IF(FULL_CHECK) WK4_INV(ny1)=0.0d0
                  SUM2=SUM2+YP(ny4,1,nx)*YP(ny4,1,nx)
                  SUM4=0.0d0
                  DO nolist2=1,NPLIST3(0) !list of heart nodes
                    np2=NPLIST3(nolist2)
                    DO nhx2=1,NHP(np2,nr0,nx)
                      nh2=NH_LOC(nhx2,nx)
                      DO nv2=1,NVHP(nh2,np2,2,nr0)
                        DO nk2=1,
     '                    MAX(NKH(nh2,np2,2,nr0)-KTYP93(2,nr0),1)
                          ny2=NYNP(nk2,nv2,nh2,np2,2,2,nr0)
C                         heart col number of GQ
                          ny3=NYNP(nk2,nv2,nh2,np2,0,2,nr0)
C                         global ny # for this variable
                          CALL SPARSE(ny1,ny2,NYT(1,2,nx),nz,
     '                      NZ_GQ_M,NZT(2,nx),ISC_GQ,ISR_GQ,KTYP24,
     '                      ERROR,*9999)
                          SUM4=SUM4+GQ(nz)*YP(ny3,1,nx) !GQhh*phih
                          IF(FULL_CHECK) WK4_INV(ny1)=
     '                      WK4_INV(ny1)-GQ(nz)*YP(ny3,1,nx)
                        ENDDO !nk2
                      ENDDO !nv2
                    ENDDO !nh2
                  ENDDO !nolist (np2)
                  SUM=SUM+SUM4*SUM4 !Finding the norm
                ENDDO !nk1
              ENDDO !nv1
            ENDDO !nh1
          ENDDO !nolist (np1)
          SUM=DSQRT(SUM)
          SUM2=DSQRT(SUM2)
          WRITE(OP_STRING,'('' Norm of dphihdn      ='',D12.4)')SUM2
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Norm of GQhh*dphihdn ='',D12.4)')SUM
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)

C***      Checking GQhb*dphibdn
          irow=0
          SUM=0.0d0
          SUM2=0.0d0
          DO nolist1=1,NPLIST3(0) !list of heart nodes
            np1=NPLIST3(nolist1)
            DO nhx1=1,NHP(np1,nr0,nx)
              nh1=NH_LOC(nhx1,nx)
              DO nv1=1,NVHP(nh1,np1,2,nr0)
                DO nk1=1,MAX(NKH(nh1,np1,2,nr0)-KTYP93(2,nr0),1)
                  ny1=NYNP(nk1,nv1,nh1,np1,1,2,nr0) !heart row # of GQ
                  irow=irow+1
                  SUM4=0.0d0
                  DO nolist2=1,NPLIST4(0) !list of body nodes
                    np2=NPLIST4(nolist2)
                    DO nhx2=1,NHP(np2,nr1,nx)
                      nh2=NH_LOC(nhx2,nx)
                      DO nv2=1,NVHP(nh2,np2,2,nr1)
                        DO nk2=1,
     '                    MAX(NKH(nh2,np2,2,nr1)-KTYP93(2,nr1),1)
                          ny2=NYNP(nk2,nv2,nh2,np2,2,2,nr1)
C                         body col number of GQ
                          ny3=NYNP(nk2,nv2,nh2,np2,0,2,nr1)
C                         global ny #
                          CALL SPARSE(ny1,ny2,NYT(1,2,nx),nz,
     '                      NZ_GQ_M,NZT(2,nx),ISC_GQ,ISR_GQ,KTYP24,
     '                      ERROR,*9999)
                          SUM4=SUM4+GQ(nz)*YP(ny3,1,nx)
                          IF(FULL_CHECK) WK4_INV(ny1)=
     '                      WK4_INV(ny1)-GQ(nz)*YP(ny3,1,nx)
                          IF(irow.EQ.1) THEN
                            SUM2=SUM2+YP(ny3,1,nx)*YP(ny3,1,nx)
                          ENDIF
                        ENDDO !nk2
                      ENDDO !nv2
                    ENDDO !nh2
                  ENDDO !nolist (np2)
                  SUM=SUM+SUM4*SUM4 !Finding the norm
                ENDDO !nk1
              ENDDO !nv1
            ENDDO !nh1
          ENDDO !nolist (np1)
          SUM=DSQRT(SUM)
          SUM2=DSQRT(SUM2)
          WRITE(OP_STRING,'('' Norm of dphibdn      ='',D12.4)')SUM2
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Norm of GQhb*dphibdn ='',D12.4)')SUM
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)

C***      Checking GKhh*phih
          SUM=0.0d0
          SUM2=0.0d0
          DO nolist1=1,NPLIST3(0) !list of heart nodes
          np1=NPLIST3(nolist1)
            DO nhx1=1,NHP(np1,nr0,nx)
              nh1=NH_LOC(nhx1,nx)
              DO nv1=1,NVHP(nh1,np1,1,nr0)
                DO nk1=1,MAX(NKH(nh1,np1,1,nr0)-KTYP93(1,nr0),1)
                  ny1=NYNP(nk1,nv1,nh1,np1,1,1,nr0) !heart row # of GK
                  ny4=NYNP(nk1,nv1,nh1,np1,0,1,nr0) !global ny# for this var.
                  SUM2=SUM2+YP(ny4,1,nx)*YP(ny4,1,nx)
                  SUM4=0.0d0
                  DO nolist2=1,NPLIST3(0) !list of heart nodes
                    np2=NPLIST3(nolist2)
                    DO nhx2=1,NHP(np2,nr0,nx)
                      nh2=NH_LOC(nhx2,nx)
                      DO nv2=1,NVHP(nh2,np2,1,nr0)
                        DO nk2=1,
     '                    MAX(NKH(nh2,np2,1,nr0)-KTYP93(1,nr0),1)
                          ny2=NYNP(nk2,nv2,nh2,np2,2,1,nr0)
                          !heart col number of GK
                          ny3=NYNP(nk2,nv2,nh2,np2,0,1,nr0)
C                         global ny # for this variable
                          CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,
     '                      NZ_GK_M,NZT(1,nx),ISC_GK,ISR_GK,KTYP24,
     '                      ERROR,*9999)
                          SUM4=SUM4+GK(nz)*YP(ny3,1,nx)
                          IF(FULL_CHECK) WK4_INV(ny1)=
     '                      WK4_INV(ny1)-GK(nz)*YP(ny3,1,nx)
                        ENDDO !nk2
                      ENDDO !nv2
                    ENDDO !nh2
                  ENDDO !nolist (np2)
                  SUM=SUM+SUM4*SUM4 !Finding the norm
                ENDDO !nk1
              ENDDO !nv1
            ENDDO !nh1
          ENDDO !nolist (np1)
          SUM=DSQRT(SUM)
          SUM2=DSQRT(SUM2)
          WRITE(OP_STRING,'('' Norm of phih         ='',D12.4)')SUM2
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Norm of GKhh*phih    ='',D12.4)')SUM
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)

C***      Checking GKhb*phib
          irow=0
          SUM=0.0d0
          SUM2=0.0d0
          DO nolist1=1,NPLIST3(0) !list of heart nodes
            np1=NPLIST3(nolist1)
            DO nhx1=1,NHP(np1,nr0,nx)
              nh1=NH_LOC(nhx1,nx)
              DO nv1=1,NVHP(nh1,np1,1,nr0)
                DO nk1=1,MAX(NKH(nh1,np1,1,nr0)-KTYP93(1,nr0),1)
                  ny1=NYNP(nk1,nv1,nh1,np1,1,1,nr0) !heart row # of GK
                  irow=irow+1
                  SUM4=0.0d0
                  DO nolist2=1,NPLIST4(0) !list of body nodes
                    np2=NPLIST4(nolist2)
                    DO nhx2=1,NHP(np2,nr1,nx)
                      nh2=NH_LOC(nhx2,nx)
                      DO nv2=1,NVHP(nh2,np2,1,nr1)
                        DO nk2=1,
     '                    MAX(NKH(nh2,np2,1,nr1)-KTYP93(1,nr1),1)
                          ny2=NYNP(nk2,nv2,nh2,np2,2,1,nr1)
C                         body col number of GK
                          ny3=NYNP(nk2,nv2,nh2,np2,0,1,nr1)
C                         global ny # for this variable
                          CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,
     '                      NZ_GK_M,NZT(1,nx),ISC_GK,ISR_GK,KTYP24,
     '                      ERROR,*9999)
                          SUM4=SUM4+GK(nz)*YP(ny3,1,nx)
                          IF(FULL_CHECK) WK4_INV(ny1)=
     '                      WK4_INV(ny1)-GK(nz)*YP(ny3,1,nx)
                          IF(irow.EQ.1) THEN
                            SUM2=SUM2+YP(ny3,1,nx)*YP(ny3,1,nx)
                          ENDIF
                        ENDDO !nk2
                      ENDDO !nv2
                    ENDDO !nh2
                  ENDDO !nolist (np2)
                  SUM=SUM+SUM4*SUM4 !Finding the norm
                ENDDO !nk1
              ENDDO !nv1
            ENDDO !nh1
          ENDDO !nolist (np1)
          SUM=DSQRT(SUM)
          SUM2=DSQRT(SUM2)
          WRITE(OP_STRING,'('' Norm of phib         ='',D12.4)')SUM2
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Norm of GKhb*phib    ='',D12.4)')SUM
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          IF(FULL_CHECK) THEN
            SUM=0.0d0
            DO nolist1=1,NPLIST3(0) !list of heart nodes
              np1=NPLIST3(nolist1)
              DO nhx1=1,NHP(np1,nr0,nx)
                nh1=NH_LOC(nhx1,nx)
                DO nv1=1,NVHP(nh1,np1,1,nr0)
                  DO nk1=1,MAX(NKH(nh1,np1,1,nr0)-KTYP93(1,nr0),1)
                    ny1=NYNP(nk1,nv1,nh1,np1,1,1,nr0)
C                   !heart row # of GK
                    SUM=WK4_INV(ny1)*WK4_INV(ny1)
                  ENDDO !nk1
                ENDDO !nv1
              ENDDO !nh1
            ENDDO !nolist (np1)
            SUM=DSQRT(SUM)
            WRITE(OP_STRING,'('' Norm of first eqtn   ='',D12.4)')SUM
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF !full_check

C***      Checking GQbh*dphihdn
          irow=0
          SUM=0.0d0
          DO nolist1=1,NPLIST4(0) !list of body nodes
            np1=NPLIST4(nolist1)
            DO nhx1=1,NHP(np1,nr1,nx)
              nh1=NH_LOC(nhx1,nx)
              DO nv1=1,NVHP(nh1,np1,2,nr1)
                DO nk1=1,MAX(NKH(nh1,np1,2,nr1)-KTYP93(2,nr1),1)
                  ny1=NYNP(nk1,nv1,nh1,np1,1,2,nr1) !body row # of GQ
                  irow=irow+1
                  IF(FULL_CHECK) WK4_INV(irow)=0.0d0
                  SUM4=0.0d0
                  DO nolist2=1,NPLIST3(0) !list of heart nodes
                    np2=NPLIST3(nolist2)
                    DO nhx2=1,NHP(np2,nr0,nx)
                      nh2=NH_LOC(nhx2,nx)
                      DO nv2=1,NVHP(nh2,np2,2,nr0)
                        DO nk2=1,
     '                    MAX(NKH(nh2,np2,2,nr0)-KTYP93(2,nr0),1)
                          ny2=NYNP(nk2,nv2,nh2,np2,2,2,nr0)
C                         heart col number of GQ
                          ny3=NYNP(nk2,nv2,nh2,np2,0,2,nr0)
C                         global ny # for this variable
                          CALL SPARSE(ny1,ny2,NYT(1,2,nx),nz,
     '                      NZ_GQ_M,NZT(2,nx),ISC_GQ,ISR_GQ,KTYP24,
     '                      ERROR,*9999)
                          SUM4=SUM4+GQ(nz)*YP(ny3,1,nx)
                          IF(FULL_CHECK) WK4_INV(irow)=
     '                      WK4_INV(irow)-GQ(nz)*YP(ny3,1,nx)
                        ENDDO !nk2
                      ENDDO !nv2
                    ENDDO !nh2
                  ENDDO !nolist (np2)
                  SUM=SUM+SUM4*SUM4 !Finding the norm
                ENDDO !nk1
              ENDDO !nv1
            ENDDO !nh1
          ENDDO !nolist (np1)
          SUM=DSQRT(SUM)
          WRITE(OP_STRING,'('' Norm of GQbh*dphihdn ='',D12.4)')SUM
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)

C***      Checking GQbb*dphibdn
          irow=0
          SUM=0.0d0
          DO nolist1=1,NPLIST4(0) !list of body nodes
            np1=NPLIST4(nolist1)
            DO nhx1=1,NHP(np1,nr1,nx)
              nh1=NH_LOC(nhx1,nx)
              DO nv1=1,NVHP(nh1,np1,2,nr1)
                DO nk1=1,MAX(NKH(nh1,np1,2,nr1)-KTYP93(2,nr1),1)
                  ny1=NYNP(nk1,nv1,nh1,np1,1,2,nr1) !body row # of GQ
                  irow=irow+1
                  SUM4=0.0d0
                  DO nolist2=1,NPLIST4(0) !list of body nodes
                    np2=NPLIST4(nolist2)
                    DO nhx2=1,NHP(np2,nr1,nx)
                      nh2=NH_LOC(nhx2,nx)
                      DO nv2=1,NVHP(nh2,np2,2,nr1)
                        DO nk2=1,
     '                    MAX(NKH(nh2,np2,2,nr1)-KTYP93(2,nr1),1)
                          ny2=NYNP(nk2,nv2,nh2,np2,2,2,nr1)
                          !body col number of GQ
                          ny3=NYNP(nk2,nv2,nh2,np2,0,2,nr1)
C                         global ny # for this variable
                          CALL SPARSE(ny1,ny2,NYT(1,2,nx),nz,
     '                      NZ_GQ_M,NZT(2,nx),ISC_GQ,ISR_GQ,KTYP24,
     '                      ERROR,*9999)
                          SUM4=SUM4+GQ(nz)*YP(ny3,1,nx)
                          IF(FULL_CHECK) WK4_INV(irow)=
     '                      WK4_INV(irow)-GQ(nz)*YP(ny3,1,nx)
                        ENDDO !nk2
                      ENDDO !nv2
                    ENDDO !nh2
                  ENDDO !nolist (np2)
                  SUM=SUM+SUM4*SUM4 !Finding the norm
                ENDDO !nk1
              ENDDO !nv1
            ENDDO !nh1
          ENDDO !nolist (np1)
          SUM=DSQRT(SUM)
          WRITE(OP_STRING,'('' Norm of GQbb*dphibdn ='',D12.4)')SUM
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)

C***      Checking GKbh*phih
          irow=0
          SUM=0.0d0
          DO nolist1=1,NPLIST4(0) !list of body nodes
            np1=NPLIST4(nolist1)
            DO nhx1=1,NHP(np1,nr1,nx)
              nh1=NH_LOC(nhx1,nx)
              DO nv1=1,NVHP(nh1,np1,1,nr1)
                DO nk1=1,MAX(NKH(nh1,np1,1,nr1)-KTYP93(1,nr1),1)
                  ny1=NYNP(nk1,nv1,nh1,np1,1,1,nr1) !body row # of GK
                  irow=irow+1
                  SUM4=0.0d0
                  DO nolist2=1,NPLIST3(0) !list of heart nodes
                    np2=NPLIST3(nolist2)
                    DO nhx2=1,NHP(np2,nr0,nx)
                      nh2=NH_LOC(nhx2,nx)
                      DO nv2=1,NVHP(nh2,np2,1,nr0)
                        DO nk2=1,
     '                    MAX(NKH(nh2,np2,1,nr0)-KTYP93(1,nr0),1)
                          ny2=NYNP(nk2,nv2,nh2,np2,2,1,nr0)
C                         heart col number of GK
                          ny3=NYNP(nk2,nv2,nh2,np2,0,1,nr0)
C                         global ny # for this variable
                          CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,
     '                      NZ_GK_M,NZT(1,nx),ISC_GK,ISR_GK,KTYP24,
     '                      ERROR,*9999)
                          SUM4=SUM4+GK(nz)*YP(ny3,1,nx)
                          IF(FULL_CHECK) WK4_INV(irow)=
     '                      WK4_INV(irow)-GK(nz)*YP(ny3,1,nx)
                        ENDDO !nk2
                      ENDDO !nv2
                    ENDDO !nh2
                  ENDDO !nolist (np2)
                  SUM=SUM+SUM4*SUM4 !Finding the norm
                ENDDO !nk1
              ENDDO !nv1
            ENDDO !nh1
          ENDDO !nolist (np1)
          SUM=DSQRT(SUM)
          WRITE(OP_STRING,'('' Norm of GKbh*phih    ='',D12.4)')SUM
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)

C***      Checking GKbb*phib
          irow=0
          SUM=0.0d0
          DO nolist1=1,NPLIST4(0) !list of body nodes
            np1=NPLIST4(nolist1)
            DO nhx1=1,NHP(np1,nr1,nx)
              nh1=NH_LOC(nhx1,nx)
              DO nv1=1,NVHP(nh1,np1,1,nr1)
                DO nk1=1,MAX(NKH(nh1,np1,1,nr1)-KTYP93(1,nr1),1)
                  ny1=NYNP(nk1,nv1,nh1,np1,1,1,nr1) !body row # of GK
                  irow=irow+1
                  SUM4=0.0d0
                  DO nolist2=1,NPLIST4(0) !list of body nodes
                    np2=NPLIST4(nolist2)
                    DO nhx2=1,NHP(np2,nr1,nx)
                      nh2=NH_LOC(nhx2,nx)
                      DO nv2=1,NVHP(nh2,np2,1,nr1)
                        DO nk2=1,
     '                    MAX(NKH(nh2,np2,1,nr1)-KTYP93(1,nr1),1)
                          ny2=NYNP(nk2,nv2,nh2,np2,2,1,nr1)
C                         body col number of GK
                          ny3=NYNP(nk2,nv2,nh2,np2,0,1,nr1)
C                         global ny # for this variable
                          CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,
     '                      NZ_GK_M,NZT(1,nx),ISC_GK,ISR_GK,KTYP24,
     '                      ERROR,*9999)
                          SUM4=SUM4+GK(nz)*YP(ny3,1,nx)
                          IF(FULL_CHECK) WK4_INV(irow)=
     '                      WK4_INV(irow)-GK(nz)*YP(ny3,1,nx)
                        ENDDO !nk2
                      ENDDO !nv2
                    ENDDO !nh2
                  ENDDO !nolist (np2)
                  SUM=SUM+SUM4*SUM4 !Finding the norm
                ENDDO !nk1
              ENDDO !nv1
            ENDDO !nh1
          ENDDO !nolist (np1)
          SUM=DSQRT(SUM)
          WRITE(OP_STRING,'('' Norm of GKbb*phib    ='',D12.4)')SUM
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          IF(FULL_CHECK) THEN
            irow=0
            SUM=0.0d0
            DO nolist1=1,NPLIST4(0) !list of body nodes
              np1=NPLIST4(nolist1)
              DO nhx1=1,NHP(np1,nr1,nx)
                nh1=NH_LOC(nhx1,nx)
                DO nv1=1,NVHP(nh1,np1,1,nr1)
                  DO nk1=1,MAX(NKH(nh1,np1,1,nr1)-KTYP93(1,nr1),1)
                    irow=irow+1
                    SUM=WK4_INV(irow)*WK4_INV(irow)
                  ENDDO !nk1
                ENDDO !nv1
              ENDDO !nh1
            ENDDO !nolist (np1)
            SUM=DSQRT(SUM)
            WRITE(OP_STRING,'('' Norm of first eqtn   ='',D12.4)')SUM
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF !full_check
        ENDIF !check

C***    Evaluate inv[GQ_hh] - store in WK1_INV
C***    Find entries in GQ matrix that comprise GQ_hh
        irow=0
        DO nolist1=1,NPLIST3(0) !list of heart nodes
          np1=NPLIST3(nolist1)
          DO nhx1=1,NHP(np1,nr0,nx)
            nh1=NH_LOC(nhx1,nx)
            DO nv1=1,NVHP(nh1,np1,2,nr0)
              DO nk1=1,MAX(NKH(nh1,np1,2,nr0)-KTYP93(2,nr0),1)
                ny1=NYNP(nk1,nv1,nh1,np1,1,2,nr0) !heart row number of GQ
                irow=irow+1
                icol=0
                DO nolist2=1,NPLIST3(0) !list of heart nodes
                  np2=NPLIST3(nolist2)
                  DO nhx2=1,NHP(np2,nr0,nx)
                    nh2=NH_LOC(nhx2,nx)
                    DO nv2=1,NVHP(nh2,np2,2,nr0)
                      DO nk2=1,
     '                  MAX(NKH(nh2,np2,2,nr0)-KTYP93(2,nr0),1)
                        ny2=NYNP(nk2,nv2,nh2,np2,2,2,nr0)
                        !heart col number of GQ
                        icol=icol+1
                        CALL SPARSE(ny1,ny2,NYT(1,2,nx),nz,
     '                    NZ_GQ_M,NZT(2,nx),ISC_GQ,ISR_GQ,KTYP24,
     '                    ERROR,*9999)
                        WK2_INV(irow,icol)=GQ(nz)
                      ENDDO !nk2
                    ENDDO !nv2
                  ENDDO !nh2
                ENDDO !nolist (np2)
              ENDDO !nk1
            ENDDO !nv1
          ENDDO !nh1
        ENDDO !nolist (np1)
        IFAIL=1

        CALL DGETRF(ISIZE_GQ_HH,ISIZE_GQ_HH,WK2_INV,NY_TRANSFER_M,
     '    IPIV,IFAIL)
        IF(IFAIL.NE.0)THEN
          WRITE(OP_STRING,*)' IFAIL=',IFAIL,' in GQ_hh inverse calc'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ERROR='>>IFAIL<>0 in DGETRF'
          GOTO 9999
        ENDIF
        CALL DGETRI(ISIZE_GQ_HH,WK2_INV,NY_TRANSFER_M,IPIV,
     '    WK4_INV,NY_TRANSFER_M,IFAIL)
        IF(IFAIL.NE.0)THEN
          WRITE(OP_STRING,*)' IFAIL=',IFAIL,' in GQ_hh inverse calc'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ERROR='>>IFAIL<>0 in DGETRI'
          GOTO 9999
        ENDIF
        DO irow=1,ISIZE_GQ_HH
          DO icol=1,ISIZE_GQ_HH
            WK1_INV(irow,icol)=WK2_INV(irow,icol)
          ENDDO
        ENDDO

        IF(OPFILE)THEN
          WRITE(OP_STRING,'(/'' GQ_hh inverse'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO irow=1,ISIZE_GQ_HH
            WRITE(OP_STRING,'('' Row '',I5,'' :'',5D12.4,'
     '        //'/:(12X,5D12.4))') irow,
     '        (WK1_INV(irow,icol),icol=1,ISIZE_GQ_HH)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDIF

C***    Form GQ_bh*inv[GQ_hh] - store in WK2_INV
        isize=0
        DO nolist2=1,NPLIST4(0) !list of body nodes
          np2=NPLIST4(nolist2)
          DO nhx2=1,NHP(np2,nr1,nx)
            nh2=NH_LOC(nhx2,nx)
            DO nv2=1,NVHP(nh2,np2,2,nr1)
              DO nk2=1,MAX(NKH(nh2,np2,2,nr1)-KTYP93(2,nr1),1)
                ny1=NYNP(nk2,nv2,nh2,np2,1,2,nr1) !body row number of GQ
                isize=isize+1
                DO icol=1,ISIZE_GQ_BB
                  SUM=0.0d0
                  irow=0
                  DO nolist1=1,NPLIST3(0) !list of heart nodes
                    np1=NPLIST3(nolist1)
                    DO nhx1=1,NHP(np1,nr0,nx)
                      nh1=NH_LOC(nhx1,nx)
                      DO nv1=1,NVHP(nh1,np1,1,nr0)
                        DO nk1=1,
     '                    MAX(NKH(nh1,np1,1,nr0)-KTYP93(1,nr0),1)
                          irow=irow+1
                          ny2=NYNP(nk1,nv1,nh1,np1,2,2,nr0)
                          !heart column number of GQ
                          CALL SPARSE(ny1,ny2,NYT(1,2,nx),nz,
     '                      NZ_GQ_M,NZT(2,nx),ISC_GQ,ISR_GQ,KTYP24,
     '                      ERROR,*9999)
                          SUM=SUM+GQ(nz)*WK1_INV(irow,icol)
                        ENDDO !nk1
                      ENDDO !nv1
                    ENDDO !nh1
                  ENDDO !np1
                  WK2_INV(isize,icol)=SUM
                ENDDO !icol
              ENDDO !nk2
            ENDDO !nv2
          ENDDO !nh2
        ENDDO !nolist2 (np2)

        IF(OPFILE)THEN
          WRITE(OP_STRING,'(/'' GQ_bh*inv[GQ_hh]'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO irow=1,ISIZE_GQ_BB
            WRITE(OP_STRING,'('' Row '',I5,'' :'',5D12.4,'
     '        //'/:(12X,5D12.4))') irow,
     '        (WK2_INV(irow,icol),icol=1,ISIZE_GQ_HH)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDIF

C***    Form {GQ_bh*inv[GQ_hh]}*GK_hb - store in WK3_INV
        DO irow=1,ISIZE_GQ_BB !rows of {GQ_BH*inv[GQ_hh]}
          icol=0
          DO nolist2=1,NPLIST4(0) !list of body surface nodes
            np2=NPLIST4(nolist2)
            DO nhx2=1,NHP(np2,nr1,nx)
              nh2=NH_LOC(nhx2,nx)
              DO nv2=1,NVHP(nh2,np2,1,nr1)
                DO nk2=1,MAX(NKH(nh2,np2,1,nr1)-KTYP93(1,nr1),1)
                ny2=NYNP(nk2,nv2,nh2,np2,2,1,nr1) !body column # of GK
                  SUM=0.0d0
                  icol=icol+1
                  icol2=0
                  DO nolist1=1,NPLIST3(0) !list of heart nodes
                    np1=NPLIST3(nolist1)
                    DO nhx1=1,NHP(np1,nr0,nx)
                      nh1=NH_LOC(nhx1,nx)
                      DO nv1=1,NVHP(nh1,np1,1,nr0)
                        DO nk1=1,
     '                    MAX(NKH(nh1,np1,1,nr0)-KTYP93(1,nr0),1)
                          icol2=icol2+1
                          ny1=NYNP(nk1,nv1,nh1,np1,1,1,nr0)
                          !heart row number of GK
                          CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,
     '                      NZ_GK_M,NZT(1,nx),ISC_GK,ISR_GK,KTYP24,
     '                      ERROR,*9999)
                          SUM=SUM+WK2_INV(irow,icol2)*GK(nz)
                        ENDDO !nk1
                      ENDDO !nv1
                    ENDDO !nh1
                  ENDDO !np1
                  WK3_INV(irow,icol)=SUM
                ENDDO !nk2
              ENDDO !nv2
            ENDDO !nh2
          ENDDO !nolist2 (np2)
        ENDDO !irow

        IF(OPFILE)THEN
          WRITE(OP_STRING,'(/'' GQ_bh*inv[GQ_hh]*GK_hb'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO irow=1,ISIZE_GQ_BB
            WRITE(OP_STRING,'('' Row '',I5,'' :'',5D12.4,'
     '        //'/:(12X,5D12.4))') irow,
     '        (WK3_INV(irow,icol),icol=1,ISIZE_GQ_BB)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDIF

C***    Form GK_bb-GQ_bh*[inv[GQ_hh]*GK_hb] - overwrite WK3_INV
        irow=0
        DO nolist2=1,NPLIST4(0) !list of body nodes
          np2=NPLIST4(nolist2)
          DO nhx2=1,NHP(np2,nr1,nx)
            nh2=NH_LOC(nhx2,nx)
            DO nv2=1,NVHP(nh2,np2,2,nr1)
              DO nk2=1,MAX(NKH(nh2,np2,2,nr1)-KTYP93(2,nr1),1)
                ny1=NYNP(nk2,nv2,nh2,np2,1,1,nr1) !body row number of GK
                irow=irow+1
                icol=0
                DO nolist1=1,NPLIST4(0) !list of body nodes
                  np1=NPLIST4(nolist1)
                  DO nhx1=1,NHP(np1,nr1,nx)
                    nh1=NH_LOC(nhx1,nx)
                    DO nv1=1,NVHP(nh1,np1,1,nr1)
                      DO nk1=1,
     '                  MAX(NKH(nh1,np1,1,nr1)-KTYP93(1,nr1),1)
                        ny2=NYNP(nk1,nv1,nh1,np1,2,1,nr1)
                        !body column number of GK
                        icol=icol+1
                        CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,
     '                    NZ_GK_M,NZT(1,nx),ISC_GK,ISR_GK,KTYP24,
     '                    ERROR,*9999)
                        WK3_INV(irow,icol)=GK(nz)-WK3_INV(irow,icol)
                      ENDDO !nk1
                    ENDDO !nv1
                  ENDDO !nh1
                ENDDO !nolist1 (np1)
              ENDDO !nk2
            ENDDO !nv2
          ENDDO !nh2
        ENDDO !nolist2 (np2)

        IF(CHECK) THEN
C***      Checking  {GK_bb-GQ_bh*inv[GQ_hh]*GK_hb}*phib
          irow=0
          SUM=0.0d0
          DO nolist2=1,NPLIST4(0) !list of body nodes
            np2=NPLIST4(nolist2)
            DO nhx2=1,NHP(np2,nr1,nx)
              nh2=NH_LOC(nhx2,nx)
              DO nv2=1,NVHP(nh2,np2,1,nr1)
                DO nk2=1,MAX(NKH(nh2,np2,1,nr1)-KTYP93(1,nr1),1)
                  irow=irow+1
                  icol=0
                  SUM4=0.0d0
                  DO nolist1=1,NPLIST4(0) !list of body nodes
                    np1=NPLIST4(nolist1)
                    DO nhx1=1,NHP(np1,nr1,nx)
                      nh1=NH_LOC(nhx1,nx)
                      DO nv1=1,NVHP(nh1,np1,1,nr1)
                        DO nk1=1,
     '                    MAX(NKH(nh1,np1,1,nr1)-KTYP93(1,nr1),1)
                          ny3=NYNP(nk1,nv1,nh1,np1,0,1,nr1)
C                         global ny # for this variable
                          icol=icol+1
                          SUM4=SUM4+WK3_INV(irow,icol)*YP(ny3,1,nx)
                        ENDDO !nk1
                      ENDDO !nv1
                    ENDDO !nh1
                  ENDDO !nolist1 (np1)
                  SUM=SUM+SUM4*SUM4 !Finding the norm
                ENDDO !nk2
              ENDDO !nv2
            ENDDO !nh2
          ENDDO !nolist2 (np2)
          SUM=DSQRT(SUM)
          WRITE(OP_STRING,'('' Norm of (..lhs..)*phib ='',D12.4)')SUM
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF !check

        IF(OPFILE)THEN
          WRITE(OP_STRING,'(/'' GK_bb-GQ_bh*inv[GQ_hh]*GK_hb'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO irow=1,ISIZE_GQ_BB
            WRITE(OP_STRING,'('' Row '',I5,'' :'',5D12.4,'
     '        //'/:(12X,5D12.4))') irow,
     '        (WK3_INV(irow,icol),icol=1,ISIZE_GQ_BB)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDIF

C***    Form inv{GK_bb-GQ_bh*inv[GQ_hh]*GK_hb} - overwrite WK1_INV
        IFAIL=1

        CALL DGETRF(ISIZE_GK_BB,ISIZE_GK_BB,WK3_INV,NY_TRANSFER_M,
     '    IPIV,IFAIL)
        IF(IFAIL.NE.0)THEN
          WRITE(OP_STRING,*)' IFAIL=',IFAIL,' in second inverse calc'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ERROR='>>IFAIL<>0 in DGETRF'
          GOTO 9999
        ENDIF
        CALL DGETRI(ISIZE_GK_BB,WK3_INV,NY_TRANSFER_M,IPIV,
     '    WK4_INV,NY_TRANSFER_M,IFAIL)
        IF(IFAIL.NE.0)THEN
          WRITE(OP_STRING,*)' IFAIL=',IFAIL,' in second inverse calc'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ERROR='>>IFAIL<>0 in DGETRI'
          GOTO 9999
        ENDIF
        DO irow=1,ISIZE_GK_BB
          DO icol=1,ISIZE_GK_BB
            WK1_INV(irow,icol)=WK3_INV(irow,icol)
          ENDDO
        ENDDO

        IF(OPFILE)THEN
          WRITE(OP_STRING,'(/''inv[GK_bb-GQ_bh*inv[GQ_hh]*GK_hb]'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO irow=1,ISIZE_GQ_BB
            WRITE(OP_STRING,'('' Row '',I5,'' :'',5D12.4,'
     '        //'/:(12X,5D12.4))') irow,
     '        (WK1_INV(irow,icol),icol=1,ISIZE_GQ_BB)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDIF

C***    Form {GQ_bh*inv[GQ_hh]}*GK_hh - store in WK3_INV
        DO irow=1,ISIZE_GQ_BB !rows of GQ_bb*inv[GK_hh]
          icol=0
          DO nolist2=1,NPLIST3(0) !list of heart surface nodes
            np2=NPLIST3(nolist2)
            DO nhx2=1,NHP(np2,nr0,nx)
              nh2=NH_LOC(nhx2,nx)
              DO nv2=1,NVHP(nh2,np2,1,nr0)
                DO nk2=1,MAX(NKH(nh2,np2,1,nr0)-KTYP93(1,nr0),1)
                  ny2=NYNP(nk2,nv2,nh2,np2,2,1,nr0) !heart column # of GK
                  icol=icol+1
                  SUM=0.0d0
                  icol2=0
                  DO nolist1=1,NPLIST3(0) !list of heart nodes
                    np1=NPLIST3(nolist1)
                    DO nhx1=1,NHP(np1,nr0,nx)
                      nh1=NH_LOC(nhx1,nx)
                      DO nv1=1,NVHP(nh1,np1,1,nr0)
                        DO nk1=1,
     '                    MAX(NKH(nh1,np1,1,nr0)-KTYP93(1,nr0),1)
                          icol2=icol2+1
                          ny1=NYNP(nk1,nv1,nh1,np1,1,1,nr0)
                          !heart row number of GK
                          CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,
     '                      NZ_GK_M,NZT(1,nx),ISC_GK,ISR_GK,KTYP24,
     '                      ERROR,*9999)
                          SUM=SUM+WK2_INV(irow,icol2)*GK(nz)
                        ENDDO !nk1
                      ENDDO !nv1
                    ENDDO !nh1
                  ENDDO !np1
                  WK3_INV(irow,icol)=SUM
                ENDDO !nk2
              ENDDO !nv2
            ENDDO !nh2
          ENDDO !nolist2 (np2)
        ENDDO !irow

        IF(OPFILE)THEN
          WRITE(OP_STRING,'(/'' GQ_bh*inv[GQ_hh]*GK_hh'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO irow=1,ISIZE_GQ_BB
            WRITE(OP_STRING,'('' Row '',I5,'' :'',5D12.4,'
     '        //'/:(12X,5D12.4))') irow,
     '        (WK3_INV(irow,icol),icol=1,ISIZE_GQ_HH)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDIF

C***    Form {GQ_bh*inv[GQ_hh]*GK_hh}-GK_bh - overwrite WK3_INV
        irow=0
        DO nolist2=1,NPLIST4(0) !list of body nodes
          np2=NPLIST4(nolist2)
          DO nhx2=1,NHP(np2,nr1,nx)
            nh2=NH_LOC(nhx2,nx)
            DO nv2=1,NVHP(nh2,np2,2,nr1)
              DO nk2=1,MAX(NKH(nh2,np2,2,nr1)-KTYP93(2,nr1),1)
                ny1=NYNP(nk2,nv2,nh2,np2,1,1,nr1) !body row number of GK
                irow=irow+1
                icol=0
                DO nolist1=1,NPLIST3(0) !list of heart nodes
                  np1=NPLIST3(nolist1)
                  DO nhx1=1,NHP(np1,nr0,nx)
                    nh1=NH_LOC(nhx1,nx)
                    DO nv1=1,NVHP(nh1,np1,1,nr0)
                      DO nk1=1,
     '                  MAX(NKH(nh1,np1,1,nr0)-KTYP93(1,nr0),1)
                        ny2=NYNP(nk1,nv1,nh1,np1,2,1,nr0)
                        !heart column number of GK
                        icol=icol+1
                        CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,
     '                    NZ_GK_M,NZT(1,nx),ISC_GK,ISR_GK,KTYP24,
     '                    ERROR,*9999)
                        WK3_INV(irow,icol)=WK3_INV(irow,icol)-GK(nz)
                      ENDDO !nk1
                    ENDDO !nv1
                  ENDDO !nh1
                ENDDO !nolist1 (np1)
              ENDDO !nk2
            ENDDO !nv2
          ENDDO !nh2
        ENDDO !nolist2 (np2)

        IF(CHECK) THEN
C***      Checking {GQ_bh*inv[GQ_hh]*GK_hh-GK_bh}*phih
          irow=0
          SUM=0.0d0
          DO nolist2=1,NPLIST4(0) !list of body nodes
            np2=NPLIST4(nolist2)
            DO nhx2=1,NHP(np2,nr1,nx)
              nh2=NH_LOC(nhx2,nx)
              DO nv2=1,NVHP(nh2,np2,1,nr1)
                DO nk2=1,MAX(NKH(nh2,np2,1,nr1)-KTYP93(1,nr1),1)
                  irow=irow+1
                  icol=0
                  SUM4=0.0d0
                  DO nolist1=1,NPLIST3(0) !list of heart nodes
                    np1=NPLIST3(nolist1)
                    DO nhx1=1,NHP(np1,nr0,nx)
                      nh1=NH_LOC(nhx1,nx)
                      DO nv1=1,NVHP(nh1,np1,1,nr0)
                        DO nk1=1,
     '                    MAX(NKH(nh1,np1,1,nr0)-KTYP93(1,nr0),1)
                          ny3=NYNP(nk1,nv1,nh1,np1,0,1,nr0)
C                         global ny # for this variable
                          icol=icol+1
                          SUM4=SUM4+WK3_INV(irow,icol)*YP(ny3,1,nx)
                        ENDDO !nk1
                      ENDDO !nv1
                    ENDDO !nh1
                  ENDDO !nolist1 (np1)
                  SUM=SUM+SUM4*SUM4 !Finding the norm
                ENDDO !nk2
              ENDDO !nv2
            ENDDO !nh2
          ENDDO !nolist2 (np2)
          SUM=DSQRT(SUM)
          WRITE(OP_STRING,'('' Norm of (..rhs...)*phih ='',D12.4)')
     '      SUM
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF !check

        IF(OPFILE)THEN
          WRITE(OP_STRING,'(/'' GQ_bh*inv[GQ_hh]*GK_hh-GK_bh'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO irow=1,ISIZE_GQ_BB
            WRITE(OP_STRING,'('' Row '',I5,'' :'',5D12.4,'
     '        //'/:(12X,5D12.4))') irow,
     '        (WK3_INV(irow,icol),icol=1,ISIZE_GQ_HH)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDIF

C***      Form final product - T_BH the transfer matrix
        DO irow=1,ISIZE_GK_BB
          DO icol=1,ISIZE_GQ_HH
            SUM=0.0d0
            DO isize=1,ISIZE_GK_BB
              SUM=SUM+WK1_INV(irow,isize)*WK3_INV(isize,icol)
            ENDDO
            T_BH(irow,icol)=SUM
            TBH_OK=(TBH_OK.OR.SUM.GT.ZERO_TOL)
          ENDDO
        ENDDO

      ENDIF !KTYP94

      IF(OPFILE)THEN
        WRITE(OP_STRING,'(/'' Transfer Matrix'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO irow=1,ISIZE_GK_BB
          WRITE(OP_STRING,'('' Row '',I5,'' :'',5D12.4,'
     '      //'/:(12X,5D12.4))') irow,
     '      (T_BH(irow,icol),icol=1,ISIZE_GQ_HH)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

      IF(CHECK) THEN
C***    Checking  T_bh*phih-phib
        irow=0
        SUM2=0.0d0
        DO nolist2=1,NPLIST4(0) !list of body nodes
          np2=NPLIST4(nolist2)
          DO nhx2=1,NHP(np2,nr1,nx)
            nh2=NH_LOC(nhx2,nx)
            DO nv2=1,NVHP(nh2,np2,1,nr1)
              DO nk2=1,MAX(NKH(nh2,np2,1,nr1)-KTYP93(1,nr1),1)
                ny3=NYNP(nk2,nv2,nh2,np2,0,1,nr1)
C               global ny# for this var
                irow=irow+1
                icol=0
                SUM=0.0d0
                DO nolist1=1,NPLIST3(0) !list of heart nodes
                  np1=NPLIST3(nolist1)
                  DO nhx1=1,NHP(np1,nr0,nx)
                    nh1=NH_LOC(nhx1,nx)
                    DO nv1=1,NVHP(nh1,np1,1,nr0)
                      DO nk1=1,
     '                  MAX(NKH(nh1,np1,1,nr0)-KTYP93(1,nr0),1)
                        ny1=NYNP(nk1,nv1,nh1,np1,0,1,nr0)
                        !heart column number of GK
                        icol=icol+1
                        SUM=SUM+T_BH(irow,icol)*YP(ny1,1,nx)
C                       t_bh*phih
                      ENDDO !nk1
                    ENDDO !nv1
                  ENDDO !nh1
                ENDDO !np1
                SUM=SUM-YP(ny3,1,nx) !-phib
                SUM2=SUM2+SUM*SUM
              ENDDO !nk2
            ENDDO !nv2
          ENDDO !nh2
        ENDDO !np2
        SUM2=DSQRT(SUM2)
        WRITE(OP_STRING,'('' Norm of Tbh*phih-phib ='',D12.4)')SUM2
       CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF !check

      IF(SVD) THEN !find singular values and condition number of T_BH
        DO irow=1,ISIZE_GK_BB
          DO icol=1,ISIZE_GQ_HH
            WK1_INV(irow,icol)=T_BH(irow,icol)
            WK2_INV(icol,irow)=T_BH(irow,icol) !transpose of T_bh
          ENDDO
        ENDDO
        IFAIL=1
        LWORK=MAX(3*MIN(ISIZE_GK_BB,ISIZE_GK_HH)+
     '    MAX(ISIZE_GK_BB,ISIZE_GK_HH),
     '    5*MIN(ISIZE_GK_BB,ISIZE_GK_HH)-4)
        CALL ASSERT(LWORK.LE.NY_TRANSFER_M*NY_TRANSFER_M,
     '    '>>Work array too small for SVD of T_BH',ERROR,*9999)

        CALL DGESVD('N','N',ISIZE_GK_BB,ISIZE_GK_HH,WK1_INV,
     '    NY_TRANSFER_M,WK4_INV,Q,1,P,1,WK3_INV,
     '    NY_TRANSFER_M*NY_TRANSFER_M,IFAIL)

        IF(IFAIL.NE.0)THEN
          WRITE(OP_STRING,*)' IFAIL=',IFAIL,' in SVD of T_BH'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ERROR='>>IFAIL<>0 in DGESVD'
          GOTO 9999
        ENDIF


        isize=MIN(ISIZE_GK_BB,ISIZE_GK_HH) !the last computed  SV
        IF(WK4_INV(isize).GT.RDELTA)THEN
          CONDITION_NUMBER=WK4_INV(1)/WK4_INV(isize)
        ELSE
          CONDITION_NUMBER=RMAX !default value
        ENDIF


        IF(OPFILE)THEN
          WRITE(OP_STRING,'(/'' Condition number of T_bh :'',D12.4)')
     '      CONDITION_NUMBER
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C LKC 4-AMY-1999 if large problems exceeds OP_STRING array
C so write out the first 25 lines of 5
          IF(INT(ISIZE_GK_HH/5).LT.25) THEN
            WRITE(OP_STRING,'('' Singular values of T_bh'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(13X,5D12.4,'
     '        //'/:(13X,5D12.4))')
     '        (WK4_INV(icol),icol=1,ISIZE_GK_HH)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE
            WRITE(OP_STRING,'('' First 125 of '',I8,
     '        '' Singular values of T_bh'')') ISIZE_GK_HH
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(13X,5D12.4,'
     '        //'/:(13X,5D12.4))')
     '        (WK4_INV(icol),icol=1,100)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
        WRITE(OP_STRING,'(/'' Condition number of T_bh :'',D12.4)')
     '    CONDITION_NUMBER
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

C LKC 4-AMY-1999 if large problems exceeds OP_STRING array
C so write out the first 25 lines of 5
        IF(INT(ISIZE_GK_HH/5).LT.25) THEN
          WRITE(OP_STRING,'('' Singular values of T_bh'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(13X,5D12.4,'
     '      //'/:(13X,5D12.4))')
     '      (WK4_INV(icol),icol=1,ISIZE_GK_HH)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ELSE
          WRITE(OP_STRING,'('' First 125 of '',I8,
     '      '' Singular values of T_bh'')') ISIZE_GK_HH
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(13X,5D12.4,'
     '      //'/:(13X,5D12.4))')
     '      (WK4_INV(icol),icol=1,100)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF

C***    Work out T_bh^t*T_bh and the condition number of
C***    this matrix (important in inverse solutions).
        DO irow=1,ISIZE_GK_HH
          DO icol=1,ISIZE_GK_HH
            SUM=0.0d0
            DO isize=1,ISIZE_GK_BB
              SUM=SUM+WK2_INV(irow,isize)*T_BH(isize,icol)
            ENDDO
            WK3_INV(irow,icol)=SUM !T_BH^t*T_BH
          ENDDO
        ENDDO
        IFAIL=1
        CALL DGESVD('N','N',ISIZE_GK_HH,ISIZE_GK_HH,WK3_INV,
     '    NY_TRANSFER_M,WK4_INV,Q,1,P,1,WK1_INV,
     '    NY_TRANSFER_M*NY_TRANSFER_M,IFAIL)

        IF(IFAIL.NE.0)THEN
          WRITE(OP_STRING,*)' IFAIL=',IFAIL,
     '      ' in SVD of T_BH^t*T_BH'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ERROR='>>IFAIL<>0 in DGESVD'
          GOTO 9999
        ENDIF

        isize=MIN(ISIZE_GK_BB,ISIZE_GK_HH) !the last computed  SV
        IF(WK4_INV(isize).GT.RDELTA)THEN
          CONDITION_NUMBER=WK4_INV(1)/WK4_INV(isize)
        ELSE
          CONDITION_NUMBER=RMAX !default value
        ENDIF


        IF(OPFILE)THEN
          WRITE(OP_STRING,'(/'' Condition number of T_bh^t*T_bh :'
     '       //''',D12.4)')CONDITION_NUMBER
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C LKC 4-AMY-1999 if large problems exceeds OP_STRING array
C so write out the first 25 lines of 5
          IF(INT(ISIZE_GK_HH/5).LT.25) THEN
            WRITE(OP_STRING,'('' Singular values of T_bh^t*T_bh'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(13X,5D12.4,'
     '        //'/:(13X,5D12.4))')
     '        (WK4_INV(icol),icol=1,ISIZE_GK_HH)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE
            WRITE(OP_STRING,'('' First 125 of '',I8,
     '        '' Singular values of T_bh^t*T_bh'')') ISIZE_GK_HH
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(13X,5D12.4,'
     '        //'/:(13X,5D12.4))')
     '        (WK4_INV(icol),icol=1,100)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
        WRITE(OP_STRING,'(/'' Condition number of T_bh^t*T_bh :'
     '     //''',D12.4)')CONDITION_NUMBER
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

C LKC 4-AMY-1999 if large problems exceeds OP_STRING array
C so write out the first 25 lines of 5
        IF(INT(ISIZE_GK_HH/5).LT.25) THEN
          WRITE(OP_STRING,'('' Singular values of T_bh^t*T_bh'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(13X,5D12.4,'
     '      //'/:(13X,5D12.4))')
     '      (WK4_INV(icol),icol=1,ISIZE_GK_HH)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ELSE
          WRITE(OP_STRING,'('' First 125 of '',I8,
     '      '' Singular values of T_bh^t*T_bh'')') ISIZE_GK_HH
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(13X,5D12.4,'
     '      //'/:(13X,5D12.4))')
     '      (WK4_INV(icol),icol=1,100)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF ! ISIZE_GK_HH

      ENDIF !SVD

      CALL EXITS('EVTRSF_DYNAM')
      RETURN
 9999 CALL ERRORS('EVTRSF_DYNAM',ERROR)
      CALL EXITS('EVTRSF_DYNAM')
      RETURN 1
      END


