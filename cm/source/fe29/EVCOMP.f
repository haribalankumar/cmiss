      SUBROUTINE EVCOMP(IBT,IDO,INP,LD,NAN,NBH,NBJ,NBJF,NDDATA,NEELEM,
     &  NEP,NFF,NFFACE,NFLIST,NHE,NHP,NKEF,NKH,NKHE,NKJE,NLL,NLNO,NMNO,
     &  NNB,NNF,NNL,NONL,NONM,NONY,NP1OPT,NPF,NP_INTERFACE,NPL,NPLIST,
     &  NPNE,NPNF,NPNODE,NPNY,NRE,NRLIST,NSB,NVHE,NVHP,NVJE,NVJF,NVJL,
     &  NVJP,NW,NXI,NXLIST,NYNE,NYNO,NYNP,NYNR,PAOPTY,AQ,CE,
     &  CELL_RCQS_VALUE,CG,CGE,CONY,CP,CURVCORRECT,DF,FEXT,PG,RG,SE,SF,
     &  WG,XA,XE,XG,XID,XIG,XIP,XN,XP,YG,YGF,YP,ZA,ZA1,ZD,ZE,ZG,ZG1,ZP,
     &  ZP1,STRING,FIX,ERROR,*)

C#### Subroutine: EVCOMP
C###  Description:
C###    <HTML>
C###    Evaluates outcomes from deformation of a compressible tissue.
C###    </HTML>

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     &  LD(NDM),NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     &  NBJF(NJM,NFM),NDDATA(0:NDM,0:NRM),NEELEM(0:NE_R_M,0:NRM),
     &  NEP(NPM),NFF(6,NEM),NFFACE(0:NF_R_M,NRM),NFLIST(0:NFM),
     &  NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),NKEF(0:4,16,6,NBFM),
     &  NKH(NHM,NPM,NCM,0:NRM),NKHE(NKM,NNM,NHM,NEM),
     &  NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),NLNO(NOPM,NXM),
     &  NMNO(1:2,0:NOPM,NXM),NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),
     &  NNL(0:4,12,NBFM),NONL(NLM,NXM),NONM(NMM,NPM,NXM),
     &  NONY(0:NOYM,NYM,NRCM,0:NRM,NXM),NP_INTERFACE(0:NPM,0:3),
     &  NP1OPT(NOPM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),NPL(5,0:3,NLM),
     &  NPLIST(0:NPM),NPNF(NNM,NBFM),NPNODE(0:NP_R_M,0:NRM),
     &  NPNY(0:6,NYM,0:NRCM,NXM),NRE(NEM),NSB(NKM,NNM,NBFM),
     &  NRLIST(0:NRM),NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     &  NVJE(NNM,NBFM,NJM,NEM),NVJF(NNM,NBFM,NJM),NVJL(4,NJM,NLM),
     &  NVJP(NJM,NPM),NW(NEM,3,NXM),NXI(-NIM:NIM,0:NEIM,0:NEM),
     &  NXLIST(0:NXM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     &  NYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),
     &  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     &  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),PAOPTY(NOPM)
      REAL*8 AQ(NMAQM,NQM),CE(NMM,NEM,NXM),CELL_RCQS_VALUE(NQRM,NQVM),
     '  CG(NMM,NGM),CGE(NMM,NGM,NEM,NXM),
     '  CONY(0:NOYM,NYM,NRCM,0:NRM,NXM),CP(NMM,NPM,NXM),
     '  CURVCORRECT(2,2,NNM,NEM),DF(NFM),FEXT(NIFEXTM,NGM,NEM),
     &  PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM,NEM),SF(NSM,NBFM),
     &  WG(NGM,NBM),XA(NAM,NJM),XE(NSM,NJM),XG(NJM,NUM),XID(NIM,NDM),
     &  XIG(NIM,NGM,NBM),XIP(NIM,NPM),XN(NJM,NGM),XP(NKM,NVM,NJM,NPM),
     &  YG(NIYGM,NGM,NEM),YGF(NIYGFM,NGFM,NFM),YP(NYM,NIYM,NXM),
     &  ZA(NAM,NHM,NCM,NEM),ZA1(NAM,NHM,NCM,NEM),ZD(NJM,NDM),
     &  ZE(NSM,NHM),ZG(NHM,NUM),ZG1(NHM,NUM),ZP(NKM,NVM,NHM,NPM,NCM),
     &  ZP1(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM)

!     Local Variables
      INTEGER IBEG,IEND,N3CO,nb,nd,ne,ng,
     &  nj_field,njj_field,nn,noelem,nonode,np,nr,
     &  nr_slave,nxc
      REAL*8 AZ,AZL(3,3),AZU(3,3),PHI(3),
     &  PST(3),RG2D,RGZ,RGZ2D,RM(3,3),TC(3,3),TG(3,3),TN(3,3),
     &  TNA,XI(3)
      CHARACTER COORDS*9,STRESSTYPE*17,TYPE*10,STRING2*255
      LOGICAL ALL_REGIONS,CBBREV,COMPL,ELASTICITY,FIX_NODES,
     &  LIMIT,RATIO,RECOIL,UNIFORM
! Functions
      INTEGER IFROMC
      REAL*8 DET

      CALL ENTERS('EVCOMP',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM evaluate compressible
C###  Parameter:        <recoil/compliance/elastance/ratio/limit>
C###    Specify which parameter will be evaluated. 'Recoil' is for the
C###    current elastic recoil pressure; 'compliance' is an
C###    approximation of the current tissue compliance; 'elastance' is 1
C###    /compliance; 'ratio' is the ratio of deformed to undeformed volume;
C###    'limit' is for setting zero displacement BCs.
C###  Parameter:        for recoil: <uniform/mean [uniform]>
C###    Specify whether the recoil pressure is evaluated using the cube
C###    root of the volume change ratio (uniform), or as the mean of the
C###    principal components (mean).
C###  Parameter:        for limit: <group_name>
C###    Specify the node group that the fixed nodes will be added to.
C###  Parameter:        <nodes/data>
C###    Specify whether evaluated at nodes or datapoints. When
C###    evaluating at nodes, the nodes are in a different region from
C###    the compressible mechanics solution.
C###  Parameter:        for nodes: <in (region) [1]>
C###    Specify which region the nodes are in where the solution will be
C###    evaluated.
C###  Parameter:        <field (#) [1]>
C###    Specify the field number for storing the evaluation.
C###  Parameter:        <region (#s/all)[1]>
C###    Specify the region number of the solution.
C###  Parameter:        <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//
     &    '<recoil/compliance/elastance/ratio [recoil]>'
        OP_STRING(3)=BLANK(1:15)//'for recoil: <uniform/mean [uniform]>'
        OP_STRING(4)=BLANK(1:15)//'<nodes/data [nodes]>'
        OP_STRING(5)=BLANK(1:15)//'for nodes: <in (region) [1]>'
        OP_STRING(6)=BLANK(1:15)//'<field (#) [1]>'
        OP_STRING(7)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(8)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','EVCOMP',ERROR,*9999)
      ELSE

        COMPL=.FALSE.
        ELASTICITY=.FALSE.
        RATIO=.FALSE.
        RECOIL=.FALSE.
        UNIFORM=.TRUE.
        LIMIT=.FALSE.
        
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        nr=NRLIST(1)

        IF(CBBREV(CO,'COMPLIANCE',4,noco+1,NTCO,N3CO)) COMPL=.TRUE.
        IF(CBBREV(CO,'ELASTICITY',3,noco+1,NTCO,N3CO)) ELASTICITY=.TRUE.
        IF(CBBREV(CO,'RATIO',3,noco+1,NTCO,N3CO))      RATIO=.TRUE.
        IF(CBBREV(CO,'RECOIL',3,noco+1,NTCO,N3CO))     RECOIL=.TRUE.
        IF(CBBREV(CO,'LIMIT',3,noco+1,NTCO,N3CO))      LIMIT=.TRUE.
        IF(CBBREV(CO,'UNIFORM',3,noco+1,NTCO,N3CO))THEN
          UNIFORM= .TRUE.
        ELSE IF(CBBREV(CO,'MEAN',3,noco+1,NTCO,N3CO))THEN
          UNIFORM= .FALSE.
        ENDIF
          IF(CBBREV(CO,'NOPENALTY',3,noco+1,NTCO,N3CO))THEN
           NO_PENALTY= .TRUE.
          ENDIF
 

        IF(.NOT.COMPL.AND..NOT.ELASTICITY.AND..NOT.RATIO.AND.
     &    .NOT.RECOIL.AND..NOT.LIMIT) RECOIL=.TRUE.

        IF(CBBREV(CO,'DATA',3,noco+1,NTCO,N3CO)) THEN
          TYPE='DATA'
          IF(CBBREV(CO,'FIELD',3,noco+1,NTCO,N3CO)) THEN
            nj_field=IFROMC(CO(N3CO+1))+NJT
          ELSE
            nj_field=1+NJT
          ENDIF
          CALL ASSERT(nj_field.GT.0,'>>Field not defined',ERROR,*9999)
        ELSE !nodes by default
c        IF(CBBREV(CO,'NODES',3,noco+1,NTCO,N3CO)) THEN
          TYPE='NODES'
          IF(CBBREV(CO,'IN',2,noco+1,NTCO,N3CO)) THEN
            nr_slave=nr
            nr=IFROMC(CO(N3CO+1))
          ELSE
            nr_slave=nr
            nr=nr_slave+1
          ENDIF
          IF(CBBREV(CO,'FIELD',3,noco+1,NTCO,N3CO)) THEN
            njj_field=IFROMC(CO(N3CO+1))
          ELSE
            njj_field=1
          ENDIF
          nj_field=NJ_LOC(NJL_FIEL,njj_field,NRLIST(1))
          CALL ASSERT(nj_field.GT.0,'>>Field not defined',ERROR,*9999)
        ENDIF
        
        IF(.NOT.LIMIT)THEN  
        IF(TYPE(1:5).EQ.'NODES')THEN
          CALL YPZP(1,NBH,NEELEM,NHE(1,nxc),NHP(1,nr,nxc),NKH(1,1,1,nr),
     '      NPNODE,nr,NVHP(1,1,1,nr),nxc,NYNE,NYNP,YP(1,1,nxc),ZA,ZP,
     '      ERROR,*9999)
          DO nonode=1,NPNODE(0,nr_slave)
            np=NPNODE(nonode,nr_slave)
            ne=NEP(np) !host element number
            IF(ne.NE.0)THEN !i.e. if inside a solution element
              nb=NBJ(1,ne)
              CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     &          nr,NVJE(1,1,1,ne),SE(1,1,ne),XA,XE,XP,ERROR,*9999)
              CALL ZPZE(NBH(1,1,ne),1,NHE(ne,nxc),NKHE(1,1,1,ne),
     '          NPF(1,1),NPNE(1,1,ne),NRE(ne),NVHE(1,1,1,ne),
     '          NW(ne,1,nxc),nxc,CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '          ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
              CALL CPCG(1,NBH(1,1,ne),NPNE(1,1,ne),nr,nxc,
     '          CE(1,ne,nxc),CG,CGE(1,1,ne,nxc),CP,PG,ERROR,*9999)
              XI(1)=XIP(1,np) !previously calculated xi location
              XI(2)=XIP(2,np)
              XI(3)=XIP(3,np)
              COORDS='Fibre'
              STRESSTYPE='Total'
              CALL ZETX50(COORDS,'Cauchy',STRESSTYPE,IBT,IDO,INP,NAN,
     &          NBH(1,1,ne),NBJ(1,ne),0,NHE(ne,nxc),NPNE(1,1,ne),nr,ne,
     &          nxc,CE(1,ne,nxc),CG,CP,FEXT(1,1,ne),PG,PHI,PST,RG(1),
     &          RG2D,RGZ,RGZ2D,RM,TC,TG,TN,TNA,XE,XG,XI,YG(1,1,ne),
     &          ZE,ZG,ERROR,*9999)
              CALL ZGMG(nb,nr,AZ,AZL,AZU,ZG,ERROR,*9999)

              IF(RECOIL)THEN !pressure w.r.t. deformed geometry
                CALL COMPRESSIBLE_HYDROSTATIC_PRESSURE(XP(1,1,nj_field,
     &            np),CG(1,1),TC,DSQRT(DET(AZL)),UNIFORM,ERROR,*9999)
              ELSE IF(RATIO)THEN !evaluate ratio of deformed to undeformed
                XP(1,1,nj_field,np)=DSQRT(DET(AZL))
              ELSE IF(COMPL)THEN !evaluate C/V0=1/(dP/dV*V0)
                CALL COMPRESSIBLE_HYDROSTATIC_ELASTICITY(XP(1,1,
     &            nj_field,np),CG(1,1),TC,DSQRT(DET(AZL)),ERROR,*9999)
              ELSE IF(ELASTICITY)THEN !evaluate E*V0=dP/dV*V0
                CALL COMPRESSIBLE_HYDROSTATIC_ELASTICITY(XP(1,1,
     &            nj_field,np),CG(1,1),TC,DSQRT(DET(AZL)),ERROR,*9999)
                XP(1,1,nj_field,np)=1.d0/XP(1,1,nj_field,np)
              ENDIF

            ELSE
              !node not internal
              XP(1,1,nj_field,np)=0.d0
            ENDIF
          ENDDO !nonode
        ELSE IF(TYPE(1:4).EQ.'DATA')THEN
          DO nd=1,NDT
            ne=LD(nd)
            IF(ne.EQ.0)THEN
              ZD(nj_field,nd)=0.d0
            ELSE
              nb=NBJ(1,ne)
              CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     &          nr,NVJE(1,1,1,ne),SE(1,1,ne),XA,XE,XP,ERROR,*9999)
              CALL ZPZE(NBH(1,1,ne),1,NHE(ne,nxc),NKHE(1,1,1,ne),
     '          NPF(1,1),NPNE(1,1,ne),NRE(ne),NVHE(1,1,1,ne),
     '          NW(ne,1,nxc),nxc,CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '          ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
              CALL CPCG(1,NBH(1,1,ne),NPNE(1,1,ne),nr,nxc,
     '          CE(1,ne,nxc),CG,CGE(1,1,ne,nxc),CP,PG,ERROR,*9999)
              XI(1)=XID(1,nd)
              XI(2)=XID(2,nd)
              XI(3)=XID(3,nd)
              COORDS='Fibre'
              STRESSTYPE='Total'

              CALL ZETX50_1(COORDS,'Cauchy',STRESSTYPE,IBT,IDO,INP,NAN,
     &          NBH(1,1,ne),NBJ(1,ne),0,NHE(ne,nxc),NPNE(1,1,ne),nr,ne,
     &          nxc,CE(1,ne,nxc),CG,CP,FEXT(1,1,ne),PG,PHI,PST,RG(1),
     &          RG2D,RGZ,RGZ2D,RM,TC,TG,TN,TNA,XE,XG,XI,YG,ZE,ZG,
     &          ERROR,*9999)

              CALL ZGMG(nb,nr,AZ,AZL,AZU,ZG,ERROR,*9999)

              IF(RECOIL)THEN !pressure w.r.t. deformed geometry
                CALL COMPRESSIBLE_HYDROSTATIC_PRESSURE(ZD(nj_field,nd),
     &            CG(1,1),TC,DSQRT(DET(AZL)),UNIFORM,ERROR,*9999)
              ELSE IF(RATIO)THEN !evaluate ratio of deformed to undeformed
c                ZD(nj_field,nd)=MAX(DSQRT(DET(AZL)),1.d0)
                ZD(nj_field,nd)=DSQRT(DET(AZL))
              ELSE IF(COMPL)THEN !evaluate C/V0=1/(dP/dV*V0)
                CALL COMPRESSIBLE_HYDROSTATIC_ELASTICITY(ZD(nj_field,
     &            nd),CG(1,1),TC,DSQRT(DET(AZL)),ERROR,*9999)
                ZD(nj_field,nd)=1.d0/ZD(nj_field,nd)
              ELSE IF(ELASTICITY)THEN !evaluate E*V0=dP/dV*V0
                CALL COMPRESSIBLE_HYDROSTATIC_ELASTICITY(ZD(nj_field,
     &            nd),CG(1,1),TC,DSQRT(DET(AZL)),ERROR,*9999)
              ENDIF
            ENDIF
          ENDDO !nonode
        ENDIF
        ELSEIF(LIMIT)THEN
           NPLIST(0) = 0

      DO noelem=1,NEELEM(0,nr)
         ne=NEELEM(noelem,nr)
         nb=NBJ(1,ne)
         FIX_NODES=.FALSE.
         CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     &     nr,NVJE(1,1,1,ne),SE(1,1,ne),XA,XE,XP,ERROR,*9999)
         CALL ZPZE(NBH(1,1,ne),1,NHE(ne,nxc),NKHE(1,1,1,ne),
     '     NPF(1,1),NPNE(1,1,ne),NRE(ne),NVHE(1,1,1,ne),
     '     NW(ne,1,nxc),nxc,CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '     ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
         CALL CPCG(1,NBH(1,1,ne),NPNE(1,1,ne),nr,nxc,
     '     CE(1,ne,nxc),CG,CGE(1,1,ne,nxc),CP,PG,ERROR,*9999)
           DO ng=1,NGT(nb)
             COORDS='Fibre'
             STRESSTYPE='Total'
c             CALL ZETX50(COORDS,'Cauchy',STRESSTYPE,IBT,IDO,INP,NAN,
c     &         NBH(1,1,ne),NBJ(1,ne),ng,NHE(ne,nxc),NPNE(1,1,ne),nr,
c     &         nxc,CE(1,ne,nxc),CG,CP,FEXT(1,ng,ne),PG,PHI,PST,
c     &         RG(ng),RG2D,RGZ,RGZ2D,RM,TC,TG,TN,TNA,XE,XG,XIG(1,ng,
c     &         nb),YG(1,ng,ne),ZE,ZG,ERROR,*9999)
              CALL ZETX50(COORDS,'Cauchy',STRESSTYPE,IBT,IDO,INP,NAN,
     &          NBH(1,1,ne),NBJ(1,ne),ng,NHE(ne,nxc),NPNE(1,1,ne),nr,ne,
     &          nxc,CE(1,ne,nxc),CG,CP,FEXT(1,ng,ne),PG,PHI,PST,RG(1),
     &          RG2D,RGZ,RGZ2D,RM,TC,TG,TN,TNA,XE,XG,XIG(1,ng,nb),
     &          YG(1,ng,ne),ZE,ZG,ERROR,*9999)
             CALL ZGMG(nb,nr,AZ,AZL,AZU,ZG,ERROR,*9999)
             YG(10,ng,ne)=DSQRT(DET(AZL)) !ratio deformed to undeformed
             IF(YG(10,ng,ne).LE.1.d0)THEN
C                write(*,*)'yg',ne,ng,YG(10,ng,ne)
                FIX_NODES=.TRUE.
             ENDIF
           ENDDO  !ng
           IF(FIX_NODES)THEN
c              WRITE(*,*) 'FIX NODES',ne
              DO nn=1,8
                 NPLIST(0)=NPLIST(0)+1
                 NPLIST(NPLIST(0))=NPNE(nn,nb,ne)
              ENDDO
           ENDIF
        ENDDO  !noelem

        IF(NPLIST(0).GT.0)THEN
          STRING2='FIXED_NODES'
C         Sort and remove duplicates from the node list
          CALL ILISTRMDUP(NPLIST(0),NPLIST(1),ERROR,*9999)
          CALL GRNODE_SUB(NPLIST,STRING2,
     &      .TRUE.,ERROR,*9999)
       ENDIF

        ENDIF
      ENDIF

      CALL EXITS('EVCOMP')
      RETURN
 9999 CALL ERRORS('EVCOMP',ERROR)
      CALL EXITS('EVCOMP')
      RETURN 1
      END



