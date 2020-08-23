      SUBROUTINE EVPRESS(IBT,IDO,INP,LD,NAN,NBH,NBJ,NBJF,NDDATA,NEELEM,
     &  NEP,NFF,NFFACE,NFLIST,NHE,NHP,NKEF,NKH,NKHE,NKJE,NLL,NLNO,NMNO,
     &  NNB,NNF,NNL,NONL,NONM,NONY,NP1OPT,NPF,NP_INTERFACE,NPL,NPLIST,
     &  NPNE,NPNF,NPNODE,NPNY,NRE,NRLIST,NSB,NVHE,NVHP,NVJE,NVJF,NVJL,
     &  NVJP,NW,NXI,NXLIST,NYNE,NYNO,NYNP,NYNR,PAOPTY,AQ,CE,
     &  CELL_RCQS_VALUE,CG,CGE,CONY,CP,CURVCORRECT,DF,FEXT,PG,RG,SE,SF,
     &  WG,XA,XE,XG,XID,XIG,XIP,XN,XP,YG,YGF,YP,ZA,ZA1,ZD,ZE,ZG,ZG1,ZP,
     &  ZP1,STRING,FIX,ERROR,*)

C#### Subroutine: EVPRESS
C###  Description:
C###    <HTML>
C###    EVPRESS evaluates effective pressures at nodes for a
C###    finite deformation elasticity problem. This has been implemented 
C###    for compressible mechanics, where an effective hydrostatic
C###    pressure is developed as the tissue changes in volume.<BR>
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
      INTEGER IBEG,IEND,i,n1,N3CO,nb,nd,ndd,ne,nef,nf,ng,
     &  nj_field,nj_vol,njj_field,nn,nne,noelem,nonode,np,nr,
     &  nr_slave,nv,nxc
      REAL*8 AZ,AZL(3,3),AZU(3,3),const,EVAL(3),EVEC(3,3),PHI(3),
     &  PST(3),RG2D,RGZ,RGZ2D,RM(3,3),TC(3,3),TC1,TC2,TG(3,3),TN(3,3),
     &  TNA,XI_SPOT,XI(3)
      CHARACTER COORDS*9,STRESSTYPE*17,TYPE*10
      LOGICAL ALL_REGIONS,CBBREV,INITIAL,SURFACES,UNIFORM,
     &  XI_POINT,XI_DATA_POINT
! Functions
      INTEGER IFROMC
      LOGICAL INLIST
      REAL*8 DET,RFROMC

      CALL ENTERS('EVPRESS',*9999)
      
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM evaluate pressures
C###  Parameter:        <surface/nodes/data/gauss>
C###    Specify whether pressures calculated on surfaces, at nodes,
C###    datapoints, or gauss points .
C###  Parameter:        <faces (#s/groups/all)[all]>
C###    Specify the faces on which to evaluate pressures.
C###  Parameter:        <uniform/averaged [averaged]>
C###    Specify the faces on which to evaluate pressures.
C###  Parameter:        <region (#s/all)[1]>
C###    Specify the region numbers to evaluate.
C###  Parameter:        <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Evaluates pressures at nodes.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<surface/nodes/data/gauss [nodes]>'
        OP_STRING(3)=BLANK(1:15)//'<faces (#s/groups/all)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<uniform/averaged [averaged]>'
        OP_STRING(5)=BLANK(1:15)//'<field (#) [1]>'
        OP_STRING(6)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(7)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','EVPRESS',ERROR,*9999)
      ELSE

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_FACES(NFFACE,NFLIST,noco,NRLIST,NTCO,CO,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        nr=NRLIST(1)

        IF(CBBREV(CO,'UNIFORM',3,noco+1,NTCO,N3CO)) THEN
          UNIFORM=.TRUE.
        ELSE
          UNIFORM = .FALSE.
        ENDIF

        XI_POINT=.TRUE. !default

        IF(CBBREV(CO,'NODES',3,noco+1,NTCO,N3CO)) THEN
          TYPE='NODES'
C          XI_POINT=.TRUE.
C          XI_DATA_POINT=.FALSE.
C          SURFACES=.FALSE.
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
        ELSE IF(CBBREV(CO,'DATA',3,noco+1,NTCO,N3CO)) THEN
          TYPE='DATA'
C          XI_POINT=.FALSE.
C          XI_DATA_POINT=.TRUE.
C          SURFACES=.FALSE.
          IF(CBBREV(CO,'FIELD',3,noco+1,NTCO,N3CO)) THEN
            nj_field=IFROMC(CO(N3CO+1))+NJT
          ELSE
            nj_field=1+NJT
          ENDIF
          IF(CBBREV(CO,'VOLUME',3,noco+1,NTCO,N3CO)) THEN
            nj_vol=NJ_LOC(NJL_FIEL,IFROMC(CO(N3CO+1)),NRLIST(1))
          ELSE
            nj_vol=1+NJT
          ENDIF
          CALL ASSERT(nj_field.GT.0,'>>Field not defined',ERROR,*9999)
          
        ELSE IF(CBBREV(CO,'GAUSS',3,noco+1,NTCO,N3CO)) THEN
          TYPE='GAUSS'
          IF(CBBREV(CO,'INITIAL',3,noco+1,NTCO,N3CO)) THEN
            INITIAL=.TRUE.
            const=RFROMC(CO(N3CO+1))
          ELSE
            INITIAL=.FALSE.
          ENDIF
        ELSE IF(CBBREV(CO,'SURFACE',3,noco+1,NTCO,N3CO)) THEN
          TYPE='SURFACE'
C          XI_POINT=.FALSE.
C          SURFACES=.TRUE.
          IF(CBBREV(CO,'AT',1,noco+1,NTCO,N3CO)) THEN
            XI_SPOT=RFROMC(CO(N3CO+1))
          ENDIF
          IF(CBBREV(CO,'FIELD',3,noco+1,NTCO,N3CO)) THEN
            njj_field=IFROMC(CO(N3CO+1))
          ELSE
            njj_field=1
          ENDIF
          nj_field=NJ_LOC(NJL_FIEL,njj_field,NRLIST(1))
          CALL ASSERT(nj_field.GT.0,'>>Field not defined',ERROR,*9999)
        ENDIF
        

        IF(TYPE(1:7).EQ.'SURFACE')THEN
          NPLIST(0)=0
          DO i=1,NPM
            NPLIST(i)=0
          ENDDO !i
          CALL YPZP(1,NBH,NEELEM,NHE(1,nxc),NHP(1,nr,nxc),NKH(1,1,1,nr),
     '      NPNODE,nr,NVHP(1,1,1,nr),nxc,NYNE,NYNP,YP(1,1,nxc),ZA,ZP,
     '      ERROR,*9999)
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            nb=NBJ(1,ne)
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     &        nr,NVJE(1,1,1,ne),SE(1,1,ne),XA,XE,XP,ERROR,*9999)
            CALL ZPZE(NBH(1,1,ne),1,NHE(ne,nxc),NKHE(1,1,1,ne),
     '        NPF(1,1),NPNE(1,1,ne),NRE(ne),NVHE(1,1,1,ne),
     '        NW(ne,1,nxc),nxc,CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '        ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
            CALL CPCG(1,NBH(1,1,ne),NPNE(1,1,ne),nr,nxc,
     '        CE(1,ne,nxc),CG,CGE(1,1,ne,nxc),CP,PG,ERROR,*9999)
            DO nef=1,6
              nf=NFF(nef,ne) !global face number
              IF(nf.NE.0)THEN
                IF(NPF(5,nf).EQ.1)THEN !in single element, boundary
                  DO nne=1,NNF(0,nef,nb)
                    nn=NNF(1+nne,nef,nb)
                    np=NPNE(nn,nb,ne)
 !check to see whether already done
                    IF(.NOT.INLIST(np,NPLIST(1),MIN(NPLIST(0),NPM),
     &                n1))THEN
                      NPLIST(0)=NPLIST(0)+1
                      IF(NPLIST(0).LE.NPM) NPLIST(NPLIST(0))=np
                      XI(1)=1.d0-XI_SPOT
                      XI(2)=1.d0-XI_SPOT
                      XI(3)=1.d0-XI_SPOT
                      IF(nn.EQ.2.OR.nn.EQ.4.OR.nn.EQ.6.OR.nn.EQ.8)
     &                  XI(1)=XI_SPOT
                      IF(nn.EQ.3.OR.nn.EQ.4.OR.nn.EQ.7.OR.nn.EQ.8)
     &                  XI(2)=XI_SPOT
                      IF(nn.GT.4) XI(3)=XI_SPOT
                      COORDS='Fibre'
                      STRESSTYPE='Total'
                      
                      CALL ZETX50(COORDS,'Cauchy',STRESSTYPE,IBT,IDO,
     &                  INP,NAN,NBH(1,1,ne),NBJ(1,ne),0,NHE(ne,nxc),
     &                  NPNE(1,1,ne),nr,ne,nxc,CE(1,ne,nxc),CG,CP,
     &                  FEXT(1,1,ne),PG,PHI,PST,RG(1),RG2D,RGZ,RGZ2D,RM,
     &                  TC,TG,TN,TNA,XE,XG,XI,YG(1,1,ne),ZE,ZG,ERROR,
     &                  *9999)
                      CALL ZGMG(nb,nr,AZ,AZL,AZU,ZG,ERROR,*9999)
C Use Cauchy stress tensor to get pressure w.r.t. deformed geometry
                      DO nv=1,NVJP(nj_field,np)
                        CALL COMPRESSIBLE_HYDROSTATIC_PRESSURE(
     &                    XP(1,nv,nj_field,np),CG(1,1),TC,
     &                    DSQRT(DET(AZL)),UNIFORM,ERROR,*9999)
                      ENDDO !nv
                    ENDIF !INLIST
                  ENDDO !nne
                ENDIF !NPF
              ENDIF !nf.NE.0
            ENDDO !nef
          ENDDO !noelem
        ELSE IF(TYPE(1:5).EQ.'NODES')THEN
          CALL YPZP(1,NBH,NEELEM,NHE(1,nxc),NHP(1,nr,nxc),NKH(1,1,1,nr),
     '      NPNODE,nr,NVHP(1,1,1,nr),nxc,NYNE,NYNP,YP(1,1,nxc),ZA,ZP,
     '      ERROR,*9999)
          DO nonode=1,NPNODE(0,nr_slave)
            np=NPNODE(nonode,nr_slave)
            ne=NEP(np)
            IF(ne.NE.0)THEN
              nb=NBJ(1,ne)
              CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     &          nr,NVJE(1,1,1,ne),SE(1,1,ne),XA,XE,XP,ERROR,*9999)
              CALL ZPZE(NBH(1,1,ne),1,NHE(ne,nxc),NKHE(1,1,1,ne),
     '          NPF(1,1),NPNE(1,1,ne),NRE(ne),NVHE(1,1,1,ne),
     '          NW(ne,1,nxc),nxc,CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '          ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
              CALL CPCG(1,NBH(1,1,ne),NPNE(1,1,ne),nr,nxc,
     '          CE(1,ne,nxc),CG,CGE(1,1,ne,nxc),CP,PG,ERROR,*9999)
              XI(1)=XIP(1,np)
              XI(2)=XIP(2,np)
              XI(3)=XIP(3,np)
              COORDS='Fibre'
              STRESSTYPE='Total'
              CALL ZETX50(COORDS,'Cauchy',STRESSTYPE,IBT,IDO,INP,NAN,
     &          NBH(1,1,ne),NBJ(1,ne),0,NHE(ne,nxc),NPNE(1,1,ne),nr,ne,
     &          nxc,CE(1,ne,nxc),CG,CP,FEXT(1,1,ne),PG,PHI,PST,RG(1),
     &          RG2D,RGZ,RGZ2D,RM,TC,TG,TN,TNA,XE,XG,XI,YG(1,1,ne),ZE,
     &          ZG,ERROR,*9999)
              CALL ZGMG(nb,nr,AZ,AZL,AZU,ZG,ERROR,*9999)
C Use Cauchy stress tensor to get pressure w.r.t. deformed geometry

              CALL COMPRESSIBLE_HYDROSTATIC_PRESSURE(XP(1,1,nj_field,
     &          np),CG(1,1),TC,DSQRT(DET(AZL)),UNIFORM,ERROR,*9999)

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
              ZD(nj_field+1,nd)=0.d0
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
              CALL ZETX50(COORDS,'Cauchy',STRESSTYPE,IBT,IDO,INP,NAN,
     &          NBH(1,1,ne),NBJ(1,ne),0,NHE(ne,nxc),NPNE(1,1,ne),nr,ne,
     &          nxc,CE(1,ne,nxc),CG,CP,FEXT(1,1,ne),PG,PHI,PST,RG(1),
     &          RG2D,RGZ,RGZ2D,RM,TC,TG,TN,TNA,XE,XG,XI,YG(1,1,ne),ZE,
     &          ZG,ERROR,*9999)
C Use Cauchy stress tensor to get pressure w.r.t. deformed geometry
              CALL ZGMG(nb,nr,AZ,AZL,AZU,ZG,ERROR,*9999)
              ZD(nj_field+1,nd)=DSQRT(DET(AZL)) !V/V0
              IF(ZD(nj_field+1,nd).LE.1.d0)THEN
                ZD(nj_field+1,nd)=1.d0
              ENDIF
              CALL COMPRESSIBLE_HYDROSTATIC_PRESSURE(ZD(nj_field,nd),
     &          CG(1,1),TC,ZD(nj_field+1,nd),UNIFORM,ERROR,*9999)
c              ZD(nj_vol,nd)=ZD(nj_field+1,nd)*ZD(nj_vol,nd)
              ZD(nj_vol+1,nd)=ZD(nj_field+1,nd)*ZD(nj_vol,nd)
            ENDIF
          ENDDO !nonode
        ELSE IF(TYPE(1:5).EQ.'GAUSS')THEN
          IF(.NOT.INITIAL)THEN
             MEAN_VV_RATIO = 0.d0
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              nb=NBJ(1,ne)
              CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     &          nr,NVJE(1,1,1,ne),SE(1,1,ne),XA,XE,XP,ERROR,*9999)
              CALL ZPZE(NBH(1,1,ne),1,NHE(ne,nxc),NKHE(1,1,1,ne),
     '          NPF(1,1),NPNE(1,1,ne),NRE(ne),NVHE(1,1,1,ne),
     '          NW(ne,1,nxc),nxc,CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '          ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
              CALL CPCG(1,NBH(1,1,ne),NPNE(1,1,ne),nr,nxc,
     '          CE(1,ne,nxc),CG,CGE(1,1,ne,nxc),CP,PG,ERROR,*9999)
              DO ng=1,NGT(nb)
                COORDS='Fibre'
                STRESSTYPE='Total'
                CALL ZETX50(COORDS,'Cauchy',STRESSTYPE,IBT,IDO,INP,NAN,
     &            NBH(1,1,ne),NBJ(1,ne),ng,NHE(ne,nxc),NPNE(1,1,ne),nr,
     &            ne,nxc,CE(1,ne,nxc),CG,CP,FEXT(1,ng,ne),PG,PHI,PST,
     &            RG(ng),RG2D,RGZ,RGZ2D,RM,TC,TG,TN,TNA,XE,XG,XIG(1,ng,
     &            nb),YG(1,ng,ne),ZE,ZG,ERROR,*9999)
C Use Cauchy stress tensor to get pressure w.r.t. deformed geometry
c                YG(11,ng,ne)=(TC(1,1)+TC(2,2)+TC(3,3))/3.d0
                CALL ZGMG(nb,nr,AZ,AZL,AZU,ZG,ERROR,*9999)
                YG(10,ng,ne)=DSQRT(DET(AZL)) !ratio deformed to undeformed
                write(*,*) ng,ne,DSQRT(DET(AZL))
                CALL COMPRESSIBLE_HYDROSTATIC_PRESSURE(YG(11,ng,ne),
     &            CG(1,1),TC,DSQRT(DET(AZL)),UNIFORM,ERROR,*9999)
                YG(12,ng,ne)=YG(11,ng,ne) 
                MEAN_VV_RATIO=MEAN_VV_RATIO+YG(10,ng,ne)

              ENDDO !ng
            ENDDO !noelem
            MEAN_VV_RATIO=MEAN_VV_RATIO/DBLE(NGT(nb)*NEELEM(0,nr))
          ELSE
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              nb=NBJ(1,ne)
              DO ng=1,NGT(nb)
                YG(10,ng,ne) = 2.d0
                YG(11,ng,ne) = const
              ENDDO !ng
            ENDDO !noelem
          ENDIF
        ENDIF
        
      ENDIF

      CALL EXITS('EVPRESS')
      RETURN
 9999 CALL ERRORS('EVPRESS',ERROR)
      CALL EXITS('EVPRESS')
      RETURN 1
      END



