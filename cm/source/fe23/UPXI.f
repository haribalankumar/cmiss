      SUBROUTINE UPXI(IBT,IDO,INP,ISDANO,ISDAPR,ISDATA,ISDATR,ISEG,FD,
     '  LD,LN,MXI,NBJ,NBJF,NBH,NBHF,NDDL,NDLT,NDP,NEELEM,NELIST,NFF,
     &  NFFACE,NFLIST,NHE,NKEF,NKHE,NKJE,NNF,NPF,NPLIST,NPNE,NPNF,NRE,
     &  NRLIST,NVHE,NVJE,NVJF,NVJP,NW,NXI,CE,CG,CGE,CP,CURVCORRECT,PG,
     &  SE,SF,SQ,WD,WDL,XA,XE,XG,XID,XIDL,XP,ZA,ZD,ZDD,ZDL,ZE,ZP,CSEG,
     &  STRING,ERROR,*)
C SMAR009 19/01/99 removed NPL,
C#### Subroutine: UPXI
C###  Description:
C###    UPXI updates Xi coordinates of data points.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn' !MEM_INIT
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc' !DP_TYPE
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISDANO(NWM,NEM),ISDAPR(NWM,NEM),ISDATA(NWM,NGRSEGM),
     &  ISDATR(NWM,NEM),ISEG(*),FD(NDM),LD(NDM),LN(0:NEM),
     &  NBHF(NHM,NCM,NFM),MXI(2,NEM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NBH(NHM,NCM,NEM),NDDL(NEM,NDEM),NDLT(NEM),
     &  NDP(NDM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),NFF(6,NEM),
     &  NFFACE(0:NF_R_M,NRM),NFLIST(0:NFM),NHE(NEM,NXM),
     &  NKEF(0:4,16,6,NBFM),NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     &  NNF(0:17,6,NBFM),NPF(9,NFM),NPLIST(0:NPM),
     &  NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),
     &  NRE(NEM),NRLIST(0:NRM),NVHE(NNM,NBFM,NHM,NEM),
     &  NVJE(NNM,NBFM,NJM,NEM),NVJF(NNM,NBFM,NJM),NVJP(NJM,NPM),
     &  NW(NEM,3,NXM),NXI(-NIM:NIM,0:NEIM,0:NEM)
C SMAR009 19/01/99 removed NPL(5,0:3,NLM),
      REAL*8 CE(NMM,NEM,NXM),CG(NMM,NGM),CGE(NMM,NGM,NEM,NXM),
     '  CP(NMM,NPM,NXM),
     '  CURVCORRECT(2,2,NNM,NEM),PG(NSM,NUM,NGM,NBM),SE(NSM,NBFM,NEM),
     '  SF(NSM,NBFM),SQ(NDM),
     '  WD(NJM,NDM),WDL(NHM,NDEM),XA(NAM,NJM),XE(NSM,NJM),XG(NJM,NUM),
     '  XID(NIM,NDM),XIDL(NIM,NDEM),XP(NKM,NVM,NJM,NPM),
     '  ZA(NAM,NHM,NCM,NEM),ZD(NJM,NDM),ZDD(NJM,NDM),ZDL(NHM,NDEM),
     '  ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER from_field(3),IBEG,IEND,IFROMC,INDX,iw,L,L1,N3CO,nj,
     &  njj_field,nolist,nr,to_field(3)
      INTEGER*4 NCLOSFACE_PTR,NDP_INV_PTR,NKJF_PTR,FDSQ_PTR
      REAL*8 RFROMC,XI_1,XI_2,XI_3,ZDLMIN
      CHARACTER TITLE*80,TYPE*11
      LOGICAL ALL_REGIONS,CBBREV,DATA_UPDATE,EDGE

      CALL ENTERS('UPXI',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM update xi
C###  Parameter:      <in (ELEMENT#s/all)[all]>
C###   Specify the element numbers within which to update the xi values
C###  Description:
C###    Update Xi coordinate projections of data points in specified
C###    elements.

        OP_STRING(1)=STRING(1:IEND)//' <in (ELEMENT#s/all)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc','UPXI',ERROR,*9999)
      ELSE
        
        IF(CBBREV(CO,'POINTS',4,noco+1,NTCO,N3CO)) THEN
          DATA_UPDATE=.FALSE.
        ELSE
          DATA_UPDATE=.TRUE.
        ENDIF
        
        IF(CBBREV(CO,'IN',1,noco+1,noco+1,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NE_R_M,NELIST(0),NELIST(1),ERROR,*9999)
        ELSE
          DO L=1,LN(0)
            nolist=IABS(LN(L))
            NELIST(nolist)=nolist
          ENDDO
          NELIST(0)=LN(0)
        ENDIF

        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,*) 'NELIST: ',
     '      (NELIST(nolist),nolist=1,NELIST(0))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        IF(DATA_UPDATE)THEN
          DO nolist=1,NELIST(0)
            L1=NELIST(nolist)
            CALL NWXID(32,L1,LD,LN,NBJ,NBJF,NDDL,NDLT,NKJE,NKEF,
     '        NNF,NPF,NPNE,NPNF,NRE,NVJE,NVJF,NXI,SE,SF,SQ,
     '        XA,XE,XID,XP,ZD,ERROR,*9999)
C     SMAR009 18/01/99 removed NPL from list
          ENDDO
          TYPE(1:11)='PROJECTIONS'
          iw=2*NJT-3
          ZDLMIN=0.d0
          CALL ACWK(iw,0,ERROR,*9999)
          INDX=0
          CALL SGDATA(INDX,IBT,IDO,INP,ISDANO,ISDAPR,ISDATA(iw,NTDATA),
     '      ISDATR,ISEG,iw,LD,MXI,NBJ,NBH,NDDL,NDLT,NDP,NEELEM,NELIST,
     '      NKHE,NKJE,NPF,NPNE,NRE,NVHE,NVJE,NW,CE,CG,
     '      CGE,CP,CSEG,CURVCORRECT,PG,SE,.TRUE.,TITLE,TYPE,WD,WDL,
     '      XA,XE,XG,XID,XIDL,XP,ZA,ZD,ZDD,ZDL,ZDLMIN,ZE,ZP,ERROR,*9999)
          CALL DAWK(iw,0,ERROR,*9999)
        ELSE
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          CALL CUPPER(CO(N3CO+1)(IBEG:IEND),STRING)
c          CALL PARSILG(NDLIST,NDM,'NODES',STRING,ERROR,*9999)
          CALL PARSILG(NPLIST,NPM,'NODES',STRING,ERROR,*9999)
          IF(CBBREV(CO,'EDGES',4,noco+1,NTCO,N3CO)) THEN
            EDGE=.TRUE.
          ELSE
            EDGE=.FALSE.
          ENDIF
          IF(CBBREV(CO,'XI1',3,noco+1,NTCO,N3CO)) THEN
            XI_1=RFROMC(CO(N3CO+1))
          ELSE
            XI_1=2.d0
          ENDIF
          IF(CBBREV(CO,'XI2',3,noco+1,NTCO,N3CO)) THEN
            XI_2=RFROMC(CO(N3CO+1))
          ELSE
            XI_2=2.d0
          ENDIF
          IF(CBBREV(CO,'XI3',3,noco+1,NTCO,N3CO)) THEN
            XI_3=RFROMC(CO(N3CO+1))
          ELSE
            XI_3=2.d0
          ENDIF
          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,
     '      ALL_REGIONS,ERROR,*9999)
          nr=NRLIST(1)
          IF(CBBREV(CO,'FROM_FIELD',4,noco+1,NTCO,N3CO)) THEN
            DO nj=1,NJT
              njj_field=IFROMC(CO(N3CO+nj))
              from_field(nj)=NJ_LOC(NJL_FIEL,njj_field,nr)
            ENDDO
          ELSE
            DO nj=1,NJT
              from_field(nj)=NJ_LOC(NJL_FIEL,nj,nr)
            ENDDO
          ENDIF
          IF(CBBREV(CO,'TO_FIELD',2,noco+1,NTCO,N3CO)) THEN
            DO nj=1,NJT
              njj_field=IFROMC(CO(N3CO+nj))
              to_field(nj)=NJ_LOC(NJL_FIEL,njj_field,nr)
            ENDDO
          ELSE
            DO nj=1,NJT
              to_field(nj)=NJ_LOC(NJL_FIEL,nj+3,nr)
            ENDDO
          ENDIF
          CALL PARSE_FACES(NFFACE,NFLIST,noco-1,NRLIST,NTCO,CO,
     '      ERROR,*9999)
          NDP_INV_PTR=0
          NCLOSFACE_PTR=0
          NKJF_PTR=0
          FDSQ_PTR=0
          CALL ALLOCATE_MEMORY(NDM,1,INTTYPE,NDP_INV_PTR,
     &      MEM_INIT,ERROR,*9999)
          CALL ALLOCATE_MEMORY(NFM+1,1,INTTYPE,NCLOSFACE_PTR,
     &      MEM_INIT,ERROR,*9999)
          CALL ALLOCATE_MEMORY(NKM*NNM*NJM,1,INTTYPE,NKJF_PTR,
     &      MEM_INIT,ERROR,*9999)
          CALL ALLOCATE_MEMORY(NFM+1,1,DPTYPE,FDSQ_PTR,
     &      MEM_INIT,ERROR,*9999)
          CALL UPXI_POINTS(FD,from_field,IBT,IDO,INP,LD,NBH,NBHF,
     &      NBJ,NBJF,%VAL(NCLOSFACE_PTR),NDP,%VAL(NDP_INV_PTR),
     &      NFF,NFLIST,NHE,NKEF,NKHE,NKJE,%VAL(NKJF_PTR),NNF,
     &      NPF,NPLIST,NPNE,NPNF,NRE,NVHE,NVJE,NVJF,NVJP,1,NXI,
     &      to_field,%VAL(FDSQ_PTR),SE,SF,SQ,XA,XE,XI_1,
     &      XI_2,XI_3,XID,XP,ZD,ZE,ZP,EDGE,ERROR,*9999)
          CALL FREE_MEMORY(NDP_INV_PTR,ERROR,*9999)
          CALL FREE_MEMORY(NCLOSFACE_PTR,ERROR,*9999)
          CALL FREE_MEMORY(NKJF_PTR,ERROR,*9999)
          CALL FREE_MEMORY(FDSQ_PTR,ERROR,*9999)

        ENDIF
      ENDIF

      CALL EXITS('UPXI')
      RETURN
 9999 CALL ERRORS('UPXI',ERROR)
      CALL EXITS('UPXI')
      RETURN 1
      END



