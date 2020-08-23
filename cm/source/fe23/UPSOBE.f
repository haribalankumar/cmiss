      SUBROUTINE UPSOBE(IBT,IDO,INP,ISIZE_MFI,
     '  ISIZE_PHI,ISIZE_TBH,LD,LDR,LGE,
     '  NAN,NBH,NBHF,NBJ,NBJF,NDDL,NDLT,NEELEM,NEL,NELIST,
     '  NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,NKHE,
     '  NKJE,NLL,NLNO,NMNO,NNB,NNF,NNL,NONL,NONM,NONY,
     '  NP_INTERFACE,NP1OPT,NPF,NPL,NPLIST3,NPNE,
     '  NPNODE,NPNY,NRE,NRLIST,NSB,NVHE,NVHP,NVJE,NVJL,
     '  NW,NXI,NXLIST,NYNE,NYNO,NYNP,NYNR,PAOPTY,Z_CONT_LIST,AQ,
     '  CE,CELL_RCQS_VALUE,CG,CGE,CONY,CP,CURVCORRECT,
     '  DL,D_RE,D_RI3,
     &  D_RP,D_TG,D_ZG,ES,FEXT,FGRAD,LAPL,LAPLSQR,MFI,PAOPTI,PBOPTI,
     '  PG,PHI,PHI_H,PMIN,PMAX,
     '  RE1,RE2,RESID,RESIDM,RESJAC,RG,SE,T_BH,WG,
     '  WK1_INV,WU,XA,
     '  XE,XG,XID,XIG,XN,XP,YG,YGF,YP,ZA,ZA1,Z_CONT,ZD,ZE,ZG,
     '  ZG1,ZP,ZP1,STRING,FIX,ERROR,*)

C#### Subroutine: UPSOBE
C###  Description:
C###    UPSOBE updates Sobolev weights and scaling factors.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'data00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'opti00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISIZE_MFI(3,NSSM),ISIZE_PHI(2),ISIZE_TBH(2),LD(NDM),LDR(0:NDM),
     '  LGE(NHM*NSM),NAN(NIM,NAM,NBFM),
     '  NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NDDL(NEM,NDEM),NDLT(NEM),NEELEM(0:NE_R_M,0:NRM),NEL(0:NELM,NLM),
     '  NELIST(0:NEM),NFF(6,NEM),NFFACE(0:NF_R_M,NRM),NGAP(NIM,NBM),
     '  NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),
     '  NKB(2,2,2,NNM,NBFM),NKEF(0:4,16,6,NBFM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),
     '  NLNO(NOPM,NXM),NMNO(1:2,0:NOPM,NXM),
     '  NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),
     '  NONL(NLM,NXM),NONM(NMM,NPM,NXM),
     '  NONY(0:NOYM,NYM,NRCM,0:NRM,NXM),NP_INTERFACE(0:NPM,0:3),
     '  NP1OPT(NOPM),NPF(9,NFM),NPL(5,0:3,NLM),NPLIST3(0:NPM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),NRE(NEM),NSB(NKM,NNM,NBFM),
     '  NRLIST(0:NRM),NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJL(4,NJM,NLM),NW(NEM,3,NXM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),NXLIST(0:NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),PAOPTY(NOPM),
     '  Z_CONT_LIST(NDM,2,7)
      REAL*8 AQ(NMAQM,NQM),CE(NMM,NEM,NXM),CELL_RCQS_VALUE(NQRM,NQVM),
     '  CG(NMM,NGM),CGE(NMM,NGM,NEM,NXM),
     '  CONY(0:NOYM,NYM,NRCM,0:NRM,NXM),CP(NMM,NPM,NXM),
     '  CURVCORRECT(2,2,NNM,NEM),DL(3,NLM),
     '  D_RE(NSM,NHM,NOPM),D_RI3(NHM*NSM),D_RP(NYM,NYM),
     '  D_TG(3,3,NHM*NSM),D_ZG(NHM,NUM,NHM*NSM),ES(NHM*NSM,NHM*NSM),
     '  FEXT(NIFEXTM,NGM,NEM),FGRAD(*),
     '  LAPL(NY_TRANSFER_M,NY_TRANSFER_M),
     '  LAPLSQR(NY_TRANSFER_M,NY_TRANSFER_M),
     '  MFI(NDM,NTSM,3,NSSM),
     '  PAOPTI(*),PBOPTI(*),PG(NSM,NUM,NGM,NBM),
     '  PHI(NY_TRANSFER_M,NTSM),PHI_H(NY_TRANSFER_M,NTSM),
     '  PMIN(*),PMAX(*),RE1(NSM,NHM),RE2(NSM,NHM),
     '  RESID(*),RESIDM(*),RESJAC(NREM,*),RG(NGM),SE(NSM,NBFM,NEM),
     '  T_BH(NY_TRANSFER_M,NY_TRANSFER_M),
     '  WG(NGM,NBM),WK1_INV(NY_TRANSFER_M,NY_TRANSFER_M),
     '  WU(0:NUM+1,NEM),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XG(NJM,NUM),XID(NIM,NDM),
     '  XIG(NIM,NGM,NBM),XN(NJM,NGM),XP(NKM,NVM,NJM,NPM),
     '  YG(NIYGM,NGM,NEM),YGF(NIYGFM,NGFM,NFM),YP(NYM,NIYM,NXM),
     '  ZA(NAM,NHM,NCM,NEM),ZA1(NAM,NHM,NCM,NEM),
     '  Z_CONT(NDM,2,67),ZD(NJM,NDM),ZE(NSM,NHM),
     '  ZG(NHM,NUM),ZG1(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM),ZP1(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER IBEG,IEND,nd,ndl,ne,noelem,nr,N3CO
      REAL*8 SOBVALUE,DATARESID,RFROMC,VALUE
      CHARACTER*11 TYPE
      LOGICAL ABBREV,CBBREV

      CALL ENTERS('UPSOBE',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM update Sobolev from [factor]
C###  Parameter:      <value #[1.0]>
C###  Specify the factor to update the Sobolev smoothing weights from
C###  Description:
C###    Updates the global Sobolev smoothing weight (WU(0,ne)) by
C###    muliplying the current global weight for each element by the
C###    the specified value.

        OP_STRING(1)=STRING(1:IEND)//' from [factor]'
     '    //' <value #[1.0]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update Sobolev from data
C###  Description:
C###    Updates the global Sobolev smoothing weight (WU(0,ne)) to
C###    be some ratio of the Sobolev objective and the data objective
C###    (Currently only implemted for data fitting by optimisation).

        OP_STRING(1)=STRING(1:IEND)//' from data'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update Sobolev from projections
C###  Description:
C###    Updates the global Sobolev smoothing weight (WU(0,ne)) such
C###    WU(0,ne) is zero for elements that contain data point
C###    projections and unchanged for elements that don't.

        OP_STRING(1)=STRING(1:IEND)//' from projections'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc','UPSOBE',ERROR,*9999)
      ELSE
        nr=1 !temporary

        CALL ASSERT(KTYP8.GT.0,'>>Define fit first',ERROR,*9999)
        CALL ASSERT(KTYP12.EQ.1.OR.KTYP12.EQ.2,
     '    '>>Sobolev smoothing not defined',ERROR,*9999)
        IF(CBBREV(CO,'FROM',2,noco+1,NTCO,N3CO)) THEN
          IF(ABBREV(CO(N3CO+1),'FACTOR',2)) THEN
            TYPE='FACTOR'
          ELSE IF(ABBREV(CO(N3CO+1),'DATA',2)) THEN
            TYPE='DATA'
          ELSE IF(ABBREV(CO(N3CO+1),'PROJECTIONS',2)) THEN
            CALL ASSERT(CALC_XI,'>>Calculate xi positions first',
     '        ERROR,*9999)
            TYPE='PROJECTIONS'
          ELSE
            TYPE='FACTOR'
          ENDIF
        ELSE
          TYPE='FACTOR'
        ENDIF
        IF(TYPE(1:6).EQ.'FACTOR') THEN
          IF(CBBREV(CO,'VALUE',2,noco+1,NTCO,N3CO)) THEN
            IF(N3CO.LT.NTCO) THEN
              VALUE=RFROMC(CO(N3CO+1))
            ELSE
              VALUE=1.0d0
            ENDIF
          ELSE
            VALUE=1.0d0
          ENDIF
        ENDIF

        IF(TYPE(1:4).EQ.'DATA') THEN
          IF(KTYP8.EQ.6) THEN ! Data fitting by optimisation
            IF(KTYP29.EQ.1) THEN
              CO(1)='FEM'
              CO(2)='EVALUATE'
              CO(3)='RESIDUALS'
              CO(4)='WRT'
              CO(5)='DATA_FITTING'
              NTCO=5
              noco=3
              CALL EVRESI(IBT,IDO,INP,ISIZE_MFI,
     '          ISIZE_PHI,ISIZE_TBH,LD,LDR,LGE,
     '          NAN,NBH,NBHF,NBJ,NBJF,NEELEM,NEL,NELIST,%VAL(0),
     '          NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,NKHE,
     '          NKJE,NLL,NLNO,NMNO,NNB,NNF,NNL,NONL,NONM,NONY,
     '          NP1OPT,NPF,NP_INTERFACE,NPL,NPLIST3,NPNE,
     '          NPNODE,NPNY,NRE,NRLIST,NSB,NVHE,NVHP,NVJE,NVJL,
     '          NW,NXI,NXLIST,NYNE,NYNO,NYNP,NYNR,PAOPTY,Z_CONT_LIST,
     '          AQ,CE,CELL_RCQS_VALUE,CG,CGE,CONY,CP,
     '          CURVCORRECT,DL,
     '          D_RE,D_RI3,D_RP,D_TG,D_ZG,ES,FEXT,FGRAD,LAPL,LAPLSQR,
     '          MFI,PAOPTI,PBOPTI,PG,PHI,PHI_H,PMIN,PMAX,RE1,RE2,RESID,
     '          RESJAC,RG,SE,T_BH,
     '          WG,WK1_INV,WU,XA,XE,XG,XID,XIG,XN,XP,YG,YGF,YP,ZA,ZA1,
     '          Z_CONT,
     '          ZD,ZE,ZG,ZG1,ZP,ZP1,STRING,FIX,ERROR,*9999)
C             DO noelem=1,NEELEM(0,nr)
C               ne=NEELEM(noelem,nr)
C               DATARESID=0.d0
C               DO ndl=1,NDLT(ne)
C                 nd=NDDL(ne,ndl)
C                 DATARESID=DATARESID+RESID(nd)**2
C             ENDDO
C             SOBVALUE=WU(10,ne)**2 !The Sobolev value for that elem
C             IF(SOBVALUE.EQ.0.d0) THEN
C               WU(0,ne)=1.d0
C             ELSE
C               WU(0,ne)=DATARESID/(DBLE(NDLT(ne))*SOBVALUE)
C             ENDIF
C             IF(DOP) THEN
C               WRITE(OP_STRING,'('' Data residual (ne='',I4,'') = '',
C     '           E12.5,''WU(10,ne) = '',E12.5,'', WU(0,ne) = '',
C     '           E12.5)') ne,DATARESID,WU(10,ne),WU(0,ne)
C               CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C             ENDIF
C             ENDDO
              DATARESID=0.d0
              DO nd=1,NDT
                DATARESID=DATARESID+RESID(nd)**2
              ENDDO
              SOBVALUE=RESID(NT_RES)**2
              IF(SOBVALUE.EQ.0.d0) THEN
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  WU(0,ne)=1.d0
                ENDDO
              ELSE
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  WU(0,ne)=DATARESID/(DBLE(NDT)*SOBVALUE)
                ENDDO
              ENDIF
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'('' Data residual = '','
     '            //'E12.5,'' Sobolev residual = '',E12.5,'
     '            //''' WU(0,1) = '',E12.5)')
     '            DATARESID,SOBVALUE,WU(0,NEELEM(1,nr))
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF

            ELSE IF(KTYP29.EQ.2) THEN
              CO(1)='FEM'
              CO(2)='EVALUATE'
              CO(3)='OBJECTIVE'
              CO(4)='WITH'
              CO(5)='OPTIMISATION'
              NTCO=5
              noco=3
              CALL EVOBJE(IBT,IDO,INP,LD,LGE,NAN,NBH,NBHF,NBJ,NBJF,
     '          NEELEM,NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,
     '          NKHE,NKJE,NLL,NLNO,NMNO,NNB,NNF,NNL,NONL,NONM,NONY,
     '          NP_INTERFACE,NP1OPT,NPF,NPL,NPNE,NPNODE,NPNY,NRE,NSB,
     '          NVHE,NVHP,NVJE,NW,NXI,NYNE,NYNO,NYNP,NYNR,PAOPTY,
     '          Z_CONT_LIST,
     '          CE,CG,CGE,CP,CURVCORRECT,DL,FEXT,
     '          FGRAD,PAOPTI,PBOPTI,PG,PMIN,PMAX,RE1,RESIDM,
     '          RESJAC,RG,SE,WG,WU,XA,XE,XG,XID,XIG,XP,YG,
     '          YGF,YP,ZA,
     '          ZA1,Z_CONT,ZD,ZE,ZG,ZP,ZP1,STRING,FIX,
     '          ERROR,*9999)
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                DATARESID=0.0d0
                DO ndl=1,NDLT(ne)
                  nd=NDDL(ne,ndl)
                  DATARESID=DATARESID+RESIDM(nd)
                ENDDO
                SOBVALUE=WU(NUM+1,ne) !The Sobolev value for that elem
                IF(SOBVALUE.EQ.0.d0) THEN
                  WU(0,ne)=1.d0
                ELSE
                  WU(0,ne)=DATARESID/(DBLE(NDLT(ne))*SOBVALUE)
                ENDIF
              ENDDO
            ENDIF
          ENDIF
        ELSE IF(TYPE(1:6).EQ.'FACTOR') THEN
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            WU(0,ne)=Wu(0,ne)*VALUE
          ENDDO !noelem (ne)
        ELSE IF(TYPE(1:11).EQ.'PROJECTIONS') THEN
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            IF(NDLT(ne).NE.0) THEN
              WU(0,ne)=0.0d0
            ENDIF
          ENDDO !noelem (ne)
        ENDIF
      ENDIF

      CALL EXITS('UPSOBE')
      RETURN
 9999 CALL ERRORS('UPSOBE',ERROR)
      CALL EXITS('UPSOBE')
      RETURN 1
      END
CC AJPE

