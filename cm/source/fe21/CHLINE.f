      SUBROUTINE CHLINE(IBT,IDO,INP,ISEG,ISELNO,ISFIBR,ISFIEL,ISLINE,
     '  ISLINO,ISL2BE,ISL3BE,ISNONO,ISN2BE,ISN3BE,MXI,NAN,NBH,NBJ,
     '  NEELEM,NEL,NGAP,NHE,NHP,NKH,NKHE,NKJ,NKJE,NLATNE,NLL,
     '  NLLIST,NPF,NPL,NPNE,NPNODE,NQNE,NQNLAT,NQS,NQXI,NRE,
     '  NVHE,NVHP,NVJE,NVJL,NW,NYNE,NYNP,CURVCORRECT,
     '  DL,SE,XA,XB,XE,XG,XP,XQ,YG,YP,YQ,ZA,ZE,ZP,
     '  CSEG,STRING,FIX,ERROR,*)

C#### Subroutine: CHLINE
C###  Description:
C###    CHLINE changes lines.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'fbgr00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISEG(*),ISELNO(NWM,NEM),ISFIBR(NWM,NEM,NGRSEGM),ISFIEL(NWM,NEM),
     '  ISLINE(NWM,2*NGRSEGM),ISLINO(NWM),ISL2BE(NLM),ISL3BE(NLM),
     '  ISNONO(NWM,NPM),ISN2BE(NLM),ISN3BE(NLM),MXI(2,NEM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NEL(0:NELM,NLM),NGAP(NIM,NBM),NHE(NEM,NXM),
     '  NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),NLATNE(NEQM+1),
     '  NLL(12,NEM),NLLIST(0:NLM),NPF(9,NFM),NPL(5,0:3,NLM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NQNE(NEQM,NQEM),NQNLAT(NEQM*NQEM),
     '  NQS(NEQM),NQXI(0:NIM,NQSCM),
     '  NRE(NEM),NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM), 
     '  NVJE(NNM,NBFM,NJM,NEM),NVJL(4,NJM,NLM),NW(NEM,3,NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),DL(3,NLM),SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM),XB(2,NJM,NLM),
     '  XE(NSM,NJM),XG(NJM,NUM),XP(NKM,NVM,NJM,NPM),XQ(NJM,NQM),
     '  YG(NIYGM,NGM,NEM),YP(NYM,NIYM,NXM),YQ(NYQM,NIQM,NAM,NXM),
     '  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,iw,IWK(6),N3CO,noiw,nolist,NTIW,nx
      LOGICAL ABBREV,CBBREV,MOUSE

      CALL ENTERS('CHLINE',*9999)

C ??? cpb 22/11/94 why is nx used here ?
      nx=1 !temporary cpb 22/11/94

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM change lines;m
C###  Parameter:      <on WS#[1]>
C###    Specify the workstation (GX window) to draw the
C###    points on.
C###  Description:
C###    Change the shape of the specified lines using Bezier control
C###    points.  Use the left mouse button to 'pick' a control point
C###    and then to relocate it.  Use the central mouse button to quit
C###    when you have finished changing lines.

        OP_STRING(1)=STRING(1:IEND)//';m <on WS#[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM change lines
C###  Parameter:      <with (LINES/all)[all]>
C###  Parameter:      <on WS#[1]>
C###  Description:
C###

        OP_STRING(1)=STRING(1:IEND)//' <with (LINES/all)[all]>'
     '    //' <on WS#[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','CHLINE',ERROR,*9999)
      ELSE
        CALL CHECKQ('M',noco,1,CO,COQU,STRING,*1)
        MOUSE=.FALSE.
        IF(ABBREV(COQU(noco,1),'P',1)) THEN
          IOTYPE=1
        ELSE IF(ABBREV(COQU(noco,1),'M',1)) THEN
          MOUSE=.TRUE.
        ENDIF
        IF(CBBREV(CO,'WITH',1,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NL_R_M,NLLIST(0),NLLIST(1),ERROR,*9999)
        ELSE
          NLLIST(0)=NLT
          DO nolist=1,NLLIST(0)
            NLLIST(nolist)=nolist
          ENDDO
        ENDIF
        IF(CBBREV(CO,'ON',1,noco+1,NTCO,N3CO)) THEN
          IWK(1)=IFROMC(CO(N3CO+1))
        ELSE
          IWK(1)=1
        ENDIF

        IF(MOUSE) THEN
          iw=IWK(1)
          IWK(2)=iw
C          AUX=.FALSE.
          IF(iw.EQ.5.OR.iw.EQ.6) THEN
            RHTRAN=iw.EQ.5
            LFTRAN=iw.EQ.6
            IF(RHTRAN.OR.LFTRAN)THEN
c              CALL WKST_WINDOW(3,0.0,511.0,0.0,511.0,ERROR,*9999)
C              DIAG=DSQRT(511.d0*512.d0*2.d0)
            ENDIF
            NTIW=2
            IF(NJT.EQ.2)THEN
              IWK(2)=1
            ELSE
              IWK(2)=3
            ENDIF
            CALL UPVIEW(IBT,IDO,INP,ISEG,ISELNO,ISFIBR,ISFIEL,ISLINE,
     '        ISLINO,ISNONO,IWK,MXI,NAN,NBH,NBJ,NEELEM,NGAP,NHE(1,nx),
     '        NHP(1,0,nx),NKH,NKHE,NKJE,NLATNE,NLL,NLLIST,NPF,NPL,NPNE,
     '        NPNODE,NQNE,NQNLAT,
     '        NQS,NQXI,NRE,NTIW,NVHE,NVHP,NVJE,NW(1,1,nx),
     '        nx,NYNE,NYNP,CURVCORRECT,DL,SE,XA,XE,XG,XP,XQ,
     '        YG,YP(1,1,nx),
     '        YQ,ZA,ZE,ZP,CSEG,FIX(1,1,nx),ERROR,*9999)
          ENDIF
          CALL CRHERM(IDO,ISEG,ISLINE,ISLINO,
     '      ISL2BE,ISL3BE,ISN2BE,ISN3BE,IWK(2),NBJ,NEL,NKJ,
     '      NLLIST,NPL,NVJL,DL,XB,XP,ZP,CSEG,ERROR,*9999)
          IF(NJT.EQ.2) THEN
            NTIW=1
            IWK(1)=1
          ELSE IF(NJT.GT.2) THEN
            IF(iw.EQ.1) THEN
              iw=2
            ELSE IF(iw.EQ.2) THEN
              iw=1
            ELSE IF(iw.EQ.5) THEN
              iw=6
            ELSE IF(iw.EQ.6) THEN
              iw=5
            ENDIF
C            AUX=.TRUE.
            IWK(1)=iw
            IWK(2)=iw
            NTIW=1
            IF(iw.EQ.5.OR.iw.EQ.6) THEN
              RHTRAN=iw.EQ.5
              LFTRAN=iw.EQ.6
              NTIW=2
              IWK(2)=3
            ENDIF
            CALL UPVIEW(IBT,IDO,INP,ISEG,ISELNO,ISFIBR,ISFIEL,ISLINE,
     '        ISLINO,ISNONO,IWK,MXI,NAN,NBH,NBJ,NEELEM,NGAP,NHE(1,nx),
     '        NHP(1,0,nx),NKH,NKHE,NKJE,NLATNE,NLL,NLLIST,NPF,NPL,
     '        NPNE,NPNODE,NQNE,NQNLAT,NQS,NQXI,NRE,NTIW,NVHE,NVHP,NVJE,
     '        NW(1,1,nx),nx,NYNE,NYNP,CURVCORRECT,
     '        DL,SE,XA,XE,XG,XP,XQ,YG,
     '        YP(1,1,nx),YQ,ZA,ZE,ZP,CSEG,FIX(1,1,nx),ERROR,*9999)
            CALL CRHERM(IDO,ISEG,ISLINE,ISLINO,
     '        ISL2BE,ISL3BE,ISN2BE,ISN3BE,IWK(2),NBJ,NEL,NKJ,
     '        NLLIST,NPL,NVJL,DL,XB,XP,ZP,CSEG,ERROR,*9999)
            NTIW=2*NJT-3
            DO noiw=1,NTIW
              IWK(noiw)=noiw
            ENDDO
            IF(FBGRAF) THEN
              NTIW=5
              IWK(4)=5
              IWK(5)=6
            ENDIF
          ENDIF
          CALL UPVIEW(IBT,IDO,INP,ISEG,ISELNO,ISFIBR,ISFIEL,ISLINE,
     '      ISLINO,ISNONO,IWK,MXI,NAN,NBH,NBJ,NEELEM,NGAP,NHE(1,nx),
     '      NHP(1,0,nx),NKH,NKHE,NKJE,NLATNE,NLL,NLLIST,NPF,NPL,NPNE,
     '      NPNODE,NQNE,NQNLAT,
     '      NQS,NQXI,NRE,NTIW,NVHE,NVHP,NVJE,NW(1,1,nx),
     '      nx,NYNE,NYNP,CURVCORRECT,DL,SE,XA,XE,XG,XP,XQ,
     '      YG,YP(1,1,nx),YQ,ZA,ZE,
     '      ZP,CSEG,FIX(1,1,nx),ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('CHLINE')
      RETURN
 9999 CALL ERRORS('CHLINE',ERROR)
      CALL EXITS('CHLINE')
      RETURN 1
      END


