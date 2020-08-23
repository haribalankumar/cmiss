      SUBROUTINE UPRESI(IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,NEELEM,
     '  NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,NKHE,NKJE,NNB,NNF,NPF,
     '  NPNE,NPNODE,NPNY,NRE,NSB,NVHE,NVHP,NVJE,NW,NXI,NXLIST,NYNE,
     '  NYNP,NYNR,Z_CONT_LIST,CE,CG,CGE,CP,CURVCORRECT,FEXT,
     '  PG,RE,RG,
     '  SE,WG,XA,XG,XP,YG,YGF,YP,ZA,ZA1,Z_CONT,ZE,ZD,
     '  ZP,ZP1,STRING,
     '  FIX,ERROR,*)

C#### Subroutine: UPRESI
C###  Description:
C###    UPRESI updates finite element node residual vector YP(ny,4) for
C###    nonlinear problems by calling ZPRP with current material
C###    parameters and initial conditions. If requested the increments
C###    are added in to ZP.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LGE(NHM*NSM),NAN(NIM,NAM,NBFM),
     '  NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NEELEM(0:NE_R_M,0:NRM),NFF(6,NEM),NFFACE(0:NF_R_M,NRM),
     '  NGAP(NIM,NBM),NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),
     '  NKB(2,2,2,NNM,NBFM),NKEF(0:4,16,6,NBFM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),NRE(NEM),NSB(NKM,NNM,NBFM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM),NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NXLIST(0:NXM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),Z_CONT_LIST(NDM,2,7)
      REAL*8 CE(NMM,NEM,NXM),CG(NMM,NGM),CGE(NMM,NGM,NEM,NXM),
     '  CP(NMM,NPM,NXM),CURVCORRECT(2,2,NNM,NEM),FEXT(NIFEXTM,NGM,NEM),
     '  PG(NSM,NUM,NGM,NBM),RE(NSM,NHM),
     '  RG(NGM),SE(NSM,NBFM,NEM),WG(NGM,NBM),
     '  XA(NAM,NJM),XG(NJM,NUM),XP(NKM,NVM,NJM,NPM),
     '  YG(NIYGM,NGM,NEM),YGF(NIYGFM,NGFM,NFM),YP(NYM,NIYM,NXM),
     '  ZA1(NAM,NHM,NCM,NEM),ZA(NAM,NHM,NCM,NEM),
     '  Z_CONT(NDM,2,67),ZE(NSM,NHM),ZD(NJM,NDM),
     '  ZP(NKM,NVM,NHM,NPM,NCM),ZP1(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,INDEX,N3CO,nc,nh,nhx,nk,nonode,
     '  np,nr,nv,nx,nxc,ny,POINT,RESIDUAL
      REAL*8 RFROMC,WEIGHT
      LOGICAL CBBREV,DATA,INCREM

      CALL ENTERS('UPRESI',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM update residuals
C###  Parameter:      <increment>
C###    Specify wiether increments are added in to ZP
C###  Parameter:      <region #[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description: updates finite element node residual vector YP(ny,4) for
C###    nonlinear problems by calling ZPRP with current material
C###    parameters and initial conditions. If requested the increments
C###    are added in to ZP.
C###

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<increment>'
        OP_STRING(3)=BLANK(1:15)//'<region #[1]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update residuals data
C###  Parameter:      <point #[1]>
C###    Specify the data point number
C###  Parameter:      <index #[1]>
C###    Specify the data point field
C###  Parameter:      <weight #[1.0]>
C###    Specify the weight to scale the data point value
C###  Parameter:      <residual #[1]>
C###    Specify the index of the residual vector to load
C###  Description: updates vector YP(residual,9) with data value

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<data>'
        OP_STRING(3)=BLANK(1:15)//'<point #[1]>'
        OP_STRING(4)=BLANK(1:15)//'<index #[1]>'
        OP_STRING(5)=BLANK(1:15)//'<weight #[1.0]>'
        OP_STRING(6)=BLANK(1:15)//'<residual #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc','UPRESI',ERROR,*9999)
      ELSE

        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)

        IF(CBBREV(CO,'REGION',1,noco+1,NTCO,N3CO)) THEN
          nr=IFROMC(CO(N3CO+1))
        ELSE
          nr=1
        ENDIF

        IF(CBBREV(CO,'INCREMENT',1,noco+1,NTCO,N3CO)) THEN
          INCREM=.TRUE.
        ELSE
          INCREM=.FALSE.
        ENDIF

        IF(CBBREV(CO,'DATA',1,noco+1,NTCO,N3CO)) THEN
          DATA=.TRUE.
        ELSE
          DATA=.FALSE.
        ENDIF

        IF(DATA) THEN
          IF(CBBREV(CO,'POINT',3,noco+1,NTCO,N3CO)) THEN
            POINT=IFROMC(CO(N3CO+1))
          ELSE
            POINT=1
          ENDIF
          IF(CBBREV(CO,'INDEX',3,noco+1,NTCO,N3CO)) THEN
            INDEX=IFROMC(CO(N3CO+1))
          ELSE
            INDEX=1
          ENDIF
C new HS/MPN 15/7/05 adding weight scale factor to data point value
          IF(CBBREV(CO,'WEIGHT',3,noco+1,NTCO,N3CO)) THEN
            WEIGHT=RFROMC(CO(N3CO+1))
          ELSE
            WEIGHT=1.0D0
          ENDIF
          IF(CBBREV(CO,'RESIDUAL',3,noco+1,NTCO,N3CO)) THEN
            RESIDUAL=IFROMC(CO(N3CO+1))
          ELSE
            RESIDUAL=1
          ENDIF
C news HS 17/5/04 changing index of YP from 8 to 9 since it interfered with the line search option in LSFUNC.f
C see also RESFUN.f line 1406, had to be adjusted for residual calculation.
          YP(RESIDUAL,9,1)=WEIGHT*ZD(INDEX,POINT)
C newe

        ELSE
C CPB 8/6/94 Adding NX_LOC
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '      ERROR,*9999)
          CALL ASSERT(ITYP6(nr,nx).EQ.2,
     '      '>>Not nonlinear problem',ERROR,*9999)

          CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),NKH(1,1,1,nr),
     '      NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,
     '      YP(1,1,nx),ZA,ZP,ERROR,*9999)
          IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN !cnst Vol
C           Put reference state for cavity from YP(ny,10) into
C           ZA1,ZP1 for ZERE55
            CALL YPZP(10,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '        NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,
     '        YP(1,1,nx),ZA1,ZP1,ERROR,*9999)
          ENDIF
          IF(INCREM) THEN !add incremental b.c.s
            nc=1 !TEMPORARY AJP 18-12-91
            ny=0
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              DO nhx=1,NHP(np,nr,nx)
                nh=NH_LOC(nhx,nx)
                DO nv=1,NVHP(nh,np,nc,nr)
                  DO nk=1,NKH(nh,np,nc,nr)
                    ny=ny+1
                    IF(FIX(ny,1,nx)) THEN
                      ZP(nk,nv,nh,np,nc)=ZP(nk,nv,nh,np,nc)+YP(ny,2,nx)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF

          CALL ZPRP(IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,NEELEM,
     '      NFF,NFFACE,NGAP,NHE(1,nx),NHP(1,nr,nx),NKB,NKEF,
     '      NKH(1,1,1,nr),NKHE,NKJE,NNB,NNF,
     '      NPF,NPNE,NPNODE,NPNY(0,1,0,
     '      nx),nr,NRE,NSB,NVHE,NVHP(1,1,1,nr),
     '      NVJE,NW(1,1,nx),nx,NXI,
     '      NYNE,NYNP,NYNR(0,0,1,nr,nx),Z_CONT_LIST,CE(1,1,nx),CG,
     '      CGE(1,1,1,nx),CP(1,1,nx),CURVCORRECT,
     '      FEXT,PG,RE,RG,SE,WG,XA,
     '      XG,XP,YG,YGF,YP(1,1,nx),ZA,ZA1,%VAL(0),Z_CONT,
     '      ZE,%VAL(0),ZP,ZP1,%VAL(0),FIX,
     '      ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('UPRESI')
      RETURN
 9999 CALL ERRORS('UPRESI',ERROR)
      CALL EXITS('UPRESI')
      RETURN 1
      END


