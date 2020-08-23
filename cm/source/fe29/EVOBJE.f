      SUBROUTINE EVOBJE(IBT,IDO,INP,LD,LGE,NAN,NBH,NBHF,NBJ,NBJF,
     '  NEELEM,NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,NKHE,NKJE,
     '  NLL,NLNO,NMNO,NNB,NNF,NNL,NONL,NONM,NONY,NP_INTERFACE,
     '  NP1OPT,NPF,NPL,NPNE,NPNODE,NPNY,NRE,NSB,
     '  NVHE,NVHP,NVJE,NW,NXI,NYNE,NYNO,NYNP,NYNR,PAOPTY,Z_CONT_LIST,
     '  CE,CG,CGE,CP,CURVCORRECT,DL,FEXT,FGRAD,PAOPTI,
     '  PBOPTI,PG,PMIN,
     '  PMAX,RE1,RESIDM,RESJAC,RG,SE,WG,WU,XA,XE,
     '  XG,XID,XIG,XP,
     '  YG,YGF,YP,ZA,ZA1,Z_CONT,ZD,ZE,ZG,ZP,ZP1,STRING,FIX,
     '  ERROR,*)

C#### Subroutine: EVOBJE
C###  Description:
C###    EVOBJE evaluates the objective function.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'opti00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LD(NDM),LGE(NHM*NSM),NAN(NIM,NAM,NBFM),
     '  NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NEELEM(0:NE_R_M,0:NRM),NFF(6,NEM),NFFACE(0:NF_R_M,NRM),
     '  NGAP(NIM,NBM),NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),
     '  NKB(2,2,2,NNM,NBFM),NKEF(0:4,16,6,NBFM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NKH(NHM,NPM,NCM,0:NRM),NLL(12,NEM),
     '  NLNO(NOPM,NXM),NMNO(1:2,0:NOPM,NXM),
     '  NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),NONL(NLM,NXM),
     '  NONM(NMM,NPM,NXM),NONY(0:NOYM,NYM,NRCM,0:NRM,NXM),
     '  NP_INTERFACE(0:NPM,0:3),NP1OPT(NOPM),NPF(9,NFM),NPL(5,0:3,NLM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),NRE(NEM),NSB(NKM,NNM,NBFM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM),NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),NYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),PAOPTY(NOPM),
     '  Z_CONT_LIST(NDM,2,7)
      REAL*8 CE(NMM,NEM,NXM),CG(NMM,NGM),CGE(NMM,NGM,NEM,NXM),
     '  CP(NMM,NPM,NXM),CURVCORRECT(2,2,NNM,NEM),DL(3,NLM),
     '  FEXT(NIFEXTM,NGM,NEM),FGRAD(NOPM),
     '  PAOPTI(NOPM),PBOPTI(NOPM),PG(NSM,NUM,NGM,NBM),
     '  PMIN(NOPM+NCOM),PMAX(NOPM+NCOM),RE1(NSM,NHM),RESIDM(NREM),
     '  RESJAC(NREM,NOPM),RG(NGM),SE(NSM,NBFM,NEM),
     '  WG(NGM,NBM),WU(0:NUM+1,NEM),XA(NAM,NJM),XE(NSM,NJM),XG(NJM,NUM),
     '  XIG(NIM,NGM,NBM),XID(NIM,NDM),XP(NKM,NVM,NJM,NPM),
     '  YG(NIYGM,NGM,NEM),YGF(NIYGFM,NGFM,NFM),YP(NYM,NIYM,NXM),
     '  ZA(NAM,NHM,NCM,NEM),ZA1(NAM,NHM,NCM,NEM),Z_CONT(NDM,2,67),
     '  ZD(NJM,NDM),ZE(NSM,NHM),ZG(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM),ZP1(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,IFROMC,
     '  n1step,n2step,n3step,
     '  N3CO,ne,no,noelem,noopti,nopara,nostep,nm,nr,NTLIST,
     '  NTSTEP,NTSTEP1,NTSTEP2,NTSTEP3,nxc,nx
      REAL*8 DELTA,DFUNC,FUNC1,FUNC2,P1STEP,P2STEP,P3STEP,
     '  TEMP(12),ZE1(NSM,NHM)
      CHARACTER CHAR1*144,CHAR2*11,CTEMP*30,FILE*(MXCH)
      LOGICAL CBBREV,DERIV,OPFILE

      CALL ENTERS('EVOBJE',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
! Note the param values require a nx for default values.
        nx=1 !temporary
        CALL STRING_TRIM(STRING,IBEG,IEND)
        WRITE(CHAR1,'(E11.4)') PAOPTI(1)
        DO noopti=2,NMNO(1,0,nx)
          CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
          WRITE(CHAR2,'(E11.4)') PAOPTI(noopti)
          CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
          CHAR1=CHAR1(IBEG1:IEND1)//','//CHAR2(IBEG2:IEND2)
        ENDDO
        CALL STRING_TRIM(CHAR1,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM evaluate objective<;FILENAME> with
C###  Parameter:         optimisation_parameters
C###
C###  Parameter:        <step NO_STEPS#[0]>
C###    Specify the number of steps to use.
C###  Parameter:        <region #[1]>
C###    Specify the region numbers to evaluate.
C###  Description:
C###    Evaluates an objective function.
C###      Used with optimisation routines.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
     '    //' with optimisation_parameters'
        OP_STRING(2)=BLANK(1:15)//'<step NO_STEPS#[0]>'
        OP_STRING(3)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM evaluate objective with material_parameters
C###  Parameter:        <region #[1]>
C###    Specify the region numbers to evaluate.
C###  Description:
C###    Evaluates an objective with respect to material
C###      parameters. Used with optimisation routines.

        OP_STRING(1)=STRING(1:IEND)
     '    //' with material_parameters'
        OP_STRING(2)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM evaluate objective with <PARAM_VALUES[0.0000E+00]>
C###  Parameter:        <region #[1]>
C###    Specify the region numbers to evaluate.
C###  Description:
C###

        OP_STRING(1)=STRING(1:IEND)
     '    //' with <PARAM_VALUES['//CHAR1(IBEG1:IEND1)//']>'
        OP_STRING(2)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM evaluate objective derivatives
C###  Parameter:        <wrt PARAMETER#[wrt 1]>
C###    Specify the parameters with evaluate wrt.
C###  Parameter:        <region #[1]>
C###    Specify the region numbers to evaluate.
C###  Description:
C###    Evaluates the derivatives of an objective.
C###      Used with optimisation routines.

        OP_STRING(1)=STRING(1:IEND)//' derivatives'
        OP_STRING(2)=BLANK(1:15)//'<wrt PARAMETER#[wrt 1]>'
        OP_STRING(3)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','EVOBJE',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN
          FILE=COQU(noco,1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          OPFILE=.TRUE.
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opobje','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
          IOFI=IOOP
        ENDIF

        CALL ASSERT(USE_NPSOL.GT.0,'>>Cannot evaluate objective'
     '    //' as USE_NPSOL is zero in parameters.inc',ERROR,*9999)

        nxc=1 !temporary
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_OPTI,ERROR,*9999)
        CALL ASSERT(nx.GT.0,
     '    '>>no nx defined for this optimisation class',ERROR,*9999)

        IF(CBBREV(CO,'REGION',1,noco+1,NTCO,N3CO)) THEN
          nr=IFROMC(CO(N3CO+1))
        ELSE
          nr=1
        ENDIF

        IF(CBBREV(CO,'WITH',2,noco+1,NTCO,N3CO)) THEN
          DERIV=.FALSE.
          IF(CBBREV(CO,'OPTIMISATION_PARAMETERS',1,N3CO+1,N3CO+1,N3CO))
     '      THEN
C            CALL ASSERT(NMNO(1,0,nx).LE.12,
C     '        '>>Increase dimension of PALIST',ERROR,*9999)
            DO noopti=1,NMNO(1,0,nx)
              PBOPTI(noopti)=PAOPTI(noopti)
            ENDDO
            IF(CBBREV(CO,'STEP',1,noco+1,NTCO,N3CO)) THEN
              NTSTEP=IFROMC(CO(N3CO+1))
            ELSE
              NTSTEP=0
            ENDIF
          ELSE IF(CBBREV(CO,'MATERIAL_PARAMETERS',1,N3CO+1,N3CO+1,
     '      N3CO)) THEN
            DO noopti=1,NMNO(1,0,nx)
              nm=NMNO(1,noopti,nx)
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                no=NONM(nm,noelem,nx)
                PBOPTI(no)=CE(nm,ne,nx)
              ENDDO
            ENDDO
          ELSE !use specified list of parameters
            CALL PARSRL(CO(N3CO+1),NOPM,NTLIST,PBOPTI,ERROR,*9999)
            NTSTEP=0
          ENDIF

        ELSE IF(CBBREV(CO,'DERIVATIVES',1,noco+1,NTCO,N3CO)) THEN
          DERIV=.TRUE.
          IF(CBBREV(CO,'WRT',2,noco+2,NTCO,N3CO)) THEN
            nopara=IFROMC(CO(N3CO+1))
          ELSE
            nopara=1
          ENDIF
          NTSTEP=0

        ELSE
          CALL STRING_TRIM(STRING,IBEG,IEND)
          CALL STRING_TRIM(CO(noco),IBEG1,IEND1)
          CTEMP=CO(noco)(IBEG1:IEND1)
          CALL STRING_TRIM(CTEMP,IBEG1,IEND1)
          STRING=STRING(1:IEND)//' '//CTEMP(IBEG1:IEND1)
          CO(noco+1)='?'
          GO TO 1
        ENDIF

        IF(KTYP26.EQ.1) THEN !Store material parameters temporarily
          CALL ASSERT(NMNO(1,0,nx).LE.12,'>>Increase dimension of TEMP',
     '      ERROR,*9999)
          DO noopti=1,NMNO(1,0,nx)
            nm=NMNO(1,noopti,nx)
            TEMP(noopti)=CE(nm,NEELEM(1,nr),nx)
            !NOTE: needs updating for varying case
          ENDDO
        ENDIF

        IF(NTSTEP.EQ.0) THEN !evaluate objective at current parameters
          IF(KTYP26.EQ.1) THEN
            DO noopti=1,NMNO(1,0,nx)
              nm=NMNO(1,noopti,nx)
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                CE(nm,ne,nx)=PBOPTI(noopti)
              ENDDO
            ENDDO
          ENDIF
          CALL OBJFUN(IBT,IDO,INP,LD,LGE,NAN,NBH,NBHF,NBJ,NBJF,
     '      NEELEM,NFF,NFFACE,NGAP,NHE(1,nx),NHP(1,nr,nx),
     '      NKB,NKEF,NKH,NKHE,NKJE,NLL,NLNO(1,nx),NNB,NNF,NNL,
     '      NONL(1,nx),NONY(0,1,1,nr,nx),NP_INTERFACE,NP1OPT,
     '      NPF,NPL,NPNE,NPNODE,NPNY(0,1,0,nx),nr,NRE,NSB,
     '      NVHE,NVHP(1,1,1,nr),NVJE,NW(1,1,nx),
     '      nx,NXI,NYNE,NYNO(0,1,1,nr,nx),NYNP,NYNR(0,0,1,nr,nx),
     '      PAOPTY,Z_CONT_LIST,CE(1,1,nx),CG,CGE(1,1,1,nx),
     '      CURVCORRECT,CP(1,1,nx),DL,FEXT,FGRAD,PAOPTI,PG,RE1,
     '      RESIDM,RESJAC,RG,SE,WG,wu,XA,XE,XG,XID,XIG,
     '      XP,YG,YGF,
     '      YP(1,1,nx),ZA,ZA1,Z_CONT,ZD,ZE,ZE1,ZG,ZP,ZP1,
     '      FIX(1,1,nx),ERROR,*9999)

          IF(DERIV) THEN !calc deriv of object func wrt param
            FUNC1=FUNC
            DELTA=1.0D-4
            nm=NMNO(1,nopara,nx)
            IF(KTYP26.EQ.1) THEN
              DO noelem=1,NEELEM(0,nr) !perturb param by delta
                ne=NEELEM(noelem,nr)
                CE(nm,ne,nx)=PBOPTI(nopara)+DELTA
              ENDDO
            ENDIF
            CALL OBJFUN(IBT,IDO,INP,LD,LGE,NAN,NBH,NBHF,NBJ,NBJF,
     '        NEELEM,NFF,NFFACE,NGAP,NHE(1,nx),NHP(1,nr,nx),
     '        NKB,NKEF,NKH,NKHE,NKJE,NLL,NLNO(1,nx),NNB,NNF,NNL,
     '        NONL(1,nx),NONY(0,1,1,nr,nx),NP_INTERFACE,NP1OPT,
     '        NPF,NPL,NPNE,NPNODE,NPNY(0,1,0,nx),nr,NRE,NSB,
     '        NVHE,NVHP(1,1,1,nr),NVJE,NW(1,1,nx),
     '        nx,NXI,NYNE,NYNO(0,1,1,nr,nx),NYNP,NYNR(0,0,1,nr,nx),
     '        PAOPTY,Z_CONT_LIST,CE(1,1,nx),CG,CGE(1,1,1,nx),
     '        CP(1,1,nx),CURVCORRECT,DL,FEXT,FGRAD,PAOPTI,PG,RE1,
     '        RESIDM,RESJAC,RG,SE,WG,WU,XA,XE,XG,XID,XIG,
     '        XP,YG,YGF,
     '        YP(1,1,nx),ZA,ZA1,Z_CONT,ZD,ZE,ZE1,ZG,ZP,ZP1,
     '        FIX(1,1,nx),ERROR,*9999)
            FUNC2=FUNC
            DFUNC=(FUNC2-FUNC1)/DELTA
            WRITE(OP_STRING,'('' Obj. func deriv wrt opt param #'','
     '        //'I1,'' = '',E11.3,'' at parameters:'',8E11.3)')
     '        nopara,DFUNC,(PBOPTI(noopti),noopti=1,NTOPTI)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          ELSE IF((ktyp26.ne.1.or.ktyp27.ne.2).AND.  !not sum of sqrd reac diffs
     '            (ktyp26.ne.2.or.ktyp27.ne.5).AND.  !not data fitting
     '            (ktyp26.ne.2.or.ktyp27.ne.6).AND.  !not fluid interface cond.
     '            (ktyp27.ne.3).AND.                 !not zero flux differences
     '            (ktyp27.ne.4))THEN                 !not hydrostatic Pr. cond.
C         ...all but KTYP27=3 use routine E04UPF (has its own print monitor)
            WRITE(OP_STRING,'('' Obj fn ='',E11.4,'
     '        //''' at params:'',12E10.3,(/10E10.3))')
     '        FUNC,(PBOPTI(noopti),noopti=1,NTOPTI)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
!PJH 9-APR-1992 added following:
          ELSE IF(.NOT.(KTYP26.EQ.2.AND.KTYP27.EQ.5)) THEN !not data fitting
            WRITE(OP_STRING,'('' Obj fn ='',E11.4,'
     '        //''' at params:'',12E10.3,(/10E10.3))')
     '        FUNC,(PBOPTI(noopti),noopti=1,NTOPTI)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF

        ELSE IF(NTSTEP.GT.0) THEN !evaluate obj. at stepped parameters
          nostep=0
          IF(NTOPTI.GE.1) THEN
            NTSTEP1=NTSTEP
            P1STEP=(PMAX(1)-PMIN(1))/DBLE(NTSTEP1-1)
          ELSE
            NTSTEP1=1
          ENDIF
          IF(NTOPTI.GE.2) THEN
            NTSTEP2=NTSTEP
            P2STEP=(PMAX(2)-PMIN(2))/DBLE(NTSTEP2-1)
          ELSE
            NTSTEP2=1
          ENDIF
          IF(NTOPTI.GE.3) THEN
            NTSTEP3=NTSTEP
            P3STEP=(PMAX(3)-PMIN(3))/DBLE(NTSTEP3-1)
          ELSE
            NTSTEP3=1
          ENDIF
          DO n1step=1,NTSTEP1
            IF(NTOPTI.GE.1) THEN
              !Update material params from optim.n params
              DO noelem=1,NEELEM(0,1)
                ne=NEELEM(noelem,1)
                CE(NMNO(1,1,nx),ne,nx)=PMIN(1)+DBLE(n1step-1)*P1STEP
              ENDDO
            ENDIF
            DO n2step=1,NTSTEP2
              IF(NTOPTI.GE.2) THEN
                !Update material params from optim.n params
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  CE(NMNO(1,2,nx),ne,nx)=PMIN(2)+DBLE(n2step-1)*P2STEP
                ENDDO
              ENDIF
              DO n3step=1,NTSTEP3
                IF(NTOPTI.GE.3) THEN
                  !Update material params from optim.n params
                  DO noelem=1,NEELEM(0,nr)
                    ne=NEELEM(noelem,nr)
                    CE(NMNO(1,3,nx),ne,nx)=PMIN(3)+DBLE(n3step-1)*P3STEP
                  ENDDO
                ENDIF
                nostep=nostep+1

                CALL OBJFUN(IBT,IDO,INP,LD,LGE,NAN,NBH,NBHF,NBJ,NBJF,
     '            NEELEM,NFF,NFFACE,NGAP,NHE(1,nx),
     '            NHP(1,nr,nx),NKB,NKEF,NKH,NKHE,NKJE,NLL,NLNO(1,nx),
     '            NNB,NNF,NNL,NONL(1,nx),NONY(0,1,1,nr,nx),NP_INTERFACE,
     '            NP1OPT,NPF,NPL,NPNE,NPNODE,NPNY(0,1,0,nx),nr,NRE,NSB,
     '            NVHE,NVHP(1,1,1,nr),NVJE,NW(1,1,nx),nx,NXI,
     '            NYNE,NYNO(0,1,1,nr,nx),NYNP,NYNR(0,0,1,nr,nx),PAOPTY,
     '            Z_CONT_LIST,CE(1,1,nx),CG,CGE(1,1,1,nx),
     '            CP(1,1,nx),
     '            CURVCORRECT,DL,FEXT,FGRAD,PAOPTI,PG,RE1,RESIDM,RESJAC,
     '            RG,SE,WG,WU,XA,XE,XG,XID,XIG,XP,YG,YGF,
     '            YP(1,1,nx),ZA,
     '            ZA1,Z_CONT,ZD,ZE,ZE1,ZG,ZP,ZP1,FIX(1,1,nx),
     '            ERROR,*9999)

                IF(NTOPTI.EQ.1) THEN
                  WRITE(OP_STRING,'('' Step '',I5,'' Obj fn = '',E11.4,'
     '              //''' at param:'',E11.3)')
     '              nostep,FUNC,CE(NMNO(1,1,nx),1,nx)
                ELSE IF(NTOPTI.EQ.2) THEN
                  WRITE(OP_STRING,'('' Step '',I5,'' Obj fn = '',E11.4,'
     '              //''' at params:'',2E11.3)')
     '              nostep,FUNC,CE(NMNO(1,1,nx),1,nx),
     '              CE(NMNO(1,2,nx),1,nx)
                ELSE IF(NTOPTI.EQ.3) THEN
                  WRITE(OP_STRING,'('' Step '',I5,'' Obj fn = '',E11.4,'
     '              //''' at params:'',3E11.3)') nostep,FUNC,
     '              CE(NMNO(1,1,nx),1,nx),CE(NMNO(1,2,nx),1,nx),
     '              CE(NMNO(1,3,nx),1,nx)
                ENDIF
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDDO !n3step
            ENDDO !n2step
          ENDDO !n1step
        ENDIF

        IF(KTYP26.EQ.1) THEN !Return material params to original values
          DO noopti=1,NMNO(1,0,nx)
            nm=NMNO(1,noopti,nx)
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              CE(nm,ne,nx)=TEMP(noopti)
            ENDDO
          ENDDO
        ENDIF

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('EVOBJE')
      RETURN
 9999 CALL ERRORS('EVOBJE',ERROR)
      CALL EXITS('EVOBJE')
      RETURN 1
      END


