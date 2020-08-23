      SUBROUTINE EVCONT(IBT,IDO,NBJ,NEEDCON,NEELEM,NEL,NKH,
     '  NLL,NLNO,NNL,NONL,NONY,NPL,NPNE,NPNODE,
     '  NPNY,NRLIST,NVHP,NVJL,NVJP,NXLIST,NYNO,NYNP,
     '  PAOPTY,CM,CJACM,CONTR,CONJAC,DL,PAOPTI,SE,XP,STRING,ERROR,*)

C#### Subroutine: EVCONT
C###  Description:
C###    EVCONT evaluates constraints.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),NBJ(NJM,NEM),
     '  NEEDCON(*),NEELEM(0:NE_R_M,0:NRM),NEL(0:NELM,NLM),
     '  NKH(NHM,NPM,NCM,0:NRM),NLL(12,NEM),NLNO(NOPM,NXM),
     '  NNL(0:4,12,NBFM),NONL(NLM,NXM),NONY(0:NOYM,NYM,NRCM,0:NRM,NXM),
     '  NPL(5,0:3,NLM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),NRLIST(0:NRM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJL(4,NJM,NLM),NVJP(NJM,NPM),NXLIST(0:NXM),
     '  NYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),PAOPTY(NOPM)
      REAL*8 CM(*),CJACM(NCOM,*),CONTR(*),CONJAC(NCOM,*),DL(3,NLM),
     '  PAOPTI(*),SE(NSM,NBFM,NEM),XP(NKM,NVM,NJM,NPM)
      CHARACTER STRING*(MXCH),ERROR*(*)
!     Local Variables
      INTEGER IBEG,IEND,N3CO,no_nrlist,nr,nx,nxc
      CHARACTER FILE*(MXCH),PARAMTYPE*20
      LOGICAL ABBREV,ALL_REGIONS,CBBREV,OPFILE

      CALL ENTERS('EVCONT',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM evaluate constraints<;FILENAME>
C###  Parameter:        <wrt PARAMTYPE[IPOPTI_params]>
C###    Specify the parameter type.
C###  Parameter:        <region (#s/all)[1]>
C###    Specify the region numbers to evaluate.
C###  Parameter:        <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Evaluates constraints for optimisations.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<wrt PARAMTYPE[IPOPTI_params]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','EVCONT',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN
          FILE=COQU(noco,1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          OPFILE=.TRUE.
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opcont','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
          IOFI=IOOP
        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)

        IF(CBBREV(CO,'WRT',2,noco+1,NTCO,N3CO)) THEN
          IF(ABBREV(CO(N3CO+1),'DATA_FITTING',1)) THEN
            PARAMTYPE='DATA_FITTING'
          ELSE IF(ABBREV(CO(N3CO+1),'IPOPTI_params',1)) THEN
            PARAMTYPE='IPOPTI_params'
          ELSE
            CO(noco+1)='?'
            GO TO 1
          ENDIF
        ELSE
          PARAMTYPE='IPOPTI_params'
        ENDIF

C CPB 2/6/94 Need to put in sparse confun for minos
C CPB 8/6/94 Adding NX_LOC
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_OPTI,ERROR,*9999)
        CALL ASSERT(nx.GT.0,
     '    '>>No nx defined for this optimisation class',
     '    ERROR,*9999)

        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          IF(PARAMTYPE(1:13).EQ.'IPOPTI_params') THEN
            CALL CONFUN(IBT,IDO,NBJ,NEEDCON,NEELEM,NEL,NKH,
     '        NLL,NLNO(1,nx),NNL,NONL(1,nx),NONY(0,1,1,nr,nx),
     '        NPL,NPNE,NPNODE,NPNY(0,1,0,nx),nr,NVHP,NVJL,NVJP,
     '        NYNO(0,1,1,nr,nx),NYNP,PAOPTY,CM,CJACM,CONTR,CONJAC,
     '        DL,PAOPTI,SE,XP,PARAMTYPE,ERROR,*9999)
          ELSE IF(PARAMTYPE(1:12).EQ.'DATA_FITTING') THEN
            CALL CONFUN(IBT,IDO,NBJ,NEEDCON,NEELEM,NEL,NKH,
     '        NLL,NLNO(1,nx),NNL,NONL(1,nx),NONY(0,1,1,nr,nx),
     '        NPL,NPNE,NPNODE,NPNY(0,1,0,nx),nr,NVHP,NVJL,NVJP,
     '        NYNO(0,1,1,nr,nx),NYNP,PAOPTY,CM,CJACM,CONTR,CONJAC,DL,
     '        PAOPTI,SE,XP,PARAMTYPE,ERROR,*9999)
          ENDIF
        ENDDO !no_nrlist

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('EVCONT')
      RETURN
 9999 CALL ERRORS('EVCONT',ERROR)
      CALL EXITS('EVCONT')
      RETURN 1
      END


