      SUBROUTINE DEINIT(IBT,IDO,INP,ITHRES,NAN,NBH,NBHF,NBJ,NBJF,NEELEM,
     '  NEL,NELIST,NENP,NENQ,NFF,NGAP,NHE,NHP,NHQ,NKEF,NKH,NKHE,NKJE,
     '  NLL,NNF,NNL,NODENVCB,NPL,NPLIST,NP_INTERFACE,NPF,NPNE,NPNF,
     '  NPNODE,NPNY,NQET,NQLIST,NQNE,NQS,NRLIST,NTIME_POINTS,NTIME_NR,
     &  NVCB,NVCNODE,NVHE,NVHF,NVHP,NVJE,NVJF,NW,NWQ,NXI,NXLIST,NYNE,
     &  NYNP,NYNQ,NYNR,TV_BC_SET,AQ,BBM,CE,CG,CGE,CP,CQ,CURVCORRECT,DF,
     &  DL,PG,RCQS,RG,SE,SF,THRES,WG,XA,XE,XG,XIG,XP,XQ,YG,YP,YQ,YQS,ZA,
     &  ZE,ZG,ZP,FIX,FIXQ,STRING,TIME_VARIABLE_NAMES,ERROR,*)

C#### Subroutine: DEINIT
C###  Description:
C###    DEINIT defines initial conditions with prompted input or from
C###    filename.ipinit (or filename.irinit for coupled problems).

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'anal00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'solv00.cmn'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ITHRES(3,NGM,NEM),NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),
     '  NBHF(NHM,NCM,NFM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NEELEM(0:NE_R_M,0:NRM),NEL(0:NELM,NLM),NELIST(0:NEM),
     '  NENP(NPM,0:NEPM,0:NRM),NENQ(0:8,NQM),NFF(6,NEM),NGAP(NIM,NBM),
     '  NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),NHQ(NRM,NXM),
     '  NKEF(0:4,16,6,NBFM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),
     '  NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),NODENVCB(NVCBM),NPF(9,NFM),
     '  NPL(5,0:3,NLM),NPLIST(0:NPM),NP_INTERFACE(0:NPM,0:3),
     '  NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),NQET(NQSCM),NQLIST(0:NQM),
     &  NQNE(NEQM,NQEM),NQS(NEQM),NRLIST(0:NRM),
     &  NTIME_POINTS(NTIMEVARSM),NTIME_NR(0:NTIMEVARSM,NRM),
     &  NVCB(-1:3,NVCBM),NVCNODE(2,NP_R_M),NVHE(NNM,NBFM,NHM,NEM),
     &  NVHF(NNM,NBFM,NHM),NVHP(NHM,NPM,NCM,0:NRM),
     &  NVJE(NNM,NBFM,NJM,NEM),NVJF(NNM,NBFM,NJM),NW(NEM,3,NXM),
     &  NWQ(8,0:NQM,NAM),NXI(-NIM:NIM,0:NEIM,0:NEM),NXLIST(0:NXM),
     &  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     &  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),NYNQ(NHM,NQM,0:NRCM,NXM),
     &  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),TV_BC_SET(0:NIQSM,0:NQM)
C SMAR009 19/01/99 removed NODENVC(NVCM),
      REAL*8 AQ(NMAQM,NQM),BBM(2,NEM),CE(NMM,NEM,NXM),CG(NMM,NGM),
     '  CGE(NMM,NGM,NEM,NXM),CP(NMM,NPM,NXM),CQ(NMM,NQM,NXM),
     '  CURVCORRECT(2,2,NNM,NEM),DF(NFM),DL(3,NLM),PG(NSM,NUM,NGM,NBM),
     '  RCQS(NQRM),RG(NGM),SE(NSM,NBFM,NEM),SF(NSM,NBFM),
     '  THRES(3,NGM,NEM),WG(NGM,NBM),XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XG(NJM,NUM),XIG(NIM,NGM,NBM),XP(NKM,NVM,NJM,NPM),XQ(NJM,NQM),
     '  YG(NIYGM,NGM,NEM),YP(NYM,NIYM,NXM),YQ(NYQM,NIQM,NAM,NXM),
     '  YQS(NIQSM,NQM),ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZG(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(*),
     '  TIME_VARIABLE_NAMES(NTIMEVARSM)*(*)
      LOGICAL FIX(NYM,NIYFIXM,NXM),FIXQ(NYQM,NIYFIXM,NXM)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,IPFILE,N3CO,
     '  no_nrlist,nr,nx,nxc
!     Functions      
      REAL*8 CONC_EVAL,RFROMC
      CHARACTER FILE*(MXCH),STATUS*3

CC AJPs - 191297 - RGB
      LOGICAL ACTIVATION,ALL_REGIONS,BC,CALCU,CBBREV,FILIO,
     '  FIRST_TIME,FIX_ZERO,FLOW_FIXED,GENER,LUNG_MESH,MOUSE,NOFIX,
     '  PRESSURE_FIXED,UPDATE_REF
CC AJPe

      CALL ENTERS('DEINIT',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define initial;d/l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Define initial and boundary conditions for the problem. These
C###    are read from or written to the file FILENAME.ipinit in the
C###    directory PATH. For a static problem, the essential and flux
C###    boundary conditions are defined. For a dynamic problem, the
C###    inital solution is also set up.
C###  Parameter:      <update_reference>
C###    Copy the current solution from YP(ny,1) into YP(ny,10) for
C###    use as the reference state for FE50 cavity problems.
C###  Parameter:      <fix/nofix[fix]>
C###    Specify whether to initialize (elasticity problems) or set
C###    (other problems) the FIX array (fix) or not (nofix).
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//';d/g/l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<update_reference>'
        OP_STRING(3)=BLANK(1:15)//'<fix/nofix[fix]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(5)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define initial;c lung
C###  Description:
C###    Define initial and boundary conditions for time dependent
C###    pulmonary problem.
C###  Parameter:      <temp_in>
C###    Specify the inspired temperature.
C###  Parameter:      <rh_in>
C###    Specify the inspired relative humidity.
C###  Parameter:      <core_temp>
C###    Specify the body core temperature.
C###  Parameter:      <mouth_temp>
C###    Specify the surrounding temperature at the model entrance.
C###  Parameter:      <sf_rate>
C###    Specify the rate of replenishment of airway surface liquid.  The
C###    default is a maximal value from the literature: 2 microL/cm^2/min.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.

        OP_STRING(1)=STRING(1:IEND)//';c lung'
        OP_STRING(2)=BLANK(1:15)//'<temp_in [30deg]>'
        OP_STRING(3)=BLANK(1:15)//'<rh_in [100%]>'
        OP_STRING(4)=BLANK(1:15)//'<core_temp [37deg]>'
        OP_STRING(5)=BLANK(1:15)//'<mouth_temp [core_temp]>'
        OP_STRING(6)=BLANK(1:15)//'<sf_rate [2microL/cm^2/min]>'
        OP_STRING(7)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define initial;g
C###  Description:
C###    Generate the initial and boundary conditions. For a static
C###    problem, the essential and flux boundary conditions are defined.
C###    For a dynamic problem, the initIal solution is also set up.
C###  Parameter:      <fix_zero>
C###    Specify both the flux and zero potential for the zero_potential
C###    point.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//';g'
        OP_STRING(2)=BLANK(1:15)//'<fix_zero>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
C#### Command: FEM define initial;g activation
C###  Description:
C###    Generate the initial and boundary conditions from an analytic
C###    activation sequence.
C###  Parameter:      <fix/nofix[fix]>
C###    Specify whether to set the FIX array (fix) or not (nofix). If
C###    not the values will still be stored in YP.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//';g activation'
        OP_STRING(2)=BLANK(1:15)//'<fix/nofix[fix]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
C---------------------------------------------------------------------


      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEINIT',ERROR,*9999)
      ELSE
        IPFILE=3 !is input file version number on 21-May-1993

        CALL PARSE_QUALIFIERS('CDGLPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)

        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

        CALL ASSERT(CALL_EQUA,'>>Problem type not defined',ERROR,*9999)
        CALL ASSERT(CALL_MATE,'>>Material parameters not defined',
     '    ERROR,*9999)

C CPB 8/6/94 Adding NX_LOC
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        IF(GENER) CALL ASSERT(CALL_ANAL,
     '    '>>Define analytic problem first',ERROR,*9999)

        IF(CBBREV(CO,'UPDATE_REFERENCE',1,noco+1,NTCO,N3CO)) THEN
          UPDATE_REF=.TRUE.
        ELSE
          UPDATE_REF=.FALSE.
        ENDIF

        IF(CBBREV(CO,'FIX_ZERO',1,noco+1,NTCO,N3CO)) THEN
          FIX_ZERO=.TRUE.
        ELSE
          FIX_ZERO=.FALSE.
        ENDIF

C LKC 12-JUL-2000 Allow read in ipinit to be nofixed
        ACTIVATION=.FALSE.
        IF(CBBREV(CO,'ACTIVATION',1,noco+1,NTCO,N3CO)) THEN
          ACTIVATION=.TRUE.
        ENDIF

        NOFIX=.FALSE.
        IF(CBBREV(CO,'NOFIX',3,noco+1,NTCO,N3CO)) THEN
          NOFIX=.TRUE.
        ENDIF

        IF(CALCU)THEN
          LUNG_MESH=.FALSE.
          IF(CBBREV(CO,'LUNG',4,noco+1,NTCO,N3CO)) THEN
c            CALL ASSERT(ITYP5(nr,nx).EQ.2.AND.ITYP2(nr,nx).EQ.11,
c     '        '>>Cannot use for this problem type',ERROR,*9999)
            LUNG_MESH=.TRUE.
            !set default values
            nr=NRLIST(1)
c            nx=NXLIST(1)
            CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '        ERROR,*9999)
C           Parse node(s) for applying boundary condition            
            CALL PARSE_NODES(NPNODE,NPLIST,N3CO,NRLIST,NTCO,CO,
     '        ERROR,*9999)
            BC=.FALSE.
            IF(CBBREV(CO,'INITIAL',4,noco+1,NTCO,N3CO))THEN
              CONC_INIT=RFROMC(CO(N3CO+1))
            ENDIF
            IF(CBBREV(CO,'FLUX',4,noco+1,NTCO,N3CO))THEN
              BC=.TRUE.
              FLUX_BC=.TRUE.
            ELSE IF(CBBREV(CO,'FIXED',4,noco+1,NTCO,N3CO))THEN
              BC=.TRUE.
              FLUX_BC=.FALSE.
C             Get fixed value(s) 
              IF(ITYP3(nr,nx).EQ.1)THEN !inert gas mixing
                IF(CBBREV(CO,'CONCENTRATION',4,noco+1,NTCO,N3CO))THEN
                  CONC_IN=RFROMC(CO(N3CO+1))
                ELSE
                  CONC_IN=1.D0
                ENDIF
                IF(.NOT.ADD) CONC_INIT=0.d0 !AJS - CONC_INIT was being overwritten when adding initial conditions
              ELSE IF(ITYP3(nr,nx).EQ.2)THEN !water and heat transfer
                IF(CBBREV(CO,'TEMPERATURE',4,noco+1,NTCO,N3CO))THEN
                  TEMP_IN=RFROMC(CO(N3CO+1))+273.15d0 !K
                ELSE
                  TEMP_IN=30.d0+273.15d0 !K
                ENDIF
                IF(CBBREV(CO,'RH',2,noco+1,NTCO,N3CO))THEN
                  RH_IN=RFROMC(CO(N3CO+1))
                ELSE 
                  RH_IN=100.d0 !%
                ENDIF
                CORE_TEMP=37.d0+273.15d0 !K
                SF_RATE=2.d0 !microL/cm^2/min
                IF(CBBREV(CO,'CORE_TEMP',4,noco+1,NTCO,N3CO))
     '            CORE_TEMP=RFROMC(CO(N3CO+1))+273.15d0 !K
                MOUTH_TEMP=CORE_TEMP !K
                IF(CBBREV(CO,'MOUTH_TEMP',3,noco+1,NTCO,N3CO))
     '            MOUTH_TEMP=RFROMC(CO(N3CO+1))+273.15d0 !K
                IF(CBBREV(CO,'SF_RATE',2,noco+1,NTCO,N3CO))
     '            SF_RATE=RFROMC(CO(N3CO+1))
                AH_IN=CONC_EVAL(TEMP_IN)*RH_IN/100.d0
              ELSE IF(ITYP3(nr,nx).EQ.4)THEN !simple P-R-F
                PRESSURE_FIXED=.FALSE.
                FLOW_FIXED=.FALSE.
                IF(CBBREV(CO,'PRESSURE',4,noco+1,NTCO,N3CO))
     &            PRESSURE_FIXED=.TRUE.
                IF(CBBREV(CO,'FLOW',4,noco+1,NTCO,N3CO))
     &            FLOW_FIXED=.TRUE.
              ENDIF
            ELSE IF(CBBREV(CO,'FIELD',3,noco+1,NTCO,N3CO))THEN
              BC=.TRUE.
              FLOW_FIXED=.TRUE.        !Used for initialising in INIT_2_11
              PRESSURE_FIXED=.FALSE.
              FLUX_BC=.FALSE.
            ENDIF
            
c            CALL ASSERT(BC,'>>Must specify boundary condition',ERROR,
c     &        *9999)
              
          ELSE
            WRITE(OP_STRING,'('' Not implemented '')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF !LUNG
        ENDIF !CALCU


        IF(FILIO) THEN
          IF(IOTYPE.EQ.3) THEN !writing ipinit
            CALL ASSERT(CALL_INIT,'>>Must define initial first',
     '        ERROR,*9999)
          ENDIF
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            IF(.NOT.GENER) CALL OPEN_FILE(ALL_REGIONS,IPFILE,
     '        FILE,'init',STATUS,ERR,ERROR,*9999)

            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              IF(ALL_REGIONS.AND..NOT.GENER) CALL PROMPT_REGION_ALL(nr,
     '          ERROR,*9999)

              IF(GENER) CALL ASSERT(ANAL_CHOICE(nr).NE.0,
     '          '>>Analytic choice is not set for this region',
     '          ERROR,*9999)

C LKC 20-JAN-98 added NQXI
C MLB 28-APR-98 removed NQXI not used or passed from fem_dynam

CC AJPs - 191297 - RGB
              CALL IPINIT(IBT,IDO,INP,ITHRES,NAN,NBH,NBHF,NBJ,NBJF,
     '          NEELEM,NEL,NENP,NENQ,NFF,NGAP,NHE(1,nx),
     '          NHP(1,nr,nx),NHQ(nr,nx),NKEF,NKH(1,1,1,nr),NKHE,NKJE,
     '          NLL,NNF,NNL,NODENVCB,NPF,NPL,NPLIST,NP_INTERFACE,NPNE,
     '          NPNF,NPNODE,NPNY,NQET,NQLIST,NQNE,NQS,nr,
     '          NTIME_NR,NVCB,NVCNODE,NVHE,NVHF,NVHP(1,1,1,nr),NVJE,
     '          NVJF,NW(1,1,nx),NWQ,nx,NXLIST,NXI,NYNE,NYNP,
     '          NYNQ(1,1,0,nx),NYNR(0,0,1,0,nx),TV_BC_SET,AQ,CE(1,1,nx),
     '          CG,CGE(1,1,1,nx),CP(1,1,nx),CQ(1,1,nx),CURVCORRECT,DF,
     '          DL,PG,RCQS,RG,SE,SF,THRES,WG,XA,XE,XG,XIG,XP,XQ,YG,
     '          YP(1,1,nx),YQ,YQS,ZA,ZE,ZG,ZP,TIME_VARIABLE_NAMES,
     '          ACTIVATION,ALL_REGIONS,FIX(1,1,nx),FIXQ,FIX_ZERO,GENER,
     '          NOFIX,UPDATE_REF,ERROR,*9999)
C  SMAR009 18/01/99 removed NODENVC,
CC AJPe
            ENDDO !no_nrlist

            IF(.NOT.GENER) CALL CLOSEF(IFILE,ERROR,*9999)
            IF(BACKUP) IOTYPE=2
          ENDDO

        ELSE IF(MOUSE) THEN
        ELSE IF(CALCU)THEN
          IF(LUNG_MESH)THEN
            !ITYP5(nr,nx)=2; ITYP2(nr,nx)=11
            nr=NRLIST(1)
            CALL INIT_2_11(NBJ,NEELEM(0,nr),NELIST,NENP,NPLIST,NPNE,
     &        NPNODE(0,nr),nr,NTIME_POINTS,NVJE,nx,NXI,NYNE,
     &        NYNP(1,1,1,1,0,1,nr),BBM,CE(1,1,nx),CP(1,1,nx),XP,
     &        YP(1,1,nx),BC,FIX(1,1,nx),FLOW_FIXED,PRESSURE_FIXED,
     &        ERROR,*9999)
          ENDIF
        ENDIF
        IF(IOTYPE.NE.3) THEN !not writing the ipinit file
          DO nr=1,NRT
            ASSEMBLE_SOLUTION(nr,nx)=.FALSE.
          ENDDO !nr
          CALL_INIT=.TRUE.
          IF(CALL_SOLV) THEN
C           May need to 'define solve' again to set mappings.
C           Define solve only needed if type of b.c.s is changed.
            WRITE(OP_STRING,'('' >>Warning: May need to '
     '        //'define solve again.'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
      ENDIF

      CALL EXITS('DEINIT')
      RETURN
 9999 CALL ERRORS('DEINIT',ERROR)
      CALL EXITS('DEINIT')
      RETURN 1
      END


