      SUBROUTINE SOLVE(CELL_ICQS_VALUE,DIPOLE_CEN_NTIME,
     '  DIPOLE_DIR_NTIME,IBT,ICQS,ICQS_SPATIAL,IDO,IICQS_SPATIAL,
     '  IRCQS_SPATIAL,INP,ISC_GD,ISC_GKK,ISC_GM,ISC_GMM,ISR_GD,ISR_GKK,
     '  ISR_GM,ISR_GMM,ISEG,ISELNO,ISFIBR,ISFIEL,ISLINE,ISLINO,ISNONO,
     '  ISTATE,IUSER,LGE,MXI,NAN,NAQ,NBH,NBHF,NBJ,NBJF,NDET,
     '  NDIPOLES,NEELEM,NELIST,NENFVC,NENP,NENQ,NEP,NFF,NFFACE,
     '  NFVC,NGAP,NHE,NHP,NHQ,NKB,NKEF,NKH,NKHE,NKJE,NLATNE,NLATNQ,
     '  NLATPNQ,NLL,NLLIST,NLQ,NNB,NNF,NODENVC,NODENVCB,NONY,NORD,NMNO,
     &  NP_INTERFACE,NPB,NPF,NPL,NPLIST,NPLIST2,NPNE,NPNODE,NPNY,NQET,
     &  NQGP,NQGP_PIVOT,NQGW,NQLIST,NQNE,NQNLAT,NQNP,NQNY,NQS,NQSCNB,
     &  NQXI,NRE,NRLIST,NRLIST2,NSB,NTIME_INTERP,NTIME_POINTS,NTIME_NR,
     &  NVCB,NVCNODE,NVHE,NVHP,NVJE,NVJP,NW,NWQ,NXI,NXLIST,NXQ,NYNE,
     &  NYNO,NYNP,NYNQ,NYNR,NYQNR,TV_BC_SET,Z_CONT_LIST,ACINUS,AQ,BBM,
     &  CE,CELL_RCQS_VALUE,CG,CGE,CONY,CP,CQ,CURVCORRECT,CYNO,DET,
     &  DIPOLE_CEN,DIPOLE_DIR,DL,DNUDXQ,DRDN,DRDNO,DXDXIQ,DXDXIQ2,ED,
     &  EIGVAL,EIGVEC,EM,ER,ES,FEXT,GCHQ,GD,GKK,GM,GMM,GR,GRR,GUQ,
     &  PAOPTI,PG,PMIN,PMAX,PROPQ,R,RAD,RCQS,RCQS_SPATIAL,RD,RE1,RESID,
     &  RESJAC,RG,RHS,SE,TIME_VALUES,USER,VC,VC_INIT,WG,XA,XAB,XC,XE,XG,
     &  XG1,XIG,XIP,XIQ,XN,XNFV,XN_GRAD,XO,XP,XQ,XR,XR_GRAD,YG,YGF,YP,
     &  YQ,YQS,ZA,ZA1,Z_CONT,ZE,ZG,ZNFV,ZP,ZP1,CSEG,STRING,FIX,FIXQ,
     &  ERROR,*)

C#### Subroutine: SOLVE
C###  Description:
C###    SOLVE calls up solution routines.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'coup00.cmn'
      INCLUDE 'gen000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'gks000.cmn'
      INCLUDE 'ipma50.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'load00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'marc00.cmn'
      INCLUDE 'nonl00.cmn'      
      INCLUDE 'ptr00.cmn'
      INCLUDE 'quas00.cmn'
      INCLUDE 'solv00.cmn'
      INCLUDE 'time02.cmn'

!     Parameter List
      INTEGER CELL_ICQS_VALUE(NQIM,NQVM),
     '  DIPOLE_CEN_NTIME(NDIPOLEM,NRM,NXM),
     '  DIPOLE_DIR_NTIME(NDIPOLEM,NRM,NXM),IBT(3,NIM,NBFM),ICQS(NQIM),
     '  ICQS_SPATIAL(NQISVM,NQM),IDO(NKM,NNM,0:NIM,NBFM),
     '  IICQS_SPATIAL(0:NQISVM,NQVM),IRCQS_SPATIAL(0:NQRSVM,NQVM),
     '  INP(NNM,NIM,NBFM),ISC_GD(NISC_GDM),ISC_GKK(NISC_GKKM,NXM),
     '  ISC_GM(NISC_GMM),ISC_GMM(NISC_GMMM),ISR_GD(NISR_GDM),
     '  ISR_GKK(NISR_GKKM,NXM),ISR_GM(NISR_GMM),ISR_GMM(NISR_GMMM),
     '  ISEG(*),ISELNO(NWM,NEM),ISFIBR(NWM,NEM,NGRSEGM),ISFIEL(NWM,NEM),
     '  ISLINE(NWM,2*NGRSEGM),ISLINO(NWM),ISNONO(NWM,NPM),ISTATE(*),
     '  IUSER(*),LGE(NHM*NSM,NRCM),MXI(2,NEM),
     '  NAN(NIM,NAM,NBFM),NAQ(NQM,NAM),NBH(NHM,NCM,NEM),
     '  NBHF(NHM,NCM,NFM),NBJ(NJM,NEM),NBJF(NJM,NFM),NDET(NBFM,0:NNM),
     '  NDIPOLES(NRM,NXM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NENFVC(0:NFVCM,NFVM),NENP(NPM,0:NEPM,0:NRM),
     '  NENQ(0:8,NQM),NEP(NPM),NFF(6,NEM),NFFACE(0:NF_R_M,NRM),
     &  NFVC(2,0:NFVCM,NVCM),NGAP(NIM,NBM),NHE(NEM,NXM),
     '  NHP(NPM,0:NRM,NXM),NHQ(NRM,NXM),NKB(2,2,2,NNM,NBFM),
     '  NKEF(0:4,16,6,NBFM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKJE(NKM,NNM,NJM,NEM),NKHE(NKM,NNM,NHM,NEM),
     '  NLATNE(NEQM+1),NLATNQ(NEQM*NQEM),NLATPNQ(NQM),
     '  NLL(12,NEM),NLLIST(0:NLM),NLQ(NQM),NNB(4,4,4,NBFM),
     &  NNF(0:17,6,NBFM),NODENVC(NVCM),NODENVCB(NVCBM),
     &  NONY(0:NOYM,NYM,NRCM,0:NRM,NXM),NORD(5,NE_R_M),
     &  NMNO(1:2,0:NOPM,NXM),NP_INTERFACE(0:NPM,0:3),NPB(0:NP_R_M,5),
     &  NPF(9,NFM),NPL(5,0:3,NLM),NPLIST(0:NPM),NPLIST2(0:NPM),
     &  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     &  NPNY(0:6,NYM,0:NRCM,NXM),NQET(NQSCM),NQGP(0:NQGM,NQM),
     &  NQGP_PIVOT(NQGM,NQM),NQLIST(0:NQM),NQNLAT(NEQM*NQEM),
     &  NQNE(NEQM,NQEM),NQNP(NPM),NQNY(2,NYQM,0:NRCM,NXM),NQS(NEQM),
     &  NQSCNB(NQSCM),NQXI(0:NIM,NQSCM),NRE(NEM),NRLIST(0:NRM),
     &  NRLIST2(0:NRM),NSB(NKM,NNM,NBFM),NTIME_INTERP(NTIMEVARSM),
     &  NTIME_POINTS(NTIMEVARSM),NTIME_NR(0:NTIMEVARSM,NRM),
     &  NVCB(-1:3,NVCBM),NVCNODE(2,NP_R_M),NVHE(NNM,NBFM,NHM,NEM),
     &  NVHP(NHM,NPM,NCM,0:NRM),NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),
     &  NW(NEM,3,NXM),
     &  NWQ(8,0:NQM,NAM),NXI(-NIM:NIM,0:NEIM,0:NEM),NXLIST(0:NXM),
     &  NXQ(-NIM:NIM,0:4,0:NQM,NAM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     &  NYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),
     &  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),NYNQ(NHM,NQM,0:NRCM,NXM),
     &  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),NYQNR(0:NYQM,0:NRCM,NCM,
     &  0:NRM,NXM),TV_BC_SET(0:NIQSM,0:NQM),Z_CONT_LIST(NDM,2,7)
      REAL*8 ACINUS(4,NEM),AQ(NMAQM,NQM),BBM(2,NEM),CE(NMM,NEM,NXM),
     '  CELL_RCQS_VALUE(NQRM,NQVM),CG(NMM,NGM),CGE(NMM,NGM,NEM,NXM),
     '  CONY(0:NOYM,NYM,NRCM,0:NRM,NXM),CP(NMM,NPM,NXM),CQ(NMM,NQM,NXM),
     '  CURVCORRECT(2,2,NNM,NEM),CYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),
     '  DET(NBFM,0:NNM,NGM,6),DIPOLE_CEN(4,0:NDIPTIMM,NDIPOLEM,NRM,NXM),
     '  DIPOLE_DIR(4,0:NDIPTIMM,NDIPOLEM,NRM,NXM),DL(3,NLM),
     '  DNUDXQ(3,3,NQM),DRDN(NGM),DRDNO(NGM,NKM),
     '  DXDXIQ(3,3,NQM),
     '  DXDXIQ2(3,3,NQM),ED(NHM*NSM,NHM*NSM),
     '  EIGVAL(NTM,2),EIGVEC(NOM,NTM,2),EM(NHM*NSM,NHM*NSM),ER(NHM*NSM),
     '  ES(NHM*NSM,NHM*NSM),FEXT(NIFEXTM,NGM,NEM),GCHQ(3,NQM),
     '  GD(NZ_GD_M),GKK(NZ_GKK_M,NXM),GM(NZ_GM_M),GMM(NZ_GMM_M),
     '  GR(NYROWM),GRR(NOM),GUQ(3,3,NQM),NQGW(NQGM,NQM),PAOPTI(*),
     '  PG(NSM,NUM,NGM,NBM),PMIN(*),PMAX(*),PROPQ(3,3,4,2,NQM,NXM),
     '  R(NOPM,*),RAD(NGM),
     '  RCQS(NQRM),RCQS_SPATIAL(NQRSVM,NQM),RD(NGM),RE1(NSM,NHM),
     '  RESID(*),RESJAC(NREM,*),RG(NGM),RHS(NQM),
     '  SE(NSM,NBFM,NEM),TIME_VALUES(2,0:NTIMEPOINTSM+1,NTIMEVARSM),
     '  USER(*),VC(0:NVCM),VC_INIT(2,NVCM),WG(NGM,NBM),XA(NAM,NJM,NEM),
     &  XAB(NORM,NEM),XC(*),XE(NSM,NJM),XG(NJM,NUM),XG1(NJM,NUM,NGM),
     &  XIG(NIM,NGM,NBM),XIP(NIM,NPM),XIQ(NIM,NQM),XN(NJM,NGM),
     &  XN_GRAD(NJM,NGM),XNFV(-(NJM+1):NJM,NFVM),XO(NOM,NXM),
     &  XP(NKM,NVM,NJM,NPM),XQ(NJM,NQM),XR(NJM,NGM),XR_GRAD(NJM,NGM),
     &  YG(NIYGM,NGM,NEM),YGF(NIYGFM,NGFM,NFM),YP(NYM,NIYM,NXM),
     &  YQ(NYQM,NIQM,NAM,NXM),YQS(NIQSM,NQM),ZA(NAM,NHM,NCM,NEM),
     &  ZA1(NAM,NHM,NCM,NEM),Z_CONT(NDM,2,67),ZE(NSM,NHM),ZG(NHM,NUM),
     '  ZNFV(NFVM),ZP(NKM,NVM,NHM,NPM,NCM),
     '  ZP1(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM),FIXQ(NYQM,NIYFIXM,NXM)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,IWK(6),N3CO,nonrlist,
     '  nrcounter,nr,nr_gkk,nrr,NTIW,nx,nx1,nx2,nx3,nxc,nzz
      INTEGER*4 WORK_PTR
      REAL*8 ALPHA,ERRMAX,RFROMC,ZE1(NSM,NHM)
      CHARACTER FILE*(MXCH)
      LOGICAL ALL_REGIONS,ASSEMBLE,CBBREV,ITER8,OPFILE,QUASI,
     '  WARM_START,GRAVITY_RESET

      nx1=0
      nx2=0
      nx3=0
      OPFILE=.FALSE.

      CALL ENTERS('SOLVE',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM solve<;FILENAME>
C###  Description:
C###    Solve boundary value problem with specified steps and
C###    iterations (if nonlinear).  The 'increment FACTOR' parameter
C###    allows you to modify previously specified dependent variable
C###    or flux/load increments.  The 'update' parameter will cause the
C###    specified workstations to be updated as the solution proceeds.
C###  Parameter:     <(separate/coupled)[separate]>
C###    Solve either a single region (separate) or a multiple region
C###    (coupled) problem.
C###  Parameter:     <(undeformed/deformed)[deformed]>
C###    Solve for either the deformed (deformed) coordinates or the undeformed
C###    (undeformed) coordinates 
C###  Parameter:     <step #[1]>
C###    Specify the number of steps taken for each increment for non
C###    linear problems.
C###  Parameter:     <iterate #[8]>
C###    Specifies the maximum number of iterations allowed to reach
C###    convergence on each step. For grid problems this invokes
C###    a bidomain iteration step.
C###  Parameter:     <increment FACTOR[1]>
C###    A highly nonlinear problem may wish to be solved in increments
C###    and the solution outputted at the end of each increment. The
C###    increments are defined as a factor of the overall solution
C###    procedure.
C###  Parameter:     <gravity_increment_reset (TRUE/FALSE)[FALSE]>
C###    TRUE resets the gravity increment to 0.0. FALSE does not change
C###    the gravity increment vector at all.
C###  Parameter:     <error LIMIT[1.0d-10]>
C###    Specifies the convergence tolerance required.
C###  Parameter:     <warm>
C###    For nonlinear FE problems which use the modified newton method
C###    the stiffness matrix is not recalculated on the first iteration.
C###  Parameter:     <update WSS[none]>
C###    Updates a specified window.
C###  Parameter:     <alpha #[0.0]>
C###    This parameter is used with the iterate option. It is
C###    alpha*grid + (1-alpha)*bem for updating grid right hand
C###    side vectors. Alpha is a real in the range [0-1].
C###  Parameter:     <region (#s/all)[1]> (separate only)
C###    Specify the region numbers to use. The "all" keyword will
C###    use all currently defined regions. The region list is set up in
C###    define couple for coupled problems
C###  Parameter:     <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<(separate/coupled)[separate]>'
        OP_STRING(3)=BLANK(1:15)//'<(undeformed/deformed)[deformed]>'
        OP_STRING(4)=BLANK(1:15)//'<step #[1]>'
        OP_STRING(5)=BLANK(1:15)//'<iterate #[8]>'
        OP_STRING(6)=BLANK(1:15)//'<increment FACTOR[1]>'
        OP_STRING(7)=BLANK(1:15)//'<gravity_increment_reset [false]>'
        OP_STRING(8)=BLANK(1:15)//'<error LIMIT[1.0d-10]>'
        OP_STRING(9)=BLANK(1:15)//'<warm>'
        OP_STRING(10)=BLANK(1:15)//'<update WSS[none]>'
        OP_STRING(11)=BLANK(1:15)//'<alpha #[0.0]>'
        OP_STRING(12)=BLANK(1:15)//'<region (#s/all)[1]>(separate only)'
        OP_STRING(13)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM solve<;FILENAME>
C###  Description:
C###    Solve time dependent problems with specified parameters.
C###  Parameter:     <(separate/coupled)[separate]>
C###    Solve either a single region (separate) or a multiple region
C###    (coupled) problem.
C###  Parameter:     <restart>
C###    Start solving again from the solution generated at the
C###    preceding time step.
C###  Parameter:     <from INITIAL_TIME#[T0]>
C###    Specifies the initial time for the solution procedure
C###  Parameter:     <to (FINAL_TIME#/for TIME_INTERVAL#)[T1]>
C###    Defines the final time or the length of the time interval.
C###  Parameter:     <delta_t STEP[DT]>
C###    Defines the time steps taken for the solution procedure.
C###  Parameter:     <at #)[0.0]>
C###    Defines the time for a static problem
C###  Parameter:     <update WSS[none]>
C###    Updates a specified window.
C###  Parameter:     <region (#s/all)[1]> (separate only)
C###    Specify the region numbers to use. The "all" keyword will
C###    use all currently defined regions. The region list is set up in
C###    define couple for coupled problems
C###  Parameter:     <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<(separate/coupled)[separate]>'
        OP_STRING(3)=BLANK(1:15)//'<restart>'
        OP_STRING(4)=BLANK(1:15)//'<from INITIAL_TIME[T0]>'
        OP_STRING(5)=BLANK(1:15)//'<to (FINAL_TIME#/for TIME_INTERVAL#)'
     '                          //'[T1]>'
        OP_STRING(6)=BLANK(1:15)//'<delta_t STEP[DT]>'
        OP_STRING(7)=BLANK(1:15)//'<at #[--]>'
        OP_STRING(8)=BLANK(1:15)//'<update WSS[none]>'
        OP_STRING(9)=BLANK(1:15)//'<reassemble>'
        OP_STRING(10)=BLANK(1:15)//'<region (#s/all)[1]>'
     '                           //' (separate only)'
        OP_STRING(11)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM solve<;FILENAME> activation
C###  Description:
C###    This command is for solving grid activation problems
C###    with more control from the command line.
C###  Parameter:     <(transmembrane/extracellular/both/none)[both]>
C###    Specify which of the bidomain equations to solve.
C###  Parameter:     <(update_vm/noupdate_vm)[update_vm]>
C###    The update Vm option is the normal one which stores the
C###    transmembrane potential at the end of each time step.
C###    The noupdate Vm is used when iterating and the solution
C###    which is to be stored for Vm at the current time has not
C###    yet been finalised.
C###  Parameter:     <static>
C###    This parameter is used for iterations and does not increase
C###    the current solution time.
C###  Parameter:     <(separate/coupled)[separate]>
C###    Specify whether the problem is a grid only problem (separate)
C###    or if it it a coupled grid/boundary element problem (coupled).
C###  Parameter:     <restart>
C###    Start solving again from the solution generated at the
C###    preceding time step.
C###  Parameter:     <to (FINAL_TIME#/for TIME_INTERVAL#)[T1]>
C###    Defines the final time or the length of the time interval.
C###  Parameter:     <cptype #[1]>
C###    Defines the final time or the length of the time interval.
C###  Parameter:     <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will
C###    use all currently defined regions.
C###  Parameter:     <class (#s)[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> activation'
        OP_STRING(2)=BLANK(1:15)//'<(transmembrane/extracellular/'
     '    //'both/none)[both]>'
        OP_STRING(3)=BLANK(1:15)//'<(update_vm/noupdate_vm)[update_vm]>'
        OP_STRING(4)=BLANK(1:15)//'<static>'
        OP_STRING(5)=BLANK(1:15)//'<(separate/coupled)[separate]>'
        OP_STRING(6)=BLANK(1:15)//'<restart>'
        OP_STRING(7)=BLANK(1:15)//'<to (FINAL_TIME#/for TIME_INTERVAL#)'
     '                          //'[T1]>'
        OP_STRING(8)=BLANK(1:15)//'<cptype #[1]>'
        OP_STRING(9)=BLANK(1:15)//'<region (#s/all)[1]> (separate only)'
        OP_STRING(10)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','SOLVE',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN
          FILE=COQU(noco,1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          OPFILE=.TRUE.
          IO4=IOFILE1
          CALL OPENF(IO4,'DISK',FILE(IBEG:IEND)//'.solve','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          IO4=IOOP
          OPFILE=.FALSE.
        ENDIF

        CALL ASSERT(CALL_MATE,'>>Need to define material',
     '    ERROR,*9999)
        CALL ASSERT(CALL_INIT,'>>Need to define initial',
     '    ERROR,*9999)
        CALL ASSERT(CALL_SOLV,'>>Need to define solve',ERROR,*9999)

        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
C CPB 8/6/94 Adding NX_LOC
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)
C news VJ 17Nov2004 Adding an option of solving for undeformed state or deformed (default)
C                   where you change the calculation of the element stiffness matrices to be 
C                   derivatives with respect to undeformed coordinates or deformed coordinates
        IF(CBBREV(CO,'UNDEFORMED',3,noco+1,NTCO,N3CO)) THEN
          KTYP5L=2
        ELSEIF(CBBREV(CO,'DEFORMED',3,noco+1,NTCO,N3CO)) THEN
          KTYP5L=1
        ELSE
          KTYP5L=1
        ENDIF
C newe VJ
        IF(CBBREV(CO,'WARM',3,noco+1,NTCO,N3CO)) THEN
          WARM_START=.TRUE.
        ELSE
          WARM_START=.FALSE.
        ENDIF

        IF(ISC_GK_PTR(nx).EQ.0) THEN
          CALL ALLOCATE_MEMORY(NISC_GKM*USE_SPARSE,0,INTTYPE,
     '      ISC_GK_PTR(nx),MEM_INIT,ERROR,*9999)
        ENDIF
        IF(ISR_GK_PTR(nx).EQ.0) THEN
          CALL ALLOCATE_MEMORY(NISR_GKM*USE_SPARSE,0,INTTYPE,
     '      ISR_GK_PTR(nx),MEM_INIT,ERROR,*9999)
        ENDIF
        IF(GK_PTR(nx).EQ.0) THEN
          CALL ALLOCATE_MEMORY(NZ_GK_M,1,DPTYPE,GK_PTR(nx),MEM_INIT,
     '      ERROR,*9999)
          WARM_START=.FALSE.
        ENDIF
        IF(ISC_GQ_PTR(nx).EQ.0) THEN
          CALL ALLOCATE_MEMORY(NISC_GQM*USE_SPARSE,0,INTTYPE,
     '      ISC_GQ_PTR(nx),MEM_INIT,ERROR,*9999)
        ENDIF
        IF(ISR_GQ_PTR(nx).EQ.0) THEN
          CALL ALLOCATE_MEMORY(NISR_GQM*USE_SPARSE,0,INTTYPE,
     '      ISR_GQ_PTR(nx),MEM_INIT,ERROR,*9999)
        ENDIF
        IF(GQ_PTR(nx).EQ.0) THEN
          CALL ALLOCATE_MEMORY(NZ_GQ_M,1,DPTYPE,GQ_PTR(nx),MEM_INIT,
     '      ERROR,*9999)
        ENDIF

        IF(CBBREV(CO,'ACTIVATION',3,noco+1,NTCO,N3CO)) THEN
          SOLVE_ACTIV=.TRUE.
          USE_SOLVE_ACTIV=.TRUE.
        ELSE
          SOLVE_ACTIV=.FALSE.
        ENDIF

        IF(.NOT.SOLVE_ACTIV) THEN
          IF(CBBREV(CO,'COUPLED',2,noco+1,NTCO,N3CO)) THEN
            CALL ASSERT(IS_COUPLED(nx),
     '        '>>Define solve for coupled problem',ERROR,*9999)
          ELSE
            CALL ASSERT(.NOT.IS_COUPLED(nx),
     '        '>>Define solve for separate problem',ERROR,*9999)
          ENDIF
        ENDIF

        IF(IS_COUPLED(nx)) THEN
          NRLIST(0)=COUP_NRLIST(0,nx)
          DO nonrlist=1,COUP_NRLIST(0,nx)
            NRLIST(nonrlist)=COUP_NRLIST(nonrlist,nx)
          ENDDO
        ELSE
          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,
     '      ERROR,*9999)
        ENDIF

        IF(CBBREV(CO,'STEP',1,noco+1,NTCO,N3CO)) THEN
          NTLOAD=IFROMC(CO(N3CO+1))
        ELSE
          NTLOAD=1
        ENDIF

        IF(CBBREV(CO,'ALPHA',3,noco+1,NTCO,N3CO)) THEN
          ALPHA=RFROMC(CO(N3CO+1))
        ELSE
          ALPHA=0.0d0
        ENDIF
C news VJ 21Jan2004: NTITER is not used for cellml problems
C                    therefore don't set default of 8 if cell
C                    class is being solved
        IF(CBBREV(CO,'ITERATE',2,noco+1,NTCO,N3CO)) THEN
          NTITER=IFROMC(CO(N3CO+1))
          ITER8=.TRUE.
        ELSE
          DO nonrlist=1,NRLIST(0)
            nrcounter=NRLIST(nonrlist)
            IF(.NOT.((ITYP5(nrcounter,nxc).EQ.2).AND.
     &        (ITYP2(nrcounter,nxc).EQ.9).AND.
     &        (ITYP19(nrcounter,nxc).EQ.1).AND.
     &        (ITYP3(nrcounter,nxc).EQ.10).AND.
     &        (KTYP33.EQ.6))) THEN !not time integration, electrical, user defined,CELLML 
              NTITER=8
            ENDIF
            ITER8=.FALSE.          
          ENDDO !nr
        ENDIF
C newe VJ
        IF(CBBREV(CO,'INCREMENT',2,noco+1,NTCO,N3CO)) THEN
          FACTOR=RFROMC(CO(N3CO+1))
        ELSE
          FACTOR=1.0d0
        ENDIF
C
C  OR 13-09-06
C     Included a flag to set the gravity increment to 0. If not reset
C     multitiple FEM solve commands would lead to linearly increasing
C     gravity "constant"
C
        IF(CBBREV(CO,'GRAVITY_INCREMENT_RESET',3,noco+1,NTCO,N3CO)) THEN
          GRAVITY_RESET=CBBREV(CO(N3CO+1),'TRUE',3,1,1,N3CO)
          IF (GRAVITY_RESET) THEN
            b_inc(1)=0.0d0
            b_inc(2)=0.0d0
            b_inc(3)=0.0d0
          ENDIF
        ELSE
          GRAVITY_RESET=.FALSE.
        ENDIF

        IF(CBBREV(CO,'RESTART',3,noco+1,NTCO,N3CO)) THEN
          RESTART=.TRUE.
        ELSE
          RESTART=.FALSE.
        ENDIF

        IF(CBBREV(CO,'REASSEMBLE',3,noco+1,NTCO,N3CO)) THEN
          REASSEMBLE=.TRUE.
        ELSE
          REASSEMBLE=.FALSE.
        ENDIF

C LKC 15-APR-2002 moving inside the 'AT' loop.
C
C        IF(CBBREV(CO,'FROM',2,noco+1,NTCO,N3CO)) THEN
C          TSTART=RFROMC(CO(N3CO+1))
C        ELSE
C          TSTART=T0
C        ENDIF
C
C        IF(CBBREV(CO,'TO',1,noco+1,NTCO,N3CO)) THEN
C          TFINISH=RFROMC(CO(N3CO+1))
C        ELSE IF(CBBREV(CO,'FOR',2,noco+1,NTCO,N3CO)) THEN
C          IF(RESTART) THEN
C            TFINISH=T_RESTART+RFROMC(CO(N3CO+1))
C          ELSE
C            TFINISH=TSTART+RFROMC(CO(N3CO+1))
C          ENDIF
C        ELSE
C          TFINISH=T1
C        ENDIF

C LKC 15-APR-2002 Sets the "start" time for a problem that is
C  solved as a static problem yet has a time dependent solution
C  Default is 0.d0. The variables are overwritten if a
C  non-static problem is solved ... ie FROM/TO

        IF(CBBREV(CO,'AT',2,noco+1,NTCO,N3CO)) THEN
          TSTART=RFROMC(CO(N3CO+1))
        ELSE
          TSTART=0.D0

          IF(CBBREV(CO,'FROM',2,noco+1,NTCO,N3CO)) THEN
            TSTART=RFROMC(CO(N3CO+1))
          ELSE
            TSTART=T0
          ENDIF

          IF(CBBREV(CO,'TO',1,noco+1,NTCO,N3CO)) THEN
            TFINISH=RFROMC(CO(N3CO+1))
          ELSE IF(CBBREV(CO,'FOR',2,noco+1,NTCO,N3CO)) THEN
            IF(RESTART) THEN
              TFINISH=T_RESTART(nx)+RFROMC(CO(N3CO+1))
            ELSE
              TFINISH=TSTART+RFROMC(CO(N3CO+1))
            ENDIF
          ELSE
            TFINISH=T1
          ENDIF
        ENDIF

        IF(CBBREV(CO,'DELTA_T',1,noco+1,NTCO,N3CO)) THEN
          DT=RFROMC(CO(N3CO+1))
C PM 26-JUL-01
c          TINCR=DT
        ELSE
C PJH 30Jun98 the default here should be the value set in define solve
!         DT=1.0D-1
C PM 26-JUL-01
          DT=TINCR
        ENDIF

        IF(CBBREV(CO,'ERROR',1,noco+1,NTCO,N3CO)) THEN
          ERRMAX=RFROMC(CO(N3CO+1))
        ELSE
C XSL NEWS 11Aug10 default ERRMAX needs to be squared for energy norm
          IF(KTYP007.EQ.4) THEN ! scalar product R.delta
            ERRMAX=1.0D-20
          ELSE ! Other convergence criteria
            ERRMAX=1.0D-10
          ENDIF
C XSL NEWE
        ENDIF

C*** 22/04/08 JHC Removed contact options below. 
C                 Contact parameters now defined through .ipcont file
CC JWF 12/3/2003 Options for contact mechanics.
CC Just checks nr=1. Since if region 1 is a contact problem
CC then whole problem is a contact problem.
C        IF(CBBREV(CO,'CONTACT_STIFFNESS',4,noco+1,NTCO,N3CO)) THEN
C          CALL ASSERT(KTYP5G(1).GE.1,'>>Define contact in ipequa',
C     '    ERROR,*9999)
C          CONT_STIFF=RFROMC(CO(N3CO+1))
C        ELSE
C          CONT_STIFF=1.0d0
C        ENDIF
CC*** 18/03/08 JHC Added frictional contact penalty parameter
C        IF(CBBREV(CO,'TANGENT_STIFFNESS',4,noco+1,NTCO,N3CO)) THEN
C          CALL ASSERT(KTYP5G(1).GE.1,'>>Define contact in ipequa',
C     '    ERROR,*9999)
C          FRIC_STIFF=RFROMC(CO(N3CO+1))
C        ELSE
C          FRIC_STIFF=1.0d0
C        ENDIF
C        IF(CBBREV(CO,'TIED_STIFFNESS',4,noco+1,NTCO,N3CO)) THEN
C          CALL ASSERT(KTYP5G(1).GE.1,'>>Define contact in ipequa',
C     '    ERROR,*9999)
C          TIED_STIFF=RFROMC(CO(N3CO+1))
C        ELSE
C          TIED_STIFF=1.0d0
C        ENDIF
C        IF(CBBREV(CO,'FRICTION_COEFFICIENT',4,noco+1,NTCO,N3CO)) THEN
C          CALL ASSERT(KTYP5G(1).GE.1,'>>Define contact in ipequa',
C     '    ERROR,*9999)
C          FRIC_COEFF=RFROMC(CO(N3CO+1))
C        ELSE
C          FRIC_COEFF=0.1d0
C        ENDIF
C        IF(CBBREV(CO,'AUGMENTATIONS',3,noco+1,NTCO,N3CO)) THEN
C          CALL ASSERT(KTYP5G(1).GE.1,'>>Define contact in ipequa',
C     '    ERROR,*9999)        
C          AUGMENT=IFROMC(CO(N3CO+1))
C        ELSE
C          AUGMENT=0
C        ENDIF
C        IF(CBBREV(CO,'GAP_TOLERANCE',3,noco+1,NTCO,N3CO)) THEN
C          CALL ASSERT(KTYP5G(1).GE.1,'>>Define contact in ipequa',
C     '    ERROR,*9999)
C          GAP_TOL=RFROMC(CO(N3CO+1))
C        ELSE
C          GAP_TOL=0.0d0
C        ENDIF                                          
CC*** 28/02/08 JHC Added option to use penalty method                                           
C        IF(CBBREV(CO,'PENALTY',3,noco+1,NTCO,N3CO)) THEN
C          CALL ASSERT(KTYP5G(1).GE.1,'>>Define contact in ipequa',
C     '    ERROR,*9999)
C          PENALTY=1
C        ELSE
C          PENALTY=0
C        ENDIF
                          
        IF(SOLVE_ACTIV) THEN
          IF(CBBREV(CO,'UPDATE_VM',2,noco+1,NTCO,N3CO)) THEN
            UPDATE_VM=.TRUE.
          ELSE IF(CBBREV(CO,'NOUPDATE_VM',4,noco+1,NTCO,N3CO)) THEN
            UPDATE_VM=.FALSE.
          ELSE
            UPDATE_VM=.TRUE.
          ENDIF
          UPVU=.FALSE.
        ELSE
C!!!    This option is out of date for most/all problems
          IF(CBBREV(CO,'UPDATE',1,noco+1,NTCO,N3CO)) THEN
            UPVU=.TRUE.
            CALL PARSIL(CO(N3CO+1),4,NTIW,IWK,ERROR,*9999)
          ELSE
            UPVU=.FALSE.
          ENDIF
        ENDIF

        IF(SOLVE_ACTIV) THEN
          IF(CBBREV(CO,'BOTH',4,noco+1,NTCO,N3CO)) THEN
            TRANSMEMBRANE=.TRUE.
            EXTRACELLULAR=.TRUE.
          ELSE IF(CBBREV(CO,'NONE',4,noco+1,NTCO,N3CO)) THEN
            TRANSMEMBRANE=.FALSE.
            EXTRACELLULAR=.FALSE.
          ELSE IF(CBBREV(CO,'TRANSMEMBRANE',3,noco+1,NTCO,N3CO)) THEN
            TRANSMEMBRANE=.TRUE.
            EXTRACELLULAR=.FALSE.
          ELSE IF(CBBREV(CO,'EXTRACELLULAR',3,noco+1,NTCO,N3CO)) THEN
            EXTRACELLULAR=.TRUE.
            TRANSMEMBRANE=.FALSE.
          ELSE
            TRANSMEMBRANE=.TRUE.
            EXTRACELLULAR=.TRUE.
          ENDIF

          IF(CBBREV(CO,'STATIC',3,noco+1,NTCO,N3CO)) THEN
            STATIC=.TRUE.
          ELSE
            STATIC=.FALSE.
          ENDIF

          IF(CBBREV(CO,'COUPLED',3,noco+1,NTCO,N3CO)) THEN
            COUPLEDGRIDBEM=.TRUE.
          ELSE IF(CBBREV(CO,'SEPARATE',3,noco+1,NTCO,N3CO)) THEN
            COUPLEDGRIDBEM=.FALSE.
          ELSE
            COUPLEDGRIDBEM=.FALSE.
          ENDIF

          IF(CBBREV(CO,'CPTYPE',5,noco+1,NTCO,N3CO)) THEN
            CPTYPE=IFROMC(CO(N3CO+1))
          ELSE
            CPTYPE=1
          ENDIF

          IF(CBBREV(CO,'FIXED_NODE',5,noco+1,NTCO,N3CO)) THEN
            SOL_ACT_FIX_NODE=IFROMC(CO(N3CO+1))
          ELSE
            SOL_ACT_FIX_NODE=0
          ENDIF
        ENDIF

        IF(CBBREV(CO,'CALC_ACTIV_TIMES',6,noco+1,NTCO,N3CO)) THEN
          CALC_ACTIV_TIMES=.TRUE.
        ELSE
          CALC_ACTIV_TIMES=.FALSE.
        ENDIF

CC CPB 8/6/94 Adding NX_LOC
C        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
C        CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
C     '    ERROR,*9999)

        QUASI=.FALSE.
        DO nonrlist=1,NRLIST(0)
          nr=NRLIST(nonrlist)
          IF(KTYP53(nr).EQ.3) CALL ASSERT(CALL_ACTI, '>>Need to define
     &         active',ERROR,*9999)
          IF(ITYP5(nr,nx).EQ.4) THEN
            QUASI=.TRUE.
          ELSE IF(ITYP5(nr,nx).EQ.5.AND.ITYP2(nr,nx).EQ.3) THEN
            CALL ASSERT(CALL_UPGAUS_EIK,
     '        '>>Need to update Gauss arrays',ERROR,*9999)
          ENDIF
        ENDDO !nonrlist

        IF(QUASI) THEN
          QUASIREGLIST(0)=NRLIST(0)
          DO nonrlist=1,NRLIST(0)
            nr=NRLIST(nonrlist)
            QUASIREGLIST(nonrlist)=NRLIST(nonrlist)
          ENDDO
          COUPLED_QUASI=IS_COUPLED(nx)
        ENDIF

C If the matrices need to be assembled then assemble then
C assemble_global is set in DEEQUA and DEMATE.

        ASSEMBLE=.FALSE.
        DO nonrlist=1,NRLIST(0)
          nr=NRLIST(nonrlist)
          IF(ITYP5(nr,nx).NE.2.AND.ITYP6(nr,nx).NE.2.AND.
     '      KTYP5G(nr).EQ.0) THEN
            !not time dep, or
            !not non-linear, or
            !not a contact problem
            !as assembly is :
            !done in a call from MARCH  routines for time-dep probs
            !done in a call from NONLIN routine for nonlinear probs
            !done in a call from NONLIN_CONT routine for contact probs
            IF(.NOT.ASSEMBLE_GLOBAL(nr,nx)) ASSEMBLE=.TRUE.
          ENDIF
        ENDDO !nr

C       ASSEMBLE=.TRUE.

        IF(SOLVE_ACTIV) ASSEMBLE=.FALSE.

        IF(ASSEMBLE) THEN

C cpb 28/10/98 Adding direct assembly of solution matrices

          nr=NRLIST(1)
          IF((ITYP5(nr,nx).EQ.1.OR.ITYP5(nr,nx).EQ.4).AND.
     '      ITYP6(nr,nx).EQ.1.AND.
     '      (ITYP4(nr,nx).EQ.1.OR.ITYP4(nr,nx).EQ.2)) THEN
C           Problem is static or quasi-static, linear and uses finite
C           elements or boundary elements.
            IF(.NOT.CALC_GLOBAL(nr,nx)) THEN
C***          Need to calculate the sparsity pattern for the global
C***          matrices. This is because you will need to keep the global
C***          matrix entries that correspond to fixed mesh variables in
C***          order to back-subsitute to find the fluxes. However at the
C***          time when the global sparsity patterns are set up (i.e.
C***          in ipequa) you do not know which mesh variables are fixed.
              IF(KTYP24.NE.0) THEN
                WORK_PTR=0
                CALL ALLOCATE_MEMORY(NYT(1,1,nx)*NYT(2,1,nx),1,CHARTYPE,
     '            WORK_PTR,MEM_INIT,ERROR,*9999)
              ENDIF
              CALL CALC_SPARSE_SOLVE(NISC_GKM,NISR_GKM,
     '          %VAL(ISC_GK_PTR(nx)),%VAL(ISR_GK_PTR(nx)),
     '          LGE,NYT(1,1,nx),NYT(2,1,nx),NBH,1,NEELEM,NHE,NPNE,
     '          NPNY(0,1,0,nx),NRLIST,NVHE,nx,NYNE,NYNP,NYNR,NZ_GK_M,
     '          KTYP24,%VAL(WORK_PTR),FIX(1,1,nx),.TRUE.,ERROR,*9999)
              IF(KTYP24.NE.0) CALL FREE_MEMORY(WORK_PTR,ERROR,*9999)
              IF(COUPLED_BEM(nx)) THEN
                IF(KTYP24.NE.0) THEN
                  WORK_PTR=0
                  CALL ALLOCATE_MEMORY(NYT(1,2,nx)*NYT(2,2,nx),1,
     '              CHARTYPE,WORK_PTR,MEM_INIT,ERROR,*9999)
                ENDIF
                CALL CALC_SPARSE_SOLVE(NISC_GQM,NISR_GQM,
     '            %VAL(ISC_GQ_PTR(nx)),%VAL(ISR_GQ_PTR(nx)),
     '            LGE,NYT(1,2,nx),NYT(2,2,nx),NBH,2,NEELEM,NHE,NPNE,
     '            NPNY(0,1,0,nx),NRLIST,NVHE,nx,NYNE,NYNP,NYNR,NZ_GQ_M,
     '            KTYP24,%VAL(WORK_PTR),FIX(1,1,nx),.TRUE.,ERROR,*9999)
                IF(KTYP24.NE.0) CALL FREE_MEMORY(WORK_PTR,ERROR,*9999)
              ENDIF

              IF(FIRSTS(nx)) THEN

                IF(IS_COUPLED(nx)) THEN
                  nrr=0
                ELSE
                  nrr=NRLIST(1)
                ENDIF
                IF(SPARSEGKK(nx).NE.0) THEN
                  WORK_PTR=0
                  CALL ALLOCATE_MEMORY(NOT(1,1,nrr,nx)*NOT(2,1,nrr,nx),
     '              1,CHARTYPE,WORK_PTR,MEM_INIT,ERROR,*9999)
                ENDIF
                CALL CALC_SPARSE_GKK(%VAL(ISC_GK_PTR(nx)),ISC_GKK(1,nx),
     '            %VAL(ISC_GQ_PTR(nx)),%VAL(ISR_GK_PTR(nx)),
     '            ISR_GKK(1,nx),%VAL(ISR_GQ_PTR(nx)),
     '            LGE,NBH,NENP,NHE(1,nx),
     '            NOT(1,1,nrr,nx),NOT(2,1,nrr,nx),NONY(0,1,1,nrr,nx),
     '            NP_INTERFACE,NPNE,NPNY(0,1,0,nx),nrr,NRE,NVHE,nx,
     '            NYNE,NYNP,NYNR(0,0,1,nrr,nx),
     '            %VAL(GK_PTR(nx)),
     '            %VAL(GQ_PTR(nx)),
     '            %VAL(WORK_PTR),IS_COUPLED(nx),.TRUE.,ERROR,*9999)
                IF(SPARSEGKK(nx).NE.0) CALL FREE_MEMORY(WORK_PTR,
     '            ERROR,*9999)

              ENDIF

C$OMP PARALLEL DO
C$OMP&  PRIVATE(nzz),
C$OMP&  SHARED(NZZT,GKK)
              DO nzz=1,NZZT(1,nr,nx)
                GKK(nzz,nx)=0.0d0
              ENDDO !nzz
C$OMP END PARALLEL DO

            ENDIF
          ENDIF

          IF(IS_COUPLED(nx)) THEN
            CALL INIT_SPARSE_MATRIX(NYNR(0,2,1,0,nx),NYNR(0,2,1,0,nx),
     '        %VAL(ISC_GK_PTR(nx)),%VAL(ISR_GK_PTR(nx)),
     '        0,0,0,NYT(1,1,nx),NZ_GK_M,NZT(1,nx),
     '        NYNR(0,1,1,0,nx),NYNR(0,1,1,0,nx),KTYP24,
     '        %VAL(GK_PTR(nx)),ERROR,*9999)
            IF(COUPLED_BEM(nx)) THEN
              CALL INIT_SPARSE_MATRIX(NYNR(0,2,2,0,nx),
     '          NYNR(0,2,2,0,nx),%VAL(ISC_GQ_PTR(nx)),
     '          %VAL(ISR_GQ_PTR(nx)),0,0,0,NYT(1,2,nx),
     '          NZ_GQ_M,NZT(2,nx),NYNR(0,1,2,0,nx),NYNR(0,1,2,0,nx),
     '          KTYP24,%VAL(GQ_PTR(nx)),ERROR,*9999)
            ENDIF
          ENDIF

          DO nonrlist=1,NRLIST(0)
            nr=NRLIST(nonrlist)
            IF(IS_COUPLED(nx)) THEN
              nr_gkk=0
            ELSE
              nr_gkk=nr
            ENDIF
            CALL GENASSEM(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,IBT,IDO,
     '        INP,ISC_GD,%VAL(ISC_GK_PTR(nx)),ISC_GKK(1,nx),ISC_GM,
     '        %VAL(ISC_GQ_PTR(nx)),ISR_GD,%VAL(ISR_GK_PTR(nx)),
     '        ISR_GKK(1,nx),ISR_GM,%VAL(ISR_GQ_PTR(nx)),LGE,NBH,NBJ,
     '        NDET,NDIPOLES,NEELEM,NENP,NGAP,NHE(1,nx),NHP(1,nr,nx),
     '        NKH(1,1,1,nr),NKHE,NKJE,NLL,NONY(0,1,1,nr_gkk,nx),NORD,
     '        NPF,NP_INTERFACE,NPB,NPNE,NPNODE,NPNY(0,1,0,nx),nr,nr_gkk,
     '        NRE,NRLIST,NVHE,NVHP(1,1,1,nr),NVJE,NW(1,1,nx),nx,NYNE,
     '        NYNP,NYNR(0,0,1,nr,nx),CE(1,1,nx),CG,CGE(1,1,1,nx),
     '        CONY(0,1,1,nr_gkk,nx),CP(1,1,nx),CURVCORRECT,DET,
     '        DIPOLE_CEN,DIPOLE_DIR,DL,DRDN,DRDNO,ED,EM,ER,ES,GD,
     '        %VAL(GK_PTR(nx)),GKK(1,nx),GM,%VAL(GQ_PTR(nx)),GR,PG,
     '        RAD,RD,RG,SE,WG,XA,XE,XG,XG1,XIG,XN,XN_GRAD,XP,XR,XR_GRAD,
     '        YG,ZA,ZE,ZG,ZP,ERROR,*9999)
          ENDDO !nonrlist

          IF(IS_COUPLED(nx)) THEN !write out global(nr=0) stiff matrices
            IF(IWRIT4(nr,nx).GE.4) THEN
              WRITE(OP_STRING,'(/'' Global region (nr=0):'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              IF(COUPLED_BEM(nx)) THEN
                WRITE(OP_STRING,'(/'' Global stiffness matrix GK:'')')
              ELSE
                WRITE(OP_STRING,
     '           '(/'' Global load vector GR & stiffness matrix GK:'')')
              ENDIF
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' NYNR(0,1,1,0)='',I5,'
     '          //''', NYNR(0,2,1,0)='',I5)') NYNR(0,1,1,0,nx),
     '          NYNR(0,2,1,0,nx)
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              IF(COUPLED_BEM(nx)) THEN
                CALL OPSTFMAT(NYNR(0,2,1,0,nx),%VAL(ISC_GK_PTR(nx)),
     '            %VAL(ISR_GK_PTR(nx)),IOOP,
     '            NYT(1,1,nx),NYT(2,1,nx),NZT(1,nx),NYNR(0,1,1,0,nx),
     '            KTYP24,
     '            %VAL(GK_PTR(nx)),
     '            GR,'GK ','GR ',.FALSE.,.TRUE.,.FALSE.,ERROR,*9999)
                WRITE(OP_STRING,'(/'' Global stiffness matrix GQ:'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' NYNR(0,1,2,0)='',I5,'
     '            //''', NYNR(0,2,1,0)='',I5)') NYNR(0,1,2,0,nx),
     '            NYNR(0,2,2,0,nx)
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                CALL OPSTFMAT(NYNR(0,2,2,0,nx),%VAL(ISC_GQ_PTR(nx)),
     '            %VAL(ISR_GQ_PTR(nx)),IOOP,
     '            NYT(1,2,nx),NYT(2,2,nx),NZT(2,nx),NYNR(0,1,2,0,nx),
     '            KTYP24,%VAL(GQ_PTR(nx)),GR,'GQ ','GR ',.FALSE.,
     '            .TRUE.,.FALSE.,ERROR,*9999)
              ELSE
                CALL OPSTFMAT(NYNR(0,2,1,0,nx),%VAL(ISC_GK_PTR(nx)),
     '            %VAL(ISR_GK_PTR(nx)),IOOP,
     '            NYT(1,1,nx),NYT(2,1,nx),NZT(1,nx),NYNR(0,1,1,0,nx),
     '            KTYP24,%VAL(GK_PTR(nx)),
     '            GR,'GK ','GR ',.FALSE.,.TRUE.,.TRUE.,ERROR,*9999)
              ENDIF
            ENDIF !iwrit4
          ENDIF !is_coupled
        ENDIF !assemble

        nr=NRLIST(1)
C SGM 26Oct2000 grid-based Finite element also
C MLT 29Nov02 Grid Finite Volume also
        IF((ITYP4(nr,nx).EQ.4.OR.ITYP4(nr,nx).EQ.6.OR.
     '      ITYP4(nr,nx).EQ.7).AND.
     '    (ITYP2(nr,nx).EQ.3).AND.
     '    (ITYP6(nr,nx).EQ.1).AND.(ITYP5(nr,nx).EQ.1).AND.
     '    SOLVE_ACTIV) THEN
          IF(CPTYPE.LE.3) THEN
            nx1=NXLIST(nx+1)
            nx2=nx
            nx3=nx
          ELSE IF(CPTYPE.EQ.4) THEN
            nx1=NXLIST(nx+1)
            nx2=NXLIST(nx+2)
            nx3=nx
          ELSE IF(CPTYPE.EQ.6) THEN
            nx1=NXLIST(2)
            nx2=NXLIST(3)
            nx3=nx
          ELSE IF(CPTYPE.EQ.7) THEN
            nx1=NXLIST(2)
            nx2=NXLIST(3)
            nx3=NXLIST(4)
          ENDIF
C SGM 26Oct2000 grid-based Finite element also
C MLT 29Nov02 Grid Finite Volume also
        ELSE IF((ITYP19(nr,nx).EQ.1).AND.
     '      (ITYP4(nr,nx).EQ.4.OR.ITYP4(nr,nx).EQ.6.OR.
     '       ITYP4(nr,nx).EQ.7).AND.
     '      (ITYP2(nr,nx).EQ.9).AND.(ITYP5(nr,nx).EQ.2).AND.
     '      SOLVE_ACTIV) THEN
          IF(CPTYPE.EQ.1) THEN
            nx1=nx
            nx2=nx
            nx3=nx
          ELSE IF(CPTYPE.EQ.2) THEN
            nx1=NXLIST(NXLIST(0))
            nx2=nx
            nx3=nx
          ELSE IF(CPTYPE.EQ.3) THEN
            nx1=NXLIST(NXLIST(0)-1)
            nx2=nx
            nx3=nx
          ELSE IF(CPTYPE.EQ.4) THEN
            nx1=NXLIST(NXLIST(0)-1)
            nx2=NXLIST(NXLIST(0))
            nx3=nx
          ELSE IF(CPTYPE.EQ.5) THEN
            nx1=NXLIST(NXLIST(0)-2)
            nx2=NXLIST(NXLIST(0)-1)
            nx3=nx
          ELSE IF(CPTYPE.EQ.6) THEN
            nx1=NXLIST(3)
            nx2=NXLIST(4)
            nx3=nx
          ELSE IF(CPTYPE.EQ.7) THEN
            nx1=NXLIST(3)
            nx2=NXLIST(4)
            nx3=NXLIST(5)
          ENDIF
        ELSE
          nx1=nx
          nx2=nx
          nx3=nx
        ENDIF

        CALL GENSOL(CELL_ICQS_VALUE,DIPOLE_CEN_NTIME,
     &    DIPOLE_DIR_NTIME,IBT,ICQS,ICQS_SPATIAL,IDO,IICQS_SPATIAL,
     '    IRCQS_SPATIAL,INP,ISC_GD,%VAL(ISC_GK_PTR(nx)),ISC_GKK,
     '    ISC_GM,ISC_GMM,%VAL(ISC_GQ_PTR(nx)),ISR_GD,
     '    %VAL(ISR_GK_PTR(nx)),ISR_GKK,ISR_GM,ISR_GMM,
     '    %VAL(ISR_GQ_PTR(nx)),ISEG,ISELNO,ISFIBR,ISFIEL,
     '    ISLINE,ISLINO,ISNONO,ISTATE,IUSER,IWK,LGE,MXI,NAN,NAQ,
     '    NBH,NBHF,NBJ,NBJF,NDET,NDIPOLES,NEELEM,NELIST,NENFVC,
     '    NENP,NENQ,NLATNE,NLATNQ,NLATPNQ,NQNLAT,NEP,NFF,NFFACE,NFVC,
     &    NGAP,NHE,NHP,NHQ,NKB,NKEF,NKH,NKHE,NKJE,NLL,NLLIST,NLQ,NMNO,
     &    NNB,NNF,NODENVC,NODENVCB,NONY,NORD,NPB,NP_INTERFACE,NPF,NPL,
     &    NPLIST,NPLIST2,NPNE,NPNODE,NPNY,NQET,NQGP,NQGP_PIVOT,NQGW,
     &    NQLIST,NQNE,NQNP,NQNY,NQS,NQSCNB,NQXI,NRLIST,NRLIST2,NRE,NSB,
     &    NTIME_INTERP,NTIME_POINTS,NTIME_NR,NTIW,NVCB,NVCNODE,NVHE,
     &    NVHP,NVJE,NVJP,NW(1,1,nx),NWQ,nx,NXI,NXLIST,NXQ,NYNE,NYNO,
     &    NYNP,NYNQ,NYNR,NYQNR,TV_BC_SET,Z_CONT_LIST,ACINUS,ALPHA,AQ,
     &    BBM,CE,CELL_RCQS_VALUE,CG,CGE,CONY,CP,CQ,CURVCORRECT,CYNO,
     &    DET,DIPOLE_CEN,DIPOLE_DIR,DL,DNUDXQ,DRDN,DRDNO,DXDXIQ,DXDXIQ2,
     &    ED,EIGVAL,EIGVEC,EM,ER,ERRMAX,ES,FEXT,GCHQ,GD,
     &    %VAL(GK_PTR(nx1)),%VAL(GK_PTR(nx2)),%VAL(GK_PTR(nx3)),GKK,GM,
     &    GMM,%VAL(GQ_PTR(nx1)),%VAL(GQ_PTR(nx2)),%VAL(GQ_PTR(nx3)),GR,
     &    GRR,GUQ,PAOPTI,PG,PMIN,PMAX,PROPQ,R,RAD,RCQS,RCQS_SPATIAL,RD,
     &    RE1,RESID,RESJAC,RG,RHS,SE,TIME_VALUES,USER,VC,
     &    VC_INIT,WARM_START,WG,XA,XAB,XC,XE,XG,XG1,XIG,XIP,XIQ,XN,XNFV,
     &    XN_GRAD,XO,XP,XQ,XR,XR_GRAD,YG,YGF,YP,YQ,YQS,ZA,ZA1,Z_CONT,ZE,
     &    ZE1,ZG,ZNFV,ZP,ZP1,CSEG,FIX,FIXQ,ITER8,
     &     ERROR,*9999)

        IF(OPFILE) THEN
          CALL CLOSEF(IO4,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('SOLVE')
      RETURN
 9999 CALL ERRORS('SOLVE',ERROR)
      IF(OPFILE) CLOSE(UNIT=IO4)
      CALL EXITS('SOLVE')
      RETURN 1
      END


