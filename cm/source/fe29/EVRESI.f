      SUBROUTINE EVRESI(IBT,IDO,INP,ISIZE_MFI,
     '  ISIZE_PHI,ISIZE_TBH,LD,LDR,LGE,
     '  NAN,NBH,NBHF,NBJ,NBJF,NEELEM,NEL,NELIST,NENP,
     '  NFF,NFFACE,NGAP,NHE,NHP,
     '  NKB,NKEF,NKH,NKHE,NKJE,NLL,NLNO,NMNO,NNB,NNF,NNL,NONL,NONM,
     '  NONY,NP1OPT,NPF,NP_INTERFACE,NPL,NPLIST3,NPNE,
     '  NPNODE,NPNY,NRE,NRLIST,NSB,NVHE,NVHP,NVJE,NVJL,
     '  NW,NXI,NXLIST,NYNE,NYNO,NYNP,NYNR,PAOPTY,Z_CONT_LIST,
     '  AQ,CE,CELL_RCQS_VALUE,CG,CGE,CONY,CP,CURVCORRECT,
     '  DL,D_RE,D_RI3,D_RP,D_TG,D_ZG,ES,FEXT,FGRAD,LAPL,LAPLSQR,MFI,
     '  PAOPTI,PBOPTI,PG,PHI,PHI_H,PMIN,PMAX,RE1,RE2,RESID,
     '  RESJAC,RG,SE,T_BH,WG,WK1_INV,WU,XA,XE,XG,XID,
     '  XIG,XN,XP,
     '  YG,YGF,YP,ZA,ZA1,Z_CONT,ZD,ZE,ZG,ZG1,ZP,ZP1,STRING,
     '  FIX,ERROR,*)

C#### Subroutine: EVRESI
C###  Description:
C###    <HTML>
C###    EVRESI evaluates residuals.<BR>
C###    WITH is used to eval resids with a specified set of parameters.
C###    <BR>
C###    PARAMTYPE is the type of parameters being used:
C###    <UL>
C###    <LI> 'GEOMETRIC_PARAMETERS' ->eval finite elast. resids for solve
C###    <LI> 'MATERIAL_PARAMETERS'  ->eval for mat. param. optim. resids
C###    <LI> 'DATA_FITTING'         ->eval data fitting resids
C###    <LI> 'TORSO_CUSTOMISE'      ->eval torso customise residuals
C###    <LI> 'MAGNETIC'             ->eval magnetic field
C###    <LI> 'MOMENTS'              ->eval second moments
C###    <LI> 'POTENTIAL'            ->eval potential residuals
C###    <LI> 'SHEAR'                ->eval shear strain residuals
C###    <LI> 'STRESS'               ->eval stress residuals
C###    <LI> 'STRAIN'               ->eval strain residuals
C###    <LI> 'PRESSURE'             ->eval strain from pressure relationship
C###    </UL>
C###    ANALYT_DERIVS/NUM_DERIVS evals the derivs of the resids wrt the
C###      PARAMTYPE params analytically/using numerical differencing
C###      (not implemented for DATA_FITTING)
C###    <BR>
C###    STEP is used for incrementing a material parameter
C###    <BR>
C###    COUPLED is used to evaluate residuals for coupled region (nr=0)
C###    <BR>
C###    VIEW/NOVIEW does/doesn't print residuals to screen
C###    </HTML>

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'coup00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'opti00.cmn'
      INCLUDE 'tol00.cmn'
      INCLUDE 'trsf00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  ISIZE_MFI(3,NSSM),INP(NNM,NIM,NBFM),
     '  ISIZE_PHI(2),ISIZE_TBH(2),LD(NDM),LDR(0:NDM),LGE(NHM*NSM,NRCM),
     '  NAN(NIM,NAM,NBFM),
     '  NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NEELEM(0:NE_R_M,0:NRM),NEL(0:NELM,NLM),NELIST(0:NEM),
     '  NENP(NPM,0:NEPM,0:NRM),
     '  NFF(6,NEM),NFFACE(0:NF_R_M,NRM),NGAP(NIM,NBM),
     '  NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),
     '  NKB(2,2,2,NNM,NBFM),NKEF(0:4,16,6,NBFM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),
     '  NLNO(NOPM,NXM),NMNO(1:2,0:NOPM,NXM),
     '  NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),
     '  NONL(NLM,NXM),NONM(NMM,NPM,NXM),NONY(0:NOYM,NYM,NRCM,0:NRM,NXM),
     '  NP_INTERFACE(0:NPM,0:3),NP1OPT(NOPM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NPL(5,0:3,NLM),NPLIST3(0:NPM),
     '  NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),NRE(NEM),NSB(NKM,NNM,NBFM),
     '  NRLIST(0:NRM),NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJL(4,NJM,NLM),
     '  NW(NEM,3,NXM),NXI(-NIM:NIM,0:NEIM,0:NEM),NXLIST(0:NXM),
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
     '  FEXT(NIFEXTM,NGM,NEM),FGRAD(NOPM),
     '  LAPL(NY_TRANSFER_M,NY_TRANSFER_M),
     '  LAPLSQR(NY_TRANSFER_M,NY_TRANSFER_M),
     '  MFI(NDM,NTSM,3,NSSM),
     '  PAOPTI(NOPM),PBOPTI(NOPM),PG(NSM,NUM,NGM,NBM),
     '  PHI(NY_TRANSFER_M,NTSM),PHI_H(NY_TRANSFER_M,NTSM),
     '  PMIN(NOPM+NCOM),PMAX(NOPM+NCOM),RE1(NSM,NHM),RE2(NSM,NHM),
     '  RESID(NREM),RESJAC(NREM,NOPM),RG(NGM),SE(NSM,NBFM,NEM),
     '  T_BH(NY_TRANSFER_M,NY_TRANSFER_M),
     '  WG(NGM,NBM),WK1_INV(NY_TRANSFER_M,NY_TRANSFER_M),
     '  WU(0:NUM+1,NEM),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XG(NJM,NUM),XID(NIM,NDM),XIG(NIM,NGM,NBM),
     '  XN(NJM,NGM),XP(NKM,NVM,NJM,NPM),
     '  YG(NIYGM,NGM,NEM),YGF(NIYGFM,NGFM,NFM),YP(NYM,NIYM,NXM),
     '  ZA(NAM,NHM,NCM,NEM),ZA1(NAM,NHM,NCM,NEM),Z_CONT(NDM,2,67),
     '  ZD(NJM,NDM),ZE(NSM,NHM),ZG(NHM,NUM),ZG1(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM),ZP1(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM)

!     Local Variables
      INTEGER ERR,i,IBEG,IBEG1,IBEG2,IEND,IEND1,
     '  IEND2,KTYPTEMP,LIST_RESID,
     '  MODE,n1step,n2step,N3CO,nd0,nd1,
     '  ne,NGLIST(0:NGM),nm,nm1,nm2,nmlist,NMMX_LOCAL,
     '  no,no1,no2,noelem,no_nglist,
     '  no_nrlist,no_nynr,no_nynr1,no_nynr2,
     '  nonode,noopti,nopara,nopara1,nopara2,nores,nostep,
     '  noy,noy1,noy2,np,nr,nr_coup,NTLIST,NTSTEP,
     '  nx,nxc,nx_opt,nx_sol,ny,ny1,ny2
      PARAMETER(NMMX_LOCAL=30)      
      REAL*8 CELL_RCQS_VALUE_LOCAL(NMMX_LOCAL),DELTA,P1STEP,P2STEP,
     '  SUMSQR,SUMSQreactions,ZE1(NSM,NHM)
      CHARACTER COORDS*9,CHAR1*144,CHAR2*11,CHAR3*3,
     '  PARAMTYPE*20,FILE*(MXCH)
      LOGICAL ABBREV,ALL_REGIONS,ANALYTIC,CBBREV,COUPLED_RES,
     '  DEFORMED,DERIV,OPFILE,PARAMS_CHANGED,VIEW
      

! Functions
      INTEGER CALC_SAMPLE_FROM_TIME,IFROMC
      REAL*8 CALC_TIME_FROM_SAMPLE,RFROMC

      CALL ENTERS('EVRESI',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        WRITE(CHAR1,'(E11.4)') PAOPTI(1)
        nxc=1!temporary
        CALL NX_LOC(NX_INQUIRE,nxc,nx_opt,NX_OPTI,ERROR,*9999)
        IF(nx_opt.GT.0) THEN
          DO noopti=2,NMNO(1,0,nx_opt)
            CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
            WRITE(CHAR2,'(E11.4)') PAOPTI(noopti)
            CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
            CHAR1=CHAR1(IBEG1:IEND1)//','//CHAR2(IBEG2:IEND2)
          ENDDO
        ELSE
          CHAR1=' '
        ENDIF
        CALL STRING_TRIM(CHAR1,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM evaluate residuals<;FILENAME> wrt mat_params
C###  Parameter:        <with PARAM_VALUES[]>
C###    Specify the parameters values
C###  Parameter:        <step NO_STEPS#[0]>
C###    Increment a material parameter
C###  Parameter:        <(view/noview)[view]>
C###    Determine whether to print out the residuals
C###  Parameter:        analyt_derivs
C###  Parameter:          <paramnum PARAM# (omit for all)[all]>
C###    Evalates the derivatives of the residuals wrt the
C###    PARAMTYPE parameters analytically.
C###  Parameter:        num_derivs
C###  Parameter:          <paramnum PARAM# (omit for all)[all]>
C###  Parameter:          <interval SIZE[1.E-6]>
C###    Evalates the derivatives of the residuals wrt the
C###    PARAMTYPE parameters using numerical differencing. Can
C###    also specify the size of interval for numerical differencing.
C###  Parameter:        <region (#s/all)[1]>
C###    Specify the region numbers to evaluate.
C###  Parameter:        <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Evaluate residuals with respect to material
C###      parameters optimisation residuals.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> wrt mat_params'
        OP_STRING(2)=BLANK(1:15)
     '    //'<with PARAM_VALUES['//CHAR1(IBEG1:IEND1)//']>'
        OP_STRING(3)=BLANK(1:15)//'<step NO_STEPS[0]>'
        OP_STRING(4)=BLANK(1:15)//'<(view/noview)[view]>'
        OP_STRING(5)=BLANK(1:15)//'analyt_derivs'
        OP_STRING(6)=BLANK(1:15)
     '    //'    <paramnum PARAM# (omit for all)[all]>'
        OP_STRING(7)=BLANK(1:15)//'num_derivs'
        OP_STRING(8)=BLANK(1:15)
     '    //'    <paramnum PARAM# (omit for all)[all]>'
        OP_STRING(9)=BLANK(1:15)//'    <interval SIZE[1.E-6]>'
        OP_STRING(10)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(11)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM evaluate residuals<;FILENAME> wrt geom_params
C###  Parameter:        <coupled>
C###    Used to evaluate residuals for coupled problem
C###  Parameter:        <(view/noview)[view]>
C###    Determine whether to print out the residuals
C###  Parameter:        analyt_derivs
C###  Parameter:          <paramnum PARAM# (omit for all)[all]>
C###    Evalates the derivatives of the residuals wrt the
C###    PARAMTYPE parameters analytically.
C###  Parameter:        num_derivs
C###  Parameter:          <paramnum PARAM# (omit for all)[all]>
C###  Parameter:          <interval SIZE[1.E-6]>
C###    Evalates the derivatives of the residuals wrt the
C###    PARAMTYPE parameters using numerical differencing.
C###    Specifying the interval for differencing.
C###  Parameter:        <region (#s/all)[1]>
C###    Specify the region numbers to evaluate.
C###  Parameter:        <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Evaluate residuals wrt to geometric parameters.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> wrt geom_params'
        OP_STRING(2)=BLANK(1:15)//'<coupled>'
C###    Used to evaluate residuals for coupled problem
        OP_STRING(3)=BLANK(1:15)//'<(view/noview)[view]>'
C###    Determine whether to print out the residuals
        OP_STRING(4)=BLANK(1:15)//'analyt_derivs'
        OP_STRING(5)=BLANK(1:15)
     '    //'    <paramnum PARAM# (omit for all)[all]>'
        OP_STRING(6)=BLANK(1:15)//'num_derivs'
        OP_STRING(7)=BLANK(1:15)
     '    //'    <paramnum PARAM# (omit for all)[all]>'
        OP_STRING(8)=BLANK(1:15)//'    <interval SIZE[1.E-6]>'
        OP_STRING(9)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(10)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM evaluate residuals<;FILENAME> wrt data_fitting
C###  Parameter:        <(function/jacobian/both)[both]>
C###    Evaluate the function/jacobian or both
C###  Parameter:        <region (#/all)[1]>
C###    Specify the region numbers to evaluate.
C###  Parameter:        <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Parameter:        <data GROUP#>
C###    Specify which data points to work on.
C###  Parameter:        <undeformed/deformed [undeformed]>
C###    Specify the nodal values to use.
C###  Parameter:        <(view/noview)[noview]>
C###    Determine whether to print out the residuals.
C###  Description:
C###    Evaluate the residuals for data fitting problems

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> wrt data_fitting'
        OP_STRING(2)=BLANK(1:15)//'<(function/jacobian/both)[both]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#/all)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(5)=BLANK(1:15)//'<data GROUP#>'
        OP_STRING(6)=BLANK(1:15)//'<undeformed/deformed [undeformed]>'
        OP_STRING(7)=BLANK(1:15)//'<(view/noview)[noview]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM evaluate residuals<;FILENAME>

C###  Parameter:        <wrt (torso_customise/moments)[]>
C###    Evaluates residuals wrt:
C###    TORSO_CUSTOMISE: Refer to Masters Thesis Joanne Crocombe.
C###    MOMENTS:  wrt to second moments
C###  Parameter:        <region (#/all)[1]>
C###    Specify the region numbers to evaluate.
C###  Parameter:        <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Evaluate residuals for a torso customisation.
C###

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<region (#/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
C#### Command: FEM evaluate residuals<;FILENAME> wrt magnetic
C###  Parameter:        <tstart (#/begining)[begining]>
C###    Specify start time
C###  Parameter:        <tend (#/end)[end]>
C###    Specify end time
C###  Parameter:        <list_residuals #>[0]
C###    Outputs additional residual information. The default value
C###    of 0 does not output any information. The higher the value
C###    the higher the output level.
C###  Description:
C###    Evaluates the residuals in magnetic fields

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> wrt magnetic'
        OP_STRING(2)=BLANK(1:15)//'<tstart (#/begining)[begining]>'
        OP_STRING(3)=BLANK(1:15)//'<tend (#/end)[end]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#/all)[1]>'
        OP_STRING(5)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(6)=BLANK(1:15)//'<list_residuals #>[0]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
C#### Command: FEM evaluate residuals<;FILENAME> wrt potential
C###  Parameter:        <(function/jacobian/both)[both]>
C###    Evaluate the function/jacobian or both
C###  Parameter:        <tstart (#/begining)[begining]>
C###    Specify start time
C###  Parameter:        <tend (#/end)[end]>
C###    Specify end time
C###  Parameter:        <class #>[1]
C###    Specify the class number (of solve type) to use.
C###  Parameter:        <list_residuals #>[0]
C###    Outputs additional residual information. The default value
C###    of 0 does not output any information. The higher the value
C###    the higher the output level.
C###  Description:
C###    Evaluates the differences between potentials in the PHI and
C###    PHI_H arrays at each time

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> wrt potential'
        OP_STRING(2)=BLANK(1:15)//'<(function/jacobian/both)[both]>'
        OP_STRING(3)=BLANK(1:15)//'<tstart (#/begining)[begining]>'
        OP_STRING(4)=BLANK(1:15)//'<tend (#/end)[end]>'
        OP_STRING(5)=BLANK(1:15)//'<class #>[1]'
        OP_STRING(6)=BLANK(1:15)//'<list_residuals #>[0]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
C#### Command: FEM evaluate residuals<;FILENAME> wrt stress
C###  Parameter:        <(fibre/wall/reference)[fibre]>
C###    Evaluate the residuals using fibre stress or strain
C###  Parameter:        <class #>[1]
C###    Specify the class number (of solve type) to use.
C###  Parameter:        <list_residuals #>[0]
C###    Outputs additional residual information. The default value
C###    of 0 does not output any information. The higher the value
C###    the higher the output level.
C###  Description:
C###    Evaluates stress residuals

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> wrt stress'
        OP_STRING(2)=BLANK(1:15)//'<(fibre/wall/reference)[fibre]>'
        OP_STRING(3)=BLANK(1:15)//'<class #>[1]'
        OP_STRING(4)=BLANK(1:15)//'<list_residuals #>[0]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
C#### Command: FEM evaluate residuals<;FILENAME> wrt strain
C###  Parameter:        <(fibre/wall/reference)[fibre]>
C###    Evaluate the residuals using fibre stress or strain
C###  Parameter:        <class #>[1]
C###    Specify the class number (of solve type) to use.
C###  Parameter:        <list_residuals #>[0]
C###    Outputs additional residual information. The default value
C###    of 0 does not output any information. The higher the value
C###    the higher the output level.
C###  Description:
C###    Evaluates strain residuals

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> wrt strain'
        OP_STRING(2)=BLANK(1:15)//'<(fibre/wall/reference)[fibre]>'
        OP_STRING(3)=BLANK(1:15)//'<class #>[1]'
        OP_STRING(4)=BLANK(1:15)//'<list_residuals #>[0]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
C#### Command: FEM evaluate residuals<;FILENAME> wrt pressure
C###  Parameter:        <(fibre/wall/reference)[fibre]>
C###    Evaluate the residuals using pressure strain relationship
C###  Parameter:        <class #>[1]
C###    Specify the class number (of solve type) to use.
C###  Parameter:        <list_residuals #>[0]
C###    Outputs additional residual information. The default value
C###    of 0 does not output any information. The higher the value
C###    the higher the output level.
C###  Description:
C###    Evaluates pressure residuals

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> wrt pressure'
        OP_STRING(2)=BLANK(1:15)//'<(fibre/wall/reference)[fibre]>'
        OP_STRING(3)=BLANK(1:15)//'<class #>[1]'
        OP_STRING(4)=BLANK(1:15)//'<list_residuals #>[0]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
C#### Command: FEM evaluate residuals<;FILENAME> wrt reaction
C###    Evaluate the residuals using reaction force data
C###  Parameter:        <class #>[1]
C###    Specify the class number (of solve type) to use.
C###  Parameter:        <list_residuals #>[0]
C###    Outputs additional residual information. The default value
C###    of 0 does not output any information. The higher the value
C###    the higher the output level.
C###  Description:
C###    Evaluates pressure residuals

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> wrt reaction'
        OP_STRING(3)=BLANK(1:15)//'<class #>[1]'
        OP_STRING(4)=BLANK(1:15)//'<list_residuals #>[0]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
C#### Command: FEM evaluate residuals data_error
C###  Parameter:        <data field#[1,2,3]>
C###    The original data
C###  Parameter:        <model field#[4,5,6]>
C###    The data fields which hold the data points
C###    calculated in the model

        OP_STRING(1)=STRING(1:IEND)//' data_error'
        OP_STRING(2)=BLANK(1:15)//'<data field#[1,2,3]>'
        OP_STRING(3)=BLANK(1:15)//'<model field#[4,5,6]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','EVRESI',ERROR,*9999)
      ELSE

C!!! LKC 18-DEC-1999 This stuff could be speed up a great deal
C!!!   especially as it is continuously called by the optimiser.
C!!!   Need to somehow by pass all this parsing after the first
C!!!   time?
C!!!   A quick way to overcome this is to move all the heavily
C!!!   used calls to the top of the if blocks.


        DEFORMED=.FALSE.
        IF(NTCOQU(noco).GT.0) THEN
          FILE=COQU(noco,1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          OPFILE=.TRUE.
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opresi','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
          IOFI=IOOP
        ENDIF

        CALL ASSERT(USE_NPSOL.GT.0,'>>Set USE_NPSOL to 1 '
     '    //'in parameters file to evaluate residuals',ERROR,*9999)

        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)

C       rgb 25-03-2000 Initialising string so doesn't crop up in ftncheck
        COORDS='INIT'

        IF(CBBREV(CO,'COUPLED',2,noco+1,NTCO,N3CO)) THEN
          CALL NX_LOC(NX_INQUIRE,nxc,nx_sol,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx_sol.GT.0,
     '      '>>No nx defined for this solve class',ERROR,*9999)
          CALL ASSERT(IS_COUPLED(nx_sol),
     '      '>>Define solve for coupled problem',ERROR,*9999)
          COUPLED_RES=.TRUE.
          NRLIST(0)=COUP_NRLIST(0,nx_sol)
          DO no_nrlist=1,COUP_NRLIST(0,nx_sol)
            NRLIST(no_nrlist)=COUP_NRLIST(no_nrlist,nx_sol)
          ENDDO
        ELSE
          COUPLED_RES=.FALSE.
          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,
     '      ERROR,*9999)
        ENDIF
C        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
C        nxc=NXLIST(1)

        IF(CBBREV(CO,'WRT',2,noco+1,NTCO,N3CO)) THEN
          IF(ABBREV(CO(N3CO+1),'GEOM_PARAMS',1)) THEN
            PARAMTYPE='GEOMETRIC_PARAMETERS'
          ELSE IF(ABBREV(CO(N3CO+1),'MAT_PARAMS',1)) THEN
            PARAMTYPE='MATERIAL_PARAMETERS'
          ELSE IF(ABBREV(CO(N3CO+1),'SHEAR',1)) THEN
            PARAMTYPE='SHEAR'
          ELSE IF(ABBREV(CO(N3CO+1),'DATA_FITTING',1)) THEN
            PARAMTYPE='DATA_FITTING'
          ELSE IF(ABBREV(CO(N3CO+1),'TORSO_CUSTOMISE',1)) THEN
            PARAMTYPE='TORSO_CUSTOMISE'
          ELSE IF(ABBREV(CO(N3CO+1),'MOMENTS',1)) THEN
            PARAMTYPE='MOMENTS'
          ELSE IF(ABBREV(CO(N3CO+1),'MAGNETIC',1)) THEN
            PARAMTYPE='MAGNETIC'
          ELSE IF(ABBREV(CO(N3CO+1),'POTENTIAL',1)) THEN
            PARAMTYPE='POTENTIAL'
          ELSE IF(ABBREV(CO(N3CO+1),'STRESS',1)) THEN
            PARAMTYPE='STRESS'
          ELSE IF(ABBREV(CO(N3CO+1),'STRAIN',1)) THEN
            PARAMTYPE='STRAIN'
          ELSE IF(ABBREV(CO(N3CO+1),'PRESSURE',1)) THEN
            PARAMTYPE='PRESSURE'
          ELSE IF(ABBREV(CO(N3CO+1),'REACTION',1)) THEN
            PARAMTYPE='REACTION'
          ELSE
            ERROR='Unknown parameter'
            GOTO 9999
          ENDIF
        ELSE
          PARAMTYPE='MATERIAL_PARAMETERS'
        ENDIF

        IF(PARAMTYPE(1:19).EQ.'MATERIAL_PARAMETERS'.OR.
     '    PARAMTYPE(1:20).EQ.'GEOMETRIC_PARAMETERS'.OR.
     '    PARAMTYPE(1:5).EQ.'SHEAR'.OR.
     '    PARAMTYPE(1:6).EQ.'STRESS'.OR.
     '    PARAMTYPE(1:6).EQ.'STRAIN'.OR.
     '    PARAMTYPE(1:8).EQ.'REACTION'.OR.
     '    PARAMTYPE(1:8).EQ.'PRESSURE') THEN
          CALL ASSERT(USE_NONLIN.GT.0,'>>Set USE_NONLIN to 1 to '
     '      //'evaluate residuals',ERROR,*9999)
        ENDIF

        MODE=2 !according to EO4UPF
        IF(PARAMTYPE(1:19).EQ.'MATERIAL_PARAMETERS') THEN

C cpb 19/3/95 To be completely general the nx classes for solve and
C and optimisation should be different.

          CALL ASSERT(.NOT.COUPLED_RES,'>>Not implemented for coupled '
     '      //'problems',ERROR,*9999)
          CALL NX_LOC(NX_INQUIRE,nxc,nx_opt,NX_OPTI,ERROR,*9999)
          CALL ASSERT(nx_opt.GT.0,
     '      '>>No nx defined for this optimisation class',
     '      ERROR,*9999)
          CALL NX_LOC(NX_INQUIRE,nxc,nx_sol,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx_sol.GT.0,
     '      '>>No nx defined for this solve class',
     '      ERROR,*9999)
C HS 26/4/04 moved here from below
          CALL ASSERT(NMNO(1,0,nx_opt).LE.NOPM,'>>NOPM too small',
     '      ERROR,*9999)
          CALL ASSERT(NTOPTI.LE.NOPM,'>>NOPM too small',
     '      ERROR,*9999)
C news HS 18-April-2004:
C added this assert for optimisation of material parameters through the
C CellML environment.
          CALL ASSERT(NMNO(1,0,nx_opt).LE.NMMX_LOCAL,'>>NMMX_LOCAL '
     &      //'too small',ERROR,*9999)
C newe
        ELSE IF(PARAMTYPE(1:20).EQ.'GEOMETRIC_PARAMETERS') THEN
          CALL ASSERT(CALL_SOLV,'>>Solution mapping arrays not set up:'
     '      //' use define solve',ERROR,*9999)
          CALL NX_LOC(NX_INQUIRE,nxc,nx_sol,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx_sol.GT.0,
     '      '>>No nx defined for this solve class',
     '      ERROR,*9999)
          nx_opt=1
        ELSE IF(PARAMTYPE(1:5).EQ.'SHEAR') THEN
          CALL NX_LOC(NX_INQUIRE,nxc,nx_opt,NX_OPTI,ERROR,*9999)
          CALL ASSERT(nx_opt.GT.0,
     '      '>>No nx defined for this optimisation class',
     '      ERROR,*9999)
          CALL NX_LOC(NX_INQUIRE,nxc,nx_sol,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx_sol.GT.0,
     '      '>>No nx defined for this solve class',
     '      ERROR,*9999)
        ELSE IF(PARAMTYPE(1:20).EQ.'GEOMETRIC_PARAMETERS') THEN
          CALL ASSERT(CALL_SOLV,'>>Solution mapping arrays not set up:'
     '      //' use define solve',ERROR,*9999)
          CALL NX_LOC(NX_INQUIRE,nxc,nx_sol,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx_sol.GT.0,
     '      '>>No nx defined for this solve class',
     '      ERROR,*9999)
          nx_opt=1
        ELSE IF(PARAMTYPE(1:12).EQ.'DATA_FITTING') THEN
          IF(CBBREV(CO,'DEFORMED',2,noco+1,NTCO,N3CO)) THEN
            DEFORMED=.TRUE.
          ELSE
            DEFORMED=.FALSE.
          ENDIF
          CALL PARSE_DATA(nd0,nd1,noco+1,NTCO,CO,ERROR,*9999)
          IF(CBBREV(CO,'FUNCTION',2,noco+1,NTCO,N3CO)) THEN
            MODE=0
          ELSEIF(CBBREV(CO,'JACOBIAN',2,noco+1,NTCO,N3CO)) THEN
            MODE=1
          ELSEIF(CBBREV(CO,'BOTH',2,noco+1,NTCO,N3CO)) THEN
            MODE=2
          ENDIF
          CALL ASSERT(.NOT.COUPLED_RES,'>>Not implemented for coupled '
     '      //'problems',ERROR,*9999)
          CALL NX_LOC(NX_INQUIRE,nxc,nx_opt,NX_OPTI,ERROR,*9999)
          CALL ASSERT(nx_opt.GT.0,
     '      '>>No nx defined for this optimisation class',
     '      ERROR,*9999)
          CALL NX_LOC(NX_INQUIRE,nxc,nx_sol,NX_SOLVE,ERROR,*9999)
        ELSE IF(PARAMTYPE(1:15).EQ.'TORSO_CUSTOMISE') THEN
          CALL ASSERT(.NOT.COUPLED_RES,'>>Not implemented for coupled '
     '      //'problems',ERROR,*9999)
          CALL NX_LOC(NX_INQUIRE,nxc,nx_opt,NX_OPTI,ERROR,*9999)
          CALL ASSERT(nx_opt.GT.0,
     '      '>>No nx defined for this optimisation class',
     '      ERROR,*9999)
        ELSE IF(PARAMTYPE(1:15).EQ.'MOMENTS') THEN
          CALL ASSERT(.NOT.COUPLED_RES,'>>Not implemented for coupled '
     '      //'problems',ERROR,*9999)
          CALL NX_LOC(NX_INQUIRE,nxc,nx_opt,NX_OPTI,ERROR,*9999)
          CALL ASSERT(nx_opt.GT.0,
     '      '>>No nx defined for this optimisation class',
     '      ERROR,*9999)
        ELSE IF(PARAMTYPE(1:9).EQ.'POTENTIAL') THEN
          CALL ASSERT(EVALUATE_TRANSFER,'>>Evaluate transfer first',
     '      ERROR,*9999)
          CALL ASSERT(EVALUATE_PHI,'>>Evaluate PHI first',
     '      ERROR,*9999)
          CALL NX_LOC(NX_INQUIRE,nxc,nx_opt,NX_OPTI,ERROR,*9999)
          CALL ASSERT(nx_opt.GT.0,
     '      '>>No nx defined for this optimisation class',
     '      ERROR,*9999)

          CALL NX_LOC(NX_INQUIRE,nxc,nx_sol,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx_sol.GT.0,
     '      '>>No nx defined for this solve class',
     '      ERROR,*9999)
          IF(CBBREV(CO,'FUNCTION',2,noco+1,NTCO,N3CO)) THEN
            MODE=0
          ELSEIF(CBBREV(CO,'JACOBIAN',2,noco+1,NTCO,N3CO)) THEN
            MODE=1
          ELSEIF(CBBREV(CO,'BOTH',2,noco+1,NTCO,N3CO)) THEN
            MODE=2
          ENDIF
          IF(CBBREV(CO,'TSTART',2,noco+1,NTCO,N3CO)) THEN
            TOPTI_START=RFROMC(CO(N3CO+1))
          ELSE
            TOPTI_START=ACTN_MIN(2)
          ENDIF
          IF(CBBREV(CO,'TEND',2,noco+1,NTCO,N3CO)) THEN
            TOPTI_END=RFROMC(CO(N3CO+1))
          ELSE
            TOPTI_END=ACTN_MAX(2)
          ENDIF
          ERR=0
          SOPTI_START=CALC_SAMPLE_FROM_TIME(TOPTI_START,ERR,ERROR)
          IF(ERR.NE.0) GOTO 9999
          SOPTI_END=CALC_SAMPLE_FROM_TIME(TOPTI_END,ERR,ERROR)
          IF(ERR.NE.0) GOTO 9999
        ELSE IF((PARAMTYPE(1:6).EQ.'STRESS').OR.
     '    (PARAMTYPE(1:6).EQ.'STRAIN')) THEN
          CALL NX_LOC(NX_INQUIRE,nxc,nx_opt,NX_OPTI,ERROR,*9999)
          CALL ASSERT(nx_opt.GT.0,
     '      '>>No nx defined for this optimisation class',
     '      ERROR,*9999)

          CALL NX_LOC(NX_INQUIRE,nxc,nx_sol,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx_sol.GT.0,
     '      '>>No nx defined for this solve class',
     '      ERROR,*9999)
          IF(CBBREV(CO,'REFERENCE',3,noco+1,NTCO,N3CO)) THEN
            COORDS='Reference'
          ELSE IF(CBBREV(CO,'WALL',3,noco+1,NTCO,N3CO)) THEN
            COORDS='Wall'
          ELSE
            COORDS='Fibre'
          ENDIF
        ENDIF

        IF(CBBREV(CO,'STEP',1,noco+1,NTCO,N3CO)) THEN
          NTSTEP=IFROMC(CO(N3CO+1))
        ELSE
          NTSTEP=0
        ENDIF

        IF(CBBREV(CO,'WITH',2,noco+1,NTCO,N3CO)) THEN
          CALL PARSRL(CO(N3CO+1),20,NTLIST,PBOPTI,ERROR,*9999)
          NTSTEP=0
        ELSE
          IF((PARAMTYPE(1:19).EQ.'MATERIAL_PARAMETERS').OR.
     '      (PARAMTYPE(1:5).EQ.'SHEAR')) THEN
            DO noopti=1,NMNO(1,0,nx_opt)
              PBOPTI(noopti)=PAOPTI(noopti)
            ENDDO
          ENDIF
        ENDIF

        IF(CBBREV(CO,'ANALYT_DERIVS',2,noco+1,NTCO,N3CO)) THEN
          DERIV=.TRUE.
          ANALYTIC=.TRUE.
        ELSE IF(CBBREV(CO,'NUM_DERIVS',1,noco+1,NTCO,N3CO)) THEN
          DERIV=.TRUE.
          ANALYTIC=.FALSE.
          IF(CBBREV(CO,'INTERVAL',1,noco+2,NTCO,N3CO)) THEN
            DELTA=RFROMC(CO(N3CO+1))
          ELSE
            DELTA=1.0d-6
          ENDIF
        ELSE
          DERIV=.FALSE.
          ANALYTIC=.FALSE.
        ENDIF

        IF(PARAMTYPE(1:12).EQ.'DATA_FITTING') THEN
          IF(CBBREV(CO,'VIEW',2,noco+1,NTCO,N3CO)) THEN
            VIEW=.TRUE.
          ELSE
            VIEW=.FALSE.
          ENDIF
        ELSE
          IF(CBBREV(CO,'NOVIEW',2,noco+1,NTCO,N3CO)) THEN
            VIEW=.FALSE.
          ELSE
            VIEW=.TRUE.
          ENDIF
        ENDIF


C LKC 18-DEC-1999 adding LIST_RESID - output residual info.
        IF(CBBREV(CO,'LIST_RESIDUALS',4,noco+1,NTCO,N3CO)) THEN
          LIST_RESID=IFROMC(CO(N3CO+1))
        ELSE
          LIST_RESID=0
        ENDIF


        IF(PARAMTYPE(1:19).EQ.'MATERIAL_PARAMETERS') THEN
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            CALL ASSERT(NOT(1,1,nr,nx_sol).LE.NREM,
     '        '>>NREM too small',ERROR,*9999)
C HS 26/4/04 moved further up
C            CALL ASSERT(NMNO(1,0,nx_opt).LE.NOPM,'>>NOPM too small',
C     '        ERROR,*9999)
C            CALL ASSERT(NTOPTI.LE.NOPM,'>>NOPM too small',
C     '        ERROR,*9999)
            IF(DERIV) THEN
              IF(CBBREV(CO,'PARAMNUM',1,noco+2,NTCO,N3CO)) THEN
                nopara1=IFROMC(CO(N3CO+1))
                nopara2=nopara1
              ELSE !calc and print derivs wrt all params
                nopara1=1
                nopara2=NMNO(1,0,nx_opt)
              ENDIF
            ENDIF

            PARAMS_CHANGED=.FALSE.
            IF(KTYP26.EQ.1)THEN
C ***         Store material parameters temporarily
C OLD AJP 18/3/96
c              CALL ASSERT(NMNO(1,0,nx_opt).LE.20,
c     '          '>>Increase dim. of TEMP',ERROR,*9999)
c              DO noopti=1,NMNO(1,0,nx_opt)
c                nm=NMNO(1,noopti,nx_opt)
c                TEMP(noopti)=CE(nm,NEELEM(1,nr),nx_sol)
cC!!!            NOTE: needs updating for varying case
c              ENDDO !noopti
C news AJP 18/3/96
              DO nmlist=1,NMNO(1,0,nx_opt)
                nm=NMNO(1,nmlist,nx_opt)
C news HS 18-April-2004:
C added this loop for optimisation of material parameters through the
C CellML environment.
                IF((KTYP54(nr).EQ.3).AND.(KTYP3B.EQ.2)) THEN !grid coupling
                  CELL_RCQS_VALUE_LOCAL(nmlist)
     &              =CELL_RCQS_VALUE(nm,1)
C newe
                ELSE
                  IF(ILP(nm,1,nr,nx_sol).EQ.1.OR. !constant spatially
     '               ILP(nm,1,nr,nx_sol).EQ.2) THEN !piecewise const
                    DO noelem=1,NEELEM(0,nr)
                      ne=NEELEM(noelem,nr)
                      CE(nm,ne,nx_opt)=CE(nm,ne,nx_sol)
                    ENDDO !noelem (ne)
                  ELSE IF(ILP(nm,1,nr,nx_sol).EQ.3) THEN !piecewise linear
                    DO nonode=1,NPNODE(0,nr)
                      np=NPNODE(nonode,nr)
                      CP(nm,np,nx_opt)=CP(nm,np,nx_sol)
                    ENDDO !nonode (np)
                  ELSE IF(ILP(nm,1,nr,nx_sol).EQ.4) THEN !Gauss points
C!!!                needs completing
                  ENDIF !ilp
                ENDIF !ktyp54
              ENDDO !nmlist
            ENDIF !ktyp26=1
C newe

            IF(NTSTEP.EQ.0) THEN !eval residuals at current parameters
              IF(KTYP26.EQ.1) THEN !material parameters
                IF(ITYP5(nr,nx_opt).EQ.2.AND.ITYP2(nr,nx_opt).EQ.9)
     '            THEN !act'ion model
                ELSE !all others
                  DO nmlist=1,NMNO(1,0,nx_opt)
                    nm=NMNO(1,nmlist,nx_opt)
C news HS 18-April-2004:
C added this loop for optimisation of material parameters through the
C CellML environment.
                    IF((KTYP54(nr).EQ.3).AND.(KTYP3B.EQ.2)) THEN !grid coupling
                      no=NONM(nm,NEELEM(1,nr),nx_opt)
                      CELL_RCQS_VALUE(nm,1)=PBOPTI(no) 
C newe
                    ELSE
                      IF(ILP(nm,1,nr,nx_sol).EQ.1.OR. !constant spatially
     '                  ILP(nm,1,nr,nx_sol).EQ.2) THEN !piecewise const
                        DO noelem=1,NEELEM(0,nr)
                          ne=NEELEM(noelem,nr)
                          no=NONM(nm,noelem,nx_opt)
                          CE(nm,ne,nx_sol)=PBOPTI(no)
                          IF(ITYP1(nr,nx_sol).EQ.5.AND. !FE50 problem
     '                      KTYP55(nr).EQ.3 !stress fn of fib/trans strain
     '                      .AND.KTYP56(nr).EQ.3) THEN !pole-zero law
C news MPN 28-Apr-93:
C                         If optimising for a pole parameter in the
C                         pole zero law then must update the shear pole
C                         parameters accordingly.
C!!!                      NOTE: if the fibre distribution model changes
C!!!                      this will need to be changed also
                            IF(nm.EQ.2.OR.nm.EQ.5.OR.nm.EQ.8) THEN
                              IF(nm.EQ.2) THEN !pole(1,1)
                                nm1=11 !update pole(1,2)
                                nm2=14 !update pole(1,3)
                              ELSE IF(nm.EQ.5) THEN !pole(2,2)
                                nm1=17 !update pole(2,1)
                                nm2=20 !update pole(2,3)
                              ELSE IF(nm.EQ.8) THEN !pole(3,3)
                                nm1=23 !update pole(3,1)
                                nm2=26 !update pole(3,2)
                              ENDIF
                              CE(nm1,ne,nx_sol)=2.0d0*CE(nm,ne,nx_sol)
     '                          /DSQRT(1.0d0+2.d00*CE(nm,ne,nx_sol))
                              CE(nm2,ne,nx_sol)=2.0d0*CE(nm,ne,nx_sol)
     '                          /DSQRT(1.0d0+2.0d0*CE(nm,ne,nx_sol))
                            ENDIF
                          ENDIF !ityp1
                        ENDDO !noelem (ne)
                      ELSE IF(ILP(nm,1,nr,nx_sol).EQ.3) THEN !p'wise lin
                        DO nonode=1,NPNODE(0,nr)
                          np=NPNODE(nonode,nr)
                          no=NONM(nm,nonode,nx_opt)
C PJH 23/6/98             CP(nm,np,nx_opt)=PBOPTI(no)
                          CP(nm,np,nx_sol)=PBOPTI(no)
                        ENDDO !nonode (np)
                      ENDIF !ILP
                    ENDIF !ktyp54
                  ENDDO !nmlist
                  PARAMS_CHANGED=.TRUE.
                ENDIF !ityp5=2 & ityp2=9
              ENDIF !ktyp26=1

              IF(.NOT.ANALYTIC) THEN
                CALL RESFUN(COORDS,DEFORMED,PARAMTYPE,IBT,IDO,INP,
     '            ISIZE_PHI,ISIZE_TBH,LD,LDR,LGE,LIST_RESID,
     '            MODE,NAN,NBH,NBHF,NBJ,NBJF,nd0,nd1,NEELEM,NEL,NELIST,
     '            NENP,NFF,NFFACE,NGAP,NGLIST,
     '            NHE(1,nx_sol),NHP(1,0,nx_sol),NKB,
     '            NKEF,NKH,NKHE,NKJE,NLL,
     '            NLNO(1,nx_opt),NNB,NNF,NNL,NONL(1,nx_opt),
     '            NONY(0,1,1,0,nx_sol),NP_INTERFACE,NP1OPT,
     '            NPF,NPL,NPLIST3,NPNE,NPNODE,NPNY(0,1,0,nx_sol),
     '            nr,NRE,NRLIST,NSB,NVHE,NVHP,NVJE,NVJL,
     '            NW(1,1,nx_sol),nx_sol,NXI,NYNE,NYNO(0,1,1,0,nx_sol),
     '            NYNP,NYNR(0,0,1,0,nx_sol),PAOPTY,Z_CONT_LIST,AQ,
     '            CE(1,1,nx_sol),CG,CGE(1,1,1,nx_sol),
     '            CONY(0,1,1,0,nx_sol),CP(1,1,nx_sol),CURVCORRECT,DL,
     '            FEXT,FGRAD,LAPL,LAPLSQR,PAOPTI,PG,PHI,PHI_H,RE1,
     '            RESID,RESJAC,RG,SE,T_BH,WG,WK1_INV,WU,
     '            XA,XE,XG,XID,
     '            XIG,XN,XP,YG,YGF,YP(1,1,nx_sol),ZA,ZA1,Z_CONT,ZD,ZE,
     '            ZE1,ZG,ZP,ZP1,COUPLED_RES,FIX,ERROR,*9999)
              ENDIF

              IF(DERIV) THEN !calc deriv of residuals wrt mat.l params
                IF(ANALYTIC) THEN !Calc derivatives analytically
                  CALL D_RESFUN(PARAMTYPE,IBT,IDO,INP,LGE,NAN,NBH,NBJ,
     '              NBJF,NEELEM,NFF,NGAP,NHE(1,nx_sol),NHP(1,nr,nx_sol),
     '              NKEF,NKH(1,1,1,nr),NKHE,NKJE,NMNO(1,0,nx_opt),
     '              NNF,NPF,NPNE,NPNY(0,1,0,nx_sol),
     '              NPNODE,nr,NRE,NVHE,NVHP(1,1,1,nr),NVJE,
     '              NW(1,1,nx_sol),nx_sol,NXI,
     '              NYNE,NYNP,NYNR(0,0,1,nr,nx_sol),
     '              CE(1,1,nx_sol),CG,CGE(1,1,1,nx_sol),
     '              CP(1,1,nx_sol),CURVCORRECT,
     '              D_RE,D_RI3,D_RP,D_TG,D_ZG,ES,FEXT,
     '              FIX(1,1,nx_sol),PG,RE1,RE2,RG,SE,WG,
     '              XA,XE,XG,XP,YG,YP(1,1,nx_sol),
     '              ZA,ZE,ZE1,ZG,ZG1,ZP,ERROR,*9999)
                ELSE !Calc derivs wrt opt params using FD's
C                 Temporarily store residuals in YP(ny,8)
                  DO no_nynr=1,NYNR(0,1,1,nr,nx_sol) !loop over rows
                    ny=NYNR(no_nynr,1,1,nr,nx_sol) !is row number
                    IF(NPNY(0,ny,0,nx_sol).EQ.1) THEN
                      np=NPNY(4,ny,0,nx_sol)
C GMH 8/1/97 Update cmgui link
                      CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                    ENDIF
                    YP(ny,8,nx_sol)=YP(ny,4,nx_sol)
                  ENDDO !no_nynr
                  DO nopara=nopara1,nopara2
                    nm=NMNO(1,nopara,nx_opt)
                    IF(KTYP26.EQ.1) THEN
                      DO noelem=1,NEELEM(0,nr) !perturb mat param
                        ne=NEELEM(noelem,nr) !..by delta
                        CE(nm,ne,nx_sol)=PBOPTI(nopara)+DELTA
                        IF(ITYP1(nr,nx_sol).EQ.5.AND. !FE50 problem
     '                    KTYP55(nr).EQ.3 !stress fn of fib/trans strain
     '                    .AND.KTYP56(nr).EQ.3) THEN !pole-zero law
C news MPN 28-Apr-93:
C                         If optimising for a pole parameter in the
C                         pole zero law then must update the shear pole
C                         parameters accordingly.
C!!!                      NOTE: if the fibre distribution model changes
C!!!                      this will need to be changed also
                          IF(nm.EQ.2.OR.nm.EQ.5.OR.nm.EQ.8) THEN
                            IF(nm.EQ.2) THEN !pole(1,1)
                              nm1=11 !update pole(1,2)
                              nm2=14 !update pole(1,3)
                            ELSE IF(nm.EQ.5) THEN !pole(2,2)
                              nm1=17 !update pole(2,1)
                              nm2=20 !update pole(2,3)
                            ELSE IF(nm.EQ.8) THEN !pole(3,3)
                              nm1=23 !update pole(3,1)
                              nm2=26 !update pole(3,2)
                            ENDIF
                            CE(nm1,ne,nx_sol)=2.d0*CE(nm,ne,nx_sol)
     '                        /DSQRT(1.d0+2.d0*CE(nm,ne,nx_sol))
                            CE(nm2,ne,nx_sol)=2.d0*CE(nm,ne,nx_sol)
     '                        /DSQRT(1.d0+2.d0*CE(nm,ne,nx_sol))
                          ENDIF
                        ENDIF !pole-zero law
C newe
                      ENDDO !noelem
                    ENDIF
                    CALL RESFUN(COORDS,DEFORMED,PARAMTYPE,IBT,IDO,INP,
     '                ISIZE_PHI,ISIZE_TBH,LD,LDR,
     '                LGE,LIST_RESID,
     '                MODE,NAN,NBH,NBHF,NBJ,NBJF,nd0,nd1,NEELEM,NEL,
     '                NELIST,
     '                NENP,NFF,NFFACE,NGAP,NGLIST,
     '                NHE(1,nx_sol),NHP(1,0,nx_sol),NKB,
     '                NKEF,NKH(1,1,1,nr),NKHE,NKJE,NLL,NLNO(1,nx_opt),
     '                NNB,NNF,NNL,NONL(1,nx_opt),NONY(0,1,1,0,nx_sol),
     '                NP_INTERFACE,NP1OPT,
     '                NPF,NPL,NPLIST3,NPNE,NPNODE,NPNY(0,1,0,nx_sol),
     '                nr,NRE,NRLIST,NSB,NVHE,NVHP,NVJE,NVJL,
     '                NW(1,1,nx_sol),nx_sol,NXI,NYNE,
     '                NYNO(0,1,1,0,nx_sol),NYNP,
     '                NYNR(0,0,1,0,nx_sol),PAOPTY,Z_CONT_LIST,
     '                AQ,CE(1,1,nx_sol),CG,CGE(1,1,1,nx_sol),
     '                CONY(0,1,1,0,nx_sol),CP(1,1,nx_sol),CURVCORRECT,
     '                DL,FEXT,FGRAD,LAPL,LAPLSQR,PAOPTI,PG,PHI,PHI_H,
     '                RE1,RESID,RESJAC,RG,SE,T_BH,WG,
     '                WK1_INV,WU,XA,XE,
     '                XG,XID,XIG,XN,XP,YG,YGF,YP(1,1,nx_sol),ZA,ZA1,
     '                Z_CONT,
     '                ZD,ZE,ZE1,ZG,ZP,ZP1,COUPLED_RES,FIX,
     '                ERROR,*9999)
                    IF(KTYP26.EQ.1) THEN
                      DO noelem=1,NEELEM(0,nr) !reset mat param
                        ne=NEELEM(noelem,nr)
                        CE(nm,ne,nx_sol)=PBOPTI(nopara)
                        IF(ITYP1(nr,nx_sol).EQ.5.AND. !FE50 problem
     '                    KTYP55(nr).EQ.3 !stress fn of fib/trans strain
     '                    .AND.KTYP56(nr).EQ.3) THEN !pole-zero law
C news MPN 28-Apr-93:
C                         If optimising for a pole parameter in the
C                         pole zero law then must update the shear pole
C                         parameters accordingly.
C!!!                      NOTE: if the fibre distribution model changes
C!!!                      this will need to be changed also
                          IF(nm.EQ.2.OR.nm.EQ.5.OR.nm.EQ.8) THEN
                            IF(nm.EQ.2) THEN !pole(1,1)
                              nm1=11 !update pole(1,2)
                              nm2=14 !update pole(1,3)
                            ELSE IF(nm.EQ.5) THEN !pole(2,2)
                              nm1=17 !update pole(2,1)
                              nm2=20 !update pole(2,3)
                            ELSE IF(nm.EQ.8) THEN !pole(3,3)
                              nm1=23 !update pole(3,1)
                              nm2=26 !update pole(3,2)
                            ENDIF
                            CE(nm1,ne,nx_sol)=2.d0*CE(nm,ne,nx_sol)
     '                        /DSQRT(1.d0+2.d0*CE(nm,ne,nx_sol))
                            CE(nm2,ne,nx_sol)=2.d0*CE(nm,ne,nx_sol)
     '                        /DSQRT(1.d0+2.d0*CE(nm,ne,nx_sol))
                          ENDIF
                        ENDIF !pole-zero law
C newe
                      ENDDO !noelem
                    ENDIF
                    DO no_nynr=1,NYNR(0,1,1,nr,nx_sol) !loop over rows
                      ny=NYNR(no_nynr,1,1,nr,nx_sol) !is row number
                      D_RP(ny,nopara)=(YP(ny,4,nx_sol)-YP(ny,8,nx_sol))
     '                  /DELTA
                    ENDDO !no_nynr
                  ENDDO !nopara
                ENDIF
                DO no_nynr=1,NYNR(0,1,1,nr,nx_sol) !loop over rows
                  ny=NYNR(no_nynr,1,1,nr,nx_sol) !is row number
                  DO noy=1,NONY(0,ny,2,nr,nx_sol)
                    no=NONY(noy,ny,2,nr,nx_sol)
                    DO noopti=1,NMNO(1,0,nx_opt)
                      RESJAC(no,noopti)=D_RP(ny,noopti)
                    ENDDO !noopti
                  ENDDO !noy
                ENDDO !no_nynr

C old MPN 20/10/93 - not needed ????
c            ELSE IF((ktyp26.ne.1.or.ktyp27.ne.2).AND. !not ssq re diffs
c     '             (ktyp26.ne.2.or.ktyp27.ne.6).AND.  !not fl int cond.
c     '             (ktyp27.ne.3).AND.                 !not 0 flux diffs
c     '             (ktyp27.ne.4))THEN                 !not hyd pr cond.
C     ...all but KTYP27=3 use routine E04UPF (has its own print monitor)
c                WRITE(OP_STRING,'('' Objective function = '',E11.3,'
c     '         //''' at parameters:'',8E11.3)')
c     '         FUNC,(PALIST(noopti),noopti=1,NTOPTI)
c                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF

              IF(DOP.OR.VIEW) THEN
                IF(DERIV) THEN !print 1st derivs of resids wrt opt pars
                  IF(ANALYTIC) THEN
                    WRITE(OP_STRING,'(/'' Analytic derivatives of '
     '                //'residuals:'')')
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  ELSE
                    WRITE(OP_STRING,'(/'' Finite Difference '
     '                //'derivatives of residuals:'')')
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  ENDIF
                  WRITE(OP_STRING,'('' ..at parameters:  '','
     '              //'5D12.4/,(19X,5D12.4))')
     '              (PBOPTI(noopti),noopti=1,NMNO(1,0,nx_opt))
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'('' ..wrt opt param nm:    '','
     '              //'5(I2,10X)/,(24X,5(I2,10X)))')
     '              (nopara,nopara=nopara1,nopara2)
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  DO no=1,NOT(1,1,nr,nx_sol)
                    WRITE(OP_STRING,'('' RESJAC('',I4,'
     '                //''',nm) = '',5D12.4/,(19X,5D12.4))')
     '                no,(RESJAC(no,nopara),nopara=nopara1,nopara2)
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  ENDDO !no
                ELSE !print residuals
                  WRITE(OP_STRING,'(/'' Resids at opt pars:'','
     '              //'5D12.5/,(20X,5D12.5))')
     '              (PBOPTI(noopti),noopti=1,NTOPTI)
C old ajp 19/3/96  (PBOPTI(noopti),noopti=1,NMNO(1,0,nx_opt))
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  SUMSQR=0.0d0
                  SUMSQreactions=0.d0
                  IF(ITYP5(nr,nx_sol).EQ.2.AND.ITYP2(nr,nx_sol)
     '              .EQ.9) THEN !activation
                    DO nores=1,NT_RES
                      SUMSQR=SUMSQR+RESID(nores)*RESID(nores)
                      WRITE(OP_STRING,'('' RESID('',I5,'') = '','
     '                  //'D12.4)') nores,RESID(nores)
                    ENDDO !nores
                    WRITE(OP_STRING,'('' Sum of squared '','
     '                //'''residuals = '',D12.4)') SUMSQR
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  ELSE
C PJH 4Mar96        DO no=1,NOT(1,1,nr,nx_opt)
                    DO no=1,NOT(1,1,nr,nx_sol)
                      SUMSQR=SUMSQR+RESID(no)*RESID(no)
                      ny=NYNO(1,no,1,nr,nx_sol)
C??? cpb 28/3/96 Check this iy=5
                      SUMSQreactions=SUMSQreactions+YP(ny,5,nx_sol)**2
                      IF(NPNY(0,ny,1,nx_sol).EQ.1) THEN !nodal variable
                        WRITE(OP_STRING,'('' RESID('',I5,'') = '','
     '                    //'D12.4,'' at ny='',I5,'' nk='',I2,'
     '                    //''' nv='',I2,'' nh='',I2,'' np='',I4,'
     '                    //''' nc='',I2)')
     '                    no,RESID(no),ny,(NPNY(i,ny,1,nx_sol),i=1,5)
                        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                      ELSE IF(NPNY(0,ny,1,nx_sol).EQ.2) THEN !aux var
                        WRITE(OP_STRING,'('' RESID('',I5,'') = '','
     '                    //'D12.4,'' at ny='',I5,'' na='',I2,'
     '                    //''' nh='',I2,'' nc='',I2,'' ne='',I4)')
     '                    no,RESID(no),ny,(NPNY(i,ny,1,nx_sol),i=1,4)
                        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                      ENDIF
                    ENDDO !no
                    WRITE(OP_STRING,'('' Sum of squared residuals='','
     '                //'D12.4,'' Sum of squared reactions='',D12.4)')
     '                SUMSQR,SUMSQreactions
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                    IF(SUMSQreactions.GT.1.0d-6) THEN
                      WRITE(OP_STRING,'('' %VAF='',F6.2)')
     '                  100.d0*(1.0d0-SUMSQR/SUMSQreactions)
                      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                    ENDIF
                  ENDIF !ITYP5/2
                ENDIF !DERIV
              ENDIF !DOP.or.VIEW

            ELSE IF(NTSTEP.GT.0) THEN !eval resids at stepped params
              nostep=0
              DO n1step=1,NTSTEP
C ***           Update material params from optimisation parameters
                P1STEP=(PMAX(1)-PMIN(1))/DBLE(NTSTEP-1)
                DO noelem=1,NEELEM(0,1)
                  ne=NEELEM(noelem,1)
                  CE(NMNO(1,1,nx_opt),ne,nx_sol)=PMIN(1)
     '              +DBLE(n1step-1)*P1STEP
                ENDDO !noelem
                DO n2step=1,NTSTEP
                  P2STEP=(PMAX(2)-PMIN(2))/DBLE(NTSTEP-1)
                  DO noelem=1,NEELEM(0,nr)
                    ne=NEELEM(noelem,nr)
                    CE(NMNO(1,2,nx_opt),ne,nx_sol)=PMIN(2)
     '                +DBLE(n2step-1)*P2STEP
                  ENDDO !noelem
                  nostep=nostep+1
                  
                  CALL RESFUN(COORDS,DEFORMED,PARAMTYPE,IBT,IDO,INP,
     '              ISIZE_PHI,ISIZE_TBH,LD,LDR,
     '              LGE,LIST_RESID,MODE,NAN,
     '              NBH,NBHF,NBJ,NBJF,nd0,nd1,NEELEM,NEL,NELIST,
     '              NPNE,NFF,NFFACE,NGAP,NGLIST,NHE(1,nx_sol),
     '              NHP(1,0,nx_sol),NKB,NKEF,
     '              NKH,NKHE,NKJE,NLL,NLNO(1,nx_opt),NNB,NNF,
     '              NNL,NONL(1,nx_opt),
     '              NONY(0,1,1,0,nx_sol),NP_INTERFACE,NP1OPT,
     '              NPF,NPL,NPLIST3,NPNE,NPNODE,NPNY(0,1,0,nx_sol),
     '              nr,NRE,NRLIST,NSB,NVHE,NVHP,NVJE,NVJL,
     '              NW(1,1,nx_sol),nx_sol,NXI,NYNE,
     '              NYNO(0,1,1,0,nx_sol),NYNP,
     '              NYNR(0,0,1,0,nx_sol),PAOPTY,Z_CONT_LIST,
     '              AQ,CE(1,1,nx_sol),CG,CGE(1,1,1,nx_sol),
     '              CONY(0,1,1,0,nx_sol),
     '              CP(1,1,nx_sol),CURVCORRECT,DL,FEXT,
     '              FGRAD,LAPL,LAPLSQR,
     '              PAOPTI,PG,PHI,PHI_H,RE1,RESID,
     '              RESJAC,RG,SE,T_BH,WG,WK1_INV,WU,
     '              XA,XE,XG,XID,XIG,XN,XP,
     '              YG,YGF,YP(1,1,nx_sol),
     '              ZA,ZA1,Z_CONT,ZD,ZE,ZE1,ZG,ZP,ZP1,
     '              COUPLED_RES,FIX,ERROR,*9999)

                  IF(DOP.OR.VIEW) THEN
                    WRITE(OP_STRING,'(/'' Step'',I5,'
     '                //''':  Parameters:'','
     '                //'8D11.3)') nostep,CE(NMNO(1,1,nx_opt),1,nx_sol),
     '                CE(NMNO(1,2,nx_opt),1,nx_sol)
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                    WRITE(OP_STRING,'('' Residuals:'')')
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                    SUMSQR=0.0d0
                    DO no=1,NOT(1,1,nr,nx_sol)
                      SUMSQR=SUMSQR+RESID(no)*RESID(no)
                      WRITE(OP_STRING,'('' RESID('',I4,'') = '','
     '                  //'D12.4)') no,RESID(no)
                      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                    ENDDO
                    WRITE(OP_STRING,'('' Sum of squares of '','
     '                //'''residuals: '',D12.4)') SUMSQR
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDDO !n2step
              ENDDO !n1step
            ENDIF

            IF(KTYP26.EQ.1) THEN !Return mat params to original values
              DO nmlist=1,NMNO(1,0,nx_opt)
                nm=NMNO(1,nmlist,nx_opt)
C news HS 18-April-2004:
C added this loop for optimisation of material parameters through the
C CellML environment.
                IF((KTYP54(nr).EQ.3).AND.(KTYP3B.EQ.2)) THEN !grid coupling
                  CELL_RCQS_VALUE(nm,1)=
     &              CELL_RCQS_VALUE_LOCAL(nmlist)
C newe
                ELSE
                  IF(ILP(nm,1,nr,nx_sol).EQ.1.OR. !constant spatially
     '              ILP(nm,1,nr,nx_sol).EQ.2) THEN !piecewise const
                    DO noelem=1,NEELEM(0,nr)
                      ne=NEELEM(noelem,nr)
                      CE(nm,ne,nx_sol)=CE(nm,ne,nx_opt)
                      IF(ITYP1(nr,nx_sol).EQ.5.AND. !FE50 problem
     '                  KTYP55(nr).EQ.3 !stress fn of fib/trans strain
     '                  .AND.KTYP56(nr).EQ.3) THEN !pole-zero law
C news MPN 28-Apr-93:
C                     If optimising for a pole parameter in the
C                     pole zero law then must update the shear pole
C                     parameters accordingly.
C!!!                  NOTE: if the fibre distribution model changes
C!!!                  this will need to be changed also
                        IF(nm.EQ.2.OR.nm.EQ.5.OR.nm.EQ.8) THEN
                          IF(nm.EQ.2) THEN !pole(1,1)
                            nm1=11 !update pole(1,2)
                            nm2=14 !update pole(1,3)
                          ELSE IF(nm.EQ.5) THEN !pole(2,2)
                            nm1=17 !update pole(2,1)
                            nm2=20 !update pole(2,3)
                          ELSE IF(nm.EQ.8) THEN !pole(3,3)
                            nm1=23 !update pole(3,1)
                            nm2=26 !update pole(3,2)
                          ENDIF
                          CE(nm1,ne,nx_sol)=2.d0*CE(nm,ne,nx_sol)
     '                      /DSQRT(1.d0+2.d0*CE(nm,ne,nx_sol))
                          CE(nm2,ne,nx_sol)=2.d0*CE(nm,ne,nx_sol)
     '                      /DSQRT(1.d0+2.d0*CE(nm,ne,nx_sol))
                        ENDIF
                      ENDIF !ityp1=5 & ktype55=3
C newe
                    ENDDO !noelem
                  ELSE IF(ILP(nm,1,nr,nx_sol).EQ.3) THEN !piecewise linear
                    DO nonode=1,NPNODE(0,nr)
                      np=NPNODE(nonode,nr)
                      CP(nm,np,nx_sol)=CP(nm,np,nx_opt)
                    ENDDO !nonode (np)
                  ELSE IF(ILP(nm,1,nr,nx_sol).EQ.4) THEN !Gauss points
C!!!                needs completing
                  ENDIF !ilp
                ENDIF !ktyp54=3 & ktyp3b=2, grid coupling
              ENDDO !nmlist
              PARAMS_CHANGED=.FALSE.
            ENDIF
          ENDDO !no_nrlist

        ELSE IF(PARAMTYPE(1:5).EQ.'SHEAR') THEN
            nr=1  
            CALL RESFUN(COORDS,DEFORMED,PARAMTYPE,IBT,IDO,INP,
     '        ISIZE_PHI,ISIZE_TBH,LD,LDR,LGE,LIST_RESID,
     '        MODE,NAN,NBH,NBHF,NBJ,NBJF,nd0,nd1,NEELEM,NEL,NELIST,
     '        NENP,NFF,NFFACE,NGAP,NGLIST,
     '        NHE(1,nx_sol),NHP(1,0,nx_sol),NKB,
     '        NKEF,NKH,NKHE,NKJE,NLL,
     '        NLNO(1,nx_opt),NNB,NNF,NNL,NONL(1,nx_opt),
     '        NONY(0,1,1,0,nx_sol),NP_INTERFACE,NP1OPT,
     '        NPF,NPL,NPLIST3,NPNE,NPNODE,NPNY(0,1,0,nx_sol),
     '        nr,NRE,NRLIST,NSB,NVHE,NVHP,NVJE,NVJL,
     '        NW(1,1,nx_sol),nx_sol,NXI,NYNE,NYNO(0,1,1,0,nx_sol),
     '        NYNP,NYNR(0,0,1,0,nx_sol),
     '        PAOPTY,Z_CONT_LIST,AQ,CE(1,1,nx_sol),CG,
     '        CGE(1,1,1,nx_sol),CONY(0,1,1,0,nx_sol),
     '        CP(1,1,nx_sol),CURVCORRECT,DL,FEXT,FGRAD,
     '        LAPL,LAPLSQR,
     '        PAOPTI,PG,PHI,PHI_H,RE1,RESID,RESJAC,RG,
     '        SE,T_BH,WG,WK1_INV,WU,XA,XE,XG,XID,XIG,XN,
     '        XP,YG,YGF,
     '        YP(1,1,nx_sol),ZA,ZA1,Z_CONT,ZD,ZE,ZE1,ZG,ZP,ZP1,
     '        COUPLED_RES,FIX,ERROR,*9999)

        ELSE IF(PARAMTYPE(1:20).EQ.'GEOMETRIC_PARAMETERS') THEN
          IF(COUPLED_RES) THEN !Eval resids for coupled prob
            CALL ASSERT(.NOT.DERIV,'>>Derivs not implemented '
     '        //'for coupled problems',ERROR,*9999)
            nr_coup=0 !region number for coupled problem
C news MPN 24Dec1999: added asserts
            CALL ASSERT(NOT(1,1,nr_coup,nx_sol).LE.NREM,
     '        '>>NREM too small',ERROR,*9999)
            CALL ASSERT(NOT(2,1,nr_coup,nx_sol).LE.NOPM,
     '        '>>NOPM too small',ERROR,*9999)
C end new
            CALL RESFUN(COORDS,DEFORMED,PARAMTYPE,IBT,IDO,INP,
     '        ISIZE_PHI,ISIZE_TBH,LD,LDR,LGE,LIST_RESID,
     '        MODE,NAN,NBH,NBHF,NBJ,NBJF,nd0,nd1,NEELEM,NEL,NELIST,
     '        NENP,NFF,NFFACE,NGAP,NGLIST,NHE(1,nx_sol),
     '        NHP(1,0,nx_sol),NKB,NKEF,
     '        NKH,NKHE,NKJE,NLL,NLNO(1,nx_opt),NNB,NNF,NNL,
     '        NONL(1,nx_opt),NONY(0,1,1,0,nx_sol),
     '        NP_INTERFACE,NP1OPT,
     '        NPF,NPL,NPLIST3,NPNE,NPNODE,NPNY(0,1,0,nx_sol),
     '        nr_coup,NRE,NRLIST,NSB,
     '        NVHE,NVHP,NVJE,NVJL,NW(1,1,nx_sol),nx_sol,NXI,
     '        NYNE,NYNO(0,1,1,0,nx_sol),
     '        NYNP,NYNR(0,0,1,0,nx_sol),PAOPTY,Z_CONT_LIST,
     '        AQ,CE(1,1,nx_sol),CG,CGE(1,1,1,nx_sol),
     '        CONY(0,1,1,0,nx_sol),
     '        CP(1,1,nx_sol),CURVCORRECT,DL,FEXT,FGRAD,
     '        LAPL,LAPLSQR,
     '        PAOPTI,PG,PHI,PHI_H,RE1,RESID,RESJAC,RG,SE,T_BH,
     '        WG,
     '        WK1_INV,WU,
     '        XA,XE,XG,XID,XIG,XN,XP,YG,YGF,YP(1,1,nx_sol),
     '        ZA,ZA1,Z_CONT,ZD,ZE,ZE1,ZG,ZP,ZP1,
     '        COUPLED_RES,FIX,ERROR,*9999)
            IF(DOP.OR.VIEW) THEN
              WRITE(OP_STRING,'(/'' Residuals for current set of '
     '          //'geometric variables:'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              SUMSQR=0.0d0
              DO no=1,NOT(1,1,nr_coup,nx_sol)
                SUMSQR=SUMSQR+RESID(no)*RESID(no)
C!!!            Note only printing indices for FIRST ny coupled to no
                ny=NYNO(1,no,1,nr_coup,nx_sol)
                IF(NPNY(0,ny,1,nx_sol).EQ.1) THEN !nodal variable
                  WRITE(OP_STRING,'('' RESID('',I5,'') = '','
     '              //'D12.4,'' at ny='',I5,'' nk='',I2,'
     '              //''' nv='',I2,'' nh='',I2,'' np='',I4,'
     '              //''' nc='',I2,'' nr='',I2)')
     '              no,RESID(no),ny,(NPNY(i,ny,1,nx_sol),i=1,6)
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ELSE IF(NPNY(0,ny,1,nx_sol).EQ.2) THEN !element var
                  WRITE(OP_STRING,'('' RESID('',I5,'') = '','
     '              //'D12.4,'' at ny='',I5,'' na='',I2,'' nh='',I2,'
     '              //''' nc='',I2,'' ne='',I4,'' nr='',I2)')
     '              no,RESID(no),ny,(NPNY(i,ny,1,nx_sol),i=1,5)
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ENDIF
              ENDDO !no
              WRITE(OP_STRING,'('' Sum of squares of residuals: '','
     '          //'D12.4)') SUMSQR
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
          ELSE !not evaluating resids for coupled problem
            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              CALL ASSERT(NOT(1,1,nr,nx_sol).LE.NREM,
     '          '>>NREM too small',ERROR,*9999)
              CALL ASSERT(NOT(2,1,nr,nx_sol).LE.NOPM,
     '          '>>NOPM too small',ERROR,*9999)
              IF(.NOT.DERIV) THEN
                CALL RESFUN(COORDS,DEFORMED,PARAMTYPE,IBT,IDO,INP,
     '            ISIZE_PHI,ISIZE_TBH,LD,LDR,LGE,LIST_RESID,
     '            MODE,NAN,NBH,NBHF,NBJ,NBJF,nd0,nd1,NEELEM,NEL,NELIST,
     '            NENP,NFF,NFFACE,NGAP,NGLIST,NHE(1,nx_sol),
     '            NHP(1,0,nx_sol),NKB,NKEF,
     '            NKH,NKHE,NKJE,NLL,NLNO(1,nx_opt),NNB,NNF,NNL,
     '            NONL(1,nx_opt),NONY(0,1,1,0,nx_sol),
     '            NP_INTERFACE,NP1OPT,
     '            NPF,NPL,NPLIST3,NPNE,NPNODE,NPNY(0,1,0,nx_sol),
     '            nr,NRE,NRLIST,NSB,
     '            NVHE,NVHP,NVJE,NVJL,NW(1,1,nx_sol),nx_sol,NXI,
     '            NYNE,NYNO(0,1,1,0,nx_sol),
     '            NYNP,NYNR(0,0,1,0,nx_sol),PAOPTY,Z_CONT_LIST,
     '            AQ,CE(1,1,nx_sol),CG,CGE(1,1,1,nx_sol),
     '            CONY(0,1,1,0,nx_sol),
     '            CP(1,1,nx_sol),CURVCORRECT,DL,FEXT,FGRAD,
     '            LAPL,LAPLSQR,
     '            PAOPTI,PG,PHI,PHI_H,RE1,RESID,RESJAC,RG,SE,T_BH,
     '            WG,
     '            WK1_INV,WU,
     '            XA,XE,XG,XID,XIG,XN,XP,
     '            YG,YGF,YP(1,1,nx_sol),
     '            ZA,ZA1,Z_CONT,ZD,ZE,ZE1,ZG,ZP,ZP1,
     '            COUPLED_RES,FIX,ERROR,*9999)

              ELSE IF(DERIV) THEN !calc deriv of residuals wrt geom vars
                IF(CBBREV(CO,'PARAMNUM',1,noco+2,NTCO,N3CO)) THEN
                  nopara1=IFROMC(CO(N3CO+1))
                  nopara2=nopara1
                ELSE !calc and print derivs wrt all params
                  nopara1=1
                  nopara2=NOT(2,1,nr,nx_sol)
                ENDIF
                KTYPTEMP=KTYP1D
                IF(ANALYTIC) THEN !Calc derivatives analytically
                  KTYP1D=1
                ELSE !Calc derivatives using F.D.'s
                  KTYP1D=2
                ENDIF
                CALL D_RESFUN(PARAMTYPE,IBT,IDO,INP,LGE,
     '            NAN,NBH,NBJ,NBJF,NEELEM,NFF,NGAP,NHE(1,nx_sol),
     '            NHP(1,nr,nx_sol),NKEF,NKH(1,1,1,nr),NKHE,NKJE,
     '            NMNO(1,0,nx_opt),NNF,
     '            NPF,NPNE,NPNY(0,1,0,nx_sol),NPNODE,nr,NRE,
     '            NVHE,NVHP(1,1,1,nr),NVJE,NW(1,1,nx_sol),nx_sol,NXI,
     '            NYNE,NYNP,NYNR(0,0,1,nr,nx_sol),
     '            CE(1,1,nx_sol),CG,CGE(1,1,1,nx_sol),CP(1,1,nx_sol),
     '            CURVCORRECT,D_RE,D_RI3,D_RP,D_TG,D_ZG,ES,
     '            FEXT,FIX(1,1,nx_sol),PG,RE1,RE2,RG,SE,
     '            WG,XA,XE,XG,XP,YG,YP(1,1,nx_sol),
     '            ZA,ZE,ZE1,ZG,ZG1,ZP,ERROR,*9999)
                KTYP1D=KTYPTEMP

                DO no_nynr1=1,NYNR(0,1,1,nr,nx_sol) !loop over rows
                  ny1=NYNR(no_nynr1,1,1,nr,nx_sol) !is row number
                  DO noy1=1,NONY(0,ny1,1,nr,nx_sol) !loop nos 4 eqn ny1
                    no1=NONY(noy1,ny1,1,nr,nx_sol) !no assc with eqn ny1
                    DO no_nynr2=1,NYNR(0,0,1,nr,nx_sol) !loop glob vars
                      ny2=NYNR(no_nynr2,0,1,nr,nx_sol) !is global var #
                      DO noy2=1,NONY(0,ny2,2,nr,nx_sol) !loop nos 4 ny2
                        no2=NONY(noy2,ny2,2,nr,nx_sol) !no assc with ny2
                        RESJAC(no1,no2)=D_RP(ny1,ny2)
                      ENDDO !noy2
                    ENDDO !no_nynr2
                  ENDDO !noy1
                ENDDO !no_nynr1
              ENDIF

              IF(DOP.OR.VIEW) THEN
                IF(DERIV) THEN !print 1st derivs of resids wrt geom vars
                  IF(nopara1.EQ.nopara2) THEN
                    WRITE(CHAR3,'(I3)')nopara1
                    CHAR1='geometric variable #'//CHAR3(1:3)
                  ELSE
                    CHAR1='all geometric variables'
                  ENDIF
                  IF(ANALYTIC) THEN
                    WRITE(OP_STRING,'(/'' Analytic derivatives of '','
     '                //'''residuals wrt '//CHAR1(1:23)//''')')
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  ELSE
                    WRITE(OP_STRING,'(/'' FD derivatives of residuals'
     '                //' wrt '//CHAR1(1:23)//''')')
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  ENDIF
                  DO no=1,NOT(1,1,nr,nx_sol)
                    WRITE(OP_STRING,'('' RESJAC('',I3,'',no=1..):'','
     '                //'5D12.4/,(20X,5D12.4))')
     '                no,(RESJAC(no,noopti),noopti=nopara1,nopara2)
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  ENDDO !no
                ELSE !print residuals
                  WRITE(OP_STRING,'(/'' Residuals for current set of '
     '              //'geometric variables:'')')
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  SUMSQR=0.0d0
                  DO no=1,NOT(1,1,nr,nx_sol)
                    SUMSQR=SUMSQR+RESID(no)*RESID(no)
                    WRITE(OP_STRING,'('' RESID('',I4,'') = '',D12.4)')
     '                no,RESID(no)
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  ENDDO !no
                  WRITE(OP_STRING,'('' Sum of squares of residuals: '','
     '              //'D12.4)') SUMSQR
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ENDIF
              ENDIF
            ENDDO !no_nrlist
          ENDIF

        ELSE IF(PARAMTYPE(1:12).EQ.'DATA_FITTING'.OR.
     '      PARAMTYPE(1:15).EQ.'TORSO_CUSTOMISE'.OR.
     '      PARAMTYPE(1:7).EQ.'MOMENTS') THEN

          IF(KTYP26.EQ.1.AND.KTYP27.EQ.5) THEN
            nx=nx_sol
          ELSE
            nx=nx_opt
          ENDIF

          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            CALL RESFUN(COORDS,DEFORMED,PARAMTYPE,IBT,IDO,INP,
     '        ISIZE_PHI,ISIZE_TBH,LD,LDR,LGE,LIST_RESID,
     '        MODE,NAN,NBH,NBHF,NBJ,NBJF,nd0,nd1,
     '        NEELEM,NEL,NELIST,NENP,
     '        NFF,NFFACE,NGAP,NGLIST,
     '        NHE(1,nx),NHP(1,0,nx),
     '        NKB,NKEF,NKH,NKHE,NKJE,NLL,
     '        NLNO(1,nx),NNB,NNF,NNL,NONL(1,nx),
     '        NONY(0,1,1,0,nx),NP_INTERFACE,NP1OPT,
     '        NPF,NPL,NPLIST3,NPNE,NPNODE,NPNY(0,1,0,nx),
     '        nr,NRE,NRLIST,NSB,
     '        NVHE,NVHP,NVJE,NVJL,NW(1,1,nx),nx,NXI,
     '        NYNE,NYNO(0,1,1,0,nx),
     '        NYNP,NYNR(0,0,1,0,nx),PAOPTY,Z_CONT_LIST,
     '        AQ,CE(1,1,nx),CG,CGE(1,1,1,nx),
     '        CONY(0,1,1,0,nx),
     '        CP(1,1,nx),CURVCORRECT,DL,FEXT,FGRAD,
     '        LAPL,LAPLSQR,
     '        PAOPTI,PG,PHI,PHI_H,RE1,RESID,RESJAC,RG,SE,T_BH,
     '        WG,
     '        WK1_INV,WU,
     '        XA,XE,XG,XID,XIG,XN,XP,YG,YGF,YP(1,1,nx),
     '        ZA,ZA1,Z_CONT,ZD,ZE,ZE1,ZG,ZP,ZP1,COUPLED_RES,FIX,
     '        ERROR,*9999)

C new MPN/VY 30Mar2009
              IF(DOP.OR.VIEW) THEN
                  WRITE(OP_STRING,'(/'' Residuals for current set of '
     '              //'variables:'')')
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  SUMSQR=0.0d0
                  DO no=1,NT_RES
                    SUMSQR=SUMSQR+RESID(no)*RESID(no)
                    WRITE(OP_STRING,'('' RESID('',I4,'') = '',D12.4)')
     '                no,RESID(no)
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  ENDDO !no
                  WRITE(OP_STRING,'('' Sum of squares of residuals: '','
     '              //'D12.4)') SUMSQR
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
C end new
          ENDDO !no_nrlist


        ELSE IF(PARAMTYPE(1:8).EQ.'MAGNETIC') THEN

          IF(CBBREV(CO,'TSTART',2,noco+1,NTCO,N3CO)) THEN
            TOPTI_START=RFROMC(CO(N3CO+1))
          ELSE
C LKC 10-JUL-2002 wrong defaults
C            TOPTI_START=ACTN_MIN(2)
            TOPTI_START=0.D0
          ENDIF
          IF(CBBREV(CO,'TEND',2,noco+1,NTCO,N3CO)) THEN
            TOPTI_END=RFROMC(CO(N3CO+1))
          ELSE
C LKC 10-JUL-2002 wrong defaults
C            TOPTI_END=ACTN_MAX(2)
            TOPTI_END=CALC_TIME_FROM_SAMPLE(ISIZE_MFI(2,1))
          ENDIF
          ERR=0

          IF(CBBREV(CO,'FREQ',3,noco+1,NTCO,N3CO)) THEN
            TRSF_FREQUENCY=RFROMC(CO(N3CO+1))
          ENDIF

          CALL RESFUN_DIPOLE(ISIZE_MFI,ISIZE_PHI,LIST_RESID,
     '      TOPTI_START,TOPTI_END,MFI,PHI,PHI_H,RESID,ERROR,*9999)


        ELSE IF(PARAMTYPE(1:9).EQ.'POTENTIAL') THEN
          nr=TRSF_NR_FIRST
          IF(ACTN_IRESFUN.EQ.1) THEN  ! Normal residual calculation
            CALL RESFUN(COORDS,DEFORMED,PARAMTYPE,IBT,IDO,INP,
     '        ISIZE_PHI,ISIZE_TBH,
     '        LD,LDR,LGE,LIST_RESID,
     '        MODE,NAN,NBH,NBHF,NBJ,NBJF,nd0,nd1,
     '        NEELEM,NEL,NELIST,NENP,
     '        NFF,NFFACE,NGAP,NGLIST,
     '        NHE(1,nx_sol),NHP(1,0,nx_sol),
     '        NKB,NKEF,NKH,NKHE,NKJE,NLL,
     '        NLNO(1,nx_opt),NNB,NNF,NNL,NONL(1,nx_opt),
     '        NONY(0,1,1,0,nx_opt),NP_INTERFACE,NP1OPT,
     '        NPF,NPL,NPLIST3,NPNE,NPNODE,NPNY(0,1,0,nx_sol),
     '        nr,NRE,NRLIST,NSB,
     '        NVHE,NVHP,NVJE,NVJL,NW(1,1,nx_sol),nx_sol,NXI,
     '        NYNE,NYNO(0,1,1,0,nx_opt),
     '        NYNP,NYNR(0,0,1,0,nx_sol),PAOPTY,Z_CONT_LIST,
     '        AQ,CE(1,1,nx_sol),CG,CGE(1,1,1,nx_sol),
     '        CONY(0,1,1,0,nx_opt),
     '        CP(1,1,nx_sol),CURVCORRECT,DL,FEXT,FGRAD,
     '        LAPL,LAPLSQR,
     '        PAOPTI,PG,PHI,PHI_H,RE1,RESID,RESJAC,RG,SE,T_BH,
     '        WG,
     '        WK1_INV,WU,
     '        XA,XE,XG,XID,XIG,XN,XP,YG,YGF,YP(1,1,nx_sol),
     '        ZA,ZA1,Z_CONT,ZD,ZE,ZE1,ZG,ZP,ZP1,COUPLED_RES,FIX,
     '        ERROR,*9999)
          ELSE IF(ACTN_IRESFUN.EQ.2) THEN ! New improved opti version
            IF(DABS(CC_OBJ_WEIGHT).GT.ZERO_TOL) THEN
              CALL RESFUN_ACTN_CC(PARAMTYPE,IBT,
     '          ISIZE_TBH,LIST_RESID,
     '          MODE,NBH,NENP,NPLIST3,NPNE,nr,nx_sol,NXI,
     '          NYNO(0,1,1,0,nx_opt),NYNP,PHI,PHI_H,RESID,RESJAC,
     '          T_BH,WK1_INV,XP,YP(1,1,nx_sol),ERROR,*9999)
            ELSE
              CALL RESFUN_ACTN(ISIZE_TBH,LIST_RESID,
     '          MODE,NPLIST3,nr,NYNO(0,1,1,0,nx_opt),
     '          NYNR(0,0,1,0,nx_sol),LAPL,LAPLSQR,
     '          PHI,PHI_H,RESID,RESJAC,T_BH,YP(1,1,nx_sol),
     '          ERROR,*9999)
            ENDIF
          ENDIF
        ELSE IF((PARAMTYPE(1:6).EQ.'STRESS').OR.
     '    (PARAMTYPE(1:6).EQ.'STRAIN')) THEN
          CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '    ERROR,*9999)
          NGLIST(0)=NGT(NBH(NH_LOC(1,nx_sol),1,NELIST(1)))
          DO no_nglist=1,NGLIST(0)
            NGLIST(no_nglist)=no_nglist
          ENDDO
          nr=1  
          CALL RESFUN(COORDS,DEFORMED,PARAMTYPE,IBT,IDO,INP,
     '      ISIZE_PHI,ISIZE_TBH,LD,LDR,LGE,LIST_RESID,
     '      MODE,NAN,NBH,NBHF,NBJ,NBJF,nd0,nd1,NEELEM,NEL,NELIST,
     '      NENP,NFF,NFFACE,NGAP,NGLIST,
     '      NHE(1,nx_sol),NHP(1,0,nx_sol),NKB,
     '      NKEF,NKH,NKHE,NKJE,NLL,
     '      NLNO(1,nx_opt),NNB,NNF,NNL,NONL(1,nx_opt),
     '      NONY(0,1,1,0,nx_sol),NP_INTERFACE,NP1OPT,
     '      NPF,NPL,NPLIST3,NPNE,NPNODE,NPNY(0,1,0,nx_opt),
     '      nr,NRE,NRLIST,NSB,NVHE,NVHP,NVJE,NVJL,
     '      NW(1,1,nx_sol),nx_sol,NXI,NYNE,NYNO(0,1,1,0,nx_opt),
     '      NYNP,NYNR(0,0,1,0,nx_sol),PAOPTY,Z_CONT_LIST,AQ,
     '      CE(1,1,nx_sol),CG,CGE(1,1,1,nx_sol),
     '      CONY(0,1,1,0,nx_sol),
     '      CP(1,1,nx_sol),CURVCORRECT,DL,FEXT,FGRAD,
     '      LAPL,LAPLSQR,
     '      PAOPTI,PG,PHI,PHI_H,RE1,RESID,RESJAC,RG,
     '      SE,T_BH,WG,WK1_INV,WU,XA,XE,XG,XID,XIG,XN,
     '      XP,YG,YGF,
     '      YP(1,1,nx_sol),ZA,ZA1,Z_CONT,ZD,ZE,ZE1,ZG,ZP,ZP1,
     '      COUPLED_RES,FIX,ERROR,*9999)

        ELSE IF(PARAMTYPE(1:8).EQ.'PRESSURE') THEN
          nr=1
          nx_sol=1
          nx_opt=1 
          CALL RESFUN(COORDS,DEFORMED,PARAMTYPE,IBT,IDO,INP,
     '      ISIZE_PHI,ISIZE_TBH,LD,LDR,LGE,LIST_RESID,
     '      MODE,NAN,NBH,NBHF,NBJ,NBJF,nd0,nd1,NEELEM,NEL,NELIST,
     '      NENP,NFF,NFFACE,NGAP,NGLIST,
     '      NHE(1,nx_sol),NHP(1,0,nx_sol),NKB,
     '      NKEF,NKH,NKHE,NKJE,NLL,
     '      NLNO(1,nx_opt),NNB,NNF,NNL,NONL(1,nx_opt),
     '      NONY(0,1,1,0,nx_sol),NP_INTERFACE,NP1OPT,
     '      NPF,NPL,NPLIST3,NPNE,NPNODE,NPNY(0,1,0,nx_opt),
     '      nr,NRE,NRLIST,NSB,NVHE,NVHP,NVJE,NVJL,
     '      NW(1,1,nx_sol),nx_sol,NXI,NYNE,NYNO(0,1,1,0,nx_opt),
     '      NYNP,NYNR(0,0,1,0,nx_sol),
     '      PAOPTY,Z_CONT_LIST,AQ,CE(1,1,nx_sol),CG,
     '      CGE(1,1,1,nx_sol),CONY(0,1,1,0,nx_sol),
     '      CP(1,1,nx_sol),CURVCORRECT,DL,FEXT,FGRAD,
     '      LAPL,LAPLSQR,
     '      PAOPTI,PG,PHI,PHI_H,RE1,RESID,RESJAC,RG,
     '      SE,T_BH,WG,WK1_INV,WU,XA,XE,XG,XID,XIG,XN,
     '      XP,YG,YGF,
     '      YP(1,1,nx_sol),ZA,ZA1,Z_CONT,ZD,ZE,ZE1,ZG,ZP,ZP1,
     '      COUPLED_RES,FIX,ERROR,*9999)

        ELSE IF(PARAMTYPE(1:8).EQ.'REACTION') THEN
          nr=1
          nx_sol=1
          nx_opt=1 
          CALL RESFUN(COORDS,DEFORMED,PARAMTYPE,IBT,IDO,INP,
     '      ISIZE_PHI,ISIZE_TBH,LD,LDR,LGE,LIST_RESID,
     '      MODE,NAN,NBH,NBHF,NBJ,NBJF,nd0,nd1,NEELEM,NEL,NELIST,
     '      NENP,NFF,NFFACE,NGAP,NGLIST,
     '      NHE(1,nx_sol),NHP(1,0,nx_sol),NKB,
     '      NKEF,NKH,NKHE,NKJE,NLL,
     '      NLNO(1,nx_opt),NNB,NNF,NNL,NONL(1,nx_opt),
     '      NONY(0,1,1,0,nx_sol),NP_INTERFACE,NP1OPT,
     '      NPF,NPL,NPLIST3,NPNE,NPNODE,NPNY(0,1,0,nx_opt),
     '      nr,NRE,NRLIST,NSB,NVHE,NVHP,NVJE,NVJL,
     '      NW(1,1,nx_sol),nx_sol,NXI,NYNE,NYNO(0,1,1,0,nx_opt),
     '      NYNP,NYNR(0,0,1,0,nx_sol),
     '      PAOPTY,Z_CONT_LIST,AQ,CE(1,1,nx_sol),CG,
     '      CGE(1,1,1,nx_sol),CONY(0,1,1,0,nx_sol),
     '      CP(1,1,nx_sol),CURVCORRECT,DL,FEXT,FGRAD,
     '      LAPL,LAPLSQR,
     '      PAOPTI,PG,PHI,PHI_H,RE1,RESID,RESJAC,RG,
     '      SE,T_BH,WG,WK1_INV,WU,XA,XE,XG,XID,XIG,XN,
     '      XP,YG,YGF,
     '      YP(1,1,nx_sol),ZA,ZA1,Z_CONT,ZD,ZE,ZE1,ZG,ZP,ZP1,
     '      COUPLED_RES,FIX,ERROR,*9999)
        ENDIF

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('EVRESI')
      RETURN
 9999 CALL ERRORS('EVRESI',ERROR)
!news MPN 28-Apr-95: reset material params for errors
C *** If an error has occurred during material param opt then reset
C *** material parameter array
      IF(PARAMTYPE(1:19).EQ.'MATERIAL_PARAMETERS') THEN
        IF(KTYP26.EQ.1) THEN
          IF(PARAMS_CHANGED) THEN !Return mat params to original values
            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              DO nmlist=1,NMNO(1,0,nx_opt)
                nm=NMNO(1,nmlist,nx_opt)
                IF(ILP(nm,1,nr,nx_sol).EQ.1.OR. !constant spatially
     '             ILP(nm,1,nr,nx_sol).EQ.2) THEN !piecewise const
                  DO noelem=1,NEELEM(0,nr)
                    ne=NEELEM(noelem,nr)
                    CE(nm,ne,nx_sol)=CE(nm,ne,nx_opt)
                    IF(ITYP1(nr,nx_sol).EQ.5.AND. !FE50 problem
     '                KTYP55(nr).EQ.3 !stress fn of fib/trans strain
     '                .AND.KTYP56(nr).EQ.3) THEN !pole-zero law
C news MPN 28-Apr-93:
C                   If optimising for a pole parameter in the
C                   pole zero law then must update the shear pole
C                   parameters accordingly.
C!!!                NOTE: if the fibre distribution model changes
C!!!                this will need to be changed also
                      IF(nm.EQ.2.OR.nm.EQ.5.OR.nm.EQ.8) THEN
                        IF(nm.EQ.2) THEN !pole(1,1)
                          nm1=11 !update pole(1,2)
                          nm2=14 !update pole(1,3)
                        ELSE IF(nm.EQ.5) THEN !pole(2,2)
                          nm1=17 !update pole(2,1)
                          nm2=20 !update pole(2,3)
                        ELSE IF(nm.EQ.8) THEN !pole(3,3)
                          nm1=23 !update pole(3,1)
                          nm2=26 !update pole(3,2)
                        ENDIF
                        CE(nm1,ne,nx_sol)=2.d0*CE(nm,ne,nx_sol)
     '                    /DSQRT(1.d0+2.d0*CE(nm,ne,nx_sol))
                        CE(nm2,ne,nx_sol)=2.d0*CE(nm,ne,nx_sol)
     '                    /DSQRT(1.d0+2.d0*CE(nm,ne,nx_sol))
                      ENDIF
                    ENDIF
                  ENDDO !noelem
                ELSE IF(ILP(nm,1,nr,nx_sol).EQ.3) THEN !piecewise lin
                  DO nonode=1,NPNODE(0,nr)
                    np=NPNODE(nonode,nr)
                    CP(nm,np,nx_sol)=CP(nm,np,nx_opt)
                  ENDDO !nonode (np)
                ENDIF !ILP
              ENDDO !nmlist
            ENDDO !no_nrlist (nr)
          ENDIF
        ENDIF
      ENDIF
!newe
      CALL EXITS('EVRESI')
      RETURN 1
      END


