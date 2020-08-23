      SUBROUTINE UPNODE(IBT,IDO,INP,NBH,NBJ,NBJF,NEELEM,NEL,NELIST,
     '  NENP,NFF,NFFACE,NHE,NHP,NKJE,NKEF,NKH,NKJ,NLF,NLL,NLLINE,
     '  NNF,NNL,NPF,NP_INTERFACE,NPL,NPLIST,NPNE,NPNF,NPNODE,
     '  NQNP,NQS,NQXI,NRE,NRLIST,NVHP,NVJE,NVJF,NVJL,NVJP,NWP,NWQ,
     '  NXLIST,NXQ,NYNE,NYNP,CE,CP,CURVCORRECT,DF,DL,DLL,PAOPTI,
     '  PBOPTI,PG,RG,SE,SF,SP,WG,XA,XE,XG,XP,XQ,YP,YQ,ZA,ZP,
     '  STRING,ERROR,*)

C#### Subroutine: UPNODE
C###  Description:
C###    UPNODE updates XP(nk,nv,nj,np) array with fitted field values
C###    or updates nodal derivatives to derivs wrt arclength.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'fit000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'nqloc00.inc'
      INCLUDE 'opti00.cmn'
      INCLUDE 'stab00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NEELEM(0:NE_R_M,0:NRM),NEL(0:NELM,NLM),NELIST(0:NEM),
     '  NENP(NPM,0:NEPM,0:NRM),
     '  NFF(6,NEM),NFFACE(0:NF_R_M,NRM),NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),
     '  NKJE(NKM,NNM,NJM,NEM),NKEF(0:4,16,6,NBFM),
     '  NKH(NHM,NPM,NCM,0:NRM),
     '  NKJ(NJM,NPM),NLF(4,NFM),NLL(12,NEM),NLLINE(0:NL_R_M,0:NRM),
     '  NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),NPF(9,NFM),
     '  NP_INTERFACE(0:NPM,0:3),NPL(5,0:3,NLM),NPLIST(0:NPM),
     '  NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),NPNODE(0:NP_R_M,0:NRM),
     '  NQNP(NPM),NQS(NEQM),NQXI(0:NIM,NQSCM),
     '  NRE(NEM),NRLIST(0:NRM),NXQ(-NIM:NIM,0:4,0:NQM,NAM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NVJE(NNM,NBFM,NJM,NEM),
     '  NVJL(4,NJM,NLM),NVJF(NNM,NBFM,NJM),
     '  NVJP(NJM,NPM),NWP(NPM,2),NWQ(8,0:NQM,NAM),
     '  NXLIST(0:NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM,NXM),CP(NMM,NPM,NXM),
     '  CURVCORRECT(2,2,NNM,NEM),
     '  DF(NFM),DL(3,NLM),DLL(3,NLM),
     '  PAOPTI(*),PBOPTI(*),PG(NSM,NUM,NGM,NBM),RG(NGM),
     '  SE(NSM,NBFM,NEM),SF(NSM,NBFM),SP(NKM,NBFM,NPM),WG(NGM,NBM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),
     '  XP(NKM,NVM,NJM,NPM),XQ(NJM,NQM),YP(NYM,NIYM,NXM),
     '  YQ(NYQM,NIQM,NAM,NXM),
     '  ZA(NAM,NHM,NCM,NEM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER I,IBEG,IDRN,IEND,IFROMC,IG(4),il,INT_TEMP(1),N,N3CO,nb,
     '  nb_cubics(0:3),nb_e,nb_f,nb_new,nb_old,ne,ne1,nf,nf_e,ng,NGA,nh,
     '  nhj,ni,niqV,nj,nj_field,nj_field1,nj_field2,
     '  nj_geom,njj,njl,NJTOT,nj_update,nk,NK1,NKLIST(3),
     '  NKLIST_old(4),NKLIST_new(4),NKMAX,nk_new,NKNI(3),nk_old,NKTOT,
     '  nl,NL_COUNT,nl_temp,NL_HOLD(2),nt,NITB,ni2,
     '  nn,nn_e,nn_f,NN_TOT,no_interface,no_nf,noelem,nolist,no_nklist,
     '  no_nl,nonode,NO_DERIV,noopti,np,npp,NP1,NP2,NP3,
     '  NPR,nq,nr,nrr,nr_bem,nr_ext,NTRL,
     '  NUM_NKS(0:3),nv,nv1,nv2,nx,nxc,nx_bem,nx_ext,ny,nyd,nq1,nq2,
     '  NOLINE(0:NJM),has_line
      REAL*8 DERIVT,DERIV_FACTOR(2),DIFF,DQ,DQT,
     '  DS,DST(100),DSDS(3),
     '  DX(3),DXDS(3),DXNDS(3,2),DX_update,FLUX(100,2),G(3),K,SUM,TEMP,
     '  NEW_THETA,NORMAL(3),
     '  TANGENT(3,2),VELOC,W,WG_LOCAL(10),
     '  X_geom,XIGG(10),XIPOS,XN_LOCAL(2,3,4),Y_geom,ZERO_TOLERANCE
      CHARACTER TYPE*9
      LOGICAL ADD_ORTHOG,ALL_REGIONS,ARCLENGTH,ATNODES,CONVERT,
     '  CHECK,CONTINUOUS,DEFORM,DERIV,EXCEPTNODES,
     '  FIT,FOUND,GEOMETRY,GRID,HANGING,
     '  HERMITE,INLIST,INTERFACE,INTERFACE2,
     '  LINEAR,POTENTIAL,REPLACE,RESCAL,
     '  MATERIAL,SC_FAC,SWAP,UPDATE_ALL,ZERO,INDIVIDUAL_VERSIONS
      PARAMETER(ZERO_TOLERANCE=1.0d-8)
      ! Functions
      LOGICAL ABBREV,CBBREV
      DATA NKNI/2,3,5/,NKLIST/2,3,5/,NUM_NKS/1,2,4,8/
      DATA XIGG/0.5000000000000D0,0.2113248654051D0,0.7886751345948D0,
     '          0.1127016653792D0,0.5000000000000D0,0.8872983346207D0,
     '          0.0694318442029D0,0.3300094782075D0,0.6699905217924D0,
     '          0.9305681557970D0/
      DATA WG_LOCAL/1.0000000000000D0,0.5000000000000D0,0.50000000000D0,
     '          0.2777777777778D0,0.4444444444444D0,0.2777777777778D0,
     '          0.1739274225687D0,0.3260725774313D0,0.3260725774313D0,
     '          0.1739274225687D0/
      DATA   IG/0,1,3,6/,NGA/4/

      CALL ENTERS('UPNODE',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM update nodes <(replace/add_linear/add_orthogonal)[replace]>
C###  Parameter:      <in (nj#/all)[all]>
C###  Specify the nj to be updated
C###  Parameter:      <convert>
C###  Specify adding a  rc field to polar geom
C###  Parameter:      <rescale>
C###  Specify wiether scaling factors are recalculated
C###  Parameter:      <(undeformed/deformed)[undeformed]>
C###  Specify wiether deformed or undeformed geometry is used
C###  Parameter:      <swap>
C###  Specify swapping  ZP and XP
C###  Parameter:      <region #s/all[1]>
C###    Specify the region numbers to update.
C###  Description:
C###    Update nodal coordinates at nj=NJ_NUMBER (1, 2, or 3).  Scaling
C###    factors are recalculated if 'rescale' is included. Use of
C###    add_orthog will treat the field as if it was an orthogonal
C###    displacement (e.g. fitting a thickness).

        OP_STRING(1)=STRING(1:IEND)//' <(replace/add_linear/add_or'
     '    //'thogonal)[replace]>'
        OP_STRING(2)=BLANK(1:15)//'<in (nj#/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<convert>'
        OP_STRING(4)=BLANK(1:15)//'<rescale>'
        OP_STRING(5)=BLANK(1:15)//'<(undeformed/deformed)[undeformed]>'
        OP_STRING(6)=BLANK(1:15)//'<swap>'
        OP_STRING(7)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update nodes derivative NUMBER#[1]
C###  Parameter:      <linear>
C###    Specify linear calculation
C###  Parameter:      <wrt (xi/arclength)[arclength]>
C###    Specify whether calculated derivatives are with res[ect to xi or
C###    arclength.  This only has an effect if `linear' is specified or
C###    all derivatives are initially zero.
C###  Parameter:      <node (#s/all)[all]>
C###    Specify the node numbers to update
C###  Parameter:      <in (nj#/all)[all]>
C###    Specify nj numbers to update
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <versions (average/individual)[average]>
C###    Specify whether to use the average of all versions of a node to
C###    calculate the derivative or to consider each version individually.
C###  Description:
C###    Rescales the nodal derivatives.  If `linear' is specified or
C###    all the derivatives are initially zero,
C###    derivatives are calculated on the basis of straight line lengths
C###    between nodes.  Otherwise all derivatives are scaled to be wrt
C###    arc length (cubic Hermite basis).

        OP_STRING(1)=STRING(1:IEND)//' derivative NUMBER#[1]'
        OP_STRING(2)=BLANK(1:15)//'<linear>'
        OP_STRING(3)=BLANK(1:15)//'<wrt (xi/arclength)[arclength]>'
        OP_STRING(4)=BLANK(1:15)//'<node (#s/all)[all]>'
        OP_STRING(5)=BLANK(1:15)//'<in (nj#/all)[all]>'
        OP_STRING(6)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(7)=BLANK(1:15)//'<versions (average/individual)'
     &    //'[average]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update nodes hanging
C###  Parameter:      <node (#s/all)[all]>
C###    Specify the node numbers to update
C###  Parameter:      <element (#s/all)[all]>
C###    Specify the elements numbers to search
C###  Parameter:      <direction (xi#)[1]>
C###    Specify the hanging direction
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Description:
C###    Sets up hanging node information, which
C###    element and face a node hangs on

        OP_STRING(1)=STRING(1:IEND)//' hanging'
        OP_STRING(2)=BLANK(1:15)//'<node (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<elements (#s/all)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<direction (xi#)[1]>'
        OP_STRING(5)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update nodes Hermite
C###  Parameter:      <in nj#[1]>
C###    Specify nj numbers to update
C###  Parameter:      <newbasis nb#[2]>
C###    Specify the new basis function number
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Description:
C###    Given a geometric variable nj described by basis NBJ(nj,ne)
C###    in element ne, check that the node derivatives for every node
C###    is appropriate for a new basis (newbasis #). If not then add
C###    the necessary new deravites to the XP list and rearrange
C###    existing derivatives.
C###    For example: if nj is described by an (Xi1,Xi2,Xi3) basis
C###    of linear/cubic-Hermite/linear (thus nj has 1 deravitive at
C###    each node) and you want to change the interpolation to
C###    an (Xi1,Xi2,Xi3) basis of bicubic-Hermite/linear, then add the
C###    two new required derivatives [dx/ds(1) and d2x/ds(1)ds(2)],
C###    move the existing derivative [dx/ds(2)] to the correct location
C###    (if required) and set the new derivs [df/ds(1) where f =
C###    x,dx/ds(2)] to DELTAf/(line length).

        OP_STRING(1)=STRING(1:IEND)//' Hermite'
        OP_STRING(2)=BLANK(1:15)//'<in nj#[1]>'
        OP_STRING(3)=BLANK(1:15)//'<newbasis nb#[2]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update nodes continuous_internal
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Description:
C###    Updates internal nodes so that continiuty between elements is
C###    preserved.  This is done by replacing multiple versions of nodal
C###    values with the arithmetic means of all the versions.

        OP_STRING(1)=STRING(1:IEND)//' continuous_internal'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
CC KAT 23Jan99: Does nothing
C
CC#### Command: FEM update nodes xi
CC###  Parameter:      <at NODES[region_2_nodes]>
CC###  Parameter:      <region (#s/all)[1]>
CC###    Specify the region numbers to update.
CC###  Description:
C
C        OP_STRING(1)=STRING(1:IEND)//' xi'
C        OP_STRING(2)=BLANK(1:15)//'<at NODES[region_2_nodes]>'
C        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
C        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
C---------------------------------------------------------------------

C#### Command: FEM update nodes interface
C###  Parameter:     <(position/increment/flux/time)[position]>
C###    Specify the nodal quantity to update
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Description: Updates the position,increment,flux or time of
C###  the nodes at the interface

        OP_STRING(1)=STRING(1:IEND)//' interface'
        OP_STRING(2)=BLANK(1:15)
     '    //'<(position/increment/flux/time)[position]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update nodes scale_factor
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Description: Updates the scale factors at the nodes
C###

        OP_STRING(1)=STRING(1:IEND)//' scale_factor'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update nodes fit
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Description:
C###

        OP_STRING(1)=STRING(1:IEND)//' fit'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update nodes grid
C###  Description: Updates the nodes from the grid variables
C###

        OP_STRING(1)=STRING(1:IEND)//' grid'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update nodes potential
C###  Description:
C###    This command is used in coupled bidomain problems to update
C###    the node potential and arc length derivatives from the
C###    extracellular potential.
C###  Parameter:   <from_class #[2]>
C###    Specify the class of the extracellular potential
C###  Parameter:   <to_class #[1]>
C###    Specify the class of the BEM problem
C###  Parameter:   <region (#s/all)[1]>
C###     Specify the regions where first is the BEM region
C###     and the second is the grid reg
C###  Parameter:   <nodes (#s)[1]>
C###     Specify the list of the nodes to update the potential at

        OP_STRING(1)=STRING(1:IEND)//' potential'
        OP_STRING(2)=BLANK(1:15)//'<from_class #[2]>'
        OP_STRING(3)=BLANK(1:15)//'<to_class #[1]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(5)=BLANK(1:15)//'<nodes (#s)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update nodes geometry from grid
C###  Description:
C###    Move boundary element nodes to match grid points in
C###    coupled problems. After a refine, the BE nodes can move
C###    away from the finite element based grid points which have
C###    not been refined. In areas of large curvature this difference
C###    can be significant for coupled grid-bem problems.
C###  Parameter:   <region (#s/all)[1]>
C###     Specify the boundary element region to be adjusted.
C###     and the second is the grid reg
C###  Parameter:   <except nodes (#s)[0]>
C###     The default is to adjust all nodes in the set region. If
C###     there is a cusp or another point where the node should not
C###     be adjusted then the except_nodes qualifier can accept a
C###     list of node numbers to be left out of the adjustment.

        OP_STRING(1)=STRING(1:IEND)//' geometry from grid'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<except nodes (#s)[0]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
C#### Command: FEM update nodes material #[1] from field #[1]
C###  Description:
C###    Updates the nodal material parameter from a field value
C###  Parameter:   <region (#s/all)[1]>
C###    Specify the region numbers to update.

        OP_STRING(1)=STRING(1:IEND)//' material #[1] from field #[1]'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc','UPNODE',ERROR,*9999)
      ELSE
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        nr=NRLIST(1) !may need updating later

        IF(CBBREV(CO,'REPLACE',3,noco+1,NTCO,N3CO)) THEN
          REPLACE=.TRUE.
          ADD_ORTHOG=.FALSE.
        ELSE IF(CBBREV(CO,'ADD_LINEAR',5,noco+1,NTCO,N3CO)) THEN
          REPLACE=.FALSE.
          ADD_ORTHOG=.FALSE.
        ELSE IF(CBBREV(CO,'ADD_ORTHOG',5,noco+1,NTCO,N3CO)) THEN
          REPLACE=.FALSE.
          ADD_ORTHOG=.TRUE.
        ELSE
          REPLACE=.TRUE.
          ADD_ORTHOG=.FALSE.
        ENDIF

        IF(CBBREV(CO,'IN',1,noco+1,NTCO,N3CO)) THEN
C KAT 13May99: User is supplying nj
C          njj=IFROMC(CO(N3CO+1))
C          nj_update=NJ_LOC(NJL_GEOM,njj,nr)
C          CALL ASSERT(nj_update.LE.NJT,'>>nj # not defined',ERROR,*9999)
          nj_update=IFROMC(CO(N3CO+1))
          CALL ASSERT(nj_update.LE.NJ_LOC(0,0,nr),'>>nj # not defined',
     '      ERROR,*9999)
          NJTOT=1
          UPDATE_ALL=.FALSE.
        ELSE
          IF(NJT.EQ.1) THEN
            nj_update=NJ_LOC(NJL_GEOM,1,nr)
          ELSE
            nj_update=0
          ENDIF
          NJTOT=NJT
          UPDATE_ALL=.TRUE.
        ENDIF

        IF(CBBREV(CO,'RESCALE',3,noco+1,NTCO,N3CO)) THEN
          RESCAL=.TRUE.
        ELSE
          RESCAL=.FALSE.
        ENDIF

        IF(CBBREV(CO,'DEFORMED',3,noco+1,NTCO,N3CO)) THEN
          DEFORM=.TRUE.
        ELSE
          DEFORM=.FALSE.
        ENDIF

        IF(CBBREV(CO,'SWAP',3,noco+1,NTCO,N3CO)) THEN
          SWAP=.TRUE.
        ELSE
          SWAP=.FALSE.
        ENDIF

        IF(CBBREV(CO,'LINEAR',3,noco+1,NTCO,N3CO)) THEN
          LINEAR=.TRUE.
        ELSE
          LINEAR=.FALSE.
        ENDIF

        IF(CBBREV(CO,'MATERIAL',3,noco+1,NTCO,N3CO)) THEN
          MATERIAL=.TRUE.
          CALL PARSIL(CO(N3CO+1),1,NTRL,INT_TEMP,ERROR,*9999)
          il=INT_TEMP(1)
          IF(CBBREV(CO,'FIELD',3,noco+1,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),1,NTRL,INT_TEMP,ERROR,*9999)
            nj_field=INT_TEMP(1)
          ELSE
            nj_field=1
          ENDIF
        ELSE
          MATERIAL=.FALSE.
        ENDIF

        IF(CBBREV(CO,'GEOMETRY',4,noco+1,NTCO,N3CO)) THEN
          GEOMETRY=.TRUE.
          IF(CBBREV(CO,'AT',2,noco+1,NTCO,N3CO)) THEN
            ATNODES =.TRUE.
            np1=IFROMC(CO(N3CO+1))
            IF(CBBREV(CO,'FROM',4,noco+1,NTCO,N3CO)) THEN
              np2=IFROMC(CO(N3CO+1))
            ENDIF
          ELSE
            ATNODES=.FALSE.
          ENDIF
        
        ELSE
          GEOMETRY=.FALSE.
        ENDIF

        IF(CBBREV(CO,'DERIVATIVE',3,noco+1,NTCO,N3CO)) THEN
          DERIV=.TRUE.
          HERMITE=.FALSE.
          NO_DERIV=IFROMC(CO(N3CO+1))
          IF(NO_DERIV.EQ.0) THEN
            NO_DERIV=1
          ENDIF
          CALL PARSE_NODES(NPNODE,NPLIST,noco,NRLIST,NTCO,CO,
     '      ERROR,*9999)
          IF(CBBREV(CO,'WRT',3,noco+1,NTCO,N3CO)) THEN
            N3CO=N3CO+1 !argument
            IF(ABBREV(CO(N3CO),'XI',2)) THEN
              ARCLENGTH=.FALSE.
            ELSEIF(ABBREV(CO(N3CO),'ARCLENGTH',3)) THEN
              ARCLENGTH=.TRUE.
            ELSE
              ERROR='Unrecognized option for WRT'
              GOTO 9999
            ENDIF
          ELSE
            ARCLENGTH=.TRUE.
          ENDIF
          IF(ARCLENGTH) THEN
            WRITE(OP_STRING,'('' Update derivs wrt s('',I1,'')'')')
     '        NO_DERIV
          ELSE
            WRITE(OP_STRING,'('' Update derivs wrt Xi('',I1,'')'')')
     '        NO_DERIV
          ENDIF
          INDIVIDUAL_VERSIONS=.FALSE.
          IF(CBBREV(CO,'VERSIONS',3,noco+1,NTCO,N3CO)) THEN
            N3CO=N3CO+1 !argument
            IF(ABBREV(CO(N3CO),'INDIVIDUAL',3)) THEN
              INDIVIDUAL_VERSIONS=.TRUE.
            ENDIF
          ENDIF
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

        ELSE IF(CBBREV(CO,'HERMITE',1,noco+1,NTCO,N3CO)) THEN
          CALL ASSERT(CALL_NODE,'>>No nodes defined',ERROR,*9999)
          CALL ASSERT(CALL_ELEM,'>>No elements defined',ERROR,*9999)
          DERIV=.TRUE.
          HERMITE=.TRUE.
          IF(nj_update.EQ.0) THEN
            nj_update=1
            NJTOT=1
            UPDATE_ALL=.FALSE.
          ENDIF
          IF(CBBREV(CO,'NEWBASIS',1,noco+1,NTCO,N3CO)) THEN
            nb_new=IFROMC(CO(N3CO+1))
          ELSE
            nb_new=2
          ENDIF
          CALL ASSERT(nb_new.LE.NBT,'>>New basis not defined',
     '      ERROR,*9999)
        ELSE
          DERIV=.FALSE.
          HERMITE=.FALSE.
        ENDIF

        IF(CBBREV(CO,'CONTINUOUS_INTERNAL',4,noco+1,NTCO,N3CO)) THEN
          CONTINUOUS=.TRUE.
        ELSE
          CONTINUOUS=.FALSE.
        ENDIF

C KAT 23Jan99: Does nothing
C        IF(CBBREV(CO,'XI',1,noco+1,NTCO,N3CO)) THEN
C          XI=.TRUE.
C        ELSE
C          XI=.FALSE.
C        ENDIF

        IF(CBBREV(CO,'FIT',2,noco+1,NTCO,N3CO)) THEN
          CALL ASSERT(KTYP8.EQ.1,
     '      '>>Update node fit is only defined for geometric fitting',
     '      ERROR,*9999)
          FIT=.TRUE.
        ELSE
          FIT=.FALSE.
        ENDIF

        IF(CBBREV(CO,'INTERFACE',3,noco+1,NTCO,N3CO)) THEN
          INTERFACE=.TRUE.
          IF(CBBREV(CO,'POSITION',1,noco+1,NTCO,N3CO)) THEN
            TYPE='POSITION'
          ELSE IF(CBBREV(CO,'INCREMENT',1,noco+1,NTCO,N3CO)) THEN
            TYPE='INCREMENT'
          ELSE IF(CBBREV(CO,'FLUX',1,noco+1,NTCO,N3CO)) THEN
            TYPE='FLUX'
          ELSE IF(CBBREV(CO,'TIME',1,noco+1,NTCO,N3CO)) THEN
            TYPE='TIME'
          ELSE
            TYPE='POSITION'
          ENDIF
        ELSE
          INTERFACE=.FALSE.
        ENDIF

        IF(CBBREV(CO,'SCALE_FACTOR',2,noco+1,NTCO,N3CO)) THEN
          SC_FAC=.TRUE.
        ELSE
          SC_FAC=.FALSE.
        ENDIF

        IF(CBBREV(CO,'GRID',2,noco+1,NTCO,N3CO)) THEN
          IF(CBBREV(CO,'BEM_REGION',5,noco+1,NTCO,N3CO)) THEN
            nr_bem=IFROMC(CO(N3CO+1))
          ELSE
            nr_bem=1
          ENDIF
          IF(CBBREV(CO,'GRID_REGION',6,noco+1,NTCO,N3CO)) THEN
            nr_ext=IFROMC(CO(N3CO+1))
          ELSE
            nr_ext=2
          ENDIF

          IF(CBBREV(CO,'EXCLUDE',3,noco+1,NTCO,N3CO)) THEN
            NRLIST(0)=2
            NRLIST(1)=nr_bem
            NRLIST(2)=nr_ext
            CALL PARSE_NODES(NPNODE,NPLIST,noco,NRLIST,NTCO,CO,
     '        ERROR,*9999)

            CALL ASSERT(NPLIST(0).GT.0,'>>No nodes entered',ERROR,*9999)

!temporarily write into CPLST until new arrays finalised.
            CPLST(0,1)=NPLIST(0)
            DO np=1,NPLIST(0)
              CPLST(np,1)=NPLIST(np)
            ENDDO
          ENDIF

          GRID=.TRUE.
        ELSE
          CPLST(0,1)=0
          GRID=.FALSE.
        ENDIF

        IF(CBBREV(CO,'HANGING',3,noco+1,NTCO,N3CO)) THEN
          HANGING=.TRUE.
          CALL PARSE_NODES(NPNODE,NPLIST,noco,NRLIST,NTCO,CO,
     '      ERROR,*9999)
          CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '        ERROR,*9999)
          IF(CBBREV(CO,'DIRECTION',3,noco+1,NTCO,N3CO)) THEN
            IDRN=IFROMC(CO(N3CO+1))
          ELSE
            IDRN=1
          ENDIF
        ELSE
          HANGING=.FALSE.
        ENDIF

        IF(CBBREV(CO,'CONVERT',1,noco+1,NTCO,N3CO)) THEN
          CONVERT=.TRUE.
        ELSE
          CONVERT=.FALSE.
        ENDIF

        IF(CBBREV(CO,'POTENTIAL',3,noco+1,NTCO,N3CO)) THEN
! extracellular potential class
          IF(CBBREV(CO,'FROM_CLASS',2,noco+1,NTCO,N3CO)) THEN
            nxc=IFROMC(CO(N3CO+1))
          ELSE
            nxc=3
          ENDIF
          CALL NX_LOC(NX_INQUIRE,nxc,nx_ext,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx_ext.GT.0,'>>No nx defined for GRID class',
     '      ERROR,*9999)

! boundary element class
          IF(CBBREV(CO,'TO_CLASS',2,noco+1,NTCO,N3CO)) THEN
            nxc=IFROMC(CO(N3CO+1))
          ELSE
            nxc=1
          ENDIF
          CALL NX_LOC(NX_INQUIRE,nxc,nx_bem,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx_bem.GT.0,'>>No nx defined for BEM class',
     '      ERROR,*9999)

! set up regions
          CALL ASSERT(NRLIST(0).GT.1,'>>Must specify 2 regions',
     '      ERROR,*9999)
          nr_bem=NRLIST(1)
          nr_ext=NRLIST(2)

! nodes to update
          IF(CBBREV(CO,'NODES',1,noco+1,NTCO,N3CO)) THEN
            CALL PARSE_NODES(NPNODE,NPLIST,noco,NRLIST,NTCO,CO,
     '        ERROR,*9999)
          ELSE !use np_interface to find nodes
            NPLIST(0)=0
            DO npp=1,NPNODE(0,nr_ext)
              np=NPNODE(npp,nr_ext)
              DO nrr=1,NP_INTERFACE(np,0)
                IF(NP_INTERFACE(np,nrr).EQ.nr_bem) THEN
                  NPLIST(0)=NPLIST(0)+1
                  NPLIST(NPLIST(0))=np
                ENDIF
              ENDDO
            ENDDO
          ENDIF

          IF(CBBREV(CO,'CHECK',2,noco+1,NTCO,N3CO)) THEN
            CHECK=.TRUE.
          ELSE
            CHECK=.FALSE.
          ENDIF

          POTENTIAL=.TRUE.
        ELSE
          POTENTIAL=.FALSE.
        ENDIF

C KAT 23Jan99: Does nothing
C        IF(XI) THEN
C          IF(CBBREV(CO,'AT',1,noco+1,NTCO,N3CO)) THEN
C            CALL PARSIL(CO(N3CO+1),NP_R_M,NPLIST(0),NPLIST(1),
C     '        ERROR,*9999)
C          ELSE
C            DO nonode=1,NPNODE(0,2)
C              NPLIST(nonode)=NPNODE(nonode,2)
C            ENDDO
C            NPLIST(0)=NPNODE(0,2)
C          ENDIF
C        ENDIF

        IF(DEFORM.OR.(INTERFACE.AND.TYPE(1:4).EQ.'FLUX').OR.
     '    ADD_ORTHOG) THEN
C         get nx for solve class
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nxc=NXLIST(1)
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          IF(nx.EQ.0) THEN
C           setup nx for solve class
            CALL NX_LOC(NX_ALLOCATE,nxc,nx,NX_SOLVE,ERROR,*9999)
          ENDIF
        ENDIF

        IF(DEFORM.AND..NOT.CALL_INIT) THEN
C ***     28-Oct-89: The following assumes an isoparametric formulation.
C ***     To generalize, need to define dependent variables
C ***     (using define initial, for example) before updating nodes.
          ITYP1(nr,nx)=5 !is finite elasticity
          ITYP6(nr,nx)=2 !is nonlinear
          CALL_INIT=.TRUE.
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            NHE(ne,nx)=NJ_LOC(NJL_GEOM,0,NRE(ne))
            DO nj=1,NJ_LOC(NJL_GEOM,0,NRE(ne))
              NBH(nj,1,ne)=NBJ(nj,ne) !AJP change 18-12-91
            ENDDO !nj
          ENDDO !noelem (ne)
        ENDIF !deform

C*** Cater for options: deriv/Hermite/interface/xi/scale_factor/fit

        IF(DERIV) THEN !Calculate nodal derivs from coord values
          ZERO=.TRUE.
          IF(LINEAR) THEN
C CS 28/12/2000 changing to selected nodes
C            DO nonode=1,NPNODE(0,nr)
C              np=NPNODE(nonode,nr)
            DO nonode=1,NPLIST(0)
              np=NPLIST(nonode)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
              IF(NJTOT.EQ.NJT) THEN
                DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                  nj=NJ_LOC(NJL_GEOM,njj,nr)
                  DO nv=1,NVJP(nj,np)
C KAT 6Dec00: Only deleted the specified derivative.
                    nk=NKNI(NO_DERIV)
                    XP(nk,nv,nj,np)=0.0d0
                  ENDDO !nv
                ENDDO !njj (nj)
              ELSE IF(NJTOT.EQ.1) THEN !just check requested geom var
                DO nv=1,NVJP(nj_update,np)
                  nk=NKNI(NO_DERIV)
                  XP(nk,nv,nj_update,np)=0.0d0
                ENDDO !nv
              ENDIF !NJTOT
            ENDDO !nonode
          ELSE
C           Determine whether derivatives are all zero.
C CS 28/12/2000 changing to selected nodes
C            DO nonode=1,NPNODE(0,nr)
C              np=NPNODE(nonode,nr)
            DO nonode=1,NPLIST(0)
              np=NPLIST(nonode)
              IF(NJTOT.EQ.NJT)THEN
                DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                  nj=NJ_LOC(NJL_GEOM,njj,nr)
                  DO nv=1,NVJP(nj,np)
                    IF(HERMITE) THEN
                      DO nk=2,NKJ(nj,np)
                        IF(DABS(XP(nk,nv,nj,np)).GT.1.0D-8) ZERO=.FALSE.
                      ENDDO
                    ELSE !not HERMITE
                      nk=NKNI(NO_DERIV)
                      IF(DABS(XP(nk,nv,nj,np)).GT.1.0D-8) ZERO=.FALSE.
                    ENDIF !HERMITE
                  ENDDO !nv
                ENDDO !njj (nj)
              ELSE IF(NJTOT.EQ.1)THEN !just check requested geom var
                DO nv=1,NVJP(nj_update,np)
                  IF(HERMITE) THEN
                    DO nk=2,NKJ(nj_update,np)
                      IF(DABS(XP(nk,nv,nj_update,np)).GT.1.0D-8)
     '                  ZERO=.FALSE.
                    ENDDO !nk
                  ELSE !not HERMITE
                    nk=NKNI(NO_DERIV)
                    IF(DABS(XP(nk,nv,nj_update,np)).GT.1.0D-8)
     '                ZERO=.FALSE.
                  ENDIF !HERMITE
                ENDDO !nv
              ENDIF !NJTOT
            ENDDO !nonode(np)
          ENDIF !linear

          IF(ZERO.OR.HERMITE) THEN
            IF(ZERO.AND..NOT.LINEAR) THEN
              WRITE(OP_STRING,'('' All derivs are currently zero'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
C news MPN 2Feb96
            IF(HERMITE) THEN
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                nb_old=NBJ(nj_update,ne)
                CALL ASSERT(NIT(nb_old).EQ.NIT(nb_new),'>>New and old'
     '            //' bases for specified nj must have same # ni.s in'
     '            //'all elems',ERROR,*9999)
              ENDDO !ne
              nb_old=NBJ(nj_update,1) !for first element (update later)
              NO_DERIV=0
              nb_cubics(0)=0
              DO ni=1,NIT(nb_new)
                IF(IBT(1,ni,nb_new).EQ.2.AND.IBT(2,ni,nb_new).EQ.1) THEN
C                 new basis is cubic Hermite in current ni dirn
                  IF(IBT(1,ni,nb_old).EQ.1.AND.IBT(2,ni,nb_old).EQ.1)
     '              THEN
C                   old basis is linear in current ni dirn
                    IF(NO_DERIV.EQ.0) THEN
                      NO_DERIV=ni
                    ELSE
                      ERROR='>>Change only one Xi dirn interp at a time'
                      GO TO 9999
                    ENDIF
                  ELSE IF(IBT(1,ni,nb_old).EQ.2.AND.
     '                IBT(2,ni,nb_old).EQ.1) THEN
C                   old basis is cubic Hermite in current ni dirn
C                   count the number of cubic Hermite dirns for
C                   basis nb_old
                    nb_cubics(0)=nb_cubics(0)+1
                    nb_cubics(nb_cubics(0))=ni
                  ELSE !old basis is not linear or cubic
                    ERROR='>>Can only update linear to cubic Hermite'
                    GO TO 9999
                  ENDIF !IBT for old basis
                ENDIF !IBT for new basis
              ENDDO !ni
              CALL ASSERT(NO_DERIV.NE.0,'>>No derivatives to update!',
     '          ERROR,*9999)
C             Increase NKJ(nj_update,np) for every node
C             and add new derivatives onto end of list
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
C GMH 8/1/97 Update cmgui link
                CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                DO nv=1,NVJP(nj_update,np)
C                 Initialise all new derivs
                  DO nk=NKJ(nj_update,np)+1,NUM_NKS(nb_cubics(0)+1)
                    XP(nk,nv,nj_update,np)=0.0d0
                  ENDDO !nk
                  IF(.NOT.ZERO) THEN
C                   move old derivs into new deriv positions
                    IF(nb_cubics(0).EQ.1) THEN
C                     old basis is cubic Hermite in only 1 Xi dirn
                      IF(NO_DERIV.LT.nb_cubics(1)) THEN
C                       deriv wrt s(1) -> deriv wrt s(2)
                        XP(3,nv,nj_update,np)=XP(2,nv,nj_update,np)
                        XP(2,nv,nj_update,np)=0.0d0
                      ENDIF
                    ELSE IF(nb_cubics(0).EQ.2) THEN
C                     old basis is cubic Hermite in 2 Xi dirns
                      IF(NO_DERIV.LT.nb_cubics(1)) THEN
C                       deriv wrt s(1)+s(2) -> deriv wrt s(2)+s(3)
                        XP(7,nv,nj_update,np)=XP(4,nv,nj_update,np)
                        XP(4,nv,nj_update,np)=0.0d0
C                       deriv wrt s(2) -> deriv wrt s(3)
                        XP(5,nv,nj_update,np)=XP(3,nv,nj_update,np)
                        XP(3,nv,nj_update,np)=0.0d0
C                       deriv wrt s(1) -> deriv wrt s(2)
                        XP(3,nv,nj_update,np)=XP(2,nv,nj_update,np)
                        XP(2,nv,nj_update,np)=0.0d0
                      ELSE IF(NO_DERIV.LT.nb_cubics(2)) THEN
C                       deriv wrt s(1)+s(2) -> deriv wrt s(3)+s(1)
                        XP(6,nv,nj_update,np)=XP(4,nv,nj_update,np)
                        XP(4,nv,nj_update,np)=0.0d0
C                       deriv wrt s(2) -> deriv wrt s(3)
                        XP(5,nv,nj_update,np)=XP(3,nv,nj_update,np)
                        XP(3,nv,nj_update,np)=0.0d0
                      ENDIF
                    ELSE
                      ERROR='>>should never get into here!!'
                      GO TO 9999
                    ENDIF
                  ENDIF !.NOT.ZERO
                ENDDO !nv
                NKJ(nj_update,np)=NUM_NKS(nb_cubics(0)+1)
              ENDDO !nonode (np)
C ***         Update NBJ(nj_update,ne), NKE(nk,nn,nb,ne),
C ***         NPNE(nn,nb,ne) and NVJE(nn,nb,nj,ne) arrays
C ***         for nb_new from nb_old
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                NBJ(nj_update,ne)=nb_new
                DO nn=1,NNT(nb_new)
                  NPNE(nn,nb_new,ne)=NPNE(nn,nb_old,ne)
                  NVJE(nn,nb_new,nj_update,ne)=
     '              NVJE(nn,nb_old,nj_update,ne)
                ENDDO !nn
              ENDDO !noelem (ne)
C KAT 8Feb01: Derivs and mapping are only changed when bases are
C             changed. i.e. HERMITE.
C            ENDIF !HERMITE
C ***         Update NKE(nk,nn,nj,ne)
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                DO njl=1,3 !geometry/fibres/field
                  DO njj=1,NJ_LOC(njl,0,nr)
                    nj=NJ_LOC(njl,njj,nr)
                    nb=NBJ(nj,ne)
                    DO nn=1,NNT(nb)
                      DO nk=1,NKT(nn,nb)
                        NKJE(nk,nn,nj,ne)=nk
                      ENDDO !nk
                    ENDDO !nn
                  ENDDO !njj
                ENDDO !njl
              ENDDO !noelem (ne)
C             Update NENP array
              CALL CALC_NENP(NBJ,NEELEM,NENP,NPNE,ERROR,*9999)
C ***         23-6-92.AJP. Added to correct NKE if required
              CALL DERIV_INFO(IDO,NBJ,NEELEM,NENP,NKJE,NKJ,
     '          NPNE,NPNODE,XP,ERROR,*9999)
C             set up basis fn, #derivs and mapping arrays for geometry
C             also update NKJ(nj,np) array (Must be done after call
C             deriv_info)
              CALL GLOBALJ(NBJ,NEELEM,NKJ,NPLIST,NPNE,NPNODE,
     '          ERROR,*9999)
C            IF(HERMITE) THEN
C KAT 26Aug98: I have moved the call to LINCAL to after the nl
C              loop so that the lengths DL are not corrupted by
C              calculations based on the new bases which don't have the
C              required derivatives set yet.
C ***         Calculate new derivatives
              DO nl=1,NLT
                ni=NPL(1,0,nl)
                IF(ni.EQ.NO_DERIV) THEN !update specified derivative
                  DS=DL(3,nl) !length of current line
                  NP1=NPL(2,1,nl)
                  NP2=NPL(3,1,nl)
                  nv1=NVJL(1,nj_update,nl)
                  nv2=NVJL(2,nj_update,nl)
C GMH 8/1/97 Update cmgui link
                  CALL NODE_CHANGE(NP1,.FALSE.,ERROR,*9999)
C GMH 8/1/97 Update cmgui link
                  CALL NODE_CHANGE(NP2,.FALSE.,ERROR,*9999)
                  DO nn=1,2
                    nl_temp=NPL(1+nn,0,nl)
                    IF(nl_temp.EQ.0) THEN
                      DERIV_FACTOR(nn)=1.0d0
                    ELSE IF(NVJL(3-nn,nj_update,nl_temp).
     '                EQ.NVJL(nn,nj_update,nl)) THEN
C                     There is another line using the same version.
                      IF(NPL(4-nn,0,nl_temp).EQ.nl) THEN
C                       There is another line in the -Xi dirn so make the
C                       new deriv the average of the two calc'ed derivs.
                        DERIV_FACTOR(nn)=0.5d0
                      ELSE
C                       The other line does not know it is connected to
C                       this line.  This may happen if a line connects to
C                       more than one line.  Use calculation from other
C                       lines for derivative at this node.
                        DERIV_FACTOR(nn)=0.0d0
                      ENDIF
                    ELSE
                      DERIV_FACTOR(nn)=1.0d0
                    ENDIF
                  ENDDO !nn
C                 Find global nk numbers for derivatives to be updated
                  IF(nb_cubics(0).EQ.0) THEN
C                   old basis is linear in all Xi directions
                    IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                    call mp_setlock()
                      WRITE(OP_STRING,'('' Old basis is linear '
     '                  //'in all Xi dirns'')')
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                    call mp_unsetlock()
                    ENDIF
                    NKLIST_old(1)=1
                    NKLIST_new(1)=2
                  ELSE IF(nb_cubics(0).EQ.1) THEN
C                   old basis is cubic Hermite in only 1 Xi dirn
                    IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                    call mp_setlock()
                      WRITE(OP_STRING,'('' Old basis is cubic Hermite '
     '                  //'in only 1 Xi dirn'')')
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                    call mp_unsetlock()
                    ENDIF
                    IF(NIT(nb_new).EQ.2) THEN
                      IF(  IBT(1,1,nb_old).EQ.1.AND.IBT(2,1,nb_old).EQ.1
     '                .AND.IBT(1,2,nb_old).EQ.2.AND.IBT(2,2,nb_old).EQ.1
     '                  ) THEN
C                       old basis is linear/cubic-Hermite in Xi1/2
                        NKLIST_old(1)=1
                        NKLIST_old(2)=3
                        NKLIST_new(1)=2
                        NKLIST_new(2)=4
                      ELSE
C                       old basis is cubic-Hermite/linear in Xi1/2
                        NKLIST_old(1)=1
                        NKLIST_old(2)=2
                        NKLIST_new(1)=3
                        NKLIST_new(2)=4
                      ENDIF
                    ELSE IF(NIT(nb_new).EQ.3) THEN
                      IF(  IBT(1,1,nb_old).EQ.1.AND.IBT(2,1,nb_old).EQ.1
     '                .AND.IBT(1,2,nb_old).EQ.2.AND.IBT(2,2,nb_old).EQ.1
     '                .AND.IBT(1,3,nb_old).EQ.1.AND.IBT(2,3,nb_old).EQ.1
     '                .AND.IBT(1,1,nb_new).EQ.2.AND.IBT(2,1,nb_new).EQ.1
     '                .AND.IBT(1,2,nb_new).EQ.2.AND.IBT(2,2,nb_new).EQ.1
     '                .AND.IBT(1,3,nb_new).EQ.1.AND.IBT(2,3,nb_new).EQ.1
     '                 .OR.IBT(1,1,nb_old).EQ.1.AND.IBT(2,1,nb_old).EQ.1
     '                .AND.IBT(1,2,nb_old).EQ.1.AND.IBT(2,2,nb_old).EQ.1
     '                .AND.IBT(1,3,nb_old).EQ.2.AND.IBT(2,3,nb_old).EQ.1
     '               .AND.(IBT(1,1,nb_new).EQ.2.AND.IBT(2,1,nb_new).EQ.1
     '                .AND.IBT(1,2,nb_new).EQ.1.AND.IBT(2,2,nb_new).EQ.1
     '                .AND.IBT(1,3,nb_new).EQ.2.AND.IBT(2,3,nb_new).EQ.1
     '                 .OR.IBT(1,1,nb_new).EQ.1.AND.IBT(2,1,nb_new).EQ.1
     '                .AND.IBT(1,2,nb_new).EQ.2.AND.IBT(2,2,nb_new).EQ.1
     '                .AND.IBT(1,3,nb_new).EQ.2.AND.IBT(2,3,nb_new).EQ.1
     '                  )) THEN
C                       old basis is lin./cubic-Hermite/lin. in Xi1/2/3
C                       and new basis is bicubic-Hermite/lin. in Xi1/2/3
C                       OR
C                       old basis is bilinear/cubic-Hermite in Xi1/2/3
C                       and new basis is cub-Her/lin/cub-Her in Xi1/2/3
C                       OR
C                       old basis is bilinear/cubic-Hermite in Xi1/2/3
C                       and new basis is lin./bicubic-Hermite in Xi1/2/3
                        NKLIST_old(1)=1
                        NKLIST_old(2)=3
                        NKLIST_new(1)=2
                        NKLIST_new(2)=4
                      ELSE
C                       old basis is cubic-Hermite/bilinear in Xi1/2/3
C                       and new basis is cub-Her/lin/cub-Her in Xi1/2/3
C                       OR
C                       old basis is cubic-Hermite/bilinear in Xi1/2/3
C                       and new basis is bicubic-Hermite/lin. in Xi1/2/3
C                       OR
C                       old basis is lin./cubic-Hermite/lin. in Xi1/2/3
C                       and new basis is lin./bicubic-Hermite in Xi1/2/3
                        NKLIST_old(1)=1
                        NKLIST_old(2)=2
                        NKLIST_new(1)=3
                        NKLIST_new(2)=4
                      ENDIF
                    ENDIF
                  ELSE IF(nb_cubics(0).EQ.2) THEN
C                   old basis is cubic Hermite in 2 Xi dirns
C                   so new basis is tri-cubic-Hermite.
                    IF(IBT(1,1,nb_old).EQ.1.AND.IBT(2,1,nb_old).EQ.1
     '                .AND.IBT(1,2,nb_old).EQ.2.AND.IBT(2,2,nb_old).EQ.1
     '                .AND.IBT(1,3,nb_old).EQ.2.AND.IBT(2,3,nb_old).EQ.1
     '                ) THEN
C                     old basis is linear/bicubic-Hermite in Xi1/2/3
                      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                      call mp_setlock()
                        WRITE(OP_STRING,'('' Old basis is '
     '                    //'linear/bicubic Hermite'')')
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                      call mp_unsetlock()
                      ENDIF
                      NKLIST_old(1)=1
                      NKLIST_old(2)=3
                      NKLIST_old(3)=5
                      NKLIST_old(4)=7
                      NKLIST_new(1)=2
                      NKLIST_new(2)=4
                      NKLIST_new(3)=6
                      NKLIST_new(4)=8
                    ELSEIF(IBT(1,1,nb_old).EQ.2.AND.IBT(2,1,nb_old).EQ.1
     '                .AND.IBT(1,2,nb_old).EQ.1.AND.IBT(2,2,nb_old).EQ.1
     '                .AND.IBT(1,3,nb_old).EQ.2.AND.IBT(2,3,nb_old).EQ.1
     '                ) THEN
C                     old basis is cub-Herm/linear/cub-Herm in Xi1/2/3
                      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                      call mp_setlock()
                        WRITE(OP_STRING,'('' Old basis is '
     '                    //'cub-Herm/linear/cub-Herm'')')
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                      call mp_unsetlock()
                      ENDIF
                      NKLIST_old(1)=1
                      NKLIST_old(2)=2
                      NKLIST_old(3)=5
                      NKLIST_old(4)=6
                      NKLIST_new(1)=3
                      NKLIST_new(2)=4
                      NKLIST_new(3)=7
                      NKLIST_new(4)=8
                    ELSEIF(IBT(1,1,nb_old).EQ.2.AND.IBT(2,1,nb_old).EQ.1
     '                .AND.IBT(1,2,nb_old).EQ.2.AND.IBT(2,2,nb_old).EQ.1
     '                .AND.IBT(1,3,nb_old).EQ.1.AND.IBT(2,3,nb_old).EQ.1
     '                ) THEN
C                     old basis is bicubic-Hermite/linear in Xi1/2/3
                      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                      call mp_setlock()
                        WRITE(OP_STRING,'('' Old basis is '
     '                    //'bicubic-Hermite/linear'')')
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                      call mp_unsetlock()
                      ENDIF
                      NKLIST_old(1)=1
                      NKLIST_old(2)=2
                      NKLIST_old(3)=3
                      NKLIST_old(4)=4
                      NKLIST_new(1)=5
                      NKLIST_new(2)=6
                      NKLIST_new(3)=7
                      NKLIST_new(4)=8
                    ENDIF !basis type
                  ENDIF !nb_cubics
                  DO no_nklist=1,NUM_NKS(nb_cubics(0))
                    nk_old=NKLIST_old(no_nklist)
C                   Calculate difference in coord/deriv  nj_update
C                   and correct for angular coords that cross
C                   the 0<->2pi discontinuity
                    DX_update=XP(nk_old,nv2,nj_update,NP2)-
     '                XP(nk_old,nv1,nj_update,NP1)
                    IF(ITYP10(nr).GT.1.AND.nk_old.EQ.1) THEN
                      IF(nj_update.EQ.2) THEN
                        IF((XP(1,nv2,nj_update,NP2).GE.3.d0*PI/2.d0)
     '                    .AND.(XP(1,nv1,nj_update,NP1).LT.PI/2.d0)
     '                    .AND.DX_update.GT.PI) THEN
                          DX_update=DX_update-2.d0*PI !acr. 0->2pi discont'y
                        ELSE IF((XP(1,nv2,nj_update,NP2).LT.PI/2.d0)
     '                      .AND.(XP(1,nv1,nj_update,NP1).GT.1.5d0*PI)
     '                      .AND.DX_update.LT.-PI) THEN
                          DX_update=DX_update+2.d0*PI !acr. 2pi->0 discont'y
                        ENDIF
                      ENDIF
                      IF(ITYP10(nr).GT.2.AND.nj_update.EQ.3) THEN
                        IF((XP(1,nv2,nj_update,NP2).GT.1.5d0*PI)
     '                    .AND.(XP(1,nv1,nj_update,NP1).LT.PI/2.d0)
     '                    .AND.DX_update.GT.PI) THEN
                          DX_update=DX_update-2.d0*PI !acr. 0->2pi discont'y
                        ELSE IF((XP(1,nv2,nj_update,NP2).LT.PI/2.d0)
     '                      .AND.(XP(1,nv1,nj_update,NP1).GT.1.5d0*PI)
     '                      .AND.DX_update.LT.-PI) THEN
                          DX_update=DX_update+2.d0*PI !acr. 2pi->0 discont'y
                        ENDIF
                      ENDIF
                    ENDIF
                    nk_new=NKLIST_new(no_nklist)
                    IF(DABS(DS).GT.ZERO_TOLERANCE) THEN !non-zero length
                      XP(nk_new,nv1,nj_update,NP1)=
     '                  XP(nk_new,nv1,nj_update,NP1)+
     '                  DERIV_FACTOR(1)*DX_update/DS
                      XP(nk_new,nv2,nj_update,NP2)=
     '                  XP(nk_new,nv2,nj_update,NP2)+
     '                  DERIV_FACTOR(2)*DX_update/DS
                    ELSE !zero line length
                      WRITE(OP_STRING,'('' >>WARNING: Zero line'
     '                  //' length -> setting deriv to zero'')')
                      CALL WRITES(IOER,OP_STRING,ERROR,*9999)
                    ENDIF
                  ENDDO !no_nklist
C!!! MPN 3Feb96: Above code could be updated later:
C    Global nks should be picked up from NPL (first
C    derivs in each direction) and NPF (for cross derivs).
C    But NPL doesn't store first derivs correctly and NPF doesn't
C    store derivs at all.
C                  NK1=NPL(4,1,nl)
C                  NK2=NPL(5,1,nl)
C!!! MPN 2Feb96: This should be storage for NPL but it isn't!
C                  IF(nj_update.EQ.1) THEN
C                    NK1=NPL(4,1,nl)
C                    NK2=NPL(5,1,nl)
C                  ELSE IF(nj_update.EQ.2) THEN
C                    NK1=NPL(4,2,nl)
C                    NK2=NPL(5,2,nl)
C                  ELSE IF(nj_update.EQ.3) THEN
C                    NK1=NPL(4,3,nl)
C                    NK2=NPL(5,3,nl)
C                  ENDIF
                ENDIF !ni.EQ.NO_DERIV
              ENDDO !nl
C             Update line & face info
              CALL LINCAL(IBT,IDO,INP,0,NBJ,NEELEM,NEL,NENP,
     '          NKJE,NLL,NLLINE,NNL,NPL,NPNE,NPNODE,NRE,NVJE,NVJL,
     '          DL,SE,XP,ERROR,*9999)
              CALL FACCAL(IBT,IDO,INP,NBJ,NBJF,NEELEM,NFF,NFFACE,
     '          NKJE,NKEF,NLF,NNF,NNL,NPF,NPL,NPNE,NPNF,
     '          NRE,NVJE,NVJF,DF,PG,RG,SE,SF,WG,XA,XE,
     '          XG,XP,ERROR,*9999)
C newe
            ELSE !not HERMITE
C             calculate derivs from straight line approx
              WRITE(OP_STRING,'('' Calculate derivs from '
     '          //'straight line approx'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C ***         Update derivatives
C GMH 23/4/96 Loop for subset of nodes
              DO nolist=1,NPLIST(0)
                np=NPLIST(nolist)

C GR    7/2/07 Added individual versions option so that only the
C       version relevant to the line is considered instead of
C       taking the average of all versions.
                IF(.NOT.INDIVIDUAL_VERSIONS) THEN

C               Find all lines that go through the node
                DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                  nj=NJ_LOC(NJL_GEOM,njj,nr)
                  DX(nj)=0.0d0
                ENDDO !nj
                NL_COUNT=0
                DO nl=1,NLT
                  ni=NPL(1,0,nl)
                  NP1=NPL(2,1,nl)
                  NP2=NPL(3,1,nl)
C new PJH 15Sep92
                  IF((NP1.EQ.np.OR.NP2.EQ.np).AND.
     '              ni.EQ.NO_DERIV) THEN !update specified derivative
                    NL_COUNT=NL_COUNT+1
                    DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                      nj=NJ_LOC(NJL_GEOM,njj,nr)
C!!! NOTE: this will need updating for different versions
                      DX(nj)=DX(nj)+XP(1,1,nj,NP2)-XP(1,1,nj,NP1)
C!!!                changed GT to GE here AAY 29Aug1995
                      IF(ITYP10(nr).GT.1.AND.nj.EQ.2.AND.
     '                  (XP(1,1,nj,NP2).GE.3.d0*PI/2.d0).AND.
     '                  (XP(1,1,nj,NP1).LT.PI/2.d0).AND.
     '                  DX(nj).GT.PI) THEN !across 0->2pi discontinuity
                        DX(nj)=DX(nj)-2.d0*PI
                      ELSE IF(ITYP10(nr).GT.1.AND.nj.EQ.2.AND.
     '                    (XP(1,1,nj,NP2).LT.PI/2.d0).AND.
     '                    (XP(1,1,nj,NP1).GT.3.d0*PI/2.d0).AND.
     '                    DX(nj).LT.-PI)THEN !across 2pi->0 discont'y
                        DX(nj)=DX(nj)+2.d0*PI
                      ELSE IF(ITYP10(nr).GT.2.AND.nj.EQ.3.AND.
     '                    (XP(1,1,nj,NP2).GT.3.d0*PI/2.d0).AND.
     '                    (XP(1,1,nj,NP1).LT.PI/2.d0).AND.
     '                    DX(nj).GT.PI)THEN !across 0->2pi discontinuity
                        DX(nj)=DX(nj)-2.d0*PI
                      ELSE IF(ITYP10(nr).GT.1.AND.nj.EQ.3.AND.
     '                    (XP(1,1,nj,NP2).LT.PI/2.d0).AND.
     '                    (XP(1,1,nj,NP1).GT.3.d0*PI/2.d0).AND.
     '                    DX(nj).LT.-PI)THEN !across 2pi->0 discont'y
                        DX(nj)=DX(nj)+2.d0*PI
                      ENDIF
                    ENDDO !njj (nj)
                  ENDIF !node,deriv
                ENDDO !nl
                SUM=0.d0
                DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                  nj=NJ_LOC(NJL_GEOM,njj,nr)
                  SUM=SUM+DX(nj)*DX(nj)
                ENDDO !nj
                IF(ARCLENGTH) THEN
                  IF(DABS(SUM).GT.RDELTA) THEN
                    DS=DSQRT(SUM)
C!!! KAT 6Dec00: This only seems designed for arc-length scale-factors.
C CPB 1/7/96 Normalise DX
                    DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                      nj=NJ_LOC(NJL_GEOM,njj,nr)
                      DX(nj)=DX(nj)/DS
                    ENDDO !nj
                  ENDIF
                ELSE !xi derivs
                  IF(NL_COUNT.GT.0) THEN
C                   Take an average of values from the lines
                    DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                      nj=NJ_LOC(NJL_GEOM,njj,nr)
                      DX(nj)=DX(nj)/NL_COUNT
                    ENDDO !nj
                  ENDIF !NL_COUNT>0
                ENDIF !ARCLENGTH

C KAT 6Dec00:   We probably could try to obtain nk from NPL (but
C               according to MPN 3Feb96 these are wrong) or from
C               inverting NUNK but we just use the simple NKNI.  This
C               only needs to be done once in this subroutine but if it
C               were done using NUNK it would be done here.
                nk=NKNI(NO_DERIV)
                IF(NJTOT.EQ.NJT) THEN
                  DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                    nj=NJ_LOC(NJL_GEOM,njj,nr)
C CS 28/12/2000 adding loop over versions
                    DO nv=1,NVJP(nj,np)
                      IF (nv.GE.2) THEN
                        WRITE(OP_STRING,'('' >>WARNING: Updating'
     '                  //' higher versions too'')')
                        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
                      ENDIF
                      XP(nk,nv,nj,np)=DX(nj)
                    ENDDO
                  ENDDO
                ELSE IF(NJTOT.EQ.1)THEN
                  DO nv=1,NVJP(nj,np)
                    IF (nv.GE.2) THEN
                      WRITE(OP_STRING,'('' >>WARNING: Updating'
     '                  //' higher versions too'')')
                      CALL WRITES(IOER,OP_STRING,ERROR,*9999)
                    ENDIF
                    XP(nk,1,nj_update,np)=DX(nj_update)
                  ENDDO
                ENDIF !NJTOT

C KAT 6Dec00:   We have calculated an average derivative at np.  We don't
C               want to store this at NP1 and NP2 for every line that
C               uses this node.
C                DO nl=1,NLT
C                  ni=NPL(1,0,nl)
C                  NP1=NPL(2,1,nl)
C                  NP2=NPL(3,1,nl)
CC GMH 8/1/97 Update cmgui link
C                  CALL NODE_CHANGE(NP1,.FALSE.,ERROR,*9999)
CC GMH 8/1/97 Update cmgui link
C                  CALL NODE_CHANGE(NP2,.FALSE.,ERROR,*9999)
C                  IF(((NP1.EQ.np).OR.(NP2.EQ.np)).AND.
C     '              (ni.EQ.NO_DERIV)) THEN !update specified derivative
CC                 find interpolation type and scaling type for line
CC KAT 6Dec00: This stuff either doesn't work or isn't used.
CC                    AVEARCLENGTH=.FALSE.
CC                    ARCLENGTH=.FALSE.
C                    NK1=0
C                    NK2=0
C                    IF(NJTOT.EQ.1) THEN !one geometric variable only
C                      IF(NPL(1,nj_update,nl).EQ.4) THEN !cubic Herm line
CC KAT 6Dec00: This stuff either doesn't work or isn't used.
CC                        DO nb=1,NBFT
CC                          IF(IBT(1,ni,nb).EQ.2.AND.
CC     '                      IBT(2,ni,nb).EQ.1) THEN
CC                            AVEARCLENGTH=(NBI(nb).EQ.6.OR.NBI(nb).EQ.7)
CC                            ARCLENGTH=(NBI(nb).EQ.5)
CC                          ENDIF
CC                        ENDDO !nb
C                        NK1=NPL(4,1,nl)
C                        NK2=NPL(5,1,nl)
C                      ENDIF
C                    ELSE !loop thru all geometric variables
C                      DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
C                        nj=NJ_LOC(NJL_GEOM,njj,nr)
C                        IF(NPL(1,nj,nl).EQ.4)THEN !cubic Hermite line
CC KAT 6Dec00: This stuff either doesn't work or isn't used.
CC                          DO nb=1,NBFT
CC                            IF(IBT(1,ni,nb).EQ.2.AND.IBT(2,ni,nb).EQ.1)
CC     '                        THEN
CC                              !cubic Hermite
CC                              AVEARCLENGTH=
CC     '                          (NBI(nb).EQ.6.OR.NBI(nb).EQ.7)
CCC!!! KAT 6Dec00: What about ARCLENGTH?
CC                            ENDIF
CC                          ENDDO !nb
C                          NK1=NPL(4,1,nl)
C                          NK2=NPL(5,1,nl)
C                        ENDIF
C                      ENDDO !njj (nj)
C                    ENDIF !NJTOT

CC KAT 6Dec00: This stuff either doesn't work or isn't used.
CC                    IF(AVEARCLENGTH.OR.ARCLENGTH)
CC     '                THEN !arc length scaling factors
CC                      IF(DL(3,nl).EQ.0.d0) THEN !no lengths calculated
CCC CPB 1/7/96 ALREADY NORMALISED
CCC                        DS=DSQRT(SUM) !is strt line arc-length btw nodes
CCC                     update scaleing factors
CCC!!! KAT 6Dec00: DS is not set
CC                        DL(1,nl)=DS !is deriv scaling factor at LH end
CC                        DL(2,nl)=DS !is deriv scaling factor at RH end
CC                        DL(3,nl)=DS !is arc-length
CC                      ELSE
CC                        DS=DL(3,nl) !is arc-length between nodes
CC                      ENDIF
CC                    ELSE !assume angle change or whatever in DL(1,nl)
CC                      DS=DL(1,nl)
CC                    ENDIF
C                    IF(NJTOT.EQ.NJT.AND.NK1.GT.0.AND.NK2.GT.0)THEN
C                      DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
C                        nj=NJ_LOC(NJL_GEOM,njj,nr)
CC                        XP(NK1,1,nj,NP1)=DX(nj)/DS
CC                        XP(NK2,1,nj,NP2)=DX(nj)/DS
C                        XP(NK1,1,nj,NP1)=DX(nj)
C                        XP(NK2,1,nj,NP2)=DX(nj)
C                      ENDDO
C                    ELSE IF(NJTOT.EQ.1.AND.NK1.GT.0.AND.NK2.GT.0)THEN
CC                      XP(NK1,1,nj_update,NP1)=DX(nj_update)/DS
CC                      XP(NK2,1,nj_update,NP2)=DX(nj_update)/DS
C                      XP(NK1,1,nj_update,NP1)=DX(nj_update)
C                      XP(NK2,1,nj_update,NP2)=DX(nj_update)
C                    ENDIF !NJTOT
C                  ENDIF ! ni.EQ.NO_DERIV
C                ENDDO !nl
C GR    7/2/07 Added individual versions option so that only the
C       version relevant to the line is considered instead of
C       taking the average of all versions.
                ELSE
                NOLINE(0) = 0
                has_line = 0
                DO nv=1,NVJP(nj,np)

C                 Initialise temp memory
                  DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                    nj=NJ_LOC(NJL_GEOM,njj,nr)
                    DX(nj)=0.0d0
                  ENDDO !njj
C                 Find all lines that go through the node and version
                  NL_COUNT=0
                  DO nl=1,NLT ! loop over all lines
                    ni=NPL(1,0,nl) 
                    NP1=NPL(2,1,nl)
                    NP2=NPL(3,1,nl)
                    ! find the version of each node that is
                    ! connected to the line
                    nv1=NVJL(1,nj,nl)
                    nv2=NVJL(2,nj,nl)

                    IF((((NP1.EQ.np).AND.(nv1.EQ.nv)).OR.
     &                  ((NP2.EQ.np).AND.(nv2.EQ.nv)))
     &                  .AND.ni.EQ.NO_DERIV) THEN
                      NL_COUNT=NL_COUNT+1
                      DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                        nj=NJ_LOC(NJL_GEOM,njj,nr)
C!!! NOTE: this will need updating for different versions
                      !DX(nj)=DX(nj)+XP(1,1,nj,NP2)-XP(1,1,nj,NP1)

                        ! update the delta using only the relevant version
                        DX(nj)=DX(nj)+XP(1,nv2,nj,NP2)-XP(1,nv1,nj,NP1)

                        IF(DOP) THEN
                        WRITE(OP_STRING,'(''nl '',I3,'' np '',I3,'
     &                  //''' nv '',I1,I3,'':'',I1,'
     &                  //'''<--->'',I3,'':'',I1)')
     &                  nl,np,nv,NP1,nv1,NP2,nv2
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                        WRITE(OP_STRING,'(''DX('',I1,'')='',F8.3)')
     &                    nj,DX(nj)
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                        ENDIF

C!!!                changed GT to GE here AAY 29Aug1995
                        IF(ITYP10(nr).GT.1.AND.nj.EQ.2.AND.
     '                      (XP(1,1,nj,NP2).GE.3.d0*PI/2.d0).AND.
     '                      (XP(1,1,nj,NP1).LT.PI/2.d0).AND.
     '                      DX(nj).GT.PI) THEN !across 0->2pi discontinuity
                          DX(nj)=DX(nj)-2.d0*PI
                        ELSE IF(ITYP10(nr).GT.1.AND.nj.EQ.2.AND.
     '                      (XP(1,1,nj,NP2).LT.PI/2.d0).AND.
     '                      (XP(1,1,nj,NP1).GT.3.d0*PI/2.d0).AND.
     '                      DX(nj).LT.-PI) THEN !across 2pi->0 discont'y
                          DX(nj)=DX(nj)+2.d0*PI
                        ELSE IF(ITYP10(nr).GT.2.AND.nj.EQ.3.AND.
     '                      (XP(1,1,nj,NP2).GT.3.d0*PI/2.d0).AND.
     '                      (XP(1,1,nj,NP1).LT.PI/2.d0).AND.
     '                      DX(nj).GT.PI)THEN !across 0->2pi discontinuity
                          DX(nj)=DX(nj)-2.d0*PI
                        ELSE IF(ITYP10(nr).GT.1.AND.nj.EQ.3.AND.
     '                        (XP(1,1,nj,NP2).LT.PI/2.d0).AND.
     '                        (XP(1,1,nj,NP1).GT.3.d0*PI/2.d0).AND.
     '                        DX(nj).LT.-PI)THEN !across 2pi->0 discont'y
                          DX(nj)=DX(nj)+2.d0*PI
                        ENDIF
                      ENDDO !njj (nj)
                    ENDIF !node,version
                  ENDDO !nl
                  SUM=0.d0
                  DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                    nj=NJ_LOC(NJL_GEOM,njj,nr)
                    SUM=SUM+DX(nj)*DX(nj)
                  ENDDO !nj
                  IF(ARCLENGTH) THEN
                    IF(DABS(SUM).GT.RDELTA) THEN
                      DS=DSQRT(SUM)
C!!! KAT 6Dec00: This only seems designed for arc-length scale-factors.
C CPB 1/7/96 Normalise DX
                      DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                        nj=NJ_LOC(NJL_GEOM,njj,nr)
                        DX(nj)=DX(nj)/DS
                      ENDDO !nj
                    ENDIF
                  ELSE !xi derivs
                    IF(NL_COUNT.GT.0) THEN
C                     Take an average of values from the lines
                      DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                        nj=NJ_LOC(NJL_GEOM,njj,nr)
                        DX(nj)=DX(nj)/NL_COUNT
                      ENDDO !nj
                    ENDIF !NL_COUNT>0
                  ENDIF !ARCLENGTH

                  ! if there are no lines attached to this version
                  ! save the version number so we can fill in the data later
                  IF(NL_COUNT.EQ.0) THEN
                    NOLINE(0) = NOLINE(0)+1
                    NOLINE(NOLINE(0)) = nv
                  ELSE
                    has_line = nv
                  ENDIF

C KAT 6Dec00:   We probably could try to obtain nk from NPL (but
C               according to MPN 3Feb96 these are wrong) or from
C               inverting NUNK but we just use the simple NKNI.  This
C               only needs to be done once in this subroutine but if it
C               were done using NUNK it would be done here.
                  nk=NKNI(NO_DERIV)
                  IF(NJTOT.EQ.NJT) THEN
                    DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                      nj=NJ_LOC(NJL_GEOM,njj,nr)
C CS 28/12/2000 adding loop over versions
C GR 15/7/06 loop no longer needed because we only want to update
C            the relevant version, stored in nv.
!                     DO nv=1,NVJP(nj,np)
!                       IF (nv.GE.2) THEN
!                         WRITE(OP_STRING,'('' >>WARNING: Updating'
!      '                  //' higher versions too'')')
!                         CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
!                       ENDIF
!                       XP(nk,nv,nj,np)=DX(nj)
!                     ENDDO
                      IF(DOP) THEN
                      WRITE(OP_STRING,'(''updating XP np '',I3,'' nv '
     &                  //''',I3,'' nk '',I1,'
     &                  //''' DX('',I1,'')='',F8.3)')
     &                  np,nv,nk,nj,DX(nj)
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      ENDIF
                      XP(nk,nv,nj,np)=DX(nj)
                    ENDDO
                  ELSE IF(NJTOT.EQ.1)THEN
!                   DO nv=1,NVJP(nj,np)
!                     IF (nv.GE.2) THEN
!                       WRITE(OP_STRING,'('' >>WARNING: Updating'
!      '                  //' higher versions too'')')
!                       CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
!                     ENDIF
!                     XP(nk,1,nj_update,np)=DX(nj_update)
!                   ENDDO
                    XP(nk,nv,nj_update,np)=DX(nj_update)
                  ENDIF !NJTOT
                ENDDO !nv
                ! If any versions had no lines attached then copy the calculated
                ! value from the last version that had a value. This could happen
                ! in certain types of collapsed elements. A better way to do this
                ! would be to detect the collapsed element and copy the correct
                ! versions, but this seems to work for common scenarios.
                DO nv=1,NOLINE(0)
                  nk=NKNI(NO_DERIV)
                  IF(NJTOT.EQ.NJT) THEN
                    DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                      nj=NJ_LOC(NJL_GEOM,njj,nr)
                    IF(DOP) THEN
                    WRITE(OP_STRING,'(''updating vers w/out line np '
     &                //''',I3,'' nv '',I3,'' nk '',I1,'
     &                //''' DX('',I1,'') to nv='',I3)')
     &                np,NOLINE(nv),nk,nj,has_line
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                    XP(nk,NOLINE(nv),nj,np) = XP(nk,has_line,nj,np)
                  ENDDO
                  ELSE IF(NJTOT.EQ.1)THEN
!                   DO nv=1,NVJP(nj,np)
!                     IF (nv.GE.2) THEN
!                       WRITE(OP_STRING,'('' >>WARNING: Updating'
!      '                  //' higher versions too'')')
!                       CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
!                     ENDIF
!                     XP(nk,1,nj_update,np)=DX(nj_update)
!                   ENDDO
                    IF(DOP) THEN
                    WRITE(OP_STRING,'(''updating vers w/out line np '
     &                //''',I3,'' nv '',I3,'' nk '',I1,'
     &                //''' DX('',I1,'') to nv='',I3)')
     &                np,NOLINE(nv),nk,nj_update,has_line
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                    XP(nk,NOLINE(nv),nj_update,np)
     &                 = XP(nk,has_line,nj_update,np)
                  ENDIF !NJTOT
                ENDDO

                ENDIF !INDIVIDUAL_VERSIONS

              ENDDO !np
              DO nb=1,NBFT
                CALL DLSE(IBT,IDO,nb,NEELEM,NLL,NNL,NPL,
     '            DL,SE,ERROR,*9999)
              ENDDO !nb
            ENDIF !HERMITE

          ELSE IF(.NOT.ZERO) THEN !update derivs to be wrt arc-length
            WRITE(OP_STRING,'('' Some derivatives are currently '','
     '        //'''nonzero'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' Update derivs to be wrt '
     '        //'arc-length'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
              nj=NJ_LOC(NJL_GEOM,njj,nr)
              DXDS(nj)=0.d0
            ENDDO !njj (nj)
C LKC 8-DEC-97 This does all nodes, NOT the selected (if any) ones
C            DO nonode=1,NPNODE(0,nr)
C              np=NPNODE(nonode,nr)
            DO nonode=1,NPLIST(0)
              np=NPLIST(nonode)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
C             G is metric tensor components for length calculation
              IF(ITYP10(nr).EQ.1) THEN !rectangular cartesian coords
                G(1)=1.d0
                G(2)=1.d0
                G(3)=1.d0
              ELSE IF(ITYP10(nr).EQ.2) THEN !cylindrical polar coords
                G(1)=1.d0
                G(2)=XP(1,1,1,np)*XP(1,1,1,np)
                G(3)=1.d0
              ELSE IF(ITYP10(nr).EQ.4) THEN !prolate spherical coords
                G(1)=(DSINH(XP(1,1,1,np))*DSINH(XP(1,1,1,np))+
     '            DSIN(XP(1,1,2,np))*DSIN(XP(1,1,2,np)))*FOCUS*FOCUS
                G(2)=G(1)
                G(3)=DSINH(XP(1,1,1,np))*DSINH(XP(1,1,1,np))*
     '            DSIN(XP(1,1,2,np))*DSIN(XP(1,1,2,np))*FOCUS*FOCUS
              ELSE !spherical or oblate coords
                ERROR='Dinner isn''t ready yet'
                GOTO 9999
              ENDIF
C             find number of first derivatives
              NKMAX=0
              DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                nj=NJ_LOC(NJL_GEOM,njj,nr)
                IF(NKJ(nj,np).GT.NKMAX) NKMAX=NKJ(nj,np)
              ENDDO !njj (nj)
c             IF(NKMAX.EQ.2) THEN
c               NKTOT=2
c             ELSE IF(NKMAX.EQ.4) THEN
c               NKTOT=3
c             ELSE
c               NKTOT=0
c             ENDIF
c             DO nk=2,NKTOT !loop over 1st derivs
C news        AAY 29Aug95 allow NKMAX=7 as well for tricubic basis
              IF(NKMAX.EQ.2) THEN
                NKTOT=1
              ELSE IF(NKMAX.EQ.4) THEN
                NKTOT=2
              ELSE IF(NKMAX.EQ.7) THEN
                NKTOT=3
              ELSE
                NKTOT=0
              ENDIF
              DO nk1=1,NKTOT !loop over 1st derivs
                nk=NKLIST(nk1)
C newe
c CPB 5/4/93 Change update so that the rescaling is based on the current
C line length not the value obtained from the nodal interpolation
C               Find a line going through the node
                nl=1
                FOUND=.FALSE.
C!!! Assumes nk=2 is xi1 deriv
                DO WHILE(.NOT.FOUND.AND.nl.LT.NLT)
                  IF((np.EQ.NPL(2,1,nl).OR.np.EQ.NPL(3,1,nl)).AND.
     '              ((nk-1).EQ.NPL(1,0,nl))) THEN
                    FOUND=.TRUE.
                  ELSE
                    nl=nl+1
                  ENDIF
                ENDDO !nl while loop
C               If a line has been found use its new length to scale
C               the line
                IF(FOUND) THEN
                  IF(DL(3,nl).EQ.0.d0) THEN
                    DSDS(nk-1)=DL(1,nl)
                  ELSE
                    DSDS(nk-1)=DL(1,nl)/DL(3,nl)
                  ENDIF
C                 dsds is the deriv of old arc. len. wrt new arc. len.
                  IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                  call mp_setlock()
                    WRITE(OP_STRING,'('' np='',I3,'' nl='',I3,'
     '                //''' dsds='',E12.5)') np,nl,DSDS(nk-1)
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                  call mp_unsetlock()
                  ENDIF
                  DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                    nj=NJ_LOC(NJL_GEOM,njj,nr)
                    XP(nk,1,nj,np)=XP(nk,1,nj,np)*DSDS(nk-1)
                  ENDDO !njj (nj)
                  WRITE(OP_STRING,'('' Node '',I3,'' updated '','
     '              //'''using arc length rescaling'')') np
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ELSE !not FOUND
C                 Otherwise used a length scaling interpolated from
C                 nodal variables.
                  SUM=0.d0
                  DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                    nj=NJ_LOC(NJL_GEOM,njj,nr)
                    DXDS(nj)=XP(nk,1,nj,np)
                    SUM=SUM+DXDS(nj)*DXDS(nj)*G(nj)
                  ENDDO
                  SUM=DSQRT(SUM) !length of "tangent" vector
                  IF(SUM.GT.1.0D-5) THEN
                    DSDS(nk-1)=1.d0/SUM
                    DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                      nj=NJ_LOC(NJL_GEOM,njj,nr)
                      XP(nk,1,nj,np)=DXDS(nj)/DSDS(nk-1) !norm'd deriv
                    ENDDO
                  ELSE
C                   derivative of old arc. len. wrt new arc. len.
                    DSDS(nk-1)=1.d0
                  ENDIF
                ENDIF !FOUND
              ENDDO !nk1 (nk)
C             cross derivative
              IF(NKMAX.EQ.4)THEN
                DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                  nj=NJ_LOC(NJL_GEOM,njj,nr)
                  XP(4,1,nj,np)=XP(4,1,nj,np)*DSDS(1)*DSDS(2)
                ENDDO !njj (nj)
              ELSE IF(NKMAX.EQ.7) THEN
C news          AAY 29Aug95 allow NKMAX=7 as well for tricubic basis
                DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                  nj=NJ_LOC(NJL_GEOM,njj,nr)
                  XP(4,1,nj,NP)=XP(4,1,nj,NP)*DSDS(1)*DSDS(2)
                  XP(6,1,nj,NP)=XP(6,1,nj,NP)*DSDS(1)*DSDS(3)
                  XP(7,1,nj,NP)=XP(7,1,nj,NP)*DSDS(2)*DSDS(3)
                  XP(8,1,nj,NP)=XP(8,1,nj,NP)*DSDS(1)*DSDS(2)*DSDS(3)
                ENDDO !njj (nj)
C newe
              ENDIF

              IF(KTYP8.EQ.1.AND.NJ_LOC(NJL_FIEL,0,nr).EQ.NJT) THEN
C               geom fitting
                DO njj=1,NJ_LOC(NJL_FIEL,0,nr) !update field var derivs
                  DO nk=2,NKJ(nj,np)
                    XP(nk,1,NJ_LOC(NJL_FIEL,njj,nr),np)=
     '                XP(nk,1,NJ_LOC(NJL_GEOM,njj,nr),np)
                  ENDDO !nk
                ENDDO !njj (nj)
              ENDIF
            ENDDO !np

C           recalculate arc lengths
            DO nl=1,NLT
              CALL ARCSCA(IDO,0,0,0,NBJ,NEL(0,nl),nl,
     '          NPL(1,0,nl),NPNE,NVJL(1,1,nl),DL,1.0D-6,XP,ERROR,*9999)
            ENDDO !nl
C           store new scaling factors
            DO nb=1,NBFT
              CALL DLSE(IBT,IDO,nb,NEELEM,NLL,NNL,NPL,
     '          DL,SE,ERROR,*9999)
            ENDDO !nb
          ENDIF !not ZERO

        ELSE IF(CONTINUOUS) THEN !ensure continuity at internal nodes
C         Update only internal nodes.
          CALL ASSERT(NJ_LOC(NJL_GEOM,0,nr).EQ.3,
     '      '>>Only implemented for 3D',ERROR,*9999)
          DO njl=1,3
            DO njj=1,NJ_LOC(njl,0,nr)
              nj=NJ_LOC(njl,njj,nr)
C             Initiallize nodes as internal.
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                NPLIST(np)=1
              ENDDO !nonode
C             Check if internal by searching for boundary faces and marking
C             their nodes as external.
              DO no_nf=1,NFFACE(0,nr)
                nf=NFFACE(no_nf,nr)
                IF(NPF(5,nf).EQ.1) THEN !external
                  ne=NPF(6,nf) !adjacent element
                  nf_e=NPF(8,nf)
                  nb_e=NBJ(nj,ne)
                  nb_f=NBJF(nj,nf)
                  DO nn_f=1,NNT(nb_f)
                    nn_e=NNF(1+nn_f,nf_e,nb_e)
                    np=NPNE(nn_e,nb_e,ne)
                    NPLIST(np)=0
                  ENDDO !nn_f
                ENDIF !external
              ENDDO !no_nf
C             Reduce internal nodes to one version
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                IF(NPLIST(np).EQ.1.AND.NENP(np,0,nr).EQ.8) THEN
C                 Internal and probably not a collapsed node
                  IF(NVJP(nj,np).GT.1) THEN !multiple versions
                    DO nk=1,NKJ(nj,np)
C!!! need to improve for angles
                      SUM=0.0d0
                      DO nv=1,NVJP(nj,np)
                        SUM=SUM+XP(nk,nv,nj,np)
                      ENDDO !nv
                      XP(nk,1,nj,np)=SUM/NVJP(nj,np)
                    ENDDO !nk
                    NVJP(nj,np)=1
                  ENDIF !multiple versions
                ENDIF !internal
              ENDDO !nonode
            ENDDO !njj
          ENDDO !njl
C         Update NVJE
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            DO njl=1,3
              DO njj=1,NJ_LOC(njl,0,nr)
                nj=NJ_LOC(njl,njj,nr)
                nb=NBJ(nj,ne)
                DO nn=1,NNT(nb)
                  np=NPNE(nn,nb,ne)
                  IF(NVJE(nn,nb,nj,ne).GT.NVJP(nj,np)) THEN
                    NVJE(nn,nb,nj,ne)=1 !internal
                  ENDIF
                ENDDO !nn
              ENDDO !njj
            ENDDO !njl
          ENDDO !noelem
C         Update NVJL
          DO no_nl=1,NLLINE(0,nr)
            nl=NLLINE(no_nl,nr)
            DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
              nj=NJ_LOC(NJL_GEOM,njj,nr)
              IF(NPL(1,1,nl).EQ.1) THEN !linear Lagrange basis
                NN_TOT=2
              ELSE IF(NPL(1,1,nl).EQ.2) THEN !quadratic Lagrange
                NN_TOT=3
              ELSE IF(NPL(1,1,nl).EQ.3) THEN !cubic Lagrange
                NN_TOT=4
              ELSE IF(NPL(1,1,nl).EQ.4) THEN !cubic Hermite
                NN_TOT=2
              ENDIF
              DO nn=1,NN_TOT
                np=NPL(1+nn,1,nl)
                IF(NVJL(nn,nj,nl).GT.NVJP(nj,np)) NVJL(nn,nj,nl)=1
              ENDDO !nn
            ENDDO !njj
          ENDDO !no_nl
C         Update scale factors
          CALL LINSCA(IBT,IDO,0,1,NBJ,NEELEM,NEL,NLL,NLLINE,
     '      NNL,NPL,NPNE,nr,NRE,NVJL,DL,SE,XP,ERROR,*9999)

C
C KAT 23Jan99: Does nothing
C        ELSE IF(XI) THEN !update Xj coords of nodes from Xi values
C
C          ERROR='>>Not implemented: Does nothing!'
C          GOTO 9999

        ELSE IF(INTERFACE) THEN !update interface nodes
          IF(TYPE(1:8).EQ.'POSITION') THEN !update interface position
C           PAOPTI has latest estimates of interface position
            DO nrr=1,2
              DO no_interface=1,NP_INTERFACE(0,nrr)
                NPR=NP_INTERFACE(no_interface,nrr)
C GMH 8/1/97 Update cmgui link
                CALL NODE_CHANGE(NPR,.FALSE.,ERROR,*9999)
                XP(1,1,2,NPR)=XP(1,1,2,NPR)+PAOPTI(2*no_interface-1)
                XP(2,1,2,NPR)=XP(2,1,2,NPR)+PAOPTI(2*no_interface)
              ENDDO
            ENDDO

          ELSE IF(TYPE(1:9).EQ.'INCREMENT') THEN !incr. interface pos
C           PAOPTI has latest estimates of interface position
C           PBOPTI has previous estimates of interface position
            DO nrr=1,2
              DO no_interface=1,NP_INTERFACE(0,nrr)
                NPR=NP_INTERFACE(no_interface,nrr)
C GMH 8/1/97 Update cmgui link
                CALL NODE_CHANGE(NPR,.FALSE.,ERROR,*9999)
                XP(1,1,2,NPR)=XP(1,1,2,NPR)+PAOPTI(2*no_interface-1)
     '            -PBOPTI(2*no_interface-1)
                XP(2,1,2,NPR)=XP(2,1,2,NPR)+PAOPTI(2*no_interface)
     '            -PBOPTI(2*no_interface)
              ENDDO
            ENDDO
C           Update old values
            DO no_interface=1,NP_INTERFACE(0,1)
              PBOPTI(2*no_interface-1)=PAOPTI(2*no_interface-1)
              PBOPTI(2*no_interface)  =PAOPTI(2*no_interface)
            ENDDO

          ELSE IF(TYPE(1:4).EQ.'FLUX') THEN !update interface
            DO nrr=1,2 !flux values
              CALL YPZP(2,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '          NKH(1,1,1,nr),NPNODE,nrr,NVHP(1,1,1,nr),nx,NYNE,NYNP,
     '          YP(1,1,nx),ZA,ZP,ERROR,*9999)
              ne=NEELEM(1,nrr)
              VELOC=CE(2,ne,nx)
              WRITE(OP_STRING,'('' velocity in region '',i1,'' ='','
     '          //'e12.3)') nrr,VELOC
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              DO no_interface=1,NP_INTERFACE(0,nrr)
                NPR=NP_INTERFACE(no_interface,nrr)
                NP3=NP_INTERFACE(no_interface,3)
                ZP(1,1,1,NP3,1)=PAOPTI(2*no_interface-1)
                ZP(2,1,1,NP3,1)=PAOPTI(2*no_interface)
                ZP(1,1,1,NPR,1)=(ZP(1,1,1,NP3,1)-ZP(1,1,2,NP3,1))
     '            /TINCR+VELOC*ZP(2,1,1,NP3,1)
                IF(nrr.EQ.2) THEN
                  ZP(1,1,1,NPR,1)=-ZP(1,1,1,NPR,1)
                ENDIF
                WRITE(OP_STRING,'('' flux at node '',I2,'' is '','
     '            //'E12.3)') NPR,ZP(1,1,1,NPR,1)
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDDO
              CALL ZPYP(2,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '          NKH(1,1,1,nrr),NPNODE,nrr,NVHP(1,1,1,nrr),nx,NYNE,NYNP,
     '          YP(1,1,nx),ZA,ZP,ERROR,*9999)
            ENDDO !nrr

          ELSE IF(TYPE(1:4).EQ.'TIME') THEN
C           update previous time step phi(i=1,2) & eta(i=3) values
            DO I=1,3
              DO no_interface=1,NP_INTERFACE(0,I)
                np=NP_INTERFACE(no_interface,I)
                ZP(1,1,2,np,1)=ZP(1,1,1,np,1)
              ENDDO
            ENDDO
          ENDIF

        ELSE IF(FIT) THEN

          DO njj=1,NUM_FIT(0)
            DO nhj=1,NUM_FIT(njj)
              nj_geom=NJ_FIT(nhj,njj) !where the variable for fit njj
              nj_field=NLH_FIT(nhj,1,njj) !the fit variable to store
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
C GMH 8/1/97 Update cmgui link
                CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                DO nv=1,NVJP(nj_field,np)
                  DO nk=1,NKJ(nj_field,np)
C GMH 15/3/96 Should be xp
C                    XP(nk,nv,nj,np)=ZP(nk,nv,nh,np,1)
                    XP(nk,nv,nj_geom,np)=XP(nk,nv,nj_field,np)
                  ENDDO !nk
                ENDDO !nv
              ENDDO !nonode
            ENDDO !nhj
          ENDDO !njj

        ELSE IF(GEOMETRY) THEN
          IF(.NOT.ATNODES)THEN
            CALL ASSERT(.NOT.UP_NQNP,' >>Update node grid first',
     '        ERROR,*9999)
            
            IF(CBBREV(CO,'EXCEPT',2,noco+1,NTCO,N3CO)) THEN
              CALL PARSE_NODES(NPNODE,NPLIST,noco,NRLIST,NTCO,CO,
     '          ERROR,*9999)
              IF(NPLIST(0).GT.0) THEN
                CALL ILISTRMDUP(NPLIST(0),NPLIST(1),ERROR,*9999)
                EXCEPTNODES=.TRUE.
              ELSE
                EXCEPTNODES=.FALSE.
              ENDIF
            ELSE
              EXCEPTNODES=.FALSE.
            ENDIF
            
            IF(EXCEPTNODES) THEN
              DO npp=1,NPNODE(0,nr)
                np=NPNODE(npp,nr)
                IF(.NOT.INLIST(np,NPLIST(1),NPLIST(0),NP1)) THEN
                  nq=NQNP(np)
                  DO nj=1,NJT
                    XP(1,1,nj,np)=XQ(nj,nq)
                  ENDDO
                ENDIF
              ENDDO
            ELSE
              DO npp=1,NPNODE(0,nr)
                np=NPNODE(npp,nr)
                nq=NQNP(np)
                DO nj=1,NJT
                  XP(1,1,nj,np)=XQ(nj,nq)
                ENDDO
              ENDDO
            ENDIF
          ELSE IF(ATNODES)THEN
            DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
              nj=NJ_LOC(NJL_GEOM,njj,nr)
              DO nv=1,NVJP(nj,np1)
                DO nk=1,NKJ(nj,np1)
                  XP(nk,nv,nj,np1)=XP(nk,nv,nj,np2)
                ENDDO !nk
              ENDDO !nv
            ENDDO !njj (nj)
          ENDIF
          
          ELSE IF(GRID) THEN

C          IF(CBBREV(CO,'BEM_REGION',5,noco+1,NTCO,N3CO)) THEN
C            nr_bem=IFROMC(CO(N3CO+1))
C          ELSE
C            nr_bem=1
C          ENDIF
C          IF(CBBREV(CO,'GRID_REGION',6,noco+1,NTCO,N3CO)) THEN
C            nr_ext=IFROMC(CO(N3CO+1))
C          ELSE
C            nr_ext=2
C          ENDIF

          DO np=1,NPM
            NQNP(np)=0
          ENDDO

          DO npp=1,NPNODE(0,nr_bem)
            np=NPNODE(npp,nr_bem)
            DQT=RMAX

            IF(CPLST(0,1).GT.0) THEN !exclude points
              FOUND=.FALSE.
              DO i=1,CPLST(0,1)
                IF(np.EQ.CPLST(i,1)) FOUND=.TRUE.
              ENDDO
              IF(FOUND) THEN
                DO nq=NQR(1,nr_ext),NQR(2,nr_ext)
                  DQ=0.0d0
                  DO nj=1,NJT
                    DQ=DQ+(DABS(XP(1,1,nj,np)-XQ(nj,nq)))**2
                  ENDDO
                  DQ=DSQRT(DQ)
                  IF(DQ.LT.DQT) THEN
                    DQT=DQ
                    NQNP(np)=nq
                  ENDIF
                ENDDO
              ELSE
                DO nq=NQR(1,nr_ext),NQR(2,nr_ext)
                  IF(NWQ(1,nq,1).GT.0) THEN !external grid points only
                    DQ=0.0d0
                    DO nj=1,NJT
                      DQ=DQ+(DABS(XP(1,1,nj,np)-XQ(nj,nq)))**2
                    ENDDO
                    DQ=DSQRT(DQ)
                    IF(DQ.LT.DQT) THEN
                      DQT=DQ
                      NQNP(np)=nq
                    ENDIF
                  ENDIF
                ENDDO
              ENDIF
            ELSE
              DO nq=NQR(1,nr_ext),NQR(2,nr_ext)
                IF(NWQ(1,nq,1).GT.0) THEN !external grid points only
                  DQ=0.0d0
                  DO nj=1,NJT
                    DQ=DQ+(DABS(XP(1,1,nj,np)-XQ(nj,nq)))**2
                  ENDDO
                  DQ=DSQRT(DQ)
                  IF(DQ.LT.DQT) THEN
                    DQT=DQ
                    NQNP(np)=nq
                  ENDIF
                ENDIF
              ENDDO
            ENDIF
            IF(DOP) THEN
              WRITE(OP_STRING,'('' Grid: ,'',I8,'' Node: '',
     '          I6,'' Distance: '',F12.6)') NQNP(np),np,DQT
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO
          UP_NQNP=.FALSE.

          !Check for duplicate mappings
          DO npp=1,NPNODE(0,nr_bem)
            np=NPNODE(npp,nr_bem)
            DO nj=1,NPNODE(0,nr_bem)
              nq=NPNODE(nj,nr_bem)
              IF(np.NE.nq) THEN
                IF(NQNP(np).EQ.NQNP(nq)) THEN
                  WRITE(OP_STRING,'('' WARNING: 2 nodes mapped to '
     '              //'1 grid point '',2I6)') np,nq
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDIF
              ENDIF
            ENDDO
          ENDDO

C This check may not be necessary
C          !Check to see that exclude points are internal grid points
C          IF(CPLST(0,1).GT.0) THEN !exclude points
C            DO npp=1,CPLST(0,1)
C              np=CPLST(npp,1)
C              IF(NWQ(1,NQNP(np),1).NE.0) THEN
C                ERROR='>>Internal point not found for exclude node'
C                GOTO 9999
C              ENDIF
C            ENDDO
C          ENDIF

C          DO nq=1,NQT
C            np=NPQ(nq)
CC GMH 8/1/97 Update cmgui link
C            CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
C            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
C              XP(1,1,nj,np)=XQ(nj,nq)
C            ENDDO !nj
C          ENDDO !nq

        ELSE IF(HANGING) THEN
          IF(JTYP2C.EQ.1) THEN
C            DO nr=1,NRT
              nr=1
              CALL HANGING_NODE_DETECT(IBT,IDO,IDRN,INP,NBJ,NEELEM,
     '          NELIST,NENP,NKJE,NPF,NPLIST,NPNE,NPNODE,nr,NVJE,
     '          NVJP,NWP,SE,SP,XA,XE,XP,ERROR,*9999)
C            ENDDO
          ENDIF

        ELSE IF(POTENTIAL) THEN

          CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqV,NIQ_V,ERROR,*9999)

          CALL ASSERT(.NOT.UP_NQNP,' >>Need to update node grid',
     '      ERROR,*9999)


C! Find grid point which lies at nodes
C          IF(UP_NQNP) THEN
C            DO npp=1,NPLIST(0)
C              np=NPLIST(npp)
C              ne=NENP(np,1,nr_ext)
C              NQNP(np)=0
C
C              NQMAXE=1
C              DO ni=1,NIT(NBJ(1,ne))
C                NQMAXE=NQMAXE*NQXI(ni,NQS(ne))
C              ENDDO
C
C              DQT=1D10
C              DO nqq=1,NQMAXE
C                nq=NQNE(ne,nqq)
C                DQ=0.0d0
C                DO nj=1,NJT
C                  DQ=DQ+(DABS(XP(1,1,nj,np)-XQ(nj,nq)))**2
C                ENDDO
C                DQ=DSQRT(DQ)
C                IF(DQ.LT.DQT) THEN
C                  DQT=DQ
C                  NQNP(np)=nq
C                ENDIF
C              ENDDO
C              CALL ASSERT(NQNP(np).GT.0,
C     '          'No grid point found for boundary node',ERROR,*9999)
C            ENDDO
C            UP_NQNP=.FALSE.
C          ENDIF

          DO npp=1,NPLIST(0)
            np=NPLIST(npp)
            ne=NENP(np,1,nr_ext)
            IF(.NOT.CHECK) THEN
! Find 2 global line numbers for arcs adjacent to grid point
              NL_HOLD(1)=0
              NL_HOLD(2)=0
              DO nl=1,4
                nl_temp=NLL(nl,ne)
                IF(nl_temp.NE.0) THEN
                  IF((NPL(1,0,nl_temp).EQ.1).AND.
     '              (NPL(2,1,nl_temp).EQ.np)) THEN
                    NL_HOLD(2)=nl_temp
                    NL_HOLD(1)=NPL(3,0,nl_temp)
                  ELSE IF((NPL(1,0,nl_temp).EQ.1)
     '              .AND.(NPL(3,1,nl_temp).EQ.np)) THEN
                    NL_HOLD(1)=nl_temp
                    NL_HOLD(2)=NPL(2,0,nl_temp)
                  ENDIF
                ENDIF
              ENDDO

! Calculate arc-length derivatives
              DST(npp)=0.0d0
              nl=NL_HOLD(1)
              IF(nl.NE.0) THEN
                SUM=0.0d0
                ne1=NEL(NEL(0,nl),nl) !grid points in highest region
                DQ=1.0d0/(DBLE(NQXI(1,NQS(ne1)))-1.0d0)

 ! Calculate array XN_LOCAL of coords at line nodes
                ni=NPL(1,0,nl)
                nt=2
                IF(NPL(1,1,nl).EQ.2) nt=3
                IF(NPL(1,1,nl).EQ.3) nt=4
                DO N=1,nt
                  np=NPL(N+1,1,nl)
                  DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                    XN_LOCAL(1,nj,N)=XP(1,1,nj,np)
                    IF((NPL(1,nj,nl).EQ.4).OR.(NPL(1,nj,nl).EQ.6)
     '                .OR.(NPL(1,nj,nl).EQ.7)) THEN
                      ne=NEL(NEL(0,nl),nl)
                      nb=NBJ(nj,ne)
                      NITB=NIT(nb)
                      CALL ASSERT(NITB.EQ.2,'ERROR2',ERROR,*9999)
                      ni2=1+MOD(ni,2)
                      nn=1
                      FOUND=.FALSE.
                      DO WHILE((nn.LE.NNT(nb)).AND.(.NOT.FOUND))
                        IF(np.EQ.NPNE(nn,nb,ne))THEN
                          FOUND=.TRUE.
                        ELSE
                          nn=nn+1
                        ENDIF
                      ENDDO
                      IF(.NOT.FOUND)THEN
                        ERROR='Could not find local node in ARCLEN'
                        GOTO 9999
                      ENDIF
                      DO nk=2,NKT(nn,nb)
                        IF(IDO(nk,nn,ni,nb).EQ.2.AND.
     '                    IDO(nk,nn,ni2,nb).EQ.1) THEN
                          XN_LOCAL(2,nj,N)=XP(nk,1,nj,np)
                        ENDIF
                      ENDDO
                    ENDIF
                  ENDDO
                ENDDO

 ! Contribution from line seg. before grid point
                DO ng=1,NGA
                  XIPOS=1.0d0-(DQ*XIGG(IG(NGA)+ng))
                  W=WG_LOCAL(IG(NGA)+ng)
                  CALL ARCDER(NPL(1,0,nl),DERIVT,DL(1,nl),XIPOS,
     '              XN_LOCAL,ERROR,*9999)
                  SUM=SUM+(W*DERIVT*DQ)
                ENDDO
                DST(npp)=DST(npp)+SUM
              ENDIF !nl.NE.0

              nl=NL_HOLD(2)
              IF(nl.NE.0) THEN
                SUM=0.0d0
                ne1=NEL(NEL(0,nl),nl)
                DQ=1.0d0/(DBLE(NQXI(1,NQS(ne1)))-1.0d0)

 ! Calculate array XN_LOCAL of coords at line nodes
                ni=NPL(1,0,nl)
                nt=2
                IF(NPL(1,1,nl).EQ.2) nt=3
                IF(NPL(1,1,nl).EQ.3) nt=4
                DO N=1,nt
                  np=NPL(N+1,1,nl)
                  DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                    XN_LOCAL(1,nj,N)=XP(1,1,nj,np)
                    IF((NPL(1,nj,nl).EQ.4).OR.(NPL(1,nj,nl).EQ.6)
     '                .OR.(NPL(1,nj,nl).EQ.7)) THEN
                      ne=NEL(NEL(0,nl),nl)
                      nb=NBJ(nj,ne)
                      NITB=NIT(nb)
                      CALL ASSERT(NITB.EQ.2,'ERROR3',ERROR,*9999)
                      ni2=1+MOD(ni,2)
                      nn=1
                      FOUND=.FALSE.
                      DO WHILE((nn.LE.NNT(nb)).AND.(.NOT.FOUND))
                        IF(np.EQ.NPNE(nn,nb,ne))THEN
                          FOUND=.TRUE.
                        ELSE
                          nn=nn+1
                        ENDIF
                      ENDDO
                      IF(.NOT.FOUND)THEN
                        ERROR='Could not find local node in ARCLEN'
                        GOTO 9999
                      ENDIF
                      DO nk=2,NKT(nn,nb)
                        IF(IDO(nk,nn,ni,nb).EQ.2.AND.
     '                    IDO(nk,nn,ni2,nb).EQ.1) THEN
                          XN_LOCAL(2,nj,N)=XP(nk,1,nj,np)
                        ENDIF
                      ENDDO
                    ENDIF
                  ENDDO
                ENDDO

 ! Contribution from line seg. after grid point
                DO ng=1,NGA
                  XIPOS=DQ*XIGG(IG(NGA)+ng)
                  W=WG_LOCAL(IG(NGA)+ng)
                  CALL ARCDER(NPL(1,0,nl),DERIVT,DL(1,nl),XIPOS,
     '              XN_LOCAL,ERROR,*9999)
                  SUM=SUM+(W*DERIVT*DQ)
                ENDDO
                DST(npp)=DST(npp)+SUM
              ENDIF !nl.NE.0
            ENDIF
          ENDDO

! Update boundary conditions
! should look at positions in NYNP to find correct places to
! update bc's
          IF(.NOT.CHECK) THEN
            DO npp=1,NPLIST(0)
              np=NPLIST(npp)
              nq=NQNP(np)
              nq1=NXQ(-1,1,nq,1)
              nq2=NXQ(1,1,nq,1)
              nv=1 !Temporary
              nh=NH_LOC(1,nx_bem)
              ny=NYNP(1,nv,nh,np,0,1,nr_bem)
              nyd=NYNP(2,nv,nh,np,0,1,nr_bem)

              YP(ny,1,nx_bem)=YQ(nq,niqV,1,nx_ext)
              YP(nyd,1,nx_bem)=(YQ(nq2,niqV,1,nx_ext)-
     '          YQ(nq1,niqV,1,nx_ext))/DST(npp)
C              YP(ny,4,nx_bem)=YP(2*np-1,1,nx_bem)
C              YP(nyd,4,nx_bem)=YP(2*np,1,nx_bem)
              YP(ny,4,nx_bem)=YP(ny,1,nx_bem)
              YP(nyd,4,nx_bem)=YP(nyd,1,nx_bem)
              IF(DOP) THEN
                WRITE(OP_STRING,'(''yp,dyp,ds'',3F12.6)')
     '            YP(ny,1,nx_bem),YP(nyd,1,nx_bem),DST(npp)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO
          ELSE !check
            DO npp=1,NPLIST(0)
              np=NPLIST(npp)
              nq=NQNP(np)
              nq1=NWQ(1,nq,1)
              nq2=NWQ(2,nq,1)
              FLUX(npp,1)=-YQ(nq,niqV,1,nx_ext)+
     '          ((4.0d0/3.0d0)*YQ(nq1,niqV,1,nx_ext))-
     '          ((1.0d0/3.0d0)*YQ(nq2,niqV,1,nx_ext))
              nv=1 !Temporary
              nh=NH_LOC(1,nx_bem)
              ny=NYNP(1,nv,nh,np,0,2,nr_bem)
              FLUX(npp,2)=YP(ny,1,nx_bem)
              YP(ny,7,nx_bem)=-FLUX(npp,1)
            ENDDO
            DIFF=0.0d0
            DO npp=1,NPLIST(0)
              DIFF=DIFF+DABS((FLUX(npp,1)-FLUX(npp,2)))
              WRITE(OP_STRING,'(''FLUX: node,grid,bem'',I5,2F12.6)')
     '          NPLIST(npp),FLUX(npp,1),FLUX(npp,2)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
            WRITE(OP_STRING,'(''Total Difference: '',F12.6)') DIFF
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

          NQHOLD=NQNP(NPLIST(1))

        ELSE IF(.NOT.DERIV.AND..NOT.INTERFACE.AND.
     '      .NOT.SC_FAC.AND..NOT.FIT.AND..NOT.MATERIAL) THEN
          IF(UPDATE_ALL) THEN !update all nodal coords
            IF(CONVERT.AND..NOT.REPLACE) THEN !add rc field to polar geom
              nj_field1=NJ_LOC(NJL_FIEL,1,nr)
              nj_field2=NJ_LOC(NJL_FIEL,2,nr)
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
C GMH 8/1/97 Update cmgui link
                CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                nk=1 !temporary PJH 23FEB96
                nv=1 !temporary
                IF(ITYP10(nr).EQ.2) THEN !geom is in cyl. polar coords
                  X_geom=XP(nk,nv,1,np)*DCOS(XP(nk,nv,2,np))
                  Y_geom=XP(nk,nv,1,np)*DSIN(XP(nk,nv,2,np))
                  X_geom=X_geom+XP(nk,nv,nj_field1,np)
                  Y_geom=Y_geom+XP(nk,nv,nj_field2,np)
                  XP(nk,nv,1,np)=DSQRT(X_geom**2+Y_geom**2)
                  IF(X_geom.GT.0.d0) THEN
                    NEW_THETA=DATAN2(Y_geom,X_geom)
                    IF(NEW_THETA.GE.0.d0) THEN
                      XP(nk,nv,2,np)=NEW_THETA
                    ELSE
                      XP(nk,nv,2,np)=NEW_THETA+2.d0*PI
                    ENDIF
                  ENDIF
                ELSE IF(ITYP10(nr).EQ.3) THEN !sph. polar coords
                ELSE IF(ITYP10(nr).EQ.4) THEN !prolate coords
                ENDIF
              ENDDO !nonode

            ELSE IF(.NOT.CONVERT) THEN
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
C GMH 8/1/97 Update cmgui link
                CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                IF(.NOT.ADD_ORTHOG) THEN
                  DO njj=1,NJ_LOC(NJL_FIEL,0,nr) !field variables
                    nj_field=NJ_LOC(NJL_FIEL,njj,nr)
                    nj=NJ_LOC(NJL_GEOM,njj,nr) !corresponding geom var
                    DO nv=1,NVJP(nj_field,np)
                      DO nk=1,NKJ(nj_field,np)
                        IF(DEFORM) THEN
C                       transfers coord to ZP
                          ZP(nk,nv,nj,np,1)=XP(nk,nv,nj_field,np)
                          WRITE(OP_STRING,'('' zp='',D12.3)')
     '                      ZP(nk,nv,nj,np,1)
                          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C news AAY 29Aug1995 added swap to swap ZP and XP
                        ELSE IF(SWAP)THEN
                          TEMP=ZP(nk,nv,nj,np,1)
                          ZP(nk,nv,nj,np,1)=XP(nk,nv,nj,np)
                          XP(nk,nv,nj,np)=TEMP
C newe
                        ENDIF !deform/swap
                        IF(REPLACE) THEN !field var replaces geom var
                          XP(nk,nv,nj,np)=XP(nk,nv,nj_field,np)
                        ELSE IF(.NOT.REPLACE) THEN !field var adds to geom
                          XP(nk,nv,nj,np)=XP(nk,nv,nj,np)
     '                      +XP(nk,nv,nj_field,np)
                        ENDIF !replace
                      ENDDO !nk
                    ENDDO !nv
                  ENDDO !njj
                ELSE !field adds orthogonally
C news AJP 26 July 1997
C Firstly find normal vector at np
                  ne=NENP(np,1,nr) !any element containing np will do
                  CALL GET_TNVECTOR(IBT,IDO,INP,NBJ,NENP(1,0,nr),
     '              NKJE,np,NPF,
     '              NP_INTERFACE,NPNE,nr,NRE,NVJE,nx,NORMAL,SE,
     '              TANGENT,XA,XE,XP,ERROR,*9999)
                  nb=NBJ(1,ne)
                  CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '              NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '              SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
                  INTERFACE2=(NP_INTERFACE(np,0).GT.1).AND.
     '              (NP_INTERFACE(np,1).NE.nr)
                  CALL CALC_CURVCORRECT(IBT(1,1,nb),IDO(1,1,0,nb),
     '              INP(1,1,nb),nb,ne,NPNE(1,1,ne),nr,CE(1,ne,nx),
     '              CURVCORRECT(1,1,1,ne),DXNDS,
     '              SE(1,1,ne),XE,XP,INTERFACE2,ERROR,*9999)
C Add field var - assume only one field variable (the thickness)
                  DO njj=1,NJ_LOC(NJL_GEOM,0,nr) !geometric vars
                    nj=NJ_LOC(NJL_GEOM,njj,nr)
                    nj_field=NJ_LOC(NJL_FIEL,1,nr) !field nj
                    DO nv=1,NVJP(nj_field,np)
                      !add field value in orthogonal direction
                      XP(1,nv,nj,np)=XP(1,nv,nj,np)+
     '                  XP(1,nv,nj_field,np)*NORMAL(nj)
                      DO nk=2,min(NKJ(nj,np),3)
                        !add derivatives (include curvature correction)
                        XP(nk,nv,nj,np)=XP(nk,nv,nj,np)+
     '                    XP(nk,nv,nj_field,np)*NORMAL(nj)+
     '                    XP(1,nv,nj_field,np)*DXNDS(nj,nk-1)
                      ENDDO !nk
                      IF(NKJ(nj,np).EQ.4) THEN
                        XP(4,nv,nj,np)=XP(4,nv,nj,np)+
     '                    XP(4,nv,nj_field,np)*NORMAL(nj)
                        !NB Cross derivative term not being calculated
                        !correctly. Curvature correction difficult to
                        !evaluate (i.e. the d2xnds1ds2 term)
                      ENDIF !nk=4
                    ENDDO !nv
                  ENDDO !nj
C newe
                ENDIF
              ENDDO !nonode
            ENDIF !convert

          ELSE IF(.NOT.UPDATE_ALL) THEN !update spec. coord nj_update
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
C!!! this nj needs checking
              nj_field=NJ_LOC(NJL_FIEL,1,nr)
              CALL ASSERT(nj_field.GT.0,'>>nj_field is zero',
     '          ERROR,*9999)
              DO nv=1,NVJP(nj_update,np)
                DO nk=1,NKJ(nj_update,np)
                  IF(DEFORM) THEN !transfer field coord to ZP
                    ZP(nk,nv,nj_update,np,1)=XP(nk,nv,nj_field,np)
                  ELSE IF(REPLACE) THEN !field var replaces geom var
                    XP(nk,nv,nj_update,np)=XP(nk,nv,nj_field,np)
                  ELSE IF(.NOT.REPLACE) THEN !field var adds to geom
                    XP(nk,nv,nj_update,np)=XP(nk,nv,nj_update,np)
     '                +XP(nk,nv,nj_field,np)
                  ENDIF
                ENDDO !nk
              ENDDO !nv
            ENDDO !nonode
          ENDIF !update_all

        ENDIF !deriv/Hermite/interface/xi/scale_factor/fit

        IF(RESCAL) THEN
          DO nl=1,NLT
            DLL(3,nl)=DL(3,nl) !are old scaling factors
          ENDDO
          IF(DEFORM) THEN
c cpb 13/5/93 Deform case not corrected so that normalised nodal
C derivatives will be achieved.
            DO nl=1,NLT
C              CALL ARCLEN(IDO,NBH,NEL(0,nl),NHP(1,nr,nx),nl,NPL(1,0,nl),
C     '          NPNE,DL,ZP,ERROR,*9999)
C              CALL ARCSCA(IDO,0,0,0,NBH,NEL(0,nl),NHP(1,nr,nx),nl,
C     '          NPL(1,0,nl),NPNE,DL,1.0D-6,ZP,ERROR,*9999)
              CALL ARCSCA(IDO,0,0,0,NBH,NEL(0,nl),nl,
     '          NPL(1,0,nl),NPNE,NVJL(1,1,nl),DL,1.0D-6,ZP,ERROR,*9999)
            ENDDO
          ELSE
c cpb 13/5/93 Changed update node rescale so that the rescaling now
C normalises the nodal derivatives. This results in a non unit scaling
C across the line ie. DL(1,nl) and DL(2,nl) are not the same. In order
C to correct this and achieve arc length scaling an optimised arc length
C needs to be found.
            DO nonode=1,NPNODE(0,nr)
C             Loop over the nodes
              np=NPNODE(nonode,nr)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
C             find number of first derivatives at that node
              NKMAX=0
              DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                nj=NJ_LOC(NJL_GEOM,njj,nr)
                IF(NKJ(nj,np).GT.NKMAX) NKMAX=NKJ(nj,np)
              ENDDO
C              IF(NKMAX.EQ.2) THEN
C                NKTOT=2
C              ELSE IF(NKMAX.EQ.4) THEN
C                NKTOT=3
C              ELSE
C                NKTOT=0
C              ENDIF
C              Loop over the first derivatives at that node
C              DO nk=2,NKTOT
C                 Find the normalising factor K
C news        AAY 29Aug95 allow NKMAX=7 as well for tricubic basis
              IF(NKMAX.EQ.2) THEN
                NKTOT=1
              ELSE IF(NKMAX.EQ.4) THEN
                NKTOT=2
              ELSE IF(NKMAX.EQ.7) THEN
                NKTOT=3
              ELSE
                NKTOT=0
              ENDIF
              DO nk1=1,NKTOT !loop ver 1st derivs
                nk=NKLIST(nk1)
C newe
                K=0.d0
                DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                  nj=NJ_LOC(NJL_GEOM,njj,nr)
                  K=K+XP(nk,1,nj,np)**2.d0
                ENDDO
                K=1.d0/DSQRT(K)
C               Normalise the nodal derivatives
                DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                  nj=NJ_LOC(NJL_GEOM,njj,nr)
                  XP(nk,1,nj,np)=XP(nk,1,nj,np)*K
C                 If bicubic Hermite adjust cross-derivative
                  IF(NKMAX.EQ.4) THEN
                    XP(4,1,nj,np)=XP(4,1,nj,np)*K
                  ENDIF
                ENDDO
C               Find any lines going through that node
                DO nl=1,NLT
                  DO N=1,2
                    IF(NPL(1+N,1,nl).eq.np.AND.
     '                NPL(3+N,1,nl).eq.nk) THEN
C                     Line goes through the node in the direction of
C                     the current derivative number so adjust the
C                     local line scaling factor
                      DL(N,nl)=DL(N,nl)/K
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
c           Find the current est of the line lengths and put it in DL(3)
            DO nl=1,NLT
              CALL ARCLEN(IDO,NBJ,NEL(0,nl),nl,NPL(1,0,nl),NPNE,
     '          NVJL(1,1,nl),DL,XP,ERROR,*9999)
            ENDDO
          ENDIF
          DO nb=1,NBFT
            CALL DLSE(IBT,IDO,nb,NEELEM,NLL,NNL,NPL,
     '        DL,SE,ERROR,*9999)
          ENDDO
C          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
C            nj=NJ_LOC(NJL_GEOM,njj,nr)
C            nb=NBJ(nj,1)
C            IF(NBI(nb).EQ.5) NBI(nb)=2 !to keep DL(1,nl) & DL(2,nl) separate
C          ENDDO
          WRITE(OP_STRING,'('' Nodal values updated from field '
     '      //'values and normalised'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Scale factors adjusted by the '
     '      //'normalising factor'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF !rescal

c cpb 13/5/93 Added update node scale_factor. This sets the local
C scaling factor at each end of the line to the current estimate of the
C line length.
        IF(SC_FAC) THEN
          IF(KTYP26.EQ.2.AND.KTYP27.EQ.5) THEN !optimise line scaling
C           Update line scale params from optimisation params
            DO noopti=1,NTOPTI
              nl=noopti
              DL(1,nl)=PAOPTI(noopti)*DL(3,nl)
              DL(2,nl)=DL(1,nl)
            ENDDO
          ELSE
            DO nl=1,NLT
              DL(1,nl)=DL(3,nl)
              DL(2,nl)=DL(3,nl)
            ENDDO
          ENDIF
          DO nb=1,NBFT
            CALL DLSE(IBT,IDO,nb,NEELEM,NLL,NNL,NPL,
     '        DL,SE,ERROR,*9999)
          ENDDO
          WRITE(OP_STRING,'('' DL(3,nl) copied to DL(1,nl) '
     '      //'and DL(2,nl)'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

        IF(DEFORM.OR.SWAP) THEN
          CALL ZPYP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),NKH(1,1,1,nr),
     '      NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,
     '      YP(1,1,nx),ZA,ZP,ERROR,*9999)
        ENDIF

        IF(MATERIAL) THEN
          nx=1 !temporary
          nj=NJ_LOC(NJL_FIEL,nj_field,nr)
          DO npp=1,NPNODE(0,nr)
            np=NPNODE(npp,nr)
            CP(il,np,nx) = XP(1,1,nj,np)
          ENDDO
          CALL_UPNODE_MAT=.TRUE.
        ENDIF
      ENDIF

      CALL EXITS('UPNODE')
      RETURN
 9999 CALL ERRORS('UPNODE',ERROR)
      CALL EXITS('UPNODE')
      RETURN 1
      END


