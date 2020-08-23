      SUBROUTINE DEMESH(IBT,IDO,INP,LD,LD_NP,NBJ,NBJF,
     '  NEELEM,NEL,NELIST,NELIST2,NENFVC,NENP,NEP,NFF,
     &  NFFACE,NFVC,NKJE,NKEF,NKJ,NLF,NLL,NLLINE,NLLIST,NNB,NNF,NNL,
     &  NODENVC,NODENVCB,NORD,NPF,NPL,NPLIST,NPLIST2,NPNE,NPNF,NPNODE,
     &  NPQ,NQLIST,NRE,NRLIST,NUNK,NVCNODE,NVJE,NVJF,NVJL,NVJP,NXI,NXQ,
     &  NP_INTERFACE,BBM,CE,DF,DL,PG,RG,SE,SF,VC,VC_INIT,WD,WG,XA,XAB,
     &  XE,XG,XID,XIP,XNFV,XP,XQ,ZA,ZD,STRING,ERROR,*)

C#### Subroutine: DEMESH
C###  Description:
C###    DEMESH defines complete mesh for specialized problems.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'genmesh.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grou00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lungex00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'pulm00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LD(NDM),LD_NP(NDM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     &  NEELEM(0:NE_R_M,0:NRM),NEL(0:NELM,NLM),
     &  NELIST(0:NEM),NELIST2(0:NEM),
     &  NENFVC(0:NFVCM,NFVM),NENP(NPM,0:NEPM,0:NRM),NEP(NPM),NFF(6,NEM),
     &  NFFACE(0:NF_R_M,NRM),NFVC(2,0:NFVCM,NVCM),NKJE(NKM,NNM,NJM,NEM),
     &  NKEF(0:4,16,6,NBFM),NKJ(NJM,NPM),NLF(4,NFM),NLL(12,NEM),
     &  NLLINE(0:NL_R_M,0:NRM),NLLIST(0:NLM),NNB(4,4,4,NBFM),
     &  NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),NODENVC(NVCM),NODENVCB(NVCBM),
     &  NORD(5,NE_R_M),NPF(9,NFM),NPL(5,0:3,NLM),NPLIST(0:NPM),
     &  NPLIST2(0:NPM),NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),
     &  NPNODE(0:NP_R_M,0:NRM),NPQ(NQM),NQLIST(0:NQM),NRE(NEM),
     &  NRLIST(0:NRM),NUNK(NKM,NJM,NPM),NVCNODE(2,NP_R_M),
     &  NVJE(NNM,NBFM,NJM,NEM),NVJF(NNM,NBFM,NJM),NVJL(4,NJM,NLM),
     &  NVJP(NJM,NPM),NXI(-NIM:NIM,0:NEIM,0:NEM),
     &  NXQ(-NIM:NIM,0:4,0:NQM,NAM),NP_INTERFACE(0:NPM,0:3),nr_coronary
      REAL*8 BBM(2,NEM),CE(NMM,NEM,NXM),DF(NFM),
     '  DL(3,NLM),PG(NSM,NUM,NGM,NBM),
     '  RG(NGM),SE(NSM,NBFM,NEM),SF(NSM,NBFM),VC(0:NVCM),
     '  VC_INIT(2,NVCM),WD(NJM,NDM),WG(NGM,NBM),XA(NAM,NJM,NEM),
     '  XAB(NORM,NEM),XE(NSM,NJM),XG(NJM,NUM),XID(NIM,NDM),XIP(NIM,NPM),
     '  XNFV(-(NJM+1):NJM,NFVM),XP(NKM,NVM,NJM,NPM),XQ(NJM,NQM),
     '  ZA(NAM,NHM,NCM,NEM),ZD(NJM,NDM)
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER ERR,i,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,IFROMC,
     &  IPFILE,MIN_ORDER,N3CO,nb_airway,nb_delaunay,nb_target,
     &  nb_voronoi,N_BDRY,N_DELTA,NE_START,njj,N_IBDRY,
     &  N_INTNL,NMAX_GENERATION,noelem,noelem_start,nonode,
     &  NOVERSIONS(0:NJM),NP_ATTACH,NP_START,nr,Nrefine,nr_host,
     &  nr_target,num_field,num_gen,nx,PATHEND,POINT_LIMIT,
     &  SCALE_DIST_LIMIT,TERMINAL_ORDER
      REAL*8 ACINUS_VOLUME,angle_max,angle_min,branch_angle,
     &  CONST_EFIELD,D_RATIO,DIAM_RATIO,D_MAX_ORDER,fraction,
     &  LDMinimum,length_limit,length_max,length_ratio,min_length,
     &  MRNA(3),OFFSET(3),radiation_diam,RATIO_LENGTHS,RFROMC,RNA(3),
     &  ROTATION_LIMIT,Spread,SV_FREQ,tumour_coord(3),tumour_diam
      CHARACTER FILE*(MXCH),GROUP_NAME*255,
     &  INPUT_PATH*100,MESH_TYPE*30,STATUS*3
      LOGICAL ABBREV,ADD_SUPER,AIRWAY_MESH,ALL_REGIONS,ARTERIES,
     &  BDM_MESH,CALCU,CALC_XI,CONTINUE,CONTINUOUS,CONVERT_MESH,
     &  CBBREV,DELAUNAY_MESH,EFIELD,EFIELD_VAL,FIELD,FILIO,
     &  FIRST_TIME,FRACTAL,GENER,LADDER_MESH,LUNG,LUNG_TOTAL,
     &  MAKE_GROUP,MATCHING_MESH,MOUSE,ORDER,RADIATION,REDUCE,
     &  REDUCING_DISTANCE_LIMIT,REGULAR_MESH,REVERSE,REPLACE_TRUNK,
     &  RESTRICT_BRANCH_ROTATION,RESTRICT_BRANCH_PLANE,
     &  SET_NJJ_ALVEOLI,
     &  SET_NJJ_RADIUS,SMOOTH_MESH_1D,SYMMETRIC,TERMINAL_MESH,
     &  LPM_FIELD,VEINS,VORONOI_MESH
      INTEGER*4 DISTAL,LDTMP1_PTR,LDTMP2_PTR,NDLIST_PTR,NHOST_PTR,
     &  NE_OLD_PTR,NE_REACTIVATE_PTR,NSTEM_PTR,NE_TEMP_PTR
!     Functions
      INTEGER LEN_TRIM

      CALL ENTERS('DEMESH',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define mesh;l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Defines mesh parameters. Mesh properties are read from or
C###    written to the file FILENAME (with extension .ipmesh) in the
C###    directory specified by PATH.
C###  Parameter:      <coronary_reg #>
C###    Specify the region of coronary trees if using coronary meshes
C###    only.
C###  Parameter:      <RNA_position x y z>
C###    Specify the RNA coordinates. For DNA mesh only
C###  Parameter:      <as (NAME)>
C###    Specify a name to group the mesh elements as.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<RNA_position x y z>'
        OP_STRING(2)=BLANK(1:15)//'<coronary_reg #>'
        OP_STRING(3)=BLANK(1:15)//'<as (NAME)>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
        
C---------------------------------------------------------------------

C#### Command: FEM define mesh;l/p/r/w<;FILENAME[$current]><;(PATH
C###  /example)[$current]> lung
C###  Description:
C###    Defines mesh parameters. Mesh properties are read from or
C###    written to the file FILENAME (with extension .ipmesh) in the
C###    directory specified by PATH. 
C###  Parameter:      <reduce>
C###    Specifies that the mesh will be reduced (for ;r; only).
C###  Parameter:      <order/generation TERMINAL_ORDER[1]>
C###    Specifies the order or generation to reduce to (for ;r; only).
C###  Parameter:      <min_order>
C###    Specifies the minimum order that the existing tree goes down to.  
C###  Parameter:      <total>
C###    Specifies that the total information for the tree mesh is read
C###    or written. i.e. coordinates, fields, elements.
C###  Parameter:      <region [1]>
C###    Specify the region number of the mesh.

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<reduce>'
        OP_STRING(3)=BLANK(1:15)//'<order/generation TERMINAL_ORDER[1]>'
        OP_STRING(4)=BLANK(1:15)//'<min_order #[1]>'
        OP_STRING(5)=BLANK(1:15)//'<total>'
        OP_STRING(6)=BLANK(1:15)//'<region [1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define mesh;c lung arteries/veins
C###  Description:
C###    Create a circulation mesh based upon an existing airway mesh
C###  Parameter:      <arteries/veins [arteries]>
C###    Specify the circulation mesh type.
C###  Parameter:      <from AIRWAY_REGION(#)[1]>
C###    Specify the region number of the existing airway tree.
C###  Parameter:      <add_supernumerary>
C###    Specifies that supernumerary branches will be added to mesh.
C###  Parameter:      <field>
C###    Sets up field structures for radius field, call UPMESH to input
C###    field values.
C###  Parameter:      <region VESSEL_REGION(#)[2]>
C###    Specify the region number of the created circulation mesh.

        OP_STRING(1)=STRING(1:IEND)//';c lung'
        OP_STRING(2)=BLANK(1:15)//'<arteries/veins [arteries]>'
        OP_STRING(3)=BLANK(1:15)//'<from AIRWAY_REGION#[1]> '
        OP_STRING(4)=BLANK(1:15)//'<add_supernumerary>'
        OP_STRING(6)=BLANK(1:15)//'<field>'
        OP_STRING(7)=BLANK(1:15)//'<region VESSEL_REGION#[2]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define mesh;c airway
C###  Description:
C###    Create an airway mesh using a specified generation scheme. If
C###    no scheme is specified then the command just sets up mesh structures
C###    (ordering, fields for radii and alveoli) for the existing elements.
C###  Parameter:      <symmetric/horsfield_delta [symmetric]>
C###    Specify whether the model will be symmetric or a Horsfield
C###    asymmetric 'Delta' model (see Horsfield (1971)).
C###  Parameter:      <delta (N_DELTA) [3]>
C###    Specify the delta value for the Horsfield mesh.
C###  Parameter:      <field>
C###    Sets up the field structures for radius and alveolar fields,
C###    call UPMESH to input field values.
C###  Parameter:      <region (#)[1]>
C###    Specify the single region number for the generated mesh.

        OP_STRING(1)=STRING(1:IEND)//';c airway'
        OP_STRING(2)=BLANK(1:15)//'<symmetric/horsfield [symmetric]>'
        OP_STRING(3)=BLANK(1:15)//'<delta (N_DELTA) [3]>'
        OP_STRING(4)=BLANK(1:15)//'<field>'
        OP_STRING(5)=BLANK(1:15)//'<region (#)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define mesh;c bifurcating_distributive
C###  Description:
C###    Create a volume filling mesh using the bifurcation-distributive
C###    algorithm in 3D only.  Details of the algorithm can be found in
C###    Tawhai et al. Ann. Biomed. Eng. 28(7):793-802 (2000).
C###  Parameter:      <region (#)[1]>
C###    Specify a single host region into which the mesh will be
C###    generated.
C###  Parameter:      <target_region (#)[region+1]>
C###    Specify the region where the mesh will stored.
C###  Parameter:      <elements (#s/all)[1]>
C###    Specify the parent elements from which generation will begin.
C###  Parameter:      <basis_function (#)[1]>
C###    Specify the basis function # for 1D tree mesh.
C###  Parameter:      <point (xi/regular/data)[xi]>
C###    Specify the method for spacing the seed points: xi option spaces
C###    the seed points regularly in xi space; regular option spaces the
C###    seed points regularly in global space; data option is for
C###    reading the points from an .ipdata file (must have already been
C###    read in).
C###  Parameter:      <density (#)[0.01]>
C###    Specify the approximate density of seed points (points/mm^3).
C###  Parameter:      <angle_maximum (#)[180]>
C###    Specify the maximum angle between a parent and daughter branch.
C###  Parameter:      <length_minimum (#)[1.2]>
C###    Specify the limiting length (terminal) for a branch.
C###  Parameter:      <fraction_branching (#)[0.4]>
C###    Specify the branching fraction = (length of generated branch)
C###    /(distance from end of parent to centre of mass of seed points).
C###  Parameter:      <diameter>
C###    Specifies the diameters of vessels in the mesh.        
C###  Parameter:      <min_barnch_length [0.0]>
C###    Allows specification of the minimum length of a terminal element

        OP_STRING(1)=STRING(1:IEND)//';c bifurcating_distributive'
        OP_STRING(2)=BLANK(1:15)//'<region (#)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<target_region (#)[region+1]>'
        OP_STRING(4)=BLANK(1:15)//'<parent_elements (#s/all)[1]>'
        OP_STRING(5)=BLANK(1:15)//'<basis_function (#)[1]>'
        OP_STRING(6)=BLANK(1:15)//'<point (xi/random/regular/data)[xi]>'
        OP_STRING(7)=BLANK(1:15)//'<density (#)[0.01]>'
        OP_STRING(8)=BLANK(1:15)//'<angle_maximum (#)[180]>'
        OP_STRING(9)=BLANK(1:15)//'<length_minimum (#)[1.2]>'
        OP_STRING(10)=BLANK(1:15)//'<fraction_branching (#)[0.4]>'
        OP_STRING(11)=BLANK(1:15)//'<diameter>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define mesh;c delaunay 
C###  Description:
C###    Create a Delaunay mesh from boundary (B), internal boundary (IB)
C###    and internal (IN) nodes defined in node groups.  If the
C###    'define node;c delaunay' command has been called previously,
C###    then node groups (boundary, internal_boundary, internal) will
C###    have created. Otherwise the B, IB, and IN nodes should be
C###    pre-grouped or else listed in this command accordingly.  All of
C###    the nodes in the region must have a specified type (B, IB, or
C###    IN). The only default valuse for this comman are if the above
C###    node groups have already been set up.  The B and IB nodes must
C###    specified, but the IN nodes are optional.
C###  Parameter:      <B_nodes (#s)[boundary]>
C###    Specify the list of boundary nodes.
C###  Parameter:      <IB_nodes (#s)[internal_boundary]>
C###    Specify the list of internal boundary nodes.
C###  Parameter:      <IN_nodes (#s)[internal]>
C###    Specify the list of internal nodes (optional).
C###  Parameter:      <basis_function (#)[1]>
C###    Specify the 2D/3D simplex basis function # for Delaunay
C###    triangles/tetrahedra.
C###  Parameter:      <region (#)[1]>
C###    Specify the region that contains only B, IB, and IN nodes.

        OP_STRING(1)=STRING(1:IEND)//';c delaunay'
        OP_STRING(2)=BLANK(1:15)//'<B_nodes (#s)[boundary]>'
        OP_STRING(3)=BLANK(1:15)//'<IB_nodes (#s)[internal_boundary]>'
        OP_STRING(4)=BLANK(1:15)//'<IN_nodes (#s)[internal]>'
        OP_STRING(5)=BLANK(1:15)//'<basis_function (#)[1]>'
        OP_STRING(6)=BLANK(1:15)//'<region (#)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------


C#### Command: FEM define mesh;c voronoi
C###  Description:
C###    Create a Voronoi mesh from boundary (B), internal boundary (IB)
C###    and internal (IN) nodes defined in node groups.  The nodes must
C###    been triangulated using the define mesh;c delaunay command, or
C###    an external routine.
C###  Parameter:      <B_nodes (#s)[boundary]>
C###    Specify the list of boundary nodes.
C###  Parameter:      <IB_nodes (#s)[internal_boundary]>
C###    Specify the list of internal boundary nodes.
C###  Parameter:      <IN_nodes (#s)[internal]>
C###    Specify the list of internal nodes (optional).
C###  Parameter:      <convert_nodes [false]>
C###    Specify that the Voronoi 'mesh' is to be converted to nodes (at
C###    the cell vertices) and elements (2D triangles covering the
C###    faces). Selecting this option will cause the Voronoi nodes and
C###    elements to overwrite the Delaunay nodes and elements.
C###  Parameter:      <basis_function (#)[1]>
C###    Specify the 2D simplex basis function # for Voronoi face
C###    elements. Used only with 'convert_mesh'.
C###  Parameter:      <region (#)[1]>
C###    Specify the region that contains only B, IB, and IN nodes.

        OP_STRING(1)=STRING(1:IEND)//';c voronoi'
        OP_STRING(2)=BLANK(1:15)//'<B_nodes (#s)[boundary]>'
        OP_STRING(3)=BLANK(1:15)//'<IB_nodes (#s)[internal_boundary]>'
        OP_STRING(4)=BLANK(1:15)//'<IN_nodes (#s)[internal]>'
        OP_STRING(5)=BLANK(1:15)//'<convert_mesh [false]>'
        OP_STRING(6)=BLANK(1:15)//'<basis_function (#)[1]>'
        OP_STRING(7)=BLANK(1:15)//'<region (#)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define mesh;c matching
C###    Create a mesh by copying the geometry of an existing mesh. This
C###    has an option to join the new mesh to the existing one to
C###    create a 'ladder' model. The basis function and region are
C###    the same as for the existing mesh.
C###  Parameter:      <element (#s/group)[all]>
C###    Specify the existing mesh elements to duplicate
C###  Parameter:      <ladder [false]>
C###    Specify whether the new mesh is to be joined to the existing
C###    mesh. Note that this happens at nodes intermediate between
C###    bifurcations.
C###  Parameter:      <reverse [false]>
C###    Specify whether the new mesh has reversed connectivity. i.e.
C###    the Xi direction is reversed.
C###  Parameter:      <region (#)[1]>
C###    Specify the region that contains the existing mesh.
 
         OP_STRING(1)=STRING(1:IEND)//';c matching'
         OP_STRING(2)=BLANK(1:15)//'<element (#s/group)[all]>'
         OP_STRING(3)=BLANK(1:15)//'<ladder [false]>'
         OP_STRING(4)=BLANK(1:15)//'<reverse [false]>'
         OP_STRING(5)=BLANK(1:15)//'<region (#)[1]>'
         CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define mesh;c tumour
C###  Description:
C###    Define a lung tumour region. Set the location and diameter.
C###    Option to specify irradiated region (diameter) and radiation
C###    dose.
C###  Parameter:      <location (x# y# z#)[0.5 0.5 0.5]>
C###    Specify x/y/z location of centre of tumour. Location given as
C###    fractional coordinate between 0.0 (minimum) and 1.0 (maximum.
C###  Parameter:      <diameter (#)[0.0]>
C###    Specify the tumour diameter. Must be defined if location is 
C###    defined.
C###  Parameter:      <irradiated_diameter [0]>
C###    Specify irradiated volume by specifying the diameter.
C###  Parameter:      <dose (#)[0.5]>
C###    Specify the radiation dose between 0.0 and 1.0.

        OP_STRING(1)=STRING(1:IEND)//';c tumour'
        OP_STRING(2)=BLANK(1:15)//'<location (x# y# z#)[0.5 0.5 0.5]>'
        OP_STRING(3)=BLANK(1:15)//'<diameter (#)[0.0]>'
        OP_STRING(4)=BLANK(1:15)//'<irradiated_diameter [0]>'
        OP_STRING(5)=BLANK(1:15)//'<dose (#)[0.5]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)


C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEMESH',ERROR,*9999)
      ELSE
C LKC 15-MAY-2000 Check bases have been defined
        CALL ASSERT(CALL_BASE,
     '    '>> Define Bases first ',ERROR,*9999)
        IPFILE=1 !is input file version number on 23-Jun 1994
        CALL PARSE_QUALIFIERS(' CDLPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
! New AJP 25-1-94
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        
        IF(CBBREV(CO,'CORONARY_REG',2,noco+1,NTCO,N3CO)) THEN
          nr_coronary=IFROMC(CO(N3CO+1))!def region for coronary mesh
        ELSE
          nr_coronary=0
        ENDIF
        
        IF(CBBREV(CO,'RNA_POSITION',3,noco+1,NTCO,N3CO)) THEN
          RNA(1)=RFROMC(CO(N3CO+1))
          RNA(2)=RFROMC(CO(N3CO+2))
          RNA(3)=RFROMC(CO(N3CO+3))
          MRNA(1)=RFROMC(CO(N3CO+4))
          MRNA(2)=RFROMC(CO(N3CO+5))
          MRNA(3)=RFROMC(CO(N3CO+6))
        ELSE
          RNA(1)=0.0d0
          RNA(2)=0.0d0
          RNA(3)=0.0d0
          MRNA(1)=0.0d0
          MRNA(2)=0.0d0
          MRNA(3)=0.0d0
        ENDIF

        IF(CBBREV(CO,'AS',2,noco+1,NTCO,N3CO)) THEN
          MAKE_GROUP=.TRUE.
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          GROUP_NAME=CO(N3CO+1)(IBEG:IEND)
        ELSE
          MAKE_GROUP=.FALSE.
        ENDIF
        
        IF(CBBREV(CO,'PARENT',2,noco+1,NTCO,N3CO)) THEN
          CALL PARSE_ELEMENTS(NEELEM,NELIST2,noco,NRLIST,NTCO,CO,
     '      ERROR,*9999) !gives an element list of parents
        ELSE
!         ne_parent=1
          NELIST2(0)=1
        ENDIF
        
        IF(CBBREV(CO,'IN',2,noco+1,NTCO,N3CO)) THEN
          nr_host=IFROMC(CO(N3CO+1))
        ELSE
          nr_host=1
        ENDIF

        IF(CBBREV(CO,'REFINE',2,noco+1,NTCO,N3CO)) THEN
          Nrefine=IFROMC(CO(N3CO+1))+1
        ELSE
          Nrefine=1
        ENDIF
c        IF(CBBREV(CO,'MINIMUM_LENGTH',3,noco+1,NTCO,N3CO)) THEN
c          MINLENGTH=RFROMC(CO(N3CO+1))
c        ELSE
c          MINLENGTH=0.d0
c        ENDIF
        IF(ADD)THEN
          IF(CBBREV(CO,'APPEND',4,noco+1,NTCO,N3CO)) THEN
            CALL PARSE_ELEMENTS(NEELEM,NELIST2,noco,NRLIST,NTCO,
     &        CO,ERROR,*9999) !list of elements to attach to
          ELSE
            NELIST2(0)=0
          ENDIF
        ENDIF !ADD
        IF(CBBREV(CO,'SPREAD',3,noco+1,NTCO,N3CO)) THEN
          Spread=RFROMC(CO(N3CO+1))
        ELSE
          Spread=0.d0
        ENDIF

C To specify set up of fields with versions to match mesh structure
        LUNG_TOTAL=.FALSE.
        SET_NJJ_RADIUS=.FALSE.
        SET_NJJ_ALVEOLI=.FALSE.
        FIELD=.FALSE.

        IF(CBBREV(CO,'FIELD',4,noco+1,NTCO,N3CO)) THEN
          nr=NRLIST(1) !only for elements from a single region
          FIELD=.TRUE.
          IF(ADD)THEN
            NUM_FIELD = NJ_LOC(NJL_FIEL,0,nr)
          ELSE
            NUM_FIELD=1 !default
          ENDIF
          IF(CBBREV(CO,'NUM_FIELD',3,noco+1,NTCO,N3CO))THEN
            NUM_FIELD=IFROMC(CO(N3CO+1))
          ENDIF
          IF(CBBREV(CO,'NOVERSIONS',3,noco+1,NTCO,N3CO))THEN
            CALL PARSIL(CO(N3CO+1),NJM,NOVERSIONS(0),
     &        NOVERSIONS(1),ERROR,*9999)
          ELSE
            NOVERSIONS(0)=0
            NOVERSIONS(1)=0
          ENDIF
          IF(CBBREV(CO,'RADIUS_FIELD',3,noco+1,NTCO,N3CO))THEN
            njj=IFROMC(CO(N3CO+1))
c            IF(ADD) SET_NJJ_RADIUS=.TRUE.
          ELSE
            njj=1
          ENDIF
          nj_radius=NJ_LOC(NJL_FIEL,njj,nr)
          IF(CBBREV(CO,'ALVEOLAR_FIELD',3,noco+1,NTCO,N3CO))THEN
            njj=IFROMC(CO(N3CO+1))
            IF(ADD) SET_NJJ_ALVEOLI=.TRUE.
          ELSE
            njj=njj+1 !next field above njj for radius
          ENDIF
          nj_alveoli=NJ_LOC(NJL_FIEL,njj,nr)
        ENDIF
C...  SET UP ELEMENT FIELDS (currently only for coupled pulmonary 
C...   perfusion models)
       EFIELD=.FALSE.
       LPM_FIELD=.FALSE.
       EFIELD_VAL=.FALSE.
       IF(CBBREV(CO,'EFIELD',4,noco+1,NTCO,N3CO)) THEN
          nr=NRLIST(1) !only for elements from a single region
          EFIELD=.TRUE.
          IF(CBBREV(CO,'LPM_FIELD',4,noco+1,NTCO,N3CO))THEN
          LPM_FIELD=.TRUE.
            nej_cap=IFROMC(CO(N3CO+1))
          ENDIF
          IF(CBBREV(CO,'VALUE',4,noco+1,NTCO,N3CO))THEN
            EFIELD_VAL=.TRUE.
            CONST_EFIELD=IFROMC(CO(N3CO+1))
          ENDIF
        ENDIF
        
 !set default
        REDUCE=.FALSE.
        TERMINAL_ORDER=1
        ORDER=.FALSE.
        IF(CBBREV(CO,'REDUCE',3,noco+1,NTCO,N3CO)) THEN
          REDUCE=.TRUE.
          ORDER=.TRUE.
          IF(ABBREV(CO(N3CO+1),'ORDER',3)) THEN
            ORDER=.TRUE.
          ELSE IF(ABBREV(CO(N3CO+1),'GENERATION',3)) THEN
            ORDER=.FALSE.
          ENDIF
          TERMINAL_ORDER=IFROMC(CO(N3CO+2))
        ENDIF

        IF(CBBREV(CO,'LUNG',2,noco+1,NTCO,N3CO))THEN
          LUNG=.TRUE.
          IF(CBBREV(CO,'TOTAL',2,noco+1,NTCO,N3CO))THEN
            LUNG_TOTAL=.TRUE.
            IF(CBBREV(CO,'AT_NODE',2,noco+1,NTCO,N3CO))THEN
              NP_ATTACH=IFROMC(CO(N3CO+1))
            ELSE
              NP_ATTACH=0
            ENDIF
          ELSE
            LUNG_TOTAL=.FALSE.
          ENDIF
        ELSE
          LUNG=.FALSE.
        ENDIF

        ARTERIES=.FALSE.
        VEINS=.FALSE.
        IF(CALCU)THEN
          AIRWAY_MESH=.FALSE.
          SMOOTH_MESH_1D=.FALSE.
          LUNG=.FALSE.
          DELAUNAY_MESH=.FALSE.
          VORONOI_MESH=.FALSE.
          ADD_SUPER=.FALSE.
          SYMMETRIC=.FALSE.
          LUNG_TOTAL=.FALSE.
          REGULAR_MESH=.FALSE.
          LUNG_TUMOUR=.FALSE.
          MATCHING_MESH=.FALSE.
          IF(CBBREV(CO,'REGULAR',3,noco+1,NTCO,N3CO))THEN
            JTYP14=1
            REGULAR_MESH=.TRUE.
            nr=NRLIST(1)
          ELSEIF(CBBREV(CO,'SMOOTH1D',3,noco+1,NTCO,N3CO))THEN
            SMOOTH_MESH_1D=.TRUE.
            IF(CBBREV(CO,'ANGLE_MINIMUM',2,noco+1,NTCO,N3CO))THEN
              angle_min=RFROMC(CO(N3CO+1))*PI/180.d0 !to radians
            ELSE
              angle_min = PI
            ENDIF
            IF(CBBREV(CO,'LENGTH_MINIMUM',2,noco+1,NTCO,N3CO))THEN
              length_limit=RFROMC(CO(N3CO+1))
            ELSE
              length_limit = 1.6d0
            ENDIF
            IF(CBBREV(CO,'MAXIMUM',3,noco+1,NTCO,N3CO))THEN
              length_max=RFROMC(CO(N3CO+1))
            ELSE
              length_max = 1.6d0
            ENDIF
            IF(CBBREV(CO,'LDMINIMUM',2,noco+1,NTCO,N3CO))THEN
              LDMinimum=RFROMC(CO(N3CO+1))
            ELSE
              LDMinimum = 0.d0
            ENDIF
          ELSEIF(CBBREV(CO,'LUNG',2,noco+1,NTCO,N3CO))THEN
            LUNG=.TRUE.
            IF(CBBREV(CO,'SYMMETRIC',3,noco+1,NTCO,N3CO))THEN
              SYMMETRIC=.TRUE.
              IF(CBBREV(CO,'NUM_GENERATIONS',5,noco+1,NTCO,N3CO)) THEN
                num_gen=IFROMC(CO(N3CO+1))
              ELSE
                num_gen=1
              ENDIF
              IF(CBBREV(CO,'BRANCH_ANGLE',9,noco+1,NTCO,N3CO)) THEN
                branch_angle=RFROMC(CO(N3CO+1))
              ELSE
                branch_angle=0.785d0
              ENDIF
            ELSEIF(CBBREV(CO,'VEINS',2,noco+1,NTCO,N3CO)) THEN
              VEINS=.TRUE.
              ARTERIES=.FALSE.
            ELSEIF(CBBREV(CO,'ARTERIES',2,noco+1,NTCO,N3CO)) THEN
              VEINS=.FALSE.
              ARTERIES=.TRUE.
            ENDIF
            IF(CBBREV(CO,'FROM',2,noco+1,NTCO,N3CO)) THEN
              nr_host=IFROMC(CO(N3CO+1))
            ELSE
              nr_host=1
            ENDIF
            IF(CBBREV(CO,'ADD_SUPERNUMERARY',5,noco+1,NTCO,N3CO)) THEN
              ADD_SUPER=.TRUE.
              D_RATIO=0.1d0 !Supernumerary to Conventional vessel diam ratio
              SV_FREQ=1.8d0 !default frequency of supernumerary vessels
              IF(CBBREV(CO,'D_RATIO',5,noco+1,NTCO,N3CO)) THEN
                D_RATIO=RFROMC(CO(N3CO+1))
              ENDIF
              IF(CBBREV(CO,'SV_FREQ',5,noco+1,NTCO,N3CO)) THEN
                SV_FREQ=RFROMC(CO(N3CO+1))
              ENDIF
              DIAM_RATIO=1.5997d0 !default - arterial values
              D_MAX_ORDER=30.d0
              IF(CBBREV(CO,'DIAM_RATIO',8,noco+1,NTCO,N3CO)) THEN
                DIAM_RATIO=RFROMC(CO(N3CO+1))
              ENDIF
              IF(CBBREV(CO,'D_MAX_ORDER',5,noco+1,NTCO,N3CO)) THEN
                D_MAX_ORDER=RFROMC(CO(N3CO+1))
              ENDIF
            ENDIF
            MIN_ORDER=0 !default
            nr=NRLIST(1) !mesh in this region
          ENDIF !LUNG
          BDM_MESH=.FALSE.
          IF(CBBREV(CO,'BIFURCATING_DISTRIBUTIVE',4,noco+1,NTCO,
     '      N3CO))THEN
            BDM_MESH=.TRUE.
            !set default values
            CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,
     &        *9999)
            nr=NRLIST(1)
            nr_target=nr+1
            nb_target=1
            POINT_LIMIT=1
            angle_max=60.d0*PI/180.d0 !to radians
            angle_min=60.d0*PI/180.d0 !to radians
            length_limit=0.5d0
            fraction=0.4d0
            min_length=0.0d0 !min branch length or terminal branches produced in BD algorithm
            SCALE_DIST_LIMIT=16
            RATIO_LENGTHS=0.75d0
            SMOOTH_MESH_1D=.FALSE.
            CALC_XI=.FALSE.
            CONTINUOUS=.FALSE.
            ROTATION_LIMIT=0.d0
            IF(CBBREV(CO,'TARGET_REGION',3,noco+1,NTCO,N3CO)) 
     '        nr_target=IFROMC(CO(N3CO+1))
            NRLIST(0)=1
            NRLIST(1)=nr_target !temporary, to get parent list
            CALL PARSE_ELEMENTS(NEELEM,NELIST2,noco,NRLIST,NTCO,CO,
     '        ERROR,*9999) !list of parent elements, else all elements
            NRLIST(1)=nr
            IF(CBBREV(CO,'BASIS_FUNCTION',2,noco+1,NTCO,N3CO))
     '        nb_target=IFROMC(CO(N3CO+1))
            IF(CBBREV(CO,'POINT_LIMIT',6,noco+1,NTCO,N3CO))
     '        POINT_LIMIT=IFROMC(CO(N3CO+1))
            IF(CBBREV(CO,'MINOR',2,noco+1,NTCO,N3CO))
     '        angle_max=RFROMC(CO(N3CO+1))*PI/180.d0 !to radians
            IF(CBBREV(CO,'MAJOR',2,noco+1,NTCO,N3CO))
     '        angle_min=RFROMC(CO(N3CO+1))*PI/180.d0 !to radians
            IF(CBBREV(CO,'ANGLE_MAXIMUM',2,noco+1,NTCO,N3CO))
     '        angle_max=RFROMC(CO(N3CO+1))*PI/180.d0 !to radians
            IF(CBBREV(CO,'LENGTH_MINIMUM',2,noco+1,NTCO,N3CO))
     '        length_limit=RFROMC(CO(N3CO+1))
            IF(CBBREV(CO,'FRACTION_BRANCHING',2,noco+1,NTCO,N3CO))
     '        fraction=RFROMC(CO(N3CO+1))
            IF(CBBREV(CO,'MIN_BRANCH_LENGTH',2,noco+1,NTCO,N3CO)) !only needed for blood vessels for flow soln.
     '        min_length=RFROMC(CO(N3CO+1))
            IF(CBBREV(CO,'SCALE_DISTANCE',3,noco+1,NTCO,N3CO))THEN
               SCALE_DIST_LIMIT=IFROMC(CO(N3CO+1))
               REDUCING_DISTANCE_LIMIT=.TRUE.
            ELSE
               SCALE_DIST_LIMIT=1000000
               REDUCING_DISTANCE_LIMIT=.FALSE.
            ENDIF
            IF(CBBREV(CO,'RATIO_LENGTHS',3,noco+1,NTCO,N3CO))
     '        RATIO_LENGTHS=RFROMC(CO(N3CO+1))
            IF(CBBREV(CO,'XI',2,noco+1,NTCO,N3CO))THEN
              CALC_XI=.TRUE.
            ENDIF
            IF(CBBREV(CO,'SMOOTH',3,noco+1,NTCO,N3CO))THEN
              SMOOTH_MESH_1D=.TRUE.
            ENDIF
            IF(CBBREV(CO,'CONTINUOUS',3,noco+1,NTCO,N3CO))THEN
              CONTINUOUS=.TRUE.
            ENDIF
            IF(CBBREV(CO,'ROTATION_LIMIT',3,noco+1,NTCO,N3CO))THEN
               ROTATION_LIMIT=RFROMC(CO(N3CO+1))
               RESTRICT_BRANCH_ROTATION=.TRUE.
               RESTRICT_BRANCH_PLANE=.TRUE.
            ELSE
               RESTRICT_BRANCH_ROTATION=.FALSE.
               RESTRICT_BRANCH_PLANE=.FALSE.
            ENDIF
            IF(CBBREV(CO,'NP0',3,noco+1,NTCO,N3CO))THEN
              NP_START=IFROMC(CO(N3CO+1))
            ELSE
              NP_START=NPT(0)
            ENDIF
            IF(CBBREV(CO,'NE0',3,noco+1,NTCO,N3CO))THEN
              NE_START=IFROMC(CO(N3CO+1))
            ELSE
              NE_START=NET(0)
            ENDIF
          ENDIF ! bifurcating-distributive method

          IF(CBBREV(CO,'DELAUNAY',2,noco+1,NTCO,N3CO)) THEN
            DELAUNAY_MESH=.TRUE.
            ! note: need to put in default to group names.
           !set default values
            nr=NRLIST(1)
            nb_delaunay=1
            N_BDRY=0
            N_IBDRY=0
            N_INTNL=0
            IF(CBBREV(CO,'B_NODES',1,noco+1,NTCO,N3CO)) THEN
              !get boundary node group name
              CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
              CALL CUPPER(CO(N3CO+1)(IBEG:IEND),STRING)
              CDATA(1)='NODES'
              !put group nodes into NPLIST
              CALL PARSILG(NPLIST,NPM,CDATA(1),STRING,ERROR,*9999)
              N_BDRY=NPLIST(0)
            ENDIF
            IF(CBBREV(CO,'IB_NODES',2,noco+1,NTCO,N3CO)) THEN
              !get boundary node group name
              CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
              CALL CUPPER(CO(N3CO+1)(IBEG:IEND),STRING)
              CDATA(1)='NODES'
              !put group nodes into NPLIST
              CALL PARSILG(NPLIST2,NPM,CDATA(1),STRING,ERROR,*9999)
              N_IBDRY=NPLIST2(0)
              !add IB node list to B node list; just using the one array
              DO nonode=1,NPLIST2(0)
                NPLIST(N_BDRY+nonode)=NPLIST2(nonode)
              ENDDO !nonode
            ENDIF
            IF(CBBREV(CO,'IN_NODES',2,noco+1,NTCO,N3CO)) THEN
              !get internal node group name
              CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
              CALL CUPPER(CO(N3CO+1)(IBEG:IEND),STRING)
              CDATA(1)='NODES'
              !put group nodes into NPLIST
              CALL PARSILG(NPLIST2,NPM,CDATA(1),STRING,ERROR,*9999)
              N_INTNL=NPLIST2(0)
              !add IN node list to B and IB node list; just using the one array
              DO nonode=1,NPLIST2(0)
                NPLIST(N_BDRY+N_IBDRY+nonode)=NPLIST2(nonode)
              ENDDO !nonode
            ENDIF
            IF(CBBREV(CO,'BASIS_FUNCTION',2,noco+1,NTCO,N3CO))
     '        nb_delaunay=IFROMC(CO(N3CO+1))
          ENDIF !DELAUNAY_MESH

          IF(CBBREV(CO,'VORONOI',2,noco+1,NTCO,N3CO)) THEN
            VORONOI_MESH=.TRUE.
           !set default values
            nr=NRLIST(1)
            nr_target=1
            nb_voronoi=1
            CONVERT_MESH=.FALSE.
            N_BDRY=0
            N_IBDRY=0
            N_INTNL=0
            NMAX_GENERATION=1
            IF(CBBREV(CO,'B_NODES',1,noco+1,NTCO,N3CO)) THEN
              !get boundary node group name
              CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
              CALL CUPPER(CO(N3CO+1)(IBEG:IEND),STRING)
              CDATA(1)='NODES'
              !put group nodes into NPLIST
              CALL PARSILG(NPLIST,NPM,CDATA(1),STRING,ERROR,*9999)
              N_BDRY=NPLIST(0)
            ENDIF
            IF(CBBREV(CO,'IB_NODES',2,noco+1,NTCO,N3CO)) THEN
              !get boundary node group name
              CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
              CALL CUPPER(CO(N3CO+1)(IBEG:IEND),STRING)
              CDATA(1)='NODES'
              !put group nodes into NPLIST
              CALL PARSILG(NPLIST2,NPM,CDATA(1),STRING,ERROR,*9999)
              N_IBDRY=NPLIST2(0)
              !add IB node list to B node list; just using the one array
              DO nonode=1,NPLIST2(0)
                NPLIST(N_BDRY+nonode)=NPLIST2(nonode)
              ENDDO !nonode
            ENDIF
            IF(CBBREV(CO,'IN_NODES',2,noco+1,NTCO,N3CO)) THEN
              !get internal node group name
              CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
              CALL CUPPER(CO(N3CO+1)(IBEG:IEND),STRING)
              CDATA(1)='NODES'
              !put group nodes into NPLIST
              CALL PARSILG(NPLIST2,NPM,CDATA(1),STRING,ERROR,*9999)
              N_INTNL=NPLIST2(0)
              !add IN node list to B and IB node list; just using the one array
              DO nonode=1,NPLIST2(0)
                NPLIST(N_BDRY+N_IBDRY+nonode)=NPLIST2(nonode)
              ENDDO !nonode
            ENDIF
            IF(CBBREV(CO,'CONVERT_MESH',4,noco+1,NTCO,N3CO))THEN
              CONVERT_MESH=.TRUE.
              IF(CBBREV(CO,'BASIS_FUNCTION',4,noco+1,NTCO,N3CO))
     '          nb_voronoi=IFROMC(CO(N3CO+1))
              IF(CBBREV(CO,'TARGET_REGION',4,noco+1,NTCO,N3CO))
     '          nr_target=IFROMC(CO(N3CO+1))
              IF(CBBREV(CO,'GENERATION',3,noco+1,NTCO,N3CO))
     '          NMAX_GENERATION=IFROMC(CO(N3CO+1))
            ENDIF !CONVERT
          ENDIF !VORONOI         
          IF(CBBREV(CO,'AIRWAY',3,noco+1,NTCO,N3CO)) THEN
            AIRWAY_MESH=.TRUE.
            IF(NRLIST(0).GT.1)THEN
              WRITE(OP_STRING,
     &          '('' WARNING: only implemented for single region'')')
              CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
            ENDIF
            nr=NRLIST(1)
            
            IF(CBBREV(CO,'BASIS_FUNCTION',3,noco+1,NTCO,N3CO))THEN
              nb_airway=IFROMC(CO(N3CO+1))
            ELSE
              nb_airway=1
            ENDIF

            IF(CBBREV(CO,'AT_NODE',2,noco+1,NTCO,N3CO))THEN
              NP_ATTACH=IFROMC(CO(N3CO+1))
            ELSE
              NP_ATTACH=0
            ENDIF

            MESH_TYPE='DEFAULT'
            IF(CBBREV(CO,'SYMMETRIC',3,noco+1,NTCO,N3CO))THEN
              MESH_TYPE='SYMMETRIC'
              IF(CBBREV(CO,'NUM_GENERATIONS',5,noco+1,NTCO,N3CO)) THEN
                NTT_GEN=IFROMC(CO(N3CO+1))
              ELSE
                NTT_GEN=1
              ENDIF
              DO i=1,NTT_GEN
                B_BRANCH_GEN(i)=1
                L_BRANCH_GEN(i)=1.d0
              ENDDO !i
            ELSE IF(CBBREV(CO,'HORSFIELD',3,noco+1,NTCO,N3CO))THEN
              MESH_TYPE='HORSFIELD'
              noco=N3CO
              IF(CBBREV(CO,'DELTA',3,noco+1,NTCO,N3CO))THEN
                N_DELTA=IFROMC(CO(N3CO+1))
              ELSE
                N_DELTA=3
              ENDIF
            ELSE IF(CBBREV(CO,'MULTI_HBW',4,noco+1,NTCO,N3CO))THEN
              MESH_TYPE='MULTI_HBW'
            ELSE IF(CBBREV(CO,'LUMPED_PARAMETER',4,noco+1,NTCO,
     &          N3CO))THEN
              MESH_TYPE='LUMPED_PARAMETER'
              IF(CBBREV(CO,'VOLUME',3,noco+1,NTCO,N3CO)) THEN
                ACINUS_VOLUME=RFROMC(CO(N3CO+1))
              ELSE
                ACINUS_VOLUME=3.0d6/DBLE(NTB)
              ENDIF
              IF(CBBREV(CO,'SURFACE_AREA',3,noco+1,NTCO,N3CO)) THEN
                AIR_BLOOD_SURFACE=RFROMC(CO(N3CO+1))
              ELSE
                AIR_BLOOD_SURFACE=100.d0/NTERMINAL !DBLE(NTB) !m2 (default = 100 m2)
              ENDIF
              IF(CBBREV(CO,'BLOOD_VOLUME',3,noco+1,NTCO,N3CO)) THEN
                CAP_BLOOD_VOL=RFROMC(CO(N3CO+1))
              ELSE
                CAP_BLOOD_VOL=0.194d0/NTERMINAL !DBLE(NTB) !litre (default = 194 ml)
              ENDIF
!               IF(CBBREV(CO,'PERFUSION',3,noco+1,NTCO,N3CO)) THEN
!                 ACINUS_BLOOD_FLOW=RFROMC(CO(N3CO+1))
!               ELSE
!                 ACINUS_BLOOD_FLOW=5.d0/60.d0/NTERMINAL !DBLE(NTB) !litre/s (default = 5 litre/min)
!               ENDIF
              IF(CBBREV(CO,'BARRIER',3,noco+1,NTCO,N3CO)) THEN
                BARRIER_THICKNESS=RFROMC(CO(N3CO+1))
              ELSE
                BARRIER_THICKNESS=1.11d-6 !m
              ENDIF
!               WRITE(OP_STRING,'('' Lumped parameter:  Volume='',D12.3,
!      &          '' mm3, Gas exchange surface='',D12.3,'' m2, Capillary' 
!      &          //' blood volume='',D12.3,'' l,  air-blood barrier '
!      &          //'thickness='',D12.3,'' um'')') ACINUS_VOLUME,
!      &          AIR_BLOOD_SURFACE,CAP_BLOOD_VOL,BARRIER_THICKNESS
!               CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
!               WRITE(OP_STRING,'('' Total lung alveolar surface area = ''
!      &          ,F8.3,'' m2, Total lung capillary blood volume = '',
!      &          D12.3,'' l'')')AIR_BLOOD_SURFACE*NTERMINAL,
!      &          CAP_BLOOD_VOL*NTERMINAL
!               CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              IF(CBBREV(CO,'SPREAD',3,noco+1,NTCO,N3CO)) THEN
                Spread=RFROMC(CO(N3CO+1))
              ELSE
                Spread=0.d0
              ENDIF
            ENDIF
            IF(CBBREV(CO,'FIELD',4,noco+1,NTCO,N3CO))THEN
              FIELD=.TRUE.
            ENDIF
c must specify 'FIELD', so will be done by stuff up higher
c            IF(MESH_TYPE.EQ.'DEFAULT'.OR.FIELD)THEN
c              IF(CBBREV(CO,'RADIUS_FIELD',3,noco+1,NTCO,N3CO))THEN
c                njj=IFROMC(CO(N3CO+1))
c                nj_radius=NJ_LOC(NJL_FIEL,njj,nr)
c                SET_NJJ_RADIUS=.FALSE.
c              ELSE
c                SET_NJJ_RADIUS=.TRUE.
c              ENDIF
c              IF(CBBREV(CO,'ALVEOLAR_FIELD',3,noco+1,NTCO,N3CO))THEN
c                njj=IFROMC(CO(N3CO+1))
c                nj_alveoli=NJ_LOC(NJL_FIEL,njj,nr)
c                SET_NJJ_ALVEOLI=.FALSE.
c              ELSE
c                SET_NJJ_ALVEOLI=.TRUE.
c              ENDIF
c            ENDIF
          ENDIF
          IF(SMOOTH_MESH_1D)THEN
            CALL MESH_SMOOTH(NBJ,NELIST,NPNE,NVJE,NXI,angle_min,
     &        LDMinimum,length_limit,LENGTH_MAX,XP,ERROR,*9999)
          ENDIF

C ARC 01/2011 Allowing creation of a matching mesh identical to one previously created and connection of the origninal and new meshes
          IF(CBBREV(CO,'MATCHING',2,noco+1,NTCO,N3CO)) THEN
            MATCHING_MESH=.TRUE.
           !set default values
            nr=NRLIST(1)
            DISTAL=0
            LADDER_MESH=.FALSE.
            REPLACE_TRUNK=.FALSE.
            REVERSE=.FALSE.
            TERMINAL_MESH=.FALSE.
            OFFSET(1)=0.d0
            OFFSET(2)=0.d0
            OFFSET(3)=0.d0
            CALL PARSE_ELEMENTS(NEELEM,NELIST2,noco,NRLIST,NTCO,CO,
     &        ERROR,*9999) !gives an element list from the original mesh to duplicate
            IF(CBBREV(CO,'LADDER_MESH',3,noco+1,NTCO,N3CO))THEN
             LADDER_MESH=.TRUE.
            ENDIF
            IF(CBBREV(CO,'REPLACE_TRUNK',3,noco+1,NTCO,N3CO))THEN
              REPLACE_TRUNK=.TRUE.
              IF(CBBREV(CO,'DISTAL',3,noco+1,NTCO,N3CO))THEN
               DISTAL=IFROMC(CO(N3CO+1))
              ELSE
               DISTAL=69
               write(*,*) 'Warning: Assuming your trunk goes to node 69'
              ENDIF
            ENDIF
            IF(CBBREV(CO,'REVERSE',3,noco+1,NTCO,N3CO)) REVERSE=.TRUE.
            IF(CBBREV(CO,'TERMINAL_MESH',3,noco+1,NTCO,N3CO))THEN
             TERMINAL_MESH=.TRUE.
            ENDIF
            IF(CBBREV(CO,'OFFSET',3,noco+1,NTCO,N3CO))THEN
              OFFSET(1)=RFROMC(CO(N3CO+1))
              OFFSET(2)=RFROMC(CO(N3CO+2))
              OFFSET(3)=RFROMC(CO(N3CO+3))
            ENDIF
          ENDIF

C AJS 11/2010 Define tumour and irradiated volume in airway mesh
          IF(CBBREV(CO,'TUMOUR',3,noco+1,NTCO,N3CO))THEN
            LUNG_TUMOUR=.TRUE.
           !set default values
            nr=NRLIST(1)
            tumour_coord(1)=0.5d0 !x
            tumour_coord(2)=0.5d0 !y
            tumour_coord(3)=0.5d0 !z
            tumour_diam=0.0d0
            RADIATION=.FALSE.
            radiation_diam=0.0d0
            radiation_dose=0.0d0
            IF(CBBREV(CO,'LOCATION',3,noco+1,NTCO,N3CO))THEN
              DO i=1,3 !for each coordinate
                tumour_coord(i)=RFROMC(CO(N3CO+i))
                CALL ASSERT((tumour_coord(i).GE.0.0d0.AND.
     &            tumour_coord(i).LE.1.0d0),
     &            '>> Coordinates must be between 0 and 1',ERROR,*9999) 
              ENDDO !i
            ENDIF
            IF(CBBREV(CO,'DIAMETER',3,noco+1,NTCO,N3CO))
     &        tumour_diam=RFROMC(CO(N3CO+1))
            CALL ASSERT((tumour_diam.GE.0.0d0),
     &            '>> Diameter must be >0',ERROR,*9999)
            IF(CBBREV(CO,'IRRADIATED_DIAMETER',3,noco+1,NTCO,N3CO))THEN
              radiation_diam=RFROMC(CO(N3CO+1))
              RADIATION=.TRUE.
              CALL ASSERT((radiation_diam.GE.0.0d0).AND.
     &          (radiation_diam.GT.tumour_diam),
     &          '>> Radiation diameter must be larger than tumour',
     &          ERROR,*9999)
              radiation_dose=0.5d0 !default
              IF(CBBREV(CO,'DOSE',3,noco+1,NTCO,N3CO))
     &          radiation_dose=RFROMC(CO(N3CO+1))
              CALL ASSERT((radiation_dose.GE.0.0d0.AND.radiation_dose
     &          .LE.1.0d0),'>> Radiation dose must be between 0 and 1'
     &          ,ERROR,*9999)
            ENDIF !RADIATION
            CALL DEFINETUMOUR(NPNODE(0,nr),radiation_diam,tumour_coord,
     &        tumour_diam,XP,RADIATION,ERROR,*9999)
          ENDIF !LUNG_TUMOUR

        ENDIF !CALCU

        nx=1 !temporary
        IF(FILIO) THEN
          DIAM_STRAHLER=.FALSE. !diameter-defined Strahler ordering
          CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          PATHEND=LEN_TRIM(FILE)
          CONTINUE=.TRUE.
          DO WHILE(CONTINUE.AND.PATHEND.GT.0)
            IF(FILE(PATHEND:PATHEND).EQ.'/') THEN
              CONTINUE=.FALSE.
            ELSE
              PATHEND=PATHEND-1
            ENDIF
          ENDDO
          IF(PATHEND.EQ.0) THEN
            INPUT_PATH='./'
          ELSE
            INPUT_PATH=FILE(:PATHEND)
          ENDIF
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
C ajp 25-1-94 ALL_REGIONS=.FALSE. !If parse_regions not called, set ALL to be false
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'mesh',
     '        STATUS,ERR,ERROR,*9999)
            
            CALL IPMESH(IBT,IDO,INP,LD,MIN_ORDER,NBJ,NBJF,NEELEM,
     &        NEL,NELIST,NELIST2,NENFVC,NENP,NEP,NFF,
     &        NFFACE,NFVC,NKJE,NKEF,NKJ,NLF,NLL,NLLINE,NNB,NNF,NNL,
     &        NODENVC,NODENVCB,NORD,NP_ATTACH,NP_INTERFACE,NPF,NPL,
     &        NPLIST,NPNE,NPNF,NPNODE,NPQ,NQLIST,nr_coronary,
     &        NRE,Nrefine,nr_host,NRLIST,NVCNODE,NVJE,NVJF,NVJL,NVJP,
     &        NXI,NXQ,TERMINAL_ORDER,BBM,CE(1,1,nx),DF,DL,PG,RG,
     &        MRNA,RNA,SE,SF,Spread,VC,VC_INIT,WG,XA,XAB,XE,XG,XID,XIP,
     &        XNFV,XP,XQ,ZA,ZD,GROUP_NAME,
     &        INPUT_PATH(:LEN_TRIM(INPUT_PATH)),ADD_SUPER,ARTERIES,
     &        CALCU,LUNG_TOTAL,MAKE_GROUP,ORDER,REDUCE,VEINS,ERROR,
     &        *9999)
            CALL CLOSEF(IFILE,ERROR,*9999)
            IF(BACKUP) IOTYPE=2
            CALL GLOBALJ(NBJ,NEELEM,NKJ,NPLIST,NPNE,NPNODE,
     '        ERROR,*9999)
            IF(FIELD) THEN
              CALL GNMESH_FIELD(NBJ,NEELEM,NENP,NPNODE,NKJ,NKJE,
     &          NOVERSIONS,NPNE,nr,NUM_FIELD,NVJE,NVJP,SE,XP,ERROR,
     &          *9999)
              IF(CBBREV(CO,'RADIUS_FIELD',3,noco+1,NTCO,N3CO))THEN
                njj=IFROMC(CO(N3CO+1))
              ELSE
                njj=1
              ENDIF
              nj_radius=NJ_LOC(NJL_FIEL,njj,nr)

            ENDIF
          ENDDO
        ELSE IF(CALCU) THEN
          nr=NRLIST(1)
          noelem_start=NEELEM(0,nr) !records highest elem # before new created

          IF(REGULAR_MESH)THEN
            CALL IPMESH1(NBJ,NEELEM,NENP,NKJE,NKJ,NP_INTERFACE,
     '        NPNE,NPNODE,nr,NRE,NVJE,NVJP,SE,XP,CALCU,ERROR,*9999)
            CALL LINCAL(IBT,IDO,INP,0,NBJ,NEELEM,NEL,NENP,NKJE,NLL,
     '        NLLINE,NNL,NPL,NPNE,NPNODE,NRE,NVJE,NVJL,DL,SE,XP,ERROR,
     '        *9999)
            CALL FACCAL(IBT,IDO,INP,NBJ,NBJF,NEELEM,NFF,NFFACE,NKJE,
     '        NKEF,NLF,NNF,NNL,NPF,NPL,NPNE,NPNF,NRE,NVJE,NVJF,DF,PG,RG,
     '        SE,SF,WG,XA,XE,XG,XP,ERROR,*9999)
            
          ELSE IF(BDM_MESH)THEN
            JTYP14=4
            LDTMP1_PTR=0
            LDTMP2_PTR=0
            NDLIST_PTR=0
            NE_OLD_PTR=0
            NE_TEMP_PTR=0
            NE_REACTIVATE_PTR=0
            NHOST_PTR=0
            NSTEM_PTR=0
            CALL ALLOCATE_MEMORY(NDM+1,0,INTTYPE,NDLIST_PTR,.TRUE.,
     &        ERROR,*9999)
            CALL ALLOCATE_MEMORY(NE_R_M,0,INTTYPE,NE_OLD_PTR,.TRUE.,
     &        ERROR,*9999)
            CALL ALLOCATE_MEMORY(NE_R_M,0,INTTYPE,NE_TEMP_PTR,.TRUE.,
     &        ERROR,*9999)
            CALL ALLOCATE_MEMORY(NDM,0,INTTYPE,LDTMP1_PTR,.TRUE.,
     &        ERROR,*9999)
            CALL ALLOCATE_MEMORY(NP_R_M,0,INTTYPE,LDTMP2_PTR,.TRUE.,
     &        ERROR,*9999)
            CALL ALLOCATE_MEMORY(NEM+1,0,INTTYPE,NE_REACTIVATE_PTR,
     &        .TRUE.,ERROR,*9999)
            CALL ALLOCATE_MEMORY(NP_R_M,0,INTTYPE,NHOST_PTR,.TRUE.,
     &        ERROR,*9999)
            CALL ALLOCATE_MEMORY(NE_R_M*2,0,INTTYPE,NSTEM_PTR,.TRUE.,
     &        ERROR,*9999)

            IF(CONTINUOUS)THEN
              CALL GNBDMESH_CONTINUOUS(LD,LD_NP,nb_target,NBJ,
     &          NEELEM,NELIST2,NENP,%VAL(NE_OLD_PTR),%VAL(NE_TEMP_PTR),
     &          NKJ,NKJE,NNB,NP_INTERFACE,NPNE,NPNODE,nr_target,NRE,
     &          NVJE,NVJP,NXI,POINT_LIMIT,SCALE_DIST_LIMIT,angle_max,
     &          fraction,length_limit,min_length,SE,WD,XP,ZD,ERROR,
     &          *9999)
            ELSE
              CALL GNBDMESH(IBT,IDO,INP,LD,LD_NP,%VAL(LDTMP1_PTR),
     &          %VAL(LDTMP2_PTR),nb_target,NBJ,%VAL(NDLIST_PTR),
     &          NEELEM,NELIST,NELIST2,NENP,%VAL(NE_OLD_PTR),NEP,
     &          %VAL(NE_REACTIVATE_PTR),NE_START,%VAL(NE_TEMP_PTR),
     &          %VAL(NHOST_PTR),NKJ,NKJE,NNB,NORD,NPF,NP_INTERFACE,
     &          NPNE,NPNODE,NP_START,nr,nr_target,NRE,%VAL(NSTEM_PTR),
     &          NVJE,NVJP,NXI,POINT_LIMIT,SCALE_DIST_LIMIT,angle_min,
     &          angle_max,fraction,length_limit,min_length,
     &          RATIO_LENGTHS,ROTATION_LIMIT,SE,WD,XA,XE,XIP,XP,ZD,
     &          CALC_XI,REDUCING_DISTANCE_LIMIT,
     &          RESTRICT_BRANCH_ROTATION,RESTRICT_BRANCH_PLANE,ERROR,
     &          *9999)
           ENDIF
            CALL FREE_MEMORY(NDLIST_PTR,ERROR,*9999)
            CALL FREE_MEMORY(NE_OLD_PTR,ERROR,*9999)
            CALL FREE_MEMORY(NE_TEMP_PTR,ERROR,*9999)
            CALL FREE_MEMORY(LDTMP1_PTR,ERROR,*9999)
            CALL FREE_MEMORY(LDTMP2_PTR,ERROR,*9999)
            CALL FREE_MEMORY(NE_REACTIVATE_PTR,ERROR,*9999)
            CALL FREE_MEMORY(NHOST_PTR,ERROR,*9999)
            CALL FREE_MEMORY(NSTEM_PTR,ERROR,*9999)
            IF(MAKE_GROUP)THEN
              STRING=GROUP_NAME
              i=0
              DO noelem=noelem_start+1,NEELEM(0,nr)
                i=i+1
                NELIST(i)=NEELEM(noelem,nr)
              ENDDO !noelem
              NELIST(0)=i
              CALL GRELEM_SUB(NELIST,STRING,.TRUE.,ERROR,
     &          *9999)
            ENDIF
          ELSEIF(MATCHING_MESH)THEN
            JTYP14=4
            IF(.NOT.REPLACE_TRUNK)THEN
             CALL GNMESH_MATCH(NBJ,NEELEM,NENP,NKJ,NKJE,NORD,
     &           NP_INTERFACE,NPNE,NPNODE,nr,NRE,NVJE,NVJP,NXI,OFFSET,
     &           SE,XAB,XP,EFIELD,LADDER_MESH,REVERSE,TERMINAL_MESH,
     &           ERROR,*9999)  
            ELSE
             CALL GNMESH_MATCH_TRUNK(DISTAL,NBJ,NEELEM,NELIST2,NENP,NKJ,
     &         NKJE,NORD,NP_INTERFACE,NPNE,NPNODE,NPLIST,NPLIST2,nr,NRE,
     &         NVJE,NVJP,NXI,OFFSET,SE,XAB,XP,EFIELD,LADDER_MESH,
     &         REVERSE,TERMINAL_MESH,ERROR,*9999)
            ENDIF

          ELSEIF(AIRWAY_MESH)THEN
            JTYP14=4
            CALL GNMESH1(nb_airway,NBJ,N_DELTA,NEELEM,NELIST2,NENP,NKJ,
     &        NKJE,NORD,NP_INTERFACE,NPNE,NPNODE,nr,NRE,Nrefine,
     &        NVJE,NVJP,NXI,ACINUS_VOLUME,BBM,SE,XP,MESH_TYPE,
     &        ERROR,*9999)
            IF(FIELD) THEN
              CALL GNMESH_FIELD(NBJ,NEELEM,NENP,NPNODE,NKJ,NKJE,
     &          NOVERSIONS,NPNE,nr,NUM_FIELD,NVJE,NVJP,SE,XP,ERROR,
     &          *9999)

            ENDIF
            IF(MAKE_GROUP)THEN
              STRING=GROUP_NAME
              IF(noelem_start.EQ.NEELEM(0,nr))THEN
C             no new elements were added, just mesh structures set up
                DO noelem=1,NEELEM(0,nr)
                  NELIST(noelem)=NEELEM(noelem,nr)
                ENDDO
                NELIST(0)=NEELEM(0,nr)
              ELSE !new elements added
                i=0
                DO noelem=noelem_start+1,NEELEM(0,nr)
                  i=i+1
                  NELIST(i)=NEELEM(noelem,nr)
                ENDDO !noelem
                NELIST(0)=i
              ENDIF
              CALL GRELEM_SUB(NELIST,STRING,.TRUE.,ERROR,*9999)
            ENDIF

            IF(REDUCE) THEN
              CALL MESH_REDUCE(NBJ,NEELEM,NENP,NNB,NORD,
     '          NP_INTERFACE,NPNE,NPNODE,nr,NXI,TERMINAL_ORDER,
     '          ORDER,ERROR,*9999)
            ENDIF
          ELSEIF(LUNG)THEN !define arterial and venous trees
C...        Check airway host mesh already defined
            JTYP14=4 !stochastic fractal tree
            IF(SYMMETRIC) THEN
              IF(CBBREV(CO,'PARENT',3,noco+1,NTCO,N3CO)) THEN
                CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '            ERROR,*9999) !hopefully this gives an element list of parents
              ELSE
!               ne_parent=1
                NELIST(0)=0
              ENDIF
              FRACTAL=.FALSE.
              IF(CBBREV(CO,'FRACTAL',3,noco+1,NTCO,N3CO)) THEN
                FRACTAL=.TRUE.
                IF(CBBREV(CO,'LENGTH',3,noco+1,NTCO,N3CO)) THEN
                  length_ratio=RFROMC(CO(N3CO+1)) !fractal length ratio
                ELSE
                  length_ratio=0.70d0
                ENDIF
              ENDIF
             CALL SYMMETRIC_TREE(NBJ,NEELEM,NELIST,NENP,NKJ,
     &          NKJE,NP_INTERFACE,NPNE,NPNODE,NORD,nr,num_gen,NVJE,
     &          NVJP,NXI,branch_angle,length_ratio,WEIBEL_LENGTH,SE,XP,
     &          FRACTAL,ERROR,*9999)
            ELSEIF(VEINS)THEN
              CALL ASSERT(NEELEM(0,nr_host).NE.0,
     &          '>> Host mesh not defined',ERROR,*9999) 
              CALL GNVEIN(NBJ,NEELEM,NENP,NKJ,NKJE,NP_INTERFACE,NPLIST,
     &          NPNE,NPNODE,nr_host,NRE,nr,NVJE,NVJP,NXI,CE,SE,XP,ERROR,
     &          *9999)
C             ELSE IF(ARTERIES)THEN
C             CALL ASSERT(NEELEM(0,nr_host).NE.0,
C     '          '>> Airway host mesh not defined',ERROR,*9999) 
C           CALL GNARTRY(MIN_ORDER,NBJ(1,NEELEM(1,nr_host)),NBJ,
C     '          NEELEM,NENP,NKJ,NKJE,NORD,NP_INTERFACE,NPNE,NPNODE,
C     '          nr_host,NRE,nr,NVJE,NVJP,NXI,CE,SE,XP,ADD_SUPER,
C     '          ARTERIES,ERROR,*9999)
C           ELSE IF(REORDER) THEN !allocates orders & generates diameters
C            IF(REORDER) THEN
C             DIAM_STRAHLER=.FALSE.
CC              DIAM_STRAHLER=.TRUE.
C              nb=NBJ(1,NEELEM(1,nr))
C              CALL BRANCH_ORD(MIN_ORDER,NEELEM,NORD,nr,NXI,ERROR,*9999)
C              IF(VEINS) THEN
C                CALL DIAM_DEF_STRAHLER(MIN_ORDER,NEELEM,NORD,nr,NXI,CE,
C     '            HORSFIELD_VEIN_DIAM,ERROR,*9999)
C              ELSE
C                CALL DIAM_DEF_STRAHLER(MIN_ORDER,NEELEM,NORD,nr,NXI,CE,
C             '            HORSFIELD_ARTERY_DIAM,ERROR,*9999)
            ENDIF
            IF(ADD_SUPER) THEN !add supernumeraries
 !put in an ASSERT to check mesh has been ordered and diameters allocated
              WRITE(OP_STRING,
     &          '($,'' Order & diameter must be defined '')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              IF(VEINS) THEN
                CALL SUPERNUMERARY(NBJ,NEELEM,NELIST,NELIST2,NENP,
     &            NKJ,NKJE,NORD,NP_INTERFACE,NPNE,NPNODE,nr,NRE,NVJE,
     &            NVJP,NXI,CE,D_RATIO,DIAM_RATIO,D_MAX_ORDER,
     &            HORSFIELD_VEIN_LENGTH,SE,SV_FREQ,XP,ZD,VEINS,ERROR,
     &            *9999)
              ELSE
                CALL SUPERNUMERARY(NBJ,NEELEM,NELIST,NELIST2,NENP,
     &            NKJ,NKJE,NORD,NP_INTERFACE,NPNE,NPNODE,nr,NRE,NVJE,
     &            NVJP,NXI,CE,D_RATIO,DIAM_RATIO,D_MAX_ORDER,
     &            HORSFIELD_ARTERY_LENGTH,SE,SV_FREQ,XP,ZD,VEINS,ERROR,
     &            *9999)
              ENDIF
            ENDIF
            IF(REDUCE) THEN
 !put in an ASSERT to check mesh has been ordered 
              WRITE(OP_STRING,
     &          '($,'' Mesh must be ordered before it is reduced'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              CALL MESH_REDUCE(NBJ,NEELEM,NENP,NNB,NORD,
     '          NP_INTERFACE,NPNE,NPNODE,nr,NXI,TERMINAL_ORDER,
     '          ORDER,ERROR,*9999)
            ENDIF
           IF(FIELD) THEN
 !sets up the field variables to store radius values in fields
             CALL GNMESH_FIELD(NBJ,NEELEM,NENP,NPNODE,NKJ,NKJE,
     &         NOVERSIONS,NPNE,nr,NUM_FIELD,NVJE,NVJP,SE,XP,ERROR,
     &          *9999)
              IF(CBBREV(CO,'RADIUS_FIELD',3,noco+1,NTCO,N3CO))THEN
                njj=IFROMC(CO(N3CO+1))
              ELSE
                njj=1
              ENDIF
              nj_radius=NJ_LOC(NJL_FIEL,njj,nr)
            ENDIF
          ELSEIF(DELAUNAY_MESH)THEN
            CALL DELAUNAY(nb_delaunay,NBJ,NEELEM,NENP,NKJE,NPLIST,
     '        NPNE,NPNODE,nr,NRE,NVJE,N_BDRY,N_IBDRY,
     '        N_INTNL,SE,XP,.FALSE.,ERROR,*9999)
            CALL LINCAL(IBT,IDO,INP,0,NBJ,NEELEM,NEL,NENP,NKJE,NLL,
     '        NLLINE,NNL,NPL,NPNE,NPNODE,NRE,NVJE,NVJL,DL,SE,XP,ERROR,
     '        *9999)
            CALL FACCAL(IBT,IDO,INP,NBJ,NBJF,NEELEM,NFF,NFFACE,NKJE,
     '        NKEF,NLF,NNF,NNL,NPF,NPL,NPNE,NPNF,NRE,NVJE,NVJF,DF,PG,RG,
     '        SE,SF,WG,XA,XE,XG,XP,ERROR,*9999)
          ELSEIF(VORONOI_MESH)THEN
C           NF_LIST_PTR=0
C           CALL ALLOCATE_MEMORY(NFM+1,0,INTTYPE,NF_LIST_PTR,.TRUE.,
C           '        ERROR,*9999)
            CALL VORONOI(IBT,nb_voronoi,N_BDRY,N_IBDRY,N_INTNL,NBJ,
     '        NEELEM,NENFVC,NENP,NFVC,NKJ,NKJE,NMAX_GENERATION,NODENVC,
     '        NODENVCB,NP_INTERFACE,NPLIST,NPNE,NPNODE,nr,NRE,nr_host,
     '        nr_target,NVCNODE,NVJE,NVJP,NXI,SE,VC,VC_INIT,XNFV,XP,ZA,
     '        CONVERT_MESH,ERROR,*9999)
C           CALL FREE_MEMORY(NF_LIST_PTR,ERROR,*9999)
            CALL LINCAL(IBT,IDO,INP,0,NBJ,NEELEM,NEL,NENP,NKJE,NLL,
     '        NLLINE,NNL,NPL,NPNE,NPNODE,NRE,NVJE,NVJL,DL,SE,XP,ERROR,
     '        *9999)
            CALL FACCAL(IBT,IDO,INP,NBJ,NBJF,NEELEM,NFF,NFFACE,NKJE,
     '        NKEF,NLF,NNF,NNL,NPF,NPL,NPNE,NPNF,NRE,NVJE,NVJF,DF,PG,RG,
     '        SE,SF,WG,XA,XE,XG,XP,ERROR,*9999)
          ELSE
            CALL ASSERT(CALL_MESH,
     '        '>>Get real! You haven''t defined a mesh yet!!',ERROR,
     '        *9999)
            
C MLB 19 Feb 1997 both FLOW and TIME_INCREMENT are involved in
C not used after these questions and are in the OMIT file. If this
C is changed please remove them.
C            WRITE(OP_STRING,'($,'' Enter flow at 1st branch ?  '')')
C            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C            READ (IOIP,*) FLOW
C            WRITE(OP_STRING,'($,'' Enter time step ?  '')')
C            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C            READ (IOIP,*) TIME_INCREMENT
          ENDIF
        ENDIF
        
        CALL_MESH=.TRUE.
        IF(jtyp14.NE.2)THEN !Eccentric spheres model
          CALL_EQUA=.FALSE.
          CALL_INIT=.FALSE. !to ensure that NHP     is redefined
          CALL_MATE=.FALSE. !to ensure that CP etc are redefined
        ENDIF
        
C LKC 29-FEB-2000 Moved code from below
        IF(ITYP5(1,1).EQ.2.AND.((ITYP2(1,1).EQ.5.AND.
     '    ITYP3(1,1).EQ.2).OR.(ITYP2(1,1).EQ.11)))THEN
          !special case for large 1D lung trees
          !(we don't care about nunk)
        ELSEIF(ADD_SUPER) THEN
          !don't call for large supernumerary tree problems
        ELSE
C         CS 14/9/98 new added NUNK
c          CALL CALC_NUNK(IDO,NBJ,NENP,NKJE,NKJ,NPNE,NPNODE,
c     '      NRLIST,NUNK,ERROR,*9999)
        ENDIF
        IF(EFIELD)THEN
          CALL GNMESH_EFIELD(NEELEM,nr,CONST_EFIELD,XAB,LPM_FIELD,
     &         EFIELD_VAL,ERROR,*9999)
        ENDIF
      ENDIF


      CALL EXITS('DEMESH')
      RETURN
 9999 CALL ERRORS('DEMESH',ERROR)
      CALL EXITS('DEMESH')
      RETURN 1
      END

