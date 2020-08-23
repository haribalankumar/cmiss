      SUBROUTINE DENODS(IBT,IDO,INP,ISEG,ISNONO,NBH,NBJ,NBJF,NEELEM,
     '  NEL,NELIST,NENP,NFF,NFFACE,NFLIST,NHE,NHP,NKJE,NKEF,NKH,NKJ,NLF,
     '  NLL,NLLINE,NNF,NNL,NPF,NPINTER,NP_INTERFACE,NPL,NPLIST,NPNE,
     '  NPNF,NPNODE,NPNY,NRE,NRLIST,NVCNP,NVHE,NVHP,NVJE,NVJF,NVJL,NVJP,
     '  NW,NXI,NXLIST,NYNE,NYNP,NYNR,DF,DL,PG,RG,SE,SF,WG,XA,XE,
     '  XG,XP,YP,ZD,ZD2,ZP,FIX,CSEG,STRING,ERROR,*)

C#### Subroutine: DENODS
C###  Description:
C###    DENODS defines nodal parameters.

C**** NPNODE(0,nr) is the total number of nodes in region nr.
C**** NPNODE(nonode,nr), nonode=1..NPNODE(0,nr) are the node numbers.

C#### Variable: NPT(0:nr)
C###  Type: INTEGER
C###  Set_up: IPNODE
C###  Description:
C###    NPT(0) is the highest node number in all regions.
C###    NPT(nr) is the highest node number in region nr.

C LKC 7-JUL-1998 adding binary option

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'parameters.inc'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISEG(*),ISNONO(NWM,NPM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NBJF(NJM,NFM),NEELEM(0:NE_R_M,0:NRM),NEL(0:NELM,NLM),
     '  NELIST(0:NEM),NENP(NPM,0:NEPM,0:NRM),NFF(6,NEM),
     '  NFFACE(0:NF_R_M,NRM),NFLIST(0:NFM),NHE(NEM,NXM),
     '  NHP(NPM,0:NRM,NXM),NKJE(NKM,NNM,NJM,NEM),NKEF(0:4,16,6,NBFM),
     '  NKH(NHM,NPM,NCM,0:NRM),NKJ(NJM,NPM),
     '  NLF(4,NFM),NLL(12,NEM),NLLINE(0:NL_R_M,0:NRM),NNF(0:17,6,NBFM),
     '  NNL(0:4,12,NBFM),NPF(9,NFM),NPINTER(3,0:20),
     '  NP_INTERFACE(0:NPM,0:3),NPL(5,0:3,NLM),NPLIST(0:NPM),
     '  NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),
     '  NRE(NEM),NRLIST(0:NRM),NVCNP(NP_R_M),NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NVJE(NNM,NBFM,NJM,NEM),
     '  NVJF(NNM,NBFM,NJM),NVJL(4,NJM,NLM),NVJP(NJM,NPM),NW(NEM,3,NXM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),NXLIST(0:NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM)
      REAL*8 DF(NFM),DL(3,NLM),PG(NSM,NUM,NGM,NBM),RG(NGM),
     '  SE(NSM,NBFM,NEM),SF(NSM,NBFM),
     '  WG(NGM,NBM),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XG(NJM,NUM),XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM),
     '  ZD(NJM,NDM),ZD2(NJM,NDM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(*)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER ERR,i,IBEG,IBEG1,IBEG2,ICHAR,IEND,IEND1,IEND2,IFROMC,IK,
     '  INDEX,INDEX_TEXT,INFO,IPFILE,INSTAT,iw,IWK(6),n,N3CO,nb,NB9,nc,
     &  nd,NDIVISION,NDLIST(7),ne,nh,nhx,nhx_MAX,ni,NIY,nj,nk,nn,NN5(5),
     &  NNTOT,NODE,NO_DERIVS,noelem,N_OFFSET,nolist,nonode,nonode_start,
     &  NOQUES,np,NP1,NP2,NP3,no_nrlist,nr,nr_target,nr_tree,NTINTER,
     &  NTIW,NTLIST,nv,nx,nxc,ny
      INTEGER*4 NFNP_PTR,XP_B_PTR,XP_IB_PTR,XPNP_PTR
      REAL*8 ARC_INCREM,ARC_LENGTH,del_radius,DIFF,PXI,RFROMC,
     '  XI(3),XI1(5),XI2(5),XID,XS(6)
      CHARACTER FILE*(MXCH),GROUP_NAME*255,STATUS*3,TYPE*9
      LOGICAL ALL_REGIONS,CALCU,CBBREV,COMMAND_IP,DELAUNAY_MESH,DELIST,
     '  EXIST,FILIO,FILEIP,FIRST_TIME,FRAME,GENER,HANGING_NODE,
     '  INTERNAL_NODES,INTERPOLATED_NODES,MAKE_GROUP,MOUSE,
     &  SET_PRESS_BASIS,REGULAR,QUAD
      DATA NN5/2,4,5,6,8/
      DATA XI1/0.5D0,0.0D0,0.5D0,1.0D0,0.5D0/,
     '     XI2/0.0D0,0.5D0,0.5D0,0.5D0,1.0D0/
C KAT 24/5/00: also set below.
C     Common blocks should be initialized in a block data.
C MLB 21/3/03: now initialized in BLK50
C      DATA NHT50/2,2,3,4,3,2, 2,2,3,4,3,2, 3,4,4,4,5,3, 3,3,4,4,4,3,
C     '  3,3,4,4,4,3, 3,3,4,4,4,3/

! New Local Variables - remove most

      CHARACTER FILEFORMAT*6
      
  
      CALL ENTERS('DENODS',*9999)
      ICHAR=999

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define nodes;l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Define the coordinate and derivative values of each node. These
C###    are read from or written to the file FILENAME.ipnode in the
C###    directory PATH.
C###  Parameter:      <(reference/dependent NIY# <nc #[1]> <class #[1]>)[reference]>
C###    Specify whether to define reference coordinate or dependent
C###    variable type NIY# values.  For the 'dependent' option, 'nc #'
C###    specifies the type of variable, and 'class #' specifies the
C###    class number (of solve type).
C###  Parameter:      <time DT[0.1]>
C###    For voronoi meshes, specify the time step.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <(ascii/binary)[ascii]>
C###    Specify whether the node file is an ascii or binary file.
C###  Parameter:      <as GROUP_NAME>
C###    Allows you to specify a group name for the nodes        
C###  Parameter:      <nodes #s/GROUP_NAME>
C###    write out only those nodes specified. Preserves the node numbers.
C###    Only applies to the w qualifier.             

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)
     '    //'<(reference/dependent NIY# <nc #[1]> '
     '    //'<class #[1]>)[reference]>'
        OP_STRING(3)=BLANK(1:15)//'<time DT#[0.1]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(5)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(6)=BLANK(1:15)//'<as GROUP_NAME>'
        OP_STRING(7)=BLANK(1:15)//'<nodes #s/GROUP_NAME>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define nodes;l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]> interpolate
C###  Description:
C###    Define a new node with coordinate and derivative values that are
C###    an average of two other nodes' values.
C###  Parameter:      <(DERIV#/all)[all]>
C###    The derivative numbers of the values to interpolate. The remainder
C###    are left uninitialised.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
     '    //' interpolate'
        OP_STRING(2)=BLANK(1:15)//'<(DERIV#/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define nodes;m<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Define nodal coordinates with the mouse. The left mouse
C###    button 'locates' a node. The central mouse button quits
C###    when you have finished defining nodes.
C###  Parameter:      <on WS#[1]>
C###    Specify the workstation.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.

        OP_STRING(1)=STRING(1:IEND)//';m'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<on WS#[1]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define nodes;c
C###  Description:
C###    Create new node or place existing nodes equally spaced along
C###    an arc described by data points.
C###  Parameter:      <at (#s/all)[all]>
C###    Specify the node numbers
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <derivs=#DERIVS[0]>
C###    Specify the number of derivatives at each node. If #DERIVS>0 the
C###    first derivative is in the direction of the arc. If #DERIVS>1
C###    and the number of coordinates is 0, the second derivative is
C###    orthogonal to the first derivative, and the third derivative is
C###    set to zero (even if it doesn't exist).
C###  Description:
C###    Creates new nodes or places existing nodes equally spaced along
C###    an arc described but data points.
CC###  Parameter:      <xi=COORD#[1]>

        OP_STRING(1)=STRING(1:IEND)//';c'
        OP_STRING(2)=BLANK(1:15)//'<at (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
C        OP_STRING(4)=BLANK(1:15)//'<xi=COORD#[1]>'
        OP_STRING(4)=BLANK(1:15)//'<derivs=#DERIVS[0]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define nodes;c hanging element #
C###  Parameter:      <Xi #,#,#[0.5,0.5,0.5]>
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Description:
C###    Create new hanging node in a particular element at
C###    specified Xi coordinates.
C###    This has been set up with the creation of hanging nodes in
C###    mind, hence the 'hanging' qualifier. However, at the moment
C###    it simply creates a new node at the specified Xi coordinates,
C###    no mapping information is set up.

        OP_STRING(1)=STRING(1:IEND)//';c hanging element #'
        OP_STRING(2)=BLANK(1:15)//'<Xi # # #[0.5 0.5 0.5]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define nodes;c nodes=ELEMENT_NODE_TOTAL
C###  Description:
C###    If ELEMENT_NODE_TOTAL=9, convert each element into quadratic
C###    elements, creating new nodes if necessary; otherwise, just
C###    recalculate lines and faces.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.

        OP_STRING(1)=STRING(1:IEND)//';c nodes=ELEMENT_NODE_TOTAL'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define nodes;c delaunay
C###  Description:
C###    Calculates Delaunay boundary (B) and  internal boundary (IB) nodes.
C###    Internal (IN) nodes are calculated upon request, either evenly spaced
C###    in Xi or else using alveolar duct centrelines.  The B and IB nodes
C###    are created over the surface of external faces only.  Node groups are
C###    automatically created as 'boundary', 'internal_boundary', and 'internal'.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <elements (#s/all)[all]>
C###    Specify the elements (list or group) over which to calculate B,
C###    IB, and IN nodes. Elements must be in specified (or default) region.
C###  Parameter:      <target_region (#s)[region]>
C###    Specify the region into which the B, IB, and IN nodes will be
C###    added.
C###  Parameter:      <divisions (#)[2]>
C###    Specify the number of divisions in Xi space along an element.
C###    i.e. 2 divisions will give a pair of B-IB nodes at Xi=0.5.
C###  Parameter:      <radius (#)[0.05]>
C###    Specify the radius of the circumcircle/circumsphere with centre
C###    the surface/edge of a host element, and that will pass through a
C###    pair of B-IB nodes that are normal to the surface.
C###  Parameter:      <internal>
C###    Specify whether to create internal nodes (default is none).
C###  Parameter:      <regular/alveolar [regular]>
C###    To be used only when 'internal' is specified. The default option
C###    is to space the internal nodes regularly in Xi, with the spacing
C###    defined by the 'divisions' parameter.  The second option is to
C###    place internal points along and surrounding the centrelines of
C###    an embedded mesh.
C###  Parameter:      <tree_region (#)[region+1]>
C###    To be used only when 'alveolar' is specified.  This parameter
C###    defines the region in which is an embedded tree mesh; this mesh
C###    is used as a template for calculating IN nodes, and must be in a
C###    separate region.

        OP_STRING(1)=STRING(1:IEND)//';c delaunay'
        OP_STRING(2)=BLANK(1:15)//'<region (#)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<elements (#s/all)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<target_region (#)[region]>'
        OP_STRING(5)=BLANK(1:15)//'<divisions (#)[2]>'
        OP_STRING(6)=BLANK(1:15)//'<radius (#)[0.05]>'
        OP_STRING(7)=BLANK(1:15)//'<internal>'
        OP_STRING(8)=BLANK(1:15)//'<regular/alveolar [regular]>'
        OP_STRING(9)=BLANK(1:15)//'<tree_region (#)[region+1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define nodes
C###  Description:
C###    Define a new node with no file output.
C###  Parameter:      <number ID#[NPT+1]>
C###    Specify the node number
C###  Parameter:      <x=XCOORD#[0.0]>
C###    Specify its x-coordinate
C###  Parameter:      <y=YCOORD#[0.0]>
C###    Specify its y-coordinate
C###  Parameter:      <z=ZCOORD#[0.0]>
C###    Specify its z-coordinate

        OP_STRING(1)=STRING(1:IEND)//' <number ID#[NPT+1]>'
        OP_STRING(2)=BLANK(1:15)//'<x=XCOORD#[0.0]>'
        OP_STRING(3)=BLANK(1:15)//'<y=YCOORD#[0.0]>'
        OP_STRING(4)=BLANK(1:15)//'<z=ZCOORD#[0.0]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DENODS',ERROR,*9999)
      ELSE 
        IPFILE=2 !is input file version number on 16-Nov-1994
        COMMAND_IP=.FALSE.
        CALL PARSE_QUALIFIERS(' CDLMPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
        IF(MOUSE) THEN
          CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)
          CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
        ENDIF
        IF(NTCOQU(noco).EQ.0) THEN !command line input
          COMMAND_IP=.TRUE.
        ENDIF

        NOQUES=0
        FILEIP=.FALSE.
        INTERPOLATED_NODES=.FALSE.

        IF(CBBREV(CO,'REFERENCE',3,noco+1,NTCO,N3CO)) THEN
          TYPE='REFERENCE'
        ELSE IF(CBBREV(CO,'DEPENDENT',3,noco+1,NTCO,N3CO)) THEN
          TYPE='DEPENDENT'
          NIY=IFROMC(CO(N3CO+1))
          IF(CBBREV(CO,'NC',2,noco+1,NTCO,N3CO)) THEN
            nc=IFROMC(CO(N3CO+1))
          ELSE
            nc=1
          ENDIF
C         Allocate and lock nx for solve class
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nxc=NXLIST(1)
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          IF(nx.EQ.0) THEN
            CALL NX_LOC(NX_ALLOCATE_AND_LOCK,nxc,nx,NX_SOLVE,
     '        ERROR,*9999)
          ENDIF
        ELSE
          TYPE='REFERENCE'
        ENDIF

        IF(CBBREV(CO,'TIME',2,noco+1,NTCO,N3CO)) THEN
          DT=RFROMC(CO(N3CO+1))
        ELSE
          DT=0.1D0
        ENDIF

        IF(CBBREV(CO,'OFFSET',3,noco+1,NTCO,N3CO)) THEN
          N_OFFSET=IFROMC(CO(N3CO+1))
        ELSE
          N_OFFSET=0
        ENDIF

        IF(CBBREV(CO,'INTERPOLATE',2,noco+1,NTCO,N3CO)) THEN
          WRITE(OP_STRING,*) ' Note: elements should be defined ',
     '      'before nodes when using this command'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          INTERPOLATED_NODES=.TRUE.
          IF(NTCO.GT.N3CO) THEN
            DELIST=.TRUE.
            CALL PARSIL(CO(N3CO+1),7,NTLIST,NDLIST,ERROR,*9999)
          ELSE
            DELIST=.FALSE.
          ENDIF
        ENDIF

        HANGING_NODE=.FALSE.
        IF(CBBREV(CO,'HANGING',2,noco+1,NTCO,N3CO)) THEN
          HANGING_NODE=.TRUE.
          IF(CBBREV(CO,'ELEMENT',3,noco+1,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),1,NTLIST,ne,ERROR,*9999)
            IF(NTLIST.NE.0) THEN
               IF(CBBREV(CO,'XI',2,noco+1,NTCO,N3CO)) THEN
                 CALL PARSRL(CO(N3CO+1),3,NTLIST,XI,ERROR,*9999)
                 IF(NTLIST.EQ.2) THEN
                   XI(3)=0.0d0
                 ELSE IF((NTLIST.LE.1).OR.(NTLIST.GT.3)) THEN
                   CALL ASSERT(.FALSE.,'>>Incorrect number '
     '             //'of Xi coordinates',ERROR,*9999)
                 ENDIF
                 DO i=1,3
                   IF((XI(i).LT.0.0d0).OR.(XI(i).GT.1.0d0)) THEN
                     CALL ASSERT(.FALSE.,'>>Xi coordinates '
     '               //'out of bounds',ERROR,*9999)
                   ENDIF
                 ENDDO
               ELSE
                 XI(1)=0.5d0
                 XI(2)=0.5d0
                 IF(NIT(ne).EQ.3) THEN
                   XI(3)=0.5d0
                 ELSE
                   XI(3)=0.0d0
                 ENDIF
              ENDIF
            ELSE
               CALL ASSERT(.FALSE.,'>>element number missing',
     '         ERROR,*9999)
            ENDIF
          ELSE
            CALL ASSERT(.FALSE.,'>>element qualifier missing',
     '        ERROR,*9999)
          ENDIF
        ENDIF
        DELAUNAY_MESH=.FALSE.
        IF(CBBREV(CO,'DELAUNAY',5,noco+1,NTCO,N3CO)) THEN
          IF(DOP)THEN
            WRITE(OP_STRING,'('' Delaunay mesh'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          DELAUNAY_MESH=.TRUE.
          !note: only for elements in a single region
          NRLIST(0)=1
          nr=NRLIST(1)
          IF(CBBREV(CO,'ELEMENT',3,noco+1,NTCO,N3CO)) THEN
            CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '        ERROR,*9999)
          ELSE !all elements in region nr
            NELIST(0)=0
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              NELIST(0)=NELIST(0)+1
              NELIST(NELIST(0))=ne
            ENDDO !noelem
          ENDIF !ELEMENT
          IF(DOP)THEN
            WRITE(OP_STRING,'(I4,'' Elements,'',I2,''..'')') NELIST(0),
     '        NELIST(1)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
C Set default values for options
          nr_target=NRLIST(1) !default region in which to put B,IB,IN nodes
          NDIVISION=2 !default number of divisions in an Xi direction
          del_radius=0.05d0 !default circumsphere radius
          INTERNAL_NODES=.FALSE. !default is to not specify internal nodes
          REGULAR=.TRUE. !only used with internal; default is regular Xi spacing
          nr_tree=nr_target+1
          IF(CBBREV(CO,'TARGET_REGION',3,noco+1,NTCO,N3CO))
     '      nr_target=IFROMC(CO(N3CO+1))
          IF(CBBREV(CO,'DIVISIONS',3,noco+1,NTCO,N3CO))
     '      NDIVISION=IFROMC(CO(N3CO+1))
          IF(CBBREV(CO,'RADIUS',3,noco+1,NTCO,N3CO))
     '      del_radius=RFROMC(CO(N3CO+1))
          IF(CBBREV(CO,'INTERNAL',3,noco+1,NTCO,N3CO))
     '      INTERNAL_NODES=.TRUE.
          IF(INTERNAL_NODES)THEN
            IF(CBBREV(CO,'ALVEOLAR',3,noco+1,NTCO,N3CO))THEN
              REGULAR=.FALSE.
              IF(CBBREV(CO,'TREE_REGION',3,noco+1,NTCO,N3CO))THEN
                nr_tree=IFROMC(CO(N3CO+1))
                CALL ASSERT((nr_tree.NE.nr).AND.(nr_tree.NE.nr_target),
     '            '>>Tree must be in separate region',ERROR,*9999)
              ENDIF !TREE_REGION
            ENDIF !ALVEOLAR
          ENDIF !INTERNAL
        ENDIF !DELAUNAY

C LKC 27-JUN-1998 Parse for fileformat
        IF(CBBREV(CO,'BINARY',2,noco+1,NTCO,N3CO)) THEN
          FILEFORMAT='BINARY'
        ELSE
          FILEFORMAT='ASCII'
        ENDIF

C KSB 28 June 2004 Group nodes 'as'
        IF(CBBREV(CO,'AS',2,noco+1,NTCO,N3CO)) THEN
          MAKE_GROUP=.TRUE.
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          GROUP_NAME=CO(N3CO+1)(IBEG:IEND)
          nonode_start=NPNODE(0,NRLIST(1))
        ELSE
          MAKE_GROUP=.FALSE.
        ENDIF

        IF(COMMAND_IP) THEN
          nr=NRLIST(1)
          IF(CBBREV(CO,'NUMBER',1,noco+1,NTCO,N3CO)) THEN
            np=IFROMC(CO(N3CO+1))
          ELSE
            np=NPT(nr)+1
          ENDIF
C GMH 8/1/97 Update cmgui link
          CALL NODE_CHANGE(np,.TRUE.,ERROR,*9999)
          NPNODE(0,nr)=NPNODE(0,nr)+1
          NPNODE(NPNODE(0,nr),nr)=np
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
            NKJ(nj,np)=1
            NVJP(nj,np)=1
          ENDDO !nj
          NPT(nr)=NPNODE(0,nr)
          IF(CBBREV(CO,'X',1,noco+1,NTCO,N3CO)) THEN
            XP(1,1,1,np)=RFROMC(CO(N3CO+1))
          ELSE
            XP(1,1,1,np)=0.d0
          ENDIF
          IF(CBBREV(CO,'Y',1,noco+1,NTCO,N3CO)) THEN
            XP(1,1,2,np)=RFROMC(CO(N3CO+1))
          ELSE
            XP(1,1,2,np)=0.d0
          ENDIF
          IF(CBBREV(CO,'Z',1,noco+1,NTCO,N3CO)) THEN
            XP(1,1,3,np)=RFROMC(CO(N3CO+1))
          ELSE
            XP(1,1,3,np)=0.d0
          ENDIF

          CALL INTERFACE(NP_INTERFACE,NPNODE,ERROR,*9999)

        ELSE IF(CALCU.AND..NOT.HANGING_NODE.AND..NOT.DELAUNAY_MESH) THEN
          IF(CBBREV(CO,'NODES',1,noco+1,NTCO,N3CO)) THEN
            NNTOT=IFROMC(CO(N3CO+1))
            INTERPOLATED_NODES=.TRUE.
          ELSE IF(CBBREV(CO,'AT',1,noco+1,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),NP_R_M,NPLIST(0),NPLIST(1),
     '        ERROR,*9999)
          ELSE
            NPLIST(0)=0
            DO nr=1,NRT
              DO nonode=NPLIST(0)+1,NPLIST(0)+NPNODE(0,nr)
                NPLIST(nonode)=NPNODE(nonode,nr)
              ENDDO
              NPLIST(0)=NPLIST(0)+NPNODE(0,nr)
            ENDDO
          ENDIF
          IF(.NOT.INTERPOLATED_NODES) THEN
            CALL ASSERT(NPLIST(0).GT.1,'>>Too few nodes in list',
     '        ERROR,*9999)
C KAT 6Nov98: Not used.
C            IF(CBBREV(CO,'XI',1,noco+1,NTCO,N3CO)) THEN
C              ni=IFROMC(CO(N3CO+1))
C            ELSE
C              ni=1
C            ENDIF
            IF(CBBREV(CO,'DERIVS',1,noco+1,NTCO,N3CO)) THEN
              NO_DERIVS=IFROMC(CO(N3CO+1))
            ELSE
              NO_DERIVS=0
            ENDIF
          ENDIF
        ENDIF

        FRAME=.FALSE.
        IF(MOUSE) THEN
          INDEX=INDEX_TEXT(0,'WIDTH1','FONT1','BLACK')
        ENDIF

C----------------------------- File io --------------------------------

        IF(FILIO) THEN
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.

            IF(FILEFORMAT(1:6).EQ.'BINARY') THEN

              IF(IOTYPE.EQ.1.OR.IOTYPE.EQ.4.OR.MOUSE.OR.CALCU) THEN
                ERROR='Binary elements only for r/w '
                GOTO 9999
              ENDIF


              IF(IOTYPE.EQ.3) THEN !write out
C*** WRITEOUT BINARY
                CALL IONODE(IOFILE1,NKJ,NPNODE,NRLIST,NVJP,XP,
     '            'WRITE',FILEFORMAT,FILE,'OPEN',ERROR,*9999)
                CALL IONODE(IOFILE1,NKJ,NPNODE,NRLIST,NVJP,XP,
     '            'WRITE',FILEFORMAT,FILE,'NODE_DATA',ERROR,*9999)
                CALL IONODE(IOFILE1,NKJ,NPNODE,NRLIST,NVJP,XP,
     '            'CLOSE',FILEFORMAT,FILE,' ',ERROR,*9999)
C*** READIN BINARY
              ELSE !read in
                CALL IONODE(IOFILE1,NKJ,NPNODE,NRLIST,NVJP,XP,
     '            'READ',FILEFORMAT,FILE,'OPEN',ERROR,*9999)
                CALL IONODE(IOFILE1,NKJ,NPNODE,NRLIST,NVJP,XP,
     '            'READ',FILEFORMAT,FILE,'NODE_DATA',ERROR,*9999)
                CALL IONODE(IOFILE1,NKJ,NPNODE,NRLIST,NVJP,XP,
     '            'CLOSE',FILEFORMAT,FILE,' ',ERROR,*9999)

C*** setup arrays associated with nodes
                NPNODE(0,0)=0
                NPT(0)=0
                DO nr=1,NRLIST(0)
                  NPNODE(0,0)=NPNODE(0,0)+NPNODE(0,nr)
                  NPT(nr)=0
                  DO nonode=1,NPNODE(0,nr)
                    IF(NPNODE(nonode,nr).GT.NPT(nr))
     '                NPT(nr)=NPNODE(nonode,nr)
                  ENDDO
                  IF(NPT(nr).GT.NPT(0)) NPT(0)=NPT(nr)
                ENDDO

C LKC 27-JUN-1999 NP_INTERFACE not setup
C*** Setup regions that interface nodes belong to
                CALL INTERFACE(NP_INTERFACE,NPNODE,ERROR,*9999)
              ENDIF !IOTYPE

            ELSE IF(FILEFORMAT(1:5).EQ.'ASCII') THEN
              CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'node',STATUS,
     '          ERR,ERROR,*9999)

              DO no_nrlist=1,NRLIST(0)
                nr=NRLIST(no_nrlist)
                IF(ALL_REGIONS) CALL PROMPT_REGION_ALL(nr,ERROR,*9999)

                IF(INTERPOLATED_NODES) THEN
                  NTINTER=0 !is initial number of interpolated nodes
 100              FORMAT='($,'' Enter node numbers of new node & 2'
     '              //'  adjacent nodes [exit]: '',3I4)'
                  CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,
     &              FILEIP,FORMAT,3,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              IDATA,IZERO,0,NPM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     &              RMAX,INFO,ERROR,*9999)
                  IF(IDATA(1).GT.0) THEN
                    NP1=IDATA(1)
                    CALL NODE_CHANGE(NP1,.TRUE.,ERROR,*9999)
                    NPNODE(0,nr)=NPNODE(0,nr)+1
                    NPNODE(NPNODE(0,nr),nr)=NP1
                    NP2=IDATA(2)
                    NP3=IDATA(3)
C                 ..update nodal arrays
                    DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                      IF(DELIST) THEN
!deriv numbers specified in DELIST (from 0)
                        DO nolist=1,NTLIST
                          nk=NDLIST(nolist)+1
                          XP(nk,1,nj,NP1)=0.5D0*(XP(nk,1,nj,NP2)+
     '                      XP(nk,1,nj,NP3))
                        ENDDO
                      ELSE !no deriv numbers specified
                        DO nk=1,NKJ(nj,NP2)
                          XP(nk,1,nj,NP1)=0.5D0*(XP(nk,1,nj,NP2)+
     '                      XP(nk,1,nj,NP3))
                        ENDDO
                      ENDIF
                      NKJ(nj,NP1)=NKJ(nj,NP2)
                    ENDDO
                    NTINTER=NTINTER+1 !is current number of inter nodes
                    CALL ASSERT(NTINTER.LE.20,'>>Incr array dimensions',
     '                ERROR,*9999)
                    NPINTER(1,NTINTER)=NP1
                    NPINTER(2,NTINTER)=NP2
                    NPINTER(3,NTINTER)=NP3
                    GO TO 100
                  ENDIF
                  NPINTER(1,0)=NTINTER !is total number of inter nodes
                  CALL ISORT(NPNODE(0,nr),NPNODE(1,nr))
                  DO nonode=1,NPNODE(0,nr)
                    IF(NPNODE(nonode,nr).GT.NPT(nr))
     '                NPT(nr)=NPNODE(nonode,nr)
                    IF(NPNODE(nonode,nr).GT.NPT(0))
     '                NPT(0)=NPNODE(nonode,nr)
                  ENDDO

                ELSE !not.INTERPOLATED_NODES
                  IF(TYPE(1:9).EQ.'REFERENCE') THEN
                    nx=1 !temporary

                    CALL IPNODE('REFERENCE',nc,NHP(1,nr,nx),
     '                NKH(1,1,1,nr),NKJ,N_OFFSET,NP_INTERFACE,NPLIST,
     '                NPNODE,nr,NRLIST,NVHP(1,1,1,nr),NVJP,
     '                XP,ZP,ERROR,*9999)

                    IF(NPT(nr).GT.NPT(0))NPT(0)=NPT(nr)

                  ELSE IF(TYPE(1:9).EQ.'DEPENDENT') THEN
                    IF(.NOT.CALL_EQUA) THEN
                      WRITE(OP_STRING,
     '                  '('' >>WARNING: equation type '
     '                  //'has not been defined:'')')
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      ITYP1(nr,nx)=5 !assume finite elasticity
                      ITYP2(nr,nx)=2
                      WRITE(OP_STRING,
     '                  '(''     - assuming finite elastic '
     '                  //'strains (ITYP1,ITYP2 set=5,2).'')')
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      ITYP6(nr,nx)=2
                      WRITE(OP_STRING,'(''     - assuming nonlinear '
     '                  //'analysis (ITYP6 set=2).'')')
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      IF(NJT.EQ.3) THEN
                        KTYP51(nr)=3 !assume 3D problem
                        WRITE(OP_STRING,
     '                    '(''     - assuming 3D problem '
     '                    //'(KTYP51 set=3)'')')
                        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      ELSE
                        KTYP51(nr)=1 !assume plane stress for 2D
                        WRITE(OP_STRING,
     '                    '(''     - assuming plane stress '
     '                    //'(KTYP51 set=1).'')')
                        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      ENDIF
                      KTYP52(nr)=2 !assume incompressible material

C MLB 14-jul-98 Commenting out to fix alpha batch job
C                      WRITE(OP_STRING,
C     '                  '(''     - assuming incompressible '//
C     '                  '''material (KTYP52 set=2).'')')
C                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

                      NCT(nr,nx)=2 !variable and reaction
C                   Set up NHT50:
C                   Other NHT50 values are const
C                   (defined in data statement above)
                    NHT50(1,1)=NJ_LOC(NJL_GEOM,0,nr)
                    NHT50(2,1)=NJ_LOC(NJL_GEOM,0,nr)
                    NHT50(1,4)=NJ_LOC(NJL_GEOM,0,nr)+1
                    NHT50(2,4)=NJ_LOC(NJL_GEOM,0,nr)
                    NHT50(1,5)=NJ_LOC(NJL_GEOM,0,nr)
                    NHT50(2,5)=NJ_LOC(NJL_GEOM,0,nr)
C                   Set up NHE and NH_LOC
                    nhx_MAX=0
                    DO noelem=1,NEELEM(0,nr)
                      ne=NEELEM(noelem,nr)
                      NHE(ne,nx)=NHT50(KTYP52(nr),KTYP51(nr))
                      IF(NHE(ne,nx).GT.nhx_MAX) nhx_MAX=NHE(ne,nx)
                    ENDDO !noelem
                    CALL CALC_NH_LOC(nhx_MAX,nx,ERROR,*9999)
C                   Set dep var basis NBH same as geom basis NBJ
                    SET_PRESS_BASIS=.FALSE.
                    DO noelem=1,NEELEM(0,nr)
                      ne=NEELEM(noelem,nr)
                      DO nhx=1,NHE(ne,nx)
                        nh=NH_LOC(nhx,nx)
                        IF(nhx.LE.NJ_LOC(NJL_GEOM,0,nr)) THEN
                          NBH(nh,1,ne)=NBJ(nh,ne)
                        ELSE !hyd press
                          NBH(nh,1,ne)=1
                          SET_PRESS_BASIS=.TRUE.
                        ENDIF
C                       copy NBH for nc=2 from nc=1. Force dep var
C                       basis is assumed to be the same as
C                       displacement dep var basis.
                          NBH(nh,2,ne)=NBH(nh,1,ne)
                        ENDDO !nhx
                      ENDDO !noelem (ne)
                      WRITE(OP_STRING,'(''     - assuming dep vars use '
     '                  //'the same bases as the geometric vars.'')')
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      IF(SET_PRESS_BASIS) THEN
                        WRITE(OP_STRING,'(''     - using basis #1 for '
     '                    //'hyd press (N/A for kinematic analysis)'')')
                        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      ENDIF
C                   Set NHP for the problem.
                    CALL CALC_NHP(NBH,NEELEM,NHE(1,nx),NHP(1,0,nx),
     '                NPNE,NPNODE,nr,nx,ERROR,*9999)
C                   Setup dependent version mapping arrays
                    CALL CALC_VERSIONS_DEP(NBH,NBJ,NEELEM,
     '                NHE(1,nx),NHP(1,nr,nx),NPNE,NPNODE,
     '                nr,NVHE,NVHP,NVJE,NVJP,nx,ERROR,*9999)
C                   Initialise all NKH in this region
                    CALL CALC_NKH(NBH,NEELEM,NHP(1,nr,nx),
     '                NKH,NPNE,NPNODE,nr,NW(1,1,nx),nx,ERROR,*9999)
C                   Calculate ny's etc
                      CALL CALC_NY_MAPS_DEP(NBH,NEELEM,
     '                  NHP(1,0,nx),NKH,NP_INTERFACE,NPNODE,
     '                  NPNY(0,1,0,nx),nr,NVHP,nx,NYNE,NYNP,
     '                  NYNR(0,0,1,0,nx),ERROR,*9999)
                      WRITE(OP_STRING,'(''   NOTE: if these options '
     '                  //'are inappropriate then "define equation"'')')
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    ENDIF !.NOT.CALL_EQUA

                    IF(IOTYPE.EQ.3) THEN !write file
                      DO nonode=1,NPNODE(0,nr)
                        np=NPNODE(nonode,nr)
                        DO nhx=1,NHP(np,nr,nx)
                          nh=NH_LOC(nhx,nx)
                          DO nv=1,NVHP(nh,np,nc,nr)
                            DO nk=1,MAX(NKH(nh,np,nc,nr)-KTYP93(1,nr),1)
                              ny=NYNP(nk,nv,nh,np,0,nc,nr)
                              ZP(nk,nv,nh,np,nc)=YP(ny,NIY,nx)
                            ENDDO !nk
                          ENDDO !nv
                        ENDDO !nh
                      ENDDO !nonode (np)
                    ENDIF !iotype=3

                    CALL IPNODE('DEPENDENT',nc,
     '                NHP(1,nr,nx),NKH(1,1,1,nr),NKJ,
     '                N_OFFSET,NP_INTERFACE,NPLIST,
     '                NPNODE,nr,NRLIST,NVHP(1,1,1,nr),NVJP,
     '                XP,ZP,ERROR,*9999)

                    IF(NPT(nr).GT.NPT(0)) NPT(0)=NPT(nr)
C                 NOTE: this doesn't call routine ZPYP as auxiliary
C                 vars & nc=2 vars have not been defined/init'ed
                    DO nonode=1,NPNODE(0,nr)
                      np=NPNODE(nonode,nr)
                      DO nhx=1,NHP(np,nr,nx)
                        nh=NH_LOC(nhx,nx)
                        DO nv=1,NVHP(nh,np,nc,nr)
                          DO nk=1,MAX(NKH(nh,np,nc,nr)-KTYP93(1,nr),1)
                            ny=NYNP(nk,nv,nh,np,0,nc,nr)
                            YP(ny,NIY,nx)=ZP(nk,nv,nh,np,nc)
                          ENDDO ! nk
                        ENDDO ! nv
                      ENDDO ! nh
                    ENDDO ! nonode (np)
                  ENDIF
                ENDIF !interpolated_nodes

                IF(NET(nr).GT.0.AND.IOTYPE.NE.3.AND.
     '            (USE_VORONOI.EQ.0)) THEN
C               ..calculate line topology if node set complete
C               ..not if using Voronoi cells - this is calculated
C                 later as part of mesh reconnection - RGB 17/4/98
                  NPLIST(0)=0
                  DO noelem=1,NEELEM(0,nr)
                    ne=NEELEM(noelem,nr)
                    nb=NBJ(1,ne)
                    DO nn=1,NNT(nb)
                      np=NPNE(nn,nb,ne)
                      DO n=1,NPLIST(0)
                        IF(NPLIST(n).eq.np) GO TO 10
                      ENDDO
                      NPLIST(0)=NPLIST(0)+1
                      NPLIST(NPLIST(0))=np
 10                   CONTINUE
                    ENDDO
                  ENDDO
                  CALL ISORT(NPLIST(0),NPLIST(1))

                  IF(NPLIST(0).EQ.NPNODE(0,nr)) THEN !node set complete

C KAT 14Sep98:  NPL contains the information needed for LINCAL
CC                   ..arc-length scaling requires NXI before SE calc'ed
C                    CALL NENXI(IBT,INP,NBJ,NEELEM,NPNE,NXI,ERROR,*9999)

C                   ..calculate line & face info if node set complete
                    CALL LINCAL(IBT,IDO,INP,0,NBJ,NEELEM,NEL,NENP,
     '                NKJE,NLL,NLLINE,NNL,NPL,NPNE,
     '                NPNODE,NRE,NVJE,NVJL,DL,SE,XP,ERROR,*9999)
                    CALL FACCAL(IBT,IDO,INP,NBJ,NBJF,NEELEM,NFF,NFFACE,
     '                NKJE,NKEF,NLF,NNF,NNL,NPF,NPL,NPNE,NPNF,NRE,
     '                NVJE,NVJF,DF,PG,RG,SE,SF,WG,XA,XE,XG,XP,
     '                ERROR,*9999)

                  ELSE
                    WRITE(OP_STRING,
     '                '('' >>No lines calculated because number of '
     '                //'element nodes does not match number of global '
     '                //'nodes'')')
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDIF
              ENDDO !noregion
              CALL CLOSEF(IFILE,ERROR,*9999)
              CALL ASSERT(NPT(nr).LE.NPM,'>>NPM too small',ERROR,
     '          *9999)
              IF(BACKUP) IOTYPE=2
            ELSE
              ERROR='Unknown FILEFORMAT'
              GOTO 9999
            ENDIF !FILEFORMAT

          ENDDO !first time | backup

C----------------------------- Mouse ip --------------------------------

        ELSE IF(MOUSE) THEN
          nx=1 !temporary
          CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'node',STATUS,
     '      ERR,ERROR,*9999)
          NPLIST(0)=0
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            IF(NPT(nr).GT.0) THEN
              !Clear any existing nodes unless doing a define;add
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                DO i=1,NTIW
                  iw=IWK(i)
                  IF(ISEG(ISNONO(iw,np)).GT.0) THEN
                    CALL ACWK(iw,1,ERROR,*9999)
                    CALL DELETE_SEGMENT(ISNONO(iw,np),ISEG,iw,ERROR,
     '                *9999)
                    CALL DAWK(iw,1,ERROR,*9999)
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
            NPT(nr)=0
            NPNODE(0,nr)=0
            np=0
            nonode=0
c cpb 14/3/94 don't calculate node list from element nodes
C            DO noelem=1,NEELEM(0,nr)
C              ne=NEELEM(noelem,nr)
C              nb=NBJ(1,ne)
C              DO nn=1,NNT(nb)
C                np=NPNE(nn,nb,ne)
C                DO nonode=1,NPLIST(0)
C                  IF(NPLIST(nonode).eq.np) GO TO 20
C                ENDDO
C                NPLIST(0)=NPLIST(0)+1
C                NPLIST(NPLIST(0))=np
C 20             CONTINUE
C              ENDDO
C            ENDDO
C            CALL ISORT(NPLIST(0),NPLIST)
C AJP       NPT_MAX=0         AJP 10-4-93
C           DO nr=1,nr
C             IF(NPT(nr).GT.NPT_MAX) NPT_MAX=NPT(nr)
C           ENDDO
C           np=NPT_MAX
C            np=NPT(0)
C            nonode=NPNODE(0,nr)
            INSTAT=1
            DO WHILE(INSTAT.EQ.1)
              iw=IWK(1)
              nonode=nonode+1
C              IF(nonode.LE.NPLIST(0)) THEN
C                np=NPLIST(nonode)
C              ELSE IF(nonode.GT.NPLIST(0)) THEN
              np=np+1
              CALL ASSERT(np.LE.NPM,'>>Increase NPM',ERROR,*9999)
C              ENDIF
              CALL DENODE(INDEX,INSTAT,ISEG,ISNONO,iw,NHP,
     '          NKH,np,NPNODE,NYNP,XP,.TRUE.,.FALSE.,FIX,FRAME,
     '          CSEG,'NEW',
     '          ERROR,*9999)
              IF(INSTAT.EQ.1) THEN
                NPNODE(nonode,nr)=np
                NPNODE(0,nr)=NPNODE(0,nr)+1
                DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                  NKJ(nj,np)=1
                ENDDO
              ENDIF
            ENDDO
            NPT(nr)=NPNODE(0,nr)
            IF(NPNODE(0,nr).GT.NPNODE(0,0)) NPNODE(0,0)=NPNODE(0,nr)
            IF(NPT(nr).GT.NPT(0)) NPT(0)=NPT(nr)
            CALL IPNODE('REFERENCE',nc,NHP(1,nr,nx),NKH(1,1,1,nr),NKJ,
     &        N_OFFSET,NP_INTERFACE,NPLIST,NPNODE,nr,NRLIST,
     &        NVHP(1,1,1,nr),NVJP,XP,ZP,ERROR,*9999)
            CALL CLOSEF(IFILE,ERROR,*9999)
            CALL INTERFACE(NP_INTERFACE,NPNODE,ERROR,*9999)
          ENDDO

C----------------------------- Calculate -------------------------------

        ELSE IF(CALCU) THEN
          nx=1 !temporary
          IF(INTERPOLATED_NODES) THEN !interp nodes at NNTOT positions
            IF(NNTOT.EQ.9) THEN !biquadratic basis
              QUAD=.FALSE.
              nb=0
              DO WHILE(nb.LT.NBFT.AND..NOT.QUAD) !find nb of biquadr
                nb=nb+1
                IF(NNT(nb).EQ.9) THEN
                  QUAD=.TRUE.
                  NB9=nb
                ENDIF
              ENDDO
              IF(DOP) THEN
                WRITE(OP_STRING,'('' nb9='',I2)') NB9
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              CALL ASSERT(QUAD,'>>no biquadratic basis defined',ERROR,
     '          *9999)
              DO no_nrlist=1,NRLIST(0)
                nr=NRLIST(no_nrlist)
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '              NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '              SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
                  nb=NBJ(1,ne)
                  NPNE(1,NB9,ne)=NPNE(1,nb,ne)
                  NPNE(3,NB9,ne)=NPNE(2,nb,ne)
                  NPNE(7,NB9,ne)=NPNE(3,nb,ne)
                  NPNE(9,NB9,ne)=NPNE(4,nb,ne)
                  DO i=1,5 !to create 5 new nodes
                    nn=NN5(i)
                    XI(1)=XI1(i)
                    XI(2)=XI2(i)
                    DO nj=1,NJT !to find coords of proposed new node position
                      nb=NBJ(nj,ne)
                      XS(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '                  nb,1,XI,XE(1,nj))
                    ENDDO
                    IF(DOP) THEN
                      WRITE(OP_STRING,FMT='('' XI: '',3F7.3)')
     '                  (XI(ni),ni=1,3)
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      WRITE(OP_STRING,FMT='('' XS: '',3F7.3)')
     '                  (XS(nj),nj=1,NJT)
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF

                    DO np=1,NPT(nr) !to search for existing nodes
                      IK=0
                      DO nj=1,NJT
                        DIFF= DABS(XS(nj)-XP(1,1,nj,np))
                        IF((DIFF.LT.(1.0D-3)* DABS(XS(nj)).OR.
     '                    DIFF.LT.1.0D-6).
     '                    OR.(ITYP10(1).eq.2.and.nj.EQ.2.
     '                    AND.( DABS(XP(1,1,1,np)).
     '                    LT.1.0D-6.OR.DABS(DIFF-2.0D0*PI).LT.1.0D-6)).
     '                    OR.(ITYP10(1).eq.3.and.nj.EQ.2.
     '                    AND.( DABS(XP(1,1,1,np)).
     '                    LT.1.0D-6.OR.DABS(DIFF-2.0D0*PI).LT.1.0D-6)).
     '                    OR.(ITYP10(1).eq.3.and.nj.EQ.3.
     '                    AND.( DABS(XP(1,1,1,np)).
     '                    LT.1.0D-6.OR.DABS(DIFF-2.0D0*PI).LT.1.0D-6)).
     '                    OR.(ITYP10(1).eq.4.and.nj.EQ.3.
     '                    AND.( DABS(XP(1,1,2,np)).
     '                    LT.1.0D-6.OR.DABS(XP(1,1,2,np)-PI).LT.1.0D-6.
     '                    OR. DABS(DIFF-2.0D0*PI).LT.1.0D-5))) THEN
                            IK=IK+1
                        ENDIF
                      ENDDO
                      IF(IK.EQ.NJT) THEN
                        EXIST=.TRUE.
                        NODE=np
                        GO TO 41
                      ELSE
                        EXIST=.FALSE.
                      ENDIF
                    ENDDO
 41                 CONTINUE

                    IF(.NOT.EXIST) THEN
                      NPT(nr)=NPT(nr)+1
                      CALL ASSERT(NPT(nr).LE.NPM,
     '                  '>>Too many nodes - increase NPM',ERROR,*9999)
                      NPNODE(0,nr)=NPNODE(0,nr)+1
                      NPNODE(NPNODE(0,nr),nr)=NPT(nr)
                      NODE=NPT(nr)
C GMH 8/1/97 Update cmgui link
                      CALL NODE_CHANGE(NODE,.FALSE.,ERROR,*9999)

                      IF(ITYP10(1).GE.2) THEN !If an angle is 2*pi set it to zero
                        IF(DABS(XS(2)-2.0D0*PI).LT.1.0D-5) XS(2)=0.0D0
                        IF(ITYP10(1).GE.3) THEN
                          IF(DABS(XS(3)-2.0D0*PI).LT.1.0D-5) XS(3)=0.0D0
                        ENDIF
                      ENDIF

                      NHP(NPT(nr),nr,nx)=
     '                  NHP(NPNODE(NPNODE(0,nr)-1,nr),nr,nx)
                      DO nj=1,NJ_LOC(0,0,nr)
                        NKJ(nj,NPT(nr))=
     '                    NKJ(nj,NPNODE(NPNODE(0,nr)-1,nr))
                      ENDDO
                      DO nhx=1,NHP(NPT(nr),nr,nx)
                        nh=NH_LOC(nhx,nx)
                        NKH(nh,NPT(nr),nc,nr)=
     '                   NKH(nh,NPNODE(NPNODE(0,nr)-1,nr),nc,nr)
                      ENDDO
                      DO nj=1,NJ_LOC(0,0,nr)
                        XP(1,1,nj,NODE)=XS(nj)
                      ENDDO
                    ENDIF
                    IF(DOP) THEN
                      WRITE(OP_STRING,'('' nn='',I1,'' EXIST='',L1,'
     '                  //''' NODE='',I4)') nn,EXIST,NODE
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                    NPNE(nn,NB9,ne)=NODE
                  ENDDO
                  DO nj=1,NJT
                    NBJ(nj,ne)=NB9
                    DO i=1,5
                      nn=NN5(i)
                      NKJE(1,nn,nj,ne)=1
                    ENDDO
                  ENDDO
C                  DO i=1,5
C                    nn=NN5(i)
C                    NKE(1,nn,NB9,ne)=1
C                  ENDDO
                  SE(1,NB9,ne)=1.0D0
                ENDDO
              ENDDO
            ENDIF

C KAT 14Sep98:  NPL contains the information needed for LINCAL
C            !Arc-length scaling requires NXI before SE calc'ed
C            CALL NENXI(IBT,INP,NBJ,NEELEM,NPNE,NXI,ERROR,*9999)

            CALL LINCAL(IBT,IDO,INP,0,NBJ,NEELEM,NEL,NENP,
     '        NKJE,NLL,NLLINE,NNL,NPL,NPNE,
     '        NPNODE,NRE,NVJE,NVJL,DL,SE,XP,ERROR,*9999)
            CALL FACCAL(IBT,IDO,INP,NBJ,NBJF,NEELEM,NFF,NFFACE,NKJE,
     '        NKEF,NLF,NNF,NNL,NPF,NPL,NPNE,NPNF,NRE,
     '        NVJE,NVJF,DF,PG,RG,SE,SF,WG,XA,XE,XG,XP,
     '        ERROR,*9999)
          ELSE IF(HANGING_NODE) THEN
C CS new 14/12/98
C MLB NOTE: This is not the best/correct way to do this.
C  something like REFINE_SETNODE should be used to create
C  the node properly including things like scale factors and
C  derivatives
            nr=NRLIST(1)
            NPNODE(0,nr)=NPNODE(0,nr)+1
            NPNODE(NPNODE(0,nr),nr)=NPT(0)+1
            NPT(0)=NPT(0)+1
            NPT(nr)=NPT(0)
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '        SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
              nb=NBJ(nj,ne)
              XS(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '          nb,1,XI,XE(1,nj))
            ENDDO
            NHP(NPT(nr),nr,nx)=
     '        NHP(NPNODE(NPNODE(0,nr)-1,nr),nr,nx)
            DO nj=1,NJ_LOC(0,0,nr)
              NKJ(nj,NPT(nr))=
     '          NKJ(nj,NPNODE(NPNODE(0,nr)-1,nr))
            ENDDO
            DO nhx=1,NHP(NPT(nr),nr,nx)
              nh=NH_LOC(nhx,nx)
              NKH(nh,NPT(nr),nc,nr)=
     '          NKH(nh,NPNODE(NPNODE(0,nr)-1,nr),nc,nr)
            ENDDO
            DO nj=1,NJ_LOC(0,0,nr)
              XP(1,1,nj,NPT(nr))=XS(nj)
              NVJP(nj,NPT(nr))=1
              !temp set derivatives to zero
              DO nk=2,NKJ(nj,NPT(nr))
                XP(nk,1,nj,NPT(nr))=0.0d0
              ENDDO
            ENDDO
            CALL INTERFACE(NP_INTERFACE,NPNODE,ERROR,*9999)
          ELSE IF(DELAUNAY_MESH)THEN
            nr=NRLIST(1) !only for elements in a single region
            XPNP_PTR=0
            XP_B_PTR=0
            XP_IB_PTR=0
c            BDRY_PTR=0
            NFNP_PTR=0
c            NXF_PTR=0
c            MEM_INIT=.FALSE.
            CALL ALLOCATE_MEMORY(3*NP_R_M,1,DPTYPE,XPNP_PTR,MEM_INIT,
     '        ERROR,*9999)
            CALL ALLOCATE_MEMORY(3*NE_R_M,1,DPTYPE,XP_B_PTR,MEM_INIT,
     '        ERROR,*9999)
            CALL ALLOCATE_MEMORY(3*NE_R_M,1,DPTYPE,XP_IB_PTR,MEM_INIT,
     '        ERROR,*9999)
c            CALL ALLOCATE_MEMORY(NE_R_M,1,INTTYPE,BDRY_PTR,MEM_INIT,
c     '        ERROR,*9999)
            CALL ALLOCATE_MEMORY(11*NP_R_M,1,INTTYPE,NFNP_PTR,MEM_INIT,
     '        ERROR,*9999)
c           CALL ALLOCATE_MEMORY(7*2*NFM,1,INTTYPE,NXF_PTR,MEM_INIT,
c     '        ERROR,*9999)
            
            CALL DELAUNAY_NODES(IBT,IDO,INP,NBJ,
     '        NDIVISION,NEELEM,NELIST,NFF,NFFACE,NFLIST,%VAL(NFNP_PTR),
     '        NKJ,NKJE,NNF,NPF,NPLIST,NPNE,NPNODE,
     '        NP_INTERFACE,nr,nr_target,nr_tree,NVCNP,NVJE,NVJP,
     '        NXI,del_radius,SE,XE,XP,%VAL(XP_B_PTR),%VAL(XP_IB_PTR),
     &        %VAL(XPNP_PTR),INTERNAL_NODES,REGULAR,ERROR,*9999)
            
            CALL FREE_MEMORY(XPNP_PTR,ERROR,*9999)
            CALL FREE_MEMORY(XP_B_PTR,ERROR,*9999)
            CALL FREE_MEMORY(XP_IB_PTR,ERROR,*9999)
c            CALL FREE_MEMORY(BDRY_PTR,ERROR,*9999)
            CALL FREE_MEMORY(NFNP_PTR,ERROR,*9999)
c            CALL FREE_MEMORY(NXF_PTR,ERROR,*9999)

          ELSE IF(.NOT.INTERPOLATED_NODES) THEN
!           !calculate nodes from data
            ARC_LENGTH=0.0d0
            DO nd=2,NDT
              ARC_INCREM=0.0d0
              DO nj=1,NJT
                ARC_INCREM=ARC_INCREM+(ZD(nj,nd)-ZD(nj,nd-1))**2
              ENDDO
              ARC_LENGTH=ARC_LENGTH+DSQRT(ARC_INCREM)
              ZD2(1,nd)=ARC_LENGTH
            ENDDO
            IF(DOP) THEN
              WRITE(OP_STRING,'('' Arc-length='',E12.4)') ARC_LENGTH
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            ARC_INCREM=ARC_LENGTH/DBLE(NPLIST(0)-1)
            ARC_LENGTH=0.0D0
            nd=1
            nr=NRLIST(1)
            DO nolist=1,NPLIST(0)
              np=NPLIST(nolist)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
              NPNODE(0,nr)=NPNODE(0,nr)+1
              NPNODE(NPNODE(0,nr),nr)=np
              IF(nolist.EQ.1) THEN
                DO nj=1,NJT
                  XP(1,1,nj,np)=ZD(nj,1)
                ENDDO
              ELSE IF(nolist.EQ.NPLIST(0)) THEN
                DO nj=1,NJT
                  XP(1,1,nj,np)=ZD(nj,NDT)
                ENDDO
              ELSE
                ARC_LENGTH=ARC_LENGTH+ARC_INCREM
                DO WHILE(ZD2(1,nd).LT.ARC_LENGTH)
                  nd=nd+1
                ENDDO
                XID=(ARC_LENGTH-ZD2(1,nd-1))/(ZD2(1,nd)-ZD2(1,nd-1))
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' np='',i4,'
     '              //''' arc_length='',e12.3,'' nd='',i6)')
     '              np,arc_length,nd
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                DO nj=1,NJT
                  XP(1,1,nj,np)=ZD(nj,nd-1)*(1.0D0-XID)+ZD(nj,nd)*XID
                  IF(DOP) THEN
                    WRITE(OP_STRING,
     '                '('' xp(1,'',i2,'',np)='',e12.3)') nj,
     '                xp(1,1,nj,np)
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDDO
              ENDIF
              DO nj=1,NJT
                NKJ(nj,np)=NO_DERIVS+1
              ENDDO
              IF(NO_DERIVS.GT.0) THEN
                IF(nolist.EQ.1) THEN
                  DO nj=1,NJT
                    XP(2,1,nj,np)=(ZD(nj,2)-ZD(nj,1))
     '                /(ZD2(1,2)-ZD2(1,1))
                  ENDDO
                ELSE IF(nolist.EQ.NPLIST(0)) THEN
                  DO nj=1,NJT
                    XP(2,1,nj,np)=(ZD(nj,NDT)-ZD(nj,NDT-1))
     '                /(ZD2(1,NDT)-ZD2(1,NDT-1))
                  ENDDO
                ELSE
                  DO nj=1,NJT
                    XP(2,1,nj,np)=(ZD(nj,nd)-ZD(nj,nd-1))
     '                /(ZD2(1,nd)-ZD2(1,nd-1))
                  ENDDO
                ENDIF
                IF(NO_DERIVS.GT.1) THEN
                  XP(3,1,2,np)=-XP(2,1,1,np)
                  XP(3,1,1,np)= XP(2,1,2,np)
                  DO nj=1,NJT
                    XP(4,1,nj,np)=0.0D0
                  ENDDO
                ENDIF
              ENDIF
            ENDDO
            NPT(nr)=NPNODE(0,nr)
            CALL ASSERT(NPT(nr).LE.NPM,'>>NPM too small',ERROR,
     '        *9999)
          ENDIF
        ENDIF !file io/mouse/calc.d

        CALL_NODE=.TRUE.

C KSB Adding node group:
        IF(MAKE_GROUP) THEN
          STRING=GROUP_NAME
          i=0
          nr=NRLIST(1)
          DO nonode=nonode_start+1,NPNODE(0,nr)
            i=i+1
            NPLIST(i)=NPNODE(nonode,nr)
          ENDDO !noelem
          NPLIST(0)=i
          CALL GRNODE_SUB(NPLIST,STRING,.TRUE.,ERROR,*9999)          
        ENDIF !MAKE_GROUP
      ENDIF

      IF(FILEFORMAT(1:6).EQ.'BINARY') THEN
        CALL IONODE(IOFILE1,NKJ,NPNODE,NRLIST,NVJP,XP,
     '    'CLOSE',FILEFORMAT,FILE,' ',ERROR,*9998)
      ENDIF

      CALL EXITS('DENODS')
      RETURN
 9999 CALL ERRORS('DENODS',ERROR)
 9998 CALL EXITS('DENODS')
      RETURN 1
      END

