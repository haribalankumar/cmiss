      SUBROUTINE DEXI(FD,IBT,IDO,INP,LD,LD_NP,LN,NBH,NBHF,NBJ,NBJF,
     '  NDDL,NDLT,NDP,NEELEM,NELIST,NENP,NEP,NFF,NFFACE,
     '  NFLIST,NHE,NKEF,NKHE,NKJE,NNB,NNF,NPF,NPNE,NPNF,NPNODE,NRE,
     '  NRLIST,NRLIST2,NVHE,NVJE,NVJF,NW,NXI,Z_CONT_LIST,
     '  CURVCORRECT,SE,SF,SQ,XA,XAB,XE,XID,XIP,XP,ZA,
     &  ZD,ZDD,ZE,ZP,STRING,ERROR,*)

C#### Subroutine: DEXI
C###  Description:
C###    DEXI defines data xi points.
C GMH 30/9/95 There were some vast tracts of commented out crap here,
C     so those regions have been relegated to *** ARCHIVE ***
C     (if you REALLY want to look at them...)
C CS  March 98 Split of some code into new routines
C     DEXI_ORTHOG,DEXI_1D,DEXI_NONLIN,DEXI_LINEAR,DEXI_EXISTING,
C     for maintainence and efficiency.
c## JWF 11/12/01 Added extra parameters NBJF,NKEF,NNF,NPNF,NVJF,SF
c## for face projections. Also added new logical CLOSEST_FACE.
      

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'data00.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'four00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ioda00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'map000.cmn'
      INCLUDE 'nonl00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER FD(NDM),IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),LD(NDM),LD_NP(NDM),LN(0:NEM),NBH(NHM,NCM,NEM),
     '  NBHF(NHM,NCM,NFM),NBJ(NJM,NEM),NBJF(NJM,NFM),NDDL(NEM,NDEM),
     '  NDLT(NEM),NDP(NDM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NENP(NPM,0:NEPM,0:NRM),NEP(NPM),
     &  NFF(6,NEM),NFFACE(0:NF_R_M,NRM),NFLIST(0:NFM),NHE(NEM,NXM),
     '  NKEF(0:4,16,6,NBFM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),
     '  NPNODE(0:NP_R_M,0:NRM),NRE(NEM),NRLIST(0:NRM),
     '  NRLIST2(0:NRM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJF(NNM,NBFM,NJM),NW(NEM,3,NXM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),Z_CONT_LIST(NDM,2,7)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),SE(NSM,NBFM,NEM),SF(NSM,NBFM),
     '  SQ(NDM),XA(NAM,NJM,NEM),XAB(NORM,NEM),XE(NSM,NJM),XID(NIM,NDM),
     '  XIP(NIM,NPM),XP(NKM,NVM,NJM,NPM),
     &  ZA(NAM,NHM,NCM,NEM),ZD(NJM,NDM),ZDD(NJM,NDM),
     &  ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(*)

!     Local Variables
      INTEGER A_LIST(NNM),CONTACT_STATUS,DEFLIST(3),direction,i,IBEG,
     &  IBEG1,IBEG2,
     '  icol,IEND,IEND1,IEND2,IFROMC,INT_TEMP(1),IT,j,n1list,N3CO,nb,
     '  nbface,NCLUSTER,ND0,ND1,ND2,nd,NDFXI(2),nd_increment,ne,nef,
     '  NELIST2(0:1000),nf,nf_global,ni,NILIST(1),NITB,nj,nj1,nj2,
     '  nj_elem,nj_num,NKJF(NKM,NNM,NJM),NKTB,nn,noelem,noface,nolist,
     '  nonode,nonr,NOXIPT(3),np,nplist,NPLIST2(0:NPM),nr,nr1,nr2,
     '  ns1,ns2,ns3,ns4,NSEED,nt,NTILDEF,NTIL,NTIMES,NTLIST,nv,nx,
     '  LDTEMP,START_ELEMENT,tempcount,nsearch,CONT_WEIGHT
      REAL*8 A1,A2,AB(3),ALFA,B1,B2,BDOTC,BETA,BMAG(3),
     '  CMAG(3),CVEC(3,3),C1,C2,coord_value,
     '  D1,D2,DELTA,DENOM1,DENOM2,DIFF,DIST,FACTOR_TOL,GAMA,LIMIT,
     &  NODE_POS(3),PXI,REAL_TEMP(1),
     '  RFROMC,SLOPE,SQMAX,SQND,THETAMIN,THETAMAX,
     '  TMIN(2),TMAX(2),TREF,X0,X1,X2,XD(3),XI(3),XI3OFF,
     '  XI_1,XI_2,XI_3,Z1(4),Z2(3),Z3(3),Z4(3)
      CHARACTER CHAR1*20,CIW*1,FILE*(MXCH),STATUS*3,TYPE*20
      LOGICAL ALL_REGIONS,CALCU,CBBREV,CENTROID,CLOSEST,
     '  CLOSEST_FACE,FRUSTUM,CUBIC,DEFORM,EXCLUDE,EXTERNAL,EXTRAPOLATE,
     '  FAST,FILIO,FINISHED,FOUND,GENER,GRPGRID,INLIST,INTERNAL,
     '  MOUSE,NEW,NODE_ASSOC,NUM_FIELD,ORDER_CLOSEST,ORTHOG,PASS2,
     '  QUADRATIC,REFERENCE,SUBDIV,SET_XI_1,SET_XI_2,SET_XI_3,SPECIFY,
     '  TAG2D,USE_LOOSE_TOL,CFORTHO,CFCROSS,NEW_CONT_PTS,BOUNDARY_CONT

      CALL ENTERS('DEXI',*9999)


C ??? Why is nx used here ???
      nx=1 ! temporary cpb 22/11/94
      nv=1 !temporary
      GRPGRID=.FALSE.

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)
        WRITE(CIW,'(I1)') 2*NJT-3+IMAP

C---------------------------------------------------------------------

C#### Command: FEM define xi;r/w<;FILENAME[$current]><;(PATH
C###  /example)[$current]> Description:
C###    Defines the xi postions of data points. Positions are read from
C###    or written to the a file FILENAME (extension .ipxi) in the
C###    directory specified by PATH with $current specifing the current
C###    default file.
C###  Parameter:      <total NIT#[NJT-1]>
C###    Specifies the number of xi directions
C###  Parameter:      <reference>
C###    Specifies the output to be labeled using the original data
C###    point numbers.
C###  Parameter:      <coupled_nodes>
C###    Specifies the output to include single nodes associated with
C###    data points. Only valid following use of data points to generate
C###    a volume-filling tree.
C###  Parameter:      <element (#s/all)[all]>
C###    Specifies the element numbers containing data points.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###

        OP_STRING(1)=STRING(1:IEND)//';r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<total NIT#[NJT-1]>'
        OP_STRING(3)=BLANK(1:15)//'<reference>'
        OP_STRING(4)=BLANK(1:15)//'<coupled_nodes>'
        OP_STRING(5)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(6)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define xi;r/w<;FILENAME[$current]><;(PATH
C###  /example)[$current] closest_face external/orthogonal/faces>
C###  Description:
C###    Defines the xi postions of data points projected on to the faces
C###    of volume elements. Positions are read from or written to the
C###    file FILENAME (extension .ipxi) in the directory specified by
C###    PATH with $current specifing the current default file. Each data
C###    point has the element no. associated with it and the local face
C###    no. to which data pt. was projected. If the orthogonal option is
C###    specified then only data points having an orthogonal projection
C###    are included. 
C###  Parameter:      <reference>
C###    Specifies the output to be labeled using the original data
C###    point numbers.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###

        OP_STRING(1)=STRING(1:IEND)//';r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']> '
     '    //'closest_face'
        OP_STRING(2)=BLANK(1:15)//'external/orthogonal/faces '
     &    //'(FACE_GROUP_NAME)'
        OP_STRING(3)=BLANK(1:15)//'<reference>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define xi;r/w<;FILENAME[$current]><;(PATH
C###  /example)[$current] nodes> Description:
C###    Defines the xi postions of data points. Positions are read from
C###    or written to the a file FILENAME (extension .ipxi) in the
C###    directory specified by PATH with $current specifing the current
C###    default file.
C###  Parameter:      <total NIT#[NJT-1]>
C###    Specifies the number of xi directions
C###  Parameter:      <element (#s/all)[all]>
C###    Specifies the element numbers containing data points.
C###  Parameter:      <of REGION# (#s)[2]>
C###    Specify the region numbers which the nodes are in.
C###  Parameter:      <in REGION# (#s)[1]>
C###    Specify the region number within which the nodes are contained.

        OP_STRING(1)=STRING(1:IEND)//';r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']> nodes'
        OP_STRING(2)=BLANK(1:15)//'<total NIT#[NJT-1]>'
        OP_STRING(3)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define xi;c
C###  Description:
C###    Calculates the xi position of a data points within an element.
C###  Parameter:      <(new/old)[new]>
C###  Parameter:      <(linear/quadratic/cubic/orthogonal/closest
C###  /closest_face/frustum)[linear]>
C###    Specifies the type of interpolation or projection used to
C###    determine the xi position.  With 'orthogonal' the error
C###    projection must be orthogonal to the element.  With 'closest'
C###    simply the nearest point is found. For 'closest_face'
C###    the projection is onto faces of 3D elements and the closest is
C###    determined by distance to the center of the face. 'frustum'
C###    uses a frustum culling technique to determine which face to
C###    project the data point onto.  If the option 'external' is
C###    specified with closest_face or frustum then only external
C###    faces are considered.
C###  Parameter:      <extrapolate>
C###    Extrapolates xi positions outside the range of 0 to 1
C###  Parameter:      <element (#s/all)[all]>
C###    Specifies the element number within which to calculate the
C###    xi positions
C###  Parameter:      <faces (#s/all)[all]>
C###    With the closest_face and frustum options only the specified
C###    faces are considered.
C###  Parameter:      <(undeformed/deformed)[undeformed] DEFLIST[xyz]>
C###  Parameter:      <max DIST#[10000.0]>
C###    Specifies the maximum distance within which to calculate a
C###    projection.
C###  Parameter:      <accept XI3_OFFSET#[0.0]>
C###    Specifies acceptable xi3 offset outside the range of [0,1]
C###  Parameter:      <xi_3 XI_3#[0.0]>
C###    Specifies the value of xi3 in xi calcalculation
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <seed_points #[1]>
C###    Specifies the number of seed points (in each xi direction)
C###     to use for each element to find the initial starting position
C###     to begin non-linear search for the closest point.
C###     One seed point has a search location of [0.5,0.5] while two seed
C###     points has searches from xi=[1/3,2/3].        
C###  Parameter:      <search_start #[2]>
C###    Specifies the number of faces to use for making an initial guess.
C###  Parameter:      <cross_boundaries>
C###    Use the element boundary crossing algorithm. Note that search_start
C###    should set to 1 for this to be most efficient.

        OP_STRING(1)=STRING(1:IEND)//';c'
        OP_STRING(2)=BLANK(1:15)//'<(new/old)[new]>'
        OP_STRING(3)=BLANK(1:15)
     &    //'<(linear/quadr/cubic/orthogonal/closest/'
        OP_STRING(4)=BLANK(1:15)
     &    //'closest_face(external/orthogonal)/frustum(external)/'
        OP_STRING(5)=BLANK(1:15)
     &    //'contain [linear]>'
        OP_STRING(6)=BLANK(1:15)   
     &    //'<(undeform/deform)[undeform] DEFLIST[xyz]>'
        OP_STRING(7)=BLANK(1:15)//'<extrapolate>'
        OP_STRING(8)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(9)=BLANK(1:15)//'<faces (#s/all)[all]>'
        OP_STRING(10)=BLANK(1:15)//'<max DIST#[10000.0]>'
        OP_STRING(11)=BLANK(1:15)//'<accept XI3_OFFSET#[0.0]>'
        OP_STRING(12)=BLANK(1:15)//'<xi_3 XI_3#[0.0]>'
        OP_STRING(13)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(14)=BLANK(1:15)//'<seed_points #[1]>'
        OP_STRING(15)=BLANK(1:15)//'<search_start #[2]>'
        OP_STRING(16)=BLANK(1:15)//'<cross_boundaries>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C GMH 30/9/95 Allowing the user to force the projection process
C       to start from a specific position.
C       The documentation below is wrong - need to show that they
C       can specify - which data point(s)
C       - start element number
C       - start xi position
C       As far as implementation goes, we want to do exactly the
C       same things as for de xi;c old, but restart from a new
C       position (user-specified)
C GMH 30/9/95 ??? Should <from first,last> be a global option?

C#### Command: FEM define xi;c specify
C###  Description:
C###    This command allows the user to force the projection process
C###    to start from a specific xi position.
C###  Parameter:      <from (DATA_POINT#/all)[all]>
C###    Specify the numbers of data points
C###  Parameter:      <element#)[1]>
C###    Specify the start element number
C###  Parameter:      <xi_1 XI_1#[0.0]>
C###    Specify the start xi1 position
C###  Parameter:      <xi_2 XI_2#[0.0]>
C###    Specify the start xi2 position
C###  Parameter:      <xi_3 XI_3#[0.0]>
C###    Specify the start xi3 position
C###

        OP_STRING(1)=STRING(1:IEND)//';c specify'
        OP_STRING(2)=BLANK(1:15)//'<from (DATA_POINT#/all)[all]>'
        OP_STRING(3)=BLANK(1:15)
     '    //'<(linear/quadr/cubic/orthogonal/closest)[linear]>'
        OP_STRING(4)=BLANK(1:15)
     '    //'closest_face(external/orthogonal)/contain [linear]'
        OP_STRING(5)=BLANK(1:15)//'<extrapolate>'
        OP_STRING(6)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(7)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(8)=BLANK(1:15)
     '    //'<(undeform/deform)[undeform] DEFLIST[xyz]>'
        OP_STRING(9)=BLANK(1:15)//'<max DIST#[10000.0]>'
        OP_STRING(10)=BLANK(1:15)//'<accept XI3_OFFSET#[0.0]>'
        OP_STRING(11)=BLANK(1:15)//'<xi_1 XI_1#[0.0]>'
        OP_STRING(12)=BLANK(1:15)//'<xi_2 XI_2#[0.0]>'
        OP_STRING(13)=BLANK(1:15)//'<xi_3 XI_3#[0.0]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define xi;c centroid
C###  Description:
C###    Data are projected onto elements by radial lines passing thru
C###    the centroid of each data scan plane (used by David Rahdert).
C###    See examples 2183 & 2184.
C###  Parameter:      <element (#s/all)[all]>
C###    Specify the element numbers to used. The "all" keyword will
C###    use all currently defined elements in the given regions.
C###  Parameter:      <cluster NCLUSTER#[72]>
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.

        OP_STRING(1)=STRING(1:IEND)//';c'
        OP_STRING(2)=BLANK(1:15)//'centroid'
        OP_STRING(3)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<cluster NCLUSTER#[72]>'
        OP_STRING(5)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define xi;c
C###  Description:
C###    Calculates the xi position of data points within an element.
C###  Parameter:      <sub DXI1#,DXI2#[10,10]>
C###  Parameter:      <element (#s/all)[all]>
C###    Specify the element numbers to used. The "all" keyword will
C###    use all currently defined elements in the given regions.

        OP_STRING(1)=STRING(1:IEND)//';c'
        OP_STRING(2)=BLANK(1:15)//'<sub DXI1#,DXI2#[10,10]>'
        OP_STRING(3)=BLANK(1:15)//'<element (#s/all)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define xi;c tag2d/contour
C###  Description:
C###    Not used? Should be archived.
C###  Parameter:      <max DIST#[10000.0]>
C###  Parameter:      <surf INDEX#[1]>
        OP_STRING(1)=STRING(1:IEND)//';c'
        OP_STRING(2)=BLANK(1:15)//'tag2d/contour'
        OP_STRING(3)=BLANK(1:15)//'<max DIST#[10000.0]>'
        OP_STRING(4)=BLANK(1:15)//'<surf INDEX#[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define xi;c subdivision
C###  Description:
C###    Calculate approxiamte xi possitions using subdivision of
C###    elements.
C###  Parameter:      <numsub nx,ny,nz[10,10,10]>
        OP_STRING(1)=STRING(1:IEND)//';c'
        OP_STRING(2)=BLANK(1:15)//'subdivision'
        OP_STRING(3)=BLANK(1:15)//'<numsub nx,ny,nz[10,10,10]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define xi;c beads
C###  Description:
C###    Calculate xi-coordinates of projections of bead positions on to
C###    xi_3 lines of 3D elements using normalized lengths in each
C###    column. It is assumed that data points 1,2,3 are epicardial
C###    beads numbered clockwise, that data points 4,5,6 are bottom
C###    beads of cols 1,2,3, and that elements numbers increase from
C###    epicardium to endocardium.
C###  Parameter:      <element (#s/all)[all]>
C###    Specify the element numbers to used. The "all" keyword will
C###    use all currently defined elements in the given regions.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.


        OP_STRING(1)=STRING(1:IEND)//';c beads'
        OP_STRING(3)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define xi;c nodes
C###  Description:
C###    Calculates the xi position nodes in REGION#[2] contained
C###    in a host mesh in REGION#[1].
C###  Parameter:      <(contain/closest)[contain]>
C###    Contain specifies that the xi position to be calculated in 
C###    contained within a host mesh. Closest finds the xi positions 
C###    for the nearest point.
C###  Parameter:      <of REGION#[2]> <in REGION#[1]>
C###    Specifies the region of the host mesh REGION#[1] and the region
C###    containing the nodes to be calculated REGION#[2].
C###  Parameter:      <loose_tolerance>
C###    Specify that a loose tolerance will be used for convergence
C###    testing of nodal projections, where the default is 10
C###    *converg_tol.
C###  Parameter:      <by FACTOR_TOL[10]>
C###    Specify a factor by which the convergence tolerance will be
C###    relaxed (must have also specified loose_tolerance).
C###  Parameter:      <order_closest>
C###    Orders the element search so that elements surrounding the
C###    closest node in host region are searched first.

        OP_STRING(1)=STRING(1:IEND)//';c nodes'
        OP_STRING(2)=BLANK(1:15)//'<(contain/closest)[contain]>'
        OP_STRING(3)=BLANK(1:15)//'<of REGION#[2]> <in REGION#[1]>'
        OP_STRING(4)=BLANK(1:15)//'<loose_tolerance>'
        OP_STRING(5)=BLANK(1:15)//'<by FACTOR_TOL[10]>'
        OP_STRING(6)=BLANK(1:15)//'<order_closest>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define xi;c contact_points
C###  Description:
C###    Calculates xi-coordinates of contact points on specified faces.
C###    These are also the integration points for integrating the
C###    contact force over the surface, so the number of integration
C###    points per xi direction must be specified.
C###    Eg: faces 1,3 points 5, will place 25 equally spaced contact points
C###    on each face 1 and 3. The 'Add' option allows the user to define further
C###    contact points in the same region of a different type (e.g. friction). 
C###  Parameter:      <faces (#s)>
C###  Parameter:      <points (#)>
C###  Parameter:      <by (#)>
C###  Parameter:      <ADD>
C###  Parameter:      <boundary>
        OP_STRING(1)=STRING(1:IEND)//';c contact_points'
        OP_STRING(2)=BLANK(1:15)//'<faces (#s)>'
        OP_STRING(3)=BLANK(1:15)//'<points (#)>'
        OP_STRING(4)=BLANK(1:15)//'<by (#)>'
        OP_STRING(5)=BLANK(1:15)//'<ADD>'
        OP_STRING(6)=BLANK(1:15)//'<direction (XMAX,YMAX,ZMAX)>'
        OP_STRING(7)=BLANK(1:15)//'<limit>'
        OP_STRING(8)=BLANK(1:15)//'<boundary>'

        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define xi;c contact_tied
C###  Description:
C###    Same as contact_points option but defines points as
C###    tied contact points. Tied contact also offers the option
C###    XI_END. If this is specified instead of faces then the
C###    external ends of muscles are covered in tie points. This 
C###    currently only works if the muscle has xi 1 in the fibre
C###    direction and the xi system follows a RH rule.
C###  Parameter:      <faces (#s)>
C###  Parameter:      <XI_END>
C###  Parameter:      <points (#)>
C###  Parameter:      <by (#)>
C###  Parameter:      <ADD>
C###  Parameter:      <boundary>
        OP_STRING(1)=STRING(1:IEND)//';c contact_tied'
        OP_STRING(2)=BLANK(1:15)//'<faces (#s)>'
        OP_STRING(3)=BLANK(1:15)//'<XI_END>'
        OP_STRING(4)=BLANK(1:15)//'<points (#)>'
        OP_STRING(5)=BLANK(1:15)//'<by (#)>'
        OP_STRING(6)=BLANK(1:15)//'<ADD>'
        OP_STRING(7)=BLANK(1:15)//'<boundary>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define xi;c contact_friction
C###  Description:
C###    Same as contact_points option but defines points as 
C###    including friction (resistance in shear direction). 
C###  Parameter:      <faces (#s)>
C###  Parameter:      <points (#)>
C###  Parameter:      <by (#)>
C###  Parameter:      <ADD>
C###  Parameter:      <boundary>
        OP_STRING(1)=STRING(1:IEND)//';c contact_friction'
        OP_STRING(2)=BLANK(1:15)//'<faces (#s)>'
        OP_STRING(3)=BLANK(1:15)//'<points (#)>'
        OP_STRING(4)=BLANK(1:15)//'<by (#)>'
        OP_STRING(5)=BLANK(1:15)//'<ADD>'
        OP_STRING(6)=BLANK(1:15)//'<direction (XMAX,YMAX,ZMAX)>'
        OP_STRING(7)=BLANK(1:15)//'<limit>'
        OP_STRING(8)=BLANK(1:15)//'<boundary>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEXI',ERROR,*9999)
      ELSE
C       nc=1 ! Temporary MPN 12-Nov-94

C!!! LKC Why to we need KTYP11 here? IF we are writing out xi positions
C then the method used to calculate the xi positions shouldn't really
C matter should it? Actually why do we need a different file format for
C projections onto faces or element????
C        
C At present if you are writing out the xi positions for a face
C projection you must specify a method which resets KTYP12 back to 2.
        
        KTYP11=1  ! setting KTYP11 for element fitting

        CALL PARSE_QUALIFIERS(' CDLMPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

C LKC 2-JUN-1998 assert
        CALL ASSERT(CALL_ELEM,'>>Define elements first',ERROR,*9999)

        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

        IF(CBBREV(CO,'BEADS',1,noco+1,NTCO,N3CO)) THEN
          TYPE='BEADS'
        ELSE IF(CBBREV(CO,'CONTACT_POINTS',1,noco+1,NTCO,N3CO)) THEN
          TYPE='CONTACT'
          CONTACT_STATUS=1
        ELSE IF(CBBREV(CO,'CONTACT_TIED',1,noco+1,NTCO,N3CO)) THEN
          TYPE='CONTACT'
          CONTACT_STATUS=2
        ELSE IF(CBBREV(CO,'CONTACT_FRICTION',1,noco+1,NTCO,N3CO)) THEN
          TYPE='CONTACT'
          CONTACT_STATUS=3    
        ELSE IF(CBBREV(CO,'NODES',1,noco+1,NTCO,N3CO)) THEN
          TYPE='NODES'
C??? GDR 8Mar05 what is the following for?
          IF(CBBREV(CO,'CONTAIN',2,noco+1,NTCO,N3CO)) THEN
            INTERNAL=.TRUE.
          ELSE
            INTERNAL=.TRUE.
          ENDIF
          
C LKC 12-MAR-1999  Generalise for multiple regions
C         IF(CBBREV(CO,'OF',2,noco+1,NTCO,N3CO)) THEN
C         nr2=IFROMC(CO(N3CO+1))
C         C LKC 16-FEB-1999
C         CALL ASSERT(NEELEM(0,nr2).GT.0,
C         '        '>> No elements defined in region',ERROR,*9999)
C         ELSE
C         nr2=2
C         ENDIF

          IF(CBBREV(CO,'LOOSE_TOLERANCE',3,noco+1,NTCO,N3CO)) THEN
            USE_LOOSE_TOL=.TRUE.
            IF(CBBREV(CO,'BY',2,noco+1,NTCO,N3CO)) THEN
              FACTOR_TOL=RFROMC(CO(N3CO+1))
            ELSE
              FACTOR_TOL=10.d0
            ENDIF
          ELSE
            USE_LOOSE_TOL=.FALSE.
            FACTOR_TOL=0
          ENDIF

          IF(CBBREV(CO,'ORDER_CLOSEST',3,noco+1,NTCO,N3CO)) THEN
            ORDER_CLOSEST=.TRUE.
          ELSE
            ORDER_CLOSEST=.FALSE.
          ENDIF

          IF(CBBREV(CO,'OF',2,noco+1,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),NRM,NRLIST2(0),NRLIST2(1),
     '        ERROR,*9999)
            CALL ASSERT(NRLIST2(0).GT.0,
     '        '>> No regions entered for OF ',ERROR,*9999)
            DO nr2=1,NRLIST2(0)
c              CALL ASSERT(NEELEM(0,NRLIST2(nr2)).GT.0,
c             '          '>> No elements defined in region',ERROR,*9999)
              IF(NEELEM(0,NRLIST2(nr2)).EQ.0) THEN
                WRITE(OP_STRING,
     &            '(''WARNING: no elements defined in region'')')
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO
          ELSE
            NRLIST2(0)=1
            NRLIST2(1)=2
            CALL ASSERT(NEELEM(0,NRLIST2(1)).GT.0,
     '        '>> No elements defined in region 2',ERROR,*9999)
          ENDIF

          IF(CBBREV(CO,'IN',2,noco+1,NTCO,N3CO)) THEN
            nr1=IFROMC(CO(N3CO+1))
            NRLIST(0)=1
            NRLIST(1)=nr1
          ELSE
            NRLIST(0)=1
            NRLIST(1)=1
C LKC 11-FEB-2003 initialise nr1 as well            
            nr1=1
          ENDIF
          IF(CBBREV(CO,'CLOSEST',2,noco+1,NTCO,N3CO)) THEN
            CLOSEST=.TRUE.
          ELSE
            CLOSEST=.FALSE.
          ENDIF
          IF(CBBREV(CO,'NEW',1,noco+1,NTCO,N3CO)) THEN
            NEW=.TRUE.
          ELSE IF(CBBREV(CO,'OLD',2,noco+1,NTCO,N3CO)) THEN
            NEW=.FALSE.
          ELSE
            NEW=.TRUE.
          ENDIF

          IF(CBBREV(CO,'STORE_FIELD',3,noco+1,NTCO,N3CO))THEN
            NUM_FIELD=.TRUE.
            nj_num=NEJ_LOC(IFROMC(CO(N3CO+1)),nr1)
            nr=NRLIST2(1)
            nj_elem=NJ_LOC(NJL_FIEL,IFROMC(CO(N3CO+2)),nr)
            CALL ASSERT(nj_num.GT.0,'>>Non-existant field',ERROR,*9999)
            CALL ASSERT(nj_elem.GT.0,'>>Non-existant field',ERROR,*9999)
            IF(CBBREV(CO,'TERMINAL_NODES',4,noco+1,NTCO,N3CO))THEN
              CALL PARSILG(NPLIST2,NPM,'NODES','TERMINAL_NODES',
     '          ERROR,*9999)
              CALL ASSERT(NPLIST2(0).GT.0,'>>Node group not found',
     '          ERROR,*9999)
            ELSE
              CALL PARSE_NODES(NPNODE,NPLIST2,noco,NRLIST,NTCO,CO,
     '          ERROR,*9999)
            ENDIF
          ELSE
            NUM_FIELD=.FALSE.
          ENDIF


          ELSE IF(CBBREV(CO,'VERTEX',2,noco+1,NTCO,N3CO)) THEN
          TYPE='VERTEX'
        ELSE
          TYPE='XIPOS'

C LKC 18-FEB-1999
          IF(NEELEM(0,NRLIST(1)).LE.0) THEN
            OP_STRING(1)='>> WARNING: No elements in this region'
            CALL WRITES(IOER,OP_STRING,ERROR,*9999)
          ENDIF

          IF(CBBREV(CO,'REFERENCE',2,noco+1,NTCO,N3CO)) THEN
            REFERENCE=.TRUE.
          ELSE
            REFERENCE=.FALSE.
          ENDIF

          IF(CBBREV(CO,'COUPLED_NODES',4,noco+1,NTCO,N3CO)) THEN
            NODE_ASSOC=.TRUE.
c            IF(IOTYPE.EQ.3) CALL ASSERT(CALL_TREE,
c     &        '>>Generate tree first',ERROR,*9999)
          ELSE
            NODE_ASSOC=.FALSE.
          ENDIF

            CUBIC=.FALSE.
            QUADRATIC=.FALSE.
            ORTHOG=.FALSE.
            CLOSEST=.FALSE.
            CLOSEST_FACE=.FALSE.
            FAST=.FALSE.
            FRUSTUM=.FALSE.
            INTERNAL=.FALSE.
          SUBDIV=.FALSE.
          IF(CBBREV(CO,'LINEAR',1,noco+1,NTCO,N3CO)) THEN
C GDR 20Mar2005 Changed order of parsing so that orthogonal can be
C   used as an option for closest_face          
          ELSE IF(CBBREV(CO,'CLOSEST_FACE',8,noco+1,NTCO,N3CO)) THEN
            CLOSEST_FACE=.TRUE.
            INTERNAL=.FALSE.
            KTYP11=2     ! face fitting
          ELSE IF(CBBREV(CO,'ORTHOGONAL',2,noco+1,NTCO,N3CO)) THEN
            ORTHOG=.TRUE.
          ELSE IF(CBBREV(CO,'CLOSEST',2,noco+1,NTCO,N3CO)) THEN
            CLOSEST=.TRUE.
          ELSE IF(CBBREV(CO,'FRUSTUM',2,noco+1,NTCO,N3CO)) THEN
            FRUSTUM=.TRUE.
            INTERNAL=.FALSE.
            KTYP11=2     ! face fitting
          ELSE IF(CBBREV(CO,'CUBIC',2,noco+1,NTCO,N3CO)) THEN
            CUBIC=.TRUE.
          ELSE IF(CBBREV(CO,'QUADRATIC',2,noco+1,NTCO,N3CO)) THEN
            QUADRATIC=.TRUE.
          ELSE IF(CBBREV(CO,'CONTAIN',2,noco+1,NTCO,N3CO)) THEN
            INTERNAL=.TRUE.
            IF(CBBREV(CO,'FAST',3,noco+1,NTCO,N3CO)) FAST=.TRUE.
          ELSE IF(CBBREV(CO,'SUBDIVISION',2,noco+1,NTCO,N3CO)) THEN
            SUBDIV=.TRUE.
            IF(CBBREV(CO,'NUMSUB',3,noco+1,NTCO,N3CO)) THEN
              CALL PARSIL(CO(N3CO+1),3,tempcount,NOXIPT,ERROR,*9999)
          ELSE
              NOXIPT(1)=10
              NOXIPT(2)=10
              NOXIPT(3)=10
          ENDIF
          ELSE
          ENDIF

          IF(CBBREV(CO,'NEW',1,noco+1,NTCO,N3CO)) THEN
            NEW=.TRUE.
          ELSE IF(CBBREV(CO,'OLD',2,noco+1,NTCO,N3CO)) THEN
            NEW=.FALSE.
          ELSE
            NEW=.TRUE.
          ENDIF

C not used
C         IF(CBBREV(CO,'NODES',2,noco+1,NTCO,N3CO)) THEN
C         NODE=.TRUE.
C         ELSE
C         NODE=.FALSE.
C         ENDIF

          IF(ADD) THEN
            !calculate xi for OLD to NDT
            ND0=NDTOLD+1
            ND1=NDT
          ELSE
            ND0=1
            ND1=NDT
          ENDIF

C GMH 30/95 Check for 'specify' - and chain through as .NOT.NEW
          IF(CBBREV(CO,'SPECIFY',2,noco+1,NTCO,N3CO)) THEN
            SPECIFY=.TRUE.
            NEW=.FALSE.
            ! Which data point to use
            IF(CBBREV(CO,'FROM',2,noco+1,NTCO,N3CO)) THEN
              !Get the data point number

C LKC pass in an array
C             CALL PARSIL(CO(N3CO+1),1,NTIL,ND0,ERROR,*9999)
C             ND1=ND0 !only do one point
              CALL PARSIL(CO(N3CO+1),1,NTIL,INT_TEMP(1),ERROR,*9999)
              ND1=INT_TEMP(1)
            ENDIF
            ! Starting element number
            IF(CBBREV(CO,'ELEMENT',2,noco+1,NTCO,N3CO)) THEN
              !Get the element number
C LKC pass in an array
C             CALL PARSIL(CO(N3CO+1),1,NTIL,START_ELEMENT,
C             '          ERROR,*9999)
              CALL PARSIL(CO(N3CO+1),1,NTIL,INT_TEMP,ERROR,*9999)
              START_ELEMENT=INT_TEMP(1)
            ELSE
              START_ELEMENT=1
            ENDIF
            ! Get the xi positions
            IF(CBBREV(CO,'XI_1',4,NOCO+1,NTCO,N3CO)) THEN
              SET_XI_1 = .TRUE.

C LKC 24-APR-1998 need to pass in an array
C             CALL PARSRL(CO(N3CO+1),1,NTIL,XI_1,ERROR,*9999)
              CALL PARSRL(CO(N3CO+1),1,NTIL,XI,ERROR,*9999)
              XI_1=XI(1)
            ELSE
              SET_XI_1 = .FALSE.
              XI_1=0.0d0
            ENDIF
            IF(CBBREV(CO,'XI_2',4,NOCO+1,NTCO,N3CO)) THEN
              SET_XI_2 = .TRUE.
C LKC 24-APR-1998 need to pass in an array
C             CALL PARSRL(CO(N3CO+1),1,NTIL,XI_2,ERROR,*9999)
              CALL PARSRL(CO(N3CO+1),1,NTIL,XI,ERROR,*9999)
              XI_2=XI(1)
            ELSE
              SET_XI_2 = .FALSE.
              XI_2=0.0d0
            ENDIF
          ELSE
            SPECIFY=.FALSE.
          ENDIF

          IF(CBBREV(CO,'DEFORM',3,NOCO+1,NTCO,N3CO)) THEN
            DEFORM=.TRUE.
            IF(NTCO.GE.N3CO+1)THEN
              CALL PARSIL(CO(N3CO+1),3,NTILDEF,DEFLIST,ERROR,*9999)
              IF(DOP)THEN
                WRITE(*,*)' DEFLIST=',DEFLIST
              ENDIF
            ELSE
              NTILDEF=3
              DO i=1,3
                DEFLIST(i)=i
              ENDDO
            ENDIF
          ELSE
            DEFORM=.FALSE.
          ENDIF

          IF(CBBREV(CO,'MAX',3,NOCO+1,NTCO,N3CO)) THEN
C MPN 7-NOV-2014: implementing MAX
C LKC 27-APR-98 Need to pass in array
C CS  11/5/98 SQMAX set but not used. LKC to look into/implement this
            CALL PARSRL(CO(N3CO+1),1,NTIL,REAL_TEMP,ERROR,*9999)
            SQMAX=REAL_TEMP(1)*REAL_TEMP(1)
C            CALL ASSERT(SQMAX.LT.(SQMAX-1.0d0),'>>Warning MAX qualifier'
C     '        //'not implemented',ERROR,*9999)
          ELSE
            SQMAX = 1.d4*1.d4
          ENDIF

          IF(CBBREV(CO,'ACCEPT',3,NOCO+1,NTCO,N3CO)) THEN
C LKC 29-APR-98 Need to pass in array
C           CALL PARSRL(CO(N3CO+1),1,NTIL,XI3OFF,ERROR,*9999)
            CALL PARSRL(CO(N3CO+1),1,NTIL,REAL_TEMP,ERROR,*9999)
            XI3OFF=REAL_TEMP(1)
          ELSE
            XI3OFF = 0.d0
          ENDIF
          XI3MIN = -XI3OFF
          XI3MAX = 1.d0 + XI3OFF

          IF(CBBREV(CO,'EXTRAPOLATE',2,noco+1,NTCO,N3CO)) THEN
            EXTRAPOLATE=.TRUE.
          ELSE
            EXTRAPOLATE=.FALSE.
          ENDIF

          IF(CBBREV(CO,'CENTROID',3,noco+1,NTCO,N3CO)) THEN
            CENTROID=.TRUE.
            IF(CBBREV(CO,'CLUSTER',3,noco+1,NTCO,N3CO)) THEN
C LKC need to pass in array
C             CALL PARSIL(CO(N3CO+1),1,NTIL,NCLUSTER,ERROR,*9999)
              CALL PARSIL(CO(N3CO+1),1,NTIL,INT_TEMP(1),ERROR,*9999)
              NCLUSTER=INT_TEMP(1)
            ELSE
              NCLUSTER=72
            ENDIF
          ELSE
            CENTROID=.FALSE.
          ENDIF

!news     AAY 20 March 95 xi_3 keyword to specify what plane of xi_3
!         to project the data onto
!         also added "add" option to def da;c xi
          IF(CBBREV(CO,'XI_3',4,NOCO+1,NTCO,N3CO)) THEN
            SET_XI_3 = .TRUE.
C LKC 27-APR-98 Need to pass in array
C           CALL PARSRL(CO(N3CO+1),1,NTIL,XI_3,ERROR,*9999)
            CALL PARSRL(CO(N3CO+1),1,NTIL,REAL_TEMP,ERROR,*9999)
            XI_3=REAL_TEMP(1)

C news MPN 7Nov97: need to set logical here
          ELSE
! PJH 23Jul96 Don't set Xi3=0 in general case
c           XI_3 = 0.0d0
            SET_XI_3 = .FALSE.
C old
C           ! PJH 23Jul96 Don't set Xi3=0 in general case
C           c          ELSE
C           c            XI_3 = 0.0d0
C           c            SET_XI_3 = .FALSE.
          ENDIF
          TAG2D=.FALSE.
!         CONTOUR=.FALSE.
        ENDIF !beads/xi
!newe


C LKC 2-MAR-2005 Adding seed point options for closest face        
        IF(CBBREV(CO,'SEED_POINTS',4,NOCO+1,NTCO,N3CO)) THEN
          NSEED=IFROMC(CO(N3CO+1))
          IF(CLOSEST_FACE.EQV..FALSE.) THEN
            ERROR='>> Only implemented for Closest face calculations'
            GOTO 9999
          ENDIF
        ELSE
          NSEED=1
        ENDIF
        
C GDR 21-MAR-2005 Adding search start option for closest face        
        IF(CBBREV(CO,'SEARCH_START',4,NOCO+1,NTCO,N3CO)) THEN
          nsearch=IFROMC(CO(N3CO+1))
          IF(CLOSEST_FACE.EQV..FALSE.) THEN
            ERROR='>> Only implemented for Closest face calculations'
            GOTO 9999
          ENDIF
        ELSE
          nsearch=2
        ENDIF
                
        CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '    ERROR,*9999)

        IF(FILIO) THEN
C CPB 8/4/94 This needs to be generalised for NJ_LOC
          CALL STRING_TRIM(FILE,IBEG,IEND)

          IF((TYPE(1:5).EQ.'XIPOS').OR.(TYPE(1:5).EQ.'NODES')) THEN
            IF(CBBREV(CO,'TOTAL',1,noco+1,NTCO,N3CO)) THEN
C LKC 27-APR-98 need to pass in array
C             CALL PARSIL(CO(N3CO+1),1,NTIL,NITB,ERROR,*9999)
              CALL PARSIL(CO(N3CO+1),1,NTIL,INT_TEMP,ERROR,*9999)
              NITB=INT_TEMP(1)
            ELSE IF(NEELEM(1,NRLIST(1)).NE.0.AND.
     '          NBJ(1,NEELEM(1,NRLIST(1))).NE.0)THEN
              NITB=NIT(NBJ(1,NEELEM(1,NRLIST(1))))
            ELSE
              NITB=NJT-1
              NITB=NITB+NJ_LOC(NJL_FIBR,0,0)+NJ_LOC(NJL_FIEL,0,0)
            ENDIF
            NXIDEF=NITB
            CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.ipxi',STATUS,
     '        'SEQUEN','FORMATTED',132,ERROR,*9999)

C!!!
C!!! LKC 12-MAR-1999 Needs generalisation for multiple regions
C!!!
C           IF(TYPE(1:5).EQ.'NODES') THEN
C           CALL IOXI_NODE(IFILE,NEP,NITB,NPNODE,XIP,2,ERROR,*9999)
            IF(TYPE(1:5).EQ.'NODES') THEN
              CALL IOXI_NODE(IFILE,NEP,NITB,NPNODE,XIP,NRLIST2(1),
     '          ERROR,*9999)
            ELSE IF(TYPE(1:5).EQ.'XIPOS') THEN
              CALL IOXI(FD,IFILE,NDP,LD,LD_NP,NITB,NRE,NRLIST,
     '          XID,NODE_ASSOC,REFERENCE,ERROR,*9999)       
C GBS 31-Oct-2000 Moved from below (not relevant for nodes)
              IF(IOTYPE.EQ.2.AND.USE_DATA.EQ.1) THEN
                IF(KTYP11.EQ.1) THEN !element fitting
                  CALL DSTATS(LD,NDDL,NDLT,ERROR,*9999)
                ELSEIF(KTYP11.EQ.2) THEN !face fitting
                  CALL DSTATS_FACE(FD,LD,NDDL,NDLT,NFF,ERROR,*9999)
C PM 14Aug02 : Added for face fitting
                  IF(CBBREV(CO,'EXTERNAL',3,noco+1,NTCO,N3CO)) THEN
                    EXTERNAL=.TRUE.
                    CALL EXT_FACEFIT_LIST(NBJ,NELIST,NFF,NFLIST,
     '                NPF,EXTERNAL,ERROR,*9999)
                  ELSE
                    EXTERNAL=.FALSE.
                  ENDIF
                  IF(CBBREV(CO,'FACES',3,noco+1,NTCO,N3CO)) THEN
                    CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,
     '                ALL_REGIONS,ERROR,*9999)
                    CALL PARSE_FACES(NFFACE,NFLIST,noco-1,NRLIST,
     '                NTCO,CO,ERROR,*9999)
                  ENDIF
                ENDIF
              ENDIF

C LKC 28-FEB-2000 Moved from below (not relevant for nodes)
C             calc #points per frame (used in time history plots)
C MLB july 2001 - this bit of code causes errors
C             TREF=XID(NITB,1)
C             nd=1
C             DO WHILE(TREF.EQ.XID(NITB,nd).and.nd.LE.NDT)
C             nd=nd+1
C             ENDDO
C             NPPF=nd-1
              TREF=XID(NITB,1)
              nd=1
              DO WHILE(DABS(TREF-XID(NITB,nd)).LT.ZERO_TOL
     '          .and.nd.LT.NDT)
                nd=nd+1
              ENDDO
              IF(DABS(TREF-XID(NITB,nd)).LT.ZERO_TOL) THEN
                NPPF=nd
              ELSE
                NPPF=nd-1
              ENDIF
            ENDIF

            CALC_XI=.TRUE.
            CALL CLOSEF(IFILE,ERROR,*9999)

C LKC 28-FEB-2000 Moved into the data sections (not relevant for nodes)
C           C           calc #points per frame (used in time history
C           plots)
C           TREF=XID(NITB,1)
C           nd=1
C           DO WHILE(TREF.EQ.XID(NITB,nd).and.nd.LE.NDT)
C           nd=nd+1
C           ENDDO
C           NPPF=nd-1

          ENDIF !xipos

          
        ELSE IF(CALCU) THEN !Calculate data position info

          CALL ASSERT(CALL_ELEM,'>>Define elements first',ERROR,*9999)

C!!! LKC 16-NOV-1998 region hardcoded to 1, need to change
C         when data have regions
C         NXIDEF=NIT(NBJ(1,NEELEM(1,1)))
C         L=0
C         DO noelem=1,NEELEM(0,1) !to init LN (recalc.d by define fit)
C         ne=NEELEM(noelem,1)
C         L=L+1
C         LN(L)=ne
C         ENDDO
C         LN(0)=L

C LKC 30-MAY-2002, the L loop is unneccessary, use the noelem variable
C         also don't need the intermediate ne variable
C         L=0
C         DO noelem=1,NEELEM(0,NRLIST(1))
C         ne=NEELEM(noelem,NRLIST(1))
C         L=L+1 !to init LN (recalc.d by define fit)
C         LN(L)=ne
C         ENDDO
C         LN(0)=L

          NXIDEF=NIT(NBJ(1,NEELEM(1,NRLIST(1))))

C PM 14Aug02:
          IF(KTYP11.NE.2) THEN ! not face fitting
            CALL ASSERT(USE_DATA.EQ.1,'>Must set USE_DATA=1 in .ippara',
     &        ERROR,*9999)            
            DO noelem=1,NEELEM(0,NRLIST(1))
              LN(noelem)=NEELEM(noelem,NRLIST(1))
            ENDDO
            LN(0)=NEELEM(0,NRLIST(1))
          ENDIF

          IF(TYPE(1:5).EQ.'XIPOS') THEN

C!!! LKC 16-NOV-1998 region hardcoded to 1, need to change
C           when data have regions
C           
C           CALL ASSERT(NET(1).GT.0,
C           '      '>>no elements defined',ERROR,*9999)
            CALL ASSERT(NET(NRLIST(1)).GT.0,
     '        '>>no elements defined',ERROR,*9999)

            IF(ITYP10(1).EQ.1) THEN !rect. cart.
              NJ1=1
              NJ2=2
              IF(CENTROID) THEN
C**             find centroids and store in ZDD(1..3,nd)
                NTIMES=NDT/NCLUSTER
                DO nt=1,NTIMES
                  Z1(1)=0.0D0
                  Z1(2)=0.0D0
                  Z1(3)=0.0D0
                  DO nd=NCLUSTER*(nt-1)+1,NCLUSTER*nt
                    DO nj=1,NJT
                      Z1(nj)=Z1(nj)+ZD(nj,nd)
                    ENDDO
                  ENDDO
C                 centroid is
                  DO nj=1,NJT
                    Z1(nj)=Z1(nj)/NCLUSTER
                  ENDDO
                  DO nd=NCLUSTER*(nt-1)+1,NCLUSTER*nt
                    DO nj=1,NJT
                      ZDD(nj,nd)=Z1(nj)
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF !centroid

            ELSE IF(ITYP10(1).GE.2) THEN
              IF(NJT.EQ.2) THEN
                NJ1=1
                NJ2=2
              ELSE IF(NJT.EQ.3) THEN
                NJ1=2
                NJ2=3
              ENDIF
            ENDIF

C GMH HUGE comment block removed

            IF(CENTROID) THEN
              IF(NEW) THEN !assign approximate element locations
                DO nd=1,NDT !to initialize arrays
                  LD(nd)=0
                  SQ(nd)=0.0D0
                ENDDO
C               do for all elements
                DO nolist=1,NELIST(0)
                  ne=NELIST(nolist)
                  nb=NBJ(NJ1,ne)
                  NITB=NIT(NBJ(NJ1,ne))
                  IF(NITB.EQ.2) THEN          !2D elements
C                   do for all unassociated data points
                    DO nd=1,NDT
                      IF(LD(nd).EQ.0) THEN
C                       get vectors from centroid for all nodes of ne
                        DO nj=1,NJT
                          XD(nj)=ZD(nj,nd)-ZDD(nj,nd) !is proj vector
                          Z1(nj)=XP(1,nv,nj,NPNE(1,nb,ne))-ZDD(nj,nd)
                          Z2(nj)=XP(1,nv,nj,NPNE(2,nb,ne))-ZDD(nj,nd)
                          Z3(nj)=XP(1,nv,nj,NPNE(3,nb,ne))-ZDD(nj,nd)
                          Z4(nj)=XP(1,nv,nj,NPNE(4,nb,ne))-ZDD(nj,nd)
                        ENDDO
C                       convert to spherical polar coords
                        CALL ZX(3,XD,XD)
                        CALL ZX(3,Z1,Z1)
                        CALL ZX(3,Z2,Z2)
                        CALL ZX(3,Z3,Z3)
                        CALL ZX(3,Z4,Z4)
C                       do for both angles
                        DO IT=2,3
                          TMIN(IT-1)=DMIN1(Z1(IT),Z2(IT),Z3(IT),Z4(IT))
                          TMAX(IT-1)=DMAX1(Z1(IT),Z2(IT),Z3(IT),Z4(IT))
C                         assume element spans smallest
C                         angle between nodes
                          IF(TMAX(IT-1)-TMIN(IT-1).GT.
     '                      TMIN(IT-1)+2.0D0*PI-TMAX(IT-1)) THEN
                            IF(Z1(IT).LT.0.0D0) Z1(IT)=Z1(IT)+2.0D0*PI
                            IF(Z2(IT).LT.0.0D0) Z2(IT)=Z2(IT)+2.0D0*PI
                            IF(Z3(IT).LT.0.0D0) Z3(IT)=Z3(IT)+2.0D0*PI
                            IF(Z4(IT).LT.0.0D0) Z4(IT)=Z4(IT)+2.0D0*PI
                            IF(XD(IT).LT.0.0D0) XD(IT)=XD(IT)+2.0D0*PI
                            TMIN(IT-1)=DMIN1(Z1(IT),Z2(IT),Z3(IT),
     '                        Z4(IT))
                            TMAX(IT-1)=DMAX1(Z1(IT),Z2(IT),Z3(IT),
     '                        Z4(IT))
                          ENDIF
                        ENDDO
C                       IF(DOP)WRITE(*,*)' Data point',nd,'element',ne
C                       IF(DOP)WRITE(*,*)tmin(1),tmax(1),tmin(2),
C                       '                    tmax(2),xd(2),xd(3)
                        IF(XD(2).GE.TMIN(1).AND.XD(2).LE.TMAX(1).AND.
     '                    XD(3).GE.TMIN(2).AND.XD(3).LE.TMAX(2)) THEN
C                         associate data point with element
                          LD(nd)=ne
                          XID(1,nd)=0.5D0
                          XID(2,nd)=0.5D0
                          IF(DOP) THEN
                            WRITE(OP_STRING,*)' Data point',nd,
     '                        ' is ','in element',ne
                            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                          ENDIF
                        ENDIF
                      ENDIF
                    ENDDO !for all data points
                  ENDIF !2D element
                ENDDO !for all elements
              ENDIF !new
C             find neighbouring elements
              CALL NENXI(NBJ,NEELEM,NENP,NNB,NPNE,NXI,ERROR,*9999)

C GMH HUGE comment block removed

C             now loop over data points
              DO nd=1,NDT
C MLB 21/3/97 initialising SQND
                SQND=1.0d0
                IF(LD(nd).EQ.0) THEN
                  LD(nd)=1
                  XID(1,nd)=0.5D0
                  XID(2,nd)=0.5D0
                ENDIF
                IF(INLIST(LD(nd),NELIST(1),NELIST(0),n1list)) THEN
                  IF(DOP) THEN
                    WRITE(OP_STRING,
     '                '(//,'' *****>Data point '',I5)') nd
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                  !initialise ne and xi to starting position
                  ne=LD(nd)
                  XI(1)=XID(1,nd)
                  XI(2)=XID(2,nd)
                  XI(3)=1.0D0
                  FINISHED=.FALSE.
                  DO WHILE(.NOT.FINISHED)
                    IF(DOP) THEN
                      WRITE(OP_STRING,
     '                  '(/'' Element '',I4,'' Xi ='',3F10.4)')
     '                  ne,XI(1),XI(2),XI(3)
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                    CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),
     '                NPF(1,1),NPNE(1,1,ne),
     '                NRE(ne),NVJE(1,1,1,ne),
     '                SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
                    NITB=NIT(NBJ(NJ1,ne))
                    IF(NITB.EQ.1) THEN !1D elements
                    ELSE IF(NITB.EQ.2) THEN !2D elements
                      IF(ITYP10(1).EQ.1) THEN !rect cart coords
                        DO nj=1,NJT
                          Z1(nj)=ZD(nj,nd)
                          Z2(nj)=ZDD(nj,nd)
                        ENDDO
                        ERROR='>> PROJ21 Moved to FE_ARCHIVE'
                        GOTO 9999
C                       CALL PROJ21(IBT,IDO,INP,IT,ITMAX,NBJ(1,ne),ne,
C                       '                    SQND,XE,XI,Z1,Z2,ERROR,
C                       *9999)
                      ENDIF
                      IF(DABS(SQND).LT.10D-1) THEN !converged
                        XID(1,nd)=XI(1)
                        XID(2,nd)=XI(2)
                        LD(nd)=ne
                        SQ(nd)=SQND
                        FINISHED=.TRUE.
                        IF(DOP) THEN
                          WRITE(OP_STRING,
     '                      '('' Convergence reached''/'
     '                      //''' nd='',I4,'' LD='',I4,'' XID='','
     '                      //'2E12.6,'//''' SQ='',E12.6,'' ZD='','
     '                      //'3(E12.6,1X))')
     '                      nd,LD(nd),(XID(ni,nd),ni=1,2),SQ(nd),
     '                      (ZD(nj,nd),nj=1,NJT)
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                        ENDIF
                      ELSE IF(XI(1).GT.0.0D0.AND.XI(1).LT.1.0D0.AND
     '                    .XI(2).GT.0.0D0.AND.XI(2).LT.1.0D0) THEN
                        !still inside element and intersection not found
                        WRITE(OP_STRING,'('' WARNING:'
     '                    //' Convergence not reached in PROJ2X''/'
     '                    //''' nd='',I4,'' LD='',I4,'' XID='',2E12.6,'
     '                    //''' SQ='',E12.6,'' ZD='',3(E12.6,1X))')
     '                    nd,LD(nd),(XID(ni,nd),ni=1,2),SQ(nd),
     '                    (ZD(nj,nd),nj=1,NJT)
                        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                        !save current, set ld(nd) to zero as flag
                        XID(1,nd)=XI(1)
                        XID(2,nd)=XI(2)
                        SQ(nd)=SQND
                        LD(nd)=0
                        FINISHED=.TRUE.
                      ENDIF
                      IF(.NOT.FINISHED) THEN !at a boundary
                        !shift to neighbour element & save solution so far
                        SQ(nd)=SQND
                        XID(1,nd)=XI(1)
                        XID(2,nd)=XI(2)
                        LD(nd)=ne
                        !find the element moved to
                        IF(XI(1).LE.0.0D0) THEN
                          XID(1,nd)=0.0D0
                          XI(1)=1.0D0
                          ne=NXI(-1,1,ne)
                        ELSE IF(XI(1).GE.1.0D0) THEN
                          XID(1,nd)=1.0D0
                          XI(1)=0.0D0
                          ne=NXI(1,1,ne)
                        ENDIF
                        IF(XI(2).LE.0.0D0) THEN
                          XID(2,nd)=0.0D0
                          XI(2)=1.0D0
                          ne=NXI(-2,1,ne)
                        ELSE IF(XI(2).GE.1.0D0) THEN
                          XID(2,nd)=1.0D0
                          XI(2)=0.0D0
                          ne=NXI(2,1,ne)
                        ENDIF
                        FINISHED=(ne.EQ.0)
                        IF(DOP) THEN
                          WRITE(OP_STRING,'(/'' Move to Element '','
     '                      //'I4,'' Xi ='',2F10.4)') ne,XI(1),XI(2)
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                        ENDIF
                      ENDIF !still in element
                    ENDIF !2d element
                  ENDDO !until finished
                ENDIF !in element list
              ENDDO  !for all data points

C GMH HUGE comment block removed

!news AAY 9 Dec 94
            ELSE IF(TAG2D) THEN
              !for all stripe data points do

C LC 25/2/97 archived section :  PJH 8Dec95

            ELSE

C LKC 30-MAY-2002 Doesn't seem necessary to initialise the arrays
C             here as LD and SQ are initialised in the DEXI_* routines
C             C CS 5/3/1998 Split into problem specific routines
C             IF(NEW) THEN
C             DO ND=ND0,ND1 !to initialize arrays
C             LD(ND)=0
C             SQ(ND)=0.d0
C             ENDDO
C             ENDIF

              NITB=NIT(NBJ(1,NELIST(1)))
              IF(ORTHOG) THEN
                CALL DEXI_ORTHOG(IBT,IDO,INP,LD,NBJ,
     '            NBH,nd,ND0,ND1,ne,NELIST,NHE,ni,NITB,nj1,
     '            NKHE,NKJE,NPF,NPNE,nolist,nr,NRE,NVHE,NVJE,NW(1,1,nx),
     '            nx,START_ELEMENT,CURVCORRECT,SE,SQ,SQND,XA,XE,XI,
     '            XI_1,XI_2,XI_3,XID,XP,ZA,ZD,ZP,DEFORM,
     '            NEW,SET_XI_1,SET_XI_2,SET_XI_3,SPECIFY,
     '            ERROR,*9999)
              ELSE IF(CLOSEST) THEN
                CALL DEXI_CLOSEST(IBT,IDO,INP,LD,NBJ,
     '            NBH,ND0,ND1,NELIST,NHE,nj1,
     '            NKHE,NKJE,NPF,NPNE,NRE,NVHE,NVJE,NW(1,1,nx),nx,NXI,
     '            START_ELEMENT,CURVCORRECT,SE,SQ,SQMAX,XA,XE,XI,
     '            XI_1,XI_2,XI_3,XID,XP,ZA,ZD,ZP,DEFORM,
     '            NEW,SET_XI_1,SET_XI_2,SET_XI_3,SPECIFY,
     '            ERROR,*9999)
              ELSE IF(SUBDIV) THEN
                CALL DEXI_SUBDIV(
     '            DEFORM,IBT,IDO,INP,LD,NBH,NBJ,ND0,ND1,
     '            NELIST,NHE,NKHE,NKJE,
     '            NOXIPT,NPF,NPNE,NRE,NVHE,NVJE,NW,nx,CURVCORRECT,SE,SQ,
     '            SQMAX,XA,XE,XID,XP,ZA,ZD,ZP,
     '            ERROR,*9999)
C Glenn Ramsey 16/2/05 added 'frustum' face projection     
              ELSE IF(CLOSEST_FACE.OR.FRUSTUM) THEN ! face proj subroutine
                CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,
     '            ERROR,*9999)
                IF(CBBREV(CO,'EXTERNAL',3,noco+1,NTCO,N3CO)) THEN
                  EXTERNAL=.TRUE.
                  CALL EXT_FACEFIT_LIST(NBJ,NELIST,NFF,NFLIST,
     '              NPF,EXTERNAL,ERROR,*9999)
                ELSE
                  EXTERNAL=.FALSE.
                ENDIF
                IF(CBBREV(CO,'FACES',3,noco+1,NTCO,N3CO)) THEN
                  CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,
     '              ALL_REGIONS,ERROR,*9999)
                  CALL PARSE_FACES(NFFACE,NFLIST,noco-1,NRLIST,NTCO,CO,
     '              ERROR,*9999)
C PM 14Aug02:
C                 NELIST(0)=0
C                 DO nf=1,NFLIST(0)
C                 NELIST(0)=NELIST(0)+1
C                 NELIST(NELIST(0))=NPF(6,NFLIST(nf))
C                 ENDDO !nf
C                 ELSE
C                 NFLIST(0)=0
                ENDIF
C GDR 20Mar2005 Adding orthogonal option for closest face        
                CFORTHO=.FALSE.
                IF(CBBREV(CO,'ORTHOGONAL',2,NOCO+1,NTCO,N3CO)) THEN
                  CFORTHO=.TRUE.
                  IF(CLOSEST_FACE.EQV..FALSE.) THEN
                    ERROR=
     &                '>> Only implemented for closest face calcs'
                    GOTO 9999
                  ENDIF
                ENDIF
C GDR 20Mar2005 Adding cross_boundaries option for closest face        
                CFCROSS=.FALSE.
                IF(CBBREV(CO,'CROSS_BOUNDARIES',
     &              4,NOCO+1,NTCO,N3CO)) THEN
                  CFCROSS=.TRUE.
                  IF(CLOSEST_FACE.EQV..FALSE.) THEN
                    ERROR=
     &                '>> Only implemented for closest face calcs'
                    GOTO 9999
                  ENDIF
                ENDIF

                IF(CLOSEST_FACE) THEN 
                  CALL DEXI_CLOSEST_FACE(FD,IBT,IDO,INP,LD,NBH,NBHF,NBJ,
     &              NBJF,ND0,ND1,NFF,NFLIST,NHE,NKEF,
     &              NKHE,NKJE,NNF,NPF,NPNE,NPNF,NRE,NSEED,nsearch,
     &              NVHE,NVJE,NVJF,nx,NXI,SE,SF,SQ,XA,XE,XI,
     &              XI_1,XI_2,XI_3,XID,XP,ZD,ZE,ZP,DEFORM,NEW,
     &              CFORTHO,SET_XI_1,SET_XI_2,SET_XI_3,SPECIFY,CFCROSS,
     &              ERROR,*9999)
     
                ELSEIF(FRUSTUM) THEN
                  CALL DEXI_IN_FRUSTUM(FD,IBT,IDO,INP,LD,NBH,NBJ,
     &              NBJF,ND0,ND1,NFF,NFLIST,NEELEM,NHE,NKEF,
     &              NKHE,NKJE,NNF,NPF,NPNE,NPNF,NRE,NVHE,NVJE,NVJF,NW,
     &              nx,NXI,CURVCORRECT,SE,SF,SQ,XA,XE,XI,
     &              XI_1,XI_2,XI_3,XID,XP,ZA,ZD,ZP,DEFORM,NEW,
     &              SET_XI_1,SET_XI_2,SET_XI_3,SPECIFY,ERROR,*9999)
                ENDIF

              ELSE
                IF(NITB.EQ.1) THEN !1D elements
                  CALL DEXI_1D(IBT,IDO,INP,LD,nb,NBJ,
     '              NBH,nd,ND0,ND1,ne,NELIST,NHE,ni,NITB,nj1,nj2,
     '              NKHE,NKJE,NKTB,NPF,NPNE,nolist,nr,NRE,ns1,ns2,NVHE,
     '              NVJE,NW(1,1,nx),nx,START_ELEMENT,A1,A2,
     '              CURVCORRECT,D1,D2,DIST,
     '              SE,SLOPE,SQ,SQND,X1,X2,XD,XI,XI_1,XI_2,XI_3,
     '              XA,XE,XID,XP,ZA,ZD,ZP,DEFORM,FOUND,NEW,ORTHOG,
     '              SET_XI_1,SET_XI_2,SET_XI_3,SPECIFY,ERROR,*9999)
                ELSE IF(NITB.GE.2) THEN !2D/3D elements
C GMH 30/10/95 Commented regions can be found in cmiss_archive.
C Not sure why this division is here for 2D elements, but I will leave
C it.  Should this just default to linear?
                  IF(CUBIC) THEN !use 1D nonlinear calculation
C GMH HUGE comment block removed
                  ELSE IF(QUADRATIC) THEN !use 2D nonlinear calculation
C GMH HUGE comment block removed
                  ELSE IF((.NOT.ORTHOG).AND.INTERNAL) THEN
C calculates the xi position within an element
C CS 7/3/98 Split into 3 passes (one linear approx and
C                   2 nonlinear) to increase speed

                    CALL DEXI_LINEAR(IBT,IDO,INP,LD,NBJ,
     '                NBH,nb,nd,ND0,ND1,ne,NELIST,NHE,ni,NITB,nj1,nj2,
     '                NKHE,NKJE,NKTB,NPF,NPNE,nolist,nr,NRE,ns1,ns2,ns3,
     '                ns4,NVHE,NVJE,NW(1,1,nx),nx,
     '                START_ELEMENT,A1,A2,ALFA,B1,B2,
     '                BETA,C1,C2,CURVCORRECT,D1,D2,DELTA,DENOM1,DENOM2,
     '                DIFF,GAMA,SE,SQ,SQND,THETAMIN,THETAMAX,X0,X1,X2,
     '                XA,XD,XI,XE,XI_1,XI_2,XI_3,XID,XP,ZA,ZD,ZP,
     '                DEFORM,EXCLUDE,EXTRAPOLATE,FOUND,NEW,ORTHOG,
     '                SET_XI_1,SET_XI_2,SET_XI_3,SPECIFY,ERROR,*9999)

                    IF(.NOT.FAST)THEN
                      PASS2=.FALSE.
                      CALL DEXI_NONLIN(IBT,IDO,INP,LD,LDTEMP,NBJ,NBH,
     '                  nd,ND0,ND1,ne,NELIST,NHE,ni,NITB,nj1,NKHE,NKJE,
     '                  NPF,NPNE,nolist,nr,NRE,NVHE,NVJE,NW(1,1,nx),nx,
     '                  CURVCORRECT,
     '                  SE,XA,XE,XI,XID,%VAL(0),XP,%VAL(0),
     '                  ZA,ZD,ZP,DEFORM,FOUND,GRPGRID,PASS2,
     '                  SPECIFY,ERROR,*9999)

                      PASS2=.TRUE.
                      CALL DEXI_NONLIN(IBT,IDO,INP,LD,LDTEMP,NBJ,NBH,
     '                  nd,ND0,ND1,ne,NELIST,NHE,ni,NITB,nj1,NKHE,NKJE,
     '                  NPF,NPNE,nolist,nr,NRE,NVHE,NVJE,NW(1,1,nx),nx,
     '                  CURVCORRECT,
     '                  SE,XA,XE,XI,XID,%VAL(0),XP,%VAL(0),
     '                  ZA,ZD,ZP,DEFORM,FOUND,GRPGRID,PASS2,
     '                  SPECIFY,ERROR,*9999)
                    ELSE !fast approach
                      CALL DEXI_NONLIN_CONTAIN(IBT,IDO,INP,LD,LDTEMP,
     '                  NBJ,NBH,nd,ND0,ND1,ne,NELIST,NHE,ni,NITB,nj1,
     '                  NKHE,NKJE,NPF,NPNE,nolist,nr,NRE,NVHE,NVJE,
     '                  NW(1,1,nx),nx,NXI,CURVCORRECT,
     '                  SE,XA,XE,XI,XID,%VAL(0),XP,%VAL(0),
     '                  ZA,ZD,ZP,DEFORM,FOUND,GRPGRID,PASS2,
     '                  SPECIFY,ERROR,*9999)
                    ENDIF

                  ELSE IF(.NOT.ORTHOG) THEN !use linear Xi calculation
                    CALL DEXI_LINEAR(IBT,IDO,INP,LD,NBJ,
     '                NBH,nb,nd,ND0,ND1,ne,NELIST,NHE,ni,NITB,nj1,nj2,
     '                NKHE,NKJE,NKTB,NPF,NPNE,nolist,nr,NRE,ns1,ns2,ns3,
     '                ns4,NVHE,NVJE,NW(1,1,nx),nx,
     '                START_ELEMENT,A1,A2,ALFA,B1,B2,
     '                BETA,C1,C2,CURVCORRECT,D1,D2,DELTA,DENOM1,DENOM2,
     '                DIFF,GAMA,SE,SQ,SQND,THETAMIN,THETAMAX,X0,X1,X2,
     '                XA,XD,XI,XE,XI_1,XI_2,XI_3,XID,XP,ZA,ZD,ZP,
     '                DEFORM,EXCLUDE,EXTRAPOLATE,FOUND,NEW,ORTHOG,
     '                SET_XI_1,SET_XI_2,SET_XI_3,SPECIFY,ERROR,*9999)
                  ENDIF !linear/orthog/cubic/internal
                ENDIF !nitb=1,2,3
              ENDIF !orthog/non-orthog

              IF(EXTRAPOLATE) THEN
                DO nd=1,NDT
                  IF(LD(nd).LT.0) THEN !data point currently not in any element
                    LD(nd)=-LD(nd)
                    XID(3,nd)=1.0D0
                  ENDIF
                ENDDO
              ENDIF
            ENDIF !orthog/centroid etc
C PM 14Aug02:
            IF(KTYP11.EQ.1) THEN
              CALL DSTATS(LD,NDDL,NDLT,ERROR,*9999)

              DO nolist=1,NELIST(0)
                ne=NELIST(nolist)
                CALL ASSERT(NDLT(ne).LE.NDEM,
     '            '>>NDEM too small',ERROR,*9999)
              ENDDO
            ELSEIF(KTYP11.EQ.2) THEN
              CALL DSTATS_FACE(FD,LD,NDDL,NDLT,NFF,ERROR,*9999)
            ENDIF
            CALC_XI=.TRUE.

          ELSE IF(TYPE(1:5).EQ.'NODES') THEN
C calculate the xi position of nodes in nr2 contained by host elements
C in nr1

C!!! Start of code for xi calculations (finished command parsing)
 !KSB 2007 Adding initialisation of array NEP
            DO nonode=1,NPNODE(0,nr2)
              np=NPNODE(nonode,nr2)
              NEP(np)=0
            ENDDO

C!!! LKC 16-NOV-1998
C!!! Need to loop over NRLIST when data has regions
C!!! Should change nr1, nr2 to use NRLIST ?
C           NITB=NIT(NBJ(1,NEELEM(1,1)))
            NITB=NIT(NBJ(1,NEELEM(1,nr1)))
            CALL ASSERT(NITB.eq.3,
     '        '>>only implemented in 3D',ERROR,*9999)
            CALL ASSERT(ITYP10(nr1).eq.1,
     '        '>>only implemented in cartesian',ERROR,*9999)

C LKC 12-MAR-1999 Adding loop over nrlist2 for region generalisation
            IF(.NOT.CLOSEST) THEN
              DO nonr=1,NRLIST2(0)
                nr2=NRLIST2(nonr)
                noelem=0
                DO nplist=1,NPNODE(0,nr2)
                  np=NPNODE(nplist,nr2)
                  IF(ORDER_CLOSEST)THEN
C MHT 18-11-02 New: routine to order search sequence of elements.
                    DO nj=1,NJT
                      XD(nj)=XP(1,1,nj,np)
                    ENDDO
                    CALL DEXI_ORDER_ELEMENTS(NEELEM(0,nr1),NELIST2,
     '                NENP(1,0,nr1),NPNODE(0,nr1),XD,XP,ERROR,*9999)
C MHT end 18-11-02
                  ENDIF
                  nolist=1
                  LDTEMP=0
                  DO WHILE ((nolist.le.NEELEM(0,nr1)).AND.(LDTEMP.EQ.0))
c old             ne=NEELEM(nolist,nr1)
                    IF(ORDER_CLOSEST)THEN
                      ne=NELIST2(nolist)
                    ELSE
                      ne=NEELEM(nolist,nr1)
                    ENDIF
                    CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),
     '                NPF(1,1),NPNE(1,1,ne),
     '                NRE(ne),NVJE(1,1,1,ne),
     '                SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
                    DO ni=1,NITB !initialising XI
                      XI(ni)=0.5d0
                      NODE_POS(ni)=XP(1,1,ni,np)
                    ENDDO
                    CALL DEXI_POINT(IBT,IDO,INP,LDTEMP,NBJ,
     '                ne,NITB,NRE(ne),FACTOR_TOL,XE,XI,XIP(1,np),
     '                NODE_POS,USE_LOOSE_TOL,ERROR,*9999)
                    IF(LDTEMP.GT.0) NEP(np)=LDTEMP
                    nolist=nolist+1
                  ENDDO
                  IF(NEP(np).EQ.0) noelem=noelem+1 !counter for lost nodes
                ENDDO !nplist

C LKC 17-APR-1999 added warning
                IF(noelem.NE.0) THEN
                  WRITE(OP_STRING,
     '              '('' >> WARNING: No elements found for'',I5,
     '              '' nodes in region'',I2)') noelem,nr2
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDIF
              ENDDO !nonr
            ELSE
              CALL DEXI_CLOSEST_NODE(IBT,IDO,INP,NBJ,
     '          NBH,NELIST,NEP,NEW,NHE,nj1,NKHE,NKJE,NPF,NPNE,
     '          NPNODE,NRE,NRLIST2,NVHE,NVJE,NW(1,1,nx),nx,
     '          NXI,START_ELEMENT,CURVCORRECT,SE,XA,XE,XI,
     '          XI_1,XI_2,XI_3,XIP,XP,ZA,ZP,
     '          ERROR,*9999)
            ENDIF !.NOT.CLOSEST
            IF(NUM_FIELD)THEN
              DO nonode=1,NPLIST2(0)
               np=NPLIST2(nonode)
               CALL ASSERT(NEP(np).NE.0,
     '      '>>Terminal node outside tissue block',ERROR,*9999)
               ne=NEP(np)
               XAB(nj_num,ne)=XAB(nj_num,ne)+1.d0
               XP(1,1,nj_elem,np)=ne
              ENDDO
            ENDIF
          ELSE IF(TYPE(1:5).EQ.'BEADS') THEN
            nb=NBJ(1,1)
            NITB=NIT(nb)
            NKTB=NKT(0,nb)
            IF(NET(1).GT.0.AND.NITB.EQ.3.AND.NJT.EQ.3.AND.NDT.GT.6) THEN
C ***         Calculate xi-coordinates of projections of bead positions
C ***         on to xi_3 lines of 3D elements using normalized
C ***         lengths in each column. It is assumed that data points
C ***         1,2,3 are epicardial beads numbered clockwise, that data
C ***         points 4,5,6 are bottom beads of cols 1,2,3, and that
C ***         elements numbers increase from epicardium to endocardium.
              DO icol=1,3
                DO ni=1,NITB
                  CVEC(ni,icol)=ZD(ni,icol+3)-ZD(ni,icol)
                ENDDO
                CMAG(icol)=0.0D0
                DO ni=1,NITB
                  CMAG(icol)=CMAG(icol)+CVEC(ni,icol)**2
                ENDDO
                CMAG(icol)=DSQRT(CMAG(icol))
              ENDDO
              XID(1,1)=0.0D0
              XID(2,1)=0.0D0
              XID(3,1)=0.0D0
              XID(1,2)=1.0D0
              XID(2,2)=0.0D0
              XID(3,2)=0.0D0
              XID(1,3)=0.5D0
              XID(2,3)=1.0D0
              XID(3,3)=0.0D0
              XID(1,4)=0.0D0
              XID(2,4)=0.0D0
              XID(3,4)=NET(1)
              XID(1,5)=1.0D0
              XID(2,5)=0.0D0
              XID(3,5)=NET(1)
              XID(1,6)=0.5D0
              XID(2,6)=1.0D0
              XID(3,6)=NET(1)
              DO nd=7,NDT
                DO icol=1,3
                  BMAG(icol)=0.0D0
                  BDOTC=0.0D0
                  DO ni=1,3
                    BMAG(icol)=BMAG(icol)+(ZD(ni,nd)-ZD(ni,icol))**2
                    BDOTC=BDOTC+(ZD(ni,nd)-ZD(ni,icol))*CVEC(ni,icol)
                  ENDDO
                  BMAG(icol)=DSQRT(BMAG(icol))
                  AB(icol)=DABS(BDOTC/(BMAG(icol)*CMAG(icol)))
                ENDDO
                IF(AB(1).GT.AB(2).AND.AB(1).GT.AB(3)) THEN
                  XID(1,nd)=0.0D0
                  XID(2,nd)=0.0D0
                  XID(3,nd)=BMAG(1)/CMAG(1)*NET(1)
                ELSE IF(AB(2).GT.AB(1).AND.AB(2).GT.AB(3)) THEN
                  XID(1,nd)=1.0D0
                  XID(2,nd)=0.0D0
                  XID(3,nd)=BMAG(2)/CMAG(2)*NET(1)
                ELSE IF(AB(3).GT.AB(1).AND.AB(3).GT.AB(2)) THEN
                  XID(1,nd)=0.5D0
                  XID(2,nd)=1.0D0
                  XID(3,nd)=BMAG(3)/CMAG(3)*NET(1)
                ENDIF
              ENDDO
              DO nd=1,NDT
                i=INT(XID(3,nd))
                IF(i.GE.0.AND.i.LE.(NET(1)-1)) THEN
                  LD(nd)=i+1
                  XID(3,nd)=XID(3,nd)-i
                ELSE IF(i.LT.0) THEN
                  LD(nd)=1
                  XID(3,nd)=0.0D0
                ELSE IF(i.GT.(NET(1)-1)) THEN
                  LD(nd)=i
                  XID(3,nd)=1.0D0
                ENDIF
                CALL DSTATS(LD,NDDL,NDLT,ERROR,*9999)
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' Bead point '',I3,'' lies in '
     '              //'element '',I3,'' at Xi coords: '',3E11.3)')
     '              nd,LD(nd),(XID(ni,nd),ni=1,NITB)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
              ENDDO
            ELSE IF(NET(1).EQ.0) THEN
              OP_STRING(1)=' >>no elements defined'
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ELSE IF(njt.NE.3.or.nitb.NE.3) THEN
              OP_STRING(1)=' >>Bead option requires 3D elem.s'
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ELSE IF(NDT.LT.6) THEN
              OP_STRING(1)=' >>At least six beads are needed'
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
C JWF 15/01/02 - Defines contact points 
          ELSE IF(TYPE(1:7).EQ.'CONTACT') THEN

            IF(CBBREV(CO,'FACES',2,noco+1,NTCO,N3CO)) THEN
              CALL PARSE_FACES(NFFACE,NFLIST,noco-1,NRLIST,NTCO,CO,
     '          ERROR,*9999)
            ENDIF

            i=0
            IF(CBBREV(CO,'XI_END',2,noco+1,NTCO,N3CO)) THEN
              DO nonr=1,NRLIST(0)
                nr=NRLIST(nonr)
                DO noface=1,NFFACE(0,nr) 
                  nf=NFFACE(noface,nr)
                  IF (NPF(5,nf).EQ.1) THEN !external
                    IF ((NPF(1,nf).EQ.2).AND.(NPF(3,nf).EQ.3)) THEN
                      i=i+1
                      NFLIST(i)=nf
                      NFLIST(0)=i                 
                    ENDIF
                  ENDIF
                ENDDO !noface
              ENDDO !nonr
            ENDIF

            IF(CBBREV(CO,'POINTS',2,noco+1,NTCO,N3CO)) THEN
              CALL PARSIL(CO(N3CO+1),1,NTLIST,NILIST,
     '          ERROR,*9999)
C*** 22/02/08 JHC Assume number of contact points in xi_1 is the same as in xi_2
              NDFXI(1)=NILIST(1) ! data point # per face
              NDFXI(2)=NILIST(1) ! data point # per face

C*** 22/02/08 JHC allow different number of contact points at different face
C                 cont_wt is not used anymore
C              CONT_WT=NDFXI(1)**2
            ELSE
              NDFXI(1)=4 ! default is 4
              NDFXI(2)=4 ! default is 4
C              CONT_WT=NDFXI(1)**2
            ENDIF
            CONT_WEIGHT=NDFXI(1)*NDFXI(1)

C*** 22/02/08 JHC Number of contact points in xi_2 different from in xi_1
            IF(CBBREV(CO,'BY',2,noco+1,NTCO,N3CO)) THEN
              CALL PARSIL(CO(N3CO+1),1,NTLIST,NILIST,
     '          ERROR,*9999)
              NDFXI(2)=NILIST(1) ! data point # per face
C              CONT_WT=NDFXI(1)*NDFXI(2)
              CONT_WEIGHT=NDFXI(1)*NDFXI(2)
            ENDIF

C*** 22/02/08 JHC option to put contact points on element boundary as well as internally
            IF(CBBREV(CO,'BOUNDARY',4,noco+1,NTCO,N3CO)) THEN
              BOUNDARY_CONT=.TRUE.
            ELSE
              BOUNDARY_CONT=.FALSE.
            ENDIF

C new JHC/MPN 19Jul2004 added check for size of Z_CONT_LIST
            WRITE(CHAR1,'(I5)') NFLIST(0)*NDFXI(1)*NDFXI(2)
            CALL STRING_TRIM(CHAR1,IBEG,IEND)
            CALL ASSERT((NFLIST(0)*NDFXI(1)*NDFXI(2)).LE.NDM,
     &        '>>Increase NDM to at least '//CHAR1(IBEG:IEND),ERROR,
     &        *9999)
 
C new JHC/MPN 29-OCT-2004 Select subset of contact points within faces
            direction=0 !initialize variable
            IF(CBBREV(CO,'XMAX',4,noco+1,NTCO,N3CO)) THEN
              direction=1
              CALL PARSRL(CO(N3CO+1),1,NTIL,XD,ERROR,*9999)
              LIMIT=XD(1)
            ELSE IF(CBBREV(CO,'YMAX',4,noco+1,NTCO,N3CO)) THEN
              direction=2
              CALL PARSRL(CO(N3CO+1),1,NTIL,XD,ERROR,*9999)
              LIMIT=XD(1)
            ELSE IF(CBBREV(CO,'ZMAX',4,noco+1,NTCO,N3CO)) THEN
              direction=3
              CALL PARSRL(CO(N3CO+1),1,NTIL,XD,ERROR,*9999)
              LIMIT=XD(1)
            ENDIF
c end new

C JWF 09/02/05: Option added so that can define contact points to be different
C types (i.e. frictionless and friction) within the same region. If 'ADD' is specified then 
C the contact points global data pt # will continue from previous function call rather than be re-initialised.
            IF(CBBREV(CO,'ADD',2,noco+1,NTCO,N3CO)) THEN
C              CALL PARSIL(CO(N3CO+1),1,NTLIST,NILIST,
C     '          ERROR,*9999)
              NEW_CONT_PTS=.TRUE. 
              nd_increment=NDT ! set to previous total value         
            ELSE
              NEW_CONT_PTS=.FALSE. 
              nd_increment=0  ! initialise global data pt #
            ENDIF

            DO nf_global=1,NFLIST(0) !loop over faces

              nf=NFLIST(nf_global) !global #
              ne=NPF(6,nf) !element #
              nef=NPF(8,nf)  !local face #
              nr=NRE(ne)  !region #
              nbface=NBJF(1,nf)
     
C             Get face information
              CALL CALC_FACE_INFORMATION_IND(NBJ(1,ne),NBJF(1,nf),
     '          nef,NKJE(1,1,1,ne),NKEF,NKJF,NNF,NPNE(1,1,ne),NPNF,
     '          nr,NVJE(1,1,1,ne),NVJF,SE(1,1,ne),SF,ERROR,*9999)
  
C new JHC/MPN 29-OCT-2004 Select subset of contact points within faces
C*** 22/02/08 JHC XPXE is needed in the case when new contact points are added 
C                 which may concide with existing contact points along element boundary
              IF(direction.NE.0.OR.BOUNDARY_CONT) THEN
                CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),
     &            NPF(1,1),NPNE(1,1,ne),NRE(ne),
     &            NVJE(1,1,1,ne),SE(1,1,ne),
     &            XA(1,1,ne),XE,XP,ERROR,*9999)
              ENDIF
c end new

C             set initial xi_1,xi_2 values for ND1,ND2 loops below
              IF(.NOT.BOUNDARY_CONT) THEN ! contact points on internal element 
                XI_1=(1.0d0/DBLE(NDFXI(1)))*0.5d0
                XI_2=(1.0d0/DBLE(NDFXI(2)))*0.5d0
              ELSE ! contact points on element boundary as well as on internal element
                XI_1=0.0d0
                XI_2=0.0d0
              ENDIF

C             Calculate third XI direction ie 0.0 or 1.0
              DO nn=1,NNT(nbface)
                A_LIST(nn)=NPNF(nn,nbface)
              ENDDO
              nb=NBJ(1,ne)
              IF(INLIST(NPNE(1,nb,ne),A_LIST,NNT(nbface),
     '        N1LIST)) THEN
                XI_3=0.0d0
              ELSE
                XI_3=1.0d0
              ENDIF
      
              DO nd1=1,NDFXI(1)
                DO nd2=1,NDFXI(2)
                  
                  nd_increment=nd_increment+1 ! global data pt #

c                 Set status of points 1-frictionless, 2-tied contact
c                 and 3-including friction.
                  Z_CONT_LIST(nd_increment,1,4)=CONTACT_STATUS
                  Z_CONT_LIST(nd_increment,2,4)=CONTACT_STATUS

C*** 22/02/08 JHC Add a check statement to see if any new contact points have been added
C                 It is important for frictional contact as their frictional forces need 
C                 to be intialised to be zero
C WARNING: it assumes that these values are zero or at least less than 1 when uninitialised
                  IF(Z_CONT_LIST(nd_increment,1,5).LT.1) THEN
                    Z_CONT_LIST(nd_increment,1,5)=1
                    Z_CONT_LIST(nd_increment,2,5)=1
                  ENDIF

C*** 22/02/08 JHC Z_CONT_LIST(nd,j,6) contains the total contact points used to calculate 
C                 integration weights
                  Z_CONT_LIST(nd_increment,1,6)=CONT_WEIGHT
                  Z_CONT_LIST(nd_increment,2,6)=CONT_WEIGHT

                  LD(nd_increment)=ne
                  FD(nd_increment)=nef
                  ! maps local xi to the global xi
                  XID(NPF(1,nf),nd_increment)=XI_1
                  XID(NPF(3,nf),nd_increment)=XI_2
                  XID(NNF(1,nef,nb),nd_increment)=XI_3
  
C new JHC/MPN 29-OCT-2004 Select subset of contact points within faces
                  IF(direction.NE.0) THEN
                    nb=NBJ(direction,ne)
                    coord_value=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                INP(1,1,nb),nb,1,XID(1,nd_increment),
     &                XE(1,direction))
                    IF(coord_value.GE.LIMIT) THEN
                      nd_increment=nd_increment-1
                    ENDIF
                  ENDIF
c end new

C*** 22/02/08 JHC when ADD comment is used, it is necessary to search for contact points that are geometrically
C                 the same so that Lag Mult values can be updated with the ones that are previously calculated 
C                 in fe50/CONTACT_RESIDUAL.f. To do so, set the flag Z_CONT_LIST(nd,i,7) to be the existing contact point number
                  IF(BOUNDARY_CONT) THEN 
                    DO j=1,NJ_LOC(NJL_GEOM,0,nr) ! x,y,z
                      nb=NBJ(j,ne)
                      ZDD(j,nd_increment)=PXI(IBT(1,1,nb),
     &                  IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     &                  XID(1,nd_increment),XE(1,j))
                    ENDDO
                    IF(Z_CONT_LIST(nd_increment,1,5).EQ.1) THEN
                      Z_CONT_LIST(nd_increment,1,7)=0
                      Z_CONT_LIST(nd_increment,2,7)=0
                      IF (NEW_CONT_PTS) THEN
                        DO j=1,NDT ! loop over all contact points defined previously
                          IF (ABS(ZDD(1,nd_increment)-ZDD(1,j))
     &                      .LE.1.0d-8.AND.
     &                      ABS(ZDD(2,nd_increment)-ZDD(2,j))
     &                      .LE.1.0d-8.AND.
     &                      ABS(ZDD(3,nd_increment)-ZDD(3,j))
     &                      .LE.1.0d-8) THEN ! found the identical data point
                            Z_CONT_LIST(nd_increment,1,7)=j ! save which data it is identical to
                            Z_CONT_LIST(nd_increment,2,7)=j
                          ENDIF
                        ENDDO
                      ENDIF
                    ENDIF
                  ENDIF

                  IF(.NOT.BOUNDARY_CONT) THEN
                    XI_2=XI_2+(1.0d0/DBLE(NDFXI(2)))
                  ELSE
                    XI_2=XI_2+1.0d0/(DBLE(NDFXI(2))-1.0d0)
                  ENDIF
C*** 22/02/08 JHC update NDP
                  NDP(nd_increment)=nd_increment

                ENDDO !nd2

                IF(.NOT.BOUNDARY_CONT) THEN
                  XI_1=XI_1+(1.0d0/DBLE(NDFXI(1)))
                  XI_2=(1.0d0/DBLE(NDFXI(2)))*0.5d0
                ELSE
                  XI_1=XI_1+1.0d0/(DBLE(NDFXI(1))-1.0d0)
                  XI_2=0.0d0
                ENDIF
              ENDDO !nd1
            ENDDO !nf_global

!           set global variables
            CALC_XI=.TRUE.
            NDT=nd_increment

          ENDIF

        ENDIF !filio/calc
        CALL_DATA=.TRUE.

C LKC 5-JUL-1999
        CALL_XI=.TRUE.
      ENDIF
      
      CALL EXITS('DEXI')
      RETURN
 9999 CALL ERRORS('DEXI',ERROR)
      CALL EXITS('DEXI')
      RETURN 1
      END

