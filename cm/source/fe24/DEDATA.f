      SUBROUTINE DEDATA(IBT,IDO,INP,ISDATA,ISEG,LD,LDR,LN,MXI,NAN,
     '  NBJ,NBH,NDADJ,NDDATA,NDDL,NDLT,NDP,NEELEM,NELIST,NENP,NHE,NHP,
     '  NKH,NKHE,NKJ,NKJE,NLL,NLS_NDATA_CONT,NLS_NDATA_IMAG,
     '  NP_INTERFACE,NPF,NPLIST,NPNE,NPNODE,NRE,NRLIST,
     '  NVHE,NVHP,NVJE,NW,NXI,NYNE,NYNP,
     '  CE,CURVCORRECT,DET,DL,DRDN,PG,RAD,RD,RG,SE,WD,WG,XA,XE,
     '  XG,XG1,XID,XIG,XN,XP,XR,XQ,YD,YP,ZA,Z_CONT,ZD,ZE,ZF,ZP,
     '  CSEG,STRING,ERROR,*)

C#### Subroutine: DEDATA
C###  Description:
C###    DEDATA defines geometric data points.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b10.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'data00.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grou00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ioda00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'map000.cmn'
      INCLUDE 'sign00.cmn'
      INCLUDE 'time00.cmn'
      INCLUDE 'tree00.cmn'
!     Parameter List
C*** Pass NRLIST
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISDATA(NWM,NGRSEGM),ISEG(*),LD(NDM),LDR(0:NDM),
     '  LN(0:NEM),MXI(2,NEM),
     '  NBJ(NJM,NEM),NAN(NIM,NAM,NBFM),
     '  NBH(NHM,NCM,NEM),NDADJ(6,NDM),NDDATA(0:NDM,0:NRM),
     '  NDDL(NEM,NDEM),NDLT(NEM),NDP(NDM),NEELEM(0:NE_R_M,0:NRM),
     '  NELIST(0:NEM),NENP(NPM,0:NEPM,0:NRM),NHE(NEM,NXM),
     '  NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),
     '  NLS_NDATA_CONT(NDM),NLS_NDATA_IMAG(NDM),NP_INTERFACE(0:NPM,0:3),
     '  NPF(9,NFM),NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NRE(NEM),NRLIST(0:NRM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM),NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM,NXM),CURVCORRECT(2,2,NNM,NEM),
     '  DET(NBFM,0:NNM,NGM,6),DL(3,NLM),DRDN(NGM),
     '  PG(NSM,NUM,NGM,NBM),RAD(NGM),RD(NGM),
     '  RG(NGM),SE(NSM,NBFM,NEM),WD(NJM,NDM),
     '  WG(NGM,NBM),XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),
     '  XG1(NJM,NUM,NGM),XID(NIM,NDM),XIG(NIM,NGM,NBM),
     &  XN(NJM,NGM),XP(NKM,NVM,NJM,NPM),XR(NJM,NGM),XQ(NJM,NQM),YD(NHM),
     '  YP(NYM,NIYM,NXM),ZA(NAM,NHM,NCM,NEM),Z_CONT(NDM,2,67),
     '  ZD(NJM,NDM),ZE(NSM,NHM),ZF(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(*)

!     Local Variables
      EXTERNAL DDOT
      INTEGER AXIS_XI, AXIS_DIV,CERROR(50),CIRC_XI,ERR,i,
     '  I2,IBEG,IBEG1,IBEG2,ICHAR,ICOORD,
     '  IEND,IEND1,IEND2,INDEX,INDEX_OLD,INDEX_POLYLINE,
     '  INDEX_POLYMARKER,INFO,INSTAT,IOSTAT,iw,IW_DEFAULT,IWK(6),
     '  J,l,MAXIT,N,nb,N3CO,NCENTRE_REG,nd,ND0,ND1,ne,nee,ng,ni,nj,NJTT,      
     '  NJTtemp,nl,nnd,noelem,noiw,nolist,NOQUES,np,NPOINTS_PER_ELEM,nq,
     '  nr,nrr,NSTEP,NSURFACE_REG,NT_FIELD,NTIL,NTIW,NUMTIMEDATA,nv,nx,
     '  START_ELEMENT
      REAL*8 ABS_RN,ABS_U,ABS_P,ALFA,ALFA_VECTOR(3),ANGLE,B(3),
     '  BETA,BETA_DASH,BETA_pu(2),BETA_VECTOR(3),C(3),CC,
     '  COS_BETA,DDOT,DEF(3),density,DIST,DSDXI,DSTEP,DX,DXI,
     '  dX_dXi1,dX_dXi2,dY_dXi1,dY_dXi2,dZ_dXi1,dZ_dXi2,F,FL,FUNC,
     '  G1_VECTOR(3),G2_VECTOR(3),GP(3),GPD(3),GW(3),
     '  MINDIST,PXI,RADIUS_OF_INTEREST,RANGE,REAL_TEMP(1),RN_VECTOR(3),
     '  RTSEC,S,SIN_BETA,SIGNALMAX(9),
     '  SIGNALMIN(9),SIGNAL_TIME,STEP,SWAP,SXI,THETA,tmp,
     &  TOLERANCE,U_VECTOR(3),UNDEF(3),P_VECTOR(3),VAL,X,X1,X2,XD(3),
     &  X_DISCRETISATION,XI(3),XID1,XID2,XL,XPFP(3),
     '  XWC1,XWC2,XXI,XXI2(1),Y,YWC1,YWC2
      CHARACTER ANGLE_TYPE*7,ASSERT_STRING*255,CHAR1*1,CIW*1,
     '  FILE*(MXCH),FILEFORM*10,FILEFORMAT*6,SPREAD_TYPE*7,STATUS*3,
     &  TITLE*80,TYPE*20,CLABEL*52
      LOGICAL ABBREV,ALL_REGIONS,BEM,CALCU,CBBREV,
     '  COMMAND_IP,CONTINUE,CONVERG,ERROR_F,
     '  ENDFILE,FILEIP,FILIO,FIRST,GENER,Group_name_defined,
     &  ISBINFILEOPEN,MOUSE,OK,SIGNAL_VALUES,UNDEFORMED,VIEW

!     Functions
      INTEGER IFROMC
      REAL*8 RFROMC

      DATA GP/ 0.11270D0, 0.50000D0, 0.88720D0 /,
     '  GW/ 0.27777D0, 0.44444D0, 0.27777D0 /


      CALL ENTERS('DEDATA',*9999)

      ICHAR=999

C ??? Why is nx used here ???
      nx=1 ! temporary cpb 22/11/94
      nv=1 !temporary

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)
        WRITE(CIW,'(I1)') 2*NJT-3+IMAP

C---------------------------------------------------------------------

C#### Command: FEM define data;l/p/r/w<;FILENAME[$current]<;(PATH/example)[$current]>
C###  Description:
C###    Define data point values. The data is read from or written
C###    to the file FILENAME (extension .ipdata) in the directory
C###    specified by PATH.  This command enters data of type geometry,
C###    fibre, field or sheet. Each line of geometry data consists of an
C###    integer data point number, followed by a value for each
C###    coordinate (2 or 3 dimensions), followed by a weight
C###    for each coordinate.  Fibre and field data has an additional
C###    value for the fibre direction angle, or other field value,
C###    which follows the coordinate values, and an additional weight
C###    corresponding to that field value, which follows the other
C###    weights.
C###  Parameter:      <(geometry/fibre/field/sheet/contact)[geometry]>
C###    Specify the type of data to define.
C###    10/10/08 JHC added new option "contact" to write out contact 
C###                 forces. Geometric coordinates of the contact data
C###                 data points are not written out as they are re-calculated.
C###  Parameter:      <as LABEL>
C###    Specify that the data is to be read into the data group of
C###    name LABEL.
C###  Parameter:      <append>
C###    Rather than clearing the data set, data points are
C###    added to the end.
C###  Parameter:      <(radians/degrees)[radians]>
C###    Specify whether or not the angle values are in degrees or
C###    radians.
C###  Parameter:      <num_fields [0]>
C###    Number of field values to read in (over-rides
C###    the number of defined fields)
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to associate data with.

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
C*** 10/10/08 JHC added new parameter for contact
C        OP_STRING(2)=BLANK(1:15)
C     '    //'<(geometry/fibre/field/sheet)[geometry]>'
        OP_STRING(2)=BLANK(1:15)
     &    //'<(geometry/fibre/field/sheet/contact)[geometry]>'
        OP_STRING(3)=BLANK(1:15)//'<as LABEL>' !new AAY 7Dec95
        OP_STRING(3)=BLANK(1:15)//'<append>'
        OP_STRING(4)=BLANK(1:15)//'<(radians/degrees)[radians]>'
        OP_STRING(5)=BLANK(1:15)//'<(num_fields)[0]>'
        OP_STRING(6)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define data;m
C###  Description:
C###    Define the data information using the mouse. The mouse is used
C###    to locate the positions of the data, and any extra values at
C###    the location is entered via the keyboard.
C###  Parameter:      <(geometry/fibre/field/sheet)[geometry]>
C###    Specify the type of data to define.
C###  Parameter:      <on WS#[1]>
C###    Specify the workstation number on which to define data.
C###  Parameter:      <(view/noview)[view]>
C###    Specify whether to display any information to screen
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to associate data with.

        OP_STRING(1)=STRING(1:IEND)//';m'
        OP_STRING(2)=BLANK(1:15)
     '    //'<(geometry/fibre/field/sheet)[geometry]>'
        OP_STRING(3)=BLANK(1:15)//'<on WS#['//CIW//']>'
        OP_STRING(4)=BLANK(1:15)//'<(view/noview)[view]>'
        OP_STRING(5)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define data;r<;FILENAME[$current]><;(PATH/example)[$current]> planes
C###  Description:
C###    Read the image plane numbers for the data points. The data is
C###    read to the file FILENAME (extension .iplada) in the directory
C###    specified by PATH.
C###  Parameter:      <data GROUP_NAME>
C###    Specify the data group name to read the data points into.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to associate data with.

        OP_STRING(1)=STRING(1:IEND)//';r'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']> planes'
        OP_STRING(2)=BLANK(1:15)//'<data GROUP_NAME>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define data;c cat
C###  Description:
C###    ???
C###  Parameter:      <element (#s/all)[all]>
C###    Specify the element numbers to used. The "all" keyword will
C###    use all currently defined elements in the given regions.
C###    Specify the element numbers to used. The "all" keyword will
C###    use all currently defined elements in the given regions.
C###  Parameter:      <step STEPSIZE#[5]>
C###    Specify the step size to use.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.

        OP_STRING(1)=STRING(1:IEND)//';c cat'
        OP_STRING(3)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<step STEPSIZE#[5]>'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define data;m geometry BEM
C###  Description:
C###    Define using the mouse the geometric location of the data
C###    point at which to calculate the boundary element domain
C###    solution.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to associate data with.

        OP_STRING(1)=STRING(1:IEND)//';m geometry BEM'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define data;r<;FILENAME[$current]><;(PATH/example)[$current]> signal
C###  Description:
C###    Read the electrode(data) positions from a signal file. The
C###    signal data is read from the file FILENAME (extension .ipsign
C###    or .binsig depending on whether or not the signal file is an
C###    ascii or binary file) in the directory specified by PATH.
C###  Parameter:     <values <time #[0.0]>>
C###    Specify that the the electode values (e.g. potentials) are also
C###    read. The time at which the values are read may also be
C###    specified.
C###  Parameter:     <(ascii/binary)[ascii]>
C###    Specify whether or not the signal file is a ascii or binary
C###    file.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to associate data with.

        OP_STRING(1)=STRING(1:IEND)//';r'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']> signal'
        OP_STRING(2)=BLANK(1:15)//'<values <time #[0.0]>>'
        OP_STRING(3)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define data;c sheet
C###  Description:
C###    ????

        OP_STRING(1)=STRING(1:IEND)//';c sheet'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define data;c stripe
C###  Description:
C###    Calculate data points along a SPAMM stripe (OBSELETE).
C###  Parameter:     <group NAME>
C###    Specify the group name to use.
C###  Parameter:     <total #DATA_POINTS>
C###    Specfiy the total number of data points to calculate.


        OP_STRING(1)=STRING(1:IEND)
     '    //';c stripe group NAME total #DATA_POINTS'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define data;r/w<;FILENAME[$current]><;(PATH/example)[$current]> electrode
C###  Description:
C###    Define electrode positions and data. OBSOLETE.
C###  Parameter:      <file_format (activation/signal/ipfile/map3d)[ipfile]>

        OP_STRING(1)=STRING(1:IEND)//';r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']> electrode'
        OP_STRING(2)=BLANK(1:15)
     '    //'<file_format (activation/emap/ipfile/map3d)[ipfile]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define data;c grid
C###  Description:
C###    Calcualte the nearest grid points to each electrode (data
C###    point).

        OP_STRING(1)=STRING(1:IEND)//';c grid'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define data;c from_grid
C###  Description:
C###    Create a data point for each grid point
C###  Parameter:      <num_fields [0]>
C###    Number of field values to define

        OP_STRING(1)=STRING(1:IEND)//';c from_grid'
        OP_STRING(2)=BLANK(1:15)//'<(num_fields)[0]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define data
C###  Description:
C###    Define data values using command line values.
C###  Parameter:      <number #[NDT+1]>
C###    Specify the number of the data point to define. Defaults to
C###    to the next number past the highest defined data point number.
C###  Parameter:      <x=XCOORD#[0.0]>
C###    Specify the x coordinate value of the data point being defined.
C###  Parameter:      <y=YCOORD#[0.0]>
C###    Specify the y coordinate value of the data point being defined.
C###  Parameter:      <z=ZCOORD#[0.0]>
C###    Specify the z coordinate value of the data point being defined.
C###  Parameter:      <fibre_angle=FIBRE_ANGLE#>
C###    Specify the value of the fibre angle field (if any) at the data
C###    point being defined.
C###  Parameter:      <x_tangent=XTANGENT#>
C###    Specify the x component of the tangent vector at the data point
C###    being defined.
C###  Parameter:      <y_tangent=YTANGENT#>
C###    Specify the y component of the tangent vector at the data point
C###    being defined.
C###  Parameter:      <z_tangent=ZTANGENT#>
C###    Specify the z component of the tangent vector at the data point
C###    being defined.
C###  Parameter:      <x_normal=XNORMAL#>
C###    Specify the x component of the normal vector at the data point
C###    being defined.
C###  Parameter:      <y_normal=YNORMAL#>
C###    Specify the y component of the normal vector at the data point
C###    being defined.
C###  Parameter:      <z_normal=ZNORMAL#>
C###    Specify the z component of the normal vector at the data point
C###    being defined.
C###  Parameter:      <(degrees/radians)[degrees]>
C###    Specify whether or not the angle values are in degrees or
C###    radians.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to associate data with.

        OP_STRING(1)=STRING(1:IEND)//' <number #[NDT+1]>'
        OP_STRING(2)=BLANK(1:15)//'<x=XCOORD#[0.0]>'
        OP_STRING(3)=BLANK(1:15)//'<y=YCOORD#[0.0]>'
        OP_STRING(4)=BLANK(1:15)//'<z=ZCOORD#[0.0]>'
        OP_STRING(5)=BLANK(1:15)//'<fibre_angle=FIBRE_ANGLE#>'
        OP_STRING(6)=BLANK(1:15)//'<x_tangent=XTANGENT#>'
        OP_STRING(7)=BLANK(1:15)//'<y_tangent=YTANGENT#>'
        OP_STRING(8)=BLANK(1:15)//'<z_tangent=ZTANGENT#>'
        OP_STRING(9)=BLANK(1:15)//'<x_normal=XNORMAL#>'
        OP_STRING(10)=BLANK(1:15)//'<y_normal=YNORMAL#>'
        OP_STRING(11)=BLANK(1:15)//'<z_normal=ZNORMAL#>'
        OP_STRING(12)=BLANK(1:15)//'<(degrees/radians)[degrees]>'
        OP_STRING(13)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define data;c from_xi
C###  Description:
C###    Calculates the geometric positions of data points given
C###    the element and their xi locations. This information is
C###    obtained by previously reading a *.ipxi file
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to associate data with.
C###  Parameter:      <(undeformed/deformed) [undeformed]>
C###    Specify the geometery to use.

        OP_STRING(1)=STRING(1:IEND)//';c from_xi'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<(undeformed/deformed) [undeformed]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define data;c from_tag
C###  Description:
C###    Calculates a deformed set of tag intersection points
C###    from the current set of data which is assumed to be
C###    from an undeformed data set with the weights being the
C###    image plane normal

        OP_STRING(1)=STRING(1:IEND)//';c from_tag'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define data;c centreline 
C###  Description:
C###    Creates a line of data points along the centre of mass of a tubular
C###    surface mesh.  The data points will be generated along the axial xi 
C###    direction beginning from the start element.  
C###    The centre of mass points are calculated by iterating around the 
C###    circumferential direction using gaussian quadrature.
C###    The result will be exact if cubic elements are use in the 
C###    circumferential xi direction.
C###    The number of axis_divisions is used to determine how many points to
C###    generate for each axial element.  The generated points
C###    are stored in the same region as the surface mesh and will 
C###    overwrite any other data points in that region.
C###  Parameter:      <axis_xi_direction [1]>
C###  Parameter:      <circ_xi_direction [2]>
C###  Parameter:      <axis_divisions [1]>
C###  Parameter:      <region (#)[1]>
C###  Parameter:      <start_element [1]>
C###    Specify the region that contains the tube shaped mesh

        OP_STRING(1)=STRING(1:IEND)//';c centreline'
        OP_STRING(2)=BLANK(1:15)//'<circ_xi_direction [1]>'
        OP_STRING(3)=BLANK(1:15)//'<axis_xi_direction [2]>'
        OP_STRING(4)=BLANK(1:15)//'<axis_divisions [1]>'
        OP_STRING(5)=BLANK(1:15)//'<region (#) [1]>'
        OP_STRING(6)=BLANK(1:15)//'<start_element [1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define data;c area
C###  Description:
C###    Creates data points with an area field along a 1D cubic hermite
C###    centre line. Each area value is the cross sectional area of the
C###    mesh at the point.  The cross section is calculated as the
C###    intersection of the mesh and a plane normal to the xi derivative
C###    at the data point. The centre mesh needs to be loaded into
C###    one region and the surface mesh into another.  Data points are
C###    generated at xi=0 if the default of 1 point per element is
C###    specified. In general for n data points the xi values are
C###    (i-1)/n where i=1:n 
C###    The generated points are stored in the same region as the
C###    surface mesh and will  overwrite any other data points in that
C###    region.  The cross sectional area is calculated by an iterative
C###    scheme which approximates the area by using more and more
C###    triangles until the difference between successive values is less
C###    than the specified tolerance.
C###    For curved surfaces it may be necessary to restrict the
C###    intersection of the plane with the mesh to be within a radius of
C###    interest.  If the surface bends back on itself it is necessary
C###    to specify a radius of interest that excludes the other elements
C###    which also intersect which we do not wish to consider at that point.
C###  Parameter:      <surface_region (#) [1]>
C###  Parameter:      <centre_region (#) [2]>
C###  Parameter:      <points_per_element [1]>
C###  Parameter:      <tolerance [0.001]>
C###  Parameter:      <radius_of_interest [10]>

C###    Specify the region that contains the tube shaped mesh

        OP_STRING(1)=STRING(1:IEND)//';c area'
        OP_STRING(2)=BLANK(1:15)//'<surface_region (#) [1]>'
        OP_STRING(3)=BLANK(1:15)//'<centre_region (#) [2]>'
        OP_STRING(4)=BLANK(1:15)//'<points_per_element [1]>'
        OP_STRING(5)=BLANK(1:15)//'<tolerance [0.001]>'
        OP_STRING(6)=BLANK(1:15)//'<radius_of_interest [10]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define data;c fill_volume
C###  Description:
C###    Fills volume elements with data points.
C###  Parameter:      <elements [all]>
C###    Specify which volume elements to fill.
C###  Parameter:      <spread regular/xi/random [regular]>
C###    Specify the type of data point spread: a regular grid in global
C###    space, a regular grid in xi space, or a random distribution.
C###  Parameter:      <density/spacing [density]>
C###    Specify the density of data points (number of points per mm^3)
C###    or else specify the linear spacing (only for a regular grid).
C###  Parameter:          <VALUE [1]>
C###    Specify the value for density or linear spacing.

        OP_STRING(1)=STRING(1:IEND)//';c fill_volume'
        OP_STRING(2)=BLANK(1:15)//'<elements [all]>'
        OP_STRING(3)=BLANK(1:15)//'<spread (regular/xi/random)'
     &    //'[regular]>'
        OP_STRING(4)=BLANK(1:15)//'<density/spacing [density]>'
        OP_STRING(5)=BLANK(1:15)//'    <VALUE [1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define data;c cylinder_tree
C###  Description:
C###    Creates data points with a random distribution over the 'surface'
C###    of a 1D mesh with associated diameters.
C###  Parameter:      <nodes [all]>
C###  Parameter:      <density [1]>

        OP_STRING(1)=STRING(1:IEND)//';c cylinder_tree'
        OP_STRING(2)=BLANK(1:15)//'<nodes [all]>'
        OP_STRING(3)=BLANK(1:15)//'<density [1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define data;c from_node
C###  Description:
C###    Creates data points located at the nodes
C###  Parameter:      <nodes [all]>

        OP_STRING(1)=STRING(1:IEND)//';c from_node'
        OP_STRING(2)=BLANK(1:15)//'<nodes [all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEDATA',ERROR,*9999)
      ELSE

C LKC 10-JAN-2000 new assert
        CALL ASSERT(USE_DATA.EQ.1,'Set USE_DATA to 1',ERROR,*9999)

        COMMAND_IP=.FALSE.
        CALL PARSE_QUALIFIERS(' CDLMPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
        IF(NTCOQU(noco).EQ.0) THEN !command line input
          COMMAND_IP=.TRUE.
        ENDIF

C LKC 4-JUL-1999 Generalise for regions
        CALL ASSERT(NRLIST(0).EQ.1,
     '    '>> Only implemented for 1 region',ERROR,*9999)
        nr=NRLIST(1)

        FILEIP=.FALSE.
        NOQUES=0
        BEM =.FALSE.
        Group_name_defined=.FALSE.
        ANGLE_TYPE='RADIANS'

        IF(CBBREV(CO,'NOVIEW',3,noco+1,NTCO,N3CO)) THEN
          VIEW=.FALSE.
        ELSE
          VIEW=.TRUE.
        ENDIF

        IF(NTCO.GE.noco+1) THEN
          IF(CBBREV(CO,'CAT',1,noco+1,NTCO,N3CO)) THEN
            TYPE='CATDAT'
            IF(CBBREV(CO,'STEP',2,noco+1,NTCO,N3CO)) THEN
C LKC Need to pass in array
C              CALL PARSRL(CO(N3CO+1),1,NTIL,STEP,ERROR,*9999)
              CALL PARSRL(CO(N3CO+1),1,NTIL,REAL_TEMP,ERROR,*9999)
              STEP=REAL_TEMP(1)
            ELSE
              STEP=5.0D0
            ENDIF
            IF(DOP) THEN
              WRITE(OP_STRING,'('' stepsize='')') STEP
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ELSE IF(CBBREV(CO,'ELECTRODE',2,noco+1,NTCO,N3CO)) THEN
            TYPE='ELECTRODE'
            IF(CBBREV(CO,'FILE_FORMAT',3,noco+1,NTCO,N3CO)) THEN
              IF(ABBREV(CO(N3CO+1),'ACTIVATION',1)) THEN
                FILEFORM='ACTIVATION'
              ELSE IF(ABBREV(CO(N3CO+1),'EMAP',1)) THEN
                FILEFORM='EMAP'
              ELSE IF(ABBREV(CO(N3CO+1),'IPFILE',1)) THEN
                FILEFORM='IPFILE'
              ELSE IF(ABBREV(CO(N3CO+1),'MAP3D',1)) THEN
                FILEFORM='MAP3D'
              ENDIF
            ELSE
              FILEFORM='IPFILE'
            ENDIF
          ELSE IF(CBBREV(CO,'FIBRES',3,noco+1,NTCO,N3CO)) THEN
            TYPE='FIBRES'
C*** 10/10/08 JHC Added option for contact
          ELSE IF(CBBREV(CO,'CONTACT',5,noco+1,NTCO,N3CO)) THEN
            TYPE='CONTACT'
          ELSE IF(CBBREV(CO,'FIELD',3,noco+1,NTCO,N3CO)) THEN
            TYPE='FIELD'
C LKC 7-MAR-2002 Not compulsory
C            CALL ASSERT(NJ_LOC(NJL_FIEL,0,0).GT.0,
C     '        ' >>Field variable not defined',
C     '        ERROR,*9999)

            IF(NJ_LOC(NJL_FIEL,0,0).LE.0) THEN
              WRITE(OP_STRING(1),
     '          '('' WARNING: Field variables not defined yet '')')
              WRITE(OP_STRING(2),
     '          '('' WARNING: Setting up the field variables now '')')

              CALL WRITES(IOER,OP_STRING,ERROR,*9999)
            ENDIF
            IF(IOTYPE.EQ.3)THEN !for writing data fields
               IF(CBBREV(CO,'NUM_FIELDS',3,noco+1,NTCO,N3CO)) THEN
                  NT_FIELD=IFROMC(CO(N3CO+1))
               ELSE
                  NT_FIELD = NJ_LOC(NJL_FIEL,0,nr)
               ENDIF
            ELSE
               NT_FIELD = NJ_LOC(NJL_FIEL,0,nr)
            ENDIF

          ELSE IF(CBBREV(CO,'GEOMETRY',2,noco+1,NTCO,N3CO)) THEN
            TYPE='GEOMETRY'
            IF(CBBREV(CO,'BEM',1,noco+1,NTCO,N3CO)) THEN
              BEM=.TRUE.
            ENDIF
          ELSE IF(CBBREV(CO,'GRID',2,noco+1,NTCO,N3CO)) THEN
            TYPE='GRID'
          ELSE IF(CBBREV(CO,'FROM_GRID',9,noco+1,NTCO,N3CO)) THEN
            TYPE='FROM_GRID'

C LKC 15-NOV-98 Unused ? - No command can access this section
C
C          ELSE IF(CBBREV(CO,'HISTORY',1,noco+1,NTCO,N3CO)) THEN
C            TYPE='HISTORY'
C            IF(CBBREV(CO,'DATA',2,noco+1,NTCO,N3CO)) THEN
CC LKC 24-APR-1998 Need to pass in array
CC              CALL PARSIL(CO(N3CO+1),1,NTIL,NDHIST,ERROR,*9999)
C              CALL PARSIL(CO(N3CO+1),1,NTIL,INT_TEMP,ERROR,*9999)
C              NDHIST=INT_TEMP(1)
C            ELSE
C              NDHIST=1
C            ENDIF
C            IF(CBBREV(CO,'COORD',2,noco+1,NTCO,N3CO)) THEN
CC LKC 24-APR-1998 Need to pass in array
CC             CALL PARSIL(CO(N3CO+1),1,NTIL,NJHIST,ERROR,*9999)
C              CALL PARSIL(CO(N3CO+1),1,NTIL,INT_TEMP,ERROR,*9999)
C              NJHIST=INT_TEMP(1)
C            ELSE
C              NJHIST=1
C            ENDIF
C            IF(DOP) THEN
C              WRITE(OP_STRING,'('' history data at ndhist,njhist '','
C     '          //'I5,I5)') ndhist,njhist
C              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C            ENDIF
C LKC 15-NOV-98 Unused ? - No command can access this section
C          ELSE IF(CBBREV(CO,'NUMBERS',1,noco+1,NTCO,N3CO)) THEN
C            TYPE='NUMBERS'
C          ELSE IF(CBBREV(CO,'PARAMETER',2,noco+1,NTCO,N3CO)) THEN
C            TYPE='PARAMETER'
C          ELSE IF(CBBREV(CO,'PROJECTIONS',2,noco+1,NTCO,N3CO)) THEN
C            TYPE='PROJECTIONS'

          ELSE IF(CBBREV(CO,'VALUES',1,noco+1,NTCO,N3CO)) THEN
            TYPE='VALUES'
          ELSE IF(CBBREV(CO,'SHEETS',2,noco+1,NTCO,N3CO)) THEN
            TYPE='SHEETS'
          ELSE IF(CBBREV(CO,'SIGNAL',2,noco+1,NTCO,N3CO)) THEN
            TYPE='SIGNAL'
            IF(.NOT.CALCU) THEN
              SIGNAL_VALUES=.FALSE.
              IF(CBBREV(CO,'VALUES',2,noco+1,NTCO,N3CO)) THEN
                SIGNAL_VALUES=.TRUE.
                IF(ABBREV(CO(N3CO+1),'TIME',2)) THEN
                  IF(N3CO+2.LE.NTCO) THEN
                    SIGNAL_TIME=RFROMC(CO(N3CO+2))
                  ELSE
                    ERROR='>>No signal time specified'
                    GOTO 9999
                  ENDIF
                ELSE
                  SIGNAL_TIME=0.0d0
                ENDIF
              ENDIF
              IF(CBBREV(CO,'BINARY',2,noco+1,NTCO,N3CO)) THEN
                FILEFORMAT='BINARY'
              ELSE
                FILEFORMAT='ASCII'
              ENDIF
            ENDIF

C LKC 15-NOV-98 Unused ? - No command can access this section
C          ELSE IF(CBBREV(CO,'SLICE',2,noco+1,NTCO,N3CO)) THEN
C            TYPE='SLICE'

          ELSE IF(CBBREV(CO,'STRIPE',2,noco+1,NTCO,N3CO)) THEN
            TYPE='STRIPE'

C LKC 15-NOV-98 Unused ? - No command can access this section
C          ELSE IF(CBBREV(CO,'TRACE',1,noco+1,NTCO,N3CO)) THEN
C            TYPE='TRACE'
          ELSE IF(CBBREV(CO,'PLANES',4,NOCO+1,NTCO,N3CO)) THEN
            TYPE='PLANES'
            CALL PARSE_DATA(nd0,nd1,noco,NTCO,CO,ERROR,*9999)
          ELSE IF(CBBREV(CO,'FROM_XI',7,NOCO+1,NTCO,N3CO)) THEN
            TYPE='FROM_XI'
            IF(CBBREV(CO,'DEFORMED',2,noco+1,NTCO,N3CO)) THEN
              UNDEFORMED=.FALSE.
            ELSE
              UNDEFORMED=.TRUE.
            ENDIF
          ELSE IF(CBBREV(CO,'FROM_TAG',7,NOCO+1,NTCO,N3CO)) THEN
             TYPE='FROM_TAG'
          ELSE IF(CBBREV(CO,'CENTRELINE',2,NOCO+1,NTCO,N3CO)) THEN
            TYPE='CENTRELINE'
            IF(CBBREV(CO,'AXIS_XI_DIRECTION',6,noco+1,
     '        NTCO,N3CO)) THEN
              AXIS_XI = IFROMC(CO(N3CO+1))
            ELSE
              AXIS_XI = 2
            ENDIF
            WRITE(*,'(''Axis direction: xi '', I1)') AXIS_XI
            IF(CBBREV(CO,'CIRC_XI_DIRECTION',6,noco+1,
     '        NTCO,N3CO)) THEN
              CIRC_XI = IFROMC(CO(N3CO+1))
            ELSE
              CIRC_XI = 1
            ENDIF
            WRITE(*,'(''Circumference direction: xi '', 
     '        I1)') CIRC_XI
            IF(CBBREV(CO,'AXIS_DIVISIONS',6,noco+1,
     '        NTCO,N3CO)) THEN
              AXIS_DIV = IFROMC(CO(N3CO+1))
            ELSE
              AXIS_DIV = 1
            ENDIF
            WRITE(*,'(''Number of axial divisions:'', 
     '        I2)') AXIS_DIV
            IF(CBBREV(CO,'REGION',1,noco+1,NTCO,N3CO)) THEN
              nr=IFROMC(CO(N3CO+1)) !set region
            ENDIF !REGION      
            WRITE(*,'(''Region: '',I1)' ) nr 
            IF(CBBREV(CO,'START_ELEMENT',1,
     '        noco+1,NTCO,N3CO)) THEN
              START_ELEMENT=IFROMC(CO(N3CO+1))
            ELSE
              START_ELEMENT = 1
            ENDIF !REGION      
            WRITE(*,'(''Start element: '',I3)') 
     '        START_ELEMENT            
          ELSE IF(CBBREV(CO,'AREA',2,NOCO+1,NTCO,N3CO)) THEN
            TYPE='AREA'
            IF(CBBREV(CO,'SURFACE_REGION',2,noco+1,
     '        NTCO,N3CO)) THEN
              NSURFACE_REG = IFROMC(CO(N3CO+1))
            ELSE
              NSURFACE_REG = 1
            ENDIF !surface_region
            WRITE(*,'(''Surface region '', I1)') NSURFACE_REG
            IF(CBBREV(CO,'CENTRE_REGION',2,noco+1,
     '        NTCO,N3CO)) THEN
              NCENTRE_REG = IFROMC(CO(N3CO+1))
            ELSE
              NCENTRE_REG = 2
            ENDIF !centre_region
            WRITE(*,'(''Centre region '',I1)') NCENTRE_REG
            IF(CBBREV(CO,'POINTS_PER_ELEMENT',2,noco+1,
     '        NTCO,N3CO)) THEN
              NPOINTS_PER_ELEM = IFROMC(CO(N3CO+1))
            ELSE
              NPOINTS_PER_ELEM = 1
            ENDIF !points_per_element
            WRITE(*,'(''Points per element'',I2)') NPOINTS_PER_ELEM
            IF(CBBREV(CO,'TOLERANCE',2,noco+1,NTCO,N3CO)) THEN
              TOLERANCE=RFROMC(CO(N3CO+1))
            ELSE
              TOLERANCE=0.001d0
            ENDIF !tolerance
            WRITE(*,'(''Tolerance: '',G12.5)' ) TOLERANCE
            IF(CBBREV(CO,'RADIUS_OF_INTEREST',2,
     '        noco+1,NTCO,N3CO)) THEN
              RADIUS_OF_INTEREST=RFROMC(CO(N3CO+1))
            ELSE
              RADIUS_OF_INTEREST = 10.0d0
            ENDIF !radius_of_interest
            WRITE(*,'(''Radius of interest: '',G12.5)') 
     '        RADIUS_OF_INTEREST  
          ELSE IF(CBBREV(CO,'CYLINDER_TREE',2,NOCO+1,NTCO,N3CO)) THEN
            TYPE='CYLINDER_TREE'
            FIRST=.TRUE.
            CALL PARSE_NODES(NPNODE,NPLIST,N3CO,NRLIST,NTCO,CO,
     '        ERROR,*9999)
            IF(CBBREV(CO,'DENSITY',3,noco+1,NTCO,N3CO)) THEN
              density = RFROMC(CO(N3CO+1))
            ELSE
              density=1.d0
            ENDIF
            IF(CBBREV(CO,'NOT_FIRST',3,noco+1,NTCO,N3CO)) FIRST=.FALSE.
          ELSE IF(CBBREV(CO,'FROM_NODE',2,NOCO+1,NTCO,N3CO)) THEN
            TYPE='FROM_NODE'
            CALL PARSE_NODES(NPNODE,NPLIST,N3CO,NRLIST,NTCO,CO,
     '        ERROR,*9999)
          ELSE IF(CBBREV(CO,'FILL_VOLUME',4,NOCO+1,NTCO,N3CO)) THEN
            TYPE='FILL_VOLUME'
            CALL PARSE_ELEMENTS(NEELEM,NELIST,N3CO,NRLIST,NTCO,CO,
     '        ERROR,*9999)
C           spread regular/xi/random [regular]
            IF(CBBREV(CO,'SPREAD',3,noco+1,NTCO,N3CO)) THEN
              IF(CBBREV(CO,'XI',2,noco+1,NTCO,N3CO)) THEN
                SPREAD_TYPE='XI'
              ELSE IF(CBBREV(CO,'RANDOM',2,noco+1,NTCO,N3CO)) THEN
                SPREAD_TYPE='RANDOM'
              ELSE
                SPREAD_TYPE='REGULAR'
              ENDIF
            ELSE
              SPREAD_TYPE='REGULAR'
            ENDIF
            IF(CBBREV(CO,'SPACING',3,noco+1,NTCO,N3CO)) THEN
              X_DISCRETISATION = RFROMC(CO(N3CO+1))
            ELSE IF(CBBREV(CO,'DENSITY',3,noco+1,NTCO,N3CO)) THEN
              X_DISCRETISATION= (1.d0/RFROMC(CO(N3CO+1)))**(1.d0/3.d0)
            ENDIF
          ELSE
            TYPE='GEOMETRY'
          ENDIF

          IF(CBBREV(CO,'AS',2,noco+1,NTCO,N3CO)) THEN
            NTGRDA=NTGRDA+1
            CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
            LAGRDA(NTGRDA)=CO(N3CO+1)(IBEG:IEND)
            Group_name_defined=.TRUE.
          ENDIF
          ADD=.FALSE.
          IF(CBBREV(CO,'APPEND',6,noco+1,NTCO,N3CO)) THEN
            ADD=.TRUE.
          ENDIF
        ELSE IF(NTCO.EQ.noco) THEN
          TYPE='GEOMETRY'
        ELSE
          TYPE='GEOMETRY'
        ENDIF

        IF(TYPE(1:6).EQ.'FIBRES'.OR.TYPE(1:6).EQ.'SHEETS') THEN
          IF(CBBREV(CO,'RADIANS',3,noco+1,NTCO,N3CO)) THEN
            ANGLE_TYPE='RADIANS'
          ELSE IF(CBBREV(CO,'DEGREES',3,noco+1,NTCO,N3CO)) THEN
            ANGLE_TYPE='DEGREES'
          ELSE
            ANGLE_TYPE='RADIANS'
          ENDIF
        ENDIF

        IF(TYPE(1:5).EQ.'SLICE') THEN
          IF(CBBREV(CO,'X',1,noco+1,NTCO,N3CO)) THEN
            nj=1
          ELSE IF(CBBREV(CO,'Y',1,noco+1,NTCO,N3CO)) THEN
            nj=2
          ELSE IF(CBBREV(CO,'Z',1,noco+1,NTCO,N3CO)) THEN
            nj=3
          ENDIF
          VAL=RFROMC(CO(N3CO+1))
          IF(CBBREV(CO,'RANGE',1,noco+1,NTCO,N3CO)) THEN
            RANGE=RFROMC(CO(N3CO+1))
          ELSE
            RANGE=0.0D0
          ENDIF
          XI3MIN=VAL-RANGE
          XI3MAX=VAL+RANGE

        ELSE IF(TYPE(1:6).EQ.'STRIPE') THEN

          WRITE(OP_STRING,'(''OLD CODE AREA - CODE REMOVED'')')

C MLB 19 Feb 1997 Unused ????
C          IF(CBBREV(CO,'GROUP',1,noco+1,NTCO,N3CO)) THEN
C            CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
C            GROUP_NAME=CO(N3CO+1)(IBEG:IEND)
C          ELSE
C            ERROR='>>>Group name needed'
C            GO TO 9999
C          ENDIF
C          IF(CBBREV(CO,'TOTAL',1,noco+1,NTCO,N3CO)) THEN
C            NT_STRIPE_DATA=IFROMC(CO(N3CO+1))
C          ELSE
C            NT_STRIPE_DATA=2000
C          ENDIF
        ENDIF

        IF(TYPE(1:8).EQ.'GEOMETRY'.OR.TYPE(1:5).EQ.'FIELD'.OR
     '    .TYPE(1:7).EQ.'NUMBERS'.OR.TYPE(1:6).EQ.'VALUES'.OR
     '    .TYPE(1:11).EQ.'PROJECTIONS'.OR.TYPE(1:6).EQ.'FIBRES'.OR
     '    .TYPE(1:5).EQ.'TRACE'.OR.TYPE(1:9).EQ.'PARAMETER'.OR
     '    .TYPE(1:6).EQ.'CATDAT'.OR
     '    .TYPE(1:6).EQ.'SHEETS'.OR.TYPE(1:7).EQ.'CONTACT') THEN
          CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '      ERROR,*9999)
        ENDIF

        IF(MOUSE) THEN
          IF(TYPE(1:8).EQ.'GEOMETRY'.OR.TYPE(1:5).EQ.'FIELD') THEN
            CALL STRING_TRIM(TYPE,IBEG,IEND)
            IF(CBBREV(CO,'DEFORMED',3,noco+1,NTCO,N3CO)) THEN
              CALL ASSERT(NJ_LOC(NJL_FIEL,0,1).EQ.NJT,
     '          '>>no deformed coords',ERROR,*9999)
              TYPE=TYPE(IBEG:IEND)//'_DEFORMED'
            ENDIF
            IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
              INDEX=INDEX_POLYMARKER(0,'PLUS','SIZE1',CO(N3CO+1))
            ELSE
              INDEX=INDEX_POLYMARKER(0,'PLUS','SIZE1','BLUE')
            ENDIF
          ELSE IF(TYPE(1:6).EQ.'FIBRES'.OR.TYPE(1:6).EQ.'SHEETS') THEN
            IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
              INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
            ELSE
              INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLUE')
            ENDIF
          ELSE IF(TYPE(1:5).EQ.'TRACE') THEN
            IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
              INDEX=INDEX_POLYMARKER(0,'CIRCLE','SIZE1',CO(N3CO+1))
            ELSE
              INDEX=INDEX_POLYMARKER(0,'CIRCLE','SIZE1','CYAN')
            ENDIF
          ELSE IF(TYPE(1:11).EQ.'PROJECTIONS') THEN
            IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
              INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
            ELSE
              INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','CYAN')
            ENDIF
          ENDIF
          IF(TYPE(1:8).EQ.'GEOMETRY'.OR.TYPE(1:5).EQ.'FIELD'.OR
     '      .TYPE(1:5).EQ.'TRACE'.OR.TYPE(1:6).EQ.'FIBRES'
     '      .OR.TYPE(1:11).EQ.'PROJECTIONS') THEN
            IF(CBBREV(CO,'FROM',2,noco+1,NTCO,N3CO)) THEN
              XI3MIN=RFROMC(CO(N3CO+1))
            ELSE
              XI3MIN=0.0D0
            ENDIF
            IF(CBBREV(CO,'TO',2,noco+1,NTCO,N3CO)) THEN
              XI3MAX=RFROMC(CO(N3CO+1))
            ELSE
              XI3MAX=1.0d0
            ENDIF
            IF(CBBREV(CO,'WHEN',2,noco+1,NTCO,N3CO)) THEN
              TIME=RFROMC(CO(N3CO+1))
            ELSE
              TIME=0.0D0
            ENDIF
          ELSE IF(TYPE(1:6).EQ.'SHEETS') THEN
            IF(CBBREV(CO,'THETA',2,noco+1,NTCO,N3CO)) THEN
              SHEET_THETA=PI*RFROMC(CO(N3CO+1))/180.0D0
            ELSE
              SHEET_THETA=0.0D0
            ENDIF
            IF(CBBREV(CO,'RANGE',1,noco+1,NTCO,N3CO)) THEN
              SHEET_RANGE=RFROMC(CO(N3CO+1))*PI/180.0D0
            ELSE
C CPB 21/3/93 Should this be PI/180.0d0 ????? NO (PJH)
              SHEET_RANGE=PI/18.0d0
            ENDIF
          ENDIF
        ENDIF

!Set default workstation numbers
        IF(MOUSE) THEN
          IW_DEFAULT=2*NJT-3+IMAP
        ENDIF

        IF(MOUSE) THEN
          CALL WS_LIST(IWK,IW_DEFAULT,NTIW,noco,NTCO,CO,ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,*)' iwk(1..NTIW)=',
     '        (IWK(noiw),noiw=1,NTIW)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF

        IF(COMMAND_IP) THEN
          nr=NRLIST(1)
          IF(CBBREV(CO,'NUMBER',1,noco+1,NTCO,N3CO)) THEN
            nd=IFROMC(CO(N3CO+1))
          ELSE
            nd=NDT+1
          ENDIF
          NDT=NDT+1
          NDP(nd)=nd
          IF(CBBREV(CO,'X',1,noco+1,NTCO,N3CO)) THEN
            ZD(1,nd)=RFROMC(CO(N3CO+1))
          ELSE
            ZD(1,nd)=0.0d0
          ENDIF
          IF(CBBREV(CO,'Y',1,noco+1,NTCO,N3CO)) THEN
            ZD(2,nd)=RFROMC(CO(N3CO+1))
          ELSE
            ZD(2,nd)=0.0d0
          ENDIF
          IF(CBBREV(CO,'Z',1,noco+1,NTCO,N3CO)) THEN
            ZD(3,nd)=RFROMC(CO(N3CO+1))
          ELSE
            ZD(3,nd)=0.0d0
          ENDIF
!News AJP 7/4/95 Enter fibre angles from front end
          IF(CBBREV(CO,'FIBRE_ANGLE',7,noco+1,NTCO,N3CO)) THEN
            NJTT=NJT+1
            ZD(NJT+1,nd)=RFROMC(CO(N3CO+1))
            IF(.NOT.CBBREV(CO,'RADIANS',1,noco+1,NTCO,N3CO)) THEN
!Default input is degrees, storage is radians
              ZD(NJT+1,nd)=ZD(NJT+1,nd)*PI/180.0d0
            ENDIF
          ELSE
            NJTT=NJT
            ZD(NJT+1,nd)=0.0d0
          ENDIF
!Newe
          DO nj=1,NJTT
            WD(nj,nd)=1.0d0
          ENDDO

        ELSE IF(FILIO) THEN
C CPB 8/4/94 This needs to be generalised for NJ_LOC
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IF(TYPE(1:6).NE.'PLANES'
     '      .AND.TYPE(1:9).NE.'ELECTRODE'
     '      .AND.TYPE(1:6).NE.'SIGNAL') THEN
C           IF(iotype.NE.3) THEN
C             NTDATA=NTDATA+1
C             DATYPE(NTDATA)=TYPE
C             NDTOLD=0
C           ENDIF
            IF(IOTYPE.GT.1) THEN
              CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.ipdata',
     '          STATUS,'SEQUEN','FORMATTED',132,ERROR,*9999)
            ENDIF

            IF(TYPE(1:8).EQ.'GEOMETRY'.OR.TYPE(1:6).EQ.'FIBRES'.OR
     '        .TYPE(1:6).EQ.'SHEETS'.OR.TYPE(1:5).EQ.'FIELD'.OR
     &        .TYPE(1:7).EQ.'CONTACT') THEN

              IF(IOTYPE.EQ.1) THEN
                IF(ITYP10(1).GT.1) THEN
                  FORMAT=
     '              '($,'' Specify whether coords entered as (1)rect.'
     '              //' cart. or (2)curvilinear [1]: '')'
                  CALL GINOUT(IOTYPE,IPINTE,IVDU,0,0,0,NOQUES,FILEIP,
     &              FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,
     &              1,2,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     '              *9999)
C MLB 19/2/97 ICOORD does not appear to be used after is is set here
C it is in the OMIT file, please remove if this code is changed.
                  ICOORD=IDATA(1)
                ELSE
                  ICOORD=1
                ENDIF
                nd=0
                CONTINUE=.TRUE.
                DO WHILE(CONTINUE)
                  nd=nd+1
                  WRITE(CHAR1,'(I1)') NJT
                  FORMAT=
     '              '(/$,'' Enter '//CHAR1(1:1)//' coords [Exit]: '')'
                  RDEFLT(1)=-1.0D6
                  CALL GINOUT(IOTYPE,IPREAL,IVDU,0,0,0,NOQUES,FILEIP,
     &              FORMAT,NJT,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     &              IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,
     &              INFO,ERROR,*9999)
                  IF(RDATA(1).GT.-1.0D5) THEN
                    NDP(nd)=nd
                    IF(ITYP10(1).GE.2) THEN
                      RDATA(2)=RDATA(2)*PI/180.0D0
                    ENDIF
                    IF(ITYP10(1).GE.3) THEN
                      RDATA(3)=RDATA(3)*PI/180.0D0
                    ENDIF
                    CALL XZ(ITYP10(1),RDATA,ZD(1,nd))
                    DO nj=1,NJT
                      WD(nj,nd)=1.0D0
                    ENDDO
                    IF(TYPE(1:6).EQ.'FIBRES') THEN
                      FORMAT=
     '                  '(/$,'' Enter fibre angle in degrees [0]: '')'
                      CALL GINOUT(IOTYPE,IPREAL,IVDU,0,0,0,NOQUES,
     &                  FILEIP,FORMAT,1,
     '                  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     '                  IMAX,LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,
     '                  ERROR,*9999)
                      ZD(NJT+1,nd)=RDATA(1)*PI/180.0D0
                      WD(NJT+1,nd)=1.0D0
                    ELSE IF(TYPE(1:6).EQ.'SHEETS') THEN
                      FORMAT=
     '                  '(/$,'' Enter sheet angle in degrees [0]: '')'
                      CALL GINOUT(IOTYPE,IPREAL,IVDU,0,0,0,NOQUES,
     &                  FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &                  IDATA,IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RZERO,
     &                  -RMAX,RMAX,INFO,ERROR,*9999)
                      ZD(NJT+3,nd)=RDATA(1)*PI/180.0D0
                      WD(NJT+3,nd)=1.0D0
                    ELSE IF(TYPE(1:5).EQ.'FIELD') THEN
                      FORMAT='(/$,'' Enter field value [0]: '')'
                      CALL GINOUT(IOTYPE,IPREAL,IVDU,0,0,0,NOQUES,
     &                  FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &                  IDATA,IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,
     &                  -RMAX,RMAX,INFO,ERROR,*9999)
                      ZD(NJT+1,nd)=RDATA(1)
                      WD(NJT+1,nd)=1.0D0
                    ENDIF
                  ELSE
                    CONTINUE=.FALSE.
                    NDT=nd-1
                  ENDIF !rdata>-1x10^5
                ENDDO !while continue
!news           AAY 7Dec95 added data groups (this will replace stripe/contour)
                IF(Group_name_defined) THEN
                  LIGRDA(1,NTGRDA)=NDTOLD+1
                  LIGRDA(2,NTGRDA)=NDT
                ENDIF
              ELSE IF(IOTYPE.EQ.2.OR.IOTYPE.EQ.4) THEN !read data file
                IF(TYPE(1:6).EQ.'FIBRES') THEN
                  NJTT=NJT+1
                  WRITE(OP_STRING,'('' >> 1 Fibre angle '
     '              //'being read'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  CALL_DATA_FIBRE=.TRUE.
                ELSE IF(TYPE(1:6).EQ.'SHEETS') THEN
                  NJTT=NJT+3
                  WRITE(OP_STRING,'('' >> 1 Sheet angle '
     '              //'being read'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  CALL_DATA_SHEET=.TRUE.
                ELSE IF(TYPE(1:5).EQ.'FIELD') THEN

C LKC 4-JUL-1999 Generalise for regions
C                  IF(NJ_LOC(NJL_FIEL,0,1).GT.0) THEN
C                    NT_FIELD=NJ_LOC(NJL_FIEL,0,1)
                  IF(NJ_LOC(NJL_FIEL,0,nr).GT.0) THEN
                    NT_FIELD=NJ_LOC(NJL_FIEL,0,nr)
                  ELSE
                    NT_FIELD=1
                  ENDIF
C PJH 23Feb96     NJTT=NJT+1

C LKC 7-MAR-2002 Ability to override the number of field values
                  IF(CBBREV(CO,'NUM_FIELDS',3,noco+1,NTCO,N3CO)) THEN
                    NT_FIELD=IFROMC(CO(N3CO+1))
                    
C                    IF(NJ_LOC(NJL_FIEL,0,nr).NE.NT_FIELD) THEN
C DMAL 17 June 2003 Doesn't reduce the number of fields if they
C already exist. And warning user of increase in the number of fields
                    IF(NJ_LOC(NJL_FIEL,0,nr).LE.NT_FIELD) THEN
                      WRITE(OP_STRING,'(''>> Warning : Increasing '
     '                  //'the number of fields.'')')
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      NJTtemp = NJ_LOC(NJL_GEOM,0,nr)+
     '                NJ_LOC(NJL_FIBR,0,nr)+NT_FIELD
                      WRITE(ASSERT_STRING,'(''>> Increase NJM '
     '                  //'to '',I1)') NJTtemp
                      CALL ASSERT(NJTtemp.LE.NJM,
     '                  ASSERT_STRING,ERROR,*9999)
                      NJ_LOC(NJL_FIEL,0,nr)=NT_FIELD
                      DO nj=1,NJ_LOC(NJL_FIEL,0,nr)
C                       NJ_LOC(NJL_FIEL,nj,nr)=nj+NJ_LOC(0,0,nr)
                        NJ_LOC(NJL_FIEL,nj,nr)=nj+
     '                    NJ_LOC(NJL_GEOM,0,nr)+NJ_LOC(NJL_FIBR,0,nr)
                      ENDDO
C LKC 17-JAN-2003 Touching up the way the number fields are calculated
C                     when data points with fields are reread in.
C                    NJ_LOC(0,0,nr)=NJ_LOC(0,0,nr)+NT_FIELD
C                    NJ_LOC(0,0,0)=NJ_LOC(0,0,0)+NT_FIELD
                      NJ_LOC(0,0,nr)=NJ_LOC(1,0,nr)+NJ_LOC(2,0,nr)
     '                  +NT_FIELD
                      NJ_LOC(0,0,0)=NJ_LOC(1,0,nr)+NJ_LOC(2,0,nr)
     '                  +NT_FIELD
                      NJ_LOC(NJL_FIEL,0,0)=NJ_LOC(NJL_FIEL,0,0)+NT_FIELD
                      NJ_LOC(NJL_FIEL,0,nr)=NT_FIELD
                    ENDIF

                  ENDIF

                  NJTT=NJT+NT_FIELD
                  WRITE(OP_STRING,'('' >> '',I1,'' Field data being '
     '              //'read'')') NT_FIELD
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  CALL_DATA_FIELD=.TRUE.
C*** 10/10/08 JHC Added contact option
                ELSE IF(TYPE(1:7).EQ.'CONTACT') THEN
                  NJTT=NJT
                  WRITE(OP_STRING,'('' >> Contact data being read'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  CALL_DATA_CONT=.TRUE.
                ELSE
                  NJTT=NJT

C LKC 24-MAY-1998 Additional information when reading in
C                  WRITE(OP_STRING,'('' >> Geometric data '
C                 '              //'being read'')')


                  IF(VIEW) THEN
                    IF(NJT.EQ.1) THEN
                      WRITE(OP_STRING,
     '                  '('' >> 1D Geometric data being read'')')
                    ELSEIF(NJT.EQ.2) THEN
                      WRITE(OP_STRING,
     '                  '('' >> 2D Geometric data being read'')')
                    ELSEIF(NJT.EQ.3) THEN
                      WRITE(OP_STRING,
     '                  '('' >> 3D Geometric data being read'')')
                    ENDIF
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ENDIF !DISPLAY
                ENDIF

C*** 10/10/08 JHC Pass in Z_CONT to IODATA.f
C                CALL IODATA('READ',ANGLE_TYPE,TYPE,IFILE,NDP,NJTT,
C     '            TITLE,WD,ZD,ERROR,*9999)
                CALL IODATA('READ',ANGLE_TYPE,TYPE,IFILE,NDP,NJTT,
     '            TITLE,WD,Z_CONT,ZD,ERROR,*9999)

                IF(VIEW) THEN
                  OP_STRING(1)=TITLE
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDIF

                IF(Group_name_defined) THEN
                  LIGRDA(1,NTGRDA)=NDTOLD+1
                  LIGRDA(2,NTGRDA)=NDT
                ENDIF

              ELSE IF(IOTYPE.EQ.3) THEN !write data file
C LKC 15-NOV-1998 This is not making use of the command parameters.
C   Changing all the logical .OR.s to .ANDS.
C
C                IF(NJ_LOC(NJL_FIBR,0,1).EQ.1.OR.TYPE(1:6).EQ.
C     '            'FIBRES') THEN
C                  NJTT=NJT+1
C                  WRITE(OP_STRING,'(''>> 1 Fibre angle being '
C     '              //'written'')')
C                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C                ELSE IF(NJ_LOC(NJL_FIBR,0,1).EQ.2.OR.TYPE(1:6).EQ.
C     '              'FIBRES') THEN
C                  NJTT=NJT+2
C                  WRITE(OP_STRING,
C     '              '(''>> 1 fibre imbrecation angle being written'')')
C                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C                ELSE IF(NJ_LOC(NJL_FIBR,0,1).EQ.3.OR.TYPE(1:6).EQ.
C     '              'SHEETS') THEN
C                  NJTT=NJT+3
C                  WRITE(OP_STRING,'(''>> 1 Sheet angle being '
C     '              //'written'')')
C                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C                ELSE IF(NJ_LOC(NJL_FIEL,0,1).GT.0) THEN
C                  NJTT=NJT+NJ_LOC(NJL_FIEL,0,1)
C                  WRITE(OP_STRING,
C     '              '(''>> '',I1,'' Field value being written'')')
C     '              NJ_LOC(NJL_FIEL,0,1)
C                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C                ELSE IF(CALL_DATA_FIELD.OR.TYPE(1:5).EQ.'FIELD') THEN
C                  NJTT=NJT+1
C                  WRITE(OP_STRING,'(''>> 1 Field data '
C     '              //'being written'')')
C                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C                ELSE
C                  NJTT=NJT
C                  WRITE(OP_STRING,'(''>> Geometric data being '
C     '              //'written'')')
C                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C                ENDIF

C!!! CS 16 MAY 2001 changing this one back to OR because thats
C!!! what it should be. Don't understand why these were changed.
                IF(NJ_LOC(NJL_FIEL,0,nr).GT.0.AND.TYPE(1:5).EQ.
     '              'FIELD') THEN

C LKC 4-JUL-1999 generialise for regions
C                  NJTT=NJT+NJ_LOC(NJL_FIEL,0,1)
C MHT 23.03.11. using NT_FIELD to control number of fields output
c                  NJTT=NJT+NJ_LOC(NJL_FIEL,0,nr)
                  NJTT = NJT+NT_FIELD
c                  WRITE(OP_STRING,
c     '              '(''>> '',I1,'' Field value being written'')')
c     '              NJ_LOC(NJL_FIEL,0,nr)
                  WRITE(OP_STRING,
     '              '(''>> '',I1,'' Field value being written'')')
     '              NT_FIELD
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ELSE IF(CALL_DATA_FIELD.AND.TYPE(1:5).EQ.'FIELD') THEN
                  NJTT=NJT+NJ_LOC(NJL_FIEL,0,nr)
                  WRITE(OP_STRING,'(''>> 1 Field data '
     '              //'being written'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ELSE IF(NJ_LOC(NJL_FIBR,0,1).EQ.1.OR.TYPE(1:6).EQ.
c                IF(NJ_LOC(NJL_FIBR,0,1).EQ.1.OR.TYPE(1:6).EQ.
     '            'FIBRES') THEN
                  NJTT=NJT+1
                  WRITE(OP_STRING,'(''>> 1 Fibre angle being '
     '              //'written'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ELSE IF(NJ_LOC(NJL_FIBR,0,1).EQ.2.AND.TYPE(1:6).EQ.
     '              'FIBRES') THEN
                  NJTT=NJT+2
                  WRITE(OP_STRING,
     '              '(''>> 1 fibre imbrecation angle being written'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ELSE IF(NJ_LOC(NJL_FIBR,0,1).EQ.3.AND.TYPE(1:6).EQ.
     '              'SHEETS') THEN
                  NJTT=NJT+3
                  WRITE(OP_STRING,'(''>> 1 Sheet angle being '
     '              //'written'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ELSE IF(NJ_LOC(NJL_FIEL,0,nr).GT.0.AND.TYPE(1:5).EQ.
     '              'FIELD') THEN

C LKC 4-JUL-1999 generialise for regions
C                  NJTT=NJT+NJ_LOC(NJL_FIEL,0,1)

C MHT 23.03.11. using NT_FIELD to control number of fields output
C                  NJTT=NJT+NJ_LOC(NJL_FIEL,0,nr)
                  NJTT = NJT+NT_FIELD
                  WRITE(OP_STRING,
     '              '(''>> '',I1,'' Field value being written'')')
C     '              NJ_LOC(NJL_FIEL,0,nr)
     '              NT_FIELD
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ELSE IF(CALL_DATA_FIELD.AND.TYPE(1:5).EQ.'FIELD') THEN
                  NJTT=NJT+1
                  WRITE(OP_STRING,'(''>> 1 Field data '
     '              //'being written'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C*** 10/10/08 JHC Added contact pressure option
                ELSE IF(CALL_CONT.AND.TYPE(1:7).EQ.'CONTACT') THEN
                  NJTT=NJT
                  WRITE(OP_STRING,'(''>> Contact data being '
     '              //'written'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ELSE
                  NJTT=NJT
                  WRITE(OP_STRING,'(''>> Geometric data being '
     '              //'written'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDIF
C*** 10/10/08 JHC Pass in Z_CONT to IODATA.f
C                CALL IODATA('WRITE',ANGLE_TYPE,TYPE,IFILE,NDP,NJTT,
C     '            TITLE,WD,ZD,ERROR,*9999)
                CALL IODATA('WRITE',ANGLE_TYPE,TYPE,IFILE,NDP,NJTT,
     '            TITLE,WD,Z_CONT,ZD,ERROR,*9999)
              ENDIF

            ENDIF
            IF(IOTYPE.GT.1) THEN
              CALL CLOSEF(IFILE,ERROR,*9999)
            ENDIF
          ELSE IF(TYPE(1:9).EQ.'ELECTRODE') THEN
            IF(FILEFORM(1:10).EQ.'ACTIVATION') THEN
C  Reads activation time along with electrode position from .IPDATA
              CALL STRING_TRIM(FILE,IBEG,IEND)
              CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.ipdata',
     '          STATUS,'SEQUEN','FORMATTED',132,ERROR,*9999)
              IF(IOTYPE.EQ.2) THEN ! read data file
                WRITE(OP_STRING,'('' >> Electrode geometric data '
     '            //'and activation times being read'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C*** 10/10/08 JHC Pass in Z_CONT to IODATA.f
C                CALL IODATA('READ',ANGLE_TYPE,TYPE,IFILE,NDP,NJT+1,
C     '            TITLE,WD,ZD,ERROR,*9999)
                CALL IODATA('READ',ANGLE_TYPE,TYPE,IFILE,NDP,NJT+1,
     '            TITLE,WD,Z_CONT,ZD,ERROR,*9999)
                OP_STRING(1)=TITLE
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ELSE IF(IOTYPE.EQ.3) THEN ! write data file
                WRITE(OP_STRING,'('' >> Electrode geometric data '
     '            //'and activation times being written'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C*** 10/10/08 JHC Pass in Z_CONT to IODATA.f
C                CALL IODATA('WRITE',ANGLE_TYPE,TYPE,IFILE,NDP,NJT+1,
C     '            TITLE,WD,ZD,ERROR,*9999)
                CALL IODATA('WRITE',ANGLE_TYPE,TYPE,IFILE,NDP,NJT+1,
     '            TITLE,WD,Z_CONT,ZD,ERROR,*9999)
              ENDIF
              CALL CLOSEF(IFILE,ERROR,*9999)

            ELSE IF(FILEFORM(1:4).EQ.'EMAP') THEN
              ERROR='>> Not yet implemented'
              GOTO 9999

            ELSE IF(FILEFORM(1:6).EQ.'IPFILE') THEN
              CALL STRING_TRIM(FILE,IBEG,IEND)
              CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.ipdata',
     '          STATUS,'SEQUEN','FORMATTED',132,ERROR,*9999)
              IF(IOTYPE.EQ.2) THEN ! read data file
                WRITE(OP_STRING,'('' >> Electrode geometric data '
     '            //'being read'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C*** 10/10/08 JHC Pass in Z_CONT to IODATA.f
C                CALL IODATA('READ',ANGLE_TYPE,TYPE,IFILE,NDP,NJT,
C     '            TITLE,WD,ZD,ERROR,*9999)
                CALL IODATA('READ',ANGLE_TYPE,TYPE,IFILE,NDP,NJT,
     '            TITLE,WD,Z_CONT,ZD,ERROR,*9999)
                OP_STRING(1)=TITLE
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ELSE IF(IOTYPE.EQ.3) THEN ! write data file
                WRITE(OP_STRING,'('' >> Electrode geometric data '
     '            //'being written'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C*** 10/10/08 JHC Pass in Z_CONT to IODATA.f
C                CALL IODATA('WRITE',ANGLE_TYPE,TYPE,IFILE,NDP,NJT,
C     '            TITLE,WD,ZD,ERROR,*9999)
                CALL IODATA('WRITE',ANGLE_TYPE,TYPE,IFILE,NDP,NJT,
     '            TITLE,WD,Z_CONT,ZD,ERROR,*9999)
              ENDIF
              CALL CLOSEF(IFILE,ERROR,*9999)

            ELSE IF(FILEFORM(1:5).EQ.'MAP3D') THEN
              IFILE=IOFILE2
              CALL STRING_TRIM(FILE,IBEG,IEND)
              IF(IOTYPE.EQ.2) THEN ! read data file
                CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.pts',
     '            STATUS,'SEQUEN','FORMATTED',132,ERROR,*9999)
                nd=0
 10             nd=nd+1
                IF(nd.GT.NDM) THEN
                  ERROR='>># Data points > NDM'
                  GOTO 9999
                ENDIF
                READ(UNIT=IFILE,FMT=*,IOSTAT=IOSTAT,END=20)
     '            (ZD(nj,nd),nj=1,NJT)
                NDP(nd)=nd
                DO nj=1,NJT
                  WD(nj,nd)=1.0D0
                ENDDO
                GOTO 10
 20             CONTINUE
                NDT=nd-1
                CALL CLOSEF(IFILE,ERROR,*9999)
              ELSE IF(IOTYPE.EQ.3) THEN ! write data file
                CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.pts',
     '            'NEW','SEQUEN','FORMATTED',132,ERROR,*9999)
                DO nd=1,NDT
                  WRITE(IFILE,*) (ZD(nj,nd),nj=1,NJT)
                ENDDO
                CALL CLOSEF(IFILE,ERROR,*9999)
              ENDIF
            ENDIF

!news     AAY 11 Dec 94 read in image plane numbers for data points
          ELSE IF(TYPE(1:6).EQ.'PLANES')THEN
            CALL ASSERT(NDT.LE.10000,'>>Error: increase array '
     '        //'dimension in NLS_NDATA_CONT',ERROR,*9999)
            CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.IPLADA',STATUS,
     '        'SEQUEN','FORMATTED',132,ERROR,*9999)
            DO nd=nd0,nd1
              READ(IFILE,*) N,NLS_NDATA_CONT(nd),NLS_NDATA_IMAG(nd)
            ENDDO
            CALL CLOSEF(IFILE,ERROR,*9999)
!newe

C cpb 8/2/97 Adding data from a signal file.
          ELSE IF(TYPE(1:6).EQ.'SIGNAL') THEN

C***        Open the signal file
            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,SIGNAL_TIME,WD,XID,ZD,'READ',FILEFORMAT,
     '        FILE(IBEG:IEND),'OPEN',ENDFILE,.TRUE.,ERROR,*9999)

C***        Read electrode geometry etc.
            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,SIGNAL_TIME,WD,XID,ZD,'READ',FILEFORMAT,
     '        FILE(IBEG:IEND),'ELECTRODE_DATA',ENDFILE,.TRUE.,
     '        ERROR,*9999)
            CALC_XI=SIGNAL_ELEMLOC(IOFILE1).EQ.1

C cpb 8/2/97 Don't implement values at the moment as IOSIGN needs to
C be modified to return the signal data at a given time.
            IF(SIGNAL_VALUES) THEN
              ERROR='>>Not implemented'
              GOTO 9999
            ENDIF

C***        Close the signal file
            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,SIGNAL_TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,
     '        FILE(IBEG:IEND),' ',ENDFILE,.TRUE.,ERROR,*9999)

          ENDIF

        ELSE IF(MOUSE) THEN

C LKC 4-JUL-1999 Generalisation for different regions
C          nr=1 !AJP 16/1/96
          nd=NDT
          NDTOLD=NDT
          iw=IWK(1)
          CALL ACWK(iw,0,ERROR,*9999)

C cpb 28/7/94 Adding data segments
          CLABEL='1' !Web error fix
          CALL OPEN_SEGMENT(ISDATA(iw,nr),ISEG,iw,CLABEL,INDEX,
     '      INDEX_OLD,nr,1,CSEG,ERROR,*9999)

          WRITE(OP_STRING,'('' >>Locate data point on '',I1)') iw
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CALL LOCATOR(INSTAT,0.0D0,XWC1,0.0D0,YWC1,
     '      ERROR,*9999)
          DO WHILE(INSTAT.EQ.1)
            nd=nd+1
            CALL ASSERT(nd.LE.NDM,'>>NDM too small',ERROR,*9999)
            NDP(nd)=nd
            IF(iw.EQ.4.AND.MAP_PROJEC(1:2).EQ.'XI') THEN
              DO nr=1,NRT
                DO noelem=1,NEELEM(0,nr)
                  NEE=NEELEM(noelem,nr)
                  XID1=DBLE(MAX_XI)*(XWC1+1.0d0)/2.0d0-
     '              DBLE(MXI(1,NEE)-1)
                  XID2=DBLE(MAX_XI)*(YWC1+1.0d0)/2.0d0-
     '              DBLE(MXI(2,NEE)-1)
                  IF(XID1.GE.0.0d0.AND.XID1.LT.1.0d0.AND.
     '              XID2.GE.0.0d0.AND.XID2.LT.1.0d0) THEN
                    XID(1,nd)=XID1
                    XID(2,nd)=XID2
                    LD(nd)=NEE
                    GO TO 400
                  ENDIF
                ENDDO
              ENDDO
 400          ne=NEE !is element no for current data point
              CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '          NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '          SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
              DO nj=1,NJT
                nb=NBJ(nj,ne)
                ZD(nj,nd)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '            nb,1,XID(1,nd),XE(1,nj))
              ENDDO
              MXI1=MXI(1,ne)
              MXI2=MXI(2,ne)
              CALL POLYMARKER(1,iw,1,XID(1,nd),ERROR,*9999)

            ELSE
              IF(NJT.EQ.2) THEN
                ZD(1,nd)=XWC1
                ZD(2,nd)=YWC1
              ELSE IF(NJT.EQ.3) THEN
                IF(iw.EQ.1) THEN
                  ZD(1,nd)=XWC1
                  ZD(3,nd)=YWC1
                ELSE IF(iw.EQ.4) THEN
                  CALL PSCOORD1(XD,XWC1,YWC1,ERROR,*9999)
                  WRITE(OP_STRING,'('' Prolate coordinates: '','
     '              //'3E11.3)')
     '              (XD(nj),nj=1,NJT)
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  CALL XZ(ITYP10(1),XD,ZD(1,nd))
                ENDIF
              ENDIF
              CALL POLYMARKER(1,iw,1,ZD(1,nd),ERROR,*9999)
            ENDIF

            IF(njt.eq.3.and.iw.NE.4) THEN
              CALL DAWK(iw,0,ERROR,*9999)
              CALL ACWK(IWK(2),0,ERROR,*9999)
              WRITE(OP_STRING,'('' >>Locate data point on '',I1)')
     '          IWK(2)
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              CALL LOCATOR(INSTAT,0.0D0,XWC2,
     '          0.0D0,YWC2,ERROR,*9999)
              ZD(2,nd)=XWC2
              CALL POLYMARKER(1,IWK(2),1,ZD(1,nd),ERROR,*9999)
              CALL DAWK(IWK(2),0,ERROR,*9999)
              CALL ACWK(3,0,ERROR,*9999)
              CALL POLYMARKER(1,3,1,ZD(1,nd),ERROR,*9999)
              CALL DAWK(3,0,ERROR,*9999)
              CALL ACWK(iw,0,ERROR,*9999)
            ENDIF
            DO nj=1,NJT
              WD(nj,nd)=1.0D0
            ENDDO
            IF(TYPE(1:8).EQ.'GEOMETRY') THEN
c             IF(NPT.GT.0) THEN
c               IF(NJT.EQ.2) THEN
c                 WRITE(IOIP,'($,'' >>Enter closest line number: '')')
c               ELSE IF(NJT.EQ.3) THEN
c                 WRITE(IOIP,'($,'' >>Enter closest face number: '')')
c               ENDIF
c               READ(IOIP,*) LD(nd)
c             ENDIF
            ELSE IF(TYPE(1:6).EQ.'FIBRES') THEN
              WRITE(OP_STRING,'($,'' >>Enter fibre angle in '
     '          //'degrees: '')')
              CALL WRITES(IOIP,OP_STRING,ERROR,*9999)
              READ(IOIP,*) ANGLE
              ZD(NJT+1,nd)=ANGLE*PI/180.0D0
              WD(NJT+1,nd)=1.0D0
            ELSE IF(TYPE(1:5).EQ.'FIELD') THEN
              IF(NJ_LOC(NJL_FIEL,0,nr).EQ.1) THEN
                FORMAT='($,'' >>Enter field value [0.0]:'',D12.4)'
                RDEFLT(1)=0.0d0
                WRITE(OP_STRING,'($,'' >>Enter field value: '')')
                CALL WRITES(IOIP,OP_STRING,ERROR,*9999)
                READ(IOIP,*) ZD(NJT+1,nd)
                WD(NJT+1,nd)=1.0D0
              ELSE IF(NJ_LOC(NJL_FIEL,0,nr).EQ.2) THEN
                WRITE(OP_STRING,'($,'' >>Enter 2 field values: '')')
                CALL WRITES(IOIP,OP_STRING,ERROR,*9999)
                READ(IOIP,*) ZD(NJT+1,nd),
     '            ZD(NJT+2,nd)
                WD(NJT+1,nd)=1.0d0
                WD(NJT+2,nd)=1.0d0
              ENDIF
            ENDIF
            IF(BEM) THEN !Calc the solution at the current data point
              DO nj=1,NJT
                XPFP(nj)=ZD(nj,nd)
              ENDDO
              ne=1 ! This will eventually identify an element in the BE region.
              nr=1 !Needs fixing
              nx=1 !temporary
              CALL EQTYPE(IBT,NBH,NEELEM,nr,NW(1,1,nx),nx,ERROR,*9999)
C MLB please leave
C              CALL DOMSOL(NBH,NBJ,NEELEM,NHE(1,nx),
C     '          NHP(1,nr,nx),NKH(1,1,1,nr),NKHE,NKJE,
C     '          NLL,NPF,NP_INTERFACE,NPNE,NPNODE,
C     '          nr,NRE,NVHE,NVHP(1,1,1,nr),NVJE,NW(1,1,nx),nx,
C     '          NYNE,NYNP,CE(1,1,nx),CURVCORRECT,DL,PG,RG,SE,WG,
C     '          XA,XE,XG1,XP,XPFP,
C     '          YD,YP(1,1,nx),ZA,ZE,ZF,ZP,ERROR,*9999)
              CALL DOMSOL2(NBH,NBJ,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '          NKH(1,1,1,nr),NKHE,NKJE,NLL,NPF,NP_INTERFACE,NPNE,
     '          NPNODE,nr,NRE,NVHE,NVHP(1,1,1,nr),NVJE,NW(1,1,nx),nx,
     '          NYNE,NYNP,CE(1,1,nx),CURVCORRECT,DET,DL,DRDN,PG,RAD,
     '          RD,RG,SE,WG,XA,XE,XG1,XN,XP,XPFP,XR,YD,YP(1,1,nx),ZA,
     '          ZE,ZF,ZP,ERROR,*9999)
              IF(.NOT.COMPLEX) THEN
                ZD(NJT+1,nd)=YD(1)
              ELSE !Complex solution.  Store magnitude of the vector
              ENDIF
              WRITE(OP_STRING,'(3(E10.4,4X))')
     '          ZD(1,nd),ZD(2,nd),ZD(3,nd)
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            WRITE(OP_STRING,'('' >>Locate data point on '',I1)') iw
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            CALL LOCATOR(INSTAT,0.0D0,XWC1,0.0D0,
     '        YWC1,ERROR,*9999)
          ENDDO

C cpb 28/7/94 Adding data segments
          CALL CLOSE_SEGMENT(ISDATA(iw,nr),iw,ERROR,*9999)

          CALL DAWK(iw,0,ERROR,*9999)
          NDT=nd

          IF(iw.EQ.4.AND.MAP_PROJEC(1:2).EQ.'XI') THEN
            CALL DSTATS(LD,NDDL,NDLT,ERROR,*9999)
            NXIDEF=2
          ENDIF

          CALL STRING_TRIM(FILE,IBEG,IEND)
          CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.ipdata',STATUS,
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
          IF(TYPE(1:8).EQ.'GEOMETRY') THEN
C*** 10/10/08 JHC Pass in Z_CONT to IODATA.f
C            CALL IODATA('WRITE',ANGLE_TYPE,TYPE,IFILE,NDP,NJT,TITLE,
C     '        WD,ZD,ERROR,*9999)
            CALL IODATA('WRITE',ANGLE_TYPE,TYPE,IFILE,NDP,NJT,TITLE,
     '        WD,Z_CONT,ZD,ERROR,*9999)
          ELSE IF(TYPE(1:6).EQ.'FIBRES') THEN
C*** 10/10/08 JHC Pass in Z_CONT to IODATA.f
C            CALL IODATA('WRITE',ANGLE_TYPE,TYPE,IFILE,NDP,NJT+1,TITLE,
C     '        WD,ZD,ERROR,*9999)
            CALL IODATA('WRITE',ANGLE_TYPE,TYPE,IFILE,NDP,NJT+1,TITLE,
     '        WD,Z_CONT,ZD,ERROR,*9999)
          ELSE IF(TYPE(1:5).EQ.'FIELD') THEN
C*** 10/10/08 JHC Pass in Z_CONT to IODATA.f
C            CALL IODATA('WRITE',ANGLE_TYPE,TYPE,IFILE,NDP,NJT+1,TITLE,
C     '        WD,ZD,ERROR,*9999)
            CALL IODATA('WRITE',ANGLE_TYPE,TYPE,IFILE,NDP,NJT+1,TITLE,
     '        WD,Z_CONT,ZD,ERROR,*9999)
          ENDIF
          CALL CLOSEF(IFILE,ERROR,*9999)

        ELSE IF(CALCU) THEN !Calculate data position info

          IF(TYPE(1:6).EQ.'STRIPE') THEN
!old MPN 1-Nov-94: CALC_STRIPE_PTS *** ARCHIVED ***
c            CALL CALC_STRIPE_PTS(NT_STRIPE_DATA,NDP,WD,ZD,GROUP_NAME,
c     '        ERROR,*9999)
          ELSE IF(TYPE(1:5).NE.'SHEET') THEN
            NXIDEF=NIT(NBJ(1,NEELEM(1,1)))
            L=0
            DO noelem=1,NEELEM(0,1) !to initialize LN (recalc.d by define fit)
              ne=NEELEM(noelem,1)
              L=L+1
              LN(L)=ne
            ENDDO
            LN(0)=L
          ENDIF

          IF(TYPE(1:5).EQ.'SHEET') THEN
            DO nd=1,NDT
              ne=LD(nd)
              IF(DOP) THEN
                WRITE(OP_STRING,'(/'' nd='',I6,'' element '',I4)')
     '            nd,ne
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              IF(ne.GT.0) THEN
                XI(1)=XID(1,nd)
                XI(2)=XID(2,nd)
                XI(3)=XID(3,nd)
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' Xi coords: '',2E12.3)')
     '              XI(1),XI(2)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '            NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '            SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)

C new CS 16/1/2002 updating
C !!! Code nees updating for new fibre/sheet angle transformations
C !!! (see subr mat_vec_ng etc)
C                CALL ASSERT(.FALSE.,'>>ERROR: Code is out of date',
C     '            ERROR,*9999)

!Lamda
C                nb=NBJ(1,ne)
C                RLAMDA     = PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
C     '            nb,1,XI,XE(1,1))
C                IF(DOP) THEN
C                  WRITE(OP_STRING,'('' RLAMDA ='',E12.3)') RLAMDA
C                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                ENDIF
C                dLAMDA_dXi1= PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
C     '            nb,2,XI,XE(1,1))
C                dLAMDA_dXi2= PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
C     '            nb,4,XI,XE(1,1))
C                sinh_LAMDA = DSINH(RLAMDA)
C                cosh_LAMDA = DCOSH(RLAMDA)
C                IF(DABS(cosh_LAMDA).GT.1.0D-6) THEN
C                  tanh_LAMDA = sinh_LAMDA/cosh_LAMDA
C                ELSE
C                  tanh_LAMDA = 1.0D8
C                ENDIF
C                IF(DABS(sinh_LAMDA).GT.1.0D-6) THEN
C                  coth_LAMDA = cosh_LAMDA/sinh_LAMDA
C                ELSE
C                  coth_LAMDA = 1.0D8
C                ENDIF
!Mu
C                nb=NBJ(2,ne)
C                RMU        = PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
C     '            nb,1,XI,XE(1,2))
C                IF(DOP) THEN
C                  WRITE(OP_STRING,'(''    RMU ='',E12.3)') RMU
C                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                ENDIF
C                dMU_dXi1   = PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
C     '            nb,2,XI,XE(1,2))
C                dMU_dXi2   = PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
C     '            nb,4,XI,XE(1,2))
C                sin_MU     = DSIN(RMU)
C                cos_MU     = DCOS(RMU)
C                IF(DABS(cos_MU).GT.1.0D-6) THEN
C                  tan_MU = sin_MU/cos_MU
C                ELSE
C                  tan_MU = 1.0D8
C                ENDIF
C                IF(DABS(sin_MU).GT.1.0D-6) THEN
C                  cot_MU     = cos_MU/sin_MU
C                ELSE
C                  cot_MU = 1.0D8
C                ENDIF
!Theta
C                nb=NBJ(3,ne)
C                THETA      = PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
C     '            nb,1,XI,XE(1,3))
C                IF(DOP) THEN
C                  WRITE(OP_STRING,'(''  THETA ='',E12.3)') THETA
C                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                ENDIF
C                dTHETA_dXi1= PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
C     '            nb,2,XI,XE(1,3))
C                dTHETA_dXi2= PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
C     '            nb,4,XI,XE(1,3))
C                sin_THETA  = DSIN(THETA)
C                cos_THETA  = DCOS(THETA)
!X
                nb=NBJ(1,ne)
                X = PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '            nb,1,XI,XE(1,1))
                dX_dXi1= PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '            nb,2,XI,XE(1,1))
                dX_dXi2= PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '            nb,4,XI,XE(1,1))
!Y
                nb=NBJ(2,ne)
                Y = PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '            nb,1,XI,XE(1,2))
                dY_dXi1= PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '            nb,2,XI,XE(1,2))
                dY_dXi2= PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '            nb,4,XI,XE(1,2))
!Z
                nb=NBJ(3,ne)
C                Z = PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
C     '            nb,1,XI,XE(1,3))
                dZ_dXi1= PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '            nb,2,XI,XE(1,3))
                dZ_dXi2= PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '            nb,4,XI,XE(1,3))
!Theta
                 THETA = DATAN(Y/X)
                 IF((X.LT.0).AND.(Y.GT.0)) THETA=PI+THETA
                 IF((X.LT.0).AND.(Y.LT.0)) THETA=PI+THETA
                 IF((X.GT.0).AND.(Y.LT.0)) THETA=2.0d0*PI+THETA
!Alfa (fibre angle)
                nb=NBJ(4,ne)
                ALFA       = PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '            nb,1,XI,XE(1,4))
                IF(DOP) THEN
                  WRITE(OP_STRING,'(''   ALFA ='',E12.3)') ALFA
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
C                sin_ALFA   = DSIN(ALFA)
C                cos_ALFA   = DCOS(ALFA)
!Beta (sheet angle)
                BETA_dash  =  PI/2.0d0-ZD(NJT+3,nd)
                IF(DOP) THEN
                  WRITE(OP_STRING,
     '              '(''   BETA-dash ='',E12.3)') BETA_dash
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                sin_BETA   = DSIN(BETA_dash)
                cos_BETA   = DCOS(BETA_dash)
!Gama
C                ARG1=coth_LAMDA*dLAMDA_dXi1+cot_MU*dMU_dXi1
C                IF(DABS(ARG1).GT.1.0D-6) THEN
C                  GAMA =DATAN2(-ARG1,dTHETA_dXi1)
C                ELSE
C                  GAMA=0.0D0
C                ENDIF
C                IF(GAMA.GT.PI/2.0D0) THEN
C                  GAMA=GAMA-PI
C                ELSE IF(GAMA.LT.-PI/2.0D0) THEN
C                  GAMA=GAMA+PI
C                ENDIF
C                IF(DOP) THEN
C                  WRITE(OP_STRING,'(''   GAMA ='',E12.3)') GAMA
C                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                ENDIF
C                sin_GAMA   = DSIN(GAMA)
C                cos_GAMA   = DCOS(GAMA)
!Delta
C                ARG1=(tan_MU*dLAMDA_dXi2+tanh_LAMDA*dMU_dXi2)*cos_GAMA
C                IF(DABS(ARG1).GT.1.0D-6) THEN
C                  DELTA=DATAN2(ARG1,(tan_MU*dMU_dXi2-tanh_LAMDA
C     '              *dLAMDA_dXi2))
C                ELSE
C                  DELTA=0.0D0
C                ENDIF
C                IF(DELTA.GT.PI/2.0D0) THEN
C                  DELTA=DELTA-PI
C                ELSE IF(DELTA.LT.-PI/2.0D0) THEN
C                  DELTA=DELTA+PI
C                ENDIF
C                IF(DOP) THEN
C                  WRITE(OP_STRING,'(''  DELTA ='',E12.3)') DELTA
C                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                ENDIF
C                sin_DELTA  = DSIN(DELTA)
C                cos_DELTA  = DCOS(DELTA)
!Alfa vector
C                ALFA_VECTOR(1)=-sin_ALFA*cos_DELTA
C                ALFA_VECTOR(2)= cos_ALFA*sin_GAMA+sin_ALFA*cos_GAMA
C     '            *sin_DELTA
C                ALFA_VECTOR(3)=-cos_ALFA*cos_GAMA+sin_ALFA*sin_GAMA
C     '            *sin_DELTA
C                CALL ROTATION(1,THETA,ALFA_VECTOR,ERROR,*9999)
C                IF(DOP) THEN
C                  WRITE(OP_STRING,'('' ALFA_vector:'',3E12.3)')
C     '              (ALFA_VECTOR(i),i=1,3)
C                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                ENDIF
!Alfa vector
                  CALL MAT_VEC_XI(IBT,IDO,INP,NAN,NBJ(1,ne),nr,
     '              ALFA_VECTOR,B,C,XE,XG,Xi,.TRUE.,ERROR,*9999)
!Beta-dash vector
C                BETA_VECTOR(1)=-sin_BETA
C                BETA_VECTOR(2)= cos_BETA
C                BETA_VECTOR(3)= 0.0D0
C                CALL ROTATION(1,THETA,BETA_VECTOR,ERROR,*9999)
C                IF(DOP) THEN
C                  WRITE(OP_STRING,'('' BETA-dash_vector:'',3E12.3)')
C     '              (BETA_VECTOR(i),i=1,3)
C                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                ENDIF
!Beta-dash vector
                BETA_VECTOR(1)=-sin_BETA
                BETA_VECTOR(2)= cos_BETA
                BETA_VECTOR(3)= THETA
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' BETA-dash_vector:'',3E12.3)')
     '              (BETA_VECTOR(i),i=1,3)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                ! convert to rc
                tmp=BETA_VECTOR(1)
                BETA_VECTOR(1)=BETA_VECTOR(2)*cos(THETA)
                BETA_VECTOR(2)=BETA_VECTOR(2)*sin(THETA)
                BETA_VECTOR(3)=tmp
!Normal vector
                CALL CROSS(ALFA_VECTOR,BETA_VECTOR,RN_VECTOR)
                ABS_RN = DSQRT(RN_VECTOR(1)**2+RN_VECTOR(2)**2
     '            +RN_VECTOR(3)**2)
                RN_VECTOR(1) = RN_VECTOR(1)/ABS_RN
                RN_VECTOR(2) = RN_VECTOR(2)/ABS_RN
                RN_VECTOR(3) = RN_VECTOR(3)/ABS_RN
                IF(DOP) THEN
                  WRITE(OP_STRING,'(''   RN_vector:'',3E12.3)')
     '              (RN_VECTOR(i),i=1,3)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
!Beta vector
                CALL CROSS(RN_VECTOR,ALFA_VECTOR,BETA_VECTOR)
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' BETA_vector:'',3E12.3)')
     '              (BETA_VECTOR(i),i=1,3)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
!Xi(1) base vector
C                G1_VECTOR(1) = sinh_LAMDA*dLAMDA_dXi1*cos_MU
C     '            -cosh_LAMDA*sin_MU*dMU_dXi1
C                G1_VECTOR(2) = cosh_LAMDA*dLAMDA_dXi1*sin_MU
C     '            +sinh_LAMDA*cos_MU*dMU_dXi1
C                G1_VECTOR(3) = sinh_LAMDA*sin_MU*dTHETA_dXi1
C                ABS_G1 = DSQRT(G1_VECTOR(1)**2+G1_VECTOR(2)**2
C     '            +G1_VECTOR(3)**2)
C                G1_VECTOR(1) = G1_VECTOR(1)/ABS_G1
C                G1_VECTOR(2) = G1_VECTOR(2)/ABS_G1
C                G1_VECTOR(3) = G1_VECTOR(3)/ABS_G1
C                CALL ROTATION(1,THETA,G1_VECTOR,ERROR,*9999)
C                IF(DOP) THEN
C                  WRITE(OP_STRING,'(''   G1_vector:'',3E12.3)')
C     '              (G1_VECTOR(i),i=1,3)
C                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                ENDIF
!Xi(2) base vector
C                G2_VECTOR(1) = sinh_LAMDA*dLAMDA_dXi2*cos_MU
C     '            -cosh_LAMDA*sin_MU*dMU_dXi2
C                G2_VECTOR(2) = cosh_LAMDA*dLAMDA_dXi2*sin_MU
C     '            +sinh_LAMDA*cos_MU*dMU_dXi2
C                G2_VECTOR(3) = sinh_LAMDA*sin_MU*dTHETA_dXi2
C                ABS_G2 = DSQRT(G2_VECTOR(1)**2+G2_VECTOR(2)**2
C     '            +G2_VECTOR(3)**2)
C                G2_VECTOR(1) = G2_VECTOR(1)/ABS_G2
C                G2_VECTOR(2) = G2_VECTOR(2)/ABS_G2
C                G2_VECTOR(3) = G2_VECTOR(3)/ABS_G2
C                CALL ROTATION(1,THETA,G2_VECTOR,ERROR,*9999)
C                IF(DOP) THEN
C                  WRITE(OP_STRING,'(''   G2_vector:'',3E12.3)')
C     '              (G2_VECTOR(i),i=1,3)
C                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                ENDIF
!Xi(1) base vector
                G1_VECTOR(1)=dX_dXi1
                G1_VECTOR(2)=dY_dXi1
                G1_VECTOR(3)=dZ_dXi1
!Xi(2) base vector
                G2_VECTOR(1)=dX_dXi2
                G2_VECTOR(2)=dY_dXi2
                G2_VECTOR(3)=dZ_dXi2
!u vector
                CALL CROSS(G1_VECTOR,G2_VECTOR,U_VECTOR)
                ABS_U = DSQRT(U_VECTOR(1)**2+U_VECTOR(2)**2+U_VECTOR(3)
     '            **2)
                U_VECTOR(1)  = U_VECTOR(1)/ABS_U
                U_VECTOR(2)  = U_VECTOR(2)/ABS_U
                U_VECTOR(3)  = U_VECTOR(3)/ABS_U
                IF(DOP) THEN
                  WRITE(OP_STRING,'(''    U_vector:'',3E12.3)')
     '              (U_VECTOR(i),i=1,3)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
!p vector
                CALL CROSS(U_VECTOR,ALFA_VECTOR,P_VECTOR)
                ABS_P = DSQRT(P_VECTOR(1)**2+P_VECTOR(2)**2+
     '            P_VECTOR(3)**2)
                P_VECTOR(1)  = P_VECTOR(1)/ABS_P
                P_VECTOR(2)  = P_VECTOR(2)/ABS_P
                P_VECTOR(3)  = P_VECTOR(3)/ABS_P
                IF(DOP) THEN
                  WRITE(OP_STRING,'(''    P_vector:'',3E12.3)')
     '              (P_VECTOR(i),i=1,3)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
!Beta (corrected sheet angle)
C                BETA_pu(1)=
C     '            BETA_VECTOR(1)*(-cos_ALFA         *cos_DELTA)
C     '            +BETA_VECTOR(2)*(-sin_ALFA*sin_GAMA*cos_THETA
C     '            -sin_ALFA*cos_GAMA          *sin_THETA
C     '            +cos_ALFA*cos_GAMA*sin_DELTA*cos_THETA
C     '            -cos_ALFA*sin_GAMA*sin_DELTA*sin_THETA)
C     '            +BETA_VECTOR(3)*(-sin_ALFA*sin_GAMA*sin_THETA
C     '            +sin_ALFA*cos_GAMA          *cos_THETA
C     '            +cos_ALFA*cos_GAMA*sin_DELTA*sin_THETA
C     '            +cos_ALFA*sin_GAMA*sin_DELTA*cos_THETA)
C                BETA_pu(2)=
C     '            BETA_VECTOR(1)*                 sin_DELTA
C     '            +BETA_VECTOR(2)*( cos_GAMA*cos_DELTA*cos_THETA
C     '            -sin_GAMA         *cos_DELTA*sin_THETA)
C     '            +BETA_VECTOR(3)*( cos_GAMA*cos_DELTA*sin_THETA
C     '            +sin_GAMA         *cos_DELTA*cos_THETA)
C                IF(DOP) THEN
C                  WRITE(OP_STRING,'('' BETA_pu_1='',E12.3,'
C     '              //''' BETA_pu_2='',E12.3)') BETA_pu(1),BETA_pu(2)
C                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                ENDIF
!Beta (corrected sheet angle)
                BETA_pu(1)=DDOT(3,BETA_VECTOR,1,U_VECTOR,1)
                BETA_pu(2)=DDOT(3,BETA_VECTOR,1,P_VECTOR,1)
                BETA=PI/2.0d0-DATAN2(BETA_pu(2),BETA_pu(1))
                IF(BETA.LT.0.0d0) THEN
                  BETA=BETA+PI
                ELSE IF(BETA.GT.PI) THEN
                  BETA=BETA-PI
                ENDIF
                ZD(NJT+3,nd)=BETA
                IF(DOP) THEN
                  WRITE(OP_STRING,'(''   BETA ='',E12.3)') BETA
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF

              ENDIF
            ENDDO
            CALC_SHEET=.TRUE.

          ELSE IF(TYPE(1:4).EQ.'GRID') THEN
C****  Calculate nearest grid points to each electrode
            CALL ASSERT(NQT.GT.0,'>>Grid points not defined',
     '        ERROR,*9999)
            CALL ASSERT(NDT.GT.0,'>>Electrode positions not defined',
     '        ERROR, *9999)
            DO nd=1,NDT
              MINDIST=1.d5
              DO nq=1,NQT
                DIST=0.0d0
                DO nj=1,NJT
                  DIST=DIST+(ZD(nj,nd)-XQ(nj,nq))*(ZD(nj,nd)-XQ(nj,nq))
                ENDDO
                IF(DIST.LT.MINDIST) THEN
                  LDR(nd)=nq
                  MINDIST=DIST
                ENDIF
              ENDDO
              IF(dop) THEN
                WRITE(OP_STRING,'('' nd ='',I5,''('',I5,'') nq ='',I7,'
     '            //''' r ='',F10.2)') nd,NDP(nd),LDR(nd),sqrt(MINDIST)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO

          ELSE IF(TYPE(1:9).EQ.'FROM_GRID') THEN
C****  Create a data point for each grid point
            CALL ASSERT(NQT.GT.0,'>>Grid points not defined',
     '        ERROR,*9999)
            CALL ASSERT(NDM.GE.NQT,'>>Increase NDM',
     '        ERROR, *9999)
C GDR 8mar05 Added the following loop to copy the grid coords to the
C            data coords. Previously the "from_grid" command did not
C            work. 
            NDT=NQT
            DO nd=1,NDT
              DO nj=1,NJT
                ZD(nj,nd)=XQ(nj,nd)
                WD(nj,nd)=1.0D0
                NDP(nd)=nd
              ENDDO
            ENDDO
            IF(CBBREV(CO,'NUM_FIELDS',3,noco+1,NTCO,N3CO)) THEN
              NT_FIELD=IFROMC(CO(N3CO+1))
              IF(NJ_LOC(NJL_FIEL,0,nr).LE.NT_FIELD) THEN
                NJTtemp = NJ_LOC(NJL_GEOM,0,nr)+
     '            NJ_LOC(NJL_FIBR,0,nr)+NT_FIELD
                WRITE(ASSERT_STRING,'(''>> Increase NJM '
     '            //'to '',I1)') NJTtemp
                CALL ASSERT(NJTtemp.LE.NJM,
     '            ASSERT_STRING,ERROR,*9999)
                NJ_LOC(NJL_FIEL,0,nr)=NT_FIELD
                DO nj=1,NJ_LOC(NJL_FIEL,0,nr)
C                 NJ_LOC(NJL_FIEL,nj,nr)=nj+NJ_LOC(0,0,nr)
                  NJ_LOC(NJL_FIEL,nj,nr)=nj+
     '              NJ_LOC(NJL_GEOM,0,nr)+NJ_LOC(NJL_FIBR,0,nr)
                ENDDO
C LKC 17-JAN-2003 Touching up the way the number fields are calculated
C               when data points with fields are reread in.
C               NJ_LOC(0,0,nr)=NJ_LOC(0,0,nr)+NT_FIELD
C               NJ_LOC(0,0,0)=NJ_LOC(0,0,0)+NT_FIELD
                NJ_LOC(0,0,nr)=NJ_LOC(1,0,nr)+NJ_LOC(2,0,nr)
     '            +NT_FIELD
                NJ_LOC(0,0,0)=NJ_LOC(1,0,nr)+NJ_LOC(2,0,nr)
     '            +NT_FIELD
                NJ_LOC(NJL_FIEL,0,0)=NJ_LOC(NJL_FIEL,0,0)+NT_FIELD
                NJ_LOC(NJL_FIEL,0,nr)=NT_FIELD
              ENDIF
            ENDIF
              

          ELSE IF(TYPE(1:6).EQ.'CATDAT') THEN
            IF(NET(1).GT.0.AND.NIT(NBJ(1,1)).EQ.1) THEN
              NBRANCH=1
              NBD(1,NBRANCH)=0
              NBD(2,NBRANCH)=0
              NBD(3,NBRANCH)=0
              nd=0
              DO nolist=1,NELIST(0)
                ne=NELIST(nolist)
                nb=NBJ(1,ne)
                np=NPNE(1,nb,ne)
                IF(NKJ(1,np).GT.2) THEN
C ***             this element is the start of a new branch
                  IF(DOP) THEN
                    WRITE(OP_STRING,*)
     '                ' New branch after',nd,' data points'
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                  NBRANCH=NBRANCH+1
                  NBD(1,NBRANCH)=0
                  NBD(2,NBRANCH)=0
                  NBD(3,NBRANCH)=0
C ***           search for parent branch & add to its list of children
                  DO nnd=1,nd
                    IF(DABS(ZD(1,nnd)-XP(1,nv,1,np)).LT.1.0D0.AND.
     '                DABS(ZD(2,nnd)-XP(1,nv,2,np)).LT.1.0D0.AND.
     '                DABS(ZD(3,nnd)-XP(1,nv,3,np)).LT.1.0D0) THEN
                      IF(DOP) THEN
                        WRITE(OP_STRING,*)
     '                    ' Parent of',NBRANCH,' is',LD(nnd)
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      ENDIF
                      IF(NBD(2,LD(nnd)).EQ.0) THEN
                        NBD(2,LD(nnd))=NBRANCH
                      ELSE IF(NBD(3,LD(nnd)).EQ.0) THEN
                        NBD(3,LD(nnd))=NBRANCH
                      ELSE
                        ERROR='Run out of children in NBD'
                        GOTO 9999
                      ENDIF
                    ENDIF
                  ENDDO
                ENDIF
                CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '            NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '            SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
                nl=NLL(1,ne)
                IF(nl.NE.0) THEN
                  NSTEP=NINT(DL(3,nl)/STEP)
                  NBD(1,NBRANCH)=NBD(1,NBRANCH)+NSTEP
                  DSTEP=DL(3,nl)/DBLE(NSTEP)
                  IF(DOP) THEN
                    WRITE(OP_STRING,*)' ne,nl,nbranch,nstep,dstep= ',
     '                ne,nl,NBRANCH,NSTEP,DSTEP
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                  XXI=0.0D0
                  DO i=1,NSTEP
                    S=DSTEP*i
                    nd=nd+1
                    IF(nd.GT.NDM) THEN
                      ERROR='NDM is too small'
                      GOTO 9999
                    ENDIF
                    LD(nd)=NBRANCH
C ***               calculate xi for this arclength - Newton iteration
                    DXI=1.0D0
                    DO WHILE (DABS(DXI).GT.0.00001D0)
C ***                 calculate arclength sxi at xi=xxi - Gaussian quad
                      DO ng=1,3
                        GPD(ng)=GP(ng)*XXI
                      ENDDO
                      SXI=0.0D0
                      DO ng=1,3
                        FUNC=0.0D0
                        DO nj=1,3
                          FUNC=FUNC+PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                      INP(1,1,nb),nb,2,GPD(ng),XE(1,nj))**2
                        ENDDO
                        SXI=SXI+DSQRT(FUNC)*GW(ng)*XXI
                      ENDDO
C ***                 now get derivative of s wrt xi at xxi
                      DSDXI=0.0D0
                      XXI2(1)=XXI
                      DO nj=1,3
                        DSDXI=DSDXI+PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                    INP(1,1,nb),nb,2,XXI2,XE(1,nj))**2
                      ENDDO
                      DSDXI=DSQRT(DSDXI)
C ***                 perform Newton step
                      DXI=(S-SXI)/DSDXI
                      XXI=XXI+DXI
                    ENDDO
                    XID(1,nd)=XXI
                    DO nj=1,NJT
                      ZD(nj,nd)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                  INP(1,1,nb),nb,1,XID(1,nd),XE(1,nj))
                      WD(nj,nd)=1.0D0
                    ENDDO
                    IF(DOP) THEN
                      WRITE(OP_STRING,*)' xxi,zd=',
     '                  xxi,(zd(i2,nd),i2=1,3)
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO
              NDT=nd
            ELSE IF(NET(1).EQ.0) THEN
              OP_STRING(1)=' >>No elements defined'
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
          ELSEIF(TYPE(1:7).EQ.'FROM_XI') THEN
            CALL ASSERT(CALC_XI,'>>Read in XI first',ERROR,*9999)
            CALL ASSERT(NDT.GT.0,'>>NDT.LE.0',ERROR,*9999)

C LKC 25-OCT-98 Need to setup NDDATA & WD & NDP
C note: nr is not setup for data yet

            DO nd=1,NDT
C KFA 25-JAN-2001
C Added changes which allow the calculation of data points using
C deformed geometery and also checks which skip the calculations
C for non-projected data points.
C NOTE : Deformed coords are in XE.
              IF (LD(nd).gt.0) THEN
                IF(.NOT.UNDEFORMED) THEN
                  CALL ZPZE(NBH(1,1,LD(nd)),1,NHE(LD(nd),nx),
     '              NKHE(1,1,1,LD(nd)),NPF(1,1),NPNE(1,1,LD(nd)),
     '              NRE(LD(nd)),NVHE(1,1,1,LD(nd)),NW(LD(nd),1,nx),
     '              nx,CURVCORRECT(1,1,1,LD(nd)),SE(1,1,LD(nd)),
     '              ZA(1,1,1,LD(nd)),XE,ZP,ERROR,*9999)
                ELSE
                  CALL XPXE(NBJ(1,LD(nd)),NKJE(1,1,1,LD(nd)),
     '              NPF(1,1),NPNE(1,1,LD(nd)),NRE(LD(nd)),
     '              NVJE(1,1,1,LD(nd)),SE(1,1,LD(nd)),
     '              XA(1,1,LD(nd)),XE,XP,ERROR,*9999)
                ENDIF
                DO ni=1,NIT(NBJ(1,LD(nd)))
                  XI(ni)=XID(ni,nd)
                ENDDO
                DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                  nb=NBJ(nj,LD(nd))
                  ZD(nj,nd)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '              INP(1,1,nb),nb,1,XI,XE(1,nj))
                  WD(nj,nd)=1.D0
                ENDDO !nj
                IF(NDP(nd).EQ.0)THEN
                  NDP(nd)=nd
                ENDIF
              ELSE
                DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
C GRC 25-JUN-2001 Do not clear coordinates here. Makes sense to leave
C them as they are - a behaviour needed for host-mesh warping of part
C of a data set, where points outside the host mesh should stay where
C they are. Additionally, the code must not rely on this command to
C initialise data point fields!
C                  ZD(nj,nd)=0.D0
                  WD(nj,nd)=0.D0
                ENDDO
              ENDIF
C LKC done latter on, need to do it here when there are regions
!              NDDATA(nd,nr)=nd
C GRC 25-JUN-2001 Do not reset NDP here. This allows us to write data
C files that have been host-mesh warped and still relate their NDP
C numbers back to the source data. Even if a data point is not in an
C element leave NDP as it was -- we want those points to remain
C unchanged in number and position when we warp those data points that
C are in elements. Clearing NDP should never be the default behaviour!
C              NDP(nd)=nd
            ENDDO
          ELSEIF(TYPE(1:8).EQ.'FROM_TAG') THEN
C KFA 2002/11/25
C The following code tags a point and an image normal
C and calculates the resulting deformed tag points
C wrt to image tagging.

            ERROR_F=.FALSE.

C C$OMP       PARALLEL DO
C C$OMP&      PRIVATE(CC,CONVERG,DEF,DX,F,FL,J,nd,nj,OK,RTSEC,
C C$OMP&        SWAP,UNDEF,X1,X2,XL)
            DO nd=1,NDT
              IF(.NOT.ERROR_F) THEN
C                IF(MOD(nd,NDT/100).EQ.0) THEN
C                  WRITE(OP_STRING,'(I5,''/'',I5)') nd,NDT
C                  CALL WRITES(IODI,OP_STRING,ERROR,*200)
C                ENDIF
! Calculate first function value
                DO nj=1,3
                  UNDEF(nj)=ZD(nj,nd)
                ENDDO

                CALL DEFORM_POINT(CURVCORRECT,IBT,IDO,INP,LD(nd),
     '            NBH,NBJ,
     '            NEELEM,NHE,NKJE,NKHE,NPF,NPNE,NRE,NVHE,NVJE,NW,
     '            SE,OK,XA,DEF,XE,XID(1,nd),XP,ZA,UNDEF,ZP,ERROR,*200)
                CONVERG=.FALSE.

                IF(OK) THEN     ! check if point is present
                ! calculate image plane constant
                  CC=ZD(1,nd)*WD(1,nd)+ZD(2,nd)*WD(2,nd)+ZD(3,nd)*
     '              WD(3,nd)
                  X2=0.0d0
                  F=-CC
                ! calcualtes dotprod of normal and deformed point
                  DO nj=1,3
                    F=F+WD(nj,nd)*DEF(nj)
                  ENDDO          !nj
                  ! Calculate second function value
                  X1=0.001d0    ! starting delta
                  FL=-CC
                  DO nj=1,3
                    UNDEF(nj)=ZD(nj,nd)+X1*WD(nj,nd)
                  ENDDO
                  CALL DEFORM_POINT(CURVCORRECT,IBT,IDO,INP,LD(nd),NBH,
     '              NBJ,NEELEM,NHE,NKJE,NKHE,NPF,NPNE,NRE,NVHE,
     '              NVJE,NW,SE,OK,XA,DEF,XE,XID(1,nd),XP,ZA,UNDEF,ZP,
     '              ERROR,*200)
                  IF(OK) THEN
                    DO nj=1,3
                      FL=FL+WD(nj,nd)*DEF(nj)
                    ENDDO
                    IF(DABS(FL).LT.DABS(F)) THEN
                      RTSEC=X1
                      XL=X2
                      SWAP=FL
                      FL=F
                      F=SWAP
                    ELSE
                      XL=X1
                      RTSEC=X2
                    ENDIF
                    J=0
                    MAXIT=50
                    DO WHILE (J.LE.MAXIT.AND..NOT.CONVERG)
                      J=J+1
                      DX=(XL-RTSEC)*F/(F-FL)
                      XL=RTSEC
                      FL=F
                      RTSEC=RTSEC+DX
                      F=-CC
                      DO nj=1,3
                        UNDEF(nj)=ZD(nj,nd)+RTSEC*WD(nj,nd)
                      ENDDO
                      CALL DEFORM_POINT(CURVCORRECT,IBT,IDO,
     '                  INP,LD(nd),NBH,NBJ,NEELEM,NHE,NKJE,NKHE,NPF,
     '                  NPNE,NRE,NVHE,NVJE,NW,SE,OK,XA,DEF,XE,
     '                  XID(1,nd),
     '                  XP,ZA,UNDEF,ZP,ERROR,*200)
                      IF(OK) THEN
                        DO nj=1,3
                          F=F+WD(nj,nd)*DEF(nj)
                        ENDDO
                      ELSE
                        J=MAXIT+1
                      ENDIF
                      IF(DABS(DX).LT.0.000001d0.OR.F.LT.0.000001d0) THEN
                      !IF(DABS(DX).LT.0.000001d0.OR.F.EQ.0.d0) THEN
                        CONVERG=.TRUE.
                      ENDIF
C                      IF(DOP) THEN
C                        WRITE(OP_STRING,'(3E11.3,'' -> '',3E11.3)')
C     '                    ZD(1,nd),ZD(2,nd),ZD(3,nd),DEF(1),DEF(2),
C     '                    DEF(3)
C                        CALL WRITES(IODI,OP_STRING,ERROR,*200)
C                      ENDIF
                    ENDDO
                  ENDIF
                ENDIF
                IF(CONVERG.EQV..TRUE.) THEN
                  DO nj=1,3
                    ZD(nj,nd)=DEF(nj)
                  ENDDO
                ELSE
C point is scrapped (doesn't project)
C                     WRITE(OP_STRING,'(I5,'' '',E12.3)') J,DX
C                     CALL WRITES(IODI,OP_STRING,ERROR,*200)
                  DO nj=1,3
                    WD(nj,nd)=0.0d0
                    ZD(nj,nd)=-1000.0d0
                  ENDDO
                ENDIF
                GO TO 202
C             This statement is designed to be skipped if no error
C             occurs. However if a error occurs within a subroutine
C             the alternate return points to line 200 to set the flag
 200            CONTINUE
C$OMP CRITICAL(DEDATA1)
                ERROR_F=.TRUE.
                WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
                CALL WRITES(IOER,OP_STRING,ERROR,*201)
                WRITE(OP_STRING,'(/'' >>An error occurred - '
     '            //'results may be unreliable!'')')
                CALL WRITES(IODI,OP_STRING,ERROR,*201)
 201            CONTINUE
C$OMP END CRITICAL(DEDATA1)
 202            CONTINUE
             ENDIF
           ENDDO
C C$OMP      END PARALLEL DO
          ELSEIF(TYPE(1:10).EQ.'CENTRELINE') THEN
            CALL CENTRE_LINE(AXIS_DIV,AXIS_XI,CIRC_XI,IBT,IDO,INP,NBJ,            
     '        NDDATA,NDP,NEELEM,NKJE,NPF,NPNE,nr,NVJE,NXI,SE,            
     '        START_ELEMENT,WD,WG,XA,XE,XIG,XP,ZD,ERROR,*9999)            
          ELSEIF(TYPE(1:4).EQ.'AREA') THEN
            CALL AREAFIELD(IBT,IDO,INP,NBJ,NCENTRE_REG,NDDATA,NDP,            
     '        NEELEM,NKJE,NPF,NPOINTS_PER_ELEM,NPNE,NSURFACE_REG,NVJE,
     '        SE,RADIUS_OF_INTEREST,TOLERANCE,WD,XA,XE,XP,ZD,
     '        ERROR,*9999)
          ELSEIF(TYPE(1:13).EQ.'CYLINDER_TREE') THEN
            CALL CYLINDER_DATA(NBJ,NDDATA,NDP,NENP,NPLIST,NPNE,nr,NVJE,
     '         density,WD,XP,ZD,FIRST,ERROR,*9999)
          ELSEIF(TYPE(1:13).EQ.'FROM_NODE') THEN
             CALL DENODE_FROM_DATA(LD,NBJ,NDDATA,NDP,NENP,NPLIST,NPNE,
     &         nr,WD,XID,XP,ZD,ERROR,*9999)
          ELSEIF(TYPE(1:11).EQ.'FILL_VOLUME')THEN
            CALL GNNEZD(IBT,IDO,INP,LD,NBJ,NDDATA,NELIST,NKJE,NPF,
     &        NPNE,nr,NVJE,PG,RG,SE,WD,WG,XA,X_DISCRETISATION,XE,XG,XID,
     &        XP,ZD,SPREAD_TYPE,ERROR,*9999)
            CALL DSTATS(LD,NDDL,NDLT,ERROR,*9999)
          ENDIF
        ENDIF

C cpb 28/10/95 Adding NDDATA
C        NDDATA(0,0)=NDT
C        NDDATA(0,1)=NDT
C        DO nr=0,1
C          DO nd=1,NDT
C            NDDATA(nd,nr)=nd
C          ENDDO
C        ENDDO !nr

C LKC 4-JUL-1999 Generalisation for different regions
!  need to change here when there are multiple data regions

        NDDATA(0,0)=NDT
        NDDATA(0,NRLIST(1))=NDT
C LKC 29-MAY-20000 should be looping over NRLIST(0)
C        DO nrr=1,NRLIST(1)
        DO nrr=1,NRLIST(0)
          nr=NRLIST(nrr)
          DO nd=1,NDDATA(0,nr)
            NDDATA(nd,nr)=nd
          ENDDO
        ENDDO !nrr (nr)

        CALL_DATA=.TRUE.
      ENDIF

      CALL EXITS('DEDATA')
      RETURN
 9999 CALL ERRORS('DEDATA',ERROR)
      IF(ISBINFILEOPEN(IOFILE1))
     '  CALL BINARYCLOSEFILE(IOFILE1,ERR,CERROR)
      CALL EXITS('DEDATA')
      RETURN 1
      END


