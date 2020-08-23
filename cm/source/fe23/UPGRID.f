      SUBROUTINE UPGRID(IBT,IDO,INP,ICQS_SPATIAL,IRCQS_SPATIAL,NAN,NAQ,
     '  NBH,NBJ,NEELEM,NELIST,NENQ,NGAP,NHE,NHP,NKH,NKHE,NKJE,NLL,
     '  NLQ,NPF,NPL,NPNE,NPNODE,NQGP,NQLIST,NQNE,NQXI,NQS,NRLIST,
     '  NVHE,NVHP,NVJE,NW,NWQ,NXLIST,NXQ,NYNE,NYNP,AQ,CE,CP,CQ,
     '  CURVCORRECT,DL,DNUDXQ,DXDXIQ,DXDXIQ2,FEXT,GCHQ,GUQ,PG,PROPQ,
     '  RCQS_SPATIAL,SE,XA,XE,XG,XIQ,XP,XQ,YG,YP,YQ,YQS,ZA,ZE,ZG,ZP,
     '  STRING,RET_ERROR,*)          

C#### Subroutine: UPGRID
C###  Description:
C###    UPGRID updates grid point parameters, including geometric
C###    positions, domain metrics, material values (conductivity 
C###    tensors), potentials and source terms for the bidomain model.

C#### Variable: XQ(nj,nq)
C###  Type: REAL*8
C###  Set_up: UPGRID,CALC_LATT_XQ
C###  Description:
C###    XQ(nj,nq) is rectangular cartesian coordinates of grid point nq.

      IMPLICIT NONE

      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='UPGRID')

      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cellml.cmn'
      INCLUDE 'error0.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'maqloc00.cmn'
      INCLUDE 'maqloc00.inc'
      INCLUDE 'nqloc00.inc'
      INCLUDE 'time02.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     &  ICQS_SPATIAL(NQISVM,NQM),IRCQS_SPATIAL(0:NQRSVM,NQVM),
     &  NAN(NIM,NAM,NBFM),NAQ(NQM,NAM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     &  INP(NNM,NIM,NBFM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     &  NENQ(0:8,NQM),NGAP(NIM,NBM),NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),
     &  NKH(NHM,NPM,NCM,0:NRM),NKHE(NKM,NNM,NHM,NEM),
     &  NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),NLQ(NQM),NPF(9,NFM),
     &  NPL(5,0:3,NLM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     &  NQGP(0:NQGM,NQM),NQLIST(0:NQM),NQNE(NEQM,NQEM),NQS(NEQM),
     &  NQXI(0:NIM,NQSCM),NRLIST(0:NRM),NVHE(NNM,NBFM,NHM,NEM),
     &  NVHP(NHM,NPM,NCM,0:NRM),NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM),
     &  NWQ(8,0:NQM,NAM),NXQ(-NIM:NIM,0:4,0:NQM,NAM),NXLIST(0:NXM),
     &  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     &  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 AQ(NMAQM,NQM),CE(NMM,NEM,NXM),CP(NMM,NPM,NXM),
     &  CQ(NMM,NQM,NXM),CURVCORRECT(2,2,NNM,NEM),DL(3,NLM),
     &  DNUDXQ(3,3,NQM),DXDXIQ(3,3,NQM),DXDXIQ2(3,3,NQM),
     &  FEXT(NIFEXTM,NGM,NEM),GCHQ(3,NQM),GUQ(3,3,NQM),
     &  PG(NSM,NUM,NGM,NBM),PROPQ(3,3,4,2,NQM,NXM),
     &  RCQS_SPATIAL(NQRSVM,NQM),SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),
     &  XE(NSM,NJM),XG(NJM,NUM),XIQ(NIM,NQM),XP(NKM,NVM,NJM,NPM),
     &  XQ(NJM,NQM),XQD(NJM,NUM),YG(NIYGM,NGM,NEM),YP(NYM,NIYM,NXM),
     '  YQ(NYQM,NIQM,NAM,NXM),YQS(NIQSM,NQM),
     '  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZG(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER RET_ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER atime,cellvarcmpt,CELL_VARIANT_NUM,CMPTLIST(0:6),
     '  field_index,i,IBEG,ICHAR,ICOOR_LOCAL,IDO1,IDO2,
     '  IDO3,IEND,IFNTYP,II,IJ,IK,iL,INFO,IFROMC,j,jb,k,l,maqx,maqy,
     '  maqz,mLL,mRR,N3CO,na,na_level,nb,nbb,NBJ_LOCAL(40),nbq,nc,ncmpt,
     '  NCW,ndx,ne,nee,neq,ngBL,ngBR,ngTL,ngTR,NG_row,nh,nhx,nicmpt,
     '  nii2,nij2,nik2,ni,nii,nij,nij1,nik,niqV,niqs,niqSAC,NITB,nj,
     '  njj1,njj2,nk,node1,node2,no_dxi,nn,noelem,NU1(0:3),no_nelist,
     '  no_nrlist,NOQUES,nq,nq1,nq2,nq_adj,nqI,nqL,nqq,nqrsv,nqv,NQ_row,
     '  n1list,nr,nr1,nr2,nrr,nTT,nu,nx,nxc,nx_d,nx2,POSITION(4),
     '  RESET_GRID_MAX,SCHEME,nnelem(NNM)

C!!! KAT: not merging from Oxford
C     '  ,ss_index
      PARAMETER (RESET_GRID_MAX=500)
      INTEGER ATXI,RESET_GRID_PT(0:RESET_GRID_MAX),
     '  NQLIST2(0:NQM),nyqs,nrcqs
      PARAMETER (NCW=35) !CW must be dimen.d the same size as CE array
      REAL*8 A_VECTOR(3),AZ,AZL(3,3),AZU(3,3),B_VECTOR(3),CHTOFF(3,3,3),
     &  Current,C_VECTOR(3),CW(NCW),DIFFUSION,DBM(3,3,3),delx,DETM,
     &  df1dxi,df2dxi,df3dxi,df4dxi,dNudXi(3,3),dS,dxdxi,dydxi,dzdxi,
     &  dXdNu(3,3),dXidNu(3,3),DXIX(3,3),Esac,Ext,DEFORMATION_VALUE,
     &  GX,GY,GL(3,3),GU(3,3),f1,f2,f3,f4,maxextn,maxextn1,maxextn2,PH3,
     &  pole1,pole2,PXI,RDUMMY(1),RFROMC,S,S0,S1,SAC_scale_factor,
     &  SAC_threshold,SIGMAMINUSE,SIGMAMINUSI,SIGMAPLUSE,SIGMAPLUSI,SUM,
     &  V,Temp_var,X_TEMP(3),X3G(4,3),XI(3),Xi1,Xi2,Xi3,XIPOS(3),
     &  XQ_TEMP(3),YIELD_TOL,const_value,FIELD_VAL(NNM)
      CHARACTER ERROR*(ERRSTRLEN),CHAR1*3,TYPE*22
      LOGICAL ADAPTIVE,ALL_LEVELS,ALL_REGIONS,CALC_ZE,CBBREV,DEFORMED,
     &     DIRN,DEFORMATION_RCQS,DEFORMATION_YQS,ERROR_FLAG,FILEIP,
     &     FOUND,FOUNDGP,INLIST,QUADBASIS,MAPG,UP_Vm,YQS_FIELD,YQS_CONST
     &     ,RCQS_FIELD,RCQS_CONST,FROM_RCQS

!     Functions
      INTEGER LEN_TRIM          
C     ,LOCAL_GAUSS_DIM
C OR 10-08-2005:
C     Changed LOCAL_GAUSS_DIM to NQEM (maximum number of grid points)
C     PARAMETER (LOCAL_GAUSS_DIM=64)
      REAL*8 PFXI,AWG(NQEM),AXIG(3,NQEM)

      DATA NU1/1,2,4,7/
      DATA YIELD_TOL/0.90d0/ !max extension ratio tolerance
C VJ Changed EXTENSION_RATIO* to DEFORMATION* for better meaning for variable
      CALL ENTERS(ROUTINENAME,*9999)
      
      nc=1 !temporary
      ICHAR=999
      ERROR_FLAG=.FALSE.

C OR 10-08-2005: Start initalize local variables to prevent warnings.
C      IFNTYP=0
C      na_level=0
C      nbq=0
C      nyqs=1
C      nrcqs=3
C      ALL_LEVELS=.FALSE.
C      CALC_ZE=.FALSE.
      DEFORMATION_YQS=.TRUE.
      DEFORMATION_RCQS=.FALSE.
C OR 10-08-2005: End ininitalizing
     
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM update grid yqs 
C###  Parameter:      <(yqs #)[1]>
C###    YQS index to update.
C###  Parameter:      <field (#)[1]>
C###    To update yqs form evaluated field
C###  Parameter:      <constant (#) [0.0]>
C###    To update yqs by constant value
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Parameter:      <from_rcqs index #>
C###    Copy values from RCQS_SPATIAL to YQS. This is used
C###    to print out the RCQS_SPATIAL values for debugging.
C###    Use fem inquire cell_variable to discover the index
C###    for a named cellml variable.
C###  Description:
C###    This command evaluates the specified field at grid
C###    points and stores the values in the specified YQS index

        OP_STRING(1)=STRING(1:IEND)//' yqs'
        OP_STRING(2)=BLANK(1:IEND)//'<(yqs index #)[1]'
        OP_STRING(3)=BLANK(1:IEND)//'<field (#)[1]>'
        OP_STRING(4)=BLANK(1:IEND)//'<constant (#)[1]>'
        OP_STRING(5)=BLANK(1:IEND)//'<region (#s/all)[1]>'
        OP_STRING(6)=BLANK(1:IEND)//'<class #[1]>'
        OP_STRING(7)=BLANK(1:IEND)//'<from_rcqs index #>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update grid rcqs
C###  Parameter:      <(rcqs index #)[1]>
C###    RQCS index to update.
C###  Parameter:      <field (#)[1]>
C###    To update rcqs form evaluated field
C###  Parameter:      <constant (#)[0.0]>
C###    To update rcqs by constant value
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    This command evaluates the specified field at grid
C###    points and stores the values in the specified YQS index
C###    It is important to note that the specific RCQS variable
C###    needs to be spatially variant (ipcell file)
        
        OP_STRING(1)=STRING(1:IEND)//' rcqs'
        OP_STRING(2)=BLANK(1:IEND)//'<(rcqs index #)[1]>'
        OP_STRING(3)=BLANK(1:IEND)//'<field (#)[1]>'
        OP_STRING(4)=BLANK(1:IEND)//'<constant (#)[1]>'        
        OP_STRING(5)=BLANK(1:IEND)//'<region (#s/all)[1]>'
        OP_STRING(6)=BLANK(1:IEND)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update grid geometry
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    This command calculates the geometric positions of the grid
C###    points and stores them in XQ. The grid is spaced evenly in
C###    material coordinates. The geometric positions are calculated
C###    by interpolation using the finite element basis functions used
C###    to describe the elements in which the grid points lie.

        OP_STRING(1)=STRING(1:IEND)//' geometry'
        OP_STRING(2)=BLANK(1:IEND)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:IEND)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update grid geometry deformed
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <from_class #[1]>
C###    Specify the class number of the problem to update from
C###  Description:
C###    This command changes the geometric positions of the grid
C###    points when the underlying finite element mesh has undergone
C###    a deformation. Because the grid points don't move in material
C###    space, they will deform with the finite elements.

        OP_STRING(1)=STRING(1:IEND)//' geometry deformed'
        OP_STRING(2)=BLANK(1:IEND)//'<from_class #[1]>'
        OP_STRING(3)=BLANK(1:IEND)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:IEND)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update grid metric
C###  Parameter:      <level #)>[all]
C###    Specify the grid level to update. Default is all levels.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    This command generates the metric tensor information
C###    associated with a grid. This information is necessary
C###    for activation solutions.

        OP_STRING(1)=STRING(1:IEND)//' metric'
        OP_STRING(2)=BLANK(1:IEND)//'<level #[all]>'
        OP_STRING(3)=BLANK(1:IEND)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:IEND)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update grid metric deformed
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <from_class #[1]>
C###    Specify the class number of the problem to update from
C###  Description:
C###    This command updates the metric tensor information
C###    associated with a grid after a deformation has occured. The
C###    metric properties are updated for the deformed fibre angles.

        OP_STRING(1)=STRING(1:IEND)//' metric deformed'
        OP_STRING(2)=BLANK(1:IEND)//'<from_class #[1]>'
        OP_STRING(3)=BLANK(1:IEND)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:IEND)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update grid material
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Update diffusivity tensor with metric information.

        OP_STRING(1)=STRING(1:IEND)//' material'
        OP_STRING(2)=BLANK(1:IEND)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:IEND)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update grid solution
C###  Parameter:      <from_region #[2]>
C###    Specify the region from which the solution is to
C###    be updated.
C###  Parameter:      <to_region #[1]>
C###    Specify the region to which the solution is to be
C###    updated.
C###  Parameter:      <class #[1]>
C###    Specify the solution class where the solution is
C###    stored.
C###  Description:
C###    Update a grid solution from one region to another.
C###    This may only be done after a call to define the coupling
C###    between the two regions.

        OP_STRING(1)=STRING(1:IEND)//' solution'
        OP_STRING(2)=BLANK(1:IEND)//'<from_region #[2]>'
        OP_STRING(3)=BLANK(1:IEND)//'<to_region #[1]>'
        OP_STRING(4)=BLANK(1:IEND)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update grid strain
C###  Parameter:      <SAC_extension THRESHOLD[1]>
C###   Specify the stretch activated channal threshold
C###  Parameter:      <SAC_scale FACTOR[1]>
C###    Specify the stretch activated channel scale factor
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Uses the fibre extension ratio computed from FEXT(1) in
C###    surrounding gauss pts to set a SAC current in YQ(7).
C###    Currently FHN & VCD models only.
C###    YQ(7) is subtracted from Iapp in subroutine GEN_INT_RHS in FE30.
C###    SAC_extension sets the threshold for the SAC
C###    SAC_scale sets the scale factor on the SAC current

        OP_STRING(1)=STRING(1:IEND)//' strain'
        OP_STRING(2)=BLANK(1:IEND)//'<SAC_extension THRESHOLD[1]>'
        OP_STRING(3)=BLANK(1:IEND)//'<SAC_scale FACTOR[1]>'
        OP_STRING(4)=BLANK(1:IEND)//'<region (#s/all)[1]>'
        OP_STRING(5)=BLANK(1:IEND)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C news MPN 28Jun2000: calc'ing values of extn ratio at grid pts
C                    This may eventually be used to supercede the
C                    STRAIN option immediately above
C                    (still under debate!)
C news VJ 8Dec2003: added option of calc'ing values of green strain at grid pts

C#### Command: FEM update grid extension_ratio/green_strain
C###  Parameter:      <component #[1]>
C###    Specifies the extension ratio component, where 1/2/3
C###    are fibre/sheet/sheet-normal components and 4/5/6 are shear
C###    extension ratio components (f-s/f-sn/s-sn).
C###  Parameter:      <yqs/rcqs #s[1]>
C###    Specifies the list of grid point variable indices for each
C###    cellular variant.Each index for each variant must be specified
C###  Parameter:      <all_variants>
C###    if option selected, only one set of indices matching the number
C###    of components for one cell variant needs to be specified. 
C###    The indices for all other variants is assumed to be the same
C###  Parameter:      <no_ze_calc>
C###    if option selected, do not reinterpolate nodal variables to element
C###    variables. It is assumed that this has already been done
C###  Parameter:      <atgrid/atxi <xi_1 VALUE[0.5]>
C###    <xi_2 VALUE[0.5]><xi_3 VALUE[0.5]>[atgrid]>
C###    Interpolate to the selected xi-point.
C###    Note this is a hack as it temporarily stores the selected xi
C###    point values in the first grid point position.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the regions to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use for extension
C###    ratio computations.
C###  Description:
C###    Computes the specified extension ratio component (with respect
C###    to fibre coordinates) at grid points and stores them in the
C###    grid point variable array, YQS or RCQS. Different indices can be
C###    used for each cellular variant defined. If the RCQS array is
C###    used, then the appropriate cellular real parameter MUST be
C###    spatially varying and the calculated values will actually be
C###    stored in the cellular spatial array RCQS_SPATIAL.

        OP_STRING(1)=STRING(1:IEND)//' extension_ratio/green_strain'
        OP_STRING(2)=BLANK(1:IEND)//'<element (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:IEND)//'<component #[1]>'
        OP_STRING(4)=BLANK(1:IEND)//'<yqs/rcqs #s[1]>'
        OP_STRING(5)=BLANK(1:IEND)//'<all_variants>'
        OP_STRING(6)=BLANK(1:IEND)//'<no_ze_calc>'
        OP_STRING(7)=BLANK(1:IEND)//'<atgrid/atxi <xi_1 VALUE[0.5]>'
     '    //'<xi_2 VALUE[0.5]><xi_3 VALUE[0.5]>[atgrid]>'
        OP_STRING(8)=BLANK(1:IEND)//'<region (#s/all)[1]>'
        OP_STRING(9)=BLANK(1:IEND)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C newe MPN 28Jun2000

C---------------------------------------------------------------------

C#### Command: FEM update grid source
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description: Updates the source information at grid points

        OP_STRING(1)=STRING(1:IEND)//' source'
        OP_STRING(2)=BLANK(1:IEND)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:IEND)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update grid adaptive_levels
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description: Updates the source information at grid points

        OP_STRING(1)=STRING(1:IEND)//' adaptive_levels'
        OP_STRING(2)=BLANK(1:IEND)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:IEND)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update grid activation_time
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description: Updates the activation times for grid points into
C###   the YQ array with an NIQ of 4 so time information at grid points
C###   may be exported to cmgui or plotted as a field

        OP_STRING(1)=STRING(1:IEND)//' activation_time'
        OP_STRING(2)=BLANK(1:IEND)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:IEND)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update grid time
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description: Update previous time step solution for multigrid
C###    transient heat eqtn

        OP_STRING(1)=STRING(1:IEND)//' time'
        OP_STRING(2)=BLANK(1:IEND)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:IEND)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update grid connectivity
C###  Parameter:      <in_region #[1]>
C###    Specify the region containing the host grid.
C###  Parameter:      <from_region #[2]>
C###    Specify the region containing elements which cleave the grid.
C###  Description:
C###    Update grid connectivity array NXQ by disconnecting grid pts
C###    which lie on either side of a cleavage plane.
C###    Note NWQ contains information about the adjacent grid points.

        OP_STRING(1)=STRING(1:IEND)//' connectivity'
        OP_STRING(2)=BLANK(1:IEND)//'<in_region #[1]>'
        OP_STRING(3)=BLANK(1:IEND)//'<from_region #[2]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update grid connectivity
C###  Parameter:      <grid (#s/name)[none]>
C###    Specify which grid points to cleave at
C###  Parameter:      <direction #[1]>
C###    Specify the xi direction in which to cleave the mesh.
C###  Description:
C###    Update grid connectivity array NXQ by disconnecting grid pts
C###    which lie on either side of a cleavage plane.
C###    Note NWQ contains information about the adjacent grid points.

        OP_STRING(1)=STRING(1:IEND)//' connectivity'
        OP_STRING(2)=BLANK(1:IEND)//'<grid (#s/name)[none]>'
        OP_STRING(3)=BLANK(1:IEND)//'<direction #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update grid stimulus
C###  Parameter:      <grid (#s/name)[none]>
C###    Specify the grid points at which the stimulus needs to be
C###    calculated. Only the list of interstitial grid points are
C###    needed.
C###  Parameter:      <diffusion #[0.0]>
C###    Specify the diffusion coefficient to multiply the potential
C###    difference by.   
C###  Parameter:      <direction #[all]>
C###    Specify a single direction in which to update the psuedo
C###    stimulus magnitude, by default all directions are updated.
C###  Parameter:      <map_grid  (#s/name[none]>
C###    Specify a list of grid points that is maping to another list of grid 
C###    points.
C###  Parameter:      <up_potential  (#s/name[none]>
C###    Specify a list of ICC grid points that is used to calculate the 
C###    potential gradient which is specify in Alievs cell model
C###  Description:
C###    Update grid stimulus calculates the potential difference
C###    based pseudo stimulus magnitude that occurs between ICC and
C###    smooth muscle cells in GI problems.

        OP_STRING(1)=STRING(1:IEND)//' stimulus'
        OP_STRING(2)=BLANK(1:IEND)//'<grid (#s/name)[none]>'
        OP_STRING(3)=BLANK(1:IEND)//'<diffusion #[0.0]>'
        OP_STRING(4)=BLANK(1:IEND)//'<direction #[all]>'
        OP_STRING(5)=BLANK(1:IEND)//'<map_grid grid (#s/name)[none]>'
        OP_STRING(6)=BLANK(1:IEND)//'<up_potential>'
     
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
C---------------------------------------------------------------------
      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc',ROUTINENAME,ERROR,*9999)
      ELSE
        IF(CBBREV(CO,'ADAPTIVE_LEVELS'  ,3,noco+1,NTCO,N3CO)) THEN
          TYPE='ADAPTIVE'
        ELSE IF(CBBREV(CO,'ACTIVATION_TIMES',3,noco+1,NTCO,N3CO)) THEN
          TYPE='ACTIVTIMES'
        ELSE IF(CBBREV(CO,'CONNECTIVITY',3,noco+1,NTCO,N3CO)) THEN
          TYPE='CONNECTIVITY'
        ELSE IF(CBBREV(CO,'GEOMETRY'    ,3,noco+1,NTCO,N3CO)) THEN
          TYPE='GEOMETRY'
        ELSE IF(CBBREV(CO,'MATERIAL'    ,3,noco+1,NTCO,N3CO)) THEN
          TYPE='MATERIAL'
        ELSE IF(CBBREV(CO,'METRIC'      ,3,noco+1,NTCO,N3CO)) THEN
          TYPE='METRIC'
        ELSE IF(CBBREV(CO,'SOLUTION'    ,3,noco+1,NTCO,N3CO)) THEN
         TYPE='SOLUTION'
        ELSE IF(CBBREV(CO,'SOURCE'      ,3,noco+1,NTCO,N3CO)) THEN
          TYPE='SOURCE'
        ELSE IF(CBBREV(CO,'STIMULUS'    ,3,noco+1,NTCO,N3CO)) THEN
          TYPE='STIMULUS'
        ELSE IF(CBBREV(CO,'STRAIN'      ,3,noco+1,NTCO,N3CO)) THEN
          TYPE='STRAIN'
C news MPN 28Jun2000: calc'ing values of extn ratio at grid pts
        ELSE IF(CBBREV(CO,'EXTENSION_RATIO',3,noco+1,NTCO,N3CO)) THEN
          TYPE='EXTENSION_RATIO'
C newe MPN 28Jun2000
C news VJ 8Dec2003: calc'ing values of green strain at grid pts
        ELSE IF(CBBREV(CO,'GREEN_STRAIN',3,noco+1,NTCO,N3CO)) THEN
          TYPE='GREEN_STRAIN'
C newe VJ 8Dec2003
        ELSE IF(CBBREV(CO,'TIME'        ,1,noco+1,NTCO,N3CO)) THEN
          TYPE='TIME'
C OR 23/04/2006
        ELSE IF(CBBREV(CO,'YQS'         ,3,noco+1,NTCO,N3CO)) THEN
          TYPE='YQS'
          nyqs=IFROMC(CO(N3CO+1))
        ELSE IF(CBBREV(CO,'RCQS'        ,3,noco+1,NTCO,N3CO)) THEN
          TYPE='RCQS'
          nrcqs=IFROMC(CO(N3CO+1))
        ENDIF
C     --type------------------------------------------------------------
C
C     OR 23/04/2006 INTRODUCE NEW FUNCTIONALITY BY CHANGING THE YQS
C          VALUES BY INTERPOLATION FROM A FIELD OR BY A CONSTANT VALUE
C
        IF(TYPE(1:3).EQ.'YQS') THEN
          FROM_RCQS=.FALSE.
          YQS_FIELD=.FALSE.
          YQS_CONST=.FALSE.
          field_index=1
          const_value=0.0d1
          IF(CBBREV(CO,'FIELD',3,noco+1,NTCO,N3CO)) THEN
            YQS_FIELD=.TRUE.
            field_index=IFROMC(CO(N3CO+1))
          ELSE IF(CBBREV(CO,'CONSTANT',2,noco+1,NTCO,N3CO)) THEN 
            YQS_CONST=.TRUE.
            const_value=RFROMC(CO(N3CO+1))
          ELSE IF(CBBREV(CO,'FROM_RCQS',3,noco+1,NTCO,N3CO)) THEN
            IF(CBBREV(CO,'INDEX',3,noco+1,NTCO,N3CO)) THEN
               FROM_RCQS=.TRUE.
               nrcqs=IFROMC(CO(N3CO+1))
            ENDIF
          ENDIF
          IF (((YQS_FIELD.EQV..TRUE.).AND.(YQS_CONST.EQV..TRUE.)).OR.
     &         ((YQS_FIELD.EQV..FALSE.).AND.(YQS_CONST.EQV..FALSE.))
     &          .AND.FROM_RCQS.EQV..FALSE.)
     &         THEN
            CALL ASSERT(.FALSE.,'>> FEM update grid yqs #'//
     '          ' from_rcqs # command needs an array index'//
     '          ' as argument',ERROR,*9999)
          ENDIF
          
          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS, ERROR,
     &         *9999)
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nxc=NXLIST(1)

          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            DO no_nelist=1,NEELEM(0,nr)
              ne=NEELEM(no_nelist,nr)

              IF (YQS_FIELD) THEN
                nj=NJ_LOC(NJL_FIEL,field_index,nr)
                nb=NBJ(nj,ne)
                NITB=NIT(nb)
              ENDIF
              SCHEME=NQS(ne)
              II=MAX(1,NQXI(1,SCHEME))
              IJ=1
              IK=1
              IF(NQXI(0,SCHEME).GT.1) IJ=MAX(1,NQXI(2,SCHEME)) 
              IF(NQXI(0,SCHEME).GT.2) IK=MAX(1,NQXI(3,SCHEME))
              
              IF(YQS_FIELD) THEN
                CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1), NPNE(1,1,ne
     &               ),nr,NVJE(1,1,1,ne), SE(1,1,ne),XA,XE,XP,ERROR,
     &               *9999)
              ENDIF
C     Loop over the grid points in each element
              DO nik=1,IK
                DO nij=1,IJ
                  DO nii=1,II
                    neq=nii+((nij-1)*NQXI(1,SCHEME)) !local grid pt #
                    IF(NQXI(0,SCHEME).GT.1) neq=neq+((nik-1)* NQXI(1
     &                   ,SCHEME)*NQXI(2,SCHEME))
                    nq=NQNE(ne,neq) !global grid pt #
C     Evaluate each grid point only once
                    IF(NENQ(1,nq).EQ.ne) THEN
C     Local Xi coords of grid point nq in element ne
                      IF (YQS_FIELD) THEN
                        IF(II.NE.1) XI(1)=DBLE(nii-1)/DBLE(II-1)
                        IF(IJ.NE.1) XI(2)=DBLE(nij-1)/DBLE(IJ-1)
                        IF(IK.NE.1) XI(3)=DBLE(nik-1)/DBLE(IK-1)
                        YQS(nyqs,nq)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     &                       INP(1,1,nb),nb,1,Xi,XE(1,nj))
                      ELSEIF (FROM_RCQS) THEN
                        YQS(nyqs,nq)=RCQS_SPATIAL(nrcqs-2,nq)
                      ELSE 
                        YQS(nyqs,nq)=const_value
                      ENDIF
                    ENDIF       !NENQ(1,nq).EQ.ne
                  ENDDO         !nii
                ENDDO           !nij
              ENDDO             !nik
            ENDDO               !ne
          ENDDO                 !nr
          
C--   type---------------------------------------------------------------
C
C     OR 23/04/2006 INTRODUCE NEW FUNCTIONALITY BY CHANGING THE RCQS
C          VALUES BY INTERPOLATION FROM A FIELD OR BY A CONSTANT VALUE
C
        ELSE IF(TYPE(1:4).EQ.'RCQS') THEN
          RCQS_FIELD=.FALSE.
          RCQS_CONST=.FALSE.
          field_index=1
          const_value=0.0d1
          IF(CBBREV(CO,'FIELD',3,noco+1,NTCO,N3CO)) THEN
            RCQS_FIELD=.TRUE.
            field_index=IFROMC(CO(N3CO+1))
          ELSE IF(CBBREV(CO,'CONST',2,noco+1,NTCO,N3CO)) THEN 
            RCQS_CONST=.TRUE.
            const_value=RFROMC(CO(N3CO+1))
          ENDIF
          IF (((RCQS_FIELD.EQV..TRUE.).AND.(RCQS_CONST.EQV..TRUE.)).OR.
     &         ((RCQS_FIELD.EQV..FALSE.).AND.(RCQS_CONST.EQV..FALSE.)))
     &         THEN
            CALL ASSERT(.FALSE.,'>> FEM update grid rcqs command'//
     &        ' needs either a field or a constant value as'//
     &        ' argument',ERROR,*9999)
          ENDIF
          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS, ERROR,
     &         *9999)
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nxc=NXLIST(1)

          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            DO no_nelist=1,NEELEM(0,nr)
              ne=NEELEM(no_nelist,nr)
              
              IF(RCQS_FIELD) THEN
                nj=NJ_LOC(NJL_FIEL,field_index,nr)                
                CALL ASSERT(nj.GT.0,
     '               'field does not exist',ERROR,*9999)
                nb=NBJ(nj,ne)
                NITB=NIT(nb)
              ENDIF
              SCHEME=NQS(ne)
              II=MAX(1,NQXI(1,SCHEME))
              IJ=1
              IK=1
              IF(NQXI(0,SCHEME).GT.1) IJ=MAX(1,NQXI(2,SCHEME)) 
              IF(NQXI(0,SCHEME).GT.2) IK=MAX(1,NQXI(3,SCHEME))
              
              IF(RCQS_FIELD) THEN
                DO nn=1,NNT(nb)
                  nnelem(nn)=NPNE(nn,nb,ne)
                  FIELD_VAL(nn)=XP(1,1,nj,nnelem(nn)) ! Using only version 1
                ENDDO           !nn
              ENDIF
C     Loop over the grid points in each element
              DO nik=1,IK
                DO nij=1,IJ
                  DO nii=1,II
                    neq=nii+((nij-1)*NQXI(1,SCHEME)) !local grid pt #
                    IF(NQXI(0,SCHEME).GT.1) neq=neq+((nik-1)* NQXI(1
     &                   ,SCHEME)*NQXI(2,SCHEME))
                    nq=NQNE(ne,neq) !global grid pt #
C     Evaluate each grid point only once
                    IF(NENQ(1,nq).EQ.ne) THEN
C     Local Xi coords of grid point nq in element ne
                      IF (RCQS_FIELD) THEN
                        XI(1)=XQ(1,nq) 
                        XI(2)=XQ(2,nq)
                        XI(3)=XQ(3,nq)  
C     E11 is always the first entry and it is assigned to PARAM(3). This
C     translates to the first index of rcqs_spatial. Hence the subtraction of 2.
                        RCQS_SPATIAL(nrcqs-2,nq)=PXI(IBT(1,1,nb),IDO(1,1
     &                       ,0,nb),INP(1,1,nb),nb,1,Xi,FIELD_VAL)
                      ELSE
                        RCQS_SPATIAL(nrcqs-2,nq)=const_value
                      ENDIF
                    ENDIF       !NENQ(1,nq).EQ.ne
                  ENDDO         !nii
                ENDDO           !nij
              ENDDO             !nik
            ENDDO               !ne
          ENDDO                 !nr
          

C--   type---------------------------------------------------------------
        ELSE IF((TYPE(1:8).EQ.'GEOMETRY')
     '    .OR.(TYPE(1:6).EQ.'METRIC')) THEN
          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,
     '      ERROR,*9999)
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nxc=NXLIST(1)

          IF(CBBREV(CO,'DEFORMED',3,noco+1,NTCO,N3CO)) THEN
            DEFORMED=.TRUE.
            IF(CBBREV(CO,'FROM_CLASS',2,noco+1,NTCO,N3CO)) THEN
              nxc=IFROMC(CO(N3CO+1))
            ELSE
              nxc=1
            ENDIF
            IF(nxc.NE.0) THEN
              CALL NX_LOC(NX_INQUIRE,nxc,nx_d,NX_SOLVE,ERROR,*9999)
              CALL ASSERT(nx_d.NE.0,
     '          '>>No nx defined for this solve class',ERROR,*9999)
            ELSE
              nx_d=0
            ENDIF
          ELSE
            DEFORMED=.FALSE.
            nx_d=0 !not used
          ENDIF

        ELSE IF(TYPE(1:8).EQ.'SOLUTION') THEN
          IF(CBBREV(CO,'CLASS',6,noco+1,NTCO,N3CO)) THEN
            nxc=IFROMC(CO(N3CO+1))
          ELSE
            nxc=1
          ENDIF
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.NE.0,'>>Invalid solve class',ERROR,*9999)

        ELSE IF(TYPE(1:12).EQ.'CONNECTIVITY') THEN
          IF(CBBREV(CO,'IN_REGION',2,noco+1,NTCO,N3CO)) THEN
            nr1=IFROMC(CO(N3CO+1))
          ELSE
            nr1=1
          ENDIF
          IF(CBBREV(CO,'FROM_REGION',2,noco+1,NTCO,N3CO)) THEN
            nr2=IFROMC(CO(N3CO+1))
          ELSE
            nr2=2
          ENDIF

C MLB - 22-01-2003 - this is unnecessary, doesn't use regions/classes
C        ELSE IF(TYPE(1:8).EQ.'ADAPTIVE') THEN
C          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,
C     '      ERROR,*9999)
C          nx=1
C
C          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,
C     '      ERROR,*9999)

        ELSE !all other types
          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,
     '      ERROR,*9999)
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nxc=NXLIST(1)
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '       ERROR,*9999)

C news MPN 28Jun2000: added if statement: don't do this for extn ratio or green strain
          IF((TYPE(1:15).NE.'EXTENSION_RATIO').AND.
     '     (TYPE(1:12).NE.'GREEN_STRAIN')) THEN
            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
C SGM 26Oct2000 grid-based Finite element also
C MLT 29Nov02 Grid Finite Volume also
              CALL ASSERT(ITYP4(nr,nx).EQ.4.OR.ITYP4(nr,nx).EQ.6.OR.
     '                    ITYP4(nr,nx).EQ.7,
     '          '>>Must be collocation or grid-based FE/FV solution'
     '          ,ERROR,*9999)
            ENDDO !no_nrlist
          ENDIF
        ENDIF !type

        IF(TYPE(1:6).EQ.'METRIC') THEN
!         Check whether grid is adaptive
          IF(NMGT.GT.1) THEN !more than one grid level
            ADAPTIVE=.TRUE.
          ELSE               !fine grid level only
            ADAPTIVE=.FALSE.
          ENDIF
          IF(CBBREV(CO,'LEVEL',1,noco+1,NTCO,N3CO)) THEN
            ALL_LEVELS=.FALSE.
            na_level=IFROMC(CO(N3CO+1))
            IF(.NOT.ADAPTIVE) na_level=1 !since only na=1 is defined
          ELSE
            ALL_LEVELS=.TRUE.
          ENDIF
        ENDIF !metric

        IF(TYPE(1:6).EQ.'STRAIN') THEN
          IF(CBBREV(CO,'SAC_EXTENSION',5,noco+1,NTCO,N3CO)) THEN
            SAC_threshold=RFROMC(CO(N3CO+1))
          ELSE
            SAC_threshold=1.d0
          ENDIF
          WRITE(OP_STRING,'('' SAC_threshold='',F7.4)') SAC_threshold
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          IF(CBBREV(CO,'SAC_SCALE',5,noco+1,NTCO,N3CO)) THEN
            SAC_scale_factor=RFROMC(CO(N3CO+1))
          ELSE
            SAC_scale_factor=1.d0
          ENDIF
          WRITE(OP_STRING,'('' SAC_scale_factor='',F7.4)')
     '      SAC_scale_factor
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF !strain

C news MPN 28Jun2000: calc'ing values of extn ratio/green strain at grid pts
        IF((TYPE(1:15).EQ.'EXTENSION_RATIO').OR.
     '    (TYPE(1:12).EQ.'GREEN_STRAIN')) THEN
C news VJ 21Jan2004: Added element option to update grid extension_ratio/green_strain
          CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,
     '      NTCO,CO,ERROR,*9999) 
C newe VJ 21Jan2004
          IF(CBBREV(CO,'COMPONENT',3,noco+1,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),6,CMPTLIST(0),CMPTLIST(1),
     '        ERROR,*9999)
          ELSE
C            CMPTLIST(0)=1
C            CMPTLIST(1)=1
C VJ used to default to 1 component, but now assert that command must be specific
            CALL ASSERT(.FALSE.,'>>Command must be used to explicit.
     '        Specify atleast 1 component to update and specify array
     '        index to be updated for each component and for cell
     '        variant',ERROR,*9999)              
          ENDIF
C news VJ 29Jan2004: Only call Assert on num of components if
C         option is extension ratio. Don't know why it was done in the
C         first place. Restriction not required for GREEN_STRAIN option
          IF ((TYPE(1:15).EQ.'EXTENSION_RATIO')) THEN
            CALL ASSERT(CMPTLIST(0).EQ.1,
     '        '>>May only specify a single component',ERROR,*9999)
          ENDIF
C newe VJ 29Jan2004
          IF(CBBREV(CO,'YQS',3,noco+1,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),NQM,NQLIST(0),NQLIST(1),ERROR,*9999)
            DEFORMATION_YQS=.TRUE.
            DEFORMATION_RCQS=.FALSE.
          ELSEIF(CBBREV(CO,'RCQS',4,noco+1,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),NQM,NQLIST(0),NQLIST(1),ERROR,*9999)
            DEFORMATION_YQS=.FALSE.
            DEFORMATION_RCQS=.TRUE.
          ELSE
            CALL ASSERT(.FALSE.,'>>Must specify either YQS or RCQS',
     '        ERROR,*9999)
          ENDIF
         
C          IF(NQLIST(0).GT.1) CALL ASSERT(NQLIST(0).NE.CELL_NUM_VARIANTS,
C     '      '>>Must specify a index for each variant or a single index',
C     '      ERROR,*9999)
C news VJ 29Jan2004: Check that component list has numbers between
C          1 and 6
          DO nicmpt=1,CMPTLIST(0)
            CALL ASSERT(CMPTLIST(nicmpt).GE.1.AND.
     '        CMPTLIST(nicmpt).LE.6,
     '        '>>Extn ratio or green strain cmpt indices
     '        must be 1<=ncmpt<=6',ERROR,*9999)
          ENDDO !nicmpt
          IF(CBBREV(CO,'ALL_VARIANTS',12,noco+1,NTCO,N3CO)) THEN
C  copy the indices given for 1 variant to all other variants
            NQLIST(0)=CMPTLIST(0)*CELL_NUM_VARIANTS
            DO nqv=2,CELL_NUM_VARIANTS
              cellvarcmpt=CMPTLIST(0)*(nqv-1)
              DO nicmpt=1,CMPTLIST(0)
                NQLIST(nicmpt+cellvarcmpt)=NQLIST(nicmpt)
              ENDDO !nicmpt
            ENDDO !nqv
          ENDIF
C          Check NQLIST has equal components for each of the cell variants
C          For a pure passive mechanics problem, all cell variants will be used 
          CALL ASSERT(NQLIST(0).EQ.(CMPTLIST(0)*CELL_NUM_VARIANTS),
     '      '>>Must have the same number of array indices for each
     '      cell variant as the number of components being updated',
     '      ERROR,*9999)

C         DPN: If we are using RCQS array, need to map a RCQS index
C         to the appropriate row in the RCQS_SPATIAL array
          IF(DEFORMATION_RCQS) THEN
            DO nqv=1,CELL_NUM_VARIANTS
C             store intermediate NQLIST index calculation
              cellvarcmpt=CMPTLIST(0)*(nqv-1)
              DO nicmpt=1,CMPTLIST(0)
                FOUND=.FALSE.
                DO nqrsv=1,IRCQS_SPATIAL(0,nqv)
                  IF(NQLIST(nicmpt+cellvarcmpt).EQ.
     '                IRCQS_SPATIAL(nqrsv,nqv)) THEN
                    NQLIST(nicmpt+cellvarcmpt)=nqrsv
                    FOUND=.TRUE.
                  ENDIF
                ENDDO !nicmpt
              ENDDO !nqrsv
              CALL ASSERT(FOUND,
     '          '>>An extension ratio/green_strain index is not ,
     '          spatially varying',ERROR,*9999)
            ENDDO !nqv
          ENDIF
C news VJ 14Jan2004 If no_ze_calc chosen, interpolation from nodal
C                    to elemental variables will not be done
C                    It is assumed that this is done already
          IF(CBBREV(CO,'NO_ZE_CALC',3,noco+1,NTCO,N3CO)) THEN
            CALC_ZE = .FALSE.
          ELSE 
            CALC_ZE = .TRUE.
          ENDIF
C newe VJ 14Jan2004
C MPN 5Aug2014 adding ability to interpolate to an xi point
          IF(CBBREV(CO,'ATXI',4,noco+1,NTCO,N3CO)) THEN
            ATXI=1
            IF(CBBREV(CO,'XI_1',4,noco+1,NTCO,N3CO)) THEN
              XIPOS(1)=RFROMC(CO(N3CO+1))
            ELSE
              XIPOS(1)=0.5D0
            ENDIF
            IF(CBBREV(CO,'XI_2',4,noco+1,NTCO,N3CO)) THEN
              XIPOS(2)=RFROMC(CO(N3CO+1))
            ELSE
              XIPOS(2)=0.5D0
            ENDIF
            IF(CBBREV(CO,'XI_3',4,noco+1,NTCO,N3CO)) THEN
              XIPOS(3)=RFROMC(CO(N3CO+1))
            ELSE
              XIPOS(3)=0.5D0
            ENDIF
          ELSE
            ATXI=0
          ENDIF
        ENDIF !extension_ratio/green_strain
C newe MPN 28Jun2000

        IF(TYPE(1:6).EQ.'SOURCE') THEN
          NOQUES=0
          FILEIP=.FALSE.
          IOTYPE=1

          FORMAT='('' Enter function to use as source term [1]:'''//
     '      '/''   (1) -2(pi)^2 sin(pi x) sin(pi y)'''//
     '      '/''   (2) Unused'''//
     '      '/''   (3) Unused'''//
     '      '/''   (4) Unused'''//
     '      '/''   (5) Unused'''//
     '      '/$,''    '',I1)'
          IDEFLT(1)=1
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,1,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IDATA(1).GT.0) THEN
            IFNTYP=IDATA(1)
          ELSE
            IFNTYP=1
          ENDIF
        ENDIF !source/bidomain
C--type---------------------------------------------------------------

        IF(TYPE(1:8).EQ.'STIMULUS') THEN

          IF(CBBREV(CO,'DIFFUSION',3,noco+1,NTCO,N3CO)) THEN
            DIFFUSION=RFROMC(CO(N3CO+1))
          ELSE
            DIFFUSION=0.0d0
          ENDIF

          IF(CBBREV(CO,'DIRECTION',3,noco+1,NTCO,N3CO)) THEN
            DIRN=.TRUE.
            nicmpt=IFROMC(CO(N3CO+1))
          ELSE
            DIRN=.FALSE.
            nicmpt=0
          ENDIF

          CALL PARSE_GRID(NQLIST,noco,NTCO,CO,ERROR,*9999)

C ASL 21/Oct/2004 to map 2 lines together
          IF(CBBREV(CO,'MAP_GRID',3,noco+1,NTCO,N3CO)) THEN
C            Get the grid group after the command "map_grid and store it in 
C            NQLIST2"             
             CDATA(1)='GRIDS'
             CALL PARSILG(NQLIST2,NQM,CDATA(1),CO(N3CO+1),ERROR,*9999)

C            set the mapping grid status to be true           
             MAPG=.TRUE.         
          ELSE
             MAPG=.FALSE.
          ENDIF
C ASL 31/Jan/2005 to use potential gradient
          IF(CBBREV(CO,'UP_POTENTIAL',3,noco+1,NTCO,N3CO)) THEN
             UP_Vm=.TRUE.         
          ELSE
             UP_Vm=.FALSE.
          ENDIF

          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,
     '         ERROR,*9999)
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nx=NXLIST(1)

C         get the niqV which can be use for finding the index in YQ
          CALL NIQ_LOC(NIQ_INQUIRE, NIQ_ION, niqV,
     '                 NIQ_V,ERROR,*9999)

          IF(DOP) THEN
            WRITE(OP_STRING,'('' Diffusion coefficient '',F12.6)')
     '        DIFFUSION
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' at grid points '')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            CALL WRITE_LONG(INTTYPE,1,1,IODI,NQLIST(0),10,10,
     '        NQLIST(1),RDUMMY,'(1X,10I6)','(1X,10I6)',ERROR,*9999)
            IF(DIRN) THEN
              WRITE(OP_STRING,'('' Updating direction '',I4)') nicmpt
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ELSE
              WRITE(OP_STRING,'('' Updating all directions '')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF

C ASL 31/Jan/2005 initialise temp varible YQS which has the size of grid(1..NQT) to be zero
          IF(NQLIST(0).GT.0) THEN
            CALL ASSERT(NIQSM.GE.3,
     '              ' >> The NIQSM needs to be >= 3',ERROR,*9999)
C$OMP       PARALLEL DO
C$OMP&      PRIVATE(nq)
C$OMP&      SHARED(NQT,YQS)
            DO nq=1,NQT
              YQS(3,nq)=0.0d0
            ENDDO !nq
C$OMP       END PARALLEL DO

            NITB=NQXI(0,NQS(NEELEM(1,nr)))

C ASL 21/Oct/2004  Making one-one grid maping
            IF (MAPG) THEN
C              make sure two list has the same number
               CALL ASSERT(NQLIST(0).EQ.NQLIST2(0),
     '              ' >> The lists are the same size',ERROR,*9999)
                      
               DO nq=1,NQLIST(0)
C                 maping the YQ value from ICC(group1) to LM(group2)
                  nqI=NQLIST(nq)
                  nqL=NQLIST2(nq)
C                 the difference of voltage both layers * the coupling_coff give the flux moves from ICC to LM
C                 however, if coupling_coff=1,V(ICC)=10, V(LM)=-10, we want to make sure only half of
C                 differenceplus to the LM and half of the difference get substracted from ICC.
                  temp_var=DIFFUSION*
     '                 (YQ(nqL,niqV,1,nx)-YQ(nqI,niqV,1,nx))
                  YQ(nqI,niqV,1,nx)=YQ(nqI,niqV,1,nx)+temp_var
                  YQ(nqL,niqV,1,nx)=YQ(nqL,niqV,1,nx)-temp_var
               ENDDO !i

C ASL modify 31/Jan/2005 for Alievs' cell model, the intracellular domain was cut and reconnect by potential gradient
            ELSEIF (DIRN) THEN
               CALL ASSERT((nicmpt.LE.NITB).AND.(nicmpt.GT.0),
     &              '>>Invalid xi direction (1<=xi<=NIT)',ERROR,*9999)

C$OMP          PARALLEL DO
C$OMP&         PRIVATE(ni,nqI,nqL,nqq)
C$OMP&         SHARED(DIFFUSION,nicmpt,NQLIST,NXQ,YQS)
C              finding which direction is the LM layer by checking the sign of ni. Calculate the last term of Aliev equation
C              and store in the YQS array which later is added into YQ array. alpha=1 in LM, alpha=-1 in ICC
               DO nqq=1,NQLIST(0)
                  nqI=NQLIST(nqq)
                  ni=-nicmpt
                  IF(NXQ(ni,1,nqI,1).LT.0) THEN
                     nqL=ABS(NXQ(ni,1,nqI,1))
                     YQS(3,nqI)=YQS(3,nqI)+DIFFUSION*
     &                    (YQS(1,nqL)-YQS(1,nqI))
                     YQS(3,nqL)=YQS(3,nqL)-DIFFUSION*
     &                    (YQS(1,nqL)-YQS(1,nqI))
                  ENDIF
                  ni=nicmpt
                  IF(NXQ(ni,1,nqI,1).LT.0) THEN
                     nqL=ABS(NXQ(ni,1,nqI,1))
                     YQS(3,nqI)=YQS(3,nqI)+DIFFUSION*
     &                    (YQS(1,nqL)-YQS(1,nqI))
                     YQS(3,nqL)=YQS(3,nqL)-DIFFUSION*
     &                    (YQS(1,nqL)-YQS(1,nqI))
                  ENDIF
               ENDDO !nqq
C$OMP         END PARALLEL DO
            ELSE
C$OMP         PARALLEL DO
C$OMP&        PRIVATE(ni,NITB,nqI,nqL,nqq)
C$OMP&        SHARED(DIFFUSION,NQLIST,NXQ,YQS)
               DO nqq=1,NQLIST(0)
                  nqI=NQLIST(nqq)
                  DO ni=-NITB,NITB
                     IF(NXQ(ni,1,nqI,1).LT.0) THEN
                        nqL=ABS(NXQ(ni,1,nqI,1))
                        YQS(3,nqI)=YQS(3,nqI)+DIFFUSION*
     &                       (YQS(1,nqL)-YQS(1,nqI))
                        YQS(3,nqL)=YQS(3,nqL)-DIFFUSION*
     &                       (YQS(1,nqL)-YQS(1,nqI))
                     ENDIF
                  ENDDO !ni
               ENDDO !nqq
C$OMP         END PARALLEL DO
            ENDIF !dirn/mapg            

C ASL 31/Jan/2005 update YQ array from YQS, so the transmembrane potential info is writen in YQ array
            IF (UP_Vm) THEN
               DO nq=1,NQT
                  YQ(nq,niqV,1,nx)=YQ(nq,niqV,1,nx)+YQS(3,nq)
               ENDDO !nq
            ENDIF !UP_Vm
          ENDIF !NQLIST>0
      ENDIF !end stimulus

C--type---------------------------------------------------------------

        IF(TYPE(1:8).EQ.'MATERIAL') THEN
C *** Compute Dij, the diffusion matrix, at each grid point in terms of
C *** xi coordinates.  This becomes a full matrix generated from the
C *** diagonal diffusion tensor in material coordinates.
C *** Similarly for the bidomain model, compute the conductivity
C *** tensors.

          !Check for 2 classes if bidomain
          IF(KTYP32.EQ.2) THEN !bidomain
            CALL ASSERT(NXLIST(0).GE.2,
     '        ' >>Specify 2 classes for bidomain problems',ERROR,*9999)
            nxc=NXLIST(2)
            CALL NX_LOC(NX_INQUIRE,nxc,nx2,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx2.GT.0,'>>No nx defined for this solve class',
     '         ERROR,*9999)
          ELSE
            nx2=nx
          ENDIF

          !Initialise local arrays
          DO j=1,3
            DO i=1,3
              dNudXi(i,j)=0.0d0
              dXidNu(i,j)=0.0d0
            ENDDO
          ENDDO

          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            NITB=NQXI(0,NQS(NEELEM(1,nr)))

            IF(ITYP5(nr,nx).EQ.1.AND.ITYP2(nr,nx).EQ.5) THEN
              !Loop over all grid points in the current region

C$OMP         PARALLEL DO
C$OMP&        PRIVATE(DETM,dNudXi,dXidNu,nj,nq,i,j,k)
C$OMP&        SHARED(CQ,DNUDXQ,DXDXIQ,NITB,NQR,nr,nx,PROPQ)
              DO nq=NQR(1,nr),NQR(2,nr)
                DO k=1,4
                  DO j=1,3
                    DO i=1,3
                      !Initialise conductivity tensor arrays
                      PROPQ(i,j,k,1,nq,nx)=0.0d0
                      PROPQ(i,j,k,2,nq,nx)=0.0d0
                    ENDDO
                  ENDDO
                ENDDO

                !Compute dnu/dxi = dnu/dx * dx/dxi
                DO j=1,NITB
                  DO i=1,NITB
                    dNudXi(i,j)=0.0d0
                    DO nj=1,NITB
                      dNudXi(i,j)=dNudXi(i,j)+
     '                  DNUDXQ(i,nj,nq)*DXDXIQ(nj,j,nq)
                    ENDDO
                  ENDDO
                ENDDO

                !Compute dxi/dnu = 1/(dnu/dxi)
                CALL INVERT(NITB,dNudXi,dXidNu,DETM)

C *** For div(k.grad(u))=f the conductivities are in CQ(2..4,nq,nx)
                IF(NITB.EQ.3) THEN
                  DO j=1,3
                    DO i=1,3
                      PROPQ(i,j,1,1,nq,nx)=
     '                   CQ(2,nq,nx)*dXidNu(i,1)*dNudXi(1,j)
     '                  +CQ(3,nq,nx)*dXidNu(i,2)*dNudXi(2,j)
     '                  +CQ(4,nq,nx)*dXidNu(i,3)*dNudXi(3,j)
                    ENDDO !j
                  ENDDO !i
                ELSE
                  DO j=1,NITB
                    DO i=1,NITB
                      PROPQ(i,j,1,1,nq,nx)=
     '                   CQ(2,nq,nx)*dXidNu(i,1)*dNudXi(1,j)
     '                  +CQ(3,nq,nx)*dXidNu(i,2)*dNudXi(2,j)
                    ENDDO !j
                  ENDDO !i
                ENDIF
              ENDDO !nq
C$OMP         END PARALLEL DO

            ELSE IF(ITYP5(nr,nx).EQ.2.AND.ITYP2(nr,nx).EQ.9) THEN
              !Loop over all grid points in the current region

C$OMP         PARALLEL DO
C$OMP&        PRIVATE(DETM,dNudXi,dXidNu,i,j,k,nj,nq)
C$OMP&        SHARED(CQ,DNUDXQ,DXDXIQ,NITB,NQR,nr,nx,nx2,PROPQ)
              DO nq=NQR(1,nr),NQR(2,nr)
                DO k=1,4
                  DO j=1,3
                    DO i=1,3
                      !Initialise conductivity tensor arrays
                      PROPQ(i,j,k,1,nq,nx)=0.0d0
                      PROPQ(i,j,k,2,nq,nx)=0.0d0
                    ENDDO
                  ENDDO
                ENDDO

C                CALL INVERT(NITB,DNUDXQ(1,1,nq),dXidNu,DETM)

                !Compute dxi/dnu = 1/(dnu/dxi)
                IF(USE_LAT.EQ.0) THEN
                  !Compute dnu/dxi = dnu/dx * dx/dxi
                  DO j=1,NITB
                    DO i=1,NITB
                      dNudXi(i,j)=0.0d0
                      DO nj=1,NITB
                        dNudXi(i,j)=dNudXi(i,j)+
     '                    DNUDXQ(i,nj,nq)*DXDXIQ(nj,j,nq)
                      ENDDO
                    ENDDO
                  ENDDO
                  CALL INVERT(NITB,dNudXi,dXidNu,DETM)
                ELSE
                  CALL INVERT(NITB,DNUDXQ(1,1,nq),dXdNu,DETM)
                ENDIF

C *** For Bidomain
C ***   Compute and store C(I)ij(nq) and C(E)ij in PROPQ
C ***     (already initialised)
C ***   Intracellular conductivities in CQ(3;4;5,nq,nx)
C ***   Extracellular conductivities in CQ(6;7;8,nq,nx)

                IF(USE_LAT.EQ.0) THEN
                  DO i=1,NITB
                    DO j=1,NITB
                      DO k=1,NITB
                        PROPQ(i,j,1,1,nq,nx)=PROPQ(i,j,1,1,nq,nx)+
     '                    CQ(k+2,nq,nx)*dXidNu(i,k)*dNudXi(k,j)
                        PROPQ(i,j,1,2,nq,nx)=PROPQ(i,j,1,2,nq,nx)+
     '                    CQ(k+5,nq,nx)*dXidNu(i,k)*dNudXi(k,j)
                      ENDDO
                    ENDDO !j
                  ENDDO !i
                ELSE
                  DO i=1,NITB
                    DO j=1,NITB
                      DO k=1,NITB
                        PROPQ(i,j,1,1,nq,nx)=PROPQ(i,j,1,1,nq,nx)+
     '                    CQ(k+2,nq,nx)*dXdNu(i,k)*DNUDXQ(k,j,nq)
                        PROPQ(i,j,1,2,nq,nx)=PROPQ(i,j,1,2,nq,nx)+
     '                    CQ(k+5,nq,nx)*dXdNu(i,k)*DNUDXQ(k,j,nq)
                      ENDDO
                    ENDDO !j
                  ENDDO !i                 
                ENDIF

                CQ(6,nq,nx2)=CQ(6,nq,nx)
                CQ(7,nq,nx2)=CQ(7,nq,nx)
                CQ(8,nq,nx2)=CQ(8,nq,nx)

              ENDDO !nq
C$OMP         END PARALLEL DO

            ELSE IF(ITYP5(nr,nx).EQ.1.AND.ITYP2(nr,nx).EQ.3) THEN
              !Loop over all grid points in the current region

              DO nq=NQR(1,nr),NQR(2,nr)
                DO k=1,4
                  DO j=1,3
                    DO i=1,3
                      !Initialise conductivity tensor arrays
                      PROPQ(i,j,k,1,nq,nx)=0.0d0
                      PROPQ(i,j,k,2,nq,nx)=0.0d0
                    ENDDO
                  ENDDO
                ENDDO

                !Compute dnu/dxi = dnu/dx * dx/dxi
                DO j=1,NITB
                  DO i=1,NITB
                    dNudXi(i,j)=0.0d0
                    DO nj=1,NITB
                      dNudXi(i,j)=dNudXi(i,j)+
     '                  DNUDXQ(i,nj,nq)*DXDXIQ(nj,j,nq)
                    ENDDO
                  ENDDO
                ENDDO

                !Compute dxi/dnu = 1/(dnu/dxi)
                CALL INVERT(NITB,dNudXi,dXidNu,DETM)

                !material for standard laplace
                CQ(1,nq,nx)=1.0d0
                CQ(2,nq,nx)=1.0d0
                CQ(3,nq,nx)=1.0d0

                DO j=1,NITB
                  DO i=1,NITB
                    DO k=1,NITB
                      PROPQ(i,j,1,1,nq,nx)=PROPQ(i,j,1,1,nq,nx)+
     '                  CQ(k,nq,nx)*dXidNu(i,k)*dNudXi(k,j)
                    ENDDO
                  ENDDO !j
                ENDDO !i

                !extracellular equivalent
                CQ(6,nq,nx)=1.0d0
                CQ(7,nq,nx)=1.0d0
                CQ(8,nq,nx)=1.0d0

              ENDDO !nq
            ENDIF !ityp4

            IF(DOP) THEN
              WRITE(OP_STRING,'('' C(I)ij, C(E)ij:'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              DO nq=NQR(1,nr),NQR(2,nr)
                WRITE(OP_STRING,'('' nq: '',I5)') nq
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(9(F12.5))')
     '            ((PROPQ(i,j,1,1,nq,nx),j=1,nitb),i=1,nitb)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(9(F12.5))')
     '            ((PROPQ(i,j,1,2,nq,nx),j=1,nitb),i=1,nitb)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDDO !nq
            ENDIF

C *** Now compute Dij,k using a first order finite difference about
C *** each grid point.  This only needs to be computed for grid points
C *** that are internal (not on external bdy).
C *** Similarly for the bidomain model, the derivatives of the
C *** conductivity tensors.

C*** Possible data dependance here. Rewrite with k as the outer loop???
C$OMP       PARALLEL DO
C$OMP&      PRIVATE(i,j,k,l,nq,SIGMAMINUSE,SIGMAMINUSI,SIGMAPLUSE,
C$OMP&        SIGMAPLUSI)
C$OMP&      SHARED(NITB,NQR,nr,NWQ,nx,NXQ,PROPQ)
            DO nq=NQR(1,nr),NQR(2,nr)
              IF(NWQ(1,nq,1).EQ.0) THEN !internal g.p.
                DO k=1,NITB
                  DO j=1,NITB
                    DO i=1,NITB
                      SIGMAPLUSI=0.0d0
                      SIGMAMINUSI=0.0d0
                      SIGMAPLUSE=0.0d0
                      SIGMAMINUSE=0.0d0

                      IF(NXQ(k,0,nq,1).GT.0) THEN
                        DO l=1,NXQ(k,0,nq,1)
                          SIGMAPLUSI=SIGMAPLUSI+
     '                      PROPQ(i,j,1,1,NXQ(k,l,nq,1),nx)
                          SIGMAPLUSE=SIGMAPLUSE+
     '                      PROPQ(i,j,1,2,NXQ(k,l,nq,1),nx)
                        ENDDO !l
                        SIGMAPLUSI=SIGMAPLUSI/DBLE(NXQ(k,0,nq,1))
                        SIGMAPLUSE=SIGMAPLUSE/DBLE(NXQ(k,0,nq,1))
                      ENDIF ! IF(NXQ(k,0,nq,1).GT.0)

                      IF(NXQ(-k,0,nq,1).GT.0) THEN
                        DO l=1,NXQ(-k,0,nq,1)
                          SIGMAMINUSI=SIGMAMINUSI+
     '                      PROPQ(i,j,1,1,NXQ(-k,l,nq,1),nx)
                          SIGMAMINUSE=SIGMAMINUSE+
     '                      PROPQ(i,j,1,2,NXQ(-k,l,nq,1),nx)
                        ENDDO !l
                        SIGMAMINUSI=SIGMAMINUSI/DBLE(NXQ(-k,0,nq,1))
                        SIGMAMINUSE=SIGMAMINUSE/DBLE(NXQ(-k,0,nq,1))
                      ENDIF ! IF(NXQ(-k,0,nq,1).GT.0)

                      PROPQ(i,j,k+1,1,nq,nx)=SIGMAPLUSI-SIGMAMINUSI
                      PROPQ(i,j,k+1,2,nq,nx)=SIGMAPLUSE-SIGMAMINUSE

                    ENDDO !j
                  ENDDO !i
                ENDDO !k
              ENDIF !nwq
            ENDDO !nq
C$OMP       END PARALLEL DO

            IF(DOP) THEN
              WRITE(OP_STRING,'('' C(I)ij,k, C(E)ij,k:'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              DO nq=NQR(1,nr),NQR(2,nr)
                WRITE(OP_STRING,'('' nq: '',I5)') nq
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(9(F12.5))')
     '            (((PROPQ(i,j,k,1,nq,nx),k=2,nitb+1),
     '            j=1,nitb),i=1,nitb)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(9(F12.5))')
     '            (((PROPQ(i,j,k,2,nq,nx),k=2,nitb+1),
     '            j=1,nitb),i=1,nitb)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDDO
            ENDIF

          ENDDO !no_nrlist

          UP_GRID_MATERIAL=.TRUE.

C--type---------------------------------------------------------------

        ELSE IF(TYPE(1:8).EQ.'SOLUTION') THEN

          CALL ASSERT(CALL_GRID,' >>Define grid first',ERROR,*9999)
          CALL ASSERT(CALL_SOLV,' >>Define solve first',ERROR,*9999)
          CALL ASSERT(CALL_COUP,' >>Define couple first',ERROR,*9999)
          CALL ASSERT(KTYP90.EQ.8,' Invalid coupling type',ERROR,*9999)
          CALL ASSERT(CPLST(0,1).LE.1000,' >>Too many coupling points',
     '      ERROR,*9999)

          CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqV,NIQ_V,ERROR,*9999)

          DO i=1,CPLST(0,1)
            nq=CPLST(i,2)
            YQ(nq,niqV,1,nx)=(YQ(nq,niqV,1,nx)+
     '        YQ(CPLST(i,1),niqV,1,nx))/2.0d0
            DO ni=1,NQXI(0,NQS(NENQ(1,nq)))
              IF(NXQ(ni,1,nq,1).GT.0) THEN
                YQ(NXQ(ni,1,nq,1),niqV,1,nx)=
     '            (YQ(NXQ(ni,1,nq,1),niqV,1,nx)+
     '            YQ(CPLST(i,1),niqV,1,nx))/2.0d0
              ENDIF
              IF(NXQ(-ni,1,nq,1).GT.0) THEN
                YQ(NXQ(-ni,1,nq,1),niqV,1,nx)=
     '            (YQ(NXQ(-ni,1,nq,1),niqV,1,nx)+
     '            YQ(CPLST(i,1),niqV,1,nx))/2.0d0
              ENDIF
            ENDDO
          ENDDO

C--type---------------------------------------------------------------

        ELSE IF(TYPE(1:6).EQ.'STRAIN') THEN

          CALL ASSERT(NIQM.GE.7,'>>Increase NIQM (>=7)',ERROR,*9999)
          CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqV,NIQ_V,ERROR,*9999)
          CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqSAC,NIQ_SAC,ERROR,*9999)

          Esac=-20.d0 !(mV) reversal potential for SAC

C Generate inward Isac current from stretch if above threshold
C Do interior collocation pts first, then bdry pts

          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)

              NQ_row=9
              NG_row=3

              DO j=2,NQ_row-1    !loop over collocation points
                DO i=2,NQ_row-1
                  nqq=i+(j-1)*NQ_row
                  nq=NQNE(ne,nqq)
C                  nq=NQGE(i+(j-1)*NQ_row,ne,nb_extended) !is coll.n pt#

                  mLL=INT(REAL(i+1)/3.0)  !m(left   Gauss pt index)
                  nBB=INT(REAL(j+1)/3.0)  !n(bottom Gauss pt index)
                  mRR=min(mLL+1,NG_row)    !m(right  Gauss pt index)
                  nTT=min(nBB+1,NG_row)    !n(top    Gauss pt index)

                  ngBL=mLL+(nBB-1)*NG_row  !ng(bottom left  Gauss pt#)
                  ngBR=mRR+(nBB-1)*NG_row  !ng(bottom right Gauss pt#)
                  ngTL=mLL+(nTT-1)*NG_row  !ng(top left     Gauss pt#)
                  ngTR=mRR+(nTT-1)*NG_row  !ng(top right    Gauss pt#)

                  iL=3*mLL-1  !locates bottom left Gauss pt index
                  jB=3*nBB-1  ! in the collocation grid

                  Xi1=DBLE(i-iL)/DBLE(NG_row) !how far ij colloc.n pt
                  Xi2=DBLE(j-jB)/DBLE(NG_row) ! is betw adj Gauss pts

!                 Fibre extension computed from surrounding Gauss pts
                  Ext=(1.d0-Xi1)*(1.d0-Xi2)*FEXT(1,ngBL,ne)
     '               +      Xi1 *(1.d0-Xi2)*FEXT(1,ngBR,ne)
     '               +(1.d0-Xi1)*      Xi2 *FEXT(1,ngTL,ne)
     '               +      Xi1 *      Xi2 *FEXT(1,ngTR,ne)

C                  IF(DOP) THEN
CC$                  call mp_setlock()
C                    WRITE(*,'(/'' Xi1='',F5.2,'' Xi2='',F5.2)') Xi1,Xi2
C                    WRITE(*,'('' Exten. ratio at nq='',I3,'
C     '               //''' :'',F6.3)') nq,Ext
CC$                  call mp_unsetlock()
C                  ENDIF

                  IF(ITYP3(nr,nx).EQ.2) THEN      !FHN
                    IF(Ext.GT.SAC_threshold) THEN !SAC current generated
                      V=YQ(nq,niqV,1,nx) !membrane potential
                      Current=-(V-Esac)*CQ(18,nq,nx)
     '                        *(Ext-SAC_threshold)*SAC_scale_factor
                    ELSE
                      Current=0.d0
                    ENDIF
                    IF(Current.LT.0.0d0) THEN !above threshold
                      YQ(nq,niqSAC,1,nx)=Current
                    ELSE
                      YQ(nq,niqSAC,1,nx)=0.0d0
                    ENDIF
                  ELSE IF(ITYP3(nr,nx).EQ.3) THEN !VCD
                    IF(Ext.GT.SAC_threshold) THEN !SAC current generated
                      V=YQ(nq,niqV,1,nx) !membrane potential
                      Current=-(V-Esac)*CQ(13,nq,nx)
     '                  *(Ext-SAC_threshold)*SAC_scale_factor
                    ELSE
                      Current=0.d0
                    ENDIF
                    IF(Current.LT.0.0d0) THEN !above threshold
                      YQ(nq,niqSAC,1,nx)=Current
                    ELSE
                      YQ(nq,niqSAC,1,nx)=0.0d0
                    ENDIF
C!!! KAT: not merging from Oxford
C                  ELSE IF((ITYP3(nr,nx).EQ.10).AND.(KTYP33.EQ.5)) THEN !User_cell5 i.e. the ph dependance of Noble98
C                    ss_index=ss
C                    YQS(SS,nq)=Ext
C                    Current=0.d0
                  ENDIF !ityp3

                  IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP               CRITICAL(UPGRID_5)
                    WRITE(OP_STRING,'(''nq:'',I6,'' Current:'',F12.6)')
     '                nq,Current
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP               END CRITICAL(UPGRID_5)
                  ENDIF

                ENDDO !i
              ENDDO !j
            ENDDO !noelem

C         Boundary collocation pts
            IF ((ITYP3(nr,nx).EQ.3).OR.(ITYP3(nr,nx).EQ.2)) THEN
C For VCD and FHN models only, the SAC currents for the Noble98 equations
C are calculated as part of the cellular equations

              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)

C         RH bdry of element ne
                DO j=1,9
                  nqq=9*j
                  nq=NQNE(ne,nqq)
C                nq=NQGE(9*j,ne,nb_extended)
                  nq1=NXQ(-1,1,nq,1)
                  nq2=NXQ( 1,1,nq,1)
C                IF(DOP) THEN
CC$                call mp_setlock()
C                  WRITE(*,'(/'' nq='',I3,'' nq1='',I3,'' nq2='',I3)')
C     '              nq,nq1,nq2
CC$                call mp_unsetlock()
C                ENDIF
                  IF(nq2.GT.0) THEN !RH point exists
                    Current=0.5d0*(YQ(nq1,niqSAC,1,nx)+
     '                YQ(nq2,niqSAC,1,nx))

                    IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP               CRITICAL(UPGRID_6)
                     WRITE(OP_STRING,'(''nq:'',I6,'' Current:'',F12.6)')
     '                  nq,Current
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP               END CRITICAL(UPGRID_6)
                    ENDIF

                    IF(Current.LT.0.0d0) THEN !above threshold
                      YQ(nq,niqSAC,1,nx)=Current
                    ELSE
                      YQ(nq,niqSAC,1,nx)=0.d0
                    ENDIF
                  ENDIF
                ENDDO !j

C         TOP bdry of element ne
                DO i=1,9
                  nqq=72+i
                  nq=NQNE(ne,nqq)
C                nq =NQGE(72+i,ne,nb_extended)
                  nq1=NXQ(-2,1,nq,1)
                  nq2=NXQ( 2,1,nq,1)
C                IF(DOP) THEN
CC$                call mp_setlock()
C                  WRITE(*,'(/'' nq='',I3,'' nq1='',I3,'' nq2='',I3)')
C     '              nq,nq1,nq2
CC$                call mp_unsetlock()
C                ENDIF
                  IF(nq2.GT.0) THEN !TOP point exists
                    Current=0.5d0*(YQ(nq1,niqSAC,1,nx)+
     '                YQ(nq2,niqSAC,1,nx))

                    IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP               CRITICAL(UPGRID_7)
                     WRITE(OP_STRING,'(''nq:'',I6,'' Current:'',F12.6)')
     '                  nq,Current
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP               END CRITICAL(UPGRID_7)
                    ENDIF

                    IF(Current.LT.0.0d0) THEN !above threshold
                      YQ(nq,niqSAC,1,nx)=Current
                    ELSE
                      YQ(nq,niqSAC,1,nx)=0.d0
                    ENDIF
                  ENDIF
                ENDDO !i
              ENDDO !noelem
            ENDIF ! VCD and FHN models
          ENDDO !no_nrlist

C--type---------------------------------------------------------------

C news MPN 28Jun2000: calc'ing values of extn ratio at grid pts
C                    This may eventually be used to supercede the
C                    STRAIN option immediately above
C                    (still under debate!)

        ELSE IF((TYPE(1:15).EQ.'EXTENSION_RATIO').OR.
     '    (TYPE(1:12).EQ.'GREEN_STRAIN')) THEN
          
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
C news VJ: Added YPZP call. ZP is meant to be a temporary array which
C          is filled up using YPZP before using it.
            IF(CALC_ZE) THEN
              CALL YPZP(1,NBH,NEELEM,NHE(1,nxc),NHP(1,nr,nxc),
     &          NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nxc,NYNE,
     &          NYNP,YP(1,1,nxc),ZA,ZP,ERROR,*9999)
            ENDIF
C newe VJ
            CALL ASSERT(KTYP53(nr).GT.1,
     '        '>>Stresses must be referred to Nu in constitutive law',
     '        ERROR,*9999)

            RESET_GRID_PT(0)=0
            DO noelem=1,NELIST(0)
              ne=NELIST(noelem)
C news VJ 14Jan2004 If no_ze_calc chosen, interpolation from nodal
C                    to elemental variables will not be done
C                    It is assumed that this is done already
              IF(CALC_ZE) THEN
                CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '            NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '            SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
                CALL ZPZE(NBH(1,1,ne),nc,
     '            NHE(ne,nx),NKHE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '            nr,NVHE(1,1,1,ne),NW(ne,1,nx),nx,
     '            CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),
     '            ZE,ZP,ERROR,*9999)
              ENDIF
C newe VJ 14Jan2004
              nb=NBH(NH_LOC(1,nx),nc,ne)
              NITB=NIT(nb)
              SCHEME=NQS(ne)
              IF(ATXI.EQ.0) THEN !compute at grid pts
                II=MAX(1,NQXI(1,SCHEME))
                IJ=1
                IK=1
                IF(NQXI(0,SCHEME).GT.1) IJ=MAX(1,NQXI(2,SCHEME))
                IF(NQXI(0,SCHEME).GT.2) IK=MAX(1,NQXI(3,SCHEME))
                DO i=1,3
                 XI(i)=0.5d0 !initialise
                ENDDO !i
              ELSE !compute at specified xi point
                II=1
                IJ=1
                IK=1
                DO i=1,3
                 XI(i)=XIPOS(i) !initialise
                ENDDO !i
              ENDIF
              IF(KTYP3B.EQ.2) THEN ! Gauss point grid scheme
C VJ added assert to check that number of gauss points does not exceed
C AWG and AXIG dimension NQEM
                CALL ASSERT(NGT(nb).LE.NQEM,
     &            'Increase size of local arrays AWG and AXIG',
     &            ERROR,*9999)
                CALL GAUSS1(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '            NGAP(1,nb),PG(1,1,1,nb),AWG,AXIG,
     '            ERROR,*9999)
              ENDIF
C             Loop over the grid points in each element
              DO nik=1,IK
                DO nij=1,IJ
                  DO nii=1,II
                    neq=nii+((nij-1)*NQXI(1,SCHEME)) !local grid pt #
                    IF(NQXI(0,SCHEME).GT.1) neq=neq+((nik-1)*
     '                NQXI(1,SCHEME)*NQXI(2,SCHEME))
                    nq=NQNE(ne,neq)                  !global grid pt #
                    
C                   Evaluate each grid point only once
                    IF(NENQ(1,nq).EQ.ne) THEN
C                     Local Xi coords of grid point nq in element ne
                      IF(ATXI.EQ.0) THEN !compute at grid pts
                        IF(KTYP3B.EQ.2) THEN ! Gauss point grid scheme
                          XI(1) = AXIG(1,neq)
                          XI(2) = AXIG(2,neq)
                          XI(3) = AXIG(3,neq)
                        ELSE
                          IF(II.NE.1) XI(1)=DBLE(nii-1)/DBLE(II-1)
                          IF(IJ.NE.1) XI(2)=DBLE(nij-1)/DBLE(IJ-1)
                          IF(IK.NE.1) XI(3)=DBLE(nik-1)/DBLE(IK-1)
                        ENDIF
                      ELSE !compute at specified xi point
                        DO i=1,3
                          XI(i)=XIPOS(i)
                        ENDDO !i
                      ENDIF
                      IF(DOP) THEN
                        WRITE(OP_STRING,'(/''***********'')')
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                        WRITE(OP_STRING,'('' nq = '',I5)') nq
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                        WRITE(OP_STRING,'(''***********'')')
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                        WRITE(OP_STRING,'('' nii, nij, nik: '',3I5)')
     '                    nii,nij,nik
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                        WRITE(OP_STRING,'('' NQXI(1..,'',I2,'')'',3I5)')
     '                    SCHEME,
     '                    NQXI(1,SCHEME),NQXI(2,SCHEME),NQXI(3,SCHEME)
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                        WRITE(OP_STRING,'('' XI: '',3F10.4)')
     '                    XI(1),XI(2),XI(3)
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      ENDIF !dop
                      IF(KTYP56(nr).EQ.3.AND.KTYP55(nr).EQ.3) THEN
C                       Interpolate matl params at XI for pole-zero law
                        CALL CPXI(1,IBT,IDO,INP,NPNE(1,1,ne),nr,nx,
     '                    CE(1,ne,nx),CP(1,1,nx),CW,XI,ERROR,*9999)
                      ENDIF
C                     Interpolate midwall geom var.s XG and derivs wrt Xi
                      CALL XEXW(0,IBT,IDO,INP,NAN,NBJ(1,ne),nr,XE,XG,XI,
     '                  ERROR,*9999)
C                     Get derivs of Xi wrt undeformed Nu (body/fibre)
C                     coords, dXidNu
                      CALL DXIDXM(NITB,nr,dXidNu,DETM,XG,'Fibre',
     '                  ERROR,*9999)
C                     Interpolate dep var.s ZG and derivs wrt Nu (JP=1)
                      CALL ZEZW(0,1,IBT,IDO,INP,NAN,NBH(1,nc,ne),
     '                  NHE(ne,nx),nr,nx,dXidNu,ZE,ZG,XI,ERROR,*9999)
C                     Calculate deformed metric tensors wrt Nu (AZL,AZU)
                      CALL ZGMG(nb,nr,AZ,AZL,AZU,ZG,ERROR,*9999)
C news VJ 29Jan2004 NQLIST now contains indices of RCQS for every cell variant
C                   to update with components in CMPTLIST. Hence need to find
C                   what cell variant grid point uses and then move to the index
C                   of NQLIST corresponding to the component updated at grid points
C                   and the cell variant the grid point models
                      DO nicmpt=1,CMPTLIST(0)
                        ncmpt=CMPTLIST(nicmpt)
C identify cell variant used for grid point
                        CELL_VARIANT_NUM=ICQS_SPATIAL(1,nq)
                        cellvarcmpt=CMPTLIST(0)*(CELL_VARIANT_NUM-1)
                        niqs=NQLIST(nicmpt+cellvarcmpt)
                        DEFORMATION_VALUE = 0.0d0
C news VJ 8Dec2003: added option of calc'ing values of green strain at grid pts
                        IF(TYPE(1:15).EQ.'EXTENSION_RATIO') THEN
                          IF(ncmpt.LE.3) THEN !axial component
                            IF(KTYP3B.EQ.2) THEN !Gauss point grid scheme
C!!!                          This relies un YG being set up
C!!!                          not using the above calcs, using com file
C!!!                          to fill YG instead for now.
                              DEFORMATION_VALUE=DSQRT(2.0d0*YG(1,
     '                          neq,ne)+1.0d0)
                            ELSE
                              DEFORMATION_VALUE=DSQRT(AZL(ncmpt,
     '                          ncmpt))
                            ENDIF
                          ELSE IF(ncmpt.EQ.4) THEN
                            DEFORMATION_VALUE=DSQRT(AZL(1,2)) !check this?
                          ELSE IF(ncmpt.EQ.5) THEN
                            DEFORMATION_VALUE=DSQRT(AZL(1,3)) !check this?
                          ELSE IF(ncmpt.EQ.6) THEN
                            DEFORMATION_VALUE=DSQRT(AZL(2,3)) !check this?
                          ENDIF
                          IF(DABS(DEFORMATION_VALUE).LT.ZERO_TOL)
     '                      THEN
C                            extn ratio was not computed for some reason,
C                           so flag the grid pt for now so that it can be
C                           set to the value of a neighbouring pt later.
C                           eg. this occurs since fibre vectors
C                           cannot be computed at apex
C                           of a prolate spheroid etc.
                            RESET_GRID_PT(0)=RESET_GRID_PT(0)+1
                            IF(RESET_GRID_PT(0).LE.RESET_GRID_MAX)
     '                        RESET_GRID_PT(RESET_GRID_PT(0))=nq
                            IF(DOP) THEN
                              WRITE(OP_STRING,
     '                          '('' >>WARNING: Bad extension ratio'','
     '                          //''' at grid point nq='',I5)') nq
                              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                              WRITE(OP_STRING,
     '                          '(''            ... computed as '','
     '                          //'D12.5)') DEFORMATION_VALUE
                              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                            ENDIF !dop
                          ENDIF !EXT_R_V.LT.ZERO_TOL
                          IF(KTYP56(nr).EQ.3.AND.KTYP55(nr).EQ.3) THEN
C                           Check that extn ratio is below yield
C                           for pole-zero law
                            IF(ncmpt.LE.3) THEN       !axial component
                              pole1=CW((ncmpt-1)*3+2) !axial yield strain
                              maxextn=DSQRT(2.d0*pole1+1.d0)
                            ELSE                      !shear component
                              pole1=CW((ncmpt-1)*3+2) !shear yield strain 1
                              pole2=CW((ncmpt-1)*3+11)!shear yield strain 2
                              maxextn1=DSQRT(2.d0*pole1+1.d0)
                              maxextn2=DSQRT(2.d0*pole2+1.d0)
                              maxextn=MIN(maxextn1,maxextn2)
                            ENDIF
                            IF(DEFORMATION_VALUE.GT.(YIELD_TOL
     '                        *maxextn)) THEN
C                             extn ratio is beyond 90% of limit
                              IF(DOP) THEN
                                WRITE(OP_STRING,'('' >>WARNING: '
     '                            //'Extension ratio exceeded '','
     '                            //'I2,''% of max at grid point nq='','
     '                            //'I5)') INT(YIELD_TOL*100.d0),nq
                                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                                WRITE(OP_STRING,
     '                            '(''            computed as '',D12.5,'
     '                            //''' (yield extension= '',D12.5,'
     '                            //''') ... but reset to '',D12.5)')
     '                            maxextn,DEFORMATION_VALUE,
     '                            YIELD_TOL*maxextn
                                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                              ENDIF !dop
                              DEFORMATION_VALUE=YIELD_TOL*maxextn
                            ENDIF
                          ENDIF !KTYPEs for pole-zero law
C news VJ 9Dec2003 added an option to allow updates of green strain components 
                        ELSE IF(TYPE(1:12).EQ.'GREEN_STRAIN') THEN
                          IF(ncmpt.LE.3) THEN !axial component
                            DEFORMATION_VALUE=
     '                        0.5d0*(AZL(ncmpt,ncmpt)-1)
                          ELSE IF(ncmpt.EQ.4) THEN
                            DEFORMATION_VALUE=0.5d0*AZL(1,2)
                          ELSE IF(ncmpt.EQ.5) THEN
                            DEFORMATION_VALUE=0.5d0*AZL(1,3)
                          ELSE IF(ncmpt.EQ.6) THEN
                            DEFORMATION_VALUE=0.5d0*AZL(2,3)
                          ENDIF
                        ENDIF !TYPE
C newe VJ 9Dec2003         
                        IF(DEFORMATION_YQS) THEN
                          YQS(niqs,nq)=DEFORMATION_VALUE
                        ELSE
                          RCQS_SPATIAL(niqs,nq)=DEFORMATION_VALUE
                        ENDIF
                      ENDDO !nicmpt (ncmpt)
                      IF(DOP) THEN
                        WRITE(OP_STRING,'(/'' YQS(1..,nq): '',6E12.5)')
     '                    (YQS(NQLIST(nicmpt),nq),nicmpt=1,CMPTLIST(0))
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      ENDIF !dop
                    ENDIF !NENQ(1,nq).EQ.ne
                  ENDDO !nii
                ENDDO !nij
              ENDDO !nik
            ENDDO !noelem (ne)
            WRITE(CHAR1(1:3),'(I3)') RESET_GRID_PT(0)
            CALL ASSERT(RESET_GRID_PT(0).LE.RESET_GRID_MAX,
     '        '>>Increase dimension of RESET_GRID_PT array to '//CHAR1,
     '        ERROR,*9999)
C           For grid points that haven't been set, use value(s) at
C           a neighbouring grid point (or 1.0 if can't find one)
C                         eg. this occurs since fibre vectors
C                         cannot be computed at apex
C                         of a prolate spheroid etc.
            DO nqq=1,RESET_GRID_PT(0)
              nq=RESET_GRID_PT(nqq)
              FOUNDGP=.FALSE.
C             Search in all +/- Xi dirns for neighbour not in the list
              ni=-4
              DO WHILE(ni.LT.3.AND..NOT.FOUNDGP)
                ni=ni+1
                IF(ni.EQ.0) ni=ni+1
                IF(NXQ(ni,0,nq,1).NE.0) THEN !neighbour(s) exist
                  nq_adj=NXQ(ni,1,nq,1)
C                 make sure neighbour not in list of non-calc'ed g.p.s
                  IF(.NOT.INLIST(nq_adj,RESET_GRID_PT(1),
     '              RESET_GRID_PT(0),n1list)) FOUNDGP=.TRUE.
                ENDIF !NXQ(ni,0,nq,1).NE.0
              ENDDO !WHILE(ni.LE.3.AND..NOT.FOUNDGP)
C             Set YQS values
              DO nicmpt=1,CMPTLIST(0)
                niqs=NQLIST(nicmpt)
                IF(DEFORMATION_YQS) THEN
                  IF(FOUNDGP.AND.DABS(YQS(niqs,nq_adj)).GT.ZERO_TOL)
     '              THEN
                    YQS(niqs,nq)=YQS(niqs,nq_adj)
                  ELSE !there are no non-zero neighbours
                    YQS(niqs,nq)=1.d0
                  ENDIF !FOUNDGP
                  IF(DOP) THEN
                    WRITE(OP_STRING,
     '                '(''  Extn ratio at grid point nq='','
     '                //'I5,'' reset to YQS(niqs,nq)= '',D12.5)')
     '                nq,YQS(niqs,nq)
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ENDIF !dop
                ELSE
                  IF(FOUNDGP.AND.DABS(RCQS_SPATIAL(niqs,
     '              nq_adj)).GT.ZERO_TOL) THEN
                    RCQS_SPATIAL(niqs,nq)=RCQS_SPATIAL(niqs,nq_adj)
                  ELSE !there are no non-zero neighbours
                    RCQS_SPATIAL(niqs,nq)=1.d0
                  ENDIF !FOUNDGP
                  IF(DOP) THEN
                    WRITE(OP_STRING,
     '                '(''  Extn ratio at grid point nq='','
     '                //'I5,'' reset to RCQS_SPATIAL(niqs,nq)= '','
     '                //'D12.5)') nq,RCQS_SPATIAL(niqs,nq)
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ENDIF !dop
                ENDIF
              ENDDO !nicmpt
            ENDDO !nqq
          ENDDO !no_nrlist

C newe MPN 28Jun2000

C--type---------------------------------------------------------------

        ELSE IF(TYPE(1:10).EQ.'ACTIVTIMES') THEN
C       Update activation time into YQ(nq,2)
          CALL MAQ_LOC(MAQ_INQUIRE,MAQ_TIME,atime,MAQ_ACTIV_TIME,
     '      ERROR,*9999)
          DO nq=1,NQT
            YQ(nq,4,1,nx)=AQ(atime,nq)
          ENDDO !nq
C--type---------------------------------------------------------------

        ELSE IF(TYPE(1:4).EQ.'TIME') THEN
C       Update previous time step solution for multigrid transient heat eqtn
          DO nq=1,NQT
            YQ(nq,6,1,nx)=YQ(nq,1,1,nx)
          ENDDO !nq

C--type---------------------------------------------------------------

        ELSE IF(TYPE(1:12).EQ.'CONNECTIVITY') THEN

          IF(CBBREV(CO,'GRID',4,noco+1,NTCO,N3CO)) THEN
            CALL PARSE_GRID(NQLIST,noco,NTCO,CO,ERROR,*9999)

            IF(CBBREV(CO,'DIRECTION',4,noco+1,NTCO,N3CO)) THEN
              ni=IFROMC(CO(N3CO+1))
            ELSE
              ni=1
            ENDIF

            IF(DOP) THEN
              WRITE(OP_STRING,'('' Xi direction to break '',I3)') ni
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' at grid points '')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              CALL WRITE_LONG(INTTYPE,1,1,IODI,NQLIST(0),10,10,
     '          NQLIST(1),RDUMMY,'(1X,10I6)','(1X,10I6)',ERROR,*9999)
            ENDIF

            na=1
            DO nqq=1,NQLIST(0)
              nq=NQLIST(nqq)
              DO i=1,NXQ(ni,0,nq,na)
C MLT 27Nov02 Use absolute NXQ value to ensure nq2 remains positive
                nq2=ABS(NXQ(ni,i,nq,na))

                DO j=1,NXQ(-ni,0,nq2,na)
C MLT 27Nov02 Use absolute NXQ values to compare with positive nq and
C in assignments to ensure correct negative sign for broken link
                  IF(ABS(NXQ(-ni,j,nq2,na)).EQ.nq) THEN
                    NXQ(-ni,j,nq2,na)=-ABS(NXQ(-ni,j,nq2,na))
                  ENDIF
                ENDDO !j

                NXQ(ni,i,nq,na)=-ABS(NXQ(ni,i,nq,na))
              ENDDO !i
            ENDDO !nqq

          ELSE
            CALL UPGRID_CONN(NBJ,NEELEM,NKJE,NPF,NPNE,NQXI,
     '        nr1,nr2,NVJE,NXQ,NQS,NQNE,
     '        SE,XA,XE,XP,XQ,ERROR,*9999)
          ENDIF

C--type---------------------------------------------------------------

        ELSE IF(TYPE(1:8).EQ.'ADAPTIVE') THEN
          CALL ASSERT(CALL_GRID,' >>Define grid first',ERROR,*9999)
          CALL CalculateNLQ(NLQ,NWQ,NXQ,AQ,YQ(1,1,1,nx),ERROR,*9999)
          CALL ConstructNLQ(NAQ,NLQ,NWQ,NXQ,ERROR,*9999)

C--type---------------------------------------------------------------

        ELSE IF(TYPE(1:6).EQ.'SOURCE') THEN
C       Calculate source term for collocation grid
          DO nq=1,NQT
            GX=XQ(1,nq)
            GY=XQ(2,nq)
            IF(IFNTYP.EQ.1) THEN
              CQ(1,nq,nx)=-2.d0*PI*PI*DSIN(PI*GX)*DSIN(PI*GY)
            ENDIF !ifntyp
          ENDDO !nq

C--type---------------------------------------------------------------
C*** Martin Buist, June 1997

        ELSE IF(TYPE(1:8).EQ.'GEOMETRY') THEN
!         Calculate geometric coordinates of grid points
          CALL ASSERT(CALL_GRID,'>>No grid defined',ERROR,*9999)

C!!! DAH 25-JUL-2002 Update grid geometry for lattice grid
C    scheme. DXDIXQ2 and XQ are made here.        

          IF(USE_LAT.EQ.1) THEN
            CALL ASSERT(.NOT.DEFORMED,
     '        '>>Lattice grid scheme not implemented for deformation',
     '        ERROR,*9999)
!           Initialise NLQ
            Do nq=1,NQM
              NLQ(nq)=1 !for fine grid
            ENDDO !nq            
            CALL CALC_LATT_XQ(IBT,IDO,INP,NAN,NBJ,NENQ,NKJE,
     '        NPF,NPNE,nr,NVJE,DXDXIQ2,SE,XA,XE,XIQ,XP,
     '        XQ,ERROR,*9999)

          ELSE !USE_LAT==0
!         Initialise DXIX
            DO i=1,3
              DO j=1,3
                DXIX(j,i)=0.0d0
              ENDDO
            ENDDO
            
C           Initialise NLQ
            DO nq=1,NQM
              NLQ(nq)=1 !for fine grid
            ENDDO !nq
            
C           Store undeformed grid coordinates in AQ so the deformation
C           gradient tensor can be calculated for deformation problems
            IF(DEFORMED) THEN
              CALL MAQ_LOC(MAQ_INQUIRE,MAQ_COORD,maqx,MAQ_X_UNDEF,
     '          ERROR,*9999)
              IF(DOP) THEN
                IF(maqx.NE.0) THEN
                  WRITE(OP_STRING,'('' MAQ undef. X coord found '')')
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
              ENDIF
              IF(maqx.EQ.0) THEN
                CALL ASSERT(MAQ_LIST(0)+1.LE.NMAQM,'>>Increase NMAQM',
     '            ERROR,*9999)
                CALL MAQ_LOC(MAQ_ALLOCATE_AND_LOCK,MAQ_COORD,maqx,
     '            MAQ_X_UNDEF,ERROR,*9999)
                DO nq=1,NQM
                  AQ(maqx,nq)=XQ(1,nq)
                ENDDO !nq
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' MAQ undef. X coord locked '')')
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
              ENDIF
              IF(NJT.GE.2) THEN
                CALL MAQ_LOC(MAQ_INQUIRE,MAQ_COORD,maqy,MAQ_Y_UNDEF,
     '            ERROR,*9999)
                IF(DOP) THEN
                  IF(maqy.NE.0) THEN
                    WRITE(OP_STRING,'('' MAQ undef. Y coord found '')')
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDIF
                IF(maqy.EQ.0) THEN
                  CALL ASSERT(MAQ_LIST(0)+1.LE.NMAQM,'>>Increase NMAQM',
     '              ERROR,*9999)
                  CALL MAQ_LOC(MAQ_ALLOCATE_AND_LOCK,MAQ_COORD,maqy,
     '              MAQ_Y_UNDEF,ERROR,*9999)
                  DO nq=1,NQM
                    AQ(maqy,nq)=XQ(2,nq)
                  ENDDO !nq
                  IF(DOP) THEN
                    WRITE(OP_STRING,'('' MAQ undef. Y coord locked '')')
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDIF
              ENDIF
              IF(NJT.GE.3) THEN
                CALL MAQ_LOC(MAQ_INQUIRE,MAQ_COORD,maqz,MAQ_Z_UNDEF,
     '            ERROR,*9999)
                IF(DOP) THEN
                  IF(maqz.NE.0) THEN
                    WRITE(OP_STRING,'('' MAQ undef. Z coord found '')')
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDIF
                IF(maqz.EQ.0) THEN
                  CALL ASSERT(MAQ_LIST(0)+1.LE.NMAQM,'>>Increase NMAQM',
     '              ERROR,*9999)
                  CALL MAQ_LOC(MAQ_ALLOCATE_AND_LOCK,MAQ_COORD,maqz,
     '              MAQ_Z_UNDEF,ERROR,*9999)
                  DO nq=1,NQM
                    AQ(maqz,nq)=XQ(3,nq)
                  ENDDO !nq
                  IF(DOP) THEN
                    WRITE(OP_STRING,'('' MAQ undef. Z coord locked '')')
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
            
            DO nrr=1,NRLIST(0)
              nr=NRLIST(nrr)
              
              IF(DOP) THEN
                WRITE(OP_STRING,'(/'' NQR(1,'',I2,'')='',I8,'
     '            //''' NQR(2,'',I2,'')='',I8)') nr,NQR(1,nr),nr,NQR(2,
     '            nr)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF !DOP

C$OMP       PARALLEL DO
C$OMP&      PRIVATE(i,II,IJ,IK,IODI,IOER,nb,ne,nee,neq,nh,nhx,nii,nij,
C$OMP&        nik,nj,njj1,njj2,nq,nu,SCHEME,X_TEMP,XE,XG,XI,XQ_TEMP,ZE,
C$OMP&        ZG,ERROR)
C$OMP&      SHARED(CURVCORRECT,DEFORMED,DXDXIQ2,DXIX,IBT,IDO,INP,NAN,nc,
C$OMP&        NBH,NBJ,NEELEM,NENQ,NHE,NKHE,NKJE,NPF,NPNE,NQNE,NQS,NQXI,
C$OMP&        nr,nx_d,NVHE,NVJE,NW,SE,XA,XP,XQ,ZA,ZP,ERROR_FLAG)

              DO nee=1,NEELEM(0,nr)
                IF(.NOT.ERROR_FLAG) THEN
                  ne=NEELEM(nee,nr)
                  SCHEME=NQS(ne)
                  
C                 Map global parameters to local element parameters
                  CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '              NPNE(1,1,ne),nr,NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,
     '              ne),XE,XP,ERROR,*103)
                  IF(DEFORMED.AND.(nx_d.GT.0)) THEN
                    CALL ZPZE(NBH(1,1,ne),nc,NHE(ne,nx_d),NKHE(1,1,1,
     '                ne),NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),
     '                NW(ne,1,nx_d),nx_d,CURVCORRECT(1,1,1,ne),SE(1,1,
     '                ne),ZA(1,1,1,ne),ZE,ZP,ERROR,*103)
                  ENDIF
                  
                  II=MAX(1,NQXI(1,SCHEME))
                  IJ=1
                  IK=1
                  IF(NQXI(0,SCHEME).GT.1) IJ=MAX(1,NQXI(2,SCHEME))
                  IF(NQXI(0,SCHEME).GT.2) IK=MAX(1,NQXI(3,SCHEME))
                  DO i=1,3
                    X_TEMP(i)=0.0d0
                    XI(i)=0.0d0
                  ENDDO !i

C MPN JHC NOV 2004 adding grid at gauss for geometry
                  IF(KTYP3B.EQ.2) THEN ! Gauss point grid scheme
C VJ added assert to check that number of gauss points does not exceed
C AWG and AXIG dimension NQEM
                    nb=NBJ(NJ_LOC(NJL_GEOM,1,nr),ne)
                    CALL ASSERT(NGT(nb).LE.NQEM,
     &                'Increase size of local arrays AWG and AXIG',
     &                ERROR,*103)
                    CALL GAUSS1(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     &                nb,NGAP(1,nb),PG(1,1,1,nb),AWG,AXIG,ERROR,*103)
                  ENDIF
C end new NOV 2004

C                 Loop over the grid points in each element
                  DO nik=1,IK
                    DO nij=1,IJ
                      DO nii=1,II
                        neq=nii+((nij-1)*NQXI(1,SCHEME))
                        IF(NQXI(0,SCHEME).GT.1) neq=neq+((nik-1)*
     '                    NQXI(1,SCHEME)*NQXI(2,SCHEME))
                        nq=NQNE(ne,neq)
C                       Evaluate each grid point only once
                        IF(NENQ(1,nq).EQ.ne) THEN
C MPN JHC NOV 2004 adding grid at gauss for geometry
C                         Local xi coords of grid point nq in element ne
                          IF(KTYP3B.EQ.2) THEN ! Gauss point grid scheme
                            XI(1) = AXIG(1,neq)
                            XI(2) = AXIG(2,neq)
                            XI(3) = AXIG(3,neq)
                          ELSE
                            IF(II.NE.1) XI(1)=DBLE(nii-1)/DBLE(II-1)
                            IF(IJ.NE.1) XI(2)=DBLE(nij-1)/DBLE(IJ-1)
                            IF(IK.NE.1) XI(3)=DBLE(nik-1)/DBLE(IK-1)
                          ENDIF
C end new NOV 2004

C                         Use basis function to get geometric position
                          IF(DEFORMED.AND.(nx_d.GT.0)) THEN
                            CALL ZEZW(0,0,IBT,IDO,INP,NAN,NBH(1,1,ne),
     '                        NHE(ne,nx_d),nr,nx_d,DXIX,ZE,ZG,XI,
     '                        ERROR,*103)
                            DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
                              nj=NJ_LOC(NJL_GEOM,nhx,nr)
                              nh=NH_LOC(nhx,nx_d)
                              DO ni=0,NIT(NBH(nh,nc,ne))
                                nu=NU1(ni)
                                XG(nj,nu)=ZG(nhx,nu)
                              ENDDO !nu
                            ENDDO !nh/nj
C                           KAT 20Dec99: Need XG set up for fibres also
                            njj1=NJL_FIBR !fibres
                            DO njj2=1,NJ_LOC(njj1,0,nr)
                              nj=NJ_LOC(njj1,njj2,nr)
                              nb=NBJ(nj,ne)
                              XG(nj,1)=PFXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                          INP(1,1,nb),NAN(1,1,nb),nb,1,XE(1,nj),
     '                          XI)
                            ENDDO !nj-fibre
                          ELSE
                            CALL XEXW(0,IBT,IDO,INP,NAN,NBJ(1,ne),nr,XE,
     '                        XG,XI,ERROR,*103)
                          ENDIF
                          
C!!!                    KAT 11Oct00:
C                         DXDXIQ2 looks wrong if ITYP10(nr) != 1.
                          DO nj=1,NJT
                            X_TEMP(nj)=XG(nj,1)
                            DXDXIQ2(nj,1,nq)=XG(nj,2)
                            CALL ASSERT(NUM.GE.4,
     &                        '>>Increase NUM, must be >= to 4',ERROR,
     &                        *103)
                            DXDXIQ2(nj,2,nq)=XG(nj,4)
                            IF(NJT.EQ.3) THEN
                              CALL ASSERT(NUM.GE.7,
     &                          '>>Increase NUM, must be >= to 7',ERROR,
     &                          *103)
                              DXDXIQ2(nj,3,nq)=XG(nj,7)
                            ENDIF
                          ENDDO !nj
                          
C                         Change coords - all grid pts in rect.cart.
                          CALL XZ(ITYP10(nr),X_TEMP,XQ_TEMP)
                          DO nj=1,NJT
                            XQ(nj,nq)=XQ_TEMP(nj)
                          ENDDO

C                         Put fibre variables into XQ
                          njj1=NJL_FIBR !fibres
                          DO njj2=1,NJ_LOC(njj1,0,nr)
                            nj=NJ_LOC(njj1,njj2,nr)
                            XQ(nj,nq)=XG(nj,1)
                          ENDDO !nj-fibre
                        ENDIF
                      ENDDO !nii
                    ENDDO !nij
                  ENDDO !nik
                  
                  GO TO 104
C               This statement is designed to be skipped if no error
C               occurs. However if an error occurs within a subroutine
C               the alternate return points to line 100 to set the
C               flag
 103              CONTINUE

                  ERROR_FLAG=.TRUE.
                  IF(ERROR.NE.' ') THEN
                    CALL FLAG_ERROR(0,ERROR(:LEN_TRIM(ERROR)))
                  ENDIF
 104              CONTINUE            
                ENDIF !not error flag
              ENDDO !element
C$OMP       END PARALLEL DO

C KAT 8Sep00: Aborting on error
              IF(ERROR_FLAG) THEN
                ERROR=' '
                GOTO 9999
              ENDIF
            ENDDO !region

C         -------------------

C PM 27-11-01:
            IF(N_DISCRET.EQ.2) THEN ! arc_length based discretisation

C discretisation of 1D cubic Hermite elements with equal spacing on the
C arc-length is performed using the Trapezoidal rule.Note that arc-
C length in ipelem is calcualted with 4-point Gaussian quadrature. Since
C the error of numerical integration from the two methods is 
C inconsistent,
C the arc-length of elements is recalculated using the Trapezoidal rule.
C Setting N_DOSCRET(default value =1) to a different value is currently
C available only under 'CORONARY' option in IPGRID for 1D network.

              nq=1
              DO nrr=1,NRLIST(0)
                nr=NRLIST(nrr)
                
                DO nee=1,NEELEM(0,nr)
                  ne=NEELEM(nee,nr)
                  IF (NPL(1,1,NLL(1,ne)).EQ.4) THEN ! 1D Cubic Hermite
                    SCHEME=NQS(ne)
                    II=MAX(1,NQXI(1,SCHEME))
                    ndx=1
                    nb=NBJ(1,ne)
                    node1=NPNE(1,nb,ne)
                    node2=NPNE(2,nb,ne)
                    
                    S0=(XP(2,1,1,node1)**2+XP(2,1,2,node1)**2+
     '                XP(2,1,3,node1)**2)**0.5d0
                    S1=(XP(2,1,1,node2)**2+XP(2,1,2,node2)**2+
     '                XP(2,1,3,node2)**2)**0.5d0
                    
                    S=0.0d0
                    DO no_dxi=1,9999,1
                      Xi3=(no_dxi*1.0d0)/(10000.0d0) 
                      df1dxi=PH3(1,1,2,Xi3) !-6.0*xi+6.0*xi**2
                      df2dxi=PH3(1,2,2,Xi3) !1.0-4.0*xi+3.0*xi**2
                      df3dxi=PH3(2,1,2,Xi3) !6.0*xi-6.0*xi**2
                      df4dxi=PH3(2,2,2,Xi3) !-2.0*xi+3.0*xi**2
                      dxdxi=df1dxi*XP(1,1,1,node1)+
     '                  df2dxi*XP(2,1,1,node1)*DL(1,NLL(1,ne))+
     '                  df3dxi*XP(1,1,1,node2)+
     '                  df4dxi*XP(2,1,1,node2)*DL(2,NLL(1,ne))
                      dydxi=df1dxi*XP(1,1,2,node1)+
     '                  df2dxi*XP(2,1,2,node1)*DL(1,NLL(1,ne))+
     '                  df3dxi*XP(1,1,2,node2)+
     '                  df4dxi*XP(2,1,2,node2)*DL(2,NLL(1,ne))
                      dzdxi=df1dxi*XP(1,1,3,node1)+
     '                  df2dxi*XP(2,1,3,node1)*DL(1,NLL(1,ne))+
     '                  df3dxi*XP(1,1,3,node2)+
     '                  df4dxi*XP(2,1,3,node2)*DL(2,NLL(1,ne))
                      dS=(dxdxi**2+dydxi**2+dzdxi**2)**0.5d0
                      S=S+dS
                    ENDDO
                    S=(0.5d0*(S0+S1)+S)*(1.0d0/10000.0d0)
                    delx=S/((II-1)*1.0d0)
                    
                    IF(nq.EQ.NQNE(ne,1)) THEN
                      IF(NENQ(1,nq).EQ.ne) THEN
                        DO nj=1,NJT
                          XQ(nj,NQNE(ne,1))=XP(1,1,nj,node1)
                        ENDDO
                        nq=nq+1
                      ENDIF
                    ENDIF
                    
                    S0=(XP(2,1,1,node1)**2+XP(2,1,2,node1)**2+
     '                XP(2,1,3,node1)**2)**0.5d0
                    S=0.0d0
                    S1=0.0d0
                    
                    DO no_dxi=1,100000,1
                      IF (no_dxi.EQ.100000) THEN
                        DO nj=1,NJT
                          XQ(nj,NQNE(ne,II))=XP(1,1,nj,node2)
                        ENDDO
                        nq=nq+1
                      ELSE
                        Xi3=(no_dxi*1.0d0)/(100000.0d0)
                        df1dxi=PH3(1,1,2,Xi3) !-6.0*xi+6.0*xi**2
                        df2dxi=PH3(1,2,2,Xi3) !1.0-4.0*xi+3.0*xi**2
                        df3dxi=PH3(2,1,2,Xi3) !6.0*xi-6.0*xi**2
                        df4dxi=PH3(2,2,2,Xi3) !-2.0*xi+3.0*xi**2
                        dxdxi=df1dxi*XP(1,1,1,node1)+
     '                    df2dxi*XP(2,1,1,node1)*DL(1,NLL(1,ne))+
     '                    df3dxi*XP(1,1,1,node2)+
     '                    df4dxi*XP(2,1,1,node2)*DL(2,NLL(1,ne))
                        dydxi=df1dxi*XP(1,1,2,node1)+
     '                    df2dxi*XP(2,1,2,node1)*DL(1,NLL(1,ne))+
     '                    df3dxi*XP(1,1,2,node2)+
     '                    df4dxi*XP(2,1,2,node2)*DL(2,NLL(1,ne))
                        dzdxi=df1dxi*XP(1,1,3,node1)+
     '                    df2dxi*XP(2,1,3,node1)*DL(1,NLL(1,ne))+
     '                    df3dxi*XP(1,1,3,node2)+
     '                    df4dxi*XP(2,1,3,node2)*DL(2,NLL(1,ne))
                        dS=(dxdxi**2+dydxi**2+dzdxi**2)**0.5d0
                        S=S+(0.5d0*dS)
                        S1=(0.5d0*S0+S)*(1.0d0/100000.0d0)
                        
                        IF (ABS((S1-delx*ndx)).LE.0.001d0) THEN
                          
                          f1=PH3(1,1,1,Xi3) !1.0-3.0*xi**2+2.0*xi**3
                          f2=PH3(1,2,1,Xi3) !xi-2.0*xi**2+xi**3
                          f3=PH3(2,1,1,Xi3) !3.0*xi**2-2.0*xi**3
                          f4=PH3(2,2,1,Xi3) !-xi**2+xi**3
                          DO nj=1,NJT
                            XQ(nj,nq)=f1*XP(1,1,nj,node1)+
     '                        f2*XP(2,1,nj,node1)*DL(1,NLL(1,ne))+
     '                        f3*XP(1,1,nj,node2)+
     '                        f4*XP(2,1,nj,node2)*DL(2,NLL(1,ne))
                          ENDDO
                          nq=nq+1
                          ndx=ndx+1
                        ELSE
                          S=S+(0.5d0*dS)
                          GOTO 105
                        ENDIF
                      ENDIF
 105                ENDDO
                  ELSE
                    CALL ASSERT(NPL(1,1,NLL(1,ne)).NE.4,
     '                '>>Not 1D cubic Hermite elements',ERROR,*9999)
                  ENDIF !NPL
                ENDDO !element
              ENDDO !region
            ENDIF ! N_DISCRET

C         -------------------

            IF(DOP) THEN
              DO nq=1,NQT
                WRITE(OP_STRING,'(''XQ'',I6,3F10.4)') nq,XQ(1,nq),
     '            XQ(2,nq),XQ(3,nq)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDDO
            ENDIF
          ENDIF !USE_LAT
        
C--type---------------------------------------------------------------
C*** Martin Buist, June 1997, Updated June 1999

        ELSE IF(TYPE(1:8).EQ.'METRIC') THEN

C         Initialise X3G
          DO j=1,3
            DO i=1,4
              X3G(i,j)=0.0d0
            ENDDO !i
          ENDDO !j

          IF(DEFORMED) THEN
C           Check this sizing as the deformed code swaps XE,ZE
            CALL ASSERT(NHM.EQ.NJM,'>>Need NHM = NJM',ERROR,*9999)
C           Make sure this doesn't cause a problem later in ZEZW
            CALL ASSERT(NH_LOC(1,nx_d).LE.40,
     '        '>>Increase length of NBJ_LOCAL in UPGRID',ERROR,*9999)
          ENDIF

          DO nrr=1,NRLIST(0)
            nr=NRLIST(nrr)
            ICOOR_LOCAL=ITYP10(nr)
            ITYP10(nr)=1
          
!         Loop over grid points
C DAH 6-08-02 Metric terms are set up differently for lattice grid scheme
            IF(USE_LAT.EQ.1) THEN
              DO nq=NQR(1,nr),NQR(2,nr)
                IF(ALL_LEVELS) THEN !all grid levels
                  na=NLQ(nq)
                ELSE !specified level only
                  na=na_level
                ENDIF
                IF(.NOT.ERROR_FLAG.AND.na.GT.0.AND.na.LE.NAM) THEN
                  ne=NENQ(1,nq)
                  nb=NBJ(1,ne)
                  NITB=NIT(nb)
                  IF(DOP) THEN
                    WRITE(OP_STRING,'(/'' nq='',I6,'' na='',I2)') nq,na 
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF !DOP
                  CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '              NPNE(1,1,ne),nr,NVJE(1,1,1,ne),SE(1,1,ne),
     '              XA(1,1,ne),XE,XP,ERROR,*9999)
                  CALL XEXW(1,IBT,IDO,INP,NAN,NBJ(1,ne),nr,XE,
     '              XQD,XIQ(1,nq),ERROR,*9999)
                  DO ni=1,NITB
                    DO nj=1,NITB
                      DXDXIQ(nj,ni,nq)=XQD(nj,NU1(ni))
                    ENDDO !nj
                  ENDDO !ni                
                  
!                 Calc dnu/dx (direction cosines of material coords)
                  IF(NITB.EQ.1) THEN
 !it makes no sense to have a rotated fibre angle
 !in 1d as we are using a 'ds' sense
                    DO nij1=1,NJT
                      DO nij2=1,NJT
                        DNUDXQ(nij1,nij2,nq)=0.0d0
                      ENDDO !nij2
                    ENDDO !nij1
                    DNUDXQ(1,1,nq)=1.0d0
                  ELSE
                    IF(CALL_FIBR) THEN !fibres defined
C!!! Deformed not implemented yet for lattice method
C                     IF(DEFORMED.AND.(nx_d.GT.0)) THEN
C                     !Maybe need to do zpze?
C                     ng_temp=0 ! to compute at xi position 
C                     CALL MAT_VEC_DEF(IBT,IDO,INP,NAN,NBH(1,nc,ne),
C                     '              NBJ_LOCAL,ng_temp,NHE(ne,nx_d),nr,
C                     nx_d,
C                     '              A_VECTOR,B_VECTOR,C_VECTOR,PG,XE,
C                     XQD,
C                     XI,
C                     '              ZE,ZG,ERROR,*106)
C                     ELSE

C*** 31-AUG-02 This has not been tested yet with fibres                    
                      CALL MAT_VEC_XI(IBT,IDO,INP,NAN,NBJ(1,ne),nr,
     '                  A_VECTOR,B_VECTOR,C_VECTOR,XE,XQD,
     '                  XIQ(1,nq),.FALSE.,ERROR,*9999)
C                     ENDIF
                      DO nj=1,NJT
                        DNUDXQ(1,nj,nq)=A_VECTOR(nj)
                        IF(NJT.GT.1) DNUDXQ(2,nj,nq)=B_VECTOR(nj)
                        IF(NJT.GT.2) DNUDXQ(3,nj,nq)=C_VECTOR(nj)
                      ENDDO !nj
                    ELSE !no fibres defined
                      DO nij1=1,NJT
                        DO nij2=1,NJT
                          IF(nij1.EQ.nij2) THEN
                            DNUDXQ(nij1,nij2,nq)=1.0d0
                          ELSE
                            DNUDXQ(nij1,nij2,nq)=0.0d0
                          ENDIF !diag
                        ENDDO !nij2
                      ENDDO !nij1
                    ENDIF !call_fibr
                  ENDIF
                ENDIF
              ENDDO

              CALL MAQ_LOC(MAQ_ALLOCATE_AND_LOCK,MAQ_COORD,maqx,
     &          MAQ_NORMAL_X,ERROR,*9999)
              CALL MAQ_LOC(MAQ_ALLOCATE_AND_LOCK,MAQ_COORD,maqy,
     &          MAQ_NORMAL_Y,ERROR,*9999)
              CALL MAQ_LOC(MAQ_ALLOCATE_AND_LOCK,MAQ_COORD,maqz,
     &          MAQ_NORMAL_Z,ERROR,*9999)

              DO nq=NQR(1,nr),NQR(2,nr)
                IF(NWQ(1,nq,1).NE.0) THEN
                  CALL NORM_LATTICE(NJT,nq,NWQ(1,nq,1),DXDXIQ(1,1,nq),
     &              DXDXIQ2(1,1,nq),XQ_TEMP,ERROR,*9999)
                  AQ(maqx,nq)=XQ_TEMP(1)
                  AQ(maqy,nq)=XQ_TEMP(2)
                  AQ(maqz,nq)=XQ_TEMP(3)
                ELSE
                  AQ(maqx,nq)=0.0d0
                  AQ(maqy,nq)=0.0d0
                  AQ(maqz,nq)=0.0d0
                ENDIF
              ENDDO
            ELSE !not lattice grid scheme

C             Find/create a quadratic basis fn to use for interpolation
              QUADBASIS=.FALSE.
              DO nb=1,NBT
                NITB=NIT(nb)
                IF(NITB.EQ.NIT(NBJ(1,NEELEM(1,nr)))) THEN
                  DO ni=1,NITB
                    IF(IBT(1,ni,nb).EQ.1) THEN
                      IF(IBT(2,ni,nb).EQ.2) THEN
                        nbq=nb
                        QUADBASIS=.TRUE.
                      ENDIF !quadratic
                    ENDIF !lagrange
                  ENDDO !ni
                ENDIF !same # xi's
              ENDDO !nb

              IF(QUADBASIS) THEN !basis already defined
                DO nj=1,40
                  NBJ_LOCAL(nj)=nbq
                ENDDO
                NITB=NIT(nbq)
                CALL ASSERT(NJ_LOC(0,0,nr).LE.40,
     '            ' >>Increase NBJ_LOCAL in UPGRID (metric)',
     '            ERROR,*9999)
              ELSE !create basis
                CALL ASSERT(NBT +1.LE.NBM ,' >>Increase NBM' ,
     &            ERROR,*9999)
                CALL ASSERT(NBFT+1.LE.NBFM,' >>Increase NBFM',
     &            ERROR,*9999)

                IF(NIT(NBJ(1,NEELEM(1,nr))).EQ.1) THEN
                  WRITE(OP_STRING,'('' Adding a quadratic basis'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ELSE IF(NIT(NBJ(1,NEELEM(1,nr))).EQ.2) THEN
                  WRITE(OP_STRING,'('' Adding a biquadratic basis'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ELSE IF(NIT(NBJ(1,NEELEM(1,nr))).EQ.3) THEN
                  WRITE(OP_STRING,'('' Adding a triquadratic basis'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDIF

                nb=NBT+1
                NBT=NBT+1
                NBFT=NBFT+1

                NBC(nb)=1
                NBI(nb)=1
                NBASEF(nb,0)=1
                NBASEF(nb,1)=nb
                NBASEF(nb,2)=nb
                NFBASE(1,nb)=nb
                NFBASE(2,nb)=nb

                NAT(nb)=0 !aux. bases
C   *** DPN 01 December 1999 - This always screws over physiome problems!!
C   ***   with nb being greater the NBFM, so adding assert!!
                IF(nb.GT.NBFM) THEN
                  WRITE(ERROR,'(''>> Increase NBFM >= '',I12)') nb
                  GOTO 9999
                ENDIF
                NAN(1,1,nb)=0
                NIT(nb)=NIT(NBJ(1,NEELEM(1,nr))) !xi dirns
                NUT(nb)=NIT(nb)*NIT(nb)+2
                CALL ASSERT(NUT(nb).LE.NUM,' >>Increase NUM',
     &            ERROR,*9999)
                DO ni=1,NIT(nb)
                  IBT(1,ni,nb)=1 !lagrange
                  IBT(2,ni,nb)=2 !quadratic
                  NGAP(ni,nb)=3 !3 (x3x3) gauss pts
                ENDDO !ni
                NGT(nb)=NGAP(1,nb)*NIT(nb)
                CALL ASSERT(NGT(nb).LE.NGM,' >>Increase NGM',
     &            ERROR,*9999)
                NNT(nb)=3**NIT(nb)
                CALL ASSERT(NNT(nb).LE.NNM,' >>Increase NNM',
     &            ERROR,*9999)

                !INP
                DO ni=1,4
                  POSITION(ni)=1
                ENDDO !ni
                DO nn=1,NNT(nb)
                  ni=1
                  DO WHILE(POSITION(ni).GT.3)
                    POSITION(ni)=1
                    ni=ni+1
                    POSITION(ni)=POSITION(ni)+1
                  ENDDO !ni
                  DO ni=1,NIT(nb)
                    INP(nn,ni,nb)=POSITION(ni)
                  ENDDO !ni
                  POSITION(1)=POSITION(1)+1
                ENDDO !nn

                !IDO
                NST(nb)=0
                NKT(0,nb)=1
                DO nn=1,NNT(nb)
                  NKT(nn,nb)=1
                  DO ni=1,NIT(nb)
                    DO nk=1,NKT(nn,nb)
                      IDO(nk,nn,ni,nb)=1
                    ENDDO !nk
                  ENDDO !ni
                  IF(NKT(nn,nb).GT.NKT(0,nb)) NKT(0,nb)=NKT(nn,nb)
                  NST(nb)=NST(nb)+NKT(nn,nb)
                ENDDO !nn
                CALL ASSERT(NKT(0,nb).LE.NKM,' >>Increase NKM',
     '            ERROR,*9999)
                CALL ASSERT(NST(nb).LE.NSM,' >>Increase NSM',
     &            ERROR,*9999)
                DO nn=1,NNT(nb)
                  DO nk=1,NKT(nn,nb)
                    IDO1=IDO(nk,nn,1,nb)
                    IDO2=1
                    IDO3=1
                    IF(NIT(nb).GE.2) IDO2=IDO(nk,nn,2,nb)
                    IF(NIT(nb).EQ.3) IDO3=IDO(nk,nn,3,nb)
                    IF(IDO1.EQ.1.AND.IDO2.EQ.1.AND.IDO3.EQ.1)
     '                IDO(nk,nn,0,nb)=1
                    IF(IDO1.EQ.2.AND.IDO2.EQ.1.AND.IDO3.EQ.1)
     '                IDO(nk,nn,0,nb)=2
                    IF(IDO1.EQ.1.AND.IDO2.EQ.2.AND.IDO3.EQ.1)
     '                IDO(nk,nn,0,nb)=4
                    IF(IDO1.EQ.2.AND.IDO2.EQ.2.AND.IDO3.EQ.1)
     '                IDO(nk,nn,0,nb)=6
                    IF(IDO1.EQ.1.AND.IDO2.EQ.1.AND.IDO3.EQ.2)
     '                IDO(nk,nn,0,nb)=7
                    IF(IDO1.EQ.2.AND.IDO2.EQ.1.AND.IDO3.EQ.2)
     '                IDO(nk,nn,0,nb)=9
                    IF(IDO1.EQ.1.AND.IDO2.EQ.2.AND.IDO3.EQ.2)
     '                IDO(nk,nn,0,nb)=10
                    IF(IDO1.EQ.2.AND.IDO2.EQ.2.AND.IDO3.EQ.2)
     '                IDO(nk,nn,0,nb)=11
                  ENDDO !nk
                ENDDO !nn

                NITB=NIT(nb)
                DO nj=1,40
                  NBJ_LOCAL(nj)=nb
                ENDDO
                CALL ASSERT(NJ_LOC(0,0,nr).LE.40,
     '            ' >>Increase NBJ_LOCAL in UPGRID (metric)',
     '            ERROR,*9999)
              ENDIF !quad basis found/created

C$OMP         PARALLEL DO
C$OMP&        PRIVATE(A_VECTOR,B_VECTOR,CHTOFF,C_VECTOR,DBM,DETM,GL,GU,
C$OMP&          na,ne,ng_row,ni,nii2,nij1,nij2,nik2,nj,nq,NQLIST,SUM,XE,
C$OMP&          XG,XI,ZE,ZG,ERROR)
C$OMP&        SHARED(ALL_LEVELS,AQ,DXDXIQ,DNUDXQ,ERROR_FLAG,GCHQ,GUQ,
C$OMP&          IBT,IDO,INP,IODI,IOER,IOOP,na_level,NAN,NBJ,NBJ_LOCAL,
C$OMP&          NEELEM,NENQ,NHM,NITB,NLQ,NQGP,NU1,nr,NWQ,NXQ,PG,
C$OMP&          DEFORMED,OP_STRING,XQ,X3G,nx_d)

              DO nq=NQR(1,nr),NQR(2,nr)
                IF(ALL_LEVELS) THEN !all grid levels
                  na=NLQ(nq)
                ELSE                !specified level only
                  na=na_level
                ENDIF
                IF(.NOT.ERROR_FLAG.AND.na.GT.0.AND.na.LE.NAM) THEN
                  ne=NENQ(1,nq)
                  IF(DOP) THEN
                    WRITE(OP_STRING,'(/'' nq='',I6,'' na='',I2)') nq,na
                    CALL WRITES(IODI,OP_STRING,ERROR,*106)
                  ENDIF !DOP

                  IF(NWQ(1,nq,na).EQ.0) THEN !internal grid point
                    CALL XQXE(NBJ(1,NEELEM(1,nr)),NENQ,nq,NQGP,NQLIST,
     &                nr,NXQ(-NIM,0,0,na),XE,XQ,DEFORMED,ERROR,*106)
                  ELSE                       !external grid point
                    CALL XQXE(NBJ(1,NEELEM(1,nr)),NENQ,NWQ(1,nq,na),
     &                NQGP,NQLIST,nr,NXQ(-NIM,0,0,na),XE,XQ,DEFORMED,
     &                ERROR,*106)
                  ENDIF !internal/external

                  CALL CALC_GRID_XI(NITB,NWQ(1,0,na),nq,
     &              NXQ(-NIM,0,0,na),XI,ERROR,*106)
                  CALL XEXW(1,IBT,IDO,INP,NAN,NBJ_LOCAL,nr,XE,XG,XI,
     &              ERROR,*106)

C                 Gi(lower) = dx(nj)/dxi(ni) in direction nj
                  DO ni=1,NITB
                    DO nj=1,NITB
                      DXDXIQ(nj,ni,nq)=XG(nj,NU1(ni))
                    ENDDO !nj
                  ENDDO !ni

                  IF(DOP) THEN
                    DO ni=1,NITB
                      DO nj=1,NITB
                        WRITE(OP_STRING,'(''dxdxiq'',3I6,F12.6)')
     '                    nj,ni,nq,DXDXIQ(nj,ni,nq)
                        CALL WRITES(IODI,OP_STRING,ERROR,*106)
                      ENDDO !nj
                    ENDDO !ni
                  ENDIF !dop

C                 Calculate Gij(lower) = dx(nj)/dxi(nii)*dx(nj)/dxi(nij)
                  DO nij2=1,3
                    DO nii2=1,3
                      GL(nii2,nij2)=0.0d0
                      GU(nii2,nij2)=0.0d0
                    ENDDO !nii2
                  ENDDO !nij2
                  DO nij2=1,NITB
                    DO nii2=1,NITB
                      DO nj=1,NITB
                        GL(nii2,nij2)=GL(nii2,nij2)+
     '                    (DXDXIQ(nj,nii2,nq)*DXDXIQ(nj,nij2,nq))
                      ENDDO !nj
                    ENDDO !nii2
                  ENDDO !nij2

C                 Calculate Gij(upper) by inversion
                  CALL INVERT(NITB,GL,GU,DETM)

                  IF(DABS(DETM).LT.ZERO_TOL) THEN !cannot invert
                    DO nij2=1,NITB
                      DO nii2=1,NITB
                        GUQ(nii2,nij2,nq)=0.0d0
                      ENDDO !nii2
                      GCHQ(nij2,nq)=0.0d0
                    ENDDO !nij2
                    DO nij2=1,NJT
                      DO nij1=1,NJT
                        DNUDXQ(nij1,nij2,nq)=0.0d0
                      ENDDO !nij2
                    ENDDO !nij1
                  ELSE !can invert

C                   Write Gij(upper) at grid point nq into GUQ
                    DO nij2=1,NITB
                      DO nii2=1,NITB
                        GUQ(nii2,nij2,nq)=GU(nii2,nij2)
                      ENDDO !nij2
                    ENDDO !nii2

                    IF(DOP) THEN
                      DO nij2=1,NITB
                        DO nii2=1,NITB
                          WRITE(OP_STRING,'(''guq'',3I6,F12.6)')
     '                      nii2,nij2,nq,GUQ(nii2,nij2,nq)
                          CALL WRITES(IODI,OP_STRING,ERROR,*106)
                        ENDDO !nij2
                      ENDDO !nii2
                    ENDIF !dop

C                   Calc the Christoffel symbol (CHTOFF)
                    CALL TOFFEL(1,NBJ_LOCAL(1),nr,CHTOFF,DBM,GU,XG,X3G,
     '                .FALSE.,ERROR,*106)

C                   Calc GCHQ(k) = CHTOFF(k,i,j)*Gij(upper)
                    DO nik2=1,NITB
                      SUM=0.0d0
                      DO nij2=1,NITB
                        DO nii2=1,NITB
                          SUM=SUM+(CHTOFF(nik2,nii2,nij2)*GU(nii2,nij2))
                        ENDDO !nij2
                      ENDDO !nii2
                      GCHQ(nik2,nq)=SUM
                    ENDDO !nik2

                    IF(DOP) THEN
                      DO nik2=1,NITB
                        WRITE(OP_STRING,'('' gchq'',2I6,F12.6)')
     '                    nik2,nq,GCHQ(nik2,nq)
                        CALL WRITES(IODI,OP_STRING,ERROR,*106)
                      ENDDO !nik2
                    ENDIF !dop

C                   Calc dnu/dx (direction cosines of material coords)
                    IF(NITB.EQ.1) THEN
                      !it makes no sense to have a rotated fibre angle
                      !in 1d as we are using a 'ds' sense
                      DO nij2=1,NJT
                        DO nij1=1,NJT
                          DNUDXQ(nij1,nij2,nq)=0.0d0
                        ENDDO !nij2
                      ENDDO !nij1
                      DNUDXQ(1,1,nq)=1.0d0
                    ELSE
                      IF(CALL_FIBR) THEN !fibres defined
C                       MLB 22-01-2003
C                       New code to correctly calculate the deformed
C                       material axes from the undeformed and deformed
C                       grid point positions. This now parallels what is
C                       done in the non-grid mechanics code elsewhere.
                        IF(DEFORMED) THEN
C                         Calculate ZE(XE) for the undeformed coordinates
                          IF(NWQ(1,nq,na).EQ.0) THEN !internal grid point
                            CALL XQXE_REF(NENQ,NITB,nq,NQLIST,nr,
     '                        NXQ(-NIM,0,0,na),NHM,AQ,ZE,ERROR,*106)
                          ELSE !external grid point
                            CALL XQXE_REF(NENQ,NITB,NWQ(1,nq,na),NQLIST,
     '                        nr,NXQ(-NIM,0,0,na),NHM,AQ,ZE,ERROR,*106)
                          ENDIF !internal/external

                          ng_row=0 ! to compute at xi position
                          CALL MAT_VEC_DEF(IBT,IDO,INP,NAN,NBJ_LOCAL,
     '                      NBJ_LOCAL,ng_row,NJ_LOC(NJL_GEOM,0,nr),nr,
     '                      nx_d,A_VECTOR,B_VECTOR,C_VECTOR,PG,ZE,XG,XI,
     '                      XE,ZG,.FALSE.,ERROR,*106)
                        ELSE
                          CALL MAT_VEC_XI(IBT,IDO,INP,NAN,NBJ_LOCAL,nr,
     '                      A_VECTOR,B_VECTOR,C_VECTOR,XE,XG,XI,.FALSE.,
     '                      ERROR,*106)
                        ENDIF
                        DO nj=1,NJT
                          DNUDXQ(1,nj,nq)=A_VECTOR(nj)
                          IF(NJT.GT.1) DNUDXQ(2,nj,nq)=B_VECTOR(nj)
                          IF(NJT.GT.2) DNUDXQ(3,nj,nq)=C_VECTOR(nj)
                        ENDDO !nj
                      ELSE !no fibres defined
                        DO nij2=1,NJT
                          DO nij1=1,NJT
                            IF(nij1.EQ.nij2) THEN
                              DNUDXQ(nij1,nij2,nq)=1.0d0
                            ELSE
                              DNUDXQ(nij1,nij2,nq)=0.0d0
                            ENDIF !diag
                          ENDDO !nij2
                        ENDDO !nij1
                      ENDIF !call_fibr
                    ENDIF !1d

                    IF(DOP) THEN
                      DO nj=1,NJT
                        WRITE(OP_STRING,'(''dnu/dxq'',3F12.6)')
     '                    DNUDXQ(1,nj,nq),DNUDXQ(2,nj,nq),
     '                    DNUDXQ(3,nj,nq)
                        CALL WRITES(IODI,OP_STRING,ERROR,*106)
                      ENDDO !nj
                    ENDIF !dop
                  ENDIF !can invert

                  GO TO 108
C                 This statement is designed to be skipped if no error
C                 occur. However if a error occurs within a subroutine
C                 the alternate return points to line 100 to set the
C                 flag
 106              CONTINUE

                  ERROR_FLAG=.TRUE.
                  IF(ERROR.NE.' ') THEN
                    CALL FLAG_ERROR(0,ERROR(:LEN_TRIM(ERROR)))
                  ENDIF
 108              CONTINUE

                ENDIF !not error_flag
              ENDDO !nq
C$OMP         END PARALLEL DO

C   KAT 8Sep00: Aborting on error
              IF(ERROR_FLAG) THEN
                ERROR=' '
                GOTO 9999
              ENDIF
            ENDIF !use_lat            
            ITYP10(nr)=ICOOR_LOCAL
          ENDDO !nr

          UP_GRID_TENSOR=.TRUE.

        ENDIF !type
      ENDIF !question mark help

      CALL EXITS(ROUTINENAME)
      RETURN
 9999 CALL ERRORS(ROUTINENAME,ERROR)
      RET_ERROR=ERROR
      CALL EXITS(ROUTINENAME)
      RETURN 1
      END
