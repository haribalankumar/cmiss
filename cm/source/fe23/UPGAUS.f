      SUBROUTINE UPGAUS(IBT,IDO,ILPIN,
     '  INP,NAN,NBH,NBHF,NBJ,NBJF,NEELEM,NELIST_INCLUDE,
     '  NELIST_EXCLUDE,NFFACE,NGLIST,NHE,NHP,NKB,NKH,NKHE,NKJE,NNF,
     '  NPF,NPNE,NPNODE,NQET,NQLIST,NQNE,NQS,NQSCNB,NQXI,NRE,NRLIST,
     '  NSB,NVHE,NVHP,NVJE,NW,NXLIST,NYNE,NYNP,
     '  AQ,CE,CG,CGE,CP,CURVCORRECT,FEXT,PG,PGNQE,RGX,SE,WG,
     '  XA,XE,XG,XIG,XP,YG,YGF,YP,YQ,YQS,ZA,ZE,ZG,ZP,STRING,RET_ERROR,*)

C#### Subroutine: UPGAUS
C###  Description:
C###    <html><pre>UPGAUS updates YG(nj,ng,ne) array from ZG array,
C###         or YG(niyg,ng,ne) array from principal stresses,
C###         or YG(niyg,ng,ne) array from activation timing,
C###         or YG(niyg,ng,ne) array from activation potentials/recovery,
C###         or YG(niyg,ng,ne) array from material parameters.
C###         or YG(niyg,ng,ne) array from diffusion of phi_m
C###         or FEXT(ni,ng,ne) array from calcium at grid points
C###         or YG(niyg,ng,ne) array from cellular/grid point fields
C###         or YG(niyg,ng,ne) array from YG(niyg2,ng,ne)
C###    </pre></html>

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'cell02.cmn'
      INCLUDE 'cellml.cmn'
      INCLUDE 'error0.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'maqloc00.inc'
      INCLUDE 'nqloc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  ILPIN(NMM,NRM,NXM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),
     '  NBJ(NJM,NEM),NBJF(NJM,NFM),NEELEM(0:NE_R_M,0:NRM),
     '  NELIST_INCLUDE(0:NEM),NELIST_EXCLUDE(0:NEM),
     '  NFFACE(0:NF_R_M,NRM),NGLIST(0:NGM),NHE(NEM,NXM),
     '  NHP(NPM,0:NRM,NXM),NKB(2,2,2,NNM,NBFM),
     '  NKH(NHM,NPM,NCM,0:NRM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NNF(0:17,6,NBFM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),NQET(NQSCM),
     '  NQLIST(0:NQM),NQNE(NEQM,NQEM),NQS(NEQM),NQSCNB(NQSCM),
     '  NQXI(0:NIM,NQSCM),NRE(NEM),NRLIST(0:NRM),NSB(NKM,NNM,NBFM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM),NXLIST(0:NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 AQ(NMAQM,NQM),CE(NMM,NEM,NXM),CG(NMM,NGM),
     '  CGE(NMM,NGM,NEM,NXM),CP(NMM,NPM,NXM),
     '  CURVCORRECT(2,2,NNM,NEM),FEXT(NIFEXTM,NGM,NEM),
     '  PG(NSM,NUM,NGM,NBM),PGNQE(NGM,NQEM,NQSCM),
     '  RGX(*),SE(NSM,NBFM,NEM),WG(NGM,NBM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XG(NJM,NUM),XIG(NIM,NGM,NBM),XP(NKM,NVM,NJM,NPM),
     '  YG(NIYGM,NGM,NEM),YGF(NIYGFM,NGFM,NFM),YP(NYM,NIYM,NXM),
     '  YQ(NYQM,NIQM,NAM,NXM),YQS(NIQSM,NQM),ZA(NAM,NHM,NCM,NEM),
     '  ZE(NSM,NHM),ZG(NHM,NUM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER RET_ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER atime,DESTINATION_FIELD,i,IBEG,ICHAR,IEND,IFNTYP,Infarct,
     '  IFROMC,INFO,j,N1LIST,N3CO,
     '  nb,NBH_temp(12,1),NBJ_temp(12),nb_stress,nc,ne,ng,ng_calc,
     '  nh,nhx,ni,nig,niqq,niqs,niqV,NITB,nj,njj1,njj2,NM1,NM2,
     '  noelem,elemadded,no_nrlist,
     '  NOQUES,nqq,nr,NU1(0:3),nx,nxc,nx_s,nyg,SCHEME,SOURCE_FIELD,
     '  YQS_INDEX
c cpb 18/10/96 Adding dynamic allocation for parallel local arrays
      INTEGER*4 XE_PTR,XG_PTR,ZE_PTR,ZG_PTR
      REAL*8 ADdotAS,AD_VECTOR(3),ALFA,ASdotCS,ASdotFD,AS_VECTOR(3),
     '  BD_VECTOR(3),BETA,BSdotFD,BS_VECTOR(3),BXdotBS,BX_VECTOR(3),
     '  Ca_factor,constant_value,CD_VECTOR(3),CS_VECTOR(3),CXdotBS,
     &  CX_VECTOR(3),
     '  DACOS_MOD,DDOT,DXIX(3,3),DZDX(3,3),EG(3,3),Ed(3,3),Ev(3,3),
     '  ENERGY,FD_VECTOR(3),GAMA,GD_VECTOR(3),GL(3,3),GU(3,3),GX,GY,
     '  HD_VECTOR(3),PHI(3),PST(3),
     '  R(3,3),RFROMC,RGX2D,RGZ,RGZ2D,RI1,RI2,RI3,RM(3,3),
     '  SUM,TC(3,3),TG(3,3),TN(3,3),TNA,U(3,3),XI(3),YG_FACTOR,
     '  YQS_FACTOR,se_gauss,se_dev,se_vol,C(3,3),detC
      CHARACTER CHAR*1,COORDS*9,TYPE*22,TYPE2*16,STRESSTYPE*17,
     '  MANIPULATION*30
      CHARACTER ERROR*(ERRSTRLEN)
      LOGICAL ALL_REGIONS,CBBREV,COLLOCATION,CONSTANT,ERROR_FLAG,
     '  FILEIP,INLIST,OP_TIMING,
     '  PRINCIPAL
      DATA NU1/1,2,4,7/
C      INTEGER ne1,noelem1
!     Functions
      REAL*8 DET

      CALL ENTERS('UPGAUS',*9999)

C CPB 18/10/96 Intialise pointers so the error free will work properly

      XE_PTR=0
      XG_PTR=0
      ZE_PTR=0
      ZG_PTR=0
      ICHAR=999

      nb_stress=0

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        WRITE(CHAR(1:1),'(I1)') NJ_LOC(NJL_FIEL,1,1)

C---------------------------------------------------------------------

C#### Command: FEM update gauss (field/dependent/grad/pressure)
C###  Parameter:      <from (nh/NJ_INDEX)[0]>
C###   Specifies the array index
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description: Updates the gauss point array from the
C###   field,dependent,grad and pressure
C###

        OP_STRING(1)=STRING(1:IEND)
     '    //' (field/dependent/grad/pressure)'
        OP_STRING(2)=BLANK(1:15)
     '    //'<from (nh/NJ_INDEX)['//CHAR(1:1)//']>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update gauss stress
C###  Parameter:      <(fibre/wall/ref)[fibre]>
C###   Specify the coordinate system the stress tensor is calculated
C###    with respect to
C###  Parameter:      <(total/passive/active/total_no_hydro)[total]>
C###    Specifies whether to use total stress, or just passive or
C###    active stress components.
C new CS 9/1/2001 for Espen - Specify whether to drop the the hydrostatic terms
C###  Parameter:      <(components/principal)[components]>
C###     Specify whether principle or tensor components are updated
C###  Parameter:      <(Cauchy/Piola/Nominal)[Cauchy]>
C###     Specify  wiether Cauchy,Piola or Nominal stress are updated
C###  Parameter:      <basis (#/geom)[geom]>
C###      Specify the basis function number used to update the stress
C###      values
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Interpolates stresses at Gauss points and stores them
C###    in the YG array.  The 'fibre/wall/ref' option specifies the
C###    coordinate system to which stresses are referred.  The
C###    'components/principal' option is used if principal stresses are
C###    required.  'Cauchy/Piola/Nominal' is the type of the stress
C###    tensor.  The'basis' option is used to select which basis is
C###    used to interpolate thestress field.

        OP_STRING(1)=STRING(1:IEND)//' stress'
        OP_STRING(2)=BLANK(1:15)//'<(fibre/wall/ref)[fibre]>'
        OP_STRING(3)=BLANK(1:15)//'<(total/passive'
     '    //'/active/total_no_hydro)[total]>'
        OP_STRING(4)=BLANK(1:15)
     '    //'<(components/principal)[components]>'
        OP_STRING(5)=BLANK(1:15)
     '    //'<(Cauchy/Piola/Nominal)[Cauchy]>'
        OP_STRING(6)=BLANK(1:15)//'<basis (#/geom)[geom]>'
        OP_STRING(7)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(8)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update gauss strain
C###  Parameter:      <(fibre/wall/ref)[fibre]>
C###   Specify the coordinate system the strain tensor is calculated
C###    with respect to
C###  Parameter:      <(components/principal/extension_ratios/invariants)[components]>
C###     Specify wiether principle, tensor components,
C###     extension_ratios or invariants are updated
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Interpolates strains at Gauss points and stores them
C###    in the YG array.  The 'fibre/wall/ref' option specifies the
C###    coordinate system to which strains are referred.  The
C###    'components/principalextension_ratios/invariants' option
C###    for the different measures of deformation.

        OP_STRING(1)=STRING(1:IEND)//' strain'
        OP_STRING(2)=BLANK(1:15)//'<(fibre/wall/ref)[fibre]>'
        OP_STRING(3)=BLANK(1:15)//'<(components/principal/'
     '    //'extension_ratios/invariants)[components]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(5)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update gauss strain_energy
C###  Description:
C###    Computes strain energy at Gauss points and stores it
C###    in the YG array.
C###    yg[1] = total strain energy density
C###    yg[2] = deviatoric component of strain energy density
C###    yg[3] = volumetric component of strain energy density
C###    Note that total = deviatoric + volumetric

        OP_STRING(1)=STRING(1:IEND)//' strain_energy'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
C---------------------------------------------------------------------

C#### Command: FEM update gauss activation
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Parameter:      <basis #[1]>
C###    Specify the basis function used to interpolate at the
C###    gauss points
C###  Description: Update date the gauss point array from activation
C###   variaables
C###

        OP_STRING(1)=STRING(1:IEND)//' activation'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(4)=BLANK(1:15)//'<basis #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update gauss potential
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <basis #[1]>
C###    Specify the basis function used to interpolate at the
C###    gauss points
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description: Update date the gauss point array from potential
C###   variaables
C###

        OP_STRING(1)=STRING(1:IEND)//' potential'
        OP_STRING(2)=BLANK(1:15)//'<basis #[1]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update gauss material
C###  Parameter:      <from MATERIAL#[1]>
C###      Specify the material number
C###  Parameter:      <to GAUSS_ARRAY#[1]>
C###      Specify the gauss point array numbers
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description: Update date the gauss point array from material
C###   variaables
C###

        OP_STRING(1)=STRING(1:IEND)//' material'
        OP_STRING(2)=BLANK(1:15)//'<from MATERIAL#[1]>'
        OP_STRING(3)=BLANK(1:15)//'<to GAUSS_ARRAY#[1]>'
        OP_STRING(4)=BLANK(1:15)//'<region #s/all[1]>'
        OP_STRING(5)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update gauss source
C###  Parameter:      <basis #[1]>
C###    Specify the basis function used to interpolate at the
C###    gauss points
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description: Updates the source variables at gauss points
C###

        OP_STRING(1)=STRING(1:IEND)//' source'
        OP_STRING(2)=BLANK(1:15)//'<basis #[1]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update gauss calcium
C###  Parameter:      <basis #[1]>
C###    Specify the basis function used to interpolate at the
C###    gauss points
C###  Parameter:      <factor #[1.0]>
C###    Specify the  [Ca]i factor
C###  Parameter:      <infarct ELEMENT#[0]>
C###    Specify the infarcted element numbers
C###  Parameter:      <region (#s/all)[1>]
C###    Specify the region numbers to update.
C###  Parameter:      <from_class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Updates [Ca]i at Gauss pts from surrounding grid pts.
C###    Used for controlling active mechanics from electrical models.

        OP_STRING(1)=STRING(1:IEND)//' calcium'
        OP_STRING(2)=BLANK(1:15)//'<basis #[1]>'
        OP_STRING(3)=BLANK(1:15)//'<factor #[1.0]>'
        OP_STRING(4)=BLANK(1:15)//'<infarct ELEMENT#[0]>'
        OP_STRING(5)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(6)=BLANK(1:15)//'<from_class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C news MPN 1Jun2000: copying/calc'ing values of gridvars at Gauss pts
C                    This may eventually be used to supercede the
C                    CALCIUM option immediately above
C                    (still under debate!)

C#### Command: FEM update gauss gridvars
C###  Parameter:      <yqs #s[1]>
C###    Specifies the list of grid point variable indices.
C###  Parameter:      <yg #s[1]>
C###    Specifies the list of gauss point location indices.
C###  Parameter:      <exclude element (GROUP/#s)[all]>
C###    Specify the element numbers to exclude (eg. infarcted elements).
C###  Parameter:      <scale_factor (factor)[1.0]>
C###    Used to scale the contribution from all grid points to the gauss
C###    point values (eg. to gradually increment the tension value)
C###  Parameter:      <region (#s/all)[all]>
C###    Specify the regions to update.
C###  Description:
C###    Sets up the Gauss point array (YG) to contain (interpolated)
C###    copies of grid point variables (YQS).  A corresponding order of
C###    the specified the grid point and Gauss point indices is assumed.
C news VJ 28Jan2004: Added include element option to update gridvars
        OP_STRING(1)=STRING(1:IEND)//' gridvars'
        OP_STRING(2)=BLANK(1:15)//'<yqs #s[1]>'
        OP_STRING(3)=BLANK(1:15)//'<yg #s[1]>'
        OP_STRING(4)=BLANK(1:15)//'<include element (GROUP/#s)[all]>'
        OP_STRING(5)=BLANK(1:15)//'<exclude element (GROUP/#s)[all]>'
        OP_STRING(6)=BLANK(1:15)//'<scale_factor (factor)[1.0]>'
        OP_STRING(7)=BLANK(1:15)//'<region (#s/all)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C newe MPN 1Jun2000

C---------------------------------------------------------------------

C *** DPN 28 April 2001 - Adding a command to copy/manipulate one
C ***   field in YG to another

C#### Command: FEM update gauss gauss_field
C###  Parameter:      <destination_field #[1]>
C###    Specifies the destination Gauss point field
C###  Parameter:      <source_field #[2]>
C###    Specifies the source Gauss point field
C###  Parameter:      <substitute/difference [substitute]>
C###    The type of manipulation to perform
C###  Parameter:      <factor #[1.0]>
C###    A factor to be used by some of the above manipulations
C###  Parameter:      <region (#s/all)[all]>
C###    Specify the regions to update.
C###  Description:
C###    Updates the destination field by the given manipulation of
C###    the source field. The manipulations implemented are simple
C###    substitution and difference, where a given fraction
C###    (specified by the factor parameter) of the difference
C###    between the source and destination field is added to the
C###    destination field.

        OP_STRING(1)=STRING(1:IEND)//' gauss_field'
        OP_STRING(2)=BLANK(1:15)//'<destination_field #[1]>'
        OP_STRING(3)=BLANK(1:15)//'<source_field #[2]>'
        OP_STRING(4)= BLANK(1:15)
     '    //'<substitute/difference [substitute]>'
        OP_STRING(5)=BLANK(1:15)//'<factor #[1.0]>'
        OP_STRING(6)=BLANK(1:15)//'<region (#s/all)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update gauss deformed_fibres
C###  Parameter:      <collocation>
C###   Specify wiether the extended basis is used
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Calculates fibre angles appropriate for the deformed mesh
C###    and puts them into the Gauss point array YG. If the
C###    COLLOCATION option is chosen the extended basis is used
C###    so that deformed fibre angles are calculated at evenly
C###    spaced points.

        OP_STRING(1)=STRING(1:IEND)//' deformed_fibres'
        OP_STRING(2)=BLANK(1:15)//'<collocation>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update gauss eikonal
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Parameter:      <basis #[1]>
C###    Specify the basis function used to interpolate at the
C###    gauss points
C###  Description:
C###    Sets up the gauss point arrays from geometry and material
C###    parameters for solution of an eikonal equation.  If 'timing' is
C###    specified then CPU and elapsed times are output.

        OP_STRING(1)=STRING(1:IEND)//' eikonal'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(4)=BLANK(1:15)//'<timing>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc','UPGAUS',ERROR,*9999)
      ELSE
        IF(CBBREV(CO,'FIELD',1,noco+1,NTCO,N3CO)) THEN
          TYPE='FIELD'
          IF(CBBREV(CO,'CONSTANT',1,noco+1,NTCO,N3CO)) THEN
             CONSTANT=.TRUE.
             constant_value=RFROMC(CO(N3CO+1))
             IF(CBBREV(CO,'YG_FIELD',2,noco+1,NTCO,N3CO)) THEN
               nyg=IFROMC(CO(N3CO+1))
             ELSE
               nyg=1
             ENDIF
          ELSE
             CONSTANT=.FALSE.
          ENDIF
        ELSE IF(CBBREV(CO,'DEPENDENT',1,noco+1,NTCO,N3CO)) THEN
          TYPE='DEPENDENT'
        ELSE IF(CBBREV(CO,'GRADIENT',3,noco+1,NTCO,N3CO)) THEN
          TYPE='GRADIENT'
        ELSE IF(CBBREV(CO,'STRESS',4,noco+1,NTCO,N3CO)) THEN
          TYPE='STRESS'
          IF(CBBREV(CO,'PIOLA',3,noco+1,NTCO,N3CO)) THEN
            TYPE2='Piola'
          ELSE IF(CBBREV(CO,'NOMINAL',3,noco+1,NTCO,N3CO)) THEN
            TYPE2='Nominal'
          ELSE
            TYPE2='Cauchy'
          ENDIF
C news MPN 24May2000: total or just passive or active stress cmpts
          IF(CBBREV(CO,'TOTAL',3,noco+1,NTCO,N3CO)) THEN
            STRESSTYPE='Total'
          ELSE IF(CBBREV(CO,'TOTAL_NO_HYDRO',8,noco+1,NTCO,N3CO)) THEN
            STRESSTYPE='Total_no_hydro'
          ELSE IF(CBBREV(CO,'PASSIVE',3,noco+1,NTCO,N3CO)) THEN
            STRESSTYPE='Passive'
          ELSE IF(CBBREV(CO,'ACTIVE',3,noco+1,NTCO,N3CO)) THEN
            STRESSTYPE='Active'
          ELSE
            STRESSTYPE='Total'
          ENDIF
C newe MPN 24May2000:
C MPN 26Feb98: Changed nb to nb_stress to avoid confusion elsewhere
          IF(CBBREV(CO,'BASIS',2,noco+1,NTCO,N3CO)) THEN
            nb_stress=IFROMC(CO(N3CO+1)) ! otherwise it is zero
          ENDIF

        ELSE IF(CBBREV(CO,'STRAIN_ENERGY',8,noco+1,NTCO,N3CO)) THEN
          TYPE='STRAIN_ENERGY'
        ELSE IF(CBBREV(CO,'STRAIN',4,noco+1,NTCO,N3CO)) THEN
          TYPE='STRAIN'
          IF(CBBREV(CO,'EXTENSION_RATIOS',3,noco+1,NTCO,N3CO)) THEN
            TYPE2='EXTENSION_RATIOS'
          ELSE IF(CBBREV(CO,'INVARIANTS',3,noco+1,NTCO,N3CO)) THEN
            TYPE2='INVARIANTS'
          ELSE
            TYPE2='COMPONENTS'
          ENDIF
        ELSE IF(CBBREV(CO,'PRESSURE',2,noco+1,NTCO,N3CO)) THEN
          TYPE='PRESSURE'
        ELSE IF(CBBREV(CO,'ACTIVATION',1,noco+1,NTCO,N3CO)) THEN
          TYPE='ACTIVATION'
        ELSE IF(CBBREV(CO,'POTENTIAL',2,noco+1,NTCO,N3CO)) THEN
          TYPE='POTENTIAL'
        ELSE IF(CBBREV(CO,'MATERIAL',1,noco+1,NTCO,N3CO)) THEN
          TYPE='MATERIAL'
        ELSE IF(CBBREV(CO,'SOURCE',2,noco+1,NTCO,N3CO)) THEN
          TYPE='SOURCE'
        ELSE IF(CBBREV(CO,'CALCIUM',1,noco+1,NTCO,N3CO)) THEN
          TYPE='CALCIUM'
        ELSE IF(CBBREV(CO,'DEFORMED_FIBRES',2,noco+1,NTCO,N3CO)) THEN
          TYPE='DEFFIBRES'
        ELSE IF(CBBREV(CO,'EIKONAL',2,noco+1,NTCO,N3CO)) THEN
          TYPE='EIKONAL'
C news MPN 1Jun2000: copying/calc'ing values of gridvars at Gauss pts
        ELSE IF(CBBREV(CO,'GRIDVARS',4,noco+1,NTCO,N3CO)) THEN
          TYPE='GRIDVARS'
C newe MPN 1Jun2000
C *** DPN 28 April 2001 - Manipulation of YG fields
        ELSE IF(CBBREV(CO,'GAUSS_FIELD',5,noco+1,NTCO,N3CO)) THEN
          TYPE='GAUSS_FIELD'
        ELSE
          CO(noco+1)='?'
          GO TO 1
        ENDIF

!Initialise DXIX
        DO i=1,3
          DO j=1,3
            DXIX(i,j)=0.0d0
          ENDDO
        ENDDO

        IF(CBBREV(CO,'REFERENCE',3,noco+1,NTCO,N3CO)) THEN
          COORDS='Reference'
        ELSE IF(CBBREV(CO,'WALL',3,noco+1,NTCO,N3CO)) THEN
          COORDS='Wall'
        ELSE
          COORDS='Fibre'
        ENDIF

        IF(CBBREV(CO,'PRINCIPAL',1,noco+1,NTCO,N3CO)) THEN
          PRINCIPAL=.TRUE.
        ELSE
          PRINCIPAL=.FALSE.
        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

        IF(TYPE(1:7).EQ.'CALCIUM') THEN
          IF(CBBREV(CO,'FROM_CLASS',2,noco+1,NTCO,N3CO)) THEN
            nxc=IFROMC(CO(N3CO+1))
C     JBD 15/6/10 Added option of updating Ca2+ from YQS index
            IF(CBBREV(CO,'YQS_INDEX',2,noco+1,NTCO,N3CO)) THEN
              YQS_INDEX=IFROMC(CO(N3CO+1))
              CALL ASSERT(YQS_INDEX.GT.0,
     &             '>>Invalid YQS index',ERROR,*9999)
            ENDIF
          ELSE
            nxc=1
          ENDIF
          CALL NX_LOC(NX_INQUIRE,nxc,nx_s,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx_s.NE.0,'>>Invalid source class',ERROR,*9999)

          IF(CBBREV(CO,'FACTOR',2,noco+1,NTCO,N3CO)) THEN
            Ca_factor=RFROMC(CO(N3CO+1))
          ELSE
            Ca_factor=1.0d0
          ENDIF

          IF(CBBREV(CO,'INFARCT',3,noco+1,NTCO,N3CO)) THEN
            Infarct=IFROMC(CO(N3CO+1)) !element# for infarct (where [Ca]=0)
          ELSE
            Infarct=0
          ENDIF
        ELSE
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nxc=NXLIST(1)
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.NE.0,'>>Invalid class number',ERROR,*9999)
        ENDIF !type=calcium

C news MPN 1Jun2000: copying/calc'ing values of gridvars at Gauss pts
        IF(TYPE(1:8).EQ.'GRIDVARS') THEN
          IF(CBBREV(CO,'YQS',3,noco+1,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),NQM,NQLIST(0),NQLIST(1),ERROR,*9999)
          ELSE
            NQLIST(0)=1
            NQLIST(1)=1
          ENDIF
          IF(CBBREV(CO,'YG',2,noco+1,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),NGM,NGLIST(0),NGLIST(1),ERROR,*9999)
          ELSE
            NGLIST(0)=1
            NGLIST(1)=1
          ENDIF
          CALL ASSERT(NQLIST(0).EQ.NGLIST(0),
     '      '>>Must be the same number of YG and YQS indices',
     '      ERROR,*9999)

C *** DPN 28 March 2001 - Check the variable indices
          DO niqq=1,NQLIST(0)
            niqs=NQLIST(niqq)
            nig=NGLIST(niqq)
            CALL ASSERT(niqs.LE.NIQSM,'>>Need to increase NIQSM',
     '        ERROR,*9999)
            CALL ASSERT(nig.LE.NIYGM,'>>Need to increase NIYGM',
     '        ERROR,*9999)
          ENDDO


C news VJ 28Jan2004: Changed NELIST to have only included elem at all times
C          IF(CBBREV(CO,'EXCLUDE',2,noco+1,NTCO,N3CO)) THEN
C            CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
C     '        ERROR,*9999) !exclude elements in list
C          ELSE
C            NELIST(0)=0 !include all elements
C          ENDIF !exclude
C news VJ 17Feb2004 Changed code to have two arrays - NELIST_INCLUDE and NELIST_EXCLUDE
C need these arrays as all of YG components must be initialised to 0 for elements excluded and included
          IF(CBBREV(CO,'EXCLUDE',2,noco+1,NTCO,N3CO)) THEN
            CALL PARSE_ELEMENTS(NEELEM,NELIST_EXCLUDE,noco,NRLIST,
     '        NTCO,CO,ERROR,*9999) !include elements or exclude (depending on command) in list            
            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              NELIST_INCLUDE(0)=NEELEM(0,nr)-NELIST_EXCLUDE(0)
              elemadded=0
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                IF(.NOT.INLIST(ne,NELIST_EXCLUDE(1),
     &            NELIST_EXCLUDE(0),N1LIST)) THEN
                  elemadded=elemadded+1
                  NELIST_INCLUDE(elemadded)=ne
                ENDIF
              ENDDO !ne
            ENDDO !nr
          ELSE !include all or specific set of elements
            CALL PARSE_ELEMENTS(NEELEM,NELIST_INCLUDE,noco,NRLIST,
     '        NTCO,CO,ERROR,*9999) !include elements
            NELIST_EXCLUDE(0)=0 !exclude option not used          
          ENDIF !exclude
C newe VJ 
          IF(CBBREV(CO,'SCALE_FACTOR',2,noco+1,NTCO,N3CO)) THEN
            YQS_FACTOR = RFROMC(CO(N3CO+1))
          ELSE
            YQS_FACTOR = 1.0d0
          ENDIF

        ENDIF !type=gridvars
C newe MPN 1Jun2000

C *** DPN 28 April 2001 - Manipulation of YG fields
        IF(TYPE(1:11).EQ.'GAUSS_FIELD') THEN
          IF(CBBREV(CO,'DESTINATION_FIELD',3,noco+1,NTCO,N3CO)) THEN
            DESTINATION_FIELD = IFROMC(CO(N3CO+1))
          ELSE
            DESTINATION_FIELD = 1
          ENDIF
          IF(CBBREV(CO,'SOURCE_FIELD',3,noco+1,NTCO,N3CO)) THEN
            SOURCE_FIELD = IFROMC(CO(N3CO+1))
          ELSE
            SOURCE_FIELD = 2
          ENDIF
          !Check the field indices
          CALL ASSERT(DESTINATION_FIELD.LE.NIYGM,
     '      '>>Need to increase NIYGM',ERROR,*9999)
          CALL ASSERT(SOURCE_FIELD.LE.NIYGM,
     '      '>>Need to increase NIYGM',ERROR,*9999)

          IF(CBBREV(CO,'FACTOR',3,noco+1,NTCO,N3CO)) THEN
            YG_FACTOR = RFROMC(CO(N3CO+1))
          ELSE
            YG_FACTOR = 1.0d0
          ENDIF

          IF(CBBREV(CO,'SUBSTITUTE',4,noco+1,NTCO,N3CO)) THEN
            MANIPULATION = 'SUBSTITUTE'
          ELSE IF(CBBREV(CO,'DIFFERENCE',4,noco+1,NTCO,
     '        N3CO)) THEN
            MANIPULATION = 'DIFFERENCE'
          ELSE
            MANIPULATION = 'SUBSTITUTE'
          ENDIF

        ENDIF !type=gauss_field

        IF(TYPE(1:8).EQ.'MATERIAL') THEN
          IF(CBBREV(CO,'FROM',1,noco+1,NTCO,N3CO)) THEN
            NM1=IFROMC(CO(N3CO+1))
          ELSE
            NM1=1
          ENDIF
          IF(CBBREV(CO,'TO',1,noco+1,NTCO,N3CO)) THEN
            NM2=IFROMC(CO(N3CO+1))
          ELSE
            NM2=1
          ENDIF

        ELSE IF(TYPE(1:6).EQ.'SOURCE') THEN
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

        ELSE IF(TYPE(1:7).EQ.'EIKONAL') THEN
          CALL ASSERT(CALL_MATE,'>>ERROR: Define material first',
     '      ERROR,*9999)
          OP_TIMING=CBBREV(CO,'TIMING',2,noco+1,NTCO,N3CO)

        ELSE
          IF(CBBREV(CO,'FROM',1,noco+1,NTCO,N3CO)) THEN
            nh=IFROMC(CO(N3CO+1))
          ELSE
            nh=NJ_LOC(NJL_FIEL,1,NRLIST(1)) !default is field variable
          ENDIF
        ENDIF

        COLLOCATION=.FALSE.
        IF(TYPE(1:9).EQ.'DEFFIBRES') THEN
          IF(CBBREV(CO,'COLLOCATION',2,noco+1,NTCO,N3CO)) THEN
            COLLOCATION=.TRUE.

C old MLB 18 August 1997
C            nb_extended=1
C            DO WHILE (nb_extended.LE.NBT.AND.NBC(nb_extended).NE.7)
C              nb_extended=nb_extended+1
C            ENDDO
C            CALL ASSERT((NBC(nb_extended).EQ.7),
C     '        '>>Extended basis function not defined',ERROR,*9999)

          ENDIF
        ENDIF

C MLB 29 August 1996 Moved inside region loop
C        IF(TYPE(1:10).EQ.'ACTIVATION'.OR.TYPE(1:9).EQ.'POTENTIAL'.OR.
C     '    TYPE(1:8).EQ.'BIDOMAIN'.OR.TYPE(1:7).EQ.'CALCIUM') THEN
C          IF(CBBREV(CO,'BASIS',2,noco+1,NTCO,N3CO)) THEN
C            nb_extended=IFROMC(CO(N3CO+1))
C          ELSE
C            nb_extended=1
C            DO WHILE (nb_extended.LE.NBT.AND.NBC(nb_extended).NE.7)
C              nb_extended=nb_extended+1
C            ENDDO
C          ENDIF
C          CALL ASSERT((NBC(nb_extended).EQ.7),
C     '      '>>Extended basis function not defined',ERROR,*9999)
C        ENDIF

c cpb 20/3/95 Maybe the region loop should be down a level or two for
C efficiency

        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)

          IF(TYPE(1:9).EQ.'DEPENDENT'.OR.TYPE(1:8).EQ.'GRADIENT'
     '      .OR.TYPE(1:8).EQ.'PRESSURE'.OR.TYPE(1:9).EQ.'DEFFIBRES'
     '      .OR.TYPE(1:6).EQ.'STRESS'.OR.TYPE(1:6).EQ.'STRAIN') THEN
            CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),NKH(1,1,1,nr),
     '        NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,
     '        ERROR,*9999)
          ENDIF

          IF(TYPE(1:10).EQ.'ACTIVATION'.OR.TYPE(1:9).EQ.'POTENTIAL'.OR.
     '      TYPE(1:8).EQ.'BIDOMAIN'.OR.TYPE(1:7).EQ.'CALCIUM') THEN
            IF(CBBREV(CO,'BASIS',2,noco+1,NTCO,N3CO)) THEN
              nb=IFROMC(CO(N3CO+1))
            ELSE
              nb=NBJ(1,NEELEM(1,nr))
            ENDIF
C old MLB 18 August 1997
C              NXIELEM=NIT(NBJ(1,NEELEM(1,nr)))
C              !Find basis function using Extended Lagrange
C              nb_extended=1
C              DO nbb=1,NBT
C                IF((NXIELEM.NE.NIT(nb_extended)).OR.
C     '            (NBC(nb_extended).NE.7)) THEN
C                  nb_extended=nb_extended+1
C                ENDIF
C              ENDDO !nbb
C            ENDIF
C            CALL ASSERT((NBC(nb_extended).EQ.7),
C     '        '>>Extended basis function not defined',ERROR,*9999)
          ENDIF

          IF(TYPE(1:9).EQ.'DEPENDENT') THEN
            nc=1 !Temporary AJP 18-12-91
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              CALL ZPZE(NBH(1,1,ne),nc,NHE(ne,nx),NKHE(1,1,1,ne),
     '          NPF(1,1),NPNE(1,1,ne),NRE(ne),NVHE(1,1,1,ne),
     '          NW(ne,1,nx),nx,CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '          ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
              DO ng=1,NGT(NBH(nh,nc,ne))
                CALL ZEZG(1,NBH(1,nc,ne),ng,NHE(ne,nx),nx,DXIX,PG,
     '            ZE,ZG,ERROR,*9999)
                YG(nh,ng,ne)=ZG(nh,1)
              ENDDO
            ENDDO

          ELSE IF(TYPE(1:8).EQ.'GRADIENT') THEN
            nc=1 !Temporary AJP 18-12-91
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '          NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '          SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
              CALL ZPZE(NBH(1,1,ne),nc,NHE(ne,nx),NKHE(1,1,1,ne),
     '          NPF(1,1),
     '          NPNE(1,1,ne),NRE(ne),NVHE(1,1,1,ne),NW(ne,1,nx),
     '          nx,CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),ZE,
     '          ZP,ERROR,*9999)
              DO ng=1,NGT(NBH(nh,nc,ne))
                CALL XEXG(NBJ(1,ne),ng,nr,PG,
     '            XE,XG,ERROR,*9999)
                CALL XGMG(0,0,NBJ(1,ne),nr,DXIX,GL,GU,
     '            RGX(1),XG,ERROR,*9999)
                CALL ZEZG(1,NBH(1,nc,ne),ng,NHE(ne,nx),nx,DXIX,PG,
     '            ZE,ZG,ERROR,*9999)
                SUM=0.0d0
                DO ni=1,NIT(NBJ(1,ne))
                  SUM=SUM+ZG(nh,NU1(ni))**2.0d0
                ENDDO
                YG(nh,ng,ne)=DSQRT(SUM) !|grad z|
              ENDDO
            ENDDO

          ELSE IF(TYPE(1:5).EQ.'FIELD') THEN
            IF(CONSTANT)THEN
              CALL PARSE_ELEMENTS(NEELEM,NELIST_INCLUDE,noco,NRLIST,
     '          NTCO,CO,ERROR,*9999) !include elements
              DO noelem=1,NELIST_INCLUDE(0)
                ne=NELIST_INCLUDE(noelem)
                DO ng=1,NGT(NBJ(nh,ne))
                  YG(nyg,ng,ne)=constant_value
                ENDDO
              ENDDO
              
            ELSE
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '            NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '            SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
C MPN 26Feb98: avoid confusion with nb in other places
C old                nb=NBJ(nh,ne)
C old                 DO ng=1,NGT(nb)
                DO ng=1,NGT(NBJ(nh,ne))
                  CALL XEXG(NBJ(1,ne),ng,nr,PG,
     '              XE,XG,ERROR,*9999)
                  YG(1,ng,ne)=XG(nh,1)
                ENDDO
              ENDDO
            ENDIF

          ELSE IF(TYPE(1:6).EQ.'STRESS'.OR.TYPE(1:6).EQ.'STRAIN') THEN
            nc=1 !Temporary AJP 18-12-91

            IF(TYPE(1:6).EQ.'STRESS'.AND.
     '        KTYP54(nr).EQ.3) THEN !Gauss point stresses (grid coupling)

              WRITE(OP_STRING,'('' >>Warning: have you executed '
     '  //'[FEM update gauss gridvars yqs 2,3,4,5,6,7 yg 1,2,3,4,5,6]?'
     '  //''')')
              CALL WRITES(IOER,OP_STRING,ERROR,*9999)

C MPN 6Aug2014: Ideally the following 3 commands should be issued here -
C then the warning above can be removed.
C However it gives a seg fault (e.g. due to memory stomping, or
C recursion?
C FEM update grid green_strain comp 1,2,3,4,5,6 RCQS 3,6,8,4,5,7 ALL_VARIANTS
C FEM solve class 2 restart to 0
C FEM update gauss gridvars yqs 2,3,4,5,6,7 yg 1,2,3,4,5,6
C (See fe23/OPSTRE.f for function calls)

            ENDIF

            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)

C CS 10/7/97 Adding ability to choose nb
C MPN 26Feb98: moved into IF block below (no good for FE50 problems)
              IF(ITYP1(nr,nx).EQ.3) THEN
                ERROR='Not implemented for PDEs'
                GO TO 9999
              ELSE IF(ITYP1(nr,nx).EQ.4) THEN  !linear elasticity
C CS 10/7/97 Adding ability to choose nb
                DO njj1=1,3 !geom/fibres/field
                  DO njj2=1,NJ_LOC(njj1,0,nr)
                    nj=NJ_LOC(njj1,njj2,nr)
                    CALL ASSERT(nj.LE.12,'>>ERROR: increase size of '
     '                //'NBJ_temp to NJM',ERROR,*9999)
                    IF(nb_stress.EQ.0) THEN
                      NBJ_temp(nj)=NBJ(nj,ne)
                    ELSE
C                     NOTE: this isn't general, ie. cannot use
C                     different nb's for different geometric coords
C                     and/or dependent vars.
                      NBJ_temp(nj)=nb_stress
                    ENDIF
                  ENDDO !njj2
                ENDDO !njj1
                DO nhx=1,NHE(ne,nx)
                  nh=NH_LOC(nhx,nx)
                  CALL ASSERT(nh.LE.12,'>>ERROR: increase size of '
     '              //'NBH_temp to NHM',ERROR,*9999)
                  IF(nb_stress.EQ.0) THEN
                    NBH_temp(nh,1)=NBH(nh,1,ne)
                  ELSE
C                   NOTE: this isn't general, ie. cannot use
C                   different nb's for different geometric coords
C                   and/or dependent vars.
                    NBH_temp(nh,1)=nb_stress
                  ENDIF
                ENDDO !nhx
                IF(nb_stress.EQ.0) nb_stress=NBH(NH_LOC(1,nx),1,ne)
C old               nb_stress=NBH(NH_LOC(1,nx),1,ne)
                IF(TYPE(1:6).EQ.'STRAIN') THEN
                  ERROR='Not implemented for Linear Elastic Strains'
                  GO TO 9999
                ENDIF
                CALL XPXE(NBJ_temp,NKJE(1,1,1,ne),NPF(1,1),
     '            NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '            SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
                CALL ZPZE(NBH_temp,nc,
     '            NHE(ne,nx),NKHE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '            NRE(ne),NVHE(1,1,1,ne),NW(ne,1,nx),nx,
     '            CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),
     '            ZE,ZP,ERROR,*9999)
                CALL CPCG(NW(ne,1,nx),nb_stress,NPNE(1,1,ne),nr,nx,
     '            CE(1,ne,nx),CG,CGE(1,1,ne,nx),
     '            CP(1,1,nx),PG,ERROR,*9999)
                CALL OPST40(nb_stress,NBH_temp(1,1),NBJ_temp,ne,
     '            NHE(ne,nx),NW(ne,1,nx),nx,CG,ENERGY,PG,TYPE,
     '            WG,XE,XG,YG(1,1,ne),ZE,ZG,PRINCIPAL,COORDS,
     '            ERROR,*9999)
              ELSE IF(ITYP1(nr,nx).EQ.5) THEN !finite elasticity
C new MPN 26Feb98: Changing the basis as for FE40 causes different
C interpolations of the geometric/dependent vars and hence different
C (ie. incorrect) stress/strain calcs.  Stresses are computed at
C different Gauss pts by passing their Xi coords (XIG) to ZETX50/ZEEX50.
                CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '            NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '            SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
                CALL ZPZE(NBH(1,1,ne),nc,
     '            NHE(ne,nx),NKHE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '            NRE(ne),NVHE(1,1,1,ne),NW(ne,1,nx),nx,
     '            CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),
     '            ZE,ZP,ERROR,*9999)
                IF(nb_stress.EQ.0) nb_stress=NBH(NH_LOC(1,nx),1,ne)
C end new
                CALL CPCG(1,nb_stress,NPNE(1,1,ne),nr,nx,
     '            CE(1,ne,nx),CG,CGE(1,1,ne,nx),
     '            CP(1,1,nx),PG,ERROR,*9999)
                DO ng=1,NGT(nb_stress)
                  IF(TYPE(1:6).EQ.'STRESS') THEN

                    CALL ZETX50(COORDS,TYPE2,STRESSTYPE,IBT,IDO,INP,NAN,
     '                NBH(1,1,ne),NBJ(1,ne),ng,NHE(ne,nx),
     '                NPNE(1,1,ne),nr,ne,nx,
     '                CE(1,ne,nx),CG,CP(1,1,nx),FEXT(1,ng,ne),
     '                PG,PHI,PST,
     '                RGX(ng),RGX2D,RGZ,RGZ2D,
     '                RM,TC,TG,TN,TNA,XE,XG,XIG(1,ng,nb_stress),
     '                YG(1,ng,ne),ZE,ZG,ERROR,*9999)
                    IF(PRINCIPAL) THEN
                      YG(1,ng,ne)=PST(1) !max principal stress
                      YG(2,ng,ne)=PST(2) !2nd principal stress
                      YG(3,ng,ne)=PST(3) !min principal stress
                      YG(4,ng,ne)=PHI(1) !1st principal angle
                      YG(5,ng,ne)=PHI(2) !2nd principal angle
                      YG(6,ng,ne)=PHI(3) !3rd principal angle
                    ELSE
                      IF(TYPE2(1:6).EQ.'Cauchy') THEN
                        CALL ASSERT(NIYGM.GE.6,
     &                    '>>Increase NIYGM to at least 6',ERROR,*9999)
                        YG(1,ng,ne)=TC(1,1)
                        YG(2,ng,ne)=TC(2,2)
                        YG(3,ng,ne)=TC(3,3)
                        YG(4,ng,ne)=TC(1,2)
                        YG(5,ng,ne)=TC(1,3)
                        YG(6,ng,ne)=TC(2,3)
                      ELSE IF(TYPE2(1:5).EQ.'Piola') THEN
                        YG(1,ng,ne)=TG(1,1)
                        YG(2,ng,ne)=TG(2,2)
                        YG(3,ng,ne)=TG(3,3)
                        YG(4,ng,ne)=TG(1,2)
                        YG(5,ng,ne)=TG(1,3)
                        YG(6,ng,ne)=TG(2,3)
                      ELSE IF(TYPE2(1:7).EQ.'Nominal') THEN
                        YG(1,ng,ne)=TN(1,1)
                        YG(2,ng,ne)=TN(2,2)
                        YG(3,ng,ne)=TN(3,3)
                        YG(4,ng,ne)=TN(1,2)
                        YG(5,ng,ne)=TN(1,3)
                        YG(6,ng,ne)=TN(2,3)
                      ENDIF
                    ENDIF
                  ELSE IF(TYPE(1:13).EQ.'STRAIN_ENERGY') THEN
                    CALL ZEEX50(COORDS,IBT,IDO,INP,NAN,
     '                NBH(1,1,ne),NBJ(1,ne),ng,NHE(ne,nx),
     '                NPNE(1,1,ne),nr,nx,
     '                DZDX,CE(1,ne,nx),CG,CP(1,1,nx),EG,PG,PHI,PST,
     '                R,RGX(ng),RI1,RI2,RI3,RM,U,
     '                XE,XG,XIG(1,ng,nb_stress),
     '                ZE,ZG,ERROR,*9999)

                    ! calculate the deviatoric component of Green's strain
                    ! C = 2*E+I, Cd = det(C)**(-1/3)*C, Ed = (1/2)*(Cd-I)
                    ! volumetric component is EG-Ed
                    DO i=1,3
                      DO j=1,3
                        C(i,j)=2*EG(i,j)
                        IF(i.EQ.j) C(i,j)=C(i,j)+1
                      ENDDO
                    ENDDO
                    detC=DET(C)**(-1.0D0/3.0D0)

                    DO i=1,3
                      DO j=1,3
                        Ed(i,j)=0.5d0*(detC*C(i,j))
                        IF(i.EQ.j) Ed(i,j)=Ed(i,j)-0.5d0
                        Ev(i,j)=EG(i,j)-Ed(i,j)
                      ENDDO
                    ENDDO

                    CALL STRAIN_ENERGY(nr,ne,ng,Ed,se_dev,ERROR,*9999)
                    CALL STRAIN_ENERGY(nr,ne,ng,Ev,se_vol,ERROR,*9999)
                    CALL STRAIN_ENERGY(nr,ne,ng,EG,se_gauss,ERROR,
     &                *9999)

                    YG(1,ng,ne)=se_gauss ! strain energy density
                    YG(2,ng,ne)=se_dev
                    YG(3,ng,ne)=se_vol
                    YG(4,ng,ne)=0
                    YG(5,ng,ne)=0
                    YG(6,ng,ne)=0

                  ELSE IF(TYPE(1:6).EQ.'STRAIN') THEN
C new MPN 26Feb98: (see above)
                    CALL ZEEX50(COORDS,IBT,IDO,INP,NAN,
     '                NBH(1,1,ne),NBJ(1,ne),ng,NHE(ne,nx),
     '                NPNE(1,1,ne),nr,nx,
     '                DZDX,CE(1,ne,nx),CG,CP(1,1,nx),EG,PG,PHI,PST,
     '                R,RGX(ng),RI1,RI2,RI3,RM,U,
     '                XE,XG,XIG(1,ng,nb_stress),
     '                ZE,ZG,ERROR,*9999)
C old
C                    CALL ZEEX50(COORDS,IBT,IDO,INP,NAN,
C     '                NBH_temp(1,1),NBJ_temp,ng,NHE(ne,nx),
C     '                NPNE(1,1,ne),nr,nx,
C     '                DZDX,CE(1,ne,nx),CG,CP(1,1,nx),EG,PG,PHI,PST,
C     '                R,RGX(ng),RI1,RI2,RI3,RM,U,
C     '                XE,XG,XIG(1,ng,nb_stress),
C     '                ZE,ZG,ERROR,*9999)
                    IF(PRINCIPAL) THEN
                      YG(1,ng,ne)=PST(1) !max principal strain
                      YG(2,ng,ne)=PST(2) !2nd principal strain
                      YG(3,ng,ne)=PST(3) !min principal strain
                      YG(4,ng,ne)=PHI(1) !1st principal angle
                      YG(5,ng,ne)=PHI(2) !2nd principal angle
                      YG(6,ng,ne)=PHI(3) !3rd principal angle
                    ELSE
                      IF(TYPE2(1:10).EQ.'COMPONENTS') THEN
                        YG(1,ng,ne)=EG(1,1)
                        YG(2,ng,ne)=EG(2,2)
                        YG(3,ng,ne)=EG(3,3)
                        YG(4,ng,ne)=EG(1,2)
                        YG(5,ng,ne)=EG(1,3)
                        YG(6,ng,ne)=EG(2,3)
                      ELSE IF(TYPE2(1:16).EQ.'EXTENSION_RATIOS') THEN
C CS 31/10/00 adding DABS because sometime shear terms become negative
                        YG(1,ng,ne)=DSQRT(DABS(2.0d0*EG(1,1)+1.0d0))
                        YG(2,ng,ne)=DSQRT(DABS(2.0d0*EG(2,2)+1.0d0))
                        YG(3,ng,ne)=DSQRT(DABS(2.0d0*EG(3,3)+1.0d0))
                        YG(4,ng,ne)=DSQRT(DABS(2.0d0*EG(1,2)+1.0d0))
                        YG(5,ng,ne)=DSQRT(DABS(2.0d0*EG(1,3)+1.0d0))
                        YG(6,ng,ne)=DSQRT(DABS(2.0d0*EG(2,3)+1.0d0))
                      ELSE IF(TYPE2(1:10).EQ.'INVARIANTS') THEN
                        YG(1,ng,ne)=RI3
                        YG(2,ng,ne)=RI2
                        YG(3,ng,ne)=RI1
                      ENDIF
                    ENDIF
                  ENDIF

                ENDDO
!newe
              ENDIF
            ENDDO

          ELSE IF(TYPE(1:8).EQ.'PRESSURE') THEN
            nc=1 !Temporary AJP 18-12-91
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              CALL ZPZE(NBH(1,1,ne),nc,NHE(ne,nx),NKHE(1,1,1,ne),
     '          NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1,nx),nx,
     '          CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,
     '          ERROR,*9999)
C MPN 26Feb98: avoid confusion with nb in other places
C old              nb=NBH(nh,nc,ne)
C old              DO ng=1,NGT(nb)
              DO ng=1,NGT(NBH(nh,nc,ne))
                CALL XEXG(NBJ(1,ne),ng,nr,PG,
     '            XE,XG,ERROR,*9999)
                CALL XGMG(1,0,NBJ(1,ne),nr,DXIX,GL,GU,
     '            RGX(1),XG,ERROR,*9999)
                CALL ZEZG(1,NBH(1,nc,ne),ng,NHE(ne,nx),nx,DXIX,PG,
     '            ZE,ZG,ERROR,*9999)
                YG(nh,ng,ne)=0.5d0*(ZG(nh,2)**2.d0+ZG(nh,4)**2.d0)
              ENDDO
            ENDDO

          ELSE IF(TYPE(1:10).EQ.'ACTIVATION') THEN
            IF(ITYP5(nr,nx).EQ.2) THEN
              IF(.NOT.CALL_GMAP) THEN
                CALL GEN_GRID_MAP(IBT,IDO,INP,NQET,NQSCNB,NQXI,PGNQE,
     '            XIG,ERROR,*9999)
                CALL_GMAP=.TRUE.
              ENDIF

              CALL MAQ_LOC(MAQ_INQUIRE,MAQ_TIME,atime,MAQ_ACTIV_TIME,
     '          ERROR,*9999)

              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)

C old MLB 18 August 1997
c               WRITE(OP_STRING,'('' Note: Gauss array YG is '
c     '           //'updated from THRES(2&3,ng,ne) only as YG(1,ng,ne) '
c     '           //'already defined'')')
c               CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C                DO ng=1,NGT(nb_extended)
C                  IF(ITYP2(nr,nx).EQ.8) THEN      !Threshold
C                    DO I=2,3
C                      YG(I,ng,ne)=THRES(I,ng,ne)
C                    ENDDO
C                  ELSE IF(ITYP2(nr,nx).EQ.9) THEN !Bidomain
C                    nq=NQGE(ng,ne,nb_extended)
C                    IF(NWQ(1,nq,1).EQ.0) THEN !Internal
C                      YG(1,ng,ne)=YQ(nq,1,3,nx) !Activation time
C                    ELSE !Boundary
C                      YG(1,ng,ne)=YQ(NWQ(1,nq,1),1,3,nx) !copy neighbour
C                    ENDIF
C                  ENDIF
C                ENDDO

                IF(ITYP2(nr,nx).EQ.8) THEN !Threshold
                  WRITE(OP_STRING,'('' Threshold updates not '
     '              //'currently supported'')')
                ENDIF

                SCHEME=NQS(ne)
C MPN 26Feb98: avoid confusion with nb in other places
C old                nb=NQSCNB(SCHEME)
C old                DO ng=1,NGT(nb)
                DO ng=1,NGT(NQSCNB(SCHEME))
                  YG(1,ng,ne)=0.0d0
                  DO nqq=1,NQET(SCHEME)
                    YG(1,ng,ne)=YG(1,ng,ne)+(AQ(atime,NQNE(ne,nqq))*
     '                PGNQE(ng,nqq,SCHEME))
                  ENDDO !niq
                ENDDO !ng
              ENDDO !noelem
            ENDIF

          ELSE IF(TYPE(1:9).EQ.'POTENTIAL') THEN
            IF(ITYP5(nr,nx).EQ.2) THEN !time-dependent
              IF(ITYP19(nr,nx).EQ.1.AND.ITYP2(nr,nx).EQ.9) THEN !Bidomain
C MLB 18 August 1997
C obselete
C                DO noelem=1,NEELEM(0,nr)
C                  ne=NEELEM(noelem,nr)
C                  DO ng=1,NGT(nb_extended)
C                    nq=NQGE(ng,ne,nb_extended)
C                    DO niq=1,MIN(NIQM,NJM)
C                      YG(niq,ng,ne)=YQ(nq,niq,1,nx)
C                    ENDDO
C                  ENDDO !ng
C                ENDDO !noelem

                IF(.NOT.CALL_GMAP) THEN
                  CALL GEN_GRID_MAP(IBT,IDO,INP,NQET,NQSCNB,NQXI,PGNQE,
     '              XIG,ERROR,*9999)
                  CALL_GMAP=.TRUE.
                ENDIF

                CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqV,NIQ_V,ERROR,*9999)

                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  SCHEME=NQS(ne)
C MPN 26Feb98: avoid confusion with nb in other places
C old                  nb=NQSCNB(SCHEME)
C old                  DO ng=1,NGT(nb)
                  DO ng=1,NGT(NQSCNB(SCHEME))
                    YG(1,ng,ne)=0.0d0
                    DO nqq=1,NQET(SCHEME)
                      YG(1,ng,ne)=YG(1,ng,ne)+
     '                  (YQ(NQNE(ne,nqq),niqV,1,nx)*
     '                  PGNQE(ng,nqq,SCHEME))
                    ENDDO !niq
                  ENDDO !ng
                ENDDO !noelem
              ENDIF !ityp2
            ENDIF !ityp5

          ELSE IF(TYPE(1:8).EQ.'MATERIAL') THEN
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              IF(ITYP1(nr,nx).EQ.3) THEN
                CALL CPCG(1,NBH(1,1,ne),NPNE(1,1,ne),nr,nx,
     '            CE(1,ne,nx),CG,CGE(1,1,ne,nx),
     '            CP(1,1,nx),PG,ERROR,*9999)
              ELSE
C CS 23/9/00 Don't quite understand this. NW seems to be wrong
C                CALL CPCG(NW(ne,1,nx),NBH(1,1,ne),NPNE(1,1,ne),nr,nx,
C     '            CE(1,ne,nx),CG,CGE(1,1,ne,nx),
C     '            CP(1,1,nx),PG,ERROR,*9999)
                CALL CPCG(1,NBH(1,1,ne),NPNE(1,1,ne),nr,nx,
     '            CE(1,ne,nx),CG,CGE(1,1,ne,nx),
     '            CP(1,1,nx),PG,ERROR,*9999)
              ENDIF
C MPN 26Feb98: avoid confusion with nb in other places
C old              nb=NBJ(1,ne)
C old              DO ng=1,NGT(nb)
              IF((ILPIN(NM1,nr,nx).EQ.4).OR.
     '          (ILPIN(NM1,nr,nx).EQ.5)) THEN !defined by Gauss/grid
                DO ng=1,NGT(NBJ(1,ne))
                  YG(NM2,ng,ne)=CGE(NM1,ng,ne,nx)
                ENDDO
              ELSE
                DO ng=1,NGT(NBJ(1,ne))
                  YG(NM2,ng,ne)=CG(NM1,ng)
                ENDDO
              ENDIF
            ENDDO

          ELSE IF(TYPE(1:6).EQ.'SOURCE') THEN
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '          NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '          SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
              DO ng=1,NGT(NBJ(1,ne))
                CALL XEXG(NBJ(1,ne),ng,nr,PG,
     '            XE,XG,ERROR,*9999)
                GX=XG(1,1)
                GY=XG(2,1)
                IF(IFNTYP.EQ.1) THEN
                  YG(1,ng,ne)=-2.d0*PI*PI*DSIN(PI*GX)*DSIN(PI*GY)
                ENDIF !ifntyp
              ENDDO !ng
            ENDDO !noelem

          ELSE IF(TYPE(1:7).EQ.'CALCIUM') THEN
C**   Compute calcium levels for mechanics at gauss points.
C**   Multiply by calcium fudge factor (temporary, hopefully!)
            IF(ITYP3(nr,nx_s).GE.2) THEN ! FHN/VCD/BR/LR
              IF(.NOT.CALL_GMAP) THEN
                CALL GEN_GRID_MAP(IBT,IDO,INP,NQET,NQSCNB,NQXI,PGNQE,
     '            XIG,ERROR,*9999)
                CALL_GMAP=.TRUE.
              ENDIF
              IF(ITYP3(nr,nx_s).EQ.2) THEN !fhn
C                CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqCa,NIQ_CALCIUM,
C     '            ERROR,*9999)
                YQS_INDEX=CELL_STATE_OFFSET(1)+3-1
              ELSEIF(ITYP3(nr,nx_s).EQ.3) THEN !vcd
C                CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqCa,NIQ_CALCIUM,
C     '            ERROR,*9999)
                YQS_INDEX=CELL_STATE_OFFSET(1)+3-1
              ELSEIF(ITYP3(nr,nx_s).EQ.4) THEN !br
C                CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqCa,NIQ_CAI,
C     '            ERROR,*9999)
                YQS_INDEX=CELL_STATE_OFFSET(1)+8-1
              ELSEIF(ITYP3(nr,nx_s).EQ.5) THEN !jrw
C                CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqCa,NIQ_CAI,
C     '            ERROR,*9999)
                CALL ASSERT(.FALSE.,' >>Not implemented',ERROR,*9999)
              ELSEIF(ITYP3(nr,nx_s).EQ.6) THEN !lr
C                CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqCa,NIQ_CAI,
C     '            ERROR,*9999)
                CALL ASSERT(.FALSE.,' >>Not implemented',ERROR,*9999)
C!!! KAT: not merging from Oxford
C              ELSEIF((ITYP3(nr,nx_s).EQ.10).AND.(KTYP33.EQ.5)) THEN !user_cell5
CC                CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqCa,NIQ_CAI,
CC     '            ERROR,*9999)
C                YQS_INDEX=Cai !need to check this with PJH, NPS May  9 2000
              ELSEIF(ITYP3(nr,nx_s).EQ.8) THEN !N98
C *** DPN 25 September 2000 - this is all bad, but need to get calcium
C ***   values into FEXT somehow when using N98 on its own - can be
C ***   fixed by using N98-HMT and sticking tension directly into YG
C *** Also, all the above are wrong - should only be using
C *** CELL_STATE_OFFSET !!!!

C *** DPN 17 January 2003 - fixed the above by changing the definitions
C               of CELL_STATE_OFFSET(variant) in the *_INIT_GRID
C               routines for the simplified models. And now assuming
C               variant 1 here for N98.

                YQS_INDEX=CELL_STATE_OFFSET(1)+4-1
              ENDIF
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                SCHEME=NQS(ne)
C MPN 26Feb98: avoid confusion with nb in other places
C old                nb=NQSCNB(SCHEME)
C old                DO ng=1,NGT(nb)
                DO ng=1,NGT(NQSCNB(SCHEME))
                  IF(ne.EQ.Infarct) THEN !Elem infarcted so Ca=0
                    FEXT(4,ng,ne)=0.0d0
                  ELSE
                    FEXT(4,ng,ne)=0.0d0
                    DO nqq=1,NQET(SCHEME)
C                      FEXT(4,ng,ne)=FEXT(4,ng,ne)+
C     '                  (YQ(NQNE(ne,nqq),niqCa,1,nx_s)
C     '                    *PGNQE(ng,nqq,SCHEME))*Ca_factor
                      FEXT(4,ng,ne)=FEXT(4,ng,ne)+
     '                  (YQS(YQS_INDEX,NQNE(ne,nqq))*
     '                  PGNQE(ng,nqq,SCHEME))*Ca_factor
                    ENDDO !nqq
                  ENDIF !infarct
                ENDDO !ng
              ENDDO !noelem

C old MLB 18 August 1997
C              DO noelem=1,NEELEM(0,nr)
C                ne=NEELEM(noelem,nr)
C                nb=NBJ(1,ne)
C loop over gauss pts for mechanics basis
C                DO ng=1,NGT(nb)
C                  nq=NQGE(ng,ne,nb) !Closest grid point to (ng,ne)
C                  FEXT(4,ng,ne)=YQ(nq,4,1,nx_s)*Ca_factor
C                ENDDO !ng
C              ENDDO !ne

            ELSE
              WRITE(OP_STRING,
     '          '('' >>>Requires a model describing calcium'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF !ityp3(nr,nx).ge.2

C news MPN 1Jun2000: copying/calc'ing values of gridvars at Gauss pts
C                    This may eventually be used to supercede the
C                    CALCIUM option immediately above
C                    (still under debate!)
          ELSE IF(TYPE(1:8).EQ.'GRIDVARS') THEN

             ! If using gridpoints at gauss we don't need mapping
             IF(KTYP3B.EQ.1) THEN
                IF(.NOT.CALL_GMAP) THEN
                   CALL GEN_GRID_MAP(IBT,IDO,INP,NQET,NQSCNB,NQXI,PGNQE,
     '                  XIG,ERROR,*9999)
                   CALL_GMAP=.TRUE.
                ENDIF
             ENDIF

            DO niqq=1,NQLIST(0)
              niqs=NQLIST(niqq)
              nig=NGLIST(niqq)
C VJ news 28Jan2004: Modified do loop to loop over NELIST_INCLUDE alone
C              DO noelem=1,NEELEM(0,nr)
C                ne=NEELEM(noelem,nr)
C CS 25/6/2002 This rearrangement seems to fix AIX optimiser problems
C                NE_EXCLUDE=.FALSE.
C                DO noelem1=1,NELIST(0)
C                  ne1=NELIST(noelem1)
C                  IF(ne1.EQ.ne) NE_EXCLUDE=.TRUE. !exclude this element
C                ENDDO !noelem1 (ne1)
C                IF(INLIST(ne,NELIST(1),NELIST(0),N1LIST)) THEN
C                  NE_EXCLUDE=.TRUE.
C                ELSE
C                  NE_EXCLUDE=.FALSE.
C                ENDIF
C KAT 15Feb00:  Initializing YG here even if EXCLUDEd because it isn't
C               set anywhere else (example b362).
C                IF(.NOT.NE_EXCLUDE) THEN
C                SCHEME=NQS(ne)
C                DO ng=1,NGT(NQSCNB(SCHEME))
C                  YG(nig,ng,ne)=0.0d0
C                  IF(.NOT.NE_EXCLUDE) THEN
C                    IF(KTYP3B.EQ.1) THEN
C                      DO nqq=1,NQET(SCHEME)
C *** DPN 26 March 2001 - Adding ability to scale the contribution from
C *** YQS to YG
c                        YG(nig,ng,ne)=YG(nig,ng,ne)+
c     '                    (YQS(niqs,NQNE(ne,nqq))*PGNQE(ng,nqq,SCHEME))
C                        YG(nig,ng,ne)=YG(nig,ng,ne)+
C     '                    (YQS(niqs,NQNE(ne,nqq))*YQS_FACTOR*
C     '                    PGNQE(ng,nqq,SCHEME))
C                      ENDDO !nqq
C                    ELSE ! Gauss point grid scheme
C                      YG(nig,ng,ne)=YQS(niqs,NQNE(ne,ng))*YQS_FACTOR
C                    ENDIF
C                  ENDIF !not ne_exclude
C                ENDDO !ng
C                ENDIF !not ne_exclude
C              ENDDO !noelem (ne)
              DO noelem=1,NELIST_INCLUDE(0)
                ne=NELIST_INCLUDE(noelem)
                SCHEME=NQS(ne)
                DO ng=1,NGT(NQSCNB(SCHEME))
                  YG(nig,ng,ne)=0.0d0
                  IF(KTYP3B.EQ.1) THEN
                    DO nqq=1,NQET(SCHEME)
C *** DPN 26 March 2001 - Adding ability to scale the contribution from
C *** YQS to YG
c                      YG(nig,ng,ne)=YG(nig,ng,ne)+
c     '                  (YQS(niqs,NQNE(ne,nqq))*PGNQE(ng,nqq,SCHEME))
                      YG(nig,ng,ne)=YG(nig,ng,ne)+
     '                  (YQS(niqs,NQNE(ne,nqq))*YQS_FACTOR*
     '                  PGNQE(ng,nqq,SCHEME))
                    ENDDO !nqq
                  ELSE ! Gauss point grid scheme
                    YG(nig,ng,ne)=YQS(niqs,NQNE(ne,ng))*YQS_FACTOR
                  ENDIF
                ENDDO !ng
              ENDDO !noelem
C set YG values to 0 for excluded elements as well. Important for malloc and malloc7 executables
              DO noelem=1,NELIST_EXCLUDE(0)
                ne=NELIST_EXCLUDE(noelem)
                SCHEME=NQS(ne)
                DO ng=1,NGT(NQSCNB(SCHEME))
                  YG(nig,ng,ne)=0.0d0
                ENDDO
              ENDDO
C newe VJ
            ENDDO !niqq (niqs)
C newe MPN 1Jun2000

C *** DPN 28 April 2001 - Manipulation of YG fields
          ELSE IF(TYPE(1:11).EQ.'GAUSS_FIELD') THEN
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              SCHEME=NQS(ne)
              DO ng=1,NGT(NQSCNB(SCHEME))
                IF(MANIPULATION(1:10).EQ.'SUBSTITUTE') THEN
                  YG(DESTINATION_FIELD,ng,ne) =
     '              YG(SOURCE_FIELD,ng,ne)*YG_FACTOR
                ELSE IF(MANIPULATION(1:22).EQ.'DIFFERENCE')
     '              THEN
                  YG(DESTINATION_FIELD,ng,ne) =
     '              YG(DESTINATION_FIELD,ng,ne) +
     '              (YG(SOURCE_FIELD,ng,ne)-YG(DESTINATION_FIELD,ng,ne))
     '              * YG_FACTOR
                ENDIF
              ENDDO ! ng=1,NGT(NQSCNB(SCHEME))
            ENDDO ! noelem=1,NEELEM(0,nr)

          ELSE IF(TYPE(1:9).EQ.'DEFFIBRES') THEN
            nc=1 !LHS
            ERROR_FLAG=.FALSE.
C cpb 18/10/96 Dynamically allocate parallel local arrays
C$OMP       PARALLEL DO
C$OMP&        PRIVATE(nb,ne,ng,ng_calc,ni,NITB,noelem,
C$OMP&          ADdotAS,AD_VECTOR,ALFA,ASdotCS,ASdotFD,AS_VECTOR,
C$OMP&          BD_VECTOR,BETA,BSdotFD,BS_VECTOR,BXdotBS,BX_VECTOR,
C$OMP&          CD_VECTOR,CS_VECTOR,CXdotBS,CX_VECTOR,FD_VECTOR,
C$OMP&          GAMA,GD_VECTOR,HD_VECTOR,IODI,MEM_INIT,XE_PTR,XG_PTR,XI,
C$OMP&          ZE_PTR,ZG_PTR,ERROR)
C$OMP&        SHARED(IBT,IDO,INP,NAN,nc,NBH,NBJ,NEELEM,NGT,NHE,NIT,
C$OMP&          NKHE,NKJE,NPF,NPNE,nr,NVHE,NVJE,NW,nx,
C$OMP&          PG,PI,SE,XA,XIG,XP,YG,ZA,ZP,COLLOCATION,ERROR_FLAG)
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              IF(.NOT.ERROR_FLAG) THEN

C CPB 18/10/96 Intialise pointers so they are zero in the parallel loop

                XE_PTR=0
                XG_PTR=0
                ZE_PTR=0
                ZG_PTR=0

c cpb 18/10/96  Dynamically allocation parallel local arrays

                CALL ALLOCATE_MEMORY(NJM*NSM,1,DPTYPE,XE_PTR,
     '            MEM_INIT,ERROR,*100)
                CALL ALLOCATE_MEMORY(NJM*NUM,1,DPTYPE,XG_PTR,
     '            MEM_INIT,ERROR,*100)
                CALL ALLOCATE_MEMORY(NSM*NHM,1,DPTYPE,ZE_PTR,
     '            MEM_INIT,ERROR,*100)
                CALL ALLOCATE_MEMORY(NHM*NUM,1,DPTYPE,ZG_PTR,
     '            MEM_INIT,ERROR,*100)

C cpb 18/10/96 Use dynamically allocated parallel local arrays
C                CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
C     '            NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
C     '            SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*100)
C                CALL ZPZE(NBH(1,1,ne),nc,NHE(ne,nx),NKHE(1,1,1,ne),
C     '            NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),
C     '            NW(ne,1,nx),nx,
C     '            CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,
C     '            ERROR,*100)
                CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '            NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '            SE(1,1,ne),XA(1,1,ne),%VAL(XE_PTR),XP,ERROR,*100)
                CALL ZPZE(NBH(1,1,ne),nc,NHE(ne,nx),NKHE(1,1,1,ne),
     '            NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),
     '            NW(ne,1,nx),nx,
     '            CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),
     '            %VAL(ZE_PTR),ZP,ERROR,*100)
                IF(COLLOCATION) THEN
                  nb=NQSCNB(NQS(ne))
                ELSE
                  nb=NBH(NH_LOC(1,nx),nc,ne)
                ENDIF
                NITB=NIT(nb)
                DO ng=1,NGT(nb)
C                  IF(.NOT.COLLOCATION.OR.
C     '              NQGE(ng,ne,nb).GT.0) THEN !nq is defined for ng,ne
                    IF(COLLOCATION) THEN
                      DO ni=1,NITB
                        XI(ni)=XIG(ni,ng,nb)
                      ENDDO
                      ng_calc=0
                    ELSE
                      ng_calc=ng
C cpb 18/10/96 Use dynamically allocated parallel local arrays
C                      CALL XEXG(NBJ(1,ne),ng,nr,PG,XE,XG,
C     '                  ERROR,*100)
                      CALL XEXG(NBJ(1,ne),ng,nr,PG,
     '                  %VAL(XE_PTR),%VAL(XG_PTR),ERROR,*100)
                    ENDIF
C                   Form orthonormal set of deformed material vectors
C                   wrt rc coords at ng or at XI(ni) for collocation.
                    CALL MAT_VEC_DEF(IBT,IDO,INP,NAN,NBH(1,nc,ne),
     '                NBJ(1,ne),ng_calc,NHE(ne,nx),nr,nx,AD_VECTOR,
     '                BD_VECTOR,CD_VECTOR,PG,%VAL(XE_PTR),%VAL(XG_PTR),
     '                XI,%VAL(ZE_PTR),%VAL(ZG_PTR),.TRUE.,ERROR,*100)
C                   Compute deformed wall ref vecs wrt rc coords at ng
                    CALL FIBRE_REF_VECS_DEF(IBT,IDO,INP,NAN,
     '                NBH(1,nc,ne),ng_calc,NHE(ne,nx),NITB,nr,nx,
     '                FD_VECTOR,GD_VECTOR,HD_VECTOR,
     '                PG,XI,%VAL(ZE_PTR),%VAL(ZG_PTR),ERROR,*100)
C                   Calc BS_VECTOR = BD_VECTOR rotated about AD_VECTOR
C                   into (Xi1,Xi2) plane with the unit normal, HD_VECTOR
                    CALL CROSS(HD_VECTOR,AD_VECTOR,BS_VECTOR)
                    CALL NORMALISE(3,BS_VECTOR,ERROR,*100)
                    CALL CROSS(AD_VECTOR,BS_VECTOR,CS_VECTOR) !req below
C                   Calc deformed sheet angle
                    BXdotBS=DDOT(3,BX_VECTOR,1,BS_VECTOR,1)
                    GAMA=DACOS_MOD(BXdotBS)
                    CXdotBS=DDOT(3,CX_VECTOR,1,BS_VECTOR,1)
                    IF(DACOS_MOD(CXdotBS).LT.(PI/2.0d0)) GAMA=-GAMA
C                   Calc AS_VECTOR = AD_VECTOR rotated BS_VECTOR
C                   into (Xi1,Xi2) plane with the unit normal, HD_VECTOR
                    CALL CROSS(BS_VECTOR,HD_VECTOR,AS_VECTOR)
                    CALL NORMALISE(3,AS_VECTOR,ERROR,*100)
C                   Calc deformed imbrication angle
                    ADdotAS=DDOT(3,AD_VECTOR,1,AS_VECTOR,1)
                    BETA=DACOS_MOD(ADdotAS)
                    ASdotCS=DDOT(3,AS_VECTOR,1,CS_VECTOR,1)
                    IF(DACOS_MOD(ASdotCS).GT.(PI/2.0d0)) BETA=-BETA
C                   Calc deformed fibre angle
                    ASdotFD=DDOT(3,AS_VECTOR,1,FD_VECTOR,1)
                    ALFA=DACOS_MOD(ASdotFD)
                    BSdotFD=DDOT(3,BS_VECTOR,1,FD_VECTOR,1)
                    IF(DACOS_MOD(BSdotFD).LT.(PI/2.0d0)) ALFA=-ALFA
C                   Store deformed material axis angles in YG
                    YG(1,ng,ne)=ALFA !fibre angle
                    YG(2,ng,ne)=BETA !imbrication angle
C MPN 2Nov2001: +PI/2 no longer necessary as sheet angles fixed.
                    YG(3,ng,ne)=GAMA !sheet angle
C                    YG(3,ng,ne)=GAMA+PI/2.0d0 !sheet angle
CC                   NOTE the +PI/2 is temporary
CC                   (see comment in MAT_VEC_ROTATE, FE02)
C                  ENDIF !.NOT.COLLOCATION or nq exists for this ng,ne
                ENDDO !ng

C cpb 18/10/96 Free dynamically allocated arrays

                CALL FREE_MEMORY(XE_PTR,ERROR,*100)
                CALL FREE_MEMORY(XG_PTR,ERROR,*100)
                CALL FREE_MEMORY(ZE_PTR,ERROR,*100)
                CALL FREE_MEMORY(ZG_PTR,ERROR,*100)

                GO TO 102
C                 This statement is designed to be skipped if no error
C                 occur. However if a error occurs within a subroutine
C                 the alternate return jumps to line 100 to set the flag
 100              CONTINUE
C KAT 14May01: mp_setlock not OPENMP.
C              Critical section is not essential.
CC$                call mp_setlock()
                  ERROR_FLAG=.TRUE.
                  WRITE(OP_STRING,'(/'' >>ERROR: An error occurred - '
     '              //'results may be unreliable!'')')
                  CALL WRITES(IODI,OP_STRING,ERROR,*101)
 101              CONTINUE
CC$                call mp_unsetlock()
 102            CONTINUE
              ENDIF !.NOT.ERROR_FLAG
            ENDDO !noelem
            CALL ASSERT(.NOT.ERROR_FLAG,'>>An error occurred during '
     '        //'update gauss deformed_fibres',ERROR,*9999)

          ELSE IF(TYPE(1:7).EQ.'EIKONAL') THEN
            IF(OP_TIMING.AND.NRLIST(0).GT.1) THEN
C             Region heading for timing information.
              WRITE(OP_STRING,'(/'' Region '',I1,'':'')') nr
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            IF(ITYP5(nr,nx).EQ.5.AND.ITYP2(nr,nx).EQ.3) THEN !eikonal equation
              CALL UPGAUS_EIK(IBT,IDO,NBH,NBHF,NBJ,NBJF,NEELEM(0,nr),
     '          NFFACE(0,nr),NHE,NKB,NKJE,NNF,NPF,NPNE,nr,NSB,NVJE,
     '          NW(1,1,nx),nx,
     '          CE(1,1,nx),CG,CGE(1,1,1,nx),CP(1,1,nx),PG,SE,WG,
     '          XA,XP,YG,YGF,OP_TIMING,ERROR,*9999)
            ELSE
              IF(OP_TIMING.OR.NRLIST(0).EQ.1) THEN
                WRITE(OP_STRING,'(/'' The equation is not eikonal'')')
              ELSE
                WRITE(OP_STRING,'('' The equation for region '',I1,'//
     '            ''' is not eikonal'')') nr
              ENDIF
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF !TYPE
        ENDDO !no_nrlist
        IF(TYPE(1:7).EQ.'EIKONAL') THEN
          CALL_UPGAUS_EIK=.TRUE.
        ENDIF
      ENDIF

      CALL EXITS('UPGAUS')
      RETURN
 9999 IF(XE_PTR.NE.0) CALL FREE_MEMORY(XE_PTR,ERROR,*1111)
      IF(XG_PTR.NE.0) CALL FREE_MEMORY(XG_PTR,ERROR,*1111)
      IF(ZE_PTR.NE.0) CALL FREE_MEMORY(ZE_PTR,ERROR,*1111)
      IF(ZG_PTR.NE.0) CALL FREE_MEMORY(ZG_PTR,ERROR,*1111)
 1111 CALL ERRORS('UPGAUS',ERROR)
      RET_ERROR=ERROR
      CALL EXITS('UPGAUS')
      RETURN 1
      END


