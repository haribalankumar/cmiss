      SUBROUTINE EVFIBR(IBT,IDO,INP,LD,NAN,NBJ,NEELEM,NENP,NKJE,NKJ,
     '  NNB,NPF,NPNE,NPNODE,NVJE,NXI,SE,STRING,WD,XA,XE,XG,XG1,
     '  XID,XP,ZD,ERROR,*)

C#### Subroutine: EVFIBR
C###  Description:
C###    <HTML>
C###    EVFIBR evaluates fibre angle at a specified point in a prolate
C###    spheriodal mesh for which lambda varies linearly with Xi3 and
C###    mu and theta vary tri-linearly.
C###    <BR>
C###    <PRE>
C###     XIPOS are the material Xi positions of the point.
C###     XPOS(2),XPOS(3) are the mu and theta coords of the point.
C###     S3 is the proportion of the way thru the wall to interpolate to.
C###     NUMELEM is a list of all the elements for which the point lies in
C###       the Xi3=0 and Xi3=1 faces (different lambda's).
C###     XG1LIST store the lamda values of the elements in the above list
C###       at the (mu,theta) position on the Xi3=0 and Xi3=1 faces.
C###     XIF1LIST/XIF2LIST store the Xi1 and Xi2 coords of the point on
C###       Xi3=0/Xi3=1 face for the elements in the list (XIFACE1/XIFACE2
C###       store these coords temporarily).
C###     INDEXELEM stores indices to elements in the list in order from
C###       endo to epi ie)INDEXELEM(1)      is most endo element,
C###                      INDEXELEM(NINDEX) is most epi  element.
C###     BASALLIST stores true for elements around the basal ring when mu
C###       is > mu on the basal face. ie) if mu>=120 then BASALLIST will be
C###       true for the elements in the list if the max mu coord for the
C###       mesh is 115 (BASAL stores this temporarily)
C###     WALL_THICK is a measure of the wall thickness at the mu,theta pt.
C###     WALL_POS is the lambda position to interpolate to.
C###     XRECT are the rectangular cart. coords of the pt, from which mu
C###       and theta coords are to be calculated (lambda is ignored)
C###     DATA files should be have the format:
C###              first  column: data point number
C###              second column: X or LAMBDA
C###              third  column: Y or MU
C###              fourth column: Z or THETA
C###              fifth  column: S3 (proportion of dist through the wall)
C###              (sixth and seventh columns unused)
C###       where (X,Y,Z) are used if the RECTANGULAR_CART option is given
C###       and (LAMBDA,MU,THETA) are used for the PROLATE_SPHEROIDAL option
C###    </PRE>

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cmis00.cmn'
      INCLUDE 'data00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'head00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LD(NDM),NAN(NIM,NAM,NBFM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),NKJE(NKM,NNM,NJM,NEM),NKJ(NJM,NPM),
     '  NNB(4,4,4,NBFM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NVJE(NNM,NBFM,NJM,NEM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 SE(NSM,NBFM,NEM),WD(NJM,NDM),XA(NAM,NJM),XE(NSM,NJM),
     '  XG(NJM,NUM),XG1(NJM,NUM,NGM),XID(NIM,NDM),XP(NKM,NVM,NJM,NPM),
     '  ZD(NJM,NDM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER MAXWALLELEM,NJMAX,NKMAX
      PARAMETER(MAXWALLELEM=10,NJMAX=6,NKMAX=8)
      INTEGER i,IBEG,IBEG1,IEND,IEND1,
     '  iface,IFROMC,INDEX,INDEXELEM(MAXWALLELEM),J,
     '  N3CO,nb,nd,NDMAX,ne,NEENDO,NE_EVFIBR,ni,NINDEX,
     '  nj,njj,NJ1,NJ2,NJ3,nk,NKTB,
     '  noelem,nr,NS1,NS2,NS3,NS4,
     '  NTLIST,nu,NU1(3),NUMELEM(0:MAXWALLELEM),WRT_XI,
     '  IOFI2  !may need to be put into cbdi02.cmn at some stage
      REAL*8 A1,A2,ALFA,ALPHA,B1,B2,BETA,C1,C2,COSHL,COSM,COST,D1,D2,
     '  DELTA,DENOM1,DENOM2,DSDXI(3),DXDXI(3),DXDXI1(3),DXDXI2(3),
     '  FIBRE_VECT(3),FOCUSDATA,FOCUSTEMP,GAMA,MU,
     '  PXI,RFROMC,RR,S3,SCAFAC(NJMAX,NKMAX),
     '  Sheet_angle_max,Sheet_angle_min,SINHL,SINM,SINT,
     '  THETA,THETAMAX,THETAMIN,TOL,WALL_POS,WALL_THICK,
     '  X1,X2,XG1MIN,XG1LIST(MAXWALLELEM,2),Xi(3),Radius_xi,
     '  XIF1LIST(MAXWALLELEM,2),XIF2LIST(MAXWALLELEM,2),
     '  XIFACE1(2),XIFACE2(2),XIPOS(3),XPOS(3),XRECT(3),XX,
     '  NORM1,NORM2,NORM3,FDOTXI1,FDOTXI2
      CHARACTER CHAR1*1,CHAR2*4,CHAR3*11,CHAR4*100,CHAR6*5,COORDS*9,
     '  CTEMP*30,FIELD(3)*5,FILE*(MXCH),FORMAT*132,INPUT_TYPE*11,
     '  MODE*7,TYPE*5
      LOGICAL BASAL,BASALLIST(MAXWALLELEM),CBBREV,CORRECT,CURVIL,
     '  DATAPTS,EXCLUDE,
     '  ONFACE0,ONFACE1,OPFILE,Single_point,TRANSFORM
      DATA FIELD/'Fibre','Imbrc','Sheet'/
      DATA NU1/2,4,7/
      DATA TOL/1.0D-5/

      CALL ENTERS('EVFIBR',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        WRITE(CHAR3(1:11),'(E11.4)') FOCUS
        CALL STRING_TRIM(CHAR3,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM evaluate fibre<;FILENAME> x_prolate_spheroidal
C###  Parameter:        <MU#[0.0],THETA#[0.0]>
C###  Parameter:        <Radius_xi=XI#[0.0]{>=0.0,<=1.0}>
C###  Parameter:        <s3=S3#[101 pts]>
C###  Parameter:        <(degrees/radians)[degrees]>
C###    Specify whether THETA is in degrees or radians.
C###  Parameter:        <region (#s/all)[1]>
C###    Specify the region numbers to evaluate.
C###  Description:
C###

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> x_prolate_spheroidal'
        OP_STRING(2)=BLANK(1:15)//'<MU#[0.0],THETA#[0.0]>'
        OP_STRING(3)=BLANK(1:15)//'<Radius_xi=XI#[0.0]{>=0.0,<=1.0}>'
        OP_STRING(4)=BLANK(1:15)//'<s3=S3#[0.0]>'
        OP_STRING(5)=BLANK(1:15)//'<(degrees/radians)[degrees]>'
        OP_STRING(6)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM evaluate fibre<;FILENAME> x_cylindrical_polar
C###  Parameter:        <X#[0.0],THETA#[0.0]>
C###  Parameter:        <s3=S3#[101 pts]>
C###  Parameter:        <(degrees/radians)[degrees]>
C###    Specify whether THETA is in degrees or radians.
C###  Parameter:        <region (#s/all)[1]>
C###    Specify the region numbers to evaluate.
C###  Description:
C###    Heart specific fibre evaluation.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> x_cylindrical_polar'
        OP_STRING(2)=BLANK(1:15)//'<X#[0.0],THETA#[0.0]>'
        OP_STRING(3)=BLANK(1:15)//'<s3=S3#[0.0]>'
        OP_STRING(4)=BLANK(1:15)//'<(degrees/radians)[degrees]>'
        OP_STRING(5)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM evaluate fibre x_rectangular_cart
C###  Parameter:        <X#[0.0],Y#[0.0],Z#[0.0]>
C###  Parameter:        <s3=S3#[0.0]>
C###  Parameter:        <focus=FOCUS#[1.0]>
C###  Parameter:        <(degrees/radians)[degrees]>
C###    Specify  degrees or radians.
C###  Parameter:        <region (#s/all)[1]>
C###    Specify the region numbers to evaluate.
C###  Description:
C###    Heart specific fibre evaluation.

        OP_STRING(1)=STRING(1:IEND)//' x_rectangular_cart'
        OP_STRING(2)=BLANK(1:15)//'<X#[0.0],Y#[0.0],Z#[0.0]>'
        OP_STRING(3)=BLANK(1:15)//'<s3=S3#[0.0]>'
        OP_STRING(4)=BLANK(1:15)
     '    //'<focus=FOCUS#['//CHAR3(IBEG1:IEND1)//']>'
        OP_STRING(5)=BLANK(1:15)//'<(degrees/radians)[degrees]>'
        OP_STRING(6)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM evaluate fibre data
C###  Parameter:        <rectangular_cart/
C###  Parameter:         prolate_spheroidal[prolate_spheroidal]>
C###    Specify Rectangular cartesian or Prolate spheroidal
C###  Parameter:        <focus=FOCUS#[1.0]>
C###  Parameter:        <(degrees/radians)[degrees]>
C###    Specify  degrees or radians.
C###  Parameter:        <region (#s/all)[1]>
C###    Specify the region numbers to evaluate.
C###  Description:
C###    Heart specific fibre evaluation.

        OP_STRING(1)=STRING(1:IEND)//' data'
        OP_STRING(2)=BLANK(1:15)//'<rectangular_cart/'
     '    //'prolate_spheroidal[prolate_spheroidal]>'
        OP_STRING(3)=BLANK(1:15)
     '    //'<focus=FOCUS#['//CHAR3(IBEG1:IEND1)//']>'
        OP_STRING(4)=BLANK(1:15)//'<(degrees/radians)[degrees]>'
        OP_STRING(5)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM evaluate fibre angle
C###  Parameter:        <from data_vector[data_vector]>
C###     Uses the data arrays ZD
C###  Parameter:        <wrt_xi #[1]>
C###    Evalutate wrt a specified xi direction
C###  Parameter:        <correct>
C###    Angle is corrected to lie between 0 and 360 Degrees
C###  Parameter:        <(degrees/radians)[degrees]>
C###    Specify  degrees or radians.
C###  Description:
C###    Evaluates the fibre angle with respect to a xi coordinate. If
C###    the fibre angle is from a data_vector then it is assumed that
C###    the fibre is given as a a vector in space in the NJT+1 - NJT+NJT
C###    positions of ZD. The fibre angle is the evaluated as the angle
C###    between this vector and the xi_one/xi_two vector. If correct
C###    is specified the angle is corrected so that it lies between
C###    0 and 2 pi radians (0 to 360 degrees) as opposed to -pi to pi
C###    radians (-180 to 180 degress). The result overwrites
C###    ZD(NJT+1,nd) with the fibre angle.

        OP_STRING(1)=STRING(1:IEND)//' angle'
        OP_STRING(2)=BLANK(1:15)//'<from data_vector[data_vector]>'
        OP_STRING(3)=BLANK(1:15)//'<wrt_xi #[1]>'
        OP_STRING(3)=BLANK(1:15)//'<correct>'
        OP_STRING(4)=BLANK(1:15)//'<(degrees/radians)[degrees]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM evaluate fibre xi
C###  Parameter:        <XI_1#[0.0]{>=0.0,<=1.0},XI_2#[0.0]{>=0.0,<=1.0},XI_3#[0.0]{>=0.0,<=1.0}>
C###  Parameter:        <in ELEMENT#[1]>
C###    Specify the element to be evaluated
C###  Parameter:        <(degrees/radians)[degrees]>
C###    Specify degrees or radians.
C###  Parameter:        <region #[1]>
C###    Specify the region number to evaluate.
C###  Description:
C###    Heart specific fibre evaluation.

        OP_STRING(1)=STRING(1:IEND)//' xi <XI_1#[0.0]{>=0.0,<=1.0},'
     '    //'XI_2#[0.0]{>=0.0,<=1.0},XI_3#[0.0]{>=0.0,<=1.0>'
        OP_STRING(2)=BLANK(1:15)//'<in ELEMENT#[1]>'
        OP_STRING(3)=BLANK(1:15)//'<(degrees/radians)[degrees]>'
        OP_STRING(4)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','EVFIBR',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN
          FILE=COQU(noco,1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          OPFILE=.TRUE.
          IOFI=IOFILE1
          IOFI2=IOFILE2
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.ipfibr','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
          IOFI=IOOP
          IOFI2=IOOP
        ENDIF

        IF(CBBREV(CO,'REGION',1,noco+1,NTCO,N3CO)) THEN
          nr=IFROMC(CO(N3CO+1))
        ELSE
          nr=1
        ENDIF

        IF(CBBREV(CO,'X_RECTANGULAR_CART',3,noco+1,NTCO,N3CO)) THEN

          CALL ASSERT(ITYP10(nr).EQ.4,
     '      '>>Valid only in prolate sph. coords',ERROR,*9999)
          CALL ASSERT(NJ_LOC(NJL_FIBR,0,nr).GT.0,
     '      '>>Fibre field not defined',ERROR,*9999)

          DO nj=1,NJT
            XIPOS(nj)=0.d0
            XPOS(nj)=0.d0
          ENDDO

C ***     Rect. cart. reference coords are specified
          CALL PARSRL(CO(N3CO+1),3,NTLIST,XRECT,ERROR,*9999)
          IF(CBBREV(CO,'S3',2,noco+2,NTCO,N3CO)) THEN
            S3=RFROMC(CO(N3CO+1))
            CALL ASSERT(S3.GE.0.d0.AND.S3.LE.1.d0,'0 <= S3 <= 1 !',
     '        ERROR,*9999)
          ELSE
            S3=0.d0
          ENDIF
          IF(CBBREV(CO,'FOCUS',1,noco+2,NTCO,N3CO)) THEN
            FOCUSDATA=RFROMC(CO(N3CO+1))
            CALL ASSERT(FOCUS.GT.0.d0,'FOCUS>0!',ERROR,*9999)
          ELSE
            FOCUSDATA=FOCUS
          ENDIF
          CURVIL=.FALSE.
          TRANSFORM=.TRUE.
          DATAPTS=.FALSE.
          NDMAX=1

        ELSE IF(CBBREV(CO,'DATA',1,noco+1,NTCO,N3CO)) THEN

          CALL ASSERT(ITYP10(nr).EQ.4,
     '      '>>Valid only in prolate sph. coords',ERROR,*9999)
          CALL ASSERT(NJ_LOC(NJL_FIBR,0,nr).GT.0,
     '      '>>Fibre field not defined',ERROR,*9999)

          DO nj=1,NJT
            XIPOS(nj)=0.d0
            XPOS(nj)=0.d0
          ENDDO

C ***     Data points have been read in from an IPDATA file with
C ***     reference coords in first 3 positions ie) ZD(1..3,nd),
C ***     and S3 coord in 4th position ie) WD(1,nd).
          CALL ASSERT(CALL_DATA,'>>no data points - use DEFINE DATA',
     '      ERROR,*9999)
          DO nd=1,NDMAX
            S3=WD(1,nd)
            CALL ASSERT(S3.GE.0.d0.AND.S3.LE.1.d0,'0 <= S3 <= 1 !',
     '        ERROR,*9999)
          ENDDO
          IF(CBBREV(CO,'RECTANGULAR_CART',1,noco+2,NTCO,N3CO)) THEN
            CURVIL=.FALSE.
          ELSE
            CURVIL=.TRUE.
          ENDIF
          IF(CBBREV(CO,'FOCUS',1,noco+2,NTCO,N3CO)) THEN
            FOCUSDATA=RFROMC(CO(N3CO+1))
            CALL ASSERT(FOCUS.GT.0.d0,'FOCUS>0!',ERROR,*9999)
          ELSE
            FOCUSDATA=FOCUS
          ENDIF
          TRANSFORM=.TRUE.
          DATAPTS=.TRUE.
          NDMAX=NDT

        ELSE IF(CBBREV(CO,'XI',2,noco+1,NTCO,N3CO)) THEN

          CALL ASSERT(ITYP10(nr).EQ.4,
     '      '>>Valid only in prolate sph. coords',ERROR,*9999)
          CALL ASSERT(NJ_LOC(NJL_FIBR,0,nr).GT.0,
     '      '>>Fibre field not defined',ERROR,*9999)

          DO nj=1,NJT
            XIPOS(nj)=0.d0
            XPOS(nj)=0.d0
          ENDDO

C ***     Xi material coords are specified
          CALL PARSRL(CO(N3CO+1),3,NTLIST,XIPOS,ERROR,*9999)
          DO ni=1,3
            WRITE(CHAR1(1:1),'(I1)') ni
            CALL ASSERT(XIPOS(ni).GE.0.d0.AND.XIPOS(ni).LE.1.d0,
     '        '0 <= Xi'//CHAR1(1:1)//' <= 1 !',ERROR,*9999)
          ENDDO
          IF(CBBREV(CO,'IN',1,noco+1,NTCO,N3CO)) THEN
            NE_EVFIBR=IFROMC(CO(N3CO+1))
          ELSE
            NE_EVFIBR=1
          ENDIF
          CURVIL    =.FALSE.
          TRANSFORM =.FALSE.
          DATAPTS   =.FALSE.
          NDMAX=1

        ELSE IF(CBBREV(CO,'X_PROLATE_SPHEROIDAL',3,noco+1,NTCO,N3CO))
     '    THEN

          CALL ASSERT(ITYP10(nr).EQ.4,
     '      '>>Valid only in prolate sph. coords',ERROR,*9999)
          CALL ASSERT(NJ_LOC(NJL_FIBR,0,nr).GT.0,
     '      '>>Fibre field not defined',ERROR,*9999)

          DO nj=1,NJT
            XIPOS(nj)=0.d0
            XPOS(nj)=0.d0
          ENDDO

C ***     Prolate spheroidal reference coords are specified
          CALL PARSRL(CO(N3CO+1),2,NTLIST,XPOS(2),ERROR,*9999)
          IF(CBBREV(CO,'RADIUS_XI',2,noco+1,NTCO,N3CO)) THEN
            Radius_xi=RFROMC(CO(N3CO+1))
          ELSE
            Radius_xi=0.d0
          ENDIF
          IF(CBBREV(CO,'S3',2,noco+1,NTCO,N3CO)) THEN
            S3=RFROMC(CO(N3CO+1))
            CALL ASSERT(S3.GE.0.d0.AND.S3.LE.1.d0,'0 <= S3 <= 1 !',
     '        ERROR,*9999)
            Single_point=.TRUE.
          ELSE
            Single_point=.FALSE.
          ENDIF
          COORDS    ='PROLATE'
          CURVIL    =.TRUE.
          TRANSFORM =.TRUE.
          DATAPTS   =.FALSE.
          IF(Single_point) THEN
            NDMAX=1
          ELSE
            NDMAX=101
          ENDIF

        ELSE IF(CBBREV(CO,'X_CYLINDRICAL_POLAR',3,noco+1,NTCO,N3CO))
     '    THEN

          CALL ASSERT(ITYP10(nr).EQ.4,
     '      '>>Valid only in prolate sph. coords',ERROR,*9999)
          CALL ASSERT(NJ_LOC(NJL_FIBR,0,nr).GT.0,
     '      '>>Fibre field not defined',ERROR,*9999)

          DO nj=1,NJT
            XIPOS(nj)=0.d0
            XPOS(nj)=0.d0
          ENDDO

C ***     Cylindrical polar reference coords are specified
          CALL PARSRL(CO(N3CO+1),2,NTLIST,XPOS(2),ERROR,*9999)
          IF(CBBREV(CO,'S3',2,noco+2,NTCO,N3CO)) THEN
            S3=RFROMC(CO(N3CO+1))
            CALL ASSERT(S3.GE.0.d0.AND.S3.LE.1.d0,'0 <= S3 <= 1 !',
     '        ERROR,*9999)
            Single_point=.TRUE.
          ELSE
            Single_point=.FALSE.
          ENDIF
          COORDS    ='CYL_POLAR'
          CURVIL    =.TRUE.
          TRANSFORM =.TRUE.
          DATAPTS   =.FALSE.
          IF(Single_point) THEN
            NDMAX=1
          ELSE
            NDMAX=101
          ENDIF

        ELSE IF(CBBREV(CO,'ANGLE',3,noco+1,NTCO,N3CO)) THEN

          TYPE='ANGLE'

          IF(CBBREV(CO,'FROM',2,noco+1,NTCO,N3CO)) THEN
            INPUT_TYPE=CO(N3CO+1)(1:11)
            CALL ASSERT(INPUT_TYPE(1:11).EQ.'DATA_VECTOR',
     '        '>>Unknown from type',ERROR,*9999)
          ELSE
            INPUT_TYPE='DATA_VECTOR'
          ENDIF
          IF(CBBREV(CO,'WRT_XI',2,noco+1,NTCO,N3CO)) THEN
            WRT_XI=IFROMC(CO(N3CO+1))
            CALL ASSERT(WRT_XI.GE.1.AND.WRT_XI.LE.2,
     '        '>>Invalid wrt_xi number',ERROR,*9999)
          ELSE
            WRT_XI=1
          ENDIF
          IF(CBBREV(CO,'CORRECT',2,noco+1,NTCO,N3CO)) THEN
            CORRECT=.TRUE.
          ELSE
            CORRECT=.FALSE.
          ENDIF

        ELSE
          CALL STRING_TRIM(STRING,IBEG,IEND)
          CALL STRING_TRIM(CO(noco),IBEG1,IEND1)
          CTEMP=CO(noco)(IBEG1:IEND1)
          CALL STRING_TRIM(CTEMP,IBEG1,IEND1)
          STRING=STRING(1:IEND)//' '//CTEMP(IBEG1:IEND1)
          CO(noco+1)='?'
          GO TO 1
        ENDIF

        IF(CBBREV(CO,'RADIANS',1,noco+1,NTCO,N3CO)) THEN
          MODE='RADIANS'
        ELSE
          MODE='DEGREES'
        ENDIF

        IF(TYPE(1:5).EQ.'ANGLE') THEN

          IF(INPUT_TYPE(1:11).EQ.'DATA_VECTOR') THEN

            CALL ASSERT(CALC_XI,'>>Define Xi positions first',
     '        ERROR,*9999)
            nr=1
            DO nd=1,NDT
              IF(LD(nd).NE.0) THEN
                ne=LD(nd)
                CALL ASSERT(NIT(NBJ(1,ne)).EQ.2,
     '            '>>Invalid number of xi directions',ERROR,*9999)
                CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '            NPNE(1,1,ne),nr,NVJE(1,1,1,ne),SE(1,1,ne),XA,
     '            XE,XP,ERROR,*9999)
                NORM1=0.0d0
                NORM2=0.0d0
                NORM3=0.0d0
                DO nj=1,NJT
                  nb=NBJ(nj,ne)
                  DXDXI1(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '              nb,2,XID(1,nd),XE(1,nj))
                  NORM1=NORM1+DXDXI1(nj)**2
                  DXDXI2(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '              nb,4,XID(1,nd),XE(1,nj))
                  NORM2=NORM2+DXDXI2(nj)**2
                  FIBRE_VECT(nj)=ZD(NJT+nj,nd)
                  NORM3=NORM3+FIBRE_VECT(nj)**2
                ENDDO !nj
                NORM1=DSQRT(NORM1)
                NORM2=DSQRT(NORM2)
                NORM3=DSQRT(NORM3)
                FDOTXI1=0.0d0
                FDOTXI2=0.0d0
                DO nj=1,NJT
                  DXDXI1(nj)=DXDXI1(nj)/NORM1
                  DXDXI2(nj)=DXDXI2(nj)/NORM2
                  FIBRE_VECT(nj)=FIBRE_VECT(nj)/NORM3
                  FDOTXI1=FDOTXI1+FIBRE_VECT(nj)*DXDXI1(nj)
                  FDOTXI2=FDOTXI2+FIBRE_VECT(nj)*DXDXI2(nj)
                ENDDO !nj
                IF(WRT_XI.EQ.1) THEN
                  ALPHA=DACOS(FDOTXI1)
                  IF(FDOTXI2.LT.0.0d0) ALPHA=-ALPHA
                ELSE
                  ALPHA=DACOS(FDOTXI2)
                  IF(FDOTXI1.LT.0.0d0) ALPHA=-ALPHA
                ENDIF
                IF(CORRECT.AND.ALPHA.LT.0.0d0) ALPHA=ALPHA+2.0d0*PI
                IF(MODE(1:7).EQ.'DEGREES') THEN
                  ALPHA=ALPHA*180.0d0/PI
                ENDIF
                ZD(NJT+1,nd)=ALPHA
              ENDIF
            ENDDO !nd

          ENDIF

        ELSE

          IF(MODE(1:7).EQ.'DEGREES') THEN !degrees -> radians
            IF(COORDS(1:7).EQ.'PROLATE') THEN
              XPOS(2)=XPOS(2)*PI/180.d0 !mu
              XPOS(3)=XPOS(3)*PI/180.d0 !theta
C GMH 13/2/97 Unused
C            Radius=Radius*PI/180.d0   !radius in mu,theta space
            ELSE IF(COORDS(1:9).EQ.'CYL_POLAR') THEN
              XPOS(3)=XPOS(3)*PI/180.d0 !theta
            ENDIF
          ENDIF

          IF(DATAPTS.AND.OPFILE) THEN
C ***     Write header lines of IPFIBR file
            WRITE(OP_STRING,FMT='(A)') 'CMISS Version '//CMISS
     '        //' IPFIBR File Version  1'
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,FMT='(A)') 'Heading: '//HEADING
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,FMT='(1X)')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            FORMAT='('' Specify whether [1]: '''//
     '        '/''   (1) fibre angle in Xi_1&2 plane'''//
     '        '/''   (2) fibre and imbrication angles'''//
     '        '/''   (3) '''//
     '        '/$,''    '',I1)'
            IF(NJ_LOC(NJL_FIBR,0,nr).EQ.2) THEN
              WRITE(OP_STRING,FMT=FORMAT) 2
            ELSE
              WRITE(OP_STRING,FMT=FORMAT) 1
            ENDIF
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            IF(NJ_LOC(NJL_FIBR,0,nr).EQ.2) THEN
              FORMAT='('' Specify how fibre angle is defined [1]: '''//
     '          '/''   (1) wrt Xi(2) coordinate'''//
     '          '/''   (2) wrt Xi(3) coordinate'''//
     '          '/$,''    '',I1)'
            ELSE
              FORMAT='('' Specify how fibre angle is defined [1]: '''//
     '          '/''   (1) wrt Xi(1) coordinate'''//
     '          '/''   (2) wrt Xi(2) coordinate'''//
     '          '/$,''    '',I1)'
            ENDIF
            WRITE(OP_STRING,FMT=FORMAT) JTYP12
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            FORMAT='('' Specify whether angles entered in [1]: '''//
     '        '/''   (1) degrees'''//
     '        '/''   (2) radians'''//
     '        '/$,''    '',I1)'
            IF(MODE.EQ.'DEGREES') THEN
              WRITE(OP_STRING,FMT=FORMAT) 1
            ELSE
              WRITE(OP_STRING,FMT=FORMAT) 2
            ENDIF
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            FORMAT='($,'' Is the basis function for the fibre angle '
     '        //'element dependent [N]? '',A)'
            WRITE(OP_STRING,FMT=FORMAT) 'Y'
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              WRITE(CHAR2,'(I1)') NBJ(NJ_LOC(NJL_FIBR,1,nr),ne)
              WRITE(CHAR6,'(I5)') ne
              FORMAT='($,'' The basis function type number for the '
     '          //'fibre angle in element '//CHAR6(1:5)//' is ['
     '          //CHAR2(1:1)//']: '','//'I1)'
              WRITE(OP_STRING,FMT=FORMAT) NBJ(NJ_LOC(NJL_FIBR,1,nr),ne)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO !noelem
            IF(NJ_LOC(NJL_FIBR,0,nr).EQ.2) THEN
              FORMAT='($,'' Is the basis function for the imbrication '
     '          //'angle element dependent [N]? '',A)'
              WRITE(OP_STRING,FMT=FORMAT) 'Y'
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                WRITE(CHAR2,'(I1)') NBJ(NJ_LOC(NJL_FIBR,2,nr),ne)
                WRITE(CHAR6,'(I5)') ne
                FORMAT='($,'' The basis function type number for the '
     '            //'imbrication angle in element '//CHAR6(1:5)
     '            //' is ['//CHAR2(1:1)//']: '',I1)'
                WRITE(OP_STRING,FMT=FORMAT)
     '            NBJ(NJ_LOC(NJL_FIBR,2,nr),ne)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDDO !noelem
            ENDIF
            WRITE(CHAR2(1:4),'(I4)') NDT
            FORMAT='($,'' The number of nodes is ['//CHAR2(1:4)//
     '        ']: '',I4)'
            WRITE(OP_STRING,FMT=FORMAT) NDT
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            FORMAT='($,'' Do you want prompting for different versions '
     '        //'of the fibre field [N]? '',A)'
            WRITE(OP_STRING,FMT=FORMAT) 'N'
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(CHAR1(1:1),'(I1)')NKJ(NJ_LOC(NJL_FIBR,1,nr),
     '        NPNODE(1,nr))-1
            FORMAT='($,'' The number of derivatives for the fibre angle'
     '        //' is ['//CHAR1(1:1)//']: '',I1)'
            WRITE(OP_STRING,FMT=FORMAT)
     '        NKJ(NJ_LOC(NJL_FIBR,1,nr),NPNODE(1,nr))-1
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            IF(NJ_LOC(NJL_FIBR,0,nr).EQ.2) THEN
              FORMAT='($,'' Do you want prompting for different '
     '          //'versions of the imbrication field [N]? '',A)'
              WRITE(OP_STRING,FMT=FORMAT) 'N'
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(CHAR1(1:1),'(I1)')NKJ(NJ_LOC(NJL_FIBR,2,nr),
     '          NPNODE(1,nr))-1
              FORMAT='($,'' The number of derivatives for the '
     '          //'imbrication angle is ['//CHAR1(1:1)//']: '',I1)'
              WRITE(OP_STRING,FMT=FORMAT)
     '          NKJ(NJ_LOC(NJL_FIBR,1,nr),NPNODE(1,nr))-1
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF

            IF(NJ_LOC(NJL_FIBR,0,nr).GE.3) THEN
C ***       Write header lines of IPSHEE file
              WRITE(OP_STRING,FMT='(A)') 'CMISS Version '//CMISS
     '          //' IPSHEE File Version  1'
              CALL WRITES(IOFI2,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,FMT='(A)') 'Heading: '//HEADING
              CALL WRITES(IOFI2,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,FMT='(1X)')
              CALL WRITES(IOFI2,OP_STRING,ERROR,*9999)
              FORMAT='('' Specify whether angles entered in [1]: '''//
     '          '/''   (1) degrees'''//
     '          '/''   (2) radians'''//
     '          '/$,''    '',I1)'
              IF(MODE.EQ.'DEGREES') THEN
                WRITE(OP_STRING,FMT=FORMAT) 1
              ELSE
                WRITE(OP_STRING,FMT=FORMAT) 2
              ENDIF
              CALL WRITES(IOFI2,OP_STRING,ERROR,*9999)
              FORMAT='($,'' Is the basis function for the sheet angle '
     '          //'element dependent [N]? '',A)'
              WRITE(OP_STRING,FMT=FORMAT) 'Y'
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                WRITE(CHAR2,'(I1)') NBJ(NJ_LOC(NJL_FIBR,3,nr),ne)
                WRITE(CHAR6,'(I5)') ne
                FORMAT='($,'' The basis function type number for the '
     '            //'sheet angle in element '//CHAR6(1:5)//' is ['
     '            //CHAR2(1:1)//']: '','//'I1)'
                WRITE(OP_STRING,FMT=FORMAT)
     '            NBJ(NJ_LOC(NJL_FIBR,3,nr),ne)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDDO !noelem
              WRITE(CHAR2(1:4),'(I4)') NDT
              FORMAT='($,'' The number of nodes is ['//CHAR2(1:4)//
     '          ']: '',I4)'
              WRITE(OP_STRING,FMT=FORMAT) NDT
              CALL WRITES(IOFI2,OP_STRING,ERROR,*9999)
              FORMAT='($,'' Do you want prompting for different '
     '          //'versions of the sheet field [N]? '',A)'
              WRITE(OP_STRING,FMT=FORMAT) 'N'
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(CHAR1(1:1),'(I1)')NKJ(NJ_LOC(NJL_FIBR,3,nr),
     '          NPNODE(1,nr))-1
              FORMAT='($,'' The number of derivatives for the sheet '
     '          //'field is ['//CHAR1(1:1)//']: '',I1)'
              WRITE(OP_STRING,FMT=FORMAT)
     '          NKJ(NJ_LOC(NJL_FIBR,3,nr),NPNODE(1,nr))-1
              CALL WRITES(IOFI2,OP_STRING,ERROR,*9999)
            ENDIF

          ENDIF !DATAPTS.AND.OPFILE

          DO nd=1,NDMAX
            IF(.NOT.DATAPTS.AND.(COORDS(1:7).EQ.'PROLATE'.OR.
     '        COORDS(1:9).EQ.'CYL_POLAR').
     '        AND..NOT.Single_point) THEN
              S3=DBLE(nd-1)/1.d2
            ENDIF
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
              WRITE(OP_STRING,'(//'' Data pt '',I3,'' S3='',F6.2)')
     '          nd,S3
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
            ENDIF

            IF(TRANSFORM) THEN
C ***       Reference coords have been given so must transform them into
C ***       Xi1 and Xi2 material coords. Note: Xi3 has been set above.

              IF(CURVIL) THEN
C ***         Prolate sph or cyl polar coords have been entered
                IF(DATAPTS) THEN
                  XPOS(2)=ZD(2,nd) !mu
                  XPOS(3)=ZD(3,nd) !theta
                  IF(MODE(1:7).EQ.'DEGREES') THEN !degrees -> radians
                    IF(COORDS(1:7).EQ.'PROLATE') THEN
                      XPOS(2)=XPOS(2)*PI/180.d0 !mu
                      XPOS(3)=XPOS(3)*PI/180.d0 !theta
C GMH 13/2/97 Unused
C                    Radius=Radius*PI/180.d0   !radius in mu,theta space
                    ELSE IF(COORDS(1:9).EQ.'CYL_POLAR') THEN
                      XPOS(3)=XPOS(3)*PI/180.d0 !theta
                    ENDIF
                  ENDIF
                ENDIF !datapts
              ELSE
C ***         Rect. cart. coords have been entered so must transform
C ***         them to prolate spheroidal coords with a focus
C ***         applicable to the entered coords.
                FOCUSTEMP=FOCUS
                FOCUS=FOCUSDATA
                IF(DATAPTS) THEN
                  CALL ZX(ITYP10(nr),ZD(1,nd),XPOS)
                ELSE
                  CALL ZX(ITYP10(nr),XRECT,XPOS)
                ENDIF
                FOCUS=FOCUSTEMP
              ENDIF !curvil

              NE_EVFIBR=0
              NEENDO=0
              NUMELEM(0)=0
              DO ne=1,MAXWALLELEM
                NUMELEM(ne)=0
                INDEXELEM(ne)=0
                DO J=1,2
                  XG1LIST(ne,J)=0.d0
                  XIF1LIST(ne,J)=0.d0
                  XIF2LIST(ne,J)=0.d0
                ENDDO
                BASALLIST(ne)=.FALSE.
              ENDDO !ne
              DO J=1,2
                XIFACE1(J)=0.d0
                XIFACE2(J)=0.d0
              ENDDO
              NINDEX=0
              XG1MIN=RMAX

C ***     Find neighbouring elements
              CALL NENXI(NBJ,NEELEM,NENP,NNB,NPNE,NXI,ERROR,*9999)

C ***     Find all elements through which const mu,theta or x,theta
C ***     line passes.
C ***     Transform theta and mu into Xi1 and Xi2 for those elements.
C ***     Calc Xi3 based upon S3 (proportion through the wall)
              NJ1=1 !lambda
              NJ2=2 !mu
              NJ3=3 !theta
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                  WRITE(OP_STRING,'(/'' Element '',I3)') ne
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
                ENDIF
                CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '            NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '            SE(1,1,ne),XA,XE,XP,ERROR,*9999)

                ONFACE0=.FALSE.
                ONFACE1=.FALSE.
                BASAL=.FALSE.
                DO iface=0,1 !loop over xi3=0, xi3=1 faces
                  nb=NBJ(NJ2,ne) !nb for mu coord
                  NKTB=NKT(0,nb)
                  NS1=1        +iface*4*NKTB
                  NS2=1+NKTB   +iface*4*NKTB
                  NS3=1+2*NKTB +iface*4*NKTB
                  NS4=1+3*NKTB +iface*4*NKTB
                  A1=XE(NS2,NJ2)-XE(NS1,NJ2)
                  B1=XE(NS3,NJ2)-XE(NS1,NJ2)
                  C1=XE(NS1,NJ2)-XE(NS2,NJ2)-XE(NS3,NJ2)+XE(NS4,NJ2)
                  nb=NBJ(NJ3,ne) !nb for theta coord
                  NKTB=NKT(0,nb)
                  NS1=1        +iface*4*NKTB
                  NS2=1+NKTB   +iface*4*NKTB
                  NS3=1+2*NKTB +iface*4*NKTB
                  NS4=1+3*NKTB +iface*4*NKTB
                  A2=XE(NS2,NJ3)-XE(NS1,NJ3)
                  B2=XE(NS3,NJ3)-XE(NS1,NJ3)
                  C2=XE(NS1,NJ3)-XE(NS2,NJ3)-XE(NS3,NJ3)+XE(NS4,NJ3)

                  IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                    WRITE(OP_STRING,'('' A1,B1,C1='',3E11.3)') A1,B1,C1
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    WRITE(OP_STRING,'('' A2,B2,C2='',3E11.3)') A2,B2,C2
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                  ENDIF
                  ALFA=A2*C1-A1*C2
                  EXCLUDE=.TRUE.
                  IF(COORDS(1:7).EQ.'PROLATE') THEN
                    MU=XPOS(2)
                  ELSE IF(COORDS(1:9).EQ.'CYL_POLAR') THEN
! calc mu from x,theta using what r?
                  ENDIF
                  THETA=XPOS(3)
                  IF(ITYP10(nr).EQ.4) THEN !check if lies within element
                    THETAMAX=DMAX1(XE(NS1,3),XE(NS3,3))
                    THETAMIN=DMIN1(XE(NS2,3),XE(NS4,3))
                    IF(THETA.LE.THETAMAX.AND.THETA.GE.THETAMIN) THEN
                      EXCLUDE=.FALSE.
                      D2=THETA-XE(1,NJ3)
                    ELSE IF((THETA+2.d0*PI).LE.THETAMAX
     '                  .AND.(THETA+2.d0*PI).GE.THETAMIN) THEN
                      EXCLUDE=.FALSE.
                      D2=THETA+2.d0*PI-XE(1,NJ3)
                    ELSE IF((THETA-2.d0*PI).LE.THETAMAX
     '                  .AND.(THETA-2.d0*PI).GE.THETAMIN) THEN
                      EXCLUDE=.FALSE.
                      D2=THETA-2.d0*PI-XE(1,NJ3)
                    ENDIF
                  ELSE
                    EXCLUDE=.FALSE.
                    D2=THETA-XE(1,NJ3)
                  ENDIF
                  IF(.NOT.EXCLUDE) THEN
                    D1=MU-XE(1,NJ2)
                    BETA=C2*D1-C1*D2+A2*B1-A1*B2
                    GAMA=B2*D1-B1*D2
                    DELTA=BETA*BETA-4.d0*ALFA*GAMA
                    IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                  call mp_setlock()
                      WRITE(OP_STRING,
     '                  '('' D1='',E11.3,'' D2='',E11.3,'
     '                  //''' ALFA='',E11.3,'' BETA='',E11.3,'
     '                  //''' GAMA='',E11.3,'' DELTA='',E11.3)')
     '                  D1,D2,ALFA,BETA,GAMA,DELTA
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                  call mp_unsetlock()
                    ENDIF !dop
                    IF(DABS(ALFA).GT.10D-6) THEN
                      IF(DELTA.GE.0.d0) THEN
                        XIPOS(1)=(-BETA+DSQRT(DELTA))/(2.d0*ALFA)
                        IF(XIPOS(1).LT.0.d0.OR.XIPOS(1).GT.1.d0) THEN
                          XIPOS(1)=(-BETA-DSQRT(DELTA))/(2.d0*ALFA)
                        ENDIF
                      ENDIF
                    ELSE IF(DABS(ALFA).LE.1.0D-6) THEN
                      IF(DABS(BETA).GT.1.0D-6) THEN
                        XIPOS(1)=-GAMA/BETA
                      ELSE
                        XIPOS(1)=-1.d0
                      ENDIF
                    ENDIF
                    IF(XIPOS(1).GE.0.d0.AND.XIPOS(1).LE.1.d0) THEN
                      DENOM1=B1+C1*XIPOS(1)
                      IF(DABS(DENOM1).GT.1.0D-6) THEN
                        XIPOS(2)=(D1-A1*XIPOS(1))/DENOM1
                      ELSE
                        DENOM2=B2+C2*XIPOS(1)
                        IF(DABS(DENOM2).GT.1.0D-6) THEN
                          XIPOS(2)=(D2-A2*XIPOS(1))/DENOM2
                        ELSE
                          WRITE(OP_STRING,'('' Xi2 cannot be defined '
     '                      //'in elem. '',I5)') ne
                          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                        ENDIF
                      ENDIF
                      IF(XIPOS(2).GT.1.d0.AND.NXI(2,1,ne).EQ.0) THEN
C ***                 If desired pt is above basal faces of mesh then
C ***                 use fibres angle on basal face (ie set Xi2=1)
                        XIPOS(2)=1.d0
                        BASAL=.TRUE. !set true to flag a message later
                      ENDIF
                      IF(XIPOS(2).GE.0.d0.AND.XIPOS(2).LE.1.d0) THEN
                        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                      call mp_setlock()
                          WRITE(OP_STRING,'('' Xi:'',3E11.3,'
     '                      //''' (linear calc.)'')') (XIPOS(ni),ni=1,3)
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                      call mp_unsetlock()
                        ENDIF !dop
                        nb=NBJ(NJ2,ne)
                        IF(NKT(0,nb).GT.1) THEN !update Xi(2)
                          XIPOS(2)=0.d0
                          X1=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '                      nb,1,XIPOS,XE(1,NJ2))
                          XIPOS(2)=1.d0
                          X2=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '                      nb,1,XIPOS,XE(1,NJ2))
                          IF(DABS(X2-X1).GT.1.0D-6) THEN !PJH 13APR91
                            XIPOS(2)=(MU-X1)/(X2-X1)
                          ELSE
                            WRITE(OP_STRING,'('' Warning!!!!'','
     '                        //'''X1=X2'')')
                            CALL WRITES(IOER,OP_STRING,ERROR,*9999)
                            XIPOS(2)=1.d0
                          ENDIF
                        ENDIF !nkt

C ***               Point lies on current face ie) 0<=xi1,xi2<=1
C ***               so interpolate to find XG for pt on the face
                        XIPOS(3)=DBLE(iface)
                        CALL XEXW(0,IBT,IDO,INP,NAN,NBJ(1,ne),nr,
     '                    XE,XG1(1,1,iface+1),XIPOS,ERROR,*9999)
                        IF(iface.EQ.0) THEN
                          ONFACE0=.TRUE. !on xi3=0 face
                          XIFACE1(1)=XIPOS(1) !store Xi1 for Xi3=0 face
                          XIFACE1(2)=XIPOS(2)
                        ELSE IF(iface.EQ.1) THEN
                          ONFACE1=.TRUE. !on xi3=1 face
                          XIFACE2(1)=XIPOS(1)
                          XIFACE2(2)=XIPOS(2)
                        ENDIF
                      ENDIF !0<xi(2)<1
                    ENDIF !0<xi(1)<1
                  ENDIF !.NOT.EXCLUDE
                ENDDO !iface

                IF(ONFACE0.AND.ONFACE1) THEN !pt in both xi3=0 & 1 faces
                  NUMELEM(0)=NUMELEM(0)+1
                  CALL ASSERT(NUMELEM(0).LE.MAXWALLELEM,
     '              '>>Increase MAXWALLELEM',ERROR,*9999)
                  NUMELEM(NUMELEM(0))=ne !add elem to list
                  DO J=1,2
                    XG1LIST(NUMELEM(0),J)=XG1(1,1,J) !store lambda@xi3=0,1
                    XIF1LIST(NUMELEM(0),J)=XIFACE1(J) !store face1 xi1,xi2
                    XIF2LIST(NUMELEM(0),J)=XIFACE2(J) !store face2 xi1,xi2
                  ENDDO
                  BASALLIST(NUMELEM(0))=BASAL !store if mu>mumax for elem
                  IF(XG1(1,1,1).LT.XG1MIN)THEN !find innermost element
                    XG1MIN=XG1(1,1,1) !by comparing lambdas on
                    NEENDO=ne !xi3=0 face
                    INDEXELEM(1)=NUMELEM(0) !index of neendo in list
                  ENDIF
                  IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                    FORMAT='('' Lies in element '',I3)'
                    WRITE(OP_STRING,FORMAT) ne
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    FORMAT='('' Element length='',E11.4)'
                    WRITE(OP_STRING,FORMAT) XG1(1,1,2)-XG1(1,1,1)
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                  ENDIF !dop

                ELSE IF(ONFACE0) THEN
                  IF(DATAPTS) THEN
                    FORMAT='('' Current data point:'',I3,'
     '                //'''; >>Mu and Theta lie in Xi3=0 face but not'
     '                //' Xi3=1 face for element '',I3)'
                    WRITE(OP_STRING,FORMAT) nd,ne
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ELSE
                    FORMAT='('' >>Mu and Theta lie in Xi3=0 face but '
     '                //'not Xi3=1 face for element '',I3)'
                    WRITE(OP_STRING,FORMAT) ne
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ENDIF
                  GO TO 100 !skip printing out any info

                ELSE IF(ONFACE1) THEN
                  IF(DATAPTS) THEN
                    FORMAT='('' Current data point:'',I3,'
     '                //'''; >>Mu and Theta lie in Xi3=1 face but not'
     '                //' Xi3=0 face for element '',I3)'
                    WRITE(OP_STRING,FORMAT) nd,ne
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ELSE
                    FORMAT='('' >>Mu and Theta lie in Xi3=1 face but '
     '                //'not Xi3=0 face for element '',I3)'
                    WRITE(OP_STRING,FORMAT) ne
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ENDIF
                  GO TO 100 !skip printing out any info
                ENDIF !onface0/1
              ENDDO !noelem

              IF(NUMELEM(0).GE.1) THEN !pt lies in an elem
C ***         Store positions of element #'s in list ordered from
C ***         endo to epi elements
                ne=NEENDO
                NINDEX=1
                DO WHILE(NXI(3,1,ne).ne.0)
                  ne=NXI(3,1,ne) !step to adjacent elem in +'ve xi3 dirn
                  NINDEX=NINDEX+1
                  DO noelem=1,NUMELEM(0)
                    IF(NUMELEM(noelem).eq.ne) INDEXELEM(NINDEX)=noelem
                  ENDDO
                ENDDO

C ***         Calc wall thickness by stepping to adjacent elems in
C ***         +'ve xi3 dirn from elem at endo
                WALL_THICK=0.d0
                DO I=1,NINDEX
                  INDEX=INDEXELEM(I)
                  WALL_THICK=WALL_THICK+XG1LIST(INDEX,2)-
     '              XG1LIST(INDEX,1)
                ENDDO

C ***         Find element number and xi3 coord for that elem
C ***         in which the specified point lies
C ***         based on wall_pos = S3*wall_thick + lamda(neendo,xi3=1)
                WALL_POS=S3*WALL_THICK+XG1LIST(INDEXELEM(1),1)
                DO I=1,NINDEX
                  INDEX=INDEXELEM(I)
                  IF((WALL_POS-XG1LIST(INDEX,1).GT.-TOL).AND.
     '              (WALL_POS-XG1LIST(INDEX,2).LT. TOL)) THEN
                    NE_EVFIBR=NUMELEM(INDEX)
C ***             Calc Xi3 as a proportion of dist through elem
                    XIPOS(3)=(WALL_POS-XG1LIST(INDEX,1))/
     '                (XG1LIST(INDEX,2)-XG1LIST(INDEX,1))
                    IF(XIPOS(3).GT.1.d0) XIPOS(3)=1.d0

C ***           Calc Xi1 (& Xi2) by LINEARLY interpolating between
C ***           the Xi1 (& Xi2) values on the Xi3=0 and Xi3=1 faces
                    XIPOS(1)=(1.d0-XIPOS(3))*XIF1LIST(INDEX,1) +
     '                XIPOS(3) *XIF2LIST(INDEX,1)
                    XIPOS(2)=(1.d0-XIPOS(3))*XIF1LIST(INDEX,2) +
     '                XIPOS(3) *XIF2LIST(INDEX,2)

                    IF(BASALLIST(INDEX))THEN !mu>mumax => on basal face
                      IF(DATAPTS) THEN
                        FORMAT='('' >>WARNING: Mu lies above basal end'
     '                    //' of mesh. Using fibre angle on basal face'
     '                    //' for data point '',I3,'' in element '',I3)'
                        WRITE(OP_STRING,FORMAT) nd,ne
                        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      ELSE
                        FORMAT='('' >>WARNING: Mu lies above basal end'
     '                    //' of mesh. Using fibre angle on basal '
     '                    //'face'')'
                        WRITE(OP_STRING,FORMAT)
                        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      ENDIF !datapts
                    ENDIF !basallist

C ***           Check that Lambda varies linearly in Xi3 so that the
C ***           linear interpolation used above for Xi1 & Xi2 is valid
                    nb=NBJ(NJ1,NE_EVFIBR)
                    CALL ASSERT(IBT(1,3,nb).EQ.1,
     '                '>>Valid only for Lagrange elements',ERROR,*9999)
                    CALL ASSERT(IBT(2,3,nb).EQ.1,
     '                '>>Lambda must vary linearly in Xi3',ERROR,*9999)
                  ENDIF !wall_pos etc
                ENDDO !i
              ENDIF !numelem(0).ge.1

              IF(NE_EVFIBR.EQ.0) THEN !point lies in no element
                IF(DATAPTS) THEN
                  FORMAT='('' Data point:'',I3,'
     '              //''', does not lie in any element.'')'
                  WRITE(OP_STRING,FORMAT) nd
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ELSE
                  FORMAT='('' Point does not lie in any element.'')'
                  WRITE(OP_STRING,FORMAT)
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDIF
                GO TO 100 !skip printing out any info
              ENDIF !ne_evfibr=0

              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
                WRITE(OP_STRING,'(/'' *****'')')
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' XPOS:'',3E11.3)')
     '            (XPOS(nj),nj=1,3)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' XIPOS:'',3E11.3)')
     '            (XIPOS(ni),ni=1,3)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
              ENDIF !dop
            ENDIF !transform

            CALL XPXE(NBJ(1,NE_EVFIBR),
     '        NKJE(1,1,1,NE_EVFIBR),NPF(1,1),NPNE(1,1,NE_EVFIBR),
     '        nr,NVJE(1,1,1,NE_EVFIBR),
     '        SE(1,1,NE_EVFIBR),XA,XE,XP,ERROR,*9999)
            CALL XEXW(0,IBT,IDO,INP,NAN,NBJ(1,NE_EVFIBR),nr,
     '        XE,XG,XIPOS,ERROR,*9999)

C ***   Calculate rate of change arc length wrt Xi to scale the
C ***   derivatives of angles wrt Xi
            SINHL=DSINH(XG(1,1))
            COSHL=DCOSH(XG(1,1))
            SINM=DSIN(XG(2,1))
            COSM=DCOS(XG(2,1))
            SINT=DSIN(XG(3,1))
            COST=DCOS(XG(3,1))
            DO ni=1,NIT(NBJ(1,NE_EVFIBR))
              DXDXI(1)=FOCUS*(SINHL*COSM*XG(1,NU1(ni))
     '          -COSHL*SINM*XG(2,NU1(ni)))
              DXDXI(2)=FOCUS*(COSHL*SINM*COST*XG(1,NU1(ni))
     '          +SINHL*COSM*COST*XG(2,NU1(ni))
     '          -SINHL*SINM*SINT*XG(3,NU1(ni)))
              DXDXI(3)=FOCUS*(COSHL*SINM*SINT*XG(1,NU1(ni))
     '          +SINHL*COSM*SINT*XG(2,NU1(ni))
     '          +SINHL*SINM*COST*XG(3,NU1(ni)))
              DSDXI(ni)=DSQRT(DXDXI(1)*DXDXI(1)+DXDXI(2)*DXDXI(2)
     '          +DXDXI(3)*DXDXI(3))
            ENDDO !ni

C ***   Determine scale factors for XG array
            CALL ASSERT(NJMAX.GE.NJM,'>>Increase NJMAX',ERROR,*9999)
            CALL ASSERT(NKMAX.GE.NKM,'>>Increase NKMAX',ERROR,*9999)
            DO nj=1,NJ_loc(0,0,nr)
              DO nk=1,NKJ(nj,NPNODE(1,nr))
                SCAFAC(nj,nk)=1.d0
                nu=IDO(nk,1,0,NBJ(nj,NE_EVFIBR))
                IF(nu.eq.2.or.nu.eq.4.or.nu.EQ.7) THEN
                  IF(nu.EQ.2) THEN
                    SCAFAC(nj,nk)=DSDXI(1)
                  ELSE IF(nu.EQ.4) THEN
                    SCAFAC(nj,nk)=DSDXI(2)
                  ELSE IF(nu.EQ.7) THEN
                    SCAFAC(nj,nk)=DSDXI(3)
                  ENDIF
                ELSE IF(nu.eq.6.or.nu.eq.9.or.nu.EQ.10) THEN
                  IF(nu.EQ.6) THEN
                    SCAFAC(nj,nk)=DSDXI(1)*DSDXI(2)
                  ELSE IF(nu.EQ.9) THEN
                    SCAFAC(nj,nk)=DSDXI(1)*DSDXI(3)
                  ELSE IF(nu.EQ.10) THEN
                    SCAFAC(nj,nk)=DSDXI(2)*DSDXI(3)
                  ENDIF !nu
                ENDIF !nu
              ENDDO !nk
            ENDDO !nj

            IF(DATAPTS) THEN
C ***     Write interpolated fibre angles in IPFIBR file format
              WRITE(CHAR2(1:4),'(I4)') nd
              FORMAT='(/$,'' Node number ['//CHAR2(1:4)//']: '',I4)'
              WRITE(OP_STRING,FMT=FORMAT) nd
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(CHAR1(1:1),'(I1)') NJ_LOC(NJL_FIBR,1,nr)
              DO nk=1,NKJ(NJ_LOC(NJL_FIBR,1,nr),NPNODE(1,nr))
                nu=IDO(nk,1,0,NBJ(NJ_LOC(NJL_FIBR,1,nr),NE_EVFIBR))
                IF(MODE(1:7).EQ.'DEGREES') THEN
                  WRITE(CHAR4,'(E12.5)') XG(NJ_LOC(NJL_FIBR,1,nr),nk)
     '              *180.d0/PI
                ELSE
                  WRITE(CHAR4,'(E12.5)') XG(NJ_LOC(NJL_FIBR,1,nr),nk)
                ENDIF
                IF(nu.EQ.1) THEN
                  FORMAT='($,'' The Xj('//CHAR1(1:1)//') coordinate'//
     '              ' is ['//CHAR4(1:12)//']: '',E12.5)'
                ELSE IF(nu.eq.2.or.nu.eq.4.or.nu.EQ.7) THEN
                  IF(nu.EQ.2) THEN
                    CHAR2='1'
                  ELSE IF(nu.EQ.4) THEN
                    CHAR2='2'
                  ELSE IF(nu.EQ.7) THEN
                    CHAR2='3'
                  ENDIF
                  FORMAT='($,'' The Xj('//CHAR1(1:1)//') derivative'//
     '              ' wrt s('//CHAR2(1:1)//') is ['//
     '              CHAR4(1:12)//']: '',E12.5)'
                ELSE IF(nu.eq.6.or.nu.eq.9.or.nu.EQ.10) THEN
                  IF(nu.EQ.6) THEN
                    CHAR2='1'
                    CHAR3='2'
                  ELSE IF(nu.EQ.9) THEN
                    CHAR2='1'
                    CHAR3='3'
                  ELSE IF(nu.EQ.10) THEN
                    CHAR2='2'
                    CHAR3='3'
                  ENDIF
                  FORMAT='($,'' The Xj('//CHAR1(1:1)//') derivative'//
     '              ' wrt s('//CHAR2(1:1)//') & s('//CHAR3(1:1)//
     '              ') is ['//CHAR4(1:12)//']: '',E12.5)'
                ENDIF !nu
                IF(MODE(1:7).EQ.'DEGREES') THEN
                  WRITE(OP_STRING,FMT=FORMAT)
     '              XG(NJ_LOC(NJL_FIBR,1,nr),nk)/
     '              SCAFAC(NJ_LOC(NJL_FIBR,1,nr),nk)*180.d0/PI
                ELSE
                  WRITE(OP_STRING,FMT=FORMAT)
     '              XG(NJ_LOC(NJL_FIBR,1,nr),nk)/
     '              SCAFAC(NJ_LOC(NJL_FIBR,1,nr),nk)
                ENDIF
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDDO !nk

              IF(NJ_LOC(NJL_FIBR,0,nr).EQ.2) THEN
                WRITE(CHAR1,'(I1)') NJ_LOC(NJL_FIBR,2,nr)
                DO nk=1,NKJ(NJ_LOC(NJL_FIBR,2,nr),NPNODE(1,nr))
                  nu=IDO(nk,1,0,NBJ(NJ_LOC(NJL_FIBR,2,nr),NE_EVFIBR))
                  IF(MODE(1:7).EQ.'DEGREES') THEN
                    WRITE(CHAR4,'(E12.5)') XG(NJ_LOC(NJL_FIBR,2,nr),nk)
     '                *180.d0/PI
                  ELSE
                    WRITE(CHAR4,'(E12.5)') XG(NJ_LOC(NJL_FIBR,2,nr),nk)
                  ENDIF
                  IF(nu.EQ.1) THEN
                    FORMAT='($,'' The Xj('//CHAR1(1:1)//') coordinate'//
     '                ' is ['//CHAR4(1:12)//']: '',E12.5)'
                  ELSE IF(nu.eq.2.or.nu.eq.4.or.nu.EQ.7) THEN
                    IF(nu.EQ.2) THEN
                      CHAR2='1'
                    ELSE IF(nu.EQ.4) THEN
                      CHAR2='2'
                    ELSE IF(nu.EQ.7) THEN
                      CHAR2='3'
                    ENDIF
                    FORMAT='($,'' The Xj('//CHAR1(1:1)//') derivative'//
     '                ' wrt s('//CHAR2(1:1)//') is ['//
     '                CHAR4(1:12)//']: '',E12.5)'
                  ELSE IF(nu.eq.6.or.nu.eq.9.or.nu.EQ.10) THEN
                    IF(nu.EQ.6) THEN
                      CHAR2='1'
                      CHAR3='2'
                    ELSE IF(nu.EQ.9) THEN
                      CHAR2='1'
                      CHAR3='3'
                    ELSE IF(nu.EQ.10) THEN
                      CHAR2='2'
                      CHAR3='3'
                    ENDIF
                    FORMAT='($,'' The Xj('//CHAR1(1:1)//') derivative'//
     '                ' wrt s('//CHAR2(1:1)//') & s('//CHAR3(1:1)//
     '                ') is ['//CHAR4(1:12)//']: '',E12.5)'
                  ENDIF !nu
                  IF(MODE(1:7).EQ.'DEGREES') THEN
                    WRITE(OP_STRING,FMT=FORMAT)
     '                XG(NJ_LOC(NJL_FIBR,2,nr),nk)/
     '                SCAFAC(NJ_LOC(NJL_FIBR,2,nr),nk)*180.d0/PI
                  ELSE
                    WRITE(OP_STRING,FMT=FORMAT)
     '                XG(NJ_LOC(NJL_FIBR,2,nr),nk)/
     '                SCAFAC(NJ_LOC(NJL_FIBR,2,nr),nk)
                  ENDIF
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ENDDO !nk
              ENDIF

              IF(NJ_LOC(NJL_FIBR,0,nr).GE.3) THEN
C ***         Write interpolated sheet angles in IPSHEE file format
                WRITE(CHAR2,'(I4)') nd
                FORMAT='(/$,'' Node number ['//CHAR2(1:4)//']: '',I4)'
                WRITE(OP_STRING,FMT=FORMAT) nd
                CALL WRITES(IOFI2,OP_STRING,ERROR,*9999)
                WRITE(CHAR1(1:1),'(I1)') NJ_LOC(NJL_FIBR,3,nr)
                DO nk=1,NKJ(NJ_LOC(NJL_FIBR,3,nr),NPNODE(1,nr))
                  nu=IDO(nk,1,0,NBJ(NJ_LOC(NJL_FIBR,3,nr),NE_EVFIBR))
                  IF(MODE(1:7).EQ.'DEGREES') THEN
                    WRITE(CHAR4,'(E12.5)') XG(NJ_LOC(NJL_FIBR,3,nr),nk)
     '                *180.d0/PI
                  ELSE
                    WRITE(CHAR4,'(E12.5)') XG(NJ_LOC(NJL_FIBR,3,nr),nk)
                  ENDIF
                  IF(nu.EQ.1) THEN
                    FORMAT='($,'' The Xj('//CHAR1(1:1)//') coordinate'//
     '                ' is ['//CHAR4(1:12)//']: '',E12.5)'

                  ELSE IF(nu.eq.2.or.nu.eq.4.or.nu.EQ.7) THEN
                    IF(nu.EQ.2) THEN
                      CHAR2='1'
                    ELSE IF(nu.EQ.4) THEN
                      CHAR2='2'
                    ELSE IF(nu.EQ.7) THEN
                      CHAR2='3'
                    ENDIF
                    FORMAT='($,'' The Xj('//CHAR1(1:1)//') derivative'//
     '                ' wrt s('//CHAR2(1:1)//') is ['//
     '                CHAR4(1:12)//']: '',E12.5)'
                  ELSE IF(nu.eq.6.or.nu.eq.9.or.nu.EQ.10) THEN
                    IF(nu.EQ.6) THEN
                      CHAR2='1'
                      CHAR3='2'
                    ELSE IF(nu.EQ.9) THEN
                      CHAR2='1'
                      CHAR3='3'
                    ELSE IF(nu.EQ.10) THEN
                      CHAR2='2'
                      CHAR3='3'
                    ENDIF
                    FORMAT='($,'' The Xj('//CHAR1(1:1)//') derivative'//
     '                ' wrt s('//CHAR2(1:1)//') & s('//CHAR3(1:1)//
     '                ') is ['//CHAR4(1:12)//']: '',E12.5)'
                  ENDIF !nu
                  IF(MODE(1:7).EQ.'DEGREES') THEN
                    WRITE(OP_STRING,FMT=FORMAT)
     '                XG(NJ_LOC(NJL_FIBR,3,nr),nk)
     '                /SCAFAC(NJ_LOC(NJL_FIBR,3,nr),nk)*180.d0/PI
                  ELSE
                    WRITE(OP_STRING,FMT=FORMAT)
     '                XG(NJ_LOC(NJL_FIBR,3,nr),nk)/
     '                SCAFAC(NJ_LOC(NJL_FIBR,3,nr),nk)
                  ENDIF
                  CALL WRITES(IOFI2,OP_STRING,ERROR,*9999)
                ENDDO !nk
              ENDIF

            ELSE !.not.DATAPTS
C ***       Output interpolated fibre (and sheet) angles
              FORMAT='(/'' ne='',I2,'' Xi(i):  '',3(E11.4,2X))'
              WRITE(OP_STRING,FMT=FORMAT) NE_EVFIBR,(XIPOS(ni),
     '          ni=1,NIT(NBJ(1,NE_EVFIBR)))
              CALL WRITES(IOFI2,OP_STRING,ERROR,*9999)

              FORMAT='( ''      '', '' X(j,0): '','
     '          //'E11.4,2(2X,E11.4,'' ('',E11.4,''deg)''))'
              WRITE(OP_STRING,FMT=FORMAT) XG(1,1),
     '          (XG(nj,1),XG(nj,1)*180.d0/PI,nj=2,NJT)
              CALL WRITES(IOFI2,OP_STRING,ERROR,*9999)

              XX=FOCUS*DCOSH(XG(1,1))*DCOS(XG(2,1))
              RR=FOCUS*DSINH(XG(1,1))*DSIN(XG(2,1))
              FORMAT='(7X,''( x= '',E11.4,'' r='',E11.4,'')'')'
              WRITE(OP_STRING,FMT=FORMAT) XX,RR
              CALL WRITES(IOFI2,OP_STRING,ERROR,*9999)

              FORMAT='(''       '',A5,'' angle & derivs: '','
     '          //'E11.4,'' ('',E11.4,''deg) '',3(E11.4,X))'
              DO njj=1,NJ_LOC(NJL_FIBR,0,nr)
                nj=NJ_LOC(NJL_FIBR,njj,nr)
                IF(NBJ(nj,NE_EVFIBR).ne.0) THEN
                  WRITE(OP_STRING,FMT=FORMAT) FIELD(njj),
     '              XG(nj,1),XG(nj,1)*180.d0/PI,
     '              ((XG(nj,nk)/SCAFAC(nj,nk)),nk=2,
     '              NKT(0,NBJ(nj,NE_EVFIBR)))
                  CALL WRITES(IOFI2,OP_STRING,ERROR,*9999)
                ENDIF
              ENDDO !njj

              WRITE(OP_STRING,FMT='('' Sheet angles within Xi radius '
     '          //'of '',F5.3)') Radius_xi
              CALL WRITES(IOFI2,OP_STRING,ERROR,*9999)
              Sheet_angle_max=-10.d0
              Sheet_angle_min= 10.d0
              DO i=1,20
                Xi(1)=XIPOS(1)+Radius_xi*DCOS((i-1)*PI/10.d0)
                Xi(2)=XIPOS(2)+Radius_xi*DSIN((i-1)*PI/10.d0)
                Xi(3)=XIPOS(3)
                nj=NJ_LOC(NJL_FIBR,3,nr) !nj index for sheet angle
                nb=NBJ(nj,ne) !basis#  for sheet angle
                GAMA=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '            nb,1,Xi,XE(1,nj))
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                  WRITE(OP_STRING,FMT='('' Xi: '',3E11.4,'
     '              //''' Sheet angle: '',E11.4)') Xi(1),Xi(2),Xi(3),
     '              GAMA
                  CALL WRITES(IOFI2,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
                ENDIF !dop
                IF(GAMA.GT.Sheet_angle_max) Sheet_angle_max=GAMA
                IF(GAMA.LT.Sheet_angle_min) Sheet_angle_min=GAMA
              ENDDO !i
              WRITE(OP_STRING,FMT='('' Sheet angle range: '',2E13.4)')
     '          Sheet_angle_min,Sheet_angle_max
              CALL WRITES(IOFI2,OP_STRING,ERROR,*9999)
            ENDIF !datapts
 100      ENDDO !nd

        ENDIF
        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IF(NJ_LOC(NJL_FIBR,0,nr).GE.3) CALL CLOSEF(IOFI2,ERROR,*9999)
          IOFI=IOOP
          IOFI2=IOOP
        ENDIF !opfile

      ENDIF

      CALL EXITS('EVFIBR')
      RETURN
 9999 CALL ERRORS('EVFIBR',ERROR)
      CALL EXITS('EVFIBR')
      RETURN 1
      END


