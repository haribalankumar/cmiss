      SUBROUTINE LIDATA(IBT,IDO,INP,LD,NBH,NBJ,NDDL,NDLT,
     '  NDP,NEELEM,NELIST,NFLIST,NHE,NHP,NKH,NKHE,NKJE,NPF,NPNE,
     '  NPNODE,NRE,NRLIST,NVHE,NVHP,NVJE,NYNE,NYNP,CURVCORRECT,
     '  EDD,SE,SQ,WD,XA,XE,XID,XP,YP,ZA,ZD,ZE,ZP,STRING,ERROR,*)

C#### Subroutine: LIDATA
C###  Description:
C###    LIDATA lists data point arrays.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'data00.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LD(NDM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NDDL(NEM,NDEM),NDLT(NEM),
     '  NDP(NDM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),NFLIST(0:NFM),
     '  NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NRE(NEM),NRLIST(0:NRM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),EDD(NDM),SE(NSM,NBFM,NEM),
     '  SQ(NDM),WD(NJM,NDM),XA(NAM,NJM,NEM),XE(NSM,NJM),XID(NIM,NDM),
     '  XP(NKM,NVM,NJM,NPM),
     '  YP(NYM,NIYM,NXM),ZA(NAM,NHM,NCM,NEM),ZD(NJM,NDM),ZE(NSM,NHM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,N3CO,NDLIST(2),nr,NTIL,NTRL,nx
      REAL*8 EDDMIN,MU(2),RADIUS,RFROMC,THETA(2),TIME
      CHARACTER FILE*100,RANGE*9
      LOGICAL ALL_REGIONS,BETWEEN,CBBREV,DATA,Dot_product,ELEMENTS,ERR,
     '  FULL,GROUPS,OPFILE,OPSTAT,STATIC,UNDEFORMED

      CALL ENTERS('LIDATA',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list data<;FILENAME>
C###  Description:
C###    Lists data in specified elements and specified range of Xi(3)
C###    to the screen or to the file FILENAME.opdata if the qualifier
C###    is present.  The
C###    'greater ERROR' parameter results in only data parameters with
C###    error greater than ERROR being listed.  The 'full' parameter
C###    results in a more complete listing.
C###  Parameter:    <element (#s/all)[all]>
C###    Specify the list of elements in which to list the data. The
C###    'all' option will list the data for all currently defined
C###    elements in the specified regions.
C###  Parameter:    <from XI_3#[0.0]{>=0.0,<=1.0}>
C###    Specify the start of the range of Xi(3) of the data points to
C###    list. Data points with Xi(3) values outside the range specified
C###    by this and the 'to XI_3#' option will not be listed.
C###  Parameter:    <to XI_3#[1.0]{>=0.0,<=1.0}>
C###    Specify the end of the range of Xi(3) of the data points to
C###    list. Data points with Xi(3) values outside the range specified
C###    by this and the 'from XI_3#' option will not be listed.
C###  Parameter:    <when (TIME#/all)[all]>
C###    Specify the time at which to list the data point values. The
C###    'all' option will list the data point values at all times.
C###  Parameter:    <greater ERROR#[0.0]>
C###    Specify that only data points with an error greater than the
C###    error specified are to be listed. This option is good for
C###    identifying data points that have large errors in their
C###    projections (i.e. have been projected into the wrong element).
C###  Parameter:    <full>
C###    Specify a more complete listing of the data point values.
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the region numbers of the elements in which to list
C###    the data points. The 'all' parameter specifies all currently
C###    defined regions.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<from XI_3#[0.0]{>=0.0,<=1.0}>'
        OP_STRING(4)=BLANK(1:15)//'<to XI_3#[1.0]{>=0.0,<=1.0}'
        OP_STRING(4)=BLANK(1:15)//'<when (TIME#/all)[all]>'
        OP_STRING(5)=BLANK(1:15)//'<greater ERROR#[0.0]>'
        OP_STRING(6)=BLANK(1:15)//'<full>'
        OP_STRING(7)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)


C---------------------------------------------------------------------

C#### Command: FEM list data statistics
C###  Parameter:    <region (#s/all)[1]>
C###  Description:
C###    List general statistics about a set of data.
C###    eg. the centroid, max/min positions.
C###  Parameter: <between DATAPOINT1#,DATAPOINT2#>
C###    Reduces the output to the subset specified by
C###    DATAPOINT1 and DATAPOINT2. Also calculates
C###    the distance between the two points.

        OP_STRING(1)=STRING(1:IEND)//' statistics'
        OP_STRING(2)=BLANK(1:15)//'<between DATAPOINT1#,DATAPOINT2#>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)


C---------------------------------------------------------------------

C#### Command: FEM list data groups
C###  Parameter:    <region (#s/all)[1]>
C###  Description:
C###    List data groups.

        OP_STRING(1)=STRING(1:IEND)//' groups'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)


C---------------------------------------------------------------------

C#### Command: FEM list data errors
C###  Description:
C###    List the data point errors and statistics. Gives the average, RMS etc. errors
C###    for whole mesh and
C###    (if use 'by_elements') for each specified element.
C###  Parameter:    <element (#s/all)[all]>
C###    Specify the list of elements in which to list the data. The
C###    'all' option will list the data for all currently defined
C###    elements in the specified regions.
C###  Parameter:    <greater ERROR#[0.0]>
C###    Specify that only data points with an error greater than the
C###    error specified are to be listed. This option is good for
C###    identifying data points that have large errors in their
C###    projections (i.e. have been projected into the wrong element).
C###  Parameter:    <by_elements>
C###    Specify that the data point errors and statistics are to be
C###    listed for each element. If this parameter is omitted then
C###    the data point errors and statistics are listed for the
C###    entire set of elements as specified by the 'elements' parameter.
C###  Parameter:    <full>
C###    Specify a more complete listing of the data point values.
C###  Parameter:    <dot>
C###    Specify that the data point errors are to be adjusted by the
C###    data point weights i.e. that the data point error is given by
C###    the dot product of the data error vector and the data weight
C###    vector.
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the region numbers of the elements in which to list
C###    the data points. The 'all' parameter specifies all currently
C###    defined regions.
C###  Parameter:    <(undeformed/deformed) [undeformed]>
C###    Specify the geometery to project onto.

        OP_STRING(1)=STRING(1:IEND)//' errors'
        OP_STRING(2)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<greater ERROR#[0.0]>'
        OP_STRING(4)=BLANK(1:15)//'<by_elements>'
        OP_STRING(5)=BLANK(1:15)//'<full>'
        OP_STRING(6)=BLANK(1:15)//'<dot>'
        OP_STRING(7)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(8)=BLANK(1:15)//'<(undeformed/deformed) [undeformed]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list data range
C###  Description:
C###    List data points within specified range of mu and theta.
C###  Parameter:    <mu=MU1#[0.0]{degrees},MU2#[0.0]{degrees}>
C###    Specify the start and finish values of the mu range in degrees.
C###  Parameter:    <theta=THETA1#,[0.0]{degrees}THETA2#[0.0]{degrees}>
C###    Specify the start and finish values of the theta range
C###    in degrees.

        OP_STRING(1)=STRING(1:IEND)//' range'
        OP_STRING(2)=BLANK(1:15)//
     '    '<mu=MU1#[0.0]{degrees},MU2#[0.0]{degrees}'
        OP_STRING(3)=BLANK(1:15)//
     '    '<theta=THETA1#[0.0]{degrees},THETA2#[0.0]{degrees}'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list data circle
C###  Description:
C###    List data points within Xi-coord circle around mu & theta.
C###  Parameter:    <mu=MU#[0.0]{degrees}>
C###    Specify the mu coordinate of the centre of the circle in
C###    degrees.
C###  Parameter:    <theta=THETA#[0.0]{degrees}>
C###    Specify the theta coordinate of the centre of the circle in
C###    degrees.
C###  Parameter:    <radius=RADIUS_in_Mu_Theta#[0.0]>
C###    Specify the radius of the circle.

        OP_STRING(1)=STRING(1:IEND)//' circle'
        OP_STRING(2)=BLANK(1:15)//'<mu=MU#[0.0]{degrees}>'
        OP_STRING(3)=BLANK(1:15)//'<theta=THETA#[0.0]{degrees}>'
        OP_STRING(4)=BLANK(1:15)//'<radius=RADIUS_in_Mu_Theta#[0.0]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIDATA',ERROR,*9999)
      ELSE

C LKC 10-JAN-2000 new assert
        CALL ASSERT(USE_DATA.EQ.1,'Set USE_DATA to 1',ERROR,*9999)

        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opdata','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL ASSERT(NRLIST(0).EQ.1,'>> Only implemented for 1 region'
     '    ,ERROR,*9999)
        nr=NRLIST(1)
        nx=1 !temporary

        GROUPS=.FALSE.
        IF(CBBREV(CO,'ERRORS',2,noco+1,NTCO,N3CO)) THEN
          ERR  =.TRUE.
          DATA =.FALSE.
          OPSTAT=.FALSE.
          BETWEEN=.FALSE.
          RANGE=' '
C CPB 12/1/92 adding list data error elements
          IF(CBBREV(CO,'BY_ELEMENTS',2,noco+1,NTCO,N3CO)) THEN
            ELEMENTS=.TRUE.
          ELSE
            ELEMENTS=.FALSE.
          ENDIF
          IF(CBBREV(CO,'DOT',1,noco+1,NTCO,N3CO)) THEN
            Dot_product=.TRUE.
          ELSE
            Dot_product=.FALSE.
          ENDIF
          IF(CBBREV(CO,'DEFORMED',1,noco+1,NTCO,N3CO)) THEN
            UNDEFORMED=.FALSE.
          ELSE
            UNDEFORMED=.TRUE.
          ENDIF
        ELSE IF(CBBREV(CO,'RANGE',2,noco+1,NTCO,N3CO)) THEN
          ERR  =.FALSE.
          DATA =.FALSE.
          OPSTAT=.FALSE.
          BETWEEN=.FALSE.
          RANGE='RECTANGLE'
          IF(CBBREV(CO,'MU',1,noco+1,NTCO,N3CO)) THEN
            CALL PARSRL(CO(N3CO+1),2,NTRL,MU,ERROR,*9999)
            MU(1)=MU(1)*PI/180.0d0
            MU(2)=MU(2)*PI/180.0d0
          ELSE
            MU(1)=0.d0
            MU(2)=0.d0
          ENDIF
          IF(CBBREV(CO,'THETA',1,noco+1,NTCO,N3CO)) THEN
            CALL PARSRL(CO(N3CO+1),2,NTRL,THETA,ERROR,*9999)
            THETA(1)=THETA(1)*PI/180.0d0
            THETA(2)=THETA(2)*PI/180.0d0
          ELSE
            THETA(1)=0.d0
            THETA(2)=0.d0
          ENDIF

        ELSE IF(CBBREV(CO,'CIRCLE',2,noco+1,NTCO,N3CO)) THEN
          ERR  =.FALSE.
          DATA =.FALSE.
          OPSTAT=.FALSE.
          BETWEEN=.FALSE.
          RANGE='CIRCLE'
          IF(CBBREV(CO,'MU',1,noco+1,NTCO,N3CO)) THEN
            CALL PARSRL(CO(N3CO+1),1,NTRL,MU,ERROR,*9999)
            MU(1)=MU(1)*PI/180.0d0
          ELSE
            MU(1)=0.d0
          ENDIF
          IF(CBBREV(CO,'THETA',1,noco+1,NTCO,N3CO)) THEN
            CALL PARSRL(CO(N3CO+1),2,NTRL,THETA,ERROR,*9999)
            THETA(1)=THETA(1)*PI/180.0d0
          ELSE
            THETA(1)=0.d0
          ENDIF
          IF(CBBREV(CO,'RADIUS',1,noco+1,NTCO,N3CO)) THEN
            RADIUS=RFROMC(CO(N3CO+1))*PI/180.0d0
          ELSE
            RADIUS=0.d0
          ENDIF

        ELSE IF(CBBREV(CO,'STATISTICS',4,noco+1,NTCO,N3CO)) THEN
          OPSTAT=.TRUE.
          DATA =.FALSE.
          ERR  =.FALSE.
          IF(CBBREV(CO,'BETWEEN',3,noco+1,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),2,NTIL,NDLIST,ERROR,*9999)
            BETWEEN=.TRUE.
          ELSE
            BETWEEN=.FALSE.
          ENDIF
        ELSE IF(CBBREV(CO,'GROUPS',4,noco+1,NTCO,N3CO)) THEN
          GROUPS=.TRUE.
        ELSE ! defaults
          ERR  =.FALSE.
          DATA =.TRUE.
          OPSTAT=.FALSE.
          BETWEEN=.FALSE.
          RANGE=' '
        ENDIF

        CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '    ERROR,*9999)

        IF(CBBREV(CO,'FROM',2,noco+1,NTCO,N3CO)) THEN
          XI3MIN=RFROMC(CO(N3CO+1))
        ELSE
          XI3MIN=0.0d0
        ENDIF
        IF(CBBREV(CO,'TO',2,noco+1,NTCO,N3CO)) THEN
          XI3MAX=RFROMC(CO(N3CO+1))
        ELSE
          XI3MAX=1.0D0
        ENDIF

        IF(CBBREV(CO,'GREATER',1,noco+1,NTCO,N3CO)) THEN
          CALL PARSRE(CO(N3CO+1),EDDMIN,ERROR,*9999)
        ELSE
          EDDMIN=0.0d0
        ENDIF

        IF(CBBREV(CO,'FULL',2,noco+1,NTCO,N3CO)) THEN
          FULL=.TRUE.
        ELSE
          FULL=.FALSE.
        ENDIF

        IF(CBBREV(CO,'WHEN',2,noco+1,NTCO,N3CO)) THEN
          TIME=RFROMC(CO(N3CO+1))
          STATIC=.FALSE.
        ELSE
          STATIC=.TRUE.
          TIME=0.0d0
        ENDIF

        CALL OPDATA(GROUPS,IBT,IDO,INP,LD,NBH,NBJ,NDDL,
     '    NDLIST,NDLT,NDP,
     '    NEELEM,NELIST,NFLIST,NHE(1,nx),NHP(1,nr,nx),
     '    NKH(1,1,1,nr),NKHE,NKJE,
     '    NPF,NPNE,NPNODE,nr,NRE,NVHE,
     '    NVHP(1,1,1,nr),NVJE,nx,NYNE,NYNP,CURVCORRECT,
     '    EDD,EDDMIN,MU,RADIUS,SE,SQ,THETA,TIME,WD,XA,XE,XID,XP,YP,
     '    ZA,ZD,ZE,ZP,BETWEEN,DATA,Dot_product,ELEMENTS,ERR,FULL,OPSTAT,
     '    UNDEFORMED,
     '    RANGE,STATIC,ERROR,*9999)

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LIDATA')
      RETURN
 9999 CALL ERRORS('LIDATA',ERROR)
      CALL EXITS('LIDATA')
      RETURN 1
      END


