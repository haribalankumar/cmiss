      SUBROUTINE LIGRID(NQLIST,NRLIST,NAQ,NEELEM,NENQ,NLATNE,NLATNQ,
     '  NLATPNQ,NLQ,NQET,NQGP,NQGP_PIVOT,NQNE,NQNLAT,NQNP,NQS,NQSCNB,
     '  NQXI,NWQ,NXQ,NXLIST,AQ,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,GCHQ,GUQ,PROPQ,
     '  XQ,YQ,YQS,STRING,ERROR,*)

C#### Subroutine: LIGRID
C###  Description:
C###    LIGRID lists finite difference grid point information

      IMPLICIT NONE
      INCLUDE 'call00.cmn'
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'nqloc00.inc'

!     Parameter List
      INTEGER NQLIST(0:NQM),NRLIST(0:NRM),NXLIST(0:NXM),
     &  NAQ(NQM,NAM),NEELEM(0:NE_R_M,0:NRM),NENQ(0:8,NQM),
     &  NLATNE(NEQM+1),NLATNQ(NEQM*NQEM),NLATPNQ(NQM),NLQ(NQM),
     &  NQET(NQSCM),NQGP(0:NQGM,NQM),NQGP_PIVOT(NQGM,NQM),
     &  NQNE(NEQM,NQEM),NQNLAT(NEQM*NQEM),NQNP(NPM),NQS(NEQM),
     &  NQSCNB(NQSCM),NQXI(0:NIM,NQSCM),NWQ(8,0:NQM,NAM),
     &  NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      REAL*8 AQ(NMAQM,NQM),CQ(NMM,NQM,NXM),DNUDXQ(3,3,NQM),
     &  DXDXIQ(3,3,NQM),DXDXIQ2(3,3,NQM),GCHQ(3,NQM),GUQ(3,3,NQM),
     &  PROPQ(3,3,4,2,NQM,NXM),XQ(NJM,NQM),YQ(NYQM,NIQM,NAM,NXM),
     &  YQS(NIQSM,NQM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER ARRAY_INDEX,IBEG,IEND,IFROMC,N3CO,na,ne,NELIST(0:NEM),
     &  no_nelist,no_nqlist,no_nrlist,nr,nx,nxc
      REAL*8 RADIUS(2),RFROMC
      CHARACTER FILE*100,TYPE*13
      LOGICAL ALL_REGIONS,ARC_LENGTH,CBBREV,FLUX,FULL,GRID_LATTICE,
     &  HISTORY,NIDS,OPFILE,POTENTIAL,RMS
!     Functions
      LOGICAL ABBREV

      CALL ENTERS('LIGRID',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list grid<;FILENAME>
C###  Description:
C###    List grid point information. Grid point properties are
C###    written to the file FILENAME (with extension .opgrid)
C###  Parameter:    <grid (#s/all)[all]>
C###    Specify which grid points are to be listed. A list of
C###    grid points can be listed or a group name can be given.
C###  Parameter:    <elements (#)[1]>
C###    Specify which element's grid points are to be listed.
C###  Parameter:    <full>
C###    List full information about grid points.
C###  Parameter:    <history>
C###    This parameter is currently unsupported.
C###  Parameter:    <lattice>
C###    List the lattice grid scheme information.       
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<grid (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<element (#)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<full>'
        OP_STRING(5)=BLANK(1:15)//'<history>'
        OP_STRING(6)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(7)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list grid<;FILENAME> geometry/material/yq/yqs
C###  Description:
C###    Lists grid point material information or solution/grid variable
C###    array information. Grid point information is
C###    written to the file FILENAME (with extension .opgrid)
C###  Parameter:    <index (#/all)[all]>
C###    Specify which array index is to be listed.  Specifying `v'
C###    indicates the potential.  If no index is specified then all of
C###    them are listed.
C###  Parameter:    <grid (#s/all)[all]>
C###    Specify which grid points are to be listed. A list of
C###    grid points can be listed or a group name can be given.
C###  Parameter:    <elements (#)[1]>
C###    Specify which element's grid points are to be listed.
C###  Parameter:    <level #[1]>
C###    For multigrid applications specify the grid level at which
C###    the information is specified.
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> material/yq/yqs'
        OP_STRING(2)=BLANK(1:15)//'<index (#/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<grid (#s/all)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<element (#)[1]>'
        OP_STRING(5)=BLANK(1:15)//'<level #[1]>'
        OP_STRING(6)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(7)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list grid<;FILENAME> connectivity
C###  Description:
C###    This command outputs the NXQ array of grid connectivity
C###    for the first multigrid level only. Grid point information
C###    is written to the file FILENAME (with extension .opgrid)
C###  Parameter:    <grid (#s/all)[all]>
C###    Specify which grid points are to be listed. A list of
C###    grid points can be listed or a group name can be given.
C###  Parameter:    <elements (#)[1]>
C###    Specify which element's grid points are to be listed.
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> connectivity'
        OP_STRING(2)=BLANK(1:15)//'<grid (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<element (#)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(5)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list grid<;FILENAME> groups
C###  Description:
C###    This command lists all grid groups. Group information
C###    is written to the file FILENAME (with extension .opgrid)

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> groups'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list grid<;FILENAME> statistics
C###  Description:
C###    This command calculates and outputs grid spacing
C###    information. The information is written to the file
C###    FILENAME (with extension .opgrid)
C###  Parameter:    <grid (#s/all)[all]>
C###    Specify which grid points are to be listed. A list of
C###    grid points can be listed or a group name can be given.
C###  Parameter:    <elements (#)[1]>
C###    Specify which element's grid points are to be listed.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> statistics'
        OP_STRING(2)=BLANK(1:15)//'<grid (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<element (#)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list grid<;FILENAME> adaptive_levels
C###  Description:
C###    This command lists grid point levels for adaptive
C###    residual calculations.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> adaptive_levels'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list grid<;FILENAME> nodes
C###  Description:
C###    This command checks that a coupled FD/BEM problem has
C###    been defined then lists grid points numbers coupled to node
C###    numbers. The information is written to the file
C###    FILENAME (with extension .opgrid)

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> nodes'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list grid<;FILENAME> flux
C###  Description:
C###    This command calculates the flux at all boundary grid points
C###    for the specified regions and class. The information is
C###    written to the file FILENAME (with extension .opgrid)
C###  Parameter:    <grid (#s/all)[all]>
C###    Specify which grid points are to be listed. A list of
C###    grid points can be listed or a group name can be given.
C###  Parameter:    <level #[1]>
C###    For multigrid applications specify the grid level at which
C###    the information is specified.
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> flux'
        OP_STRING(2)=BLANK(1:15)//'<grid (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<level #[1]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(5)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list grid<;FILENAME> errors
C###  Description:
C###    This command is used for solve8 problems in which Laplace's
C###    equation is solved and checked against a known analytic
C###    solution. The command gives a summary of the differences
C###    between the analytic and computed solutions. The information is
C###    written to the file FILENAME (with extension .opgrid)
C###  Parameter:    <grid (#s/all)[all]>
C###    Specify which grid points are to be listed. A list of
C###    grid points can be listed or a group name can be given.
C###  Parameter:    <level #[1]>
C###    For multigrid applications specify the grid level at which
C###    the information is specified.
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Parameter:    <RMS>
C###    A measure of the absolute error.
C###  Parameter:    <NIDS>
C###    A measure of the relative error.
C###  Parameter:    <potential>
C###    Give errors in potentials.
C###  Parameter:    <flux>
C###    Give errors in boundary fluxes.
C###  Parameter:    <arc_length_derivative>
C###    Give errors in arc length derivative values, note this is
C###    only valid for circular domains.
C###  Parameter:    <inner_radius #[1]>
C###    For use with arc length derivative terms, gives radius of
C###    the inside of the domain.
C###  Parameter:    <outer_radius #[2]>
C###    For use with arc length derivative terms, gives radius of
C###    the inside of the domain.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> errors'
        OP_STRING(2)=BLANK(1:15)//'<grid (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<level #[1]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(5)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(6)=BLANK(1:15)//'<RMS>'
        OP_STRING(7)=BLANK(1:15)//'<NIDS>'
        OP_STRING(8)=BLANK(1:15)//'<potential>'
        OP_STRING(9)=BLANK(1:15)//'<flux>'
        OP_STRING(10)=BLANK(1:15)//'<arc_length_derivative>'
        OP_STRING(11)=BLANK(1:15)//'<inner_radius #[1]>'
        OP_STRING(12)=BLANK(1:15)//'<outer_radius #[2]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list grid<;FILENAME> activation_time
C###  Description:
C###    This command is used for listing the activation times
C###    for activation problems based on the steepest Vm gradient.
C###    The information is
C###    written to the file FILENAME (with extension .opgrid)
C###  Parameter:    <grid (#s/all)[all]>
C###    Specify which grid points are to be listed. A list of
C###    grid points can be listed or a group name can be given.
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> activation_time'
        OP_STRING(2)=BLANK(1:15)//'<grid (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list grid<;FILENAME> auxiliary
C###  Description:
C###    This command is used for listing the auxiliary parameters defined
C###    for all grid points
C###  Parameter:    <grid (#s/all)[all]>
C###    Specify which grid points are to be listed. A list of
C###    grid points can be listed or a group name can be given.
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> auxiliary'
        OP_STRING(2)=BLANK(1:15)//'<grid (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list grid<;FILENAME> lattice
C###  Description:
C###    This command is used for listing the lattice grid
C###    points and their mapping to the global grid scheme.
C###    The information is written to the file
C###    FILENAME (with extension .opgrid)
C###  Parameter:    <by_grid>
C###    List the global grid points and the lattice grid points
C###    that they are mapped to.
C###  Parameter:    <grid (#s/all)[all]>
C###    Specify which grid points are to be listed. A list of
C###    grid points can be listed or a group name can be given.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> lattice'
        OP_STRING(2)=BLANK(1:15)//'<by_grid>'
        OP_STRING(3)=BLANK(1:15)//'<grid (#s/all)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

        
      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIGRID',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opgrid','NEW',
     '      'SEQUEN','FORMATTED',255,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        IF(USE_LAT.EQ.0) THEN
          CALL ASSERT(NQT.GT.0,'>>No grid points defined',ERROR,*9999)
        ELSE
          CALL ASSERT(CALL_GRID,' Must define grid first',
     '      ERROR,*9999)
        ENDIF
        IF(CBBREV(CO,'GEOMETRY',2,noco+1,NTCO,N3CO)) THEN
          TYPE='GEOMETRY'
        ELSE IF(CBBREV(CO,'MATERIAL',1,noco+1,NTCO,N3CO)) THEN
          TYPE='MATERIAL'
          CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
     '      ERROR,*9999)
        ELSE IF(CBBREV(CO,'CONNECTIVITY',3,noco+1,NTCO,N3CO)) THEN
          TYPE='CONNECTIVITY'
        ELSE IF(CBBREV(CO,'STATISTICS',3,noco+1,NTCO,N3CO)) THEN
          TYPE='STATISTICS'
        ELSE IF(CBBREV(CO,'ADAPTIVE_LEVELS',2,noco+1,NTCO,N3CO)) THEN
          TYPE='ADAPTIVE'
        ELSE IF(CBBREV(CO,'ACTIVATION_TIMES',2,noco+1,NTCO,N3CO)) THEN
          TYPE='ACTIVTIMES'
        ELSE IF(CBBREV(CO,'AUXILIARY',2,noco+1,NTCO,N3CO)) THEN
          TYPE='AUXILIARY'
        ELSE IF(CBBREV(CO,'NODES',1,noco+1,NTCO,N3CO)) THEN
          TYPE='NODES'
          CALL ASSERT(.NOT.UP_NQNP,' >>Update node grid first',
     '      ERROR,*9999)
        ELSE IF(CBBREV(CO,'YQS',3,noco+1,NTCO,N3CO)) THEN
          TYPE='YQS'
        ELSE IF(CBBREV(CO,'YQ',2,noco+1,NTCO,N3CO)) THEN
          TYPE='YQ'
          CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
     '      ERROR,*9999)
        ELSE IF(CBBREV(CO,'GROUPS',1,noco+1,NTCO,N3CO)) THEN
          TYPE='GROUPS'
C        ELSE IF(CBBREV(CO,'FLUX',2,noco+1,NTCO,N3CO)) THEN
C          TYPE='FLUX'
C          CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
C     '      ERROR,*9999)
        ELSE IF(CBBREV(CO,'ERRORS',3,noco+1,NTCO,N3CO)) THEN
          TYPE='ERRORS'
          IF(CBBREV(CO,'RMS',3,noco+1,NTCO,N3CO)) THEN
            RMS=.TRUE.
          ELSE
            RMS=.FALSE.
          ENDIF
          IF(CBBREV(CO,'NIDS',4,noco+1,NTCO,N3CO)) THEN
            NIDS=.TRUE.
          ELSE
            NIDS=.FALSE.
          ENDIF
          IF(CBBREV(CO,'POTENTIAL',3,noco+1,NTCO,N3CO)) THEN
            POTENTIAL=.TRUE.
          ELSE
            POTENTIAL=.FALSE.
          ENDIF
          IF(CBBREV(CO,'FLUX',4,noco+1,NTCO,N3CO)) THEN
            FLUX=.TRUE.
          ELSE
            FLUX=.FALSE.
          ENDIF
          IF(CBBREV(CO,'ARC_LENGTH_DERIVATIVE',3,noco+1,NTCO,N3CO)) THEN
            ARC_LENGTH=.TRUE.
            IF(CBBREV(CO,'INNER_RADIUS',3,noco+1,NTCO,N3CO)) THEN
              RADIUS(1)=RFROMC(CO(N3CO+1))
            ELSE
              RADIUS(1)=1.0d0
            ENDIF
            IF(CBBREV(CO,'OUTER_RADIUS',3,noco+1,NTCO,N3CO)) THEN
              RADIUS(2)=RFROMC(CO(N3CO+1))
            ELSE
              RADIUS(2)=1.0d0
            ENDIF
          ELSE
            ARC_LENGTH=.FALSE.
          ENDIF
        ELSE IF(CBBREV(CO,'FLUX',2,noco+1,NTCO,N3CO)) THEN
          TYPE='FLUX'
          CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
     '      ERROR,*9999)
        ELSE IF(CBBREV(CO,'LATTICE',2,noco+1,NTCO,N3CO)) THEN
          CALL ASSERT(USE_LAT.EQ.1,
     '      '>>Choose lattice grid points in ipgrid',ERROR,*9999)
          TYPE='LATTICE'
        ELSE
          TYPE='CURRENT'
C PJH 1Jun99 'list grid' should work without eqtn defined
C          CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
C     '      ERROR,*9999)
          IF(CBBREV(CO,'FULL',2,noco+1,NTCO,N3CO)) THEN
            FULL=.TRUE.
          ELSE
            FULL=.FALSE.
          ENDIF
          IF(CBBREV(CO,'HISTORY',1,noco+1,NTCO,N3CO)) THEN
            HISTORY=.TRUE.
          ELSE
            HISTORY=.FALSE.
          ENDIF
        ENDIF !type

        IF(TYPE(1:8).EQ.'MATERIAL'.OR.
     '    TYPE(1:3).EQ.'YQS'.OR.
     '    TYPE(1:2).EQ.'YQ') THEN
          IF(CBBREV(CO,'INDEX',3,noco+1,NTCO,N3CO)) THEN
            IF(TYPE.EQ.'YQ'.AND.ABBREV(CO(N3CO+1),'V',1)) THEN !potential
              CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,ARRAY_INDEX,NIQ_V,
     '          ERROR,*9999)
            ELSE
              ARRAY_INDEX=IFROMC(CO(N3CO+1))
            ENDIF
          ELSE
            ARRAY_INDEX=0
          ENDIF
        ENDIF

        IF(nx.EQ.0) nx=1

        IF(TYPE(1:6).NE.'GROUPS') THEN
          IF(CBBREV(CO,'LATTICE',3,noco+1,NTCO,N3CO)) THEN
            IF(CBBREV(CO,'BY_GRID',2,noco+1,NTCO,N3CO)) THEN
              GRID_LATTICE=.TRUE.
              CALL PARSE_GRID(NQLIST,noco,NTCO,CO,ERROR,*9999)
            ELSE
              GRID_LATTICE=.FALSE.
            ENDIF
          ELSEIF(CBBREV(CO,'GRID',1,noco+1,NTCO,N3CO)) THEN
            CALL PARSE_GRID(NQLIST,noco,NTCO,CO,ERROR,*9999)
          ELSEIF(CBBREV(CO,'ELEMENT',1,noco+1,NTCO,N3CO)) THEN
C *** DPN 30 November 1999 - adding option to list grid points within
C ***   specified elements only
            CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     &        ERROR,*9999)
            DO no_nelist=1,NELIST(0)
              ne=NELIST(no_nelist)
              CALL ASSERT(ne.GT.0,'Invalid element number',ERROR,*9999)
C***          add the grid points from the specified element into NQLIST
              DO no_nqlist=1,NQET(NQS(ne))
                NQLIST(NQLIST(0)+no_nqlist)=NQNE(ne,no_nqlist)
              ENDDO !no_nqlist
              NQLIST(0)=NQLIST(0)+NQET(NQS(ne))
            ENDDO
          ELSEIF(CBBREV(CO,'REGION',2,noco+1,NTCO,N3CO)) THEN
            NQLIST(0)=0
            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              DO no_nqlist=NQR(1,nr),NQR(2,nr)
                NQLIST(0)=NQLIST(0)+1
                NQLIST(NQLIST(0))=no_nqlist
              ENDDO !no_nqlist
            ENDDO !no_nrlist
          ELSE
            NQLIST(0)=NQT
            DO no_nqlist=1,NQLIST(0)
              NQLIST(no_nqlist)=no_nqlist
            ENDDO
          ENDIF
        ENDIF

        IF(TYPE(1:8).EQ.'ADAPTIVE') THEN
          na=1
        ELSE
          IF(CBBREV(CO,'LEVEL',1,noco+1,NTCO,N3CO)) THEN
            na=IFROMC(CO(N3CO+1))
          ELSE
            na=1
          ENDIF
        ENDIF !levels

        IF(TYPE(1:6).EQ.'GROUPS') THEN
          CALL OPGRIDG(ERROR,*9999)
        ELSE
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            CALL OPGRID(TYPE,FULL,GRID_LATTICE,HISTORY,ARRAY_INDEX,
     '        na,NQLIST,nr,NAQ,NEELEM,NENQ,NLATNE,NLATNQ,NLATPNQ,
     '        NLQ,NQGP,NQGP_PIVOT,NQNLAT,NQNP,NQS,NQSCNB,NQXI,NWQ,nx,
     '        nxc,NXQ,AQ,CQ(1,1,nx),DNUDXQ,DXDXIQ,DXDXIQ2,GCHQ,GUQ,
     '        PROPQ(1,1,1,1,1,nx),RADIUS,XQ,YQ(1,1,na,nx),YQS,
     '        ARC_LENGTH,FLUX,NIDS,POTENTIAL,RMS,ERROR,*9999)
          ENDDO !no_nrlist
        ENDIF

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF

      ENDIF

      CALL EXITS('LIGRID')
      RETURN
 9999 CALL ERRORS('LIGRID',ERROR)
      CALL EXITS('LIGRID')
      RETURN 1
      END


