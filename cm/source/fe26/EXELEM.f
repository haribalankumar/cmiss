      SUBROUTINE EXELEM(IBT,IDO,INP,NBH,NBJ,NBJF,NEELEM,NELIST,
     '  NFF,NFLIST,NHE,NHP,NKB,NKHE,NKJ,NKJE,NLATNE,NLF,NLL,NLLIST,NNB,
     '  NNF,NPNE,NQET,NQNE,NQNLAT,NQS,NQXI,NRLIST,NSB,NVHE,NVJE,NVJP,NW,
     '  NXLIST,NYNQ,SE,YQ,YQS,CELL_ICQS_VALUE,CELL_ICQS_SPATIAL,
     '  CELL_ICQS_NAMES,CELL_RCQS_VALUE,CELL_RCQS_SPATIAL,
     '  CELL_RCQS_NAMES,CELL_YQS_VALUE,CELL_YQS_SPATIAL ,CELL_YQS_NAMES,
     '  ICQS_SPATIAL,IICQS_SPATIAL,IRCQS_SPATIAL,RCQS_SPATIAL,
     '  STRING,ERROR,*)
      
C#### Subroutine: EXELEM
C###  Description:
C###    EXELEM exports element data from finite element data base.

C****   ELEM_TYPE  is 'geometry'
C****   ELEM_NAME  is '1..#elements' or group name
C****   ELEM_TOTAL is total #elements in exported list
C****   NPNE(nn,nb,ne) are global node numbers
C****   SE(ns,nb,ne) are element scaling factors
C**** NOTE: for socket connection CONNID2 is defined as a parameter
C**** specifying the data transfer socket number
C     SPECIAL_BASIS_FLAG is zero except where the basis is not a tensor
C     product basis
C       1 is Hermite simplex with apex at node 1
C       2 is Hermite simplex with apex at node 3
C       3 is sector element
C       4 other bases with cross derivs set to zero

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cell02.cmn'
      INCLUDE 'cellml.cmn'
      INCLUDE 'cmis00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'fsklib.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     &  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     &  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),NFF(6,NEM),NHE(NEM,NXM),
     &  NFLIST(0:NFM),
     &  NHP(NPM,0:NRM,NXM),NKB(2,2,2,NNM,NBFM),NKHE(NKM,NNM,NHM,NEM),
     &  NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),NLATNE(NEQM+1),NLF(4,NFM),
     &  NLL(12,NEM),NLLIST(0:NLM),NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),
     &  NPNE(NNM,NBFM,NEM),NQET(NQSCM),NQNE(NEQM,NQEM),
     &  NQNLAT(NEQM*NQEM),NQS(NEQM),NQXI(0:NIM,NQSCM),NRLIST(0:NRM),
     &  NSB(NKM,NNM,NBFM),NVHE(NNM,NBFM,NHM,NEM),NVJE(NNM,NBFM,NJM,NEM),
     &  NVJP(NJM,NPM),NW(NEM,3,NXM),NXLIST(0:NXM),
     &  NYNQ(NHM,NQM,0:NRCM,NXM),CELL_ICQS_SPATIAL(NQIM,NQVM),
     &  CELL_ICQS_VALUE(NQIM,NQVM),CELL_RCQS_SPATIAL(NQRM,NQVM),
     &  CELL_YQS_SPATIAL(NIQSM,NQVM),ICQS_SPATIAL(NQISVM,NQM),
     &  IICQS_SPATIAL(0:NQISVM,NQVM),IRCQS_SPATIAL(0:NQRSVM,NQVM)
      REAL*8 SE(NSM,NBFM,NEM),YQ(NYQM,NIQM,NAM,NXM),YQS(NIQSM,NQM),
     '  CELL_RCQS_VALUE(NQRM,NQVM),CELL_YQS_VALUE(NIQSM,NQVM),
     '  RCQS_SPATIAL(NQRSVM,NQM)
      CHARACTER ERROR*(*),
     '  CELL_ICQS_NAMES(NQIM,NQVM)*(CELL_NAME_LENGTH),
     '  CELL_RCQS_NAMES(NQRM,NQVM)*(CELL_NAME_LENGTH),
     '  CELL_YQS_NAMES(NIQSM,NQVM)*(CELL_NAME_LENGTH),
     '  STRING*(MXCH)
!     Local Variables
      INTEGER MAX_ELEMENT_NODES,SIMPLEX_BASIS_DIRNS(0:3)
      PARAMETER (MAX_ELEMENT_NODES=64)
      INTEGER CLEN,COLLAPSEDNODE,COUNT,DATA_TYPE,ELEMENT_DIMENSION,
     '  ELEMENT_NODE_LIST(0:MAX_ELEMENT_NODES),
     '  ELEM_NODES_PREV,FIELD_BASE_TYPE,
     '  FIELD_NODE_LIST(0:MAX_ELEMENT_NODES),
     '  I,IBEG,IBEG1,IBEG2,IBEG3,IBEG4,IDOI(3),IDOTEMP(8,3),IEND,IEND1,
     '  IEND2,IEND3,IEND4,IFROMC,INTSTR(1024),ISOCKET(10),iy,j,k,MK,
     '  N3CO,na,NAE,NAF,nb,NBLIST(0:99),nc,ne,NE_FIELD,
     '  NE_SHAPE,nf,nf_e,NFIELDT,ni,ni1,ni2,nicollapse,nj,
     '  njh,NJH_COUNT,NJH_FIELD_BASE_LIST(0:500),NJH_LIST(99),NJHT,
     '  njj,njj1,njj2,njTOT_NE,njTOT_NE_FIELD,
     '  nk,nkk,nkk2,NKLIST(8),NKLIST1(8),NK_TOT,nl,nlat,NLFF(4),nn,nnn,
     '  NN_TOT,NO_NELIST,no_njh,np,NP_ELEMENT,nqsc,nqsc_old,nqsv,
     '  nr,nrr,nqesc,NS_TOT,nss,NST_LIST(NBFM),NUMNODES(3),num_nn,
     '  NUMCOLLAPSE,nv_offset,nx,nxc,
     '  offset_elem,offset_node,offset_line,
     '  offset_face,POSITION(4),SFLIST(64),SPECIAL_BASIS_FLAG
      REAL*8 SE_LIST(NKM)
      CHARACTER BASES*54,CHAR1*5,CHAR2*5,CHAR3*5,CHAR4*20,
     '  ELEM_NAME*50,ELEM_TYPE*50,FIELD_EX_TYPE*8,FIELD_NAME*50,
     '  FILE*200,OUTPUT*11,SHAPE_STRING*30,SIMPLEX_BASIS_LIST*6,
     '  GRIDVARIABLE*3
      LOGICAL ABBREV,ALL_REGIONS,ATCOLLAPSE,AUTONAME,
     '  CBBREV,CELL_FIELD,COLLAPSED,DATAFILE,
     '  FIELDS_CHANGED,FIRST_ELEMENT,FOUND,GRID_NUMBERS,LARGE_FORMAT,
     '  NJLIST(3),SET_FIELD_NAME,SHAPE_CHANGED,EXPORT_CELL,CELL_TEST
      PARAMETER(DATA_TYPE = 9)  !Element data

      DATA IDOI/1,1,1/

      CALL ENTERS('EXELEM',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        WRITE(CHAR2,'(I5)') NEELEM(1,1)
        WRITE(CHAR3,'(I5)') NEELEM(NEELEM(0,1),1)
        CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
        CALL STRING_TRIM(CHAR3,IBEG3,IEND3)
        CHAR4=CHAR2(IBEG2:IEND2)//'..'//CHAR3(IBEG3:IEND3)
        CALL STRING_TRIM(CHAR4,IBEG4,IEND4)

C---------------------------------------------------------------------

C#### Command: FEM export elements<;FILENAME[default]>
C###  Parameter:      <(geometry/field/material)[geometry]>
C###    Specify the type of elements to export.
C###    `geometry' includes fibres and fields if defined.
C###    `field' exports dependent variable field only.
C###    `material' exports grid point material parameters only.
C###  Parameter:    <elements (GROUP/#s/all)[all]>
C###    Specify either element group, element numbers or all elements
C###    to be exported.
C###  Parameter:    <region (#s/all)[1]>
C###    Limit to elements of specified regions.
C###  Parameter:      <as NAME[0...0]>
C###    Name the element group with a character name.
C###  Parameter:      <cell>
C###    If present in the command, exports spatially varying real
C###    cellular material parameters - Must be used with the `field'
C###    qualifier
C###  Parameter:      <name NAME[fieldname]>
C###    Override the default field name.
C###  Parameter:      <to (datafile/Motif)[datafile]>
C###    Specify the destination file format (`datafile' for exporting to
C###    CMGUI).
C###  Parameter:      <autoname>
C###    Create a CMGUI link compatible file name automatically.
C###  Parameter:      <offset_elem OFFSET[0]>
C####   Add OFFSET to element numbers.
C###  Parameter:      <offset_node OFFSET[offset_elem]>
C###    Add OFFSET to node numbers
C###  Parameter:      <offset_line OFFSET[offset_elem]>
C####   Add OFFSET to line numbers.
C###  Parameter:      <(yq/yqs)[yq]>
C###    For `field' problems using collocation grids specify whether
C###    the YQ or YQS arrays is to be used for the export.
C###  Parameter:      <iy #[1]>
C###    For `field' specify the index in the YP, YQ or YQS dependent
C###    variable array for the desired aspect of the solution.
C###  Parameter:      <nc #[1]>
C###    For `field' specify the type of dependent variable.
C###  Parameter:      <using (fit/solve)[solve]>
C###    For `field' specify the problem type.
C###  Parameter:      <class #[1]>
C###    For `field' specify the class number.
C###  Parameter:      <grid_numbers>
C###    If a grid problem is defined, export grid point numbers as
C###    an element based field with the element geometry. The default
C###    is to not export the grid point numbers. This option does not
C###    work when exporting field values.
C###  Description:
C###    Export elements from CMISS, normally to CMGUI.
CC###    Field information can also be exported. The field is then
CC###    plotted by selecting the plot scalar field option in the
CC###    Graphical element editor.

        OP_STRING(1)=STRING(1:IEND)
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
        OP_STRING(2)=BLANK(1:15)
     '    //'<(geometry/field/material)[geometry]>'
        OP_STRING(3)=BLANK(1:15)//'<element (GROUP/#s)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(5)=BLANK(1:15)//'<as NAME['//CHAR4(IBEG4:IEND4)//']>'
        OP_STRING(6)=BLANK(1:15)//'<to (datafile/Motif)[datafile]>'
        OP_STRING(7)=BLANK(1:15)//'<autoname>'
        OP_STRING(8)=BLANK(1:15)//'<offset_elem OFFSET[0]>'
        OP_STRING(9)=BLANK(1:15)//'<offset_node OFFSET[offset_elem]>'
        OP_STRING(10)=BLANK(1:15)//'<offset_line OFFSET[offset_elem]>'
        OP_STRING(11)=BLANK(1:15)//'<offset_face OFFSET[offset_elem]>'
        OP_STRING(12)=BLANK(1:15)//'<(yq/yqs)[yq]>'
        OP_STRING(13)=BLANK(1:15)//'<iy #[1]>'
        OP_STRING(14)=BLANK(1:15)//'<nc #[1]>'
        OP_STRING(15)=BLANK(1:15)//'<using (fit/solve)[solve]>'
        OP_STRING(16)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(17)=BLANK(1:15)//'<grid_numbers>'
        OP_STRING(18)=BLANK(1:15)//'<cell>'
        OP_STRING(19)=BLANK(1:15)//'<name NAME[fieldname]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe26','doc','EXELEM',ERROR,*9999)
      ELSE
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        nr=NRLIST(1)
        CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '    ERROR,*9999)

        GRID_NUMBERS=.FALSE.

C LKC 13-JUL-1999 adding field naming option
        SET_FIELD_NAME=.FALSE.

C LKC 27-SEP-1999 Initialise
        ELEM_NODES_PREV=0

C RGB initialising string
        SHAPE_STRING = ' '

C LKC 1-JUL-1999 Not initialised - only used for exnode (field headings)
        FIELD_EX_TYPE=' '

C LKC 26-JUN-1999 check there are elements to export
        CALL ASSERT(NELIST(0).GT.0,'>> No Elements to export ',
     '    ERROR,*9999)

C GMH 9/2/97 Allow automatic generation of cmgui link naming convention
        IF(CBBREV(CO,'AUTONAME',4,noco+1,NTCO,N3CO)) THEN
          AUTONAME=.TRUE.
        ELSE
          AUTONAME=.FALSE.
        ENDIF

C GMH 9/2/97 Allowing user to specify IY
        IF(CBBREV(CO,'IY',2,noco+1,NTCO,N3CO)) THEN
          IY=IFROMC(CO(N3CO+1))
        ELSE
          IY=1
        ENDIF

        IF(CBBREV(CO,'YQS',3,noco+1,NTCO,N3CO)) THEN
          GRIDVARIABLE='YQS'
        ELSE
          GRIDVARIABLE='YQ'
        ENDIF

        IF(CBBREV(CO,'GEOMETRY',2,noco+1,NTCO,N3CO)) THEN
          nx=1
          ELEM_TYPE='GEOMETRY'

          IF(CBBREV(CO,'GRID_NUMBERS',4,noco+1,NTCO,N3CO)) THEN
            GRID_NUMBERS=.TRUE.
          ELSE
            GRID_NUMBERS=.FALSE.
          ENDIF

        ELSE IF(CBBREV(CO,'FIELD',2,noco+1,NTCO,N3CO)) THEN
          ELEM_TYPE='FIELD'
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nxc=NXLIST(1)
          IF(CBBREV(CO,'USING',2,noco+1,NTCO,N3CO)) THEN
            IF(ABBREV(CO(N3CO+1),'FIT',2)) THEN
              CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_FIT,ERROR,*9999)
              CALL ASSERT(nx.GT.0,'>>No nx defined for this fit class',
     '          ERROR,*9999)
            ELSE
              CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
              CALL ASSERT(nx.GT.0,'>>No nx defined for this solve '
     '          //'class',ERROR,*9999)
              IF(.NOT.CALL_INIT) THEN
                WRITE(OP_STRING,'('' >>WARNING: Initial conditions have'
     '            //' not been set up!'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
            ENDIF
          ELSE
            CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '        ERROR,*9999)
            IF(.NOT.CALL_INIT) THEN
              WRITE(OP_STRING,'('' >>WARNING: Initial conditions have'
     '          //' not been set up!'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF
          IF(CBBREV(CO,'NAME',3,noco+1,NTCO,N3CO)) THEN
            SET_FIELD_NAME=.TRUE.
            CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
            FIELD_NAME=CO(N3CO+1)(IBEG1:IEND1)
          ELSE
            FIELD_NAME='-'
          ENDIF

        ELSE IF(CBBREV(CO,'MATERIAL',2,noco+1,NTCO,N3CO)) THEN
          nx=1
          ELEM_TYPE='MATERIAL'

        ELSE
          nx=1
          ELEM_TYPE='GEOMETRY'

          IF(CBBREV(CO,'GRID_NUMBERS',4,noco+1,NTCO,N3CO)) THEN
            GRID_NUMBERS=.TRUE.
          ELSE
            GRID_NUMBERS=.FALSE.
          ENDIF

        ENDIF !element type

        IF(CBBREV(CO,'CELL',4,noco+1,NTCO,N3CO)) THEN
          CALL ASSERT(ELEM_TYPE(1:5).EQ.'FIELD',
     '      'Need to also use the field qualifier',ERROR,*9999)
          CALL ASSERT(CALL_CELL,'>>Must define cell first',
     '      ERROR,*9999)
          CALL ASSERT(CALL_CELL_MATE,
     '      '>>Must define cell materials first',ERROR,*9999)
          EXPORT_CELL=.TRUE.
        ELSE
          EXPORT_CELL=.FALSE.
        ENDIF

        IF(CBBREV(CO,'NC',2,noco+1,NTCO,N3CO)) THEN
          nc=IFROMC(CO(N3CO+1))
        ELSE
          nc=1
        ENDIF

C news MPN 9-11-95: OFFSET option to add const to all elems/faces/lines
c cpb 14/8/95 Extending offset with different #s for ne/np/nl and nf
        IF(CBBREV(CO,'OFFSET_ELEM',8,noco+1,NTCO,N3CO)) THEN
          offset_elem=IFROMC(CO(N3CO+1))
        ELSE
          offset_elem=0
        ENDIF
        IF(CBBREV(CO,'OFFSET_NODE',8,noco+1,NTCO,N3CO)) THEN
          offset_node=IFROMC(CO(N3CO+1))
        ELSE
          offset_node=offset_elem
        ENDIF
        IF(CBBREV(CO,'OFFSET_LINE',8,noco+1,NTCO,N3CO)) THEN
          offset_line=IFROMC(CO(N3CO+1))
        ELSE
          offset_line=offset_elem
        ENDIF
        IF(CBBREV(CO,'OFFSET_FACE',8,noco+1,NTCO,N3CO)) THEN
          offset_face=IFROMC(CO(N3CO+1))
        ELSE
          offset_face=offset_elem
        ENDIF

C news MPN 25-Jul-95: NAME option for exelem file
        IF(CBBREV(CO,'AS',1,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          ELEM_NAME=CO(N3CO+1)(IBEG1:IEND1)
        ELSE
          IF(AUTONAME) THEN
            WRITE(CHAR1,'(I5)') nr
            CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
            ELEM_NAME='region_'//CHAR1(IBEG1:IEND1)
          ELSE !autoname
            WRITE(CHAR1,'(I5)') NELIST(1)+offset_elem
            WRITE(CHAR2,'(I5)') NELIST(NELIST(0))+offset_elem
            CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
            CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
            ELEM_NAME=CHAR1(IBEG1:IEND1)//'..'//CHAR2(IBEG2:IEND2)
          ENDIF !autoname
        ENDIF

        IF(CBBREV(CO,'TO',1,noco+1,NTCO,N3CO)) THEN
          IF(CBBREV(CO,'DATAFILE',2,noco+1,NTCO,N3CO)) THEN
            OUTPUT='DATAFILE'
            DATAFILE=.TRUE.
C KAT 15Feb99: SPREADSHEET not available:
C          ELSE IF(CBBREV(CO,'SPREADSHEET',2,noco+1,NTCO,N3CO)) THEN
C            OUTPUT='SPREADSHEET'
C            DATAFILE=.FALSE.
          ELSE IF(CBBREV(CO,'MOTIF',2,noco+1,NTCO,N3CO)) THEN
            OUTPUT='MOTIF'
            DATAFILE=.FALSE.
          ELSE IF(CBBREV(CO,'EXPLORER',2,noco+1,NTCO,N3CO)) THEN
            OUTPUT='EXPLORER'
            DATAFILE=.FALSE.
          ELSE
            OUTPUT='DATAFILE'
            DATAFILE=.TRUE.
          ENDIF
        ELSE
          OUTPUT='DATAFILE'
          DATAFILE=.TRUE.
        ENDIF

C GMH 26/12/96 Explicitly intitialise FIELD_NODE_LIST
        DO np=0,MAX_ELEMENT_NODES
          FIELD_NODE_LIST(np)=0
        ENDDO !np

        na=1 !temporary

        IF((OUTPUT(1:8).EQ.'DATAFILE').OR.(OUTPUT(1:5).EQ.'MOTIF')) THEN
          CALL CHECKF(1,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          IF(NELIST(0).GT.0) THEN
            IF(DATAFILE) THEN
              CALL STRING_TRIM(FILE,IBEG,IEND)
              CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.exelem','NEW',
     '          'SEQUEN','FORMATTED',132,ERROR,*9999)
C**             write the group name
              CALL STRING_TRIM(ELEM_NAME,IBEG,IEND)
              WRITE(IFILE,'( '' Group name: '',A)') ELEM_NAME(IBEG:IEND)
            ELSE
              IF(USE_SOCKET) THEN
C**             send the data type identifier
                IF(FSKWRITE(DATA_TYPE,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
C**             send the group name
                CALL STRING_TRIM(ELEM_NAME,IBEG,IEND)
                CLEN=FSKLEN(ELEM_NAME(IBEG:IEND))
                CALL FSKF2C(ELEM_NAME(IBEG:IEND),CLEN,INTSTR)
                IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
                IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
     '            GOTO 9999
              ELSE
                ERROR=
     '            '>Can only export to MOTIF if sockets are being used'
                GOTO 9999
              ENDIF
            ENDIF !datafile/socket

C**         zero the array for storing the faces used
            DO nf=0,NFT
              NFLIST(nf)=0
            ENDDO

C**         zero the array for storing the lines used
            DO nl=0,NLT
              NLLIST(nl)=0
            ENDDO

C**         determine the faces and lines used
            DO no_nelist=1,NELIST(0)
              NE=NELIST(no_nelist)
              IF(NJ_LOC(NJL_GEOM,0,nr).GT.0) THEN
                nb=NBJ(1,NE)
                IF(NIT(nb).GE.3) THEN
C**               add faces to the faces list
                  DO nf_e=1,NFE(nb)
                    nf=NFF(nf_e,NE)
                    IF(nf.NE.0) THEN
                      IF(NFLIST(nf).EQ.0) THEN
                        NFLIST(nf)=1
                        NFLIST(0)=NFLIST(0)+1
                      ENDIF
                    ENDIF
                  ENDDO !nf_e
                ELSE IF(NIT(nb).GE.2) THEN
C**               add lines to the lines list
                  DO NAE=1,NLE(nb)
                    nl=NLL(NAE,NE)
                    IF(nl.GT.0) THEN
                      IF(NLLIST(nl).EQ.0) THEN
                        NLLIST(nl)=1
                        NLLIST(0)=NLLIST(0)+1
                      ENDIF
                    ENDIF !nl.GT.0
                  ENDDO !nae
                ENDIF
              ENDIF
            ENDDO !no_nelist (ne)

            IF(NFLIST(0).GT.0) THEN
              DO nf=1,NFT
                IF(NFLIST(nf).NE.0) THEN
                  nb=NBJF(1,NF)
C**               add lines to the lines list
                  DO NAF=1,NLE(nb)
                    nl=NLF(NAF,NF)
                    IF(nl.GT.0) THEN
                      IF(NLLIST(nl).EQ.0) THEN
                        NLLIST(nl)=1
                        NLLIST(0)=NLLIST(0)+1
                      ENDIF
                    ENDIF !nl.GT.0
                  ENDDO !naf
                ENDIF
              ENDDO !nf
            ENDIF !NFLIST(0)>0

C**         write the lines
            IF(DATAFILE) THEN
              IF(NLLIST(0).GT.0) THEN
                WRITE(IFILE,'( '' Shape.  Dimension=1'' )')
                DO nl=1,NLT

                  IF(NLLIST(nl).NE.0) THEN
C LKC 26-AUG-2002 Need to check for more than 100000 lines
C                   WRITE(IFILE,'( '' Element: 0 0 '',I5)')
C                   '                nl+offset_line
                    IF(nl+offset_line.LT.100000) THEN
                      FORMAT='( '' Element: 0 0 '',I5)'
                    ELSE
                      FORMAT='( '' Element: 0 0 '',I9)'
                    ENDIF
                    WRITE(IFILE,FORMAT) nl+offset_line
                  ENDIF

                ENDDO !nl
              ENDIF
            ELSE
              IF(NLLIST(0).GT.0) THEN
                IF(FSKWRITE(ICHAR('S'),SK_CHAR,1,CONNID2).EQ.-1)
     '            GOTO 9999
                ISOCKET(1)=1
                IF(FSKWRITE(ISOCKET,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
                ISOCKET(1)=3
                DO nl=1,NLT
                  IF(NLLIST(nl).NE.0) THEN
                    IF(FSKWRITE(ICHAR('E'),SK_CHAR,1,CONNID2).EQ.-1)
     '                GOTO 9999
                    ISOCKET(2)=nl
                    IF(FSKWRITE(ISOCKET,SK_LONG_INT,2,CONNID2).EQ.-1)
     '                GOTO 9999
                  ENDIF
                ENDDO !nl
              ENDIF !NLLIST(0)>0
            ENDIF !datafile/socket

C**         write the faces
            IF(DATAFILE) THEN
              IF(NFLIST(0).GT.0) THEN
                IF(IBT(1,1,nb).EQ.3.AND.IBT(3,1,nb).EQ.2) THEN
                  WRITE(IFILE,'( '' Shape.  Dimension=2, '
     '              //'simplex(2)*simplex'' )')
                ELSE
                  WRITE(IFILE,'( '' Shape.  Dimension=2'' )')
                ENDIF
                DO nf=1,NFT
                  IF(NFLIST(nf).NE.0) THEN
                    WRITE(IFILE,'( '' Element: 0 '',I5,'' 0'' )')
     '                nf+offset_face
                    WRITE(IFILE,'( ''   Faces:'' )')
                    nb=NBJF(1,nf)
C**                 write the lines for the face
c cpb 3/8/95 Adding sector export
C                    CALL ASSERT(NLE(nb).EQ.4,
C     '                ' Incorrect #lines for face',ERROR,*9999)
                    CALL ASSERT(NLE(nb).GT.2,
     '                ' Incorrect #lines for face',ERROR,*9999)

                    IF(NLE(nb).EQ.4) THEN

                      LARGE_FORMAT=.FALSE.

C KAT 2Dec98: no offset for no line
                      DO naf=1,4
                        nl=NLF(naf,nf)
                        IF(nl.NE.0) THEN
                          NLFF(naf)=nl+offset_line
                          IF(NLFF(naf).GT.99999) LARGE_FORMAT=.TRUE.
                        ELSE
                          NLFF(naf)=0
                        ENDIF
                      ENDDO !naf
C LKC 26-AUG-2002 Need to check for more then 100000 faces
C                      WRITE(IFILE,'( ''   0 0 '',I5)') NLFF(3)
C                      WRITE(IFILE,'( ''   0 0 '',I5)') NLFF(4)
C                      WRITE(IFILE,'( ''   0 0 '',I5)') NLFF(1)
C                      WRITE(IFILE,'( ''   0 0 '',I5)') NLFF(2)

                      IF(LARGE_FORMAT) THEN
                        FORMAT='( ''   0 0 '',I9)'
                      ELSE
                        FORMAT='( ''   0 0 '',I5)'
                      ENDIF

                      WRITE(IFILE,FORMAT) NLFF(3)
                      WRITE(IFILE,FORMAT) NLFF(4)
                      WRITE(IFILE,FORMAT) NLFF(1)
                      WRITE(IFILE,FORMAT) NLFF(2)

                    ELSE
C DB 4Jun01 nb is already the basis for the face, so shouldn't use npf
C                      ni1=NPF(1,nf)
C                      ni2=NPF(3,nf)
                      ni1=1
                      ni2=2

C LKC 26-AUG-2002 Need to check for more than 100000 lines
                      LARGE_FORMAT=.FALSE.
                      DO naf=1,3
                        IF(NLF(naf,nf).GE.100000) LARGE_FORMAT=.TRUE.
                      ENDDO

                      IF(LARGE_FORMAT) THEN
                        FORMAT='( ''   0 0 '',I9)'
                      ELSE
                        FORMAT='( ''   0 0 '',I5)'
                      ENDIF

                      IF(IBT(1,ni1,nb).EQ.5) THEN
                        WRITE(IFILE,FORMAT) NLF(2,nf)+offset_line
                        WRITE(IFILE,FORMAT) NLF(3,nf)+offset_line
                        WRITE(IFILE,'( ''   0 0     0'')')
                        WRITE(IFILE,FORMAT) NLF(1,nf)+offset_line
                      ELSE IF(IBT(1,ni1,nb).EQ.6) THEN
                        WRITE(IFILE,FORMAT) NLF(2,nf)+offset_line
                        WRITE(IFILE,FORMAT) NLF(3,nf)+offset_line
                        WRITE(IFILE,FORMAT) NLF(1,nf)+offset_line
                        WRITE(IFILE,'( ''   0 0     0'')')
                      ELSE IF(IBT(1,ni2,nb).EQ.5) THEN
                        WRITE(IFILE,'( ''   0 0     0'')')
                        WRITE(IFILE,FORMAT) NLF(3,nf)+offset_line
                        WRITE(IFILE,FORMAT) NLF(1,nf)+offset_line
                        WRITE(IFILE,FORMAT) NLF(2,nf)+offset_line
                      ELSE IF(IBT(1,ni2,nb).EQ.6) THEN
                        WRITE(IFILE,FORMAT) NLF(3,nf)+offset_line
                        WRITE(IFILE,'( ''   0 0     0'')')
                        WRITE(IFILE,FORMAT) NLF(1,nf)+offset_line
                        WRITE(IFILE,FORMAT) NLF(2,nf)+offset_line
                      ELSE IF(IBT(1,1,nb).EQ.3.AND.IBT(2,1,nb).EQ.4)
     '                    THEN !Hermite simplex
                        IF(NKT(1,nb).EQ.1) THEN !apex node 1
                          WRITE(IFILE,FORMAT) NLF(1,nf)+offset_line
                          WRITE(IFILE,FORMAT) NLF(2,nf)+offset_line
                          WRITE(IFILE,'( ''   0 0 0'')')
                          WRITE(IFILE,FORMAT) NLF(3,nf)+offset_line
                        ELSE IF(NKT(3,nb).EQ.1) THEN !apex node 3
                          WRITE(IFILE,FORMAT) NLF(2,nf)+offset_line
                          WRITE(IFILE,FORMAT) NLF(3,nf)+offset_line
                          WRITE(IFILE,FORMAT) NLF(1,nf)+offset_line
                          WRITE(IFILE,'( ''   0 0 0'')')
                        ENDIF
                      ELSE IF(IBT(1,ni1,nb).EQ.3.AND.IBT(3,1,nb).EQ.2)
     '                    THEN
                        WRITE(IFILE,FORMAT) NLF(1,nf)+offset_line
                        WRITE(IFILE,FORMAT) NLF(2,nf)+offset_line
                        WRITE(IFILE,FORMAT) NLF(3,nf)+offset_line
                      ELSE
                        ERROR='>>Unknown sector type'
                        GOTO 9999
                      ENDIF
                    ENDIF !NLE(nb)
                  ENDIF !NFLIST(nf)/=0
                ENDDO !nf
              ENDIF !NFLIST(0)>0

            ELSE !socket
              IF(NFLIST(0).GT.0) THEN
                IF(FSKWRITE(ICHAR('S'),SK_CHAR,1,CONNID2).EQ.-1)
     '            GOTO 9999
                ISOCKET(1)=2
                IF(FSKWRITE(ISOCKET,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
                ISOCKET(1)=2
                ISOCKET(3)=3
                ISOCKET(5)=3
                ISOCKET(7)=3
                ISOCKET(9)=3
                DO nf=1,NFT
                  IF(NFLIST(nf).NE.0) THEN
                    IF(FSKWRITE(ICHAR('E'),SK_CHAR,1,CONNID2).EQ.-1)
     '                GOTO 9999
                    ISOCKET(2)=nf
                    nb=NBJF(1,nf)
C**                 write the lines for the face
c cpb 3/8/95 Adding sector export
C                    CALL ASSERT(NLE(nb).EQ.4,
C     '                ' Incorrect #lines for face',ERROR,*9999)
                    CALL ASSERT(NLE(nb).GT.2,
     '                ' Incorrect #lines for face',ERROR,*9999)
                    IF(NLE(nb).EQ.4) THEN
                      ISOCKET(4)=NLF(3,nf)
                      ISOCKET(6)=NLF(4,nf)
                      ISOCKET(8)=NLF(1,nf)
                      ISOCKET(10)=NLF(2,nf)
                    ELSE
C DB 4Jun01 nb is already the basis for the face, so shouldn't use npf
C                      ni1=NPF(1,nf)
C                      ni2=NPF(3,nf)
                      ni1=1
                      ni2=2
                      IF(IBT(1,ni1,nb).EQ.5) THEN
                        ISOCKET(4)=NLF(3,nf)
                        ISOCKET(6)=NLF(4,nf)
                        ISOCKET(8)=0
                        ISOCKET(10)=NLF(2,nf)
                      ELSE IF(IBT(1,ni1,nb).EQ.6) THEN
                        ISOCKET(4)=NLF(3,nf)
                        ISOCKET(6)=NLF(4,nf)
                        ISOCKET(8)=NLF(1,nf)
                        ISOCKET(10)=0
                      ELSE IF(IBT(1,ni2,nb).EQ.5) THEN
                        ISOCKET(4)=0
                        ISOCKET(6)=NLF(4,nf)
                        ISOCKET(8)=NLF(1,nf)
                        ISOCKET(10)=NLF(2,nf)
                      ELSE IF(IBT(1,ni2,nb).EQ.6) THEN
                        ISOCKET(4)=NLF(3,nf)
                        ISOCKET(6)=0
                        ISOCKET(8)=NLF(1,nf)
                        ISOCKET(10)=NLF(2,nf)
                      ELSE IF(IBT(1,1,nb).EQ.3.AND.IBT(1,2,nb).EQ.4)
     '                    THEN !Hermite simplex
                        IF(NKT(1,nb).EQ.1) THEN !Apex node 1
                          ISOCKET(4)=NLF(1,nf)
                          ISOCKET(6)=NLF(2,nf)
                          ISOCKET(8)=0
                          ISOCKET(10)=NLF(3,nf)
                        ELSE IF(NKT(3,nb).EQ.1) THEN !Apex node 3
                          ISOCKET(4)=NLF(2,nf)
                          ISOCKET(6)=NLF(3,nf)
                          ISOCKET(8)=NLF(1,nf)
                          ISOCKET(10)=0
                        ENDIF
                      ELSE IF(IBT(1,ni1,nb).EQ.3.AND.IBT(3,1,nb).EQ.2)
     '                    THEN
                          ISOCKET(4)=NLF(1,nf)
                          ISOCKET(6)=NLF(2,nf)
                          ISOCKET(8)=NLF(3,nf)
                      ELSE
                        ERROR='>>Unknown sector type'
                        GOTO 9999
                      ENDIF
                    ENDIF
                    IF(FSKWRITE(ISOCKET,SK_LONG_INT,10,CONNID2).EQ.-1)
     '                GOTO 9999
                  ENDIF !NFLIST(nf)/=0
                ENDDO !nf
              ENDIF !NFLIST(0)>0
            ENDIF !datafile/socket

C**         write the elements
            FIRST_ELEMENT=.TRUE.
            NE_FIELD=NELIST(1)
            NE_SHAPE=NELIST(1)

            nqsc_old=-1

C---------------------------------------------------------------------

            IF((ELEM_TYPE(1:8).EQ.'GEOMETRY').OR.
     '        (ELEM_TYPE(1:5).EQ.'FIELD')) THEN
              DO no_nelist=1,NELIST(0)
                NE=NELIST(no_nelist)
                IF(ELEM_TYPE(1:8).EQ.'GEOMETRY') THEN
                  NJHT=NJ_LOC(NJL_GEOM,0,nr)+NJ_LOC(NJL_FIBR,0,nr)+
     '              NJ_LOC(NJL_FIEL,0,nr)
                  NJH_COUNT=0
                  DO njh=1,NJ_LOC(NJL_GEOM,0,nr)
                    NJH_COUNT=NJH_COUNT+1
                    NJH_LIST(NJH_COUNT)=NJ_LOC(NJL_GEOM,njh,nr)
                  ENDDO
                  DO njh=1,NJ_LOC(NJL_FIBR,0,nr)
                    NJH_COUNT=NJH_COUNT+1
                    NJH_LIST(NJH_COUNT)=NJ_LOC(NJL_FIBR,njh,nr)
                  ENDDO
                  DO njh=1,NJ_LOC(NJL_FIEL,0,nr)
                    NJH_COUNT=NJH_COUNT+1
                    NJH_LIST(NJH_COUNT)=NJ_LOC(NJL_FIEL,njh,nr)
                  ENDDO
C maps njh to nj and nj to njh NPS 21/11/97
                  NJHT=NJH_COUNT
                ELSE !field

C ***
C *** DPN - NHE(ne,nx) is the number of dependent variables defined in
C ***       element ne for problem type nx
C ***
C ???     - NJHT is the number of fields ??? i.e. the number of
C ???       parameters to export ?????
C ***

                  IF(EXPORT_CELL) THEN
                    NJHT=1
c                    NJHT=IRCQS_SPATIAL(0)
c                    DO njh=1,NJHT
c                      NJH_LIST(njh)=IRCQS_SPATIAL(njh)
c                    ENDDO
                  ELSE
                    NJHT=NHE(ne,nx)
C SGM 26Oct2000 grid-based Finite element also
C MLT 29Nov02 Grid Finite Volume also
                    IF((ITYP4(nr,nx).EQ.4.OR.ITYP4(nr,nx).EQ.6.OR.
     '                  ITYP4(nr,nx).EQ.7).AND.
     '                (NJHT.EQ.0)) NJHT=1
                    DO njh=1,NJHT
                      NJH_LIST(njh)=njh
                    ENDDO
                  ENDIF
                ENDIF !geometry/field

                IF(NJHT.GT.0) THEN
C**               check if the shape has changed
                  SHAPE_CHANGED=.FALSE.
                  FIELDS_CHANGED=.FALSE.
                  IF(ELEM_TYPE(1:8).EQ.'GEOMETRY') THEN
                    IF(FIRST_ELEMENT.OR.
     '                (NIT(NBJ(1,NE_SHAPE)).NE.NIT(NBJ(1,NE))))
     '                THEN
                      SHAPE_CHANGED=.TRUE.
                      FIELDS_CHANGED=.TRUE.
                    ENDIF
                  ELSE !ELEM_TYPE='FIELD'
                    IF(ITYP19(nr,nx).EQ.2.AND.ITYP3(nr,nx).EQ.4) THEN !Infarct
                      IF(.NOT.EXPORT_CELL.AND.FIRST_ELEMENT) THEN
                        SHAPE_CHANGED=.TRUE.
                        FIELDS_CHANGED=.TRUE.
                      ELSEIF(EXPORT_CELL.AND.FIRST_ELEMENT) THEN
                        SHAPE_CHANGED=.TRUE.
                        FIELDS_CHANGED=.TRUE.
                      ENDIF
                    ELSE
                      IF(.NOT.EXPORT_CELL.AND.(FIRST_ELEMENT.OR.
     '                  (NIT(NBH(NH_LOC(1,nx),1,NE_SHAPE)).NE.
     '                  NIT(NBH(NH_LOC(1,nx),1,NE))))) THEN
                        SHAPE_CHANGED=.TRUE.
                        FIELDS_CHANGED=.TRUE.
                      ELSEIF(EXPORT_CELL.AND.FIRST_ELEMENT) THEN
                        SHAPE_CHANGED=.TRUE.
                        FIELDS_CHANGED=.TRUE.
                      ENDIF
                    ENDIF !Infarct
                  ENDIF !ELEM_TYPE
C**               calculate the element node list and check if the
C**               fields have changed
                  DO nb=0,NBFT
                    NBLIST(nb)=0
                  ENDDO
                  ELEMENT_NODE_LIST(0)=0
                  DO no_njh=1,NJHT
                    njh=NJH_LIST(no_njh)
                    IF(ELEM_TYPE(1:8).EQ.'GEOMETRY') THEN
C**                     reorder nj's

                      IF((1.LE.no_njh).AND.
     '                  (no_njh.LE.NJ_LOC(NJL_GEOM,0,nr))) THEN
                        NJ=NJ_LOC(NJL_GEOM,no_njh,nr)
                      ELSE IF((NJ_LOC(NJL_GEOM,0,nr).LT.no_njh).AND.
     '                    (no_njh.LE.NJ_LOC(NJL_GEOM,0,nr)+
     '                    NJ_LOC(NJL_FIBR,0,nr))) THEN
                        NJ=
     '                    NJ_LOC(NJL_FIBR,no_njh-
     '                    NJ_LOC(NJL_GEOM,0,nr),nr)
                      ELSE IF((NJ_LOC(NJL_GEOM,0,nr)+
     '                    NJ_LOC(NJL_FIBR,0,nr)
     '                    .LT.no_njh).AND.(no_njh.LE.
     '                    NJ_LOC(NJL_GEOM,0,nr)+
     '                    NJ_LOC(NJL_FIBR,0,nr)+
     '                    NJ_LOC(NJL_FIEL,0,nr))) THEN
                        NJ=NJ_LOC(NJL_FIEL,no_njh-NJ_LOC(NJL_FIBR,0,nr)-
     '                    NJ_LOC(NJL_GEOM,0,nr),nr)
                      ELSE
                        ERROR='>>Invalid nj'
                        GOTO 9999
                      ENDIF

                      nb=NBJ(nj,NE)
                      IF(nb.EQ.NBJ(njh,NE_FIELD)) THEN
                        nn=NNT(nb)
                        njTOT_NE=0
                        njTOT_NE_FIELD=0
                        DO njj1=1,3
                          DO njj2=1,NJ_LOC(njj1,0,nr)
                            nj=NJ_LOC(njj1,njj2,nr)
                            IF(NVJP(nj,NPNE(nn,nb,NE)).GT.0)
     '                        njTOT_NE=njTOT_NE+1
                            IF(NVJP(nj,NPNE(nn,nb,NE_FIELD)).GT.0)
     '                        njTOT_NE_FIELD=njTOT_NE_FIELD+1
                          ENDDO !njj2
                        ENDDO !njj1
                        DO WHILE((nn.GT.1).AND.(NVJE(nn,nb,njh,ne).EQ.
     '                    NVJE(nn,nb,njh,NE_FIELD).AND.
     '                    (njTOT_NE.EQ.njTOT_NE_FIELD)))
                          nn=nn-1
                          njTOT_NE=0
                          njTOT_NE_FIELD=0
                          DO njj1=1,3
                            DO njj2=1,NJ_LOC(njj1,0,nr)
                              nj=NJ_LOC(njj1,njj2,nr)
                              IF(NVJP(nj,NPNE(nn,nb,NE)).GT.0)
     '                          njTOT_NE=njTOT_NE+1
                              IF(NVJP(nj,NPNE(nn,nb,NE_FIELD)).GT.0)
     '                          njTOT_NE_FIELD=njTOT_NE_FIELD+1
                            ENDDO !njj2
                          ENDDO !njj1
                        ENDDO !nn
                        IF(NVJE(nn,nb,njh,ne).NE.
     '                    NVJE(nn,nb,njh,NE_FIELD).OR.
     '                    njTOT_NE.NE.njTOT_NE_FIELD) THEN
                          FIELDS_CHANGED=.TRUE.
                        ENDIF
                      ELSE
                        FIELDS_CHANGED=.TRUE.
                      ENDIF
                    ELSE
                      IF(EXPORT_CELL) THEN
                        ! ???? DPN ?????
                        nb=NBJ(1,ne)
                      ELSE
                        nb=NBH(NH_LOC(njh,nx),1,NE)
                      ENDIF
C SGM 26Oct2000 grid-based Finite element also
C MLT 29Nov02 Grid Finite Volume also
                      IF(ITYP4(nr,nx).EQ.4.OR.ITYP4(nr,nx).EQ.6.OR.
     '                  ITYP4(nr,nx).EQ.7) nb=NBJ(1,ne) !collocation
                      IF(EXPORT_CELL) THEN
                        CELL_TEST=.FALSE.
                      ELSE
                        CELL_TEST=nb.EQ.NBH(NH_LOC(njh,nx),1,NE_FIELD)
                      ENDIF
                      IF(CELL_TEST) THEN

                        IF(NBC(nb).NE.0) THEN !Stops auxilliary basis export
                          nn=NNT(nb)
                          DO WHILE(nn.GT.1.AND.
     '                      NVHE(nn,nb,NH_LOC(njh,nx),ne).EQ.
     '                      NVHE(nn,nb,NH_LOC(njh,nx),NE_FIELD).AND.
     '                      NHP(NPNE(nn,nb,NE),nr,nx).EQ.
     '                    NHP(NPNE(nn,nb,NE_FIELD),nr,nx))
                            nn=nn-1
                          ENDDO !nn
                          IF(NVHE(nn,nb,NH_LOC(njh,nx),ne).NE.
     '                      NVHE(nn,nb,NH_LOC(njh,nx),NE_FIELD).OR.
     '                      NHP(NPNE(nn,nb,NE),nr,nx).NE.
     '                      NHP(NPNE(nn,nb,NE_FIELD),nr,nx)) THEN
                            FIELDS_CHANGED=.TRUE.
                          ENDIF
                        ENDIF
                      ELSE !basis has changed for nj
                        FIELDS_CHANGED=.TRUE.
                      ENDIF
                    ENDIF
C!!! DB.  Temporary.  Stops auxilliary basis export
                    IF(NBC(nb).NE.0) THEN
                      IF(NBLIST(nb).EQ.0) THEN
C**                       this basis not previously found
                        NBLIST(nb)=1
                        NBLIST(0)=NBLIST(0)+1
                      ENDIF
                      DO nn=1,NNT(nb)
                        NP=NPNE(nn,nb,NE)
                        NP_ELEMENT=ELEMENT_NODE_LIST(0)
                        DO WHILE ((NP_ELEMENT.GT.0).AND.
     '                    (NP.NE.ELEMENT_NODE_LIST(NP_ELEMENT)))
                          NP_ELEMENT=NP_ELEMENT-1
                        ENDDO !while (NP_ELEMENT.GT.0)
                        IF(NP_ELEMENT.LE.0) THEN
                          CALL ASSERT(ELEMENT_NODE_LIST(0)+1.LE.
     '                      MAX_ELEMENT_NODES,
     '                      '>>Too many local nodes',
     '                      ERROR,*9999)
                          ELEMENT_NODE_LIST(0)=ELEMENT_NODE_LIST(0)+1
                          ELEMENT_NODE_LIST(ELEMENT_NODE_LIST(0))=NP
                          NP_ELEMENT=ELEMENT_NODE_LIST(0)
                        ENDIF
C!!! LKC 19-MAY-1999 is this wierd ?
                        IF(NPNE(nn,nb,NE_FIELD).NE.
     '                    FIELD_NODE_LIST(NP_ELEMENT)) THEN
                          FIELDS_CHANGED=.TRUE.
                        ENDIF
                      ENDDO !nn

C LKC 19-MAY-1999 fix for repeated nodes in an element
C eg for the case with a three noded bilinear element
                      IF(ELEMENT_NODE_LIST(0).NE.ELEM_NODES_PREV)
     '                  THEN
                        FIELDS_CHANGED=.TRUE.
                      ENDIF
                      ELEM_NODES_PREV=ELEMENT_NODE_LIST(0)
                    ENDIF
                  ENDDO !njh

                  IF(ELEM_TYPE(1:5).EQ.'FIELD'.AND.
     '              (ITYP2(nr,nx).EQ.2.OR.
     '              ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3)) THEN
C                   For Finite elasticity or fluid mechanics/
C                   const vol constraint, check if deformed fibre
C                   fields have changed
                    DO njj=1,NJ_LOC(NJL_FIBR,0,nr)
                      njh=NJ_LOC(NJL_FIBR,njj,nr)
                      nb=NBJ(njh,NE)
                      IF(NBLIST(nb).EQ.0) THEN
C**                     this basis not previously found
                        NBLIST(nb)=1
                        NBLIST(0)=NBLIST(0)+1
                      ENDIF
                      IF(nb.EQ.NBJ(njh,NE_FIELD)) THEN
                        nn=NNT(nb)
                        njTOT_NE=0
                        njTOT_NE_FIELD=0
                        DO njj1=1,3
                          DO njj2=1,NJ_LOC(njj1,0,nr)
                            nj=NJ_LOC(njj1,njj2,nr)
                            IF(NVJP(nj,NPNE(nn,nb,NE)).GT.0)
     '                        njTOT_NE=njTOT_NE+1
                            IF(NVJP(nj,NPNE(nn,nb,NE_FIELD)).GT.0)
     '                        njTOT_NE_FIELD=njTOT_NE_FIELD+1
                          ENDDO !njj2
                        ENDDO !njj1
                        DO WHILE((nn.GT.1).AND.
     '                    (NVJE(nn,nb,njh,ne).EQ.
     '                    NVJE(nn,nb,njh,NE_FIELD).AND.
     '                    (njTOT_NE.EQ.njTOT_NE_FIELD)))
                          nn=nn-1
                          njTOT_NE=0
                          njTOT_NE_FIELD=0
                          DO njj1=1,3
                            DO njj2=1,NJ_LOC(njj1,0,nr)
                              nj=NJ_LOC(njj1,njj2,nr)
                              IF(NVJP(nj,NPNE(nn,nb,NE)).GT.0)
     '                          njTOT_NE=njTOT_NE+1
                              IF(NVJP(nj,NPNE(nn,nb,NE_FIELD)).GT.0)
     '                          njTOT_NE_FIELD=njTOT_NE_FIELD+1
                            ENDDO !njj2
                          ENDDO !njj1
                        ENDDO !nn
                        IF(NVJE(nn,nb,njh,ne).NE.
     '                    NVJE(nn,nb,njh,NE_FIELD).OR.
     '                    njTOT_NE.NE.njTOT_NE_FIELD) THEN
                          FIELDS_CHANGED=.TRUE.
                        ENDIF
                      ELSE !basis has changed for nj
                        FIELDS_CHANGED=.TRUE.
                      ENDIF
                    ENDDO !njj
                  ENDIF !FIELD and finite elasticity

                  IF(SHAPE_CHANGED) THEN
                    IF(ELEM_TYPE(1:8).EQ.'GEOMETRY') THEN
                      ELEMENT_DIMENSION=NIT(NBJ(1,NE))
                    ELSE
                      IF(EXPORT_CELL) THEN
                        ELEMENT_DIMENSION=NIT(NBJ(1,NE))
                      ELSE
C!!! KAT 2007-01-30: If field is grid based then there may not be a
C!!!                 dependent variable basis function, and a geometry
C!!!                 basis function would be suitable.
                        ELEMENT_DIMENSION=NIT(NBH(NH_LOC(1,nx),1,NE))
                      ENDIF
                    ENDIF
                    IF(DATAFILE) THEN
C New Simplex elements CS 12/11/97
                      SIMPLEX_BASIS_DIRNS(0)=0
                      SIMPLEX_BASIS_DIRNS(1)=0
                      DO ni=1,ELEMENT_DIMENSION
                        IF(IBT(1,ni,nb).EQ.3) THEN
C LKC 20-NOV-97 not for hermite simplex
                          IF(IBT(2,ni,nb).NE.4) THEN
                            SIMPLEX_BASIS_DIRNS(0)=
     '                        SIMPLEX_BASIS_DIRNS(0)+1
                            SIMPLEX_BASIS_DIRNS(SIMPLEX_BASIS_DIRNS(0))
     '                        =ni
                          ENDIF
                        ENDIF
                      ENDDO
                      IF(SIMPLEX_BASIS_DIRNS(0).NE.0) THEN
                         SIMPLEX_BASIS_LIST='('
                         DO i=1,SIMPLEX_BASIS_DIRNS(0) - 1
                           WRITE(CHAR1,'(I1)') SIMPLEX_BASIS_DIRNS(i+1)
                           CALL STRING_TRIM(SIMPLEX_BASIS_LIST,IBEG,
     '                       IEND)
                           SIMPLEX_BASIS_LIST=SIMPLEX_BASIS_LIST(
     '                       IBEG:IEND)//CHAR1
                           CALL STRING_TRIM(SIMPLEX_BASIS_LIST,IBEG,
     '                       IEND)
                          SIMPLEX_BASIS_LIST=SIMPLEX_BASIS_LIST(
     '                       IBEG:IEND)//';'
                        ENDDO
                        CALL STRING_TRIM(SIMPLEX_BASIS_LIST,IBEG2,IEND2)
                        SIMPLEX_BASIS_LIST=SIMPLEX_BASIS_LIST(
     '                    IBEG2:IEND2-1)//')'
                        DO ni=1,ELEMENT_DIMENSION
                          IF(IBT(1,ni,nb).EQ.3) THEN
                            CALL STRING_TRIM(SHAPE_STRING,IBEG,IEND)
                            SHAPE_STRING=
     '                        SHAPE_STRING(IBEG:IEND)//'simplex'
                          ELSE
                            CALL STRING_TRIM(SHAPE_STRING,IBEG,IEND)
                            SHAPE_STRING=SHAPE_STRING(IBEG:IEND)//'line'
                          ENDIF
                          IF(ni.EQ.SIMPLEX_BASIS_DIRNS(1)) THEN
                            CALL STRING_TRIM(SHAPE_STRING,IBEG,IEND)
                            SHAPE_STRING=SHAPE_STRING(IBEG:IEND)
     '                        //SIMPLEX_BASIS_LIST (IBEG2:IEND2)
                          ENDIF
                          IF(ni.LT.ELEMENT_DIMENSION) THEN
                            CALL STRING_TRIM(SHAPE_STRING,IBEG,IEND)
                            SHAPE_STRING=SHAPE_STRING(IBEG:IEND)//'*'
                          ENDIF
                        ENDDO
                      ENDIF
                      IF(SIMPLEX_BASIS_DIRNS(0).NE.0) THEN
                        CALL STRING_TRIM(SHAPE_STRING,IBEG,IEND)
                        WRITE(IFILE,'( '' Shape.  Dimension='',I1,'', '
     '                    //''',A)')
     '                    ELEMENT_DIMENSION,SHAPE_STRING(IBEG:IEND)
                      ELSE
                        WRITE(IFILE,'( '' Shape.  Dimension='',I1)')
     '                    ELEMENT_DIMENSION
                      ENDIF
                    ELSE
                      IF(FSKWRITE(ICHAR('S'),SK_CHAR,1,CONNID2).EQ.-1)
     '                  GOTO 9999
                      IF(FSKWRITE(ELEMENT_DIMENSION,SK_LONG_INT,1,
     '                  CONNID2).EQ.-1) GOTO 9999
                    ENDIF
                    NE_SHAPE=NE
                  ENDIF
C CPB 26/11/98 Check if the grid scheme (if any) have changed.
                  IF(ELEM_TYPE(1:5).EQ.'FIELD') THEN
C SGM 26Oct2000 grid-based Finite element also
C MLT 29Nov02 Grid Finite Volume also
                    IF(ITYP4(nr,nx).EQ.4.OR.ITYP4(nr,nx).EQ.6.OR.
     '                 ITYP4(nr,nx).EQ.7) THEN !collocation or grid based FE/FV
C *** DPN 30 September 1999 - nqsc_old is always -1, which is wrong!!
C ***   FIELDS_CHANGED should not always be .TRUE.
c                      IF(nqsc_old.EQ.-1) THEN
c                        FIELDS_CHANGED=.TRUE.
c                      ELSE
c                        nqsc=NQS(ne)
c                        IF(nqsc.NE.nqsc_old) THEN
c                          nqsc_old=nqsc
c                          FIELDS_CHANGED=.TRUE.
c                        ENDIF
c                      ENDIF
                      IF(nqsc_old.EQ.-1) THEN
                        FIELDS_CHANGED=.TRUE.
                        nqsc_old=NQS(ne)
                      ELSE
                        nqsc=NQS(ne)
                        IF(nqsc.NE.nqsc_old) THEN
                          nqsc_old=nqsc
                          FIELDS_CHANGED=.TRUE.
C DPN 01 February 2000 - This is bad, only sets FIELDS_CHANGED to TRUE
C if the grid scheme has changed but other things can also change
C i.e. Number of nodes per element
C                        ELSE
C                          FIELDS_CHANGED=.FALSE.
                        ENDIF
                      ENDIF
                    ENDIF
                  ELSEIF(GRID_NUMBERS) THEN

C LKC 6-NOV-2003 Ensure grids have been defined
                    CALL ASSERT(CALL_GRID,'>> Define Grid first',
     '                ERROR,*9999)
                    IF(nqsc_old.EQ.-1) THEN
                      FIELDS_CHANGED=.TRUE.
                      nqsc_old=NQS(ne)
                    ELSE
                      nqsc=NQS(ne)
                      IF(nqsc.NE.nqsc_old) THEN
                        nqsc_old=nqsc
                        FIELDS_CHANGED=.TRUE.
                      ENDIF
                    ENDIF
                  ENDIF
                  IF(FIELDS_CHANGED) THEN
                    NE_FIELD=NE
                    DO NP_ELEMENT=0,ELEMENT_NODE_LIST(0)
                      FIELD_NODE_LIST(NP_ELEMENT)=
     '                  ELEMENT_NODE_LIST(NP_ELEMENT)
                    ENDDO
C**                 write the scale factor information
                    IF(DATAFILE) THEN
                      WRITE(IFILE,'( '' #Scale factor sets='',I2)')
     '                  NBLIST(0)
                    ELSE
                      IF(FSKWRITE(ICHAR('#'),SK_CHAR,1,CONNID2).EQ.-1)
     '                  GOTO 9999
                      ISOCKET(1)=NBLIST(0)
                      IF(FSKWRITE(ISOCKET,SK_LONG_INT,1,CONNID2).EQ.
     '                  -1) GOTO 9999
                    ENDIF
                    DO nb=1,NBFT
                      IF(NBLIST(nb).NE.0) THEN
                        CALL CALCULATE_BASIS_STRING(IBT,nb,DATAFILE,
     '                    BASES,SPECIAL_BASIS_FLAG,ERROR,*9999)
                        IF(SPECIAL_BASIS_FLAG.GE.0) THEN
C CPB 3/8/95 Sector elements are now the same as the other elements
C for the ns calculation
C                          IF(SPECIAL_BASIS_FLAG.EQ.3) THEN
C                            IF(2.EQ.NIT(nb)) THEN
C                              NN_LAY=NNT(nb)
C                            ELSE
C                              IF(IBT(1,3,NB).EQ.2) THEN
C                                NN_LAY=NNT(nb)/2
C                              ELSE
C                                NN_LAY=NNT(nb)/(IBT(2,3,NB)+1)
C                              ENDIF
C                            ENDIF
C                            DO nn=1,NNT(nb)
C                              IF(MOD(nn,NN_LAY).EQ.1) THEN
C                                NS_TOT=NS_TOT+(IBT(2,1,NB)+1)*NKT(nn,NB)
C                              ELSE
C                                NS_TOT=NS_TOT+NKT(nn,NB)
C                              ENDIF
C                            ENDDO
                          IF(SPECIAL_BASIS_FLAG.GE.3) THEN
C KAT 25Jan00:              Cmgui can only have one # scale factors per
C                           basis type, but does not know the difference
C                           between tensor product basis functions and
C                           those that are already collapsed.  Therefore
C                           the full number of scale factors for a
C                           tensor product basis function are always
C                           output.
                            NS_TOT=1
                            DO ni=1,NIT(nb)
                              IF(IBT(1,ni,nb).EQ.1) THEN !Lagrange
                                NS_TOT=NS_TOT*(IBT(2,ni,nb)+1)
                              ELSE IF(IBT(1,ni,nb).LE.4) THEN
                                !Hermite/Simplex/Seredipity
                                IF(IBT(2,ni,nb).EQ.1) THEN !cubic Hermite
                                  NS_TOT=NS_TOT*4
                                ELSE !Quadratic
                                  NS_TOT=NS_TOT*3
                                ENDIF !IBT(2,ni,nb)
                              ELSE IF(IBT(1,ni,nb).LE.8) THEN !Sector
                                IF(IBT(2,ni,nb).EQ.4) THEN !Hermite
                                  NS_TOT=NS_TOT*4
                                ELSE !Lagrange
                                  NS_TOT=NS_TOT*(IBT(2,ni,nb)+1)
                                ENDIF !IBT(2,ni,nb)
                              ELSE !Fourier
                                NS_TOT=NS_TOT*IBT(2,ni,nb)
                              ENDIF !IBT(1,ni,nb)
                            ENDDO !ni
                            IF(NS_TOT.GT.NSM) THEN
                              IEND=0
                              CALL APPENDC(IEND,' >>Increase NSM to ',
     '                          ERROR)
                              CALL APPENDI(IEND,NS_TOT,ERROR)
                              GOTO 9999
                            ENDIF
                          ELSE
                            NS_TOT=0
                            DO nn=1,NNT(nb)
                              NS_TOT=NS_TOT+NKT(nn,nb)
                            ENDDO
                          ENDIF
                          IF(DATAFILE) THEN
                            CALL STRING_TRIM(BASES,IBEG,IEND)
                            WRITE(IFILE,
     '                        '(3X,A,'', #Scale factors='',I2)')
     '                        BASES(IBEG:IEND),NS_TOT
                          ELSE
                            IF(FSKWRITE(NS_TOT,SK_LONG_INT,1,
     '                        CONNID2).EQ.-1) GOTO 9999
                          ENDIF
                          NST_LIST(nb)=NS_TOT
                        ENDIF
                      ENDIF
                    ENDDO !nb
C**                 write the node information
                    IF(DATAFILE) THEN
                      WRITE(IFILE,'( '' #Nodes='',I12)')
     '                  ELEMENT_NODE_LIST(0)
                    ELSE
                      IF(FSKWRITE(ELEMENT_NODE_LIST(0),SK_LONG_INT,1,
     '                  CONNID2).EQ.-1) GOTO 9999
                    ENDIF

C LKC 20-DEC-97 new method for determining num of fields for use
C  in export_field_heading

C*** Creating  a NJLIST
C    ! contains the types of existing fields at a given node
                    DO njj1=1,3 !geometry, fibre, general
                      NJLIST(njj1)=.FALSE.
                    ENDDO

                    DO nrr=1,NRLIST(0) !all reg. exporting
                      DO njj1=1,3 !each field type
                        IF(NJ_LOC(njj1,0,NRLIST(nrr)).NE.0) THEN
                          NJLIST(njj1)=.TRUE.
                        ENDIF
                      ENDDO
                    ENDDO

C*** Write the field information
                    NJH_FIELD_BASE_LIST(0)=0

C ***
C *** DPN - this loop should be fine, just need to add to
C ***       EXPORT_FIELD_HEADING
C ***
C ***     - DON'T want to go into this loop ???
C ***
                    IF(EXPORT_CELL) NJHT=-NJHT
                    DO no_njh=1,NJHT
                      njh=NJH_LIST(no_njh)

                      IF(ELEM_TYPE(1:8).EQ.'GEOMETRY') THEN
C**                     reorder nj's
                        IF((1.LE.no_njh).AND.
     '                    (no_njh.LE.NJ_LOC(NJL_GEOM,0,nr))) THEN
                          NJ=NJ_LOC(NJL_GEOM,no_njh,nr)
                        ELSE IF((NJ_LOC(NJL_GEOM,0,nr).LT.no_njh).AND.
     '                      (no_njh.LE.NJ_LOC(NJL_GEOM,0,nr)+
     '                      NJ_LOC(NJL_FIBR,0,nr))) THEN
                          NJ=
     '                      NJ_LOC(NJL_FIBR,no_njh-
     '                      NJ_LOC(NJL_GEOM,0,nr),nr)
                        ELSE IF((NJ_LOC(NJL_GEOM,0,nr)+
     '                      NJ_LOC(NJL_FIBR,0,nr)
     '                      .LT.no_njh).AND.(no_njh.LE.
     '                      NJ_LOC(NJL_GEOM,0,nr)+
     '                      NJ_LOC(NJL_FIBR,0,nr)+
     '                      NJ_LOC(NJL_FIEL,0,nr))) THEN
                          NJ=NJ_LOC(NJL_FIEL,no_njh-
     '                      NJ_LOC(NJL_FIBR,0,nr)-
     '                      NJ_LOC(NJL_GEOM,0,nr),nr)
                        ELSE
                          ERROR='>>Invalid nj'
                          GOTO 9999
                        ENDIF
C!!! KAT: But this nj has a different meaning!
                        nb=NBJ(nj,NE)
                        CALL EXPORT_FIELD_HEADING(ELEMENT_DIMENSION,
     '                    FIELD_BASE_TYPE,IBT,IFILE,0,iy,nb,NBH,
     '                    nc,0,ne,NFIELDT,0,NHE(1,nx),NHP(1,nr,nx),
     '                    nj,NJLIST,0,nr,1,NW(1,1,nx),nx,
     '                    SPECIAL_BASIS_FLAG,AUTONAME,
     '                    DATAFILE,.FALSE.,FIELD_EX_TYPE,FIELD_NAME,
     '                    SET_FIELD_NAME,GRID_NUMBERS,ERROR,*9999)
                      ELSE !field
C MLB 220803            nb=NBH(njh,1,NE)
                        nb=NBH(NH_LOC(njh,nx),1,NE)
C SGM 26Oct2000 grid-based Finite element also
C MLT 29Nov02 Grid Finite Volume also
                        IF(ITYP4(nr,nx).EQ.4.OR.ITYP4(nr,nx).EQ.6.OR.
     '                     ITYP4(nr,nx).EQ.7) THEN !Collocation or grid based FE/FV
                          nb=NBJ(1,ne)
C KAT 23/3/00: Transferred from EXPORT_FIELD_HEADING_CELL
                          NJH_FIELD_BASE_LIST(0)=
     '                      NJH_FIELD_BASE_LIST(0)+1
                          NJH_FIELD_BASE_LIST(NJH_FIELD_BASE_LIST(0))=
     '                      njh
                        ENDIF
                        CALL EXPORT_FIELD_HEADING(ELEMENT_DIMENSION,
     '                    FIELD_BASE_TYPE,IBT,IFILE,0,iy,nb,NBH,
     '                    nc,0,ne,NFIELDT,njh,NHE(1,nx),NHP(1,nr,nx),0,
     '                    NJLIST,0,nr,1,NW(1,1,nx),nx,
     '                    SPECIAL_BASIS_FLAG,AUTONAME,
     '                    DATAFILE,.FALSE.,FIELD_EX_TYPE,FIELD_NAME,
     '                    SET_FIELD_NAME,GRID_NUMBERS,ERROR,*9999)
                      ENDIF !geometry/field
                      IF(FIELD_BASE_TYPE.EQ.1) THEN
C!!! DB.  Temporary.  Stops auxilliary basis export
                        IF(NBC(nb).EQ.0) THEN
                          SPECIAL_BASIS_FLAG=-1
                        ENDIF
                        IF(SPECIAL_BASIS_FLAG.GE.0) THEN
                          NS_TOT=0
                          DO I=1,nb-1
                            IF(NBLIST(I).NE.0) THEN
C KAT 25Jan00: Exporting scale factors corresponding to a full tensor product.
                              NS_TOT=NS_TOT+NST_LIST(I)
C                              DO nn=1,NNT(I)
C                                NS_TOT=NS_TOT+NKT(nn,I)
C                              ENDDO
                            ENDIF !NBLIST(I).NE.0
                          ENDDO !I
                          IF(SPECIAL_BASIS_FLAG.EQ.0) THEN
C**                       not Hermite simplex
                            IF(DATAFILE) THEN
                              WRITE(IFILE,'(5X,''#Nodes='',I2)') NNT(nb)
                            ELSE
                              IF(FSKWRITE(NNT(nb),SK_LONG_INT,1,
     '                          CONNID2).EQ.-1) GOTO 9999
                            ENDIF
                            DO nn=1,NNT(nb)
                              NP=NPNE(nn,nb,NE)
                              NP_ELEMENT=ELEMENT_NODE_LIST(0)
                              DO WHILE ((NP_ELEMENT.GT.0).AND.
     '                          (NP.NE.ELEMENT_NODE_LIST(NP_ELEMENT)))
                                NP_ELEMENT=NP_ELEMENT-1
                              ENDDO
                              IF(DATAFILE) THEN
                                WRITE(IFILE,'(5X,I2,'
     '                            //'''.  #Values='',I1)')
     '                            NP_ELEMENT,NKT(nn,nb)
                              ELSE
                                IF(FSKWRITE(NP_ELEMENT,SK_LONG_INT,1,
     '                            CONNID2).EQ.-1) GOTO 9999
                                IF(FSKWRITE(NKT(nn,NB),SK_LONG_INT,1,
     '                            CONNID2).EQ.-1) GOTO 9999
                              ENDIF
                              IF(ELEM_TYPE(1:8).EQ.'GEOMETRY') THEN
                                nv_offset=NVJE(nn,nb,njh,ne)
                                nv_offset=(nv_offset-1)*NKJ(njh,NP)
                                DO mk=1,NKT(nn,nb)
                                  NKLIST(mk)=NKJE(mk,nn,njh,ne)
                                ENDDO !mk
                              ELSE
                                nv_offset=NVHE(nn,nb,NH_LOC(njh,nx),ne)
C!!! KAT 9Nov98: Not consistent with node file
                                nv_offset=(nv_offset-1)*NKT(nn,NB)
                                DO mk=1,NKT(nn,nb)
                                  NKLIST(mk)=
     '                              NKHE(mk,nn,NH_LOC(njh,nx),ne)
                                ENDDO !mk
                              ENDIF
C KAT 9Nov98: Not consistent wtih node file
C                            nv_offset=(nv_offset-1)*NKT(nn,NB)
                              IF(DATAFILE) THEN
                                WRITE(IFILE,'(7X,''Value indices:  '','
     '                            //'12(1X,I3))')
     '                            (nv_offset+NKLIST(mk),mk=1,NKT(nn,nb))
                                WRITE(IFILE,
     '                            '(7X,''Scale factor indices:'','
     '                            //'12(1X,I3))')
     '                            (NS_TOT+NK,NK=1,NKT(nn,NB))
                              ELSE
                                DO mk=1,NKT(nn,NB)
                                  IF(FSKWRITE(nv_offset+
     '                              NKLIST(mk),
     '                              SK_LONG_INT,1,CONNID2).EQ.-1)
     '                              GOTO 9999
                                ENDDO !mk
                                IF(FSKWRITE(NS_TOT+1,SK_LONG_INT,
     '                            NKT(nn,NB),CONNID2).EQ.-1) GOTO 9999
                              ENDIF
                              NS_TOT=NS_TOT+NKT(nn,NB)
                            ENDDO !nn
                          ELSE IF(SPECIAL_BASIS_FLAG.EQ.1) THEN
C**                       apex at node 1
                            IF(DATAFILE) THEN
                              WRITE(IFILE,'(5X,''#Nodes=4'' )')
                            ELSE
                              IF(FSKWRITE(4,SK_LONG_INT,1,
     '                          CONNID2).EQ.-1)
     '                          GOTO 9999
                            ENDIF
C**                       node 1
                            NP=NPNE(1,nb,NE)
                            NP_ELEMENT=ELEMENT_NODE_LIST(0)
                            DO WHILE ((NP_ELEMENT.GT.0).AND.
     '                        (NP.NE.ELEMENT_NODE_LIST(NP_ELEMENT)))
                              NP_ELEMENT=NP_ELEMENT-1
                            ENDDO
                            IF(DATAFILE) THEN
                              WRITE(IFILE,'(5X,I2,''.  #Values=2'')')
     '                          NP_ELEMENT
                            ELSE
                              ISOCKET(1)=NP_ELEMENT
                              ISOCKET(2)=2
                            ENDIF
                            IF(ELEM_TYPE(1:8).EQ.'GEOMETRY') THEN
                              nv_offset=NVJE(1,nb,njh,ne)
                              NKLIST(1)=NKJE(1,1,njh,ne)
                            ELSE
                              nv_offset=NVHE(1,nb,NH_LOC(njh,nx),ne)
                              NKLIST(1)=NKHE(1,1,NH_LOC(njh,nx),ne)
                            ENDIF
                            nv_offset=(nv_offset-1)*NKT(1,NB)
                            IF(DATAFILE) THEN
                              WRITE(IFILE,'(7X,''Value indices:  '','
     '                          //'12(1X,I3))') nv_offset+
     '                          NKLIST(1),0
                              WRITE(IFILE,'(7X,'
     '                          //'''Scale factor indices:'','
     '                          //'12(1X,I3))') NS_TOT+1,0
                            ELSE
                              ISOCKET(3)=nv_offset+NKLIST(1)
                              ISOCKET(4)=0
                              ISOCKET(5)=NS_TOT+1
                              ISOCKET(6)=0
                              IF(FSKWRITE(ISOCKET,SK_LONG_INT,6,
     '                          CONNID2).EQ.-1) GOTO 9999
                            ENDIF
C                         node 2 (repeat node 1)
                            IF(DATAFILE) THEN
                              WRITE(IFILE,'(5X,I2,''.  #Values=2'')')
     '                          NP_ELEMENT
                              WRITE(IFILE,'(7X,''Value indices:  '','
     '                          //'12(1X,I3))') nv_offset+
     '                          NKLIST(1),0
                              WRITE(IFILE,'(7X,'
     '                          //'''Scale factor indices:'','
     '                          //'12(1X,I3))') NS_TOT+1,0
                            ELSE
                              IF(FSKWRITE(ISOCKET,SK_LONG_INT,6,
     '                          CONNID2).EQ.-1) GOTO 9999
                            ENDIF
                            NS_TOT=NS_TOT+NKT(1,NB)
C                         node 3
                            NP=NPNE(2,nb,NE)
                            NP_ELEMENT=ELEMENT_NODE_LIST(0)
                            DO WHILE ((NP_ELEMENT.GT.0).AND.
     '                        (NP.NE.ELEMENT_NODE_LIST(NP_ELEMENT)))
                              NP_ELEMENT=NP_ELEMENT-1
                            ENDDO
                            IF(DATAFILE) THEN
                              WRITE(IFILE,'(5X,I2,''.  #Values='',I1)')
     '                          NP_ELEMENT,NKT(2,NB)
                            ELSE
                              IF(FSKWRITE(NP_ELEMENT,SK_LONG_INT,1,
     '                          CONNID2).EQ.-1) GOTO 9999
                              IF(FSKWRITE(NKT(2,NB),SK_LONG_INT,1,
     '                          CONNID2).EQ.-1) GOTO 9999
                            ENDIF
                            IF(ELEM_TYPE(1:8).EQ.'GEOMETRY') THEN
                              nv_offset=NVJE(2,nb,njh,ne)
                              DO mk=1,NKT(2,nb)
                                NKLIST(mk)=NKJE(mk,2,njh,ne)
                              ENDDO !mk
                            ELSE
                              nv_offset=NVHE(2,nb,NH_LOC(njh,nx),ne)
                              DO mk=1,NKT(2,nb)
                                NKLIST(mk)=
     '                            NKHE(mk,2,NH_LOC(njh,nx),ne)
                              ENDDO !mk
                            ENDIF
                            nv_offset=(nv_offset-1)*NKT(2,NB)
                            IF(DATAFILE) THEN
                              WRITE(IFILE,'(7X,''Value indices:  '','
     '                          //'12(1X,I3))')
     '                          (nv_offset+NKLIST(mk),
     '                          MK=1,NKT(2,NB))
                              WRITE(IFILE,'(7X,'
     '                          //'''Scale factor indices:'','
     '                          //'12(1X,I3))')
     '                          (NS_TOT+NK,NK=1,NKT(2,NB))
                            ELSE
                              DO mk=1,NKT(2,NB)
                                IF(FSKWRITE(nv_offset+NKLIST(mk),
     '                            SK_LONG_INT,1,CONNID2).EQ.-1)
     '                            GOTO 9999
                              ENDDO
                              IF(FSKWRITE(NS_TOT+1,SK_LONG_INT,
     '                          NKT(2,NB),
     '                          CONNID2).EQ.-1) GOTO 9999
                            ENDIF
                            NS_TOT=NS_TOT+NKT(2,NB)
C                         node 4
                            NP=NPNE(3,nb,NE)
                            NP_ELEMENT=ELEMENT_NODE_LIST(0)
                            DO WHILE ((NP_ELEMENT.GT.0).AND.
     '                        (NP.NE.ELEMENT_NODE_LIST(NP_ELEMENT)))
                              NP_ELEMENT=NP_ELEMENT-1
                            ENDDO
                            IF(DATAFILE) THEN
                              WRITE(IFILE,'(5X,I2,''.  #Values='',I1)')
     '                          NP_ELEMENT,NKT(3,NB)
                            ELSE
                              IF(FSKWRITE(NP_ELEMENT,SK_LONG_INT,1,
     '                          CONNID2).EQ.-1) GOTO 9999
                              IF(FSKWRITE(NKT(3,NB),SK_LONG_INT,1,
     '                          CONNID2).EQ.-1) GOTO 9999
                            ENDIF
                            IF(ELEM_TYPE(1:8).EQ.'GEOMETRY') THEN
                              nv_offset=NVJE(3,nb,njh,ne)
                              DO mk=1,NKT(3,nb)
                                NKLIST(mk)=NKJE(mk,3,njh,ne)
                              ENDDO !mk
                            ELSE
                              nv_offset=NVHE(3,nb,NH_LOC(njh,nx),ne)
                              DO mk=1,NKT(3,nb)
                                NKLIST(mk)=
     '                            NKHE(mk,3,NH_LOC(njh,nx),ne)
                              ENDDO !mk
                            ENDIF
                            nv_offset=(nv_offset-1)*NKT(3,NB)
                            IF(DATAFILE) THEN
                              WRITE(IFILE,'(7X,''Value indices:  '','
     '                          //'12(1X,I3))')
     '                          (nv_offset+NKLIST(mk),MK=1,
     '                          NKT(3,NB))
                              WRITE(IFILE,'(7X,'
     '                          //'''Scale factor indices:'','
     '                          //'12(1X,I3))')
     '                          (NS_TOT+NK,NK=1,NKT(3,NB))
                            ELSE
                              DO mk=1,NKT(3,NB)
                                IF(FSKWRITE(nv_offset+NKLIST(mk),
     '                            SK_LONG_INT,1,CONNID2).EQ.-1)
     '                            GOTO 9999
                              ENDDO
                              IF(FSKWRITE(NS_TOT+1,SK_LONG_INT,
     '                          NKT(3,NB),
     '                          CONNID2).EQ.-1) GOTO 9999
                            ENDIF
                          ELSE IF(SPECIAL_BASIS_FLAG.EQ.2) THEN
C**                       apex at node 3
                            IF(DATAFILE) THEN
                              WRITE(IFILE,'(5X,''#Nodes=4'' )')
                            ELSE
                              IF(FSKWRITE(4,SK_LONG_INT,1,
     '                          CONNID2).EQ.-1)
     '                          GOTO 9999
                            ENDIF
C                         node 1
                            NP=NPNE(1,nb,NE)
                            NP_ELEMENT=ELEMENT_NODE_LIST(0)
                            DO WHILE ((NP_ELEMENT.GT.0).AND.
     '                        (NP.NE.ELEMENT_NODE_LIST(NP_ELEMENT)))
                              NP_ELEMENT=NP_ELEMENT-1
                            ENDDO
                            IF(DATAFILE) THEN
                              WRITE(IFILE,'(5X,I2,''.  #Values='',I1)')
     '                          NP_ELEMENT,NKT(1,NB)
                            ELSE
                              IF(FSKWRITE(NP_ELEMENT,SK_LONG_INT,1,
     '                          CONNID2).EQ.-1) GOTO 9999
                              IF(FSKWRITE(NKT(1,NB),SK_LONG_INT,1,
     '                          CONNID2).EQ.-1) GOTO 9999
                            ENDIF
                            IF(ELEM_TYPE(1:8).EQ.'GEOMETRY') THEN
                              nv_offset=NVJE(1,nb,njh,ne)
                              DO mk=1,NKT(1,nb)
                                NKLIST(mk)=NKJE(mk,1,njh,ne)
                              ENDDO !mk
                            ELSE
                              nv_offset=NVHE(1,nb,NH_LOC(njh,nx),ne)
                              DO mk=1,NKT(1,nb)
                                NKLIST(mk)=
     '                            NKHE(mk,1,NH_LOC(njh,nx),ne)
                              ENDDO !mk
                            ENDIF
                            nv_offset=(nv_offset-1)*NKT(1,NB)
                            IF(DATAFILE) THEN
                              WRITE(IFILE,'(7X,''Value indices:  '','
     '                          //'12(1X,I3))')
     '                          (nv_offset+NKLIST(mk),mk=1,
     '                          NKT(1,NB))
                              WRITE(IFILE,'(7X,'//
     '                          '''Scale factor indices:'',12(1X,I3))')
     '                          (NS_TOT+NK,NK=1,NKT(1,NB))
                            ELSE
                              DO mk=1,NKT(1,NB)
                                IF(FSKWRITE(nv_offset+NKLIST(mk),
     '                            SK_LONG_INT,1,CONNID2).EQ.-1)
     '                            GOTO 9999
                              ENDDO
                              IF(FSKWRITE(NS_TOT+1,SK_LONG_INT,
     '                          NKT(1,NB),
     '                          CONNID2).EQ.-1) GOTO 9999
                            ENDIF
                            NS_TOT=NS_TOT+NKT(1,NB)
C                         node 2
                            NP=NPNE(2,nb,NE)
                            NP_ELEMENT=ELEMENT_NODE_LIST(0)
                            DO WHILE ((NP_ELEMENT.GT.0).AND.
     '                        (NP.NE.ELEMENT_NODE_LIST(NP_ELEMENT)))
                              NP_ELEMENT=NP_ELEMENT-1
                            ENDDO
                            IF(DATAFILE) THEN
                              WRITE(IFILE,'(5X,I2,''.  #Values='',I1)')
     '                          NP_ELEMENT,NKT(2,NB)
                            ELSE
                              IF(FSKWRITE(NP_ELEMENT,SK_LONG_INT,1,
     '                          CONNID2).EQ.-1) GOTO 9999
                              IF(FSKWRITE(NKT(2,NB),SK_LONG_INT,1,
     '                          CONNID2).EQ.-1) GOTO 9999
                            ENDIF
                            IF(ELEM_TYPE(1:8).EQ.'GEOMETRY') THEN
                              nv_offset=NVJE(2,nb,njh,ne)
                              DO mk=1,NKT(2,nb)
                                NKLIST(mk)=NKJE(mk,2,njh,ne)
                              ENDDO !mk
                            ELSE
                              nv_offset=NVHE(2,nb,NH_LOC(njh,nx),ne)
                              DO mk=1,NKT(2,nb)
                                NKLIST(mk)=
     '                            NKHE(mk,2,NH_LOC(njh,nx),ne)
                              ENDDO !mk
                            ENDIF
                            nv_offset=(nv_offset-1)*NKT(2,NB)
                            IF(DATAFILE) THEN
                              WRITE(IFILE,'(7X,''Value indices:  '','
     '                          //'12(1X,I3))')
     '                          (nv_offset+NKLIST(mk),mk=1,
     '                          NKT(2,NB))
                              WRITE(IFILE,'(7X,'
     '                          //'''Scale factor indices:'','
     '                          //'12(1X,I3))')
     '                          (NS_TOT+NK,NK=1,NKT(2,NB))
                            ELSE
                              DO mk=1,NKT(2,NB)
                                IF(FSKWRITE(nv_offset+NKLIST(mk),
     '                            SK_LONG_INT,1,CONNID2).EQ.-1)
     '                            GOTO 9999
                              ENDDO
                              IF(FSKWRITE(NS_TOT+1,SK_LONG_INT,
     '                          NKT(2,NB),
     '                          CONNID2).EQ.-1) GOTO 9999
                            ENDIF
                            NS_TOT=NS_TOT+NKT(2,NB)
C                         node 3
                            NP=NPNE(3,nb,NE)
                            NP_ELEMENT=ELEMENT_NODE_LIST(0)
                            DO WHILE ((NP_ELEMENT.GT.0).AND.
     '                        (NP.NE.ELEMENT_NODE_LIST(NP_ELEMENT)))
                              NP_ELEMENT=NP_ELEMENT-1
                            ENDDO
                            IF(DATAFILE) THEN
                              WRITE(IFILE,'(5X,I2,''.  #Values=2'')')
     '                          NP_ELEMENT
                            ELSE
                              ISOCKET(1)=NP_ELEMENT
                              ISOCKET(2)=2
                            ENDIF
                            IF(ELEM_TYPE(1:8).EQ.'GEOMETRY') THEN
                              nv_offset=NVJE(3,nb,njh,ne)
                              NKLIST(1)=NKJE(1,3,njh,ne)
                            ELSE
                              nv_offset=NVHE(3,nb,NH_LOC(njh,nx),ne)
                              NKLIST(1)=NKHE(1,3,NH_LOC(njh,nx),ne)
                            ENDIF
                            nv_offset=(nv_offset-1)*NKT(3,NB)
                            IF(DATAFILE) THEN
                              WRITE(IFILE,'(7X,''Value indices:  '','
     '                          //'12(1X,I3))') nv_offset+
     '                          NKLIST(1),0
                              WRITE(IFILE,'(7X,'
     '                          //'''Scale factor indices:'','
     '                          //'12(1X,I3))') NS_TOT+1,0
                            ELSE
                              ISOCKET(3)=nv_offset+NKLIST(1)
                              ISOCKET(4)=0
                              ISOCKET(5)=NS_TOT+1
                              ISOCKET(6)=0
                              IF(FSKWRITE(ISOCKET,SK_LONG_INT,6,
     '                          CONNID2).EQ.-1) GOTO 9999
                            ENDIF
C                         node 4 (repeat of node 3)
                            IF(DATAFILE) THEN
                              WRITE(IFILE,'(5X,I2,''.  #Values=2'' )')
     '                          NP_ELEMENT
                              WRITE(IFILE,'(7X,''Value indices:  '','
     '                          //'12(1X,I3))') nv_offset+
     '                          NKLIST(1),0
                              WRITE(IFILE,'(7X,'
     '                          //'''Scale factor indices:'','
     '                          //'12(1X,I3))') NS_TOT+1,0
                            ELSE
                              IF(FSKWRITE(ISOCKET,SK_LONG_INT,6,
     '                          CONNID2).EQ.-1) GOTO 9999
                            ENDIF
                          ELSE IF(SPECIAL_BASIS_FLAG.GE.3) THEN
C**                       sector element

C 25/2/97 LC archived section : extending sectors

                            NN_TOT=1
                            DO ni=1,NIT(nb)
                              IF(IBT(1,ni,nb).EQ.1) THEN
                                NUM_NN=IBT(2,ni,nb)+1
                              ELSE IF(IBT(1,ni,nb).EQ.2) THEN
                                NUM_NN=2
                              ELSE IF(IBT(1,ni,nb).EQ.5.OR.
     '                            IBT(1,ni,nb).EQ.6) THEN
                                IF(IBT(2,ni,nb).EQ.4) THEN
                                  NUM_NN=2
                                ELSE
                                  NUM_NN=IBT(2,ni,nb)+1
                                ENDIF
                              ENDIF
                              NUMNODES(ni)=NUM_NN
                              NN_TOT=NN_TOT*NUM_NN
                              POSITION(ni)=1
                            ENDDO !ni
                            IF(DATAFILE) THEN
                              WRITE(IFILE,'(5X,''#Nodes='',I2)') NN_TOT
                            ELSE
                              IF(FSKWRITE(NN_TOT,SK_LONG_INT,1,
     '                          CONNID2).EQ.-1) GOTO 9999
                            ENDIF
                            NP_ELEMENT=0
                            DO i=1,NN_TOT
                              NUMCOLLAPSE=0
                              ATCOLLAPSE=.FALSE.
                              DO ni=1,NIT(nb)
                                IF(IBT(1,ni,nb).EQ.5) THEN
                                  IF(POSITION(IBT(3,ni,nb)).EQ.1) THEN
                                    ATCOLLAPSE=.TRUE.
                                    NUMCOLLAPSE=NUMCOLLAPSE+1
                                    nicollapse=ni
                                  ENDIF
                                ELSE IF(IBT(1,ni,nb).EQ.6) THEN
                                  IF(POSITION(IBT(3,ni,nb)).EQ.
     '                              NUMNODES(IBT(3,ni,nb))) THEN
                                    ATCOLLAPSE=.TRUE.
                                    NUMCOLLAPSE=NUMCOLLAPSE+1
                                    nicollapse=ni
                                  ENDIF
                                ENDIF
                              ENDDO !ni
                              IF(ATCOLLAPSE) THEN
                                IF(NUMCOLLAPSE.EQ.2) THEN
                                  nn=1
                                  FOUND=.FALSE.
                                  ni=IBT(3,nicollapse,nb)
                                  DO WHILE(.NOT.FOUND.AND.nn.LE.NNT(nb))
                                    IF(IBT(1,nicollapse,nb).EQ.5) THEN
                                      FOUND=POSITION(ni).EQ.1
                                    ELSE
                                      FOUND=POSITION(ni).EQ.NUMNODES(ni)
                                    ENDIF
                                    IF(.NOT.FOUND) nn=nn+1
                                  ENDDO
                                ELSE
                                  nn=1
                                  FOUND=.FALSE.
                                  ni2=IBT(3,nicollapse,nb)
                                  DO WHILE(.NOT.FOUND.AND.nn.LE.NNT(nb))
                                    COUNT=0
                                    DO ni=1,NIT(nb)
                                      IF(ni.NE.nicollapse) THEN
                                        IF(IBT(1,nicollapse,nb).EQ.5)
     '                                    THEN
                                          IF(INP(nn,ni,nb).EQ.
     '                                      POSITION(ni).AND.
     '                                      INP(nn,ni2,nb).EQ.1)
     '                                      COUNT=COUNT+1
                                        ELSE
                                          IF(INP(nn,ni,nb).EQ.
     '                                      POSITION(ni).AND.
     '                                      INP(nn,ni2,nb).EQ.
     '                                      NUMNODES(ni)) COUNT=COUNT+1
                                        ENDIF
                                      ENDIF
                                    ENDDO
                                    IF(COUNT.EQ.(NIT(nb)-1)) THEN
                                      FOUND=.TRUE.
                                    ELSE
                                      nn=nn+1
                                    ENDIF
                                  ENDDO
                                ENDIF
                              ELSE
                                nn=1
                                FOUND=.FALSE.
                                DO WHILE(.NOT.FOUND.AND.nn.LE.NNT(nb))
                                  COUNT=0
                                  DO ni=1,NIT(nb)
                                    IF(INP(nn,ni,nb).EQ.POSITION(ni))
     '                                COUNT=COUNT+1
                                  ENDDO
                                  IF(COUNT.EQ.NIT(nb)) THEN
                                    FOUND=.TRUE.
                                  ELSE
                                    nn=nn+1
                                  ENDIF
                                ENDDO
                              ENDIF
                              CALL ASSERT(FOUND,
     '                          '>>Could not find local node',
     '                          ERROR,*9999)
                              NP_ELEMENT=nn
                              NK_TOT=1
                              DO ni=1,NIT(nb)
                                IF((IBT(1,ni,nb).EQ.2.AND.
     '                            (IBT(2,ni,nb).EQ.1.OR.
     '                            POSITION(ni).NE.(IBT(2,ni,nb)-1))).OR.
     '                            ((IBT(1,ni,nb).EQ.5.OR.
     '                            IBT(1,ni,nb).EQ.6).AND.
     '                            IBT(2,ni,nb).EQ.4)) THEN
                                  DO nk=1,NK_TOT
                                    IDOTEMP(nk,ni)=1
                                    DO ni2=1,ni-1
                                      IDOTEMP(nk+NK_TOT,ni2)=
     '                                  IDOTEMP(nk,ni2)
                                    ENDDO !ni2
                                  ENDDO !nk
                                  DO nk=NK_TOT+1,NK_TOT*2
                                    IDOTEMP(nk,ni)=2
                                  ENDDO !nk
                                  NK_TOT=NK_TOT*2
                                ELSE
                                  DO nk=1,NK_TOT
                                    IDOTEMP(nk,ni)=1
                                  ENDDO !nk
                                ENDIF
                              ENDDO !ni
                              IF(DATAFILE) THEN
                                WRITE(IFILE,'(5X,I2,''.  #Values='','
     '                            //'I1)') NP_ELEMENT,NK_TOT
                              ELSE
                                IF(FSKWRITE(NP_ELEMENT,SK_LONG_INT,1,
     '                            CONNID2).EQ.-1) GOTO 9999
                                IF(FSKWRITE(NK_TOT,SK_LONG_INT,1,
     '                            CONNID2).EQ.-1) GOTO 9999
                              ENDIF

C LKC 15-SEP-1999 This doesn't work if exporting bicubic geometry
C  which has had a bilinear field fitted to it. (see example e84)
C  We do not want to reset NS_TOT as the Scale factor indices are wrong.
C  This is probably the case for fibres as well (not sure). Onlyt reset
C  NS_TOT for the geometry scale factor indices.
C                              NS_TOT=0

C KAT 25Jan00: Exporting scale factors corresponding to a full tensor product.
CC If the basis function changes between the previous geometry/fibre/field
CC  then don't reset the Scale factor number NS_TOT
C                              IF(((
C     '                          (no_njh.EQ.NJ_LOC(NJL_FIBR,1,nr)).OR.
C     '                          (no_njh.EQ.NJ_LOC(NJL_FIEL,1,nr))).AND.
C     '                          i.EQ.1  ! first node of element
C     '                          )) THEN ! If first component of a field type
C                                NS_TOT_PREV=NS_TOT
C                              ENDIF
C                              NS_TOT=0
C                              DO nn=1,NP_ELEMENT-1
C                                NS_TOT=NS_TOT+NKT(nn,nb)
C                              ENDDO
CC If sectors and exporting fields then change the offset
CC KAT 25Jan00: Already know SPECIAL_BASIS_FLAG is 3
C                              IF(no_njh.GT.NJ_LOC(NJL_GEOM,0,nr)) THEN
CC                              IF(SPECIAL_BASIS_FLAG.EQ.3.AND.
CC     '                          no_njh.GT.NJ_LOC(NJL_GEOM,0,nr)) THEN
C                                NS_TOT=NS_TOT+NS_TOT_PREV
C                              ENDIF

                              IF(ELEM_TYPE(1:8).EQ.'GEOMETRY') THEN
                                nv_offset=NVJE(NP_ELEMENT,nb,njh,ne)
                                DO mk=1,NKT(NP_ELEMENT,nb)
                                  NKLIST1(mk)=NKJE(mk,NP_ELEMENT,njh,ne)
                                ENDDO !mk
                              ELSE
                                nv_offset=
     '                            NVHE(NP_ELEMENT,nb,NH_LOC(njh,nx),ne)
                                DO mk=1,NKT(NP_ELEMENT,nb)
                                  NKLIST1(mk)=NKHE(mk,NP_ELEMENT,
     '                              NH_LOC(njh,nx),ne)
                                ENDDO !mk
                              ENDIF
                              nv_offset=(nv_offset-1)*NKT(NP_ELEMENT,NB)
                              IF(NK_TOT.EQ.NKT(NP_ELEMENT,nb)) THEN
C                               # derivatives at the node is the same
C                               as the maximum # so use default case.
                                DO nk=1,NK_TOT
                                  NKLIST(nk)=nv_offset+NKLIST1(nk)
                                  SFLIST(nk)=NS_TOT+nk
                                ENDDO
                              ELSE
C                               Must find which are the real derivatives
C                               and which derivatives will have zero
C                               for the scale and value indicies
                                DO nk=1,NK_TOT
                                  nkk=1
                                  FOUND=.FALSE.
                                  DO WHILE(.NOT.FOUND.AND.nkk.LE.
     '                              NKT(NP_ELEMENT,nb))
                                    COUNT=0
                                    DO ni=1,NIT(nb)
                                      IF(IDO(nkk,NP_ELEMENT,ni,nb).EQ.
     '                                  IDOTEMP(nk,ni)) COUNT=COUNT+1
                                    ENDDO !ni
                                    IF(COUNT.EQ.NIT(nb)) THEN
                                      FOUND=.TRUE.
                                    ELSE
                                      nkk=nkk+1
                                    ENDIF
                                  ENDDO
                                  IF(FOUND) THEN
                                    NKLIST(nk)=nv_offset+NKLIST1(nkk)
                                    SFLIST(nk)=NS_TOT+nkk
                                  ELSE
                                    NKLIST(nk)=0
                                    SFLIST(nk)=0
                                  ENDIF
                                ENDDO !nk
                              ENDIF
                              IF(DATAFILE) THEN
                                WRITE(IFILE,'(7X,''Value indices:  '','
     '                            //'12(1X,I3))') (NKLIST(nk),nk=1,
     '                            NK_TOT)
                                WRITE(IFILE,'(7X,''Scale factor '
     '                            //'indices:'',12(1X,I3))')
     '                            (SFLIST(nk),nk=1,NK_TOT)
                              ELSE
                                DO nk=1,NK_TOT
                                  IF(FSKWRITE(NKLIST(nk),SK_LONG_INT,1,
     '                              CONNID2).EQ.-1) GOTO 9999
                                ENDDO
                                DO nk=1,NK_TOT
                                  IF(FSKWRITE(SFLIST(nk),SK_LONG_INT,1,
     '                              CONNID2).EQ.-1) GOTO 9999
                                ENDDO
                              ENDIF
                              POSITION(1)=POSITION(1)+1
                              DO ni=1,NIT(nb)
                                IF(POSITION(ni).GT.NUMNODES(ni)) THEN
                                  POSITION(ni)=1
                                  POSITION(ni+1)=POSITION(ni+1)+1
                                ENDIF
                              ENDDO !ni
                              NS_TOT=NS_TOT+NK_TOT
                            ENDDO !i
                          ENDIF
                        ENDIF
                      ELSE IF(FIELD_BASE_TYPE.EQ.2) THEN !Grid Based
                        IF(DATAFILE) THEN
                          nqsc=NQS(ne)
                          IF(NIT(nb).EQ.1) THEN
                            WRITE(IFILE,'(3X,''#xi1='',I5)')
     '                        NQXI(1,nqsc)-1
                          ELSE IF(NIT(nb).EQ.2) THEN
                            WRITE(IFILE,'(3X,''#xi1='',I5,'
     '                        //''', #xi2='',I5)') NQXI(1,nqsc)-1,
     '                        NQXI(2,nqsc)-1
                          ELSE IF(NIT(nb).EQ.3) THEN
                            WRITE(IFILE,'(3X,''#xi1='',I5,'
     '                        //''', #xi2='',I5,'', #xi3='',I5)')
     '                        NQXI(1,nqsc)-1,NQXI(2,nqsc)-1,
     '                        NQXI(3,nqsc)-1
                          ENDIF
                        ELSE
                          ERROR='>>Not implemented'
                          GOTO 9999
                        ENDIF
                      ELSE
                        ERROR='>>Invalid field base type'
                        GOTO 9999
                      ENDIF
                    ENDDO !njh

                    IF(EXPORT_CELL) NJHT=-NJHT

                    IF(GRID_NUMBERS) THEN
C DB 30-JUL-99 Exports grid point number with geometry
                      IF((ELEM_TYPE(1:8).EQ.'GEOMETRY').AND.(NQT.GT.0))
     '                  THEN
C *** DPN 20 January 2000 - can now use integer fields
c                        WRITE(IFILE,
c     '                    '( I5,'') grid_point_number, coordinate,'//
c     '                    ' rectangular cartesian, #Components=1'' )')
c     '                    NFIELDT
                        WRITE(IFILE,
     '                    '( I5,'') grid_point_number, field,'//
     '                    ' integer, #Components=1'' )')
     '                    NFIELDT
C MLB 27/8/99 only works in 3d!
C                      WRITE(IFILE,'( '' grid_point_number.'
C     '                  //' l.Lagrange*l.Lagrange*l.Lagrange,'
C     '                  //' no modify, grid based.'' )')
                        nqsc=NQS(ne)
C MLT 27/11/02 replace use of NIT(nb) to determine number of xi
C coordinates with NQXI(0,nqsc)
                        IF(NQXI(0,nqsc).EQ.1) THEN
                          WRITE(IFILE,'( '' grid_point_number.'
     '                      //' l.Lagrange,'
     '                      //' no modify, grid based.'' )')
                          WRITE(IFILE,'(3X,''#xi1='',I5)')
     '                      NQXI(1,nqsc)-1
                        ELSE IF(NQXI(0,nqsc).EQ.2) THEN
                          WRITE(IFILE,'( '' grid_point_number.'
     '                      //' l.Lagrange*l.Lagrange,'
     '                      //' no modify, grid based.'' )')
                          WRITE(IFILE,'(3X,''#xi1='',I5,'
     '                      //''', #xi2='',I5)') NQXI(1,nqsc)-1,
     '                      NQXI(2,nqsc)-1
                        ELSE IF(NQXI(0,nqsc).EQ.3) THEN
                          WRITE(IFILE,'( '' grid_point_number.'
     '                      //' l.Lagrange*l.Lagrange*l.Lagrange,'
     '                      //' no modify, grid based.'' )')
                          WRITE(IFILE,'(3X,''#xi1='',I5,'
     '                      //''', #xi2='',I5,'', #xi3='',I5)')
     '                      NQXI(1,nqsc)-1,NQXI(2,nqsc)-1,
     '                      NQXI(3,nqsc)-1
                        ENDIF
                      ENDIF
                    ENDIF

                    IF(ELEM_TYPE(1:5).EQ.'FIELD'.AND.
     '                (ITYP2(nr,nx).EQ.2.OR.
     '                ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3)) THEN
C                     For Finite elasticity or fluid mechanics/
C                     const vol constraint, write out the deformed
C                     fibre fields.
                      DO njj=1,NJ_LOC(NJL_FIBR,0,nr)
                        nj=NJ_LOC(NJL_FIBR,njj,nr)
                        nb=NBJ(nj,NE)
                        njh=njj
                        CALL EXPORT_FIELD_HEADING(ELEMENT_DIMENSION,
     '                    FIELD_BASE_TYPE,IBT,IFILE,0,iy,nb,NBH,
     '                    nc,0,ne,NFIELDT,njh,NHE(1,nx),NHP(1,nr,nx),0,
     '                    NJLIST,0,nr,1,NW(1,1,nx),nx,
     '                    SPECIAL_BASIS_FLAG,AUTONAME,
     '                    DATAFILE,.TRUE.,FIELD_EX_TYPE,FIELD_NAME,
     '                    .FALSE.,GRID_NUMBERS,ERROR,*9999)
                        NS_TOT=0
                        DO I=1,NB-1
                          IF(NBLIST(I).NE.0) THEN
                            DO nn=1,NNT(I)
                              NS_TOT=NS_TOT+NKT(nn,I)
                            ENDDO
                          ENDIF
                        ENDDO !I
                        IF(DATAFILE) THEN
                          WRITE(IFILE,'(5X,''#Nodes='',I2)') NNT(nb)
                        ELSE
                          IF(FSKWRITE(NNT(nb),SK_LONG_INT,1,
     '                      CONNID2).EQ.-1) GOTO 9999
                        ENDIF
                        DO nn=1,NNT(nb)
                          NP=NPNE(nn,nb,NE)
                          NP_ELEMENT=ELEMENT_NODE_LIST(0)
                          DO WHILE ((NP_ELEMENT.GT.0).AND.
     '                      (NP.NE.ELEMENT_NODE_LIST(NP_ELEMENT)))
                            NP_ELEMENT=NP_ELEMENT-1
                          ENDDO
                          IF(DATAFILE) THEN
                            WRITE(IFILE,'(5X,I2,''.  #Values='',I1)')
     '                        NP_ELEMENT,NKT(nn,NB)
                          ELSE
                            IF(FSKWRITE(NP_ELEMENT,SK_LONG_INT,1,
     '                        CONNID2).EQ.-1) GOTO 9999
                            IF(FSKWRITE(NKT(nn,NB),SK_LONG_INT,1,
     '                        CONNID2).EQ.-1) GOTO 9999
                          ENDIF
                          nv_offset=(NVJE(nn,nb,nj,ne)-1)*NKT(nn,nb)
                          IF(DATAFILE) THEN
                            WRITE(IFILE,'(7X,''Value indices:  '','
     '                        //'12(1X,I3))')
     '                        (nv_offset+
     '                        NKJE(MK,nn,nj,NE),MK=1,NKT(nn,NB))
                            WRITE(IFILE,
     '                        '(7X,''Scale factor indices:'','
     '                        //'12(1X,I3))')
     '                        (NS_TOT+NK,NK=1,NKT(nn,NB))
                          ELSE
                            DO mk=1,NKT(nn,NB)
                              IF(FSKWRITE(nv_offset+NKJE(mk,nn,nj,NE),
     '                          SK_LONG_INT,1,CONNID2).EQ.-1)
     '                          GOTO 9999
                            ENDDO !mk
                            IF(FSKWRITE(NS_TOT+1,SK_LONG_INT,
     '                        NKT(nn,NB),CONNID2).EQ.-1) GOTO 9999
                          ENDIF
                          NS_TOT=NS_TOT+NKT(nn,NB)
                        ENDDO !nn
                      ENDDO !njj
                    ELSEIF(EXPORT_CELL) THEN
                      nb=NBJ(1,ne)
                      CALL EXPORT_FIELD_HEADING_CELL(IFILE,nb,ne,
     '                  NFIELDT,NQS,NQXI,
     '                  CELL_ICQS_VALUE,CELL_ICQS_SPATIAL,
     '                  CELL_ICQS_NAMES,CELL_RCQS_VALUE,
     '                  CELL_RCQS_SPATIAL,CELL_RCQS_NAMES,
     '                  CELL_YQS_VALUE,CELL_YQS_SPATIAL,CELL_YQS_NAMES,
     '                  ICQS_SPATIAL,IICQS_SPATIAL,IRCQS_SPATIAL,
     '                  RCQS_SPATIAL,ERROR,*9999)
C KAT 23/3/00: Transferred from EXPORT_FIELD_HEADING_CELL
C *** NJH_FIELD_BASE_LIST(0) is the number of fields
                      NJH_FIELD_BASE_LIST(0)=NFIELDT
C *** Loop through the fields (is this an appropriate loop var name?)
                      DO no_njh=1,NFIELDT
                        !used to loop through the variables to be exported
                        NJH_FIELD_BASE_LIST(no_njh)=no_njh
                      ENDDO !j
                    ENDIF !FIELD and finite elasticity, EXPORT_CELL
                  ENDIF !FIELDS_CHANGED
                  IF(ELEM_TYPE(1:8).EQ.'GEOMETRY') THEN
                    nb=NBJ(1,NE)
                  ELSE
C SGM 26Oct2000 grid-based Finite element also
C MLT 29Nov02 Grid Finite Volume also
                    IF(ITYP4(nr,nx).EQ.4.OR.ITYP4(nr,nx).EQ.6.OR.
     '                 ITYP4(nr,nx).EQ.7) THEN
                      nb=NBJ(1,ne) !collocation or grid based FE/FV
                    ELSE
                      nb=NBH(NH_LOC(1,nx),1,NE)
                    ENDIF
                  ENDIF
C**               write the element
                  IF(DATAFILE) THEN
                    WRITE(IFILE,'(1X,''Element: '',I12,'' 0 0'' )')
     '                NE+offset_elem
C**                 write the faces
                    IF(NFE(nb).GT.0) WRITE(IFILE,'(3X,''Faces: '' )')
                  ELSE
                    IF(FSKWRITE(ICHAR('E'),SK_CHAR,1,CONNID2).EQ.-1)
     '                GOTO 9999
                    ISOCKET(1)=1
                    ISOCKET(2)=NE
                    IF(FSKWRITE(ISOCKET,SK_LONG_INT,2,CONNID2).EQ.-1)
     '                GOTO 9999
                  ENDIF
                  CALL CALCULATE_BASIS_STRING(IBT,nb,.TRUE.,BASES,
     '              SPECIAL_BASIS_FLAG,ERROR,*9999)
                  IF(SPECIAL_BASIS_FLAG.EQ.0.OR.
     '              SPECIAL_BASIS_FLAG.EQ.4) THEN
C**                 not Hermite simplex
                    IF(NIT(nb).GE.3) THEN
                      IF(DATAFILE) THEN
                        DO nf_e=1,NFE(nb)
C MPN 26Feb98: Don't add offset for collapsed faces (ie. NFF=0)
                          IF(NFF(nf_e,NE).NE.0) THEN
                            WRITE(IFILE,'(5X,''0 '',I5,'' 0'' )')
     '                        NFF(nf_e,NE)+offset_face
                          ELSE
                            WRITE(IFILE,'(5X,''0 '',I5,'' 0'' )')
     '                        NFF(nf_e,NE)
                          ENDIF
C old                          WRITE(IFILE,'(5X,''0 '',I5,'' 0'' )')
C old     '                      NFF(nf_e,NE)+offset_face
                        ENDDO !nf_e
                      ELSE
                        ISOCKET(1)=2
                        DO nf_e=1,NFE(nb)
                          ISOCKET(2)=NFF(nf_e,NE)
                          IF(FSKWRITE(ISOCKET,SK_LONG_INT,2,
     '                      CONNID2).EQ.-1) GOTO 9999
                        ENDDO !nf_e
                      ENDIF
                    ELSE IF(NIT(nb).GE.2) THEN
C?????????????????????DB.  Numbering order for faces and lines
C                       inconsistent
C                      DO NAE=1,NLE(nb)
C                        WRITE(IFILE,'(5X,''0 0 '',I5)')
C     '                    NLL(NAE,NE)+offset_line
C                      ENDDO

C New simplex elements 12/11/97
C                      CALL ASSERT(NLE(nb).EQ.4,
C     '                  ' Incorrect #lines for element',ERROR,*9999)
                      IF(DATAFILE) THEN


C KAT 2Dec98: no offset for no line
                        LARGE_FORMAT=.FALSE.
                        DO nae=1,NLE(nb)
                          nl=NLL(nae,ne)
                          IF(nl.NE.0) THEN
                            NLFF(nae)=nl+offset_line
                            IF(NLFF(nae).GT.99999) LARGE_FORMAT=.TRUE.
                          ELSE
                            NLFF(nae)=0
                          ENDIF
                        ENDDO !nae

C LKC 26-AUG-2002 Need to check for more then 100000 faces
                        IF(LARGE_FORMAT) THEN
                          FORMAT='( ''   0 0 '',I9)'
                        ELSE
                          FORMAT='( ''   0 0 '',I5)'
                        ENDIF

                        IF(NLE(nb).EQ.3) THEN
                          WRITE(IFILE,FORMAT) NLFF(2)
                          WRITE(IFILE,FORMAT) NLFF(1)
                          WRITE(IFILE,FORMAT) NLFF(3)
                        ELSE IF(NLE(nb).EQ.4) THEN
                          WRITE(IFILE,FORMAT) NLFF(3)
                          WRITE(IFILE,FORMAT) NLFF(4)
                          WRITE(IFILE,FORMAT) NLFF(1)
                          WRITE(IFILE,FORMAT) NLFF(2)
                        ENDIF
                      ELSE
                        IF(NLE(nb).EQ.3) THEN
                          ISOCKET(1)=3
                          ISOCKET(2)=NLL(2,NE)
                          ISOCKET(3)=3
                          ISOCKET(4)=NLL(1,NE)
                          ISOCKET(5)=3
                          ISOCKET(6)=NLL(3,NE)
                          IF(FSKWRITE(ISOCKET,SK_LONG_INT,6,
     '                      CONNID2).EQ.-1) GOTO 9999
                        ELSE IF(NLE(nb).EQ.4) THEN
                          ISOCKET(1)=3
                          ISOCKET(2)=NLL(3,NE)
                          ISOCKET(3)=3
                          ISOCKET(4)=NLL(4,NE)
                          ISOCKET(5)=3
                          ISOCKET(6)=NLL(1,NE)
                          ISOCKET(7)=3
                          ISOCKET(8)=NLL(2,NE)
                          IF(FSKWRITE(ISOCKET,SK_LONG_INT,8,
     '                      CONNID2).EQ.-1) GOTO 9999
                        ENDIF
                      ENDIF
                    ENDIF
                  ELSE IF(SPECIAL_BASIS_FLAG.EQ.1) THEN
C**                 apex at node 1
                    CALL ASSERT(NLE(nb).EQ.3,
     '                ' Incorrect #lines for element',ERROR,*9999)
                    IF(DATAFILE) THEN
                      WRITE(IFILE,'( ''   0 0 '',I5)') NLL(1,NE)+
     '                  offset_line
                      WRITE(IFILE,'( ''   0 0 '',I5)') NLL(2,NE)+
     '                  offset_line
                      WRITE(IFILE,'( ''   0 0 0'')')
                      WRITE(IFILE,'( ''   0 0 '',I5)') NLL(3,NE)+
     '                  offset_line
                    ELSE
                      ISOCKET(1)=3
                      ISOCKET(2)=NLL(1,NE)
                      ISOCKET(3)=3
                      ISOCKET(4)=NLL(2,NE)
                      ISOCKET(5)=3
                      ISOCKET(6)=0
                      ISOCKET(7)=3
                      ISOCKET(8)=NLL(3,NE)
                      IF(FSKWRITE(ISOCKET,SK_LONG_INT,8,
     '                  CONNID2).EQ.-1) GOTO 9999
                    ENDIF
                  ELSE IF(SPECIAL_BASIS_FLAG.EQ.2) THEN
C**                 apex at node 3
                    CALL ASSERT(NLE(nb).EQ.3,
     '                ' Incorrect #lines for element',ERROR,*9999)
                    IF(DATAFILE) THEN
                      WRITE(IFILE,'( ''   0 0 '',I5)') NLL(2,NE)+
     '                  offset_line
                      WRITE(IFILE,'( ''   0 0 '',I5)') NLL(3,NE)+
     '                  offset_line
                      WRITE(IFILE,'( ''   0 0 '',I5)') NLL(1,NE)+
     '                  offset_line
                      WRITE(IFILE,'( ''   0 0 0'')')
                    ELSE
                      ISOCKET(1)=3
                      ISOCKET(2)=NLL(2,NE)
                      ISOCKET(3)=3
                      ISOCKET(4)=NLL(3,NE)
                      ISOCKET(5)=3
                      ISOCKET(6)=NLL(1,NE)
                      ISOCKET(7)=3
                      ISOCKET(8)=0
                      IF(FSKWRITE(ISOCKET,SK_LONG_INT,8,
     '                  CONNID2).EQ.-1) GOTO 9999
                    ENDIF
                  ELSE IF(SPECIAL_BASIS_FLAG.EQ.3) THEN
C**                 sector
                    IF(NIT(nb).GE.3) THEN
c cpb 4/8/95 Generalising sector export
C                      CALL ASSERT(NFE(nb).EQ.5,
C     '                  ' Incorrect #faces for element',ERROR,*9999)
C                      IF(DATAFILE) THEN
C                        WRITE(IFILE,'(5X,''0 '',I5,'' 0'' )')
C     '                    NFF(1,NE)+offset_face
C                        WRITE(IFILE,'(5X,''0 '',I5,'' 0'' )')
C     '                    NFF(2,NE)+offset_face
C                        WRITE(IFILE,'( ''   0 0 0'')')
C                        WRITE(IFILE,'(5X,''0 '',I5,'' 0'' )')
C     '                    NFF(3,NE)+offset_face
C                        WRITE(IFILE,'(5X,''0 '',I5,'' 0'' )')
C     '                    NFF(4,NE)+offset_face
C                        WRITE(IFILE,'(5X,''0 '',I5,'' 0'' )')
C     '                    NFF(5,NE)+offset_face
C                      ELSE
C                        ISOCKET(1)=2
C                        ISOCKET(2)=NFF(1,NE)
C                        ISOCKET(3)=2
C                        ISOCKET(4)=NFF(2,NE)
C                        ISOCKET(5)=2
C                        ISOCKET(6)=NFF(3,NE)
C                        ISOCKET(7)=2
C                        ISOCKET(8)=NFF(4,NE)
C                        ISOCKET(9)=2
C                        ISOCKET(10)=NFF(5,NE)
C                        IF(FSKWRITE(ISOCKET,SK_LONG_INT,10,
C     '                      CONNID2).EQ.-1) GOTO 9999
C                      ENDIF
                      DO ni=1,3
                        IF(IBT(1,ni,nb).EQ.5.OR.IBT(1,ni,nb).EQ.6) THEN
                          ni2=IBT(3,ni,nb)
                          IF(IBT(1,ni,nb).EQ.5) THEN
                            COLLAPSEDNODE=1
                          ELSE
                            IF(IBT(1,ni2,nb).EQ.1) THEN
                              COLLAPSEDNODE=IBT(2,ni2,nb)+1
                            ELSE IF(IBT(1,ni2,nb).EQ.2) THEN
                              COLLAPSEDNODE=2
                            ELSE IF(IBT(1,ni2,nb).EQ.5.OR.
     '                        IBT(1,ni2,nb).EQ.6) THEN
                              IF(IBT(2,ni2,nb).EQ.4) THEN
                                COLLAPSEDNODE=2
                              ELSE
                                COLLAPSEDNODE=IBT(2,ni2,nb)+1
                              ENDIF
                            ENDIF
                          ENDIF
                          COLLAPSED=.TRUE.
                        ELSE
                          COLLAPSED=.FALSE.
                        ENDIF
                        DO i=1,2
                          nf=1
                          FOUND=.FALSE.
                          DO WHILE(.NOT.FOUND.AND.nf.LE.NFE(nb))
                            IF(NNF(1,nf,nb).EQ.ni) THEN
                              COUNT=0
                              DO nnn=1,NNF(0,nf,nb)
                                nn=NNF(1+nnn,nf,nb)
                                IF(COLLAPSED) THEN
                                  IF(I.EQ.1.AND.INP(nn,ni,nb).EQ.1.OR.
     '                              I.EQ.2.AND.(INP(nn,ni,nb).NE.1.OR.
     '                              INP(nn,ni2,nb).EQ.COLLAPSEDNODE))
     '                              COUNT=COUNT+1
                                ELSE
                                  IF(I.EQ.1.AND.INP(nn,ni,nb).EQ.1.OR.
     '                              I.EQ.2.AND.INP(nn,ni,nb).NE.1)
     '                              COUNT=COUNT+1
                                ENDIF
                              ENDDO
                              IF(COUNT.EQ.NNF(0,nf,nb)) FOUND=.TRUE.
                            ENDIF
                            IF(.NOT.FOUND) nf=nf+1
                          ENDDO
                          IF(FOUND) THEN
                            IF(DATAFILE) THEN
C MPN 26Feb98: Don't add offset for collapsed faces (ie. NFF=0)
                              IF(NFF(nf,NE).NE.0) THEN
                                WRITE(IFILE,'(5X,''0 '',I5,'' 0'' )')
     '                            NFF(nf,NE)+offset_face
                              ELSE
                                WRITE(IFILE,'(5X,''0 '',I5,'' 0'' )')
     '                            NFF(nf,NE)
                              ENDIF
C old                              WRITE(IFILE,'(5X,''0 '',I5,'' 0'' )')
C old     '                          NFF(nf,NE)+offset_face
                            ELSE
                              ISOCKET(1)=2
                              ISOCKET(2)=NFF(nf,NE)
                              IF(FSKWRITE(ISOCKET,SK_LONG_INT,2,
     '                          CONNID2).EQ.-1) GOTO 9999
                            ENDIF
                          ELSE
                            IF(DATAFILE) THEN
                              WRITE(IFILE,'(5X,''0     0 0'' )')
                            ELSE
                              ISOCKET(1)=2
                              ISOCKET(2)=0
                              IF(FSKWRITE(ISOCKET,SK_LONG_INT,2,
     '                          CONNID2).EQ.-1) GOTO 9999
                            ENDIF
                          ENDIF
                        ENDDO !i
                      ENDDO !ni
                    ELSE IF(NIT(nb).GE.2) THEN
c cpb 4/8/95 Generalising sector export
C                      CALL ASSERT(NLE(nb).EQ.3,
C     '                  ' Incorrect #lines for element',ERROR,*9999)
C                      IF(DATAFILE) THEN
C                        WRITE(IFILE,'( ''   0 0 0'')')
C                        WRITE(IFILE,'( ''   0 0 '',I5)')
C     '                    NLL(3,NE)+offset_line
C                        WRITE(IFILE,'( ''   0 0 '',I5)')
C     '                    NLL(1,NE)+offset_line
C                        WRITE(IFILE,'( ''   0 0 '',I5)')
C     '                    NLL(2,NE)+offset_line
C                      ELSE
C                        ISOCKET(1)=3
C                        ISOCKET(2)=0
C                        ISOCKET(3)=3
C                        ISOCKET(4)=NLL(3,NE)
C                        ISOCKET(5)=3
C                        ISOCKET(6)=NLL(1,NE)
C                        ISOCKET(7)=3
C                        ISOCKET(8)=NLL(2,NE)
C                        IF(FSKWRITE(ISOCKET,SK_LONG_INT,8,
C     '                      CONNID2).EQ.-1) GOTO 9999
C                     ENDIF


C LKC 26-AUG-2002 Checking that the format size is sufficient.
                      LARGE_FORMAT=.FALSE.
                      DO nae=1,3
                        IF(NLL(nae,NE)+offset_line.GT.99999)
     '                    LARGE_FORMAT=.TRUE.
                      ENDDO
                      IF(LARGE_FORMAT) THEN
                        FORMAT='(''   0 0 '',I9)'
                      ELSE
                        FORMAT='(''   0 0 '',I5)'
                      ENDIF

                      IF(IBT(1,1,nb).EQ.5) THEN
                        IF(DATAFILE) THEN
                          WRITE(IFILE,FORMAT) NLL(2,NE)+offset_line
                          WRITE(IFILE,FORMAT) NLL(3,NE)+offset_line
                          WRITE(IFILE,'( ''   0 0     0'')')
                          WRITE(IFILE,FORMAT) NLL(1,NE)+offset_line
                        ELSE
                          ISOCKET(1)=3
                          ISOCKET(2)=NLL(2,NE)
                          ISOCKET(3)=3
                          ISOCKET(4)=NLL(3,NE)
                          ISOCKET(5)=3
                          ISOCKET(6)=0
                          ISOCKET(7)=3
                          ISOCKET(8)=NLL(1,NE)
                          IF(FSKWRITE(ISOCKET,SK_LONG_INT,8,
     '                      CONNID2).EQ.-1) GOTO 9999
                        ENDIF
                      ELSE IF(IBT(1,1,nb).EQ.6) THEN
                        IF(DATAFILE) THEN
                          WRITE(IFILE,FORMAT) NLL(2,NE)+offset_line
                          WRITE(IFILE,FORMAT) NLL(3,NE)+offset_line
                          WRITE(IFILE,FORMAT) NLL(1,NE)+offset_line
                          WRITE(IFILE,'( ''   0 0     0'')')
                        ELSE
                          ISOCKET(1)=3
                          ISOCKET(2)=NLL(2,NE)
                          ISOCKET(3)=3
                          ISOCKET(4)=NLL(3,NE)
                          ISOCKET(5)=3
                          ISOCKET(6)=NLL(1,NE)
                          ISOCKET(7)=3
                          ISOCKET(8)=0
                          IF(FSKWRITE(ISOCKET,SK_LONG_INT,8,
     '                      CONNID2).EQ.-1) GOTO 9999
                        ENDIF
                      ELSE IF(IBT(1,2,nb).EQ.5) THEN
                        IF(DATAFILE) THEN
                          WRITE(IFILE,'( ''   0 0     0'')')
                          WRITE(IFILE,FORMAT) NLL(3,NE)+offset_line
                          WRITE(IFILE,FORMAT) NLL(1,NE)+offset_line
                          WRITE(IFILE,FORMAT) NLL(2,NE)+offset_line
                        ELSE
                          ISOCKET(1)=3
                          ISOCKET(2)=0
                          ISOCKET(3)=3
                          ISOCKET(4)=NLL(3,NE)
                          ISOCKET(5)=3
                          ISOCKET(6)=NLL(1,NE)
                          ISOCKET(7)=3
                          ISOCKET(8)=NLL(2,NE)
                          IF(FSKWRITE(ISOCKET,SK_LONG_INT,8,
     '                      CONNID2).EQ.-1) GOTO 9999
                        ENDIF
                      ELSE
                        IF(DATAFILE) THEN
                          WRITE(IFILE,FORMAT)
     '                      NLL(3,NE)+offset_line
                          WRITE(IFILE,'( ''   0 0     0'')')
                          WRITE(IFILE,FORMAT) NLL(1,NE)+offset_line
                          WRITE(IFILE,FORMAT) NLL(2,NE)+offset_line
                        ELSE
                          ISOCKET(1)=3
                          ISOCKET(2)=NLL(3,NE)
                          ISOCKET(3)=3
                          ISOCKET(4)=0
                          ISOCKET(5)=3
                          ISOCKET(6)=NLL(1,NE)
                          ISOCKET(7)=3
                          ISOCKET(8)=NLL(2,NE)
                          IF(FSKWRITE(ISOCKET,SK_LONG_INT,8,
     '                      CONNID2).EQ.-1) GOTO 9999
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDIF
C**               write out the grid values
C DB 30-JUL-99 Exports grid point number with geometry
                  IF((NJH_FIELD_BASE_LIST(0).GT.0).OR.
     '              ((ELEM_TYPE(1:8).EQ.'GEOMETRY').AND.(NQT.GT.0)
     '              .AND.GRID_NUMBERS)) THEN
C                  IF((NJH_FIELD_BASE_LIST(0).GT.0)) THEN
                    IF(DATAFILE) THEN
                      WRITE(IFILE,'(3X,''Values:'' )')
                    ELSE
                      ERROR='>>Not implemented'
                      GOTO 9999
                    ENDIF
                    DO no_njh=1,NJH_FIELD_BASE_LIST(0)
                      njh=NJH_FIELD_BASE_LIST(no_njh)
                      IF(ELEM_TYPE(1:8).EQ.'GEOMETRY') THEN
                        ERROR='>>Not implemented'
                        GOTO 9999
                      ELSE IF(ELEM_TYPE(1:5).EQ.'FIELD') THEN
                        nqsc=NQS(ne)
C!!! cpb 25/1/99 THIS IS WRONG BUT NH_LOC IS NOT USED CORRECTLY FOR
C!!! GRID PROBLEMS. I'M LEAVING THIS AS IT IS UNTIL THIS IS FIXED
C                        WRITE(IFILE,'(4X,5(1X,E12.5),:/(4X,5(1X,E12.5)))')
C     '                    (YQ(NYNQ(NH_LOC(njh,nx),NQNE(ne,nqesc),0,nx),
C     '                    iy,na,nx),nqesc=1,NQET(nqsc))

C ***
C *** DPN - Just need to stick a real switch in here, and write out
C ***       the spatially varying parameters
C ***
C ???     - the only change required after NJHT is set-up ???
C ???     - just need to do everything right in EXPORT_FIELD_HEADING ??
C ***
C ***     - something like RCQS_SPATIAL(njh,nq),nq=NQR(1,nr),NQR(2,nr)
C ***

                        IF(EXPORT_CELL) THEN
C *** DPN 17 November 1999 - Write out the grid point values for all
C ***   Spatially varying integer and real parameters and all state
C ***   variables.

C KAT 2006-06-15 something is wrong with k or nqsv here:  see IPMAT3_CELL1.

                          IF(njh.EQ.1) THEN
                            ! model ID - do nothing
                          ELSEIF(njh.EQ.2) THEN
C ***                       The cell type field (variant numbers)
                            WRITE(IFILE,'(4X,5(1X,I12),:/(4X,5(1X,'
     '                        //'I12)))') (ICQS_SPATIAL(1,NQNE(ne,
     '                        nqesc)),nqesc=1,NQET(nqsc))
                          ELSEIF(njh.LT.NQIT+2) THEN
C ***                       The integer parameters
                            j=njh !skip the variant number
                            CELL_FIELD=.FALSE.
                            DO k=1,IICQS_SPATIAL(0,1)
                              IF(IICQS_SPATIAL(k,1).EQ.j) THEN
                                CELL_FIELD=.TRUE.
                                nqsv=k
                              ENDIF
                            ENDDO !k
                            IF(CELL_FIELD) THEN
                              !normal element field
                              WRITE(IFILE,'(4X,5(1X,I12),:/(4X,5(1X,'
     '                          //'I12)))') (ICQS_SPATIAL(nqsv,NQNE(ne,
     '                          nqesc)),nqesc=1,NQET(nqsc))
                            ELSE
                              !indexed field
                              !** do nothing **
                            ENDIF
                          ELSEIF(njh.LT.NQRT+NQIT+2) THEN
C ***                       The real parameters
                            j=njh-1-NQIT
                            CELL_FIELD=.FALSE.
                            DO k=1,IRCQS_SPATIAL(0,1)
                              IF(IRCQS_SPATIAL(k,1).EQ.j) THEN
                                CELL_FIELD=.TRUE.
                                nqsv=k
                              ENDIF
                            ENDDO !k
                            IF(CELL_FIELD) THEN
                              !normal element field
                              WRITE(IFILE,'(4X,5(1X,E12.5),:/(4X,5(1X,'
     '                          //'E12.5)))') (RCQS_SPATIAL(nqsv,
     '                          NQNE(ne,nqesc)),nqesc=1,NQET(nqsc))
                            ELSE
                              !indexed field
                              ! ** do nothing **
                            ENDIF
                          ELSE
C ***                       State variables
                            j=njh-1-NQIT-NQRT
                            WRITE(IFILE,'(4X,5(1X,E12.5),:/(4X,5(1X,'
     '                        //'E12.5)))') (YQS(CELL_STATE_OFFSET(1)-1
     '                        +j,NQNE(ne,nqesc)),nqesc=1,NQET(nqsc))
                          ENDIF
                        ELSE
                          IF(USE_LAT.EQ.0) THEN
                            IF(GRIDVARIABLE(1:3).EQ.'YQS') THEN
                              WRITE(IFILE,'(4X,5(1X,E12.5),:/(4X,5(1X,'
     '                          //'E12.5)))') (YQS(iy,NYNQ(njh,NQNE(ne,
     '                          nqesc),0,nx)),nqesc=1,NQET(nqsc))
                            ELSE
                              WRITE(IFILE,'(4X,5(1X,E12.5),:/(4X,5(1X,'
     '                          //'E12.5)))') (YQ(NYNQ(njh,NQNE(ne,
     '                          nqesc),0,nx),iy,na,nx),nqesc=1,
     '                          NQET(nqsc))
                            ENDIF
                          ELSE !lattice grid method
                            IF(GRIDVARIABLE(1:3).EQ.'YQS') THEN
                              WRITE(IFILE,'(4X,5(1X,E12.5),:/(4X,5(1X,'
     '                          //'E12.5)))') (YQS(iy,NYNQ(njh,
     '                          NQNLAT(nlat),0,nx)),
     '                          nlat=NLATNE(ne),(NLATNE(ne+1)-1))
                            ELSE
                              WRITE(IFILE,'(4X,5(1X,E12.5),:/(4X,5(1X,'
     '                          //'E12.5)))') (YQ(NYNQ(njh,
     '                          NQNLAT(nlat),0,nx),iy,na,nx),
     '                          nlat=NLATNE(ne),(NLATNE(ne+1)-1))
                            ENDIF
                          ENDIF
                        ENDIF
                      ENDIF
                    ENDDO !no_njh
C DB 30-JUL-99 Exports grid point number with geometry
                    IF((ELEM_TYPE(1:8).EQ.'GEOMETRY').AND.(NQT.GT.0)
     '                .AND.GRID_NUMBERS) THEN
                      IF(USE_LAT.EQ.0) THEN
                        nqsc=NQS(ne)
                        WRITE(IFILE,'(8(1X,I9))')
     '                    (NQNE(ne,nqesc),nqesc=1,NQET(nqsc))
                      ELSE !lattice grid method
                        WRITE(IFILE,'(8(1X,I9))') (NQNLAT(nlat),
     '                    nlat=NLATNE(ne),(NLATNE(ne+1)-1))
                      ENDIF
                    ENDIF
                  ENDIF
C**               write the nodes
                  IF(DATAFILE) THEN
                    WRITE(IFILE,'(3X,''Nodes:'' )')
                    WRITE(IFILE,'(4X,16(1X,I12))')
     '                (ELEMENT_NODE_LIST(I)+offset_node,
     '                I=1,ELEMENT_NODE_LIST(0))
C**                 write the scale factors
                    WRITE(IFILE,'(3X,''Scale factors:'' )')
                  ELSE
                    IF(FSKWRITE(ELEMENT_NODE_LIST(1),SK_LONG_INT,
     '                ELEMENT_NODE_LIST(0),CONNID2).EQ.-1) GOTO 9999
                  ENDIF
                  DO nb=1,NBFT
                    IF(NBLIST(nb).NE.0) THEN
                      CALL CALCULATE_BASIS_STRING(IBT,nb,.TRUE.,BASES,
     '                  SPECIAL_BASIS_FLAG,ERROR,*9999)

C 25/2/97 LC archived section :
C         CPB 3/8/95 Sector elements are now the same as the other elements
C         for the ns calculation

C KAT 25Jan00: Exporting scale factors corresponding to a full tensor product.
                      IF(SPECIAL_BASIS_FLAG.GE.3) THEN
                        NN_TOT=1
                        DO ni=1,NIT(nb)
                          IF(IBT(1,ni,nb).EQ.1) THEN !Lagrange
                            NUM_NN=IBT(2,ni,nb)+1
                          ELSE IF(IBT(1,ni,nb).LE.4) THEN
                            !Hermite/Simplex/Seredipity
                            NUM_NN=2
                          ELSE IF(IBT(1,ni,nb).LE.8) THEN !Sector
                            IF(IBT(2,ni,nb).EQ.4) THEN !Hermite
                              NUM_NN=2
                            ELSE !Lagrange
                              NUM_NN=IBT(2,ni,nb)+1
                            ENDIF
                          ELSE !Fourier
                            NUM_NN=IBT(2,ni,nb) !is this right?
                          ENDIF !IBT(1,ni,nb)
                          NUMNODES(ni)=NUM_NN
                          NN_TOT=NN_TOT*NUM_NN
                        ENDDO !ni
                        DO ni=1,3
                          POSITION(ni)=1
                        ENDDO !ni
                        DO nn=1,NN_TOT !nodes in full tensor product
                          nnn= !node in nb
     '                      NNB(POSITION(1),POSITION(2),POSITION(3),nb)
                          CALL ASSERT(nnn.NE.0,
     '                      '>>Could not find local node',ERROR,*9999)
C                         Calculate IDO map of full tensor product
                          NK_TOT=1
                          DO ni=1,NIT(nb)
                            IF((IBT(1,ni,nb).EQ.2.AND. !Hermite
     '                        (IBT(2,ni,nb).EQ.1.OR.
     '                        POSITION(ni).NE.(IBT(2,ni,nb)-1)) !deriv
     '                        ).OR.(
     '                        (IBT(1,ni,nb).EQ.5.OR.IBT(1,ni,nb).EQ.6) !sector
     '                        .AND.IBT(2,ni,nb).EQ.4) !Hermite
     '                        ) THEN
                              DO nk=1,NK_TOT
                                IDOTEMP(nk,ni)=1
                                DO ni2=1,ni-1
                                  IDOTEMP(nk+NK_TOT,ni2)=
     '                              IDOTEMP(nk,ni2)
                                ENDDO !ni2
                              ENDDO !nk
                              DO nk=NK_TOT+1,NK_TOT*2
                                IDOTEMP(nk,ni)=2
                              ENDDO !nk
                              NK_TOT=NK_TOT*2
                            ELSE
                              DO nk=1,NK_TOT
                                IDOTEMP(nk,ni)=1
                              ENDDO !nk
                            ENDIF
                          ENDDO !ni
                          DO ni=NIT(nb)+1,3
                            DO nk=1,NK_TOT
                              IDOTEMP(nk,ni)=1
                            ENDDO !nk
                          ENDDO !ni
                          DO nk=1,NK_TOT !derivs in full tensor product
                            nkk=NKB(IDOTEMP(nk,1),IDOTEMP(nk,2),
     '                        IDOTEMP(nk,3),nnn,nb) !deriv in nb
                            IF(nkk.NE.0) THEN
                              nss=NSB(nkk,nnn,nb) !dof in nb
                              SE_LIST(nk)=SE(nss,nb,ne)
                            ELSE
C                             This dof of the full tensor product is not
C                             used by this basis in cmgui but other
C                             bases may need it.  So, try to calculate
C                             scale factor from lower derivatives.
                              SE_LIST(nk)=1.0d0
                              DO ni=1,NIT(nb)
                                IF(IDOTEMP(nk,ni).EQ.2) THEN
                                  IDOI(ni)=2
                                  nkk2=
     '                              NKB(IDOI(1),IDOI(2),IDOI(3),nnn,nb)
                                  IDOI(ni)=1
                                  IF(nkk2.NE.0) THEN
                                    nss=NSB(nkk2,nnn,nb)
                                    SE_LIST(nk)=
     '                                SE_LIST(nk)*SE(nss,nb,ne)
                                  ENDIF !nkk2
                                ENDIF !IDOTEMP
                              ENDDO !ni
                            ENDIF !nkk
                          ENDDO !nk
                          POSITION(1)=POSITION(1)+1
                          DO ni=1,NIT(nb)
                            IF(POSITION(ni).GT.NUMNODES(ni)) THEN
                              POSITION(ni)=1
                              POSITION(ni+1)=POSITION(ni+1)+1
                            ENDIF
                          ENDDO !ni
                          IF(DATAFILE) THEN
                            WRITE(IFILE,'(4X,4(1X,E24.16))')
     '                        (SE_LIST(nk),nk=1,NK_TOT)
                          ELSE
                            IF(FSKWRITE(SE_LIST,SK_DOUBLE_FLOAT,
     '                        NK_TOT,CONNID2).EQ.-1) GOTO 9999
                          ENDIF
                        ENDDO !nn
                      ELSE
                        IF(DATAFILE) THEN
                          WRITE(IFILE,'(4X,5(1X,E24.16))')
     '                      (SE(nk,nb,NE),nk=1,NST_LIST(nb))
                        ELSE
                          IF(FSKWRITE(SE(1,nb,NE),SK_DOUBLE_FLOAT,
     '                      NST_LIST(nb),CONNID2).EQ.-1) GOTO 9999
                        ENDIF
C                      NS_TOT=0
C                      DO nn=1,NNT(nb)
C                        NS_TOT=NS_TOT+NKT(nn,NB)
C                      ENDDO
C                      IF(DATAFILE) THEN
C                        WRITE(IFILE,'(4X,5(1X,E24.16))')
C     '                    (SE(nk,nb,NE),nk=1,NS_TOT(nb))
C                      ELSE
C                        IF(FSKWRITE(SE(1,nb,NE),SK_DOUBLE_FLOAT,
C     '                    NS_TOT,CONNID2).EQ.-1) GOTO 9999
C                      ENDIF
                      ENDIF
                    ENDIF
                  ENDDO !nb
                  FIRST_ELEMENT=.FALSE.
                ENDIF !NJHT>0

              ENDDO !no_nelist (ne)

C---------------------------------------------------------------------

            ELSE IF(ELEM_TYPE(1:8).EQ.'MATERIAL') THEN

C ***         Write out the variant number(s)
              IF(CELL_VARIANTS_USED) THEN
                !need to write out elem based field
                WRITE(IFILE,'(/'' Variant'')')
                WRITE(IFILE,'(5(1X,I5))') (ICQS_SPATIAL(1,j),j=1,NQT)
              ELSE
                !variant is a constant field
                WRITE(IFILE,'(/A,I5)') ' Variant',
     '            ICQS_SPATIAL(1,1) !use the value for the first g pt
              ENDIF

C *** DPN 21 July 1999 - stuck with only being able to write out one
C *** name when exporting spatially varying parameters, but this is
C *** not a bad thing ???

C ***         Write out the state variables
              !always a element based field ???
              WRITE(IFILE,'(/'' State variables'')')
              DO j=1,CELL_NUM_STATE(1)
                WRITE(IFILE,'(I5,3X,A)') j,
     '            CELL_YQS_NAMES(CELL_STATE_OFFSET(1)-1+j,1)
                WRITE(IFILE,'(5X,5E15.5)') (CELL_YQS_VALUE(
     '            CELL_STATE_OFFSET(1)-1+j,ICQS_SPATIAL(1,k)),k=1,NQT)
              ENDDO !j
C ***         Write out the model parameters
              WRITE(IFILE,'(/'' Model parameters'')')
              DO j=1,CELL_NUM_MODEL(1)
                CELL_FIELD=.FALSE.
                DO k=1,IICQS_SPATIAL(0,1)
                  IF(IICQS_SPATIAL(k,1).EQ.j+CELL_MODEL_OFFSET(1)-1)
     '              THEN
                    CELL_FIELD=.TRUE.
                    nqsv=k
                  ENDIF
                ENDDO !k
                IF(CELL_FIELD) THEN
                  WRITE(IFILE,'(I5,3X,A)') j,
     '              CELL_ICQS_NAMES(CELL_MODEL_OFFSET(1)-1+j,1)
                  WRITE(IFILE,'(5X,5I5)') (ICQS_SPATIAL(nqsv,k),k=1,NQT)
                ELSE
                  DO k=1,CELL_NUM_VARIANTS
                    WRITE(IFILE,'(I5,3X,A,''(variant '',I5,'')'',I5)')
     '                j,CELL_ICQS_NAMES(CELL_MODEL_OFFSET(1)-1+j,k),k,
     '                CELL_ICQS_VALUE(CELL_MODEL_OFFSET(1)-1+j,k)
                  ENDDO !k
                ENDIF !cell_field
              ENDDO !j
C ***         Write out the control parameters
              WRITE(IFILE,'(/'' Control parameters'')')
              DO j=1,CELL_NUM_CONTROL(1)
                CELL_FIELD=.FALSE.
                DO k=1,IICQS_SPATIAL(0,1)
                  IF(IICQS_SPATIAL(k,1).EQ.j+CELL_CONTROL_OFFSET(1)-1)
     '              THEN
                    CELL_FIELD=.TRUE.
                    nqsv=k
                  ENDIF
                ENDDO !k
                IF(CELL_FIELD) THEN
                  WRITE(IFILE,'(I5,3X,A)') j,
     '              CELL_ICQS_NAMES(CELL_CONTROL_OFFSET(1)-1+j,1)
                  WRITE(IFILE,'(5X,5I5)') (ICQS_SPATIAL(nqsv,k),k=1,NQT)
                ELSE
                  DO k=1,CELL_NUM_VARIANTS
                    WRITE(IFILE,'(I5,3X,A,''(variant '',I5,'')'',I5)')
     '                j,CELL_ICQS_NAMES(CELL_CONTROL_OFFSET(1)-1+j,k),k,
     '                CELL_ICQS_VALUE(CELL_CONTROL_OFFSET(1)-1+j,k)
                  ENDDO !k
                ENDIF !cell_field
              ENDDO !j
C ***         Write out the material parameters
              WRITE(IFILE,'(/'' Material parameters'')')
              DO j=1,CELL_NUM_PARAMETERS(1)
                CELL_FIELD=.FALSE.
                DO k=1,IRCQS_SPATIAL(0,1)
                  IF(IRCQS_SPATIAL(k,1).EQ.j+
     '              CELL_PARAMETERS_OFFSET(1)-1) THEN
                    CELL_FIELD=.TRUE.
                    nqsv=k
                  ENDIF
                ENDDO !k
                IF(CELL_FIELD) THEN
                  WRITE(IFILE,'(I5,3X,A)') j,
     '              CELL_RCQS_NAMES(CELL_PARAMETERS_OFFSET(1)-1+j,1)
                  WRITE(IFILE,'(5X,5E15.5)')
     '              (RCQS_SPATIAL(nqsv,k),k=1,NQT)
                ELSE
                  DO k=1,CELL_NUM_VARIANTS
                    WRITE(IFILE,
     '                '(I5,3X,A,''(variant '',I5,'')'',E15.5)') j,
     '                CELL_RCQS_NAMES(CELL_PARAMETERS_OFFSET(1)-1+j,k),
     '                k,CELL_RCQS_VALUE(CELL_PARAMETERS_OFFSET(1)-1+j,k)
                  ENDDO !k
                ENDIF !cell_field
              ENDDO !j
C ***         Write out the protocol parameters
              WRITE(IFILE,'(/'' Protocol parameters'')')
              DO j=1,CELL_NUM_PROTOCOL(1)
                CELL_FIELD=.FALSE.
                DO k=1,IRCQS_SPATIAL(0,1)
                  IF(IRCQS_SPATIAL(k,1).EQ.j+
     '              CELL_PROTOCOL_OFFSET(1)-1) THEN
                    CELL_FIELD=.TRUE.
                    nqsv=k
                  ENDIF
                ENDDO !k
                IF(CELL_FIELD) THEN
                  WRITE(IFILE,'(I5,3X,A)') j,
     '              CELL_RCQS_NAMES(CELL_PROTOCOL_OFFSET(1)-1+j,1)
                  WRITE(IFILE,'(5X,5E15.5)')
     '              (RCQS_SPATIAL(nqsv,k),k=1,NQT)
                ELSE
                  DO k=1,CELL_NUM_VARIANTS
                    WRITE(IFILE,
     '                '(I5,3X,A,''(variant '',I5,'')'',E15.5)') j,
     '                CELL_RCQS_NAMES(CELL_PROTOCOL_OFFSET(1)-1+j,k),k,
     '                CELL_RCQS_VALUE(CELL_PROTOCOL_OFFSET(1)-1+j,k)
                  ENDDO !k
                ENDIF !cell_field
              ENDDO !j
C ***         Write out the AII parameters
              WRITE(IFILE,'(/'' Additional integer input parameters'')')
              DO j=1,CELL_NUM_AII(1)
                CELL_FIELD=.FALSE.
                DO k=1,IICQS_SPATIAL(0,1)
                  IF(IICQS_SPATIAL(k,1).EQ.j+CELL_AII_OFFSET(1)-1) THEN
                    CELL_FIELD=.TRUE.
                    nqsv=k
                  ENDIF
                ENDDO !k
                IF(CELL_FIELD) THEN
                  WRITE(IFILE,'(I5,3X,A)') j,
     '              CELL_ICQS_NAMES(CELL_AII_OFFSET(1)-1+j,1)
                  WRITE(IFILE,'(5X,5I5)') (ICQS_SPATIAL(nqsv,k),k=1,NQT)
                ELSE
                  DO k=1,CELL_NUM_VARIANTS
                    WRITE(IFILE,'(I5,3X,A,''(variant '',I5,'')'',I5)')
     '                j,CELL_ICQS_NAMES(CELL_AII_OFFSET(1)-1+j,k),k,
     '                CELL_ICQS_VALUE(CELL_AII_OFFSET(1)-1+j,k)
                  ENDDO !k
                ENDIF !cell_field
              ENDDO !j
C ***         Write out the ARI parameters
              WRITE(IFILE,'(/'' Additional real input parameters'')')
              DO j=1,CELL_NUM_ARI(1)
                CELL_FIELD=.FALSE.
                DO k=1,IRCQS_SPATIAL(0,1)
                  IF(IRCQS_SPATIAL(k,1).EQ.j+
     '              CELL_ARI_OFFSET(1)-1) THEN
                    CELL_FIELD=.TRUE.
                    nqsv=k
                  ENDIF
                ENDDO !k
                IF(CELL_FIELD) THEN
                  WRITE(IFILE,'(I5,3X,A)') j,
     '              CELL_RCQS_NAMES(CELL_ARI_OFFSET(1)-1+j,1)
                  WRITE(IFILE,'(5X,5E15.5)')
     '              (RCQS_SPATIAL(nqsv,k),k=1,NQT)
                ELSE
                  DO k=1,CELL_NUM_VARIANTS
                    WRITE(IFILE,
     '                '(I5,3X,A,''(variant '',I5,'')'',E15.5)') j,
     '                CELL_RCQS_NAMES(CELL_ARI_OFFSET(1)-1+j,k),k,
     '                CELL_RCQS_VALUE(CELL_ARI_OFFSET(1)-1+j,k)
                  ENDDO !k
                ENDIF !cell_field
              ENDDO !j

C---------------------------------------------------------------------

            ENDIF !ELEM_TYPE=geometry/field/material

            IF(DATAFILE) THEN
              CALL CLOSEF(IFILE,ERROR,*9999)
            ELSE
              IF(FSKWRITE(0,SK_CHAR,1,CONNID2).EQ.-1) GOTO 9999
            ENDIF
          ENDIF !NELIST(0)>0

        ELSE IF(OUTPUT(1:8).EQ.'EXPLORER') THEN

        ENDIF !Output=datafile/Motif/Explorer

      ENDIF !??

      CALL EXITS('EXELEM')
      RETURN
 9999 CALL ERRORS('EXELEM',ERROR)
      CALL EXITS('EXELEM')
      RETURN 1
      END


