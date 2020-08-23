      SUBROUTINE GRGRID(IBT,IDO,INP,LD,NBJ,NBH,NEELEM,NELIST,NENQ,NHE,
     &  NKHE,NKJE,NLATNE,NLATPNQ,NPF,NPNE,NQET,NQLIST,NQNE,NQNLAT,NQS,
     '  NQXI,NRE,NRLIST,NVHE,NVJE,NW,NWQ,NXQ,CURVCORRECT,SE,SQ,XA,XE,
     '  XID,XIQ,XP,XQ,ZA,ZD,ZP,STRING,ERROR,*)

C#### Subroutine: GRGRID
C###  Description:
C###    GRGRID groups grid points.
C**** NTGRGR is number of grid groups currently defined.
C**** LAGRGR(nogrgr) is label given to group number NOGRGR.
C**** LIGRGR(0,nogrgr) is number in list for group number NOGRGR.
C**** LIGRGR(1..,nogrgr) is list for group number NOGRGR.

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'grou00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     &  LD(NDM),NBJ(NJM,NEM),NBH(NHM,NCM,NEM),NEELEM(0:NE_R_M,0:NRM),
     &  NELIST(0:NEM),NENQ(0:8,NQM),NHE(NEM,NXM),NKHE(NKM,NNM,NHM,NEM),
     &  NKJE(NKM,NNM,NJM,NEM),NLATNE(NEQM+1),NLATPNQ(NQM),NPF(9,NFM),
     &  NPNE(NNM,NBFM,NEM),NQET(NQSCM),NQLIST(0:NQM),NQNE(NEQM,NQEM),
     &  NQNLAT(NEQM*NQEM),NQS(NEQM),NQXI(0:NIM,NQSCM),NRE(NEM),
     &  NRLIST(0:NRM),NVHE(NNM,NBFM,NHM,NEM),NVJE(NNM,NBFM,NJM,NEM),
     &  NW(NEM,3,NXM),NWQ(8,0:NQM,NAM),NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),SE(NSM,NBFM,NEM),SQ(NDM),
     &  XA(NAM,NJM,NEM),XE(NSM,NJM),XID(NIM,NDM),XIQ(NIM,NQM),
     &  XP(NKM,NVM,NJM,NPM),XQ(NJM,NQM),ZA(NAM,NHM,NCM,NEM),ZD(NJM,NDM),
     &  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER dirn,ERR,i,IBEG,IBEG1,IBEG2,ICHAR,IEND,IEND1,IEND2,ij,ik,
     &  INFO,IPFILE,j,k,LDTEMP,N3CO,IFROMC,mq,na,ND0,ND1,ne,noelem,
     &  nogr,nogrgr,nolist,NOQUES,nq,nqq,nq_end,nq_start,N1GRGR,ni,nii,
     &  nij,nik,nitb,NJ1,nlat,NOG,nqsc,nr,nrr,ntemp,NTGRGR_READ,
     &  NTR,NWQ_MAP(3,26),pt,nx,Xi_coord,xi_dirn
      REAL*8 distance,origin(3),radius(1),XI(3)
      CHARACTER CHAR*30,CHAR2*2,FILE*100,LABEL*30,STATUS*3
      LOGICAL ALL_REGIONS,CALCU,CBBREV,DEFORM,FILIO,FILEIP,FIRST_TIME,
     &  FOUND,GENER,GRPGRID,INLIST,MOUSE,ONEOFF,PASS2,SPECIFY,UNSORTED
!     Functions
      INTEGER ILISTMBR
      LOGICAL ABBREV,ZEROED

      DATA NWQ_MAP/1,0,0,
     &  -1,0,0,
     &  0,1,0,
     &  0,-1,0,
     &  0,0,1,
     &  0,0,-1,
     &  1,1,0,
     &  1,-1,0,
     &  -1,1,0,
     &  -1,-1,0,
     &  1,0,1,
     &  1,0,-1,
     &  -1,0,1,
     &  -1,0,-1,
     &  0,1,1,
     &  0,1,-1,
     &  0,-1,1,
     &  0,-1,-1,
     &  1,1,1,
     &  1,1,-1,
     &  1,-1,1,
     &  1,-1,-1,
     &  -1,1,1,
     &  -1,1,-1,
     &  -1,-1,1,
     &  -1,-1,-1/     

C!!! DAH 20-08-02 Not all the option implemented for lattice method
C                 yet.      
      
      CALL ENTERS('GRGRID',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
C        CHAR2=CFROMI(NTGRGR+1,'(I2)')
        WRITE(CHAR2,'(I2)') NTGRGR+1
        CALL STRING_TRIM(CHAR2,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM group grids grid GRID#s/GRID_GROUP
C###  Parameter:     <as LABEL[grid_1]>
C###    Specifies the name of the grid group
C###  Parameter:     <(sorted/unsorted)[sorted]>
C###    The sorted command sorts the grid point numbers in ascending
C###    order and removes duplicates.
C###  Parameter:     <region (#s/all)[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped.
C###  Description:
C###    Group grid numbers/groups into a group.  <Unsorted> leaves
C###    specified order and duplicates in list.

        OP_STRING(1)=STRING(1:IEND)//' grid GRID#s/GRID_GROUP'
        OP_STRING(2)=BLANK(1:15)//
     &    '<as LABEL[grid'//CHAR2(IBEG2:IEND2)//']>'
        OP_STRING(3)=BLANK(1:15)//'<(sorted/unsorted)[sorted]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM group grids;l/r/w<;<PATH/>FILENAME>
C###  Parameter:     <region (#s/all)[1]>
C###    Specify the region in which the grid points lie
C###  Description:
C###    Group grid numbers/groups into a group.

        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        OP_STRING(1)=STRING(1:IEND)//';l/r/w'//
     &    '<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
C     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM group grids within elements ELEMENT#s/ELEMENT_GROUP
C###  Parameter:     <as LABEL[grid_1]>
C###    Specifies the name of the grid group
C###  Parameter:     <(sorted/unsorted)[sorted]>
C###    The sorted command sorts the grid point numbers in ascending
C###    order and removes duplicates.
C###  Parameter:     <region (#s/all)[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped.
C###  Description:
C###    Group grids within elements, groups grid points which lie within the
C###    specified elements. <Unsorted> leaves
C###    specified order and duplicates in list.

        OP_STRING(1)=STRING(1:IEND)//' within elements '//
     &    'ELEMENT#s/ELEMENT_GROUP'
        OP_STRING(2)=BLANK(1:15)//
     &    '<as LABEL[grid'//CHAR2(IBEG2:IEND2)//']>'
        OP_STRING(3)=BLANK(1:15)//'<(sorted/unsorted)[sorted]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM group grids
C###  Parameter:     <xi(1/2/3) (high/low)>
C###    Identifies the grid points on a boundary where an xi
C###    coordinate is at its upper or lower limit.  e.g. `xi2 low' indicates
C###    the boundary where the domain is in the positive xi2 direction.
C###  Parameter:     <as LABEL[grid_1]>
C###    Specifies the name of the group grid.
C###  Parameter:     <oneoff>
C###    On global boundary, "oneoff" means one point inside boundary.
C###  Parameter:     <(sorted/unsorted)[sorted]>
C###    The sorted command sorts the grid point numbers in ascending
C###     order and removes duplicates.
C###  Parameter:     <element (#s/all)[all]>
C###    Specify the list of elements to use
C###  Parameter:     <region (#s/all)[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped.
C###  Description:
C###    Group grid numbers/groups into a group.
C###

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<xi(1/2/3) (high/low)>'
        OP_STRING(3)=BLANK(1:15)//
     &    '<as LABEL[grid_'//CHAR2(IBEG2:IEND2)//']>'
        OP_STRING(4)=BLANK(1:15)//'<oneoff>'
        OP_STRING(5)=BLANK(1:15)//'<(sorted/unsorted)[sorted]>'
        OP_STRING(6)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(7)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM group grids external
C###  Parameter:     <as LABEL[grid_1]>
C###    Specifies the name of the grid group.
C###  Parameter:     <oneoff>
C###    On global boundary, "oneoff" means one point inside boundary.
C###  Parameter:     <(sorted/unsorted)[sorted]>
C###    The sorted command sorts the grid point numbers in ascending
C###    order and removes duplicates.
C###  Parameter:     <element (#s/all)[all]>
C###    Specify the list of elements to use
C###  Parameter:     <region (#s/all)[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped.
C###  Description:
C###    Groups all grid points on boundary of specified elements.

        OP_STRING(1)=STRING(1:IEND)//' external'
        OP_STRING(2)=BLANK(1:15)//
     &    '<as LABEL>[grid_'//CHAR2(IBEG2:IEND2)//']'
        OP_STRING(3)=BLANK(1:15)//'<oneoff>'
        OP_STRING(4)=BLANK(1:15)//'<(sorted/unsorted)[sorted]>'
        OP_STRING(5)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(6)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM group grids element NELIST
C###  Parameter:     <as LABEL[grid_1]>
C###    Specifies the name of the grid group.
C###  Parameter:     <(sorted/unsorted)[sorted]>
C###    The sorted command sorts the grid point numbers in ascending
C###    order and removes duplicates.
C###  Parameter:     <region (#s/all)[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped.
C###  Description:
C###    Group all points in element(s) NELIST

        OP_STRING(1)=STRING(1:IEND)//' element NELIST'
        OP_STRING(2)=BLANK(1:15)//
     &    '<as LABEL[grid_'//CHAR2(IBEG2:IEND2)//']>'
        OP_STRING(3)=BLANK(1:15)//'<(sorted/unsorted)[sorted]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM group grids square/cube #
C###  Parameter:     <as LABEL[grid_1]>
C###    Specifies the name of the grid group.
C###  Parameter:     <(sorted/unsorted)[sorted]>
C###    The sorted command sorts the grid point numbers in ascending
C###    order and removes duplicates.
C###  Parameter:     <region (#s/all)[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped.
C###  Description:
C###    Group all points surrounding a single specified grid point
C###    (#) to form a square in 2D or a cube in 3D.

        OP_STRING(1)=STRING(1:IEND)//' square/cube #'
        OP_STRING(2)=BLANK(1:15)//
     &    '<as LABEL[grid'//CHAR2(IBEG2:IEND2)//']>'
        OP_STRING(3)=BLANK(1:15)//'<(sorted/unsorted)[sorted]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM group grids line #
C###  Parameter:     <as LABEL[grid_1]>
C###    Specifies the name of the grid group.
C###  Parameter:     <xidirn #[1]{>=1,<=3}>
C###    Specifies the Xi direction to create the group of grid points
C###    for.
C###  Parameter:     <(positive/negative)[positive]>
C###    Specifies whether the grid points are recorded in either the
C###    positive or negative Xi direction.
C###  Parameter:     <(sorted/unsorted)[sorted]>
C###    The sorted command sorts the grid point numbers in ascending
C###    order and removes duplicates.
C###  Parameter:     <region (#s/all)[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped.
C###  Description:
C###    Start at grid #, and follow xi direction specified until either
C###    returning to same point or reaching boundary.  Gives list of
C###    points in a straight line.

        OP_STRING(1)=STRING(1:IEND)//' line #'
        OP_STRING(2)=BLANK(1:15)//
     &    '<as LABEL>[grid'//CHAR2(IBEG2:IEND2)//']'
        OP_STRING(3)=BLANK(1:15)//'<xidirn #[1]{>=1,<=3}>'
        OP_STRING(4)=BLANK(1:15)//'<(positive/negative)[positive]>'
        OP_STRING(5)=BLANK(1:15)//'<(sorted/unsorted)[sorted]>'
        OP_STRING(6)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM group grids function sphere
C###  Parameter:     <origin X#[0.0],Y#[0.0],Z#[0.0]>
C###   Specify the origin of the sphere within to group the grid points
C###  Parameter:     <radius #[1.0]>
C###   Specify the radius of the sphere within to group the grid points
C###  Parameter:     <as LABEL[grid_1]>
C###    Specifies the name of the group grid.
C###  Parameter:     <(sorted/unsorted)[sorted]>
C###    The sorted command sorts the grid point numbers in ascending
C###     order and removes duplicates.
C###  Parameter:     <region (#s/all)[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped.
C###  Description:
C###    Group grid numbers/groups within a user specified sphere
C###    into a group.

        OP_STRING(1)=STRING(1:IEND)//' function sphere'
        OP_STRING(2)=BLANK(1:15)//'<origin X#[0.0],Y#[0.0],Z#[0.0]>'
        OP_STRING(3)=BLANK(1:15)//'<radius #[1.0]>'
        OP_STRING(4)=BLANK(1:15)//
     &    '<as LABEL[grid_'//CHAR2(IBEG2:IEND2)//']>'
        OP_STRING(5)=BLANK(1:15)//'<(sorted/unsorted)[sorted]>'
        OP_STRING(6)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM group grids face
C###  Parameter:     <xi(1/2/3) (high/low)>
C###    Identifies the grid points on a face where an xi
C###    coordinate is at its upper or lower limit.
C###    e.g. `xi2 low' indicates the boundary where the domain is in
C###    the positive xi2 direction.
C###  Parameter:     <element (#s/all)[all]>
C###    Specify the list of elements to use
C###  Parameter:     <as LABEL[grid_1]>
C###    Specifies the name of the group grid.
C###  Parameter:     <oneoff>
C###    "oneoff" means one point inside boundary.
C###  Parameter:     <(sorted/unsorted)[sorted]>
C###    The sorted command sorts the grid point numbers in ascending
C###     order and removes duplicates.
C###  Parameter:     <region (#s/all)[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###  Description:
C###    Group grid points on a specific face of a list of elements
C###    regardless of whether they are internal or external to the
C###    mesh.

        OP_STRING(1)=STRING(1:IEND)//' face'
        OP_STRING(2)=BLANK(1:15)//'<xi(1/2/3) (high/low)>'
        OP_STRING(3)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(4)=BLANK(1:15)//
     &    '<as LABEL[grid_'//CHAR2(IBEG2:IEND2)//']>'
        OP_STRING(5)=BLANK(1:15)//'<(sorted/unsorted)[sorted]>'
        OP_STRING(6)=BLANK(1:15)//'<oneoff[]>'
        OP_STRING(7)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM group grid ungrouped
C###  Parameter:     <as LABEL[grid_1]>
C###    Specifies the name of the grid group
C###  Parameter:     <region #[1]>
C###    Specify the region in which to group the grid points
C###  Description:
C###    Groups grid points that are not already in a group
C###    into a group.

        OP_STRING(1)=STRING(1:IEND)//' ungrouped'
        OP_STRING(2)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//
     '    '<as LABEL[grid_'//CHAR2(IBEG2:IEND2)//']>'
        OP_STRING(4)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','GRGRID',ERROR,*9999)
      ELSE

C DAH 06-03-03 Added in assert for CALL_GRID
        CALL ASSERT(CALL_GRID,' Must define grid first',
     &    ERROR,*9999)

        IPFILE=1 !is input file version number on 21 Feb 2000
        CALL PARSE_QUALIFIERS(' PLRW',noco,1,CO,COQU,
     &    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

C new CS 20 Feb 2000
        IF(CBBREV(CO,'UNSORTED',1,noco+1,NTCO,N3CO)) THEN
          UNSORTED=.TRUE.
        ELSE
          UNSORTED=.FALSE.
        ENDIF
        Xi_coord=1

C new CS 20 Feb 2000
        IF(FILIO) THEN
          FILEIP=.FALSE.
          NOQUES=0
          CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            ERR=0
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'grgp',
     &        STATUS,ERR,ERROR,*9999)

            IF(IOTYPE.EQ.3) THEN
              IDATA(1)=NTGRGR
            ENDIF
            FORMAT='($,'' The number of grid groups is [1]: '',I5)'
            ICHAR=999 ! Web error unless set
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,
     &        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,NQM,
     &        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              NTGRGR_READ=IDATA(1)
            ELSE
              NTGRGR_READ=NTGRGR
            ENDIF

            DO nogr=1,NTGRGR_READ
              IF(IOTYPE.NE.3) THEN
                CDEFLT(1)='grid1'
                ICHAR=LEN('grid1')
              ELSE
                CDATA(1)=LAGRGR(nogr)
                CALL STRING_TRIM(CDATA(1),IBEG,IEND)
                ICHAR=LEN(CDATA(1)(IBEG:IEND))
              ENDIF
              FORMAT='(/$,'' Grid point group name [grid1]: '',A)'
              CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,1,0,NOQUES,
     &          FILEIP,FORMAT,1,
     &          ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     &          LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
                NTGRGR=NTGRGR+1 !is new total
                CALL ASSERT(NTGRGR.LE.GRGR_MAXGRP,
     &            '>>Can''t create any more groups. '//
     &            'Consider reusing an existing group.',ERROR,*9999)
                CALL STRING_TRIM(CDATA(1),IBEG,IEND)
                LAGRGR(NTGRGR)=CDATA(1)(IBEG:IEND)
              ENDIF

              IF(IOTYPE.EQ.3) THEN
                IDATA(1)=NLIGRGR(nogr)
              ENDIF
              FORMAT='($,'' The number of grid points is [1]: '',I5)'
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &          FORMAT,1,
     &          ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,NQM,
     &          LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
                NLIGRGR(NTGRGR)=IDATA(1)
                NQLIST(0)=NLIGRGR(NTGRGR)
              ENDIF

              IF(IOTYPE.EQ.3) THEN
                CALL ILIST_COPY(NLIGRGR(nogr),
     &            %VAL(LIGRGR_PTR(nogr)),NQLIST(1))
                NQLIST(0)=NLIGRGR(nogr)
              ENDIF

              DO pt=1,NQLIST(0)
                IF(IOTYPE.EQ.3) THEN
                  IDATA(1)=NQLIST(pt)
                ENDIF
                FORMAT='($,'' '',I10)'
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &            FORMAT,1,
     &          ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,NQM,
     &          LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) THEN
                  NQLIST(pt)=IDATA(1)
                ENDIF
              ENDDO

              IF(IOTYPE.NE.3) THEN
                CALL ALLOCATE_MEMORY(NQLIST(0),0,
     &            INTTYPE,LIGRGR_PTR(NTGRGR),MEM_INIT,ERROR,*9999)
                CALL ILIST_COPY(NQLIST(0),NQLIST(1),
     &            %VAL(LIGRGR_PTR(NTGRGR)))
              ENDIF

            ENDDO
            CALL CLOSEF(IFILE,ERROR,*9999)
            IF(BACKUP) IOTYPE=2
          ENDDO

C Points on the faces of a list of elements
C MLB for ryas002, 20-1-03
C MLB adding 1,2d 12-3-03
        ELSE IF(CBBREV(CO,'FACE',2,noco+1,NTCO,N3CO)) THEN

          IF(CBBREV(CO,'XI1',3,noco+1,NTCO,N3CO)
     &      .OR.CBBREV(CO,'XI2',3,noco+1,NTCO,N3CO)
     &      .OR.CBBREV(CO,'XI3',3,noco+1,NTCO,N3CO)) THEN
            Xi_coord=IFROMC(CO(N3CO)(3:3))
            N3CO=N3CO+1
            IF(ABBREV(CO(N3CO),'HIGH',2)) THEN
            ELSE IF(ABBREV(CO(N3CO),'LOW',2)) THEN
              Xi_coord=-Xi_coord
            ELSE
              CALL FLAG_ERROR(1,'Unrecognized xi limit')
              ERROR=' '
              GO TO 9999
            ENDIF
          ELSE
            CALL FLAG_ERROR(1,'Unrecognized xi direction')
            ERROR=' '
            GO TO 9999
          ENDIF

          ONEOFF=CBBREV(CO,'ONEOFF',3,noco+1,NTCO,N3CO)
          IF(ONEOFF.AND.((Xi_coord.EQ.3).OR.(Xi_coord.EQ.-3))) THEN
            CALL ASSERT(.FALSE.,'>>oneoff is not available in 3D',
     &        ERROR,*9999)
          ENDIF

          na=1 !Temporary
          NQLIST(0)=0
          CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,ERROR,
     &      *9999)

          DO noelem=1,NELIST(0) !elements
            ne=NELIST(noelem)

            IF(NIT(NBJ(1,ne)).EQ.3) THEN
              IF(Xi_coord.EQ.-3) THEN
                nq_start=1
                nq_end=NQXI(1,NQS(ne))*NQXI(2,NQS(ne))
                DO nqq=nq_start,nq_end
                  NQLIST(0)=NQLIST(0)+1
                  NQLIST(NQLIST(0))=NQNE(ne,nqq)
                ENDDO !nqq
              ELSEIF(Xi_coord.EQ.-2) THEN
                nq_start=1
                nq_end=NQXI(1,NQS(ne))*NQXI(3,NQS(ne))
                pt=nq_start
                mq=0
                DO nqq=nq_start,nq_end
                  mq=mq+1
                  NQLIST(0)=NQLIST(0)+1
                  NQLIST(NQLIST(0))=NQNE(ne,pt)
                  IF(mq.EQ.NQXI(1,NQS(ne))) THEN
                    mq=0
                    pt=pt+(NQXI(1,NQS(ne))*(NQXI(2,NQS(ne))-1))
                  ENDIF
                  pt=pt+1
                ENDDO !nqq
              ELSEIF(Xi_coord.EQ.-1) THEN
                nq_start=1
                nq_end=(NQXI(1,NQS(ne))*NQXI(2,NQS(ne))*
     &            NQXI(3,NQS(ne)))-(NQXI(1,NQS(ne))-1)
                DO nqq=nq_start,nq_end,NQXI(1,NQS(ne))
                  NQLIST(0)=NQLIST(0)+1
                  NQLIST(NQLIST(0))=NQNE(ne,nqq)
                ENDDO !nqq
              ELSEIF(Xi_coord.EQ.1) THEN
                nq_start=NQXI(1,NQS(ne))
                nq_end=NQXI(1,NQS(ne))*NQXI(2,NQS(ne))*NQXI(3,NQS(ne))
                DO nqq=nq_start,nq_end,NQXI(1,NQS(ne))
                  NQLIST(0)=NQLIST(0)+1
                  NQLIST(NQLIST(0))=NQNE(ne,nqq)
                ENDDO !nqq
              ELSEIF(Xi_coord.EQ.2) THEN
                nq_start=(NQXI(1,NQS(ne))*(NQXI(2,NQS(ne))-1))+1
                nq_end=NQXI(1,NQS(ne))*NQXI(3,NQS(ne))
                pt=nq_start
                mq=0
                DO nqq=1,nq_end
                  mq=mq+1
                  NQLIST(0)=NQLIST(0)+1
                  NQLIST(NQLIST(0))=NQNE(ne,pt)
                  IF(mq.EQ.NQXI(1,NQS(ne))) THEN
                    mq=0
                    pt=pt+NQXI(1,NQS(ne))*(NQXI(2,NQS(ne))-1)
                  ENDIF
                  pt=pt+1
                ENDDO !nqq
              ELSEIF(Xi_coord.EQ.3) THEN
                nq_start=(NQXI(1,NQS(ne))*NQXI(2,NQS(ne))*
     &            (NQXI(3,NQS(ne))-1))+1
                nq_end=NQXI(1,NQS(ne))*NQXI(2,NQS(ne))*NQXI(3,NQS(ne))
                DO nqq=nq_start,nq_end
                  NQLIST(0)=NQLIST(0)+1
                  NQLIST(NQLIST(0))=NQNE(ne,nqq)
                ENDDO !nqq
              ENDIF

            ELSE IF(NIT(NBJ(1,ne)).EQ.2) THEN

              IF(Xi_coord.EQ.-2) THEN
                IF(ONEOFF) THEN
                  nq_start=NQXI(1,NQS(ne))+1
                  nq_end=NQXI(1,NQS(ne))*2
                ELSE
                  nq_start=1
                  nq_end=NQXI(1,NQS(ne))
                ENDIF
                DO nqq=nq_start,nq_end
                  NQLIST(0)=NQLIST(0)+1
                  NQLIST(NQLIST(0))=NQNE(ne,nqq)
                ENDDO !nqq
              ELSEIF(Xi_coord.EQ.-1) THEN
                nq_start=1
                nq_end=NQXI(2,NQS(ne))
                IF(ONEOFF) THEN
                  DO nqq=nq_start,nq_end
                    mq=(nqq-1)*NQXI(1,NQS(ne))+2
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=NQNE(ne,mq)
                  ENDDO !nqq
                ELSE
                  DO nqq=nq_start,nq_end
                    mq=(nqq-1)*NQXI(1,NQS(ne))+1
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=NQNE(ne,mq)
                  ENDDO !nqq
                ENDIF
              ELSEIF(Xi_coord.EQ.1) THEN
                nq_start=1
                nq_end=NQXI(2,NQS(ne))
                IF(ONEOFF) THEN
                  DO nqq=1,nq_end
                    mq=(nqq*NQXI(1,NQS(ne)))-1
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=NQNE(ne,mq)
                  ENDDO !nqq
                ELSE
                  DO nqq=1,nq_end
                    mq=nqq*NQXI(1,NQS(ne))
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=NQNE(ne,mq)
                  ENDDO !nqq
                ENDIF
              ELSEIF(Xi_coord.EQ.2) THEN
                IF(ONEOFF) THEN
                  nq_start=((NQXI(2,NQS(ne))-2)*NQXI(1,NQS(ne)))+1
                  nq_end=((NQXI(2,NQS(ne))-1)*NQXI(1,NQS(ne)))
                ELSE
                  nq_start=((NQXI(2,NQS(ne))-1)*NQXI(1,NQS(ne)))+1
                  nq_end=NQXI(1,NQS(ne))*NQXI(2,NQS(ne))
                ENDIF
                DO nqq=nq_start,nq_end
                  NQLIST(0)=NQLIST(0)+1
                  NQLIST(NQLIST(0))=NQNE(ne,nqq)
                ENDDO !nqq
              ENDIF
   
            ELSE IF(NIT(NBJ(1,ne)).EQ.1) THEN

              IF(Xi_coord.EQ.-1) THEN
                NQLIST(0)=NQLIST(0)+1
                IF(ONEOFF) THEN
                  NQLIST(NQLIST(0))=NQNE(ne,2)
                ELSE
                  NQLIST(NQLIST(0))=NQNE(ne,1)
                ENDIF
              ELSEIF(Xi_coord.EQ.1) THEN
                NQLIST(0)=NQLIST(0)+1
                IF(ONEOFF) THEN
                  NQLIST(NQLIST(0))=NQNE(ne,(NQXI(1,NQS(ne))-1))
                ELSE
                  NQLIST(NQLIST(0))=NQNE(ne,NQXI(1,NQS(ne)))
                ENDIF
              ENDIF
            ENDIF

          ENDDO !noelem

          WRITE(OP_STRING(1),
     &      '(''>>Increase NQM >= '',I12)') NQLIST(0)
          CALL STRING_TRIM(OP_STRING(1),IBEG,IEND)
          CALL ASSERT(NQLIST(0).LE.NQM,OP_STRING(1)(IBEG:IEND),ERROR,
     &      *9999)

C new CS 18 Feb 2000
        ELSE IF(CBBREV(CO,'WITHIN',2,noco+1,NTCO,N3CO)) THEN
          CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     &      ERROR,*9999)
          CALL ASSERT(NET(NRLIST(1)).GT.0,
     &      '>>no elements defined',ERROR,*9999)
          nr=NRLIST(1)
          NJ1=1

          GRPGRID=.TRUE.
          NITB=NIT(NBJ(1,NELIST(1)))
          SPECIFY=.FALSE.
          DEFORM=.FALSE.


          DO nq=1,NQT
            ND0=nq
            ND1=nq
            LD(1)=0
            SQ(1)=0.d0
            nx=1
            PASS2=.TRUE.
            CALL DEXI_NONLIN(IBT,IDO,INP,LD,LDTEMP,NBJ,NBH,
     &        nq,ND0,ND1,ne,NELIST,NHE,ni,NITB,nj1,NKHE,NKJE,NPF,
     &        NPNE,nolist,nr,NRE,NVHE,NVJE,NW(1,1,nx),nx,CURVCORRECT,
     &        SE,XA,XE,XI,XID,XIQ,XP,XQ,ZA,ZD,ZP,DEFORM,
     &        FOUND,GRPGRID,PASS2,SPECIFY,ERROR,*9999)

            IF(LD(1).NE.0) THEN
              NQLIST(0)=NQLIST(0)+1
              NQLIST(NQLIST(0))=nq
            ENDIF

          ENDDO !nq

          WRITE(OP_STRING(1),
     &      '(''>>Increase NQM >= '',I12)') NQLIST(0)
          CALL STRING_TRIM(OP_STRING(1),IBEG,IEND)
          CALL ASSERT(NQLIST(0).LE.NQM,OP_STRING(1)(IBEG:IEND),ERROR,
     &      *9999)

C  Grid points on Xi# boundary
        ELSE IF(CBBREV(CO,'XI1',3,noco+1,NTCO,N3CO)
     &    .OR.CBBREV(CO,'XI2',3,noco+1,NTCO,N3CO)
     &    .OR.CBBREV(CO,'XI3',3,noco+1,NTCO,N3CO)) THEN
          Xi_coord=IFROMC(CO(N3CO)(3:3))
          N3CO=N3CO+1
C KAT 14Sep00: 0,1 are obselete as they are misleading
          IF(CO(N3CO).EQ.'1'.OR.ABBREV(CO(N3CO),'HIGH',2)) THEN
          ELSEIF(CO(N3CO).EQ.'0'.OR.ABBREV(CO(N3CO),'LOW',2)) THEN
            Xi_coord=-Xi_coord
          ELSE
            CALL FLAG_ERROR(1,'Unrecognized xi limit')
            ERROR=' '
            GO TO 9999
          ENDIF
          ONEOFF=CBBREV(CO,'ONEOFF',3,noco+1,NTCO,N3CO)
          na=1 !Temporary
          NQLIST(0)=0
          IF(CBBREV(CO,'ELEMENT',2,noco+1,NTCO,N3CO)) THEN
            IF(USE_LAT.EQ.1) THEN
              WRITE(OP_STRING,'(''>>Warning: The lattice grid scheme'//
     &          ' defines element grid point hosting differently'//
     &          ' from the standard collocation scheme.'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     &        ERROR,*9999)
            DO noelem=1,NELIST(0) !elements
              ne=NELIST(noelem)
              IF(USE_LAT.EQ.0) THEN
                DO nqq=1,NQET(NQS(ne))
                  nq=NQNE(ne,nqq)
                  IF(NXQ(Xi_coord,1,nq,na).EQ.0) THEN
                    NQLIST(0)=NQLIST(0)+1
                    IF(NQLIST(0).LE.NQM) THEN
                      IF(ONEOFF) THEN
                        NQLIST(NQLIST(0))=NWQ(1,nq,na) !Point one in from bdy
                      ELSE
                        NQLIST(NQLIST(0))=nq
                      ENDIF
                    ENDIF
                  ENDIF !NXQ
                ENDDO !nqq
              ELSE !lattice
                nr=NRE(ne)
                DO nq=NQR(1,nr),NQR(2,nr)
                  IF(NWQ(1,nq,na).NE.0) THEN                  
                    IF(NENQ(1,nq).EQ.ne) THEN
                      IF(ONEOFF) THEN
                        IF(Xi_coord.EQ.1) THEN
                          IF(NWQ_MAP(1,NWQ(1,nq,na)).EQ.-1) THEN
                            NQLIST(0)=NQLIST(0)+1                          
                            nlat=NLATPNQ(nq)
                            CALL FIND_LATT_NEIJK(i,j,k,ne,nlat,NLATNE,
     &                        NQXI,NQS,nqsc,ERROR,*9999)
                            NQLIST(NQLIST(0))= NQNLAT(NLATNE(ne)+
     &                        ((k-1+NWQ_MAP(3,NWQ(1,nq,na)))*NQXI(2,
     &                        nqsc)+j-1+NWQ_MAP(2,NWQ(1,nq,na)))*NQXI(1,
     &                        nqsc)+i-1+NWQ_MAP(1,NWQ(1,nq,na)))
                          ENDIF                        
                        ELSEIF(Xi_coord.EQ.-1) THEN
                          IF(NWQ_MAP(1,NWQ(1,nq,na)).EQ.1) THEN
                            NQLIST(0)=NQLIST(0)+1
                            nlat=NLATPNQ(nq)
                            CALL FIND_LATT_NEIJK(i,j,k,ne,nlat,NLATNE,
     &                        NQXI,NQS,nqsc,ERROR,*9999)
                            NQLIST(NQLIST(0))= NQNLAT(NLATNE(ne)+
     &                        ((k-1+NWQ_MAP(3,NWQ(1,nq,na)))*NQXI(2,
     &                        nqsc)+j-1+NWQ_MAP(2,NWQ(1,nq,na)))*NQXI(1,
     &                        nqsc)+i-1+NWQ_MAP(1,NWQ(1,nq,na)))
                          ENDIF
                        ELSEIF(Xi_coord.EQ.2) THEN
                          IF(NWQ_MAP(2,NWQ(1,nq,na)).EQ.-1) THEN
                            NQLIST(0)=NQLIST(0)+1
                            nlat=NLATPNQ(nq)
                            CALL FIND_LATT_NEIJK(i,j,k,ne,nlat,NLATNE,
     &                        NQXI,NQS,nqsc,ERROR,*9999)
                            NQLIST(NQLIST(0))= NQNLAT(NLATNE(ne)+
     &                        ((k-1+NWQ_MAP(3,NWQ(1,nq,na)))*NQXI(2,
     &                        nqsc)+j-1+NWQ_MAP(2,NWQ(1,nq,na)))*NQXI(1,
     &                        nqsc)+i-1+NWQ_MAP(1,NWQ(1,nq,na)))
                          ENDIF                           
                        ELSEIF(Xi_coord.EQ.-2) THEN
                          IF(NWQ_MAP(2,NWQ(1,nq,na)).EQ.1) THEN
                            NQLIST(0)=NQLIST(0)+1
                            nlat=NLATPNQ(nq)
                            CALL FIND_LATT_NEIJK(i,j,k,ne,nlat,NLATNE,
     &                        NQXI,NQS,nqsc,ERROR,*9999)
                            NQLIST(NQLIST(0))= NQNLAT(NLATNE(ne)+
     &                        ((k-1+NWQ_MAP(3,NWQ(1,nq,na)))*NQXI(2,
     &                        nqsc)+j-1+NWQ_MAP(2,NWQ(1,nq,na)))*NQXI(1,
     &                        nqsc)+i-1+NWQ_MAP(1,NWQ(1,nq,na)))
                          ENDIF                         
                        ELSEIF(Xi_coord.EQ.3) THEN
                          IF(NWQ_MAP(3,NWQ(1,nq,na)).EQ.-1) THEN
                            NQLIST(0)=NQLIST(0)+1
                            nlat=NLATPNQ(nq)
                            CALL FIND_LATT_NEIJK(i,j,k,ne,nlat,NLATNE,
     &                        NQXI,NQS,nqsc,ERROR,*9999)
                            NQLIST(NQLIST(0))= NQNLAT(NLATNE(ne)+
     &                        ((k-1+NWQ_MAP(3,NWQ(1,nq,na)))*NQXI(2,
     &                        nqsc)+j-1+NWQ_MAP(2,NWQ(1,nq,na)))*NQXI(1,
     &                        nqsc)+i-1+NWQ_MAP(1,NWQ(1,nq,na)))
                          ENDIF                           
                        ELSE
                          IF(NWQ_MAP(3,NWQ(1,nq,na)).EQ.1) THEN
                            NQLIST(0)=NQLIST(0)+1
                            nlat=NLATPNQ(nq)
                            CALL FIND_LATT_NEIJK(i,j,k,ne,nlat,NLATNE,
     &                        NQXI,NQS,nqsc,ERROR,*9999)
                            NQLIST(NQLIST(0))= NQNLAT(NLATNE(ne)+
     &                        ((k-1+NWQ_MAP(3,NWQ(1,nq,na)))*NQXI(2,
     &                        nqsc)+j-1+NWQ_MAP(2,NWQ(1,nq,na)))*NQXI(1,
     &                        nqsc)+i-1+NWQ_MAP(1,NWQ(1,nq,na)))
                          ENDIF
                        ENDIF
                      ELSE
                        IF(Xi_coord.EQ.1) THEN
                          IF(NWQ_MAP(1,NWQ(1,nq,na)).EQ.-1) THEN
                            NQLIST(0)=NQLIST(0)+1 
                            NQLIST(NQLIST(0))=nq
                          ENDIF                        
                        ELSEIF(Xi_coord.EQ.-1) THEN
                          IF(NWQ_MAP(1,NWQ(1,nq,na)).EQ.1) THEN
                            NQLIST(0)=NQLIST(0)+1
                            NQLIST(NQLIST(0))=nq
                          ENDIF
                        ELSEIF(Xi_coord.EQ.2) THEN
                          IF(NWQ_MAP(2,NWQ(1,nq,na)).EQ.-1) THEN
                            NQLIST(0)=NQLIST(0)+1 
                            NQLIST(NQLIST(0))=nq
                          ENDIF                           
                        ELSEIF(Xi_coord.EQ.-2) THEN
                          IF(NWQ_MAP(2,NWQ(1,nq,na)).EQ.1) THEN
                            NQLIST(0)=NQLIST(0)+1
                            NQLIST(NQLIST(0))=nq                         
                          ENDIF                         
                        ELSEIF(Xi_coord.EQ.3) THEN
                          IF(NWQ_MAP(3,NWQ(1,nq,na)).EQ.-1) THEN
                            NQLIST(0)=NQLIST(0)+1                
                            NQLIST(NQLIST(0))=nq                         
                          ENDIF                           
                        ELSE
                          IF(NWQ_MAP(3,NWQ(1,nq,na)).EQ.1) THEN
                            NQLIST(0)=NQLIST(0)+1                  
                            NQLIST(NQLIST(0))=nq                         
                          ENDIF                         
                        ENDIF                   
                      ENDIF
                    ENDIF
                  ENDIF
                ENDDO
              ENDIF
            ENDDO !noelem
          ELSE
            DO nrr=1,NRLIST(0)
              nr=NRLIST(nrr)
              DO nq=NQR(1,nr),NQR(2,nr)
                IF(USE_LAT.EQ.0) THEN
                  IF(NXQ(Xi_coord,1,nq,na).EQ.0) THEN
                    NQLIST(0)=NQLIST(0)+1
                    IF(NQLIST(0).LE.NQM) THEN
                      IF(ONEOFF) THEN
                        NQLIST(NQLIST(0))=NWQ(1,nq,na) 
!Point one in from bdy
                      ELSE
                        NQLIST(NQLIST(0))=nq
                      ENDIF
                    ENDIF
                  ENDIF
                ELSE !use lattice
                  IF(ONEOFF) THEN
                    IF(NWQ(1,nq,na).NE.0) THEN
                      IF(Xi_coord.EQ.1) THEN
                        IF(NWQ_MAP(1,NWQ(1,nq,na)).EQ.-1) THEN
                          NQLIST(0)=NQLIST(0)+1                          
                          nlat=NLATPNQ(nq)
                          CALL FIND_LATT_NEIJK(i,j,k,ne,nlat,NLATNE,
     &                      NQXI,NQS,nqsc,ERROR,*9999)
                          NQLIST(NQLIST(0))= NQNLAT(NLATNE(ne)+
     &                      ((k-1+NWQ_MAP(3,NWQ(1,nq,na)))*NQXI(2,nqsc)+
     &                      j-1+NWQ_MAP(2,NWQ(1,nq,na)))*NQXI(1,nqsc)+
     &                      i-1+NWQ_MAP(1,NWQ(1,nq,na)))
                        ENDIF                        
                      ELSEIF(Xi_coord.EQ.-1) THEN
                        IF(NWQ_MAP(1,NWQ(1,nq,na)).EQ.1) THEN
                          NQLIST(0)=NQLIST(0)+1
                          nlat=NLATPNQ(nq)
                          CALL FIND_LATT_NEIJK(i,j,k,ne,nlat,NLATNE,
     &                      NQXI,NQS,nqsc,ERROR,*9999)
                          NQLIST(NQLIST(0))= NQNLAT(NLATNE(ne)+
     &                      ((k-1+NWQ_MAP(3,NWQ(1,nq,na)))*NQXI(2,nqsc)+
     &                      j-1+NWQ_MAP(2,NWQ(1,nq,na)))*NQXI(1,nqsc)+
     &                      i-1+NWQ_MAP(1,NWQ(1,nq,na)))
                        ENDIF
                      ELSEIF(Xi_coord.EQ.2) THEN
                        IF(NWQ_MAP(2,NWQ(1,nq,na)).EQ.-1) THEN
                          NQLIST(0)=NQLIST(0)+1
                          nlat=NLATPNQ(nq)
                          CALL FIND_LATT_NEIJK(i,j,k,ne,nlat,NLATNE,
     &                      NQXI,NQS,nqsc,ERROR,*9999)
                          NQLIST(NQLIST(0))= NQNLAT(NLATNE(ne)+
     &                      ((k-1+NWQ_MAP(3,NWQ(1,nq,na)))*NQXI(2,nqsc)+
     &                      j-1+NWQ_MAP(2,NWQ(1,nq,na)))*NQXI(1,nqsc)+
     &                      i-1+NWQ_MAP(1,NWQ(1,nq,na)))
                        ENDIF                           
                      ELSEIF(Xi_coord.EQ.-2) THEN
                        IF(NWQ_MAP(2,NWQ(1,nq,na)).EQ.1) THEN
                          NQLIST(0)=NQLIST(0)+1
                          nlat=NLATPNQ(nq)
                          CALL FIND_LATT_NEIJK(i,j,k,ne,nlat,NLATNE,
     &                      NQXI,NQS,nqsc,ERROR,*9999)
                          NQLIST(NQLIST(0))= NQNLAT(NLATNE(ne)+
     &                      ((k-1+NWQ_MAP(3,NWQ(1,nq,na)))*NQXI(2,nqsc)+
     &                      j-1+NWQ_MAP(2,NWQ(1,nq,na)))*NQXI(1,nqsc)+
     &                      i-1+NWQ_MAP(1,NWQ(1,nq,na)))
                        ENDIF                         
                      ELSEIF(Xi_coord.EQ.3) THEN
                        IF(NWQ_MAP(3,NWQ(1,nq,na)).EQ.-1) THEN
                          NQLIST(0)=NQLIST(0)+1
                          nlat=NLATPNQ(nq)
                          CALL FIND_LATT_NEIJK(i,j,k,ne,nlat,NLATNE,
     &                      NQXI,NQS,nqsc,ERROR,*9999)
                          NQLIST(NQLIST(0))= NQNLAT(NLATNE(ne)+
     &                      ((k-1+NWQ_MAP(3,NWQ(1,nq,na)))*NQXI(2,nqsc)+
     &                      j-1+NWQ_MAP(2,NWQ(1,nq,na)))*NQXI(1,nqsc)+
     &                      i-1+NWQ_MAP(1,NWQ(1,nq,na)))
                        ENDIF                           
                      ELSE
                        IF(NWQ_MAP(3,NWQ(1,nq,na)).EQ.1) THEN
                          NQLIST(0)=NQLIST(0)+1
                          nlat=NLATPNQ(nq)
                          CALL FIND_LATT_NEIJK(i,j,k,ne,nlat,NLATNE,
     &                      NQXI,NQS,nqsc,ERROR,*9999)
                          NQLIST(NQLIST(0))= NQNLAT(NLATNE(ne)+
     &                      ((k-1+NWQ_MAP(3,NWQ(1,nq,na)))*NQXI(2,nqsc)+
     &                      j-1+NWQ_MAP(2,NWQ(1,nq,na)))*NQXI(1,nqsc)+
     &                      i-1+NWQ_MAP(1,NWQ(1,nq,na)))
                        ENDIF                         
                      ENDIF
                    ENDIF
                  ELSE
                    IF(NWQ(1,nq,na).NE.0) THEN
                      IF(Xi_coord.EQ.1) THEN
                        IF(NWQ_MAP(1,NWQ(1,nq,na)).EQ.-1) THEN
                          NQLIST(0)=NQLIST(0)+1 
                          NQLIST(NQLIST(0))=nq
                        ENDIF                        
                      ELSEIF(Xi_coord.EQ.-1) THEN
                        IF(NWQ_MAP(1,NWQ(1,nq,na)).EQ.1) THEN
                          NQLIST(0)=NQLIST(0)+1
                          NQLIST(NQLIST(0))=nq
                        ENDIF
                      ELSEIF(Xi_coord.EQ.2) THEN
                        IF(NWQ_MAP(2,NWQ(1,nq,na)).EQ.-1) THEN
                          NQLIST(0)=NQLIST(0)+1 
                          NQLIST(NQLIST(0))=nq
                        ENDIF                           
                      ELSEIF(Xi_coord.EQ.-2) THEN
                        IF(NWQ_MAP(2,NWQ(1,nq,na)).EQ.1) THEN
                          NQLIST(0)=NQLIST(0)+1
                          NQLIST(NQLIST(0))=nq                         
                        ENDIF                         
                      ELSEIF(Xi_coord.EQ.3) THEN
                        IF(NWQ_MAP(3,NWQ(1,nq,na)).EQ.-1) THEN
                          NQLIST(0)=NQLIST(0)+1                
                          NQLIST(NQLIST(0))=nq                         
                        ENDIF                           
                      ELSE
                        IF(NWQ_MAP(3,NWQ(1,nq,na)).EQ.1) THEN
                          NQLIST(0)=NQLIST(0)+1                  
                          NQLIST(NQLIST(0))=nq                         
                        ENDIF                         
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF !use_lat
              ENDDO !nq
            ENDDO !nr
          ENDIF !element
          WRITE(OP_STRING(1),
     &      '(''>>Increase NQM >= '',I12)') NQLIST(0)
          CALL STRING_TRIM(OP_STRING(1),IBEG,IEND)
          CALL ASSERT(NQLIST(0).LE.NQM,OP_STRING(1)(IBEG:IEND),ERROR,
     &      *9999)

C  All external grid points
        ELSE IF(CBBREV(CO,'EXTERNAL',2,noco+1,NTCO,N3CO)) THEN
          ONEOFF=.FALSE.
          na=1 !Temporary
          IF(CBBREV(CO,'ONEOFF',3,noco+1,NTCO,N3CO)) THEN
            ONEOFF=.TRUE.
          ENDIF
          NQLIST(0)=0
          IF(CBBREV(CO,'ELEMENT',2,noco+1,NTCO,N3CO)) THEN
            CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     &        ERROR,*9999)
            DO noelem=1,NELIST(0) !elements
              ne=NELIST(noelem)
              DO nqq=1,NQET(NQS(ne))
                nq=NQNE(ne,nqq)
                IF(NWQ(1,nq,1).NE.0) THEN
                  NQLIST(0)=NQLIST(0)+1
                  IF(NQLIST(0).LE.NQM) THEN
                    IF(ONEOFF) THEN
                      IF(USE_LAT.EQ.0) THEN
                        NQLIST(NQLIST(0))=NWQ(1,nq,na) !Point one in from bdy
                      ELSE
                        nlat=NLATPNQ(nq)
                        CALL FIND_LATT_NEIJK(i,j,k,ne,nlat,NLATNE,NQXI,
     &                    NQS,nqsc,ERROR,*9999)
                        NQLIST(NQLIST(0))= NQNLAT(NLATNE(ne)+
     &                    ((k-1+NWQ_MAP(3,NWQ(1,nq,na)))*NQXI(2,nqsc)+
     &                    j-1+NWQ_MAP(2,NWQ(1,nq,na)))*NQXI(1,nqsc)+
     &                    i-1+NWQ_MAP(1,NWQ(1,nq,na)))
                      ENDIF
                    ELSE
                      NQLIST(NQLIST(0))=nq
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO !nqq
            ENDDO !noelem
          ELSE
            DO nrr=1,NRLIST(0)
              nr=NRLIST(nrr)
              DO nq=NQR(1,nr),NQR(2,nr)
                IF(NWQ(1,nq,1).NE.0) THEN
                  NQLIST(0)=NQLIST(0)+1
                  IF(NQLIST(0).LE.NQM) THEN
                    IF(ONEOFF) THEN
                      IF(USE_LAT.EQ.0) THEN
                        NQLIST(NQLIST(0))=NWQ(1,nq,na) !Point one in from bdy
                      ELSE
                        nlat=NLATPNQ(nq)
                        CALL FIND_LATT_NEIJK(i,j,k,ne,nlat,NLATNE,NQXI,
     &                    NQS,nqsc,ERROR,*9999)                        
                        NQLIST(NQLIST(0))= NQNLAT(NLATNE(ne)+
     &                    ((k-1+NWQ_MAP(3,NWQ(1,nq,na)))*NQXI(2,nqsc)+
     &                    j-1+NWQ_MAP(2,NWQ(1,nq,na)))*NQXI(1,nqsc)+
     &                    i-1+NWQ_MAP(1,NWQ(1,nq,na)))
                      ENDIF
                    ELSE
                      NQLIST(NQLIST(0))=nq
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO !nq
            ENDDO !nr
          ENDIF
          WRITE(OP_STRING(1),
     &      '(''>>Increase NQM >= '',I12)') NQLIST(0)
          CALL STRING_TRIM(OP_STRING(1),IBEG,IEND)
          CALL ASSERT(NQLIST(0).LE.NQM,OP_STRING(1)(IBEG:IEND),ERROR,
     &      *9999)

C Points contained within specified element(s)
        ELSE IF(CBBREV(CO,'ELEMENT',2,noco+1,NTCO,N3CO)) THEN
          CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     &      ERROR,*9999)

C MLB 18 August 1997
C old
C          IF(CBBREV(CO,'BASIS',2,noco+1,NTCO,N3CO)) THEN
C            nb_extended=IFROMC(CO(N3CO+1))
C          ELSE
C            nb_extended=1
C            DO WHILE (nb_extended.LE.NBT.AND.NBC(nb_extended).NE.7)
C              nb_extended=nb_extended+1
C            ENDDO
C          ENDIF
C          CALL ASSERT((NBC(nb_extended).EQ.7),
C     '      'Extended basis function not defined',ERROR,*9999)
C          NQLIST(0)=0
C          DO noelem=1,NELIST(0) !element list
C            ne=NELIST(noelem)
C            DO ng=1,NGT(nb_extended)
C              NQLIST(0)=NQLIST(0)+1
C              NQLIST(NQLIST(0))=NQGE(ng,ne,nb_extended)
C            ENDDO
C          ENDDO
          NQLIST(0)=0
          IF(USE_LAT.EQ.0) THEN
            DO noelem=1,NELIST(0) !element list
              ne=NELIST(noelem)
              DO nqq=1,NQET(NQS(ne))
                nq=NQNE(ne,nqq)
                NQLIST(0)=NQLIST(0)+1
                NQLIST(NQLIST(0))=nq
              ENDDO
            ENDDO
          ELSE
            DO noelem=1,NELIST(0) !element list
              ne=NELIST(noelem)
              DO nqq=NLATNE(ne),NLATNE(ne+1)-1
                nq=NQNLAT(nqq)
                NQLIST(0)=NQLIST(0)+1
                NQLIST(NQLIST(0))=nq
              ENDDO
            ENDDO
          ENDIF

          WRITE(OP_STRING(1),
     &      '(''>>Increase NQM >= '',I12)') NQLIST(0)
          CALL STRING_TRIM(OP_STRING(1),IBEG,IEND)
          CALL ASSERT(NQLIST(0).LE.NQM,OP_STRING(1)(IBEG:IEND),ERROR,
     &      *9999)

C Surrounding points in a square/cube
        ELSE IF(CBBREV(CO,'SQUARE',1,noco+1,NTCO,N3CO).OR.
     &     CBBREV(CO,'CUBE',1,noco+1,NTCO,N3CO)) THEN
          nq=IFROMC(CO(N3CO+1))
          IF(CO(N3CO)(1:1).EQ.'S'.OR.CO(N3CO)(1:1).EQ.'s') THEN
            NITB=2
          ELSE IF(CO(N3CO)(1:1).EQ.'C'.OR.CO(N3CO)(1:1).EQ.'c') THEN
            NITB=3
          ENDIF
          ik=MAX(0,NITB-2) !zero for 1,2-D, one for 3-D
          ij=MIN(NITB-1,1) !zero for 1-D, one for 2,3-D
          NQLIST(0)=0
          na=1 !Temporary
          DO nik=-ik,ik
            DO nij=-ij,ij
              DO nii=-1,1
                mq=NXQ(nii,1,NXQ(nij*2,1,NXQ(nik*3,1,nq,na),na),na)
                IF(mq.GT.0) THEN
                  NQLIST(0)=NQLIST(0)+1
                  IF(NQLIST(0).LE.NQM) NQLIST(NQLIST(0))=mq
                ENDIF
              ENDDO !nii
            ENDDO !nij
          ENDDO !nik
          WRITE(OP_STRING(1),
     &      '(''>>Increase NQM >= '',I12)') NQLIST(0)
          CALL STRING_TRIM(OP_STRING(1),IBEG,IEND)
          CALL ASSERT(NQLIST(0).LE.NQM,OP_STRING(1)(IBEG:IEND),ERROR,
     &      *9999)

C Use a function to group points
C new CS June 10 1999, just doing a sphere for now
        ELSE IF(CBBREV(CO,'FUNCTION',2,noco+1,NTCO,N3CO)) THEN
          IF(CBBREV(CO,'SPHERE',1,noco+1,NTCO,N3CO)) THEN
            IF(CBBREV(CO,'RADIUS',1,noco+1,NTCO,N3CO)) THEN
              CALL PARSRL(CO(N3CO+1),1,NTR,radius,ERROR,*9999)
            ELSE
              radius(1)=1.0d0
            ENDIF
            IF(CBBREV(CO,'ORIGIN',2,noco+1,NTCO,N3CO)) THEN
              CALL PARSRL(CO(N3CO+1),3,NTR,origin,ERROR,*9999)
            ELSE
              origin(1)=0.d0
              origin(2)=0.d0
              origin(3)=0.d0
            ENDIF
            NQLIST(0)=0
            na=1 !Temporary
            IF(NJT.EQ.2)THEN
              DO nq=1,NQT
                distance=SQRT(ABS(origin(1)-XQ(1,nq))**2+
     &            ABS(origin(2)-XQ(2,nq))**2)
                IF(distance.LE.radius(1)) THEN
                  NQLIST(0)=NQLIST(0)+1
                  NQLIST(NQLIST(0))=nq
                ENDIF
              ENDDO
            ELSE
              DO nq=1,NQT
                distance=SQRT(ABS(origin(1)-XQ(1,nq))**2+
     &            ABS(origin(2)-XQ(2,nq))**2+
     &            ABS(origin(3)-XQ(3,nq))**2)
                IF(distance.LE.radius(1)) THEN
                  NQLIST(0)=NQLIST(0)+1
                  NQLIST(NQLIST(0))=nq
                ENDIF
              ENDDO
            ENDIF
          ENDIF
          WRITE(OP_STRING(1),
     &      '(''>>Increase NQM >= '',I12)') NQLIST(0)
          CALL STRING_TRIM(OP_STRING(1),IBEG,IEND)
          CALL ASSERT(NQLIST(0).LE.NQM,OP_STRING(1)(IBEG:IEND),ERROR,
     &      *9999)

C Points in a line from first point
        ELSE IF(CBBREV(CO,'LINE',1,noco+1,NTCO,N3CO)) THEN
          na=1 !Temporary
          nq_start=IFROMC(CO(N3CO+1))
          IF(CBBREV(CO,'XIDIRN',3,noco+1,NTCO,N3CO)) THEN
            xi_dirn=IFROMC(CO(N3CO+1))
          ELSE
            xi_dirn=1
          ENDIF
          IF(CBBREV(CO,'POSITIVE',3,noco+1,NTCO,N3CO)) THEN
            DIRN=1
          ELSE IF(CBBREV(CO,'NEGATIVE',3,noco+1,NTCO,N3CO)) THEN
            DIRN=-1
          ELSE
            DIRN=1
          ENDIF
          xi_dirn=xi_dirn*DIRN
          nq=nq_start
          NQLIST(0)=0
          DO WHILE (nq.NE.0)
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=nq
            nq=NXQ(xi_dirn,1,nq,na)
            IF(nq.EQ.nq_start) nq=0
          ENDDO
          WRITE(OP_STRING(1),
     &      '(''>>Increase NQM >= '',I12)') NQLIST(0)
          CALL STRING_TRIM(OP_STRING(1),IBEG,IEND)
          CALL ASSERT(NQLIST(0).LE.NQM,OP_STRING(1)(IBEG:IEND),ERROR,
     &      *9999)

C  List of point numbers/groups

        ELSE IF(CBBREV(CO,'UNGROUPED',2,noco,NTCO,N3CO)) THEN
          DO noelem=1,NELIST(0) !element list
            ne=NELIST(noelem)
            DO nqq=1,NQET(NQS(ne)) !local grid point number in ne
              nq=NQNE(ne,nqq) !global grid point number
              NQLIST(0)=NQLIST(0)+1
              NQLIST(NQLIST(0))=nq
            ENDDO !nqq
          ENDDO !noelem
          ZEROED=.FALSE.
          DO nogrgr=1,NTGRGR !for each grid group
            DO nog=1,NLIGRGR(nogrgr) !number of grid points in group
              nq=ILISTMBR(%VAL(LIGRGR_PTR(nogrgr)),nog)
              IF(INLIST(nq,NQLIST(1),NQLIST(0),NTEMP))THEN
                NQLIST(NTEMP)=0 !remove from the list
                ZEROED=.TRUE.
              ENDIF
            ENDDO !nog
          ENDDO !nogrgr
          IF(ZEROED)THEN
            NTEMP=NQLIST(0)+1
            CALL ILISTRMDUP(NTEMP,NQLIST(0),ERROR,*9999)
          ELSE
            NTEMP=NQLIST(0)
            CALL ILISTRMDUP(NTEMP,NQLIST(1),ERROR,*9999)
          ENDIF
          NQLIST(0)=NTEMP !contains number of unique elements
        ELSE
          CALL PARSE_GRID(NQLIST,noco,NTCO,CO,ERROR,*9999)
        ENDIF

        IF(.NOT.UNSORTED) THEN
C         Sort and remove duplicates from the grid pt list
          CALL ILISTRMDUP(NQLIST(0),NQLIST(1),ERROR,*9999)
        ENDIF

        IF(CBBREV(CO,'AS',1,noco+1,NTCO,N3CO)) THEN
          !Check whether group name already exists
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          CALL CUPPER(CO(N3CO+1)(IBEG:IEND),CHAR)
          N1GRGR=0
          DO nogrgr=1,NTGRGR
            CALL CUPPER(LAGRGR(nogrgr),LABEL)
            CALL STRING_TRIM(LABEL,IBEG2,IEND2)
            IF(CHAR(IBEG:IEND).EQ.LABEL(IBEG2:IEND2)) THEN
              N1GRGR=nogrgr !is existing group label ID
              GO TO 100
            ENDIF
          ENDDO
 100      IF(N1GRGR.EQ.0) THEN !need new group label
            NTGRGR=NTGRGR+1 !is new total #groups
            CALL ASSERT(NTGRGR.LE.GRGR_MAXGRP,
     &        '>>Can''t create any more groups. '//
     &        'Consider reusing an existing group.',ERROR,*9999)
            N1GRGR=NTGRGR
            CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
            LAGRGR(N1GRGR)=CO(N3CO+1)(IBEG:IEND) !is new group label
          ENDIF

        ELSE
          IF(.NOT.FILIO) THEN
            NTGRGR=NTGRGR+1 !is new total
            CALL ASSERT(NTGRGR.LE.GRGR_MAXGRP,
     &        '>>Can''t create any more groups. '//
     &        'Consider reusing an existing group.',ERROR,*9999)
            N1GRGR=NTGRGR
C          CHAR2=CFROMI(N1GRGR,'(I2)')
            WRITE(CHAR2,'(I2)') N1GRGR
            CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
            LAGRGR(N1GRGR)='grid_'//CHAR2(IBEG2:IEND2) !new group label
          ENDIF !.NOT.FILIO
        ENDIF

        IF(.NOT.FILIO) THEN
C KAT 11May99: dynamic groups
C        LIGRGR(0,N1GRGR)=NQLIST(0)
C        CALL ASSERT(NQLIST(0).LE.GRGR_MAXDATA,
C     '    '>>Increase maximum grid group size, GRGR_MAXDATA in '
C     '    //'grou00.cmn',ERROR,*9999)
C        DO nolist=1,NQLIST(0)
C          LIGRGR(nolist,N1GRGR)=NQLIST(nolist)
C        ENDDO
          NLIGRGR(N1GRGR)=NQLIST(0)
          CALL ALLOCATE_MEMORY(NQLIST(0),0,INTTYPE,LIGRGR_PTR(N1GRGR),
     &      MEM_INIT,ERROR,*9999)
          CALL ILIST_COPY(NQLIST(0),NQLIST(1),
     &      %VAL(LIGRGR_PTR(N1GRGR)))

C        IF(DOP) WRITE(*,'(10(1X,I6))')
C     '    (LIGRGR(nolist,N1GRGR),nolist=1,LIGRGR(0,N1GRGR))
        ENDIF !.NOT.FILIO

      ENDIF

      CALL EXITS('GRGRID')
      RETURN
 9999 CALL ERRORS('GRGRID',ERROR)
      CALL EXITS('GRGRID')
      RETURN 1
      END


