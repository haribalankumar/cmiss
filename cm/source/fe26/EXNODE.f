      SUBROUTINE EXNODE(IBT,NBH,NEELEM,NHE,NHP,NHQ,
     '  NKH,NKJ,NP_INTERFACE,NPLIST,NPNODE,NPNY,NQNY,
     '  NRLIST,NRLIST2,NVHP,NVJP,NXLIST,NYNE,NYNP,NYNR,
     '  NYQNR,NW,XP,YP,YQ,YQS,ZA,ZA1,ZP,ZP1,STRING,ERROR,*)

C#### Subroutine: EXNODE
C###  Description:
C###    <HTML>
C###    EXNODE exports node data from finite element data base.
C###    <BR><BR>
C###    <UL>
C###    <LI> NODE_TYPE  is 'geometry' or 'field'
C###    <LI> NODE_NAME  is '1..#nodes' or group name
C###    <LI> NODE_TOTAL is total # nodes in exported list
C###    <LI> NKJ(nj,np) are #derivatives for variables nj of node np
C###    <LI> XP(nk,nv,nj,np) are geometric coordinates & derivatives
C###    </UL>
C### NOTE: for socket connection CONNID2 is defined as a parameter
C### specifying the data transfer socket number
C###  </HTML>

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cmis00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'fsklib.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),
     '  NBH(NHM,NCM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),
     '  NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),NHQ(NRM,NXM),
     '  NKH(NHM,NPM,NCM,0:NRM),NKJ(NJM,NPM),
     '  NP_INTERFACE(0:NPM,0:3),NPLIST(0:NPM),
     '  NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM,NXM),
     '  NQNY(2,NYQM,0:NRCM,NXM),NRLIST(0:NRM),NRLIST2(0:NRM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NVJP(NJM,NPM),
     '  NW(NEM,3,NXM),NXLIST(0:NXM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),
     '  NYQNR(0:NYQM,0:NRCM,NCM,0:NRM,NXM)
      REAL*8 XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM),YQ(NYQM,NIQM,NAM,NXM),
     '  YQS(NIQSM,NQM),ZA(NAM,NHM,NCM,NEM),ZA1(NAM,NHM,NCM,NEM),
     '  ZP(NKM,NVM,NHM,NPM,NCM),ZP1(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)

!     Local Variables
      INTEGER MAX_NH
      PARAMETER (MAX_NH=12)
      INTEGER CLEN,DATA_TYPE,FIELD_BASE_TYPE,IBEG,IBEG1,IBEG2,
     '  IBEG3,IBEG4,
     '  IDUMMY,IEND,IEND1,IEND2,IEND3,IEND4,IFROMC,INTSTR(1024),IY,
     '  na,N3CO,nc,NFIELDT,nh,NHPMAX,NHPMAX_PREV,nhx,
     '  NIQLIST(0:1),NIQSLIST(0:1),NIYLIST(0:16),
     '  nj,njj,NJJ1,NJJ2,
     '  nk,NKHMAXLIST_NP(MAX_NH),NKHMAXLIST_NPFIELD(MAX_NH),NKHMAX_NR,
     '  nohist,no_interface,NOLIST,nonr,nonrlist,nonylist,np,NP_FIELD,
     '  nr,nr1,nrr,NUM_FIELD,NUM_FIELDS,NUM_NP_FIELDS,NUMTIMEDATA,
     '  nv,NVHPMAXLIST_NP(MAX_NH),NVHPMAXLIST_NPFIELD(MAX_NH),
     '  NVHPMAX_NR,nx,nxc,ny,OFFSET,VALUE_INDEX,VERSION_NUM
      REAL*8 FREQUENCY,RELATIVE,RFROMC,TEND,TIME,TOTMAX,TOTMIN,TSTART,
     '  YPMAX(16),YPMIN(16)
      CHARACTER CHAR1*5,CHAR2*5,CHAR3*5,CHAR4*10,FIELD_EX_TYPE*8,
     '  FIELD_NAME*50,FILE*200,FILEFORMAT*6,FILENAME*200,
     '  HISTORYFNAME*200,OUTPUT*11,NODE_NAME*50,NODE_TYPE*50
      !SMAR009 22/12/98 ,ERROR1*100
      LOGICAL ABBREV,ALL_REGIONS,AUTONAME,
     '  CBBREV,DATAFILE,
     '  ENDFILE,FIELDS_CHANGED,FIN,FIRST_NODE,HISTORY,
     '  NJLIST(3),ONEVERSION,SET_FIELD_NAME,SET_FREQUENCY,
     '  YQDATA,YQSDATA,YPDATA

      PARAMETER(DATA_TYPE = 7) !Node list data - must correspond to GUI

!     Functions
      LOGICAL INLIST

      CALL ENTERS('EXNODE',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        WRITE(CHAR2,'(I5)') NPNODE(1,1)
        WRITE(CHAR3,'(I5)') NPNODE(NPNODE(0,1),1)
        CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
        CALL STRING_TRIM(CHAR3,IBEG3,IEND3)
        CHAR4=CHAR2(IBEG2:IEND2)//'..'//CHAR3(IBEG3:IEND3)
        CALL STRING_TRIM(CHAR4,IBEG4,IEND4)

C---------------------------------------------------------------------

C#### Command: FEM export nodes<;FILENAME[default]> [geometry]
C###  Parameter:      <node (GROUP/#s)[all]>
C###    Specify the nodes to export.
C###  Parameter:      <as NAME[0..0]>
C###    Label the file with a character name.
C###  Parameter:      <to (datafile/Motif)[datafile]>
C###    Specify the destination file format (`datafile' for exporting to
C###    CMGUI).
C###  Parameter:      <autoname>
C###    Create a CMGUI link compatible file name automatically.
C###  Parameter:      <offset OFFSET[0]>
C###    Add OFFSET to node numbers.
C###  Parameter:      <region (#s/all)[1]>
C###    Limit to nodes in the specified regions.
C###  Parameter:      <version #[all]>
C###    Limit to only the specified version.
C###  Description:
C###    Exports nodes from CMISS, normally to CMGUI.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME['//FILE00(IBEG1:IEND1)
     '    //']> [geometry]'
        OP_STRING(2)=BLANK(1:15)//'<node (GROUP/#s)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<as NAME['//CHAR4(IBEG4:IEND4)//']>'
        OP_STRING(4)=BLANK(1:15)//'<to (datafile/Motif)[datafile]>'
        OP_STRING(5)=BLANK(1:15)//'<autoname>'
        OP_STRING(6)=BLANK(1:15)//'<offset OFFSET[0]>'
        OP_STRING(7)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(8)=BLANK(1:15)//'<version #[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM export nodes<;FILENAME[default]> field <name FIELDNAME>
C###  Parameter:      <node (GROUP/#s/all)[all]>
C###    Specify the nodes to export.
C###  Parameter:      <as NAME[0..0]>
C###    Label the file with a character name.
C###  Parameter:      <to (datafile/Motif)[datafile]>
C###    Specify the destination file format (`datafile' for exporting to
C###    CMGUI).
C###  Parameter:      <autoname>
C###    Create a CMGUI link compatible file name automatically.
C###  Parameter:      <offset OFFSET[0]>
C###    Add OFFSET to node numbers.
C###  Parameter:      <frequency FREQUENCY[1.0]>
C###    Specify a frequency to convert between stored values and real-time
C###  Parameter:      <iy #[1]>
C###    Specify the index in the YP dependent variable array
C###    for the desired aspect of the solution.
C###  Parameter:      <nc #[1]>
C###    Specify the type of dependent variable.
C###  Parameter:      <region (#s/all)[1]>
C###    Limit to nodes in the specified regions.
C###  Parameter:      <using (fit/solve)[solve]>
C###    Specify the problem type.
C###  Parameter:      <class #[1]>
C###    Specify the class number.
C###  Parameter:      <version #[1]>
C###    Limit to only the specified version.
C###  Description:
C###    Exports nodes with field variables from CMISS,
C###    normally to CMGUI.  If a field name is not specified, then
C###    a default field name will be generated from the problem type.
C###    If a frequency is specified, it is used to convert the stored results
C###    to physical results (e.g. time-step results to real-time).  The
C###    frequency number is the number of time-steps per unit of real time.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME['//FILE00(IBEG1:IEND1)
     '    //']> field <fieldname>'
        OP_STRING(2)=BLANK(1:15)//'<node (GROUP/#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<as NAME['//CHAR4(IBEG4:IEND4)//']>'
        OP_STRING(4)=BLANK(1:15)//'<to (datafile/Motif)[datafile]>'
        OP_STRING(5)=BLANK(1:15)//'<autoname>'
        OP_STRING(6)=BLANK(1:15)//'<offset OFFSET[0]>'
        OP_STRING(7)=BLANK(1:15)//'<frequency FREQUENCY[1.0]>'
        OP_STRING(8)=BLANK(1:15)//'<iy #[1]>'
        OP_STRING(9)=BLANK(1:15)//'<nc #[1]>'
        OP_STRING(10)=BLANK(1:15)//'<region #[1]>'
        OP_STRING(11)=BLANK(1:15)//'<using (fit/solve)[solve]>'
        OP_STRING(12)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(13)=BLANK(1:15)//'<version #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM export nodes<;FILENAME[default]> history <FILENAME>
C###  Parameter:      <(ascii/binary)[ascii]>
C###    Specify whether the file is stored as a binary or
C###    ascii file.
C###  Parameter:      <node (GROUP/#s/all)[all]>
C###    Specify the nodes to export.
C###  Parameter:      <as NAME[0..0]>
C###    Label the file with a character name.
C###  Parameter:      <autoname>
C###    Create a CMGUI link compatible file name automatically.
C###  Parameter:      <offset OFFSET[0]>
C###    Add OFFSET to node numbers.
C###  Parameter:      <tstart TSTART[beginning]>
C###    Limit to times > TSTART.
C###  Parameter:      <tend TEND[end]>
C###    Limit to times <= TEND.
C###  Parameter:      <nc #[1]>
C###    Specify the type of dependent variable.
C###  Parameter:      <iy #[1]>
C###    Specify the index in the YP dependent variable array
C###    for the desired aspect of the solution.
C###  Parameter:      <region (#s/all)[1]>
C###    Limit to nodes in the specified regions.
C###  Parameter:      <intoregion (#s/same)[same]>
C###    Specifies the regions of the exported nodes
C###  Parameter:      <using (fit/solve)[solve]>
C###    Specify the problem type.
C###  Parameter:      <class #[1]>
C###    Specify the class number to use.
C###  Parameter:      <version #[1]>
C###    Limit to only the specified version.
C###  Description:
C###    Exports nodes from CMISS, normally to CMGUI (ie datafile).

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME['//FILE00(IBEG1:IEND1)
     '    //']> history <FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(3)=BLANK(1:15)//'<node GROUP/#s[all]>'
        OP_STRING(4)=BLANK(1:15)//'<as NAME['//CHAR4(IBEG4:IEND4)//']>'
        OP_STRING(5)=BLANK(1:15)//'<autoname>'
        OP_STRING(6)=BLANK(1:15)//'<offset OFFSET[0]>'
        OP_STRING(7)=BLANK(1:15)//'<tstart TIME[beginning]>'
        OP_STRING(8)=BLANK(1:15)//'<tend TIME[end]>'
        OP_STRING(9)=BLANK(1:15)//'<nc #[1]>'
        OP_STRING(10)=BLANK(1:15)//'<iy #[1]>'
        OP_STRING(11)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(12)=BLANK(1:15)//'<intoregion (#s/same)[same]>'
        OP_STRING(13)=BLANK(1:15)//'<using (fit/solve)[solve]>'
        OP_STRING(14)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(15)=BLANK(1:15)//'<version #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe26','doc','EXNODE',ERROR,*9999)
      ELSE

C LKC 6-FEB-1998 Moved initialisation from below
        HISTORY=.FALSE. !unless history specified and nx good
        SET_FIELD_NAME=.FALSE.

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
! AJP 18/8/99
!        nr=NRLIST(1)
        CALL PARSE_NODES(NPNODE,NPLIST,noco,NRLIST,NTCO,CO,
     '    ERROR,*9999)

        IF(CBBREV(CO,'GEOMETRY',2,noco+1,NTCO,N3CO)) THEN
          NODE_TYPE='GEOMETRY'
        ELSE IF(CBBREV(CO,'FIELD',2,noco+1,NTCO,N3CO)) THEN
          NODE_TYPE='FIELD'
C AJPs 31/5/99 Adding field name.
          IF(CBBREV(CO,'NAME',2,noco+1,NTCO,N3CO)) THEN
            CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
            SET_FIELD_NAME=.TRUE.
            FIELD_NAME=CO(N3CO+1)(IBEG1:IEND1)
          ELSE
             FIELD_NAME='-'
          ENDIF
          
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
          CALL ASSERT(NHM.LE.MAX_NH,
     '      '>>Size of NKHMAXLIST and NVHPMAXLIST '
     '      //'arrays must be at least NHM',ERROR,*9999)
        ELSE
          nx=1
          NODE_TYPE='GEOMETRY'
        ENDIF


C GMH 26/12/96 Explicitly initialise NKHMAXLIST_NPFIELD
        DO nh=1,MAX_NH
          NKHMAXLIST_NPFIELD(nh)=0
          NVHPMAXLIST_NPFIELD(nh)=0
          NKHMAXLIST_NP(nh)=0
          NVHPMAXLIST_NP(nh)=0
        ENDDO !nh

        IF(CBBREV(CO,'RADIUS',2,noco+1,NTCO,N3CO)) THEN
          FIELD_EX_TYPE='RADIUS'
        ELSE IF(CBBREV(CO,'MEAN',2,noco+1,NTCO,N3CO)) THEN
          FIELD_EX_TYPE='MEAN' !mean concentration
        ELSE IF(CBBREV(CO,'PRESSURE',6,noco+1,NTCO,N3CO)) THEN
          FIELD_EX_TYPE='PRESSURE' !capillary pressure
        ELSE
          FIELD_EX_TYPE='-'
        ENDIF

C GMH 15/11/95 Allowing user to specify IY
        IF(CBBREV(CO,'IY',2,noco+1,NTCO,N3CO)) THEN
          IY=IFROMC(CO(N3CO+1))
          DO nonylist=0,16
            NIYLIST(nonylist)=0
          ENDDO
          NIYLIST(0)=2
          NIYLIST(1)=1
          NIYLIST(2)=IY
        ELSE
          DO nonylist=0,16
            NIYLIST(nonylist)=0
          ENDDO
          IY=1
          NIYLIST(0)=1
          NIYLIST(1)=1
        ENDIF

C GMH 9/2/97 Allow automatic generation of cmgui link naming convention
        IF(CBBREV(CO,'AUTONAME',4,noco+1,NTCO,N3CO)) THEN
          AUTONAME=.TRUE.
        ELSE
          AUTONAME=.FALSE.
        ENDIF

C news MPN 9-11-95: OFFSET option to add const to all nodes
        IF(CBBREV(CO,'OFFSET',2,noco+1,NTCO,N3CO)) THEN
          OFFSET=IFROMC(CO(N3CO+1))
        ELSE
          OFFSET=0
        ENDIF

C news MPN 25-Jul-95: NAME option for exnode file
        IF(CBBREV(CO,'AS',2,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          NODE_NAME=CO(N3CO+1)(IBEG1:IEND1)
        ELSE
          IF(AUTONAME) THEN
! AJP 18/8/99
!            WRITE(CHAR1,'(I5)') nr
            WRITE(CHAR1,'(I5)') NRLIST(1)
            CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
            NODE_NAME='region_'//CHAR1(IBEG1:IEND1)
          ELSE !autoname
            WRITE(CHAR1,'(I5)') NPLIST(1)+OFFSET
            WRITE(CHAR2,'(I5)') NPLIST(NPLIST(0))+OFFSET
            CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
            CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
            NODE_NAME=CHAR1(IBEG1:IEND1)//'..'//CHAR2(IBEG2:IEND2)
          ENDIF !autoname
        ENDIF

C news AJP 13/9/99 Adding frequency dependency
        IF(CBBREV(CO,'FREQUENCY',2,noco+1,NTCO,N3CO)) THEN
          FREQUENCY=RFROMC(CO(N3CO+1))
          SET_FREQUENCY=.TRUE.
        ELSE
          FREQUENCY=1.0d0
          SET_FREQUENCY=.FALSE.
        ENDIF

        IF(CBBREV(CO,'TO',1,noco+1,NTCO,N3CO)) THEN
          IF(CBBREV(CO,'DATAFILE',2,noco+1,NTCO,N3CO)) THEN
            OUTPUT='DATAFILE'
            DATAFILE=.TRUE.
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

        IF(CBBREV(CO,'NC',2,noco+1,NTCO,N3CO)) THEN
          nc=IFROMC(CO(N3CO+1))
        ELSE
          nc=1
        ENDIF

        IF(CBBREV(CO,'VERSION',2,noco+1,NTCO,N3CO)) THEN
          ONEVERSION=.TRUE.
          VERSION_NUM=IFROMC(CO(N3CO+1))
        ELSE
          ONEVERSION=.FALSE.
        ENDIF

C LKC 6-FEB-1998 Moved initialisation up
C        HISTORY=.FALSE. !unless history specified and nx good
        IF(CBBREV(CO,'HISTORY',2,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          HISTORYFNAME=CO(N3CO+1)(IBEG:IEND)
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
            ENDIF
          ELSE
            CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '        ERROR,*9999)
          ENDIF

          IF(CBBREV(CO,'BINARY',2,noco+1,NTCO,N3CO)) THEN
            FILEFORMAT='BINARY'
          ELSE
            FILEFORMAT='ASCII'
          ENDIF

          IF(CBBREV(CO,'TSTART',2,noco+1,NTCO,N3CO)) THEN
            TSTART=RFROMC(CO(N3CO+1))
          ELSE
            TSTART=-RMAX
          ENDIF

          IF(CBBREV(CO,'TEND',2,noco+1,NTCO,N3CO)) THEN
            TEND=RFROMC(CO(N3CO+1))
          ELSE
            TEND=RMAX
          ENDIF

C          IF(CBBREV(CO,'COUPLED',1,noco+1,NTCO,N3CO)) THEN
C            CALL ASSERT(IS_COUPLED(nx),'>>Define coupled solve first',
C     '        ERROR,*9999)
C            COUPLED=.TRUE.
C          ELSE
C            COUPLED=.FALSE.
C          ENDIF
CC GMH 27/12/96 COUPLED is not used
C          COUPLED=COUPLED

          IF(CBBREV(CO,'INTOREGION',2,noco+1,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),NRM,NRLIST2(0),NRLIST2(1),
     '        ERROR,*9999)
          ELSE
            NRLIST2(0)=NRLIST(0)
            DO nonr=1,NRLIST(0)
              NRLIST2(nonr)=NRLIST(nonr)
            ENDDO
          ENDIF
          HISTORY=.TRUE.

        ELSE
          IF(NODE_TYPE.NE.'FIELD') nx=1
        ENDIF

C*** Exporting a HISTORY file
        IF(HISTORY) THEN
          CALL CHECKF(1,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

C LKC 25-AUG-1999 NJTLIST needs to be initialised
          DO njj1=1,3
            NJLIST(njj1)=.FALSE.
          ENDDO

          NIQLIST(0)=0
          NIQSLIST(0)=0
          YPDATA=.TRUE.
          YQDATA=.FALSE.
          YQSDATA=.FALSE.
          na=1
          CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '      NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,NUMTIMEDATA,
     '      nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,YP(1,1,nx),
     '      YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'READ',FILEFORMAT,HISTORYFNAME,
     '      'OPEN',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,ERROR,*9999)
C AJP 18/8/99
          IF(FILEFORMAT.EQ.'ASCII') THEN !need to find NUMTIMEDATA
            NUMTIMEDATA=0
            ENDFILE=.FALSE.
            DO WHILE(.NOT.ENDFILE)
              CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '          NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,
     '          NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),
     '          NYQNR(0,0,1,0,nx),TIME,YP(1,1,nx),YPMAX,YPMIN,
     '          YQ(1,1,1,nx),YQS,'READ',FILEFORMAT,HISTORYFNAME,
     '          'TIME_DATA',
     '          ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,ERROR,*9999)
              IF(.NOT.ENDFILE) NUMTIMEDATA=NUMTIMEDATA+1
            ENDDO
C*** RESET the input file
            CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '        NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,
     '        NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),
     '        TIME,YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'READ',
     '        FILEFORMAT,HISTORYFNAME,
     '        'RESET',ENDFILE,.TRUE.,YPDATA,YQDATA,
     '        YQSDATA,ERROR,*9999)
          ENDIF
          CALL ASSERT(NUMTIMEDATA.GT.0,
     '      '>>No signal, NUMTIMEDATA.LE.0',ERROR,*9999)
C AJPe 18/8/99
          TOTMAX=-RMAX
          TOTMIN=RMAX
          DO nohist=1,NUMTIMEDATA

            CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '        NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,NUMTIMEDATA,
     '        nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,YP(1,1,nx),
     '        YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'READ',FILEFORMAT,
     '        HISTORYFNAME,'TIME_DATA',ENDFILE,.TRUE.,YPDATA,YQDATA,
     '        YQSDATA,ERROR,*9999)

            IF(TIME.GT.TSTART.AND.TIME.LE.TEND) THEN
              DO nonrlist=1,NRLIST(0)
                nr=NRLIST(nonrlist)
                CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '            NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,
     '            YP(1,IY,nx),ZA,ZP,ERROR,*9999)
              ENDDO !nr
              CALL STRING_TRIM(FILE,IBEG,IEND)
              WRITE(CHAR4,'(I4)') nohist+1000
              CALL STRING_TRIM(CHAR4,IBEG1,IEND1)
              FILENAME=FILE(IBEG:IEND)//'_'//CHAR4(IBEG1:IEND1)
     '          //'.exnode'
              WRITE(OP_STRING,'('' Time = '',D12.4,' //
     '          '''     Creating '',A)')
     '          TIME,FILENAME
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

C LKC 13-JUN-1999 The min and max values calculated in IOHIST correspond
C  to all the data stored. If a subset of nodes are supplied for exporting
C  then the appropriate values are not output. Will calculate the correct
C  values further down.
C              WRITE(OP_STRING,'(''    Min. Value='',D12.4,'
C     '          //''', Max. Value='',D12.4)') YPMIN(1),YPMAX(1)
C              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

                YPMIN(1)=RMAX
                YPMAX(1)=-RMAX
                DO NOLIST=1,NPLIST(0)
                  np=NPLIST(NOLIST)
                  DO nonr=1,NP_INTERFACE(np,0)
                    nr=NP_INTERFACE(np,nonr)
                    DO nhx=1,NHP(np,nr,nx)
                      nh=NH_LOC(nhx,nx)
                      DO nk=1,NKH(nh,np,nc,nr)
                        DO nv=1,NVHP(nh,np,nc,nr)
                          IF (YPMAX(1).LT.ZP(nk,nv,nh,np,nc)) THEN
                            YPMAX(1)=ZP(nk,nv,nh,np,nc)
                          ENDIF
                          IF (YPMIN(1).GT.ZP(nk,nv,nh,np,nc)) THEN
                            YPMIN(1)=ZP(nk,nv,nh,np,nc)
                          ENDIF
                        ENDDO !nv
                      ENDDO ! nk
                    ENDDO !nhx (nh)
                  ENDDO !nonr
                ENDDO !NOLIST (np)
                WRITE(OP_STRING,'(''    Min. Value='',D12.4,'
     '            //''', Max. Value='',D12.4)') YPMIN(1),YPMAX(1)
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)


              CALL OPENF(IFILE,'DISK',FILENAME,'NEW','SEQUEN',
     '          'FORMATTED',132,ERROR,*9999)
C**           write the group name
              CALL STRING_TRIM(NODE_NAME,IBEG,IEND)
              WRITE(IFILE,'( '' Group name: '',A)')
     '          NODE_NAME(IBEG:IEND)
              FIRST_NODE=.TRUE.
              NP_FIELD=NPLIST(1)
              NHPMAX_PREV=0

              DO NOLIST=1,NPLIST(0)
                np=NPLIST(NOLIST)

                CALL ASSERT(NP_INTERFACE(np,0).GT.0,
     '            '>>Node does not belong to any regions',ERROR,*9999)

C LKC 14-JUL-1999 This is wrong
C   This may pick up the wrong region for interface nodes if
C   we are actully trying to export a specific region. Not
C   sure what to do when exporting multiple regions
C   with different dependent variables
C
CC GMH 5/9/96 Get a valid region for this node
C                nr=NP_INTERFACE(np,1)

C*** calculate the maximum NHP's for a np in all regions
                NHPMAX=0
                DO nrr=1,NP_INTERFACE(np,0)
                  nr=NP_INTERFACE(np,nrr)
                  IF( NHPMAX.LT.NHP(NP,nr,nx).AND.
     '              INLIST(nr,NRLIST(1),NRLIST(0),IBEG) ) THEN
                    NHPMAX=NHP(NP,nr,nx)
                  ENDIF
                ENDDO


                IF(NHPMAX.EQ.0) THEN
                  IEND=0
                  CALL APPENDC(IEND,
     '              'WARNING: No Dependent variable for node ',
     '              OP_STRING(1))
                  CALL APPENDI(IEND,np,OP_STRING(1))
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ENDIF

C**             Find max NKH and NVHP across all regions that
C**             contain NP for each nh

C LKC 14-JUL-1999 is not region specific
C                DO nhx=1,NHP(NP,nr,nx)
                DO nhx=1,NHPMAX
                  nh=NH_LOC(nhx,nx)
                  NKHMAX_NR=0
                  NVHPMAX_NR=0
                  DO no_interface=1,NP_INTERFACE(NP,0)
                    nr1=NP_INTERFACE(NP,no_interface)
                    IF(NKHMAX_NR.LT.NKH(nh,NP,nc,nr1))
     '                NKHMAX_NR=NKH(nh,NP,nc,nr1)
                    IF(NVHPMAX_NR.LT.NVHP(nh,NP,nc,nr1))
     '                NVHPMAX_NR=NVHP(nh,NP,nc,nr1)
                  ENDDO !nr1
                  IF(ONEVERSION) THEN
                    NVHPMAX_NR=1
                  ENDIF
                  NKHMAXLIST_NP(nh)=NKHMAX_NR
                  NVHPMAXLIST_NP(nh)=NVHPMAX_NR
                ENDDO !nhx (nh)
C**             check if the fields have changed
                FIELDS_CHANGED=.FALSE.

C LKC 14-JUL-1999 This is not region specific
C                IF(FIRST_NODE.OR.
C     '            (NHP(NP_FIELD,nr,nx).NE.NHP(np,nr,nx))) THEN

                IF(FIRST_NODE.OR.NHPMAX.NE.NHPMAX_PREV) THEN
                  FIELDS_CHANGED=.TRUE.
                ELSE
!AJP 18/8/99
                  nhx=NHP(NP_FIELD,NRLIST(1),nx) !assumes same NH for all regions!!
!                  nhx=NHP(NP_FIELD,nr,nx)
                  nh=NH_LOC(nhx,nx)
                  DO WHILE(nhx.GT.1.AND.
     '              NKHMAXLIST_NP(nh).EQ.NKHMAXLIST_NPFIELD(nh).AND.
     '              NVHPMAXLIST_NP(nh).EQ.NVHPMAXLIST_NPFIELD(nh))
                    nhx=nhx-1
                    nh=NH_LOC(nhx,nx)
                  ENDDO
                  IF(nhx.GT.0.AND.
     '              (NKHMAXLIST_NP(nh).NE.NKHMAXLIST_NPFIELD(nh).OR.
     '              NVHPMAXLIST_NP(nh).NE.NVHPMAXLIST_NPFIELD(nh))) THEN
                    FIELDS_CHANGED=.TRUE.
                  ENDIF
                ENDIF
                IF(FIELDS_CHANGED) THEN
                  NP_FIELD=NP
                  NHPMAX_PREV=NHPMAX
                  VALUE_INDEX=1
!AJP 18/8/99
!                  DO nhx=1,NHP(NP_FIELD,nr,nx)
                  DO nhx=1,NHP(NP_FIELD,NRLIST(1),nx) !assumes same NH for each nr!!!
                    nh=NH_LOC(nhx,nx)
                    NKHMAXLIST_NPFIELD(nh)=NKHMAXLIST_NP(nh)
                    NVHPMAXLIST_NPFIELD(nh)=NVHPMAXLIST_NP(nh)
                    CALL EXPORT_FIELD_HEADING(0,FIELD_BASE_TYPE,
     '                IBT,IFILE,VALUE_INDEX,iy,0,
     '                NBH,nc,NKHMAXLIST_NPFIELD(nh)-1,
     '                0,NFIELDT,nhx,NHE(1,nx),NHP(1,NRLIST(1),nx),
     '                0,NJLIST,NP_FIELD,
     '                NRLIST(1),NVHPMAXLIST_NPFIELD(nh),
     '                NW(1,1,nx),nx,IDUMMY,AUTONAME,DATAFILE,.FALSE.,
     '                FIELD_EX_TYPE,FIELD_NAME,
     '                SET_FIELD_NAME,.FALSE.,ERROR,*9999)
                    VALUE_INDEX=VALUE_INDEX+
     '                NKHMAXLIST_NPFIELD(nh)*NVHPMAXLIST_NPFIELD(nh)
                  ENDDO !nhx (nh)
                ENDIF

C**             write the node
                WRITE(IFILE,'(1X,''Node: '',I12)') NP+OFFSET

                IF(ONEVERSION) THEN
!AJP 18/8/99
!                  DO nhx=1,NHP(np,nr,nx)
                  DO nhx=1,NHP(np,NRLIST(1),nx) !assumes same NH for every nr
                    nh=NH_LOC(nhx,nx)
                    IF(NKHMAXLIST_NP(nh).GT.0.AND.
     '                NVHPMAXLIST_NP(nh).GT.0)
     '                WRITE(IFILE,'(2X,5(1X,E24.16))')
     '                (ZP(nk,VERSION_NUM,nh,np,nc),
     '                nk=1,NKHMAXLIST_NP(nh))
                  ENDDO !nhx (nh)
                ELSE !.not.oneversion
!AJP 18/8/99
!                  DO nhx=1,NHP(np,nr,nx)
                  DO nhx=1,NHP(np,NRLIST(1),nx) !assumes same NH for every nr
                    nh=NH_LOC(nhx,nx)
                    IF(NKHMAXLIST_NP(nh).GT.0.AND.
     '                NVHPMAXLIST_NP(nh).GT.0)
     '                WRITE(IFILE,'(2X,5(1X,E24.16))')
     '                ((ZP(nk,nv,nh,np,nc),nk=1,NKHMAXLIST_NP(nh)),
     '                nv=1,NVHPMAXLIST_NP(nh))
                  ENDDO !nhx (nh)
                ENDIF
                FIRST_NODE=.FALSE.
              ENDDO !nolist (np)

              CALL CLOSEF(IFILE,ERROR,*9999)
              IF(YPMAX(1).GT.TOTMAX) TOTMAX=YPMAX(1)
              IF(YPMIN(1).LT.TOTMIN) TOTMIN=YPMIN(1)
            ENDIF
          ENDDO !nohist
          WRITE(OP_STRING,'('' Total Min. ='',D12.4,'', Total Max. ='','
     '      //'D12.4)') TOTMIN,TOTMAX
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

C*** Close history file
          CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '      NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,NUMTIMEDATA,
     '      nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,YP(1,1,nx),
     '      YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'CLOSE',FILEFORMAT,FILE,' ',
     '      ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,ERROR,*9999)

C*** Exporting DATAFILE or MOTIF
        ELSE
          IF((OUTPUT(1:8).EQ.'DATAFILE').OR.(OUTPUT(1:5).EQ.'MOTIF'))
     '      THEN
            CALL CHECKF(1,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
            IF(NPLIST(0).GT.0) THEN
              IF(DATAFILE) THEN
                CALL STRING_TRIM(FILE,IBEG,IEND)
                CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.exnode',
     '            'NEW','SEQUEN','FORMATTED',132,ERROR,*9999)
C**             write the group name
                CALL STRING_TRIM(NODE_NAME,IBEG,IEND)
                WRITE(IFILE,'( '' Group name: '',A)')
     '            NODE_NAME(IBEG:IEND)
              ELSE
                IF(USE_SOCKET) THEN
C**               send the data type identifier
                  IF(FSKWRITE(DATA_TYPE,SK_LONG_INT,1,CONNID2).EQ.-1)
     '              GOTO 9999
C**               send the group name
                  CALL STRING_TRIM(NODE_NAME,IBEG,IEND)
                  CLEN=FSKLEN(NODE_NAME(IBEG:IEND))
                  CALL FSKF2C(NODE_NAME(IBEG:IEND),CLEN,INTSTR)
                  IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
     '              GOTO 9999
                  IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
     '              GOTO 9999
                ELSE
                  ERROR='>>Can only export to MOTIF if sockets are '
     '              //'being used'
                  GOTO 9999
                ENDIF
              ENDIF
              FIRST_NODE=.TRUE.
              NP_FIELD=NPLIST(1)
              NUM_NP_FIELDS=0
C*** Exporting Geometry
              IF(NODE_TYPE(1:8).EQ.'GEOMETRY') THEN
                DO nolist=1,NPLIST(0)
                  NP=NPLIST(nolist)
C*** Check if the fields have changed
                  FIELDS_CHANGED=.FALSE.

C LKC 5-FEB-1999 Appears to be unnecessary and if no equations
C  defined then NHP is uninitialised
C                  IF(FIRST_NODE.OR.(NHP(NP_FIELD,nr,nx).NE.
C     '              NHP(np,nr,nx))) THEN
                  IF(FIRST_NODE) THEN
                    FIELDS_CHANGED=.TRUE.
                  ELSE
                    NUM_FIELDS=0
                    DO njj1=1,3 !geometry, fibre, general
                      NUM_FIELD=0
                      DO nonr=1,NP_INTERFACE(NP,0)
                        nr=NP_INTERFACE(NP,nonr)
                        IF(NJ_LOC(njj1,0,nr).GT.NUM_FIELD) THEN
                          NUM_FIELD=NJ_LOC(njj1,0,nr)
                          nrr=nr
                        ENDIF
                      ENDDO
                      njj2=NJ_LOC(njj1,0,nrr)
                      NUM_FIELDS=NUM_FIELDS+NUM_FIELD
                      IF(njj2.GT.0) THEN
                        nj=NJ_LOC(njj1,njj2,nrr)
                        FIN=.FALSE.
                        DO WHILE(.NOT.FIN)
                          IF(njj2.LE.1) THEN
                            FIN=.TRUE.
                          ELSE
                            IF(NKJ(nj,NP_FIELD).NE.NKJ(nj,NP)) THEN
                              FIN=.TRUE.
                            ELSE
                              IF(NVJP(nj,NP_FIELD).NE.NVJP(nj,NP))
     '                          THEN
                                FIN=.TRUE.
                              ELSE
                                njj2=njj2-1
                                nj=NJ_LOC(njj1,njj2,nrr)
                              ENDIF
                            ENDIF
                          ENDIF
                        ENDDO
                        IF(njj2.GT.0.AND.
     '                    (NKJ(nj,NP_FIELD).NE.NKJ(nj,NP).OR.
     '                    NVJP(nj,NP_FIELD).NE.NVJP(nj,NP))) THEN
                          FIELDS_CHANGED=.TRUE.
                        ENDIF
                      ENDIF
                    ENDDO !njj1
                    IF(NUM_FIELDS.NE.NUM_NP_FIELDS) THEN
                      FIELDS_CHANGED=.TRUE.
                    ENDIF
                  ENDIF !FIRST_NODE
                  IF(FIELDS_CHANGED) THEN
                    IF(.NOT.DATAFILE) THEN
                      IF(FSKWRITE(ICHAR('#'),SK_CHAR,1,CONNID2).EQ.-1)
     '                  GOTO 9999
                    ENDIF

C LKC 15-NOV-97 - new section to calc the number of field for a node

C*** Creating  a NJLIST
C    ! contains the types of existing fields at a given node
                    DO njj1=1,3 !geometry, fibre, general
                      NJLIST(njj1)=.FALSE.
                    ENDDO

                    DO nonr=1,NP_INTERFACE(NP,0) !interface
C 25-JUN-1999  LKC adding nrr
C                      DO nr=1,NRLIST(0) !all reg. exporting
C                      DO nrr=1,NRLIST(0) !all reg. exporting
C                        IF(NP_INTERFACE(NP,nonr).EQ.NRLIST(nr)) THEN
C LKC 1-JUL-1999 Another modification

                      DO nrr=1,NRLIST(0) !all reg. exporting
                        nr=NRLIST(nrr)
                        IF(NP_INTERFACE(NP,nonr).EQ.NRLIST(nrr)) THEN
                          DO njj1=1,3 !each field type
C KAT 25Jul99: Needing to export something even if no versions exist.
C KAT 21May99: Checking versions exist
C                            DO njj2=1,NJ_LOC(njj1,0,nr)
C                              nj=NJ_LOC(njj1,njj2,nr)
C                              IF(NVJP(nj,NP).NE.0) THEN
C                                NJLIST(njj1)=.TRUE.
C                              ENDIF
C                            ENDDO
                            IF(NJ_LOC(njj1,0,NRLIST(nrr)).NE.0)
     '                        NJLIST(njj1)=.TRUE.
                          ENDDO
                        ENDIF
                      ENDDO
                    ENDDO

C*** Write the field information
                    NP_FIELD=NP
                    VALUE_INDEX=1
                    NUM_NP_FIELDS=0
                    DO njj1=1,3 !geometry, fibre, general
                      NUM_FIELD=0
                      DO nonr=1,NP_INTERFACE(NP_FIELD,0)
                        nr=NP_INTERFACE(NP_FIELD,nonr)
                        IF(NJ_LOC(njj1,0,nr).GT.NUM_FIELD) THEN
                          NUM_FIELD=NJ_LOC(njj1,0,nr)
                          nrr=nr
                        ENDIF
                      ENDDO
                      IF(NUM_FIELD.GT.0) THEN
                        DO njj2=1,NJ_LOC(njj1,0,nrr)
                          nj=NJ_LOC(njj1,njj2,nrr)
C KAT 25Jul99: Needing to export something even if no versions exist.
C                          IF(NVJP(nj,NP_FIELD).NE.0) THEN
                          IF(ONEVERSION) THEN
                            CALL EXPORT_FIELD_HEADING(0,
     '                        FIELD_BASE_TYPE,IBT,IFILE,
     '                        VALUE_INDEX,iy,0,NBH,
     '                        nc,MAX(NKJ(nj,NP_FIELD)-1,0),0,NFIELDT,0,
     '                        NHE(1,nx),NHP(1,nrr,nx),nj,
     '                        NJLIST,NP_FIELD,nrr,
     '                        1,NW(1,1,nx),nx,IDUMMY,AUTONAME,DATAFILE,
     '                        .FALSE.,FIELD_EX_TYPE,FIELD_NAME,
     '                        SET_FIELD_NAME,.FALSE.,ERROR,*9999)
                            VALUE_INDEX=VALUE_INDEX+
     '                        MAX(NKJ(nj,NP_FIELD),1)
                          ELSE !.not.oneversion
                            CALL EXPORT_FIELD_HEADING(0,
     '                        FIELD_BASE_TYPE,IBT,IFILE,
     '                        VALUE_INDEX,iy,0,NBH,
     '                        nc,MAX(NKJ(nj,NP_FIELD)-1,0),0,NFIELDT,0,
     '                        NHE(1,nx),NHP(1,nrr,nx),nj,
     '                        NJLIST,NP_FIELD,nrr,NVJP(nj,NP_FIELD),
     '                        NW(1,1,nx),nx,IDUMMY,AUTONAME,
     '                        DATAFILE,.FALSE.,FIELD_EX_TYPE,
     '                        FIELD_NAME,SET_FIELD_NAME,.FALSE.,
     '                        ERROR,*9999)
                            VALUE_INDEX=VALUE_INDEX+
     '                        MAX(NKJ(nj,NP_FIELD)*NVJP(nj,NP_FIELD),1)
                          ENDIF
                          NUM_NP_FIELDS=NUM_NP_FIELDS+1
C                          ENDIF !NVJP!=0
                         ENDDO !njj2
C                        NUM_NP_FIELDS=NUM_NP_FIELDS+NUM_FIELD
                      ENDIF !NUMFIELD>0
                    ENDDO !njj1
                  ENDIF !FIELDS_CHANGED
C***              write the node
                  IF(DATAFILE) THEN
                    WRITE(IFILE,'(1X,''Node: '',I12)') NP+OFFSET
                  ELSE
                    IF(FSKWRITE(ICHAR('N'),SK_CHAR,1,CONNID2).EQ.-1)
     '                GOTO 9999
                    IF(FSKWRITE(NP+OFFSET,SK_LONG_INT,1,CONNID2).EQ.-1)
     '                GOTO 9999
                  ENDIF

                  DO njj1=1,3 !geometry, fibre, general
                    NUM_FIELD=0
                    DO nonr=1,NP_INTERFACE(NP,0)
                      nr=NP_INTERFACE(NP,nonr)
                      IF(NJ_LOC(njj1,0,nr).GT.NUM_FIELD) THEN
                        NUM_FIELD=NJ_LOC(njj1,0,nr)
                        nrr=nr
                      ENDIF
                    ENDDO
C KAT 25Jul99: Check fields exist for this nj type.
                    IF(NUM_FIELD.GT.0) THEN
                      DO njj2=1,NJ_LOC(njj1,0,nrr)
                        nj=NJ_LOC(njj1,njj2,nrr)
                        IF(DATAFILE) THEN
C KAT 25Jul99: Needing to export something even if no versions exist.
                          IF(NKJ(nj,NP).GT.0.AND.NVJP(nj,NP).GT.0) THEN
                            IF(ONEVERSION) THEN
                              WRITE(IFILE,'(2X,5(1X,E24.16))')
     '                          (XP(nk,VERSION_NUM,nj,NP),
     '                          nk=1,NKJ(nj,NP))
                            ELSE !.not.oneversion
                              WRITE(IFILE,'(2X,5(1X,E24.16))')
     '                          ((XP(nk,nv,nj,NP),nk=1,NKJ(nj,NP)),
     '                          nv=1,NVJP(nj,NP))
                            ENDIF
                          ELSE
                            WRITE(IFILE,'(3X,I1)') 0
                          ENDIF
                        ELSE
                          IF(ONEVERSION) THEN
                            IF(FSKWRITE(XP(1,VERSION_NUM,nj,NP),
     '                        SK_DOUBLE_FLOAT,NKJ(nj,NP),CONNID2).EQ.-1)
     '                        GOTO 9999
                          ELSE !.not.oneversion
                            DO nv=1,NVJP(nj,NP)
                              IF(FSKWRITE(XP(1,nv,nj,NP),
     '                          SK_DOUBLE_FLOAT,NKJ(nj,NP),CONNID2)
     '                          .EQ.-1) GOTO 9999
                            ENDDO !nv
                          ENDIF
                        ENDIF
                      ENDDO !njj2
                    ENDIF
                  ENDDO !njj1
                  FIRST_NODE=.FALSE.
                ENDDO !nolist (np)

C*** Exporting a Field
              ELSE IF(NODE_TYPE(1:5).EQ.'FIELD') THEN


C LKC 26-JAN-2007 Modify to allow exporting fields for multiple regions                
C                nr=NRLIST(1) !temporary - needs fixing AJP 18/8/99

                DO nrr=1,NRLIST(0)
                  nr=NRLIST(nrr)

                  DO njj1=1,3 !geometry, fibre, general
                    NJLIST(njj1)=.FALSE.
                    IF(NJ_LOC(njj1,0,nr).NE.0) NJLIST(njj1)=.TRUE.
                  ENDDO
                
                  CALL YPZP(IY,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '              NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,
     '              YP(1,1,nx),ZA,ZP,ERROR,*9999)

                
                  IF((ITYP2(nr,nx).EQ.2.OR. !Finite elasticity
     '              ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) !const vol fluid
     '              .AND.NJ_LOC(NJL_FIBR,0,nr).GT.0) THEN !fibres defined
                    CALL YPZP(5,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '                NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,
     &                NYNP,YP(1,1,nx),ZA1,ZP1,ERROR,*9999)
                  ENDIF
                  
                ENDDO !nrr
                
                DO NOLIST=1,NPLIST(0)
                  NP=NPLIST(NOLIST)
C GMH 5/9/96 Get a valid region for this node
                  CALL ASSERT(NP_INTERFACE(np,0).GT.0,
     '              '>>Node does not belong to any regions',ERROR,*9999)
                  nr=NP_INTERFACE(np,1) !assumes same NH for each NR
C**               Find max NKH and NVHP across all regions that
C**               contain NP for each nh
                  DO nhx=1,NHP(NP,nr,nx)
                    nh=NH_LOC(nhx,nx)
                    NKHMAX_NR=0
                    NVHPMAX_NR=0
                    DO no_interface=1,NP_INTERFACE(NP,0)
                      nr1=NP_INTERFACE(NP,no_interface)
                      IF(NKHMAX_NR.LT.NKH(nh,NP,nc,nr1)) THEN
C LKC 24-NOV-2005 add check if cross derivatives set to 0 in ipequa
C file and not in the ipbase file - then export 1 less varible.
C I'm guess this is assuming that KTYP93 can only be 1 for a bicubic field.                        
C     '                  NKHMAX_NR=NKH(nh,NP,nc,nr1)
                        IF(NKH(nh,NP,nc,nr1).EQ.1) THEN
                          NKHMAX_NR=1 !top/bottom of the mesh
                        ELSE
                          NKHMAX_NR=NKH(nh,NP,nc,nr1)-KTYP93(nc,nr)
                        ENDIF
                      ENDIF
                      IF(NVHPMAX_NR.LT.NVHP(nh,NP,nc,nr1))
     '                  NVHPMAX_NR=NVHP(nh,NP,nc,nr1)
                    ENDDO !nr1
                    IF(ONEVERSION) THEN
                      NVHPMAX_NR=1
                    ENDIF
                    NKHMAXLIST_NP(nh)=NKHMAX_NR
                    NVHPMAXLIST_NP(nh)=NVHPMAX_NR
                  ENDDO !nhx (nh)
C**               check if the fields have changed
                  FIELDS_CHANGED=.FALSE.
                  IF(FIRST_NODE.OR.
     '              (NHP(NP_FIELD,nr,nx).NE.NHP(np,nr,nx))) THEN
                    FIELDS_CHANGED=.TRUE.
                  ELSE
                    nhx=NHP(NP_FIELD,nr,nx)
                    nh=NH_LOC(nhx,nx)
                    DO WHILE(nhx.GT.1.AND.
     '                NKHMAXLIST_NP(nh).EQ.NKHMAXLIST_NPFIELD(nh).AND.
     '                NVHPMAXLIST_NP(nh).EQ.NVHPMAXLIST_NPFIELD(nh))
                      nhx=nhx-1
                      nh=NH_LOC(nhx,nx)
                    ENDDO
                    IF(nhx.GT.0.AND.
     '                (NKHMAXLIST_NP(nh).NE.NKHMAXLIST_NPFIELD(nh).OR.
     '                NVHPMAXLIST_NP(nh).NE.NVHPMAXLIST_NPFIELD(nh)))
     '                THEN
                      FIELDS_CHANGED=.TRUE.
                    ENDIF
                    IF(ITYP2(nr,nx).EQ.2.OR.
     '                ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN
C                     Check if deformed fibre fields have changed for
C                     Finite elasticity or fluid mechanics/
C                     const vol constraint
                      njj2=NJ_LOC(NJL_FIBR,0,nr)
                      IF (njj2.GE.1) THEN
                        nj=NJ_LOC(NJL_FIBR,njj2,nr)
                        DO WHILE (njj2.GT.1.AND.
     '                    NKJ(nj,NP_FIELD).EQ.NKJ(nj,NP).AND.
     '                    NVJP(nj,NP_FIELD).EQ.NVJP(nj,NP))
                          njj2=njj2-1
                          IF (njj2.NE.0) nj=NJ_LOC(NJL_FIBR,njj2,nr)
                        ENDDO
                      ENDIF
                      IF(njj2.GT.0.AND.
     '                  (NKJ(nj,NP_FIELD).NE.NKJ(nj,NP).OR.
     '                  NVJP(nj,NP_FIELD).NE.NVJP(nj,NP))) THEN
                        FIELDS_CHANGED=.TRUE.
                      ENDIF
                      IF(njj2.EQ.0) THEN
                        FIELDS_CHANGED=.TRUE.
                        NJLIST(2)=.FALSE.
                      ENDIF
                    ENDIF !Finite elasticity
                  ENDIF !FIRST_NODE...
                  IF(FIELDS_CHANGED) THEN
                    IF(.NOT.DATAFILE) THEN
                      IF(FSKWRITE(ICHAR('#'),SK_CHAR,1,CONNID2).EQ.-1)
     '                  GOTO 9999
                    ENDIF
                    NP_FIELD=NP
                    VALUE_INDEX=1
                    DO nhx=1,NHP(NP_FIELD,nr,nx)
                      nh=NH_LOC(nhx,nx)
                      NKHMAXLIST_NPFIELD(nh)=NKHMAXLIST_NP(nh)
                      NVHPMAXLIST_NPFIELD(nh)=NVHPMAXLIST_NP(nh)
                      CALL EXPORT_FIELD_HEADING(0,FIELD_BASE_TYPE,
     '                  IBT,IFILE,VALUE_INDEX,
     '                  iy,0,NBH,nc,NKHMAXLIST_NPFIELD(nh)-1,
     '                  0,NFIELDT,nhx,NHE(1,nx),NHP(1,nr,nx),0,
     '                  NJLIST,
     '                  NP_FIELD,nr,
     '                  NVHPMAXLIST_NPFIELD(nh),
     '                  NW(1,1,nx),nx,IDUMMY,AUTONAME,
     '                  DATAFILE,.FALSE.,FIELD_EX_TYPE,
     '                  FIELD_NAME,SET_FIELD_NAME,.FALSE.,ERROR,*9999)
                      VALUE_INDEX=VALUE_INDEX+
     '                  NKHMAXLIST_NPFIELD(nh)*NVHPMAXLIST_NPFIELD(nh)
                    ENDDO !nhx (nh)
                    IF(ITYP2(nr,nx).EQ.2.OR.
     '                ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN
C                     Finite elasticity or fluid mechanics/
C                     const vol constraint
                      DO njj=1,NJ_LOC(NJL_FIBR,0,nr)
                        nj=NJ_LOC(NJL_FIBR,njj,nr)
                        nh=njj
                        CALL EXPORT_FIELD_HEADING(0,FIELD_BASE_TYPE,
     '                    IBT,IFILE,
     '                    VALUE_INDEX,iy,0,NBH,
     '                    nc,NKJ(nj,NP_FIELD)-1,
     '                    0,NFIELDT,nh,NHE(1,nx),NHP(1,nr,nx),0,
     '                    NJLIST,NP_FIELD,nr,
     '                    NVJP(nj,NP_FIELD),
     '                    NW(1,1,nx),nx,IDUMMY,AUTONAME,
     '                    DATAFILE,.TRUE.,FIELD_EX_TYPE,
     '                    FIELD_NAME,SET_FIELD_NAME,.FALSE.,ERROR,*9999)
                        VALUE_INDEX=VALUE_INDEX+
     '                    NKJ(nj,NP_FIELD)*NVJP(nj,NP_FIELD)
                      ENDDO !njj
                    ENDIF !Finite elasticity
                  ENDIF !FIELDS_CHANGED
C**               write the node
                  IF(DATAFILE) THEN
                    WRITE(IFILE,'(1X,''Node: '',I12)') NP+OFFSET
                  ELSE
                    IF(FSKWRITE(ICHAR('N'),SK_CHAR,1,CONNID2).EQ.-1)
     '                GOTO 9999
                    IF(FSKWRITE(NP+OFFSET,SK_LONG_INT,1,CONNID2).EQ.-1)
     '                GOTO 9999
                  ENDIF
                  DO nhx=1,NHP(np,nr,nx)
                    nh=NH_LOC(nhx,nx)
                    IF(DATAFILE) THEN
! news AJP 13/9/99 Using Frequency
                      IF(SET_FREQUENCY) THEN
                        IF(FREQUENCY.GT.RDELTA) THEN
                          DO nk=1,NKHMAXLIST_NP(nh)
                            DO nv=1,NVHPMAXLIST_NP(nh)
                              ZP(nk,nv,nh,np,nc)=
     '                          ZP(nk,nv,nh,np,nc)/FREQUENCY
                            ENDDO !nv
                          ENDDO !nk
                        ELSE
                          CALL ASSERT(1.LT.0,
     '                      '>>Frequency zero or negative??',
     '                      ERROR,*9999)
                        ENDIF !FREQUENCY > 0
                      ENDIF !SET_FREQUENCY
                      IF(NKHMAXLIST_NP(nh).GT.0.AND.
     '                  NVHPMAXLIST_NP(nh).GT.0) THEN
                        IF(ONEVERSION) THEN
                          IF(FIELD_EX_TYPE.EQ.'RADIUS')THEN
                            nj=NJ_LOC(3,1,nr)
                            WRITE(IFILE,'(2X,5(1X,E24.16))')
     '                        (XP(nk,1,nj,np),
     '                        nk=1,NKHMAXLIST_NP(nh))
                          ELSE IF(FIELD_EX_TYPE.EQ.'MEAN')THEN
                            nk=NKHMAXLIST_NP(1)
                            nv=1
                            ny=NYNP(nk,nv,nh,np,1,1,nr)
                            WRITE(IFILE,'(2X,5(1X,E24.16))')YP(ny,8,1)
                          ELSE IF(FIELD_EX_TYPE.EQ.'PRESSURE')THEN
                            WRITE(IFILE,'(2X,5(1X,E24.16))')
     '                        ZP(nk,VERSION_NUM,nh,NP,nc)
                          ELSE
                            WRITE(IFILE,'(2X,5(1X,E24.16))')
     '                        (ZP(nk,VERSION_NUM,nh,NP,nc),
     '                        nk=1,NKHMAXLIST_NP(nh))
                          ENDIF
                        ELSE !.not.oneversion
                          IF(FIELD_EX_TYPE.EQ.'RADIUS')THEN
                            nj=NJ_LOC(3,1,nr)
                            WRITE(IFILE,'(2X,5(1X,E24.16))')
     '                        (XP(nk,1,nj,np),
     '                        nk=1,NKHMAXLIST_NP(nh))
                          ELSE IF(FIELD_EX_TYPE.EQ.'MEAN')THEN
                            nk=NKHMAXLIST_NP(1)
                            nv=1
                            ny=NYNP(nk,nv,nh,np,1,1,nr)
                            WRITE(IFILE,'(2X,5(1X,E24.16))')YP(ny,8,1)
                          ELSE IF(FIELD_EX_TYPE.EQ.'RH')THEN
                            IF(nhx.EQ.1)THEN
                              RELATIVE=ZP(1,1,nh+1,NP,nc)/
     '                          0.4055384d0*DEXP(4.97d3/
     '                          ZP(1,1,nh,NP,nc))
                              WRITE(IFILE,'(2X,5(1X,E24.16))')
     '                          RELATIVE
                            ENDIF
                          ELSE IF(FIELD_EX_TYPE.EQ.'PRESSURE')THEN
                            WRITE(IFILE,'(2X,5(1X,E24.16))')
     '                        ZP(1,1,nh,NP,nc)
                          ELSE
                            WRITE(IFILE,'(2X,5(1X,E24.16))')
     '                        ((ZP(nk,nv,nh,NP,nc),
     '                        nk=1,NKHMAXLIST_NP(nh)),
     '                        nv=1,NVHPMAXLIST_NP(nh))
                          ENDIF
                        ENDIF
                      ENDIF
                    ELSE
                      IF(ONEVERSION) THEN
                        IF(FSKWRITE(ZP(1,VERSION_NUM,nh,NP,nc),
     '                    SK_DOUBLE_FLOAT,NKHMAXLIST_NP(nh),
     '                    CONNID2).EQ.-1) GOTO 9999
                      ELSE
                        DO nv=1,NVHPMAXLIST_NP(nh)
                          IF(FSKWRITE(ZP(1,nv,nh,NP,nc),SK_DOUBLE_FLOAT,
     '                      NKHMAXLIST_NP(nh),CONNID2).EQ.-1) GOTO 9999
                        ENDDO !nv
                      ENDIF
                    ENDIF
                  ENDDO
                  IF(ITYP2(nr,nx).EQ.2.OR.
     '              ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN
C                   Finite elasticity or fluid mechanics/const vol
                    DO njj=1,NJ_LOC(NJL_FIBR,0,nr)
                      nj=NJ_LOC(NJL_FIBR,njj,nr)
                      nh=NH_LOC(njj,nx)
                      IF(DATAFILE) THEN
                        IF(ONEVERSION) THEN
                          WRITE(IFILE,'(2X,5(1X,E24.16))')
     '                      (ZP1(nk,VERSION_NUM,nh,NP,1),nk=1,NKJ(nj,
     '                      NP))
                        ELSE !.not.oneversion
                          WRITE(IFILE,'(2X,5(1X,E24.16))')
     '                      ((ZP1(nk,nv,nh,NP,1),nk=1,NKJ(nj,NP)),
     '                      nv=1,NVJP(nj,NP))
                        ENDIF
                      ELSE
                        IF(ONEVERSION) THEN
                          IF(FSKWRITE(ZP1(1,VERSION_NUM,nh,NP,1),
     '                      SK_DOUBLE_FLOAT,NKJ(nj,NP),CONNID2).EQ.-1)
     '                      GOTO 9999
                        ELSE !.not.oneversion
                          DO nv=1,NVJP(nj,NP)
                            IF(FSKWRITE(ZP1(1,nv,nh,NP,1),
     '                        SK_DOUBLE_FLOAT,NKJ(nj,NP),CONNID2).EQ.-1)
     '                        GOTO 9999
                          ENDDO !nv
                        ENDIF
                      ENDIF
                    ENDDO !njj
                  ENDIF !Finite elasticity
                  FIRST_NODE=.FALSE.
                ENDDO !nolist (np)
              ENDIF !NODE_TYPE
              IF(DATAFILE) THEN
                CALL CLOSEF(IFILE,ERROR,*9999)
              ELSE
                IF(FSKWRITE(0,SK_CHAR,1,CONNID2).EQ.-1) GOTO 9999
              ENDIF
            ENDIF

          ELSE IF(OUTPUT(1:8).EQ.'EXPLORER') THEN
C LKC 25-AUG-1999 added error message
            ERROR='Export to EXPLORER not implemented'
            GOTO 9999
          ENDIF

        ENDIF

      ENDIF

      CALL EXITS('EXNODE')
      RETURN
 9999 IF(HISTORY) THEN
        CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '    NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,NUMTIMEDATA,
     '    nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,YP(1,1,nx),
     '    YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'CLOSE',FILEFORMAT,FILE,' ',
     '    ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,ERROR,*1111)
      ENDIF
 1111 CALL ERRORS('EXNODE',ERROR)
      CALL EXITS('EXNODE')
      RETURN 1
      END


