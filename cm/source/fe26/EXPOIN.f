      SUBROUTINE EXPOIN(LIST,NBH,NBJ,NHE,NKHE,NKJE,NPF,NPNE,
     '  NPNODE,NRE,NVHE,NVJE,NW,NXLIST,NYNQ,AQ,CE,CG,CGE,
     '  CP,CURVCORRECT,PG,
     '  SE,XA,XE,XG,XP,XQ,YQ,ZA,ZE,ZG,ZP,
     '  YQS,
     '  CELL_ICQS_VALUE,CELL_ICQS_SPATIAL,CELL_ICQS_NAMES,
     '  CELL_RCQS_VALUE,CELL_RCQS_SPATIAL,CELL_RCQS_NAMES,
     '  CELL_YQS_VALUE ,CELL_YQS_SPATIAL ,CELL_YQS_NAMES ,
     '  ICQS_SPATIAL,IICQS_SPATIAL,IRCQS_SPATIAL,RCQS_SPATIAL,
     '  STRING,ERROR,*)

C#### Subroutine: EXPOIN
C###  Description:
C###    EXPOIN exports point data from finite element data base.

C**** NOTE: for socket connection CONNID2 is defined as a parameter
C**** specifying the data transfer socket number
C**** 22/6/94 Now conforms to .EXNODE file structure (GBS)
C?????DB.  This needs re-thinking.  Maybe we need several numbers for
C?????  each cmgui node (like the element, face and line numbers for
C?????  the cmgui elements
C**** 22/04/96 GBS Added export of potentials and time at gridpoints
C *** DPN 13 September 1999 - added cellular materials at gridpoints

      IMPLICIT none
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
      INCLUDE 'grou00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'loc00.cmn'
      INCLUDE 'maqloc00.inc'
      INCLUDE 'nqloc00.inc'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER LIST(0:NLISTM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NHE(NEM,NXM),NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),NRE(NEM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM),
     '  NXLIST(0:NXM),NYNQ(NHM,NQM,0:NRCM,NXM),
     '  CELL_ICQS_SPATIAL(NQIM,NQVM),CELL_ICQS_VALUE(NQIM,NQVM),
     '  CELL_RCQS_SPATIAL(NQRM,NQVM),CELL_YQS_SPATIAL(NIQSM,NQVM),
     '  ICQS_SPATIAL(NQISVM,NQM),
     '  IICQS_SPATIAL(0:NQISVM,NQVM),
     '  IRCQS_SPATIAL(0:NQRSVM,NQVM)
      REAL*8 AQ(NMAQM,NQM),CE(NMM,NEM,NXM),CG(NMM,NGM),
     '  CGE(NMM,NGM,NEM,NXM),CP(NMM,NPM,NXM),
     '  CURVCORRECT(2,2,NNM,NEM),PG(NSM,NUM,NGM,NBM),SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),
     '  XP(NKM,NVM,NJM,NPM),XQ(NJM,NQM),YQ(NYQM,NIQM,NAM,NXM),
     '  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZG(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM),YQS(NIQSM,NQM),
     '  CELL_RCQS_VALUE(NQRM,NQVM),CELL_YQS_VALUE(NIQSM,NQVM),
     '  RCQS_SPATIAL(NQRSVM,NQM)
      CHARACTER ERROR*(*),
     '  CELL_ICQS_NAMES(NQIM,NQVM)*(CELL_NAME_LENGTH),
     '  CELL_RCQS_NAMES(NQRM,NQVM)*(CELL_NAME_LENGTH),
     '  CELL_YQS_NAMES(NIQSM,NQVM)*(CELL_NAME_LENGTH),
     '  STRING*(MXCH)
!     Local Variables
      INTEGER atime,CLEN,DATA_TYPE,IBEG,IBEG1,IBEG2,IBEG3,IBEG4,IEND,
     '  IEND1,IEND2,IEND3,IEND4,IFROMC,ILISTMBR,index,
     '  INTSTR(1024),ISOCKET(3),N1GRGA,N3CO,nb,
     '  nc,ne,ne_num_points,ng,niqV,nj,nogrga,
     '  NOLIST,nong,NONODE,nq,nr,NUMFIELDS,
     '  NUM_FIELD_LIST,nx,nxc,ny_p,ny_r,ny_v,OFFSET,VAL_INDEX,nqsv,j,k
      REAL*8 ACT_FIELD,DXINU(3,3),EG(3,3),FIELD_LIST(99),GL(3,3),
     '  GU(3,3),PATH_FIELD,POT_FIELD,PRESS_FIELD,RADIUS_FIELD,RG,
     '  RT(3,3),TG(3,3),VEL_FIELD
      CHARACTER CHAR*30,CHAR1*5,CHAR2*5,CHAR3*5,
     '  COMPONENT_NAME*80,FIELD_NAME*80,
     '  FILE*100,LABEL*30,OUTPUT*11,POINT_NAME*50,POINT_TYPE*50
      LOGICAL ACTIVATION_TIME,CBBREV,CELL_FIELD,DATAFILE,PATH,GROUPED,
     '  MATERIAL,POTENTIAL,PRESSURE,RADIUS,STRESS,VELOCITY,NXCHECK
      PARAMETER(DATA_TYPE = 4)  !Point Data

      CALL ENTERS('EXPOIN',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)

C GMH 22/10/95 Repeat for each case to name groups correctly
        WRITE(CHAR2,'(I5)') 1+10000
        WRITE(CHAR3,'(I5)') NDT+10000
        CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
        CALL STRING_TRIM(CHAR3,IBEG3,IEND3)
        POINT_NAME=CHAR2(IBEG2:IEND2)//'..'//CHAR3(IBEG3:IEND3)
        CALL STRING_TRIM(POINT_NAME,IBEG4,IEND4)

C---------------------------------------------------------------------

C#### Command: FEM export points<;FILENAME[default]> gauss
C###  Parameter:      <group NAME>
C###    Specify Gauss point group to export.
C###  Parameter:      <3dstress>
C###    Export stresses with the points.
C###  Parameter:      <as NAME[10001..10000]>
C###    Label the file with a character name.
C###  Parameter:      <element LIST[all]>
C###    Limit to Gauss points in the specified elements.
C###  Parameter:      <region #[1]>
C###    Specify the region number.
C###  Parameter:      <to (datafile/Motif)[datafile]>
C###    Specify the destination file format (`datafile' for exporting to
C###    CMGUI).
C###  Parameter:      <offset OFFSET[10000]>
C###    Start labelling points from OFFSET+1.
C###  Description:
C###    Exports points for each Gauss point, normally to CMGUI.
C###  Parameter"      <nocheck>
C###    Exports grid points without checking if a problem type has been
C###    defined (nx).

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME['//FILE00(IBEG1:IEND1)
     '    //']> gauss'
C GMH 27/10/95 Hard coded for 3D so far...
        OP_STRING(2)=BLANK(1:15)//'<group NAME>'
        OP_STRING(3)=BLANK(1:15)//'<3dstress>'
        OP_STRING(4)=BLANK(1:15)//'<as NAME['//POINT_NAME(IBEG4:IEND4)
     '    //']>'
        OP_STRING(5)=BLANK(1:15)//'<element LIST[all]>'
        OP_STRING(6)=BLANK(1:15)//'<region #[1]>'
        OP_STRING(7)=BLANK(1:15)//'<to (datafile/Motif)[datafile]>'
        OP_STRING(8)=BLANK(1:15)//'<offset OFFSET[10000]>'
        OP_STRING(9)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(10)=BLANK(1:15)//'<nocheck>'        
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C GMH 22/10/95 Repeat for each case to name groups correctly
        WRITE(CHAR2,'(I5)') 1+10000
        WRITE(CHAR3,'(I5)') NQT+10000
        CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
        CALL STRING_TRIM(CHAR3,IBEG3,IEND3)
        POINT_NAME=CHAR2(IBEG2:IEND2)//'..'//CHAR3(IBEG3:IEND3)
        CALL STRING_TRIM(POINT_NAME,IBEG4,IEND4)

C---------------------------------------------------------------------

C#### Command: FEM export points<;FILENAME[default]> grid
C###  Description:
C###    Export grid points as nodes to CMGUI.
C###  Parameter:      <as NAME[10001..10000]>
C###    Label the file with a character name.
C###  Parameter:      <number LIST[all]>
C###    Specify the points to export.
C###  Parameter:      <region #[1]>
C###    Specify the region number.
C###  Parameter:      <to (datafile/Motif)[datafile]>
C###    Specify the destination file format (`datafile' for exporting to
C###    CMGUI).
C###  Parameter:      <offset OFFSET[10000]>
C###    Add OFFSET to point numbers.
C###  Parameter:      <class>
C###    Specify the class number (of solve type) to use.
C###  Parameter:      <activation_time>
C###    Export the activation time field for activation problems
C###  Parameter:      <material>
C###    Export the cellular material parameters
C###  Parameter:      <potential>
C###    Export the potential field for activation problems
C###  Parameter:      <pressure>
C###    Export the pressure field for coronary problems
C###  Parameter:      <radius>
C###    Export the radius field
C###  Parameter:      <velocity>
C###    Export the velocity field for coronary problems
C###  Parameter:      <path>
C###    Export the path field for coronary problems
C###  Parameter"      <nocheck>
C###    Exports grid points without checking if a problem type has been
C###    defined (nx).

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME['//FILE00(IBEG1:IEND1)
     '    //']> grid'
        OP_STRING(2)=BLANK(1:15)//'<as NAME['//POINT_NAME(IBEG4:IEND4)
     '    //']>'
        OP_STRING(3)=BLANK(1:15)//'<number LIST[all]>'
        OP_STRING(4)=BLANK(1:15)//'<region #[1]>'
        OP_STRING(5)=BLANK(1:15)//'<to (datafile/Motif)[datafile]>'
        OP_STRING(6)=BLANK(1:15)//'<offset OFFSET[10000]>'
        OP_STRING(7)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(8)=BLANK(1:15)//'<potential>'
        OP_STRING(9)=BLANK(1:15)//'<activation_time>'
        OP_STRING(10)=BLANK(1:15)//'<pressure>'
        OP_STRING(11)=BLANK(1:15)//'<radius>'
        OP_STRING(12)=BLANK(1:15)//'<velocity>'
        OP_STRING(13)=BLANK(1:15)//'<path>'
        OP_STRING(14)=BLANK(1:15)//'<material>'
        OP_STRING(15)=BLANK(1:15)//'<nocheck>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe26','doc','EXPOIN',ERROR,*9999)
      ELSE
        nc=1 !Temporary MPN 12-Nov-94
        STRESS=.FALSE.

        
        IF(CBBREV(CO,'DATA',2,noco+1,NTCO,N3CO)) THEN
          POINT_TYPE='DATA'
          CALL ASSERT(.FALSE.,'>>Use export data',ERROR,*9999)
        ELSE IF(CBBREV(CO,'GAUSS',2,noco+1,NTCO,N3CO)) THEN
          POINT_TYPE='GAUSS'
        ELSE IF(CBBREV(CO,'GRID',2,noco+1,NTCO,N3CO)) THEN
          POINT_TYPE='GRID'
        ELSE
          CO(noco+1)='?'
          GO TO 1
        ENDIF
        
        IF(CBBREV(CO,'NOCHECK',2,noco+1,NTCO,N3CO)) THEN
          NXCHECK=.FALSE.
          nx=1
        ELSE
          NXCHECK=.TRUE.
        ENDIF

C new GMH 27/10/95 Adding export stress to point gauss
        IF(POINT_TYPE(1:5).EQ.'GAUSS') THEN
          IF(CBBREV(CO,'3DSTRESS',1,noco+1,NTCO,N3CO)) THEN
            STRESS=.TRUE.
          ELSE
            STRESS=.FALSE.
          ENDIF
        ENDIF

C GMH 13/7/96 If we are exporting data, then we do not have any nx
C dependence.
C KAT 15Feb99: point_type != data
C        IF(POINT_TYPE(1:4).NE.'DATA') THEN
C CS 11Jun99 if not exporting Gauss stress then do not need nx depenence
C OR 19-06-07 That seems still to be wrong. To implement the comments from
C     above for GAUSS and GRID points, it is necessary to first initialise 
C     the variables STRESS to FALSE and nx=1. Further, it is only necessary 
C     to check if nx is defined if gauss points and stress is outputted. 
C     STRESS in all other cases set to FALSE. Hence removed the lines below
C     and included them above after Stress=.TRUE.
C     
C        IF((POINT_TYPE(1:5).NE.'GAUSS').OR.STRESS) THEN
C          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
C          nxc=NXLIST(1)
C          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
C          CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
C     '      ERROR,*9999)
C        ENDIF
C        ELSE
CC         Ensure it is not used - give seg fault
C          nx=0
C        ENDIF !not data

C news MPN 9-11-95: OFFSET option to add const to all points
        IF(CBBREV(CO,'OFFSET',2,noco+1,NTCO,N3CO)) THEN
          OFFSET=IFROMC(CO(N3CO+1))
        ELSE
          OFFSET=10000
        ENDIF
        WRITE(CHAR2,'(I5)') 1+OFFSET
        WRITE(CHAR3,'(I5)') NQT+OFFSET
        CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
        CALL STRING_TRIM(CHAR3,IBEG3,IEND3)
        POINT_NAME=CHAR2(IBEG2:IEND2)//'..'//CHAR3(IBEG3:IEND3)
        CALL STRING_TRIM(POINT_NAME,IBEG4,IEND4)

C new GMH 19/9/95 Adding export error to data point
C KAT 15Feb99: point_type != data
C        IF(POINT_TYPE(1:4).EQ.'DATA') THEN
C          IF(CBBREV(CO,'ERROR',1,noco+1,NTCO,N3CO)) THEN
C            ERR=.TRUE.
C          ELSE
C            ERR=.FALSE.
C          ENDIF
C        ENDIF

        IF(CBBREV(CO,'REGION',1,noco+1,NTCO,N3CO)) THEN
          nr=IFROMC(CO(N3CO+1))
        ELSE
          nr=1
        ENDIF

        IF(POINT_TYPE(1:4).EQ.'GRID') THEN
          IF(CBBREV(CO,'ACTIVATION_TIME',3,noco+1,NTCO,N3CO)) THEN
            ACTIVATION_TIME=.TRUE.
          ELSE
            ACTIVATION_TIME=.FALSE.
          ENDIF
          IF(CBBREV(CO,'POTENTIAL',3,noco+1,NTCO,N3CO)) THEN
            POTENTIAL=.TRUE.
          ELSE
            POTENTIAL=.FALSE.
          ENDIF
          IF(CBBREV(CO,'MATERIAL',3,noco+1,NTCO,N3CO)) THEN
            MATERIAL=.TRUE.
            CALL ASSERT(CALL_CELL,'>>Must define cell first',
     '        ERROR,*9999)
            CALL ASSERT(CALL_CELL_MATE,
     '        '>>Must define cell materials first',ERROR,*9999)
          ELSE
            MATERIAL=.FALSE.
          ENDIF
          IF(CBBREV(CO,'PRESSURE',3,noco+1,NTCO,N3CO)) THEN
            PRESSURE=.TRUE.
          ELSE
            PRESSURE=.FALSE.
          ENDIF
          IF(CBBREV(CO,'RADIUS',3,noco+1,NTCO,N3CO)) THEN
            RADIUS=.TRUE.
          ELSE
            RADIUS=.FALSE.
          ENDIF
          IF(CBBREV(CO,'VELOCITY',3,noco+1,NTCO,N3CO)) THEN
            VELOCITY=.TRUE.
          ELSE
            VELOCITY=.FALSE.
          ENDIF
          IF(CBBREV(CO,'PATH',3,noco+1,NTCO,N3CO)) THEN
            PATH=.TRUE.
          ELSE
            PATH=.FALSE.
          ENDIF
C
C     OR 21-06-07 It is not essential to have a NX_SOLVE class defined
C     when writing out grid points, e.g. if non of the other BOOLEAN
C     variables is specified as TRUE. 
C          
          IF (NXCHECK) THEN 
            CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
            nxc=NXLIST(1)
            CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     &           ERROR,*9999)
          ENDIF
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

        IF(CBBREV(CO,'NUMBER' ,2,noco+1,NTCO,N3CO).OR
     '    .CBBREV(CO,'ELEMENT',2,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),100,LIST(0),LIST(1),ERROR,*9999)
          WRITE(CHAR1,'(I5)') LIST(1)+OFFSET
          WRITE(CHAR2,'(I5)') LIST(LIST(0))+OFFSET
          CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
          CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
          POINT_NAME=CHAR1(IBEG1:IEND1)//'..'//CHAR2(IBEG2:IEND2)
!        ELSE IF(CBBREV(CO,'GROUP',2,noco+1,NTCO,N3CO)) THEN
!          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
!          POINT_NAME=CO(N3CO+1)(IBEG:IEND)
        ELSE
          IF(POINT_TYPE(1:5).EQ.'NODES') THEN
            LIST(0)=NPNODE(0,nr)
            DO NONODE=1,NPNODE(0,nr)
              LIST(NONODE)=NPNODE(NONODE,nr)
            ENDDO
C KAT 15Feb99: point_type != data
C          ELSE IF(POINT_TYPE(1:4).EQ.'DATA') THEN
C            LIST(0)=NDT
C            DO nd=1,NDT
C              LIST(nd)=nd
C            ENDDO
          ELSE IF(POINT_TYPE(1:5).EQ.'GAUSS') THEN
            IF(CBBREV(CO,'GROUP',1,noco+1,NTCO,N3CO)) THEN
              CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
              CALL CUPPER(CO(N3CO+1)(IBEG:IEND),CHAR)
              N1GRGA=0
              DO nogrga=1,NTGRGA
                CALL CUPPER(LAGRGA(nogrga),LABEL)
                CALL STRING_TRIM(LABEL,IBEG2,IEND2)
                IF(CHAR(IBEG:IEND).EQ.LABEL(IBEG2:IEND2)) THEN
                  N1GRGA=nogrga !is existing group label ID
                ENDIF
              ENDDO
              nogrga=N1GRGA
              LIST(0)=ILISTMBR(%VAL(LIGRGA_PTR(nogrga)),1)
              GROUPED=.TRUE.
            ELSE
              LIST(0)=NET(nr)
              DO ne=1,NET(nr)
                LIST(ne)=ne
              ENDDO
              GROUPED=.FALSE.
            ENDIF
          ELSE IF(POINT_TYPE(1:5).EQ.'GRID') THEN
            LIST(0)=NQT
            DO nq=1,NQT
              LIST(nq)=nq
            ENDDO
          ENDIF
C news GMH 22/10/95 NAME option for exnode file
          IF(CBBREV(CO,'AS',1,noco+1,NTCO,N3CO)) THEN
            CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
            POINT_NAME=CO(N3CO+1)(IBEG1:IEND1)
          ELSE
            WRITE(CHAR1,'(I5)') LIST(1)+OFFSET
            WRITE(CHAR2,'(I5)') LIST(LIST(0))+OFFSET
            CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
            CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
            POINT_NAME=CHAR1(IBEG1:IEND1)//'..'//CHAR2(IBEG2:IEND2)
          ENDIF
        ENDIF

        IF((OUTPUT(1:8).EQ.'DATAFILE').OR.(OUTPUT(1:5).EQ.'MOTIF')) THEN
          CALL CHECKF(1,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          IF(LIST(0).GT.0) THEN
            IF(DATAFILE) THEN
              CALL STRING_TRIM(FILE,IBEG,IEND)
              CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.exnode','NEW',
     '        'SEQUEN','FORMATTED',132,ERROR,*9999)
C**           write the group name
              CALL STRING_TRIM(POINT_NAME,IBEG,IEND)
              WRITE(IFILE,'( '' Group name: '',A)')
     '          POINT_NAME(IBEG:IEND)
            ELSE
              IF(USE_SOCKET) THEN
C**             send the data type identifier
                IF(FSKWRITE(DATA_TYPE,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
C**             send the group name
                CALL STRING_TRIM(POINT_NAME,IBEG,IEND)
                CLEN=FSKLEN(POINT_NAME(IBEG:IEND))
                CALL FSKF2C(POINT_NAME(IBEG:IEND),CLEN,INTSTR)
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
C KAT 15Feb99: point_type != data
C            IF(POINT_TYPE(1:4).EQ.'DATA') THEN
C              IF(ERR) THEN
C                IF(DATAFILE) THEN
C                  WRITE(IFILE,'(1X,''#Fields=2'')')
C                ELSE
C                  ISOCKET(1)=2
C                  IF(FSKWRITE(ISOCKET,SK_LONG_INT,1,CONNID2).EQ.-1)
C     '              GOTO 9999
C                ENDIF
C              ELSE
C                IF(DATAFILE) THEN
C                  WRITE(IFILE,'(1X,''#Fields=1'')')
C                ELSE
C                  IF(FSKWRITE(ICHAR('#'),SK_CHAR,1,CONNID2).EQ.-1)
C     '              GOTO 9999
C                ENDIF
C              ENDIF
C            ELSE
            IF(POINT_TYPE(1:5).EQ.'GAUSS') THEN
              IF(STRESS) THEN
                IF(DATAFILE) THEN
                  WRITE(IFILE,'(1X,''#Fields=4'')')
                ELSE
                  ISOCKET(1)=4
                  IF(FSKWRITE(ISOCKET,SK_LONG_INT,1,CONNID2).EQ.-1)
     '              GOTO 9999
                ENDIF
              ELSE
                IF(DATAFILE) THEN
                  WRITE(IFILE,'(1X,''#Fields=1'')')
                ELSE
                  IF(FSKWRITE(ICHAR('#'),SK_CHAR,1,CONNID2).EQ.-1)
     '              GOTO 9999
                ENDIF
              ENDIF
            ELSE IF(POINT_TYPE(1:4).EQ.'GRID') THEN
              IF(DATAFILE) THEN
                NUMFIELDS=1
                IF(ACTIVATION_TIME) NUMFIELDS=NUMFIELDS+1
                IF(POTENTIAL) NUMFIELDS=NUMFIELDS+1
                IF(PRESSURE) NUMFIELDS=NUMFIELDS+1
                IF(RADIUS) NUMFIELDS=NUMFIELDS+1
                IF(VELOCITY) NUMFIELDS=NUMFIELDS+1
                IF(PATH) NUMFIELDS=NUMFIELDS+1
                IF(MATERIAL) NUMFIELDS=NUMFIELDS+1+1+ !cell_type+modelID
     '            CELL_NUM_STATE(1)+(NQIT-1)+NQRT
C PM 26-JUL-01 : flow through elastic tube
                IF((ITYP3(nr,nx).EQ.1).AND.(ITYP2(nr,nx).EQ.5).AND.
     '            (.NOT.PRESSURE).AND.(.NOT.RADIUS).AND.
     '            (.NOT.VELOCITY)) NUMFIELDS=4

C *** DPN 13 September 1999 - need to be more general!
c                CALL ASSERT(NUMFIELDS.LT.1000,' >>Too many fields',
c     '            ERROR,*9999)  !fields written into I1
                WRITE(IFILE,'(1X,''#Fields='',I12)') NUMFIELDS
              ELSE
                ISOCKET(1)=4
                IF(FSKWRITE(ISOCKET,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
              ENDIF
            ENDIF

            IF(DATAFILE) THEN
              IF(POINT_TYPE(1:4).EQ.'GRID') THEN
                WRITE(IFILE,'(1X,''1) coordinates, coordinate, '
     '            //'rectangular cartesian, #Components='',I1)')
     '            NJ_LOC(NJL_GEOM,0,nr)
                DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                  IF(nj.EQ.1) THEN
                    WRITE(IFILE,'(1X,''  x.  Value index= 1, '
     '                //'#Derivatives=0'')')
                  ELSEIF(nj.EQ.2) THEN
                    WRITE(IFILE,'(1X,''  y.  Value index= 2, '
     '                //'#Derivatives=0'')')
                  ELSEIF(nj.EQ.3) THEN
                    WRITE(IFILE,'(1X,''  z.  Value index= 3, '
     '                //'#Derivatives=0'')')
                  ENDIF
                ENDDO

              ELSE
                WRITE(IFILE,'(1X,''1) coordinates, coordinate, '
     '            //'rectangular cartesian, #Components=3'')')
                DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                  IF(nj.EQ.1) THEN
                    WRITE(IFILE,'(1X,''  x.  Value index= 1, '
     '                //'#Derivatives=0'')')
                  ELSEIF(nj.EQ.2) THEN
                    WRITE(IFILE,'(1X,''  y.  Value index= 2, '
     '                //'#Derivatives=0'')')
                  ELSEIF(nj.EQ.3) THEN
                    WRITE(IFILE,'(1X,''  z.  Value index= 3, '
     '                //'#Derivatives=0'')')
                  ENDIF
                ENDDO
              ENDIF
            ELSE
              FIELD_NAME='coordinates'
              CALL STRING_TRIM(FIELD_NAME,IBEG,IEND)
              CLEN=FSKLEN(FIELD_NAME(IBEG:IEND))
              CALL FSKF2C(FIELD_NAME(IBEG:IEND),CLEN,INTSTR)
              IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
     '          GOTO 9999
              IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
     '          GOTO 9999
              ISOCKET(1)=1
              ISOCKET(2)=1
              ISOCKET(3)=3
              IF(FSKWRITE(ISOCKET,SK_LONG_INT,3,CONNID2).EQ.-1)
     '          GOTO 9999
              COMPONENT_NAME='x'
              CALL STRING_TRIM(COMPONENT_NAME,IBEG,IEND)
              CLEN=FSKLEN(COMPONENT_NAME(IBEG:IEND))
              CALL FSKF2C(COMPONENT_NAME(IBEG:IEND),CLEN,INTSTR)
              IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
     '          GOTO 9999
              IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
     '          GOTO 9999
              ISOCKET(1)=1
              ISOCKET(2)=0
              IF(FSKWRITE(ISOCKET,SK_LONG_INT,2,CONNID2).EQ.-1)
     '          GOTO 9999
              COMPONENT_NAME='y'
              CALL STRING_TRIM(COMPONENT_NAME,IBEG,IEND)
              CLEN=FSKLEN(COMPONENT_NAME(IBEG:IEND))
              CALL FSKF2C(COMPONENT_NAME(IBEG:IEND),CLEN,INTSTR)
              IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
     '          GOTO 9999
              IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
     '          GOTO 9999
              ISOCKET(1)=2
              ISOCKET(2)=0
              IF(FSKWRITE(ISOCKET,SK_LONG_INT,2,CONNID2).EQ.-1)
     '          GOTO 9999
              COMPONENT_NAME='z'
              CALL STRING_TRIM(COMPONENT_NAME,IBEG,IEND)
              CLEN=FSKLEN(COMPONENT_NAME(IBEG:IEND))
              CALL FSKF2C(COMPONENT_NAME(IBEG:IEND),CLEN,INTSTR)
              IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
              IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
     '          GOTO 9999
              ISOCKET(1)=3
              ISOCKET(2)=0
              IF(FSKWRITE(ISOCKET,SK_LONG_INT,2,CONNID2).EQ.-1)
     '          GOTO 9999
            ENDIF
C KAT 15Feb99: point_type != data
C            IF((POINT_TYPE(1:4).EQ.'DATA').AND.ERR) THEN
C              IF(DATAFILE) THEN
C                WRITE(IFILE,'(1X,''2) error, field, '
C     '            //'rectangular cartesian, #Components=3'')')
C                WRITE(IFILE,'(1X,''  x.  Value index= 4, '
C     '            //'#Derivatives=0'')')
C                WRITE(IFILE,'(1X,''  y.  Value index= 5, '
C     '            //'#Derivatives=0'')')
C                WRITE(IFILE,'(1X,''  z.  Value index= 6, '
C     '            //'#Derivatives=0'')')
C              ELSE
C                FIELD_NAME='error'
C                CALL STRING_TRIM(FIELD_NAME,IBEG,IEND)
C                CLEN=FSKLEN(FIELD_NAME(IBEG:IEND))
C                CALL FSKF2C(FIELD_NAME(IBEG:IEND),CLEN,INTSTR)
C                IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
C     '            GOTO 9999
C                IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
C     '            GOTO 9999
C                ISOCKET(1)=3 !field
C                ISOCKET(2)=1 !rectangular cartesian
C                ISOCKET(3)=3 !components
C                IF(FSKWRITE(ISOCKET,SK_LONG_INT,3,CONNID2).EQ.-1)
C     '            GOTO 9999
C                COMPONENT_NAME='x'
C                CALL STRING_TRIM(COMPONENT_NAME,IBEG,IEND)
C                CLEN=FSKLEN(COMPONENT_NAME(IBEG:IEND))
C                CALL FSKF2C(COMPONENT_NAME(IBEG:IEND),CLEN,INTSTR)
C                IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
C     '            GOTO 9999
C                IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
C     '            GOTO 9999
C                ISOCKET(1)=4
C                ISOCKET(2)=0
C                IF(FSKWRITE(ISOCKET,SK_LONG_INT,2,CONNID2).EQ.-1)
C     '            GOTO 9999
C                COMPONENT_NAME='y'
C                CALL STRING_TRIM(COMPONENT_NAME,IBEG,IEND)
C                CLEN=FSKLEN(COMPONENT_NAME(IBEG:IEND))
C                CALL FSKF2C(COMPONENT_NAME(IBEG:IEND),CLEN,INTSTR)
C                IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
C     '            GOTO 9999
C                IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
C     '            GOTO 9999
C                ISOCKET(1)=5
C                ISOCKET(2)=0
C                IF(FSKWRITE(ISOCKET,SK_LONG_INT,2,CONNID2).EQ.-1)
C     '            GOTO 9999
C                COMPONENT_NAME='z'
C                CALL STRING_TRIM(COMPONENT_NAME,IBEG,IEND)
C                CLEN=FSKLEN(COMPONENT_NAME(IBEG:IEND))
C                CALL FSKF2C(COMPONENT_NAME(IBEG:IEND),CLEN,INTSTR)
C                IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
C     '            GOTO 9999
C                IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
C     '            GOTO 9999
C                ISOCKET(1)=6
C                ISOCKET(2)=0
C                IF(FSKWRITE(ISOCKET,SK_LONG_INT,2,CONNID2).EQ.-1)
C     '            GOTO 9999
C              ENDIF
C            ENDIF
            IF((POINT_TYPE(1:5).EQ.'GAUSS').AND.STRESS) THEN
              IF(DATAFILE) THEN
                WRITE(IFILE,'(1X,''2) stress1, field, '
     '            //'rectangular cartesian, #Components=3'')')
                WRITE(IFILE,'(1X,''  x.  Value index= 4, '
     '            //'#Derivatives=0'')')
                WRITE(IFILE,'(1X,''  y.  Value index= 5, '
     '            //'#Derivatives=0'')')
                WRITE(IFILE,'(1X,''  z.  Value index= 6, '
     '            //'#Derivatives=0'')')
                WRITE(IFILE,'(1X,''3) stress2, field, '
     '            //'rectangular cartesian, #Components=3'')')
                WRITE(IFILE,'(1X,''  x.  Value index= 7, '
     '            //'#Derivatives=0'')')
                WRITE(IFILE,'(1X,''  y.  Value index= 8, '
     '            //'#Derivatives=0'')')
                WRITE(IFILE,'(1X,''  z.  Value index= 9, '
     '            //'#Derivatives=0'')')
                WRITE(IFILE,'(1X,''4) stress3, field, '
     '            //'rectangular cartesian, #Components=3'')')
                WRITE(IFILE,'(1X,''  x.  Value index= 10, '
     '            //'#Derivatives=0'')')
                WRITE(IFILE,'(1X,''  y.  Value index= 11, '
     '            //'#Derivatives=0'')')
                WRITE(IFILE,'(1X,''  z.  Value index= 12, '
     '            //'#Derivatives=0'')')
              ELSE
                FIELD_NAME='stress1'
                CALL STRING_TRIM(FIELD_NAME,IBEG,IEND)
                CLEN=FSKLEN(FIELD_NAME(IBEG:IEND))
                CALL FSKF2C(FIELD_NAME(IBEG:IEND),CLEN,INTSTR)
                IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
                IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
     '            GOTO 9999
                ISOCKET(1)=3 !field
                ISOCKET(2)=1 !rectangular cartesian
                ISOCKET(3)=3 !components
                IF(FSKWRITE(ISOCKET,SK_LONG_INT,3,CONNID2).EQ.-1)
     '            GOTO 9999
                COMPONENT_NAME='x'
                CALL STRING_TRIM(COMPONENT_NAME,IBEG,IEND)
                CLEN=FSKLEN(COMPONENT_NAME(IBEG:IEND))
                CALL FSKF2C(COMPONENT_NAME(IBEG:IEND),CLEN,INTSTR)
                IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
                IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
     '            GOTO 9999
                ISOCKET(1)=4
                ISOCKET(2)=0
                IF(FSKWRITE(ISOCKET,SK_LONG_INT,2,CONNID2).EQ.-1)
     '            GOTO 9999
                COMPONENT_NAME='y'
                CALL STRING_TRIM(COMPONENT_NAME,IBEG,IEND)
                CLEN=FSKLEN(COMPONENT_NAME(IBEG:IEND))
                CALL FSKF2C(COMPONENT_NAME(IBEG:IEND),CLEN,INTSTR)
                IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
                IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
     '            GOTO 9999
                ISOCKET(1)=5
                ISOCKET(2)=0
                IF(FSKWRITE(ISOCKET,SK_LONG_INT,2,CONNID2).EQ.-1)
     '            GOTO 9999
                COMPONENT_NAME='z'
                CALL STRING_TRIM(COMPONENT_NAME,IBEG,IEND)
                CLEN=FSKLEN(COMPONENT_NAME(IBEG:IEND))
                CALL FSKF2C(COMPONENT_NAME(IBEG:IEND),CLEN,INTSTR)
                IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
                IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
     '            GOTO 9999
                ISOCKET(1)=6
                ISOCKET(2)=0
                IF(FSKWRITE(ISOCKET,SK_LONG_INT,2,CONNID2).EQ.-1)
     '            GOTO 9999
                FIELD_NAME='stress2'
                CALL STRING_TRIM(FIELD_NAME,IBEG,IEND)
                CLEN=FSKLEN(FIELD_NAME(IBEG:IEND))
                CALL FSKF2C(FIELD_NAME(IBEG:IEND),CLEN,INTSTR)
                IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
                IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
     '            GOTO 9999
                ISOCKET(1)=3 !field
                ISOCKET(2)=1 !rectangular cartesian
                ISOCKET(3)=3 !components
                IF(FSKWRITE(ISOCKET,SK_LONG_INT,3,CONNID2).EQ.-1)
     '            GOTO 9999
                COMPONENT_NAME='x'
                CALL STRING_TRIM(COMPONENT_NAME,IBEG,IEND)
                CLEN=FSKLEN(COMPONENT_NAME(IBEG:IEND))
                CALL FSKF2C(COMPONENT_NAME(IBEG:IEND),CLEN,INTSTR)
                IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
                IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
     '            GOTO 9999
                ISOCKET(1)=7
                ISOCKET(2)=0
                IF(FSKWRITE(ISOCKET,SK_LONG_INT,2,CONNID2).EQ.-1)
     '            GOTO 9999
                COMPONENT_NAME='y'
                CALL STRING_TRIM(COMPONENT_NAME,IBEG,IEND)
                CLEN=FSKLEN(COMPONENT_NAME(IBEG:IEND))
                CALL FSKF2C(COMPONENT_NAME(IBEG:IEND),CLEN,INTSTR)
                IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
                IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
     '            GOTO 9999
                ISOCKET(1)=8
                ISOCKET(2)=0
                IF(FSKWRITE(ISOCKET,SK_LONG_INT,2,CONNID2).EQ.-1)
     '            GOTO 9999
                COMPONENT_NAME='z'
                CALL STRING_TRIM(COMPONENT_NAME,IBEG,IEND)
                CLEN=FSKLEN(COMPONENT_NAME(IBEG:IEND))
                CALL FSKF2C(COMPONENT_NAME(IBEG:IEND),CLEN,INTSTR)
                IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
                IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
     '            GOTO 9999
                ISOCKET(1)=9
                ISOCKET(2)=0
                IF(FSKWRITE(ISOCKET,SK_LONG_INT,2,CONNID2).EQ.-1)
     '            GOTO 9999
                FIELD_NAME='stress3'
                CALL STRING_TRIM(FIELD_NAME,IBEG,IEND)
                CLEN=FSKLEN(FIELD_NAME(IBEG:IEND))
                CALL FSKF2C(FIELD_NAME(IBEG:IEND),CLEN,INTSTR)
                IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
                IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
     '            GOTO 9999
                ISOCKET(1)=3 !field
                ISOCKET(2)=1 !rectangular cartesian
                ISOCKET(3)=3 !components
                IF(FSKWRITE(ISOCKET,SK_LONG_INT,3,CONNID2).EQ.-1)
     '            GOTO 9999
                COMPONENT_NAME='x'
                CALL STRING_TRIM(COMPONENT_NAME,IBEG,IEND)
                CLEN=FSKLEN(COMPONENT_NAME(IBEG:IEND))
                CALL FSKF2C(COMPONENT_NAME(IBEG:IEND),CLEN,INTSTR)
                IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
                IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
     '            GOTO 9999
                ISOCKET(1)=10
                ISOCKET(2)=0
                IF(FSKWRITE(ISOCKET,SK_LONG_INT,2,CONNID2).EQ.-1)
     '            GOTO 9999
                COMPONENT_NAME='y'
                CALL STRING_TRIM(COMPONENT_NAME,IBEG,IEND)
                CLEN=FSKLEN(COMPONENT_NAME(IBEG:IEND))
                CALL FSKF2C(COMPONENT_NAME(IBEG:IEND),CLEN,INTSTR)
                IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
                IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
     '            GOTO 9999
                ISOCKET(1)=11
                ISOCKET(2)=0
                IF(FSKWRITE(ISOCKET,SK_LONG_INT,2,CONNID2).EQ.-1)
     '            GOTO 9999
                COMPONENT_NAME='z'
                CALL STRING_TRIM(COMPONENT_NAME,IBEG,IEND)
                CLEN=FSKLEN(COMPONENT_NAME(IBEG:IEND))
                CALL FSKF2C(COMPONENT_NAME(IBEG:IEND),CLEN,INTSTR)
                IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
                IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
     '            GOTO 9999
                ISOCKET(1)=12
                ISOCKET(2)=0
                IF(FSKWRITE(ISOCKET,SK_LONG_INT,2,CONNID2).EQ.-1)
     '            GOTO 9999
              ENDIF
            ENDIF
            IF(POINT_TYPE(1:4).EQ.'GRID') THEN
              IF(DATAFILE) THEN
                NUMFIELDS=2
                VAL_INDEX=NJ_LOC(NJL_GEOM,0,nr)+1
                IF(ACTIVATION_TIME) THEN
                  WRITE(IFILE,'(1X,I1,'') activation_time, field, '
     '              //'rectangular cartesian, #Components=1'')')
     '              NUMFIELDS
                  WRITE(IFILE,'(1X,''  activation_time.  '
     '              //'Value index='',I2,'', #Derivatives=0'')')
     '              VAL_INDEX
                  NUMFIELDS=NUMFIELDS+1
                  VAL_INDEX=VAL_INDEX+1
                ENDIF
                IF(POTENTIAL) THEN
                  WRITE(IFILE,'(1X,I1,'') potential, field, '
     '              //'rectangular cartesian, #Components=1'')')
     '              NUMFIELDS
                  WRITE(IFILE,'(1X,''  potential.  '
     '              //'Value index='',I2,'', #Derivatives=0'')')
     '              VAL_INDEX
                  NUMFIELDS=NUMFIELDS+1
                  VAL_INDEX=VAL_INDEX+1
                ENDIF
                IF(PRESSURE) THEN
                  WRITE(IFILE,'(1X,I1,'') pressure, field, '
     '              //'rectangular cartesian, #Components=1'')')
     '              NUMFIELDS
                  WRITE(IFILE,'(1X,''  pressure.  '
     '              //'Value index='',I2,'', #Derivatives=0'')')
     '              VAL_INDEX
                  NUMFIELDS=NUMFIELDS+1
                  VAL_INDEX=VAL_INDEX+1
                ENDIF
                IF(RADIUS) THEN
                  WRITE(IFILE,'(1X,I1,'') radius, field, '
     '              //'rectangular cartesian, #Components=1'')')
     '              NUMFIELDS
                  WRITE(IFILE,'(1X,''  radius.  '
     '              //'Value index='',I2,'', #Derivatives=0'')')
     '              VAL_INDEX
                  NUMFIELDS=NUMFIELDS+1
                  VAL_INDEX=VAL_INDEX+1
                ENDIF
                IF(VELOCITY) THEN
                  WRITE(IFILE,'(1X,I1,'') velocity, field, '
     '              //'rectangular cartesian, #Components=1'')')
     '              NUMFIELDS
                  WRITE(IFILE,'(1X,''  velocity.  '
     '              //'Value index='',I2,'', #Derivatives=0'')')
     '              VAL_INDEX
                  NUMFIELDS=NUMFIELDS+1
                  VAL_INDEX=VAL_INDEX+1
                ENDIF
                IF(PATH) THEN
                  WRITE(IFILE,'(1X,I1,'') path, field, '
     '              //'rectangular cartesian, #Components=1'')')
     '              NUMFIELDS
                  WRITE(IFILE,'(1X,''  path.  '
     '              //'Value index='',I2,'', #Derivatives=0'')')
     '              VAL_INDEX
                  NUMFIELDS=NUMFIELDS+1
                  VAL_INDEX=VAL_INDEX+1
                ENDIF

C PM 26-JUL-01 : fluid in elastic tube
                IF((ITYP3(nr,nx).EQ.1).AND.(ITYP2(nr,nx).EQ.5).AND.
     '            (.NOT.PRESSURE).AND.(.NOT.RADIUS).AND.
     '            (.NOT.VELOCITY)) THEN
                  WRITE(IFILE,'(1X,I1,'') pressure, field, '
     '              //'rectangular cartesian, #Components=1'')')
     '              NUMFIELDS
                  WRITE(IFILE,'(1X,''  pressure.  '
     '              //'Value index='',I2,'', #Derivatives=0'')')
     '              VAL_INDEX
                  NUMFIELDS=NUMFIELDS+1
                  VAL_INDEX=VAL_INDEX+1
                  WRITE(IFILE,'(1X,I1,'') radius, field, '
     '              //'rectangular cartesian, #Components=1'')')
     '              NUMFIELDS
                  WRITE(IFILE,'(1X,''  radius.  '
     '              //'Value index='',I2,'', #Derivatives=0'')')
     '              VAL_INDEX
                  NUMFIELDS=NUMFIELDS+1
                  VAL_INDEX=VAL_INDEX+1
                  WRITE(IFILE,'(1X,I1,'') velocity, field, '
     '              //'rectangular cartesian, #Components=1'')')
     '              NUMFIELDS
                  WRITE(IFILE,'(1X,''  velocity.  '
     '              //'Value index='',I2,'', #Derivatives=0'')')
     '              VAL_INDEX
                ENDIF

                IF(MATERIAL) THEN

C *** DPN 13 September 1999 - write out the headers for the cellular
C *** parameters:
C ***  - Model ID is a constant string field
C ***  - The cell type is a node based field (==variant number)
C ***  - State variables are always fields
C ***  - Integer and real parameters are indexed fields unless they are
C ***     spatially varying, in which case they are normal node fields

                  !First write out the header for the model name
C *** The model name is given by ITYP19_ITYP3_KTYP33, the individual
C     numbers can be obtained from the string and used to set the model
C     type in Cell.
                  WRITE(IFILE,'(1X,I1,'') model_id, field, constant, '
     '              //'string, #Components=1'')')
     '              NUMFIELDS
                  WRITE(IFILE,'(1X,''  ID.'')')
                  NUMFIELDS=NUMFIELDS+1
                  !The header for the cell type field
                  WRITE(IFILE,'(1X,I1,'') cell_type, field, integer, '
     '              //'#Components=1'')')
     '              NUMFIELDS
                  WRITE(IFILE,'(1X,''  number.  '
     '              //'Value index= 1, #Derivatives=0'')')
                  NUMFIELDS=NUMFIELDS+1
                  !The headers for the state variables
                  DO j=1,CELL_NUM_STATE(1)
                    CALL STRING_TRIM(CELL_YQS_NAMES(
     '                CELL_STATE_OFFSET(1)-1+j,1),IBEG,IEND)
                    WRITE(CHAR,'(I12)') NUMFIELDS
                    CALL STRING_TRIM(CHAR,IBEG1,IEND1)
                    WRITE(IFILE,'(1X,A,'') '',A,'', field, real, '
     '                //'#Components=1'')')
     '                CHAR(IBEG1:IEND1),
     '                CELL_YQS_NAMES(CELL_STATE_OFFSET(1)-1+j,
     '                1)(IBEG:IEND)
                    WRITE(IFILE,'(1X,''  number.  '
     '                //'Value index= 1, #Derivatives=0'')')
                    NUMFIELDS=NUMFIELDS+1
                  ENDDO !j
                  !The headers for the integer parameters
                  DO j=2,NQIT !skip variant number
                    CALL STRING_TRIM(CELL_ICQS_NAMES(j,1),IBEG,IEND)
                    WRITE(CHAR,'(I12)') NUMFIELDS
                    CALL STRING_TRIM(CHAR,IBEG1,IEND1)
                    CELL_FIELD=.FALSE.
                    DO k=1,IICQS_SPATIAL(0,1)
                      IF(IICQS_SPATIAL(k,1).EQ.j) THEN
                        CELL_FIELD=.TRUE.
                      ENDIF
                    ENDDO !k
                    IF(CELL_FIELD) THEN
                      !normal node field
                      WRITE(IFILE,'(1X,A,'') '',A,'', field, integer, '
     '                  //'#Components=1'')') CHAR(IBEG1:IEND1),
     '                  CELL_ICQS_NAMES(j,1)(IBEG:IEND)
                      WRITE(IFILE,'(1X,''  number.  '
     '                  //'Value index= 1, #Derivatives=0'')')
                    ELSE
                      !indexed field
                      WRITE(IFILE,'(1X,A,'') '',A,'', field, indexed, '
     '                  //'Index_field=cell_type, #Values='',I12,'', '
     '                  //'integer, #Components=1'')')
     '                  CHAR(IBEG1:IEND1),
     '                  CELL_ICQS_NAMES(j,1)(IBEG:IEND),
     '                  CELL_NUM_VARIANTS
                      WRITE(IFILE,'(1X,''  number.'')')
                    ENDIF
                    NUMFIELDS=NUMFIELDS+1
                  ENDDO !j
                  !The headers for the real parameters
                  DO j=1,NQRT
                    CALL STRING_TRIM(CELL_RCQS_NAMES(j,1),IBEG,IEND)
                    WRITE(CHAR,'(I12)') NUMFIELDS
                    CALL STRING_TRIM(CHAR,IBEG1,IEND1)
                    CELL_FIELD=.FALSE.
                    DO k=1,IRCQS_SPATIAL(0,1)
                      IF(IRCQS_SPATIAL(k,1).EQ.j) THEN
                        CELL_FIELD=.TRUE.
                      ENDIF
                    ENDDO !k
                    IF(CELL_FIELD) THEN
                      !normal node field
                      WRITE(IFILE,'(1X,A,'') '',A,'', field, real, '
     '                  //'#Components=1'')') CHAR(IBEG1:IEND1),
     '                  CELL_RCQS_NAMES(j,1)(IBEG:IEND)
                      WRITE(IFILE,'(1X,''  number.  '
     '                  //'Value index= 1, #Derivatives=0'')')
                    ELSE
                      !indexed field
                      WRITE(IFILE,'(1X,A,'') '',A,'', field, indexed, '
     '                  //'Index_field=cell_type, #Values='',I12,'', '
     '                  //'real, #Components=1'')')
     '                  CHAR(IBEG1:IEND1),
     '                  CELL_RCQS_NAMES(j,1)(IBEG:IEND),
     '                  CELL_NUM_VARIANTS
                      WRITE(IFILE,'(1X,''  number.'')')
                    ENDIF
                    NUMFIELDS=NUMFIELDS+1
                  ENDDO !j
                  !Now write out the values for the indexed fields
                  WRITE(IFILE,'(1X,''Values:'')')
                  !The model ID
                  WRITE(CHAR1,'(I5)') ITYP19(nr,nx)
                  CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
                  WRITE(CHAR2,'(I5)') ITYP3(nr,nx)
                  CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
                  WRITE(CHAR3,'(I5)') KTYP33
                  CALL STRING_TRIM(CHAR3,IBEG3,IEND3)
                  WRITE(CHAR,'(A,''_'',A,''_'',A)')
     '              CHAR1(IBEG1:IEND1),CHAR2(IBEG2:IEND2),
     '              CHAR3(IBEG3:IEND3)
                  CALL STRING_TRIM(CHAR,IBEG,IEND)
                  WRITE(IFILE,'(1X,A)') CHAR(IBEG:IEND)
                  !The integer parameters
                  DO j=2,NQIT !skip variant number
                    CELL_FIELD=.FALSE.
                    DO k=1,IICQS_SPATIAL(0,1)
                      IF(IICQS_SPATIAL(k,1).EQ.j) THEN
                        CELL_FIELD=.TRUE.
                      ENDIF
                    ENDDO !k
                    IF(CELL_FIELD) THEN
                      !normal node field
                      ! ** do nothing **
                    ELSE
                      !indexed field - write out the variant values
                      WRITE(IFILE,'(I12)')
     '                  (CELL_ICQS_VALUE(j,k),k=1,CELL_NUM_VARIANTS)
                    ENDIF
                  ENDDO !j
                  !The real parameters
                  DO j=1,NQRT
                    CELL_FIELD=.FALSE.
                    DO k=1,IRCQS_SPATIAL(0,1)
                      IF(IRCQS_SPATIAL(k,1).EQ.j) THEN
                        CELL_FIELD=.TRUE.
                      ENDIF
                    ENDDO !k
                    IF(CELL_FIELD) THEN
                      !normal node field
                      ! ** do nothing **
                    ELSE
                      !indexed field
                      WRITE(IFILE,'(9E13.5)')
     '                  (CELL_RCQS_VALUE(j,k),k=1,CELL_NUM_VARIANTS)
                    ENDIF
                  ENDDO !j
                ENDIF !MATERIAL
              ELSE
                FIELD_NAME='transmembrane_potential'
                CALL STRING_TRIM(FIELD_NAME,IBEG,IEND)
                CLEN=FSKLEN(FIELD_NAME(IBEG:IEND))
                CALL FSKF2C(FIELD_NAME(IBEG:IEND),CLEN,INTSTR)
                IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
                IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
     '            GOTO 9999
                ISOCKET(1)=3
                ISOCKET(2)=1
                ISOCKET(3)=1
                IF(FSKWRITE(ISOCKET,SK_LONG_INT,3,CONNID2).EQ.-1)
     '            GOTO 9999
                CALL STRING_TRIM(FIELD_NAME,IBEG,IEND)
                CLEN=FSKLEN(FIELD_NAME(IBEG:IEND))
                CALL FSKF2C(FIELD_NAME(IBEG:IEND),CLEN,INTSTR)
                IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
                IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
     '            GOTO 9999
                ISOCKET(1)=4
                ISOCKET(2)=0
                IF(FSKWRITE(ISOCKET,SK_LONG_INT,2,CONNID2).EQ.-1)
     '            GOTO 9999
                FIELD_NAME='extracellular_potential'
                CALL STRING_TRIM(FIELD_NAME,IBEG,IEND)
                CLEN=FSKLEN(FIELD_NAME(IBEG:IEND))
                CALL FSKF2C(FIELD_NAME(IBEG:IEND),CLEN,INTSTR)
                IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
                IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
     '            GOTO 9999
                ISOCKET(1)=3
                ISOCKET(2)=1
                ISOCKET(3)=1
                IF(FSKWRITE(ISOCKET,SK_LONG_INT,3,CONNID2).EQ.-1)
     '            GOTO 9999
                CALL STRING_TRIM(FIELD_NAME,IBEG,IEND)
                CLEN=FSKLEN(FIELD_NAME(IBEG:IEND))
                CALL FSKF2C(FIELD_NAME(IBEG:IEND),CLEN,INTSTR)
                IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
                IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
     '            GOTO 9999
                ISOCKET(1)=5
                ISOCKET(2)=0
                IF(FSKWRITE(ISOCKET,SK_LONG_INT,2,CONNID2).EQ.-1)
     '            GOTO 9999
                FIELD_NAME='recovery'
                CALL STRING_TRIM(FIELD_NAME,IBEG,IEND)
                CLEN=FSKLEN(FIELD_NAME(IBEG:IEND))
                CALL FSKF2C(FIELD_NAME(IBEG:IEND),CLEN,INTSTR)
                IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
                IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
     '            GOTO 9999
                ISOCKET(1)=3
                ISOCKET(2)=1
                ISOCKET(3)=1
                IF(FSKWRITE(ISOCKET,SK_LONG_INT,3,CONNID2).EQ.-1)
     '            GOTO 9999
                CALL STRING_TRIM(FIELD_NAME,IBEG,IEND)
                CLEN=FSKLEN(FIELD_NAME(IBEG:IEND))
                CALL FSKF2C(FIELD_NAME(IBEG:IEND),CLEN,INTSTR)
                IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
                IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
     '            GOTO 9999
                ISOCKET(1)=6
                ISOCKET(2)=0
                IF(FSKWRITE(ISOCKET,SK_LONG_INT,2,CONNID2).EQ.-1)
     '            GOTO 9999
                FIELD_NAME='activation_time'
                CALL STRING_TRIM(FIELD_NAME,IBEG,IEND)
                CLEN=FSKLEN(FIELD_NAME(IBEG:IEND))
                CALL FSKF2C(FIELD_NAME(IBEG:IEND),CLEN,INTSTR)
                IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
                IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
     '            GOTO 9999
                ISOCKET(1)=3
                ISOCKET(2)=1
                ISOCKET(3)=1
                IF(FSKWRITE(ISOCKET,SK_LONG_INT,3,CONNID2).EQ.-1)
     '            GOTO 9999
                CALL STRING_TRIM(FIELD_NAME,IBEG,IEND)
                CLEN=FSKLEN(FIELD_NAME(IBEG:IEND))
                CALL FSKF2C(FIELD_NAME(IBEG:IEND),CLEN,INTSTR)
                IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
                IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
     '            GOTO 9999
                ISOCKET(1)=7
                ISOCKET(2)=0
                IF(FSKWRITE(ISOCKET,SK_LONG_INT,2,CONNID2).EQ.-1)
     '            GOTO 9999
              ENDIF
            ENDIF
            index=2
            ne_num_points=-1
            DO NOLIST=1,LIST(0)
C KAT 15Feb99: point_type != data
C              IF(POINT_TYPE(1:4).EQ.'DATA') THEN
C                nd=LIST(NOLIST)
C                IF(DATAFILE) THEN
C                  WRITE(IFILE,'(1X,''Node: '',I7)') nd+OFFSET
C                  WRITE(IFILE,'(1X,3E13.5)')  (ZD(nj,nd),nj=1,3)
C                ELSE
C                  IF(FSKWRITE(ICHAR('N'),SK_CHAR,1,CONNID2).EQ.-1)
C     '              GOTO 9999
C                  ISOCKET(1)=nd+10000
C                  IF(FSKWRITE(ISOCKET,SK_LONG_INT,1,CONNID2).EQ.-1)
C     '              GOTO 9999
C                  IF(FSKWRITE(ZD(1,nd),SK_DOUBLE_FLOAT,3,CONNID2).EQ.-1)
C     '              GOTO 9999
C                ENDIF
C                IF(ERR) THEN    !write residuals
C                  !ERROR code taken from OPDATA
C                  ne=LD(nd)     !element number for data point nd
C                  IF(ne.GT.0) THEN
C                    CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF,NPNE(1,1,ne),
C     '                NRE(ne),NVJE(1,1,1,ne),
C             '                SE(1,1,ne),XA(1,1,ne),XE,XP,
C             ERROR,*9999)
C                    IF((KTYP8.LE.1.OR.KTYP8.EQ.6).AND.ITYP6(1,1).LE.1)
C     '                THEN
C                      DO nj=1,NJ_LOC(NJL_GEOM,0,NRE(ne))
C                        nb=NBJ(nj,ne)
C                        X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
C     '                    nb,1,XID(1,nd),XE(1,nj))
C                      ENDDO
C                      CALL XZ(ITYP10(NRE(ne)),X,Z) !trans coords to r.c.
C                      DO nj=1,NJ_LOC(NJL_GEOM,0,NRE(ne))
C                        Z(nj)=Z(nj)-ZD(nj,nd)
C                      ENDDO
C                    ELSE
C                      ERROR='>>Not implemented'
C                      GOTO 9999
C                    ENDIF
C                  ELSE
C                    DO nj=1,3
C                      Z(nj)=0.0D0
C                    ENDDO
C                  ENDIF
C                  IF(DATAFILE) THEN
C                    WRITE(IFILE,'(1X,3E13.5)')  (Z(nj),nj=1,3)
C                  ELSE
C                    IF(FSKWRITE(Z,SK_DOUBLE_FLOAT,3,
C     '                CONNID2).EQ.-1)
C     '                GOTO 9999
C                  ENDIF
C                ENDIF
C              ELSE
              IF(POINT_TYPE(1:5).EQ.'GAUSS') THEN
                IF(GROUPED) THEN
                  index=index+ne_num_points+1
                  ne=ILISTMBR(%VAL(LIGRGA_PTR(nogrga)),index)
                  index=index+1
                  ne_num_points=ILISTMBR(%VAL(LIGRGA_PTR(nogrga)),index)
                ELSE
                  ne=LIST(NOLIST)
                  nb=NBJ(1,ne)
                  ne_num_points=NGT(nb)
                ENDIF
                CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '            NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '            SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
                nb=NBJ(1,ne)
C GMH 27/10/95 All this is taken from upgaus and opst40
                IF(STRESS) THEN
                  CALL ZPZE(NBH(1,1,ne),nc,NHE(ne,nx),NKHE(1,1,1,ne),
     '              NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),
     '              NW(ne,1,nx),nx,CURVCORRECT(1,1,1,ne),
     '              SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
                  CALL CPCG(NW(ne,1,nx),nb,NPNE(1,1,ne),nr,nx,
     '              CE(1,ne,nx),CG,CGE(1,1,ne,nx),
     '              CP(1,1,nx),PG,ERROR,*9999)
                ENDIF
                DO nong=1,ne_num_points
                  IF(GROUPED) THEN
                    ng=ILISTMBR(%VAL(LIGRGA_PTR(nogrga)),index+nong)
                  ELSE
                    ng=nong
                  ENDIF
                  CALL XEXG(NBJ(1,ne),ng,NRE(ne),PG,
     '              XE,XG,ERROR,*9999)
                  IF(DATAFILE) THEN
                    WRITE(IFILE,'(1X,''Node:'',I7)')
     '                ng+(NOLIST-1)*NGT(nb)+OFFSET
                    WRITE(IFILE,'(1X,3E13.5)') (XG(nj,1),nj=1,NJT)
                  ELSE
                    IF(FSKWRITE(ICHAR('N'),SK_CHAR,1,CONNID2).EQ.-1)
     '                GOTO 9999
                    ISOCKET(1)=ng+(NOLIST-1)*NGT(nb)+10000
                    IF(FSKWRITE(ISOCKET,SK_LONG_INT,1,CONNID2).EQ.-1)
     '                GOTO 9999
                    IF(FSKWRITE(XG(1,1),SK_DOUBLE_FLOAT,3,
     '                CONNID2).EQ.-1) GOTO 9999
                  ENDIF
C GMH 27/10/95 All this is taken from upgaus and opst40
                  IF(STRESS) THEN
                    CALL XGMG(1,NIT(NBJ(1,ne)),NBJ(1,ne),NRE(ne),
     '                DXINU,GL,GU,RG,XG,
     '                ERROR,*9999)
                    CALL ZEZG(1,NBH(1,nc,ne),ng,NHE(ne,nx),nx,DXINU,
     '                PG,ZE,ZG,ERROR,*9999)
                    CALL CGS9(NRE(ne),NW(ne,1,nx),nx,CG(1,ng),EG,RT,
     '                TG,ZG,ERROR,*9999)
                    ! now actually write TG
                    IF(DATAFILE) THEN
                      WRITE(IFILE,'(1X,3E13.5)')  (TG(1,nj),nj=1,3)
                      WRITE(IFILE,'(1X,3E13.5)')  (TG(2,nj),nj=1,3)
                      WRITE(IFILE,'(1X,3E13.5)')  (TG(3,nj),nj=1,3)
                    ELSE
                      IF(FSKWRITE(TG(1,1),SK_DOUBLE_FLOAT,3,
     '                  CONNID2).EQ.-1)
     '                  GOTO 9999
                      IF(FSKWRITE(TG(2,1),SK_DOUBLE_FLOAT,3,
     '                  CONNID2).EQ.-1)
     '                  GOTO 9999
                      IF(FSKWRITE(TG(3,1),SK_DOUBLE_FLOAT,3,
     '                  CONNID2).EQ.-1)
     '                  GOTO 9999
                    ENDIF
                  ENDIF
                ENDDO

              ELSE IF(POINT_TYPE(1:4).EQ.'GRID') THEN
                nq=LIST(NOLIST)
                IF(DATAFILE) THEN
                  IF((ITYP3(nr,nx).EQ.1).AND.(ITYP2(nr,nx).EQ.5)) THEN
                    !flow in elastic tubes
                    ny_p=NYNQ(1,nq,0,nx)
                    ny_r=NYNQ(2,nq,0,nx)
                    ny_v=NYNQ(3,nq,0,nx)
                    PRESS_FIELD=YQ(ny_p,1,1,nx)
                    RADIUS_FIELD=YQ(ny_r,1,1,nx)

C PM 26-JUL-01 : to accommodate for negative velocities
C                    VEL_FIELD=MAX(YQ(ny_v,1,1,nx),LOOSE_TOL)
                    VEL_FIELD=YQ(ny_v,1,1,nx)
                    IF(PATH) THEN
                      PATH_FIELD=MAX(YQ(ny_r,6,1,nx),LOOSE_TOL)
                    ELSE
                      PATH_FIELD=0.0d0
                    ENDIF
                    POT_FIELD=0.0d0
                    ACT_FIELD=0.0d0
                  ELSE
                    CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqV,NIQ_V,
     '                ERROR,*9999)
                    CALL MAQ_LOC(MAQ_INQUIRE,MAQ_TIME,atime,
     '                MAQ_ACTIV_TIME,ERROR,*9999)
                    IF(niqV.EQ.0) niqV=1
                    IF(atime.EQ.0) atime=1
                    ACT_FIELD=AQ(atime,nq)
                    POT_FIELD=YQ(nq,niqV,1,nx)
                    PRESS_FIELD=0.0d0
                    RADIUS_FIELD=0.5d0
                    VEL_FIELD=0.0d0
                    PATH_FIELD=0.0d0
                  ENDIF

                  WRITE(IFILE,'(1X,''Node: '',I7)') nq+OFFSET
                  NUM_FIELD_LIST=0
                  DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                    NUM_FIELD_LIST=NUM_FIELD_LIST+1
                    FIELD_LIST(NUM_FIELD_LIST)=XQ(nj,nq)
                  ENDDO
                  IF(ACTIVATION_TIME) THEN
                    NUM_FIELD_LIST=NUM_FIELD_LIST+1
                    FIELD_LIST(NUM_FIELD_LIST)=ACT_FIELD
                  ENDIF
                  IF(POTENTIAL) THEN
                    NUM_FIELD_LIST=NUM_FIELD_LIST+1
                    FIELD_LIST(NUM_FIELD_LIST)=POT_FIELD
                  ENDIF
                  IF(PRESSURE) THEN
                    NUM_FIELD_LIST=NUM_FIELD_LIST+1
                    FIELD_LIST(NUM_FIELD_LIST)=PRESS_FIELD
                  ENDIF
                  IF(RADIUS) THEN
                    NUM_FIELD_LIST=NUM_FIELD_LIST+1
                    FIELD_LIST(NUM_FIELD_LIST)=RADIUS_FIELD
                  ENDIF
                  IF(VELOCITY) THEN
                    NUM_FIELD_LIST=NUM_FIELD_LIST+1
                    FIELD_LIST(NUM_FIELD_LIST)=VEL_FIELD
                  ENDIF
                  IF(PATH) THEN
                    NUM_FIELD_LIST=NUM_FIELD_LIST+1
                    FIELD_LIST(NUM_FIELD_LIST)=PATH_FIELD
                  ENDIF

C PM 26-JUL-01 : fluid in elastic tube
                  IF((ITYP3(nr,nx).EQ.1).AND.(ITYP2(nr,nx).EQ.5).AND.
     '               (.NOT.PRESSURE).AND.(.NOT.RADIUS).AND.
     '               (.NOT.VELOCITY)) THEN
                    NUM_FIELD_LIST=NUM_FIELD_LIST+3
                    FIELD_LIST(NUM_FIELD_LIST-2)=PRESS_FIELD
                    FIELD_LIST(NUM_FIELD_LIST-1)=RADIUS_FIELD
                    FIELD_LIST(NUM_FIELD_LIST)=VEL_FIELD
                  ENDIF

                  WRITE(IFILE,'(1X,9E13.5)')
     '              (FIELD_LIST(nj),nj=1,NUM_FIELD_LIST)

                  IF(MATERIAL) THEN
C *** DPN 13 September 1999 - Adding cellular material parameters
                    !The cell type field
                    WRITE(IFILE,'(1X,I2)') ICQS_SPATIAL(1,nq)
                    !State variables
                    DO j=1,CELL_NUM_STATE(1)
                      WRITE(IFILE,'(1X,9E13.5)')
     '                  YQS(CELL_STATE_OFFSET(1)-1+j,nq)
                    ENDDO !j
                    !The integer parameters
                    DO j=2,NQIT !skip variant number
                      CELL_FIELD=.FALSE.
                      DO k=1,IICQS_SPATIAL(0,1)
                        IF(IICQS_SPATIAL(k,1).EQ.j) THEN
                          CELL_FIELD=.TRUE.
                          nqsv=k
                        ENDIF
                      ENDDO !k
                      IF(CELL_FIELD) THEN
                        !normal node field
                        WRITE(CHAR,'(I12)') ICQS_SPATIAL(nqsv,nq)
                        CALL STRING_TRIM(CHAR,IBEG,IEND)
                        WRITE(IFILE,'(1X,A)') CHAR(IBEG:IEND)
                      ELSE
                        !indexed field
                        ! ** do nothing **
                      ENDIF
                    ENDDO !j
                    !The real parameters
                    DO j=1,NQRT
                      CELL_FIELD=.FALSE.
                      DO k=1,IRCQS_SPATIAL(0,1)
                        IF(IRCQS_SPATIAL(k,1).EQ.j) THEN
                          CELL_FIELD=.TRUE.
                          nqsv=k
                        ENDIF
                      ENDDO !k
                      IF(CELL_FIELD) THEN
                        !normal node field
                        WRITE(IFILE,'(1X,9E13.5)') RCQS_SPATIAL(nqsv,nq)
                      ELSE
                        !indexed field
                        ! ** do nothing **
                      ENDIF
                    ENDDO !j
                  ENDIF !MATERIAL

                ELSE
                  IF(FSKWRITE(ICHAR('N'),SK_CHAR,1,CONNID2).EQ.-1)
     '              GOTO 9999
                  ISOCKET(1)=nq+10000
                  IF(FSKWRITE(ISOCKET,SK_LONG_INT,1,CONNID2).EQ.-1)
     '              GOTO 9999
                  IF(FSKWRITE(XQ(1,nq),SK_DOUBLE_FLOAT,3,CONNID2).EQ.-1)
     '              GOTO 9999
                  IF(FSKWRITE(YQ(nq,1,1,nx),SK_DOUBLE_FLOAT,1,
     '              CONNID2).EQ.-1)
     '              GOTO 9999
                  IF(FSKWRITE(YQ(nq,2,1,nx),SK_DOUBLE_FLOAT,1,
     '              CONNID2).EQ.-1)
     '              GOTO 9999
                  IF(FSKWRITE(YQ(nq,3,1,nx),SK_DOUBLE_FLOAT,1,
     '              CONNID2).EQ.-1)
     '              GOTO 9999
                  IF(FSKWRITE(YQ(nq,1,3,nx),SK_DOUBLE_FLOAT,1,
     '              CONNID2).EQ.-1)
     '              GOTO 9999
                ENDIF
              ENDIF
            ENDDO
            IF(DATAFILE) THEN
              CALL CLOSEF(IFILE,ERROR,*9999)
            ELSE
              IF(FSKWRITE(0,SK_CHAR,1,CONNID2).EQ.-1) GOTO 9999
            ENDIF
          ENDIF

        ELSE IF(OUTPUT(1:8).EQ.'EXPLORER') THEN

        ENDIF

      ENDIF

      CALL EXITS('EXPOIN')
      RETURN
 9999 CALL ERRORS('EXPOIN',ERROR)
      CALL EXITS('EXPOIN')
      RETURN 1
      END


