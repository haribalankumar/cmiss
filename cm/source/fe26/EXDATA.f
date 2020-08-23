      SUBROUTINE EXDATA(IBT,IDO,INP,LD,LIST,NBJ,NDDATA,NDP,NENQ,NKJE,
     &  NPF,NPNE,NRE,NRLIST,NVJE,NWQ,NXQ,AQ,DXDXIQ,DXDXIQ2,SE,XA,XE,
     &  XID,XIQ,XP,Z_CONT,ZD,STRING,ERROR,*)

C#### Subroutine: EXDATA
C###  Description:
C###    EXDATA exports data points from finite element data base.

C**** NOTE: for socket connection CONNID2 is defined as a parameter
C**** specifying the data transfer socket number
C**** 22/6/94 Now conforms to .EXNODE file structure (GBS)
C?????DB.  This needs re-thinking.  Maybe we need several numbers for
C?????  each cmgui node (like the element, face and line numbers for
C?????  the cmgui elements
C**** 22/04/96 GBS Added export of potentials and time at gridpoints
C**** CS 14/7/98 made this EXDATA routine from EXPOIN

      IMPLICIT none
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'      
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cmis00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'fsklib.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'loc00.cmn'
      INCLUDE 'maqloc00.inc'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LD(NDM),LIST(0:NLISTM),NBJ(NJM,NEM),NDDATA(0:NDM,0:NRM),
     '  NDP(NDM),NENQ(0:8,NQM),
     '  NKJE(NKM,NNM,NJM,NEM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NRE(NEM),NRLIST(0:NRM),NVJE(NNM,NBFM,NJM,NEM),NWQ(8,0:NQM,NAM),
     '  NXQ(-NIM:NIM,0:4,0:NQM)
      REAL*8 AQ(NMAQM,NQM),DXDXIQ(3,3,NQM),DXDXIQ2(3,3,NQM),
     &  SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),XE(NSM,NJM),XID(NIM,NDM),
     &  XIQ(NIM,NQM),XP(NKM,NVM,NJM,NPM),Z_CONT(NDM,2,67),ZD(NJM,NDM)
      CHARACTER ERROR*(*),STRING*(MXCH)

!     Local Variables
      INTEGER MX_FIELDS
      PARAMETER (MX_FIELDS=9)

      INTEGER CLEN,DATA_TYPE,IBEG,IBEG1,IBEG2,IBEG3,
     '  IBEG4,IEND,IEND1,IEND2,IEND3,IEND4,IFROMC,INTSTR(1024),
     '  ISOCKET(3),maq,N3CO,nb,nd,ne,NFIELDS(0:MX_FIELDS),ni,
     '  NITB,nj,NOLIST,nq,nr,OFFSET,WARN_NDP
      REAL*8 PXI,X(6),XNLOCAL(3),Z(6)
      CHARACTER CHAR1*5,CHAR2*5,CHAR3*5,CMISS_NUMBER*3,
     '  COMPONENT_NAME*80,FIELD_NAME*80,
     '  FILE*200,OUTPUT*11,POINT_NAME*50
      LOGICAL ALL_REGIONS,CBBREV,CONTACT,DATAFILE,ERR,LATTICE,
     '  MASTER,NORMALS
      PARAMETER(DATA_TYPE = 4)  !Point Data


      CALL ENTERS('EXDATA',*9999)


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

C#### Command: FEM export data<;FILENAME[default]>
C###  Parameter:      <as NAME[10001..10000]>
C###    Name the data set.
C###  Parameter:      <number LIST[all]>
C###    Defines the data to include in the exported data set.
C new CS 15/11/99 added export using nd for label
C###  Parameter:      <cmiss_number (NDP/nd)[NDP]>
C###    Numbers the data points using the original data point
C###    numbers NDP, or the data point loop counter nd. The later is
C###    useful if there are multiple data points with the same
C###    original number.
C###  Parameter:      <error>
C###    Include with the data a field named `error' containing the
C###    projections to elements.
C###  Parameter:      <contact>
C###    Include with the data fields associated with contact mechanics:
C###    frictionless and frictional contact pressures, contact gap,
C###    slip functions
C###  Parameter:      <master>
C###    Include with the data fields of normal and tangential vectors
C###    defined at the master face for contact mechanics problems
C###  Parameter:      <to (datafile/Motif)[datafile]>
C###    Specify the destination file format (`datafile' for exporting to
C###    CMGUI).
C###  Parameter:      <num_fields #)[0]>
C###    Export field components associated with the data
C###  Parameter:      <field_name NAME)[UnknownField]>
C###    Name of the field to be exported
C###  Parameter:      <offset OFFSET[10000]>
C###    Add OFFSET to data points' reference numbers when
C###    calculating their CMGUI node numbers.
C###  Parameter:      <lattice>
C###    Export the lattice grid points.         
C###  Parameter:      <normals>
C###    Export the surface normals at grid points as vectors.        
C###  Parameter:      <region (#s/all)[1]>
C###    Limit to data points of specified regions.
C###  Description:
C###    Exports data to the CMGUI.
CC###    However to export data so it can be displayed by the graphics
CC###    interface use the ipdata_to_exnode.awk file to display the data
CC###    points as nodes. Go to the cmiss utils page from the programmers
CC###    help to get details on how to use this file.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME['//FILE00(IBEG1:IEND1)
     '    //']>'
        OP_STRING(2)=BLANK(1:15)//'<as NAME['//POINT_NAME(IBEG4:IEND4)
     '    //']>'
        OP_STRING(3)=BLANK(1:15)
     '    //'<number LIST[all]>'
        OP_STRING(4)=BLANK(1:15)
     '    //'<cmiss_number (NDP/nd)[NDP]>'
        OP_STRING(5)=BLANK(1:15)//'<error>'
        OP_STRING(6)=BLANK(1:15)//'<contact>'
        OP_STRING(7)=BLANK(1:15)//'<master>'
        OP_STRING(8)=BLANK(1:15)
     '    //'<to (datafile/Motif)[datafile]>'
        OP_STRING(9)=BLANK(1:15)//'<field_nums #s)[0]>'
        OP_STRING(10)=BLANK(1:15)//'<field_name NAME)[UnknownField]>'
        OP_STRING(11)=BLANK(1:15)//'<offset OFFSET[10000]>'
        OP_STRING(12)=BLANK(1:15)//'<lattice>'
        OP_STRING(13)=BLANK(1:15)//'<normals>'
        OP_STRING(14)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C KAT 15Feb99: already done
CC GMH 22/10/95 Repeat for each case to name groups correctly
C        WRITE(CHAR2,'(I5)') 1+10000
C        WRITE(CHAR3,'(I5)') NET(1)+10000
C        CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
C        CALL STRING_TRIM(CHAR3,IBEG3,IEND3)
C        POINT_NAME=CHAR2(IBEG2:IEND2)//'..'//CHAR3(IBEG3:IEND3)
C        CALL STRING_TRIM(POINT_NAME,IBEG4,IEND4)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe26','doc','EXPOIN',ERROR,*9999)
      ELSE
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        nr=NRLIST(1)

C LKC 18-AUG-1998
C     Note : Remove offset once CMGUI data is upto date
C news MPN 9-11-95: OFFSET option to add const to all points

        IF(CBBREV(CO,'OFFSET',2,noco+1,NTCO,N3CO)) THEN
          OFFSET=IFROMC(CO(N3CO+1))
        ELSE
          OFFSET=10000
        ENDIF


C KAT 15Feb99: Done later anyway
C        WRITE(CHAR2,'(I5)') 1+OFFSET
C        WRITE(CHAR3,'(I5)') NQT+OFFSET
C        CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
C        CALL STRING_TRIM(CHAR3,IBEG3,IEND3)
C        POINT_NAME=CHAR2(IBEG2:IEND2)//'..'//CHAR3(IBEG3:IEND3)
C        CALL STRING_TRIM(POINT_NAME,IBEG4,IEND4)

C new CS 15/11/99 added export using nd for label
        IF(CBBREV(CO,'CMISS_NUMBER',2,noco+1,NTCO,N3CO)) THEN
          IF(CBBREV(CO,'NDP',3,noco+1,NTCO,N3CO)) THEN
            CMISS_NUMBER='NDP'
          ELSE IF(CBBREV(CO,'ND',2,noco+1,NTCO,N3CO)) THEN
            CMISS_NUMBER='ND'
          ENDIF
        ELSE
          CMISS_NUMBER='NDP'
        ENDIF

C new GMH 19/9/95 Adding export error to data point
        IF(CBBREV(CO,'ERROR',3,noco+1,NTCO,N3CO)) THEN
          ERR=.TRUE.
        ELSE
          ERR=.FALSE.
        ENDIF

C NEWS JHC 19/02/08 Adding export contact pressure to data point
        IF(CBBREV(CO,'CONTACT',3,noco+1,NTCO,N3CO)) THEN
          CONTACT=.TRUE.
        ELSE
          CONTACT=.FALSE.
        ENDIF
C NEWE
C*** 27/05/08 JHC Exporting normal and tangent vectors defined at master face
        IF(CBBREV(CO,'MASTER',4,noco+1,NTCO,N3CO)) THEN
          MASTER=.TRUE.
        ELSE
          MASTER=.FALSE.
        ENDIF
        IF(CBBREV(CO,'FIELD_NUMS',7,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),MX_FIELDS,NFIELDS(0),NFIELDS(1),
     '      ERROR,*9999)
          DO nolist=1,NFIELDS(0)
            nj=NFIELDS(nolist)
            IF(nj.GT.NJ_LOC(NJL_FIEL,0,nr)) THEN
              WRITE(ERROR,
     '          '(''Field '',I1,'' not defined for export'')') nj
              GOTO 9999
            ENDIF
          ENDDO !nonlist
        ELSE
          NFIELDS(0)=0
          NFIELDS(1)=0
        ENDIF


        IF(CBBREV(CO,'TO',2,noco+1,NTCO,N3CO)) THEN
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

C DAH 24-MAY-2002
C     Adding in lattice grid point export

        IF(CBBREV(CO,'LATTICE',4,noco+1,NTCO,N3CO)) THEN
          LATTICE=.TRUE.
        ELSE
          LATTICE=.FALSE.
        ENDIF

C DAH 21-03-02 Adding in export of surface normals at gridpoints. 

        IF(CBBREV(CO,'NORMALS',4,noco+1,NTCO,N3CO)) THEN
          NORMALS=.TRUE.
        ELSE
          NORMALS=.FALSE.
        ENDIF

C Check that only valid options are specified for NORMALS and LATTICE
        IF(LATTICE.AND.NORMALS) THEN
          ERROR=
     '      '>> Normals and lattice points cannot be exported together'
          GOTO 9999
        ENDIF
        
        IF(LATTICE.OR.NORMALS) THEN
          CALL ASSERT(.NOT.ERR,
     '      '>> ERROR is not valid for NORMALS or LATTICE',ERROR,*9999)
        ENDIF

C KAT 15Feb99: ELEMENT not used
        IF(CBBREV(CO,'NUMBER' ,2,noco+1,NTCO,N3CO)) THEN
C        IF(CBBREV(CO,'NUMBER' ,2,noco+1,NTCO,N3CO).OR
C     '    .CBBREV(CO,'ELEMENT',2,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NLISTM,LIST(0),LIST(1),ERROR,*9999)
C KAT 15Feb99: GROUP not used. Use AS
C          WRITE(CHAR1,'(I5)') LIST(1)+OFFSET
C          WRITE(CHAR2,'(I5)') LIST(LIST(0))+OFFSET
C          CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
C          CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
C          POINT_NAME=CHAR1(IBEG1:IEND1)//'..'//CHAR2(IBEG2:IEND2)
C        ELSE IF(CBBREV(CO,'GROUP',2,noco+1,NTCO,N3CO)) THEN
C          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
C         POINT_NAME=CO(N3CO+1)(IBEG:IEND)

        ELSEIF(.NOT.LATTICE.AND..NOT.NORMALS) THEN !need this check

C LKC 21-MAY-1999 better to use NDDATA instead of NDT
C          LIST(0)=NDT
C          DO nd=1,NDT
C            LIST(nd)=nd
C
C awaiting data region dependency

C KFA 06-DEC-2001 why not output empty files if thats whats asked for?
C          CALL ASSERT(NDDATA(0,nr).GT.0,
C     '      '>> No Data points for this region',ERROR,*9999)
          LIST(0)=NDDATA(0,nr)
          DO nd=1,NDDATA(0,nr)
            LIST(nd)=NDDATA(nd,nr)
          ENDDO


C KAT 15Feb99: GROUP no longer used. Use AS
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
C KAT 15Feb99: GROUP no longer used. Use AS
C        ENDIF

        IF(CBBREV(CO,'FIELD_NAME',7,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          FIELD_NAME=CO(N3CO+1)(IBEG1:IEND1)
        ELSE
          FIELD_NAME='UnknownField'
        ENDIF

        
        IF((OUTPUT(1:8).EQ.'DATAFILE').OR.(OUTPUT(1:5).EQ.'MOTIF')) THEN
          CALL CHECKF(1,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

C LKC 20-JUL-2005 should  this have been checked in with the lattice code?
          IF(.NOT.LATTICE.AND..NOT.NORMALS) THEN
C          IF(LIST(0).GT.0) THEN

            IF(DATAFILE) THEN
              CALL STRING_TRIM(FILE,IBEG,IEND)
              CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.exdata','NEW',
     '        'SEQUEN','FORMATTED',132,ERROR,*9999)
C  ####################################################################
            IF(FILE(IBEG:12).EQ.'./DATAerror1') THEN
              open (unit=12,
     $          file=FILE(IBEG:IEND)//'.bin',
     $          form='unformatted',
     $          access='direct',
     $          recl=NJ_LOC(NJL_GEOM,0,nr)*8,
     $          ERR=9999)
            END IF
C  ####################################################################
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
            ENDIF

            IF(ERR) THEN
              IF(DATAFILE) THEN
                WRITE(IFILE,'(1X,''#Fields=2'')')
              ELSE
                ISOCKET(1)=2
                IF(FSKWRITE(ISOCKET,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
              ENDIF
C NEWS JHC 19/02/08 added for contact pressure
            ELSEIF(CONTACT) THEN
              IF(DATAFILE) THEN
                WRITE(IFILE,'(1X,''#Fields=6'')')
              ENDIF
C NEWE JHC
C*** 27/05/08 normal and tangent vectors calculated on master face
            ELSEIF(MASTER) THEN
              IF(DATAFILE) THEN
                WRITE(IFILE,'(1X,''#Fields=4'')')
              ENDIF
            ELSEIF(NFIELDS(0).GT.0) THEN
              IF(DATAFILE) THEN
                WRITE(IFILE,'(1X,''#Fields=2'')')
              ELSE
                ISOCKET(1)=2
                IF(FSKWRITE(ISOCKET,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
              ENDIF
            ELSE
              IF(DATAFILE) THEN
                WRITE(IFILE,'(1X,''#Fields=1'')')
              ELSE
                IF(FSKWRITE(ICHAR('#'),SK_CHAR,1,CONNID2).EQ.-1)
     '            GOTO 9999
              ENDIF
            ENDIF
            IF(DATAFILE) THEN

C LKC 4-SEP-2000 Generalise for 1d/2d/3d
              IF(NJ_LOC(NJL_GEOM,0,nr).EQ.1) THEN
                WRITE(IFILE,'(1X,''1) coordinates, coordinate, '
     '            //'rectangular cartesian, #Components=1'')')
                WRITE(IFILE,'(1X,''  x.  Value index= 1, '
     '          //'#Derivatives=0'')')
              ELSEIF(NJ_LOC(NJL_GEOM,0,nr).EQ.2) THEN
                WRITE(IFILE,'(1X,''1) coordinates, coordinate, '
     '            //'rectangular cartesian, #Components=2'')')
                WRITE(IFILE,'(1X,''  x.  Value index= 1, '
     '          //'#Derivatives=0'')')
                WRITE(IFILE,'(1X,''  y.  Value index= 2, '
     '            //'#Derivatives=0'')')
              ELSE
                WRITE(IFILE,'(1X,''1) coordinates, coordinate, '
     '            //'rectangular cartesian, #Components=3'')')
                WRITE(IFILE,'(1X,''  x.  Value index= 1, '
     '          //'#Derivatives=0'')')
                WRITE(IFILE,'(1X,''  y.  Value index= 2, '
     '            //'#Derivatives=0'')')
                WRITE(IFILE,'(1X,''  z.  Value index= 3, '
     '            //'#Derivatives=0'')')
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

            IF(ERR) THEN
              IF(DATAFILE) THEN
C LKC 4-SEP-2000 Generalise for 1d/2d/3d
                IF(NJ_LOC(NJL_GEOM,0,nr).EQ.1) THEN
                  WRITE(IFILE,'(1X,''2) error, field, '
     '              //'rectangular cartesian, #Components=1'')')
                  WRITE(IFILE,'(1X,''  x.  Value index= 2, '
     '              //'#Derivatives=0'')')
                ELSEIF(NJ_LOC(NJL_GEOM,0,nr).EQ.2) THEN
                  WRITE(IFILE,'(1X,''2) error, field, '
     '              //'rectangular cartesian, #Components=2'')')
                  WRITE(IFILE,'(1X,''  x.  Value index= 3, '
     '              //'#Derivatives=0'')')
                  WRITE(IFILE,'(1X,''  y.  Value index= 4, '
     '              //'#Derivatives=0'')')
                ELSE
                  WRITE(IFILE,'(1X,''2) error, field, '
     '              //'rectangular cartesian, #Components=3'')')
                  WRITE(IFILE,'(1X,''  x.  Value index= 4, '
     '            //'#Derivatives=0'')')
                  WRITE(IFILE,'(1X,''  y.  Value index= 5, '
     '              //'#Derivatives=0'')')
                  WRITE(IFILE,'(1X,''  z.  Value index= 6, '
     '              //'#Derivatives=0'')')
                ENDIF
              ELSE
                FIELD_NAME='error'
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
              ENDIF
C NEWS JHC 15/04/08 Added for contact pressure
C            Note that only works for 3D element
            ELSEIF(CONTACT) THEN
              IF(DATAFILE) THEN
                WRITE(IFILE,'(1X,''2) ContactGap, field, '
     &            //'rectangular cartesian, #Components=1'')')
                WRITE(IFILE,'(1X,''  ContactGap.  Value index= 4, '
     &            //'#Derivatives=0'')')
                WRITE(IFILE,'(1X,''3) ContactPressure, field, '
     &            //'rectangular cartesian, #Components=1'')')
                WRITE(IFILE,'(1X,''  ContactPressure.  Value index= 5, '
     &            //'#Derivatives=0'')')
                WRITE(IFILE,'(1X,''4) Frictional_1, field, '
     &            //'rectangular cartesian, #Components=1'')')
                WRITE(IFILE,'(1X,''  Frictional_1.  Value index= 6, '
     &            //'#Derivatives=0'')')
                WRITE(IFILE,'(1X,''5) Frictional_2, field, '
     &            //'rectangular cartesian, #Components=1'')')
                WRITE(IFILE,'(1X,''  Frictional_2.  Value index= 7, '
     &            //'#Derivatives=0'')')
                WRITE(IFILE,'(1X,''6) Slip_Func, field, '
     &            //'rectangular cartesian, #Components=1'')')
                WRITE(IFILE,'(1X,''  Slip_Func.  Value index= 8, '
     &            //'#Derivatives=0'')')
              ENDIF
C NEWE JHC
C*** 27/05/08 JHC normal and tangent vectors at master face
C            Note that only works for 3D element
            ELSEIF(MASTER) THEN
              IF(DATAFILE) THEN
                WRITE(IFILE,'(1X,''2) normal, field, '
     &            //'rectangular cartesian, #Components=3'')')
                WRITE(IFILE,'(1X,''  x.  Value index= 4, '
     &          //'#Derivatives=0'')')
                WRITE(IFILE,'(1X,''  y.  Value index= 5, '
     &            //'#Derivatives=0'')')
                WRITE(IFILE,'(1X,''  z.  Value index= 6, '
     &            //'#Derivatives=0'')')
                WRITE(IFILE,'(1X,''3) tangent_1, field, '
     &            //'rectangular cartesian, #Components=3'')')
                WRITE(IFILE,'(1X,''  x.  Value index= 7, '
     &          //'#Derivatives=0'')')
                WRITE(IFILE,'(1X,''  y.  Value index= 8, '
     &            //'#Derivatives=0'')')
                WRITE(IFILE,'(1X,''  z.  Value index= 9, '
     &            //'#Derivatives=0'')')
                WRITE(IFILE,'(1X,''4) tangent_2, field, '
     &            //'rectangular cartesian, #Components=3'')')
                WRITE(IFILE,'(1X,''  x.  Value index= 10, '
     &          //'#Derivatives=0'')')
                WRITE(IFILE,'(1X,''  y.  Value index= 11, '
     &            //'#Derivatives=0'')')
                WRITE(IFILE,'(1X,''  z.  Value index= 12, '
     &            //'#Derivatives=0'')')
              ENDIF
            ELSEIF(NFIELDS(0).GT.0) THEN

C LKC 7-MAR-2002
C Exporting 1 field with #components=NFIELDS along the NJT
C geometric positions for the data points. The components
C are INTS because they do not necessarily correspond to x/y/z

              CALL STRING_TRIM(FIELD_NAME,IBEG1,IEND1)
              WRITE(IFILE,'(1X,''2) '',A,'', field, '
     '          //'rectangular cartesian, #Components='',I1)')
     '          FIELD_NAME(IBEG1:IEND1),NFIELDS(0)
              DO nj=1,NFIELDS(0)
                WRITE(IFILE,'(1X,''  '',I1,''.  Value index= '',I1,
     '            '', #Derivatives=0'')') nj,NJ_LOC(NJL_FIEL,nj,nr)
              ENDDO

            ENDIF


            WARN_NDP=0
            DO NOLIST=1,LIST(0)

C LKC 18-AUG-1998 Change data numbering here to use NDP
C LKC 13-MAY-1999 Modification if NDP not setup
              nd=LIST(NOLIST)
              IF(DATAFILE) THEN
C A weak check to guess if NDP(nd) is unitialised.
                IF(NDP(nd).LE.0.OR.NDP(nd).GT.100000000) THEN
                  WARN_NDP=WARN_NDP+1
                  NDP(nd)=nolist
                ENDIF

                IF(CMISS_NUMBER.EQ."ND") THEN
                  WRITE(IFILE,'(1X,''Node: '',I9)') nd+OFFSET
                ELSE
                  WRITE(IFILE,'(1X,''Node: '',I9)') NDP(nd)+OFFSET
                ENDIF

C LKC 4-SEP-2000 Generalise for 1/2/3d
C                  WRITE(IFILE,'(1X,3E13.5)')  (ZD(nj,nd),nj=1,3)
                  WRITE(IFILE,'(1X,3E13.5)')
     '            (ZD(nj,nd),nj=1,NJ_LOC(NJL_GEOM,0,nr))
              ELSE
                IF(FSKWRITE(ICHAR('N'),SK_CHAR,1,CONNID2).EQ.-1)
     '            GOTO 9999
                ISOCKET(1)=nd+10000
                IF(FSKWRITE(ISOCKET,SK_LONG_INT,1,CONNID2).EQ.-1)
     '            GOTO 9999
                IF(FSKWRITE(ZD(1,nd),SK_DOUBLE_FLOAT,3,CONNID2).EQ.-1)
     '            GOTO 9999
              ENDIF

              IF(ERR) THEN !write residuals
C                           ERROR code taken from OPDATA
                ne=LD(nd) !element number for data point nd
                IF(ne.GT.0) THEN
                  CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF,NPNE(1,1,ne),
     '              NRE(ne),NVJE(1,1,1,ne),
     '              SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
C JHC (MPN) 23Apr2004: need this for nonlinear problems
                  IF(KTYP8.LE.1.OR.KTYP8.EQ.6.OR.KTYP8.EQ.4) THEN
C                  IF((KTYP8.LE.1.OR.KTYP8.EQ.6.OR.KTYP8.EQ.4)
C     '              .AND.ITYP6(1,1).LE.1)
C     '              THEN
                    DO nj=1,NJ_LOC(NJL_GEOM,0,NRE(ne))
                      nb=NBJ(nj,ne)
                      X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '                  nb,1,XID(1,nd),XE(1,nj))
                    ENDDO
                    CALL XZ(ITYP10(NRE(ne)),X,Z) !trans coords to r.c.
                    DO nj=1,NJ_LOC(NJL_GEOM,0,NRE(ne))
                      Z(nj)=Z(nj)-ZD(nj,nd)
                    ENDDO
                  ELSE
                    ERROR='>>This error export not implemented'
                    GOTO 9999
                  ENDIF
                ELSE
C LKC 4-SEP-2000 Generalise for 1d/2d/3d
C                  DO nj=1,3
C                    Z(nj)=0.0D0
C                  ENDDO
                  DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                    Z(nj)=0.0D0
                  ENDDO

                ENDIF
                IF(DATAFILE) THEN
C 23-MAR-2001 Generalise for 1d/2d/3d
C                  WRITE(IFILE,'(1X,3E13.5)')  (Z(nj),nj=1,3)
                  WRITE(IFILE,'(1X,3E25.15)')
     '              (Z(nj),nj=1,NJ_LOC(NJL_GEOM,0,nr))
C  ####################################################################
                    IF(FILE(IBEG:12).EQ.'./DATAerror1') THEN
                      IF (NJ_LOC(NJL_GEOM,0,nr).EQ.2) THEN
                        write(12,rec=NOLIST,ERR=9999)Z(1),Z(2)
                      ELSE
                        write(12,rec=NOLIST,ERR=9999)Z(1),Z(2),Z(3)
                      ENDIF
                    ENDIF
C  ####################################################################
                ELSE
                  IF(FSKWRITE(Z,SK_DOUBLE_FLOAT,3,
     '              CONNID2).EQ.-1)
     '              GOTO 9999
                ENDIF
C*** 12/03/08 JHC Added for contact pressure
              ELSEIF(CONTACT) THEN !write residuals
                IF(DATAFILE) THEN
                  WRITE(IFILE,'(1X,3E13.5)') Z_CONT(nd,1,13) ! contact gap
                  WRITE(IFILE,'(1X,3E13.5)') Z_CONT(nd,1,18) ! contact pressure
                  WRITE(IFILE,'(1X,3E13.5)') Z_CONT(nd,1,20) ! frictional (1) contact force
                  WRITE(IFILE,'(1X,3E13.5)') Z_CONT(nd,1,22) ! frictional (2) contact force
                  WRITE(IFILE,'(1X,3E13.5)') Z_CONT(nd,1,45) ! Phi_trial
                ELSE
                  IF(FSKWRITE(Z,SK_DOUBLE_FLOAT,3,
     '              CONNID2).EQ.-1)
     '              GOTO 9999
                ENDIF
C*** 27/05/08 JHC Added for normal and tangent vectors at master face
              ELSEIF(MASTER) THEN !write normal and tangent vectors
                IF(DATAFILE) THEN
                  WRITE(IFILE,'(1X,3E13.5)') (Z_CONT(nd,1,3+nj),nj=1,3) ! normal vector
                  WRITE(IFILE,'(1X,3E13.5)') (Z_CONT(nd,1,6+nj),nj=1,3) ! tangent_1 vector
                  WRITE(IFILE,'(1X,3E13.5)') (Z_CONT(nd,1,9+nj),nj=1,3) ! tangent_2 vector
                ELSE
                  IF(FSKWRITE(Z,SK_DOUBLE_FLOAT,3,
     '              CONNID2).EQ.-1)
     '              GOTO 9999
                ENDIF
              ELSEIF(NFIELDS(0).GT.0) THEN
                WRITE(IFILE,'(1X,3E13.5)')
     '            (ZD(NJ_LOC(NJL_FIEL,NFIELDS(nj),nr),nd),
     '            nj=1,NFIELDS(0))
              ENDIF
            ENDDO !NOLIST
C ###################################################
            WRITE(*,*) FILE(IBEG:14)
            WRITE(*,*) './DATAerror1'
            IF(FILE(IBEG:12).EQ.'./DATAerror1') THEN
                close(unit=12)
            ENDIF
C ###################################################
C LKC 22-MAY-1999 new warning
            IF(WARN_NDP.NE.0) THEN
              WRITE(OP_STRING,'('' >> Warning: NDP not setup, '','
     '          //'''altering '',I7,'' of '',I7,'
     '          //''' reference data num.'')') WARN_NDP,LIST(0)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF


            IF(DATAFILE) THEN
              CALL CLOSEF(IFILE,ERROR,*9999)
            ELSE
              IF(FSKWRITE(0,SK_CHAR,1,CONNID2).EQ.-1) GOTO 9999
            ENDIF
          ENDIF !not LATTICE/NORMAL


          IF(LATTICE) THEN
            CALL ASSERT(CALL_GRID,' Must define grid first',
     '        ERROR,*9999)
            CALL ASSERT(USE_LAT.EQ.1,
     '        ' Choose lattice grid points in ipgrid',ERROR,*9999) 
            IF(DATAFILE) THEN
              CALL STRING_TRIM(FILE,IBEG,IEND)              
              CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.exdata','NEW',
     '          'SEQUEN','FORMATTED',132,ERROR,*9999)
C LKC 6-MAY-2003 Not Used              
C              CONTINUE=.TRUE.
              nq=0
              !write the file header info
              CALL STRING_TRIM(POINT_NAME,IBEG,IEND)
              WRITE(IFILE,'( '' Group name: '',A)')
     '          POINT_NAME(IBEG:IEND)
              WRITE(IFILE,'(1X,''#Fields=1'')')
              WRITE(IFILE,'(1X,''1) element_xi, field, '
     '          //'element_xi, #Components=1'')')
              WRITE(IFILE,'(1X,''  1.  Value index= 1, '
     '          //'#Derivatives=0'')')     
              DO nq=1,NQT
C               Get ne for nq
                ne=NENQ(1,nq)
                nb=NBJ(1,ne)
                NITB=NIT(nb)
                WRITE(IFILE,'(1X,''Node: '',I9)') nq+OFFSET
                IF(NITB.EQ.1) THEN !1d
                  WRITE(IFILE,
     '              '(10X,''E '',I3,'' 1 '',1E13.5)') ne,XIQ(1,nq)
                ELSEIF(NITB.EQ.2) THEN !2d
                  WRITE(IFILE,
     '              '(10X,''E '',I3,'' 2 '',2E13.5)') 
     '              ne,(XIQ(ni,nq),ni=1,2)
                ELSE !3d
                  WRITE(IFILE,
     '              '(10X,''E '',I3,'' 3 '',3E13.5)') 
     '              ne,(XIQ(ni,nq),ni=1,3)
                ENDIF
              ENDDO
              CALL CLOSEF(IFILE,ERROR,*9999)              
            ELSE
              ERROR='Lattice data export only set up for datafiles'
              GOTO 9999
            ENDIF
          ENDIF

          IF(NORMALS) THEN
            CALL ASSERT(CALL_GRID,' Must define grid first',
     '        ERROR,*9999)
            IF(DATAFILE) THEN
              CALL STRING_TRIM(FILE,IBEG,IEND)              
              CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.exdata','NEW',
     '          'SEQUEN','FORMATTED',132,ERROR,*9999)
              CALL STRING_TRIM(POINT_NAME,IBEG,IEND)             
              WRITE(IFILE,'( '' Group name: '',A)')
     '          POINT_NAME(IBEG:IEND)
              WRITE(IFILE,'(1X,''#Fields=2'')')
              WRITE(IFILE,'(1X,''1) element_xi, field, '
     '          //'element_xi, #Components=1'')')
              WRITE(IFILE,'(1X,''  1.  Value index= 1, '
     '          //'#Derivatives=0'')')
              WRITE(IFILE,'(1X,''2) axis1, field, '
     '          //'rectangular cartesian, #Components=3'')')
              WRITE(IFILE,'(1X,''  x.  Value index= 2, '
     '          //'#Derivatives=0'')')
              WRITE(IFILE,'(1X,''  y.  Value index= 3, '
     '          //'#Derivatives=0'')')
              WRITE(IFILE,'(1X,''  z.  Value index= 4, '
     '          //'#Derivatives=0'')')              
              DO nq=1,NQT
                IF(NWQ(1,nq,1).NE.0) THEN
                  ne=NENQ(1,nq)
                  WRITE(IFILE,'(1X,''Node: '',I9)') nq+OFFSET
                  IF(NITB.EQ.1) THEN !1d
                    WRITE(OP_STRING,
     '                '(''>>1d normal export not implemented'')')
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)                   
                  ELSEIF(NITB.EQ.2) THEN !2d
                    WRITE(OP_STRING,
     '                '(''>>2d normal export not implemented'')')
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)                    
                  ELSE !3d
                    WRITE(IFILE,
     '                '(10X,''E '',I3,'' 3 '',3E13.5)') 
     '                ne,(XIQ(ni,nq),ni=1,3)
                    IF(USE_LAT.EQ.1) THEN !normals from lattice grid points 
C                      CALL NORM_LATTICE(3,nq,NWQ(1,nq,1),DXDXIQ,DXDXIQ2,
C     '                  XNLOCAL,ERROR,*9999)
                      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_COORD,maq,
     &                  MAQ_NORMAL_X,ERROR,*9999)
                      XNLOCAL(1)=AQ(maq,nq)
                      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_COORD,maq,
     &                  MAQ_NORMAL_Y,ERROR,*9999)
                      XNLOCAL(2)=AQ(maq,nq)
                      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_COORD,maq,
     &                  MAQ_NORMAL_Z,ERROR,*9999)
                      XNLOCAL(3)=AQ(maq,nq)
                      WRITE(IFILE,'(10X,3E13.5)') (XNLOCAL(ni),ni=1,3)
                    ELSE !normals from collocation grid points
                      CALL NORM31(3,nq,NXQ,DXDXIQ,DXDXIQ2,
     '                  XNLOCAL,ERROR,*9999)
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO
              CALL CLOSEF(IFILE,ERROR,*9999)
            ELSE
              ERROR='Lattice data export only set up for datafiles'
              GOTO 9999                   
            ENDIF
          ENDIF
          
        ELSE IF(OUTPUT(1:8).EQ.'EXPLORER') THEN
          ERROR='Explorer not implemented'
          GOTO 9999
        ENDIF
      ENDIF

      CALL EXITS('EXDATA')
      RETURN
 9999 CALL ERRORS('EXDATA',ERROR)
      CALL EXITS('EXDATA')
      RETURN 1
      END


