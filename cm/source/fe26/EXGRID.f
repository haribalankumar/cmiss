      SUBROUTINE EXGRID(NEELEM,NLATNE,NQLIST,NQNLAT,NQS,NQSCNB,NQXI,
     &  NRLIST,NWQ,NXQ,XQ,YQ,STRING,ERROR,*)

C#### Subroutine: EXGRID
C###  Description:
C###    EXGRID creates an exelem file using grid point
C###    connectivity where the nodes are created using
C###    export points grid. Elements are exported as (2/4/8)
C###    noded elements with linear basis functions.
C***  Written by Martin Buist, June 1998
C***  Adding FieldML - Greg Sands, Sept 2003

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'

!     Parameter List
      INTEGER NEELEM(0:NE_R_M,0:NRM),NLATNE(NEQM+1),NQLIST(0:NQM),
     &  NQNLAT(NEQM*NQEM),NQS(NEQM),NQSCNB(NQSCM),NQXI(0:NIM,NQSCM),
     &  NRLIST(0:NRM),NWQ(8,0:NQM,NAM),NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      REAL*8 XQ(NJM,NQM),YQ(NYQM,NIQM,NAM,NXM)
      CHARACTER ERROR*(*),STRING*(MXCH)

!     Local Variables
      INTEGER ELEMENT,FACE(24),IBEG,IBEG1,IBEG2,IBEG3,IBEG4,IEND,IEND1,
     &  IEND2,IEND3,IEND4,IFROMC,mq,mq2,mq3,mq4,mq5,mq6,mq7,mq8,N3CO,nb,
     &  ne,nf,ni,nj,nk,nl,nn,no_nrlist,nq,nqq,nr,NUMFIELDS,
     &  offset_elem,offset_node,offset_line,offset_face,SCHEME

      REAL*8 SCALEFACTOR
      CHARACTER BASES*54,CHAR1*5,CHAR2*5,CHAR3*5,CHAR4*20,ELEM_NAME*50,
     &  FILE*100
      LOGICAL ACTIVATION_TIME,ALL_REGIONS,CBBREV,ELASTIC_TUBE,FACE1,
     &  FACE2,FACE3,FIELDML,LIBMESH,PATH,POTENTIAL,PRESSURE,RADIUS,
     &  SURFACE,VELOCITY,XI1XI2,XI1XI3,XI2XI3

      DATA FACE/3,7,2,10,5,8,4,11,1,9,3,5,
     &          6,12,7,8,2,4,1,6,10,11,9,12/

      CALL ENTERS('EXGRID',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        WRITE(CHAR2,'(I5)') NQT
        WRITE(CHAR3,'(I5)') NQT
        CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
        CALL STRING_TRIM(CHAR3,IBEG3,IEND3)
        CHAR4=CHAR2(IBEG2:IEND2)//'..'//CHAR3(IBEG3:IEND3)
        CALL STRING_TRIM(CHAR4,IBEG4,IEND4)

C---------------------------------------------------------------------

C#### Command: FEM export grid<;FILENAME[default]>
C###  Description:
C###    This command calculates linear elements between grid points
C###    and outputs the information in an exelem file to be used
C###    along with an exnode file generated from export points in
C###    Cmgui.
C###  Parameter:      <as NAME[0...0]>
C###    Name the field.
C###  Parameter:      <fieldml>
C###    Export as a FieldML file.
C###  Parameter:      <offset_elem OFFSET[0]>
C####   Add OFFSET to element numbers.
C###  Parameter:      <offset_node OFFSET[offset_elem]>
C###    Add OFFSET to node numbers
C###  Parameter:      <offset_line OFFSET[offset_elem]>
C####   Add OFFSET to line numbers.
C###  Parameter:      <offset_face OFFSET[offset_elem]>
C####   Add OFFSET to face numbers.
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the regions.
C###  Parameter:      <surface>
C###    In 3D only make bilinear surface elements, not trilinear
C###    volume elements
C###  Parameter:      <activation_time>
C###    Export the activation time field for activation problems
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
C###  Parameter:      <libmesh>
C###    Only export highest-order elements for libmesh

        OP_STRING(1)=STRING(1:IEND)
     &    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
        OP_STRING(2)=BLANK(1:15)//'<as NAME['//CHAR4(IBEG4:IEND4)//']>'
        OP_STRING(3)=BLANK(1:15)//'<fieldml>'
        OP_STRING(4)=BLANK(1:15)//'<offset_elem OFFSET[0]>'
        OP_STRING(5)=BLANK(1:15)//'<offset_node OFFSET[offset_elem]>'
        OP_STRING(6)=BLANK(1:15)//'<offset_line OFFSET[offset_elem]>'
        OP_STRING(7)=BLANK(1:15)//'<offset_face OFFSET[offset_elem]>'
        OP_STRING(8)=BLANK(1:15)//'<region #[1]>'
        OP_STRING(9)=BLANK(1:15)//'<surface>'
        OP_STRING(10)=BLANK(1:15)//'<activation_time>'
        OP_STRING(11)=BLANK(1:15)//'<potential>'
        OP_STRING(12)=BLANK(1:15)//'<pressure>'
        OP_STRING(13)=BLANK(1:15)//'<radius>'
        OP_STRING(14)=BLANK(1:15)//'<velocity>'
        OP_STRING(15)=BLANK(1:15)//'<path>'
        OP_STRING(16)=BLANK(1:15)//'<libmesh>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe26','doc','EXGRID',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRID.EQ.1,' Set USE_GRID to 1 in ippara',
     &    ERROR,*9999)
        CALL ASSERT(CALL_GRID,' Must define grid first',ERROR,*9999)

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

        IF(CBBREV(CO,'FIELDML',1,noco+1,NTCO,N3CO)) THEN
          CALL ASSERT(USE_LAT.EQ.1,'FieldML export only coded for '
     &      //'lattice-based grids',ERROR,*9999)
          FIELDML=.TRUE.
        ELSE
          FIELDML=.FALSE.
        ENDIF

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

C PM 26-JUL-01
        IF(CBBREV(CO,'ELASTIC_TUBE',3,noco+1,NTCO,N3CO)) THEN
          ELASTIC_TUBE=.TRUE.
        ELSE
          ELASTIC_TUBE=.FALSE.
        ENDIF
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

        IF(CBBREV(CO,'SURFACE',3,noco+1,NTCO,N3CO)) THEN
          SURFACE=.TRUE.
        ELSE
          SURFACE=.FALSE.
        ENDIF

        IF(CBBREV(CO,'LIBMESH',3,noco+1,NTCO,N3CO)) THEN
          LIBMESH=.TRUE.
        ELSE
          LIBMESH=.FALSE.
        ENDIF

        IF(CBBREV(CO,'AS',1,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          ELEM_NAME=CO(N3CO+1)(IBEG1:IEND1)
        ELSE
          WRITE(CHAR1,'(I5)') offset_node+1
          WRITE(CHAR2,'(I5)') offset_node+NQT
          CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
          CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
          ELEM_NAME=CHAR1(IBEG1:IEND1)//'..'//CHAR2(IBEG2:IEND2)
        ENDIF

        CALL CHECKF(1,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
        IF(FIELDML) THEN
          CALL STRING_TRIM(FILE,IBEG,IEND)
          CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.fml','NEW',
     &      'SEQUEN','FORMATTED',132,ERROR,*9999)
          WRITE(IFILE,'(''<?xml version="1.0" '
     &      //'encoding="iso-8859-1"?>'')')
          WRITE(IFILE,'(''<!-- Generated by Cmiss(cm) -->'')')
          WRITE(IFILE,'(''<fieldml xmlns:fieldml='
     &      //'"http://www.physiome.org.nz/fieldml/0.1#"'
     &      //' xmlns="http://www.physiome.org.nz/fieldml/0.1#">'')')
          CALL EXGRID_FIELDML(NLATNE,NQNLAT,NQS,NQSCNB,NQXI,NRLIST,
     &      offset_elem,offset_node,XQ,YQ,
     &      ACTIVATION_TIME,ELASTIC_TUBE,PATH,POTENTIAL,PRESSURE,RADIUS,
     &      VELOCITY,ERROR,*9999)
          WRITE(IFILE,'(''</fieldml>'')')
          CALL CLOSEF(IFILE,ERROR,*9999)
        ELSEIF(NQT.GT.0) THEN !grid points defined
          CALL STRING_TRIM(FILE,IBEG,IEND)
          CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.exelem','NEW',
     &      'SEQUEN','FORMATTED',132,ERROR,*9999)

          CALL STRING_TRIM(ELEM_NAME,IBEG,IEND)
          WRITE(IFILE,'( '' Group name: '',A)') ELEM_NAME(IBEG:IEND)

          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            SCHEME=NQS(NEELEM(1,nr))
            nb=NQSCNB(SCHEME)

            IF(NQXI(0,SCHEME).EQ.1) THEN

              WRITE(BASES,'(A)') 'l.Lagrange'
              CALL STRING_TRIM(BASES,IBEG,IEND)
              WRITE(IFILE,'( '' Shape.  Dimension=1'' )')
              WRITE(IFILE,'( '' #Scale factor sets= 1'')')
              WRITE(IFILE,'(3X,A,'', #Scale factors='',I2)')
     &          BASES(IBEG:IEND),NNT(nb)

              WRITE(IFILE,'( '' #Nodes='',I2)') NNT(nb)
              NUMFIELDS=1
              IF(ACTIVATION_TIME) NUMFIELDS=NUMFIELDS+1
              IF(POTENTIAL) NUMFIELDS=NUMFIELDS+1
              IF(PRESSURE) NUMFIELDS=NUMFIELDS+1
              IF(RADIUS) NUMFIELDS=NUMFIELDS+1
              IF(VELOCITY) NUMFIELDS=NUMFIELDS+1
              IF(PATH) NUMFIELDS=NUMFIELDS+1
C PM 26-JUL-01
              IF(ELASTIC_TUBE) NUMFIELDS=4
              CALL ASSERT(NUMFIELDS.LT.10,' >>Too many fields',
     &          ERROR,*9999)  !fields written into I1
              WRITE(IFILE,'(1X,''#Fields='',I1)') NUMFIELDS

              WRITE(IFILE,'(1X,''1) coordinates, coordinate, '
     &          //'rectangular cartesian, #Components='',I1)')
     &          NJ_LOC(NJL_GEOM,0,nr)

              DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                IF(nj.EQ.1) THEN
                  WRITE(CHAR1,'(A)') 'x'
                  CALL STRING_TRIM(CHAR1,IBEG2,IEND2)
                ELSE IF(nj.EQ.2) THEN
                  WRITE(CHAR1,'(A)') 'y'
                  CALL STRING_TRIM(CHAR1,IBEG2,IEND2)
                ELSE IF(nj.EQ.3) THEN
                  WRITE(CHAR1,'(A)') 'z'
                  CALL STRING_TRIM(CHAR1,IBEG2,IEND2)
                ELSE
                  CALL ASSERT(.FALSE.,' Incorrect number on njs',
     &              ERROR,*9999)
                ENDIF

                WRITE(IFILE,'(3X,A,''.  '',A,'', no modify, '
     &            //'standard node based.'')') CHAR1(IBEG2:IEND2),
     &            BASES(IBEG:IEND)

                WRITE(IFILE,'(5X,''#Nodes='',I2)') NNT(nb)

                DO nn=1,NNT(nb)
                  WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                  WRITE(IFILE,'(7X,''Value indices:     1'')')
                  WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                ENDDO
              ENDDO

              NUMFIELDS=2
              IF(ACTIVATION_TIME) THEN
                WRITE(IFILE,'(1X,I1,'') activation_time, field, '
     &            //'rectangular cartesian, #Components=1'')')
     &            NUMFIELDS
                WRITE(IFILE,'(3X,''activation_time.  '',A,'', no '
     &            //'modify, standard node based.'')') BASES(IBEG:IEND)
                WRITE(IFILE,'(5X,''#Nodes='',I2)') NNT(nb)
                DO nn=1,NNT(nb)
                  WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                  WRITE(IFILE,'(7X,''Value indices:     1'')')
                  WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                ENDDO
                NUMFIELDS=NUMFIELDS+1
              ENDIF
              IF(POTENTIAL) THEN
                WRITE(IFILE,'(1X,I1,'') potential, field, '
     &            //'rectangular cartesian, #Components=1'')')
     &            NUMFIELDS
                WRITE(IFILE,'(3X,''potential.  '',A,'', no '
     &            //'modify, standard node based.'')') BASES(IBEG:IEND)
                WRITE(IFILE,'(5X,''#Nodes='',I2)') NNT(nb)
                DO nn=1,NNT(nb)
                  WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                  WRITE(IFILE,'(7X,''Value indices:     1'')')
                  WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                ENDDO
                NUMFIELDS=NUMFIELDS+1
              ENDIF
              IF(PRESSURE) THEN
                WRITE(IFILE,'(1X,I1,'') pressure, field, '
     &            //'rectangular cartesian, #Components=1'')')
     &            NUMFIELDS
                WRITE(IFILE,'(3X,''pressure.  '',A,'', no '
     &            //'modify, standard node based.'')') BASES(IBEG:IEND)
                WRITE(IFILE,'(5X,''#Nodes='',I2)') NNT(nb)
                DO nn=1,NNT(nb)
                  WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                  WRITE(IFILE,'(7X,''Value indices:     1'')')
                  WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                ENDDO
                NUMFIELDS=NUMFIELDS+1
              ENDIF
              IF(RADIUS) THEN
                WRITE(IFILE,'(1X,I1,'') radius, field, '
     &            //'rectangular cartesian, #Components=1'')')
     &            NUMFIELDS
                WRITE(IFILE,'(3X,''radius.  '',A,'', no '
     &            //'modify, standard node based.'')') BASES(IBEG:IEND)
                WRITE(IFILE,'(5X,''#Nodes='',I2)') NNT(nb)
                DO nn=1,NNT(nb)
                  WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                  WRITE(IFILE,'(7X,''Value indices:     1'')')
                  WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                ENDDO
                NUMFIELDS=NUMFIELDS+1
              ENDIF
              IF(VELOCITY) THEN
                WRITE(IFILE,'(1X,I1,'') velocity, field, '
     &            //'rectangular cartesian, #Components=1'')')
     &            NUMFIELDS
                WRITE(IFILE,'(3X,''velocity.  '',A,'', no '
     &            //'modify, standard node based.'')') BASES(IBEG:IEND)
                WRITE(IFILE,'(5X,''#Nodes='',I2)') NNT(nb)
                DO nn=1,NNT(nb)
                  WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                  WRITE(IFILE,'(7X,''Value indices:     1'')')
                  WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                ENDDO
                NUMFIELDS=NUMFIELDS+1
              ENDIF
              IF(PATH) THEN
                WRITE(IFILE,'(1X,I1,'') path, field, '
     &            //'rectangular cartesian, #Components=1'')')
     &            NUMFIELDS
                WRITE(IFILE,'(3X,''path.  '',A,'', no '
     &            //'modify, standard node based.'')') BASES(IBEG:IEND)
                WRITE(IFILE,'(5X,''#Nodes='',I2)') NNT(nb)
                DO nn=1,NNT(nb)
                  WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                  WRITE(IFILE,'(7X,''Value indices:     1'')')
                  WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                ENDDO
                NUMFIELDS=NUMFIELDS+1
              ENDIF

C PM 26-JUL-01 : flow in elastic tube
              IF(ELASTIC_TUBE) THEN
                WRITE(IFILE,'(1X,I1,'') pressure, field, '
     &            //'rectangular cartesian, #Components=1'')')
     &            NUMFIELDS
                WRITE(IFILE,'(3X,''pressure.  '',A,'', no '
     &            //'modify, standard node based.'')') BASES(IBEG:IEND)
                WRITE(IFILE,'(5X,''#Nodes='',I2)') NNT(nb)
                DO nn=1,NNT(nb)
                  WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                  WRITE(IFILE,'(7X,''Value indices:     1'')')
                  WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                ENDDO
                WRITE(IFILE,'(1X,I1,'') radius, field, '
     &            //'rectangular cartesian, #Components=1'')')
     &            NUMFIELDS+1
                WRITE(IFILE,'(3X,''radius.  '',A,'', no '
     &            //'modify, standard node based.'')') BASES(IBEG:IEND)
                WRITE(IFILE,'(5X,''#Nodes='',I2)') NNT(nb)
                DO nn=1,NNT(nb)
                  WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                  WRITE(IFILE,'(7X,''Value indices:     1'')')
                  WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                ENDDO
                WRITE(IFILE,'(1X,I1,'') velocity, field, '
     &            //'rectangular cartesian, #Components=1'')')
     &            NUMFIELDS+2
                WRITE(IFILE,'(3X,''velocity.  '',A,'', no '
     &            //'modify, standard node based.'')') BASES(IBEG:IEND)
                WRITE(IFILE,'(5X,''#Nodes='',I2)') NNT(nb)
                DO nn=1,NNT(nb)
                  WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                  WRITE(IFILE,'(7X,''Value indices:     1'')')
                  WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                ENDDO
C                NUMFIELDS=NUMFIELDS+1
              ENDIF

              ELEMENT=1
              SCALEFACTOR=1.0d0
              DO nq=NQR(1,nr),NQR(2,nr)
                IF(NXQ(1,0,nq,1).GT.0) THEN
                  DO nqq=1,NXQ(1,0,nq,1)
                    mq=NXQ(1,nqq,nq,1)
                    WRITE(IFILE,'(1X,''Element: '',I6,'' 0 0'' )')
     &                ELEMENT+offset_elem
                    ELEMENT=ELEMENT+1
                    WRITE(IFILE,'(3X,''Nodes:'' )')
                    WRITE(IFILE,'(4X,2(1X,I6))') nq+offset_node,
     &                mq+offset_node
                    WRITE(IFILE,'(3X,''Scale factors:'' )')
                    WRITE(IFILE,'(4X,5(1X,E24.16))') SCALEFACTOR,
     &                SCALEFACTOR
                  ENDDO
                ENDIF
              ENDDO

            ELSE IF(NQXI(0,SCHEME).EQ.2) THEN

              NQLIST(0)=0
              DO nq=NQR(1,nr),NQR(2,nr)
                IF((NXQ(1,0,nq,1).GT.0).AND.(NXQ(2,0,nq,1).GT.0)) THEN
                  IF(NXQ(1,0,NXQ(2,1,nq,1),1).GT.0) THEN
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq
                  ENDIF
                ENDIF
              ENDDO
              CALL ASSERT(NQLIST(0).GT.0,' No grid elements created',
     &          ERROR,*9999)

              IF( .NOT.LIBMESH ) THEN

                WRITE(IFILE,'( '' Shape.  Dimension=1'' )')
                DO nl=1,NLE(nb)*NQLIST(0)
                  WRITE(IFILE,'( '' Element: 0 0 '',I6)')
     &                 nl+offset_line
                ENDDO

              ENDIF

              WRITE(BASES,'(A)') 'l.Lagrange*l.Lagrange'
              CALL STRING_TRIM(BASES,IBEG,IEND)
              WRITE(IFILE,'( '' Shape.  Dimension=2'' )')
              WRITE(IFILE,'( '' #Scale factor sets= 1'')')
              WRITE(IFILE,'(3X,A,'', #Scale factors='',I2)')
     &          BASES(IBEG:IEND),NNT(nb)
              WRITE(IFILE,'( '' #Nodes='',I2)') NNT(nb)
              NUMFIELDS=1
              IF(ACTIVATION_TIME) NUMFIELDS=NUMFIELDS+1
              IF(POTENTIAL) NUMFIELDS=NUMFIELDS+1
              IF(PRESSURE) NUMFIELDS=NUMFIELDS+1
              IF(RADIUS) NUMFIELDS=NUMFIELDS+1
              IF(VELOCITY) NUMFIELDS=NUMFIELDS+1
              IF(PATH) NUMFIELDS=NUMFIELDS+1

              CALL ASSERT(NUMFIELDS.LT.10,' >>Too many fields',
     &          ERROR,*9999)  !fields written into I1
              WRITE(IFILE,'(1X,''#Fields='',I1)') NUMFIELDS
              WRITE(IFILE,'(1X,''1) coordinates, coordinate, '
     &          //'rectangular cartesian, #Components='',I1)')
     &          NJ_LOC(NJL_GEOM,0,nr)

              DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                IF(nj.EQ.1) THEN
                  WRITE(CHAR1,'(A)') 'x'
                  CALL STRING_TRIM(CHAR1,IBEG2,IEND2)
                ELSEIF(nj.EQ.2) THEN
                  WRITE(CHAR1,'(A)') 'y'
                  CALL STRING_TRIM(CHAR1,IBEG2,IEND2)
                ELSEIF(nj.EQ.3) THEN
                  WRITE(CHAR1,'(A)') 'z'
                  CALL STRING_TRIM(CHAR1,IBEG2,IEND2)
                ELSE
                  CALL ASSERT(.FALSE.,' Incorrect number on njs',
     &              ERROR,*9999)
                ENDIF

                WRITE(IFILE,'(3X,A,''.  '',A,'', no modify, '
     &            //'standard node based.'')') CHAR1(IBEG2:IEND2),
     &            BASES(IBEG:IEND)

                WRITE(IFILE,'(5X,''#Nodes='',I2)') NNT(nb)
                DO nn=1,NNT(nb)
                  WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                  WRITE(IFILE,'(7X,''Value indices:     1'')')
                  WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                ENDDO
              ENDDO

              NUMFIELDS=2
              IF(ACTIVATION_TIME) THEN
                WRITE(IFILE,'(1X,I1,'') activation_time, field, '
     &            //'rectangular cartesian, #Components=1'')')
     &            NUMFIELDS
                WRITE(IFILE,'(3X,''activation_time.  '',A,'', no '
     &            //'modify, standard node based.'')') BASES(IBEG:IEND)
                WRITE(IFILE,'(5X,''#Nodes='',I2)') NNT(nb)
                DO nn=1,NNT(nb)
                  WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                  WRITE(IFILE,'(7X,''Value indices:     1'')')
                  WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                ENDDO
                NUMFIELDS=NUMFIELDS+1
              ENDIF
              IF(POTENTIAL) THEN
                WRITE(IFILE,'(1X,I1,'') potential, field, '
     &            //'rectangular cartesian, #Components=1'')')
     &            NUMFIELDS
                WRITE(IFILE,'(3X,''potential.  '',A,'', no '
     &            //'modify, standard node based.'')') BASES(IBEG:IEND)
                WRITE(IFILE,'(5X,''#Nodes='',I2)') NNT(nb)
                DO nn=1,NNT(nb)
                  WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                  WRITE(IFILE,'(7X,''Value indices:     1'')')
                  WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                ENDDO
                NUMFIELDS=NUMFIELDS+1
              ENDIF
              IF(PRESSURE) THEN
                WRITE(IFILE,'(1X,I1,'') pressure, field, '
     &            //'rectangular cartesian, #Components=1'')')
     &            NUMFIELDS
                WRITE(IFILE,'(3X,''pressure.  '',A,'', no '
     &            //'modify, standard node based.'')') BASES(IBEG:IEND)
                WRITE(IFILE,'(5X,''#Nodes='',I2)') NNT(nb)
                DO nn=1,NNT(nb)
                  WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                  WRITE(IFILE,'(7X,''Value indices:     1'')')
                  WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                ENDDO
                NUMFIELDS=NUMFIELDS+1
              ENDIF
              IF(RADIUS) THEN
                WRITE(IFILE,'(1X,I1,'') radius, field, '
     &            //'rectangular cartesian, #Components=1'')')
     &            NUMFIELDS
                WRITE(IFILE,'(3X,''radius.  '',A,'', no '
     &            //'modify, standard node based.'')') BASES(IBEG:IEND)
                WRITE(IFILE,'(5X,''#Nodes='',I2)') NNT(nb)
                DO nn=1,NNT(nb)
                  WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                  WRITE(IFILE,'(7X,''Value indices:     1'')')
                  WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                ENDDO
                NUMFIELDS=NUMFIELDS+1
              ENDIF
              IF(VELOCITY) THEN
                WRITE(IFILE,'(1X,I1,'') velocity, field, '
     &            //'rectangular cartesian, #Components=1'')')
     &            NUMFIELDS
                WRITE(IFILE,'(3X,''velocity.  '',A,'', no '
     &            //'modify, standard node based.'')') BASES(IBEG:IEND)
                WRITE(IFILE,'(5X,''#Nodes='',I2)') NNT(nb)
                DO nn=1,NNT(nb)
                  WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                  WRITE(IFILE,'(7X,''Value indices:     1'')')
                  WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                ENDDO
                NUMFIELDS=NUMFIELDS+1
              ENDIF
              IF(PATH) THEN
                WRITE(IFILE,'(1X,I1,'') path, field, '
     &            //'rectangular cartesian, #Components=1'')')
     &            NUMFIELDS
                WRITE(IFILE,'(3X,''path.  '',A,'', no '
     &            //'modify, standard node based.'')') BASES(IBEG:IEND)
                WRITE(IFILE,'(5X,''#Nodes='',I2)') NNT(nb)
                DO nn=1,NNT(nb)
                  WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                  WRITE(IFILE,'(7X,''Value indices:     1'')')
                  WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                ENDDO
                NUMFIELDS=NUMFIELDS+1
              ENDIF

              ELEMENT=1
              SCALEFACTOR=1.0d0
              DO nq=1,NQLIST(0)
                mq=NQLIST(nq)
                mq2=NXQ(1,1,mq,1)
                mq3=NXQ(2,1,mq,1)
                mq4=NXQ(1,1,NXQ(2,1,mq,1),1)
                WRITE(IFILE,'(1X,''Element: '',I6,'' 0 0'' )')
     &            ELEMENT+offset_elem
                WRITE(IFILE,'(3X,''Faces:'' )')
                DO nl=1,NLE(nb)
                  WRITE(IFILE,'(3X,''0 0 '',I6)')
     &              nl+(NLE(nb)*(ELEMENT-1))+offset_line
                ENDDO
                WRITE(IFILE,'(3X,''Nodes:'' )')
                WRITE(IFILE,'(4X,4(1X,I6))') mq+offset_node,
     &            mq2+offset_node,mq3+offset_node,mq4+offset_node
                WRITE(IFILE,'(3X,''Scale factors:'' )')
                WRITE(IFILE,'(4X,5(1X,E24.16))') SCALEFACTOR,
     &            SCALEFACTOR,SCALEFACTOR,SCALEFACTOR
                ELEMENT=ELEMENT+1
              ENDDO

            ELSEIF(NQXI(0,SCHEME).EQ.3) THEN

              IF(SURFACE) THEN
                NQLIST(0)=0
                DO nq=NQR(1,nr),NQR(2,nr)
                  IF(NWQ(1,nq,1).NE.0) THEN
                    IF((NXQ(1,0,nq,1).GT.0).AND.(NXQ(2,0,nq,1).GT.0))
     &                THEN
                      IF(NXQ(1,0,NXQ(2,1,nq,1),1).GT.0) THEN
                        IF((NWQ(1,NXQ(1,1,nq,1),1).NE.0).AND.
     &                    (NWQ(1,NXQ(2,1,nq,1),1).NE.0).AND.
     &                    (NWQ(1,NXQ(1,1,NXQ(2,1,nq,1),1),1).NE.0)) THEN
                          NQLIST(0)=NQLIST(0)+1
                          NQLIST(NQLIST(0))=nq
                        ENDIF
                      ENDIF
                    ENDIF
                    IF((NXQ(1,0,nq,1).GT.0).AND.(NXQ(3,0,nq,1).GT.0))
     &                THEN
                      IF(NXQ(1,0,NXQ(3,1,nq,1),1).GT.0) THEN
                        IF((NWQ(1,NXQ(1,1,nq,1),1).NE.0).AND.
     &                    (NWQ(1,NXQ(3,1,nq,1),1).NE.0).AND.
     &                    (NWQ(1,NXQ(1,1,NXQ(3,1,nq,1),1),1).NE.0)) THEN
                          NQLIST(0)=NQLIST(0)+1
                          NQLIST(NQLIST(0))=nq
                        ENDIF
                      ENDIF
                    ENDIF
                    IF((NXQ(2,0,nq,1).GT.0).AND.(NXQ(3,0,nq,1).GT.0))
     &                THEN
                      IF(NXQ(2,0,NXQ(3,1,nq,1),1).GT.0) THEN
                        IF((NWQ(1,NXQ(2,1,nq,1),1).NE.0).AND.
     &                    (NWQ(1,NXQ(3,1,nq,1),1).NE.0).AND.
     &                    (NWQ(1,NXQ(2,1,NXQ(3,1,nq,1),1),1).NE.0)) THEN
                          NQLIST(0)=NQLIST(0)+1
                          NQLIST(NQLIST(0))=nq
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDIF
                ENDDO
                CALL ASSERT(NQLIST(0).GT.0,' No grid elements created',
     &            ERROR,*9999)

                WRITE(IFILE,'( '' Shape.  Dimension=1'' )')
                DO nl=1,4*NQLIST(0)
                  WRITE(IFILE,'( '' Element: 0 0 '',I6)')
     &              nl+offset_line
                ENDDO

                WRITE(BASES,'(A)') 'l.Lagrange*l.Lagrange'
                CALL STRING_TRIM(BASES,IBEG,IEND)
                WRITE(IFILE,'( '' Shape.  Dimension=2'' )')
                WRITE(IFILE,'( '' #Scale factor sets= 1'')')
                WRITE(IFILE,'(3X,A,'', #Scale factors= 4'')')
     &            BASES(IBEG:IEND)
                WRITE(IFILE,'( '' #Nodes= 4'')')
                NUMFIELDS=1
                IF(ACTIVATION_TIME) NUMFIELDS=NUMFIELDS+1
                IF(POTENTIAL) NUMFIELDS=NUMFIELDS+1
                IF(PRESSURE) NUMFIELDS=NUMFIELDS+1
                IF(RADIUS) NUMFIELDS=NUMFIELDS+1
                IF(VELOCITY) NUMFIELDS=NUMFIELDS+1
                IF(PATH) NUMFIELDS=NUMFIELDS+1

                CALL ASSERT(NUMFIELDS.LT.10,' >>Too many fields',
     &            ERROR,*9999)  !fields written into I1
                WRITE(IFILE,'(1X,''#Fields='',I1)') NUMFIELDS
                WRITE(IFILE,'(1X,''1) coordinates, coordinate, '
     &            //'rectangular cartesian, #Components='',I1)')
     &            NJ_LOC(NJL_GEOM,0,nr)

                DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                  IF(nj.EQ.1) THEN
                    WRITE(CHAR1,'(A)') 'x'
                    CALL STRING_TRIM(CHAR1,IBEG2,IEND2)
                  ELSEIF(nj.EQ.2) THEN
                    WRITE(CHAR1,'(A)') 'y'
                    CALL STRING_TRIM(CHAR1,IBEG2,IEND2)
                  ELSEIF(nj.EQ.3) THEN
                    WRITE(CHAR1,'(A)') 'z'
                    CALL STRING_TRIM(CHAR1,IBEG2,IEND2)
                  ELSE
                    CALL ASSERT(.FALSE.,' Incorrect number on njs',
     &                ERROR,*9999)
                  ENDIF

                  WRITE(IFILE,'(3X,A,''.  '',A,'', no modify, '
     &              //'standard node based.'')') CHAR1(IBEG2:IEND2),
     &              BASES(IBEG:IEND)

                  WRITE(IFILE,'(5X,''#Nodes= 4'')')
                  DO nn=1,4
                    WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                    WRITE(IFILE,'(7X,''Value indices:     1'')')
                    WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                  ENDDO
                ENDDO

                NUMFIELDS=2
                IF(ACTIVATION_TIME) THEN
                  WRITE(IFILE,'(1X,I1,'') activation_time, field, '
     &              //'rectangular cartesian, #Components=1'')')
     &              NUMFIELDS
                  WRITE(IFILE,'(3X,''activation_time.  '',A,'', no '
     &              //'modify, standard node based.'')')
     &            BASES(IBEG:IEND)
                  WRITE(IFILE,'(5X,''#Nodes= 4'')')
                  DO nn=1,4
                    WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                    WRITE(IFILE,'(7X,''Value indices:     1'')')
                    WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                  ENDDO
                  NUMFIELDS=NUMFIELDS+1
                ENDIF
                IF(POTENTIAL) THEN
                  WRITE(IFILE,'(1X,I1,'') potential, field, '
     &              //'rectangular cartesian, #Components=1'')')
     &              NUMFIELDS
                  WRITE(IFILE,'(3X,''potential.  '',A,'', no '
     &              //'modify, standard node based.'')')
     &            BASES(IBEG:IEND)
                  WRITE(IFILE,'(5X,''#Nodes= 4'')')
                  DO nn=1,4
                    WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                    WRITE(IFILE,'(7X,''Value indices:     1'')')
                    WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                  ENDDO
                  NUMFIELDS=NUMFIELDS+1
                ENDIF
                IF(PRESSURE) THEN
                  WRITE(IFILE,'(1X,I1,'') pressure, field, '
     &              //'rectangular cartesian, #Components=1'')')
     &              NUMFIELDS
                  WRITE(IFILE,'(3X,''pressure.  '',A,'', no '
     &              //'modify, standard node based.'')')
     &            BASES(IBEG:IEND)
                  WRITE(IFILE,'(5X,''#Nodes= 4'')')
                  DO nn=1,4
                    WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                    WRITE(IFILE,'(7X,''Value indices:     1'')')
                    WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                  ENDDO
                  NUMFIELDS=NUMFIELDS+1
                ENDIF
                IF(RADIUS) THEN
                  WRITE(IFILE,'(1X,I1,'') radius, field, '
     &              //'rectangular cartesian, #Components=1'')')
     &              NUMFIELDS
                  WRITE(IFILE,'(3X,''radius.  '',A,'', no '
     &              //'modify, standard node based.'')')
     &            BASES(IBEG:IEND)
                  WRITE(IFILE,'(5X,''#Nodes= 4'')')
                  DO nn=1,4
                    WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                    WRITE(IFILE,'(7X,''Value indices:     1'')')
                    WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                  ENDDO
                  NUMFIELDS=NUMFIELDS+1
                ENDIF
                IF(VELOCITY) THEN
                  WRITE(IFILE,'(1X,I1,'') velocity, field, '
     &              //'rectangular cartesian, #Components=1'')')
     &              NUMFIELDS
                  WRITE(IFILE,'(3X,''velocity.  '',A,'', no '
     &              //'modify, standard node based.'')')
     &            BASES(IBEG:IEND)
                  WRITE(IFILE,'(5X,''#Nodes= 4'')')
                  DO nn=1,4
                    WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                    WRITE(IFILE,'(7X,''Value indices:     1'')')
                    WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                  ENDDO
                  NUMFIELDS=NUMFIELDS+1
                ENDIF
                IF(PATH) THEN
                  WRITE(IFILE,'(1X,I1,'') path, field, '
     &              //'rectangular cartesian, #Components=1'')')
     &              NUMFIELDS
                  WRITE(IFILE,'(3X,''path.  '',A,'', no '
     &              //'modify, standard node based.'')')
     &            BASES(IBEG:IEND)
                  WRITE(IFILE,'(5X,''#Nodes= 4'')')
                  DO nn=1,4
                    WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                    WRITE(IFILE,'(7X,''Value indices:     1'')')
                    WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                  ENDDO
                  NUMFIELDS=NUMFIELDS+1
                ENDIF

                ELEMENT=1
                SCALEFACTOR=1.0d0
                DO nq=1,NQLIST(0)
                  XI1XI2=.FALSE.
                  XI1XI3=.FALSE.
                  XI2XI3=.FALSE.
                  IF(NQLIST(nq-1).NE.NQLIST(nq)) THEN
                    FACE1=.FALSE.
                    FACE2=.FALSE.
                    FACE3=.FALSE.
                  ENDIF
                  mq=NQLIST(nq)

                  IF(NWQ(1,mq,1).NE.0) THEN
                    IF((NXQ(1,0,mq,1).GT.0).AND.(NXQ(2,0,mq,1).GT.0))
     &                THEN
                      IF(NXQ(1,0,NXQ(2,1,mq,1),1).GT.0) THEN
                        IF((NWQ(1,NXQ(1,1,mq,1),1).NE.0).AND.
     &                    (NWQ(1,NXQ(2,1,mq,1),1).NE.0).AND.
     &                    (NWQ(1,NXQ(1,1,NXQ(2,1,mq,1),1),1).NE.0)) THEN
                          IF(.NOT.FACE1) THEN
                            XI1XI2=.TRUE.
                            FACE1=.TRUE.
                            GOTO 100
                          ENDIF
                        ENDIF
                      ENDIF
                    ENDIF
                    IF((NXQ(1,0,mq,1).GT.0).AND.(NXQ(3,0,mq,1).GT.0))
     &                THEN
                      IF(NXQ(1,0,NXQ(3,1,mq,1),1).GT.0) THEN
                        IF((NWQ(1,NXQ(1,1,mq,1),1).NE.0).AND.
     &                    (NWQ(1,NXQ(3,1,mq,1),1).NE.0).AND.
     &                    (NWQ(1,NXQ(1,1,NXQ(3,1,mq,1),1),1).NE.0)) THEN
                          IF(.NOT.FACE2) THEN
                            XI1XI3=.TRUE.
                            FACE2=.TRUE.
                            GOTO 100
                          ENDIF
                        ENDIF
                      ENDIF
                    ENDIF
                    IF((NXQ(2,0,mq,1).GT.0).AND.(NXQ(3,0,mq,1).GT.0))
     &                THEN
                      IF(NXQ(2,0,NXQ(3,1,mq,1),1).GT.0) THEN
                        IF((NWQ(1,NXQ(2,1,mq,1),1).NE.0).AND.
     &                    (NWQ(1,NXQ(3,1,mq,1),1).NE.0).AND.
     &                    (NWQ(1,NXQ(2,1,NXQ(3,1,mq,1),1),1).NE.0)) THEN
                          IF(.NOT.FACE3) THEN
                            XI2XI3=.TRUE.
                            FACE3=.TRUE.
                            GOTO 100
                          ENDIF
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDIF

 100              IF(XI1XI2) THEN
                    mq2=NXQ(1,1,mq,1)
                    mq3=NXQ(2,1,mq,1)
                    mq4=NXQ(1,1,NXQ(2,1,mq,1),1)
                  ELSEIF(XI1XI3) THEN
                    mq2=NXQ(1,1,mq,1)
                    mq3=NXQ(3,1,mq,1)
                    mq4=NXQ(1,1,NXQ(3,1,mq,1),1)
                  ELSEIF(XI2XI3) THEN
                    mq2=NXQ(2,1,mq,1)
                    mq3=NXQ(3,1,mq,1)
                    mq4=NXQ(2,1,NXQ(3,1,mq,1),1)
                  ENDIF

                  WRITE(IFILE,'(1X,''Element: '',I6,'' 0 0'' )')
     &              ELEMENT+offset_elem
                  WRITE(IFILE,'(3X,''Faces:'' )')
                  DO nl=1,4
                    WRITE(IFILE,'(3X,''0 0 '',I6)')
     &                nl+(4*(ELEMENT-1))+offset_line
                  ENDDO
                  WRITE(IFILE,'(3X,''Nodes:'' )')
                  WRITE(IFILE,'(4X,4(1X,I6))') mq+offset_node,
     &              mq2+offset_node,mq3+offset_node,mq4+offset_node
                  WRITE(IFILE,'(3X,''Scale factors:'' )')
                  WRITE(IFILE,'(4X,5(1X,E24.16))') SCALEFACTOR,
     &              SCALEFACTOR,SCALEFACTOR,SCALEFACTOR
                  ELEMENT=ELEMENT+1
                ENDDO

              ELSE

                NQLIST(0)=0
                IF(USE_LAT.EQ.0) THEN
                  DO nq=NQR(1,nr),NQR(2,nr)
                    IF((NXQ(1,0,nq,1).GT.0).AND.(NXQ(2,0,nq,1).GT.0)
     &                .AND.(NXQ(3,0,nq,1).GT.0)) THEN
                      IF((NXQ(1,0,NXQ(2,1,nq,1),1).GT.0)
     &                  .AND.(NXQ(1,0,NXQ(3,1,nq,1),1).GT.0)
     &                  .AND.(NXQ(2,0,NXQ(3,1,nq,1),1).GT.0)) THEN
                        IF(NXQ(2,0,NXQ(1,1,NXQ(3,1,nq,1),1),1).GT.0)THEN
                          NQLIST(0)=NQLIST(0)+1
                          NQLIST(NQLIST(0))=nq
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDDO
                ELSE !Lattice-based grid
                  DO ne=1,NEQM
C***                Loop over subelements (if any)
                    DO nk=0,MAX(NQXI(3,SCHEME)-2,0)
                      DO nj=0,MAX(NQXI(2,SCHEME)-2,0)
                        DO ni=0,MAX(NQXI(1,SCHEME)-2,0)
                          NQLIST(0)=NQLIST(0)+1
                          NQLIST(NQLIST(0))=NLATNE(ne)+ni+
     &                      nj*NQXI(1,SCHEME)+
     &                      nk*NQXI(1,SCHEME)*NQXI(2,SCHEME)
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                ENDIF
                CALL ASSERT(NQLIST(0).GT.0,' No grid elements created',
     &            ERROR,*9999)

                IF( .NOT.LIBMESH ) THEN

                  WRITE(IFILE,'( '' Shape.  Dimension=1'' )')
                  DO nl=1,NLE(nb)*NQLIST(0)
                    WRITE(IFILE,'( '' Element: 0 0 '',I6)')
     &                   nl+offset_line
                  ENDDO

                  WRITE(IFILE,'( '' Shape.  Dimension=2'' )')
                  DO nq=1,NQLIST(0)
                    DO nf=1,NFE(nb)
                      WRITE(IFILE,'( '' Element: 0 '',I6,'' 0'' )')
     &                     nf+((nq-1)*NFE(nb))+offset_face
                      WRITE(IFILE,'( ''   Faces:'' )')
                      WRITE(IFILE,'( ''   0 0 '',I6)')(FACE(1+4*(nf-1))+
     &                     (NLE(nb)*(nq-1)))+offset_line
                      WRITE(IFILE,'( ''   0 0 '',I6)')(FACE(2+4*(nf-1))+
     &                     (NLE(nb)*(nq-1)))+offset_line
                      WRITE(IFILE,'( ''   0 0 '',I6)')(FACE(3+4*(nf-1))+
     &                     (NLE(nb)*(nq-1)))+offset_line
                      WRITE(IFILE,'( ''   0 0 '',I6)')(FACE(4+4*(nf-1))+
     &                     (NLE(nb)*(nq-1)))+offset_line
                    ENDDO
                  ENDDO

                ENDIF

                WRITE(BASES,'(A)') 'l.Lagrange*l.Lagrange*l.Lagrange'
                CALL STRING_TRIM(BASES,IBEG,IEND)
                WRITE(IFILE,'( '' Shape.  Dimension=3'' )')
                WRITE(IFILE,'( '' #Scale factor sets= 1'')')
                WRITE(IFILE,'(3X,A,'', #Scale factors='',I2)')
     &            BASES(IBEG:IEND),NNT(nb)
                WRITE(IFILE,'( '' #Nodes='',I2)') NNT(nb)
                NUMFIELDS=1
                IF(ACTIVATION_TIME) NUMFIELDS=NUMFIELDS+1
                IF(POTENTIAL) NUMFIELDS=NUMFIELDS+1
                IF(PRESSURE) NUMFIELDS=NUMFIELDS+1
                IF(RADIUS) NUMFIELDS=NUMFIELDS+1
                IF(VELOCITY) NUMFIELDS=NUMFIELDS+1
                IF(PATH) NUMFIELDS=NUMFIELDS+1

                CALL ASSERT(NUMFIELDS.LT.10,' >>Too many fields',
     &            ERROR,*9999)  !fields written into I1
                WRITE(IFILE,'(1X,''#Fields='',I1)') NUMFIELDS
                WRITE(IFILE,'(1X,''1) coordinates, coordinate, '
     &            //'rectangular cartesian, #Components='',I1)')
     &            NJ_LOC(NJL_GEOM,0,nr)

                DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                  IF(nj.EQ.1) THEN
                    WRITE(CHAR1,'(A)') 'x'
                    CALL STRING_TRIM(CHAR1,IBEG2,IEND2)
                  ELSEIF(nj.EQ.2) THEN
                    WRITE(CHAR1,'(A)') 'y'
                    CALL STRING_TRIM(CHAR1,IBEG2,IEND2)
                  ELSEIF(nj.EQ.3) THEN
                    WRITE(CHAR1,'(A)') 'z'
                    CALL STRING_TRIM(CHAR1,IBEG2,IEND2)
                  ELSE
                    CALL ASSERT(.FALSE.,' Incorrect number on njs',
     &                ERROR,*9999)
                  ENDIF

                  WRITE(IFILE,'(3X,A,''.  '',A,'', no modify, '
     &              //'standard node based.'')') CHAR1(IBEG2:IEND2),
     &              BASES(IBEG:IEND)

                  WRITE(IFILE,'(5X,''#Nodes='',I2)') NNT(nb)
                  DO nn=1,NNT(nb)
                    WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                    WRITE(IFILE,'(7X,''Value indices:     1'')')
                    WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                  ENDDO
                ENDDO

                NUMFIELDS=2
                IF(ACTIVATION_TIME) THEN
                  WRITE(IFILE,'(1X,I1,'') activation_time, field, '
     &              //'rectangular cartesian, #Components=1'')')
     &              NUMFIELDS
                  WRITE(IFILE,'(3X,''activation_time.  '',A,'', no '
     &              //'modify, standard node based.'')')
     &            BASES(IBEG:IEND)
                  WRITE(IFILE,'(5X,''#Nodes='',I2)') NNT(nb)
                  DO nn=1,NNT(nb)
                    WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                    WRITE(IFILE,'(7X,''Value indices:     1'')')
                    WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                  ENDDO
                  NUMFIELDS=NUMFIELDS+1
                ENDIF
                IF(POTENTIAL) THEN
                  WRITE(IFILE,'(1X,I1,'') potential, field, '
     &              //'rectangular cartesian, #Components=1'')')
     &              NUMFIELDS
                  WRITE(IFILE,'(3X,''potential.  '',A,'', no '
     &              //'modify, standard node based.'')')
     &            BASES(IBEG:IEND)
                  WRITE(IFILE,'(5X,''#Nodes='',I2)') NNT(nb)
                  DO nn=1,NNT(nb)
                    WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                    WRITE(IFILE,'(7X,''Value indices:     1'')')
                    WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                  ENDDO
                  NUMFIELDS=NUMFIELDS+1
                ENDIF
                IF(PRESSURE) THEN
                  WRITE(IFILE,'(1X,I1,'') pressure, field, '
     &              //'rectangular cartesian, #Components=1'')')
     &              NUMFIELDS
                  WRITE(IFILE,'(3X,''pressure.  '',A,'', no '
     &              //'modify, standard node based.'')')
     &            BASES(IBEG:IEND)
                  WRITE(IFILE,'(5X,''#Nodes='',I2)') NNT(nb)
                  DO nn=1,NNT(nb)
                    WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                    WRITE(IFILE,'(7X,''Value indices:     1'')')
                    WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                  ENDDO
                  NUMFIELDS=NUMFIELDS+1
                ENDIF
                IF(RADIUS) THEN
                  WRITE(IFILE,'(1X,I1,'') radius, field, '
     &              //'rectangular cartesian, #Components=1'')')
     &              NUMFIELDS
                  WRITE(IFILE,'(3X,''radius.  '',A,'', no '
     &              //'modify, standard node based.'')')
     &            BASES(IBEG:IEND)
                  WRITE(IFILE,'(5X,''#Nodes='',I2)') NNT(nb)
                  DO nn=1,NNT(nb)
                    WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                    WRITE(IFILE,'(7X,''Value indices:     1'')')
                    WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                  ENDDO
                  NUMFIELDS=NUMFIELDS+1
                ENDIF
                IF(VELOCITY) THEN
                  WRITE(IFILE,'(1X,I1,'') velocity, field, '
     &              //'rectangular cartesian, #Components=1'')')
     &              NUMFIELDS
                  WRITE(IFILE,'(3X,''velocity.  '',A,'', no '
     &              //'modify, standard node based.'')')
     &            BASES(IBEG:IEND)
                  WRITE(IFILE,'(5X,''#Nodes='',I2)') NNT(nb)
                  DO nn=1,NNT(nb)
                    WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                    WRITE(IFILE,'(7X,''Value indices:     1'')')
                    WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                  ENDDO
                  NUMFIELDS=NUMFIELDS+1
                ENDIF
                IF(PATH) THEN
                  WRITE(IFILE,'(1X,I1,'') path, field, '
     &              //'rectangular cartesian, #Components=1'')')
     &              NUMFIELDS
                  WRITE(IFILE,'(3X,''path.  '',A,'', no '
     &              //'modify, standard node based.'')')
     &            BASES(IBEG:IEND)
                  WRITE(IFILE,'(5X,''#Nodes='',I2)') NNT(nb)
                  DO nn=1,NNT(nb)
                    WRITE(IFILE,'(5X,I2,''.  #Values=1'')') nn
                    WRITE(IFILE,'(7X,''Value indices:     1'')')
                    WRITE(IFILE,'(7X,''Scale factor indices: '',I3)') nn
                  ENDDO
                  NUMFIELDS=NUMFIELDS+1
                ENDIF

                ELEMENT=1
                SCALEFACTOR=1.0d0
                DO nq=1,NQLIST(0)
                  IF(USE_LAT.EQ.0.AND.(.NOT.LIBMESH)) THEN
                    mq=NQLIST(nq)
                    mq2=NXQ(1,1,mq,1)
                    mq3=NXQ(2,1,mq,1)
                    mq4=NXQ(1,1,NXQ(2,1,mq,1),1)
                    mq5=NXQ(3,1,mq,1)
                    mq6=NXQ(1,1,NXQ(3,1,mq,1),1)
                    mq7=NXQ(2,1,NXQ(3,1,mq,1),1)
                    mq8=NXQ(1,1,NXQ(2,1,NXQ(3,1,mq,1),1),1)
                  ELSEIF(USE_LAT.EQ.0.AND.LIBMESH) THEN
                    mq=NQLIST(nq)
                    mq2=NXQ(1,1,mq,1)
                    mq3=NXQ(1,1,NXQ(2,1,mq,1),1)
                    mq4=NXQ(2,1,mq,1)
                    mq5=NXQ(3,1,mq,1)
                    mq6=NXQ(1,1,NXQ(3,1,mq,1),1)
                    mq7=NXQ(1,1,NXQ(2,1,NXQ(3,1,mq,1),1),1)
                    mq8=NXQ(2,1,NXQ(3,1,mq,1),1)
                  ELSE !Lattice-based grid
                    mq=NQNLAT(NQLIST(nq))
                    mq2=NQNLAT(NQLIST(nq)+1)
                    mq3=NQNLAT(NQLIST(nq)+NQXI(1,SCHEME))
                    mq4=NQNLAT(NQLIST(nq)+1+NQXI(1,SCHEME))
                    mq5=NQNLAT(NQLIST(nq)+
     &                NQXI(1,SCHEME)*NQXI(2,SCHEME))
                    mq6=NQNLAT(NQLIST(nq)+1+
     &                NQXI(1,SCHEME)*NQXI(2,SCHEME))
                    mq7=NQNLAT(NQLIST(nq)+NQXI(1,SCHEME)+
     &                NQXI(1,SCHEME)*NQXI(2,SCHEME))
                    mq8=NQNLAT(NQLIST(nq)+1+NQXI(1,SCHEME)+
     &                NQXI(1,SCHEME)*NQXI(2,SCHEME))
                  ENDIF
                  WRITE(IFILE,'(1X,''Element: '',I6,'' 0 0'' )')
     &              ELEMENT+offset_elem
                  WRITE(IFILE,'(3X,''Faces:'' )')
                  DO nf=1,NFE(nb)
                    WRITE(IFILE,'(5X,''0 '',I6,'' 0'' )')
     &                nf+(NFE(nb)*(ELEMENT-1))+offset_face
                  ENDDO
                  WRITE(IFILE,'(3X,''Nodes:'' )')
                  WRITE(IFILE,'(4X,8(1X,I6))') mq+offset_node,
     &              mq2+offset_node,mq3+offset_node,mq4+offset_node,
     &              mq5+offset_node,mq6+offset_node,mq7+offset_node,
     &              mq8+offset_node
                  WRITE(IFILE,'(3X,''Scale factors:'' )')
                  WRITE(IFILE,'(4X,5(1X,E24.16))') SCALEFACTOR,
     &              SCALEFACTOR,SCALEFACTOR,SCALEFACTOR,SCALEFACTOR
                  WRITE(IFILE,'(2X,5(1X,E24.16))') SCALEFACTOR,
     &              SCALEFACTOR,SCALEFACTOR
                  ELEMENT=ELEMENT+1
                ENDDO

              ENDIF
            ELSE
              CALL ASSERT(.FALSE.,' Incorrect number of xi directions',
     &          ERROR,*9999)
            ENDIF
          ENDDO
          CALL CLOSEF(IFILE,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('EXGRID')
      RETURN
 9999 CALL ERRORS('EXGRID',ERROR)
      CALL EXITS('EXGRID')
      RETURN 1
      END
