      SUBROUTINE DEMAP(ISEG,ISMAP,MXI,NBJ,NEELEM,NELIST,
     '  NENP,NNB,NPNE,NRLIST,NXI,CSEG,STRING,ERROR,*)

C#### Subroutine: DEMAP
C###  Description:
C###    DEMAP defines map in Cylindrical, Hammer, Polar, Rectangular or
C###    Xi projection.
C###    DEMAP also defines the mapping for parameters in a model.
C###    For example, a material parameter in an equation can be mapped
C###    to a field such the value material parameter is subsituted with
C###    the field value during the solve. This mapping can be defined
C###    for other parameters such as those for CellML models.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'graf00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'map000.cmn'
!     Parameter List
      INTEGER ISEG(*),ISMAP(NGRSEGM),
     '  MXI(2,NEM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NELIST(0:NEM),NENP(NPM,0:NEPM,0:NRM),NNB(4,4,4,NBFM),
     '  NPNE(NNM,NBFM,NEM),NRLIST(0:NRM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),INDEX,INDEX_POLYLINE,N3CO,NBREAKCO,
     '  ne,ne1,NE_BOT_LEFT,NE_BOTTOM,NESTRT,noiw,nr,NTIW
      REAL*8 RFROMC
      CHARACTER DESTINATION*16,DESTINATION_CHILD*255,
     &  DESTINATION_PARENT*255,SOURCE*16,SOURCE_CHILD*255,
     &  SOURCE_PARENT*255
      LOGICAL ALL_REGIONS,CBBREV,CONTINUE

      CALL ENTERS('DEMAP',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM define map xi
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <element (#s/all)[all]>
C###    Specify the element numbers to used. The "all" keyword will
C###    use all currently defined elements in the given regions.
C###  Parameter:      <rgb=RGB[cyan]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Description:
C###

        OP_STRING(1)=STRING(1:IEND)//' xi'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<rgb=RGB[cyan]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define map cylindrical/hammer/polar
C###  Parameter:      <on WS#[4]>
C###    Specify the workstation (GX window) to draw the
C###    (object) on.
C###  Parameter:      <mu=DEGREES#[0.0]{degrees}>
C###  Parameter:      <rotate DEGREES#[0.0]{degrees}>
C###  Parameter:      <label>
C###  Parameter:      <rgb=RGB[black]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Description:
C###    This command creates a new workstation as a projection onto 2D
C###    of the 3D coordinates.  For a Hammer map, the option label will
C###    display the labels LAD,PDA,LV which are relevant only for a
C###    projection of cardiac geometry.  The projection can be rotated
C###    in the theta direction by a user specified number of degrees.

        OP_STRING(1)=STRING(1:IEND)//' cylindrical/hammer/polar'
        OP_STRING(2)=BLANK(1:15)//'<on WS#[4]>'
        OP_STRING(3)=BLANK(1:15)//'<mu=DEGREES#[0.0]{degrees}>'
        OP_STRING(4)=BLANK(1:15)//'<rotate DEGREES#[0.0]{degrees}>'
        OP_STRING(5)=BLANK(1:15)//'<label>'
        OP_STRING(6)=BLANK(1:15)//'<rgb=RGB[black]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define map rectangular
C###  Parameter:      <on WS#[4]>
C###    Specify the workstation (GX window) to draw the
C###    (object) on.
C###  Parameter:      <x_size SIZE#[0.7]>
C###  Parameter:      <y_size SIZE#[0.7]>
C###  Parameter:      <rgb=RGB[black]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Description:
C###    Define a rectangular map on specified workstation with
C###    specified dimensions.

        OP_STRING(1)=STRING(1:IEND)//' rectangular'
        OP_STRING(2)=BLANK(1:15)//'<on WS#[4]>'
        OP_STRING(3)=BLANK(1:15)//'<x_size SIZE#[0.7]>'
        OP_STRING(4)=BLANK(1:15)//'<y_size SIZE#[0.7]>'
        OP_STRING(5)=BLANK(1:15)//'<rgb=RGB[black]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define map
C###  Parameter:         (equation NAME input NAME)
C###    Specify an equation input variable to being mapped.
C###  Parameter:      to
C###    All settings after "to" relate to the destination of the mapping.
C###  Parameter:         (initial_value)
C###    Resets the mapping to the initial value.
C###  Parameter:         (constant NAME component NAME)
C###    Specify the mapping to a constant component.
C###  Parameter:         (field NAME component NAME)
C###    Specify the mapping to a field component.
C###  Parameter:         (maths NAME output NAME)
C###    Specify the mapping to a maths output.
C###  Parameter:         (equation NAME output NAME)
C###    Specify the mapping to an equation output.
C###  Description:
C###    Define a mapping for two parameters in a model.
C###    For example, a material parameter in an equation can be mapped
C###    to a field such the value material parameter is subsituted with
C###    the field value during the solve. This mapping can be defined
C###    for other parameters such as those for CellML models.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'   (maths NAME variable NAME)'
        OP_STRING(3)=BLANK(1:15)//'   (equation NAME variable NAME)'
        OP_STRING(4)=BLANK(1:15)//'to'
        OP_STRING(5)=BLANK(1:15)//'   (initial_value)'
        OP_STRING(6)=BLANK(1:15)//'   (constant NAME component NAME)'
        OP_STRING(7)=BLANK(1:15)//'   (field NAME component NAME)'
        OP_STRING(8)=BLANK(1:15)//'   (maths NAME variable NAME)'
        OP_STRING(9)=BLANK(1:15)//'   (equation NAME variable NAME)'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEMAP',ERROR,*9999)
      ELSE
        nr=1 !needs updating
        
        IF(CBBREV(CO,'TO',2,noco+1,NTCO,NBREAKCO))THEN

          IF(CBBREV(CO,'MATHS',4,noco+1,NBREAKCO,N3CO))THEN
            SOURCE='MATHS'
            SOURCE_PARENT=CO(N3CO+1)
            IF(CBBREV(CO,'VARIABLE',3,noco+1,NBREAKCO,N3CO))THEN
              SOURCE_CHILD=CO(N3CO+1)
            ELSE
              CALL ASSERT(.TRUE.,'>>Maths input variable not defined',
     &          ERROR,*9999)
            ENDIF
          ELSEIF(CBBREV(CO,'EQUATION',4,noco+1,NBREAKCO,N3CO))THEN
            SOURCE='EQUATION'
            SOURCE_PARENT=CO(N3CO+1)
            IF(CBBREV(CO,'VARIABLE',3,noco+1,NBREAKCO,N3CO))THEN
              SOURCE_CHILD=CO(N3CO+1)
            ELSE
              CALL ASSERT(.TRUE.,
     &          '>>Equation input variable not defined',ERROR,*9999)
            ENDIF
          ELSE
            CALL ASSERT(.TRUE.,'>>Source type not defined',
     &        ERROR,*9999)
          ENDIF

          IF(CBBREV(CO,'INITIAL_VALUE',4,NBREAKCO+1,NTCO,N3CO))THEN
            DESTINATION='INITIAL_VALUE'
            DESTINATION_PARENT=' '
            DESTINATION_CHILD=' '
          ELSEIF(CBBREV(CO,'CONSTANT',5,NBREAKCO+1,NTCO,N3CO))THEN
            DESTINATION='CONSTANT'
            DESTINATION_PARENT=CO(N3CO+1)
            IF(CBBREV(CO,'COMPONENT',4,NBREAKCO+1,NTCO,N3CO))THEN
              DESTINATION_CHILD=CO(N3CO+1)
            ELSE
              CALL ASSERT(.TRUE.,'>>Constant component not defined',
     &          ERROR,*9999)
            ENDIF
          ELSEIF(CBBREV(CO,'FIELD',5,NBREAKCO+1,NTCO,N3CO))THEN
            DESTINATION='FIELD'
            DESTINATION_PARENT=CO(N3CO+1)
            IF(CBBREV(CO,'COMPONENT',4,NBREAKCO+1,NTCO,N3CO))THEN
              DESTINATION_CHILD=CO(N3CO+1)
            ELSE
              CALL ASSERT(.TRUE.,'>>Field component not defined',
     &          ERROR,*9999)
            ENDIF
          ELSEIF(CBBREV(CO,'MATHS',4,NBREAKCO+1,NTCO,N3CO))THEN
            DESTINATION='MATHS'
            DESTINATION_PARENT=CO(N3CO+1)
            IF(CBBREV(CO,'VARIABLE',3,NBREAKCO+1,NTCO,N3CO))THEN
              DESTINATION_CHILD=CO(N3CO+1)
            ELSE
              CALL ASSERT(.TRUE.,'>>Maths output variable not defined',
     &          ERROR,*9999)
            ENDIF
          ELSEIF(CBBREV(CO,'EQUATION',4,NBREAKCO+1,NTCO,N3CO))THEN
            DESTINATION='EQUATION'
            DESTINATION_PARENT=CO(N3CO+1)
            IF(CBBREV(CO,'VARIABLE',3,NBREAKCO+1,NTCO,N3CO))THEN
              DESTINATION_CHILD=CO(N3CO+1)
            ELSE
              CALL ASSERT(.TRUE.,
     &          '>>Equation output variable not defined',ERROR,*9999)
            ENDIF
          ELSE
            CALL ASSERT(.TRUE.,'>>Destination type not defined',
     &        ERROR,*9999)
          ENDIF
          
          IF(SOURCE.EQ.'MATHS')THEN
            CALL SET_MATHS_INPUT_MAP(SOURCE_PARENT,SOURCE_CHILD,
     &        DESTINATION,DESTINATION_PARENT,DESTINATION_CHILD,
     &        ERROR,*9999)
          ELSEIF(SOURCE.EQ.'EQUATION')THEN
            CALL SET_EQUATION_INPUT_MAP(SOURCE_PARENT,SOURCE_CHILD,
     &        DESTINATION,DESTINATION_PARENT,DESTINATION_CHILD,
     &        ERROR,*9999)
          ENDIF
            
        ELSE

          IF(CBBREV(CO,'CYLINDRICAL',1,noco+1,NTCO,N3CO)) THEN
            MAP_PROJEC='CYLINDRICAL'
            CALL ASSERT(NJT.EQ.3,'>>Map only for 3D',ERROR,*9999)
          ELSE IF(CBBREV(CO,'HAMMER',1,noco+1,NTCO,N3CO)) THEN
            MAP_PROJEC='HAMMER'
            CALL ASSERT(NJT.EQ.3,'>>Map only for 3D',ERROR,*9999)
          ELSE IF(CBBREV(CO,'POLAR',1,noco+1,NTCO,N3CO)) THEN
            MAP_PROJEC='POLAR'
            CALL ASSERT(NJT.EQ.3,'>>Map only for 3D',ERROR,*9999)
          ELSE IF(CBBREV(CO,'RECTANGULAR',1,noco+1,NTCO,N3CO)) THEN
            MAP_PROJEC='RECTANGULAR'
            CALL ASSERT(ITYP10(1).EQ.1,
     '        '>>Rectangular map only for Cartesian',ERROR,*9999)
          ELSE IF(CBBREV(CO,'XI',1,noco+1,NTCO,N3CO)) THEN
            MAP_PROJEC='XI'
          ELSE
            MAP_PROJEC='XI'
          ENDIF
          CALL WS_LIST(IWK,4,NTIW,noco,NTCO,CO,ERROR,*9999)

          IF(MAP_PROJEC(1:2).EQ.'XI') THEN
            CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,
     '        ERROR,*9999)
            CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '        ERROR,*9999)
            IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
              INDEX=INDEX_POLYLINE(0,'DOTTED','WIDTH1',CO(N3CO+1))
            ELSE
              INDEX=INDEX_POLYLINE(0,'DOTTED','WIDTH1','CYAN')
            ENDIF
          ELSE
            IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
              INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
            ELSE
              INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLACK')
            ENDIF
          ENDIF

          DO noiw=1,NTIW
            iw=IWK(noiw)
            CALL ASSERT(iw.EQ.4.OR.iw.EQ.15,'>>Map only for 4 or 15',
     '        ERROR,*9999)
          ENDDO

          IF(MAP_PROJEC(1:11).EQ.'CYLINDRICAL') THEN

          ELSE IF(MAP_PROJEC(1:6).EQ.'HAMMER') THEN
            IF(CBBREV(CO,'MU',1,noco+1,NTCO,N3CO)) THEN
              RMU=RFROMC(CO(N3CO+1))
            ELSE
              RMU=-1.0D0
            ENDIF
            IF(CBBREV(CO,'ROTATION',1,noco+1,NTCO,N3CO)) THEN
              RMAP=RFROMC(CO(N3CO+1))
            ELSE
              RMAP=0.0D0
            ENDIF
            IF(CBBREV(CO,'LABEL',1,noco+1,NTCO,N3CO)) THEN
              LABEL=.TRUE.
            ELSE
              LABEL=.FALSE.
            ENDIF

          ELSE IF(MAP_PROJEC(1:5).EQ.'POLAR') THEN
            RMAP=DSQRT(DBLE(YMAX)**2+DBLE(ZMAX)**2)

          ELSE IF(MAP_PROJEC(1:11).EQ.'RECTANGULAR') THEN
            IF(CBBREV(CO,'X_SIZE',1,noco+1,NTCO,N3CO)) THEN
              XI1MAP=RFROMC(CO(N3CO+1))
            ELSE
              XI1MAP=0.7D0
            ENDIF
            IF(CBBREV(CO,'Y_SIZE',1,noco+1,NTCO,N3CO)) THEN
              XI2MAP=RFROMC(CO(N3CO+1))
            ELSE
              XI2MAP=0.7D0
            ENDIF

          ELSE IF(MAP_PROJEC(1:2).EQ.'XI') THEN
            CALL NENXI(NBJ,NEELEM,NENP,NNB,NPNE,NXI,ERROR,*9999)
            ne=NEELEM(1,nr)
            ne1=ne
            NESTRT=0
            DO WHILE(NXI(-1,1,ne).gt.0.and.NE.NE.NESTRT) !left of 1st row
              NESTRT=ne1
              ne=NXI(-1,1,ne)
              WRITE(OP_STRING,'('' ne='',i3)') ne
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDDO
            DO WHILE(NXI(-2,1,ne).GT.0)   !bottom of 1st col
              ne=NXI(-2,1,ne)
            ENDDO
            NE_BOT_LEFT=ne              !is bottom left element
            IF(DOP) THEN
              WRITE(OP_STRING,'('' NE_BOT_LEFT='',I3)') NE_BOT_LEFT
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            NE_BOTTOM=NE_BOT_LEFT
            NET_XI1=1
            CONTINUE=.TRUE.
            DO WHILE(CONTINUE)    !loops over columns
              IF(DOP) THEN
                WRITE(OP_STRING,'('' ne_bottom='',I3)') NE_BOTTOM
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              ne=NE_BOTTOM
              NET_XI2=1
              MXI(1,ne)=NET_XI1
              MXI(2,ne)=NET_XI2
              IF(DOP) THEN
                WRITE(OP_STRING,'('' ne='',I3,'' MXI(1..,ne)='',2I3)')
     '            ne,MXI(1,ne),MXI(2,ne)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              DO WHILE(NXI(2,1,ne).GT.0)  !go to top of current col
                ne=NXI(2,1,ne)
                NET_XI2=NET_XI2+1 !..counting no of elements in Xi2 dir.n
                MXI(1,ne)=NET_XI1 !..recording relative Xi1 position of ne
                MXI(2,ne)=NET_XI2 !..recording relative Xi2 position of ne
                IF(DOP) THEN
                  WRITE(OP_STRING,
     '              '('' ne='',I3,'' MXI(1..,ne)='',2I3)')
     '              ne,MXI(1,ne),MXI(2,ne)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
              ENDDO
              IF(NXI(1,1,NE_BOTTOM).GT.0) THEN
                NET_XI1=NET_XI1+1  !..counting no of elements in Xi1 dir.n
                NE_BOTTOM=NXI(1,1,NE_BOTTOM) !for next bottom to the right
              ELSE
                CONTINUE=.FALSE.
              ENDIF
            ENDDO
            IF(DOP) THEN
              WRITE(OP_STRING,'('' NET_XI1='',I3,'' NET_XI2='',I3)')
     '          NET_XI1,NET_XI2
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            MAX_XI=MAX(NET_XI1,NET_XI2)
          ENDIF

          DO noiw=1,NTIW
            iw=IWK(noiw)
            CALL ACWK(iw,0,ERROR,*9999)
            IMAP=1
c PJH 4-Jan-1991 NTMAP=NTMAP+1
            NTMAP=1
            CALL ASSERT(NTMAP.LE.NRM,'>>NRM too small',ERROR,*9999)
            CALL SGMAP(INDEX,ISEG,ISMAP(NTMAP),iw,CSEG,ERROR,*9999)
            CALL DAWK(iw,0,ERROR,*9999)
          ENDDO
        ENDIF
      ENDIF

      CALL EXITS('DEMAP')
      RETURN
 9999 CALL ERRORS('DEMAP',ERROR)
      CALL EXITS('DEMAP')
      RETURN 1
      END


