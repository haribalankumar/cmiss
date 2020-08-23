      SUBROUTINE DRSURF(IBT,IDO,INP,ISEG,ISSURF,NAN,NBH,NBJ,
     '  NEELEM,NELIST,NHE,NHP,NKH,NKHE,NKJE,NPF,NPNE,NPNODE,NRE,
     '  NVHE,NVHP,NVJE,NW,NYNE,NYNP,
     '  CURVCORRECT,PG,RG,SE,XA,XE,XG,XP,YP,ZA,ZE,ZG,ZP,CSEG,
     '  STRING,ERROR,*)

C#### Subroutine: DRSURF
C###  Description:
C###    DRSURF draws element surface.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grou00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'moti00.cmn'
      INCLUDE 'ntsg00.cmn'
      INCLUDE 'scal00.cmn'
      INCLUDE 'surf00.cmn'
      INCLUDE 'time00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISEG(*),ISSURF(NWM,NGRSEGM),NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),NHE(NEM,NXM),
     '  NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKJE(NKM,NNM,NJM,NEM),NKHE(NKM,NNM,NHM,NEM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),NRE(NEM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),PG(NSM,NUM,NGM,NBM),RG(NGM),
     '  SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM),ZA(NAM,NHM,NCM,NEM),
     '  ZE(NSM,NHM),ZG(NHM,NUM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,INDEX,INDEX_POLYLINE,iw,IWK(6),n1grel,N3CO,
     '  noelem,nogrel,noiw,NOXIPT(2),NO_SURF(2),np,nr,nt,NTIW,NTRL,
     '  NTXIPT,nx
      REAL*8 RL(2),XI3
      CHARACTER CHAR*30,LABEL*30
      LOGICAL CBBREV,COLOUR,DEFORM,STATIC

      CALL ENTERS('DRSURF',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM draw surface
C###  Parameter:    <(geometry/field/strain)[geometry]>
C###    Specifies wiether the surface is drawn from geometry,
C###    field, or strain information.
C###  Parameter:    <(P1/P2/PA/E11/E22/E33/E12/E13/E23)[P1]>
C###    Specifies the tensor component with which to draw thw
C###    surface with
C###  Parameter:    <rgb=RGB[black]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Parameter:    <(greyscale/colour)[greyscale]>
C###    Specify whether the surface is to be colour or greyscale.
C###  Parameter:    <(dots/lines/pattern)[pattern]>
C###     Draw the surface as dots, lines or a pattern
C###  Parameter:    <(undeformed/deformed)[undeformed]>
C###     Draw undeformed/deformed geometry
C###  Parameter:    <in (all/ELEMENT#s/GROUP_LABEL#)[all]>
C###    Specify the element numbers or group labels within which to
C###    draw the surface in.
C###  Parameter:    <at XI_3#[0.0]>
C###    Distance along Xi_3 direction where surface will be drawn.
C###  Parameter:    <when (static/TIME)[static]>
C###    Specify wiether the surface is static or the time
C###  Parameter:    <by NO_PT#s[10,10]>
C###    Specifies the number of points
C###  Parameter:    <on WS_ID#[3]>
C###    Specify the workstation (GX window) to draw the
C###    surface on.
C###  Parameter:    <surface_index (1/2/3/4/5)[4]>
C###    Specify the surface index
C###  Parameter:    <range MIN#[0.0],MAX#[1.0]>
C###    Specifies the range of the surface
C###  Description:
C###    Draw element surfaces.
CC###  Parameter:    <depthcue (on/off)[off]>
CC###    Choose to use depth cue
CC###  Parameter:    <edges (on/off)[on]>
CC###    Choose to use edges
CC###  Parameter:    <cull (front/back/none)[none]>

        OP_STRING( 1)=STRING(1:IEND)
     '    //' <(geometry/field/strain)[geometry]>'
        OP_STRING( 2)=BLANK(1:15)
     '    //'<(P1/P2/PA/E11/E22/E33/E12/E13/E23)[P1]>'
        OP_STRING( 3)=BLANK(1:15)//'<rgb=RGB[black]>'
        OP_STRING( 4)=BLANK(1:15)//'<(greyscale/colour)[greyscale]>'
        OP_STRING( 5)=BLANK(1:15)//'<(dots/lines/pattern)[pattern]>'
        OP_STRING( 6)=BLANK(1:15)//'<(undeformed/deformed)[undeformed]>'
        OP_STRING( 7)=BLANK(1:15)
     '  //'<in (all/ELEMENT#s/of GROUP_LABEL#)[all]>'
        OP_STRING( 8)=BLANK(1:15)//'<at Xi_3#[0.0]>'
        OP_STRING( 9)=BLANK(1:15)//'<when (static/TIME)[static]>'
        OP_STRING(10)=BLANK(1:15)//'<by NO_PT#s[10,10]>'
        OP_STRING(11)=BLANK(1:15)//'<on WS_ID#[3]>'
        OP_STRING(12)=BLANK(1:15)//'<surface_index (1/2/3/4/5)[4]>'
C        OP_STRING(13)=BLANK(1:15)//'<depthcue (on/off)[off]>'
C        OP_STRING(14)=BLANK(1:15)//'<edges (on/off)[on]>'
C        OP_STRING(15)=BLANK(1:15)//'<cull (front/back/none)[none]>'
        OP_STRING(13)=BLANK(1:15)//'<range MIN#[0.0],MAX#[1.0]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRSURF',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        nx=1 ! Temporary cpb 22/11/94
        CALL WS_LIST(IWK,3,NTIW,noco,NTCO,CO,ERROR,*9999)
        IF(CBBREV(CO,'GEOMETRY',2,noco+1,NTCO,N3CO)) THEN
          VARIABLE_TYPE='GEOMETRY'
        ELSE IF(CBBREV(CO,'FIELD',2,noco+1,NTCO,N3CO)) THEN
          VARIABLE_TYPE='FIELD'
        ELSE IF(CBBREV(CO,'STRAIN',2,noco+1,NTCO,N3CO)) THEN
          VARIABLE_TYPE='STRAIN'
          IF(CBBREV(CO,'P1',2,noco+1,NTCO,N3CO)) THEN
            VARIABLE_TYPE(7:8)='P1'
          ELSE IF(CBBREV(CO,'P2',2,noco+1,NTCO,N3CO)) THEN
            VARIABLE_TYPE(7:8)='P2'
          ELSE IF(CBBREV(CO,'PA',2,noco+1,NTCO,N3CO)) THEN
            VARIABLE_TYPE(7:8)='PA'
          ELSE IF(CBBREV(CO,'E11',3,noco+1,NTCO,N3CO)) THEN
            VARIABLE_TYPE(7:9)='E11'
          ELSE IF(CBBREV(CO,'E22',3,noco+1,NTCO,N3CO)) THEN
            VARIABLE_TYPE(7:9)='E22'
          ELSE IF(CBBREV(CO,'E33',3,noco+1,NTCO,N3CO)) THEN
            VARIABLE_TYPE(7:9)='E33'
          ELSE IF(CBBREV(CO,'E12',3,noco+1,NTCO,N3CO)) THEN
            VARIABLE_TYPE(7:9)='E12'
          ELSE IF(CBBREV(CO,'E13',3,noco+1,NTCO,N3CO)) THEN
            VARIABLE_TYPE(7:9)='E13'
          ELSE IF(CBBREV(CO,'E23',3,noco+1,NTCO,N3CO)) THEN
            VARIABLE_TYPE(7:9)='E23'
          ELSE
            VARIABLE_TYPE(7:8)='P1'
          ENDIF
        ELSE
          VARIABLE_TYPE='GEOMETRY'
        ENDIF

        IF(CBBREV(CO,'COLOUR',1,noco+1,NTCO,N3CO)) THEN
          COLOUR=.TRUE.
        ELSE
          COLOUR=.FALSE.
        ENDIF
        IF(CBBREV(CO,'DOTS',3,noco+1,NTCO,N3CO)) THEN
          SURFACE_TYPE='DOTS'
        ELSE IF(CBBREV(CO,'LINES',3,noco+1,NTCO,N3CO)) THEN
          SURFACE_TYPE='LINES'
        ELSE IF(CBBREV(CO,'PATTERN',3,noco+1,NTCO,N3CO)) THEN
          SURFACE_TYPE='PATTERN'
        ELSE
          SURFACE_TYPE='PATTERN'
        ENDIF
        IF(CBBREV(CO,'UNDEFORMED',3,noco+1,NTCO,N3CO)) THEN
          DEFORM=.FALSE.
        ELSE IF(CBBREV(CO,'DEFORMED',3,noco+1,NTCO,N3CO)) THEN
          DEFORM=.TRUE.
        ELSE
          DEFORM=.FALSE.
        ENDIF
        IF(CBBREV(CO,'IN',2,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NE_R_M,NELIST(0),NELIST(1),ERROR,*9999)
        ELSE IF(CBBREV(CO,'OF',2,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
C          CHAR=CUPPER(CO(N3CO+1)(IBEG:IEND))
          CALL CUPPER(CO(N3CO+1)(IBEG:IEND),CHAR)
          n1grel=0
          DO nogrel=1,NTGREL
C            LABEL=CUPPER(LAGREL(nogrel))
            CALL CUPPER(LAGREL(nogrel),LABEL)
            IF(CHAR(IBEG:IEND).EQ.LABEL(IBEG:IEND)) THEN
              n1grel=nogrel
              GO TO 100
            ENDIF
          ENDDO
 100      IF(n1grel.GT.0) THEN
C KAT 4Nov99 Updating for dynamic groups
            NELIST(0)=NLIGREL(n1grel)
            CALL ILIST_COPY(NELIST(0),
     '        %VAL(LIGREL_PTR(nogrel)),NELIST(1))
C            NELIST(0)=LIGREL(0,n1grel)
C            DO nolist=1,NELIST(0)
C              NELIST(nolist)=LIGREL(nolist,n1grel)
C            ENDDO
          ELSE
            WRITE(OP_STRING,'('' >>No such label defined'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ELSE
          NELIST(0)=0
          DO nr=1,NRT
            DO noelem=NELIST(0)+1,NELIST(0)+NEELEM(0,nr)
              NELIST(noelem)=NEELEM(noelem,nr)
            ENDDO
            NELIST(0)=NELIST(0)+NEELEM(0,nr)
          ENDDO
        ENDIF
        IF(CBBREV(CO,'AT',2,noco+1,NTCO,N3CO)) THEN
          CALL PARSRE(CO(N3CO+1),XI3,ERROR,*9999)
        ELSE
          XI3=0.0D0
        ENDIF
        IF(CBBREV(CO,'WHEN',2,noco+1,NTCO,N3CO)) THEN
          CALL PARSRE(CO(N3CO+1),TIME,ERROR,*9999)
          STATIC=.FALSE.
        ELSE
          STATIC=.TRUE.
          TIME=0.0D0
        ENDIF
        IF(CBBREV(CO,'BY',2,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),2,NTXIPT,NOXIPT,ERROR,*9999)
        ELSE
          NOXIPT(1)=10
          NOXIPT(2)=10
        ENDIF
        IF(CBBREV(CO,'SURFACE_INDEX',3,noco+1,NTCO,N3CO)) THEN
          NO_SURF(1)=INDEX_SURF
          CALL PARSIL(CO(N3CO+1),1,nt,NO_SURF(1),ERROR,*9999)
        ELSE
          IF(VARIABLE_TYPE(1:8).EQ.'GEOMETRY') THEN
            INDEX_SURF=4
          ELSE IF(VARIABLE_TYPE(1:5).EQ.'FIELD') THEN
            INDEX_SURF=2
          ELSE IF(VARIABLE_TYPE(1:6).EQ.'STRAIN') THEN
            INDEX_SURF=2
          ENDIF
        ENDIF
C        IF(CBBREV(CO,'DEPTHCUE',3,noco+1,NTCO,N3CO)) THEN
C          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
C          DEPTHCUE=CO(N3CO+1)(IBEG:IEND)
C        ELSE
C          DEPTHCUE='OFF'
C        ENDIF
C        IF(CBBREV(CO,'EDGES',3,noco+1,NTCO,N3CO)) THEN
C          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
C          EDGES=CO(N3CO+1)(IBEG:IEND)
C        ELSE
C          IF(INDEX_SURF.EQ.2.AND.VARIABLE_TYPE(1:8).EQ.'GEOMETRY') THEN
C            EDGES='ON'
C          ELSE
C            EDGES='OFF'
C          ENDIF
C        ENDIF
C        IF(CBBREV(CO,'CULL',3,noco+1,NTCO,N3CO)) THEN
C          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
C          CULL=CO(N3CO+1)(IBEG:IEND)
C        ELSE
C          CULL='NONE'
C        ENDIF
        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
        ELSE
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLACK')
        ENDIF

        IF(VARIABLE_TYPE(1:5).EQ.'FIELD') THEN
          ZMINI=XP(1,1,NJ_LOC(NJL_FIEL,1,nr),1)
          ZMAXI=XP(1,1,NJ_LOC(NJL_FIEL,1,nr),1)
          DO np=2,NPT(1)
            IF(XP(1,1,NJ_LOC(NJL_FIEL,1,nr),np).LT.ZMINI) ZMINI=
     '        XP(1,1,NJ_LOC(NJL_FIEL,1,nr),np)
            IF(XP(1,1,NJ_LOC(NJL_FIEL,1,nr),np).GT.ZMAXI) ZMAXI=
     '        XP(1,1,NJ_LOC(NJL_FIEL,1,nr),np)
          ENDDO
          ZDIFF=ZMAXI-ZMINI
          IF(DABS(ZDIFF).LT.1.0D-6) ZDIFF=1.0D0
          IF(DOP) THEN
            WRITE(OP_STRING,'('' ZMINI='',E12.3,'' ZMAXI='',E12.3)')
     '        ZMINI,ZMAXI
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ELSE IF(VARIABLE_TYPE(1:6).EQ.'STRAIN') THEN
          IF(CBBREV(CO,'RANGE',3,noco+1,NTCO,N3CO)) THEN
            CALL PARSRL(CO(N3CO+1),2,NTRL,RL,ERROR,*9999)
            ZMINI=RL(1)
            ZMAXI=RL(2)
            ZDIFF=ZMAXI-ZMINI
            IF(DOP) THEN
              WRITE(OP_STRING,
     '          '('' ZMINI='',E12.3,'' ZMAXI='',E12.3)')
     '          ZMINI,ZMAXI
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ELSE
            ZMINI=0.0D0
            ZMAXI=1.0D0
            ZDIFF=1.0D0
          ENDIF
        ENDIF

        IF(ADD) THEN
          NTSURF=NTSURF+1
        ELSE IF(NTSURF.EQ.0) THEN
          NTSURF=1
        ENDIF
        CALL ASSERT(NTSURF.LE.NRM,'>>NRM too small',ERROR,*9999)
        IF(DEFORM) THEN
          nr=1 !Needs fixing
          CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '      NKH(1,1,1,nr),NPNODE,
     '      nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)
        ENDIF
        IF(.NOT.STATIC) THEN
          CALL YPZP(MOTION_IY,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '      NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,
     '      YP(1,1,nx),ZA,ZP,ERROR,*9999)
        ENDIF
        DO noiw=1,NTIW
          IW=IWK(noiw)
          CALL ACWK(iw,1,ERROR,*9999)
          CALL SGSURF(INDEX,IBT,IDO,INP,ISEG,ISSURF(iw,NTSURF),iw,NAN,
     '      NBH,NBJ,NELIST,NHE(1,nx),NKHE,NKJE,NTSURF,NOXIPT,NPF,
     '      NPNE,NRE,NVHE,NVJE,NW(1,1,nx),nx,
     '      CURVCORRECT,PG,RG,SE,XA,XE,XG,XI3,XP,ZA,ZE,ZG,ZP,
     '      COLOUR,DEFORM,STATIC,CSEG,ERROR,*9999)
          IF(VARIABLE_TYPE(1:5).EQ.'FIELD') THEN
            NTSCAL=NTSCAL+1
C           CALL SGSCAL(INDEX,ISEG,ISSCAL(iw,NTSCAL),iw,NTSCAL,COLOUR,
C    '        CSEG,ERROR,*9999)
          ENDIF
          CALL DAWK(iw,1,ERROR,*9999)
        ENDDO
      ENDIF

      CALL EXITS('DRSURF')
      RETURN
 9999 CALL ERRORS('DRSURF',ERROR)
      CALL EXITS('DRSURF')
      RETURN 1
      END


