      SUBROUTINE DRSTRA(IBT,IDO,INP,ISEG,ISSTRA,LD,NAN,NBH,NBJ,
     '  NEELEM,NELIST,NHE,NHP,NKH,NKHE,NKJE,NPF,NPNE,NPNODE,
     '  NRE,NRLIST,NVHE,NVHP,NVJE,NW,NXLIST,NYNE,NYNP,
     '  CURVCORRECT,PG,RGX,SE,XA,XE,XG,XID,XIG,XP,YP,
     '  ZA,ZE,ZG,ZP,CSEG,STRING,ERROR,*)

C#### Subroutine: DRSTRA
C###  Description:
C###    DRSTRA draws principal strains for vector plots.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbst02.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'moti00.cmn'
      INCLUDE 'ntsg00.cmn'
      INCLUDE 'time00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISEG(*),ISSTRA(NEM,NGRSEGM),LD(NDM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKJE(NKM,NNM,NJM,NEM),NKHE(NKM,NNM,NHM,NEM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NRE(NEM),NRLIST(0:NRM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM),NXLIST(0:NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),PG(NSM,NUM,NGM,NBM),RGX(NGM),
     '  SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XG(NJM,NUM),XID(NIM,NDM),
     '  XIG(NIM,NGM,NBM),XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM),
     '  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZG(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,INDEX,INDEX_POLYLINE,IPOINTTYP,
     '  ISSTRA_index,iw,IWK(6),N3CO,ne,noelem,noiw,no_nelist,
     '  no_nrlist,nostra,nr,nr1,NTIW,nx,nxc
      REAL*8 RFROMC
      CHARACTER TYPE*8
      LOGICAL ALL_REGIONS,CBBREV,STATIC

      CALL ENTERS('DRSTRA',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM draw strain vector
C###  Parameter:    <(gauss/data)[gauss]>
C###  Parameter:    <in (all/ELEMENT#s)[all]>
C###    Specify the element(s) in which to draw the strain vectors.
C###    The "all" keyword specifies all currently defined elements.
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the worksation(s). The all keyword specifies all
C###    currently defined workstations.
C###  Parameter:    <scale FACTOR#[1.0]>
C###    Specifes the scale factor for drawing the strain vectors.
C###  Parameter:    <rgb=RGB[black]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Draws principal strains for vector plots.

        OP_STRING(1)=STRING(1:IEND)//' vector'
        OP_STRING(2)=BLANK(1:15)//'<(gauss/data)[gauss]>'
        OP_STRING(3)=BLANK(1:15)//'<in (all/ELEMENT#s)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<on (all/WS#s)[all]>'
        OP_STRING(5)=BLANK(1:15)//'<scale FACTOR#[1.0]>'
        OP_STRING(6)=BLANK(1:15)//'<rgb=RGB[black]>'
        OP_STRING(7)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(8)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM draw strain
C###  Parameter:    <number (next/N#)[next]>
C###  Parameter:    <in (all/ELEMENT#s)[all]>
C###    Specify the element(s) in which to draw the strain vectors.
C###    The "all" keyword specifies all the currentlky defined
C###    elements.
C###  Parameter:    <at XI3#[0]{>=0.0,<=1.0}>
C###  Parameter:    <when (static/TIME)[static]>
C###  Parameter:    <wrt TIME_REF#[0.0]>
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the worksation(s). The all keyword specifies all
C###    currently defined workstations.
C###  Parameter:    <scale FACTOR#[1.0]>
C###    Specify the scale factor for drawing the strain vectors.
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Draws principal strain vectors.

        OP_STRING(1)=STRING(1:IEND)//' <number (next/N#)[next]>'
        OP_STRING(2)=BLANK(1:15)//'<in (all/ELEMENT#s)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<at XI3#[0]{>=0.0,<=1.0}>'
        OP_STRING(4)=BLANK(1:15)//'<when (static/TIME)[static]>'
        OP_STRING(5)=BLANK(1:15)//'<wrt TIME_REF#[0.0]>'
        OP_STRING(6)=BLANK(1:15)//'<on (all/WS#s)[all]>'
        OP_STRING(7)=BLANK(1:15)//'<scale FACTOR#>[1.0]>'
        OP_STRING(8)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(9)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRSTRA',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        IF(CBBREV(CO,'VECTOR',1,noco+1,NTCO,N3CO)) THEN
          TYPE='VECTOR'
        ELSE IF(CBBREV(CO,'SURFACE',2,noco+1,NTCO,N3CO)) THEN
          TYPE='SURFACE'
        ELSE
          TYPE='SURFACE'
        ENDIF

        IF(CBBREV(CO,'DATA',1,noco+1,NTCO,N3CO)) THEN
          IPOINTTYP=2 !data points
        ELSE
          IPOINTTYP=1 !Gauss Points
        ENDIF

        IF(CBBREV(CO,'IN',1,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NE_R_M,NELIST(0),NELIST(1),ERROR,*9999)
        ELSE
          NELIST(0)=0
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            DO noelem=1,NEELEM(0,nr)
              NELIST(noelem+NELIST(0))=NEELEM(noelem,nr)
            ENDDO !noelem
            NELIST(0)=NELIST(0)+NEELEM(0,nr)
          ENDDO !no_nrlist (nr)
        ENDIF

        IF(CBBREV(CO,'SCALE',2,noco+1,NTCO,N3CO)) THEN
          SCALE=RFROMC(CO(N3CO+1))
        ELSE
          SCALE=1.0D0
        ENDIF

        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
        ELSE
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLACK')
        ENDIF

        IF(TYPE(1:6).EQ.'VECTOR') THEN
          STATIC=.TRUE.
          IF(ADD) THEN
            NTSTRA=NTSTRA+1
          ELSE IF(NTSTRA.EQ.0) THEN
            NTSTRA=1
          ENDIF
          CALL ASSERT(NTSTRA.LE.NRM,'>>NRM too small',ERROR,*9999)

        ELSE
          IF(CBBREV(CO,'WHEN',2,noco+1,NTCO,N3CO)) THEN
            TIME=RFROMC(CO(N3CO+1))
            STATIC=.FALSE.
          ELSE
            STATIC=.TRUE.
          ENDIF
          IF(CBBREV(CO,'WRT',2,noco+1,NTCO,N3CO)) THEN
            TIME_REF=RFROMC(CO(N3CO+1))
          ELSE
            TIME_REF=0.0D0
          ENDIF
          IF(CBBREV(CO,'NUMBER',2,noco+1,NTCO,N3CO)) THEN
            nostra=IFROMC(CO(N3CO+1))
            NTSTRA=MAX(NTSTRA,nostra)
          ELSE
            NTSTRA=NTSTRA+1
            nostra=NTSTRA
          ENDIF
          CALL ASSERT(NTSTRA.LE.NRM,'>>NRM too small',ERROR,*9999)
        ENDIF

        nr1=0
        DO no_nelist=1,NELIST(0)
          ne=NELIST(no_nelist)
          nr=NRE(ne)
          IF(nr.NE.nr1) THEN
            IF(TYPE(1:6).EQ.'VECTOR') THEN
              CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '          NKH(1,1,1,nr),NPNODE,
     '          nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,
     '          ERROR,*9999)
              ISSTRA_index=NTSTRA
            ELSE
              IF(.NOT.STATIC) THEN !motion field in YP
                CALL YPZP(MOTION_IY,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '            NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,
     '            YP(1,1,nx),ZA,ZP,ERROR,*9999)
              ENDIF
              ISSTRA_index=nostra
            ENDIF
            nr1=nr
          ENDIF !nr.NE.nr1
          DO noiw=1,NTIW
            iw=IWK(noiw)
            CALL ACWK(iw,1,ERROR,*9999)
            CALL SGSTRA(INDEX,IBT,IDO,INP,IPOINTTYP,ISEG,
     '        ISSTRA(1,ISSTRA_index),
     '        iw,LD,NAN,NBH(1,1,ne),NBJ(1,ne),ne,NHE(ne,nx),
     '        NKHE(1,1,1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),nr,
     '        NVHE(1,1,1,ne),NVJE(1,1,1,ne),NW(ne,1,nx),nx,
     '        CURVCORRECT(1,1,1,ne),PG,RGX,SE(1,1,ne),STATIC,
     '        XA,XE,XG,XID,XIG,XP,ZA(1,1,1,ne),ZE,ZG,ZP,CSEG,
     '        ERROR,*9999)
            CALL DAWK(iw,1,ERROR,*9999)
          ENDDO !noiw (iw)
        ENDDO !no_nelist (ne)
      ENDIF

      CALL EXITS('DRSTRA')
      RETURN
 9999 CALL ERRORS('DRSTRA',ERROR)
      CALL EXITS('DRSTRA')
      RETURN 1
      END



