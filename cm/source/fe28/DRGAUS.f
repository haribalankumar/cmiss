      SUBROUTINE DRGAUS(ISEG,ISGAUS,MXI,NBH,NBJ,
     '  NEELEM,NELIST,NGLIST,NHE,NHP,
     '  NKH,NKHE,NKJE,NPF,NPNE,NPNODE,NRE,NRLIST,
     '  NVHE,NVHP,NVJE,NW,NXLIST,NYNE,NYNP,
     '  CURVCORRECT,SE,PG,XA,XE,XG,XIG,XP,YP,
     '  ZA,ZE,ZG,ZP,CSEG,STRING,ERROR,*)

C#### Subroutine: DRGAUS
C###  Description:
C###    DRGAUS draws Gauss point positions in segments ISGAUS(iw,ng,ne).

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER ISEG(*),ISGAUS(NWM,NGM,NEM),MXI(2,NEM),NBH(NHM,NCM,NEM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NGLIST(0:NGM),NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),
     '  NKH(NHM,NPM,NCM,0:NRM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NRE(NEM),NRLIST(0:NRM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM),NXLIST(0:NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),PG(NSM,NUM,NGM,NBM),
     '  SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),XIG(NIM,NGM,NBM),
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM),ZA(NAM,NHM,NCM,NEM),
     '  ZE(NSM,NHM),ZG(NHM,NUM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER i,IBEG,IBEG1,IEND,IEND1,INDEX,
     '  INDEX_POLYLINE,iw,IWK(6),j,
     '  n1list,N3CO,nb,nc,ne,ng,nj,noiw,nolist,no_nrlist,nr,
     '  NTIW,nx,nxc
      REAL*8 DXIX(3,3),X(3)
      CHARACTER TYPE*6
      LOGICAL ALL_GAUSS_POINTS,ALL_REGIONS,CBBREV,DEFORM

      CALL ENTERS('DRGAUS',*9999)
      nc=1 ! Temporary MPN 12-Nov-94
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM draw gauss
C###  Parameter:    <(point/axes)[point]>
C###    Draw a point or draw a set of axes at each gauss point.
C###  Parameter:    <number (all/GAUSS#)[all]>
C###  Parameter:    <elements (all/#s/GROUP)[all]>
C###     Define what elements are to be included.
C###  Parameter:    <(undeformed/deformed)[undeformed]>
C###     Draw undeformed/deformed geometry
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the workstation (GX window) to draw the
C###    points on.
C###  Parameter:    <rgb=RGB[green]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Parameter:    <region #[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Draws Gauss point positions, as either points or  axes, in the
C###    specified elements. Gauss points are drawn to the specified
C###    workstation in the specified colour.

        OP_STRING(1)=STRING(1:IEND)//' <(point/axes)[point]>'
        OP_STRING(2)=BLANK(1:15)//'<number (all/GAUSS#)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<elements (all/#s/GROUP)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<(undeformed/deformed)[undeformed]>'
        OP_STRING(5)=BLANK(1:15)//'<on (all/WS#s)[all]>'
        OP_STRING(6)=BLANK(1:15)//'<rgb=RGB[green]>'
        OP_STRING(7)=BLANK(1:15)//'<region #[1]>'
        OP_STRING(8)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRGAUS',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '    ERROR,*9999)

        IF(CBBREV(CO,'POINT',1,noco+1,NTCO,N3CO)) THEN
          TYPE='POINT'
        ELSE IF(CBBREV(CO,'AXES',1,noco+1,NTCO,N3CO)) THEN
          TYPE='AXES'
        ELSE
          TYPE='POINT'
        ENDIF

        IF(CBBREV(CO,'NUMBER',1,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NGM,NGLIST(0),NGLIST(1),ERROR,*9999)
          ALL_GAUSS_POINTS=.FALSE.
        ELSE
          ALL_GAUSS_POINTS=.TRUE.
        ENDIF

        IF(CBBREV(CO,'UNDEFORMED',1,noco+1,NTCO,N3CO)) THEN
          DEFORM=.FALSE.
        ELSE IF(CBBREV(CO,'DEFORMED',1,noco+1,NTCO,N3CO)) THEN
          DEFORM=.TRUE.
        ELSE
          DEFORM=.FALSE.
        ENDIF

        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
        ELSE
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','GREEN')
        ENDIF

        DO i=1,3
          DO j=1,3
            DXIX(i,j)=0.0d0
          ENDDO
        ENDDO

        IF(DEFORM) THEN
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nxc=NXLIST(1)
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '      ERROR,*9999)
          nr=NRLIST(1)
          CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),NKH(1,1,1,nr),
     '      NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,
     '      ERROR,*9999)
        ENDIF

        DO noiw=1,NTIW
          iw=IWK(noiw)
          CALL ACWK(iw,1,ERROR,*9999)
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            DO nolist=1,NELIST(0)
              ne=NELIST(nolist)
              IF(DEFORM) THEN
                CALL ZPZE(NBH(1,1,ne),nc,NHE(ne,nx),NKHE(1,1,1,ne),
     '            NPF(1,1),NPNE(1,1,ne),NRE(ne),NVHE(1,1,1,ne),
     '            NW(ne,1,nx),nx,CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '            ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
              ELSE
                CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '            NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '            SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
              ENDIF
              nb=NBJ(1,ne)
              IF(ALL_GAUSS_POINTS) THEN
                NGLIST(0)=NGT(nb)
                DO ng=1,NGT(nb)
                  NGLIST(ng)=ng
                ENDDO
              ENDIF
              DO n1list=1,NGLIST(0)
                ng=NGLIST(n1list)
                IF(DEFORM) THEN
C MPN 4Apr97 DXIX not used in ZEZG since IP=0.
C GMH 2/9/95 ??? DXIX should be set before call.
                  CALL ZEZG(0,NBH(1,1,ne),ng,NHE(1,nx),nx,DXIX,PG,ZE,ZG,
     '              ERROR,*9999)
                  DO nj=1,NJT
                    X(nj)=ZG(nj,1)
                  ENDDO
                ELSE
                  CALL XEXG(NBJ(1,ne),ng,NRE(ne),PG,
     '              XE,XG,ERROR,*9999)
                  DO nj=1,NJT
                    X(nj)=XG(nj,1)
                  ENDDO
                ENDIF
                CALL SGGAUS(INDEX,ISEG,ISGAUS(iw,ng,ne),iw,MXI(1,ne),
     '            ne,ng,NIT(nb),nr,
     '            CSEG,TYPE,X,XG,XIG(1,ng,nb),DEFORM,ERROR,*9999)
              ENDDO !n1list (ng)
            ENDDO !nolist (ne)
          ENDDO !no_nrlist (nr)
          CALL DAWK(iw,1,ERROR,*9999)
        ENDDO !noiw
      ENDIF

      CALL EXITS('DRGAUS')
      RETURN
 9999 CALL ERRORS('DRGAUS',ERROR)
      CALL EXITS('DRGAUS')
      RETURN 1
      END


