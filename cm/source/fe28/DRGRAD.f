      SUBROUTINE DRGRAD(ISEG,ISGRAD,NBH,NBJ,NEELEM,NELIST,NHE,NHP,
     '  NKH,NKHE,NKJE,NPF,NPNE,NPNODE,NRE,NRLIST,NVHE,NVHP,NVJE,NW,
     '  NXLIST,NYNE,NYNP,CE,CG,CGE,CP,CURVCORRECT,PG,RG,SE,XA,XE,XG,
     '  XP,YP,ZA,ZE,ZG,ZP,CSEG,STRING,ERROR,*)

C#### Subroutine: DRGRAD
C###  Description:
C###    DRGRAD draws gradient vector (eg heat flux).

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'grad00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISGRAD(NEM,NGRSEGM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),NHE(NEM,NXM),
     '  NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NRE(NEM),NRLIST(0:NRM),
     '  NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NVJE(NNM,NBFM,NJM,NEM),
     '  NW(NEM,3,NXM),NXLIST(0:NXM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM,NXM),CG(NMM,NGM),CGE(NMM,NGM,NEM,NXM),
     '  CP(NMM,NPM,NXM),
     '  CURVCORRECT(2,2,NNM,NEM),PG(NSM,NUM,NGM,NBM),RG(NGM),
     '  SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),XP(NKM,NVM,NJM,NPM),
     '  YP(NYM,NIYM,NXM),ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZG(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,INDEX,INDEX_POLYLINE,iw,IWK(6),
     '  N3CO,ne,nolist,no_nrlist,nr,noiw,NTIW,nx,nxc
      REAL*8 RFROMC
      LOGICAL ALL_REGIONS,CBBREV

      CALL ENTERS('DRGRAD',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM draw gradient
C###  Parameter:    <elements (all/#s/GROUP)[all]>
C###    Define what elements are to be included.
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the workstation (GX window) to draw the
C###    gradient on.
C###  Parameter:    <scale FACTOR#[1]>
C###    Define scale of gradient lines.
C###  Parameter:    <rgb=RGB[black]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Parameter:    <regions (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Draws gradient vectors (i.e., heat flux).

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<elements (all/#s/GROUP)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<on (all/WSS)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<scale FACTOR[1]>'
        OP_STRING(5)=BLANK(1:15)//'<rgb=RGB[black]>'
        OP_STRING(6)=BLANK(1:15)//'<regions (#s/all)[1]>'
        OP_STRING(7)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRGRAD',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '    ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        IF(CBBREV(CO,'SCALE',1,noco+1,NTCO,N3CO)) THEN
          SCALE=RFROMC(CO(N3CO+1))
        ELSE
          SCALE=1.0D0
        ENDIF

        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
        ELSE
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLACK')
        ENDIF

        IF(ADD) THEN
          NTGRAD=NTGRAD+1
        ELSE IF(NTGRAD.EQ.0) THEN
          NTGRAD=1
        ENDIF
        CALL ASSERT(NTGRAD.LE.NRM,'>>NRM too small',ERROR,*9999)

        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '      NKH(1,1,1,nr),NPNODE,
     '      nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,1,nx),
     '      ZA,ZP,ERROR,*9999)
          DO nolist=1,NELIST(0)
            ne=NELIST(nolist)
            IF(nr.EQ.NRE(ne)) THEN
              DO noiw=1,NTIW
                IW=IWK(noiw)
                CALL ACWK(iw,1,ERROR,*9999)
                CALL SGGRAD(INDEX,ISEG,ISGRAD(ne,NTGRAD),iw,NBH(1,1,ne),
     '            NBJ(1,ne),ne,NHE(ne,nx),NKHE(1,1,1,ne),NKJE(1,1,1,ne),
     '            NPF(1,1),NPNE(1,1,ne),nr,
     '            NVHE(1,1,1,ne),NVJE(1,1,1,ne),NW(ne,1,nx),nx,
     '            CE(1,ne,nx),CG,CGE,CP(1,1,nx),CSEG,
     '            CURVCORRECT(1,1,1,ne),PG,RG,SE(1,1,ne),XA,XE,XG,XP,
     '            ZA(1,1,1,ne),ZE,ZG,ZP,ERROR,*9999)
                CALL DAWK(iw,1,ERROR,*9999)
              ENDDO !noiw
            ENDIF !nr
          ENDDO !nolist (ne)
        ENDDO !no_nrlist (nr)
      ENDIF

      CALL EXITS('DRGRAD')
      RETURN
 9999 CALL ERRORS('DRGRAD',ERROR)
      CALL EXITS('DRGRAD')
      RETURN 1
      END


