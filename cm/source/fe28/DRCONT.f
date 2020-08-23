      SUBROUTINE DRCONT(IBT,IDO,INP,ISCONO,ISCONT,ISEG,MXI,
     '  NAN,NBJ,NBH,NEELEM,NELIST,NHE,NHP,NKH,NKHE,NKJE,NPF,NPNE,NPNODE,
     '  NRE,NRLIST,NTCOVA,NVHE,NVHP,NVJE,NW,NXLIST,NYNE,NYNP,
     '  COVA,CURVCORRECT,PG,SE,XA,XE,XG,XP,YG,YP,ZA,ZE,ZG,ZP,CSEG,
     '  STRING,ERROR,*)

C#### Subroutine: DRCONT
C###  Description:
C###    DRCONT draws contours and contour numbers in one workstation
C###    viewport.  Contour values are held in COVA(ne,nocova),
C###    nocova=1,NTCOVA(ne) where NTCOVA(ne) is number of contours in
C###    element ne (segment ISCONO(i,ne) identifies number).  Maximum
C###    of 50 contours per element may be drawn.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'moti00.cmn'
      INCLUDE 'ntsg00.cmn'
      INCLUDE 'scal00.cmn'
      INCLUDE 'time00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISCONO(NHM,NEM),ISCONT(NHM,NEM,NGRSEGM),ISEG(*),MXI(2,NEM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),NHE(NEM,NXM),
     '  NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKJE(NKM,NNM,NJM,NEM),NKHE(NKM,NNM,NHM,NEM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NRE(NEM),NRLIST(0:NRM),NTCOVA(NEM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM),NXLIST(0:NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 COVA(NEM,50),CURVCORRECT(2,2,NNM,NEM),PG(NSM,NUM,NGM,NBM),
     '  SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),
     '  XP(NKM,NVM,NJM,NPM),YG(NIYGM,NGM,NEM),YP(NYM,NIYM,NXM),
     '  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZG(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER i,IBEG,IBEG2,IEND,IEND2,IFROMC,INDEX,INDEX_POLYLINE,iw,
     '  IWK(6),N3CO,nb,nc,ne,njj,nh,nocova,noiw,nolist,no_nrlist,norl,
     '  nr,ns,NTIW,NTRL,nx,nxc
      REAL*8 RFROMC,RL(100),TOL,XICONT,XISTEP
      CHARACTER CHAR2*9,CONTYP*10
      LOGICAL ALL_REGIONS,CBBREV,LABEL,STATIC

      CALL ENTERS('DRCONT',*9999)

      nc=1 ! Temporary MPN 12-Nov-94

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        nr=1 !temporary for test below
        nx=1 !temporary for test below
        IF(ITYP2(nr,nx).GT.0) THEN !dependent variables drawn
          CHAR2='dependent'
        ELSE
          CHAR2='field'
        ENDIF
        CALL STRING_TRIM(CHAR2,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM draw contour
C###  Parameter:    <(field/dependent/strain/Cauchy/vector)[field]>
C###    Define what the contour is of.
C###  Parameter:    <value (full/LIST_OF_VALUE#s)[full]>
C###  Parameter:    <of VARIABLE#[1]>
C###  Parameter:    <elements (all/#s/GROUP)[all]>
C###    Define what elements are to be included.
C###  Parameter:    <on WS_ID#[1]>
C###    Specify the window number contour will be drawn on (default #1).
C###  Parameter:    <Xi_3 #[0.0]{>=0.0,<=1.0}>
C###    Distance along Xi_3 direction where contours will be drawn.
C###  Parameter:    <tol VALUE#[0.01]>
C###  Parameter:    <step SIZE#[0.02]>
C###    Step size between contours
C###  Parameter:    <when (static/TIME#)[static]>
C###  Parameter:    <wrt TIME_REF#[0.0]>
C###  Parameter:    <rgb=RGB[black]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Parameter:    <label (on/off)[on]>
C###  Parameter:    <regions (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    This command draws contours of the solution variable specified
C###    on the appropriate workstation.


        OP_STRING( 1)=STRING(1:IEND)
     '    //' <(field/dependent/strain/cauchy/vector)['
     '    //CHAR2(IBEG2:IEND2)//']>'
        OP_STRING( 2)=BLANK(1:15)
     '    //' <value (full/LIST_OF_VALUE#s)[full]>'
        OP_STRING( 3)=BLANK(1:15)//'<of VARIABLE#[1]>'
        OP_STRING( 4)=BLANK(1:15)//'<elements (all/#s/GROUP)[all]>'
        OP_STRING( 5)=BLANK(1:15)//'<on WS_ID#[1]>'
        OP_STRING( 6)=BLANK(1:15)//'<Xi_3 #[0.0]{>=0.0,<=1.0}>'
        OP_STRING( 7)=BLANK(1:15)//'<tol VALUE#[0.01]>'
        OP_STRING( 8)=BLANK(1:15)//'<step SIZE#>[0.02]'
        OP_STRING( 9)=BLANK(1:15)//'<when (static/TIME#)[static]>'
        OP_STRING(10)=BLANK(1:15)//'<wrt TIME_REF#[0.0]>'
        OP_STRING(11)=BLANK(1:15)//'<rgb=RGB[black]>'
        OP_STRING(12)=BLANK(1:15)//'<label (on/off)[on]>'
        OP_STRING(13)=BLANK(1:15)//'<regions (#s/all)[1]>'
        OP_STRING(14)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRCONT',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '    ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

! GBS 28/01/96
        NTCONT=1
!        NTCONT=NTCONT+1
!        CALL ASSERT(NTCONT.LE.NRM,'>>NRM too small',ERROR,*9999)
        IF(CBBREV(CO,'FIELD',1,noco+1,NTCO,N3CO)) THEN
          CONTYP='FIELD'
        ELSE IF(CBBREV(CO,'DEPENDENT',1,noco+1,NTCO,N3CO)) THEN
          CONTYP='DEPENDENT'
        ELSE IF(CBBREV(CO,'STRAIN',1,noco+1,NTCO,N3CO)) THEN
          CONTYP='STRAIN'
        ELSE IF(CBBREV(CO,'CAUCHY',1,noco+1,NTCO,N3CO)) THEN
          CONTYP='CAUCHY'
        ELSE IF(CBBREV(CO,'VECTOR',1,noco+1,NTCO,N3CO)) THEN
          CONTYP='VECTOR'
        ELSE
          nr=NRLIST(1)
          nx=1 !temporary for test below
          IF(ITYP2(nr,nx).GT.0) THEN !dependent variables drawn
            CONTYP='DEPENDENT'
          ELSE
            CONTYP='FIELD'
          ENDIF
        ENDIF

        IF(CONTYP(1:5).NE.'FIELD') THEN
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '      ERROR,*9999)
        ELSE
          nx=1
        ENDIF

        IF(CBBREV(CO,'OF',2,noco+2,NTCO,N3CO)) THEN
          nh=IFROMC(CO(N3CO+1))
        ELSE
          nh=1
        ENDIF
        njj=nh !for fields

        IF(CBBREV(CO,'VALUE',1,noco+1,NTCO,N3CO)) THEN
          CALL PARSRL(CO(N3CO+1),100,NTRL,RL,ERROR,*9999)
        ELSE
          nr=NRLIST(1)
          CALL CM_RANGE(nh,0,NBH,NBJ,NEELEM,NHE(1,nx),NHP(1,0,nx),
     '      NKH,NKHE,NKJE,NPF,NPNE,NPNODE,NRLIST,NVHE,NVHP,NVJE,
     '      NW(1,1,nx),nx,NYNE,NYNP,
     '      CURVCORRECT,PG,SE,CONTYP,XA,XE,XG,XP,YG,YP(1,1,nx),ZA,ZE,
     '      ZG,ZP,ERROR,*9999)
          NTRL=10
          DO norl=1,NTRL
            RL(norl)=ZMINI+DBLE(norl-1)/DBLE(NTRL-1)*(ZMAXI-ZMINI)
          ENDDO
        ENDIF
C        ZVAL_MIN=RL(1)
C        ZVAL_MAX=RL(NTRL)

        IF(CBBREV(CO,'XI_3',1,noco+2,NTCO,N3CO)) THEN
          XICONT=RFROMC(CO(N3CO+1))
        ELSE
          XICONT=0.0D0
        ENDIF

        IF(CBBREV(CO,'TOL',1,noco+2,NTCO,N3CO)) THEN
          TOL=RFROMC(CO(N3CO+1))
        ELSE
          TOL=0.01D0
        ENDIF

        IF(CBBREV(CO,'STEP',3,noco+2,NTCO,N3CO)) THEN
          XISTEP=RFROMC(CO(N3CO+1))
        ELSE
          XISTEP=0.02D0
        ENDIF

        IF(CBBREV(CO,'WHEN',3,noco+2,NTCO,N3CO)) THEN
          TIME=RFROMC(CO(N3CO+1))
          STATIC=.FALSE.
        ELSE
          STATIC=.TRUE.
        ENDIF

        IF(CBBREV(CO,'WRT',3,noco+2,NTCO,N3CO)) THEN
          TIME_REF=RFROMC(CO(N3CO+1))
        ELSE
          TIME_REF=0.0D0
        ENDIF

        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
        ELSE
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLACK')
        ENDIF

        IF(CBBREV(CO,'LABEL',1,noco+1,NTCO,N3CO)) THEN
          IF(CBBREV(CO,'ON',2,N3CO+1,N3CO+1,N3CO)) THEN
            LABEL=.TRUE.
          ELSE IF(CBBREV(CO,'OFF',2,N3CO+1,N3CO+1,N3CO)) THEN
            LABEL=.FALSE.
          ELSE
            LABEL=.TRUE.
          ENDIF
        ELSE
          LABEL=.TRUE.
        ENDIF

        IF(CONTYP(1:5).EQ.'FIELD') THEN
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            CALL ASSERT(NBJ(NJ_LOC(NJL_FIEL,njj,nr),NELIST(1)).GT.0,
     '        '>>Field not defined',ERROR,*9999)
          ENDDO
        ELSE IF(CONTYP(1:9).EQ.'DEPENDENT') THEN
          CALL ASSERT(NBH(nh,1,NELIST(1)).GT.0,
     '      '>>Field not defined',ERROR,*9999)
        ENDIF

C GMH 4/11/96 Use code for getting IW number properly (from DRFIEL)
        DO noiw=1,NTIW
          IW=IWK(noiw)
C       IW=IWK(1)
          IF(DOP) THEN
            WRITE(OP_STRING,'('' NTCONT='',I3,'' CONTYP='',A10,'//
     '        ''' IW='',I3,'' XICONT='',F5.2)') NTCONT,
     '        CONTYP(1:10),iw,XICONT
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,
     '        '('' NELIST='',(40I3))') (NELIST(i),i=0,NELIST(0))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' NTRL  ='',I3,'' RL='',(12F8.2))')
     '        NTRL,(RL(i),i=1,NTRL)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

C CPB 10/9/92 Changed workstation update mode
          CALL ACWK(iw,1,ERROR,*9999)
C          CALL ACWK(iw,0,ERROR,*9999)
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            IF(CONTYP(1:5).NE.'FIELD') THEN
              IF(STATIC) THEN
                CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '            NKH(1,1,1,nr),NPNODE,
     '            nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,
     '            ERROR,*9999)
              ELSE IF(.NOT.STATIC) THEN
                CALL YPZP(MOTION_IY,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '            NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,
     '            YP(1,1,nx),ZA,ZP,ERROR,*9999)
              ENDIF
            ENDIF
            DO nolist=1,NELIST(0)
              ne=NELIST(nolist)
              IF(nr.EQ.NRE(ne)) THEN
                IF(DOP) THEN
                  WRITE(OP_STRING,
     '              '('' Element '',I3,'' variable '',I1)') ne,nh
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '            NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '            SE(1,1,ne),XA(1,1,ne),XE,XP,
     '            ERROR,*9999)
C             IF(DOP) THEN
C               DO nj=1,NJ_LOC(NJL_GEOM,0,nr)+JTYP9+JTYP11
C                 WRITE(IOOP,'('' XE(ns,'',I1,''): '',8E11.3,:,(/8E11.3))')
C    '              nj,(XE(ns,nj),ns=1,NST(NBJ(nj,ne)))
C               ENDDO
C             ENDIF
                IF(CONTYP(1:5).EQ.'FIELD') THEN
                  IF(DOP) THEN
                    WRITE(OP_STRING,
     '                '('' Transfer field from XE to ZE'')')
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                  nb=NBJ(NJ_LOC(NJL_FIEL,njj,nr),ne)
                  DO ns=1,NST(nb)
                    ZE(ns,1)=XE(ns,NJ_LOC(NJL_FIEL,njj,nr))
                  ENDDO
                ELSE
                  CALL ZPZE(NBH(1,1,ne),nc,NHE(ne,nx),NKHE(1,1,1,ne),
     '              NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),
     '              NW(ne,1,nx),nx,CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '              ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
                  nb=NBH(nh,1,ne)
                ENDIF
C             IF(DOP) THEN
C              WRITE(IOOP,'('' ZE(ns,'',I1,''): '',8E11.3,:,(/8E11.3))')
C    '            nh,(ZE(ns,nh),ns=1,NST(nb))
C             ENDIF
                NTCOVA(ne)=NTRL
                DO nocova=1,NTCOVA(ne)
                  COVA(ne,nocova)=RL(nocova)
                ENDDO
                IF(CONTYP(1:5).EQ.'FIELD') THEN !NBH etc not defined
                  CALL SGCONT(INDEX,IBT,IDO,INP,ISCONO(1,ne),
     '              ISCONT(1,ne,NTCONT),ISEG,iw,MXI(1,ne),NAN,
     '              NBJ(NJ_LOC(NJL_FIEL,njj,nr),ne),NBJ(1,ne),ne,1,1,
     '              nr,NTCOVA(ne),nx,CONTYP,COVA,CSEG,LABEL,PG,
     '              STATIC,TOL,
     '              XE,XG,XICONT,XISTEP,ZE,ZG,ERROR,*9999)
                ELSE
                  CALL SGCONT(INDEX,IBT,IDO,INP,ISCONO(nh,ne),
     '              ISCONT(nh,ne,NTCONT),ISEG,iw,MXI(1,ne),NAN,
     '              NBH(1,1,ne),NBJ(1,ne),ne,nh,NHE(ne,nx),
     '              nr,NTCOVA(ne),nx,CONTYP,COVA,CSEG,LABEL,PG,
     '              STATIC,TOL,
     '              XE,XG,XICONT,XISTEP,ZE,ZG,ERROR,*9999)
                ENDIF
              ENDIF
            ENDDO !nolist (ne)
          ENDDO !no_nrlist (nr)
          CALL DAWK(iw,1,ERROR,*9999)
C        CALL DAWK(iw,0,ERROR,*9999)
        ENDDO !noiw
      ENDIF

      CALL EXITS('DRCONT')
      RETURN
 9999 CALL ERRORS('DRCONT',ERROR)
      CALL EXITS('DRCONT')
      RETURN 1
      END


