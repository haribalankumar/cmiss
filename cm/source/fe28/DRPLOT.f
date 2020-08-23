      SUBROUTINE DRPLOT(ISEG,ISPLOT,NBJ,NBH,NEELEM,NELIST,
     '  NHE,NHP,NKH,NKHE,NKJE,NPF,NPL,NPNE,NPNODE,NRE,
     '  NVHE,NVHP,NVJE,NW,NYNE,NYNP,CURVCORRECT,
     '  DL,SE,XA,XE,XP,YP,ZA,ZE,ZP,CSEG,STRING,ERROR,*)

C#### Subroutine: DRPLOT
C###  Description:
C###    DRPLOT draws 3D Phigs plot on workstation 9.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'graf00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISPLOT(NHM,0:NEM,NGRSEGM),NBH(NHM,NCM,NEM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),
     '  NPL(5,0:3,NLM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NRE(NEM),NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),DL(3,NLM),SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM),ZA(NAM,NHM,NCM,NEM),
     '  ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,INDEX,INDEX_OLD,INDEX_POLYLINE,
     '  iw,N3CO,nb,nc,ne,nh,nj,nl,nodx,noelem,nolist,nr,ns,NTDX,nx
      REAL*8 PXL(2,20),X_AXIS(3,2),XL(20,3)
      LOGICAL CBBREV

      CALL ENTERS('DRPLOT',*9999)
      nx=1 ! temporary cpb 22/11/94
      nc=1 ! Temporary MPN 12-Nov-94
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM draw plot
C###  Parameter:    <of VARIABLE#[1]>
C###    Specifies the variable number
C###  Parameter:    <in (all/ELEMENT#s)[all]>
C###    Specifies the element numbers
C###  Parameter:    <xi_3 VALUE#[0.0]{>=0.0,<=1.0}>
C###    Specifies the xi3 values
C###  Parameter:    <rgb=RGB[black]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Description:
C###    draws 3D Phigs plot

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<of VARIABLE#[1]>'
        OP_STRING(3)=BLANK(1:15)//'<in (all/ELEMENT#s)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<Xi_3 VALUE#[0.0]{>=0.0,<=1.0}>'
        OP_STRING(5)=BLANK(1:15)//'<rgb=RGB[black]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRPLOT',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        IW=9
        NTPLOT=NTPLOT+1
        CALL ASSERT(NTPLOT.LE.NRM,'>>NRM too small',ERROR,*9999)
        IF(CBBREV(CO,'OF',1,noco+1,NTCO,N3CO)) THEN
          nh=IFROMC(CO(N3CO+1))
        ELSE
          nh=1
        ENDIF
        IF(CBBREV(CO,'IN',1,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NE_R_M,NELIST(0),NELIST(1),ERROR,*9999)
        ELSE
          NELIST(0)=0
          DO nr=1,NRT
            DO noelem=NELIST(0)+1,NELIST(0)+NEELEM(0,nr)
              NELIST(noelem)=NEELEM(noelem,nr)
            ENDDO
            NELIST(0)=NELIST(0)+NEELEM(0,nr)
          ENDDO
        ENDIF
        IF(CBBREV(CO,'XI_3',1,noco+1,NTCO,N3CO)) THEN
C GMH 2/9/95 Unused          XIPLOT=RFROMC(CO(N3CO+1))
        ELSE
C GMH 2/9/95 Unused          XIPLOT=0.0D0
        ENDIF
        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
        ELSE
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLACK')
        ENDIF

        nr=1 !Needs fixing
        IF(ITYP1(nr,nx).GT.1.AND.ITYP2(nr,nx).NE.11) THEN
!new MPN 6-Jan-95: current soln now stored in YP(ny,1) for nonlin probs
          CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '      NKH(1,1,1,nr),NPNODE,
     '      nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)
!old
c          IF(ITYP1(nr,nx).LE.4) THEN
c            CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),NKH(1,1,1,nr),NPNODE,
c     '        nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)
c          ELSE IF(ITYP1(nr,nx).EQ.5) THEN
c            CALL YPZP(4,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),NKH(1,1,1,nr),NPNODE,
c     '        nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)
c          ENDIF
        ENDIF

        CALL ACWK(iw,1,ERROR,*9999)
        CALL OPEN_SEGMENT(ISPLOT(nh,0,NTPLOT),ISEG,iw,'PLOT',INDEX,
     '    INDEX_OLD,0,1,CSEG,ERROR,*9999)

C ***   Draw axes and mesh lines at z=0 in parent structure
        DO nj=1,3
          X_AXIS(nj,1)=0.0D0
          X_AXIS(nj,2)=0.0D0
        ENDDO
        X_AXIS(1,2)=0.8D0*DBLE(XMAX)
        CALL POLYLINE(1,iw,2,X_AXIS,ERROR,*9999)
        CALL TEXT(1,iw,'X',X_AXIS(1,2),ERROR,*9999)
        X_AXIS(1,2)=0.0D0
        X_AXIS(2,2)=0.8D0*DBLE(YMAX)
        CALL POLYLINE(1,iw,2,X_AXIS,ERROR,*9999)
        CALL TEXT(1,iw,'Y',X_AXIS(1,2),ERROR,*9999)
        X_AXIS(2,2)=0.0D0
        X_AXIS(3,2)=1.1D0*DBLE(ZMAX)
        CALL POLYLINE(1,iw,2,X_AXIS,ERROR,*9999)
        CALL TEXT(1,iw,'Z',X_AXIS(1,2),ERROR,*9999)

        DO nl=1,NLT
          IF(NJT.EQ.2) THEN
C GMG 2/9/95 Initialise NTDX
            NTDX=0
            CALL XPXL('UNDEFORMED',1,NPL(1,0,nl),nr,NTDX,DL(1,nl),XL,XP,
     '        ERROR,*9999)
            DO nodx=1,NTDX
              PXL(1,nodx)=XL(nodx,1)
              PXL(2,nodx)=XL(nodx,2)
            ENDDO
            CALL POLYLINE(1,iw,NTDX,PXL,ERROR,*9999)
          ENDIF
        ENDDO
        CALL CLOSE_SEGMENT(ISPLOT(nh,0,NTPLOT),iw,ERROR,*9999)
        CALL DAWK(iw,1,ERROR,*9999)

C ***   Draw fill area sets at z=dependent variable nh with one
C ***   structure for each element
        DO nolist=1,NELIST(0)
          CALL ACWK(iw,1,ERROR,*9999)
          ne=NELIST(nolist)
          IF(DOP) THEN
            WRITE(OP_STRING,
     '        '('' Element '',I3,'' variable '',I1)') ne,nh
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '      NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '      SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
          IF(DOP) THEN
C CPB 8/4/94 Adding NJ_LOC
C            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)+JTYP9+JTYP11
           DO nj=1,NJ_LOC(0,0,nr)
              WRITE(OP_STRING,
     '          '('' XE(ns,'',I1,''): '',8E11.3,:,(/8E11.3))')
     '          nj,(XE(ns,nj),ns=1,NST(NBJ(nj,ne)))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
           ENDDO
          ENDIF
          IF(ITYP1(nr,nx).EQ.1.OR.ITYP2(nr,nx).EQ.11) THEN
            IF(DOP) THEN
              WRITE(OP_STRING,'('' Transfer field from XE to ZE'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            nb=NBJ(NJ_LOC(NJL_FIEL,1,nr),ne)
            DO ns=1,NST(nb)
              ZE(ns,nh)=XE(ns,NJ_LOC(NJL_FIEL,1,nr))
            ENDDO
          ELSE
            nr=1 !may need generalizing
            CALL ZPZE(NBH(1,1,ne),nc,NHE(ne,nx),NKHE(1,1,1,ne),
     '        NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1,nx),
     '        nx,CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),
     '        ZE,ZP,ERROR,*9999)
            nb=NBH(nh,1,ne)
          ENDIF
          IF(DOP) THEN
            WRITE(OP_STRING,
     '        '('' ZE(ns,'',I1,''): '',8E11.3,:,(/8E11.3))')
     '        nh,(ZE(ns,nh),ns=1,NST(nb))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CALL SGPLOT(INDEX,ISEG,ISPLOT(nh,ne,NTPLOT),iw,nb,ne,nh,
     '      XE,ZE,CSEG,ERROR,*9999)
          CALL DAWK(iw,1,ERROR,*9999)
        ENDDO
      ENDIF

      CALL EXITS('DRPLOT')
      RETURN
 9999 CALL ERRORS('DRPLOT',ERROR)
      CALL EXITS('DRPLOT')
      RETURN 1
      END


