      SUBROUTINE DRSHEE(IBT,IDO,INP,ISEG,ISSHEE,
     '  NAN,NBH,NBJ,NEELEM,NELIST,NENP,NHE,NHP,NKJE,NKH,NNB,
     '  NPF,NPNE,NPNODE,NRE,NVHP,NVJE,NXI,NYNE,NYNP,
     '  SE,XA,XE,XG,XP,YP,ZA,ZP,CSEG,STRING,ERROR,*)

C#### Subroutine: DRSHEE
C###  Description:
C###    DRSHEE draws element sheet orientations.  Constant vectors of
C###    Xi-coordinate length DXIF are drawn on plane Xi(3)=XIF at
C###    increments of DXI1 and DXI2 in Xi(1) and Xi(2) direction.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISEG(*),ISSHEE(NWM,NEM,NGRSEGM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),NENP(NPM,0:NEPM,0:NRM),
     '  NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),NKJE(NKM,NNM,NJM,NEM),
     '  NKH(NHM,NPM,NCM,0:NRM),NNB(4,4,4,NBFM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NRE(NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM),ZA(NAM,NHM,NCM,NEM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG1,IEND,IEND1,IFROMC,INDEX,INDEX_POLYLINE,iw,
     '  N3CO,ne,noelem,nolist,nr,NTDXI,nx
      REAL*8 DXI(3),DXI2,DXI3,RFROMC,THETA
      LOGICAL CBBREV

      CALL ENTERS('DRSHEE',*9999)
      nx=1 ! temporary cpb 22/11/94
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM draw sheets
C###  Parameter:    <theta=ANGLE#[0.0]{degrees}>
C###  Parameter:    <dxi VALUE#s[0.2,0.13]>
C###    Increments int the Xi(1) and Xi(2) directions.
C###  Parameter:    <region #[1]>
C###    Specify the region number to use.
C###  Parameter:    <in (all/ELEMENT#s[all]>
C###    Specify the element(s) in which to draw the sheets. The "all"
C###    keyword specifies all currently defined elements.
C###  Parameter:    <rgb=RGB[blue]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Description:
C###    Draws element sheet orientations.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<theta=ANGLE#[0.0]{degrees}>'
        OP_STRING(3)=BLANK(1:15)//'<dxi VALUE#s[0.2,0.13]>'
        OP_STRING(4)=BLANK(1:15)//'<region #[1]>'
        OP_STRING(5)=BLANK(1:15)//'<in ELEMENT#s[all]>'
        OP_STRING(6)=BLANK(1:15)//'<rgb=RGB[blue]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRSHEE',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        IF(CBBREV(CO,'REGION',1,noco+1,NTCO,N3CO)) THEN
          nr=IFROMC(CO(N3CO+1))
        ELSE
          nr=1
        ENDIF

        CALL ASSERT(NJ_LOC(njl_fibr,0,nr).EQ.3,
     '    '>>No sheets defined',ERROR,*9999)
        IF(CBBREV(CO,'THETA',2,noco+1,NTCO,N3CO)) THEN
          THETA=PI*RFROMC(CO(N3CO+1))/180.0D0
        ELSE
          THETA=0.0D0
        ENDIF
        IF(CBBREV(CO,'DXI',1,noco+1,NTCO,N3CO)) THEN
          CALL PARSRL(CO(N3CO+1),3,NTDXI,DXI,ERROR,*9999)
          DXI2=DXI(1)
          DXI3=DXI(2)
C          IF(NTDXI.GT.2) THEN
C            DXIS=DXI(3)
C          ELSE
C            DXIS=0.13D0
C          ENDIF
        ELSE
          DXI2=0.20D0
          DXI3=0.20D0
C          DXIS=0.13D0
        ENDIF
        IF(CBBREV(CO,'IN',1,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NE_R_M,NELIST(0),NELIST(1),ERROR,*9999)
        ELSE
          NELIST(0)=NEELEM(0,nr)
          DO noelem=1,NEELEM(0,nr)
            NELIST(noelem)=NEELEM(noelem,nr)
          ENDDO
        ENDIF
        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
        ELSE
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLUE')
        ENDIF

        IF(ADD) THEN
          NTSHEE=NTSHEE+1
        ELSE IF(NTSHEE.EQ.0) THEN
          NTSHEE=1
        ENDIF
        CALL ASSERT(NTSHEE.LE.NRM,'>>NRM too small',ERROR,*9999)
        IF(ITYP2(nr,nx).EQ.14.OR.ITYP2(nr,nx).EQ.15) THEN
          CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '      NKH(1,1,1,nr),NPNODE,
     '      nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)
        ENDIF
        CALL NENXI(NBJ,NEELEM,NENP,NNB,NPNE,NXI,ERROR,*9999)
        IW=13
C CPB 10/9/92 Changed workstation update mode
        CALL ACWK(iw,1,ERROR,*9999)
        DO nolist=1,NELIST(0)
          ne=NELIST(nolist)
C GMH 2/9/95 Unused          nb=NBJ(1,ne)
          IF(NXI(-3,1,ne).EQ.0.AND.
     '      (NET(1).LE.30.OR.(NET(1).gt.30.and.NE.GT.30))) THEN
!                              !restrict to endo temporary!!!!!
            CALL SGSHEE(INDEX,IBT,IDO,INP,ISEG,ISSHEE(iw,ne,NTSHEE),iw,
     '        NAN,NBJ,ne,NEELEM,NKJE,NPF,NPNE,NRE,NVJE,NXI,
     '        DXI2,DXI3,SE,THETA,XA,XE,XG,XP,CSEG,ERROR,*9999)
          ENDIF
        ENDDO
        CALL DAWK(iw,1,ERROR,*9999)
      ENDIF

      CALL EXITS('DRSHEE')
      RETURN
 9999 CALL ERRORS('DRSHEE',ERROR)
      CALL EXITS('DRSHEE')
      RETURN 1
      END


