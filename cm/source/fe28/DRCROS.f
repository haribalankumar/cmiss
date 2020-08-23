      SUBROUTINE DRCROS(IBT,IDO,INP,ISCROS,ISEG,NAN,NBJ,NEELEM,
     '  NENP,NHE,NKHE,NNB,NPF,NPNE,NRE,NVJE,NXI,NXLIST,
     '  PG,SE,XA,XE,XG,XP,ZD,ZE,ZG,CSEG,STRING,ERROR,*)

C#### Subroutine: DRCROS
C###  Description:
C###    DRCROS draws cross-section in segment ISCROS(iw,nocros).

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'ntsg00.cmn'
      INCLUDE 'post00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISCROS(NWM,NGRSEGM),ISEG(*),NAN(NIM,NAM,NBFM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NENP(NPM,0:NEPM,0:NRM),
     '  NHE(NEM,NXM),NKHE(NKM,NNM,NHM,NEM),NNB(4,4,4,NBFM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NRE(NEM),NVJE(NNM,NBFM,NJM,NEM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),NXLIST(0:NXM)
      REAL*8 PG(NSM,NUM,NGM,NBM),SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),XP(NKM,NVM,NJM,NPM),
     '  ZD(NJM,NDM),ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,INDEX,INDEX_POLYLINE,iw,j,N3CO,N4CO,nj,nx,nxc
      REAL*8 RFROMC,ZDD_LOCAL(3),ZVAL
      CHARACTER SHEET_TYPE*6
      LOGICAL CBBREV,DATA

      CALL ENTERS('DRCROS',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM draw cross-section
C###  Parameter:    <x/y/z (POSITION#/data)>
C###  Parameter:    <post>
C###  Parameter:    <rgb=RGB[red]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    This command draws a cross-section of the element boundaries at
C###    a specified value of x,y or z in the y-z, x-z, or x-y planes
C###    respectively.  Alternatively, the position of the cross-section
C###    can be defined from a data file containing the coordinates of
C###    the single point at which the cross-sections are desired.

        OP_STRING(1)=STRING(1:IEND)//' x/y/z POSITION#/data'
        OP_STRING(2)=BLANK(1:15)//'<post>'
        OP_STRING(3)=BLANK(1:15)//'<rgb=RGB[red]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM draw cross-section radial/normal
C###  Parameter:    <theta=ANGLE#[0.0]{degrees}>
C###  Parameter:    <rgb=RGB[red]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###

        OP_STRING(1)=STRING(1:IEND)//' radial/normal'
        OP_STRING(2)=BLANK(1:15)//'<theta=ANGLE#[0.0]{degrees}>'
        OP_STRING(3)=BLANK(1:15)//'<rgb=RGB[red]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRCROS',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        IF(CBBREV(CO,'POST',1,noco+1,NTCO,N3CO)) THEN
          POSTSCRIPT=.TRUE.
        ENDIF

        IF(CBBREV(CO,'X',1,noco+1,NTCO,N3CO)) THEN
          nj=1
          IW=2
        ELSE IF(CBBREV(CO,'Y',1,noco+1,NTCO,N3CO)) THEN
          nj=2
          IW=1
        ELSE IF(CBBREV(CO,'Z',1,noco+1,NTCO,N3CO)) THEN
          nj=3
          IW=3
        ELSE IF(CBBREV(CO,'RADIAL',1,noco+1,NTCO,N3CO)) THEN
          SHEET_TYPE='RADIAL'
          IF(ITYP10(1).EQ.2) THEN      !cyl. polar coords
            nj=2
          ELSE IF(ITYP10(1).EQ.4) THEN !prolate coords
            nj=3
          ENDIF
          iw=13
        ELSE IF(CBBREV(CO,'NORMAL',1,noco+1,NTCO,N3CO)) THEN
          SHEET_TYPE='NORMAL'
          IF(ITYP10(1).EQ.2) THEN      !cyl. polar coords
            nj=2
          ELSE IF(ITYP10(1).EQ.4) THEN !prolate coords
            nj=3
          ENDIF
          IW=13
        ELSE
          CO(noco+1)='?'
          GO TO 1
        ENDIF

        IF(CBBREV(CO,'DATA',1,noco+1,NTCO,N4CO)) THEN
          DATA=.TRUE.
          ZVAL=ZD(nj,1)
          DO j=1,NJT
            ZDD_LOCAL(j)=ZD(j,1)
          ENDDO
        ELSE IF(CBBREV(CO,'THETA',1,noco+1,NTCO,N3CO)) THEN
          DATA=.FALSE.
          ZVAL=RFROMC(CO(N3CO+1))
          IF(IW.EQ.13) ZVAL=ZVAL*PI/180.0D0
        ELSE
          ZVAL=RFROMC(CO(N3CO+1))
        ENDIF

        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
        ELSE
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','RED')
        ENDIF

        NTCROS=1
        CALL ASSERT(NTCROS.LE.NRM,'>>NRM too small',ERROR,*9999)

C CPB 28/3/94 Changed workstation update mode
            CALL ACWK(iw,1,ERROR,*9999)
c            CALL ACWK(iw,0,ERROR,*9999)
        CALL SGCROS(INDEX,IBT,IDO,INP,ISCROS(iw,NTCROS),ISEG,iw,NAN,
     '    NBJ,NEELEM,NENP,NHE(1,nx),nj,NKHE,NNB,NTCROS,
     '    NPF,NPNE,NRE,NVJE,NXI,nx,CSEG,DATA,PG,SE,
     '    SHEET_TYPE,XA,XE,XG,XP,ZDD_LOCAL,ZE,ZG,ZVAL,ERROR,*9999)
        CALL DAWK(iw,1,ERROR,*9999)
C        CALL DAWK(iw,0,ERROR,*9999)
      ENDIF

      CALL EXITS('DRCROS')
      RETURN
 9999 CALL ERRORS('DRCROS',ERROR)
      CALL EXITS('DRCROS')
      RETURN 1
      END


