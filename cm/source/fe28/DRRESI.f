      SUBROUTINE DRRESI(ISEG,ISRESI,NBH,NEELEM,NHE,NHP,NKH,
     '  NPNODE,NPNY,NVHP,NXLIST,NYNE,NYNO,NYNP,
     '  RESID,YP,ZA,ZP,CSEG,STRING,ERROR,*)

C#### Subroutine: DRRESI
C###  Description:
C###    DRRESI draws residuals.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'reac00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISRESI(NWM),NBH(NHM,NCM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),
     '  NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),
     '  NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM,NXM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NXLIST(0:NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 RESID(*),YP(NYM,NIYM,NXM),ZA(NAM,NHM,NCM,NEM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,INDEX,INDEX_POLYLINE,iw,IWK(6),
     '  N3CO,noiw,nr,NTIW,nx,nxc
      REAL*8 RFROMC
      LOGICAL CBBREV

      CALL ENTERS('DRRESI',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM draw residuals
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the worksation(s). The all keyword specifies all
C###    currently defined workstations.
C###  Parameter:    <scale FACTOR#[1.0]>
C###    Scales the magnitude of the residual vectors for drawing.
C###  Parameter:    <rgb=RGB[black]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Draws the residual vectors on the specified workstation, in
C###    the specified colour. The scale can be adjusted to enable
C###    viewing.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<on (all/WS#s[all]>'
        OP_STRING(3)=BLANK(1:15)//'<scale FACTOR#[1.0]>'
        OP_STRING(4)=BLANK(1:15)//'<rgb=RGB[black]>'
        OP_STRING(5)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRRESI',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_OPTI,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this opti class',
     '    ERROR,*9999)

        IF(CBBREV(CO,'SCALE',1,noco+1,NTCO,N3CO)) THEN
          FACTOR=RFROMC(CO(N3CO+1))
        ELSE
          FACTOR=1.0d0
        ENDIF

        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
        ELSE
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLACK')
        ENDIF

        nr=1 !Needs fixing

        CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '    NKH(1,1,1,nr),NPNODE,
     '    nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)

        DO noiw=1,NTIW
          IW=IWK(noiw)
          CALL ACWK(iw,1,ERROR,*9999)
          CALL SGRESI(INDEX,ISEG,ISRESI(iw),iw,
     '      NPNY(0,1,0,nx),nr,nx,NYNO(0,1,1,nr,nx),RESID,ZP,
     '      CSEG,ERROR,*9999)
          CALL DAWK(iw,1,ERROR,*9999)
        ENDDO
      ENDIF

      CALL EXITS('DRRESI')
      RETURN
 9999 CALL ERRORS('DRRESI',ERROR)
      CALL EXITS('DRRESI')
      RETURN 1
      END


