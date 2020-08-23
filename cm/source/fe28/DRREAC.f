      SUBROUTINE DRREAC(ISEG,ISREAC,NBH,NEELEM,NHE,NHP,NKH,
     '  NPL,NPNODE,NRLIST,NVHP,NXLIST,NYNE,NYNP,NYNR,
     '  XP,YP,ZA,ZP,CSEG,STRING,FIX,ERROR,*)

C#### Subroutine: DRREAC
C###  Description:
C###    DRREAC draws reactions.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'reac00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISREAC(NWM),NBH(NHM,NCM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),
     '  NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),
     '  NPL(5,0:3,NLM),
     '  NPNODE(0:NP_R_M,0:NRM),NRLIST(0:NRM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NXLIST(0:NXM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM)
      REAL*8 XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM),ZA(NAM,NHM,NCM,NEM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,INDEX,INDEX_POLYLINE,iw,IWK(6),iy,
     '  nc,N3CO,noiw,no_nrlist,nr,NTIW,nx,nxc
      REAL*8 RFROMC
      LOGICAL ALL_REGIONS,CBBREV

      CALL ENTERS('DRREAC',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM draw reactions
C###  Parameter:    <iy=IY#[4]> <nc=NC#[1]> {finite elasticity problems}
C###  Parameter:    <iy=IY#[1]> <nc=NC#[2]> {all other problem types}
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the worksation(s). The all keyword specifies all
C###    currently defined workstations.
C###  Parameter:    <scale FACTOR#[1.0]>
C###    Scales the magnitude of the reaction forces for drawing.
C###  Parameter:    <rgb=RGB[black]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Parameter:    <(all/constrained)[all]>
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Draws the reaction forces from the current solution on the
C###    specified workstation, in the specified colour. The scale can be
C###    adjusted to enable viewing.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<iy=IY#[4]> <nc=NC#[1]>'
     '    //'{finite elasticity problems}'
        OP_STRING(3)=BLANK(1:15)//'<iy=IY#[1]> <nc=NC#[2]>'
     '    //'{all other problem types}'
        OP_STRING(4)=BLANK(1:15)//'<on (all/WS#s)[all]>'
        OP_STRING(5)=BLANK(1:15)//'<scale FACTOR#[1.0]>'
        OP_STRING(6)=BLANK(1:15)//'<rgb=RGB[black]>'
        OP_STRING(7)=BLANK(1:15)//'<(all/constrained)[all]>'
        OP_STRING(8)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(9)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRREAC',ERROR,*9999)
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

        IF(CBBREV(CO,'CONSTRAINED',1,noco+1,NTCO,N3CO)) THEN
          CONSTR=.TRUE.
        ELSE IF(CBBREV(CO,'ALL',1,noco+1,NTCO,N3CO)) THEN
          CONSTR=.FALSE.
        ELSE
          CONSTR=.FALSE.
        ENDIF

        IF(CBBREV(CO,'IY',1,noco+1,NTCO,N3CO)) THEN
          iy=IFROMC(CO(N3CO+1))
        ELSE
          IF(ITYP1(NRLIST(1),nx).EQ.5) THEN !finite elasticity problems
            iy=4
          ELSE  !all other problem types
            iy=1
          ENDIF
        ENDIF

        IF(CBBREV(CO,'NC',1,noco+1,NTCO,N3CO)) THEN
          nc=IFROMC(CO(N3CO+1))
        ELSE
          IF(ITYP1(NRLIST(1),nx).EQ.5) THEN !finite elasticity problems
            nc=1
          ELSE  !all other problem types
            nc=2
          ENDIF
        ENDIF

        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
        ELSE
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLACK')
        ENDIF

        IF(CBBREV(CO,'SCALE',1,noco+1,NTCO,N3CO)) THEN
          FACTOR=RFROMC(CO(N3CO+1))
        ELSE
          FACTOR=1.0D0
        ENDIF

        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
C         Put deformed solution into ZP (iy=1 in YPZP)
          CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),NKH(1,1,1,nr),
     '      NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,
     '      YP(1,1,nx),ZA,ZP,ERROR,*9999)
          DO noiw=1,NTIW
            IW=IWK(noiw)
            CALL ACWK(iw,1,ERROR,*9999)
            IF(DOP) THEN
              WRITE(*,'('' nc='',I2,'' nr='',I2,'' nx='',I2)') nc,nr,nx
            ENDIF
            CALL SGREAC(INDEX,ISEG,ISREAC(iw),iw,iy,nc,NHP(1,nr,nx),
     '        NKH(1,1,1,nr),NPL,NPNODE,nr,nx,NYNR(0,0,nc,nr,nx),
     '        CSEG,FIX(1,1,nx),XP,YP(1,1,nx),ZP,ERROR,*9999)
            CALL DAWK(iw,1,ERROR,*9999)
          ENDDO !noiw

        ENDDO !no_nrlist
      ENDIF

      CALL EXITS('DRREAC')
      RETURN
 9999 CALL ERRORS('DRREAC',ERROR)
      CALL EXITS('DRREAC')
      RETURN 1
      END

