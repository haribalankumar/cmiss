      SUBROUTINE DRNODS(ISEG,ISNONO,NBH,NEELEM,NHE,NHP,NKH,
     '  NPNODE,NRLIST,NVHP,NXLIST,NYNE,NYNP,
     '  XP,YP,ZA,ZP,CSEG,STRING,FIX,ERROR,*)

C#### Subroutine: DRNODS
C###  Description:
C###    DRNODS draws nodal parameters.
C**** NPNODE(0,nr) is the total number of nodes in region nr.
C**** NPNODE(nonode,nr), nonode=1..NPNODE(0,nr) are the node numbers.
C**** NPT(nr) is the highest node number.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER ISEG(*),ISNONO(NWM,NPM),NBH(NHM,NCM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),
     '  NKH(NHM,NPM,NCM,0:NRM),NPNODE(0:NP_R_M,0:NRM),NRLIST(0:NRM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NXLIST(0:NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
      REAL*8 XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM),
     '  ZA(NAM,NHM,NCM,NEM),ZP(NKM,NVM,NHM,NPM,NCM)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER IBEG,IBEG1,IEND,IEND1,INDEX,INDEX_TEXT,iw,IWK(6),
     '  N3CO,nj,njj,nh,nhx,noiw,nolist,nonode,np,nr,NTIW,nxc,nx
      REAL*8 X(3),Z(3)
      LOGICAL ALL_REGIONS,CBBREV,DEFORM

      CALL ENTERS('DRNODS',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM draw nodes
C###  Parameter:    <(undeformed/deformed)[undeformed]>
C###     Draw undeformed/deformed geometry
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the workstation (GX window) to draw the
C###    nodes on.
C###  Parameter:    <rgb=RGB[black]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Parameter:    <region (#s/all)[all]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Draws the defined nodes on the specified workstation, in the
C###    specified colour.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<(undeformed/deformed)[undeformed]>'
        OP_STRING(3)=BLANK(1:15)//'<on (all/WS#s)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<rgb=RGB[black]>'
        OP_STRING(5)=BLANK(1:15)//'<region (#s/all)[all]>'
        OP_STRING(6)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRNODS',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        ALL_REGIONS=.FALSE.

        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

C new MPN 7Apr97
        IF(CBBREV(CO,'UNDEFORMED',1,noco+1,NTCO,N3CO)) THEN
          DEFORM=.FALSE.
        ELSE IF(CBBREV(CO,'DEFORMED',1,noco+1,NTCO,N3CO)) THEN
          DEFORM=.TRUE.
        ELSE
          DEFORM=.FALSE.
        ENDIF
        IF(DEFORM) THEN
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nxc=NXLIST(1)
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.GT.0,'>>No nx defined for this class',
     '      ERROR,*9999)
        ELSE
          nx=1
        ENDIF

        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
          INDEX=INDEX_TEXT(0,'WIDTH1','FONT1',CO(N3CO+1))
        ELSE
          INDEX=INDEX_TEXT(0,'WIDTH1','FONT1','BLACK')
        ENDIF

        DO nj=1,3
          X(nj)=0.0d0
        ENDDO

C LC 25/2/97 archived section : Is this code used? PJH 7Oct94

C GMH 2/9/95 Unused        FRAME=.FALSE.
C        DO noiw=1,NTIW
C          IW=IWK(noiw)
C GMH 2/9/95 Unused          IF(IW.GE.5) FRAME=.TRUE.
C        ENDDO

        DO noiw=1,NTIW
          IW=IWK(noiw)
          CALL ACWK(iw,1,ERROR,*9999)
          DO nolist=1,NRLIST(0)
            nr=NRLIST(nolist)
            CALL ASSERT(NPT(nr).GT.0,'>>no nodes defined',ERROR,*9999)
            IF(DEFORM) THEN
              CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '          NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,
     '          YP(1,1,nx),ZA,ZP,ERROR,*9999)
            ENDIF
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              IF(DEFORM) THEN
                DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
                  nh=NH_LOC(nhx,nx)
                  X(nhx)=ZP(1,1,nh,np,1)
                ENDDO !nh
                CALL XZ(ITYP11(nr),X,Z)
              ELSE
                DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                  nj=NJ_LOC(NJL_GEOM,njj,nr)
                  X(nj)=XP(1,1,nj,np)
                ENDDO !nj
                CALL XZ(ITYP10(nr),X,Z)
              ENDIF
              CALL SGNODE(INDEX,ISEG,ISNONO(iw,np),iw,NHP(1,nr,nx),
     '          NKH(1,1,1,nr),np,nr,nx,NYNP,FIX(1,1,nx),Z(1),Z(2),Z(3),
     '          CSEG,ERROR,*9999)
            ENDDO !nonode
          ENDDO !nolist
          CALL DAWK(iw,1,ERROR,*9999)
        ENDDO !noiw (iw)
      ENDIF

      CALL EXITS('DRNODS')
      RETURN
 9999 CALL ERRORS('DRNODS',ERROR)
      CALL EXITS('DRNODS')
      RETURN 1
      END


