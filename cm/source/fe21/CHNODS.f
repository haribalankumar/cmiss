      SUBROUTINE CHNODS(IBT,IDO,INP,ISEG,ISELNO,ISFIBR,ISFIEL,
     '  ISLINE,ISLINO,ISNONO,MXI,NAN,NBH,NBJ,NEELEM,NEL,NELIST,
     '  NENP,NGAP,NHE,NHP,NKH,NKHE,NKJ,NKJE,NLATNE,
     '  NLL,NLLINE,NLLIST,
     '  NNL,NPF,NPL,NPLIST,NPNE,NPNODE,NQNE,NQNLAT,
     '  NQS,NQXI,NRE,NRLIST,
     '  NVHE,NVHP,NVJE,NVJL,NVJP,NW,NWP,NXLIST,NYNE,NYNP,CURVCORRECT,DL,
     '  PAOPTI,SE,XA,XE,XG,XP,XQ,YG,YP,YQ,ZA,ZC,ZD,ZE,ZP,CSEG,STRING,
     '  FIX,ERROR,*)

C#### Subroutine: CHNODS
C###  Description:
C###    CHNODS allows user to redefine position of nodes.
C###    np is node number of deleted and  then recreated node.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'map000.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISEG(*),ISELNO(NWM,NEM),ISFIBR(NWM,NEM,NGRSEGM),ISFIEL(NWM,NEM),
     '  ISLINE(NWM,2*NGRSEGM),ISLINO(NWM),ISNONO(NWM,NPM),MXI(2,NEM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NEL(0:NELM,NLM),NELIST(0:NEM),
     '  NENP(NPM,0:NEPM,0:NRM),NGAP(NIM,NBM),NHE(NEM,NXM),
     '  NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),NLATNE(NEQM+1),
     '  NLL(12,NEM),NLLINE(0:NL_R_M,0:NRM),NLLIST(0:NLM),
     '  NNL(0:4,12,NBFM),NPF(9,NFM),NPL(5,0:3,NLM),NPLIST(0:NPM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),NQNE(NEQM,NQEM),
     '  NQNLAT(NEQM*NQEM),NQS(NEQM),NQXI(0:NIM,NQSCM),NRE(NEM),
     '  NRLIST(0:NRM),NW(NEM,3,NXM),NWP(NPM,2),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJL(4,NJM,NLM),NVJP(NJM,NPM),
     '  NXLIST(0:NXM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),DL(3,NLM),PAOPTI(*),
     '  SE(NSM,NBFM,NEM),XA(NAM,NJM),
     '  XE(NSM,NJM),XG(NJM,NUM),XP(NKM,NVM,NJM,NPM),
     '  XQ(NJM,NQM),YG(NIYGM,NGM,NEM),
     '  YP(NYM,NIYM,NXM),YQ(NYQM,NIQM,NAM,NXM),
     '  ZA(NAM,NHM,NCM,NEM),ZC(NJM,NEM),ZD(NJM,NDM),
     '  ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM)

!     Local Variables
      INTEGER DIRECTION,I,IBEG,ICHAR,IEND,IFROMC,IIW,
     '  INFO,INSTAT,INT_TEMP(1),
     '  IPICK,IPLANE,
     '  ISEGM,iw,IWK(6),i1,i2,J,LD,MIN_ELEM,NITB,
     '  NOQUES,N3CO,nb,ne,ne2,nh,ni,nj,njj1,njj2,nk,nl,nn,nn2,noelem,
     '  noiw,nolist,
     '  nolist2,no_nelist,nonode,nonr,np,NP1,nr,
     '  nr_list,nrr,NTITLV(1),NTIW,NTLV,NTRL,
     '  nu,num_deriv,nv,nx,nxc,MAX_VERSIONS,pos_inlist
      REAL*8 AMOUNT,ANGLE,AXIS(3),CEN(3),DATD,DATW,DISTANCE,
     '  ELEM_VECT(3,3),MESHD,MESHW,MIN_DIST,MIN_XI(3),NEW_XI,
     '  PXI,RFROMC,RL(1),RL1(3),RL2(3),RL3(3),SUM,TOTAL(3),
     '  TRANS(3,4),TRANS2(3,4),X(11,20),XD(3),XI(3),XX(3),YY(3,NVM),Z(3)
      CHARACTER CHAR*80,FILE*100
      LOGICAL ABBREV,ALL_REGIONS,BOTH,CBBREV,CENTRE,DEFORM,FILEIP,
     '  FOUND,FRAME,INLIST,MOUSE,ORTHOG,PICK1,
     &  PROJECT,PROMPT,REFL,ROTATE,SCALE,TRANSL,UNQUAL
      LOGICAL XFORM             ! Kumar
      LOGICAL TMATRIX
      REAL*8 T(12)
      
      CALL ENTERS('CHNODS',*9999)

      nx=1 !temporary cpb 22/11/94
    
      NOQUES=0
      FILEIP=.FALSE.
      ICHAR=9999
      XFORM=.FALSE.                                                ! Kumar
C     OR> 08/05/06 Initialise variables
      TMATRIX=.FALSE.
      TRANSL=.FALSE.
      SCALE=.FALSE.
      ROTATE=.FALSE.
      REFL=.FALSE.
      PROJECT=.FALSE.
      PICK1=.FALSE.
      ORTHOG=.FALSE.
      CENTRE=.FALSE.
      DIRECTION=0
      NEW_XI=0.0d0
      

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM change nodes;p
C###  Description:
C###    Change nodal coordinates using prompts.

        OP_STRING(1)=STRING(1:IEND)//';p'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM change nodes;m
C###  Parameter:      <deform>
C###    Used to pick the node and then to relocate it. No node number has
C###    to specified
C###  Parameter:      <node (NODE#/pick)[pick]>
C###    Used to pick the node and then to relocate it.
C###    No node number has to specified
C###  Parameter:      <on WS#[1]>
C###    Specify the workstation (GX window) to draw the
C###    points on.

C###  Description:
C###    Change the coordinates of a specified or 'picked' node with the
C###    mouse.  Use the left mouse button to 'pick' a node and then to
C###    relocate it.  Use the central mouse button to quit when you
C###    have finished changing nodes.

        OP_STRING(1)=STRING(1:IEND)//';m'
        OP_STRING(2)=BLANK(1:15)//'<deform>'
        OP_STRING(3)=BLANK(1:15)//'<node (NODE#/pick)[pick]>'
        OP_STRING(4)=BLANK(1:15)//'<on WS#[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM change nodes deformed
C###  Parameter:      <node (#s/all)[all]>
C###    Specifies the node number(s) to change to the deformed coordinates
C###  Parameter:      <region (#s/all)[1]>
C###    Specifies the region number to change to the deformed coordinates
C###  Description:
C###    Change nodal coordinates to deformed coordinates.
C**** CS 25/6/1997

        OP_STRING(1)=STRING(1:IEND)//'deformed'
        OP_STRING(2)=BLANK(1:15)//'<node (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM change nodes translate by DX#,DY#,DZ#
C###  Parameter:      <node (#s/all)[all]>
C###  Parameter:      <region (#s/all)[1]>
C###  Description:
C###    Translate the specified nodal coordinates by the increments
C###    DX,DY,DZ.

        OP_STRING(1)=STRING(1:IEND)//' translate by DX#,DY#,DZ#'
        OP_STRING(2)=BLANK(1:15)//'<node (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM change nodes scale by SX#,SY#,SZ#
C###  Parameter:      <node (#s/all)[all]>
C###  Parameter:      <region (#s/all)[1]>
C###  Description:
C###    Scale the specified nodal coordinates by the ratio SX,SY,SZ.

        OP_STRING(1)=STRING(1:IEND)//' scale by SX#,SY#,SZ#'
        OP_STRING(2)=BLANK(1:15)//'<node (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM change nodes scale to data
C###  Parameter:    <region (#s/all)[1]>
C###  Description:
C###    Scales the mesh in x and y directions so the mesh size is
C###    closer to the data set size.

        OP_STRING(1)=STRING(1:IEND)//' scale to data'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C GMH 31/10/95 Adding orthogonal - use for generating 3d volume from
C surface.

C#### Command: FEM change nodes orthogonal
C###  Parameter:      <node (#s/all)[all]>
C###    Specify the nodes to be moved
C###  Parameter:      <distance R#[0.0]>
C###    Specify the distance the nodes are to be moved
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the regions number to update
C###  Description:
C###    Moves the nodes by the specified distance in a direction
C###    orthogonal to the derivatives specified at the nodes.

        OP_STRING(1)=STRING(1:IEND)//' orthogonal'
        OP_STRING(2)=BLANK(1:15)//'<node (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<distance R#[0.0]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM change nodes reflect in PLANE[1]
C###  Parameter:      <node (#s/all)[all]>
C###    Specify the nodes to reflect
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the regions to apply the reflection to
C###  Description:
C###    Reflect the specified nodes in the specified plane.

        OP_STRING(1)=STRING(1:IEND)//' reflect in PLANE[1]'
        OP_STRING(2)=BLANK(1:15)//'<node (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM change nodes rotate by ANGLE
C###  Parameter:      <node (#s/all)[all]>
C###    Specify the the nodes to move
C###  Parameter:      <about AX#[0.0],AY#[0.0],AZ#[0.0]>
C###    Specify ta point though which the axis of rotation passes
C###  Parameter:      <axis BX#[0.0],BY#[0.0],BZ#[1.0]>
C###    Specify the direction vector of the axis of rotation
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the regions to rotate
C###  Description:
C###    Rotate coordinates of specified nodes by ANGLE about axis
C###    through point (AX,AY,AZ) with direction vector (BX,BY,BZ).

        OP_STRING(1)=STRING(1:IEND)//' rotate by ANGLE'
        OP_STRING(2)=BLANK(1:15)//'<node (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<about AX#[0.0],AY#[0.0],AZ#[0.0]>'
        OP_STRING(4)=BLANK(1:15)//'<axis BX#[0.0],BY#[0.0],BZ#[1.0]>'
        OP_STRING(5)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM change nodes rotate from optimiser
C###  Description:
C###    Rotate nodes in mesh by angle given by an optimisation routine

        OP_STRING(1)=STRING(1:IEND)//' rotate from optimiser'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C GMH 31/10/95 Adding projection - ie move the node to the
C       nearest point on a specified element, and make all
C       derivatives lie in the plane of that element

C#### Command: FEM change nodes project
C###  Parameter:      <node (#s/all)[all]>
C####   Specify the nodes to move
C###  Parameter:      <element (#s/all)[all]>
C###    Specify the elements to project the nodes onto.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the regions to apply change to
C###  Description:
C###    Moves the node to the nearest point on the specified elements,
C###    constrain all derivatives to lie in the plane of the element.

        OP_STRING(1)=STRING(1:IEND)//' project'
        OP_STRING(2)=BLANK(1:15)//'<node (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM change nodes xform FILENAME
C###  Description:
C###

        OP_STRING(1)=STRING(1:IEND)//' xform FILENAME'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C     OR> 08/05/06 Adding command that the transformation matrix also
C     can be provided from the command line and not only from a file.
        
C#### Command: FEM change nodes tmatrix <3x4 Matrix entries>
C###  Parameter:      <node (#s/all)[all]>
C###    Specify nodes that are modified by a rigid body translation
C###  Parameter:      <region (#s/all)[1]>
C###  Description: 
C###        Applys to the specified nodes in the specified region
C###        the 3D transformation matrix. The entries for the 
C###        transformation matrix are in the following orde:
C###        r11, r12, r13, t1, r21, r22, r23, t2, r31, r32, r33, t3
C###        where r stands for the entries of the rotation matrix
C###        and t for the translation entries
        
        OP_STRING(1)=STRING(1:IEND)//' tmatrix'
        OP_STRING(2)=BLANK(1:15)//'<node (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<matrix (12 comma-separated '
     &          //'numbers) >'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------


C#### Command: FEM change nodes hanging
C###  Parameter:      <node (#s/search)[search]>
C###    Specify the nodes to shift or search for nodes that are
C###    extremely close to elements that they are not connected to
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the regions to change.
C###  Description:
C###    Moves the hanging node to the edge of the element it hangs in

C CS  Added 14/5/98

        OP_STRING(1)=STRING(1:IEND)//' hanging'
        OP_STRING(2)=BLANK(1:15)//'<node (#s/search)[search]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','CHNODS',ERROR,*9999)
      ELSE
        CALL CHECKQ(' PM',noco,1,CO,COQU,STRING,*1)
        PROMPT=.FALSE.
        MOUSE =.FALSE.
        UNQUAL=.TRUE.

        IF(ABBREV(COQU(noco,1),'P',1)) THEN
          PROMPT=.TRUE.
          UNQUAL=.FALSE.
        ELSE IF(ABBREV(COQU(noco,1),'M',1)) THEN
          MOUSE =.TRUE.
          UNQUAL=.FALSE.
!news   AAY 20 Jun transformation read from file
        ELSE IF(CBBREV(CO,'XFORM',1,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          FILE=CO(N3CO+1)(IBEG:IEND)
          CALL OPENF(IFILE,'DISK',FILE(1:IEND-IBEG+1)//'.TRN','OLD',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
          CALL GETRAN(IFILE,TRANS)  !this is for nodal values
          CLOSE(UNIT=IFILE)                                              
          CALL OPENF(IFILE,'DISK',FILE(1:IEND-IBEG+1)//'.TRN','OLD',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
          CALL GETRAN(IFILE,TRANS2)
          CLOSE(UNIT=IFILE)                                        
          DO I=1,3
            TRANS2(I,4)=0.d0  !this is for the derivatives
          ENDDO
C ******************************************************************* Kumar 
          XFORM=.TRUE.
C     DO I=1,3
C     PRINT*, TRANS(I,1),TRANS(I,2),TRANS(I,3),TRANS(I,4)
C     ENDDO
C     DO I=1,3
C     PRINT*, TRANS2(I,1),TRANS2(I,2),TRANS2(I,3),TRANS2(I,4)
C     ENDDO
C     PRINT*, XFORM
C ******************************************************************* Kumar
!newe
        ENDIF


        IF(MOUSE) THEN
          IF(CBBREV(CO,'DEFORM',1,noco+1,NTCO,N3CO)) THEN
            DEFORM=.TRUE.
          ELSE
            DEFORM=.FALSE.
          ENDIF
          IF(CBBREV(CO,'ON',1,noco+1,NTCO,N3CO)) THEN
            iw=IFROMC(CO(N3CO+1))
          ELSE
            iw=1
          ENDIF
          IF(CBBREV(CO,'NODE',1,noco+1,NTCO,N3CO)) THEN
            np=IFROMC(CO(N3CO+1))
            PICK1=.FALSE.
          ELSE
            PICK1=.TRUE.
          ENDIF
        ENDIF

        IF(UNQUAL.AND.NTCO.EQ.noco) THEN
          CO(noco+1)='?'
          GO TO 1
        ENDIF

        IF(PROMPT) THEN
 100      FORMAT='($,'' >>Enter node,coordinate & derivative'//
     '      ' numbers [exit]: '',3I4)'
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      3,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,NPT(1),
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
c         CALL IINOUT(1,IOIP,0,FORMAT,3,IDATA,IZERO,0,NPT(1),INFO,
c    '      ERROR,*9999)
          IF(IDATA(1).GT.0) THEN
            np=IDATA(1)
            nj=IDATA(2)
            nk=IDATA(3)+1
C GMH 8/1/97 Update cmgui link
            CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
            FORMAT='($,'' >>Enter new value: '')'
            WRITE(OP_STRING,FORMAT)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            READ(IOIP,'((A))') CHAR
            CALL PARSTR(CHAR,1,NTLV,NTITLV,1,RL,ERROR,*9999)
            XP(nk,1,nj,np)=RL(1)
            GO TO 100
          ENDIF

        ELSE IF(MOUSE) THEN
          IF(PICK1) THEN
            NP1=NPNODE(1,1)
            CALL ASSERT(ISNONO(iw,NP1).GT.0,
     '        '>>Node segments not defined',ERROR,*9999)
            IIW=iw
            FRAME=.FALSE.
            IF(iw.EQ.5.OR.iw.EQ.6) THEN
              FRAME=.TRUE.
              IF(NJT.EQ.2)THEN
                IIW=1
              ELSE
                IIW=3
              ENDIF
            ENDIF
            DO nrr=1,NRT
              DO nonode=1,NPNODE(0,nrr)
                np=NPNODE(nonode,nrr)
                CALL DETECT(iw,ISEG,ISNONO(IIW,np),'DETECTABLE',
     '            ERROR,*9999)
              ENDDO
            ENDDO
            INSTAT=1
            IF(NJT.EQ.2.OR.iw.EQ.4) THEN
              BOTH=.FALSE.
            ELSE
              BOTH=.TRUE.
            ENDIF
            DO WHILE(INSTAT.EQ.1)
              WRITE(OP_STRING,'('' >>Pick node on '',I1,'':'')') IIW
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              CALL ACWK(IIW,0,ERROR,*9999)
              CALL PICK(IIW,'REQUEST',INSTAT,ISEGM,IPICK,
     '          ERROR,*9999)
              CALL DAWK(IIW,0,ERROR,*9999)
              IF(INSTAT.EQ.1) THEN
               np=IFROMC(CSEG(ISEGM)(53:57))
                IF(DEFORM) THEN
                  CALL CHNODE(ISEG,ISELNO,ISLINE,ISLINO,ISNONO,
     '              iw,MXI,NBJ,NEELEM,NHP,NKH,NLL,NLLIST,np,
     '              NPL,NPNE,NPNODE,nx,NYNP,DL,XP,ZC,ZP,CSEG,BOTH,
     '              DEFORM,FIX,FRAME,ERROR,*9999)
                ELSE IF(.NOT.DEFORM) THEN
                  CALL CHNODE(ISEG,ISELNO,ISLINE,ISLINO,ISNONO,
     '              iw,MXI,NBJ,NEELEM,NHP,NKH,NLL,NLLIST,np,
     '              NPL,NPNE,NPNODE,nx,NYNP,DL,XP,ZC,ZP,CSEG,BOTH,
     '              DEFORM,FIX,FRAME,ERROR,*9999)
                ENDIF
                CALL DETECT(iw,ISEG,ISNONO(IIW,np),'DETECTABLE',
     '            ERROR,*9999)
              ENDIF
            ENDDO
            DO nrr=1,NRT
              DO nonode=1,NPNODE(0,nrr)
                np=NPNODE(nonode,nrr)
                CALL DETECT(iw,ISEG,ISNONO(IIW,np),'UNDETECTABLE',
     '            ERROR,*9999)
              ENDDO
            ENDDO
          ELSE IF(.NOT.PICK1) THEN
            BOTH=.FALSE.
            IF(iw.EQ.5.OR.iw.EQ.6) THEN
              FRAME=.TRUE.
              BOTH=.TRUE.
            ELSE
              FRAME=.FALSE.
            ENDIF
            IF(DEFORM) THEN
              CALL CHNODE(ISEG,ISELNO,ISLINE,ISLINO,ISNONO,iw,MXI,
     '          NBJ,NEELEM,NHP,NKH,NLL,NLLIST,
     '          np,NPL,NPNE,NPNODE,nx,NYNP,
     '          DL,XP,ZC,ZP,CSEG,BOTH,DEFORM,FIX,FRAME,ERROR,*9999)
            ELSE IF(.NOT.DEFORM) THEN
              CALL CHNODE(ISEG,ISELNO,ISLINE,ISLINO,ISNONO,iw,MXI,
     '          NBJ,NEELEM,NHP,NKH,NLL,NLLIST,
     '          np,NPL,NPNE,NPNODE,nx,NYNP,
     '          DL,XP,ZC,ZP,CSEG,BOTH,DEFORM,FIX,FRAME,ERROR,*9999)
            ENDIF
          ENDIF

        ELSE IF(UNQUAL.AND.ABBREV(CO(noco+1),'DEFORMED',1).
     &         OR.ABBREV(CO(noco+1),'TRANSLATE'  ,2).
     &         OR.ABBREV(CO(noco+1),'SCALE'      ,1).
     &         OR.ABBREV(CO(noco+1),'ORTHOGONAL' ,2).
     &         OR.ABBREV(CO(noco+1),'REFLECT'    ,2).
     &         OR.ABBREV(CO(noco+1),'ROTATE'     ,2).
     &         OR.ABBREV(CO(noco+1),'PROJECT'    ,2).
     &         OR.ABBREV(CO(noco+1),'HANGING'    ,2).
     &         OR.ABBREV(CO(noco+1),'CENTRE'     ,2).
     &         OR.ABBREV(CO(noco+1),'XFORM'      ,2). ! Kumar
     &         OR.ABBREV(CO(noco+1),'TMATRIX'    ,2)) THEN
          
          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,
     '         *9999)
          CALL PARSE_NODES(NPNODE,NPLIST,noco,NRLIST,NTCO,CO,
     '         ERROR,*9999)
     
C     ************************************************************ Kumar
C     CALL RESET(TRANS)
C     CALL RESET(TRANS2)
C     DEFORM=.FALSE.
C     TRANSL=.FALSE.
C     SCALE =.FALSE.
C     ROTATE=.FALSE.
C     REFL =.FALSE.
C     PROJECT =.FALSE.
C     ORTHOG =.FALSE.
          IF(.NOT.XFORM) THEN
            CALL RESET(TRANS)
            CALL RESET(TRANS2)
            DEFORM=.FALSE.
            TRANSL=.FALSE.
            SCALE =.FALSE.
            ROTATE=.FALSE.
            REFL =.FALSE.
            PROJECT =.FALSE.
            ORTHOG =.FALSE.
          ENDIF  
! ***************************************************************** Kumar     

          IF(ABBREV(CO(noco+1),'DEFORMED',1)) THEN
            DEFORM=.TRUE.
            CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
            nxc=NXLIST(1)
            CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '        ERROR,*9999)
            DO nr_list=1,NRLIST(0)
              nr=NRLIST(nr_list)
              CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '          NKH(1,1,1,nr),NPNODE,
     '          nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,
     '          ERROR,*9999)
              DO nolist=1,NPLIST(0)
                np=NPLIST(nolist)
                DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
                  nj=NJ_LOC(NJL_GEOM,njj2,nr)
                  nh=NH_LOC(njj2,nx)
                  DO nv=1,NVJP(nj,np)
C GMH 8/1/97 Update cmgui link
                    CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                    XP(1,nv,nj,np)=XP(1,nv,nj,np)+ZP(1,nv,nh,np,1)
                  ENDDO  !nv
                ENDDO  !nj
              ENDDO  !nolist
            ENDDO !nrlist
          ELSE IF(ABBREV(CO(noco+1),'TRANSLATE',2)) THEN
            TRANSL=.TRUE.
            IF(CBBREV(CO,'BY',1,noco+1,NTCO,N3CO)) THEN
              CALL PARSRL(CO(N3CO+1),3,NTRL,RL1,ERROR,*9999)
            ENDIF
C MPN 1Apr2003 checking the number of increments is valid
            CALL ASSERT(NTRL.EQ.NJT,'>> # increments entered must '
     '        //'equal # global coordinates',ERROR,*9999)
            AMOUNT=0.d0
            DO nj=1,NJT
              AMOUNT=AMOUNT+RL1(nj)**2
            ENDDO !nj
C MPN 1Apr2003 initialise remaining RL1 terms for use in SHIFT
            DO nj=NJT+1,3
              RL1(nj)=0.d0
            ENDDO !nj
            AMOUNT=DSQRT(AMOUNT)
            CALL SHIFT(RL1,AMOUNT,TRANS,ERROR,*9999)
          ELSE IF(ABBREV(CO(noco+1),'REFLECT',2)) THEN
            REFL=.TRUE.
            IF(CBBREV(CO,'IN',2,noco+1,NTCO,N3CO)) THEN
C LKC 29-APR-98 Need to pass in array
              CALL PARSIL(CO(N3CO+1),1,NTRL,INT_TEMP,ERROR,*9999)
              IPLANE=INT_TEMP(1)
            ELSE
              IPLANE=1
            ENDIF
            TRANS(IPLANE,IPLANE)=-1.d0
            TRANS2(IPLANE,IPLANE)=-1.d0
          ELSE IF(ABBREV(CO(noco+1),'SCALE',1)) THEN
            SCALE=.TRUE.
            IF(CBBREV(CO,'BY',1,noco+1,NTCO,N3CO)) THEN
              CALL PARSRL(CO(N3CO+1),3,NTRL,RL1,ERROR,*9999)
            ENDIF
            IF(CBBREV(CO,'TO',2,noco+1,NTCO,N3CO)) THEN
              CALL MESHXY(IBT,IDO,INP,NBJ,1,NEELEM,NPNE,MESHD,SE,
     '               XP,MESHW,0.0d0,ERROR,*9999)
              CALL DATXY(2,DATD,DATW,0.0d0,ZD,ERROR,*9999)
              RL1(1)=DATW/MESHW
              RL1(2)=DATD/MESHD
              RL1(3)=1.d0
            ENDIF
            CALL SCALE1(RL1(1),TRANS)
            CALL SCALE2(RL1(2),TRANS)
            IF(NJT.EQ.3) THEN
              CALL SCALE3(RL1(3),TRANS)
            ENDIF
          ELSE IF(ABBREV(CO(noco+1),'ROTATE',2)) THEN
            ROTATE=.TRUE.
            IF(CBBREV(CO,'BY',1,noco+1,NTCO,N3CO)) THEN
              ANGLE=RFROMC(CO(N3CO+1))*PI/180.d0
            ELSE IF(CBBREV(CO,'FROM',2,noco+1,NTCO,N3CO)) THEN
              ANGLE=PAOPTI(1)*PI/180.d0
            ELSE
              ANGLE=0.d0
            ENDIF
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' Angle='',D12.4)') ANGLE
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
            IF(CBBREV(CO,'ABOUT',2,noco+1,NTCO,N3CO)) THEN
              CALL PARSRL(CO(N3CO+1),3,NTRL,RL2,ERROR,*9999)
            ELSE
              RL2(1)=0.d0
              RL2(2)=0.d0
              RL2(3)=0.d0
            ENDIF
            IF(CBBREV(CO,'AXIS',2,noco+1,NTCO,N3CO)) THEN
              CALL PARSRL(CO(N3CO+1),3,NTRL,RL3,ERROR,*9999)
            ELSE
              RL3(1)=0.d0
              RL3(2)=0.d0
              RL3(3)=1.d0
            ENDIF
            AMOUNT=DSQRT(RL2(1)**2+RL2(2)**2+RL2(3)**2)
            IF(AMOUNT.GT.1.0D-6) THEN
              CALL SHIFT(RL2,-AMOUNT,TRANS,ERROR,*9999)
            ENDIF
C cpb 4/4/98 Don't see why this is here
C            CALL ZZ(RL3,AXIS,TRANS)
            AXIS(1)=RL3(1)
            AXIS(2)=RL3(2)
            AXIS(3)=RL3(3)
            CALL TWIST(AXIS,ANGLE,TRANS,ERROR,*9999)
            CALL TWIST(AXIS,ANGLE,TRANS2,ERROR,*9999)
            IF(AMOUNT.GT.1.0D-6) THEN
              CALL SHIFT(RL2,AMOUNT,TRANS,ERROR,*9999)
            ENDIF
          ELSE IF(ABBREV(CO(noco+1),'PROJECT',2)) THEN
            PROJECT=.TRUE.
            CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '        ERROR,*9999)
          ELSE IF(ABBREV(CO(noco+1),'ORTHOGONAL',2)) THEN
            ORTHOG=.TRUE.
            IF(CBBREV(CO,'DISTANCE',1,noco+1,NTCO,N3CO)) THEN
              DISTANCE=RFROMC(CO(N3CO+1))
            ELSE
              DISTANCE=0.d0
            ENDIF
          ELSE IF(ABBREV(CO(noco+1),'HANGING',2)) THEN
C****       CS new 14/5/98
C           If the global node np is in any element other
C           than the ones it is conected to it is 'hanging'
            nv=1
            DO nonr=1,NRLIST(0)
              nr=NRLIST(nonr)
              DO nonode=1,NPLIST(0)
                np=NPLIST(nonode)
                NWP(np,1)=0
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  nb=NBJ(1,ne)
                  CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '              NPNE(1,1,ne),nr,NVJE(1,1,1,ne),SE(1,1,ne),
     '              XA,XE,XP,ERROR,*9999)
                  NITB=NIT(NBJ(1,ne))
                  XD(1)=XP(1,1,1,np)
                  XD(2)=XP(1,1,2,np)
                  XD(3)=XP(1,1,3,np)
                  DO ni=1,NITB
                    XI(ni)=0.5d0
                  ENDDO
                  CALL DEXI_POINT(IBT,IDO,INP,LD,NBJ,ne,NITB,nr,
     '              0.d0,XE,XI,XI,XD,.FALSE.,ERROR,*9999)
                  DO no_nelist=1,NENP(np,0,nr)
                    NELIST(no_nelist)=NENP(np,no_nelist,nr)
                  ENDDO
                  IF((LD.NE.0).AND..NOT.
     '              (INLIST(ne,NELIST(1),NENP(np,0,nr),
     '              pos_inlist))) THEN
                    NWP(np,1)=ne
C                   Move the hanging node
                    DISTANCE=0.5d0
                    DO ni=1,NIT(NBJ(1,ne))
                      IF(XI(ni).GT.0.5d0) THEN
                        IF((1.0d0-XI(ni)).LT.DISTANCE) THEN
                          DISTANCE=1.0d0-XI(ni)
                          DIRECTION=ni
                          NEW_XI=1.0d0
                        ENDIF
                      ELSE
                        IF(XI(ni).LT.DISTANCE) THEN
                          DISTANCE=XI(ni)
                          DIRECTION=ni
                          NEW_XI=0.0d0
                        ENDIF
                      ENDIF
                    ENDDO ! ni
                    XI(DIRECTION)=NEW_XI

                    CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                    nv=1
                    DO njj1=1,NJ_LOC(NJL_GEOM,0,nr)
                      nj=NJ_LOC(NJL_GEOM,njj1,nr)
                      DO nu=1,NUT(nb)
                        X(nu,nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                    INP(1,1,nb),nb,nu,XI,XE(1,nj))
                      ENDDO !nu
                    ENDDO !njj1

                    DO njj1=1,NJ_LOC(NJL_GEOM,0,nr)
                      nj=NJ_LOC(NJL_GEOM,njj1,nr)
                      XP(1,nv,nj,np)=X(1,nj)
                      IF((NBI(nb).GE.5).AND.(NBI(nb).LE.7)) THEN
                        DO nk=2,NKT(0,nb)
                          IF(nk.EQ.2) THEN
                            XP(2,nv,nj,np)=X(2,nj)
                          ELSE IF(nk.EQ.3) THEN
                            XP(3,nv,nj,np)=X(4,nj)
                          ELSE IF(nk.EQ.4) THEN
                            XP(4,nv,nj,np)=X(6,nj)
                          ELSE IF(nk.EQ.5) THEN
                            XP(5,nv,nj,np)=X(7,nj)
                          ELSE IF(nk.EQ.6) THEN
                            XP(6,nv,nj,np)=X(9,nj)
                          ELSE IF(nk.EQ.7) THEN
                            XP(7,nv,nj,np)=X(10,nj)
                          ELSE IF(nk.EQ.8) THEN
                            XP(8,nv,nj,np)=X(11,nj)
                          ENDIF
                        ENDDO !nk
                      ELSE
                        CALL ASSERT(.FALSE.,'>>Not implemented for '
     '                    //'basis type',ERROR,*9999)
                      ENDIF
                    ENDDO !njj1
                  ENDIF ! LD
                ENDDO !noelem
              ENDDO !nonode
            ENDDO !nonr
          ELSEIF(CBBREV(CO,'CENTRE',3,noco+1,NTCO,N3CO)) THEN
            CENTRE=.TRUE.
            DO nj=1,3
              CEN(nj)=0.d0
              TOTAL(nj)=0.d0
              AXIS(nj)=0.d0
            ENDDO
            IF(CBBREV(CO,'AT',1,noco+1,NTCO,N3CO)) THEN
              CALL PARSRL(CO(N3CO+1),3,NTRL,CEN,ERROR,*9999)
            ENDIF
            DO nolist=1,NPLIST(0)
              np=NPLIST(nolist)
              DO nj=1,NJT
                TOTAL(nj)=XP(1,1,nj,np)+TOTAL(nj)
              ENDDO
            ENDDO
              AMOUNT=0.d0
            DO nj=1,NJT
              TOTAL(nj)=TOTAL(nj)/NPLIST(0)
              AXIS(nj)=CEN(nj)-TOTAL(nj)
              AMOUNT=AMOUNT+AXIS(nj)**2
            ENDDO
            AMOUNT=DSQRT(AMOUNT)
            CALL RESET(TRANS)
            CALL SHIFT(AXIS,AMOUNT,TRANS,ERROR,*9999)
          ELSEIF(CBBREV(CO,'TMATRIX',2,noco+1,NTCO,N3CO)) THEN
            TMATRIX=.TRUE.
            DO I=1,12 
              T(I)=0.0d0 
            ENDDO
            IF(CBBREV(CO,'MATRIX',1,noco+1,NTCO,N3CO)) THEN
              CALL PARSRL(CO(N3CO+1),12,NTRL,T,ERROR,*9999)
              IF(NTRL.NE.12) THEN
                WRITE(OP_STRING,'('' >>WARNING: Less then 12 entries ' 
     &               //'for the transformation matrix are specified! ' 
     '               //'The missing numbers get initialised by 0.'')')
                CALL WRITES(IOER,OP_STRING,ERROR,*9999)
              ENDIF
              DO I=1,3
                DO J=1,3
                  TRANS(I,J)=T(4*(I-1)+J)
                  TRANS2(I,J)=T(4*(I-1)+J)
                ENDDO
                TRANS(I,J)=T(4*(I-1)+4)
              ENDDO 
            ELSE
              CALL RESET(TRANS)
              CALL RESET(TRANS2)
              WRITE(OP_STRING,'('' >>WARNING: No transformation '
     &             // 'matrix specified, reset transformation '
     &             // 'matrix to identity. '')')
              CALL WRITES(IOER,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF

! ***************************************************************** Kumar
!          IF(TRANSL.OR.SCALE.OR.ROTATE.OR.REFL.OR.CENTRE) THEN
          IF(TRANSL.OR.SCALE.OR.ROTATE.OR.REFL.OR.CENTRE.OR.XFORM.OR
     &         .TMATRIX) THEN

! ***************************************************************** Kumar 
            IF(DOP) THEN                                  
C     KAT 14May01: Can't branch out of critical section.
C     Critical section is not essential.
CC$           call mp_setlock()
              DO I=1,3
                WRITE(OP_STRING,'('' TRANS('',I1,'',j):'',4E12.3)')
     &               I,(TRANS(I,J),J=1,4)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDDO
CC$           call mp_unsetlock()
            ENDIF
            DO nv=1,NVM
              YY(1,nv)=0.d0
              YY(2,nv)=0.d0
              YY(3,nv)=0.d0
            ENDDO
C AJP 3/12/98 Find the maximum number of versions
            MAX_VERSIONS=1
            DO nolist=1,NPLIST(0)
              np=NPLIST(nolist)
              DO nj=1,NJT
                IF(NVJP(nj,np).GT.MAX_VERSIONS)MAX_VERSIONS=NVJP(nj,np)
              ENDDO !nj
            ENDDO !np

            DO nolist=1,NPLIST(0)
              np=NPLIST(nolist)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
C             transform nodal values

              DO nj=1,NJT-1
                IF(NVJP(nj,np).NE.NVJP(nj+1,np)) THEN
                  WRITE(OP_STRING,'(''Warning: version numbers '
     '              //'inconsistent for transforms'')')
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
              ENDDO

C AJP 3/12/98
              CALL ASSERT(MAX_VERSIONS.LE.NVM,
     '          'Increase dim. of YY in chnods',ERROR,*9999)
C old             CALL ASSERT(NVJP(1,np).LE.12,
C old     '          'Increase dim. of YY in chnods',ERROR,*9999)
C old              DO nv=1,NVJP(1,np)
              DO nv=1,MAX_VERSIONS ! ajp 3/12/98
                DO nj=1,NJT
                  YY(nj,nv)=XP(1,nv,nj,np)
                ENDDO
                DO nj=1,NJT
                  IF(nv.LE.NVJP(nj,np)) THEN !an nv version of XP exists
                    YY(nj,nv)=XP(1,nv,nj,np)
                  ELSE !an nv version of XP does not exist..
                    YY(nj,nv)=XP(1,1,nj,np) !so use 1st version
                  ENDIF
                ENDDO !nj
                Z(1)=YY(1,nv)
                Z(2)=YY(2,nv)
                Z(3)=YY(3,nv)
                DO i2=1,3
                  YY(i2,nv)=0.0D0
                  DO i1=1,3
                    YY(i2,nv)=YY(i2,nv)+Z(i1)*TRANS(i2,i1)
                  ENDDO
                  YY(i2,nv)=YY(i2,nv)+TRANS(i2,4)
                ENDDO
C AJPs 3/12/98 Move outside loop
C                DO nj=1,NJT
C                  XP(1,nv,nj,np)=YY(nj,nv)
C                ENDDO
              ENDDO !nv
              DO nv=1,MAX_VERSIONS
                DO nj=1,NJT
                  XP(1,nv,nj,np)=YY(nj,nv)
                ENDDO
              ENDDO !nv
C AJPe 3/12/98

              IF(NKJ(1,np).GT.1) THEN
C               assume basis type is bicubic H and transform derivs
                DO nk=2,NKJ(1,np)
C                  DO nv=1,NVJP(1,np) !ajp 3/12/98
                  DO nv=1,MAX_VERSIONS
                    DO nj=1,NJT
                      IF(nv.LE.NVJP(nj,np)) THEN !an nv version of XP(nk) exists
                        YY(nj,nv)=XP(nk,nv,nj,np)
                      ELSE !an nv version of XP does not exist for this nk..
                        YY(nj,nv)=XP(nk,1,nj,np) !so use 1st version
                      ENDIF
                    ENDDO !nj
                    Z(1)=YY(1,nv)
                    Z(2)=YY(2,nv)
                    Z(3)=YY(3,nv)
                    DO i2=1,3
                      YY(i2,nv)=0.0D0
                      DO i1=1,3
                        YY(i2,nv)=YY(i2,nv)+Z(i1)*TRANS2(i2,i1)
                      ENDDO
                      YY(i2,nv)=YY(i2,nv)+TRANS2(i2,4)
                    ENDDO !i2
C AJPs 3/12/98 Move outside nv loop
C                    DO nj=1,NJT
C                      XP(nk,nv,nj,np)=YY(nj,nv)
C                    ENDDO
                  ENDDO !nv
                  DO nv=1,MAX_VERSIONS
                    DO nj=1,NJT
                      XP(nk,nv,nj,np)=YY(nj,nv)
                    ENDDO
                  ENDDO !nv
C AJPe
                ENDDO !nk
              ENDIF !nkj(1,np)>1
            ENDDO !np

            IF(SCALE) THEN
              DO nb=1,NBFT
C cpb 5/12/96 readding arc length scaling
c cpb swapping over nbi = 4/5
C                IF(NBI(nb).EQ.3.OR.NBI(nb).EQ.4) THEN
                IF(NBI(nb).EQ.3.OR.(NBI(nb).GE.5.AND.NBI(nb).LE.7)) THEN
                  DO nl=1,NLT
                    CALL ARCSCA(IDO,0,0,0,NBJ,NEL(0,nl),nl,NPL(1,0,nl),
     '                NPNE,NVJL(1,1,nl),DL,1.0D-6,XP,ERROR,*9999)
                  ENDDO
                ENDIF
              ENDDO
            ENDIF
          ELSE IF(ORTHOG) THEN
            DO nolist=1,NPLIST(0)
              np=NPLIST(nolist)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
C             find current location
              DO nj=1,NJT
                XX(nj)=XP(1,1,nj,np)
              ENDDO
C               find an element using this node
              CALL ASSERT(NENP(np,0,0).GT.0,
     '          '>>Must have element defined at node',ERROR,*9999)
              ne=NENP(np,1,0)
              nb=NBJ(1,ne)
C               find the nn number
              nn=0
              DO nn2=1,NNT(nb)
                IF (NPNE(nn2,nb,ne).EQ.np) nn=nn2
              ENDDO
              CALL ASSERT(nn.GT.0,
     '          '>>Could not find global node',ERROR,*9999)
C             find a orthogonal vector
              num_deriv=0
              DO nk=1,NKJ(1,np)
C KAT 26Feb01: It seems that the inverse of NKJE should be here.
C                nu=IDO(NKE(nk,nn,nb,ne),nn,0,nb)
                nu=IDO(nk,nn,0,nb)
                IF((nu.EQ.2).OR.(nu.EQ.4).OR.(nu.EQ.7)) THEN
                  num_deriv=num_deriv+1
                  CALL ASSERT(num_deriv.LE.3,
     '              '>>Too many derivatives at node',ERROR,*9999)
                  DO nj=1,NJT
                    ELEM_VECT(nj,num_deriv)=XP(nk,1,nj,np)
                  ENDDO
                ENDIF
              ENDDO !nk
              IF(num_deriv.EQ.0) THEN
C               No vectors - zero orthogonal
                DO nj=1,NJT
                  ELEM_VECT(nj,3)=0.0d0
                ENDDO
              ELSE IF(num_deriv.EQ.2) THEN
C               Get the cross productes
                CALL CROSS(ELEM_VECT(1,1),ELEM_VECT(1,2),
     '            ELEM_VECT(1,3)) !result in elem_vect(nj,3)
                SUM=0.0D0
                DO nj=1,NJT
                  SUM=SUM+ELEM_VECT(nj,3)**2
                ENDDO
                IF(SUM.GT.1.0D-6) THEN
                  SUM=DSQRT(SUM)
                  DO nj=1,NJT
                    ELEM_VECT(nj,3)=ELEM_VECT(nj,3)/SUM
                  ENDDO
                ELSE
                  WRITE(OP_STRING,'('' Zero normal np='',I5,'
     '              //''' nk='',I2)') np,nk
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  DO nj=1,NJT
                    IF(nj.EQ.1) THEN
                      ELEM_VECT(nj,3)=1.0D0
                    ELSE
                      ELEM_VECT(nj,3)=0.0D0
                    ENDIF
                  ENDDO
                ENDIF !zero normal
              ELSE
                ERROR='>>Invalid number of derivatives'
                GO TO 9999
              ENDIF
C             Move the node in the normal direction
              DO nj=1,NJT
                XP(1,1,nj,np)=XP(1,1,nj,np)+DISTANCE*ELEM_VECT(nj,3)
              ENDDO
            ENDDO !nodes
          ELSE IF(PROJECT) THEN
C GMH 31/10/95 Find the closest point on the given elements to the
C   specified node.
            DO nolist=1,NPLIST(0)
              np=NPLIST(nolist)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
C             find current location
              DO nj=1,NJT
                XX(nj)=XP(1,1,nj,np)
              ENDDO
C             loop over elements finding closest projection
              MIN_DIST=0.0D0
              MIN_ELEM=0
              DO nolist2=1,NELIST(0)
                ne=NELIST(nolist2)
                NITB=NIT(NBJ(1,ne))
                CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),
     '            NPF(1,1),NPNE(1,1,ne),
     '            NRE(ne),NVJE(1,1,1,ne),
     '            SE(1,1,ne),XA,XE,XP,ERROR,*9999)
                DO ni=1,NITB !Start in the center of the element
                  XI(ni)=0.5D0
                ENDDO
                FOUND=.FALSE.
                CALL PROJ_ORTHOG(IBT,IDO,INP,NBJ(1,ne),DISTANCE,
     '            XE,XI,XX,FOUND,ERROR,*9999)
                IF(FOUND) THEN
                  IF(MIN_ELEM.EQ.0.OR.DISTANCE.LT.MIN_DIST) THEN
                    DO ni=1,NITB !save xi values
                      MIN_XI(ni)=XI(ni)
                    ENDDO
                    MIN_ELEM=ne
                    MIN_DIST=DISTANCE
                  ENDIF
                ENDIF
              ENDDO !ne
              IF(MIN_ELEM.NE.0) THEN
C               update nodal coord and derivatives
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                  WRITE(OP_STRING,'('' Proj for node '',I5,'
     '              //''' is element '',I5,'
     '              //''' Distance is '',E12.3)') np,MIN_ELEM,MIN_DIST
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                ENDIF
                ne=MIN_ELEM
                NITB=NIT(NBJ(1,ne))
                CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),
     '            NPF(1,1),NPNE(1,1,ne),
     '            NRE(ne),NVJE(1,1,1,ne),
     '            SE(1,1,ne),XA,XE,XP,ERROR,*9999)
C               calculate the xi directions at this point
                DO nj=1,NJT
                  nb=NBJ(nj,ne)
                  DO ni=1,NITB
                    IF(ni.EQ.1) nu=2
                    IF(ni.EQ.2) nu=4
                    IF(ni.EQ.3) nu=7
                    ELEM_VECT(nj,ni)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                INP(1,1,nb),nb,nu,MIN_XI,XE(1,nj))
                  ENDDO
                ENDDO
C                 find an element using this node
                CALL ASSERT(NENP(np,0,0).GT.0,
     '            '>>Must have element defined at node',ERROR,*9999)
                ne2=NENP(np,1,0)
                nb=NBJ(1,ne2)
C                 find the nn number
                nn=0
                DO nn2=1,NNT(nb)
                  IF (NPNE(nn2,nb,ne2).EQ.np) nn=nn2
                ENDDO
                CALL ASSERT(nn.GT.0,
     '            '>>Could not find global node',ERROR,*9999)
C               update nodal values
                DO nk=1,NKJ(1,np)
C KAT 26Feb01: It seems that the inverse of NKJE should be here.
C                  nu=IDO(NKE(nk,nn,nb,ne2),nn,0,nb)
                  nu=IDO(nk,nn,0,nb)
                  IF(nu.EQ.1) THEN !nodal value
                    DO nj=1,NJT
                      nb=NBJ(nj,ne)
                      XP(1,1,nj,np)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                  INP(1,1,nb),nb,1,MIN_XI,XE(1,nj)) !nu=1
                    ENDDO

                  ELSE IF((nu.EQ.2).OR.(nu.EQ.4).OR.(nu.EQ.7)) THEN
C                   How we treat the vector depends on the dimension of
C                   the element
                    IF(NITB.EQ.1) THEN
C                     the derivative is just s1
                      DO nj=1,NJT
                        XX(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                    INP(1,1,nb),nb,2,MIN_XI,XE(1,nj)) !nu=2
                      ENDDO
                    ELSE IF(NITB.EQ.2) THEN
                      DO nj=1,NJT
                        XX(nj)=XP(nk,1,nj,np)
                      ENDDO
C                     if 2D, then the vectors can stay as they are, else
C                     take the normal (s1 cross s2), subtract component
C                     in the normal direction, what we are left with is
C                     in the plane of the element.
                      IF(NJT.EQ.3) THEN
                        CALL CROSS(ELEM_VECT(1,1),ELEM_VECT(1,2),
     '                    ELEM_VECT(1,3)) !result in elem_vect(nj,3)
C                       normalise the normal
                        SUM=0.0D0
                        DO nj=1,NJT
                          SUM=SUM+ELEM_VECT(nj,3)**2
                        ENDDO
                        IF(SUM.GT.1.0D-6) THEN
                          SUM=DSQRT(SUM)
                          DO nj=1,NJT
                            ELEM_VECT(nj,3)=ELEM_VECT(nj,3)/SUM
                          ENDDO
                        ELSE
                          WRITE(OP_STRING,'('' Zero normal np='',I5,'
     '                      //''' nk='',I2)') np,nk
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                          DO nj=1,NJT
                            IF(nj.EQ.1) THEN
                              ELEM_VECT(nj,3)=1.0D0
                            ELSE
                              ELEM_VECT(nj,3)=0.0D0
                            ENDIF
                          ENDDO
                        ENDIF
C                       dot product of normal and vector
                        SUM=0.0D0
                        DO nj=1,NJT
                          SUM=SUM+XX(nj)*ELEM_VECT(nj,3)
                        ENDDO
                        DO nj=1,NJT
                          XX(nj)=XX(nj)-SUM*ELEM_VECT(nj,3)
                        ENDDO
                      ENDIF
                    ELSE IF(NITB.EQ.3) THEN
                      ERROR='>> Not implemented'
                      GOTO 9999
                    ENDIF
C                   scale the vector to unit length
                    SUM=0.0D0
                    DO nj=1,NJT
                      SUM=SUM+XX(nj)**2
                    ENDDO
                    IF(SUM.GT.1.0D-6) THEN
                      SUM=DSQRT(SUM)
                      DO nj=1,NJT
                        XP(nk,1,nj,np)=XX(nj)/SUM
                      ENDDO
                    ELSE
                      WRITE(OP_STRING,'('' Zero vector np='',I5,'
     '                  //''' nk='',I2)') np,nk
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      DO nj=1,NJT
                        IF(nj.EQ.1) THEN
                          XP(nk,1,nj,np)=1.0D0
                        ELSE
                          XP(nk,1,nj,np)=0.0D0
                        ENDIF
                      ENDDO
                    ENDIF

                  ELSE !double derivative
C                   for the moment, set the cross to zero
                    DO nj=1,NJT
                      XP(nk,1,nj,np)=0.0D0
                    ENDDO

                  ENDIF
                ENDDO !nk
              ELSE
                WRITE(OP_STRING,'('' No min proj for node '',I5)') np
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO !nodes
          ENDIF !project
        ENDIF !unqualified

        IF(MOUSE) THEN
          IF(NET(1).GT.0) THEN
            CALL LINCAL(IBT,IDO,INP,0,NBJ,NEELEM,NEL,NENP,
     '        NKJE,NLL,NLLINE,NNL,NPL,NPNE,NPNODE,NRE,NVJE,NVJL,
     '        DL,SE,XP,ERROR,*9999)
          ENDIF
C         DO nb=1,NBFT
C           IF(NBI(nb).EQ.3.OR.NBI(nb).EQ.5) THEN
C             DO nl=1,NLT
C               CALL ARCSCA(IDO,0,0,0,NBJ,NEL(0,nl),nl
C    '            NPL(1,0,nl),NPNE,DL,1.0D-6,XP,ERROR,*9999)
C             ENDDO
C           ENDIF
C         ENDDO
          NTIW=2*NJT-3+IMAP
          DO noiw=1,NTIW
            IWK(noiw)=noiw
          ENDDO
          CALL UPVIEW(IBT,IDO,INP,ISEG,ISELNO,ISFIBR,ISFIEL,ISLINE,
     '      ISLINO,ISNONO,IWK,MXI,NAN,NBH,NBJ,NEELEM,NGAP,NHE(1,nx),
     '      NHP(1,0,nx),NKH,NKHE,NKJE,NLATNE,NLL,NLLIST,NPF,NPL,NPNE,
     '      NPNODE,NQNE,NQNLAT,
     '      NQS,NQXI,NRE,NTIW,NVHE,NVHP,NVJE,NW(1,1,nx),
     '      nx,NYNE,NYNP,CURVCORRECT,DL,SE,XA,XE,XG,XP,XQ,
     '      YG,YP(1,1,nx),YQ,ZA,ZE,
     '      ZP,CSEG,FIX(1,1,nx),ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('CHNODS')
      RETURN
 9999 CALL ERRORS('CHNODS',ERROR)
      CALL EXITS('CHNODS')
      RETURN 1
      END


