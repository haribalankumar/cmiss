      SUBROUTINE UPINIT(IBT,IDO,INP,LD_NP,NBH,NBJ,NEELEM,NELIST,NENP,
     &  NENQ,NFFACE,NFLIST,NHE,NHP,NKH,NKHE,NKJ,NKJE,NNF,NONY,NPF,
     &  NPLIST1,NPLIST3,NPNE,NPNODE,NP_INTERFACE,NQNP,NQS,NQXI,NRE,
     &  NRLIST,NVHE,NVHP,NVJE,NW,NXLIST,NXQ,NYNE,NYNP,AQ,CQ,CURVCORRECT,
     &  DF,DNUDXQ,DXDXIQ,DXDXIQ2,SE,XA,XAB,XE,XP,XQ,YP,YQ,ZA,ZCROSSING,
     &  ZD,ZE,ZP,STRING,FIX,FIXQ,ERROR,*)
C      SUBROUTINE UPINIT(NBJ,NEELEM,NELIST,NENQ,NHP,NKH,NKJ,NONY,NPLIST1,
C     &  NPLIST3,NPNE,NPNODE,NP_INTERFACE,NQNP,NQS,NQXI,NRLIST,NVHP,
C     &  NXLIST,NXQ,NYNE,NYNP,AQ,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,XAB,XP,XQ,YP,
C     '  YQ,ZCROSSING,STRING,FIX,FIXQ,ERROR,*)

C#### Subroutine: UPINIT
C###  Description:
C###    UPINIT updates the initial conditions.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'pulm00.cmn'
      INCLUDE 'nqloc00.inc'
      INCLUDE 'tol00.cmn'
      INCLUDE 'trsf00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     &  LD_NP(NDM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     &  NELIST(0:NEM),NENP(NPM,0:NEPM,0:NRM),NENQ(0:8,NQM),
     &  NFFACE(0:NF_R_M,NRM),NFLIST(0:NFM),NHE(NEM,NXM),
     &  NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),NKHE(NKM,NNM,NHM,NEM),
     &  NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),NNF(0:17,6,NBFM),
     &  NONY(0:NOYM,NYM,NRCM,0:NRM,NXM),NPF(9,NFM),NPLIST1(0:NPM),
     &  NPLIST3(0:NPM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     &  NP_INTERFACE(0:NPM,0:3),NQNP(NPM),NQS(NEQM),NQXI(0:NIM,NQSCM),
     &  NRE(NEM),NRLIST(0:NRM),NVHE(NNM,NBFM,NHM,NEM),
     &  NVHP(NHM,NPM,NCM,0:NRM),NVJE(NNM,NBFM,NJM,NEM),NXLIST(0:NXM),
     &  NW(NEM,3,NXM),NXQ(-NIM:NIM,0:4,0:NQM,NAM),
     &  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     &  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 AQ(NMAQM,NQM),CURVCORRECT(2,2,NNM,NEM),CQ(NMM,NQM,NXM),
     &  DF(NFM),DNUDXQ(3,3,NQM),DXDXIQ(3,3,NQM),DXDXIQ2(3,3,NQM),
     &  SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),XAB(NORM,NEM),XE(NSM,NJM),
     &  XP(NKM,NVM,NJM,NPM),XQ(NJM,NQM),YP(NYM,NIYM,NXM),
     &  YQ(NYQM,NIQM,NAM,NXM),ZA(NAM,NHM,NCM,NEM),ZD(NJM,NDM),
     &  ZCROSSING(NY_TRANSFER_M,NTSM),ZE(NSM,NHM),
     &  ZP(NKM,NVM,NHM,NPM,NCM)
C      REAL*8 AQ(NMAQM,NQM),CQ(NMM,NQM,NXM),DNUDXQ(3,3,NQM),
C     &  DXDXIQ(3,3,NQM),DXDXIQ2(3,3,NQM),XAB(NORM,NEM),
C     '  XP(NKM,NVM,NJM,NPM),XQ(NJM,NQM),YP(NYM,NIYM,NXM),
C     '  YQ(NYQM,NIQM,NAM,NXM),ZCROSSING(NY_TRANSFER_M,NTSM)
      CHARACTER ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM),FIXQ(NYQM,NIYFIXM,NXM)
!     Local Variables
      INTEGER depvar,IBEG,IEND,in1,in2,irow,n1,N3CO,nb,nc,nd,ne,nef,nf,
     &  NFNP(0:10,NP_R_M),nh,nhx,niqV,niy,nj,nj_field,njj_field,nk,nn,
     &  nne,noelem,noface,nonode,no_nrlist,np,npp,np2,nq,nrr,
     &  nr_bem,nr_bem2,nts,nv,nr,nx,nxc,nx_bem,nx_ext,nx_upd,
     &  ny,nyd
      REAL*8 A,ALPHA,AREA,B,CALC_TIME_FROM_SAMPLE,constant_value,DIST1,
     &  DIST2,DPHIDS,D2PHIDNDS,DX(3,3),FLUX,ONEMINUSALPHA,PHI,RFROMC,
     &  SCALE_F,SUM_XNORM(3),X(3),X05(3),XI(3),XI05(3),XNORM(3),
     '  grav,maxdemens(6)
      CHARACTER TYPE*9
      LOGICAL ABBREV,ADD_TO,ALL_REGIONS,AT_ELEMENTS,AT_NODES,BC,CBBREV,
     &  DERIV,DISPLACEMENT,ERROR_FLAG,FOUND,INCR,INTERPOLATE,
     &  NORMAL_FORCE,UNDEFORMED
!     Functions
      INTEGER IFROMC
      REAL*8 PXI
      LOGICAL INLIST
      
      CALL ENTERS('UPINIT',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM update initial
C###  Parameter:      <from (field/grid/iterate)[field]>
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <class #s[1]>
C###    Specify the class number (of solve type) to use.
C###  Parameter:      <alpha #[1.0]>
C###    This parameter is used in conjunction with the grid
C###    qualifier. It is used as alpha*grid + (1-alpha)*bem
C###    is the new bem initial condition. Alpha can be any real
C###    number in the [0-1] range.
C###  Description:
C###    If the GRID option is specified there are 3 regions
C###    and 3 classes expected.
C###
C###    The first region is the region which contains the
C###      boundary conditions to be updated.
C###    The second region is the boundary element region which
C###      lies beneath the grid mesh.
C###    The third region is the region which contains the grid
C###      mesh.
C###
C###    The first class is the extracellular potential class for
C###      the grid solution.
C###    The second class is the boundary element class to which
C###      the boundary conditions are updated.
C###    The third class is the extracellular update class.
C###
C###    If FIX(ny,niy,nx_bem) is set for a potential boundary
C###    condition then a potential will be calculated and set.
C###    If FIX(ny,niy,nx_bem) is set for a flux boundary
C###    condition then a flux will be calculated and set.
C###    If cubic hermite elements are used then arc length
C###    derivatives will also be calculated and set.

        OP_STRING(1)=STRING(1:IEND)//' <from (field/grid/iterate)'
     '    //'[field]>'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(4)=BLANK(1:15)//'<alpha #[1.0]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update initial
C###  Parameter:      <from zcrossing>
C###    Specify the use of the zcrossing array.
C###  Parameter:      <interpolate>
C###    Interpolate around zero-crossing to get more accurate estimate.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Updates the YP array using the initial estimates of the
C###    activation time as given by the zcrossing array.

        OP_STRING(1)=STRING(1:IEND)//' <from zcrossing>'
        OP_STRING(2)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update initial
C###  Description:
C###    Update initial conditions defined at nodes or elements (in first
C###    specified region) from either a node-based or data-based field
C###    (in second specified region).
C###  Parameter:      <elements/nodes (#s/all)[elements/all]>
C###    Specify whether nodal or elemental boundary conditions used, and
C###    list the element or node numbers/groups.
C###  Parameter:      <region (#)[1]>
C###    Specify the region number that contains the nodes or elements.        
C###  Parameter:      <field/datafield (#) [1]>
C###    Specify whether the field is nodal (field) or data (datafield)
C###    and supply the field number.        
C###  Parameter:      <displacement/force [displacement]>
C###    Specify whether the update is for displacement or force BCs.        
C###  Parameter:      <region (#)[1]>
C###    Specify the region number containing the field.
C###  Parameter:      <class #s[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//' <elem/nodes (#s/all)[elem/all]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#) [1]>'
        OP_STRING(2)=BLANK(1:15)//'<field/datafield (#) [1]>'
        OP_STRING(2)=BLANK(1:15)//'<displacement/force [displacement]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#) [1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #s[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
C#### Command: FEM update initial flow
C###  Description:
C###    Update boundary flow conditions defined at elements from fields
C###    stored at data points, mapped to the nodes at the ends of
C###    terminal elements in the specified element list or group.
C###  Parameter:      <element (#s/all)[all]>
C###  Parameter:      <region (#)[1]>
C###    Specify the region number that contains the elements.        
C###  Parameter:      <field (#) [1]>
C###    Specify the data field number.        
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region number containing the field.
C###  Parameter:      <class #s[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//' <normal>'
        OP_STRING(2)=BLANK(1:15)//'<force/displacement [displacement]>'
        OP_STRING(3)=BLANK(1:15)//'<faces  #s/all [all]>'
        OP_STRING(4)=BLANK(1:15)//'<region #/all [all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc','UPINIT',ERROR,*9999)
      ELSE
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        
C MHT for update initial nodes/elements from field        
        AT_NODES=.FALSE.
        AT_ELEMENTS=.FALSE.
        DISPLACEMENT=.TRUE.
        BC=.FALSE. !line causing the problem
        ADD_TO=.FALSE.
        INCR=.FALSE.
        IF(CBBREV(CO,'FROM',2,noco+1,NTCO,N3CO)) THEN
          IF(ABBREV(CO(N3CO+1),'FIELD',2)) THEN
            TYPE='FIELD'
            nr=NRLIST(1)
            njj_field=IFROMC(CO(N3CO+2))
            IF(CBBREV(CO,'NODES',4,noco+1,NTCO,N3CO)) THEN
              CALL PARSE_NODES(NPNODE,NPLIST1,noco,NRLIST,NTCO,CO,ERROR,
     &          *9999)
              IF(NPLIST1(0).GT.0)THEN
                AT_NODES=.TRUE.
              ENDIF
            ELSE IF(CBBREV(CO,'ELEMENTS',4,noco+1,NTCO,N3CO)) THEN
              CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     &          ERROR,*9999)
              IF(NELIST(0).GT.0) AT_ELEMENTS=.TRUE.
            ENDIF
            IF(AT_NODES)THEN
              nj_field=NJ_LOC(NJL_FIEL,njj_field,nr)
              CALL ASSERT(nj_field.GT.0,'>>Field not defined',ERROR,
     &          *9999)
            ELSE IF(AT_ELEMENTS)THEN
              nj_field=NEJ_LOC(njj_field,nr)
              CALL ASSERT(nj_field.GT.0,'>>Field not defined',ERROR,
     '          *9999)
            ENDIF

            IF(CBBREV(CO,'FORCE',4,noco+1,NTCO,N3CO))THEN
               DISPLACEMENT=.FALSE.
               IF(CBBREV(CO,'BC',2,noco+1,NTCO,N3CO))THEN
                  BC=.TRUE.
               ENDIF
            ENDIF
            IF(CBBREV(CO,'DEPVAR',3,noco+1,NTCO,N3CO)) THEN
              depvar=IFROMC(CO(N3CO+1))
            ELSE
              depvar=1
            ENDIF
            SCALE_F=1.d0
            IF(CBBREV(CO,'SCALE',3,noco+1,NTCO,N3CO)) SCALE_F=
     &        RFROMC(CO(N3CO+1))
            IF(CBBREV(CO,'ADD_TO',4,noco+1,NTCO,N3CO))
     &        ADD_TO=.TRUE.
            IF(CBBREV(CO,'INCR',4,noco+1,NTCO,N3CO))
     &        INCR=.TRUE.
            
          ELSE IF(ABBREV(CO(N3CO+1),'GRID',2)) THEN
            TYPE='GRID'
          ELSE IF(ABBREV(CO(N3CO+1),'ITERATE',2)) THEN
            TYPE='ITERATE'
          ELSE IF(ABBREV(CO(N3CO+1),'ZCROSSING',2)) THEN
            TYPE='ZCROSSING'
          ELSE
            TYPE='FIELD'
          ENDIF
        ELSE IF(CBBREV(CO,'CONSTANT',4,noco+1,NTCO,N3CO)) THEN
          AT_NODES=.FALSE.
          AT_ELEMENTS=.FALSE.
          TYPE = 'CONSTANT'
          constant_value=RFROMC(CO(N3CO+1))
          noco=3
          IF(CBBREV(CO,'NODES',4,noco+1,NTCO,N3CO)) THEN
            CALL PARSE_NODES(NPNODE,NPLIST1,noco,NRLIST,NTCO,CO,ERROR,
     &        *9999)
            IF(NPLIST1(0).GT.0)THEN
              AT_NODES=.TRUE.
            ENDIF
C.. KSB 01/2008: Below used if want to set the inlet reference node for gravity calculations
            IF(CBBREV(CO,'INLET_REF_NODE',5,noco+1,NTCO,N3CO)) THEN
              np_in=IFROMC(CO(N3CO+1))
            ELSE
              np_in=0 !must define this if want to use it: write error later
            ENDIF
          ELSE IF(CBBREV(CO,'ELEMENTS',4,noco+1,NTCO,N3CO)) THEN
            CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     &        ERROR,*9999)
            IF(NELIST(0).GT.0) AT_ELEMENTS=.TRUE.
          ENDIF
        ELSE
          AT_NODES=.FALSE.
          AT_ELEMENTS=.FALSE.
          TYPE='FIELD'
          IF(CBBREV(CO,'NORMAL',4,noco+1,NTCO,N3CO)) THEN
            TYPE = 'NORMAL'
            IF(CBBREV(CO,'NODES',4,noco+1,NTCO,N3CO)) THEN
              CALL PARSE_NODES(NPNODE,NPLIST1,noco,NRLIST,NTCO,CO,ERROR,
     &          *9999)
              IF(NPLIST1(0).GT.0)THEN
                AT_NODES=.TRUE.
              ENDIF
            ENDIF
            CALL PARSE_FACES(NFFACE,NFLIST,noco-1,NRLIST,NTCO,CO,
     '        ERROR,*9999)
            NORMAL_FORCE=.FALSE.
            IF(CBBREV(CO,'FORCE',4,noco+1,NTCO,N3CO)) NORMAL_FORCE=
     &        .TRUE.
            UNDEFORMED=.TRUE.
            IF(CBBREV(CO,'DEFORMED',3,noco+1,NTCO,N3CO)) UNDEFORMED=
     &        .FALSE.
            SCALE_F=1.d0
            IF(CBBREV(CO,'SCALE',3,noco+1,NTCO,N3CO)) SCALE_F=
     &        RFROMC(CO(N3CO+1))
          ENDIF
        ENDIF
        IF(CBBREV(CO,'GRAVITY',4,noco+1,NTCO,N3CO)) THEN
          GRAVFACT=RFROMC(CO(N3CO+1))
C... KSB: Adding gravity vector direction          
          DO nj=1,NJT
            grav_vect(nj)=0.d0 !initialise
          ENDDO
          IF(CBBREV(CO,'UPRIGHT',4,noco+1,NTCO,N3CO)) THEN
            POSITION='UPRIGHT'
            grav_vect(3)=-1.d0
          ELSEIF(CBBREV(CO,'INVERTED',4,noco+1,NTCO,N3CO)) THEN
            POSITION='INVERTED'
            grav_vect(3)=1.d0
          ELSEIF(CBBREV(CO,'PRONE',4,noco+1,NTCO,N3CO)) THEN
            POSITION='PRONE'
            grav_vect(2)=1.d0
          ELSEIF(CBBREV(CO,'SUPINE',4,noco+1,NTCO,N3CO)) THEN
            POSITION='SUPINE'
            grav_vect(2)=-1.d0
          ENDIF
        ELSE
          GRAVFACT=0
          POSITION='UPRIGHT'
        ENDIF
      ENDIF
      
      IF(TYPE(1:5).EQ.'FIELD') THEN
        IF(.NOT.AT_NODES.AND.(.NOT.AT_ELEMENTS))THEN
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_FIT,ERROR,*9999)
          CALL ASSERT(nx.GT.0,'>>No nx defined for this fit class',
     '      ERROR,*9999)
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
C CPB Should loop over number of field variables here
              DO nhx=1,NHP(np,nr,nx)
                nh=NH_LOC(nhx,nx)
                DO nv=1,NVHP(nh,np,1,nr)
                  DO nk=1,NKJ(NJ_LOC(NJL_FIEL,1,nr),np)
                    ny=NYNP(nk,nv,nh,np,0,1,nr)
                    IF(FIX(ny,1,nx)) THEN !ny is a fixed b.c.
                      YP(ny,1,nx)=XP(nk,nv,NJ_LOC(NJL_FIEL,1,nr),np)
                    ENDIF
                  ENDDO !nk loop
                ENDDO !nv loop
              ENDDO !nh loop
            ENDDO !nonode loop
          ENDDO !no_nrlist loop
        ELSE IF(AT_NODES)THEN
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.GT.0,'>>No nx defined for this fit class',
     '      ERROR,*9999)
          IF(DISPLACEMENT)THEN
            nc=1
            niy=2
          ELSE
            nc=2
C MHT 8/11/10 for force initials
            IF(INCR)THEN
              niy=2 !adds to force increments (BCs)
            ELSE
              niy=3 !adds to force initial
 ! note: also do for niy=1, because when written out to file IPINIT_ELAS
              ! transfers the current solution from YP to ZP using YPZP
            ENDIF

C            niy=2 !adds to force increments
          ENDIF
          maxdemens(1)=1000.d0
          maxdemens(2)=-1000.d0
          maxdemens(1)=1000.d0
          maxdemens(2)=-1000.d0
          maxdemens(3)=1000.d0
          maxdemens(4)=-1000.d0
          maxdemens(5)=1000.d0
          maxdemens(6)=-1000.d0
          DO nonode=1,NPLIST1(0)
            np=NPLIST1(nonode)
            maxdemens(1)=min(maxdemens(1),xp(1,1,1,np))
            maxdemens(2)=max(maxdemens(2),xp(1,1,1,np))
            maxdemens(3)=min(maxdemens(3),xp(1,1,2,np))
            maxdemens(4)=max(maxdemens(4),xp(1,1,2,np))
            maxdemens(5)=min(maxdemens(5),xp(1,1,3,np))
            maxdemens(6)=max(maxdemens(6),xp(1,1,3,np))
          ENDDO !np
          DO nonode=1,NPLIST1(0)
            np=NPLIST1(nonode)
            IF(POSITION(1:7).EQ.'UPRIGHT')THEN
              grav=PULMAT(1)*ABS(XP(1,1,3,np)-maxdemens(6))*9810.d0
     &          *gravfact
            ELSE IF(POSITION(1:8).EQ.'INVERTED')THEN
              grav=PULMAT(1)*ABS(XP(1,1,3,np)-maxdemens(5))*9810.d0
     &          *gravfact
            ELSE IF(POSITION(1:5).EQ.'PRONE')THEN
              grav=PULMAT(1)*ABS(XP(1,1,2,np)-maxdemens(3))*9810.d0
     &          *gravfact
            ELSE IF(POSITION(1:6).EQ.'SUPINE')THEN
              grav=PULMAT(1)*ABS(XP(1,1,2,np)-maxdemens(4))*9810.d0
     &          *gravfact
            ELSE
             grav=0.d0
            ENDIF
            nhx=depvar
            nh=NH_LOC(nhx,nx)
            DO nv=1,NVHP(nh,np,1,nr)
              DO nk=1,NKJ(nj_field,np)
                ny=NYNP(nk,nv,nh,np,0,nc,nr)
                IF(ADD_TO)THEN
                  YP(ny,niy,nx)=YP(ny,niy,nx)+XP(nk,nv,nj_field,np)
     &              *SCALE_F
                  IF(DISPLACEMENT)THEN
                    ZP(nk,nv,nh,np,nc)=ZP(nk,nv,nh,np,nc)+
     &                XP(nk,nv,nj_field,np)*SCALE_F
                  ENDIF
                  IF(.NOT.DISPLACEMENT)THEN
                    YP(ny,1,nx)=YP(ny,1,nx)+XP(nk,nv,nj_field,np)
     &                *SCALE_F
                    ZP(nk,nv,nh,np,nc)=ZP(nk,nv,nh,np,nc)+
     &                XP(nk,nv,nj_field,np)*SCALE_F
                  ENDIF
                ELSE
                  YP(ny,niy,nx)=XP(nk,nv,nj_field,np)*SCALE_F+grav
                  IF(.NOT.DISPLACEMENT)THEN !force
                     IF(.NOT.BC)THEN
                        YP(ny,1,nx)=XP(nk,nv,nj_field,np)*SCALE_F
                        ZP(nk,nv,nh,np,nc)=XP(nk,nv,nj_field,np)*SCALE_F
                     ENDIF
                  ENDIF
                ENDIF
              ENDDO !nk loop
            ENDDO !nv loop
          ENDDO !nonode loop
        ELSE IF(AT_ELEMENTS)THEN
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.GT.0,'>>No nx defined for this fit class',
     '      ERROR,*9999)
          DO noelem=1,NELIST(0)
            ne=NELIST(noelem)
            nb=NBJ(nj_field,ne)
            np=NPNE(2,nb,ne)
            DO nhx=1,NHP(np,nr,nx) !nx is solve class
              nh=NH_LOC(nhx,nx)
              ny=NYNE(1,nh,0,1,ne)
              YP(ny,1,nx)=XAB(nj_field,ne)
            ENDDO !nh loop
          ENDDO !noelem loop
        ENDIF

      ELSE IF(TYPE(1:4).EQ.'GRID') THEN
          
        IF(CBBREV(CO,'ALPHA',3,noco+1,NTCO,N3CO)) THEN
          ALPHA=RFROMC(CO(N3CO+1))
        ELSE
          ALPHA=1.0d0
        ENDIF
        
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx_ext,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx_ext.GT.0,'>>No nx defined for GRID class',
     '    ERROR,*9999)

        !boundary element class
        nxc=NXLIST(2)
        CALL NX_LOC(NX_INQUIRE,nxc,nx_bem,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx_bem.GT.0,'>>No nx defined for BEM class',
     '    ERROR,*9999)

        !extracellular update class
        nxc=NXLIST(3)
        CALL NX_LOC(NX_INQUIRE,nxc,nx_upd,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx_upd.GT.0,'>>No nx defined for GRID class',
     '    ERROR,*9999)

        !set up regions
        CALL ASSERT(NRLIST(0).GT.1,'>>Must specify 2 regions',
     '    ERROR,*9999)
        !torso boundary element region
        nr_bem=NRLIST(1)
        !heart boundary element region
        nr_bem2=NRLIST(2)

        !calculate the nodes which will be updated
        NPLIST1(0)=0
        DO npp=1,NPNODE(0,nr_bem2)
          np=NPNODE(npp,nr_bem2)
          DO nrr=1,NP_INTERFACE(np,0)
            IF(NP_INTERFACE(np,nrr).EQ.nr_bem) THEN
              NPLIST1(0)=NPLIST1(0)+1
              NPLIST1(NPLIST1(0))=np
            ENDIF
          ENDDO
        ENDDO

        !find the niq potential index
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqV,NIQ_V,ERROR,*9999)
        IF(niqV.EQ.0) niqV=1

        CALL ASSERT(.NOT.UP_NQNP,' >>Need to update node grid',
     '    ERROR,*9999)
        CALL ASSERT(.NOT.UP_NENQ,' >>Need to update element grid',
     '    ERROR,*9999)

        !Check that potential boundary conditions
        IF(CPLST(0,1).GT.0) THEN
          DO npp=1,CPLST(0,1)
            np=CPLST(npp,1)
            DO np2=1,NPNODE(0,nr_bem)
              IF(np.EQ.NPNODE(np2,nr_bem)) THEN
                nv=1
                nh=NH_LOC(1,nx_bem)
                nc=1 !look for a potential boundary condition
                ny=NYNP(1,nv,nh,np,0,nc,nr_bem)
                IF(ny.GT.0) THEN !may not be in region nr_bem!
                  CALL ASSERT(FIX(ny,1,nx_bem),
     '              '>>Must set potential bc at excluded nodes',
     '              ERROR,*9999)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDIF

        !loop over nodes and update boundary conditions
        ERROR_FLAG=.FALSE.
        CALL ASSERT(ALPHA.GT.-ZERO_TOL,'Alpha must be [0-1]',
     '    ERROR,*9999)
        ONEMINUSALPHA=1.0d0-ALPHA
        CALL ASSERT(ONEMINUSALPHA.GT.-ZERO_TOL,'Alpha must be [0-1]',
     '    ERROR,*9999)

C$OMP   PARALLEL DO
C$OMP&  PRIVATE(DERIV,D2PHIDNDS,DPHIDS,FLUX,nc,nh,np,npp,nq,nv,
C$OMP&    ny,nyd,PHI)
C$OMP&  SHARED(ALPHA,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,ERROR,ERROR_FLAG,
C$OMP&    FIX,FIXQ,IODI,NENQ,NH_LOC,niqV,NKH,
C$OMP&    NPLIST1,NQNP,NQS,NQXI,nr_bem,nx_bem,nx_ext,nx_upd,NXQ,
C$OMP&    NYNP,ONEMINUSALPHA,XQ,YP,YQ)
        DO npp=1,NPLIST1(0)
          IF(.NOT.ERROR_FLAG) THEN
            np=NPLIST1(npp)
            nq=NQNP(np)

C LKC 10-JUL-2000 Looks like this whole block is unused
C
C            !Do a check on BEM normal reversals
CC MHT 24-03-00 REVERSED set but not used
CC            REVERSED=.FALSE.
C            DO nep=1,NENP(np,0,nr_bem)
C              ne=NENP(np,nep,nr_bem)
CC              IF(NW(ne,3).EQ.1) REVERSED=.TRUE.
C            ENDDO !nep

            !Update boundary conditions
            nv=1
            nh=NH_LOC(1,nx_bem)
            nc=1 !look for a potential boundary condition
            ny=NYNP(1,nv,nh,np,0,nc,nr_bem)

            IF(NKH(nh,np,nc,nr_bem).GT.1) THEN
              DERIV=.TRUE.
            ELSE
              DERIV=.FALSE.
            ENDIF
            IF(DERIV) nyd=NYNP(2,nv,nh,np,0,nc,nr_bem)

            IF(FIX(ny,1,nx_bem)) THEN
              PHI=YQ(nq,niqV,1,nx_ext)
              IF(DERIV) THEN
                CALL NQDS(2,1,NENQ,1,niqV,nq,NQS,NQXI,NXQ,AQ,
     &            CQ(6,nq,nx_ext),DNUDXQ,DPHIDS,DXDXIQ,DXDXIQ2,XQ,
     &            YQ(1,1,1,nx_ext),.FALSE.,ERROR,*100)
              ENDIF

              IF(FIXQ(nq,3,nx_upd)) THEN
                !dphi/ds and dphi/dn are reversed
                YP(ny,1,nx_bem)=(ONEMINUSALPHA*YP(ny,1,nx_bem))+
     '            (ALPHA*PHI)
                IF(DERIV) YP(nyd,1,nx_bem)=
     '            (ONEMINUSALPHA*YP(nyd,1,nx_bem))-(ALPHA*DPHIDS)
              ELSE
                !dphi/dn and d2phi/dnds are reversed
                YP(ny,1,nx_bem)=(ONEMINUSALPHA*YP(ny,1,nx_bem))+
     '            (ALPHA*PHI)
                IF(DERIV) YP(nyd,1,nx_bem)=
     '            (ONEMINUSALPHA*YP(nyd,1,nx_bem))+(ALPHA*DPHIDS)
              ENDIF
            ELSE
              nc=2 !look for a flux boundary condition
              ny=NYNP(1,nv,nh,np,0,nc,nr_bem)
              IF(NKH(nh,np,nc,nr_bem).GT.1) THEN
                DERIV=.TRUE.
              ELSE
                DERIV=.FALSE.
              ENDIF
              IF(DERIV) nyd=NYNP(2,nv,nh,np,0,nc,nr_bem)

              IF(FIX(ny,1,nx_bem)) THEN
                CALL GGRADPHIQDN(NENQ,niqV,nq,NQS,NQXI,NXQ,AQ,
     &            CQ(6,nq,nx_ext),DNUDXQ,DXDXIQ,DXDXIQ2,FLUX,
     &            YQ(1,1,1,nx_ext),ERROR,*100)
                IF(DERIV) THEN
                  CALL NQDS(2,1,NENQ,1,niqV,nq,NQS,NQXI,NXQ,AQ,
     &              CQ(6,nq,nx_ext),DNUDXQ,D2PHIDNDS,DXDXIQ,DXDXIQ2,XQ,
     &              YQ(1,1,1,nx_ext),.TRUE.,ERROR,*100)
                ENDIF

                IF(FIXQ(nq,3,nx_upd)) THEN
                  !dphi/ds and dphi/dn are reversed
                  YP(ny,1,nx_bem)=(ONEMINUSALPHA*YP(ny,1,nx_bem))-
     '              (ALPHA*FLUX)
                  IF(DERIV) YP(nyd,1,nx_bem)=
     '              (ONEMINUSALPHA*YP(nyd,1,nx_bem))+(ALPHA*D2PHIDNDS)
                ELSE
                  !dphi/dn and d2phi/dnds are reversed
                  YP(ny,1,nx_bem)=(ONEMINUSALPHA*YP(ny,1,nx_bem))-
     '              (ALPHA*FLUX)
                  IF(DERIV) YP(nyd,1,nx_bem)=
     '              (ONEMINUSALPHA*YP(nyd,1,nx_bem))-(ALPHA*D2PHIDNDS)
                ENDIF
              ELSE
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(UPINIT_1)
                WRITE(OP_STRING,'('' >>ERROR no boundary condition '
     '            //'set'')')
                CALL WRITES(IODI,OP_STRING,ERROR,*100)
CC$OMP END CRITICAL(UPINIT_1)
                GOTO 100
              ENDIF
            ENDIF
            GO TO 102
            !this statement is designed to be skipped if no error
            !occurs. However if a error occurs within a subroutine
            !the alternate return jumps to line 100 to set the flag
 100        CONTINUE
C$OMP CRITICAL(UPINIT_2)
            ERROR_FLAG=.TRUE.
            WRITE(OP_STRING,'(/'' >>ERROR: An error occurred - '
     '        //'results may be unreliable'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*101)
 101        CONTINUE
C$OMP END CRITICAL(UPINIT_2)
 102        CONTINUE
          ENDIF !error_flag
        ENDDO
C$OMP END PARALLEL DO

      ELSE IF(TYPE(1:5).EQ.'ITERATE') THEN
        IF(KTYP20.EQ.1) THEN
          IF(KTYP21.EQ.1) THEN !Forward problem calculation
C Update initial conditions from fitted potential field value based
C on the fit variables to deterimine which type nodal conditions to
C transfer; i.e. if a particular ny is in the fit then that ny is
C transfered from the XP array into the initial condition arrays.
C Other initial condition values are assumed to be already set from
C a previous call to DEINIT.
C NOTE : it is assumed that the potential values are sitting at
C nj=NJ_LOC(NJL_FIEL,1,nr) in the XP array.

C CPB 11/2/94 I'm not sure if the region stuff is being handled
C correctly

C CPB 8/6/94 Adding NX_LOC
            CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_FIT,ERROR,*9999)
            CALL ASSERT(nx.GT.0,'>>No nx defined for this fit class',
     '        ERROR,*9999)

            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
C GMH 8/1/97 Update cmgui link
                CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                DO nhx=1,NHP(np,nr,nx)
                  nh=NH_LOC(nhx,nx)
                  DO nv=1,NVHP(nh,np,1,nr)
                    DO nk=1,NKJ(NJ_LOC(NJL_FIEL,1,nr),np)
                      ny=NYNP(nk,nv,nh,np,0,1,nr)
                      IF(NONY(0,ny,2,nr,nx).GT.0) THEN !ny is in the fit
                        IF(FIX(ny,1,nx)) THEN      !ny is a fixed b.c.
                          YP(ny,1,nx)=XP(nk,nv,NJ_LOC(NJL_FIEL,1,nr),np)
                        ELSE
                          WRITE(OP_STRING,'('' >>WARNING: Var ny='','
     '                      //'I4,'' was in the fit but is not a '','
     '                      //'''fixed boundary condition'')') ny
                          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                        ENDIF
                      ENDIF
                    ENDDO !nk loop
                  ENDDO !nv loop
                ENDDO !nh loop
              ENDDO !nonode loop
            ENDDO !no_nrlist loop
          ELSE IF(KTYP21.EQ.2) THEN
          ENDIF
        ELSE IF(KTYP20.EQ.2) THEN
          IF(KTYP21.EQ.1) THEN
          ELSE IF(KTYP21.EQ.2) THEN
          ENDIF
        ENDIF
      ELSEIF(TYPE(1:9).EQ.'ZCROSSING') THEN
        IF(CBBREV(CO,'INTERPOLATE',2,noco+1,NTCO,N3CO)) THEN
          INTERPOLATE=.TRUE.
        ELSE
          INTERPOLATE=.FALSE.
        ENDIF

        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)
        CALL ASSERT(NPLIST3(0).GT.0,'>>Call define transfer first',
     '    ERROR,*9999)
        CALL ASSERT(EVALUATE_TRANSFER,'>>Evaluate transfer first',
     '    ERROR,*9999)

        irow=0
        DO nonode=1,NPLIST3(0) !list of heart nodes
          np=NPLIST3(nonode)
          DO nhx=1,NHP(np,TRSF_NR_FIRST,nx)
            nh=NH_LOC(nhx,nx)
            DO nv=1,NVHP(nh,np,1,TRSF_NR_FIRST)
C             !Currently we only have position information
C             DO nk1=1,MAX(NKH(nh1,np1,1,TRSF_NR_FIRST)-
C     '       KTYP93(1,TRSF_NR_FIRST),1)
              nk=1
              ny=NYNP(nk,nv,nh,np,0,1,TRSF_NR_FIRST)
              irow=irow+1 !Appropriate row of ZCROSSING

              nts=1
              FOUND=.FALSE.
              DO WHILE(.NOT.FOUND.AND.nts.LT.NTST)
! AJP 9/9/99      IF(ZCROSSING(irow,nts).LT.0) THEN
                IF(ZCROSSING(irow,nts).LE.0.0d0) THEN
                  nts=nts+1
                ELSE
                  FOUND=.TRUE.
                ENDIF
              ENDDO

              CALL ASSERT(nts.LT.NTST,
     '          '>>Cannot find zero-crossing time',ERROR,*9999)
              YP(ny,1,nx)=CALC_TIME_FROM_SAMPLE(nts)
C GBS 16-Nov-2001 Added calculation to interpolate exact zcrossing time
              IF(INTERPOLATE) THEN
                A=ZCROSSING(irow,nts-1)
                B=ZCROSSING(irow,nts)
                YP(ny,1,nx)=YP(ny,1,nx)-(B/(B-A))/TRSF_FREQUENCY
              ENDIF
              FIX(ny,1,nx)=.FALSE.
              YP(ny,3,nx)=YP(ny,1,nx)  ! initial zcrossing soln retained

C*** YP(ny,1,nx) contains initial estimate of activation time at np
C*** Critical point times are fixed in "def opti"
            ENDDO ! nv
          ENDDO ! nhx (nh)
        ENDDO ! nonode (np)
      ELSE IF(TYPE(1:9).EQ.'DATAFIELD')THEN
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this fit class',
     '    ERROR,*9999)

        DO nd=1,NDT
          np=LD_NP(nd) !coupling from data point to node
          IF(AT_NODES)THEN
            DO nhx=1,NHP(np,nr,nx) !nx is solve class
              nh=NH_LOC(nhx,nx)
              DO nv=1,NVHP(nh,np,1,nr)
                nk=1
                ny=NYNP(nk,nv,nh,np,0,1,nr)
                YP(ny,1,nx)=ZD(nj_field,nd)
              ENDDO !nv loop
            ENDDO !nh loop
          ELSE IF(AT_ELEMENTS)THEN
            IF(NENP(np,0,nr).NE.1)THEN
              WRITE(OP_STRING,'('' >>WARNING: error in mapping  '')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            ne=NENP(np,1,nr)
            nb=NBJ(nj_field,ne)
            DO nhx=1,NHP(np,nr,nx) !nx is solve class
              nh=NH_LOC(nhx,nx)
              DO nv=1,NVHP(nh,np,1,nr)
                ny=NYNE(1,nh,0,1,ne)
                YP(ny,1,nx)=ZD(nj_field,nd)
              ENDDO !nv loop
            ENDDO !nh loop
          ENDIF
        ENDDO !nd
        
      ELSE IF(TYPE(1:8).EQ.'CONSTANT') THEN
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)
        IF(AT_NODES)THEN
c          nx=1 !fudge for Kerry at moment (7 May 2004)
          nr=NRLIST(1)
          IF(COMPLIANCE.EQ.3.OR.COMPLIANCE.EQ.6) THEN !Use inlet node co-ords as reference height for grav
            CALL ASSERT(np_in.NE.0,
     &        '>>Inlet node not defined: use INLET_REF_NODE [np#]',
     &        ERROR,*9999)
            maxdemens(1)=XP(1,1,1,np_in)
            maxdemens(2)=XP(1,1,1,np_in)
            maxdemens(3)=XP(1,1,2,np_in)
            maxdemens(4)=XP(1,1,2,np_in)         
            maxdemens(5)=XP(1,1,3,np_in)
            maxdemens(6)=XP(1,1,3,np_in)
          ELSE
C Set up max/min dimensions
            maxdemens(1)=1000.d0
            maxdemens(2)=-1000.d0
            maxdemens(1)=1000.d0
            maxdemens(2)=-1000.d0
            maxdemens(3)=1000.d0
            maxdemens(4)=-1000.d0
            maxdemens(5)=1000.d0
            maxdemens(6)=-1000.d0
            DO nonode=1,NPLIST1(0)
              np=NPLIST1(nonode)
              maxdemens(1)=min(maxdemens(1),xp(1,1,1,np))
              maxdemens(2)=max(maxdemens(2),xp(1,1,1,np))
              maxdemens(3)=min(maxdemens(3),xp(1,1,2,np))
              maxdemens(4)=max(maxdemens(4),xp(1,1,2,np))
              maxdemens(5)=min(maxdemens(5),xp(1,1,3,np))
              maxdemens(6)=max(maxdemens(6),xp(1,1,3,np))
            ENDDO !np
          ENDIF !COMPLIANCE.EQ.3.OR.COMPLIANCE.EQ.6: for blood flow problems
          DO nonode=1,NPLIST1(0)
            np=NPLIST1(nonode)
             IF(COMPLIANCE.EQ.3.OR.COMPLIANCE.EQ.6) THEN
              grav=0.d0
             IF(COUPLE_VIA_LPM.EQ.'Y')THEN!the outlet node is the venous trunk and we don't want to adjust for gravity
              ELSE
              DO nj=1,NJT
                grav=grav+(PULMAT(1)*grav_vect(nj)*9810.d0*gravfact
     &            *(XP(1,1,nj,np)-XP(1,1,nj,np_in)))
              ENDDO !nj
             ENDIF!coupled via LPM
            ELSE
              IF(POSITION(1:7).EQ.'UPRIGHT')THEN
                grav=PULMAT(1)*ABS(XP(1,1,3,np)-maxdemens(6))*9810.d0
     &            *gravfact
              ELSE IF(POSITION(1:8).EQ.'INVERTED')THEN
                grav=PULMAT(1)*ABS(XP(1,1,3,np)-maxdemens(5))*9810.d0
     &            *gravfact
              ELSE IF(POSITION(1:5).EQ.'PRONE')THEN
                grav=PULMAT(1)*ABS(XP(1,1,2,np)-maxdemens(3))*9810.d0
     &            *gravfact
              ELSE IF(POSITION(1:6).EQ.'SUPINE')THEN
                grav=PULMAT(1)*ABS(XP(1,1,2,np)-maxdemens(4))*9810.d0
     &            *gravfact
              ELSE
                grav=0.d0
              ENDIF
            ENDIF !COMPLIANCE.EQ.3.OR.COMPLIANCE.EQ.6
            DO nhx=1,NHP(np,nr,nx) !nx is solve class
              nh=NH_LOC(nhx,nx)
              DO nv=1,NVHP(nh,np,1,nr)
                ny=NYNP(1,nv,nh,np,0,1,nr)
                YP(ny,1,nx)=constant_value+grav
              ENDDO !nv loop
            ENDDO !nh loop
          ENDDO !nonode loop
        ELSE IF(AT_ELEMENTS)THEN
          nr=NRLIST(1)
          nh=NH_LOC(1,nx)
          DO noelem=1,NELIST(0)
            ne=NELIST(noelem)
            ny=NYNE(1,nh,1,1,ne)
            YP(ny,1,nx)=constant_value
          ENDDO
c          WRITE(OP_STRING,'('' >> Not implemented'')')
c          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

      ELSE IF(TYPE(1:6).EQ.'NORMAL')THEN
        nr=NRLIST(1)
        IF(NORMAL_FORCE)THEN
          nc=2 !forces
        ELSE
          nc=1 !displacements
        ENDIF
        NPLIST1(0)=0
        DO np=1,NP_R_M
          NFNP(0,np)=0 !initialising
        ENDDO !np
C Get a list of nodes in the face list. 
C Loop over each face in the list of faces        
        DO noface=1,NFLIST(0) 
          nf=NFLIST(noface)
          nef=NPF(8,nf) !local face number
          ne=NPF(6,nf) !host element
          nb=NBJ(1,ne)
C Loop over the face nodes          
          DO nne=1,NNF(0,nef,nb)
            nn=NNF(1+nne,nef,nb) !nn'th element node
            np=NPNE(nn,nb,ne) !global node at nn'th position in element
            NFNP(0,np)=NFNP(0,np)+1
            NFNP(NFNP(0,np),np)=nf
C Check that the face node has not already been done            
            IF(.NOT.INLIST(np,NPLIST1(1),MIN(NPLIST1(0),NPM),n1))THEN
              NPLIST1(0)=NPLIST1(0)+1
              IF(NPLIST1(0).LE.NPM) NPLIST1(NPLIST1(0))=np
            ENDIF !INLIST
          ENDDO !nne
        ENDDO !noface

C For list of nodes, average the normal calculated at XP for each face.
        DO nonode=1,NPLIST1(0)
          DO nj=1,3
            SUM_XNORM(nj)=0.d0
          ENDDO !nj
          AREA=0.d0
          np=NPLIST1(nonode)
          DO noface=1,NFNP(0,np) !for each surrounding face
            nf=NFNP(noface,np) !global face number
            ne=NPF(6,nf) !host element
C Get nodal information for host element
            IF(UNDEFORMED)THEN
              CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     &          nr,NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,
     &          *9999)
            ELSE
              CALL ZPZE(NBH(1,1,ne),1,NHE(ne,nxc),NKHE(1,1,1,ne),
     '          NPF(1,1),NPNE(1,1,ne),NRE(ne),NVHE(1,1,1,ne),
     '          NW(ne,1,nxc),nxc,CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '          ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
            ENDIF
C Get the face directions, and 0 or 1 for current face              
c            XI1=NNF(1,nef,nb) !normal to face
            in1=NPF(1,nf) !first Xi direction for face nf
            in2=NPF(3,nf) !second Xi direction for face nf
            nn=1
            DO WHILE(NPNE(nn,nb,ne).NE.np)
              nn=nn+1
            ENDDO !while
C Find XI location, must be at a corner            
            IF(nn.EQ.2.OR.nn.EQ.4.OR.nn.EQ.6.OR.nn.EQ.8) XI(1)=1.d0
            IF(nn.EQ.3.OR.nn.EQ.4.OR.nn.EQ.7.OR.nn.EQ.8) XI(2)=1.d0
            IF(nn.EQ.5.OR.nn.EQ.6.OR.nn.EQ.7.OR.nn.EQ.8) XI(3)=1.d0

            IF(UNDEFORMED)THEN
              DO nj=1,3
                XI(3)=0.d0 !initialise XI location of node
                XI05(nj)=0.5d0 !center of the 3D element
                X(nj)=XP(1,1,nj,np) !nodal location
              ENDDO !nj
              DO nj=1,3 !note that 3D elements assumed throughout
                X05(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '            XI05,XE(1,nj)) !global coordinates of element mid-point
                DX(nj,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '            nb,2,XI,XE(1,nj)) !global Xi1 derivatives
                DX(nj,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '            nb,4,XI,XE(1,nj)) !global Xi2 derivatives
                DX(nj,3)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '            nb,7,XI,XE(1,nj)) !global Xi3 derivatives
              ENDDO !nj
            ELSE
              DO nj=1,3
                XI(3)=0.d0 !initialise XI location of node
                XI05(nj)=0.5d0 !center of the 3D element
                X(nj)=ZP(1,1,nj,np,1) !nodal location
              ENDDO !nj
              DO nj=1,3 !note that 3D elements assumed throughout
                X05(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '            XI05,ZE(1,nj)) !global coordinates of element mid-point
                DX(nj,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '            nb,2,XI,ZE(1,nj)) !global Xi1 derivatives
                DX(nj,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '            nb,4,XI,ZE(1,nj)) !global Xi2 derivatives
                DX(nj,3)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '            nb,7,XI,ZE(1,nj)) !global Xi3 derivatives
              ENDDO !nj
            ENDIF !UNDEFORMED
            
C Calculate unit outward normal to boundary surface
            XNORM(1)=(DX(2,in1)*DX(3,in2)-DX(3,in1)*DX(2,in2))
            XNORM(2)=(DX(3,in1)*DX(1,in2)-DX(1,in1)*DX(3,in2))
            XNORM(3)=(DX(1,in1)*DX(2,in2)-DX(2,in1)*DX(1,in2))
            CALL NORMALISE(3,XNORM,ERROR,*9999) !normalise XNORM
            DIST1=0.d0
            DIST2=0.d0
            DO nj=1,3
              DIST1=DIST1+(X(nj)+XNORM(nj)-X05(nj))**2.d0
              DIST2=DIST2+(X(nj)-XNORM(nj)-X05(nj))**2.d0
            ENDDO !nj
            DIST1=DSQRT(DIST1)
            DIST2=DSQRT(DIST2)
            IF(DIST1.LT.DIST2)THEN !determine outwards direction
              DO nj=1,3
                XNORM(nj)=-XNORM(nj)
              ENDDO !nj
            ENDIF
C Add contribution to outwards normal from each face            
            DO nj=1,3
              SUM_XNORM(nj)=SUM_XNORM(nj)+XNORM(nj)
            ENDDO !nj
            AREA=AREA+DF(nf) !sum the area of the adjoining faces
          ENDDO !noface

          CALL NORMALISE(3,SUM_XNORM,ERROR,*9999) !normalise XNORM
          
C Put normal into ny for displacement or force
          DO nhx=1,NJ_LOC(NJL_GEOM,0,nr) !note, only geometry
            nh=NH_LOC(nhx,nxc)
            DO nv=1,NVHP(nh,np,1,nr)
              ny=NYNP(1,nv,nh,np,0,nc,nr) !ny for displacement or force B.C.
              IF(FIX(ny,1,nxc))THEN
                IF(NORMAL_FORCE)THEN
                  YP(ny,2,nxc)=SUM_XNORM(nhx)*AREA*SCALE_F
C                    YP(ny,2,nxc)=SUM_XNORM(nhx)
                ELSE
                  YP(ny,2,nxc)=SUM_XNORM(nhx)*SCALE_F
                ENDIF
              ENDIF !FIX
            ENDDO !nv
          ENDDO !nhx
        ENDDO !nonode
      ENDIF

      CALL_INIT=.TRUE.
      CALL EXITS('UPINIT')
      RETURN
 9999 CALL ERRORS('UPINIT',ERROR)
      CALL EXITS('UPINIT')
      RETURN 1
      END
