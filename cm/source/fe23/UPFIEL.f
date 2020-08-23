      SUBROUTINE UPFIEL(IBT,IDO,INP,NBH,NBHF,NBJ,NBJF,
     '  NEELEM,NEL,NENP,NFF,NFFACE,
     '  NFLIST,NKB,NKEF,NKH,NKHE,NKJ,NKJE,NLF,NLL,NLLINE,NLLIST,
     '  NNF,NNL,NPF,NPL,NPLIST,NPNE,NPNF,
     '  NPNODE,NPNY,NRE,NRLIST,NVHE,NVHP,NVJE,NVJF,
     '  NVJP,NWP,NXI,NXLIST,NYNO,NYNE,NYNP,
     '  CP,DF,PAOPTI,PG,RG,SE,SF,WG,XA,XAB,XE,XG,
     '  XP,YP,ZD,STRING,FIX,ERROR,*)

C#### Subroutine: UPFIEL
C###  Description:
C###    Updates field variables from geometry or solution variables.

C CPB 20/5/94 Havn't put in much in the way of mapping multiple njs
C on to multiple field njs

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'opti00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),
     '  NBJ(NJM,NEM),NBJF(NJM,NFM),NEELEM(0:NE_R_M,0:NRM),
     '  NEL(0:NELM,NLM),NENP(NPM,0:NEPM,0:NRM),NFF(6,NEM),
     '  NFFACE(0:NF_R_M,NRM),NFLIST(0:NFM),NKB(2,2,2,NNM,NBFM),
     '  NKEF(0:4,16,6,NBFM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),
     '  NLF(4,NFM),NLL(12,NEM),NLLINE(0:NL_R_M,0:NRM),NLLIST(0:NLM),
     '  NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),NPF(9,NFM),NPL(5,0:3,NLM),
     '  NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),
     '  NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM,NXM),NRE(NEM),
     &  NRLIST(0:NRM),NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJF(NNM,NBFM,NJM),
     '  NVJP(NJM,NPM),NWP(NPM,2),NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NXLIST(0:NXM),NYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),NYNE(NAM,NHM,
     '  0:NRCM,NCM,NEM),NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CP(NMM,NPM,NXM),DF(NFM),PAOPTI(NOPM),
     '  RG(NGM),PG(NSM,NUM,NGM,NBM),
     '  SE(NSM,NBFM,NEM),SF(NSM,NBFM),WG(NGM,NBM),
     '  XA(NAM,NJM,NEM),XAB(NORM,NEM),XE(NSM,NJM),XG(NJM,NUM),
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM),
     '  ZD(NJM,NDM)
      CHARACTER ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER FREELIST(NJ_LOC_MX),IBEG,IEND,IFROMC,iy,LD,
     '  N3CO,nb,nc,nd,ne,ne1,
     '  nh,ni,NITB,
     '  nj,nj1,nj2,nj_elem,njj,NJF,nk,nm,nn,
     '  no,no_nrlist,noelem,
     '  nonode,np,nr,nr1,nu(1),NUM_FIELD,numfree,nv,
     '  nx,nxc,nx_opt,ny,PART2
      INTEGER*4 njdir,np_max,np_min,XPT_PTR
      REAL*8 avrad,DIST,grad,href,hmax,hmin,mean,MIN_DIST,PXI,RFROMC,
     '  XD(3),XI(3)
      LOGICAL ABBREV,ALL_REGIONS,CBBREV,ELEM_FIELD,HANGING
      CHARACTER OPERATION*16,TYPE*9,UPDATE*16

C SEN 20/01/03 nx not initialised?
      DATA nx / 0 /

      CALL ENTERS('UPFIEL',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

        UPDATE='HELPFIELD'
        CALL UPFG(UPDATE,%VAL(0),%VAL(0),%VAL(0),%VAL(0),
     '    %VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),
     '    %VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),
     &    STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update field from geometry/solution
C###  Parameter:      <YP_index #[1]>
C###   Specify the YP index to update the field from
C###  Parameter:      <nh nh#[1]>
C###   Specify the dependent variable to update the field from
C###  Parameter:      <to FIELD_VAR#[1]>
C###   Specify the field variable to update
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description: Updates the field variables from geometry or
C###   solution variables

        OP_STRING(1)=STRING(1:IEND)//' from geometry/solution'
        OP_STRING(2)=BLANK(1:15)//'<YP_index #[1]>'
        OP_STRING(3)=BLANK(1:15)//'<nh nh#[1]>'
        OP_STRING(4)=BLANK(1:15)//'<to FIELD_VAR#[1]>'
        OP_STRING(5)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(6)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update field from material
C###  Parameter:      <CP_index #[1]>
C###   Specify the CP index to update the field from
C###  Parameter:      <to FIELD_VAR#[1]>
C###   Specify the field variable to update
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <using (fit/optimisation/solve)[solve]>
C###    Specify the method used to update the field variables
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Update field values with node based material parameter

        OP_STRING(1)=STRING(1:IEND)//' from material'
        OP_STRING(2)=BLANK(1:15)//'<CP_index #[1]>'
        OP_STRING(3)=BLANK(1:15)//'<to FIELD_VAR#[1]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(5)=BLANK(1:15)
     '    //'<using (fit/optimisation/solve)[solve]>'
        OP_STRING(6)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
C#### Command: FEM update field distance
C###  Parameter:      <to FIELD_VAR#[1]>
C###   Specify the field variable to update
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Description:
C###    Updates field values with distance from the nearest data point

C**** Created by Carey Stevens 19/6/98

        OP_STRING(1)=STRING(1:IEND)//' distance'
        OP_STRING(2)=BLANK(1:15)//'<to FIELD_VAR#[1]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update field eikonal
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Parameter:      <basis #[1]>
C###    Specify the basis number
C###  Description:
C###    Sets up the field variable 1 from dependent variable versions
C###    and boundary conditions for solution of an eikonal equation.

        OP_STRING(1)=STRING(1:IEND)//' eikonal'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update field velocity
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Sets up the velocity field variable for lung fluid flow
C###    problems. Calculates velocity from flow and radius fields.

        OP_STRING(1)=STRING(1:IEND)//' velocity'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
C---------------------------------------------------------------------

C#### Command: FEM update field FIELD_VAR#[1] linear
C###  Parameter:      <mean #>
C###    Specify the mean field value.
C###  Parameter:      <gradient #>
C###    Specify the field gradient.
C###  Parameter:      <direction #(nj=1/2/3) [1]>
C###    Specify the coordinate direction to apply the linear gradient.
C###  Description:
C###    Sets up the linear field variable. Takes minimum and maximum 
C###    coordinate values in the nj direction specified and applies
C###    a linear distribution using the mean and gradient.

        OP_STRING(1)=STRING(1:IEND)//' linear'
        OP_STRING(2)=BLANK(1:15)//'<mean #>'
        OP_STRING(3)=BLANK(1:15)//'<gradient #>'
        OP_STRING(4)=BLANK(1:15)//'<direction #(nj=1/2/3) [1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update field FIELD_VAR#[1] summation
C###  Description:
C###    Sums up the field from terminal nodes to inlet node. 
C###    For use with 1D tree meshes (e.g. pulmonary airways).
C###    Performed on all nodes/elements.

        OP_STRING(1)=STRING(1:IEND)//' summation'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
 
      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc','UPFIEL',ERROR,*9999)
      ELSE
        
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)

        IF(CBBREV(CO,'SUBSTITUTE',2,noco+1,NTCO,n3co)) THEN
          TYPE='OPERATION'
          OPERATION='SUBSTITUTE'
          PART2=N3CO
        ELSEIF(CBBREV(CO,'ADD',2,noco+1,NTCO,n3co)) THEN
          TYPE='OPERATION'
          OPERATION='ADD'
          PART2=N3CO
        ELSEIF(CBBREV(CO,'SUBTRACT',2,noco+1,NTCO,n3co)) THEN
          TYPE='OPERATION'
          OPERATION='SUBTRACT'
          PART2=N3CO
        ELSEIF(CBBREV(CO,'MULTIPLY',2,noco+1,NTCO,n3co)) THEN
          TYPE='OPERATION'
          OPERATION='MULTIPLY'
          PART2=N3CO
        ELSEIF(CBBREV(CO,'DIVIDE',2,noco+1,NTCO,n3co)) THEN
          TYPE='OPERATION'
          OPERATION='DIVIDE'
          PART2=N3CO
        ELSEIF(CBBREV(CO,'FROM',2,noco+1,NTCO,n3co)) THEN
          IF(ABBREV(CO(n3co+1),'GEOMETRY',1)) THEN
            TYPE='GEOMETRY'
          ELSE IF(ABBREV(CO(n3co+1),'SOLUTION',1)) THEN
            TYPE='SOLUTION'
          ELSE IF(ABBREV(CO(n3co+1),'MATERIAL',1)) THEN
            TYPE='MATERIAL'
          ELSE IF(ABBREV(CO(n3co+1),'RADIUS',1)) THEN
            TYPE='RADIUS'
          ELSE
            CO(noco+1)='?'
            GO TO 1
          ENDIF
        ELSE IF(CBBREV(CO,'DISTANCE',2,noco+1,NTCO,n3co)) THEN
          TYPE='DISTANCE'
        ELSE IF(CBBREV(CO,'EIKONAL',2,noco+1,NTCO,n3co)) THEN
          TYPE='EIKONAL'
        ELSE IF(CBBREV(CO,'OPTIMISE',2,noco+1,NTCO,n3co)) THEN
          TYPE='OPTIMISE'
        ELSE IF(CBBREV(CO,'VELOCITY',2,noco+1,NTCO,n3co)) THEN
          TYPE='VELOCITY'
        ELSE IF(CBBREV(CO,'LINEAR',2,noco+1,NTCO,n3co)) THEN
          TYPE='LINEAR'
        ELSE IF(CBBREV(CO,'SUMMATION',3,noco+1,NTCO,n3co)) THEN
          TYPE='SUMMATION'
        ELSE
          CO(noco+1)='?'
          GO TO 1
        ENDIF

        IF(CBBREV(CO,'ELEMENTS',2,noco+1,NTCO,n3co)) THEN
          ELEM_FIELD=.TRUE.
        ELSE
          ELEM_FIELD=.FALSE.
        ENDIF
        
        IF((TYPE(1:8).NE.'GEOMETRY').AND.
     '    (TYPE(1:8).NE.'OPTIMISE'))THEN
          IF(.NOT.ELEM_FIELD)THEN
            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              CALL ASSERT(NJ_LOC(NJL_FIEL,0,nr).GT.0,
     '          '>>Define field first',ERROR,*9999)
            ENDDO !nr
          ELSE
            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              CALL ASSERT(NEJ_LOC(0,nr).GT.0,'>>Define field first',
     '          ERROR,*9999)
            ENDDO
          ENDIF
        ENDIF !not geometry

        IF(TYPE(1:9).EQ.'OPERATION') THEN

          UPDATE='FIELD'
          CALL UPFG(UPDATE,CP,%VAL(0),NEELEM,NKH,NKJ,NPLIST,NPNODE,
     '      NRLIST,NVHP,NVJP,NYNP,OPERATION,PART2,XAB,XP,%VAL(0),YP,ZD,
     '      STRING,ERROR,*9999)

        ELSEIF(TYPE(1:8).EQ.'GEOMETRY') THEN
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            IF(NJ_LOC(NJL_FIEL,0,nr).EQ.0) THEN
C             Create new fields.
              NUM_FIELD=NJ_LOC(NJL_GEOM,0,nr)
C             Set up NJ_LOC(nj,njj).
C             Store free nj's in FREELIST & check enough room.
              nj=0
              numfree=0
              DO WHILE((numfree.LT.NUM_FIELD).AND.(nj.LE.NJM))
                nj=nj+1
                CALL ASSERT(nj.LE.3*NJ_LOC_MX,'>>Increase NJ_LOC_MX in '
     '            //'loc00.cmn',ERROR,*9999)
                IF(NJ_TYPE(nj,1).EQ.0) THEN
C                 Empty space in that nj location
                  numfree=numfree+1
                  CALL ASSERT(numfree.LE.NJ_LOC_MX,
     '              '>>Increase NJ_LOC_MX in '//'loc00.cmn',ERROR,*9999)
                  FREELIST(numfree)=nj
                ENDIF
              ENDDO
              CALL ASSERT(nj.LE.NJM,' >>Increase NJM',ERROR,*9999)
C             Store field in free space
              DO numfree=1,NUM_FIELD
                nj=FREELIST(numfree)
                NJ_LOC(NJL_FIEL,numfree,nr)=nj
                NJ_TYPE(nj,1)=NJL_FIEL
                NJ_TYPE(nj,2)=numfree
                IF(nj.GT.NJ_LOC(0,0,nr)) NJ_LOC(0,0,nr)=nj
                IF(nj.GT.NJ_LOC(NJL_FIEL,0,0)) NJ_LOC(NJL_FIEL,0,0)=nj
              ENDDO
              NJ_LOC(NJL_FIEL,0,nr)=NUM_FIELD
              IF(NJ_LOC(0,0,nr).GT.NJ_LOC(0,0,0))
     '          NJ_LOC(0,0,0)=NJ_LOC(0,0,nr)
C             Copy geometry info
              DO njj=1,NJ_LOC(NJL_FIEL,0,nr)
                nj1=NJ_LOC(NJL_GEOM,njj,nr)
                nj2=NJ_LOC(NJL_FIEL,njj,nr)
                DO nonode=1,NPNODE(0,nr)
                  np=NPNODE(nonode,nr)
                  NVJP(nj2,np)=NVJP(nj1,np)
                  NKJ(nj2,np)=NKJ(nj1,np)
                  DO nv=1,NVJP(nj1,np)
                    DO nk=1,NKJ(nj1,np)
                      XP(nk,nv,nj2,np)=XP(nk,nv,nj1,np)
                    ENDDO !nk
                  ENDDO !nv
                ENDDO !nonode (np)
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  nb=NBJ(nj1,ne)
                  NBJ(nj2,ne)=nb
                  DO nn=1,NNT(nb)
                    NVJE(nn,nb,nj2,ne)=NVJE(nn,nb,nj1,ne)
                    DO nk=1,NKT(nn,nb)
                      NKJE(nk,nn,nj2,ne)=NKJE(nk,nn,nj1,ne)
                    ENDDO !nk
                  ENDDO !nn
                ENDDO !noelem (ne)
              ENDDO !njj

C new CS 13/1/2002 Setting up faces for new fields
              CALL FACCAL(IBT,IDO,INP,NBJ,NBJF,NEELEM,NFF,NFFACE,
     '          NKJE,NKEF,NLF,NNF,NNL,NPF,NPL,NPNE,NPNF,
     '          NRE,NVJE,NVJF,DF,PG,RG,SE,SF,WG,XA,
     '          XE,XG,XP,ERROR,*9999)

              CALL_FIEL=.TRUE.
              CALL_ELFD=.TRUE.
            ELSE
              CALL ASSERT(NJ_LOC(NJL_GEOM,0,nr).LE.NJ_LOC(NJL_FIEL,0,
     '          nr),'>>Increase number of field variables',ERROR,*9999)
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                  nj1=NJ_LOC(NJL_GEOM,njj,nr)
                  nj2=NJ_LOC(NJL_FIEL,njj,nr)
                  DO nv=1,NVJP(nj1,np)
                    DO nk=1,NKJ(nj1,np)
                      XP(nk,nv,nj2,np)=XP(nk,nv,nj1,np)
                    ENDDO !nk
                  ENDDO !nv
                ENDDO !njj
              ENDDO !nonode (np)
            ENDIF !NJ_LOC(NJL_FIEL,0,nr)
          ENDDO !nr

        ELSE IF(TYPE(1:8).EQ.'SOLUTION') THEN
          IF(CBBREV(CO,'YP_INDEX',2,noco+1,NTCO,N3CO)) THEN
            iy=IFROMC(CO(N3CO+1))
          ELSE
            iy=1
          ENDIF
          IF(CBBREV(CO,'NH',2,noco+1,NTCO,N3CO)) THEN
            nh=IFROMC(CO(N3CO+1))
          ELSE
            nh=1
          ENDIF
          IF(CBBREV(CO,'TO',2,noco+1,NTCO,N3CO)) THEN
            NJF=IFROMC(CO(N3CO+1))
          ELSE
            NJF=1
          ENDIF
C??? CS why? iy=5..7 should be OK?
C          CALL ASSERT(iy.GT.0.AND.iy.LE.5,'>>IY out of range',
C     '      ERROR,*9999)
          CALL ASSERT(iy.GT.0.AND.iy.LE.7,'>>IY out of range',
     '      ERROR,*9999)
          CALL ASSERT((nh.GT.0).OR.(nh.LT.NH_LOC(0,nx)),
     '       '>>dependent variable out of range',ERROR,*9999)
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '      ERROR,*9999)

          IF(ELEM_FIELD)THEN
            nj=NEJ_LOC(NJF,nr)
            DO noelem=1,neelem(0,nr)
              ne=neelem(noelem,nr)
              ny=NYNE(1,nh,1,1,ne)
              XAB(nj,ne)=YP(ny,iy,nx)
            ENDDO
          ELSE
C KAT 2Jan99: Updating field interpolation as well as values so that
C             field stores solution exactly.
            nc=1 !should be specified on command line
            DO no_nrlist=1,NRLIST(0) !loop over regions
              nr=NRLIST(no_nrlist)
              CALL ASSERT(NJF.GT.0.AND.NJF.LE.NJ_LOC(NJL_FIEL,0,nr),
     '          '>>Field variable number out of range',ERROR,*9999)
              nj=NJ_LOC(NJL_FIEL,NJF,nr)
C              DO no_nynr=1,NYNR(0,0,1,nr,nx) !loop over global vars
C                ny=NYNR(no_nynr,0,1,nr,nx) !is global var number
C                nk=NPNY(1,ny,0,nx)
C                nv=NPNY(2,ny,0,nx)
C                np=NPNY(4,ny,0,nx)
C                CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
C                XP(nk,nv,nj,np)=YP(ny,iy,nx)
C              ENDDO !no_nynr
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                NVJP(nj,np)=NVHP(nh,np,nc,nr)
                NKJ(nj,np)=NKH(nh,np,nc,nr)
                DO nv=1,NVJP(nj,np)
                  DO nk=1,NKJ(nj,np)-KTYP93(nc,nr)
                    ny=NYNP(nk,nv,nh,np,0,nc,nr)
                    XP(nk,nv,nj,np)=YP(ny,iy,nx)
                  ENDDO !nk
                ENDDO !nv
              ENDDO !nonode
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                nb=NBH(nh,nc,ne)
                NBJ(nj,ne)=nb
                DO nn=1,NNT(nb)
                  NVJE(nn,nb,nj,ne)=NVHE(nn,nb,nh,ne)
                  DO nk=1,NKT(nn,nb)
                    NKJE(nk,nn,nj,ne)=NKHE(nk,nn,nh,ne)
                  ENDDO !nk
                ENDDO !nn
              ENDDO !noelem
            ENDDO !no_nrlist
          ENDIF !ELEM_FIELD
          
        ELSE IF(TYPE(1:8).EQ.'MATERIAL') THEN
          IF(CBBREV(CO,'CP_INDEX',2,noco+1,NTCO,N3CO)) THEN
            nm=IFROMC(CO(N3CO+1))
          ELSE
            nm=1
          ENDIF
          CALL ASSERT(nm.GT.0.AND.nm.LE.NMM,'>>nm out of range',
     '      ERROR,*9999)

          IF(CBBREV(CO,'TO',2,noco+1,NTCO,N3CO)) THEN
            NJF=IFROMC(CO(N3CO+1))
          ELSE
            NJF=1
          ENDIF

          IF(CBBREV(CO,'USING',2,noco+1,NTCO,N3CO)) THEN
            IF(ABBREV(CO(N3CO+1),'SOLVE',2)) THEN
              CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
              CALL ASSERT(nx.GT.0,
     '          '>>No nx defined for this solve class',ERROR,*9999)
            ELSE IF(ABBREV(CO(N3CO+1),'FIT',2)) THEN
              CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_FIT,ERROR,*9999)
              CALL ASSERT(nx.GT.0,
     '          '>>No nx defined for this fit class',ERROR,*9999)
            ELSE IF(ABBREV(CO(N3CO+1),'OPTIMISE',2)) THEN
              CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_OPTI,ERROR,*9999)
              CALL ASSERT(nx.GT.0,
     '          '>>No nx defined for this optimisation class',
     '          ERROR,*9999)
            ELSE
              CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
              CALL ASSERT(nx.GT.0,
     '          '>>No nx defined for this solve class',ERROR,*9999)
            ENDIF
          ELSE
            CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx.GT.0,
     '        '>>No nx defined for this solve class',ERROR,*9999)
          ENDIF

          DO no_nrlist=1,NRLIST(0) !loop over regions
            nr=NRLIST(no_nrlist)
            CALL ASSERT(NJF.GT.0.AND.NJF.LE.NJ_LOC(NJL_FIEL,0,nr),
     '        '>>Field variable number out of range',ERROR,*9999)
            nj=NJ_LOC(NJL_FIEL,NJF,nr)
            DO nonode=1,NPNODE(0,nr) !loop over nodes
              np=NPNODE(nonode,nr)
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
              XP(1,1,nj,np)=CP(nm,np,nx)
            ENDDO !no_np
          ENDDO !no_nrlist
 
        ELSE IF(TYPE(1:8).EQ.'DISTANCE') THEN
          IF(CBBREV(CO,'TO',2,noco+1,NTCO,N3CO)) THEN
            NJF=IFROMC(CO(N3CO+1))
          ELSE
            NJF=1
          ENDIF
          HANGING=.FALSE.
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              CALL ASSERT(NJF.GT.0.AND.NJF.LE.NJ_LOC(NJL_FIEL,0,nr),
     '          '>>Field variable number out of range',ERROR,*9999)
              IF(NWP(np,1).EQ.0) THEN ! is not a hanging node
                nj=NJ_LOC(NJL_FIEL,NJF,nr)
                CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                DO nv=1,NVJP(nj,np)
                  MIN_DIST=RMAX
                  DO nd=1,NDT
                    DIST=DSQRT((XP(1,nv,1,np)-ZD(1,nd))**2+
     '                (XP(1,nv,2,np)-ZD(2,nd))**2+
     '                (XP(1,nv,3,np)-ZD(3,nd))**2)
                    IF(DIST.LT.MIN_DIST) THEN
                      MIN_DIST=DIST
                    ENDIF
                  ENDDO
                  XP(1,nv,nj,np)=MIN_DIST
                ENDDO !nv
              ELSE
                HANGING=.TRUE.
              ENDIF !NWP
            ENDDO !nonode (np)
          ENDDO !nr
          IF(HANGING) THEN
            nu(1)=1
            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              IF(NWP(np,1).NE.0) THEN ! is a hanging node
                nj=NJ_LOC(NJL_FIEL,NJF,nr)
                ne=NWP(np,1)

                CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '            NPNE(1,1,ne),nr,NVJE(1,1,1,ne),SE(1,1,ne),
     '            XA(1,1,ne),XE,XP,ERROR,*9999)
                NITB=NIT(NBJ(1,ne))
                XD(1)=XP(1,1,1,np)
                XD(2)=XP(1,1,2,np)
                XD(3)=XP(1,1,3,np)
                DO ni=1,NITB
                  XI(ni)=0.5d0
                ENDDO
                CALL DEXI_POINT(IBT,IDO,INP,LD,NBJ,ne,NITB,nr,
     '            0.d0,XE,XI,XI,XD,.FALSE.,ERROR,*9999)
                CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                DO nv=1,NVJP(nj,np)
                    nb=NBJ(nj,ne)
                    XP(1,nv,nj,np)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                INP(1,1,nb),nb,nu(1),Xi,XE(1,nj))
                ENDDO !nv
              ENDIF !NWP
            ENDDO !nonode (np)
          ENDDO !nr

          ENDIF

        ELSE IF(TYPE(1:6).EQ.'RADIUS')THEN
          nj=NEJ_LOC(IFROMC(CO(noco+1)),nr)
          nb=1   !NBJ(nj_radius,ne)
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            avrad=0.d0
            DO nn=1,2
              np=NPNE(nn,nb,ne)
              nv=NVJE(nn,nb,nj_radius,ne)
              avrad=avrad+XP(1,nv,nj_radius,np)
            ENDDO !nn
            XAB(nej_strtrad,ne)=avrad/2.d0
          ENDDO !nonode (np)
        ELSE IF(TYPE(1:7).EQ.'EIKONAL') THEN
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '      ERROR,*9999)
          CALL ASSERT(CALL_INIT,'>>Define initial conditions first',
     '      ERROR,*9999)
          XPT_PTR=0
          CALL ALLOCATE_MEMORY(NKM,1,DPTYPE,XPT_PTR,
     '      MEM_INIT,ERROR,*9999)
          CALL UPFIEL_EIK(INP,NBH,NBHF,NBJ,NBJF,NEELEM,NEL,NFF,NFFACE,
     '      NFLIST,NKB,NKHE,NKJ,NKJE,NLF,NLL,NLLINE,NLLIST,NNF,NNL,NPF,
     '      NPNE,NPNODE,NRLIST,NVHE,NVHP,NVJE,NVJP,nx,NYNP,XP,
     '      %VAL(XPT_PTR),FIX(1,1,nx),ERROR,*9991)
          CALL FREE_MEMORY(XPT_PTR,ERROR,*9999)
        ELSE IF(TYPE(1:8).EQ.'OPTIMISE') THEN
          CALL NX_LOC(NX_INQUIRE,nxc,nx_opt,NX_OPTI,ERROR,*9999)
C!!! CS 14/4/2000 Looping over multiple regions here,
C!!! but optimising routines don't seem to handle multiple
C!!! regions at the moment though
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            DO no=1,NTOPTI
              ny=NYNO(1,no,1,nr,nx_opt)
              nk=NPNY(1,ny,0,nx_opt)
              nv=NPNY(2,ny,0,nx_opt)
              nj=NPNY(3,ny,0,nx_opt)
              np=NPNY(4,ny,0,nx_opt)
              IF(nk.EQ.1) THEN
                XP(nk,nv,nj,np)=PAOPTI(no)
                DOP=.TRUE.
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' Updated '
     '              //'XP(..,np) '',I4,'' = '',D10.3)')
     '              np,XP(nk,nv,nj,np)*180.0d0/PI
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                 ENDIF
                DOP=.FALSE.
              ELSE
                XP(nk,nv,nj,np)=PAOPTI(no)
              ENDIF
            ENDDO !no
          ENDDO !no_nrlist
        ELSE IF(TYPE(1:8).EQ.'VELOCITY')THEN !AJS Sep 2007
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            CALL ASSERT(nej_vel.NE.0.0d0,
     '        '>>Define element velocity field first',ERROR,*9999)
            CALL ASSERT(nej_strtrad.NE.0.0d0,
     '        '>>Define element radius field first',ERROR,*9999)
            CALL ASSERT(nej_flow.NE.0.0d0,
     '        '>>Define element flow field first',ERROR,*9999)
            IF(XAB(nej_strtrad,ne).EQ.0.0d0)THEN
              XAB(nej_vel,ne)=0.d0
            ELSE
              XAB(nej_vel,ne)=XAB(nej_flow,ne)/
     '          (PI*(XAB(nej_strtrad,ne)**2.d0))
            ENDIF
          ENDDO !noelem
        ELSE IF(TYPE(1:6).EQ.'LINEAR') THEN
C AJS: Add option to apply linear gradient to field
          IF(CBBREV(CO,'FIELD',3,1,NTCO,N3CO)) THEN
            NJF=IFROMC(CO(N3CO+1))
          ELSE
            NJF=1
          ENDIF
          nj=NJ_LOC(NJL_FIEL,NJF,nr)
          CALL PARSE_NODES(NPNODE,NPLIST,noco,NRLIST,NTCO,CO,ERROR,
     '      *9999)
          CALL ASSERT(NPLIST(0).GT.0,'>>Nodes not found',ERROR,
     '      *9999)
          IF(CBBREV(CO,'MEAN',3,noco+1,NTCO,N3CO))THEN
            mean=RFROMC(CO(N3CO+1))
          ELSE
            CALL ASSERT(.FALSE.,'>>Define mean value.',ERROR,*9999)
          ENDIF
          IF(CBBREV(CO,'GRADIENT',3,noco+1,NTCO,N3CO))THEN
            grad=RFROMC(CO(N3CO+1))
          ELSE
            CALL ASSERT(.FALSE.,'>>Define linear gradient.',ERROR,*9999)
          ENDIF
          IF(CBBREV(CO,'NJ_DIRECTION',3,noco+1,NTCO,N3CO))THEN
            njdir=IFROMC(CO(N3CO+1))
            CALL ASSERT((njdir.EQ.1.OR.njdir.EQ.2.OR.njdir.EQ.3),
     &        '>>nj direction should be 1/2/3.',ERROR,*9999)
          ELSE
            njdir=3 !default = z-direction
          ENDIF

C Set reference height at halfway between min and max in njdir direction
          hmax=XP(1,1,njdir,NPLIST(1))
          hmin=XP(1,1,njdir,NPLIST(1))
          DO nonode=1,NPLIST(0)
            np=NPLIST(nonode)
            IF(XP(1,1,njdir,np).GE.hmax)THEN
              hmax=XP(1,1,njdir,np)
              np_max=np
            ENDIF
            IF(XP(1,1,njdir,np).LE.hmin)THEN
              hmin=XP(1,1,njdir,np)
              np_min=np
            ENDIF
          ENDDO !nonode
          href=0.5d0*(hmax+hmin) 
!           write(*,*) 'nj field = ',nj,' height=',hmax-hmin
C Set field value as field = mean + gradient*(distance from ref height)
          DO nonode=1,NPLIST(0)
            np=NPLIST(nonode)
            ne=NENP(np,1,nr)
            nb=NBJ(nj,ne)
            !nv=NVJE(2,nb,nj,ne)
            DO nv=1,NVJP(nj,np)!NVJE(2,nb,nj,ne)
              XP(1,nv,nj,np)=mean+grad*(XP(1,1,njdir,np)-href)
!               IF(np.EQ.np_max)write(*,*)' At npmax',np,' field=',
!      &          XP(1,nv,nj,np)
!               IF(np.EQ.np_min)write(*,*)' At npmin',np,' field=',
!      &          XP(1,nv,nj,np)
            ENDDO !nv
          ENDDO !nonode
        ELSE IF(TYPE(1:9).EQ.'SUMMATION') THEN
C AJS: Add option to perform summation up tree structure (i.e. flow conservation)
          IF(CBBREV(CO,'FIELD',3,1,NTCO,N3CO)) THEN !specifies field number
            NJF=IFROMC(CO(N3CO+1))
          ELSE
            NJF=1
          ENDIF
          nj=NJ_LOC(NJL_FIEL,NJF,nr)
          CALL SumFlows(NBJ,NEELEM(0,1),nj,NPNE,NXI,NVJE,XP,
     &      ERROR,*9999)
          IF(CBBREV(CO,'NORMALISE',3,1,NTCO,N3CO))THEN
            CALL NORMALISEGENERAL(NBJ,NEELEM(0,1),nj,NPNE,NXI,NVJE,XP,
     &        ERROR,*9999)
          ENDIF !NORMALISE
        ENDIF !type
      ENDIF

      CALL EXITS('UPFIEL')
      RETURN
 9991 IF(XPT_PTR.NE.0) CALL FREE_MEMORY(XPT_PTR,ERROR,*9999)
 9999 CALL ERRORS('UPFIEL',ERROR)
      CALL EXITS('UPFIEL')
      RETURN 1
      END



