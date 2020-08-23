      SUBROUTINE EVREAC(IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,NEELEM,
     '  NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,NKHE,NKJE,NNB,NNF,
     '  NONY,NPLIST,NPF,NPNE,NPNODE,NPNY,
     '  NRE,NRLIST,NSB,NVHE,NVHP,NVJE,
     '  NW,NXI,NXLIST,NYNE,NYNO,NYNP,NYNR,Z_CONT_LIST,CE,CG,CGE,
     '  CONY,
     '  CP,CURVCORRECT,FEXT,FIX,GRR,PG,RE,RG,SE,WG,XA,
     '  XG,XP,YG,YGF,YP,
     '  ZA,ZA1,Z_CONT,ZE,ZP,ZP1,STRING,ERROR,*)

C#### Subroutine: EVREAC
C###  Description:
C###    EVREAC lists finite element node residuals for nonlinear
C###    problems by calling ZPRP (NOTE: uses current material params).
C###    FULL_BCS is used to evaluate reactions with the full set of
C###      specified bcs instead of the current bc increments
C###    COUPLED is used to evaluate reactions for coupled region (nr=0)
C###    VIEW/NOVIEW does/doesn't print residuals to screen

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'coup00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LGE(NHM*NSM),NAN(NIM,NAM,NBFM),
     '  NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NEELEM(0:NE_R_M,0:NRM),NFF(6,NEM),NFFACE(0:NF_R_M,NRM),
     '  NGAP(NIM,NBM),NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),
     '  NKB(2,2,2,NNM,NBFM),NKEF(0:4,16,6,NBFM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),
     '  NONY(0:NOYM,NYM,NRCM,0:NRM,NXM),
     '  NPLIST(0:NPM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM,NXM),
     '  NRE(NEM),NRLIST(0:NRM),NSB(NKM,NNM,NBFM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM),NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NXLIST(0:NXM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),
     '  Z_CONT_LIST(NDM,2,7)
      REAL*8 CE(NMM,NEM,NXM),CG(NMM,NGM),CGE(NMM,NGM,NEM,NXM),
     '  CONY(0:NOYM,NYM,NRCM,0:NRM,NXM),
     '  CP(NMM,NPM,NXM),CURVCORRECT(2,2,NNM,NEM),FEXT(NIFEXTM,NGM,NEM),
     '  GRR(NOM),PG(NSM,NUM,NGM,NBM),
     '  RE(NSM,NHM),RG(NGM),SE(NSM,NBFM,NEM),
     '  WG(NGM,NBM),XA(NAM,NJM,NEM),XG(NJM,NUM),XP(NKM,NVM,NJM,NPM),
     '  YG(NIYGM,NGM,NEM),YGF(NIYGFM,NGFM,NFM),YP(NYM,NIYM,NXM),
     '  ZA(NAM,NHM,NCM,NEM),ZA1(NAM,NHM,NCM,NEM),Z_CONT(NDM,2,67),
     '  ZE(NSM,NHM),
     '  ZP(NKM,NVM,NHM,NPM,NCM),ZP1(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER IN_DIRECTION,GETNYR,i,IBEG,IEND,IFROMC,N1LIST,
     '  N3CO,na,nc,ne,nh,nhx,nk,NH_FIELD,NJF,nj,
     '  no,noelem,nonode,no_nrlist,no_nynr,noy,np,nr,nr_coup,
     '  nv,nx,nxc,ny,ny1,ny_num,RESIDUAL_NUM
      REAL*8 co1,ERRMAX,RATIO(0:6),RFROMC,SUMfixed,
     &  SUMfree,SUMreaction,WEIGHT
      CHARACTER FILE*(MXCH)
      LOGICAL ALL_NODES,ALL_REGIONS,CBBREV,
     '  COUPLED_REAC,DIRECTION,FULL_BCS,
     '  INLIST,OPFILE,RESIDUAL,SINGLE_NUM,VIEW,NJF_CHECK,
     &     NH_CHECK

      CALL ENTERS('EVREAC',*9999)
      
C XSL NEWS 13Aug2010 Changed to be consistant with default value in SOLVE
C Hacked because ERRMAX could be a user input para defined in SOLVE
C and passed into CALC_CONV_RATIO to determine when the simulation should terminate or continue (controlled in NONLIN)
C when the residual is close to ZERO_TOl
C For the purpose of output into the OPREAC file, this value does not play a critial role
      IF(KTYP007.EQ.4) THEN ! scalar product R.delta
        ERRMAX=1.0D-20
      ELSE ! Other convergence criteria
        ERRMAX=1.0D-10
      ENDIF
C XSL NEWE

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM evaluate reactions<;FILENAME>
C###  Parameter:        <number (ID#s/all)[all]>
C###    Specify which reactions to evaluate
C###  Parameter:        <full_bcs>
C###    Used to evaluate reactions with the full set of
C###      specified bcs instead of the current bc increments
C###  Parameter:        <coupled>
C###    Specify a coupled problem
C###  Parameter:        <(view/noview)[view]>
C###    Specify to display the residuals
C###  Parameter:        <nodes (#s/GROUP/all)[all]>
C###    Specify the node numbers to evaulate
C###  Parameter:        <direction (all/1/2/3)[all]>
C###    Specify the direction of the reaction
C###  Parameter:        <weight #[1.0]>
C###    Specify the weight to scale the reaction value
C###  Parameter:        <residual (#/1)[1]>
C###    Specify the residual number to store in 
C###    (for optimisation problems)
C###  Para,eter:        <nh (1/2/3)>
C###    If parameter 'to' is specified, the reacion values
C###    of dependent variable nh are stored in field 'to'.
C###    If parameter 'to' is not specified, nothing happens.
C###  Parameter:        <to (fields #)>
C###    Specifies the field number to store the computed reaction
C###    values of dependent variable 'nh'. Parameter 'to' can only
C###    be specified in conjunction with parameter 'nh'.
C###  Parameter:        <region (#s/all)[1]>
C###    Specify the region numbers to evaluate.
C###  Parameter:        <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    EVREAC lists finite element node residuals for nonlinear
C###      problems

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> '
     '       //'<number (ID#s/all)[all]>'
        OP_STRING(2)=BLANK(1:15)//'<full_bcs>'
        OP_STRING(3)=BLANK(1:15)//'<coupled>'
        OP_STRING(4)=BLANK(1:15)//'<(view/noview)[view]>'
        OP_STRING(5)=BLANK(1:15)//'<nodes (#s/GROUP/all)[all]>'
        OP_STRING(6)=BLANK(1:15)//'<direction (all/1/2/3)[all]>'
        OP_STRING(7)=BLANK(1:15)//'<weight #[1.0]>'
        OP_STRING(9)=BLANK(1:15)//'<residual (#/1)[1]>'
        OP_STRING(10)=BLANK(1:15)//'<nh (1/2/3)>'
        OP_STRING(11)=BLANK(1:15)//'<to (field #)>'
        OP_STRING(12)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(13)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
        
      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','EVREAC',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN
          FILE=COQU(noco,1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          OPFILE=.TRUE.
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opreac','NEW',
     &         'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
          IOFI=IOOP
        ENDIF
        
        CALL ASSERT(USE_NONLIN.GT.0,'>>Set USE_NONLIN to 1 '
     &       //'in parameters file to evaluate reactions',ERROR,*9999)
        
        CALL ASSERT(CALL_SOLV,'>>Solution mapping arrays not set up:'
     &       //' use define solve',ERROR,*9999)
        
        IF(CBBREV(CO,'NUMBER',1,noco+1,NTCO,N3CO)) THEN
          SINGLE_NUM=.TRUE.
          ny_num=IFROMC(CO(N3CO+1))
        ELSE
          SINGLE_NUM=.FALSE.
        ENDIF
        
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     &       ERROR,*9999)
        
        IF(CBBREV(CO,'COUPLED',2,noco+1,NTCO,N3CO)) THEN
          CALL ASSERT(IS_COUPLED(nx),
     &         '>>Define solve for coupled problem',ERROR,*9999)
          CALL ASSERT(.NOT.SINGLE_NUM,'>>Cannot print a single ny '
     &         //'reaction for the coupled region',ERROR,*9999)
          COUPLED_REAC=.TRUE.
          NRLIST(0)=COUP_NRLIST(0,nx)
          DO no_nrlist=1,COUP_NRLIST(0,nx)
            NRLIST(no_nrlist)=COUP_NRLIST(no_nrlist,nx)
          ENDDO
        ELSE
          COUPLED_REAC=.FALSE.
          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,
     &         ERROR,*9999)
        ENDIF
C     CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
C        nxc=NXLIST(1)
C        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
C        CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
C     &    ERROR,*9999)

        IF(CBBREV(CO,'FULL_BCS',1,noco+1,NTCO,N3CO)) THEN
          FULL_BCS=.TRUE.
        ELSE
          FULL_BCS=.FALSE.
        ENDIF
        
        IF(CBBREV(CO,'NOVIEW',3,noco+1,NTCO,N3CO)) THEN
          VIEW=.FALSE.
        ELSE
          VIEW=.TRUE.
        ENDIF
        
        IF(CBBREV(CO,'NODES',3,noco+1,NTCO,N3CO)) THEN
          CALL PARSE_NODES(NPNODE,NPLIST,noco,NRLIST,NTCO,CO,
     &         ERROR,*9999)
          ALL_NODES=.FALSE.
        ELSE                    !list all nodes
          ALL_NODES=.TRUE.
        ENDIF
        
        IF(CBBREV(CO,'DIRECTION',3,noco+1,NTCO,N3CO)) THEN
          IN_DIRECTION=IFROMC(CO(N3CO+1))
          DIRECTION=.TRUE.
        ELSE
          IN_DIRECTION=0
          DIRECTION=.FALSE.
        ENDIF
        
C     new HS/MPN 15/7/05 adding weight scale factor to reaction value
        IF(CBBREV(CO,'WEIGHT',3,noco+1,NTCO,N3CO)) THEN
          WEIGHT=RFROMC(CO(N3CO+1))
        ELSE
          WEIGHT=1.0D0
        ENDIF
        
        IF(CBBREV(CO,'RESIDUAL',3,noco+1,NTCO,N3CO)) THEN
          RESIDUAL_NUM=IFROMC(CO(N3CO+1))
          RESIDUAL=.TRUE.
        ELSE
          RESIDUAL_NUM=0
          RESIDUAL=.FALSE.
        ENDIF
        
C     OR> 10/03/06
C     
C     INCLUDE PARAMETERS TO STORE COMPUTED REACTION VALUES IN A FIELD

        IF(CBBREV(CO,'NH',2,noco+1,NTCO,N3CO)) THEN
          NH_CHECK=.TRUE.
          NH_FIELD=IFROMC(CO(N3CO+1))
        ELSE
          NH_CHECK=.FALSE.
          NH_FIELD=0
        ENDIF
        
        IF(CBBREV(CO,'TO',2,noco+1,NTCO,N3CO)) THEN
          IF(NH_CHECK) THEN
            NJF_CHECK=.TRUE.
            NJF=IFROMC(CO(N3CO+1))
          ELSE
            NJF_CHECK=.FALSE.
            NJF=0
            CALL ASSERT(NH_CHECK,'>>Need to define PARAMETER nh '
     &           //'variable in conjunction with writing '
     &           //'reaction values to a field',ERROR,*9999)
          ENDIF
        ELSE
          NJF_CHECK=.FALSE.
          NJF=0
        ENDIF
        
        
        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          
          CALL ASSERT(nr.GT.0,'>>Evaluate reactions separately '
     &         //'for coupled regions',ERROR,*9999)
          CALL ASSERT(ITYP6(nr,nx).EQ.2,
     &         '>>No nonlinear problem defined',ERROR,*9999)

          CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),NKH(1,1,1,nr),
     &         NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,
     &         YP(1,1,nx),ZA,ZP,ERROR,*9999)
          IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN !cnst Vol
C     Put reference state for cavity from YP(ny,10) into
C     ZA1,ZP1 for ZERE55
            CALL YPZP(10,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     &           NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,
     &           YP(1,1,nx),ZA1,ZP1,ERROR,*9999)
          ENDIF

          IF(FULL_BCS) THEN
C           Apply full set of bcs to ZP and ZA. NOTE: don't change YP
C           as only want incr.s added temporarily for this calculation.
C ***       Add incremental b.c.'s to ZP
            DO nc=1,NCT(nr,nx)
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                DO nhx=1,NHP(np,nr,nx)
                  nh=NH_LOC(nhx,nx)
                  DO nv=1,NVHP(nh,np,nc,nr)
                    DO nk=1,NKH(nh,np,nc,nr)
                      ny=NYNP(nk,nv,nh,np,0,nc,nr)
                      IF(FIX(ny,1,nx).AND.FIX(ny,2,nx)) THEN
                        ZP(nk,nv,nh,np,nc)=ZP(nk,nv,nh,np,nc)
     &                       +YP(ny,2,nx)
                      ENDIF
                    ENDDO       !nk
                  ENDDO         !nv
                ENDDO           !nh
              ENDDO             !np
            ENDDO               !nc
C ***       Add incremental b.c.'s to ZA
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              DO nhx=1,NHE(ne,nx)
                nh=NH_LOC(nhx,nx)
                DO na=1,NAT(NBH(nh,1,ne))
                  ny=NYNE(na,nh,0,1,ne)
                  IF(FIX(ny,1,nx).AND.FIX(ny,2,nx)) THEN
                    ZA(na,nh,1,ne)=ZA(na,nh,1,ne)+YP(ny,2,nx)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDIF

          CALL ZPRP(IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,NEELEM,
     &      NFF,NFFACE,NGAP,NHE(1,nx),NHP(1,nr,nx),
     &      NKB,NKEF,NKH(1,1,1,nr),NKHE,NKJE,NNB,NNF,NPF,NPNE,
     &      NPNODE,NPNY(0,1,0,nx),nr,NRE,NSB,NVHE,NVHP(1,1,1,nr),
     &      NVJE,NW(1,1,nx),nx,NXI,NYNE,NYNP,NYNR(0,0,1,nr,nx),
     &      Z_CONT_LIST,
     &      CE(1,1,nx),CG,CGE(1,1,1,nx),CP(1,1,nx),
     &      CURVCORRECT,FEXT,PG,RE,RG,SE,WG,XA,XG,XP,
     &      YG,YGF,
     &      YP(1,1,nx),ZA,ZA1,%VAL(0),Z_CONT,ZE,%VAL(0),ZP,ZP1,
     &      %VAL(0),FIX,
     &      ERROR,*9999)

          IF(.NOT.COUPLED_REAC.AND..NOT.SINGLE_NUM) THEN
            SUMfree=0.0d0
            SUMfixed=0.0d0
            SUMreaction=0.0d0
            DO no_nynr=1,NYNR(0,0,1,nr,nx) !Loop over global vars
              ny=NYNR(no_nynr,0,1,nr,nx) !is the global var #
              ny1=GETNYR(1,NPNY(0,1,0,nx),nr,1,0,ny,NYNE,NYNP) !row#           

              IF(NPNY(0,ny,0,nx).EQ.1) THEN !nodally based variable
                nk=NPNY(1,ny,0,nx)
                nv=NPNY(2,ny,0,nx)
                nh=NPNY(3,ny,0,nx)
                np=NPNY(4,ny,0,nx)
                
                IF(ALL_NODES.OR.
     &               INLIST(np,NPLIST(1),NPLIST(0),N1LIST)) THEN
                  IF((.NOT.DIRECTION).OR.
     &                 ((nh.EQ.IN_DIRECTION).AND.nk.EQ.1)) THEN
                    
                    IF((nh.eq.1.and.nk.EQ.1).OR.RESIDUAL) THEN
                      IF(FIX(ny,1,nx)) THEN
                        SUMfixed=SUMfixed+DABS(YP(ny1,4,nx))
                        SUMreaction=SUMreaction+YP(ny1,4,nx)
                        IF(VIEW) THEN !print reactions for each region
                          FORMAT='('' np='',I7,'' nh='',I2,'
     &                         //''' nv='',I2,'' nk='',I2,'' ny='',I4,'
     &                         //''' Reaction= '',D12.4,'' (fixed)'')'
                          WRITE(OP_STRING,FORMAT)
     &                         np,nh,nv,nk,ny,YP(ny1,4,nx)
                          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                        ENDIF
                        IF(NJF_CHECK.AND.(NH.EQ.NH_FIELD)) THEN
                          nj=NJ_LOC(NJL_FIEL,NJF,nr)
                          XP(nk,nv,nj,np)=YP(ny1,4,nx)
                        ENDIF 
                      ELSE  
                        SUMfree=SUMfree+DABS(YP(ny1,4,nx)) 
                        SUMreaction=SUMreaction+YP(ny1,4,nx)
                        IF(VIEW) THEN !print reactions for each region
                          FORMAT='('' np='',I7,'' nh='',I2,'
     &                         //''' nv='',I2,'' nk='',I2,'' ny='',I4,'
     &                         //''' Reaction= '',D12.4)'
                          WRITE(OP_STRING,FORMAT)
     &                         np,nh,nv,nk,ny,YP(ny1,4,nx)
                          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                        ENDIF
                        IF(NJF_CHECK.AND.(NH.EQ.NH_FIELD)) THEN
                          nj=NJ_LOC(NJL_FIEL,NJF,nr)
                          XP(nk,nv,nj,np)=YP(ny1,4,nx)
                        ENDIF 
                      ENDIF     !FIX
                    ELSE IF(nh.gt.1.or.nk.GT.1) THEN
                      IF(FIX(ny,1,nx)) THEN
                        SUMfixed=SUMfixed+DABS(YP(ny1,4,nx))
                        IF(VIEW) THEN !print reactions for each region
                          FORMAT='(11X,'' nh='',I2,'' nv='',I2,'
     &                         //''' nk='',I2,'' ny='',I4,'
     &                         //''' Reaction= '',D12.4,'' (fixed)'')'
                          WRITE(OP_STRING,FORMAT)
     &                         nh,nv,nk,ny,YP(ny1,4,nx)
                          CALL WRITES(IOFI,OP_STRING,ERROR,*9999) 
                        ENDIF
                        IF(NJF_CHECK.AND.(NH.EQ.NH_FIELD)) THEN
                          nj=NJ_LOC(NJL_FIEL,NJF,nr)
                          XP(nk,nv,nj,np)=YP(ny1,4,nx)
                        ENDIF 
                      ELSE
                        SUMfree=SUMfree+DABS(YP(ny1,4,nx))
                        IF(VIEW) THEN !print reactions for each region
                          FORMAT='(11X,'' nh='',I2,'' nv='',I2,'
     &                         //''' nk='',I2,'' ny='',I4,'
     &                         //''' Reaction= '',D12.4)'
                          WRITE(OP_STRING,FORMAT)
     &                         nh,nv,nk,ny,YP(ny1,4,nx)
                          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                        ENDIF
                        IF(NJF_CHECK.AND.(NH.EQ.NH_FIELD)) THEN
                          nj=NJ_LOC(NJL_FIEL,NJF,nr)
                          XP(nk,nv,nj,np)=YP(ny1,4,nx)
                        ENDIF 
                      ENDIF     !FIX
                    ENDIF       !nh/nk
                  ENDIF         ! DIRECTION
                ENDIF           ! NPLIST
                
              ELSE IF((NPNY(0,ny,0,nx).EQ.2).AND.
     &               (.NOT.RESIDUAL)) THEN !element based var
                na=NPNY(1,ny,0,nx)
                nh=NPNY(2,ny,0,nx)
                ne=NPNY(4,ny,0,nx)
                IF(na.EQ.1) THEN
                  IF(FIX(ny,1,nx)) THEN
                    SUMfixed=SUMfixed+DABS(YP(ny1,4,nx))
                    IF(VIEW) THEN !print reactions for each region
                      FORMAT='('' ne='',I3,'' nh='',I2,'' na='',I2,'
     &                     //''' ny='',I7,'
     &                     //''' Reaction= '',D12.4,'' (fixed)'')'
                      WRITE(OP_STRING,FORMAT) ne,nh,na,ny,
     &                     YP(ny1,4,nx)
                      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                    ENDIF
                    IF(NJF_CHECK) THEN
C     OR> DOES ELEMENT BASED VARIABLES MAKE SENSE
                    ENDIF 
                  ELSE
                    SUMfree=SUMfree+DABS(YP(ny1,4,nx))
                    IF(VIEW) THEN !print reactions for each region
                      FORMAT='('' ne='',I3,'' nh='',I2,'' na='',I2,'
     &                     //''' ny='',I7,'' Reaction= '',D12.4)'
                      WRITE(OP_STRING,FORMAT) ne,nh,na,ny,
     &                     YP(ny1,4,nx)
                      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                    ENDIF
                    IF(NJF_CHECK) THEN
C     OR> DOES ELEMENT BASED VARIABLES MAKE SENSE
                    ENDIF 
                  ENDIF         !FIX
                ELSE IF(na.GT.1) THEN
                  IF(FIX(ny,1,nx)) THEN
                    SUMfixed=SUMfixed+DABS(YP(ny1,4,nx))
                    IF(VIEW) THEN !print reactions for each region
                      FORMAT='(7X,'' nh='',I2,'' na='',I2,'' ny='','
     &                     //'I7,'' Reaction= '',D12.4,'' (fixed)'')'
                      WRITE(OP_STRING,FORMAT) nh,na,ny,YP(ny1,4,nx)
                      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                    ENDIF
                    IF(NJF_CHECK) THEN
C     OR> DOES ELEMENT BASED VARIABLES MAKE SENSE
                    ENDIF 
                  ELSE
                    SUMfree=SUMfree+DABS(YP(ny1,4,nx))
                    IF(VIEW) THEN !print reactions for each region
                      FORMAT='(7X,'' nh='',I2,'' na='',I2,'
     &                     //''' ny='',I7,'' Reaction= '',D12.4)'
                      WRITE(OP_STRING,FORMAT) nh,na,ny,YP(ny1,4,nx)
                      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                    ENDIF
                    IF(NJF_CHECK) THEN
C     OR> DOES ELEMENT BASED VARIABLES MAKE SENSE
                    ENDIF 
                  ENDIF         !FIX
                ENDIF           !na.EQ.1
              ENDIF             !NPNY - node or elem based var
            ENDDO               !no_nynr (ny)
            WRITE(OP_STRING,'('' Sum of absolute constrained '
     &           //'reactions: '',D12.4)') SUMfixed
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' Sum of absolute unconstrained '
     &           //'reactions: '',D12.4)') SUMfree
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            IF (RESIDUAL) THEN
              WRITE(OP_STRING,'('' Sum of residual '
     &             //'reactions: '',D12.4)') SUMreaction
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              YP(RESIDUAL_NUM,7,1)=WEIGHT*SUMreaction                  
            ENDIF
C*** 08/08/07 XSL Calculate and output convergence ratio
C 25/11/08 XSL Hacked, local var (in NONLIN) CONTACT was set to false
C to provide the correct logic when evaluating energy norm.
C Note: The residual evaluated here will be different from the terminal output
C for energy norm because YP5 has been updated,
C similarly unconstrained R will be different for contact problem as ZPRP was 
C performed with updated deformed coord.
            CALL CALC_CONV_RATIO(1,IOFI,NBH,ne,NEELEM,
     &        NONY(0,1,1,0,nx),NPNY(0,1,0,nx),nr,nr,nx,
     &        NYNO(0,1,1,0,nx),NYNR(0,0,1,0,nx),
     &        CONY(0,1,1,0,nx),ERRMAX,GRR,RATIO,YG,YP(1,1,nx),.FALSE.,
     &        .TRUE.,
     &        ERROR,*9999)
          ENDIF                 !.NOT.COUPLED_REAC.AND..NOT.SINGLE_NUM
        ENDDO                   !nrlist
        
        IF(SINGLE_NUM) THEN
          FORMAT='('' Reaction at ny ='',I7,'' is '',D12.4)'
          WRITE(OP_STRING,FORMAT) ny_num,YP(ny_num,4,nx)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE IF(COUPLED_REAC) THEN !print react.s for coupled reg (nr=0)
          nr_coup=0             !global coupled region #
C         initialise GRR (temporary storage for reactions)
          DO no=1,NOT(1,1,nr_coup,nx)
            GRR(no)=0.0d0
          ENDDO                 !no
          DO no_nynr=1,NYNR(0,1,1,nr_coup,nx) !loop rows of glob reg 0
            ny1=NYNR(no_nynr,1,1,nr_coup,nx) !is row number
            DO noy=1,NONY(0,ny1,1,nr_coup,nx) !soln rows assc with ny1
              no=NONY(noy,ny1,1,nr_coup,nx) !is row number for ny1
              co1=CONY(noy,ny1,1,nr_coup,nx) !is coupling coeff for ny1
              GRR(no)=GRR(no)+YP(ny1,4,nx)*co1
            ENDDO               !noy
          ENDDO                 !no_nynr (ny)
          IF(VIEW) THEN
            WRITE(OP_STRING,'(/'' Reactions associated with free '
     &           //'variables of the global region (nr=0):'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            SUMfree=0.0d0
            DO no=1,NOT(1,1,nr_coup,nx)
              SUMfree=SUMfree+DABS(GRR(no))
C!!!          Note only printing indices for FIRST ny coupled to no
              ny=NYNO(1,no,1,nr_coup,nx)
              IF(NPNY(0,ny,1,nx).EQ.1) THEN !nodal variable
                WRITE(OP_STRING,'('' Reaction = '','
     &               //'D12.4,''  at ny='',I5,'' nk='',I2,'
     &               //''' nv='',I2,'' nh='',I2,'' np='',I7,'
     &               //''' nc='',I2,'' nr='',I2)')
     &               GRR(no),ny,(NPNY(i,ny,1,nx),i=1,6)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ELSE IF(NPNY(0,ny,1,nx).EQ.2) THEN !element var
                WRITE(OP_STRING,'('' Reaction = '','
     &               //'D12.4,''  at ny='',I5,'' na='',I2,'' nh='',I2,'
     &               //''' nc='',I2,'' ne='',I4,'' nr='',I2)')
     &               GRR(no),ny,(NPNY(i,ny,1,nx),i=1,5)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO               !no
            WRITE(OP_STRING,'('' Sum of absolute unconstrained '
     &           //'reactions: '',D12.4)') SUMfree
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C*** 08/08/07 XSL Calculate and output convergence ratio
C 25/11/08 XSL Hacked, local var (in NONLIN) CONTACT was set to false
C to provide the correct logic when evaluating energy norm.
C Note: The residual evaluated here will be different from the terminal output
C for energy norm because YP5 has been updated,
C similarly unconstrained R will be different for contact problem as ZPRP was 
C performed with updated deformed coord.
            CALL CALC_CONV_RATIO(1,IOFI,NBH,ne,NEELEM,
     &        NONY(0,1,1,0,nx),NPNY(0,1,0,nx),nr,nr_coup,nx,
     &        NYNO(0,1,1,0,nx),NYNR(0,0,1,0,nx),
     &        CONY(0,1,1,0,nx),ERRMAX,GRR,RATIO,YG,YP(1,1,nx),.FALSE.,
     &        .TRUE.,
     &        ERROR,*9999)
          ENDIF                 !VIEW
        ENDIF                   !COUPLED_REAC
        
        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF
      
      CALL EXITS('EVREAC')
      RETURN
 9999 CALL ERRORS('EVREAC',ERROR)
      CALL EXITS('EVREAC')
      RETURN 1
      END



