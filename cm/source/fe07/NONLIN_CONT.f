      SUBROUTINE NONLIN_CONT(IBT,IDO,INP,ISC_GK,
     '  ISC_GKK,ISC_GQ,ISR_GK,ISR_GKK,ISR_GQ,LGE,NAN,NBH,
     '  NBHF,NBJ,NBJF,NEELEM,NENP,NFF,NFFACE,NGAP,
     '  NHE,NHP,NKB,NKEF,NKH,NKHE,NKJE,NMNO,NNB,NNF,NONY,
     '  NORD,NPF,NP_INTERFACE,NPNE,NPNODE,NPNY,nr_solve,
     '  NRE,NRLIST,NSB,NVHE,NVHP,NVJE,NW,nx,NXI,NYNE,NYNO,
     '  NYNP,NYNR,Z_CONT_LIST,CE,CG,CGE,CONY,CP,
     '  CURVCORRECT,CYNO,ERRMAX,FEXT,
     '  GK,GKK,GQ,GR,GRR,PG,RE1,RG,SE,WG,XA,
     '  XG,XO,XP,YG,YGF,YP,ZA,ZA1,Z_CONT,ZE,ZE1,ZP,ZP1,
     '  FIX,ERROR,*)

C#### Subroutine: NONLIN_CONT
C###  Description:
C###    NONLIN_CONT controls the solution of contact mechanics
C###    problems using an augmented Newton-Raphson procedure. 
C###    Bodies involved may be described through linear or 
C###    finite elasticity.

C****  YP(ny,1)  current equilibrium solution - T+dt
C****  YP(ny,2)  prescribed dep var/force increms set by FIX(ny,1)
C****  YP(ny,3)  prescribed initial equilibrium solution
C****  YP(ny,4)  current set of equilibrium equation residuals
C****  YP(ny,5)  temporary storage 1 used to store solution
C****            increments to add to YP(ny,1)
C****  YP(ny,6)  storage of contact surface pressure
C****  YP(ny,8)  previous solution at iteration i-1 
C****  YP(ny,11) velocity at iteration i

C***   CONT_IT - Iteration counter
C***   LOAD_IT - Load step counter
C***   AUG_IT - Augmentation counter
C***   CONV_IT - Convergence check counter

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'acti01.cmn'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'disp00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'fsklib.inc'
      INCLUDE 'gen000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'gks000.cmn'
      INCLUDE 'host00.cmn'
      INCLUDE 'host00.inc'
      INCLUDE 'ipma50.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'load00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'ktyp70.cmn'
      INCLUDE 'nonl00.cmn'
      INCLUDE 'ofst00.cmn'
      INCLUDE 'press00.cmn'
      INCLUDE 'time02.cmn'
      INCLUDE 'ptr00.cmn'
      INCLUDE 'solv00.cmn'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),ISC_GK(NISC_GKM),
     '  ISC_GKK(NISC_GKKM),ISC_GQ(NISC_GQM),ISR_GK(NISR_GKM),
     '  ISR_GKK(NISR_GKKM),ISR_GQ(NISR_GQM),LGE(NHM*NSM,NRCM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),
     '  NBJ(NJM,NEM),NBJF(NJM,NFM),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),NFF(6,NEM),
     '  NFFACE(0:NF_R_M,NRM),NGAP(NIM,NBM),NHE(NEM),
     '  NHP(NPM,0:NRM),NKB(2,2,2,NNM,NBFM),NKEF(0:4,16,6,NBFM),
     '  NKH(NHM,NPM,NCM,0:NRM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),
     '  NONY(0:NOYM,NYM,NRCM,0:NRM),NORD(5,NE_R_M),NMNO(1:2,0:NOPM),
     '  NPF(9,NFM),NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM),nr_solve,
     '  NRE(NEM),NRLIST(0:NRM),NSB(NKM,NNM,NBFM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx,
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNO(0:NYOM,NOOPM,NRCM,0:NRM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM),
     '  Z_CONT_LIST(NDM,2,7)
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),CGE(NMM,NGM,NEM),
     '  CONY(0:NOYM,NYM,NRCM,0:NRM),CP(NMM,NPM),
     '  CURVCORRECT(2,2,NNM,NEM),CYNO(0:NYOM,NOOPM,NRCM,0:NRM),
     '  ERRMAX,
     '  FEXT(NIFEXTM,NGM,NEM),GK(NZ_GK_M),GKK(NZ_GKK_M),
     '  GQ(NZ_GQ_M),GR(NYROWM),GRR(NOM),PG(NSM,NUM,NGM,NBM),
     '  RE1(NSM,NHM),RG(NGM),SE(NSM,NBFM,NEM),
     '  WG(NGM,NBM),XA(NAM,NJM,NEM),
     '  XG(NJM,NUM),XO(NOM),XP(NKM,NVM,NJM,NPM),YG(NIYGM,NGM,NEM),
     '  YGF(NIYGFM,NGFM,NFM),YP(NYM,NIYM),ZA(NAM,NHM,NCM,NEM),
     '  ZA1(NAM,NHM,NCM,NEM),Z_CONT(NDM,2,67),ZE(NSM,NHM),
     '  ZE1(NSM,NHM),
     '  ZP(NKM,NVM,NHM,NPM,NCM),ZP1(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM),ERROR_FLAG,GQ_ASSEM

!     Local Variables
      INTEGER ERR,GETNYR,IP,nd,j,nc,no,noload,nonrlist,no_nynr,
     '  noy,np,nonode,nhx,nj,nh,nv,nk,nr,nr1,NWRIT,ny,ny1,ny2,nyo,
     '  METHOD

      REAL*8 co,CONT_PTS_CONT,CONT_PTS_TIED,GAP_CHECK_CONT,
     '  GAP_CHECK_TIED,GAP_SUM_CONT,GAP_SUM_TIED,RATIO,
     '  RSUM_CONSTRAINED,RSUM_SOLINCR,RSUM_UNCONSTRAIN,ZERO_TOL,
     '  ZAA(NAM,NHM,NCM,NEM),ZPA(NKM,NVM,NHM,NPM,NCM)

      REAL ELAPSED_TIME,TIME_START1(1),TIME_STOP(1) 

C ***  ZPA,ZAA store acceleration information 

      LOGICAL CONVERGED,OUTPUT,REITER,SOLVE,UPDATE_MATRIX

      DATA ZERO_TOL /1.0d-8/

      CALL ENTERS('NONLIN_CONT',*9999)

      CALL ASSERT(NIYM.GE.5,
     '    '>>Increase NIYM to 5 for general NONLIN problems',
     '    ERROR,*9999)

      CALL ASSERT(NIYFIXM.GE.1,
     '    '>>Increase NIYFIXM to 1 for general NONLIN problems',
     '    ERROR,*9999)

C *** NJM must be 67 or greater to store contact point field values.
C *** in WD(NJM,nd) and ZD(NJM,nd) updated in UPDATA (fe23)
      CALL ASSERT(NJM.GE.67,
     '    '>>Increase NJM to 67 for general CONTACT problems',
     '    ERROR,*9999)

      DO nonrlist=1,NRLIST(0)
        nr=NRLIST(nonrlist)
        IF(ITYP1(nr,nx).EQ.5) THEN
          CALL ASSERT(USE_NONLIN.GT.0,'>>Cannot solve NONLIN problems '
     '      //'as USE_NONLIN is zero in parameters.inc',ERROR,*9999)
          CALL ASSERT(NIFEXTM.GE.7,
     '      '>>Increase NIFEXTM to 7 for finite deformation problems',
     '      ERROR,*9999)
        ENDIF
      ENDDO !nr

      TIME=TSTART
      DTIME=DT
                
      nr1=NRLIST(1)
      METHOD=ITYP9(nr1,nx)
C *** Check that solution methods are the same in all regions
      DO nonrlist=2,NRLIST(0)
        nr=NRLIST(nonrlist)
        CALL ASSERT(ITYP9(nr,nx).EQ.METHOD,
     '    '>>Iteration method not same in all regions',ERROR,*9999)
        CALL ASSERT(ITYP5(nr,nx).EQ.ITYP5(nr1,nx),
     '    '>>Problem type not same in all regions',ERROR,*9999)
      ENDDO

C *** Loop over all regions in list
      SOLVE=.FALSE.
      DO nonrlist=1,NRLIST(0)
        nr=NRLIST(nonrlist)
C ***   Solve equations (if at least one dof is not fixed)
        DO no_nynr=1,NYNR(0,0,1,nr) !loop over global variables
          ny=NYNR(no_nynr,0,1,nr) !is global variable number
          IF(.NOT.FIX(ny,1)) SOLVE=.TRUE.
        ENDDO !no_nynr
      ENDDO !nonrlist (nr)

      IF(NTLOAD.EQ.0) THEN
        REITER=.TRUE.
        NTLOAD=1
      ELSE
        REITER=.FALSE.
      ENDIF

      OUTPUT=.FALSE.
      NWRIT=1

      DO noload=1,NTLOAD
      
        NOSTEP=noload !to pass via /NONL/ (nonl00.cmn) to other routines
        RATIO=0.0d0
        IF(IWRIT1(nr_solve,nx).GT.0) THEN
          IF(noload.EQ.1.OR.IWRIT2(nr_solve,nx).EQ.2.
     '      OR.NWRIT.EQ.IWRIT1(nr_solve,nx)) THEN
            OUTPUT=.TRUE. !output initial residual
            NWRIT=1
          ELSE
            OUTPUT=.FALSE.
            NWRIT=NWRIT+1
          ENDIF
        ENDIF

        DO nonrlist=1,NRLIST(0)

C ***     Only add loads once
          nr=NRLIST(nonrlist)
                    
          IF(.NOT.REITER) THEN
                        
            IF(nonrlist.EQ.1) THEN
              WRITE(OP_STRING,'(/'' Load step'',I3/,1X,12(''=''))')
     '          noload
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF !nonrlist
                    
C ***       Setup previous solution and velocity information at t=0.0
            IF (CONT_IT.EQ.1) THEN
              DO nc=1,NCT(nr,nx) !loop over RHS(displ) and LHS(force) vars
                DO no_nynr=1,NYNR(0,0,nc,nr) !loop over global variables
                  ny=NYNR(no_nynr,0,nc,nr) !is global variable number
                  IF(NPNY(0,ny,0).EQ.1) THEN
                    np=NPNY(4,ny,0)
                    CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                  ENDIF
                  IF (KTYP5I(nr).EQ.1) THEN ! inertia
                    YP(ny,12)=0.0d0 ! Previous velocity=0.0 at time=0.0
                    YP(ny,8)=YP(ny,1) ! Previous solution at time=0.0
                  ENDIF 
                ENDDO !no_nynr (ny)
              ENDDO !nc
            ENDIF !CONT_IT
   
            IF (LOAD_IT.EQ.1) THEN
C ***         Check which points are in contact
C ***         Only do once at beginning of each load step
              DO nd=1,NDT
C ***           Check signed normal gap Z_CONT(nd,j,13) for slave (j=1)
                IF (Z_CONT_LIST(nd,1,4).NE.2) THEN ! contact
                  IF (Z_CONT(nd,1,13).GE.GAP_TOL) THEN ! then in contact
                    Z_CONT_LIST(nd,1,3)=1
                    Z_CONT_LIST(nd,2,3)=1
                  ELSE ! not in contact
C ***             Check if problem starts with active contraints
                    IF (KTYP5J(nr).GE.1) THEN ! active
                      Z_CONT_LIST(nd,1,3)=1
                      Z_CONT_LIST(nd,2,3)=1                                                                                  
                    ELSE ! not active                                   
                      Z_CONT_LIST(nd,1,3)=0
                      Z_CONT_LIST(nd,2,3)=0                                                                                    
                    ENDIF ! KTYP5J
                  ENDIF
                ELSE ! tied contact -always in contact
                  Z_CONT_LIST(nd,1,3)=1
                  Z_CONT_LIST(nd,2,3)=1   
                ENDIF 
              ENDDO !nd 
              KTYP5J(nr)=0 ! turn of initial active constraint                                           
            ENDIF !LOAD_IT
  
            IF (LOAD_IT.EQ.1) THEN           
              DO nc=1,NCT(nr,nx) !loop over RHS(displ) and LHS(force) vars
                DO no_nynr=1,NYNR(0,0,nc,nr) !loop over global variables
                  ny=NYNR(no_nynr,0,nc,nr) !is global variable number
                  IF(NPNY(0,ny,0).EQ.1) THEN
                    np=NPNY(4,ny,0)
                    CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                  ENDIF
C ***             Apply increments to displacement and force b.c.s
                  IF(FIX(ny,1)) THEN !ny has incremented ess bdry cond
                    YP(ny,1)=YP(ny,1)+YP(ny,2)*FACTOR
                  ENDIF  
                  IF (KTYP5I(nr).EQ.1) THEN ! inertia
                    YP(ny,11)=(YP(ny,1)-YP(ny,8))/T_inc ! velocity 
                    YP(ny,14)=(YP(ny,11)-YP(ny,12))/T_inc ! acceleration 
                  ENDIF
                  YP(ny,8)=YP(ny,1) ! Update prescribed boundary conditions  
                ENDDO !no_nynr
              ENDDO !nc
            ENDIF !LOAD_IT

            IF (CONT_IT.EQ.1) THEN
              IF (KTYP5I(nr).EQ.1) THEN ! inertia
C ***           Get ZPA,ZAA terms from YP(ny,11) : velocity 
                CALL YPZP(11,NBH,NEELEM,NHE,NHP(1,nr),NKH(1,1,1,nr),
     '          NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZAA,ZPA,
     '          ERROR,*9999)
              ENDIF !KTYP5I
            ENDIF !CONT_IT

            IF (CONT_IT.EQ.1) THEN
              IF(ITYP2(1,nx).EQ.1) THEN !linear elasticity
C ***           For contact problems initially set ZP to XP
C ***           Only do this once (ie first iteration)
                nc=1 ! dependent variable
C ***           force (nc=2) RHS variables already present in ZP at this point.
                DO nonode=1,NPNODE(0,nr)
                  np=NPNODE(nonode,nr)
                  DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
                    nj=NJ_LOC(NJL_GEOM,nhx,nr)
                    nh=NH_LOC(nhx,nx)
                    DO nv=1,NVHP(nh,np,nc,nr)
                      DO nk=1,NKH(nh,np,nc,nr)
                        ZP(nk,nv,nh,np,nc)=XP(nk,nv,nj,np)
                      ENDDO !nk
                    ENDDO !nv
                  ENDDO !nh
                ENDDO !nonode (np)
C ***           Copy ZP into YP(ny,3) initial conditions
C ***           Only do this at begining
                CALL ZPYP(3,NBH,NEELEM,NHE,NHP,NKH,
     '            NPNODE,nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
C ***           Copy ZP into YP(ny,1) current solution
                CALL ZPYP(1,NBH,NEELEM,NHE,NHP,NKH,
     '            NPNODE,nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
              ENDIF ! ITYP2
            ENDIF !CONT_IT
                   
          ENDIF !REITER
                  
C ***     Get residuals
          IF(ITYP2(1,nx).EQ.1) THEN !linear elasticity
                
C ***       Initialise residuals
            DO no_nynr=1,NYNR(0,1,1,nr) !loop over rows
              ny=NYNR(no_nynr,1,1,nr) !is row number
              YP(ny,4)=0.0d0
            ENDDO

C ***       Add in global loads
            DO no_nynr=1,NYNR(0,1,1,nr) !loop over rows
              ny=NYNR(no_nynr,1,1,nr) !row number
              IF(NPNY(0,ny,0).EQ.1) THEN
                np=NPNY(4,ny,0)
                CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
              ENDIF
              ny1=GETNYR(2,NPNY,nr,0,1,ny,NYNE,NYNP) !is RHS var #
              YP(ny,4)=YP(ny,4)-YP(ny1,1)
            ENDDO !no_nynr

          ELSE !finite elasticity

C ***       Get residual for finite elasticity
            CALL YPZP(1,NBH,NEELEM,NHE,NHP(1,nr),NKH(1,1,1,nr),NPNODE,
     '        nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)

            IF (KTYP5I(nr).EQ.1) THEN ! inertia
C ***         Get ZPA,ZAA terms from YP(ny,11) : velocity 
              CALL YPZP(11,NBH,NEELEM,NHE,NHP(1,nr),NKH(1,1,1,nr),
     '          NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZAA,ZPA,
     '          ERROR,*9999)
            ENDIF

            CALL ZPRP(IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,NEELEM,
     '        NFF,NFFACE,NGAP,NHE,NHP(1,nr),NKB,NKEF,NKH(1,1,1,nr),
     '        NKHE,NKJE,NNB,NNF,NPF,NPNE,NPNODE,NPNY,nr,NRE,NSB,
     '        NVHE,NVHP(1,1,1,nr),NVJE,NW,nx,NXI,NYNE,NYNP,
     '        NYNR(0,0,1,nr),Z_CONT_LIST,CE,CG,CGE,CP,
     '        CURVCORRECT,FEXT,PG,RE1,RG,
     '        SE,WG,XA,XG,XP,YG,YGF,YP,ZA,ZA1,ZAA,Z_CONT,
     '        ZE,ZE1,ZP,ZP1,ZPA,FIX,
     '        ERROR,*9999)

          ENDIF !ITYP2   

          ERROR_FLAG=.FALSE.
                  
C  ***    Initialise YP(ny,6) - contact forces
          DO no_nynr=1,NYNR(0,0,1,nr) !loop over global vars
            ny=NYNR(no_nynr,0,1,nr) !is global variable number
            IF(NPNY(0,ny,0).EQ.1) THEN
              np=NPNY(4,ny,0)
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
            ENDIF
            YP(ny,6)=0.0d0      
          ENDDO !no_nynr

C ***     Modify global residual vector YP(ny,4) with contact contribution
          IF(.NOT.ERROR_FLAG) THEN
            CALL CPU_TIMER(CPU_USER,TIME_START1)
            CALL CONTACT_RESIDUAL(IBT,IDO,INP,NBH,NBJ,NBJF,NFF,NHE,
     '        NKEF,NKHE,NKJE,NNF,NPF,NPNE,nr,NRE,NVHE,NVJE,NW,nx,NYNP,
     '        Z_CONT_LIST,CURVCORRECT,SE,YP,XA,XP,Z_CONT,ZA,ZE,ZP,FIX,
     '        ERROR,*300)
            GO TO 301
C ***       This statement is designed to be skipped if no error
C ***       occurs.  However if a error occurs within a subroutine
C ***       the alternate return points to line 300 to set the flag

 300        CONTINUE
            ERROR_FLAG=.TRUE.
            WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
            CALL WRITES(IOER,OP_STRING,ERROR,*301)
 301        CONTINUE
          ENDIF !.NOT.ERROR_FLAG
          CALL ASSERT(.NOT.ERROR_FLAG,'>>An error occurred during '
     '      //'contact residual calculations',ERROR,*9999)
          CALL CPU_TIMER(CPU_USER,TIME_STOP)
          ELAPSED_TIME=TIME_STOP(1)-TIME_START1(1)
          IF(IWRIT4(nr,nx).GE.1) THEN
            WRITE(OP_STRING,
     '      '('' CPU time of 1 thread for global residual contact '
     '       //'calcs: '',D11.4,'' s'')') ELAPSED_TIME
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

C ***     Output residual stats
          IF(OUTPUT.AND.IWRIT3(nr,nx).GT.0) THEN
            IF(REITER) THEN
              IP=1
            ELSE
              IP=2
            ENDIF
            CALL ZPOP(IP,NBH,1,NEELEM,NHP(1,nr),NKH(1,1,1,nr),NPNODE,
     '        nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,4),ZA,ZP,FIX(1,1),
     '        ERROR,*9999)
          ENDIF
                 
        ENDDO !nonrlist (nr)

        RSUM_CONSTRAINED=0.0D0
        RSUM_UNCONSTRAIN=0.0D0

C ***   Calc sum of constrained and unconstrained residuals
C ***   using GRR as temp storage for constraint reduced eqn resids
        DO no=1,NOT(1,1,nr_solve,nx) !loop over global soln rows
          GRR(no)=0.0d0
        ENDDO !no

        DO no_nynr=1,NYNR(0,1,1,nr_solve) !loop over rows
          ny1=NYNR(no_nynr,1,1,nr_solve) !is row number
          IF(NONY(0,ny1,1,nr_solve).GT.0) THEN !free dependent variable
            ny2=GETNYR(2,NPNY,nr_solve,0,1,ny1,NYNE,NYNP) !is RHS var #
            DO noy=1,NONY(0,ny1,1,nr_solve) !loop on rows assoc with ny1
              no=NONY(noy,ny1,1,nr_solve) !is row number for ny1
              co=CONY(noy,ny1,1,nr_solve) !is coupling coeff for ny1
              GRR(no)=GRR(no)+(YP(ny1,4)-YP(ny2,1))*co
            ENDDO !noy
          ELSE !bdry cond applied to dependent variable
            RSUM_CONSTRAINED=RSUM_CONSTRAINED+DABS(YP(ny1,4))
          ENDIF
        ENDDO !no_nynr

        DO no=1,NOT(1,1,nr_solve,nx) !loop over global soln rows
          DO nyo=1,NYNO(0,no,1,nr_solve)
            ny1=NYNO(nyo,no,1,nr_solve) !is row number
            ny2=GETNYR(2,NPNY,nr_solve,0,1,ny1,NYNE,NYNP) !is RHS var #
C ***       Only add in residual if no force bc is applied
            IF(.NOT.FIX(ny2,1))
     '        RSUM_UNCONSTRAIN=RSUM_UNCONSTRAIN+DABS(GRR(no))
          ENDDO !nyo
        ENDDO !no
            
C ***   Calc sum of absolute solution vector increments
        RSUM_SOLINCR=0.0d0
        IF (CONT_IT.EQ.1) THEN
          RSUM_SOLINCR=0.0d0 ! first iteration
        ELSE
          DO no=1,NOT(2,1,nr_solve,nx) !loop over global soln variables
            IF(NYNO(0,no,2,nr_solve).GT.0) THEN
              ny=NYNO(1,no,2,nr_solve) !is first global variable number
              RSUM_SOLINCR=RSUM_SOLINCR+DABS(YP(ny,5)) !coupled to no
            ENDIF !NYNO
          ENDDO !no
        ENDIF !CONT_IT

        IF(RSUM_CONSTRAINED.LT.ZERO_TOL) THEN
          IF(RSUM_UNCONSTRAIN.LT.ZERO_TOL) THEN
            RATIO=ZERO_TOL
          ELSE
            RATIO=1.0d0/ZERO_TOL
          ENDIF
        ELSE
          RATIO=RSUM_UNCONSTRAIN/RSUM_CONSTRAINED
        ENDIF
        WRITE(OP_STRING,
     '    '(/'' Sum of unconstrained residuals='',D11.4,'
     '    //'/'' Sum of constrained residuals  ='',D11.4,'
     '    //'/'' Their ratio (unconstr/constr) ='',D11.4,'
     '    //'/'' Sum of soln vector increments ='',D11.4)')
     '    RSUM_UNCONSTRAIN,RSUM_CONSTRAINED,RATIO,RSUM_SOLINCR
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

        CONVERGED=.FALSE.
       IF((RSUM_UNCONSTRAIN.LE.ERRMAX).OR.(CONV_IT.GT.NTITER)) THEN
c        IF(RSUM_UNCONSTRAIN.LE.ERRMAX) THEN
          IF (CONV_IT.EQ.1) THEN
            CONVERGED=.FALSE.
            CALL SET_USER_INTEGER('CONVERGED',0,ERR)
          ELSE
            CONVERGED=.TRUE.
            CALL SET_USER_INTEGER('CONVERGED',1,ERR)
          ENDIF
        ENDIF

        IF(CONVERGED) THEN

C          WRITE(OP_STRING,'(/'' Convergence achieved after '','
C     '      //'I2,'' iterations'')') LOAD_IT

C NEWS JHC 10-03-05 outputs iterations to 3-significant figures     
          WRITE(OP_STRING,'(/'' Convergence achieved after '','
     '      //'I3,'' iterations'')') LOAD_IT
C NEWE
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C ***     Check if gap convergence is achieved.
          GAP_SUM_TIED=0.0d0
          CONT_PTS_TIED=0.0d0
          GAP_SUM_CONT=0.0d0
          CONT_PTS_CONT=0.0d0
                     
          DO nd=1,NDT
C ***       Check signed normal gap Z_CONT(nd,j,13) for slave (j=1)
            IF (Z_CONT_LIST(nd,1,3).EQ.1) THEN ! contact point
              IF (Z_CONT_LIST(nd,1,4).NE.2) THEN ! contact (normal or friction)
                IF (Z_CONT(nd,1,13).GE.0.0d0) THEN ! then in contact
                  GAP_SUM_CONT=GAP_SUM_CONT+(Z_CONT(nd,1,13)**2)
                  CONT_PTS_CONT=CONT_PTS_CONT+1.0d0
                ENDIF      
              ELSE ! tied contact point                 
                GAP_SUM_TIED=GAP_SUM_TIED+(Z_CONT(nd,1,13)**2)
                CONT_PTS_TIED=CONT_PTS_TIED+1.0d0
              ENDIF
            ENDIF
          ENDDO !nd

          IF (CONT_PTS_CONT.EQ.0.0d0) THEN
            GAP_CHECK_CONT=0.0d0
          ELSE
            GAP_CHECK_CONT=DSQRT(GAP_SUM_CONT/CONT_PTS_CONT)
          ENDIF
          IF (CONT_PTS_TIED.EQ.0.0d0) THEN
            GAP_CHECK_TIED=0.0d0
          ELSE
            GAP_CHECK_TIED=DSQRT(GAP_SUM_TIED/CONT_PTS_TIED)
          ENDIF
    
          WRITE(OP_STRING,'(/'' Contact Gap '',D11.4,
     '      '''')') GAP_CHECK_CONT
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
  
          WRITE(OP_STRING,'(/'' Tied Gap '',D11.4,
     '      '''')') GAP_CHECK_TIED
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            
          AUG_IT=AUG_IT+1
          IF (AUG_IT.LE.AUGMENT) THEN
          
C ***       Reset CONVERGED status
            CONVERGED=.FALSE.
            CALL SET_USER_INTEGER('CONVERGED',0,ERR)
            WRITE(OP_STRING,'(/'' Augmentation '',I2,
     '        '''')') AUG_IT
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
    
C ***       Perform augmentation
            DO j=1,2
              DO nd=1,NDT
                Z_CONT(nd,j,17)=Z_CONT(nd,j,18) ! Update normal contact force estimate
                Z_CONT(nd,j,19)=Z_CONT(nd,j,20) ! Update tangential1 contact force estimate
                Z_CONT(nd,j,21)=Z_CONT(nd,j,22) ! Update tangential2 contact force estimate
              ENDDO !nd
            ENDDO !j
    
C ***       Set YP(ny,1) back to previous solution
            DO nonrlist=1,NRLIST(0) ! update initial guess
              nr=NRLIST(nonrlist)
              DO no_nynr=1,NYNR(0,0,1,nr) !loop over global vars
                ny=NYNR(no_nynr,0,1,nr) !is global variable number
                IF(NPNY(0,ny,0).EQ.1) THEN
                  np=NPNY(4,ny,0)
                  CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                ENDIF
                YP(ny,1)=YP(ny,8) ! restore previous solution 
                YP(ny,11)=YP(ny,12) ! restore previous velocity 
                IF (KTYP5I(nr).EQ.1) THEN ! inertia
                  YP(ny,11)=(YP(ny,1)-YP(ny,8))/T_inc ! velocity
                  YP(ny,14)=(YP(ny,11)-YP(ny,12))/T_inc ! acceleration 
                ENDIF            
              ENDDO !no_nynr
            ENDDO !nonrlist (nr)  
            CONV_IT=1 ! convergence counter         
            GOTO 9990

          ENDIF !AUG_IT

          IF (KTYP5I(nr).EQ.1) THEN ! inertia
            DO nonrlist=1,NRLIST(0) ! update initial guess
              nr=NRLIST(nonrlist)
              DO no_nynr=1,NYNR(0,0,1,nr) !loop over global vars
                ny=NYNR(no_nynr,0,1,nr) !is global variable number
                IF(NPNY(0,ny,0).EQ.1) THEN
                  np=NPNY(4,ny,0)
                  CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                ENDIF
                YP(ny,8)=YP(ny,1) ! store previous solution
                YP(ny,12)=YP(ny,11) ! store previous velocity     
              ENDDO !no_nynr
            ENDDO !nonrlist (nr)
          ENDIF

          CONV_IT=1 ! convergence counter
          LOAD_IT=1 ! reset load step counter
          AUG_IT=0 ! reset augmentation count
          
          GOTO 9990
 
        ELSE

          WRITE(OP_STRING,
     '      '(/'' Convergence has not been reached after '','
     '      //'I3,'' iterations'')') LOAD_IT
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CALL SET_USER_INTEGER('CONVERGED',0,ERR)

        ENDIF !Converged

        IF(NTITER.EQ.0) GOTO 9990

        IF(.NOT.SOLVE) THEN
          WRITE(OP_STRING,'(/'' no free variables'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          GO TO 9998
        ENDIF

C ***   Assemble stiffness matrix
        IF(ITYP2(1,nx).EQ.1) THEN !linear elasticity

          DO nonrlist=1,NRLIST(0)
            nr=NRLIST(nonrlist)
            IF(IWRIT4(nr,nx).GE.1) THEN
              WRITE(OP_STRING,'(/'' Region '',I1,'':'')') nr
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            GQ_ASSEM=.FALSE.
            TIME=TSTART
            GQ_ASSEM=.TRUE.
            CALL ASSEMBLE1(IBT,IDO,INP,ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,
     '        ISR_GKK,ISR_GQ,NBH,NBJ,NEELEM,NHE,NHP(1,nr),NKJE,
     '        NONY(0,1,1,nr),NORD,NPF,NP_INTERFACE,NPNE,NPNY,nr,
     '        nr_solve,NRE,NVHE,NVJE,NW,nx,NYNE,NYNP,NYNR(0,0,1,nr),CE,
     '        CGE,CONY(0,1,1,nr),CP,CURVCORRECT,GK,GKK,GQ,GR,PG,SE,WG,
     '        XA,XP,YG,GQ_ASSEM,.TRUE.,.TRUE.,ERROR,*9999)
            ASSEMBLE_GLOBAL(nr,nx)=.TRUE.
          ENDDO !nonrlist

        ELSE !finite elasticiy

C ***     Assemble for finite elasticy
          CALL ASSEMBLE5(IBT,IDO,INP,ISC_GK,ISR_GK,NAN,NBH,NBHF,
     '      NBJ,NBJF,NEELEM,NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,
     '      NKH,NKHE,NKJE,NMNO,NNB,NNF,NPF,NPNE,NPNODE,NPNY,NRE,
     '      NRLIST,nr_solve,NSB,NVHE,NVHP,NVJE,NW,nx,NXI,NYNE,NYNP,
     '      NYNR,CE,CGE,CP,CURVCORRECT,
     '      FEXT,GK,PG,SE,WG,XA,XP,YG,YGF,YP,ZA,ZA1,
     '      ZAA,ZP,ZP1,ZPA,FIX,ERROR,*9999)

        ENDIF !ITYP2

C ***   Modify global stiffness matrix GK with contact contribution
        IF(.NOT.ERROR_FLAG) THEN
          CALL CPU_TIMER(CPU_USER,TIME_START1)
          CALL CONTACT_STIFFNESS(IBT,IDO,INP,ISC_GK,ISR_GK,NBH,NBHF,
     '      NFF,NHE,NKEF,NKHE,NNF,NPF,NPNE,NRE,NVHE,nx,NYNP,
     '      Z_CONT_LIST,GK,SE,Z_CONT,FIX,
     '      ERROR,*400)
          GO TO 401
C ***     This statement is designed to be skipped if no error
C ***     occurs.  However if a error occurs within a subroutine
C ***     the alternate return points to line 400 to set the flag

 400      CONTINUE
          ERROR_FLAG=.TRUE.
          WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
          CALL WRITES(IOER,OP_STRING,ERROR,*401)
 401      CONTINUE
        ENDIF !.NOT.ERROR_FLAG
        CALL ASSERT(.NOT.ERROR_FLAG,'>>An error occurred during '
     '    //'contact stiffness calculations',ERROR,*9999)
        CALL CPU_TIMER(CPU_USER,TIME_STOP)
          ELAPSED_TIME=TIME_STOP(1)-TIME_START1(1)
        IF(IWRIT4(nr,nx).GE.1) THEN
          WRITE(OP_STRING,
     '    '(/'' For contact stiffness calcs: CPU time of 1 thread:'','
     '    //'D11.4,'' s'')') ELAPSED_TIME
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

C ***   New matrix
        FIRSTS(nx)=.TRUE.
        UPDATE_MATRIX=.TRUE.

C ***   Solve global system of equations
        CALL SOLVE5(ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,ISR_GKK,ISR_GQ,
     '    LGE,NBH,NENP,NHE,NONY(0,1,1,nr_solve),NP_INTERFACE,NPNE,
     '    NPNY,nr_solve,NRE,NVHE,nx,NYNE,
     '    NYNO(0,1,1,nr_solve),NYNP,NYNR(0,0,1,nr_solve),
     '    CONY(0,1,1,nr_solve),CYNO(0,1,1,nr_solve),GK,GKK,GQ,
     '    GR,GRR,XO,YP,FIRSTS(nx),UPDATE_MATRIX,ERROR,*9999)

        DO nonrlist=1,NRLIST(0) ! update initial guess
          nr=NRLIST(nonrlist)
C ***     Update current solution by adding increments
          DO no_nynr=1,NYNR(0,0,1,nr) !loop over global vars
            ny=NYNR(no_nynr,0,1,nr) !is global variable number
            IF(NPNY(0,ny,0).EQ.1) THEN
              np=NPNY(4,ny,0)
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
            ENDIF
            YP(ny,1)=YP(ny,1)+YP(ny,5)
            IF (KTYP5I(nr).EQ.1) THEN ! inertia
              YP(ny,11)=(YP(ny,1)-YP(ny,8))/T_inc ! velocity
              YP(ny,14)=(YP(ny,11)-YP(ny,12))/T_inc ! acceleration
            ENDIF
          ENDDO !no_nynr
        ENDDO !nonrlist (nr)

        IF (KTYP5I(nr).EQ.1) THEN ! inertia
C ***     Copy YP(ny,11) (current solution) into ZAA,ZPA
          DO nonrlist=1,NRLIST(0)
            nr=NRLIST(nonrlist)
            CALL YPZP(11,NBH,NEELEM,NHE,NHP(1,nr),NKH(1,1,1,nr),
     '        NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZAA,ZPA,
     '        ERROR,*9999)
          ENDDO
        ENDIF

        IF(ITYP2(1,nx).EQ.1) THEN !linear elasticity
C ***      Copy YP(ny,1) (current solution) into ZP
           DO nonrlist=1,NRLIST(0)
             nr=NRLIST(nonrlist)
             CALL YPZP(1,NBH,NEELEM,NHE,NHP(1,nr),NKH(1,1,1,nr),
     '         NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZA,ZP,
     '         ERROR,*9999)
           ENDDO
        ENDIF ! ITYP2        

        CONT_IT=CONT_IT+1
        LOAD_IT=LOAD_IT+1
        CONV_IT=CONV_IT+1

 9990   CONTINUE

      ENDDO !noload (load step)

 9998 CALL EXITS('NONLIN_CONT')
      RETURN

 9999 CALL ERRORS('NONLIN_CONT',ERROR)
      CALL EXITS('NONLIN_CONT')
      RETURN 1
      END


