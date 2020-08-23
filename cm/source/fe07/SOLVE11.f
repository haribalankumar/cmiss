      SUBROUTINE SOLVE11(IBT,IDO,INP,ISC_GD,ISR_GD,ISC_GK,ISC_GKK,
     &  ISR_GK,ISR_GKK,nb,NBH,NBJ,NEELEM,NELIST,NENP,NHE,NHST,NKJE,
     &  NO_NE,NONY,NORD,NPF,NPLIST2,NPNE,NPNODE,NPNY,nr,NVJE,NVJP,nx,
     &  NXI,NYNE,NYNO,NYNP,NYNR,NY_OFST,NZ_ESED,SPARSE_S,BBM,CE,CG,CGE,
     &  CONY,CYNO,CP,ED,EM,ER,ES,GD,GK,GKK,GK2,GR,GRR,PG,
     &  RG,SE,STACK_ED,STACK_EM,STACK_ES,WG,XA,XAB,XG,XO,XP,
     &  YG,YP,ZE,ZG,DYNAM1,FIRST_A,FIX,LINEAR,SMOOTHING,UPDATE_MATRIX,
     &  UPDATE_VECTOR,UPDATE_VELOCITY,FIRST_NX,RET_ERROR,*)

C#### Subroutine: SOLVE11
C###  Description:
C###    SOLVE11 calls ASSEMBLE11 to construct the unreduced system of
C###    matrices for pulmonary transport problems.  Boundary conditions
C###    are applied to reduce the system, and SOLVE-SYSTEM is called to
C###    produce an initial solution.  For iterative problems (coupled
C###    water-heat, coupled pressure-flow in capillaries) the system is
C###    updated and the solution is recalculated until convergence is
C###    reached.
      
      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b08.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'error0.inc'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'loc00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'pulm00.cmn'
      INCLUDE 'ptr00.cmn'
      INCLUDE 'solv00.cmn'
      INCLUDE 'time02.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISC_GD(NISC_GDM),ISR_GD(NISR_GDM),ISC_GK(NISC_GKM),
     '  ISC_GKK(NISC_GKKM),ISR_GK(NISR_GKM),ISR_GKK(NISR_GKKM),nb,
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NEELEM(0:NE_R_M),NELIST(0:NEM),
     '  NENP(NPM,0:NEPM),NHE(NEM),NHST,NKJE(NKM,NNM,NJM,NEM),
     &  NO_NE(NEM),NORD(5,NE_R_M),NPF(9,NFM),NPLIST2(0:NPM),
     &  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M),NPNY(0:6,NYM,0:NRCM),
     '  nr,NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),nx,
     &  NXI(-NIM:NIM,0:NEIM,0:NEM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     &  NYNO(0:NYOM,NOOPM,NRCM),NY_OFST(NYM),NONY(0:NOYM,NYM,NRCM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM),NYNR(0:NY_R_M,0:NRCM,NCM),
     '  NZ_ESED(0:12,NE_R_M),SPARSE_S
      REAL*8 BBM(2,NEM),CE(NMM,NEM),CG(NMM,NGM),CGE(NMM,NGM,NEM),
     '  CONY(0:NOYM,NYM,NRCM),CYNO(0:NYOM,NOOPM,NRCM),
     '  CP(NMM,NPM),FLOW,GD(NZ_GD_M),GK(NZ_GK_M),GKK(NZ_GKK_M),
     &  GK2(NZ_GK_M),GR(NYROWM),GRR(NOM),PG(NSM,NUM,NGM,NBM),
     &  SE(NSM,NBFM,NEM),STACK_ED(12,NE_R_M),STACK_EM(12,NE_R_M),
     &  STACK_ES(12,NE_R_M),WG(NGM,NBM),XA(NAM,NJM,NEM),XO(NOM),
     &  XP(NKM,NVM,NJM,NPM),YG(NIYGM,NGM,NEM),YP(NYM,NIYM),ZE(NSM,NHM),
     &  ZG(NHM,NUM),ED(NHM*NSM,NHM*NSM),EM(NHM*NSM,NHM*NSM),ER(NHM*NSM),
     &  ES(NHM*NSM,NHM*NSM),RG(NGM),XAB(NORM,NEM),XG(NJM,NUM)
      LOGICAL DYNAM1,FIRST_NX,FIRST_A,FIX(NYM,NIYFIXM),
     &  LINEAR,SMOOTHING,UPDATE_MATRIX,UPDATE_VECTOR,
     &  UPDATE_VELOCITY
      CHARACTER RET_ERROR*(*)
!     Local Variables
      INTEGER i,KOUNT,ne,ne0,N_ERR,N_MIN_ERR,nh,nhx,nhs_cap,nj,
     &  noelem,noelem2,nonode,no_nynr,no_nynr1,no_nynr2,no1,no2,noy1,np,
     &  np1,np2,nyo2,ny,ny1,ny2,ny3,nz,nzr,nzz,ost1,ost2
c KOUNT_NE(NEM),
      REAL*8 AA,BB,co1,co2,ERR,ERR2,ERR_RAD,HEIGHT(NJT),L_in,L_out,
     &  LPM_FUNC,MIN_ERR,P1,P2,Ppl,Q01,R_in,R_out,stretch,x_cap,y_cap,
     &  z_cap,k_factor

      REAL TIME_START,TIME_STOP

      LOGICAL CONVERGED,FIRST_IT,UPDATE_LOCAL
      CHARACTER CHAR*1,ERROR*(ERRSTRLEN),STRING*255
!     External function
      INTEGER GETNYR_NONR,IDIGITS

      CALL ENTERS('SOLVE11',*9999)
      UPDATE_LOCAL=.TRUE.
      FIRST_IT=.TRUE. !first iteration
      MIN_ERR=1.d10
      N_MIN_ERR=0
      nzr=0
      IF(ITYP3(nr,nx).LE.2.OR.ITYP3(nr,nx).EQ.5)THEN
        A1=DT*THETA(1)
      ELSE
        A1=1.d0 !capillary, currently no time-dependence
      ENDIF !ITYP3
c      IF(ITYP3(nr,nx).EQ.3.OR.ITYP3(nr,nx).EQ.6) THEN
c        DO noelem=1,NEM
c          KOUNT_NE(noelem)=0 !records # iterations a vessel has zero diameter
c        ENDDO
c      ENDIF
      
      IF(FIRST_NX)THEN

C***  Store the previous solution
        IF(7.GT.NIYM) THEN
          WRITE(CHAR,'(I1)') IDIGITS(7)
          WRITE(ERROR,'(''>>Increase NIYM to '',I'//CHAR//')') 7
          GO TO 9999
        ENDIF
        
        DO no_nynr=1,NYNR(0,0,1) !Loop over global vars of GK
          ny=NYNR(no_nynr,0,1) !global variable #
          YP(ny,2)=YP(ny,1)
          YP(ny,4)=0.d0
          YP(ny,7)=0.d0
        ENDDO !no_nynr
        CALL CALC_SPARSE_GKK_1DTREE(NISC_GKKM,NISR_GKKM,ISC_GK,ISC_GKK,
     &    ISR_GK,ISR_GKK,NYT(1,1,nx),NYT(2,1,nx),NBJ,NDIAG(nx),NEELEM,
     &    NENP,NHST,NPLIST2,NPNE,nr,nx,NXI,NYNE,NYNP,NYNR,NZZT(1,nr,nx),
     &    1,FIX,ERROR,*9999)
      
        FIRST_NX=.FALSE.
        FIRST_INSP=.FALSE.
      ENDIF !FIRST

      IF(NYNR(NYNR(0,0,1),0,1).GT.NOM) THEN
        WRITE(CHAR,'(I1)') IDIGITS(NYNR(NYNR(0,0,1),0,1))
        WRITE(ERROR,'(''>>Increase NOM to '',I'//CHAR//')')
     &    NYNR(NYNR(0,0,1),0,1)
        GO TO 9999
      ENDIF
      DO no_nynr=1,NYNR(0,0,1) !Loop over global vars of GK
        ny=NYNR(no_nynr,0,1) !global variable #
        YP(ny,3)=YP(ny,1) !storage of solution at kth iteration
        YP(ny,6)=YP(ny,2) !storage of wall solution at kth iteration
        XO(ny)=YP(ny,1)
      ENDDO !no_nynr
      CONVERGED=.FALSE.
      KOUNT=0
      DO WHILE(.NOT.CONVERGED)
        KOUNT=KOUNT+1
C***  Calculate the sparsity pattern for the reduced matrices
        IF(UPDATE_MATRIX)THEN
          DO nzz=1,NZZT(1,nr,nx)
            GKK(nzz)=0.0d0
          ENDDO !nzz
          IF(ITYP3(nr,nx).LE.2.OR.(ITYP3(nr,nx).GT.2.
     &      AND.KOUNT.EQ.1).OR.ITYP3(nr,nx).EQ.5)THEN
            CALL ASSEMBLE11(IBT,IDO,INP,ISC_GD,ISC_GK,ISR_GD,ISR_GK,nb,
     &        NBH,NBJ,NEELEM,NENP,NHE,NKJE,NO_NE,NORD,NPF,NPLIST2,NPNE,
     &        nr,NVJE,nx,NYNE,NYNP,NYNR,NZ_ESED,nzr,CE,CG,CGE,CP,ED,EM,
     '        ER,ES,GD,GK,GK2,GR,PG,RG,SE,STACK_ED,STACK_EM,
     &        STACK_ES,WG,XA,XAB,XG,XP,YG,ZE,ZG,UPDATE_MATRIX,
     &        UPDATE_VECTOR,UPDATE_VELOCITY,ERROR,*9999)
          ELSE !Poiselle flow, update resistance values only
            DO noelem=1,NEELEM(0) !update for all ne
              nz=NZ_ESED(1,noelem)
              ne=NEELEM(noelem)
              IF(ITYP3(nr,nx).EQ.3.OR.ITYP3(nr,nx).EQ.6)THEN !Kamm capillary blood flow
                GK(nz)=-CE(nm_Rseg,ne) !contribution from Q(1-2)
              ELSEIF(ITYP3(nr,nx).EQ.4)THEN !Simple P-R-airflow solution
                IF(COUPLE_VIA_LPM.EQ.'Y'.OR.
     &           COMPLIANCE.EQ.3.OR.COMPLIANCE.EQ.6) THEN
                   GK(nz)=-XAB(nej_resis,ne)
                ELSE
                   GK(nz)=XAB(nej_resis,ne)
                ENDIF
              ENDIF
            ENDDO !noelem
          ENDIF
          IF(NOT(2,1,nr,nx).EQ.0) THEN
            ERROR=' >>The number of unknowns is zero'
            GOTO 9999
          ENDIF !not
          CALL ASSERT(DYNAM1,'>>DYNAM1 is not set',ERROR,*9999)
          IF(KOUNT.EQ.1)THEN !only set up the first time through
            ost1=0
            DO no_nynr1=1,NYNR(0,1,1) !Loop over global rows of GK
              ny1=NYNR(no_nynr1,1,1) !is row #
              IF(FIX(ny1,1))THEN !new line
                ost1=ost1+1 !count # of fixed BC
              ENDIF
              NY_OFST(ny1)=ost1
            ENDDO !no_nynr1
          ENDIF
        ELSE
          IF(ITYP3(nr,nx).LE.2)THEN !quick updating for inert gas/wht
            DO nzz=1,NZZT(1,nr,nx)
              GKK(nzz)=0.0d0
            ENDDO !nzz
            DO nz=1,NZT(1,nx)
              GK(nz)=0.d0
              GK2(nz)=0.d0
              GD(nz)=0.d0
            ENDDO !nz
            DO noelem2=1,NELIST(0) !only for ne with changing geometry
              ne=NELIST(noelem2)
              noelem=NO_NE(ne)
              nhs_cap=0
              np=NPNE(2,1,ne)
              CALL ASSEMBLE11_DYNAM(IBT,IDO,INP,ISC_GK,ISR_GK,nb,
     &          NBH(1,1,ne),NBJ(1,ne),ne,NENP,NHE(ne),nhs_cap,
     &          NKJE(1,1,1,ne),NORD,NPF(1,1),NPLIST2,NPNE,nr,
     &          NVJE(1,1,1,ne),nx,NYNE,NYNP,NZ_ESED(0,noelem),nzr,
     '          CE(1,ne),CG,CGE(1,1,ne),CP,ED,EM,ER,ES,GK,GR,PG,RG,
     &          SE(1,1,ne),STACK_ED(1,noelem),STACK_EM(1,noelem),
     &          STACK_ES(1,noelem),WG,XA,XAB,XG,XP,YG(1,1,ne),ZE,ZG,
     &          .TRUE.,UPDATE_VECTOR,ERROR,*9999)
            ENDDO !noelem 
            DO noelem=1,NEELEM(0) !update for all ne
              DO i=1,NZ_ESED(0,noelem)
                nz=NZ_ESED(i,noelem)
                GK(nz)=GK(nz)+STACK_ES(i,noelem)
                GD(nz)=GD(nz)+STACK_ED(i,noelem)
                GK2(nz)=GK2(nz)+STACK_EM(i,noelem)
              ENDDO !i
            ENDDO !noelem
          ENDIF
        ENDIF !UPDATE_MATRIX
        
        DO no1=1,NOT(1,1,nr,nx)
          GRR(no1)=0.0d0
        ENDDO !no1
C***  Apply any flux boundary conditions
        IF(ITYP3(nr,nx).EQ.1)THEN
          DO no_nynr1=1,NYNR(0,0,2) !loop over global vars for nc=2
            ny1=NYNR(no_nynr1,0,2) !is global (flux) variable number
            IF(FIX(ny1,1)) THEN !flux is set as a b.c.
              ny2=GETNYR_NONR(2,NPNY,1,0,ny1,NYNE,NYNP) !flux row #
              DO noy1=1,NONY(0,ny2,1)
                no1=NONY(noy1,ny2,1)
                co1=CONY(noy1,ny2,1)
                GRR(no1)=YP(ny1,3)*co1*DT !GRR plus appl flux bc
              ENDDO !noy1
            ENDIF
          ENDDO !no_nynr1 (ny1)
        ENDIF !ITYP3(nr,nx).EQ.1
        
C***  Set a local logical for updating the solution matrix. i.e. for
C***  inert gas mixing or water and heat transfer the updating is
C***  always done, but for flow models it is only doen if UPDATE_MATRIX
C***  is .true.
c        IF(UPDATE_MATRIX)THEN
        UPDATE_LOCAL=.TRUE.
C***  Assemble the reduced system of matrices
        DO no_nynr1=1,NYNR(0,1,1) !Loop over global rows of GK
          ny1=NYNR(no_nynr1,1,1) !is row #
          IF(.NOT.FIX(ny1,1).OR.(ITYP3(nr,nx).EQ.3.OR.ITYP3(nr,nx)
     &      .EQ.4.OR.ITYP3(nr,nx).EQ.6))THEN
            ost1=0
            no1=NONY(1,ny1,1) !is no# attached to row ny1
            BB=GR(ny1) !get reduced R.H.S.vector
            DO no_nynr2=ISR_GK(ny1),ISR_GK(ny1+1)-1 !each row entry
              ny2=ISC_GK(no_nynr2) !ny2 for column #
              ny3=ny2
c                ny3=GETNYR_NONR(1,NPNY,2,0,ny2,NYNE,NYNP) !local GK var #
c                CALL SPARSE(ny1,ny3,NYT(1,1,nx),nz,NZ_GK_M,NZT(1,nx),
c     '            ISC_GK,ISR_GK,6,ERROR,*100)
              CALL SPARSE(ny1,ny3,NYT(1,1,nx),nz,NZ_GK_M,NZT(1,nx),
     '          ISC_GK,ISR_GK,6,ERROR,*9999)
              IF(ITYP3(nr,nx).EQ.1)THEN !inert gas mixing
                BB=BB-(GK(nz)+GK2(nz))*YP(ny2,3) !-K.c^n
              ELSE IF(ITYP3(nr,nx).EQ.2)THEN !thermofluid dynamics
                BB=BB-GK(nz)*YP(ny2,3)+GK2(nz)*YP(ny2,2)
              ELSEIF(ITYP3(nr,nx).EQ.5)THEN ! MCT
                BB=BB-(GK(nz)+GK2(nz))*YP(ny2,3) ! R - K.c^n
              ELSEIF(ITYP3(nr,nx).GE.3.AND.FIX(ny3,1)) THEN !simple flow
                BB=BB-GK(nz)*YP(ny3,3)
              ENDIF !ITYP3
              IF(UPDATE_LOCAL)THEN
                IF(ITYP3(nr,nx).EQ.1) THEN !inert gas mixing
                  AA=GD(nz)+A1*(GK(nz)+GK2(nz)) !M+K*dt*theta
                ELSE IF(ITYP3(nr,nx).EQ.2)THEN !thermofluid dynamics
                  AA=GD(nz)+A1*GK(nz) !M+K*dt*theta
                ELSEIF(ITYP3(nr,nx).EQ.5)THEN ! MCT
                  AA=GD(nz)+A1*(GK(nz)+GK2(nz)) ! M + K*dt*theta
                ELSEIF(ITYP3(nr,nx).GE.3) THEN !capillaries
                  AA=A1*GK(nz) !K*dt*theta
                ENDIF
                IF(.NOT.FIX(ny2,1))THEN
                  IF(ITYP3(nr,nx).EQ.3.OR.ITYP3(nr,nx).EQ.4.OR.
     &              ITYP3(nr,nx).EQ.6)THEN
C... for capillaries only remove columns, not rows, therefore ost1=0
                    ost1=0
                    ost2=NY_OFST(ny2)
                  ELSE
                    ost1=NY_OFST(ny1)
                    ost2=NY_OFST(ny2)  
                  ENDIF
                  CALL SPARSE(ny1-ost1,ny3-ost2,NYT(1,1,nx),nzz,
     '              NZ_GKK_M,NZZT(1,nr,nx),ISC_GKK,ISR_GKK,SPARSE_S,
     '              ERROR,*9999)
                  WRITE(STRING,'(''>>nzz=0 for ny1,ny2 = '',I5,I5)')
     &              ny1,ny2
                  CALL ASSERT(nzz.GT.0,STRING,ERROR,*9999)
                  GKK(nzz)=GKK(nzz)+AA
                ENDIF !.NOT.FIX
              ENDIF !UPDATE_LOCAL
            ENDDO !no_nynr2
            CALL ASSERT(no1.LE.NOM,'>>Increase NOM',ERROR,*9999)
            GRR(no1)=GRR(no1)+BB !cmht
          ENDIF !FIX
        ENDDO !no_nynr1

c        IF(ITYP3(nr,nx).EQ.1)THEN
c          CALL REPLACE_FIXED_CONCS(NBJ,NEELEM,NPNE,NVJE,nx,
c     '      NXI,NYNP,XP,YP(1,1),FIX,ERROR,*9999)
c        ENDIF

C***  Solve reduced system
        CALL SOLVE_SYSTEM(ISC_GKK,ISR_GKK,NDIAG(nx),NDIAG(nx),NDIAG(nx),
     &    NZZT(1,nr,nx),IWRIT4(nr,nx),PRECON_CODE(1),SOLVEROPTION(1),
     &    SPARSEGKK(1),GKK,GRR,XO,FIRST_A,.TRUE.,.FALSE.,nx,
     &    ERROR,*9999)
        ERR=0.d0
        N_ERR=0
C***  Put solution values into YP
        DO no2=1,NOT(2,1,nr,nx)
          DO nyo2=1,NYNO(0,no2,2)
            ny2=NYNO(nyo2,no2,2)
            co2=CYNO(nyo2,no2,2)
            IF(ITYP3(nr,nx).LE.2.OR.ITYP3(nr,nx).EQ.5)THEN
              IF(DABS(YP(ny2,4)).GT.ZERO_TOL)THEN
                ERR=ERR+(DT*(XO(no2)*co2-YP(ny2,4))/
     '            (YP(ny2,3)+DT*YP(ny2,4)))**2.d0
              ENDIF
              IF(.NOT.FIX(ny2,1))THEN
                YP(ny2,4)=XO(no2)*co2 !dc/dt
                YP(ny2,1)=YP(ny2,3)+DT*XO(no2)*co2 !c^(n+1)=c^(n)+dc
             ENDIF
             IF(DABS(YP(ny2,1)).LT.ZERO_TOL) YP(ny2,1)=ZERO_TOL              
             IF(DABS(YP(ny2,1)).GT.1.d0) YP(ny2,1)=1.d0          

            ELSE IF(ITYP3(nr,nx).GE.3.and.ITYP3(nr,nx).ne.4)THEN
              YP(ny2,5)=YP(ny2,1) !temp storage of previous solution
              YP(ny2,1)=XO(no2)*co2 !new pressure & flow solutions
              IF(DABS(YP(ny2,1)).GT.ZERO_TOL)THEN
                YP(ny2,4)=(YP(ny2,5)-YP(ny2,1))**2.d0/YP(ny2,1)**2
                ERR=ERR+YP(ny2,4) !((YP(ny2,5)-YP(ny2,1))/YP(ny2,1))**2
              ENDIF
            ELSEIF(ITYP3(nr,nx).EQ.4)THEN
              YP(ny2,5)=YP(ny2,1) !temp storage of previous solution
              YP(ny2,1)=XO(no2)*co2 !new pressure & flow solutions
              !IF(KOUNT.EQ.2)THEN
              !YP(ny2,1)=(YP(ny2,1)+YP(ny2,5))/2.d0 !FIRST IT ALWAYS OVERESTIMATE THE PVR SO TRYING TO REIN IT IN A BIT 
                !write(*,*) YP(ny2,1),YP(ny2,5)
                !Use WRITE(OP_STRING...) and CALL WRITES for this - see above (in comments)
              !ENDIF
              IF(DABS(YP(ny2,1)).GT.ZERO_TOL)THEN
                YP(ny2,4)=(YP(ny2,5)-YP(ny2,1))**2.d0/YP(ny2,1)**2
                ERR=ERR+YP(ny2,4) !((YP(ny2,5)-YP(ny2,1))/YP(ny2,1))**2
              ENDIF
            ENDIF !ITYP3
          ENDDO !nyo2
        ENDDO !no2

        IF(SMOOTHING)THEN
          IF(FLOW.LT.0.d0)THEN
            DO nhx=1,NHST
              nh=NH_LOC(nhx,nx)
              CALL SMOOTH_EXPIRATION(NBJ,NEELEM,nh,NORD,NPNE,NVJE,NXI,
     '          NYNP,XP,YP(1,1),ERROR,*9999)
            ENDDO
          ENDIF
        ENDIF
        IF(N_ERR.NE.0) ERR=ERR/N_ERR
C***  Iterate for non-linear systems
        IF(.NOT.LINEAR)THEN
          IF(ITYP3(nr,nx).EQ.2)THEN
            CALL RADIAL_SOL(NBJ,NENP,NPNE,NPNODE,nr,nx,NXI,NYNP,CE,
     &        CP,ERR_RAD,XP,YP,ERROR,*9999)
c            IF(ERR.LE.CONVERG_TOL*1.d3.AND.(.NOT.FIRST_IT))THEN
            IF(ERR.LE.1.d-5.AND.(.NOT.FIRST_IT))THEN
              CONVERGED=.TRUE.
              DO nonode=1,NPNODE(0) !put radial solution into CP
                np=NPNODE(nonode)
                DO i=2,7
                  CP(i,np)=CP(i+6,np)
                ENDDO !i
                IF(ITYP12(nr,nx).EQ.2)THEN
                  DO i=1,17
                    XP(1,1,nj_ode+i,np)=XP(1,2,nj_ode+i,np)
                  ENDDO
                    
                ENDIF
              ENDDO !nonode
            ELSE !if error not converged
              IF(ERR.GE.MIN_ERR)THEN
                N_MIN_ERR=N_MIN_ERR+1
              ELSE
                MIN_ERR=MIN(MIN_ERR,ERR)
              ENDIF
              UPDATE_MATRIX=.FALSE.
              UPDATE_VECTOR=.FALSE.
            ENDIF !ERR
          ELSE IF(ITYP3(nr,nx).EQ.3.OR.ITYP3(nr,nx).EQ.6)THEN !capillary blood flow
C          ELSE IF(ITYP3(nr,nx).EQ.3)THEN !capillary blood flow
            ERR2=0.d0 !initialise
c            DO noelem=1,NEELEM(0) !This is actually done in CALC_HEMATOCRIT
c              ne=NEELEM(noelem) !element #
c              ny=NYNE(1,1,1,1,ne)!ny#,assumes na=1 for flow,make general
c              CE(nm_flow,ne)=YP(ny,1) !stores flow solution for export
c            ENDDO !noelem
            IF(KOUNT.EQ.1) THEN
C... to determine inlet & outlet vessels (only do 1st time through)
              INLETS(0)=0
              OUTLETS(0)=0
              DO noelem=1,NEELEM(0)
                ne=NEELEM(noelem)
                IF(NENP(NPNE(1,nb,ne),0).EQ.1.OR.
     '            NENP(NPNE(2,nb,ne),0).EQ.1) THEN !feed vessel
                  IF(NENP(NPNE(1,nb,ne),0).EQ.1) THEN
                    np1=NPNE(1,nb,ne)
                    np2=NPNE(2,nb,ne)
                    P1=YP(NYNP(1,1,1,np1,0,1),1) !find out direction of flow
                    P2=YP(NYNP(1,1,1,np2,0,1),1)
                  ELSE IF(NENP(NPNE(2,nb,ne),0).EQ.1) THEN
                    np1=NPNE(2,nb,ne)
                    np2=NPNE(1,nb,ne)
                    P1=YP(NYNP(1,1,1,np1,0,1),1)
                    P2=YP(NYNP(1,1,1,np2,0,1),1)
                  ENDIF
                  IF(P2.LT.P1) THEN !inlet vessel
                    INLETS(0)=INLETS(0)+1
                    INLETS(INLETS(0))=ne
                  ELSE IF(P1.LT.P2) THEN !outlet vessel
                    OUTLETS(0)=OUTLETS(0)+1
                    CALL ASSERT(OUTLETS(0).LT.200000,
     &                '>>outlet array not'//
     &                ' large enough, must increase in pulm00.cmn',
     &                ERROR,*9999)
                    OUTLETS(OUTLETS(0))=ne
                  ENDIF
                ENDIF !feed vessel
              ENDDO !noelem
            ENDIF !KOUNT.EQ.1
C... calculates hematocrit in each element
            CALL CALC_HEMATOCRIT(nb,NEELEM,NENP,NPNE,NORD,NYNE,CE,
     &        ERR2,YP,ERROR,*9999)
C... CAP_DIMENSION includes diameter model for alveolar capillary
C... vessels, and for larger arterioles/venules.
C... For alveolar vessels - it calculates: a0 & C0, from this - a*, b*,
C... Dh. For arterioles/venules: diameter(D)=D0+compliance*Pressure
c            CALL CAP_DIMENSION(ITYP3(nr,nx),KOUNT,KOUNT_NE,nb,NEELEM,
c     &        NENP,NPNE,NVJE,NYNP,CE,XP,YP,ERROR,*9999)
C... recalculates Refd & Rfactor using a*, b* approximations
            CALL CALC_HEMODYNAMICS(ITYP3(nr,nx),nb,NEELEM,NENP,NPNE,
     &        NVJE,CE,XP,ERROR,*9999)              
            ERR=ERR/NYNR(0,1,1)   !sum of the error divided by the
            ERR2=ERR2/NYNR(0,1,1) !total number of variables
C!!! why convergence errors so large? hematocrit error largest, b/c some
C!!! segments have Hd changing to 0.d0, if flow too small.
C!!! solution error (err) important for convergence not err2 (Hd error)
           IF(ERR.LT.1.0d-6) THEN !soln converged
              CONVERGED=.TRUE.
            ELSE !if error not converged
              IF(ERR.GE.MIN_ERR) THEN !.OR.ERR2.GE.MIN_ERR)
                N_MIN_ERR=N_MIN_ERR+1
              ELSE
                MIN_ERR=ERR
              ENDIF
            ENDIF !ERR
          ELSE IF(ITYP3(nr,nx).EQ.4) THEN
            IF (COMPLIANCE.ne.2)THEN
              CALL CALC_PRESS_AREA(FIRST_IT,NBJ,NEELEM,NORD,NPNE,NVJE,
     '          NVJP,NYNP,XAB,XP,YP,ERROR,*9999)
            ENDIF
            CALL CALC_RESIS_FLOW(KOUNT,NBJ,NEELEM,NENP,NORD,NPNE,NYNE,
     &        NVJE,
     '        XAB,XP,YP,ERROR,*9999)
C...  ARC 07-02-2011 ADDING COUPLED PERFUSION MODELS
            IF((COMPLIANCE.EQ.3.OR.COMPLIANCE.EQ.6).AND.
     &         COUPLE_VIA_LPM.EQ.'Y') THEN !Coupled perfusion - must calculate micro-vessel resistance
               DO noelem=1,NEELEM(0)
                 ne=NEELEM(noelem)
                 IF(XAB(nej_cap,ne).EQ.1.d0)THEN !check its a capillary
                   nh=NH_LOC(1,nx)
!!     NB/ Check na,nh values always=1 should set up more generally
!< !                     sum_capillaries=sum_capillaries+1
                   ne0=NXI(-1,1,ne)!upstream element number
                   P1=YP(NYNP(1,1,nh,NPNE(1,nb,ne),0,1),1) !pressure at start node of capillary element
                   P2=YP(NYNP(1,1,nh,NPNE(2,nb,ne),0,1),1)!pressure at end node of capillary element
                   Q01=YP(NYNE(1,nh,1,1,ne0),1) !flow in element upstream of capillary element !mm^3/s
                   IF(COMPLIANCE.EQ.3) THEN
                     Ppl=XP(1,1,7,NPNE(1,nb,ne))*98.06d0 !cm H2O -!Pa temporarily hard coded! 15.02.08
                     stretch=1.d0 !No mechanics therefore stretch factor=1.d0
                   ELSEIF(COMPLIANCE.EQ.6) THEN !Coupled to mechanics pleural & vessel stretch fields
                     Ppl=XP(1,1,nj_pleural,NPNE(1,nb,ne))*98.07d0 !cmH2O->Pa
                     IF(nj_stretch.EQ.0) THEN
                       stretch=1.d0 !Set to 1 i.e. no effect if undefined
                     ELSE
                       stretch=XP(1,1,nj_stretch,NPNE(1,nb,ne))
                     ENDIF
                   ENDIF
                   R_in=XP(1,1,nj_radius_R0,NPNE(1,nb,ne))*
     &                 sqrt(1.d0/stretch) !Including mechanics stretch factor (as is CALC_PRESS_AREA)
                   R_out=XP(1,1,nj_radius_R0,NPNE(2,nb,ne))*
     &                 sqrt(1.d0/stretch)
                   IF(nj_hypoxia.NE.0) THEN
                     k_factor=XP(1,1,nj_hypoxia,NPNE(1,nb,ne)) !Arterial constriction factor
                     IF(k_factor.GT.1.d0) THEN
                       !Use WRITE(OP_STRING...) and CALL WRITES for this
                       write(*,*) "CHANGE CODE k_factor",
     '                   NPNE(1,nb,ne),k_factor
                       !k_factor=1.d0 !If undefined (=2) then have equal to 1
                     ENDIF
                   ELSE
                     k_factor=1.d0 !If undefined Then have equal to 1 (no constriction)
                   ENDIF

                   !Length of up and downstream elts   
                   np1=NPNE(1,nb,ne0)
                   np2=NPNE(2,nb,ne0)
                   L_in=((XP(1,1,1,np2)-XP(1,1,1,np1))**2.d0+
     &                 (XP(1,1,2,np2)-XP(1,1,2,np1))**2.d0+
     &                 (XP(1,1,3,np2)-XP(1,1,3,np1))**2.d0)**0.5d0
                   L_out=L_in !Trees are currently always the same at terminals
                   x_cap=XP(1,1,1,NPNE(1,nb,ne))
                   y_cap=XP(1,1,2,NPNE(1,nb,ne))
                   z_cap=XP(1,1,3,NPNE(1,nb,ne))
                   DO nj=1,NJT
                      HEIGHT(nj)=XP(1,1,nj,NPNE(1,nb,ne))-
     &                    XP(1,1,nj,np_in) !gravitational head - same for art & vein now!
                   ENDDO
                   CALL CAP_FLOW_PARAM(ne,L_in,L_out,Ppl,
     &                R_in,R_out,stretch,ERROR,*9999)
                   IF(LADDER.EQ.1) THEN !Connections only at terminals
!                     CALL CAP_FLOW_SIMPLE(ne,HEIGHT,LPM_FUNC,P1,P2,
!     &                 Ppl,Q01,R_in,R_out,x_cap,y_cap,z_cap,ERROR,*9999)
                   ELSEIF(LADDER.EQ.2) THEN !Include the ladder model
                    CALL CAP_FLOW_LADDER(ne,HEIGHT,LPM_FUNC,P1,P2,
     &                   Ppl,Q01,R_in,R_out,x_cap,y_cap,z_cap,k_factor,
     &                   .FALSE.,ERROR,*9999)
                   ELSEIF(LADDER.EQ.3)THEN !multibranching acinus
                     CALL CAP_FLOW_MBA(ne,NORD,LPM_FUNC,P1,P2,Ppl,Q01,      
     &                   x_cap,y_cap,z_cap,.FALSE.,ERROR,*9999)
                   ENDIF 
                 XAB(nej_resis,ne)=LPM_FUNC
                 
                 ENDIF!Capillary element
              ENDDO !loop elements            
            ENDIF !Coupled perfusion
            ERR=ERR/NYNR(0,1,1) !sum of error divided by no of variables
            IF(ERR.LE.1.d-6.AND.(.NOT.FIRST_IT))THEN
c            IF(ERR.LE.CONVERG_TOL*1.d3.AND.(.NOT.FIRST_IT))THEN
              CONVERGED=.TRUE.
              WRITE(OP_STRING,
     &          '('' Convergence achieved after '',I5,'' iterations'')')
     &          KOUNT
              CALL WRITES(IOER,OP_STRING,ERROR,*9999)  
            ELSE !if error not converged
              IF(ERR.GE.MIN_ERR) THEN !.OR.ERR2.GE.MIN_ERR)
                N_MIN_ERR=N_MIN_ERR+1
              ELSE
                MIN_ERR=ERR
              ENDIF
              WRITE(OP_STRING,
     &          '('' Not converged, error = '',D14.6,
     '     D14.6)')
     '          ERR
              CALL WRITES(IOER,OP_STRING,ERROR,*9999)
            ENDIF !ERR
          ELSE !should not be non-linear
            CONVERGED=.TRUE.
          ENDIF !ITYP3(nr,nx)
          IF(ITYP3(nr,nx).LE.2) THEN
            IF(N_MIN_ERR.GT.10)THEN
              WRITE(OP_STRING,'('' Not converging, error is'',
     '          D14.6)') MAX(ERR,ERR_RAD)
              CALL WRITES(IOER,OP_STRING,ERROR,*9999)
              DO nonode=1,NPNODE(0) !put radial solution into CP
                np=NPNODE(nonode)
                DO i=2,7
                  CP(i,np)=CP(i+6,np)
                ENDDO !i
              ENDDO !nonode
              CONVERGED=.TRUE.
            ENDIF
          ELSE IF(ITYP3(nr,nx).GE.3) THEN
            IF(N_MIN_ERR.GT.40) THEN !takes longer to converge
              WRITE(OP_STRING,'('' Not converging, error is'',D14.6, 
     '     D14.6)')
     '          ERR
              CALL WRITES(IOER,OP_STRING,ERROR,*9999)
             IF(N_MIN_ERR.GT.60.OR.KOUNT.GT.70) THEN !not converged
                CONVERGED=.TRUE.
              ENDIF
            ENDIF
          ENDIF !ITYP3
        ELSE IF(LINEAR)THEN
          CONVERGED=.TRUE.
        ENDIF !LINEAR
        FIRST_IT=.FALSE.
      ENDDO !WHILE .NOT.CONVERGED
      UPDATE_VECTOR=.FALSE.
      UPDATE_MATRIX=.FALSE.
      IF(ITYP3(nr,nx).EQ.3.OR.
     &   (ITYP3(nr,nx).EQ.4.AND.COMPLIANCE.EQ.3)) THEN 

        WRITE(OP_STRING,'('' Error in solution: '' ,2(D14.6))')ERR,ERR2
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' # iterations '' ,I5)') KOUNT
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
      ENDIF
      CALL EXITS('SOLVE11')
      RETURN
 9999 CALL ERRORS('SOLVE11',ERROR)
      RET_ERROR=ERROR
      CALL EXITS('SOLVE11')
      RETURN 1
      END



