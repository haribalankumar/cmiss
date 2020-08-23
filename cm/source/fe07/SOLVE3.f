      SUBROUTINE SOLVE3(ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,ISR_GKK,ISR_GQ,
     '  LGE,NBH,NENP,NHE,NONY,NP_INTERFACE,NPNE,NPNY,nr,NRE,NVHE,nx,
     '  NYNE,NYNO,NYNP,NYNR,TIME_STEP,CONY,CYNO,GD,GK,GKK,GM,GQ,GR,
     '  GRR,XO,YP,DYNAM1,DYNAM2,FIRST_A,FIX,UPDATE_MATRIX,ERROR,*)

C#### Subroutine: SOLVE3
C###  Description:
C###    SOLVE3 solves time dependent advection diffusion equations.
C###    It returns the increment to add to the solution vector in
C###    YP(ny,1). Hence the solution at time T+dT is the current
C###    solution + the increment.

C**** Matrices are assembled in ASSEMBLE3 and just reorganised here.
C**** TIME_STEP controls the type of time integration i.e.
C****   TIME_STEP=0 : Static solution of initial accelerations.
C****       "     1 : first  order (Crank-Nicholson) time integration.
C****       "     2 : second order (Newmark; Wood-Zienkiewicz etc)
C****       "     3 : third  order (Houbolt; Hilbert-Hughes etc)
C**** On entry YP(ny,1) contains b.c.s, defined by FIX(ny,1) or the
C**** current estimate of soln for nonlinear case. On exit YP(ny,1)
C**** contains the incremental solution.

      IMPLICIT NONE
      INCLUDE 'b08.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'solv00.cmn'
      INCLUDE 'time02.cmn'
!     Parameter List
      INTEGER ISC_GK(NISC_GKM),ISC_GKK(NISC_GKKM),ISC_GQ(NISC_GQM),
     '  ISR_GK(NISR_GKM),ISR_GKK(NISR_GKKM),ISR_GQ(NISR_GQM),
     '  LGE(NHM*NSM,NRCM),NBH(NHM,NCM,NEM),NENP(NPM,0:NEPM,0:NRM),
     '  NHE(NEM),NONY(0:NOYM,NYM,NRCM),NP_INTERFACE(0:NPM,0:3),
     '  NPNE(NNM,NBFM,NEM),
     '  NPNY(0:6,NYM,0:NRCM),nr,NRE(NEM),NVHE(NNM,NBFM,NHM,NEM),nx,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),NYNO(0:NYOM,NOOPM,NRCM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),NYNR(0:NY_R_M,0:NRCM,NCM),
     '  TIME_STEP
      REAL*8 CONY(0:NOYM,NYM,NRCM),CYNO(0:NYOM,NOOPM,NRCM),GD(NZ_GD_M),
     '  GK(NZ_GK_M),GKK(NZ_GKK_M),GM(NZ_GM_M),GQ(NZ_GQ_M),GR(NYROWM),
     '  GRR(NOM),XO(NOM),YP(NYM,NIYM)
      LOGICAL DYNAM1,DYNAM2,FIRST_A,FIX(NYM,NIYFIXM),UPDATE_MATRIX
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER GETNYR,no1,no2,no_nynr1,no_nynr2,
     '  noy1,noy2,np,ny,ny1,ny2,ny3,nyo2,nyr,nz,nzz
      INTEGER*4 WORK_PTR
      REAL ELAPSED_TIME,TIME_START1(1),TIME_START2(1),TIME_STOP(1)
      REAL*8 AA,BB,co1,co2,SUM
      LOGICAL MEM_INIT,X_INIT

      CALL ENTERS('SOLVE3',*9999)

C cpb 28/3/96 Rewritting SOLVE3. The old routine is in ARCHIVE.

      IF(NOT(2,1,nr,nx).EQ.0) THEN
        ERROR=' >>The number of unknowns is zero'
        GOTO 9999
      ENDIF

      IF(TIME_STEP.EQ.0) THEN
        CALL ASSERT(DYNAM2,'>>DYNAM2 is not set',ERROR,*9999)
      ELSE IF(TIME_STEP.EQ.1) THEN
        CALL ASSERT(DYNAM1,'>>DYNAM1 is not set',ERROR,*9999)
      ELSE IF(TIME_STEP.EQ.2) THEN
        CALL ASSERT(DYNAM1.AND.DYNAM2,'>>DYNAM1 and DYNAM2 are not '
     '    //'both set',ERROR,*9999)
      ELSE IF(TIME_STEP.EQ.3) THEN
        CALL ASSERT(DYNAM1.AND.DYNAM2,'>>DYNAM1 and DYNAM2 are not '
     '    //'both set',ERROR,*9999)
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_START1)
      CALL CPU_TIMER(CPU_USER,TIME_START2)

C*** Setup and initialise arrays

      IF(UPDATE_MATRIX) THEN
        IF(FIRST_A) THEN
          IF(SPARSEGKK(nx).NE.0) THEN
            WORK_PTR=0
            MEM_INIT=.FALSE.
            CALL ALLOCATE_MEMORY(NOT(1,1,nr,nx)*NOT(2,1,nr,nx),1,
     '        CHARTYPE,WORK_PTR,MEM_INIT,ERROR,*9999)
          ENDIF
          CALL CALC_SPARSE_GKK(ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,ISR_GKK,
     '      ISR_GQ,LGE,NBH,NENP,NHE,NOT(1,1,nr,nx),NOT(2,1,nr,nx),
     '      NONY,NP_INTERFACE,NPNE,NPNY,nr,NRE,NVHE,nx,NYNE,NYNP,
     '      NYNR,GK,GQ,%VAL(WORK_PTR),.FALSE.,.TRUE.,ERROR,*9999)
          IF(SPARSEGKK(nx).NE.0) CALL FREE_MEMORY(WORK_PTR,ERROR,*9999)
        ENDIF
        DO nzz=1,NZZT(1,nr,nx)
          GKK(nzz)=0.0d0
        ENDDO !nzz
      ENDIF
      DO no1=1,NOT(1,1,nr,nx)
        GRR(no1)=0.0d0
      ENDDO !no1

C*** Generate reduced system

C*** Update RHS vector GRR from flux b.c.s

      DO no_nynr1=1,NYNR(0,0,2) !loop over global variables for nc=2
        ny1=NYNR(no_nynr1,0,2) !is global (flux) variable number
        IF(FIX(ny1,1)) THEN !flux is set as a b.c.
          ny2=GETNYR(2,NPNY,nr,1,0,ny1,NYNE,NYNP) !glob (flux) row #
          DO noy1=1,NONY(0,ny2,1)
            no1=NONY(noy1,ny2,1)
            co1=CONY(noy1,ny2,1)
            GRR(no1)=GRR(no1)+YP(ny1,1)*co1 !GRR val to appl flux bc
          ENDDO !noy1
        ENDIF
      ENDDO !no_nynr1 (ny1)

      DO no_nynr1=1,NYNR(0,1,1) !Loop over global rows of GK
        ny1=NYNR(no_nynr1,1,1) !is row #
        DO noy1=1,NONY(0,ny1,1) !loop over #no's attached to row ny1
          no1=NONY(noy1,ny1,1) !is no# attached to row ny1
          co1=CONY(noy1,ny1,1) !is coupling coeff for row mapping
C                               ie row_no1=a*row_ny1+b*row_ny2
          BB=GR(ny1) !get reduced R.H.S.vector
          DO no_nynr2=1,NYNR(0,0,1) !loop over the #cols of GK
            ny2=NYNR(no_nynr2,0,1) !is global variable #
            ny3=GETNYR(1,NPNY,nr,2,0,ny2,NYNE,NYNP) !local GK var #
            CALL SPARSE(ny1,ny3,NYT(1,1,nx),nz,NZ_GK_M,
     '        NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
            IF(nz.NE.0) THEN
              IF(TIME_STEP.GE.1) BB=BB-GK(nz)*YP(ny2,8)
              IF(TIME_STEP.GE.2) BB=BB-GD(nz)*YP(ny2,12)
              IF(TIME_STEP.EQ.3) BB=BB-GM(nz)*YP(ny2,15)
              IF(.NOT.DYNAM1.AND..NOT.DYNAM2) THEN
                AA=GK(nz)
              ELSE
                IF(TIME_STEP.EQ.0) THEN
                  AA=GM(nz)
                ELSE IF(TIME_STEP.EQ.1) THEN
                  AA=GD(nz)+A1*GK(nz)
                ELSE IF(TIME_STEP.EQ.2) THEN
                  AA=GM(nz)+A1*GD(nz)+A2*GK(nz)
                ELSE IF(TIME_STEP.EQ.3) THEN
                  AA=A1*GM(nz)+A2*GD(nz)+A3*GK(nz)
                ENDIF
              ENDIF
              IF(FIX(ny2,1).AND..NOT.FIX(ny2,3)) BB=BB-AA*YP(ny2,9)
              DO noy2=1,NONY(0,ny2,2) !loop over #no's for var ny2
                no2=NONY(noy2,ny2,2) !no# attached to ny2
                co2=CONY(noy2,ny2,2) !coup coeff for the col mappng
C                                     i.e. var_no1=a*var_ny1+b*var_ny2
                CALL SPARSE(no1,no2,NOT(1,1,nr,nx),nzz,NZ_GKK_M,
     '            NZZT(1,nr,nx),ISC_GKK,ISR_GKK,SPARSEGKK(nx),
     '            ERROR,*9999)
                IF(nzz.NE.0.AND.UPDATE_MATRIX) GKK(nzz)=GKK(nzz)+
     '            AA*co1*co2
              ENDDO !noy2
            ENDIF
C cpb 1/4/96 No additive constant implemented yet.
          ENDDO !no_nynr2
          GRR(no1)=GRR(no1)+BB*co1
        ENDDO !noy1
      ENDDO !no_nynr1

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
      IF(IWRIT4(nr,nx).GE.1) THEN
        WRITE(OP_STRING,'(/'' CPU time for solution matrix '
     '    //'initialisation and assembly: '',D11.4,'' s'')')
     '    ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_START2)
      IF(KTYP4.NE.0.AND.FIRST_A) THEN !Output global matrices
        CALL WRITE_SOL_MATRIX(ISC_GKK,ISR_GKK,nr,nx,GKK,GRR,
     '    ERROR,*9999)
      ENDIF

      IF(IWRIT4(nr,nx).GE.3) THEN
        IF(UPDATE_MATRIX) THEN
          WRITE(OP_STRING,'(/'' Global load vector GRR & stiffness '
     '      //'matrix GKK:'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NOT(1,1,nr,nx)='',I5,'
     '      //''', NOT(2,1,nr,nx)='',I5)') NOT(1,1,nr,nx),
     '      NOT(2,1,nr,nx)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ELSE
          WRITE(OP_STRING,'(/'' Global load vector GRR:'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NOT(1,1,nr,nx)='',I5)') NOT(1,1,nr,nx)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        CALL OPSTFMAT(NYNR(0,1,1),ISC_GKK,ISR_GKK,IOOP,NOT(1,1,nr,nx),
     '    NOT(2,1,nr,nx),NZZT(1,nr,nx),NYNR(0,2,1),SPARSEGKK(nx),GKK,
     '    GRR,'GKK','GRR',.TRUE.,UPDATE_MATRIX,.TRUE.,ERROR,
     '    *9999)
      ENDIF

       CALL CPU_TIMER(CPU_USER,TIME_STOP)
       ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
       IF(IWRIT4(nr,nx).GE.1) THEN
        IF(IWRIT4(nr,nx).GE.4.OR.KTYP4.NE.0) THEN
          WRITE(OP_STRING,'(/'' CPU time for solution matrix '
     '      //'output: '',D11.4,'' s'')') ELAPSED_TIME
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF

C*** Solve reduced system

      CALL CPU_TIMER(CPU_USER,TIME_START2)

      ! Use zero as an initial solution guess for the first iteration.
      ! Otherwise use the solution from the previous time.
      ! A separate variable to FIRST is required because FIRST is
      ! modified in SOLVE_SYSTEM.
      X_INIT=FIRST_A

      CALL SOLVE_SYSTEM(ISC_GKK,ISR_GKK,NOT(1,1,nr,nx),NOT(1,1,nr,nx),
     '  NOT(2,1,nr,nx),NZZT(1,nr,nx),IWRIT4(nr,nx),PRECON_CODE(nx),
     '  SOLVEROPTION(nx),SPARSEGKK(nx),GKK,GRR,XO,FIRST_A,UPDATE_MATRIX,
     '  X_INIT,nx,ERROR,*9999)

C*** Put solution values into YP(ny,1)

      DO no2=1,NOT(2,1,nr,nx)
        DO nyo2=1,NYNO(0,no2,2)
          ny2=NYNO(nyo2,no2,2)
          co2=CYNO(nyo2,no2,2)
          IF(NPNY(0,ny2,0).EQ.1) THEN
            np=NPNY(4,ny2,0)
C GMH 8/1/97 Update cmgui link
            CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
          ENDIF
          YP(ny2,1)=XO(no2)*co2
        ENDDO !nyo2
      ENDDO !no2

C*** Backsubstitute to find unknown fluxes

      DO no_nynr1=1,NYNR(0,0,1) !Loop over the global variables of nc=1
        ny=NYNR(no_nynr1,0,1) !is global variable number
        IF(FIX(ny,1)) THEN !global variable was set as a bc
          ny1=NYNR(no_nynr1,1,1) !row corresponding the the variable
          SUM=0.0d0
          DO no_nynr2=1,NYNR(0,2,1) !loop over local columns
            ny2=NYNR(no_nynr2,2,1) !is local column number
            ny3=GETNYR(1,NPNY,nr,0,2,ny2,NYNE,NYNP)
            !is global variable number for local column ny2
            CALL SPARSE(ny1,ny3,NYT(1,1,nx),nz,NZ_GK_M,
     '        NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
            IF(nz.NE.0) SUM=SUM+GK(nz)*YP(ny3,1)
          ENDDO
          nyr=GETNYR(2,NPNY,nr,0,0,ny,NYNE,NYNP) !correspond. rhs var. #
          IF(NPNY(0,nyr,0).EQ.1) THEN
            np=NPNY(4,nyr,0)
C GMH 8/1/97 Update cmgui link
            CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
          ENDIF
          YP(nyr,1)=SUM
        ENDIF
      ENDDO !no_nynr

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
      IF(IWRIT4(nr,nx).GE.1) THEN
        WRITE(OP_STRING,'(/'' CPU time for storage and back '
     '    //'substitution: '',D11.4,'' s'')') ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START1(1)
      IF(IWRIT4(nr,nx).GE.1) THEN
        WRITE(OP_STRING,'(/'' Total CPU time for solution: '','
     '    //'D11.4,'' s'')') ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('SOLVE3')
      RETURN
 9999 CALL ERRORS('SOLVE3',ERROR)
      CALL EXITS('SOLVE3')
      RETURN 1
      END


