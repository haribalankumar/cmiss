      SUBROUTINE POST_SOL(nb,NBJ,NEELEM,NELIST,NELIST2,NENP,NORD,NPNE,
     &  NPNODE,nr,NVJE,nx,NXI,NYNE,NYNP,ACINUS,BBM,CE,CP,T,XAB,XP,
     &  YP,GAS,UPDATE_MATRIX,ERROR,*)

C#### Subroutine: POST_SOL
C###  Description:
C###    POST_SOL does post-solution updates for pulmonary transport
C###    problems.

      IMPLICIT NONE
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'loc00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'lungex00.cmn'
      INCLUDE 'pulm00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M),NELIST(0:NEM),
     '  NENP(NPM,0:NEPM),NORD(5,NE_R_M),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M),nr,NVJE(NNM,NBFM,NJM,NEM),nx,
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM)
      REAL*8 ACINUS(4,NEM),BBM(2,NEM),CE(NMM,NEM),CP(NMM,NPM),MASS,T,
     '  XAB(NORM,NEM),XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM)
      LOGICAL GAS,UPDATE_MATRIX
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,ne,nh,ny,gen,
     '  np2,i,j,nextgen,NELIST2(0:NEM),enum
      INTEGER*4 PATHWAY_PTR,PLOT_DATA_PTR,TIME_PTR
      REAL*8 genresis(24),bresis
      LOGICAL BRANCH

      CALL ENTERS('POST_SOL',*9999)

      IF(ITYP3(nr,nx).LE.2) THEN !airways
        CALL CALCMASS(NBJ,NEELEM,NH_LOC(1,nx),NORD,NPNE,NVJE,NXI,
     &    NYNP,BBM,MASS,XAB,XP,YP(1,1),ERROR,*9999)
        MASS_CURRENT(nx)=MASS

        IF(ITYP7(nr,nx).EQ.1)THEN !not gas exchange
          IDEAL_MASS(nx)=IDEAL_MASS(nx)+INLET_FLOW(nx)*DT*YP(1,1) !change better Y
        ELSE !gas exchange
!         write(*,*) 'IDEAL_MASS=',IDEAL_MASS(nx),' O2_EXCHANGED=',
!      &    O2_EXCHANGED
          IDEAL_MASS(nx)=IDEAL_MASS(nx)+INLET_FLOW(nx)*DT*YP(1,1)+
     &      O2_EXCHANGED !mass at start of time step + inspired mass + exchanged mass
        ENDIF

        IF(ITYP7(nr,nx).GT.1) CALL GAS_EXCHANGE_METRICS(nb,NEELEM,NPNE,
     &    NPNODE,nx,NXI,NVJE,CE,ACINUS,BBM,T,XP,ERROR,*9999)

        IF(INLET_FLOW(nx).LT.0.d0) SUM_EXP_VOL=SUM_EXP_VOL
     &    +DABS(INLET_FLOW(nx)*DT)
        CALL OPLUNG(nb,NEELEM,NPNE,NPNODE,nr,nx,NVJE,NYNP,BBM,CE,CP,
     &    T,XAB,XP,YP,GAS,ERROR,*9999)
        IF(N_SOLN.EQ.IWRIT6)THEN
          N_SOLN=0
        ENDIF
      ELSE IF(ITYP3(nr,nx).EQ.3) THEN !pulmonary capillary blood flow
        PATHWAY_PTR=0
        PLOT_DATA_PTR=0
        TIME_PTR=0
        CALL ALLOCATE_MEMORY((MAX_PATH+1)*FACTORS,1,DPTYPE,PATHWAY_PTR,
     '    MEM_INIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(4500*FACTORS,1,DPTYPE,PLOT_DATA_PTR,
     '    MEM_INIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(NEM,1,DPTYPE,TIME_PTR,MEM_INIT,ERROR,
     &    *9999)        
        CALL CALC_TRANSIT_TIME(nb,NEELEM,NENP,NPNE,
     '    NYNE,NYNP,CE,%VAL(PATHWAY_PTR),%VAL(PLOT_DATA_PTR),
     &    %VAL(TIME_PTR),YP(1,1),ERROR,*9999)
         !calculates RBC and neutrophil transit times
        CALL FREE_MEMORY(PATHWAY_PTR,ERROR,*9999)
        CALL FREE_MEMORY(PLOT_DATA_PTR,ERROR,*9999)
        CALL FREE_MEMORY(TIME_PTR,ERROR,*9999)
C... outputs model solution data, geometric, pressures and flow
C... results for model analysis.        
        CALL OPPCAP(nb,NEELEM,NENP,NPNE,NYNE,NYNP,CE,YP,ERROR,
     '    *9999)
        UPDATE_MATRIX=.TRUE.
      ELSE IF(ITYP3(nr,nx).EQ.5)THEN !MCT
        WRITE(OP_STRING,'('' Time='',F6.3)') T
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
      ELSE
        UPDATE_MATRIX=.TRUE.
      ENDIF !ITYP3
     
      CALL EXITS('POST_SOL')
      RETURN

 9999 CALL ERRORS('POST_SOL',ERROR)
      CALL EXITS('POST_SOL')
      RETURN 1
      END



