      SUBROUTINE INIT_PULM(nb,NBJ,nbreath,NEELEM,NHP,NHST,NKJE,
     '  NORD,NPNE,NPNODE,nr,NVJE,NVJP,NYNP,nx,NXI,N_TM,SPARSE_P,
     &  SPARSE_S,ACINUS,BBM,CE,MASS,SUM_N2,SUM_WATER,T,
     &  water_loss,XAB,XP,YP,DYNAM1,FIRST_A,FIX,GAS,LINEAR,
     &  SCALE,SMOOTHING,UPDATE_NY,UPDATE_VECTOR,ERROR,*)

C#### Subroutine: INIT_PULM
C###  Description:
C###    INIT_PULM sets up the initial parameters for gas mixing
C###    or temperature/water vapour distribution modelling in lung
C###    models.
C###    Acinus parameters stored in ACINUS array for gas exchange problems:
C###      ACINUS(nm_time,ne)=time elapsed for RBC in capillaries
C###      ACINUS(nm_avPAO2,ne)=breath-averaged alveolar air PO2
C###      ACINUS(nm_avPcO2,ne)=breath-averaged end-capillary PO2
C###      ACINUS(nm_PcO2,ne)=end-capillary PO2
C***  Created by Merryn Howatson Tawhai, April 1998

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grou00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'loc00.cmn'
      INCLUDE 'lungex00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'marc00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'moti00.cmn'
      INCLUDE 'pulm00.cmn'
      INCLUDE 'solv00.cmn'
      INCLUDE 'time02.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER nb,NBJ(NJM,NEM),nbreath,NEELEM(0:NE_R_M),NHP(NPM),NHST,
     &  NKJE(NKM,NNM,NJM,NEM),NORD(5,NE_R_M),
     &  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M),nr,NVJE(NNM,NBFM,NJM,NEM),
     &  NVJP(NJM,NPM),nx,NXI(-NIM:NIM,0:NEIM,0:NEM),N_TM,
     &  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM),SPARSE_P,SPARSE_S
      REAL*8 ACINUS(4,NEM),BBM(2,NEM),CE(NMM,NEM),
     &  MASS,SUM_N2,SUM_WATER,T,water_loss,XAB(NORM,NEM),
     &  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM)
      CHARACTER ERROR*(*)
      LOGICAL DYNAM1,FIRST_A,FIX(NYM,NIYFIXM),GAS,LINEAR,SCALE,
     &  SMOOTHING,UPDATE_NY,UPDATE_VECTOR
!     Local Variables
      INTEGER count_terms,IBEG,IEND,j,ne,nn,nonode,np,nv,noelem,
     &  VALUES(1,10)
      REAL*8 avtt,avtt2,sumQ,sumVc
      CHARACTER ERRMSG*50

      CALL ENTERS('INIT_PULM',*9999)
      SMOOTHING=.FALSE.
      UPDATE_NY=.TRUE.
      N_TM=1 !initialise # of time points for boundary condition
      DYNAM1=.TRUE. !Use GD
      N_SOLN=0
      UPDATE_VECTOR=.TRUE.
      scale=.false.
      NHST=NHP(NPNODE(1)) !# of dependent variables (1 or 2)
      SOLV_MAXIT(2)=SOLV_MAXIT(1)
      SOLV_TOL(2)=SOLV_TOL(1)
      BBM_MASS_CHANGE=0.0d0
      nm_time=1 !ACINUS(nm_time,ne)=time elapsed for RBC in capillaries
      nm_avPAO2=2 !ACINUS(nm_avPAO2,ne)=breath-averaged alveolar air PO2
      nm_avPcO2=3 !ACINUS(nm_avPcO2,ne)=breath-averaged end-capillary PO2
      nm_PcO2=4 !ACINUS(nm_PcO2,ne)=end-capillary PO2 at time t
      IF(ITYP3(nr,nx).EQ.1)THEN
        GAS=.TRUE.
        LINEAR=.TRUE.
        nbreath=1
        SUM_N2=0.d0
        CALL ASSERT(XP(1,1,nj_flow,NPNODE(1)).EQ.1.d0,
     &    '>>Flows must be fractional (nj_flow field). Use "> fem '//
     &    'update field # summation normalise" to define flow values'//
     &    ' as a fraction of the inlet flow.',ERROR,*9999)

C AJS 02/11 Add gas exchange initialisations
        IF(ITYP7(nr,nx).GT.1)THEN
!           BBM_MASS_CHANGE=0.0d0
!           InspExpIterations=0
          IF(nx.EQ.1)THEN !initialise at the start of an inspiration
            BreathIterations=0
            DO noelem=1,NEELEM(0)
              ne=NEELEM(noelem)
              IF(NXI(1,0,ne).EQ.0)THEN
                ACINUS(nm_avPAO2,ne)=CONC_INIT*(O2molVol*1.d3)*
     &            (P_ATM-p_water) !breath-averaged PAO2 for ne
                ACINUS(nm_avPcO2,ne)=MIXED_VENOUS_PO2!breath-averaged Pc'O2 for ne
                PAO2_BREATH=0.0d0 !spatial (over all acini) and temporal (over breath) average 
                PartO2_BREATH=0.0d0 !temporal (over breath) average PaO2
              ENDIF !terminal
            ENDDO!noelem
          ENDIF !nx
        ENDIF !ITYP7

      ELSE IF(ITYP3(nr,nx).EQ.2)THEN
        GAS=.FALSE.
        LINEAR=.FALSE.
        water_loss=0.d0
        SUM_WATER=0.d0
        TMIN=MIN(CORE_TEMP,MOUTH_TEMP,TEMP_IN)
        TMAX=MAX(CORE_TEMP,MOUTH_TEMP,TEMP_IN)
c        SMOOTHING=.TRUE.
          WHT(1)=0.5d0 !gamma inspiration (water)
          WHT(2)=2.d0 ! beta inspiration (heat)
          WHT(3)=0.5d0 !gamma expiration (water) !makes no difference!
          WHT(4)=2.d0 ! beta expiration (heat)
          WHT(5)=1.d0 !division inspiration
          WHT(6)=1.d0 !division expiration
      ELSE IF(ITYP3(nr,nx).EQ.4)THEN
        GAS=.FALSE.
        LINEAR=.FALSE.
      ELSE
        GAS=.FALSE.
        LINEAR=.TRUE.
      ENDIF

C SEN The Cmiss solvers all use Compressed-Row storage
      SPARSE_P=1
      SPARSE_S=1
      IF(.NOT.RESTART)THEN

C Blood flow and RBC transit time for each acinus (req'd in GAS_EXCH_ACINUS)
        IF(ITYP3(nr,nx).EQ.1.AND.ITYP7(nr,nx).GT.1)THEN
          PO2_arterial=MIXED_VENOUS_PO2
          avtt=0.0d0
          avtt2=0.0d0 !Q-weighted
          Q_ACINUS=0.0d0
          V_ACINUS=0.0d0
          sumVc=0.0d0
          count_terms=0
          DO noelem=1,NEELEM(0)
            ne=NEELEM(noelem)
            nb=NBJ(nj_radius,ne)
            DO nn=1,2
              np=NPNE(nn,nb,ne) 
              IF(NTB.GT.0.AND.NXI(1,0,ne).EQ.0)THEN ! terminal
                nv=NVJE(nn,nb,nj_flow,ne)
                IF(XP(1,nv,nj_flow,np).LT.0.0d0) WRITE(IOOP,
     &            '('' >>Negative ventilation for node '',I6)')np
                nv=NVJE(nn,nb,nj_Qdot,ne)
                IF(XP(1,nv,nj_Qdot,np).LT.0.0d0) WRITE(IOOP,
     &            '('' >>Negative perfusion for node '',I6)')np
              ELSE
                nv=NVJE(nn,nb,nj_flow,ne)
                IF(XP(1,nv,nj_flow,np).LT.0.0d0)THEN 
                  WRITE(ERRMSG,'('' >>Negative ventilation for node '',
     &              I6,'' (not a terminal node!!)'')')np
                  CALL ASSERT(.TRUE.,ERRMSG,ERROR,*9999)
                ENDIF
              ENDIF
            ENDDO !nn
          ENDDO!noelem
          DO noelem=1,NEELEM(0)
            ne=NEELEM(noelem)
            nb=NBJ(nj_flow,ne)
            np=NPNE(2,nb,ne) !terminal node
            IF(NXI(1,0,ne).EQ.0)THEN ! terminal
              nv=NVJE(2,nb,nj_Qdot,ne)
              XP(1,nv,nj_pob,np)=MIXED_VENOUS_PO2 !Initialise blood PO2
!               IF(SET_Q_DISTRIBUTION)THEN
!                 XP(1,nv,nj_tt,np)=ACINUS_VBLOOD/(XP(1,nv,nj_Qdot,np)*1.d-6) ! ttRBC = capillary blood volume/Q = litre/(litre/s)
!                 XP(1,1,nj_Vc,np)=ACINUS_VBLOOD !capillary blood volume
!               ENDIF
              avtt=avtt+XP(1,1,nj_tt,np)
              avtt2=avtt2+XP(1,1,nj_tt,np)*XP(1,nv,nj_Qdot,np)
              count_terms=count_terms+1
              Q_ACINUS=Q_ACINUS+XP(1,nv,nj_Qdot,np)
              V_ACINUS=V_ACINUS+XP(1,nv,nj_flow,np)
              sumVc=sumVc+XP(1,1,nj_Vc,np)
              ACINUS(nm_time,ne)=0.0d0!+DT
              ACINUS(nm_PcO2,ne)=MIXED_VENOUS_PO2
            ENDIF
          ENDDO!noelem
          avtt=avtt/count_terms
          avtt2=avtt2/Q_ACINUS
          Q_ACINUS=Q_ACINUS*CARDIAC_OUTPUT !fractional perfusion*cardiac output mm3/s
          V_ACINUS=V_ACINUS*INLET_FLOW(1) !fractional ventilation*inspiratory flow mm3/s
          write(*,*) 'INLET_FLOW(1)=',INLET_FLOW(1)

C Total acinar ventilation, perfusion and blood volume
C NB: Calculate V total assuming inspiration duration = expiration duration
C     therefore multiply by 30 seconds (half of one minute)
         WRITE(OP_STRING,'('' Total ventilation= '',E12.3,
     &      '' L/min,  Total perfusion= '',E12.3,'' L/min,  Total Vc= ''
     &      ,E12.3,''ml,  Mixed venous blood PO2='',F8.3,'' mmHg'')')
     &      V_ACINUS*1.d-6*30.d0,Q_ACINUS*1.d-6*60.d0,sumVc*1.d-3,
     &      MIXED_VENOUS_PO2
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C Average acinus values
          Q_ACINUS=Q_ACINUS/DBLE(count_terms)
          CALL SET_USER_DOUBLE('Q_ACINUS',Q_ACINUS,ERROR)
          V_ACINUS=V_ACINUS/DBLE(count_terms)
          CALL SET_USER_DOUBLE('V_ACINUS',V_ACINUS,ERROR)

          WRITE(OP_STRING,'('' Average acinar ventilation= '',E10.2,
     &      '' mm3/s  Average acinar perfusion= '',E10.2,'' mm3/s '//
     &      ' Average RBC transit time: (1) over '',I6,
     &      '' units= '',F6.3,''s  (2) Q-weighted= '',F6.3,'' s'')')
     &      V_ACINUS,Q_ACINUS,count_terms,avtt,avtt2
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

        ENDIF !ityp
      ENDIF ! not restart

      IF(RESTART)THEN
C        T=T_LAST
C        FIRST_A=.FALSE.
         T=TSTART
         FIRST_A=.TRUE.
        IF(IWRIT1(nr,nx).NE.0)THEN !have defined frequency of output
C     Open output files
          CALL STRING_TRIM(FILE02,IBEG,IEND)
          IF(GAS.AND..NOT.OPENF_HISTORY)THEN
            OPENF_HISTORY=.TRUE.
            IF(LUNG_OP.GE.1)THEN
              CALL OPENF(IOFILE2,'DISK',FILE02(IBEG:IEND)
     &          //'.history','NEW','SEQUEN','FORMATTED',160,ERROR,*9999)
            ENDIF
c            IF(GAS)THEN
c              IF(LUNG_OP.GE.2)THEN
c                CALL OPENF(IOFILE3,'DISK',FILE02(IBEG:IEND)//'.snhist',
c     '            'NEW','SEQUEN','FORMATTED',160,ERROR,*9999)
c              ENDIF
c              IF(LUNG_OP.GE.3)THEN
c                CALL OPENF(IOFILE4,'DISK',FILE02(IBEG:IEND)//'.acinus',
c     '            'NEW','SEQUEN','FORMATTED',160,ERROR,*9999)
c              ENDIF !GAS
c            ENDIF
          ENDIF

           IF(GAS.AND.ITYP7(nr,nx).GT.1)THEN !gas exchange
              IF(LUNG_OP.GE.1)THEN
                CALL OPENF(IOFILE2,'DISK',
     '            FILE02(IBEG:IEND)//'.history','OLD',
     '            'APPEND','FORMATTED',160,ERROR,*9999)
                OPENF_HISTORY=.TRUE. !lung history file has been opened
              ENDIF
              IF(LUNG_OP.GE.2)THEN
                CALL OPENF(IOFILE3,'DISK',FILE02(IBEG:IEND)//'.error',
     '            'OLD','APPEND','FORMATTED',160,ERROR,*9999)
              ENDIF
              IF(LUNG_OP.GE.3)THEN
                CALL OPENF(IOFILE4,'DISK',FILE02(IBEG:IEND)//'.nodes',
     '            'OLD','APPEND','FORMATTED',160,ERROR,*9999)
              ENDIF
            ENDIF

        ENDIF !IWRIT1
      ELSE !not restart

        FIRST_A=.TRUE.
        T=TSTART
        MASS=0.d0
        CUMULATIVE_VOLUME=0.d0
        SUM_EXP_VOL=0.d0 !initialise the cumulative expired volume

        IF(ITYP3(nr,nx).LE.2)THEN
          IF(INITIAL_VOLUME.LT.1.d6)THEN !small problems
            WRITE(OP_STRING,'('' Initial volume is'',E12.3,'' ml'')')
     &        INITIAL_VOLUME/1.d3
          ELSE
            WRITE(OP_STRING,'('' Initial volume is'',E12.3,'' L'')')
     &        INITIAL_VOLUME/1.d6
          ENDIF
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)


C...Format for screen output (gas mixing + gas exchange)
          IF(GAS.AND.ITYP7(nr,nx).GT.1)THEN
            WRITE(OP_STRING,'(''  Time   VolErr%   MassErr%   '//
     &        'inletPO2     avePAO2       PaO2      '//
     &        'DLO2(ml/s/mmHg)   O2uptake(ml/min)'')') 
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF

          DV_TOTAL(nx)=0.d0 !AJS 11/2010 Add solve class for DV_TOTAL
          BSUM=0.d0
          CUM=0.d0 !cumulative expired volume
          OPENF_HISTORY=.FALSE.
          IF(IWRIT1(nr,nx).NE.0)THEN !have defined frequency of output
C     Open output files
            CALL STRING_TRIM(FILE02,IBEG,IEND)

            IF(GAS.AND.ITYP7(nr,nx).GT.1)THEN !gas exchange
              IF(LUNG_OP.GE.1)THEN
                CALL OPENF(IOFILE2,'DISK',
     '            FILE02(IBEG:IEND)//'.history','NEW',
     '            'SEQUEN','FORMATTED',160,ERROR,*9999)
                WRITE(IOFILE2,
     &            '(''     Time    PAO2      PaO2     Inlet_PO2 '')') 
                OPENF_HISTORY=.TRUE. !lung history file has been opened
              ENDIF
              IF(LUNG_OP.GE.2)THEN
                CALL OPENF(IOFILE3,'DISK',FILE02(IBEG:IEND)//'.error',
     '            'NEW','SEQUEN','FORMATTED',160,ERROR,*9999)
              ENDIF
              IF(LUNG_OP.GE.3)THEN
                CALL OPENF(IOFILE4,'DISK',FILE02(IBEG:IEND)//'.nodes',
     '            'NEW','SEQUEN','FORMATTED',160,ERROR,*9999)
                DO nonode=1,NP_OP_LIST(0)
                  VALUES(1,nonode)=NP_OP_LIST(nonode)
                ENDDO !nonode
                WRITE(IOFILE4,'(''   Time'',30(I8))') (VALUES(1,j),j=1,
     &          NP_OP_LIST(0))
              ENDIF

            ELSE 

              CALL STRING_TRIM(FILE02,IBEG,IEND)
              CALL OPENF(IOFILE3,'DISK',FILE02(IBEG:IEND)//'.snhist',
     '          'NEW','SEQUEN','FORMATTED',160,ERROR,*9999)
              WRITE(IOFILE3,'('' Mass error '')') 

              IF(LUNG_OP.GE.1)THEN
                CALL OPENF(IOFILE2,'DISK',
     '          FILE02(IBEG:IEND)//'.history','NEW',
     '          'SEQUEN','FORMATTED',160,ERROR,*9999)
                WRITE(IOFILE2,'(''   Cumul_exp_vol  Inlet_conc'')') 
                OPENF_HISTORY=.TRUE. !lung history file has been opened
              ENDIF
              IF(GAS)THEN
                IF(LUNG_OP.GE.2)THEN
                  CALL OPENF(IOFILE3,'DISK',
     '              FILE02(IBEG:IEND)//'.snhist','NEW','SEQUEN',
     '              'FORMATTED',160,ERROR,*9999)
                ENDIF
                IF(LUNG_OP.GE.3)THEN
                  CALL OPENF(IOFILE4,'DISK',
     '              FILE02(IBEG:IEND)//'.acinus','NEW','SEQUEN',
     '              'FORMATTED',160,ERROR,*9999)
                ENDIF
              ENDIF !GAS
            ENDIF
          ENDIF !IWRIT1
c          IF(.NOT.GAS)THEN !water and heat transfer
c          IF(ne_trachea.EQ.1) SMOOTHING=.TRUE.
c            DO nlpm=1,NTB !scale BBM volumes and calculate ratio to total
c              XAB(3,nlpm)=XAB(2,nlpm) !reinitialise the volume
c              XAB(4,nlpm)=XAB(2,nlpm) !volume at previous time-step
c              XAB(5,nlpm)=CONC_INIT !concentration of washin gas
c              XAB(6,nlpm)=CONC_INIT*XAB(3,nlpm) !mass in XAB at st insp
c            ENDDO !nlpm
c          ENDIF
          MASS=0.d0
          nb=NBJ(1,NEELEM(1))

        ENDIF
      ENDIF !RESTART


      IF(ITYP3(nr,nx).LE.2)THEN
C     Calculate the inital mass of washin gas
        CALL CALCMASS(NBJ,NEELEM,NH_LOC(1,nx),NORD,NPNE,NVJE,NXI,
     '    NYNP(1,1,1,1,0,1),BBM,MASS,XAB,XP,YP(1,1),ERROR,*9999)
      ENDIF
      IF(.NOT.RESTART) MASS_0=MASS

      CALL EXITS('INIT_PULM')
      RETURN
 9999 CALL ERRORS('INIT_PULM',ERROR)
      CALL EXITS('INIT_PULM')
      RETURN 1
      END


