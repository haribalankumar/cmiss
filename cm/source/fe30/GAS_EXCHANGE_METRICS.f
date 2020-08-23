      SUBROUTINE GAS_EXCHANGE_METRICS(nb,NEELEM,NPNE,NPNODE,nx,NXI,NVJE,
     &  CE,ACINUS,BBM,T,XP,ERROR,*)

C#### Subroutine: GAS_EXCHANGE_METRICS
C###  Description:
C###    GAS_EXCHANGE_METRICS 
C###    Calculates gas exchange outputs
C***  Created by AJS, Feb 2011

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'loc00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lungex00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'time02.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER nb,NEELEM(0:NE_R_M),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M),
     &  nx,NXI(-NIM:NIM,0:NEIM,0:NEM),NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 ACINUS(4,NEM),CE(NMM,NEM),BBM(2,NEM),T,XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ne,noelem,nonode,np,nv
      REAL*8 BBMtotal,BBMvol(10),Vtotal

      CALL ENTERS('GAS_EXCHANGE_METRICS',*9999)

C Initialise
      BBMtotal=0.0d0!initialise
      PAO2_AVERAGE=0.0d0
!       PAO2_EXPIRED=0.0d0
      Vtotal=0.0d0

C Find minimum, maximum, and mean PAO2 over inspiration/expiration
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        nv=NVJE(2,nb,nj_poa,ne)
        np=NPNE(2,nb,ne) !terminal node
        IF(NXI(1,0,ne).EQ.0.AND.XP(1,1,nj_alveoli,np).EQ.1.d0)THEN ! terminal
          ACINUS(nm_avPAO2,ne)=ACINUS(nm_avPAO2,ne)+XP(1,nv,nj_poa,np)
          ACINUS(nm_avPcO2,ne)=ACINUS(nm_avPcO2,ne)+ACINUS(nm_PcO2,ne)
          BBMtotal=BBMtotal+BBM(1,ne)
          Vtotal=Vtotal+XP(1,nv,nj_flow,np)
          PAO2_AVERAGE=PAO2_AVERAGE+BBM(1,ne)*XP(1,nv,nj_poa,np) !overall average
!           PAO2_EXPIRED=PAO2_EXPIRED+XP(1,nv,nj_flow,np)*
!      &      XP(1,nv,nj_poa,np) !expired average
        ENDIF
      ENDDO !noelem
      PAO2_AVERAGE=PAO2_AVERAGE/BBMtotal !volume-weighted average at time T
!       PAO2_EXPIRED=PAO2_EXPIRED/Vtotal !ventilation-weighted average at time T
      PAO2_BREATH=PAO2_BREATH+PAO2_AVERAGE !spatial & temporal average
      PartO2_BREATH=PartO2_BREATH+PO2_arterial !temporally-averaged PaO2

      IF(nx.EQ.2.AND.(TFINISH-T).LT.DT)THEN !end of expiration
        PAO2_BREATH=PAO2_BREATH/BreathIterations
        PartO2_BREATH=PartO2_BREATH/BreathIterations
        DO noelem=1,NEELEM(0)
          ne=NEELEM(noelem)
          nv=NVJE(2,nb,nj_poa,ne)
          np=NPNE(2,nb,ne) !terminal node
          IF(NXI(1,0,ne).EQ.0.AND.XP(1,1,nj_alveoli,np).EQ.1.d0)THEN ! terminal
            ACINUS(nm_avPAO2,ne)=ACINUS(nm_avPAO2,ne)/BreathIterations !breath-averaged PAO2 for ne
            ACINUS(nm_avPcO2,ne)=ACINUS(nm_avPcO2,ne)/BreathIterations !breath-averaged Pc'O2 for ne
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('GAS_EXCHANGE_METRICS')
      RETURN

 9999 CALL ERRORS('GAS_EXCHANGE_METRICS',ERROR)
      CALL EXITS('GAS_EXCHANGE_METRICS')
      RETURN 1
      END



