      SUBROUTINE PRE_SOL(nb,NBJ,NEELEM,NENP,NORD,NPNE,nr,NVJE,nx,
     &  NXI,NYNE,NYNP,BBM,CE,XAB,XP,YP,T,FIRST_A,GAS,
     &  SCALE,UPDATE_MATRIX,UPDATE_VECTOR,update_ny,ERROR,*)

C#### Subroutine: PRE_SOL
C###  Description:
C###    PRE_SOL does pre-solution updates for pulmonary transport
C###    problems.
      
      IMPLICIT NONE
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'loc00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'pulm00.cmn'
!     Parameter List
      INTEGER nb,NBJ(NJM,NEM),NEELEM(0:NE_R_M),NENP(NPM,0:NEPM),
     &  NORD(5,NE_R_M),NPNE(NNM,NBFM,NEM),nr,
     '  NVJE(NNM,NBFM,NJM,NEM),nx,NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM)
      REAL*8 BBM(2,NEM),CE(NMM,NEM),XAB(NORM,NEM),
     '  XP(NKM,NVM,NJM,NPM),YP(NYM)
      LOGICAL FIRST_A,GAS,SCALE,UPDATE_MATRIX,
     &  UPDATE_VECTOR,update_ny
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 NUM_BLOCKED,T,TOTAL_NE

      CALL ENTERS('PRE_SOL',*9999)

      IF(T.GT.0.05d0) SCALE=.TRUE.

      T_cycle=T-0.5d0*DT-BSUM !do calculations at mid time-step

      IF(ITYP3(nr,nx).LE.2) THEN !AdvDiff
        UPDATE_VECTOR=.TRUE.
        UPDATE_NY=.TRUE.

      ELSE IF(ITYP3(nr,nx).EQ.2) THEN ! coupled
        UPDATE_VECTOR=.TRUE.
        UPDATE_NY=.TRUE.
C...AJS Particle stuff not committed to global cmiss (NJW's code...)
!         IF(ITYP19(nr,nx).EQ.2)THEN !particles
! c        CALL MESH_FLOW(NBJ,NEELEM,NPNE,NVJE,nx,NXI,XP,ERROR,*9999)
! c call subroutine to update L1 and apparent diffusion coeff
! C store ADC in CE(1,ne). store L1 in CE(2,ne).
!           CALL UPDATE_PARTICLES(NBJ,NEELEM,NORD,NPNE,NVJE,nx,NXI,CE,XP,
!      &      ERROR,*9999)
!         ENDIF

      ELSE IF(ITYP3(nr,nx).GE.3) THEN !capillaries
        UPDATE_MATRIX=.TRUE.  !iterating with hematocrit
        FIRST_A=.TRUE.
C... Allocates elements blocked with white blood cells (WBCs)
C... Boundary conditions already put into YP in IPINI3
C... Initialises hemodynamic properties
        IF(ITYP3(nr,nx).EQ.3)THEN
          CALL WBC_BLOCK(nb,NEELEM,NPNE,NYNP,CE,YP,ERROR,*9999)
          CALL CALC_HEMODYNAMICS(ITYP3(nr,nx),nb,NEELEM,NENP,NPNE,NVJE,
     &      CE,XP,ERROR,*9999)
          NUM_BLOCKED=NE_BLOCK(0)
          TOTAL_NE=NEELEM(0)
          NUM_BLOCKED=(NUM_BLOCKED/TOTAL_NE)*100.d0
          WRITE(OP_STRING,'('' Percentage of blocked segments '//
     '    'in whole mesh:'',F12.6)') NUM_BLOCKED
          CALL WRITES(IOER,OP_STRING,ERROR,*9999)
        ELSE IF(ITYP3(NR,NX).EQ.4)THEN
          CALL CALC_RESIS_FLOW(1,NBJ,NEELEM,NENP,NORD,NPNE,NYNE,NVJE,
     '      XAB,XP,YP,ERROR,*9999)
        ENDIF
      ELSE IF(ITYP3(nr,nx).EQ.5)THEN
        UPDATE_VECTOR=.TRUE.
        UPDATE_NY=.TRUE.
      ENDIF !ITYP3

      CALL EXITS('PRE_SOL')
      RETURN

 9999 CALL ERRORS('PRE_SOL',ERROR)
      CALL EXITS('PRE_SOL')
      RETURN 1
      END



