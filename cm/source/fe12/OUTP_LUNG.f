      SUBROUTINE OUTP_LUNG(NBJ,NELIST,NENP,NPNE,nr,nx,NYNE,NYNP,
     &  CE,YP,ERROR,*)

C#### Subroutine:  OUTP_LUNG
C###  Description:
C###    OUTP_LUNG outputs solution variables for pulmonary transport
C###    problems.

C*** Written by KSB, 20th July, 2004.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'pulm00.cmn'
      INCLUDE 'lung00.cmn'
!     Parameter list
      INTEGER NBJ(NJM,NEM),NELIST(0:NEM),NENP(NPM,0:NEPM),
     &  NPNE(NNM,NBFM,NEM),nr,nx,
     &  NYNE(NAM,NHM,0:NRCM,NCM,NEM),NYNP(NKM,NVM,NHM,NPM)
      REAL*8 CE(NMM,NEM),YP(NYM,NIYM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER nb
      INTEGER*4 PATHWAY_PTR,PLOT_DATA_PTR,TIME_PTR

      CALL ENTERS('OUTP_LUNG',*9999)
      
      IF(ITYP3(nr,nx).EQ.3.OR.ITYP3(nr,nx).EQ.5.OR.
     '  ITYP3(nr,nx).EQ.6)THEN
 !pulmonary capillary blood flow
        PATHWAY_PTR=0
        PLOT_DATA_PTR=0
        TIME_PTR=0
        CALL ALLOCATE_MEMORY((MAX_PATH+1)*FACTORS,1,DPTYPE,PATHWAY_PTR,
     '    MEM_INIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(4500*FACTORS,1,DPTYPE,PLOT_DATA_PTR,
     '    MEM_INIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(NEM,1,DPTYPE,TIME_PTR,MEM_INIT,ERROR,
     &    *9999)
        nb=nbj(1,NELIST(1)) !basis function for 1st element
        CALL CALC_TRANSIT_TIME(nb,NELIST,NENP,NPNE,
     '    NYNE,NYNP,CE,%VAL(PATHWAY_PTR),%VAL(PLOT_DATA_PTR),
     &    %VAL(TIME_PTR),YP(1,1),ERROR,*9999)
         !calculates RBC and neutrophil transit times
        CALL FREE_MEMORY(PATHWAY_PTR,ERROR,*9999)
        CALL FREE_MEMORY(PLOT_DATA_PTR,ERROR,*9999)
        CALL FREE_MEMORY(TIME_PTR,ERROR,*9999)
C... outputs model solution data, geometric, pressures and flow
C... results for model analysis.        
        CALL OPPCAP(nb,NELIST,NENP,NPNE,NYNE,NYNP,CE,YP,ERROR,
     '    *9999)
      ENDIF
      
      CALL EXITS('OUTP_LUNG')
      RETURN
 9999 CALL ERRORS('OUTP_LUNG',ERROR)
      CALL EXITS('OUTP_LUNG')
      RETURN 1
      END
      
