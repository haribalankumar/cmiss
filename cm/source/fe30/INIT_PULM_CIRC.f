      SUBROUTINE INIT_PULM_CIRC(ISC_GK,ISR_GK,ISC_GKK,ISR_GKK,nb,NBJ,
     '  NEELEM,NENP,NHP,NHST,NPNE,NPNODE,nr,NYNE,
     '  NYNP,NYNR,nx,NXI,nxl,SPARSE_P,SPARSE_S,T,DYNAM1,
     '  FIX,LINEAR,UPDATE_VECTOR,ERROR,*)

C####  Subroutine: INIT_PULM_CIRC
C###   Description:
C###     INIT_PULM_CIRC sets up any arrays required for calculations
C###     in pulmonary capillary system.

C***   Created by KSB, October 2001.


      IMPLICIT NONE
      INCLUDE 'b12.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'marc00.cmn'
      INCLUDE 'solv00.cmn'

!     Parameter List
      INTEGER ISC_GK(NISC_GKM),ISR_GK(NISR_GKM),ISC_GKK(NISC_GKKM),
     &  ISR_GKK(NISR_GKKM),NEELEM(0:NE_R_M),nb,NBJ(NJM,NEM),
     &  NENP(NPM,0:NEPM),NHP(NPM),NHST,
     &  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M),nr,
     &  NYNE(NAM,NHM,0:NRCM,NCM,NEM),NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM),
     &  NYNR(0:NY_R_M,0:NRCM,NCM),nx,NXI(-NIM:NIM,0:NEIM,0:NEM),nxl,
     &  SPARSE_P,SPARSE_S
      REAL*8 T
      CHARACTER ERROR*(*)
      LOGICAL DYNAM1,FIX(NYM,NIYFIXM),LINEAR,UPDATE_VECTOR
!     Local Variables
c      INTEGER ny

      CALL ENTERS('INIT_PULM_CIRC',*9999)

      DYNAM1=.TRUE. !don't use GD ?
      LINEAR=.FALSE. !non-linear, iterate w.r.t hematocrit distribution
      UPDATE_VECTOR=.TRUE.
      NHST=NHP(NPNODE(1))
      nxl=nx !temporary, b/c Merryn's stuff uses nxl
      N_SOLN=0
      nb=NBJ(1,NEELEM(1))
C SEN The Cmiss solvers all use Compressed-Row storage
      SPARSE_P=1
      SPARSE_S=1
      IF(RESTART)THEN
        T=TSTART-DT
      ELSE
        T=TSTART
      ENDIF
      CALL CALC_SPARSE_GKK_1DTREE(NISC_GKKM,NISR_GKKM,ISC_GK,ISC_GKK,
     &  ISR_GK,ISR_GKK,NYT(1,1,nx),NYT(2,1,nx),NBJ,NDIAG(nx),NEELEM,
     &  NENP,NHST,NPNE,nr,nx,NXI,NYNE,NYNP,NYNR,NZZT(1,nr,nx),
     &  SPARSE_P,FIX,ERROR,*9999)
      
C... Sets up NONY, NYNO, CONY, CYNO, NOT correctly (only rows removed)         
c      CALL GLOBAL_LUNG(NONY,NPNY,nr,nxl,NYNE,NYNO,NYNP,NYNR,CONY,CYNO,
c     &  FIX,ERROR,*9999)
      
      CALL EXITS('INIT_PULM_CIRC ')
      RETURN
 9999 CALL ERRORS('INIT_PULM_CIRC ',ERROR)
      CALL EXITS('INIT_PULM_CIRC ')
      RETURN 1
      END

