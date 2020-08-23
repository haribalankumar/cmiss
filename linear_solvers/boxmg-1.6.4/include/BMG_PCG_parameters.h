C =======================================================================
C   
C   Include file BMG_PCG_parameters.h    
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   This include file provides a single resource for defining
C   commonly used parameters in the pcg code used with the BOXMG
C   family of codes.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   2003/07/01  - Written (T.Austin)
C
C =======================================================================

C =======================================================================
C ----------------------------------------------------
C     Parameter Indexing
C ---------------------------------

C ---------------------------------
C     INTEGER Parameters
C ---------------------------------

      INTEGER   id_BMG_PCG_NMG_CYCLES,
     &          id_BMG_PCG_STOP_TEST,
     &          id_BMG_PCG_PRECON,
     &          id_BMG_PCG_BMG_SETUP,
     &          id_BMG_PCG_MAX_ITERS,
     &          id_BMG_PCG_NUM_ITERS

      PARAMETER ( id_BMG_PCG_NMG_CYCLES   = 1,
     &            id_BMG_PCG_STOP_TEST    = 2,
     &            id_BMG_PCG_PRECON       = 3,
     &            id_BMG_PCG_BMG_SETUP    = 4,
     &            id_BMG_PCG_MAX_ITERS    = 5,
     &            id_BMG_PCG_NUM_ITERS    = 6  )

      INTEGER   NBMG_PCG_iPARMS 
      PARAMETER ( NBMG_PCG_iPARMS = 6 )


C -------------------------------
C     REAL Parameters
C -------------------------------

      INTEGER   id_BMG_PCG_STOP_TOL,
     &          id_BMG_PCG_FINAL_TOL  
      PARAMETER ( id_BMG_PCG_STOP_TOL = 1,
     &            id_BMG_PCG_FINAL_TOL = 2 )

      INTEGER   NBMG_PCG_rPARMS 
      PARAMETER ( NBMG_PCG_rPARMS = 2 )

C --------------------------------
C     Stopping Critieria
C --------------------------------

      INTEGER BMG_PCG_STOP_ABS_RES_L2,
     &        BMG_PCG_STOP_REL_RES_L2,
     &        BMG_PCG_STOP_ABS_RES_M2,
     &        BMG_PCG_STOP_REL_RES_M2

      PARAMETER( BMG_PCG_STOP_ABS_RES_L2 = 0,
     &           BMG_PCG_STOP_REL_RES_L2 = 1,
     &           BMG_PCG_STOP_ABS_RES_M2 = 2,
     &           BMG_PCG_STOP_REL_RES_M2 = 3   ) 

C --------------------------------
C     Preconditioners
C --------------------------------

      INTEGER BMG_PCG_PRECON_NONE,
     &        BMG_PCG_PRECON_DIAG,
     &        BMG_PCG_PRECON_BMG

      PARAMETER( BMG_PCG_PRECON_NONE  = 1,
     &           BMG_PCG_PRECON_DIAG  = 2,
     &           BMG_PCG_PRECON_BMG   = 3  ) 

C -------------------------------
C     BMG-PCG Setup options:
C -------------------------------

      INTEGER BMG_PCG_BMG_SETUP_all,
     &        BMG_PCG_BMG_SETUP_only,
     &        BMG_PCG_BMG_SETUP_none

      PARAMETER ( BMG_PCG_BMG_SETUP_all  = 0,
     &            BMG_PCG_BMG_SETUP_only = 1,
     &            BMG_PCG_BMG_SETUP_none = 2 )

C =======================================================================

