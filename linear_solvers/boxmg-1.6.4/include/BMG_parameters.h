C =======================================================================
C   
C   Include file constants.h    
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   This include file provides a single resource for defining
C   commonly used parameters in the BOXMG family of codes.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   2000/02/22  - Written (M.Berndt)
C   2000/03/06  - added 2D relaxation and stopping criteria parameters
C   2000/03/07  - added IO/Debugging parameters (JDM)
C   2000/12/09  - added Cycle Parameters        (JDM)
C   2002/04/26  - added pointer parameters      (JDM)
C
C =======================================================================

C =======================================================================
C ----------------------------------------------------
C     Parameter Indexing
C ---------------------------------

C ---------------------------------
C     INTEGER Parameters
C ---------------------------------

      INTEGER   id_BMG2_DIM_NOG,
     &          id_BMG2_DIM_NF,
     &          id_BMG2_DIM_NC,
     &          id_BMG2_DIM_NSO,
     &          id_BMG2_DIM_NSOR,
     &          id_BMG2_DIM_NCI,
     &          id_BMG2_DIM_NCBW,
     &          id_BMG2_DIM_NCU, 
     &          id_BMG2_POINTERS,
     &          id_BMG2_STENCIL,
     &          id_BMG2_BC,
     &          id_BMG2_SETUP,
     &          id_BMG2_INTERP,
     &          id_BMG2_RELAX,
     &          id_BMG2_RELAX_SYM,
     &          id_BMG2_NRELAX_DOWN,
     &          id_BMG2_NRELAX_UP,
     &          id_BMG2_NRELAX_FG,
     &          id_BMG2_CYCLE_CLASS,
     &          id_BMG2_NCYCLE_TYPE,
     &          id_BMG2_FMG_NNCYCLE,
     &          id_BMG2_MAX_ITERS,
     &          id_BMG2_NUM_ITERS,
     &          id_BMG2_STOP_TEST,
     &          id_BMG2_MIN_NOG,
     &          id_BMG2_CG_MIN_DIM,
     &          id_BMG2_CG_TYPE,
     &          id_BMG2_CG_CONSTRUCT,
     &          id_BMG2_Err_Code,
     &          id_BMG2_Ext_Err_Code,
     &          id_BMG3_DIM_NOG,
     &          id_BMG3_DIM_NF,
     &          id_BMG3_DIM_NC,
     &          id_BMG3_DIM_NSO,
     &          id_BMG3_DIM_NSOR,
     &          id_BMG3_DIM_NCI,
     &          id_BMG3_DIM_NCBW,
     &          id_BMG3_DIM_NCU, 
     &          id_BMG3_DIM_NSO3, 
     &          id_BMG3_DIM_NSR3, 
     &          id_BMG3_DIM_NCI3, 
     &          id_BMG3_DIM_NCPL, 
     &          id_BMG3_POINTERS,
     &          id_BMG3_STENCIL,
     &          id_BMG3_BC,
     &          id_BMG3_SETUP,
     &          id_BMG3_INTERP,
     &          id_BMG3_RELAX,
     &          id_BMG3_RELAX_SYM,
     &          id_BMG3_NRELAX_DOWN,
     &          id_BMG3_NRELAX_UP,
     &          id_BMG3_NRELAX_FG,
     &          id_BMG3_ADAPT_TYPE,
     &          id_BMG3_CYCLE_CLASS,
     &          id_BMG3_NCYCLE_TYPE,
     &          id_BMG3_FMG_NNCYCLE,
     &          id_BMG3_MAX_ITERS,
     &          id_BMG3_NUM_ITERS,
     &          id_BMG3_STOP_TEST,
     &          id_BMG3_MIN_NOG,
     &          id_BMG3_CG_MIN_DIM,
     &          id_BMG3_CG_TYPE,
     &          id_BMG3_CG_CONSTRUCT,
     &          id_BMG3_Err_Code,
     &          id_BMG3_Ext_Err_Code,
     &          NBMG2_iPARMS, NBMG3_iPARMS, NBMG_iPARMS 

      PARAMETER ( id_BMG2_DIM_NOG      =  1,
     &            id_BMG2_DIM_NF       =  2,
     &            id_BMG2_DIM_NC       =  3,
     &            id_BMG2_DIM_NSO      =  4,
     &            id_BMG2_DIM_NSOR     =  5,
     &            id_BMG2_DIM_NCI      =  6,
     &            id_BMG2_DIM_NCBW     =  7,
     &            id_BMG2_DIM_NCU      =  8, 
     &            id_BMG2_POINTERS     =  9,
     &            id_BMG2_STENCIL      =  10,
     &            id_BMG2_BC           =  11, 
     &            id_BMG2_SETUP        =  12,
     &            id_BMG2_INTERP       =  13,
     &            id_BMG2_RELAX        =  14,
     &            id_BMG2_RELAX_SYM    =  15,
     &            id_BMG2_NRELAX_DOWN  =  16,
     &            id_BMG2_NRELAX_UP    =  17,
     &            id_BMG2_NRELAX_FG    =  18,
     &            id_BMG2_CYCLE_CLASS  =  19,
     &            id_BMG2_NCYCLE_TYPE  =  20,
     &            id_BMG2_FMG_NNCYCLE  =  21,
     &            id_BMG2_MAX_ITERS    =  22,
     &            id_BMG2_NUM_ITERS    =  23,
     &            id_BMG2_STOP_TEST    =  24,    
     &            id_BMG2_MIN_NOG      =  25,
     &            id_BMG2_CG_MIN_DIM   =  26,
     &            id_BMG2_CG_CONSTRUCT =  27,
     &            id_BMG2_CG_TYPE      =  28,
     &            id_BMG2_Err_Code     =  29,
     &            id_BMG2_Ext_Err_Code =  30,
     &            id_BMG3_DIM_NOG      =  31,
     &            id_BMG3_DIM_NF       =  32,
     &            id_BMG3_DIM_NC       =  33,
     &            id_BMG3_DIM_NSO      =  34,
     &            id_BMG3_DIM_NSOR     =  35,
     &            id_BMG3_DIM_NCI      =  36,
     &            id_BMG3_DIM_NCBW     =  37,
     &            id_BMG3_DIM_NCU      =  38, 
     &            id_BMG3_DIM_NSO3     =  39, 
     &            id_BMG3_DIM_NSR3     =  40, 
     &            id_BMG3_DIM_NCI3     =  41, 
     &            id_BMG3_DIM_NCPL     =  42, 
     &            id_BMG3_POINTERS     =  43,
     &            id_BMG3_STENCIL      =  44,
     &            id_BMG3_BC           =  45,
     &            id_BMG3_SETUP        =  46,
     &            id_BMG3_INTERP       =  47,
     &            id_BMG3_RELAX        =  48,
     &            id_BMG3_RELAX_SYM    =  49,
     &            id_BMG3_NRELAX_DOWN  =  50,
     &            id_BMG3_NRELAX_UP    =  51,
     &            id_BMG3_NRELAX_FG    =  52,
     &            id_BMG3_ADAPT_TYPE   =  53,
     &            id_BMG3_CYCLE_CLASS  =  54,
     &            id_BMG3_NCYCLE_TYPE  =  55,
     &            id_BMG3_FMG_NNCYCLE  =  56,
     &            id_BMG3_MAX_ITERS    =  57,
     &            id_BMG3_NUM_ITERS    =  58,
     &            id_BMG3_STOP_TEST    =  59,
     &            id_BMG3_MIN_NOG      =  60,
     &            id_BMG3_CG_MIN_DIM   =  61,
     &            id_BMG3_CG_TYPE      =  62,
     &            id_BMG3_CG_CONSTRUCT =  63,
     &            id_BMG3_Err_Code     =  64,
     &            id_BMG3_Ext_Err_Code =  65,
     &            NBMG2_iPARMS = 30,      ! Number of 2D parameters
     &            NBMG3_iPARMS = 35,      ! Number of 3D parameters
     &            NBMG_iPARMS  = 65       ! Dimension of BMG_iPARM
     &            )

C -------------------------------
C     REAL Parameters
C -------------------------------

      INTEGER   id_BMG2_STOP_TOL,
     &          id_BMG2_FINAL_TOL,
     &          id_BMG2_ANORM,
     &          id_BMG3_STOP_TOL,
     &          id_BMG3_FINAL_TOL,
     &          id_BMG3_ANORM,
     &          NBMG_rPARMS 

      PARAMETER ( id_BMG2_STOP_TOL  = 1,
     &            id_BMG2_FINAL_TOL = 2,
     &            id_BMG2_ANORM     = 3,
     &            id_BMG3_STOP_TOL  = 4,
     &            id_BMG3_FINAL_TOL = 5,
     &            id_BMG3_ANORM     = 6,
     &            NBMG_rPARMS = 6       )


C ----------------------------------------------------
C ============================================================================
C ----------------------------------------------------
C     Parameter Values
C -------------------------------

C -------------------------------
C     Stencil Type
C -------------------------------

      INTEGER BMG_STENCIL_5pt,
     &        BMG_STENCIL_9pt,
     &        BMG_STENCIL_7pt,
     &        BMG_STENCIL_27pt

      PARAMETER ( BMG_STENCIL_5pt  = 1,
     &            BMG_STENCIL_9pt  = 2,
     &            BMG_STENCIL_7pt  = 1,
     &            BMG_STENCIL_27pt = 2  )


C -------------------------------
C     Periodicity
C -------------------------------

      INTEGER  BMG_BCs_definite,
     &         BMG_BCs_def_per_x,
     &         BMG_BCs_def_per_y,
     &         BMG_BCs_def_per_z,
     &         BMG_BCs_def_per_xy,
     &         BMG_BCs_indef_per_x,
     &         BMG_BCs_indef_per_y,
     &         BMG_BCs_indef_per_z,
     &         BMG_BCs_indef_per_xy,
     &         BMG_BCs_indef_nonper

      PARAMETER ( BMG_BCs_definite     =  0,
     &            BMG_BCs_def_per_x    =  1,
     &            BMG_BCs_def_per_y    =  2,
     &            BMG_BCs_def_per_z    =  3,
     &            BMG_BCs_def_per_xy   =  4,
     &            BMG_BCs_indef_per_x  = -1,
     &            BMG_BCs_indef_per_y  = -2,
     &            BMG_BCs_indef_per_z  = -3,
     &            BMG_BCs_indef_per_xy = -4,
     &            BMG_BCs_indef_nonper = -5  )

C -------------------------------
C     Setup options:
C -------------------------------

      INTEGER BMG_SETUP_only,
     &        BMG_SETUP_none,
     &        BMG_SETUP_opers, 
     &        BMG_SETUP_ptrs_opers 

      PARAMETER ( BMG_SETUP_only  = 3,
     &            BMG_SETUP_none  = 2,
     &            BMG_SETUP_opers = 1,
     &            BMG_SETUP_ptrs_opers = 0 )

C -------------------------------
C     Memory Allocation:
C -------------------------------

      INTEGER BMG_USE_pointers,
     &        BMG_NO_pointers

      PARAMETER ( BMG_USE_pointers = 1,
     &            BMG_NO_pointers = 0   )

C -------------------------------
C     Interpolation:
C -------------------------------

      INTEGER  BMG_INTERP_OI,
     &         BMG_INTERP_GMG

      PARAMETER( BMG_INTERP_OI   = 1,
     &           BMG_INTERP_GMG  = 2 )

C -------------------------------
C     Relaxation:
C -------------------------------
      
      INTEGER  BMG_GS_RB_point, 
     &         BMG_GS_RB_x_lines,
     &         BMG_GS_RB_y_lines,
     &         BMG_GS_RB_x_y_lines,
     &         BMG_GS_RB_planes_xy_yz_xz,
     &         BMG_GS_adaptive,
     &         BMG_CG_relaxation

      PARAMETER( BMG_GS_RB_point           = 1,
     &           BMG_GS_RB_x_lines         = 2,
     &           BMG_GS_RB_y_lines         = 3,
     &           BMG_GS_RB_x_y_lines       = 4,
     &           BMG_GS_RB_planes_xy_yz_xz = 5,
     &           BMG_GS_adaptive           = 6,
     &           BMG_CG_relaxation         = 7 )


C --------------------------------
C     Symmetry of the MG n-cycle:
C --------------------------------

      INTEGER  BMG_RELAX_NONSYM,
     &         BMG_RELAX_SYM 
      PARAMETER( BMG_RELAX_NONSYM = 0, 
     &           BMG_RELAX_SYM = 1    )

      INTEGER  BMG_DOWN, 
     &         BMG_UP
      PARAMETER(  BMG_DOWN = 0, 
     &            BMG_UP = 1   )

C --------------------------------
C     Stopping Critieria
C --------------------------------

      INTEGER BMG_STOP_ABS_RES_L2,
     &        BMG_STOP_REL_RES_L2

      PARAMETER( BMG_STOP_ABS_RES_L2 = 0, 
     &           BMG_STOP_REL_RES_L2 = 1 ) 

C --------------------------------
C     Cycle Class and Type
C --------------------------------
      
      INTEGER   BMG_FMG_CYCLE, 
     &          BMG_N_CYCLE, 
     &          BMG_V_CYCLE, 
     &          BMG_W_CYCLE

      PARAMETER ( BMG_FMG_CYCLE = 0,
     &            BMG_N_CYCLE   = 1,
     &            BMG_V_CYCLE   = 1,
     &            BMG_W_CYCLE   = 2 )


C --------------------------------
C     Coarse-Grid Operator Type
C --------------------------------

      INTEGER   BMG_CG_ITLI_IzIyIx,
     &          BMG_CG_ITLI,
     &          BMG_CG_USER

      PARAMETER ( BMG_CG_ITLI_IzIyIx = 1,
     &            BMG_CG_ITLI        = 2,
     &            BMG_CG_USER        = 3 )


C --------------------------------
C     I^{T} L I Construction
C --------------------------------

      INTEGER   BMG_CG_CONS_explicit,
     &          BMG_CG_CONS_block

      PARAMETER ( BMG_CG_CONS_explicit = 1,
     &            BMG_CG_CONS_block    = 2  )
                

C ----------------------------------------------------
C ============================================================================
C ----------------------------------------------------
C     IOFLAG Indexing
C ---------------------------------

      INTEGER  iBMG2_BUG_STENCIL_FG,
     &         iBMG2_BUG_STENCIL_CG,
     &         iBMG2_BUG_STENCIL_CG1,
     &         iBMG2_BUG_RESTRICT,
     &         iBMG2_BUG_INTERP,
     &         iBMG2_BUG_RES_CG_SOLVE,         
     &         iBMG2_BUG_RES_INTERP,
     &         iBMG2_BUG_RES_RELAX,
     &         iBMG2_BUG_RES_RESTRICT,
     &         iBMG2_BUG_PARAMETERS,
     &         iBMG2_OUT_ITERATIONS,
     &         iBMG2_OUT_STENCIL_TTY,
     &         iBMG2_OUT_RESTRICT_TTY,
     &         iBMG2_OUT_INTERP_TTY,
     &         iBMG2_OUT_TIME_CYCLING,
     &         iBMG2_OUT_TIME_SETUP,
     &         iBMG2_OUT_TIME_TOTAL,
     &         iBMG2_OUT_WSPACE_SIZE,
     &         iBMG2_OUT_WSPACE_POINT,
     &         iBMG2_WARN_ZERO_RESIDUAL,
     &         iBMG2_WARN_ZERO_RHS,
     &         iBMG2_OUT_STOP_ERROR,
     &         iBMG3_BUG_STENCIL_FG,
     &         iBMG3_BUG_STENCIL_CG,
     &         iBMG3_BUG_STENCIL_CG1,
     &         iBMG3_BUG_RESTRICT,
     &         iBMG3_BUG_INTERP,
     &         iBMG3_BUG_RES_CG_SOLVE,         
     &         iBMG3_BUG_RES_INTERP,
     &         iBMG3_BUG_RES_RELAX,         
     &         iBMG3_BUG_RES_RESTRICT,
     &         iBMG3_BUG_PARAMETERS,
     &         iBMG3_OUT_ITERATIONS,
     &         iBMG3_OUT_STENCIL_TTY,
     &         iBMG3_OUT_RESTRICT_TTY,
     &         iBMG3_OUT_INTERP_TTY,
     &         iBMG3_OUT_TIME_CYCLING,
     &         iBMG3_OUT_TIME_SETUP,
     &         iBMG3_OUT_TIME_TOTAL,
     &         iBMG3_OUT_WSPACE_SIZE,
     &         iBMG3_OUT_WSPACE_POINT,
     &         iBMG3_WARN_ZERO_RESIDUAL,
     &         iBMG3_WARN_ZERO_RHS,
     &         iBMG3_OUT_STOP_ERROR,
     &         iBMG_OUT_RHS,
     &         iBMG_OUT_SOLUTION,
     &         NBMG2_IOFLAG, NBMG3_IOFLAG, NBMG_IOFLAG

      PARAMETER( iBMG2_BUG_STENCIL_FG     = 1,
     &           iBMG2_BUG_STENCIL_CG     = 2, 
     &           iBMG2_BUG_STENCIL_CG1    = 3,
     &           iBMG2_BUG_RESTRICT       = 4,
     &           iBMG2_BUG_INTERP         = 5,
     &           iBMG2_BUG_RES_INTERP     = 6,
     &           iBMG2_BUG_RES_RESTRICT   = 7,
     &           iBMG2_BUG_RES_RELAX      = 8,
     &           iBMG2_BUG_RES_CG_SOLVE   = 9,
     &           iBMG2_BUG_PARAMETERS     = 10,
     &           iBMG2_OUT_WSPACE_SIZE    = 11,
     &           iBMG2_OUT_WSPACE_POINT   = 12,
     &           iBMG2_OUT_TIME_SETUP     = 13,
     &           iBMG2_OUT_TIME_CYCLING   = 14,
     &           iBMG2_OUT_TIME_TOTAL     = 15,
     &           iBMG2_OUT_ITERATIONS     = 16,
     &           iBMG2_OUT_STENCIL_TTY    = 17,
     &           iBMG2_OUT_RESTRICT_TTY   = 18,
     &           iBMG2_OUT_INTERP_TTY     = 19,
     &           iBMG2_WARN_ZERO_RESIDUAL = 20,
     &           iBMG2_WARN_ZERO_RHS      = 21,
     &           iBMG2_OUT_STOP_ERROR     = 22,
     &           iBMG3_BUG_STENCIL_FG     = 23,
     &           iBMG3_BUG_STENCIL_CG     = 24, 
     &           iBMG3_BUG_STENCIL_CG1    = 25,
     &           iBMG3_BUG_RESTRICT       = 26,
     &           iBMG3_BUG_INTERP         = 27,
     &           iBMG3_BUG_RES_INTERP     = 28,
     &           iBMG3_BUG_RES_RESTRICT   = 29,
     &           iBMG3_BUG_RES_RELAX      = 30,
     &           iBMG3_BUG_RES_CG_SOLVE   = 31,
     &           iBMG3_BUG_PARAMETERS     = 32,
     &           iBMG3_OUT_WSPACE_SIZE    = 33,
     &           iBMG3_OUT_WSPACE_POINT   = 34,
     &           iBMG3_OUT_TIME_SETUP     = 35,
     &           iBMG3_OUT_TIME_CYCLING   = 36,
     &           iBMG3_OUT_TIME_TOTAL     = 37,
     &           iBMG3_OUT_ITERATIONS     = 38,
     &           iBMG3_OUT_STENCIL_TTY    = 39,
     &           iBMG3_OUT_RESTRICT_TTY   = 40,
     &           iBMG3_OUT_INTERP_TTY     = 41,
     &           iBMG3_WARN_ZERO_RESIDUAL = 42,
     &           iBMG3_WARN_ZERO_RHS      = 23,
     &           iBMG3_OUT_STOP_ERROR     = 44,
     &           iBMG_OUT_SOLUTION        = 45, 
     &           iBMG_OUT_RHS             = 46,
     &           NBMG2_IOFLAG = 23,       ! Number of 2D I/O FLAGS
     &           NBMG3_IOFLAG = 23,       ! Number of 3D I/O FLAGS
     &           NBMG_IOFLAG  = 46        ! Dimension of BMG_IOFLAG
     &           )

C ----------------------------------------------------
C =========================================================================



