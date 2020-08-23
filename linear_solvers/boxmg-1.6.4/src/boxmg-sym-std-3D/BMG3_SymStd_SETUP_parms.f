      SUBROUTINE BMG3_SymStd_SETUP_parms( 
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAGS
     &                )

C ==========================================================================
C
C   BMG3_SymStd_SETUP_parms.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_SETUP_parms.f sets up default cycle parameters
C   which should work for most problems.  This includes setting
C   up the default parameters for plane relaxation.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2000/12/28 (JDM)
C
C ==================================================================
C   INPUT:
C ========================
C
C
C
C ==================================================================
C   OUTPUT:
C ===========================
C
C
C
C ==================================================================
C   LOCAL:
C ========================
C
C
C
C ==========================================================================

      IMPLICIT NONE

C -----------------------------
C     Includes
C
      INCLUDE 'BMG_constants.h'
      INCLUDE 'BMG_parameters.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER  BMG_iPARMS(NBMG_iPARMS)
      REAL*8   BMG_rPARMS(NBMG_rPARMS)
      LOGICAL  BMG_IOFLAGS(NBMG_IOFLAG)

C ----------------------------
C     Local Declarations
C
      INTEGER i

C ==========================================================================

C --------------------------------------------------------------------------
C     ********** SETUP 3D CYCLE
C --------------------------------------------------------------------------

C -------------------------------------
C    Memory allocation:
C -------------------------------------

      BMG_iPARMS(id_BMG3_POINTERS) = BMG_NO_pointers

C -------------------------------------
C    Stencil:
C -------------------------------------

      BMG_iPARMS(id_BMG3_STENCIL)  = BMG_STENCIL_7pt

C -------------------------------------
C     Setup:
C -------------------------------------

      BMG_iPARMS(id_BMG3_SETUP) = BMG_SETUP_Ptrs_Opers

C -------------------------------------
C     Interp:
C -------------------------------------

      BMG_iPARMS(id_BMG3_INTERP) = BMG_INTERP_OI

C -------------------------------------
C     Relaxation:
C -------------------------------------

      BMG_iPARMS(id_BMG3_RELAX)     
     &          = BMG_GS_RB_planes_xy_yz_xz

      BMG_iPARMS(id_BMG3_RELAX_SYM)    = BMG_RELAX_NONSYM

      BMG_iPARMS(id_BMG3_NRELAX_DOWN) = 1
      BMG_iPARMS(id_BMG3_NRELAX_UP)   = 1
      BMG_iPARMS(id_BMG3_NRELAX_FG)   = 0   ! not used at the moment

C -------------------------------------
C     Cycle Class and Type
C -------------------------------------

      BMG_iPARMS(id_BMG3_CYCLE_CLASS) = BMG_N_CYCLE
      BMG_iPARMS(id_BMG3_NCYCLE_TYPE) = BMG_W_CYCLE
      BMG_iPARMS(id_BMG3_FMG_NNCYCLE) = 1

C -------------------------------------
C     Stopping Criteria
C -------------------------------------

      BMG_iPARMS(id_BMG3_MAX_ITERS) = 10
      BMG_iPARMS(id_BMG3_STOP_TEST) = BMG_STOP_REL_RES_L2
      BMG_rPARMS(id_BMG3_STOP_TOL)  = 1D-6

C -------------------------------------
C     Coarsening limit
C -------------------------------------

      BMG_iPARMS(id_BMG3_CG_MIN_DIM) = 3

C -------------------------------------
C     Coarse-Grid Operator TYPE
C -------------------------------------

      BMG_iPARMS(id_BMG3_CG_TYPE) = BMG_CG_ITLI

C -------------------------------------
C     Coarse-Grid Operator CONSTRUCTION
C -------------------------------------

      BMG_iPARMS(id_BMG3_CG_CONSTRUCT) = BMG_CG_CONS_explicit

C -------------------------------------
C     Minimum number of grids
C -------------------------------------

      BMG_iPARMS(id_BMG3_MIN_NOG) = 2

C -------------------------------------
C     Set error parameters to zero
C -------------------------------------

      BMG_iPARMS(id_BMG3_Err_Code) = 0
      BMG_iPARMS(id_BMG3_Ext_Err_Code) = 0


C ==========================================================================

C --------------------------------------------------------------------------
C     ******** SETUP PLANE RELAXATION 
C --------------------------------------------------------------------------

C -------------------------------------
C     Setup:
C -------------------------------------

      BMG_iPARMS(id_BMG2_SETUP) = BMG_SETUP_none

C -------------------------------------
C     Interpolation:
C -------------------------------------

      BMG_iPARMS(id_BMG2_INTERP) = BMG_INTERP_OI

C -------------------------------------
C     Relaxation:
C -------------------------------------

      BMG_iPARMS(id_BMG2_RELAX)      = BMG_GS_RB_x_y_lines
      BMG_iPARMS(id_BMG2_RELAX_SYM ) = BMG_RELAX_NONSYM

      BMG_iPARMS(id_BMG2_NRELAX_DOWN ) = 1
      BMG_iPARMS(id_BMG2_NRELAX_UP )   = 1
      BMG_iPARMS(id_BMG2_NRELAX_FG )   = 0   ! not used at the moment

C -------------------------------------
C     Cycle Class and Type
C -------------------------------------

      BMG_iPARMS(id_BMG2_CYCLE_CLASS) = BMG_N_CYCLE
      BMG_iPARMS(id_BMG2_NCYCLE_TYPE) = BMG_V_CYCLE
      BMG_iPARMS(id_BMG2_FMG_NNCYCLE) = 1

C -------------------------------------
C     Stopping Criteria
C -------------------------------------

      BMG_iPARMS(id_BMG2_MAX_ITERS) = 1
      BMG_iPARMS(id_BMG2_STOP_TEST) = BMG_STOP_REL_RES_L2
      BMG_rPARMS(id_BMG2_STOP_TOL)  = rZERO

C -------------------------------------
C     Coarsening limit
C -------------------------------------

      BMG_iPARMS(id_BMG2_CG_MIN_DIM) = 3

C -------------------------------------
C     Coarse-Grid Operator TYPE
C -------------------------------------

      BMG_iPARMS(id_BMG2_CG_TYPE) = BMG_CG_ITLI

C -------------------------------------
C     Coarse-Grid Operator CONSTRUCTION
C -------------------------------------

      BMG_iPARMS(id_BMG2_CG_CONSTRUCT) = BMG_CG_CONS_explicit

C -------------------------------------
C     Minimum number of grids
C -------------------------------------

      BMG_iPARMS(id_BMG2_MIN_NOG) = 2

C ==========================================================================

C -------------------------------------
C     I/O Parameters
C -------------------------------------

      DO i=1, NBMG_IOFLAG
         BMG_IOFLAGS(i)=.FALSE.
      ENDDO

C ==========================================================================


C =====================================

      RETURN
      END
