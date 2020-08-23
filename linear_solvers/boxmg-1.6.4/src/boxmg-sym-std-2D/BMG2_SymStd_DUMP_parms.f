      SUBROUTINE BMG2_SymStd_DUMP_parms( 
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG 
     &                )

C ==========================================================================
C
C   BMG2_SymStd_DUMP_parms.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_SymStd_DUMP_parms.f outputs the parameter arrays.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2003/01/30 (JDM)
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
      INCLUDE 'BMG_workspace.h'

C ---------------------------
C     Argument Declarations:
C
      INTEGER BMG_iPARMS(NBMG_iPARMS)
      REAL*8  BMG_rPARMS(NBMG_rPARMS)
      LOGICAL BMG_IOFLAG(NBMG_IOFLAG)

C --------------------------
C     Local Declarations:
C

C ==========================================================================

      WRITE(*,*) '****  Parameters: 2D '
      WRITE(*,*) 
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_POINTERS)    = ',
     &           BMG_iPARMS(id_BMG2_POINTERS)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_STENCIL)     = ', 
     &           BMG_iPARMS(id_BMG2_STENCIL)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_SETUP)       = ',
     &           BMG_iPARMS(id_BMG2_SETUP)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_RELAX)       = ',
     &           BMG_iPARMS(id_BMG2_RELAX)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_RELAX_SYM)   = ',
     &           BMG_iPARMS(id_BMG2_RELAX_SYM )
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_NRELAX_DOWN) = ',
     &           BMG_iPARMS(id_BMG2_NRELAX_DOWN)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_NRELAX_UP)   = ',
     &           BMG_iPARMS(id_BMG2_NRELAX_UP)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_NRELAX_FG)   = ',
     &           BMG_iPARMS(id_BMG2_NRELAX_FG)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_CYCLE_CLASS) = ',
     &           BMG_iPARMS(id_BMG2_CYCLE_CLASS)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_NCYCLE_TYPE) = ',
     &           BMG_iPARMS(id_BMG2_NCYCLE_TYPE)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_FMG_NNCYCLE) = ',
     &           BMG_iPARMS(id_BMG2_FMG_NNCYCLE)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_MAX_ITERS)   = ', 
     &           BMG_iPARMS(id_BMG2_MAX_ITERS)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_STOP_TEST)   = ',
     &           BMG_iPARMS(id_BMG2_STOP_TEST)
      WRITE(*,*) 'BMG_rPARMS(id_BMG2_STOP_TOL)    = ',
     &           BMG_rPARMS(id_BMG2_STOP_TOL)
      WRITE(*,*) 'BMG_iPARMS(id_BMG2_CG_MIN_DIM)  = ', 
     &           BMG_iPARMS(id_BMG2_CG_MIN_DIM)


      RETURN
      END
