      SUBROUTINE BMG_PCG_SET_parms(
     &     NDIM,MAXIT,PRECON,NPRECON,IFD,RELAX,SYMRLX,TOL, 
     &     MG_CYCLE,MG_OUTPUT,PERIODIC1,PERIODIC2,PERIODIC3,INDEF,
     &     BMG_iPARMS,BMG_rPARMS,BMG_IOFLAG,BMG_PCG_OUTPUT,
     &     BMG_PCG_iPARMS,BMG_PCG_rPARMS
     &     )

      IMPLICIT NONE

C -----------------------------
C     Includes
C -----------------------------

      INCLUDE 'BMG_constants.h'
      INCLUDE 'BMG_parameters.h'
      INCLUDE 'BMG_PCG_parameters.h'

C ----------------------------
C     Argument Declarations
C ----------------------------

      INTEGER IFD, MAXIT, PRECON, NPRECON, NDIM, NX, 
     &        RELAX, MG_CYCLE, MG_OUTPUT, SYMRLX
      LOGICAL PERIODIC1, PERIODIC2, PERIODIC3, INDEF
      REAL*8  TOL
      !
      ! BoxMG Cycle and I/O Parameters
      !
      INTEGER  BMG_iPARMS(NBMG_iPARMS), BMG_PCG_iPARMS(NBMG_PCG_iPARMS)
      REAL*8   BMG_rPARMS(NBMG_rPARMS), BMG_PCG_rPARMS(NBMG_PCG_rPARMS)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG), BMG_PCG_OUTPUT
      !
      ! Local Variables
      !
      INTEGER I

C ==========================================================================

      DO I=1,NBMG_iPARMS
         BMG_iPARMS(I) = 0
      ENDDO

      DO I=1,NBMG_rPARMS
         BMG_rPARMS(I) = 0.D0
      ENDDO

      DO I=1,NBMG_IOFLAG
         BMG_IOFLAG(I) = .FALSE.
      ENDDO

      DO I=1,NBMG_PCG_iPARMS
         BMG_PCG_iPARMS(I) = 0
      ENDDO
      
      DO I=1,NBMG_PCG_rPARMS
         BMG_PCG_rPARMS(I) = 0.D0
      ENDDO

      BMG_PCG_OUTPUT = .FALSE.

      IF( NDIM.EQ.2 ) THEN

         !
         !  Default BOXMG parameters
         !
         CALL BMG2_SymStd_SETUP_parms( 
     &        BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG 
     &        ) 

         !
         !  Set BOXMG preconditioner parameters 
         !
         IF( IFD.NE.1 ) THEN
            BMG_iPARMS(id_BMG2_STENCIL) = BMG_STENCIL_9pt
         ELSE
            BMG_iPARMS(id_BMG2_STENCIL) = BMG_STENCIL_5pt
         ENDIF
         BMG_iPARMS(id_BMG2_SETUP)        = BMG_SETUP_ptrs_opers
         BMG_iPARMS(id_BMG2_CG_CONSTRUCT) = BMG_CG_CONS_block
         !
         ! Set BC
         !
         IF( INDEF ) THEN
           IF( PERIODIC1 ) THEN
             BMG_iPARMS(id_BMG2_BC) = BMG_BCs_indef_per_x
           ELSEIF( PERIODIC2 ) THEN
             BMG_iPARMS(id_BMG2_BC) = BMG_BCs_indef_per_y
           ELSE
             BMG_iPARMS(id_BMG2_BC) = BMG_BCs_indef_nonper
           ENDIF
         ELSE
           IF( PERIODIC1 ) THEN
             BMG_iPARMS(id_BMG2_BC) = BMG_BCs_def_per_x
           ELSEIF( PERIODIC2 ) THEN
             BMG_iPARMS(id_BMG2_BC) = BMG_BCs_def_per_y
           ELSE
             BMG_iPARMS(id_BMG2_BC) = BMG_BCs_definite
           ENDIF
         ENDIF

         IF( PRECON.EQ.3 ) THEN
            BMG_iPARMS(id_BMG2_RELAX)       = RELAX
            BMG_iPARMS(id_BMG2_RELAX_SYM)   = BMG_RELAX_SYM
            BMG_iPARMS(id_BMG2_NCYCLE_TYPE) = MG_CYCLE
            BMG_iPARMS(id_BMG2_MAX_ITERS)   = NPRECON
         ENDIF

         IF( MG_OUTPUT.GE.1 ) THEN
            BMG_PCG_OUTPUT  =  .TRUE.
         ENDIF

         IF( MG_OUTPUT.GE.2 ) THEN
            BMG_IOFLAG(iBMG2_OUT_ITERATIONS)   = .TRUE.
         ENDIF

         IF( MG_OUTPUT.GE.3 ) THEN
            BMG_IOFLAG(iBMG2_BUG_RES_CG_SOLVE) = .TRUE.
            BMG_IOFLAG(iBMG2_BUG_RES_RELAX)    = .TRUE.
         ENDIF

         IF( MG_OUTPUT.GE.4 ) THEN
            BMG_IOFLAG(iBMG2_OUT_TIME_CYCLING) = .TRUE.
         ENDIF

      ELSEIF( NDIM.EQ.3 ) THEN

         !
         !  Default BOXMG parameters
         !
         CALL BMG3_SymStd_SETUP_parms( 
     &        BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG 
     &        )  

         !
         !  Set BOXMG preconditioner parameters
         !
         IF( IFD.NE.1 ) THEN
            BMG_iPARMS(id_BMG3_STENCIL) = BMG_STENCIL_27pt
         ELSE
            BMG_iPARMS(id_BMG3_STENCIL) = BMG_STENCIL_7pt
         ENDIF
         !
         BMG_iPARMS(id_BMG3_SETUP)        =  BMG_SETUP_ptrs_opers
         !
         BMG_iPARMS(id_BMG3_CG_TYPE)      =  BMG_CG_ITLI
         !
         BMG_iPARMS(id_BMG3_CG_CONSTRUCT) =  BMG_CG_CONS_block
         !
         ! Set Default BC
         !
         IF( INDEF ) THEN
           IF( PERIODIC1 ) THEN
             BMG_iPARMS(id_BMG3_BC) = BMG_BCs_indef_per_x
           ELSEIF( PERIODIC2 ) THEN
             BMG_iPARMS(id_BMG3_BC) = BMG_BCs_indef_per_y
           ELSEIF( PERIODIC3 ) THEN
             BMG_iPARMS(id_BMG3_BC) = BMG_BCs_indef_per_z
           ELSE
             BMG_iPARMS(id_BMG3_BC) = BMG_BCs_indef_nonper
           ENDIF
         ELSE
           IF( PERIODIC1 ) THEN
             BMG_iPARMS(id_BMG2_BC) = BMG_BCs_def_per_x
           ELSEIF( PERIODIC2 ) THEN
             BMG_iPARMS(id_BMG2_BC) = BMG_BCs_def_per_y
           ELSEIF( PERIODIC3 ) THEN
             BMG_iPARMS(id_BMG3_BC) = BMG_BCs_def_per_z
           ELSE
             BMG_iPARMS(id_BMG2_BC) = BMG_BCs_definite
           ENDIF
         ENDIF
         !
         IF( PRECON.EQ.3 ) THEN
            BMG_iPARMS(id_BMG3_RELAX)       = RELAX
            BMG_iPARMS(id_BMG3_RELAX_SYM)   = BMG_RELAX_SYM 
            BMG_iPARMS(id_BMG3_NCYCLE_TYPE) = BMG_V_CYCLE
            BMG_iPARMS(id_BMG3_MAX_ITERS)   = NPRECON
         ENDIF

         IF( MG_OUTPUT.GE.1 ) THEN
            BMG_PCG_OUTPUT  =  .TRUE.
         ENDIF

         IF( MG_OUTPUT.GE.2 ) THEN
            BMG_IOFLAG(iBMG3_OUT_ITERATIONS)   = .TRUE.
         ENDIF

         IF( MG_OUTPUT.GE.3 ) THEN
            BMG_IOFLAG(iBMG3_BUG_RES_CG_SOLVE) = .TRUE.
            BMG_IOFLAG(iBMG3_BUG_RES_RELAX)    = .TRUE.
         ENDIF

         IF( MG_OUTPUT.GE.4 ) THEN
            BMG_IOFLAG(iBMG3_OUT_TIME_CYCLING) = .TRUE.
         ENDIF
         
      ENDIF

      !
      !  Set PCG parameters
      !
      BMG_PCG_iPARMS(id_BMG_PCG_PRECON)     = PRECON

      BMG_PCG_iPARMS(id_BMG_PCG_STOP_TEST)  = BMG_PCG_STOP_REL_RES_L2
      BMG_PCG_iPARMS(id_BMG_PCG_MAX_ITERS)  = MAXIT
      BMG_PCG_rPARMS(id_BMG_PCG_STOP_TOL)   = TOL
         

      RETURN
      END
