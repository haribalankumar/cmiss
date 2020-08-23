      SUBROUTINE BMG_SET_parms(
     &     ND,MAXIT,IFD,RELAX,SYMRLX,TOL,MG_CYCLE,MG_OUTPUT, 
     &     PERIODIC1,PERIODIC2,PERIODIC3,INDEF,
     &     BMG_iPARMS,BMG_rPARMS,BMG_IOFLAG
     &     )

      IMPLICIT NONE

C -----------------------------
C     Includes
C -----------------------------

      INCLUDE 'BMG_constants.h'
      INCLUDE 'BMG_parameters.h'

C ----------------------------
C     Argument Declarations
C ----------------------------

      INTEGER IFD, MAXIT, ND, NX, RELAX, MG_CYCLE, MG_OUTPUT, SYMRLX
      LOGICAL INDEF, PERIODIC1, PERIODIC2, PERIODIC3
      REAL*8  TOL
      !
      ! BoxMG Cycle and I/O Parameters
      !
      INTEGER  BMG_iPARMS(NBMG_iPARMS)
      REAL*8   BMG_rPARMS(NBMG_rPARMS)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)
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
            
      IF( ND.EQ.2 ) THEN

         !
         !  Default BOXMG parameters
         !
         CALL BMG2_SymStd_SETUP_parms( 
     &        BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG 
     &        ) 

         !
         ! * CMISS DEVELOPER override default parameters *
         !
         IF( IFD.NE.1 ) THEN
            BMG_iPARMS(id_BMG2_STENCIL) = BMG_STENCIL_9pt
         ELSE
            BMG_iPARMS(id_BMG2_STENCIL) = BMG_STENCIL_5pt
         ENDIF
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
         !
         BMG_iPARMS(id_BMG2_STOP_TEST) = BMG_STOP_REL_RES_L2
         !
         ! * CMISS USER override default parameters *
         !
         BMG_iPARMS(id_BMG2_RELAX)       = RELAX
         BMG_iPARMS(id_BMG2_RELAX_SYM)   = SYMRLX

         IF( MG_CYCLE.EQ.1 ) THEN
            BMG_iPARMS(id_BMG2_CYCLE_CLASS) = BMG_N_CYCLE
            BMG_iPARMS(id_BMG2_NCYCLE_TYPE) = BMG_V_CYCLE
         ELSEIF( MG_CYCLE.EQ.3 ) THEN
            BMG_iPARMS(id_BMG2_CYCLE_CLASS) = BMG_FMG_CYCLE
            BMG_iPARMS(id_BMG2_NCYCLE_TYPE) = BMG_V_CYCLE
         ELSEIF( MG_CYCLE.EQ.4 ) THEN
            BMG_iPARMS(id_BMG2_CYCLE_CLASS) = BMG_N_CYCLE
            BMG_iPARMS(id_BMG2_NCYCLE_TYPE) = BMG_W_CYCLE
         ENDIF
         
         BMG_iPARMS(id_BMG2_MAX_ITERS) = MAXIT
         BMG_rPARMS(id_BMG2_STOP_TOL)  = TOL

         !
         ! * CMISS USER override default I/O parameters *
         !
         IF( MG_OUTPUT.EQ.1 ) THEN
            BMG_IOFLAG(iBMG2_OUT_ITERATIONS)   = .TRUE.
         ELSEIF( MG_OUTPUT.EQ.2 ) THEN
            BMG_IOFLAG(iBMG2_OUT_ITERATIONS)   = .TRUE.
            BMG_IOFLAG(iBMG2_BUG_RES_CG_SOLVE) = .TRUE.
            BMG_IOFLAG(iBMG2_BUG_RES_RELAX)    = .TRUE.
         ELSEIF( MG_OUTPUT.EQ.3 ) THEN
            BMG_IOFLAG(iBMG2_OUT_ITERATIONS)   = .TRUE.
            BMG_IOFLAG(iBMG2_BUG_RES_CG_SOLVE) = .TRUE.
            BMG_IOFLAG(iBMG2_BUG_RES_RELAX)    = .TRUE.
            BMG_IOFLAG(iBMG2_OUT_TIME_CYCLING) = .TRUE.
         ELSEIF( MG_OUTPUT.EQ.4 ) THEN
            BMG_IOFLAG(iBMG2_OUT_ITERATIONS)   = .TRUE.
            BMG_IOFLAG(iBMG2_BUG_RES_CG_SOLVE) = .TRUE.
            BMG_IOFLAG(iBMG2_BUG_RES_RELAX)    = .TRUE.
            BMG_IOFLAG(iBMG2_OUT_TIME_CYCLING) = .TRUE.
            BMG_IOFLAG(iBMG2_OUT_TIME_SETUP)   = .TRUE.
            BMG_IOFLAG(iBMG2_OUT_WSPACE_POINT) = .TRUE.
            BMG_IOFLAG(iBMG2_OUT_WSPACE_SIZE)  = .TRUE.
         ENDIF

      ELSEIF( ND.EQ.3 ) THEN

         !
         !  Default BOXMG parameters
         !
         CALL BMG3_SymStd_SETUP_parms( 
     &        BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG 
     &        )   

         !
         ! * CMISS DEVELOPER override default parameters *
         !
         IF( IFD.NE.1 ) THEN
            BMG_iPARMS(id_BMG3_STENCIL) = BMG_STENCIL_27pt
         ELSE
            BMG_iPARMS(id_BMG3_STENCIL) = BMG_STENCIL_7pt
         ENDIF
         !
         BMG_iPARMS(id_BMG3_CG_TYPE) =  BMG_CG_ITLI
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
         BMG_iPARMS(id_BMG3_STOP_TEST) = BMG_STOP_REL_RES_L2
         !
         BMG_iPARMS(id_BMG3_NRELAX_DOWN) = 1
         BMG_iPARMS(id_BMG3_NRELAX_UP)   = 1
         BMG_iPARMS(id_BMG2_MAX_ITERS)   = 1
         BMG_rPARMS(id_BMG2_STOP_TOL)    = 1e-4
         
         BMG_iPARMS(id_BMG3_CG_MIN_DIM) = 3

         !
         ! * CMISS USER override default parameters *
         !
         BMG_iPARMS(id_BMG3_RELAX)        = RELAX
         BMG_iPARMS(id_BMG3_RELAX_SYM)    = SYMRLX
         BMG_iPARMS(id_BMG2_RELAX_SYM)    = SYMRLX

         IF( MG_CYCLE.EQ.1 ) THEN
            BMG_iPARMS(id_BMG3_CYCLE_CLASS) = BMG_N_CYCLE
            BMG_iPARMS(id_BMG3_NCYCLE_TYPE) = BMG_V_CYCLE
         ELSEIF( MG_CYCLE.EQ.3 ) THEN
            BMG_iPARMS(id_BMG3_CYCLE_CLASS) = BMG_FMG_CYCLE
            BMG_iPARMS(id_BMG3_NCYCLE_TYPE) = BMG_V_CYCLE
         ELSEIF( MG_CYCLE.EQ.4 ) THEN
            BMG_iPARMS(id_BMG3_CYCLE_CLASS) = BMG_N_CYCLE
            BMG_iPARMS(id_BMG3_NCYCLE_TYPE) = BMG_W_CYCLE
         ENDIF

         BMG_iPARMS(id_BMG3_MAX_ITERS) = MAXIT
         BMG_rPARMS(id_BMG3_STOP_TOL)  = TOL


         !
         ! * CMISS USER override default I/O parameters *
         !
         IF( MG_OUTPUT.EQ.1 ) THEN
            BMG_IOFLAG(iBMG3_OUT_ITERATIONS)   = .TRUE.
         ELSEIF( MG_OUTPUT.EQ.2 ) THEN
            BMG_IOFLAG(iBMG3_OUT_ITERATIONS)   = .TRUE.
            BMG_IOFLAG(iBMG3_BUG_RES_CG_SOLVE) = .TRUE.
            BMG_IOFLAG(iBMG3_BUG_RES_RELAX)    = .TRUE.
         ELSEIF( MG_OUTPUT.EQ.3 ) THEN
            BMG_IOFLAG(iBMG3_OUT_ITERATIONS)   = .TRUE.
            BMG_IOFLAG(iBMG3_BUG_RES_CG_SOLVE) = .TRUE.
            BMG_IOFLAG(iBMG3_BUG_RES_RELAX)    = .TRUE.
            BMG_IOFLAG(iBMG3_OUT_TIME_CYCLING) = .TRUE.
         ELSEIF( MG_OUTPUT.EQ.4 ) THEN
            BMG_IOFLAG(iBMG3_OUT_ITERATIONS)   = .TRUE.
            BMG_IOFLAG(iBMG3_BUG_RES_CG_SOLVE) = .TRUE.
            BMG_IOFLAG(iBMG3_BUG_RES_RELAX)    = .TRUE.
            BMG_IOFLAG(iBMG3_OUT_TIME_CYCLING) = .TRUE.
            BMG_IOFLAG(iBMG3_OUT_TIME_SETUP)   = .TRUE.
            BMG_IOFLAG(iBMG3_OUT_WSPACE_POINT) = .TRUE.
            BMG_IOFLAG(iBMG3_OUT_WSPACE_SIZE)  = .TRUE.
         ENDIF

      ENDIF

      RETURN
      END
