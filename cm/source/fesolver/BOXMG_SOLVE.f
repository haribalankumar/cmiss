      SUBROUTINE BOXMG_SOLVE(
     &              METHOD,ISC_A,ISR_A,LDA,M,N,NZA,
     &              OUTPUTCODE,SOLVER,SPARSE_A,A,B,X,FIRST_A,
     &              UPDATE_A,X_INIT,NEELEM,NQNE,NQS,NQXI,NXI,
     &              NR,NX,ERROR,*)

C#### Subroutine: BOXMG_SOLVE
C###  Description:
C###    BOXMG_SOLVE solves a linear system of equations Ax=b that is
C###    generated from a logically rectangular (2D or 3D) problem. In
C###    the context of CMISS we are concerned with problems consisting
C###    of one element and a collocation, finite volume, or finite element
C###    discretization.  The solver used is determined by SOLVER and the 
C###    sparsity used is determined by SPARSENESS.  FIRST_A and UPDATE_A 
C###    should be .TRUE. for the first call.  On return FIRST_A will be 
C###    set to .FALSE..   On subsequent calls UPDATE_A should be 
C###    .TRUE./.FALSE. if refactorization is necessary/unnessary, and FIRST_A 
C###    should be set to .TRUE. again if the sparsity pattern changes.

      IMPLICIT NONE
 
      INCLUDE 'BMG_parameters.h'
      INCLUDE 'BMG_PCG_parameters.h'
      INCLUDE 'BMG_workspace.h'

      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'bmg00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'ptr00.cmn'
      INCLUDE 'solver.inc'
      INCLUDE 'solv00.cmn'
      INCLUDE 'time02.cmn'
      INCLUDE 'grid00.cmn'

!     Parameter List
      INTEGER LDA,M,N,NR,NX,NZA
      INTEGER SPARSE_A,SOLVER,OUTPUTCODE,METHOD
      INTEGER ISC_A(*),ISR_A(*)
      INTEGER NQXI(0:NIM,NQSCM),NXI(-NIM:NIM,0:NEIM,0:NEM),
     &        NQNE(NEQM,NQEM),NQS(NEQM),NEELEM(0:NE_R_M,0:NRM)
      REAL *8 A(*)
      REAL *8 B(N),X(N)
      CHARACTER ERROR*(*)
      LOGICAL FIRST_A,UPDATE_A,X_INIT

!     Local Variables

      INTEGER I,ITER,NDIM,NEXI1,NEXI2,NEXI3,NXI1,NXI2,NXI3
      REAL TIME1,TIME2,TIME3,FACT_TIME,SOLN_TIME,TOTL_TIME
      REAL*8 RESID,TOL
      LOGICAL DENSE,FACTOR,INIT,PCG_BOXMG,MG_BOXMG,
     &        PERIODIC1,PERIODIC2,PERIODIC3,INDEF

!     BMG Local Variables
      INTEGER BMG_BC,BMG_NOGm,BMG_NFm,BMG_NSOm,
     &        IFD,MAPSZ,MAXIT,MGOUT,NPRECON,NStncl,PRECON,
     &        PSYCLE,RELAX,SYMRLX
      INTEGER pSI,pSR

      INTEGER  BMG_IJK_MAP_SIZE
      EXTERNAL BMG_IJK_MAP_SIZE  

C ===============================================
C     BEGIN EXECUTION.
C ===============================================

      CALL ENTERS('BOXMG_SOLVE',*9999)

      CALL ASSERT(nx.LE.9,'>>Can only have nx.le.9  at present',
     &  ERROR,*9999)

      DENSE   = SPARSE_A.EQ.0  

      PCG_BOXMG = SOLVER.EQ.SOLV_PCG_BMG
      MG_BOXMG  = SOLVER.EQ.SOLV_BMG

      FACTOR  = UPDATE_A.OR.FIRST_A

      MGOUT   = MG_OUTPUT(NX)
      RELAX   = MG_RELAX_TYPE(NX)
      PRECON  = PRECON_CODE(NX)
      NPRECON = SOLV_NPRECON(NX)
      PSYCLE  = MG_CYCLE(NX)
      SYMRLX  = MG_RELAX_SYM(NX)

      TOL     = SOLV_TOL(NX)
      MAXIT   = SOLV_MAXIT(NX)

      !****************
      ! CAUSING ERRORS
      !****************     
      NDIM = NQXI(0,1)

      IF(SPARSE_A.NE.0 .AND. SPARSE_A.NE.1) THEN
        ERROR='>>Only able to solve Dense and Compressed-Row systems'
        GOTO 9999
      ENDIF

C     Optionally print out array stats
      IF(FACTOR .AND. OUTPUTCODE.GE.2) THEN
        CALL CHECK_MATRIX(A,LDA,N,SPARSE_A,ISC_A,ISR_A,ERROR,*9999)
      ENDIF

C     Optionally fix the system in case we have some rows that are zero
      IF(SOLV_FIXIT(NX)) THEN
        CALL FIX_MATRIX(A,LDA,N,SPARSE_A,ISC_A,ISR_A,ERROR,*9999)
      ENDIF

C ===============================================
C     BOXMG Specific Setup and Memory Allocation
C ===============================================

      CALL CPU_TIMER(CPU_USER,TIME1)

      IF(FIRST_A.OR.(.NOT.BMG_ALLOCATED(NX)) ) THEN

C        Set parameters

         IF(NDIM.EQ.2) THEN
            IF(METHOD.EQ.7) THEN
               IFD=1
               NStncl=3
            ELSE
               IFD=2
               NStncl=5
            ENDIF
         ELSEIF(NDIM.EQ.3) THEN
             IF(METHOD.EQ.7) THEN
                IFD=1
                NStncl=4
             ELSE
                IFD=2
                NStncl=14
             ENDIF
         ENDIF

         !
         ! Setup CM_MESH_TO_BMG and in the process get NXI1, NXI2, NXI3
         !
         CALL SETUP_CM_MESH_TO_BMG(
     &        NDIM,NXI1,NXI2,NXI3,NEELEM,NEXI1,NEXI2,NEXI3,
     &        NQNE,NQS,NQXI,NR,NX,NXI,PERIODIC1,PERIODIC2,
     &        PERIODIC3,ERROR,*9999
     &        )    

         BMG_DIMS(1,NX) = NXI1
         BMG_DIMS(2,NX) = NXI2
         BMG_DIMS(3,NX) = NXI3

         CALL CHECK_INDEFINITE(INDEF,A,ISR_A,ISC_A,N)

         IF( MG_BOXMG ) THEN
            !
            ! Standard BOXMG Solve
            !
            CALL BMG_SET_parms(
     &          NDIM,MAXIT,IFD,RELAX,SYMRLX,TOL,PSYCLE,MGOUT,
     &          PERIODIC1,PERIODIC2,PERIODIC3,INDEF,BMG_iPARMS(1,NX),
     &          BMG_rPARMS(1,NX),BMG_IOFLAG(1,NX)
     &          )

         ELSEIF( PCG_BOXMG ) THEN
            !
            ! Preconditioned CG using BOXMG preconditioner
            !
            CALL BMG_PCG_SET_parms(
     &          NDIM,MAXIT,PRECON,NPRECON,IFD,RELAX,SYMRLX,TOL,
     &          PSYCLE,MGOUT,PERIODIC1,PERIODIC2,PERIODIC3,INDEF,
     &          BMG_iPARMS(1,NX),BMG_rPARMS(1,NX),
     &          BMG_IOFLAG(1,NX),BMG_PCG_OUTPUT(NX),
     &          BMG_PCG_iPARMS(1,NX),BMG_PCG_rPARMS(1,NX)
     &          )

         ENDIF

C        Set what to store in workspace array 
      
         DO i=1, NBMG_InWORK
            BMG_InWORK(i,NX) = .FALSE.
         ENDDO

C        Store SO, U, Q, and RES all in BMG_rWORK array
         BMG_InWORK(i_InWORK_SO,NX)    = .TRUE.        
         BMG_InWORK(i_InWORK_U,NX)     = .TRUE. 
         BMG_InWORK(i_InWORK_Q,NX)     = .TRUE. 
         BMG_InWORK(i_InWORK_RES,NX)   = .TRUE. 

         IF( PCG_BOXMG ) THEN
            BMG_InWORK(i_InWORK_PCG_P,NX) = .TRUE. 
            BMG_InWORK(i_InWORK_PCG_R,NX) = .TRUE. 
            BMG_InWORK(i_InWORK_PCG_Z,NX) = .TRUE. 
         ENDIF

C        Later put all PCG arrays into rWORK when using PCG 

         IF(NDIM.EQ.2) THEN
            BMG_IOFLAG(iBMG2_OUT_STOP_ERROR,NX) = .TRUE. 
            BMG_iPARMS(id_BMG2_POINTERS,NX)     = BMG_USE_pointers
         ELSEIF(NDIM.EQ.3) THEN
            BMG_IOFLAG(iBMG3_OUT_STOP_ERROR,NX) = .TRUE.
            BMG_iPARMS(id_BMG3_POINTERS,NX)     = BMG_USE_pointers
         ENDIF
      
         pSI=1
         pSR=1

         BMG_NOGm = 1
         BMG_NFm  = 1
         BMG_NSOm = 1

         CALL BMG_SymStd_SETUP_PtrWork( 
     &           NDIM,BMG_DIMS(1,NX),
     &           BMG_DIMS(2,NX),BMG_DIMS(3,NX),
     &           BMG_iPARMS(1,NX), 
     &           BMG_NOGm,BMG_NFm,BMG_NSOm, 
     &           NBMG_iWORK(NX),NBMG_rWORK(NX),
     &           NBMG_iWORK_PL(NX),NBMG_rWORK_PL(NX),
     &           BMG_pWORK(1,NX),BMG_InWORK(1,NX), 
     &           pSR,pSI,BMG_IOFLAG(1,NX)
     &           ) 

C        Free/Allocate BOXMG memory as necessary
         CALL BOXMG_FREE(NX,FIRST_A,ERROR,*9999) 
         CALL BOXMG_ALLOC(NX,ERROR,*9999)

         IF( PCG_BOXMG ) THEN
            IF( BMG_PCG_RES(NX).NE.0 ) THEN 
               CALL FREE_MEMORY(BMG_PCG_RES(NX),ERROR,*9999)
            ENDIF
            CALL ALLOCATE_MEMORY(
     &              MAXIT,1,DPTYPE,BMG_PCG_RES(NX),INIT,ERROR,*9999
     &              )
         ENDIF

      ENDIF

      IF( FACTOR.OR.(.NOT.BMG_ALLOCATED(NX)) ) THEN 
        
         CALL BMG_UTILS_iV_zero(%val(BMG_iWORK_PTR(NX)),NBMG_iWORK(NX))
         CALL BMG_UTILS_rV_zero(%val(BMG_rWORK_PTR(NX)),NBMG_rWORK(NX))

         CALL BMG_UTILS_iV_zero(
     &           %val(BMG_iWORK_PL_PTR(NX)),NBMG_iWORK_PL(NX)
     &           )
         CALL BMG_UTILS_rV_zero(
     &           %val(BMG_rWORK_PL_PTR(NX)),NBMG_rWORK_PL(NX)
     &           )

         INIT = .TRUE.
         IF(.NOT.BMG_ALLOCATED(NX)) THEN

           CALL ALLOCATE_MEMORY(
     &          BMG_DIMS(1,NX)*BMG_DIMS(2,NX)*BMG_DIMS(3,NX),
     &          1,INTTYPE,BMG_BDY_PTR(NX),INIT,ERROR,*9999
     &          )
         
           CALL BMG_UTILS_iV_zero(
     &          %val(BMG_BDY_PTR(NX)),
     &          BMG_DIMS(1,NX)*BMG_DIMS(2,NX)*BMG_DIMS(3,NX)
     &          )
         
           IF(NDIM.EQ.2) THEN
             
             CALL ALLOCATE_MEMORY(
     &            2*BMG_DIMS(1,NX)*BMG_DIMS(2,NX),1,INTTYPE,
     &            CM_IJK_PTR(NX),INIT,ERROR,*9999
     &            )      
             CALL BMG_UTILS_iV_zero(
     &            %val(CM_IJK_PTR(NX)),
     &            2*BMG_DIMS(1,NX)*BMG_DIMS(2,NX)
     &            )

           ELSEIF(NDIM.EQ.3) THEN
             
             CALL ALLOCATE_MEMORY(
     &            3*BMG_DIMS(1,NX)*BMG_DIMS(2,NX)*BMG_DIMS(3,NX),
     &            1,INTTYPE,CM_IJK_PTR(NX),INIT,ERROR,*9999
     &            )
             CALL BMG_UTILS_iV_zero(
     &            %val(CM_IJK_PTR(NX)),
     &            3*BMG_DIMS(1,NX)*BMG_DIMS(2,NX)*BMG_DIMS(3,NX)
     &            )
             
           ENDIF

         ENDIF

         CALL BMG_CONSTRUCT_CM_IJK_MAP(
     &           NDIM,BMG_DIMS(1,NX),BMG_DIMS(2,NX),BMG_DIMS(3,NX),
     &           %val(CM_IJK_PTR(NX))
     &           )

         CALL BMG_CONVERT_csr_matrix( 
     &           NDIM,
     &           BMG_DIMS(1,NX),BMG_DIMS(2,NX),BMG_DIMS(3,NX),
     &           BMG_iPARMS(1,NX),%val(BMG_rWORK_PTR(NX)), 
     &           BMG_pWORK(1,NX),%val(BMG_BDY_PTR(NX)),
     &           %val(BMG_TO_CM_PTR(NX)),
     &           A,ISR_A,ISC_A,NStncl
     &           )
         
           
         MAPSZ = BMG_IJK_MAP_SIZE(NDIM,
     &          BMG_DIMS(1,NX)+2,
     &          BMG_DIMS(2,NX)+2,
     &          BMG_DIMS(3,NX)+2
     &          )

         IF(.NOT.BMG_ALLOCATED(NX)) THEN
           !
           CALL ALLOCATE_MEMORY(
     &          MAPSZ,1,INTTYPE,BMG_IJK_PTR(NX),INIT,ERROR,*9999
     &          )   
           !
         ENDIF

         CALL BMG_UTILS_iV_zero(%val(BMG_IJK_PTR(NX)),MAPSZ)

         IF( MG_BOXMG ) THEN

            IF(NDIM.EQ.2) THEN
               BMG_iPARMS(id_BMG2_SETUP,NX) = BMG_SETUP_only
            ELSEIF(NDIM.EQ.3) THEN
               BMG_iPARMS(id_BMG3_SETUP,NX) = BMG_SETUP_only
            ENDIF
            
            !
            ! Create prolongation, restriction, and coarse-grid operators
            !
            CALL BMG_SETUP(
     &              NDIM,BMG_DIMS(1,NX),BMG_DIMS(2,NX),BMG_DIMS(3,NX),
     &              BMG_iPARMS(1,NX),BMG_rPARMS(1,NX),BMG_IOFLAG(1,NX),
     &              %val(BMG_rWORK_PTR(NX)),%val(BMG_iWORK_PTR(NX)),
     &              BMG_pWORK(1,NX),%val(BMG_iWORK_PL_PTR(NX)),
     &              NBMG_iWORK_PL(NX),%val(BMG_rWORK_PL_PTR(NX)),
     &              NBMG_rWORK_PL(NX),%val(BMG_IJK_PTR(NX))
     &              )

         ELSEIF( PCG_BOXMG ) THEN
            
            BMG_PCG_iPARMS(id_BMG_PCG_BMG_SETUP,NX) = 
     &                                    BMG_PCG_BMG_SETUP_only
            
            !
            ! Create prolongation, restriction, and coarse-grid operators
            !
            CALL BMG_PCG_SETUP(
     &           NDIM,BMG_DIMS(1,NX),BMG_DIMS(2,NX),BMG_DIMS(3,NX),
     &           BMG_PCG_iPARMS(1,NX), BMG_PCG_rPARMS(1,NX), 
     &           BMG_PCG_OUTPUT(NX), %val(BMG_PCG_RES(NX)), MAXIT, 
     &           BMG_iPARMS(1,NX),BMG_rPARMS(1,NX),BMG_IOFLAG(1,NX),
     &           %val(BMG_rWORK_PTR(NX)),%val(BMG_iWORK_PTR(NX)),
     &           NBMG_rWORK(NX),NBMG_iWORK(NX),BMG_pWORK(1,NX),
     &           %val(BMG_iWORK_PL_PTR(NX)),NBMG_iWORK_PL(NX),
     &           %val(BMG_rWORK_PL_PTR(NX)),NBMG_rWORK_PL(NX),
     &           %val(BMG_IJK_PTR(NX))
     &           )
           
         ENDIF

         BMG_ALLOCATED(NX) = .TRUE.

      ENDIF

c     Write CMISS rhs into BOXMG rhs
      IF(NDIM.EQ.2) THEN
         BMG_BC = BMG_iPARMS(id_BMG2_BC,NX) 
      ELSEIF(NDIM.EQ.3) THEN
         BMG_BC = BMG_iPARMS(id_BMG3_BC,NX) 
      ENDIF

      CALL BMG_CONVERT_rhs(
     &        NDIM,BMG_DIMS(1,NX),BMG_DIMS(2,NX),BMG_DIMS(3,NX), 
     &        %val(BMG_rWORK_PTR(NX)),BMG_pWORK(1,NX),
     &        BMG_iPARMS(1,NX),%val(BMG_BDY_PTR(NX)),
     &        %val(CM_IJK_PTR(NX)),%val(BMG_TO_CM_PTR(NX)),
     &        A,ISR_A,ISC_A,B
     &        )

c     Set boundary conditions for solution vector
      CALL SOLV_FIX_BOUNDS(A,N,N,X,B,SPARSE_A,ISC_A,ISR_A,ERROR,*9999)

c     Write CMISS initial guess into BOXMG initial guess
      CALL BMG_CONVERT_soln(
     &        0,NDIM,BMG_DIMS(1,NX),BMG_DIMS(2,NX),BMG_DIMS(3,NX),
     &        %val(CM_IJK_PTR(NX)),%val(BMG_TO_CM_PTR(NX)),
     &        %val(BMG_rWORK_PTR(NX)),BMG_pWORK(1,NX),X
     &        )

      CALL CPU_TIMER(CPU_USER,TIME2)

      IF(NDIM.EQ.2) THEN
         BMG_iPARMS(id_BMG2_SETUP,NX) = BMG_SETUP_none
      ELSEIF(NDIM.EQ.3) THEN
         BMG_iPARMS(id_BMG3_SETUP,NX) = BMG_SETUP_none
      ENDIF

      IF( MG_BOXMG ) THEN

         IF(NDIM.EQ.2) THEN 
            BMG_iPARMS(id_BMG2_NUM_ITERS,NX) = 0
            BMG_rPARMS(id_BMG2_FINAL_TOL,NX) = 0.0D0
         ELSEIF(NDIM.EQ.3) THEN 
            BMG_iPARMS(id_BMG3_NUM_ITERS,NX) = 0
            BMG_rPARMS(id_BMG3_FINAL_TOL,NX) = 0.0D0 
         ENDIF

         CALL BMG_SOLVE(
     &           NDIM,BMG_DIMS(1,NX),BMG_DIMS(2,NX),BMG_DIMS(3,NX),
     &           BMG_iPARMS(1,NX),BMG_rPARMS(1,NX),BMG_IOFLAG(1,NX),
     &           %val(BMG_rWORK_PTR(NX)),%val(BMG_iWORK_PTR(NX)),
     &           BMG_pWORK(1,NX),%val(BMG_iWORK_PL_PTR(NX)),
     &           NBMG_iWORK_PL(NX),%val(BMG_rWORK_PL_PTR(NX)),
     &           NBMG_rWORK_PL(NX),%val(BMG_IJK_PTR(NX))
     &           )

         IF(NDIM.EQ.2) THEN 
            ITER  = BMG_iPARMS(id_BMG2_NUM_ITERS,NX)
            RESID = BMG_rPARMS(id_BMG2_FINAL_TOL,NX) 
         ELSEIF(NDIM.EQ.3) THEN 
            ITER  = BMG_iPARMS(id_BMG3_NUM_ITERS,NX) 
            RESID = BMG_rPARMS(id_BMG3_FINAL_TOL,NX) 
         ENDIF

      ELSEIF( PCG_BOXMG ) THEN

         BMG_PCG_iPARMS(id_BMG_PCG_NUM_ITERS,NX) = 0
         BMG_PCG_rPARMS(id_BMG_PCG_FINAL_TOL,NX) = 0.0D0

         CALL BMG_PCG_SOLVE(
     &           NDIM,BMG_DIMS(1,NX),BMG_DIMS(2,NX),BMG_DIMS(3,NX),
     &           BMG_PCG_iPARMS(1,NX),BMG_PCG_rPARMS(1,NX),
     &           BMG_PCG_OUTPUT(NX),%val(BMG_PCG_RES(NX)),MAXIT, 
     &           BMG_iPARMS(1,NX),BMG_rPARMS(1,NX),BMG_IOFLAG(1,NX),
     &           %val(BMG_rWORK_PTR(NX)),%val(BMG_iWORK_PTR(NX)),
     &           NBMG_rWORK(NX),NBMG_iWORK(NX),BMG_pWORK(1,NX),
     &           %val(BMG_iWORK_PL_PTR(NX)),NBMG_iWORK_PL(NX),
     &           %val(BMG_rWORK_PL_PTR(NX)),NBMG_rWORK_PL(NX),
     &           %val(BMG_IJK_PTR(NX))
     &           )

         ITER  = BMG_PCG_iPARMS(id_BMG_PCG_NUM_ITERS,NX)
         RESID = BMG_PCG_rPARMS(id_BMG_PCG_FINAL_TOL,NX)

      ENDIF 

      CALL BMG_CONVERT_soln( 
     &        1,NDIM,BMG_DIMS(1,NX),BMG_DIMS(2,NX),BMG_DIMS(3,NX),
     &        %val(CM_IJK_PTR(NX)),%val(BMG_TO_CM_PTR(NX)),
     &        %val(BMG_rWORK_PTR(NX)),BMG_pWORK(1,NX),X 
     &        )


      CALL CPU_TIMER(CPU_USER,TIME3) 

      IF(OUTPUTCODE.GE.1) THEN
         WRITE(OP_STRING,7000) ' Number of iterations: ',ITER
         CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
         WRITE(OP_STRING,7100) ' Residual of solution: ',RESID
         CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

C     Calculate and print the residual, and array/RHS norms etc
      IF(OUTPUTCODE.GE.1) THEN
         CALL CHECK_SOLN(A,LDA,N,B,X,SPARSE_A,ISC_A,ISR_A,ERROR,*9999)
      ENDIF

C     Timing information
      IF(OUTPUTCODE.GE.1) THEN
         FACT_TIME=TIME2-TIME1
         SOLN_TIME=TIME3-TIME2
         TOTL_TIME=TIME3-TIME1

         IF(FACTOR) THEN
            WRITE(OP_STRING,7200) ' CPU time to factorise system:     ', 
     &      FACT_TIME,'s'
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
         ENDIF 

         WRITE(OP_STRING,7200) ' CPU time to solve system:         ', 
     &        SOLN_TIME,'s'
         CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

         IF(FACTOR) THEN
            WRITE(OP_STRING,7200) ' Sum of factor and solution times: ', 
     &           TOTL_TIME,'s' 
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999) 
         ENDIF 
      ENDIF 

      CALL EXITS('BOXMG_SOLVE') 
      RETURN 

 405  FORMAT (1X,A,1X,1P,E16.9,/)
 406  FORMAT (/,1X,A,1X,1P,E16.9)
 410  FORMAT (/,1X,A,1X,/)
 2000 FORMAT(I7,2X,I7)
 2010 FORMAT(I7,2X,E15.9)
 7000 FORMAT(A,I10) 
 7100 FORMAT(A,G12.6) 
 7200 FORMAT(A,1P,D11.4,1X,A) 

 9999 CALL ERRORS('BOXMG_SOLVE',ERROR) 
      CALL EXITS('BOXMG_SOLVE') 
      RETURN 1 
      END 


      SUBROUTINE SETUP_CM_MESH_TO_BMG(
     &        NDIM,NXI1,NXI2,NXI3,NEELEM,NEXI1,NEXI2,NEXI3,
     &        NQNE,NQS,NQXI,NR,NX,NXI,PERIODIC1,PERIODIC2,
     E        PERIODIC3,ERROR,*)

C#### Subroutine: SETUP_CM_MESH_TO_BMG
C###  Description:
C###    SETUP_CM_MESH_TO_BMG sets up a map from a tensor product mesh
C###    in cmiss to an nxi1 X nxi2 X nxi3 logically rectangular grid for
C###    BOXMG.

      IMPLICIT NONE

      INCLUDE 'BMG_parameters.h'
      INCLUDE 'BMG_PCG_parameters.h'
      INCLUDE 'BMG_workspace.h'

      INCLUDE 'bmg00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'grid00.cmn'

!     Parameter List
      INTEGER NDIM,NEELEM(0:NE_R_M,0:NRM),NEXI1,NEXI2,NEXI3,
     &        NQNE(NEQM,NQEM),NQS(NEQM),NQXI(0:NIM,NQSCM),NR,
     &        NXI(-NIM:NIM,0:NEIM,0:NEM),NX,NXI1,NXI2,NXI3
      CHARACTER ERROR*(*)

!     Local Variables
      INTEGER CELEM
      INTEGER NTE,NEXTELEM

      LOGICAL PERIODIC1, PERIODIC2, PERIODIC3
      
C ===============================================
C     BEGIN EXECUTION.
C ===============================================

      CALL ENTERS('SETUP_CM_MESH_TO_BMG',*9999)

      NTE = NEELEM(0,nr)  ! Number Total Elements

      IF( NTE.EQ.1 ) THEN

         NDIM = NQXI(0,1)

         NEXI1 = 1
         NEXI2 = 1
         NEXI3 = 1

         NXI1 = NQXI(1,1)
         NXI2 = NQXI(2,1)
         NXI3 = NQXI(3,1)

         PERIODIC1 = .FALSE.
         PERIODIC2 = .FALSE.
         PERIODIC3 = .FALSE.

      ELSE

         PERIODIC1 = .FALSE.
         PERIODIC2 = .FALSE.
         PERIODIC3 = .FALSE.

         NDIM = NQXI(0,1)

         ! 
         ! Determine Number of Elements in each direction
         ! 
          
         ! xi1 direction
         
         NEXI1 = 1
         
         IF( NXI(1,1,NEELEM(1,nr)).NE.0 ) THEN

            CELEM    = 1 
            NEXTELEM = NXI(1,1,NEELEM(1,nr))
            DO WHILE( NEXTELEM.NE.0 .AND. NEXTELEM.NE.CELEM ) 

               NEXI1 = NEXI1 + 1
               NEXTELEM = NXI(1,1,NEXTELEM)

               IF( NEXTELEM.EQ.CELEM ) PERIODIC1 = .TRUE.

            ENDDO

         ENDIF

         IF( NXI(-1,1,NEELEM(1,nr)).NE.0 .AND. .NOT.PERIODIC1 ) THEN

            NEXTELEM = NXI(-1,1,NEELEM(1,nr))
            DO WHILE( NEXTELEM.NE.0 .AND. .NOT.PERIODIC1 ) 

               NEXI1 = NEXI1 + 1
               NEXTELEM = NXI(-1,1,NEXTELEM)

            ENDDO

         ENDIF

         ! xi2 direction

         NEXI2 = 1
         IF( NXI(2,1,NEELEM(1,nr)).NE.0 ) THEN

            CELEM    = 1
            NEXTELEM = NXI(2,1,NEELEM(1,nr))
            DO WHILE( NEXTELEM.NE.0 .AND. NEXTELEM.NE.CELEM ) 

               NEXI2 = NEXI2 + 1
               NEXTELEM = NXI(2,1,NEXTELEM)

               IF( NEXTELEM.EQ.CELEM ) PERIODIC2 = .TRUE.

            ENDDO

         ENDIF

         IF( NXI(-2,1,NEELEM(1,nr)).NE.0 .AND. .NOT.PERIODIC2 ) THEN

            NEXTELEM = NXI(-2,1,NEELEM(1,nr))
            DO WHILE( NEXTELEM.NE.0 ) 

               NEXI2 = NEXI2 + 1
               NEXTELEM = NXI(-2,1,NEXTELEM)

            ENDDO

         ENDIF

         ! xi3 direction

         NEXI3 = 1
         IF( NXI(3,1,NEELEM(1,nr)).NE.0 ) THEN

            CELEM    = 1
            NEXTELEM = NXI(3,1,NEELEM(1,nr))
            DO WHILE( NEXTELEM.NE.0 .AND. NEXTELEM.NE.CELEM ) 

               NEXI3 = NEXI3 + 1
               NEXTELEM = NXI(3,1,NEXTELEM)

               IF( NEXTELEM.EQ.CELEM ) PERIODIC3 = .TRUE.

            ENDDO

         ENDIF

         IF( NXI(-3,1,NEELEM(1,nr)).NE.0 .AND. .NOT.PERIODIC3 ) THEN

            NEXTELEM = NXI(-3,1,NEELEM(1,nr))
            DO WHILE( NEXTELEM.NE.0 ) 

               NEXI3 = NEXI3 + 1
               NEXTELEM = NXI(-3,1,NEXTELEM)

            ENDDO

         ENDIF

      ENDIF

      !
      ! Create map of elements
      ! 
      CALL ALLOCATE_MEMORY(
     &     NEXI1*NEXI2*NEXI3,1,INTTYPE,CM_MESH_DIST_PTR(NX),
     &     .TRUE.,ERROR,*9999
     &     )      

      ! 
      ! Setup CM_MESH_DIST
      !
      CALL SETUP_CM_MESH_DIST(
     &     %val(CM_MESH_DIST_PTR(NX)),NEELEM,NEXI1,NEXI2,NEXI3,NTE,
     &     NQXI,NXI,PERIODIC1,PERIODIC2,PERIODIC3
     &     )

      !
      ! Get sizes in three principal directions
      !     
      CALL GET_TOTAL_MESH_SIZES(
     &     %val(CM_MESH_DIST_PTR(NX)),NEXI1,NEXI2,NEXI3,NQS,NTE,NQXI,
     &     NXI1,NXI2,NXI3,PERIODIC1,PERIODIC2,PERIODIC3
     &     )

      CALL ALLOCATE_MEMORY(
     &     NXI1*NXI2*NXI3,1,INTTYPE,BMG_TO_CM_PTR(NX),
     &     .TRUE.,ERROR,*9999
     &     )      

      ! 
      ! Setup CM_MESH_TO_BMG
      !
      CALL FILLUP_CM_MESH_TO_BMG(
     &     %val(CM_MESH_DIST_PTR(NX)),%val(BMG_TO_CM_PTR(NX)),
     &     NEXI1,NEXI2,NEXI3,NQNE,NQS,NQXI,NXI1,NXI2,NXI3,NXI
     &     )
      

      CALL EXITS('SETUP_CM_MESH_TO_BMG') 
      RETURN 
      
 9999 CALL ERRORS('SETUP_CM_MESH_TO_BMG',ERROR) 
      CALL EXITS('SETUP_CM_MESH_TO_BMG') 
      RETURN 1 
      END 

      SUBROUTINE FILLUP_CM_MESH_TO_BMG(
     &     CM_MESH_DIST,CM_MESH_TO_BMG,NEXI1,NEXI2,NEXI3,
     &     NQNE,NQS,NQXI,NXI1,NXI2,NXI3,NXI
     &     )

      IMPLICIT NONE

      INCLUDE 'bmg00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'

!     Parameter List
      INTEGER NEXI1,NEXI2,NEXI3,NXI1,NXI2,NXI3,NQNE(NEQM,NQEM),
     &        NQS(NEQM),NQXI(0:NIM,NQSCM),NXI(-NIM:NIM,0:NEIM,0:NEM)

      INTEGER CM_MESH_TO_BMG(NXI1,NXI2,NXI3)
      INTEGER CM_MESH_DIST(NEXI1,NEXI2,NEXI3)
      
!     Local Variables
      INTEGER ELEM,I1,I2,J1,J2,K1,K2
      INTEGER EINDEX,INDEX1,INDEX2,INDEX3,NNXI1,NNXI2,NNXI3
      INTEGER ISTART,JSTART,KSTART

C ===============================================
C     BEGIN EXECUTION.
C ===============================================

      DO K2=1,NEXI3

         IF( K2.EQ.1 ) THEN
            KSTART = 1 
         ELSE
            KSTART = KSTART + NQXI(3,NQS(CM_MESH_DIST(1,1,K2-1))) - 1
         ENDIF

         DO J2=1,NEXI2

            IF( J2.EQ.1 ) THEN
               JSTART = 1 
            ELSE
               JSTART = JSTART + NQXI(2,NQS(CM_MESH_DIST(1,J2-1,1))) - 1
            ENDIF

            DO I2=1,NEXI1

               IF( I2.EQ.1 ) THEN
                  ISTART = 1 
               ELSE
                  ISTART = ISTART+NQXI(1,NQS(CM_MESH_DIST(I2-1,1,1)))-1
               ENDIF
            
               ELEM = CM_MESH_DIST(I2,J2,K2)

               NNXI1 = NQXI(1,NQS(ELEM))
               NNXI2 = NQXI(2,NQS(ELEM))
               NNXI3 = NQXI(3,NQS(ELEM))

               DO K1=1,NNXI3
               DO J1=1,NNXI2
               DO I1=1,NNXI1

                  EINDEX = (K1-1)*(NNXI1*NNXI2) + (J1-1)*NNXI1 + I1

                  INDEX1 = ISTART + I1 - 1 
                  INDEX2 = JSTART + J1 - 1
                  INDEX3 = KSTART + K1 - 1

                  CM_MESH_TO_BMG(INDEX1,INDEX2,INDEX3) =  
     &                 NQNE(ELEM,EINDEX)

               ENDDO
               ENDDO
               ENDDO
            
            ENDDO
         ENDDO
      ENDDO

      RETURN 
      END 


      SUBROUTINE SETUP_CM_MESH_DIST(
     &     CM_MESH_DIST,NEELEM,NEXI1,NEXI2,NEXI3,NTE,NQXI,NXI,
     &     PERIODIC1,PERIODIC2,PERIODIC3
     &     )

      IMPLICIT NONE

      INCLUDE 'bmg00.cmn'
      INCLUDE 'geom00.cmn'

!     Parameter List
      INTEGER NEELEM(0:NE_R_M,0:NRM),NEXI1,NEXI2,NEXI3,
     &        NQXI(0:NIM,NQSCM),NXI(-NIM:NIM,0:NEIM,0:NEM), NTE

      LOGICAL PERIODIC1,PERIODIC2,PERIODIC3

      INTEGER CM_MESH_DIST(NEXI1,NEXI2,NEXI3)
      
!     Local Variables
      INTEGER ELM,NXTELM,I,J,K,INITIAL1,INITIAL2,INITIAL3,
     &        INDEX1,INDEX2,INDEX3,NXTELM1,NXTELM2,NXTELM3

      LOGICAL COND1,COND2,COND3,ELMFOUND

C ===============================================
C     BEGIN EXECUTION.
C ===============================================

      !
      ! Get first element in "corner"
      !
      ELMFOUND = .FALSE.
      !
      ELM = 0
      !
      DO WHILE ( .NOT.ELMFOUND  ) 
         !
         ELM = ELM + 1
         !
         IF( (NXI(-1,1,NEELEM(ELM,1)).EQ.0 .OR. PERIODIC1) .AND.
     &       (NXI(-2,1,NEELEM(ELM,1)).EQ.0 .OR. PERIODIC2) .AND.
     &       (NXI(-3,1,NEELEM(ELM,1)).EQ.0 .OR. PERIODIC3)) THEN
               
           CM_MESH_DIST(1,1,1) = NEELEM(ELM,1)
           ELMFOUND = .TRUE.

        ENDIF

      ENDDO 

      COND3  = .TRUE.
      INDEX3 = 1

      INITIAL3 = CM_MESH_DIST(1,1,1)      
      NXTELM3 = INITIAL3

      DO WHILE ( NXTELM3.NE.0 .AND. 
     &          (NXTELM3.NE.INITIAL3.OR.COND3) )
           
         COND2  = .TRUE.
         INDEX2 = 1

         INITIAL2 = NXTELM3
         NXTELM2 = INITIAL2

         DO WHILE ( NXTELM2.NE.0 .AND. 
     &             (NXTELM2.NE.INITIAL2.OR.COND2) )

            COND1  = .TRUE.
            INDEX1 = 1

            INITIAL1 = NXTELM2
            NXTELM1 = INITIAL1

            DO WHILE ( NXTELM1.NE.0 .AND. 
     &                (NXTELM1.NE.INITIAL1.OR.COND1) )
               !
               CM_MESH_DIST(INDEX1,INDEX2,INDEX3) = NXTELM1
               !
               NXTELM1 = NXI(1,1,NXTELM1)
               COND1   = .FALSE.
               INDEX1  = INDEX1 + 1
               !
            ENDDO

            NXTELM2 = NXI(2,1,NXTELM2)
            INDEX2 = INDEX2 + 1

            COND2 = .FALSE.

         ENDDO

         NXTELM3 = NXI(3,1,NXTELM3)
         INDEX3 = INDEX3 + 1

         COND3 = .FALSE.
               
      ENDDO

      RETURN 
      END 

      SUBROUTINE GET_TOTAL_MESH_SIZES(
     &     CM_MESH_DIST,NEXI1,NEXI2,NEXI3,NQS,NTE,NQXI,
     &     NXI1,NXI2,NXI3,PERIODIC1,PERIODIC2,PERIODIC3
     &     )

      IMPLICIT NONE

      INCLUDE 'bmg00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'

!     Parameter List
      INTEGER NEXI1,NEXI2,NEXI3,NQS(NEQM),NQXI(0:NIM,NQSCM),NTE,
     &        NXI1,NXI2,NXI3

      INTEGER CM_MESH_DIST(NEXI1,NEXI2,NEXI3)
      LOGICAL PERIODIC1,PERIODIC2,PERIODIC3 
      
      INTEGER I,J,K

C ===============================================
C     BEGIN EXECUTION.
C ===============================================

      NXI1 = 0
      DO I=1,NEXI1
         NXI1 = NXI1 + NQXI(1,NQS(CM_MESH_DIST(I,1,1)))
      ENDDO
      NXI1 = NXI1 - NEXI1 
      IF( .NOT.PERIODIC1 ) NXI1 = NXI1 + 1

      NXI2 = 0
      DO J=1,NEXI2
         NXI2 = NXI2 + NQXI(2,NQS(CM_MESH_DIST(1,J,1)))
      ENDDO
      NXI2 = NXI2 - NEXI2
      IF( .NOT.PERIODIC2 ) NXI2 = NXI2 + 1

      NXI3 = 0
      DO K=1,NEXI3
         NXI3 = NXI3 + NQXI(3,NQS(CM_MESH_DIST(1,1,K)))
      ENDDO
      NXI3 = NXI3 - NEXI3 
      IF( .NOT.PERIODIC3 ) NXI3 = NXI3 + 1

      RETURN 
      END 

      SUBROUTINE CHECK_INDEFINITE(INDEF,A,IA,JA,N)

      IMPLICIT NONE

      INCLUDE 'bmg00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'

!     Parameter List      
      INTEGER N
      INTEGER IA(*),JA(*)
      REAL*8  A(*)
      LOGICAL INDEF

!     Local Variables
      INTEGER I,J
      REAL*8  SUM

C ===============================================
C     BEGIN EXECUTION.
C ===============================================

      INDEF = .TRUE.
      !
      DO I=1,N
         !
         SUM = 0.D0
         !
         DO J=IA(I),IA(I+1)-1
            SUM = SUM + A(J)
         ENDDO
         !
         IF( DABS(SUM).GT.1.D-10 ) THEN
             INDEF = .FALSE.
             RETURN
         ENDIF
         !
      ENDDO

      END
