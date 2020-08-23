      SUBROUTINE SOLVE_SYSTEM(ISC_A,ISR_A,LDA,M,N,NZA,OUTPUTCODE,PRECON,
     '  SOLVER,SPARSE_A,A,B,X,FIRST_A,UPDATE_A,X_INIT,NX,ERROR,*)

C#### Subroutine: SOLVE_SYSTEM
C###  Description:
C###    SOLVE_SYSTEM solves a linear system of equations Ax=b. The
C###    solver used is determined by SOLVER and the sparsity used
C###    is determined by SPARSENESS.
C###    FIRST_A and UPDATE_A should be .TRUE. for the first call.
C###    On return FIRST_A will be set to .FALSE..
C###    On subsequent calls UPDATE_A should be .TRUE./.FALSE. if
C###    refactorization is necessary/unnessary, and FIRST_A should
C###    be set to .TRUE. again if the sparsity pattern changes.

      IMPLICIT NONE

      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'ptr00.cmn'
      INCLUDE 'solv00.cmn'
      INCLUDE 'time02.cmn'
      INCLUDE 'solver.inc'

!     Parameter List
      INTEGER LDA,M,N,NZA,SPARSE_A,ISC_A(*),ISR_A(*),SOLVER,PRECON,
     '  OUTPUTCODE,NX
      REAL*8 A(*),B(N),X(N)
      CHARACTER ERROR*(*)
      LOGICAL FIRST_A,UPDATE_A,X_INIT
!     Local Variables
      INTEGER ITER
      REAL TIME1,TIME2,TIME3,FACT_TIME,SOLN_TIME,TOTL_TIME
      REAL*8 RESID
      LOGICAL DENSE,DIRECT,FACTOR,SAVE_SYSTEM,AMG1R6
C     Only needed when we have OpenMP
C$    INTEGER NCPU
C$    INTEGER OMP_GET_MAX_THREADS


      CALL ENTERS('SOLVE_SYSTEM',*9999)

       write(*,*) 'calling solve system'

      CALL ASSERT(nx.LE.9,'>>Can only have nx.le.9  at present',
     '  ERROR,*9999)

      DENSE=SPARSE_A.EQ.0
      DIRECT=(SOLVER.EQ.SOLV_LU .OR. SOLVER.EQ.SOLV_SVD .OR.
     '  SOLVER.EQ.SOLV_LSQ .OR. SOLVER.EQ.SOLV_CHOLESKY .OR.
     '  SOLVER.EQ.SOLV_LDL .OR. SOLVER.EQ.SOLV_SUPERLU  .OR.
     '  SOLVER.EQ.SOLV_UMFPACK2 .OR. SOLVER.EQ.SOLV_UMFPACK4)
      FACTOR=(UPDATE_A.OR.FIRST_A)
      SAVE_SYSTEM=DIRECT.AND.DENSE.AND.OUTPUTCODE.GE.1

      AMG1R6 = (SOLVER.EQ.SOLV_AMG)

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

C     Free/Allocate memory as necissary
      IF(FIRST_A) THEN
        CALL FREE_SOLVER(PRECON,SOLVER,SPARSE_A,NX,FIRST_A,ERROR,*9999)
        CALL ALLOC_SOLVER(LDA,M,N,NZA,SAVE_SYSTEM,PRECON,SOLVER,
     '    SPARSE_A,NX,ERROR,*9999)
      ELSE IF(UPDATE_A) THEN
        CALL FREE_SOLVER(PRECON,SOLVER,SPARSE_A,NX,FIRST_A,ERROR,*9999)
      ENDIF

C     Direct solvers: LU, SVD, ...
      IF(DIRECT) THEN

C       Dense solvers: use LAPACK
        IF(DENSE) THEN

          CALL CPU_TIMER(CPU_USER,TIME1)

C         Factorise the system
          IF(FACTOR) THEN
            IF(OUTPUTCODE.GE.1) THEN
              CALL PDCOPY(LDA*N,A,1,%VAL(RSOLV3_PTR(NX)),1)
            ENDIF

            CALL DIR_FACTOR(A,LDA,M,N,%VAL(ISOLV1_PTR(NX)),SOLVER,
     '        %VAL(RSOLV1_PTR(NX)),%VAL(RSOLV2_PTR(NX)),RSOLV1_LEN(NX),
     '        SOLV_ANORM(NX),ERROR,*9999)
            FIRST_A=.FALSE.
          ENDIF
          CALL CPU_TIMER(CPU_USER,TIME2)

C         Solve the problem
          CALL DIR_SOLVE(A,LDA,LDA,M,N,X,B,%VAL(ISOLV1_PTR(NX)),SOLVER,
     '      %VAL(RSOLV1_PTR(NX)),%VAL(RSOLV2_PTR(NX)),RSOLV2_LEN(NX),
     '      SOLV_ANORM(NX),OUTPUTCODE,ERROR,*9999)
          CALL CPU_TIMER(CPU_USER,TIME3)

C         Calculate and print the residual, and array/RHS norms etc
          IF(OUTPUTCODE.GE.1) THEN
            CALL CHECK_SOLN(%VAL(RSOLV3_PTR(NX)),LDA,N,B,X,SPARSE_A,
     '        ISC_A,ISR_A,ERROR,*9999)
          ENDIF

C       Sparse arrays: use Umfpack/SuperLU
        ELSE IF(SOLVER.EQ.SOLV_UMFPACK4) THEN

          CALL CPU_TIMER(CPU_USER,TIME1)

C         Factorise the system
          IF(FIRST_A) THEN
            CALL UMFPACK4_FACTOR(A,N,NZA,ISC_A,ISR_A,ISOLV1_PTR(NX),
     '        RSOLV1_PTR(NX),SOLV_UMF_PARAM(1,NX),
     '        SOLV_UMF_CONTROL(1,NX),SOLV_UMF_INFO(1,NX),
     '        SOLV_ANORM(NX),ERROR,*9999)
            FIRST_A=.FALSE.
          ELSE IF(UPDATE_A) THEN
            CALL UMFPACK4_REFACTOR(A,N,NZA,ISC_A,ISR_A,ISOLV1_PTR(NX),
     '        RSOLV1_PTR(NX),SOLV_UMF_CONTROL(1,NX),SOLV_UMF_INFO(1,NX),
     '        SOLV_ANORM(NX),ERROR,*9999)
          ENDIF
          CALL CPU_TIMER(CPU_USER,TIME2)

C         Solve the problem
          CALL UMFPACK4_SOLVE(A,N,NZA,ISC_A,ISR_A,B,X,RSOLV1_PTR(NX),
     '      SOLV_UMF_CONTROL(1,NX),SOLV_UMF_INFO(1,NX),SOLV_ANORM(NX),
     '      OUTPUTCODE,ERROR,*9999)
          CALL CPU_TIMER(CPU_USER,TIME3)

C         Calculate and print the residual, and array/RHS norms etc
          IF(OUTPUTCODE.GE.1) THEN
           CALL CHECK_SOLN(A,LDA,N,B,X,SPARSE_A,ISC_A,ISR_A,ERROR,*9999)
          ENDIF

C       The sparse, direct, SuperLU solver
        ELSE IF(SOLVER.EQ.SOLV_SUPERLU) THEN

          CALL CPU_TIMER(CPU_USER,TIME1)

C         Ensure SuperLU knows how many CPU's we are using. The following
C         commented code gets compiled in by OpenMP compilers.
C$        NCPU=OMP_GET_MAX_THREADS()
C$        CALL SUPERLU_SETPARAM('NUMBER OF PROCS',NCPU,
C$   '      SOLV_SLU_PARAM(1,NX),ERROR,*9999)

C         Factorise the system
          IF(FIRST_A) THEN
            CALL SUPERLU_FACTOR(A,N,NZA,ISC_A,ISR_A,RSOLV1_PTR(NX),
     '        SOLV_SLU_PARAM(1,NX),SOLV_ANORM(NX),ERROR,*9999)
            FIRST_A=.FALSE.
          ELSE IF(UPDATE_A) THEN
            CALL SUPERLU_REFACTOR(A,N,NZA,ISC_A,ISR_A,RSOLV1_PTR(NX),
     '        SOLV_SLU_PARAM(1,NX),SOLV_ANORM(NX),ERROR,*9999)
          ENDIF
          CALL CPU_TIMER(CPU_USER,TIME2)

C         Solve the problem
          CALL SUPERLU_SOLVE(B,X,N,NZA,RSOLV1_PTR(NX),
     '      SOLV_SLU_PARAM(1,NX),SOLV_ANORM(NX),OUTPUTCODE,ERROR,*9999)
          CALL CPU_TIMER(CPU_USER,TIME3)
 
C         Calculate and print the residual, and array/RHS norms etc
          IF(OUTPUTCODE.GE.1) THEN
           CALL CHECK_SOLN(A,LDA,N,B,X,SPARSE_A,ISC_A,ISR_A,ERROR,*9999)
          ENDIF 
        ELSE
          ERROR='>>Unknown sparse direct solver'
          GOTO 9999
        ENDIF

C     Iterative solvers: AMG, Jacobi, Conjugate Gradient, ...
      ELSE  

        IF(AMG1R6) THEN  ! John Ruge's Serial AMG

          CALL CPU_TIMER(CPU_USER,TIME1)
  
C         Create prolongation, restriction, and coarse-grid operators
          IF(FIRST_A) THEN 
            CALL AMG_SET_PARMS(AMG_CYCLE(NX),AMG_NDIMS(1,NX),
     '        AMG_IPARMS(1,NX),AMG_RPARMS(1,NX),SOLV_MAXIT(NX), 
     '        N,NX,NZA,SOLV_TOL(NX),ERROR,*9999)
            CALL AMG_COPY_MATRIX(N,AMG_NDIMS(1,NX),AMG_NDIMS(4,NX),
     '        NZA,ISC_A,ISR_A,A,%val(AMG_AA_PTR(NX)),
     '        %val(AMG_IA_PTR(NX)),%val(AMG_JA_PTR(NX)),
     '        %val(AMG_BDY_PTR(NX)),SPARSE_A,ERROR,*9999)
            CALL AMG_WRAPPER(%val(AMG_AA_PTR(NX)), 
     '        %val(AMG_IA_PTR(NX)),%val(AMG_JA_PTR(NX)), 
     '        %val(AMG_U_PTR(NX)),%val(AMG_F_PTR(NX)),
     '        %val(AMG_IG_PTR(NX)),A,ISR_A,ISC_A,X,B,
     '        %val(AMG_BDY_PTR(NX)),AMG_NDIMS(1,NX),
     '        AMG_IPARMS(1,NX),AMG_RPARMS(1,NX),
     '        ERROR,*9999)    
          ELSEIF(UPDATE_A) THEN     
            CALL AMG_COPY_MATRIX(N,AMG_NDIMS(1,NX),AMG_NDIMS(4,NX),
     '        NZA,ISC_A,ISR_A,A,%val(AMG_AA_PTR(NX)),
     '        %val(AMG_IA_PTR(NX)),%val(AMG_JA_PTR(NX)),  
     '        %val(AMG_BDY_PTR(NX)),SPARSE_A,ERROR,*9999)  
            CALL AMG_WRAPPER(%val(AMG_AA_PTR(NX)), 
     '        %val(AMG_IA_PTR(NX)),%val(AMG_JA_PTR(NX)), 
     '        %val(AMG_U_PTR(NX)),%val(AMG_F_PTR(NX)),  
     '        %val(AMG_IG_PTR(NX)),A,ISR_A,ISC_A,X,B,
     '        %val(AMG_BDY_PTR(NX)),AMG_NDIMS(1,NX), 
     '        AMG_IPARMS(1,NX),AMG_RPARMS(1,NX),
     '        ERROR,*9999)     
          ENDIF 
          CALL CPU_TIMER(CPU_USER,TIME2)   

C         Solve the problem  
          AMG_IPARMS(2,NX) = 2    
          CALL AMG_WRAPPER(%val(AMG_AA_PTR(NX)),
     '      %val(AMG_IA_PTR(NX)),%val(AMG_JA_PTR(NX)),
     '      %val(AMG_U_PTR(NX)),%val(AMG_F_PTR(NX)),
     '      %val(AMG_IG_PTR(NX)),A,ISR_A,ISC_A,X,B,
     '      %val(AMG_BDY_PTR(NX)),AMG_NDIMS(1,NX), 
     '      AMG_IPARMS(1,NX),AMG_RPARMS(1,NX),
     '      ERROR,*9999)  
          CALL CPU_TIMER(CPU_USER,TIME3)
                     
        ELSE             ! Stuart Norris iterative solvers 

          CALL CPU_TIMER(CPU_USER,TIME1)  
 
C         Factorise the system 
          IF(FACTOR) THEN
            CALL ITER_FACTOR(A,LDA,N,%VAL(RSOLV1_PTR(NX)),
     '        SOLV_ANORM(NX),SPARSE_A,ISC_A,ISR_A,SOLVER,PRECON,
     '        OUTPUTCODE,ERROR,*9999)
            FIRST_A=.FALSE.
          ENDIF
          CALL CPU_TIMER(CPU_USER,TIME2) 

C         Solve the problem
          IF(X_INIT) CALL PDLOAD(N,0.0D0,X,1)
          ITER=SOLV_MAXIT(NX)
          RESID=SOLV_TOL(NX)
          CALL ITER_SOLVE(A,LDA,N,X,B,%VAL(RSOLV1_PTR(NX)),SPARSE_A,
     '      ISC_A,ISR_A,SOLV_ANORM(NX),RESID,SOLV_OMEGA(NX),ITER,
     '      SOLV_NRES(NX),SOLV_NPRECON(NX),SOLVER,PRECON,OUTPUTCODE,
     '      ERROR,*9999)
          CALL CPU_TIMER(CPU_USER,TIME3)

C         Number of Iterations/Residual
          IF(OUTPUTCODE.GE.1) THEN
            WRITE(OP_STRING,7000) ' Number of iterations: ',ITER
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,7100) ' Residual of solution: ',RESID
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF

C         Calculate and print the residual, and array/RHS norms etc
          IF(OUTPUTCODE.GE.1) THEN
           CALL CHECK_SOLN(A,LDA,N,B,X,SPARSE_A,ISC_A,ISR_A,ERROR,*9999)
          ENDIF

        ENDIF

      ENDIF

C     Timing information
      IF(OUTPUTCODE.GE.1) THEN
        FACT_TIME=TIME2-TIME1
        SOLN_TIME=TIME3-TIME2
        TOTL_TIME=TIME3-TIME1

        IF(FACTOR) THEN
          WRITE(OP_STRING,7200) ' CPU time to factorise system:     ',
     '      FACT_TIME,'s'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

        WRITE(OP_STRING,7200) ' CPU time to solve system:         ',
     '    SOLN_TIME,'s'
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

        IF(FACTOR) THEN
          WRITE(OP_STRING,7200) ' Sum of factor and solution times: ',
     '      TOTL_TIME,'s'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('SOLVE_SYSTEM')
      RETURN

 7000 FORMAT(A,I10)
 7100 FORMAT(A,G12.6)
 7200 FORMAT(A,1P,D11.4,1X,A)

 9999 CALL ERRORS('SOLVE_SYSTEM',ERROR)
      CALL EXITS('SOLVE_SYSTEM')
      RETURN 1
      END
