      SUBROUTINE ALLOC_SOLVER(LDA,M,N,NZA,SAVE_SYSTEM,PRECON,SOLVER,
     '  SPARSE_A,NX,ERROR,*)

C#### Subroutine: ALLOC_SOLVER
C###  Description:
C###    ALLOC_SOLVER allocates memory for the linear solvers in SOLVE_SYSTEM.
C###    On completion the SOLV_ALLOC variable is set to true to flag that
C###    memory is allocated.

      IMPLICIT NONE

      INCLUDE 'cbfe01.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'ptr00.cmn'
      INCLUDE 'solv00.cmn'
      INCLUDE 'solver.inc'

!     Parameter List
      INTEGER LDA,M,N,NZA,SPARSE_A,SOLVER,PRECON,NX
      CHARACTER ERROR*(*)
      LOGICAL SAVE_SYSTEM
!     Local Variables
      LOGICAL DENSE,DIRECT,AMG1R6


      CALL ENTERS('ALLOC_SOLVER',*9999)

      CALL ASSERT(nx.LE.9,'>>Can only have nx.le.9  at present',
     '  ERROR,*9999)

      DENSE=SPARSE_A.EQ.0
      DIRECT=(SOLVER.EQ.SOLV_LU  .OR. SOLVER.EQ.SOLV_UMFPACK2 .OR.
     '        SOLVER.EQ.SOLV_LSQ .OR. SOLVER.EQ.SOLV_CHOLESKY .OR.
     '        SOLVER.EQ.SOLV_LDL .OR. SOLVER.EQ.SOLV_SUPERLU  .OR.
     '        SOLVER.EQ.SOLV_SVD .OR. SOLVER.EQ.SOLV_UMFPACK4)

      AMG1R6=(SOLVER.EQ.SOLV_AMG)

C     Direct solvers: LU, SVD, ...
      IF(DIRECT) THEN

C       Dense solvers: use LAPACK
        IF(DENSE) THEN
          CALL DIR_ALLOC(LDA,M,N,SOLVER,RSOLV1_PTR(NX),RSOLV2_PTR(NX),
     '      ISOLV1_PTR(NX),RSOLV1_LEN(NX),RSOLV2_LEN(NX),ERROR,*9999)

C         Array to store original system, so we can calculate the residual
          IF(SAVE_SYSTEM) THEN
            CALL ALLOCATE_MEMORY(LDA*N,1,DPTYPE,RSOLV3_PTR(NX),MEM_INIT,
     '        ERROR,*9999)
          ENDIF

C       Sparse arrays: use Umfpack/SuperLU
        ELSE IF(SOLVER.EQ.SOLV_UMFPACK4) THEN

C       The sparse, direct, SuperLU solver
        ELSE IF(SOLVER.EQ.SOLV_SUPERLU) THEN

        ELSE
          ERROR='>>Unknown sparse direct solver'
          GOTO 9999
        ENDIF

C     Iterative solvers: AMG1r6, Jacobi, Conjugate Gradient, ...
      ELSE 

        IF(AMG1R6) THEN ! Algebraic Multigrid (Ruge's Serial Version) 
          CALL AMG_ALLOC(LDA,N,NZA,AMG_U_PTR(NX),AMG_F_PTR(NX),
     '      AMG_IG_PTR(NX),AMG_AA_PTR(NX),AMG_IA_PTR(NX),AMG_JA_PTR(NX), 
     '      AMG_BDY_PTR(NX),ERROR,*9999)
        ELSE            ! Solvers designed by Stuart Norris
          CALL ITER_ALLOC(LDA,N,NZA,SOLVER,PRECON,SPARSE_A,
     '      RSOLV1_PTR(nx),ERROR,*9999)
        ENDIF  

      ENDIF

C     Flag memory as allocated
      SOLV_ALLOC(NX)=.TRUE.

      CALL EXITS('ALLOC_SOLVER')
      RETURN

 9999 CALL ERRORS('ALLOC_SOLVER',ERROR)
      CALL EXITS('ALLOC_SOLVER')
      RETURN 1
      END


