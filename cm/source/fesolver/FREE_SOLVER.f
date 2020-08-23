      SUBROUTINE FREE_SOLVER(PRECON,SOLVER,SPARSE_A,NX,FIRST_A,ERROR,*)

C#### Subroutine: FREE_SOLVER
C###  Description:
C###    FREE_SOLVER frees memory used by the linear solvers in SOLVE_SYSTEM.
C###    The memory is flagged as allocated by the logical variable SOLV_ALLOC
C###    held in common and described in the ptr00.cmn file.

      IMPLICIT NONE

      INCLUDE 'ptr00.cmn'
      INCLUDE 'solv00.cmn'
      INCLUDE 'solver.inc'

!     Parameter List
      INTEGER SOLVER,PRECON,SPARSE_A,NX
      LOGICAL FIRST_A
      CHARACTER ERROR*(*)
!     Local Variables
      LOGICAL DENSE,DIRECT,AMG1R6,FREE_I,FREE_R
!     Block data we wish to load
      EXTERNAL SOLVER_MEMORY_FLAGS


      CALL ENTERS('FREE_SOLVER',*9999)

      CALL ASSERT(nx.LE.9,'>>Can only have nx.le.9  at present',
     '  ERROR,*9999)

      DENSE=SPARSE_A.EQ.0
      DIRECT=(SOLVER.EQ.SOLV_LU .OR. SOLVER.EQ.SOLV_SVD .OR.
     '  SOLVER.EQ.SOLV_LSQ .OR. SOLVER.EQ.SOLV_CHOLESKY .OR.
     '  SOLVER.EQ.SOLV_LDL .OR. SOLVER.EQ.SOLV_SUPERLU  .OR.
     '  SOLVER.EQ.SOLV_UMFPACK2 .OR. SOLVER.EQ.SOLV_UMFPACK4)

      AMG1R6=(SOLVER.EQ.SOLV_AMG)

      IF(SOLV_ALLOC(NX)) THEN

C       Direct solvers: LU, SVD, ...
        IF(DIRECT) THEN

C         Dense solvers: use LAPACK
          IF(DENSE) THEN
            IF(FIRST_A) THEN
              CALL DIR_FREE(SOLVER,RSOLV1_PTR(NX),RSOLV2_PTR(NX),
     '          ISOLV1_PTR(NX),ERROR,*9999)

              IF(RSOLV3_PTR(NX).NE.0) THEN
                CALL FREE_MEMORY(RSOLV3_PTR(NX),ERROR,*9999)
              ENDIF
            ENDIF

C         Sparse arrays: use Umfpack/SuperLU
          ELSE IF(SOLVER.EQ.SOLV_UMFPACK4) THEN
            FREE_R=SOLV_ALLOC(NX)
            FREE_I=FIRST_A
            CALL UMFPACK4_FREE(ISOLV1_PTR(NX),FREE_I,RSOLV1_PTR(NX),
     '        FREE_R,ERROR,*9999)

C         The sparse, direct, SuperLU solver
          ELSE IF(SOLVER.EQ.SOLV_SUPERLU) THEN
            IF(FIRST_A) THEN
              CALL SUPERLU_FREE(RSOLV1_PTR(NX),ERROR,*9999)
            ENDIF

          ELSE
            ERROR='>>Unknown sparse direct solver'
            GOTO 9999
          ENDIF

C       Iterative solvers: Jacobi, Conjugate Gradient, ...
        ELSE

          IF(FIRST_A) THEN
            IF(AMG1R6) THEN
              CALL AMG_FREE(AMG_U_PTR(nx),AMG_F_PTR(nx),AMG_IG_PTR(nx),
     '          AMG_AA_PTR(nx),AMG_IA_PTR(nx),AMG_JA_PTR(nx),
     '          AMG_BDY_PTR(nx),ERROR,*9999)
            ELSE
              CALL ITER_FREE(SOLVER,PRECON,RSOLV1_PTR(nx),ERROR,*9999)
            ENDIF
          ENDIF

        ENDIF

      ENDIF

C     Flag memory as free'ed
      IF(FIRST_A) SOLV_ALLOC(NX)=.FALSE.

      CALL EXITS('FREE_SOLVER')
      RETURN

 9999 CALL ERRORS('FREE_SOLVER',ERROR)
      CALL EXITS('FREE_SOLVER')
      RETURN 1
      END


