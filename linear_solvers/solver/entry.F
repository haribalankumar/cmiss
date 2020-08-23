
      SUBROUTINE PRECON_ALLOC(LDA,N,NZA,PRECON,SPARSE_A,D_PTR,ERROR,*)

C#### Subroutine: PRECON_ALLOC
C###  Description:
C###    PRECON_ALLOC allocates memory for the factorised matrix
C###    used by the preconditioner specified by PRECON.
C###  Written by Stuart Norris 28/05/02

      IMPLICIT NONE
      INCLUDE 'mach00.inc'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'solver.inc'
!     Parameter List
      INTEGER LDA,N,NZA,PRECON,SPARSE_A
      CHARACTER ERROR*(*)
      POINTER D_PTR
!     Local Variables
      INTEGER NDIM


C     No preconditioning
      IF(PRECON.EQ.PRE_NONE) THEN
        D_PTR=0

C     Stationary methods
      ELSE IF(PRECON.EQ.PRE_POINT .OR. PRECON.EQ.PRE_JACOBI .OR. 
     '    PRECON.EQ.PRE_SOR .OR. PRECON.EQ.PRE_SCALE) THEN
        CALL ALLOCATE_MEMORY(N,1,DPTYPE,D_PTR,NOINIT,ERROR,*9999)

C     Level 0 Incomplete factorisations
      ELSE IF(PRECON.EQ.PRE_ILU0 .OR. PRECON.EQ.PRE_IC0 .OR.
     '    PRECON.EQ.PRE_ILDL0) THEN
        CALL ALLOCATE_MEMORY(N,1,DPTYPE,D_PTR,NOINIT,ERROR,*9999)

C     Level 1 Incomplete factorisations
      ELSE IF(PRECON.EQ.PRE_ILU1 .OR. PRECON.EQ.PRE_IC1 .OR. 
     '    PRECON.EQ.PRE_ILDL1) THEN
        IF(SPARSE_A.EQ.0) THEN
          NDIM=LDA*N
        ELSE IF(SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4) THEN
          NDIM=NZA
        ELSE
          ERROR='>>Sparsity type not implemented'
          GOTO 9999
        ENDIF
        CALL ALLOCATE_MEMORY(NDIM,1,DPTYPE,D_PTR,NOINIT,ERROR,*9999)

C     Bzzzzzzt, Error
      ELSE
        ERROR='>>No Such Preconditioner'
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('PRECON_ALLOC',ERROR)
      RETURN 1
      END


      SUBROUTINE PRECON_FREE(PRECON,D_PTR,ERROR,*)

C#### Subroutine: PRECON_FREE
C###  Description:
C###    PRECON_FREE frees memory allocated in PRECON_ALLOC.
C###  Written by Stuart Norris 04/06/02

      IMPLICIT NONE
      INCLUDE 'solver.inc'
!     Parameter List
      INTEGER PRECON
      CHARACTER ERROR*(*)
      POINTER D_PTR
!     Local Variables


C     No preconditioning
      IF(PRECON.EQ.PRE_NONE) THEN

C     Stationary methods
      ELSE IF(PRECON.EQ.PRE_POINT .OR. PRECON.EQ.PRE_JACOBI .OR. 
     '    PRECON.EQ.PRE_SOR .OR. PRECON.EQ.PRE_SCALE) THEN
        IF(D_PTR.NE.0) CALL FREE_MEMORY(D_PTR,ERROR,*9999)

C     Level 0 Incomplete factorisations
      ELSE IF(PRECON.EQ.PRE_ILU0 .OR. PRECON.EQ.PRE_ILDL0 .OR.
     '    PRECON.EQ.PRE_IC0) THEN
        IF(D_PTR.NE.0) CALL FREE_MEMORY(D_PTR,ERROR,*9999)

C     Level 1 Incomplete factorisations
      ELSE IF(PRECON.EQ.PRE_ILU1 .OR. PRECON.EQ.PRE_ILDL1 .OR. 
     '    PRECON.EQ.PRE_IC1) THEN
        IF(D_PTR.NE.0) CALL FREE_MEMORY(D_PTR,ERROR,*9999)

C     Bzzzzzzt, Error
      ELSE
        ERROR='>>No Such Preconditioner'
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('PRECON_FREE',ERROR)
      RETURN 1
      END


      SUBROUTINE PRECON_FACTOR(A,LDA,N,D,PRECON,SPARSE_A,ISC_A,ISR_A,
     '  OUTPUTCODE,ERROR,*)

C#### Subroutine: PRECON_FACTOR
C###  Description:
C###    PRECON_FACTOR factorises a system that is to be used in
C###    a preconditioned conjugate gradient solver.
C###  Written by Stuart Norris 26/02/02

      IMPLICIT NONE
      INCLUDE 'mach00.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'solver.inc'
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*)
      INTEGER LDA,N,PRECON,SPARSE_A,OUTPUTCODE
      REAL*8 A(*),D(*)
      CHARACTER ERROR*(*)
!     Local Variables


C     No preconditioning
      IF(PRECON.EQ.PRE_NONE) THEN
        WRITE(OP_STRING,'('' No Preconditioning'')')
        IF(OUTPUTCODE.GE.2) CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

C     Stationary methods
      ELSE IF(PRECON.EQ.PRE_POINT .OR. PRECON.EQ.PRE_JACOBI .OR. 
     '    PRECON.EQ.PRE_SOR .OR. PRECON.EQ.PRE_SCALE) THEN
        IF(PRECON.EQ.PRE_POINT) THEN
          WRITE(OP_STRING,'('' Point Jacobi Preconditioning'')')
          IF(OUTPUTCODE.GE.2) CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CALL POINT_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,ERROR,*9999)
        ELSE IF(PRECON.EQ.PRE_SCALE) THEN
          WRITE(OP_STRING,'('' Row-Scale Point Jacobi '//
     '      'Preconditioning'')')
          IF(OUTPUTCODE.GE.2) CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CALL ROWSCALE_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,ERROR,
     '      *9999)
        ELSE IF(PRECON.EQ.PRE_JACOBI) THEN
          WRITE(OP_STRING,'('' Jacobi Preconditioning'')')
          IF(OUTPUTCODE.GE.2) CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CALL JACOBI_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,ERROR,*9999)
        ELSE IF(PRECON.EQ.PRE_SOR) THEN
          WRITE(OP_STRING,'('' SSOR Preconditioning'')')
          IF(OUTPUTCODE.GE.2) CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CALL JACOBI_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,ERROR,*9999)
        ENDIF

C     Level 0 Incomplete factorisations
      ELSE IF(PRECON.EQ.PRE_IC0) THEN
        WRITE(OP_STRING,'('' Incomplete Cholesky(0) Preconditioning'')')
        IF(OUTPUTCODE.GE.2) CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        CALL INCOMP_CHOL0_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,ERROR,
     '    *9999)
      ELSE IF(PRECON.EQ.PRE_ILU0) THEN
        WRITE(OP_STRING,'('' Incomplete LU(0) Preconditioning'')')
        IF(OUTPUTCODE.GE.2) CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        CALL INCOMP_LU0_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,ERROR,
     '    *9999)
      ELSE IF(PRECON.EQ.PRE_ILDL0) THEN
        WRITE(OP_STRING,'('' Incomplete LDL(0) Preconditioning'')')
        IF(OUTPUTCODE.GE.2) CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        CALL INCOMP_LDL0_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,ERROR,
     '    *9999)

C     Level 1 Incomplete factorisations
      ELSE IF(PRECON.EQ.PRE_IC1) THEN
        WRITE(OP_STRING,'('' Incomplete Cholesky(1) Preconditioning'')')
        IF(OUTPUTCODE.GE.2) CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        CALL INCOMP_CHOL1_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,ERROR,
     '    *9999)
      ELSE IF(PRECON.EQ.PRE_ILU1) THEN
        WRITE(OP_STRING,'('' Incomplete LU(1) Preconditioning'')')
        IF(OUTPUTCODE.GE.2) CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        CALL INCOMP_LU1_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,ERROR,
     '    *9999)
      ELSE IF(PRECON.EQ.PRE_ILDL1) THEN
        WRITE(OP_STRING,'('' Incomplete LDL(1) Preconditioning'')')
        IF(OUTPUTCODE.GE.2) CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        CALL INCOMP_LDL1_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,ERROR,
     '    *9999)

C     Bzzzzzzt, Error
      ELSE
        ERROR='>>No Such Preconditioner'
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('PRECON_FACTOR',ERROR)
      RETURN 1
      END


      SUBROUTINE PRECON_SOLVE(A,LDA,N,X,B,D,E,PRECON,SPARSE_A,
     '  ISC_A,ISR_A,OMEGA,MAXIT,ERROR,*)

C#### Subroutine: PRECON_SOLVE
C###  Description:
C###    PRECON_SOLVE applies ITER iterations of the PRECON preconditioner
C###    to the system A x = b.
C###  Written by Stuart Norris 26/02/02

      IMPLICIT NONE
      INCLUDE 'solver.inc'
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*),LDA,N,PRECON,SPARSE_A,MAXIT
      REAL*8 A(*),X(*),B(*),D(*),E(*),OMEGA
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I
      EXTERNAL PDCOPY


C SEN Do we need to do this?
C$OMP PARALLEL DO PRIVATE(I)
      DO I=1,N
        X(I)=0.0D0
      ENDDO
C$OMP END PARALLEL DO

C     No preconditioning
      IF(PRECON.EQ.PRE_NONE) THEN
        CALL PDCOPY(N,B,1,X,1)

C     Stationary methods
      ELSE IF(PRECON.EQ.PRE_POINT .OR. PRECON.EQ.PRE_SCALE) THEN
C$OMP   PARALLEL DO PRIVATE(I)
        DO I=1,N
          X(I)=B(I)*D(I)
        ENDDO
C$OMP   END PARALLEL DO
      ELSE IF(PRECON.EQ.PRE_JACOBI) THEN
        CALL PRECON_JACOBI(A,LDA,N,X,B,D,E,SPARSE_A,ISC_A,ISR_A,OMEGA,
     '    MAXIT,ERROR,*9999)
      ELSE IF(PRECON.EQ.PRE_SOR) THEN
        CALL PRECON_SSOR(A,LDA,N,X,B,D,SPARSE_A,ISC_A,ISR_A,OMEGA,
     '    MAXIT,ERROR,*9999)

C     Level 0 Incomplete factorisations
      ELSE IF(PRECON.EQ.PRE_IC0) THEN
        CALL PRECON_INCOMP_CHOL0(A,LDA,N,X,B,D,E,SPARSE_A,ISC_A,ISR_A,
     '    MAXIT,ERROR,*9999)
      ELSE IF(PRECON.EQ.PRE_ILU0) THEN
        CALL PRECON_INCOMP_LU0(A,LDA,N,X,B,D,E,SPARSE_A,ISC_A,ISR_A,
     '    MAXIT,ERROR,*9999)
      ELSE IF(PRECON.EQ.PRE_ILDL0) THEN
        CALL PRECON_INCOMP_LDL0(A,LDA,N,X,B,D,E,SPARSE_A,ISC_A,ISR_A,
     '    MAXIT,ERROR,*9999)

C     Level 1 Incomplete factorisations
      ELSE IF(PRECON.EQ.PRE_IC1) THEN
        CALL PRECON_INCOMP_CHOL1(A,LDA,N,X,B,D,E,SPARSE_A,ISC_A,ISR_A,
     '    MAXIT,ERROR,*9999)
      ELSE IF(PRECON.EQ.PRE_ILU1) THEN
        CALL PRECON_INCOMP_LU1(A,LDA,N,X,B,D,E,SPARSE_A,ISC_A,ISR_A,
     '    MAXIT,ERROR,*9999)
      ELSE IF(PRECON.EQ.PRE_ILDL1) THEN
        CALL PRECON_INCOMP_LDL1(A,LDA,N,X,B,D,E,SPARSE_A,ISC_A,ISR_A,
     '    MAXIT,ERROR,*9999)

C     Bzzzzzzt, Error
      ELSE
        ERROR='>>No Such Preconditioner'
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('PRECON_SOLVE',ERROR)
      RETURN 1
      END


      SUBROUTINE ITER_ALLOC(LDA,N,NZA,SOLVER,PRECON,SPARSE_A,D_PTR,
     '  ERROR,*)

C#### Subroutine: ITER_ALLOC
C###  Description:
C###    ITER_ALLOC allocates memory for the factorised matrix
C###    used by the solver specified by PRECON.
C###  Written by Stuart Norris 28/05/02

      IMPLICIT NONE
      INCLUDE 'mach00.inc'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'solver.inc'
!     Parameter List
      INTEGER LDA,N,NZA,SOLVER,PRECON,SPARSE_A
      CHARACTER ERROR*(*)
      POINTER D_PTR
!     Local Variables
      INTEGER NDIM


C     Stationary methods
      IF(SOLVER.EQ.SOLV_JACOBI .OR. SOLVER.EQ.SOLV_SOR) THEN
        NDIM=N
        CALL ALLOCATE_MEMORY(NDIM,1,DPTYPE,D_PTR,NOINIT,ERROR,*9999)

C     Level 0 Incomplete factorisations
      ELSE IF(SOLVER.EQ.SOLV_IC0 .OR. SOLVER.EQ.SOLV_ILU0 .OR. 
     '    SOLVER.EQ.SOLV_ILDL0) THEN
        NDIM=N
        CALL ALLOCATE_MEMORY(NDIM,1,DPTYPE,D_PTR,NOINIT,ERROR,*9999)

C     Level 1 Incomplete factorisations
      ELSE IF(SOLVER.EQ.SOLV_IC1 .OR. SOLVER.EQ.SOLV_ILU1 .OR. 
     '    SOLVER.EQ.SOLV_ILDL1) THEN
        IF(SPARSE_A.EQ.0) THEN
          NDIM=LDA*N
        ELSE IF(SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4) THEN
          NDIM=NZA
        ELSE
          ERROR='>>Sparsity type not implemented'
          GOTO 9999
        ENDIF
        CALL ALLOCATE_MEMORY(NDIM,1,DPTYPE,D_PTR,NOINIT,ERROR,*9999)

C     Krylov space methods -- preconditioners matricies
      ELSE IF(SOLVER.EQ.SOLV_CG .OR. SOLVER.EQ.SOLV_BICGSTAB .OR. 
     '    SOLVER.EQ.SOLV_GMRES) THEN
        CALL PRECON_ALLOC(LDA,N,NZA,PRECON,SPARSE_A,D_PTR,ERROR,*9999)

C     Bzzzzzzt, Error
      ELSE
        ERROR='>>No Such Solver'
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('ITER_ALLOC',ERROR)
      RETURN 1
      END


      SUBROUTINE ITER_FREE(SOLVER,PRECON,D_PTR,ERROR,*)

C#### Subroutine: ITER_FREE
C###  Description:
C###    ITER_FREE frees memory allocated with ITER_ALLOC.
C###  Written by Stuart Norris 04/06/02

      IMPLICIT NONE
      INCLUDE 'solver.inc'
!     Parameter List
      INTEGER SOLVER,PRECON
      CHARACTER ERROR*(*)
      POINTER D_PTR
!     Local Variables


C     Stationary methods
      IF(SOLVER.EQ.SOLV_JACOBI .OR. SOLVER.EQ.SOLV_SOR) THEN
        IF(D_PTR.NE.0) CALL FREE_MEMORY(D_PTR,ERROR,*9999)

C     Level 0 Incomplete factorisations
      ELSE IF(SOLVER.EQ.SOLV_IC0 .OR. SOLVER.EQ.SOLV_ILU0 .OR. 
     '    SOLVER.EQ.SOLV_ILDL0) THEN
        IF(D_PTR.NE.0) CALL FREE_MEMORY(D_PTR,ERROR,*9999)

C     Level 1 Incomplete factorisations
      ELSE IF(SOLVER.EQ.SOLV_IC1 .OR. SOLVER.EQ.SOLV_ILU1 .OR. 
     '    SOLVER.EQ.SOLV_ILDL1) THEN
        IF(D_PTR.NE.0) CALL FREE_MEMORY(D_PTR,ERROR,*9999)

C     Krylov space methods -- preconditioners matricies
      ELSE IF(SOLVER.EQ.SOLV_CG .OR. SOLVER.EQ.SOLV_BICGSTAB .OR. 
     '    SOLVER.EQ.SOLV_GMRES) THEN
        CALL PRECON_FREE(PRECON,D_PTR,ERROR,*9999)

C     Bzzzzzzt, Error
      ELSE
        ERROR='>>No Such Solver'
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('ITER_FREE',ERROR)
      RETURN 1
      END


      SUBROUTINE ITER_FACTOR(A,LDA,N,D,ANORM,SPARSE_A,ISC_A,ISR_A,
     '  SOLVER,PRECON,OUTPUTCODE,ERROR,*)

C#### Subroutine: ITER_FACTOR
C###  Description:
C###    ITER_FACTOR factorises a system that is to be used in
C###    one of the iterative solvers
C###  Written by Stuart Norris 26/03/02

      IMPLICIT NONE
      INCLUDE 'mach00.inc'
      INCLUDE 'solver.inc'
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*),LDA,N,SOLVER,PRECON,SPARSE_A,OUTPUTCODE
      REAL*8 A(*),D(*),ANORM
      CHARACTER ERROR*(*)
!     Local Constants
      INTEGER NORMTYPE
      PARAMETER(NORMTYPE=2)
!     Local variables


C     Factorise the system:
C     Stationary methods
      IF(SOLVER.EQ.SOLV_JACOBI .OR. SOLVER.EQ.SOLV_SOR) THEN
        CALL JACOBI_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,ERROR,*9999)

C     Level 0 Incomplete factorisations
      ELSE IF(SOLVER.EQ.SOLV_IC0) THEN
        CALL INCOMP_CHOL0_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,ERROR,
     '    *9999)
      ELSE IF(SOLVER.EQ.SOLV_ILU0) THEN
        CALL INCOMP_LU0_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,ERROR,
     '    *9999)
      ELSE IF(SOLVER.EQ.SOLV_ILDL0) THEN
        CALL INCOMP_LDL0_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,ERROR,
     '    *9999)

C     Level 1 Incomplete factorisations
      ELSE IF(SOLVER.EQ.SOLV_IC1) THEN
        CALL INCOMP_CHOL1_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,ERROR,
     '    *9999)
      ELSE IF(SOLVER.EQ.SOLV_ILU1) THEN
        CALL INCOMP_LU1_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,ERROR,
     '    *9999)
      ELSE IF(SOLVER.EQ.SOLV_ILDL1) THEN
        CALL INCOMP_LDL1_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,ERROR,
     '    *9999)

C     Krylov space methods -- we might want to precondition
      ELSE IF(SOLVER.EQ.SOLV_CG .OR. SOLVER.EQ.SOLV_BICGSTAB .OR. 
     '    SOLVER.EQ.SOLV_GMRES) THEN
        CALL PRECON_FACTOR(A,LDA,N,D,PRECON,SPARSE_A,ISC_A,ISR_A,
     '    OUTPUTCODE,ERROR,*9999)

      ELSE
        ERROR='>>No such iterative solver'
        GOTO 9999
      ENDIF

C     Get the norm of A
      CALL CALC_MAT_NORM(ANORM,A,LDA,N,SPARSE_A,ISC_A,ISR_A,NORMTYPE,
     '  ERROR,*9999)

      RETURN

 9999 CALL ERRORS('ITER_FACTOR',ERROR)
      RETURN 1
      END


      SUBROUTINE ITER_SOLVE(A,LDA,N,X,B,D,SPARSE_A,ISC_A,ISR_A,ANORM,
     '  RESID,OMEGA,ITER,NRES,NPRECON,SOLVER,PRECON,OUTPUTCODE,ERROR,*)

C#### Subroutine: ITER_SOLVE
C###  Description:
C###    ITER_SOLV is a wrapper around the iterative solvers, 
C###    allowing a more intelligable naming of arrays.
C###  Written by Stuart Norris 28/02/02

      IMPLICIT NONE
      INCLUDE 'solver.inc'
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*),LDA,N,SPARSE_A
      INTEGER ITER,NRES,NPRECON,SOLVER,PRECON,OUTPUTCODE
      REAL*8 A(*),X(*),B(*),D(*),ANORM,OMEGA,RESID
      CHARACTER ERROR*(*)
!     Local constants
      LOGICAL FACTORISE
      PARAMETER(FACTORISE=.FALSE.)
!     Local Variables
      INTEGER NRESL


C     Fix boundary conditions
      CALL SOLV_FIX_BOUNDS(A,LDA,N,X,B,SPARSE_A,ISC_A,ISR_A,ERROR,*9999)

C     Stationary methods
      IF(SOLVER.EQ.SOLV_JACOBI) THEN
        CALL JACOBI(A,LDA,N,X,B,D,SPARSE_A,ISC_A,ISR_A,ANORM,RESID,
     '    OMEGA,ITER,FACTORISE,OUTPUTCODE,ERROR,*9999)
      ELSE IF(SOLVER.EQ.SOLV_SOR) THEN
        CALL SOR(A,LDA,N,X,B,D,SPARSE_A,ISC_A,ISR_A,ANORM,RESID,
     '    OMEGA,ITER,FACTORISE,OUTPUTCODE,ERROR,*9999)

C     Level 0 Incomplete factorisations
      ELSE IF(SOLVER.EQ.SOLV_IC0) THEN
        CALL INCOMP_CHOL0(A,LDA,N,X,B,D,SPARSE_A,ISC_A,ISR_A,ANORM,
     '    RESID,ITER,FACTORISE,OUTPUTCODE,ERROR,*9999)
      ELSE IF(SOLVER.EQ.SOLV_ILU0) THEN
        CALL INCOMP_LU0(A,LDA,N,X,B,D,SPARSE_A,ISC_A,ISR_A,ANORM,
     '    RESID,ITER,FACTORISE,OUTPUTCODE,ERROR,*9999)
      ELSE IF(SOLVER.EQ.SOLV_ILDL0) THEN
        CALL INCOMP_LDL0(A,LDA,N,X,B,D,SPARSE_A,ISC_A,ISR_A,ANORM,
     '    RESID,ITER,FACTORISE,OUTPUTCODE,ERROR,*9999)

C     Level 1 Incomplete factorisations
      ELSE IF(SOLVER.EQ.SOLV_IC1) THEN
        CALL INCOMP_CHOL1(A,LDA,N,X,B,D,SPARSE_A,ISC_A,ISR_A,ANORM,
     '    RESID,ITER,FACTORISE,OUTPUTCODE,ERROR,*9999)
      ELSE IF(SOLVER.EQ.SOLV_ILU1) THEN
        CALL INCOMP_LU1(A,LDA,N,X,B,D,SPARSE_A,ISC_A,ISR_A,ANORM,
     '    RESID,ITER,FACTORISE,OUTPUTCODE,ERROR,*9999)
      ELSE IF(SOLVER.EQ.SOLV_ILDL1) THEN
        CALL INCOMP_LDL1(A,LDA,N,X,B,D,SPARSE_A,ISC_A,ISR_A,ANORM,
     '    RESID,ITER,FACTORISE,OUTPUTCODE,ERROR,*9999)

C     Krylov space methods
      ELSE IF(SOLVER.EQ.SOLV_CG) THEN
        NRESL=NRES
        IF(NRES.LE.0) NRESL=ITER+1
        CALL CONJGRAD(A,LDA,N,X,B,D,SPARSE_A,ISC_A,ISR_A,ANORM,RESID,
     '    OMEGA,ITER,NRESL,NPRECON,PRECON,FACTORISE,OUTPUTCODE,
     '    ERROR,*9999)
      ELSE IF(SOLVER.EQ.SOLV_BICGSTAB) THEN
        NRESL=NRES
        IF(NRES.LE.0) NRESL=ITER+1
        CALL BICGSTAB(A,LDA,N,X,B,D,SPARSE_A,ISC_A,ISR_A,ANORM,RESID,
     '    OMEGA,ITER,NRESL,NPRECON,PRECON,FACTORISE,OUTPUTCODE,
     '    ERROR,*9999)
      ELSE IF(SOLVER.EQ.SOLV_GMRES) THEN
        NRESL=NRES
        IF(NRES.LE.0) NRESL=10
        CALL GMRES(A,LDA,N,X,B,D,SPARSE_A,ISC_A,ISR_A,ANORM,RESID,
     '    OMEGA,ITER,NRESL,NPRECON,PRECON,FACTORISE,OUTPUTCODE,
     '    ERROR,*9999)

      ELSE
        ERROR='>>Unknown iterative solver'
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('ITER_SOLVE',ERROR)
      RETURN 1
      END


      SUBROUTINE AMG_ALLOC(LDA,N,NZA,U_PTR,F_PTR,IG_PTR,
     &  A_PTR,IA_PTR,JA_PTR,BDY_PTR,ERROR,*)

C#### Subroutine: AMG_ALLOC
C###  Description:
C###    AMG_ALLOC allocates memory needed to call amg1r6, the John 
C###    Ruge serial version of amg with some OpenMP directives added.
C###  Written by Travis Austin 02/12/04

      IMPLICIT NONE
      INCLUDE 'solver.inc'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER LDA,N,NZA
      CHARACTER ERROR*(*)
      POINTER U_PTR,F_PTR,IG_PTR,A_PTR,IA_PTR,JA_PTR,BDY_PTR
!     Local Variables
      INTEGER NDA, NDU, NDG

      NDA = 4*NZA+5*N
      NDU = int(2.8*N)
      NDG = int(5.4*N)

!     Allocate REAL arrays
      CALL ALLOCATE_MEMORY(NDA,1,DPTYPE,A_PTR,NOINIT,ERROR,*9999) 
      CALL ALLOCATE_MEMORY(NDU,1,DPTYPE,U_PTR,NOINIT,ERROR,*9999) 
      CALL ALLOCATE_MEMORY(NDU,1,DPTYPE,F_PTR,NOINIT,ERROR,*9999) 

!     Allocate INTEGER arrays
      CALL ALLOCATE_MEMORY(NDU,1,INTTYPE, IA_PTR,NOINIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(NDA,1,INTTYPE, JA_PTR,NOINIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(NDG,1,INTTYPE, IG_PTR,NOINIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY( N ,1,INTTYPE,BDY_PTR,NOINIT,ERROR,*9999)

      RETURN

 9999 CALL ERRORS('AMG_ALLOC',ERROR)
      RETURN 1
      END

      SUBROUTINE AMG_FREE(U_PTR,F_PTR,IG_PTR,
     '  A_PTR,IA_PTR,JA_PTR,BDY_PTR,ERROR,*)

C#### Subroutine: AMG_FREE
C###  Description:
C###    AMG_FREE frees memory allocated with AMG_ALLOC.
C###  Written by Travis Austin 02/12/04
 
      IMPLICIT NONE
!     Parameter List
      INTEGER LDA,M,N,SOLVER,WORK1_LEN,WORK2_LEN
      CHARACTER ERROR*(*)
      POINTER U_PTR,F_PTR,IG_PTR,A_PTR,IA_PTR,JA_PTR,BDY_PTR  

      IF(U_PTR.NE.0) CALL FREE_MEMORY(U_PTR,ERROR,*9999)
      IF(F_PTR.NE.0) CALL FREE_MEMORY(F_PTR,ERROR,*9999)
      IF(A_PTR.NE.0) CALL FREE_MEMORY(A_PTR,ERROR,*9999)

      IF(IA_PTR.NE.0) CALL FREE_MEMORY(IA_PTR,ERROR,*9999)
      IF(JA_PTR.NE.0) CALL FREE_MEMORY(JA_PTR,ERROR,*9999)
      IF(IG_PTR.NE.0) CALL FREE_MEMORY(IG_PTR,ERROR,*9999)

      IF(BDY_PTR.NE.0) CALL FREE_MEMORY(BDY_PTR,ERROR,*9999)
      
      RETURN

 9999 CALL ERRORS('AMG_FREE',ERROR)
      RETURN 1
      END

      SUBROUTINE DIR_ALLOC(LDA,M,N,SOLVER,WORK1_PTR,WORK2_PTR,
     '  IWORK_PTR,WORK1_LEN,WORK2_LEN,ERROR,*)

C#### Subroutine: DIR_ALLOC
C###  Description:
C###    DIR_ALLOC allocates memory for the factorised matrix used by 
C###    the LAPACK dense direct solvers.
C###  Written by Stuart Norris 28/05/02

      IMPLICIT NONE
      INCLUDE 'mach00.inc'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'solver.inc'
!     Parameter List
      INTEGER LDA,M,N,SOLVER,WORK1_LEN,WORK2_LEN
      CHARACTER ERROR*(*)
      POINTER WORK1_PTR,WORK2_PTR,IWORK_PTR
!     Local Variables
      INTEGER NB
!     Functions
      INTEGER ILAENV
      EXTERNAL ILAENV


C     LU Decomposition
      IF(SOLVER.EQ.SOLV_LU) THEN
        WORK1_LEN=0
        WORK2_LEN=0
        WORK1_PTR=0
        WORK2_PTR=0
        CALL ALLOCATE_MEMORY(N,1,INTTYPE,IWORK_PTR,NOINIT,ERROR,*9999)

C     SV Decomposition.
      ELSE IF(SOLVER.EQ.SOLV_SVD) THEN
        WORK1_LEN=MIN(M,N)
        WORK2_LEN=N**2
        CALL ALLOCATE_MEMORY(WORK1_LEN,1,DPTYPE,WORK1_PTR,NOINIT,
     '    ERROR,*9999)
        CALL ALLOCATE_MEMORY(WORK2_LEN,1,DPTYPE,WORK2_PTR,NOINIT,
     '    ERROR,*9999)
        IWORK_PTR=0

C     LSQ solve -- find a QR factorisation.
      ELSE IF(SOLVER.EQ.SOLV_LSQ) THEN
        IF (M.GE.N) THEN
          NB=ILAENV(1,'DGEQRF',' ',M,N,-1,-1)
          WORK2_LEN=NB*N
        ELSE
          NB=ILAENV(1,'DGELQF',' ',M,N,-1,-1)
          WORK2_LEN=NB*M
        ENDIF
        WORK1_LEN=MIN(M,N)
        CALL ALLOCATE_MEMORY(WORK1_LEN,1,DPTYPE,WORK1_PTR,NOINIT,
     '    ERROR,*9999)
        CALL ALLOCATE_MEMORY(WORK2_LEN,1,DPTYPE,WORK2_PTR,NOINIT,
     '    ERROR,*9999)
        IWORK_PTR=0

C     Cholesky factorisation
      ELSE IF(SOLVER.EQ.SOLV_CHOLESKY) THEN
        WORK1_LEN=0
        WORK2_LEN=0
        WORK1_PTR=0
        WORK2_PTR=0
        IWORK_PTR=0

C     LDL factorisation
      ELSE IF(SOLVER.EQ.SOLV_LDL) THEN
        WORK1_LEN=0
        WORK2_LEN=0
        WORK1_PTR=0
        WORK2_PTR=0
        IWORK_PTR=0

C     Bzzzzzzt, Error
      ELSE
        ERROR='>>No Such Solver'
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('DIR_ALLOC',ERROR)
      RETURN 1
      END


      SUBROUTINE DIR_FREE(SOLVER,WORK1_PTR,WORK2_PTR,IWORK_PTR,ERROR,*)

C#### Subroutine: DIR_FREE
C###  Description:
C###    DIR_FREE frees memory allocated in DIR_ALLOC
C###  Written by Stuart Norris 04/06/02

      IMPLICIT NONE
      INCLUDE 'solver.inc'
!     Parameter List
      INTEGER SOLVER
      CHARACTER ERROR*(*)
      POINTER WORK1_PTR,WORK2_PTR,IWORK_PTR
!     Local Variables


C     LU Decomposition
      IF(SOLVER.EQ.SOLV_LU) THEN
        IF(IWORK_PTR.NE.0) CALL FREE_MEMORY(IWORK_PTR,ERROR,*9999)

C     SV Decomposition.
      ELSE IF(SOLVER.EQ.SOLV_SVD) THEN
        IF(WORK2_PTR.NE.0) CALL FREE_MEMORY(WORK2_PTR,ERROR,*9999)
        IF(WORK1_PTR.NE.0) CALL FREE_MEMORY(WORK1_PTR,ERROR,*9999)

C     LSQ solve -- find a QR factorisation.
      ELSE IF(SOLVER.EQ.SOLV_LSQ) THEN
        IF(WORK2_PTR.NE.0) CALL FREE_MEMORY(WORK2_PTR,ERROR,*9999)
        IF(WORK1_PTR.NE.0) CALL FREE_MEMORY(WORK1_PTR,ERROR,*9999)

C     Cholesky factorisation
      ELSE IF(SOLVER.EQ.SOLV_CHOLESKY) THEN

C     LDL factorisation
      ELSE IF(SOLVER.EQ.SOLV_LDL) THEN

C     Bzzzzzzt, Error
      ELSE
        ERROR='>>No Such Solver'
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('DIR_FREE',ERROR)
      RETURN 1
      END


      SUBROUTINE DIR_FACTOR(A,LDA,M,N,IWORK_A,SOLVER,WORK1,WORK2,
     '  WORK2_LEN,ANORM,ERROR,*)

C#### Subroutine: DIR_FACTOR
C###  Description:
C###    DIR_FACTOR factorises a system that is to be used in
C###    one of the dense direct solvers
C###  Written by Stuart Norris 04/04/02

      IMPLICIT NONE
      INCLUDE 'mach00.inc'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'solver.inc'
!     Parameter List
      INTEGER LDA,M,N,IWORK_A(*),SOLVER,WORK2_LEN
      REAL*8 A(*),WORK1(*),WORK2(*),ANORM
      CHARACTER ERROR*(*)
!     Local Constants
      INTEGER NORMTYPE,SPARSE_A
      PARAMETER(NORMTYPE=2,SPARSE_A=0)
!     Local variables
      INTEGER INFO,WORK_LEN,LDVT,LDU,ISC_A(1),ISR_A(1)
      POINTER WORK_PTR
      REAL*8 U(1)
      CHARACTER*(*) NAME*6


      INFO=0
      NAME=' '

C     Get the norm of A
      CALL CALC_MAT_NORM(ANORM,A,LDA,N,SPARSE_A,ISC_A,ISR_A,NORMTYPE,
     '  ERROR,*9999)

C     Factorise the system:

C     LU Decomposition
      IF(SOLVER.EQ.SOLV_LU) THEN
        CALL ASSERT(LDA.GE.MAX(1,M),
     '    '>> Incorrect setup for DGETRF LDA >= max(1,M)',
     '    ERROR,*9999)
        CALL ASSERT(LDA.GE.1 .AND. LDA.GE.N,
     '    '>> Incorrect setup for DGETRF LDA >= max(1,N)',
     '    ERROR,*9999)

        CALL DGETRF(M,N,A,LDA,IWORK_A,INFO)
        NAME='DGETRF'

C     SV Decomposition.
      ELSE IF(SOLVER.EQ.SOLV_SVD) THEN
        CALL ASSERT(LDA.GE.MAX(1,M),
     '    '>> Incorrect setup for DGESVD LDA >= max(1,M)',
     '    ERROR,*9999)

        !LWORK=MAX(3*MIN(M,N)+MAX(M,N),5*MIN(M,N)-4)
        WORK_LEN=3*MIN(M,N)+MAX(M,N)+5*MIN(M,N)-4
        LDVT=N
        LDU=M
        WORK_PTR=0
        CALL ALLOCATE_MEMORY(WORK_LEN,1,DPTYPE,WORK_PTR,NOINIT,ERROR,
     '    *9999)

        CALL DGESVD('O','A',M,N,A,LDA,WORK1,U,LDU,WORK2,LDVT,
     '    %VAL(WORK_PTR),WORK_LEN,INFO)
        NAME='DGESVD'

        IF(WORK_PTR.NE.0) CALL FREE_MEMORY(WORK_PTR,ERROR,*9999)

C     LSQ solve -- find a QR factorisation.
      ELSE IF(SOLVER.EQ.SOLV_LSQ) THEN
        IF (M.GE.N) THEN
          CALL DGEQRF(M,N,A,LDA,WORK1,WORK2,WORK2_LEN,INFO)
          NAME='DGEQRF'
        ELSE
          CALL DGELQF(M,N,A,LDA,WORK1,WORK2,WORK2_LEN,INFO)
          NAME='DGELQF'
        ENDIF

C     Cholskey factorisation
      ELSE IF(SOLVER.EQ.SOLV_CHOLESKY) THEN !Cholskey
        CALL ASSERT(LDA.GE.N .AND. N.GE.1,
     '    '>> Incorrect setup for DPOTRF LDA >= max(1,N)',
     '    ERROR,*9999)

        CALL DPOTRF('Upper',N,A,LDA,INFO)
        NAME='DPOTRF'

C     LDL factorisation
      ELSE IF(SOLVER.EQ.SOLV_LDL) THEN !LDL
        CALL ASSERT(LDA.GE.N .AND. N.GE.1,
     '    '>> Incorrect setup for LDLFACTOR LDA >= max(1,N)',
     '    ERROR,*9999)

        CALL LDL_FACTOR(A,LDA,N,ERROR,*9999)

C     Bzzzzzt. Unknown Solver
      ELSE
        WRITE(ERROR,'(''>>Unknown dense solver'')')
        GOTO 9999
      ENDIF

      IF(INFO.NE.0) THEN
        WRITE(ERROR,'('' >>INFO='',I6,'' in  '',A,'''')') 
     '    INFO,NAME
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('DIR_FACTOR',ERROR)
      RETURN 1
      END


      SUBROUTINE DIR_SOLVE(A,LDA,MMAX,M,N,X,B,IWORK_A,SOLVER,
     '  WORK1,WORK2,WORK2_LEN,ANORM,OUTPUTCODE,ERROR,*)

C#### Subroutine: DIR_SOLVE
C###  Description:
C###    DIR_SOLV is a wrapper around the direct solvers, 
C###    allowing a more intelligable naming of arrays etc.
C###  Written by Stuart Norris 04/04/02

      IMPLICIT NONE
      INCLUDE 'mach00.inc'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'solver.inc'
!     Parameter List
      INTEGER IWORK_A(*),LDA,MMAX,M,N,SOLVER,WORK2_LEN,OUTPUTCODE
      REAL*8 A(*),X(*),B(*),WORK1(*),WORK2(*),ANORM
      CHARACTER ERROR*(*)
!     Local parameters
      INTEGER I,INFO,NRHS,WORK_LEN,LDVT
      POINTER WORK_PTR
      REAL*8 SVDTOL,BNORM
      CHARACTER*(*) TRANS*(1),NAME*(6)
!     Functions
      REAL*8 DNRM2
      EXTERNAL DCOPY,DNRM2


      TRANS='N'
      NRHS=1
      LDVT=N
      SVDTOL=1.0D-7
      INFO=0
      NAME=' '

C     LU Decomposition
      IF(SOLVER.EQ.SOLV_LU) THEN
        IF(OUTPUTCODE.GE.2) THEN
          BNORM=DNRM2(N,B,1)
          CALL WRITE_DENSE_SETTINGS(ANORM,BNORM,N,
     '      'Lapack LU Decomposition',ERROR,*9999)
        ENDIF

        CALL DCOPY(M,B,1,X,1)
        CALL DGETRS(TRANS,N,NRHS,A,LDA,IWORK_A,X,MMAX,INFO)
        NAME='DGETRS'

C     SV Decomposition
      ELSE IF(SOLVER.EQ.SOLV_SVD) THEN
        IF(OUTPUTCODE.GE.2) THEN
          BNORM=DNRM2(N,B,1)
          CALL WRITE_DENSE_SETTINGS(ANORM,BNORM,N,'Lapack SVD',
     '      ERROR,*9999)
        ENDIF

        WORK_PTR=0
        WORK_LEN=N
        CALL ALLOCATE_MEMORY(WORK_LEN,1,DPTYPE,WORK_PTR,NOINIT,ERROR,
     '    *9999)

        CALL SVD_SOLVE(M,N,SVDTOL,WORK1,A,LDA,WORK2,LDVT,X,B,
     '    %VAL(WORK_PTR),WORK_LEN,ERROR,*9999)
        NAME='SVD_BA'

        IF(WORK_PTR.NE.0) CALL FREE_MEMORY(WORK_PTR,ERROR,*9999)

C     Least Squares
      ELSE IF(SOLVER.EQ.SOLV_LSQ) THEN
        IF(OUTPUTCODE.GE.2) THEN
          BNORM=DNRM2(N,B,1)
          CALL WRITE_DENSE_SETTINGS(ANORM,BNORM,N,'Lapack LSQ',
     '      ERROR,*9999)
        ENDIF

        IF (M.GE.N) THEN
          WORK_PTR=0
          WORK_LEN=M
          CALL ALLOCATE_MEMORY(WORK_LEN,1,DPTYPE,WORK_PTR,NOINIT,ERROR,
     '      *9999)
          CALL DCOPY(M,B,1,%VAL(WORK_PTR),1)

          CALL DORMQR('Left','Transpose',M,1,N,A,LDA,WORK1,
     '      %VAL(WORK_PTR),M,WORK2,WORK2_LEN,INFO)
          NAME='DORMQR'
          IF(INFO.EQ.0) THEN
            CALL DCOPY(N,%VAL(WORK_PTR),1,X,1)
            CALL DTRSM('Left','Upper','No transpose','Non-unit',N,1,
     '        1.0D0,A,LDA,X,M)
          ENDIF

          IF(WORK_PTR.NE.0) CALL FREE_MEMORY(WORK_PTR,ERROR,*9999)
        ELSE
          CALL DCOPY(M,B,1,X,1)
          CALL DTRSM('Left','Lower','No transpose','Non-unit',M,1,
     '      1.0D0,A,LDA,X,M)
          DO I=M+1,N
            X(I)=0.0D0
          ENDDO
          CALL DORMLQ('Left','Transpose',N,1,M,A,LDA,WORK1,X,N,WORK2,
     '      WORK2_LEN,INFO)
          NAME='DORMLQ'
        ENDIF

C     Cholesky Factorisation
      ELSE IF(SOLVER.EQ.SOLV_CHOLESKY) THEN
        IF(OUTPUTCODE.GE.2) THEN
          BNORM=DNRM2(N,B,1)
          CALL WRITE_DENSE_SETTINGS(ANORM,BNORM,N,
     '      'Lapack Cholesky Factorisation',ERROR,*9999)
        ENDIF

        CALL DCOPY(M,B,1,X,1)
        CALL DPOTRS('Upper',N,NRHS,A,LDA,X,MMAX,INFO)
        NAME='DPOTRS'

C     LDL Decomposition
      ELSE IF(SOLVER.EQ.SOLV_LDL) THEN
        IF(OUTPUTCODE.GE.2) THEN
          BNORM=DNRM2(N,B,1)
          CALL WRITE_DENSE_SETTINGS(ANORM,BNORM,N,
     '      'Dense LDL Decomposition',ERROR,*9999)
        ENDIF

        CALL LDL_SOLVE(A,LDA,N,X,B,ERROR,*9999)

C     Bzzzzzt. Unknown Solver
      ELSE
        WRITE(ERROR,'(''>>Unknown direct solver'')')
        GOTO 9999
      ENDIF

C     Catch errors from LAPACK routines
      IF(INFO.NE.0) THEN
        WRITE(ERROR,'('' >>INFO='',I6,'' in  '',A,'''')') INFO,NAME
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('DIR_SOLVE',ERROR)
      RETURN 1
      END


      SUBROUTINE WRITE_DENSE_SETTINGS(ANORM,BNORM,N,NAME,ERROR,*)

C#### Subroutine: WRITE_DENSE_SETTINGS
C###  Description:
C###    Write out the solver settings for the LAPACK solvers.
C###  Written by Stuart Norris 07/08/02

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER N
      REAL*8 ANORM,BNORM
      CHARACTER NAME*(*),ERROR*(*)
!     Local Constants
      INTEGER SPARSE_A
      PARAMETER(SPARSE_A=0)
!     Local variables
      INTEGER NZA
!     Functions
      INTEGER LEN_TRIM
      EXTERNAL LEN_TRIM


      NZA=N**2

      WRITE(OP_STRING,'(A)') NAME(:LEN_TRIM(NAME))
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'()')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      CALL WRITE_MATRIX_DIMS(SPARSE_A,N,NZA,ANORM,BNORM,ERROR,*9999)

      RETURN

 9999 CALL ERRORS('WRITE_DENSE_SETTINGS',ERROR)
      RETURN 1
      END
