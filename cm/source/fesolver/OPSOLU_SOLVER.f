      SUBROUTINE OPSOLU_SOLVER(nx,ERROR,*)

C#### Subroutine: OPSOLU_SOLVER
C###  Description:
C###    Outputs solution parameters for solvers in libsolver.

      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='OPSOLU_SOLVER')
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'solv00.cmn'
      INCLUDE 'solver.inc'
!     Parameter List
      INTEGER nx
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER SLU_CO
      CHARACTER SOLTITLE(24)*28,PRECONTITLE(7)*30,
     '  SPARSETITLE(0:5)*32,SUPERLUCOLORDERTITLE(0:3)*34

      DATA SOLTITLE /'LU Decomposition            ',
     '                'Single Value Decomposition  ',
     '                'Least Squares               ',
     '                'Cholesky Decomposition      ',
     '                'Jacobi Iteration            ',
     '                'Succesive Over Relaxation   ',
     '                'Incomplete LU(0)            ',
     '                'Incomplete LU(1)            ',
     '                'Conjugate Gradient          ',
     '                'Biconjugate Gradient Stab   ',
     '                'Generalised Minimum Residual',
     '                ' ',
     '                ' ',
     '                ' ',
     '                ' ',
     '                ' ',
     '                ' ',
     '                ' ',
     '                ' ',
     '                ' ',
     '                'SuperLU Sparse LU',
     '                'Umfpack 4.0 Sparse LU',
     '                'Umfpack 2.2 Sparse LU',
     '                'Harwell MA28 Sparse LU'/

      DATA SUPERLUCOLORDERTITLE /'Natural ordering                  ',
     '                            'Minimum degree ordering from A^T.A',
     '                            'Minimum degree ordering from A+A^T',
     '                            'Column Minimum degree ordering'/

      DATA PRECONTITLE /'None                          ',
     '                   'Point Jacobi                  ',
     '                   'Jacobi Iteration              ',
     '                   'Symmetric succesive relaxation',
     '                   'Incomplete LU Decomposition(0)',
     '                   'Incomplete LU Decomposition(1)',
     '                   'Row scale                     '/

      DATA SPARSETITLE/'fully populated matrix          ',
     '                  'compressed row sparse matrix    ',
     '                  'row-column sparse matrix        ',
     '                  'compressed column sparse matrix ',
     '                  'sorted row-column sparse matrix ',
     '                  'umfpack row-column sparse matrix'/

      CALL ENTERS(ROUTINENAME,*9999)

      WRITE(OP_STRING,'(''  Solution procedure is '',A)')
     '  SOLTITLE(SOLVEROPTION(nx))
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''  Solution matrix is stored as a '',A)')
     '  SPARSETITLE(SPARSEGKK(nx))
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      IF(SOLVEROPTION(nx).EQ.SOLV_SUPERLU) THEN !LU-sparse
        CALL SUPERLU_GETPARAM('COLUMN ORDERING',SLU_CO,
     '    SOLV_SLU_PARAM(1,NX),ERROR,*9999)
        WRITE(OP_STRING,'(''  Column ordering is '',A)')
     '    SUPERLUCOLORDERTITLE(SLU_CO)
      ELSE IF(SOLVEROPTION(nx).EQ.SOLV_BICGSTAB.OR.
     '    SOLVEROPTION(nx).EQ.SOLV_CG .OR.
     '    SOLVEROPTION(nx).EQ.SOLV_GMRES) THEN
        WRITE(OP_STRING,'(''  Maximum number of iterations '
     '    //'for the iterative solver is '',I6)') SOLV_MAXIT(nx)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''  Krylov refactorisation frequency '
     '    //' is '',I3)') SOLV_NRES(nx)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''  Iterative solver tolerance is '','
     '    //'D11.4)') SOLV_TOL(nx)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''  Type of preconditioning is '','
     '    //'A)') PRECONTITLE(PRECON_CODE(nx))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS(ROUTINENAME)
      RETURN
 9999 CALL ERRORS(ROUTINENAME,ERROR)
      CALL EXITS(ROUTINENAME)
      RETURN 1
      END


