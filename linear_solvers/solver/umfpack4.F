
      SUBROUTINE UMFPACK4_RESET(PARAM,ERROR,*)

C#### Subroutine: UMFPACK4_RESET
C###  Description:
C###    Initialise the Umfpack superLU parameter array to the default values
C###  Written by Stuart Norris 31/07/02

      IMPLICIT NONE
!     Parameter List
      REAL*8 PARAM(40)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I


      DO I=1,40
        PARAM(I)=-1.0D0
      ENDDO

      RETURN

 9999 CALL ERRORS('UMFPACK4_RESET',ERROR)
      RETURN 1
      END


      SUBROUTINE UMFPACK4_FACTOR(A,N,NZA,ISC_A,ISR_A,SYMBOLIC_PTR,
     '  NUMERIC_PTR,UMFDEF,CONTROL,INFO,ANORM,ERROR,*)

C#### Subroutine: UMFPACK4_FACTOR
C###  Description:
C###    Factorise a system.
C###  Written by Stuart Norris 11/07/02

      IMPLICIT NONE
      INCLUDE 'mach00.inc'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'solver.inc'
!     Parameter List
      INTEGER N,NZA,ISC_A(*),ISR_A(*)
      REAL*8 A(*),UMFDEF(*),CONTROL(20),INFO(90),ANORM
      CHARACTER ERROR*(*)
      POINTER NUMERIC_PTR,SYMBOLIC_PTR
!     Local Constants
      INTEGER NORMTYPE,SPARSE_A
      PARAMETER(NORMTYPE=2,SPARSE_A=1)
!     Local Variables
      INTEGER IERR
      CHARACTER*32 NAME,BUFF


      IF(SPARSE_A.NE.1) THEN
        ERROR='>>Only able to solve compressed row matricies'
        GOTO 9999
      ENDIF

C     Get the norm of A
      CALL CALC_MAT_NORM(ANORM,A,N,N,SPARSE_A,ISC_A,ISR_A,NORMTYPE,
     '  ERROR,*9999)

C     Factor the system
      CALL UMFPACK4_FACTOR_C(A,N,NZA,ISC_A,ISR_A,SYMBOLIC_PTR,
     '  NUMERIC_PTR,UMFDEF,CONTROL,INFO,NAME,IERR)

      IF(IERR.NE.0) THEN
        CALL UMFPACK4_ERROR(IERR,BUFF)
        ERROR='>>'//BUFF
        GOTO 9999        
      ENDIF

      RETURN

 9999 CALL ERRORS('UMFPACK4_FACTOR',ERROR)
      RETURN 1
      END


      SUBROUTINE UMFPACK4_REFACTOR(A,N,NZA,ISC_A,ISR_A,SYMBOLIC_PTR,
     '  NUMERIC_PTR,CONTROL,INFO,ANORM,ERROR,*)

C#### Subroutine: UMFPACK4_REFACTOR
C###  Description:
C###    Factorise a system using a previously supplied sparsity pattern.
C###  Written by Stuart Norris 11/07/02

      IMPLICIT NONE
      INCLUDE 'mach00.inc'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'solver.inc'
!     Parameter List
      INTEGER N,NZA,ISC_A(*),ISR_A(*)
      REAL*8 A(*),CONTROL(20),INFO(90),ANORM
      CHARACTER ERROR*(*)
      POINTER NUMERIC_PTR,SYMBOLIC_PTR
!     Local Constants
      INTEGER NORMTYPE,SPARSE_A
      PARAMETER(NORMTYPE=2,SPARSE_A=1)
!     Local Variables
      INTEGER IERR
      CHARACTER*32 NAME,BUFF


      IF(SPARSE_A.NE.1) THEN
        ERROR='>>Only able to solve compressed row matricies'
        GOTO 9999
      ENDIF

C     Get the norm of A
      CALL CALC_MAT_NORM(ANORM,A,N,N,SPARSE_A,ISC_A,ISR_A,NORMTYPE,
     '  ERROR,*9999)

C     Factor the system using the previous sparsity pattern
      CALL UMFPACK4_REFACTOR_C(A,N,NZA,ISC_A,ISR_A,SYMBOLIC_PTR,
     '  NUMERIC_PTR,CONTROL,INFO,NAME,IERR)

      IF(IERR.NE.0) THEN
        CALL UMFPACK4_ERROR(IERR,BUFF)
        ERROR='>>'//BUFF
        GOTO 9999        
      ENDIF

      RETURN

 9999 CALL ERRORS('UMFPACK4_REFACTOR',ERROR)
      RETURN 1
      END


      SUBROUTINE UMFPACK4_SOLVE(A,N,NZA,ISC_A,ISR_A,B,X,NUMERIC_PTR,
     '  CONTROL,INFO,ANORM,OUTPUTCODE,ERROR,*)

C#### Subroutine: UMFPACK4_SOLVE
C###  Description:
C###    Solve a system that has been factorised using UMFPACK4_FACTOR
C###    or UMFPACK4_REFACTOR.
C###  Written by Stuart Norris 11/07/02

      IMPLICIT NONE
!     Parameter List
      INTEGER N,NZA,ISC_A(*),ISR_A(*),OUTPUTCODE
      REAL*8 A(*),B(*),X(*),CONTROL(20),INFO(90),ANORM
      POINTER NUMERIC_PTR
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IERR
      REAL*8 BNORM
      CHARACTER*80 BUFF
!     External functions
      REAL*8 DNRM2
      EXTERNAL DNRM2


      IF(OUTPUTCODE.GE.2) THEN
        BNORM=DNRM2(N,B,1)

        CALL WRITE_UMFPACK_SETTINGS(N,NZA,ANORM,BNORM,CONTROL,
     '    ERROR,*9999)
      ENDIF

      CALL UMFPACK4_SOLVE_C(A,N,NZA,ISC_A,ISR_A,B,X,NUMERIC_PTR,
     '  CONTROL,INFO,IERR)

      IF(IERR.NE.0) THEN
        CALL UMFPACK4_ERROR(IERR,BUFF)
        ERROR='>>'//BUFF
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('UMFPACK4_SOLVE',ERROR)
      RETURN 1
      END


      SUBROUTINE UMFPACK4_FREE(SYMBOLIC_PTR,FREE_SYMBOLIC,
     '  NUMERIC_PTR,FREE_NUMERIC,ERROR,*)

C#### Subroutine: UMFPACK4_FREE
C###  Description:
C###    Free memory used by the Umfpack 4.0 solver.
C###  Written by Stuart Norris 11/07/02

      IMPLICIT NONE
!     Parameter List
      LOGICAL FREE_SYMBOLIC,FREE_NUMERIC
      POINTER SYMBOLIC_PTR,NUMERIC_PTR
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IFREE_SYMBOLIC,IFREE_NUMERIC


      IF(FREE_SYMBOLIC) THEN
        IFREE_SYMBOLIC=1
      ELSE
        IFREE_SYMBOLIC=0
      ENDIF
      IF(FREE_NUMERIC) THEN
        IFREE_NUMERIC=1
      ELSE
        IFREE_NUMERIC=0
      ENDIF
      CALL UMFPACK4_FREE_C(SYMBOLIC_PTR,IFREE_SYMBOLIC,
     '  NUMERIC_PTR,IFREE_NUMERIC)

      RETURN

 9999 CALL ERRORS('UMFPACK4_FREE',ERROR)
      RETURN 1
      END


      SUBROUTINE WRITE_UMFPACK_SETTINGS(N,NZA,ANORM,BNORM,
     '  CONTROL,ERROR,*)

C#### Subroutine: WRITE_UMFPACK_SETTINGS
C###  Description:
C###    Write out the solver settings for the Umfpack solver.
C###  Written by Stuart Norris 07/08/02

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER N,NZA
      REAL*8 ANORM,BNORM,CONTROL(*)
      CHARACTER ERROR*(*)
!     Local Constants
      INTEGER SPARSE_A
      PARAMETER(SPARSE_A=1)


      WRITE(OP_STRING,'('' Umfpack 4.0'')')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'()')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,7000) ' Pivot Tolerence    : ',CONTROL(4)
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'()')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      CALL WRITE_MATRIX_DIMS(SPARSE_A,N,NZA,ANORM,BNORM,ERROR,*9999)

      RETURN

 7000 FORMAT(A,1P,E11.4)

 9999 CALL ERRORS('WRITE_UMFPACK_SETTINGS',ERROR)
      RETURN 1
      END
