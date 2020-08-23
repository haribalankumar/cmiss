
      SUBROUTINE SUPERLU_RESET(PARAM,ERROR,*)

C#### Subroutine: SUPERLU_RESET
C###  Description:
C###    Initialise the superLU parameter array to the default values
C###  Written by Stuart Norris 04/06/02

      IMPLICIT NONE
!     Parameter List
      INTEGER PARAM(10)
      CHARACTER ERROR*(*)
!     Local Variables


      CALL SUPERLU_RESETP(PARAM)

      RETURN

 9999 CALL ERRORS('SUPERLU_RESET',ERROR)
      RETURN 1
      END


      SUBROUTINE SUPERLU_SETPARAM(LABEL,VALUE,PARAM,ERROR,*)

C#### Subroutine: SUPERLU_SETPARAM
C###  Description:
C###    Set superLU parameters
C###  Written by Stuart Norris 04/06/02

      IMPLICIT NONE
!     Parameter List
      INTEGER VALUE,PARAM(10)
      CHARACTER LABEL*(*),ERROR*(*)
!     Local Variables


      IF(LABEL.EQ.'PANEL SIZE') THEN
        PARAM(1)=VALUE
      ELSE IF(LABEL.EQ.'RELAX PARAMETER') THEN
        PARAM(2)=VALUE
      ELSE IF(LABEL.EQ.'SUPERNODE SIZE') THEN
        PARAM(3)=VALUE
      ELSE IF(LABEL.EQ.'MIN ROW DIM') THEN
        PARAM(4)=VALUE
      ELSE IF(LABEL.EQ.'MIN COL DIM') THEN
        PARAM(5)=VALUE
      ELSE IF(LABEL.EQ.'ESTIMATED FILL') THEN
        PARAM(6)=VALUE
        PARAM(7)=VALUE
        PARAM(8)=VALUE
      ELSE IF(LABEL.EQ.'ESTIMATED L FILL') THEN
        PARAM(6)=VALUE
      ELSE IF(LABEL.EQ.'ESTIMATED U FILL') THEN
        PARAM(7)=VALUE
      ELSE IF(LABEL.EQ.'ESTIMATED L SUB') THEN
        PARAM(8)=VALUE
      ELSE IF(LABEL.EQ.'COLUMN ORDERING') THEN
        PARAM(9)=VALUE
      ELSE IF(LABEL.EQ.'NUMBER OF PROCS') THEN
        PARAM(10)=VALUE
      ELSE
        ERROR='>>No such parameter '''//LABEL//''''
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('SUPERLU_SETPARAM',ERROR)
      RETURN 1
      END


      SUBROUTINE SUPERLU_GETPARAM(LABEL,VALUE,PARAM,ERROR,*)

C#### Subroutine: SUPERLU_GETPARAM
C###  Description:
C###    Read superLU parameters
C###  Written by Stuart Norris 31/07/02

      IMPLICIT NONE
!     Parameter List
      INTEGER VALUE,PARAM(10)
      CHARACTER LABEL*(*),ERROR*(*)
!     Local Variables


      IF(LABEL.EQ.'PANEL SIZE') THEN
        VALUE=PARAM(1)
      ELSE IF(LABEL.EQ.'RELAX PARAMETER') THEN
        VALUE=PARAM(2)
      ELSE IF(LABEL.EQ.'SUPERNODE SIZE') THEN
        VALUE=PARAM(3)
      ELSE IF(LABEL.EQ.'MIN ROW DIM') THEN
        VALUE=PARAM(4)
      ELSE IF(LABEL.EQ.'MIN COL DIM') THEN
        VALUE=PARAM(5)
      ELSE IF(LABEL.EQ.'ESTIMATED FILL') THEN
        VALUE=PARAM(6)
      ELSE IF(LABEL.EQ.'ESTIMATED L FILL') THEN
        VALUE=PARAM(6)
      ELSE IF(LABEL.EQ.'ESTIMATED U FILL') THEN
        VALUE=PARAM(7)
      ELSE IF(LABEL.EQ.'ESTIMATED L SUB') THEN
        VALUE=PARAM(8)
      ELSE IF(LABEL.EQ.'COLUMN ORDERING') THEN
        VALUE=PARAM(9)
      ELSE IF(LABEL.EQ.'NUMBER OF PROCS') THEN
        VALUE=PARAM(10)
      ELSE
        ERROR='>>No such parameter '''//LABEL//''''
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('SUPERLU_GETPARAM',ERROR)
      RETURN 1
      END


      SUBROUTINE SUPERLU_FACTOR(A,N,NZA,ISC_A,ISR_A,SYMBOLIC_PTR,
     '  PARAM,ANORM,ERROR,*)

C#### Subroutine: SUPERLU_FACTOR
C###  Description:
C###    Factorise a sparse array using superLU
C###  Written by Stuart Norris 04/06/02

      IMPLICIT NONE
!     Parameter List
      INTEGER N,NZA,ISC_A(*),ISR_A(*),PARAM(10)
      REAL*8 A(*),ANORM
      CHARACTER ERROR*(*)
      POINTER SYMBOLIC_PTR
!     Local Constants
      INTEGER NORMTYPE,SPARSE_A
      PARAMETER(NORMTYPE=2,SPARSE_A=1)
!     Local Variables
      INTEGER INFO


C     Get the norm of A
      CALL CALC_MAT_NORM(ANORM,A,N,N,SPARSE_A,ISC_A,ISR_A,NORMTYPE,
     '  ERROR,*9999)

C     Factorise
      INFO=0
      CALL SUPERLU_DGSTRF(N,NZA,A,ISR_A,ISC_A,SYMBOLIC_PTR,PARAM,INFO)
      CALL SUPERLU_ERROR(INFO,N,ERROR,*9999)

      RETURN

 9999 CALL ERRORS('SUPERLU_FACTOR',ERROR)
      RETURN 1
      END


      SUBROUTINE SUPERLU_REFACTOR(A,N,NZA,ISC_A,ISR_A,SYMBOLIC_PTR,
     '  PARAM,ANORM,ERROR,*)

C#### Subroutine: SUPERLU_FACTOR
C###  Description:
C###    Refactorise a sparse array using superLU
C###  Written by Stuart Norris 04/06/02

      IMPLICIT NONE
!     Parameter List
      INTEGER N,NZA,ISC_A(*),ISR_A(*),PARAM(10)
      REAL*8 A(*),ANORM
      CHARACTER ERROR*(*)
      POINTER SYMBOLIC_PTR
!     Local Constants
      INTEGER NORMTYPE,SPARSE_A
      PARAMETER(NORMTYPE=2,SPARSE_A=1)
!     Local Variables
      INTEGER INFO


C     Get the norm of A
      CALL CALC_MAT_NORM(ANORM,A,N,N,SPARSE_A,ISC_A,ISR_A,NORMTYPE,
     '  ERROR,*9999)

C     Re-Factorise
      INFO=0
      CALL SUPERLU_DGSTRF_REFACT(N,NZA,A,ISR_A,ISC_A,SYMBOLIC_PTR,
     '  PARAM,INFO)
      CALL SUPERLU_ERROR(INFO,N,ERROR,*9999)

      RETURN

 9999 CALL ERRORS('SUPERLU_REFACTOR',ERROR)
      RETURN 1
      END


      SUBROUTINE SUPERLU_ERROR(INFO,N,ERROR,*)

C#### Subroutine: SUPERLU_ERROR
C###  Description:
C###    Interpret error codes from the superLU solver
C###  Written by Stuart Norris 04/06/02

      IMPLICIT NONE
!     Parameter List
      INTEGER INFO,N
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER M


      IF(INFO.LT.0) THEN
        M=ABS(INFO)
        WRITE(ERROR,7000) '>>The ',M,'th argument had an illegal value'
        GOTO 9999
      ELSE IF(INFO.GT.0 .AND. INFO.LT.N) THEN
        WRITE(ERROR,7000) '>>The matrix is singular at the ',INFO,
     '    'th column'
        GOTO 9999
      ELSE IF(INFO.GT.N) THEN
        M=INFO-N
        WRITE(ERROR,7000) '>>Malloc failure with ',M,' bytes allocated'
        GOTO 9999
      ENDIF

 7000  FORMAT(A,I10,A)

      RETURN

 9999 CALL ERRORS('SUPERLU_ERROR',ERROR)
      RETURN 1
      END


      SUBROUTINE SUPERLU_SOLVE(B,X,N,NZA,SYMBOLIC_PTR,PARAM,ANORM,
     '  OUTPUTCODE,ERROR,*)

C#### Subroutine: SUPERLU_SOLVE
C###  Description:
C###    Solve a sparse system using superLU
C###  Written by Stuart Norris 04/06/02

      IMPLICIT NONE
!     Parameter List
      INTEGER N,NZA,PARAM(10),OUTPUTCODE
      REAL*8 B(*),X(*),ANORM
      CHARACTER ERROR*(*)
      POINTER SYMBOLIC_PTR
!     Local Variables
      INTEGER INFO
      REAL*8 BNORM
!     External functions
      REAL*8 DNRM2
      EXTERNAL DNRM2


      IF(OUTPUTCODE.GE.2) THEN
        BNORM=DNRM2(N,B,1)

        CALL WRITE_SUPERLU_SETTINGS(N,NZA,ANORM,BNORM,PARAM,ERROR,*9999)
      ENDIF

      CALL DCOPY(N,B,1,X,1)
      INFO=0
      CALL SUPERLU_DGSTRS('T',N,X,1,N,SYMBOLIC_PTR,PARAM,INFO)
      CALL SUPERLU_ERROR(INFO,N,ERROR,*9999)

      RETURN

 9999 CALL ERRORS('SUPERLU_SOLVE',ERROR)
      RETURN 1
      END


      SUBROUTINE SUPERLU_FREE(SYMBOLIC_PTR,ERROR,*)

C#### Subroutine: SUPERLU_FREE
C###  Description:
C###    Free memory used by the superLU solver factorisation
C###  Written by Stuart Norris 04/06/02

      IMPLICIT NONE
!     Parameter List
      POINTER SYMBOLIC_PTR
      CHARACTER ERROR*(*)
!     Local Variables


      CALL SUPERLU_DESTROY(SYMBOLIC_PTR)
        
      RETURN

 9999 CALL ERRORS('SUPERLU_FREE',ERROR)
      RETURN 1
      END


      SUBROUTINE WRITE_SUPERLU_SETTINGS(N,NZA,ANORM,BNORM,PARAM,ERROR,*)

C#### Subroutine: WRITE_SUPERLU_SETTINGS
C###  Description:
C###    Write out the solver settings for the SuperLU solver.
C###  Written by Stuart Norris 07/08/02

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER N,NZA,PARAM(10)
      REAL*8 ANORM,BNORM
      CHARACTER ERROR*(*)
!     Local Constants
      INTEGER SPARSE_A
      PARAMETER(SPARSE_A=1)
!     Local Variables
      INTEGER VAL,VERSION
      LOGICAL MULTI


      MULTI=.FALSE.
C$    MULTI=.TRUE.
      VERSION=2
C$    VERSION=1

      IF(VERSION.EQ.2) THEN
        WRITE(OP_STRING,'('' SuperLU 2.0 Solver'')')
      ELSE
        WRITE(OP_STRING,'('' SuperLU 1.0 Solver'')')
      ENDIF
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'()')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      CALL SUPERLU_GETPARAM('COLUMN ORDERING',VAL,PARAM,ERROR,*9999)
      WRITE(OP_STRING,7000) ' Column ordering    : ',VAL
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      CALL SUPERLU_GETPARAM('PANEL SIZE',VAL,PARAM,ERROR,*9999)
      WRITE(OP_STRING,7000) ' Panel Size         : ',VAL
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      CALL SUPERLU_GETPARAM('RELAX PARAMETER',VAL,PARAM,ERROR,*9999)
      WRITE(OP_STRING,7000) ' Relaxation factor  : ',VAL
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      CALL SUPERLU_GETPARAM('SUPERNODE SIZE',VAL,PARAM,ERROR,*9999)
      WRITE(OP_STRING,7000) ' Max supernode size : ',VAL
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      CALL SUPERLU_GETPARAM('MIN ROW DIM',VAL,PARAM,ERROR,*9999)
      WRITE(OP_STRING,7000) ' Minimum row dim    : ',VAL
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      CALL SUPERLU_GETPARAM('MIN COL DIM',VAL,PARAM,ERROR,*9999)
      WRITE(OP_STRING,7000) ' Minimum col dim    : ',VAL
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      IF(.NOT.MULTI) THEN
        CALL SUPERLU_GETPARAM('ESTIMATED FILL',VAL,PARAM,ERROR,*9999)
        WRITE(OP_STRING,7000) ' Estimated Fill     : ',VAL
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ELSE
        CALL SUPERLU_GETPARAM('ESTIMATED L FILL',VAL,PARAM,ERROR,*9999)
        WRITE(OP_STRING,7000) ' Estimated L Fill   : ',VAL
       CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
       CALL SUPERLU_GETPARAM('ESTIMATED U FILL',VAL,PARAM,ERROR,*9999)
       WRITE(OP_STRING,7000) ' Estimated U Fill   : ',VAL
       CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
       CALL SUPERLU_GETPARAM('ESTIMATED L SUB',VAL,PARAM,ERROR,*9999)
       WRITE(OP_STRING,7000) ' Estimated L Sub    : ',VAL
       CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF
      WRITE(OP_STRING,'()')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      CALL WRITE_MATRIX_DIMS(SPARSE_A,N,NZA,ANORM,BNORM,ERROR,*9999)

      RETURN

 7000 FORMAT(A,5X,I5)

 9999 CALL ERRORS('WRITE_SUPERLU_SETTINGS',ERROR)
      RETURN 1
      END
