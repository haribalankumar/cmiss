
      SUBROUTINE MA28AD_ALLOC(N,NZA,NICN,NIRN,LICN,LIRN,D_PTR,WORK_PTR,
     '  IKEEP_PTR,ISC2_PTR,ERROR,*)

C#### Subroutine: MA28AD_ALLOC
C###  Description:
C###    MA28AD_ALLOC allocates memory for the Harwell MA28AD() sparse
C###    direct linear solver.
C###  Written by Stuart Norris 03/06/02

      IMPLICIT NONE
      INCLUDE 'mach00.inc'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'solver.inc'
!     Parameter List
      INTEGER N,NZA,NICN,NIRN,LICN,LIRN
      POINTER D_PTR,WORK_PTR,IKEEP_PTR,ISC2_PTR
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER LKEEP


      IF(MAXINT/NICN.LT.NZA) THEN
        LICN=MAXINT
      ELSE
        LICN=NICN*NZA
      ENDIF
      IF(MAXINT/NIRN.LT.NZA) THEN
        LIRN=MAXINT
      ELSE
        LIRN=NIRN*NZA
      ENDIF
      IF(MAXINT/5.LT.N) THEN
        WRITE(ERROR,7000) '>>Cannot alloc 5*',N,' integers'
        GOTO 9999
      ELSE
        LKEEP=5*N
      ENDIF

      CALL ALLOCATE_MEMORY(    N,1, DPTYPE, WORK_PTR,NOINIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(LKEEP,1,INTTYPE,IKEEP_PTR,NOINIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY( LICN,1,INTTYPE, ISC2_PTR,NOINIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY( LICN,1, DPTYPE,    D_PTR,NOINIT,ERROR,*9999)

      RETURN

 7000 FORMAT(A,I20,A)

 9999 CALL ERRORS('MA28AD_ALLOC',ERROR)
      RETURN 1
      END


      SUBROUTINE MA28AD_FREE(D_PTR,WORK_PTR,IKEEP_PTR,ISC2_PTR,ERROR,*)

C#### Subroutine: MA28AD_FREE
C###  Description:
C###    MA28AD_FREE frees memory allocated in MA28AD_ALLOC.
C###  Written by Stuart Norris 04/06/02

      IMPLICIT NONE
      POINTER D_PTR,WORK_PTR,IKEEP_PTR,ISC2_PTR
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables


      IF(D_PTR.NE.0)     CALL FREE_MEMORY(    D_PTR,ERROR,*9999)
      IF(WORK_PTR.NE.0)  CALL FREE_MEMORY( WORK_PTR,ERROR,*9999)
      IF(IKEEP_PTR.NE.0) CALL FREE_MEMORY(IKEEP_PTR,ERROR,*9999)
      IF(ISC2_PTR.NE.0)  CALL FREE_MEMORY( ISC2_PTR,ERROR,*9999)

      RETURN

 9999 CALL ERRORS('MA28AD_FREE',ERROR)
      RETURN 1
      END


      SUBROUTINE MA28AD_FACTOR(A,D,N,NZA,LICN,LIRN,ISC_A,ISR_A,USET,
     '  IKEEP,ISC2_A,WORK,ERROR,*)

C#### Subroutine: MA28AD_FACTOR
C###  Description:
C###    Factorise a sparse system using Harwell's MA28AD LU code
C###  Written by Stuart Norris 08/05/02

      IMPLICIT NONE
      INCLUDE 'mach00.inc'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'solver.inc'
!     Parameter List
      INTEGER N,NZA,LICN,LIRN,ISR_A(*),ISC_A(*),IKEEP(*),ISC2_A(*)
      REAL*8 A(*),D(*),WORK(*),USET
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IERR,LWORK
      POINTER IWORK_PTR,ISR2_PTR


      IF(MAXINT/8.LT.N) THEN
        WRITE(ERROR,7000) '>>Cannot alloc 8*',N,' integers'
        GOTO 9999
      ELSE
        LWORK=8*N
      ENDIF

C     Malloc and copy arrays
      CALL ALLOCATE_MEMORY(LWORK,1,INTTYPE,IWORK_PTR,NOINIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY( LIRN,1,INTTYPE, ISR2_PTR,NOINIT,ERROR,*9999)

C     CALL ICOPY(NZA,ISC_A,1,ISC2_A,1)
C     CALL ICOPY(NZA,ISR_A,1,%VAL(ISR2_PTR),1)
C     CALL DCOPY(NZA,A,1,D,1)
      CALL MA28AD_SETUP(A,D,N,NZA,ISC_A,ISR_A,ISC2_A,%VAL(ISR2_PTR),
     '  ERROR,*9999)

C     Factorise
      CALL MA28AD(N,NZA,D,LICN,%VAL(ISR2_PTR),LIRN,ISC2_A,USET,IKEEP,
     '  %VAL(IWORK_PTR),WORK,IERR)
      CALL MA28AD_ERR(IERR,LIRN,LICN,ERROR,*9999)
      IF(IERR.GT.0) THEN
        WRITE(OP_STRING,'(A)') ERROR
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ELSE IF(IERR.LT.0) THEN
        GOTO 9999
      ENDIF

C     Clear memory
      IF(ISR2_PTR.NE.0)  CALL FREE_MEMORY(ISR2_PTR,ERROR,*9999)
      IF(IWORK_PTR.NE.0) CALL FREE_MEMORY(IWORK_PTR,ERROR,*9999)

      RETURN

 7000 FORMAT(A,I20,A)

 9999 CALL ERRORS('MA28AD_FACTOR',ERROR)
      RETURN 1
      RETURN
      END


      SUBROUTINE MA28AD_SOLVE(A,B,X,N,NZA,LICN,IKEEP,ISC2_A,WORK,
     '  ERROR,*)

C#### Subroutine: MA28AD_SOLVE
C###  Description:
C###    Solve a system that has been factorised with Harwell's MA28AD
C###  Written by Stuart Norris 08/05/02

      IMPLICIT NONE
!     Parameter List
      INTEGER N,NZA,LICN,IKEEP(*),ISC2_A(*)
      REAL*8 A(*),B(*),X(*),WORK(*)
      CHARACTER ERROR*(*)
!     Local Variables


      CALL DCOPY(N,B,1,X,1)
      CALL MA28CD(N,A,LICN,ISC2_A,IKEEP,X,WORK,1)

      RETURN

 9999 CALL ERRORS('MA28AD_SOLVE',ERROR)
      RETURN 1
      RETURN
      END


      SUBROUTINE MA28AD_SETUP(A,D,N,NZA,ISC_A,ISR_A,ISC2_A,ISR2_A,
     '  ERROR,*)

C#### Subroutine: MA28AD_SETUP
C###  Description:
C###    Set up data arrays for the Harwell MA28AD() sparse LU solver.
C###  Written by Stuart Norris 02/05/02

      IMPLICIT NONE
!     Parameter List
      INTEGER N,NZA,ISC_A(*),ISR_A(*),ISC2_A(*),ISR2_A(*)
      REAL*8 A(*),D(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,IJ


      CALL DCOPY(NZA,A,1,D,1)
      CALL ICOPY(NZA,ISC_A,1,ISC2_A,1)

      DO I=1,N
        DO IJ=ISR_A(I),ISR_A(I+1)-1
          ISR2_A(IJ)=I
        ENDDO
      ENDDO

      RETURN

 9999 CALL ERRORS('MA28AD_SETUP',ERROR)
      RETURN 1
      END


      SUBROUTINE MA28AD_ERR(IERR,LIRN,LICN,ERR,*)

C#### Subroutine: MA28AD_ERR
C###  Description:
C###    Print errors from Harwell's MA28AD solver
C###  Written by Stuart Norris 08/05/02

      IMPLICIT NONE
!     Parameter List
      INTEGER IERR,LIRN,LICN
      CHARACTER ERR*(*)
!     Local Variables
      REAL*8 EPS,RMIN,RESID
      INTEGER IRNCP,ICNCP,MINIRN,MINICN,IRANK
      LOGICAL ABORT1,ABORT2
      COMMON /MA28FD/ EPS,RMIN,RESID,IRNCP,ICNCP,MINIRN,MINICN,IRANK,
     '  ABORT1,ABORT2


C     OK decomposition
      IF(IERR.EQ.0) THEN
        ERR='Successful decomposition'
      ELSE IF(IERR.EQ.1) THEN
        ERR='Successful decomposition on a structurally singular matrix'
      ELSE IF(IERR.EQ.2) THEN
        ERR='Successful decomposition on a numerically singular matrix'

C     Problems
      ELSE IF(IERR.EQ.-1) THEN
        ERR='Matrix structurally singular'
      ELSE IF(IERR.EQ.-2) THEN
        ERR='Matrix numerically singular'
      ELSE IF(IERR.EQ.-3) THEN
        WRITE(ERR,7000) 'LIRN too small....increase to ',MINIRN,' from '
     '    ,LIRN
      ELSE IF(IERR.EQ.-4) THEN
        ERR='LICN far too small'
      ELSE IF(IERR.EQ.-5) THEN
        WRITE(ERR,7000) 'LICN too small....increase to ',MINICN,' from '
     '    ,LICN
      ELSE IF(IERR.EQ.-6) THEN
        WRITE(ERR,7000) 'LIRN and LICN too small....increase to',MINIRN,
     '    ' and ',MINICN,' from ',LIRN,' and ',LICN
      ELSE IF(IERR.EQ.-7) THEN
        ERR='Error encountered during block triangularization phase'
      ELSE IF(IERR.EQ.-8) THEN
        ERR='LIRN < NZ'
      ELSE IF(IERR.EQ.-9) THEN
        ERR='LICN < NZ'
      ELSE IF(IERR.EQ.-10) THEN
        ERR='NZ <= 0'
      ELSE IF(IERR.EQ.-11) THEN
        ERR='N not in range 1 <= N <= MAXINT'
      ELSE IF(IERR.EQ.-12) THEN
        ERR='Row or column index out-of-range'
      ELSE IF(IERR.EQ.-13) THEN
        ERR='Nonzero was not present during the previous call'
      ELSE IF(IERR.EQ.-14) THEN
        ERR='More than one nonzero in same position in matrix'
      ELSE IF(IERR.EQ.-15) THEN
        ERR='A call to MA28BD follows a call where entries were dropped'
      ELSE IF(IERR.EQ.-16) THEN
        ERR='More than MAXIT iterations required'
      ELSE IF(IERR.EQ.-17) THEN
        ERR='Convergence of iterative refinement too slow'
      ELSE
        WRITE(ERR,7000) 'Unknown Error ',IERR
      ENDIF

      RETURN

 7000 FORMAT(4(A,I10))

 9999 CALL ERRORS('MA28AD_ERR',ERR)
      RETURN 1
      END
