
      SUBROUTINE UMD2FA_ALLOC(N,NZA,NUMF,NINDX,NDARR,WORK_PTR,INDEX_PTR,
     '  ERROR,*)

C#### Subroutine: UMD2FA_ALLOC
C###  Description:
C###    UMD2FA_ALLOC allocates memory for the UMFPACK LU sparse direct
C###    linear solver, UMD2FA().
C###  Written by Stuart Norris 03/06/02

      IMPLICIT NONE
      INCLUDE 'mach00.inc'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'solver.inc'
      INTEGER N,NZA,NUMF,NINDX,NDARR
      POINTER WORK_PTR,INDEX_PTR
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables


      IF(MAXINT/(2*NUMF).LT.NZA) THEN
        WRITE(0,*) 'UMD2FA_ALLOC: Warning, restricting malloc'
        NINDX=MAXINT
        NDARR=(MAXINT/2) + 1
      ELSE
        NINDX=NZA*NUMF*2
        NDARR=NZA*NUMF
      ENDIF
      CALL ALLOCATE_MEMORY(NINDX,1,INTTYPE,INDEX_PTR,NOINIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(NDARR,1, DPTYPE, WORK_PTR,NOINIT,ERROR,*9999)

      RETURN

 9999 CALL ERRORS('UMD2FA_ALLOC',ERROR)
      RETURN 1
      END


      SUBROUTINE UMD2FA_FREE(WORK_PTR,INDEX_PTR,ERROR,*)

C#### Subroutine: UMD2FA_FREE
C###  Description:
C###    UMD2FA_FREE frees memory allocated in UMD2FA_ALLOC.
C###  Written by Stuart Norris 04/06/02

      IMPLICIT NONE
      POINTER WORK_PTR,INDEX_PTR
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables


      IF(WORK_PTR.NE.0) CALL FREE_MEMORY(WORK_PTR,ERROR,*9999)
      IF(INDEX_PTR.NE.0) CALL FREE_MEMORY(INDEX_PTR,ERROR,*9999)

      RETURN

 9999 CALL ERRORS('UMD2FA_FREE',ERROR)
      RETURN 1
      END


      SUBROUTINE UMD2FA_FACTOR(A,N,NZA,ISC_A,ISR_A,WORK,INDEX,NINDX,
     '  NDARR,JOB,KEEP,ICNTL,INFO,CNTL,RINFO,UMFDEF,OUTPUTCODE,ERROR,*)

C#### Subroutine: UMD2FA_FACTOR
C###  Description:
C###    Factorise a sparse system using UMPACK's UMD2FA() LU code
C###  Written by Stuart Norris 08/05/02

      IMPLICIT NONE
!     Parameter List
      INTEGER N,NZA,ISR_A(*),ISC_A(*),INDEX(*),OUTPUTCODE
      INTEGER NINDX,NDARR,JOB,KEEP(20),ICNTL(20),INFO(40)
      REAL*8 A(*),WORK(*),CNTL(10),RINFO(20),UMFDEF(*)
      CHARACTER ERROR*(*)
!     Local Variables


C     Set solver defaults, and copy arrays
      CALL UMD2FA_SETUP(N,NZA,A,WORK,ISC_A,ISR_A,INDEX,ERROR,*9999)
      CALL UMD21I(KEEP,CNTL,ICNTL)
      IF(OUTPUTCODE.GT.2) ICNTL(3)=4
C     ICNTL(6)=1
C     IF(UMF_ITER) JOB=1
      IF(UMFDEF(1).GE.0.0D0)  CNTL(1)=UMFDEF(1)
      IF(UMFDEF(2).GE.0.0D0) ICNTL(7)=NINT(UMFDEF(2))
      IF(UMFDEF(3).GE.0.0D0) ICNTL(4)=NINT(UMFDEF(3))
      IF(UMFDEF(4).GE.0.0D0) ICNTL(5)=NINT(UMFDEF(4))
      IF(UMFDEF(5).GE.0.0D0) THEN
        ICNTL(6)=NINT(UMFDEF(5))
      ELSE
        ICNTL(6)=1
      ENDIF
      IF(UMFDEF(6).GE.0.0D0)  CNTL(2)=UMFDEF(6)
      IF(UMFDEF(7).GE.0.0D0)  KEEP(7)=NINT(UMFDEF(7))
      IF(UMFDEF(8).GE.0.0D0)  KEEP(8)=NINT(UMFDEF(8))

C     Factorise
      CALL UMD2FA(N,NZA,JOB,.FALSE.,NDARR,NINDX,WORK,INDEX,KEEP,CNTL,
     '  ICNTL,INFO,RINFO)
      CALL UMD2FA_ERR(INFO(1),ERROR,*9999)

      RETURN

 9999 CALL ERRORS('UMD2FA_FACTOR',ERROR)
      RETURN 1
      RETURN
      END


      SUBROUTINE UMD2FA_SOLVE(A,B,X,N,NZA,WORK,INDEX,NINDX,NDARR,KEEP,
     '  ICNTL,INFO,CNTL,RINFO,ERROR,*)

C#### Subroutine: UMD2FA_SOLVE
C###  Description:
C###    Solve a system that has been factorised with UMFPACKS's UMD2FA
C###  Written by Stuart Norris 03/06/02

      IMPLICIT NONE
      INCLUDE 'mach00.inc'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'solver.inc'
!     Parameter List
      INTEGER N,NZA,INDEX(*)
      INTEGER NINDX,NDARR,JOB,KEEP(20),ICNTL(20),INFO(40)
      REAL*8 A(*),WORK(*),B(*),X(*),CNTL(10),RINFO(20)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER LWORK
      POINTER WORK_PTR


C     Allocate workspace
      IF(MAXINT/4.LT.N) THEN
        WRITE(ERROR,7000) '>>Cannot alloc 4*',N,' doubles'
        GOTO 9999
      ENDIF
      LWORK=4*N
      CALL ALLOCATE_MEMORY(LWORK,1,DPTYPE,WORK_PTR,NOINIT,ERROR,*9999)

      CALL UMD2SO(N,0,.FALSE.,NDARR,NINDX,WORK,INDEX,KEEP,B,X,
     '  %VAL(WORK_PTR),CNTL,ICNTL,INFO,RINFO)
      IF(INFO(1).LT.0) THEN
        WRITE(ERROR,'(A,I10)') '>>UMD2SO returns ',INFO(1)
        GOTO 9999
      ENDIF

C     Free memory
      IF(WORK_PTR.NE.0) CALL FREE_MEMORY(WORK_PTR,ERROR,*9999)

      RETURN

 7000 FORMAT(A,I20,A)

 9999 CALL ERRORS('UMD2FA_SOLVE',ERROR)
      RETURN 1
      END


      SUBROUTINE UMD2FA_SETUP(N,NZA,A,WORK,ISC_A,ISR_A,INDEX,ERROR,*)

C#### Subroutine: UMD2FA_SETUP
C###  Description:
C###    Set up data arrays for the Umfpack sparse LU solver.
C###  Written by Stuart Norris 09/05/02

      IMPLICIT NONE
!     Parameter List
      INTEGER N,NZA,ISC_A(*),ISR_A(*),INDEX(*)
      REAL*8 A(*),WORK(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,IJ


      CALL DCOPY(NZA,A,1,WORK,1)
      CALL ICOPY(NZA,ISC_A,1,INDEX(NZA+1),1)

      DO I=1,N
        DO IJ=ISR_A(I),ISR_A(I+1)-1
          INDEX(IJ)=I
        ENDDO
      ENDDO

      RETURN

 9999 CALL ERRORS('UMD2FA_SETUP',ERROR)
      RETURN 1
      END


      SUBROUTINE UMD2FA_ERR(IERR,ERROR,*)

C#### Subroutine: UMD2FA_ERR
C###  Description:
C###    Print errors from UMFPACKS's UMD2FA solver
C###  Written by Stuart Norris 03/06/02

      IMPLICIT NONE
!     Parameter List
      INTEGER IERR
      CHARACTER ERROR*(*)
!     Local Variables


C     Errors
      IF(IERR.LT.0) THEN
        IF(IERR.EQ.-1) THEN
          ERROR='>>Matrix order not positive in UMFPACK'
          GOTO 9999
        ELSE IF(IERR.EQ.-2) THEN
          ERROR='>>Number of entries not positive in UMFPACK'
          GOTO 9999
        ELSE IF(IERR.EQ.-3) THEN
          ERROR='>>LINDEX too small in UMFPACK'
          GOTO 9999
        ELSEIF(IERR.EQ.-4) THEN
          ERROR='>>LVALUE too small in UMFPACK'
          GOTO 9999
        ELSE IF(IERR.EQ.-5) THEN
          ERROR='>>LINDEX and LVALUE too small in UMFPACK'
          GOTO 9999
        ELSE
          WRITE(ERROR,'(''>>UMFPACK UMD2FA returns '',I10)') IERR
          GOTO 9999
        ENDIF

C     OK decomposition
      ELSE IF(IERR.EQ.0) THEN
        ERROR='Successful decomposition'

C     Positive error message
      ELSE IF(IERR.GT.1) THEN
        WRITE(ERROR,'(''>>UMFPACK UMD2FA returns '',I10)') IERR

      ENDIF

      RETURN

 9999 CALL ERRORS('UMD2FA_ERR',ERROR)
      RETURN 1
      END

