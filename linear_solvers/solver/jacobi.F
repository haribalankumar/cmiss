
      SUBROUTINE POINT_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,ERROR,*)

C#### Subroutine: POINT_FACTOR
C###  Description:
C###    POINT_FACTOR factorises a system for use with the point Jacobi 
C###    preconditioner. Unlike Jacobi iteration it is kludged to allow
C###    systems with zero elements on the diagonal.
C###
C###    I don't like this code, but it is given for backwards compatability
C###    with what is in CMISS.
C###  Written by Stuart Norris 09/80/02

      IMPLICIT NONE
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*),LDA,N,SPARSE_A
      REAL*8 A(*),D(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J


C     Dense matrix.
      IF(SPARSE_A.EQ.0) THEN

C$OMP   PARALLEL DO PRIVATE(I)
        DO I=1,N
          IF(A(I+N*(I-1)).EQ.0.0D0) THEN
            D(I)=1.0D0 ! Kludge. We should die here
          ELSE
            D(I)=1.0D0/A(I+LDA*(I-1))
          ENDIF
        ENDDO
C$OMP   END PARALLEL DO

C     Sparse matricies: CMISS compressed row
      ELSE IF(SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4) THEN

C$OMP   PARALLEL DO PRIVATE(I,J)
        DO I=1,N
          DO J=ISR_A(I),ISR_A(I+1)-1
            IF(ISC_A(J).EQ.I) THEN
              IF(A(J).EQ.0.0D0) THEN
                D(I)=1.0D0 ! Kludge. We should die here
              ELSE
                D(I)=1.0D0/A(J)
              ENDIF
              GOTO 100
            ENDIF
          ENDDO
          D(I)=1.0D0 ! Kludge. We should die here
 100      CONTINUE
        ENDDO
C$OMP   END PARALLEL DO

      ELSE
        ERROR='>>Sparsity type not implemented'
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('POINT_FACTOR',ERROR)
      RETURN 1
      END


      SUBROUTINE ROWPIV_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,IPIVOT_A,
     '  ERROR,*)

C#### Subroutine: ROWPIV_FACTOR
C###  Description:
C###    ROWPIV_FACTOR factorises a system for use with the row pivoted
C###    point Jacobi preconditioner as present in CMISS.
C###
C###    This code matches "Kludge for Kludge" the code in CMISS. It is
C###    beyond the authors skill to see how the original CMISS code was
C###    expected to do anything.
C###  Written by Stuart Norris 09/80/02

      IMPLICIT NONE
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*),LDA,N,SPARSE_A,IPIVOT_A(*)
      REAL*8 A(*),D(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J,IJ
      REAL*8 MAXVAL
      LOGICAL FOUND


C     Dense matrix.
      IF(SPARSE_A.EQ.0) THEN

C$OMP   PARALLEL DO PRIVATE(I)
        DO I=1,N
          IF(A(I+N*(I-1)).EQ.0.0D0) THEN
            D(I)=1.0D0 ! Kludge
          ELSE
            D(I)=1.0D0/A(I+LDA*(I-1))
          ENDIF
        ENDDO
C$OMP   END PARALLEL DO

C     Sparse matricies: CMISS compressed row
      ELSE IF(SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4) THEN

        DO I=1,N
          IPIVOT_A(I)=I
C         Find the diagonal
          DO J=ISR_A(I),ISR_A(I+1)-1
            IF(ISC_A(J).EQ.I) THEN
              IF(A(J).EQ.0.0D0) THEN
                D(I)=1.0D0
                FOUND=.FALSE.
              ELSE
                D(I)=1.0D0/A(J)
                FOUND=.TRUE.
              ENDIF
              GOTO 100
            ENDIF
          ENDDO
          D(I)=1.0D0
          FOUND=.FALSE.
 100      CONTINUE

C         If we havn't found a non-zero diagonal element, then search for 
C         a suitable pivot partner
          IF(.NOT.FOUND) THEN
            DO J=1,I-1
              IF(IPIVOT_A(J).EQ.J) THEN
                MAXVAL=0.0D0
                DO IJ=ISR_A(I),ISR_A(I+1)-1
                  IF(ABS(A(IJ)).GT.ABS(MAXVAL)) THEN
                    MAXVAL=A(IJ)
                  ENDIF
                ENDDO
                IF(ABS(MAXVAL).GT.0.0D0) THEN
                  IPIVOT_A(I)=J
                  IPIVOT_A(J)=I
                  D(I)=1.0D0/D(J)
                  D(J)=1.0D0/MAXVAL
                  GOTO 200
                ENDIF
              ENDIF
            ENDDO

C           DO IJ=ISR_A(I),ISR_A(I+1)-1
C             K=ISC_A(IJ)
C             IF(K.NE.I.AND.IPIVOT_A(K).EQ.K .AND. A(IJ).NE.0.0D0) THEN
C               DO KJ=ISR_A(K),ISR_A(K+1)-1
C                 IF(ISC_A(KJ).EQ.I .AND. A(KJ).NE.0.0D0) THEN
C                   IPIVOT_A(K)=I
C                   IPIVOT_A(I)=K
C                   D(K)=1.0D0/A(IJ)
C                   D(I)=1.0D0/A(KJ)
C                   GOTO 200
C                 ENDIF
C               ENDDO
C             ENDIF
C           ENDDO

          ENDIF
 200      CONTINUE

        ENDDO

      ELSE
        ERROR='>>Sparsity type not implemented'
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('ROWPIV_FACTOR',ERROR)
      RETURN 1
      END


      SUBROUTINE ROWSCALE_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,ERROR,*)

C#### Subroutine: ROWSCALE_FACTOR
C###  Description:
C###    ROWSCALE_FACTOR factorises a system for use with the row scale
C###    point Jacobi preconditioner -- the system is scaled according to 
C###    magnitude of the row.
C###  Written by Stuart Norris 09/80/02

      IMPLICIT NONE
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*),LDA,N,SPARSE_A
      REAL*8 A(*),D(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,IJ
      REAL*8 SUM
!     External Functions
      REAL*8 DASUM
      EXTERNAL DASUM


C     Dense matrix.
      IF(SPARSE_A.EQ.0) THEN

C$OMP   PARALLEL DO PRIVATE(I,SUM), FIRSTPRIVATE(N,LDA)
        DO I=1,N
          SUM=DASUM(N,A(I),LDA)

          IF(SUM.EQ.0.0D0) THEN
            D(I)=1.0D0
          ELSE
C           D(I)=DBLE(N)/SUM
            D(I)=2.0D0/SUM
          ENDIF
        ENDDO
C$OMP   END PARALLEL DO

C     Sparse matricies: CMISS compressed row
      ELSE IF(SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4) THEN

C$OMP   PARALLEL DO PRIVATE(I,SUM)
        DO I=1,N
          SUM=0.0D0
          DO IJ=ISR_A(I),ISR_A(I+1)-1
            SUM=SUM+ABS(A(IJ))
          ENDDO

          IF(SUM.EQ.0.0D0) THEN
            D(I)=1.0D0
          ELSE
C           D(I)=DBLE(ISR_A(I+1)-ISR_A(I))/SUM
            D(I)=2.0/SUM
          ENDIF
        ENDDO
C$OMP   END PARALLEL DO

      ELSE
        ERROR='>>Sparsity type not implemented'
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('ROWSCALE_FACTOR',ERROR)
      RETURN 1
      END


      SUBROUTINE JACOBI_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,ERROR,*)

C#### Subroutine: JACOBI_FACTOR
C###  Description:
C###    JACOBI_FACTOR factorises a system for use in Jacobian iteration,
C###    either the preconditioner or the stand alone solver. All it entails 
C###    is the identification of the diagonal, and the inversion of this
C###    value since divides are expensive (or used to be. Must check to
C###    see if they still are).
C###    May also be used for SOR, SSOR and Point Jacobi.
C###  Written by Stuart Norris 26/02/02

      IMPLICIT NONE
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*),LDA,N,SPARSE_A
      REAL*8 A(*),D(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J,ID0,ID1


C     Dense matrix.
      IF(SPARSE_A.EQ.0) THEN

        ID0=N+1
C$OMP   PARALLEL DO PRIVATE(I), REDUCTION(MIN:ID0)
        DO I=1,N
C         Rather than just check for cases where A(i,i)=0, we do this
C         minimum reduction, so that OpenMP will not complain about
C         jumping out of a parallel loop
          IF(A(I+N*(I-1)).EQ.0.0D0) THEN
            ID0=MIN(ID0,I)
C           D(I)=1.0D0 ! Kludge
          ELSE
            D(I)=1.0D0/A(I+LDA*(I-1))
          ENDIF
        ENDDO
C$OMP   END PARALLEL DO
        IF(ID0.LE.N) THEN
          WRITE(ERROR,'(A,2(1X,I9))') '>>Zero on diagonal: ',ID0,N
          GOTO 9999
        ENDIF

C     Sparse matricies: CMISS compressed row
      ELSE IF(SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4) THEN

        ID0=N+1
        ID1=N+1
C$OMP   PARALLEL DO PRIVATE(I,J), REDUCTION(MIN:ID0,ID1)
        DO I=1,N
          DO J=ISR_A(I),ISR_A(I+1)-1
            IF(ISC_A(J).EQ.I) THEN
C             Rather than just check for cases where A(i,i)=0, we do this
C             minimum reduction, so that OpenMP will not complain about 
C             jumping out of a parallel loop
              IF(A(J).EQ.0.0D0) THEN
                ID0=MIN(ID0,I)
C               D(I)=1.0D0 ! Kludge
              ELSE
                D(I)=1.0D0/A(J)
              ENDIF
              GOTO 100
            ENDIF
          ENDDO
          ID1=MIN(ID1,I)
 100      CONTINUE
        ENDDO
C$OMP   END PARALLEL DO
        IF(ID0.LE.N) THEN
          WRITE(ERROR,'(A,2(1X,I9))') '>>Zero on diagonal: ',ID0,N
          GOTO 9999
        ENDIF
        IF(ID1.LE.N) THEN
          WRITE(ERROR,'(A,2(1X,I9))') '>>No diagonal element: ',ID1,N
          GOTO 9999
        ENDIF

      ELSE
        ERROR='>>Sparsity type not implemented'
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('JACOBI_FACTOR',ERROR)
      RETURN 1
      END


      SUBROUTINE PRECON_JACOBI(A,LDA,N,X,B,D,DX,SPARSE_A,ISC_A,ISR_A,
     '  OMEGA,MAXIT,ERROR,*)

C#### Subroutine: PRECON_JACOBI
C###  Description:
C###    PRECON_JACOBI applies MAXIT iterations of the Jacobi preconditioner
C###    to the system 
C###           A x = b
C###  Written by Stuart Norris 26/02/02

      IMPLICIT NONE
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*),LDA,N,SPARSE_A,MAXIT
      REAL*8 A(*),X(*),B(*),D(*),DX(*),OMEGA
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,IJ,ITER
      REAL*8 SUM
!     Functions
      REAL*8 DDOT
      EXTERNAL DDOT,PDAXPY


C     Apply preconditioning
C$OMP PARALLEL DO PRIVATE(I)
      DO I=1,N
        X(I)=OMEGA*B(I)*D(I)
      ENDDO
C$OMP END PARALLEL DO

      DO ITER=2,MAXIT

C       Dense arrays
        IF(SPARSE_A.EQ.0) THEN
C$OMP     PARALLEL DO PRIVATE(I,SUM), FIRSTPRIVATE(N,LDA)
          DO I=1,N
            DX(I)=D(I)*(B(I)-DDOT(N,A(I),LDA,X,1))
          ENDDO
C$OMP     END PARALLEL DO

C       Sparse arrays
        ELSE IF(SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4) THEN
C$OMP     PARALLEL DO PRIVATE(I,IJ,SUM)
          DO I=1,N
            SUM=0.0D0
            DO IJ=ISR_A(I),ISR_A(I+1)-1
              SUM=SUM+A(IJ)*X(ISC_A(IJ))
            ENDDO
            DX(I)=D(I)*(B(I)-SUM)
          ENDDO
C$OMP     END PARALLEL DO

        ELSE
          ERROR='>>Sparsity type not implemented'
          GOTO 9999
        ENDIF

C       Update solution vector
        CALL PDAXPY(N,OMEGA,DX,1,X,1)

      ENDDO

      RETURN

 9999 CALL ERRORS('PRECON_JACOBI',ERROR)
      RETURN 1
      END


      SUBROUTINE PRECON_SSOR(A,LDA,N,X,B,D,SPARSE_A,ISC_A,ISR_A,OMEGA,
     '  MAXIT,ERROR,*)

C#### Subroutine: PRECON_SSOR
C###  Description:
C###    PRECON_SSOR applies MAXIT iterations of the SSOR preconditioner
C###    to the system 
C###           A x = b
C###    The OpenMP directives have been deleted to ensure that the serial
C###    and parallel libraries give the same results.
C###  Written by Stuart Norris 26/02/02

      IMPLICIT NONE
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*),LDA,N,SPARSE_A,MAXIT
      REAL*8 A(*),X(*),B(*),D(*),OMEGA
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,IJ,ITER
      REAL*8 DX,SUM
!     Functions
      REAL*8 PDDOT
      EXTERNAL PDDOT,DAXPY


C     Apply preconditioning
      DO ITER=1,(MAXIT+1)/2

C       Dense arrays
        IF(SPARSE_A.EQ.0) THEN

          DO I=1,N
            DX=D(I)*(B(I)-PDDOT(N,A(I),LDA,X,1))

            X(I)=X(I)+OMEGA*DX
          ENDDO

          DO I=N,1,-1
            DX=D(I)*(B(I)-PDDOT(N,A(I),LDA,X,1))

            X(I)=X(I)+OMEGA*DX
          ENDDO

C       Sparse arrays
        ELSE IF(SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4) THEN

          DO I=1,N
            SUM=0.0D0
            DO IJ=ISR_A(I),ISR_A(I+1)-1
              SUM=SUM+A(IJ)*X(ISC_A(IJ))
            ENDDO
            DX=D(I)*(B(I)-SUM)

            X(I)=X(I)+OMEGA*DX
          ENDDO

          DO I=N,1,-1
            SUM=0.0D0
            DO IJ=ISR_A(I),ISR_A(I+1)-1
              SUM=SUM+A(IJ)*X(ISC_A(IJ))
            ENDDO
            DX=D(I)*(B(I)-SUM)

            X(I)=X(I)+OMEGA*DX
          ENDDO

        ELSE
          ERROR='>>Sparsity type not implemented'
          GOTO 9999
        ENDIF

      ENDDO

      RETURN

 9999 CALL ERRORS('PRECON_SSOR',ERROR)
      RETURN 1
      END


      SUBROUTINE JACOBI(A,LDA,N,X,B,D,SPARSE_A,ISC_A,ISR_A,ANORM,RESID,
     '  OMEGA,ITER,FACTORISE,OUTPUTCODE,ERROR,*)

C#### Subroutine: JACOBI
C###  Description:
C###    JACOBI solves a system of linear equations using Jacobi's method,
C###    possibly the simplest iterative scheme. For each iteration x is 
C###    updated by
C###          x <= x + omega dx
C###    where dx is calculated from 
C###          dx <= (b - Ax + Dx)/D - x
C###    D being the diagonal of A.
C###  Written by Stuart Norris 26/02/02

      IMPLICIT NONE
      INCLUDE 'mach00.inc'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'solver.inc'
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*),LDA,N,SPARSE_A,ITER,OUTPUTCODE
      REAL*8 A(*),X(*),B(*),D(*),ANORM,RESID,OMEGA
      CHARACTER ERROR*(*)
      LOGICAL FACTORISE
!     Local Variables
      POINTER DX_PTR
      DATA DX_PTR / 0 /


C     Check input
      IF(ITER.LE.0) THEN
        ERROR='>>Maxit less than 0'
        GOTO 9999
      ENDIF
      IF(RESID.LT.0.0D0) THEN
        ERROR='>>Tolerance less than 0.0'
        GOTO 9999
      ENDIF

C     Factorise the system
      IF(FACTORISE) THEN
        CALL JACOBI_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,ERROR,*9999)
      ENDIF

C     Solve the system
      CALL ALLOCATE_MEMORY(N,1,DPTYPE,DX_PTR,NOINIT,ERROR,*9999)

      CALL JACOBI_SUB(A,LDA,N,X,B,D,%VAL(DX_PTR),SPARSE_A,ISC_A,ISR_A,
     '  ANORM,RESID,OMEGA,ITER,OUTPUTCODE,ERROR,*9999)

      IF(DX_PTR.NE.0) CALL FREE_MEMORY(DX_PTR,ERROR,*9999)
      RETURN

 9999 CALL ERRORS('JACOBI',ERROR)
      RETURN 1
      END


      SUBROUTINE JACOBI_SUB(A,LDA,N,X,B,D,DX,SPARSE_A,ISC_A,ISR_A,ANORM,
     '  RESID,OMEGA,ITER,OUTPUTCODE,ERROR,*)

C#### Subroutine: JACOBI_SUB
C###  Description:
C###    JACOBI solves a system of linear equations using Jacobi's method,
C###    possibly the simplest iterative scheme. For each iteration x is 
C###    updated by
C###          x <= x + omega dx
C###    where dx is calculated from 
C###          dx <= (b - Ax + Dx)/D - x
C###    D being the diagonal of A.
C###  Written by Stuart Norris 26/02/02

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*),LDA,N,SPARSE_A,ITER,OUTPUTCODE
      REAL*8 A(*),X(*),B(*),D(*),DX(*),ANORM,RESID,OMEGA
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,IJ,MAXIT
      REAL*8 TOL,BNORM,SCALED_TOL,SUM,RSUM
!     Functions
      REAL*8 DDOT,PDNRM2
      EXTERNAL DDOT,PDNRM2,PDAXPY


      TOL=RESID
      MAXIT=ITER
      BNORM=PDNRM2(N,B,1)
      SCALED_TOL=TOL*BNORM
C$OMP PARALLEL DO PRIVATE(I)
      DO I=1,N
        DX(I)=0.0D0
      ENDDO
C$OMP END PARALLEL DO

      IF(OUTPUTCODE.GE.2) THEN
        CALL WRITE_JACOBI_SETTINGS(N,SPARSE_A,ISR_A,ANORM,BNORM,
     '    TOL,OMEGA,MAXIT,'Jacobi Solver',ERROR,*9999)
      ENDIF

C     Solve the system
      DO ITER=1,MAXIT

C       Calculate DX, the solution update.

C       Dense arrays
        IF(SPARSE_A.EQ.0) THEN
          RSUM=0.0D0
C$OMP     PARALLEL DO PRIVATE(I,SUM), FIRSTPRIVATE(N,LDA),
C$OMP'      REDUCTION(+:RSUM)
          DO I=1,N
            SUM=B(I)-DDOT(N,A(I),LDA,X,1)

            RSUM=RSUM+SUM**2
            DX(I)=D(I)*SUM
          ENDDO
C$OMP     END PARALLEL DO

C       Sparse arrays
        ELSE IF(SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4) THEN

          RSUM=0.0D0
C$OMP     PARALLEL DO PRIVATE(I,IJ,SUM), REDUCTION(+:RSUM)
          DO I=1,N
            SUM=B(I)
            DO IJ=ISR_A(I),ISR_A(I+1)-1
              SUM=SUM-A(IJ)*X(ISC_A(IJ))
            ENDDO

            RSUM=RSUM+SUM**2
            DX(I)=D(I)*SUM
          ENDDO
C$OMP     END PARALLEL DO

        ELSE
          ERROR='>>Sparsity type not implemented'
          GOTO 9999
        ENDIF

C       Update the solution
        RESID=SQRT(RSUM)
        CALL PDAXPY(N,OMEGA,DX,1,X,1)

C       On the first pass, estimate the scaling factor for the residual.
        IF(ITER.EQ.1) SCALED_TOL=TOL*(ANORM*PDNRM2(N,X,1) + BNORM)

C       Print out solution, and check for convergance
        IF(OUTPUTCODE.GE.2) THEN
          WRITE(OP_STRING,'(1P,I5,1x,G12.6)') ITER,RESID
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        IF(RESID.LE.SCALED_TOL) GOTO 100

      ENDDO
      ITER=MAXIT

 100  CONTINUE

      RETURN

 9999 CALL ERRORS('JACOBI_SUB',ERROR)
      RETURN 1
      END


      SUBROUTINE SOR(A,LDA,N,X,B,D,SPARSE_A,ISC_A,ISR_A,ANORM,RESID,
     '  OMEGA,ITER,FACTORISE,OUTPUTCODE,ERROR,*)

C#### Subroutine: SOR
C###  Description:
C###    SOR solves a system of linear equations using Sucessive Over
C###    relaxation. The OpenMP directives have been deleted to ensure
C###    that the serial and parallel libraries give the same results.
C###  Written by Stuart Norris 26/02/02

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*),LDA,N,SPARSE_A,ITER,OUTPUTCODE
      REAL*8 A(*),X(*),B(*),D(*),ANORM,RESID,OMEGA
      CHARACTER ERROR*(*)
      LOGICAL FACTORISE
!     Local Variables
      INTEGER I,IJ,MAXIT
      REAL*8 TOL,DX,RSUM,SUM,BNORM,SCALED_TOL
!     Functions
      REAL*8 PDDOT,PDNRM2
      EXTERNAL PDDOT,PDNRM2


      TOL=RESID
      MAXIT=ITER
      BNORM=PDNRM2(N,B,1)
      SCALED_TOL=TOL*BNORM
      IF(OUTPUTCODE.GE.2) THEN
        CALL WRITE_JACOBI_SETTINGS(N,SPARSE_A,ISR_A,ANORM,BNORM,
     '    TOL,OMEGA,MAXIT,'SOR Solver',ERROR,*9999)
      ENDIF

C     Check input
      IF(MAXIT.LE.0) THEN
        ERROR='>>Maxit less than 0'
        GOTO 9999
      ENDIF
      IF(TOL.LT.0.0d0) THEN
        ERROR='>>Tolerance less than 0.0'
        GOTO 9999
      ENDIF

C     Factorise the system
      IF(FACTORISE) THEN
        CALL JACOBI_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,ERROR,*9999)
      ENDIF

C     Solve the system
      DO ITER=1,MAXIT

C       Calculate DX, and update the solution

C       Dense arrays
        IF(SPARSE_A.EQ.0) THEN

          RSUM=0.0D0
          DO I=1,N
            SUM=B(I)-PDDOT(N,A(I),LDA,X,1)
            DX=D(I)*SUM

            RSUM=RSUM+SUM**2
            X(I)=X(I)+OMEGA*DX
          ENDDO

C       Sparse arrays
        ELSE IF(SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4) THEN

          RSUM=0.0D0
          DO I=1,N
            SUM=B(I)
            DO IJ=ISR_A(I),ISR_A(I+1)-1
              SUM=SUM-A(IJ)*X(ISC_A(IJ))
            ENDDO
            DX=D(I)*SUM

            RSUM=RSUM+SUM**2
            X(I)=X(I)+OMEGA*DX
          ENDDO

        ELSE
          ERROR='>>Sparsity type not implemented'
          GOTO 9999
        ENDIF

        RESID=SQRT(RSUM)
        IF(ITER.EQ.1) SCALED_TOL=TOL*(ANORM*PDNRM2(N,X,1) + BNORM)

C       Print out solution, and check for convergance
        IF(OUTPUTCODE.GE.2) THEN
          WRITE(OP_STRING,'(1P,I5,1x,G12.6)') ITER,RESID
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        IF(RESID.LE.SCALED_TOL) GOTO 100

      ENDDO
      ITER=MAXIT

 100  CONTINUE

      RETURN

 9999 CALL ERRORS('SOR',ERROR)
      RETURN 1
      END


      SUBROUTINE WRITE_JACOBI_SETTINGS(N,SPARSE_A,ISR_A,ANORM,BNORM,
     '  TOL,OMEGA,MAXIT,NAME,ERROR,*)

C#### Subroutine: WRITE_JACOBI_SETTINGS
C###  Description:
C###    Write out the solver settings for the Jacobi/SOR solvers.
C###  Written by Stuart Norris 07/08/02

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER N,SPARSE_A,ISR_A(*),MAXIT
      REAL*8 ANORM,BNORM,TOL,OMEGA
      CHARACTER NAME*(*),ERROR*(*)
!     Local variables
      INTEGER NZA
!     Functions
      INTEGER LEN_TRIM
      EXTERNAL LEN_TRIM


      IF(SPARSE_A.EQ.0) THEN
        NZA=N**2
      ELSE
        NZA=ISR_A(N+1)-1
      ENDIF

      WRITE(OP_STRING,'(1X,A)') NAME(:LEN_TRIM(NAME))
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'()')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Relaxation Factor  :'',1P,E11.4)') OMEGA
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Solution Tolerence :'',1P,E11.4)') TOL
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Maximum Iterations :'',I11)') MAXIT
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'()')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      CALL WRITE_MATRIX_DIMS(SPARSE_A,N,NZA,ANORM,BNORM,ERROR,*9999)

      RETURN

 9999 CALL ERRORS('WRITE_JACOBI_SETTINGS',ERROR)
      RETURN 1
      END
