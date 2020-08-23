
      SUBROUTINE INCOMP_CHOL0_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,
     '  ERROR,*)

C#### Subroutine: INCOMP_CHOL0_FACTOR
C###  Description:
C###    INCOMP_CHOL0_FACTOR factorises a system of equations for solution
C###    with the Incomplete Cholesky(0) iterative solver or preconditioner.
C###  Written by Stuart Norris 26/02/02

      IMPLICIT NONE
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*),LDA,N,SPARSE_A
      REAL*8 A(*),D(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J
      REAL*8 DSUM,RDUM


      DO I=1,N
        D(I)=0.0D0
      ENDDO

C     Dense matrix.
      IF(SPARSE_A.EQ.0) THEN

        DO I=1,N
          IF(A(I+LDA*(I-1)).EQ.0.0D0) THEN
            ERROR='>>Zero on diagonal'
            GOTO 9999
          ENDIF

          RDUM=A(I+LDA*(I-1))
          DSUM=0.0D0
C$OMP     PARALLEL DO PRIVATE(J), REDUCTION(+:DSUM)
          DO J=1,I-1
            DSUM=DSUM+(A(I+LDA*(J-1))*D(J))**2
          ENDDO
C$OMP     END PARALLEL DO

          IF(RDUM.LE.DSUM) THEN
            IF(RDUM.EQ.DSUM) THEN
              ERROR='>>Zero decomposition element'
            ELSE
              ERROR='>>System is not positive definate'
            ENDIF
            GOTO 9999
          ENDIF
          D(I)=1.0D0/SQRT(RDUM-DSUM)
        ENDDO

C     Sparse matricies: CMISS compressed row
      ELSE IF(SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4) THEN

        DO I=1,N
          RDUM=0.0D0
          DSUM=0.0D0
          DO J=ISR_A(I),ISR_A(I+1)-1
            IF(ISC_A(J).GE.I) THEN
              IF(ISC_A(J).EQ.I) THEN
                IF(A(J).EQ.0.0D0) THEN
                  ERROR='>>Zero on diagonal'
                  GOTO 9999
                ENDIF
                RDUM=A(J)
              ENDIF
              GOTO 100
            ENDIF
            DSUM=DSUM+(A(J)*D(ISC_A(J)))**2
          ENDDO
 100      CONTINUE

          IF(RDUM.LE.DSUM) THEN
            IF(RDUM.EQ.DSUM) THEN
              ERROR='>>Zero decomposition element'
            ELSE
              ERROR='>>System is not positive definate'
            ENDIF
            GOTO 9999
          ENDIF
          D(I)=1.0D0/SQRT(RDUM-DSUM)
        ENDDO

      ELSE
        ERROR='>>Sparsity type not implemented'
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('INCOMP_CHOL0_FACTOR',ERROR)
      RETURN 1
      END


      SUBROUTINE PRECON_INCOMP_CHOL0(A,LDA,N,X,B,D,R,SPARSE_A,ISC_A,
     '  ISR_A,MAXIT,ERROR,*)

C#### Subroutine: PRECON_INCOMP_CHOL0
C###  Description:
C###    PRECON_INCOMP_CHOL0 applies MAXIT iterations of the Incomplete 
C###    Cholskey(0) solver to precondition the system
C###           A x = b
C###  Written by Stuart Norris 27/02/02

      IMPLICIT NONE
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*),LDA,N,SPARSE_A,MAXIT
      REAL*8 A(*),X(*),B(*),D(*),R(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J,ITER
      REAL*8 DSUM,RTMP
!     Functions
      REAL*8 PDDOT
      EXTERNAL PDDOT


C     Apply preconditioning
      DO I=1,N
        R(I)=0.0D0
      ENDDO

      DO ITER=1,MAXIT

C       Dense arrays
        IF(SPARSE_A.EQ.0) THEN

C         Calculate the residual, and take the forward sweep
          DO I=1,N
            RTMP=B(I)-PDDOT(N,A(I),LDA,X,1)

            DSUM=0.0D0
C$OMP       PARALLEL DO PRIVATE(J), REDUCTION(+:DSUM)
            DO J=1,I-1
              DSUM=DSUM+A(I+LDA*(J-1))*R(J)*D(J)
            ENDDO
C$OMP       END PARALLEL DO
            R(I)=D(I)*(RTMP-DSUM)
          ENDDO

C         Backwards sweep, and update the solution
          DO I=N,1,-1
            DSUM=PDDOT(N-I,A(I+LDA*I),LDA,R(I+1),1)
            R(I)=D(I)*( R(I) - D(I)*DSUM )

            X(I)=X(I)+R(I)
          ENDDO

C       Sparse arrays
        ELSE IF(SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4) THEN

C         Calculate the residual, and take the forward sweep
          DO I=1,N

            RTMP=B(I)
            DO J=ISR_A(I),ISR_A(I+1)-1
              RTMP=RTMP-A(J)*X(ISC_A(J))
            ENDDO

            DSUM=0.0D0
            DO J=ISR_A(I),ISR_A(I+1)-1
              IF(ISC_A(J).GE.I) GOTO 100
              DSUM=DSUM+A(J)*R(ISC_A(J))*D(ISC_A(J))
            ENDDO
 100        CONTINUE
            R(I)=D(I)*(RTMP-DSUM)
          ENDDO

C         Backwards sweep, and update the solution
          DO I=N,1,-1

            DSUM=0.0D0
            DO J=ISR_A(I+1)-1,ISR_A(I),-1
              IF(ISC_A(J).LE.I) GOTO 200
              DSUM=DSUM+A(J)*R(ISC_A(J))
            ENDDO
 200        CONTINUE
            R(I)=D(I)*( R(I) - D(I)*DSUM )

            X(I)=X(I)+R(I)
          ENDDO

        ELSE
          ERROR='>>Sparsity type not implemented'
          GOTO 9999
        ENDIF

      ENDDO

      RETURN

 9999 CALL ERRORS('PRECON_INCOMP_CHOL0',ERROR)
      RETURN 1
      END


      SUBROUTINE INCOMP_CHOL0(A,LDA,N,X,B,D,SPARSE_A,ISC_A,ISR_A,ANORM,
     '  RESID,ITER,FACTORISE,OUTPUTCODE,ERROR,*)

C#### Subroutine: INCOMP_CHOL0
C###  Description:
C###    INCOMP_CHOL0 solves a system of linear equations using Incomplete
C###    Cholskey(0) decomposition. The A matrix is decomposed using no infill 
C###    (or pivoting), and then the a sequence of forward and backward
C###    relaxation passes are made.
C###
C###    INCOMP_CHOL0() is a wrapper around the INCOMP_CHOL0_SUB() routine,
C###    to allow memory allocation etc.
C###  Written by Stuart Norris 26/02/02

      IMPLICIT NONE
      INCLUDE 'mach00.inc'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'solver.inc'
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*),LDA,N,SPARSE_A,ITER,OUTPUTCODE
      REAL*8 A(*),X(*),B(*),D(*),ANORM,RESID
      CHARACTER ERROR*(*)
      LOGICAL FACTORISE
!     Local Parameters
      INTEGER QUIET
      PARAMETER(QUIET=1)
!     Local Variables
      INTEGER NZA
      POINTER R_PTR
      DATA R_PTR / 0 /


      IF(OUTPUTCODE.GE.2) THEN
        WRITE(OP_STRING,'('' Incomplete Cholesky(0) Factorisation'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF
      NZA=1
      IF(SPARSE_A.NE.0) NZA=ISR_A(N+1)-1

C     Check input
      IF(ITER.LE.0) THEN
        ERROR='>>Maxit less than 0'
        GOTO 9999
      ENDIF
      IF(RESID.LT.0.0d0) THEN
        ERROR='>>Tolerance less than 0.0'
        GOTO 9999
      ENDIF
      IF(OUTPUTCODE.GE.3) THEN
        CALL CHECK_POS_DEF(A,LDA,N,SPARSE_A,ISC_A,ISR_A,QUIET,ERROR,
     '    *9999)
        CALL CHECK_SYMM(A,LDA,N,NZA,SPARSE_A,ISC_A,ISR_A,QUIET,ERROR,
     '    *9999)
      ENDIF

C     Factorise the system
      IF(FACTORISE) THEN
        CALL INCOMP_CHOL0_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,ERROR,
     '    *9999)
      ENDIF

C     Solve the system
      CALL ALLOCATE_MEMORY(N,1,DPTYPE,R_PTR,NOINIT,ERROR,*9999)

      CALL INCOMP_CHOL0_SUB(A,LDA,N,X,B,D,%VAL(R_PTR),SPARSE_A,ISC_A,
     '  ISR_A,ANORM,RESID,ITER,OUTPUTCODE,ERROR,*9999)

      IF(R_PTR.NE.0) CALL FREE_MEMORY(R_PTR,ERROR,*9999)

C     Bye!
      RETURN

 9999 CALL ERRORS('INCOMP_CHOL0',ERROR)
      RETURN 1
      END


      SUBROUTINE INCOMP_CHOL0_SUB(A,LDA,N,X,B,D,R,SPARSE_A,ISC_A,ISR_A,
     '  ANORM,RESID,ITER,OUTPUTCODE,ERROR,*)

C#### Subroutine: INCOMP_CHOL0_SUB
C###  Description:
C###    INCOMP_CHOL0 solves a system of linear equations using Incomplete
C###    Cholskey(0) decomposition. The A matrix is decomposed using no infill 
C###    (or pivoting), and then the a sequence of forward and backward
C###    relaxation passes are made.
C###
C###    INCOMP_CHOL0_SUB() is not meant to be called directly. Instead
C###    you are meant to call INCOMP_CHOL0() which handes memory issues.
C###  Written by Stuart Norris 26/02/02

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*),LDA,N,SPARSE_A,ITER,OUTPUTCODE
      REAL*8 A(*),X(*),B(*),D(*),R(*),ANORM,RESID
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J,MAXIT
      REAL*8 TOL,DSUM,RTMP,BNORM,SCALED_TOL
!     Functions
      REAL*8 PDDOT,PDNRM2
      EXTERNAL PDDOT,PDNRM2


      TOL=RESID
      MAXIT=ITER
      BNORM=PDNRM2(N,B,1)
      SCALED_TOL=TOL*BNORM
      DO I=1,N
        R(I)=0.0D0
      ENDDO

C     Solve the system
      DO ITER=1,MAXIT

C       Dense arrays
        IF(SPARSE_A.EQ.0) THEN

C         Calculate the residual, and take the forward sweep
          RESID=0.0D0
          DO I=1,N
            RTMP=B(I)-PDDOT(N,A(I),LDA,X,1)
            RESID=RESID+RTMP**2

            DSUM=0.0D0
C$OMP       PARALLEL DO PRIVATE(J), REDUCTION(+:DSUM)
            DO J=1,I-1
              DSUM=DSUM+A(I+LDA*(J-1))*R(J)*D(J)
            ENDDO
C$OMP       END PARALLEL DO
            R(I)=D(I)*(RTMP-DSUM)
          ENDDO
          RESID=SQRT(RESID)

C         Backwards sweep, and update the solution
          DO I=N,1,-1
            DSUM=PDDOT(N-I,A(I+LDA*I),LDA,R(I+1),1)
            R(I)=D(I)*( R(I) - D(I)*DSUM )

            X(I)=X(I)+R(I)
          ENDDO

C       Sparse arrays
        ELSE IF(SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4) THEN

C         Calculate the residual, and take the forward sweep
          RESID=0.0D0
          DO I=1,N
            RTMP=B(I)
            DO J=ISR_A(I),ISR_A(I+1)-1
              RTMP=RTMP-A(J)*X(ISC_A(J))
            ENDDO
            RESID=RESID+RTMP**2

            DSUM=0.0D0
            DO J=ISR_A(I),ISR_A(I+1)-1
              IF(ISC_A(J).GE.I) GOTO 100
              DSUM=DSUM+A(J)*R(ISC_A(J))*D(ISC_A(J))
            ENDDO
 100        CONTINUE
            R(I)=D(I)*(RTMP-DSUM)
          ENDDO
          RESID=SQRT(RESID)

C         Backwards sweep, and update the solution
          DO I=N,1,-1

            DSUM=0.0D0
            DO J=ISR_A(I+1)-1,ISR_A(I),-1
              IF(ISC_A(J).LE.I) GOTO 200
              DSUM=DSUM+A(J)*R(ISC_A(J))
            ENDDO
 200        CONTINUE
            R(I)=D(I)*( R(I) - D(I)*DSUM )

            X(I)=X(I)+R(I)
          ENDDO

        ELSE
          ERROR='>>Sparsity type not implemented'
          GOTO 9999
        ENDIF

        IF(ITER.EQ.1) SCALED_TOL=TOL*(ANORM*PDNRM2(N,X,1) + BNORM)

C       Print out solution, and check for convergance
        IF(OUTPUTCODE.GE.2) THEN
          WRITE(OP_STRING,'(1P,I5,1x,G12.6)') ITER,RESID
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        IF(RESID.LE.SCALED_TOL) GOTO 300

      ENDDO
      ITER=MAXIT

 300  CONTINUE

      RETURN

 9999 CALL ERRORS('INCOMP_CHOL0_SUB',ERROR)
      RETURN 1
      END
