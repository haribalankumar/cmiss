
      SUBROUTINE INCOMP_LDL0_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,
     '  ERROR,*)

C#### Subroutine: INCOMP_LDL0_FACTOR
C###  Description:
C###    INCOMP_LDL0_FACTOR factorises a system of equations for solution
C###    with the Incomplete LDL(0) iterative solver or preconditioner.
C###  Written by Stuart Norris 10/01/03

      IMPLICIT NONE
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*),LDA,N,SPARSE_A
      REAL*8 A(*),D(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J,II,IJ
      REAL*8 DSUM,RDUM


      DO I=1,N
        D(I)=0.0D0
      ENDDO

C     Dense matrix.
      IF(SPARSE_A.EQ.0) THEN

        DO I=1,N
          II=I+LDA*(I-1)
          IF(A(II).EQ.0.0D0) THEN
            ERROR='>>Zero on diagonal'
            GOTO 9999
          ENDIF
          RDUM=A(II)

          DSUM=0.0D0
          DO J=1,I-1
            IJ=I+LDA*(J-1)
            DSUM=DSUM+D(J)*A(IJ)**2
          ENDDO

          IF(RDUM.LE.DSUM) THEN
            IF(RDUM.EQ.DSUM) THEN
              ERROR='>>Zero decomposition element'
            ELSE
              ERROR='>>System is not positive definate'
            ENDIF
            GOTO 9999
          ENDIF
          D(I)=1.0D0/(A(II)-DSUM)
        ENDDO

C     Sparse matricies: CMISS compressed row
      ELSE IF(SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4) THEN

        DO I=1,N

          RDUM=0.0D0
          DSUM=0.0D0
          DO IJ=ISR_A(I),ISR_A(I+1)-1
            J=ISC_A(IJ)
            IF(J.GE.I) THEN
              IF(J.NE.I) THEN
                ERROR='>>Missing diagonal element'
                GOTO 9999
              ELSE IF(A(IJ).EQ.0.0D0) THEN
                ERROR='>>Zero on diagonal'
                GOTO 9999
              ENDIF
              RDUM=A(IJ)
              GOTO 200
            ENDIF
            DSUM=DSUM+D(J)*A(IJ)**2
          ENDDO
 200      CONTINUE

          IF(RDUM.LE.DSUM) THEN
            IF(RDUM.EQ.DSUM) THEN
              ERROR='>>Zero decomposition element'
            ELSE
              ERROR='>>System is not positive definate'
            ENDIF
            GOTO 9999
          ENDIF
          D(I)=1.0D0/(RDUM-DSUM)
        ENDDO

      ELSE
        ERROR='>>Sparsity type not implemented'
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('INCOMP_LDL0_FACTOR',ERROR)
      RETURN 1
      END


      SUBROUTINE PRECON_INCOMP_LDL0(A,LDA,N,X,B,D,R,SPARSE_A,ISC_A,
     '  ISR_A,MAXIT,ERROR,*)

C#### Subroutine: PRECON_INCOMP_LDL0
C###  Description:
C###    PRECON_INCOMP_LDL0 applies MAXIT iterations of the Incomplete
C###    LDL(0) factorisation solver to precondition the system
C###           A x = b
C###  Written by Stuart Norris 10/01/03

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

            DSUM=PDDOT(I-1,A(I),LDA,R,1)
            R(I)=D(I)*(RTMP-DSUM)
          ENDDO

C         Backwards sweep, and update the solution
          DO I=N,1,-1
            DSUM=PDDOT(N-I,A(I+LDA*I),LDA,R(I+1),1)
            R(I)=R(I)-D(I)*DSUM

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
              DSUM=DSUM+A(J)*R(ISC_A(J))
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

            R(I)=R(I)-D(I)*DSUM
            X(I)=X(I)+R(I)
          ENDDO

        ELSE
          ERROR='>>Sparsity type not implemented'
          GOTO 9999
        ENDIF

      ENDDO

      RETURN

 9999 CALL ERRORS('PRECON_INCOMP_LDL0',ERROR)
      RETURN 1
      END



      SUBROUTINE INCOMP_LDL0(A,LDA,N,X,B,D,SPARSE_A,ISC_A,ISR_A,ANORM,
     '  RESID,ITER,FACTORISE,OUTPUTCODE,ERROR,*)

C#### Subroutine: INCOMP_LDL0
C###  Description:
C###    INCOMP_LDL0 solves a system of linear equations using Incomplete
C###    LDL(0) factorisation. The A matrix is factorised using no infill 
C###    (or pivoting), and a sequence of forward and backward relaxation 
C###    passes are made to solve the system.
C###
C###    The level 0 version of the solver uses a very sparse factorisation
C###    where only the diagonal is stored.
C###  Written by Stuart Norris 10/01/03

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
!     Local Variables
      POINTER R_PTR
      DATA R_PTR / 0 /


      IF(OUTPUTCODE.GE.2) THEN
        WRITE(OP_STRING,'('' Incomplete LDL(0) Decomposition'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

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
        CALL INCOMP_LDL0_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,ERROR,
     '    *9999)
      ENDIF

C     Solve the system
      CALL ALLOCATE_MEMORY(N,1,DPTYPE,R_PTR,NOINIT,ERROR,*9999)

      CALL INCOMP_LDL0_SUB(A,LDA,N,X,B,D,%VAL(R_PTR),SPARSE_A,ISC_A,
     '  ISR_A,ANORM,RESID,ITER,OUTPUTCODE,ERROR,*9999)

      IF(R_PTR.NE.0) CALL FREE_MEMORY(R_PTR,ERROR,*9999)

C     Bye!
      RETURN

 9999 CALL ERRORS('INCOMP_LDL0',ERROR)
      RETURN 1
      END


      SUBROUTINE INCOMP_LDL0_SUB(A,LDA,N,X,B,D,R,SPARSE_A,ISC_A,ISR_A,
     '  ANORM,RESID,ITER,OUTPUTCODE,ERROR,*)

C#### Subroutine: INCOMP_LDL0_SUB
C###  Description:
C###    INCOMP_LDL0_SUB solves a system of linear equations using Incomplete
C###    LDL(0) factorisation. The A matrix is factorised using no infill 
C###    (or pivoting), and a sequence of forward and backward relaxation 
C###    passes are made to solve the system.
C###
C###    The level 0 version of the solver uses a very sparse factorisation
C###    where only the diagonal is stored.
C###  Written by Stuart Norris 10/01/03

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*),LDA,N,SPARSE_A,ITER,OUTPUTCODE
      REAL*8 A(*),X(*),B(*),D(*),R(*),ANORM,RESID
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,IJ,MAXIT
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

            DSUM=PDDOT(I-1,A(I),LDA,R,1)
            R(I)=(RTMP-DSUM)*D(I)
          ENDDO
          RESID=SQRT(RESID)

C         Backwards sweep, and update the solution
          DO I=N,1,-1
            DSUM=PDDOT(N-I,A(I+LDA*I),LDA,R(I+1),1)
            R(I)=R(I)-DSUM*D(I)

            X(I)=X(I)+R(I)
          ENDDO

C       Sparse arrays
        ELSE IF(SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4) THEN

C         Calculate the residual, and take the forward sweep
          RESID=0.0D0
          DO I=1,N
            RTMP=B(I)
            DO IJ=ISR_A(I),ISR_A(I+1)-1
              RTMP=RTMP-A(IJ)*X(ISC_A(IJ))
            ENDDO
            RESID=RESID+RTMP**2

            DSUM=0.0D0
            DO IJ=ISR_A(I),ISR_A(I+1)-1
              IF(ISC_A(IJ).GE.I) GOTO 100
              DSUM=DSUM+A(IJ)*R(ISC_A(IJ))
            ENDDO
 100        CONTINUE
            R(I)=D(I)*(RTMP-DSUM)
          ENDDO
          RESID=SQRT(RESID)

C         Backwards sweep, and update the solution
          DO I=N,1,-1
            DSUM=0.0D0
            DO IJ=ISR_A(I+1)-1,ISR_A(I),-1
              IF(ISC_A(IJ).LE.I) GOTO 200
              DSUM=DSUM+A(IJ)*R(ISC_A(IJ))
            ENDDO
 200        CONTINUE

            R(I)=R(I)-D(I)*DSUM
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

 9999 CALL ERRORS('INCOMP_LDL0_SUB',ERROR)
      RETURN 1
      END
