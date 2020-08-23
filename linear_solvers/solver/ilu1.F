
      SUBROUTINE INCOMP_LU1_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,
     '  ERROR,*)

C#### Subroutine: INCOMP_LU1_FACTOR
C###  Description:
C###    INCOMP_LU1_FACTOR Factorises a system of equations for solution
C###    with Incomplete LU (1) factorisation.
C###
C###    Note that there is not much point using this routine for a dense
C###    matrix (it uses the same number of operations as dense LU), but
C###    it is included to test and debug the sparse matrix operations.
C###  Written by Stuart Norris 21/05/02

      IMPLICIT NONE
      INCLUDE 'mach00.inc'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'solver.inc'
!     Parameter List
      INTEGER LDA,N,SPARSE_A
      INTEGER ISC_A(*),ISR_A(*)
      REAL*8 A(*),D(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER NZA,IDUM(1),JDUM(1),KDUM(1),LDUM(1)
      POINTER IDX_PTR,ITC_PTR,ITR_PTR,IWK_PTR
      DATA IDX_PTR,ITC_PTR,ITR_PTR,IWK_PTR / 4*0 /


C     Dense matrix.
      IF(SPARSE_A.EQ.0) THEN
        CALL INCOMP_LU1_FACTOR_SUB(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,
     '    IDUM,JDUM,KDUM,LDUM,ERROR,*9999)

C     Sparse matricies: CMISS compressed row
      ELSE IF(SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4) THEN
        NZA=ISR_A(N+1)-1

        CALL ALLOCATE_MEMORY(NZA,1,INTTYPE,IDX_PTR,NOINIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(N+1,1,INTTYPE,ITC_PTR,NOINIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(NZA,1,INTTYPE,ITR_PTR,NOINIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(  N,1,INTTYPE,IWK_PTR,NOINIT,ERROR,*9999)

        CALL INCOMP_LU1_FACTOR_SUB(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,
     '    %VAL(ITC_PTR),%VAL(ITR_PTR),%VAL(IDX_PTR),%VAL(IWK_PTR),
     '    ERROR,*9999)

        IF(IDX_PTR.NE.0) CALL FREE_MEMORY(IDX_PTR,ERROR,*9999)
        IF(ITC_PTR.NE.0) CALL FREE_MEMORY(ITC_PTR,ERROR,*9999)
        IF(ITR_PTR.NE.0) CALL FREE_MEMORY(ITR_PTR,ERROR,*9999)
        IF(IWK_PTR.NE.0) CALL FREE_MEMORY(IWK_PTR,ERROR,*9999)

      ELSE
        ERROR='>>Sparsity type not implemented'
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('INCOMP_LU1_FACTOR',ERROR)
      RETURN 1
      END


      SUBROUTINE INCOMP_LU1_FACTOR_SUB(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,
     '  ITC_A,ITR_A,IDX_A,IWORK,ERROR,*)

C#### Subroutine: INCOMP_LU1_FACTOR_SUB
C###  Description:
C###    INCOMP_LU1_FACTOR_SUB Factorises a system of equations for solution
C###    with Incomplete LU (1) factorisation.
C###
C###    Note that there is not much point using this routine for a dense
C###    matrix (it uses the same number of operations as dense LU), but
C###    it is included to test and debug the sparse matrix operations.
C###  Written by Stuart Norris 21/05/02

      IMPLICIT NONE
!     Parameter List
      INTEGER LDA,N,SPARSE_A,IWORK(*)
      INTEGER ISC_A(*),ISR_A(*),ITC_A(*),ITR_A(*),IDX_A(*)
      REAL*8 A(*),D(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J,IJ,IK,KJ,JJ,J1,KK,K1,ICOL,NZA
      REAL*8 SUM
!     Functions
      REAL*8 PDDOT
      EXTERNAL PDDOT


C     Dense matrix
      IF(SPARSE_A.EQ.0) THEN

        CALL DCOPY(LDA*N,A,1,D,1)

        DO J=1,N
          J1=1+(J-1)*LDA
          IJ=J1
          DO I=1,J
            IJ=I+(J-1)*LDA
            IF(A(IJ).NE.0.0D0) THEN
C             SUM=0.0D0
C             DO K=1,I-1
C               SUM=SUM+D(I,K)*D(K,J)
C             ENDDO
              SUM=PDDOT(I-1,D(I),LDA,D(J1),1)

              D(IJ)=D(IJ)-SUM
            ENDIF
          ENDDO
          JJ=IJ

          IF(D(JJ).EQ.0.0D0) THEN
            WRITE(ERROR,'(A,2(1X,I9))') '>>Zero on diagonal:',J,N
            GOTO 9999
          ENDIF
          D(JJ)=1.0D0/D(JJ)

          DO I=J+1,N
            IJ=I+(J-1)*LDA
            IF(A(IJ).NE.0.0D0) THEN
C             SUM=0.0D0
C             DO K=1,J-1
C               SUM=SUM+D(I,K)*D(K,J)
C             ENDDO
              SUM=PDDOT(J-1,D(I),LDA,D(J1),1)

              D(IJ)=(D(IJ)-SUM)*D(JJ)
            ENDIF
          ENDDO
        ENDDO

C     Sparse matricies: CMISS compressed row
      ELSE IF(SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4) THEN

C       Copy A onto D, and generate a compressed column index for 
C       use in the factorisation.
        NZA=ISR_A(N+1)-1
        CALL DCOPY(NZA,A,1,D,1)
        CALL TRANS_INDEX(N,NZA,ISC_A,ISR_A,ITC_A,ITR_A,IDX_A,IWORK,
     '    ERROR,*9999)

        DO J=1,N

C         Set the column index
          IF(J.EQ.1) THEN
            DO I=1,N
              IWORK(I)=0
            ENDDO
          ELSE
            DO I=ITC_A(J-1),ITC_A(J)-1
              IWORK(ITR_A(I))=0
            ENDDO
          ENDIF
          DO I=ITC_A(J),ITC_A(J+1)-1
            IWORK(ITR_A(I))=IDX_A(I)
          ENDDO

C         Factor
          DO KK=ITC_A(J),ITC_A(J+1)-1
            I=ITR_A(KK)
            IF(I.GT.J) GOTO 200
            IJ=IDX_A(KK)

            SUM=0.0D0
C           DO K=1,I-1
            DO IK=ISR_A(I),ISR_A(I+1)-1
              IF(ISC_A(IK).GE.I) GOTO 100
              ICOL=ISC_A(IK)
              KJ=IWORK(ICOL)
              IF(KJ.NE.0) SUM=SUM+D(IK)*D(KJ)
            ENDDO
 100        CONTINUE
            D(IJ)=D(IJ)-SUM
          ENDDO
 200      CONTINUE
          K1=KK

          JJ=IWORK(J)
          IF(JJ.EQ.0) THEN
            WRITE(ERROR,'(A,2(1X,I9))') '>>No diagonal element:',J,N
            GOTO 9999
          ELSE IF(D(JJ).EQ.0.0D0) THEN
            WRITE(ERROR,'(A,2(1X,I9))') '>>Zero on diagonal:',J,N
            GOTO 9999
          ENDIF
          D(JJ)=1.0D0/D(JJ)

          DO KK=K1,ITC_A(J+1)-1
            IF(ITR_A(KK).GT.J) THEN
              IJ=IDX_A(KK)
              I=ITR_A(KK)

              SUM=0.0D0
C             DO K=1,J-1
              DO IK=ISR_A(I),ISR_A(I+1)-1
                IF(ISC_A(IK).GE.J) GOTO 300
                ICOL=ISC_A(IK)
                KJ=IWORK(ICOL)
                IF(KJ.NE.0) SUM=SUM+D(IK)*D(KJ)
              ENDDO
 300          CONTINUE
              D(IJ)=(D(IJ)-SUM)*D(JJ)
            ENDIF
          ENDDO

        ENDDO

      ELSE
        ERROR='>>Sparsity type not implemented'
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('INCOMP_LU1_FACTOR_SUB',ERROR)
      RETURN 1
      END


      SUBROUTINE PRECON_INCOMP_LU1(A,LDA,N,X,B,D,R,SPARSE_A,ISC_A,ISR_A,
     '  MAXIT,ERROR,*)

C#### Subroutine: PRECON_INCOMP_LU1
C###  Description:
C###    PRECON_INCOMP_LU1 applies MAXIT iterations of the Incomplete 
C###    LU (1) solver to precondition the system
C###           A x = b
C###  Written by Stuart Norris 21/05/02

      IMPLICIT NONE
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*),LDA,N,SPARSE_A,MAXIT
      REAL*8 A(*),X(*),B(*),D(*),R(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,II,IJ,IP,ITER
      REAL*8 DSUM,DIAG,RTMP
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
            R(I)=B(I)-PDDOT(N,A(I),LDA,X,1)

            DSUM=PDDOT(I-1,D(I),LDA,R(1),1)

            R(I)=R(I)-DSUM
          ENDDO

C         Backwards sweep, and update the solution
          DO I=N,1,-1
            II=I+(I-1)*LDA
            IP=I+I*LDA

            DSUM=PDDOT(N-I,D(IP),LDA,R(I+1),1)
            DIAG=D(II)

            R(I)=(R(I)-DSUM)*DIAG
            X(I)=X(I)+R(I)
          ENDDO

C       Sparse arrays
        ELSE IF(SPARSE_A.EQ.1.OR.SPARSE_A.EQ.4) THEN

C         Calculate the residual, and take the forward sweep
          DO I=1,N
            RTMP=0.0D0
            DO IJ=ISR_A(I),ISR_A(I+1)-1
              RTMP=RTMP+A(IJ)*X(ISC_A(IJ))
            ENDDO
            RTMP=B(I)-RTMP

            DSUM=0.0D0
            IF(I.GT.1) THEN
              DO IJ=ISR_A(I),ISR_A(I+1)-1
                IF(ISC_A(IJ).GE.I) GOTO 100
                DSUM=DSUM+D(IJ)*R(ISC_A(IJ))
              ENDDO
 100          CONTINUE
            ENDIF

            R(I)=RTMP-DSUM
          ENDDO

C         Backwards sweep, and update the solution
          DO I=N,1,-1
            DSUM=0.0D0
            DIAG=1.0D0
            DO IJ=ISR_A(I+1)-1,ISR_A(I),-1
              IF(ISC_A(IJ).LE.I) THEN
                IF(ISC_A(IJ).EQ.I) DIAG=D(IJ)
                GOTO 200
              ENDIF
              DSUM=DSUM+D(IJ)*R(ISC_A(IJ))
            ENDDO
 200        CONTINUE

            R(I)=(R(I)-DSUM)*DIAG
            X(I)=X(I)+R(I)
          ENDDO

        ELSE
          ERROR='>>Sparsity type not implemented'
          GOTO 9999
        ENDIF

      ENDDO

      RETURN

 9999 CALL ERRORS('PRECON_INCOMP_LU1',ERROR)
      RETURN 1
      END


      SUBROUTINE INCOMP_LU1(A,LDA,N,X,B,D,SPARSE_A,ISC_A,ISR_A,ANORM,
     '  RESID,ITER,FACTORISE,OUTPUTCODE,ERROR,*)

C#### Subroutine: INCOMP_LU1
C###  Description:
C###    INCOMP_LU1 solves a system of linear equations using Incomplete
C###    LU (1) factorisation The A matrix is decomposed using no infill 
C###    (or pivoting), and then the a sequence of forward and backward
C###    relaxation passes are made.
C###  Written by Stuart Norris 21/03/02

      IMPLICIT NONE
      INCLUDE 'mach00.inc'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'solver.inc'
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*)
      INTEGER LDA,N,SPARSE_A,ITER,OUTPUTCODE
      REAL*8 A(*),X(*),B(*),D(*),ANORM,RESID
      CHARACTER ERROR*(*)
      LOGICAL FACTORISE
!     Local Variables
      POINTER R_PTR
      DATA R_PTR / 0 /


      IF(OUTPUTCODE.GE.2) THEN
        WRITE(OP_STRING,'('' Incomplete LU(1) Factorisation'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

C     Check input
      IF(ITER.LE.0) THEN
        ERROR='>>Maxit less than 0'
        GOTO 9999
      ENDIF
      IF(RESID.LT.0.0d0) THEN
        ERROR='>>Tolerance less than 0.0'
        GOTO 9999
      ENDIF

C     Factorise the system
      IF(FACTORISE) THEN
        CALL INCOMP_LU1_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,ERROR,
     '    *9999)
      ENDIF

C     Solve the system
      CALL ALLOCATE_MEMORY(N,1,DPTYPE,R_PTR,NOINIT,ERROR,*9999)

      CALL INCOMP_LU1_SUB(A,LDA,N,X,B,D,%VAL(R_PTR),SPARSE_A,
     '  ISC_A,ISR_A,ANORM,RESID,ITER,OUTPUTCODE,ERROR,*9999)

      IF(R_PTR.NE.0) CALL FREE_MEMORY(R_PTR,ERROR,*9999)

C     Bye!
      RETURN

 9999 CALL ERRORS('INCOMP_LU1',ERROR)
      RETURN 1
      END


      SUBROUTINE INCOMP_LU1_SUB(A,LDA,N,X,B,D,R,SPARSE_A,ISC_A,ISR_A,
     '  ANORM,RESID,ITER,OUTPUTCODE,ERROR,*)

C#### Subroutine: INCOMP_LU1
C###  Description:
C###    INCOMP_LU1 solves a system of linear equations using Incomplete
C###    LU (1) factorisation The A matrix is decomposed using no infill 
C###    (or pivoting), and then the a sequence of forward and backward
C###    relaxation passes are made.
C###  Written by Stuart Norris 21/03/02

      IMPLICIT NONE
      INCLUDE 'mach00.inc'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*)
      INTEGER LDA,N,SPARSE_A,ITER,OUTPUTCODE
      REAL*8 A(*),X(*),B(*),D(*),R(*),ANORM,RESID
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,II,IJ,IP,MAXIT
      REAL*8 TOL,DSUM,DIAG,RTMP,BNORM,SCALED_TOL
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

C         Calculate the residual and take the forward sweep
          DO I=1,N
            R(I)=B(I)-PDDOT(N,A(I),LDA,X,1)
          ENDDO
          RESID=PDNRM2(N,R,1)

          DO I=1,N
            DSUM=PDDOT(I-1,D(I),LDA,R,1)

            R(I)=R(I)-DSUM
          ENDDO

C         Backwards sweep, and update the solution
          DO I=N,1,-1
            II=I+(I-1)*LDA
            IP=I+I*LDA

            DSUM=PDDOT(N-I,D(IP),LDA,R(I+1),1)
            DIAG=D(II)

            R(I)=(R(I)-DSUM)*DIAG
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
            IF(I.GT.1) THEN
              DO IJ=ISR_A(I),ISR_A(I+1)-1
                IF(ISC_A(IJ).GE.I) GOTO 100
                DSUM=DSUM+D(IJ)*R(ISC_A(IJ))
              ENDDO
 100          CONTINUE
            ENDIF
            R(I)=RTMP-DSUM
          ENDDO
          RESID=SQRT(RESID)

C         Backwards sweep, and update the solution
          DO I=N,1,-1
            DSUM=0.0D0
            DIAG=1.0D0
            DO IJ=ISR_A(I+1)-1,ISR_A(I),-1
              IF(ISC_A(IJ).LE.I) THEN
                IF(ISC_A(IJ).EQ.I) DIAG=D(IJ)
                GOTO 200
              ENDIF
              DSUM=DSUM+D(IJ)*R(ISC_A(IJ))
            ENDDO
 200        CONTINUE

            R(I)=(R(I)-DSUM)*DIAG
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

 9999 CALL ERRORS('INCOMP_LU1_SUB',ERROR)
      RETURN 1
      END
