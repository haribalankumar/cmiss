
      SUBROUTINE INCOMP_LDL1_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,
     '  ERROR,*)

C#### Subroutine: INCOMP_LDL1_FACTOR
C###  Description:
C###    INCOMP_LDL1_FACTOR Factorises a system of equations for solution
C###    with Incomplete LDL(1) factorisation
C###  Written by Stuart Norris 28/03/02

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
      INTEGER NZA,IDUM(1),JDUM(1)
      POINTER IDX_PTR,IWK_PTR
      DATA IDX_PTR,IWK_PTR / 2*0 /


C     Dense matrix.
      IF(SPARSE_A.EQ.0) THEN
        CALL INCOMP_LDL1_FACTOR_SUB(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,
     '    IDUM,JDUM,ERROR,*9999)

C     Sparse matricies: CMISS compressed row.
      ELSE IF(SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4) THEN
        NZA=ISR_A(N+1)-1
        CALL ALLOCATE_MEMORY(NZA,1,INTTYPE,IDX_PTR,NOINIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(  N,1,INTTYPE,IWK_PTR,NOINIT,ERROR,*9999)

        CALL INCOMP_LDL1_FACTOR_SUB(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,
     '    %VAL(IDX_PTR),%VAL(IWK_PTR),ERROR,*9999)

        IF(IDX_PTR.NE.0) CALL FREE_MEMORY(IDX_PTR,ERROR,*9999)
        IF(IWK_PTR.NE.0) CALL FREE_MEMORY(IWK_PTR,ERROR,*9999)

      ELSE
        ERROR='>>Sparsity type not implemented'
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('INCOMP_LDL1_FACTOR',ERROR)
      RETURN 1
      END


      SUBROUTINE INCOMP_LDL1_FACTOR_SUB(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,
     '  IDX_A,IWORK,ERROR,*)

C#### Subroutine: INCOMP_LDL1_FACTOR_SUB
C###  Description:
C###    INCOMP_LDL1_FACTOR_SUB Factorises a system of equations for solution
C###    with Incomplete LDL(1) factorisation
C###  Written by Stuart Norris 28/03/02

      IMPLICIT NONE
!     Parameter List
      INTEGER LDA,N,SPARSE_A
      INTEGER ISC_A(*),ISR_A(*),IDX_A(*),IWORK(*)
      REAL*8 A(*),D(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J,K,IJ,IK,JI,JJ,JK,KJ,NZA
      REAL*8 SUM,DJJ
!     Functions
      REAL*8 DDOT
      EXTERNAL DDOT


C     Dense matrix.
      IF(SPARSE_A.EQ.0) THEN

        CALL DCOPY(LDA*N,A,1,D,1)

        DO J=1,N
          SUM=0.0D0
          DO K=1,J-1
            JK=J+(K-1)*LDA
            KJ=K+(J-1)*LDA
            SUM=SUM+D(JK)*D(KJ)
          ENDDO
          JJ=J+(J-1)*LDA
          DJJ=D(JJ)-SUM

          IF(A(JJ).EQ.0.0D0) THEN
            ERROR='>>No diagonal element'
            GOTO 9999
          ELSE IF(DJJ.EQ.0.0D0) THEN
            ERROR='>>Zero on diagonal'
            GOTO 9999
          ENDIF
          D(JJ)=1.0D0/DJJ

          DO I=J+1,N
            IJ=I+(J-1)*LDA
            JI=J+(I-1)*LDA
            IF(A(JI).NE.0.0D0) THEN
              SUM=0.0D0
              DO K=1,J-1
                IK=I+(K-1)*LDA
                KJ=K+(J-1)*LDA
                SUM=SUM+D(IK)*D(KJ)
              ENDDO

              D(JI)=(D(JI)-SUM)*D(JJ)
              D(IJ)=D(JI)*DJJ
            ENDIF
          ENDDO

        ENDDO

C     Sparse matricies: CMISS compressed row.
      ELSE IF(SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4) THEN

C       Copy A onto D, and generate a compressed column index for 
C       use in the factorisation.

        NZA=ISR_A(N+1)-1
        CALL DCOPY(NZA,A,1,D,1)
        CALL TRANS_INDEX_SYMM(N,NZA,ISC_A,ISR_A,IDX_A,IWORK,ERROR,*9999)

        DO J=1,N

          IF(J.EQ.1) THEN
            DO I=1,N
              IWORK(I)=0
            ENDDO
          ELSE
            DO I=ISR_A(J-1),ISR_A(J)-1
              IWORK(ISC_A(I))=0
            ENDDO
          ENDIF
          DO I=ISR_A(J),ISR_A(J+1)-1
            IWORK(ISC_A(I))=IDX_A(I)
          ENDDO

          SUM=0.0D0
          DO JK=ISR_A(J),ISR_A(J+1)-1
            IF(ISC_A(JK).GE.J) GOTO 100
            K=ISC_A(JK)
            KJ=IWORK(K)
            IF(KJ.NE.0) SUM=SUM+D(JK)*D(KJ)
          ENDDO
 100      CONTINUE
          JJ=IWORK(J)
          DJJ=D(JJ)-SUM

          IF(JJ.EQ.0.0D0) THEN
            ERROR='>>No diagonal element'
            GOTO 9999
          ELSE IF(DJJ.EQ.0.0D0) THEN
            ERROR='>>Zero on diagonal'
            GOTO 9999
          ENDIF
          D(JJ)=1.0D0/DJJ

          DO JI=ISR_A(J+1)-1,ISR_A(J),-1
            IF(ISC_A(JI).LE.J) GOTO 300
            I=ISC_A(JI)
            IJ=IWORK(I)
            SUM=0.0D0
            DO IK=ISR_A(I),ISR_A(I+1)-1
              IF(ISC_A(IK).GE.J) GOTO 200
              K=ISC_A(IK)
              KJ=IWORK(K)
              IF(KJ.NE.0) SUM=SUM+D(IK)*D(KJ)
            ENDDO
 200        CONTINUE

            D(JI)=(D(JI)-SUM)*D(JJ)
            D(IJ)=D(JI)*DJJ
          ENDDO
 300      CONTINUE

        ENDDO

      ELSE
        ERROR='>>Sparsity type not implemented'
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('INCOMP_LDL1_FACTOR_SUB',ERROR)
      RETURN 1
      END


      SUBROUTINE PRECON_INCOMP_LDL1(A,LDA,N,X,B,D,R,SPARSE_A,ISC_A,
     '  ISR_A,MAXIT,ERROR,*)

C#### Subroutine: PRECON_INCOMP_LDL1
C###  Description:
C###    PRECON_INCOMP_LDL1 applies MAXIT iterations of the Incomplete 
C###    LDL(1) solver to precondition the system
C###           A x = b
C###  Written by Stuart Norris 08/04/02

      IMPLICIT NONE
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*),LDA,N,SPARSE_A,MAXIT
      REAL*8 A(*),X(*),B(*),D(*),R(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,IJ,ITER
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
            RTMP=B(I)-PDDOT(N,A(I),LDA,X,1)

            DSUM=PDDOT(I-1,D(I),LDA,R,1)

            R(I)=(RTMP-DSUM)*D(I+(I-1)*LDA)
          ENDDO

C         Backwards sweep, and update the solution
          DO I=N,1,-1
            DSUM=PDDOT(N-I,D(I+I*LDA),LDA,R(I+1),1)
            R(I)=R(I)-DSUM

            X(I)=X(I)+R(I)
          ENDDO

C       Sparse arrays
        ELSE IF(SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4) THEN

C         Calculate the residual, and take the forward sweep
          DO I=1,N
            RTMP=B(I)
            DO IJ=ISR_A(I),ISR_A(I+1)-1
              RTMP=RTMP-A(IJ)*X(ISC_A(IJ))
            ENDDO

            DSUM=0.0D0
            DIAG=1.0D0
            DO IJ=ISR_A(I),ISR_A(I+1)-1
              IF(ISC_A(IJ).GE.I) THEN
                DIAG=D(IJ)
                GOTO 100
              ENDIF
              DSUM=DSUM+D(IJ)*R(ISC_A(IJ))
            ENDDO
 100        CONTINUE

            R(I)=(RTMP-DSUM)*DIAG
          ENDDO

C         Backwards sweep, and update the solution
          DO I=N,1,-1
            DSUM=0.0D0
            DO IJ=ISR_A(I+1)-1,ISR_A(I),-1
              IF(ISC_A(IJ).LE.I) GOTO 200
              DSUM=DSUM+D(IJ)*R(ISC_A(IJ))
            ENDDO
 200        CONTINUE

            R(I)=R(I)-DSUM
            X(I)=X(I)+R(I)
          ENDDO

        ELSE
          ERROR='>>Sparsity type not implemented'
          GOTO 9999
        ENDIF

      ENDDO

      RETURN

 9999 CALL ERRORS('PRECON_INCOMP_LDL1',ERROR)
      RETURN 1
      END



      SUBROUTINE INCOMP_LDL1(A,LDA,N,X,B,D,SPARSE_A,ISC_A,ISR_A,ANORM,
     '  RESID,ITER,FACTORISE,OUTPUTCODE,ERROR,*)

C#### Subroutine: INCOMP_LDL1
C###  Description:
C###    INCOMP_LDL1 solves a system of linear equations using Incomplete
C###    LDL(1) factorisation The A matrix is decomposed using no infill 
C###    (or pivoting), and then the a sequence of forward and backward
C###    relaxation passes are made.
C###  Written by Stuart Norris 21/03/02

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
        WRITE(OP_STRING,'('' Incomplete LDL(1) Factorisation'')')
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
        CALL INCOMP_LDL1_FACTOR(A,LDA,N,D,SPARSE_A,ISC_A,ISR_A,ERROR,
     '    *9999)
      ENDIF

C     Solve the system
      CALL ALLOCATE_MEMORY(N,1,DPTYPE,R_PTR,NOINIT,ERROR,*9999)

      CALL INCOMP_LDL1_SUB(A,LDA,N,X,B,D,%VAL(R_PTR),SPARSE_A,ISC_A,
     '  ISR_A,ANORM,RESID,ITER,OUTPUTCODE,ERROR,*9999)

      IF(R_PTR.NE.0) CALL FREE_MEMORY(R_PTR,ERROR,*9999)

C     Bye!
      RETURN

 9999 CALL ERRORS('INCOMP_LDL1',ERROR)
      RETURN 1
      END


      SUBROUTINE INCOMP_LDL1_SUB(A,LDA,N,X,B,D,R,SPARSE_A,ISC_A,ISR_A,
     '  ANORM,RESID,ITER,OUTPUTCODE,ERROR,*)

C#### Subroutine: INCOMP_LDL1_SUB
C###  Description:
C###    INCOMP_LDL1_SUB solves a system of linear equations using Incomplete
C###    LDL(1) factorisation The A matrix is decomposed using no infill 
C###    (or pivoting), and then the a sequence of forward and backward
C###    relaxation passes are made.
C###  Written by Stuart Norris 21/03/02

      IMPLICIT NONE
      INCLUDE 'mach00.inc'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*),LDA,N,SPARSE_A,ITER,OUTPUTCODE
      REAL*8 A(*),X(*),B(*),D(*),R(*),ANORM,RESID
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,II,IJ,MAXIT
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

C         Calculate the residual, and take the forward sweep
          RESID=0.0D0
          DO I=1,N
            RTMP=B(I)-PDDOT(N,A(I),LDA,X,1)
            RESID=RESID+RTMP**2

C           DSUM=0.0D0
C           DO J=1,I-1
C             IJ=I+(J-1)*LDA
C             DSUM=DSUM+D(IJ)*R(J)
C           ENDDO
            DSUM=PDDOT(I-1,D(I),LDA,R(1),1)
            II=I+(I-1)*LDA

            R(I)=(RTMP-DSUM)*D(II)
          ENDDO
          RESID=SQRT(RESID)

C         Backwards sweep, and update the solution
          DO I=N,1,-1
C           DSUM=0.0D0
C           DO J=I+1,N
C             IJ=I+(J-1)*LDA
C             DSUM=DSUM+D(IJ)*R(J)
C           ENDDO
            DSUM=PDDOT(N-I,D(I+I*LDA),LDA,R(I+1),1)

            R(I)=R(I)-DSUM

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
            DIAG=1.0D0
            DO IJ=ISR_A(I),ISR_A(I+1)-1
              IF(ISC_A(IJ).GE.I) THEN
                DIAG=D(IJ)
                GOTO 100
              ENDIF
              DSUM=DSUM+D(IJ)*R(ISC_A(IJ))
            ENDDO
 100        CONTINUE

            R(I)=(RTMP-DSUM)*DIAG
          ENDDO
          RESID=SQRT(RESID)

C         Backwards sweep, and update the solution
          DO I=N,1,-1
            DSUM=0.0D0
            DO IJ=ISR_A(I+1)-1,ISR_A(I),-1
              IF(ISC_A(IJ).LE.I) GOTO 200
              DSUM=DSUM+D(IJ)*R(ISC_A(IJ))
            ENDDO
 200        CONTINUE

            R(I)=R(I)-DSUM
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

 9999 CALL ERRORS('INCOMP_LDL1_SUB',ERROR)
      RETURN 1
      END
