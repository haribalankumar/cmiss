
      SUBROUTINE CHOL_FACTOR(A,LDA,N,ERROR,*)

C#### Subroutine: CHOL_FACTOR
C###  Description:
C###    CHOL_FACTOR produces a Cholesky factorisation of a system of equations.
C###    The factorised system is stored at L^T in the upper triangular portion
C###    of A.
C###  Written by Stuart Norris 05/04/02

      IMPLICIT NONE
!     Parameter List
      INTEGER LDA,N
      REAL*8 A(LDA,*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J,K
      REAL*8 SUM
!     Functions
      REAL*8 DDOT
      EXTERNAL DDOT


      DO I=1,N
        SUM=0.0D0
        DO K=1,I-1
          SUM=SUM+A(K,I)**2
        ENDDO
        A(I,I)=A(I,I)-SUM

        IF(A(I,I).EQ.0.0D0) THEN
          ERROR='>>Zero on diagonal'
          GOTO 9999
        ELSE IF(A(I,I).LT.0.0D0) THEN
          ERROR='>>Matrix not positive definate'
          GOTO 9999
        ENDIF
        A(I,I)=1.0D0/SQRT(A(I,I))

        DO J=I+1,N
          SUM=0.0D0
          DO K=1,I-1
            SUM=SUM+A(K,J)*A(K,I)
          ENDDO
          A(I,J)=(A(I,J)-SUM)*A(I,I)
        ENDDO
      ENDDO

      RETURN

 9999 CALL ERRORS('CHOL_FACTOR',ERROR)
      RETURN 1
      END


      SUBROUTINE CHOL_SOLVE(A,LDA,N,X,B,ERROR,*)

C#### Subroutine: CHOL_SOLVE
C###  Description:
C###    CHOL_SOLVE solves a system that has been factorised using CHOL_FACTOR.
C###    The factorised system is stored at L^T in the upper triangular portion
C###    of A.
C###  Written by Stuart Norris 05/04/02

      IMPLICIT NONE
!     Parameter List
      INTEGER LDA,N
      REAL*8 A(LDA,*),X(*),B(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J
      REAL*8 SUM
!     Functions
      REAL*8 DDOT
      EXTERNAL DDOT


      DO I=1,N
        SUM=0.0D0
        DO J=1,I-1
          SUM=SUM+A(J,I)*X(J)
        ENDDO

        X(I)=(B(I)-SUM)*A(I,I)
      ENDDO

      DO I=N,1,-1
        SUM=0.0D0
        DO J=I+1,N
          SUM=SUM+A(I,J)*X(J)
        ENDDO

        X(I)=(X(I)-SUM)*A(I,I)
      ENDDO

      RETURN

 9999 CALL ERRORS('CHOL_SOLVE',ERROR)
      RETURN 1
      END


      SUBROUTINE LDL_FACTOR(A,LDA,N,ERROR,*)

C#### Subroutine: LDL_FACTOR
C###  Description:
C###    Factorise a symmetric system using LDL factorisation. The
C###    resulting system should be solved using LDL_SOLVE.
C###  Written by Stuart Norris 04/04/02

      IMPLICIT NONE
!     Parameter List
      INTEGER LDA,N
      REAL*8 A(LDA,*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER J
      REAL*8 AJJ
C     INTEGER I,K
C     REAL*8 SUM
!     Functions
      REAL*8 DDOT
      EXTERNAL DDOT


C     DO I=1,N
C       SUM=0.0D0
C       DO K=1,I-1
C         SUM=SUM+A(K,K)*A(K,I)**2
C       ENDDO
C       A(I,I)=A(I,I)-SUM
C
C       IF(ABS(A(I,I)).EQ.0.0D0) THEN
C         ERROR='>>Zero on diagonal'
C         GOTO 9999
C       ELSE IF(A(I,I).LT.0.0D0) THEN
C         ERROR='>>Matrix not positive definate'
C         GOTO 9999
C       ENDIF
C       A(I,I)=1.0D0/A(I,I)
C
C       DO J=I+1,N
C         SUM=0.0D0
C         DO K=1,I-1
C           SUM=SUM+A(K,J)*A(K,I)*A(K,K)
C         ENDDO
C         A(I,J)=A(I,J)-SUM
C       ENDDO
C     ENDDO


      AJJ=A(1,1)
      IF(AJJ.EQ.0.0D0) THEN
        ERROR='>>Zero on diagonal'
        GOTO 9999
      ENDIF
      A(1,1)=1.0D0/AJJ
      CALL DSCAL(N-1,A(1,1),A(2,1),1)

      DO J=2,N
C       SUM=0.0D0
C       DO K=1,J-1
C         SUM=SUM+A(J,K)*A(K,J)
C       ENDDO
C       AJJ=A(J,J)-SUM
        AJJ=A(J,J)-DDOT(J-1,A(J,1),LDA,A(1,J),1)
        IF(AJJ.EQ.0.0D0) THEN
          ERROR='>>Zero on diagonal'
          GOTO 9999
        ENDIF
        A(J,J)=1.0D0/AJJ

C       DO I=J+1,N
C         SUM=0.0D0
C         DO K=1,J-1
C           SUM=SUM+A(I,K)*A(K,J)
C         ENDDO
C         A(I,J)=(A(I,J)-SUM)*A(J,J)
C       ENDDO
        IF(J.LT.N) THEN
          CALL DGEMV('N',N-J,J-1,-A(J,J),A(J+1,1),LDA,
     '      A(1,J),1,A(J,J),A(J+1,J),1)
        ENDIF

C       DO I=J+1,N
C         A(J,I)=A(I,J)*AJJ
C       ENDDO
        IF(J.LT.N) THEN
          CALL DCOPY(N-J,A(J+1,J),1,A(J,J+1),LDA)
          CALL DSCAL(N-J,AJJ,A(J,J+1),LDA)
        ENDIF
      ENDDO

      RETURN

 9999 CALL ERRORS('LDL_FACTOR',ERROR)
      RETURN 1
      END


      SUBROUTINE LDL_SOLVE(A,LDA,N,X,B,ERROR,*)

C#### Subroutine: LDL_SOLVE
C###  Description:
C###    Solve a symmetric system of equations using LDL decomposition.
C###    The system must be symmetric, and should have been factorised
C###    using LDL_FACTOR. The factorised system is stored as L^T in
C###    the matrix A.
C###  Written by Stuart Norris 04/04/02

      IMPLICIT NONE
!     Parameter List
      INTEGER LDA,N
      REAL*8 A(LDA,*),X(*),B(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I
C     INTEGER J
C     REAL*8 SUM
!     Functions
      REAL*8 DDOT
      EXTERNAL DDOT


      X(1)=B(1)*A(1,1)
      DO I=2,N
C       SUM=0.0D0
C       DO J=1,I-1
C         SUM=SUM+A(J,I)*X(J)
C       ENDDO
C       X(I)=(B(I)-SUM)*A(I,I)
        X(I)=(B(I)-DDOT(I-1,A(1,I),1,X(1),1))*A(I,I)
      ENDDO

CC    DO I=N,1,-1
CC      SUM=0.0D0
CC      DO J=I+1,N
CC        SUM=SUM+A(I,J)*X(J)
CC      ENDDO
CC      X(I)=X(I)-SUM*A(I,I)
CC    ENDDO

      DO I=N-1,1,-1
C       SUM=0.0D0
C       DO J=I+1,N
C         SUM=SUM+A(J,I)*X(J)
C       ENDDO
C       X(I)=X(I)-SUM
        X(I)=X(I)-DDOT(N-I,A(I+1,I),1,X(I+1),1)
      ENDDO

      RETURN

 9999 CALL ERRORS('LDL_SOLVE',ERROR)
      RETURN 1
      END


      SUBROUTINE LU_FACTOR(A,LDA,N,ERROR,*)

C#### Subroutine: LU_FACTOR
C###  Description:
C###    Factorise a linear system using LU factorisation. The
C###    resulting system should be solved using LU_SOLVE.
C###  Written by Stuart Norris 05/04/02

      IMPLICIT NONE
!     Parameter List
      INTEGER LDA,N
      REAL*8 A(LDA,*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J,K
      REAL*8 SUM


      DO J=1,N
        DO I=1,J
          SUM=0.0D0
          DO K=1,I-1
            SUM=SUM+A(I,K)*A(K,J)
          ENDDO

          A(I,J)=A(I,J)-SUM
        ENDDO

        IF(ABS(A(J,J)).EQ.0.0D0) THEN
          ERROR='>>Zero on diagonal'
          GOTO 9999
        ENDIF
        A(J,J)=1.0D0/A(J,J)

        DO I=J+1,N
          SUM=0.0D0
          DO K=1,J-1
            SUM=SUM+A(I,K)*A(K,J)
          ENDDO

          A(I,J)=(A(I,J)-SUM)*A(J,J)
        ENDDO
      ENDDO

      RETURN

 9999 CALL ERRORS('LU_FACTOR',ERROR)
      RETURN 1
      END


      SUBROUTINE LU_SOLVE(A,LDA,N,X,B,ERROR,*)

C#### Subroutine: LU_SOLVE
C###  Description:
C###    Solve a symmetric system of equations using LU decomposition.
C###    The system should have been factorised using LU_FACTOR. The 
C###    factorised system is stored as L+U in the matrix A.
C###  Written by Stuart Norris 05/04/02

      IMPLICIT NONE
!     Parameter List
      INTEGER LDA,N
      REAL*8 A(LDA,*),X(*),B(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J
      REAL*8 SUM
!     Functions
      REAL*8 DDOT
      EXTERNAL DDOT


      DO I=1,N
        SUM=0.0D0
        DO J=1,I-1
          SUM=SUM+A(I,J)*X(J)
        ENDDO

        X(I)=B(I)-SUM
      ENDDO
      
      DO I=N,1,-1
        SUM=0.0D0
        DO J=I+1,N
          SUM=SUM+A(I,J)*X(J)
        ENDDO

        X(I)=(X(I)-SUM)*A(I,I)
      ENDDO

      RETURN

 9999 CALL ERRORS('LU_SOLVE',ERROR)
      RETURN 1
      END


      SUBROUTINE SVD_SOLVE(M,N,SVDTOL,W,U,LDU,VT,LDVT,X,B,WORK,LWORK,
     '  ERROR,*)

C#### Subroutine: SVD_SOLVE
C###  Description:
C###    SVD_SOLVE solves a dense linear system that has been factorised
C###    using the LAPACK DGESVD routine. The original loops have been 
C###    replaced by BLAS.
C###  Written by Stuart Norris 19/12/01

      IMPLICIT NONE

!     Parameter List
      INTEGER M,N,LDU,LDVT,LWORK
      REAL*8 SVDTOL,U(LDU,N),VT(LDVT,N),W(*),X(N),B(M),WORK(N)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER J
!     Functions
      REAL*8 DDOT


      CALL ASSERT(LWORK.GE.N,'>> must have LWORK >= N',ERROR,*9999)
      CALL ASSERT(SVDTOL.GT.0.D0,'>> must have SVDTOL >= 0',ERROR,*9999)

C     Forward substitute.

C$OMP PARALLEL DO PRIVATE(J,M)
      DO J=1,N
        IF(J.GT.M) THEN
          WORK(J)=0.0D0
        ELSE IF(W(J).LT.ABS(SVDTOL)) THEN
          WORK(J)=0.0D0
        ELSE
C         WORK(J)=0.0D0
C         DO I=1,M
C           WORK(J)=WORK(J)+U(I,J)*B(I)
C         ENDDO
C         WORK(J)=WORK(J)/W(J)
          WORK(J)=DDOT(M,U(1,J),1,B,1)/W(J)
        ENDIF
      ENDDO
C$OMP END PARALLEL DO

C     Backward substitute.

C     DO J=1,N
C       X(J)=0.0D0
C       DO I=1,N
C         X(J)=X(J)+VT(I,J)*WORK(I)
C       ENDDO
C     ENDDO
      CALL DGEMV('T',N,N,1.0D0,VT,LDVT,WORK,1,0.0D0,X,1)

      RETURN
 9999 CALL ERRORS('SVD_SOLVE',ERROR)
      RETURN 1
      END
