      SUBROUTINE MP_PSEUDO_INVERSE(A,M,N,equations,used_N,
     '  used_equations,ERROR,*)

C#### Subroutine: MP_PSEUDO_INVERSE
C###  Description:
C###    MP_PSEUDO_INVERSE performs the Moore-Penrose pseudo
C###    inverse upon an (N,equations) matrix A and returns
C###    this inverse in the (equations,N) matrix M.
C###    The two variables used_N and used_equations allow
C###    the used space of A to be less than the maximum
C###    dimensions.      

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'tol00.cmn'
      
      REAL*8 SIGVAL
      PARAMETER(SIGVAL=1.d-6)
!     Parameter List
      INTEGER N,equations,used_N,used_equations
      REAL*8 A(N,equations),M(equations,N)
      CHARACTER ERROR*(*)

!     Local Variables
      INTEGER i,j,rank,pr,pc,k,
     '  PSr(N),QSc(equations),index,PI,QI,QJ,PK,PJ,QK,I_M,J_M,K_M
      REAL*8 row_mult,UUT,LTL,CLUUT(equations,equations),
     '  CLLTL(equations,equations),
     '  B(equations) 
      LOGICAL CONTINUE      
      
      CALL ENTERS('MP_PSEUDO_INVERSE',*9999)

C     Permutation vectors
      DO i=1,used_N
        PSr(i)=i
      ENDDO
      DO i=1,used_equations
        QSc(i)=i
      ENDDO
C     Initialise rank
      rank=0
C     Perform Gauss elimination on A
      DO i=1,MIN(used_N,used_equations)
C       Pivoting. The pivoting strategy used here looks for a
C       significant value, with preference given to searching
C       down the pivot column. It may be better to use a
C       strategy comparing all absolute values.
        pr=i
        pc=i
        IF (ABS(A(PSr(pr),QSc(pc))).LT.SIGVAL) THEN
          CONTINUE=.TRUE.
          j=i
          DO WHILE(j.LE.used_equations.AND.CONTINUE)
            k=i
            DO WHILE(k.LE.used_N.AND.CONTINUE)
              IF (ABS(A(PSr(k),QSc(j))).GT.SIGVAL) THEN
                pr=k
                pc=j
                CONTINUE=.FALSE.
              ENDIF
              k=k+1
            ENDDO !k
            j=j+1
          ENDDO !j
        ENDIF
C       Store the pivot information
        index=PSr(i)
        PSr(i)=PSr(pr)
        PSr(pr)=index
        index=QSc(i)
        QSc(i)=QSc(pc)
        QSc(pc)=index
        PI=PSr(i)
        QI=QSc(i)
C       Perform row operations
        IF (ABS(A(PI,QI)).GT.ZERO_TOL) THEN
          DO j=i+1,used_N
            PJ=PSr(j)
            row_mult=A(PJ,QI)/A(PI,QI)
            DO k=i+1,used_equations
              QK=QSc(k)
              A(PJ,QK)=A(PJ,QK)-row_mult*A(PI,QK)
            ENDDO !k
            A(PJ,QI)=row_mult
          ENDDO !j
          rank=rank+1
        ELSE !no more row reductions possible
          GOTO 100
        ENDIF
      ENDDO !generalised LU
C     Form UUT and LTL and factorise with Cholesky
 100  DO i=1,rank
        PI=PSr(i)
        QI=QSc(i)
        DO j=i,rank
          PJ=PSR(j)
          QJ=QSc(j)
          UUT=A(PI,QJ)*A(PJ,QJ)
          IF(i.EQ.j) THEN
            LTL=1
          ELSE
            LTL=A(PJ,QI)
          ENDIF
          DO k=j+1,MIN(used_N,used_equations)
            PK=PSr(k)
            QK=QSc(k)
            UUT=UUT+A(PI,QK)*A(PJ,QK)
            LTL=LTL+A(PK,QI)*A(PK,QJ)
          ENDDO !k
          DO k=MIN(used_N,used_equations)+1,used_N
            PK=PSr(k)
            LTL=LTL+A(PK,QI)*A(PK,QJ)
          ENDDO !k
          DO k=MIN(used_N,used_equations)+1,used_equations
            QK=QSc(k)
            UUT=UUT+A(PI,QK)*A(PJ,QK)
          ENDDO !k
C         Calculate the Cholesky factorisation
          CLUUT(j,i)=UUT
          CLLTL(j,i)=LTL
          DO k=1,i-1
            CLUUT(j,i)=CLUUT(j,i)-CLUUT(i,k)*CLUUT(j,k)
            CLLTL(j,i)=CLLTL(j,i)-CLLTL(i,k)*CLLTL(j,k)
          ENDDO
          IF(j.EQ.i) THEN
            CLUUT(i,i)=SQRT(CLUUT(i,i)) 
            CLLTL(i,i)=SQRT(CLLTL(i,i)) 
          ELSE
            CALL ASSERT(CLUUT(i,i).NE.0,
     '        '>> Divide by zero',ERROR,*9999)
            CLUUT(j,i)=CLUUT(j,i)/CLUUT(i,i) 
            CLLTL(j,i)=CLLTL(j,i)/CLLTL(i,i) 
          ENDIF
        ENDDO !j
      ENDDO !i
      
C     Initialise the first column of M, the weights for nq
C      DO i=1,used_equations
C        M(i,1)=0
C      ENDDO

C     Perform forward and back substitutions and products to form M         

C     Forward substitution
      DO i=1,used_N
        DO j=1,MIN(i,rank)
          IF(i.EQ.j) THEN
            B(j)=1
          ELSE
            B(j)=A(PSr(i),QSc(j))
          ENDIF
          DO k=1,j-1
            B(j)=B(j)-CLLTL(j,k)*B(k)
          ENDDO !k
          B(j)=B(j)/CLLTL(j,j)
        ENDDO !j
        DO j=min(i,rank)+1,rank
          B(j)=0
          DO k=1,j-1
            B(j)=B(j)-CLLTL(j,k)*B(k)
          ENDDO !k
          B(j)=B(j)/CLLTL(j,j)
        ENDDO !j
C       Back substitution
        DO j=rank,1,-1
          DO k=j+1,rank
            B(j)=B(j)-CLLTL(k,j)*B(k)
          ENDDO !k
          B(j)=B(j)/CLLTL(j,j)
        ENDDO !j
C       Forward substitution
        DO j=1,rank
          DO k=1,j-1
            B(j)=B(j)-CLUUT(j,k)*B(k)
          ENDDO !k
          B(j)=B(j)/CLUUT(j,j)
        ENDDO !j
C       Back substitute
        DO j=rank,1,-1
          DO k=j+1,rank
            B(j)=B(j)-CLUUT(k,j)*B(k)
          ENDDO !k
          B(j)=B(j)/CLUUT(j,j)
        ENDDO !j
        
C       Calculate M
        DO j=1,used_equations
          I_M=PSr(i)
          J_M=QSc(j)
          M(J_M,I_M)=0
          DO k=1,MIN(j,rank)
            K_M=PSr(k)
            M(J_M,I_M)=M(J_M,I_M)+A(K_M,J_M)*B(k)
          ENDDO !k
        ENDDO !j
      ENDDO !i

      CALL EXITS('MP_PSEUDO_INVERSE')
      RETURN
 9999 CALL ERRORS('MP_PSEUDO_INVERSE',ERROR)
      CALL EXITS('MP_PSEUDO_INVERSE')
      RETURN 1      
      END      

