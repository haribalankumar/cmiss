
      SUBROUTINE CONJGRAD(A,LDA,N,X,B,D,SPARSE_A,ISC_A,ISR_A,ANORM,
     '  RESID,OMEGA,ITER,NRES,NPRECON,PRECON,FACTORISE,OUTPUTCODE,
     '  ERROR,*)

C#### Subroutine: CONJGRAD
C###  Description:
C###    CONJGRAD solves a system of linear equations using a preconditioned
C###    Conjugate Gradient scheme.
C###  Written by Stuart Norris 26/02/02

      IMPLICIT NONE
      INCLUDE 'mach00.inc'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'solver.inc'
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*)
      INTEGER LDA,N,SPARSE_A,ITER,NRES,NPRECON,PRECON,OUTPUTCODE
      REAL*8 A(*),X(*),B(*),D(*),ANORM,RESID,OMEGA
      CHARACTER ERROR*(*)
      LOGICAL FACTORISE
!     Local Variables
      POINTER E_PTR,P_PTR,Q_PTR,R_PTR,Z_PTR
      DATA E_PTR,P_PTR,Q_PTR,R_PTR,Z_PTR / 5*0 /


      CALL ALLOCATE_MEMORY(N,1,DPTYPE,E_PTR,NOINIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(N,1,DPTYPE,P_PTR,NOINIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(N,1,DPTYPE,Q_PTR,NOINIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(N,1,DPTYPE,R_PTR,NOINIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(N,1,DPTYPE,Z_PTR,NOINIT,ERROR,*9999)

      CALL CONJGRAD_SUB(A,LDA,N,X,B,D,%VAL(E_PTR),%VAL(P_PTR),
     '  %VAL(Q_PTR),%VAL(R_PTR),%VAL(Z_PTR),SPARSE_A,ISC_A,ISR_A,
     '  ANORM,RESID,OMEGA,ITER,NRES,NPRECON,PRECON,FACTORISE,
     '  OUTPUTCODE,ERROR,*9999)

      IF(Z_PTR.NE.0) CALL FREE_MEMORY(Z_PTR,ERROR,*9999)
      IF(R_PTR.NE.0) CALL FREE_MEMORY(R_PTR,ERROR,*9999)
      IF(Q_PTR.NE.0) CALL FREE_MEMORY(Q_PTR,ERROR,*9999)
      IF(P_PTR.NE.0) CALL FREE_MEMORY(P_PTR,ERROR,*9999)
      IF(E_PTR.NE.0) CALL FREE_MEMORY(E_PTR,ERROR,*9999)

      RETURN

 9999 CALL ERRORS('CONJGRAD',ERROR)
      RETURN 1
      END


      SUBROUTINE CONJGRAD_SUB(A,LDA,N,X,B,D,E,P,Q,R,Z,SPARSE_A,ISC_A,
     '  ISR_A,ANORM,RESID,OMEGA,ITER,NRES,NPRECON,PRECON,FACTORISE,
     '  OUTPUTCODE,ERROR,*)

C#### Subroutine: CONJGRAD
C###  Description:
C###    CONJGRAD solves a system of linear equations using a preconditioned
C###    Conjugate Gradient scheme.
C###  Written by Stuart Norris 26/02/02

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*)
      INTEGER LDA,N,SPARSE_A,ITER,NRES,NPRECON,PRECON,OUTPUTCODE
      REAL*8 A(*),X(N),B(N),D(N),E(N),P(N),Q(N),R(N),Z(N)
      REAL*8 ANORM,RESID,OMEGA
      CHARACTER ERROR*(*)
      LOGICAL FACTORISE
!     Local Parameters
      INTEGER QUIET
      PARAMETER(QUIET=1)
!     Local Variables
      INTEGER I,MAXIT,NZA
      REAL*8 ALPHA,BETA,RHO,RHO_OLD,TOL,BNORM,SCALED_TOL,EPS
!     Functions
      REAL*8 PDDOT,PDNRM2,DLAMCH
      EXTERNAL PDDOT,PDNRM2,PDCOPY,PDAXPY,DLAMCH


      TOL=RESID
      MAXIT=ITER
      RHO=1.0D0
      BETA=1.0D0
      BNORM=PDNRM2(N,B,1)
      EPS=DLAMCH('EPS')**2
      SCALED_TOL=TOL*BNORM
      NZA=1

      IF(SPARSE_A.NE.0) NZA=ISR_A(N+1)-1
C$OMP PARALLEL DO PRIVATE(I)
      DO I=1,N
        E(I)=0.0D0
        P(I)=0.0D0
        Q(I)=0.0D0
        R(I)=0.0D0
        Z(I)=0.0D0
      ENDDO
C$OMP END PARALLEL DO

      IF(OUTPUTCODE.GE.2) THEN
        CALL WRITE_KRYLOV_SETTINGS(N,SPARSE_A,ISR_A,ANORM,BNORM,TOL,
     '    MAXIT,NRES,NPRECON,'Conjugate Gradient Solver',ERROR,*9999)
      ENDIF

C     Check input
      IF(MAXIT.LE.0) THEN
        ERROR='>>Maxit less than 0'
        GOTO 9999
      ENDIF
      IF(TOL.LT.0.0D0) THEN
        ERROR='>>Tolerance less than 0.0'
        GOTO 9999
      ENDIF

C     If first time, factorise equations for use with precondioner
      IF(FACTORISE) THEN
        CALL PRECON_FACTOR(A,LDA,N,D,PRECON,SPARSE_A,ISC_A,ISR_A,
     '    ERROR,*9999)
      ENDIF

C     Check for null system: if we have a zero RHS we exit.
      IF(BNORM.LE.0.0D0) THEN
        ITER=0
        RESID=0.0D0
C$OMP   PARALLEL DO PRIVATE(I)
        DO I=1,N
          X(I)=0.0D0
        ENDDO
C$OMP   END PARALLEL DO
        IF(OUTPUTCODE.GE.2) THEN
          WRITE(OP_STRING,7000) ITER,RESID
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        RETURN
      ENDIF

C     Calculate the residual
      CALL CALC_RESIDUAL(A,LDA,N,X,B,R,SPARSE_A,ISC_A,ISR_A,ERROR,*9999)

C     Solve the system
      DO ITER=1,MAXIT

C       Precondition
        CALL PRECON_SOLVE(A,LDA,N,Z,R,D,E,PRECON,SPARSE_A,ISC_A,ISR_A,
     '    OMEGA,NPRECON,ERROR,*9999)

C       Calculate the search direction P
        RHO_OLD=RHO
        RHO=PDDOT(N,R,1,Z,1)

        IF(ITER.EQ.1) THEN
          CALL PDCOPY(N,Z,1,P,1)
        ELSE
          IF(ABS(RHO_OLD).LE.ABS(RHO*EPS)) THEN
            WRITE(ERROR,7100) '>>Beta = Inf (',RHO,'/',RHO_OLD,')'
            GOTO 9999
          ENDIF
          BETA=RHO/RHO_OLD
          CALL PDAXPY(N,BETA,P,1,Z,1)
          CALL PDCOPY(N,Z,1,P,1)
        ENDIF

C       Calculate the search distance
        CALL MATRIX_MULT(A,LDA,N,P,Q,SPARSE_A,ISC_A,ISR_A,ERROR,*9999)

        ALPHA=PDDOT(N,P,1,Q,1)
        IF(ABS(ALPHA).LE.ABS(RHO*EPS)) THEN
          WRITE(ERROR,7100) '>>Alpha = Inf (',RHO,'/',ALPHA,')'
          GOTO 9999
        ENDIF
        ALPHA=RHO/ALPHA

C       Update the solution. Every NRES iterations update the residual
C       estimate (otherwise it tends to wander off away from the true
C       value).
        CALL PDAXPY(N,ALPHA,P,1,X,1)
        IF(MOD(ITER,NRES).EQ.0) THEN
          CALL CALC_RESIDUAL(A,LDA,N,X,B,R,SPARSE_A,ISC_A,ISR_A,ERROR,
     '      *9999)
        ELSE
          ALPHA=-ALPHA
          CALL PDAXPY(N,ALPHA,Q,1,R,1)
        ENDIF
        RESID=PDNRM2(N,R,1)
        IF(ITER.EQ.1) SCALED_TOL=TOL*(ANORM*PDNRM2(N,X,1) + BNORM)

C       Print out solution, and check for convergance
        IF(OUTPUTCODE.GE.2) THEN
          WRITE(OP_STRING,7000) ITER,RESID
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        IF(RESID.LE.SCALED_TOL) GOTO 100

      ENDDO
      ITER=MAXIT

 100  CONTINUE

      RETURN

 7000 FORMAT(1P,I5,1x,G12.6)
 7100 FORMAT(3(A,E12.5))

 9999 CALL ERRORS('CONJGRAD_SUB',ERROR)
      RETURN 1
      END


      SUBROUTINE BICGSTAB(A,LDA,N,X,B,D,SPARSE_A,ISC_A,ISR_A,ANORM,
     '  RESID,OMEGA,ITER,NRES,NPRECON,PRECON,FACTORISE,OUTPUTCODE,
     '  ERROR,*)

C#### Subroutine: BICGSTAB
C###  Description:
C###    BiCGSTAB solves a system of linear equations using the preconditioned
C###    preconditioned Bi-Conjugate Gradient STABilised scheme.
C###  Written by Stuart Norris 28/02/02

      IMPLICIT NONE
      INCLUDE 'mach00.inc'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'solver.inc'
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*)
      INTEGER LDA,N,SPARSE_A,ITER,NRES,NPRECON,PRECON,OUTPUTCODE
      REAL*8 A(*),X(*),B(*),D(*),ANORM,RESID,OMEGA
      CHARACTER ERROR*(*)
      LOGICAL FACTORISE
!     Local Variables
      POINTER E_PTR,P_PTR,Q_PTR,R_PTR,S_PTR,T_PTR,U_PTR,V_PTR,W_PTR
      DATA E_PTR,P_PTR,Q_PTR,R_PTR,S_PTR,T_PTR,U_PTR,V_PTR,W_PTR / 9*0 /


      CALL ALLOCATE_MEMORY(N,1,DPTYPE,E_PTR,NOINIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(N,1,DPTYPE,P_PTR,NOINIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(N,1,DPTYPE,Q_PTR,NOINIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(N,1,DPTYPE,R_PTR,NOINIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(N,1,DPTYPE,S_PTR,NOINIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(N,1,DPTYPE,T_PTR,NOINIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(N,1,DPTYPE,U_PTR,NOINIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(N,1,DPTYPE,V_PTR,NOINIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(N,1,DPTYPE,W_PTR,NOINIT,ERROR,*9999)

      CALL BICGSTAB_SUB(A,LDA,N,X,B,D,%VAL(E_PTR),%VAL(P_PTR),
     '  %VAL(Q_PTR),%VAL(R_PTR),%VAL(S_PTR),%VAL(T_PTR),%VAL(U_PTR),
     '  %VAL(V_PTR),%VAL(W_PTR),SPARSE_A,ISC_A,ISR_A,ANORM,RESID,OMEGA,
     '  ITER,NRES,NPRECON,PRECON,FACTORISE,OUTPUTCODE,ERROR,*9999)

      IF(W_PTR.NE.0) CALL FREE_MEMORY(W_PTR,ERROR,*9999)
      IF(V_PTR.NE.0) CALL FREE_MEMORY(V_PTR,ERROR,*9999)
      IF(U_PTR.NE.0) CALL FREE_MEMORY(U_PTR,ERROR,*9999)
      IF(T_PTR.NE.0) CALL FREE_MEMORY(T_PTR,ERROR,*9999)
      IF(S_PTR.NE.0) CALL FREE_MEMORY(S_PTR,ERROR,*9999)
      IF(R_PTR.NE.0) CALL FREE_MEMORY(R_PTR,ERROR,*9999)
      IF(Q_PTR.NE.0) CALL FREE_MEMORY(Q_PTR,ERROR,*9999)
      IF(P_PTR.NE.0) CALL FREE_MEMORY(P_PTR,ERROR,*9999)
      IF(E_PTR.NE.0) CALL FREE_MEMORY(E_PTR,ERROR,*9999)

      RETURN

 9999 CALL ERRORS('BICGSTAB',ERROR)
      RETURN 1
      END


      SUBROUTINE BICGSTAB_SUB(A,LDA,N,X,B,D,E,P,Q,R,S,T,U,V,W,SPARSE_A,
     '  ISC_A,ISR_A,ANORM,RESID,OMEGA,ITER,NRES,NPRECON,PRECON,
     '  FACTORISE,OUTPUTCODE,ERROR,*)

C#### Subroutine: BICGSTAB_SUB
C###  Description:
C###    BiCGSTAB solves a system of linear equations using the preconditioned
C###    preconditioned Bi-Conjugate Gradient STABilised scheme.
C###  Written by Stuart Norris 28/02/02

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*)
      INTEGER LDA,N,SPARSE_A,ITER,NRES,NPRECON,PRECON,OUTPUTCODE
      REAL*8 A(*),X(N),B(N),D(N),E(N),P(N),Q(N),R(N),S(N),T(N),U(N),V(N)
      REAL*8 W(N),ANORM,RESID,OMEGA
      CHARACTER ERROR*(*)
      LOGICAL FACTORISE
!     Local Variables
      INTEGER I,MAXIT
      REAL*8 ALPHA,BETA,GAMMA,RHO,RHO_OLD,SUM1,SUM2,TOL,BNORM,SCALED_TOL
      REAL*8 EPS
!     Functions
      REAL*8 PDDOT,PDNRM2,DLAMCH
      EXTERNAL PDDOT,PDNRM2,PDCOPY,PDAXPY,DLAMCH


      TOL=RESID
      MAXIT=ITER
      RHO=1.0D0
      BETA=1.0D0
      ALPHA=1.0D0
      GAMMA=1.0D0
      BNORM=PDNRM2(N,B,1)
      EPS=DLAMCH('EPS')**2
      SCALED_TOL=TOL*BNORM
C$OMP PARALLEL DO PRIVATE(I)
      DO I=1,N
        E(I)=0.0D0
        P(I)=0.0D0
        Q(I)=0.0D0
        R(I)=0.0D0
        S(I)=0.0D0
        T(I)=0.0D0
        U(I)=0.0D0
        V(I)=0.0D0
        W(I)=0.0D0
      ENDDO
C$OMP END PARALLEL DO

      IF(OUTPUTCODE.GE.2) THEN
        CALL WRITE_KRYLOV_SETTINGS(N,SPARSE_A,ISR_A,ANORM,BNORM,TOL,
     '    MAXIT,NRES,NPRECON,'BiCGSTAB Solver',ERROR,*9999)
      ENDIF

C     Check input
      IF(MAXIT.LE.0) THEN
        ERROR='>>Maxit less than 0'
        GOTO 9999
      ENDIF
      IF(TOL.LT.0.0D0) THEN
        ERROR='>>Tolerance less than 0.0'
        GOTO 9999
      ENDIF

C     If first time, factorise equations for use with precondioner
      IF(FACTORISE) THEN
        CALL PRECON_FACTOR(A,LDA,N,D,PRECON,SPARSE_A,ISC_A,ISR_A,
     '    ERROR,*9999)
      ENDIF

C     Check for null system: if we have a zero RHS we exit.
      IF(BNORM.LE.0.0D0) THEN
        ITER=0
        RESID=0.0D0
C$OMP   PARALLEL DO PRIVATE(I)
        DO I=1,N
          X(I)=0.0D0
        ENDDO
C$OMP   END PARALLEL DO
        IF(OUTPUTCODE.GE.2) THEN
          WRITE(OP_STRING,7000) ITER,RESID
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        RETURN
      ENDIF

C     Calculate the residual
      CALL CALC_RESIDUAL(A,LDA,N,X,B,R,SPARSE_A,ISC_A,ISR_A,ERROR,*9999)
      CALL PDCOPY(N,R,1,U,1)

C     Solve the system
      DO ITER=1,MAXIT

        RHO_OLD=RHO
        RHO=PDDOT(N,R,1,U,1)
        IF(ABS(RHO).LT.EPS) THEN
          ERROR='>>Breakdown: R and R~^T have become orthogonal'
          GOTO 9999
        ENDIF

        IF(ITER.EQ.1) THEN
          CALL PDCOPY(N,R,1,P,1)
        ELSE
          IF(ABS(RHO_OLD).LE.ABS(RHO*EPS)) THEN
            WRITE(ERROR,7100) '>>Beta = Inf (',RHO,'/',RHO_OLD,')'
            GOTO 9999
          ENDIF
          IF(ABS(GAMMA).LE.ABS(ALPHA*EPS)) THEN
            WRITE(ERROR,7100) '>>Beta = Inf (',ALPHA,'/',GAMMA,')'
            GOTO 9999
          ENDIF

          BETA=(RHO/RHO_OLD)*(ALPHA/GAMMA)
C$OMP     PARALLEL DO PRIVATE(I)
          DO I=1,N
            P(I)=R(I)+BETA*(P(I)-GAMMA*V(I))
          ENDDO
C$OMP     END PARALLEL DO
        ENDIF

C       Precondition
        CALL PRECON_SOLVE(A,LDA,N,Q,P,D,E,PRECON,SPARSE_A,ISC_A,ISR_A,
     '    OMEGA,NPRECON,ERROR,*9999)
        CALL MATRIX_MULT(A,LDA,N,Q,V,SPARSE_A,ISC_A,ISR_A,ERROR,*9999)

        SUM1=PDDOT(N,U,1,V,1)
        IF(ABS(SUM1).LE.ABS(RHO*EPS)) THEN
          WRITE(ERROR,7100) '>>Alpha = Inf (',RHO,'/',SUM1,')'
          GOTO 9999
        ENDIF
        ALPHA=RHO/SUM1

        CALL PDCOPY(N,R,1,S,1)
        CALL PDAXPY(N,-ALPHA,V,1,S,1)

        RESID=PDNRM2(N,S,1)
        IF(ITER.EQ.1) SCALED_TOL=TOL*(ANORM*PDNRM2(N,X,1) + BNORM)
        IF(RESID.LT.SCALED_TOL) THEN
          CALL PDAXPY(N,ALPHA,Q,1,X,1)
          IF(OUTPUTCODE.GE.2) THEN
            WRITE(OP_STRING,7000) ITER,RESID
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
          GOTO 100
        ENDIF

C       Precondition
        CALL PRECON_SOLVE(A,LDA,N,W,S,D,E,PRECON,SPARSE_A,ISC_A,ISR_A,
     '    OMEGA,NPRECON,ERROR,*9999)
        CALL MATRIX_MULT(A,LDA,N,W,T,SPARSE_A,ISC_A,ISR_A,ERROR,*9999)

        SUM1=PDNRM2(N,T,1)**2
        SUM2=PDDOT(N,T,1,S,1)
        IF(ABS(SUM1).LE.ABS(SUM2*EPS)) THEN
          WRITE(ERROR,7100) '>>Gamma = Inf (',SUM2,'/',SUM1,')'
          GOTO 9999
        ENDIF
        GAMMA=SUM2/SUM1

C       Update the solution. Every NRES iterations update the residual
C       estimate (otherwise it tends to wander off away from the true
C       value).
        CALL PDAXPY(N,ALPHA,Q,1,X,1)
        CALL PDAXPY(N,GAMMA,W,1,X,1)
        IF(MOD(ITER,NRES).EQ.0) THEN
          CALL CALC_RESIDUAL(A,LDA,N,X,B,R,SPARSE_A,ISC_A,ISR_A,ERROR,
     '      *9999)
        ELSE
          CALL PDCOPY(N,S,1,R,1)
          CALL PDAXPY(N,-GAMMA,T,1,R,1)
        ENDIF
        RESID=PDNRM2(N,R,1)
        IF(ITER.EQ.1) SCALED_TOL=TOL*(ANORM*PDNRM2(N,X,1) + BNORM)

C       Print out solution, and check for convergance
        IF(OUTPUTCODE.GE.2) THEN
          WRITE(OP_STRING,7000) ITER,RESID
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        IF(RESID.LE.SCALED_TOL) GOTO 100

C       Check for breakdown
        IF(ABS(GAMMA).LT.EPS) THEN
          ERROR='>>Breakdown: S and T have become orthogonal to T''.T'
          GOTO 9999
        ENDIF

      ENDDO
      ITER=MAXIT

 100  CONTINUE

      RETURN

 7000 FORMAT(1P,I5,1x,G12.6)
 7100 FORMAT(3(A,E12.5))

 9999 CALL ERRORS('BICGSTAB_SUB',ERROR)
      RETURN 1
      END


      SUBROUTINE GMRES(A,LDA,N,X,B,D,SPARSE_A,ISC_A,ISR_A,ANORM,RESID,
     '  OMEGA,ITER,NRES,NPRECON,PRECON,FACTORISE,OUTPUTCODE,ERROR,*)
C
C#### Subroutine: GMRES
C###  Description:
C###    GRMRES solves a system of linear equations using the preconditioned
C###    preconditioned Generalised Minimum Residual scheme.
C###  Written by Stuart Norris 01/03/02

      IMPLICIT NONE
      INCLUDE 'mach00.inc'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'solver.inc'
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*)
      INTEGER LDA,N,SPARSE_A,ITER,NRES,NPRECON,PRECON,OUTPUTCODE
      REAL*8 A(*),X(N),B(N),D(N),ANORM,RESID,OMEGA
      CHARACTER ERROR*(*)
      LOGICAL FACTORISE
!     Local Variables
      INTEGER NDIM1,NDIM2,NDIM3,LNRES
      POINTER E_PTR,R_PTR,V_PTR,S_PTR,T_PTR,U_PTR,Y_PTR,H_PTR
      DATA E_PTR,R_PTR,V_PTR,S_PTR,T_PTR,U_PTR,Y_PTR,H_PTR / 8*0 /


      IF(NRES.LE.0) THEN
        LNRES=10
      ELSE
        LNRES=NRES
      ENDIF
      NDIM1=N*(LNRES+1)
      NDIM2=LNRES+1
      NDIM3=LNRES*(LNRES+1)

      CALL ALLOCATE_MEMORY(N    ,1,DPTYPE,E_PTR,NOINIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(N    ,1,DPTYPE,R_PTR,NOINIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(NDIM1,1,DPTYPE,V_PTR,NOINIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(NDIM2,1,DPTYPE,S_PTR,NOINIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(NDIM2,1,DPTYPE,T_PTR,NOINIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(NDIM2,1,DPTYPE,U_PTR,NOINIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(NDIM2,1,DPTYPE,Y_PTR,NOINIT,ERROR,*9999)
      CALL ALLOCATE_MEMORY(NDIM3,1,DPTYPE,H_PTR,NOINIT,ERROR,*9999)

      CALL GMRES_SUB(A,LDA,N,X,B,D,%VAL(E_PTR),%VAL(R_PTR),%VAL(V_PTR),
     '  %VAL(S_PTR),%VAL(T_PTR),%VAL(U_PTR),%VAL(Y_PTR),%VAL(H_PTR),
     '  SPARSE_A,ISC_A,ISR_A,ANORM,RESID,OMEGA,ITER,LNRES,NPRECON,
     '  PRECON,FACTORISE,OUTPUTCODE,ERROR,*9999)

      IF(H_PTR.NE.0) CALL FREE_MEMORY(H_PTR,ERROR,*9999)
      IF(Y_PTR.NE.0) CALL FREE_MEMORY(Y_PTR,ERROR,*9999)
      IF(U_PTR.NE.0) CALL FREE_MEMORY(U_PTR,ERROR,*9999)
      IF(T_PTR.NE.0) CALL FREE_MEMORY(T_PTR,ERROR,*9999)
      IF(S_PTR.NE.0) CALL FREE_MEMORY(S_PTR,ERROR,*9999)
      IF(V_PTR.NE.0) CALL FREE_MEMORY(V_PTR,ERROR,*9999)
      IF(R_PTR.NE.0) CALL FREE_MEMORY(R_PTR,ERROR,*9999)
      IF(E_PTR.NE.0) CALL FREE_MEMORY(E_PTR,ERROR,*9999)

      RETURN

 9999 CALL ERRORS('GMRES',ERROR)
      RETURN 1
      END


      SUBROUTINE GMRES_SUB(A,LDA,N,X,B,D,E,R,V,S,T,U,Y,H,SPARSE_A,
     '  ISC_A,ISR_A,ANORM,RESID,OMEGA,ITER,NRES,NPRECON,PRECON,
     '  FACTORISE,OUTPUTCODE,ERROR,*)
C
C#### Subroutine: GMRES_SUB
C###  Description:
C###    GRMRES solves a system of linear equations using the preconditioned
C###    preconditioned Generalised Minimum Residual scheme.
C###  Written by Stuart Norris 01/03/02

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER ISC_A(*),ISR_A(*)
      INTEGER LDA,N,SPARSE_A,ITER,NRES,NPRECON,PRECON,OUTPUTCODE
      REAL*8 A(*),X(N),B(N),D(N),E(N),R(N),V(N,NRES+1)
      REAL*8 S(NRES+1),T(NRES+1),U(NRES+1),Y(NRES+1),H(NRES+1,NRES)
      REAL*8 ANORM,RESID,OMEGA
      CHARACTER ERROR*(*)
      LOGICAL FACTORISE
!     Local Variables
      INTEGER I,J,L,MAXIT
      REAL*8 TOL,ABSV1,ALPHA,BETA,SIGMA,TMP,BNORM,SCALED_TOL,EPS
!     Functions
      REAL*8 PDDOT,PDNRM2,DLAMCH
      EXTERNAL PDDOT,PDNRM2,PDCOPY,PDAXPY,DLAMCH


      TOL=RESID
      MAXIT=ITER
      BNORM=PDNRM2(N,B,1)
      EPS=DLAMCH('EPS')**2
      SCALED_TOL=TOL*BNORM
      SIGMA=1.0D0
C$OMP PARALLEL DO PRIVATE(I)
      DO I=1,N
        E(I)=0.0D0
        R(I)=0.0D0
      ENDDO
C$OMP END PARALLEL DO
      DO J=1,NRES+1
C$OMP   PARALLEL DO PRIVATE(I)
        DO I=1,N
          V(I,J)=0.0D0
        ENDDO
C$OMP   END PARALLEL DO
      ENDDO

      IF(OUTPUTCODE.GE.2) THEN
        CALL WRITE_KRYLOV_SETTINGS(N,SPARSE_A,ISR_A,ANORM,BNORM,TOL,
     '    MAXIT,NRES,NPRECON,'GMRES Solver',ERROR,*9999)
      ENDIF

C     Check input
      IF(MAXIT.LE.0) THEN
        ERROR='>>Maxit less than 0'
        GOTO 9999
      ENDIF
      IF(TOL.LT.0.0D0) THEN
        ERROR='>>Tolerance less than 0.0'
        GOTO 9999
      ENDIF

C     If first time, factorise equations for use with precondioner
      IF(FACTORISE) THEN
        CALL PRECON_FACTOR(A,LDA,N,D,PRECON,SPARSE_A,ISC_A,ISR_A,
     '    ERROR,*9999)
      ENDIF

C     Check for null system: if we have a zero RHS we exit.
      IF(BNORM.LE.0.0D0) THEN
        ITER=0
        RESID=0.0D0
C$OMP   PARALLEL DO PRIVATE(I)
        DO I=1,N
          X(I)=0.0D0
        ENDDO
C$OMP   END PARALLEL DO
        IF(OUTPUTCODE.GE.2) THEN
          WRITE(OP_STRING,7000) ITER,RESID
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        RETURN
      ENDIF

C     Solve the system
      ITER=1
      DO WHILE(ITER.LE.MAXIT)

C       Calculate the residual
        CALL CALC_RESIDUAL(A,LDA,N,X,B,R,SPARSE_A,ISC_A,ISR_A,ERROR,
     '    *9999)
        RESID=PDNRM2(N,R,1)

C       Precondition
        CALL PRECON_SOLVE(A,LDA,N,V(1,1),R,D,E,PRECON,SPARSE_A,ISC_A,
     '    ISR_A,OMEGA,NPRECON,ERROR,*9999)

        ABSV1=PDNRM2(N,V(1,1),1)
        IF(ABS(ABSV1).LE.ABS(RESID*EPS)) THEN
          WRITE(ERROR,7100) '>>Sigma = Inf (',RESID,'/',ABSV1,')'
          GOTO 9999
        ENDIF
        IF(ITER.EQ.1) THEN
          SIGMA=RESID/ABSV1
          SCALED_TOL=TOL*(ANORM*PDNRM2(N,X,1) + BNORM)
        ENDIF

C       Print out solution, and check for convergance
        IF(RESID.LT.SCALED_TOL) THEN
          IF(OUTPUTCODE.GE.2) THEN
            WRITE(OP_STRING,7000) ITER,RESID
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
          GOTO 200
        ENDIF

        S(1)=ABSV1
        TMP=1.0D0/ABSV1
C       DO K=1,N
C         V(K,1)=V(K,1)*TMP
C       ENDDO
        CALL PDSCAL(N,TMP,V,1)

C       Inner loop
        DO I=1,NRES

          CALL MATRIX_MULT(A,LDA,N,V(1,I),R,SPARSE_A,ISC_A,ISR_A,ERROR,
     '      *9999)

C         Precondition
          CALL PRECON_SOLVE(A,LDA,N,V(1,I+1),R,D,E,PRECON,SPARSE_A,
     '      ISC_A,ISR_A,OMEGA,NPRECON,ERROR,*9999)

          DO J=1,I
            H(J,I)=PDDOT(N,V(1,I+1),1,V(1,J),1)
            CALL PDAXPY(N,-H(J,I),V(1,J),1,V(1,I+1),1)
          ENDDO
          H(I+1,I)=PDNRM2(N,V(1,I+1),1)
          IF(ABS(H(I+1,I)).LE.EPS) THEN
            WRITE(ERROR,7100) '>>H_i+1,i = 0 (',H(I+1,I),')'
            GOTO 9999
          ENDIF
          TMP=1.0D0/H(I+1,I)
C         DO K=1,N
C           V(K,I+1)=V(K,I+1)*TMP
C         ENDDO
          CALL PDSCAL(N,TMP,V(1,I+1),1)

          DO J=1,I-1
            ALPHA=T(J)*H(  J,I)+U(J)*H(J+1,I)
            BETA =T(J)*H(J+1,I)-U(J)*H(  J,I)
            H(J,I)=ALPHA
            H(J+1,I)=BETA
          ENDDO
          TMP=SQRT(H(I,I)**2+H(I+1,I)**2)
          IF(ABS(TMP).LE.EPS) THEN
            WRITE(ERROR,7100) '>>H_i,i^2 + H_i+1,i^2 = 0 (',TMP,')'
            GOTO 9999
          ENDIF
          T(I)=H(I,I)/TMP
          U(I)=H(I+1,I)/TMP
          H(I,I)=TMP
          H(I+1,I)=0.0D0

          S(I+1)=-(U(I)*S(I))
          S(I)=S(I)*T(I)
          RESID=SIGMA*ABS(S(I+1))

C         Print out solution, and check for convergance
          IF(OUTPUTCODE.GE.2) THEN
            WRITE(OP_STRING,7000) ITER,RESID
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
          IF(RESID.LT.SCALED_TOL.OR.ITER.GE.MAXIT) GOTO 100

          ITER=ITER+1
        ENDDO
        I=I-1
 100    CONTINUE

        DO J=I,1,-1
          IF(ABS(H(J,J)).LE.ABS(S(J)*EPS)) THEN
            WRITE(ERROR,7100) '>>Y(J) = Inf (',S(J),'/',H(J,J),')'
            GOTO 9999
          ENDIF
          Y(J)=S(J)/H(J,J)
          DO L=J-1,1,-1
            S(L)=S(L)-Y(J)*H(L,J)
          ENDDO
        ENDDO

        DO J=I,1,-1
          CALL PDAXPY(N,Y(J),V(1,J),1,X,1)
        ENDDO
        IF(ITER.EQ.1) SCALED_TOL=TOL*(ANORM*PDNRM2(N,X,1) + BNORM)
        IF(RESID.LT.SCALED_TOL.OR.ITER.GE.MAXIT) GOTO 200

      ENDDO
      ITER=MAXIT

 200  CONTINUE

      RETURN

 7000 FORMAT(1P,I5,1x,G12.6)
 7100 FORMAT(3(A,E12.5))

 9999 CALL ERRORS('GMRES_SUB',ERROR)
      RETURN 1
      END


      SUBROUTINE WRITE_KRYLOV_SETTINGS(N,SPARSE_A,ISR_A,ANORM,BNORM,TOL,
     '  MAXIT,NRES,NPRECON,NAME,ERROR,*)

C#### Subroutine: WRITE_KRYLOV_SETTINGS
C###  Description:
C###    Write out the solver settings for the Krylov space solvers.
C###  Written by Stuart Norris 07/08/02

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER N,SPARSE_A,ISR_A(*),MAXIT,NRES,NPRECON
      REAL*8 ANORM,BNORM,TOL
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

      WRITE(OP_STRING,'(2A)') ' ',NAME(:LEN_TRIM(NAME))
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'()')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Solution Tolerence :'',1P,E11.4)') TOL
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Maximum Iterations :'',I11)') MAXIT
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Restart Frequency  :'',I11)') NRES
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Precon. Iterations :'',I11)') NPRECON
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'()')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      CALL WRITE_MATRIX_DIMS(SPARSE_A,N,NZA,ANORM,BNORM,ERROR,*9999)

      RETURN

 9999 CALL ERRORS('WRITE_KRYLOV_SETTINGS',ERROR)
      RETURN 1
      END
