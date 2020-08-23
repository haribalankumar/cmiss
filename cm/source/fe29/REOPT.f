      SUBROUTINE REOPT(CONSTR,M,N,NRHS,UA,LDUA,SM,LDSM,VTA,LDVTA,B,
     '  LDB,SB,VTB,LDVTB,X,LDX,REG_PARAMETER,WORK,LWORK,ERROR,*)

C#### Subroutine: REOPT
C###  Description:
C###    Computes the regularisation parameters for the optimimum
C###  relative error (RE).
C###  Note:
C###    LWORK >= MIN(M,N)*(2 + N) + NTST*(1 + N) + N
CC JMB 13-OCT-2000

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER LDB, LDSM, LDUA, LDVTA, LDVTB, LDX, LWORK, M, N, NRHS
      REAL*8 B(LDB,*), REG_PARAMETER(*), SB(*), SM(LDSM,*), UA(LDUA,*),
     '  VTA(LDVTA,*), VTB(LDVTB,*), X(LDX,*), WORK(*)
      CHARACTER CONSTR, ERROR*(*)
!     Local Variables
      INTEGER i, j, k, l, MN
      REAL*8 AX, BX, MINRE, ONE, RE, SCALE, SSQ, ZERO
      PARAMETER (ONE = 1.0d0, ZERO = 0.0d0)
!     Functions
      REAL*8 DNRM2, FMIN, REFUN
      LOGICAL LSAME
      EXTERNAL REFUN

      CALL ENTERS('REOPT',*9999)

      ! Initialisation
      MN = MIN(M,N)
      IF( LWORK.LT.(MN*(2 + N) + NTST*(1 + N) + N) ) GOTO 9999
      IF( .NOT.LSAME(CONSTR, 'I') ) GOTO 9999

      IF( ISTABILISE.EQ.2 ) THEN
        ! TGSVD
        DO j = 1,NRHS
          CALL FOURIERCOEFFS(CONSTR, M, MN, UA, LDUA, B(1,j),
     '      WORK(1), ERROR, *9999)
          DO i = 1,MN
            WORK(i) = WORK(i)/SM(i,1)
          ENDDO
          MINRE = 1.D+6
          REG_PARAMETER(j) = ONE
          DO i = 1,MN
            CALL DGEMV('T', i, N, ONE, VTA, LDVTA, WORK(1), 1,
     '        ZERO, WORK(MN + 1), 1)
            IF( ICOUPLING.EQ.1 ) THEN
               CALL DAXPY(N, -ONE, X(1,j), 1, WORK(MN + 1), 1)
               RE = DNRM2(N, WORK(MN + 1), 1)
            ELSE
              SCALE = ZERO
              SSQ = ONE
              DO k = 1,NTST
                DO l = 1,N
                  WORK(MN + N + l) = SB(j)*VTB(j,k)* WORK(MN + l)
     '              - X(l,k)
                ENDDO
                CALL DLASSQ(N, WORK(MN + N + 1), 1, SCALE, SSQ)
              ENDDO
              RE = SCALE*DSQRT(SSQ)
            ENDIF
            IF( RE.LT.MINRE ) THEN
              MINRE = RE
              REG_PARAMETER(j) = DBLE(i)
            ENDIF
          ENDDO
        ENDDO
      ELSEIF( ISTABILISE.EQ.3 ) THEN
        ! Tikhonov
        CALL DLACPY('F', MN, N, VTA, LDVTA, WORK(MN + 1), MN)
        IF ( ICOUPLING.EQ.2 ) THEN
          CALL DLACPY('F', N, NTST, X, LDX, WORK(MN*(N + 1) + 1), N)
        ENDIF
        AX = MAX(SM(MN,1),DSQRT(LOOSE_TOL))
        BX = SM(1,1)
        DO j = 1,NRHS
          IF( ICOUPLING.EQ.1 ) THEN
            CALL DCOPY(N, X(1,j), 1, WORK(MN*(N + 1) + 1), 1)
          ELSE
            CALL DCOPY(NTST, VTB(j,1), LDVTB,
     '        WORK(MN*(N + 1) + N*NTST + 1), 1)
            CALL DSCAL(NTST, SB(j), WORK(MN*(N + 1) + N*NTST + 1), 1)
          ENDIF
          CALL FOURIERCOEFFS(CONSTR, M, MN, UA, LDUA, B(1,j), WORK(1),
     '      ERROR, *9999)
          REG_PARAMETER(j) = FMIN(AX, BX, REFUN, M, N, WORK(1), SM(1,1),
     '      WORK(MN + 1))
        ENDDO
      ENDIF

      CALL EXITS('REOPT')
      RETURN
 9999 CALL ERRORS('REOPT',ERROR)
      CALL EXITS('REOPT')
      RETURN 1
      END

C---------------------------------------------------------------------
