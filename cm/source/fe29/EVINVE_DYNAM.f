      SUBROUTINE EVINVE_DYNAM(ISEG,ISIZE_TBH,ISPLOTXY,PHI,PHI_H,
     '  PHI_H_EXACT,REG_PARAMETER,SIGMA_PHI,SM_T_BH,T_BH,U_PHI,U_T_BH,
     '  VT_PHI,VT_T_BH,WORK,WORK_PHI,WORK_T_BH,CSEG,ERROR,*)

C#### Subroutine: EVINVE_DYNAM
C###  Description:
C###    EVINVE_DYNAM produces the inverse signals PHI_H.
C**** The inverse solution can be found by eiter solving the equations
C**** in a least square sense, Tikhonov regularisation or Trucated SVD
C**** regularisation.
CC JMB 13-OCT-2000

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER ISEG(*),ISIZE_TBH(*),ISPLOTXY(*)
      REAL*8 PHI(NY_TRANSFER_M,*), PHI_H(NY_TRANSFER_M,*),
     '  PHI_H_EXACT(NY_TRANSFER_M,*), REG_PARAMETER(0:NTSM),
     '  SIGMA_PHI(*), SM_T_BH(NY_TRANSFER_M,*),
     '  T_BH(NY_TRANSFER_M,NY_TRANSFER_M),
     '  U_PHI(NY_TRANSFER_M,*), U_T_BH(NY_TRANSFER_M,*), VT_PHI(NTSM,*),
     '  VT_T_BH(NY_TRANSFER_M,*), WORK(NY_TRANSFER_M,*),
     '  WORK_PHI(NY_TRANSFER_M,*), WORK_T_BH(NY_TRANSFER_M,*)
      CHARACTER CSEG(*)*(*), ERROR*(*)
!     Local Variables
      INTEGER i, IFAIL, j, K, LWORK, M, MT, N, NRHS
      INTEGER*4 IWORK_PTR
      REAL*8 ONE, ZERO
      CHARACTER CONSTR*1
      PARAMETER (K = 3, ONE = 1.0d0, ZERO = 0.0d0)

      CALL ENTERS('EVINVE_DYNAM',*9999)

      ! Initialisation
      M = ISIZE_TBH(1)
      N = ISIZE_TBH(2)
      LWORK = NY_TRANSFER_M*NY_TRANSFER_M

      ! T_BH gets overwritten by DGELS and DGESVD
      CALL DLACPY('F', M, N, T_BH, NY_TRANSFER_M, WORK_T_BH,
     '  NY_TRANSFER_M)

      IF( ISTABILISE.EQ.1 ) THEN
        ! No regularisation
        CALL DLACPY('F', M, NTST, PHI, NY_TRANSFER_M, WORK_PHI,
     '    NY_TRANSFER_M)
        IFAIL = 1
        CALL DGELS('N', M, N, NTST, WORK_T_BH, NY_TRANSFER_M,
     '    WORK_PHI, NY_TRANSFER_M, WORK, LWORK, IFAIL)
        IF( IFAIL.NE.0 ) THEN
          WRITE(OP_STRING,*)' IFAIL=',IFAIL,' in LS for T_BH'
          CALL WRITES(IOOP, OP_STRING ,ERROR, *9999)
          ERROR = '>>IFAIL<>0 in DGELS'
          GOTO 9999
        ENDIF
        CALL DLACPY('F', N, NTST, WORK_PHI, NY_TRANSFER_M, PHI_H,
     '    NY_TRANSFER_M)
      ELSE
        IF( ICONSTRAINT.EQ.1 ) THEN
          ! Standard CSVD. Zero-order Tikhonov or TSVD
          CONSTR = 'I'
          IFAIL = 1
          CALL DGESVD('S', 'S', M, N, WORK_T_BH, NY_TRANSFER_M,
     '      SM_T_BH(1,1), U_T_BH, NY_TRANSFER_M, VT_T_BH, NY_TRANSFER_M,
     '      WORK, LWORK, IFAIL)
          IF( IFAIL.NE.0 ) THEN
            WRITE(OP_STRING,*)' IFAIL=',IFAIL,' in SVD for T_BH'
            CALL WRITES(IOOP, OP_STRING, ERROR, *9999)
            ERROR = '>>IFAIL<>0 in DGESVD'
            GOTO 9999
          ENDIF

        ELSE
          ! Standard CGSVD. First or second-order Tikhonov or TGSVD
        ENDIF

        IF( ICOUPLING.EQ.2 ) THEN
          ! Greensite
          CALL DLACPY('F', M, NTST, PHI, NY_TRANSFER_M, WORK_PHI,
     '      NY_TRANSFER_M)
          IFAIL = 1
          CALL DGESVD('S', 'S', M, NTST, WORK_PHI, NY_TRANSFER_M,
     '      SIGMA_PHI, U_PHI, NY_TRANSFER_M, VT_PHI, NTSM, WORK, LWORK,
     '      IFAIL)
          IF( IFAIL.NE.0 ) THEN
            WRITE(OP_STRING,*)' IFAIL=',IFAIL,' in SVD for PHI'
            CALL WRITES(IOOP, OP_STRING, ERROR, *9999)
            ERROR = '>>IFAIL<>0 in DGESVD'
            GOTO 9999
          ENDIF
          IF( SVD_CUTOFF_RATIO.GT.ZERO ) THEN
            ! Manual inter-equation truncation rank
            NRHS = INT(SVD_CUTOFF_RATIO)
          ELSE
            MT = MIN(M,NTST)
            IWORK_PTR = 0
            CALL ALLOCATE_MEMORY(MT + K + 1, 0, DPTYPE, IWORK_PTR,
     '        MEM_INIT, ERROR, *9999)
            CALL INTEREQRANK(MT, K, NRHS, SIGMA_PHI, %VAL(IWORK_PTR),
     '        WORK, LWORK, ERROR, *9999)
            CALL FREE_MEMORY(IWORK_PTR, ERROR, *9999)
          ENDIF
          CALL DLACPY('F', M, NRHS, U_PHI, NY_TRANSFER_M, WORK_PHI,
     '      NY_TRANSFER_M)
        ELSE
          NRHS = NTST
          CALL DLACPY('F', M, NRHS, PHI, NY_TRANSFER_M, WORK_PHI,
     '      NY_TRANSFER_M)
        ENDIF

        ! Modularity
        REG_PARAMETER(0) = DBLE(NRHS)

        IF( .NOT.EVALUATE_INVERSE ) THEN
          IF( IREGULARISE.EQ.1 ) THEN
            ! GCV
            CALL GCV(CONSTR, M, N, NRHS, U_T_BH, NY_TRANSFER_M,
     '        SM_T_BH, NY_TRANSFER_M, WORK_PHI, NY_TRANSFER_M,
     '        REG_PARAMETER(1), WORK, LWORK, ERROR, *9999)
          ELSE IF( IREGULARISE.EQ.2 ) THEN
            ! L-curve
            CALL LCURVE(CONSTR, M, N, NRHS, U_T_BH, NY_TRANSFER_M,
     '        SM_T_BH, NY_TRANSFER_M, WORK_PHI, NY_TRANSFER_M,
     '        REG_PARAMETER(1), WORK, LWORK, ERROR, *9999)
          ELSE IF( IREGULARISE.EQ.3 ) THEN
            ! Picard
            CALL PICARD(CONSTR, M, N, NRHS, U_T_BH, NY_TRANSFER_M,
     '        SM_T_BH, NY_TRANSFER_M, WORK_PHI, NY_TRANSFER_M,
     '        REG_PARAMETER(1), WORK, LWORK, ISEG, ISPLOTXY, CSEG,
     '        ERROR, *9999)
          ELSE IF( IREGULARISE.EQ.4 ) THEN
            ! Quasiopt
            CALL QUASIOPT(CONSTR, M, N,NRHS,  U_T_BH, NY_TRANSFER_M,
     '        SM_T_BH, NY_TRANSFER_M, WORK_PHI, NY_TRANSFER_M,
     '        REG_PARAMETER(1), WORK, LWORK, ERROR, *9999)
          ELSE IF( IREGULARISE.EQ.5 ) THEN
            ! Optimal
            CALL REOPT(CONSTR, M, N, NRHS, U_T_BH, NY_TRANSFER_M,
     '        SM_T_BH, NY_TRANSFER_M, VT_T_BH, NY_TRANSFER_M,
     '        WORK_PHI, NY_TRANSFER_M, SIGMA_PHI, VT_PHI, NTSM,
     '        PHI_H_EXACT, NY_TRANSFER_M, REG_PARAMETER(1), WORK, LWORK,
     '        ERROR, *9999)
          ELSE IF( IREGULARISE.EQ.6 ) THEN
            ! CRESO
            CALL CRESO(CONSTR, M, N, NRHS, U_T_BH, NY_TRANSFER_M,
     '        SM_T_BH, NY_TRANSFER_M, WORK_PHI, NY_TRANSFER_M,
     '        REG_PARAMETER(1), WORK, LWORK, ERROR, *9999)
          ELSE IF( IREGULARISE.EQ.7 ) THEN
            ! Zero-crossing
            CALL ZEROCROSSING(CONSTR, M, N, NRHS, U_T_BH, NY_TRANSFER_M,
     '        SM_T_BH, NY_TRANSFER_M, WORK_PHI, NY_TRANSFER_M,
     '        REG_PARAMETER(1), WORK, LWORK, ERROR, *9999)
          ENDIF
        ENDIF
        IF( ISTABILISE.EQ.2 ) THEN
          ! TGSVD pseudo inverse
          DO j = 1,NRHS
            CALL TGSVD(CONSTR, M, N, U_T_BH, NY_TRANSFER_M, SM_T_BH,
     '        NY_TRANSFER_M, VT_T_BH, NY_TRANSFER_M, PHI_H(1,j),
     '        WORK_PHI(1,j), INT(REG_PARAMETER(j)), WORK, LWORK,
     '        ERROR, *9999)
          ENDDO
        ELSE IF( ISTABILISE.EQ.3 ) THEN
          ! Tikhonov pseudo inverse
          DO j = 1,NRHS
            CALL TIKHONOV(CONSTR, M, N, U_T_BH, NY_TRANSFER_M,
     '        SM_T_BH, NY_TRANSFER_M, VT_T_BH, NY_TRANSFER_M,
     '        PHI_H(1,j), WORK_PHI(1,j), REG_PARAMETER(j), WORK, LWORK,
     '        ERROR, *9999)
          ENDDO
        ENDIF
        IF( ICOUPLING.EQ.2 ) THEN
          ! Greensite
          DO j = 1,NRHS
            DO i = 1,N
              WORK_PHI(i,j) = SIGMA_PHI(j)*PHI_H(i,j)
            ENDDO
          ENDDO
          CALL DGEMM('N', 'N', N, NTST, NRHS, ONE, WORK_PHI,
     '      NY_TRANSFER_M, VT_PHI, NTSM, ZERO, PHI_H, NY_TRANSFER_M)
        ENDIF
      ENDIF

      CALL EXITS('EVINVE_DYNAM')
      RETURN
 9999 CALL ERRORS('EVINVE_DYNAM',ERROR)
      CALL EXITS('EVINVE_DYNAM')
      RETURN 1
      END


