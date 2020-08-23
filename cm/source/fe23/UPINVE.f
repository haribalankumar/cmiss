      SUBROUTINE UPINVE(ISIZE_PHI,nBody,nHeart,NHP,NKH,NPLIST3,NVHP,
     '  NXLIST,NYNP,NYNR,LAPL,LAPLSQR,PHI,T_BH,YP,STRING,ERROR,*)

C#### Subroutine: UPINVE
C###  Description:
C###    UPINVE updates inverse of transfer matrix (with regularisation).
C###    Presently coded for activation times using an iterative method.

C     ISIZE_TBH passed into separate variables nBody and nHeart for
C     dimensioning local arrays
C     ISIZE_TBH(1) is nBody  = number of torso points
C     ISIZE_TBH(2) is nHeart = number of heart points

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'time02.cmn'
      INCLUDE 'tol00.cmn'
      INCLUDE 'trsf00.cmn'

!     Parameter List
      INTEGER ISIZE_PHI(2),nBody,nHeart,
     '  NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),NPLIST3(0:NPM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NXLIST(0:NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM)
      REAL*8 LAPL(NY_TRANSFER_M,NY_TRANSFER_M),
     '  LAPLSQR(NY_TRANSFER_M,NY_TRANSFER_M),PHI(NY_TRANSFER_M,NTSM),
     '  T_BH(NY_TRANSFER_M,NY_TRANSFER_M),YP(NYM,NIYM,NXM)
      CHARACTER STRING*(MXCH),ERROR*(*)

!     Parameters
C     For spline fitting: K=order, NSP=num spline points, NUMLC=num L-curves
      INTEGER K, NSP, NUMLC
      PARAMETER (K = 4, NSP = 2000, NUMLC = 100)
      REAL*8 R,C,TOL
      PARAMETER (R = 0.618033989d0, C = 1.0d0-R, TOL = 1.0d-3)

!     Local Variables
      INTEGER ERR,i,iB,IBEG,IEND,iH,iHH,iK,iLC,iter,IT_NUM,j,N3CO,
     '  nh,nk,nLC,no_nynr,no_nynr1,np,nr,nts,nv,ny,ny1,
     '  S_END,S_START,nx,nxc
      REAL T0(1),T1(1),T2(1)
      REAL*8 dPHIdTAU(nBody,nHeart,NTSM),dPHIdTAU2(nHeart,nHeart),
     '  d_PHIDIFF(nHeart),d_TAU(nBody,NTSM),
     '  DX,F1,F2,K1,K2,KAPPA,KAPPAMAX,L,LAMBDA,LC(NUMLC,4),
     '  LFACTOR,LNORM,LNORM_L,NORM_MIN,PHIDIFF(nBody,NTSM),
     '  RNORM,RNORM_L,SPLC(0:NSP,0:2,3),TEND,TIME,TSTART,X,X0,X1,X2,X3
      CHARACTER FILE*100
      LOGICAL OPFILE,FOUND,LOGSCALE

!     Functions
      INTEGER CALC_SAMPLE_FROM_TIME,IFROMC
      REAL*8 BVALUE,CALC_FORWARD_ACTN,CALC_LNORM,CALC_RNORM,
     '  CALC_TIME_FROM_SAMPLE,RFROMC
      LOGICAL CBBREV

      CALL ENTERS('UPINVE',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM update inverse<;FILENAME> activation
C###  Parameter:      <tstart #[0.0]>
C###    Defines the lower time limit for which the inverse is calculated.
C###  Parameter:      <tend #[end]>
C###    Defines the upper time limit for which the inverse is calculated.
C###  Parameter:      <iterate #[1]>
C###    Number of iterations to perform.
C###  Parameter:      <lambda_override #>
C###    Override the constant regularisation parameter.
C###  Parameter:      <factor_lambda #[1.0d4]>
C###    Search factor for regularisation parameter (relative to default
C###    value).
C###  Parameter:      <logscale>
C###    Defines whether l-curve corner is found in log space (default no).
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Performs one or more iterations of the iterative method for
C###    converging on an activation inverse solution. IF FILENAME is
C###    specified, L-curve values are written to FILENAME.opinve.
C###    Questions in the ipinve file control the behaviour of the
C###    updates.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<tstart #[0.0]>'
        OP_STRING(3)=BLANK(1:15)//'<tend #[end]>'
        OP_STRING(4)=BLANK(1:15)//'<iterate #[1]>'
        OP_STRING(5)=BLANK(1:15)//'<lambda_override #>'
        OP_STRING(6)=BLANK(1:15)//'<factor_lambda #[1.0d4]>'
        OP_STRING(7)=BLANK(1:15)//'<logscale>'
        OP_STRING(8)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc','UPINVE',ERROR,*9999)
      ELSE

        IF(NTCOQU(noco).GT.0) THEN
          OPFILE = .TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI = IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opinve','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE = .FALSE.
          IOFI = IOOP
        ENDIF

        CALL ASSERT(EVALUATE_TRANSFER,'>>Evaluate transfer first',
     '    ERROR,*9999)

        IF(CBBREV(CO,'TSTART',2,noco+1,NTCO,N3CO)) THEN
          TSTART = RFROMC(CO(N3CO+1))
        ELSE
          TSTART = 0.0d0
        ENDIF
        IF(CBBREV(CO,'TEND',2,noco+1,NTCO,N3CO)) THEN
          TEND = RFROMC(CO(N3CO+1))
        ELSE
          TEND = CALC_TIME_FROM_SAMPLE(ISIZE_PHI(2))
        ENDIF
        ERR=0
        S_START = CALC_SAMPLE_FROM_TIME(TSTART,ERR,ERROR)
          IF(ERR.NE.0) GOTO 9999
        S_END = CALC_SAMPLE_FROM_TIME(TEND,ERR,ERROR)
          IF(ERR.NE.0) GOTO 9999
        IF(CBBREV(CO,'ITERATE',1,noco+1,NTCO,N3CO)) THEN
          IT_NUM = IFROMC(CO(N3CO+1))
        ELSE
          IT_NUM = 1
        ENDIF

        IF(CBBREV(CO,'LAMBDA_OVERRIDE',2,noco+1,NTCO,N3CO)) THEN
          REG_PARAM_LAPLACE = RFROMC(CO(N3CO+1))
        ENDIF
        IF(CBBREV(CO,'FACTOR_LAMBDA',1,noco+1,NTCO,N3CO)) THEN
          LFACTOR = RFROMC(CO(N3CO+1))
        ELSE
          LFACTOR = 1.0d4 ! initial search factor for Lambda
        ENDIF
        IF(CBBREV(CO,'LOGSCALE',2,noco+1,NTCO,N3CO)) THEN
          LOGSCALE = .TRUE.
        ELSE
          LOGSCALE = .FALSE.
        ENDIF

        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc = NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        ERR = 0

        WRITE(OP_STRING,'(A6,I5,A6,I5)') 'Body',nBody,'Heart',nHeart
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

        nr = TRSF_NR_FIRST  ! transfer region for first surface (heart)
        nh = 1
        nv = 1
        nk = 1
        WRITE(OP_STRING(1),'(A48,A24)') 'Results','Timings'
        WRITE(OP_STRING(2),'(4A12,3A8)') 'Iteration','|Res|',
     '    '|Sol|','Lambda','Setup','Lambda','Total'
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

        LAMBDA = REG_PARAM_LAPLACE  !default value

        DO iter = 1,IT_NUM
          CALL CPU_TIMER(CPU_USER,T0(1))

C***      Zero work arrays
          DO iB = 1,nBody
            DO nts = S_START,S_END
              PHIDIFF(iB,nts) = 0.0d0       ! stores (Phi-Phi^)|k
              d_TAU(iB,nts) = 0.0d0         ! stores (dPhi/dTau)^T.Tau|k
            ENDDO
            DO iH = 1,nHeart
              DO nts = S_START,S_END
                dPHIdTAU(iB,iH,nts) = 0.0d0 ! stores (dPhi/dTau)
              ENDDO
            ENDDO
          ENDDO
          DO iH = 1,nHeart
            d_PHIDIFF(iH) = 0.0d0           ! stores (dPhi/dTau)^T.(Phi-Phi^)
            DO iHH = 1,nHeart
              dPHIdTAU2(iH,iHH) = 0.0d0     ! stores (dPhi/dTau)^T.(dPhi/dTau)
            ENDDO
          ENDDO

C***      Store initial activation times
C***      -       tau  in YP(ny,8,nx)
C***      -  L^tL(tau) in YP(ny,5,nx)
C***      If surface Laplacian regularisation
C***        L is the Laplacian -> premultiply initial activation times
C***      otherwise
C***        L is the identity -> YP(ny,5,nx) is the original activation time
          DO no_nynr = 1,NYNR(0,0,1,nr,nx)
            ny = NYNR(no_nynr,0,1,nr,nx)
            YP(ny,8,nx) = YP(ny,1,nx)
            IF(IREGULARISE.EQ.1) THEN       !no regularisation
              YP(ny,5,nx) = YP(ny,1,nx)
            ELSE IF(IREGULARISE.EQ.2) THEN  !surface Laplacian
              YP(ny,5,nx) = 0.0d0
              DO no_nynr1 = 1,NYNR(0,0,1,nr,nx)
                ny1 = NYNR(no_nynr1,0,1,nr,nx)
                YP(ny,5,nx) = YP(ny,5,nx) +
     '            LAPLSQR(no_nynr,no_nynr1)*YP(ny1,1,nx)
              ENDDO
            ENDIF !IREGULARISE
          ENDDO
C***      Compute solution norm = norm of L(tau)
          LNORM = CALC_LNORM(NYNR(0,0,1,nr,nx),LAPL,YP(1,1,nx))

          DO nts = S_START,S_END
            TIME = CALC_TIME_FROM_SAMPLE(nts)
C***        Compute Phi-Phi^
            DO iB = 1,nBody
              PHIDIFF(iB,nts) = PHI(iB,nts) -
     '          CALC_FORWARD_ACTN(0,0,0,0,iB,NHP(1,0,nx),NKH,NPLIST3,
     '          NVHP,nx,NYNP,0.5d0/TRSF_FREQUENCY,T_BH,TIME,YP(1,1,nx),
     '          TRSF_ACTN_WAVE_INTERPOLATE,ERR,ERROR)
                IF (ERR.NE.0) GOTO 9999
C***          Compute dPhi/dTau
              DO iH = 1,nHeart
                np = NPLIST3(iH)
                ny = NYNP(nk,nv,nh,np,0,1,nr)
                dPHIdTAU(iB,iH,nts) = CALC_FORWARD_ACTN(2,0,ny,iH,iB,
     '            NHP(1,0,nx),NKH,NPLIST3,NVHP,nx,NYNP,
     '            0.5d0/TRSF_FREQUENCY,T_BH,TIME,YP(1,1,nx),
     '            TRSF_ACTN_WAVE_INTERPOLATE,ERR,ERROR)
                  IF (ERR.NE.0) GOTO 9999
                d_TAU(iB,nts) = d_TAU(iB,nts) +
     '            dPHIdTAU(iB,iH,nts)*YP(ny,1,nx)
              ENDDO
            ENDDO

C***        Sum over time (dPhi/dTau)^T.(Phi-Phi_h) -> d_PHIDIFF
            CALL DGEMV('T',nBody,nHeart,1.0d0,dPHIdTAU(1,1,nts),
     '        nBody,PHIDIFF(1,nts),1,1.0d0,d_PHIDIFF,1)
C***        Sum over time (dPhi/dTau)^T.(dPhi/dTau) -> dPHIdTAU2
            CALL DGEMM('T','N',nHeart,nHeart,nBody,1.0d0,
     '        dPHIdTAU(1,1,nts),nBody,dPHIdTAU(1,1,nts),
     '        nBody,1.0d0,dPHIdTAU2,nHeart)
          ENDDO !nts

C***      Compute residual norm

C LKC 10-JAN-2006 Seems T_BH isn't supposed to be passed in
C          RNORM = CALC_RNORM(nBody,nHeart,NYNR(0,0,1,nr,nx),
C     '      S_START,S_END,dPHIdTAU,d_TAU,PHIDIFF,T_BH,YP(1,1,nx))
          RNORM = CALC_RNORM(nBody,nHeart,NYNR(0,0,1,nr,nx),
     '      S_START,S_END,dPHIdTAU,d_TAU,PHIDIFF,YP(1,1,nx))
          CALL CPU_TIMER(CPU_USER,T1(1))

          IF(IREGULARISE.EQ.1) THEN       !no regularisation
            LAMBDA = 0.0d0

          ELSE IF(IREGULARISE.EQ.2) THEN  !surface Laplacian regularisation
            IF (ICONSTRAINT.EQ.1) THEN    !Constant Lambda
              LAMBDA = REG_PARAM_LAPLACE
              IF(OPFILE) THEN
C***            Compute an L-Curve
                WRITE(IOFI,'(/3A15)') 'L','Lnorm','Rnorm'
                L = LAMBDA/LFACTOR
                DO WHILE (L.LT.LAMBDA*LFACTOR)
                  CALL CALC_INVERSE_SOLN(nBody,nHeart,NYNR(0,0,1,nr,nx),
     '              S_START,S_END,dPHIdTAU,dPHIdTAU2,d_PHIDIFF,d_TAU,L,
     '              LAPL,LAPLSQR,LNORM_L,PHIDIFF,RNORM_L,
     '              YP(1,1,nx),OPFILE,ERROR,*9999)
                  L = L*SQRT(SQRT(10.0d0))
                ENDDO
              ENDIF

            ELSE IF(ICONSTRAINT.EQ.2) THEN  !Ratio
              IF(LNORM.GT.ZERO_TOL) THEN
                LAMBDA = RNORM/LNORM
              ELSE
                LAMBDA = REG_PARAM_LAPLACE
              ENDIF

            ELSE IF(ICONSTRAINT.EQ.3) THEN !L-curve zero
C***        Choose Lambda to balance the solution and residual norms
C***        => looking for
C***            Rnorm - L*Lnorm = 0
              IF(OPFILE) THEN
                WRITE(IOFI,'(/3A15)') 'L','Lnorm','Rnorm'
              ENDIF
C             Bracket zero
              FOUND = .FALSE.
              L = LAMBDA/LFACTOR
              DO WHILE (L.LT.LAMBDA*LFACTOR.AND..NOT.FOUND)
                CALL CALC_INVERSE_SOLN(nBody,nHeart,NYNR(0,0,1,nr,nx),
     '            S_START,S_END,dPHIdTAU,dPHIdTAU2,d_PHIDIFF,d_TAU,L,
     '            LAPL,LAPLSQR,LNORM_L,PHIDIFF,RNORM_L,YP(1,1,nx),
     '            OPFILE,ERROR,*9999)
                IF(.NOT.FOUND.AND.RNORM_L.LT.(L*LNORM_L)) THEN
                  LAMBDA = L
                  FOUND  = .TRUE.
                ENDIF
                L = L*10.0d0
              ENDDO
              IF(FOUND) THEN  !bisection search for zero
                X2 = LAMBDA
                X1 = LAMBDA/10.0d0
                DX = X2 - X1
                F2 = RNORM_L - LAMBDA*LNORM_L
                DO WHILE(DX.GT.TOL*ABS(X2).AND.ABS(F2).GT.ZERO_TOL)
                  DX = DX*0.5d0
                  X2 = X1 + DX
                  LAMBDA = X2
                  CALL CALC_INVERSE_SOLN(nBody,nHeart,NYNR(0,0,1,nr,nx),
     '              S_START,S_END,dPHIdTAU,dPHIdTAU2,d_PHIDIFF,d_TAU,
     '              LAMBDA,LAPL,LAPLSQR,LNORM_L,PHIDIFF,RNORM_L,
     '              YP(1,1,nx),OPFILE,ERROR,*9999)
                  F2 = RNORM_L - LAMBDA*LNORM_L
                  IF(F2.GT.0.0d0) X1 = X2
                ENDDO
                LFACTOR = 1.0d2 ! reduce search for lambda
              ELSE  !if no zero, then use default value
                LAMBDA = REG_PARAM_LAPLACE
              ENDIF

            ELSE IF( ICONSTRAINT.EQ.4 ) THEN  !Maximum curvature of L-curve
C***          Compute an L-Curve
              L = LAMBDA/LFACTOR
              nLC = 0
              DO WHILE ( L.LT.LAMBDA*LFACTOR .AND. nLC.LT.NUMLC )
                nLC = nLC+1
                CALL CALC_INVERSE_SOLN(nBody,nHeart,NYNR(0,0,1,nr,nx),
     '            S_START,S_END,dPHIdTAU,dPHIdTAU2,d_PHIDIFF,d_TAU,L,
     '            LAPL,LAPLSQR,LNORM_L,PHIDIFF,RNORM_L,YP(1,1,nx),
     '            .FALSE.,ERROR,*9999)
                IF( LOGSCALE ) THEN
                  LC(nLC,1) = DLOG(RNORM_L)
                  LC(nLC,2) = DLOG(LNORM_L)
                ELSE
                  LC(nLC,1) = RNORM_L
                  LC(nLC,2) = LNORM_L
                ENDIF
                LC(nLC,3) = L
                L = L * SQRT(SQRT(10.0d0))
              ENDDO
              ! Fit a 2-D spline curve to the L-curve
              DO i = 1,nLC+K
                LC(i,4) = DBLE(i)
              ENDDO
              ! Evaluate spline between K .. nLC+1
              DO i = 0,NSP
                 X = DBLE(K) + DBLE(i) * DBLE(nLC+1-K)/DBLE(NSP)
                DO j = 0,2     !Value, first and second derivative
                  DO iLC = 1,3 !splines for Rnorm, Lnorm and Lambda
                    SPLC(i,j,iLC) = BVALUE(LC(1,4),LC(1,iLC),nLC,K,X,j)
                  ENDDO
                ENDDO
              ENDDO
              ! Compute the corner of the discretised spline
              KAPPAMAX = 0.0d0
              iK = 0
              IF( OPFILE ) THEN
                ! Write out values and curvature
                WRITE(IOFI,'(/4A15)') 'L', 'Lnorm', 'Rnorm', 'Kappa'
              ENDIF
              DO i = 0,NSP
                K1 = SPLC(i,1,1)*SPLC(i,2,2) - SPLC(i,2,1)*SPLC(i,1,2)
                K2 = (SPLC(i,1,1)**2.0d0 + SPLC(i,1,2)**2.0d0)**(1.5d0)
                IF( DABS(K2).GT.ZERO_TOL ) THEN
                  KAPPA = K1/K2
                  IF( KAPPA.GT.KAPPAMAX ) THEN
                    KAPPAMAX = KAPPA
                    iK = i
                  ENDIF
                  IF( OPFILE .AND. MOD(i,10).EQ.0 ) THEN
                    WRITE(IOFI,'(4F15.8)') SPLC(i,0,3), SPLC(i,0,2),
     '                SPLC(i,0,1), KAPPA
                  ENDIF
                ENDIF
              ENDDO
              IF( KAPPAMAX.GT.ZERO_TOL ) THEN
                LAMBDA = SPLC(iK,0,3)
                LFACTOR = 1.0d2 !reduce search for lambda
              ELSE
                LAMBDA = REG_PARAM_LAPLACE
              ENDIF

            ELSE IF( ICONSTRAINT.EQ.5 ) THEN  !Minimum product of norms
              IF(OPFILE) THEN
                WRITE(IOFI,'(/3A15)') 'L','Lnorm','Rnorm'
              ENDIF
C             Bracket minimum
              NORM_MIN = 1.0d12
              L = LAMBDA/LFACTOR
              FOUND = .FALSE.
              DO WHILE (L.LT.LAMBDA*LFACTOR.AND..NOT.FOUND)
                CALL CALC_INVERSE_SOLN(nBody,nHeart,NYNR(0,0,1,nr,nx),
     '            S_START,S_END,dPHIdTAU,dPHIdTAU2,d_PHIDIFF,d_TAU,L,
     '            LAPL,LAPLSQR,LNORM_L,PHIDIFF,RNORM_L,YP(1,1,nx),
     '            OPFILE,ERROR,*9999)
                IF(RNORM_L*LNORM_L.LT.NORM_MIN) THEN
                  LAMBDA = L
                  NORM_MIN = RNORM_L*LNORM_L
                ELSE
                  FOUND = .TRUE.
                ENDIF
                L = L*10.0d0
              ENDDO
              IF(FOUND) THEN
C               Golden section search for minimum
                X0 = LAMBDA/10.0d0
                X1 = LAMBDA
                F1 = NORM_MIN
                X3 = LAMBDA*10.0d0
                X2 = R*X1 + C*X3
                CALL CALC_INVERSE_SOLN(nBody,nHeart,NYNR(0,0,1,nr,nx),
     '            S_START,S_END,dPHIdTAU,dPHIdTAU2,d_PHIDIFF,d_TAU,X2,
     '            LAPL,LAPLSQR,LNORM_L,PHIDIFF,RNORM_L,YP(1,1,nx),
     '            OPFILE,ERROR,*9999)
                F2 = RNORM_L*LNORM_L
                DO WHILE (ABS(X3-X0).GT.TOL*(ABS(X1)+ABS(X2)))
                  IF(F2.LT.F1) THEN
                    X0 = X1
                    X1 = X2
                    X2 = R*X1 + C*X3
                    F1 = F2
                    CALL CALC_INVERSE_SOLN(nBody,nHeart,
     '                NYNR(0,0,1,nr,nx),S_START,S_END,dPHIdTAU,
     '                dPHIdTAU2,d_PHIDIFF,d_TAU,X2,LAPL,LAPLSQR,LNORM_L,
     '                PHIDIFF,RNORM_L,YP(1,1,nx),OPFILE,ERROR,*9999)
                    F2 = RNORM_L*LNORM_L
                  ELSE
                    X3 = X2
                    X2 = X1
                    X1 = R*X2 + C*X0
                    F2 = F1
                    CALL CALC_INVERSE_SOLN(nBody,nHeart,
     '                NYNR(0,0,1,nr,nx),S_START,S_END,dPHIdTAU,
     '                dPHIdTAU2,d_PHIDIFF,d_TAU,X1,LAPL,LAPLSQR,LNORM_L,
     '                PHIDIFF,RNORM_L,YP(1,1,nx),OPFILE,ERROR,*9999)
                    F1 = RNORM_L*LNORM_L
                  ENDIF
                ENDDO
                IF(F1.LT.F2) THEN
                  LAMBDA = X1
                ELSE
                  LAMBDA = X2
                ENDIF
                LFACTOR = 1.0d2 !reduce search for lambda
              ELSE ! Minimum not found
                LAMBDA = REG_PARAM_LAPLACE
              ENDIF
            ENDIF !ICONSTRAINT
          ENDIF  !IREGULARISE (?regularisation)

          CALL CALC_INVERSE_SOLN(nBody,nHeart,NYNR(0,0,1,nr,nx),
     '      S_START,S_END,dPHIdTAU,dPHIdTAU2,d_PHIDIFF,d_TAU,LAMBDA,
     '      LAPL,LAPLSQR,LNORM_L,PHIDIFF,RNORM_L,YP(1,1,nx),
     '      .FALSE.,ERROR,*9999)

          CALL CPU_TIMER(CPU_USER,T2(1))

          WRITE(OP_STRING,'(I6,6X,3E12.4,3F8.2)')
     '      iter,RNORM_L,LNORM_L,LAMBDA,
     '      T1(1)-T0(1),T2(1)-T1(1),T2(1)-T0(1)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDDO !iter

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI = IOOP
        ENDIF
      ENDIF

      CALL EXITS('UPINVE')
      RETURN
 9999 CALL ERRORS('UPINVE',ERROR)
      CALL EXITS('UPINVE')
      RETURN 1
      END


