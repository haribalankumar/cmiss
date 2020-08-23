      SUBROUTINE NAGMINA(ISTATE,PAOPTI,PMIN,PMAX,R,RESID,RESJAC,XC,
     '  IUSER,USER,OUTPTPRN,ERROR,*)

C#### Subroutine: NAGMINA
C###  Description:
C###    NAGMINA is a minimisation routine using MINLSSQP.


C**** KTYP26=1 KTYP27=2: (Sum of squared reaction differences)
C****   Starting point is PAOPTI in NTOPTI-dimensional space.
C****   Uses objective function given by subr FUNCT1 -old: in ARCHIVE
C****   Returns location of min in PAOPTI and function value in FRET.
C**** KTYP26=1 KTYP27=6: Sum of squared react. diff. unknown geom
C**** KTYP26=2 KTYP27=1: (Minimize squared area of trapezoids)
C**** KTYP26=2 KTYP27=2: (Minimize curvature of stripes)
C**** KTYP26=2 KTYP27=3: (Zero flux differences)
C****   Starting point is PAOPTI in NTOPTI-dimensional space.
C****   Uses objective function given by FUNCT3 -old: in ARCHIVE
C**** KTYP26=2 KTYP27=5: (Geometric data fitting)
C**** KTYP26=2 KTYP27=6: (Fluid interface condition)
C****   Starting point is PAOPTI in NTOPTI-dimensional space.
C****   Uses objective function given by subr FUNCT1 -old: in ARCHIVE
C****   Returns location of min in PAOPTI and function value in FRET.
C**** KTYP26=2 KTYP27=7: (Aerofoil wake pressure difference)
C****   PAOPTI contains initial y positions of wake nodes
C**** KTYP26=2 KTYP27=8: (Aerofoil lift)
C****   PAOPTI contains initial y positions of aerofoil nodes
C**** KTYP26=2 KTYP27=9: (Boundary layer thickness condition)
C**** KTYP26=2 KTYP27=10: Torso mesh parameter optimisation
C**** KTYP26=2 KTYP27=11: Difference in second moments optimisation
C**** KTYP26=2 KTYP27=12: (Inverse activation times)
C**** KTYP26=2 KTYP27=13: (Dipole Inverses)
C**** KTYP26=3 : Micro-structure, fibre stress

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'ofst00.cmn'
      INCLUDE 'opti00.cmn'
      INCLUDE 'time02.cmn'
!     Parameter List
      INTEGER ISTATE(*),IUSER(*)
      REAL*8 PAOPTI(*),PMIN(*),PMAX(*),R(NOPM,*),RESID(*),
     '  RESJAC(NREM,*),USER(*),XC(*)
      CHARACTER ERROR*(*)
      LOGICAL OUTPTPRN
!     Local Variables
      INTEGER endcol,ERR,i,j,k,LDFJ,LDR,NCNLN,noopti,nr,nxc,
     '  nx_opt,IFAIL,ITER,IEVAL
      REAL ELAPSED_TIME,TIME_START(1),TIME_STOP(1)
      REAL*8 CORREL(NOPM,NOPM),COVAR(NOPM,NOPM),OBJF,OPTION(56)
      CHARACTER STRING*(MXCH),CHAR8*10
      LOGICAL CALCCORR,END
      EXTERNAL FUNCT2

      CALL ENTERS('NAGMINA',*9999)

C CPB 19/3/94 Need to look at nx classes here

      nxc=1 ! temporary

      CALL NX_LOC(NX_INQUIRE,nxc,nx_opt,NX_OPTI,ERROR,*9999)
      CALL ASSERT(nx_opt.GT.0,'>>no nx defined for optimisation',
     '  ERROR,*9999)

      LDFJ=NREM  !1st dimension of RESJAC
      LDR=NOPM   !1st dimension of R


C!!! LKC 22-MAR-1999 Removing the loop as it does not appear appropriate
C      DO nr=1,NRT

      IFAIL=0 ! for cases when MINLSSQP is not called

C     mat params, sqd reaction diffs or geometric sqd
      IF(KTYP26.EQ.1.AND.(KTYP27.EQ.2.OR.KTYP27.EQ.5.OR.
     '  KTYP27.EQ.6)) THEN 
C        nr=1 !in future may need to put into nr loop
C        rgb 7/4/98. nr=1 now inside loop DO nr=1,NRT

C!!! LKC 22-MAR-1999 re-hard coding nr=1 and adding warning
        nr=1
        WRITE(OP_STRING,
     '    '(/'' >>WARNING: Region number is fixed to 1'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        IF(ITYP5(nr,1).EQ.2.AND.ITYP19(nr,1).EQ.1.AND.
     '    ITYP2(nr,1).EQ.9) THEN !Activation model

          CALL ASSERT(CALL_SOLV,'>>Solution mapping arrays not defined',
     '      ERROR,*9999)
          CO(1)='FEM'
           CO(2)='DEFINE'
          CO(3)='INITIAL'
          NTCO=3
          noco=1
          COQU(3,1)='W'
          COQU(3,2)='OPTIMISE'
          NTCOQU(1)=0
          NTCOQU(2)=0
          NTCOQU(3)=2
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)
          CO(1)='FEM'
           CO(2)='DEFINE'
          CO(3)='SOLVE'
          NTCO=3
          noco=1
          COQU(3,1)='W'
          COQU(3,2)='OPTIMISE'
          NTCOQU(1)=0
          NTCOQU(2)=0
          NTCOQU(3)=2
          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
     '      ERROR,*9999)
C!!! NT_RES should be set before entering NAGMIN
          NT_RES=NDT
        ENDIF

        CALL ASSERT(NT_RES.LE.NREM,'>>NREM too small',ERROR,*9999)
        CALL ASSERT(NTOPTI.LE.NOPM,'>>NOPM too small',ERROR,*9999)
        NCNLN=0

        CALL MINLSSQPOPTION('Default',OPTION)
        CALL MINLSSQPOPTION('Function Tolerance = 1.0D-12',OPTION)
        CALL MINLSSQPOPTIONR('Optimality Tolerance',OPTTOL,OPTION)
        CALL MINLSSQPOPTIONL('Verbose',(IPPLEV.GE.10),OPTION)
        CALL MINLSSQPOPTIONL('QP Verb',(MOD(IPPLEV,20).EQ.0),OPTION)
        CALL MINLSSQPOPTIONL('Line Verb',(MOD(IPPLEV,30).EQ.0),OPTION)
        CALL MINLSSQPOPTIONR('Linesearch Lamda',LNSTOL,OPTION)
        CALL MINLSSQPOPTIONR('Line Stepmax',STEPLIM,OPTION)
        CALL MINLSSQPOPTIONI('Derivative level',IPDLEV,OPTION)
        CALL MINLSSQPOPTIONI('Verify level',IPVLEV,OPTION)
        CALL MINLSSQPOPTIONI('Start Verify at Variable',ISROCV,OPTION)
        CALL MINLSSQPOPTIONI('Stop Verify at Variable',ISPOCV,OPTION)
        IF(KTYP27.EQ.5.OR.KTYP27.EQ.6) THEN
          CALL MINLSSQPOPTIONL('Central Differences',.FALSE.,OPTION)
        ELSE
          CALL MINLSSQPOPTIONL('Central Differences',.TRUE.,OPTION)
        ENDIF
        IF(DIFFER.GT.0.0D0) THEN !difference has been set by user
          CALL MINLSSQPOPTIONR('Difference DX',DIFFER,OPTION)
        ENDIF
        CALL MINLSSQPOPTIONL('Warm Start',WARMST,OPTION)

        DO noopti=1,NTOPTI
          XC(noopti)=PAOPTI(noopti)
          if(dop) then
            WRITE(OP_STRING,'('' Initial Value '',I6,'': '',D10.3)')
     '        noopti,XC(noopti)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          endif
        ENDDO

        CALL CPU_TIMER(CPU_USER,TIME_START)
        CALL MINLSSQP(OBJF,XC,RESID,RESJAC,LDFJ,NT_RES,NTOPTI,ITER,
     '    IEVAL,R,LDR,PMIN,PMAX,ISTATE,FUNCT2,IUSER,USER,OPTION,IFAIL)
        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)
        DO noopti=1,NTOPTI
          PAOPTI(noopti)=XC(noopti)
          WRITE(OP_STRING,'('' Parameter '',I2,'':'',D25.15)')
     '    noopti,XC(noopti)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDDO

        IF(OUTPTPRN) THEN
          WRITE(OP_STRING,'(/'' ****** Material Parameter '
     '      //'Optimisation ******'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Optimisation time:'',D11.4,'' s'')')
     '      ELAPSED_TIME
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          IF(IPPLEV.EQ.0) THEN
            IF(ifail.eq.0) THEN
              WRITE(OP_STRING,'('' Objective function = '',D12.4)') OBJF
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' Took '',I3,'' iterations.'')') iter
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF
C VYW/MPN 11/03/2011 Fixing long-time bug in covariance calculation
C ***     Calc Cholesky factor of Hessian matrix R
          CALL DPOTRF('U',NTOPTI,R,NOPM,ERR)
C ***     Calc inverse of Hessian matrix where R is its Cholesky factor
C         CALL F07FJF('U',NTOPTI,R,NOPM,ERR)
          CALL DPOTRI('U',NTOPTI,R,NOPM,ERR)
          IF(ERR.EQ.0) THEN
C KFA 2003/01/17 removed upper limit of 20 NTOPTI, and added
C  paging to the output to print more columns
            IF(NTOPTI.LT.NT_RES.AND.NTOPTI.LT.NOPM) THEN
C ***         Calculate and output parameter covariance matrix
              WRITE(OP_STRING,'(/'' Parameter Covariances:'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              DO i=1,NTOPTI
                DO j=i,NTOPTI
                  COVAR(i,j)=2.0D0*OBJF*R(i,j)/DBLE(NT_RES-NTOPTI)
                ENDDO
C                WRITE(CHAR8,FMT='(I3)') (i-1)*12+1
C                WRITE(OP_STRING,'('//CHAR8//'X,10D12.4)')
C     '            (COVAR(i,j),j=i,NTOPTI)
C                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDDO
              DO k=1,NTOPTI,5  !This is the number of pages
                endcol = k+4
                IF(endcol.GT.NTOPTI) endcol=NTOPTI
                WRITE(OP_STRING,'(''noopti='',I3,'','',I3)')
     '            k,endcol
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                DO i=1,k  !This is the number of complete rows
                  WRITE(OP_STRING,'(5D12.4)')
     '              (COVAR(i,j),j=k,endcol)
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDDO
                DO i=k+1,endcol !This is the incomplete rows
                  WRITE(CHAR8,FMT='(I3)') (i-k)*12
                  WRITE(OP_STRING,'('//CHAR8//'X,5D12.4)')
     '              (COVAR(i,j),j=i,endcol)
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDDO
              ENDDO
              CALCCORR=.TRUE.
              DO i=1,NTOPTI
                IF(COVAR(i,i).LT.1.0D-10) CALCCORR=.FALSE.
              ENDDO
              IF(CALCCORR) THEN
C ***           Calculate and output parameter correlation matrix
                WRITE(OP_STRING,'('' Parameter Correlations:'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                DO i=1,NTOPTI
                  DO j=i,NTOPTI
                    CORREL(i,j)=COVAR(i,j)/DSQRT(COVAR(i,i)*COVAR(j,j))
                  ENDDO
C                  WRITE(CHAR8,FMT='(I3)') (i-1)*12+1
C                  WRITE(OP_STRING,'('//CHAR8//'X,10D12.4)')
C     '              (CORREL(i,j),j=i,NTOPTI)
C                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDDO
                DO k=1,NTOPTI,5  !This is the number of pages
                  endcol = k+4
                  IF(endcol.GT.NTOPTI) endcol=NTOPTI
                  WRITE(OP_STRING,'(''noopti='',I3,'','',I3)')
     '              k,endcol
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  DO i=1,k  !This is the number of complete rows
                    WRITE(OP_STRING,'(5D12.4)')
     '                (CORREL(i,j),j=k,endcol)
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ENDDO
                  DO i=k+1,endcol !This is the incomplete rows
                    WRITE(CHAR8,FMT='(I3)') (i-k)*12
                    WRITE(OP_STRING,'('//CHAR8//'X,5D12.4)')
     '                (CORREL(i,j),j=i,endcol)
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ENDDO
                ENDDO
              ELSE  !CALCCORR=false
                WRITE(OP_STRING,'('' >>Cannot calculate '
     '            //'correlations of params with zero covariance(s)'')')
                CALL WRITES(IOER,OP_STRING,ERROR,*9999)
              ENDIF
            ELSE IF(NTOPTI.LT.NT_RES) THEN  !NTOPTI > 20
              WRITE(OP_STRING,'('' >>Increase dim of CORREL & COVAR'')')
              CALL WRITES(IOER,OP_STRING,ERROR,*9999)
            ELSE  !NTOPTI > NTRES
              WRITE(OP_STRING,'('' >>Cannot calculate '
     '          //'covariances unless #opt pars < #resids '')')
              CALL WRITES(IOER,OP_STRING,ERROR,*9999)
            ENDIF
          ELSE  !ERR <> 0
            WRITE(OP_STRING,'('' Could not find the inverse of '
     '        //'the Hessian: ERR='',I5,'' in DPOTRI'')') ERR
            CALL WRITES(IOER,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF

      ELSE IF((KTYP26.EQ.2.AND.KTYP27.EQ.6).OR.   !fluid interface
     '        (KTYP26.EQ.2.AND.KTYP27.EQ.1)) THEN !min squared areas

        WRITE(OP_STRING,'('' NOT IN OPERATION'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      ELSE IF(KTYP26.EQ.2.AND.KTYP27.EQ.2) THEN !Min curvatures

      ELSE IF(KTYP26.EQ.2.AND.KTYP27.EQ.3) THEN !Zero flux differences

        WRITE(OP_STRING,'('' NOT IN OPERATION'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      ELSE IF(KTYP26.EQ.2.AND.KTYP27.EQ.5) THEN !Data fitting
        CALL ASSERT(NREM.GE.NT_RES,'>>NREM too small',ERROR,*9999)
        CALL ASSERT(NOPM.GE.NTOPTI,'>>NOPM too small',ERROR,*9999)
        CALL ASSERT(NOPM.GE.NTCNTR,'>>NCOM too small',ERROR,*9999)
        DO noopti=1,NTOPTI
          XC(noopti)=PAOPTI(noopti)
        ENDDO
        NCNLN=NTCNTR
        CALL ASSERT(NCNLN.EQ.0,'>>MinSQP doesn''t have non-bound '//
     '    'constraints',ERROR,*9999)

C        IF(KTYP12.EQ.1) THEN !Sobolev smoothing
CC CPB 13/4/94 Calculate the scaling factor for the Sobolev smoothing
CC as the ratio of the Sobolev value for the mesh to the sum of squares
CC the data errors for the element
C         CO(1)='FEM'
C          CO(2)='UPDATE'
C          CO(3)='SOBOLEV'
C          NTCO=3
C          noco=1
C          CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,
C    '      ERROR,*9999)
C        ENDIF

        CALL MINLSSQPOPTION('Default',OPTION)
        CALL MINLSSQPOPTION('Function Tolerance = 1.0D-12',OPTION)
        CALL MINLSSQPOPTIONR('Optimality Tolerance',OPTTOL,OPTION)
        CALL MINLSSQPOPTIONL('Verbose',(IPPLEV.GE.10),OPTION)
        CALL MINLSSQPOPTIONL('QP Verb',(MOD(IPPLEV,20).EQ.0),OPTION)
        CALL MINLSSQPOPTIONL('Line Verb',(MOD(IPPLEV,30).EQ.0),OPTION)
        CALL MINLSSQPOPTIONR('Linesearch Lamda',LNSTOL,OPTION)
        CALL MINLSSQPOPTIONR('Line Stepmax',STEPLIM,OPTION)
        CALL MINLSSQPOPTIONI('Derivative level',IPDLEV,OPTION)
        CALL MINLSSQPOPTIONI('Verify level',IPVLEV,OPTION)
        CALL MINLSSQPOPTIONI('Start Verify at Variable',ISROCV,OPTION)
        CALL MINLSSQPOPTIONI('Stop Verify at Variable',ISPOCV,OPTION)
        CALL MINLSSQPOPTIONL('Central Differences',.TRUE.,OPTION)
        IF(DIFFER.GT.0.0D0) THEN !difference has been set by user
          CALL MINLSSQPOPTIONR('Difference DX',DIFFER,OPTION)
        ENDIF
        CALL MINLSSQPOPTIONL('Warm Start',WARMST,OPTION)
C       CALL MINLSSQPOPTIONR('Nonlinear Feasibility Tolerance ',NLFTOL,OPTION)
C       Hessian is constant so there should be no need to reset it.
        CALL MINLSSQPOPTIONI('Reset Frequency',0,OPTION)
        IF(KTYP29B.EQ.2) CALL MINLSSQPOPTION('Unit Initial Hessian',
     '    OPTION)

        CALL CPU_TIMER(CPU_USER,TIME_START)
        CALL MINLSSQP(OBJF,XC,RESID,RESJAC,LDFJ,NT_RES,NTOPTI,ITER,
     '    IEVAL,R,LDR,PMIN,PMAX,ISTATE,FUNCT2,IUSER,USER,OPTION,IFAIL)
        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)

        DO noopti=1,NTOPTI
          PAOPTI(noopti)=XC(noopti)
        ENDDO

        IF(ifail.eq.0) THEN
          WRITE(OP_STRING,'(/'' Objective function = '',D12.4)') OBJF
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Took '',I4,'' iterations.'')') iter
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C XSL 23Aug2010 IWRIT5 is defined in iwrit00.cmn as IWRIT5(nr,nx)
C But for some cases nr, nx_opt were not initialised in IPOPTI
C Hence rename IWRIT5 to IWRIT5_OPTI
          IF(IWRIT5_OPTI.GT.1) THEN
            WRITE(OP_STRING,'('' Solution time:'',D11.4,'' seconds'')')
     '        ELAPSED_TIME
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF !IWRIT5
        ENDIF

      ELSE IF(KTYP26.EQ.2.AND.KTYP27.EQ.7) THEN !aerofoil wake & stress
        CALL ASSERT(NREM.GE.NT_RES,'>>NREM too small',ERROR,*9999)
        CALL ASSERT(NOPM.GE.NTOPTI,'>>NOPM too small',ERROR,*9999)
        CALL ASSERT(NOPM.GE.NTCNTR,'>>NCOM too small',ERROR,*9999)

        DO noopti=1,NTOPTI
          XC(noopti)=PAOPTI(noopti)
        ENDDO
        NCNLN=NTCNTR
        CALL ASSERT(NCNLN.EQ.0,'>>MinSQP doesn''t have non-bound '//
     '    'constraints',ERROR,*9999)

        CALL MINLSSQPOPTION('Default',OPTION)
        CALL MINLSSQPOPTION('Function Tolerance = 1.0D-12',OPTION)
        CALL MINLSSQPOPTIONR('Optimality Tolerance',OPTTOL,OPTION)
        CALL MINLSSQPOPTIONL('Verbose',(IPPLEV.GE.10),OPTION)
        CALL MINLSSQPOPTIONL('QP Verb',(MOD(IPPLEV,20).EQ.0),OPTION)
        CALL MINLSSQPOPTIONL('Line Verb',(MOD(IPPLEV,30).EQ.0),OPTION)
        CALL MINLSSQPOPTIONR('Linesearch Lamda',LNSTOL,OPTION)
        CALL MINLSSQPOPTIONR('Line Stepmax',STEPLIM,OPTION)
        CALL MINLSSQPOPTIONI('Derivative level',IPDLEV,OPTION)
        CALL MINLSSQPOPTIONI('Verify level',IPVLEV,OPTION)
        CALL MINLSSQPOPTIONI('Start Verify at Variable',ISROCV,OPTION)
        CALL MINLSSQPOPTIONI('Stop Verify at Variable',ISPOCV,OPTION)
        CALL MINLSSQPOPTIONL('Central Differences',.TRUE.,OPTION)
        IF(DIFFER.GT.0.0D0) THEN !difference has been set by user
          CALL MINLSSQPOPTIONR('Difference DX',DIFFER,OPTION)
        ENDIF
        CALL MINLSSQPOPTIONL('Warm Start',WARMST,OPTION)
C       CALL MINLSSQPOPTIONR('Nonlinear Feasibility Tolerance ',NLFTOL,OPTION)

        CALL CPU_TIMER(CPU_USER,TIME_START)
        CALL MINLSSQP(OBJF,XC,RESID,RESJAC,LDFJ,NT_RES,NTOPTI,ITER,
     '    IEVAL,R,LDR,PMIN,PMAX,ISTATE,FUNCT2,IUSER,USER,OPTION,IFAIL)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)

        DO noopti=1,NTOPTI
          PAOPTI(noopti)=XC(noopti)
        ENDDO

        IF(ifail.eq.0) THEN
          WRITE(OP_STRING,'('' Objective function = '',D12.4)') OBJF
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Took '',I3,'' iterations.'')') iter
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

      ELSE IF(KTYP26.EQ.2.AND.KTYP27.EQ.8) THEN !aerofoil lift & wake
        CALL ASSERT(NREM.GE.NT_RES,'>>NREM too small',ERROR,*9999)
        CALL ASSERT(NOPM.GE.NTOPTI,'>>NOPM too small',ERROR,*9999)
        CALL ASSERT(NOPM.GE.NTCNTR,'>>NCOM too small',ERROR,*9999)
        DO noopti=1,NTOPTI
          XC(noopti)=PAOPTI(noopti)
        ENDDO
        NCNLN=NTCNTR
        CALL ASSERT(NCNLN.EQ.0,'>>MinSQP doesn''t have non-bound '//
     '    'constraints',ERROR,*9999)

        CALL MINLSSQPOPTION('Default',OPTION)
        CALL MINLSSQPOPTION('Function Tolerance = 1.0D-12',OPTION)
        CALL MINLSSQPOPTIONR('Optimality Tolerance',OPTTOL,OPTION)
        CALL MINLSSQPOPTIONL('Verbose',(IPPLEV.GE.10),OPTION)
        CALL MINLSSQPOPTIONL('QP Verb',(MOD(IPPLEV,20).EQ.0),OPTION)
        CALL MINLSSQPOPTIONL('Line Verb',(MOD(IPPLEV,30).EQ.0),OPTION)
        CALL MINLSSQPOPTIONR('Linesearch Lamda',LNSTOL,OPTION)
        CALL MINLSSQPOPTIONR('Line Stepmax',STEPLIM,OPTION)
        CALL MINLSSQPOPTIONI('Derivative level',IPDLEV,OPTION)
        CALL MINLSSQPOPTIONI('Verify level',IPVLEV,OPTION)
        CALL MINLSSQPOPTIONI('Start Verify at Variable',ISROCV,OPTION)
        CALL MINLSSQPOPTIONI('Stop Verify at Variable',ISPOCV,OPTION)
        CALL MINLSSQPOPTIONL('Central Differences',.TRUE.,OPTION)
        IF(DIFFER.GT.0.0D0) THEN !difference has been set by user
          CALL MINLSSQPOPTIONR('Difference DX',DIFFER,OPTION)
        ENDIF
        CALL MINLSSQPOPTIONL('Warm Start',WARMST,OPTION)
C       CALL MINLSSQPOPTIONR('Nonlinear Feasibility Tolerance ',NLFTOL,OPTION)

        CALL CPU_TIMER(CPU_USER,TIME_START)
        CALL MINLSSQP(OBJF,XC,RESID,RESJAC,LDFJ,NT_RES,NTOPTI,ITER,
     '    IEVAL,R,LDR,PMIN,PMAX,ISTATE,FUNCT2,IUSER,USER,OPTION,IFAIL)
        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)

        DO noopti=1,NTOPTI
          PAOPTI(noopti)=XC(noopti)
        ENDDO

        IF(ifail.eq.0) THEN
          WRITE(OP_STRING,'(/'' Objective function = '',D12.4)') OBJF
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Took '',I3,'' iterations.'')') iter
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          IF(IWRIT5_OPTI.GT.1) THEN
            WRITE(OP_STRING,'('' Solution time:'',D11.4,'' seconds'')')
     '        ELAPSED_TIME
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF !IWRIT5
        ENDIF
      ELSE IF(KTYP26.EQ.2.AND.KTYP27.EQ.10) THEN !torso mesh
        CALL ASSERT(NREM.GE.NT_RES,'>>NREM too small',ERROR,*9999)
        CALL ASSERT(NOPM.GE.NTOPTI,'>>NOPM too small',ERROR,*9999)
        DO noopti=1,NTOPTI
          XC(noopti)=PAOPTI(noopti)
        ENDDO
        NCNLN=0

        CALL MINLSSQPOPTION('Default',OPTION)
        CALL MINLSSQPOPTION('Function Tolerance = 1.0D-12',OPTION)
        CALL MINLSSQPOPTIONR('Optimality Tolerance',OPTTOL,OPTION)
        CALL MINLSSQPOPTIONL('Verbose',(IPPLEV.GE.10),OPTION)
        CALL MINLSSQPOPTIONL('QP Verb',(MOD(IPPLEV,20).EQ.0),OPTION)
        CALL MINLSSQPOPTIONL('Line Verb',(MOD(IPPLEV,30).EQ.0),OPTION)
        CALL MINLSSQPOPTIONR('Linesearch Lamda',LNSTOL,OPTION)
        CALL MINLSSQPOPTIONR('Line Stepmax',STEPLIM,OPTION)
        CALL MINLSSQPOPTIONI('Derivative level',IPDLEV,OPTION)
        CALL MINLSSQPOPTIONI('Verify level',IPVLEV,OPTION)
        CALL MINLSSQPOPTIONI('Start Verify at Variable',ISROCV,OPTION)
        CALL MINLSSQPOPTIONI('Stop Verify at Variable',ISPOCV,OPTION)
        CALL MINLSSQPOPTIONL('Central Differences',.TRUE.,OPTION)
        IF(DIFFER.GT.0.0D0) THEN !difference has been set by user
          CALL MINLSSQPOPTIONR('Difference DX',DIFFER,OPTION)
        ENDIF
        CALL MINLSSQPOPTIONL('Warm Start',WARMST,OPTION)

        CALL CPU_TIMER(CPU_USER,TIME_START)
        CALL MINLSSQP(OBJF,XC,RESID,RESJAC,LDFJ,NT_RES,NTOPTI,ITER,
     '    IEVAL,R,LDR,PMIN,PMAX,ISTATE,FUNCT2,IUSER,USER,OPTION,IFAIL)
        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)

        DO noopti=1,NTOPTI
          PAOPTI(noopti)=XC(noopti)
        ENDDO

        IF(IFAIL.EQ.0) THEN
          WRITE(OP_STRING,'(/'' Objective function = '',D12.4)') OBJF
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Took '',I3,'' iterations.'')') iter
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          IF(IWRIT5_OPTI.GT.1) THEN
            WRITE(OP_STRING,'('' Solution time:'',D11.4,'' seconds'')')
     '        ELAPSED_TIME
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF !IWRIT5
        ENDIF
      ELSE IF(KTYP26.EQ.2.AND.KTYP27.EQ.11) THEN !2nd moments
        CALL ASSERT(NREM.GE.NT_RES,'>>NREM too small',ERROR,*9999)
        CALL ASSERT(NOPM.GE.NTOPTI,'>>NOPM too small',ERROR,*9999)
        DO noopti=1,NTOPTI
          XC(noopti)=PAOPTI(noopti)
        ENDDO
        NCNLN=0

        CALL MINLSSQPOPTION('Default',OPTION)
        CALL MINLSSQPOPTION('Function Tolerance = 1.0D-12',OPTION)
        CALL MINLSSQPOPTIONR('Optimality Tolerance',OPTTOL,OPTION)
        CALL MINLSSQPOPTIONL('Verbose',(IPPLEV.GE.10),OPTION)
        CALL MINLSSQPOPTIONL('QP Verb',(MOD(IPPLEV,20).EQ.0),OPTION)
        CALL MINLSSQPOPTIONL('Line Verb',(MOD(IPPLEV,30).EQ.0),OPTION)
        CALL MINLSSQPOPTIONR('Linesearch Lamda',LNSTOL,OPTION)
        CALL MINLSSQPOPTIONR('Line Stepmax',STEPLIM,OPTION)
        CALL MINLSSQPOPTIONI('Derivative level',IPDLEV,OPTION)
        CALL MINLSSQPOPTIONI('Verify level',IPVLEV,OPTION)
        CALL MINLSSQPOPTIONI('Start Verify at Variable',ISROCV,OPTION)
        CALL MINLSSQPOPTIONI('Stop Verify at Variable',ISPOCV,OPTION)
        CALL MINLSSQPOPTIONL('Central Differences',.TRUE.,OPTION)
        IF(DIFFER.GT.0.0D0) THEN !difference has been set by user
          CALL MINLSSQPOPTIONR('Difference DX',DIFFER,OPTION)
        ENDIF
        CALL MINLSSQPOPTIONL('Warm Start',WARMST,OPTION)

        CALL CPU_TIMER(CPU_USER,TIME_START)
        CALL MINLSSQP(OBJF,XC,RESID,RESJAC,LDFJ,NT_RES,NTOPTI,ITER,
     '    IEVAL,R,LDR,PMIN,PMAX,ISTATE,FUNCT2,IUSER,USER,OPTION,IFAIL)
        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)

        DO noopti=1,NTOPTI
          PAOPTI(noopti)=XC(noopti)
        ENDDO

        IF(IFAIL.EQ.0) THEN
          WRITE(OP_STRING,'(/'' Objective function = '',D12.4)') OBJF
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Took '',I3,'' iterations.'')') iter
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          IF(IWRIT5_OPTI.GT.1) THEN
            WRITE(OP_STRING,'('' Solution time:'',D11.4,'' seconds'')')
     '        ELAPSED_TIME
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF !IWRIT5
        ENDIF
      ELSE IF(KTYP26.EQ.2.AND.KTYP27.EQ.12) THEN !Activation time
        CALL ASSERT(NREM.GE.NT_RES,'>>NREM too small',ERROR,*9999)
        CALL ASSERT(NOPM.GE.NTOPTI,'>>NOPM too small',ERROR,*9999)
        DO noopti=1,NTOPTI
          XC(noopti)=PAOPTI(noopti)
        ENDDO
        NCNLN=0

        CALL MINLSSQPOPTION('Default',OPTION)
        CALL MINLSSQPOPTIONR('Line Stepmax',STEPLIM,OPTION)
        CALL MINLSSQPOPTIONL('Verbose',(IPPLEV.GE.10),OPTION)
        CALL MINLSSQPOPTIONL('QP Verb',(IPPLEV.GE.20),OPTION)
        CALL MINLSSQPOPTIONL('Line Verb',(IPPLEV.GE.30),OPTION)
c       CALL MINLSSQPOPTIONI('QP Output File',6,OPTION)
        CALL MINLSSQPOPTIONR('Linesearch Lamda',LNSTOL,OPTION)
        CALL MINLSSQPOPTIONR('Optimality Tolerance',OPTTOL,OPTION)
        CALL MINLSSQPOPTIONR('Function Tolerance',FUNPREC,OPTION)
        CALL MINLSSQPOPTIONI('Maximum Iterations',MAX_MAJOR_ITER,OPTION)
        CALL MINLSSQPOPTIONI('QP Maximum Iterations',MAX_MINOR_ITER,
     '    OPTION)
        CALL MINLSSQPOPTIONI('Derivative level',IPDLEV,OPTION)
        CALL MINLSSQPOPTIONI('Verify level',IPVLEV,OPTION)
        CALL MINLSSQPOPTIONI('Start Verify at Variable',ISROCV,OPTION)
        CALL MINLSSQPOPTIONI('Stop Verify at Variable',ISPOCV,OPTION)
        CALL MINLSSQPOPTIONL('Central Differences',.FALSE.,OPTION)
        IF(DIFFER.GT.0.0D0) THEN !difference has been set by user
          CALL MINLSSQPOPTIONR('Difference DX',DIFFER,OPTION)
        ENDIF
        CALL MINLSSQPOPTIONL('Warm Start',WARMST,OPTION)
C KAT: For backward compatibility.  The Hessian is not J^t J but
C      sometimes setting the Hessian approximation to this improves
C      performance because the real Hessian is changing so fast.  A
C      better solution would be to reformulate the residuals so that
C      they are not squares.  The Hessian is then more constant and the
C      optimizer can better estimate the Hessian.
        CALL MINLSSQPOPTIONI('Reset Frequency',2,OPTION)
C SEN: Can make a big difference in some cases
C      Need to be made options.
C       CALL MINLSSQPOPTIONL('Unit Hessian',.TRUE.,OPTION)
C       CALL MINLSSQPOPTIONI('Reset Frequency',0,OPTION)

        CALL CPU_TIMER(CPU_USER,TIME_START)
        CALL MINLSSQP(OBJF,XC,RESID,RESJAC,LDFJ,NT_RES,NTOPTI,ITER,
     '    IEVAL,R,LDR,PMIN,PMAX,ISTATE,FUNCT2,IUSER,USER,OPTION,IFAIL)
        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)

        DO noopti=1,NTOPTI
          PAOPTI(noopti)=XC(noopti)
        ENDDO


        WRITE(OP_STRING,'(/'' Objective function = '',D12.4)') OBJF
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Took '',I4,'' iterations'')') iter
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' And  '',I4,'' evaluations.'')') ieval
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        IF(IWRIT5_OPTI.GT.1) THEN
          WRITE(OP_STRING,'('' Solution time:'',D11.4,'' seconds'')')
     '      ELAPSED_TIME
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF !IWRIT5


      ELSE IF(KTYP27.EQ.13) THEN !Dipole inverses

        CALL ASSERT(NREM.GE.NT_RES,'>>NREM too small',ERROR,*9999)
        CALL ASSERT(NOPM.GE.NTOPTI,'>>NOPM too small',ERROR,*9999)
        DO noopti=1,NTOPTI
          XC(noopti)=PAOPTI(noopti)
        ENDDO
        NCNLN=0

        CALL MINLSSQPOPTION('Default',OPTION)
        CALL MINLSSQPOPTIONR('Line Stepmax',STEPLIM,OPTION)
        CALL MINLSSQPOPTIONL('Verbose',(IPPLEV.GE.10),OPTION)
        CALL MINLSSQPOPTIONL('QP Verb',(IPPLEV.GE.20),OPTION)
        CALL MINLSSQPOPTIONL('Line Verb',(IPPLEV.GE.30),OPTION)
        CALL MINLSSQPOPTIONR('Linesearch Lamda',LNSTOL,OPTION)
        CALL MINLSSQPOPTIONR('Optimality Tolerance',OPTTOL,OPTION)
        CALL MINLSSQPOPTIONR('Function Tolerance',FUNPREC,OPTION)
        CALL MINLSSQPOPTIONI('Maximum Iterations',MAX_MAJOR_ITER,OPTION)
        CALL MINLSSQPOPTIONI('QP Maximum Iterations',MAX_MINOR_ITER,
     '    OPTION)
        CALL MINLSSQPOPTIONI('Derivative level',IPDLEV,OPTION)
        CALL MINLSSQPOPTIONI('Verify level',IPVLEV,OPTION)
        CALL MINLSSQPOPTIONI('Start Verify at Variable',ISROCV,OPTION)
        CALL MINLSSQPOPTIONI('Stop Verify at Variable',ISPOCV,OPTION)
        CALL MINLSSQPOPTIONL('Central Differences',.FALSE.,OPTION)
        IF(DIFFER.GT.0.0D0) THEN !difference has been set by user
          CALL MINLSSQPOPTIONR('Difference DX',DIFFER,OPTION)
        ENDIF
        CALL MINLSSQPOPTIONL('Warm Start',WARMST,OPTION)
C KAT: For backward compatibility with the old minlssqp default.  The
C       Hessian is not J^t*J but sometimes setting the Hessian approximation
C       to this improves performance because the real Hessian is changing so
C       fast.  A better solution would be to reformulate the residuals so that
C       they are not squares.  The Hessian is then more constant and the
C       optimizer can better estimate the Hessian.
        CALL MINLSSQPOPTIONI('Reset Frequency',2,OPTION)
C SEN: Can make a big difference in some cases
C      Need to be made options.
C       CALL MINLSSQPOPTIONL('Unit Hessian',.TRUE.,OPTION)
C       CALL MINLSSQPOPTIONI('Reset Frequency',0,OPTION)

        CALL CPU_TIMER(CPU_USER,TIME_START)
        CALL MINLSSQP(OBJF,XC,RESID,RESJAC,LDFJ,NT_RES,NTOPTI,ITER,
     '    IEVAL,R,LDR,PMIN,PMAX,ISTATE,FUNCT2,IUSER,USER,OPTION,IFAIL)
        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)

        DO noopti=1,NTOPTI
          PAOPTI(noopti)=XC(noopti)
        ENDDO

        WRITE(OP_STRING,'(/'' Objective function = '',D12.4)') OBJF
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Took '',I4,'' iterations'')') iter
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' And  '',I4,'' evaluations.'')') ieval
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        IF(IWRIT5_OPTI.GT.1) THEN
          WRITE(OP_STRING,'('' Solution time:'',D11.4,'' seconds'')')
     '      ELAPSED_TIME
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF !IWRIT5


      ELSE IF(KTYP26.EQ.3) THEN !Micro-structure
        nr=1
        WRITE(OP_STRING,
     '    '(/'' >>WARNING: Region number is fixed to 1'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        NCNLN=0
        CALL ASSERT(NT_RES.LE.NREM,'>>NREM too small',ERROR,*9999)
        CALL ASSERT(NTOPTI.LE.NOPM,'>>NOPM too small',ERROR,*9999)

        CALL MINLSSQPOPTION('Default',OPTION)
        CALL MINLSSQPOPTION('Function Tolerance = 1.0D-12',OPTION)
        CALL MINLSSQPOPTIONR('Optimality Tolerance',OPTTOL,OPTION)
        CALL MINLSSQPOPTIONL('Verbose',(IPPLEV.GE.10),OPTION)
        CALL MINLSSQPOPTIONL('QP Verb',(MOD(IPPLEV,20).EQ.0),OPTION)
        CALL MINLSSQPOPTIONL('Line Verb',(MOD(IPPLEV,30).EQ.0),OPTION)
        CALL MINLSSQPOPTIONR('Linesearch Lamda',LNSTOL,OPTION)
        CALL MINLSSQPOPTIONR('Line Stepmax',STEPLIM,OPTION)
        CALL MINLSSQPOPTIONI('Derivative level',IPDLEV,OPTION)
        CALL MINLSSQPOPTIONI('Verify level',IPVLEV,OPTION)
        CALL MINLSSQPOPTIONI('Start Verify at Variable',ISROCV,OPTION)
        CALL MINLSSQPOPTIONI('Stop Verify at Variable',ISPOCV,OPTION)
        CALL MINLSSQPOPTIONL('Central Differences',.TRUE.,OPTION)
        IF(DIFFER.GT.0.0D0) THEN !difference has been set by user
          CALL MINLSSQPOPTIONR('Difference DX',DIFFER,OPTION)
        ENDIF
        CALL MINLSSQPOPTIONL('Warm Start',WARMST,OPTION)

        DO noopti=1,NTOPTI
          XC(noopti)=PAOPTI(noopti)
        ENDDO

        CALL CPU_TIMER(CPU_USER,TIME_START)
        CALL MINLSSQP(OBJF,XC,RESID,RESJAC,LDFJ,NT_RES,NTOPTI,ITER,
     '    IEVAL,R,LDR,PMIN,PMAX,ISTATE,FUNCT2,IUSER,USER,OPTION,IFAIL)
        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)
        DO noopti=1,NTOPTI
          PAOPTI(noopti)=XC(noopti)
        ENDDO
        WRITE(OP_STRING,'('' Optimisation time:'',D11.4,'' s'')')
     '    ELAPSED_TIME

      ELSE IF(KTYP26.EQ.4) THEN !Holmes constants
        nr=1
        WRITE(OP_STRING,
     '    '(/'' >>WARNING: Region number is fixed to 1'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        NCNLN=0
        CALL ASSERT(NT_RES.LE.NREM,'>>NREM too small',ERROR,*9999)
        CALL ASSERT(NTOPTI.LE.NOPM,'>>NOPM too small',ERROR,*9999)

        CALL MINLSSQPOPTION('Default',OPTION)
        CALL MINLSSQPOPTION('Function Tolerance = 1.0D-12',OPTION)
        CALL MINLSSQPOPTIONR('Optimality Tolerance',OPTTOL,OPTION)
        CALL MINLSSQPOPTIONL('Verbose',(IPPLEV.GE.10),OPTION)
C       CALL MINLSSQPOPTIONI('Output file',11,OPTION)!So it isnt mixed with Careys crap
        CALL MINLSSQPOPTIONL('QP Verb',(MOD(IPPLEV,20).EQ.0),OPTION)
        CALL MINLSSQPOPTIONL('Line Verb',(MOD(IPPLEV,30).EQ.0),OPTION)
        CALL MINLSSQPOPTIONR('Linesearch Lamda',LNSTOL,OPTION)
        CALL MINLSSQPOPTIONR('Line Stepmax',STEPLIM,OPTION)
        CALL MINLSSQPOPTIONI('Derivative level',IPDLEV,OPTION)
        CALL MINLSSQPOPTIONI('Verify level',IPVLEV,OPTION)
        CALL MINLSSQPOPTIONI('Start Verify at Variable',ISROCV,OPTION)
        CALL MINLSSQPOPTIONI('Stop Verify at Variable',ISPOCV,OPTION)
        CALL MINLSSQPOPTIONL('Central Differences',.TRUE.,OPTION)
        IF(DIFFER.GT.0.0D0) THEN !difference has been set by user
        CALL MINLSSQPOPTIONR('Difference DX',DIFFER,OPTION)
        ENDIF
        CALL MINLSSQPOPTIONL('Warm Start',WARMST,OPTION)

        DO noopti=1,NTOPTI
          XC(noopti)=PAOPTI(noopti)
        ENDDO

        CALL CPU_TIMER(CPU_USER,TIME_START)
        CALL MINLSSQP(OBJF,XC,RESID,RESJAC,LDFJ,NT_RES,NTOPTI,ITER,
     '    IEVAL,R,LDR,PMIN,PMAX,ISTATE,FUNCT2,IUSER,USER,OPTION,IFAIL)
        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)
        DO noopti=1,NTOPTI
          PAOPTI(noopti)=XC(noopti)
        ENDDO
        WRITE(OP_STRING,'('' Optimisation time:'',D11.4,'' s'')')
     '    ELAPSED_TIME
      ELSE
        ERROR='Unknown optimisation problem'
        GOTO 9999
      ENDIF

      IF((KTYP27.LE.4.OR.KTYP27.EQ.7.OR.KTYP27.EQ.8).AND.OUTPTPRN) THEN
        WRITE(OP_STRING,'('' Parameters: '',/(1X,10D12.4))')
     '    (PAOPTI(noopti),noopti=1,NTOPTI)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      IF( IFAIL.NE.0 ) THEN
        CALL FLAG_ERROR(0,' ')
        CALL WRITE_CHAR(IOER,'MINLSSQP failed: IFAIL=',ERR)
        CALL WRITE_INT(IOER,IFAIL,ERR)
        CALL WRITE_CHAR(IOER,NEWLINE,ERR)
        ERROR=' '
        GOTO 9999
      ENDIF
C!!! LKC removing nr=1,NRT loop
C      ENDDO

      CALL EXITS('NAGMINA')
      RETURN
 9999 CALL ERRORS('NAGMINA',ERROR)
      CALL EXITS('NAGMINA')
      RETURN 1
      END


