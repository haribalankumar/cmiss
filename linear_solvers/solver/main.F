
      PROGRAM TEST_SOLVERS

      IMPLICIT NONE
      INCLUDE 'mach00.inc'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'solver.inc'
!     Parameters
      INTEGER M,N,LDA,NZA          ! Dimension of the equation set
      INTEGER NMAX,MMAX            ! Dimension of auxillary arrays
      INTEGER E_LEN,F_LEN
!     System arrays
      INTEGER SPARSE_A
      POINTER ISR_PTR
      POINTER ISC_PTR
      POINTER A_PTR
      POINTER B_PTR
      POINTER X_PTR
      POINTER D_PTR
      POINTER E_PTR
      POINTER F_PTR
      POINTER I_PTR
      POINTER J_PTR

!     Harwell stuff
C     INTEGER LICN,LIRN
      INTEGER NICN,NIRN
C     REAL*8 USET

!     Umfpack stuff
C     INTEGER KEEP(20),ICNTL(20),INFO(40),JOB,NINDX,NDARR
C     REAL*8 CNTL(10),RINFO(20)
      INTEGER NUMF
      LOGICAL UMF_ITER
      REAL*8 UMFDEF(40),UMF_CONTROL(20),UMF_INFO(90)

!     SuperLU stuff
      INTEGER COL_ORDER,SLU_FILL,PARAM(10),SLU_RELAX

!     Factorised crap
      INTEGER ITER,NRES,NPRECON,ITER_TMP
      INTEGER SOLVER,PRECON
      INTEGER OUTPUTCODE
      REAL*8 ANORM,RESID,OMEGA,RESID_TMP

!     Local Variables
      INTEGER I,N_FACT,N_RFAC,N_SOLV,NX,NY,NZ,IERR,ITMP,NCPU
      INTEGER PTYPE,MTYPE,KLUDGE,NNX,NNY,NNZ,MAXTIME,SETZERO
      LOGICAL DENSE,DUMP_FILE,WRITE_FILE,DUMP_FACT,LCHECK_SYMM_STRUCT
      LOGICAL LCHECK_IS_SYMM,LCHECK_FILL,LCHECK_DIAGDOM,ZERO_B,SOLN2D
      LOGICAL LMIX_ARRAYS,LRAND_WEIGHTS,LEXAMINE,LTEST,TRANSIENT
      LOGICAL LCHECK_MATRIX,WRITE_BOEING_FILE
      REAL*8 TIME(3,2,2),RDUM
      CHARACTER*128 BUFF,DUFF,FILENAME,ERROR,DUMPNAME,FACTNAME,WRITENAME
      CHARACTER*128 SOLNNAME,BOEINGNAME

!     Functions
      INTEGER IARGC,LEN_TRIM,ILAENV,ALARM
      REAL*8 TIMER_READ_WALL,TIMER_READ_BUSY,DDOT
      EXTERNAL TIMER_READ_WALL,TIMER_READ_BUSY
      EXTERNAL IARGC,LEN_TRIM,ILAENV,GETARG,DDOT,ALARM_WAKE,ALARM

      COMMON / ALARM_COMMON / MAXTIME

      DATA I_PTR,J_PTR,ISR_PTR,ISC_PTR         / 4*0 /
      DATA A_PTR,B_PTR,X_PTR                   / 3*0 /
      DATA D_PTR,E_PTR,F_PTR                   / 3*0 /
      DATA NICN,NIRN                           / 4,1 /
C     DATA NUMF,JOB                            / 5,0 /
      DATA COL_ORDER,KLUDGE,SLU_FILL,SLU_RELAX / 1,0,-20,0 /
      DATA LCHECK_IS_SYMM,LCHECK_FILL          / 2*.FALSE. /
      DATA LMIX_ARRAYS,LCHECK_DIAGDOM          / 2*.FALSE. /
      DATA LRAND_WEIGHTS,LEXAMINE              / 2*.FALSE. /
      DATA LCHECK_SYMM_STRUCT,LTEST            / 2*.FALSE. /
      DATA WRITE_BOEING_FILE,LCHECK_MATRIX     / 2*.FALSE. /
      DATA UMFDEF                              / 40*-1.0D0 /


      MTYPE=1
      PTYPE=0
      DUMP_FILE=.FALSE.
      DUMP_FACT=.FALSE.
      WRITE_FILE=.FALSE.
      DENSE=.FALSE.
      UMF_ITER=.FALSE.
      TRANSIENT=.TRUE.
      ZERO_B=.FALSE.
      SOLN2D=.FALSE.
      MAXTIME=-1

      SOLVER=0
      PRECON=1
      SETZERO=0

      OUTPUTCODE=2
      RESID=1.0D-9
      OMEGA=0.667
      ITER=1000
      NRES=0
      NPRECON=2
      NCPU=1
      FILENAME='input.dat'
      DUMPNAME='output.dat'
      FACTNAME='factor.dat'
      SOLNNAME='soln2d.gnu'
      WRITENAME='output.tri'

      N_FACT=1
      N_RFAC=0
      N_SOLV=1

      I=1
      DO WHILE(I.LE.IARGC())
        CALL GETARG(I,BUFF)
        IF(BUFF.EQ.'-file') THEN
          I=I+1
          CALL GETARG(I,FILENAME)
          PTYPE=10
        ELSE IF(BUFF.EQ.'-triplet') THEN
          I=I+1
          CALL GETARG(I,FILENAME)
          PTYPE=11
        ELSE IF(BUFF.EQ.'-soln2d_file') THEN
          I=I+1
          CALL GETARG(I,SOLNNAME)
          SOLN2D=.TRUE.

        ELSE IF(BUFF(1:2).EQ.'-v') THEN
          OUTPUTCODE=OUTPUTCODE+1

        ELSE IF(BUFF.EQ.'-dense') THEN
          DENSE=.TRUE.

        ELSE IF(BUFF.EQ.'-test') THEN
          LTEST=.TRUE.
        ELSE IF(BUFF.EQ.'-check_symm_struct') THEN
          LCHECK_SYMM_STRUCT=.TRUE.
        ELSE IF(BUFF.EQ.'-check_symm') THEN
          LCHECK_IS_SYMM=.TRUE.
        ELSE IF(BUFF.EQ.'-check_fill') THEN
          LCHECK_FILL=.TRUE.
        ELSE IF(BUFF.EQ.'-check_diagdom') THEN
          LCHECK_DIAGDOM=.TRUE.
        ELSE IF(BUFF.EQ.'-check_matrix') THEN
          LCHECK_MATRIX=.TRUE.
        ELSE IF(BUFF.EQ.'-mix_arrays') THEN
          LMIX_ARRAYS=.TRUE.
        ELSE IF(BUFF.EQ.'-rand_weights') THEN
          LRAND_WEIGHTS=.TRUE.
        ELSE IF(BUFF.EQ.'-examine') THEN
          LEXAMINE=.TRUE.
        ELSE IF(BUFF.EQ.'-kludge') THEN
          KLUDGE=KLUDGE+1
        ELSE IF(BUFF.EQ.'-nosolve') THEN
          SOLVER=0

        ELSE IF(BUFF.EQ.'-transient' .OR. BUFF.EQ.'-trans') THEN
          TRANSIENT=.TRUE.
        ELSE IF(BUFF.EQ.'-steady' .OR. BUFF.EQ.'-notrans') THEN
          TRANSIENT=.FALSE.

        ELSE IF(BUFF.EQ.'-dump_file') THEN
          DUMP_FILE=.TRUE.
          I=I+1
          CALL GETARG(I,DUMPNAME)

        ELSE IF(BUFF.EQ.'-write_file') THEN
          WRITE_FILE=.TRUE.
          I=I+1
          CALL GETARG(I,WRITENAME)

        ELSE IF(BUFF.EQ.'-write_file') THEN
          WRITE_FILE=.TRUE.
          I=I+1
          CALL GETARG(I,WRITENAME)

        ELSE IF(BUFF.EQ.'-write_boeing') THEN
          WRITE_BOEING_FILE=.TRUE.
          I=I+1
          CALL GETARG(I,BOEINGNAME)

        ELSE IF(BUFF(1:2).EQ.'-q') THEN
          OUTPUTCODE=0

        ELSE IF(BUFF(1:8).EQ.'-zerob') THEN
          ZERO_B=.TRUE.
        ELSE IF(BUFF(1:8).EQ.'-setzero') THEN
          SETZERO=1
          IF(BUFF(9:9).EQ.':') THEN
            READ(BUFF(10:),*,IOSTAT=IERR) SETZERO
            IF(IERR.NE.0 .OR. SETZERO.LT.0) THEN
              ERROR='>>Error reading setzero row number'
              GOTO 9999
            ENDIF
          ENDIF

        ELSE IF(BUFF.EQ.'-har_narr') THEN
          I=I+1
          CALL GETARG(I,BUFF)
          READ(BUFF,*,IOSTAT=IERR) NICN,NIRN
          IF(IERR.NE.0 .OR. NICN.LT.1 .OR. NIRN.LT.1) THEN
            ERROR='>>Error reading harwell array scale factors'
            GOTO 9999
          ENDIF

        ELSE IF(BUFF.EQ.'-maxtime' .OR. BUFF.EQ.'-alarm') THEN
          I=I+1
          CALL GETARG(I,BUFF)
          READ(BUFF,*,IOSTAT=IERR) MAXTIME
          IF(IERR.NE.0) THEN
            ERROR='>>Error reading maximum time'
            GOTO 9999
          ENDIF

C       UMF's many and varied parameters
        ELSE IF(BUFF.EQ.'-umf_narr') THEN
          I=I+1
          CALL GETARG(I,BUFF)
          READ(BUFF,*,IOSTAT=IERR) NUMF
          IF(IERR.NE.0 .OR. NUMF.LT.1) THEN
            ERROR='>>Error reading umfpack array scale factor'
            GOTO 9999
          ENDIF

        ELSE IF(BUFF.EQ.'-umf_iter') THEN
          UMFDEF(20)=1.0D0
        ELSE IF(BUFF.EQ.'-umf_noiter') THEN
          UMFDEF(20)=0.0D0

        ELSE IF(BUFF.EQ.'-umf_btri') THEN
          UMFDEF(3)=1.0D0
        ELSE IF(BUFF.EQ.'-umf_nobtri') THEN
          UMFDEF(3)=0.0D0

        ELSE IF(BUFF.EQ.'-umf_symm') THEN
          UMFDEF(5)=0.0D0
        ELSE IF(BUFF.EQ.'-umf_nosymm') THEN
          UMFDEF(5)=1.0D0

        ELSE IF(BUFF(1:5).EQ.'-umf_') THEN
          I=I+1
          CALL GETARG(I,DUFF)
          READ(DUFF,*,IOSTAT=IERR) RDUM
          IF(IERR.NE.0) THEN
            ERROR='>>Error reading '//BUFF
            GOTO 9999
          ENDIF

          IF(BUFF.EQ.'-umf_ptol') THEN
            UMFDEF(1)=RDUM
          ELSE IF(BUFF.EQ.'-umf_block') THEN
            UMFDEF(2)=RDUM

          ELSE IF(BUFF.EQ.'-umf_ncol') THEN
            UMFDEF(4)=RDUM
          ELSE IF(BUFF.EQ.'-umf_growth') THEN
            UMFDEF(6)=RDUM
          ELSE IF(BUFF.EQ.'-umf_d1') THEN
            UMFDEF(7)=RDUM
          ELSE IF(BUFF.EQ.'-umf_d2') THEN
            UMFDEF(8)=RDUM

          ELSE IF(BUFF.EQ.'-umf_drow') THEN
            UMFDEF(11)=RDUM
          ELSE IF(BUFF.EQ.'-umf_dcol') THEN
            UMFDEF(12)=RDUM
          ELSE IF(BUFF.EQ.'-umf_relax1') THEN
            UMFDEF(13)=RDUM
          ELSE IF(BUFF.EQ.'-umf_relax2') THEN
            UMFDEF(14)=RDUM
          ELSE IF(BUFF.EQ.'-umf_relax3') THEN
            UMFDEF(15)=RDUM

          ELSE
            ERROR='>>Unknown option '//BUFF
            GOTO 9999
          ENDIF

        ELSE IF(BUFF.EQ.'-slu_co') THEN
          I=I+1
          CALL GETARG(I,BUFF)
          READ(BUFF,*,IOSTAT=IERR) COL_ORDER
          IF(IERR.NE.0 .OR. COL_ORDER.LT.1) THEN
            ERROR='>>Error reading superLU column ordering'
            GOTO 9999
          ENDIF

        ELSE IF(BUFF.EQ.'-slu_fill') THEN
          I=I+1
          CALL GETARG(I,BUFF)
          READ(BUFF,*,IOSTAT=IERR) SLU_FILL
          IF(IERR.NE.0) THEN
            ERROR='>>Error reading superLU fill'
            GOTO 9999
          ENDIF

        ELSE IF(BUFF.EQ.'-slu_relax') THEN
          I=I+1
          CALL GETARG(I,BUFF)
          READ(BUFF,*,IOSTAT=IERR) SLU_RELAX
          IF(IERR.NE.0) THEN
            ERROR='>>Error reading superLU relax'
            GOTO 9999
          ENDIF

        ELSE IF(BUFF.EQ.'-finite_difference' .OR. BUFF.EQ.'-fd'
     '      .OR. BUFF.EQ.'-fd3D' .OR. BUFF.EQ.'-fd3d') THEN
          PTYPE=2
          I=I+1
          CALL GETARG(I,BUFF)
          READ(BUFF,*,IOSTAT=IERR) NX,NY,NZ
          IF(IERR.NE.0 .OR. NX.LT.3 .OR. NY.LT.3 .OR. NZ.LT.3) THEN
            ERROR='>>Error reading dimensions of finite difference mesh'
            GOTO 9999
          ENDIF

        ELSE IF(BUFF.EQ.'-fd2D' .OR. BUFF.EQ.'-fd2d') THEN
          PTYPE=3
          I=I+1
          CALL GETARG(I,BUFF)
          READ(BUFF,*,IOSTAT=IERR) NX,NY
          IF(IERR.NE.0 .OR. NX.LT.3 .OR. NY.LT.3) THEN
            ERROR='>>Error reading dimensions of finite difference mesh'
            GOTO 9999
          ENDIF

        ELSE IF(BUFF.EQ.'-finite_volume' .OR. BUFF.EQ.'-fv'
     '      .OR. BUFF.EQ.'-fv3D' .OR. BUFF.EQ.'-fv3d') THEN
          PTYPE=4
          I=I+1
          CALL GETARG(I,BUFF)
          READ(BUFF,*,IOSTAT=IERR) NX,NY,NZ
          IF(IERR.NE.0 .OR. NX.LT.3 .OR. NY.LT.3 .OR. NZ.LT.3) THEN
            ERROR='>>Error reading dimensions of finite volume mesh'
            GOTO 9999
          ENDIF

        ELSE IF(BUFF.EQ.'-fv2D' .OR. BUFF.EQ.'-fv2d') THEN
          PTYPE=5
          I=I+1
          CALL GETARG(I,BUFF)
          READ(BUFF,*,IOSTAT=IERR) NX,NY
          IF(IERR.NE.0 .OR. NX.LT.3 .OR. NY.LT.3) THEN
            ERROR='>>Error reading dimensions of finite volume mesh'
            GOTO 9999
          ENDIF

        ELSE IF(BUFF.EQ.'-fv2Db' .OR. BUFF.EQ.'-fv2db') THEN
          PTYPE=6
          I=I+1
          CALL GETARG(I,BUFF)
          READ(BUFF,*,IOSTAT=IERR) NX,NY
          IF(IERR.NE.0 .OR. NX.LT.3 .OR. NY.LT.3) THEN
            ERROR='>>Error reading dimensions of finite volume mesh'
            GOTO 9999
          ENDIF

        ELSE IF(BUFF.EQ.'-sinmesh') THEN
          MTYPE=2

        ELSE IF(BUFF(1:6).EQ.'-niter') THEN
          I=I+1
          CALL GETARG(I,BUFF)
          READ(BUFF,*,IOSTAT=IERR) ITER
          IF(IERR.NE.0 .OR. ITER.LT.0) THEN
            ERROR='>>Error reading ITER'
            GOTO 9999
          ENDIF

        ELSE IF(BUFF(1:6).EQ.'-nres') THEN
          I=I+1
          CALL GETARG(I,BUFF)
          READ(BUFF,*,IOSTAT=IERR) NRES
          IF(IERR.NE.0) THEN
            ERROR='>>Error reading NRES'
            GOTO 9999
          ENDIF

        ELSE IF(BUFF(1:8).EQ.'-nprecon') THEN
          I=I+1
          CALL GETARG(I,BUFF)
          READ(BUFF,*,IOSTAT=IERR) NPRECON
          IF(IERR.NE.0 .OR. NPRECON.LE.0) THEN
            ERROR='>>Error reading NPRECON'
            GOTO 9999
          ENDIF

        ELSE IF(BUFF(1:4).EQ.'-tol') THEN
          I=I+1
          CALL GETARG(I,BUFF)
          READ(BUFF,*,IOSTAT=IERR) RESID
          IF(IERR.NE.0 .OR. RESID.LE.0.0D0) THEN
            ERROR='>>Error reading RESID'
            GOTO 9999
          ENDIF

        ELSE IF(BUFF(1:6).EQ.'-omega') THEN
          I=I+1
          CALL GETARG(I,BUFF)
          READ(BUFF,*,IOSTAT=IERR) OMEGA
          IF(IERR.NE.0 .OR. OMEGA.LE.0.0D0) THEN
            ERROR='>>Error reading OMEGA'
            GOTO 9999
          ENDIF

        ELSE IF(BUFF(1:6).EQ.'-nfact') THEN
          I=I+1
          CALL GETARG(I,BUFF)
          READ(BUFF,*,IOSTAT=IERR) N_FACT
          IF(IERR.NE.0 .OR. N_FACT.LE.0) THEN
            ERROR='>>Error reading N_FACT'
            GOTO 9999
          ENDIF

        ELSE IF(BUFF(1:8).EQ.'-nrefact') THEN
          I=I+1
          CALL GETARG(I,BUFF)
          READ(BUFF,*,IOSTAT=IERR) N_RFAC
          IF(IERR.NE.0 .OR. N_RFAC.LT.0) THEN
            ERROR='>>Error reading N_RFACT'
            GOTO 9999
          ENDIF

        ELSE IF(BUFF(1:6).EQ.'-nsolv') THEN
          I=I+1
          CALL GETARG(I,BUFF)
          READ(BUFF,*,IOSTAT=IERR) N_SOLV
          IF(IERR.NE.0 .OR. N_SOLV.LT.0) THEN
            ERROR='>>Error reading N_SOLV'
            GOTO 9999
          ENDIF

        ELSE IF(BUFF(1:6).EQ.'-ncpu') THEN
          I=I+1
          CALL GETARG(I,BUFF)
          READ(BUFF,*,IOSTAT=IERR) NCPU
          IF(IERR.NE.0 .OR. NCPU.LE.0) THEN
            ERROR='>>Error reading NCPU'
            GOTO 9999
          ENDIF

        ELSE IF(BUFF.EQ.'-solver') THEN
          I=I+1
          CALL GETARG(I,BUFF)
          IF(BUFF.EQ.'none') THEN
            SOLVER=0
          ELSE IF(BUFF.EQ.'matmul') THEN
            SOLVER=30

          ELSE IF(BUFF.EQ.'jacobi') THEN
            SOLVER=SOLV_JACOBI
          ELSE IF(BUFF.EQ.'sor') THEN
            SOLVER=SOLV_SOR

          ELSE IF(BUFF.EQ.'ic0') THEN
            SOLVER=SOLV_IC0
          ELSE IF(BUFF.EQ.'ilu0') THEN
            SOLVER=SOLV_ILU0
          ELSE IF(BUFF.EQ.'ildl0') THEN
            SOLVER=SOLV_ILDL0

          ELSE IF(BUFF.EQ.'ic1') THEN
            SOLVER=SOLV_IC1
          ELSE IF(BUFF.EQ.'ilu1') THEN
            SOLVER=SOLV_ILU1
          ELSE IF(BUFF.EQ.'ildl1') THEN
            SOLVER=SOLV_ILDL1

          ELSE IF(BUFF.EQ.'cg') THEN
            SOLVER=SOLV_CG
          ELSE IF(BUFF.EQ.'bicgstab') THEN
            SOLVER=SOLV_BICGSTAB
          ELSE IF(BUFF.EQ.'gmres') THEN
            SOLVER=SOLV_GMRES

          ELSE IF(BUFF.EQ.'lu') THEN
            SOLVER=SOLV_LU
            DENSE=.TRUE.
          ELSE IF(BUFF.EQ.'cholesky') THEN
            SOLVER=SOLV_CHOLESKY
            DENSE=.TRUE.
          ELSE IF(BUFF.EQ.'ldl') THEN
            SOLVER=SOLV_LDL
            DENSE=.TRUE.
          ELSE IF(BUFF.EQ.'svd') THEN
            SOLVER=SOLV_SVD
            DENSE=.TRUE.
          ELSE IF(BUFF.EQ.'lsq') THEN
            SOLVER=SOLV_LSQ
            DENSE=.TRUE.
 
          ELSE IF(BUFF.EQ.'harwell') THEN
            SOLVER=SOLV_HARWELL
            DENSE=.FALSE.
          ELSE IF(BUFF.EQ.'superlu') THEN
            SOLVER=SOLV_SUPERLU
            DENSE=.FALSE.
          ELSE IF(BUFF.EQ.'umfpack' .OR. BUFF.EQ.'umfpack4') THEN
            SOLVER=SOLV_UMFPACK4
            DENSE=.FALSE.
          ELSE IF(BUFF.EQ.'umfpack2') THEN
            SOLVER=SOLV_UMFPACK2
            DENSE=.FALSE.

          ELSE
            ERROR='>>Unknown solver type '''
     '        //BUFF(:LEN_TRIM(BUFF))//''''
            GOTO 9999
          ENDIF

        ELSE IF(BUFF.EQ.'-precon') THEN
          I=I+1
          CALL GETARG(I,BUFF)
          IF(BUFF.EQ.'none') THEN
            PRECON=PRE_NONE
          ELSE IF(BUFF.EQ.'pj' .OR. BUFF.EQ.'point-jacobi') THEN
            PRECON=PRE_POINT
          ELSE IF(BUFF.EQ.'rowscale') THEN
            PRECON=PRE_SCALE

          ELSE IF(BUFF.EQ.'jacobi') THEN
            PRECON=PRE_JACOBI
          ELSE IF(BUFF.EQ.'ssor' .OR. BUFF.EQ.'sor') THEN
            PRECON=PRE_SOR

          ELSE IF(BUFF.EQ.'ic0') THEN
            PRECON=PRE_IC0
          ELSE IF(BUFF.EQ.'ilu0') THEN
            PRECON=PRE_ILU0
          ELSE IF(BUFF.EQ.'ildl0') THEN
            PRECON=PRE_ILDL0

          ELSE IF(BUFF.EQ.'ic1') THEN
            PRECON=PRE_IC1
          ELSE IF(BUFF.EQ.'ilu1') THEN
            PRECON=PRE_ILU1
          ELSE IF(BUFF.EQ.'ildl1') THEN
            PRECON=PRE_ILDL1

          ELSE
            ERROR='>>Unknown preconditioner type '''
     '        //BUFF(:LEN_TRIM(BUFF))//''''
            GOTO 9999
          ENDIF

        ELSE IF(BUFF(1:2).EQ.'-h' .OR. BUFF(1:2).EQ.'-?') THEN
          CALL USAGE()
          STOP

        ELSE
          ERROR='>>Unknown option '//BUFF
          GOTO 9999
        ENDIF

        I=I+1
      ENDDO

C     Check we have an array to work on

      IF(PTYPE.EQ.0) THEN
        ERROR='>>No problem defined'
        GOTO 9999
      ENDIF

C     Set the number of threads
C$    CALL OMP_SET_NUM_THREADS(NCPU)

C     Set up the problem
C     Read equations from vsm file
      IF(PTYPE.EQ.10) THEN
        WRITE(*,7400) 'Reading equations from ''',
     '    FILENAME(:LEN_TRIM(FILENAME)),''''
        CALL READ_SPARSE_DIM(N,NZA,FILENAME,ERROR,*9999)
        CALL ALLOCATE_MEMORY(NZA,1,INTTYPE,ISC_PTR,NOINIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(N+1,1,INTTYPE,ISR_PTR,NOINIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(NZA,1, DPTYPE,  A_PTR,NOINIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(  N,1, DPTYPE,  B_PTR,NOINIT,ERROR,*9999)
        NMAX=N
        MMAX=NZA
        CALL READ_SPARSE(%VAL(A_PTR),N,NMAX,MMAX,%VAL(B_PTR),
     '    %VAL(ISC_PTR),%VAL(ISR_PTR),FILENAME,ERROR,*9999)

C     Read equations from triplet file
      ELSE IF(PTYPE.EQ.11) THEN
        WRITE(*,7400) 'Reading equations from ''',
     '    FILENAME(:LEN_TRIM(FILENAME)),''''
        CALL READ_TRIPLET_DIM(N,NZA,FILENAME,ERROR,*9999)
        CALL ALLOCATE_MEMORY(NZA,1,INTTYPE,ISC_PTR,NOINIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(N+1,1,INTTYPE,ISR_PTR,NOINIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(NZA,1, DPTYPE,  A_PTR,NOINIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(  N,1, DPTYPE,  B_PTR,NOINIT,ERROR,*9999)
        NMAX=N
        MMAX=NZA
        CALL READ_TRIPLET(%VAL(A_PTR),N,NMAX,MMAX,%VAL(ISC_PTR),
     '    %VAL(ISR_PTR),FILENAME,ERROR,*9999)
        CALL PDLOAD(N,1.0D0,%VAL(B_PTR),1)

C     3D Finite Differences
      ELSE IF(PTYPE.EQ.2) THEN
        WRITE(*,7300) 'Equations generated on a ',NX,' x ',NY,
     '    ' x ',NZ,' FDiff mesh'
        N=(NX-2)*(NY-2)*(NZ-2)
        !NZA=N*7
        NNX=NX-4
        NNY=NY-4
        NNZ=NZ-4
        NZA=8*4 + 4*(NNX+NNY+NNZ)*5 + 2*(NNX*NNY+NNY*NNZ+NNZ*NNX)*6
     '    + NNX*NNY*NNZ*7
        CALL ALLOCATE_MEMORY(NZA,1,INTTYPE,ISC_PTR,NOINIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(N+1,1,INTTYPE,ISR_PTR,NOINIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(NZA,N,DPTYPE,A_PTR,NOINIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(N,N,DPTYPE,B_PTR,NOINIT,ERROR,*9999)

        CALL SET_FD3D_EQUATIONS(%VAL(A_PTR),LDA,N,%VAL(B_PTR),NX,NY,NZ,
     '    MTYPE,%VAL(ISC_PTR),%VAL(ISR_PTR),NZA,TRANSIENT,ERROR,*9999)

C     2D Finite Differences
      ELSE IF(PTYPE.EQ.3) THEN
        WRITE(*,7300) 'Equations generated on a ',NX,' x ',NY,
     '    ' FDiff mesh'
        N=(NX-2)*(NY-2)
        !NZA=N*5
        NNX=NX-4
        NNY=NY-4
        NZA=4*3 + 2*(NNX+NNY)*4 + NNX*NNY*5
        CALL ALLOCATE_MEMORY(NZA,1,INTTYPE,ISC_PTR,NOINIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(N+1,1,INTTYPE,ISR_PTR,NOINIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(NZA,N,DPTYPE,A_PTR,NOINIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(N,N,DPTYPE,B_PTR,NOINIT,ERROR,*9999)

        CALL SET_FD2D_EQUATIONS(%VAL(A_PTR),LDA,N,%VAL(B_PTR),NX,NY,
     '    MTYPE,%VAL(ISC_PTR),%VAL(ISR_PTR),NZA,TRANSIENT,ERROR,*9999)

C     3D Finite Volumes
      ELSE IF(PTYPE.EQ.4) THEN
        WRITE(*,7300) 'Equations generated on a ',NX,' x ',NY,
     '    ' x ',NZ,' FVol mesh'
        N=(NX-2)*(NY-2)*(NZ-2)
        !NZA=N*27
        NNX=NX-4
        NNY=NY-4
        NNZ=NZ-4
        NZA = 8*8 + 4*(NNX+NNY+NNZ)*12 + 2*(NNX*NNY+NNY*NNZ+NNZ*NNX)*18
     '    + NNX*NNY*NNZ*27
        CALL ALLOCATE_MEMORY(NZA,1,INTTYPE,ISC_PTR,NOINIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(N+1,1,INTTYPE,ISR_PTR,NOINIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(NZA,N,DPTYPE,A_PTR,NOINIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(N,N,DPTYPE,B_PTR,NOINIT,ERROR,*9999)

        CALL SET_FV3D_EQUATIONS(%VAL(A_PTR),LDA,N,%VAL(B_PTR),NX,NY,NZ,
     '    MTYPE,%VAL(ISC_PTR),%VAL(ISR_PTR),NZA,TRANSIENT,ERROR,*9999)

C     2D Finite Volumes
      ELSE IF(PTYPE.EQ.5) THEN
        WRITE(*,7300) 'Equations generated on a ',NX,' x ',NY,
     '    ' FVol mesh'
        N=(NX-2)*(NY-2)
        !NZA=N*9
        NNX=NX-4
        NNY=NY-4
        NZA=4*4 + 2*(NNX+NNY)*6 + NNX*NNY*9
        CALL ALLOCATE_MEMORY(NZA,1,INTTYPE,ISC_PTR,NOINIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(N+1,1,INTTYPE,ISR_PTR,NOINIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(NZA,N,DPTYPE,A_PTR,NOINIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(N,N,DPTYPE,B_PTR,NOINIT,ERROR,*9999)

        CALL SET_FV2D_EQUATIONS(%VAL(A_PTR),LDA,N,%VAL(B_PTR),NX,NY,
     '    MTYPE,%VAL(ISC_PTR),%VAL(ISR_PTR),NZA,TRANSIENT,ERROR,*9999)

C     2D Finite Volumes with boundarys
      ELSE IF(PTYPE.EQ.6) THEN
        WRITE(*,7300) 'Equations generated on a ',NX,' x ',NY,
     '    ' FVol(b) mesh'
        N=NX*NY
        NZA=N*9
        CALL ALLOCATE_MEMORY(NZA,1,INTTYPE,ISC_PTR,NOINIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(N+1,1,INTTYPE,ISR_PTR,NOINIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(NZA,N,DPTYPE,A_PTR,NOINIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(N,N,DPTYPE,B_PTR,NOINIT,ERROR,*9999)

        CALL SET_FV2D_EQUATIONS_BOUNDS(%VAL(A_PTR),LDA,N,%VAL(B_PTR),
     '    NX,NY,MTYPE,%VAL(ISC_PTR),%VAL(ISR_PTR),NZA,TRANSIENT,ERROR,
     '    *9999)
      ELSE
        ERROR='>>Unknown problem type'
        GOTO 9999
      ENDIF
      LDA=N
      SPARSE_A=1
      CALL ALLOCATE_MEMORY( N,1,DPTYPE,X_PTR,NOINIT,ERROR,*9999)
      CALL PDLOAD(N,0.0D0,%VAL(X_PTR),1)
      WRITE(*,'()')

C     If we are to zero the system, do so
      IF(SETZERO.GT.0) THEN
        WRITE(*,7500) 'Zeroing system'
        CALL SET_ZERO(%VAL(A_PTR),LDA,N,%VAL(B_PTR),SPARSE_A,
     '    %VAL(ISC_PTR),%VAL(ISR_PTR),SETZERO,ERROR,*9999)
      ENDIF

C     Possibly zero the Right Hand Side
      IF(ZERO_B) THEN
        WRITE(*,7500) 'Zeroing the RHS'
        CALL PDLOAD(N,1.0D0,%VAL(B_PTR),1)
      ENDIF

C     Convert to a dense array
      IF(DENSE) THEN
        CALL ALLOCATE_MEMORY(N*N,1,DPTYPE,D_PTR,NOINIT,ERROR,*9999)

        CALL SPARSE_TO_DENSE(%VAL(A_PTR),%VAL(D_PTR),LDA,N,SPARSE_A,
     '    %VAL(ISC_PTR),%VAL(ISR_PTR),ERROR,*9999)

        IF(  A_PTR.NE.0) CALL FREE_MEMORY(A_PTR,ERROR,*9999)
c       IF(ISC_PTR.NE.0) CALL FREE_MEMORY(ISC_PTR,ERROR,*9999)
c       IF(ISR_PTR.NE.0) CALL FREE_MEMORY(ISR_PTR,ERROR,*9999)
        A_PTR=D_PTR
        D_PTR=0
        NZA=N**2

        SPARSE_A=0
      ENDIF

C     Randomly weight the equations
      IF(LRAND_WEIGHTS) THEN
        CALL RAND_WEIGHT(%VAL(A_PTR),LDA,N,%VAL(B_PTR),SPARSE_A,
     '    %VAL(ISC_PTR),%VAL(ISR_PTR),ERROR,*9999)
      ENDIF

C     Mix up array ordering
      IF(LMIX_ARRAYS) THEN
        CALL MIX_ARRAYS(%VAL(A_PTR),LDA,N,NZA,%VAL(B_PTR),SPARSE_A,
     '    %VAL(ISC_PTR),%VAL(ISR_PTR),ERROR,*9999)
      ENDIF

C     Experimental reordering schemes.
      IF(KLUDGE.GT.0) THEN
        CALL KLUDGE_PIV(%VAL(A_PTR),LDA,N,NZA,%VAL(B_PTR),SPARSE_A,
     '    %VAL(ISC_PTR),%VAL(ISR_PTR),KLUDGE,ERROR,*9999)
      ENDIF

C     Check input data
      IF(LEXAMINE) THEN
C       CALL KLUDGE_INDEX(%VAL(A_PTR),LDA,N,%VAL(B_PTR),SPARSE_A,
C    '    %VAL(ISC_PTR),%VAL(ISR_PTR),ERROR,*9999)
        CALL Examine(%VAL(A_PTR),LDA,N,%VAL(B_PTR),SPARSE_A,
     '    %VAL(ISC_PTR),%VAL(ISR_PTR),ERROR,*9999)
      ENDIF

C     Check diagonal dominance
      IF(LCHECK_DIAGDOM) THEN
        CALL CHECK_DIAG_DOM(%VAL(A_PTR),LDA,N,SPARSE_A,%VAL(ISC_PTR),
     '    %VAL(ISR_PTR),OUTPUTCODE,ERROR,*9999)
      ENDIF

C     Check symmetry
      IF(LCHECK_IS_SYMM) THEN
        CALL CHECK_SYMM(%VAL(A_PTR),LDA,N,NZA,SPARSE_A,
     '    %VAL(ISC_PTR),%VAL(ISR_PTR),OUTPUTCODE,ERROR,*9999)
      ENDIF

C     Check symmetric structure
      IF(LCHECK_SYMM_STRUCT) THEN
        CALL CHECK_SYMM_FILL(%VAL(A_PTR),LDA,N,NZA,SPARSE_A,
     '    %VAL(ISC_PTR),%VAL(ISR_PTR),OUTPUTCODE,ERROR,*9999)
      ENDIF

C     Check fill pattern
      IF(LCHECK_FILL) THEN
        CALL CHECK_FILL(%VAL(A_PTR),LDA,N,SPARSE_A,
     '     %VAL(ISC_PTR),%VAL(ISR_PTR),OUTPUTCODE,ERROR,*9999)
      ENDIF

C     Check fill pattern
      IF(LCHECK_MATRIX) THEN
        CALL CHECK_MATRIX(%VAL(A_PTR),LDA,N,SPARSE_A,
     '     %VAL(ISC_PTR),%VAL(ISR_PTR),ERROR,*9999)
      ENDIF

C     Dump output file
      IF(DUMP_FILE) THEN
        WRITE(*,7400) 'Dumping equations to ''',
     '    DUMPNAME(:LEN_TRIM(DUMPNAME)),''''
        CALL WRITE_SPARSE(%VAL(A_PTR),LDA,N,NZA,%VAL(B_PTR),SPARSE_A,
     '    %VAL(ISC_PTR),%VAL(ISR_PTR),DUMPNAME,ERROR,*9999)
        STOP
      ENDIF

C     Write triplet file
      IF(WRITE_FILE) THEN
        WRITE(*,7400) 'Writeing equations to triplet file ''',
     '    WRITENAME(:LEN_TRIM(WRITENAME)),''''
        CALL WRITE_TRIPLET(%VAL(A_PTR),LDA,N,NZA,%VAL(B_PTR),SPARSE_A,
     '    %VAL(ISC_PTR),%VAL(ISR_PTR),WRITENAME,ERROR,*9999)
        STOP
      ENDIF

C     Write Boeing file
      IF(WRITE_BOEING_FILE) THEN
        WRITE(*,7400) 'Writeing equations to Boeing file ''',
     '    BOEINGNAME(:LEN_TRIM(BOEINGNAME)),''''
         CALL WRITE_BOEING(%VAL(A_PTR),N,NZA,SPARSE_A,%VAL(ISC_PTR),
     '     %VAL(ISR_PTR),BOEINGNAME,ERROR,*9999)
         STOP
      ENDIF

C     Set the maximum time
      IERR=ALARM(MAX(0,MAXTIME),ALARM_WAKE)

C     Solvers

C     Test the basic operations
      IF(SOLVER.EQ.30) THEN

        CALL PDLOAD(N,0.0D0,%VAL(B_PTR),1)

C       Matrix Multiply
        CALL TIMER_START(1)
        CALL TIMER_START(2)
        DO I=1,N_FACT
          CALL PDLOAD(N,1.0D0,%VAL(X_PTR),1)
          CALL MATRIX_MULT(%VAL(A_PTR),LDA,N,%VAL(X_PTR),%VAL(B_PTR),
     '      SPARSE_A,%VAL(ISC_PTR),%VAL(ISR_PTR),ERROR,*9999)
        ENDDO
        CALL TIMER_STOP(2)

        ITMP=SPARSE_A
        IF(LTEST) ITMP=0

        CALL TIMER_START(3)
        DO I=1,N_SOLV
          RDUM=DDOT(N,%VAL(X_PTR),1,%VAL(B_PTR),1)
        ENDDO
        CALL TIMER_STOP(3)
        CALL TIMER_STOP(1)

C     Dense code
      ELSE IF(DENSE .AND. (SOLVER.EQ.SOLV_LU .OR. SOLVER.EQ.SOLV_SVD  
     '    .OR. SOLVER.EQ.SOLV_LSQ .OR. SOLVER.EQ.SOLV_CHOLESKY .OR. 
     '    SOLVER.EQ.SOLV_LDL)) THEN

C       Alloc memory
        M=N
        CALL ALLOCATE_MEMORY(LDA*N,1,DPTYPE,D_PTR,NOINIT,ERROR,*9999)
        CALL DIR_ALLOC(LDA,M,N,SOLVER,E_PTR,F_PTR,I_PTR,E_LEN,
     '    F_LEN,ERROR,*9999)

C       Factorise
        CALL TIMER_START(1)
        CALL TIMER_START(2)
        DO I=1,N_FACT
          CALL DCOPY(LDA*N,%VAL(A_PTR),1,%VAL(D_PTR),1)
          CALL DIR_FACTOR(%VAL(D_PTR),LDA,M,N,%VAL(I_PTR),SOLVER,
     '      %VAL(E_PTR),%VAL(F_PTR),F_LEN,ANORM,ERROR,*9999)
        ENDDO
        CALL TIMER_STOP(2)

C       Solve the problem
        CALL TIMER_START(3)
        DO I=1,N_SOLV
          CALL PDLOAD(N,0.0D0,%VAL(X_PTR),1)
          CALL DIR_SOLVE(%VAL(D_PTR),LDA,LDA,M,N,%VAL(X_PTR),
     '      %VAL(B_PTR),%VAL(I_PTR),SOLVER,%VAL(E_PTR),
     '      %VAL(F_PTR),F_LEN,ANORM,OUTPUTCODE,ERROR,*9999)
        ENDDO
        CALL TIMER_STOP(3)
        CALL TIMER_STOP(1)

C     Sparse direct solvers
C     Harwell MA28AD
      ELSE IF(SOLVER.EQ.SOLV_HARWELL) THEN

        IF(OUTPUTCODE.GT.0) THEN
          WRITE(*,'(A)') ' Harwell MA28AD Sparse LU Decomposition'
        ENDIF
        ERROR='>>Harwell MA28AD solver no longer linked to'
        GOTO 9999

C       Malloc
C       USET=0.1D0
C       CALL MA28AD_ALLOC(N,NZA,NICN,NIRN,LICN,LIRN,D_PTR,E_PTR,
C    '    I_PTR,J_PTR,ERROR,*9999)

C       Factorise
C       CALL TIMER_START(1)
C       CALL TIMER_START(2)

C       DO I=1,N_FACT
C         CALL MA28AD_FACTOR(%VAL(A_PTR),%VAL(D_PTR),N,NZA,LICN,LIRN,
C    '      %VAL(ISC_PTR),%VAL(ISR_PTR),USET,%VAL(I_PTR),%VAL(J_PTR),
C    '      %VAL(E_PTR),ERROR,*9999)
C       ENDDO
C       CALL TIMER_STOP(2)

C       If flagged to do so, dump out the factorised system for a look.
C       IF(DUMP_FACT) THEN
C         CALL WRITE_MA28AD_FACTORED(%VAL(D_PTR),LDA,N,%VAL(B_PTR),
C    '      %VAL(ISR_PTR),%VAL(IC_LU_PTR),%VAL(IR_LU_PTR),
C    '      %VAL(IKEEP_PTR),FACTNAME,ERROR,*9999)
C       ENDIF

C       Solve the problem
C       CALL TIMER_START(3)
C       DO I=1,N_SOLV
C         CALL MA28AD_SOLVE(%VAL(D_PTR),%VAL(B_PTR),%VAL(X_PTR),N,NZA,
C    '      LICN,%VAL(I_PTR),%VAL(J_PTR),%VAL(E_PTR),ERROR,*9999)
C       ENDDO
C       CALL TIMER_STOP(3)
C       CALL TIMER_STOP(1)

C       Free up memory
C       CALL MA28AD_FREE(D_PTR,E_PTR,I_PTR,J_PTR,ERROR,*9999)

C     Sparse direct solvers
C     UMFPACK/MA38AD
      ELSE IF(SOLVER.EQ.SOLV_UMFPACK2) THEN

        IF(OUTPUTCODE.GT.0) THEN
          WRITE(*,'(A)') ' UMFPACK-2.2 Sparse LU Decomposition'
        ENDIF
        ERROR='>>UMFPACK-2.2 solver no longer linked to'
        GOTO 9999

C       Malloc and initialise
C       CALL UMD2FA_ALLOC(N,NZA,NUMF,NINDX,NDARR,E_PTR,I_PTR,
C    '    ERROR,*9999)
C       IF(UMF_ITER) JOB=1

C       Factorise
C       CALL TIMER_START(1)
C       CALL TIMER_START(2)

C       DO I=1,N_FACT
C         CALL UMD2FA_FACTOR(%VAL(A_PTR),N,NZA,%VAL(ISC_PTR),
C    '      %VAL(ISR_PTR),%VAL(E_PTR),%VAL(I_PTR),NINDX,NDARR,JOB,
C    '      KEEP,ICNTL,INFO,CNTL,RINFO,UMFDEF,OUTPUTCODE,ERROR,*9999)
C       ENDDO
        CALL TIMER_STOP(2)

C       Solve the problem
C       CALL TIMER_START(3)
C       DO I=1,N_SOLV
C         CALL UMD2FA_SOLVE(%VAL(A_PTR),%VAL(B_PTR),%VAL(X_PTR),N,NZA,
C    '      %VAL(E_PTR),%VAL(I_PTR),NINDX,NDARR,KEEP,ICNTL,INFO,
C    '      CNTL,RINFO,ERROR,*9999)
C       ENDDO
C       CALL TIMER_STOP(3)
C       CALL TIMER_STOP(1)

C     Sparse direct solvers
C     UMFPACK-4.0
      ELSE IF(SOLVER.EQ.SOLV_UMFPACK4) THEN

        IF(OUTPUTCODE.GT.0) THEN
          WRITE(*,'(A)') ' UMFPACK-4.0 Sparse LU Decomposition'
        ENDIF

C       Factorise
        CALL TIMER_START(1)
        CALL TIMER_START(2)
        DO I=1,N_FACT
          IF(OUTPUTCODE.GT.0) THEN
            WRITE(*,'(A)') ' Factor'
          ENDIF
          IF(I.NE.1) THEN
            CALL UMFPACK4_FREE(I_PTR,.TRUE.,D_PTR,.TRUE.,ERROR,*9999)
          ENDIF
          CALL UMFPACK4_FACTOR(%VAL(A_PTR),N,NZA,%VAL(ISC_PTR),
     '      %VAL(ISR_PTR),I_PTR,D_PTR,UMFDEF,UMF_CONTROL,UMF_INFO,ANORM,
     '      ERROR,*9999)
        ENDDO

C       Refactorise
        DO I=1,N_RFAC
          IF(OUTPUTCODE.GT.0) THEN
            WRITE(*,'(A)') ' Re-Factor'
          ENDIF

          CALL UMFPACK4_FREE(I_PTR,.FALSE.,D_PTR,.TRUE.,ERROR,*9999)
          CALL UMFPACK4_REFACTOR(%VAL(A_PTR),N,NZA,%VAL(ISC_PTR),
     '      %VAL(ISR_PTR),I_PTR,D_PTR,UMF_CONTROL,UMF_INFO,ANORM,
     '      ERROR,*9999)
        ENDDO
        CALL TIMER_STOP(2)

C       Solve the problem
        CALL TIMER_START(3)
        DO I=1,N_SOLV
          IF(OUTPUTCODE.GT.0) THEN
            WRITE(*,'(A)') ' Solve'
          ENDIF
          CALL UMFPACK4_SOLVE(%VAL(A_PTR),N,NZA,%VAL(ISC_PTR),
     '      %VAL(ISR_PTR),%VAL(B_PTR),%VAL(X_PTR),D_PTR,UMF_CONTROL,
     '      UMF_INFO,ANORM,OUTPUTCODE,ERROR,*9999)
        ENDDO
        CALL TIMER_STOP(3)
        CALL TIMER_STOP(1)

        CALL UMFPACK4_FREE(I_PTR,.TRUE.,D_PTR,.TRUE.,ERROR,*9999)

C     Sparse direct solvers
C     SuperLU
      ELSE IF(SOLVER.EQ.SOLV_SUPERLU) THEN

        IF(OUTPUTCODE.GT.0) THEN
          WRITE(*,'(A)') ' SuperLU Sparse LU Decomposition'
        ENDIF
C       Malloc and initialise
        CALL SUPERLU_RESET(PARAM,ERROR,*9999)
        CALL SUPERLU_SETPARAM('COLUMN ORDERING',COL_ORDER,PARAM,ERROR,
     '    *9999)
        CALL SUPERLU_SETPARAM('NUMBER OF PROCS',NCPU,PARAM,ERROR,*9999)
        CALL SUPERLU_SETPARAM('ESTIMATED FILL',SLU_FILL,PARAM,
     '    ERROR,*9999)
        IF(SLU_RELAX.NE.0) THEN
          CALL SUPERLU_SETPARAM('RELAX PARAMETER',SLU_RELAX,PARAM,
     '      ERROR,*9999)
        ENDIF
c       CALL SUPERLU_SETPARAM('ESTIMATED U FILL',-260,PARAM,
c    '    ERROR,*9999)
c       CALL SUPERLU_SETPARAM('ESTIMATED L SUB',-260,PARAM,
c    '    ERROR,*9999)

C       Factorise
        CALL TIMER_START(1)
        CALL TIMER_START(2)
        DO I=1,N_FACT
          IF(OUTPUTCODE.GT.0) THEN
            WRITE(*,'(A)') ' Factor'
          ENDIF
          IF(I.NE.1) THEN
            CALL SUPERLU_FREE(D_PTR,ERROR,*9999)
          ENDIF
          CALL SUPERLU_FACTOR(%VAL(A_PTR),N,NZA,%VAL(ISC_PTR),
     '      %VAL(ISR_PTR),D_PTR,PARAM,ANORM,ERROR,*9999)
        ENDDO

C       Refactorise
        DO I=1,N_RFAC
          IF(OUTPUTCODE.GT.0) THEN
            WRITE(*,'(A)') ' Re-Factor'
          ENDIF

          CALL SUPERLU_REFACTOR(%VAL(A_PTR),N,NZA,%VAL(ISC_PTR),
     '      %VAL(ISR_PTR),D_PTR,PARAM,ANORM,ERROR,*9999)
        ENDDO
        CALL TIMER_STOP(2)

C       Solve the problem
        CALL TIMER_START(3)
        DO I=1,N_SOLV
          IF(OUTPUTCODE.GT.0) THEN
            WRITE(*,'(A)') ' Solve'
          ENDIF
          CALL SUPERLU_SOLVE(%VAL(B_PTR),%VAL(X_PTR),N,NZA,D_PTR,
     '      PARAM,ANORM,OUTPUTCODE,ERROR,*9999)
        ENDDO
        CALL TIMER_STOP(3)
        CALL TIMER_STOP(1)

        CALL SUPERLU_FREE(D_PTR,ERROR,*9999)

C     Iterative methods
      ELSE IF(SOLVER.NE.0) THEN

C       Allocate arrays to hold the factorisation
        CALL ITER_ALLOC(LDA,N,NZA,SOLVER,PRECON,SPARSE_A,D_PTR,
     '    ERROR,*9999)

C       Factorise
        CALL TIMER_START(1)
        CALL TIMER_START(2)
        DO I=1,N_FACT
          CALL ITER_FACTOR(%VAL(A_PTR),LDA,N,%VAL(D_PTR),ANORM,SPARSE_A,
     '      %VAL(ISC_PTR),%VAL(ISR_PTR),SOLVER,PRECON,OUTPUTCODE,ERROR,
     '      *9999)
        ENDDO
        CALL TIMER_STOP(2)

C       Solve the problem
        ITMP=OUTPUTCODE
        CALL TIMER_START(3)
        DO I=1,N_SOLV
          CALL PDLOAD(N,0.0D0,%VAL(X_PTR),1)
          ITER_TMP=ITER
          RESID_TMP=RESID
          CALL ITER_SOLVE(%VAL(A_PTR),LDA,N,%VAL(X_PTR),%VAL(B_PTR),
     '      %VAL(D_PTR),SPARSE_A,%VAL(ISC_PTR),%VAL(ISR_PTR),ANORM,
     '      RESID_TMP,OMEGA,ITER_TMP,NRES,NPRECON,SOLVER,PRECON,
     '      OUTPUTCODE,ERROR,*9999)
          OUTPUTCODE=0
        ENDDO
        CALL TIMER_STOP(3)
        CALL TIMER_STOP(1)
        OUTPUTCODE=ITMP
        ITER=ITER_TMP
        RESID=RESID_TMP

        IF(OUTPUTCODE.EQ.0) WRITE(*,'(A,/,I5,1x,1P,G12.6)') 
     '    'Iterations/Residual',ITER,RESID

      ENDIF

C     Check the solution
      CALL CHECK_SOLN(%VAL(A_PTR),LDA,N,%VAL(B_PTR),%VAL(X_PTR),
     '  SPARSE_A,%VAL(ISC_PTR),%VAL(ISR_PTR),ERROR,*9999)

C     Print the times
      DO I=1,3
        TIME(I,1,1)=TIMER_READ_BUSY(I)
        TIME(I,2,1)=TIMER_READ_WALL(I)
      ENDDO

      DO I=1,2
        TIME(2,I,2)=TIME(2,I,1)/DBLE(N_FACT)
        TIME(3,I,2)=TIME(3,I,1)/DBLE(N_SOLV)
        TIME(1,I,2)=TIME(2,I,2)+TIME(3,I,2)
      ENDDO

      WRITE(*,7200)   ' Number of Equations : ',N
      IF(N_FACT.GT.1 .OR. N_SOLV.GT.1) THEN
        WRITE(*,7100) ' Factorisation time  :',TIME(2,1,1),TIME(2,2,1),
     '    ' :',TIME(2,1,2),TIME(2,2,2)
        WRITE(*,7100) ' Solution time       :',TIME(3,1,1),TIME(3,2,1),
     '    ' :',TIME(3,1,2),TIME(3,2,2)
        WRITE(*,7100) ' Total time          :',TIME(1,1,1),TIME(1,2,1),
     '    ' :',TIME(1,1,2),TIME(1,2,2)
      ELSE
        WRITE(*,7000) ' Factorisation time  :',TIME(2,1,1),TIME(2,2,1)
        WRITE(*,7000) ' Solution time       :',TIME(3,1,1),TIME(3,2,1)
        WRITE(*,7000) ' Total time          :',TIME(1,1,1),TIME(1,2,1)
      ENDIF

C     Possibly write out a solution
      IF(SOLN2D .AND. (PTYPE.EQ.5 .OR. PTYPE.EQ.3)) THEN
        CALL WRITE_SOLN2D(%VAL(X_PTR),N,NX,NY,SOLNNAME,ERROR,*9999)
      ENDIF

C     Free up memory (in case we routinise the code)
      IF(A_PTR.NE.0) CALL FREE_MEMORY(A_PTR,ERROR,*9999)
      IF(B_PTR.NE.0) CALL FREE_MEMORY(B_PTR,ERROR,*9999)
      IF(D_PTR.NE.0) CALL FREE_MEMORY(D_PTR,ERROR,*9999)
      IF(E_PTR.NE.0) CALL FREE_MEMORY(E_PTR,ERROR,*9999)
      IF(F_PTR.NE.0) CALL FREE_MEMORY(F_PTR,ERROR,*9999)
      IF(X_PTR.NE.0) CALL FREE_MEMORY(X_PTR,ERROR,*9999)

      IF(I_PTR.NE.0) CALL FREE_MEMORY(I_PTR,ERROR,*9999)
      IF(J_PTR.NE.0) CALL FREE_MEMORY(J_PTR,ERROR,*9999)

      CALL EXIT(0)
      STOP

 7000 FORMAT(A,2(1X,F12.6))
 7100 FORMAT(2(A,2(1X,F12.6)))
 7200 FORMAT(A,I7)
 7300 FORMAT(A,3(I5,A))
 7400 FORMAT(A,A,A)
 7500 FORMAT(A)

 9999 CALL ERRORS('MAIN',ERROR)
      CALL EXIT(-1)
      END


      SUBROUTINE ALARM_WAKE()

      IMPLICIT NONE
      INTEGER MAXTIME

      COMMON / ALARM_COMMON / MAXTIME

      WRITE(0,7100) ' Maximum Time exceeded'
      WRITE(0,7000) ' alarm_wake() called after ',MAXTIME,' seconds'
      WRITE(0,7100) ' soltest exiting ...'

      CALL EXIT(-1)
 7000 FORMAT(A,I9,A)
 7100 FORMAT(A)
      STOP
      END


      SUBROUTINE USAGE()

      IMPLICIT NONE

      WRITE(*,'()')
      WRITE(*,*) 'Usage: soltest [-file <input datafile>|'
      WRITE(*,*) '       |-f[dv]2D <NX,NY>|-f[dv] <NX,NY,NZ>|]'
      WRITE(*,*) '       |-triplet <input triplet file>]'
      WRITE(*,*) '       [-trans|-steady|-notrans]'
      WRITE(*,*) '       [-solver <solver>] [-precon <precon>]'
      WRITE(*,*) '       [-nfact <N>] [-nsolv <N>] [-verbose|-quiet]'
      WRITE(*,*) '       [-niter <N>] [-tol <R>] [-omega <R>]'
      WRITE(*,*) '       [-dump_file <output datafile>] [-sinmesh]'
      WRITE(*,*) '       [-dump_fact <factored datafile>] [-nosolve]'
      WRITE(*,*) '       [-write_file <output triplet ile>]'
      WRITE(*,*) '       [-soln2d_file <output 2D gnuplot file>]'
      WRITE(*,*) '       [-har_narr <NICN,NIRN>] [-slu_co <col_order>]'
      WRITE(*,*) '       [-slu_fill <fill>]'
      WRITE(*,*) '       [-umf_narr <NUMF>] [-umf_iter|-umf_noiter]'
      WRITE(*,*) '       [-check_symm] [-check_fill] [-check_diagdom]'
      WRITE(*,*) '       [-mix_arrays] [-rand_weights]'
      WRITE(*,*) '       [-examine] [-check_symm_struct]'
      WRITE(*,*) '       [-nprecon <N>] [-nres <N>] [-maxtime <Time>]'
      WRITE(*,*) '       [-setzero[:N]] [-zerob]'
      WRITE(*,'()')
      WRITE(*,*) '       solvers: jacobi sor ic[0|1] ilu[0|1] ildl[0|1]'
      WRITE(*,*) '                cg bicgstab gmres harwell umfpack'
      WRITE(*,*) '                superlu none matmul'
      WRITE(*,*)
      WRITE(*,*) '       preconditioners: none point-jacobi jacobi ssor'
      WRITE(*,*) '                ic[0|1] ilu[0|1] ildl[0|1] rowscale'
      WRITE(*,'()')

      CALL EXIT(1)
      STOP
      END
