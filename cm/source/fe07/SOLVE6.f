      SUBROUTINE SOLVE6(ISC_GK,ISC_GKK,ISC_GMM,ISC_GQ,ISR_GK,ISR_GKK,
     '  ISR_GMM,ISR_GQ,LGE,NBH,NENP,NHE,NONY,NP_INTERFACE,NPNE,NPNY,
     '  nr,NRE,NVHE,nx,NYNE,NYNP,NYNR,CONY,EIGVAL,EIGVEC,GK,GKK,GM,
     '  GMM,GQ,GRR,DYNAM2,FIRST_A,UPDATE_MATRIX,ERROR,*)

C#### Subroutine: SOLVE6
C###  Description:
C###    SOLVE6 solves modal analysis problems.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'eige00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'solv00.cmn'
      INCLUDE 'time02.cmn'
!     Parameter List
      INTEGER ISC_GK(NISC_GKM),ISC_GKK(NISC_GKKM),ISC_GMM(NISC_GMMM),
     '  ISC_GQ(NISC_GQM),ISR_GK(NISR_GKM),ISR_GKK(NISR_GKKM),
     '  ISR_GMM(NISR_GMMM),ISR_GQ(NISR_GQM),LGE(NHM*NSM,NRCM),
     '  NBH(NHM,NCM,NEM),NP_INTERFACE(0:NPM,0:3),
     '  NENP(NPM,0:NEPM,0:NRM),NHE(NEM),NONY(0:NOYM,NYM,NRCM),
     '  NPNE(NNM,NBFM,NEM),NPNY(0:6,NYM,0:NRCM),
     '  nr,NRE(NEM),NVHE(NNM,NBFM,NHM,NEM),nx,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM)
      REAL*8 CONY(0:NOYM,NYM,NRCM),
     '  EIGVAL(NTM,2),EIGVEC(NOM,NTM,2),GK(NZ_GK_M),GKK(NZ_GKK_M),
     '  GM(NZ_GM_M),GMM(NZ_GMM_M),GQ(NZ_GQ_M),GRR(NOM)
      CHARACTER ERROR*(*)
      LOGICAL DYNAM2,FIRST_A,UPDATE_MATRIX
!     Local Variables
      INTEGER GETNYR,no1,no2,no_nynr1,no_nynr2,noy1,
     '  noy2,nt,ny1,ny2,ny3,nz,nzz,PROBLEM_TYPE,SOLVER_TYPE
      INTEGER*4 WORK_PTR
      REAL ELAPSED_TIME,TIME_START1(1),TIME_START2(1),TIME_STOP(1)
      REAL*8 AA,BB,co1,co2

      CALL ENTERS('SOLVE6',*9999)

      PROBLEM_TYPE=2 !Temporary
      SOLVER_TYPE=1 !Temporary

      IF(NOT(2,1,nr,nx).EQ.0) THEN
        ERROR=' >>The number of unknowns is zero'
        GOTO 9999
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_START1)
      CALL CPU_TIMER(CPU_USER,TIME_START2)

C*** Setup and initialise arrays

      IF(UPDATE_MATRIX) THEN
        IF(FIRST_A) THEN
          SPARSEGKK(nx)=0
          WORK_PTR=0
          CALL CALC_SPARSE_GKK(ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,ISR_GKK,
     '      ISR_GQ,LGE,NBH,NENP,NHE,NOT(1,1,nr,nx),NOT(2,1,nr,nx),
     '      NONY,NP_INTERFACE,NPNE,NPNY,nr,NRE,NVHE,nx,NYNE,NYNP,
     '      NYNR,GK,GQ,%VAL(WORK_PTR),.FALSE.,.TRUE.,ERROR,*9999)
          NOT(1,4,nr,nx)=NOT(1,1,nr,nx)
          NOT(2,4,nr,nx)=NOT(2,1,nr,nx)
          NZZT(4,nr,nx)=NZZT(1,nr,nx)
        ENDIF
        CALL ASSERT(NZZT(1,nr,nx).LE.NZ_GKK_M,'>>Increase NZ_GKK_M',
     '    ERROR,*9999)
        CALL ASSERT(NZZT(4,nr,nx).LE.NZ_GMM_M,'>>Increase NZ_GMM_M',
     '    ERROR,*9999)
        DO nzz=1,NZZT(1,nr,nx)
          GKK(nzz)=0.0d0
          GMM(nzz)=0.0d0
        ENDDO !nzz
      ENDIF

C*** Generate reduced system

      DO no_nynr1=1,NYNR(0,1,1) !Loop over global rows of GK
        ny1=NYNR(no_nynr1,1,1) !is row #
        DO noy1=1,NONY(0,ny1,1) !loop over #no's attached to row ny1
          no1=NONY(noy1,ny1,1) !is no# attached to row ny1
          co1=CONY(noy1,ny1,1) !is coupling coeff for row mapping
C                               ie row_no1=a*row_ny1+b*row_ny2
          DO no_nynr2=1,NYNR(0,0,1) !loop over the #cols of GK
            ny2=NYNR(no_nynr2,0,1) !is global variable #
            ny3=GETNYR(1,NPNY,nr,2,0,ny2,NYNE,NYNP) !local GK var #
            CALL SPARSE(ny1,ny3,NYT(1,1,nx),nz,NZ_GK_M,
     '        NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
            IF(nz.NE.0) THEN
              AA=GK(nz)
              IF(DYNAM2) THEN
                BB=GM(nz)
              ELSE
                IF(ny1.EQ.ny3) THEN
                  BB=1.0d0
                ELSE
                  BB=0.0d0
                ENDIF
              ENDIF
              DO noy2=1,NONY(0,ny2,2) !loop over #no's for var ny2
                no2=NONY(noy2,ny2,2) !no# attached to ny2
                co2=CONY(noy2,ny2,2) !coup coeff for the col mappng
C                                     i.e. var_no1=a*var_ny1+b*var_ny2
                CALL SPARSE(no1,no2,NOT(1,1,nr,nx),nzz,NZ_GKK_M,
     '            NZZT(1,nr,nx),ISC_GKK,ISR_GKK,SPARSEGKK(nx),
     '            ERROR,*9999)
                IF(nzz.NE.0.AND.UPDATE_MATRIX) THEN
                  GKK(nzz)=GKK(nzz)+AA*co1*co2
                  GMM(nzz)=GMM(nzz)+BB*co1*co2
                ENDIF
              ENDDO !noy2
            ENDIF
          ENDDO !no_nynr2
        ENDDO !noy1
      ENDDO !no_nynr1

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
      IF(IWRIT4(nr,nx).GE.1) THEN
        WRITE(OP_STRING,'(/'' CPU time for solution matrix '
     '    //'initialisation and assembly: '',D11.4,'' s'')')
     '    ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_START2)
C      IF(KTYP4.NE.0.AND.FIRST_A) THEN !Output global matrices
C        CALL WRITE_SOL_MATRIX(ISC_GKK,ISR_GKK,nr,nx,GKK,GRR,
C     '    ERROR,*9999)
C      ENDIF

      IF(IWRIT4(nr,nx).GE.3.AND.UPDATE_MATRIX) THEN
        WRITE(OP_STRING,'(/'' Global stiffness matrix GKK:'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' NOT(1,1,nr,nx)='',I5,'
     '    //''', NOT(2,1,nr,nx)='',I5)') NOT(1,1,nr,nx),
     '    NOT(2,1,nr,nx)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        CALL OPSTFMAT(NYNR(0,1,1),ISC_GKK,ISR_GKK,IOOP,NOT(1,1,nr,nx),
     '    NOT(2,1,nr,nx),NZZT(1,nr,nx),NYNR(0,2,1),SPARSEGKK(nx),GKK,
     '    GRR,'GKK','GRR',.TRUE.,UPDATE_MATRIX,.FALSE.,ERROR,
     '    *9999)
        WRITE(OP_STRING,'(/'' Global mass matrix GMM:'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' NOT(1,4,nr,nx)='',I5,'
     '    //''', NOT(2,4,nr,nx)='',I5)') NOT(1,4,nr,nx),
     '    NOT(2,4,nr,nx)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        CALL OPSTFMAT(NYNR(0,1,1),ISC_GMM,ISR_GMM,IOOP,NOT(1,4,nr,nx),
     '    NOT(2,4,nr,nx),NZZT(4,nr,nx),NYNR(0,2,1),SPARSEGKK(nx),GMM,
     '    GRR,'GMM','GRR',.TRUE.,UPDATE_MATRIX,.FALSE.,ERROR,
     '    *9999)
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
      IF(IWRIT4(nr,nx).GE.4) THEN
        WRITE(OP_STRING,'(/'' CPU time for solution matrix '
     '    //'output: '',D11.4,'' s'')') ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

C***  Allocate space for work arrays

C      CALL GETEIGENWORK(NOM,PROBLEM_TYPE,SOLVER_TYPE,
C     '  IWORK_PTR,IWORK1_EIGEN_PTR,WORK2_PTR,
C     '  WORK1_EIGEN_PTR,
C     '  WORK2_EIGEN_PTR,WORK3_EIGEN_PTR,WORK4_EIGEN_PTR,
C     '  .TRUE.,FIRST_A,ERROR,*9999)

C***  Solve generalised eigenproblem

      CALL EIGENPROBLEM(NOT(1,1,nr,nx),NOT(2,1,nr,nx),NOT(1,4,nr,nx),
     '  NOT(2,4,nr,nx),NOM,NOT(1,1,nr,nx),NTM,NUMEIGEN,KTYP17,
     '  IWRIT4(nr,nx),PROBLEM_TYPE,SOLVER_TYPE,GKK,GMM,EIGVAL(1,1),
     '  EIGVEC(1,1,1),ERROR,*9999)

      DO nt=1,NUMEIGEN
        EIGVAL(nt,2)=0.0d0 !Real eigenvalues only
        DO no1=1,NOT(1,1,nr,nx)
          EIGVEC(no1,nt,2)=0.0d0
        ENDDO !no
      ENDDO !nt

C***  Free space for work arrays

C      CALL GETEIGENWORK(NOM,PROBLEM_TYPE,SOLVER_TYPE,
C     '  IWORK_PTR,IWORK1_EIGEN_PTR,WORK2_PTR,
C     '  WORK1_EIGEN_PTR,
C     '  WORK2_EIGEN_PTR,WORK3_EIGEN_PTR,WORK4_EIGEN_PTR,
C     '  .FALSE.,FIRST_A,ERROR,*9999)

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START1(1)
      IF(IWRIT4(nr,nx).GE.1) THEN
        WRITE(OP_STRING,'(/'' Total CPU time for solution: '','
     '    //'D11.4,'' s'')') ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('SOLVE6')
      RETURN
C 1112 CALL GETEIGENWORK(NOM,PROBLEM_TYPE,SOLVER_TYPE,
C     '  IWORK_PTR,IWORK1_EIGEN_PTR,WORK2_PTR,
C     '  WORK1_EIGEN_PTR,
C     '  WORK2_EIGEN_PTR,WORK3_EIGEN_PTR,WORK4_EIGEN_PTR,
C     '  .FALSE.,FIRST_A,ERROR,*9999)
 9999 CALL ERRORS('SOLVE6',ERROR)
      CALL EXITS('SOLVE6')
      RETURN 1
      END


