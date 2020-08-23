      SUBROUTINE SOLVE5(ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,ISR_GKK,ISR_GQ,
     '  LGE,NBH,NENP,NHE,NONY,NP_INTERFACE,NPNE,NPNY,nr,NRE,NVHE,
     '  nx,NYNE,NYNO,NYNP,NYNR,CONY,CYNO,GK,GKK,GQ,GR,GRR,
     '  XO,YP,FIRSTS,UPDATE_MATRIX,ERROR,*)

C#### Subroutine: SOLVE5
C###  Description:
C###    Solves the resulting system of linear equations for the
C###    increments of the nonlinear solver in NONLIN.

      IMPLICIT NONE
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'solv00.cmn'
      INCLUDE 'time02.cmn'
!     Parameter List
      INTEGER ISC_GK(NISC_GKM),ISC_GKK(NISC_GKKM),ISC_GQ(NISC_GQM),
     '  ISR_GK(NISR_GKM),ISR_GKK(NISR_GKKM),ISR_GQ(NISR_GQM),
     '  LGE(NHM*NSM,NRCM),NBH(NHM,NCM,NEM),NENP(NPM,0:NEPM,0:NRM),
     '  NHE(NEM),NONY(0:NOYM,NYM,NRCM),NP_INTERFACE(0:NPM,0:3),
     '  NPNE(NNM,NBFM,NEM),NPNY(0:6,NYM,0:NRCM),nr,NRE(NEM),
     '  NVHE(NNM,NBFM,NHM,NEM),nx,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),NYNO(0:NYOM,NOOPM,NRCM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),NYNR(0:NY_R_M,0:NRCM,NCM)
      REAL*8 CONY(0:NOYM,NYM,NRCM),CYNO(0:NYOM,NOOPM,NRCM),
     '  GK(NZ_GK_M),GKK(NZ_GKK_M),GQ(NZ_GQ_M),GR(NYROWM),GRR(NOM),
     '  XO(NOM),YP(NYM,NIYM)
      LOGICAL FIRSTS,UPDATE_MATRIX
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER no1,no2,no_nynr1,
     '  no_nynr2,noy1,noy2,
     '  ny1,ny2,ny3,nyo2,nz,nzz
      INTEGER*4 WORK_PTR
      REAL ELAPSED_TIME,TIME_START(1),TIME_STOP(1)
      REAL*8 co1,co2,DIRDERIV
      LOGICAL ERROR_FLAG
      CHARACTER*80 PRIVATE_ERROR
!     External Functions
      INTEGER GETNYR,LEN_TRIM

      CALL ENTERS('SOLVE5',*9999)

      CALL CPU_TIMER(CPU_USER,TIME_START)


C*** Initialise arrays
c cpb 23/9/95 New solution system
      IF(UPDATE_MATRIX) THEN
        IF(FIRSTS) THEN
C CPB 6/11/95 Temporary work array allocation
          IF(SPARSEGKK(nx).NE.0) THEN
            WORK_PTR=0
            CALL ALLOCATE_MEMORY(NOT(1,1,nr,nx)*NOT(2,1,nr,nx),1,
     '        CHARTYPE,WORK_PTR,MEM_INIT,ERROR,*9999)
          ENDIF

          CALL CALC_SPARSE_GKK(ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,ISR_GKK,
     '      ISR_GQ,LGE,NBH,NENP,NHE,NOT(1,1,nr,nx),NOT(2,1,nr,nx),
     '      NONY,NP_INTERFACE,NPNE,NPNY,nr,NRE,NVHE,nx,NYNE,NYNP,
     '      NYNR,GK,GQ,%VAL(WORK_PTR),.FALSE.,.FALSE.,ERROR,*9999)
          IF(SPARSEGKK(nx).NE.0) CALL FREE_MEMORY(WORK_PTR,ERROR,*9999)

        ENDIF
        DO nzz=1,NZZT(1,nr,nx)
          GKK(nzz)=0.0d0
        ENDDO !nzz
      ENDIF
      DO no1=1,NOT(1,1,nr,nx)
        GRR(no1)=0.0d0
      ENDDO !no1

C *** Calculate global RHS vector

      DO no_nynr1=1,NYNR(0,1,1) !loop over rows
        ny1=NYNR(no_nynr1,1,1) !is row number
        GR(ny1)=YP(ny1,4)
      ENDDO

      ERROR_FLAG=.FALSE.

C KAT 2003-04-22: YP(ny3,5) is not initialized, and the current solution
C was adjusted for increments in NONLIN before the residual was
C calculated, so YP(1,4) already contains the appropriate adjustments
C C *** Calculate global RHS vector using current soln increments from YP
C 
C C$OMP PARALLEL DO DEFAULT(NONE)
C C$OMP&  PRIVATE(no_nynr1,no_nynr2,ny1,ny2,ny3,nz,PRIVATE_ERROR)
C C$OMP&  SHARED(ERROR_FLAG,FIX,NYNR,nx,ISC_GK,ISR_GK,GK,GR,KTYP24,
C C$OMP&  NYT,NZ_GK_M,NZT,YP)
C       DO no_nynr2=1,NYNR(0,2,1) !loop over local columns of GK
C         ny3=NYNR(no_nynr2,0,1) !global var # for local var ny2
C         IF(FIX(ny3,1).AND.YP(ny3,5).NE.0.0d0) THEN
C           ny2=NYNR(no_nynr2,2,1) !local column #
C           DO no_nynr1=1,NYNR(0,1,1) !loop over rows
C             ny1=NYNR(no_nynr1,1,1) !is row number
C             CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,NZ_GK_M,
C      '        NZT(1,nx),ISC_GK,ISR_GK,KTYP24,PRIVATE_ERROR,*100)
C             IF(nz.NE.0) THEN
C C$OMP         ATOMIC
C               GR(ny1)=GR(ny1)+GK(nz)*YP(ny3,5)
C             ENDIF
C           ENDDO
C         ENDIF
C         GO TO 105
C  100      CONTINUE
C           ERROR_FLAG=.TRUE.
C           IF(PRIVATE_ERROR.NE.' ') THEN
C             CALL FLAG_ERROR(0,
C      '        PRIVATE_ERROR(:LEN_TRIM(PRIVATE_ERROR)))
C           ENDIF
C  105    CONTINUE
C       ENDDO
C C$OMP END PARALLEL DO
C       IF(ERROR_FLAG) THEN
C         ERROR=' '
C         GOTO 9999
C       ENDIF

C cpb 26/1/97 Adding GR output for MPN
      IF(IWRIT4(nr,nx).GE.4) THEN
        WRITE(OP_STRING,'(/'' Global rhs vector GR:'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' NYNR(0,1,1)='',I5)') NYNR(0,1,1)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        CALL OPSTFMAT(NYNR(0,2,1),ISC_GK,ISR_GK,IOOP,NYT(1,1,nx),
     '    NYT(2,1,nx),NZT(1,nx),NYNR(0,1,1),KTYP24,GK,GR,'GK ','GR ',
     '    .FALSE.,.FALSE.,.TRUE.,ERROR,*9999)
      ENDIF

C *** Generate constrained-reduced solution matrix GKK and RHS vector
C$OMP PARALLEL DO DEFAULT(NONE)
C$OMP&  PRIVATE(co1,co2,no_nynr1,no_nynr2,no1,no2,noy1,noy2,ny1,ny2,ny3,
C$OMP&  nz,nzz, PRIVATE_ERROR)
C$OMP&  SHARED(CONY,ERROR_FLAG,NONY,nr,NPNY,NYNE,NYNP,NYNR,nx,
C$OMP&  ISC_GK,ISC_GKK,ISR_GK,ISR_GKK,GK,GKK,GR,GRR,KTYP24,NOT,NYT,
C$OMP&  NZ_GK_M,NZ_GKK_M,NZT,NZZT,SPARSEGKK,UPDATE_MATRIX,YP)
      DO no_nynr1=1,NYNR(0,1,1) !loop over rows of GK
        ny1=NYNR(no_nynr1,1,1) !is row number
        IF(KTYP24.EQ.1) THEN !Compressed Row
          IF(UPDATE_MATRIX) THEN
            DO nz=ISR_GK(ny1),ISR_GK(ny1+1)-1
              ny3=ISC_GK(nz) !local column #
              ny2=GETNYR(1,NPNY,nr,0,2,ny3,NYNE,NYNP) !global column #
              DO noy1=1,NONY(0,ny1,1)
                no1=NONY(noy1,ny1,1) !solution row # attached to ny1
                co1=CONY(noy1,ny1,1) !row coupling coefficient
                DO noy2=1,NONY(0,ny2,2)
                  no2=NONY(noy2,ny2,2) !solution var # attached to ny2
                  co2=CONY(noy2,ny2,2) !row coupling coefficient
                  CALL SPARSE(no1,no2,NOT(1,1,nr,nx),nzz,NZ_GKK_M,
     '              NZZT(1,nr,nx),ISC_GKK,ISR_GKK,SPARSEGKK(nx),
     '              PRIVATE_ERROR,*200)
                  IF(nzz.NE.0) THEN
C$OMP               ATOMIC
                    GKK(nzz)=GKK(nzz)+GK(nz)*co1*co2
                  ENDIF
                ENDDO !noy2
              ENDDO !noy1
            ENDDO !nz
          ENDIF !UPDATE_MATRIX
          ny1=NYNR(no_nynr1,1,1) !is row number
          DO noy1=1,NONY(0,ny1,1)
            no1=NONY(noy1,ny1,1) !solution row # attached to ny1
            co1=CONY(noy1,ny1,1) !row coupling coefficient
            GRR(no1)=GRR(no1)+GR(ny1)*co1
          ENDDO !noy1
        ELSE
          DO noy1=1,NONY(0,ny1,1)
            no1=NONY(noy1,ny1,1) !solution row # attached to ny1
            co1=CONY(noy1,ny1,1) !row coupling coefficient
            IF(UPDATE_MATRIX) THEN
              DO no_nynr2=1,NYNR(0,2,1) !loop over local columns of GK
                ny2=NYNR(no_nynr2,0,1) !global column #
                ny3=NYNR(no_nynr2,2,1) !local column #
                CALL SPARSE(ny1,ny3,NYT(1,1,nx),nz,NZ_GK_M,
     '            NZT(1,nx),ISC_GK,ISR_GK,KTYP24,PRIVATE_ERROR,*200)
                IF(nz.NE.0) THEN
                  DO noy2=1,NONY(0,ny2,2)
                    no2=NONY(noy2,ny2,2) !solution var # attached to ny2
                    co2=CONY(noy2,ny2,2) !row coupling coefficient
                    CALL SPARSE(no1,no2,NOT(1,1,nr,nx),nzz,NZ_GKK_M,
     '                NZZT(1,nr,nx),ISC_GKK,ISR_GKK,SPARSEGKK(nx),
     '                PRIVATE_ERROR,*200)
                    IF(nzz.NE.0) THEN
C$OMP                 ATOMIC
                      GKK(nzz)=GKK(nzz)+GK(nz)*co1*co2
                    ENDIF
C                 reduced LHS matrix
                  ENDDO !noy2
                ENDIF
              ENDDO !no_nynr2
            ENDIF
            GRR(no1)=GRR(no1)+GR(ny1)*co1
          ENDDO !noy1
        ENDIF
        GO TO 205
 200      CONTINUE
          ERROR_FLAG=.TRUE.
          IF(PRIVATE_ERROR.NE.' ') THEN
            CALL FLAG_ERROR(0,
     '        PRIVATE_ERROR(:LEN_TRIM(PRIVATE_ERROR)))
          ENDIF
 205    CONTINUE
      ENDDO !no_nynr1
C$OMP END PARALLEL DO
      IF(ERROR_FLAG) THEN
        ERROR=' '
        GOTO 9999
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)
      IF(IWRIT4(nr,nx).GE.1) THEN
        WRITE(OP_STRING,
     '    '(/'' Global boundary condition reduction time:         '','
     '    //'D11.4,'' s'')') ELAPSED_TIME
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      CALL ASSERT(NZZT(1,nr,nx).GT.0,'>>NZZT(1,nr,nx) is zero',
     '  ERROR,*9999)

      IF(KTYP4.NE.0) THEN !Output global matrices
        CALL WRITE_SOL_MATRIX(ISC_GKK,ISR_GKK,nr,nx,GKK,GRR,ERROR,*9999)
      ENDIF

      IF(IWRIT4(nr,nx).GE.3) THEN
        WRITE(OP_STRING,
     '    '(/'' Global load vector GRR & stiffness matrix GKK:'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' NOT(1,1,nr,nx)='',I5,'
     '    //''', NOT(2,1,nr,nx)='',I5)') NOT(1,1,nr,nx),NOT(2,1,nr,nx)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
c cpb 24/9/95 Adding generic stiffness matrix output
        CALL OPSTFMAT(NYNR(0,1,1),ISC_GKK,ISR_GKK,IOOP,NOT(1,1,nr,nx),
     '    NOT(2,1,nr,nx),NZZT(1,nr,nx),NYNR(0,2,1),SPARSEGKK(nx),
     '    GKK,GRR,'GKK','GRR',.TRUE.,.TRUE.,.TRUE.,ERROR,*9999)
      ENDIF

C JWF 23/11/02  Modify solution system for contact mechanics involving linear elasticity
C     Dependent variable is now geometry and solution is solved via Newton Raphson
C     since problem is now non-linear.

      IF ((KTYP5G(1).GE.1).AND.(ITYP2(1,nx).EQ.1)) THEN
        DO no1=1,NOT(1,1,nr,nx)
          DO no2=1,NOT(2,1,nr,nx)
            CALL SPARSE(no1,no2,NOT(1,1,nr,nx),nzz,NZ_GKK_M,
     '        NZZT(1,nr,nx),ISC_GKK,ISR_GKK,SPARSEGKK(nx),
     '        PRIVATE_ERROR,*9999)
            IF(nzz.NE.0) THEN
              ny2=NYNO(1,no2,2)
              GRR(no1)=GRR(no1)+GKK(nzz)*YP(ny2,1)-GKK(nzz)*YP(ny2,3)
            ENDIF
          ENDDO !no2
        ENDDO !no1
      ENDIF

C*** Solve reduced system

C CPB 23/9/95 New solution system

      CALL SOLVE_SYSTEM(ISC_GKK,ISR_GKK,NOT(1,1,nr,nx),NOT(1,1,nr,nx),
     '  NOT(2,1,nr,nx),NZZT(1,nr,nx),IWRIT4(nr,nx),PRECON_CODE(nx),
     '  SOLVEROPTION(nx),SPARSEGKK(nx),GKK,GRR,XO,FIRSTS,UPDATE_MATRIX,
     '  .TRUE.,nx,ERROR,*9999)

C NEWS 10/08/07 XSL Check the descent condition
C   In order for the residual to decrease for each iteration,
C   it makes sense to have the derivative at the current solution
C   to be negative. Hence the residual will decrease for the 
C   next Newton step. This derivative is calculated by taking the dot
C   product of the gradient vector (GRR) and the direction (XO).

      DIRDERIV=0.0d0
      DO no1=1,NOT(1,1,nr,nx) !loop over global soln rows
          DIRDERIV=DIRDERIV-GRR(no1)*XO(no1)
      ENDDO !no1

C Uncomment to output DIRDERIV
C      WRITE(OP_STRING,
C     &  '(/'' dirDerivative ='',E11.4)')
C     &  DIRDERIV
C      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      IF (DIRDERIV.GT.0) THEN
        WRITE(OP_STRING,'('' >>Warning: Search direciton is not '
     &    //'a descent direction'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)        
      ENDIF
C NEWE XSL

C??? CS8/3/2001 This is how it this was. I can't see why
C??? the original code had been changed. I have put it back
C??? so that multiple ny->no no->ny mappings work
C??? ie so hanging nodes work.
CC *** Store solution increments in YP(ny,5)
C      DO no2=1,NOT(2,1,nr,nx)
C        co1=CYNO(0,no2,2) !additive constant
C        DO nyo2=1,NYNO(0,no2,2)
C          ny2=NYNO(nyo2,no2,2)
C          IF(NPNY(0,ny2,0).EQ.1) THEN
C            np=NPNY(4,ny2,0)
CC GMH 8/1/97 Update cmgui link
C            CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
C          ENDIF
CC          YP(ny2,5)=-co1
CC        ENDDO !nyo2
CC      ENDDO !no2
CC      DO no2=1,NOT(2,1,nr,nx)
CC        DO nyo2=1,NYNO(0,no2,2) !must be a sep loop since adding to YP
CC          ny2=NYNO(nyo2,no2,2)
C          co2=CYNO(nyo2,no2,2)
CC          YP(ny2,5)=YP(ny2,5)-XO(no2)*co2
C          YP(ny2,5)=-co1-XO(no2)*co2
C        ENDDO !nyo2
C      ENDDO !no2

C *** Store solution increments in YP(ny,5)
      DO no_nynr2=1,NYNR(0,2,1) !loop over local columns of GK
        ny2=NYNR(no_nynr2,0,1) !global variable #
C ??? KAT 2003-04-22: How do we choose an no to get the constant addend
C for YP from CYNO(0,no,2)? Which one do we use?  What if there are
C none? Using CONY(0,ny,nrc) might make more sense or a different array
C indexed by ny.  But for now let's just use 0.0 because nothing else
C makes much sense.
        YP(ny2,5)=0.0d0
      ENDDO
C This doesn't initialize YP if there are no no's for the ny.
C       DO no2=1,NOT(2,1,nr,nx)
C         co1=CYNO(0,no2,2) !additive constant
C         DO nyo2=1,NYNO(0,no2,2)
C           ny2=NYNO(nyo2,no2,2)
C           IF(NPNY(0,ny2,0).EQ.1) THEN
C             np=NPNY(4,ny2,0)
C C GMH 8/1/97 Update cmgui link
C             CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
C           ENDIF
C           YP(ny2,5)=-co1
C         ENDDO !nyo2
C       ENDDO !no2
      DO no2=1,NOT(2,1,nr,nx)
        DO nyo2=1,NYNO(0,no2,2) !must be a sep loop since adding to YP
          ny2=NYNO(nyo2,no2,2)
          co2=CYNO(nyo2,no2,2)
          YP(ny2,5)=YP(ny2,5)-XO(no2)*co2
        ENDDO !nyo2
      ENDDO !no2

      CALL EXITS('SOLVE5')
      RETURN
 9999 CALL ERRORS('SOLVE5',ERROR)
      CALL EXITS('SOLVE5')
      RETURN 1
      END


