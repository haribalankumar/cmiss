      SUBROUTINE SOLVE1(ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,ISR_GKK,ISR_GQ,
     '  LGE,NBH,NENP,NHE,NONY,NP_INTERFACE,NPNE,NPNY,nr,NRE,NVHE,nx,
     '  NYNE,NYNO,NYNP,NYNR,CONY,CYNO,GK,GKK,GQ,GR,GRR,XO,YP,
     '  FIRST_A,FIX,UPDATE_MATRIX,UPDATE_VECTOR,ERROR,*)

C#### Subroutine: SOLVE1
C###  Description:
C###    SOLVE1 finds solution of a system of linear equations
C###    no=1,NOT(1,1,nr,nx).  YP(ny,1) contains ics and bcs on entry,
C###    solution values on exit.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'solv00.cmn'
      INCLUDE 'time02.cmn'
!     Parameter List
      INTEGER ISC_GK(NISC_GKM),ISC_GKK(NISC_GKKM),ISC_GQ(NISC_GQM),
     '  ISR_GK(NISR_GKM),ISR_GKK(NISR_GKKM),ISR_GQ(NISR_GQM),
     '  LGE(NHM*NSM,NRCM),NBH(NHM,NCM,NEM),NENP(NPM,0:NEPM,0:NRM),
     '  NHE(NEM),NONY(0:NOYM,NYM,NRCM),NP_INTERFACE(0:NPM,0:3),
     '  NPNE(NNM,NBFM,NEM),
     '  NPNY(0:6,NYM,0:NRCM),nr,NRE(NEM),NVHE(NNM,NBFM,NHM,NEM),nx,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),NYNO(0:NYOM,NOOPM,NRCM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),NYNR(0:NY_R_M,0:NRCM,NCM)
      REAL*8 CONY(0:NOYM,NYM,NRCM),CYNO(0:NYOM,NOOPM,NRCM),GK(NZ_GK_M),
     '  GKK(NZ_GKK_M),GQ(NZ_GQ_M),GR(NYROWM),GRR(NOM),XO(NOM),
     '  YP(NYM,NIYM)
      CHARACTER ERROR*(*)
      LOGICAL FIRST_A,FIX(NYM,NIYFIXM),UPDATE_MATRIX,UPDATE_VECTOR
!     Local Variables
      INTEGER GETNYR,no1,no2,noy1,noy2,no_nynr,no_nynr1,no_nynr2,np,ny3,
     '  nyo1,ny1c,ny1r,ny1v,ny2r,ny2v,nz,nzz
      INTEGER*4 WORK_PTR
      REAL ELAPSED_TIME,TIME_START1(1),TIME_START2(1),TIME_STOP(1)
      REAL*8 co1,co2,co3,SUM

      CALL ENTERS('SOLVE1',*9999)

      IF(NOT(2,1,nr,nx).EQ.0) THEN
        ERROR=' >>The number of unknowns is zero'
        GOTO 9999
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_START1)
      CALL CPU_TIMER(CPU_USER,TIME_START2)

C*** Setup and initialise arrays
C*** Direct assembly of solution matrices

      IF(CALC_GLOBAL(nr,nx)) THEN
        IF(UPDATE_MATRIX) THEN
          IF(FIRST_A) THEN !Temporary work array allocation
            IF(SPARSEGKK(nx).NE.0) THEN
              WORK_PTR=0
              CALL ALLOCATE_MEMORY(NOT(1,1,nr,nx)*NOT(2,1,nr,nx),1,
     '          CHARTYPE,WORK_PTR,MEM_INIT,ERROR,*9999)
            ENDIF
            CALL CALC_SPARSE_GKK(ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,ISR_GKK,
     '        ISR_GQ,LGE,NBH,NENP,NHE,NOT(1,1,nr,nx),NOT(2,1,nr,nx),
     '        NONY,NP_INTERFACE,NPNE,NPNY,nr,NRE,NVHE,nx,NYNE,NYNP,
     '        NYNR,GK,GQ,%VAL(WORK_PTR),.FALSE.,.TRUE.,ERROR,*9999)
            IF(SPARSEGKK(nx).NE.0) CALL FREE_MEMORY(WORK_PTR,
     '        ERROR,*9999)
          ENDIF !FIRST_A
C$OMP PARALLEL DO
C$OMP&    PRIVATE(nzz),
C$OMP&    SHARED(NZZT,GKK)
          DO nzz=1,NZZT(1,nr,nx)
            GKK(nzz)=0.0d0
          ENDDO !nzz
C$OMP END PARALLEL DO
        ENDIF !UPDATE_MATRIX
      ENDIF !CALC_GLOBAL
      IF(UPDATE_VECTOR) THEN
C$OMP PARALLEL DO
C$OMP&  PRIVATE(no1),
C$OMP&  SHARED(NOT,GRR)
        DO no1=1,NOT(1,1,nr,nx)
          GRR(no1)=0.0d0
        ENDDO !no1
C$OMP END PARALLEL DO
      ENDIF

C*** Update GRR from essential & flux bcs and GK from mixed bcs
C Essential bcs
      DO no_nynr=1,NYNR(0,0,1) !loop over global variables for nc=1
        ny1v=NYNR(no_nynr,0,1) !global var#
        ny1c=GETNYR(1,NPNY,nr,2,0,ny1v,NYNE,NYNP) !local column#
        ny2v=GETNYR(2,NPNY,nr,0,0,ny1v,NYNE,NYNP) !equiv global flux var#
        IF(FIX(ny1v,1)) THEN !global var set as bc
          DO no_nynr1=1,NYNR(0,1,1) !loop over the rows of nc=1
            ny1r=NYNR(no_nynr1,1,1) !row#
            DO noy1=1,NONY(0,ny1r,1)
              no1=NONY(noy1,ny1r,1)
              co1=CONY(noy1,ny1r,1)
              CALL SPARSE(ny1r,ny1c,NYT(1,1,nx),nz,NZ_GK_M,
     '          NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
              IF(nz.NE.0) GRR(no1)=GRR(no1)-GK(nz)*YP(ny1v,1)*co1
            ENDDO !noy1
          ENDDO !no_nynr1
        ENDIF !fix
      ENDDO !no_nynr
C Flux bcs
      DO no_nynr=1,NYNR(0,0,2) !loop over global variables for nc=2
        ny1v=NYNR(no_nynr,0,1) !global var#
        ny2v=NYNR(no_nynr,0,2) !global flux#
        IF(FIX(ny2v,1).AND..NOT.FIX(ny2v,2)) THEN   !flux is set as a bc
          ny2r=GETNYR(2,NPNY,nr,1,0,ny2v,NYNE,NYNP) !global flux row#
          DO noy1=1,NONY(0,ny2r,1)
            no1= NONY(noy1,ny2r,1)
            co1= CONY(noy1,ny2r,1)
            GRR(no1)=GRR(no1)+YP(ny2v,1)*co1 !GRR from applied flux bc
          ENDDO !noy1
        ENDIF !FIX
      ENDDO !no_nynr
C Mixed bcs
      DO no_nynr=1,NYNR(0,0,1) !loop over global variables for nc=1
        ny1v=NYNR(no_nynr,0,1) !global var#
        ny2v=NYNR(no_nynr,0,2) !global flux#
        ny1c=GETNYR(1,NPNY,nr,2,0,ny1v,NYNE,NYNP) !local column#
        IF(FIX(ny2v,1).AND.FIX(ny2v,2)) THEN !flux=au
          DO no_nynr1=1,NYNR(0,1,1) !loop over the rows of nc=1
            ny1r=NYNR(no_nynr1,1,1) !row#
            IF(ny1r.EQ.ny1c) THEN
              CALL SPARSE(ny1r,ny1c,NYT(1,1,nx),nz,NZ_GK_M,
     '          NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
              IF(nz.NE.0) GK(nz)=GK(nz)-YP(ny1v,2) !!!!needs coupling coeff?
              IF(DOP.AND.nz.NE.0) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP           CRITICAL(SOLVE1_1)
                WRITE(OP_STRING,'(/'' GK updated by mixed bc at nz='','
     '            //'I8,'' ny2v='',I8,'' to GK='',D11.4)')
     '            nz,ny2v,GK(nz)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP           END CRITICAL(SOLVE1_1)
              ENDIF !DOP
            ENDIF !ny1r=ny1c
          ENDDO !no_nynr1
        ENDIF !fix
      ENDDO !no_nynr

C*** Generate reduced system

      DO no_nynr1=1,NYNR(0,1,1) !loop over global rows of GK
        ny1r=NYNR(no_nynr1,1,1)  !row#
        DO noy1=1,NONY(0,ny1r,1) !loop over #no's attached to row ny1r
          no1=NONY(noy1,ny1r,1)  !no# attached to row ny1r
          co1=CONY(noy1,ny1r,1)  !coupling coeff for row mapping
C                                ie row_no1=a*row_ny1r+b*row_ny..
          GRR(no1)=GRR(no1)+GR(ny1r)*co1 !get reduced RHS vector
          IF(CALC_GLOBAL(nr,nx)) THEN
            DO no_nynr2=1,NYNR(0,0,1) !loop over the #cols of GK
              ny1v=NYNR(no_nynr2,0,1) !global var#
              ny1c=GETNYR(1,NPNY,nr,2,0,ny1v,NYNE,NYNP) !local col#
              CALL SPARSE(ny1r,ny1c,NYT(1,1,nx),nz,NZ_GK_M,
     '          NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
              IF(nz.NE.0) THEN
                IF(UPDATE_MATRIX) THEN
                  DO noy2=1,NONY(0,ny1v,2) !loop over #no's for var ny1v
                    no2=NONY(noy2,ny1v,2)  !no# attached to ny1v
                    co2=CONY(noy2,ny1v,2)  !coup coeff for the col mapping
C                                          i.e. var_no1=a*var_ny1r+b*var_ny..
                    CALL SPARSE(no1,no2,NOT(1,1,nr,nx),nzz,NZ_GKK_M,
     '                NZZT(1,nr,nx),ISC_GKK,ISR_GKK,SPARSEGKK(nx),
     '                ERROR,*9999)
                    IF(nzz.NE.0) GKK(nzz)=GKK(nzz)+GK(nz)*co1*co2
                  ENDDO !noy2
                ENDIF
C!!! There is only one constant defined for each variable therefore
C!!! it is implied that this constant is applied to the last ny mapped
C!!! to the no. ie. no=a*ny1r+b*ny1v+(c*ny1c+d) that is a*ny1r=b*ny1v and
C!!! a*ny1r=c*ny1c+d and not a*ny1r=b*ny1v+d etc.
                co3=CONY(0,ny1v,2) !is additive constant applied to vars
                GRR(no1)=GRR(no1)-GK(nz)*co3 !additive const. in RHS vec
              ENDIF !nz.NE.0
            ENDDO !no_nynr2
          ENDIF !CALC_GLOBAL
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
      IF(UPDATE_MATRIX) THEN
        IF(KTYP4.NE.0) THEN !Output global matrices
          CALL WRITE_SOL_MATRIX(ISC_GKK,ISR_GKK,nr,nx,GKK,GRR,
     '      ERROR,*9999)
        ENDIF

        IF(IWRIT4(nr,nx).GE.3) THEN
          WRITE(OP_STRING,
     '      '(/'' Global load vector GRR & stiffness matrix GKK:'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NOT(1,1,nr,nx)='',I5,'
     '      //''', NOT(2,1,nr,nx)='',I5)') NOT(1,1,nr,nx),NOT(2,1,nr,nx)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
c cpb 24/9/95 Adding generic stiffness matrix output
          CALL OPSTFMAT(NYNR(0,1,1),ISC_GKK,ISR_GKK,IOOP,NOT(1,1,nr,nx),
     '      NOT(2,1,nr,nx),NZZT(1,nr,nx),NYNR(0,2,1),SPARSEGKK(nx),GKK,
     '      GRR,'GKK','GRR',.TRUE.,.TRUE.,.TRUE.,ERROR,*9999)
        ENDIF
      ENDIF !UPDATE_MATRIX

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
      IF(UPDATE_MATRIX) THEN
        IF(IWRIT4(nr,nx).GE.1) THEN
          IF(IWRIT4(nr,nx).GE.4.OR.KTYP4.NE.0) THEN
            WRITE(OP_STRING,'(/'' CPU time for solution matrix '
     '        //'output: '',D11.4,'' s'')') ELAPSED_TIME
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
      ENDIF

C*** Solve reduced system
      CALL CPU_TIMER(CPU_USER,TIME_START2)

C CPB 23/9/95 New solution system

      CALL SOLVE_SYSTEM(ISC_GKK,ISR_GKK,NOT(1,1,nr,nx),NOT(1,1,nr,nx),
     '  NOT(2,1,nr,nx),NZZT(1,nr,nx),IWRIT4(nr,nx),PRECON_CODE(nx),
     '  SOLVEROPTION(nx),SPARSEGKK(nx),GKK,GRR,XO,FIRST_A,UPDATE_MATRIX,
     '  .TRUE.,nx,ERROR,*9999)

C*** Put solution values into YP(ny1c,1)

      DO no1=1,NOT(2,1,nr,nx)     !initialise YP for solution dofs
        DO nyo1=1,NYNO(0,no1,2)
          ny1c=NYNO(nyo1,no1,2)
          IF(NPNY(0,ny1c,0).EQ.1) THEN
            np=NPNY(4,ny1c,0)
            CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
          ENDIF
          IF(.NOT.FIX(ny1c,1)) THEN
            YP(ny1c,1)=0.0d0
          ENDIF
        ENDDO !nyo1
      ENDDO !no1
      DO no1=1,NOT(2,1,nr,nx)
        DO nyo1=1,NYNO(0,no1,2)
          ny1c=NYNO(nyo1,no1,2)
          co1 =CYNO(nyo1,no1,2)
          YP(ny1c,1)=YP(ny1c,1)+XO(no1)*co1
        ENDDO !nyo1
        IF(NYNO(0,no1,2).GT.0) THEN !additive constant
          co1=CYNO(0,no1,2)
          YP(ny1c,1)=YP(ny1c,1)+co1
        ENDIF
      ENDDO !no1

C*** Backsubstitute to find unknown fluxes

      DO no_nynr=1,NYNR(0,0,1) !loop over the global variables of nc=1
        ny1v=NYNR(no_nynr,0,1) !is global var#
        ny2v=NYNR(no_nynr,0,2) !is global flux#
        IF(FIX(ny1v,1)) THEN !global var set as bc
          ny1r=NYNR(no_nynr,1,1) !row#
          SUM=0.0d0
          DO no_nynr2=1,NYNR(0,2,1) !loop over local columns
            ny1c=NYNR(no_nynr2,2,1) !is local column#
            ny3=GETNYR(1,NPNY,nr,0,2,ny1c,NYNE,NYNP)
            !is global variable number for local column ny1c
            CALL SPARSE(ny1r,ny1c,NYT(1,1,nx),nz,NZ_GK_M,
     '        NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
            IF(nz.NE.0) SUM=SUM+GK(nz)*YP(ny3,1)
          ENDDO
          ny2r=GETNYR(2,NPNY,nr,0,0,ny1v,NYNE,NYNP) !correspond. RHS var#
          IF(NPNY(0,ny2r,0).EQ.1) THEN
            np=NPNY(4,ny2r,0)
            CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
          ENDIF
          YP(ny2r,1)=SUM
        ELSE IF(FIX(ny2v,1).AND.FIX(ny2v,2)) THEN !flux=au bc
          YP(ny2v,1)=YP(ny1v,2)*YP(ny1v,1) !=coeff x solution var
        ENDIF
      ENDDO !no_nynr

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
      IF(IWRIT4(nr,nx).GE.1) THEN
        WRITE(OP_STRING,'(/'' CPU time for storage and back '
     '    //'substitution: '',D11.4,'' s'')') ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START1(1)
      IF(IWRIT4(nr,nx).GE.1) THEN
        WRITE(OP_STRING,'(/'' Total CPU time for solution: '','
     '    //'D11.4,'' s'')') ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('SOLVE1')
      RETURN
 9999 CALL ERRORS('SOLVE1',ERROR)
      CALL EXITS('SOLVE1')
      RETURN 1
      END


