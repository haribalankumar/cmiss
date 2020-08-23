      SUBROUTINE SOLVE2(ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,ISR_GKK,ISR_GQ,
     '  LGE,NBH,NENP,NHE,NDIPOLES,NONY,NP_INTERFACE,NPNE,NPNY,nr,
     '  NRE,NVHE,nx,NYNE,NYNO,NYNP,NYNR,CONY,CYNO,GD,GK,GKK,GQ,GRR,
     '  XO,YP,FIRST,FIX,UPDATE_MATRIX,UPDATE_SOURCE,UPDATE_VECTOR,
     '  ERROR,*)

C#### Subroutine: SOLVE2
C###  Description:
C###    SOLVE2 finds solution of an unsymmetric system of BE
C###    linear equations no=1,NOT(1,1,nr,nx).

C**** Equations are assembled as GK*U + GD = GQ*dUdN
C**** and reordered as GKK*XO=GRR.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'solv00.cmn'
      INCLUDE 'time02.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER ISC_GK(NISC_GKM),ISC_GKK(NISC_GKKM),ISC_GQ(NISC_GQM),
     '  ISR_GK(NISR_GKM),ISR_GKK(NISR_GKKM),ISR_GQ(NISR_GQM),
     '  LGE(NHM*NSM,NRCM),NBH(NHM,NCM,NEM),NENP(NPM,0:NEPM,0:NRM),
     '  NHE(NEM),NDIPOLES(NRM,NXM),NONY(0:NOYM,NYM,NRCM),
     '  NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),NPNY(0:6,NYM,0:NRCM),
     '  nr,NRE(NEM),NVHE(NNM,NBFM,NHM,NEM),nx,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),NYNO(0:NYOM,NOOPM,NRCM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),NYNR(0:NY_R_M,0:NRCM,NCM)
      REAL*8 CONY(0:NOYM,NYM,NRCM),CYNO(0:NYOM,NOOPM,NRCM),GD(NZ_GD_M),
     '  GK(NZ_GK_M),GKK(NZ_GKK_M),GQ(NZ_GQ_M),GRR(NOM),XO(NOM),
     '  YP(NYM,NIYM)
      CHARACTER ERROR*(*)
      LOGICAL FIRST,FIX(NYM,NIYFIXM),UPDATE_MATRIX,
     '  UPDATE_SOURCE,UPDATE_VECTOR
!     Local Variables
      INTEGER GETNYR,nc,nk,no1,no2,no_nynr,no_nynr1,
     '  no_nynr2,noy1,noy2,np,ny,ny1,ny2,ny3,nyo2,nz,nz1,nzz
      INTEGER*4 WORK_PTR
      REAL ELAPSED_TIME,TIME_START1(1),TIME_START2(1),TIME_STOP(1)
      REAL*8 ALPHA, co1,co2,SUM1,SUM2,P1
      LOGICAL ERROR_FLAG,UPDATE_MATRIX_DUMMY,X_INIT

      CALL ENTERS('SOLVE2',*9999)

      IF(NOT(2,1,nr,nx).EQ.0) THEN
        ERROR=' >>The number of unknowns is zero'
        GOTO 9999
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_START1)
      CALL CPU_TIMER(CPU_USER,TIME_START2)

C*** Initialise arrays

C CPB 29/10/98 Adding direct assembly of solution matrices

      IF(CALC_GLOBAL(nr,nx)) THEN
c cpb 23/9/95 New solution system
        IF(UPDATE_MATRIX) THEN
          IF(FIRST) THEN
C CPB 6/11/95 Temporary work array allocation

            IF(SPARSEGKK(nx).NE.0) THEN
              WORK_PTR=0
              CALL ALLOCATE_MEMORY(NOT(1,1,nr,nx)*NOT(2,1,nr,nx),1,
     '          CHARTYPE,WORK_PTR,MEM_INIT,ERROR,*9999)
            ENDIF
            CALL CALC_SPARSE_GKK(ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,ISR_GKK,
     '        ISR_GQ,LGE,NBH,NENP,NHE,NOT(1,1,nr,nx),NOT(2,1,nr,nx),
     '        NONY,NP_INTERFACE,NPNE,NPNY,nr,NRE,NVHE,nx,NYNE,NYNP,
     '        NYNR,GK,GQ,%VAL(WORK_PTR),.FALSE.,.TRUE.,ERROR,*9999)
            IF(SPARSEGKK(nx).NE.0) CALL FREE_MEMORY(WORK_PTR,ERROR,
     '        *9999)
          ENDIF
C$OMP PARALLEL DO
C$OMP&  PRIVATE(nzz),
C$OMP&  SHARED(NZZT,GKK)
          DO nzz=1,NZZT(1,nr,nx)
            GKK(nzz)=0.0d0
          ENDDO !nzz
C$OMP END PARALLEL DO
        ENDIF
      ENDIF
      IF(UPDATE_VECTOR) THEN
C$OMP PARALLEL DO
C$OMP&  PRIVATE(no1),
C$OMP&  SHARED(NOT,GRR)
        DO no1=1,NOT(1,1,nr,nx)
          GRR(no1)=0.0d0
        ENDDO !no1
C$OMP END PARALLEL DO
      ENDIF

C*** Generate reordered system

      IF(UPDATE_MATRIX.OR.UPDATE_VECTOR) THEN
        IF(CALC_GLOBAL(nr,nx)) THEN
          ERROR_FLAG=.FALSE.
C$OMP PARALLEL DO
C$OMP&  PRIVATE(no_nynr1,ny1,noy1,no1,co1,no_nynr2,ny2,ny3,
C$OMP&          nz1,noy2,no2,co2,nzz,nc),
C$OMP&  SHARED(KTYP24,nr,NZ_GK_M,NZ_GKK_M,NZ_GQ_M,ERROR_FLAG)
          DO no_nynr1=1,NYNR(0,1,1) !Loop over the # of rows of GK
            IF(.NOT.ERROR_FLAG) THEN
              ny1=NYNR(no_nynr1,1,1) !is the row #
              DO noy1=1,NONY(0,ny1,1) !loop over the # of no's attached
C                                      to the row ny1
                no1=NONY(noy1,ny1,1) !is the row # attached to ny1
                co1=CONY(noy1,ny1,1) !is the coupling coeff for the row
C                                     mapping ie row_no1=a*row_ny1+
C                                     b*row_ny2
                nc=1 !GK
                DO no_nynr2=1,NYNR(0,0,nc) !loop over # of cols of GK
                  ny2=NYNR(no_nynr2,0,nc) !is the global GK variable #
                  ny3=GETNYR(nc,NPNY,nr,2,0,ny2,NYNE,NYNP)
C                 is local GK variable #
                  CALL SPARSE(ny1,ny3,NYT(1,nc,nx),nz1,NZ_GK_M,
     '              NZT(nc,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*100)
C                 Check if fixed, else there exists one or more no
                  IF(FIX(ny2,1)) THEN
C                   Move over to the RHS
                    IF(UPDATE_VECTOR) THEN
                      GRR(no1)=GRR(no1)-GK(nz1)*YP(ny2,1)*co1 !*co2??
                    ENDIF
                  ELSE
C                   Loop over any solution variables attached to the GK
C                   variable
                    DO noy2=1,NONY(0,ny2,2) !loop over the # of no's
C                                            attached to the solution
C                                            variable ny2
                      no2=NONY(noy2,ny2,2) !is column # attached to ny2
                      co2=CONY(noy2,ny2,2) !is the coup coeff for the
C                                           column mapping i.e. var_no1=
C                                           a*var_ny1+b*var_ny2
                      IF(UPDATE_MATRIX) THEN
                        CALL SPARSE(no1,no2,NOT(1,1,nr,nx),nzz,
     '                    NZ_GKK_M,NZZT(1,nr,nx),ISC_GKK,ISR_GKK,
     '                    SPARSEGKK(nx),ERROR,*100)
                        IF(nzz.NE.0) THEN
                          GKK(nzz)=GKK(nzz)+GK(nz1)*co1*co2
                        ENDIF
                      ENDIF
                    ENDDO !noy2
                  ENDIF !FIX
                ENDDO !no_nynr2
                nc=2 !GQ
                DO no_nynr2=1,NYNR(0,0,nc) !loop over # of cols of GQ
                  ny2=NYNR(no_nynr2,0,nc) !is the global GQ variable #
                  ny3=GETNYR(nc,NPNY,nr,2,0,ny2,NYNE,NYNP)
C                 is local GQ variable #
                  CALL SPARSE(ny1,ny3,NYT(1,nc,nx),nz1,NZ_GQ_M,
     '              NZT(nc,nx),ISC_GQ,ISR_GQ,KTYP24,ERROR,*100)
C                 Check if fixed, else there exists one or more no
                  IF(FIX(ny2,1)) THEN
C                   Move over to the RHS
                    IF(UPDATE_VECTOR) THEN
                      GRR(no1)=GRR(no1)+GQ(nz1)*YP(ny2,1)*co1 !*co2??
                    ENDIF
                  ELSE
C                   Loop over any solution variables attached to the GK
C                   variable
                    DO noy2=1,NONY(0,ny2,2) !loop over the # of no's
C                                            attached to the solution
C                                          variable ny2
                      no2=NONY(noy2,ny2,2) !is column # attached to ny2
                      co2=CONY(noy2,ny2,2) !is the coup coeff for the
C                                           column mapping i.e. var_no1=
C                                           a*var_ny1+b*var_ny2
                      IF(UPDATE_MATRIX) THEN
                        CALL SPARSE(no1,no2,NOT(1,1,nr,nx),nzz,
     '                    NZ_GKK_M,NZZT(1,nr,nx),ISC_GKK,ISR_GKK,
     '                    SPARSEGKK(nx),ERROR,*100)
                        IF(nzz.NE.0) THEN
                          GKK(nzz)=GKK(nzz)-GQ(nz1)*co1*co2
                        ENDIF
                      ENDIF
                    ENDDO !noy2
                  ENDIF !FIX
                ENDDO !no_nynr2
C!!! Additive constants not implemented for BEM
              ENDDO !noy1
              GOTO 102
C             This statement is designed to be skipped if no error
C             occur. However if a error occurs within a subroutine
C             the alternate return points to line 100 to set the flag
 100          CONTINUE
C$OMP CRITICAL(SOLVE2_1)
              ERROR_FLAG=.TRUE.
              WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
              CALL WRITES(IOER,OP_STRING,ERROR,*101)
              WRITE(OP_STRING,'(/'' >>An error occurred - '
     '          //'results may be unreliable!'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*101)
 101          CONTINUE
C$OMP END CRITICAL(SOLVE2_1)
 102          CONTINUE
            ENDIF
          ENDDO !no_nynr1
C$OMP END PARALLEL DO
          CALL ASSERT(.NOT.ERROR_FLAG,'>>An error occurred during '
     '      //'global stiffness matrix assembly',ERROR,*9999)
        ELSE IF(UPDATE_VECTOR) THEN
C$OMP PARALLEL
C$OMP&  PRIVATE(no_nynr,ny1,ny2,no_nynr1,ny,noy1,no1,co1,nz),
C$OMP&  SHARED(KTYP24,nr,NZ_GK_M,NZ_GQ_M,GRR,ERROR_FLAG)
          ERROR_FLAG=.FALSE.
C$OMP DO
          DO no_nynr=1,NYNR(0,0,1) !Loop over the glob variables of nc=1
            IF(.NOT.ERROR_FLAG) THEN
              ny1=NYNR(no_nynr,0,1) !is the global variable#
              ny2=GETNYR(1,NPNY,nr,2,0,ny1,NYNE,NYNP) !is the local col#
              IF(FIX(ny1,1)) THEN !global variable is set as bc
                DO no_nynr1=1,NYNR(0,1,1) !Loop over the rows of nc=1
                  ny=NYNR(no_nynr1,1,1) !is the row#
                  DO noy1=1,NONY(0,ny,1)
                    no1=NONY(noy1,ny,1)
                    co1=CONY(noy1,ny,1)
                    CALL SPARSE(ny,ny2,NYT(1,1,nx),nz,NZ_GK_M,
     '                NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*200)
                    IF(nz.NE.0) GRR(no1)=GRR(no1)-GK(nz)*YP(ny1,1)*co1
                  ENDDO !noy1
                ENDDO !no_nynr1
              ENDIF !fix
              GOTO 202
C             This statement is designed to be skipped if no error
C             occur. However if a error occurs within a subroutine
C             the alternate return points to line 100 to set the flag
 200          CONTINUE
C$OMP CRITICAL(SOLVE2_2)
              ERROR_FLAG=.TRUE.
              WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
              CALL WRITES(IOER,OP_STRING,ERROR,*201)
              WRITE(OP_STRING,'(/'' >>An error occurred - '
     '          //'results may be unreliable!'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*201)
 201          CONTINUE
C$OMP END CRITICAL(SOLVE2_2)
 202          CONTINUE
            ENDIF
          ENDDO !no_nynr
C$OMP END DO NOWAIT
C$OMP DO
          DO no_nynr=1,NYNR(0,0,2) !Loop over the global vars of nc=2
            IF(.NOT.ERROR_FLAG) THEN
              ny1=NYNR(no_nynr,0,2) !is the global variable#
              ny2=GETNYR(2,NPNY,nr,2,0,ny1,NYNE,NYNP) !is the loc col #
              IF(FIX(ny1,1)) THEN !global variable is set as bc.
                DO no_nynr1=1,NYNR(0,1,2) !Loop over the rows of nc=2
                  ny=NYNR(no_nynr1,1,2) !is the row#
                  DO noy1=1,NONY(0,ny,1)
                    no1=NONY(noy1,ny,1)
                    co1=CONY(noy1,ny,1)
                    CALL SPARSE(ny,ny2,NYT(1,2,nx),nz,NZ_GQ_M,
     '                NZT(2,nx),ISC_GQ,ISR_GQ,KTYP24,ERROR,*300)
                    IF(nz.NE.0) GRR(no1)=GRR(no1)+GQ(nz)*YP(ny1,1)*co1
                  ENDDO !noy1
                ENDDO !no_nynr1
              ENDIF !fix
              GOTO 302
C             This statement is designed to be skipped if no error
C             occur. However if a error occurs within a subroutine
C             the alternate return points to line 100 to set the flag
 300          CONTINUE
C$OMP CRITICAL(SOLVE2_3)
              ERROR_FLAG=.TRUE.
              WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
              CALL WRITES(IOER,OP_STRING,ERROR,*301)
              WRITE(OP_STRING,'(/'' >>An error occurred - '
     '          //'results may be unreliable!'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*301)
 301          CONTINUE
C$OMP END CRITICAL(SOLVE2_3)
 302          CONTINUE
            ENDIF
          ENDDO !no_nynr
C$OMP END DO
C$OMP END PARALLEL
          CALL ASSERT(.NOT.ERROR_FLAG,'>>An error occurred during '
     '      //'global stiffness matrix assembly',ERROR,*9999)
        ENDIF
      ENDIF


      IF(UPDATE_SOURCE.AND.USE_DIPOLE.EQ.1) THEN

C GMH 6/11/95 IF(FIRST.AND.ITYP4(nr,nx).EQ.4) THEN
C GMH 16/11/95 should be ITYP5
        IF(ITYP5(nr,nx).EQ.4) THEN
C*** Quasistatic solution and first time. Save GRR in the YP(..,5)
C*** storage location so that you can increment GRR if necessary later
C*** on.
          IF(FIRST) THEN
            CALL DCOPY(NOT(1,1,0,nx),GRR,1,YP(1,5),1)
          ENDIF
          IF(IGREN(nr).GE.9.OR.NDIPOLES(nr,nx).GT.0) THEN

C LKC 25-JAN-2006 ERROR_FLAG did not appear to be initialised - results
C in the IF section sometimes not being executed! Also added ASSERT at the end to flag
C if the ERROR_FLAG was activated.
            ERROR_FLAG=.FALSE.

C cpb 28/11/97 Adding multiprocessing. This will not work in the
C case of mulitple rows mapped to one row
C$OMP PARALLEL DO
C$OMP&  PRIVATE(no_nynr1,np,ny1,noy1,no1,co1),
C$OMP&  SHARED(NONY,NYNR,CONY,GRR,YP,GD,ERROR_FLAG)

            DO no_nynr1=1,NYNR(0,1,1) !Loop over the # of rows of GK
              ny1=NYNR(no_nynr1,1,1) !is the row #
              IF(.NOT.ERROR_FLAG) THEN
                IF(NPNY(0,ny1,0).EQ.1) THEN
                  np=NPNY(4,ny1,0)
C GMH 8/1/97 Update cmgui link
                  CALL NODE_CHANGE(np,.FALSE.,ERROR,*400)
                ENDIF
                
                DO noy1=1,NONY(0,ny1,1) !loop over the # of no attached
C                                        to the row ny1
                  no1=NONY(noy1,ny1,1) !is the row # attached to ny1
                  co1=CONY(noy1,ny1,1) !is the coupling coeff for row
C                                 mapping ie row_no1=a*row_ny1+b*row_ny2
C CPB 4/6/95 Only change GRR by increment in source vector (assumes
C rhs vector without any source terms is temporarily stored in YP(ny,5))
                  GRR(no1)=YP(no1,5)-GD(ny1)*co1
                ENDDO !noy1
                GO TO 402
C               This statement is designed to be skipped if no error
C               occur. However if a error occurs within a subroutine
C               the alternate return points to line 200 to set the flag
 400            CONTINUE
C$OMP CRITICAL(SOLVE2_4)
                ERROR_FLAG=.TRUE.
                WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
                CALL WRITES(IOER,OP_STRING,ERROR,*401)
                WRITE(OP_STRING,'(/'' >>An error occurred - '
     '            //'results may be unreliable!'')')
                CALL WRITES(IOER,OP_STRING,ERROR,*401)
 401            CONTINUE
C$OMP END CRITICAL(SOLVE2_4)
 402            CONTINUE
              ENDIF !ERROR_FLAG
            ENDDO !no_nynr1 (ny1)
C$OMP END PARALLEL DO
            CALL ASSERT(.NOT.ERROR_FLAG,'>>An error occurred  '
     '        //'in multiprocessing updating the source',ERROR,*9999)
          ENDIF ! IGREN

C*** Static solutions etc.
        ELSE
C GMH 6/11/95 Copied from above, GRR substituted for YP(no1,5)
          IF(IGREN(nr).GE.9.OR.NDIPOLES(nr,nx).GT.0) THEN
            DO no_nynr1=1,NYNR(0,1,1) !Loop over the # of rows of GK
              ny1=NYNR(no_nynr1,1,1) !is the row #
              DO noy1=1,NONY(0,ny1,1) !loop over the # of no's attached
C                                    to the row ny1
                no1=NONY(noy1,ny1,1) !is the row # attached to ny1
                co1=CONY(noy1,ny1,1) !is the coupling coeff for the row
C                                 mapping ie row_no1=a*row_ny1+b*row_ny2
C GM 6/11/95 Add source vector to GRR
                GRR(no1)=GRR(no1)-GD(ny1)*co1
              ENDDO !noy1
            ENDDO !no_nynr1 (ny1)
          ENDIF
        ENDIF
      ENDIF !Update Source

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
      IF(IWRIT4(nr,nx).GE.1) THEN
        WRITE(OP_STRING,'(/'' CPU time for solution matrix '
     '    //'initialisation and assembly: '',D11.4,'' s'')')
     '    ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF
      CALL CPU_TIMER(CPU_USER,TIME_START2)

      IF(UPDATE_MATRIX.OR.UPDATE_VECTOR) THEN

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
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)

      IF(UPDATE_MATRIX.OR.UPDATE_VECTOR) THEN
        IF(IWRIT4(nr,nx).GE.1) THEN
          IF(IWRIT4(nr,nx).GE.4.OR.KTYP4.NE.0) THEN
            WRITE(OP_STRING,'(/'' CPU time for solution matrix '
     '        //'output: '',D11.4,'' s'')') ELAPSED_TIME
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
      ENDIF

      ! Use zero as an initial solution guess for the first iteration.
      ! Otherwise use the solution from the previous time.
      ! A separate variable to FIRST is required because FIRST is
      ! modified in SOLVE_SYSTEM.
      X_INIT=FIRST

      UPDATE_MATRIX_DUMMY=UPDATE_MATRIX !As we need two solves without
C                                        refactorisation
      IF(SALU_CONSISTENCY(nx).AND.FIRST) THEN
C       Initialise the rhs with 1's to find phi_j^1
C$OMP PARALLEL DO
C$OMP&  PRIVATE(no1),
C$OMP&  SHARED(YP)
        DO no1=1,NOT(1,1,nr,nx)
          YP(no1,6)=1.0d0
        ENDDO
C$OMP END PARALLEL DO
C       Solve to find phi_j^1
        CALL SOLVE_SYSTEM(ISC_GKK,ISR_GKK,NOT(1,1,nr,nx),NOT(1,1,nr,nx),
     '    NOT(2,1,nr,nx),NZZT(1,nr,nx),IWRIT4(nr,nx),PRECON_CODE(nx),
     '    SOLVEROPTION(nx),SPARSEGKK(nx),GKK,YP(1,6),XO,FIRST,
     '    UPDATE_MATRIX_DUMMY,X_INIT,nx,ERROR,*9999)
        UPDATE_MATRIX_DUMMY=.FALSE.
        CALL DCOPY(NOT(1,1,nr,nx),XO,1,YP(1,6),1)
      ENDIF

C*** Solve reordered system

C CPB 23/9/95 New solution system

      CALL SOLVE_SYSTEM(ISC_GKK,ISR_GKK,NOT(1,1,nr,nx),NOT(1,1,nr,nx),
     '  NOT(2,1,nr,nx),NZZT(1,nr,nx),IWRIT4(nr,nx),PRECON_CODE(nx),
     '  SOLVEROPTION(nx),SPARSEGKK(nx),GKK,GRR,XO,FIRST,
     '  UPDATE_MATRIX_DUMMY,X_INIT,nx,ERROR,*9999)

      IF(SALU_CONSISTENCY(nx)) THEN
C       Find Alpha
        SUM1=0.0d0
        SUM2=0.0d0
        ny1=GETNYR(1,NPNY,nr,1,0,SALU_NY(nx),NYNE,NYNP)
        P1=0.0d0
        DO no_nynr2=1,NYNR(0,0,1)
          ny2=NYNR(no_nynr2,0,1)
          ny3=GETNYR(1,NPNY,nr,2,0,ny2,NYNE,NYNP)
          CALL SPARSE(ny1,ny3,NYT(1,1,nx),nz1,NZ_GK_M,
     '      NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
          IF(FIX(ny2,1)) THEN
            P1=P1-GK(nz1)*YP(ny2,1)
          ELSE
            DO noy2=1,NONY(0,ny2,2)
              no2=NONY(noy2,ny2,2)
              co2=CONY(noy2,ny2,2)
              SUM1=SUM1+GK(nz1)*XO(no2)*co2
              SUM2=SUM2+GK(nz1)*YP(no2,6)*co2
            ENDDO !noy2
          ENDIF !FIX
        ENDDO !no_nynr2
        DO no_nynr2=1,NYNR(0,0,2)
          ny2=NYNR(no_nynr2,0,2)
          ny3=GETNYR(2,NPNY,nr,2,0,ny2,NYNE,NYNP)
          CALL SPARSE(ny1,ny3,NYT(1,2,nx),nz1,NZ_GQ_M,
     '      NZT(2,nx),ISC_GQ,ISR_GQ,KTYP24,ERROR,*9999)
          IF(FIX(ny2,1)) THEN
            P1=P1+GQ(nz1)*YP(ny2,1)
          ELSE
            DO noy2=1,NONY(0,ny2,2)
              no2=NONY(noy2,ny2,2)
              co2=CONY(noy2,ny2,2)
              SUM1=SUM1-GQ(nz1)*XO(no2)*co2
              SUM2=SUM2-GQ(nz1)*YP(no2,6)*co2
            ENDDO !noy2
          ENDIF !FIX
        ENDDO !no_nynr2

C LKC 17-AUG-1999 Incorrect if we are not using dipoles
C
C        IF(IGREN(nr).GE.9.OR.NDIPOLES(nr,nx).GT.0) THEN
C          P1=P1-GD(ny1)
C        ENDIF
        IF(USE_DIPOLE.EQ.1) THEN
          IF(IGREN(nr).GE.9.OR.NDIPOLES(nr,nx).GT.0) THEN
            P1=P1-GD(ny1)
          ENDIF
        ELSE
          IF(IGREN(nr).GE.9) THEN
            P1=P1-GD(ny1)
          ENDIF
        ENDIF ! USE_DIPOLES=1


C cpb 20/8/97 Modifing the salu consistency criteria to allow for
C the fixing of a potential at an arbitrary value. i.e.
C        ALPHA=(SUM1-P1+AC)/(1.0d0-SUM2)
C Note: that AC is not acutally added as it is already subtracted
C from P1 above
        IF((DABS(SUM1-P1).LT.ZERO_TOL).OR.
     '    (DABS(1.0d0-SUM2).LT.ZERO_TOL)) THEN
          ALPHA=0.0d0
        ELSE
          ALPHA=(SUM1-P1)/(1.0d0-SUM2)
        ENDIF
        IF(DOP) THEN
          WRITE(OP_STRING,'(/'' Salu consistency parameters:'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' SUM1  ='',D12.5)') SUM1
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' SUM2  ='',D12.5)') SUM2
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' P1    ='',D12.5)') P1
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' ALPHA ='',D12.5)') ALPHA
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
C       Update the solution vector
        DO no1=1,NOT(1,1,nr,nx)
          DO nyo2=1,NYNO(0,no1,2)
            ny2=NYNO(nyo2,no1,2)
            nk=NPNY(1,ny2,0)
            IF(nk.EQ.1) THEN
              XO(no1)=XO(no1)+ALPHA*YP(no1,6)
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_START2)

C*** Transfer solution vector to XO and YP
C CS new 16/12/98 initialise YP first
      DO no1=1,NOT(1,1,nr,nx)
        DO nyo2=1,NYNO(0,no1,2)
          ny2=NYNO(nyo2,no1,2)
          YP(ny2,1)=0.0d0
        ENDDO !nyo2
      ENDDO !no1

      DO no1=1,NOT(1,1,nr,nx)
        co1=CYNO(0,no1,2)
        DO nyo2=1,NYNO(0,no1,2)
          ny2=NYNO(nyo2,no1,2)
          co2=CYNO(nyo2,no1,2)
          IF(NPNY(0,ny2,0).EQ.1) THEN
            np=NPNY(4,ny2,0)
C GMH 8/1/97 Update cmgui link
            CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
          ENDIF
          YP(ny2,1)=XO(no1)*co2+co1+YP(ny2,1)
        ENDDO !nyo2
      ENDDO !no1

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
      IF(IWRIT4(nr,nx).GE.1) THEN
        WRITE(OP_STRING,'(/'' CPU time for storage: '',D11.4,'' s'')')
     '    ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START1(1)
      IF(IWRIT4(nr,nx).GE.1) THEN
        WRITE(OP_STRING,'(/'' Total CPU time for solution: '','
     '    //'D11.4,'' s'')') ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('SOLVE2')
      RETURN
 9999 CALL ERRORS('SOLVE2',ERROR)
      CALL EXITS('SOLVE2')
      RETURN 1
      END


