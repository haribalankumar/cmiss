      SUBROUTINE CALC_SPARSE_GKK(ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,ISR_GKK,
     '  ISR_GQ,LGE,NBH,NENP,NHE,M,N,NONY,NP_INTERFACE,NPNE,NPNY,nr,
     '  NRE,NVHE,nx,NYNE,NYNP,NYNR,GK,GQ,WORK_ARRAY,COUPLEDFEBE,
     '  STRIPZEROS,ERROR,*)

C#### Subroutine: CALC_SPARSE_GKK
C###  Description:
C###    CALC_SPARSE_GKK calculates sparsity patterns for global
C###    solution matrix GKK for region nr and problem type nx.
C###    If STRIPZEROS is .TRUE. then all elements of GK that are close
C###    to zero are not included in GKK; otherwise the sparsity pattern
C###    of GKK is just based on that of GK.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'coup00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'solv00.cmn'
      INCLUDE 'time02.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER ISC_GK(NISC_GKM),ISC_GKK(NISC_GKKM),ISC_GQ(NISC_GQM),
     '  ISR_GK(NISR_GKM),ISR_GKK(NISR_GKKM),ISR_GQ(NISR_GQM),
     '  LGE(NHM*NSM,NRCM),M,N,NBH(NHM,NCM,NEM),
     '  NENP(NPM,0:NEPM,0:NRM),NHE(NEM),NONY(0:NOYM,NYM,NRCM),
     '  NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),
     '  NPNY(0:6,NYM,0:NRCM),nr,NRE(NEM),
     '  NVHE(NNM,NBFM,NHM,NEM),nx,NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),NYNR(0:NY_R_M,0:NRCM,NCM)
      REAL*8 GK(NZ_GK_M),GQ(NZ_GQ_M)
      LOGICAL*1 WORK_ARRAY(N,M)
      CHARACTER ERROR*(*)
      LOGICAL COUPLEDFEBE,STRIPZEROS
!     Local Variables
      INTEGER GETNYR,ne,nhs1,nhs2,NHST(2),np,no1,no2,none,nonp,
     '  no_nynr1,no_nynr2,noy1,noy2,nr2,nrr,ntmp,ny1,ny2,ny3,ny4,ny5,
     '  nyy1,nyy2,nyy3,nzgk,nzgq,ROWLIST(0:1)
      REAL ELAPSED_TIME,TIME_START(1),TIME_STOP(1)
      REAL*8 A(1)
      LOGICAL BEM_REGION,ERROR_FLAG,INLIST,INROW,VALID_REGION

      CALL ENTERS('CALC_SPARSE_GKK',*9999)

      CALL CPU_TIMER(CPU_USER,TIME_START)

      IF(SPARSEGKK(nx).EQ.0) THEN !No sparsity
C***    Do nothing
        NZZT(1,nr,nx)=M*N
      ELSE

C***    Initialise work array
C$OMP PARALLEL DO
C$OMP&  PRIVATE(no1,no2),
C$OMP&  SHARED(M,N,WORK_ARRAY)
        DO no1=1,M
          DO no2=1,N
            WORK_ARRAY(no2,no1)=.FALSE.
          ENDDO !no1
        ENDDO !no2
C$OMP END PARALLEL DO
        IF(nr.EQ.0) THEN
          BEM_REGION=COUPLEDFEBE
        ELSE
          BEM_REGION=(ITYP4(nr,nx).EQ.2)
        ENDIF
C***    Calculate sparsity
        IF(CALC_GLOBAL(nr,nx)) THEN
          ERROR_FLAG=.FALSE.
C$OMP     PARALLEL DO
C$OMP&      PRIVATE(no_nynr1,ny1,noy1,no1,no_nynr2,ny2,ny3,ny4,ny5,
C$OMP&      nzgk,nzgq,noy2,no2),
C$OMP&      SHARED(ERROR_FLAG,KTYP24,nr,NZ_GK_M,NZ_GQ_M,WORK_ARRAY)
          DO no_nynr1=1,NYNR(0,1,1) !Loop over the rows of GK
            IF(.NOT.ERROR_FLAG) THEN
              ny1=NYNR(no_nynr1,1,1) !is the row #
              IF(KTYP24.EQ.1.AND..NOT.BEM_REGION) THEN
C               Compressed Row GK matrix and not BEM.
C               Dont need to loop over all columns.
                DO nzgk=ISR_GK(ny1),ISR_GK(ny1+1)-1 !loop over nz cols
C                                                    of GK
                  ny3=ISC_GK(nzgk) !local column #
                  ny2=GETNYR(1,NPNY,nr,0,2,ny3,NYNE,NYNP) !glob column #
                  DO noy1=1,NONY(0,ny1,1) !Loop over solu dofs
C                                          attached to ny1
                    no1=NONY(noy1,ny1,1) !solution row # attached to ny1
                    IF(.NOT.STRIPZEROS.OR.
     '                DABS(GK(nzgk)).GE.ZERO_TOL) THEN
                      DO noy2=1,NONY(0,ny2,2) !Loop sol dofs from ny2
                        no2=NONY(noy2,ny2,2) !sol var # attached to ny2
                        WORK_ARRAY(no2,no1)=.TRUE.
                      ENDDO !noy2
                    ENDIF
                  ENDDO !noy1
                ENDDO !nz
              ELSE
                DO noy1=1,NONY(0,ny1,1) !Loop over solu dofs attached
C                                        to ny1
                  no1=NONY(noy1,ny1,1) !is the solu row attached to ny1
                  DO no_nynr2=1,NYNR(0,0,1) !Loop over the columns of GK
                    ny2=NYNR(no_nynr2,0,1) !global var number
                    ny3=GETNYR(1,NPNY,nr,2,0,ny2,NYNE,NYNP) !is local #
                    CALL SPARSE(ny1,ny3,NYT(1,1,nx),nzgk,NZ_GK_M,
     '                NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*100)
                    IF(nzgk.NE.0) THEN
                      IF(.NOT.STRIPZEROS.OR.
     '                  DABS(GK(nzgk)).GE.ZERO_TOL) THEN
                        DO noy2=1,NONY(0,ny2,2) !Loop sol dofs from ny2
                          no2=NONY(noy2,ny2,2)
                          WORK_ARRAY(no2,no1)=.TRUE.
                        ENDDO !noy2
                      ENDIF
                    ENDIF
                    IF(BEM_REGION) THEN
                      ny4=GETNYR(2,NPNY,nr,0,0,ny2,NYNE,NYNP) !glo rhs #
                      ny5=GETNYR(2,NPNY,nr,2,0,ny4,NYNE,NYNP) !loc rhs #
                      CALL SPARSE(ny1,ny5,NYT(1,2,nx),nzgq,NZ_GQ_M,
     '                  NZT(2,nx),ISC_GQ,ISR_GQ,KTYP24,ERROR,*100)
                      IF(nzgq.NE.0) THEN
                        IF(.NOT.STRIPZEROS.OR
     '                    .DABS(GQ(nzgq)).GE.ZERO_TOL) THEN
                          DO noy2=1,NONY(0,ny4,2)
                            no2=NONY(noy2,ny4,2)
                            WORK_ARRAY(no2,no1)=.TRUE.
                          ENDDO !noy2
                        ENDIF
                      ENDIF
                    ENDIF !BEM_REGION
                  ENDDO !no_nynr2 (ny2)
                ENDDO !noy1
              ENDIF
              GOTO 102
C             This statement is designed to be skipped if no error
C             occurs. However if a error occurs within a subroutine
C             the alternate return points to line 100 to set the flag
 100          CONTINUE
C$OMP         CRITICAL(CALC_SPARSE_GKK_1)
              ERROR_FLAG=.TRUE.
              WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
              CALL WRITES(IOER,OP_STRING,ERROR,*101)
              WRITE(OP_STRING,'(/'' >>An error occurred - '
     '          //'results may be unreliable!'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*101)
 101          CONTINUE
C$OMP         END CRITICAL(CALC_SPARSE_GKK_1)
 102          CONTINUE
            ENDIF
          ENDDO !no_nynr1 (ny1)
C$OMP END PARALLEL DO
          IF(ERROR_FLAG) GOTO 9999
        ELSE
          DO no_nynr1=1,NYNR(0,1,1) !Loop over the rows of GK
            ny1=NYNR(no_nynr1,1,1) !row #
            DO noy1=1,NONY(0,ny1,1) !Loop over sol dofs attach to ny1
              no1=NONY(noy1,ny1,1)
              DO no_nynr2=1,NYNR(0,0,1) !Loop over the cols of GK
                ny2=NYNR(no_nynr2,0,1) !var #
C new MPN 22Sep1999: fixing bad programming
                IF(NPNY(0,ny2,0).EQ.1) THEN
                  nrr=NPNY(6,ny2,0) !node based ny
                ELSE IF(NPNY(0,ny2,0).EQ.2) THEN
                  nrr=NPNY(5,ny2,0) !element based ny
                ENDIF
C old            nrr=NPNY(6,ny2,0) !assuming node base ny's!!!!!
                IF(ITYP4(nrr,nx).EQ.1) THEN !FEM variable
                  INROW=.FALSE.
                  np=NPNY(4,ny1,1)
                  none=1
                  DO WHILE(none.LE.NENP(np,0,0).AND..NOT.INROW)
                    ne=NENP(np,none,0)
                    IF(IS_COUPLED(nx)) THEN
                      IF(INLIST(NRE(ne),COUP_NRLIST(1,nx),
     &                  COUP_NRLIST(0,nx),ntmp)) THEN
                        VALID_REGION=.TRUE.
                      ELSE
                        VALID_REGION=.FALSE.
                      ENDIF
                    ELSE
                      VALID_REGION=.TRUE.
                    ENDIF
                    IF(VALID_REGION) THEN
                      CALL MELGE(LGE,NBH(1,1,ne),1,ne,NHE(ne),NHST,
     '                  NPNE(1,1,ne),NRE(ne),NVHE(1,1,1,ne),nx,NYNE,
     '                  NYNP,ERROR,*9999)
                      nhs1=1
                      DO WHILE(nhs1.LE.NHST(1).AND..NOT.INROW)
                        nyy1=IABS(LGE(nhs1,1))
                        IF(nyy1.EQ.ny1) THEN
                          nhs2=1
                          DO WHILE(nhs2.LE.NHST(2).AND..NOT.INROW)
                            nyy2=IABS(LGE(nhs2,2))
                            nyy3=GETNYR(1,NPNY,nr,0,2,nyy2,NYNE,NYNP)
                            IF(nyy3.EQ.ny2) INROW=.TRUE.
                            nhs2=nhs2+1
                          ENDDO !nhs2
                        ENDIF
                        nhs1=nhs1+1
                      ENDDO !nhs1
                    ENDIF !coupled region
                    none=none+1
                  ENDDO !none (ne)
                  IF(INROW) THEN
                    DO noy2=1,NONY(0,ny2,2)
                      no2=NONY(noy2,ny2,2)
                      WORK_ARRAY(no2,no1)=.TRUE.
                    ENDDO !noy2 (no2)
                  ENDIF
                ELSE IF(ITYP4(nrr,nx).EQ.2) THEN !BEM variable
                  INROW=.FALSE.
                  nr2=NPNY(6,ny1,1) !assuming node base ny's!!!!!
                  np=NPNY(4,ny2,0)
                  DO nonp=1,NP_INTERFACE(np,0)
                    IF(NP_INTERFACE(np,nonp).EQ.nr2) INROW=.TRUE.
                  ENDDO
                  IF(INROW) THEN
                    DO noy2=1,NONY(0,ny2,2)
                      no2=NONY(noy2,ny2,2)
                      WORK_ARRAY(no2,no1)=.TRUE.
                    ENDDO !noy2 (no2)
                    ny4=GETNYR(2,NPNY,nr,0,0,ny2,NYNE,NYNP) !rhs var#
                    DO noy2=1,NONY(0,ny4,2)
                      no2=NONY(noy2,ny4,2)
                      WORK_ARRAY(no2,no1)=.TRUE.
                    ENDDO !noy2
                  ENDIF
                ENDIF
              ENDDO
            ENDDO !noy1 (no1)
          ENDDO !no_nynr1 (ny1)
        ENDIF

        CALL CALC_SPARSE(NISC_GKKM,NISR_GKKM,ISC_GKK,ISR_GKK,M,N,
     '    NZZT(1,nr,nx),SPARSEGKK(nx),WORK_ARRAY,ERROR,*9999)

      ENDIF

      WRITE(OP_STRING(1),'(''>>Increase NZ_GKK_M >= '',I12)') NZZT(1,nr,
     &  nx)
      CALL ASSERT(NZZT(1,nr,nx).LE.NZ_GKK_M,OP_STRING(1),ERROR,*9999)

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)
      IF(DOP) THEN
        WRITE(OP_STRING,'('' CPU time to calculate sparsity pattern '
     '    //'for region '',I1,'' is '',D11.4,'' s'')') nr,ELAPSED_TIME
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(/'' Total solution matrix size, NZZT='',I10)')
     '    NZZT(1,nr,nx)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        IF(SPARSEGKK(nx).EQ.0) THEN
          WRITE(OP_STRING,'(/'' Solution matrix is stored as a '
     '      //'fully populated matrix'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ELSE IF(SPARSEGKK(nx).EQ.1) THEN
          WRITE(OP_STRING,'(/'' Solution matrix is stored as a '
     '      //'compressed row sparse matrix'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ELSE IF(SPARSEGKK(nx).EQ.2) THEN
          WRITE(OP_STRING,'(/'' Solution matrix is stored as a '
     '      //'row-column sparse matrix'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ELSE IF(SPARSEGKK(nx).EQ.3) THEN
          WRITE(OP_STRING,'(/'' Solution matrix is stored as a '
     '      //'compressed column sparse matrix'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ELSE IF(SPARSEGKK(nx).EQ.4) THEN
          WRITE(OP_STRING,'(/'' Solution matrix is stored as a '
     '      //'sorted row-column sparse matrix'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ELSE IF(SPARSEGKK(nx).EQ.5) THEN
          WRITE(OP_STRING,'(/'' Solution matrix is stored as a '
     '      //'Umfpack row-column sparse matrix'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        CALL IO_MATRIX('WRITE','TERM',ISC_GKK,ISR_GKK,IODI,
     '    NOT(1,1,nr,nx),NOT(2,1,nr,nx),NOM,NZZT(1,nr,nx),NOM,
     '    NOM,NZ_GKK_M,ROWLIST,
     '    SPARSEGKK(nx),A,'ASCII ','SPARSITY ',.TRUE.,ERROR,*9999)
      ENDIF

      CALL EXITS('CALC_SPARSE_GKK')
      RETURN
 9999 CALL ERRORS('CALC_SPARSE_GKK',ERROR)
      CALL EXITS('CALC_SPARSE_GKK')
      RETURN 1
      END


