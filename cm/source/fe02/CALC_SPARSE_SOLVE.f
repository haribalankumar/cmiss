      SUBROUTINE CALC_SPARSE_SOLVE(ISCMAX,ISRMAX,ISC,ISR,LGE,M,N,NBH,
     '  nc,NEELEM,NHE,NPNE,NPNY,NRLIST,NVHE,nx,NYNE,NYNP,NYNR,NZMAX,
     '  SPARSENESS,WORK_ARRAY,FIX,CALC_SPARSITY,ERROR,*)

C#### Subroutine: CALC_SPARSE_SOLVE
C###  Description:
C###    CALC_SPARSE_SOLVE calculates sparsity patterns for global
C###    matrices for region nr and solution problem type nx.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'solv00.cmn'
      INCLUDE 'time02.cmn'
!     Parameter List
      INTEGER ISCMAX,ISRMAX,ISC(ISCMAX),ISR(ISRMAX),LGE(NHM*NSM,NRCM),
     '  M,N,NBH(NHM,NCM,NEM),nc,NEELEM(0:NE_R_M,0:NRM),NHE(NEM),
     '  NPNE(NNM,NBFM,NEM),NPNY(0:6,NYM,0:NRCM),NRLIST(0:NRM),
     '  NVHE(NNM,NBFM,NHM,NEM),nx,NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM),NZMAX,SPARSENESS
      LOGICAL*1 WORK_ARRAY(N,M)
      LOGICAL FIX(NYM,NIYFIXM)
      CHARACTER ERROR*(*)
      LOGICAL CALC_SPARSITY
!     Local Variables
      INTEGER GETNYR,ne,nhs1,nhs2,NHST(2),noelem,no_nrlist,
     '  no_nynr1,no_nynr2,nr,ny1,ny2,ny3
      REAL ELAPSED_TIME,TIME_START(1),TIME_STOP(1)
      CHARACTER CHAR*1
!     EXTERNAL FUNCTION
      INTEGER IDIGITS

      CALL ENTERS('CALC_SPARSE_SOLVE',*9999)

      CALL CPU_TIMER(CPU_USER,TIME_START)

      IF(SPARSENESS.EQ.0) THEN !No sparsity
        NZT(nc,nx)=M*N
      ELSE IF(CALC_SPARSITY) THEN

        CALL ASSERT(USE_SPARSE.EQ.1,
     '    '>>Set USE_SPARSE to 1 to use sparsity arrays',ERROR,*9999)

C***    Initialise work array
C$OMP PARALLEL DO
C$OMP&  PRIVATE(ny2,ny1),
C$OMP&  SHARED(M,N,WORK_ARRAY)
        DO ny1=1,M
          DO ny2=1,N
            WORK_ARRAY(ny2,ny1)=.FALSE.
          ENDDO !ny2
        ENDDO !ny1
C$OMP END PARALLEL DO
C***    Calculate sparsity pattern
        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          IF(CALC_GLOBAL(nr,nx)) THEN
C***        Calculate sparsity pattern for global matrices (i.e. GK,
C***        GQ etc.)
            IF(ITYP4(nr,nx).EQ.1) THEN !FEM
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                CALL MELGE(LGE,NBH(1,1,ne),nc,ne,NHE(ne),NHST,
     '            NPNE(1,1,ne),nr,NVHE(1,1,1,ne),nx,NYNE,NYNP,
     '            ERROR,*9999)
                DO nhs1=1,NHST(1)
                  ny1=IABS(LGE(nhs1,1))
                  DO nhs2=1,NHST(2)
                    ny2=IABS(LGE(nhs2,2))
                    WORK_ARRAY(ny2,ny1)=.TRUE.
                  ENDDO !nhs2
                ENDDO !nhs1
              ENDDO !noelem
            ELSE IF(ITYP4(nr,nx).EQ.2) THEN !BEM
C cpb 20/5/95 This could be greatly improved as all rows in the same
C BEM region have the same sparsity pattern ie. the BEM sparsity is
C a block sparsity rather than generally sparse. Hence you do not
C have to calculate the sparsity like the FEM case rather you can just
C set it.
C$OMP PARALLEL DO
C$OMP&  PRIVATE(no_nynr1,ny1,no_nynr2,ny2),
C$OMP&  SHARED(NYNR,WORK_ARRAY)
              DO no_nynr1=1,NYNR(0,1,nc,nr)
                ny1=NYNR(no_nynr1,1,nc,nr)
                DO no_nynr2=1,NYNR(0,2,nc,nr)
                  ny2=NYNR(no_nynr2,2,nc,nr)
                  WORK_ARRAY(ny2,ny1)=.TRUE.
                ENDDO !no_nynr2
              ENDDO !no_nynr1
C$OMP END PARALLEL DO
            ELSE
              ERROR='>>Invalid ITYP4'
              GOTO 9999
            ENDIF
          ELSE
C***        Don't calculate sparsity pattern for global matrices i.e.
C***        assemble problem directly to GKK. For this case you still
C***        need to store the sparsity pattern for the fixed variables
C***        as you will need the GK, GQ values etc. to calculate fluxes
C***        etc.
            IF(ITYP4(nr,nx).EQ.1) THEN !FEM
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                CALL MELGE(LGE,NBH(1,1,ne),nc,ne,NHE(ne),NHST,
     '            NPNE(1,1,ne),nr,NVHE(1,1,1,ne),nx,NYNE,NYNP,
     '            ERROR,*9999)
                DO nhs1=1,NHST(1)
                  ny1=IABS(LGE(nhs1,1)) !global row #
                  ny2=GETNYR(1,NPNY,nr,0,1,ny1,NYNE,NYNP) !global var #
                  IF(FIX(ny2,1)) THEN
                    DO nhs2=1,NHST(2)
                      ny2=IABS(LGE(nhs2,2))
                      WORK_ARRAY(ny2,ny1)=.TRUE.
                    ENDDO !nhs2
                  ELSE
                    DO nhs2=1,NHST(2)
                      ny2=IABS(LGE(nhs2,2))
                      IF(FIX(ny2,1)) WORK_ARRAY(ny2,ny1)=.TRUE.
                    ENDDO !nhs2
                  ENDIF
                ENDDO !nhs1
              ENDDO !noelem
            ELSE IF(ITYP4(nr,nx).EQ.2) THEN !BEM
C$OMP         PARALLEL DO
C$OMP&          DEFAULT(NONE),
C$OMP&          PRIVATE(no_nynr1,no_nynr2,ny1,ny2,ny3),
C$OMP&          SHARED(nc,nr,NPNY,NYNE,NYNP,NYNR,FIX,WORK_ARRAY)
              DO no_nynr1=1,NYNR(0,1,nc,nr)
                ny1=NYNR(no_nynr1,1,nc,nr) !global row #
                ny2=GETNYR(nc,NPNY,nr,0,1,ny1,NYNE,NYNP) !global var #
                IF(FIX(ny2,1)) THEN
                  DO no_nynr2=1,NYNR(0,2,nc,nr)
                    ny2=NYNR(no_nynr2,2,nc,nr)
                    WORK_ARRAY(ny2,ny1)=.TRUE.
                  ENDDO !no_nynr2
                ELSE
                  DO no_nynr2=1,NYNR(0,2,nc,nr)
                    ny2=NYNR(no_nynr2,2,nc,nr) ! local var #
                    ny3=GETNYR(nc,NPNY,nr,0,2,ny2,NYNE,NYNP) !glob var #
                    IF(FIX(ny3,1)) WORK_ARRAY(ny2,ny1)=.TRUE.
                  ENDDO !nonynr2
                ENDIF
              ENDDO !no_nynr1
C$OMP         END PARALLEL DO
            ELSE
              ERROR='>>Invalid ITYP4'
              GOTO 9999
            ENDIF
          ENDIF
        ENDDO !no_nrlist (nr)

        CALL CALC_SPARSE(ISCMAX,ISRMAX,ISC,ISR,M,N,NZT(nc,nx),
     '    SPARSENESS,WORK_ARRAY,ERROR,*9999)

      ENDIF

C LKC 21-DEC-1998 Change assert so that it outputs the current NZMAX
C      CALL ASSERT(NZT(nc,nx).LE.NZMAX,'>>Increase NZMAX',ERROR,*9999)
      IF(NZT(nc,nx).GT.NZMAX) THEN
        WRITE(CHAR,'(I1)') IDIGITS(NZT(nc,nx))
        WRITE(ERROR,'(''Increase NZ_GK_M, try '',I'//CHAR//')')
     &    NZT(nc,nx)
        GOTO 9999
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)
      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(CALC_SPARSE_SOLVE_1)
        WRITE(OP_STRING,'('' CPU time to calculate sparsity pattern '
     '    //'for region '',I1,'' is '',D11.4,'' s'')') nr,ELAPSED_TIME
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(CALC_SPARSE_SOLVE_1)
      ENDIF

      CALL EXITS('CALC_SPARSE_SOLVE')
      RETURN
 9999 CALL ERRORS('CALC_SPARSE_SOLVE',ERROR)
      CALL EXITS('CALC_SPARSE_SOLVE')
      RETURN 1
      END


