      SUBROUTINE CALC_SPARSE_FIT(ISCMAX,ISRMAX,ISC,ISR,LGE,M,N,NBH,
     '  NEELEM,NPNE,NRLIST,NVHE,nx,NYNE,NYNP,SPARSENESS,WORK_ARRAY,
     '  CALC_SPARSITY,ERROR,*)

C#### Subroutine: CALC_SPARSE_FIT
C###  Description:
C###    CALC_SPARSE_FIT calculates sparsity patterns for global
C###    matrices for region nr and fitting problem type nx.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'fit000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'time02.cmn'
!     Parameter List
      INTEGER ISCMAX,ISRMAX,ISC(ISCMAX),ISR(ISRMAX),LGE(NHM*NSM,NRCM),
     '  M,N,NBH(NHM,NCM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NPNE(NNM,NBFM,NEM),NRLIST(0:NRM),NVHE(NNM,NBFM,NHM,NEM),nx,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),SPARSENESS
      LOGICAL*1 WORK_ARRAY(N,M)
      CHARACTER ERROR*(*)
      LOGICAL CALC_SPARSITY
!     Local Variables
      INTEGER IEND,ne,nhs1,nhs2,NHST(2),njj,
     '  noelem,no_nrlist,nr,ny1,ny2
      REAL ELAPSED_TIME,TIME_START(1),TIME_STOP(1)

      CALL ENTERS('CALC_SPARSE_FIT',*9999)

      CALL CPU_TIMER(CPU_USER,TIME_START)

      IF(SPARSENESS.EQ.0) THEN !No sparsity
        NZT(1,nx)=M*N
      ELSE IF(CALC_SPARSITY) THEN

        CALL ASSERT(USE_SPARSE.EQ.1,
     '    '>>Set USE_SPARSE to 1 to use sparsity arrays',ERROR,*9999)

C***    Initialise work array
        DO ny1=1,M
          DO ny2=1,N
            WORK_ARRAY(ny2,ny1)=.FALSE.
          ENDDO !no2
        ENDDO !no1
C***    Calculate sparsity pattern
        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            DO njj=1,NUM_FIT(0)
              CALL MELGEF(LGE,NBH(1,1,ne),ne,NHST,njj,NPNE(1,1,ne),nr,
     '          NVHE(1,1,1,ne),nx,NYNE,NYNP,ERROR,*9999)
              DO nhs1=1,NHST(1)
                ny1=IABS(LGE(nhs1,1))
                DO nhs2=1,NHST(2)
                  ny2=IABS(LGE(nhs2,2))
                  CALL ASSERT(ny2.LE.N,'>>exceeded working array ' // '
     '                 bounds N ',ERROR,*9999)
                  CALL ASSERT(ny1.LE.M,'>>exceeded working array ' // '
     '                 bounds M',ERROR,*9999)
                  WORK_ARRAY(ny2,ny1)=.TRUE.
                ENDDO !nhs2
              ENDDO !nh1
            ENDDO !njj
          ENDDO !noelem
        ENDDO !no_nrlist (nr)

        CALL CALC_SPARSE(ISCMAX,ISRMAX,ISC,ISR,M,N,NZT(1,nx),
     '    SPARSENESS,WORK_ARRAY,ERROR,*9999)

      ENDIF

C KAT 14Jan99:
      IF(NZT(1,nx).GT.NZ_GK_M) THEN
        IEND=0
        CALL APPENDC(IEND,' >>Increase NZ_GK_M to ',ERROR)
        CALL APPENDI(IEND,NZT(1,nx),ERROR)
        GOTO 9999
      ENDIF
C      CALL ASSERT(NZT(1,nx).LE.NZ_GK_M,'>>Increase NZ_GK_M',
C     '  ERROR,*9999)

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' CPU time to calculate sparsity pattern '
     '    //'for region '',I1,'' is '',D11.4,'' s'')') nr,ELAPSED_TIME
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('CALC_SPARSE_FIT')
      RETURN
 9999 CALL ERRORS('CALC_SPARSE_FIT',ERROR)
      CALL EXITS('CALC_SPARSE_FIT')
      RETURN 1
      END


