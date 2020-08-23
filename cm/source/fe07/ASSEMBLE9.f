      SUBROUTINE ASSEMBLE9(BC_POINTS,BRANCH,CQ,CONECT,ER,ES,GKK,
     '  GRR,HALF_TIME_STEP,ISC_GKK,ISR_GKK,LGE,M,N,nr,nx,NHQ,
     '  NQ_START,N_TERM_Q,NXQ,NYNQ,NYQNR,UPDATE_MATRIX,
     '  UPDATE_VECTOR,WORK_ARRAY,TIME,TOT_BITIME,TOT_BOTIME,
     '  TOT_GRTIME,XQ,YQ,TIME_VALUES,NTIME_POINTS,NTIME_NR,
     &  N_VENOUS_GEOM,VENOUS_NETWORK,ERROR,*)

C#### Subroutine: ASSEMBLE9
C###  Description:
C###    ASSEMBLE9 assembles the global reduced matrice GKK, and RHS
C###    vector GRR for implicit/explicit finite difference problems.
C###    NOTE: At present because there is no need to form the non
C###    sparse form of GK, which for an explicit finite difference
C###    scheme (the only kind implimented at present) would be
C###    prohibitivly large, GKK is assembled directed in a sparse
C###    form. This may have to change for coupled problems or implicit
C###    schemes. (NPS 30/10/96)

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'solv00.cmn'
      INCLUDE 'time02.cmn'

!     Parameter List
      INTEGER BC_POINTS(3,3,0:NQM),CONECT(-1:1,0:2,NQM),
     '  ISC_GKK(NISC_GKKM),ISR_GKK(NISR_GKKM),LGE(NHM*NSM,NRCM),
     '  NHQ(NRM),N_TERM_Q(-1:1,0:NEM),nr,nx,
     '  NXQ(-NIM:NIM,0:4,0:NQM,NAM),
     '  NYQNR(0:NYQM,0:NRCM,NCM,0:NRM),NYNQ(NHM,NQM,0:NRCM),
     '  NTIME_POINTS(NTIMEVARSM),NTIME_NR(0:NTIMEVARSM,NRM),
     &  N_VENOUS_GEOM
      REAL*8 CQ(NMM,NQM),ER(NHM*NSM),ES(NHM*NSM,NHM*NSM),
     '  GKK(NZ_GKK_M),GRR(NOM),TIME,XQ(NJM,NQM),YQ(NYQM,NIQM),
     '  TIME_VALUES(2,0:NTIMEPOINTSM+1,NTIMEVARSM)
      REAL TOT_BITIME,TOT_BOTIME,TOT_GRTIME
      CHARACTER ERROR*(*),VENOUS_NETWORK
      LOGICAL HALF_TIME_STEP,UPDATE_MATRIX,UPDATE_VECTOR
!     Local Variables
      INTEGER  ADJACENT,ii,M,N,nhs1,nhs2,NHST(2),
     '  no_nq,no_nynr,nq,nqq,NQ_START,ny1,ny2,nz,nq_term
      REAL  ELAPSED_TIME,TIME_START1(1),TIME_STOP(1)
      LOGICAL BRANCH(NQM),ERROR_FLAG,TERMINAL(NQM)
      LOGICAL*1 WORK_ARRAY(N,M)

      CALL ENTERS('ASSEMBLE9',*9999)


      IF(ITYP16(nr,nx).NE.4) THEN ! GKK,GRR not used for large explicit Lax
        IF(UPDATE_MATRIX) THEN    ! Wendroff problems
          IF(ITYP16(nr,nx).EQ.1) THEN
C set work array equal to identity for sparsity pattern for explict schemes
            DO ny1=1,M
              DO ny2=1,N
                IF(ny1.ne.ny2) THEN
                  WORK_ARRAY(ny2,ny1)=.FALSE.
                ELSE
                  WORK_ARRAY(ny2,ny1)=.TRUE.
                ENDIF
              ENDDO !ny2
            ENDDO !ny1
            CALL CALC_SPARSE(NISC_GKKM,NISR_GKKM,ISC_GKK,ISR_GKK,
     '        M,N,NZT(1,nx),SPARSEGKK(nx),WORK_ARRAY,ERROR,*9999)
          ELSE
C set work array equal to identity for sparsity pattern for
C implit schemes (not implimented)
          ENDIF !ITYP16

          DO nz=1,NZT(1,nx) !initalise GKK before updating
            GKK(nz)=0.0d0
          ENDDO !nz
        ENDIF !UPDATE_MATRIX
        IF(UPDATE_VECTOR) THEN
          DO no_nynr=1,NYQNR(0,1,1,nr)
            ny1=NYQNR(no_nynr,1,1,nr)
            GRR(ny1)=0.0d0 !initalise GR before updating
          ENDDO
        ENDIF !UPDATE_VECTOR

        IF(UPDATE_MATRIX.OR.UPDATE_VECTOR) THEN
          DO nq=NQR(1,nr),NQR(2,nr)
            CALL MELGEG(LGE,nq,NHST,NHQ,nr,NYNQ,ERROR,*9999)
            CALL XPFD30(NXQ,N_TERM_Q,CONECT,CQ,ES,HALF_TIME_STEP,
     '        NHQ,nq,nr,nx,NYNQ,UPDATE_MATRIX,UPDATE_VECTOR,TIME,XQ,
     '        YQ,ERROR,*9999)

            DO nhs1=1,NHST(1)
              ny1=IABS(LGE(nhs1,1))
              IF (UPDATE_MATRIX) THEN
                DO nhs2=1,NHST(2)
                  ny2=IABS(LGE(nhs2,2))
                  CALL SPARSE(ny1,ny2,MAX(M,N),nz,NZ_GK_M,
     '              NZT(1,nx),ISC_GKK,ISR_GKK,SPARSEGKK(nx),ERROR,*9999)
                  IF(nz.NE.0) THEN
                    GKK(nz)=GKK(nz)+ES(nhs1,nhs2)
                  ENDIF
                ENDDO !nhs2
              ENDIF
              IF (UPDATE_VECTOR) THEN
                GRR(ny1)=GRR(ny1)+ER(nhs1)
              ENDIF
            ENDDO !nhs1
          ENDDO
        ENDIF !UPDATE_MATRIX or UPDATE_VECTOR

        UPDATE_MATRIX=.FALSE.
        UPDATE_VECTOR=.FALSE.

      ELSE IF(ITYP16(nr,nx).EQ.4) THEN !Lax-Wendroff
        IF(ITYP3(nr,nx).EQ.1) THEN !flow in elastic tube
          CALL CPU_TIMER(CPU_USER,TIME_START1)
          ERROR_FLAG=.FALSE.

C PM 29-NOV-01 : Terminal points where transition occurs from arteries to veins
C Initialise variable
          DO nq=NQR(1,nr),NQR(2,nr)
            TERMINAL(nq)=.FALSE.
          ENDDO

          IF(((VENOUS_NETWORK.EQ.'Y').OR.(VENOUS_NETWORK.EQ.'y'))
     '      .AND.(N_VENOUS_GEOM.EQ.2)) THEN
            DO nq=NQR(1,nr),NQR(2,nr)
              DO nq_term=1,N_TERM_Q(0,0)
                IF(nq.EQ.N_TERM_Q(0,nq_term)) THEN
                  TERMINAL(nq)=.TRUE.
                ENDIF
              ENDDO
            ENDDO
          ENDIF

C$OMP   PARALLEL DO
C$OMP&  PRIVATE(nq,ADJACENT,ii,nqq),
C$OMP&  SHARED(NQR,NXQ,NQ_START,HALF_TIME_STEP,CONECT,
C$OMP&           CQ,ER,ES,NHQ,nr,nx,NYNQ,UPDATE_MATRIX,
C$OMP&           UPDATE_VECTOR,TERMINAL,
C$OMP&           TIME,XQ,YQ,BC_POINTS,BRANCH),
C$OMP&  SCHEDULE(GUIDED)
          DO nq=NQR(1,nr),NQR(2,nr) !ipinit NPS 4/2/97

            IF(.NOT.ERROR_FLAG) THEN

              ADJACENT=0
              DO ii=-1,1,2
                DO nqq=1,NXQ(ii,0,nq,1)
                  ADJACENT=ADJACENT+1
                ENDDO
              ENDDO
              IF(((nq.EQ.NQ_START).AND.
     '          HALF_TIME_STEP).OR.(ADJACENT.EQ.2)) THEN
C the grid point nq is not an end point or a bifrucation.
C grid points around a bifuraction are corrected latter

                IF(.NOT.TERMINAL(nq)) THEN
                  CALL XPFD30(NXQ,N_TERM_Q,CONECT,CQ,ES,
     '              HALF_TIME_STEP,NHQ,nq,nr,nx,NYNQ,UPDATE_MATRIX,
     '              UPDATE_VECTOR,TIME,XQ,YQ,ERROR,*100)
                ENDIF
              ENDIF

              IF(TERMINAL(nq)) THEN
                CALL BRANCH3(nq,nr,NXQ,NYNQ,CQ,XQ,YQ,HALF_TIME_STEP,
     '            ERROR,*100)
              ENDIF

              GO TO 102
C               This statement is designed to be skipped if no error
C               occur. However if a error occurs within a subroutine
C               the alternate return points to line 100 to set the
C               flag
 100            CONTINUE
C$OMP CRITICAL(ASSEMBLE9_1)
                ERROR_FLAG=.TRUE.
                WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
                CALL WRITES(IOER,OP_STRING,ERROR,*101)
                WRITE(OP_STRING,'(/'' >>An error occurred - '
     '            //'results may be unreliable!'')')
                CALL WRITES(IOER,OP_STRING,ERROR,*101)
 101            CONTINUE
C$OMP END CRITICAL(ASSEMBLE9_1)
 102          CONTINUE
            ENDIF !.NOT.ERROR_FLAG
          ENDDO !nq
C$OMP END PARALLEL DO

          CALL CPU_TIMER(CPU_USER,TIME_STOP)
          ELAPSED_TIME=TIME_STOP(1)-TIME_START1(1)
          CALL ASSERT(.NOT.ERROR_FLAG,'>>An error occurred during '
     '      //'Lax Wendroff grid calculations',ERROR,*9999)
          TOT_GRTIME=TOT_GRTIME+ELAPSED_TIME

          IF(.NOT.HALF_TIME_STEP) THEN
            CALL CPU_TIMER(CPU_USER,TIME_START1)

C$OMP   PARALLEL DO
C$OMP&  PRIVATE(no_nq),
C$OMP&  SHARED(nr,TIME),
C$OMP&  SCHEDULE(GUIDED)
            DO no_nq=1,BC_POINTS(1,1,0) !no_nq is bifurcation index
              IF(.NOT.ERROR_FLAG) THEN
                CALL BRANCH1(BC_POINTS(1,1,no_nq),BRANCH(no_nq),
     '            CONECT,CQ,ITYP12(nr,nx),nr,NTIME_POINTS,NTIME_NR,NYNQ,
     &            TIME,XQ,YQ,TIME_VALUES,ERROR,*110)

C SMAR009 18/01/99 removed  no_bc_points,
                GO TO 112
C                 This statement is designed to be skipped if no error
C                 occur. However if a error occurs within a subroutine
C                 the alternate return points to line 100 to set the flag
 110              CONTINUE
C$OMP CRITICAL(ASSEMBLE9_2)
                  ERROR_FLAG=.TRUE.
                  WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
                  CALL WRITES(IOER,OP_STRING,ERROR,*111)
                  WRITE(OP_STRING,'(/'' >>An error occurred - '
     '              //'results may be unreliable!'')')
                  CALL WRITES(IOER,OP_STRING,ERROR,*111)
 111            CONTINUE
C$OMP END CRITICAL(ASSEMBLE9_2)
 112            CONTINUE
              ENDIF !.NOT.ERROR_FLAG
            ENDDO !no_nq

C$OMP END PARALLEL DO

            CALL CPU_TIMER(CPU_USER,TIME_STOP)
            ELAPSED_TIME=TIME_STOP(1)-TIME_START1(1)
            TOT_BITIME=TOT_BITIME+ELAPSED_TIME
            CALL CPU_TIMER(CPU_USER,TIME_START1)

C$OMP PARALLEL DO
C$OMP&  PRIVATE(no_nq),
C$OMP&  SHARED(nr,TIME),
C$OMP&  SCHEDULE(GUIDED)
            DO no_nq=1,BC_POINTS(1,1,0)
              IF(.NOT.ERROR_FLAG) THEN
                CALL BRANCH2(BC_POINTS(1,1,no_nq),BRANCH,CONECT,
     &            ITYP12(nr,nx),CQ,no_nq,nr,NTIME_POINTS,NTIME_NR,NYNQ,
     &            TIME,TIME_VALUES,XQ,YQ,ERROR,*120)

C SMAR009 18/01/99 removed  no_bc_points,

                GO TO 122
C                 This statement is designed to be skipped if no error
C                 occurs. However if a error occurs within a subroutine
C                 the alternate return points to line 100 to set the
C                 flag
 120              CONTINUE
C$OMP CRITICAL(ASSEMBLE9_3)
                  ERROR_FLAG=.TRUE.
                  WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
                  CALL WRITES(IOER,OP_STRING,ERROR,*121)
                  WRITE(OP_STRING,'(/'' >>An error occurred - '
     '              //'results may be unreliable!'')')
                  CALL WRITES(IOER,OP_STRING,ERROR,*121)
 121              CONTINUE
C$OMP END CRITICAL(ASSEMBLE9_3)
 122            CONTINUE
              ENDIF !.NOT.ERROR_FLAG
            ENDDO !no_nq

C$OMP END PARALLEL DO

            CALL CPU_TIMER(CPU_USER,TIME_STOP)
            ELAPSED_TIME=TIME_STOP(1)-TIME_START1(1)
            TOT_BOTIME=TOT_BOTIME+ELAPSED_TIME

          ENDIF !.NOT.HALF_TIME_STEP

        ENDIF !ITYP3(nr,nx)=1

      ENDIF !ITYP16(nr,nx)

      CALL EXITS('ASSEMBLE9')
      RETURN
 9999 CALL ERRORS('ASSEMBLE9',ERROR)
      CALL EXITS('ASSEMBLE9')
      RETURN 1
      END



