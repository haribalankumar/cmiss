      SUBROUTINE ASSEMBLE10(ISC_GKK,ISR_GKK,NEELEM,NENQ,NLATNE,NLATNQ,
     &  NLATPNQ,NLQ,NQGP,NQGP_PIVOT,NQNLAT,NQS,NQXI,NRLIST,NWQ,
     &  nx_ext,nx_trans,nx_upd,NXQ,AQ,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,GCHQ,
     &  GUQ,GKK,NQGW,PROPQ,XQ,BIDOMAIN,COUPBID,FIRST_A,FIXQ,
     &  IMPLICIT,UPDATE,UPDATEDT,UPDATE_MATRIX,SOLVEEIGHTPROBLEM,
     &  ERROR,*)

C#### Subroutine: ASSEMBLE10
C###  Description:
C###    ASSEMBLE10 writes directly into the compressed row storage
C###    arrays. It generates the constraint matrix for the implicit
C###    solution of grid activation problems.
C***  Created by Martin Buist, May 1997
C***  GBS August 2003 - archived all code for obsolete sparse solvers

      IMPLICIT NONE

      INCLUDE 'b01.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'solv00.cmn'
      INCLUDE 'time02.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter list
      INTEGER ISC_GKK(NISC_GKKM,NXM),ISR_GKK(NISR_GKKM,NXM),
     &  NEELEM(0:NE_R_M,0:NRM),NENQ(0:8,NQM),NLATNE(NEQM+1),
     &  NLATNQ(NEQM*NQEM),NLATPNQ(NQM),NLQ(NQM),NQGP(0:NQGM,NQM),
     &  NQGP_PIVOT(NQGM,NQM),NQNLAT(NEQM*NQEM),NQS(NEQM),
     &  NQXI(0:NIM,NQSCM),NRLIST(0:NRM),NWQ(8,0:NQM),nx_ext,nx_trans,
     &  nx_upd,NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      REAL*8 AQ(NMAQM,NQM),CQ(NMM,NQM),DNUDXQ(3,3,NQM),DXDXIQ(3,3,NQM),
     &  DXDXIQ2(3,3,NQM),GCHQ(3,NQM),GUQ(3,3,NQM),GKK(NZ_GKK_M,NXM),
     &  NQGW(NQGM,NQM),PROPQ(3,3,4,2,NQM),XQ(NJM,NQM)
      CHARACTER ERROR*(*)
      LOGICAL BIDOMAIN,COUPBID,FIRST_A,FIXQ(NYQM,NIYFIXM,NXM),IMPLICIT,
     &  SOLVEEIGHTPROBLEM,UPDATE,UPDATEDT,UPDATE_MATRIX

!     Local variables
      INTEGER DUMMY_LIST(0:1),NITB,nq,nr,nrr,nzero,nzz,
     &  nzz_full,PLACEINT,PLACEEXT
      REAL ELAPSED_TIME,TIME_START(1),TIME_STOP(1)
      REAL*8 COEFFSEXT(NQGM),DTCMAMTHETA,M(9,NQGM),THETACOEF
      LOGICAL NOTIMP

      CALL ENTERS('ASSEMBLE10',*9999)

      CALL ASSERT(IMPLICIT.OR.BIDOMAIN.OR.COUPBID,
     &  '>>Must specify either IMPLICIT, BIDOMAIN or COUPBID',
     &  ERROR,*9999)

      CALL CPU_TIMER(CPU_USER,TIME_START)      

C What is this statement for ?
      IF(THETA(1).LT.ZERO_TOL) THEN
        THETACOEF=1.0d0
      ELSE
        THETACOEF=THETA(1)
      ENDIF
    
      IF(.NOT.UPDATE) THEN
        IF(.NOT.COUPBID) THEN
C GBS 26-Aug-2003 Generate supporting points from NXI or lattice arrays
          IF(USE_LAT.EQ.0) THEN
            CALL GET_FD_POINTS(NEELEM,NENQ,NLQ,NQGP,NQGP_PIVOT,NQS,NQXI,
     &        NRLIST,NWQ,NXQ,ERROR,*9999)
          ELSE
            CALL GSUPPORT(NLATNE,NLATNQ,NLATPNQ,NQGP,NQGP_PIVOT,NQNLAT,
     &        NQS,NQXI,NRLIST,NWQ,ERROR,*9999)
          ENDIF
        ENDIF

C        CALL ASSERT(NOM.GE.NQT,'>>Increase NOM>=NQT',ERROR,*9999)
        PLACEINT=0
        PLACEEXT=0
        IF(IMPLICIT) THEN
          NOT(1,1,NRLIST(1),nx_trans)=NQT
          NOT(2,1,NRLIST(1),nx_trans)=NQT
          IF(SPARSEGKK(nx_trans).EQ.0) THEN
            CALL ASSEMBLE10_CHECK_NO_SPARSE(nzz_full,ERROR,*9999)
            NZZT(1,NRLIST(1),nx_trans)=nzz_full
C$OMP PARALLEL DO
C$OMP&  PRIVATE(nq),
C$OMP&  SHARED(GKK,nx_trans)
            DO nq=1,nzz_full
              GKK(nq,nx_trans)=0.0d0
            ENDDO !nq
C$OMP END PARALLEL DO
          ELSE
            CALL ASSEMBLE10_CHECK_SPARSE(NQGP,NRLIST,nzz,ERROR,*9999)
            NZZT(1,NRLIST(1),nx_trans)=nzz
            ISR_GKK(1,nx_trans)=1
            DO nrr=1,NRLIST(0)
              nr=NRLIST(nrr)
              DO nq=NQR(1,nr),NQR(2,nr)
                DO nzero=1,NQGP(0,nq)
                  PLACEINT=PLACEINT+1
                  ISC_GKK(PLACEINT,nx_trans)=NQGP(nzero,nq)
                ENDDO
                ISR_GKK(nq+1,nx_trans)=ISR_GKK(nq,nx_trans)+NQGP(0,nq)
              ENDDO !nq
            ENDDO !nr
          ENDIF
        ENDIF
        IF(BIDOMAIN) THEN
          NOT(1,1,NRLIST(1),nx_ext)=NQT
          NOT(2,1,NRLIST(1),nx_ext)=NQT
          IF(SPARSEGKK(nx_ext).EQ.0) THEN
            CALL ASSEMBLE10_CHECK_NO_SPARSE(nzz_full,ERROR,*9999)
            NZZT(1,NRLIST(1),nx_ext)=nzz_full
C$OMP PARALLEL DO
C$OMP&  PRIVATE(nq),
C$OMP&  SHARED(GKK,nx_ext)
            DO nq=1,NQT*NQT
              GKK(nq,nx_ext)=0.0d0
            ENDDO !nq
C$OMP END PARALLEL DO
          ELSE
            CALL ASSEMBLE10_CHECK_SPARSE(NQGP,NRLIST,nzz,ERROR,*9999)
            NZZT(1,NRLIST(1),nx_ext)=nzz
            ISR_GKK(1,nx_ext)=1
            DO nrr=1,NRLIST(0)
              nr=NRLIST(nrr)
              DO nq=NQR(1,nr),NQR(2,nr)
                DO nzero=1,NQGP(0,nq)
                  PLACEEXT=PLACEEXT+1
                  ISC_GKK(PLACEEXT,nx_ext)=NQGP(nzero,nq)
                ENDDO
                ISR_GKK(nq+1,nx_ext)=ISR_GKK(nq,nx_ext)+NQGP(0,nq)
              ENDDO !nq
            ENDDO !nr
          ENDIF
        ENDIF
        IF(COUPBID) THEN
          CALL ASSERT(USE_LAT.EQ.0,'>>Lattice method not implemented '
     &      //'for coupled bidomain',ERROR,*9999)
          NOT(1,1,NRLIST(1),nx_upd)=NQT
          NOT(2,1,NRLIST(1),nx_upd)=NQT
          IF(SPARSEGKK(nx_upd).EQ.0) THEN
            CALL ASSEMBLE10_CHECK_NO_SPARSE(nzz_full,ERROR,*9999)
            NZZT(1,NRLIST(1),nx_upd)=nzz_full
C$OMP PARALLEL DO
C$OMP&  PRIVATE(nq),
C$OMP&  SHARED(GKK,nx_upd)
            DO nq=1,NQT*NQT
              GKK(nq,nx_upd)=0.0d0
            ENDDO !nq
C$OMP END PARALLEL DO
          ELSE
            CALL ASSEMBLE10_CHECK_SPARSE(NQGP,NRLIST,nzz,ERROR,*9999)
            NZZT(1,NRLIST(1),nx_upd)=nzz
            ISR_GKK(1,nx_upd)=1
            DO nrr=1,NRLIST(0)
              nr=NRLIST(nrr)
              DO nq=NQR(1,nr),NQR(2,nr)
                DO nzero=1,NQGP(0,nq)
                  PLACEEXT=PLACEEXT+1
                  ISC_GKK(PLACEEXT,nx_ext)=NQGP(nzero,nq)
                ENDDO
                ISR_GKK(nq+1,nx_upd)=ISR_GKK(nq,nx_upd)+NQGP(0,nq)
              ENDDO !nq
            ENDDO !nr
          ENDIF
        ENDIF

      ENDIF !.NOT.UPDATE

      IF(UPDATEDT) THEN
        IF(IMPLICIT) THEN
          !Main loop
          DO nrr=1,NRLIST(0)
            nr=NRLIST(nrr)

C$OMP PARALLEL DO
C$OMP&  PRIVATE(DTCMAMTHETA,nq,nzero,PLACEINT),
C$OMP&  SHARED(ISR_GKK,GKK,NQGP,NQR,NQT,nr,NWQ,nx_trans,
C$OMP&    SPARSEGKK,CQ,DT)
            DO nq=NQR(1,nr),NQR(2,nr) !create line for each nq
              DTCMAMTHETA=DT/(CQ(1,nq)*CQ(2,nq))*THETACOEF
              IF(NWQ(1,nq).EQ.0) THEN !internal
                IF(SPARSEGKK(nx_trans).EQ.0) THEN !no sparsity
                  DO nzero=1,NQGP(0,nq)
                    GKK(((NQGP(nzero,nq)-1)*NQT)+nq,nx_trans)=
     &                -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTCMAMTHETA
                    IF(NQGP(nzero,nq).EQ.nq)
     &                GKK(((nq-1)*NQT)+nq,nx_trans)=
     &                GKK(((nq-1)*NQT)+nq,nx_trans)+1.0d0
                  ENDDO
                ELSE !compressed row
                  PLACEINT=ISR_GKK(nq,nx_trans)
                  DO nzero=1,NQGP(0,nq)
                    GKK(PLACEINT,nx_trans)=
     &                -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTCMAMTHETA
                    IF(NQGP(nzero,nq).EQ.nq) GKK(PLACEINT,nx_trans)=
     &                GKK(PLACEINT,nx_trans)+1.0d0
                    PLACEINT=PLACEINT+1
                  ENDDO
                ENDIF
              ENDIF
            ENDDO !nq
C$OMP END PARALLEL DO
          ENDDO !nr
        ENDIF
      ELSE
        NOTIMP=(.NOT.IMPLICIT)
        PLACEINT=0
        PLACEEXT=0

        !Main loop
        DO nrr=1,NRLIST(0)
          nr=NRLIST(nrr)
          NITB=NQXI(0,NQS(NEELEM(1,nr)))
          IF(IMPLICIT) THEN
            DO nq=NQR(1,nr),NQR(2,nr) !create line for each nq
c *** DPN 14 March 2000 - fixing up units
c              DTCMAMTHETA=DT/((CQ(1,nq)*1.0D-6)*CQ(2,nq))*THETACOEF
c              !10^-6 correction for Cm
              DTCMAMTHETA=DT/(CQ(1,nq)*CQ(2,nq))*THETACOEF
              IF(USE_LAT.EQ.0) THEN
                CALL CALC_GRID_COEF(NENQ,NITB,nq,NQGP,NQS,NQXI,
     &            NWQ(1,nq),NXQ(-NIM,0,0,1),nx_ext,nx_trans,COEFFSEXT,
     &            CQ,DNUDXQ,DXDXIQ,DXDXIQ2,GCHQ(1,nq),GUQ(1,1,nq),
     &            NQGW(1,nq),PROPQ(1,1,1,1,nq),.FALSE.,FIXQ,IMPLICIT,
     &            SOLVEEIGHTPROBLEM,ERROR,*9999)
              ELSE !lattice method
                CALL CALC_LATTICE_WEIGHTS(nq,NQGP(0,nq),M,XQ,
     &            ERROR,*9999)
                CALL CALC_LATT_COEF_TRANS(NITB,nq,NWQ(1,nq),NQGP(0,nq),
     &            AQ(1,nq),M,NQGW(1,nq),PROPQ,FIXQ(1,1,nx_trans),
     &            ERROR,*9999)
              ENDIF !USE_LAT.EQ.0
              DO nzero=1,NQGP(0,nq)
                IF(SPARSEGKK(nx_trans).EQ.0) THEN !no sparsity
                  PLACEINT=(NQGP(nzero,nq)-1)*NQT+nq
                ELSE
                  PLACEINT=PLACEINT+1
                ENDIF
                IF(NWQ(1,nq).EQ.0) THEN !internal
                  IF(USE_LAT.EQ.0) THEN
                    GKK(PLACEINT,nx_trans)=
     &                -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTCMAMTHETA
                  ELSE
                    GKK(PLACEINT,nx_trans)=-NQGW(nzero,nq)*DTCMAMTHETA
                  ENDIF
                  IF(NQGP(nzero,nq).EQ.nq) GKK(PLACEINT,nx_trans)=
     &              GKK(PLACEINT,nx_trans)+1.0d0
                ELSE !external
                  IF(USE_LAT.EQ.0) THEN
                    GKK(PLACEINT,nx_trans)=NQGW(NQGP_PIVOT(nzero,nq),nq)
                  ELSE
                    GKK(PLACEINT,nx_trans)=NQGW(nzero,nq)
                  ENDIF
                ENDIF
              ENDDO
            ENDDO !nq
          ENDIF !implicit

          IF(BIDOMAIN) THEN
            DO nq=NQR(1,nr),NQR(2,nr) !create line for each nq
              IF(USE_LAT.EQ.0) THEN
                CALL CALC_GRID_COEF(NENQ,NITB,nq,NQGP,NQS,NQXI,
     &            NWQ(1,nq),NXQ(-NIM,0,0,1),nx_ext,nx_trans,COEFFSEXT,
     &            CQ,DNUDXQ,DXDXIQ,DXDXIQ2,GCHQ(1,nq),GUQ(1,1,nq),
     &            NQGW(1,nq),PROPQ(1,1,1,1,nq),BIDOMAIN,FIXQ,NOTIMP,
     &            SOLVEEIGHTPROBLEM,ERROR,*9999)
              ELSE !lattice method
                CALL CALC_LATTICE_WEIGHTS(nq,NQGP(0,nq),M,XQ,
     &            ERROR,*9999)
                CALL CALC_LATT_COEF_EXT(NITB,nq,NWQ(1,nq),NQGP(0,nq),
     &            AQ(1,nq),COEFFSEXT,M,PROPQ,FIXQ(1,1,nx_ext),
     &            SOLVEEIGHTPROBLEM,ERROR,*9999)
              ENDIF
              DO nzero=1,NQGP(0,nq)
                IF(SPARSEGKK(nx_ext).EQ.0) THEN !no sparsity
                  PLACEEXT=(NQGP(nzero,nq)-1)*NQT+nq
                ELSE
                  PLACEEXT=PLACEEXT+1
                ENDIF
                IF(USE_LAT.EQ.0) THEN
                  GKK(PLACEEXT,nx_ext)=COEFFSEXT(NQGP_PIVOT(nzero,nq))
                ELSE
                  GKK(PLACEEXT,nx_ext)=COEFFSEXT(nzero)
                ENDIF
              ENDDO
            ENDDO !nq
          ENDIF

          IF(COUPBID) THEN
            DO nq=NQR(1,nr),NQR(2,nr) !create line for each nq
              CALL CALC_GRID_COEF(NENQ,NITB,nq,NQGP,NQS,NQXI,NWQ(1,nq),
     &          NXQ(-NIM,0,0,1),nx_upd,nx_trans,COEFFSEXT,CQ,DNUDXQ,
     &          DXDXIQ,DXDXIQ2,GCHQ(1,nq),GUQ(1,1,nq),NQGW(1,nq),
     &          PROPQ(1,1,1,1,nq),.TRUE.,FIXQ,IMPLICIT,
     &          SOLVEEIGHTPROBLEM,ERROR,*9999)
              DO nzero=1,NQGP(0,nq)
                IF(SPARSEGKK(nx_upd).EQ.0) THEN !no sparsity
                  PLACEEXT=(NQGP(nzero,nq)-1)*NQT+nq
                ELSE
                  PLACEEXT=PLACEEXT+1
                ENDIF
                IF(USE_LAT.EQ.0) THEN
                  GKK(PLACEEXT,nx_upd)=COEFFSEXT(NQGP_PIVOT(nzero,nq))
                ELSE
                  GKK(PLACEEXT,nx_upd)=COEFFSEXT(nzero)
                ENDIF
              ENDDO
            ENDDO !nq
          ENDIF
        ENDDO !nr

      ENDIF !updatedt

      IF(.NOT.COUPBID) THEN
        UP_GRID_MATERIAL=.FALSE.
        UP_GRID_TENSOR=.FALSE.
      ENDIF
      FIRST_A=.TRUE.
      UPDATE_MATRIX=.TRUE.

      IF(DOP) THEN
        nr=NRLIST(1)
        WRITE(OP_STRING,
     &    '(/'' Global stiffness matrix GKK - transmembrane:'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' NOT(1,1,nr,nx)='',I5,'
     &    //''', NOT(2,1,nr,nx)='',I5)') NOT(1,1,nr,nx_trans),
     &    NOT(2,1,nr,nx_trans)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        CALL OPSTFMAT(DUMMY_LIST,ISC_GKK,ISR_GKK,IOOP,
     &    NOT(1,1,nr,nx_trans),NOT(2,1,nr,nx_trans),
     &    NZZT(1,nr,nx_trans),DUMMY_LIST,SPARSEGKK(nx_trans),
     &    GKK,GKK,'GKK','GKK',.TRUE.,.TRUE.,.FALSE.,
     &    ERROR,*9999)
        IF(BIDOMAIN) THEN
          WRITE(OP_STRING,
     &      '(/'' Global stiffness matrix GKK - external:'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NOT(1,1,nr,nx)='',I5,'
     &      //''', NOT(2,1,nr,nx)='',I5)') NOT(1,1,nr,nx_ext),
     &      NOT(2,1,nr,nx_ext)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CALL OPSTFMAT(DUMMY_LIST,ISC_GKK,ISR_GKK,IOOP,
     &      NOT(1,1,nr,nx_ext),NOT(2,1,nr,nx_ext),
     &      NZZT(1,nr,nx_ext),DUMMY_LIST,SPARSEGKK(nx_ext),
     &      GKK,GKK,'GKK','GKK',.TRUE.,.TRUE.,.FALSE.,
     &      ERROR,*9999)
        ENDIF
      ENDIF
      
      IF(IWRIT5(NRLIST(1),nx_trans).GE.1) THEN
        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)
        WRITE(OP_STRING,'(1X,''Time for matrix assembly '',F8.2,'
     &    //'''s cpu'')') ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('ASSEMBLE10')
      RETURN
 9999 CALL ERRORS('ASSEMBLE10',ERROR)
      CALL EXITS('ASSEMBLE10')
      RETURN 1
      END


      SUBROUTINE ASSEMBLE10_CHECK_NO_SPARSE(nzz,ERROR,*)      

      IMPLICIT NONE

      INCLUDE 'geom00.cmn'

!     Parameter list
      INTEGER nzz
      CHARACTER ERROR*(*)

      CALL ENTERS('ASSEMBLE10_CHECK_NO_SPARSE',*9999)

      nzz=NQT*NQT
      IF(NZ_GKK_M.LT.nzz) THEN
        WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
        GOTO 9999
      ENDIF

      CALL EXITS('ASSEMBLE10_CHECK_NO_SPARSE')
      RETURN
 9999 CALL ERRORS('ASSEMBLE10_CHECK_NO_SPARSE',ERROR)
      CALL EXITS('ASSEMBLE10_CHECK_NO_SPARSE')
      RETURN 1
      END


      SUBROUTINE ASSEMBLE10_CHECK_SPARSE(NQGP,NRLIST,nzz,ERROR,*)      

      IMPLICIT NONE

      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'

!     Parameter list
      INTEGER NQGP(0:NQGM,NQM),NRLIST(0:NRM),nzz
      CHARACTER ERROR*(*)

!     Local variables
      INTEGER maxrow,nq,nr,nrr
      CALL ENTERS('ASSEMBLE10_CHECK_SPARSE',*9999)

      nzz=0
      maxrow=0
      DO nrr=1,NRLIST(0)
        nr=NRLIST(nrr)
        DO nq=NQR(1,nr),NQR(2,nr)
          nzz=nzz+NQGP(0,nq)
        ENDDO !nq
        IF(NQR(2,nr)+1.GT.maxrow) maxrow=NQR(2,nr)+1
      ENDDO !nr
      IF(NZ_GKK_M.LT.nzz) THEN
        WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
        GOTO 9999
      ELSE IF(NISC_GKKM.LT.nzz) THEN
        WRITE(ERROR,'(''>>Increase ISC_GKK_M to >= '',I12)') nzz
        GOTO 9999
      ELSE IF(NISR_GKKM.LT.maxrow) THEN
        WRITE(ERROR,'(''>>Increase ISR_GKK_M to >= '',I12)') maxrow
        GOTO 9999
      ENDIF

      CALL EXITS('ASSEMBLE10_CHECK_SPARSE')
      RETURN
 9999 CALL ERRORS('ASSEMBLE10_CHECK_SPARSE',ERROR)
      CALL EXITS('ASSEMBLE10_CHECK_SPARSE')
      RETURN 1
      END
