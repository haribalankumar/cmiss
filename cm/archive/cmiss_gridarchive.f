C#### Module: CMISS_GRIDARCHIVE
C###  Description:
C###    Old grid related cmiss routines


C FE07
C=====
C###  Routine: ADAMS_ASSEMBLE10 cpb merge
C###  Routine: ADAMS_MARCH8 cpb merge
C###  Routine: ASSEMBLE10     - pre-removal of non-compressed-row sparsity
C###  Routine: ASSEMBLE10_FE  - pre-removal of non-compressed-row sparsity
C###  Routine: MARCH3    cardiac activ.n problems (fixed time step)
C###  Routine: MARCH5    cardiac activ.n problems (bidomain)
C###  Routine: MARCH7    implicit grid solution

C FE19
C=====
C###  Routine: BR_ION        Beeler-Reuter ionic equation
C###  Routine: BR_RATES  Beeler-Reuter rate constants  
C###  Routine: DN_ION    (fn)diFrancesco-Noble ionic equation  
C###  Routine: FHN_ION   (fn)FitzFugh-Nagumo ionic equation
C###  Routine: LR_ION    (fn)Luo-Rudy ionic equation
C###  Routine: LR_RATES  Luo-Rudy rate constants  
C###  Routine: VCD_ION   (fn)van Capelle-Durrer ionic equation  

C FE30
C=====
C###  Routine: BFRONT    activ. pattern with fixed t step-Bidomain
C###  Routine: CALC_ADAMS_DTAR cpb merge 
C###  Routine: CALC_ADAMS_EXT_RHS cpb merge
C###  Routine: CALC_ADAMS_GRID_COEF cpb merge
C###  Routine: CALC_ADAMS_PSEUDO_STIMULUS cpb merge
C###  Routine: CALC_GRID_BOUND_COEF old version
C###  Routine: GEN_INT_RHS rhs vector for Vm implicit grid
C###  Routine: GET_EXT_CONTRIB phi-e component of Vm equation
C###  Routine: IONIC_CURRENT calculates ionic current info for grids
C###  Routine: TFRONT    activ. pattern with fixed t step-threshold


C FE07
C=====

      SUBROUTINE ADAMS_ASSEMBLE10(ISC_GKK,ISR_GKK,NEELEM,NQGP,
     '  NQGP_PIVOT,NQS,NQXI,NRLIST,NWQ,nx_ext,nx_trans,
     '  nx_upd,NXQ,CQ,GCHQ,GUQ,GKK,NEWDT,NQGW,OLDDT,PROPQ,BIDOMAIN,
     '  FIRST_A,FIXQ,IMPLICIT,UPDATE,UPDATEDT,UPDATE_MATRIX,ERROR,*)

C#### Subroutine: ADAMS_ASSEMBLE10
C###  Description:
C###    ADAMS_ASSEMBLE10 creates the stiffness matrices for the 
C###    adaptive time stepping ionic current grid problems that
C###    use the Adams-Moulton integrator.

      IMPLICIT NONE

      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b12.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:grid00.cmn'
      INCLUDE 'cmiss$reference:iwrit00.cmn'
      INCLUDE 'cmiss$reference:solv00.cmn'
      INCLUDE 'cmiss$reference:time02.cmn'
      INCLUDE 'cmiss$reference:tol00.cmn'

!     Parameter list
      INTEGER ISC_GKK(NISC_GKKM,NXM),ISR_GKK(NISR_GKKM,NXM),
     '  NEELEM(0:NE_R_M,0:NRM),NQGP(0:19,NQM),NQGP_PIVOT(19,NQM),
     '  NQS(NEM),NQXI(0:NIM,NQSCM),NRLIST(0:NRM),NWQ(6,0:NQM),nx_ext,
     '  nx_trans,nx_upd,NXQ(-NIM:NIM,0:4,0:NQM)
      REAL*8 CQ(NMM,NQM),GCHQ(3,NQM),GUQ(3,3,NQM),
     '  GKK(NZ_GKK_M,NXM),NEWDT,NQGW(19,NQM),OLDDT,PROPQ(3,3,4,2,NQM)
      CHARACTER ERROR*(*)
      LOGICAL BIDOMAIN,FIRST_A,FIXQ(NYQM,NIYFIXM,NXM),IMPLICIT,UPDATE,
     '  UPDATEDT,UPDATE_MATRIX
!     Local variables
      INTEGER COUNTINT,COUNTEXT,DUMMY_LIST(0:1),IJ,IK,maxrow,
     '  maxcol,nii,nij,nik,NITB,nq,nzero2,nq1,nq2,nr,nrr,
     '  nzero,nzz,PLACEINT,PLACEEXT,PLACE_PIVOT(19),TEMPPLACE(19)
      REAL*8 COEFFSEXT(19),DTCMAMTHETA,TEMPCOEFF(19),THETACOEF
      REAL ELAPSED_TIME,TIME_START(1),TIME_STOP(1)
      LOGICAL NOTIMP

      CALL ENTERS('ADAMS_ASSEMBLE10',*9999)

      CALL CPU_TIMER(CPU_USER,TIME_START)
      
      IF(.NOT.UPDATE) THEN
        CALL GET_FD_POINTS(NEELEM,NQGP,NQGP_PIVOT,NQS,NQXI,NRLIST,NWQ,
     '    NXQ,ERROR,*9999)
      ENDIF

      IF(THETA(1).LT.ZERO_TOL) THEN
        THETACOEF=1.0d0
      ELSE
        THETACOEF=THETA(1)
      ENDIF
      
C     Check array sizes
      IF(.NOT.UPDATE) THEN
        IF(IMPLICIT) THEN
          NOT(1,1,NRLIST(1),nx_trans)=NQT
          NOT(2,1,NRLIST(1),nx_trans)=NQT
          IF(SPARSEGKK(nx_trans).EQ.0) THEN !no sparsity
            nzz=NQT*NQT
            NZZT(1,NRLIST(1),nx_trans)=nzz
            IF(NZ_GKK_M.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ENDIF
          ELSE IF(SPARSEGKK(nx_trans).EQ.1) THEN !compressed row
            nzz=0
            maxrow=0
            DO nrr=1,NRLIST(0)
              nr=NRLIST(nrr)
              DO nq=NQR(1,nr),NQR(2,nr)
                DO nzero=1,NQGP(0,nq)
                  nzz=nzz+1
                ENDDO !nzero
              ENDDO !nq
              IF(NQR(2,nr)+1.GT.maxrow) maxrow=NQR(2,nr)+1
              NZZT(1,nr,nx_trans)=nzz
            ENDDO !nr
            IF(NZ_GKK_M.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ELSE IF(NISC_GKKM.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase ISC_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ELSE IF(NISR_GKKM.LT.maxrow) THEN
              WRITE(ERROR,'(''>>Increase ISR_GKK_M to >= '',I12)') 
     '          maxrow
              GOTO 9999
            ENDIF
            ISR_GKK(1,nx_trans)=1
          ELSE IF(SPARSEGKK(nx_trans).EQ.2) THEN !row/col
            nzz=0
            DO nrr=1,NRLIST(0)
              nr=NRLIST(nrr)
              DO nq=NQR(1,nr),NQR(2,nr)
                DO nzero=1,NQGP(0,nq)
                  nzz=nzz+1
                ENDDO !nzero
              ENDDO !nq
              NZZT(1,nr,nx_trans)=nzz
            ENDDO !nr
            IF(NZ_GKK_M.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ELSE IF(NISR_GKKM.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase ISR_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ELSE IF(NISC_GKKM.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase ISC_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ENDIF
          ELSE IF(SPARSEGKK(nx_trans).EQ.3) THEN !compressed column
            nzz=0
            maxcol=0
            DO nrr=1,NRLIST(0)
              nr=NRLIST(nrr)
              DO nq=NQR(1,nr),NQR(2,nr)
                DO nzero=1,NQGP(0,nq)
                  nzz=nzz+1
                ENDDO !nzero
              ENDDO !nq
              IF(NQR(2,nr)+1.GT.maxcol) maxcol=NQR(2,nr)+1
              NZZT(1,nr,nx_trans)=nzz
            ENDDO !nr
            IF(NZ_GKK_M.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ELSE IF(NISR_GKKM.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase ISR_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ELSE IF(NISC_GKKM.LT.maxcol) THEN
              WRITE(ERROR,'(''>>Increase ISC_GKK_M to >= '',I12)') 
     '          maxcol
              GOTO 9999
            ENDIF
          ELSE IF(SPARSEGKK(nx_trans).EQ.4) THEN
            ERROR='>>Not implemented'
            GOTO 9999
          ELSE IF(SPARSEGKK(nx_trans).EQ.5) THEN !row/col #2
            nzz=0
            DO nrr=1,NRLIST(0)
              nr=NRLIST(nrr)
              DO nq=NQR(1,nr),NQR(2,nr)
                DO nzero=1,NQGP(0,nq)
                  nzz=nzz+1
                ENDDO !nzero
              ENDDO !nq
              NZZT(1,nr,nx_trans)=nzz
            ENDDO !nr
            IF(NZ_GKK_M.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ELSE IF(NISC_GKKM.LT.2*nzz) THEN
              WRITE(ERROR,'(''>>Increase ISC_GKK_M to >= '',I12)') 2*nzz
              GOTO 9999
            ENDIF
          ELSE
            ERROR='>>Unknown sparsity type for GKK'
            GOTO 9999
          ENDIF
        ENDIF
        IF(BIDOMAIN) THEN
          NOT(1,1,NRLIST(1),nx_ext)=NQT
          NOT(2,1,NRLIST(1),nx_ext)=NQT
          IF(SPARSEGKK(nx_ext).EQ.0) THEN !no sparsity
            nzz=NQT*NQT
            NZZT(1,NRLIST(1),nx_ext)=nzz
            IF(NZ_GKK_M.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ENDIF
          ELSE IF(SPARSEGKK(nx_ext).EQ.1) THEN !compressed row
            nzz=0
            maxrow=0
            DO nrr=1,NRLIST(0)
              nr=NRLIST(nrr)
              DO nq=NQR(1,nr),NQR(2,nr)
                DO nzero=1,NQGP(0,nq)
                  nzz=nzz+1
                ENDDO !nzero
              ENDDO !nq
              NZZT(1,nr,nx_ext)=nzz
              IF(NQR(2,nr)+1.GT.maxrow) maxrow=NQR(2,nr)+1
            ENDDO !nr
            IF(NZ_GKK_M.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ELSE IF(NISC_GKKM.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase ISC_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ELSE IF(NISR_GKKM.LT.maxrow) THEN
              WRITE(ERROR,'(''>>Increase ISR_GKK_M to >= '',I12)') 
     '          maxrow
              GOTO 9999
            ENDIF
          ELSE IF(SPARSEGKK(nx_ext).EQ.2) THEN !row/col
            nzz=0
            DO nrr=1,NRLIST(0)
              nr=NRLIST(nrr)
              DO nq=NQR(1,nr),NQR(2,nr)
                DO nzero=1,NQGP(0,nq)
                  nzz=nzz+1
                ENDDO !nzero
              ENDDO !nq
              NZZT(1,nr,nx_ext)=nzz
            ENDDO !nr
            IF(NZ_GKK_M.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ELSE IF(NISR_GKKM.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase ISR_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ELSE IF(NISC_GKKM.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase ISC_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ENDIF
          ELSE IF(SPARSEGKK(nx_ext).EQ.3) THEN !compressed column
            nzz=0
            maxcol=0
            DO nrr=1,NRLIST(0)
              nr=NRLIST(nrr)
              DO nq=NQR(1,nr),NQR(2,nr)
                DO nzero=1,NQGP(0,nq)
                  nzz=nzz+1
                ENDDO !nzero
              ENDDO !nq
              IF(NQR(2,nr)+1.GT.maxcol) maxcol=NQR(2,nr)+1
              NZZT(1,nr,nx_ext)=nzz
            ENDDO !nr
            IF(NZ_GKK_M.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ELSE IF(NISR_GKKM.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase ISR_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ELSE IF(NISC_GKKM.LT.maxcol) THEN
              WRITE(ERROR,'(''>>Increase ISC_GKK_M to >= '',I12)') 
     '          maxcol
              GOTO 9999
            ENDIF
          ELSE IF(SPARSEGKK(nx_ext).EQ.4) THEN
            ERROR='>>Not implemented'
            GOTO 9999
          ELSE IF(SPARSEGKK(nx_ext).EQ.5) THEN !row/col #2
            nzz=0
            DO nrr=1,NRLIST(0)
              nr=NRLIST(nrr)
              DO nq=NQR(1,nr),NQR(2,nr)
                DO nzero=1,NQGP(0,nq)
                  nzz=nzz+1
                ENDDO !nzero
              ENDDO !nq
              NZZT(1,nr,nx_ext)=nzz
            ENDDO !nr
            IF(NZ_GKK_M.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ELSE IF(NISC_GKKM.LT.2*nzz) THEN
              WRITE(ERROR,'(''>>Increase ISC_GKK_M to >= '',I12)') 2*nzz
              GOTO 9999
            ENDIF
          ELSE
            ERROR='>>Unknown sparsity type for GKK'
            GOTO 9999
          ENDIF
        ENDIF
        CALL ASSERT(NOM.GE.NQT,'>>Increase NOM, min. NQT',ERROR,*9999)
      ENDIF

C     Initialising
      IF(.NOT.UPDATE) THEN
        IF(IMPLICIT) THEN
          IF(SPARSEGKK(nx_trans).EQ.0) THEN
C$OMP PARALLEL DO
C$&     PRIVATE(nq),
C$&     SHARED(GKK,nx_trans)
            DO nq=1,NQT*NQT
              GKK(nq,nx_trans)=0.0d0
            ENDDO !nq
C$OMP END PARALLEL DO
          ELSE IF(SPARSEGKK(nx_trans).EQ.1) THEN
            ISR_GKK(1,nx_trans)=1
          ELSE IF(SPARSEGKK(nx_trans).EQ.3) THEN
            ISC_GKK(1,nx_trans)=1
          ENDIF
        ENDIF
        IF(BIDOMAIN) THEN
          IF(SPARSEGKK(nx_ext).EQ.0) THEN
C$OMP PARALLEL DO
C$&     PRIVATE(nq),
C$&     SHARED(GKK,nx_ext)
            DO nq=1,NQT*NQT
              GKK(nq,nx_ext)=0.0d0
            ENDDO !nq
C$OMP END PARALLEL DO
          ELSE IF(SPARSEGKK(nx_ext).EQ.1) THEN
            ISR_GKK(1,nx_ext)=1
          ELSE IF(SPARSEGKK(nx_ext).EQ.3) THEN
            ISC_GKK(1,nx_ext)=1
          ENDIF
        ENDIF
      ENDIF

      IF(UPDATEDT) THEN
        IF(IMPLICIT) THEN
C         Main loop
          DO nrr=1,NRLIST(0)
            nr=NRLIST(nrr)

C$OMP PARALLEL DO
C$&     PRIVATE(COUNTINT,DTCMAMTHETA,IJ,IK,nii,nij,nik,NITB,nq,nq1,nq2,
C$&       nzero,PLACEINT,TEMPCOEFF,TEMPPLACE),
C$&     SHARED(ISC_GKK,ISR_GKK,GKK,NQGP,NQR,NQT,nr,NWQ,nx_trans,
C$&       SPARSEGKK,CQ,DT)
            DO nq=NQR(1,nr),NQR(2,nr) !create line for each nq
              DTCMAMTHETA=DT/(CQ(1,nq)*CQ(2,nq))*THETACOEF
              IF(SPARSEGKK(nx_trans).EQ.0) THEN !no sparsity
                IF(NWQ(1,nq).EQ.0) THEN !internal
                  DO nzero=1,NQGP(0,nq)
                    GKK(((NQGP(nzero,nq)-1)*NQT)+nq,nx_trans)=
     '                -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTCMAMTHETA
                    IF(NQGP(nzero,nq).EQ.nq) 
     '                GKK(((nq-1)*NQT)+nq,nx_trans)=
     '                GKK(((nq-1)*NQT)+nq,nx_trans)+1.0d0
                  ENDDO
                ENDIF
              ELSE IF(SPARSEGKK(nx_trans).EQ.1) THEN !compressed row
                IF(NWQ(1,nq).EQ.0) THEN !internal
                  PLACEINT=ISR_GKK(nq,nx_trans)
                  DO nzero=1,NQGP(0,nq)
                    GKK(PLACEINT,nx_trans)=
     '                -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTCMAMTHETA
                    IF(NQGP(nzero,nq).EQ.nq) GKK(PLACEINT,nx_trans)=
     '                GKK(PLACEINT,nx_trans)+1.0d0
                    PLACEINT=PLACEINT+1
                  ENDDO
                ENDIF
              ELSE IF(SPARSEGKK(nx_trans).EQ.2) THEN !row/col
                IF(NWQ(1,nq).EQ.0) THEN !internal
                  PLACEINT=ISC_GKK(nq,nx_trans)
                  DO nzero=1,NQGP(0,nq)
                    GKK(PLACEINT,nx_trans)=
     '                -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTCMAMTHETA
                    IF(NQGP(nzero,nq).EQ.nq) GKK(PLACEINT,nx_trans)=
     '                GKK(PLACEINT,nx_trans)+1.0d0
                    PLACEINT=PLACEINT+1
                  ENDDO !nzero
                ENDIF
              ELSE IF(SPARSEGKK(nx_trans).EQ.3) THEN !compressed column

C!!! THIS NEEDS REDOING FOR MP

                DTCMAMTHETA=NEWDT/(CQ(1,nq)*CQ(2,nq))*THETACOEF !10^-6 corr.

                COUNTINT=0
                NITB=NQXI(0,NQS(NEELEM(1,nr)))
                IK=MAX(0,NITB-2) !zero for 1,2-D, one for 3-D
                IJ=MIN(NITB-1,1) !zero for 1-D, one for 2,3-D
                
                DO nik=-IK,IK
                  DO nij=-IJ,IJ
                    DO nii=-1,1
                      nq1=NXQ(nii,1,NXQ(nij*2,1,NXQ(nik*3,1,nq)))
                      IF(nq1.GT.0) THEN
                        DO nzero=1,NQGP(0,nq1)
                          IF(NQGP(nzero,nq1).EQ.nq) THEN  
                            COUNTINT=COUNTINT+1
                            IF(NWQ(1,nq1).EQ.0) THEN !internal
                              TEMPCOEFF(COUNTINT)=
     '                          NQGW(NQGP_PIVOT(nzero,nq1),nq1)*
     '                          DTCMAMTHETA+1.0d0
                            ELSE !external
                              TEMPCOEFF(COUNTINT)=
     '                          NQGW(NQGP_PIVOT(nzero,nq1),nq1)
                            ENDIF
                            TEMPPLACE(COUNTINT)=nq1
                            nq2=NXQ(nii,1,NXQ(nij*2,1,NXQ(nik*3,1,nq1)))
                            IF((nq2.GT.0).AND.(nq2.NE.nq1)) THEN
                              DO nzero2=1,NQGP(0,nq2)
                                IF(NQGP(nzero2,nq2).EQ.nq) THEN  
                                  COUNTINT=COUNTINT+1
                                  IF(NWQ(1,nq2).EQ.0) THEN !internal
                                    TEMPCOEFF(COUNTINT)=
     '                                NQGW(NQGP_PIVOT(nzero2,nq2),nq2)*
     '                                DTCMAMTHETA
                                  ELSE !external
                                    TEMPCOEFF(COUNTINT)=
     '                                NQGW(NQGP_PIVOT(nzero2,nq2),nq2)
                                  ENDIF
                                  TEMPPLACE(COUNTINT)=nq2
                                ENDIF
                              ENDDO
                            ENDIF
                          ENDIF
                        ENDDO
                      ENDIF
                    ENDDO !nii
                  ENDDO !nij
                ENDDO !nik
                
                CALL ISORTP(COUNTINT,TEMPPLACE,PLACE_PIVOT)
                
                DO nzero=1,COUNTINT
                  PLACEINT=PLACEINT+1
                  GKK(PLACEINT,nx_trans)=TEMPCOEFF(PLACE_PIVOT(nzero))
                ENDDO
                
              ELSE IF(SPARSEGKK(nx_trans).EQ.5) THEN !row/col #2
                PLACEINT=ISC_GKK(nq,nx_trans)
                IF(NWQ(1,nq).EQ.0) THEN !internal
                  DO nzero=1,NQGP(0,nq)
                    GKK(PLACEINT,nx_trans)=
     '                -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTCMAMTHETA
                    IF(NQGP(nzero,nq).EQ.nq) GKK(PLACEINT,nx_trans)=
     '                GKK(PLACEINT,nx_trans)+1.0d0
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
C       Main loop
        DO nrr=1,NRLIST(0)
          nr=NRLIST(nrr)
          NITB=NQXI(0,NQS(NEELEM(1,nr)))
          IF(IMPLICIT) THEN
            DO nq=NQR(1,nr),NQR(2,nr) !create line for each nq
              DTCMAMTHETA=DT/(CQ(1,nq)*CQ(2,nq))*THETACOEF
              IF(SPARSEGKK(nx_trans).EQ.0) THEN !no sparsity
                CALL CALC_ADAMS_GRID_COEF(nq,NWQ(1,nq),NXQ(-NIM,0,nq),
     '            nx_ext,nx_trans,NITB,COEFFSEXT,GCHQ(1,nq),
     '            GUQ(1,1,nq),NQGW(1,nq),PROPQ(1,1,1,1,nq),.FALSE.,FIXQ,
     '            IMPLICIT,ERROR,*9999)
                IF(NWQ(1,nq).EQ.0) THEN !internal
                  DO nzero=1,NQGP(0,nq)
                    GKK(((NQGP(nzero,nq)-1)*NQT)+nq,nx_trans)=
     '                -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTCMAMTHETA
                    IF(NQGP(nzero,nq).EQ.nq) 
     '                GKK(((nq-1)*NQT)+nq,nx_trans)=
     '                GKK(((nq-1)*NQT)+nq,nx_trans)+1.0d0
                  ENDDO
                ELSE
                  DO nzero=1,NQGP(0,nq)
                    GKK(((NQGP(nzero,nq)-1)*NQT)+nq,nx_trans)=
     '                NQGW(NQGP_PIVOT(nzero,nq),nq)
                  ENDDO
                ENDIF
              ELSE IF(SPARSEGKK(nx_trans).EQ.1) THEN !compressed row
                CALL CALC_ADAMS_GRID_COEF(nq,NWQ(1,nq),NXQ(-NIM,0,nq),
     '            nx_ext,nx_trans,NITB,COEFFSEXT,GCHQ(1,nq),
     '            GUQ(1,1,nq),NQGW(1,nq),PROPQ(1,1,1,1,nq),.FALSE.,FIXQ,
     '            IMPLICIT,ERROR,*9999)
                COUNTINT=0
                IF(NWQ(1,nq).EQ.0) THEN !internal
                  DO nzero=1,NQGP(0,nq)
                    PLACEINT=PLACEINT+1
                    COUNTINT=COUNTINT+1
                    IF(.NOT.UPDATE) 
     '                ISC_GKK(PLACEINT,nx_trans)=NQGP(nzero,nq)
                    GKK(PLACEINT,nx_trans)=
     '                -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTCMAMTHETA
                    IF(NQGP(nzero,nq).EQ.nq) GKK(PLACEINT,nx_trans)=
     '                GKK(PLACEINT,nx_trans)+1.0d0
                  ENDDO
                ELSE !external
                  DO nzero=1,NQGP(0,nq)
                    PLACEINT=PLACEINT+1
                    COUNTINT=COUNTINT+1
                    IF(.NOT.UPDATE) 
     '                ISC_GKK(PLACEINT,nx_trans)=NQGP(nzero,nq)
                    GKK(PLACEINT,nx_trans)=
     '                NQGW(NQGP_PIVOT(nzero,nq),nq)
                  ENDDO
                ENDIF
                IF(.NOT.UPDATE) 
     '            ISR_GKK(nq+1,nx_trans)=ISR_GKK(nq,nx_trans)+COUNTINT
              ELSE IF(SPARSEGKK(nx_trans).EQ.2) THEN !row/col
                CALL CALC_ADAMS_GRID_COEF(nq,NWQ(1,nq),NXQ(-NIM,0,nq),
     '            nx_ext,nx_trans,NITB,COEFFSEXT,GCHQ(1,nq),
     '            GUQ(1,1,nq),NQGW(1,nq),PROPQ(1,1,1,1,nq),.FALSE.,FIXQ,
     '            IMPLICIT,ERROR,*9999)
                IF(NWQ(1,nq).EQ.0) THEN !internal
                  DO nzero=1,NQGP(0,nq)
                    PLACEINT=PLACEINT+1
                    IF(.NOT.UPDATE) THEN
                      ISR_GKK(PLACEINT,nx_trans)=nq
                      ISC_GKK(PLACEINT,nx_trans)=NQGP(nzero,nq)
                    ENDIF
                    GKK(PLACEINT,nx_trans)=
     '                -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTCMAMTHETA
                    IF(NQGP(nzero,nq).EQ.nq) GKK(PLACEINT,nx_trans)=
     '                GKK(PLACEINT,nx_trans)+1.0d0
                  ENDDO !nzero
                ELSE !external
                  DO nzero=1,NQGP(0,nq)
                    PLACEINT=PLACEINT+1
                    IF(.NOT.UPDATE) THEN
                      ISR_GKK(PLACEINT,nx_trans)=nq
                      ISC_GKK(PLACEINT,nx_trans)=NQGP(nzero,nq)
                    ENDIF
                    GKK(PLACEINT,nx_trans)=NQGW(NQGP_PIVOT(nzero,nq),nq)
                  ENDDO !nzero
                ENDIF
              ELSE IF(SPARSEGKK(nx_trans).EQ.3) THEN !compressed column
                COUNTINT=0
                IK=MAX(0,NITB-2) !zero for 1,2-D, one for 3-D
                IJ=MIN(NITB-1,1) !zero for 1-D, one for 2,3-D
                
                DO nik=-IK,IK
                  DO nij=-IJ,IJ
                    DO nii=-1,1
                      nq1=NXQ(nii,1,NXQ(nij*2,1,NXQ(nik*3,1,nq)))
                      IF(nq1.GT.0) THEN
                        DO nzero=1,NQGP(0,nq1)
                          IF(NQGP(nzero,nq1).EQ.nq) THEN  
                            CALL CALC_ADAMS_GRID_COEF(nq1,NWQ(1,nq1),
     '                        NXQ(-NIM,0,nq1),nx_ext,nx_trans,NITB,
     '                        COEFFSEXT,GCHQ(1,nq1),GUQ(1,1,nq1),
     '                        NQGW(1,nq1),PROPQ(1,1,1,1,nq1),.FALSE.,
     '                        FIXQ,IMPLICIT,ERROR,*9999)
                            COUNTINT=COUNTINT+1
                            IF(NWQ(1,nq1).EQ.0) THEN !internal
                              TEMPCOEFF(COUNTINT)=
     '                          NQGW(NQGP_PIVOT(nzero,nq1),nq1)*
     '                          DTCMAMTHETA+1.0d0
                            ELSE !external
                              TEMPCOEFF(COUNTINT)=
     '                          NQGW(NQGP_PIVOT(nzero,nq1),nq1)
                            ENDIF
                            TEMPPLACE(COUNTINT)=nq1
                            nq2=NXQ(nii,1,NXQ(nij*2,1,NXQ(nik*3,1,nq1)))
                            IF((nq2.GT.0).AND.(nq2.NE.nq1)) THEN
                              DO nzero2=1,NQGP(0,nq2)
                                IF(NQGP(nzero2,nq2).EQ.nq) THEN  
                                  CALL CALC_ADAMS_GRID_COEF(nq2,
     '                              NWQ(1,nq2),NXQ(-NIM,0,nq2),nx_ext,
     '                              nx_trans,NITB,COEFFSEXT,GCHQ(1,nq2),
     '                              GUQ(1,1,nq2),NQGW(1,nq2),
     '                              PROPQ(1,1,1,1,nq2),.FALSE.,FIXQ,
     '                              IMPLICIT,ERROR,*9999)
                                  COUNTINT=COUNTINT+1
                                  IF(NWQ(1,nq2).EQ.0) THEN !internal
                                    TEMPCOEFF(COUNTINT)=
     '                                -NQGW(NQGP_PIVOT(nzero2,nq2),nq2)*
     '                                DTCMAMTHETA
                                  ELSE !external
                                    TEMPCOEFF(COUNTINT)=
     '                                NQGW(NQGP_PIVOT(nzero2,nq2),nq2)
                                  ENDIF
                                  TEMPPLACE(COUNTINT)=nq2
                                ENDIF
                              ENDDO
                            ENDIF
                          ENDIF
                        ENDDO
                      ENDIF
                    ENDDO !nii
                  ENDDO !nij
                ENDDO !nik
                
                CALL ISORTP(COUNTINT,TEMPPLACE,PLACE_PIVOT)
                
                DO nzero=1,COUNTINT
                  PLACEINT=PLACEINT+1
                  IF(.NOT.UPDATE) 
     '              ISR_GKK(PLACEINT,nx_trans)=TEMPPLACE(nzero)
                  GKK(PLACEINT,nx_trans)=TEMPCOEFF(PLACE_PIVOT(nzero))
                ENDDO
                IF(.NOT.UPDATE) 
     '            ISC_GKK(nq+1,nx_trans)=ISC_GKK(nq,nx_trans)+COUNTINT
              ELSE IF(SPARSEGKK(nx_trans).EQ.5) THEN !row/col #2
                CALL CALC_ADAMS_GRID_COEF(nq,NWQ(1,nq),NXQ(-NIM,0,nq),
     '            nx_ext,nx_trans,NITB,COEFFSEXT,GCHQ(1,nq),
     '            GUQ(1,1,nq),NQGW(1,nq),PROPQ(1,1,1,1,nq),.FALSE.,
     '            FIXQ,IMPLICIT,ERROR,*9999)
                IF(NWQ(1,nq).EQ.0) THEN !internal
                  DO nzero=1,NQGP(0,nq)
                    PLACEINT=PLACEINT+1
                    IF(.NOT.UPDATE) ISC_GKK(PLACEINT,nx_trans)=nq
                    GKK(PLACEINT,nx_trans)=
     '                -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTCMAMTHETA
                    IF(NQGP(nzero,nq).EQ.nq) GKK(PLACEINT,nx_trans)=
     '                GKK(PLACEINT,nx_trans)+1.0d0
                  ENDDO
                ELSE !external
                  DO nzero=1,NQGP(0,nq)
                    PLACEINT=PLACEINT+1
                    IF(.NOT.UPDATE) ISC_GKK(PLACEINT,nx_trans)=nq
                    GKK(PLACEINT,nx_trans)=NQGW(NQGP_PIVOT(nzero,nq),nq)
                  ENDDO
                ENDIF
              ENDIF
            ENDDO !nq
          ENDIF
          IF(BIDOMAIN) THEN
            DO nq=NQR(1,nr),NQR(2,nr) !create line for each nq
              IF(SPARSEGKK(nx_ext).EQ.0) THEN !no sparsity
                CALL CALC_ADAMS_GRID_COEF(nq,NWQ(1,nq),NXQ(-NIM,0,nq),
     '            nx_ext,nx_trans,NITB,COEFFSEXT,GCHQ(1,nq),
     '            GUQ(1,1,nq),NQGW(1,nq),PROPQ(1,1,1,1,nq),BIDOMAIN,
     '            FIXQ,NOTIMP,ERROR,*9999)
                DO nzero=1,NQGP(0,nq)
                  GKK(((NQGP(nzero,nq)-1)*NQT)+nq,nx_ext)=
     '              COEFFSEXT(NQGP_PIVOT(nzero,nq))
                ENDDO
              ELSE IF(SPARSEGKK(nx_ext).EQ.1) THEN !compressed row
                CALL CALC_ADAMS_GRID_COEF(nq,NWQ(1,nq),NXQ(-NIM,0,nq),
     '            nx_ext,nx_trans,NITB,COEFFSEXT,GCHQ(1,nq),
     '            GUQ(1,1,nq),NQGW(1,nq),PROPQ(1,1,1,1,nq),BIDOMAIN,
     '            FIXQ,NOTIMP,ERROR,*9999)
                COUNTEXT=0
                DO nzero=1,NQGP(0,nq)
                  PLACEEXT=PLACEEXT+1
                  COUNTEXT=COUNTEXT+1
                  IF(.NOT.UPDATE) 
     '              ISC_GKK(PLACEEXT,nx_ext)=NQGP(nzero,nq)
                  GKK(PLACEEXT,nx_ext)=COEFFSEXT(NQGP_PIVOT(nzero,nq))
                ENDDO
                IF(.NOT.UPDATE) 
     '            ISR_GKK(nq+1,nx_ext)=ISR_GKK(nq,nx_ext)+COUNTEXT
              ELSE IF(SPARSEGKK(nx_ext).EQ.2) THEN !row/col
                CALL CALC_ADAMS_GRID_COEF(nq,NWQ(1,nq),NXQ(-NIM,0,nq),
     '            nx_ext,nx_trans,NITB,COEFFSEXT,GCHQ(1,nq),
     '            GUQ(1,1,nq),NQGW(1,nq),PROPQ(1,1,1,1,nq),BIDOMAIN,
     '            FIXQ,NOTIMP,ERROR,*9999)
                DO nzero=1,NQGP(0,nq)
                  PLACEEXT=PLACEEXT+1
                  IF(.NOT.UPDATE) THEN
                    ISR_GKK(PLACEEXT,nx_ext)=nq
                    ISC_GKK(PLACEEXT,nx_ext)=NQGP(nzero,nq)
                  ENDIF
                  GKK(PLACEEXT,nx_ext)=COEFFSEXT(NQGP_PIVOT(nzero,nq))
                ENDDO
              ELSEIF(SPARSEGKK(nx_ext).EQ.3) THEN !compressed column
                COUNTEXT=0
                IK=MAX(0,NITB-2) !zero for 1,2-D, one for 3-D
                IJ=MIN(NITB-1,1) !zero for 1-D, one for 2,3-D
                
                DO nik=-IK,IK
                  DO nij=-IJ,IJ
                    DO nii=-1,1
                      nq1=NXQ(nii,1,NXQ(nij*2,1,NXQ(nik*3,1,nq)))
                      IF(nq1.GT.0) THEN
                        DO nzero=1,NQGP(0,nq1)
                          IF(NQGP(nzero,nq1).EQ.nq) THEN  
                            CALL CALC_ADAMS_GRID_COEF(nq1,NWQ(1,nq1),
     '                        NXQ(-NIM,0,nq1),nx_ext,nx_trans,NITB,
     '                        COEFFSEXT,GCHQ(1,nq1),GUQ(1,1,nq1),
     '                        NQGW(1,nq1),PROPQ(1,1,1,1,nq1),BIDOMAIN,
     '                        FIXQ,NOTIMP,ERROR,*9999)
                            COUNTEXT=COUNTEXT+1
                            TEMPCOEFF(COUNTEXT)=
     '                        COEFFSEXT(NQGP_PIVOT(nzero,nq1))
                            TEMPPLACE(COUNTEXT)=nq1
                            nq2=NXQ(nii,1,NXQ(nij*2,1,NXQ(nik*3,1,nq1)))
                            IF((nq2.GT.0).AND.(nq2.NE.nq1)) THEN
                              DO nzero2=1,NQGP(0,nq2)
                                IF(NQGP(nzero2,nq2).EQ.nq) THEN  
                                  CALL CALC_ADAMS_GRID_COEF(nq2,
     '                              NWQ(1,nq2),NXQ(-NIM,0,nq2),nx_ext,
     '                              nx_trans,NITB,COEFFSEXT,GCHQ(1,nq2),
     '                              GUQ(1,1,nq2),NQGW(1,nq2),
     '                              PROPQ(1,1,1,1,nq2),BIDOMAIN,FIXQ,
     '                              NOTIMP,ERROR,*9999)
                                  COUNTEXT=COUNTEXT+1
                                  TEMPCOEFF(COUNTEXT)=
     '                              COEFFSEXT(NQGP_PIVOT(nzero2,nq2))
                                  TEMPPLACE(COUNTEXT)=nq2
                                ENDIF
                              ENDDO
                            ENDIF
                          ENDIF
                        ENDDO
                      ENDIF
                    ENDDO !nii
                  ENDDO !nij
                ENDDO !nik
                
                CALL ISORTP(COUNTEXT,TEMPPLACE,PLACE_PIVOT)
                
                DO nzero=1,COUNTEXT
                  PLACEEXT=PLACEEXT+1
                  IF(.NOT.UPDATE) 
     '              ISR_GKK(PLACEEXT,nx_ext)=TEMPPLACE(nzero)
                  GKK(PLACEEXT,nx_ext)=TEMPCOEFF(PLACE_PIVOT(nzero))
                ENDDO !nzero
                IF(.NOT.UPDATE) 
     '            ISC_GKK(nq+1,nx_ext)=ISC_GKK(nq,nx_ext)+COUNTEXT
              ELSE IF(SPARSEGKK(nx_ext).EQ.5) THEN !row/col #2
                CALL CALC_ADAMS_GRID_COEF(nq,NWQ(1,nq),NXQ(-NIM,0,nq),
     '            nx_ext,nx_trans,NITB,COEFFSEXT,GCHQ(1,nq),
     '            GUQ(1,1,nq),NQGW(1,nq),PROPQ(1,1,1,1,nq),
     '            BIDOMAIN,FIXQ,NOTIMP,ERROR,*9999)
                DO nzero=1,NQGP(0,nq)
                  PLACEEXT=PLACEEXT+1
                  IF(.NOT.UPDATE) ISC_GKK(PLACEEXT,nx_ext)=nq
                  GKK(PLACEEXT,nx_ext)=COEFFSEXT(NQGP_PIVOT(nzero,nq))
                ENDDO !nzero
              ENDIF
            ENDDO !nq
          ENDIF
        ENDDO !nr
        
        IF(.NOT.UPDATE) THEN
          IF(IMPLICIT) THEN
            IF(SPARSEGKK(nx_trans).EQ.5) THEN
              COUNTINT=0
              DO nrr=1,NRLIST(0)
                nr=NRLIST(nrr)
                DO nq=NQR(1,nr),NQR(2,nr)
                  DO nzero=1,NQGP(0,nq)                  
                    COUNTINT=COUNTINT+1
                    ISC_GKK(PLACEINT+COUNTINT,nx_trans)=NQGP(nzero,nq)
                  ENDDO !nzero
                ENDDO !nq
              ENDDO !nrr
            ENDIF
          ENDIF
          IF(BIDOMAIN) THEN
            IF(SPARSEGKK(nx_ext).EQ.5) THEN
              COUNTEXT=0
              DO nrr=1,NRLIST(0)
                nr=NRLIST(nrr)
                DO nq=NQR(1,nr),NQR(2,nr)
                  DO nzero=1,NQGP(0,nq)
                    COUNTEXT=COUNTEXT+1
                    ISC_GKK(PLACEEXT+COUNTEXT,nx_ext)=NQGP(nzero,nq)
                  ENDDO !nzero
                ENDDO !nq
              ENDDO !nrr
            ENDIF
          ENDIF  
        ENDIF
      ENDIF !updatedt
    
      UP_GRID_MATERIAL=.FALSE.
      UP_GRID_TENSOR=.FALSE.      
      FIRST_A=.TRUE.
      UPDATE_MATRIX=.TRUE.

      IF(.NOT.UPDATE) THEN
        IF(nx_upd.GT.0) THEN
          SPARSEGKK(nx_upd)=SPARSEGKK(nx_ext)
          NZZT(1,NRLIST(1),nx_upd)=NZZT(1,NRLIST(1),nx_ext)
          IF(SPARSEGKK(nx_upd).EQ.0) THEN
C$OMP PARALLEL DO
C$&     PRIVATE(nq),
C$&     SHARED(GKK,NRLIST,nx_ext,nx_upd)
            DO nq=1,NZZT(1,NRLIST(1),nx_ext)
              GKK(nq,nx_upd)=GKK(nq,nx_ext)
            ENDDO
C$OMP END PARALLEL DO
          ELSEIF(SPARSEGKK(nx_upd).EQ.1) THEN
C$OMP PARALLEL DO
C$&     PRIVATE(nq),
C$&     SHARED(ISC_GKK,GKK,NRLIST,nx_ext,nx_upd)
            DO nq=1,NZZT(1,NRLIST(1),nx_ext)
              GKK(nq,nx_upd)=GKK(nq,nx_ext)
              ISC_GKK(nq,nx_upd)=ISC_GKK(nq,nx_ext)
            ENDDO
C$OMP END PARALLEL DO
            DO nq=1,NQT+1
              ISR_GKK(nq,nx_upd)=ISR_GKK(nq,nx_ext)
            ENDDO
          ELSEIF((SPARSEGKK(nx_upd).EQ.2).OR.
     '      (SPARSEGKK(nx_upd).EQ.4)) THEN
C$OMP PARALLEL DO
C$&     PRIVATE(nq),
C$&     SHARED(ISC_GKK,ISR_GKK,GKK,NRLIST,nx_ext,nx_upd)
            DO nq=1,NZZT(1,NRLIST(1),nx_ext)
              GKK(nq,nx_upd)=GKK(nq,nx_ext)
              ISR_GKK(nq,nx_upd)=ISR_GKK(nq,nx_ext)
              ISC_GKK(nq,nx_upd)=ISC_GKK(nq,nx_ext)
            ENDDO
C$OMP END PARALLEL DO
          ELSEIF(SPARSEGKK(nx_upd).EQ.3) THEN
C$OMP PARALLEL DO
C$&     PRIVATE(nq),
C$&     SHARED(ISR_GKK,GKK,NRLIST,nx_ext,nx_upd)
            DO nq=1,NZZT(1,NRLIST(1),nx_ext)
              GKK(nq,nx_upd)=GKK(nq,nx_ext)
              ISR_GKK(nq,nx_upd)=ISR_GKK(nq,nx_ext)
            ENDDO
C$OMP END PARALLEL DO
            DO nq=1,NQT+1
              ISC_GKK(nq,nx_upd)=ISC_GKK(nq,nx_ext)
            ENDDO
          ELSEIF(SPARSEGKK(nx_upd).EQ.5) THEN
C$OMP PARALLEL DO
C$&     PRIVATE(nq),
C$&     SHARED(ISC_GKK,GKK,NRLIST,nx_ext,nx_upd)
            DO nq=1,NZZT(1,NRLIST(1),nx_ext)
              GKK(nq,nx_upd)=GKK(nq,nx_ext)
              ISC_GKK(nq,nx_upd)=ISC_GKK(nq,nx_ext)
              ISC_GKK(nq+NZZT(1,NRLIST(1),nx_ext),nx_upd)=
     '          ISC_GKK(nq+NZZT(1,NRLIST(1),nx_ext),nx_ext)
            ENDDO
C$OMP END PARALLEL DO
          ENDIF
        ENDIF
      ENDIF

      IF(DOP) THEN
        nr=NRLIST(1)
        WRITE(OP_STRING,
     '    '(/'' Global stiffness matrix GKK - transmembrane:'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' NOT(1,1,nr,nx)='',I5,'
     '    //''', NOT(2,1,nr,nx)='',I5)') NOT(1,1,nr,nx_trans),
     '    NOT(2,1,nr,nx_trans)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        CALL OPSTFMAT(DUMMY_LIST,ISC_GKK,ISR_GKK,IOOP,
     '    NOT(1,1,nr,nx_trans),NOT(2,1,nr,nx_trans),
     '    NZZT(1,nr,nx_trans),DUMMY_LIST,SPARSEGKK(nx_trans),
     '    GKK,GKK,'GKK','GKK',.TRUE.,.TRUE.,.FALSE.,
     '    ERROR,*9999)
        IF(BIDOMAIN) THEN
          WRITE(OP_STRING,
     '      '(/'' Global stiffness matrix GKK - external:'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NOT(1,1,nr,nx)='',I5,'
     '      //''', NOT(2,1,nr,nx)='',I5)') NOT(1,1,nr,nx_ext),
     '      NOT(2,1,nr,nx_ext)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CALL OPSTFMAT(DUMMY_LIST,ISC_GKK,ISR_GKK,IOOP,
     '      NOT(1,1,nr,nx_ext),NOT(2,1,nr,nx_ext),
     '      NZZT(1,nr,nx_ext),DUMMY_LIST,SPARSEGKK(nx_ext),
     '      GKK,GKK,'GKK','GKK',.TRUE.,.TRUE.,.FALSE.,
     '      ERROR,*9999)
        ENDIF
      ENDIF

      IF(IWRIT5(NRLIST(1),nx_trans).GE.1) THEN
        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)
        WRITE(OP_STRING,'(1X,''Time for matrix assembly '',F8.2,'
     '    //'''s cpu'')') ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('ADAMS_ASSEMBLE10')
      RETURN
 9999 CALL ERRORS('ADAMS_ASSEMBLE10',ERROR)
      CALL EXITS('ADAMS_ASSEMBLE10')
      RETURN 1
      END


      SUBROUTINE ADAMS_MARCH8(ADAMS_IWORK,IBT,IDO,INP,ISC_GKK,ISR_GKK,
     '  NBH,NBJ,NEELEM,NENQ,NHE,NHP,NKE,NKH,NPF,NP_INTERFACE,NPNE,
     '  NPNODE,NQGP,NQGP_PIVOT,NQGW,NQNE,NQS,NQXI,NRLIST,NVHE,NVHP,
     '  NW,NWQ,NXLIST,NXQ,NYNE,NYNP,ADAMS_WORK,AQ,CQ,CURVCORRECT,
     '  GCHQ,GKK,GUQ,PROPQ,SE,XQ,YP,YQ,ZA,ZE,ZP,FIXQ,ITER8,ERROR,*)

C#### Subroutine: ADAMS_MARCH8
C###  Description:
C###    ADAMS_MARCH8 is the adaptive time stepping routine for the
C###    general grid ionic current activation problem when using the
C###    adaptive Adams-Moulton integrator.

      IMPLICIT NONE

      INCLUDE 'cmiss$reference:adam00.cmn'
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b12.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:grid00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:iwrit00.cmn'
      INCLUDE 'cmiss$reference:ktyp30.cmn'
      INCLUDE 'cmiss$reference:loc00.inc'
      INCLUDE 'cmiss$reference:maqloc00.inc'
      INCLUDE 'cmiss$reference:marc00.cmn'
      INCLUDE 'cmiss$reference:nqloc00.inc'
      INCLUDE 'cmiss$reference:solv00.cmn'
      INCLUDE 'cmiss$reference:time02.cmn'
      INCLUDE 'cmiss$reference:tol00.cmn'

!     Parameter list
      INTEGER ADAMS_IWORK(ADAMS_LIWORK,NQM,NXM),IBT(3,NIM,NBFM),
     '  IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISC_GKK(NISC_GKKM,NXM),ISR_GKK(NISR_GKKM,NXM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NENQ(0:8,NQM),NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),
     '  NKE(NKM,NNM,NBFM,NEM),NKH(NHM,NPM,NCM,0:NRM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NP_INTERFACE(0:NPM,0:3),
     '  NPNODE(0:NP_R_M,0:NRM),NQGP(0:19,NQM),NQGP_PIVOT(19,NQM),
     '  NQNE(NEM,NQEM),NQS(NEM),NQXI(0:NIM,NQSCM),NRLIST(0:NRM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),NW(NEM,3),NWQ(6,0:NQM),
     '  NXLIST(0:NXM),NXQ(-NIM:NIM,0:4,0:NQM)
      REAL*8 ADAMS_WORK(ADAMS_LWORK,NQM,NXM),AQ(NMAQM,NQM),
     '  CQ(NMM,NQM),CURVCORRECT(2,2,NNM,NEM),
     '  GCHQ(3,NQM),GKK(NZ_GKK_M,NXM),GUQ(3,3,NQM),NQGW(19,NQM),
     '  PROPQ(3,3,4,2,NQM),SE(NSM,NBFM,NEM),XQ(NJM,NQM),
     '  YP(NYM,NIYM,NXM),YQ(NYQM,NIQM,NAM,NXM),ZA(NAM,NHM,NCM,NEM),
     '  ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
      LOGICAL FIXQ(NYQM,NIYFIXM,NXM),ITER8
!     Local variables
      INTEGER maqdt,maqp1t0,maqp1t1,maqp1i,maqp2t0,maqp2t1,maqp2i,
     '  niqBNDRY,niqDV,niq_old,niq_old2,niqV,NITB,nnq,NNQ_MIN,
     '  nonrlist,nq,nr,NRTEMP(0:1),NSOL,NSTEP,nxc,nx_ext,nx,
     '  nx_torso,nx_upd,YLIST(0:99)
      REAL*8 COEFFSEXT(19),DIFFUSION,OLDDT,PHI1,PHI2,
     '  PSTIMULUS,RHS(NQM),T
      REAL ELAPSED_TIME,TIME_START(1),TIME_START2(1),TIME_STOP(1),
     '  TIME_STOP2(1)
      LOGICAL BIDOMAIN,CHMTRIX,CONTINU,FIRSTADAMSDTAR,FIRSTITER,FIRST_A,
     '  IMPLICIT,UPDATE_MATRIX,UP_GRID_DELTAT,X_INIT,ERROR_FLAG

      INTEGER L_AII,L_AIO,L_ARI,L_ARO,L_CONTROL,L_MODEL,L_PARAMETERS
      PARAMETER(L_AII=1,L_AIO=1,L_ARI=1,L_ARO=1,L_CONTROL=1,
     '  L_MODEL=2,L_PARAMETERS=99)
      INTEGER AII(L_AII),AIO(L_AIO),CONTROL(L_CONTROL),ERR_CODE,
     '  MODEL(L_MODEL)
      REAL*8 ARI(L_ARI),ARO(L_ARO),PARAMETERS(L_PARAMETERS)

      SAVE FIRSTITER,NNQ_MIN,NSTEP,TIME_START
      
      SAVE MODEL,PARAMETERS,OLDDT

      CALL ENTERS('ADAMS_MARCH8',*9999)
     
      CALL ASSERT(KTYP37.EQ.5,'>>Must have KTYP37=5 in ADAMS_MARCH8',
     '  ERROR,*9999)
      
      DT=TINCR
      X_INIT=.FALSE.
      CHMTRIX=.FALSE.

C     Check for Bidomain solution
      IF(KTYP32.EQ.1) THEN !monodomain
        BIDOMAIN=.FALSE.
      ELSE IF(KTYP32.EQ.2) THEN !bidomain
        BIDOMAIN=.TRUE.

      ELSE
        CALL ASSERT(.FALSE.,'>>Invalid type in KTYP32 (mono/bidomain)',
     '    ERROR,*9999)
      ENDIF     

C     Get indices for AQ auxiliary parameters
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE1,maqp1t0,MAQ_START,
     '  ERROR,*9999)
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE1,maqp1t1,MAQ_STOP,
     '  ERROR,*9999)
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE1,maqp1i,MAQ_CURRENT,
     '  ERROR,*9999)
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE2,maqp2t0,MAQ_START,
     '  ERROR,*9999)
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE2,maqp2t1,MAQ_STOP,
     '  ERROR,*9999)
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE2,maqp2i,MAQ_CURRENT,
     '  ERROR,*9999)
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_TIME,maqdt,MAQ_DT,
     '  ERROR,*9999)

! Get class information
      nxc=NXLIST(1)
      CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
      CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '  ERROR,*9999)
      IF(BIDOMAIN) THEN
        CALL ASSERT(NXLIST(0).GE.2,'>>Error - bidomain needs 2 classes',
     '    ERROR,*9999)
        nxc=NXLIST(2)
        CALL NX_LOC(NX_INQUIRE,nxc,nx_ext,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx_ext.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)
      ELSE
        nx_ext=0
      ENDIF
      IF(NXLIST(0).GE.3) THEN
        nxc=NXLIST(3)
        CALL NX_LOC(NX_INQUIRE,nxc,nx_upd,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx_upd.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)
      ELSE
        nx_upd=0
      ENDIF
      IF(NXLIST(0).GE.4) THEN
        nxc=NXLIST(4)
        CALL NX_LOC(NX_INQUIRE,nxc,nx_torso,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx_torso.GT.0,
     '    '>>No nx defined for this solve class',ERROR,*9999)
      ELSE
        nx_torso=0
      ENDIF

C     Check for implicit/explicit solution method
      IF(THETA(1).GT.ZERO_TOL) THEN
        IMPLICIT=.TRUE.
      ELSE
        IMPLICIT=.FALSE.
      ENDIF

C     Get YQ-Y mappings
      CALL GETYQYMAP(NRLIST(1),nx,YLIST,ERROR,*9999)

      IF(RESTART) THEN
        T=T_RESTART
        FIRSTADAMSDTAR=.FALSE.
      ELSE
        IF(.NOT.ITER8) THEN
          T=TSTART
          NSTEP=0
          FIRSTITER=.TRUE.
          FIRSTADAMSDTAR=.TRUE.
          NNQ_MIN=1
          IF(ITYP3(NRLIST(1),nx).EQ.5.AND.KTYP33.EQ.2) THEN
            MODEL(1)=2
            MODEL(2)=2
            CALL JRWP_INIT_PARAM(L_PARAMETERS,PARAMETERS,ERR_CODE)
            IF(ERR_CODE.NE.0) THEN
              ERROR='>>L_PARAMETERS too small'
              GOTO 9999
            ENDIF
          ENDIF

          DO nq=1,NQT
            ADAMS_IWORK(1,nq,nx)=0
            ADAMS_WORK(1,nq,nx)=DT
          ENDDO

          IF(KTYP36.EQ.1) THEN !using DTAR
C         Initialise dynamic tracking
C$OMP PARALLEL DO
C$&     PRIVATE(nq),
C$&     SHARED(NWQ)
            DO nq=1,NQT
              NWQ(5,nq)=0
              NWQ(4,nq)=0
            ENDDO
C$OMP END PARALLEL DO
            NWQ(5,0)=0
          ELSE
            NWQ(5,0)=0
            DO nonrlist=1,NRLIST(0)
              nr=NRLIST(nonrlist)
              DO nq=NQR(1,nr),NQR(2,nr)
                IF(NWQ(1,nq).EQ.0) THEN !internal
                  NWQ(5,0)=NWQ(5,0)+1
                  NWQ(5,NWQ(5,0))=nq
                ENDIF
              ENDDO 
            ENDDO
          ENDIF

        ENDIF
      ENDIF

      OLDDT=DT
      UP_GRID_DELTAT=DT.NE.OLDDT

C     Create/update matricies if implicit
      IF(IMPLICIT.OR.BIDOMAIN) THEN
        IF(RESTART) THEN
          UPDATE_MATRIX=.FALSE.
          FIRST_A=.FALSE.
          IF(UP_GRID_MATERIAL.OR.UP_GRID_TENSOR.OR.UP_GRID_DELTAT) THEN
            CALL ADAMS_ASSEMBLE10(ISC_GKK,ISR_GKK,NEELEM,NQGP,
     '        NQGP_PIVOT,NQS,NQXI,NRLIST,NWQ,nx_ext,nx,
     '        nx_upd,NXQ,CQ,GCHQ,GUQ,GKK,DT,NQGW,OLDDT,
     '        PROPQ,BIDOMAIN,FIRST_A,FIXQ,IMPLICIT,.TRUE.,
     '        UP_GRID_DELTAT,UPDATE_MATRIX,ERROR,*9999)
            CHMTRIX=.TRUE.
          ENDIF
        ELSE
          IF(.NOT.ITER8) THEN
            CALL ADAMS_ASSEMBLE10(ISC_GKK,ISR_GKK,NEELEM,NQGP,
     '        NQGP_PIVOT,NQS,NQXI,NRLIST,NWQ,nx_ext,nx,
     '        nx_upd,NXQ,CQ,GCHQ,GUQ,GKK,DT,NQGW,OLDDT,PROPQ,
     '        BIDOMAIN,FIRST_A,FIXQ,IMPLICIT,.FALSE.,
     '        UP_GRID_DELTAT,UPDATE_MATRIX,ERROR,*9999)
            CHMTRIX=.TRUE.
          ENDIF
        ENDIF
      ELSE
        IF(.NOT.RESTART) THEN
C         Calculuate NQGP and NQGW
          CALL GET_FD_POINTS(NEELEM,NQGP,NQGP_PIVOT,NQS,NQXI,NRLIST,
     '      NWQ,NXQ,ERROR,*9999)
          DO nonrlist=1,NRLIST(0)
            nr=NRLIST(nonrlist)
            NITB=NQXI(0,NQS(NEELEM(1,nr)))
            ERROR_FLAG=.FALSE.
C$OMP PARALLEL DO
C$&     PRIVATE(nq,COEFFSEXT),
C$&     SHARED(ERROR_FLAG,GCHQ,GUQ,NITB,NQGW,NQR,NWQ,NXQ,nx_ext,nx,
C$&            PROPQ,FIXQ,IMPLICIT)
            DO nq=NQR(1,nr),NQR(2,nr)
              IF(.NOT.ERROR_FLAG) THEN
                CALL CALC_ADAMS_GRID_COEF(nq,NWQ(1,nq),NXQ(-NIM,0,nq),
     '            nx_ext,nx,NITB,COEFFSEXT,GCHQ(1,nq),GUQ(1,1,nq),
     '            NQGW(1,nq),PROPQ(1,1,1,1,nq),.FALSE.,FIXQ,IMPLICIT,
     '            ERROR,*100)
                GOTO 102
 100            CONTINUE
C$OMP CRITICAL(ADAMS_MARCH8_1)
                ERROR_FLAG=.TRUE.
                WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
                CALL WRITES(IOER,OP_STRING,ERROR,*101)
                WRITE(OP_STRING,'(/'' >>An error occurred - '
     '            //'results may be unreliable!'')')
                CALL WRITES(IODI,OP_STRING,ERROR,*101)
 101            CONTINUE
C$OMP END CRITICAL(ADAMS_MARCH8_1)
 102            CONTINUE
              ENDIF !.NOT.ERROR_FLAG  
            ENDDO !nq
C$OMP END PARALLEL DO
          ENDDO !nr
        ENDIF
      ENDIF

      CONTINU=.TRUE.
      CALL CPU_TIMER(CPU_USER,TIME_START)
      CALL CPU_TIMER(CPU_USER,TIME_START2)
      CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niq_old,NIQ_OLDSOLN,ERROR,*9999)
      CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niq_old2,NIQ_OLDSOLN2,
     '  ERROR,*9999)
      CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqV,NIQ_V,ERROR,*9999)
      CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqBNDRY,NIQ_BNDRY,ERROR,*9999)
      IF((ITYP3(NRLIST(1),nx).EQ.3.AND.KTYP33.EQ.2)) THEN !VCDC
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqDV,NIQ_DV,ERROR,*9999)
      ENDIF

      IF(.NOT.ITER8) THEN
C       Main time integration loop
        DO WHILE(CONTINU)
          NSTEP=NSTEP+1
          NSOL=0
                    
          UP_GRID_DELTAT=DT.NE.OLDDT

          IF((IMPLICIT.OR.BIDOMAIN).AND.UP_GRID_DELTAT) THEN
            CALL ADAMS_ASSEMBLE10(ISC_GKK,ISR_GKK,NEELEM,NQGP,
     '        NQGP_PIVOT,NQS,NQXI,NRLIST,NWQ,nx_ext,nx,nx_upd,
     '        NXQ,CQ,GCHQ,GUQ,GKK,DT,NQGW,OLDDT,PROPQ,BIDOMAIN,
     '        FIRST_A,FIXQ,IMPLICIT,.TRUE.,UP_GRID_DELTAT,
     '        UPDATE_MATRIX,ERROR,*9999)
            UPDATE_MATRIX=.TRUE.
            CHMTRIX=.TRUE.
          ENDIF
          
          IF(KTYP36.EQ.1) THEN
            CALL CALC_ADAMS_DTAR(maqp1t0,maqp1t1,maqp1i,maqp2t0,
     '        maqp2t1,maqp2i,niq_old,niqV,nnq_min,NQXI,NRLIST,
     '        NSOL,NWQ,nx,NXQ,ADAMS_WORK(1,1,nx),AQ,CQ,T,
     '        YQ(1,1,1,nx),FIRSTADAMSDTAR,ERROR,*9999)
          ENDIF
          
C         Store solution as basis for next time step
          IF(BIDOMAIN) THEN
            DO nonrlist=1,NRLIST(0)
              nr=NRLIST(nonrlist)
C$OMP PARALLEL DO
C$&     PRIVATE(nq),
C$&     SHARED(niq_old,niq_old2,niqV,NQR,nr,nx,nx_ext,YQ)
              DO nq=NQR(1,nr),NQR(2,nr)
C              YQ(nq,niq_old2,1,nx)=YQ(nq,niq_old,1,nx)
C              YQ(nq,niq_old2,1,nx_ext)=YQ(nq,niq_old,1,nx_ext)
                YQ(nq,niq_old,1,nx)=YQ(nq,niqV,1,nx)
                YQ(nq,niq_old,1,nx_ext)=YQ(nq,niqV,1,nx_ext)
              ENDDO
C$OMP END PARALLEL DO
            ENDDO
          ELSE
            DO nonrlist=1,NRLIST(0)
              nr=NRLIST(nonrlist)
C$OMP PARALLEL DO
C$&     PRIVATE(nq),
C$&     SHARED(niq_old,niq_old2,niqV,NQR,nr,nx,YQ)
              DO nq=NQR(1,nr),NQR(2,nr)
C              YQ(nq,niq_old2,1,nx)=YQ(nq,niq_old,1,nx)
                YQ(nq,niq_old,1,nx)=YQ(nq,niqV,1,nx)
              ENDDO
C$OMP END PARALLEL DO
            ENDDO
          ENDIF
          
          IF(IMPLICIT.AND.KTYP36.EQ.1) THEN
            DO nonrlist=1,NRLIST(0)
              nr=NRLIST(nonrlist)
C$OMP PARALLEL DO
C$&     PRIVATE(nq),
C$&     SHARED(niqV,NQR,nr,NWQ,nx,RHS,YQ)
              DO nq=NQR(1,nr),NQR(2,nr) !make sure RHS is set correctly
                IF(NWQ(1,nq).EQ.0) THEN !internal
                  RHS(nq)=YQ(nq,niqV,1,nx)
                ELSE !boundary
                  RHS(nq)=0.0d0
                ENDIF
              ENDDO
C$OMP END PARALLEL DO
            ENDDO
          ENDIF
            
C         Loop over regions
          DO nonrlist=1,NRLIST(0)
            nr=NRLIST(nonrlist)
            ERROR_FLAG=.FALSE.
C           Loop over grid points
C$OMP PARALLEL DO
C$&     PRIVATE(nnq,nq,PSTIMULUS),
C$&     FIRSTPRIVATE(PARAMETERS),
C$&     SHARED(BIDOMAIN,IMPLICIT,maqdt,maqp1t0,maqp1t1,
C$&            maqp1i,maqp2t0,maqp2t1,maqp2i,NENQ,niqV,
C$&            nnq_min,nr,NQS,NQXI,NWQ,nx,nx_ext,NXQ,AQ,
C$&            CQ,GCHQ,GUQ,PROPQ,RHS,T,YLIST,YQ,
C$&            ERROR_FLAG,IOER,IODI)
C$&     REDUCTION(+:NSOL)
            DO nnq=1,NWQ(5,0)
              IF(.NOT.ERROR_FLAG) THEN
                nq=NWQ(5,nnq)
                NSOL=NSOL+1
                IF(NWQ(1,nq).EQ.0) THEN !internal
                  CALL CALC_ADAMS_PSEUDO_STIMULUS(niqV,nq,NQGP(0,nq),
     '              NQGP_PIVOT(1,nq),nx_ext,nx,NQGW(1,nq),PSTIMULUS,
     '              THETA(1),YQ,BIDOMAIN,ERROR,*200)
                  PARAMETERS(2)=CQ(1,nq)
                  PARAMETERS(97)=PSTIMULUS
                  PARAMETERS(98)=CQ(2,nq)
                  
                  CALL INTEGRATOR(ADAMS_IWORK(1,nq,nx),AII,AIO,CONTROL,
     '              L_AII,L_AIO,L_ARI,L_ARO,L_CONTROL,L_MODEL,
     '              L_PARAMETERS, maqdt,maqp1t0,maqp1t1,maqp1i,
     '              maqp2t0,maqp2t1,maqp2i,MODEL,YLIST(0),nq,
     '              nr,nx,YLIST,ADAMS_WORK(1,nq,nx),AQ,ARI,ARO,CQ,
     '              DIFFUSION,PARAMETERS,T,XQ,YQ,ERROR,*200)
                  
                  IF(IMPLICIT) RHS(nq)=YQ(nq,niqV,1,nx)
                ENDIF
                
                GOTO 202
 200            CONTINUE
C$OMP CRITICAL(ADAMS_MARCH8_2)
                ERROR_FLAG=.TRUE.
                WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
                CALL WRITES(IOER,OP_STRING,ERROR,*201)
                WRITE(OP_STRING,'(/'' >>An error occurred - '
     '            //'results may be unreliable!'')')
                CALL WRITES(IODI,OP_STRING,ERROR,*201)
 201            CONTINUE
C$OMP END CRITICAL(ADAMS_MARCH8_2)
 202            CONTINUE
              ENDIF !.NOT.ERROR_FLAG  
            ENDDO !grid points
C$OMP END PARALLEL DO
            IF(ERROR_FLAG) GOTO 9999

            
C!!! Note this has problems for three grid point wide geometries.
C$OMP PARALLEL DO
C$&     PRIVATE(nq,PHI1,PHI2),
C$&     SHARED(IMPLICIT,niqV,NQR,NWQ,nx,RHS,YQ)
            DO nq=NQR(1,nr),NQR(2,nr)
              IF(NWQ(1,nq).NE.0) THEN !boundary grid point
                IF(IMPLICIT) THEN
                  RHS(nq)=0.0d0
                ELSE
                  PHI1=YQ(NWQ(1,nq),niqV,1,nx)
                  PHI2=YQ(NWQ(2,nq),niqV,1,nx)
                  YQ(nq,niqV,1,nx)=(4.0d0*PHI1-PHI2)/3.0d0
                ENDIF
              ENDIF
            ENDDO !grid points          
          ENDDO !region
C$OMP END PARALLEL DO
          
          IF(IMPLICIT) THEN
C           Solve the system for the transmembrane potential
            CALL SOLVE_SYSTEM(ISC_GKK(1,nx),ISR_GKK(1,nx),1,NQT,
     '        NQT,NQT,NQT,NZZT(1,NRLIST(1),nx),NZ_GKK_M,NISC_GKKM,
     '        NISR_GKKM,IWRIT4(NRLIST(1),nx),PRECON_CODE(nx),
     '        SOLVEROPTION(nx),SPARSEGKK(nx),SPARSEGKK(nx),1,
     '        GKK(1,nx),RHS,YQ(1,niqV,1,nx),FIRST_A,UPDATE_MATRIX,
     '          X_INIT,nx,ERROR,*9999)
            UPDATE_MATRIX=.FALSE.
          ENDIF
          
          IF(BIDOMAIN) THEN
            IF(CHMTRIX) THEN
              IF(.NOT.RESTART) FIRST_A=.TRUE.
              UPDATE_MATRIX=.TRUE.
              CHMTRIX=.FALSE.
            ENDIF
            
C           Generate the RHS vector for extracellular potential solution
            CALL CALC_ADAMS_EXT_RHS(niqV,niqBNDRY,NQGP,NQGP_PIVOT,NQGW,
     '        NRLIST,NWQ,nx_ext,nx,DT,PROPQ,RHS,YQ,CQ,
     '        FIXQ(1,1,nx_ext),ERROR,*9999)
            
C           Give better guess (maybe do this for the TM potentials
C           as well?)
C$OMP PARALLEL DO
C$&     PRIVATE(nq),
C$&     SHARED(niq_old,niq_old2,niqV,nx_ext,YQ)
            DO nq=1,NQT
              YQ(nq,niqV,1,nx_ext)=YQ(nq,niq_old,1,nx_ext)+0.05d0*
     '          (YQ(nq,niq_old,1,nx_ext)-YQ(nq,niq_old2,1,nx_ext))
            ENDDO
C$OMP END PARALLEL DO
            
C           Solve the system for the extracellular potential
            CALL SOLVE_SYSTEM(ISC_GKK(1,nx_ext),ISR_GKK(1,nx_ext),
     '        1,NQT,NQT,NQT,NQT,NZZT(1,NRLIST(1),nx_ext),NZ_GKK_M,
     '        NISC_GKKM,NISR_GKKM,IWRIT4(NRLIST(1),nx_ext),
     '        PRECON_CODE(nx_ext),SOLVEROPTION(nx_ext),
     '        SPARSEGKK(nx_ext),SPARSEGKK(nx_ext),1,GKK(1,nx_ext),RHS,
     '        YQ(1,niqV,1,nx_ext),FIRST_A,UPDATE_MATRIX,X_INIT,nx_ext,
     '        ERROR,*9999)
            UPDATE_MATRIX=.FALSE.
          ENDIF !bidomain        
            
          IF(DOP) THEN
            WRITE(OP_STRING,'('' Time: '',F12.6)') T
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            DO nq=1,NQT
              WRITE(OP_STRING,'('' nq, Vm:'',I6,F12.6)') nq,
     '          YQ(nq,niqV,1,nx)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDIF
            
          IF(ITYP3(NRLIST(1),nx).EQ.3.AND.KTYP33.EQ.2) THEN !VCDC
C$OMP PARALLEL DO
C$&     PRIVATE(nq),
C$&     SHARED(DT,niqDV,niq_old,niqV,nx,YQ)
            DO nq=1,NQT
              YQ(nq,niqDV,1,nx)=(YQ(nq,niqV,1,nx)-
     '          YQ(nq,niq_old,1,nx))/DT
            ENDDO
C$OMP END PARALLEL DO
          ENDIF
                                  
          T=T+DT          

C         Find the new smallest time step
          OLDDT=DT
          IF(NWQ(5,0).GT.0) THEN
            DT=RMAX
C$OMP PARALLEL DO
C$&     PRIVATE(nq,nnq),
C$&     SHARED(ADAMS_WORK,NWQ,DT,nx)
            DO nnq=1,NWQ(5,0)
              nq=NWQ(5,nnq)
              IF(ADAMS_WORK(1,nq,nx).LT.DT) THEN
C$OMP CRITICAL(ADAMS_MARCH8_3)
                DT=ADAMS_WORK(1,nq,nx)
C$OMP END CRITICAL(ADAMS_MARCH8_3)
              ENDIF
            ENDDO
C$OMP END PARALLEL DO
          ENDIF
          
          IF(IWRIT5(NRLIST(1),nx).EQ.2) THEN
            IF(MOD(NSTEP,IWRIT1(NRLIST(1),nx)).EQ.0) THEN
              CALL CPU_TIMER(CPU_USER,TIME_STOP2)
              ELAPSED_TIME=TIME_STOP2(1)-TIME_START2(1)
              WRITE(OP_STRING,'(/,1X,''Time (ms): '',F9.3,'
     '          //''', DT (us): '',F8.1,'' Iters: '',I8,'
     '          //''' Active: '',I6,'' Cpu (s): '',F6.2)')
     '          T,DT*1000.0d0,NSTEP,NWQ(5,0),ELAPSED_TIME
      	      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              CALL CPU_TIMER(CPU_USER,TIME_START2)
            ENDIF
          ENDIF
          
          IF(T.GE.TFINISH) CONTINU=.FALSE.
        ENDDO !time
C       End on main time integration loop

        T_RESTART=T
        TINCR=DT

        IF(IWRIT5(NRLIST(1),nx).GE.1) THEN
          CALL CPU_TIMER(CPU_USER,TIME_STOP)
          ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)
          WRITE(OP_STRING,'(1X,I6,'' iterations : '',F8.2,'
     '      //'''s cpu'')') NSTEP,ELAPSED_TIME
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

      ELSE !iterate

        CALL ASSERT(NXLIST(0).GE.4,
     '    '>>You need 4 classes for coupled bidomain',ERROR,*9999)
        
        IF(FIRSTITER) THEN
          FIRSTITER=.FALSE.
          UPDATE_MATRIX=.TRUE.
          FIRST_A=.TRUE.
          DO nq=1,NQT
            FIXQ(nq,3,nx_upd)=.TRUE.
          ENDDO
        ELSE
          UPDATE_MATRIX=.FALSE.
          FIRST_A=.FALSE.
        ENDIF

        NRTEMP(0)=1
        NRTEMP(1)=NRLIST(1)

C       Generate the RHS vector for extracellular potential solution
        CALL GEN_EXT_RHS(niqV,niqBNDRY,NQGP,NQGP_PIVOT,NQGW,NRTEMP,
     '    NWQ,nx_ext,nx,DT,PROPQ,RHS,YQ,CQ,FIXQ(1,1,nx_upd),
     '    ERROR,*9999)

C       Update RHS from YP fluxes
        CALL GEN_GRID_BEM_RHS(IBT,IDO,INP,ISC_GKK(1,nx_upd),
     '    ISR_GKK(1,nx_upd),NBH,NBJ,NEELEM,NENQ,NHE,NHP,NKE,NKH,NPF,
     '    NP_INTERFACE,NPNE,NPNODE,NQNE,NQS,NQXI,NVHE,NVHP,NYNE,NYNP,
     '    NRLIST,NW,NWQ,nx_torso,nx_upd,NXQ,CURVCORRECT,GKK(1,nx_upd),
     '    RHS,SE,YP,ZA,ZE,ZP,UPDATE_MATRIX,ERROR,*9999)
        
        CALL SOLVE_SYSTEM(ISC_GKK(1,nx_upd),ISR_GKK(1,nx_upd),1,NQT,NQT,
     '    NQT,NQT,NZZT(1,NRLIST(1),nx_upd),NZ_GKK_M,NISC_GKKM,NISR_GKKM,
     '    IWRIT4(NRLIST(1),nx_upd),PRECON_CODE(nx_upd),
     '    SOLVEROPTION(nx_upd),SPARSEGKK(nx_upd),SPARSEGKK(nx_upd),1,
     '    GKK(1,nx_upd),RHS,YQ(1,niqV,1,nx_ext),FIRST_A,UPDATE_MATRIX,
     '    X_INIT,nx_upd,ERROR,*9999)
      ENDIF
      
      CALL EXITS('ADAMS_MARCH8')
      RETURN
 9999 CALL ERRORS('ADAMS_MARCH8',ERROR)
      CALL EXITS('ADAMS_MARCH8')
      RETURN 1
      END


      SUBROUTINE ASSEMBLE10(ISC_GKK,ISR_GKK,NEELEM,NENQ,NLATNE,NLATNQ,
     &  NLATPNQ,NQNLAT,NLQ,NQGP,NQGP_PIVOT,NQS,NQXI,NRLIST,NWQ,
     &  nx_ext,nx_trans,nx_upd,NXQ,AQ,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,GCHQ,
     &  GUQ,GKK,NEWDT,NQGW,PROPQ,XQ,BIDOMAIN,COUPBID,FIRST_A,FIXQ,
     &  IMPLICIT,UPDATE,UPDATEDT,UPDATE_MATRIX,SOLVEEIGHTPROBLEM,
     &  ERROR,*)

C#### Subroutine: ASSEMBLE10
C###  Description:
C###    ASSEMBLE10 writes directly into the compressed row storage
C###    arrays. It generates the constraint matrix for the implicit
C###    solution of grid activation problems.
C***  Created by Martin Buist, May 1997

      IMPLICIT NONE

      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b12.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:grid00.cmn'
      INCLUDE 'cmiss$reference:iwrit00.cmn'
      INCLUDE 'cmiss$reference:solv00.cmn'
      INCLUDE 'cmiss$reference:time02.cmn'
      INCLUDE 'cmiss$reference:tol00.cmn'

!     Parameter list
      INTEGER ISC_GKK(NISC_GKKM,NXM),ISR_GKK(NISR_GKKM,NXM),
     &  NEELEM(0:NE_R_M,0:NRM),NENQ(0:8,NQM),NLATNE(NEQM+1),
     &  NLATNQ(NEQM*NQEM),NLATPNQ(NQM),NQNLAT(NEQM*NQEM),NLQ(NQM),
     &  NQGP(0:NQGM,NQM),NQGP_PIVOT(NQGM,NQM),NQS(NEQM),
     &  NQXI(0:NIM,NQSCM),NRLIST(0:NRM),NWQ(8,0:NQM),nx_ext,nx_trans,
     &  nx_upd,NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      REAL*8 AQ(NMAQM,NQM),CQ(NMM,NQM),DNUDXQ(3,3,NQM),DXDXIQ(3,3,NQM),
     &  DXDXIQ2(3,3,NQM),GCHQ(3,NQM),GUQ(3,3,NQM),GKK(NZ_GKK_M,NXM),
     &  NEWDT,NQGW(NQGM,NQM),PROPQ(3,3,4,2,NQM),XQ(NJM,NQM)
      CHARACTER ERROR*(*)
      LOGICAL BIDOMAIN,COUPBID,FIRST_A,FIXQ(NYQM,NIYFIXM,NXM),IMPLICIT,
     &  SOLVEEIGHTPROBLEM,UPDATE,UPDATEDT,UPDATE_MATRIX
!     Local variables
      INTEGER COUNTINT,COUNTEXT,DUMMY_LIST(0:1),IJ,IK,
     &  maxrow,i,j,maxcol,nii,nij,nik,NITB,nq,nq1,nq2,nr,nrr,nzero,
     &  nzero2,nzz,PLACEINT,PLACEEXT,PLACE_PIVOT(22),TEMPPLACE(22)
      REAL ELAPSED_TIME,TIME_START(1),TIME_STOP(1)
      REAL*8 COEFFSEXT(NQGM),DTCMAMTHETA,LATTCOEFS(NQGM),M(9,NQGM),
     &  TEMPCOEFF(22),THETACOEF
      LOGICAL NOTIMP

      CALL ENTERS('ASSEMBLE10',*9999)

      CALL CPU_TIMER(CPU_USER,TIME_START)      
      
      IF(.NOT.UPDATE) THEN
        IF(.NOT.COUPBID) THEN
          IF(USE_LAT.EQ.0) THEN !not needed for lattice scheme
            CALL GET_FD_POINTS(NEELEM,NENQ,NLQ,NQGP,NQGP_PIVOT,NQS,NQXI,
     '        NRLIST,NWQ,NXQ,ERROR,*9999)
          ENDIF
        ENDIF
      ENDIF

C What is this statement for ?
      IF(THETA(1).LT.ZERO_TOL) THEN
        THETACOEF=1.0d0
      ELSE
        THETACOEF=THETA(1)
      ENDIF
    
      !Check array sizes
      IF(.NOT.UPDATE) THEN
        IF(IMPLICIT) THEN
          NOT(1,1,NRLIST(1),nx_trans)=NQT
          NOT(2,1,NRLIST(1),nx_trans)=NQT
          IF(SPARSEGKK(nx_trans).EQ.0) THEN !no sparsity
            CALL ASSERT(USE_LAT.EQ.0,'>>Lattice method not implemented '
     '        //'for no sparsity',ERROR,*9999)
            nzz=NQT*NQT
            NZZT(1,NRLIST(1),nx_trans)=nzz
            IF(NZ_GKK_M.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ENDIF
          ELSE IF(SPARSEGKK(nx_trans).EQ.1) THEN !compressed row
            IF(USE_LAT.EQ.0) THEN !these checks done later
              nzz=0
              maxrow=0
              DO nrr=1,NRLIST(0)
                nr=NRLIST(nrr)
                DO nq=NQR(1,nr),NQR(2,nr)
                  DO nzero=1,NQGP(0,nq)
                    nzz=nzz+1
                  ENDDO !nzero
                ENDDO !nq
                IF(NQR(2,nr)+1.GT.maxrow) maxrow=NQR(2,nr)+1
                NZZT(1,nr,nx_trans)=nzz
              ENDDO !nr
              IF(NZ_GKK_M.LT.nzz) THEN
                WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
                GOTO 9999
              ELSE IF(NISC_GKKM.LT.nzz) THEN
                WRITE(ERROR,'(''>>Increase ISC_GKK_M to >= '',I12)') nzz
                GOTO 9999
              ELSE IF(NISR_GKKM.LT.maxrow) THEN
                WRITE(ERROR,'(''>>Increase ISR_GKK_M to >= '',I12)') 
     '            maxrow
                GOTO 9999
              ENDIF
            ENDIF
            ISR_GKK(1,nx_trans)=1
          ELSE IF(SPARSEGKK(nx_trans).EQ.2) THEN !row/col
            IF(USE_LAT.EQ.0) THEN !these checks done later
              nzz=0
              DO nrr=1,NRLIST(0)
                nr=NRLIST(nrr)
                DO nq=NQR(1,nr),NQR(2,nr)
                  DO nzero=1,NQGP(0,nq)
                    nzz=nzz+1
                  ENDDO !nzero
                ENDDO !nq
                NZZT(1,nr,nx_trans)=nzz
              ENDDO !nr
              IF(NZ_GKK_M.LT.nzz) THEN
                WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
                GOTO 9999
              ELSE IF(NISR_GKKM.LT.nzz) THEN
                WRITE(ERROR,'(''>>Increase ISR_GKK_M to >= '',I12)') nzz
                GOTO 9999
              ELSE IF(NISC_GKKM.LT.nzz) THEN
                WRITE(ERROR,'(''>>Increase ISC_GKK_M to >= '',I12)') nzz
                GOTO 9999
              ENDIF
            ENDIF
          ELSE IF(SPARSEGKK(nx_trans).EQ.3) THEN !compressed column
            CALL ASSERT(USE_LAT.EQ.0,'>>Lattice method not implemented '
     &        //'for compressed column sparsity',ERROR,*9999)
            nzz=0
            maxcol=0
            DO nrr=1,NRLIST(0)
              nr=NRLIST(nrr)
              DO nq=NQR(1,nr),NQR(2,nr)
                DO nzero=1,NQGP(0,nq)
                  nzz=nzz+1
                ENDDO !nzero
              ENDDO !nq
              IF(NQR(2,nr)+1.GT.maxcol) maxcol=NQR(2,nr)+1
              NZZT(1,nr,nx_trans)=nzz
            ENDDO !nr
            IF(NZ_GKK_M.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ELSE IF(NISR_GKKM.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase ISR_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ELSE IF(NISC_GKKM.LT.maxcol) THEN
              WRITE(ERROR,'(''>>Increase ISC_GKK_M to >= '',I12)')
     '          maxcol
              GOTO 9999
            ENDIF
          ELSE IF(SPARSEGKK(nx_trans).EQ.4) THEN
            ERROR='>>Not implemented'
            GOTO 9999
          ELSE IF(SPARSEGKK(nx_trans).EQ.5) THEN !row/col #2
            CALL ASSERT(USE_LAT.EQ.0,'>>Lattice method not implemented '
     &        //'for compressed row/column sparsity',ERROR,*9999)
            nzz=0
            DO nrr=1,NRLIST(0)
              nr=NRLIST(nrr)
              DO nq=NQR(1,nr),NQR(2,nr)
                DO nzero=1,NQGP(0,nq)
                  nzz=nzz+1
                ENDDO !nzero
              ENDDO !nq
              NZZT(1,nr,nx_trans)=nzz
            ENDDO !nr
            IF(NZ_GKK_M.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ELSE IF(NISC_GKKM.LT.2*nzz) THEN
              WRITE(ERROR,'(''>>Increase ISC_GKK_M to >= '',I12)') 2*nzz
              GOTO 9999
            ENDIF
          ENDIF
        ENDIF
        IF(BIDOMAIN) THEN
          NOT(1,1,NRLIST(1),nx_ext)=NQT
          NOT(2,1,NRLIST(1),nx_ext)=NQT
          IF(SPARSEGKK(nx_ext).EQ.0) THEN !no sparsity
            CALL ASSERT(USE_LAT.EQ.0,'>>Lattice method not implemented '
     '        //'for no sparsity',ERROR,*9999)
            nzz=NQT*NQT
            NZZT(1,NRLIST(1),nx_ext)=nzz
            IF(NZ_GKK_M.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ENDIF
          ELSE IF(SPARSEGKK(nx_ext).EQ.1) THEN !compressed row
            IF(USE_LAT.EQ.0) THEN !these checks done later
              nzz=0
              maxrow=0
              DO nrr=1,NRLIST(0)
                nr=NRLIST(nrr)
                DO nq=NQR(1,nr),NQR(2,nr)
                  DO nzero=1,NQGP(0,nq)
                    nzz=nzz+1
                  ENDDO !nzero
                ENDDO !nq
                NZZT(1,nr,nx_ext)=nzz
                IF(NQR(2,nr)+1.GT.maxrow) maxrow=NQR(2,nr)+1
              ENDDO !nr
              IF(NZ_GKK_M.LT.nzz) THEN
                WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
                GOTO 9999
              ELSE IF(NISC_GKKM.LT.nzz) THEN
                WRITE(ERROR,'(''>>Increase ISC_GKK_M to >= '',I12)') nzz
                GOTO 9999
              ELSE IF(NISR_GKKM.LT.maxrow) THEN
                WRITE(ERROR,'(''>>Increase ISR_GKK_M to >= '',I12)') 
     '            maxrow
                GOTO 9999
              ENDIF
            ENDIF
          ELSE IF(SPARSEGKK(nx_ext).EQ.2) THEN !row/col
            IF(USE_LAT.EQ.0) THEN !these checks done later
              nzz=0
              DO nrr=1,NRLIST(0)
                nr=NRLIST(nrr)
                DO nq=NQR(1,nr),NQR(2,nr)
                  DO nzero=1,NQGP(0,nq)
                    nzz=nzz+1
                  ENDDO !nzero
                ENDDO !nq
                NZZT(1,nr,nx_ext)=nzz
              ENDDO !nr
              IF(NZ_GKK_M.LT.nzz) THEN
                WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
                GOTO 9999
              ELSE IF(NISR_GKKM.LT.nzz) THEN
                WRITE(ERROR,'(''>>Increase ISR_GKK_M to >= '',I12)') nzz
                GOTO 9999
              ELSE IF(NISC_GKKM.LT.nzz) THEN
                WRITE(ERROR,'(''>>Increase ISC_GKK_M to >= '',I12)') nzz
                GOTO 9999
              ENDIF
            ENDIF
          ELSE IF(SPARSEGKK(nx_ext).EQ.3) THEN !compressed column
            CALL ASSERT(USE_LAT.EQ.0,'>>Lattice method not implemented '
     &        //'for compressed column sparsity',ERROR,*9999)
            nzz=0
            maxcol=0
            DO nrr=1,NRLIST(0)
              nr=NRLIST(nrr)
              DO nq=NQR(1,nr),NQR(2,nr)
                DO nzero=1,NQGP(0,nq)
                  nzz=nzz+1
                ENDDO !nzero
              ENDDO !nq
              IF(NQR(2,nr)+1.GT.maxcol) maxcol=NQR(2,nr)+1
              NZZT(1,nr,nx_ext)=nzz
            ENDDO !nr
            IF(NZ_GKK_M.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ELSE IF(NISR_GKKM.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase ISR_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ELSE IF(NISC_GKKM.LT.maxcol) THEN
              WRITE(ERROR,'(''>>Increase ISC_GKK_M to >= '',I12)')
     '          maxcol
              GOTO 9999
            ENDIF
          ELSE IF(SPARSEGKK(nx_ext).EQ.4) THEN
            ERROR='>>Not implemented'
            GOTO 9999
          ELSE IF(SPARSEGKK(nx_ext).EQ.5) THEN !row/col #2
            CALL ASSERT(USE_LAT.EQ.0,'>>Lattice method not implemented '
     &        //'for compressed row/column sparsity',ERROR,*9999)
            nzz=0
            DO nrr=1,NRLIST(0)
              nr=NRLIST(nrr)
              DO nq=NQR(1,nr),NQR(2,nr)
                DO nzero=1,NQGP(0,nq)
                  nzz=nzz+1
                ENDDO !nzero
              ENDDO !nq
              NZZT(1,nr,nx_ext)=nzz
            ENDDO !nr
            IF(NZ_GKK_M.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ELSE IF(NISC_GKKM.LT.2*nzz) THEN
              WRITE(ERROR,'(''>>Increase ISC_GKK_M to >= '',I12)') 2
     '          *nzz
              GOTO 9999
            ENDIF
          ENDIF
        ENDIF
        IF(COUPBID) THEN
          CALL ASSERT(USE_LAT.EQ.0,'>>Lattice method not implemented '
     &      //'for coupled bidomain',ERROR,*9999)
          NOT(1,1,NRLIST(1),nx_upd)=NQT
          NOT(2,1,NRLIST(1),nx_upd)=NQT
          IF(SPARSEGKK(nx_upd).EQ.0) THEN !no sparsity
            nzz=NQT*NQT
            NZZT(1,NRLIST(1),nx_upd)=nzz
            IF(NZ_GKK_M.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ENDIF
          ELSE IF(SPARSEGKK(nx_upd).EQ.1) THEN !compressed row
            nzz=0
            maxrow=0
            DO nrr=1,NRLIST(0)
              nr=NRLIST(nrr)
              DO nq=NQR(1,nr),NQR(2,nr)
                DO nzero=1,NQGP(0,nq)
                  nzz=nzz+1
                ENDDO !nzero
              ENDDO !nq
              NZZT(1,nr,nx_upd)=nzz
              IF(NQR(2,nr)+1.GT.maxrow) maxrow=NQR(2,nr)+1
            ENDDO !nr
            IF(NZ_GKK_M.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ELSE IF(NISC_GKKM.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase ISC_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ELSE IF(NISR_GKKM.LT.maxrow) THEN
              WRITE(ERROR,'(''>>Increase ISR_GKK_M to >= '',I12)')
     '          maxrow
              GOTO 9999
            ENDIF
          ELSE IF(SPARSEGKK(nx_upd).EQ.2) THEN !row/col
            nzz=0
            DO nrr=1,NRLIST(0)
              nr=NRLIST(nrr)
              DO nq=NQR(1,nr),NQR(2,nr)
                DO nzero=1,NQGP(0,nq)
                  nzz=nzz+1
                ENDDO !nzero
              ENDDO !nq
              NZZT(1,nr,nx_upd)=nzz
            ENDDO !nr
            IF(NZ_GKK_M.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ELSE IF(NISR_GKKM.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase ISR_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ELSE IF(NISC_GKKM.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase ISC_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ENDIF
          ELSE IF(SPARSEGKK(nx_upd).EQ.3) THEN !compressed column
            nzz=0
            maxcol=0
            DO nrr=1,NRLIST(0)
              nr=NRLIST(nrr)
              DO nq=NQR(1,nr),NQR(2,nr)
                DO nzero=1,NQGP(0,nq)
                  nzz=nzz+1
                ENDDO !nzero
              ENDDO !nq
              IF(NQR(2,nr)+1.GT.maxcol) maxcol=NQR(2,nr)+1
              NZZT(1,nr,nx_upd)=nzz
            ENDDO !nr
            IF(NZ_GKK_M.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ELSE IF(NISR_GKKM.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase ISR_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ELSE IF(NISC_GKKM.LT.maxcol) THEN
              WRITE(ERROR,'(''>>Increase ISC_GKK_M to >= '',I12)')
     '          maxcol
              GOTO 9999
            ENDIF
          ELSE IF(SPARSEGKK(nx_upd).EQ.4) THEN
            ERROR='>>Not implemented'
            GOTO 9999
          ELSE IF(SPARSEGKK(nx_upd).EQ.5) THEN !row/col #2
            nzz=0
            DO nrr=1,NRLIST(0)
              nr=NRLIST(nrr)
              DO nq=NQR(1,nr),NQR(2,nr)
                DO nzero=1,NQGP(0,nq)
                  nzz=nzz+1
                ENDDO !nzero
              ENDDO !nq
              NZZT(1,nr,nx_upd)=nzz
            ENDDO !nr
            IF(NZ_GKK_M.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ELSE IF(NISC_GKKM.LT.2*nzz) THEN
              WRITE(ERROR,'(''>>Increase ISC_GKK_M to >= '',I12)') 2*nzz
              GOTO 9999
            ENDIF
          ELSE
            ERROR='>>Unknown sparsity type for GKK'
            GOTO 9999
          ENDIF
        ENDIF
C        CALL ASSERT(NOM.GE.NQT,'>>Increase NOM, min. NQT',ERROR,*9999)
      ENDIF

      !Initialising
      IF(.NOT.UPDATE) THEN
        IF(IMPLICIT) THEN
          IF(SPARSEGKK(nx_trans).EQ.0) THEN
C$OMP PARALLEL DO
C$OMP&  PRIVATE(nq),
C$OMP&  SHARED(GKK,nx_trans)
            DO nq=1,NQT*NQT
              GKK(nq,nx_trans)=0.0d0
            ENDDO !nq
C$OMP END PARALLEL DO
          ELSE IF(SPARSEGKK(nx_trans).EQ.1) THEN
            ISR_GKK(1,nx_trans)=1
          ELSE IF(SPARSEGKK(nx_trans).EQ.3) THEN
            ISC_GKK(1,nx_trans)=1
          ENDIF
        ENDIF
        IF(BIDOMAIN) THEN
          IF(SPARSEGKK(nx_ext).EQ.0) THEN
C$OMP PARALLEL DO
C$OMP&  PRIVATE(nq),
C$OMP&  SHARED(GKK,nx_ext)
            DO nq=1,NQT*NQT
              GKK(nq,nx_ext)=0.0d0
            ENDDO !nq
C$OMP END PARALLEL DO
          ELSE IF(SPARSEGKK(nx_ext).EQ.1) THEN
            ISR_GKK(1,nx_ext)=1
          ELSE IF(SPARSEGKK(nx_ext).EQ.3) THEN
            ISC_GKK(1,nx_ext)=1
          ENDIF
        ENDIF
        IF(COUPBID) THEN
          IF(SPARSEGKK(nx_upd).EQ.0) THEN
C$OMP PARALLEL DO
C$OMP&  PRIVATE(nq),
C$OMP&  SHARED(GKK,nx_upd)
            DO nq=1,NQT*NQT
              GKK(nq,nx_upd)=0.0d0
            ENDDO !nq
C$OMP END PARALLEL DO
          ELSE IF(SPARSEGKK(nx_upd).EQ.1) THEN
            ISR_GKK(1,nx_upd)=1
          ELSE IF(SPARSEGKK(nx_upd).EQ.3) THEN
            ISC_GKK(1,nx_upd)=1
          ENDIF
        ENDIF
      ENDIF

      IF(UPDATEDT) THEN
        IF(IMPLICIT) THEN
          !Main loop
          DO nrr=1,NRLIST(0)
            nr=NRLIST(nrr)

C$OMP PARALLEL DO
C$OMP&  PRIVATE(COUNTINT,DTCMAMTHETA,IJ,IK,nii,nij,nik,NITB,nq,nq1,nq2,
C$OMP&    nzero,PLACEINT,TEMPCOEFF,TEMPPLACE),
C$OMP&  SHARED(ISC_GKK,ISR_GKK,GKK,NQGP,NQR,NQT,nr,NWQ,nx_trans,
C$OMP&    SPARSEGKK,CQ,DT)
            DO nq=NQR(1,nr),NQR(2,nr) !create line for each nq
              DTCMAMTHETA=DT/(CQ(1,nq)*CQ(2,nq))*THETACOEF
              IF(SPARSEGKK(nx_trans).EQ.0) THEN !no sparsity
                IF(NWQ(1,nq).EQ.0) THEN !internal
                  DO nzero=1,NQGP(0,nq)
                    GKK(((NQGP(nzero,nq)-1)*NQT)+nq,nx_trans)=
     '                -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTCMAMTHETA
                    IF(NQGP(nzero,nq).EQ.nq)
     '                GKK(((nq-1)*NQT)+nq,nx_trans)=
     '                GKK(((nq-1)*NQT)+nq,nx_trans)+1.0d0
                  ENDDO
                ENDIF
              ELSE IF(SPARSEGKK(nx_trans).EQ.1) THEN !compressed row
                IF(NWQ(1,nq).EQ.0) THEN !internal
                  PLACEINT=ISR_GKK(nq,nx_trans)
                  DO nzero=1,NQGP(0,nq)
                    GKK(PLACEINT,nx_trans)=
     '                -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTCMAMTHETA
                    IF(NQGP(nzero,nq).EQ.nq) GKK(PLACEINT,nx_trans)=
     '                GKK(PLACEINT,nx_trans)+1.0d0
                    PLACEINT=PLACEINT+1
                  ENDDO
                ENDIF
              ELSE IF(SPARSEGKK(nx_trans).EQ.2) THEN !row/col
                IF(NWQ(1,nq).EQ.0) THEN !internal
                  PLACEINT=ISC_GKK(nq,nx_trans)
                  DO nzero=1,NQGP(0,nq)
                    GKK(PLACEINT,nx_trans)=
     '                -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTCMAMTHETA
                    IF(NQGP(nzero,nq).EQ.nq) GKK(PLACEINT,nx_trans)=
     '                GKK(PLACEINT,nx_trans)+1.0d0
                    PLACEINT=PLACEINT+1
                  ENDDO !nzero
                ENDIF
              ELSE IF(SPARSEGKK(nx_trans).EQ.3) THEN !compressed column

C!!! THIS NEEDS REDOING FOR MP

                DTCMAMTHETA=NEWDT/(CQ(1,nq)*CQ(2,nq))*THETACOEF

                COUNTINT=0
                NITB=NQXI(0,NQS(NEELEM(1,nr)))
                IK=MAX(0,NITB-2) !zero for 1,2-D, one for 3-D
                IJ=MIN(NITB-1,1) !zero for 1-D, one for 2,3-D

                DO nik=-IK,IK
                  DO nij=-IJ,IJ
                    DO nii=-1,1
                      nq1=NXQ(nii,1,NXQ(nij*2,1,NXQ(nik*3,1,nq,1),1),1)
                      IF(nq1.GT.0) THEN
                        DO nzero=1,NQGP(0,nq1)
                          IF(NQGP(nzero,nq1).EQ.nq) THEN
                            COUNTINT=COUNTINT+1
                            IF(NWQ(1,nq1).EQ.0) THEN !internal
                              TEMPCOEFF(COUNTINT)=
     '                          NQGW(NQGP_PIVOT(nzero,nq1),nq1)*
     '                          DTCMAMTHETA+1.0d0
                            ELSE !external
                              TEMPCOEFF(COUNTINT)=
     '                          NQGW(NQGP_PIVOT(nzero,nq1),nq1)
                            ENDIF
                            TEMPPLACE(COUNTINT)=nq1
                            nq2=NXQ(nii,1,NXQ(nij*2,1,
     '                                    NXQ(nik*3,1,nq1,1),1),1)
                            IF((nq2.GT.0).AND.(nq2.NE.nq1)) THEN
                              DO nzero2=1,NQGP(0,nq2)
                                IF(NQGP(nzero2,nq2).EQ.nq) THEN
                                  COUNTINT=COUNTINT+1
                                  IF(NWQ(1,nq2).EQ.0) THEN !internal
                                    TEMPCOEFF(COUNTINT)=
     '                                NQGW(NQGP_PIVOT(nzero2,nq2),nq2)*
     '                                DTCMAMTHETA
                                  ELSE !external
                                    TEMPCOEFF(COUNTINT)=
     '                                NQGW(NQGP_PIVOT(nzero2,nq2),nq2)
                                  ENDIF
                                  TEMPPLACE(COUNTINT)=nq2
                                ENDIF
                              ENDDO
                            ENDIF
                          ENDIF
                        ENDDO
                      ENDIF
                    ENDDO !nii
                  ENDDO !nij
                ENDDO !nik

                CALL ISORTP(COUNTINT,TEMPPLACE,PLACE_PIVOT)

                DO nzero=1,COUNTINT
                  PLACEINT=PLACEINT+1
                  GKK(PLACEINT,nx_trans)=TEMPCOEFF(PLACE_PIVOT(nzero))
                ENDDO

              ELSE IF(SPARSEGKK(nx_trans).EQ.5) THEN !row/col #2
                PLACEINT=ISC_GKK(nq,nx_trans)
                IF(NWQ(1,nq).EQ.0) THEN !internal
                  DO nzero=1,NQGP(0,nq)
                    GKK(PLACEINT,nx_trans)=
     '                -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTCMAMTHETA
                    IF(NQGP(nzero,nq).EQ.nq) GKK(PLACEINT,nx_trans)=
     '                GKK(PLACEINT,nx_trans)+1.0d0
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
              IF(SPARSEGKK(nx_trans).EQ.0) THEN !no sparsity
                CALL CALC_GRID_COEF(NENQ,NITB,nq,NQGP,NQS,
     '            NQXI,NWQ(1,nq),NXQ(-NIM,0,0,1),nx_ext,nx_trans,
     '            COEFFSEXT,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,
     '            GCHQ(1,nq),GUQ(1,1,nq),
     '            NQGW(1,nq),PROPQ(1,1,1,1,nq),.FALSE.,FIXQ,
     '            IMPLICIT,SOLVEEIGHTPROBLEM,ERROR,*9999)
                IF(NWQ(1,nq).EQ.0) THEN !internal
                  DO nzero=1,NQGP(0,nq)
                    GKK(((NQGP(nzero,nq)-1)*NQT)+nq,nx_trans)=
     '                -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTCMAMTHETA
                    IF(NQGP(nzero,nq).EQ.nq)
     '                GKK(((nq-1)*NQT)+nq,nx_trans)=
     '                GKK(((nq-1)*NQT)+nq,nx_trans)+1.0d0
                  ENDDO
                ELSE
                  DO nzero=1,NQGP(0,nq)
                    GKK(((NQGP(nzero,nq)-1)*NQT)+nq,nx_trans)=
     '                NQGW(NQGP_PIVOT(nzero,nq),nq)
                  ENDDO
                ENDIF
              ELSE IF(SPARSEGKK(nx_trans).EQ.1) THEN !compressed row
                IF(USE_LAT.EQ.0) THEN
                  CALL CALC_GRID_COEF(NENQ,NITB,nq,NQGP,NQS,
     '              NQXI,NWQ(1,nq),NXQ(-NIM,0,0,1),nx_ext,nx_trans,
     '              COEFFSEXT,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,
     '              GCHQ(1,nq),GUQ(1,1,nq),
     '              NQGW(1,nq),PROPQ(1,1,1,1,nq),.FALSE.,FIXQ,
     '              IMPLICIT,SOLVEEIGHTPROBLEM,ERROR,*9999)
                  COUNTINT=0
                  IF(NWQ(1,nq).EQ.0) THEN !internal
                    DO nzero=1,NQGP(0,nq)
                      PLACEINT=PLACEINT+1
                      COUNTINT=COUNTINT+1
                      IF(.NOT.UPDATE) 
     '                  ISC_GKK(PLACEINT,nx_trans)=NQGP(nzero,nq)
                      GKK(PLACEINT,nx_trans)=
     '                  -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTCMAMTHETA
                      IF(NQGP(nzero,nq).EQ.nq) GKK(PLACEINT,nx_trans)=
     '                  GKK(PLACEINT,nx_trans)+1.0d0
                    ENDDO
                  ELSE !external
                    DO nzero=1,NQGP(0,nq)
                      PLACEINT=PLACEINT+1
                      COUNTINT=COUNTINT+1
                      IF(.NOT.UPDATE) 
     '                  ISC_GKK(PLACEINT,nx_trans)=NQGP(nzero,nq)
                      GKK(PLACEINT,nx_trans)=
     '                  NQGW(NQGP_PIVOT(nzero,nq),nq)
                    ENDDO
                  ENDIF
                  IF(.NOT.UPDATE) 
     '              ISR_GKK(nq+1,nx_trans)=ISR_GKK(nq,nx_trans)+COUNTINT
                ELSE !lattice method
                  IF(.NOT.UPDATE) THEN
                    CALL GSUPPORT(NLATNE,NLATNQ,NLATPNQ,nq,NQNLAT,NQS,
     &                NQXI,NWQ,NQGP(0,nq),ERROR,*9999)
                    CALL ASSERT(nq+1.LE.NISR_GKKM,
     &                '>> Increase NISR_GKKM',ERROR,*9999)
                    CALL ASSERT(PLACEINT+NQGP(0,nq).LE.NISC_GKKM,
     &                '>> Increase NISC_GKKM',ERROR,*9999)
                    ISR_GKK(nq+1,nx_trans)=ISR_GKK(nq,nx_trans)+
     &                NQGP(0,nq)
                    DO nzero=1,NQGP(0,nq)
                      PLACEINT=PLACEINT+1
                      ISC_GKK(PLACEINT,nx_trans)=NQGP(nzero,nq)
                    ENDDO
                  ENDIF
                  CALL CALC_LATTICE_WEIGHTS(nq,NQGP(0,nq),M,XQ,
     &              ERROR,*9999)
                  CALL CALC_LATT_COEF_TRANS(NITB,nq,NWQ(1,nq),nx_trans,
     &              NQGP(0,nq),AQ(1,nq),M,NQGW(1,nq),PROPQ,FIXQ,
     &              ERROR,*9999)
                  CALL ASSERT(PLACEINT.LE.NZ_GKK_M,
     &              '>>Increase NZ_GKK_M',ERROR,*9999)
                  DO i=1,NQGP(0,nq)
                    j=PLACEINT-NQGP(0,nq)+i
                    IF(NWQ(1,nq).EQ.0) THEN !internal
                      GKK(j,nx_trans)=-1.0d0*NQGW(i,nq)*DTCMAMTHETA
                      IF(NQGP(i,nq).EQ.nq) THEN
                        GKK(j,nx_trans)=GKK(j,nx_trans)+1.0d0
                      ENDIF
                    ELSE
                      GKK(j,nx_trans)=NQGW(i,nq)
                    ENDIF
                  ENDDO                         
                ENDIF
              ELSE IF(SPARSEGKK(nx_trans).EQ.2) THEN !row/col
                IF(USE_LAT.EQ.0) THEN
                  CALL CALC_GRID_COEF(NENQ,NITB,nq,NQGP,NQS,NQXI,
     &              NWQ(1,nq),NXQ(-NIM,0,0,1),nx_ext,nx_trans,COEFFSEXT,
     &              CQ,DNUDXQ,DXDXIQ,DXDXIQ2,GCHQ(1,nq),GUQ(1,1,nq),
     &              NQGW(1,nq),PROPQ(1,1,1,1,nq),.FALSE.,FIXQ,IMPLICIT,
     &              SOLVEEIGHTPROBLEM,ERROR,*9999)
                  IF(NWQ(1,nq).EQ.0) THEN !internal
                    DO nzero=1,NQGP(0,nq)
                      PLACEINT=PLACEINT+1
                      IF(.NOT.UPDATE) THEN
                        ISR_GKK(PLACEINT,nx_trans)=nq
                        ISC_GKK(PLACEINT,nx_trans)=NQGP(nzero,nq)
                      ENDIF
                      GKK(PLACEINT,nx_trans)=
     '                  -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTCMAMTHETA
                      IF(NQGP(nzero,nq).EQ.nq) GKK(PLACEINT,nx_trans)=
     '                  GKK(PLACEINT,nx_trans)+1.0d0
                    ENDDO !nzero
                  ELSE !external
                    DO nzero=1,NQGP(0,nq)
                      PLACEINT=PLACEINT+1
                      IF(.NOT.UPDATE) THEN
                        ISR_GKK(PLACEINT,nx_trans)=nq
                        ISC_GKK(PLACEINT,nx_trans)=NQGP(nzero,nq)
                      ENDIF
                      GKK(PLACEINT,nx_trans)=NQGW(NQGP_PIVOT(nzero,nq),
     '                  nq)
                    ENDDO !nzero
                  ENDIF
                ELSE !use lattice method
                  IF(.NOT.UPDATE) THEN
                    CALL GSUPPORT(NLATNE,NLATNQ,NLATPNQ,nq,NQNLAT,NQS,
     &                NQXI,NWQ,NQGP(0,nq),ERROR,*9999)
                    CALL ASSERT(PLACEINT+NQGP(0,nq).LE.NISR_GKKM,
     &                '>> Increase NISR_GKKM',ERROR,*9999)
                    CALL ASSERT(PLACEINT+NQGP(0,nq).LE.NISC_GKKM,
     &                '>> Increase NISC_GKKM',ERROR,*9999)
                    DO nzero=1,NQGP(0,nq)
                      PLACEINT=PLACEINT+1
                      ISR_GKK(PLACEINT,nx_trans)=nq
                      ISC_GKK(PLACEINT,nx_trans)=NQGP(nzero,nq)
                    ENDDO
                  ENDIF
                  CALL CALC_LATTICE_WEIGHTS(nq,NQGP(0,nq),M,XQ,
     &              ERROR,*9999)
                  CALL CALC_LATT_COEF_TRANS(NITB,nq,NWQ(1,nq),nx_trans,
     &              NQGP(0,nq),AQ(1,nq),M,NQGW(1,nq),PROPQ,FIXQ,
     &              ERROR,*9999)
                  CALL ASSERT(PLACEINT.LE.NZ_GKK_M,
     &              '>>Increase NZ_GKK_M',ERROR,*9999)
                  DO i=1,NQGP(0,nq)
                    j=PLACEINT-NQGP(0,nq)+i
                    IF(NWQ(1,nq).EQ.0) THEN !internal
                      GKK(j,nx_trans)=-1.0d0*NQGW(i,nq)*DTCMAMTHETA
                      IF(NQGP(i,nq).EQ.nq) THEN
                        GKK(j,nx_trans)=GKK(j,nx_trans)+1.0d0
                      ENDIF
                    ELSE
                      GKK(j,nx_trans)=NQGW(i,nq)
                    ENDIF
                  ENDDO                         
                ENDIF !lattice 
              ELSE IF(SPARSEGKK(nx_trans).EQ.3) THEN !compressed column
                COUNTINT=0
                IK=MAX(0,NITB-2) !zero for 1,2-D, one for 3-D
                IJ=MIN(NITB-1,1) !zero for 1-D, one for 2,3-D

                DO nik=-IK,IK
                  DO nij=-IJ,IJ
                    DO nii=-1,1
                      nq1=NXQ(nii,1,NXQ(nij*2,1,NXQ(nik*3,1,nq,1),1),1)
                      IF(nq1.GT.0) THEN
                        DO nzero=1,NQGP(0,nq1)
                          IF(NQGP(nzero,nq1).EQ.nq) THEN
                            CALL CALC_GRID_COEF(NENQ,NITB,nq,NQGP,
     '                        NQS,NQXI,NWQ(1,nq),
     '                        NXQ(-NIM,0,0,1),nx_ext,nx_trans,
     '                        COEFFSEXT,CQ,DNUDXQ,
     '                        DXDXIQ,DXDXIQ2,
     '                        GCHQ(1,nq),GUQ(1,1,nq),NQGW(1,nq),
     '                        PROPQ(1,1,1,1,nq),.FALSE.,FIXQ,
     '                        IMPLICIT,SOLVEEIGHTPROBLEM,ERROR,*9999)
                            COUNTINT=COUNTINT+1
                            IF(NWQ(1,nq1).EQ.0) THEN !internal
                              TEMPCOEFF(COUNTINT)=
     '                          NQGW(NQGP_PIVOT(nzero,nq1),nq1)*
     '                          DTCMAMTHETA+1.0d0
                            ELSE !external
                              TEMPCOEFF(COUNTINT)=
     '                          NQGW(NQGP_PIVOT(nzero,nq1),nq1)
                            ENDIF
                            TEMPPLACE(COUNTINT)=nq1
                            nq2=NXQ(nii,1,NXQ(nij*2,1,
     '                                    NXQ(nik*3,1,nq1,1),1),1)
                            IF((nq2.GT.0).AND.(nq2.NE.nq1)) THEN
                              DO nzero2=1,NQGP(0,nq2)
                                IF(NQGP(nzero2,nq2).EQ.nq) THEN
                                  CALL CALC_GRID_COEF(NENQ,NITB,nq,
     '                              NQGP,NQS,NQXI,NWQ(1,nq),
     '                              NXQ(-NIM,0,0,1),nx_ext,nx_trans,
     '                              COEFFSEXT,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,
     '                              GCHQ(1,nq),GUQ(1,1,nq),NQGW(1,nq),
     '                              PROPQ(1,1,1,1,nq),.FALSE.,
     '                              FIXQ,IMPLICIT,SOLVEEIGHTPROBLEM,
     '                              ERROR,*9999)
                                  COUNTINT=COUNTINT+1
                                  IF(NWQ(1,nq2).EQ.0) THEN !internal
                                    TEMPCOEFF(COUNTINT)=
     '                                -NQGW(NQGP_PIVOT(nzero2,nq2),nq2)*
     '                                DTCMAMTHETA
                                  ELSE !external
                                    TEMPCOEFF(COUNTINT)=
     '                                NQGW(NQGP_PIVOT(nzero2,nq2),nq2)
                                  ENDIF
                                  TEMPPLACE(COUNTINT)=nq2
                                ENDIF !NQGP
                              ENDDO !nzero2
                            ENDIF !nq2
                          ENDIF !NQGP
                        ENDDO !nzero
                      ENDIF !nq1
                    ENDDO !nii
                  ENDDO !nij
                ENDDO !nik

                CALL ISORTP(COUNTINT,TEMPPLACE,PLACE_PIVOT)

                DO nzero=1,COUNTINT
                  PLACEINT=PLACEINT+1
                  IF(.NOT.UPDATE)
     '              ISR_GKK(PLACEINT,nx_trans)=TEMPPLACE(nzero)
                  GKK(PLACEINT,nx_trans)=TEMPCOEFF(PLACE_PIVOT(nzero))
                ENDDO
                IF(.NOT.UPDATE)
     '            ISC_GKK(nq+1,nx_trans)=ISC_GKK(nq,nx_trans)+COUNTINT
              ELSE IF(SPARSEGKK(nx_trans).EQ.5) THEN !row/col #2
                CALL ASSERT(USE_LAT.EQ.0,
     '            '>>SPARSENESS 5 not implemented for implict schemes',
     '            ERROR,*9999)
                CALL CALC_GRID_COEF(NENQ,NITB,nq,NQGP,NQS,
     '            NQXI,NWQ(1,nq),NXQ(-NIM,0,0,1),nx_ext,nx_trans,
     '            COEFFSEXT,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,
     '            GCHQ(1,nq),GUQ(1,1,nq),
     '            NQGW(1,nq),PROPQ(1,1,1,1,nq),
     '            .FALSE.,FIXQ,IMPLICIT,SOLVEEIGHTPROBLEM,ERROR,*9999)
                IF(NWQ(1,nq).EQ.0) THEN !internal
                  DO nzero=1,NQGP(0,nq)
                    PLACEINT=PLACEINT+1
                    IF(.NOT.UPDATE) ISC_GKK(PLACEINT,nx_trans)=nq
                    GKK(PLACEINT,nx_trans)=
     '                -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTCMAMTHETA
                    IF(NQGP(nzero,nq).EQ.nq) GKK(PLACEINT,nx_trans)=
     '                GKK(PLACEINT,nx_trans)+1.0d0
                  ENDDO
                ELSE !external
                  DO nzero=1,NQGP(0,nq)
                    PLACEINT=PLACEINT+1
                    IF(.NOT.UPDATE) ISC_GKK(PLACEINT,nx_trans)=nq
                    GKK(PLACEINT,nx_trans)=NQGW(NQGP_PIVOT(nzero,nq),nq)
                  ENDDO
                ENDIF
              ENDIF !SPARSEGKK
            ENDDO !nq
          ENDIF !implicit


          IF(BIDOMAIN) THEN
            DO nq=NQR(1,nr),NQR(2,nr) !create line for each nq
              IF(SPARSEGKK(nx_ext).EQ.0) THEN !no sparsity
                CALL CALC_GRID_COEF(NENQ,NITB,nq,NQGP,NQS,
     '            NQXI,NWQ(1,nq),NXQ(-NIM,0,0,1),nx_ext,nx_trans,
     '            COEFFSEXT,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,
     '            GCHQ(1,nq),GUQ(1,1,nq),
     '            NQGW(1,nq),PROPQ(1,1,1,1,nq),
     '            BIDOMAIN,FIXQ,NOTIMP,SOLVEEIGHTPROBLEM,ERROR,*9999)
                DO nzero=1,NQGP(0,nq)
                  GKK(((NQGP(nzero,nq)-1)*NQT)+nq,nx_ext)=
     '              COEFFSEXT(NQGP_PIVOT(nzero,nq))
                ENDDO
              ELSE IF(SPARSEGKK(nx_ext).EQ.1) THEN !compressed row
                IF(USE_LAT.EQ.0) THEN
                  CALL CALC_GRID_COEF(NENQ,NITB,nq,NQGP,NQS,
     '              NQXI,NWQ(1,nq),NXQ(-NIM,0,0,1),nx_ext,nx_trans,
     '              COEFFSEXT,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,
     '              GCHQ(1,nq),GUQ(1,1,nq),
     '              NQGW(1,nq),PROPQ(1,1,1,1,nq),
     '              BIDOMAIN,FIXQ,NOTIMP,SOLVEEIGHTPROBLEM,ERROR,*9999)
                  COUNTEXT=0
                  DO nzero=1,NQGP(0,nq)
                    PLACEEXT=PLACEEXT+1
                    COUNTEXT=COUNTEXT+1
                    IF(.NOT.UPDATE) 
     '                ISC_GKK(PLACEEXT,nx_ext)=NQGP(nzero,nq)
                    GKK(PLACEEXT,nx_ext)=COEFFSEXT(NQGP_PIVOT(nzero,nq))
                  ENDDO
                  IF(.NOT.UPDATE) 
     '              ISR_GKK(nq+1,nx_ext)=ISR_GKK(nq,nx_ext)+COUNTEXT
                ELSE !lattice method
                  IF(.NOT.UPDATE) THEN
                    CALL GSUPPORT(NLATNE,NLATNQ,NLATPNQ,nq,NQNLAT,NQS,
     &                NQXI,NWQ,NQGP(0,nq),ERROR,*9999)
                    CALL ASSERT(PLACEEXT+NQGP(0,nq).LE.NISR_GKKM,
     &                '>> Increase NISR_GKKM',ERROR,*9999)
                    CALL ASSERT(PLACEEXT+NQGP(0,nq).LE.NISC_GKKM,
     &                '>> Increase NISC_GKKM',ERROR,*9999)
                    ISR_GKK(nq+1,nx_ext)=ISR_GKK(nq,nx_ext)+NQGP(0,nq)
                    DO nzero=1,NQGP(0,nq)
                      PLACEEXT=PLACEEXT+1
                      ISC_GKK(PLACEEXT,nx_ext)=NQGP(nzero,nq)
                    ENDDO
                  ENDIF
                  CALL CALC_LATTICE_WEIGHTS(nq,NQGP(0,nq),M,XQ,
     &              ERROR,*9999)
                  CALL CALC_LATT_COEF_EXT(NITB,nq,NWQ(1,nq),nx_ext,
     &              NQGP(0,nq),LATTCOEFS,AQ(1,nq),M,PROPQ,FIXQ,
     &              SOLVEEIGHTPROBLEM,ERROR,*9999)
                  CALL ASSERT(PLACEEXT.LE.NZ_GKK_M,
     &              '>>Increase NZ_GKK_M',ERROR,*9999)
                  DO i=1,NQGP(0,nq)
                    GKK(PLACEEXT-NQGP(0,nq)+i,nx_ext)=LATTCOEFS(i)
                  ENDDO                         
                ENDIF
              ELSE IF(SPARSEGKK(nx_ext).EQ.2) THEN !row/col
                IF(USE_LAT.EQ.0) THEN
                  CALL CALC_GRID_COEF(NENQ,NITB,nq,NQGP,NQS,
     '              NQXI,NWQ(1,nq),NXQ(-NIM,0,0,1),nx_ext,nx_trans,
     '              COEFFSEXT,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,
     '              GCHQ(1,nq),GUQ(1,1,nq),
     '              NQGW(1,nq),PROPQ(1,1,1,1,nq),
     '              BIDOMAIN,FIXQ,NOTIMP,SOLVEEIGHTPROBLEM,ERROR,*9999)
                  DO nzero=1,NQGP(0,nq)
                    PLACEEXT=PLACEEXT+1
                    IF(.NOT.UPDATE) THEN
                      ISR_GKK(PLACEEXT,nx_ext)=nq
                      ISC_GKK(PLACEEXT,nx_ext)=NQGP(nzero,nq)
                    ENDIF
                    GKK(PLACEEXT,nx_ext)=COEFFSEXT(NQGP_PIVOT(nzero,nq))
                  ENDDO
                ELSE !generate the matrix using the lattice scheme
                  IF(.NOT.UPDATE) THEN
                    CALL GSUPPORT(NLATNE,NLATNQ,NLATPNQ,nq,NQNLAT,NQS,
     &                NQXI,NWQ,NQGP(0,nq),ERROR,*9999)
                    CALL ASSERT(PLACEEXT+NQGP(0,nq).LE.NISR_GKKM,
     &                '>> Increase NISR_GKKM',ERROR,*9999)
                    CALL ASSERT(PLACEEXT+NQGP(0,nq).LE.NISC_GKKM,
     &                '>> Increase NISC_GKKM',ERROR,*9999)
                    DO nzero=1,NQGP(0,nq)
                      PLACEEXT=PLACEEXT+1
                      ISR_GKK(PLACEEXT,nx_ext)=nq
                      ISC_GKK(PLACEEXT,nx_ext)=NQGP(nzero,nq)
                    ENDDO
                  ENDIF
                  CALL CALC_LATTICE_WEIGHTS(nq,NQGP(0,nq),M,XQ,
     &              ERROR,*9999)
                  CALL CALC_LATT_COEF_EXT(NITB,nq,NWQ(1,nq),nx_ext,
     &              NQGP(0,nq),LATTCOEFS,AQ(1,nq),M,PROPQ,FIXQ,
     &              SOLVEEIGHTPROBLEM,ERROR,*9999)
                  CALL ASSERT(PLACEEXT.LE.NZ_GKK_M,
     &              '>>Increase NZ_GKK_M',ERROR,*9999)
                  DO i=1,NQGP(0,nq)
                    GKK(PLACEEXT-NQGP(0,nq)+i,nx_ext)=LATTCOEFS(i)
                  ENDDO
                ENDIF !lattice
              ELSE IF(SPARSEGKK(nx_ext).EQ.3) THEN !compressed column
                COUNTEXT=0
                IK=MAX(0,NITB-2) !zero for 1,2-D, one for 3-D
                IJ=MIN(NITB-1,1) !zero for 1-D, one for 2,3-D

                DO nik=-IK,IK
                  DO nij=-IJ,IJ
                    DO nii=-1,1
                      nq1=NXQ(nii,1,NXQ(nij*2,1,NXQ(nik*3,1,nq,1),1),1)
                      IF(nq1.GT.0) THEN
                        DO nzero=1,NQGP(0,nq1)
                          IF(NQGP(nzero,nq1).EQ.nq) THEN
                            CALL CALC_GRID_COEF(NENQ,NITB,nq,NQGP,
     '                        NQS,NQXI,NWQ(1,nq),
     '                        NXQ(-NIM,0,0,1),nx_ext,nx_trans,
     '                        COEFFSEXT,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,
     '                        GCHQ(1,nq),
     '                        GUQ(1,1,nq),
     '                        NQGW(1,nq),PROPQ(1,1,1,1,nq),
     '                        BIDOMAIN,FIXQ,NOTIMP,SOLVEEIGHTPROBLEM,
     '                        ERROR,*9999)
                            COUNTEXT=COUNTEXT+1
                            TEMPCOEFF(COUNTEXT)=
     '                        COEFFSEXT(NQGP_PIVOT(nzero,nq1))
                            TEMPPLACE(COUNTEXT)=nq1
                            nq2=NXQ(nii,1,NXQ(nij*2,1,
     '                                    NXQ(nik*3,1,nq1,1),1),1)
                            IF((nq2.GT.0).AND.(nq2.NE.nq1)) THEN
                              DO nzero2=1,NQGP(0,nq2)
                                IF(NQGP(nzero2,nq2).EQ.nq) THEN
                                  CALL CALC_GRID_COEF(NENQ,NITB,nq,
     '                              NQGP,NQS,NQXI,NWQ(1,nq),
     '                              NXQ(-NIM,0,0,1),nx_ext,nx_trans,
     '                              COEFFSEXT,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,
     '                              GCHQ(1,nq),GUQ(1,1,nq),NQGW(1,nq),
     '                              PROPQ(1,1,1,1,nq),
     '                              BIDOMAIN,FIXQ,NOTIMP,
     '                              SOLVEEIGHTPROBLEM,ERROR,*9999)
                                  COUNTEXT=COUNTEXT+1
                                  TEMPCOEFF(COUNTEXT)=
     '                              COEFFSEXT(NQGP_PIVOT(nzero2,nq2))
                                  TEMPPLACE(COUNTEXT)=nq2
                                ENDIF !NQGP
                              ENDDO !nzero2
                            ENDIF !nq2
                          ENDIF !NQGP
                        ENDDO !nzero
                      ENDIF !nq1
                    ENDDO !nii
                  ENDDO !nij
                ENDDO !nik

                CALL ISORTP(COUNTEXT,TEMPPLACE,PLACE_PIVOT)

                DO nzero=1,COUNTEXT
                  PLACEEXT=PLACEEXT+1
                  IF(.NOT.UPDATE)
     '              ISR_GKK(PLACEEXT,nx_ext)=TEMPPLACE(nzero)
                  GKK(PLACEEXT,nx_ext)=TEMPCOEFF(PLACE_PIVOT(nzero))
                ENDDO !nzero
                IF(.NOT.UPDATE)
     '            ISC_GKK(nq+1,nx_ext)=ISC_GKK(nq,nx_ext)+COUNTEXT
              ELSE IF(SPARSEGKK(nx_ext).EQ.5) THEN !row/col #2
                IF(USE_LAT.EQ.0) THEN
                  CALL CALC_GRID_COEF(NENQ,NITB,nq,NQGP,NQS,
     '              NQXI,NWQ(1,nq),NXQ(-NIM,0,0,1),nx_ext,nx_trans,
     '              COEFFSEXT,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,
     '              GCHQ(1,nq),GUQ(1,1,nq),
     '              NQGW(1,nq),PROPQ(1,1,1,1,nq),
     '              BIDOMAIN,FIXQ,NOTIMP,SOLVEEIGHTPROBLEM,ERROR,*9999)
                  DO nzero=1,NQGP(0,nq)
                    PLACEEXT=PLACEEXT+1
                    IF(.NOT.UPDATE) ISC_GKK(PLACEEXT,nx_ext)=nq
                    GKK(PLACEEXT,nx_ext)=COEFFSEXT(NQGP_PIVOT(nzero,nq))
                  ENDDO !nzero
                ELSE !use lattice grid scheme
                  IF(.NOT.UPDATE) THEN
                    CALL GSUPPORT(NLATNE,NLATNQ,NLATPNQ,nq,NQNLAT,NQS,
     &                NQXI,NWQ,NQGP(0,nq),ERROR,*9999)
                    CALL ASSERT(PLACEEXT+NQGP(0,nq).LE.NISR_GKKM,
     &                '>> Increase NISR_GKKM',ERROR,*9999)
                    CALL ASSERT(PLACEEXT+NQGP(0,nq).LE.NISC_GKKM,
     &                '>> Increase NISC_GKKM',ERROR,*9999)
                    DO nzero=1,NQGP(0,nq)
                      PLACEEXT=PLACEEXT+1
                      ISR_GKK(PLACEEXT,nx_ext)=nq
                      ISC_GKK(PLACEEXT,nx_ext)=NQGP(nzero,nq)
                    ENDDO
                  ENDIF
                  CALL CALC_LATTICE_WEIGHTS(nq,NQGP(0,nq),M,XQ,
     &              ERROR,*9999)
                  CALL CALC_LATT_COEF_EXT(NITB,nq,NWQ(1,nq),nx_ext,
     &              NQGP(0,nq),LATTCOEFS,AQ(1,nq),M,PROPQ,FIXQ,
     &              SOLVEEIGHTPROBLEM,ERROR,*9999)
                  CALL ASSERT(PLACEEXT.LE.NZ_GKK_M,
     &              '>>Increase NZ_GKK_M',ERROR,*9999)
                  DO i=1,NQGP(0,nq)
                    GKK(PLACEEXT-NQGP(0,nq)+i,nx_ext)=LATTCOEFS(i)
                  ENDDO
                ENDIF !USE_LAT
              ENDIF
            ENDDO !nq
          ENDIF

          IF(COUPBID) THEN
            DO nq=NQR(1,nr),NQR(2,nr) !create line for each nq
              IF(SPARSEGKK(nx_upd).EQ.0) THEN !no sparsity
                CALL CALC_GRID_COEF(NENQ,NITB,nq,NQGP,NQS,
     '            NQXI,NWQ(1,nq),NXQ(-NIM,0,0,1),nx_upd,nx_trans,
     '            COEFFSEXT,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,
     '            GCHQ(1,nq),GUQ(1,1,nq),
     '            NQGW(1,nq),PROPQ(1,1,1,1,nq),
     '            .TRUE.,FIXQ,IMPLICIT,SOLVEEIGHTPROBLEM,ERROR,*9999)
                DO nzero=1,NQGP(0,nq)
                  GKK(((NQGP(nzero,nq)-1)*NQT)+nq,nx_upd)=
     '              COEFFSEXT(NQGP_PIVOT(nzero,nq))
                ENDDO
              ELSE IF(SPARSEGKK(nx_upd).EQ.1) THEN !compressed row
                CALL CALC_GRID_COEF(NENQ,NITB,nq,NQGP,NQS,
     '            NQXI,NWQ(1,nq),NXQ(-NIM,0,0,1),nx_upd,nx_trans,
     '            COEFFSEXT,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,
     '            GCHQ(1,nq),GUQ(1,1,nq),
     '            NQGW(1,nq),PROPQ(1,1,1,1,nq),
     '            .TRUE.,FIXQ,IMPLICIT,SOLVEEIGHTPROBLEM,ERROR,*9999)
                COUNTEXT=0
                DO nzero=1,NQGP(0,nq)
                  PLACEEXT=PLACEEXT+1
                  COUNTEXT=COUNTEXT+1
                  IF(.NOT.UPDATE)
     '              ISC_GKK(PLACEEXT,nx_upd)=NQGP(nzero,nq)
                  GKK(PLACEEXT,nx_upd)=COEFFSEXT(NQGP_PIVOT(nzero,nq))
                ENDDO
                IF(.NOT.UPDATE)
     '            ISR_GKK(nq+1,nx_upd)=ISR_GKK(nq,nx_upd)+COUNTEXT
              ELSE IF(SPARSEGKK(nx_upd).EQ.2) THEN !row/col
                CALL CALC_GRID_COEF(NENQ,NITB,nq,NQGP,NQS,
     '            NQXI,NWQ(1,nq),NXQ(-NIM,0,0,1),nx_upd,nx_trans,
     '            COEFFSEXT,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,
     '            GCHQ(1,nq),GUQ(1,1,nq),
     '            NQGW(1,nq),PROPQ(1,1,1,1,nq),
     '            .TRUE.,FIXQ,IMPLICIT,SOLVEEIGHTPROBLEM,ERROR,*9999)
                DO nzero=1,NQGP(0,nq)
                  PLACEEXT=PLACEEXT+1
                  IF(.NOT.UPDATE) THEN
                    ISR_GKK(PLACEEXT,nx_upd)=nq
                    ISC_GKK(PLACEEXT,nx_upd)=NQGP(nzero,nq)
                  ENDIF
                  GKK(PLACEEXT,nx_upd)=COEFFSEXT(NQGP_PIVOT(nzero,nq))
                ENDDO
              ELSE IF(SPARSEGKK(nx_upd).EQ.3) THEN !compressed column
                COUNTEXT=0
                IK=MAX(0,NITB-2) !zero for 1,2-D, one for 3-D
                IJ=MIN(NITB-1,1) !zero for 1-D, one for 2,3-D

                DO nik=-IK,IK
                  DO nij=-IJ,IJ
                    DO nii=-1,1
                      nq1=NXQ(nii,1,NXQ(nij*2,1,NXQ(nik*3,1,nq,1),1),1)
                      IF(nq1.GT.0) THEN
                        DO nzero=1,NQGP(0,nq1)
                          IF(NQGP(nzero,nq1).EQ.nq) THEN
                            CALL CALC_GRID_COEF(NENQ,NITB,nq,NQGP,
     '                        NQS,NQXI,NWQ(1,nq),
     '                        NXQ(-NIM,0,0,1),nx_upd,nx_trans,
     '                        COEFFSEXT,CQ,DNUDXQ,
     '                        DXDXIQ,DXDXIQ2,
     '                        GCHQ(1,nq),GUQ(1,1,nq),NQGW(1,nq),
     '                        PROPQ(1,1,1,1,nq),.TRUE.,FIXQ,
     '                        IMPLICIT,SOLVEEIGHTPROBLEM,ERROR,*9999)
                            COUNTEXT=COUNTEXT+1
                            TEMPCOEFF(COUNTEXT)=
     '                        COEFFSEXT(NQGP_PIVOT(nzero,nq1))
                            TEMPPLACE(COUNTEXT)=nq1
                            nq2=NXQ(nii,1,NXQ(nij*2,1,
     '                                    NXQ(nik*3,1,nq1,1),1),1)
                            IF((nq2.GT.0).AND.(nq2.NE.nq1)) THEN
                              DO nzero2=1,NQGP(0,nq2)
                                IF(NQGP(nzero2,nq2).EQ.nq) THEN
                                  CALL CALC_GRID_COEF(NENQ,NITB,nq,
     '                              NQGP,NQS,NQXI,NWQ(1,nq),
     '                              NXQ(-NIM,0,0,1),nx_upd,nx_trans,
     '                              COEFFSEXT,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,
     '                              GCHQ(1,nq),GUQ(1,1,nq),
     '                              NQGW(1,nq),PROPQ(1,1,1,1,nq),
     '                              .TRUE.,FIXQ,IMPLICIT,
     '                              SOLVEEIGHTPROBLEM,ERROR,*9999)
                                  COUNTEXT=COUNTEXT+1
                                  TEMPCOEFF(COUNTEXT)=
     '                              COEFFSEXT(NQGP_PIVOT(nzero2,nq2))
                                  TEMPPLACE(COUNTEXT)=nq2
                                ENDIF
                              ENDDO
                            ENDIF
                          ENDIF
                        ENDDO
                      ENDIF
                    ENDDO !nii
                  ENDDO !nij
                ENDDO !nik

                CALL ISORTP(COUNTEXT,TEMPPLACE,PLACE_PIVOT)

                DO nzero=1,COUNTEXT
                  PLACEEXT=PLACEEXT+1
                  IF(.NOT.UPDATE)
     '              ISR_GKK(PLACEEXT,nx_upd)=TEMPPLACE(nzero)
                  GKK(PLACEEXT,nx_upd)=TEMPCOEFF(PLACE_PIVOT(nzero))
                ENDDO !nzero
                IF(.NOT.UPDATE)
     '            ISC_GKK(nq+1,nx_upd)=ISC_GKK(nq,nx_upd)+COUNTEXT
              ELSE IF(SPARSEGKK(nx_upd).EQ.5) THEN !row/col #2
                CALL CALC_GRID_COEF(NENQ,NITB,nq,NQGP,NQS,
     '            NQXI,NWQ(1,nq),NXQ(-NIM,0,0,1),nx_upd,nx_trans,
     '            COEFFSEXT,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,
     '            GCHQ(1,nq),GUQ(1,1,nq),
     '            NQGW(1,nq),PROPQ(1,1,1,1,nq),
     '            .TRUE.,FIXQ,IMPLICIT,SOLVEEIGHTPROBLEM,ERROR,*9999)
                DO nzero=1,NQGP(0,nq)
                  PLACEEXT=PLACEEXT+1
                  IF(.NOT.UPDATE) ISC_GKK(PLACEEXT,nx_upd)=nq
                  GKK(PLACEEXT,nx_upd)=COEFFSEXT(NQGP_PIVOT(nzero,nq))
                ENDDO !nzero
              ENDIF
            ENDDO !nq
          ENDIF
        ENDDO !nr


        IF(.NOT.UPDATE) THEN
          IF(IMPLICIT) THEN
            NZZT(1,NRLIST(1),nx_trans)=PLACEINT
            IF(SPARSEGKK(nx_trans).EQ.5) THEN
              COUNTINT=1
              DO nrr=1,NRLIST(0)
                nr=NRLIST(nrr)
                DO nq=NQR(1,nr),NQR(2,nr)
                  DO nzero=1,NQGP(0,nq)
                    ISC_GKK(PLACEINT+COUNTINT,nx_trans)=NQGP(nzero,nq)
                    COUNTINT=COUNTINT+1
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDIF
          IF(BIDOMAIN) THEN
            NZZT(1,NRLIST(1),nx_ext)=PLACEEXT
            IF(SPARSEGKK(nx_ext).EQ.5) THEN
              IF(USE_LAT.EQ.0) THEN                
                COUNTEXT=1
                DO nrr=1,NRLIST(0)
                  nr=NRLIST(nrr)
                  DO nq=NQR(1,nr),NQR(2,nr)
                    DO nzero=1,NQGP(0,nq)
                      ISC_GKK(PLACEEXT+COUNTEXT,nx_ext)=NQGP(nzero,nq)
                      COUNTEXT=COUNTEXT+1
                    ENDDO
                  ENDDO
                ENDDO
              ELSE !lattice method
C               The row/col sparsity pattern (2) is rearranged into
C               the row/col2 pattern (5).
                DO COUNTEXT=1,PLACEEXT
                  ISC_GKK(PLACEEXT+COUNTEXT,nx_ext)=
     '              ISC_GKK(COUNTEXT,nx_ext)
                  ISC_GKK(COUNTEXT,nx_ext)=ISR_GKK(COUNTEXT,nx_ext)
                ENDDO
              ENDIF
            ENDIF
          ENDIF
          IF(COUPBID) THEN
            NZZT(1,NRLIST(1),nx_upd)=PLACEEXT
            IF(SPARSEGKK(nx_upd).EQ.5) THEN
              COUNTEXT=1
              DO nrr=1,NRLIST(0)
                nr=NRLIST(nrr)
                DO nq=NQR(1,nr),NQR(2,nr)
                  DO nzero=1,NQGP(0,nq)
                    ISC_GKK(PLACEEXT+COUNTEXT,nx_upd)=NQGP(nzero,nq)
                    COUNTEXT=COUNTEXT+1
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDIF
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
     '    '(/'' Global stiffness matrix GKK - transmembrane:'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' NOT(1,1,nr,nx)='',I5,'
     '    //''', NOT(2,1,nr,nx)='',I5)') NOT(1,1,nr,nx_trans),
     '    NOT(2,1,nr,nx_trans)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        CALL OPSTFMAT(DUMMY_LIST,ISC_GKK,ISR_GKK,IOOP,
     '    NOT(1,1,nr,nx_trans),NOT(2,1,nr,nx_trans),
     '    NZZT(1,nr,nx_trans),DUMMY_LIST,SPARSEGKK(nx_trans),
     '    GKK,GKK,'GKK','GKK',.TRUE.,.TRUE.,.FALSE.,
     '    ERROR,*9999)
        IF(BIDOMAIN) THEN
          WRITE(OP_STRING,
     '      '(/'' Global stiffness matrix GKK - external:'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NOT(1,1,nr,nx)='',I5,'
     '      //''', NOT(2,1,nr,nx)='',I5)') NOT(1,1,nr,nx_ext),
     '      NOT(2,1,nr,nx_ext)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CALL OPSTFMAT(DUMMY_LIST,ISC_GKK,ISR_GKK,IOOP,
     '      NOT(1,1,nr,nx_ext),NOT(2,1,nr,nx_ext),
     '      NZZT(1,nr,nx_ext),DUMMY_LIST,SPARSEGKK(nx_ext),
     '      GKK,GKK,'GKK','GKK',.TRUE.,.TRUE.,.FALSE.,
     '      ERROR,*9999)
        ENDIF
      ENDIF

      
      IF(IWRIT5(NRLIST(1),nx_trans).GE.1) THEN
        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)
        WRITE(OP_STRING,'(1X,''Time for matrix assembly '',F8.2,'
     '    //'''s cpu'')') ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF


      CALL EXITS('ASSEMBLE10')
      RETURN
 9999 CALL ERRORS('ASSEMBLE10',ERROR)
      CALL EXITS('ASSEMBLE10')
      RETURN 1
      END


      SUBROUTINE ASSEMBLE10_FE(ISC_GKK,ISR_GKK,NBH,NEELEM,NQGP,
     '  NQGP_PIVOT,NQS,NQXI,NRLIST,nx_ext,nx_trans,nx_upd,NXQ,CQ,
     '  GKK,GM,NQGW,PG,WG,XQ,PROPQ,BIDOMAIN,COUPBID,FIRST_A,FIXQ,
     '  IMPLICIT,UPDATE,UPDATEDT,UPDATE_MATRIX,SOLVEEIGHTPROBLEM,
     '  ERROR,*)

C#### Subroutine: ASSEMBLE10_FE
C###  Description:
C###    ASSEMBLE10_FE writes directly into the compressed row storage
C###    arrays. It generates the constraint matrix for the implicit
C###    solution of activation problems, using the Grid-based Finite
C###    Element method.
C***  Created by Scott Marsden November 2000

      IMPLICIT NONE

      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b12.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:grid00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:iwrit00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
      INCLUDE 'cmiss$reference:solv00.cmn'
      INCLUDE 'cmiss$reference:time02.cmn'
      INCLUDE 'cmiss$reference:tol00.cmn'

!     Parameter list
      INTEGER ISC_GKK(NISC_GKKM,NXM),ISR_GKK(NISR_GKKM,NXM),
     '  NBH(NHM,NCM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NQGP(0:NQGM,NQM),NQGP_PIVOT(NQGM,NQM),
     '  NQS(NEQM),NQXI(0:NIM,NQSCM),NRLIST(0:NRM),nx_ext,
     '  nx_trans,nx_upd,NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      REAL*8 CQ(NMM,NQM),GKK(NZ_GKK_M,NXM),GM(NZ_GM_M),
     '  NQGW(NQGM,NQM),PG(NSM,NUM,NGM,NBM),WG(NGM,NBM),XQ(NJM,NQM),
     '  PROPQ(3,3,4,2,NQM)
      CHARACTER ERROR*(*)
      LOGICAL BIDOMAIN,COUPBID,FIRST_A,FIXQ(NYQM,NIYFIXM,NXM),
     '  IMPLICIT,SOLVEEIGHTPROBLEM,UPDATE,UPDATEDT,UPDATE_MATRIX
!     Local variables
      INTEGER COUNTINT,COUNTEXT,DUMMY_LIST(0:1),IEND,Inc_NZ_GKK_M,
     '  Inc_NISR_GKK_M,Inc_NISC_GKK_M,maxrow,maxcol,nb,NITB,
     '  nq,nr,nrr,nx,nzero,nzz,nzz_trans,nzz_ext,nzz_upd,
     '  PLACEINT,PLACEEXT
      REAL*8 COEFFSEXT(NQGM),DTTHETA
      REAL ELAPSED_TIME,TIME_START(1),TIME_STOP(1)
      LOGICAL ERROR_FLAG
C SGM 29 November 2000 not used
C      LOGICAL NOTIMP

      CALL ENTERS('ASSEMBLE10_FE',*9999)

      ERROR_FLAG=.FALSE.

      CALL CPU_TIMER(CPU_USER,TIME_START)

      IF(nx_trans.NE.0) THEN
        nx=nx_trans
      ELSE
        nx=nx_ext
      ENDIF

C SGM 17Jan01 check array sizes for cases which don't use NQGP
C     and initialise arrays.
C     (NQGP is calculated in CALC_FE_GRID_COEF later on).
      IF(.NOT.UPDATE) THEN
        IF(IMPLICIT) THEN
          NOT(1,1,NRLIST(1),nx_trans)=NQT
          NOT(2,1,NRLIST(1),nx_trans)=NQT
          IF(SPARSEGKK(nx_trans).EQ.0) THEN !no sparsity
            nzz=NQT*NQT
            NZZT(1,NRLIST(1),nx_trans)=nzz
            IF(NZ_GKK_M.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ENDIF
CC$OMP PARALLEL DO
CC$OMP&  PRIVATE(nq),
CC$OMP&  SHARED(GKK,nx_trans)
            !initialise
            DO nq=1,NQT*NQT
              GKK(nq,nx_trans)=0.0d0
            ENDDO !nq
CC$OMP END PARALLEL DO
          ELSE IF(SPARSEGKK(nx_trans).EQ.1) THEN !compressed row
            maxrow=0
            DO nrr=1,NRLIST(0)
              nr=NRLIST(nrr)
              IF(NQR(2,nr)+1.GT.maxrow) maxrow=NQR(2,nr)+1
            ENDDO !nr
            IF(NISR_GKKM.LT.maxrow) THEN
              WRITE(ERROR,'(''>>Increase ISR_GKK_M to >= '',I12)')
     '          maxrow
              GOTO 9999
            ENDIF
            !initialise
            ISR_GKK(1,nx_trans)=1
          ELSE IF(SPARSEGKK(nx_trans).EQ.2) THEN !row/col
            !nothing to do
          ELSE IF(SPARSEGKK(nx_trans).EQ.3) THEN !compressed column
            maxcol=0
            DO nrr=1,NRLIST(0)
              nr=NRLIST(nrr)
              IF(NQR(2,nr)+1.GT.maxcol) maxcol=NQR(2,nr)+1
            ENDDO !nr
            IF(NISC_GKKM.LT.maxcol) THEN
              WRITE(ERROR,'(''>>Increase ISC_GKK_M to >= '',I12)')
     '          maxcol
              GOTO 9999
            ENDIF
            !initialise
            ISC_GKK(1,nx_trans)=1
          ELSE IF(SPARSEGKK(nx_trans).EQ.4) THEN
            ERROR='>>Not implemented'
            GOTO 9999
          ELSE IF(SPARSEGKK(nx_trans).EQ.5) THEN !row/col #2
            !nothing to do
          ELSE
            ERROR='>>Unknown sparsity type for GKK'
            GOTO 9999
          ENDIF
        ENDIF !implicit
        IF(BIDOMAIN) THEN
          NOT(1,1,NRLIST(1),nx_ext)=NQT
          NOT(2,1,NRLIST(1),nx_ext)=NQT
          IF(SPARSEGKK(nx_ext).EQ.0) THEN !no sparsity
            nzz=NQT*NQT
            NZZT(1,NRLIST(1),nx_ext)=nzz
            IF(NZ_GKK_M.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ENDIF
CC$OMP PARALLEL DO
CC$OMP&  PRIVATE(nq),
CC$OMP&  SHARED(GKK,nx_ext)
            !initialise
            DO nq=1,NQT*NQT
              GKK(nq,nx_ext)=0.0d0
            ENDDO !nq
CC$OMP END PARALLEL DO
          ELSE IF(SPARSEGKK(nx_ext).EQ.1) THEN !compressed row
            maxrow=0
            DO nrr=1,NRLIST(0)
              nr=NRLIST(nrr)
              IF(NQR(2,nr)+1.GT.maxrow) maxrow=NQR(2,nr)+1
            ENDDO !nr
            IF(NISR_GKKM.LT.maxrow) THEN
              WRITE(ERROR,'(''>>Increase ISR_GKK_M to >= '',I12)')
     '          maxrow
              GOTO 9999
            ENDIF
            !initialise
            ISR_GKK(1,nx_ext)=1
          ELSE IF(SPARSEGKK(nx_ext).EQ.2) THEN !row/col
            !nothing to do
          ELSE IF(SPARSEGKK(nx_ext).EQ.3) THEN !compressed column
            maxcol=0
            DO nrr=1,NRLIST(0)
              nr=NRLIST(nrr)
              IF(NQR(2,nr)+1.GT.maxcol) maxcol=NQR(2,nr)+1
            ENDDO !nr
            IF(NISC_GKKM.LT.maxcol) THEN
              WRITE(ERROR,'(''>>Increase ISC_GKK_M to >= '',I12)')
     '          maxcol
              GOTO 9999
            ENDIF
            !initialise
            ISC_GKK(1,nx_ext)=1
          ELSE IF(SPARSEGKK(nx_ext).EQ.4) THEN
            ERROR='>>Not implemented'
            GOTO 9999
          ELSE IF(SPARSEGKK(nx_ext).EQ.5) THEN !row/col #2
            !nothing to do
          ELSE
            ERROR='>>Unknown sparsity type for GKK'
            GOTO 9999
          ENDIF
        ENDIF !bidomain
        IF(COUPBID) THEN
          NOT(1,1,NRLIST(1),nx_upd)=NQT
          NOT(2,1,NRLIST(1),nx_upd)=NQT
          IF(SPARSEGKK(nx_upd).EQ.0) THEN !no sparsity
            nzz=NQT*NQT
            NZZT(1,NRLIST(1),nx_upd)=nzz
            IF(NZ_GKK_M.LT.nzz) THEN
              WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
              GOTO 9999
            ENDIF
CC$OMP PARALLEL DO
CC$OMP&  PRIVATE(nq),
CC$OMP&  SHARED(GKK,nx_upd)
            !initialise
            DO nq=1,NQT*NQT
              GKK(nq,nx_upd)=0.0d0
            ENDDO !nq
CC$OMP END PARALLEL DO
          ELSE IF(SPARSEGKK(nx_upd).EQ.1) THEN !compressed row
            maxrow=0
            DO nrr=1,NRLIST(0)
              nr=NRLIST(nrr)
              IF(NQR(2,nr)+1.GT.maxrow) maxrow=NQR(2,nr)+1
            ENDDO !nr
            IF(NISR_GKKM.LT.maxrow) THEN
              WRITE(ERROR,'(''>>Increase ISR_GKK_M to >= '',I12)')
     '          maxrow
              GOTO 9999
            ENDIF
            !initialise
            ISR_GKK(1,nx_upd)=1
          ELSE IF(SPARSEGKK(nx_upd).EQ.2) THEN !row/col
            !nothing to do
          ELSE IF(SPARSEGKK(nx_upd).EQ.3) THEN !compressed column
            maxcol=0
            DO nrr=1,NRLIST(0)
              nr=NRLIST(nrr)
              IF(NQR(2,nr)+1.GT.maxcol) maxcol=NQR(2,nr)+1
            ENDDO !nr
            IF(NISC_GKKM.LT.maxcol) THEN
              WRITE(ERROR,'(''>>Increase ISC_GKK_M to >= '',I12)')
     '          maxcol
              GOTO 9999
            ENDIF
            !initialise
            ISC_GKK(1,nx_upd)=1
          ELSE IF(SPARSEGKK(nx_upd).EQ.4) THEN
            ERROR='>>Not implemented'
            GOTO 9999
          ELSE IF(SPARSEGKK(nx_upd).EQ.5) THEN !row/col #2
            !nothing to do
          ELSE
            ERROR='>>Unknown sparsity type for GKK'
            GOTO 9999
          ENDIF
        ENDIF !coupid
      ENDIF !not update


      PLACEINT=0
      PLACEEXT=0

      nzz_trans=0 ! nzz counter for intracellular non-zeros
      nzz_ext=0   ! nzz counter for extracellular non-zeros
      nzz_upd=0   ! COUPID

      Inc_NZ_GKK_M=0   ! value to increace NZ_GKK_M to (if any)
      Inc_NISC_GKK_M=0 ! value to increace NISC_GKK_M to (if any)
      Inc_NISR_GKK_M=0 ! value to increace NISR_GKK_M to (if any)
C nr and nq loops now outside conditionals.
      DO nrr=1,NRLIST(0)

        nr=NRLIST(nrr)
        nb=NBH(NH_LOC(1,nx),1,NEELEM(1,nr))
        NITB=NQXI(0,NQS(NEELEM(1,nr)))

C make sure the number of element  Xi coordinates equals the number of
C       Xi coordiates for the dependent variable before call to
C     CALC_FE_GRID_COEF.
        CALL ASSERT(NITB.EQ.NIT(nb),
     '    '>>#Xi-coords inconsistent between element and dependent var',
     '    ERROR,*9999)

C SGM 01Feb01 Moved following two checks inside nr loop.
        !Check array sizes
C MLT 30/11/02 Modified to account for different support for 
C grid FE and grid FV
        IF(ITYP4(nr,nx).EQ.6.AND.NQGM.LT.3**NITB) THEN ! Grid FE
          IEND=0
          CALL APPENDC(IEND,'For Grid FE NQGM needs to be at least ',
     '                 ERROR)
          CALL APPENDI(IEND,3**NITB,ERROR)
          GOTO 9999 
        ELSE IF(ITYP4(nr,nx).EQ.7.AND.NQGM.LT.2*NITB+1) THEN ! Grid FV
          IEND=0
          CALL APPENDC(IEND,'For Grid FV NQGM needs to be at least ',
     '                 ERROR)
          CALL APPENDI(IEND,2*NITB+1,ERROR)
          GOTO 9999 
        ENDIF

        IF(NZ_GM_M.LT.NQT*NQGM) THEN
          IEND=0
          CALL APPENDC(IEND,'NZ_GM_M needs to be at least ',ERROR)
          CALL APPENDI(IEND,NQT*NQGM,ERROR)
          GOTO 9999
        ENDIF

        DO nq=NQR(1,nr),NQR(2,nr) !create line for each nq

C MLT 29Nov02 Adding grid finite volumes
          IF(ITYP4(nr,nx).EQ.6) THEN ! Grid-based FE
            CALL CALC_FE_GRID_COEF(nb,NITB,nq,NQGP,NQGP_PIVOT,nr,
     '        nx_ext,NXQ(-NIM,0,0,1),COEFFSEXT,CQ,GM,NQGW(1,nq),PG,WG,
     '        XQ,BIDOMAIN,FIXQ,SOLVEEIGHTPROBLEM,ERROR,*9999)
          ELSE ! Grid FV
            CALL CALC_FV_GRID_COEF(NITB,nq,NQGP,NQGP_PIVOT,nr,
     '        nx_ext,NXQ(-NIM,0,0,1),COEFFSEXT,CQ,GM,NQGW(1,nq),
     '        XQ,PROPQ,BIDOMAIN,FIXQ,SOLVEEIGHTPROBLEM,ERROR,*9999)
          ENDIF

C Check array sizes (NZ_GKK_M,ISR_GKK,ISC_GKK)
          IF(.NOT.UPDATE) THEN
            IF(IMPLICIT) THEN
              IF(SPARSEGKK(nx_trans).EQ.0) THEN !no sparsity
                !nothing to do
              ELSE IF(SPARSEGKK(nx_trans).EQ.1) THEN !compressed row
                nzz_trans = nzz_trans + NQGP(0,nq)
                NZZT(1,nr,nx_trans)=nzz_trans
                IF(NZ_GKK_M.LT.nzz_trans) THEN
                  IF(Inc_NZ_GKK_M.LT.nzz_trans)Inc_NZ_GKK_M = nzz_trans
                ELSE IF(NISC_GKKM.LT.nzz_trans) THEN
                  IF(Inc_NISC_GKK_M.LT.nzz_trans)Inc_NISC_GKK_M =
     '              nzz_trans
                ENDIF
              ELSE IF(SPARSEGKK(nx_trans).EQ.2) THEN !row/col
                nzz_trans = nzz_trans + NQGP(0,nq)
                NZZT(1,nr,nx_trans)=nzz_trans
                IF(NZ_GKK_M.LT.nzz_trans) THEN
                  IF(Inc_NZ_GKK_M.LT.nzz_trans)Inc_NZ_GKK_M = nzz_trans
                ELSE IF(NISR_GKKM.LT.nzz_trans) THEN
                  IF(Inc_NISR_GKK_M.LT.nzz_trans)Inc_NISR_GKK_M =
     '              nzz_trans
                ELSE IF(NISC_GKKM.LT.nzz_trans) THEN
                  IF(Inc_NISC_GKK_M.LT.nzz_trans)Inc_NISC_GKK_M =
     '              nzz_trans
                ENDIF
              ELSE IF(SPARSEGKK(nx_trans).EQ.3) THEN !compressed column
                nzz_trans = nzz_trans + NQGP(0,nq)
                NZZT(1,nr,nx_trans)=nzz_trans
                IF(NZ_GKK_M.LT.nzz_trans) THEN
                  IF(Inc_NZ_GKK_M.LT.nzz_trans)Inc_NZ_GKK_M = nzz_trans
                ELSE IF(NISR_GKKM.LT.nzz_trans) THEN
                  IF(Inc_NISR_GKK_M.LT.nzz_trans)Inc_NISR_GKK_M =
     '              nzz_trans
                ENDIF
              ELSE IF(SPARSEGKK(nx_trans).EQ.4) THEN
                ERROR='>>Not implemented'
                GOTO 9999
              ELSE IF(SPARSEGKK(nx_trans).EQ.5) THEN !row/col #2
                nzz_trans = nzz_trans + NQGP(0,nq)
                NZZT(1,nr,nx_trans)=nzz_trans
                IF(NZ_GKK_M.LT.nzz_trans) THEN
                  IF(Inc_NZ_GKK_M.LT.nzz_trans)Inc_NZ_GKK_M = nzz_trans
                ELSE IF(NISC_GKKM.LT.2*nzz_trans) THEN
                  IF(Inc_NISC_GKK_M.LT.2*nzz_trans)
     '              Inc_NISC_GKK_M = 2*nzz_trans
C *** DPN 31 July 2001 - Not sure exactly what this GOTO is doing here
C                        but should at least write out an error message
                  ERROR='>>Not implemented?'
                  GOTO 9999
                ENDIF
              ELSE
                ERROR='>>Unknown sparsity type for GKK'
                GOTO 9999
              ENDIF
            ENDIF !implicit
            IF(BIDOMAIN) THEN
              IF(SPARSEGKK(nx_ext).EQ.0) THEN !no sparsity
                !nothin to do
              ELSE IF(SPARSEGKK(nx_ext).EQ.1) THEN !compressed row
                nzz_ext = nzz_ext + NQGP(0,nq)
                NZZT(1,nr,nx_ext)=nzz_ext
                IF(NZ_GKK_M.LT.nzz_ext) THEN
                  IF(Inc_NZ_GKK_M.LT.nzz_ext)Inc_NZ_GKK_M = nzz_ext
                ELSE IF(NISC_GKKM.LT.nzz_ext) THEN
                  IF(Inc_NISC_GKK_M.LT.nzz_ext)Inc_NISC_GKK_M = nzz_ext
                ENDIF
              ELSE IF(SPARSEGKK(nx_ext).EQ.2) THEN !row/col
                nzz_ext = nzz_ext + NQGP(0,nq)
                NZZT(1,nr,nx_ext)=nzz_ext
                IF(NZ_GKK_M.LT.nzz_ext) THEN
                  IF(Inc_NZ_GKK_M.LT.nzz_ext)Inc_NZ_GKK_M = nzz_ext
                ELSE IF(NISR_GKKM.LT.nzz_ext) THEN
                  IF(Inc_NISR_GKK_M.LT.nzz_ext)Inc_NISR_GKK_M = nzz_ext
                ELSE IF(NISC_GKKM.LT.nzz_ext) THEN
                  IF(Inc_NISC_GKK_M.LT.nzz_ext)Inc_NISC_GKK_M = nzz_ext
                ENDIF
              ELSE IF(SPARSEGKK(nx_ext).EQ.3) THEN !compressed column
                nzz_ext = nzz_ext + NQGP(0,nq)
                NZZT(1,nr,nx_ext)=nzz_ext
                IF(NZ_GKK_M.LT.nzz_ext) THEN
                  IF(Inc_NZ_GKK_M.LT.nzz_ext)Inc_NZ_GKK_M = nzz_ext
                ELSE IF(NISR_GKKM.LT.nzz_ext) THEN
                  IF(Inc_NISR_GKK_M.LT.nzz_ext)Inc_NISR_GKK_M = nzz_ext
                ENDIF
              ELSE IF(SPARSEGKK(nx_ext).EQ.4) THEN
                ERROR='>>Not implemented'
                GOTO 9999
              ELSE IF(SPARSEGKK(nx_ext).EQ.5) THEN !row/col #2
                nzz_ext = nzz_ext + NQGP(0,nq)
                NZZT(1,nr,nx_ext)=nzz_ext
                IF(NZ_GKK_M.LT.nzz_ext) THEN
                  IF(Inc_NZ_GKK_M.LT.nzz_ext)Inc_NZ_GKK_M = nzz_ext
                ELSE IF(NISC_GKKM.LT.2*nzz_ext) THEN
                  IF(Inc_NISC_GKK_M.LT.2*nzz_ext)Inc_NISC_GKK_M = 2
     '              *nzz_ext
                ENDIF
              ELSE
                ERROR='>>Unknown sparsity type for GKK'
                GOTO 9999
              ENDIF
            ENDIF !bidomain
            IF(COUPBID) THEN
              IF(SPARSEGKK(nx_upd).EQ.0) THEN !no sparsity
                !nothing to do
              ELSE IF(SPARSEGKK(nx_upd).EQ.1) THEN !compressed row
                nzz_upd = nzz_upd + NQGP(0,nq)
                NZZT(1,nr,nx_upd)=nzz_upd
                IF(NZ_GKK_M.LT.nzz_upd) THEN
                  IF(Inc_NZ_GKK_M.LT.nzz_upd)Inc_NZ_GKK_M = nzz_upd
                ELSE IF(NISC_GKKM.LT.nzz_upd) THEN
                  IF(Inc_NISC_GKK_M.LT.nzz_upd)Inc_NISC_GKK_M = nzz_upd
                ENDIF
              ELSE IF(SPARSEGKK(nx_upd).EQ.2) THEN !row/col
                nzz_upd = nzz_upd + NQGP(0,nq)
                NZZT(1,nr,nx_upd)=nzz_upd
                IF(NZ_GKK_M.LT.nzz_upd) THEN
                  IF(Inc_NZ_GKK_M.LT.nzz_upd)Inc_NZ_GKK_M = nzz_upd
                ELSE IF(NISR_GKKM.LT.nzz_upd) THEN
                  IF(Inc_NISR_GKK_M.LT.nzz_upd)Inc_NISR_GKK_M = nzz_upd
                ELSE IF(NISC_GKKM.LT.nzz_upd) THEN
                  IF(Inc_NISC_GKK_M.LT.nzz_upd)Inc_NISC_GKK_M = nzz_upd
                ENDIF
              ELSE IF(SPARSEGKK(nx_upd).EQ.3) THEN !compressed column
                nzz_upd = nzz_upd + NQGP(0,nq)
                NZZT(1,nr,nx_upd)=nzz_upd
                IF(NZ_GKK_M.LT.nzz_upd) THEN
                  IF(Inc_NZ_GKK_M.LT.nzz_upd)Inc_NZ_GKK_M = nzz_upd
                ELSE IF(NISR_GKKM.LT.nzz_upd) THEN
                  IF(Inc_NISR_GKK_M.LT.nzz_upd)Inc_NISR_GKK_M = nzz_upd
                ENDIF
              ELSE IF(SPARSEGKK(nx_upd).EQ.4) THEN
                ERROR='>>Not implemented'
                GOTO 9999
              ELSE IF(SPARSEGKK(nx_upd).EQ.5) THEN !row/col #2
                nzz_upd = nzz_upd + NQGP(0,nq)
                NZZT(1,nr,nx_upd)=nzz_upd
                IF(NZ_GKK_M.LT.nzz_upd) THEN
                  IF(Inc_NZ_GKK_M.LT.nzz_upd)Inc_NZ_GKK_M = nzz_upd
                ELSE IF(NISC_GKKM.LT.2*nzz_upd) THEN
                  IF(Inc_NISC_GKK_M.LT.2*nzz_upd)Inc_NISC_GKK_M = 2
     '              *nzz_upd
                ENDIF
              ELSE
                ERROR='>>Unknown sparsity type for GKK'
                GOTO 9999
              ENDIF
            ENDIF !coupid
          ENDIF !not update

          IF(Inc_NZ_GKK_M.EQ.0.AND.Inc_NISR_GKK_M.EQ.0
     '      .AND.Inc_NISC_GKK_M.EQ.0) THEN !don't need to increace array sizes
            IF(UPDATEDT) THEN
              IF(IMPLICIT) THEN
                DTTHETA=DT*THETA(1)
                IF(SPARSEGKK(nx_trans).EQ.0) THEN !no sparsity
                  DO nzero=1,NQGP(0,nq)
                    GKK(((NQGP(nzero,nq)-1)*NQT)+nq,nx_trans)=
     '                -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTTHETA +
     '                GM(NQGP_PIVOT(nzero,nq)+NQGM*(nq-1))
                  ENDDO
                ELSE IF(SPARSEGKK(nx_trans).EQ.1) THEN !compressed row
                  PLACEINT=ISR_GKK(nq,nx_trans)
                  DO nzero=1,NQGP(0,nq)
                    GKK(PLACEINT,nx_trans)=
     '                -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTTHETA +
     '                GM(NQGP_PIVOT(nzero,nq)+NQGM*(nq-1))
                    PLACEINT=PLACEINT+1
                  ENDDO !nzero
                ELSE IF(SPARSEGKK(nx_trans).EQ.2) THEN !row/col
                  PLACEINT=ISC_GKK(nq,nx_trans)
                  DO nzero=1,NQGP(0,nq)
                    GKK(PLACEINT,nx_trans)=
     '                -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTTHETA +
     '                GM(NQGP_PIVOT(nzero,nq)+NQGM*(nq-1))
                    PLACEINT=PLACEINT+1
                  ENDDO !nzero
                ELSE IF(SPARSEGKK(nx_trans).EQ.3) THEN !compressed column
C no compressed column for grid-based FE yet
                  IF(.NOT.ERROR_FLAG) THEN
                    ERROR_FLAG=.TRUE.
                    CALL FLAG_ERROR(-1,'>>compressed column storage '
     '                //'not yet implemented for the Grid-based '
     '                //'finite element method')
                  ENDIF
                ELSE IF(SPARSEGKK(nx_trans).EQ.5) THEN !row/col #2
                  PLACEINT=ISC_GKK(nq,nx_trans)
                  DO nzero=1,NQGP(0,nq)
                    GKK(PLACEINT,nx_trans)=
     '                -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTTHETA +
     '                GM(NQGP_PIVOT(nzero,nq)+NQGM*(nq-1))
                  ENDDO
                ENDIF
                IF(ERROR_FLAG) GOTO 9998
              ENDIF !implicit
            ELSE !not updatedt
              IF(IMPLICIT) THEN
                DTTHETA=DT*THETA(1)
                IF(SPARSEGKK(nx_trans).EQ.0) THEN !no sparsity
                  DO nzero=1,NQGP(0,nq)
                    GKK(((NQGP(nzero,nq)-1)*NQT)+nq,nx_trans)=
     '                -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTTHETA +
     '                GM(NQGP_PIVOT(nzero,nq)+NQGM*(nq-1))
                  ENDDO
                ELSE IF(SPARSEGKK(nx_trans).EQ.1) THEN !compressed row
                  COUNTINT=0
                  DO nzero=1,NQGP(0,nq)
                    PLACEINT=PLACEINT+1
                    COUNTINT=COUNTINT+1
                    IF(.NOT.UPDATE)
     '                ISC_GKK(PLACEINT,nx_trans)=NQGP(nzero,nq)
                    GKK(PLACEINT,nx_trans)=
     '                -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTTHETA +
     '                GM(NQGP_PIVOT(nzero,nq)+NQGM*(nq-1))
                  ENDDO
                  IF(.NOT.UPDATE)
     '              ISR_GKK(nq+1,nx_trans)=ISR_GKK(nq,nx_trans)+COUNTINT
                ELSE IF(SPARSEGKK(nx_trans).EQ.2) THEN !row/col
                  DO nzero=1,NQGP(0,nq)
                    PLACEINT=PLACEINT+1
                    IF(.NOT.UPDATE) THEN
                      ISR_GKK(PLACEINT,nx_trans)=nq
                      ISC_GKK(PLACEINT,nx_trans)=NQGP(nzero,nq)
                    ENDIF
                    GKK(PLACEINT,nx_trans)=
     '                -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTTHETA +
     '                GM(NQGP_PIVOT(nzero,nq)+NQGM*(nq-1))
                  ENDDO !nzero
                ELSE IF(SPARSEGKK(nx_trans).EQ.3) THEN !compressed column
C no compressed column for grid-based FE yet
                  CALL ASSERT(ITYP4(nr,nx_trans).NE.6,
     '              '>>compressed column '//
     '              'storage not yet implemented for the Grid-based '//
     '              'finite element method',ERROR,*9999)

                ELSE IF(SPARSEGKK(nx_trans).EQ.5) THEN !row/col #2
                  DO nzero=1,NQGP(0,nq)
                    PLACEINT=PLACEINT+1
                    IF(.NOT.UPDATE) ISC_GKK(PLACEINT,nx_trans)=nq
                    GKK(PLACEINT,nx_trans)=
     '                -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTTHETA +
     '                GM(NQGP_PIVOT(nzero,nq)+NQGM*(nq-1))
                  ENDDO
                ENDIF !SPARSEGKK
              ENDIF !implicit

              IF(BIDOMAIN) THEN
                IF(SPARSEGKK(nx_ext).EQ.0) THEN !no sparsity
                  DO nzero=1,NQGP(0,nq)
                    GKK(((NQGP(nzero,nq)-1)*NQT)+nq,nx_ext)=
     '                COEFFSEXT(NQGP_PIVOT(nzero,nq))
                  ENDDO
                ELSE IF(SPARSEGKK(nx_ext).EQ.1) THEN !compressed row
                  COUNTEXT=0
                  DO nzero=1,NQGP(0,nq)
                    PLACEEXT=PLACEEXT+1
                    COUNTEXT=COUNTEXT+1
                    IF(.NOT.UPDATE)
     '                ISC_GKK(PLACEEXT,nx_ext)=NQGP(nzero,nq)
                    GKK(PLACEEXT,nx_ext)=COEFFSEXT(NQGP_PIVOT(nzero,nq))
                  ENDDO
                  IF(.NOT.UPDATE)
     '              ISR_GKK(nq+1,nx_ext)=ISR_GKK(nq,nx_ext)+COUNTEXT
                ELSE IF(SPARSEGKK(nx_ext).EQ.2) THEN !row/col
                  DO nzero=1,NQGP(0,nq)
                    PLACEEXT=PLACEEXT+1
                    IF(.NOT.UPDATE) THEN
                      ISR_GKK(PLACEEXT,nx_ext)=nq
                      ISC_GKK(PLACEEXT,nx_ext)=NQGP(nzero,nq)
                    ENDIF
                    GKK(PLACEEXT,nx_ext)=COEFFSEXT(NQGP_PIVOT(nzero,nq))
                  ENDDO
                ELSE IF(SPARSEGKK(nx_ext).EQ.3) THEN !compressed column

C no compressed column for grid-based FE yet
                  CALL ASSERT(ITYP4(nr,nx_ext).NE.6,
     '              '>>compressed column '//
     '              'storage not yet implemented for the Grid-based '//
     '              'finite element method',ERROR,*9999)

                ELSE IF(SPARSEGKK(nx_ext).EQ.5) THEN !row/col #2
                  DO nzero=1,NQGP(0,nq)
                    PLACEEXT=PLACEEXT+1
                    IF(.NOT.UPDATE) ISC_GKK(PLACEEXT,nx_ext)=nq
                    GKK(PLACEEXT,nx_ext)=COEFFSEXT(NQGP_PIVOT(nzero,nq))
                  ENDDO !nzero
                ENDIF
              ENDIF !bidomain

              IF(COUPBID) THEN
                IF(SPARSEGKK(nx_upd).EQ.0) THEN !no sparsity
                  DO nzero=1,NQGP(0,nq)
                    GKK(((NQGP(nzero,nq)-1)*NQT)+nq,nx_upd)=
     '                COEFFSEXT(NQGP_PIVOT(nzero,nq))
                  ENDDO
                ELSE IF(SPARSEGKK(nx_upd).EQ.1) THEN !compressed row
                  COUNTEXT=0
                  DO nzero=1,NQGP(0,nq)
                    PLACEEXT=PLACEEXT+1
                    COUNTEXT=COUNTEXT+1
                    IF(.NOT.UPDATE)
     '                ISC_GKK(PLACEEXT,nx_upd)=NQGP(nzero,nq)
                    GKK(PLACEEXT,nx_upd)=COEFFSEXT(NQGP_PIVOT(nzero,nq))
                  ENDDO
                  IF(.NOT.UPDATE)
     '              ISR_GKK(nq+1,nx_upd)=ISR_GKK(nq,nx_upd)+COUNTEXT
                ELSE IF(SPARSEGKK(nx_upd).EQ.2) THEN !row/col
                  DO nzero=1,NQGP(0,nq)
                    PLACEEXT=PLACEEXT+1
                    IF(.NOT.UPDATE) THEN
                      ISR_GKK(PLACEEXT,nx_upd)=nq
                      ISC_GKK(PLACEEXT,nx_upd)=NQGP(nzero,nq)
                    ENDIF
                    GKK(PLACEEXT,nx_upd)=COEFFSEXT(NQGP_PIVOT(nzero,nq))
                  ENDDO
                ELSE IF(SPARSEGKK(nx_upd).EQ.3) THEN !compressed column
C no compressed column for grid-based FE yet
                  CALL ASSERT(ITYP4(nr,nx_upd).NE.6,
     '              '>>compressed column '//
     '              'storage not yet implemented for the Grid-based '//
     '              'finite element method',ERROR,*9999)

                ELSE IF(SPARSEGKK(nx_upd).EQ.5) THEN !row/col #2
                  DO nzero=1,NQGP(0,nq)
                    PLACEEXT=PLACEEXT+1
                    IF(.NOT.UPDATE) ISC_GKK(PLACEEXT,nx_upd)=nq
                    GKK(PLACEEXT,nx_upd)=COEFFSEXT(NQGP_PIVOT(nzero,nq))
                  ENDDO !nzero
                ENDIF
              ENDIF !coupid
            ENDIF !updatedt
          ENDIF !don't need to increace array sizes
        ENDDO !nq
      ENDDO !nr

c      write(*,'("GKK entries:")') 
c      write(*,'("nq: 451, transmembrane and extracellular")') 
c      PLACEINT=ISR_GKK(451,nx_trans)
c      write(*,'(7(F14.5,1X))') 
c     ' (GKK(PLACEINT+nzero,nx_trans),nzero=0,NQGP(0,451)-1)
c      write(*,'(7(F14.5,1X))') 
c     ' (GKK(PLACEINT+nzero,nx_ext),nzero=0,NQGP(0,451)-1)
c      write(*,'("nq: 452, transmembrane and extracellular")') 
c      PLACEINT=ISR_GKK(452,nx_trans)
c      write(*,'(7(F14.5,1X))') 
c     ' (GKK(PLACEINT+nzero,nx_trans),nzero=0,NQGP(0,452)-1)
c      write(*,'(7(F14.5,1X))') 
c     ' (GKK(PLACEINT+nzero,nx_ext),nzero=0,NQGP(0,452)-1)
c      write(*,'("nq: 453, transmembrane and extracellular")') 
c      PLACEINT=ISR_GKK(453,nx_trans)
c      write(*,'(7(F14.5,1X))') 
c     ' (GKK(PLACEINT+nzero,nx_trans),nzero=0,NQGP(0,453)-1)
c      write(*,'(7(F14.5,1X))') 
c     ' (GKK(PLACEINT+nzero,nx_ext),nzero=0,NQGP(0,453)-1)


C SGM 18Jan01 print out error statement if need to increace array sizes
      IF(Inc_NZ_GKK_M.GT.0) THEN
        WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)')Inc_NZ_GKK_M
        GOTO 9999
      ELSE IF(Inc_NISR_GKK_M.GT.0) THEN
        WRITE(ERROR,'(''>>Increase ISR_GKK_M to >= '',I12)'
     '    )Inc_NISR_GKK_M
        GOTO 9999
      ELSE IF(Inc_NISC_GKK_M.GT.0) THEN
        WRITE(ERROR,'(''>>Increase ISC_GKK_M to >= '',I12)'
     '    )Inc_NISC_GKK_M
        GOTO 9999
      ENDIF

      IF(.NOT.UPDATEDT) THEN
        IF(.NOT.UPDATE) THEN
          IF(IMPLICIT) THEN
            NZZT(1,NRLIST(1),nx_trans)=PLACEINT
            IF(SPARSEGKK(nx_trans).EQ.5) THEN
              COUNTINT=1
              DO nrr=1,NRLIST(0)
                nr=NRLIST(nrr)
                DO nq=NQR(1,nr),NQR(2,nr)
                  DO nzero=1,NQGP(0,nq)
                    ISC_GKK(PLACEINT+COUNTINT,nx_trans)=NQGP(nzero,nq)
                    COUNTINT=COUNTINT+1
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDIF
          IF(BIDOMAIN) THEN
            NZZT(1,NRLIST(1),nx_ext)=PLACEEXT
            IF(SPARSEGKK(nx_ext).EQ.5) THEN
              COUNTEXT=1
              DO nrr=1,NRLIST(0)
                nr=NRLIST(nrr)
                DO nq=NQR(1,nr),NQR(2,nr)
                  DO nzero=1,NQGP(0,nq)
                    ISC_GKK(PLACEEXT+COUNTEXT,nx_ext)=NQGP(nzero,nq)
                    COUNTEXT=COUNTEXT+1
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDIF
          IF(COUPBID) THEN
            NZZT(1,NRLIST(1),nx_upd)=PLACEEXT
            IF(SPARSEGKK(nx_upd).EQ.5) THEN
              COUNTEXT=1
              DO nrr=1,NRLIST(0)
                nr=NRLIST(nrr)
                DO nq=NQR(1,nr),NQR(2,nr)
                  DO nzero=1,NQGP(0,nq)
                    ISC_GKK(PLACEEXT+COUNTEXT,nx_upd)=NQGP(nzero,nq)
                    COUNTEXT=COUNTEXT+1
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDIF
      ENDIF ! not updatedt

      IF(.NOT.COUPBID) THEN
        UP_GRID_MATERIAL=.FALSE.
        UP_GRID_TENSOR=.FALSE.
      ENDIF
      FIRST_A=.TRUE.
      UPDATE_MATRIX=.TRUE.

      IF(DOP) THEN
        nr=NRLIST(1)
        WRITE(OP_STRING,
     '    '(/'' Global stiffness matrix GKK - transmembrane:'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' NOT(1,1,nr,nx)='',I5,'
     '    //''', NOT(2,1,nr,nx)='',I5)') NOT(1,1,nr,nx_trans),
     '    NOT(2,1,nr,nx_trans)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        CALL OPSTFMAT(DUMMY_LIST,ISC_GKK,ISR_GKK,IOOP,
     '    NOT(1,1,nr,nx_trans),NOT(2,1,nr,nx_trans),
     '    NZZT(1,nr,nx_trans),DUMMY_LIST,SPARSEGKK(nx_trans),
     '    GKK,GKK,'GKK','GKK',.TRUE.,.TRUE.,.FALSE.,
     '    ERROR,*9999)
        IF(BIDOMAIN) THEN
          WRITE(OP_STRING,
     '      '(/'' Global stiffness matrix GKK - external:'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NOT(1,1,nr,nx)='',I5,'
     '      //''', NOT(2,1,nr,nx)='',I5)') NOT(1,1,nr,nx_ext),
     '      NOT(2,1,nr,nx_ext)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CALL OPSTFMAT(DUMMY_LIST,ISC_GKK,ISR_GKK,IOOP,
     '      NOT(1,1,nr,nx_ext),NOT(2,1,nr,nx_ext),
     '      NZZT(1,nr,nx_ext),DUMMY_LIST,SPARSEGKK(nx_ext),
     '      GKK,GKK,'GKK','GKK',.TRUE.,.TRUE.,.FALSE.,
     '      ERROR,*9999)
        ENDIF
      ENDIF

      IF(IWRIT5(NRLIST(1),nx).GE.1) THEN
        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)
        WRITE(OP_STRING,'(1X,''Time for matrix assembly '',F8.2,'
     '    //'''s cpu'')') ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('ASSEMBLE10_FE')
      RETURN
 9998 ERROR=' '
 9999 CALL ERRORS('ASSEMBLE10_FE',ERROR)
      CALL EXITS('ASSEMBLE10_FE')
      RETURN 1
      END


      SUBROUTINE MARCH3(IBT,INP,ITHRES,
     '  NBJ,NEELEM,NHP,NKE,NKH,NPF,NPNE,NPNODE,
     '  nr,NRE,NVJE,NWQ,nx,NXI,NYNP,
     '  CE,CG,CP,PG,SE,THRES,XA,XE,XG,XIG,XP,YG,FIX,ERROR,*)

C#### Subroutine: MARCH3
C###  Description:
C###    MARCH3 performs time integration of cardiac activation equations
C###    with fixed time step algorithm.

C**** YP(ny,1) is solution vector at time T+DT
C**** YP(ny,3) is incremental boundary conditions
C****       4  "  solution vector at time T
C****       5  "  reaction vector at time T+DT
C**** FIX(ny,5) is .true. for input from a file (FILE07)

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b12.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:iwrit00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
      INCLUDE 'cmiss$reference:loc00.inc'
      INCLUDE 'cmiss$reference:marc00.cmn'
      INCLUDE 'cmiss$reference:time01.cmn'      
      INCLUDE 'cmiss$reference:time02.cmn'      
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ITHRES(3,NGM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NHP(NPM),
     '  NKE(NKM,NNM,NBFM,NEM),NKH(NHM,NPM,NCM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  nr,NRE(NEM),NVJE(NNM,NBFM,NJM,NEM),
     '  NWQ(6,0:NQM,NAM),nx,NXI(-NIM:NIM,0:4,0:NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),CP(NMM,NPM),PG(NSM,NUM,NGM,NBM),
     '  SE(NSM,NBFM,NEM),THRES(3,NGM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),XIG(NIM,NGM,NBM),
     '  XP(NKM,NVM,NJM,NPM),YG(NGM,NJM,NEM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM)
!     Local Variables
      INTEGER IBEG,IBEG1,IBEG5,IEND,IEND1,IEND5,
     '  MAXIT,nb,nc,ND_COL(9),ne,ng,nh,nhx,NITB,nk,
     '  no_col,noelem,nonode,no_unit,np,nq,NSE1,NSE2,NSE3,
     '  NSTEP,NT_COL,NT_UNIT,NUAT,ny,NY_UNIT(20)
      REAl ELAPSED_TIME,TIME_START(1),TIME_STOP(1)
      REAL*8 T
      CHARACTER CHAR4*4,CHAR*5,FILE*100,FORMAT*100,
     '  NAME_COL(9)*20
      LOGICAL CONTINUE,FILEIP
      SAVE NSTEP

      CALL ENTERS('MARCH3',*9999)

      nc=1 !temporary cpb 22/11/94
      DT=TINCR
      FILEIP=.FALSE.
      nb=NBJ(1,1)
      NITB=NIT(nb)
                 
      IF(RESTART) THEN
        T=T_RESTART
      	IF(NODE_HISTORY_OP) THEN
      	  CALL TRIM(FILE02,IBEG,IEND)
      	  CALL OPENF(IOFILE2,'DISK',FILE02(IBEG:IEND)//'.ipdata','OLD',
     '	    'APPEND','FORMATTED',132,ERROR,*9999)
          IF(NITB.EQ.2) THEN
      	    CALL OPENF(IOFILE3,'DISK',FILE02(IBEG:IEND)//'.history',
     '        'OLD','APPEND','FORMATTED',132,ERROR,*9999)
          ELSE IF(NITB.EQ.3) THEN
      	    CALL OPENF(IOFILE3,'DISK',FILE02(IBEG:IEND)//'.history3',
     '        'OLD','APPEND','FORMATTED',132,ERROR,*9999)
      	  ENDIF
        ENDIF
      ELSE IF(.NOT.RESTART) THEN !perform initial tasks
      	T=TSTART
      	no_unit=0
      	NSTEP=0

C *** 	Initialize ITHRES array
      	DO noelem=1,NEELEM(0,nr)
      	  ne=NEELEM(noelem,nr)
      	  nb=NBJ(1,ne)
      	  DO ng=1,NGT(nb)
      	    IF(ITHRES(1,ng,ne).EQ.0) THEN !ng currently not active
      	      IF(T.GE.YG(ng,1,ne)) ITHRES(1,ng,ne)=1
      	    ENDIF
      	  ENDDO
      	ENDDO
      	DO nonode=1,NPNODE(0,nr)
      	  np=NPNODE(nonode,nr)
          WRITE(CHAR4,FMT='(I4)') np
      	  CALL TRIM(CHAR4,IBEG1,IEND1)
      	  DO nhx=1,NHP(np)
            nh=NH_LOC(nhx,nx)
     	    DO nk=1,NKH(nh,np,nc)
      	      ny=NYNP(nk,1,nh,np,0,nc,nr)
      	      IF(nk.EQ.1.AND.FIX(ny,5)) THEN
      		no_unit=no_unit+1
      		NY_UNIT(no_unit)=20+no_unit
      		CALL TRIM(FILE07,IBEG5,IEND5)
      		FILE=FILE07(IBEG5:IEND5)//'.iphist_'
     '		  //CHAR4(IBEG1:IEND1)
      		FILEIP=.TRUE.
      		CALL TRIM(FILE,IBEG,IEND)
      		CALL OPENF(NY_UNIT(no_unit),'DISK',
     '		  FILE(IBEG:IEND),'OLD','SEQUEN','FORMATTED',
     '		  132,ERROR,*9999)
      		IF(DOP) THEN
C$        call mp_setlock()
      		  WRITE(OP_STRING,'('' FILE='',A)') FILE(IBEG:IEND)
      		  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C$        call mp_unsetlock()
      		ENDIF
      	      ENDIF
      	    ENDDO
      	  ENDDO
      	ENDDO
      	NT_UNIT=no_unit

      	IF(NODE_HISTORY_OP) THEN
      	  CALL TRIM(FILE02,IBEG,IEND)
      	  CALL OPENF(IOFILE2,'DISK',FILE02(IBEG:IEND)//'.ipdata','NEW',
     '	    'SEQUEN','FORMATTED',132,ERROR,*9999)
          IF(NITB.EQ.2) THEN
      	    CALL OPENF(IOFILE3,'DISK',FILE02(IBEG:IEND)//'.history',
     '        'NEW','SEQUEN','FORMATTED',132,ERROR,*9999)
          ELSE IF(NITB.EQ.3) THEN
      	    CALL OPENF(IOFILE3,'DISK',FILE02(IBEG:IEND)//'.history3',
     '        'NEW','SEQUEN','FORMATTED',132,ERROR,*9999)
          ENDIF
      	  NT_COL=NODE_HISTORY(0)+1
      	  DO no_col=1,NT_COL
      	    ND_COL(no_col)=no_col
      	  ENDDO
      	  WRITE(IOFILE2,'(10I2)') NT_COL,
     '      (ND_COL(no_col),no_col=1,NT_COL)
      	  DO no_col=1,NT_COL
      	    IF(no_col.EQ.1) THEN
      	      NAME_COL(no_col)='Time'
      	    ELSE
      	      np=NODE_HISTORY(no_col-1)
              WRITE(CHAR,FMT='(I5)') np
      	      CALL TRIM(CHAR,IBEG,IEND)
      	      NAME_COL(no_col)='Node '//CHAR(IBEG:IEND)
      	    ENDIF
      	    WRITE(IOFILE2,'(A)') NAME_COL(no_col)
      	  ENDDO
      	ENDIF  !NODE_HISTORY_OP

      	IF(NGT(nb).EQ.3**NITB) THEN
c     	  NMAX=3
      	ELSE IF(NGT(nb).EQ.5**NITB) THEN
c     	  NMAX=5
      	ELSE
      	  ERROR='>>>Incorrect # of Gauss points'
      	  GO TO 9999
      	ENDIF
      	CALL CPCG(1,NBJ(1,1),NPNE(1,1,1),nr,nx,CE(1,1),CG,CP,PG,
     '    ERROR,*9999)
      	NSE1=NINT(CG(5,1))  !# of Gauss points along a search in Xi1
      	NSE2=NINT(CG(6,1))  !# of Gauss points along a search in Xi2
      	NSE3=NINT(CG(7,1))  !# of Gauss points along a search in Xi3
      	IF(DOP) THEN
C$       call mp_setlock()
      	  WRITE(OP_STRING,'('' NSE1='',I3,'' NSE2='',I3,'
     '	    //''' NSE3='',I3)') NSE1,NSE2,NSE3
      	  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C$       call mp_unsetlock()
      	ENDIF
      	MAXIT=0
      	IF(IWRIT4(nr,nx).GT.0.
     '    AND.ITYP3(nr,nx).ne.3.AND.ITYP3(nr,nx).ne.5) THEN
      	  WRITE(OP_STRING,
     '	    '(//'' ** Isochrone #1 ** (initial conditions)''/)')
      	  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      	ENDIF
      	DO noelem=1,NEELEM(0,1)
      	  ne=NEELEM(noelem,1)
      	  nb=NBJ(1,ne)
      	  MAXIT=MAXIT+NGT(nb)
      	  IF(IWRIT4(nr,nx).GT.0) THEN
      	    IF(IWRIT4(nr,nx).GT.1.
     '        AND.ITYP3(nr,nx).ne.3.AND.ITYP3(nr,nx).ne.5) THEN
              WRITE(CHAR4,FMT='(I4)') ne
      	      CALL TRIM(CHAR4,IBEG,IEND)
      	      FORMAT='('' Element #'//CHAR4(IBEG:IEND)//':'')'
      	      WRITE(OP_STRING,FORMAT)
      	      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      	      WRITE(OP_STRING,'('' ITHRES: '',125I1)') 
     '		(ITHRES(1,ng,ne),ng=1,NGT(nb))
      	      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      	    ENDIF
      	  ENDIF
      	ENDDO
      	MAXIT=MAXIT+50

      ENDIF !.not.restart

      CALL CPU_TIMER(CPU_USER,TIME_START)

      CONTINUE=.TRUE.
      DO WHILE(CONTINUE)
        NSTEP=NSTEP+1
        WRITE(CHAR4,FMT='(I4)') NSTEP+1
      	CALL TRIM(CHAR4,IBEG,IEND)
      	WRITE(OP_STRING,
     '	  '('' ** Isochrone #'//CHAR4(IBEG:IEND)//'**'',(at time '','
     '	  //'D11.4,'')'')') T+DT
      	CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      	CALL TFRONT(IBT,INP,ITHRES,
     '    NBJ,NEELEM,NKE,NPF,NPNE,
     '    NRE,NSE1,NSE2,NSE3,NVJE,NXI,
     '	  CE,CG,CP,PG,SE,T,THRES,XA,XE,XG,XIG,XP,YG,
     '	  ERROR,*9999)
C Subroutine FRONT Archived GBS 27-OCT-1994
c        ELSE IF(ITYP3(nr,nx).EQ.3) THEN  !FitzHugh-Nagumo model
c          CALL FRONT(NBJ,NGAP,NKE,NQGE,NWQ(1,0,1),NXQ(-NIM,0,0,1),CQ,
c     '      DNUDXQ,DXDXIQ,GUQ,PROPQ,T,XA,XQ,ZA,ERROR,*9999)

C ***   Output Gauss point values

        IF(IWRIT4(nr,nx).GT.1) THEN
      	  DO noelem=1,NEELEM(0,1)
      	    ne=NEELEM(noelem,1)
      	    nb=NBJ(1,ne)
            WRITE(CHAR4,FMT='(I4)') ne
      	    CALL TRIM(CHAR4,IBEG,IEND)
      	    FORMAT='('' Element #'//CHAR4(IBEG:IEND)//':'')'
      	    WRITE(OP_STRING,FORMAT)
      	    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      	    WRITE(OP_STRING,'('' ITHRES: '',125I1 )') 
     '	      (ITHRES(1,ng,ne),ng=1,NGT(nb))
      	    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      	    WRITE(OP_STRING,'('' YG(ng,1,ne): '',25F8.1)')
     '	       (YG(ng,1,ne),ng=1,NGT(nb))
      	    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      	    WRITE(OP_STRING,'('' THRES(1): '',25F8.1)')
     '	       (THRES(1,ng,ne),ng=1,NGT(nb))
      	    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      	  ENDDO
        ENDIF

c       IF(UPVU) THEN
c         IF(ISOAREA) THEN !isochrone areas defined
C ***       Update colour of isochrone areas on output windows
c           DO noiw=1,NTIW
c             IW=IWK(noiw)
c             IF(IWKS(IW).GT.0) THEN
c               IF(ISISOC(IW).GT.0) THEN
c                 CALL ACWK(IW,1,ERROR,*9999)
c                 DO noelem=1,NEELEM(0,1)
c                   ne=NEELEM(noelem,1)
c                   IF(DOP) THEN
c                     WRITE(OP_STRING,'('' element '',I4)') ne
c                     CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c                   ENDIF
c                   nb=NBJ(1,ne)
c                   DO ng=1,NGT(nb)
c                     IFILL=ng+(ne-1)*NGT(nb) !is fill area index
c                     IF(COLOUR_WS) THEN !use colours ICOL=17..249
c                       ICOL=249-NINT((THRES(2,ng,ne)-ZMINI)
c    '                    /ZDIFF*232.0D0)
c                       IF(DOP) THEN
c                         WRITE(OP_STRING,'('' IFILL='',I4,'
c    '                      //''' ICOL='',I4)') IFILL,ICOL
c                         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c                       ENDIF
c                       CALL SET_FILL_REP(1,IFILL,'SOLID',1,ICOL,ERROR,
c    '                    *9999)
c                     ELSE    !monochrome patterns IPAT=1..16
c                       IPAT=16-NINT((THRES(2,ng,ne)-ZMINI)
c    '                    /ZDIFF*15.0D0)
c                       IF(DOP) THEN
c                         WRITE(OP_STRING,'('' IFILL='',I4,'
c    '                      //''' IPAT='',I4)') IFILL,IPAT
c                         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c                       ENDIF
c                       CALL SET_FILL_REP(1,IFILL,'PATTERN',
c    '                    12+(17-IPAT),1,ERROR,*9999)
c                     ENDIF
c                   ENDDO
c                 ENDDO
c                 CALL DAWK(IW,1,ERROR,*9999)
c               ENDIF
c             ENDIF
c           ENDDO         
c
c         ELSE IF(.NOT.ISOAREA) THEN !isochrone points defined
C ***       Update visibility of Gauss points on output windows
c           DO noiw=1,NTIW
c             IW=IWK(noiw)
c             IF(IWKS(IW).GT.0) THEN
c               CALL ACWK(IW,1,ERROR,*9999)
c               IF(ITYP3(nr,nx).EQ.3.OR.ITYP3(nr,nx).EQ.5) THEN
c                 DO nq=1,NQT
c                   IF(COLOUR_WS) THEN
c                     IVAL=INT(ZA(1,1,1,nq)*215.0D0)+1
c                     IF(IVAL.LT.1) IVAL=1
c                     IF(IVAL.GT.216) IVAL=216
c                     INDEX=INDEX_POLYMARKER(IVAL,'PLUS','SIZE1',
c    '                  'BLACK')
c!!! This call (and routine?) need to use XQ and YQ
c                      CALL SGGRID(INDEX,ISEG,ISGRID(IW),IW,nq,nq,
c     '                  CSEG,XA(1,1,nq),ERROR,*9999)
c                   ELSE
c                     IF(NWQ(4,nq,1).EQ.1) THEN
c                       CALL VISIB(IW,ISEG,ISGRID(IW),'VISIBLE',
c    '                    ERROR,*9999)
c                     ELSE
c                       CALL VISIB(IW,ISEG,ISGRID(IW),'INVISIBLE',
c    '                    ERROR,*9999)
c                     ENDIF
c                   ENDIF
c                 ENDDO
c               ELSE
c                 DO noelem=1,NEELEM(0,1)
c                   ne=NEELEM(noelem,1)
c                   nb=NBJ(1,ne)
c                   DO ng=1,NGT(nb)
c                     IF(ISGAUS(IW,ng,ne).GT.0) THEN
c                       IF(ITHRES(1,ng,ne).EQ.1) THEN
c                         CALL VISIB(IW,ISEG,ISGAUS(IW,ng,ne),'VISIBLE',
c    '                      ERROR,*9999)
c                       ELSE IF(ITHRES(1,ng,ne).ne.1) THEN
c                         CALL VISIB(IW,ISEG,ISGAUS(IW,ng,ne),
c    '                      'INVISIBLE',ERROR,*9999)
c                       ENDIF
c                     ENDIF
c                   ENDDO
c                 ENDDO
c               ENDIF
c               CALL DAWK(IW,1,ERROR,*9999)
c             ENDIF
c           ENDDO
c         ENDIF
c       ENDIF !upvu

C ***   Check whether any Gauss pts are partially active
        NUAT=0
        IF(ITYP2(nr,nx).EQ.9) THEN
          DO nq=1,NQT
            IF(NWQ(4,nq,1).GE.1) NUAT=NUAT+1
          ENDDO
        ELSE
          DO noelem=1,NEELEM(0,1)
            ne=NEELEM(noelem,1)
            nb=NBJ(1,ne)
            DO ng=1,NGT(nb)
C             IF(ITHRES(1,ng,ne).EQ.1) NUAT=NUAT+1
              IF(ITHRES(1,ng,ne).LE.1) NUAT=NUAT+1
            ENDDO
          ENDDO
        ENDIF

C ***   Check for termination
        T=T+DT
        IF(NUAT.EQ.0.OR.T.GE.TFINISH) THEN
          CONTINUE=.FALSE.
        ENDIF
      	IF(ITYP2(nr,nx).ne.9.AND.NSTEP.GE.MAXIT) THEN
      	  WRITE(OP_STRING,'('' Exceeded max iterations ='',I8)')MAXIT
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CONTINUE=.FALSE.
      	ENDIF

        !Check for numeric keypad entry interrupt
        CALL GETSTR2(ERROR,*9998)

      ENDDO
      CALL CLOSEF(IOFILE2,ERROR,*9999)
      IF(NODE_HISTORY_OP) CALL CLOSEF(IOFILE3,ERROR,*9999)
      
C *** Timing information
      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)
      WRITE(OP_STRING,'(//'' Solution took '',I8,'
     '  //''' iterations and '',D11.4,'' s of cpu time'')') 
     '  NSTEP,ELAPSED_TIME
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

C  Commented out GBS 29/7/93
C  This file doesn't seem to be needed any longer
C      IF(ITYP3(nr,nx).EQ.3.OR.ITYP3(nr,nx).EQ.5) THEN  !FHN/bidomain
C      	CALL OPENF(IOFILE2,'DISK',FILE02(IBEG:IEND)//'_time.ipdata',
C     '	  'NEW','SEQUEN','FORMATTED',132,ERROR,*9999)
C        WRITE(IOFILE2,'('' FHN grid point potential values at time'','
C     '	  //'F10.2)') T
C      	DO nq=1,NQT
C      	  WRITE(IOFILE2,'(I5,8D12.5)') nq,(XQ(nj,nq),nj=1,3),
C     '	    ZA(1,1,1,nq),1.0D0,1.0D0,1.0D0,1.0D0
C      	ENDDO
C        CALL CLOSEF(IOFILE2,ERROR,*9999)
C      ENDIF
      
C     IF(IWRIT4(nr,nx).GT.0.AND.NITB.EQ.3) THEN
C       CALL TRIM(FILE02,IBEG,IEND)
C       CALL OPENF(9,'DISK',FILE02(IBEG:IEND)//'.ipdat','NEW',
C    '    'SEQUEN','FORMATTED',132,ERROR,*9999)
C       WRITE(9,'('' Model epicardial Gauss point activation times'')')
C       nd=0
C       DO ne=13,24
C         CALL XPXE(NBJ(1,ne),NKE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
C    '      NQGE(1,ne),nr,NVJE(1,1,1,ne),
C    '      SE(1,1,ne),XA,XE,XP,ERROR,*9999)
C         DO ng=19,27
C           CALL XEXG(NBJ(1,ne),ng,nr,PG,XE,XG,ERROR,*9999)
C           XG(4,1)=THRES(1,ng,ne)
C           nd=nd+1
C           WRITE(9,'(I3,4D12.4,'' 1 1 1 1'')') nd,(XG(nj,1),nj=1,4)
C         ENDDO
C         CALL CLOSEF(9,ERROR,*9999)
C       ENDDO
C     ENDIF
      IF(FILEIP) THEN
        DO no_unit=1,NT_UNIT
          CALL CLOSEF(NY_UNIT(no_unit),ERROR,*9999)
        ENDDO
      ENDIF

      
      T_RESTART=T   !Just in case we come back
      CALL EXITS('MARCH3')
      RETURN

 9998 T_RESTART=T
      CALL EXITS('MARCH3')
      RETURN

 9999 IF(FILEIP) THEN
        DO no_unit=1,NT_UNIT
          CLOSE(UNIT=no_unit)
        ENDDO 
      ENDIF
      CALL ERRORS('MARCH3',ERROR)
      CALL EXITS('MARCH3')
      RETURN 1
      END


      SUBROUTINE MARCH5(NBJ,NEELEM,NRLIST,NWQ,nx,NXQ,AQ,
     '  CQ,GCHQ,GUQ,PROPQ,XQ,YQ,ERROR,*)

C#### Subroutine: MARCH5
C###  Description:
C###    MARCH5 performs collocation solution of bidomain/monodomain 
C###    equations.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b12.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:iwrit00.cmn'
      INCLUDE 'cmiss$reference:ktyp30.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
      INCLUDE 'cmiss$reference:loc00.inc'
      INCLUDE 'cmiss$reference:marc00.cmn'
      INCLUDE 'cmiss$reference:nqloc00.inc'
      INCLUDE 'cmiss$reference:time01.cmn'
      INCLUDE 'cmiss$reference:time02.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NRLIST(0:NRM),
     '  NWQ(6,0:NQM,NAM),
     '  nx,NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      REAL*8 AQ(NMAQM,NQM),CQ(NMM,NQM),GCHQ(3,NQM),GUQ(3,3,NQM),
     '  PROPQ(3,3,4,2,NQM),XQ(NJM,NQM),YQ(NYQM,NIQM,NAM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER COUNT,HIST_COUNT,IBEG,IEND,niq,niqCAI,niqCALCIUM,niqD,
     '  niqDRECOV,niqDV,niqF,niqH,niqJ,niqMM,niqOLDSOLN,niqRECOV,niqSAC,
     '  niqV,niqX,NITB,nj,njj2,no_col,nq,NQFIRST2,NQFIRST3,NQSTART,
     '  nr,NSTEP,NT_COL,NUMTIMES
      REAL ELAPSED_TIME,TIME_START(1),TIME_STOP(1)
      REAL*8 T
      CHARACTER CHAR10*10,NAME_COL(9)*20
      LOGICAL CONTINUE
      SAVE HIST_COUNT,NSTEP,NQSTART

      CALL ENTERS('MARCH5',*9999)

      DT=TINCR

C old MLB 3/9/97
C! Find basis function using Extended Lagrange 
C      nb=1
C      DO WHILE (nb.LE.NBT.AND.NBC(nb).NE.7)
C        nb=nb+1
C      ENDDO
C      CALL ASSERT((NBC(nb).EQ.7),
C     '  'Extended basis function not defined',ERROR,*9999)
C      NITB=NIT(nb)

! Temp until think of a better way
      NITB=NJT

      nr=NRLIST(1)
      IF(ITYP3(nr,nx).EQ.1) THEN
C       Cubic w/o recovery
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqV,NIQ_V,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqOLDSOLN,NIQ_OLDSOLN,
     '    ERROR,*9999)
      ELSE IF(ITYP3(nr,nx).EQ.2) THEN
C       FitzHugh-Nagumo
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqV,NIQ_V,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqRECOV,NIQ_RECOV,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqCALCIUM,NIQ_CALCIUM,
     '    ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqSAC,NIQ_SAC,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqOLDSOLN,NIQ_OLDSOLN,
     '    ERROR,*9999)
      ELSE IF(ITYP3(nr,nx).EQ.3) THEN
C       van Capelle-Durrer
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqV,NIQ_V,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqRECOV,NIQ_RECOV,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqCALCIUM,NIQ_CALCIUM,
     '    ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqSAC,NIQ_SAC,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqDV,NIQ_DV,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqDRECOV,NIQ_DRECOV,
     '    ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqOLDSOLN,NIQ_OLDSOLN,
     '    ERROR,*9999)
      ELSE IF(ITYP3(nr,nx).EQ.4) THEN
C       Beeler-Reuter
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqV,NIQ_V,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqX,NIQ_X,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqCAI,NIQ_CAI,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqMM,NIQ_M,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqH,NIQ_H,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqJ,NIQ_J,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqD,NIQ_D,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqF,NIQ_F,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqOLDSOLN,NIQ_OLDSOLN,
     '    ERROR,*9999)
      ELSE IF(ITYP3(nr,nx).EQ.5) THEN
C       Jafri-Rice-Winslow
        CALL ASSERT(.FALSE.,'>>Not implemented',ERROR,*9999)
      ELSE IF(ITYP3(nr,nx).EQ.6) THEN
C       Luo-Rudy
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqV,NIQ_V,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqX,NIQ_X,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqCAI,NIQ_CAI,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqM,NIQ_M,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqH,NIQ_H,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqJ,NIQ_J,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqD,NIQ_D,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqF,NIQ_F,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqOLDSOLN,NIQ_OLDSOLN,
     '    ERROR,*9999)            
      ENDIF !Model type

      IF(RESTART) THEN
        T=T_RESTART
      	IF(NODE_HISTORY_OP) THEN
      	  CALL TRIM(FILE02,IBEG,IEND)
      	  CALL OPENF(IOFILE2,'DISK',FILE02(IBEG:IEND)//'.ipdata','OLD',
     '	    'APPEND','FORMATTED',132,ERROR,*9999)
      	  CALL OPENF(IOFILE5,'DISK',FILE02(IBEG:IEND)//'.ipdat2','OLD',
     '	    'APPEND','FORMATTED',132,ERROR,*9999)
          IF(NITB.EQ.2) THEN
      	    CALL OPENF(IOFILE3,'DISK',FILE02(IBEG:IEND)//'.history',
     '        'OLD','APPEND','FORMATTED',132,ERROR,*9999)
          ELSE IF(NITB.EQ.3) THEN
      	    CALL OPENF(IOFILE3,'DISK',FILE02(IBEG:IEND)//'.history3',
     '        'OLD','APPEND','FORMATTED',132,ERROR,*9999)
      	  ENDIF
        ENDIF
      ELSE IF(.NOT.RESTART) THEN ! perform initial tasks
      	T=TSTART
      	NSTEP=0
        HIST_COUNT=1

      	IF(NODE_HISTORY_OP) THEN
      	  CALL TRIM(FILE02,IBEG,IEND)
      	  CALL OPENF(IOFILE2,'DISK',FILE02(IBEG:IEND)//'.ipdata','NEW',
     '	    'SEQUEN','FORMATTED',132,ERROR,*9999)
      	  CALL OPENF(IOFILE5,'DISK',FILE02(IBEG:IEND)//'.ipdat2','NEW',
     '	    'SEQUEN','FORMATTED',132,ERROR,*9999)
          IF(NITB.EQ.2) THEN
      	    CALL OPENF(IOFILE3,'DISK',FILE02(IBEG:IEND)//'.history',
     '        'NEW','SEQUEN','FORMATTED',132,ERROR,*9999)
          ELSE IF(NITB.EQ.3) THEN
      	    CALL OPENF(IOFILE3,'DISK',FILE02(IBEG:IEND)//'.history3',
     '        'NEW','SEQUEN','FORMATTED',132,ERROR,*9999)
          ENDIF
      	  NT_COL=NODE_HISTORY(0)+1
      	  WRITE(IOFILE2,'(10I2)') NT_COL,(no_col,no_col=1,NT_COL)
      	  WRITE(IOFILE5,'(12I3)') NIQM+1,(no_col,no_col=1,NIQM+1)
      	  DO no_col=1,NT_COL
      	    IF(no_col.EQ.1) THEN
      	      NAME_COL(no_col)='Time'
      	    ELSE
      	      nq=NODE_HISTORY(no_col-1)
              WRITE(CHAR10,FMT='(I6)') nq
      	      CALL TRIM(CHAR10,IBEG,IEND)
      	      NAME_COL(no_col)='Node '//CHAR10(IBEG:IEND)
      	    ENDIF
      	    WRITE(IOFILE2,'(A)') NAME_COL(no_col)
      	  ENDDO
          DO no_col=1,NIQM
            WRITE(IOFILE5,'(''YQ('',I6,'','',I2,'',1)'')') 
     '        NODE_HISTORY(1),no_col
          ENDDO
      	  WRITE(IOFILE2,'(I5,9D12.4)') NSTEP,T,
     '	    (YQ(NODE_HISTORY(no_col),1,1),no_col=1,NODE_HISTORY(0))
      	  IF(NITB.EQ.2) THEN  ! history file for reading into Explorer
      	    nq=1
      	    DO nj=1,2    ! Find grid point in bottom left corner
      	      DO WHILE(NXQ(-nj,1,nq,1).ne.0)
      		nq=NXQ(-nj,1,nq,1)
      	      ENDDO
      	    ENDDO
      	    NQSTART=nq
      	    COUNT=1   ! Find number of grid points in x dir'n
      	    DO WHILE(NXQ(1,1,nq,1).ne.0)
      	      nq=NXQ(1,1,nq,1)
      	      COUNT=COUNT+1
      	    ENDDO
      	    WRITE(IOFILE3,'(I4)') COUNT  ! and write to file
      	    nq=NQSTART
      	    COUNT=1   ! now in y dir'n
      	    DO WHILE(NXQ(2,1,nq,1).ne.0)
      	      nq=NXQ(2,1,nq,1)
      	      COUNT=COUNT+1
      	    ENDDO
      	    WRITE(IOFILE3,'(I4)') COUNT
c      	    IF(NITB.EQ.3) THEN
c      	      nq=NQSTART
c      	      COUNT=1   ! now in z dir'n
c      	      DO WHILE(NXQ(3,1,nq,1).ne.0)
c      		nq=NXQ(3,1,nq,1)
c      		COUNT=COUNT+1
c      	      ENDDO
c      	      WRITE(IOFILE3,'(I4)') COUNT
c      	    ENDIF
C* Total number of time steps written to file
      	    IF(IWRIT1(NRLIST(1),nx).EQ.0) THEN
      	      NUMTIMES=0
      	    ELSE
      	      NUMTIMES=INT((T1-T0)/DT/IWRIT1(NRLIST(1),nx))
      	    ENDIF
      	    WRITE(IOFILE3,'(I4)') NUMTIMES
      	    nq=NQSTART
      	    NQFIRST2=NQSTART
      	    NQFIRST3=NQSTART
      	    DO WHILE(nq.ne.0)
      	      DO WHILE(nq.ne.0)
      		DO WHILE(nq.ne.0)
      		  WRITE(IOFILE3,'(F10.4)') YQ(nq,1,1)
      		  nq=NXQ(1,1,nq,1)
      		ENDDO
      		NQFIRST2=NXQ(2,1,NQFIRST2,1)
      		nq=NQFIRST2
      	      ENDDO
      	      NQFIRST3=NXQ(3,1,NQFIRST3,1)
      	      NQFIRST2=NQFIRST3
      	      nq=NQFIRST3
      	    ENDDO
      	  ELSE IF(NITB.EQ.3) THEN   ! history3 file
C* File format for 3D point viewing program :
C* - number of points
C* - number of timesteps (inc. t=0)
C* - point coords
C* - potential values looped by time, point
      	    WRITE(IOFILE3,'(I6)') NQT
      	    IF(IWRIT1(NRLIST(1),nx).EQ.0) THEN
      	      NUMTIMES=0
      	    ELSE
      	      NUMTIMES=INT((T1-T0)/DT/IWRIT1(NRLIST(1),nx))
      	    ENDIF
      	    WRITE(IOFILE3,'(I6)') NUMTIMES+1
      	    DO nq=1,NQT
      	      WRITE(IOFILE3,'(3F12.5)') (XQ(nj,nq), nj=1,3)
      	    ENDDO
      	    WRITE(IOFILE3,'(8F10.4)') (YQ(nq,1,1), nq=1,NQT)
      	  ENDIF
      	  HIST_COUNT=1
      	ENDIF   ! NODE_HISTORY_OP

      	IF(KTYP36.EQ.0) THEN   ! Not using DTAR
      	  DO nq=1,NQT
      	    NWQ(5,nq,1)=nq
      	    NWQ(4,nq,1)=2
      	  ENDDO
      	  NWQ(5,0,1)=NQT
      	ELSE                   ! Using DTAR
      	  DO nq=1,NQT
      	    NWQ(5,nq,1)=0
      	    NWQ(4,nq,1)=0
      	  ENDDO
      	  NWQ(5,0,1)=0
      	ENDIF                  ! Not using DTAR
      ENDIF ! .not.restart

      CALL CPU_TIMER(CPU_USER,TIME_START)

      CONTINUE=.TRUE.
      DO WHILE(CONTINUE)
        NSTEP=NSTEP+1

      	CALL BFRONT(NBJ,NEELEM,niqCAI,niqCALCIUM,niqD,niqDRECOV,
     '    niqDV,niqF,niqH,niqJ,niqMM,niqOLDSOLN,niqRECOV,niqSAC,niqV,
     '    niqX,NRLIST,NWQ(1,0,1),nx,NXQ(-NIM,0,0,1),
     '    AQ,CQ,GCHQ,GUQ,PROPQ,T,YQ,ERROR,*9999)

C ***   Output Gauss point values

      	IF(NODE_HISTORY_OP) THEN
      	  WRITE(IOFILE2,'(I6,9D12.4)') NSTEP,T+DT,
     '	    (YQ(NODE_HISTORY(no_col),1,1),no_col=1,NODE_HISTORY(0))
      	  WRITE(IOFILE5,'(I5,6D13.5,/,5X,6D13.5)') NSTEP,T+DT,
     '	    (YQ(NODE_HISTORY(1),1,1)-CQ(9,NODE_HISTORY(1)))/
     '      (CQ(10,NODE_HISTORY(1))-CQ(9,NODE_HISTORY(1))),  ! 0-1 pot
     '      (YQ(NODE_HISTORY(1),niq,1),niq=2,NIQM)
        ENDIF
      	IF(HIST_COUNT.EQ.IWRIT1(NRLIST(1),nx)) THEN
      	  HIST_COUNT=1
          IF(NODE_HISTORY_OP) THEN
      	    IF(NITB.EQ.2) THEN  ! history file
      	      nq=NQSTART
      	      NQFIRST2=NQSTART
      	      NQFIRST3=NQSTART
      	      DO WHILE(nq.ne.0)
      		DO WHILE(nq.ne.0)
      		  DO WHILE(nq.ne.0)
      		    WRITE(IOFILE3,'(F10.4)') YQ(nq,1,1)
      		    nq=NXQ(1,1,nq,1)
      		  ENDDO
      		  NQFIRST2=NXQ(2,1,NQFIRST2,1)
      		  nq=NQFIRST2
      		ENDDO
      		NQFIRST3=NXQ(3,1,NQFIRST3,1)
      		NQFIRST2=NQFIRST3
      		nq=NQFIRST3
      	      ENDDO
      	    ELSE IF(NITB.EQ.3) THEN  ! history3 file
      	      WRITE(IOFILE3,'(8F10.4)') (YQ(nq,1,1), nq=1,NQT)
      	    ENDIF
          ENDIF
          CALL CPU_TIMER(CPU_USER,TIME_STOP)
          ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)
      	  WRITE(OP_STRING,'(1X,I6,'' iterations : '',F8.2,'
     '	    //'''s cpu :'',''active:'',I6,''/'',I6)') 
     '	    NSTEP,ELAPSED_TIME,NWQ(4,0,1),NWQ(5,0,1)
      	  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      	ELSE
      	  HIST_COUNT=HIST_COUNT+1
      	ENDIF

C ***   Check for termination
        T=T+DT
        IF(T.GE.TFINISH.OR.NWQ(4,0,1).EQ.0) THEN
      	  CONTINUE=.FALSE.
        ENDIF

        ! Check for numeric keypad entry interrupt
        CALL GETSTR2(ERROR,*9998)

      ENDDO  ! CONTINUE
      
C *** Timing information
      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)
      WRITE(OP_STRING,'('' >>>Solution took '',I6,'
     '  //''' iterations and '',D11.4,'' s of cpu time'')') 
     '  NSTEP,ELAPSED_TIME
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      IF(IWRIT4(NRLIST(1),nx).GE.1) THEN
        CALL ASSERT(.FALSE.,'>>Obselete in march5',ERROR,*9999)
C**  Write activation times to .act file
      	CALL TRIM(FILE02,IBEG,IEND)
      	CALL OPENF(IOFILE4,'DISK',FILE02(IBEG:IEND)//'.act','NEW',
     '	  'SEQUEN','FORMATTED',132,ERROR,*9999)
        DO nq=1,NQT
      	  WRITE(IOFILE4,'(I7,<NJT+1>F10.4)')
     '	    nq,(XQ(NJ_LOC(NJL_GEOM,njj2,NRLIST(1)),nq),
     '      njj2=1,NJ_LOC(NJL_GEOM,0,NRLIST(1))),YQ(nq,1,3)
        ENDDO
        CALL CLOSEF(IOFILE4,ERROR,*9999)
      	IF(IWRIT4(NRLIST(1),nx).EQ.2) THEN
      	  WRITE(OP_STRING,'('' Transmembrane Potential:'')')
      	  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      	  WRITE(OP_STRING,'(1X,12F10.5)') (YQ(nq,1,1),nq=1,NQT)
      	  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      	  WRITE(OP_STRING,'('' Activation Time:'')')
      	  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      	  WRITE(OP_STRING,'(1X,12F10.5)') (YQ(nq,1,3),nq=1,NQT)
      	  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      	ENDIF
      ENDIF

 9998 CALL CLOSEF(IOFILE2,ERROR,*9999)
      IF(NODE_HISTORY_OP) THEN
        CALL CLOSEF(IOFILE3,ERROR,*9999)
        CALL CLOSEF(IOFILE5,ERROR,*9999)
      ENDIF
      T_RESTART=T   !Just in case we come back

      CALL EXITS('MARCH5')
      RETURN
 9999 CALL ERRORS('MARCH5',ERROR)
      CALL EXITS('MARCH5')
      RETURN 1
      END


      SUBROUTINE MARCH7(IBT,IDO,INP,ISC_GKK,ISR_GKK,
     '  NBH,NBJ,NEELEM,NENQ,NHE,NHP,NKE,NKH,NPF,NP_INTERFACE,NPNE,
     '  NPNODE,NQGP,NQGP_PIVOT,NQGW,NQNE,NQS,NQXI,NRLIST,NVHE,NVHP,NYNE,
     '  NYNP,NW,NWQ,NXLIST,NXQ,AQ,CQ,CURVCORRECT,GCHQ,GKK,GUQ,PROPQ,
     '  SE,YP,YQ,ZA,ZE,ZP,FIXQ,ITER8,ERROR,*)

C#### Subroutine: MARCH7
C###  Description:
C###    MARCH7 controls the implicit solution of grid activation 
C###    problems. It can integrate both the monodomain and bidomain
C###    equations.  It is also used in the iteration on epicardial
C###    fluxes and potentials to match them across the epicardial
C###    surface.
C***  Created by Martin Buist, May 1997

      IMPLICIT NONE

      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b12.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:grid00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:iwrit00.cmn'
      INCLUDE 'cmiss$reference:ktyp30.cmn'
      INCLUDE 'cmiss$reference:loc00.inc'
      INCLUDE 'cmiss$reference:maqloc00.inc'
      INCLUDE 'cmiss$reference:marc00.cmn'
      INCLUDE 'cmiss$reference:nqloc00.inc'
      INCLUDE 'cmiss$reference:solv00.cmn'
      INCLUDE 'cmiss$reference:time02.cmn'

!     Parameter list
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),ISC_GKK(NISC_GKKM,NXM),ISR_GKK(NISR_GKKM,NXM),
     '  NBH(NHM,NCM,NEM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NENQ(0:8,NQM),NHE(NEM,NXM),
     '  NHP(NPM,0:NRM,NXM),NKE(NKM,NNM,NBFM,NEM),
     '  NKH(NHM,NPM,NCM,0:NRM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NP_INTERFACE(0:NPM,0:3),
     '  NPNODE(0:NP_R_M,0:NRM),NQGP(0:19,NQM),NQGP_PIVOT(19,NQM),
     '  NQNE(NEM,NQEM),NQS(NEM),
     '  NQXI(0:NIM,NQSCM),NRLIST(0:NRM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),NW(NEM,3),NWQ(6,0:NQM),
     '  NXLIST(0:NXM),NXQ(-NIM:NIM,0:4,0:NQM)
      REAL*8 AQ(NMAQM,NQM),CQ(NMM,NQM),CURVCORRECT(2,2,NNM,NEM),
     '  GCHQ(3,NQM),GKK(NZ_GKK_M,NXM),GUQ(3,3,NQM),NQGW(NQM,19),
     '  PROPQ(3,3,4,2,NQM),SE(NSM,NBFM,NEM),YP(NYM,NIYM,NXM),
     '  YQ(NYQM,NIQM,NAM,NXM),ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
      LOGICAL FIXQ(NYQM,NIYFIXM,NXM),ITER8
!     Local variables
      INTEGER maqp1t0,maqp1t1,maqp1i,maqp2t0,maqp2t1,maqp2i,niqCAI,
     '  niqCALCIUM,niqD,niqDRECOV,niqDV,niqF,niqH,niqJ,niqMM,niqOLDSOLN,
     '  niqRECOV,niqSAC,niqV,niqX,niqBNDRY,NRTEMP(0:1),NSTEP,nxc,nx_ext,
     '  nx_torso,nx_trans,nx_upd,nq,nr
      REAL*8 IION(NQM),RHS(NQM),SOLNINT(NQM),SOLNEXT(NQM),T
      REAL ELAPSED_TIME,TIME_START(1),TIME_STOP(1)
      LOGICAL BIDOMAIN,CHMTRIX,CONTINU,FIRST_A,FIRSTITER,
     '  UPDATE_MATRIX,X_INIT
      SAVE FIRSTITER,NSTEP,TIME_START
      
      CALL ENTERS('MARCH7',*9999)
     
      DT=TINCR
      X_INIT=.FALSE.

! Check for Bidomain solution
      IF(KTYP32.EQ.1) THEN !monodomain
        BIDOMAIN=.FALSE.
      ELSE IF(KTYP32.EQ.2) THEN !bidomain
        BIDOMAIN=.TRUE.
      ELSE
        CALL ASSERT(.FALSE.,'>>Invalid type in KTYP32 (mono/bidomain)',
     '    ERROR,*9999)
      ENDIF

! Get class information
      nxc=NXLIST(1)
      CALL NX_LOC(NX_INQUIRE,nxc,nx_trans,NX_SOLVE,ERROR,*9999)
      CALL ASSERT(nx_trans.GT.0,'>>No nx defined for this solve class',
     '  ERROR,*9999)
      IF(BIDOMAIN) THEN
        CALL ASSERT(NXLIST(0).GE.2,'>>Error - bidomain needs 2 classes',
     '    ERROR,*9999)
        nxc=NXLIST(2)
        CALL NX_LOC(NX_INQUIRE,nxc,nx_ext,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx_ext.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)
      ELSE
        nx_ext=0
      ENDIF
      IF(NXLIST(0).GE.3) THEN
        nxc=NXLIST(3)
        CALL NX_LOC(NX_INQUIRE,nxc,nx_upd,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx_upd.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)
      ELSE
        nx_upd=0
      ENDIF
      IF(NXLIST(0).GE.4) THEN
        nxc=NXLIST(4)
        CALL NX_LOC(NX_INQUIRE,nxc,nx_torso,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx_torso.GT.0,
     '    '>>No nx defined for this solve class',ERROR,*9999)
      ELSE
        nx_torso=0
      ENDIF

! Get indices for AQ auxiliary parameters
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE1,maqp1t0,MAQ_START,
     '  ERROR,*9999)
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE1,maqp1t1,MAQ_STOP,
     '  ERROR,*9999)
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE1,maqp1i,MAQ_CURRENT,
     '  ERROR,*9999)
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE2,maqp2t0,MAQ_START,
     '  ERROR,*9999)
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE2,maqp2t1,MAQ_STOP,
     '  ERROR,*9999)
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE2,maqp2i,MAQ_CURRENT,
     '  ERROR,*9999)

! Get indicies for YQ (niq's)
      nr=NRLIST(1)
      IF(ITYP3(nr,nx_trans).EQ.1) THEN
C       Cubic w/o recovery
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqV,NIQ_V,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqBNDRY,NIQ_BNDRY,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqOLDSOLN,NIQ_OLDSOLN,
     '    ERROR,*9999)
      ELSE IF(ITYP3(nr,nx_trans).EQ.2) THEN
C       FitzHugh-Nagumo
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqV,NIQ_V,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqBNDRY,NIQ_BNDRY,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqRECOV,NIQ_RECOV,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqCALCIUM,NIQ_CALCIUM,
     '    ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqSAC,NIQ_SAC,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqOLDSOLN,NIQ_OLDSOLN,
     '    ERROR,*9999)
      ELSE IF(ITYP3(nr,nx_trans).EQ.3) THEN
C       van Capelle-Durrer
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqV,NIQ_V,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqBNDRY,NIQ_BNDRY,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqRECOV,NIQ_RECOV,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqCALCIUM,NIQ_CALCIUM,
     '    ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqSAC,NIQ_SAC,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqDV,NIQ_DV,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqDRECOV,NIQ_DRECOV,
     '    ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqOLDSOLN,NIQ_OLDSOLN,
     '    ERROR,*9999)
      ELSE IF(ITYP3(nr,nx_trans).EQ.4) THEN
C       Beeler-Reuter
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqV,NIQ_V,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqBNDRY,NIQ_BNDRY,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqX,NIQ_X,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqCAI,NIQ_CAI,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqMM,NIQ_M,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqH,NIQ_H,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqJ,NIQ_J,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqD,NIQ_D,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqF,NIQ_F,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqOLDSOLN,NIQ_OLDSOLN,
     '    ERROR,*9999)
      ELSE IF(ITYP3(nr,nx_trans).EQ.5) THEN
C       Jafri-Rice-Winslow
        CALL ASSERT(.FALSE.,'>>Not implemented',ERROR,*9999)
      ELSE IF(ITYP3(nr,nx_trans).EQ.6) THEN
C       Luo-Rudy
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqV,NIQ_V,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqBNDRY,NIQ_BNDRY,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqX,NIQ_X,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqCAI,NIQ_CAI,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqM,NIQ_M,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqH,NIQ_H,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqJ,NIQ_J,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqD,NIQ_D,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqF,NIQ_F,ERROR,*9999)
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqOLDSOLN,NIQ_OLDSOLN,
     '    ERROR,*9999)            
      ENDIF !Model type

! Create/update matricies if necessary and set start time
      IF(RESTART) THEN
        UPDATE_MATRIX=.FALSE.
        FIRST_A=.FALSE.
        IF(.NOT.ITER8) THEN
          T=T_RESTART
          IF(UP_GRID_MATERIAL.OR.UP_GRID_TENSOR) THEN
            CALL ASSEMBLE10(ISC_GKK,ISR_GKK,NEELEM,NQGP,NQGP_PIVOT,
     '        NQGW,NQS,NQXI,NRLIST,NWQ,nx_ext,nx_trans,nx_upd,
     '        NXQ,CQ,GCHQ,GUQ,GKK,PROPQ,BIDOMAIN,FIRST_A,FIXQ,.TRUE.,
     '        .TRUE.,UPDATE_MATRIX,ERROR,*9999)
            CHMTRIX=.TRUE.
          ELSE
            CHMTRIX=.FALSE.
          ENDIF
        ENDIF
      ELSE
        IF(.NOT.ITER8) THEN
          T=TSTART
          NSTEP=0
          CHMTRIX=.TRUE.
          FIRSTITER=.TRUE.
          CALL ASSEMBLE10(ISC_GKK,ISR_GKK,NEELEM,NQGP,NQGP_PIVOT,
     '      NQGW,NQS,NQXI,NRLIST,NWQ,nx_ext,nx_trans,nx_upd,
     '      NXQ,CQ,GCHQ,GUQ,GKK,PROPQ,BIDOMAIN,FIRST_A,FIXQ,.TRUE.,
     '      .FALSE.,UPDATE_MATRIX,ERROR,*9999)
        ENDIF
      ENDIF

      IF(.NOT.ITER8) THEN
        CONTINU=.TRUE.
        CALL CPU_TIMER(CPU_USER,TIME_START)

! Main time integration loop
        DO WHILE(CONTINU)
          NSTEP=NSTEP+1

! Update the ionic current calculations        
          CALL IONIC_CURRENT(CQ,IION,niqCAI,niqCALCIUM,niqD,niqDRECOV,
     '      niqDV,niqF,niqH,niqJ,niqMM,niqOLDSOLN,niqRECOV,niqSAC,niqV,
     '      niqX,NRLIST,nx_trans,YQ(1,1,1,nx_trans),ERROR,*9999)

! Generate the RHS vector for transmembrane potential solution
          CALL GEN_INT_RHS(AQ,CQ,IION,maqp1t0,maqp1t1,maqp1i,maqp2t0,
     '      maqp2t1,maqp2i,niqBNDRY,niqSAC,niqV,NQGP,NQGP_PIVOT,
     '      NQGW,NRLIST,NWQ,nx_ext,nx_trans,RHS,T,YQ,FIXQ,ERROR,*9999)

! Solve the system for the transmembrane potential
          CALL SOLVE_SYSTEM(ISC_GKK(1,nx_trans),ISR_GKK(1,nx_trans),
     '      1,NQT,NQT,NQT,NQT,NZZT(1,NRLIST(1),nx_trans),NZ_GKK_M,
     '      NISC_GKKM,NISR_GKKM,IWRIT4(1,nx_trans),
     '      PRECON_CODE(nx_trans),SOLVEROPTION(nx_trans),
     '      SPARSEGKK(nx_trans),SPARSEGKK(nx_trans),1,GKK(1,nx_trans),
     '      RHS,SOLNINT,FIRST_A,UPDATE_MATRIX,X_INIT,nx_trans,
     '      ERROR,*9999)

          IF(ITYP3(NRLIST(1),nx_trans).EQ.3.AND.KTYP33.EQ.2) THEN !VCDC
            DO nq=1,NQT
              YQ(nq,niqDV,1,nx_trans)=SOLNINT(nq)-
     '          YQ(nq,niqV,1,nx_trans)
            ENDDO
          ENDIF !VCDC

          CALL DCOPY(NQT,SOLNINT,1,YQ(1,niqV,1,nx_trans),1)

          IF(ITYP3(NRLIST(1),nx_trans).EQ.3.AND.KTYP33.EQ.2) THEN !VCDC
            IF(DOP) THEN
              WRITE(OP_STRING,
     '          '('' Warning: Potential maybe cutoff at 150mV'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF !dop
            DO nq=1,NQT
              IF(YQ(nq,niqV,1,nx_trans).GT.50.d0) THEN
                YQ(nq,niqV,1,nx_trans)=50.d0 !mV cutoff
              ENDIF
            ENDDO !nq
          ENDIF !VCDC


          IF(BIDOMAIN) THEN
            IF(CHMTRIX) THEN
              IF(.NOT.RESTART) FIRST_A=.TRUE.
              UPDATE_MATRIX=.TRUE.
              CHMTRIX=.FALSE.
            ENDIF

! Generate the RHS vector for extracellular potential solution
            CALL GEN_EXT_RHS(niqV,niqBNDRY,NQGP,NQGP_PIVOT,NQGW,NRLIST,
     '        NWQ,nx_ext,nx_trans,PROPQ,RHS,YQ,CQ,FIXQ(1,1,nx_ext),
     '        ERROR,*9999)

! Solve the system for the extracellular potential
            CALL SOLVE_SYSTEM(ISC_GKK(1,nx_ext),ISR_GKK(1,nx_ext),
     '        1,NQT,NQT,NQT,NQT,NZZT(1,NRLIST(1),nx_ext),NZ_GKK_M,
     '        NISC_GKKM,NISR_GKKM,IWRIT4(1,nx_ext),PRECON_CODE(nx_ext),
     '        SOLVEROPTION(nx_ext),SPARSEGKK(nx_ext),SPARSEGKK(nx_ext),
     '        1,GKK(1,nx_ext),RHS,SOLNEXT,FIRST_A,UPDATE_MATRIX,
     '        X_INIT,nx_ext,ERROR,*9999)

            CALL DCOPY(NQT,SOLNEXT,1,YQ(1,niqV,1,nx_ext),1)

          ENDIF !bidomain

          UPDATE_MATRIX=.FALSE.
          T=T+DT
          IF(T.GE.TFINISH) CONTINU=.FALSE.
          IF(DOP) THEN
            WRITE(OP_STRING,'(''Solution: '',I8)') NSTEP
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            DO nq=1,NQT
              WRITE(OP_STRING,'(I6,2F12.6)') nq,SOLNINT(nq),SOLNEXT(nq)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDIF
        ENDDO !time
! End of main time integration loop

        T_RESTART=T
        WRITE(OP_STRING,'('' Current time is '',F13.6,''ms'')') T
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Current step is '',F13.6,''ms'')') DT
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)
        WRITE(OP_STRING,'(1X,I6,'' total iterations : '',F8.2,'
     '    //'''s cpu'')') NSTEP,ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      ELSE !iterate

        CALL ASSERT(NXLIST(0).GE.4,
     '    '>>You need 4 classes for coupled bidomain',ERROR,*9999)

        IF(FIRSTITER) THEN
          FIRSTITER=.FALSE.
          UPDATE_MATRIX=.TRUE.
          FIRST_A=.TRUE.
          DO nq=1,NQT
            FIXQ(nq,3,nx_upd)=.TRUE.
          ENDDO
        ELSE
          UPDATE_MATRIX=.FALSE.
          FIRST_A=.FALSE.
        ENDIF

        NRTEMP(0)=1
        NRTEMP(1)=NRLIST(1)

! Generate the RHS vector for extracellular potential solution
        CALL GEN_EXT_RHS(niqV,niqBNDRY,NQGP,NQGP_PIVOT,NQGW,NRTEMP,NWQ,
     '    nx_ext,nx_trans,PROPQ,RHS,YQ,CQ,FIXQ(1,1,nx_upd),ERROR,*9999)

! Update RHS from YP fluxes
        CALL GEN_GRID_BEM_RHS(IBT,IDO,INP,ISC_GKK(1,nx_upd),
     '    ISR_GKK(1,nx_upd),NBH,NBJ,NEELEM,NENQ,NHE,NHP,NKE,NKH,NPF,
     '    NP_INTERFACE,NPNE,NPNODE,NQNE,
     '    NQS,NQXI,NVHE,NVHP,NYNE,NYNP,
     '    NRLIST,NW,NWQ,nx_torso,nx_upd,NXQ,CURVCORRECT,
     '    GKK(1,nx_upd),RHS,
     '    SE,YP,ZA,ZE,ZP,UPDATE_MATRIX,ERROR,*9999)

        CALL SOLVE_SYSTEM(ISC_GKK(1,nx_upd),ISR_GKK(1,nx_upd),1,NQT,NQT,
     '    NQT,NQT,NZZT(1,NRLIST(1),nx_upd),NZ_GKK_M,NISC_GKKM,NISR_GKKM,
     '    IWRIT4(1,nx_upd),PRECON_CODE(nx_upd),SOLVEROPTION(nx_upd),
     '    SPARSEGKK(nx_upd),SPARSEGKK(nx_upd),1,GKK(1,nx_upd),RHS,
     '    SOLNEXT,FIRST_A,UPDATE_MATRIX,X_INIT,nx_upd,ERROR,*9999)

        CALL DCOPY(NQT,SOLNEXT,1,YQ(1,niqV,1,nx_ext),1)

      ENDIF

      CALL EXITS('MARCH7')
      RETURN
 9999 CALL ERRORS('MARCH7',ERROR)
      CALL EXITS('MARCH7')
      RETURN 1
      END


C FE19
C=====

      SUBROUTINE BR_ION(Y,IONIC_PARAMS,F)

C#### Subroutine: BR_ION
C###  Description:
C###    BR_ION is the Beeler-Reuter model, with the option of sodium
C###    kinetics determined by Ebihara-Johnson or Drouhard-Roberge
C###    according to KTYP33
                    
      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:ktyp30.cmn'
!     Parameter List
      REAL*8 F(*),IONIC_PARAMS(12),Y(*)
!     Local Variables
      REAL*8 GNa,GNaC,Gs,INa,Is,Ix1,IK1,Vs,VNa,PHIM,m,h,j,d,f1,x1,Cai

C     IONIC_PARAMS(1) !resting potential
C     IONIC_PARAMS(2) !sodium equilibrium (reversal) potential
C     IONIC_PARAMS(3) !sodium conductance
C     IONIC_PARAMS(4) !steady-state sodium conductance
C     IONIC_PARAMS(5) !slow current conductance

      PHIM=Y(1)
      m=Y(2)
      h=Y(3)
      j=Y(4)
      d=Y(5)
      f1=Y(6)
      x1=Y(7)
      Cai=Y(8)

      VNa=IONIC_PARAMS(2) 
      GNa=IONIC_PARAMS(3)
      GNaC=IONIC_PARAMS(4) 
      Gs=IONIC_PARAMS(5) 

      INa=0.d0
      IF(KTYP33.EQ.1) THEN !Beeler-Reuter sodium model
        INa = GNa*m*m*m*h*j + GNaC
      ELSE IF(KTYP33.EQ.2.OR.KTYP33.EQ.3) THEN 
C       Ebihara-Johnson or Drouhard-Roberge sodium model
        INa = GNa*m*m*m*h
      ENDIF !sodium model
      INa = INa*(PHIM-VNa)

      Vs = -82.3d0 - 13.0287d0*DLOG(Cai)
      Is = Gs*d*f1*(PHIM-Vs)

C     1.d-2 correction factor in the next two equations is for the 
C     change from cm^2 to mm^2

      Ix1 = 0.8d0*x1*(DEXP(0.04d0*(PHIM+77.d0))-1.d0)/
     '  (DEXP(0.04d0*(PHIM+35.d0)))!*1.d-2

      IK1 = 0.35d0 * (4.d0*(DEXP(0.04d0*(PHIM+85.d0))-1.d0) /
     '  (DEXP(0.08d0*(PHIM+53.d0))+DEXP(0.04d0*(PHIM+53.d0))) +
     '  0.2d0*(PHIM+23.d0)/(1.d0-DEXP(-0.04d0*(PHIM+23.d0))))!*1.d-2
      
      if(dop) then
        write(*,'('' BR: INa, Is, Ix1, IK1'',4E15.5)') INa,Is,Ix1,IK1
      endif
      
      F(1)=(INa+Is+Ix1+IK1)*1.0d-6 !correction factor for [mV] & [mA]

      RETURN
      END          


      SUBROUTINE BR_RATES(Y,IONIC_PARAMS,F)

C#### Subroutine: BR_RATES
C###  Description:
C###    Computes the rate constants for the Beeler-Reuter model
                    
      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b12.cmn'
      INCLUDE 'cmiss$reference:ktyp30.cmn'
!     Parameter List
      REAL*8 F(*),IONIC_PARAMS(12),Y(*)
!     Local Variables
      REAL*8 Gs,Is,Vs,V,m,h,j,d,f1,x1,Cai,
     '  Am,Bm,Ah,Bh,Aj,Bj,Ad,Bd,Af,Bf,Ax1,Bx1,Dm,Dh,Dj,Dd,Df,DCai,Dx1

C     IONIC_PARAMS(1) !resting potential
C     IONIC_PARAMS(2) !sodium equilibrium (reversal) potential
C     IONIC_PARAMS(3) !sodium conductance
C     IONIC_PARAMS(4) !steady-state sodium conductance
C     IONIC_PARAMS(5) !slow current conductance

      V=Y(1)
      m=Y(2)
      h=Y(3)
      j=Y(4)
      d=Y(5)
      f1=Y(6)
      x1=Y(7)
      Cai=Y(8)

      Gs=IONIC_PARAMS(5) 

      IF(KTYP33.EQ.1) THEN !Beeler-Reuter sodium model
        Am = -1.d0*(V+47.d0)/(DEXP(-0.1d0*(V+47.d0))-1.d0)
        Bm = 40.d0*DEXP(-0.056d0*(V+72.d0))
        Ah = 0.126d0*DEXP(-0.25d0*(V+77.d0))
        Bh = 1.7d0/(DEXP(-.082d0*(V+22.5d0))+1.d0)
        Aj = 0.055d0*DEXP(-0.25d0*(V+78.d0))/
     '    (DEXP(-0.2d0*(V+78.d0))+1.d0)
        Bj = 0.3d0/(DEXP(-0.1d0*(V+32.d0))+1.d0)
      ELSE IF(KTYP33.EQ.2) THEN !Ebihara-Johnson sodium model
        Am = (0.32d0*(V+47.13d0))/(1.d0-DEXP(-V-47.13d0))
        Bm = 0.08d0*DEXP(-V/11.d0)
        IF(V.GE.-40.d0) THEN
          Ah = 0.d0
          Bh = 1.d0/(0.13d0*(DEXP((V+10.66d0)/(-11.1d0))+1.d0))
        ELSE
          Ah = 0.135d0*DEXP((-80.d0-V)/6.8d0)
          Bh = 3.56d0*DEXP(0.079d0*V)+3.1d5*DEXP(0.35d0*V)
        ENDIF
        Aj = 0.d0
        Bj = 0.d0
      ELSE IF(KTYP33.EQ.3) THEN !Drouhard-Roberge sodium model
        Am = 0.9d0*(V+42.65d0)/(1.d0-DEXP(-0.22d0*(V+42.65d0)))
        Bm = 1.437d0*DEXP(-0.085d0*(V+39.75d0))
        Ah = 0.1d0*DEXP(-0.193d0*(V+79.65d0))
        Bh = 1.7d0/(1.d0+DEXP(-0.095d0*(V+20.5d0)))
        Aj = 0.d0
        Bj = 0.d0
      ENDIF !sodium model
      Dm=Am*(1.d0-m)-Bm*m
      Dh=Ah*(1.d0-h)-Bh*h
      Dj=Aj*(1.d0-j)-Bj*j

      Vs = -82.3d0 - 13.0287d0*DLOG(Cai)
      Is = Gs*d*f1*(V-Vs)
      DCai = -1.d-7*Is + 0.07d0*(1.d-7-Cai)

      Ad = 0.095d0*DEXP(-0.01d0*(V-5.d0))/
     '  (1.d0+DEXP(-0.072d0*(V-5.d0)))
      Bd = 0.07d0*DEXP(-(V+44.d0)/59.d0)/
     '  (1.d0+DEXP(0.05d0*(V+44.d0)))
      Af = 0.012d0*DEXP(-0.008d0*(V+28.d0))/
     '  (1.d0+DEXP(0.15d0*(V+28.d0)))
      Bf = 0.0065d0*DEXP(-0.02d0*(V+30.d0))/
     '  (1.d0+DEXP(-0.2d0*(V+30.d0)))
      Dd=Ad*(1.d0-d)-Bd*d
      Df=Af*(1.d0-f1)-Bf*f1

      Ax1 = 5.d-4*DEXP(-(V+50.d0)/12.1d0)/(1.d0+DEXP((V+50.d0)/17.5d0))
      Bx1 = 0.0013d0*DEXP(-0.06d0*(V+20.d0))/
     '  (1.d0+DEXP(-0.04d0*(V+20.d0)))
      Dx1=Ax1*(1.d0-x1)-Bx1*x1

C      m=m+DT*Dm*1.d-1 !corrections of [ms]
C      h=h+DT*Dh*1.d-1
C      j=j+DT*Dj*1.d-1
C      d=d+DT*Dd*1.d-1
C      f=f+DT*Df*1.d-1
C      x1=x1+DT*Dx1*1.d-1
C      Cai=Cai+DT*DCai*1.d-1

      F(2)=Dm*1.d-1
      F(3)=Dh*1.d-1
      F(4)=Dj*1.d-1
      F(5)=Dd*1.d-1
      F(6)=Df*1.d-1
      F(7)=Dx1*1.d-1
      F(8)=DCai*1.d-1

      RETURN      
      END


      REAL*8 FUNCTION DN_ION(PHIM,m,h,IONIC_PARAMS)

C#### Function: DN_ION
C###  Description:
C###    DN_ION is the diFrancesco-Noble model for ionic current.
                    
      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
!     Parameter List
      REAL*8 h,IONIC_PARAMS(12),m,PHIM
!     Local Variables
      REAL*8 GL,GNa,IL,INa,VL,VNa

      VL= IONIC_PARAMS(1) !equal to resting potential
      GL= IONIC_PARAMS(2) !equal to reciprocal of membrane resistance
      GNa=IONIC_PARAMS(3) !sodium conductance
      VNa=IONIC_PARAMS(4) !sodium equilibrium (reversal) potential
      
C **  Corrections applied for GNa and GL
      INa=GNa*m*m*m*h*(PHIM-VNa)*1.0d-3
      IL= GL*(PHIM-VL)*1.0d-6           

      if(dop) then
        write(*,'(''DN: INa, IL'',2E15.5)') INa, IL
      endif
      
      DN_ION=INa+IL

      RETURN
      END


      REAL*8 FUNCTION FHN_ION(PHIM,RECOV,IONIC_PARAMS)

C#### Function: FHN_ION
C###  Type: REAL*8
C###  Description:
C###    FHN_ION is the FitzHugh-Nagumo ionic equation for 
C###    bidomain solution.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:ktyp30.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
!     Parameter List
      REAL*8 IONIC_PARAMS(12),PHIM,RECOV
!     Local Variables
      REAL*8 DIFF,THRESHOLD,PHI,PLATEAU,REST,RATE,DECAY,FHN
      CHARACTER ERROR*20

C *** IONIC_PARAMS(1) is rest potential
C *** IONIC_PARAMS(2) is plateau potential
C *** IONIC_PARAMS(3) is threshold potential
C *** IONIC_PARAMS(4) is excitation rate constant
C *** IONIC_PARAMS(5) is excitation decay constant
C ***  (6...) are recovery constants

      REST=IONIC_PARAMS(1)
      PLATEAU=IONIC_PARAMS(2)
      THRESHOLD=IONIC_PARAMS(3)
      RATE=IONIC_PARAMS(4)
      DECAY=IONIC_PARAMS(5)
      
      IF(DOP) THEN
        WRITE(OP_STRING,'('' REST='',D12.4)') REST
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' PLATEAU='',D12.4)') PLATEAU
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' THRESHOLD='',D12.4)') THRESHOLD
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' RATE='',D12.4)') RATE
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' DECAY='',D12.4)') DECAY
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

C     Normalise wrt resting and plateau potentials
      DIFF = PLATEAU-REST

      PHI = (PHIM-REST)/DIFF
      THRESHOLD=(THRESHOLD-REST)/DIFF

c      FHN = RATE/DIFF/DIFF * (PHIM-REST) * (PHIM-THRESHOLD) * 
c     '  (PHIM-PLATEAU)
      FHN = RATE*PHI*(PHI-THRESHOLD)*(PHI-1.d0)

      IF(DOP) THEN
        WRITE(OP_STRING,'('' DIFF='',D12.4)') DIFF
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' THRESHOLD='',D12.4)') THRESHOLD
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' FHN='',D12.4)') FHN
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      
      IF(KTYP33.EQ.1) THEN !Standard FHN
        FHN=FHN+DECAY*RECOV
      ELSE !Roger's FHN or Panfilov FHN
C        FHN=FHN+DECAY*(PHIM-REST)*RECOV
        FHN=FHN+DECAY*PHI*RECOV
      ENDIF

      FHN_ION = FHN

      IF(DOP) THEN
        WRITE(OP_STRING,'('' FHN_ION='',D12.4)') FHN_ION
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

 9999 RETURN
      END


      REAL*8 FUNCTION LR_ION(PHIM,m,h,j,d,f,x,Cai,IONIC_PARAMS)

C#### Function: LR_ION
C###  Description:
C###    LR_ION is the (first) Luo-Rudy model

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
!     Parameter List
      REAL*8 IONIC_PARAMS(12),PHIM,m,h,j,d,f,x,Cai
!     Local Variables
      REAL*8 aK1,bK1,GK,GK1,GNa,Gs,INa,Is,IK,IK1,IKp,Ib,K1_inf,Kp,
     '  VK,VK1,Vs,VNa,Xi

C     IONIC_PARAMS(1)  !resting potential
C     IONIC_PARAMS(2)  !sodium equilibrium (reversal) potential
C     IONIC_PARAMS(3)  !sodium conductance
C     IONIC_PARAMS(4)  !slow current conductance
C     IONIC_PARAMS(9)  !potassium conductance
C     IONIC_PARAMS(10) !time-indep potassium conductance
C     IONIC_PARAMS(11) !potassium reversal potential
C     IONIC_PARAMS(12) !time-indep potassium reversal potential

      VNa=IONIC_PARAMS(2) 
      GNa=IONIC_PARAMS(3) 
      Gs=IONIC_PARAMS(4)  
      GK=IONIC_PARAMS(9) 
      GK1=IONIC_PARAMS(10)
      VK=IONIC_PARAMS(11) 
      VK1=IONIC_PARAMS(12)

      INa = GNa*m*m*m*h*j*(PHIM-VNa)

      Vs = -82.3d0 - 13.0287d0*DLOG(Cai)
      Is = Gs*d*f*(PHIM-Vs)

      IF(PHIM.GT.-100.d0) THEN
        Xi = 2.837d0*(DEXP(0.04d0*(PHIM+77.d0))-1.d0)/
     '    ((PHIM+77.d0)*DEXP(0.04d0*(PHIM+35.d0)))
      ELSE
        Xi=1.d0
      ENDIF
      IK = gK*x*Xi*(PHIM-VK)

      aK1 = 1.02d0/(1.d0+DEXP(0.2385d0*(PHIM-VK1-59.215d0)))
      bK1 = (0.49124d0*DEXP(0.08032d0*(PHIM-VK1+5.476d0))
     '  +DEXP(0.06175d0*(PHIM-VK1-594.31d0)))/
     '  (1.d0+DEXP(-0.5143d0*(PHIM-VK1+4.753d0)))
      K1_inf=aK1/(aK1+bK1)
      IK1 = gK1*K1_inf*(PHIM-VK1)

      Kp = 1.d0/(1.d0+DEXP((7.488d0-PHIM)/5.98d0))
      IKp = 0.0183d0*Kp*(PHIM-VK1)

      Ib = 0.03921d0*(PHIM+59.87d0)
      
      if(dop) then
        write(*,'('' LR: INa,Is,IK,IK1,IKp,Ib'',/6E10.5)') 
     '    INa,Is,IK,IK1,IKp,Ib
      endif
      
      LR_ION=(INa+Is+IK+IK1+IKp+Ib)*1.0d-6 !corr. factor for [mV] & [mA]

      RETURN
      END


      SUBROUTINE LR_RATES(V,m,h,j,d,f,x,Cai,IONIC_PARAMS)

C#### Subroutine: LR_RATES
C###  Description:
C###    Computes the rate constants for the Luo-Rudy model
                    
      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b12.cmn'
!     Parameter List
      REAL*8 IONIC_PARAMS(12),V,m,h,j,d,f,x,Cai
!     Local Variables
      REAL*8 Gs,Is,Vs,
     '  Am,Bm,Ah,Bh,Aj,Bj,Ad,Bd,Af,Bf,Ax,Bx,Dm,Dh,Dj,Dd,Df,DCai,Dx

C      IONIC_PARAMS(1) !resting potential
C      IONIC_PARAMS(2) !sodium equilibrium (reversal) potential
C      IONIC_PARAMS(3) !sodium conductance
C      IONIC_PARAMS(4) !slow current conductance

      Gs=IONIC_PARAMS(4) 

      Am = (0.32d0*(V+47.13d0))/(1.d0-DEXP(-V-47.13d0))
      Bm = 0.08d0*DEXP(-V/11.d0)
      IF(V.GE.-40.d0) THEN
        Ah = 0.d0
        Aj = 0.d0
        Bh = 1.d0/(0.13d0*(DEXP((V+10.66d0)/(-11.1d0))+1.d0))
        Bj = 0.3d0*DEXP(-2.535d-7*V)/(1.d0+DEXP(-0.1d0*(V+32.d0)))
       ELSE
        Ah = 0.135d0*DEXP((-80.d0-V)/6.8d0)
        Aj = (-1.2714d5*DEXP(0.2444d0*V)
     '    -3.474d-5*DEXP(-0.04391d0*V))*(V+37.78d0)/
     '    (1.d0+DEXP(0.311d0*(V+79.23d0)))
        Bh = 3.56d0*DEXP(0.079d0*V)+3.1d5*DEXP(0.35d0*V)
        Bj = 0.1212d0*DEXP(-0.01052d0*V)/
     '    (1.d0+DEXP(-0.1378d0*(V+40.14d0)))
      ENDIF
      Dm=Am*(1.d0-m)-Bm*m
      Dh=Ah*(1.d0-h)-Bh*h
      Dj=Aj*(1.d0-j)-Bj*j

      Vs = -82.3d0 - 13.0287d0*DLOG(Cai)
      Is = Gs*d*f*(V-Vs)
      DCai = -1.d-7*Is + 0.07d0*(1.d-7-Cai)

      Ad = 0.095d0*DEXP(-0.01d0*(V-5.d0))/
     '  (1.d0+DEXP(-0.072d0*(V-5.d0)))
      Bd = 0.07d0*DEXP(-(V+44.d0)/59.d0)/
     '  (1.d0+DEXP(0.05d0*(V+44.d0)))
      Af = 0.012d0*DEXP(-0.008d0*(V+28.d0))/
     '  (1.d0+DEXP(0.15d0*(V+28.d0)))
      Bf = 0.0065d0*DEXP(-0.02d0*(V+30.d0))/
     '  (1.d0+DEXP(-0.2d0*(V+30.d0)))
      Dd=Ad*(1.d0-d)-Bd*d
      Df=Af*(1.d0-f)-Bf*f

      Ax = 5.d-4*DEXP(-(V+50.d0)/12.1d0)/(1.d0+DEXP((V+50.d0)/17.5d0))
      Bx = 0.0013d0*DEXP(-0.06d0*(V+20.d0))/
     '  (1.d0+DEXP(-0.04d0*(V+20.d0)))
      Dx=Ax*(1.d0-x)-Bx*x

      m=m+DT*Dm*1.d-1 !corrections of [ms]
      h=h+DT*Dh*1.d-1
      j=j+DT*Dj*1.d-1
      d=d+DT*Dd*1.d-1
      f=f+DT*Df*1.d-1
      x=x+DT*Dx*1.d-1
      Cai=Cai+DT*DCai*1.d-1
      
      RETURN      
      END


      REAL*8 FUNCTION VCD_ION(PHIM,RECOV,D_PHIM)

C#### Function: VCD_ION
C###  Type: REAL*8
C###  Description:
C###    VCD_ION is the VanCapelle-Durrer ionic equation for bidomain 
C###    solution, with modifications from UCLA provided by Alan 
C###    Garfinkel. Original reference in fax dated March 8, 1994.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:ktyp30.cmn'
!     Parameter List
      REAL*8 D_PHIM,PHIM,RECOV
!     Local Variables
      REAL*8 F,I0,I1

      IF(PHIM.LT.-70.0d0) THEN
        I1=5.0d0+0.5d0*(PHIM+70.0d0)
      ELSE IF(phim.gt.0.0d0) THEN
        I1=6.0d0+0.425d0*PHIM     
      ELSE
        I1=5.0d0+((PHIM+70.0d0)/70.0d0)
      ENDIF
      IF(PHIM.LT.-74.3d0) THEN
        F=7.84d0+2.0d0*(PHIM+74.3d0)
      ELSE IF(PHIM.GT.-27.8d0) THEN
        F=-98.84d0+1.71d0*(PHIM+27.8d0)
      ELSE            
        F=((3.837854d-3*PHIM +0.584649d0)*PHIM +25.31834d0)*PHIM
     '    +235.6256d0
      ENDIF

      IF(KTYP33.EQ.2) THEN  !Calif. mods to VCD
        F=4.0d0*F
        IF(KTYP34.EQ.2) THEN    !"Ischemic" APD (112ms)
          IF(D_PHIM.LT.0) I1=2.0d0*I1
        ENDIF
      ENDIF
      
      I0=I1+F

      VCD_ION=(RECOV*I1+(1.0d0-RECOV)*I0) * 1.d-6

      RETURN
      END


C FE30
C=====

      SUBROUTINE BFRONT(NBJ,NEELEM,niqCAI,niqCALCIUM,niqD,niqDRECOV,
     '  niqDV,niqF,niqH,niqJ,niqMM,niqOLDSOLN,niqRECOV,niqSAC,niqV,
     '  niqX,NRLIST,NWQ,nx,NXQ,AQ,CQ,GCHQ,GUQ,PROPQ,T,YQ,ERROR,*)

C#### Subroutine: BFRONT
C###  Description:
C###    BFRONT calculates ionic current based cardiac activation 
C###    equations. It calculates activation pattern by using Bidomain 
C###    equations. The equation is modelled explicitly, with explicit 
C###    finite differences for diffusion calculations.

C**** ITYP2(nr,nx)=9 -- Monodomain/Bidomain model
C**** KTYP32 is 1,2 for monodomain/bidomain
C**** KTYP36 = 1 if using DTAR (Dynamic Tracking of the Active Region)

C**** NWQ defined in DEGRID:
C**** na dropped in this subroutine 
C**** NWQ(1,nq) is   0 if nq is internal
C***      "        mq1 if nq on external bdy (for calc no-flux b.c.)
C***  NWQ(2,nq) is   0 if nq is internal
C***      "        mq2 if nq on external bdy
C***  NWQ(4,nq) is -1 if grid point nq was active and is now inactive
C***                0  "    "     "    is not active     
C***     "          1  "    "     "    "  "   active
C***     "          2  "    "     "    "  surrounded by active points
C***  NWQ(4,0) contains number of _currently_ active points
C***  i.e. sum of NWQ(4,1..nqt)
C***  NWQ(5,nq) contains a list of active points for DTAR, with the
C***     total number of active points in NWQ(5,0)
C**** PHIMQ and PHIEQ contain the value of phi(m) and phi(e) at the
C**** grid points of the local quadratic element about nq.
C**** DPHIDX   contains phim,i+phie,i and similarly
C**** D2PHIDX2 contains phim,ij+phie,ij

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b12.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:grid00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:ktyp30.cmn'
      INCLUDE 'cmiss$reference:maqloc00.inc'
      INCLUDE 'cmiss$reference:time02.cmn'
      INCLUDE 'cmiss$reference:tol00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),niqCAI,niqCALCIUM,
     '  niqD,niqDRECOV,niqDV,niqF,niqH,niqJ,niqMM,niqOLDSOLN,niqRECOV,
     '  niqSAC,niqV,niqX,NRLIST(0:NRM),NWQ(6,0:NQM),
     '  NXQ(-NIM:NIM,0:4,0:NQM)
      REAL*8 AQ(NMAQM,NQM),CQ(NMM,NQM),GCHQ(3,NQM),GUQ(3,3,NQM),
     '  PROPQ(3,3,4,2,NQM),T,YQ(NYQM,NIQM,NAM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER atime,dpot,ib,IJ,IK,jb,maqp1t0,maqp1t1,maqp1i,maqp2t0,
     '  maqp2t1,maqp2i,mq,mqi,nbri,ni,nii,nij,nik,niq,NITB,nj,nk,
     '  nnq,nnq_min,nq,nr,nrr,NSOL_EXT,NSOL_INT,nx
      REAL*8 Am,CM,CUBIC_ION,DVMST,D_RECOV,DIFF,
     '  DPHIM,DPHIDX(3),D2PHIDX2(3,3),
     '  EPS,FHN_ION,
     '  IAPP,IION,IM,IPROP,LR_ION,NORM_PHI,PHIM,PHIEQ(-1:1,-1:1,-1:1),
     '  PHIM1,PHIM2,PHIMP,PHIMQ(-1:1,-1:1,-1:1),RATE,RECOV,SUM1,SUM2,
     '  T_CONST,V,VCD_ION,VSTAR,VTHRESH
      LOGICAL ERROR_FLAG,SAC
      SAVE NNQ_MIN

      CALL ENTERS('BFRONT',*9999)

      IF(ABS(T-T0).LT.(DT/10.d0)) THEN !First time thru
        nnq_min=1  !Start at beginning of array
      ENDIF

      SAC=.FALSE.
      IF(ITYP3(NRLIST(1),nx).EQ.2) THEN      !FHN
        IF(DABS(CQ(18,1)).GT.1.D-8) SAC=.TRUE.
      ELSEIF(ITYP3(NRLIST(1),nx).EQ.3) THEN  !VCD
        IF(DABS(CQ(13,1)).GT.1.D-8) SAC=.TRUE.
      ENDIF

      !Initialise local arrays
      DO ni=-1,1
        DO nj=-1,1
          DO nk=-1,1
            PHIMQ(nk,nj,ni)=0.d0
            PHIEQ(nk,nj,ni)=0.d0
          ENDDO
          D2PHIDX2(nj+2,ni+2)=0.d0
        ENDDO
        DPHIDX(ni+2)=0.d0
      ENDDO

      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE1,maqp1t0,MAQ_START,
     '  ERROR,*9999)
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE1,maqp1t1,MAQ_STOP,
     '  ERROR,*9999)
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE1,maqp1i,MAQ_CURRENT,
     '  ERROR,*9999)
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE2,maqp2t0,MAQ_START,
     '  ERROR,*9999)
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE2,maqp2t1,MAQ_STOP,
     '  ERROR,*9999)
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE2,maqp2i,MAQ_CURRENT,
     '  ERROR,*9999)
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_TIME,atime,MAQ_ACTIV_TIME,
     '  ERROR,*9999)
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_TIME,dpot,MAQ_M_DPOT,
     '  ERROR,*9999)

      IF(DOP) THEN
C$      CALL mp_setlock()
        WRITE(OP_STRING,'('' Allocated maqs'',6I2)') maqp1t0,maqp1t1,
     '    maqp1i,maqp2t0,maqp2t1,maqp2i
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C$      CALL mp_unsetlock()
      ENDIF      

      NSOL_INT=0
      NSOL_EXT=0
      
C***  For each g.p. nq, put the value of phi(m) and phi(e) to be used
C***  for this time step into YQ(nq,(2:3),3) - temporary storage
      DO nrr=1,NRLIST(0)
        nr=NRLIST(nrr) 
        NITB=NIT(NBJ(1,NEELEM(1,nr)))
        DO nq=NQR(1,nr),NQR(2,nr)          
          YQ(nq,niqOLDSOLN,1)=YQ(nq,niqV,1)
Cb          YQ(nq,3,3)=YQ(nq,2,1)
C***  Get applied current if a current has been defined for this grid
C***  point, and if the time is within the application times for that
C***  current pulse.  If so, then make sure it is in active list.
          IF(NWQ(4,nq).EQ.0) THEN !Not already active or prev. active
            IAPP=0.0d0

            IF((AQ(maqp1i,nq).GT.ZERO_TOL).AND.(T.GE.AQ(maqp1t0,nq))
     '        .AND.(T.LT.AQ(maqp1t1,nq))) THEN
              IAPP=AQ(maqp1i,nq)
            ELSE IF((AQ(maqp2i,nq).GT.ZERO_TOL).AND.(T.GE.
     '        AQ(maqp2t0,nq)).AND.(T.LT.AQ(maqp2t1,nq))) THEN
              IAPP=AQ(maqp2i,nq)
            ENDIF

            IF(SAC) THEN
              IF(YQ(nq,niqSAC,1).GT.1.d-8) THEN !inward SAC current
                IAPP=IAPP+YQ(nq,niqSAC,1)
                IF(DOP) THEN
C$                call mp_setlock()
                  write(*,'('' Inward SAC current at nq='',I5,'
     '              //''' is'',E12.3)') nq,IAPP
C$                call mp_unsetlock()
                ENDIF
              ENDIF
            ENDIF
            IF(IAPP.GT.1.D-8) THEN
              NWQ(4,nq)=1
              NWQ(5,0)=NWQ(5,0)+1
              NWQ(5,NWQ(5,0))=nq
              IF(DOP) THEN
C$              call mp_setlock()
                WRITE(OP_STRING,'(''  Adding nq : '',I6)') nq
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C$              call mp_unsetlock()
              ENDIF
            ENDIF
          ENDIF
        ENDDO !nq
      ENDDO !nr

      DO nrr=1,NRLIST(0) !Loop over regions
        nr=NRLIST(nrr)

C old MLB 18 August 1997
C        nxielem=NIT(NBJ(1,NEELEM(1,nr)))
C!Find basis function using Extended Lagrange 
C        nb_extended=1
C        DO nbb=1,NBT
C          IF((nxielem.NE.NIT(nb_extended)).OR.(NBC(nb_extended).NE.7))
C     '      THEN
C            nb_extended=nb_extended+1
C          ENDIF
C        ENDDO
C        NITB=NIT(nb_extended)

        NITB=NIT(NBJ(1,NEELEM(1,nr)))

C     Dynamic Tracking of the Active Region (DTAR)
C     Loop over all active points and make neighbours active if Vm is
C     greater than (resting potential + 1mV)
        IK=MAX(0,NITB-2)          !zero for 1,2-D, one for 3-D
        IJ=MIN(NITB-1,1)          !zero for 1-D, one for 2,3-D

        DO nnq=1,NWQ(5,0)
          nq=NWQ(5,nnq)
          VTHRESH=CQ(9,NQ)+1.d0 !CQ(9,nq) is the resting potential
          IF(NWQ(4,nq).EQ.1.AND.YQ(nq,niqOLDSOLN,1).GT.VTHRESH) THEN
            if(dop) then
              WRITE(OP_STRING,'(''  Scanning nq : '',I6)') nq
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            endif
C     Activate all neighbouring grid points
            DO nik=-IK,IK
              DO nij=-IJ,IJ
                DO nii=-1,1
                  IF(NXQ(nii,0,nq).LE.1) THEN
                    mq=NXQ(nii,1,NXQ(nij*2,1,NXQ(nik*3,1,nq)))
                    IF(mq.GT.0.AND.NWQ(4,mq).EQ.0) THEN 
                      !Not already active
                      if(dop) then
                        WRITE(OP_STRING,'(''    Adding mq : '',I6)') mq
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      endif
                      NWQ(4,mq)=1
                      NWQ(5,0)=NWQ(5,0)+1
                      NWQ(5,NWQ(5,0))=mq
                    ENDIF
                  ELSE
                    DO nbri=1,NXQ(nii,0,nq) !branches in xi1 from nq
                      mqi=NXQ(nii,nbri,NXQ(nij*2,1,NXQ(nik*3,1,nq)))
                      IF(mqi.GT.0.AND.NWQ(4,mqi).EQ.0) THEN
                        !if point exists & not already active
                        if(dop) then
                          WRITE(OP_STRING,
     '                      '(''    Adding mqi : '',I6)') mqi
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                        endif
                        NWQ(4,mqi)=1
                        NWQ(5,0)=NWQ(5,0)+1
                        NWQ(5,NWQ(5,0))=mqi
                        IF(NWQ(6,mqi).GT.0) THEN
                          NWQ(4,NWQ(6,mqi))=1
                          NWQ(5,0)=NWQ(5,0)+1
                          NWQ(5,NWQ(5,0))=NWQ(6,mqi)
                        ENDIF !coupled nodes
                      ENDIF
                    ENDDO !nbri
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
            NWQ(4,nq)=2
          ENDIF
        ENDDO

C       Loop over all active grid points & solve
        ERROR_FLAG=.FALSE.
C$DOACROSS local(nnq,mq,mqi,nbri,nq,nik,nij,nii,PHIMQ,PHIEQ,DPHIDX,
C$&              ib,jb,
C$&              D2PHIDX2,IPROP,ni,nj,SUM1,SUM2,IAPP,AM,IM,PHIM,IION,
C$&              CM,DPHIM,PHIMP,V,RECOV,DIFF,NORM_PHI,D_RECOV,EPS,RATE,
C$&              VSTAR,T_CONST,niq),
C$&        share(SAC,NWQ,NQR,NXQ,YQ,ERROR_FLAG),
C$&        reduction(NSOL_INT)
        DO nnq=nnq_min,NWQ(5,0)
          IF(ERROR_FLAG) GOTO 200
          nq=NWQ(5,nnq)
          IF((nq.GE.NQR(1,nr)).AND.(nq.LE.NQR(2,nr))) THEN
C***  If we have an internal g.p. continue (and still active)
            IF(NWQ(1,nq).EQ.0.AND.NWQ(4,nq).GT.0) THEN !Not on boundary
C***  Formulate a local quadratic element about nq defined by UMQ,
C***  storing the value of phi(m) at each node point mq.
              NSOL_INT=NSOL_INT+1

C MLB 28 August 1996 altered to the following 2 loops

C Initialise Phi's to zero.  Must do it here inside nq loop.
            DO nik=-IK,IK
              DO nij=-IJ,IJ
                DO nii=-1,1
                  PHIMQ(nii,nij,nik)=0.0d0
                  PHIEQ(nii,nij,nik)=0.0d0
                ENDDO !nii
              ENDDO !nij
            ENDDO !nik

            DO nik=-IK,IK                                
              DO nij=-IJ,IJ
                DO nii=-1,1
                  IF(NXQ(nii,0,nq).LE.1) THEN
                    mq=NXQ(nii,1,NXQ(nij*2,1,NXQ(nik*3,1,nq)))
C***  Value of phi(m) and phi(e) at MQ
                    IF(mq.GT.0) THEN
                      PHIMQ(nii,nij,nik)=YQ(mq,niqOLDSOLN,1)
Cb                      PHIEQ(nii,nij,nik)=YQ(mq,3,3)
                      PHIEQ(nii,nij,nik)=0.0d0
                    ELSE
                      PHIMQ(nii,nij,nik)=0.0d0
                      PHIEQ(nii,nij,nik)=0.0d0
                    ENDIF
                  ELSE
                    DO nbri=1,NXQ(nii,0,nq)
                      mqi=NXQ(nii,nbri,NXQ(nij*2,1,NXQ(nik*3,1,nq)))
                      IF(mqi.GT.0) THEN
                        IF((NXQ(nij,0,nq).GT.1).AND.(nii.EQ.0)) THEN
                          !if vertical branches
                          PHIMQ(nii,nij,nik)=PHIMQ(nii,nij,nik)+
     '                      (YQ(mqi,niqOLDSOLN,1)/DBLE(NXQ(nij,0,nq)))
Cb                          PHIEQ(nii,nij,nik)=PHIEQ(nii,nij,nik)+
Cb     '                      (YQ(mqi,3,3)/DBLE(NXQ(nij,0,nq))) 
                          PHIEQ(nii,nij,nik)=0.0d0
                        ELSE IF((NXQ(nii,0,nq).GT.1).AND.(nij.EQ.0))
     '                    THEN  !if horizontal branches
                          PHIMQ(nii,nij,nik)=PHIMQ(nii,nij,nik)+
     '                      (YQ(mqi,niqOLDSOLN,1)/DBLE(NXQ(nii,0,nq)))
Cb                          PHIEQ(nii,nij,nik)=PHIEQ(nii,nij,nik)+
Cb     '                      (YQ(mqi,3,3)/DBLE(NXQ(nii,0,nq)))    
                          PHIEQ(nii,nij,nik)=0.0d0
                        ELSE !if no branches
                          PHIMQ(nii,nij,nik)=PHIMQ(nii,nij,nik)+
     '                      YQ(mqi,niqOLDSOLN,1)
Cb                          PHIEQ(nii,nij,nik)=PHIEQ(nii,nij,nik)+
Cb     '                      YQ(mqi,3,3)
                          PHIEQ(nii,nij,nik)=0.0d0
                        ENDIF      
                      ENDIF
                    ENDDO !nbri
                  ENDIF
                ENDDO
              ENDDO
            ENDDO

              if(dop) then
C$              call mp_setlock()
                WRITE(OP_STRING,'('' nq: '',I5)') nq
                CALL WRITES(IODI,OP_STRING,ERROR,*100)
C$              call mp_unsetlock()
              endif

C     IF(DOP) THEN
C     WRITE(OP_STRING,'('' Phi(m) around nq'')')
C     CALL WRITES(IODI,OP_STRING,ERROR,*100)
C     DO NIK=-IK,IK
C     DO NIJ=-IJ,IJ
C     WRITE(OP_STRING,'(3F12.6)')
C     '      (PHIMQ(NII,NIJ,NIK),NII=-1,1)
C     CALL WRITES(IODI,OP_STRING,ERROR,*100)
C     ENDDO
C     ENDDO
C     WRITE(OP_STRING,'('' Phi(e) around nq'')')
C     CALL WRITES(IODI,OP_STRING,ERROR,*100)
C     DO NIK=-IK,IK
C     DO NIJ=-IJ,IJ
C     WRITE(OP_STRING,'(3F12.6)')
C     '      (PHIEQ(NII,NIJ,NIK),NII=-1,1)
C     CALL WRITES(IODI,OP_STRING,ERROR,*100)
C     ENDDO
C     ENDDO
C     ENDIF
C     *** Compute PHI,k by first order finite differences about nq.
              DPHIDX(1)=PHIMQ(1,0,0)-PHIMQ(-1,0,0)
     '          +PHIEQ(1,0,0)-PHIEQ(-1,0,0) 
              DPHIDX(2)=PHIMQ(0,1,0)-PHIMQ(0,-1,0)
     '          +PHIEQ(0,1,0)-PHIEQ(0,-1,0) 
              DPHIDX(3)=PHIMQ(0,0,1)-PHIMQ(0,0,-1)
     '          +PHIEQ(0,0,1)-PHIEQ(0,0,-1)
              IF(DOP) THEN
C$              call mp_setlock()
                WRITE(OP_STRING,'('' dphi/dxi = '',3F10.5)')
     '            (dphidx(ni),ni=1,3)
                CALL WRITES(IODI,OP_STRING,ERROR,*100)
C$              call mp_unsetlock()
              ENDIF
C***  Compute PHI,ij = d2PHI/dxi(i)dxi(j) by taking second order finite
C***  differences about nq.
              DO ib=1,NITB
                DO jb=1,NITB
                  IF(ib.EQ.jb) THEN
                    IF(ib.EQ.1) THEN !PHI,11
                      D2PHIDX2(ib,jb)=(PHIMQ(1,0,0)-2.0d0*PHIMQ(0,0,0)+
     '                  PHIMQ(-1,0,0)
     '                  +PHIEQ(1,0,0)-2.0d0*PHIEQ(0,0,0)+PHIEQ(-1,0,0)
     '                  )*4.0d0
                    ELSE IF(ib.EQ.2) THEN !PHI,22
                      D2PHIDX2(ib,jb)=(PHIMQ(0,1,0)-2.0d0*PHIMQ(0,0,0)+
     '                  PHIMQ(0,-1,0)
     '                  +PHIEQ(0,1,0)-2.0d0*PHIEQ(0,0,0)+PHIEQ(0,-1,0)
     '                  )*4.0d0
                    ELSE IF(ib.EQ.3) THEN !PHI,33
                      D2PHIDX2(ib,jb)=(PHIMQ(0,0,1)-2.0d0*PHIMQ(0,0,0)+
     '                  PHIMQ(0,0,-1)
     '                  +PHIEQ(0,0,1)-2.0d0*PHIEQ(0,0,0)+PHIEQ(0,0,-1)
     '                  )*4.0d0
                    ENDIF               
                  ELSE
                    IF(ib+jb.EQ.3) THEN !PHI,12 and PHI,21
                      D2PHIDX2(ib,jb)=PHIMQ(1,1,0)-PHIMQ(1,-1,0)-
     '                  PHIMQ(-1,1,0)+PHIMQ(-1,-1,0)
     '                  +PHIEQ(1,1,0)-
     '                  PHIEQ(1,-1,0)-PHIEQ(-1,1,0)+PHIEQ(-1,-1,0)
                    ELSE IF(ib+jb.EQ.4) THEN !PHI,13 and PHI,31
                      D2PHIDX2(ib,jb)=PHIMQ(1,0,1)-PHIMQ(1,0,-1)-
     '                  PHIMQ(-1,0,1)+PHIMQ(-1,0,-1)
     '                  +PHIEQ(1,0,1)-
     '                  PHIEQ(1,0,-1)-PHIEQ(-1,0,1)+PHIEQ(-1,0,-1)
                    ELSE IF(ib+jb.EQ.5) THEN !PHI,23 and PHI,32
                      D2PHIDX2(ib,jb)=PHIMQ(0,1,1)-PHIMQ(0,-1,1)-
     '                  PHIMQ(0,1,-1)+PHIMQ(0,-1,-1)
     '                  +PHIEQ(0,1,1)-
     '                  PHIEQ(0,-1,1)-PHIEQ(0,1,-1)+PHIEQ(0,-1,-1)
                    ENDIF
                  ENDIF
                ENDDO
              ENDDO
              IF(DOP) THEN
C$              call mp_setlock()
                WRITE(OP_STRING,'('' d2phi/dxi2'')')
                CALL WRITES(IODI,OP_STRING,ERROR,*100)
                DO nj=1,nitb
                  WRITE(OP_STRING,'(3F10.5)') (d2phidx2(ni,nj),
     '              ni=1,nitb)
                  CALL WRITES(IODI,OP_STRING,ERROR,*100)
                ENDDO
C$              call mp_unsetlock()
              ENDIF
          
C***  Compute Iprop = propagated current 
              IPROP=0.0d0
              DO nii=1,NITB
                DO nij=1,NITB           
                  SUM1=0.0d0
                  SUM2=0.0d0
                  DO nik=1,NITB
                    SUM1=SUM1+PROPQ(nik,nii,nij+1,1,nq)*DPHIDX(nik)
                    SUM2=SUM2+PROPQ(nik,nii,1,1,nq)*D2PHIDX2(nij,nik)
                  ENDDO
                  IPROP=IPROP+(SUM1+SUM2)*GUQ(nii,nij,nq)
                ENDDO !nij
              ENDDO !nii
              DO nii=1,NITB
                SUM1=0.0d0
                DO nij=1,NITB
                  SUM1=SUM1+PROPQ(nij,nii,1,1,nq)*DPHIDX(nij)
                ENDDO
                IPROP=IPROP-SUM1*GCHQ(nii,nq)
              ENDDO !nii
              IF(DOP) THEN
C$              call mp_setlock()
                WRITE(OP_STRING,'('' Iprop = '',F15.5,'' uA/mm^3'')')
     '            iprop*1.d3
                CALL WRITES(IODI,OP_STRING,ERROR,*100)
C$              call mp_unsetlock()
              ENDIF

              IAPP=0.0d0
              IF((AQ(maqp1i,nq).GT.ZERO_TOL).AND.(T.GE.AQ(maqp1t0,nq))
     '          .AND.(T.LT.AQ(maqp1t1,nq))) THEN
                IAPP=AQ(maqp1i,nq)
              ELSE IF((AQ(maqp2i,nq).GT.ZERO_TOL).AND.(T.GE.
     '          AQ(maqp2t0,nq)).AND.(T.LT.AQ(maqp2t1,nq))) THEN
                IAPP=AQ(maqp2i,nq)
              ENDIF

              IF(SAC) THEN
                IF(YQ(nq,niqSAC,1).GT.1.d-8) THEN !inward SAC current
                  IAPP=IAPP+YQ(nq,niqSAC,1)
                  IF(DOP) THEN
C$                  call mp_setlock()
                    write(*,'('' Inward SAC current at nq='',I5,'
     '                //''' is'',E12.3)') nq,IAPP
C$                  call mp_unsetlock()
                  ENDIF
                ENDIF
              ENDIF

              IF(DOP) THEN
C$              call mp_setlock()
                WRITE(OP_STRING,'('' Iapp  = '',F15.5,'' uA/mm^3'')')
     '            iapp*1.d3
                CALL WRITES(IODI,OP_STRING,ERROR,*100)
C$              call mp_unsetlock()
              ENDIF
          
C***  Compute transmembrane current.
              AM=CQ(2,nq)
              IM=(IPROP+IAPP)/AM
              IF(DOP) THEN
C$              call mp_setlock()
                WRITE(OP_STRING,'('' Im    = '',F15.5,'' uA/mm^2'')')
     '            IM*1.d3
                CALL WRITES(IODI,OP_STRING,ERROR,*100)
C$              call mp_unsetlock()
              ENDIF

C*** Using Euler Predictor-Corrector method to solve for new phim
              PHIM=YQ(nq,niqOLDSOLN,1)
              IF(ITYP3(nr,nx).EQ.1) THEN !Cubic w/o recovery
                IION=CUBIC_ION(PHIM,CQ(9,nq))
              ELSE IF(ITYP3(nr,nx).EQ.2) THEN !FHN
                IION=FHN_ION(PHIM,YQ(nq,niqRECOV,1),CQ(9,nq))
              ELSE IF(ITYP3(nr,nx).EQ.3) THEN !VCD
                IION=VCD_ION(PHIM,YQ(nq,niqRECOV,1),YQ(nq,niqDV,1))
              ELSE IF(ITYP3(nr,nx).EQ.4) THEN !BR
                !BR_ION(Vm,m,h,j,d,f,x1,Cai,...)
                CALL ASSERT(.FALSE.,'>>You must use march8',ERROR,*100)
C                IION=BR_ION(PHIM,YQ(nq,niqMM,1),YQ(nq,niqH,1),
C     '            YQ(nq,niqJ,1),YQ(nq,niqD,1),YQ(nq,niqF,1),
C     '            YQ(nq,niqX,1),YQ(nq,niqCAI,1),CQ(9,nq))
              ELSE IF(ITYP3(nr,nx).EQ.5) THEN !unused
              ELSE IF(ITYP3(nr,nx).EQ.6) THEN !Luo-Rudy
                IION=LR_ION(PHIM,YQ(nq,niqMM,1),YQ(nq,niqH,1),
     '            YQ(nq,niqJ,1),YQ(nq,niqD,1),YQ(nq,niqF,1),
     '            YQ(nq,niqX,1),YQ(nq,niqCAI,1),CQ(9,nq))
              ELSE IF(ITYP3(nr,nx).EQ.7) THEN !diFrancesco-Noble
C                IION=DN_ION(PHIM,YQ(nq,3,1),YQ(nq,4,1),CQ(9,nq))
              ENDIF
              IF(DOP) THEN
C$              call mp_setlock()
                WRITE(OP_STRING,'('' Iion  = '',F15.5,'' uA/mm^2'')')
     '            IION*1.d3
                CALL WRITES(IODI,OP_STRING,ERROR,*100)
C$              call mp_unsetlock()
              ENDIF
              CM=CQ(1,nq)
C     10^-6 correction value for Cm      
              DPHIM=(IM-IION)/CM*1.0d6
              PHIMP=PHIM+DPHIM*DT   !Predicted value of Phi(m)
              
CC              DVMST(nq)=DPHIM*DT
              
c     IF(DOP) THEN
c     WRITE(OP_STRING,'('' dphim,phimp: '',4F15.6)')
c     '        DPHIM,PHIMP
c     CALL WRITES(IODI,OP_STRING,ERROR,*100)
c     ENDIF
c     IF(ITYP3(nr,nx).EQ.1) THEN  !Cubic w/o recovery
c     IION=CUBIC_ION(PHIMP,CQ(9,nq))
c     ELSE IF(ITYP3(nr,nx).EQ.2) THEN  !FHN
c     IION=FHN_ION(PHIMP,YQ(nq,3,1),CQ(9,nq))
c     ELSE IF(ITYP3(nr,nx).EQ.3) THEN  !VCD
c     IION=VCD_ION(PHIMP,YQ(nq,3,1),YQ(nq,4,1))
c     ELSE IF(ITYP3(nr,nx).EQ.4) THEN  !EBJ
c     IION=EBJ_ION(PHIMP,YQ(nq,3,1),YQ(nq,4,1),CQ(9,nq))
c     ENDIF
c     IF(DOP) THEN
c     WRITE(OP_STRING,'('' Iion = '',F15.5)') IION
c     CALL WRITES(IODI,OP_STRING,ERROR,*100)
c     ENDIF
C**** Update the transmembrane potential
c     YQ(nq,1,1)=PHIM+DT*0.5d0*(DPHIM+(IM-IION)/CM) !New phi(m)
              YQ(nq,niqV,1)=PHIMP
              V=PHIMP
              IF(ITYP3(nr,nx).EQ.3.AND.KTYP33.EQ.2) THEN !VCDC mods
                YQ(nq,niqDV,1)=DPHIM
              ENDIF
C**** Store activation time for optimisation      
              IF(DABS(DPHIM).GT.AQ(dpot,nq)) THEN
                AQ(dpot,nq)=DABS(DPHIM)!Largest abs. change in potential
                AQ(atime,nq)=T         !Activation time
              ENDIF

C**** Update the recovery variable if necessary (1 step)
              IF(ITYP3(nr,nx).EQ.2) THEN !FHN
                RECOV=YQ(nq,niqRECOV,1)
 !CQ( 9,nq) is the rest potential
 !CQ(10,nq) is the plateau potential      
 !CQ(11,nq) is the threshold potential      
 !NORM_PHI  is a 0-1 normalised potential
                DIFF = CQ(10,nq)-CQ(9,nq)
                NORM_PHI=(V-CQ(9,nq))/DIFF
                IF(KTYP33.LT.3) THEN !Standard or Roger's FHN
 !CQ(14,nq) is the recovery rate constant
 !CQ(15,nq) is the recovery decay constant
                  D_RECOV=CQ(14,nq)*(NORM_PHI-CQ(15,nq)*RECOV)
                ELSE !Panfilov FHN
 !CQ(12,nq) is rate constant (same because of identical parabola)
 !CQ(14,nq) is epsilon0
 !CQ(15,nq) is mu1
 !CQ(16,nq) is mu2
 !CQ(17,nq) is tau (time const for force)
                  EPS=CQ(14,nq)+CQ(15,nq)*RECOV/(NORM_PHI+CQ(16,nq))
                  RATE=CQ(12,nq)/DIFF/DIFF
                  VSTAR=V+CQ(9,nq)-CQ(10,nq)-CQ(11,nq)
                  D_RECOV=EPS*(-RECOV-RATE*VSTAR*(V-CQ(9,nq)))
                ENDIF
                YQ(nq,niqRECOV,1)=RECOV+DT*D_RECOV
!Compute calcium level (simple equation from Panfilov)
                IF(DABS(CQ(17,nq)).GT.1d-12) THEN
                  IF(NORM_PHI.LT.0.d0) NORM_PHI=0.d0
                  YQ(nq,niqCALCIUM,1)=YQ(nq,niqCALCIUM,1)+
     '              DT*(NORM_PHI-YQ(nq,niqCALCIUM,1))/CQ(17,nq)
                ELSE
                  YQ(nq,niqCALCIUM,1)=0.d0
                ENDIF

              ELSE IF(ITYP3(nr,nx).EQ.3) THEN !vanCapelle-Durrer
                RECOV=YQ(nq,niqRECOV,1)
                IF(KTYP33.EQ.1) THEN !Original VCD
                  T_CONST=CQ(10,nq)
                ELSE IF(KTYP33.EQ.2) THEN !VCDC mods
C     From modifications provided by Alan Garfinkel
                  IF(YQ(nq,niqDRECOV,1).ge.0.0d0) THEN
                    IF(KTYP34.EQ.1) THEN !"Normal" APD (209ms)
                      T_CONST=0.5d0
                    ELSE IF(KTYP34.EQ.2) THEN !"Ischemic" APD (112ms)
                      T_CONST=0.33d0
                    ENDIF
                  ELSE IF(recov.gt.0.85d0) THEN
                    IF(KTYP34.EQ.1) THEN !"Normal" APD (209ms)
                      T_CONST=0.1d0
                    ELSE IF(KTYP34.EQ.2) THEN !"Ischemic" APD (112ms)
                      T_CONST=0.066d0
                    ENDIF
                  ELSE
                    IF(KTYP34.EQ.1) THEN !"Normal" APD (209ms)
                      T_CONST=3.0d0
                    ELSE IF(KTYP34.EQ.2) THEN !"Ischemic" APD (112ms)
                      T_CONST=3.31d0
                    ENDIF
                  ENDIF
C     Scale Factor for time constant
                  T_CONST=T_CONST*CQ(10,nq)
                ENDIF
                IF(V.LT.-80.0d0) THEN
                  D_RECOV=-RECOV/T_CONST
                ELSE IF(V.GT.-60.0d0) THEN
                  D_RECOV=(1.0d0-RECOV)/T_CONST
                ELSE
                  D_RECOV=((V+80.0d0)/20.0d0-RECOV)/T_CONST
                ENDIF
                YQ(nq,niqRECOV,1)=RECOV+DT*D_RECOV
                IF(KTYP33.EQ.2) THEN !VCDC mods
                  YQ(nq,niqDRECOV,1)=D_RECOV
                ENDIF
!Compute calcium level (simple equation from Panfilov)
                IF(DABS(CQ(12,nq)).GT.1d-12) THEN
                    NORM_PHI=(V-CQ(9,nq))/100.d0
                  IF(NORM_PHI.LT.0.d0) NORM_PHI=0.d0
                  YQ(nq,niqCALCIUM,1)=YQ(nq,niqCALCIUM,1)+
     '              DT*(NORM_PHI-YQ(nq,niqCALCIUM,1))/CQ(12,nq)
                ELSE
                  YQ(nq,niqCALCIUM,1)=0.d0
                ENDIF

              ELSE IF(ITYP3(nr,nx).EQ.4) THEN !Beeler-Reuter
                  !BR_RATES(Vm,m,h,j,d,f,x1,Cai,...)
C                CALL BR_RATES(V,YQ(nq,niqMM,1),YQ(nq,niqH,1),
C     '            YQ(nq,niqJ,1),YQ(nq,niqD,1),YQ(nq,niqF,1),
C     '            YQ(nq,niqX,1),YQ(nq,niqCAI,1),CQ(9,nq))

              ELSE IF(ITYP3(nr,nx).EQ.5) THEN !Unused

              ELSE IF(ITYP3(nr,nx).EQ.6) THEN !Luo-Rudy
                CALL LR_RATES(V,YQ(nq,niqMM,1),YQ(nq,niqH,1),
     '            YQ(nq,niqJ,1),YQ(nq,niqD,1),YQ(nq,niqF,1),
     '            YQ(nq,niqX,1),YQ(nq,niqCAI,1),CQ(9,nq))

              ELSE IF(ITYP3(nr,nx).EQ.7) THEN !diFrancesco-Noble

              ENDIF !ityp3

C     Removal from DTAR solution - Needs more thought!
              IF(KTYP36.EQ.1) THEN
C     For cubic - switch off if Phi is within 0.01 mV of plateau
                IF(ITYP3(nr,nx).EQ.1) THEN
                  IF(CQ(10,nq)-YQ(nq,niqV,1).LT.0.01) THEN
                    NWQ(4,nq)=-1
C MLB trying to increase speed
C                    NWQ(5,nq)=-1
                    if(dop) then
C$                    call mp_setlock()
                      WRITE(OP_STRING,'(''  Removing nq : '',I6)') nq
                      CALL WRITES(IODI,OP_STRING,ERROR,*100)
C$                    call mp_unsetlock()
                    endif
C MLB 12/5/97 
                  ELSE IF(YQ(nq,niqV,1).GT.CQ(10,nq)) THEN
                    NWQ(4,nq)=-1
C MLB trying to increase speed
C                    NWQ(5,nq)=-1
                    if(dop) then
C$                    call mp_setlock()
                      WRITE(OP_STRING,'(''  Removing nq : '',I6)') nq
                      CALL WRITES(IODI,OP_STRING,ERROR,*100)
C$                    call mp_unsetlock()
                    endif
                  ENDIF
C     Switch off if dPhim is less than 10^-5 for 100 time steps
                ELSE IF(ABS(DPHIM).lt.1.d-5) THEN
                  NWQ(4,nq)=NWQ(4,nq)+1
                  IF(NWQ(4,nq).EQ.102) THEN
                    NWQ(4,nq)=-1
C MLB trying to increase speed
C                    NWQ(5,nq)=-1
                    if(dop) then
C$                    call mp_setlock()
                      WRITE(OP_STRING,'(''  Removing nq : '',I6)') nq
                      CALL WRITES(IODI,OP_STRING,ERROR,*100)
C$                    call mp_unsetlock()
                    endif
                  ENDIF
                ENDIF
              ENDIF
          
              IF(DOP) THEN
C$              call mp_setlock()
                WRITE(OP_STRING,
     '            '('' Old phim: '',E14.5,'' New phim: '',E14.5)')
     '            PHIM,YQ(nq,niqV,1)                          
                CALL WRITES(IODI,OP_STRING,ERROR,*100)
                IF(ITYP3(nr,nx).EQ.2.OR.ITYP3(nr,nx).EQ.3) THEN
                  !FHN or VCD
                  WRITE(OP_STRING,'('' Recovery: '',2E14.5)')
     '              RECOV,YQ(nq,niqRECOV,1)
                  CALL WRITES(IODI,OP_STRING,ERROR,*100)
                ELSE IF(ITYP3(nr,nx).EQ.4) THEN !BR
                  WRITE(OP_STRING,'('' m,h,j,d,f,x1,Ca(i)'
     '              //': '',7E14.5)')
     '              (YQ(nq,niq,1),niq=5,9),(YQ(nq,niq,1),niq=3,4)
                  CALL WRITES(IODI,OP_STRING,ERROR,*100)
                ENDIF
                WRITE(OP_STRING,'('' ============='')')
                CALL WRITES(IODI,OP_STRING,ERROR,*100)
C$              call mp_unsetlock()
              ENDIF
            ENDIF !Internal grid point
          ENDIF !In the current region
          GOTO 200
C          IF(ERROR_FLAG) THEN
C         This statement is designed to be skipped if no error 
C         occur. However if a error occurs within a subroutine 
C         the alternate return points to line 100 to set the flag
 100        CONTINUE
C$          call mp_setlock()
            ERROR_FLAG=.TRUE.
            WRITE(OP_STRING,'(/'' >>ERROR: An error occurred - '
     '        //'results may be unreliable!'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*101)
 101        CONTINUE
C$          call mp_unsetlock()
C          ENDIF
 200      CONTINUE
        ENDDO                     !Active points
        CALL ASSERT(.NOT.ERROR_FLAG,'>>An error occurred during '
     '    //'grid point calculations',ERROR,*9999)

      ENDDO !nrr - region loop

C     Update bdy grid points
      IF(DOP) THEN
C$      call mp_setlock()
        WRITE(OP_STRING,'('' Updating boundary conditions'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C$      call mp_unsetlock()
      ENDIF
C     Loop over all active grid points & solve
        ERROR_FLAG=.FALSE.
C$DOACROSS local(nq,PHIM1,PHIM2,PHIMP,DPHIM),
C$&        share(ERROR_FLAG,NQT,NWQ,YQ,T),
C$&        reduction(NSOL_EXT)
      DO nq=1,NQT
C     Generic no-flux bdy cond.n applied at external bdy point
C     only if internal points are active
        IF(NWQ(1,nq).GT.0) THEN !On external bdy
          IF(DOP) THEN
C$          call mp_setlock()
            WRITE(OP_STRING,'('' Boundary GP'',I4)') nq
            CALL WRITES(IODI,OP_STRING,ERROR,*300)
C$          call mp_unsetlock()
          ENDIF
          IF(NWQ(4,NWQ(1,nq)).GT.0.OR.NWQ(4,NWQ(2,nq)).GT.0) THEN
            NSOL_EXT=NSOL_EXT+1
            PHIM1=YQ(NWQ(1,nq),niqV,1) !Pick up values for g.p. as
            PHIM2=YQ(NWQ(2,nq),niqV,1) !defined in DEGRID
            PHIMP=(4.0d0*PHIM1-PHIM2)/3.0d0
            DPHIM=PHIMP-YQ(nq,niqV,1)
CC            DVMST(nq)=DPHIM
            YQ(nq,niqV,1)=PHIMP
C**** Store activation time
            IF(DABS(DPHIM).GT.AQ(dpot,nq)) THEN
              AQ(dpot,nq)=DABS(DPHIM) !Largest abs change in potential
              AQ(atime,nq)=T !Activation time
            ENDIF

            IF(DOP) THEN
C$            call mp_setlock()
              WRITE(OP_STRING,
     '          '('' nq: '', I5,'' phim1,phim2,phim,dphim: '',4F12.5)')
     '          nq,phim1,phim2,phimp,dphim
              CALL WRITES(IODI,OP_STRING,ERROR,*300)
C$            call mp_unsetlock()
            ENDIF
          ELSE IF(NWQ(4,NWQ(1,nq)).LT.0.OR.NWQ(4,NWQ(2,nq)).LT.0) THEN
C           Both internal points are removed, therefore remove boundary
C           point from active solution (if DTAR)
            IF(KTYP36.EQ.1) THEN
              NWQ(4,nq)=-1
              if(dop) then
C$              call mp_setlock()
                WRITE(OP_STRING,'(''    Removing nq : '',I6)') nq
                CALL WRITES(IODI,OP_STRING,ERROR,*300)
C$              call mp_unsetlock()
              endif
            ENDIF
          ENDIF
        ENDIF
        GOTO 400
C        IF(ERROR_FLAG) THEN
C         This statement is designed to be skipped if no error 
C         occur. However if a error occurs within a subroutine 
C         the alternate return points to line 300 to set the flag
 300      CONTINUE
C$        call mp_setlock()
          ERROR_FLAG=.TRUE.
          WRITE(OP_STRING,'(/'' >>ERROR: An error occurred - '
     '      //'results may be unreliable!'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*301)
 301      CONTINUE
C$        call mp_unsetlock()
C        ENDIF
 400    CONTINUE
      ENDDO                     !Active points

      DO nnq=1,CPLST(0,1)
        nq=CPLST(nnq,1)
C        IF((NWQ(6,nq).GT.0).AND.(nq.GT.NWQ(6,nq))) THEN
          DVMST=YQ(nq,niqV,1)-YQ(NWQ(6,nq),niqV,1)
          IF((DVMST.LT.ZERO_TOL).AND.
     '      (YQ(NWQ(6,nq),niqV,1).LE.CQ(10,nq))) THEN
            YQ(nq,niqV,1)=YQ(NWQ(6,nq),niqV,1)
          ENDIF
          IF((DVMST.GT.ZERO_TOL).AND.
     '      (YQ(nq,niqV,1).LE.CQ(10,nq))) THEN
            YQ(NWQ(6,nq),niqV,1)=YQ(nq,niqV,1)
          ENDIF
C        ENDIF
      ENDDO

      CALL ASSERT(.NOT.ERROR_FLAG,'>>An error occurred during '
     '  //'boundary point calculations',ERROR,*9999)

c     WRITE(OP_STRING,'('' ----- Iteration time:'',F8.2,'
c     '  //''' s.  Active (int:bdy/DTAR): '',I6,'':'',I6,''/'',I6)')
c     '  cpu,nsol_int,nsol_ext,nwq(5,0)
c     CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

C     Total number of active points
      NWQ(4,0)=NSOL_INT+NSOL_EXT

C     Determine the minimum index of NWQ(,) containing an active 
C     data point for following iteration
      nnq=nnq_min
      DO WHILE ((NWQ(4,NWQ(5,nnq)).LE.0).AND.(nnq.LE.NQT))
        nnq=nnq+1
      ENDDO
      nnq_min=nnq
      
      CALL EXITS('BFRONT')
      RETURN
 9999 CALL ERRORS('BFRONT',ERROR)
      CALL EXITS('BFRONT')
      RETURN 1
      END


      SUBROUTINE CALC_ADAMS_DTAR(maqp1t0,maqp1t1,maqp1i,maqp2t0,
     '  maqp2t1,maqp2i,niq_old,niqV,nnq_min,NQXI,NRLIST,NSOL,NWQ,
     '  nx,NXQ,ADAMS_WORK,AQ,CQ,T,YQ,FIRST,ERROR,*)

C#### Subroutine: CALC_ADAMS_DTAR
C###  Description:
C###    CALC_ADAMS_DTAR calculates the list of "active" grid points
C###    that need to be integrated. Points need to be integrated if
C###    1) A stimulus current is applied to that point at the current
C###    time; 2) A wavefront is approaching such that the membrane
C###    potential increases past a threshold amount; 3) The time is
C###    past the next integration time for the adaptive Adams-Moulton
C###    integrator i.e. the current time is greater than the last
C###    integration output time + the last adaptive step size.
C###  See-Also: CALC_DTAR

      IMPLICIT NONE

      INCLUDE 'cmiss$reference:adam00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:grid00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:nqloc00.inc'
      INCLUDE 'cmiss$reference:tol00.cmn'

!     Parameter list
      INTEGER maqp1t0,maqp1t1,maqp1i,maqp2t0,maqp2t1,maqp2i,niq_old,
     '  niqV,nnq_min,NQXI(0:NIM,NQSCM),NRLIST(0:NRM),NSOL,NWQ(6,0:NQM),
     '  nx,NXQ(-NIM:NIM,0:4,0:NQM)
      REAL*8 ADAMS_WORK(ADAMS_LWORK,NQM),AQ(NMAQM,NQM),CQ(NMM,NQM),T,
     '  YQ(NYQM,NIQM,NAM)
      LOGICAL FIRST
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER IJ,IK,mq,mqi,nbri,nii,nij,nik,NITB,niqSAC,nnq,nq,nr,nrr
      REAL*8 IAPP,VTHRESH
      LOGICAL ERROR_FLAG,SAC

      CALL ENTERS('CALC_ADAMS_DTAR',*9999)

      NITB=NQXI(0,1)
      IK=MAX(0,NITB-2) !zero for 1,2-D, one for 3-D
      IJ=MIN(NITB-1,1) !zero for 1-D, one for 2,3-D
      NWQ(5,0)=0
      DO nrr=1,NRLIST(0)
        nr=NRLIST(nrr)
        ERROR_FLAG=.FALSE.          
C$OMP PARALLEL DO
C$&   PRIVATE(nik,nij,nii,mq,mqi,nbri,nq,IAPP,SAC,niqSAC)
C$&   SHARE(ADAMS_WORK,AQ,CQ,maqp1i,maqp1t0,maqp1t1,maqp2i,maqp2t0,
C$&     maqp2t1,niq_old,niqV,nr,nx,NWQ,T,YQ,ERROR_FLAG,IODI,IOER)
        DO nq=NQR(1,nr),NQR(2,nr)
          IF(.NOT.ERROR_FLAG) THEN
            IF(NWQ(1,nq).EQ.0) THEN !internal
              IAPP=0.0d0
              IF((AQ(maqp1i,nq).GT.ZERO_TOL).AND.(T.GE.AQ(maqp1t0,nq))
     '          .AND.(T.LT.AQ(maqp1t1,nq))) THEN
                IAPP=AQ(maqp1i,nq)
              ELSE IF((AQ(maqp2i,nq).GT.ZERO_TOL).AND.(T.GE.
     '            AQ(maqp2t0,nq)).AND.(T.LT.AQ(maqp2t1,nq))) THEN
                IAPP=AQ(maqp2i,nq)
              ENDIF
              SAC=.FALSE.
              IF(ITYP3(nr,nx).EQ.2) THEN !FHN
                IF(DABS(CQ(18,nq)).GT.ZERO_TOL) SAC=.TRUE.
              ELSEIF(ITYP3(nr,nx).EQ.3) THEN !VCD
                IF(DABS(CQ(13,nq)).GT.ZERO_TOL) SAC=.TRUE.
              ENDIF              
              IF(SAC) THEN
                CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqSAC,NIQ_SAC,
     '            ERROR,*100)
                IAPP=IAPP-YQ(nq,niqSAC,1)
              ENDIF
              
C             Add to list of points to integrate if
              IF(FIRST) THEN
C               First time around
                NWQ(4,nq)=1
C$OMP CRITICAL(CALC_ADAMS_DTAR_1)
                NWQ(5,0)=NWQ(5,0)+1
C$OMP END CRITICAL(CALC_ADAMS_DTAR_1)
                NWQ(5,NWQ(5,0))=nq
              ELSE IF(IAPP.GT.ZERO_TOL) THEN
C               There is applied current
                NWQ(4,nq)=1
C$OMP CRITICAL(CALC_ADAMS_DTAR_2)
                NWQ(5,0)=NWQ(5,0)+1
C$OMP END CRITICAL(CALC_ADAMS_DTAR_2)
                NWQ(5,NWQ(5,0))=nq
              ELSE IF(T.GE.ADAMS_WORK(4,nq)+ADAMS_WORK(2,nq)) THEN
C               or is time to integrate
                NWQ(4,nq)=1
C$OMP CRITICAL(CALC_ADAMS_DTAR_3)
                NWQ(5,0)=NWQ(5,0)+1
C$OMP END CRITICAL(CALC_ADAMS_DTAR_3)
                NWQ(5,NWQ(5,0))=nq
              ELSE IF(DABS(YQ(nq,niqV,1)-YQ(nq,niq_old,1)).GT.
     '            0.5d0) THEN
C                 or has changed by 1 mV since last update. In this case
C                 add all neighbouring points as well.
                NWQ(4,nq)=1
C$OMP CRITICAL(CALC_ADAMS_DTAR_4)
                NWQ(5,0)=NWQ(5,0)+1
                NWQ(5,NWQ(5,0))=nq
C$OMP END CRITICAL(CALC_ADAMS_DTAR_4)
                DO nik=-IK,IK
                  DO nij=-IJ,IJ
                    DO nii=-1,1
                      IF(NXQ(nii,0,nq).LE.1) THEN
                        mq=NXQ(nii,1,NXQ(nij*2,1,NXQ(nik*3,1,nq)))
                        IF(mq.GT.0.AND.NWQ(4,mq).EQ.0) THEN 
C                         Not already active
                          NWQ(4,mq)=1
C$OMP CRITICAL(CALC_ADAMS_DTAR_5)
                          NWQ(5,0)=NWQ(5,0)+1
                          NWQ(5,NWQ(5,0))=mq
C$OMP END CRITICAL(CALC_ADAMS_DTAR_5)
                        ENDIF
                      ELSE
                        DO nbri=1,NXQ(nii,0,nq) !branch in xi1 from nq
                          mqi=NXQ(nii,nbri,NXQ(nij*2,1,
     '                      NXQ(nik*3,1,nq)))
                          IF(mqi.GT.0.AND.NWQ(4,mqi).EQ.0) THEN
C                             if point exists & not already active
                            NWQ(4,mqi)=1
C$OMP CRITICAL(CALC_ADAMS_DTAR_6)
                            NWQ(5,0)=NWQ(5,0)+1
                            NWQ(5,NWQ(5,0))=mqi
C$OMP END CRITICAL(CALC_ADAMS_DTAR_6)
                            IF(NWQ(6,mqi).GT.0) THEN
C$OMP CRITICAL(CALC_ADAMS_DTAR_7)
                              NWQ(4,NWQ(6,mqi))=1
                              NWQ(5,0)=NWQ(5,0)+1
                              NWQ(5,NWQ(5,0))=NWQ(6,mqi)
C$OMP END CRITICAL(CALC_ADAMS_DTAR_7)
                            ENDIF !coupled nodes
                          ENDIF
                        ENDDO !nbri
                      ENDIF
                    ENDDO
                  ENDDO
                ENDDO
              ELSE
C               inactivate point
                NWQ(4,nq)=0
              ENDIF
            ENDIF
            
            GOTO 102
 100        CONTINUE
C$OMP CRITICAL(CALC_ADAMS_DTAR_8)
            ERROR_FLAG=.TRUE.
            WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
            CALL WRITES(IOER,OP_STRING,ERROR,*101)
            WRITE(OP_STRING,'(/'' >>An error occurred - '
     '          //'results may be unreliable!'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*101)
 101        CONTINUE
C$OMP END CRITICAL(CALC_ADAMS_DTAR_8)
 102        CONTINUE
          ENDIF !.NOT.ERROR_FLAG             
        ENDDO !nq
C$OMP END PARALLEL DO
      ENDDO !nr
      IF(FIRST) FIRST=.FALSE.

      CALL EXITS('CALC_ADAMS_DTAR')
      RETURN
 9999 CALL ERRORS('CALC_ADAMS_DTAR',ERROR)
      CALL EXITS('CALC_ADAMS_DTAR')
      RETURN 1
      END


      SUBROUTINE CALC_ADAMS_EXT_RHS(niqV,niqBNDRY,NQGP,NQGP_PIVOT,
     '  NQGW,NRLIST,NWQ,nx_ext,nx_trans,DT,PROPQ,RHS,YQ,CQ,FIXQ,ERROR,*)

C#### Subroutine: CALC_ADAMS_EXT_RHS
C###  Description:
C###    CALC_ADAMS_EXT_RHS is used to calcualte the right hand side
C###    vectors at each time step for the extracellular domain in 
C###    the solution of grid activation problems.

      IMPLICIT NONE

      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:grid00.cmn'

!     Parameter list
      INTEGER niqV,niqBNDRY,NQGP(0:19,NQM),NQGP_PIVOT(19,NQM),
     '  NRLIST(0:NRM),NWQ(6,0:NQM),nx_ext,nx_trans
      REAL*8 DT,NQGW(19,NQM),PROPQ(3,3,4,2,NQM),RHS(NQM),
     '  YQ(NYQM,NIQM,NAM,NXM),CQ(NMM,NQM)
      LOGICAL FIXQ(NYQM,NIYFIXM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER ERRORCODE,i,nq,nr,nrr
      REAL*8 CMAMDT

      CALL ENTERS('CALC_ADAMS_EXT_RHS',*9999)

      ERRORCODE=0
      DO nrr=1,NRLIST(0)
        nr=NRLIST(nrr)
C$OMP PARALLEL DO
C$&   PRIVATE(i,nq)
C$&   SHARED(CQ,ERRORCODE,FIXQ,niqBNDRY,niqV,NQGP,NQGP_PIVOT,NQGW,NWQ,
C$&     nx_ext,nx_trans,PROPQ,RHS,YQ)
        DO nq=NQR(1,nr),NQR(2,nr)
          IF(NWQ(1,nq).EQ.0) THEN !internal
            RHS(nq)=0.0d0
            DO i=1,NQGP(0,nq)
              IF(NQGP(i,nq).GT.0) RHS(nq)=RHS(nq)-
     '          (NQGW(NQGP_PIVOT(i,nq),nq)*
     '          YQ(NQGP(i,nq),niqV,1,nx_trans))
            ENDDO !i
          ELSE !boundary grid point
            IF(FIXQ(nq,1)) THEN
              RHS(nq)=YQ(nq,niqBNDRY,1,nx_ext)
            ELSEIF(FIXQ(nq,2)) THEN
              RHS(nq)=YQ(nq,niqBNDRY,1,nx_ext)
            ELSEIF(FIXQ(nq,3)) THEN
              RHS(nq)=PROPQ(1,1,1,1,nq)/(PROPQ(1,1,1,1,nq)+
     '          PROPQ(1,1,1,2,nq))
              RHS(nq)=RHS(nq)*(YQ(nq,niqV,1,nx_trans)-CQ(9,nq))
            ELSE
              ERRORCODE=nq
            ENDIF
          ENDIF
        ENDDO !nq
C$OMP END PARALLEL DO
      ENDDO !nr

      IF(ERRORCODE.NE.0) THEN
C$OMP CRITICAL(CALC_ADAMS_EXT_RHS_1)
        WRITE(OP_STRING,'('' >>Error: No boundary condition applied at'
     '    //' point '',I8)') nq
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
C$OMP END CRITICAL(CALC_ADAMS_EXT_RHS_1)
      ENDIF

      IF(DOP) THEN
C$OMP CRITICAL(CALC_ADAMS_EXT_RHS_2)
        WRITE(OP_STRING,'('' EXT_RHS:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nq=1,NQT
          WRITE(OP_STRING,'('' RHS('',I8'')='',F12.6'' uA/mm^3'')') 
     '      nq,RHS(nq)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !nq
C$OMP END CRITICAL(CALC_ADAMS_EXT_RHS_2)
      ENDIF

      CALL EXITS('CALC_ADAMS_EXT_RHS')
      RETURN
 9999 CALL ERRORS('CALC_ADAMS_EXT_RHS',ERROR)
      CALL EXITS('CALC_ADAMS_EXT_RHS')
      RETURN 1
      END


      SUBROUTINE CALC_ADAMS_GRID_COEF(nq,NWQ,NXQ,nx_ext,nx_trans,
     '  NITB,COEFFSEXT,GCHQ,GUQ,NQGW,PROPQ,BIDOMAIN,FIXQ,IMPLICIT,
     '  ERROR,*)

C#### Subroutine: CALC_ADAMS_GRID_COEF
C###  Description:
C###    CALC_ADAMS_GRID_COEF finds the coefficients of the diffusion
C###    operator for ionic activation problems i.e. it calculates 
C###    NQGW such that \diverg{\sigma_i\grad{xxx}}=\sum{i=1}{n}
C###    NQGW(n).xxx(NQGP(n)). If BIDOMAIN is .TRUE. the routine 
C###    also calculates the coefficients of the extracellular matrix
C###    i.e. for the operator \diverg{(\sigma_i+\sigma_e)\grad{\phi_e}}

      IMPLICIT NONE

      INCLUDE 'cmiss$reference:b12.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'

!     Parameter list
      INTEGER NITB,nq,NWQ(6),NXQ(-NIM:NIM,0:4),nx_ext,nx_trans
      REAL*8 COEFFSEXT(19),GCHQ(3),GUQ(3,3),NQGW(19),PROPQ(3,3,4,2)
      CHARACTER ERROR*(*)
      LOGICAL BIDOMAIN,FIXQ(NYQM,NIYFIXM,NXM),IMPLICIT
!     Local variables
      INTEGER i,n,ni,nii,nij
      REAL*8 SUM(15)

      CALL ENTERS('CALC_ADAMS_GRID_COEF',*9999)
      
      n=0
      IF(BIDOMAIN) THEN
        DO i=1,19
          COEFFSEXT(i)=0.0d0
        ENDDO
      ENDIF

      IF(NWQ(1).EQ.0) THEN !internal gp
        IF(NITB.EQ.1) THEN
          DO i=1,NXQ(-1,0)
            n=n+1
            NQGW(n)=-GUQ(1,1)*PROPQ(1,1,2,1)+
     '        GCHQ(1)*PROPQ(1,1,1,1)+
     '        4.0d0*GUQ(1,1)*PROPQ(1,1,1,1)
            NQGW(n)=NQGW(n)/DBLE(NXQ(-1,0))
            IF(BIDOMAIN) THEN
              COEFFSEXT(n)=-GUQ(1,1)*(PROPQ(1,1,2,1)+
     '          PROPQ(1,1,2,2))+
     '          GCHQ(1)*(PROPQ(1,1,1,1)+PROPQ(1,1,1,2))+
     '          4.0d0*GUQ(1,1)*(PROPQ(1,1,1,2)+
     '          PROPQ(1,1,1,1))
              COEFFSEXT(n)=COEFFSEXT(n)/DBLE(NXQ(-1,0))
            ENDIF
          ENDDO

          n=n+1
          NQGW(n)=-8.0d0*GUQ(1,1)*PROPQ(1,1,1,1)
          IF(BIDOMAIN) THEN
            COEFFSEXT(n)=-8.0d0*GUQ(1,1)*(PROPQ(1,1,1,2)+
     '        PROPQ(1,1,1,1))
          ENDIF

          DO i=1,NXQ(1,0)
            n=n+1
            NQGW(n)=-GCHQ(1)*PROPQ(1,1,1,1)+
     '        GUQ(1,1)*PROPQ(1,1,2,1)+
     '        4.0d0*GUQ(1,1)*PROPQ(1,1,1,1)
            NQGW(n)=NQGW(n)/DBLE(NXQ(1,0))
            IF(BIDOMAIN) THEN
              COEFFSEXT(n)=-GCHQ(1)*(PROPQ(1,1,1,1)+
     '          PROPQ(1,1,1,2))+
     '          GUQ(1,1)*(PROPQ(1,1,2,1)+PROPQ(1,1,2,2))+
     '          4.0d0*GUQ(1,1)*(PROPQ(1,1,1,2)+
     '          PROPQ(1,1,1,1))
              COEFFSEXT(n)=COEFFSEXT(n)/DBLE(NXQ(1,0))
            ENDIF
          ENDDO
        ELSE IF(NITB.EQ.2) THEN
          DO ni=1,8
            SUM(ni)=0.0d0
          ENDDO !ni
          DO nii=1,NITB
C           C(1j)*G1j(upper)
            SUM(1)=SUM(1)+GUQ(1,nii)*PROPQ(nii,1,1,1)
C           C(2j)*G2j(upper)
            SUM(2)=SUM(2)+GUQ(2,nii)*PROPQ(nii,2,1,1)
C           C(1l)*[Christoffel*Gij(upper)](l)
            SUM(3)=SUM(3)+GCHQ(nii)*PROPQ(nii,1,1,1)
C           C(2l)*[Christoffel*Gij(upper)](l)
            SUM(4)=SUM(4)+GCHQ(nii)*PROPQ(nii,2,1,1)
C           C(1j)*G2j(upper)
            SUM(5)=SUM(5)+GUQ(2,nii)*PROPQ(nii,1,1,1)
C           C(2j)*G1j(upper)
            SUM(6)=SUM(6)+GUQ(1,nii)*PROPQ(nii,2,1,1) 
            DO nij=1,NITB
C             C(1k),j*Gij(upper)
              SUM(7)=SUM(7)+GUQ(nii,nij)*PROPQ(nii,1,1+nij,1)
C             C(2k),j*Gij(upper)
              SUM(8)=SUM(8)+GUQ(nii,nij)*PROPQ(nii,2,1+nij,1)
            ENDDO !nij
          ENDDO !nii
          NQGW(1)=-SUM(5)-SUM(6)
          NQGW(2)=-SUM(8)+SUM(4)+4.0d0*SUM(2)
          NQGW(3)=SUM(5)+SUM(6)
          NQGW(4)=-SUM(7)+SUM(3)+4.0d0*SUM(1)
          NQGW(5)=-8.0d0*(SUM(1)+SUM(2))
          NQGW(6)=SUM(7)-SUM(3)+4.0d0*SUM(1)
          NQGW(7)=SUM(5)+SUM(6)
          NQGW(8)=SUM(8)-SUM(4)+4.0d0*SUM(2)
          NQGW(9)=-SUM(5)-SUM(6)
          IF(BIDOMAIN) THEN
            DO ni=1,8
              SUM(ni)=0.0d0
            ENDDO !ni
            DO nii=1,NITB
C             (C(1j)(int)+C(1j)(ext))*G1j(upper)
              SUM(1)=SUM(1)+GUQ(1,nii)*
     '          (PROPQ(nii,1,1,1)+PROPQ(nii,1,1,2))
C             (C(2j)(int)+C(2j)(ext))*G2j(upper)
              SUM(2)=SUM(2)+GUQ(2,nii)*
     '          (PROPQ(nii,2,1,1)+PROPQ(nii,2,1,2))
C             (C(1l)(int)+C(1l)(ext))*[Christoffel*Gij(upper)](l)
              SUM(3)=SUM(3)+GCHQ(nii)*
     '          (PROPQ(nii,1,1,1)+PROPQ(nii,1,1,2))
C             (C(2l)(int)+C(2l)(ext))*[Christoffel*Gij(upper)](l)
              SUM(4)=SUM(4)+GCHQ(nii)*
     '          (PROPQ(nii,2,1,1)+PROPQ(nii,2,1,2))
C             (C(1j)(int)+C(1j)(ext))*G2j(upper)
              SUM(5)=SUM(5)+GUQ(2,nii)*
     '          (PROPQ(nii,1,1,1)+PROPQ(nii,1,1,2))
C             (C(2j)(int)+C(2j)(ext))*G1j (upper)
              SUM(6)=SUM(6)+GUQ(1,nii)*
     '          (PROPQ(nii,2,1,1)+PROPQ(nii,2,1,2)) 
              DO nij=1,NITB
C               (C(1k),j(int)+C(1k),j(ext))*Gij(upper)
                SUM(7)=SUM(7)+GUQ(nii,nij)*
     '            (PROPQ(nii,1,1+nij,1)+PROPQ(nii,1,1+nij,2))
C               (C(2k),j(int)+C(2k),j(ext))*Gij(upper)
                SUM(8)=SUM(8)+GUQ(nii,nij)*
     '            (PROPQ(nii,2,1+nij,1)+PROPQ(nii,2,1+nij,2))
              ENDDO !nij
            ENDDO !nii
            COEFFSEXT(1)=-SUM(5)-SUM(6)
            COEFFSEXT(2)=-SUM(8)+SUM(4)+4.0d0*SUM(2)
            COEFFSEXT(3)=SUM(5)+SUM(6)
            COEFFSEXT(4)=-SUM(7)+SUM(3)+4.0d0*SUM(1)
            COEFFSEXT(5)=-8.0d0*(SUM(1)+SUM(2))
            COEFFSEXT(6)=SUM(7)-SUM(3)+4.0d0*SUM(1)
            COEFFSEXT(7)=SUM(5)+SUM(6)
            COEFFSEXT(8)=SUM(8)-SUM(4)+4.0d0*SUM(2)
            COEFFSEXT(9)=-SUM(5)-SUM(6)
          ENDIF
        ELSE IF(NITB.EQ.3) THEN
          DO ni=1,15
            SUM(ni)=0.0d0
          ENDDO !ni
          DO nii=1,NITB
C           C(1j)*G1j(upper)
            SUM(1)=SUM(1)+GUQ(1,nii)*PROPQ(nii,1,1,1)
C           C(2j)*G2j(upper)
            SUM(2)=SUM(2)+GUQ(2,nii)*PROPQ(nii,2,1,1)
C           C(3j)*G3j(upper)
            SUM(3)=SUM(3)+GUQ(3,nii)*PROPQ(nii,3,1,1)
C           C(1l)*[Christoffel*Gij(upper)](l)
            SUM(4)=SUM(4)+GCHQ(nii)*PROPQ(nii,1,1,1)
C           C(2l)*[Christoffel*Gij(upper)](l)
            SUM(5)=SUM(5)+GCHQ(nii)*PROPQ(nii,2,1,1)
C           C(3l)*[Christoffel*Gij(upper)](l)
            SUM(6)=SUM(6)+GCHQ(nii)*PROPQ(nii,3,1,1)
            DO nij=1,NITB
C             C(1k),j*Gij(upper)
              SUM(7)=SUM(7)+GUQ(nii,nij)*PROPQ(nii,1,1+nij,1)
C             C(2k),j*Gij(upper)
              SUM(8)=SUM(8)+GUQ(nii,nij)*PROPQ(nii,2,1+nij,1)
C             C(3k),j*Gij(upper)
              SUM(9)=SUM(9)+GUQ(nii,nij)*PROPQ(nii,3,1+nij,1)
            ENDDO !nij
C           C(1j)*G2j (upper)
            SUM(10)=SUM(10)+GUQ(2,nii)*PROPQ(nii,1,1,1)
C           C(2j)*G1j(upper)
            SUM(11)=SUM(11)+GUQ(1,nii)*PROPQ(nii,2,1,1)
C           C(1j)*G3j(upper)
            SUM(12)=SUM(12)+GUQ(3,nii)*PROPQ(nii,1,1,1)
C           C(3j)*G1j(upper)
            SUM(13)=SUM(13)+GUQ(1,nii)*PROPQ(nii,3,1,1)
C           C(2j)*G3j(upper)
            SUM(14)=SUM(14)+GUQ(3,nii)*PROPQ(nii,2,1,1)
C           C(3j)*G2j(upper)
            SUM(15)=SUM(15)+GUQ(2,nii)*PROPQ(nii,3,1,1)
          ENDDO !nii
          NQGW(1)=SUM(14)+SUM(15)
          NQGW(2)=SUM(12)+SUM(13)
          NQGW(3)=-SUM(9)+SUM(6)+4.0d0*SUM(3)
          NQGW(4)=-SUM(12)-SUM(13)
          NQGW(5)=-SUM(14)-SUM(15)
          NQGW(6)=SUM(10)+SUM(11)
          NQGW(7)=-SUM(8)+SUM(5)+4.0d0*SUM(2)
          NQGW(8)=-SUM(10)-SUM(11)
          NQGW(9)=-SUM(7)+SUM(4)+4.0d0*SUM(1)
          NQGW(10)=-8.0d0*(SUM(1)+SUM(2)+SUM(3))
          NQGW(11)=SUM(7)-SUM(4)+4.0d0*SUM(1)
          NQGW(12)=-SUM(10)-SUM(11)
          NQGW(13)=SUM(8)-SUM(5)+4.0d0*SUM(2)
          NQGW(14)=SUM(10)+SUM(11)
          NQGW(15)=-SUM(14)-SUM(15)
          NQGW(16)=-SUM(12)-SUM(13)
          NQGW(17)=SUM(9)-SUM(6)+4.0d0*SUM(3)
          NQGW(18)=SUM(12)+SUM(13)
          NQGW(19)=SUM(14)+SUM(15)
          IF(BIDOMAIN) THEN
            DO ni=1,15
              SUM(ni)=0.0d0
            ENDDO !ni
            DO nii=1,NITB
C             (C(1j)(int)+C(1j)(ext))*G1j(upper)
              SUM(1)=SUM(1)+GUQ(1,nii)*
     '          (PROPQ(nii,1,1,1)+PROPQ(nii,1,1,2))
C             (C(2j)(int)+C(2j)(ext))*G2j(upper)
              SUM(2)=SUM(2)+GUQ(2,nii)*(PROPQ(nii,2,1,1)+
     '          PROPQ(nii,2,1,2))
C             (C(3j)(int)+C(3j)(ext))*G3j(upper)
              SUM(3)=SUM(3)+GUQ(3,nii)*
     '          (PROPQ(nii,3,1,1)+PROPQ(nii,3,1,2))
C             (C(1l)(int)+C(1l)(ext))*[Christoffel*Gij(upper)](l)
              SUM(4)=SUM(4)+GCHQ(nii)*
     '          (PROPQ(nii,1,1,1)+PROPQ(nii,1,1,2))
C             (C(2l)(int)+C(2l)(ext))*[Christoffel*Gij(upper)](l)
              SUM(5)=SUM(5)+GCHQ(nii)*
     '          (PROPQ(nii,2,1,1)+PROPQ(nii,2,1,2))
C             (C(3l)(int)+C(3l)(ext))*[Christoffel*Gij(upper)](l)
              SUM(6)=SUM(6)+GCHQ(nii)*
     '          (PROPQ(nii,3,1,1)+PROPQ(nii,3,1,2))
              DO nij=1,NITB
C               (C(1k),j(int)+C(1k),j(ext))*Gij(upper)
                SUM(7)=SUM(7)+GUQ(nii,nij)*
     '            (PROPQ(nii,1,1+nij,1)+PROPQ(nii,1,1+nij,2))
C               (C(2k),j(int)+C(2k),j(ext))*Gij(upper)
                SUM(8)=SUM(8)+GUQ(nii,nij)*
     '            (PROPQ(nii,2,1+nij,1)+PROPQ(nii,2,1+nij,2))
C               (C(3k),j(int)+C(3k),j(ext))*Gij(upper)
                SUM(9)=SUM(9)+GUQ(nii,nij)*
     '            (PROPQ(nii,3,1+nij,1)+PROPQ(nii,3,1+nij,2))
              ENDDO !nij
C             (C(1j)(int)+C(1j)(ext))*G2j(upper)
              SUM(10)=SUM(10)+GUQ(2,nii)*
     '          (PROPQ(nii,1,1,1)+PROPQ(nii,1,1,2))
C             (C(2j)(int)+C(2j)(ext))*G1j(upper)
              SUM(11)=SUM(11)+GUQ(1,nii)*
     '          (PROPQ(nii,2,1,1)+PROPQ(nii,2,1,2))
C             (C(1j)(int)+C(1j)(ext))*G3j(upper)
              SUM(12)=SUM(12)+GUQ(3,nii)*
     '          (PROPQ(nii,1,1,1)+PROPQ(nii,1,1,2))
C             (C(3j)(int)+C(3j)(ext))*G1j(upper)
              SUM(13)=SUM(13)+GUQ(1,nii)*
     '          (PROPQ(nii,3,1,1)+PROPQ(nii,3,1,2))
C             (C(2j)*int)+C(2j)(ext)*G3j(upper)
              SUM(14)=SUM(14)+GUQ(3,nii)*
     '          (PROPQ(nii,2,1,1)+ PROPQ(nii,2,1,2))
C             (C(3j)(int)+C(3j)(ext))*G2j(upper)
              SUM(15)=SUM(15)+GUQ(2,nii)*
     '          (PROPQ(nii,3,1,1)+PROPQ(nii,3,1,2))
            ENDDO !nii
C??? Should these have there signs reversed ???
            COEFFSEXT(1)=-SUM(14)-SUM(15)
            COEFFSEXT(2)=-SUM(12)-SUM(13)
            COEFFSEXT(3)=SUM(9)-SUM(6)-4.0d0*SUM(3)
            COEFFSEXT(4)=SUM(12)+SUM(13)
            COEFFSEXT(5)=SUM(14)+SUM(15)
            COEFFSEXT(6)=-SUM(10)-SUM(11)
            COEFFSEXT(7)=SUM(8)-SUM(5)-4.0d0*SUM(2)
            COEFFSEXT(8)=SUM(10)+SUM(11)
            COEFFSEXT(9)=SUM(7)-SUM(4)-4.0d0*SUM(1)
            COEFFSEXT(10)=8.0d0*(SUM(1)+SUM(2)+SUM(3))
            COEFFSEXT(11)=-SUM(7)+SUM(4)-4.0d0*SUM(1)
            COEFFSEXT(12)=SUM(10)+SUM(11)
            COEFFSEXT(13)=-SUM(8)+SUM(5)-4.0d0*SUM(2)
            COEFFSEXT(14)=-SUM(10)-SUM(11)
            COEFFSEXT(15)=SUM(14)+SUM(15)
            COEFFSEXT(16)=SUM(12)+SUM(13)
            COEFFSEXT(17)=-SUM(9)+SUM(6)-4.0d0*SUM(3)
            COEFFSEXT(18)=-SUM(12)-SUM(13)
            COEFFSEXT(19)=-SUM(14)-SUM(15)
          ENDIF
        ENDIF
      ELSE !external grid point
        IF(IMPLICIT) THEN
          IF(FIXQ(nq,1,nx_trans)) THEN
            NQGW(1)=1.0d0
            NQGW(2)=0.0d0
            NQGW(3)=0.0d0         
          ELSEIF(FIXQ(nq,2,nx_trans)) THEN
            NQGW(1)=1.0d0
            NQGW(2)=-4.0d0/3.0d0
            NQGW(3)=1.0d0/3.0d0
          ELSE
            WRITE(ERROR,'(''>>No boundary condition set at point'
     '        //' (trans.) '',I8)') nq
            GOTO 9999
          ENDIF
        ENDIF
        IF(BIDOMAIN) THEN
          IF(FIXQ(nq,1,nx_ext)) THEN
            COEFFSEXT(1)=1.0d0
            COEFFSEXT(2)=0.0d0
            COEFFSEXT(3)=0.0d0         
          ELSEIF(FIXQ(nq,2,nx_ext)) THEN
            COEFFSEXT(1)=1.0d0
            COEFFSEXT(2)=-4.0d0/3.0d0
            COEFFSEXT(3)=1.0d0/3.0d0
          ELSEIF(FIXQ(nq,3,nx_ext)) THEN
            COEFFSEXT(1)=1.0d0
            COEFFSEXT(2)=0.0d0
            COEFFSEXT(3)=0.0d0
          ELSE
            WRITE(ERROR,'(''>>No boundary condition set at point'
     '        //' (ext.) '',I8)') nq
            GOTO 9999
          ENDIF
        ENDIF
      ENDIF
      
      CALL EXITS('CALC_ADAMS_GRID_COEF')
      RETURN
 9999 CALL ERRORS('CALC_ADAMS_GRID_COEF',ERROR)
      CALL EXITS('CALC_ADAMS_GRID_COEF')
      RETURN 1
      END      


      SUBROUTINE CALC_ADAMS_PSEUDO_STIMULUS(niqV,nq,NQGP,NQGP_PIVOT,
     '  nx_ext,nx_trans,NQGW,PSTIMULUS,THETA,YQ,BIDOMAIN,ERROR,*)

C#### Subroutine: CALC_ADAMS_PSEUDO_STIMULUS
C###  Description:
C###    CALC_ADAMS_PSEUDO_STIMULUS calculates the pseudo stimulus
C###    currents to be applied to a cell. This pseudo stimulus
C###    currents is made up of the explicit component of the 
C###    intracellular diffusive current and the contribution
C###    from the bidomain extracellular diffusive current if they
C###    exist.

      IMPLICIT NONE

      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      
!     Parameter list
      INTEGER niqV,nq,NQGP(0:19),NQGP_PIVOT(19),nx_ext,nx_trans
      REAL*8 NQGW(19),PSTIMULUS,THETA,YQ(NYQM,NIQM,NAM,NXM)
      LOGICAL BIDOMAIN
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER i

      CALL ENTERS('CALC_ADAMS_PSEUDO_STIMULUS',*9999)

      PSTIMULUS=0.0d0
      IF(THETA.LT.1.0d0) THEN
        DO i=1,NQGP(0)
          PSTIMULUS=PSTIMULUS+(NQGW(NQGP_PIVOT(i))*
     '      YQ(NQGP(i),niqV,1,nx_trans))
        ENDDO !i
        PSTIMULUS=PSTIMULUS*(1.0d0-THETA)
      ENDIF
      IF(BIDOMAIN) THEN
        DO i=1,NQGP(0)
          PSTIMULUS=PSTIMULUS+(NQGW(NQGP_PIVOT(i))*
     '      YQ(NQGP(i),niqV,1,nx_ext))
        ENDDO !i        
      ENDIF
      
      IF(DOP) THEN
C$OMP CRITICAL(CALC_ADAMS_PSEUDO_STIMULUS_1)
        WRITE(OP_STRING,'('' Pseudo stimulus = '',F11.6,'' uA/mm^3'')')
     '    PSTIMULUS
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C$OMP END CRITICAL(CALC_ADAMS_PSEUDO_STIMULUS_1)
      ENDIF

      CALL EXITS('CALC_ADAMS_PSEUDO_STIMULUS')
      RETURN
 9999 CALL ERRORS('CALC_ADAMS_PSEUDO_STIMULUS',ERROR)
      CALL EXITS('CALC_ADAMS_PSEUDO_STIMULUS')
      RETURN 1
      END


      SUBROUTINE CALC_GRID_BOUND_COEF(APPROX,NENQ,nq,NQGP,NQGP_PIVOT,
     '  NQS,NQXI,NXQ,COEFF,CQ,DNUDXQ,DXDXIQ,ERROR,*)

C#### Subroutine: CALC_GRID_BOUND_COEF
C###  Description:
C###    CALC_GRID_BOUND_COEF calculates the flux coefficients of Phi
C###    and the grid points involved in calculating the flux 
C###    coefficients. The type of approximation which is used for
C###    calculating the coefficients is determined by the APPROX
C###    variable which chooses between 2,3,4 point one sided and 2,4
C###    point two sided grid point schemes.
C***  Created by Martin Buist, January 1999

      IMPLICIT NONE

      INCLUDE 'cmiss$reference:geom00.cmn'

!     Parameter list
      INTEGER APPROX,NENQ(0:8,NQM),nq,NQGP(0:19,NQM),NQGP_PIVOT(19,NQM),
     '  NQS(NEM),NQXI(0:NIM,NQSCM),NXQ(-NIM:NIM,0:4,0:NQM)
      REAL*8 COEFF(19),CQ(NMM),DNUDXQ(3,3,NQM),DXDXIQ(3,3,NQM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER i,ne,ni,nj,nj1,nj2,nq1,nq2,nq3,nq4,SCHEME
      REAL*8 DET,DNUDX(3,3),DXDNU(3,3),DXDXI(3,3),DXI,DXIDX(3,3),
     '  G(3),XNLOCAL(3)
      LOGICAL FOUND

      CALL ENTERS('CALC_GRID_BOUND_COEF',*9999)
      
      ne=NENQ(1,nq)
      SCHEME=NQS(ne)
      NQGP(0,nq)=0
      DO i=1,19
        COEFF(i)=0.0d0
      ENDDO !i

      !calculate dxi/dx
      DO nj=1,NJT
        DO ni=1,NQXI(0,SCHEME)
          DXDXI(nj,ni)=DXDXIQ(nj,ni,nq)
        ENDDO !ni
      ENDDO !nj
      CALL ASSERT(NJT.EQ.NQXI(0,SCHEME),' >>NJT.NE.NIT - cannot invert',
     '  ERROR,*9999)
      CALL INVERT(NJT,DXDXI,DXIDX,DET)

      !calculate g(nj) conductivities
      DO nj2=1,NJT
        DO nj1=1,NJT
          DNUDX(nj1,nj2)=DNUDXQ(nj1,nj2,nq)
        ENDDO !nj1
      ENDDO !nj2
      CALL INVERT(NJT,DNUDX,DXDNU,DET)
      DO nj2=1,NJT
        G(nj2)=0.0d0
        DO nj1=1,NJT
          G(nj2)=G(nj2)+(DNUDX(nj1,nj2)*DXDNU(nj2,nj1)*CQ(nj1))
        ENDDO !nj1
        G(nj2)=G(nj2)*1000.0d0 !change to mS from Siemens
      ENDDO !nj2

      !calculate a normal vector at nq
      CALL NORM31(nq,NXQ,DXDXIQ,XNLOCAL,ERROR,*9999)

      !calculate coefficients and grid points
      DXI=0.5d0
      DO ni=1,NQXI(0,SCHEME)
        IF(NXQ(-ni,1,nq).EQ.0) THEN
          !no grid point in -ni
          !use one sided difference
          IF(APPROX.EQ.1) THEN
            nq1=NXQ(ni,1,nq)

            NQGP(0,nq)=NQGP(0,nq)+1
            NQGP(NQGP(0,nq),nq)=nq
            DO nj=1,NJT
              COEFF(NQGP(0,nq))=COEFF(NQGP(0,nq))-((XNLOCAL(nj)*G(nj)*
     '          DXIDX(ni,nj))/DXI)
            ENDDO !nj
            NQGP(0,nq)=NQGP(0,nq)+1
            NQGP(NQGP(0,nq),nq)=nq1
            DO nj=1,NJT
              COEFF(NQGP(0,nq))=COEFF(NQGP(0,nq))+((XNLOCAL(nj)*G(nj)*
     '          DXIDX(ni,nj))/DXI)
            ENDDO !nj
C            DPHIDXI(ni)=-(YQ(nq,niqV)-YQ(nq1,niqV))/DXI
          ELSEIF(APPROX.EQ.2) THEN
            nq1=NXQ(ni,1,nq)
            nq2=NXQ(ni,1,nq1)

            NQGP(0,nq)=NQGP(0,nq)+1
            NQGP(NQGP(0,nq),nq)=nq
            DO nj=1,NJT
              COEFF(NQGP(0,nq))=COEFF(NQGP(0,nq))-((3.0d0*XNLOCAL(nj)*
     '          G(nj)*DXIDX(ni,nj))/(2.0d0*DXI))
            ENDDO !nj
            NQGP(0,nq)=NQGP(0,nq)+1
            NQGP(NQGP(0,nq),nq)=nq1
            DO nj=1,NJT
              COEFF(NQGP(0,nq))=COEFF(NQGP(0,nq))+((4.0d0*XNLOCAL(nj)*
     '          G(nj)*DXIDX(ni,nj))/(2.0d0*DXI))
            ENDDO !nj
            NQGP(0,nq)=NQGP(0,nq)+1
            NQGP(NQGP(0,nq),nq)=nq2
            DO nj=1,NJT
              COEFF(NQGP(0,nq))=COEFF(NQGP(0,nq))-((1.0d0*XNLOCAL(nj)*
     '          G(nj)*DXIDX(ni,nj))/(2.0d0*DXI))
            ENDDO !nj

C            DPHIDXI(ni)=-((3.0d0*YQ(nq,niqV))-(4.0d0*YQ(nq1,niqV))+
C     '        (1.0d0*YQ(nq2,niqV)))/(2.0d0*DXI)
          ELSEIF(APPROX.EQ.3) THEN
            nq1=NXQ(ni,1,nq)
            nq2=NXQ(ni,1,nq1)
            nq3=NXQ(ni,1,nq2)

            NQGP(0,nq)=NQGP(0,nq)+1
            NQGP(NQGP(0,nq),nq)=nq
            DO nj=1,NJT
              COEFF(NQGP(0,nq))=COEFF(NQGP(0,nq))-((11.0d0*XNLOCAL(nj)*
     '          G(nj)*DXIDX(ni,nj))/(6.0d0*DXI))
            ENDDO !nj
            NQGP(0,nq)=NQGP(0,nq)+1
            NQGP(NQGP(0,nq),nq)=nq1
            DO nj=1,NJT
              COEFF(NQGP(0,nq))=COEFF(NQGP(0,nq))+((18.0d0*XNLOCAL(nj)*
     '          G(nj)*DXIDX(ni,nj))/(6.0d0*DXI))
            ENDDO !nj
            NQGP(0,nq)=NQGP(0,nq)+1
            NQGP(NQGP(0,nq),nq)=nq2
            DO nj=1,NJT
              COEFF(NQGP(0,nq))=COEFF(NQGP(0,nq))-((9.0d0*XNLOCAL(nj)*
     '          G(nj)*DXIDX(ni,nj))/(6.0d0*DXI))
            ENDDO !nj
            NQGP(0,nq)=NQGP(0,nq)+1
            NQGP(NQGP(0,nq),nq)=nq3
            DO nj=1,NJT
              COEFF(NQGP(0,nq))=COEFF(NQGP(0,nq))+((2.0d0*XNLOCAL(nj)*
     '          G(nj)*DXIDX(ni,nj))/(6.0d0*DXI))
            ENDDO !nj

C            DPHIDXI(ni)=-((11.0d0*YQ(nq,niqV))-(18.0d0*YQ(nq1,niqV))+
C     '        (9.0d0*YQ(nq2,niqV))-(2.0d0*YQ(nq3,niqV)))/
C     '        (6.0d0*DXI)
          ENDIF !approx
        ELSEIF(NXQ(ni,1,nq).EQ.0) THEN
          !no grid point in +ni
          !use one sided difference
          IF(APPROX.EQ.1) THEN
            nq1=NXQ(-ni,1,nq)

            NQGP(0,nq)=NQGP(0,nq)+1
            NQGP(NQGP(0,nq),nq)=nq
            DO nj=1,NJT
              COEFF(NQGP(0,nq))=COEFF(NQGP(0,nq))+((XNLOCAL(nj)*G(nj)*
     '          DXIDX(ni,nj))/DXI)
            ENDDO !nj
            NQGP(0,nq)=NQGP(0,nq)+1
            NQGP(NQGP(0,nq),nq)=nq1
            DO nj=1,NJT
              COEFF(NQGP(0,nq))=COEFF(NQGP(0,nq))-((XNLOCAL(nj)*G(nj)*
     '          DXIDX(ni,nj))/DXI)
            ENDDO !nj
C            DPHIDXI(ni)=-(YQ(nq,niqV)-YQ(nq1,niqV))/DXI
          ELSEIF(APPROX.EQ.2) THEN
            nq1=NXQ(-ni,1,nq)
            nq2=NXQ(-ni,1,nq1)

            NQGP(0,nq)=NQGP(0,nq)+1
            NQGP(NQGP(0,nq),nq)=nq
            DO nj=1,NJT
              COEFF(NQGP(0,nq))=COEFF(NQGP(0,nq))+((3.0d0*XNLOCAL(nj)*
     '          G(nj)*DXIDX(ni,nj))/(2.0d0*DXI))
            ENDDO !nj
            NQGP(0,nq)=NQGP(0,nq)+1
            NQGP(NQGP(0,nq),nq)=nq1
            DO nj=1,NJT
              COEFF(NQGP(0,nq))=COEFF(NQGP(0,nq))-((4.0d0*XNLOCAL(nj)*
     '          G(nj)*DXIDX(ni,nj))/(2.0d0*DXI))
            ENDDO !nj
            NQGP(0,nq)=NQGP(0,nq)+1
            NQGP(NQGP(0,nq),nq)=nq2
            DO nj=1,NJT
              COEFF(NQGP(0,nq))=COEFF(NQGP(0,nq))+((1.0d0*XNLOCAL(nj)*
     '          G(nj)*DXIDX(ni,nj))/(2.0d0*DXI))
            ENDDO !nj
C            DPHIDXI(ni)=-((3.0d0*YQ(nq,niqV))-(4.0d0*YQ(nq1,niqV))+
C     '        (1.0d0*YQ(nq2,niqV)))/(2.0d0*DXI)
          ELSEIF(APPROX.EQ.3) THEN
            nq1=NXQ(-ni,1,nq)
            nq2=NXQ(-ni,1,nq1)
            nq3=NXQ(-ni,1,nq2)

            NQGP(0,nq)=NQGP(0,nq)+1
            NQGP(NQGP(0,nq),nq)=nq
            DO nj=1,NJT
              COEFF(NQGP(0,nq))=COEFF(NQGP(0,nq))+((11.0d0*XNLOCAL(nj)*
     '          G(nj)*DXIDX(ni,nj))/(6.0d0*DXI))
            ENDDO !nj
            NQGP(0,nq)=NQGP(0,nq)+1
            NQGP(NQGP(0,nq),nq)=nq1
            DO nj=1,NJT
              COEFF(NQGP(0,nq))=COEFF(NQGP(0,nq))-((18.0d0*XNLOCAL(nj)*
     '          G(nj)*DXIDX(ni,nj))/(6.0d0*DXI))
            ENDDO !nj
            NQGP(0,nq)=NQGP(0,nq)+1
            NQGP(NQGP(0,nq),nq)=nq2
            DO nj=1,NJT
              COEFF(NQGP(0,nq))=COEFF(NQGP(0,nq))+((9.0d0*XNLOCAL(nj)*
     '          G(nj)*DXIDX(ni,nj))/(6.0d0*DXI))
            ENDDO !nj
            NQGP(0,nq)=NQGP(0,nq)+1
            NQGP(NQGP(0,nq),nq)=nq3
            DO nj=1,NJT
              COEFF(NQGP(0,nq))=COEFF(NQGP(0,nq))-((2.0d0*XNLOCAL(nj)*
     '          G(nj)*DXIDX(ni,nj))/(6.0d0*DXI))
            ENDDO !nj
C            DPHIDXI(ni)=-((11.0d0*YQ(nq,niqV))-(18.0d0*YQ(nq1,niqV))+
C     '        (9.0d0*YQ(nq2,niqV))-(2.0d0*YQ(nq3,niqV)))/
C     '        (6.0d0*DXI)
          ENDIF !approx
        ELSE
          !points in both +ni and -ni
          !we can use a two sided difference
          IF((APPROX.EQ.1).OR.(APPROX.EQ.2)) THEN
            nq1=NXQ(-ni,1,nq)
            nq2=NXQ(ni,1,nq)

            NQGP(0,nq)=NQGP(0,nq)+1
            NQGP(NQGP(0,nq),nq)=nq1
            DO nj=1,NJT
              COEFF(NQGP(0,nq))=COEFF(NQGP(0,nq))-((XNLOCAL(nj)*G(nj)*
     '          DXIDX(ni,nj))/(2.0d0*DXI))
            ENDDO !nj
            NQGP(0,nq)=NQGP(0,nq)+1
            NQGP(NQGP(0,nq),nq)=nq2
            DO nj=1,NJT
              COEFF(NQGP(0,nq))=COEFF(NQGP(0,nq))+((XNLOCAL(nj)*G(nj)*
     '          DXIDX(ni,nj))/(2.0d0*DXI))
            ENDDO !nj
C            DPHIDXI(ni)=(YQ(nq2,niqV)-YQ(nq1,niqV))/(2.0d0*DXI)
          ELSEIF(APPROX.EQ.3) THEN
            nq1=NXQ(-ni,1,nq)
            nq2=NXQ(ni,1,nq)
            nq3=NXQ(-ni,1,nq1)
            IF((nq3.LE.0).OR.(nq3.GT.NQT)) nq3=nq1
            nq4=NXQ(ni,1,nq2)
            IF((nq4.LE.0).OR.(nq4.GT.NQT)) nq4=nq1

            NQGP(0,nq)=NQGP(0,nq)+1
            NQGP(NQGP(0,nq),nq)=nq1
            DO nj=1,NJT
              COEFF(NQGP(0,nq))=COEFF(NQGP(0,nq))-((8.0d0*XNLOCAL(nj)*
     '          G(nj)*DXIDX(ni,nj))/(12.0d0*DXI))
            ENDDO !nj
            NQGP(0,nq)=NQGP(0,nq)+1
            NQGP(NQGP(0,nq),nq)=nq2
            DO nj=1,NJT
              COEFF(NQGP(0,nq))=COEFF(NQGP(0,nq))+((8.0d0*XNLOCAL(nj)*
     '          G(nj)*DXIDX(ni,nj))/(12.0d0*DXI))
            ENDDO !nj
            NQGP(0,nq)=NQGP(0,nq)+1
            NQGP(NQGP(0,nq),nq)=nq3
            DO nj=1,NJT
              COEFF(NQGP(0,nq))=COEFF(NQGP(0,nq))+((1.0d0*XNLOCAL(nj)*
     '          G(nj)*DXIDX(ni,nj))/(12.0d0*DXI))
            ENDDO !nj
            NQGP(0,nq)=NQGP(0,nq)+1
            NQGP(NQGP(0,nq),nq)=nq4
            DO nj=1,NJT
              COEFF(NQGP(0,nq))=COEFF(NQGP(0,nq))-((1.0d0*XNLOCAL(nj)*
     '          G(nj)*DXIDX(ni,nj))/(12.0d0*DXI))
            ENDDO !nj
C            DPHIDXI(ni)=(YQ(nq3,niqV)-(8.0d0*YQ(nq1,niqV))+(8.0d0*
C     '        YQ(nq2,niqV))-YQ(nq4,niqV))/(12.0d0*DXI)
          ENDIF !approx
        ENDIF !NXQ non-zero
      ENDDO !ni

      FOUND=.FALSE.
      DO nq1=1,NQGP(0,nq)
        DO nq2=1,NQGP(0,nq)
          IF((nq1.NE.nq2).AND.(NQGP(nq1,nq).EQ.NQGP(nq2,nq))) THEN
            !duplicate grid point entry
            nq3=nq1
            nq4=nq2
            FOUND=.TRUE.
          ENDIF !duplicate
        ENDDO !nq2
      ENDDO !nq1

      IF(FOUND) THEN
        COEFF(nq4)=COEFF(nq4)+COEFF(nq3)
        NQGP(0,nq)=NQGP(0,nq)-1
        DO nq1=nq3,NQGP(0,nq)
          COEFF(nq1)=COEFF(nq1+1)
          NQGP(nq1,nq)=NQGP(nq1+1,nq)
        ENDDO !nq1
      ENDIF !found

      CALL ISORTP(NQGP(0,nq),NQGP(1,nq),NQGP_PIVOT(1,nq))

      CALL EXITS('CALC_GRID_BOUND_COEF')
      RETURN
 9999 CALL ERRORS('CALC_GRID_BOUND_COEF',ERROR)
      CALL EXITS('CALC_GRID_BOUND_COEF')
      RETURN 1
      END


      SUBROUTINE GEN_INT_RHS(AQ,CQ,IION,maqp1t0,maqp1t1,maqp1i,maqp2t0,
     '  maqp2t1,maqp2i,niqBNDRY,niqSAC,niqV,NQGP,NQGP_PIVOT,NQGW,NRLIST,
     '  NWQ,nx_ext,nx_trans,RHS,T,YQ,FIXQ,ERROR,*)

C#### Subroutine: GEN_INT_RHS
C###  Description:
C###    GEN_INT_RHS is used to generate the right hand side vectors
C###    at each time step for the implicit solution of grid
C###    activation problems. (intracellular for bidomain, or monodomain)
C***  Created by Martin Buist, May 1997

      IMPLICIT NONE

      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b12.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:grid00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:ktyp30.cmn'
      INCLUDE 'cmiss$reference:tol00.cmn'

!     Parameter list
      INTEGER maqp1t0,maqp1t1,maqp1i,maqp2t0,maqp2t1,maqp2i,
     '  niqBNDRY,niqSAC,niqV,NQGP(0:19,NQM),NQGP_PIVOT(19,NQM),
     '  NRLIST(0:NRM),NWQ(6,0:NQM),nx_ext,nx_trans
      REAL*8 AQ(NMAQM,NQM),CQ(NMM,NQM),NQGW(NQM,19),RHS(NQM),T,
     '  YQ(NYQM,NIQM,NAM,NXM),IION(NQM)
      LOGICAL FIXQ(NYQM,NIYFIXM,NXM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER nq,nr,nrr
      REAL*8 DTCMAM,DTCM,IAPP
      LOGICAL SAC

      CALL ENTERS('GEN_INT_RHS',*9999)

! Set up extracellular contribution to RHS      
      IF(KTYP32.EQ.2) THEN !bidomain
        CALL GET_EXT_CONTRIB(niqV,NQGP,NQGP_PIVOT,NQGW,NRLIST,NWQ,
     '    nx_ext,RHS,YQ,ERROR,*9999)
      ENDIF

      DO nrr=1,NRLIST(0)
        nr=NRLIST(nrr)
        DO nq=NQR(1,nr),NQR(2,nr)

          SAC=.FALSE.
          IF(ITYP3(nr,nx_trans).EQ.2) THEN      !FHN
            IF(DABS(CQ(18,nq)).GT.ZERO_TOL) SAC=.TRUE.
          ELSEIF(ITYP3(nr,nx_trans).EQ.3) THEN  !VCD
            IF(DABS(CQ(13,nq)).GT.ZERO_TOL) SAC=.TRUE.
          ENDIF

! Set up current applied (intracellular)
          IAPP=0.0d0
          IF((AQ(maqp1i,nq).GT.ZERO_TOL).AND.(T.GE.AQ(maqp1t0,nq))
     '      .AND.(T.LT.AQ(maqp1t1,nq))) THEN
            IAPP=AQ(maqp1i,nq)
          ELSE IF((AQ(maqp2i,nq).GT.ZERO_TOL).AND.(T.GE.
     '      AQ(maqp2t0,nq)).AND.(T.LT.AQ(maqp2t1,nq))) THEN
            IAPP=AQ(maqp2i,nq)
          ENDIF

          IF(SAC) THEN
            CALL ASSERT(NIQM.GE.7,'>>Increase NIQM (>=7)',ERROR,*9999)
!            IF(YQ(nq,7,1,nx_trans).GT.ZERO_TOL) THEN !inward SAC current
              IAPP=IAPP-YQ(nq,niqSAC,1,nx_trans)
!            ENDIF
          ENDIF

          DTCM=DT/(CQ(1,nq)*1.0D-6) !10^-6 correction for Cm
          DTCMAM=DTCM/CQ(2,nq)

          IF(NWQ(1,nq).EQ.0) THEN !internal
            IF(KTYP32.EQ.2) THEN !bidomain
              RHS(nq)=(DTCMAM*RHS(nq))+YQ(nq,niqV,1,nx_trans)+
     '          (DTCMAM*IAPP)-(DTCM*IION(nq))
            ELSE !monodomain
              RHS(nq)=YQ(nq,niqV,1,nx_trans)+(DTCMAM*IAPP)-
     '          (DTCM*IION(nq))
            ENDIF
          ELSE
            IF(FIXQ(nq,1,nx_trans)) THEN
              RHS(nq)=YQ(nq,niqBNDRY,1,nx_trans)
            ELSEIF(FIXQ(nq,2,nx_trans)) THEN  
              RHS(nq)=YQ(nq,niqBNDRY,1,nx_trans)
            ELSE
              WRITE(OP_STRING,'(''>>No boundary condition set at'
     '          //' point '',I8)') nq
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              GOTO 9999
            ENDIF
          ENDIF

        ENDDO !nq
      ENDDO !nr

      IF(DOP) THEN
        WRITE(OP_STRING,'(''INT_RHS:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nq=1,NQT
          WRITE(OP_STRING,'(I8,F12.6)') nq,RHS(nq)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

      CALL EXITS('GEN_INT_RHS')
      RETURN
 9999 CALL ERRORS('GEN_INT_RHS',ERROR)
      CALL EXITS('GEN_INT_RHS')
      RETURN 1
      END


      SUBROUTINE GET_EXT_CONTRIB(niqV,NQGP,NQGP_PIVOT,NQGW,NRLIST,NWQ,
     '  nx_ext,RHS,YQ,ERROR,*)

C#### Subroutine: GET_EXT_CONTRIB
C###  Description:
C###    GET_EXT_CONTRIB returns the contribution of the extracellular
C###    potential to the transmembrane potential equation. This is
C###    only used in bidomain simulations, the extracellular 
C###    contribution is zero in a monodomain simulation.
C***  Created by Martin Buist, August 1997

      IMPLICIT NONE

      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:grid00.cmn'

!     Parameter list
      INTEGER niqV,NQGP(0:19,NQM),NQGP_PIVOT(19,NQM),
     '  NRLIST(0:NRM),NWQ(6,0:NQM),nx_ext
      REAL*8 NQGW(NQM,19),RHS(NQM),YQ(NYQM,NIQM,NAM,NXM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER i,nq,nr,nrr

      CALL ENTERS('GET_EXT_CONTRIB',*9999)

      DO nrr=1,NRLIST(0)
        nr=NRLIST(nrr)
!start of parallel region
C$DOACROSS local(nq,i), share(NWQ,RHS,NQGP,NQGW,YQ)
        DO nq=NQR(1,nr),NQR(2,nr)
          IF(NWQ(1,nq).EQ.0) THEN !internal
            RHS(nq)=0.0d0
            DO i=1,NQGP(0,nq)
              RHS(nq)=RHS(nq)-(NQGW(nq,NQGP_PIVOT(i,nq))*
     '          YQ(NQGP(i,nq),niqV,1,nx_ext))
            ENDDO !i
          ELSE
            RHS(nq)=0.0d0
          ENDIF
        ENDDO !nq
!end of parallel region
      ENDDO !nr 

      CALL EXITS('GET_EXT_CONTRIB')
      RETURN
 9999 CALL ERRORS('GET_EXT_CONTRIB',ERROR)
      CALL EXITS('GET_EXT_CONTRIB')
      RETURN 1
      END


      SUBROUTINE IONIC_CURRENT(CQ,IION,niqCAI,niqCALCIUM,niqD,niqDRECOV,
     '  niqDV,niqF,niqH,niqJ,niqMM,niqOLDSOLN,niqRECOV,niqSAC,niqV,
     '  niqX,NRLIST,nx,YQ,ERROR,*)

C#### Subroutine: IONIC_CURRENT
C###  Description:
C###    IONIC_CURRENT returns the value of the ionic current
C###    from the transmembrane potential for a grid point nq.
C###    recovery variables are updated if they are required.
C***  Created by Martin Buist 1 August 1997

      IMPLICIT NONE

      INCLUDE 'cmiss$reference:b12.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:grid00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:ktyp30.cmn'

!     Parameter list
      INTEGER niqCAI,niqCALCIUM,niqD,niqDRECOV,niqDV,niqF,niqH,niqJ,
     '  niqMM,niqOLDSOLN,niqRECOV,niqSAC,niqV,niqX,NRLIST(0:NRM),nx
      REAL*8 CQ(NMM,NQM),IION(NQM),YQ(NYQM,NIQM,NAM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER nq,nr,nrr
      REAL*8 CUBIC_ION,D_RECOV,DIFF,EPS,FHN_ION,LR_ION,
     '  NORM_PHI,PHIM,RATE,RECOV,T_CONST,V,VCD_ION,VSTAR

      CALL ENTERS('IONIC_CURRENT',*9999)

      niqOLDSOLN=niqOLDSOLN
      niqSAC=niqSAC

      DO nrr=1,NRLIST(0)
        nr=NRLIST(nrr)
!start of parallel region
C$DOACROSS local(nq,PHIM,V,RECOV,DIFF,NORM_PHI,D_RECOV,EPS,RATE,
C$& VSTAR,T_CONST), share(IION,YQ,CQ,nr,nx)
        DO nq=NQR(1,nr),NQR(2,nr)
          V=YQ(nq,niqV,1)
! Calculate the rate constants and/or update the recovery variables

          IF(ITYP3(nr,nx).EQ.2) THEN !FHN
            RECOV=YQ(nq,niqRECOV,1)
! CQ( 9,nq) is the rest potential
! CQ(10,nq) is the plateau potential      
! CQ(11,nq) is the threshold potential      
! NORM_PHI  is a 0-1 normalised potential
            DIFF = CQ(10,nq)-CQ(9,nq)
            NORM_PHI=(V-CQ(9,nq))/DIFF
            IF(KTYP33.LT.3) THEN !Standard or Roger's FHN
! CQ(14,nq) is the recovery rate constant
! CQ(15,nq) is the recovery decay constant
              D_RECOV=CQ(14,nq)*(NORM_PHI-CQ(15,nq)*RECOV)
            ELSE !Panfilov FHN
! CQ(12,nq) is rate constant (same because of identical parabola)
! CQ(14,nq) is epsilon0
! CQ(15,nq) is mu1
! CQ(16,nq) is mu2
! CQ(17,nq) is tau (time const for force)
              EPS=CQ(14,nq)+CQ(15,nq)*RECOV/(NORM_PHI+CQ(16,nq))
              RATE=CQ(12,nq)/DIFF/DIFF
              VSTAR=V+CQ(9,nq)-CQ(10,nq)-CQ(11,nq)
              D_RECOV=EPS*(-RECOV-RATE*VSTAR*(V-CQ(9,nq)))
            ENDIF
            YQ(nq,niqRECOV,1)=RECOV+DT*D_RECOV
! Compute calcium level (simple equation from Panfilov)
            IF(DABS(CQ(17,nq)).GT.1d-12) THEN
              IF(NORM_PHI.LT.0.d0) NORM_PHI=0.d0
              YQ(nq,niqCALCIUM,1)=YQ(nq,niqCALCIUM,1)+DT*
     '          (NORM_PHI-YQ(nq,niqCALCIUM,1))/CQ(17,nq)
            ELSE
              YQ(nq,niqCALCIUM,1)=0.d0
            ENDIF

          ELSE IF(ITYP3(nr,nx).EQ.3) THEN !vanCapelle-Durrer
            RECOV=YQ(nq,niqRECOV,1)
            IF(KTYP33.EQ.1) THEN !Original VCD
              T_CONST=CQ(10,nq)
            ELSE IF(KTYP33.EQ.2) THEN !VCDC mods
! From modifications provided by Alan Garfinkel
              IF(YQ(nq,niqDRECOV,1).ge.0.0d0) THEN  
                !recovery var Y is increasing
                IF(KTYP34.EQ.1) THEN !"Normal" APD (209ms)
                  T_CONST=0.5d0
                ELSE IF(KTYP34.EQ.2) THEN !"Ischemic" APD (112ms)
                  T_CONST=0.33d0
                ENDIF
              ELSE IF(recov.gt.0.85d0) THEN !Y decreasing & Y>0.85
                IF(KTYP34.EQ.1) THEN !"Normal" APD (209ms)
                  T_CONST=0.1d0
                ELSE IF(KTYP34.EQ.2) THEN !"Ischemic" APD (112ms)
                  T_CONST=0.066d0
                ENDIF
              ELSE                          !Y decreasing & Y<0.85
                IF(KTYP34.EQ.1) THEN !"Normal" APD (209ms)
                  T_CONST=3.0d0
                ELSE IF(KTYP34.EQ.2) THEN !"Ischemic" APD (112ms)
!PJH 5/7/98       T_CONST=3.31d0
                  T_CONST=CQ(14,nq)
                ENDIF
              ENDIF
! Scale Factor for time constant
              T_CONST=T_CONST*CQ(10,nq) !CQ(10,nq) is time constant
            ENDIF
            IF(V.LT.-80.0d0) THEN
              D_RECOV=-RECOV/T_CONST
            ELSE IF(V.GT.-60.0d0) THEN
              D_RECOV=(1.0d0-RECOV)/T_CONST
            ELSE
              D_RECOV=((V+80.0d0)/20.0d0-RECOV)/T_CONST
            ENDIF
            YQ(nq,niqRECOV,1)=RECOV+DT*D_RECOV
            IF(KTYP33.EQ.2) THEN !VCDC mods
              YQ(nq,niqDRECOV,1)=D_RECOV
            ENDIF
! Compute calcium level (simple equation from Panfilov)
            IF(DABS(CQ(12,nq)).GT.1d-12) THEN
              NORM_PHI=(V-CQ(9,nq))/100.d0
              IF(NORM_PHI.LT.0.d0) NORM_PHI=0.d0
              YQ(nq,niqCALCIUM,1)=YQ(nq,niqCALCIUM,1)+DT*
     '          (NORM_PHI-YQ(nq,niqCALCIUM,1))/CQ(12,nq)
            ELSE
              YQ(nq,niqCALCIUM,1)=0.d0
            ENDIF

          ELSE IF(ITYP3(nr,nx).EQ.4) THEN !Beeler-Reuter
              !BR_RATES(Vm,m,h,j,d,f,x1,Cai,...)
C            CALL ASSERT(.FALSE.,'>>You must use march8',ERROR,*9999)
C            CALL BR_RATES(V,YQ(nq,niqMM,1),YQ(nq,niqH,1),
C     '        YQ(nq,niqJ,1),YQ(nq,niqD,1),YQ(nq,niqF,1),
C     '        YQ(nq,niqX,1),YQ(nq,niqCAI,1),CQ(9,nq))

          ELSE IF(ITYP3(nr,nx).EQ.5) THEN !Unused

          ELSE IF(ITYP3(nr,nx).EQ.6) THEN !Luo-Rudy
            CALL LR_RATES(V,YQ(nq,niqMM,1),YQ(nq,niqH,1),
     '        YQ(nq,niqJ,1),YQ(nq,niqD,1),YQ(nq,niqF,1),
     '        YQ(nq,niqX,1),YQ(nq,niqCAI,1),CQ(9,nq))

          ELSE IF(ITYP3(nr,nx).EQ.7) THEN !diFrancesco-Noble

C            CALL DN_RATES(eventually...)

          ENDIF !ityp3

! Integrate the differential equations for the gating variables
C          CALL INTEGRATE_GATE(eventually...)

! Compute ionic current
          IION(nq)=0.0d0
          PHIM=YQ(nq,niqV,1)
          IF(ITYP3(nr,nx).EQ.1) THEN !Cubic w/o recovery
            IION(nq)=CUBIC_ION(PHIM,CQ(9,nq))
          ELSE IF(ITYP3(nr,nx).EQ.2) THEN !FHN
            IION(nq)=FHN_ION(PHIM,YQ(nq,niqRECOV,1),CQ(9,nq))
          ELSE IF(ITYP3(nr,nx).EQ.3) THEN !VCD
            IION(nq)=VCD_ION(PHIM,YQ(nq,niqRECOV,1),YQ(nq,niqDV,1))
          ELSE IF(ITYP3(nr,nx).EQ.4) THEN !BR
            !BR_ION(Vm,m,h,j,d,f,x1,Cai,...)
C            IION(nq)=BR_ION(PHIM,YQ(nq,niqMM,1),YQ(nq,niqH,1),
C     '        YQ(nq,niqJ,1),YQ(nq,niqD,1),YQ(nq,niqF,1),
C     '        YQ(nq,niqX,1),YQ(nq,niqCAI,1),CQ(9,nq))
          ELSE IF(ITYP3(nr,nx).EQ.5) THEN !unused
          ELSE IF(ITYP3(nr,nx).EQ.6) THEN !Luo-Rudy
            IION(nq)=LR_ION(PHIM,YQ(nq,niqMM,1),YQ(nq,niqH,1),
     '        YQ(nq,niqJ,1),YQ(nq,niqD,1),YQ(nq,niqF,1),
     '        YQ(nq,niqX,1),YQ(nq,niqCAI,1),CQ(9,nq))
          ELSE IF(ITYP3(nr,nx).EQ.7) THEN !diFrancesco-Noble
C            IION(nq)=DN_ION(PHIM,YQ(nq,3,1),YQ(nq,4,1),CQ(9,nq))
          ENDIF

        ENDDO !nq
!end of parallel region
      ENDDO !nr

      CALL EXITS('IONIC_CURRENT')
      RETURN
 9999 CALL ERRORS('IONIC_CURRENT',ERROR)
      CALL EXITS('IONIC_CURRENT')
      RETURN 1
      END


      SUBROUTINE TFRONT(IBT,INP,ITHRES,NBJ,NEELEM,NKE,NPF,NPNE,
     '  NRE,NSE1,NSE2,NSE3,NVJE,NXI,CE,CG,CP,PG,SE,T,THRES,XA,
     '  XE,XG,XIG,XP,YG,ERROR,*)

C#### Subroutine: TFRONT
C###  Description:
C###    TFRONT calculates activation pattern by threshold modelling
C###    with fixed time step.

C**** ITYP3(nr,nx)=1 constant wavespeed threshold model
C**** NGP3(i,j,k),i,j,k=1,3 are the Gauss point numbers ng for 3*3*3
C**** NGP5(i,j,k),i,j,k=1,3 are the Gauss point numbers ng for 5*5*5
C**** MNE(i,j,k),i,j,k=-1,1 are the element numbers surrounding ne.

C**** YG(ng,1,ne) is absolute time of Gauss pt ng becoming active
C**** NSE1 is # of Gauss points along a search in Xi1,2 dir
C**** NSE3 is # of Gauss points along a search in Xi3 dir
C**** KTYP31 = 1 for forward tracking of model
C****   "    = 2 for backward   "     "    "

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b12.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:ktyp30.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),INP(NNM,NIM,NBFM),ITHRES(3,NGM,NEM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NKE(NKM,NNM,NBFM,NEM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NRE(NEM),NSE1,NSE2,NSE3,NVJE(NNM,NBFM,NJM,NEM),
     '  NXI(-NIM:NIM,0:4,0:NEM)
      REAL*8 CE(NMM,NEM),CG1TEMP,CG2TEMP,CG3TEMP,CG(NMM,NGM),
     '  CP(NMM,NPM),PG(NSM,NUM,NGM,NBM),SE(NSM,NBFM,NEM),T,
     '  THRES(3,NGM,NEM),XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XG(NJM,NUM),XIG(NIM,NGM,NBM),XP(NKM,NVM,NJM,NPM),YG(NGM,NJM,NEM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,i1,i2,i3,II1,ISUM,j,j1,j2,j3,
     '  JJ1,k,k1,k2,k3,KK1,ME,MG,MNE(-1:1,-1:1,-1:1),
     '  nb,ne,ng,NGIMAX,NGJMAX,NGKMAX,NGP3(3,3,3),NGP5(5,5,5),
     '  NGTB,NITB,
     '  nj,nj1,NK1,NMAX,noelem,nr,nx
      REAL*8 A_VECTOR(3),B_VECTOR(3),G_VECTOR(3),C_VECTOR(3),
     '  CG1,CG2,CG3,COS_ALFA1,COS_ALFA2,COS_ALFA3,
     '  DX,DX1,DX2,DX3,DXI1,DXI2,DXI3,R,SCALAR,SQ,THREST,TIME,
     '  X(3),Z_ng(3),Z_mg(3)

      DATA NGP3/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,
     '  19,20,21,22,23,24,25,26,27/
      DATA NGP5/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,
     '  19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
     '  37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,
     '  55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,
     '  73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,
     '  91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,
     '  107,108,109,110,111,112,113,114,115,116,117,118,119,120,
     '  121,122,123,124,125/

      CALL ENTERS('TFRONT',*9999)
      CALL NENXI(IBT,INP,NBJ,NEELEM,NPNE,NXI,ERROR,*9999) 
      !to return elements surrounding ne in NXI

      nr=1 !Temporary
      nx=1 !Temporary

C *** Calc no of Gauss pts in each direction
      nb=NBJ(1,1)
      NITB=NIT(nb)                             
      NGTB=NGT(nb)
      IF(NGTB.EQ.3**NITB) THEN
        NMAX=3
      ELSE IF(NGTB.EQ.5**NITB) THEN
        NMAX=5
      ELSE
        ERROR='>>Incorrect no of Gauss points'
        GO TO 9999
      ENDIF
      NGIMAX=NMAX               !is number of Gauss points in Xi(1) dir.n
      NGJMAX=NMAX               !is number of Gauss points in Xi(2) dir.n
      IF(NITB.EQ.2) THEN
        NGKMAX=1                !is number of Gauss points in Xi(3) dir.n
      ELSE IF(NITB.EQ.3) THEN
        NGKMAX=NMAX
      ENDIF
                               
      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
C MLB 18 August 1997
C NQGE is now gone, replace with 
C NQE is now gone, CS 15/9/97
C        CALL XPXE(NBJ(1,ne),NKE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
C     '    NQGE(1,ne),NRE(ne),NVJE(1,1,1,ne),
C     '    SE(1,1,ne),XA,XE,XP,ERROR,*9999)

        CALL XPXE(NBJ(1,ne),NKE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '    NRE(ne),NVJE(1,1,1,ne),
     '    SE(1,1,ne),XA,XE,XP,ERROR,*9999)

        CALL CPCG(1,NBJ(1,ne),NPNE(1,1,ne),NRE(ne),nx,CE(1,ne),CG,CP,PG,
     '    ERROR,*9999)
        nb=NBJ(1,ne)
        NITB=NIT(nb)
        MNE(-1,-1,-1)=NXI(-1,1,NXI(-2,1,NXI(-3,1,ne)))
        MNE( 0,-1,-1)=NXI( 0,1,NXI(-2,1,NXI(-3,1,ne)))
        MNE( 1,-1,-1)=NXI( 1,1,NXI(-2,1,NXI(-3,1,ne)))
        MNE(-1, 0,-1)=NXI(-1,1,NXI( 0,1,NXI(-3,1,ne)))
        MNE( 0, 0,-1)=NXI( 0,1,NXI( 0,1,NXI(-3,1,ne)))
        MNE( 1, 0,-1)=NXI( 1,1,NXI( 0,1,NXI(-3,1,ne)))
        MNE(-1, 1,-1)=NXI(-1,1,NXI( 2,1,NXI(-3,1,ne)))
        MNE( 0, 1,-1)=NXI( 0,1,NXI( 2,1,NXI(-3,1,ne)))
        MNE( 1, 1,-1)=NXI( 1,1,NXI( 2,1,NXI(-3,1,ne)))
        MNE(-1,-1, 0)=NXI(-1,1,NXI(-2,1,NXI( 0,1,ne)))
        MNE( 0,-1, 0)=NXI( 0,1,NXI(-2,1,NXI( 0,1,ne)))
        MNE( 1,-1, 0)=NXI( 1,1,NXI(-2,1,NXI( 0,1,ne)))
        MNE(-1, 0, 0)=NXI(-1,1,NXI( 0,1,NXI( 0,1,ne)))
        MNE( 0, 0, 0)=NXI( 0,1,NXI( 0,1,NXI( 0,1,ne)))
        MNE( 1, 0, 0)=NXI( 1,1,NXI( 0,1,NXI( 0,1,ne)))
        MNE(-1, 1, 0)=NXI(-1,1,NXI( 2,1,NXI( 0,1,ne)))
        MNE( 0, 1, 0)=NXI( 0,1,NXI( 2,1,NXI( 0,1,ne)))
        MNE( 1, 1, 0)=NXI( 1,1,NXI( 2,1,NXI( 0,1,ne)))
        MNE(-1,-1, 1)=NXI(-1,1,NXI(-2,1,NXI( 3,1,ne)))
        MNE( 0,-1, 1)=NXI( 0,1,NXI(-2,1,NXI( 3,1,ne)))
        MNE( 1,-1, 1)=NXI( 1,1,NXI(-2,1,NXI( 3,1,ne)))
        MNE(-1, 0, 1)=NXI(-1,1,NXI( 0,1,NXI( 3,1,ne)))
        MNE( 0, 0, 1)=NXI( 0,1,NXI( 0,1,NXI( 3,1,ne)))
        MNE( 1, 0, 1)=NXI( 1,1,NXI( 0,1,NXI( 3,1,ne)))
        MNE(-1, 1, 1)=NXI(-1,1,NXI( 2,1,NXI( 3,1,ne)))
        MNE( 0, 1, 1)=NXI( 0,1,NXI( 2,1,NXI( 3,1,ne)))
        MNE( 1, 1, 1)=NXI( 1,1,NXI( 2,1,NXI( 3,1,ne)))
        IF(DOP) THEN
C$        call mp_setlock()
          WRITE(OP_STRING,'('' Element '',I3,'' MNE:'',27I4)')
     '      ne,(((MNE(I,J,K),I=-1,1),J=-1,1),K=-1,1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C$        call mp_unsetlock()
        ENDIF

        DO 400 k=1,NGKMAX
        DO 400 j=1,NGJMAX
        DO 400 i=1,NGIMAX
          IF(NMAX.EQ.3) THEN
            ng=NGP3(i,j,k)
          ELSE IF(NMAX.EQ.5) THEN
            ng=NGP5(i,j,k)
          ENDIF
          IF(DOP) THEN
C$          call mp_setlock()
            WRITE(OP_STRING,
     '        '('' ne='',I3,'' ng='',I4,'' I='',I2,'' J='',I2,'
     '        //''' K='',I2)') ne,ng,i,j,k
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C$          call mp_unsetlock()
          ENDIF
C ***     Have Primary Gauss point ng.
          IF(ITHRES(1,ng,ne).EQ.1) THEN !ng is in active front
! Calculate spatial coordinates at ng 
            CALL XEXG(NBJ(1,ne),ng,NRE(ne),PG,XE,XG,ERROR,*9999)
c PJH 28-04-92 CALL XGMG(0,NIT(1),1,NRE(ne),DXIX,GL,GU,RG,XG,
c     '         ERROR,*9999)
            DO nj=1,NJT
              X(nj)=XG(nj,1)
            ENDDO
            CALL XZ(ITYP10(1),X,Z_ng)
! Calculate material axis vectors at ng 
            CALL MAT_VEC_NG(NITB,NRE(ne),A_VECTOR,B_VECTOR,C_VECTOR,XG,
     '        ERROR,*9999)
            TIME=THRES(1,ng,ne)+DABS(TINCR) !is time from ng to ellipse
C PJH old   ETA=XG(NJ_LOC(NJL_FIBR,1,nr),1) !is fibre angle to Xi(1) coordinate
C PJH old   COSETA=DCOS(ETA)
C PJH old   SINETA=DSIN(ETA)
            CG1=CG(1,ng)*TIME   !is now dist along 1st (fibre) axis
            CG2=CG(2,ng)*TIME   !is now dist along 2nd (sheet) axis
            IF(NJT.EQ.3) THEN   !3D
              CG3=CG(3,ng)*TIME !is now dist along 3rd (transverse) axis
            ENDIF
            IF(ITHRES(2,ng,ne).EQ.1) THEN !speed up for Purkinje tissue
              CG1=CG1*CG(4,ng)
              CG2=CG2*CG(4,ng)
              IF(NJT.EQ.3) CG3=CG3*CG(4,ng) !3D
            ENDIF   
            IF(DOP) THEN
C$            call mp_setlock()
              CG1TEMP=CG1/TIME
              CG2TEMP=CG2/TIME
              CG3TEMP=CG3/TIME
              WRITE(OP_STRING,
     '        '('' C1 ='',E10.3,'' C2='',E10.3,'' C3='',E10.3)')
     '        CG1TEMP,CG2TEMP,CG3TEMP
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C$            call mp_unsetlock()
            ENDIF
            IF(DOP) THEN
C$            call mp_setlock()
              WRITE(OP_STRING,
     '        '(''  **Active: Time since activation='',E10.3)') TIME
C GBS old    '        //''' Eta='',E10.3)') TIME,ETA
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C$            call mp_unsetlock()
            ENDIF
            DO 300 k1=-NSE3,NSE3
              KK1=k+k1
              k2=MOD(KK1,NMAX)
              IF(k2.LE.0) k2=k2+NMAX
              IF(KK1.LT.1) THEN
                k3=-1
              ELSE IF(KK1.GT.NMAX) THEN
                k3=1     
              ELSE
                k3=0
              ENDIF
            DO 300 j1=-NSE2,NSE2
              JJ1=j+j1
              j2=MOD(JJ1,NMAX)
              IF(j2.LE.0) j2 = j2+NMAX
              IF(JJ1.LT.1) THEN
                j3=-1
              ELSE IF(JJ1.GT.NMAX) THEN
                j3=1        
              ELSE
                j3=0
              ENDIF
            DO 300 i1=-NSE1,NSE1
              II1=i+i1
              i2=MOD(II1,NMAX)
              IF(i2.LE.0) i2 = i2+NMAX
              IF(II1.LT.1)THEN
                i3=-1
              ELSE IF(II1.GT.NMAX) THEN
                i3=1
              ELSE
                i3=0
              ENDIF
              IF(NMAX.EQ.3) THEN
                MG=NGP3(i2,j2,k2)
              ELSE IF(NMAX.EQ.5) THEN
                MG=NGP5(i2,j2,k2)
              ENDIF
              ME=MNE(i3,j3,k3)
C ***         Have secondary Gauss point mg.
              IF(DOP) THEN                 
C$              call mp_setlock()
                WRITE(OP_STRING,'(''  # me='',I3,'' mg='',I3,'
     '          //'''   I1='',I2,'
     '          //''' J1='',I2,'' K1='',I2,'' I2='',I2,'' J2='',I2,'
     '          //''' K2='',I2)') ME,MG,i1,j1,k1,i2,j2,k2
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C$              call mp_unsetlock()
              ENDIF
              IF(ME.GT.0) THEN
                IF(ITHRES(1,MG,ME).EQ.0.
     '            AND.(.NOT.(mg.eq.ng.and.me.eq.ne))) THEN
                  IF(i3.EQ.0) THEN
                    DXI1=(XIG(1,MG,nb)-XIG(1,ng,nb))
                  ELSE
                    IF(XIG(1,MG,nb).GT.XIG(1,ng,nb)) THEN
                      DXI1=-1.0d0+XIG(1,MG,nb)-XIG(1,ng,nb)
                    ELSE
                      DXI1=1.0d0-XIG(1,ng,nb)+XIG(1,MG,nb)
                    ENDIF
                  ENDIF
                  IF(j3.EQ.0) THEN
                    DXI2=(XIG(2,MG,nb)-XIG(2,ng,nb))
                  ELSE
                    IF(XIG(2,MG,nb).GT.XIG(2,ng,nb)) THEN
                      DXI2=-1.0d0+XIG(2,MG,nb)-XIG(2,ng,nb)
                    ELSE
                      DXI2=1.0d0-XIG(2,ng,nb)+XIG(2,MG,nb)
                    ENDIF
                  ENDIF
                  IF(NITB.EQ.2) THEN
                    DXI3=0.0d0
                  ELSE IF(NITB.EQ.3) THEN
                    IF(K3.EQ.0) THEN
                      DXI3=(XIG(3,MG,nb)-XIG(3,ng,nb))
                    ELSE
                      IF(XIG(3,MG,nb).GT.XIG(3,ng,nb)) THEN
                        DXI3=-1.0d0+XIG(3,MG,nb)-XIG(3,ng,nb)
                      ELSE
                        DXI3=1.0d0-XIG(3,ng,nb)+XIG(3,MG,nb)
                      ENDIF       
                    ENDIF
                  ENDIF

! XGRC contains derivs of rect. cart. coords wrt Xi1,2,3
c PJH 28-04-92    CALL XZ_DERIV(ITYP10(1),2,XG,XGRC)
c PJH 28-04-92    CALL XZ_DERIV(ITYP10(1),4,XG,XGRC)
c PJH 28-04-92    CALL XZ_DERIV(ITYP10(1),7,XG,XGRC)
! Calculate spatial coordinates at mg
 
C MLB 18 August 1997
C NQGE is now gone, replace in call with NQE
C NQE is now gone, CS 15/9/97
C                  CALL XPXE(NBJ(1,me),NKE(1,1,1,me),NPF(1,1),
C     '              NPNE(1,1,me),NQGE(1,me),NRE(ne),NVJE(1,1,1,ne),
C     '              SE(1,1,me),XA,XE,XP,
C     '              ERROR,*9999)
                  CALL XPXE(NBJ(1,me),NKE(1,1,1,me),NPF(1,1),
     '              NPNE(1,1,me),NRE(ne),NVJE(1,1,1,ne),
     '              SE(1,1,me),XA,XE,XP,
     '              ERROR,*9999)
                  CALL XEXG(NBJ(1,me),mg,NRE(ne),PG,XE,XG,ERROR,*9999)
                  DO nj=1,NJT
                    X(nj)=XG(nj,1)
                  ENDDO
                  CALL XZ(ITYP10(1),X,Z_mg)
! DX1,2,3 are increments along x1,2,3 (rect. cart.) 
                  DX1=Z_mg(1)-Z_ng(1)      
                  DX2=Z_mg(2)-Z_ng(2)
                  DX3=Z_mg(3)-Z_ng(3)
                  IF(DOP) THEN
C$                  call mp_setlock()
                    WRITE(OP_STRING,'(''   DXI1='',E10.3,'' DXI2='','
     '                //'E10.3,'' DXI3='',E10.3,'' DX1='','
     '                //'E10.3,'' DX2='',E10.3,'' DX3='',E10.3)')
     '                DXI1,DXI2,DXI3,DX1,DX2,DX3
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C$                  call mp_unsetlock()
                  ENDIF
                  DX=DSQRT(DX1*DX1+DX2*DX2+DX3*DX3) !dist betw ng & mg
                  IF(DABS(DX).gt.1.d-8) THEN
! Calculate vector from ng to mg in rect. cart. coords
                    G_VECTOR(1)=DX1
                    G_VECTOR(2)=DX2
                    G_VECTOR(3)=DX3
                    DO nj=1,NJT
                      G_VECTOR(nj)=G_VECTOR(nj)/DX
                    ENDDO
! Calc. components of g_vector wrt material coords             
                    COS_ALFA1=SCALAR(NJT,G_VECTOR,A_VECTOR)
                    COS_ALFA2=SCALAR(NJT,G_VECTOR,B_VECTOR)
                    COS_ALFA3=SCALAR(NJT,G_VECTOR,C_VECTOR)
                    IF(DOP) THEN
C$                    call mp_setlock()
                      WRITE(OP_STRING,'('' cos_ALFA1='',E12.3,'
     '                  //''' cos_ALFA2='',E12.3,'
     '                  //''' cos_ALFA3='',E12.3)') 
     '                  COS_ALFA1,COS_ALFA2,COS_ALFA3
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C$                    call mp_unsetlock()
                    ENDIF
                    IF(NJT.EQ.2) THEN      !2D
                      SQ=DSQRT((COS_ALFA1/CG1)**2+
     '                  (COS_ALFA2/CG2)**2)
                    ELSE IF(NJT.EQ.3) THEN !3D
                      SQ=DSQRT((COS_ALFA1/CG1)**2
     '                  +(COS_ALFA2/CG2)**2+(COS_ALFA3/CG3)**2)
                    ENDIF                  
                    IF(sq.gt.1.d-6) THEN
                      R=1.0d0/SQ !is distance thru mg to wavefront
                    ELSE
                      WRITE(OP_STRING,*) ' SQ is zero'
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    ENDIF
                    IF(R.GE.DX) THEN !wavefront from ng is beyond mg
                      THREST=(1.0d0-DX/R)*(THRES(1,ng,ne)+DABS(TINCR))
                      IF(DOP) THEN
C$                      call mp_setlock()
                        WRITE(OP_STRING,'(''   DX='',E10.3,'' R='','
     '                    //'E10.3,'' THREST='',E10.3,'
     '                    //''' ng= '',I4,'' ne= '',I3)')
     '                  DX,R,THREST,ng,ne
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C$                      call mp_unsetlock()
                      ENDIF
                      IF(THREST.GT.THRES(1,MG,me)) THEN !mg activated from ng
                        THRES(1,MG,me)=THREST !is time since activation of mg
C                       YG(MG,1,me)=YG(ng,1,ne)+DABS(TINCR)-THREST
                        YG(MG,1,me)=T+DABS(TINCR)-THREST
                        IF(DOP) THEN
C$                        call mp_setlock()
                          WRITE(OP_STRING,
     '                      '(''   **Update: THRES(1,mg,me)='','
     '                      //'E10.3,'' TA= '',E10.3,'' ng= '',I4,'
     '                      //''' ne= '',I3)') 
     '                      THRES(1,MG,me),(T+TINCR-THRES(1,MG,me)),
     '                      ng,ne
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C$                        call mp_unsetlock()
                        ENDIF
                      ENDIF  
                    ELSE
                      IF(DOP) THEN
C$                      call mp_setlock()
                        WRITE(OP_STRING,
     '                    '('' R.LT.DX: not activated this step'')')
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C$                      call mp_unsetlock()
                      ENDIF
                    ENDIF
                  ELSE
                    WRITE(OP_STRING,*) ' DX is zero'
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDIF
              ENDIF
 300        CONTINUE
          ENDIF
 400    CONTINUE
      ENDDO

C *** If ITHRES(1,ng,ne)>0 update THRES(1,ng,ne) and
C *** set ITHRES(1,ng,ne)=1 if time since activation,THRES(1,ng,ne),
C *** is >0
      IF(DOP.AND.KTYP31.EQ.1) THEN
        WRITE(OP_STRING,'(''  Points activated this timestep T='','
     '    //'E10.3,'' TINCR was: '',E10.3)') 
     '    (T+TINCR),TINCR
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF 
      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        nb=NBJ(1,ne)
        DO ng=1,NGT(nb)
C         IF(THRES(1,ng,ne).LT.0.) THEN 
!          !Gauss pt was set to nonzero initial
C           THRES(1,ng,ne)=THRES(1,ng,ne)+DABS(TINCR)
C         ENDIF  
          IF(ITHRES(1,ng,ne).GT.0) THEN      !ng is active
            THRES(1,ng,ne)=THRES(1,ng,ne)+DABS(TINCR)
          ELSE IF(ITHRES(1,ng,ne).EQ.0) THEN !ng currently not active
            IF(T+TINCR.GE.YG(ng,1,ne)) THEN
              ITHRES(1,ng,ne)=1
              IF(DOP.AND.KTYP31.EQ.1) THEN
                WRITE(OP_STRING,'('' ng='',I4,'' ne='',I3,'
     '            //'''   Time since activated ='','
     '            //'E10.3,''  TA='',E10.3)') ng,ne,THRES(1,ng,ne),
     '            (T+TINCR-THRES(1,ng,ne))
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ENDDO

C *** Set ITHRES(1,ng,ne)=2 if all points surrounding ng are activated
      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        nb=NBJ(1,ne)
        NITB=NIT(nb)
        MNE(-1,-1,-1)=NXI(-1,1,NXI(-2,1,NXI(-3,1,ne)))
        MNE( 0,-1,-1)=NXI( 0,1,NXI(-2,1,NXI(-3,1,ne)))
        MNE( 1,-1,-1)=NXI( 1,1,NXI(-2,1,NXI(-3,1,ne)))
        MNE(-1, 0,-1)=NXI(-1,1,NXI( 0,1,NXI(-3,1,ne)))
        MNE( 0, 0,-1)=NXI( 0,1,NXI( 0,1,NXI(-3,1,ne)))
        MNE( 1, 0,-1)=NXI( 1,1,NXI( 0,1,NXI(-3,1,ne)))
        MNE(-1, 1,-1)=NXI(-1,1,NXI( 2,1,NXI(-3,1,ne)))
        MNE( 0, 1,-1)=NXI( 0,1,NXI( 2,1,NXI(-3,1,ne)))
        MNE( 1, 1,-1)=NXI( 1,1,NXI( 2,1,NXI(-3,1,ne)))
        MNE(-1,-1, 0)=NXI(-1,1,NXI(-2,1,NXI( 0,1,ne)))
        MNE( 0,-1, 0)=NXI( 0,1,NXI(-2,1,NXI( 0,1,ne)))
        MNE( 1,-1, 0)=NXI( 1,1,NXI(-2,1,NXI( 0,1,ne)))
        MNE(-1, 0, 0)=NXI(-1,1,NXI( 0,1,NXI( 0,1,ne)))
        MNE( 0, 0, 0)=NXI( 0,1,NXI( 0,1,NXI( 0,1,ne)))
        MNE( 1, 0, 0)=NXI( 1,1,NXI( 0,1,NXI( 0,1,ne)))
        MNE(-1, 1, 0)=NXI(-1,1,NXI( 2,1,NXI( 0,1,ne)))
        MNE( 0, 1, 0)=NXI( 0,1,NXI( 2,1,NXI( 0,1,ne)))
        MNE( 1, 1, 0)=NXI( 1,1,NXI( 2,1,NXI( 0,1,ne)))
        MNE(-1,-1, 1)=NXI(-1,1,NXI(-2,1,NXI( 3,1,ne)))
        MNE( 0,-1, 1)=NXI( 0,1,NXI(-2,1,NXI( 3,1,ne)))
        MNE( 1,-1, 1)=NXI( 1,1,NXI(-2,1,NXI( 3,1,ne)))
        MNE(-1, 0, 1)=NXI(-1,1,NXI( 0,1,NXI( 3,1,ne)))
        MNE( 0, 0, 1)=NXI( 0,1,NXI( 0,1,NXI( 3,1,ne)))
        MNE( 1, 0, 1)=NXI( 1,1,NXI( 0,1,NXI( 3,1,ne)))
        MNE(-1, 1, 1)=NXI(-1,1,NXI( 2,1,NXI( 3,1,ne)))
        MNE( 0, 1, 1)=NXI( 0,1,NXI( 2,1,NXI( 3,1,ne)))
        MNE( 1, 1, 1)=NXI( 1,1,NXI( 2,1,NXI( 3,1,ne)))
        DO 650 k=1,NGKMAX
        DO 650 j=1,NGJMAX
        DO 650 i=1,NGIMAX
          IF(NMAX.EQ.3)THEN
            ng=NGP3(i,j,k)
          ELSE
            ng=NGP5(i,j,k)
          ENDIF
          ISUM=0
          nj1=1
          IF(NITB.LE.2)THEN
            NK1=0
          ELSE
            NK1=1
          ENDIF
          DO k1=-NK1,NK1
            KK1=k+k1
            k2=MOD(KK1,NMAX)
            IF(k2.LE.0) k2=k2+NMAX
            IF(Kk1.LT.1) THEN
            k3=-1
          ELSE IF(KK1.GT.NMAX) THEN
            k3=1
          ELSE
            k3=0
          ENDIF
            DO j1=-nj1,nj1
              JJ1=j+j1
              j2=MOD(JJ1,NMAX)
              IF(j2.LE.0) j2 = j2+NMAX
              IF(JJ1.LT.1)THEN
                j3=-1
              ELSE IF(JJ1.GT.NMAX)THEN
                j3=1
              ELSE
                j3=0
              ENDIF
              DO i1=-nj1,nj1
                II1=i+i1
                i2=MOD(II1,NMAX)
                IF(i2.LE.0)i2 = i2+NMAX
                IF(II1.LT.1)THEN
                  i3=-1
                ELSE IF(II1.GT.NMAX)THEN
                  i3=1
                ELSE
                  i3=0
                ENDIF
                IF(NMAX.EQ.3)THEN
                  MG=NGP3(i2,j2,k2)
                ELSE
                  MG=NGP5(i2,j2,k2)
                ENDIF
                me=MNE(i3,j3,k3)
C ***           Have primary and secondary point.
                IF(me.GT.0.AND.(.NOT.(mg.eq.ng.and.me.eq.ne))) THEN
                  IF(ITHRES(1,MG,me).GE.1) ISUM=ISUM+1
                ELSE
                  IF(me.EQ.0) ISUM=ISUM+1 !for boundaries
                ENDIF
              ENDDO
            ENDDO
          ENDDO
! Check # of active Gauss points surrounding an active Gauss point
! and set ithres(1) to 2 if surrounded by active points
          IF(NITB.EQ.2) THEN
            IF(ITHRES(1,ng,ne).ne.0) THEN
              IF(ISUM.EQ.8) ITHRES(1,ng,ne)=2
            ENDIF
          ELSE IF(NITB.EQ.3) THEN
            IF(ITHRES(1,ng,ne).ne.0) THEN
              IF(ISUM.EQ.26) ITHRES(1,ng,ne)=2
            ENDIF
          ENDIF
 650    CONTINUE
      ENDDO

      CALL EXITS('TFRONT')
      RETURN
 9999 CALL ERRORS('TFRONT',ERROR)
      CALL EXITS('TFRONT')
      RETURN 1
      END
