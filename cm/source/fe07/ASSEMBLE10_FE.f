      SUBROUTINE ASSEMBLE10_FE(ISC_GKK,ISR_GKK,NEELEM,NLATNE,
     &  NQGP,NQGP_PIVOT,NQNLAT,NQS,NQSCNB,NQXI,NRLIST,nx_ext,nx_trans,
     &  nx_upd,NXQ,CQ,GKK,GM,NQGW,PG,WG,XQ,PROPQ,BIDOMAIN,COUPBID,
     &  FIRST_A,FIXQ,IMPLICIT,UPDATE,UPDATEDT,UPDATE_MATRIX,
     &  SOLVEEIGHTPROBLEM,ERROR,*)

C#### Subroutine: ASSEMBLE10_FE
C###  Description:
C###    ASSEMBLE10_FE writes directly into the compressed row storage
C###    arrays. It generates the constraint matrix for the implicit
C###    solution of activation problems, using the Grid-based Finite
C###    Element method.
C***  Created by Scott Marsden November 2000
C***  GBS August 2003 - archived all code for obsolete sparse solvers
C***  GBS Sept 2003   - added lattice-based elements to Grid-FEM

      IMPLICIT NONE

      INCLUDE 'b01.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'solv00.cmn'
      INCLUDE 'time02.cmn'
      INCLUDE 'tol00.cmn'
      INCLUDE 'mach00.inc'
      
!     Parameter list
      INTEGER ISC_GKK(NISC_GKKM,NXM),ISR_GKK(NISR_GKKM,NXM),
     &  NEELEM(0:NE_R_M,0:NRM),NLATNE(NEQM+1),NQGP(0:NQGM,NQM),
     &  NQGP_PIVOT(NQGM,NQM),NQNLAT(NEQM*NQEM),NQS(NEQM),NQSCNB(NQSCM),
     &  NQXI(0:NIM,NQSCM),NRLIST(0:NRM),nx_ext,nx_trans,nx_upd,
     &  NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      REAL*8 CQ(NMM,NQM),GKK(NZ_GKK_M,NXM),GM(NZ_GM_M),NQGW(NQGM,NQM),
     &  PG(NSM,NUM,NGM,NBM),WG(NGM,NBM),XQ(NJM,NQM),PROPQ(3,3,4,2,NQM)
      CHARACTER ERROR*(*)
      LOGICAL BIDOMAIN,COUPBID,FIRST_A,FIXQ(NYQM,NIYFIXM,NXM),IMPLICIT,
     &  SOLVEEIGHTPROBLEM,UPDATE,UPDATEDT,UPDATE_MATRIX

!     Local variables
      INTEGER DUMMY_LIST(0:1),IEND,Inc_NZ_GKK_M,Inc_NISR_GKK_M,
     &  Inc_NISC_GKK_M,maxrow,nb,NITB,nq,nq_dummy,nr,nrr,nx,nzero,nzz,
     &  nzz_trans,nzz_ext,nzz_upd,PLACEINT,PLACEEXT,SCHEME
      INTEGER*4 COEFFSEXT_PTR
      REAL*8 DTTHETA
      REAL ELAPSED_TIME,TIME_START(1),TIME_STOP(1)
      LOGICAL ERROR_FLAG
C SGM 29 November 2000 not used
C      LOGICAL NOTIMP

!     Functions
      REAL*8 GET_COEFFSEXT

      ! initialise pointer now for error handling at 9999
      COEFFSEXT_PTR=0

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
     &          maxrow
              GOTO 9999
            ENDIF
            !initialise
            ISR_GKK(1,nx_trans)=1
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
     &          maxrow
              GOTO 9999
            ENDIF
            !initialise
            ISR_GKK(1,nx_ext)=1
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
     &          maxrow
              GOTO 9999
            ENDIF
            !initialise
            ISR_GKK(1,nx_upd)=1
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
        SCHEME=NQS(NEELEM(1,nr))
        nb=NQSCNB(SCHEME)
        NITB=NQXI(0,SCHEME)
        
C make sure the number of element  Xi coordinates equals the number of
C       Xi coordiates for the dependent variable before call to
C     CALC_FE_GRID_COEF.
        CALL ASSERT(NITB.EQ.NIT(nb),
     &    '>>#Xi-coords inconsistent between element and dependent var',
     &    ERROR,*9999)
        
C SGM 01Feb01 Moved following two checks inside nr loop.
        !Check array sizes
C MLT 30/11/02 Modified to account for different support for 
C grid FE and grid FV
        IF(ITYP4(nr,nx).EQ.6.AND.NQGM.LT.3**NITB) THEN ! Grid FE
          IEND=0
          CALL APPENDC(IEND,'For Grid FE NQGM needs to be at least ',
     &      ERROR)
          CALL APPENDI(IEND,3**NITB,ERROR)
          GOTO 9999 
        ELSE IF(ITYP4(nr,nx).EQ.7.AND.NQGM.LT.2*NITB+1) THEN ! Grid FV
          IEND=0
          CALL APPENDC(IEND,'For Grid FV NQGM needs to be at least ',
     &      ERROR)
          CALL APPENDI(IEND,2*NITB+1,ERROR)
          GOTO 9999 
        ENDIF
        
        IF(NZ_GM_M.LT.NQT*NQGM) THEN
          IEND=0
          CALL APPENDC(IEND,'NZ_GM_M needs to be at least ',ERROR)
          CALL APPENDI(IEND,NQT*NQGM,ERROR)
          GOTO 9999
        ENDIF

        
C LKC 12-AUG-2005 Want to allocate COEFFSEXT because we don't need it
C have a second nq variable for the non-lattice code, and this local
C variable can be too large to be stored on the stack for large problems
C
C NOTE/CAUTION: COEFFSEXT is (NQGM,NQM) for the lattice code but only
C        COEFFSEXT (NGGM) for the non-lattice code
C
        IF(USE_LAT.EQ.1) THEN
          CALL ALLOCATE_MEMORY(NQGM*NQM,1,DPTYPE,COEFFSEXT_PTR,
     '      MEM_INIT,ERROR,*9997)
          CALL CALC_FE_LATT_COEF(nb,NITB,NLATNE,NQGP,NQGP_PIVOT,NQNLAT,
     &      NQS,NQXI,nr,nx_ext,%VAL(COEFFSEXT_PTR),CQ,GM,NQGW,PG,WG,XQ,
     &      BIDOMAIN,FIXQ,SOLVEEIGHTPROBLEM,ERROR,*9999)
        ELSE
          CALL ALLOCATE_MEMORY(NQGM,1,DPTYPE,COEFFSEXT_PTR,
     '      MEM_INIT,ERROR,*9997)
        ENDIF
        
        
        DO nq=NQR(1,nr),NQR(2,nr) !create line for each nq

C MLT 29Nov02 Adding grid finite volumes

C LKC 14-JUN-2004 Removing the second index (nq) on the local array COEFFSEXT
C as it is not used any more (assembly is done grid point by grid point) 
C and the large local sometimes crashes on a linux pc.
C The routine CALC_FE_LATT_COEF has also been modified.

          IF(USE_LAT.EQ.1) THEN
            nq_dummy=nq
          ELSE
            nq_dummy=1
          ENDIF

          
          IF(USE_LAT.EQ.0) THEN
            IF(ITYP4(nr,nx).EQ.6) THEN ! Grid-based FE              
              CALL CALC_FE_GRID_COEF(nb,NITB,nq,NQGP,NQGP_PIVOT,nr,
     &          nx_ext,NXQ(-NIM,0,0,1),%VAL(COEFFSEXT_PTR),CQ,GM,
     &          NQGW(1,nq),PG,WG,XQ,BIDOMAIN,FIXQ,SOLVEEIGHTPROBLEM,
     &          ERROR,*9999)
            ELSE ! Grid FV
              CALL CALC_FV_GRID_COEF(NITB,nq,NQGP,NQGP_PIVOT,nr,nx_ext,
     &          NXQ(-NIM,0,0,1),%VAL(COEFFSEXT_PTR),CQ,GM,NQGW(1,nq),XQ,
     &          PROPQ,BIDOMAIN,FIXQ,SOLVEEIGHTPROBLEM,ERROR,*9999)
            ENDIF            
          ENDIF !USE_LAT

          
C Check array sizes (NZ_GKK_M,ISR_GKK,ISC_GKK)
          IF(.NOT.UPDATE) THEN
            IF(IMPLICIT) THEN
              IF(SPARSEGKK(nx_trans).EQ.1) THEN !compressed row
                nzz_trans = nzz_trans + NQGP(0,nq)
                NZZT(1,nr,nx_trans)=nzz_trans
                IF(NZ_GKK_M.LT.nzz_trans) THEN
                  IF(Inc_NZ_GKK_M.LT.nzz_trans)Inc_NZ_GKK_M = nzz_trans
                ELSE IF(NISC_GKKM.LT.nzz_trans) THEN
                  IF(Inc_NISC_GKK_M.LT.nzz_trans)Inc_NISC_GKK_M =
     &              nzz_trans
                ENDIF
              ENDIF
            ENDIF !implicit
            IF(BIDOMAIN) THEN
              IF(SPARSEGKK(nx_ext).EQ.1) THEN !compressed row
                nzz_ext = nzz_ext + NQGP(0,nq)
                NZZT(1,nr,nx_ext)=nzz_ext
                IF(NZ_GKK_M.LT.nzz_ext) THEN
                  IF(Inc_NZ_GKK_M.LT.nzz_ext)Inc_NZ_GKK_M = nzz_ext
                ELSE IF(NISC_GKKM.LT.nzz_ext) THEN
                  IF(Inc_NISC_GKK_M.LT.nzz_ext)Inc_NISC_GKK_M = nzz_ext
                ENDIF
              ENDIF
            ENDIF !bidomain
            IF(COUPBID) THEN
              IF(SPARSEGKK(nx_upd).EQ.1) THEN !compressed row
                nzz_upd = nzz_upd + NQGP(0,nq)
                NZZT(1,nr,nx_upd)=nzz_upd
                IF(NZ_GKK_M.LT.nzz_upd) THEN
                  IF(Inc_NZ_GKK_M.LT.nzz_upd)Inc_NZ_GKK_M = nzz_upd
                ELSE IF(NISC_GKKM.LT.nzz_upd) THEN
                  IF(Inc_NISC_GKK_M.LT.nzz_upd)Inc_NISC_GKK_M = nzz_upd
                ENDIF
              ENDIF
            ENDIF !coupid
          ENDIF !not update

          IF(Inc_NZ_GKK_M.EQ.0.AND.Inc_NISR_GKK_M.EQ.0
     &      .AND.Inc_NISC_GKK_M.EQ.0) THEN !don't need to increase array sizes
            IF(UPDATEDT) THEN
              IF(IMPLICIT) THEN
                DTTHETA=DT*THETA(1)
                IF(SPARSEGKK(nx_trans).EQ.0) THEN !no sparsity
                  DO nzero=1,NQGP(0,nq)
                    GKK(((NQGP(nzero,nq)-1)*NQT)+nq,nx_trans)=
     &                -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTTHETA +
     &                GM(NQGP_PIVOT(nzero,nq)+NQGM*(nq-1))
                  ENDDO
                ELSE IF(SPARSEGKK(nx_trans).EQ.1) THEN !compressed row
                  PLACEINT=ISR_GKK(nq,nx_trans)
                  DO nzero=1,NQGP(0,nq)
                    GKK(PLACEINT,nx_trans)=
     &                -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTTHETA +
     &                GM(NQGP_PIVOT(nzero,nq)+NQGM*(nq-1))
                    PLACEINT=PLACEINT+1
                  ENDDO !nzero
                ENDIF
                IF(ERROR_FLAG) GOTO 9998
              ENDIF !implicit
            ELSE !not updatedt
              IF(IMPLICIT) THEN
                DTTHETA=DT*THETA(1)
                IF(SPARSEGKK(nx_trans).EQ.0) THEN !no sparsity
                  DO nzero=1,NQGP(0,nq)
                    GKK(((NQGP(nzero,nq)-1)*NQT)+nq,nx_trans)=
     &                -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTTHETA +
     &                GM(NQGP_PIVOT(nzero,nq)+NQGM*(nq-1))
                  ENDDO
                ELSE IF(SPARSEGKK(nx_trans).EQ.1) THEN !compressed row
                  DO nzero=1,NQGP(0,nq)
                    PLACEINT=PLACEINT+1
                    IF(.NOT.UPDATE)
     &                ISC_GKK(PLACEINT,nx_trans)=NQGP(nzero,nq)
                    GKK(PLACEINT,nx_trans)=
     &                -NQGW(NQGP_PIVOT(nzero,nq),nq)*DTTHETA +
     &                GM(NQGP_PIVOT(nzero,nq)+NQGM*(nq-1))
                  ENDDO
                  IF(.NOT.UPDATE) ISR_GKK(nq+1,nx_trans) =
     &              ISR_GKK(nq,nx_trans)+NQGP(0,nq)
                ENDIF !SPARSEGKK
              ENDIF !implicit

              IF(BIDOMAIN) THEN
                IF(SPARSEGKK(nx_ext).EQ.0) THEN !no sparsity
                  DO nzero=1,NQGP(0,nq)
C LKC 12-AUG-2005 Get values from dynamic array                    
C                    GKK(((NQGP(nzero,nq)-1)*NQT)+nq,nx_ext)=
C     &                COEFFSEXT(NQGP_PIVOT(nzero,nq))
                    GKK(((NQGP(nzero,nq)-1)*NQT)+nq,nx_ext)=
     &                GET_COEFFSEXT(%VAL(COEFFSEXT_PTR),
     &                NQGP_PIVOT(nzero,nq),nq_dummy)
                    
                  ENDDO
                ELSE IF(SPARSEGKK(nx_ext).EQ.1) THEN !compressed row
                  DO nzero=1,NQGP(0,nq)
                    PLACEEXT=PLACEEXT+1
                    IF(.NOT.UPDATE)
     &                ISC_GKK(PLACEEXT,nx_ext)=NQGP(nzero,nq)
C LKC 12-AUG-2005 Get values from dynamic array                    
C                    GKK(PLACEEXT,nx_ext)=
C     &                COEFFSEXT(NQGP_PIVOT(nzero,nq))
                    GKK(PLACEEXT,nx_ext)=
     &                GET_COEFFSEXT(%VAL(COEFFSEXT_PTR),
     &                NQGP_PIVOT(nzero,nq),nq_dummy)
                  ENDDO
                  IF(.NOT.UPDATE)
     &              ISR_GKK(nq+1,nx_ext)=ISR_GKK(nq,nx_ext)+NQGP(0,nq)
                ENDIF
              ENDIF !bidomain

              IF(COUPBID) THEN
                IF(SPARSEGKK(nx_upd).EQ.0) THEN !no sparsity
                  DO nzero=1,NQGP(0,nq)
C LKC 12-AUG-2005 Get values from dynamic array
C                    GKK(((NQGP(nzero,nq)-1)*NQT)+nq,nx_upd)=
C     &                COEFFSEXT(NQGP_PIVOT(nzero,nq))
                    GKK(((NQGP(nzero,nq)-1)*NQT)+nq,nx_upd)=
     &                GET_COEFFSEXT(%VAL(COEFFSEXT_PTR),
     &                NQGP_PIVOT(nzero,nq),nq_dummy)
                  ENDDO
                ELSE IF(SPARSEGKK(nx_upd).EQ.1) THEN !compressed row
                  DO nzero=1,NQGP(0,nq)
                    PLACEEXT=PLACEEXT+1
                    IF(.NOT.UPDATE)
     &                ISC_GKK(PLACEEXT,nx_upd)=NQGP(nzero,nq)
C LKC 12-AUG-2005 Get values from dynamic array                    
C                    GKK(PLACEEXT,nx_upd)=
C     &                COEFFSEXT(NQGP_PIVOT(nzero,nq))
                    GKK(PLACEEXT,nx_upd)=
     &                GET_COEFFSEXT(%VAL(COEFFSEXT_PTR),
     &                NQGP_PIVOT(nzero,nq),nq_dummy)
                  ENDDO
                  IF(.NOT.UPDATE)
     &              ISR_GKK(nq+1,nx_upd)=ISR_GKK(nq,nx_upd)+NQGP(0,nq)
                ENDIF
              ENDIF !coupid
            ENDIF !updatedt
          ENDIF !don't need to increase array sizes
        ENDDO !nq
      ENDDO !nr

      IF(Inc_NZ_GKK_M.GT.0) THEN
        WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)')Inc_NZ_GKK_M
        GOTO 9999
      ELSE IF(Inc_NISR_GKK_M.GT.0) THEN
        WRITE(ERROR,'(''>>Increase ISR_GKK_M to >= '',I12)'
     &    )Inc_NISR_GKK_M
        GOTO 9999
      ELSE IF(Inc_NISC_GKK_M.GT.0) THEN
        WRITE(ERROR,'(''>>Increase ISC_GKK_M to >= '',I12)'
     &    )Inc_NISC_GKK_M
        GOTO 9999
      ENDIF

      IF(.NOT.UPDATEDT) THEN
        IF(.NOT.UPDATE) THEN
          IF(IMPLICIT) THEN
            NZZT(1,NRLIST(1),nx_trans)=PLACEINT
          ENDIF
          IF(BIDOMAIN) THEN
            NZZT(1,NRLIST(1),nx_ext)=PLACEEXT
          ENDIF
          IF(COUPBID) THEN
            NZZT(1,NRLIST(1),nx_upd)=PLACEEXT
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

      IF(IWRIT5(NRLIST(1),nx).GE.1) THEN
        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)
        WRITE(OP_STRING,'(1X,''Time for matrix assembly '',F8.2,'
     &    //'''s cpu'')') ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      CALL FREE_MEMORY(COEFFSEXT_PTR,ERROR,*9999)      
      CALL EXITS('ASSEMBLE10_FE')
      RETURN
 9998 ERROR=' '      
 9999 IF(COEFFSEXT_PTR.NE.0) CALL FREE_MEMORY(COEFFSEXT_PTR,ERROR,*9997)
 9997 CALL ERRORS('ASSEMBLE10_FE',ERROR)
      CALL EXITS('ASSEMBLE10_FE')
      RETURN 1
      END



      REAL*8 FUNCTION GET_COEFFSEXT(COEFFSEXT,i,j)


C#### Function: GET_REALVAL
C###  Type: REAL*8
C###  Description:
C###    Returns a real value from a dynamically allocated array (COEFFSEXT)

      IMPLICIT NONE
      
      INCLUDE 'geom00.cmn'
      
      REAL*8 COEFFSEXT(NQGM,*)
      INTEGER i,j
      
      GET_COEFFSEXT=COEFFSEXT(i,j)
      
      RETURN
      END

      
