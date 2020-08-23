      SUBROUTINE ASSEMBLE2(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,IBT,IDO,
     '  INP,ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,ISR_GKK,ISR_GQ,
     '  NBH,NBJ,NDET,NDIPOLES,NEELEM,NENP,NGAP,NHE,NHP,NKH,NKHE,NKJE,
     '  NLL,NONY,NPF,NP_INTERFACE,NPNE,NPNODE,NPNY,nr,nr_gkk,
     '  NRE,NVHP,NVJE,NW,nx,NYNE,NYNP,NYNR,CE,CONY,CURVCORRECT,
     '  DET,DIPOLE_CEN,DIPOLE_DIR,DL,GD,GK,GKK,GQ,PG,
     '  SE,TIME,WG,XA,XE,XG,XIG,XP,
     '  UPDATE_MATRIX,UPDATE_SOURCE,RET_ERROR,*)

C#### Subroutine: ASSEMBLE2
C###  Description:
C###    ASSEMBLE2 assembles the global unreduced matrices GK,
C###    GM, etc. for Static linear BEM problems.

      IMPLICIT NONE
      INCLUDE 'b10.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'error0.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'iwrit00.cmn'
C      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'mxch.inc'
      INCLUDE 'time02.cmn'
!     Parameter List
      INTEGER DIPOLE_CEN_NTIME(NDIPOLEM,NRM,NXM),
     '  DIPOLE_DIR_NTIME(NDIPOLEM,NRM,NXM),IBT(3,NIM,NBFM),
     '  IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISC_GK(NISC_GKM),ISC_GKK(NISC_GKKM),ISC_GQ(NISC_GQM),
     '  ISR_GK(NISR_GKM),ISR_GKK(NISR_GKKM),
     '  ISR_GQ(NISR_GQM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NDET(NBFM,0:NNM),NDIPOLES(NRM,NXM),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),NGAP(NIM,NBM),NHE(NEM),NHP(NPM),
     '  NKH(NHM,NPM,NCM),NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NLL(12,NEM),NONY(0:NOYM,NYM,NRCM),
     '  NP_INTERFACE(0:NPM,0:3),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM),nr,nr_gkk,NRE(NEM),
     '  NVHP(NHM,NPM,NCM),NVJE(NNM,NBFM,NJM,NEM),
     '  NW(NEM,3),nx,NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM)
      REAL*8 CE(NMM,NEM),CONY(0:NOYM,NYM,NRCM),CURVCORRECT(2,2,NNM,NEM),
     '  DET(NBFM,0:NNM,NGM,6),DIPOLE_CEN(4,0:NDIPTIMM,NDIPOLEM,NRM,NXM),
     '  DIPOLE_DIR(4,0:NDIPTIMM,NDIPOLEM,NRM,NXM),DL(3,NLM),
     '  GD(NZ_GD_M),GK(NZ_GK_M),GKK(NZ_GKK_M),
     '  GQ(NZ_GQ_M),PG(NSM,NUM,NGM,NBM),
     '  SE(NSM,NBFM,NEM),
     '  TIME,WG(NGM,NBM),XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),
     '  XIG(NIM,NGM,NBM),
     '  XP(NKM,NVM,NJM,NPM)
      CHARACTER RET_ERROR*(*)
      LOGICAL UPDATE_MATRIX,UPDATE_SOURCE
!     Local Variables
      INTEGER ne,nj,njj2,noelem,nonode,no_nynr1,
     '  NO_TIMES,np,ny1

      INTEGER*4 LGKE_PTR,LGQE_PTR,NPB_PTR,
     '  DET_ADAPT_PTR,DRDN_PTR,DRDNO_PTR,GKES_PTR,GQES_PTR,PG_J_PTR,
     '  PG_Q_PTR,PG_U_PTR,RAD_PTR,RD_PTR,RG_PTR,XA_PTR,XE_PTR,XG_PTR,
     '  XG1_PTR,XIG_J_PTR,XIG_Q_PTR,XIG_U_PTR,XN_PTR,XN_GRAD_PTR,
     '  XR_PTR,XR_GRAD_PTR
      REAL  AVETIME,ELAPSED_TIME,TIME_START1(1),TIME_START2(1),
     '  TIME_START3(1),TIME_STOP(1)
      REAL*8 XPFP(3)
      CHARACTER ERROR*(ERRSTRLEN),ERROR_DUMMY*255
      LOGICAL ERROR_FLAG,FULL_EQUATIONS,INTERFACE

      CALL ENTERS('ASSEMBLE2',*9999)

C CPB 18/12/97 Adding BE stiffness matrices
C cpb 23/1/97 Fixing multiprocessing adaptive integration
C CPB 18/10/96 Intialise pointers so the error free will work properly

      LGKE_PTR=0
      LGQE_PTR=0
      NPB_PTR=0
      DET_ADAPT_PTR=0
      DRDN_PTR=0
      DRDNO_PTR=0
      GKES_PTR=0
      GQES_PTR=0
      PG_J_PTR=0
      PG_Q_PTR=0
      PG_U_PTR=0
      RAD_PTR=0
      RD_PTR=0
      RG_PTR=0
      XA_PTR=0
      XE_PTR=0
      XG_PTR=0
      XG1_PTR=0
      XIG_J_PTR=0
      XIG_Q_PTR=0
      XIG_U_PTR=0
      XN_PTR=0
      XN_GRAD_PTR=0
      XR_PTR=0
      XR_GRAD_PTR=0

C***  Initialise matrices for this region
      CALL CPU_TIMER(CPU_USER,TIME_START1)
      CALL CPU_TIMER(CPU_USER,TIME_START2)

      IF(UPDATE_MATRIX) THEN
        CALL INIT_SPARSE_MATRIX(NYNR(0,2,1),NYNR(0,2,1),ISC_GK,ISR_GK,
     '    0,0,1,NYT(1,1,nx),NZ_GK_M,NZT(1,nx),NYNR(0,1,1),NYNR(0,1,1),
     '    KTYP24,GK,ERROR,*9999)
        CALL INIT_SPARSE_MATRIX(NYNR(0,2,2),NYNR(0,2,2),ISC_GQ,ISR_GQ,
     '    0,0,1,NYT(1,2,nx),NZ_GQ_M,NZT(2,nx),NYNR(0,1,2),NYNR(0,1,2),
     '    KTYP24,GQ,ERROR,*9999)
      ENDIF


C*** LKC 28-NOV-1999 Special source poisson problem
C***    should not go into this loop - USE_DIPOLE=0 or the
C***    number of dipoles will be 0

      IF(UPDATE_SOURCE.AND.(USE_DIPOLE.EQ.1)) THEN
        IF(IGREN(nr).GE.9.OR.NDIPOLES(nr,nx).GT.0) THEN
C         Poisson or/with dipole
          NYT(1,3,nx)=NYT(1,1,nx)
          NYT(2,3,nx)=1
          CALL ASSERT(NYT(1,3,nx).LE.NZ_GD_M,'>>Increase NZ_GD_M',
     '      ERROR,*9999)
C$OMP     PARALLEL DO
C$OMP&      PRIVATE(no_nynr1,ny1),
C$OMP&      SHARED(NYNR,GD)
          DO no_nynr1=1,NYNR(0,1,1)
            ny1=NYNR(no_nynr1,1,1)
            GD(ny1)=0.0d0
          ENDDO !no_nynr1
C$OMP     END PARALLEL DO

C LKC 16-OCT-2000 If there is no dipole in the region, then GD
C  is not being initialised and then it gets assembled into
C  the RHS in solve2/9. (If there is a dipole it gets setup in XPGD)
C
C        ENDIF
C      ENDIF ! Update source

        ELSE
C$OMP PARALLEL DO
C$OMP&  PRIVATE(no_nynr1,ny1),
C$OMP&  SHARED(NYNR,GD)
          DO no_nynr1=1,NYNR(0,1,1)
            ny1=NYNR(no_nynr1,1,1)
            GD(ny1)=0.0d0
          ENDDO !no_nynr1
C$OMP END PARALLEL DO


        ENDIF ! if dipole and poisson
      ENDIF ! Update source

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
      IF(IWRIT4(nr,nx).GE.1) THEN
        WRITE(OP_STRING,'(/'' CPU time for setup and '
     '    //'initialisation: '',D11.4,'' s'')') ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      IF(UPDATE_MATRIX) THEN

C CPB 19/10/98 Checking hard coded dimensions in XEGKGQ_3DL and XPGKGQ_3DL
      IF(OPTI3DLAPLACE) THEN
        CALL ASSERT(NBFM.LE.99,'>>Increase NBFM hard coded '
     '    //'dimension in XEGKGQ_3DL and XPGKGQ_3DL',ERROR,*9999)
        CALL ASSERT(NHM.LE.9,'>>Increase NHM hard coded dimension '
     '    //'in XEGKGQ_3DL and XPGKGQ_3DL',ERROR,*9999)
        CALL ASSERT(NIM.LE.3,'>>Increase NIM hard coded dimension '
     '    //'in XEGKGQ_3DL and XPGKGQ_3DL',ERROR,*9999)
        CALL ASSERT(NJM.LE.9,'>>Increase NJM hard coded dimension '
     '    //'in XEGKGQ_3DL and XPGKGQ_3DL',ERROR,*9999)
        CALL ASSERT(NGM.LE.99,'>>Increase NGM hard coded dimension '
     '    //'in XEGKGQ_3DL and XPGKGQ_3DL',ERROR,*9999)
        CALL ASSERT(NKM.LE.8,'>>Increase NKM hard coded dimension '
     '    //'in XEGKGQ_3DL and XPGKGQ_3DL',ERROR,*9999)
        CALL ASSERT(NNM.LE.64,'>>Increase NNM hard coded dimension '
     '    //'in XEGKGQ_3DL and XPGKGQ_3DL',ERROR,*9999)
        CALL ASSERT(NSM.LE.64,'>>Increase NSM hard coded dimension '
     '    //'in XEGKGQ_3DL and XPGKGQ_3DL',ERROR,*9999)
        CALL ASSERT(NUM.LE.11,'>>Increase NUM hard coded dimension '
     '    //'in XEGKGQ_3DL and XPGKGQ_3DL',ERROR,*9999)
      ENDIF

C***  Find element stiffness matrices

C Check to see if we are in a region which doesn't require the extra
C equations at the interface.
        FULL_EQUATIONS=.FALSE.
        IF(KTYP91.EQ.1)THEN
          FULL_EQUATIONS=.TRUE.
        ENDIF !End of KTYP91 loop

C cpb 28/6/96 Adding parallel directives and node based outer loops
        CALL CPU_TIMER(CPU_USER,TIME_START2)
        AVETIME=0.0
        NO_TIMES=0
        ERROR_FLAG=.FALSE.
        IF(BEMLOOPTYPE.EQ.1) THEN !Element outer loop
C$OMP     PARALLEL DO
C$OMP&    PRIVATE(noelem,ne,np,INTERFACE,TIME_START3,TIME_STOP,
C$OMP&      ELAPSED_TIME,LGKE_PTR,LGQE_PTR,NPB_PTR,DET_ADAPT_PTR,
C$OMP&      DRDN_PTR,DRDNO_PTR,GKES_PTR,GQES_PTR,PG_J_PTR,
C$OMP&      PG_Q_PTR,PG_U_PTR,RAD_PTR,RD_PTR,RG_PTR,XA_PTR,XE_PTR,
C$OMP&      XG_PTR,XG1_PTR,XIG_J_PTR,XIG_Q_PTR,XIG_U_PTR,XN_PTR,
C$OMP&      XN_GRAD_PTR,XR_PTR,XR_GRAD_PTR,ERROR),
C$OMP&    SHARED(FULL_EQUATIONS,IOOP,MEM_INIT,NEELEM,NPNE,
C$OMP&      NP_INTERFACE,NGM,nr,nr_gkk,nx,ERROR_FLAG),
C$OMP&    REDUCTION(+:AVETIME,NO_TIMES)
          DO noelem=1,NEELEM(0,nr)
            CALL CPU_TIMER(CPU_USER,TIME_START3)
            ne=NEELEM(noelem,nr)
            IF(.NOT.ERROR_FLAG) THEN

C             INTERFACE is true if element ne is on the interface
C             between regions in a coupled problem and the region to
C             which ne belongs is NOT the region with the smallest
C             region number.
              np=NPNE(1,NBJ(1,ne),ne)
              INTERFACE=(NP_INTERFACE(np,0).GT.1).AND.
     '          (NP_INTERFACE(np,1).NE.nr)

              IF(OPTI3DLAPLACE) THEN

C cpb 23/9/98 Where possible don't dynamically allocate local
C parallel arrays but hard code the dimensions in XEGKGQ_3DL

                NPB_PTR=0
                CALL ALLOCATE_MEMORY((NP_R_M+1)*5,1,INTTYPE,NPB_PTR,
     '            MEM_INIT,ERROR,*100)

                CALL XEGKGQ_3DL(IBT,IDO,INP,ISC_GK,ISC_GKK,ISC_GQ,
     '            ISR_GK,ISR_GKK,ISR_GQ,NBH,NBJ,NDET,ne,NENP,NGAP,
     '            NKH,NKHE,NKJE,NLL,NONY,%VAL(NPB_PTR),NP_INTERFACE,NPF,
     '            NPNE,NPNODE(0,nr),NPNY,nr,nr_gkk,NRE,NVJE,NW,nx,NYNE,
     '            NYNP,CE,CONY,CURVCORRECT,DET,DL,GK,GKK,GQ,PG,
     '            SE,WG,XIG,XP,INTERFACE,FULL_EQUATIONS,ERROR,*100)

                CALL FREE_MEMORY(NPB_PTR,ERROR,*100)

              ELSE

C CPB 18/10/96 Intialise pointers so they are zero in the parallel loop

                LGKE_PTR=0
                LGQE_PTR=0
                NPB_PTR=0
                DET_ADAPT_PTR=0
                DRDN_PTR=0
                DRDNO_PTR=0
                GKES_PTR=0
                GQES_PTR=0
                PG_J_PTR=0
                PG_Q_PTR=0
                PG_U_PTR=0
                RAD_PTR=0
                RD_PTR=0
                RG_PTR=0
                XA_PTR=0
                XE_PTR=0
                XG_PTR=0
                XG1_PTR=0
                XIG_J_PTR=0
                XIG_Q_PTR=0
                XIG_U_PTR=0
                XN_PTR=0
                XN_GRAD_PTR=0
                XR_PTR=0
                XR_GRAD_PTR=0

c cpb 15/10/96 Dynamically allocate the local arrays for the parallel
C loop
                CALL ALLOCATE_MEMORY((NKM+NHM*NSM+1)*NKM*2,1,INTTYPE,
     '            LGKE_PTR,MEM_INIT,ERROR,*100)
                CALL ALLOCATE_MEMORY((NKM+NHM*NSM+1)*NKM*2,1,INTTYPE,
     '            LGQE_PTR,MEM_INIT,ERROR,*100)
                CALL ALLOCATE_MEMORY((NP_R_M+1)*5,1,INTTYPE,NPB_PTR,
     '            MEM_INIT,ERROR,*100)
                CALL ALLOCATE_MEMORY(NGM,1,DPTYPE,DRDN_PTR,
     '            MEM_INIT,ERROR,*100)
                CALL ALLOCATE_MEMORY(NGM*NKM,1,DPTYPE,DRDNO_PTR,
     '            MEM_INIT,ERROR,*100)
                CALL ALLOCATE_MEMORY((NKM+NHM*NSM+1)*NKM,1,DPTYPE,
     '            GKES_PTR,MEM_INIT,ERROR,*100)
                CALL ALLOCATE_MEMORY((NKM+NHM*NSM+1)*NKM,1,DPTYPE,
     '            GQES_PTR,MEM_INIT,ERROR,*100)
                CALL ALLOCATE_MEMORY(NGM,1,DPTYPE,RAD_PTR,
     '            MEM_INIT,ERROR,*100)
                CALL ALLOCATE_MEMORY(NGM,1,DPTYPE,RD_PTR,
     '            MEM_INIT,ERROR,*100)
                CALL ALLOCATE_MEMORY(NGM,1,DPTYPE,RG_PTR,
     '            MEM_INIT,ERROR,*100)
                CALL ALLOCATE_MEMORY(NAM*NJM*NQM,1,DPTYPE,XA_PTR,
     '            MEM_INIT,ERROR,*100)
                CALL ALLOCATE_MEMORY(NSM*NJM,1,DPTYPE,XE_PTR,
     '            MEM_INIT,ERROR,*100)
                CALL ALLOCATE_MEMORY(NJM*NUM,1,DPTYPE,XG_PTR,
     '            MEM_INIT,ERROR,*100)
                CALL ALLOCATE_MEMORY(NJM*NUM*NGM,1,DPTYPE,XG1_PTR,
     '            MEM_INIT,ERROR,*100)
                CALL ALLOCATE_MEMORY(NJM*NGM,1,DPTYPE,XN_PTR,
     '            MEM_INIT,ERROR,*100)
                CALL ALLOCATE_MEMORY(NJM*NGM,1,DPTYPE,XN_GRAD_PTR,
     '            MEM_INIT,ERROR,*100)
                CALL ALLOCATE_MEMORY(NJM*NGM,1,DPTYPE,XR_PTR,
     '            MEM_INIT,ERROR,*100)
                CALL ALLOCATE_MEMORY(NJM*NGM,1,DPTYPE,XR_GRAD_PTR,
     '            MEM_INIT,ERROR,*100)
C!!! cpb 23/1/97 Don't allocate the PG_X etc. arrays if not using
C!!! adaptive integration. This will mean an access to address 0 if
C!!! these arrays are used inside the integration routines but it saves
C!!! memory.
                IF(ADAPINT) THEN
                  CALL ALLOCATE_MEMORY(NBFM*(NNM+1)*NGM,1,DPTYPE,
     '              DET_ADAPT_PTR,MEM_INIT,ERROR,*100)
                  CALL ALLOCATE_MEMORY(NSM*NUM*NGM,1,DPTYPE,PG_J_PTR,
     '              MEM_INIT,ERROR,*100)
                  CALL ALLOCATE_MEMORY(NSM*NUM*NGM,1,DPTYPE,PG_Q_PTR,
     '              MEM_INIT,ERROR,*100)
                  CALL ALLOCATE_MEMORY(NSM*NUM*NGM,1,DPTYPE,PG_U_PTR,
     '              MEM_INIT,ERROR,*100)
                  CALL ALLOCATE_MEMORY(NIM*NGM,1,DPTYPE,XIG_J_PTR,
     '              MEM_INIT,ERROR,*100)
                  CALL ALLOCATE_MEMORY(NIM*NGM,1,DPTYPE,XIG_Q_PTR,
     '              MEM_INIT,ERROR,*100)
                  CALL ALLOCATE_MEMORY(NIM*NGM,1,DPTYPE,XIG_U_PTR,
     '              MEM_INIT,ERROR,*100)
                ENDIF

                CALL XEGKGQ(IBT,IDO,INP,ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,
     '            ISR_GKK,ISR_GQ,%VAL(LGKE_PTR),%VAL(LGQE_PTR),NBH,
     '            NBJ,NDET,ne,NENP,NGAP,NHE,NHP,NKH,NKHE,NKJE,NLL,NONY,
     '            %VAL(NPB_PTR),NP_INTERFACE,NPF,NPNE,NPNODE(0,nr),
     '            NPNY,nr,nr_gkk,NRE,NVJE,NW,nx,NYNE,NYNP,CE,CONY,
     '            CURVCORRECT,DET,%VAL(DET_ADAPT_PTR),DL,%VAL(DRDN_PTR),
     '            %VAL(DRDNO_PTR),GK,GKK,%VAL(GKES_PTR),
     '            GQ,%VAL(GQES_PTR),PG,%VAL(PG_J_PTR),%VAL(PG_Q_PTR),
     '            %VAL(PG_U_PTR),%VAL(RAD_PTR),%VAL(RD_PTR),
     '            %VAL(RG_PTR),SE,WG,%VAL(XA_PTR),%VAL(XE_PTR),
     '            %VAL(XG_PTR),%VAL(XG1_PTR),XIG,%VAL(XIG_J_PTR),
     '            %VAL(XIG_Q_PTR),%VAL(XIG_U_PTR),%VAL(XN_PTR),
     '            %VAL(XN_GRAD_PTR),XP,%VAL(XR_PTR),%VAL(XR_GRAD_PTR),
     '            INTERFACE,FULL_EQUATIONS,ERROR,*100)

C cpb 15/10/96 Free dynamically allocated arrays

                CALL FREE_MEMORY(LGKE_PTR,ERROR,*100)
                CALL FREE_MEMORY(LGQE_PTR,ERROR,*100)
                CALL FREE_MEMORY(NPB_PTR,ERROR,*100)
                CALL FREE_MEMORY(DRDN_PTR,ERROR,*100)
                CALL FREE_MEMORY(DRDNO_PTR,ERROR,*100)
                CALL FREE_MEMORY(GKES_PTR,ERROR,*100)
                CALL FREE_MEMORY(GQES_PTR,ERROR,*100)
                CALL FREE_MEMORY(RAD_PTR,ERROR,*100)
                CALL FREE_MEMORY(RD_PTR,ERROR,*100)
                CALL FREE_MEMORY(RG_PTR,ERROR,*100)
                CALL FREE_MEMORY(XA_PTR,ERROR,*100)
                CALL FREE_MEMORY(XE_PTR,ERROR,*100)
                CALL FREE_MEMORY(XG_PTR,ERROR,*100)
                CALL FREE_MEMORY(XG1_PTR,ERROR,*100)
                CALL FREE_MEMORY(XN_PTR,ERROR,*100)
                CALL FREE_MEMORY(XN_GRAD_PTR,ERROR,*100)
                CALL FREE_MEMORY(XR_PTR,ERROR,*100)
                CALL FREE_MEMORY(XR_GRAD_PTR,ERROR,*100)
                IF(ADAPINT) THEN
                  CALL FREE_MEMORY(DET_ADAPT_PTR,ERROR,*100)
                  CALL FREE_MEMORY(PG_J_PTR,ERROR,*100)
                  CALL FREE_MEMORY(PG_Q_PTR,ERROR,*100)
                  CALL FREE_MEMORY(PG_U_PTR,ERROR,*100)
                  CALL FREE_MEMORY(XIG_J_PTR,ERROR,*100)
                  CALL FREE_MEMORY(XIG_Q_PTR,ERROR,*100)
                  CALL FREE_MEMORY(XIG_U_PTR,ERROR,*100)
                ENDIF

              ENDIF

            CALL CPU_TIMER(CPU_USER,TIME_STOP)
            ELAPSED_TIME=TIME_STOP(1)-TIME_START3(1)
            AVETIME=AVETIME+ELAPSED_TIME
            NO_TIMES=NO_TIMES+1
            IF(IWRIT4(nr,nx).GE.1) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP         CRITICAL(ASSEMBLE2_2)

C LKC 27-JAN-1998 Changing output format
C             WRITE(OP_STRING,'(/'' CPU time for element '',I5,'
C             '          //''' assembly: '',D11.4,'' s'')') ne,
C             ELAPSED_TIME
C             CALL WRITES(IOOP,OP_STRING,ERROR,*100)

              WRITE(OP_STRING,'('' CPU time for element '',I5,'
     '          //''' assembly: '',D11.4,'' s'')') ne,ELAPSED_TIME
              CALL WRITES(IOOP,OP_STRING,ERROR,*100)
CC$OMP         END CRITICAL(ASSEMBLE2_2)
            ENDIF

            GO TO 102
C               This statement is designed to be skipped if no error
C               occur. However if a error occurs within a subroutine
C               the alternate return points to line 100 to set the flag
 100            CONTINUE
C$OMP CRITICAL(ASSEMBLE2_1)
                ERROR_FLAG=.TRUE.
                WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
                CALL WRITES(IOER,OP_STRING,ERROR,*101)
                WRITE(OP_STRING,'(/'' >>An error occurred - '
     '            //'results may be unreliable!'')')
                CALL WRITES(IOER,OP_STRING,ERROR,*101)
 101            CONTINUE
C$OMP END CRITICAL(ASSEMBLE2_1)
 102          CONTINUE
            ENDIF

          ENDDO !ne
C$OMP END PARALLEL DO
        ELSE IF(BEMLOOPTYPE.EQ.2) THEN !Node outer loop
C$OMP     PARALLEL DO
C$OMP&      PRIVATE(nonode,np,INTERFACE,TIME_START3,TIME_STOP,
C$OMP&        ELAPSED_TIME,LGKE_PTR,LGQE_PTR,NPB_PTR,DET_ADAPT_PTR,
C$OMP&        DRDN_PTR,DRDNO_PTR,GKES_PTR,GQES_PTR,PG_J_PTR,
C$OMP&        PG_Q_PTR,PG_U_PTR,RAD_PTR,RD_PTR,RG_PTR,XA_PTR,XE_PTR,
C$OMP&        XG_PTR,XG1_PTR,XIG_J_PTR,XIG_Q_PTR,XIG_U_PTR,XN_PTR,
C$OMP&        XN_GRAD_PTR,XR_PTR,XR_GRAD_PTR,ERROR),
C$OMP&      SHARED(FULL_EQUATIONS,IOOP,MEM_INIT,NGM,nr,nr_gkk,nx,
C$OMP&        ERROR_FLAG),
C$OMP&      REDUCTION(+:AVETIME,NO_TIMES)
          DO nonode=1,NPNODE(0,nr)
            CALL CPU_TIMER(CPU_USER,TIME_START3)
            np=NPNODE(nonode,nr)
            IF(.NOT.ERROR_FLAG) THEN

C             INTERFACE is true if node np is on the interface
C             between regions in a coupled problem and the region to
C             which np belongs is NOT the region with the smallest
C             region number.
              INTERFACE=(NP_INTERFACE(np,0).GT.1).AND.
     '          (NP_INTERFACE(np,1).NE.nr)
              IF(OPTI3DLAPLACE) THEN

C cpb 23/9/98 Don't dynamically allocate local parallel arrays but
C hard code the dimensions in XPGKGQ_3DL

                CALL XPGKGQ_3DL(IBT,IDO,INP,ISC_GK,ISC_GKK,ISC_GQ,
     '            ISR_GK,ISR_GKK,ISR_GQ,NBH,NBJ,NDET,NEELEM,NENP,
     '            NGAP,NHP,NKH,NKHE,NKJE,NLL,NONY,np,NP_INTERFACE,NPF,
     '            NPNE,NPNY,nr,nr_gkk,NRE,NVJE,NW,nx,NYNE,NYNP,CE,CONY,
     '            CURVCORRECT,DET,DL,GK,GKK,GQ,PG,SE,WG,XIG,XP,
     '            INTERFACE,FULL_EQUATIONS,ERROR,*200)

              ELSE

C CPB 18/10/96 Intialise pointers so they are zero in the parallel loop

                LGKE_PTR=0
                LGQE_PTR=0
                DET_ADAPT_PTR=0
                DRDN_PTR=0
                DRDNO_PTR=0
                GKES_PTR=0
                GQES_PTR=0
                PG_J_PTR=0
                PG_Q_PTR=0
                PG_U_PTR=0
                RAD_PTR=0
                RD_PTR=0
                RG_PTR=0
                XA_PTR=0
                XE_PTR=0
                XG_PTR=0
                XG1_PTR=0
                XIG_J_PTR=0
                XIG_Q_PTR=0
                XIG_U_PTR=0
                XN_PTR=0
                XN_GRAD_PTR=0
                XR_PTR=0
                XR_GRAD_PTR=0

c cpb 15/10/96 Dynamically allocate the local arrays for the parallel
C loop
                CALL ALLOCATE_MEMORY((NKM+NHM*NSM+1)*NKM*2,1,INTTYPE,
     '            LGKE_PTR,MEM_INIT,ERROR,*200)
                CALL ALLOCATE_MEMORY((NKM+NHM*NSM+1)*NKM*2,1,INTTYPE,
     '            LGQE_PTR,MEM_INIT,ERROR,*200)
                CALL ALLOCATE_MEMORY(NGM,1,DPTYPE,DRDN_PTR,
     '            MEM_INIT,ERROR,*200)
                CALL ALLOCATE_MEMORY(NGM*NKM,1,DPTYPE,DRDNO_PTR,
     '            MEM_INIT,ERROR,*200)
                CALL ALLOCATE_MEMORY((NKM+NHM*NSM+1)*NKM,1,DPTYPE,
     '            GKES_PTR,MEM_INIT,ERROR,*200)
                CALL ALLOCATE_MEMORY((NKM+NHM*NSM+1)*NKM,1,DPTYPE,
     '            GQES_PTR,MEM_INIT,ERROR,*200)
                CALL ALLOCATE_MEMORY(NGM,1,DPTYPE,RAD_PTR,
     '            MEM_INIT,ERROR,*200)
                CALL ALLOCATE_MEMORY(NGM,1,DPTYPE,RD_PTR,
     '            MEM_INIT,ERROR,*200)
                CALL ALLOCATE_MEMORY(NGM,1,DPTYPE,RG_PTR,
     '            MEM_INIT,ERROR,*200)
                CALL ALLOCATE_MEMORY(NAM*NJM*NQM,1,DPTYPE,XA_PTR,
     '            MEM_INIT,ERROR,*200)
                CALL ALLOCATE_MEMORY(NSM*NJM,1,DPTYPE,XE_PTR,
     '            MEM_INIT,ERROR,*200)
                CALL ALLOCATE_MEMORY(NJM*NUM,1,DPTYPE,XG_PTR,
     '            MEM_INIT,ERROR,*200)
                CALL ALLOCATE_MEMORY(NJM*NUM*NGM,1,DPTYPE,XG1_PTR,
     '            MEM_INIT,ERROR,*200)
                CALL ALLOCATE_MEMORY(NJM*NGM,1,DPTYPE,XN_PTR,
     '            MEM_INIT,ERROR,*200)
                CALL ALLOCATE_MEMORY(NJM*NGM,1,DPTYPE,XN_GRAD_PTR,
     '            MEM_INIT,ERROR,*200)
                CALL ALLOCATE_MEMORY(NJM*NGM,1,DPTYPE,XR_PTR,
     '            MEM_INIT,ERROR,*200)
                CALL ALLOCATE_MEMORY(NJM*NGM,1,DPTYPE,XR_GRAD_PTR,
     '            MEM_INIT,ERROR,*200)
C!!! cpb 23/1/97 Don't allocate the PG_X and XIG_X arrays if not using
C!!! adaptive integration. This will mean an access to address 0 if
C!!! these arrays are used inside the integration routines but it saves
C!!! memory.
                IF(ADAPINT) THEN
                  CALL ALLOCATE_MEMORY(NBFM*(NNM+1)*NGM,1,DPTYPE,
     '              DET_ADAPT_PTR,MEM_INIT,ERROR,*200)
                  CALL ALLOCATE_MEMORY(NSM*NUM*NGM,1,DPTYPE,PG_J_PTR,
     '              MEM_INIT,ERROR,*200)
                  CALL ALLOCATE_MEMORY(NSM*NUM*NGM,1,DPTYPE,PG_Q_PTR,
     '              MEM_INIT,ERROR,*200)
                  CALL ALLOCATE_MEMORY(NSM*NUM*NGM,1,DPTYPE,PG_U_PTR,
     '              MEM_INIT,ERROR,*200)
                  CALL ALLOCATE_MEMORY(NIM*NGM,1,DPTYPE,XIG_J_PTR,
     '              MEM_INIT,ERROR,*200)
                  CALL ALLOCATE_MEMORY(NIM*NGM,1,DPTYPE,XIG_Q_PTR,
     '              MEM_INIT,ERROR,*200)
                  CALL ALLOCATE_MEMORY(NIM*NGM,1,DPTYPE,XIG_U_PTR,
     '              MEM_INIT,ERROR,*200)
                ENDIF

                CALL XPGKGQ(IBT,IDO,INP,ISC_GK,ISC_GKK,ISC_GQ,
     '            ISR_GK,ISR_GKK,ISR_GQ,%VAL(LGKE_PTR),%VAL(LGQE_PTR),
     '            NBH,NBJ,NDET,NEELEM,NENP,NGAP,NHE,NHP,NKH,NKHE,NKJE,
     '            NLL,NONY,np,NP_INTERFACE,NPF,NPNE,
     '            NPNY,nr,nr_gkk,NRE,NVJE,NW,nx,NYNE,NYNP,CE,CONY,
     '            CURVCORRECT,DET,%VAL(DET_ADAPT_PTR),DL,%VAL(DRDN_PTR),
     '            %VAL(DRDNO_PTR),GK,GKK,%VAL(GKES_PTR),GQ,
     '            %VAL(GQES_PTR),PG,%VAL(PG_J_PTR),%VAL(PG_Q_PTR),
     '            %VAL(PG_U_PTR),%VAL(RAD_PTR),%VAL(RD_PTR),
     '            %VAL(RG_PTR),SE,WG,%VAL(XA_PTR),%VAL(XE_PTR),
     '            %VAL(XG_PTR),%VAL(XG1_PTR),XIG,%VAL(XIG_J_PTR),
     '            %VAL(XIG_Q_PTR),%VAL(XIG_U_PTR),
     '            %VAL(XN_PTR),%VAL(XN_GRAD_PTR),XP,
     '            %VAL(XR_PTR),%VAL(XR_GRAD_PTR),INTERFACE,
     '            FULL_EQUATIONS,ERROR,*200)

C cpb 15/10/96 Free dynamically allocated arrays

                CALL FREE_MEMORY(LGKE_PTR,ERROR,*200)
                CALL FREE_MEMORY(LGQE_PTR,ERROR,*200)
                CALL FREE_MEMORY(DRDN_PTR,ERROR,*200)
                CALL FREE_MEMORY(DRDNO_PTR,ERROR,*200)
                CALL FREE_MEMORY(GKES_PTR,ERROR,*200)
                CALL FREE_MEMORY(GQES_PTR,ERROR,*200)
                CALL FREE_MEMORY(RAD_PTR,ERROR,*200)
                CALL FREE_MEMORY(RD_PTR,ERROR,*200)
                CALL FREE_MEMORY(RG_PTR,ERROR,*200)
                CALL FREE_MEMORY(XA_PTR,ERROR,*200)
                CALL FREE_MEMORY(XE_PTR,ERROR,*200)
                CALL FREE_MEMORY(XG_PTR,ERROR,*200)
                CALL FREE_MEMORY(XG1_PTR,ERROR,*200)
                CALL FREE_MEMORY(XN_PTR,ERROR,*200)
                CALL FREE_MEMORY(XN_GRAD_PTR,ERROR,*200)
                CALL FREE_MEMORY(XR_PTR,ERROR,*200)
                CALL FREE_MEMORY(XR_GRAD_PTR,ERROR,*200)
                IF(ADAPINT) THEN
                  CALL FREE_MEMORY(DET_ADAPT_PTR,ERROR,*200)
                  CALL FREE_MEMORY(PG_J_PTR,ERROR,*200)
                  CALL FREE_MEMORY(PG_Q_PTR,ERROR,*200)
                  CALL FREE_MEMORY(PG_U_PTR,ERROR,*200)
                  CALL FREE_MEMORY(XIG_J_PTR,ERROR,*200)
                  CALL FREE_MEMORY(XIG_Q_PTR,ERROR,*200)
                  CALL FREE_MEMORY(XIG_U_PTR,ERROR,*200)
                ENDIF

              ENDIF

              CALL CPU_TIMER(CPU_USER,TIME_STOP)
              ELAPSED_TIME=TIME_STOP(1)-TIME_START3(1)
              AVETIME=AVETIME+ELAPSED_TIME
              NO_TIMES=NO_TIMES+1
              IF(IWRIT4(nr,nx).GE.1) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP           CRITICAL(ASSEMBLE2_4)
                WRITE(OP_STRING,'(/'' CPU time for node '',I5,'
     '            //''' assembly: '',D11.4,'' s'')') np,ELAPSED_TIME
                CALL WRITES(IOOP,OP_STRING,ERROR,*200)
CC$OMP           END CRITICAL(ASSEMBLE2_4)
              ENDIF

              GO TO 202
C               This statement is designed to be skipped if no error
C               occur. However if a error occurs within a subroutine
C               the alternate return points to line 100 to set the flag
 200            CONTINUE
C$OMP CRITICAL(ASSEMBLE2_3)
                ERROR_FLAG=.TRUE.
                WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
                CALL WRITES(IOER,OP_STRING,ERROR,*201)
                WRITE(OP_STRING,'(/'' >>An error occurred - '
     '            //'results may be unreliable!'')')
                CALL WRITES(IODI,OP_STRING,ERROR,*201)
 201            CONTINUE
C$OMP END CRITICAL(ASSEMBLE2_3)
 202          CONTINUE
            ENDIF

          ENDDO !ne
C$OMP END PARALLEL DO
        ELSE
          ERROR='>>Invalid BEMLOOPTYPE'
          GOTO 9999
        ENDIF
        CALL ASSERT(.NOT.ERROR_FLAG,'>>An error occurred during '
     '    //'element stiffness calculations',ERROR,*9999)

        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
        IF(IWRIT4(nr,nx).GE.1) THEN
          WRITE(OP_STRING,'(/'' CPU time for stiffness matrix assembly '
     '      //': '',D11.4,'' s'')') ELAPSED_TIME
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          IF(BEMLOOPTYPE.EQ.1) THEN
            WRITE(OP_STRING,'('' Average CPU time per element for '
     '        //'assembly:'',D11.4,'' s'')') AVETIME/NO_TIMES
          ELSE
            WRITE(OP_STRING,'('' Average CPU time per node for '
     '        //'assembly:'',D11.4,'' s'')') AVETIME/NO_TIMES
          ENDIF
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

C Need to calculate the coefficient arising from the integration of the
C singularity of the Fundamental solution when the source is located at
C the surface nodes. For a smooth surface a "Laplace-type  equation"
C e.g. Helmholtz, Poisson, Yukawa) this term is just 1/2. For an
C arbitrary cartesian surface and a Laplace-type equation this term
C becomes: (interior angle between left and right tangents)/(2*pi) [2d]
C (internal solid angle)/(4*pi) [3d]

        CALL CPU_TIMER(CPU_USER,TIME_START2)
        CALL COEFF(INP,ISC_GK,ISC_GKK,ISR_GK,ISR_GKK,NBJ,NHP,NKJE,NONY,
     '    NPF,NPNE,NPNODE(0,nr),nr,nr_gkk,NRE,NVHP,NVJE,NW,nx,NYNP,
     '    CONY,GK,GKK,PG,SE,XA,XE,XG,XP,ERROR,*9999)

        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
        IF(IWRIT4(nr,nx).GE.1) THEN
          WRITE(OP_STRING,'(/'' CPU time for coefficient '
     '      //'calculation: '',D11.4,'' s'')') ELAPSED_TIME
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

      ENDIF


C***  Calculate GD vector (domain integrals or source terms).

C*** LKC 28-NOV-1999 Special source poisson problem
C***    should not go into this loop - USE_DIPOLE=0 or the
C***    number of dipoles will be 0

      IF(UPDATE_SOURCE.AND.(USE_DIPOLE.EQ.1)) THEN
        IF(IGREN(nr).GE.9.OR.NDIPOLES(nr,nx).GT.0) THEN
          !Poisson or/with dipole
          CALL CPU_TIMER(CPU_USER,TIME_START2)
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
              nj=NJ_LOC(NJL_GEOM,njj2,nr)
              XPFP(njj2)=XP(1,1,nj,np)
            ENDDO !njj2
            CALL XPGD(DIPOLE_CEN_NTIME(1,1,nx),
     '        DIPOLE_DIR_NTIME(1,1,nx),
     '        NBH,NDIPOLES(1,nx),NENP,NKH,NP_INTERFACE,np,nr,NW,
     '        nx,NYNP,CE,DIPOLE_CEN(1,0,1,1,nx),
     '        DIPOLE_DIR(1,0,1,1,nx),GD,TIME,XG,XP,XPFP,
     '        ERROR,*9999)
          ENDDO !nonode
          CALL CPU_TIMER(CPU_USER,TIME_STOP)
          ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
          IF(IWRIT4(nr,nx).GE.1) THEN
            WRITE(OP_STRING,'(/'' CPU time for source term '
     '        //'calculation: '',D11.4,'' s'')') ELAPSED_TIME
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF !IGREN, NDIPOLES
      ENDIF !UPDATE_SOURCE


      IF(UPDATE_MATRIX) THEN
        IF(IWRIT4(nr,nx).GE.4) THEN
          CALL CPU_TIMER(CPU_USER,TIME_START2)
          WRITE(OP_STRING,'(/'' Global stiffness matrix GK:'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NYNR(0,1,1)='',I5,'', NYNR(0,2,1)='','
     '      //'I5)') NYNR(0,1,1),NYNR(0,2,1)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
c cpb 24/9/95 Adding generic stiffness matrix output
          CALL OPSTFMAT(NYNR(0,2,1),ISC_GK,ISR_GK,IOOP,NYT(1,1,nx),
     '      NYT(2,1,nx),NZT(1,nx),NYNR(0,1,1),KTYP24,GK,GK,'GK ','   ',
     '      .FALSE.,.TRUE.,.FALSE.,ERROR,*9999)
          WRITE(OP_STRING,'(/'' Global stiffness matrix GQ:'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NYNR(0,1,2)='',I5,'', NYNR(0,2,1)='','
     '      //'I5)') NYNR(0,1,2),NYNR(0,2,2)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
c cpb 24/9/95 Adding generic stiffness matrix output
          CALL OPSTFMAT(NYNR(0,2,2),ISC_GQ,ISR_GQ,IOOP,NYT(1,2,nx),
     '      NYT(2,2,nx),NZT(2,nx),NYNR(0,1,2),KTYP24,GQ,GQ,'GQ ','   ',
     '      .FALSE.,.TRUE.,.FALSE.,ERROR,*9999)

          CALL CPU_TIMER(CPU_USER,TIME_STOP)
          ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
          WRITE(OP_STRING,'(/'' CPU time for stiffness matrices '
     '      //'output: '',D11.4,'' s'')') ELAPSED_TIME
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START1(1)
      IF(IWRIT4(nr,nx).GE.1) THEN
        WRITE(OP_STRING,'(/'' Total CPU time for assembly: '','
     '    //'D11.4,'' s'')') ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('ASSEMBLE2')
      RETURN
 9999 IF(NPB_PTR.NE.0) CALL FREE_MEMORY(NPB_PTR,ERROR_DUMMY,*1114)
      IF(DET_ADAPT_PTR.NE.0) CALL FREE_MEMORY(DET_ADAPT_PTR,
     '  ERROR_DUMMY,*1114)
      IF(DRDN_PTR.NE.0) CALL FREE_MEMORY(DRDN_PTR,ERROR_DUMMY,*1114)
      IF(DRDNO_PTR.NE.0) CALL FREE_MEMORY(DRDNO_PTR,ERROR_DUMMY,*1114)
      IF(PG_J_PTR.NE.0) CALL FREE_MEMORY(PG_J_PTR,ERROR_DUMMY,*1114)
      IF(PG_Q_PTR.NE.0) CALL FREE_MEMORY(PG_Q_PTR,ERROR_DUMMY,*1114)
      IF(PG_U_PTR.NE.0) CALL FREE_MEMORY(PG_U_PTR,ERROR_DUMMY,*1114)
      IF(RAD_PTR.NE.0) CALL FREE_MEMORY(RAD_PTR,ERROR_DUMMY,*1114)
      IF(RD_PTR.NE.0) CALL FREE_MEMORY(RD_PTR,ERROR_DUMMY,*1114)
      IF(RG_PTR.NE.0) CALL FREE_MEMORY(RG_PTR,ERROR_DUMMY,*1114)
      IF(XA_PTR.NE.0) CALL FREE_MEMORY(XA_PTR,ERROR_DUMMY,*1114)
      IF(XE_PTR.NE.0) CALL FREE_MEMORY(XE_PTR,ERROR_DUMMY,*1114)
      IF(XG_PTR.NE.0) CALL FREE_MEMORY(XG_PTR,ERROR_DUMMY,*1114)
      IF(XG1_PTR.NE.0) CALL FREE_MEMORY(XG1_PTR,ERROR_DUMMY,*1114)
      IF(XIG_J_PTR.NE.0) CALL FREE_MEMORY(XIG_J_PTR,ERROR_DUMMY,*1114)
      IF(XIG_Q_PTR.NE.0) CALL FREE_MEMORY(XIG_Q_PTR,ERROR_DUMMY,*1114)
      IF(XIG_U_PTR.NE.0) CALL FREE_MEMORY(XIG_U_PTR,ERROR_DUMMY,*1114)
      IF(XN_PTR.NE.0) CALL FREE_MEMORY(XN_PTR,ERROR_DUMMY,*1114)
      IF(XN_GRAD_PTR.NE.0) CALL FREE_MEMORY(XN_GRAD_PTR,ERROR_DUMMY,
     '  *1114)
      IF(XR_PTR.NE.0) CALL FREE_MEMORY(XR_PTR,ERROR_DUMMY,*1114)
      IF(XR_GRAD_PTR.NE.0) CALL FREE_MEMORY(XR_GRAD_PTR,ERROR_DUMMY,
     '  *1114)

 1114 CALL ERRORS('ASSEMBLE2',ERROR)
      RET_ERROR=ERROR
      CALL EXITS('ASSEMBLE2')
      RETURN 1
      END


