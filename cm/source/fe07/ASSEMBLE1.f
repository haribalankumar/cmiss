      SUBROUTINE ASSEMBLE1(IBT,IDO,INP,ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,
     '  ISR_GKK,ISR_GQ,NBH,NBJ,NEELEM,NHE,NHP,NKJE,NONY,NORD,NPF,
     '  NP_INTERFACE,NPNE,NPNY,nr,nr_gkk,NRE,NVHE,NVJE,NW,nx,NYNE,NYNP,
     '  NYNR,CE,CGE,CONY,CP,CURVCORRECT,GK,GKK,GQ,GR,PG,SE,WG,XA,
     '  XP,YG,GQ_ASSEM,UPDATE_MATRIX,UPDATE_VECTOR,RET_ERROR,*)

C#### Subroutine: ASSEMBLE1
C###  Description:
C###    ASSEMBLE1 assembles the global unreduced matrices GK,
C###    GM, etc. for Static linear FEM problems.

C**** If GQ_ASSEM is .TRUE. then the GQ matrix is assembled for the
C**** nodes that are shared between appropriate regions.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'error0.inc'
      INCLUDE 'fsklib.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'host00.inc'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'mxch.inc'
      INCLUDE 'solv00.cmn'
      INCLUDE 'time02.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),ISC_GK(NISC_GKM),ISC_GKK(NISC_GKKM),
     '  ISC_GQ(NISC_GQM),ISR_GK(NISR_GKM),ISR_GKK(NISR_GKKM),
     '  ISR_GQ(NISR_GQM),NBH(NHM,NCM,NEM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NHE(NEM),NHP(NPM),
     '  NKJE(NKM,NNM,NJM,NEM),NONY(0:NOYM,NYM,NRCM),NORD(5,NE_R_M),
     '  NP_INTERFACE(0:NPM,0:3),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NPNY(0:6,NYM,0:NRCM),nr,nr_gkk,NRE(NEM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVJE(NNM,NBFM,NJM,NEM),
     '  NW(NEM,3),nx,NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM)
      REAL*8 CE(NMM,NEM),CGE(NMM,NGM,NEM),
     '  CONY(0:NOYM,NYM,NRCM),CP(NMM,NPM),
     '  CURVCORRECT(2,2,NNM,NEM),GK(NZ_GK_M),GKK(NZ_GKK_M),
     '  GQ(NZ_GQ_M),GR(NYROWM),PG(NSM,NUM,NGM,NBM),
     '  SE(NSM,NBFM,NEM),WG(NGM,NBM),XA(NAM,NJM,NEM),
     '  XP(NKM,NVM,NJM,NPM),YG(NIYGM,NGM,NEM)
      LOGICAL GQ_ASSEM,UPDATE_MATRIX,UPDATE_VECTOR
      CHARACTER RET_ERROR*(*)
!     Local Variables
      INTEGER ne,noelem,no_nynr1,NO_TIMES,ny1
C CS 6/5/97 Adding dynamic allocation for parallel local arrays
      INTEGER*4 LGE_PTR,CG_PTR,ED_PTR,EM_PTR,ER_PTR,ES_PTR,
     '  RG_PTR,XE_PTR,XG_PTR,XG1_PTR,ZE_PTR,ZG_PTR
      !SMAR009 22/12/98 ,WORK_PTR
      REAL ELAPSED_TIME,TIME_START1(1),TIME_START2(1),TIME_START3(1),
     '  TIME_STOP(1)
      CHARACTER ERROR*(ERRSTRLEN)
      LOGICAL ERROR_FLAG

      CALL ENTERS('ASSEMBLE1',*9999)

C CS 6/5/97 Intialise pointers so the error free will work properly

      LGE_PTR=0
      CG_PTR=0
      ED_PTR=0
      EM_PTR=0
      ER_PTR=0
      ES_PTR=0
      RG_PTR=0
      XE_PTR=0
      XG_PTR=0
      XG1_PTR=0
      ZE_PTR=0
      ZG_PTR=0

C cpb 27/4/95 Adding sparse global matrices

       CALL CPU_TIMER(CPU_USER,TIME_START1)

      IF(UPDATE_MATRIX.OR.UPDATE_VECTOR) THEN

C***    Initialise matrices for this region

        CALL CPU_TIMER(CPU_USER,TIME_START2)

        IF(UPDATE_MATRIX) THEN

          CALL INIT_SPARSE_MATRIX(NYNR(0,2,1),NYNR(0,2,1),ISC_GK,
     '      ISR_GK,0,0,1,NYT(1,1,nx),NZ_GK_M,NZT(1,nx),NYNR(0,1,1),
     '      NYNR(0,1,1),KTYP24,GK,ERROR,*9999)
          IF(GQ_ASSEM) THEN
            CALL INIT_SPARSE_MATRIX(NYNR(0,2,2),NYNR(0,2,2),ISC_GQ,
     '        ISR_GQ,0,0,1,NYT(1,2,nx),NZ_GQ_M,NZT(2,nx),NYNR(0,1,2),
     '        NYNR(0,1,2),KTYP24,GQ,ERROR,*9999)
          ENDIF

          IF(.NOT.CALC_GLOBAL(nr,nx)) THEN
            CALL INIT_SPARSE_MATRIX(NYNR(0,2,1),NONY(0,1,2),ISC_GKK,
     '        ISR_GKK,NOYM,NOYM,2,NOT(1,1,nr,nx),NZ_GKK_M,
     '        NZZT(1,nr,nx),NYNR(0,1,1),NONY(0,1,1),SPARSEGKK(nx),
     '        GKK,ERROR,*9999)
            IF(GQ_ASSEM) THEN
              CALL INIT_SPARSE_MATRIX(NYNR(0,2,2),NONY(0,1,2),ISC_GQ,
     '          ISR_GQ,NOYM,NOYM,2,NOT(1,1,nr,nx),NZ_GKK_M,
     '          NZZT(1,nr,nx),NYNR(0,1,2),NONY(0,1,1),SPARSEGKK(nx),
     '          GKK,ERROR,*9999)
            ENDIF
          ENDIF

        ENDIF
C        IF(UPDATE_VECTOR.AND..NOT.GQ_ASSEM) THEN
        IF(UPDATE_VECTOR) THEN
          CALL ASSERT(NYT(1,1,nx).LE.NYROWM,'>>Increase NYROWM',
     '      ERROR,*9999)
C$OMP     PARALLEL DO
C$OMP&      PRIVATE(no_nynr1,ny1),
C$OMP&      SHARED(NYNR,GR)
          DO no_nynr1=1,NYNR(0,1,1)
            ny1=NYNR(no_nynr1,1,1)
            GR(ny1)=0.0d0
          ENDDO !no_nynr (ny1)
C$OMP     END PARALLEL DO
        ENDIF

        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
        IF(IWRIT4(nr,nx).GE.1) THEN
          WRITE(OP_STRING,'(/'' CPU time for setup and '
     '      //'initialisation: '',D11.4,'' s'')') ELAPSED_TIME
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

C***  Find element stiffness matrices

        CALL CPU_TIMER(CPU_USER,TIME_START2)
        NO_TIMES=0

        ERROR_FLAG=.FALSE.
C CS 6/5/97 Dynamically allocation parallel local arrays
C$OMP   PARALLEL DO
C$OMP&    PRIVATE(ELAPSED_TIME,
C$OMP&            ne,noelem,LGE_PTR,CG_PTR,ED_PTR,EM_PTR,
C$OMP&            ES_PTR,ER_PTR,RG_PTR,XE_PTR,XG_PTR,XG1_PTR,
C$OMP&            ZE_PTR,ZG_PTR,TIME_START3,TIME_STOP,ERROR),
C$OMP&    SHARED(GQ_ASSEM,MEM_INIT,NGM,nr,nr_gkk,nx,ERROR_FLAG,
C$OMP&           UPDATE_MATRIX,UPDATE_VECTOR),
C$OMP&    REDUCTION(+:NO_TIMES)
        DO noelem=1,NEELEM(0,nr)
          CALL CPU_TIMER(CPU_USER,TIME_START3)
          ne=NEELEM(noelem,nr)
          IF(.NOT.ERROR_FLAG) THEN

C CS 6/5/97 Intialise pointers so they are zero in the parallel loop

            LGE_PTR=0
            CG_PTR=0
            ED_PTR=0
            EM_PTR=0
            ER_PTR=0
            ES_PTR=0
            RG_PTR=0
            XE_PTR=0
            XG_PTR=0
            XG1_PTR=0
            ZE_PTR=0
            ZG_PTR=0

C CS 6/5/97  Dynamically allocation parallel local arrays

            CALL ALLOCATE_MEMORY(NHM*NSM*NRCM,1,INTTYPE,LGE_PTR,
     '        MEM_INIT,ERROR,*100)
            CALL ALLOCATE_MEMORY(NMM*NGM,1,DPTYPE,CG_PTR,
     '        MEM_INIT,ERROR,*100)
            CALL ALLOCATE_MEMORY(NHM*NSM*NHM*NSM,1,DPTYPE,ED_PTR,
     '        MEM_INIT,ERROR,*100)
            CALL ALLOCATE_MEMORY(NHM*NSM*NHM*NSM,1,DPTYPE,EM_PTR,
     '        MEM_INIT,ERROR,*100)
            CALL ALLOCATE_MEMORY(NHM*NSM,1,DPTYPE,ER_PTR,
     '        MEM_INIT,ERROR,*100)
            CALL ALLOCATE_MEMORY(NHM*NSM*NHM*NSM,1,DPTYPE,ES_PTR,
     '        MEM_INIT,ERROR,*100)
            CALL ALLOCATE_MEMORY(NGM,1,DPTYPE,RG_PTR,
     '        MEM_INIT,ERROR,*100)
            CALL ALLOCATE_MEMORY(NSM*NJM,1,DPTYPE,XE_PTR,
     '        MEM_INIT,ERROR,*100)
            CALL ALLOCATE_MEMORY(NJM*NUM,1,DPTYPE,XG_PTR,
     '        MEM_INIT,ERROR,*100)
            CALL ALLOCATE_MEMORY(NJM*NUM*NGM,1,DPTYPE,XG1_PTR,
     '        MEM_INIT,ERROR,*100)
            CALL ALLOCATE_MEMORY(NSM*NHM,1,DPTYPE,ZE_PTR,
     '        MEM_INIT,ERROR,*100)
            CALL ALLOCATE_MEMORY(NHM*NUM,1,DPTYPE,ZG_PTR,
     '        MEM_INIT,ERROR,*100)

C CS 6/5/97 Must create a dynam subroutine as dynamically allocated
C arrays are used within this subroutine level. This can be fixed on
C the move to F90.

            CALL ASSEMBLE1_DYNAM(IBT,IDO,INP,ISC_GK,ISC_GKK,ISC_GQ,
     '        ISR_GK,ISR_GKK,ISR_GQ,%VAL(LGE_PTR),NBH,NBJ,ne,
     '        NEELEM,NHE,NHP,NKJE,NONY,NORD,NP_INTERFACE,NPF,NPNE,NPNY,
     '        nr,nr_gkk,NRE,NVHE,NVJE,NW,nx,NYNE,NYNP,CE,%VAL(CG_PTR),
     '        CGE,CONY,CP,CURVCORRECT,%VAL(ED_PTR),%VAL(EM_PTR),
     '        %VAL(ER_PTR),%VAL(ES_PTR),GK,GKK,GQ,GR,PG,
     '        %VAL(RG_PTR),SE,WG,XA,%VAL(XE_PTR),
     '        %VAL(XG_PTR),%VAL(XG1_PTR),XP,YG,%VAL(ZE_PTR),
     '        %VAL(ZG_PTR),GQ_ASSEM,TIME_START3,
     '        UPDATE_MATRIX,UPDATE_VECTOR,ERROR,*100)

C CS 6/5/97 Free dynamically allocated arrays

            CALL FREE_MEMORY(LGE_PTR,ERROR,*100)
            CALL FREE_MEMORY(CG_PTR,ERROR,*100)
            CALL FREE_MEMORY(ED_PTR,ERROR,*100)
            CALL FREE_MEMORY(EM_PTR,ERROR,*100)
            CALL FREE_MEMORY(ER_PTR,ERROR,*100)
            CALL FREE_MEMORY(ES_PTR,ERROR,*100)
            CALL FREE_MEMORY(RG_PTR,ERROR,*100)
            CALL FREE_MEMORY(XE_PTR,ERROR,*100)
            CALL FREE_MEMORY(XG_PTR,ERROR,*100)
            CALL FREE_MEMORY(XG1_PTR,ERROR,*100)
            CALL FREE_MEMORY(ZE_PTR,ERROR,*100)
            CALL FREE_MEMORY(ZG_PTR,ERROR,*100)

            NO_TIMES=NO_TIMES+1

            GO TO 102
C             This statement is designed to be skipped if no error
C             occurs. However if a error occurs within a subroutine
C             the alternate return points to line 100 to set the flag
 100          CONTINUE
C$OMP CRITICAL(ASSEMBLE1_1)
              ERROR_FLAG=.TRUE.
              WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
              CALL WRITES(IOER,OP_STRING,ERROR,*101)
              WRITE(OP_STRING,'(/'' >>An error occurred - '
     '          //'results may be unreliable!'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*101)
 101          CONTINUE
C$OMP END CRITICAL(ASSEMBLE1_1)
 102        CONTINUE
          ENDIF !.NOT.ERROR_FLAG
        ENDDO !noelem (ne)
C$OMP END PARALLEL DO

        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
        IF(IWRIT4(nr,nx).GE.1) THEN
          WRITE(OP_STRING,'(/'' CPU time for stiffness matrix assembly '
     '      //': '',D11.4,'' s'')') ELAPSED_TIME
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Average CPU time per element for '
     '      //'assembly:'',D11.4,'' s'')') ELAPSED_TIME/NO_TIMES
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

        IF(IWRIT4(nr,nx).GE.4) THEN
          CALL CPU_TIMER(CPU_USER,TIME_START2)
          IF(.NOT.GQ_ASSEM) THEN
            WRITE(OP_STRING,
     '        '(/'' Global load vector GR & stiffness matrix GK:'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' NYNR(0,1,1)='',I5,'', NYNR(0,2,1)='','
     '        //'I5)') NYNR(0,1,1),NYNR(0,2,1)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
c cpb 24/9/95 Adding generic stiffness matrix output
             CALL OPSTFMAT(NYNR(0,2,1),ISC_GK,ISR_GK,IOOP,NYT(1,1,nx),
     '        NYT(2,1,nx),NZT(1,nx),NYNR(0,1,1),KTYP24,GK,GR,'GK ',
     '        'GR ',.FALSE.,.TRUE.,.TRUE.,ERROR,*9999)
          ELSE
            WRITE(OP_STRING,'(/'' Global stiffness matrix GK:'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' NYNR(0,1,1)='',I5,'', NYNR(0,2,1)='','
     '        //'I5)') NYNR(0,1,1),NYNR(0,2,1)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
c cpb 24/9/95 Adding generic stiffness matrix output
            CALL OPSTFMAT(NYNR(0,2,1),ISC_GK,ISR_GK,IOOP,NYT(1,1,nx),
     '        NYT(2,1,nx),NZT(1,nx),NYNR(0,1,1),KTYP24,GK,GR,'GK ',
     '        'GR ',.FALSE.,.TRUE.,.FALSE.,ERROR,*9999)
            WRITE(OP_STRING,'(/'' Global stiffness matrix GQ:'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' NYNR(0,1,2)='',I5,'', NYNR(0,2,2)='','
     '        //'I5)') NYNR(0,1,2),NYNR(0,2,2)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
c cpb 24/9/95 Adding generic stiffness matrix output
            CALL OPSTFMAT(NYNR(0,2,2),ISC_GQ,ISR_GQ,IOOP,NYT(1,2,nx),
     '        NYT(2,2,nx),NZT(2,nx),NYNR(0,1,2),KTYP24,GQ,GR,'GQ ',
     '        'GR ',.FALSE.,.TRUE.,.FALSE.,ERROR,*9999)
          ENDIF

          CALL CPU_TIMER(CPU_USER,TIME_STOP)
          ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
          WRITE(OP_STRING,'(/'' CPU time for stiffness matrices '
     '      //'output: '',D11.4,'' s'')') ELAPSED_TIME
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

        ENDIF
      ENDIF !update

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START1(1)
      IF(IWRIT4(nr,nx).GE.1) THEN
        WRITE(OP_STRING,'(/'' Total CPU time for assembly: '','
     '    //'D11.4,'' s'')') ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('ASSEMBLE1')
      RETURN

 9999 IF(LGE_PTR.NE.0) CALL FREE_MEMORY(LGE_PTR,ERROR,*1114)
      IF(CG_PTR.NE.0) CALL FREE_MEMORY(CG_PTR,ERROR,*1114)
      IF(ES_PTR.NE.0) CALL FREE_MEMORY(ED_PTR,ERROR,*1114)
      IF(ES_PTR.NE.0) CALL FREE_MEMORY(EM_PTR,ERROR,*1114)
      IF(ES_PTR.NE.0) CALL FREE_MEMORY(ER_PTR,ERROR,*1114)
      IF(ES_PTR.NE.0) CALL FREE_MEMORY(ES_PTR,ERROR,*1114)
      IF(RG_PTR.NE.0) CALL FREE_MEMORY(RG_PTR,ERROR,*1114)
      IF(XE_PTR.NE.0) CALL FREE_MEMORY(XE_PTR,ERROR,*1114)
      IF(XG_PTR.NE.0) CALL FREE_MEMORY(XG_PTR,ERROR,*1114)
      IF(XG_PTR.NE.0) CALL FREE_MEMORY(XG1_PTR,ERROR,*1114)
      IF(ZE_PTR.NE.0) CALL FREE_MEMORY(ZE_PTR,ERROR,*1114)
      IF(ZG_PTR.NE.0) CALL FREE_MEMORY(ZG_PTR,ERROR,*1114)
 1114 CALL ERRORS('ASSEMBLE1',ERROR)
      RET_ERROR=ERROR
      CALL EXITS('ASSEMBLE1')
      RETURN 1
      END


