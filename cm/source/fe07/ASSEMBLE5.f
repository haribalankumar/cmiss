      SUBROUTINE ASSEMBLE5(IBT,IDO,INP,ISC_GK,ISR_GK,NAN,NBH,NBHF,
     '  NBJ,NBJF,NEELEM,NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,
     '  NKHE,NKJE,NMNO,NNB,NNF,NPF,NPNE,NPNODE,NPNY,NRE,NRLIST,
     '  nr_solve,NSB,NVHE,NVHP,NVJE,NW,nx,NXI,NYNE,NYNP,NYNR,
     '  CE,CGE,CP,CURVCORRECT,FEXT,GK,PG,SE,WG,XA,XP,YG,YGF,YP,
     '  ZA,ZA1,ZAA,ZP,ZP1,ZPA,FIX,RET_ERROR,*)

C#### Subroutine: ASSEMBLE5
C###  Description:
C###    ASSEMBLE5 generates element stiffness matrix ES for nonlinear
C###    problems by ZEES, and assembles into the global stiffness
C###    matrix GK.

C**** If KTYP1A=2 element stiffness matrices are computed in parallel
C**** using slave processes running on hosts set up in 'define solve'.
C**** If KTYP1A=1 elem stiff matrices are computed in series, locally.
C**** LGE(nhs,nrc) is location in global system of elem var nhs
C NEW 06/07/05 JHC parallel element stiffness matrix is no longer implemented

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'error0.inc'
      INCLUDE 'fsklib.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'host00.cmn'
      INCLUDE 'host00.inc'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'mxch.inc'
      INCLUDE 'time02.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISC_GK(NISC_GKM),ISR_GK(NISR_GKM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),
     '  NBJ(NJM,NEM),NBJF(NJM,NFM),NEELEM(0:NE_R_M,0:NRM),
     '  NFF(6,NEM),NFFACE(0:NF_R_M,NRM),NGAP(NIM,NBM),
     '  NHE(NEM),NHP(NPM,0:NRM),NKB(2,2,2,NNM,NBFM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NKEF(0:4,16,6,NBFM),
     '  NKH(NHM,NPM,NCM,0:NRM),NMNO(1:2,0:NOPM),NNB(4,4,4,NBFM),
     '  NNF(0:17,6,NBFM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM),
     '  NRE(NEM),NRLIST(0:NRM),nr_solve,NSB(NKM,NNM,NBFM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx,
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM)
      REAL*8 CE(NMM,NEM),CGE(NMM,NGM,NEM),CP(NMM,NPM),
     '  CURVCORRECT(2,2,NNM,NEM),
     '  FEXT(NIFEXTM,NGM,NEM),GK(NZ_GK_M),
     '  PG(NSM,NUM,NGM,NBM),
     '  SE(NSM,NBFM,NEM),WG(NGM,NBM),XA(NAM,NJM,NEM),
     '  XP(NKM,NVM,NJM,NPM),
     '  YG(NIYGM,NGM,NEM),YGF(NIYGFM,NGFM,NFM),YP(NYM,NIYM),
     '  ZA(NAM,NHM,NCM,NEM),ZA1(NAM,NHM,NCM,NEM),ZAA(NAM,NHM,NCM,NEM),
     '  ZP(NKM,NVM,NHM,NPM,NCM),ZP1(NKM,NVM,NHM,NPM,NCM),
     '  ZPA(NKM,NVM,NHM,NPM,NCM)
      CHARACTER RET_ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM)
!     Local Variables
      INTEGER ne,nf,noelem,noface,nonrlist,nr
c cpb 18/10/96 Adding dynamic allocation for parallel local arrays
      INTEGER*4 LGE_PTR,IDOXFT_PTR,MYMS_PTR,NSFE_PTR,NYNS_PTR,
     '  CG_PTR,D_RE_PTR,D_RI3_PTR,D_TG_PTR,D_ZG_PTR,ES_PTR,FS_PTR,
     '  RE1_PTR,RE2_PTR,RDF1_PTR,RDF2_PTR,RG_PTR,SM_PTR,SN_PTR,
     '  XDF_PTR,XE_PTR,XG_PTR,
     '  ZDF_PTR,ZE_PTR,ZE1_PTR,ZG_PTR,ZG1_PTR
      REAL ELAPSED_TIME,TIME_START1(1),
     '  TIME_START1_REAL(1),TIME_START2(1),TIME_STOP(1)
      CHARACTER ERROR*(ERRSTRLEN)
      LOGICAL ERROR_FLAG

      CALL ENTERS('ASSEMBLE5',*9999)

C CPB 18/10/96 Intialise pointers so the error free will work properly

      LGE_PTR=0
      CG_PTR=0
      D_RE_PTR=0
      D_RI3_PTR=0
      D_TG_PTR=0
      D_ZG_PTR=0
      ES_PTR=0
      RE1_PTR=0
      RE2_PTR=0
      RG_PTR=0
      XE_PTR=0
      XG_PTR=0
      ZE_PTR=0
      ZE1_PTR=0
      ZG_PTR=0
      ZG1_PTR=0

      IDOXFT_PTR=0
      MYMS_PTR=0
      NSFE_PTR=0
      NYNS_PTR=0
      FS_PTR=0
      RDF1_PTR=0
      RDF2_PTR=0
      SM_PTR=0
      SN_PTR=0
      XDF_PTR=0
      ZDF_PTR=0

C *** Initialise GK matrix for solution region nr_solve
      CALL INIT_SPARSE_MATRIX(NYNR(0,2,1,nr_solve),
     '  NYNR(0,2,1,nr_solve),ISC_GK,ISR_GK,0,0,1,
     '  NYT(1,1,nx),NZ_GK_M,NZT(1,nx),NYNR(0,1,1,nr_solve),
     '  NYNR(0,1,1,nr_solve),KTYP24,GK,ERROR,*9999)

C NEWS 06/07/05 JHC as parallel stiffness matrix calculation is no longer implemented,
C AVETIME is no longer needed
C      AVETIME=0.0
      CALL CPU_TIMER(CPU_USER,TIME_START1)
      CALL REAL_TIMER(REAL_TOTAL,TIME_START1_REAL)
      CALL REAL_TIMER(REAL_TOTAL,TIME_START2)

C *** Calculate Jacobian (global tangent stiffness matrix)
      IF(KTYP1A.EQ.1) THEN !calc elem stiff matrices in series locally

        DO nonrlist=1,NRLIST(0)
          nr=NRLIST(nonrlist)
          CALL YPZP(1,NBH,NEELEM,NHE,NHP(1,nr),NKH(1,1,1,nr),NPNODE,nr,
     '      NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
          IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN !cnst Vol
C           Put reference state for cavity from YP(ny,10) into
C           ZA1,ZP1 for ZERE55
            CALL YPZP(10,NBH,NEELEM,NHE,NHP(1,nr),NKH(1,1,1,nr),
     '        NPNODE,nr,
     '        NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZA1,ZP1,ERROR,*9999)
          ENDIF

          ERROR_FLAG=.FALSE.
C new MPN 1Feb2000: OMP parallel proc directive
C old
CC$DOACROSS local (ne,noelem,LGE_PTR,CG_PTR,D_RE_PTR,D_RI3_PTR,D_TG_PTR,
CC$&               D_ZG_PTR,ES_PTR,RE1_PTR,RE2_PTR,RG_PTR,XE_PTR,XG_PTR,
CC$&               ZE_PTR,ZE1_PTR,ZG_PTR,ZG1_PTR,ERROR)
CC$&        share (nr,nx,ERROR_FLAG)
C end old
C$OMP     PARALLEL DO
C$OMP&      PRIVATE(ne,noelem,LGE_PTR,CG_PTR,D_RE_PTR,D_RI3_PTR,
C$OMP&              D_TG_PTR,D_ZG_PTR,ES_PTR,RE1_PTR,RE2_PTR,RG_PTR,
C$OMP&              XE_PTR,XG_PTR,ZE_PTR,ZE1_PTR,ZG_PTR,ZG1_PTR,ERROR),
C$OMP&      SHARED(MEM_INIT,NGM,nr,nx,ERROR_FLAG)
          DO noelem=1,NEELEM(0,nr) !is main element loop
            ne=NEELEM(noelem,nr)
            IF(.NOT.ERROR_FLAG) THEN

C CPB 18/10/96 Intialise pointers so they are zero in the parallel loop

              LGE_PTR=0
              CG_PTR=0
              D_RE_PTR=0
              D_RI3_PTR=0
              D_TG_PTR=0
              D_ZG_PTR=0
              ES_PTR=0
              RG_PTR=0
              XE_PTR=0
              XG_PTR=0
              ZE_PTR=0
              ZE1_PTR=0
              ZG_PTR=0
              ZG1_PTR=0

c cpb 18/10/96  Dynamically allocation parallel local arrays

              CALL ALLOCATE_MEMORY(NHM*NSM*NRCM,1,INTTYPE,LGE_PTR,
     '          MEM_INIT,ERROR,*100)
              CALL ALLOCATE_MEMORY(NMM*NGM,1,DPTYPE,CG_PTR,
     '          MEM_INIT,ERROR,*100)
              CALL ALLOCATE_MEMORY(NSM*NHM*NOPM,1,DPTYPE,D_RE_PTR,
     '          MEM_INIT,ERROR,*100)
              CALL ALLOCATE_MEMORY(NHM*NSM,1,DPTYPE,D_RI3_PTR,
     '          MEM_INIT,ERROR,*100)
              CALL ALLOCATE_MEMORY(3*3*NHM*NSM,1,DPTYPE,D_TG_PTR,
     '          MEM_INIT,ERROR,*100)
              CALL ALLOCATE_MEMORY(NHM*NUM*NHM*NSM,1,DPTYPE,D_ZG_PTR,
     '          MEM_INIT,ERROR,*100)
              CALL ALLOCATE_MEMORY(NHM*NSM*NHM*NSM,1,DPTYPE,ES_PTR,
     '          MEM_INIT,ERROR,*100)
              CALL ALLOCATE_MEMORY(NGM,1,DPTYPE,RG_PTR,
     '          MEM_INIT,ERROR,*100)
              CALL ALLOCATE_MEMORY(NSM*NJM,1,DPTYPE,XE_PTR,
     '          MEM_INIT,ERROR,*100)
              CALL ALLOCATE_MEMORY(NJM*NUM,1,DPTYPE,XG_PTR,
     '          MEM_INIT,ERROR,*100)
              CALL ALLOCATE_MEMORY(NSM*NHM,1,DPTYPE,ZE_PTR,
     '          MEM_INIT,ERROR,*100)
              CALL ALLOCATE_MEMORY(NSM*NHM,1,DPTYPE,ZE1_PTR,
     '          MEM_INIT,ERROR,*100)
              CALL ALLOCATE_MEMORY(NHM*NUM,1,DPTYPE,ZG_PTR,
     '          MEM_INIT,ERROR,*100)
              CALL ALLOCATE_MEMORY(NHM*NUM,1,DPTYPE,ZG1_PTR,
     '          MEM_INIT,ERROR,*100)

C             These arrays are only necessary for finite differences
              IF(KTYP1D.NE.1) THEN
C              IF(.true..or.KTYP1D.NE.1) THEN
                RE1_PTR=0
                CALL ALLOCATE_MEMORY(NSM*NHM,1,DPTYPE,RE1_PTR,
     '            MEM_INIT,ERROR,*100)
                RE2_PTR=0
                CALL ALLOCATE_MEMORY(NSM*NHM,1,DPTYPE,RE2_PTR,
     '            MEM_INIT,ERROR,*100)
              ENDIF

C cpb 18/10/96 Must create a dynam subroutine as dynamically allocated
C arrays are used within this subroutine level. This can be fixed on
C the move to F90.

              CALL ASSEMBLE5_DYNAM(IBT,IDO,INP,ISC_GK,ISR_GK,
     '          %VAL(LGE_PTR),NAN,NBH,NBJ,NBJF,ne,NFF,NGAP,NHE,
     '          NKEF,NKHE,NKJE,NMNO,NNF,NPF,NPNE,NPNY,nr,NRE,NVHE,NVJE,
     '          NW,nx,NXI,NYNE,NYNP,CE,%VAL(CG_PTR),CGE,CP,CURVCORRECT,
     '          %VAL(D_RE_PTR),%VAL(D_RI3_PTR),%VAL(D_TG_PTR),
     '          %VAL(D_ZG_PTR),%VAL(ES_PTR),FEXT,GK,PG,%VAL(RE1_PTR),
     '          %VAL(RE2_PTR),%VAL(RG_PTR),SE,WG,XA,%VAL(XE_PTR),
     '          %VAL(XG_PTR),XP,YG,ZA,ZA1,ZAA,%VAL(ZE_PTR),
     '          %VAL(ZE1_PTR),%VAL(ZG_PTR),%VAL(ZG1_PTR),ZP,ZP1,ZPA,
     '          FIX,ERROR,*100)

C cpb 18/10/96 Free dynamically allocated arrays

              CALL FREE_MEMORY(LGE_PTR,ERROR,*100)
              CALL FREE_MEMORY(CG_PTR,ERROR,*100)
              CALL FREE_MEMORY(D_RE_PTR,ERROR,*100)
              CALL FREE_MEMORY(D_RI3_PTR,ERROR,*100)
              CALL FREE_MEMORY(D_TG_PTR,ERROR,*100)
              CALL FREE_MEMORY(D_ZG_PTR,ERROR,*100)
              CALL FREE_MEMORY(ES_PTR,ERROR,*100)
              CALL FREE_MEMORY(RG_PTR,ERROR,*100)
              CALL FREE_MEMORY(XE_PTR,ERROR,*100)
              CALL FREE_MEMORY(XG_PTR,ERROR,*100)
              CALL FREE_MEMORY(ZE_PTR,ERROR,*100)
              CALL FREE_MEMORY(ZE1_PTR,ERROR,*100)
              CALL FREE_MEMORY(ZG_PTR,ERROR,*100)
              CALL FREE_MEMORY(ZG1_PTR,ERROR,*100)
              IF(KTYP1D.NE.1) THEN
C              IF(.true..or.KTYP1D.NE.1) THEN
                CALL FREE_MEMORY(RE1_PTR,ERROR,*100)
                CALL FREE_MEMORY(RE2_PTR,ERROR,*100)
              ENDIF

              GO TO 102
C               This statement is designed to be skipped if no error
C               occur. However if a error occurs within a subroutine
C               the alternate return points to line 100 to set the flag
 100            CONTINUE
C$OMP           CRITICAL(ASSEMBLE5_1)
                ERROR_FLAG=.TRUE.
                WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
                CALL WRITES(IOER,OP_STRING,ERROR,*101)
                WRITE(OP_STRING,'(/'' >>An error occurred - '
     '            //'results may be unreliable!'')')
                CALL WRITES(IODI,OP_STRING,ERROR,*101)
 101            CONTINUE
C$OMP           END CRITICAL(ASSEMBLE5_1)
 102          CONTINUE
            ENDIF !.NOT.ERROR_FLAG
          ENDDO !noelem (ne)
C$OMP     END PARALLEL DO
        ENDDO !nonrlist (nr)
        CALL ASSERT(.NOT.ERROR_FLAG,'>>An error occurred during '
     '    //'element stiffness calculations',ERROR,*9999)

        IF(IWRIT4(nr_solve,nx).GE.1) THEN
          CALL CPU_TIMER(CPU_USER,TIME_STOP)
          ELAPSED_TIME=TIME_STOP(1)-TIME_START1(1)
          WRITE(OP_STRING,
     '      '(/'' For element stiffness calcs: CPU time of 1 thread:'','
     '      //'D11.4,'' s'')') ELAPSED_TIME
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          CALL REAL_TIMER(REAL_TOTAL,TIME_STOP)
          ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
          WRITE(OP_STRING,
     '      '( ''                              Elapsed (wall) time: '','
     '      //'D11.4,'' s'')') ELAPSED_TIME
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF


      ELSE IF(KTYP1A.EQ.2) THEN !calc elem stiff matrices in parallel
C       Form complete list of elements to assemble across all regions
C NEW 04/07/05 JHC removed parallel stuff. 
        CALL ASSERT(.FALSE.,'>>Parallel element stiffness computations '
     &    //'are no longer implemented',ERROR,*9999)

      ENDIF

      IF(ITYP5(nr,nx).EQ.5.AND.ITYP15(nr,nx).NE.0) THEN
C       Wavefront path analysis .AND. upwinding

        CALL CPU_TIMER(CPU_USER,TIME_START1)
        CALL REAL_TIMER(REAL_TOTAL,TIME_START2)

        DO nonrlist=1,NRLIST(0)
          nr=NRLIST(nonrlist)
          CALL ASSERT(NJ_LOC(NJL_GEOM,0,nr).EQ.3,
     '      '>>Only 3D upwinding is implemented',ERROR,*9999)

          ERROR_FLAG=.FALSE.
C new MPN 1Feb2000: OMP parallel proc directive
C old
CC$DOACROSS local (nf,noface,IDOXFT_PTR,MYMS_PTR,NSFE_PTR,
CC$&               NYNS_PTR,CG_PTR,FS_PTR,RDF1_PTR,RDF2_PTR,SM_PTR,
CC$&               SN_PTR,XDF_PTR,ZDF_PTR,ZE_PTR,ERROR)
CC$&        share (NPF,nr,nx,ERROR_FLAG)
C end old
C$OMP     PARALLEL DO
C$OMP&      PRIVATE(nf,noface,IDOXFT_PTR,MYMS_PTR,NSFE_PTR,
C$OMP&              NYNS_PTR,CG_PTR,FS_PTR,RDF1_PTR,RDF2_PTR,SM_PTR,
C$OMP&              SN_PTR,XDF_PTR,ZDF_PTR,ZE_PTR,ERROR),
C$OMP&      SHARED(MEM_INIT,NGM,NHM,NPF,nr,nx,ERROR_FLAG)
          DO noface=1,NFFACE(0,nr) !is main face loop
            nf=NFFACE(noface,nr)
            IF(.NOT.ERROR_FLAG) THEN
C     '        .AND.(ITYP15(nr,nx).EQ.3.EQV.NPF(5,nf).EQ.2)) THEN
CC             correct number of adjacent elements.

c             Dynamic allocation of parallel local arrays.
C             Intialise pointers so they are zero in the parallel loop.

              IDOXFT_PTR=0
              CALL ALLOCATE_MEMORY(NHM,1,INTTYPE,IDOXFT_PTR,
     '          MEM_INIT,ERROR,*200)
              MYMS_PTR=0
              CALL ALLOCATE_MEMORY((NSFM+1)*2*2*NHM,1,INTTYPE,MYMS_PTR,
     '          MEM_INIT,ERROR,*200)
              NSFE_PTR=0
              CALL ALLOCATE_MEMORY(NSM*2*2*NHM,1,INTTYPE,NSFE_PTR,
     '          MEM_INIT,ERROR,*200)
              NYNS_PTR=0
              CALL ALLOCATE_MEMORY((NSM+1)*2*2*NHM,1,INTTYPE,NYNS_PTR,
     '          MEM_INIT,ERROR,*200)
              CG_PTR=0
              CALL ALLOCATE_MEMORY(NMM*NGM,1,DPTYPE,CG_PTR,
     '          MEM_INIT,ERROR,*200)
              FS_PTR=0
              CALL ALLOCATE_MEMORY(NSFM*2*NSFM*2*2*NHM,1,DPTYPE,FS_PTR,
     '          MEM_INIT,ERROR,*200)
              SM_PTR=0
              CALL ALLOCATE_MEMORY(NSFM*2*2*NHM,1,DPTYPE,SM_PTR,
     '          MEM_INIT,ERROR,*200)
              SN_PTR=0
              CALL ALLOCATE_MEMORY(NSM*2*2*NHM,1,DPTYPE,SN_PTR,
     '          MEM_INIT,ERROR,*200)
              XDF_PTR=0
              CALL ALLOCATE_MEMORY(NSFM*2*2*NJM,1,DPTYPE,XDF_PTR,
     '          MEM_INIT,ERROR,*200)
              ZDF_PTR=0
              CALL ALLOCATE_MEMORY(NSFM*2*2*NHM,1,DPTYPE,ZDF_PTR,
     '          MEM_INIT,ERROR,*200)

C             These arrays are only necessary for finite differences
              RDF1_PTR=0
              RDF2_PTR=0
              IF(KTYP1D.NE.1) THEN
C              IF(.true..or.KTYP1D.NE.1) THEN
                CALL ALLOCATE_MEMORY(NSFM*2*NHM,1,DPTYPE,RDF1_PTR,
     '            MEM_INIT,ERROR,*200)
                CALL ALLOCATE_MEMORY(NSFM*2*NHM,1,DPTYPE,RDF2_PTR,
     '            MEM_INIT,ERROR,*200)
              ENDIF

C             Must use a dynam subroutine as dynamically allocated
C             arrays are used within this subroutine level.

              CALL ASSEMBLE5_FACE(IBT,IDO,%VAL(IDOXFT_PTR),INP,
     '          ISC_GK,ISR_GK,%VAL(MYMS_PTR),NBH,NBHF(1,1,nf),
     '          NBJ,NBJF(1,nf),nf,NHE,NKB,NKHE,NKJE,NNB,NNF,
     '          NPF(1,nf),NPNE,nr,NSB,%VAL(NSFE_PTR),
     '          NVHE,NVJE,nx,NYNP(1,1,1,1,0,1,nr),%VAL(NYNS_PTR),
     '          CE,%VAL(CG_PTR),CP,%VAL(FS_PTR),GK,PG,
     '          %VAL(RDF1_PTR),%VAL(RDF2_PTR),SE,
     '          %VAL(SM_PTR),%VAL(SN_PTR),%VAL(XDF_PTR),XP,
     '          YGF(1,1,nf),YP,%VAL(ZDF_PTR),ERROR,*200)

C             Free dynamically allocated arrays

              CALL FREE_MEMORY(IDOXFT_PTR,ERROR,*200)
              CALL FREE_MEMORY(MYMS_PTR,ERROR,*200)
              CALL FREE_MEMORY(NSFE_PTR,ERROR,*200)
              CALL FREE_MEMORY(NYNS_PTR,ERROR,*200)
              CALL FREE_MEMORY(CG_PTR,ERROR,*200)
              CALL FREE_MEMORY(FS_PTR,ERROR,*200)
              CALL FREE_MEMORY(SM_PTR,ERROR,*200)
              CALL FREE_MEMORY(SN_PTR,ERROR,*200)
              CALL FREE_MEMORY(XDF_PTR,ERROR,*200)
              CALL FREE_MEMORY(ZDF_PTR,ERROR,*200)
              IF(KTYP1D.NE.1) THEN
C              IF(.true..or.KTYP1D.NE.1) THEN
                CALL FREE_MEMORY(RDF1_PTR,ERROR,*200)
                CALL FREE_MEMORY(RDF2_PTR,ERROR,*200)
              ENDIF

              GO TO 202
C               This statement is designed to be skipped if no error
C               occurs.  However if a error occurs within a subroutine
C               the alternate return points to line 200 to set the flag
 200            CONTINUE
C$OMP           CRITICAL(ASSEMBLE5_2)
                ERROR_FLAG=.TRUE.
                WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
                CALL WRITES(IOER,OP_STRING,ERROR,*201)
 201            CONTINUE
C$OMP           END CRITICAL(ASSEMBLE5_2)
                IF(IDOXFT_PTR.NE.0)
     '            CALL FREE_MEMORY(IDOXFT_PTR,ERROR,*202)
                IF(MYMS_PTR.NE.0) CALL FREE_MEMORY(MYMS_PTR,ERROR,*202)
                IF(NSFE_PTR.NE.0) CALL FREE_MEMORY(NSFE_PTR,ERROR,*202)
                IF(NYNS_PTR.NE.0) CALL FREE_MEMORY(NYNS_PTR,ERROR,*202)
                IF(CG_PTR.NE.0) CALL FREE_MEMORY(CG_PTR,ERROR,*202)
                IF(FS_PTR.NE.0) CALL FREE_MEMORY(FS_PTR,ERROR,*202)
                IF(SM_PTR.NE.0) CALL FREE_MEMORY(SM_PTR,ERROR,*202)
                IF(SN_PTR.NE.0) CALL FREE_MEMORY(SN_PTR,ERROR,*202)
                IF(XDF_PTR.NE.0) CALL FREE_MEMORY(XDF_PTR,ERROR,*202)
                IF(ZDF_PTR.NE.0) CALL FREE_MEMORY(ZDF_PTR,ERROR,*202)
                IF(RDF1_PTR.NE.0) CALL FREE_MEMORY(RDF1_PTR,ERROR,*202)
                IF(RDF2_PTR.NE.0) CALL FREE_MEMORY(RDF2_PTR,ERROR,*202)
 202          CONTINUE
            ENDIF !.NOT.ERROR_FLAG
          ENDDO !noface (nf)
C$OMP     END PARALLEL DO
        ENDDO !nonrlist (nr)

        CALL ASSERT(.NOT.ERROR_FLAG,'>>An error occurred during '
     '    //'face stiffness calculations',ERROR,*9999)

        IF(IWRIT4(nr_solve,nx).GE.1) THEN
          CALL CPU_TIMER(CPU_USER,TIME_STOP)
          ELAPSED_TIME=TIME_STOP(1)-TIME_START1(1)
          WRITE(OP_STRING,
     '      '('' For face stiffness calcs:    CPU time of 1 thread:'','
     '      //'D11.4,'' s'')') ELAPSED_TIME
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          CALL REAL_TIMER(REAL_TOTAL,TIME_STOP)
          ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
          WRITE(OP_STRING,
     '      '( ''                              Elapsed (wall) time: '','
     '      //'D11.4,'' s'')') ELAPSED_TIME
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

      ENDIF

      IF(IWRIT4(nr_solve,nx).GE.4) THEN
        WRITE(OP_STRING,
     '    '(/'' Global stiffness matrix GK:'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' NYNR(0,1,1,nr_solve)='',I5,'
     '    //''', NYNR(0,2,1,nr_solve)='',I5)')
     '    NYNR(0,1,1,nr_solve),NYNR(0,2,1,nr_solve)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
c cpb 24/9/95 Adding generic stiffness matrix output
        CALL OPSTFMAT(NYNR(0,2,1,nr_solve),ISC_GK,ISR_GK,IOOP,
     '    NYT(1,1,nx),NYT(2,1,nx),NZT(1,nx),NYNR(0,1,1,nr_solve),
     '    KTYP24,GK,GK,'GK ','GR ',.FALSE.,.TRUE.,.FALSE.,ERROR,*9999)
      ENDIF

      CALL EXITS('ASSEMBLE5')
      RETURN

C NEW 04/07/05 JHC removed parallel stiffness matrix computations

C cpb 18/10/96 Free dynamically allocated arrays

 9999 IF(LGE_PTR.NE.0) CALL FREE_MEMORY(LGE_PTR,ERROR,*1113)
      IF(CG_PTR.NE.0) CALL FREE_MEMORY(CG_PTR,ERROR,*1113)
      IF(D_RE_PTR.NE.0) CALL FREE_MEMORY(D_RE_PTR,ERROR,*1113)
      IF(D_RI3_PTR.NE.0) CALL FREE_MEMORY(D_RI3_PTR,ERROR,*1113)
      IF(D_TG_PTR.NE.0) CALL FREE_MEMORY(D_TG_PTR,ERROR,*1113)
      IF(D_ZG_PTR.NE.0) CALL FREE_MEMORY(D_ZG_PTR,ERROR,*1113)
      IF(ES_PTR.NE.0) CALL FREE_MEMORY(ES_PTR,ERROR,*1113)
      IF(RE1_PTR.NE.0) CALL FREE_MEMORY(RE1_PTR,ERROR,*1113)
      IF(RE2_PTR.NE.0) CALL FREE_MEMORY(RE2_PTR,ERROR,*1113)
      IF(RG_PTR.NE.0) CALL FREE_MEMORY(RG_PTR,ERROR,*1113)
      IF(XE_PTR.NE.0) CALL FREE_MEMORY(XE_PTR,ERROR,*1113)
      IF(XG_PTR.NE.0) CALL FREE_MEMORY(XG_PTR,ERROR,*1113)
      IF(ZE_PTR.NE.0) CALL FREE_MEMORY(ZE_PTR,ERROR,*1113)
      IF(ZE1_PTR.NE.0) CALL FREE_MEMORY(ZE1_PTR,ERROR,*1113)
      IF(ZG_PTR.NE.0) CALL FREE_MEMORY(ZG_PTR,ERROR,*1113)
      IF(ZG1_PTR.NE.0) CALL FREE_MEMORY(ZG1_PTR,ERROR,*1113)

 1113 CALL ERRORS('ASSEMBLE5',ERROR)
      RET_ERROR=ERROR
      CALL EXITS('ASSEMBLE5')
      RETURN 1
      END


