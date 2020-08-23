      SUBROUTINE UPGAUS_EIK(IBT,IDO,NBH,NBHF,NBJ,NBJF,NEELEM,NFFACE,NHE,
     '  NKB,NKJE,NNF,NPF,NPNE,nr,NSB,NVJE,NW,nx,CE,CG,CGE,CP,PG,SE,WG,
     '  XA,XP,YG,YGF,OP_TIMING,RET_ERROR,*)

C#### Subroutine: UPGAUS_EIK
C###  Description:
C###    Sets up the gauss point arrays from geometry and material
C###    parameters for solution of an eikonal equation.

C#### Variable: YGF(niygf,ng,nf)
C###  Type: REAL*8
C###  Set_up: UPGAUS_EIK
C###  Description:
C###    YG(niyg,ng,ne) is the face Gauss point array.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'error0.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'mxch.inc'
      INCLUDE 'time02.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NEELEM(0:NE_R_M),NFFACE(0:NF_R_M),NHE(NEM),
     '  NKB(2,2,2,NNM,NBFM),NKJE(NKM,NNM,NJM,NEM),NNF(0:17,6,NBFM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),nr,NSB(NKM,NNM,NBFM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),CGE(NMM,NGM,NEM),CP(NMM,NPM),
     '  PG(NSM,NUM,NGM,NBM),SE(NSM,NBFM,NEM),WG(NGM,NBM),
     '  XA(NAM,NJM,NEM),XP(NKM,NVM,NJM,NPM),
     '  YG(NIYGM,NGM,NEM),YGF(NIYGFM,NGFM,NFM)
      LOGICAL OP_TIMING
      CHARACTER RET_ERROR*(*)
!     Local Variables
      INTEGER nb,ne,NEA(2),NEAFST,NEAXT,nf,NF_EA(2),ng,NHF,NIEF(0:2),
     '  NITB,nj,njj,noelem,noface
      INTEGER*4 XDF_PTR,XE_PTR,XG_PTR,XGF_PTR
      REAL ELAPSED_TIME,TIME_START1(1),
     '  TIME_START2(1),TIME_STOP(1)
      CHARACTER ERROR*(ERRSTRLEN)
      LOGICAL ERROR_FLAG,FLUX,INTEGRATE,MODDIFF

      CALL ENTERS('UPGAUS_EIK',*9999)

      IF(OP_TIMING) THEN
        CALL CPU_TIMER(CPU_USER,TIME_START1)
        CALL REAL_TIMER(REAL_TOTAL,TIME_START2)
      ENDIF
C      DO il=1,ILT(1,nr,nx)
C        CALL ASSERT(ABS(ILP(il,1,nr,nx)).EQ.1,
C     '    '>>Material parameters must be constant',ERROR,*9999)
C      ENDDO !il
      ERROR_FLAG=.FALSE.

C$OMP PARALLEL DO
C$OMP&PRIVATE (MODDIFF,nb,ne,nf,NITB,nj,njj,noelem,XE_PTR,XG_PTR,ERROR)
C$OMP&    SHARED (MEM_INIT,nr,nx,YG,ERROR_FLAG)
      DO noelem=1,NEELEM(0) !is main element loop
        ne=NEELEM(noelem)
        IF(NW(ne,1).GE.0.AND..NOT.ERROR_FLAG) THEN

C         Dynamic allocation of parallel local arrays.
C         Intialise pointers so they are zero in the parallel loop.
          XE_PTR=0
          CALL ALLOCATE_MEMORY(NSM*NJM,1,DPTYPE,XE_PTR,
     '      MEM_INIT,ERROR,*100)
          XG_PTR=0
          CALL ALLOCATE_MEMORY(NJM*NUM,1,DPTYPE,XG_PTR,
     '      MEM_INIT,ERROR,*100)

C         Transfer geometric variables to local element array.
          CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '      NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '      SE(1,1,ne),XA(1,1,ne),%VAL(XE_PTR),XP,ERROR,*100)

          NITB=NIT(NBJ(1,ne))
          nb=NBH(NH_LOC(1,nx),1,ne)
C         The last field may contain a diffusion modifier that
C         approaches zero at singular initiation points.
          njj=NJ_LOC(NJL_FIEL,0,nr)
          MODDIFF=njj.EQ.3
          IF(MODDIFF) THEN
            nj=NJ_LOC(NJL_FIEL,njj,nr)
            MODDIFF=NBJ(nj,ne).NE.0
          ENDIF

C         Material parameters
          CALL CPCG(1,nb,NPNE(1,1,ne),nr,nx,
     '      CE(1,ne),CG,CGE(1,1,ne),CP,PG,ERROR,*100)

C         Loop over Gauss points
          DO ng=1,NGT(nb)
C           Evalulate geometric variables at Gauss point.
            CALL XEXG(NBJ(1,ne),ng,nr,PG,%VAL(XE_PTR),%VAL(XG_PTR),
     '        ERROR,*100)
C           Calculate necessary terms from geometry and material
C           params.
            CALL XGYG(NITB,nr,nx,CG(1,ng),WG(ng,nb),
     '        %VAL(XG_PTR),YG(1,ng,ne),MODDIFF,ERROR,*100)
          ENDDO !ng

C         Free dynamically allocated arrays
          CALL FREE_MEMORY(XE_PTR,ERROR,*100)
          CALL FREE_MEMORY(XG_PTR,ERROR,*100)

          GO TO 102
C           This statement is designed to be skipped if no error
C           occur. However if a error occurs within a subroutine
C           the alternate return points to line 100 to set the flag
 100        CONTINUE
C$OMP       CRITICAL(UPGAUS_EIK_CR)
            ERROR_FLAG=.TRUE.
            WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
            CALL WRITES(IOER,OP_STRING,ERROR,*101)
 101        CONTINUE
C$OMP       END CRITICAL(UPGAUS_EIK_CR)
            IF(XE_PTR.NE.0) CALL FREE_MEMORY(XE_PTR,ERROR,*102)
            IF(XG_PTR.NE.0) CALL FREE_MEMORY(XG_PTR,ERROR,*102)
 102      CONTINUE

        ENDIF !ERROR_FLAG
      ENDDO !noelem (ne)

      CALL ASSERT(.NOT.ERROR_FLAG,'>>An error occurred during '
     '  //'element Gauss array calculations',ERROR,*9999)

      IF(OP_TIMING) THEN
        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START1(1)
        WRITE(OP_STRING,
     '    '(/'' For element Gauss calcs: CPU time of 1 thread:'','
     '    //'D11.4,'' s'')') ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        CALL REAL_TIMER(REAL_TOTAL,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
        WRITE(OP_STRING,
     '    '( ''                          Elapsed (wall) time: '','
     '      //'D11.4,'' s'')') ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      IF(ITYP15(nr,nx).NE.0) THEN !upwinding

        IF(OP_TIMING) THEN
          CALL CPU_TIMER(CPU_USER,TIME_START1)
          CALL REAL_TIMER(REAL_TOTAL,TIME_START2)
        ENDIF
C        DO il=1,ILT(1,nr,nx)
C          CALL ASSERT(ABS(ILP(il,1,nr,nx)).EQ.1,
C     '      '>>Material parameters must be constant',ERROR,*9999)
C        ENDDO !il
        CALL ASSERT(NJ_LOC(NJL_GEOM,0,nr).EQ.3,
     '    '>>Only 3D upwinding is implemented',ERROR,*9999)
        ERROR_FLAG=.FALSE.

C$OMP   PARALLEL DO
C$OMP&    PRIVATE (FLUX,INTEGRATE,nb,ne,nf,NEAFST,NEAXT,NHF,NITB,noelem,
C$OMP&    XDF_PTR,XGF_PTR,ERROR)
C$OMP&  SHARED (MEM_INIT,nr,nx,YGF,ERROR_FLAG)
        DO noface=1,NFFACE(0) !is main element loop
          nf=NFFACE(noface)
          IF(.NOT.ERROR_FLAG) THEN
C     '      .AND.(ITYP15(nr,nx).NE.3.OR.NPF(5,nf).EQ.2)) THEN
CC           correct number of adjacent elements.

C           Find connectivity information.
            CALL FACE_INT_PREP(NEA,NEAXT,NF_EA,NHE,NHF,NIEF,
     '        NPF(1,nf),nr,nx,INTEGRATE,ERROR,*200)

            IF(INTEGRATE) THEN !false if integral is known to be zero

C             Dynamic allocation of parallel local arrays.
C             Intialise pointers so they are zero in the parallel loop.
              XDF_PTR=0
              CALL ALLOCATE_MEMORY(NSFM*2*2*NJM,1,DPTYPE,XDF_PTR,
     '          MEM_INIT,ERROR,*200)
              XGF_PTR=0
              CALL ALLOCATE_MEMORY(NJM*NUM*2,1,DPTYPE,XGF_PTR,
     '          MEM_INIT,ERROR,*200)

C             Transfer geometric variables to local face array.
              CALL XPXF(IBT,IDO,NBJ,NBJF(1,nf),NEA,NEAXT,NF_EA,NIEF,
     '          NKB,NKJE,NNF,NPNE,nr,NSB,NVJE,SE,%VAL(XDF_PTR),XP,
     '          ERROR,*200)

              nb=NBHF(NH_LOC(1,nx),1,nf)
              NITB=NIT(nb)
              ne=NEA(1)

C!!!    should check on order of derivatives
              FLUX=ITYP15(nr,nx).NE.3.OR..TRUE.
              IF(FLUX) THEN
                NEAFST=NEAXT
              ELSE
                NEAFST=1 !derivs are the same in each neighbouring element
              ENDIF

C             Material parameters
              CALL CPCGF(NBJF(1,nf),NEA,NF_EA,NNF,NPNE,nr,nx,CE,CG,CP,
     '          PG,ERROR,*200)

C             Loop over Gauss points
              DO ng=1,NGT(nb)
C               Evalulate geometric variables at Gauss point.
                CALL XFXG(NBJF(1,nf),NEAFST,NIEF,ng,nr,PG,%VAL(XDF_PTR),
     '            %VAL(XGF_PTR),ERROR,*200)
C               Calculate necessary terms from geometry and material params.
                CALL XGYGF(NEAXT,NIEF,NITB,nr,nx,CG(1,ng),WG(ng,nb),
     '            %VAL(XGF_PTR),YGF(1,ng,nf),ERROR,*200)
              ENDDO !ng

C             Free dynamically allocated arrays
              CALL FREE_MEMORY(XDF_PTR,ERROR,*200)
              CALL FREE_MEMORY(XGF_PTR,ERROR,*200)

            ENDIF !INTEGRATE

            GO TO 202
C             This statement is designed to be skipped if no error
C             occurs.  However if a error occurs within a subroutine
C             the alternate return points to line 200 to set the flag
 200          CONTINUE
C$OMP         CRITICAL(UPGAUS_EIK_CR)
              ERROR_FLAG=.TRUE.
              WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
              CALL WRITES(IOER,OP_STRING,ERROR,*201)
 201          CONTINUE
C$OMP         END CRITICAL(UPGAUS_EIK_CR)
              IF(XG_PTR.NE.0) CALL FREE_MEMORY(XGF_PTR,ERROR,*202)
              IF(XDF_PTR.NE.0) CALL FREE_MEMORY(XDF_PTR,ERROR,*202)
 202        CONTINUE

          ENDIF !.NOT.ERROR_FLAG
        ENDDO !noelem

        CALL ASSERT(.NOT.ERROR_FLAG,'>>An error occurred during '
     '    //'face Gauss array calculations',ERROR,*9999)

        IF(OP_TIMING) THEN
          CALL CPU_TIMER(CPU_USER,TIME_STOP)
          ELAPSED_TIME=TIME_STOP(1)-TIME_START1(1)
          WRITE(OP_STRING,
     '      '('' For face Gauss calcs:    CPU time of 1 thread:'','
     '      //'D11.4,'' s'')') ELAPSED_TIME
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CALL REAL_TIMER(REAL_TOTAL,TIME_STOP)
          ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
          WRITE(OP_STRING,
     '      '(''                          Elapsed (wall) time: '','
     '      //'D11.4,'' s'')') ELAPSED_TIME
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

      ENDIF !upwinding

      CALL EXITS('UPGAUS_EIK')
      RETURN

 9999 CALL ERRORS('UPGAUS_EIK',ERROR)
      RET_ERROR=ERROR
      CALL EXITS('UPGAUS_EIK')
      RETURN 1
      END


