      SUBROUTINE ASSEMBLE3(IBT,IDO,INP,ISC_GD,ISC_GK,ISC_GM,ISR_GD,
     '  ISR_GK,ISR_GM,LGE,NBH,NBJ,NEELEM,NHE,NKHE,NKJE,NORD,NPF,
     '  NPNE,nr,NVHE,NVJE,NW,nx,NYNE,NYNP,NYNR,
     '  CE,CG,CGE,CP,ED,EM,ER,ES,GD,GK,GM,GR,PG,RG,SE,WG,XA,XE,XG,
     '  XP,YG,ZA,ZE,ZG,ZP,DYNAM1,DYNAM2,UPDATE_MATRIX,ERROR,*)

C#### Subroutine: ASSEMBLE3
C###  Description:
C###    ASSEMBLE3 assembles the global unreduced matrices GK, GD, GM,
C###    etc. for time dependent FEM problems.

C**** If DYNAM1 is true the problem will use the GD matrix
C**** If DYNAM2 is true the problem will use the GM matrix
C**** FIRST indicates whether or not the routine is being called for
C****   the first time. If FIRST is true the (temporary) matrix sparsity
C****   patterns are calculated. If the routine is not being called
C****   for the first time FIRST should be false.
C**** If UPDATE_MATRIX is true the GK (and GD, GM matrices) will be
C****   calculated otherwise just the R.H.S. vector (GR) will be
C****   calculated.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'coup00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'time02.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISC_GD(NISC_GDM),ISC_GK(NISC_GKM),ISC_GM(NISC_GMM),
     '  ISR_GD(NISR_GDM),ISR_GK(NISR_GKM),ISR_GM(NISR_GMM),
     '  LGE(NHM*NSM,NRCM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NHE(NEM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NORD(5,NE_R_M),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),nr,NVHE(NNM,NBFM,NHM,NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),NYNR(0:NY_R_M,0:NRCM,NCM)
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),CGE(NMM,NGM,NEM),CP(NMM,NPM),
     '  ED(NHM*NSM,NHM*NSM),EM(NHM*NSM,NHM*NSM),ER(NHM*NSM),
     '  ES(NHM*NSM,NHM*NSM),GD(NZ_GD_M),GK(NZ_GK_M),GM(NZ_GM_M),
     '  GR(NYROWM),PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM,NEM),
     '  WG(NGM,NBM),XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),
     '  XP(NKM,NVM,NJM,NPM),YG(NIYGM,NGM,NEM),ZA(NAM,NHM,NCM,NEM),
     '  ZE(NSM,NHM),ZG(NHM,NUM),ZP(NKM,NVM,NHM,NPM,NCM)
      LOGICAL DYNAM1,DYNAM2,UPDATE_MATRIX
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ne,nhs1,nhs2,NHST(2),noelem,no_nynr1,NO_TIMES,ny1,ny2,nz
      REAL AVETIME,ELAPSED_TIME,TIME_START1(1),TIME_START2(1),
     '  TIME_START3(1),TIME_STOP(1)

      CALL ENTERS('ASSEMBLE3',*9999)

      CALL CPU_TIMER(CPU_USER,TIME_START1)
      CALL CPU_TIMER(CPU_USER,TIME_START2)
      IF(UPDATE_MATRIX) THEN
C***  Initialise matrices for this region
        CALL INIT_SPARSE_MATRIX(NYNR(0,2,1),NYNR(0,2,1),ISC_GK,ISR_GK,
     '    0,0,1,NYT(1,1,nx),NZ_GK_M,NZT(1,nx),NYNR(0,1,1),NYNR(0,1,1),
     '    KTYP24,GK,ERROR,*9999)
        IF(DYNAM1) THEN
C cpb 31/3/96 At the moment the sparsity for the GD and GM matrices
C are not setup anywhere. They should be set up in ipequa but for
C things like cubic time stepping you might need the GM matrix for
C an equation that does not use it. Hence will set up the sparsity
C patterns for the GD and GM matrices here by copying the sparsity
C pattern from the GK matrix. This will have to be looked at in the
C future
          NYT(1,3,nx)=NYT(1,1,nx)
          NYT(2,3,nx)=NYT(2,1,nx)
          NZT(3,nx)=NZT(1,nx)
          CALL ASSERT(NISC_GDM.GE.NISC_GKM,'>>Increase NISC_GDM',
     '      ERROR,*9999)
          CALL ASSERT(NISR_GDM.GE.NISR_GKM,'>>Increase NISR_GDM',
     '      ERROR,*9999)
          DO nz=1,NISC_GKM
            ISC_GD(nz)=ISC_GK(nz)
          ENDDO !nz
          DO nz=1,NISR_GKM
            ISR_GD(nz)=ISR_GK(nz)
          ENDDO !nz
          CALL INIT_SPARSE_MATRIX(NYNR(0,2,1),NYNR(0,2,1),ISC_GD,ISR_GD,
     '      0,0,1,NYT(1,3,nx),NZ_GD_M,NZT(3,nx),NYNR(0,1,1),NYNR(0,1,1),
     '      KTYP24,GD,ERROR,*9999)
        ENDIF
        IF(DYNAM2) THEN
          NYT(1,4,nx)=NYT(1,1,nx)
          NYT(2,4,nx)=NYT(2,1,nx)
          NZT(4,nx)=NZT(1,nx)
          CALL ASSERT(NISC_GMM.GE.NISC_GKM,'>>Increase NISC_GMM',
     '      ERROR,*9999)
          CALL ASSERT(NISR_GMM.GE.NISR_GKM,'>>Increase NISR_GMM',
     '      ERROR,*9999)
          DO nz=1,NISC_GKM
            ISC_GM(nz)=ISC_GK(nz)
          ENDDO !nz
          DO nz=1,NISR_GKM
            ISR_GM(nz)=ISR_GK(nz)
          ENDDO !nz
          CALL INIT_SPARSE_MATRIX(NYNR(0,2,1),NYNR(0,2,1),ISC_GM,ISR_GM,
     '      0,0,1,NYT(1,4,nx),NZ_GD_M,NZT(4,nx),NYNR(0,1,1),NYNR(0,1,1),
     '      KTYP24,GM,ERROR,*9999)
        ENDIF
      ENDIF
      DO nz=1,NZT(1,nx)
        GD(nz)=0.0d0
        GK(nz)=0.0d0
      ENDDO
      DO no_nynr1=1,NYNR(0,1,1)
        ny1=NYNR(no_nynr1,1,1)
        GR(ny1)=0.0d0
      ENDDO !no_nynr (ny1)

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
      IF(IWRIT4(nr,nx).GE.1) THEN
        WRITE(OP_STRING,'(/'' CPU time for setup and '
     '    //'initialisation: '',D11.4,'' s'')') ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

C***  Find element stiffness matrices
      CALL CPU_TIMER(CPU_USER,TIME_START2)
      AVETIME=0.0
      NO_TIMES=0
      DO noelem=1,NEELEM(0,nr)
        CALL CPU_TIMER(CPU_USER,TIME_START3)
        ne=NEELEM(noelem,nr)
        IF(NW(ne,1).GE.1.AND.NW(ne,1).LT.20) THEN
          CALL MELGE(LGE,NBH(1,1,ne),1,ne,NHE(ne),NHST,
     '      NPNE(1,1,ne),nr,NVHE(1,1,1,ne),nx,NYNE,NYNP,ERROR,*9999)
          IF(IWRIT4(nr,nx).GE.5) THEN
            FORMAT='(/'' Element'',I5,'', Number of variables: '','
     '        //'''NHST(1)='',I3,'', NHST(2)='',I3)'
            WRITE(OP_STRING,FORMAT) ne,NHST(1),NHST(2)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            FORMAT='('' LGE(1..,1): '',8(1X,I5),/:(13X,8(1X,I5)))'
            WRITE(OP_STRING,FORMAT) (LGE(nhs1,1),nhs1=1,NHST(1))
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            FORMAT='('' LGE(1..,2): '',8(1X,I5),/:(13X,8(x,I5)))'
            WRITE(OP_STRING,FORMAT) (LGE(nhs2,2),nhs2=1,NHST(2))
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF

C*** Initialize element arrays

          IF(UPDATE_MATRIX) THEN
            DO nhs1=1,NHST(1)
              DO nhs2=1,NHST(2)
                ES(nhs1,nhs2)=0.0d0
              ENDDO !nhs2
              IF(DYNAM1) THEN
                DO nhs2=1,NHST(2)
                  ED(nhs1,nhs2)=0.0d0
                ENDDO !nhs2
              ENDIF
              IF(DYNAM2) THEN
                DO nhs2=1,NHST(2)
                  EM(nhs1,nhs2)=0.0d0
                ENDDO !nhs2
              ENDIF
            ENDDO !nhs1
          ENDIF
          DO nhs1=1,NHST(1)
            ER(nhs1)=0.0d0
          ENDDO !nhs1

          IF(ITYP1(nr,nx).EQ.3) THEN !partial differential equation
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '        SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
            CALL XPES30(IBT,IDO,INP,NBH(1,1,ne),NBJ(1,ne),ne,
     '        NHE(ne),NORD(5,ne),NPNE(1,1,ne),nr,nx,
     '        CE(1,ne),CG,CGE(1,1,ne),CP,ED,EM,ER,ES,PG,RG,SE(1,1,ne),
     '        WG,XE,XG,YG(1,1,ne),ZE,ZG,
     '        UPDATE_MATRIX,.TRUE.,ERROR,*9999)
          ELSE IF(ITYP1(nr,nx).EQ.4) THEN !linear elasticity
            CALL XPES40(NBH(1,1,ne),NBJ(1,ne),
     '        NHE(ne),NKJE(1,1,1,ne),NPF,NPNE(1,1,ne),
     '        nr,NVJE(1,1,1,ne),NW(ne,1),nx,
     '        CE(1,ne),CG,CGE(1,1,ne),CP,ED,EM,ER,ES,PG,RG,SE(1,1,ne),
     '        WG,XA(1,1,ne),XE,XG,XP,YG(1,1,ne),
     '        UPDATE_MATRIX,ERROR,*9999)
          ENDIF
          IF(IWRIT4(nr,nx).GE.5) THEN
            IF(UPDATE_MATRIX) THEN
              WRITE(OP_STRING,'(/'' Element load vector ER & '
     '          //'stiffness matrix ES:'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ELSE
              WRITE(OP_STRING,'(/'' Element load vector ER:'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            CALL OPESTFMAT(NHST,IOOP,ES,ER,'ES','ER',UPDATE_MATRIX,
     '        .TRUE.,ERROR,*9999)
            IF(UPDATE_MATRIX) THEN
              IF(DYNAM1) THEN
                WRITE(OP_STRING,'(/'' Element damping matrix ED:'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                CALL OPESTFMAT(NHST,IOOP,ED,ER,'ED','ER',.TRUE.,
     '            .FALSE.,ERROR,*9999)
              ENDIF
              IF(DYNAM2) THEN
                WRITE(OP_STRING,'(/'' Element mass matrix EM:'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                CALL OPESTFMAT(NHST,IOOP,EM,ER,'EM','ER',.TRUE.,
     '            .FALSE.,ERROR,*9999)
              ENDIF
            ENDIF
          ENDIF

C*** Assemble element matrices into global matrices

          IF(UPDATE_MATRIX) THEN
            DO nhs1=1,NHST(1)
              ny1=IABS(LGE(nhs1,1))
              DO nhs2=1,NHST(2)
                ny2=IABS(LGE(nhs2,2))
                CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,NZ_GK_M,
     '            NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
                GK(nz)=GK(nz)+ES(nhs1,nhs2)
              ENDDO !nhs2
              IF(DYNAM1) THEN
                DO nhs2=1,NHST(2)
                  ny2=IABS(LGE(nhs2,2))
                  CALL SPARSE(ny1,ny2,NYT(1,3,nx),nz,NZ_GD_M,
     '              NZT(3,nx),ISC_GD,ISR_GD,KTYP24,ERROR,*9999)
                  GD(nz)=GD(nz)+ED(nhs1,nhs2)
                ENDDO !nhs2
              ENDIF
              IF(DYNAM2) THEN
                DO nhs2=1,NHST(2)
                  ny2=IABS(LGE(nhs2,2))
                  CALL SPARSE(ny1,ny2,NYT(1,4,nx),nz,NZ_GM_M,
     '              NZT(4,nx),ISC_GM,ISR_GM,KTYP24,ERROR,*9999)
                  GM(nz)=GM(nz)+EM(nhs1,nhs2)
                ENDDO !nhs2
              ENDIF
            ENDDO !nhs1
          ENDIF
          DO nhs1=1,NHST(1)
            ny1=IABS(LGE(nhs1,1))
            GR(ny1)=GR(ny1)+ER(nhs1)
          ENDDO !nhs1
        ELSE
          ERROR='>>Invalid NW number'
          GOTO 9999
        ENDIF
        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START3(1)
        AVETIME=AVETIME+ELAPSED_TIME
        NO_TIMES=NO_TIMES+1
C MHT temporary...remove comments if you want
C       IF(IWRIT4(nr,nx).GE.1) THEN
C         WRITE(OP_STRING,'(/'' CPU time for element '',I5,'
C    '      //''' assembly: '',D11.4,'' s'')') ne,ELAPSED_TIME
C         CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C       ENDIF
      ENDDO !noelem (ne)

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
      IF(IWRIT4(nr,nx).GE.1) THEN
        WRITE(OP_STRING,'(/'' CPU time for stiffness matrix assembly '
     '    //': '',D11.4,'' s'')') ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Average CPU time per element for '
     '    //'assembly:'',D11.4,'' s'')') AVETIME/NO_TIMES
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      IF(IWRIT4(nr,nx).GE.4) THEN
        CALL CPU_TIMER(CPU_USER,TIME_START2)
        IF(UPDATE_MATRIX) THEN
          WRITE(OP_STRING,'(/'' Global load vector GR & stiffness '
     '      //'matrix GK:'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NYNR(0,1,1)='',I5,'', NYNR(0,2,1)='','
     '      //'I5)') NYNR(0,1,1),NYNR(0,2,1)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ELSE
          WRITE(OP_STRING,'(/'' Global load vector GR:'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NYNR(0,1,1)='',I5)') NYNR(0,1,1)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        CALL OPSTFMAT(NYNR(0,2,1),ISC_GK,ISR_GK,IOOP,NYT(1,1,nx),
     '    NYT(2,1,nx),NZT(1,nx),NYNR(0,1,1),KTYP24,GK,GR,'GK ',
     '    'GR ',.FALSE.,UPDATE_MATRIX,.TRUE.,ERROR,*9999)
        IF(UPDATE_MATRIX) THEN
          IF(DYNAM1) THEN
            WRITE(OP_STRING,'(/'' Global damping matrix GD:'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' NYNR(0,1,1)='',I5,'', NYNR(0,2,1)='','
     '        //'I5)') NYNR(0,1,1),NYNR(0,2,1)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            CALL OPSTFMAT(NYNR(0,2,1),ISC_GD,ISR_GD,IOOP,NYT(1,3,nx),
     '        NYT(2,3,nx),NZT(3,nx),NYNR(0,1,1),KTYP24,GD,GR,'GD ',
     '        'GR ',.FALSE.,.TRUE.,.FALSE.,ERROR,*9999)
          ENDIF
          IF(DYNAM2) THEN
            WRITE(OP_STRING,'(/'' Global mass matrix GM:'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' NYNR(0,1,1)='',I5,'', NYNR(0,2,1)='','
     '        //'I5)') NYNR(0,1,1),NYNR(0,2,1)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            CALL OPSTFMAT(NYNR(0,2,1),ISC_GM,ISR_GM,IOOP,NYT(1,4,nx),
     '        NYT(2,4,nx),NZT(4,nx),NYNR(0,1,1),KTYP24,GM,GR,'GM ',
     '        'GR ',.FALSE.,.TRUE.,.FALSE.,ERROR,*9999)
          ENDIF
        ENDIF

        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
        WRITE(OP_STRING,'(/'' CPU time for stiffness matrices '
     '    //'output: '',D11.4,'' s'')') ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START1(1)
      IF(IWRIT4(nr,nx).GE.1) THEN
        WRITE(OP_STRING,'(/'' Total CPU time for assembly: '','
     '    //'D11.4,'' s'')') ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('ASSEMBLE3')
      RETURN
 9999 CALL ERRORS('ASSEMBLE3',ERROR)
      CALL EXITS('ASSEMBLE3')
      RETURN 1
      END


