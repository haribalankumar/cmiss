      SUBROUTINE ASSEMBLE1_DYNAM(IBT,IDO,INP,ISC_GK,ISC_GKK,ISC_GQ,
     '  ISR_GK,ISR_GKK,ISR_GQ,LGE,NBH,NBJ,ne,NEELEM,NHE,NHP,NKJE,
     '  NONY,NORD,NP_INTERFACE,NPF,NPNE,NPNY,nr,nr_gkk,NRE,NVHE,NVJE,
     '  NW,nx,NYNE,NYNP,CE,CG,CGE,
     '  CONY,CP,CURVCORRECT,ED,EM,ER,ES,GK,GKK,
     '  GQ,GR,PG,RG,SE,WG,XA,XE,XG,XG1,XP,YG,ZE,ZG,GQ_ASSEM,
     '  TIME_START3,UPDATE_MATRIX,UPDATE_VECTOR,ERROR,*)

C#### Subroutine: ASSEMBLE1_DYNAM
C###  Description:
C###    ASSEMBLE1_DYNAM is the dynamic subroutine associted with
C###    ASSEMBLE1.


      IMPLICIT NONE
      INCLUDE 'bem000.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'mxch.inc'
      INCLUDE 'solv00.cmn'
      INCLUDE 'time02.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),ISC_GK(NISC_GKM),ISC_GKK(NISC_GKKM),
     '  ISC_GQ(NISC_GQM),ISR_GK(NISR_GKM),ISR_GKK(NISR_GKKM),
     '  ISR_GQ(NISR_GQM),LGE(NHM*NSM,NRCM),NBH(NHM,NCM,NEM),
     '  NBJ(NJM,NEM),ne,NEELEM(0:NE_R_M,0:NRM),NHE(NEM),NHP(NPM),
     '  NKJE(NKM,NNM,NJM,NEM),NONY(0:NOYM,NYM,NRCM),NORD(5,NE_R_M),
     '  NP_INTERFACE(0:NPM,0:3),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NPNY(0:6,NYM,0:NRCM),nr,nr_gkk,NRE(NEM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),
     '  NW(NEM,3),nx,NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL TIME_START3(1)
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),CGE(NMM,NGM,NEM),
     '  CONY(0:NOYM,NYM,NRCM),
     '  CP(NMM,NPM),CURVCORRECT(2,2,NNM,NEM),ED(NHM*NSM,NHM*NSM),
     '  EM(NHM*NSM,NHM*NSM),ER(NHM*NSM),ES(NHM*NSM,NHM*NSM),
     '  GK(NZ_GK_M),GKK(NZ_GKK_M),GQ(NZ_GQ_M),GR(NYROWM),
     '  PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM,NEM),
     '  WG(NGM,NBM),XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XG(NJM,NUM),XG1(NJM,NUM,NGM),XP(NKM,NVM,NJM,NPM),
     '  YG(NIYGM,NGM,NEM),
     '  ZE(NSM,NHM),ZG(NHM,NUM)
      LOGICAL GQ_ASSEM,UPDATE_MATRIX,UPDATE_VECTOR
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER GETNYR,NHST(2),nhs,nhs1,nhs2,no1,no2,noy1,noy2,ny1,ny2,
     '  ny3,ny4,nz,nzz
      REAL ELAPSED_TIME,TIME_STOP(1)
      REAL*8 co1,co2
      CHARACTER FORMAT*500

      CALL ENTERS('ASSEMBLE1_DYNAM',*9999)

      IF(NW(ne,1).GT.0) THEN
        CALL MELGE(LGE,NBH(1,1,ne),1,ne,NHE(ne),NHST,
     '    NPNE(1,1,ne),nr,NVHE(1,1,1,ne),nx,NYNE,NYNP,ERROR,*9999)
        IF(IWRIT4(nr,nx).GE.5) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(ASSEMBLE1_DYNAM_1)
          FORMAT='(/'' Element'',I5,'', Number of variables: '','
     '      //'''NHST(1)='',I3,'', NHST(2)='',I3)'
          WRITE(OP_STRING,FORMAT) ne,NHST(1),NHST(2)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          FORMAT='('' LGE(1..,1): '',8(1X,I5),/:(13X,8(1X,I5)))'
          WRITE(OP_STRING,FORMAT) (LGE(nhs,1),nhs=1,NHST(1))
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          FORMAT='('' LGE(1..,2): '',8(1X,I5),/:(13X,8(x,I5)))'
          WRITE(OP_STRING,FORMAT) (LGE(nhs,2),nhs=1,NHST(2))
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(ASSEMBLE1_DYNAM_1)
        ENDIF

        DO nhs1=1,NHST(1)
          ER(nhs1)=0.0d0
          DO nhs2=1,NHST(2)
            ES(nhs1,nhs2)=0.0d0
          ENDDO !nhs2
        ENDDO !nhs1
        IF(NW(ne,1).LE.20) THEN !finite element
          IF(ITYP1(nr,nx).EQ.3) THEN !partial differential equation
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '        SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
            IF(ITYP2(nr,nx).EQ.7) THEN !Maxwells
              CALL XPES37(NBH(1,1,ne),NBJ(1,ne),ne,NHE(ne),NPNE(1,1,ne),
     &          nr,nx,CE(1,ne),CG,CGE(1,1,ne),CP,ED,EM,ER,ES,PG,RG,WG,
     &          XE,XG,UPDATE_MATRIX,UPDATE_VECTOR,ERROR,*9999)
            ELSE
              CALL XPES30(IBT,IDO,INP,NBH(1,1,ne),NBJ(1,ne),ne,
     '          NHE(ne),NORD(5,ne),NPNE(1,1,ne),nr,nx,
     '          CE(1,ne),CG,CGE(1,1,ne),CP,ED,EM,ER,ES,PG,RG,SE(1,1,ne),
     '          WG,XE,XG,YG(1,1,ne),ZE,ZG,
     '          UPDATE_MATRIX,UPDATE_VECTOR,ERROR,*9999)
            ENDIF
          ELSE IF(ITYP1(nr,nx).EQ.4) THEN !linear elasticity
            CALL XPES40(NBH(1,1,ne),NBJ(1,ne),
     '        NHE(ne),NKJE(1,1,1,ne),NPF,NPNE(1,1,ne),
     '        nr,NVJE(1,1,1,ne),NW(ne,1),nx,
     '        CE(1,ne),CG,CGE(1,1,ne),CP,ED,EM,ER,ES,PG,RG,SE(1,1,ne),
     '        WG,XA(1,1,ne),XE,XG,XP,YG(1,1,ne),UPDATE_MATRIX,
     '        ERROR,*9999)
          ENDIF
        ELSE !boundary element
          ERROR='>>Boundary elements should use ASSEMBLE2'
          GO TO 9999
        ENDIF
        IF(IWRIT4(nr,nx).GE.5) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(ASSEMBLE1_DYNAM_2)
          WRITE(OP_STRING,
     '      '(/'' Element load vector ER & stiffness matrix ES:'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
c cpb 24/9/95 Use generic element stiffness matrix output
          CALL OPESTFMAT(NHST,IOOP,ES,ER,'ES','ER',.TRUE.,.TRUE.,
     '      ERROR,*9999)
CC$OMP END CRITICAL(ASSEMBLE1_DYNAM_2)
        ENDIF

C*** Assemble element stiffness matrix into global system.


        IF(CALC_GLOBAL(nr,nx)) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Atomic directive is probably better anyway.
CC$OMP CRITICAL(ASSEMBLE1_DYNAM_3)
          DO nhs1=1,NHST(1)
            ny1=IABS(LGE(nhs1,1))
C$OMP       ATOMIC
            GR(ny1)=GR(ny1)+ER(nhs1)
            DO nhs2=1,NHST(2)
              ny2=IABS(LGE(nhs2,2))
              CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,NZ_GK_M,
     '          NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C$OMP         ATOMIC
              GK(nz)=GK(nz)+ES(nhs1,nhs2)
            ENDDO !nhs2
          ENDDO !nhs1
CC$OMP END CRITICAL(ASSEMBLE1_DYNAM_3)
        ELSE
C KAT 20Mar01: Can't branch out of critical section.
C              Atomic directive is probably better anyway.
CC$OMP CRITICAL(ASSEMBLE1_DYNAM_4)
          DO nhs1=1,NHST(1)
            ny1=IABS(LGE(nhs1,1))
            GR(ny1)=GR(ny1)+ER(nhs1)
            ny2=GETNYR(1,NPNY,nr,0,1,ny1,NYNE,NYNP)
            DO noy1=1,NONY(0,ny2,1)
              no1=NONY(noy1,ny2,1)
              co1=CONY(noy1,ny2,1)
              DO nhs2=1,NHST(2)
                ny3=IABS(LGE(nhs2,2))
                ny4=GETNYR(1,NPNY,nr,0,2,ny3,NYNE,NYNP)
                DO noy2=1,NONY(0,ny4,2)
                  no2=NONY(noy2,ny4,2)
                  co2=CONY(noy2,ny4,2)
                  CALL SPARSE(no1,no2,NOT(1,1,nr_gkk,nx),nzz,NZ_GKK_M,
     '              NZZT(1,nr_gkk,nx),ISC_GKK,ISR_GKK,SPARSEGKK(nx),
     '              ERROR,*9999)
                  IF(nzz.NE.0) THEN
C$OMP               ATOMIC
                    GKK(nzz)=GKK(nzz)+ES(nhs1,nhs2)*co1*co2
                  ENDIF
                ENDDO !noy2
                IF(NONY(0,ny4,2).EQ.0) THEN
                  CALL SPARSE(ny1,ny3,NYT(1,1,nx),nz,NZ_GK_M,
     '              NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C$OMP             ATOMIC
                  GK(nz)=GK(nz)+ES(nhs1,nhs2)
                ENDIF
              ENDDO !nhs2
            ENDDO !noy1
            IF(NONY(0,ny2,1).EQ.0) THEN
              DO nhs2=1,NHST(2)
                ny3=IABS(LGE(nhs2,2))
                CALL SPARSE(ny1,ny3,NYT(1,1,nx),nz,NZ_GK_M,
     '            NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C$OMP           ATOMIC
                GK(nz)=GK(nz)+ES(nhs1,nhs2)
              ENDDO !nhs2
            ENDIF
          ENDDO !nhs1
CC$OMP END CRITICAL(ASSEMBLE1_DYNAM_4)
        ENDIF

        IF(GQ_ASSEM) THEN
          CALL MELGE(LGE,NBH(1,1,ne),2,ne,NHE(ne),NHST,
     '      NPNE(1,1,ne),nr,NVHE(1,1,1,ne),nx,NYNE,NYNP,ERROR,*9999)
          IF(IWRIT4(nr,nx).GE.5) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(ASSEMBLE1_DYNAM_5)
            FORMAT='(/'' Element'',I5,'', Number of RHS vars: '','
     '        //'''NHST(1)='',I3,'', NHST(2)='',I3)'
            WRITE(OP_STRING,FORMAT) ne,NHST(1),NHST(2)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            FORMAT='('' LGE(1..,1): '',8(1X,I5),/:(13X,8(1X,I5)))'
            WRITE(OP_STRING,FORMAT) (LGE(nhs,1),nhs=1,NHST(1))
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            FORMAT='('' LGE(1..,2): '',8(1X,I5),/:(13X,8(1X,I5)))'
            WRITE(OP_STRING,FORMAT) (LGE(nhs,2),nhs=1,NHST(2))
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(ASSEMBLE1_DYNAM_5)
          ENDIF
C cpb 23/1/97 adding curvature corrections. EM stores curvature
C corrections
          DO nhs1=1,NHST(1)
            DO nhs2=1,NHST(2)
              ES(nhs1,nhs2)=0.0d0 !Should use EQ ????
              EM(nhs1,nhs2)=0.0d0
             ENDDO !nhs2
          ENDDO !nhs1
          CALL XPEQ30(ISC_GQ,ISR_GQ,LGE,NBH,NBJ,ne,NEELEM,NHP,NHST,
     '      NKJE,NPF,NP_INTERFACE,NPNE,nr,NRE,NVHE,NVJE,
     '      nx,NYNP,CURVCORRECT,ES,EM,GQ,PG,RG,SE,WG,XA,XE,XG1,XP,
     '      ERROR,*9999)

          IF(IWRIT4(nr,nx).GE.5) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(ASSEMBLE1_DYNAM_6)
            WRITE(OP_STRING,
     '        '(/'' Element flux matrix EQ:'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
c cpb 24/9/95 Use generic element stiffness matrix output
            CALL OPESTFMAT(NHST,IOOP,ES,ER,'EQ','ER',.TRUE.,.FALSE.,
     '        ERROR,*9999)
            IF(BEMCURVATURECORRECTION) THEN
              WRITE(OP_STRING,
     '          '(/'' Element stiffness matrix ES:'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              CALL OPESTFMAT(NHST,IOOP,EM,ER,'ES','ER',.TRUE.,
     '          .FALSE.,ERROR,*9999)
            ENDIF
CC$OMP END CRITICAL(ASSEMBLE1_DYNAM_6)
          ENDIF

C*** Assemble element flux matrix into global system.

          IF(CALC_GLOBAL(nr,nx)) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Atomic directive is probably better anyway.
CC$OMP CRITICAL(ASSEMBLE1_DYNAM_7)
            DO nhs1=1,NHST(1)
              ny1=IABS(LGE(nhs1,1))
              DO nhs2=1,NHST(2)
                ny2=IABS(LGE(nhs2,2))
                CALL SPARSE(ny1,ny2,NYT(1,2,nx),nz,NZ_GQ_M,
     '            NZT(2,nx),ISC_GQ,ISR_GQ,KTYP24,ERROR,*9999)
C$OMP           ATOMIC
                GQ(nz)=GQ(nz)+ES(nhs1,nhs2)
                IF(BEMCURVATURECORRECTION) THEN
                  CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,NZ_GK_M,
     '              NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C$OMP             ATOMIC
                  GK(nz)=GK(nz)-EM(nhs1,nhs2)
                ENDIF
              ENDDO !nhs2
            ENDDO !nhs1
CC$OMP END CRITICAL(ASSEMBLE1_DYNAM_7)
          ELSE
C KAT 20Mar01: Can't branch out of critical section.
C              Atomic directive is probably better anyway.
CC$OMP CRITICAL(ASSEMBLE1_DYNAM_8)
            DO nhs1=1,NHST(1)
              ny1=IABS(LGE(nhs1,1))
              ny2=GETNYR(1,NPNY,nr,0,1,ny1,NYNE,NYNP)
              DO noy1=1,NONY(0,ny2,1)
                no1=NONY(noy1,ny2,1)
                co1=CONY(noy1,ny2,1)
                DO nhs2=1,NHST(2)
                  ny3=IABS(LGE(nhs2,2))
                  ny4=GETNYR(1,NPNY,nr,0,2,ny3,NYNE,NYNP)
                  DO noy2=1,NONY(0,ny4,2)
                    no2=NONY(noy2,ny4,2)
                    co2=CONY(noy2,ny4,2)
                    CALL SPARSE(no1,no2,NOT(1,1,nr_gkk,nx),nzz,NZ_GKK_M,
     '                NZZT(1,nr_gkk,nx),ISC_GKK,ISR_GKK,SPARSEGKK(nx),
     '                ERROR,*9999)
C$OMP               ATOMIC
                    GKK(nzz)=GKK(nzz)-ES(nhs1,nhs2)*co1*co2
                    IF(BEMCURVATURECORRECTION) THEN
                      CALL SPARSE(no1,no2,NOT(1,1,nr_gkk,nx),nzz,
     '                  NZ_GKK_M,NZZT(1,nr,nx),ISC_GKK,ISR_GKK,
     '                  SPARSEGKK(nx),ERROR,*9999)
C$OMP                 ATOMIC
                      GKK(nzz)=GKK(nzz)-EM(nhs1,nhs2)*co1*co2
                    ENDIF
                  ENDDO !noy2
                  IF(NONY(0,ny4,2).EQ.0) THEN
                    CALL SPARSE(ny1,ny3,NYT(1,2,nx),nz,NZ_GQ_M,
     '                NZT(2,nx),ISC_GQ,ISR_GQ,KTYP24,ERROR,*9999)
C$OMP               ATOMIC
                    GQ(nz)=GQ(nz)+ES(nhs1,nhs2)
                    IF(BEMCURVATURECORRECTION) THEN
                      CALL SPARSE(ny1,ny3,NYT(1,1,nx),nz,NZ_GK_M,
     '                  NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C$OMP                 ATOMIC
                      GK(nz)=GK(nz)-EM(nhs1,nhs2)
                    ENDIF
                  ENDIF
                ENDDO !nhs2
              ENDDO !noy1
              IF(NONY(0,ny2,1).EQ.0) THEN
                DO nhs2=1,NHST(2)
                  ny3=IABS(LGE(nhs2,2))
                  CALL SPARSE(ny1,ny3,NYT(1,2,nx),nz,NZ_GQ_M,
     '              NZT(2,nx),ISC_GQ,ISR_GQ,KTYP24,ERROR,*9999)
C$OMP             ATOMIC
                  GQ(nz)=GQ(nz)+ES(nhs1,nhs2)
                  IF(BEMCURVATURECORRECTION) THEN
                    CALL SPARSE(ny1,ny3,NYT(1,1,nx),nz,NZ_GK_M,
     '                NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C$OMP               ATOMIC
                    GK(nz)=GK(nz)-EM(nhs1,nhs2)
                  ENDIF
                ENDDO !nhs2
              ENDIF
            ENDDO !nhs1
CC$OMP END CRITICAL(ASSEMBLE1_DYNAM_8)
          ENDIF
        ENDIF !GQ_ASSEM
      ENDIF !NW(ne,1).gt.0

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START3(1)
      IF(IWRIT4(nr,nx).GE.1) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(ASSEMBLE1_DYNAM_9)
        WRITE(OP_STRING,'(/'' CPU time for element '',I5,'
     '    //''' assembly: '',D11.4,'' s'')') ne,ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(ASSEMBLE1_DYNAM_9)
      ENDIF

      CALL EXITS('ASSEMBLE1_DYNAM')
      RETURN
 9999 CALL ERRORS('ASSEMBLE1_DYNAM',ERROR)
      CALL EXITS('ASSEMBLE1_DYNAM')
      RETURN 1
      END


