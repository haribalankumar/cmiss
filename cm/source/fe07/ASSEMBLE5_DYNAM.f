      SUBROUTINE ASSEMBLE5_DYNAM(IBT,IDO,INP,ISC_GK,ISR_GK,LGE,NAN,NBH,
     '  NBJ,NBJF,ne,NFF,NGAP,NHE,NKEF,NKHE,NKJE,NMNO,NNF,NPF,NPNE,NPNY,
     '  nr,NRE,NVHE,NVJE,NW,nx,NXI,NYNE,NYNP,CE,CG,CGE,CP,CURVCORRECT,
     '  D_RE,D_RI3,D_TG,D_ZG,ES,FEXT,GK,PG,RE1,RE2,RG,SE,WG,XA,XE,XG,
     '  XP,YG,ZA,ZA1,ZAA,ZE,ZE1,ZG,ZG1,ZP,ZP1,ZPA,FIX,ERROR,*)

C#### Subroutine: ASSEMBLE5_DYNAM
C###  Description:
C###    ASSEMBLE5_DYNAM is the dynamic subroutine associted with
C###    ASSEMBLE5.


      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'mxch.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISC_GK(NISC_GKM),ISR_GK(NISR_GKM),LGE(NHM*NSM,NRCM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  ne,NFF(6,NEM),NGAP(NIM,NBM),NHE(NEM),NKEF(0:4,16,6,NBFM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NMNO(1:2,0:NOPM),NNF(0:17,6,NBFM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NPNY(0:6,NYM,0:NRCM),
     '  nr,NRE(NEM),NVHE(NNM,NBFM,NHM,NEM),NVJE(NNM,NBFM,NJM,NEM),
     '  NW(NEM,3),nx,NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),CGE(NMM,NGM,NEM),CP(NMM,NPM),
     '  CURVCORRECT(2,2,NNM,NEM),D_RE(NSM,NHM,NOPM),D_RI3(NHM*NSM),
     '  D_TG(3,3,NHM*NSM),D_ZG(NHM,NUM,NHM*NSM),ES(NHM*NSM,NHM*NSM),
     '  FEXT(NIFEXTM,NGM,NEM),GK(NZ_GK_M),PG(NSM,NUM,NGM,NBM),
     '  RE1(NSM,NHM),RE2(NSM,NHM),RG(NGM),SE(NSM,NBFM,NEM),
     '  WG(NGM,NBM),XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XG(NJM,NUM),XP(NKM,NVM,NJM,NPM),YG(NIYGM,NGM,NEM),
     '  ZA(NAM,NHM,NCM,NEM),ZA1(NAM,NHM,NCM,NEM),ZAA(NAM,NHM,NCM,NEM),
     '  ZE(NSM,NHM),ZE1(NSM,NHM),ZG(NHM,NUM),ZG1(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM),ZP1(NKM,NVM,NHM,NPM,NCM),
     '  ZPA(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM)
!     Local Variables
      INTEGER GETNYR,NHST(2),nhs,nhs1,nhs2,ny,ny1,ny2,nz
      REAL*8 ZEA(NSM,NHM)
      CHARACTER FORMAT*500
      LOGICAL ELEM

      CALL ENTERS('ASSEMBLE5_DYNAM',*9999)

      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP   CRITICAL(ASSEMBLE5_DYNAM_1)
        WRITE(OP_STRING,'(/'' ============='')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'( '' Element '',I5)') ne
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'( '' ============='')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP   END CRITICAL(ASSEMBLE5_DYNAM_1)
      ENDIF

      IF(NW(ne,1).GT.0) THEN

C ***   READ THE NOTE BELOW ***
        CALL MELGE(LGE,NBH(1,1,ne),1,ne,NHE(ne),NHST,
     '    NPNE(1,1,ne),nr,NVHE(1,1,1,ne),nx,NYNE,NYNP,
     '    ERROR,*9999)
C ***   NOTE: If anyone changes the parameter list of MELGE then either
C ***         update fe_parallel.f by:
C ***         1. passing the new arrays etc to slave process thru
C ***            sockets;
C ***     and 2. changing call to MELGE  (or 3. see MPN)

        IF(IWRIT4(nr,nx).GE.5) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP     CRITICAL(ASSEMBLE5_DYNAM_2)
          FORMAT='(/'' Element'',I5,'', Number of vars: '','
     '      //'''NHST(1)='',I3,'', NHST(2)='',I3)'
          WRITE(OP_STRING,FORMAT) ne,NHST(1),NHST(2)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          FORMAT='('' LGE(1..,1): '',10I5,:/(13X,10I5))'
          WRITE(OP_STRING,FORMAT) (LGE(nhs,1),nhs=1,NHST(1))
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          FORMAT='('' LGE(1..,2): '',10I5,:/(13X,10I5))'
          WRITE(OP_STRING,FORMAT) (LGE(nhs,2),nhs=1,NHST(2))
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
CC$OMP     END CRITICAL(ASSEMBLE5_DYNAM_2)
        ENDIF

        ELEM=.FALSE.
        DO nhs2=1,NHST(2) !loop over local variables
          ny=IABS(LGE(nhs2,2)) !local variable number
          ny1=GETNYR(1,NPNY,nr,0,2,ny,NYNE,NYNP) !global var #
          IF(.NOT.FIX(ny1,1)) ELEM=.TRUE.
          DO nhs1=1,NHST(1) !loop over rows
            ES(nhs1,nhs2)=0.0d0
          ENDDO !nhs1
        ENDDO !nhs2

        IF(ELEM) THEN
          CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '      NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '      SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
          CALL ZPZE(NBH(1,1,ne),1,NHE(ne),NKHE(1,1,1,ne),
     '      NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1),nx,
     '      CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),
     '      ZE,ZP,ERROR,*9999)
          IF (KTYP5I(nr).EQ.1) THEN !inertia
            CALL ZPZE(NBH(1,1,ne),1,NHE(ne),NKHE(1,1,1,ne),
     '        NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1),nx,
     '        CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZAA(1,1,1,ne),
     '        ZEA,ZPA,ERROR,*9999) 
          ENDIF         
          IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN
C           Put reference state for cavity from ZA1,ZP1 into ZE1
C           for ZERE55 (constant volume constraint)
            CALL ZPZE(NBH(1,1,ne),1,NHE(ne),NKHE(1,1,1,ne),
     '        NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),
     '        NW(ne,1),nx,CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '        ZA1(1,1,1,ne),ZE1,ZP1,ERROR,*9999)
          ENDIF

C ***     READ THE NOTE BELOW ***
          CALL ZEES(IBT,IDO,INP,LGE,NAN,NBH(1,1,ne),NBJ(1,ne),
     '      NBJF,ne,NFF(1,ne),NGAP,NHE(ne),NKEF,NMNO,
     '      NNF,NPNE,NPNY,nr,NRE,NW(ne,1),nx,NXI,NYNE,NYNP,
     '      CE(1,ne),CG,CGE(1,1,ne),
     '      CP,D_RE,D_RI3,D_TG,D_ZG,ES,FEXT(1,1,ne),
     '      PG,RE1,RE2,RG,SE,WG,XE,XG,YG(1,1,ne),ZE,ZE1,ZEA,ZG,ZG1,
     '      FIX,ERROR,*9999)
C ***     NOTE: If anyone changes the parameter list of ZEES
C ***           then either update fe_parallel.f by:
C ***           1. passing the new arrays etc to slave process
C ***              thru sockets;
C ***       and 2. changing call to ZEES  (or 3. see MPN)

          IF(IWRIT4(nr,nx).GE.5) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP       CRITICAL(ASSEMBLE5_DYNAM_3)
            WRITE(OP_STRING,
     '        '(/'' Element '',I5,'' stiffness matrix ES:'')') ne
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            CALL OPESTFMAT(NHST,IOOP,ES,ES,'ES','ER',
     '        .TRUE.,.FALSE.,ERROR,*9999)
CC$OMP       END CRITICAL(ASSEMBLE5_DYNAM_3)
          ENDIF

C***      Assemble element stiffness matrix into global system.
C KAT 20Mar01: Can't branch out of critical section.
C              Atomic directive is probably better anyway.
CC$OMP     CRITICAL(ASSEMBLE5_DYNAM_4)
          DO nhs1=1,NHST(1)
            ny1=IABS(LGE(nhs1,1))
            DO nhs2=1,NHST(2)
              ny2=IABS(LGE(nhs2,2))
C ***         Note ES contains glob derivs (no correction by SE)
              CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,NZ_GK_M,
     '          NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C$OMP         ATOMIC
              GK(nz)=GK(nz)+ES(nhs1,nhs2)
            ENDDO !nhs2
          ENDDO !nhs1
CC$OMP     END CRITICAL(ASSEMBLE5_DYNAM_4)
        ENDIF
      ENDIF

      CALL EXITS('ASSEMBLE5_DYNAM')
      RETURN
 9999 CALL ERRORS('ASSEMBLE5_DYNAM',ERROR)
      CALL EXITS('ASSEMBLE5_DYNAM')
      RETURN 1
      END


