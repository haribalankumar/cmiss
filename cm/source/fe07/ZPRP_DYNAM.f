      SUBROUTINE ZPRP_DYNAM(IBT,IDO,INP,LGE,NAN,NBH,NBJ,NBJF,nc,ne,
     '  NFF,NGAP,NHE,NKEF,NKHE,NKJE,NNF,NPF,NPNE,NPNY,nr,NRE,NVHE,
     '  NVJE,NW,nx,NXI,NYNE,NYNP,CE,CG,CGE,CP,CURVCORRECT,FEXT,PG,RE,RG,
     '  SE,WG,XA,XE,XG,XP,YG,YP,ZA,ZA1,ZAA,ZE,ZE1,ZG,ZP,ZP1,ZPA,ERROR,*)

C#### Subroutine: ZPRP_DYNAM
C###  Description:
C###    ZPRP_DYNAM is the dynamic subroutine associted with ZPRP.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LGE(NHM*NSM,NRCM),NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),
     '  NBJ(NJM,NEM),NBJF(NJM,NFM),nc,ne,NFF(6,NEM),NGAP(NIM,NBM),
     '  NHE(NEM),NKEF(0:4,16,6,NBFM),NKJE(NKM,NNM,NJM,NEM),
     '  NKHE(NKM,NNM,NHM,NEM),NNF(0:17,6,NBFM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NPNY(0:6,NYM,0:NRCM),nr,NRE(NEM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),nx,NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),CGE(NMM,NGM,NEM),CP(NMM,NPM),
     '  CURVCORRECT(2,2,NNM,NEM),FEXT(NIFEXTM,NGM,NEM),
     '  PG(NSM,NUM,NGM,NBM),RE(NSM,NHM),RG(NGM),
     '  SE(NSM,NBFM,NEM),WG(NGM,NBM),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XG(NJM,NUM),XP(NKM,NVM,NJM,NPM),
     '  YG(NIYGM,NGM,NEM),YP(NYM,NIYM),
     '  ZA(NAM,NHM,NCM,NEM),ZA1(NAM,NHM,NCM,NEM),ZAA(NAM,NHM,NCM,NEM),
     '  ZE(NSM,NHM),ZE1(NSM,NHM),ZG(NHM,NUM),ZP(NKM,NVM,NHM,NPM,NCM),
     '  ZP1(NKM,NVM,NHM,NPM,NCM),ZPA(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,nh,NHST(2),nhs,nhx,np,ns,ny
      REAL*8 ZEA(NSM,NHM)
      CALL ENTERS('ZPRP_DYNAM',*9999)

      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP   CRITICAL(ZPRP_DYNAM_1)
        WRITE(OP_STRING,'(/'' ============='')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'( '' Element '',I5)') ne
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'( '' ============='')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP   END CRITICAL(ZPRP_DYNAM_1)
      ENDIF
      IF(NW(ne,1).GE.0) THEN
        CALL MELGE(LGE,NBH(1,1,ne),nc,ne,NHE(ne),NHST,
     '    NPNE(1,1,ne),nr,NVHE(1,1,1,ne),nx,NYNE,NYNP,ERROR,*9999)
        CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '    NPNE(1,1,ne),nr,NVJE(1,1,1,ne),SE(1,1,ne),
     '    XA(1,1,ne),XE,XP,ERROR,*9999)
        CALL ZPZE(NBH(1,1,ne),1,NHE(ne),NKHE(1,1,1,ne),
     '    NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),
     '    NW(ne,1),nx,CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '    ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
        IF (KTYP5I(nr).EQ.1) THEN !inertia
          CALL ZPZE(NBH(1,1,ne),1,NHE(ne),NKHE(1,1,1,ne),
     '      NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),
     '      NW(ne,1),nx,CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '      ZAA(1,1,1,ne),ZEA,ZPA,ERROR,*9999)
        ENDIF
        IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN !const vol
C         Put reference state for cavity from ZA1,ZP1 into ZE1
C         for ZERE55
          CALL ZPZE(NBH(1,1,ne),1,NHE(ne),NKHE(1,1,1,ne),
     '      NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),
     '      NW(ne,1),nx,CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '      ZA1(1,1,1,ne),ZE1,ZP1,ERROR,*9999)
        ENDIF

        nb=NBH(NH_LOC(1,nx),nc,ne)
        IF(ITYP1(nr,nx).EQ.3) THEN
          CALL CPCG(1,nb,NPNE(1,1,ne),nr,nx,
     '      CE(1,ne),CG,CGE(1,1,ne),CP,PG,ERROR,*9999)
          CALL ZERE30(NBH(1,nc,ne),NBJ(1,ne),NHE(ne),nr,nx,
     '      CG,PG,RE,WG,XE,XG,YG(1,1,ne),ZE,ZG,ERROR,*9999)
C         Apply scale factors if necessary
          DO nhx=1,NHE(ne)
            nh=nh_loc(nhx,nx)
            nb=NBH(nh,nc,ne)
            IF(NBI(nb).NE.1) THEN !scale factors not unit
              DO ns=1,NST(nb)
                RE(ns,nh)=RE(ns,nh)*SE(ns,nb,ne)
              ENDDO
            ENDIF
          ENDDO
        ELSE IF(ITYP1(nr,nx).EQ.4) THEN
          CALL CPCG(NW(ne,1),nb,NPNE(1,1,ne),nr,nx,
     '      CE(1,ne),CG,CGE(1,1,ne),CP,PG,ERROR,*9999)
          CALL ZERE40(NBH(1,nc,ne),NBJ(1,ne),ne,NHE(ne),
     '      NW(ne,1),nx,CE(1,ne),CG,PG,RE,SE,WG,XE,XG,
     '      ZE,ZG,ERROR,*9999)
        ELSE IF(ITYP1(nr,nx).EQ.5) THEN
          CALL CPCG(1,nb,NPNE(1,1,ne),nr,nx,
     '      CE(1,ne),CG,CGE(1,1,ne),CP,PG,ERROR,*9999)
          IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN !cnst vol
            CALL ZERE55(INP,NBH(1,nc,ne),ne,NHE(ne),nr,nx,
     '        CG,PG,RE,WG,ZE,ZE1,ZG,ERROR,*9999)
          ELSE
            CALL ZERE50(IBT,IDO,INP,NAN,
     '        NBH(1,nc,ne),NBJ(1,ne),NBJF,ne,NFF(1,ne),NGAP,NHE(ne),
     '        NKEF,NNF,NPNE,nr,NRE,NW(ne,1),nx,NXI,
     '        CE(1,ne),CG,CP,FEXT(1,1,ne),PG,RE,RG,SE,
     '        WG,XE,XG,YG(1,1,ne),ZE,ZEA,ZG,ERROR,*9999)
          ENDIF
        ENDIF

C KAT 20Mar01: Can't branch out of critical section.
C              Atomic directive is probably better anyway.
CC$OMP   CRITICAL(ZPRP_DYNAM_2)
        nhs=0
        DO nhx=1,NH_LOC(0,nx)
          nh=NH_LOC(nhx,nx)
          DO ns=1,NST(NBH(nh,nc,ne))+NAT(NBH(nh,nc,ne))
            nhs=nhs+1
            ny=IABS(LGE(nhs,1)) !row number
            IF(NPNY(0,ny,0).EQ.1) THEN
              np=NPNY(4,ny,0)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
            ENDIF
C$OMP       ATOMIC
            YP(ny,4)=YP(ny,4)+RE(ns,nh)
          ENDDO !ns
        ENDDO !nhx
CC$OMP   END CRITICAL(ZPRP_DYNAM_2)
      ENDIF !NW>0

      CALL EXITS('ZPRP_DYNAM')
      RETURN
 9999 CALL ERRORS('ZPRP_DYNAM',ERROR)
      CALL EXITS('ZPRP_DYNAM')
      RETURN 1
      END


