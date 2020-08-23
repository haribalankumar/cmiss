      SUBROUTINE SGSTRE(INDEX,IBT,IDO,INP,ISEG,ISSTRE,iw,
     '  NAN,NBH,NBJ,ne,NHE,NKHE,NKJE,NPF,NPNE,nr,
     '  NVHE,NVJE,NW,nx,
     '  CE,CG,CGE,CP,CSEG,CURVCORRECT,FEXT,PG,RG,SE,
     '  XA,XE,XG,XIG,XP,YG,ZA,ZE,ZG,ZP,ERROR,*)

C#### Subroutine: SGSTRE
C###  Description:
C###    SGSTRE creates element stress segments ISSTRE.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'grow00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp60.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INDEX,
     '  INP(NNM,NIM,NBFM),ISEG(*),ISSTRE(NEM,NGRSEGM),
     '  iw,NAN(NIM,NAM,NBFM),
     '  NBH(NHM,NCM),NBJ(NJM),ne,NHE,
     '  NKHE(NKM,NNM,NHM),NKJE(NKM,NNM,NJM),NPF(9),NPNE(NNM,NBFM),nr,
     '  NVHE(NNM,NBFM,NHM),NVJE(NNM,NBFM,NJM),NW,nx
      REAL*8 CE(NMM),CG(NMM,NGM),CGE(NMM,NGM),
     '  CP(NMM,NPM),CURVCORRECT(2,2,NNM),
     '  FEXT(NIFEXTM,NGM),
     '  PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),XIG(NIM,NGM,NBM),
     '  XP(NKM,NVM,NJM,NPM),YG(NIYGM,NGM),ZA(NAM,NHM),ZE(NSM,NHM),
     '  ZG(NHM,NUM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER INDEX_OLD,nc,ng

      CALL ENTERS('SGSTRE',*9999)
      nc=1   !Temporary

      CALL OPEN_SEGMENT(ISSTRE(ne,NTSTRE),ISEG,iw,'STRE',INDEX,
     '  INDEX_OLD,ne,1,CSEG,ERROR,*9999)

      CALL XPXE(NBJ,NKJE,NPF,NPNE,nr,NVJE,
     '  SE,XA(1,1,ne),XE,XP,ERROR,*9999)
      CALL ZPZE(NBH,nc,NHE,NKHE,NPF,NPNE,nr,NVHE,NW,nx,
     '  CURVCORRECT,SE,ZA,ZE,ZP,ERROR,*9999)
      IF(ITYP1(nr,nx).EQ.5) THEN !finite elasticity
        CALL CPCG(1,NBH(NH_LOC(1,nx),nc),NPNE,nr,nx,CE,CG,CGE,CP,PG,
     '    ERROR,*9999)
      ELSE !all other problems
        CALL CPCG(NW,NBH(NH_LOC(1,nx),nc),NPNE,nr,nx,CE,CG,CGE,CP,PG,
     '    ERROR,*9999)
        IF(KTYP60.EQ.1) THEN !Growth law: so get density from YG(5)
          DO ng=1,NGT(NBJ(1))
            CG(5,ng)=YG(5,ng) !density
            CG(1,ng)=CG(1,ng)*CG(5,ng)**2 !Young's mod from density
          ENDDO
        ELSE IF(KTYP60.EQ.2) THEN !Growth law: so get density from YG(5)
          DO ng=1,NGT(NBJ(1))
            CG(5,ng)=YG(5,ng) !density
            CG(1,ng)=CG(1,ng)*(CG(5,ng)/GROW1)**2 !Young's mod from density
          ENDDO !ng
        ENDIF
      ENDIF
      CALL STRESS1(INDEX,IBT,IDO,INP,iw,NAN,NBH,NBJ,ne,NHE,NPNE,nr,
     '  NW,nx,CE,CG,CP,FEXT,PG,RG,XE,XG,XIG,YG,ZE,ZG,ERROR,*9999)

      CALL CLOSE_SEGMENT(ISSTRE(ne,NTSTRE),iw,ERROR,*9999)

      CALL EXITS('SGSTRE')
      RETURN
 9999 CALL ERRORS('SGSTRE',ERROR)
      CALL EXITS('SGSTRE')
      RETURN 1
      END


