      SUBROUTINE SGFIBR(INDEX,IBT,IDO,INP,ISEG,ISFIBR,iw,MXI,NAN,NBJ,ne,
     '  NKHE,NKJE,NPF,NPNE,nr,NVHE,NVJE,NW,nx,
     '  CSEG,CURVCORRECT,SE,TYPE,XA,XE,XG,XI,XP,ZA,ZE,ZP,DEFORMED,
     '  ERROR,*)

C#### Subroutine: SGFIBR
C###  Description:
C###    SGFIBR creates element fibre segment ISFIBR.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'map000.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INDEX,
     '  INP(NNM,NIM,NBFM),ISEG(*),ISFIBR,iw,MXI(2),
     '  NAN(NIM,NAM,NBFM),NBJ(NJM),ne,
     '  NKHE(NKM,NNM,NHM),NKJE(NKM,NNM,NJM),NPF(9),NPNE(NNM,NBFM),nr,
     '  NVHE(NNM,NBFM,NHM),NVJE(NNM,NBFM,NJM),
     '  NW,nx
      REAL*8 CURVCORRECT(2,2,NNM),SE(NSM,NBFM),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XG(NJM,NUM),XI(*),
     '  XP(NKM,NVM,NJM,NPM),ZA(NAM,NHM),ZE(NSM,NHM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*),TYPE*(*)
      LOGICAL DEFORMED
!     Local Variables
      INTEGER INDEX_OLD,nb,nc,nj,nk,nn,ns

      CALL ENTERS('SGFIBR',*9999)
      nc=1 !Temporary MPN 12-Nov-94
      CALL OPEN_SEGMENT(ISFIBR,ISEG,iw,'FIBR',INDEX,INDEX_OLD,
     '  ne,1,CSEG,ERROR,*9999)

      IF(iw.EQ.4) THEN
        PROJEC=MAP_PROJEC
        IF(PROJEC(1:2).EQ.'XI') THEN
          MXI1=MXI(1) !bottom left coords
          MXI2=MXI(2) !..for Xi map projection
        ENDIF
      ENDIF

      CALL XPXE(NBJ,NKJE,NPF,NPNE,nr,NVJE,
     '  SE,XA(1,1,ne),XE,XP,ERROR,*9999)
      IF(ITYP2(nr,nx).EQ.14.OR.ITYP2(nr,nx).EQ.15.OR.DEFORMED) THEN
        CALL ZPZE(NBJ,nc,NJ_LOC(NJL_GEOM,0,nr),NKHE,NPF,NPNE,
     '    nr,NVHE,NW,nx,CURVCORRECT,SE,ZA,ZE,ZP,ERROR,*9999)
        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
          nb=NBJ(nj)
          DO ns=1,NST(nb)
            XE(ns,nj)=ZE(ns,nj)
          ENDDO
        ENDDO
      ENDIF

      IF(iw.EQ.4.AND.ITYP10(nr).EQ.4.AND.PROJEC(1:6).EQ.'HAMMER') THEN
        nb=NBJ(1)                  !hammer proj with prolate coords
        DO nn=1,NNT(nb)               !set lamda=1 and zero derivs
          XE(1+(nn-1)*NKT(0,nb),1)=1.0D0
          DO nk=2,NKT(0,nb)
            XE(nk+(nn-1)*NKT(0,nb),1)=0.0D0
          ENDDO
        ENDDO
      ENDIF

      IF(TYPE(1:4).EQ.'AXIS') THEN
        CALL FIBRE2(INDEX,IBT,IDO,INP,iw,NAN,NBJ,nr,XE,XG,
     '    ERROR,*9999)
      ELSE IF(TYPE(1:5).EQ.'SHEET') THEN
        CALL FIBRE3(INDEX,IBT,IDO,INP,iw,NAN,NBJ,nr,XE,XG,ERROR,*9999)
      ELSE IF(TYPE(1:5).EQ.'HELIX') THEN
        CALL FIBRE4(INDEX,IBT,IDO,INP,iw,NAN,NBJ,nr,XE,XG,XI,
     '    ERROR,*9999)
      ENDIF

      CALL CLOSE_SEGMENT(ISFIBR,iw,ERROR,*9999)

      CALL EXITS('SGFIBR')
      RETURN
 9999 CALL ERRORS('SGFIBR',ERROR)
      CALL EXITS('SGFIBR')
      RETURN 1
      END


