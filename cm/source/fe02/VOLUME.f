      SUBROUTINE VOLUME(NBJ,NELIST,NKJE,NPF,NPNE,nr,NVJE,NW,
     '  PG,RG,SE,VOL,WG,XA,XE,XG,XN,XP,ERROR,*)

C#### Subroutine: VOLUME
C###  Description:
C###    VOLUME calculates volume enclosed by a closed surface using
C###    Gaussian quad.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NELIST(0:NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  nr,NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3)
      REAL*8 PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM,NEM),
     '  VOL,WG(NGM,NBM),XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XG(NJM,NUM),XN(NJM,NGM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,ne,ng,nj,noelem
      REAL*8 DXIX(3,3),GL(3,3),GU(3,3),XGXN
      LOGICAL INTERFACE

      CALL ENTERS('VOLUME',*9999)

      VOL=0.0d0
      INTERFACE=.FALSE.
! AJP 30/8/99
!      DO noelem=1,NEELEM(0,nr)
!        ne=NEELEM(noelem,nr)
      DO noelem=1,NELIST(0)
        ne=NELIST(noelem)
        nb=NBJ(1,ne)
        CALL ASSERT(NIT(nb).EQ.2,'>> Only implemented for 2d elements',
     '    ERROR,*9999)
        CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '    nr,NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
        DO ng=1,NGT(nb)
          CALL XEXG(NBJ(1,ne),ng,nr,PG,XE,XG,ERROR,*9999)
          CALL NORMAL(ne,nr,NW,XG,XN(1,ng),INTERFACE,ERROR,*9999)
          CALL XGMG(0,3,nb,nr,DXIX,GL,GU,RG(ng),XG,ERROR,*9999)
          XGXN=0.0d0
          DO nj=1,NJT
            XGXN=XGXN+XG(nj,1)*XN(nj,ng)
          ENDDO !nj
          VOL=VOL+XGXN*RG(ng)*WG(ng,nb)
        ENDDO !ng
      ENDDO !ne

C LKC 19-SEP-2000 Will give negative volumes if the unit normal
C  vector is not pointing out of the volume.
C      VOL=VOL/3

      VOL=DABS(VOL/3.D0)

      CALL EXITS('VOLUME')
      RETURN
 9999 CALL ERRORS('VOLUME',ERROR)
      CALL EXITS('VOLUME')
      RETURN 1
      END


