      SUBROUTINE ZPZE20(NBH,NBJ,NHE,NPNE,NVHE,nx,SE,ZE,ZP,
     '  ERROR,*)

C#### Subroutine: ZPZE20
C###  Description:
C###    ZPZE20 transfers global node parameters ZP(nk,nv,nh,np,nc) to
C###    element node parameters ZE(ns,nhx).

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM),NBJ(NJM),NHE,NPNE(NNM,NBFM),NVHE(NNM,NBFM,NHM),
     '  nx
      REAL*8 SE(NSM,NBFM),ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,nc,nh,nhx,nk,nn,NNK,NNS,np,ns,nv
      REAL*8 SEN

      CALL ENTERS('ZPZE20',*9999)
      nc=1 !Temporary AJP 18-12-91
      DO nhx=1,NHE
        nh=NH_LOC(nhx,nx) !GMH 20-Apr-95
        nb=NBH(nh)
        ns=0
        DO nn=1,NNT(nb)
          np=NPNE(nn,nb)
          nv=NVHE(nn,nb,nh)
          IF(nh.LE.2) THEN
            NNK=nh+1
            NNS=(nn-1)*NKT(0,NBJ(1))+NNK
            SEN=SE(NNS,NBJ(1))
          ELSE
            SEN=1.d0
          ENDIF
          DO nk=1,NKT(0,nb)
            ns=NS+1
            ZE(ns,nhx)=ZP(nk,nv,nh,np,nc)*SE(ns,nb)*SEN
          ENDDO
        ENDDO
      ENDDO

      CALL EXITS('ZPZE20')
      RETURN
 9999 CALL ERRORS('ZPZE20',ERROR)
      CALL EXITS('ZPZE20')
      RETURN 1
      END

