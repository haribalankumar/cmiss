      SUBROUTINE ZPZE_CAVITY(NBH,nc,NKJE,NPNE,NPNODE,
     '  nr,NVJE,
     '  SE,XP,ZE,ZP,ERROR,*)

C#### Subroutine: ZPZE_CAVITY
C###  Description:
C###    ZPZE_CAVITY transfers global node parameters ZP and XP
C###    to element node parameters ZE(ns,nhx). This avoids having
C###    to define equations for the cavity meshes when we are only
C###    interested in the deformed geometry that can be obtained
C###    from the interface nodes.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM,NCM),nc,NKJE(NKM,NNM,NJM),
     '  NPNE(NNM,NBFM),NPNODE(0:NP_R_M,0:NRM),
     '  nr,NVJE(NNM,NBFM,NJM)
      REAL*8 SE(NSM,NBFM),XP(NKM,NVM,NJM,NPM),
     '  ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER MK,nb,nh,nhx,nk,nn,np,ns,nv

      CALL ENTERS('ZPZE_CAVITY',*9999)

      DO nh=1,3
        nb=NBH(nh,1)
        ns=0
        DO nn=1,NNT(nb)
          np=NPNE(nn,nb)
          nv=NVJE(nn,nb,nh)
          DO mk=1,NKT(nn,nb)
            nk=NKJE(mk,nn,nh)
            ns=ns+1
            IF((mk.EQ.NKT(0,nb)).AND.
     '        ((KTYP93(nc,nr).EQ.1).AND.(NKT(0,nb).GT.3))) THEN
              ZE(ns,nhx)=0.0d0
            ELSE
              IF(np.LE.NPNODE(0,1)) THEN
                ZE(ns,nh)=ZP(nk,nv,nh,np,nc)*SE(ns,nb)
              ELSE
                ZE(ns,nh)=XP(nk,nv,nh,np)*SE(ns,nb)
              ENDIF
            ENDIF
          ENDDO !mk
        ENDDO !nn
      ENDDO !nhx

      CALL EXITS('ZPZE_CAVITY')
      RETURN
 9999 CALL ERRORS('ZPZE_CAVITY',ERROR)
      CALL EXITS('ZPZE_CAVITY')
      RETURN 1
      END


