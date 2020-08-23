      SUBROUTINE XPXE_POINTS(from_field,NBJ,NKJE,NPF,NPNE,nr,NVJE,SE,
     &   XA,XE,XP,ERROR,*)

C#### Subroutine: XPXE_POINTS
C###  Description:
C###    XPXE_POINTS transfers field values to XE in place of geometry.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER from_field(3),NBJ(NJM),NKJE(NKM,NNM,NJM),NPF(9),
     '  NPNE(NNM,NBFM),nr,NVJE(NNM,NBFM,NJM)
      REAL*8 SE(NSM,NBFM),XA(NAM,NJM),XE(NSM,NJM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER mk,nb,nj,nj_field,nk,nn,np,ns,nv

      CALL ENTERS('XPXE_POINTS',*9999)

      DO nj=1,NJT 
        nj_field=from_field(nj)
        nb=NBJ(nj_field)
        IF(nb.GT.0) THEN
          ns=0
          IF(NNT(nb).GT.0) THEN
!           write(*,*) 'nnt(nb)=',NNT(nb)
            DO nn=1,NNT(nb)
              np=NPNE(nn,nb)
              nv=NVJE(nn,nb,nj_field)
!              write(*,*) 'nkt=',NKT(nn,nb)
              DO mk=1,NKT(nn,nb)
                ns=ns+1
                nk=NKJE(mk,nn,nj_field)
                XE(ns,nj)=XP(nk,nv,nj_field,np)*SE(ns,nb)
!           write(*,*) 'xP=',XP(nk,nv,nj_field,np),'np=',np
!           write(*,*) 'ns,nj,xe=',ns,nj,XE(ns,nj),'np=',np
              ENDDO
            ENDDO
          ENDIF
        ENDIF !nb.GT.0
      ENDDO !nj

      CALL EXITS('XPXE_POINTS')
      RETURN
 9999 CALL ERRORS('XPXE_POINTS',ERROR)
      CALL EXITS('XPXE_POINTS')
      RETURN 1
      END


