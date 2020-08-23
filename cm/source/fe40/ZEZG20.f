      SUBROUTINE ZEZG20(NBH,ng,NHE,nx,PG,ZE,ZGG,ERROR,*)

C#### Subroutine: ZEZG20
C###  Description:
C###    ZEZG20 evaluates the Gauss point array ZGG(nu,nh) from the
C###    element array ZE(ns,nhx) at the current Gauss point ng,
C###    where nu=1,NUT(nb).

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM),ng,NHE,nx
      REAL*8 PG(NSM,NUM,NGM,NBM),ZE(NSM,NHM),ZGG(10,*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,nh,nhx,ns,nu
      REAL*8 SUM

      CALL ENTERS('ZEZG20',*9999)
      DO nhx=1,NHE
        nh=NH_LOC(nhx,nx) !GMH 20-Apr-95
        nb=NBH(nh)
        DO nu=1,NUT(nb)
          SUM=0.d0
          DO ns=1,NST(nb)
            SUM=SUM+PG(ns,nu,ng,nb)*ZE(ns,nhx)
          ENDDO
          ZGG(nu,nh)=SUM
        ENDDO
      ENDDO

      CALL EXITS('ZEZG20')
      RETURN
 9999 CALL ERRORS('ZEZG20',ERROR)
      CALL EXITS('ZEZG20')
      RETURN 1
      END


