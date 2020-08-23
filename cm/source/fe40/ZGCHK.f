      SUBROUTINE ZGCHK(NBH,NHE,nx,ZGG,ERROR,*)

C#### Subroutine: ZGCHK
C###  Description:
C###    ZGCHK checks for ZGG(nu,nh) values to make sure they are not too
C###    small to cause UNDERFLOW errors in SUBROUTINE GGRRM.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM),NHE,nx
      REAL*8 ZGG(10,*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,nh,nhx,nu

      CALL ENTERS('ZGCHK',*9999)
      DO nhx=1,NHE
        nh=NH_LOC(nhx,nx) !GMH 20-Apr-95
        nb=NBH(nh)
        DO nu=1,NUT(nb)
          IF(DABS(ZGG(nu,nh)).LT.1.D-20) ZGG(nu,nh)=0.d0
        ENDDO
      ENDDO

      CALL EXITS('ZGCHK')
      RETURN
 9999 CALL ERRORS('ZGCHK',ERROR)
      CALL EXITS('ZGCHK')
      RETURN 1
      END


