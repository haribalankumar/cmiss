      SUBROUTINE X3XG(NBJ,D3PG,XE,X3G,ERROR,*)

C#### Subroutine: X3XG
C###  Description:
C###    X3XG evaluates Gauss point array X3G of 3rd derivs of X wrt Xi
C###    coordinates from element node array XE.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBJ(NJM)
      REAL*8 D3PG(16,4,25,*),XE(NSM,NJM),X3G(4,3,*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,nd,ng,nj,ns
      REAL*8 SUM

      CALL ENTERS('X3XG',*9999)
      DO nj=1,NJ_LOC(NJL_GEOM,0,1) ! 1 is temporary RGB
        nb=NBJ(nj)
        DO nd=1,4
          DO ng=1,NGT(nb)
            SUM=0.d0
            DO ns=1,NST(nb)
              SUM=SUM+D3PG(ns,nd,ng,nb)*XE(ns,nj)
            ENDDO
            X3G(nd,nj,ng)=SUM
          ENDDO
        ENDDO
      ENDDO

      CALL EXITS('X3XG')
      RETURN
 9999 CALL ERRORS('X3XG',ERROR)
      CALL EXITS('X3XG')
      RETURN 1
      END


