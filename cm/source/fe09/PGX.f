      REAL*8 FUNCTION PGX(nb,nj,ns,DXIX,PG)

C#### Function: PGX
C###  Type: REAL*8
C###  Description:
C###    PGX calculates the first partial derivative of the ns-th
C###    term basis function PG with respect to the Xj-th coordinate.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER nb,nj,ns
      REAL*8 DXIX(3,3),PG(NSM,NUM)
!     Local Variables
      INTEGER ni,nu

      PGX=0.0D0
      DO ni=1,NIT(nb)
        nu=1+ni*(1+ni)/2
        PGX=PGX+PG(ns,nu)*DXIX(ni,nj)
      ENDDO

      RETURN
      END


