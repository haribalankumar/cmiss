      REAL*8 FUNCTION NODE_DIST(nk1,nk2,np1,np2,nv1,nv2,XP)

C#### Function: NODE_DIST
C###  Type: REAL*8
C###  Description:
C###    Determines the distance between 2 node points from
C###    CMISS data array XP
C###  See-Also: DATA_DIST

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'

!     Parameter List
      INTEGER nk1,nk2,np1,np2,nv1,nv2
      REAL*8 XP(NKM,NVM,NJM,NPM)

!     Local Variables
      INTEGER nj

      NODE_DIST=0.0d0
      DO nj=1,NJT
        NODE_DIST=NODE_DIST+ (XP(nk2,nv2,nj,np2)-XP(nk1,nv1,nj,np1))**2
      ENDDO
      NODE_DIST=DSQRT(NODE_DIST)

      RETURN
      END


